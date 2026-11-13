# The AGS solver of NLopt 2.11.0 (GN_AGS, src/algs/ags: ags.cc, solver.cc,
# evolvent.cc, by V. Sovrasov): the index method of Strongin and Sergeyev on
# the Peano type evolvent of the search box, with inequality constraints
# handled by the index scheme. The algorithm parameters are the global C
# variables of ags.cc (ags_r = 3, eps_res = 0.001, evolvent_density = 12,
# ags_eps = 0, ags_refine_loc = 0), which the NLopt API does not expose; the
# local refinement (HookeJeevesOptimizer) is off by those defaults.
#
# The C++ is compiled with floating point contraction on: the transform of
# the evolvent to the search box is fused (.fma); the other contracted
# products are exact (by -2, 2 and 0.5). The queue is the std::priority_queue
# of libc++ (.nlopt_heap_push and .nlopt_heap_pop of
# R/optimization-nlopt-stogo.R), whose order on ties of the characteristic R
# decides the next interval.
#
# Copyright (C) 2018 Sovrasov V., released under the MIT license.

# pow of the C library: R's `^` takes x * x for the exponent 2, which can
# differ from pow(x, 2.0) in the last bit; pweibull evaluates -pow(x, y)
# for x > 0 with scale 1.
.nlopt_libm_pow <- function(x, y) {
  if (is.nan(x) || is.nan(y)) return(x + y)
  if (y == 0) return(1)
  if (x > 0) return(-stats::pweibull(x, y, 1, lower.tail = FALSE, log.p = TRUE))
  x^y
}

# mapd of evolvent.cc (key 1): the point of the evolvent of density m of the
# cube [-1/2, 1/2]^n at x in [0, 1]
.nlopt_ags_mapd <- function(x, m, n) {
  n1 <- n - 1L
  nexp <- 2L^n
  d <- x
  r <- 0.5
  it <- 0L
  l <- 0L
  iw <- rep(1, n)
  y <- numeric(n)
  iu <- integer(n); iv <- integer(n)
  for (j in seq_len(m)) {
    iq <- 0L
    if (x == 1) {
      is <- nexp - 1L
      d <- 0
    } else {
      d <- d * nexp
      is <- as.integer(trunc(d))
      d <- d - is
    }
    # node
    if (is == 0L) {
      l <- n1
      iu[] <- -1L; iv[] <- -1L
    } else if (is == nexp - 1L) {
      l <- n1
      iu[] <- -1L; iv[] <- -1L
      iu[1L] <- 1L; iv[1L] <- 1L
      iv[n] <- 1L
    } else {
      iff <- nexp
      k1 <- -1L
      ii <- is
      for (i in 0:n1) {
        iff <- iff %/% 2L
        if (ii >= iff) {
          if (ii == iff && ii != 1L) { l <- i; iq <- -1L }
          ii <- ii - iff
          k2 <- 1L
        } else {
          k2 <- -1L
          if (ii == iff - 1L && ii != 0L) { l <- i; iq <- 1L }
        }
        jj <- -k1 * k2
        iv[i + 1L] <- jj
        iu[i + 1L] <- jj
        k1 <- k2
      }
      iv[l + 1L] <- iv[l + 1L] * iq
      iv[n] <- -iv[n]
    }
    t <- iu[1L]; iu[1L] <- iu[it + 1L]; iu[it + 1L] <- t
    t <- iv[1L]; iv[1L] <- iv[it + 1L]; iv[it + 1L] <- t
    if (l == 0L) l <- it else if (l == it) l <- 0L
    r <- r * 0.5
    it <- l
    iu <- iu * iw
    iw <- -iv * iw
    y <- r * iu + y
  }
  y
}

# ags_minimize with NLPSolver::Solve: `fobj` is the objective, `fcons` the
# list of scalar constraints; returns list(ret, x, minf)
.nlopt_ags <- function(fobj, fcons, lb, ub, stop, r = 3, eps = 0, eps_r = 0.001, density = 12L) {
  n <- length(lb)
  m <- length(fcons)
  iters_limit <- if (stop$maxeval != 0) stop$maxeval else .Machine$integer.max
  if (!(max(ub - lb, 0) > 0)) return(list(ret = "FAILURE", x = lb, minf = Inf))
  rho <- ub - lb
  shift <- 0.5 * (lb + ub)
  get_image <- function(x) {
    y <- if (n != 1L) .nlopt_ags_mapd(x, density, n) else x - 0.5
    .fma(rho, y, shift)
  }
  # the trials: x, y, g (constraints then objective), idx (-1 for the ends)
  cap <- 64L
  TX <- numeric(cap); TY <- matrix(0, cap, n); TG <- matrix(0, cap, m + 1L); TI <- integer(cap)
  nt <- 0L
  add_trial <- function(x, y, g, idx) {
    if (nt == length(TX)) {
      k <- length(TX)
      TX <<- c(TX, numeric(k)); TY <<- rbind(TY, matrix(0, k, n))
      TG <<- rbind(TG, matrix(0, k, m + 1L)); TI <<- c(TI, integer(k))
    }
    nt <<- nt + 1L
    TX[nt] <<- x; TY[nt, ] <<- y; TG[nt, ] <<- g; TI[nt] <<- idx
    nt
  }
  # the intervals: left and right trial, delta, characteristic R
  IL <- integer(cap); IR <- integer(cap); DL <- numeric(cap); RV <- numeric(cap)
  ni <- 0L
  add_interval <- function(pl, pr) {
    if (ni == length(IL)) {
      k <- length(IL)
      IL <<- c(IL, integer(k)); IR <<- c(IR, integer(k)); DL <<- c(DL, numeric(k)); RV <<- c(RV, numeric(k))
    }
    ni <<- ni + 1L
    IL[ni] <<- pl; IR[ni] <<- pr
    ni
  }
  sx <- integer(0) # the std::set ordered by pl.x
  hq <- integer(0) # the heap of the std::priority_queue
  H <- rep(1, m + 1L)
  Z <- rep(.Machine$double.xmax, m + 1L)
  cnt <- numeric(m + 1L)
  iterations <- 0
  need_refill <- FALSE
  need_stop <- FALSE
  min_delta <- .Machine$double.xmax
  max_idx <- -1L
  opt_y <- numeric(n); opt_g <- numeric(m + 1L); opt_idx <- -1L
  ret_code <- "SUCCESS"
  # the next point (numPoints = 1) and the interval it splits
  np_x <- 0; np_y <- numeric(n); np_g <- numeric(m + 1L); np_idx <- 0L
  next_interval <- 0L

  q_push <- function(id) hq <<- .nlopt_heap_push(hq, id, RV)
  q_pop <- function() {
    p <- .nlopt_heap_pop(hq, RV)
    hq <<- p$h
    p$top
  }
  calculate_r <- function(id) {
    pl <- IL[id]; pr <- IR[id]
    ipl <- TI[pl]; ipr <- TI[pr]
    if (ipl == ipr) {
      v <- ipr + 1L
      rh <- r * H[v]
      d <- DL[id]
      q <- (TG[pr, v] - TG[pl, v]) / rh
      d + q * q / d - 2 * (TG[pr, v] + TG[pl, v] - 2 * Z[v]) / rh
    } else if (ipl < ipr) {
      v <- ipr + 1L
      2 * DL[id] - 4 * (TG[pr, v] - Z[v]) / (r * H[v])
    } else {
      v <- ipl + 1L
      2 * DL[id] - 4 * (TG[pl, v] - Z[v]) / (r * H[v])
    }
  }
  next_coordinate <- function(id) {
    pl <- IL[id]; pr <- IR[id]
    if (TI[pr] == TI[pl]) {
      v <- TI[pr] + 1L
      dg <- TG[pr, v] - TG[pl, v]
      sg <- if (isTRUE(dg > 0)) 1 else -1
      0.5 * (TX[pr] + TX[pl]) - 0.5 * sg * .nlopt_libm_pow(abs(dg) / H[v], n) / r
    } else {
      0.5 * (TX[pr] + TX[pl])
    }
  }
  update_h <- function(value, index) {
    v <- index + 1L
    if (isTRUE(value > H[v]) || (H[v] == 1 && isTRUE(value > 1e-12))) {
      H[v] <<- value
      need_refill <<- TRUE
    }
  }
  update_all_h <- function(pos) {
    id <- sx[pos]
    pl <- IL[id]
    ipl <- TI[pl]
    if (ipl < 0L) return(invisible())
    if (ipl == TI[IR[id]]) {
      update_h(abs(TG[IR[id], ipl + 1L] - TG[pl, ipl + 1L]) / DL[id], ipl)
    } else {
      N <- length(sx)
      rp <- pos + 1L
      while (rp <= N && TI[IL[sx[rp]]] < ipl) rp <- rp + 1L
      if (rp <= N && TI[IL[sx[rp]]] >= ipl) {
        t <- IL[sx[rp]]
        update_h(abs(TG[t, ipl + 1L] - TG[pl, ipl + 1L]) / (TX[t] - TX[pl])^(1 / n), ipl)
      }
      lp <- pos - 1L
      while (lp != 1L && TI[IL[sx[lp]]] < ipl) lp <- lp - 1L
      if (lp != 1L && TI[IL[sx[lp]]] >= ipl) {
        t <- IL[sx[lp]]
        update_h(abs(TG[t, ipl + 1L] - TG[pl, ipl + 1L]) / (TX[pl] - TX[t])^(1 / n), ipl)
      }
    }
  }
  # std::set::insert; returns the position, or 0 when pl.x is present
  set_insert <- function(id) {
    xs <- TX[IL[sx]]
    x <- TX[IL[id]]
    k <- sum(xs < x)
    if (k < length(xs) && xs[k + 1L] == x) return(0L)
    sx <<- append(sx, id, after = k)
    k + 1L
  }
  refill_queue <- function() {
    hq <<- integer(0)
    for (id in sx) {
      RV[id] <<- calculate_r(id)
      q_push(id)
    }
    need_refill <<- FALSE
  }
  make_trials <- function() {
    idx <- 0L
    while (idx < m) {
      np_idx <<- idx
      val <- fcons[[idx + 1L]](np_y)
      cnt[idx + 1L] <<- cnt[idx + 1L] + 1
      np_g[idx + 1L] <<- val
      if (isTRUE(val > 0)) break
      idx <- idx + 1L
    }
    if (idx > max_idx) {
      max_idx <<- idx
      for (j in seq_len(max_idx)) Z[j] <<- -eps_r * H[j]
      need_refill <<- TRUE
    }
    if (idx == m) {
      cnt[idx + 1L] <<- cnt[idx + 1L] + 1
      np_idx <<- idx
      np_g[idx + 1L] <<- fobj(np_y)
    }
    if (np_idx == max_idx && isTRUE(np_g[max_idx + 1L] < Z[max_idx + 1L])) {
      Z[max_idx + 1L] <<- np_g[max_idx + 1L]
      need_refill <<- TRUE
    }
  }
  estimate_optimum <- function() {
    if (opt_idx < np_idx || (opt_idx == np_idx && isTRUE(opt_g[opt_idx + 1L] > np_g[np_idx + 1L]))) {
      opt_y <<- np_y; opt_g <<- np_g; opt_idx <<- np_idx
      need_refill <<- TRUE
      if (opt_idx == m && isTRUE(opt_g[opt_idx + 1L] < stop$minf_max)) need_stop <<- TRUE
    }
  }
  calculate_next_points <- function() {
    next_interval <<- q_pop()
    np_x <<- next_coordinate(next_interval)
    if (isTRUE(np_x >= TX[IR[next_interval]]) || isTRUE(np_x <= TX[IL[next_interval]])) need_stop <<- TRUE
    np_y <<- get_image(np_x)
  }
  external_stop <- function() {
    if (.nlopt_stop_time(stop)) {
      ret_code <<- "MAXTIME_REACHED"
      TRUE
    } else if (.nlopt_stop_forced(stop)) {
      ret_code <<- "FORCED_STOP"
      TRUE
    } else FALSE
  }

  failed <- FALSE
  # FirstIteration
  left <- add_trial(0, numeric(n), numeric(m + 1L), -1L)
  right <- add_trial(1, numeric(n), numeric(m + 1L), -1L)
  np_x <- 0.5
  np_y <- get_image(np_x)
  make_trials()
  estimate_optimum()
  tn <- add_trial(np_x, np_y, np_g, np_idx)
  for (pr in list(c(left, tn), c(tn, right))) {
    id <- add_interval(pr[1L], pr[2L])
    DL[id] <- (TX[pr[2L]] - TX[pr[1L]])^(1 / n)
    min_delta <- min(min_delta, DL[id])
    pos <- set_insert(id)
    update_all_h(pos)
  }
  refill_queue()
  calculate_next_points()
  make_trials()
  iterations <- iterations + 2
  repeat {
    # InsertIntervals
    tn <- add_trial(np_x, np_y, np_g, np_idx)
    old <- next_interval
    id <- add_interval(tn, IR[old])
    IR[old] <- tn
    DL[old] <- (TX[IR[old]] - TX[IL[old]])^(1 / n)
    DL[id] <- (TX[IR[id]] - TX[IL[id]])^(1 / n)
    min_delta <- min(min_delta, DL[id])
    min_delta <- min(min_delta, DL[old])
    pos <- set_insert(id)
    if (pos == 0L) { failed <- TRUE; break }
    update_all_h(pos)
    update_all_h(pos - 1L)
    if (!need_refill) {
      RV[id] <- calculate_r(id)
      RV[old] <- calculate_r(old)
      q_push(id)
      q_push(old)
    }
    estimate_optimum()
    if (need_refill || length(hq) < 1L) refill_queue()
    calculate_next_points()
    make_trials()
    need_stop <- need_stop || min_delta < eps || external_stop()
    iterations <- iterations + 1
    if (!(iterations < iters_limit && !need_stop)) break
  }
  if (failed) return(list(ret = "FAILURE", x = lb, minf = Inf))
  if (opt_idx != m) return(list(ret = "FAILURE", x = lb, minf = Inf))
  ret <- if (cnt[1L] >= iters_limit) "MAXEVAL_REACHED" else ret_code
  list(ret = ret, x = opt_y, minf = opt_g[opt_idx + 1L])
}

# the GN_AGS case of nlopt_optimize_
.nlopt_ags_dispatch <- function(o, x, stop) {
  n <- length(x)
  if (!.nlopt_finite_domain(o$lb, o$ub)) return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  if (n > 10L) return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  fc <- o$fc
  if (length(fc) > 10L || any(vapply(fc, function(cc) cc$m, 0L) != 1L))
    return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  fcons <- lapply(fc, function(cc) function(xx) cc$f(xx)[1L])
  fobj <- function(xx) {
    stop$nevals <- stop$nevals + 1L
    o$f(xx)
  }
  r <- .nlopt_ags(fobj, fcons, o$lb, o$ub, stop)
  if (r$ret == "FAILURE") return(list(ret = "FAILURE", x = x, minf = Inf))
  r
}
