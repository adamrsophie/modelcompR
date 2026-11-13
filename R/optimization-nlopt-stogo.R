# StoGO of NLopt 2.11.0 (GD_STOGO and GD_STOGO_RAND, src/algs/stogo:
# stogo.cc, global.cc, local.cc, tools.cc, linalg.cc, by S. Gudmundsson and
# K. Madsen): branch and bound on boxes, each box sampled by a regular
# pattern (2n + 1 points, or the midpoint and 2n random points for
# GD_STOGO_RAND) from which a trust region BFGS local search runs (the INI2
# initial radius of stogo_config.h).
#
# The C++ is compiled with floating point contraction on (and vectorized at
# -O3), which this port reproduces operation by operation: axpy, ger and the
# 'T' gemv are fused; norm2 and dot round the products of the first
# floor(n / 4) * 4 (norm2) or floor(n / 8) * 8 (dot) elements and add them
# in order, then fuse the rest; the 'N' gemv does the same with blocks of 8;
# the regular sampling, the dogleg radius and the lower bound of a box each
# fuse one product. The candidate and garbage sets are std::priority_queue
# of libc++, whose order on ties of the box minimum (all the boxes without
# a stationary point tie at DBL_MAX) decides the search.
#
# Copyright (c) 1998 Stefan Gudmundsson and Kaj Madsen; NLopt modifications
# Copyright (c) 2007-2014 Massachusetts Institute of Technology, released
# under the MIT license (see src/algs/stogo/COPYRIGHT).

# std::priority_queue of libc++ on `h` (ids): push_heap and the Floyd
# pop_heap, with the order key[a] < key[b] (the top has the largest key)
.nlopt_heap_sift_up <- function(h, len, key) {
  if (len > 1L) {
    L <- (len - 2L) %/% 2L
    last <- len - 1L
    if (key[h[L + 1L]] < key[h[last + 1L]]) {
      t <- h[last + 1L]
      repeat {
        h[last + 1L] <- h[L + 1L]
        last <- L
        if (L == 0L) break
        L <- (L - 1L) %/% 2L
        if (!(key[h[L + 1L]] < key[t])) break
      }
      h[last + 1L] <- t
    }
  }
  h
}
.nlopt_heap_push <- function(h, id, key) {
  force(id) # (`id` may extend `key`)
  .nlopt_heap_sift_up(c(h, id), length(h) + 1L, key)
}
# returns list(top, h)
.nlopt_heap_pop <- function(h, key) {
  len <- length(h)
  top <- h[1L]
  if (len > 1L) {
    hole <- 0L; child <- 0L
    lim <- (len - 2L) %/% 2L
    repeat {
      child <- 2L * child + 1L
      if (child + 1L < len && key[h[child + 1L]] < key[h[child + 2L]]) child <- child + 1L
      h[hole + 1L] <- h[child + 1L]
      hole <- child
      if (child > lim) break
    }
    last <- len - 1L
    if (hole == last) {
      h[hole + 1L] <- top
    } else {
      h[hole + 1L] <- h[last + 1L]
      h[last + 1L] <- top
      h <- .nlopt_heap_sift_up(h, hole + 1L, key)
    }
  }
  list(top = top, h = h[-len])
}

# the sum of a * b of linalg.cc: the first floor(n / blk) * blk products
# rounded and added in order (the vectorized loop), the rest fused
.nlopt_stogo_dot <- function(a, b, blk = 8L) {
  n <- length(a)
  k <- (n %/% blk) * blk
  s <- 0
  for (i in seq_len(k)) s <- s + a[i] * b[i]
  for (i in seq_len(n - k) + k) s <- .fma(a[i], b[i], s)
  s
}
.nlopt_stogo_norm2 <- function(a) sqrt(.nlopt_stogo_dot(a, a, 4L))
# gemv(trans, alpha, A, x, beta, y)
.nlopt_stogo_gemv <- function(trans, alpha, A, x, beta, y) {
  n <- length(x)
  s <- numeric(n)
  if (trans == "N") {
    k <- (n %/% 8L) * 8L
    for (j in seq_len(k)) s <- s + (A[, j] * x[j]) * alpha
    for (j in seq_len(n - k) + k) s <- .fma(A[, j] * x[j], alpha, s)
  } else {
    for (j in seq_len(n)) s <- .fma(A[j, ] * x[j], alpha, s)
  }
  .fma(y, beta, s)
}
# ger(alpha, x, y, A): A + alpha x y'
.nlopt_stogo_ger <- function(alpha, x, y, A) {
  n <- length(x)
  .fma(matrix(alpha * x, n, n), matrix(y, n, n, byrow = TRUE), A)
}
.nlopt_stogo_norminf <- function(x) {
  t <- .Machine$double.xmin
  for (v in abs(x)) if (isTRUE(t < v)) t <- v
  t
}

# stogo_minimize with Global::Search: `fg` returns list(val, grad); returns
# list(ret, x, minf)
.nlopt_stogo <- function(fg, lb, ub, stop, nrandom) {
  n <- length(lb)
  DBL_MAX <- .Machine$double.xmax
  eps_cl <- 0.1; rshift <- 0.3; mu <- 1e-4
  rnd_pnts <- nrandom
  det_pnts <- min(2L * n + 1L - nrandom, 2L * n + 1L)
  mac_eps <- .Machine$double.eps
  dom_lb <- lb; dom_ub <- ub
  best_f <- DBL_MAX; best_x <- numeric(n)
  fbound <- DBL_MAX
  in_time <- function() !(.nlopt_stop_evals(stop) || .nlopt_stop_time(stop))
  # MyGlobal::ObjectiveGradient
  obj_grad <- function(x) {
    r <- fg(x)
    if (!any(x < dom_lb | x > dom_ub, na.rm = TRUE) && isTRUE(r$val < best_f)) {
      best_f <<- r$val
      best_x <<- x
    }
    r
  }
  # a box: lb, ub, minf and its trials (rows of X, values F)
  new_box <- function(blb, bub) list(lb = blb, ub = bub, minf = DBL_MAX, X = matrix(0, 0, n), F = numeric(0))
  add_trial <- function(b, x, f) {
    b$X <- rbind(b$X, x, deparse.level = 0)
    b$F <- c(b$F, f)
    if (isTRUE(f < b$minf)) b$minf <- f
    b
  }
  # TBox::OutsideBox: 0 inside the box, 1 outside the box in the domain, 2
  # outside the domain
  outside_box <- function(b, x) {
    ins_box <- TRUE
    for (i in seq_len(n)) {
      if (x[i] < b$lb[i] || x[i] > b$ub[i]) ins_box <- FALSE
      if (x[i] < dom_lb[i] || x[i] > dom_ub[i]) return(2L)
    }
    if (ins_box) 0L else 1L
  }
  inside_box <- function(b, x) !any(x < b$lb | x > b$ub, na.rm = TRUE)
  # TBox::CloseToMin: the first trial within eps_cl, or NULL
  close_to_min <- function(b, v) {
    for (k in seq_along(b$F)) {
      x <- b$X[k, ]
      if (.nlopt_stogo_norm2(v - x) <= eps_cl) return(list(x = x, f = b$F[k]))
    }
    NULL
  }
  shortest_side <- function(b) {
    t <- b$ub[1L] - b$lb[1L]
    for (i in seq_len(n)[-1L]) if (b$ub[i] - b$lb[i] < t) t <- b$ub[i] - b$lb[i]
    t
  }
  longest_side <- function(b) {
    t <- b$ub[1L] - b$lb[1L]; j <- 1L
    for (i in seq_len(n)[-1L]) if (b$ub[i] - b$lb[i] > t) { t <- b$ub[i] - b$lb[i]; j <- i }
    j
  }
  closest_side <- function(b, x) {
    dist <- DBL_MAX
    for (i in seq_len(n)) {
      a1 <- x[i] - b$lb[i]; a2 <- b$ub[i] - x[i]
      t <- if (isTRUE(a2 < a1)) a2 else a1
      if (isTRUE(t < dist)) dist <- t
    }
    dist
  }
  lower_bound <- function(b, maxgrad) {
    lbound <- b$minf
    k <- length(b$F)
    for (i1 in seq_len(k - 1L)) for (i2 in (i1 + 1L):k) {
      est <- 0.5 * .fma(-maxgrad, .nlopt_stogo_norm2(b$X[i1, ] - b$X[i2, ]), b$F[i1] + b$F[i2])
      if (isTRUE(est < lbound)) lbound <- est
    }
    lbound
  }
  split_box <- function(b) {
    B1 <- new_box(b$lb, b$ub); B2 <- new_box(b$lb, b$ub)
    i <- longest_side(b)
    ns <- length(b$F)
    if (ns <= 1L) {
      w <- b$ub[i] - b$lb[i]
      m <- b$lb[i] + w / 2
      B1$ub[i] <- m; B2$lb[i] <- m
    } else {
      center <- numeric(n)
      for (k in seq_len(ns)) center <- center + b$X[k, ]
      center <- (1 / ns) * center
      dispers <- numeric(n)
      for (k in seq_len(ns)) {
        d <- center - b$X[k, ]
        dispers <- dispers + d * d
      }
      dispers <- (1 / ns) * dispers
      i <- 1L; t <- dispers[1L]
      for (k in seq_len(n)[-1L]) if (isTRUE(dispers[k] > t)) { t <- dispers[k]; i <- k }
      B1$ub[i] <- center[i]; B2$lb[i] <- center[i]
    }
    fm1 <- DBL_MAX; fm2 <- DBL_MAX
    for (k in seq_len(ns)) {
      x <- b$X[k, ]; f <- b$F[k]
      if (inside_box(B1, x)) {
        if (isTRUE(f < fm1)) fm1 <- f
        B1 <- add_trial(B1, x, f)
      } else {
        B2 <- add_trial(B2, x, f)
        if (isTRUE(f < fm2)) fm2 <- f
      }
    }
    B1$minf <- fm1; B2$minf <- fm2
    list(B1, B2)
  }

  # local(): the trust region BFGS search from x in the box b; returns
  # list(info, x, f, maxgrad); info 0 LS_Unstable, 1 LS_MaxIter, 2 LS_Old,
  # 3 LS_New, 4 LS_Out, 5 LS_MaxEvalTime
  local_search <- function(x, b, mgr) {
    ctm <- close_to_min(b, x)
    if (!is.null(ctm)) return(list(info = 2L, x = x, f = ctm$f, maxgrad = mgr))
    B <- diag(n); H <- diag(n)
    z <- numeric(n); w <- numeric(n); h_n <- numeric(n)
    r <- obj_grad(x)
    f <- r$val; g <- r$grad
    stop$nevals <- stop$nevals + 1L
    if (!in_time()) return(list(info = 5L))
    delta <- closest_side(b, x) * 0.5
    if (delta < mac_eps) delta <- (1 / 7) * shortest_side(b)
    k <- 0L; good <- FALSE; info <- 3L; outside <- 0L
    maxgrad <- mgr
    while (!good) {
      k <- k + 1L
      if (k > 50L * n) { info <- 1L; break }
      ni <- .nlopt_stogo_norminf(g)
      if (isTRUE(maxgrad < ni)) maxgrad <- ni
      h_sd <- -g
      nrm_sd <- .nlopt_stogo_norm2(h_sd)
      if (isTRUE(nrm_sd < 1e-4)) { good <- TRUE; break }
      h_n <- .nlopt_stogo_gemv("N", -1, H, g, 0, h_n)
      nrm_hn <- .nlopt_stogo_norm2(h_n)
      if (isTRUE(nrm_hn < delta)) {
        h_dl <- h_n
      } else {
        z <- .nlopt_stogo_gemv("N", 1, B, g, 0, z)
        tmp <- .nlopt_stogo_dot(g, z)
        if (isTRUE(tmp == 0)) { info <- 0L; break }
        alpha <- (nrm_sd * nrm_sd) / tmp
        h_sd <- alpha * h_sd
        nrm_sd <- abs(alpha) * nrm_sd
        if (isTRUE(nrm_sd >= delta)) {
          h_dl <- (delta / nrm_sd) * h_sd
        } else {
          d2 <- delta * delta
          s2 <- nrm_sd * nrm_sd
          nom <- d2 - s2
          snrm_hn <- nrm_hn * nrm_hn
          tmp <- .nlopt_stogo_dot(h_n, h_sd)
          den <- tmp - s2 + sqrt(.fma(tmp - d2, tmp - d2, (snrm_hn - d2) * (d2 - s2)))
          if (isTRUE(den == 0)) { info <- 0L; break }
          beta <- nom / den
          h_dl <- .fma(1 - beta, h_sd, beta * h_n)
        }
      }
      nrm_dl <- .nlopt_stogo_norm2(h_dl)
      x_new <- x + h_dl
      it <- outside_box(b, x_new)
      if (it == 1L) {
        outside <- outside + 1L
        if (outside > 1L) break
      } else if (it == 2L) {
        info <- 4L
        break
      } else {
        outside <- 0L
      }
      r <- obj_grad(x_new)
      f_new <- r$val; g_new <- r$grad
      stop$nevals <- stop$nevals + 1L
      if (!in_time()) return(list(info = 5L))
      z <- .nlopt_stogo_gemv("N", 0.5, B, h_dl, 0, z)
      ro <- (f_new - f) / (.nlopt_stogo_dot(g, h_dl) + .nlopt_stogo_dot(h_dl, z))
      if (isTRUE(ro > 0.75)) delta <- delta * 2
      if (isTRUE(ro < 0.25)) delta <- delta / 3
      if (isTRUE(ro > 0)) {
        y <- g_new - g
        alpha <- .nlopt_stogo_dot(y, h_dl)
        if (!isTRUE(alpha <= sqrt(mac_eps) * nrm_dl * .nlopt_stogo_norm2(y))) {
          z <- .nlopt_stogo_gemv("N", 1, B, h_dl, 0, z)
          beta <- -1 / .nlopt_stogo_dot(h_dl, z)
          B <- .nlopt_stogo_ger(1 / alpha, y, y, B)
          B <- .nlopt_stogo_ger(beta, z, z, B)
          z <- .nlopt_stogo_gemv("N", 1, H, y, 0, z)
          w <- .nlopt_stogo_gemv("T", 1, H, y, 0, w)
          beta <- .nlopt_stogo_dot(y, z)
          beta <- (1 + beta / alpha) / alpha
          H <- .nlopt_stogo_ger(beta, h_dl, h_dl, H)
          H <- .nlopt_stogo_ger(-1 / alpha, z, h_dl, H)
          H <- .nlopt_stogo_ger(-1 / alpha, h_dl, w, H)
        }
        if (isTRUE(nrm_dl < .nlopt_stogo_norm2(x) * 1e-4)) good <- TRUE
        ctm <- close_to_min(b, x_new)
        if (!is.null(ctm)) {
          x_new <- ctm$x; f_new <- ctm$f
          info <- 2L
          good <- TRUE
        }
        x <- x_new; g <- g_new; f <- f_new
      }
    }
    if (outside_box(b, x) != 0L) { info <- 4L; f <- DBL_MAX }
    if (info == 0L) stop("StoGO: the local search became unstable (NLopt exits the process here)")
    list(info = if (outside > 0L) 4L else info, x = x, f = f, maxgrad = maxgrad)
  }

  # Global::NewtonTest: samples the box b and runs the local searches;
  # returns list(box, maxgrad, nout)
  newton_test <- function(b) {
    nout <- 0L
    maxgrad <- 0
    samples <- list()
    # FillRandom
    for (i in seq_len(rnd_pnts)) {
      x <- numeric(n)
      for (dir in seq_len(n)) x[dir] <- .nlopt_urand(b$lb[dir], b$ub[dir])
      samples[[length(samples) + 1L]] <- x
    }
    # FillRegular
    if (det_pnts > 0L) {
      m <- abs(b$ub - b$lb) / 2 + b$lb
      x <- m
      flag <- 1L; dir <- 1L
      for (i in seq_len(det_pnts - 1L)) {
        w <- b$ub[dir] - b$lb[dir]
        x[dir] <- .fma(flag * rshift, w, m[dir])
        samples[[length(samples) + 1L]] <- x
        flag <- -flag
        if (flag == 1L && dir <= n) {
          x[dir] <- m[dir]
          dir <- dir + 1L
        }
      }
      samples[[length(samples) + 1L]] <- m
    }
    while (length(samples)) {
      x <- samples[[length(samples)]]
      samples[[length(samples)]] <- NULL
      ls <- local_search(x, b, maxgrad)
      info <- ls$info
      if (info != 5L) maxgrad <- ls$maxgrad
      if (info == 4L) {
        nout <- nout + 1L
      } else if (info == 3L) {
        b <- add_trial(b, ls$x, ls$f)
        if (isTRUE(ls$f <= fbound + mu) && isTRUE(ls$f <= b$minf + mu) && isTRUE(ls$f < stop$minf_max)) break
      }
      if (!in_time() || info == 5L) break
    }
    list(box = b, maxgrad = maxgrad, nout = nout)
  }

  # the candidate and garbage sets: heaps of box ids, keyed by -minf
  boxes <- list(); key <- numeric(0)
  cand <- integer(0); garbage <- integer(0)
  store <- function(b) {
    boxes[[length(boxes) + 1L]] <<- b
    key[length(boxes)] <<- -b$minf
    length(boxes)
  }
  cand <- .nlopt_heap_push(cand, store(new_box(lb, ub)), key)
  done <- FALSE
  while (!done) {
    while (length(cand)) {
      p <- .nlopt_heap_pop(cand, key)
      cand <- p$h
      b <- boxes[[p$top]]
      boxes[p$top] <- list(NULL)
      # ReduceOrSubdivide
      nt <- newton_test(b)
      b <- nt$box
      ns <- length(b$F)
      if (ns == 0L || (ns == 1L && nt$nout == 0L) ||
          (ns > 1L && isTRUE(lower_bound(b, nt$maxgrad) > fbound))) {
        garbage <- .nlopt_heap_push(garbage, store(b), key)
      } else {
        bb <- split_box(b)
        cand <- .nlopt_heap_push(cand, store(bb[[1L]]), key)
        cand <- .nlopt_heap_push(cand, store(bb[[2L]]), key)
      }
      if (isTRUE(b$minf < fbound)) fbound <- b$minf
      if (best_f < DBL_MAX && isTRUE(best_f < stop$minf_max)) { done <- TRUE; break }
      if (!in_time()) { done <- TRUE; break }
    }
    if (in_time()) {
      while (length(garbage)) {
        p <- .nlopt_heap_pop(garbage, key)
        garbage <- p$h
        bb <- split_box(boxes[[p$top]])
        boxes[p$top] <- list(NULL)
        cand <- .nlopt_heap_push(cand, store(bb[[1L]]), key)
        cand <- .nlopt_heap_push(cand, store(bb[[2L]]), key)
      }
    }
  }
  if (!(best_f < DBL_MAX)) return(list(ret = "FAILURE", x = lb, minf = Inf))
  list(ret = "SUCCESS", x = best_x, minf = best_f)
}

# the GD_STOGO cases of nlopt_optimize_
.nlopt_stogo_dispatch <- function(o, x, stop) {
  n <- length(x)
  if (!.nlopt_finite_domain(o$lb, o$ub)) return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  pop <- as.integer(o$population %||% 0L)
  nrandom <- if (o$algorithm == "GD_STOGO") 0L else if (pop > 0L) pop else 2L * n
  r <- .nlopt_stogo(o$fg, o$lb, o$ub, stop, nrandom)
  if (r$ret == "FAILURE") return(list(ret = "FAILURE", x = x, minf = Inf))
  r
}
