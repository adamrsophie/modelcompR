# The DIRECT and DIRECT-L algorithms of NLopt 2.11.0 (src/algs/cdirect/
# cdirect.c, Steven G. Johnson's re-implementation of Jones' DIRECT and
# Gablonsky's DIRECT-L), ported to R: GN_DIRECT, GN_DIRECT_L,
# GN_DIRECT_L_RAND and their unscaled variants (GN_DIRECT_NOSCAL,
# GN_DIRECT_L_NOSCAL, GN_DIRECT_L_RAND_NOSCAL).
#
# The red black tree of cdirect.c (util/redblack.c) keeps the hyper
# rectangles sorted by (diameter, f, age) with unique ages, so its in order
# traversal is that of an ordered set: here an integer vector `ord` of
# rectangle ids, kept sorted. The convex hull holds rectangle ids and reads
# their current diameter and f, as cdirect.c reads them through pointers to
# rectangles that the divisions of the same pass modify.
#
# Copyright (c) 2007-2014 Massachusetts Institute of Technology, released
# under the MIT license (see R/optimization-nlopt-util.R).

# (float) x: a double rounded to single precision
.nlopt_float <- function(x) readBin(writeBin(as.numeric(x), raw(), size = 4L), "double", size = 4L)

.nlopt_cd_third <- 0x1.5555555555555p-2

# cdirect_unscaled: `f` takes a point of [lb, ub]; returns list(ret, x,
# minf). `xtol_abs` is the absolute tolerance on the widths (NULL for none).
.nlopt_cdirect_unscaled <- function(f, lb, ub, x, stop, magic_eps, which_alg, xtol_abs = stop$xtol_abs) {
  n <- length(x)
  which_diam <- which_alg %% 3L
  which_div <- (which_alg %/% 3L) %% 3L
  which_opt <- (which_alg %/% 9L) %% 3L
  third <- .nlopt_cd_third
  # the state of the search, updated by the closures below with <<-
  minf <- Inf; xmin <- x; age <- 0
  cap <- 256L
  C <- matrix(0, cap, n); W <- matrix(0, cap, n)
  D <- numeric(cap); Fv <- numeric(cap); A <- numeric(cap)
  nr <- 0L
  ord <- integer(0)
  halt <- function(ret) stop(structure(class = c("nlopt_cd_halt", "condition"), list(message = ret, call = NULL)))

  diameter <- function(w) {
    if (which_diam == 0L) {
      s <- 0
      for (i in seq_len(n)) s <- s + w[i] * w[i]
      .nlopt_float(sqrt(s) * 0.5)
    } else {
      maxw <- 0
      for (i in seq_len(n)) if (w[i] > maxw) maxw <- w[i]
      .nlopt_float(maxw * 0.5)
    }
  }
  function_eval <- function(xx) {
    v <- f(xx)
    if (isTRUE(v < minf)) { minf <<- v; xmin <<- xx }
    stop$nevals <- stop$nevals + 1L
    v
  }
  FUNCTION_EVAL <- function(xx) {
    v <- function_eval(xx)
    if (.nlopt_stop_forced(stop)) halt("FORCED_STOP")
    if (minf < stop$minf_max) halt("STOPVAL_REACHED")
    if (.nlopt_stop_evals(stop)) halt("MAXEVAL_REACHED")
    if (.nlopt_stop_time(stop)) halt("MAXTIME_REACHED")
    v
  }
  new_rect <- function(cc, ww, d, fval, age) {
    if (nr == nrow(C)) {
      k <- nrow(C)
      C <<- rbind(C, matrix(0, k, n)); W <<- rbind(W, matrix(0, k, n))
      D <<- c(D, numeric(k)); Fv <<- c(Fv, numeric(k)); A <<- c(A, numeric(k))
    }
    id <- nr <<- nr + 1L
    C[id, ] <<- cc; W[id, ] <<- ww; D[id] <<- d; Fv[id] <<- fval; A[id] <<- age
    id
  }
  # the ordered set: insert by the (d, f, age) order, resort after a key
  # change
  tree_insert <- function(id) {
    o <- ord
    d <- D[id]; fv <- Fv[id]; a <- A[id]
    Do <- D[o]; Fo <- Fv[o]
    less <- Do < d | (Do == d & (Fo < fv | (Fo == fv & A[o] < a)))
    ord <<- append(o, id, after = sum(less))
  }
  tree_resort <- function(id) {
    ord <<- ord[ord != id]
    tree_insert(id)
  }
  # the first position whose diameter is >= k0 (nlopt_rb_tree_find_gt of
  # the key (k0, -HUGE_VAL))
  find_gt <- function(Dv, k0) sum(Dv < k0) + 1L

  divide_rect <- function(id) {
    cc <- C[id, ]; w <- W[id, ]
    wmax <- w[1L]; imax <- 1L
    if (n > 1L) for (i in 2L:n) if (w[i] > wmax) { imax <- i; wmax <- w[i] }
    longest <- wmax - w <= wmax * 5e-2
    nlongest <- sum(longest)
    if (which_div == 1L || (which_div == 0L && nlongest == n)) {
      fv <- numeric(2L * n)
      for (i in seq_len(n)) {
        if (longest[i]) {
          csave <- cc[i]
          cc[i] <- csave - w[i] * third
          fv[2L * i - 1L] <- FUNCTION_EVAL(cc)
          cc[i] <- csave + w[i] * third
          fv[2L * i] <- FUNCTION_EVAL(cc)
          cc[i] <- csave
        } else {
          fv[2L * i - 1L] <- fv[2L * i] <- Inf
        }
      }
      fmin <- vapply(seq_len(n), function(i) {
        a <- fv[2L * i - 1L]; b <- fv[2L * i]
        if (a < b) a else b
      }, 0)
      isort <- .nlopt_qsort(seq_len(n), function(a, b) if (fmin[a] < fmin[b]) -1 else if (fmin[a] > fmin[b]) 1 else 0)
      for (i in seq_len(nlongest)) {
        j <- isort[i]
        w[j] <- w[j] * third
        W[id, ] <<- w
        D[id] <<- diameter(w)
        A[id] <<- age; age <<- age + 1
        tree_resort(id)
        for (k in 0:1) {
          fk <- fv[2L * j - 1L + k]
          if (!is.finite(fk)) return("FAILURE")
          cn <- cc
          cn[j] <- cn[j] + w[j] * (2 * k - 1)
          nid <- new_rect(cn, w, D[id], fk, age)
          age <<- age + 1
          tree_insert(nid)
        }
      }
    } else {
      if (nlongest > 1L && which_div == 2L) {
        i <- .nlopt_iurand(nlongest)
        for (k in seq_len(n)) {
          if (longest[k]) {
            if (i == 0L) { i <- k; break }
            i <- i - 1L
          }
        }
      } else i <- imax
      w[i] <- w[i] * third
      W[id, ] <<- w
      D[id] <<- diameter(w)
      A[id] <<- age; age <<- age + 1
      tree_resort(id)
      for (k in 0:1) {
        cn <- cc
        cn[i] <- cn[i] + w[i] * (2 * k - 1)
        fk <- FUNCTION_EVAL(cn)
        nid <- new_rect(cn, w, D[id], fk, age)
        age <<- age + 1
        if (!is.finite(fk)) return("FAILURE")
        tree_insert(nid)
      }
    }
    "SUCCESS"
  }

  convex_hull <- function(allow_dups) {
    o <- ord
    N <- length(o)
    if (N == 0L) return(integer(0))
    Dv <- D[o]; Fo <- Fv[o]
    D <- D; Fv <- Fv
    hull <- integer(N + 4L); nhull <- 0L
    pos <- 1L
    xmin <- Dv[1L]; yminmin <- Fo[1L]; xmax <- Dv[N]
    if (allow_dups) {
      repeat {
        nhull <- nhull + 1L; hull[nhull] <- o[pos]
        pos <- pos + 1L
        if (!(pos <= N && Dv[pos] == xmin && Fo[pos] == yminmin)) break
      }
    } else {
      nhull <- 1L; hull[1L] <- o[1L]
    }
    if (xmin == xmax) return(hull[seq_len(nhull)])
    nmax <- find_gt(Dv, xmax * (1 - 1e-13))
    ymaxmin <- Fo[nmax]
    minslope <- (ymaxmin - yminmin) / (xmax - xmin)
    pos <- find_gt(Dv, xmin * (1 + 1e-13))
    while (pos != nmax) {
      k0 <- Dv[pos]; k1 <- Fo[pos]
      if (isTRUE(k1 > yminmin + (k0 - xmin) * minslope)) { pos <- pos + 1L; next }
      if (nhull > 0L && k0 == D[hull[nhull]]) {
        if (k1 > Fv[hull[nhull]]) {
          pos <- find_gt(Dv, k0 * (1 + 1e-13))
          next
        } else {
          if (allow_dups) { nhull <- nhull + 1L; hull[nhull] <- o[pos] }
          pos <- pos + 1L
          next
        }
      }
      while (nhull > 1L) {
        t1 <- hull[nhull]
        it2 <- nhull - 1L
        repeat {
          t2 <- hull[it2]
          it2 <- it2 - 1L
          if (!(it2 >= 1L && D[t2] == D[t1] && Fv[t2] == Fv[t1])) break
        }
        if (it2 < 1L) break
        if (isTRUE((D[t1] - D[t2]) * (k1 - Fv[t2]) - (Fv[t1] - Fv[t2]) * (k0 - D[t2]) >= 0)) break
        nhull <- nhull - 1L
      }
      nhull <- nhull + 1L; hull[nhull] <- o[pos]
      pos <- pos + 1L
    }
    if (allow_dups) {
      repeat {
        nhull <- nhull + 1L; hull[nhull] <- o[nmax]
        nmax <- nmax + 1L
        if (!(nmax <= N && Dv[nmax] == xmax && Fo[nmax] == ymaxmin)) break
      }
    } else {
      nhull <- nhull + 1L; hull[nhull] <- o[nmax]
    }
    hull[seq_len(nhull)]
  }

  xtol_rel <- stop$xtol_rel
  small <- function(w) {
    for (i in seq_len(n)) {
      if (w[i] > (if (is.null(xtol_abs)) 0 else xtol_abs[i]) && w[i] > (ub[i] - lb[i]) * xtol_rel) return(FALSE)
    }
    TRUE
  }

  divide_good_rects <- function() {
    hull <- convex_hull(which_opt != 1L)
    nhull <- length(hull)
    meps <- magic_eps
    xtol_reached <- TRUE; divided_some <- FALSE
    repeat {
      i <- 1L
      while (i <= nhull) {
        K1 <- -Inf; K2 <- -Inf
        h <- hull[i]
        im <- i - 1L
        while (im >= 1L && D[hull[im]] == D[h]) im <- im - 1L
        ip <- i + 1L
        while (ip <= nhull && D[hull[ip]] == D[h]) ip <- ip + 1L
        if (im >= 1L) K1 <- (Fv[h] - Fv[hull[im]]) / (D[h] - D[hull[im]])
        if (ip <= nhull) K2 <- (Fv[h] - Fv[hull[ip]]) / (D[h] - D[hull[ip]])
        K <- if (isTRUE(K1 > K2)) K1 else K2
        # (a comparison with NaN, from an infinite minf, is false as in C)
        if (isTRUE(Fv[h] - K * D[h] <= minf - meps * abs(minf)) || ip == nhull + 1L) {
          ret <- divide_rect(h)
          divided_some <- TRUE
          if (ret != "SUCCESS") return(ret)
          xtol_reached <- xtol_reached && small(W[h, ])
        }
        if (which_opt == 1L) i <- ip - 1L
        else if (which_opt == 2L) i <- i + .nlopt_iurand(ip - i)
        i <- i + 1L
      }
      if (!divided_some) {
        if (meps != 0) { meps <- 0; next }
        o <- ord
        pos <- length(o)
        wmax <- D[o[pos]]
        while (pos > 1L && D[o[pos - 1L]] == wmax) pos <- pos - 1L
        return(divide_rect(o[pos]))
      }
      return(if (xtol_reached) "XTOL_REACHED" else "SUCCESS")
    }
  }

  cc <- 0.5 * (lb + ub)
  ww <- ub - lb
  ret <- tryCatch({
    d0 <- diameter(ww)
    f0 <- function_eval(cc)
    id <- new_rect(cc, ww, d0, f0, age)
    age <- age + 1
    tree_insert(id)
    ret <- divide_rect(id)
    if (ret == "SUCCESS") {
      repeat {
        minf0 <- minf
        ret <- divide_good_rects()
        if (ret != "SUCCESS") break
        if (isTRUE(minf < minf0) && .nlopt_stop_f(stop, minf, minf0)) { ret <- "FTOL_REACHED"; break }
      }
    }
    ret
  }, nlopt_cd_halt = function(e) conditionMessage(e))
  list(ret = ret, x = xmin, minf = minf)
}

# cdirect: DIRECT on the unit hypercube, the point rescaled to [lb, ub] for
# each evaluation
.nlopt_cdirect <- function(f, lb, ub, x, stop, magic_eps, which_alg) {
  n <- length(x)
  wd <- ub - lb
  xu <- (x - lb) / wd
  xa <- if (!is.null(stop$xtol_abs)) stop$xtol_abs / wd
  fu <- function(z) f(lb + z * wd)
  res <- .nlopt_cdirect_unscaled(fu, rep(0, n), rep(1, n), xu, stop, magic_eps, which_alg, xtol_abs = xa)
  res$x <- lb + res$x * wd
  res
}

# the which_alg code of the cdirect cases of nlopt_optimize_. The
# _NOSCAL cases test the algorithm against NLOPT_GN_DIRECT and
# NLOPT_GN_DIRECT_L_RAND (not their _NOSCAL names), so all three unscaled
# variants run DIRECT-L (which_alg 13), as in NLopt.
.nlopt_cdirect_which <- function(alg) {
  plain <- alg == "GN_DIRECT"
  rnd <- alg == "GN_DIRECT_L_RAND"
  as.integer((!plain) + 3L * (if (rnd) 2L else !plain) + 9L * (if (rnd) 1L else !plain))
}

.nlopt_cdirect_dispatch <- function(o, x, stop) {
  if (!.nlopt_finite_domain(o$lb, o$ub)) return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  run <- if (grepl("_NOSCAL$", o$algorithm)) .nlopt_cdirect_unscaled else .nlopt_cdirect
  run(o$f, o$lb, o$ub, x, stop, .nlopt_get_param(o, "magic_eps", 0), .nlopt_cdirect_which(o$algorithm))
}
