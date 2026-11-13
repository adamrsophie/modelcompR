# The Nelder Mead simplex (LN_NELDERMEAD, src/algs/neldermead/nldrmd.c) and
# Rowan's subplex method (LN_SBPLX, sbplx.c) of NLopt 2.11.0, ported to R.
#
# Copyright (c) 2007-2014 Massachusetts Institute of Technology, released
# under the MIT license (see R/optimization-nlopt-util.R).
#
# The red black tree of nldrmd.c, which keeps the simplex sorted by f(x)
# with ties broken by the position of the vertex in memory, is a vector of
# vertex numbers ordered by (f, vertex number).

.nlopt_nm_close <- function(a, b) abs(a - b) <= 1e-13 * (abs(a) + abs(b))

# reflectpt: xnew = c + scale (c - xold) pinned to the bounds; ok is FALSE
# when xnew coincides with c or xold
.nlopt_nm_reflect <- function(c, scale, xold, lb, ub) {
  xnew <- c + scale * (c - xold)
  xnew <- ifelse(xnew < lb, lb, xnew)
  xnew <- ifelse(xnew > ub, ub, xnew)
  equalc <- all(.nlopt_nm_close(xnew, c))
  equalold <- all(.nlopt_nm_close(xnew, xold))
  list(x = xnew, ok = !(equalc || equalold))
}

# nldrmd_minimize_: `minf` holds f(x) on entry; psi > 0 replaces the x and
# f tolerances by a reduction of the simplex diameter (for subplex).
# Returns list(ret, x, minf, fdiff).
.nlopt_nldrmd_minimize_ <- function(n, f, lb, ub, x, minf, xstep, stop, psi) {
  pts <- matrix(0, n + 1L, n)         # row i: vertex i
  fv <- numeric(n + 1L)
  fdiff <- Inf
  init_diam <- 0
  done <- function(ret) list(ret = ret, x = x, minf = minf, fdiff = fdiff)
  halt <- function(ret) stop(structure(class = c("nlopt_nm_done", "condition"), list(message = ret, call = NULL)))
  # CHECK_EVAL
  check <- function(xc, fc) {
    stop$nevals <- stop$nevals + 1L
    if (.nlopt_stop_forced(stop)) halt("FORCED_STOP")
    if (fc <= minf) {
      minf <<- fc; x <<- xc
      if (minf < stop$minf_max) halt("STOPVAL_REACHED")
    }
    if (.nlopt_stop_evals(stop)) halt("MAXEVAL_REACHED")
    if (.nlopt_stop_time(stop)) halt("MAXTIME_REACHED")
  }
  run <- function() {
    pts[1L, ] <<- x
    fv[1L] <<- minf
    if (minf < stop$minf_max) return(done("STOPVAL_REACHED"))
    for (i in seq_len(n)) {
      # each vertex starts from the current best x, which check() may move
      x0 <- x
      pt <- x0
      pt[i] <- pt[i] + xstep[i]
      if (pt[i] > ub[i]) {
        pt[i] <- if (ub[i] - x0[i] > abs(xstep[i]) * 0.1) ub[i] else x0[i] - abs(xstep[i])
      }
      if (pt[i] < lb[i]) {
        if (x0[i] - lb[i] > abs(xstep[i]) * 0.1) pt[i] <- lb[i]
        else {
          pt[i] <- x0[i] + abs(xstep[i])
          if (pt[i] > ub[i]) pt[i] <- 0.5 * ((if (ub[i] - x0[i] > x0[i] - lb[i]) ub[i] else lb[i]) + x0[i])
        }
      }
      if (.nlopt_nm_close(pt[i], x0[i])) {
        stop$msg <- sprintf("starting step size led to simplex that was too small in dimension %d: %g is too close to x[%d]=%g",
                            i - 1L, pt[i], i - 1L, x0[i])
        return(done("FAILURE"))
      }
      pts[i + 1L, ] <<- pt
      fv[i + 1L] <<- f(pt)
      check(pt, fv[i + 1L])
    }
    repeat { # restart
      repeat {
        ord <- order(fv, seq_len(n + 1L))
        lo <- ord[1L]; hi <- ord[n + 1L]
        fl <- fv[lo]; xl <- pts[lo, ]
        fh <- fv[hi]; xh <- pts[hi, ]
        fdiff <<- fh - fl
        if (init_diam == 0) init_diam <<- sum(abs(xl - xh))
        if (psi <= 0 && .nlopt_stop_ftol(stop, fl, fh)) return(done("FTOL_REACHED"))
        # centroid of the vertices other than the highest
        c <- numeric(n)
        for (i in seq_len(n + 1L)) if (i != hi) c <- c + pts[i, ]
        c <- c * (1 / n)
        # x convergence: xcur = centroid plus the largest distance from it
        xcur <- numeric(n)
        for (i in seq_len(n + 1L)) {
          dx <- abs(pts[i, ] - c)
          xcur <- ifelse(dx > xcur, dx, xcur)
        }
        xcur <- xcur + c
        if (psi > 0) {
          if (sum(abs(xl - xh)) < psi * init_diam) return(done("XTOL_REACHED"))
        } else if (.nlopt_stop_x(stop, c, xcur)) return(done("XTOL_REACHED"))
        # reflection
        r <- .nlopt_nm_reflect(c, 1, xh, lb, ub)
        if (!r$ok) return(done("XTOL_REACHED"))
        xcur <- r$x
        fr <- f(xcur)
        check(xcur, fr)
        shrink <- FALSE
        if (fr < fl) { # new best point, expand
          r <- .nlopt_nm_reflect(c, 2, xh, lb, ub)
          if (!r$ok) return(done("XTOL_REACHED"))
          xh <- r$x
          fh <- f(xh)
          check(xh, fh)
          if (fh >= fr) { fh <- fr; xh <- xcur }
        } else if (fr < fv[ord[n]]) { # accept the new point
          xh <- xcur; fh <- fr
        } else { # new worst point, contract
          r <- .nlopt_nm_reflect(c, if (fh <= fr) -0.5 else 0.5, xh, lb, ub)
          if (!r$ok) return(done("XTOL_REACHED"))
          xcur <- r$x
          fc <- f(xcur)
          check(xcur, fc)
          if (fc < fr && fc < fh) { xh <- xcur; fh <- fc }
          else shrink <- TRUE
        }
        if (shrink) break
        pts[hi, ] <<- xh
        fv[hi] <<- fh
      }
      # failed contraction: shrink the simplex toward the lowest vertex
      for (i in seq_len(n + 1L)) {
        if (i != lo) {
          r <- .nlopt_nm_reflect(xl, -0.5, pts[i, ], lb, ub)
          if (!r$ok) return(done("XTOL_REACHED"))
          pts[i, ] <<- r$x
          fv[i] <<- f(r$x)
          check(r$x, fv[i])
        }
      }
    }
  }
  tryCatch(run(), nlopt_nm_done = function(e) done(conditionMessage(e)))
}

# nldrmd_minimize
.nlopt_nldrmd_minimize <- function(f, lb, ub, x, xstep, stop) {
  n <- length(x)
  minf <- f(x)
  stop$nevals <- stop$nevals + 1L
  if (.nlopt_stop_forced(stop)) return(list(ret = "FORCED_STOP", x = x, minf = minf))
  if (minf < stop$minf_max) return(list(ret = "STOPVAL_REACHED", x = x, minf = minf))
  if (.nlopt_stop_evals(stop)) return(list(ret = "MAXEVAL_REACHED", x = x, minf = minf))
  if (.nlopt_stop_time(stop)) return(list(ret = "MAXTIME_REACHED", x = x, minf = minf))
  .nlopt_nldrmd_minimize_(n, f, lb, ub, x, minf, xstep, stop, 0)[c("ret", "x", "minf")]
}

# sbplx_minimize
.nlopt_sbplx_minimize <- function(f, lb, ub, x, xstep0, stop) {
  psi <- 0.25; omega <- 0.1; nsmin <- 2L; nsmax <- 5L
  n <- length(x)
  minf <- f(x)
  stop$nevals <- stop$nevals + 1L
  if (.nlopt_stop_forced(stop)) return(list(ret = "FORCED_STOP", x = x, minf = minf))
  if (minf < stop$minf_max) return(list(ret = "STOPVAL_REACHED", x = x, minf = minf))
  if (.nlopt_stop_evals(stop)) return(list(ret = "MAXEVAL_REACHED", x = x, minf = minf))
  if (.nlopt_stop_time(stop)) return(list(ret = "MAXTIME_REACHED", x = x, minf = minf))
  xstep <- xstep0
  dx <- numeric(n)
  # Nelder Mead on the subspace p[is + 1:ns] of x
  sub <- function(is, ns) {
    idx <- p[is + seq_len(ns)]
    sf <- function(xs) { xx <- x; xx[idx] <- xs; x <<- xx; f(xx) }
    r <- .nlopt_nldrmd_minimize_(ns, sf, lb[idx], ub[idx], x[idx], minf, xstep[idx], stop, psi)
    x[idx] <<- r$x
    minf <<- r$minf
    r
  }
  repeat {
    normi <- 0
    nsubs <- 0L
    fdiff_max <- 0
    xprev <- x
    # indices sorted by decreasing |dx|
    adx <- abs(dx)
    p <- .nlopt_qsort(seq_len(n), function(i, j) if (adx[i] > adx[j]) -1 else if (adx[i] < adx[j]) 1 else 0)
    normdx <- sum(abs(dx))
    i <- 0L
    while (i + nsmin < n) {
      ns_goodness <- -Inf
      norm <- normi
      nk <- if (i + nsmax > n) n else i + nsmax
      k <- i
      while (k < i + nsmin - 1L) { norm <- norm + abs(dx[p[k + 1L]]); k <- k + 1L }
      ns <- nsmin
      k <- i + nsmin - 1L
      while (k < nk) {
        norm <- norm + abs(dx[p[k + 1L]])
        if (!(n - (k + 1L) < nsmin)) {
          goodness <- if (k + 1L < n) norm / (k + 1L) - (normdx - norm) / (n - (k + 1L)) else normdx / n
          if (goodness > ns_goodness) { ns_goodness <- goodness; ns <- (k + 1L) - i }
        }
        k <- k + 1L
      }
      for (k in i + seq_len(ns)) normi <- normi + abs(dx[p[k]])
      nsubs <- nsubs + 1L
      r <- sub(i, ns)
      if (r$fdiff > fdiff_max) fdiff_max <- r$fdiff
      if (r$ret == "FAILURE") return(list(ret = "XTOL_REACHED", x = x, minf = minf))
      if (r$ret != "XTOL_REACHED") return(list(ret = r$ret, x = x, minf = minf))
      i <- i + ns
    }
    # the last subspace
    ns <- n - i
    nsubs <- nsubs + 1L
    r <- sub(i, ns)
    if (r$fdiff > fdiff_max) fdiff_max <- r$fdiff
    if (r$ret == "FAILURE") return(list(ret = "XTOL_REACHED", x = x, minf = minf))
    if (r$ret != "XTOL_REACHED") return(list(ret = r$ret, x = x, minf = minf))
    # termination tests
    if (.nlopt_stop_ftol(stop, minf, minf + fdiff_max)) return(list(ret = "FTOL_REACHED", x = x, minf = minf))
    if (.nlopt_stop_x(stop, x, xprev)) {
      xa <- if (is.null(stop$xtol_abs)) rep(0, n) else stop$xtol_abs
      if (!any(abs(xstep) * psi > xa & abs(xstep) * psi > stop$xtol_rel * abs(x)))
        return(list(ret = "XTOL_REACHED", x = x, minf = minf))
    }
    dx <- x - xprev
    scale <- if (nsubs == 1L) psi else {
      sc <- sum(abs(dx)) / sum(abs(xstep))
      min(max(sc, omega), 1 / omega)
    }
    xstep <- ifelse(dx == 0, -(xstep * scale), abs(xstep * scale) * ifelse(dx < 0, -1, 1))
  }
}
