# BOBYQA of PRIMA (libprima 0.7.1, fortran/bobyqa), the solver that
# OptimizationPRIMA.jl calls through PRIMA.jl: Powell's method for bound
# constrained minimization without derivatives. Translated from the Fortran
# reference with the arithmetic of the libprima build: every a * b + c of the
# source is one fused multiply add (.pw_fma), inner and matrix products are
# fused chains, and MAX and MIN ignore a NaN operand. The common routines of
# powalg (calvlag, calden, updateh, shiftbase, hess_mul, quadinc) are those of
# R/optimization-prima-powell.R.

# Fortran SIGN(a, b) elementwise (the sign bit of b, so that -0 counts as
# negative)
.bq_copysign <- function(a, b) ifelse(b < 0 | (b == 0 & 1 / b < 0), -abs(a), abs(a))

# elementwise MAX and MIN that ignore a NaN operand (fmaxnm and fminnm)
.bq_pmax <- function(a, b) {
  m <- if (length(a) && length(b)) max(length(a), length(b)) else 0L
  a <- rep_len(a, m)
  b <- rep_len(b, m)
  r <- pmax(a, b)
  r[is.nan(a)] <- b[is.nan(a)]
  r[is.nan(b)] <- a[is.nan(b)]
  r
}
.bq_pmin <- function(a, b) {
  m <- if (length(a) && length(b)) max(length(a), length(b)) else 0L
  a <- rep_len(a, m)
  b <- rep_len(b, m)
  r <- pmin(a, b)
  r[is.nan(a)] <- b[is.nan(a)]
  r[is.nan(b)] <- a[is.nan(b)]
  r
}
# first index of the smallest element among those with mask TRUE (minloc with
# a mask), 0 when the mask is empty
.bq_minloc_mask <- function(v, mask) {
  k <- which(mask & !is.na(v))
  if (!length(k)) return(if (any(mask)) which(mask)[1L] else 0L)
  k[which.min(v[k])]
}

# xinbd: XBASE + STEP inside the bounds, at the bound when STEP reaches SL or SU
.bq_xinbd <- function(xbase, step, xl, xu, sl, su) {
  s <- .bq_pmax(sl, .bq_pmin(su, step))
  x <- .bq_pmax(xl, .bq_pmin(xu, xbase + s))
  lo <- which(s <= sl)
  x[lo] <- xl[lo]
  hi <- which(s >= su)
  x[hi] <- xu[hi]
  x
}

# trrad of BOBYQA
.bq_trrad <- function(delta_in, dnorm, eta1, eta2, gamma1, gamma2, ratio) {
  if (ratio <= eta1) .pw_fminnm(gamma1 * delta_in, dnorm)
  else if (ratio <= eta2) .pw_fmaxnm(gamma1 * delta_in, dnorm)
  else .pw_fmaxnm(gamma1 * delta_in, gamma2 * dnorm)
}

# A + alpha (x y' + y x') with the lower triangle computed and mirrored
# (r2update of linalg)
.bq_r2sym <- function(A, alpha, x, y) {
  n <- length(x)
  ax <- rep(alpha * x, times = n)
  ay <- rep(alpha * y, times = n)
  B <- matrix(.pw_fma(ay, rep(x, each = n), .pw_fma(ax, rep(y, each = n), A)), n, n)
  B[upper.tri(B)] <- t(B)[upper.tri(B)]
  B
}

# initxf: the initial interpolation set; sl and su are the bounds relative
# to xbase, at least rhobeg away from xbase or zero
.bq_initxf <- function(fun, iprint, maxfun, ftarget, rhobeg, xl, xu, x0, npt) {
  n <- length(x0)
  sl <- xl - x0
  su <- xu - x0
  neg <- sl < 0
  sl[neg] <- .bq_pmin(sl[neg], -rhobeg)
  x0[!neg] <- xl[!neg]
  sl[!neg] <- 0
  su[!neg] <- xu[!neg] - xl[!neg]
  pos <- su > 0
  su[pos] <- .bq_pmax(su[pos], rhobeg)
  x0[!pos] <- xu[!pos]
  sl[!pos] <- xl[!pos] - xu[!pos]
  su[!pos] <- 0
  xbase <- x0
  evaluated <- logical(npt)
  fval <- rep(.pr_realmax, npt)
  xpt <- matrix(0, n, npt)
  for (k in seq_len(n)) xpt[k, k + 1L] <- if (su[k] <= 0) -rhobeg else rhobeg
  for (k in seq_len(min(npt - n - 1L, n))) {
    xpt[k, k + n + 1L] <- -rhobeg
    if (sl[k] >= 0) xpt[k, k + n + 1L] <- .pw_fminnm(2 * rhobeg, su[k])
    if (su[k] <= 0) xpt[k, k + n + 1L] <- .pw_fmaxnm(-2 * rhobeg, sl[k])
  }
  info <- 0L
  evalk <- function(k) {
    x <- .bq_xinbd(xbase, xpt[, k], xl, xu, sl, su)
    f <- .pw_evaluatef(fun, x)
    .pw_fmsg(iprint, "Initialization", k, rhobeg, f, x)
    evaluated[k] <<- TRUE
    fval[k] <<- f
    .pw_checkexit(maxfun, k, f, ftarget, x)
  }
  for (k in seq_len(min(npt, 2L * n + 1L))) {
    info <- evalk(k)
    if (info != 0L) break
  }
  for (k in seq_len(min(npt - n, n + 1L))[-1L]) {
    if (isTRUE(xpt[k - 1L, k] * xpt[k - 1L, k + n] < 0) && isTRUE(fval[k + n] < fval[k])) {
      fval[c(k, k + n)] <- fval[c(k + n, k)]
      xpt[, c(k, k + n)] <- xpt[, c(k + n, k)]
    }
  }
  ij <- .nw_setij(n, npt)
  if (npt > 2L * n + 1L) {
    for (k in (2L * n + 2L):npt) xpt[, k] <- xpt[, ij[1L, k - 2L * n - 1L] + 1L] + xpt[, ij[2L, k - 2L * n - 1L] + 1L]
    if (info == 0L) for (k in (2L * n + 2L):npt) {
      info <- evalk(k)
      if (info != 0L) break
    }
  }
  kopt <- .bq_minloc_mask(fval, evaluated)
  list(ij = ij, kopt = kopt, nf = sum(evaluated), fval = fval, sl = sl, su = su, xbase = xbase, xpt = xpt,
       info = info)
}

# initq: the initial quadratic model; returns list(gopt, hq, pq)
.bq_initq <- function(ij, fval, xpt) {
  n <- nrow(xpt)
  npt <- ncol(xpt)
  fbase <- fval[1L]
  dg <- function(cols) xpt[cbind(seq_along(cols), cols)]
  gopt <- (fval[seq_len(n) + 1L] - fbase) / dg(seq_len(n) + 1L)
  ndiag <- min(n, npt - n - 1L)
  kd <- seq_len(ndiag)
  xa <- dg(kd + 1L)
  xb <- dg(n + kd + 1L)
  if (ndiag) gopt[kd] <- .pw_fma(gopt[kd], xb, -(((fval[n + kd + 1L] - fbase) / xb) * xa)) / (xb - xa)
  hq <- matrix(0, n, n)
  for (k in kd) hq[k, k] <- (2 * ((fval[k + 1L] - fbase) / xa[k] - (fval[n + k + 1L] - fbase) / xb[k])) / (xa[k] - xb[k])
  for (k in seq_len(max(0L, npt - 2L * n - 1L))) {
    i <- ij[1L, k]
    j <- ij[2L, k]
    hq[i, j] <- (((fbase - fval[i + 1L]) - fval[j + 1L]) + fval[k + 2L * n + 1L]) /
      (xpt[i, k + 2L * n + 1L] * xpt[j, k + 2L * n + 1L])
    hq[j, i] <- hq[i, j]
  }
  kopt <- .pw_minloc(fval)
  if (kopt != 1L) gopt <- gopt + .pw_matprod21(hq, xpt[, kopt])
  list(gopt = gopt, hq = hq, pq = numeric(npt))
}

# inith: the initial [BMAT, ZMAT] (H = ZMAT ZMAT' in BOBYQA)
.bq_inith <- function(ij, xpt) {
  n <- nrow(xpt)
  npt <- ncol(xpt)
  rhobeg <- max(abs(xpt[, 2L]))
  rhosq <- rhobeg * rhobeg
  ndiag <- min(n, npt - n - 1L)
  kd <- seq_len(ndiag)
  xa <- xpt[cbind(kd, kd + 1L)]
  xb <- xpt[cbind(kd, n + kd + 1L)]
  bmat <- matrix(0, n, npt + n)
  bmat[kd, 1L] <- -(xa + xb) / (xa * xb)
  for (k in kd) {
    bmat[k, k + n + 1L] <- -0.5 / xpt[k, k + 1L]
    bmat[k, k + 1L] <- -bmat[k, 1L] - bmat[k, k + n + 1L]
  }
  for (k in seq_len(n)[seq_len(n) > ndiag]) {
    bmat[k, 1L] <- -1 / xpt[k, k + 1L]
    bmat[k, k + 1L] <- -bmat[k, 1L]
    bmat[k, npt + k] <- -0.5 * rhosq
  }
  zmat <- matrix(0, npt, npt - n - 1L)
  zmat[1L, kd] <- sqrt(2) / (xa * xb)
  for (k in kd) {
    zmat[k + 1L, k] <- -zmat[1L, k] - sqrt(0.5) / rhosq
    zmat[k + n + 1L, k] <- sqrt(0.5) / rhosq
  }
  for (k in seq_len(npt - n - 1L)[seq_len(npt - n - 1L) > ndiag]) {
    zmat[1L, k] <- 1 / rhosq
    zmat[k + n + 1L, k] <- 1 / rhosq
    zmat[ij[, k - n] + 1L, k] <- -1 / rhosq
  }
  list(bmat = bmat, zmat = zmat)
}

# interval_max of univar: the maximizer of fun on [lb, ub] by a grid and a
# parabola through the neighbors of the best grid point
.bq_interval_max <- function(fun, lb, ub, args, grid_size) {
  if (ub <= lb) return(lb)
  xgrid <- .pw_linspace(lb, ub, grid_size)
  fgrid <- vapply(xgrid, fun, 0, args = args)
  if (all(is.nan(fgrid))) return(lb)
  kopt <- .pw_maxloc(fgrid)
  fopt <- fgrid[kopt]
  if (kopt == 1L) return(lb)
  if (kopt == grid_size) return(ub)
  fprev <- fgrid[kopt - 1L]
  fnext <- fgrid[kopt + 1L]
  step <- 0
  if (isTRUE(abs(fprev - fnext) > 0)) step <- 0.5 * ((fnext - fprev) / (((fopt + fopt) - fprev) - fnext))
  if (is.finite(step) && abs(step) > 0) lb + ((ub - lb) * ((kopt - 1) + step)) / (grid_size - 1) else xgrid[kopt]
}
.bq_interval_fun_trsbox <- function(hangt, args) {
  if (!(abs(hangt) > 0)) return(0)
  sth <- (hangt + hangt) / .pw_fma(hangt, hangt, 1)
  f <- .pw_fma(hangt, .pw_fma(hangt, args[2L], -args[3L]) - args[3L], args[1L])
  sth * .pw_fma(-(0.5 * sth), f, .pw_fma(hangt, args[4L], -args[5L]))
}

# trsbox: the truncated conjugate gradient method in the box, continued by
# rotations in two dimensional subspaces; returns list(d, crvmin)
.bq_trsbox <- function(delta, gopt_in, hq_in, pq_in, sl, su, xopt, xpt) {
  n <- length(gopt_in)
  ctest <- 0.01
  gmax <- max(abs(gopt_in))
  if (isTRUE(gmax > 999999995904)) {
    modscal <- .pw_fmaxnm(2 * .Machine$double.xmin, 1 / gmax)
    gopt <- gopt_in * modscal
    pq <- pq_in * modscal
    hq <- hq_in * modscal
    scaled <- TRUE
  } else {
    modscal <- 1
    gopt <- gopt_in
    pq <- pq_in
    hq <- hq_in
    scaled <- FALSE
  }
  sumsq <- function(v, m) .pw_inprod(v[m], v[m])
  xbdi <- integer(n)
  xbdi[xopt >= su & gopt <= 0] <- 1L
  xbdi[xopt <= sl & gopt >= 0] <- -1L
  nact <- sum(xbdi != 0L)
  d <- numeric(n)
  s <- numeric(n)
  crvmin <- -.pr_realmax
  gnew <- gopt
  gredsq <- sumsq(gnew, xbdi == 0L)
  delsq <- delta * delta
  qred <- 0
  beta <- 0
  ggsav <- 0
  itercg <- 0L
  twod_search <- FALSE
  maxiter <- min(10000L, (n - nact)^2)
  iter <- 0L
  while (iter < maxiter) {
    iter <- iter + 1L
    free <- xbdi == 0L
    resid <- delsq - sumsq(d, free)
    if (resid <= 0) {
      twod_search <- TRUE
      break
    }
    s <- if (itercg == 0L) -gnew else .pw_fma(beta, s, -gnew)
    s[!free] <- 0
    stepsq <- .pw_inprod(s, s)
    ds <- .pw_inprod(d[free], s[free])
    if (!(isTRUE(stepsq > .pr_eps * delsq) && isTRUE(gredsq * delsq > (ctest * qred) * (ctest * qred)) && !is.nan(ds)))
      break
    sqrtd <- .uo_maxval(c(sqrt(.pw_fma(stepsq, resid, ds * ds)), sqrt(stepsq * resid), abs(ds)))
    bstep <- if (ds >= 0) resid / (sqrtd + ds) else (sqrtd - ds) / stepsq
    if (!isTRUE(bstep > 0) || !is.finite(bstep)) break
    hs <- .pw_hess_mul(s, xpt, pq, hq)
    shs <- .pw_inprod(s[free], hs[free])
    stplen <- bstep
    if (isTRUE(shs > 0)) stplen <- .pw_fminnm(bstep, gredsq / shs)
    xnew <- xopt + d
    xtest <- .pw_fma(stplen, s, xnew)
    sbound <- rep(stplen, n)
    up <- which(s > 0 & xtest > su)
    sbound[up] <- (su[up] - xnew[up]) / s[up]
    lo <- which(s < 0 & xtest < sl)
    sbound[lo] <- (sl[lo] - xnew[lo]) / s[lo]
    sbound[is.nan(sbound)] <- stplen
    iact <- 0L
    if (any(sbound < stplen)) {
      iact <- .pw_minloc(sbound)
      stplen <- sbound[iact]
    }
    sdec <- 0
    if (isTRUE(stplen > 0)) {
      itercg <- itercg + 1L
      rayleighq <- shs / stepsq
      if (iact == 0L && isTRUE(rayleighq > 0)) {
        crvmin <- if (crvmin <= -.pr_realmax) rayleighq else .pw_fminnm(crvmin, rayleighq)
      }
      ggsav <- gredsq
      gnew <- .pw_fma(stplen, hs, gnew)
      gredsq <- sumsq(gnew, xbdi == 0L)
      d <- .pw_fma(stplen, s, d)
      sdec <- .pw_fmaxnm(stplen * .pw_fma(-0.5 * stplen, shs, ggsav), 0)
      qred <- qred + sdec
    }
    if (iact > 0L) {
      nact <- nact + 1L
      xbdi[iact] <- as.integer(.bq_copysign(1, s[iact]))
      if (nact >= n) break
      delsq <- .pw_fma(-d[iact], d[iact], delsq)
      if (delsq <= 0) {
        twod_search <- TRUE
        break
      }
      beta <- 0
      itercg <- 0L
      gredsq <- sumsq(gnew, xbdi == 0L)
    } else if (isTRUE(stplen < bstep)) {
      if (itercg >= n - nact || !isTRUE(sdec > ctest * qred) || is.nan(sdec) || is.nan(qred)) break
      beta <- gredsq / ggsav
    } else {
      twod_search <- TRUE
      break
    }
  }
  if (twod_search) {
    crvmin <- 0
    maxiter <- 10L * (n - nact)
  } else {
    maxiter <- 0L
  }
  nactsav <- nact - 1L
  hdred <- numeric(n)
  dredsq <- 0
  for (iter in seq_len(maxiter)) {
    xnew <- xopt + d
    xbdi[xbdi == 0L & xnew >= su] <- 1L
    xbdi[xbdi == 0L & xnew <= sl] <- -1L
    nact <- sum(xbdi != 0L)
    if (nact >= n - 1L) break
    free <- xbdi == 0L
    gredsq <- sumsq(gnew, free)
    dredg <- .pw_inprod(d[free], gnew[free])
    if (iter == 1L || nact > nactsav) {
      dredsq <- sumsq(d, free)
      dred <- d
      dred[!free] <- 0
      hdred <- .pw_hess_mul(dred, xpt, pq, hq)
      nactsav <- nact
    }
    temp <- .pw_fma(gredsq, dredsq, -(dredg * dredg))
    if (!isTRUE(temp > ((ctest * ctest) * qred) * qred) || is.nan(temp) || is.nan(qred)) break
    temp <- sqrt(temp)
    s <- .pw_fma(dredg, d, -(dredsq * gnew)) / temp
    s[!free] <- 0
    sredg <- -temp
    ssq <- .pw_fma(d, d, s * s)
    tanbd <- rep(1, n)
    sqdscr <- rep(-.pr_realmax, n)
    m <- which(free & xopt - sl < sqrt(ssq))
    sqdscr[m] <- sqrt(.bq_pmax(0, .pw_fma(-(xopt[m] - sl[m]), xopt[m] - sl[m], ssq[m])))
    m <- which(sqdscr - s > 0)
    tanbd[m] <- .bq_pmin(tanbd[m], (xnew[m] - sl[m]) / (sqdscr[m] - s[m]))
    sqdscr <- rep(-.pr_realmax, n)
    m <- which(free & su - xopt < sqrt(ssq))
    sqdscr[m] <- sqrt(.bq_pmax(0, .pw_fma(-(su[m] - xopt[m]), su[m] - xopt[m], ssq[m])))
    m <- which(sqdscr + s > 0)
    tanbd[m] <- .bq_pmin(tanbd[m], (su[m] - xnew[m]) / (sqdscr[m] + s[m]))
    tanbd[is.nan(tanbd)] <- 0
    iact <- 0L
    hangt_bd <- 1
    if (any(tanbd < 1)) {
      iact <- .pw_minloc(tanbd)
      hangt_bd <- tanbd[iact]
    }
    if (hangt_bd <= 0) break
    hs <- .pw_hess_mul(s, xpt, pq, hq)
    shs <- .pw_inprod(s[free], hs[free])
    dhs <- .pw_inprod(d[free], hs[free])
    dhd <- .pw_inprod(d[free], hdred[free])
    args <- c(shs, dhd, dhs, dredg, sredg)
    if (any(is.nan(args))) break
    grid_size <- 2L * as.integer(floor(.pw_fma(17, hangt_bd, 4.1) + 0.5))
    hangt <- .bq_interval_max(.bq_interval_fun_trsbox, 0, hangt_bd, args, grid_size)
    sdec <- .bq_interval_fun_trsbox(hangt, args)
    if (!isTRUE(sdec > 0)) break
    cth <- .pw_fma(-hangt, hangt, 1) / .pw_fma(hangt, hangt, 1)
    sth <- (hangt + hangt) / .pw_fma(hangt, hangt, 1)
    gnew <- .pw_fma(sth, hs, .pw_fma(cth - 1, hdred, gnew))
    d[free] <- .pw_fma(cth, d[free], sth * s[free])
    hdred <- .pw_fma(cth, hdred, sth * hs)
    qred <- qred + sdec
    if (iact >= 1L && iact <= n && hangt >= hangt_bd) {
      xbdi[iact] <- as.integer(.bq_copysign(1, .pw_fma(-0.5, sl[iact] + su[iact], xopt[iact] + d[iact])))
    } else if (!isTRUE(sdec > ctest * qred)) {
      break
    }
  }
  xnew <- .bq_pmax(sl, .bq_pmin(su, xopt + d))
  xnew[xbdi == -1L] <- sl[xbdi == -1L]
  xnew[xbdi == 1L] <- su[xbdi == 1L]
  d <- xnew - xopt
  if (crvmin <= -.pr_realmax || is.nan(crvmin)) crvmin <- 0
  if (scaled && crvmin > 0) crvmin <- crvmin / modscal
  list(d = d, crvmin = crvmin)
}

# errbd: the bound on the error of the model used to accept a short step
.bq_errbd <- function(crvmin, d, gopt, hq, moderr_rec, pq, rho, sl, su, xopt, xpt) {
  n <- nrow(xpt)
  xnew <- xopt + d
  gnew <- gopt + .pw_hess_mul(d, xpt, pq, hq)
  bfirst <- rep(.uo_maxval(abs(moderr_rec)), n)
  lo <- which(xnew <= sl)
  bfirst[lo] <- gnew[lo] * rho
  hi <- which(xnew >= su)
  bfirst[hi] <- -gnew[hi] * rho
  bsecond <- (0.5 * (diag(hq) + .pw_matprod21(xpt * xpt, pq))) * (rho * rho)
  ebound <- .uo_minval(.bq_pmax(bfirst, bfirst + bsecond))
  if (crvmin > 0) ebound <- .pw_fminnm(ebound, (0.125 * crvmin) * (rho * rho))
  ebound
}

# setdrop_tr of BOBYQA: the point to drop for the trust region step
.bq_setdrop_tr <- function(kopt, ximproved, bmat, d, delta, rho, xpt, zmat) {
  x <- if (ximproved) xpt[, kopt] + d else xpt[, kopt]
  distsq <- .pw_distsq(xpt, x)
  w <- .pw_vfmaxnm(distsq / (rho * rho), 1)
  w2 <- w * w
  score <- (w2 * w2) * .pw_calden(kopt, bmat, d, xpt, zmat, 1L)
  if (!ximproved) score[kopt] <- -1
  if (any(score > 1, na.rm = TRUE) || (ximproved && any(score > 0, na.rm = TRUE))) return(.pw_maxloc(score))
  if (ximproved) return(.pw_maxloc(distsq))
  0L
}

# geostep of BOBYQA: a step along a line through XOPT and an interpolation
# point, or a Cauchy step, that improves the geometry when XPT(:, knew) is
# replaced
.bq_geostep <- function(knew, kopt, bmat, delbar, sl, su, xpt, zmat) {
  n <- nrow(xpt)
  npt <- ncol(xpt)
  pqlag <- .pw_omega_col(1L, zmat, knew)
  alpha <- pqlag[knew]
  xopt <- xpt[, kopt]
  glag <- bmat[, knew] + .pw_hess_mul(xopt, xpt, pqlag)
  if (is.nan(sum(abs(glag)))) {
    d <- xpt[, knew] - xopt
    return(.pw_fminnm(0.5, delbar / .pw_norm(d)) * d)
  }
  dderiv <- .pw_matprod12(glag, xpt) - .pw_inprod(glag, xopt)
  distsq <- .pw_distsq(xpt, xopt)
  stplen <- matrix(0, 3L, npt)
  isbd <- matrix(0L, 3L, npt)
  for (k in seq_len(npt)) {
    if (k == kopt || is.nan(dderiv[k])) {
      dderiv[k] <- 0
      next
    }
    subd <- delbar / sqrt(distsq[k])
    slbd <- -subd
    ilbd <- 0L
    iubd <- 0L
    sumin <- .pw_fminnm(1, subd)
    xdiff <- xpt[, k] - xopt
    lfrac <- .bq_copysign(subd, -xdiff)
    m <- which(sl - xopt > -abs(xdiff) * subd)
    lfrac[m] <- (sl[m] - xopt[m]) / xdiff[m]
    ufrac <- .bq_copysign(subd, xdiff)
    m <- which(su - xopt < abs(xdiff) * subd)
    ufrac[m] <- (su[m] - xopt[m]) / xdiff[m]
    slbd_test <- rep(slbd, n)
    slbd_test[which(xdiff > 0)] <- lfrac[which(xdiff > 0)]
    slbd_test[which(xdiff < 0)] <- ufrac[which(xdiff < 0)]
    if (any(slbd_test > slbd, na.rm = TRUE)) {
      ilbd <- .pw_maxloc(slbd_test)
      slbd <- slbd_test[ilbd]
      ilbd <- -ilbd * as.integer(.bq_copysign(1, xdiff[ilbd]))
    }
    subd_test <- rep(subd, n)
    subd_test[which(xdiff > 0)] <- ufrac[which(xdiff > 0)]
    subd_test[which(xdiff < 0)] <- lfrac[which(xdiff < 0)]
    if (any(subd_test < subd, na.rm = TRUE)) {
      iubd <- .pw_minloc(subd_test)
      subd <- .pw_fmaxnm(sumin, subd_test[iubd])
      iubd <- iubd * as.integer(.bq_copysign(1, xdiff[iubd]))
    }
    stpm <- 0.5
    if (k == knew) {
      stpm <- slbd
      if (isTRUE(abs(1 - dderiv[k]) > 0)) stpm <- -((0.5 * dderiv[k]) / (1 - dderiv[k]))
    }
    stpm <- .pw_fmaxnm(slbd, .pw_fminnm(subd, stpm))
    stplen[, k] <- c(slbd, subd, stpm)
    isbd[, k] <- c(ilbd, iubd, 0L)
  }
  dd <- rep(dderiv, each = 3L)
  vlag <- matrix((stplen * (1 - stplen)) * dd, 3L, npt)
  vlag[, knew] <- stplen[, knew] * .pw_fma(stplen[, knew], 1 - dderiv[knew], dderiv[knew])
  vlag[is.nan(vlag)] <- 0
  t1 <- (stplen * (1 - stplen)) * rep(distsq, each = 3L)
  betabd <- 0.5 * (t1 * t1)
  predsq <- (vlag * vlag) * .pw_fma(vlag, vlag, alpha * betabd)
  predsq[is.nan(predsq)] <- 0
  ksqs <- apply(predsq, 1L, .pw_maxloc)
  isq <- .pw_maxloc(predsq[cbind(1:3, ksqs)])
  ksq <- ksqs[isq]
  stpsiz <- stplen[isq, ksq]
  ibd <- isbd[isq, ksq]
  xline <- .bq_pmax(sl, .bq_pmin(su, .pw_fma(stpsiz, xpt[, ksq] - xopt, xopt)))
  if (ibd < 0L) xline[-ibd] <- sl[-ibd]
  if (ibd > 0L) xline[ibd] <- su[ibd]
  d <- xline - xopt
  den_line <- .pw_calden(kopt, bmat, d, xpt, zmat, 1L)
  if (delbar > 0.009999999776482582) return(d)
  bigstp <- delbar + delbar
  xcauchy <- xopt
  vlagsq_cauchy <- 0
  x <- numeric(n)
  for (uphill in 0:1) {
    if (uphill == 1L) glag <- -glag
    s <- numeric(n)
    mask_free <- .bq_pmin(xopt - sl, glag) > 0 | .bq_pmax(xopt - su, glag) < 0
    s[which(mask_free)] <- bigstp
    ggfree <- .pw_inprod(glag[which(mask_free)], glag[which(mask_free)])
    if (!isTRUE(ggfree > 0)) next
    sfixsq <- 0
    grdstp <- 0
    for (k in seq_len(n)) {
      resis <- delbar * delbar - sfixsq
      if (resis <= 0) break
      ssqsav <- sfixsq
      grdstp <- sqrt(resis / ggfree)
      xtemp <- .pw_fma(-grdstp, glag, xopt)
      fixl <- s >= bigstp & xtemp <= sl
      fixu <- s >= bigstp & xtemp >= su
      mask_free <- s >= bigstp & !(fixl | fixu)
      s[which(fixl)] <- sl[which(fixl)] - xopt[which(fixl)]
      s[which(fixu)] <- su[which(fixu)] - xopt[which(fixu)]
      fx <- which(fixl | fixu)
      sfixsq <- sfixsq + .pw_inprod(s[fx], s[fx])
      ggfree <- .pw_inprod(glag[which(mask_free)], glag[which(mask_free)])
      if (!(isTRUE(sfixsq > ssqsav) && isTRUE(ggfree > 0))) break
    }
    x[which(glag > 0)] <- sl[which(glag > 0)]
    x[which(glag <= 0)] <- su[which(glag <= 0)]
    x[which(abs(s) <= 0)] <- xopt[which(abs(s) <= 0)]
    xtemp <- .bq_pmax(sl, .bq_pmin(su, .pw_fma(-grdstp, glag, xopt)))
    big <- which(s >= bigstp)
    x[big] <- xtemp[big]
    s[big] <- -grdstp * glag[big]
    gs <- .pw_inprod(glag, s)
    sxpt <- .pw_matprod12(s, xpt)
    curv <- .pw_inprod(sxpt, pqlag * sxpt)
    if (uphill == 1L) curv <- -curv
    if (isTRUE(curv > -gs) && isTRUE(curv < -2.414213562373095 * gs)) {
      scaling <- -gs / curv
      x <- .bq_pmax(sl, .bq_pmin(su, .pw_fma(scaling, s, xopt)))
      v <- (0.5 * gs) * scaling
      vlagsq <- v * v
    } else {
      v <- .pw_fma(0.5, curv, gs)
      vlagsq <- v * v
    }
    if (isTRUE(vlagsq > vlagsq_cauchy)) {
      xcauchy <- x
      vlagsq_cauchy <- vlagsq
    }
  }
  s <- xcauchy - xopt
  den_cauchy <- .pw_calden(kopt, bmat, s, xpt, zmat, 1L)
  if (isTRUE(den_cauchy[knew] > .pw_fmaxnm(den_line[knew], 0)) || is.nan(den_line[knew])) d <- s
  d
}

# tryqalt of BOBYQA: the alternative model of least Frobenius norm, compared
# through the projected gradients
.bq_tryqalt <- function(bmat, fval, ratio, sl, su, xopt, xpt, zmat, itest, gopt, hq, pq) {
  npt <- ncol(xpt)
  proj <- function(g) {
    up <- which(xopt >= su)
    g[up] <- .bq_pmax(0, g[up])
    lo <- which(xopt <= sl)
    g[lo] <- .bq_pmin(0, g[lo])
    g
  }
  pgopt <- proj(gopt)
  pqalt <- .pw_omega_mul(1L, zmat, fval)
  galt <- .pw_matprod21(bmat[, seq_len(npt), drop = FALSE], fval) + .pw_hess_mul(xopt, xpt, pqalt)
  pgalt <- proj(galt)
  if (isTRUE(ratio > 0.1) || isTRUE(.pw_inprod(pgopt, pgopt) < 10 * .pw_inprod(pgalt, pgalt))) {
    itest <- 0L
  } else {
    itest <- itest + 1L
  }
  if (itest >= 3L) {
    gopt <- galt
    pq <- pqalt
    hq[] <- 0
    itest <- 0L
  }
  list(itest = itest, gopt = gopt, hq = hq, pq = pq)
}

# updateh_rsc of rescue: the update of [BMAT, ZMAT] with given BETA and VLAG
.bq_updateh_rsc <- function(knew, beta, vlag, bmat, zmat) {
  npt <- nrow(zmat)
  n <- nrow(bmat)
  tau <- vlag[knew]
  zk <- zmat[knew, ]
  denom <- .pw_fma(.pw_inprod(zk, zk), beta, tau * tau)
  if (!(is.finite(sum(abs(vlag)) + abs(beta)) && isTRUE(denom > 0))) return(list(bmat = bmat, zmat = zmat))
  vlag[knew] <- vlag[knew] - 1
  ztest <- 1e-20 * max(abs(zmat))
  for (j in seq_len(npt - n - 1L)[-1L]) {
    if (isTRUE(abs(zmat[knew, j]) > ztest)) {
      grot <- .pw_planerot(c(zmat[knew, 1L], zmat[knew, j]))
      zmat[, c(1L, j)] <- .pw_matprod22(zmat[, c(1L, j), drop = FALSE], t(grot))
    }
    zmat[knew, j] <- 0
  }
  hcol <- c(zmat[knew, 1L] * zmat[, 1L], bmat[, knew])
  sqrtdn <- sqrt(denom)
  vk <- vlag[seq_len(npt)]
  zmat[, 1L] <- .pw_fma(tau / sqrtdn, zmat[, 1L], -((zmat[knew, 1L] / sqrtdn) * vk))
  alpha <- hcol[knew]
  jj <- npt + seq_len(n)
  v1 <- .pw_fma(alpha, vlag[jj], -(tau * hcol[jj])) / denom
  v2 <- .pw_fma(-beta, hcol[jj], -(tau * vlag[jj])) / denom
  bmat <- (bmat + outer(v1, vlag)) + outer(v2, hcol)
  B <- bmat[, jj, drop = FALSE]
  B[upper.tri(B)] <- t(B)[upper.tri(B)]
  bmat[, jj] <- B
  list(bmat = bmat, zmat = zmat)
}

# rescue: rebuilds the interpolation set and [BMAT, ZMAT] when rounding has
# spoiled the denominators (Powell's RESCUE); st holds the state of bobyqb
.bq_rescue <- function(fun, iprint, maxfun, delta, ftarget, xl, xu, st) {
  xpt <- st$xpt
  n <- nrow(xpt)
  npt <- ncol(xpt)
  kopt <- st$kopt
  st$info <- 0L
  if (st$nf >= maxfun) {
    st$bmat[] <- 0
    st$zmat[] <- 0
    st$info <- 3L
    return(st)
  }
  xopt <- xpt[, kopt]
  sl <- .bq_pmin(st$sl - xopt, 0)
  su <- .bq_pmax(st$su - xopt, 0)
  xbase <- .bq_pmin(.bq_pmax(xl, st$xbase + xopt), xu)
  xpt <- xpt - xopt
  xpt[, kopt] <- 0
  pq <- st$pq
  spq <- 0
  for (v in pq) spq <- spq + v
  v <- .pw_fma(0.5 * spq, xopt, .pw_matprod21(xpt, pq))
  hq <- .bq_r2sym(st$hq, 1, xopt, v)
  gopt <- st$gopt
  fval <- st$fval
  p1 <- .bq_pmin(delta, su)
  p2 <- .bq_pmax(-delta, sl)
  sw <- which(p1 + p2 < 0)
  t1 <- p1[sw]
  p1[sw] <- p2[sw]
  p2[sw] <- t1
  m <- which(abs(p2) < 0.5 * abs(p1))
  p2[m] <- 0.5 * p1[m]
  sfrac <- 0.5 / (n + 1)
  ptsid <- numeric(npt)
  ptsid[1L] <- sfrac
  bmat <- matrix(0, n, npt + n)
  zmat <- matrix(0, npt, npt - n - 1L)
  for (k in seq_len(n)) {
    ptsid[k + 1L] <- k + sfrac
    if (k <= npt - n - 1L) {
      ptsid[k + n + 1L] <- k / (n + 1) + sfrac
      temp <- 1 / (p1[k] - p2[k])
      bmat[k, k + 1L] <- -temp + 1 / p1[k]
      bmat[k, k + n + 1L] <- temp + 1 / p2[k]
      bmat[k, 1L] <- -bmat[k, k + 1L] - bmat[k, k + n + 1L]
      zmat[1L, k] <- sqrt(2) / abs(p1[k] * p2[k])
      zmat[k + 1L, k] <- (zmat[1L, k] * p2[k]) * temp
      zmat[k + n + 1L, k] <- (-zmat[1L, k] * p1[k]) * temp
    } else {
      bmat[k, 1L] <- -1 / p1[k]
      bmat[k, k + 1L] <- 1 / p1[k]
      bmat[k, k + npt] <- -0.5 * (p1[k] * p1[k])
    }
  }
  ij <- .nw_setij(n, npt)
  for (k in seq_len(npt)[seq_len(npt) >= 2L * n + 2L]) {
    ip <- ij[1L, k - 2L * n - 1L]
    iq <- ij[2L, k - 2L * n - 1L]
    ptsid[k] <- (ip + iq / (n + 1)) + sfrac
    temp <- 1 / (p1[ip] * p1[iq])
    zmat[c(1L, k), k - n - 1L] <- temp
    zmat[c(ip + 1L, iq + 1L), k - n - 1L] <- -temp
  }
  if (kopt != 1L) {
    bmat[, c(1L, kopt)] <- bmat[, c(kopt, 1L)]
    zmat[c(1L, kopt), ] <- zmat[c(kopt, 1L), ]
  }
  ptsid[1L] <- ptsid[kopt]
  ptsid[kopt] <- 0
  score <- sqrt(.pw_distsq(xpt, numeric(n)))
  score[kopt] <- 0
  scoreinc <- max(score)
  nprov <- npt - 1L
  # the indices ip and iq of a provisional point from ptsid
  pidx <- function(p) {
    ip <- floor(p)
    c(ip, floor(.pw_fma(n + 1, p, -((n + 1) * ip))))
  }
  wmv <- numeric(npt + n)
  hdiag <- numeric(npt)
  for (iter in seq_len(npt * npt)) {
    if (all(score <= 0) || nprov <= 1L) break
    korig <- .bq_minloc_mask(score, score > 0)
    for (k in seq_len(npt)) {
      if (k == kopt) {
        w <- 0
      } else if (ptsid[k] <= 0) {
        w <- .pw_inprod(xpt[, korig], xpt[, k])
      } else {
        pq2 <- pidx(ptsid[k])
        ip <- pq2[1L]
        iq <- pq2[2L]
        w <- if (ip > 0 && iq > 0) .pw_fma(xpt[ip, korig], p1[ip], xpt[iq, korig] * p1[iq])
          else if (ip > 0) xpt[ip, korig] * p1[ip]
          else if (iq > 0) xpt[iq, korig] * p2[iq]
          else 0
      }
      wmv[k] <- (0.5 * w) * w
    }
    wmv[npt + seq_len(n)] <- xpt[, korig]
    w1 <- wmv[seq_len(npt)]
    w1z <- .pw_matprod12(w1, zmat)
    vlag <- numeric(npt + n)
    vlag[seq_len(npt)] <- .pw_matprod21(zmat, w1z) + .pw_matprod12(wmv[npt + seq_len(n)], bmat[, seq_len(npt), drop = FALSE])
    vlag[npt + seq_len(n)] <- .pw_matprod21(bmat, wmv)
    bsum <- .pw_inprod(wmv[seq_len(n)], .pw_matprod21(bmat[, seq_len(npt), drop = FALSE], w1) + .pw_matprod21(bmat, wmv))
    xs <- .pw_inprod(xpt[, korig], xpt[, korig])
    beta <- .pw_fma(0.5, xs * xs, -.pw_inprod(w1z, w1z)) - bsum
    vlag[kopt] <- vlag[kopt] + 1
    den <- numeric(npt)
    pr <- which(ptsid > 0)
    for (r in pr) hdiag[r] <- .pw_inprod(zmat[r, ], zmat[r, ])
    den[pr] <- .pw_fma(hdiag[pr], beta, vlag[pr] * vlag[pr])
    vk <- vlag[seq_len(npt)]
    if (!(is.finite(sum(abs(vlag))) && any(den > 0.05 * .uo_maxval(vk * vk), na.rm = TRUE))) {
      score[korig] <- -score[korig] - scoreinc
      next
    }
    kprov <- .pw_maxloc(den)
    if (kprov != korig) {
      bmat[, c(kprov, korig)] <- bmat[, c(korig, kprov)]
      zmat[c(kprov, korig), ] <- zmat[c(korig, kprov), ]
      vlag[c(kprov, korig)] <- vlag[c(korig, kprov)]
    }
    ptsid[kprov] <- ptsid[korig]
    ptsid[korig] <- 0
    score[korig] <- 0
    score <- abs(score)
    up <- .bq_updateh_rsc(korig, beta, vlag, bmat, zmat)
    bmat <- up$bmat
    zmat <- up$zmat
    nprov <- nprov - 1L
  }
  kbase <- kopt
  fbase <- fval[kopt]
  nf <- st$nf
  if (nprov > 0L) for (kpt in seq_len(npt)) {
    if (ptsid[kpt] <= 0) next
    hq <- .pw_r1sym(hq, pq[kpt], xpt[, kpt])
    pq[kpt] <- 0
    pq2 <- pidx(ptsid[kpt])
    ip <- pq2[1L]
    iq <- pq2[2L]
    xp <- 0
    xq <- 0
    xpt[, kpt] <- 0
    if (ip > 0 && iq > 0) {
      xp <- p1[ip]
      xpt[ip, kpt] <- xp
      xq <- p1[iq]
      xpt[iq, kpt] <- xq
    } else if (ip > 0) {
      xp <- p1[ip]
      xpt[ip, kpt] <- xp
    } else if (iq > 0) {
      xq <- p2[iq]
      xpt[iq, kpt] <- xq
    }
    x <- .bq_xinbd(xbase, xpt[, kpt], xl, xu, sl, su)
    f <- .pw_evaluatef(fun, x)
    nf <- nf + 1L
    .pw_fmsg(iprint, "Rescue", nf, delta, f, x)
    fval[kpt] <- f
    if (isTRUE(f < fval[kopt])) kopt <- kpt
    sub <- .pw_checkexit(maxfun, nf, f, ftarget, x)
    if (sub != 0L) {
      st$info <- sub
      break
    }
    vquad <- fbase
    xxpt <- numeric(npt)
    if (ip > 0 && iq > 0) {
      vquad <- .pw_fma(xp, .pw_fma(0.5 * xp, hq[ip, ip], gopt[ip]), vquad)
      vquad <- .pw_fma(xq, .pw_fma(0.5 * xq, hq[iq, iq], gopt[iq]), vquad)
      vquad <- .pw_fma(xp * xq, hq[ip, iq], vquad)
      xxpt <- .pw_fma(xp, xpt[ip, ], xq * xpt[iq, ])
    } else if (ip > 0) {
      vquad <- .pw_fma(xp, .pw_fma(0.5 * xp, hq[ip, ip], gopt[ip]), vquad)
      xxpt <- xp * xpt[ip, ]
    } else if (iq > 0) {
      vquad <- .pw_fma(xq, .pw_fma(0.5 * xq, hq[iq, iq], gopt[iq]), vquad)
      xxpt <- xq * xpt[iq, ]
    }
    vquad <- .pw_fma(0.5, .pw_inprod(xxpt, pq * xxpt), vquad)
    moderr <- f - vquad
    gopt <- .pw_fma(moderr, bmat[, kpt], gopt)
    pqinc <- moderr * .pw_omega_col(1L, zmat, kpt)
    m <- which(ptsid <= 0)
    pq[m] <- pq[m] + pqinc[m]
    for (k in which(ptsid > 0)) {
      pq2 <- pidx(ptsid[k])
      ip <- pq2[1L]
      iq <- pq2[2L]
      if (ip > 0 && iq > 0) {
        hq[ip, ip] <- .pw_fma(pqinc[k], p1[ip] * p1[ip], hq[ip, ip])
        hq[iq, iq] <- .pw_fma(pqinc[k], p1[iq] * p1[iq], hq[iq, iq])
        hq[ip, iq] <- .pw_fma(pqinc[k] * p1[ip], p1[iq], hq[ip, iq])
        hq[iq, ip] <- hq[ip, iq]
      } else if (ip > 0) {
        hq[ip, ip] <- .pw_fma(pqinc[k], p1[ip] * p1[ip], hq[ip, ip])
      } else if (iq > 0) {
        hq[iq, iq] <- .pw_fma(pqinc[k], p2[iq] * p2[iq], hq[iq, iq])
      }
    }
    ptsid[kpt] <- 0
  }
  if (kopt != kbase) gopt <- gopt + .pw_hess_mul(xpt[, kopt], xpt, pq, hq)
  st[c("xpt", "kopt", "nf", "fval", "gopt", "hq", "pq", "sl", "su", "xbase", "bmat", "zmat")] <-
    list(xpt, kopt, nf, fval, gopt, hq, pq, sl, su, xbase, bmat, zmat)
  st
}

# bobyqb: the main loop of BOBYQA
.bq_bobyqb <- function(fun, iprint, maxfun, npt, ftarget, rhobeg, rhoend, xl, xu, x0) {
  eta1 <- 0.1
  eta2 <- 0.7
  gamma1 <- 0.5
  gamma2 <- 2
  n <- length(x0)
  ini <- .bq_initxf(fun, iprint, maxfun, ftarget, rhobeg, xl, xu, x0, npt)
  st <- ini[c("kopt", "nf", "fval", "sl", "su", "xbase", "xpt")]
  x <- .bq_xinbd(st$xbase, st$xpt[, st$kopt], xl, xu, st$sl, st$su)
  f <- st$fval[st$kopt]
  if (ini$info != 0L) {
    .pw_retmsg("BOBYQA", iprint, ini$info, st$nf, f, x)
    return(list(x = x, f = f, nf = st$nf, info = ini$info))
  }
  st[c("bmat", "zmat")] <- .bq_inith(ini$ij, st$xpt)
  st[c("gopt", "hq", "pq")] <- .bq_initq(ini$ij, st$fval, st$xpt)
  rho <- rhobeg
  delta <- rho
  ebound <- 0
  shortd <- FALSE
  trfail <- FALSE
  ratio <- -1
  dnorm_rec <- c(.pr_realmax, .pr_realmax)
  moderr_rec <- dnorm_rec
  knew_tr <- 0L
  itest <- 0L
  gamma3 <- .pw_fmaxnm(1, .pw_fminnm(0.75 * gamma2, 1.5))
  maxtr <- if (maxfun > .Machine$integer.max %/% 2L) maxfun else 2L * maxfun  # 2 maxfun overflows in Fortran
  info <- 20L
  d <- numeric(n)
  # replaces XPT(:, knew) by XOPT + d in [BMAT, ZMAT], XPT, FVAL and the model
  swap_point <- function(knew, ximproved, moderr) {
    xdrop <- st$xpt[, knew]
    xosav <- st$xpt[, st$kopt]
    h <- .pw_updateh(knew, st$kopt, d, st$xpt, 1L, st$bmat, st$zmat, posdef = TRUE)
    st$bmat <<- h$bmat
    st$zmat <<- h$zmat
    st$xpt[, knew] <<- .bq_pmax(st$sl, .bq_pmin(st$su, xosav + d))
    st$fval[knew] <<- f
    if (ximproved) st$kopt <<- knew
    q <- .nw_updateq(1L, knew, ximproved, st$bmat, d, moderr, xdrop, xosav, st$xpt, st$zmat, st$gopt, st$hq, st$pq)
    st[c("gopt", "hq", "pq")] <<- q
  }
  rescue <- function() {
    st <<- .bq_rescue(fun, iprint, maxfun, delta, ftarget, xl, xu, st)
    st$nf
  }
  evalx <- function(state, rad) {
    x <<- .bq_xinbd(st$xbase, st$xpt[, st$kopt] + d, xl, xu, st$sl, st$su)
    f <<- .pw_evaluatef(fun, x)
    st$nf <<- st$nf + 1L
    .pw_fmsg(iprint, state, st$nf, rad, f, x)
    .pw_checkexit(maxfun, st$nf, f, ftarget, x)
  }
  for (tr in seq_len(maxtr)) {
    kopt <- st$kopt
    ts <- .bq_trsbox(delta, st$gopt, st$hq, st$pq, st$sl, st$su, st$xpt[, kopt], st$xpt)
    d <- ts$d
    crvmin <- ts$crvmin
    dnorm <- .pw_fminnm(delta, .pw_norm(d))
    shortd <- isTRUE(dnorm < 0.5 * rho)
    qred <- -.pw_quadinc(d, st$xpt, st$gopt, st$pq, st$hq)
    trfail <- !isTRUE(qred > 1e-5 * (rho * rho))
    if (shortd || trfail) {
      delta <- 0.1 * delta
      if (delta <= gamma3 * rho) delta <- rho
      ebound <- .bq_errbd(crvmin, d, st$gopt, st$hq, moderr_rec, st$pq, rho, st$sl, st$su, st$xpt[, kopt], st$xpt)
    } else {
      sub <- evalx("Trust region", delta)
      if (sub != 0L) {
        info <- sub
        break
      }
      dnorm_rec <- c(dnorm_rec[2L], dnorm)
      fopt <- st$fval[kopt]
      moderr <- (f - fopt) + qred
      moderr_rec <- c(moderr_rec[2L], moderr)
      ratio <- .pr_redrat(fopt - f, qred, eta1)
      delta <- .bq_trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
      if (delta <= gamma3 * rho) delta <- rho
      ximproved <- isTRUE(f < fopt)
      vlag <- .pw_calvlag(kopt, st$bmat, d, st$xpt, st$zmat, 1L)
      den <- .pw_calden(kopt, st$bmat, d, st$xpt, st$zmat, 1L)
      vk <- vlag[seq_len(npt)]
      if (ximproved && !(is.finite(sum(abs(vlag))) && any(den > .uo_maxval(vk * vk), na.rm = TRUE))) {
        rescue()
        if (st$info != 0L) {
          info <- st$info
          break
        }
        dnorm_rec <- c(.pr_realmax, .pr_realmax)
        moderr_rec <- dnorm_rec
        d <- .bq_pmax(st$sl, .bq_pmin(st$su, d)) - st$xpt[, st$kopt]
        moderr <- (f - st$fval[st$kopt]) - .pw_quadinc(d, st$xpt, st$gopt, st$pq, st$hq)
        ximproved <- isTRUE(f < st$fval[st$kopt])
      }
      knew_tr <- .bq_setdrop_tr(st$kopt, ximproved, st$bmat, d, delta, rho, st$xpt, st$zmat)
      if (knew_tr > 0L) {
        swap_point(knew_tr, ximproved, moderr)
        qa <- .bq_tryqalt(st$bmat, st$fval - st$fval[st$kopt], ratio, st$sl, st$su, st$xpt[, st$kopt], st$xpt,
                          st$zmat, itest, st$gopt, st$hq, st$pq)
        itest <- qa$itest
        st[c("gopt", "hq", "pq")] <- qa[c("gopt", "hq", "pq")]
      }
    }
    kopt <- st$kopt
    accurate_mod <- all(abs(moderr_rec) <= ebound) && all(dnorm_rec <= rho)
    distsq <- .pw_distsq(st$xpt, st$xpt[, kopt])
    tr2 <- 10 * rho
    close_itpset <- all(distsq <= .pw_fmaxnm(delta * delta, tr2 * tr2))
    adequate_geo <- isTRUE((shortd && accurate_mod) || close_itpset)
    small_trrad <- isTRUE(.pw_fmaxnm(delta, dnorm) <= rho)
    improve_geo <- (shortd || trfail || !isTRUE(ratio > eta1) || knew_tr == 0L) && !adequate_geo
    reduce_rho <- (shortd || trfail || !isTRUE(ratio > 0) || knew_tr == 0L) && adequate_geo && small_trrad
    if (improve_geo) {
      knew_geo <- .pw_maxloc(distsq)
      delbar <- .pw_fmaxnm(.pw_fminnm(0.1 * sqrt(.uo_maxval(distsq)), delta), rho)
      d <- .bq_geostep(knew_geo, kopt, st$bmat, delbar, st$sl, st$su, st$xpt, st$zmat)
      vlag <- .pw_calvlag(kopt, st$bmat, d, st$xpt, st$zmat, 1L)
      den <- .pw_calden(kopt, st$bmat, d, st$xpt, st$zmat, 1L)
      if (!(is.finite(sum(abs(vlag))) && isTRUE(den[knew_geo] > 0.5 * (vlag[knew_geo] * vlag[knew_geo])))) {
        rescue()
        if (st$info != 0L) {
          info <- st$info
          break
        }
        dnorm_rec <- c(.pr_realmax, .pr_realmax)
        moderr_rec <- dnorm_rec
      } else {
        sub <- evalx("Geometry", delbar)
        if (sub != 0L) {
          info <- sub
          break
        }
        dnorm <- .pw_fminnm(delbar, .pw_norm(d))
        dnorm_rec <- c(dnorm_rec[2L], dnorm)
        moderr <- (f - st$fval[kopt]) - .pw_quadinc(d, st$xpt, st$gopt, st$pq, st$hq)
        moderr_rec <- c(moderr_rec[2L], moderr)
        swap_point(knew_geo, isTRUE(f < st$fval[kopt]), moderr)
      }
    }
    if (reduce_rho) {
      if (rho <= rhoend) {
        info <- 0L
        break
      }
      rhonew <- .pr_redrho(rho, rhoend)
      delta <- .pw_fmaxnm(0.5 * rho, rhonew)
      rho <- rhonew
      .pw_rhomsg(iprint, st$nf, delta, st$fval[st$kopt], rho, st$xbase + st$xpt[, st$kopt])
      dnorm_rec <- c(.pr_realmax, .pr_realmax)
      moderr_rec <- dnorm_rec
    }
    xopt <- st$xpt[, st$kopt]
    if (isTRUE(.pw_inprod(xopt, xopt) >= 1000 * (delta * delta))) {
      st$sl <- .bq_pmin(st$sl - xopt, 0)
      st$su <- .bq_pmax(st$su - xopt, 0)
      sb <- .pw_shiftbase(st$kopt, st$xbase, st$xpt, st$zmat, st$bmat, st$pq, st$hq, 1L)
      st$xpt <- sb$xpt
      st$bmat <- sb$bmat
      st$hq <- sb$hq
      st$xbase <- .bq_pmax(xl, .bq_pmin(xu, sb$xbase))
    }
  }
  if (info == 0L && shortd && st$nf < maxfun) evalx("Trust region", rho)
  if (isTRUE(st$fval[st$kopt] < f) || is.nan(f)) {
    x <- .bq_xinbd(st$xbase, st$xpt[, st$kopt], xl, xu, st$sl, st$su)
    f <- st$fval[st$kopt]
  }
  .pw_retmsg("BOBYQA", iprint, info, st$nf, f, x)
  list(x = x, f = f, nf = st$nf, info = info)
}

# bobyqa: the bounds and options as libprima prepares them, then bobyqb
.bq_bobyqa <- function(fun, x0, xl = NULL, xu = NULL, rhobeg = 1, rhoend = 1e-6 * rhobeg, ftarget = -Inf,
                       maxfun = 500L * length(x0), npt = 2L * length(x0) + 1L, iprint = 0L) {
  n <- length(x0)
  boundmax <- 0.25 * .pr_realmax
  xl <- if (is.null(xl)) rep(-boundmax, n) else as.numeric(xl)
  xu <- if (is.null(xu)) rep(boundmax, n) else as.numeric(xu)
  xl[is.nan(xl) | xl < -boundmax] <- -boundmax
  xu[is.nan(xu) | xu > boundmax] <- boundmax
  if (any(xu - xl < 2 * .pr_eps)) {
    warning("BOBYQA: There is no space between the lower and upper bounds of variable ",
            which(xu - xl < 2 * .pr_eps)[1L], ". The solver cannot continue", call. = FALSE)
    return(list(x = as.numeric(x0), f = NaN, nf = 0L, info = 6L))
  }
  x0 <- .bq_pmax(xl, .bq_pmin(xu, .pw_moderatex(as.numeric(x0))))
  op <- .pw_preproc("BOBYQA", n, maxfun, ftarget, rhobeg, rhoend, npt, xl = xl, xu = xu, x0 = x0)
  .bq_bobyqb(fun, as.integer(iprint), op$maxfun, op$npt, op$ftarget, op$rhobeg, op$rhoend, xl, xu, op$x0)
}
