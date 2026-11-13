# NEWUOA of PRIMA (Zaikun Zhang's modern reference implementation of
# Powell's NEWUOA), reconstructed from the compiled library that PRIMA.jl
# calls and checked against it bit for bit (see R/optimization-prima-powell.R
# for the method). Function names follow the Fortran routines of
# src/newuoa (newuob, initxf, initq, inith, trsapp, setdrop_tr, geostep,
# biglag, bigden, updateq, tryqalt).

# initxf: the initial interpolation set and its function values
.nw_setij <- function(n, npt) {
  m <- max(0L, npt - 2L * n - 1L)
  k <- seq_len(m)
  i <- (k - 1L) %% n + 1L
  matrix(c(rbind(i, (i + (k - 1L) %/% n) %% n + 1L)), 2L, m)
}
.nw_initxf <- function(fun, maxfun, ftarget, rhobeg, x0, npt, iprint = 0L) {
  n <- length(x0)
  xpt <- matrix(0, n, npt)
  for (k in seq_len(min(npt - 1L, n))) xpt[k, k + 1L] <- rhobeg
  for (k in seq_len(min(npt - n - 1L, n))) xpt[k, k + n + 1L] <- -rhobeg
  fval <- rep(.pr_realmax, npt)
  evaluated <- logical(npt)
  info <- 0L
  for (k in seq_len(min(npt, 2L * n + 1L))) {
    x <- x0 + xpt[, k]
    f <- .pw_evaluatef(fun, x)
    evaluated[k] <- TRUE
    fval[k] <- f
    .pw_fmsg(iprint, "Initialization", k, rhobeg, f, x)
    info <- .pw_checkexit(maxfun, k, f, ftarget, x)
    if (info != 0L) break
  }
  ij <- .nw_setij(n, npt)
  if (info == 0L && npt > 2L * n + 1L) {
    for (k in seq_len(ncol(ij))) for (r in 1:2) {
      if (fval[ij[r, k] + n + 1L] < fval[ij[r, k] + 1L]) ij[r, k] <- ij[r, k] + n
    }
    for (k in (2L * n + 2L):npt) xpt[, k] <- xpt[, ij[1L, k - 2L * n - 1L] + 1L] + xpt[, ij[2L, k - 2L * n - 1L] + 1L]
    for (k in (2L * n + 2L):npt) {
      x <- x0 + xpt[, k]
      f <- .pw_evaluatef(fun, x)
      evaluated[k] <- TRUE
      fval[k] <- f
      .pw_fmsg(iprint, "Initialization", k, rhobeg, f, x)
      info <- .pw_checkexit(maxfun, k, f, ftarget, x)
      if (info != 0L) break
    }
  }
  kopt <- which.min(ifelse(evaluated, fval, Inf))
  list(ij = ij, kopt = kopt, nf = sum(evaluated), fval = fval, xbase = x0, xpt = xpt, info = info)
}

# initq: the initial quadratic model; returns list(gopt, hq, pq)
.nw_initq <- function(ij, fval, xpt) {
  n <- nrow(xpt)
  npt <- ncol(xpt)
  rhobeg <- max(abs(xpt[, 2L]))
  fbase <- fval[1L]
  gq <- (fval[seq_len(n) + 1L] - fbase) / rhobeg
  hq <- matrix(0, n, n)
  ndiag <- min(npt - n - 1L, n)
  for (k in seq_len(ndiag)) {
    hq[k, k] <- ((fval[k + 1L] - fbase) / rhobeg - (fbase - fval[k + n + 1L]) / rhobeg) / rhobeg
    gq[k] <- 0.5 * (gq[k] + (fbase - fval[k + n + 1L]) / rhobeg)
  }
  if (npt > 2L * n + 1L) for (k in (2L * n + 2L):npt) {
    ij1 <- ij[1L, k - 2L * n - 1L]
    ij2 <- ij[2L, k - 2L * n - 1L]
    i <- (ij1 - 1L) %% n + 1L
    j <- (ij2 - 1L) %% n + 1L
    hq[i, j] <- (((fbase - fval[ij1 + 1L]) - fval[ij2 + 1L]) + fval[k]) / (xpt[i, k] * xpt[j, k])
    hq[j, i] <- hq[i, j]
  }
  kopt <- which.min(fval)
  gopt <- if (kopt != 1L) gq + .pw_matprod21(hq, xpt[, kopt]) else gq
  list(gopt = gopt, hq = hq, pq = numeric(npt))
}

# inith: the initial [BMAT, ZMAT, IDZ]
.nw_inith <- function(ij, xpt) {
  n <- nrow(xpt)
  npt <- ncol(xpt)
  rhobeg <- max(abs(xpt[, 2L]))
  rhosq <- rhobeg^2
  recip <- 1 / rhosq
  reciq <- sqrt(0.5) / rhosq
  bmat <- matrix(0, n, npt + n)
  zmat <- matrix(0, npt, npt - n - 1L)
  for (k in seq_len(n)) {
    if (k + n + 1L <= npt) {
      bmat[k, k + 1L] <- 0.5 / rhobeg
      bmat[k, k + n + 1L] <- -0.5 / rhobeg
    } else {
      bmat[k, 1L] <- -1 / rhobeg
      bmat[k, k + 1L] <- 1 / rhobeg
      bmat[k, npt + k] <- -0.5 * rhosq
    }
  }
  for (k in seq_len(min(npt - n - 1L, n))) {
    zmat[1L, k] <- -reciq - reciq
    zmat[k + 1L, k] <- reciq
    zmat[k + n + 1L, k] <- reciq
  }
  if (npt - n - 1L > n) for (k in (n + 1L):(npt - n - 1L)) {
    zmat[1L, k] <- recip
    zmat[k + n + 1L, k] <- recip
    zmat[ij[1L, k - n] + 1L, k] <- -recip
    zmat[ij[2L, k - n] + 1L, k] <- -recip
  }
  list(idz = 1L, bmat = bmat, zmat = zmat)
}

# trsapp: truncated conjugate gradient for the trust region subproblem,
# continued by searches on the boundary; returns list(d, crvmin)
.nw_circle_trsapp <- function(theta, args) {
  cth <- cos(theta)
  sth <- sin(theta)
  .pw_fma(.pw_fma(args[2L], cth, args[1L]), cth, .pw_fma(args[4L], cth, args[3L]) * sth)
}
.nw_trsapp <- function(delta, gopt, hq, pq, tol, xpt) {
  n <- length(gopt)
  maxabs <- max(abs(gopt))
  scaling <- 1
  scaled <- FALSE
  if (isTRUE(maxabs > 999999995904)) {
    scaling <- .pw_fmaxnm(1 / maxabs, 2^-1021)
    gopt <- gopt * scaling
    pq <- pq * scaling
    hq <- hq * scaling
    scaled <- TRUE
  }
  tolsq <- tol * tol
  s <- numeric(n)
  crvmin <- 0
  gg <- .pw_inprod(gopt, gopt)
  gg0 <- gg
  d <- -gopt
  hs <- numeric(n)
  ds <- 0
  ss <- 0
  delsq <- delta * delta
  qred <- 0
  dd <- gg
  twod <- FALSE
  iter <- 1L
  while (iter <= n) {
    if (is.nan(gg) || isTRUE(tolsq * gg0 >= gg)) break
    if (iter == 1L) {
      bstep <- delta / sqrt(dd)
    } else {
      resid <- delsq - ss
      if (isTRUE(resid <= 0)) {
        twod <- TRUE
        break
      }
      if (isTRUE(delsq * .pr_eps >= dd) || is.nan(ds)) break
      sqrtd <- max(sqrt(.pw_fma(ds, ds, resid * dd)), abs(ds), sqrt(resid * dd))
      bstep <- if (ds <= 0) (sqrtd - ds) / dd else resid / (ds + sqrtd)
    }
    if (isTRUE(bstep <= 0) || !is.finite(bstep)) break
    hd <- .pw_hess_mul(d, xpt, pq, hq)
    dhd <- .pw_inprod(d, hd)
    alpha <- bstep
    if (isTRUE(dhd > 0)) {
      alpha <- .pw_fminnm(gg / dhd, bstep)
      crvmin <- if (iter == 1L) dhd / dd else .pw_fminnm(dhd / dd, crvmin)
    }
    qadd <- .pw_fma(-(alpha * 0.5), dhd, gg) * alpha
    qred <- qred + qadd
    sold <- s
    s <- .pw_fma(d, alpha, s)
    ss <- .pw_inprod(s, s)
    hs <- .pw_fma(hd, alpha, hs)
    gh <- gopt + hs
    ggnew <- .pw_inprod(gh, gh)
    if (!is.finite(.pw_sumabs(s))) {
      s <- sold
      gg <- ggnew
      break
    }
    if (isTRUE(alpha >= bstep) || isTRUE(ss >= delsq)) {
      crvmin <- 0
      twod <- n != 1L && isTRUE(tolsq * gg0 < ggnew)
      gg <- ggnew
      break
    }
    if (isTRUE(tol * qred >= qadd)) {
      gg <- ggnew
      break
    }
    beta <- ggnew / gg
    d <- .pw_fma(d, beta, -gopt) - hs
    dd <- .pw_inprod(d, d)
    ds <- .pw_inprod(d, s)
    if (isTRUE(ds <= 0)) {
      gg <- ggnew
      break
    }
    iter <- iter + 1L
    gg <- ggnew
  }
  if (isTRUE(ss > 0) && twod) {
    for (k in seq_len(max(n - iter, 1L))) {
      if (is.nan(gg) || isTRUE(tolsq * gg0 >= gg)) break
      sg <- .pw_inprod(s, gopt)
      shs <- .pw_inprod(s, hs)
      gh <- gopt + hs
      d <- gh - .pw_project1(gh, s)
      if (isTRUE(.pw_norm(d) <= sqrt(tol) * sqrt(gg))) break
      d <- (.pw_norm(s) / .pw_norm(d)) * d
      if (isTRUE(abs(.pw_inprod(d, s)) >= (.pw_norm(d) * 0.1) * .pw_norm(s))) break
      if (isTRUE(.pw_norm(d) >= delta + delta)) break
      hd <- .pw_hess_mul(d, xpt, pq, hq)
      dg <- .pw_inprod(d, gopt)
      dhd <- .pw_inprod(hd, d)
      dhs <- .pw_inprod(hd, s)
      args <- c(sg, (shs - dhd) * 0.5, dg, dhs)
      theta <- .pw_circle_pick(.nw_circle_trsapp, args, 50L, FALSE)
      cth <- cos(theta)
      sth <- sin(theta)
      sold <- s
      s <- .pw_fma(s, cth, sth * d)
      if (!is.finite(.pw_sumabs(s))) {
        s <- sold
        break
      }
      f0 <- .pw_fma(args[3L] + args[4L], 0, args[1L] + args[2L])
      ft <- .pw_fma(cth, .pw_fma(cth, args[2L], args[1L]), .pw_fma(cth, args[4L], args[3L]) * sth)
      reduc <- f0 - ft
      qred <- qred + reduc
      if (isTRUE(reduc / qred <= tol)) break
      hs <- .pw_fma(hs, cth, sth * hd)
      gh <- hs + gopt
      gg <- .pw_inprod(gh, gh)
    }
  }
  if (is.nan(crvmin)) {
    crvmin <- 0
  } else if (scaled && crvmin > 0) {
    crvmin <- crvmin / scaling
  }
  list(d = s, crvmin = crvmin)
}

# setdrop_tr: the interpolation point to replace by a trust region step
# (0 for none)
.nw_setdrop_tr <- function(idz, kopt, ximproved, bmat, d, delta, rho, xpt, zmat) {
  x <- if (ximproved) xpt[, kopt] + d else xpt[, kopt]
  distsq <- .pw_distsq(xpt, x)
  den2 <- .pw_fmaxnm(delta * 0.1, rho)
  den2 <- den2 * den2
  w <- .pw_vfmaxnm(distsq / den2, 1)
  weight <- (w * w) * w
  score <- abs(.pw_calden(kopt, bmat, d, xpt, zmat, idz)) * weight
  if (!ximproved) score[kopt] <- -1
  if (any(score > 1, na.rm = TRUE) || (ximproved && any(score > 0, na.rm = TRUE))) return(.pw_maxloc(score))
  if (ximproved) return(.pw_maxloc(distsq))
  0L
}

# updateq: the model after XPT(:, knew) is replaced; returns list(gopt, hq, pq)
.nw_updateq <- function(idz, knew, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq) {
  hq <- .pw_r1sym(hq, pq[knew], xdrop)
  pq[knew] <- 0
  pqinc <- .pw_omega_col(idz, zmat, knew) * moderr
  pq <- pq + pqinc
  gopt <- .pw_hess_mul(xosav, xpt, pqinc) + .pw_fma(moderr, bmat[, knew], gopt)
  if (ximproved) gopt <- gopt + .pw_hess_mul(d, xpt, pq, hq)
  list(gopt = gopt, hq = hq, pq = pq)
}

# tryqalt: the alternative model of least Frobenius norm when the model has
# been poor three times in a row; returns list(itest, gopt, hq, pq)
.nw_tryqalt <- function(idz, bmat, fval, ratio, xopt, xpt, zmat, itest, gopt, hq, pq) {
  npt <- ncol(xpt)
  pqalt <- .pw_omega_mul(idz, zmat, fval)
  galt <- .pw_matprod21(bmat[, seq_len(npt), drop = FALSE], fval) + .pw_hess_mul(xopt, xpt, pqalt)
  if (isTRUE(ratio > 0.1) || isTRUE(.pw_inprod(gopt, gopt) < 10 * .pw_inprod(galt, galt))) {
    itest <- 0L
  } else {
    itest <- itest + 1L
    if (itest >= 3L) {
      gopt <- galt
      pq <- pqalt
      hq[] <- 0
      itest <- 0L
    }
  }
  list(itest = itest, gopt = gopt, hq = hq, pq = pq)
}

# geostep: a step of length delbar that improves the geometry of the
# interpolation set when XPT(:, knew) is replaced; biglag maximizes the
# Lagrange function, bigden the denominator of the update of H
.nw_circle_biglag <- function(theta, a) {
  cth <- cos(theta)
  sth <- sin(theta)
  .pw_fma(.pw_fma(a[5L], cth, a[3L]), sth, .pw_fma(.pw_fma(a[4L], cth, a[2L]), cth, a[1L]))
}
.nw_biglag <- function(idz, knew, kopt, bmat, delbar, xpt, zmat) {
  n <- nrow(xpt)
  hcol <- .pw_omega_col(idz, zmat, knew)
  xopt <- xpt[, kopt]
  d <- xpt[, knew] - xopt
  dd <- .pw_inprod(d, d)
  hd <- .pw_hess_mul(d, xpt, hcol)
  gopt <- .pw_hess_mul(xopt, xpt, hcol) + bmat[, knew]
  gg <- .pw_inprod(gopt, gopt)
  gd <- .pw_inprod(d, gopt)
  dhd <- .pw_inprod(d, hd)
  scaling <- delbar / sqrt(dd)
  if (isTRUE(gd * dhd < 0)) scaling <- -scaling
  temp <- if (isTRUE(gd * gd > (dd * 0.99) * gg)) 1 else 0
  tau <- .pw_fma(scaling * 0.5, abs(dhd), abs(gd)) * scaling
  if (isTRUE((delbar * delbar) * gg < (tau * tau) * 0.01)) temp <- 1
  if (!is.finite(.pw_sumabs(scaling * d))) return(d)
  d <- scaling * d
  hd <- scaling * hd
  s <- .pw_fma(hd, temp, gopt)
  for (iter in seq_len(n)) {
    ss <- .pw_inprod(s, s)
    s <- s - .pw_project1(s, d)
    if (isTRUE(.pw_norm(s) <= 2^-13 * sqrt(ss))) break
    s <- (.pw_norm(d) / .pw_norm(s)) * s
    if (isTRUE(abs(.pw_inprod(d, s)) >= (.pw_norm(d) * 0.1) * .pw_norm(s))) break
    if (isTRUE(.pw_norm(s) >= delbar + delbar)) break
    hs <- .pw_hess_mul(s, xpt, hcol)
    a1 <- .pw_inprod(s, hs) * 0.5
    a2 <- .pw_inprod(d, gopt)
    a3 <- .pw_inprod(s, gopt)
    a4 <- .pw_fma(.pw_inprod(d, hd), 0.5, -a1)
    a5 <- .pw_inprod(s, hd)
    angle <- .pw_circle_pick(.nw_circle_biglag, c(a1, a2, a3, a4, a5), 50L, TRUE)
    cth <- cos(angle)
    sth <- sin(angle)
    dold <- d
    d <- .pw_fma(cth, d, sth * s)
    if (!is.finite(.pw_sumabs(d))) {
      d <- dold
      break
    }
    fth <- .pw_fma(sth, .pw_fma(cth, a5, a3), .pw_fma(cth, .pw_fma(cth, a4, a2), a1))
    f0 <- .pw_fma(a3 + a5, 0, (a2 + a4) + a1)
    if (isTRUE(abs(fth) <= 1.1 * abs(f0))) break
    hd <- .pw_fma(cth, hd, sth * hs)
    s <- gopt + hd
  }
  d
}
.nw_circle_bigden <- function(theta, a) {
  t2 <- theta + theta
  t3 <- theta * 3
  t4 <- theta * 4
  v <- list(1, cos(theta), sin(theta), cos(t2), sin(t2), cos(t3), sin(t3), cos(t4), sin(t4))
  z <- 0
  for (i in 1:9) z <- .pw_fma(a[i], v[[i]], z)
  z
}
.nw_bigden <- function(idz, knew, kopt, bmat, d0, xpt, zmat) {
  n <- nrow(xpt)
  npt <- ncol(xpt)
  kk <- seq_len(npt)
  xopt <- xpt[, kopt]
  dnorm <- .pw_norm(d0)
  hcol <- .pw_omega_col(idz, zmat, knew)
  alpha <- hcol[knew]
  d <- d0
  dd <- .pw_inprod(d, d)
  s <- xpt[, knew] - xopt
  ds <- .pw_inprod(d, s)
  ss <- .pw_inprod(s, s)
  xoptsq <- .pw_inprod(xopt, xopt)
  if (!isTRUE(ds * ds <= (dd * 0.99) * ss)) {
    dtest <- (ds * ds) / ss
    xptemp <- xpt - xopt
    dstemp <- .pw_matprod12(d, xptemp)
    sstemp <- .pw_distsq(xptemp, numeric(n))
    dstemp[kopt] <- .pw_fma(ds, 2, 1)
    sstemp[kopt] <- ss
    k <- .pw_minloc((dstemp * dstemp) / sstemp)
    if (!isTRUE((dstemp[k] * dstemp[k]) / sstemp[k] >= dtest) && k != kopt) s <- xpt[, k] - xopt
  }
  densav <- 0
  for (iter in seq_len(n)) {
    ss <- .pw_inprod(s, s)
    s <- s - .pw_project1(s, d)
    if (isTRUE(.pw_norm(s) <= 2^-13 * sqrt(ss))) break
    s <- (s / .pw_norm(s)) * .pw_norm(d)
    if (isTRUE(abs(.pw_inprod(d, s)) >= (.pw_norm(d) * 0.1) * .pw_norm(s))) break
    if (isTRUE(.pw_norm(s) >= dnorm + dnorm)) break
    a <- .pw_inprod(xopt, d)
    b <- .pw_inprod(xopt, s)
    dd <- .pw_inprod(d, d)
    den <- numeric(9L)
    den[1L] <- (b * 0.5) * b + .pw_fma(a, a * 0.5, .pw_fma(dd, 0.5, xoptsq) * dd)
    den[2L] <- (a + a) * dd
    den[3L] <- (b + b) * dd
    den[4L] <- .pw_fma(a, a * 0.5, -((b * 0.5) * b))
    den[5L] <- a * b
    ta <- .pw_matprod12(d, xpt)
    tb <- .pw_matprod12(s, xpt)
    tc <- .pw_matprod12(xopt, xpt)
    W <- matrix(0, npt + n, 5L)
    W[kk, 1L] <- .pw_fma(ta, ta, tb * tb) * 0.25
    W[kk, 2L] <- ta * tc
    W[kk, 3L] <- tb * tc
    W[kk, 4L] <- .pw_fma(ta, ta, -(tb * tb)) * 0.25
    W[kk, 5L] <- (ta * 0.5) * tb
    W[npt + seq_len(n), 2L] <- d
    W[npt + seq_len(n), 3L] <- s
    P <- matrix(0, npt + n, 5L)
    for (jc in 1:5) {
      P[kk, jc] <- .pw_omega_mul(idz, zmat, W[kk, jc])
      if (jc == 2L || jc == 3L) P[kk, jc] <- P[kk, jc] + .pw_matprod12(W[npt + seq_len(n), jc], bmat[, kk, drop = FALSE])
      P[npt + seq_len(n), jc] <- .pw_matprod21(bmat, W[, jc])
    }
    P1 <- P[, 1L]; P2 <- P[, 2L]; P3 <- P[, 3L]; P4 <- P[, 4L]; P5 <- P[, 5L]
    W1 <- W[, 1L]; W2 <- W[, 2L]; W3 <- W[, 3L]; W4 <- W[, 4L]; W5 <- W[, 5L]
    # every den(j) is a chain over k: den = fma(A1, B1, den - pre), then
    # den(4) and den(8) take a second fma, then den(1) loses sm
    sm <- (P1 * 0.5) * W1
    sm <- .pw_fma(P2 * 0.5, W2, sm)
    sm <- .pw_fma(P3 * 0.5, W3, sm)
    sm <- .pw_fma(P4 * 0.5, W4, sm)
    sm <- .pw_fma(P5 * 0.5, W5, sm)
    ta2 <- .pw_fma(P1, W2, W1 * P2); tb2 <- .pw_fma(P2, W4, W2 * P4); tc2 <- .pw_fma(P3, W5, W3 * P5)
    ta3 <- .pw_fma(W3, P1, P3 * W1); tb3 <- .pw_fma(P2, W5, W2 * P5); tc3 <- .pw_fma(P3, W4, W3 * P4)
    ta4 <- .pw_fma(P1, W4, W1 * P4)
    ta5 <- .pw_fma(P1, W5, W1 * P5); tb5 <- .pw_fma(W3, P2, P3 * W2)
    ta9 <- .pw_fma(P4, W5, W4 * P5)
    zero <- numeric(npt + n)
    half <- rep(0.5, npt + n)
    A1 <- cbind(-(P1 * 0.5), -(tb2 + tc2), -(tb3 - tc3), -(P2 * 0.5), -tb5, -(tb2 - tc2), -(tb3 + tc3), -(P4 * 0.5), -ta9)
    B1 <- cbind(W1, half, half, W2, half, half, half, W4, half)
    pre <- cbind(zero, ta2, ta3, ta4, ta5, zero, zero, zero, zero)
    A2 <- cbind(P3 * 0.5, P5 * 0.5)
    B2 <- cbind(W3, W5)
    for (k in seq_len(npt + n)) {
      den <- .pw_fma(A1[k, ], B1[k, ], den - pre[k, ])
      den[c(4L, 8L)] <- .pw_fma(A2[k, ], B2[k, ], den[c(4L, 8L)])
      den[1L] <- den[1L] - sm[k]
    }
    Q <- P[knew, ]
    par <- (Q * Q) * 0.5
    psum <- ((((0 + par[1L]) + par[2L]) + par[3L]) + par[4L]) + par[5L]
    tp1 <- Q[1L] + Q[1L]
    denex <- numeric(9L)
    denex[1L] <- .pw_fma(den[1L], alpha, par[1L]) + psum
    denex[2L] <- .pw_fma(Q[3L], Q[5L], .pw_fma(Q[2L], Q[4L], .pw_fma(tp1, Q[2L], den[2L] * alpha)))
    denex[3L] <- .pw_fma(-Q[4L], Q[3L], .pw_fma(Q[2L], Q[5L], .pw_fma(tp1, Q[3L], den[3L] * alpha)))
    denex[4L] <- (.pw_fma(tp1, Q[4L], den[4L] * alpha) + par[2L]) - par[3L]
    denex[5L] <- .pw_fma(Q[2L], Q[3L], .pw_fma(tp1, Q[5L], den[5L] * alpha))
    denex[6L] <- .pw_fma(-Q[3L], Q[5L], .pw_fma(Q[2L], Q[4L], alpha * den[6L]))
    denex[7L] <- .pw_fma(Q[4L], Q[3L], .pw_fma(Q[2L], Q[5L], alpha * den[7L]))
    denex[8L] <- .pw_fma(alpha, den[8L], par[4L]) - par[5L]
    denex[9L] <- .pw_fma(alpha, den[9L], Q[4L] * Q[5L])
    angle <- .pw_circle_pick(.nw_circle_bigden, denex, 50L, TRUE)
    cth <- cos(angle)
    sth <- sin(angle)
    dold <- d
    d <- .pw_fma(cth, d, sth * s)
    if (!is.finite(.pw_sumabs(d))) {
      d <- dold
      break
    }
    if (iter > 1L) densav <- .pw_fmaxnm(.nw_circle_bigden(0, denex), densav)
    denmax <- .nw_circle_bigden(angle, denex)
    if (isTRUE(abs(denmax) <= 1.1 * abs(densav))) break
    densav <- denmax
    vlag <- .pw_matprod21(P, c(1, cth, sth, cos(angle + angle), sin(angle + angle)))
    tau <- vlag[knew]
    xnew <- d + xopt
    dxn <- .pw_inprod(xnew, d)
    xnsq <- .pw_inprod(xnew, xnew)
    vk <- vlag[kk]
    w <- .pw_matprod12(xnew, xpt) * ifelse(.pw_isvec(npt, 2), .pw_fma(-alpha, vk, tau * hcol),
                                           .pw_fma(tau, hcol, -(alpha * vk)))
    s <- .pw_fma(tau, bmat[, knew], alpha * (.pw_fma(dxn, xopt, xnsq * d) - vlag[npt + seq_len(n)]))
    s <- s + .pw_matprod21(xpt, w)
  }
  d
}
.nw_geostep <- function(idz, knew, kopt, bmat, delbar, xpt, zmat) {
  d <- .nw_biglag(idz, knew, kopt, bmat, delbar, xpt, zmat)
  alpha <- .pw_omega_col(idz, zmat, knew)[knew]
  vlag <- .pw_calvlag(kopt, bmat, d, xpt, zmat, idz)
  beta <- .pw_calbeta(kopt, bmat, d, xpt, zmat, idz)
  ab <- beta * alpha
  v2 <- vlag[knew] * vlag[knew]
  denom <- ab + v2
  if (isTRUE(v2 > 0) && !is.nan(beta) && !isTRUE(abs(ab / v2 + 1) <= 0.8)) return(d)
  dbig <- .nw_bigden(idz, knew, kopt, bmat, d, xpt, zmat)
  vlag <- .pw_calvlag(kopt, bmat, dbig, xpt, zmat, idz)
  beta <- .pw_calbeta(kopt, bmat, dbig, xpt, zmat, idz)
  dn <- .pw_fma(vlag[knew], vlag[knew], beta * alpha)
  if (isTRUE(abs(dn) >= abs(denom)) || is.nan(denom)) return(dbig)
  d
}

# newuob: the main loop; returns list(x, f, nf, info)
.nw_newuob <- function(fun, maxfun, npt, ftarget, rhobeg, rhoend, x0, iprint = 0L) {
  eta1 <- 0.1
  eta2 <- 0.7
  gamma1 <- 0.5
  gamma2 <- 2
  n <- length(x0)
  ini <- .nw_initxf(fun, maxfun, ftarget, rhobeg, x0, npt, iprint)
  kopt <- ini$kopt
  nf <- ini$nf
  fval <- ini$fval
  xbase <- ini$xbase
  xpt <- ini$xpt
  x <- xbase + xpt[, kopt]
  f <- fval[kopt]
  if (ini$info != 0L) {
    .pw_retmsg("NEWUOA", iprint, ini$info, nf, f, x)
    return(list(x = x, f = f, nf = nf, info = ini$info))
  }
  hm <- .nw_inith(ini$ij, xpt)
  idz <- hm$idz
  bmat <- hm$bmat
  zmat <- hm$zmat
  qm <- .nw_initq(ini$ij, fval, xpt)
  gopt <- qm$gopt
  hq <- qm$hq
  pq <- qm$pq
  ratio <- -1
  rho <- rhobeg
  delta <- rho
  dnormsav <- c(.pr_realmax, .pr_realmax)
  moderrsav <- c(.pr_realmax, .pr_realmax)
  itest <- 0L
  knew_tr <- 0L
  gamma3 <- .pw_fmaxnm(.pw_fminnm(0.75 * gamma2, 1.5), 1)
  maxtr <- if (maxfun > .Machine$integer.max %/% 2L) maxfun else 2L * maxfun  # 2 maxfun overflows in Fortran
  info <- 20L
  dnorm <- 0
  d <- numeric(n)
  # updates the model, the H matrix and the interpolation set when
  # XPT(:, knew) moves to XOPT + d (updateh, updatexf, updateq of newuob)
  swap_point <- function(knew, ximproved, moderr) {
    xdrop <- xpt[, knew]
    xosav <- xpt[, kopt]
    h <- .pw_updateh(knew, kopt, d, xpt, idz, bmat, zmat)
    idz <<- h$idz
    bmat <<- h$bmat
    zmat <<- h$zmat
    xpt[, knew] <<- xosav + d
    fval[knew] <<- f
    if (ximproved) kopt <<- knew
    q <- .nw_updateq(idz, knew, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
    gopt <<- q$gopt
    hq <<- q$hq
    pq <<- q$pq
  }
  for (tr in seq_len(maxtr)) {
    ts <- .nw_trsapp(delta, gopt, hq, pq, 0.01, xpt)
    d <- ts$d
    crvmin <- ts$crvmin
    dnorm <- .pw_fminnm(.pw_norm(d), delta)
    shortd <- isTRUE(0.5 * rho > dnorm)
    qinc <- .pw_quadinc(d, xpt, gopt, pq, hq)
    qred <- -qinc
    rhosq <- rho * rho
    trfail <- !isTRUE(qred > rhosq * 0x1.4f8b58p-17)
    if (shortd || trfail) {
      delta <- delta * 0.1
      if (delta <= gamma3 * rho) delta <- rho
    } else {
      x <- (xpt[, kopt] + d) + xbase
      f <- .pw_evaluatef(fun, x)
      nf <- nf + 1L
      .pw_fmsg(iprint, "Trust region", nf, delta, f, x)
      sub <- .pw_checkexit(maxfun, nf, f, ftarget, x)
      if (sub != 0L) {
        info <- sub
        break
      }
      fopt <- fval[kopt]
      moderr <- (f - fopt) - qinc
      ratio <- .pr_redrat(fopt - f, qred, eta1)
      delta <- .pr_trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
      if (delta <= gamma3 * rho) delta <- rho
      ximproved <- f < fopt
      knew_tr <- .nw_setdrop_tr(idz, kopt, ximproved, bmat, d, delta, rho, xpt, zmat)
      dnormsav <- c(dnormsav[2L], dnorm)
      moderrsav <- c(moderrsav[2L], moderr)
      if (knew_tr > 0L) {
        swap_point(knew_tr, ximproved, moderr)
        qa <- .nw_tryqalt(idz, bmat, fval - fval[kopt], ratio, xpt[, kopt], xpt, zmat, itest, gopt, hq, pq)
        itest <- qa$itest
        gopt <- qa$gopt
        hq <- qa$hq
        pq <- qa$pq
      }
    }
    accurate_mod <- all(abs(moderrsav) <= (crvmin * 0.125) * rhosq) && all(dnormsav <= rho)
    distsq <- .pw_distsq(xpt, xpt[, kopt])
    close_itpset <- all(distsq <= (delta * delta) * 4)
    adequate_geo <- isTRUE((shortd && accurate_mod) || close_itpset)
    small_trrad <- isTRUE(.pw_fmaxnm(delta, dnorm) <= rho)
    improve_geo <- (shortd || trfail || isTRUE(ratio <= eta1) || knew_tr == 0L) && !adequate_geo
    reduce_rho <- (shortd || trfail || isTRUE(ratio <= 0) || knew_tr == 0L) && adequate_geo && small_trrad
    if (improve_geo) {
      knew_geo <- .pw_maxloc(distsq)
      delbar <- .pw_fmaxnm(.pw_fminnm(sqrt(max(distsq, na.rm = TRUE)) * 0.1, delta * 0.5), rho)
      d <- .nw_geostep(idz, knew_geo, kopt, bmat, delbar, xpt, zmat)
      x <- (xpt[, kopt] + d) + xbase
      f <- .pw_evaluatef(fun, x)
      nf <- nf + 1L
      .pw_fmsg(iprint, "Geometry", nf, delbar, f, x)
      sub <- .pw_checkexit(maxfun, nf, f, ftarget, x)
      if (sub != 0L) {
        info <- sub
        break
      }
      dnorm <- .pw_fminnm(delbar, .pw_norm(d))
      fopt <- fval[kopt]
      moderr <- (f - fopt) - .pw_quadinc(d, xpt, gopt, pq, hq)
      swap_point(knew_geo, f < fopt, moderr)
      dnormsav <- c(dnormsav[2L], dnorm)
      moderrsav <- c(moderrsav[2L], moderr)
    }
    if (reduce_rho) {
      if (rho <= rhoend) {
        info <- 0L
        break
      }
      rhonew <- .pr_redrho(rho, rhoend)
      delta <- .pw_fmaxnm(0.5 * rho, rhonew)
      rho <- rhonew
      .pw_rhomsg(iprint, nf, delta, fval[kopt], rho, xbase + xpt[, kopt])
      dnormsav <- c(.pr_realmax, .pr_realmax)
      moderrsav <- c(.pr_realmax, .pr_realmax)
    }
    xopt <- xpt[, kopt]
    if (isTRUE(.pw_inprod(xopt, xopt) >= (delta * delta) * 100)) {
      sb <- .pw_shiftbase(kopt, xbase, xpt, zmat, bmat, pq, hq, idz)
      xbase <- sb$xbase
      xpt <- sb$xpt
      bmat <- sb$bmat
      hq <- sb$hq
    }
  }
  if (info == 0L && isTRUE(0.5 * rho > dnorm) && nf < maxfun) {
    x <- (xpt[, kopt] + d) + xbase
    f <- .pw_evaluatef(fun, x)
    nf <- nf + 1L
    .pw_fmsg(iprint, "Trust region", nf, delta, f, x)
  }
  if (isTRUE(fval[kopt] < f) || is.nan(f)) {
    x <- xbase + xpt[, kopt]
    f <- fval[kopt]
  }
  .pw_retmsg("NEWUOA", iprint, info, nf, f, x)
  list(x = x, f = f, nf = nf, info = info)
}

# newuoa: preprocessing and the call of newuob
.nw_newuoa <- function(fun, x0, rhobeg = 1, rhoend = 1e-6 * rhobeg, ftarget = -Inf,
                       maxfun = 500L * length(x0), npt = 2L * length(x0) + 1L, iprint = 0L) {
  n <- length(x0)
  x0 <- .pw_moderatex(as.numeric(x0))
  op <- .pw_preproc("NEWUOA", n, maxfun, ftarget, rhobeg, rhoend, npt)
  .nw_newuob(fun, op$maxfun, op$npt, op$ftarget, op$rhobeg, op$rhoend, x0, as.integer(iprint))
}
