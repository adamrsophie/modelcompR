# COBYLA of PRIMA (Zaikun Zhang's modern reference implementation of
# Powell's COBYLA): a port of the Python translation pyprima that SciPy
# 1.18.1 vendors (scipy/_external/pyprima/cobyla and common), which mirrors
# the Fortran code that PRIMA.jl calls. Two places differ between the two
# and follow the `fortran` flag:
#   * pyprima solves the least squares problems of TRSTLP with
#     numpy.linalg.lstsq, the Fortran code by back substitution with the
#     factorization kept by TRSTLP (lsqr in linalg.f90);
#   * pyprima compares 0 based constraint indices with the number of linear
#     constraints in TRSTLP, so the objective row counts as a constraint
#     there; the Fortran code excludes it.
# The Fortran variant is inferred from pyprima: no Fortran source of PRIMA is
# at hand, and libprima 0.7.1 (PRIMA_jll of PRIMA.jl) differs from the newer
# PRIMA that pyprima translates in some step, so its iterations split from
# this port after the first few steps while both reach the same points.
# Indices below are 1 based; comments name the pyprima function.

.pr_realmax <- .Machine$double.xmax
.pr_eps <- .Machine$double.eps
.pr_funcmax <- 1e30
.pr_constrmax <- 1e30
.pr_boundmax <- .Machine$double.xmax / 4

.pr_isminor <- function(x, ref) {
  refa <- abs(ref) + 0.1 * abs(x)
  refb <- abs(ref) + 0.2 * abs(x)
  abs(ref) >= refa | refa >= refb
}

.pr_hypot <- function(a, b) Mod(complex(real = a, imaginary = b))

# Givens rotation G with G %*% x = (r, 0) (planerot)
.pr_planerot <- function(x) {
  realmin <- .Machine$double.xmin
  if (any(is.nan(x))) {
    c <- 1; s <- 0
  } else if (all(is.infinite(x))) {
    c <- sign(x[1L]) / sqrt(2); s <- sign(x[2L]) / sqrt(2)
  } else if (abs(x[1L]) <= 0 && abs(x[2L]) <= 0) {
    c <- 1; s <- 0
  } else if (abs(x[2L]) <= .pr_eps * abs(x[1L])) {
    c <- sign(x[1L]); s <- 0
  } else if (abs(x[1L]) <= .pr_eps * abs(x[2L])) {
    c <- 0; s <- sign(x[2L])
  } else if (all(sqrt(realmin) < abs(x) & abs(x) < sqrt(.pr_realmax / 2.1))) {
    r <- sqrt(sum(x * x))
    c <- x[1L] / r; s <- x[2L] / r
  } else if (abs(x[1L]) > abs(x[2L])) {
    t <- x[2L] / x[1L]
    u <- max(1, abs(t), sqrt(1 + t * t)) * sign(x[1L])
    c <- 1 / u; s <- t / u
  } else {
    t <- x[1L] / x[2L]
    u <- max(1, abs(t), sqrt(1 + t * t)) * sign(x[2L])
    c <- t / u; s <- 1 / u
  }
  matrix(c(c, -s, s, c), 2L, 2L)
}

# add the column cvec to the QR factorization kept as Q and the diagonal of R
# (qradd_Rdiag); n is the number of columns before the addition
.pr_qradd <- function(cvec, Q, Rdiag, n) {
  m <- ncol(Q)
  cq <- as.numeric(cvec %*% Q)
  cqa <- as.numeric(abs(cvec) %*% abs(Q))
  cq[.pr_isminor(cq, cqa)] <- 0
  if (m - 1L >= n + 1L) {
    for (k in (m - 1L):(n + 1L)) {
      if (abs(cq[k + 1L]) > 0) {
        G <- .pr_planerot(cq[k:(k + 1L)])
        Q[, k:(k + 1L)] <- Q[, k:(k + 1L)] %*% t(G)
        cq[k] <- .pr_hypot(cq[k], cq[k + 1L])
      }
    }
  }
  if (n < m && abs(cq[n + 1L]) > .pr_eps^2 && !.pr_isminor(cq[n + 1L], cqa[n + 1L])) n <- n + 1L
  if (n >= 1L && n <= m) Rdiag[n] <- cq[n]
  list(Q = Q, Rdiag = Rdiag, n = n)
}

# move column i of A to the end of the factorization (qrexc_Rdiag)
.pr_qrexc <- function(A, Q, Rdiag, i) {
  n <- ncol(A)
  if (i < 1L || i > n) return(list(Q = Q, Rdiag = Rdiag))
  if (i <= n - 1L) {
    for (k in i:(n - 1L)) {
      G <- .pr_planerot(c(Rdiag[k + 1L], sum(Q[, k] * A[, k + 1L])))
      Q[, k:(k + 1L)] <- Q[, c(k + 1L, k)] %*% t(G)
    }
    Rdiag[i:(n - 1L)] <- vapply(i:(n - 1L), function(k) sum(Q[, k] * A[, k + 1L]), 0)
  }
  Rdiag[n] <- sum(Q[, n] * A[, i])
  list(Q = Q, Rdiag = Rdiag)
}

# least squares solution of A x = b (lsqr)
.pr_lsqr <- function(A, b, Q, Rdiag, fortran) {
  nc <- ncol(A)
  if (nc == 0L) return(numeric())
  if (!fortran) return(as.numeric(qr.coef(qr(A, LAPACK = TRUE), b)))
  x <- numeric(nc)
  y <- b
  for (i in nc:1L) {
    yq <- sum(y * Q[, i])
    yqa <- sum(abs(y) * abs(Q[, i]))
    if (.pr_isminor(yq, yqa)) {
      x[i] <- 0
    } else {
      x[i] <- yq / Rdiag[i]
      y <- y - x[i] * A[, i]
    }
  }
  x
}

# one stage of the trust region linear programming solver (trstlp_sub)
.pr_trstlp_sub <- function(iact, nact, stage, A, b, delta, d, vmultc, z, fortran) {
  num_vars <- nrow(A)
  mcon <- ncol(A)
  zdasav <- numeric(ncol(z))
  vmultd <- numeric(length(vmultc))
  zdota <- numeric(ncol(z))
  out <- function() list(iact = iact, nact = nact, d = d, vmultc = vmultc, z = z)
  sq <- seq_len
  if (stage == 1L) {
    iact <- seq_len(mcon)
    nact <- 0L
    d <- numeric(num_vars)
    cviol <- max(0, -b)
    vmultc <- cviol + b
    z <- diag(num_vars)
    if (mcon == 0L || cviol <= 0) return(out())
    if (all(is.nan(b))) return(out())
    icon <- which.max(-b)
    num_constraints <- mcon
    sdirn <- numeric(num_vars)
  } else {
    if (sum(d * d) >= delta * delta) return(out())
    iact[mcon] <- mcon
    vmultc[mcon] <- 0
    num_constraints <- mcon - 1L
    icon <- mcon
    cviol <- max(0, as.numeric(d %*% A[, sq(num_constraints), drop = FALSE]) - b[sq(num_constraints)])
  }
  # the objective row of stage 2 is excluded from the multiplier tests by the
  # Fortran code only (see the header)
  lincon <- function(idx) if (fortran) idx <= num_constraints else idx - 1L <= num_constraints
  for (k in sq(nact)) zdota[k] <- sum(z[, k] * A[, iact[k]])
  optold <- .pr_realmax
  nactold <- nact
  nfail <- 0L
  maxiter <- min(10000L, 100L * max(num_constraints, num_vars))
  for (iter in sq(maxiter)) {
    optnew <- if (stage == 1L) cviol else sum(d * A[, mcon])
    if (isTRUE(optnew < optold) || nact > nactold) {
      nactold <- nact
      nfail <- 0L
    } else {
      nfail <- nfail + 1L
    }
    optold <- min(optold, optnew)
    if (nfail == 3L) break
    if (icon > nact) {
      zdasav[sq(nact)] <- zdota[sq(nact)]
      nactsav <- nact
      r <- .pr_qradd(A[, iact[icon]], z, zdota, nact)
      z <- r$Q; zdota <- r$Rdiag; nact <- r$n
      if (nact == nactsav + 1L) {
        if (nact != icon) {
          vmultc[c(icon, nact)] <- c(vmultc[nact], 0)
          iact[c(icon, nact)] <- iact[c(nact, icon)]
        } else {
          vmultc[nact] <- 0
        }
      } else {
        vmultd[sq(nact)] <- .pr_lsqr(A[, iact[sq(nact)], drop = FALSE], A[, iact[icon]], z[, sq(nact), drop = FALSE],
                                     zdasav[sq(nact)], fortran)
        if (!any(vmultd[sq(nact)] > 0 & lincon(iact[sq(nact)]))) break
        if (nact < mcon) vmultd[(nact + 1L):mcon] <- -1
        ok <- vmultd[sq(nact)] > 0 & lincon(iact[sq(nact)])
        fracmult <- ifelse(ok, vmultc[sq(nact)] / vmultd[sq(nact)], .pr_realmax)
        frac <- min(fracmult)
        vmultc[sq(nact)] <- pmax(0, vmultc[sq(nact)] - frac * vmultd[sq(nact)])
        if (is.nan(zdota[nact]) || abs(zdota[nact]) <= .pr_eps^2) break
        vmultc[c(icon, nact)] <- c(0, frac)
        iact[c(icon, nact)] <- iact[c(nact, icon)]
      }
      if (stage == 2L && iact[nact] != mcon) {
        if (nact <= 1L) break
        r <- .pr_qrexc(A[, iact[sq(nact)], drop = FALSE], z, zdota[sq(nact)], nact - 1L)
        z <- r$Q; zdota[sq(nact)] <- r$Rdiag
        iact[c(nact - 1L, nact)] <- iact[c(nact, nact - 1L)]
        vmultc[c(nact - 1L, nact)] <- vmultc[c(nact, nact - 1L)]
      }
      if (is.nan(zdota[nact]) || abs(zdota[nact]) <= .pr_eps^2) break
      if (stage == 1L) {
        sdirn <- sdirn - ((sum(sdirn * A[, iact[nact]]) + 1) / zdota[nact]) * z[, nact]
      } else {
        sdirn <- -1 / zdota[nact] * z[, nact]
      }
    } else {
      r <- .pr_qrexc(A[, iact[sq(nact)], drop = FALSE], z, zdota[sq(nact)], icon)
      z <- r$Q; zdota[sq(nact)] <- r$Rdiag
      iact[icon:nact] <- c(iact[seq.int(icon + 1L, length.out = nact - icon)], iact[icon])
      vmultc[icon:nact] <- c(vmultc[seq.int(icon + 1L, length.out = nact - icon)], vmultc[icon])
      nact <- nact - 1L
      if (nact > 0L && (is.nan(zdota[nact]) || abs(zdota[nact]) <= .pr_eps^2)) break
      if (stage == 1L) {
        sdirn <- sdirn - sum(sdirn * z[, nact + 1L]) * z[, nact + 1L]
      } else {
        # a 0 based index of -1 wraps to the last column in pyprima
        k <- if (nact == 0L) ncol(z) else nact
        sdirn <- -1 / zdota[k] * z[, k]
      }
    }
    dd <- delta * delta - sum(d * d)
    ss <- sum(sdirn * sdirn)
    sd <- sum(sdirn * d)
    if (dd <= 0 || ss <= .pr_eps * delta * delta || is.nan(sd)) break
    sqrtd <- max(sqrt(ss * dd + sd * sd), abs(sd), sqrt(ss * dd))
    step <- if (sd > 0) dd / (sqrtd + sd) else (sqrtd - sd) / ss
    if (step <= 0 || !is.finite(step)) break
    if (stage == 1L) {
      if (.pr_isminor(cviol, step)) break
      step <- min(step, cviol)
    }
    dnew <- d + step * sdirn
    act <- iact[sq(nact)]
    if (stage == 1L) cviol <- max(0, as.numeric(dnew %*% A[, act, drop = FALSE]) - b[act])
    vmultd[sq(nact)] <- -.pr_lsqr(A[, act, drop = FALSE], dnew, z[, sq(nact), drop = FALSE], zdota[sq(nact)], fortran)
    if (stage == 2L && nact > 0L) vmultd[nact] <- max(0, vmultd[nact])
    Ai <- A[, iact, drop = FALSE]
    cvshift <- cviol - (as.numeric(dnew %*% Ai) - b[iact])
    cvsabs <- as.numeric(abs(dnew) %*% abs(Ai)) + abs(b[iact]) + cviol
    cvshift[.pr_isminor(cvshift, cvsabs)] <- 0
    if (nact < mcon) vmultd[(nact + 1L):mcon] <- cvshift[(nact + 1L):mcon]
    fracmult <- ifelse(vmultd < 0, vmultc / (vmultc - vmultd), .pr_realmax)
    cand <- c(1, fracmult)
    icon <- which.min(cand) - 1L
    frac <- min(cand)
    dold <- d
    d <- (1 - frac) * d + frac * dnew
    vmultc <- pmax(0, (1 - frac) * vmultc + frac * vmultd)
    if (!(is.finite(sum(abs(d))) && is.finite(sum(abs(vmultc))))) {
      d <- dold
      break
    }
    if (stage == 1L) cviol <- max(0, as.numeric(d %*% A) - b)
    if (icon < 1L || icon > mcon) break
  }
  out()
}

# minimize g'd subject to A'd <= b (least violation first) and |d| <= delta
# (trstlp)
.pr_trstlp <- function(A, b, delta, g, fortran) {
  nc <- ncol(A)
  num_vars <- nrow(A)
  vmultc <- numeric(nc + 1L)
  iact <- integer(nc + 1L)
  A_aug <- cbind(A, g)
  b_aug <- c(b, 0)
  for (i in seq_len(nc + 1L)) {
    maxval <- max(abs(A_aug[, i]))
    if (maxval > 1e12) {
      modscal <- max(2 * .Machine$double.xmin, 1 / maxval)
      A_aug[, i] <- A_aug[, i] * modscal
      b_aug[i] <- b_aug[i] * modscal
    }
  }
  s1 <- .pr_trstlp_sub(iact[seq_len(nc)], 0L, 1L, A_aug[, seq_len(nc), drop = FALSE], b_aug[seq_len(nc)], delta,
                       numeric(num_vars), vmultc[seq_len(nc)], matrix(0, num_vars, num_vars), fortran)
  iact[seq_len(nc)] <- s1$iact
  vmultc[seq_len(nc)] <- s1$vmultc
  s2 <- .pr_trstlp_sub(iact, s1$nact, 2L, A_aug, b_aug, delta, s1$d, vmultc, s1$z, fortran)
  s2$d
}

.pr_trrad <- function(delta_in, dnorm, eta1, eta2, gamma1, gamma2, ratio) {
  if (ratio <= eta1) gamma1 * dnorm
  else if (ratio <= eta2) max(gamma1 * delta_in, dnorm)
  else max(gamma1 * delta_in, gamma2 * dnorm)
}

.pr_redrat <- function(ared, pred, rshrink) {
  if (is.nan(ared)) -.pr_realmax
  else if (is.nan(pred) || pred <= 0) (if (ared > 0) rshrink / 2 else -.pr_realmax)
  else if (pred == Inf && ared == Inf) 1
  else if (pred == Inf && ared == -Inf) -.pr_realmax
  else ared / pred
}

.pr_redrho <- function(rho_in, rhoend) {
  rho_ratio <- rho_in / rhoend
  if (rho_ratio > 250) 0.1 * rho_in else if (rho_ratio <= 16) rhoend else sqrt(rho_ratio) * rhoend
}

# exit flags of PRIMA (common/infos.py)
.pr_info <- c(SMALL_TR_RADIUS = 0L, FTARGET_ACHIEVED = 1L, TRSUBP_FAILED = 2L, MAXFUN_REACHED = 3L,
              MAXTR_REACHED = 20L, NAN_INF_X = -1L, NAN_INF_F = -2L, NAN_INF_MODEL = -3L,
              NO_SPACE_BETWEEN_BOUNDS = 6L, DAMAGING_ROUNDING = 7L, ZERO_LINEAR_CONSTRAINT = 8L,
              CALLBACK_TERMINATE = 30L, INVALID_INPUT = 100L)

.pr_checkbreak_con <- function(maxfun, nf, cstrv, ctol, f, ftarget, x) {
  info <- 0L
  if (any(is.nan(x)) || any(is.infinite(x))) info <- -1L
  if (is.nan(f) || f == Inf || is.nan(cstrv) || cstrv == Inf) info <- -2L
  if (isTRUE(cstrv <= ctol && f <= ftarget)) info <- 1L
  if (nf >= maxfun) info <- 3L
  info
}

.pr_moderatef <- function(f) {
  f <- if (is.nan(f)) .pr_funcmax else f
  min(max(f, -.pr_realmax), .pr_funcmax)
}

.pr_moderatec <- function(c) {
  c[is.nan(c)] <- .pr_constrmax
  pmin(pmax(c, -.pr_constrmax), .pr_constrmax)
}

# f and the constraints (linear first, then nonlinear) at x (evaluate)
.pr_evaluate <- function(calcfc, x, m_nlcon, amat, bvec) {
  m_lcon <- length(bvec)
  constr <- numeric(m_lcon + m_nlcon)
  if (m_lcon) constr[seq_len(m_lcon)] <- as.numeric(amat %*% x) - bvec
  if (any(is.nan(x))) {
    f <- sum(x)
    constr <- rep(f, m_nlcon)
  } else {
    r <- calcfc(pmin(pmax(x, -.pr_realmax), .pr_realmax))
    f <- .pr_moderatef(r$f)
    if (m_nlcon) constr[m_lcon + seq_len(m_nlcon)] <- .pr_moderatec(r$c)
  }
  list(f = f, constr = constr)
}

.pr_isbetter <- function(f1, c1, f2, c2, ctol) {
  bad <- function(v) any(is.nan(v) | v == Inf)
  ib <- bad(c(f1, c1)) && !bad(c(f2, c2))
  ib <- ib || isTRUE(f1 < f2 && c1 <= c2)
  ib <- ib || isTRUE(f1 <= f2 && c1 < c2)
  cref <- 10 * max(.pr_eps, min(ctol, 1e-2 * .pr_constrmax))
  ib || isTRUE(f1 < .pr_realmax && c1 <= ctol && (is.nan(c2) || c2 > max(ctol, cref)))
}

# keep the evaluated points in a filter of nondominated pairs (savefilt)
.pr_savefilt <- function(filt, cstrv, ctol, cweight, f, x, constr) {
  nf <- filt$nfilt
  idx <- seq_len(nf)
  for (i in idx) if (.pr_isbetter(filt$ffilt[i], filt$cfilt[i], f, cstrv, ctol)) return(filt)
  if (nf && any(filt$ffilt[idx] <= f & filt$cfilt[idx] <= cstrv, na.rm = TRUE)) return(filt)
  keep <- !vapply(idx, function(i) .pr_isbetter(f, cstrv, filt$ffilt[i], filt$cfilt[i], ctol), TRUE)
  maxfilt <- length(filt$ffilt)
  if (sum(keep) == maxfilt) {
    cfs <- pmax(filt$cfilt - ctol, 0)
    if (cweight <= 0) {
      phi <- filt$ffilt
    } else if (cweight == Inf) {
      phi <- cfs
    } else {
      phi <- pmax(filt$ffilt, -.pr_realmax)
      phi[is.nan(phi)] <- -.pr_realmax
      phi[phi == Inf] <- .pr_realmax
      phi <- phi + cweight * cfs
    }
    phimax <- max(phi)
    cref <- max(cfs[which(phi >= phimax)])
    fref <- max(filt$ffilt[which(cfs >= cref)])
    cand <- which(!(filt$ffilt > fref))
    kworst <- cand[which.max(filt$cfilt[cand])]
    if (!length(kworst)) kworst <- 1L
    keep[kworst] <- FALSE
  }
  kept <- which(keep)
  nk <- length(kept)
  if (nk) {
    filt$xfilt[, seq_len(nk)] <- filt$xfilt[, kept, drop = FALSE]
    filt$ffilt[seq_len(nk)] <- filt$ffilt[kept]
    filt$cfilt[seq_len(nk)] <- filt$cfilt[kept]
    filt$confilt[, seq_len(nk)] <- filt$confilt[, kept, drop = FALSE]
  }
  filt$xfilt[, nk + 1L] <- x
  filt$ffilt[nk + 1L] <- f
  filt$cfilt[nk + 1L] <- cstrv
  filt$confilt[, nk + 1L] <- constr
  filt$nfilt <- nk + 1L
  filt
}

# the best point of a history or filter (selectx)
.pr_selectx <- function(fh, ch, cweight, ctol) {
  nh <- length(fh)
  anyw <- function(v) any(v, na.rm = TRUE)
  if (anyw(fh < .pr_funcmax & ch < .pr_constrmax)) {
    fref <- .pr_funcmax; cref <- .pr_constrmax
  } else if (anyw(fh < .pr_realmax & ch < .pr_constrmax)) {
    fref <- .pr_realmax; cref <- .pr_constrmax
  } else if (anyw(fh < .pr_funcmax & ch < .pr_realmax)) {
    fref <- .pr_funcmax; cref <- .pr_realmax
  } else {
    fref <- .pr_realmax; cref <- .pr_realmax
  }
  if (!anyw(fh < fref & ch < cref)) return(nh)
  chs <- pmax(ch - ctol, 0)
  cmin <- min(chs[which(fh < fref)])
  cref <- max(.pr_eps, 2 * cmin)
  phi <- if (cweight <= 0) fh else if (cweight == Inf) chs else pmax(fh, -.pr_realmax) + cweight * chs
  phimin <- min(phi[which(fh < fref & chs <= cref)])
  cref <- min(chs[which(fh < fref & phi <= phimin)])
  fref <- min(fh[which(chs <= cref)])
  cand <- which(!(fh > fref))
  cand[which.min(ch[cand])]
}

# ratio of the change in f to the change in the constraints (fcratio)
.pr_fcratio <- function(conmat, fval) {
  if (!nrow(conmat)) return(0)
  cmin <- apply(-conmat, 1L, min)
  cmax <- apply(-conmat, 1L, max)
  sel <- cmin < 0.5 * cmax
  if (any(sel) && min(fval) < max(fval)) (max(fval) - min(fval)) / min(pmax(cmax, 0)[sel] - cmin[sel]) else 0
}

.pr_findpole <- function(cpen, cval, fval) {
  jopt <- length(fval)
  phi <- fval + cpen * cval
  phimin <- min(phi)
  if (phimin < phi[jopt] || any(cval < cval[jopt] & phi <= phi[jopt])) {
    cand <- which(!(phi > phimin))
    jopt <- cand[which.min(cval[cand])]
  }
  jopt
}

.pr_inv <- function(M) tryCatch(solve(M), error = function(e) matrix(NaN, nrow(M), ncol(M)))

.pr_erri <- function(simi, S) max(abs(simi %*% S - diag(nrow(S))))

# make the best vertex the pole of the simplex (updatepole)
.pr_updatepole <- function(cpen, conmat, cval, fval, sim, simi) {
  n <- nrow(sim)
  info <- 0L
  jopt <- .pr_findpole(cpen, cval, fval)
  sim_old <- sim
  simi_old <- simi
  if (jopt <= n) {
    sim[, n + 1L] <- sim[, n + 1L] + sim[, jopt]
    sj <- sim[, jopt]
    sim[, jopt] <- 0
    sim[, seq_len(n)] <- sim[, seq_len(n)] - sj
    simi[jopt, ] <- -colSums(simi)
  }
  S <- sim[, seq_len(n), drop = FALSE]
  erri <- .pr_erri(simi, S)
  if (!isTRUE(erri <= 0.1)) {
    simi_test <- .pr_inv(S)
    erri_test <- .pr_erri(simi_test, S)
    if (isTRUE(erri_test < erri) || (is.nan(erri) && !is.nan(erri_test))) {
      simi <- simi_test
      erri <- erri_test
    }
  }
  if (isTRUE(erri <= 1)) {
    if (jopt <= n) {
      fval[c(jopt, n + 1L)] <- fval[c(n + 1L, jopt)]
      conmat[, c(jopt, n + 1L)] <- conmat[, c(n + 1L, jopt)]
      cval[c(jopt, n + 1L)] <- cval[c(n + 1L, jopt)]
    }
  } else {
    info <- 7L
    sim <- sim_old
    simi <- simi_old
  }
  list(conmat = conmat, cval = cval, fval = fval, sim = sim, simi = simi, info = info)
}

# replace vertex jdrop by the new point (updatexfc)
.pr_updatexfc <- function(jdrop, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi) {
  n <- nrow(sim)
  if (is.null(jdrop))
    return(list(conmat = conmat, cval = cval, fval = fval, sim = sim, simi = simi, info = 0L))
  if (jdrop <= n) {
    sim[, jdrop] <- d
    sj <- simi[jdrop, ] / sum(simi[jdrop, ] * d)
    simi <- simi - outer(as.numeric(simi %*% d), sj)
    simi[jdrop, ] <- sj
  } else {
    sim[, n + 1L] <- sim[, n + 1L] + d
    sim[, seq_len(n)] <- sim[, seq_len(n)] - d
    simid <- as.numeric(simi %*% d)
    simi <- simi + outer(simid, colSums(simi) / (1 - sum(simid)))
  }
  # pyprima keeps no copies here: a failed update leaves the modified arrays
  simi_mod <- simi
  S <- sim[, seq_len(n), drop = FALSE]
  erri <- .pr_erri(simi, S)
  if (!isTRUE(erri <= 0.1)) {
    simi_test <- .pr_inv(S)
    erri_test <- .pr_erri(simi_test, S)
    if (isTRUE(erri_test < erri) || (is.nan(erri) && !is.nan(erri_test))) {
      simi <- simi_test
      erri <- erri_test
    }
  }
  if (isTRUE(erri <= 1)) {
    fval[jdrop] <- f
    conmat[, jdrop] <- constr
    cval[jdrop] <- cstrv
    return(.pr_updatepole(cpen, conmat, cval, fval, sim, simi))
  }
  list(conmat = conmat, cval = cval, fval = fval, sim = sim, simi = simi_mod, info = 7L)
}

# vertex to drop after a trust region step (setdrop_tr)
.pr_setdrop_tr <- function(ximproved, d, delta, rho, sim, simi) {
  n <- nrow(sim)
  S <- sim[, seq_len(n), drop = FALSE]
  if (ximproved) {
    distsq <- c(colSums((S - d)^2), sum(d * d))
  } else {
    distsq <- c(colSums(S^2), 0)
  }
  weight <- pmax(1, distsq / max(rho, delta / 10)^2)
  simid <- as.numeric(simi %*% d)
  score <- weight * abs(c(simid, 1 - sum(simid)))
  if (!ximproved) score[n + 1L] <- -1
  score[is.nan(score)] <- -1
  jdrop <- NULL
  if (any(score > 0)) jdrop <- which.max(score)
  if (ximproved && is.null(jdrop)) jdrop <- which.max(distsq)
  jdrop
}

# gradients of the linear models of the constraints (columns)
.pr_conA <- function(amat, conmat, simi, m_lcon) {
  n <- nrow(simi)
  m <- nrow(conmat)
  A <- matrix(0, n, m)
  if (m_lcon) A[, seq_len(m_lcon)] <- t(amat)
  if (m > m_lcon) {
    idx <- (m_lcon + 1L):m
    A[, idx] <- t((conmat[idx, seq_len(n), drop = FALSE] - conmat[idx, n + 1L]) %*% simi)
  }
  A
}

# geometry step that keeps the simplex well conditioned (geostep)
.pr_geostep <- function(jdrop, amat, bvec, conmat, cpen, cval, delbar, fval, simi) {
  n <- nrow(simi)
  d <- simi[jdrop, ]
  d <- delbar * (d / sqrt(sum(d * d)))
  g <- as.numeric((fval[seq_len(n)] - fval[n + 1L]) %*% simi)
  A <- .pr_conA(amat, conmat, simi, length(bvec))
  dA <- as.numeric(d %*% A)
  cvpd <- max(0, conmat[, n + 1L] + dA)
  cvnd <- max(0, conmat[, n + 1L] - dA)
  if (-sum(d * g) + cpen * cvnd < sum(d * g) + cpen * cvpd) d <- -d
  d
}

# penalty parameter that makes the predicted merit reduction positive (getcpen)
.pr_getcpen <- function(amat, bvec, conmat, cpen, cval, delta, fval, rho, sim, simi, fortran) {
  n <- nrow(sim)
  for (iter in 0:n) {
    r <- .pr_updatepole(cpen, conmat, cval, fval, sim, simi)
    conmat <- r$conmat; cval <- r$cval; fval <- r$fval; sim <- r$sim; simi <- r$simi
    if (r$info == 7L) break
    g <- as.numeric((fval[seq_len(n)] - fval[n + 1L]) %*% simi)
    A <- .pr_conA(amat, conmat, simi, length(bvec))
    d <- .pr_trstlp(A, -conmat[, n + 1L], delta, g, fortran)
    preref <- -sum(d * g)
    prerec <- cval[n + 1L] - max(0, conmat[, n + 1L] + as.numeric(d %*% A))
    if (!isTRUE(prerec > 0 && preref < 0)) break
    cpen <- max(cpen, min(-2 * preref / prerec, .pr_realmax))
    if (.pr_findpole(cpen, cval, fval) == n + 1L) break
  }
  cpen
}

# the main loop (cobylb, with initxfc and initfilt)
.pr_cobylb <- function(calcfc, maxfilt, maxfun, amat, bvec, ctol, cweight, eta1, eta2, ftarget, gamma1, gamma2,
                       rhobeg, rhoend, constr, f, x, fortran) {
  n <- length(x)
  m <- length(constr)
  m_lcon <- length(bvec)
  m_nlcon <- m - m_lcon
  np1 <- n + 1L
  cpenmin <- .pr_eps
  ev <- function(x) .pr_evaluate(calcfc, x, m_nlcon, amat, bvec)
  # initxfc
  sim <- cbind(diag(rhobeg, n), x)
  dimnames(sim) <- NULL
  simi <- diag(1 / rhobeg, n)
  evaluated <- logical(np1)
  fval <- rep(.pr_realmax, np1)
  cval <- rep(.pr_realmax, np1)
  conmat <- matrix(.pr_realmax, m, np1)
  info <- 0L
  for (k in 0:n) {
    xk <- sim[, np1]
    if (k == 0L) {
      j <- np1
      fk <- f
      ck <- constr
    } else {
      j <- k
      xk[j] <- xk[j] + rhobeg
      e <- ev(xk)
      fk <- e$f
      ck <- e$constr
    }
    cstrv <- max(0, ck)
    evaluated[j] <- TRUE
    fval[j] <- fk
    conmat[, j] <- ck
    cval[j] <- cstrv
    subinfo <- .pr_checkbreak_con(maxfun, k, cstrv, ctol, fk, ftarget, xk)
    if (subinfo != 0L) {
      info <- subinfo
      break
    }
    if (j < np1 && fval[j] < fval[np1]) {
      fval[c(j, np1)] <- fval[c(np1, j)]
      cval[c(j, np1)] <- cval[c(np1, j)]
      conmat[, c(j, np1)] <- conmat[, c(np1, j)]
      sim[, np1] <- xk
      sim[j, seq_len(j)] <- -rhobeg
    }
  }
  nf <- sum(evaluated)
  if (all(evaluated)) simi <- .pr_inv(sim[, seq_len(n), drop = FALSE])
  # initfilt
  nfl <- min(max(maxfilt, 1L), maxfun)
  filt <- list(nfilt = 0L, cfilt = numeric(nfl), ffilt = numeric(nfl), xfilt = matrix(0, n, nfl),
               confilt = matrix(0, m, nfl))
  for (i in seq_len(np1)) {
    if (evaluated[i]) {
      xi <- if (i < np1) sim[, i] + sim[, np1] else sim[, i]
      filt <- .pr_savefilt(filt, cval[i], ctol, cweight, fval[i], xi, conmat[, i])
    }
  }
  finish <- function(info, cweight) {
    kf <- seq_len(filt$nfilt)
    kopt <- .pr_selectx(filt$ffilt[kf], filt$cfilt[kf], cweight, ctol)
    list(x = filt$xfilt[, kopt], f = filt$ffilt[kopt], constr = filt$confilt[, kopt], cstrv = filt$cfilt[kopt],
         nf = nf, info = info)
  }
  if (info != 0L) return(finish(info, cweight))
  rho <- rhobeg
  delta <- rhobeg
  cpen <- max(cpenmin, min(1e3, .pr_fcratio(conmat, fval)))
  shortd <- FALSE
  ratio <- -1
  jdrop_tr <- 1L
  gamma3 <- max(1, min(0.75 * gamma2, 1.5))
  maxtr <- 10 * maxfun
  info <- 20L
  d <- numeric(n)
  setpole <- function(r) {
    conmat <<- r$conmat; cval <<- r$cval; fval <<- r$fval; sim <<- r$sim; simi <<- r$simi
    r$info
  }
  # evaluate at x unless a vertex is that close (shared by both steps)
  newpoint <- function(x) {
    distsq <- c(colSums((x - (sim[, np1] + sim[, seq_len(n), drop = FALSE]))^2), sum((x - sim[, np1])^2))
    j <- which.min(distsq)
    if (distsq[j] <= (1e-4 * rhoend)^2) return(list(f = fval[j], constr = conmat[, j], cstrv = cval[j]))
    e <- ev(x)
    cs <- max(0, e$constr)
    nf <<- nf + 1L
    filt <<- .pr_savefilt(filt, cs, ctol, cweight, e$f, x, e$constr)
    list(f = e$f, constr = e$constr, cstrv = cs)
  }
  for (tr in seq_len(maxtr)) {
    cpen <- .pr_getcpen(amat, bvec, conmat, cpen, cval, delta, fval, rho, sim, simi, fortran)
    if (setpole(.pr_updatepole(cpen, conmat, cval, fval, sim, simi)) == 7L) {
      info <- 7L
      break
    }
    adequate_geo <- all(colSums(sim[, seq_len(n), drop = FALSE]^2) <= 4 * delta^2)
    g <- as.numeric((fval[seq_len(n)] - fval[np1]) %*% simi)
    A <- .pr_conA(amat, conmat, simi, m_lcon)
    d <- .pr_trstlp(A, -conmat[, np1], delta, g, fortran)
    dnorm <- min(delta, sqrt(sum(d * d)))
    shortd <- dnorm <= 0.1 * rho
    preref <- -sum(d * g)
    prerec <- cval[np1] - max(0, conmat[, np1] + as.numeric(d %*% A))
    prerem <- preref + cpen * prerec
    trfail <- !isTRUE(prerem > 1e-6 * min(cpen, 1) * rho)
    if (shortd || trfail) {
      delta <- delta * 0.1
      if (delta <= gamma3 * rho) delta <- rho
    } else {
      x <- sim[, np1] + d
      p <- newpoint(x)
      actrem <- (fval[np1] + cpen * cval[np1]) - (p$f + cpen * p$cstrv)
      ratio <- .pr_redrat(actrem, prerem, eta1)
      delta <- .pr_trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
      if (delta <= gamma3 * rho) delta <- rho
      ximproved <- isTRUE(actrem > 0)
      jdrop_tr <- .pr_setdrop_tr(ximproved, d, delta, rho, sim, simi)
      if (setpole(.pr_updatexfc(jdrop_tr, p$constr, cpen, p$cstrv, d, p$f, conmat, cval, fval, sim, simi)) == 7L) {
        info <- 7L
        break
      }
      subinfo <- .pr_checkbreak_con(maxfun, nf, p$cstrv, ctol, p$f, ftarget, x)
      if (subinfo != 0L) {
        info <- subinfo
        break
      }
    }
    bad_trstep <- shortd || trfail || ratio <= 0 || is.null(jdrop_tr)
    improve_geo <- bad_trstep && !adequate_geo
    reduce_rho <- bad_trstep && adequate_geo && max(delta, dnorm) <= rho
    colsq <- colSums(sim[, seq_len(n), drop = FALSE]^2)
    if (improve_geo && !all(colsq <= 4 * delta^2)) {
      jdrop_geo <- which.max(colsq)
      d <- .pr_geostep(jdrop_geo, amat, bvec, conmat, cpen, cval, delta / 2, fval, simi)
      x <- sim[, np1] + d
      p <- newpoint(x)
      if (setpole(.pr_updatexfc(jdrop_geo, p$constr, cpen, p$cstrv, d, p$f, conmat, cval, fval, sim, simi)) == 7L) {
        info <- 7L
        break
      }
      subinfo <- .pr_checkbreak_con(maxfun, nf, p$cstrv, ctol, p$f, ftarget, x)
      if (subinfo != 0L) {
        info <- subinfo
        break
      }
    }
    if (reduce_rho) {
      if (rho <= rhoend) {
        info <- 0L
        break
      }
      delta <- max(0.5 * rho, .pr_redrho(rho, rhoend))
      rho <- .pr_redrho(rho, rhoend)
      cpen <- max(cpenmin, min(cpen, .pr_fcratio(conmat, fval)))
      if (setpole(.pr_updatepole(cpen, conmat, cval, fval, sim, simi)) == 7L) {
        info <- 7L
        break
      }
    }
  }
  x <- sim[, np1] + d
  if (info == 0L && shortd && sqrt(sum((x - sim[, np1])^2)) > 1e-3 * rhoend && nf < maxfun) {
    e <- ev(x)
    nf <- nf + 1L
    filt <- .pr_savefilt(filt, max(0, e$constr), ctol, cweight, e$f, x, e$constr)
  }
  finish(info, max(cpen, cweight))
}

# cobyla(): constraints Aineq x <= bineq, Aeq x = beq, xl <= x <= xu and
# nlcon(x) <= 0, where calcfc(x) returns list(f, c = nlcon(x))
.pr_cobyla <- function(calcfc, m_nlcon, x, Aineq = NULL, bineq = NULL, Aeq = NULL, beq = NULL, xl = NULL, xu = NULL,
                       f0 = NULL, nlconstr0 = NULL, rhobeg = NULL, rhoend = NULL, ftarget = -.pr_realmax,
                       ctol = sqrt(.pr_eps), cweight = 1e8, maxfun = NULL, maxfilt = 2000L, eta1 = NULL,
                       eta2 = NULL, gamma1 = 0.5, gamma2 = 2, fortran = FALSE) {
  n <- length(x)
  if (!is.null(xl)) xl <- pmax(ifelse(is.nan(xl), -.pr_boundmax, xl), -.pr_boundmax)
  if (!is.null(xu)) xu <- pmin(ifelse(is.nan(xu), .pr_boundmax, xu), .pr_boundmax)
  # get_lincon: bounds and linear equalities become inequalities
  I <- diag(n)
  rows <- list(); rhs <- list()
  if (!is.null(xl)) { k <- which(xl > -.pr_boundmax); rows$l <- -I[k, , drop = FALSE]; rhs$l <- -xl[k] }
  if (!is.null(xu)) { k <- which(xu < .pr_boundmax); rows$u <- I[k, , drop = FALSE]; rhs$u <- xu[k] }
  if (!is.null(Aeq)) { rows$e1 <- -Aeq; rhs$e1 <- -beq; rows$e2 <- Aeq; rhs$e2 <- beq }
  if (!is.null(Aineq)) { rows$i <- Aineq; rhs$i <- bineq }
  amat <- do.call(rbind, c(list(matrix(0, 0L, n)), rows))
  bvec <- as.numeric(unlist(rhs))
  m_lcon <- length(bvec)
  mmm <- m_lcon + m_nlcon
  constr <- numeric(mmm)
  if (!is.null(f0) && !is.null(nlconstr0) && all(is.finite(x))) {
    f <- .pr_moderatef(f0)
    if (m_lcon) constr[seq_len(m_lcon)] <- .pr_moderatec(as.numeric(amat %*% x) - bvec)
    if (m_nlcon) constr[m_lcon + seq_len(m_nlcon)] <- .pr_moderatec(nlconstr0)
  } else {
    x[is.nan(x)] <- 0
    x <- pmin(pmax(x, -.pr_realmax), .pr_realmax)
    e <- .pr_evaluate(calcfc, x, m_nlcon, amat, bvec)
    f <- e$f
    constr <- e$constr
    if (m_lcon) constr[seq_len(m_lcon)] <- .pr_moderatec(constr[seq_len(m_lcon)])
  }
  # defaults and preproc
  if (is.null(rhobeg)) rhobeg <- if (!is.null(rhoend) && is.finite(rhoend) && rhoend > 0) max(10 * rhoend, 1) else 1
  if (is.null(rhoend)) rhoend <- if (rhobeg > 0) max(.pr_eps, min(1e-6 * rhobeg, 1e-6)) else 1e-6
  if (is.null(maxfun)) maxfun <- 500L * n
  if (is.null(eta1)) eta1 <- if (!is.null(eta2) && eta2 > 0 && eta2 < 1) max(.pr_eps, eta2 / 7) else 0.1
  if (is.null(eta2)) eta2 <- if (eta1 > 0 && eta1 < 1) (eta1 + 2) / 3 else 0.7
  maxfun <- max(maxfun, n + 2L)
  if (is.nan(ftarget)) ftarget <- -.pr_realmax
  maxfilt <- if (maxfilt < 1) 2000L else max(200L, maxfilt)
  maxfilt <- min(maxfun, max(200L, maxfilt))
  if (abs(eta1 - eta2) < 1e2 * .pr_eps * max(abs(eta1), 1)) eta2 <- eta1
  if (is.nan(gamma1) || gamma1 <= 0 || gamma1 >= 1) gamma1 <- 0.5
  if (is.nan(gamma2) || gamma2 < 1 || is.infinite(gamma2)) gamma2 <- 2
  if (abs(rhobeg - rhoend) < 1e2 * .pr_eps * max(abs(rhobeg), 1)) rhoend <- rhobeg
  if (rhobeg <= 0 || !is.finite(rhobeg)) rhobeg <- if (is.finite(rhoend) && rhoend > 0) max(10 * rhoend, 1) else 1
  if (rhoend <= 0 || rhobeg < rhoend || !is.finite(rhoend)) rhoend <- max(.pr_eps, min(0.1 * rhobeg, 1e-6))
  if (is.nan(ctol) || ctol < 0) ctol <- sqrt(.pr_eps)
  if (is.nan(cweight) || cweight < 0) cweight <- 1e8
  .pr_cobylb(calcfc, as.integer(maxfilt), maxfun, if (m_lcon) amat, if (m_lcon) bvec, ctol, cweight, eta1, eta2,
             ftarget, gamma1, gamma2, rhobeg, rhoend, constr, f, x, fortran)
}

.pr_cobyla_message <- function(info) {
  switch(as.character(info),
         "0" = "The lower bound for the trust region radius is reached",
         "1" = "The target function value is reached",
         "2" = "A trust region step failed to reduce the quadratic model",
         "3" = "The maximum number of function evaluations is reached",
         "20" = "The maximum number of trust region iterations is reached",
         "-1" = "NaN or Inf occurs in x",
         "-2" = "The objective or constraint function returns NaN or nearly infinite values",
         "-3" = "NaN or Inf occurs in the models",
         "6" = "There is no space between the lower and upper bounds",
         "7" = "Rounding errors are becoming damaging",
         "8" = "One of the linear constraints has a zero gradient",
         "30" = "The callback requested termination",
         paste("Return code", info))
}

# minimize(method = "COBYLA") of SciPy (_minimize_cobyla with the pyprima
# interface pyprima.minimize)
.sp_cobyla <- function(fun, x0, rhobeg = 1.0, tol = 1e-4, maxiter = 1000L, disp = 0L, catol = NULL,
                       f_target = -Inf, callback = NULL, constraints = list(), bounds = NULL) {
  ctol <- catol %||% sqrt(.pr_eps)
  n0 <- length(x0)
  lb <- if (is.null(bounds)) rep(-Inf, n0) else rep_len(bounds$lb, n0)
  ub <- if (is.null(bounds)) rep(Inf, n0) else rep_len(bounds$ub, n0)
  # nonlinear constraints lb <= c(x) <= ub as c(x) <= 0 (process_nl_constraints)
  nlc <- if (length(constraints)) function(x) {
    unlist(lapply(constraints, function(cn) {
      v <- as.numeric(cn$fun(x))
      l <- rep_len(cn$lb, length(v)); u <- rep_len(cn$ub, length(v))
      c(l[l > -Inf] - v[l > -Inf], v[u < Inf] - u[u < Inf])
    }))
  }
  # variables fixed by the bounds are removed; pyprima removes them from the
  # objective only and would pass the shortened vector to the constraints,
  # here the constraints see the full vector too
  finite_max <- function(v) max(c(1, abs(v[is.finite(v)])))
  tolb <- 10 * .pr_eps * max(n0, 1) * max(finite_max(lb), finite_max(ub))
  fixed <- lb <= ub & abs(lb - ub) < tolb
  full <- function(x) x
  if (any(fixed)) {
    fv <- pmin(pmax(0.5 * (lb[fixed] + ub[fixed]), lb[fixed]), ub[fixed])
    full <- function(x) {
      z <- numeric(n0)
      z[fixed] <- fv
      z[!fixed] <- x
      z
    }
    x0 <- x0[!fixed]; lb <- lb[!fixed]; ub <- ub[!fixed]
  }
  if (is.null(nlc)) x0 <- pmin(pmax(x0, lb, na.rm = TRUE), ub, na.rm = TRUE)
  calcfc <- function(x) {
    z <- full(x)
    list(f = as.numeric(fun(z)), c = if (!is.null(nlc)) nlc(z) else numeric())
  }
  r0 <- calcfc(x0)
  res <- .pr_cobyla(calcfc, length(r0$c), x0, xl = lb, xu = ub, f0 = r0$f, nlconstr0 = r0$c, rhobeg = rhobeg,
                    rhoend = tol, maxfun = as.integer(maxiter), ctol = ctol, ftarget = f_target)
  success <- res$cstrv <= ctol && res$info %in% c(0L, 1L)
  msg <- if (res$cstrv > ctol) {
    "Did not converge to a solution satisfying the constraints. See `maxcv` for the magnitude of the violation."
  } else {
    .pr_cobyla_message(res$info)
  }
  list(x = full(res$x), status = res$info, success = success, message = msg, nfev = res$nf, fun = res$f,
       maxcv = res$cstrv)
}
