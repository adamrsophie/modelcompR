# UOBYQA of PRIMA (Zaikun Zhang's modern reference implementation of
# Powell's UOBYQA), reconstructed from the compiled library that PRIMA.jl
# calls and checked against it bit for bit (see R/optimization-prima-powell.R
# for the method). Function names follow the Fortran routines of
# src/uobyqa (uobyqb, initxf, initl, initq, trstep, geostep, setdrop_tr,
# update) and of the linalg and powalg modules (hessenberg_hhd_trid,
# eigmin_sym_trid, smat_mul_vec, vec2smat, quadinc_ghv, calvlag_qint,
# shiftbase_qint). A quadratic model and the Lagrange functions are kept as
# a gradient followed by the upper triangle of the Hessian packed by
# columns, (i, j) with i <= j at n + j (j - 1) / 2 + i.

# fma chain z + x[1] y[1] + x[2] y[2] + ... in this order
.uo_chain <- function(x, y, z) {
  for (i in seq_along(x)) z <- .pw_fma(x[i], y[i], z)
  z
}
# maxval and minval of gfortran: NaN skipped, NaN when all are NaN
.uo_maxval <- function(v) {
  v <- v[!is.na(v)]
  if (length(v)) max(v) else NaN
}
.uo_minval <- function(v) {
  v <- v[!is.na(v)]
  if (length(v)) min(v) else NaN
}
# sign(a, b) of Fortran (the sign bit of b, so that -0 counts as negative)
.uo_copysign <- function(a, b) if (isTRUE(b < 0) || (isTRUE(b == 0) && 1 / b < 0)) -abs(a) else abs(a)

.uo_vec2smat <- function(v, n) {
  S <- matrix(0, n, n)
  S[upper.tri(S, diag = TRUE)] <- v
  S[lower.tri(S)] <- t(S)[lower.tri(S)]
  S
}
.uo_smat_mul_vec <- function(v, x) {
  n <- length(x)
  y <- numeric(n)
  for (j in seq_len(n)) {
    o <- (j * (j - 1L)) %/% 2L
    y[j] <- .uo_chain(x[seq_len(j)], v[o + seq_len(j)], 0)
    if (j > 1L) {
      i <- seq_len(j - 1L)
      y[i] <- .pw_fma(v[o + i], x[j], y[i])
    }
  }
  y
}
# the vector w with Q(x + d) - Q(x) = <w, ghv> for Q of packed form ghv
.uo_wvec <- function(d, x, m) {
  n <- length(x)
  s <- x + d
  w <- numeric(m)
  w[seq_len(n)] <- d
  for (j in seq_len(n)) {
    o <- n + (j * (j - 1L)) %/% 2L
    i <- seq_len(j)
    w[o + i] <- .pw_fma(d[i], s[j], d[j] * x[i])
    w[o + j] <- w[o + j] * 0.5
  }
  w
}
.uo_quadinc <- function(ghv, d, x) .pw_inprod(.uo_wvec(d, x, length(ghv)), ghv)
# values at xref + d of the Lagrange functions (calvlag_qint)
.uo_calvlag <- function(pl, d, xref, kref) {
  v <- .pw_matprod12(.uo_wvec(d, xref, nrow(pl)), pl)
  v[kref] <- v[kref] + 1
  v
}

# Householder tridiagonalization A = Q T Q'; returns the diagonal td, the
# subdiagonal tn and the Householder vectors in the lower part of A
.uo_hessenberg <- function(A) {
  n <- nrow(A)
  smax <- .uo_maxval(abs(A))
  scaled <- isTRUE(smax > 0) && (smax > 1e8 || smax < 9.999999747378752e-05)
  if (scaled) A <- A / smax
  td <- diag(A)
  tn <- numeric(max(n - 1L, 0L))
  w <- numeric(n)
  z <- numeric(n)
  for (k in seq_len(n - 1L)) {
    lo <- if (k + 2L <= n) (k + 2L):n else integer()
    ksq <- .uo_chain(A[lo, k], A[lo, k], 0)
    if (isTRUE(ksq > 0)) {
      temp <- A[k + 1L, k]
      tn[k] <- .uo_copysign(sqrt(.pw_fma(temp, temp, ksq)), temp)
      r <- ksq / (temp + tn[k])
      A[k + 1L, k] <- -r
      s <- sqrt(2 / .pw_fma(r, r, ksq))
      idx <- (k + 1L):n
      w[idx] <- A[idx, k] * s
      A[idx, k] <- w[idx]
      z[idx] <- td[idx] * w[idx]
      for (j in seq_len(n - 1L - k) + k) {
        ii <- (j + 1L):n
        z[ii] <- .pw_fma(A[ii, j], w[j], z[ii])
        z[j] <- .uo_chain(A[ii, j], w[ii], z[j])
      }
      wz <- .pw_inprod(w[idx], z[idx])
      lane <- .pw_isvec(n - k, 2)
      a <- ifelse(lane, w[idx] * wz - 2 * z[idx], .pw_fma(wz, w[idx], -2 * z[idx]))
      td[idx] <- .pw_fma(w[idx], a, td[idx])
      for (j in seq_len(n - 1L - k) + k) {
        ii <- (j + 1L):n
        a <- .pw_fma(-w[ii], z[j], A[ii, j])
        b <- .pw_fma(-w[ii], wz, z[ii])
        A[ii, j] <- .pw_fma(-b, w[j], a)
      }
    } else {
      tn[k] <- A[k + 1L, k]
      A[k + 1L, k] <- 0
    }
  }
  if (scaled) {
    td <- td * smax
    tn <- tn * smax
  }
  list(A = A, td = td, tn = tn)
}

# the least eigenvalue of the tridiagonal matrix [td, tn] by bisection
.uo_pivots <- function(td, tn, shift) {
  n <- length(td)
  piv <- rep(-1, n)
  piv[1L] <- td[1L] - shift
  for (k in seq_len(n - 1L)) {
    if (!isTRUE(piv[k] > 0)) break
    piv[k + 1L] <- (td[k + 1L] - shift) - (tn[k] * tn[k]) / piv[k]
  }
  piv
}
.uo_eigmin <- function(td, tn, tol) {
  n <- length(td)
  piv <- rep(-1, n)
  piv[1L] <- td[1L]
  for (k in seq_len(n - 1L)) {
    if (!isTRUE(piv[k] > 0)) break
    piv[k + 1L] <- td[k + 1L] - (tn[k] * tn[k]) / piv[k]
  }
  if (isTRUE(all(piv >= 0))) {
    eigu <- .uo_minval(piv)
    eigl <- 0
  } else {
    eigu <- .uo_minval(td)
    eigl <- -.uo_maxval((abs(td) + abs(c(0, tn))) + abs(c(tn, 0)))
  }
  pivsav <- numeric(n)
  ksave <- 0L
  pivksv <- 0
  for (iter in 1:100) {
    if (isTRUE((eigu - eigl) <= .pw_fmaxnm(abs(eigl), abs(eigu)) * tol)) break
    eig <- (eigu + eigl) * 0.5
    piv <- .uo_pivots(td, tn, eig)
    if (isTRUE(all(piv > 0))) {
      pivsav <- piv
      eigl <- eig
    } else {
      k <- which(!(piv > 0) | is.na(piv))[1L]
      pivsav[seq_len(k - 1L)] <- piv[seq_len(k - 1L)]
      if (ksave == k && isTRUE(pivksv < 0) && isTRUE((pivsav[k] - piv[k]) >= (piv[k] - pivksv))) {
        eigu <- .pw_fma(pivsav[k], eig, -(piv[k] * eigl)) / (pivsav[k] - piv[k])
        pivksv <- 0
      } else {
        eigu <- eig
        ksave <- k
        pivksv <- piv[k]
      }
    }
  }
  eigl
}

# trstep: the trust region subproblem min g'd + d'hd/2, |d| <= delta, by
# the tridiagonal method of More and Sorensen as in Powell's TRSTEP;
# returns list(d, crvmin)
.uo_trstep <- function(delta, g, h, tol) {
  n <- length(g)
  d <- numeric(n)
  crvmin <- 0
  gmax <- .uo_maxval(abs(g))
  scaled <- isTRUE(gmax > 1e8)
  scaling <- 1
  gg <- g
  h1 <- h[1L, 1L]
  if (scaled) {
    scaling <- .pw_fmaxnm(1 / gmax, 2 * .Machine$double.xmin)
    gg <- g * scaling
    h <- h * scaling
  }
  gsq <- .pw_inprod(gg, gg)
  if (is.na(gsq)) return(list(d = d, crvmin = crvmin))
  if (!isTRUE(any(abs(h) > 0))) {
    gnorm <- sqrt(gsq)
    if (gnorm > 0) d <- -(gg * (delta / gnorm))
    return(list(d = d, crvmin = crvmin))
  }
  if (n == 1L) {
    # the one dimensional case works on the unscaled g and h
    d <- .uo_copysign(delta, -g[1L])
    if (isTRUE(h1 > 0)) {
      dt <- -(g[1L] / h1)
      if (isTRUE(abs(dt) <= delta)) {
        d <- dt
        crvmin <- h1
      }
    }
    return(list(d = d, crvmin = crvmin))
  }
  hh <- .uo_hessenberg(h)
  h <- hh$A
  td <- hh$td
  tn <- hh$tn
  for (k in seq_len(n - 1L)) {
    ii <- (k + 1L):n
    s <- .pw_inprod(gg[ii], h[ii, k])
    gg[ii] <- .pw_fma(-s, h[ii, k], gg[ii])
  }
  tot <- .pw_sumabs(tn) + ((.pw_sumabs(gg) + .pw_sumabs(h)) + .pw_sumabs(td))
  if (!is.finite(tot)) return(list(d = numeric(n), crvmin = 0))
  hnorm <- .uo_maxval((abs(c(0, tn)) + abs(td)) + abs(c(tn, 0)))
  gnorm <- sqrt(gsq)
  delsq <- delta * delta
  tdmin <- .uo_minval(td)
  t0 <- if (isTRUE(tdmin < 0)) -tdmin else 0
  t1 <- gnorm / delta - hnorm
  parl <- if (isTRUE(t0 < t1)) t1 else t0
  parlest <- parl
  par <- parl
  paru <- 0
  paruest <- 0
  phil <- 0
  phiu <- 0
  posdef <- FALSE
  maxiter <- min(100L * n, 1000L)
  dold <- numeric(n)
  tnz <- c(tn, 0)
  for (iter in seq_len(maxiter)) {
    if (!is.finite(.pw_sumabs(d))) {
      d <- dold
      break
    }
    dold <- d
    piv <- numeric(n)
    piv[1L] <- td[1L] + par
    for (k in seq_len(n - 1L)) {
      if (isTRUE(piv[k] > 0)) {
        piv[k + 1L] <- (par + td[k + 1L]) - (tn[k] * tn[k]) / piv[k]
      } else if (isTRUE(abs(piv[k]) + abs(tn[k]) <= 0)) {
        piv[k + 1L] <- par + td[k + 1L]
      } else {
        break
      }
    }
    if (anyNA(piv)) break
    negc <- piv < 0 | (piv <= 0 & abs(tnz) > 0)
    negcrv <- any(negc)
    if (!negcrv) {
      if (gsq <= 0) {
        if (par <= 0) break
        paru <- par
        paruest <- par
      }
      k <- max(c(0L, which(abs(piv) + abs(tnz) <= 0)))
    } else {
      k <- min(c(n, which(negc)))
    }
    test190 <- FALSE
    if (negcrv || k > 0L) {
      d <- numeric(n)
      d[k] <- 1
      dhd <- piv[k]
      if (k < n && isTRUE(abs(tn[k]) > abs(piv[k]))) {
        temp <- par + td[k + 1L]
        piv[k + 1L] <- temp
        if (isTRUE(abs(piv[k]) >= temp)) {
          d[k + 1L] <- .uo_copysign(1, -tn[k])
          dhd <- .pw_fma(-abs(tn[k]), 2, temp + piv[k])
        } else {
          q <- tn[k] / temp
          dhd <- .pw_fma(-tn[k], q, piv[k])
          d[k + 1L] <- -q
        }
      }
      for (i in rev(seq_len(k - 1L))) {
        if (isTRUE(abs(tn[i]) > 0)) {
          d[i] <- -((tn[i] * d[i + 1L]) / piv[i])
        } else {
          d[seq_len(i)] <- 0
          break
        }
      }
      dsq <- .pw_inprod(d, d)
      parlest <- par - dhd / dsq
      parl <- par
      test190 <- TRUE
    } else if (gsq <= 0) {
      test190 <- TRUE
    }
    if (test190) {
      temp <- if (gsq <= 0) (1 - tol) * paruest else paruest
      if (isTRUE(paruest > 0) && isTRUE(temp <= parlest)) {
        dsq <- .pw_inprod(d, d)
        s <- delta / sqrt(dsq)
        d <- if (isTRUE(.pw_inprod(d, gg) > 0)) -(d * s) else d * s
        break
      }
    } else {
      d[1L] <- -(gg[1L] / piv[1L])
      for (i in seq_len(n - 1L) + 1L) d[i] <- -(.pw_fma(tn[i - 1L], d[i - 1L], gg[i]) / piv[i])
      wsq <- .pw_inprod(piv, d * d)
      for (i in rev(seq_len(n - 1L))) d[i] <- d[i] - (d[i + 1L] * tn[i]) / piv[i]
      if (!is.finite(.pw_sumabs(d))) {
        d <- dold
        break
      }
      dsq <- .pw_inprod(d, d)
      if (par <= 0 && isTRUE(delsq >= dsq)) {
        crvmin <- .uo_eigmin(td, tn, 0.01)
        break
      }
      dnorm <- sqrt(dsq)
      phi <- 1 / dnorm - 1 / delta
      temp <- .pw_fma(tol, (dsq * par) / wsq + 1, -((phi * dsq) * phi))
      if (isTRUE(temp >= 0)) {
        d <- d * (delta / dnorm)
        break
      }
      if (iter >= 2L && isTRUE(par <= parl)) break
      if (paru > 0 && isTRUE(par >= paru)) break
      if (isTRUE(phi < 0)) {
        if (posdef) {
          if (isTRUE(phi <= phil)) break
          slope <- (phi - phil) / (par - parl)
          parlest <- par - phi / slope
        } else {
          parlest <- par
        }
        slope <- if (paru > 0) (phiu - phi) / (paru - par) else 1 / gnorm
        temp <- par - phi / slope
        paruest <- if (paruest > 0) .pw_fminnm(paruest, temp) else temp
        posdef <- TRUE
        parl <- par
        phil <- phi
      } else if (posdef) {
        if (paru > 0) {
          if (isTRUE(phi >= phiu)) break
          slope <- (phiu - phi) / (paru - par)
        } else {
          slope <- 1 / gnorm
        }
        parlest <- .pw_fmaxnm(parlest, par - phi / slope)
        slope <- (phi - phil) / (par - parl)
        paruest <- par - phi / slope
        paru <- par
        phiu <- phi
      } else {
        w <- numeric(n)
        w[1L] <- 1 / piv[1L]
        for (i in seq_len(n - 1L) + 1L) {
          t <- w[i - 1L] * tn[i - 1L]
          w[i] <- if (isTRUE(t > 0)) -((t + 1) / piv[i]) else (1 - t) / piv[i]
        }
        wwsq <- .pw_inprod(piv, w * w)
        z <- w
        for (i in rev(seq_len(n - 1L))) z[i] <- w[i] - (z[i + 1L] * tn[i]) / piv[i]
        zsq <- .pw_inprod(z, z)
        dtz <- .pw_inprod(d, z)
        tempa <- abs(delsq - dsq)
        if (isTRUE(abs(dtz) > 0)) {
          tempb <- sqrt(.pw_fma(zsq, tempa, dtz * dtz))
          gam <- tempa / (.uo_copysign(tempb, dtz) + dtz)
        } else {
          gam <- sqrt(tempa / zsq)
        }
        temp <- .pw_fma(tol, .pw_fma(delsq, par, wsq), -((gam * gam) * wwsq))
        if (isTRUE(temp >= 0)) {
          d <- .pw_fma(gam, z, d)
          break
        }
        plc <- .pw_fmaxnm(par - wwsq / zsq, parlest)
        if (paru > 0) {
          if (isTRUE(phi >= phiu)) break
          slope <- (phiu - phi) / (paru - par)
        } else {
          slope <- 1 / gnorm
        }
        paru <- par
        paruest <- par
        phiu <- phi
        parlest <- .pw_fmaxnm(par - phi / slope, plc)
      }
    }
    par <- if (paru > 0) .pw_fmaxnm((parl + paru) * 0.5, parlest) else .pw_fma(parlest, 2, gnorm / delta)
    if (paruest > 0) par <- .pw_fminnm(par, paruest)
  }
  for (k in rev(seq_len(n - 1L))) {
    ii <- (k + 1L):n
    s <- .pw_inprod(d[ii], h[ii, k])
    d[ii] <- .pw_fma(-s, h[ii, k], d[ii])
  }
  if (isTRUE(delta < .pw_norm(d))) d <- d * (delta / .pw_norm(d))
  if (is.na(crvmin)) {
    crvmin <- 0
  } else if (scaled && crvmin > 0) {
    crvmin <- crvmin / scaling
  }
  list(d = d, crvmin = crvmin)
}

# initxf: the (n + 1)(n + 2)/2 initial interpolation points and their values
.uo_initxf <- function(fun, maxfun, ftarget, rhobeg, x0, iprint = 0L) {
  n <- length(x0)
  npt <- ((n + 1L) * (n + 2L)) %/% 2L
  xpt <- matrix(0, n, npt)
  for (i in seq_len(n)) xpt[i, 2L * i] <- rhobeg
  fval <- rep(.pr_realmax, npt)
  evaluated <- logical(npt)
  info <- 0L
  for (k in seq_len(2L * n + 1L)) {
    x <- x0 + xpt[, k]
    f <- .pw_evaluatef(fun, x)
    evaluated[k] <- TRUE
    fval[k] <- f
    .pw_fmsg(iprint, "Initialization", k, rhobeg, f, x)
    if (k %% 2L == 0L) xpt[, k + 1L] <- if (f < fval[1L]) xpt[, k] + xpt[, k] else -xpt[, k]
    info <- .pw_checkexit(maxfun, k, f, ftarget, x)
    if (info != 0L) break
  }
  if (info == 0L && n > 1L) {
    step <- ifelse(fval[2L * seq_len(n)] < fval[1L], rhobeg, -rhobeg)
    k <- 2L * n + 1L
    for (q in 2:n) {
      for (p in seq_len(q - 1L)) {
        k <- k + 1L
        xpt[p, k] <- step[p]
        xpt[q, k] <- step[q]
        x <- x0 + xpt[, k]
        f <- .pw_evaluatef(fun, x)
        evaluated[k] <- TRUE
        fval[k] <- f
        .pw_fmsg(iprint, "Initialization", k, rhobeg, f, x)
        info <- .pw_checkexit(maxfun, k, f, ftarget, x)
        if (info != 0L) break
      }
      if (info != 0L) break
    }
  }
  kopt <- which.min(ifelse(evaluated, fval, Inf))
  list(kopt = kopt, nf = sum(evaluated), fval = fval, xbase = x0, xpt = xpt, info = info)
}

# initl: the Lagrange functions of the initial interpolation set
.uo_initl <- function(xpt) {
  n <- nrow(xpt)
  npt <- ncol(xpt)
  rho <- max(abs(xpt[, 2L]))
  rr <- rho * rho
  pl <- matrix(0, npt - 1L, npt)
  a <- 0.5 / rho
  b <- -(2 / rr)
  cc <- 1 / rr
  for (i in seq_len(n)) {
    ih <- n + (i * (i + 1L)) %/% 2L
    if (xpt[i, 2L * i + 1L] > 0) {
      pl[i, 1L] <- -(1.5 / rho)
      pl[ih, 1L] <- cc
      pl[i, 2L * i] <- 2 / rho
      pl[ih, 2L * i] <- b
    } else {
      pl[ih, 1L] <- b
      pl[i, 2L * i] <- a
      pl[ih, 2L * i] <- cc
    }
    pl[i, 2L * i + 1L] <- -a
    pl[ih, 2L * i + 1L] <- cc
  }
  k <- 2L * n + 1L
  for (q in seq_len(n - 1L) + 1L) {
    for (p in seq_len(q - 1L)) {
      k <- k + 1L
      ih <- n + (q * (q - 1L)) %/% 2L + p
      t <- 1 / (xpt[q, k] * xpt[p, k])
      pl[ih, 1L] <- t
      pl[ih, k] <- t
      pl[ih, if (xpt[p, k] < 0) 2L * p + 1L else 2L * p] <- -t
      pl[ih, if (xpt[q, k] < 0) 2L * q + 1L else 2L * q] <- -t
    }
  }
  pl
}

# initq: the quadratic interpolant of the initial values
.uo_initq <- function(fval, xpt) {
  n <- nrow(xpt)
  npt <- ncol(xpt)
  rho <- max(abs(xpt[, 2L]))
  rr <- rho * rho
  tr <- rho + rho
  fb <- fval[1L]
  pq <- numeric(npt - 1L)
  hd <- numeric(n)
  for (i in seq_len(n)) {
    ih <- n + (i * (i + 1L)) %/% 2L
    f2 <- fval[2L * i]
    f3 <- fval[2L * i + 1L]
    if (xpt[i, 2L * i + 1L] > 0) {
      g <- (.pw_fma(f2, 4, -(fb * 3)) - f3) / tr
      hh <- .pw_fma(-f2, 2, fb + f3) / rr
    } else {
      g <- (f2 - f3) / tr
      hh <- .pw_fma(-fb, 2, f3 + f2) / rr
    }
    pq[i] <- g
    hd[i] <- hh
    pq[ih] <- hh
  }
  hr <- rr * 0.5
  k <- 2L * n + 1L
  for (q in seq_len(n - 1L) + 1L) {
    for (p in seq_len(q - 1L)) {
      k <- k + 1L
      ih <- n + (q * (q - 1L)) %/% 2L + p
      t <- .pw_fma(-pq[p], xpt[p, k], fval[k] - fb)
      s <- hd[q] + hd[p]
      t <- .pw_fma(-pq[q], xpt[q, k], t)
      t <- .pw_fma(-hr, s, t)
      pq[ih] <- t / (xpt[q, k] * xpt[p, k])
    }
  }
  pq
}

# update: PL, PQ and the interpolation set after XPT(:, knew) moves to
# XOPT + d; returns list(kopt, fval, pl, pq, xpt)
.uo_update <- function(knew, d, f, moderr, kopt, fval, pl, pq, xpt) {
  vlag <- .uo_calvlag(pl, d, xpt[, kopt], kopt)
  pk <- pl[, knew] / vlag[knew]
  pl <- pl - outer(pk, vlag)
  pl[, knew] <- pk
  pq <- .pw_fma(pk, moderr, pq)
  xpt[, knew] <- xpt[, kopt] + d
  fval[knew] <- f
  if (isTRUE(fval[knew] < fval[kopt])) kopt <- knew
  list(kopt = kopt, fval = fval, pl = pl, pq = pq, xpt = xpt)
}

# setdrop_tr: the point that the trust region step replaces (0 for none)
.uo_setdrop_tr <- function(kopt, ximproved, d, pl, rho, xpt) {
  x <- if (ximproved) xpt[, kopt] + d else xpt[, kopt]
  distsq <- .pw_distsq(xpt, x)
  w <- .pw_vfmaxnm(distsq / (rho * rho), 1)
  w <- w * w
  score <- abs(.uo_calvlag(pl, d, xpt[, kopt], kopt)) * (w * w)
  if (!ximproved) score[kopt] <- -1
  if (isTRUE(any(score > 1)) || (ximproved && isTRUE(any(score > 0)))) {
    k <- which.max(score)
    return(if (length(k)) k else 0L)
  }
  if (ximproved) return(.pw_maxloc(distsq))
  0L
}

# shiftbase_qint: moves the base point to XOPT
.uo_shiftbase <- function(kopt, pl, pq, xbase, xpt) {
  n <- nrow(xpt)
  m <- nrow(pl)
  xopt <- xpt[, kopt]
  xbase <- xbase + xopt
  xpt <- xpt - xopt
  hi <- (n + 1L):m
  pq[seq_len(n)] <- pq[seq_len(n)] + .uo_smat_mul_vec(pq[hi], xopt)
  for (k in seq_len(ncol(pl))) pl[seq_len(n), k] <- pl[seq_len(n), k] + .uo_smat_mul_vec(pl[hi, k], xopt)
  list(pl = pl, pq = pq, xbase = xbase, xpt = xpt)
}

# geostep: a step d of length delbar that makes |L_knew(XOPT + d)| large
# (Powell's LAGMAX), with the Cauchy step as a fallback
.uo_geostep <- function(knew, kopt, delbar, pl, xpt) {
  n <- nrow(xpt)
  hi <- (n + 1L):nrow(pl)
  xopt <- xpt[, kopt]
  g <- pl[seq_len(n), knew] + .uo_smat_mul_vec(pl[hi, knew], xopt)
  h <- .uo_vec2smat(pl[hi, knew], n)
  gg <- .pw_inprod(g, g)
  ghg <- .pw_inprod(g, .pw_matprod21(h, g))
  if (isTRUE(gg > 0) && is.finite(gg)) {
    dc <- g * (delbar / sqrt(gg))
    if (isTRUE(ghg < 0)) dc <- -dc
  } else {
    dc <- xpt[, knew] - xopt
    dc <- dc * .pw_fminnm(delbar / .pw_norm(dc), 0.5)
    if (isTRUE(.pw_inprod(g, dc) * .pw_inprod(dc, .pw_matprod21(h, dc)) < 0)) dc <- -dc
  }
  if (!isTRUE(any(abs(h) > 0)) || is.na(.pw_sumabs(h) + .pw_sumabs(g))) return(dc)
  if (n == 1L) return(if (isTRUE(g * h[1L, 1L] > 0)) delbar else -delbar)
  k <- .pw_maxloc(vapply(seq_len(n), function(j) .pw_inprod(h[, j], h[, j]), 0))
  v <- h[, k]
  v <- v / .pw_norm(v)
  vsq <- .pw_inprod(v, v)
  d <- .pw_matprod21(h, v)
  vhv <- .pw_inprod(v, d)
  dsq <- .pw_inprod(d, d)
  if (isTRUE(vhv * vhv <= (dsq * 0.9999) * vsq)) {
    d <- .pw_fma(-(vhv / vsq), v, d)
    wsq <- .pw_inprod(d, d)
    ratio <- sqrt(wsq / vsq)
    whw <- .pw_inprod(d, .pw_matprod21(h, d))
    v <- v * ratio
    rr <- ratio * ratio
    vhw <- ratio * wsq
    s <- .pw_fma(rr, vhv, whw)
    t <- .pw_fma(-rr, vhv, whw) * 0.5
    t <- if (isTRUE(s < 0)) t - sqrt(.pw_fma(vhw, vhw, t * t)) else t + sqrt(.pw_fma(t, t, vhw * vhw))
    d <- .pw_fma(vhw, v, t * d)
  }
  dd <- .pw_inprod(d, d)
  gd <- .pw_inprod(g, d)
  dhd <- .pw_inprod(d, .pw_matprod21(h, d))
  if (!(isTRUE(gg > 0) && isTRUE(dd > 0))) return(dc)
  v <- .pw_fma(-(gd / gg), g, d)
  vv <- .pw_inprod(v, v)
  scaling <- delbar / sqrt(dd)
  if (isTRUE(gd * dhd < 0)) scaling <- -scaling
  d <- d * scaling
  gnorm <- sqrt(gg)
  if (!(isTRUE((delbar * 0.005) * abs(dhd) < dd * gnorm) && isTRUE(vv > dd * 1e-4))) {
    return(if (isTRUE(.pw_sumabs(d) > 0)) d else dc)
  }
  hv <- .pw_matprod21(h, v)
  vhg <- .pw_inprod(g, hv)
  vhv <- .pw_inprod(v, hv)
  vnorm <- sqrt(vv)
  vhv <- vhv / vv
  ghg <- ghg / gg
  vhg <- vhg / (vnorm * gnorm)
  if (isTRUE(abs(vhg) <= .pw_fmaxnm(abs(ghg), abs(vhv)) * 0.01)) {
    vmu <- ghg - vhv
    wcos <- 1
    wsin <- 0
    wg <- gnorm
  } else {
    t <- (ghg - vhv) * 0.5
    vhg2 <- vhg * vhg
    r <- sqrt(.pw_fma(t, t, vhg2))
    vmu <- if (isTRUE(t < 0)) t - r else t + r
    r <- sqrt(.pw_fma(vmu, vmu, vhg2))
    wcos <- vmu / r
    wsin <- vhg / r
    wg <- wcos * gnorm
  }
  ta <- wcos / gnorm
  tc <- wcos / vnorm
  tdd <- wsin / gnorm
  tb <- wsin / vnorm
  d <- .pw_fma(ta, g, tb * v)
  v <- ifelse(.pw_isvec(n, 3), .pw_fma(-tdd, g, tc * v), .pw_fma(tc, v, -(tdd * g)))
  halfrt <- 0.7071067811865476
  dlin <- wg / delbar
  vl <- (wsin * gnorm) / delbar
  vlin <- -vl
  tempa <- .pw_fma(abs(vmu + vhv), 0.5, abs(dlin))
  tempb <- .pw_fma(abs(ghg - vmu), 0.5, abs(vl))
  tempc <- .pw_fma(abs(dlin) + abs(vl), halfrt, abs(ghg + vhv) * 0.25)
  if (isTRUE(tempa >= tempb) && isTRUE(tempa >= tempc)) {
    tempd <- if (isTRUE((vmu + vhv) * dlin < 0)) -delbar else delbar
    tempv <- 0
  } else if (isTRUE(tempb >= tempc)) {
    tempd <- 0
    tempv <- if (isTRUE((ghg - vmu) * vlin < 0)) -delbar else delbar
  } else {
    hd <- delbar * halfrt
    tempd <- if (isTRUE((ghg + vhv) * dlin < 0)) -hd else hd
    tempv <- if (isTRUE((ghg + vhv) * vlin < 0)) -hd else hd
  }
  d <- .pw_fma(tempd, d, tempv * v)
  if (is.na(.pw_sumabs(d))) dc else d
}

.uo_uobyqb <- function(fun, maxfun, ftarget, rhobeg, rhoend, x0, iprint = 0L) {
  eta1 <- 0.1
  eta2 <- 0.7
  gamma1 <- 0.5
  gamma2 <- 2
  n <- length(x0)
  npt <- ((n + 1L) * (n + 2L)) %/% 2L
  hi <- (n + 1L):(npt - 1L)
  ini <- .uo_initxf(fun, maxfun, ftarget, rhobeg, x0, iprint)
  kopt <- ini$kopt
  nf <- ini$nf
  fval <- ini$fval
  xbase <- ini$xbase
  xpt <- ini$xpt
  x <- xbase + xpt[, kopt]
  f <- fval[kopt]
  if (ini$info != 0L) {
    .pw_retmsg("UOBYQA", iprint, ini$info, nf, f, x)
    return(list(x = x, f = f, nf = nf, info = ini$info))
  }
  pl <- .uo_initl(xpt)
  pq <- .uo_initq(fval, xpt)
  rho <- rhobeg
  delta <- rho
  ratio <- -1
  dnormsav <- c(.pr_realmax, .pr_realmax)
  moderrsav <- c(.pr_realmax, .pr_realmax)
  gamma3 <- .pw_fmaxnm(.pw_fminnm(0.75 * gamma2, 1.5), 1)
  maxtr <- if (maxfun > .Machine$integer.max %/% 2L) maxfun else 2L * maxfun  # 2 maxfun overflows in Fortran
  info <- 20L
  dnorm <- 0
  d <- numeric(n)
  knew_tr <- 0L
  apply_update <- function(knew, moderr) {
    u <- .uo_update(knew, d, f, moderr, kopt, fval, pl, pq, xpt)
    kopt <<- u$kopt
    fval <<- u$fval
    pl <<- u$pl
    pq <<- u$pq
    xpt <<- u$xpt
  }
  for (tr in seq_len(maxtr)) {
    xopt <- xpt[, kopt]
    g <- pq[seq_len(n)] + .uo_smat_mul_vec(pq[hi], xopt)
    ts <- .uo_trstep(delta, g, .uo_vec2smat(pq[hi], n), 0.01)
    d <- ts$d
    crvmin <- ts$crvmin
    dnorm <- .pw_fminnm(.pw_norm(d), delta)
    shortd <- isTRUE(0.5 * rho > dnorm)
    qinc <- .uo_quadinc(pq, d, xopt)
    qred <- -qinc
    trfail <- !isTRUE(qred > (rho * rho) * .pr_eps)
    ddmove <- 0
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
      knew_tr <- .uo_setdrop_tr(kopt, ximproved, d, pl, rho, xpt)
      if (knew_tr > 0L) {
        dv <- xpt[, knew_tr] - xpt[, kopt]
        ddmove <- .pw_inprod(dv, dv)
        apply_update(knew_tr, moderr)
      }
      dnormsav <- c(dnormsav[2L], dnorm)
      moderrsav <- c(moderrsav[2L], moderr)
    }
    accurate_mod <- isTRUE(all(abs(moderrsav) <= (crvmin * 0.125) * (rho * rho))) && isTRUE(all(dnormsav <= rho))
    distsq <- .pw_distsq(xpt, xpt[, kopt])
    close_itpset <- isTRUE(all(distsq <= (delta * delta) * 4))
    adequate_geo <- (shortd && accurate_mod) || close_itpset
    small_trrad <- isTRUE(.pw_fmaxnm(dnorm, delta) <= rho)
    improve_geo <- (shortd || trfail || knew_tr == 0L ||
                      (isTRUE(ratio <= eta1) && isTRUE(ddmove <= (delta * delta) * 4))) && !adequate_geo
    reduce_rho <- (shortd || trfail || knew_tr == 0L || isTRUE(ratio <= 0)) && adequate_geo && small_trrad
    if (improve_geo) {
      knew_geo <- .pw_maxloc(distsq)
      delbar <- rho
      d <- .uo_geostep(knew_geo, kopt, delbar, pl, xpt)
      x <- (xpt[, kopt] + d) + xbase
      f <- .pw_evaluatef(fun, x)
      nf <- nf + 1L
      .pw_fmsg(iprint, "Geometry", nf, delbar, f, x)
      sub <- .pw_checkexit(maxfun, nf, f, ftarget, x)
      if (sub != 0L) {
        info <- sub
        break
      }
      dnorm <- .pw_fminnm(.pw_norm(d), delbar)
      moderr <- (f - fval[kopt]) - .uo_quadinc(pq, d, xpt[, kopt])
      apply_update(knew_geo, moderr)
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
    if (isTRUE(.pw_inprod(xopt, xopt) >= (delta * delta) * 1000)) {
      sb <- .uo_shiftbase(kopt, pl, pq, xbase, xpt)
      pl <- sb$pl
      pq <- sb$pq
      xbase <- sb$xbase
      xpt <- sb$xpt
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
  .pw_retmsg("UOBYQA", iprint, info, nf, f, x)
  list(x = x, f = f, nf = nf, info = info)
}

# uobyqa: preprocessing and the call of uobyqb
.uo_uobyqa <- function(fun, x0, rhobeg = 1, rhoend = 1e-6 * rhobeg, ftarget = -Inf,
                       maxfun = 500L * length(x0), iprint = 0L) {
  n <- length(x0)
  x0 <- .pw_moderatex(as.numeric(x0))
  op <- .pw_preproc("UOBYQA", n, maxfun, ftarget, rhobeg, rhoend)
  .uo_uobyqb(fun, op$maxfun, op$ftarget, op$rhobeg, op$rhoend, x0, as.integer(iprint))
}
