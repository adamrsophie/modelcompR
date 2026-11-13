# MINPACK (Garbow, Hillstrom and More, 1980) as translated to C in SciPy
# (scipy/optimize/src/minpack.c): the Powell hybrid method HYBRD and HYBRJ
# behind scipy.optimize.root(method = "hybr") and the Levenberg Marquardt
# method LMDIF and LMDER behind root(method = "lm") and
# least_squares(method = "lm"). Matrices are column major as in the C code;
# the packed triangle of HYBRD keeps the 1 based offsets of the Fortran.

.spmp_epsmch <- 2.220446049250313e-16

# sum in the order of the C loops (sum() accumulates in long double)
.spmp_ssum <- function(v) {
  s <- 0
  for (t in v) s <- s + t
  s
}

# enorm: the Euclidean norm guarded against overflow and underflow
.spmp_enorm <- function(x) {
  rgiant <- 2^63.5
  rdwarf <- 0.5 / rgiant
  agiant <- rgiant / length(x)
  ax <- abs(x)
  if (all(ax > rdwarf & ax < agiant)) return(sqrt(.spmp_ssum(ax * ax)))
  s1 <- 0; s2 <- 0; s3 <- 0; x1max <- 0; x3max <- 0
  for (xabs in ax) {
    if (xabs <= rdwarf || xabs >= agiant) {
      if (xabs > rdwarf) {
        if (xabs > x1max) {
          s1 <- 1 + s1 * (x1max / xabs)^2
          x1max <- xabs
        } else {
          s1 <- s1 + (xabs / x1max)^2
        }
      } else if (xabs > x3max) {
        s3 <- 1 + s3 * (x3max / xabs)^2
        x3max <- xabs
      } else if (xabs != 0) {
        s3 <- s3 + (xabs / x3max)^2
      }
    } else {
      s2 <- s2 + xabs * xabs
    }
  }
  if (s1 != 0) return(x1max * sqrt(s1 + (s2 / x1max) / x1max))
  if (s2 != 0) {
    if (s2 >= x3max) return(sqrt(s2 * (1 + (x3max / s2) * (x3max * s3))))
    return(sqrt(x3max * ((s2 / x3max) + (x3max * s3))))
  }
  x3max * sqrt(s3)
}

# fdjac1: forward difference Jacobian of n functions in n unknowns, banded
# when ml + mu + 1 < n
.spmp_fdjac1 <- function(fcn, x, fvec, ml, mu, epsfcn) {
  n <- length(x)
  eps <- sqrt(max(epsfcn, .spmp_epsmch))
  fjac <- matrix(0, n, n)
  msum <- ml + mu + 1L
  if (msum >= n) {
    for (j in seq_len(n)) {
      temp <- x[j]
      h <- eps * abs(temp)
      if (h == 0) h <- eps
      x[j] <- temp + h
      wa1 <- fcn(x)
      x[j] <- temp
      fjac[, j] <- (wa1 - fvec) / h
    }
    return(fjac)
  }
  i <- seq_len(n)
  for (k in seq_len(msum)) {
    js <- seq(k, n, by = msum)
    wa2 <- x[js]
    h <- eps * abs(wa2)
    h[h == 0] <- eps
    x[js] <- wa2 + h
    wa1 <- fcn(x)
    x[js] <- wa2
    for (q in seq_along(js)) {
      j <- js[q]
      fjac[, j] <- ifelse(i >= j - mu & i <= j + ml, (wa1 - fvec) / h[q], 0)
    }
  }
  fjac
}

# fdjac2: forward difference Jacobian of m functions in n unknowns
.spmp_fdjac2 <- function(fcn, x, fvec, epsfcn) {
  eps <- sqrt(max(epsfcn, .spmp_epsmch))
  fjac <- matrix(0, length(fvec), length(x))
  for (j in seq_along(x)) {
    temp <- x[j]
    h <- eps * abs(temp)
    if (h == 0) h <- eps
    x[j] <- temp + h
    wa <- fcn(x)
    x[j] <- temp
    fjac[, j] <- (wa - fvec) / h
  }
  fjac
}

# qrfac: Householder QR with optional column pivoting; R sits above the
# diagonal of a, its diagonal in rdiag, and the Householder vectors below
.spmp_qrfac <- function(a, pivot) {
  m <- nrow(a); n <- ncol(a)
  acnorm <- vapply(seq_len(n), function(j) .spmp_enorm(a[, j]), 0)
  rdiag <- acnorm
  wa <- rdiag
  ipvt <- seq_len(n)
  for (j in seq_len(min(m, n))) {
    if (pivot) {
      kmax <- j
      for (k in j:n) if (rdiag[k] > rdiag[kmax]) kmax <- k
      if (kmax != j) {
        temp <- a[, j]; a[, j] <- a[, kmax]; a[, kmax] <- temp
        rdiag[kmax] <- rdiag[j]
        wa[kmax] <- wa[j]
        k <- ipvt[j]; ipvt[j] <- ipvt[kmax]; ipvt[kmax] <- k
      }
    }
    rows <- j:m
    ajnorm <- .spmp_enorm(a[rows, j])
    if (ajnorm == 0) {
      rdiag[j] <- -ajnorm
      next
    }
    if (a[j, j] < 0) ajnorm <- -ajnorm
    a[rows, j] <- a[rows, j] / ajnorm
    a[j, j] <- a[j, j] + 1
    if (j < n) {
      for (k in (j + 1L):n) {
        temp <- .spmp_ssum(a[rows, j] * a[rows, k]) / a[j, j]
        a[rows, k] <- a[rows, k] - temp * a[rows, j]
        if (pivot && rdiag[k] != 0) {
          temp <- a[j, k] / rdiag[k]
          rdiag[k] <- rdiag[k] * sqrt(max(0, 1 - temp * temp))
          if (0.05 * (rdiag[k] / wa[k])^2 <= .spmp_epsmch) {
            rdiag[k] <- if (j < m) .spmp_enorm(a[(j + 1L):m, k]) else 0
            wa[k] <- rdiag[k]
          }
        }
      }
    }
    rdiag[j] <- -ajnorm
  }
  list(a = a, ipvt = ipvt, rdiag = rdiag, acnorm = acnorm)
}

# qform: the orthogonal matrix Q of the square factorization of qrfac
.spmp_qform <- function(q) {
  m <- nrow(q)
  if (m > 1L) for (j in 2:m) q[seq_len(j - 1L), j] <- 0
  for (k in m:1) {
    rows <- k:m
    wa <- q[rows, k]
    q[rows, k] <- 0
    q[k, k] <- 1
    if (wa[1L] != 0) {
      for (j in k:m) {
        temp <- .spmp_ssum(q[rows, j] * wa) / wa[1L]
        q[rows, j] <- q[rows, j] - temp * wa
      }
    }
  }
  q
}

# dogleg: the dogleg step in the trust region of radius delta; r is the
# upper triangle packed by rows
.spmp_dogleg <- function(n, r, diag, qtb, delta) {
  x <- numeric(n)
  jj <- (n * (n + 1L)) %/% 2L
  for (k in 0:(n - 1L)) {
    j <- n - k
    jj <- jj - (k + 1L)
    l <- jj + 2L
    ssum <- 0
    if (j < n) {
      for (i in (j + 1L):n) {
        ssum <- ssum + r[l] * x[i]
        l <- l + 1L
      }
    }
    temp <- r[jj + 1L]
    if (temp == 0) {
      l <- j
      for (i in seq_len(j)) {
        temp <- max(temp, abs(r[l]))
        l <- l + n - i
      }
      temp <- .spmp_epsmch * temp
      if (temp == 0) temp <- .spmp_epsmch
    }
    x[j] <- (qtb[j] - ssum) / temp
  }
  wa1 <- numeric(n)
  qnorm <- .spmp_enorm(diag * x)
  if (qnorm <= delta) return(x)
  l <- 1L
  for (j in seq_len(n)) {
    temp <- qtb[j]
    for (i in j:n) {
      wa1[i] <- wa1[i] + r[l] * temp
      l <- l + 1L
    }
    wa1[j] <- wa1[j] / diag[j]
  }
  gnorm <- .spmp_enorm(wa1)
  sgnorm <- 0
  alpha <- delta / qnorm
  if (gnorm != 0) {
    wa1 <- (wa1 / gnorm) / diag
    wa2 <- numeric(n)
    l <- 1L
    for (j in seq_len(n)) {
      ssum <- 0
      for (i in j:n) {
        ssum <- ssum + r[l] * wa1[i]
        l <- l + 1L
      }
      wa2[j] <- ssum
    }
    temp <- .spmp_enorm(wa2)
    sgnorm <- (gnorm / temp) / temp
    alpha <- 0
    if (sgnorm < delta) {
      bnorm <- .spmp_enorm(qtb)
      temp <- (bnorm / gnorm) * (bnorm / qnorm) * (sgnorm / delta)
      temp <- temp - (delta / qnorm) * (sgnorm / delta)^2 +
        sqrt((temp - (delta / qnorm))^2 + ((1 - (delta / qnorm)^2) * (1 - (sgnorm / delta)^2)))
      alpha <- ((delta / qnorm) * (1 - (sgnorm / delta)^2)) / temp
    }
  }
  temp <- (1 - alpha) * min(sgnorm, delta)
  temp * wa1 + alpha * x
}

# r1updt: the rank one update of the packed lower trapezoid s (stored by
# columns) through Givens rotations; returns s, v, w and sing
.spmp_r1updt <- function(m, n, s, u, v) {
  giant <- .Machine$double.xmax
  w <- numeric(m)
  jj <- (n * (2L * m - n + 1L)) %/% 2L - (m - n)
  w[n:m] <- s[jj + 0:(m - n)]
  if (n > 1L) {
    for (nmj in seq_len(n - 1L)) {
      j <- n - nmj
      jj <- jj - (m - j + 1L)
      w[j] <- 0
      if (v[j] != 0) {
        if (abs(v[n]) < abs(v[j])) {
          cotan <- v[n] / v[j]
          sn <- 0.5 / sqrt(0.25 + 0.25 * (cotan * cotan))
          cs <- sn * cotan
          tau <- 1
          if (abs(cs) * giant > 1) tau <- 1 / cs
        } else {
          tn <- v[j] / v[n]
          cs <- 0.5 / sqrt(0.25 + 0.25 * (tn * tn))
          sn <- cs * tn
          tau <- sn
        }
        v[n] <- sn * v[j] + cs * v[n]
        v[j] <- tau
        idx <- jj + 0:(m - j)
        temp <- cs * s[idx] - sn * w[j:m]
        w[j:m] <- sn * s[idx] + cs * w[j:m]
        s[idx] <- temp
      }
    }
  }
  w <- w + v[n] * u
  sing <- FALSE
  if (n > 1L) {
    for (j in seq_len(n - 1L)) {
      if (w[j] != 0) {
        if (abs(s[jj]) < abs(w[j])) {
          cotan <- s[jj] / w[j]
          sn <- 0.5 / sqrt(0.25 + 0.25 * cotan^2)
          cs <- sn * cotan
          tau <- 1
          if (abs(cs) * giant > 1) tau <- 1 / cs
        } else {
          tn <- w[j] / s[jj]
          cs <- 0.5 / sqrt(0.25 + 0.25 * tn^2)
          sn <- cs * tn
          tau <- sn
        }
        idx <- jj + 0:(m - j)
        temp <- cs * s[idx] + sn * w[j:m]
        w[j:m] <- -sn * s[idx] + cs * w[j:m]
        s[idx] <- temp
        w[j] <- tau
      }
      if (s[jj] == 0) sing <- TRUE
      jj <- jj + (m - j + 1L)
    }
  }
  s[jj + 0:(m - n)] <- w[n:m]
  if (s[jj] == 0) sing <- TRUE
  list(s = s, v = v, w = w, sing = sing)
}

# r1mpyq: a times the product of the Givens rotations stored in v and w
.spmp_r1mpyq <- function(a, v, w) {
  n <- ncol(a)
  if (n < 2L) return(a)
  rot <- function(t) if (abs(t) > 1) c(1 / t, sqrt(1 - (1 / t) * (1 / t))) else c(sqrt(1 - t * t), t)
  for (j in (n - 1L):1) {
    cs <- rot(v[j])
    temp <- cs[1L] * a[, j] - cs[2L] * a[, n]
    a[, n] <- cs[2L] * a[, j] + cs[1L] * a[, n]
    a[, j] <- temp
  }
  for (j in seq_len(n - 1L)) {
    cs <- rot(w[j])
    temp <- cs[1L] * a[, j] + cs[2L] * a[, n]
    a[, n] <- -cs[2L] * a[, j] + cs[1L] * a[, n]
    a[, j] <- temp
  }
  a
}

# qrsolv: solves the least squares problem of R augmented by diag(d)
.spmp_qrsolv <- function(n, r, ipvt, d, qtb) {
  for (j in seq_len(n - 1L)) r[(j + 1L):n, j] <- r[j, (j + 1L):n]
  x <- r[cbind(seq_len(n), seq_len(n))]
  wa <- qtb[seq_len(n)]
  sdiag <- numeric(n)
  for (j in seq_len(n)) {
    l <- ipvt[j]
    if (d[l] != 0) {
      sdiag[j:n] <- 0
      sdiag[j] <- d[l]
      qtbpj <- 0
      for (k in j:n) {
        if (sdiag[k] != 0) {
          if (abs(r[k, k]) < abs(sdiag[k])) {
            cotan <- r[k, k] / sdiag[k]
            sn <- 0.5 / sqrt(0.25 + 0.25 * cotan^2)
            cs <- sn * cotan
          } else {
            tn <- sdiag[k] / r[k, k]
            cs <- 0.5 / sqrt(0.25 + 0.25 * tn^2)
            sn <- cs * tn
          }
          r[k, k] <- r[k, k] * cs
          r[k, k] <- r[k, k] + sn * sdiag[k]
          temp <- cs * wa[k] + sn * qtbpj
          qtbpj <- -sn * wa[k] + cs * qtbpj
          wa[k] <- temp
          if (k < n) {
            i <- (k + 1L):n
            temp <- cs * r[i, k] + sn * sdiag[i]
            sdiag[i] <- -sn * r[i, k] + cs * sdiag[i]
            r[i, k] <- temp
          }
        }
      }
    }
    sdiag[j] <- r[j, j]
    r[j, j] <- x[j]
  }
  nsing <- n
  for (j in seq_len(n)) {
    if (sdiag[j] == 0 && nsing == n) nsing <- j - 1L
    if (nsing < n) wa[j] <- 0
  }
  for (k in seq_len(nsing)) {
    j <- nsing - k + 1L
    ssum <- if (j < nsing) .spmp_ssum(r[(j + 1L):nsing, j] * wa[(j + 1L):nsing]) else 0
    wa[j] <- (wa[j] - ssum) / sdiag[j]
  }
  out <- numeric(n)
  out[ipvt] <- wa
  list(r = r, x = out, sdiag = sdiag)
}

# lmpar: the Levenberg Marquardt parameter for the trust region delta
.spmp_lmpar <- function(n, r, ipvt, d, qtb, delta, par) {
  dwarf <- 0x1.0000000000000p-1022
  wa1 <- qtb[seq_len(n)]
  nsing <- n
  for (j in seq_len(n)) {
    if (r[j, j] == 0 && nsing == n) nsing <- j - 1L
    if (nsing < n) wa1[j] <- 0
  }
  for (k in seq_len(nsing)) {
    j <- nsing - k + 1L
    wa1[j] <- wa1[j] / r[j, j]
    temp <- wa1[j]
    if (j > 1L) {
      i <- seq_len(j - 1L)
      wa1[i] <- wa1[i] - r[i, j] * temp
    }
  }
  x <- numeric(n)
  x[ipvt] <- wa1
  iter <- 0L
  wa2 <- d * x
  dxnorm <- .spmp_enorm(wa2)
  fp <- dxnorm - delta
  if (fp <= 0.1 * delta) return(list(par = 0, x = x))
  parl <- 0
  if (nsing == n) {
    wa1 <- d[ipvt] * (wa2[ipvt] / dxnorm)
    for (j in seq_len(n)) {
      ssum <- if (j > 1L) .spmp_ssum(r[seq_len(j - 1L), j] * wa1[seq_len(j - 1L)]) else 0
      wa1[j] <- (wa1[j] - ssum) / r[j, j]
    }
    temp <- .spmp_enorm(wa1)
    parl <- ((fp / delta) / temp) / temp
  }
  for (j in seq_len(n)) wa1[j] <- .spmp_ssum(r[seq_len(j), j] * qtb[seq_len(j)]) / d[ipvt[j]]
  gnorm <- .spmp_enorm(wa1)
  paru <- gnorm / delta
  if (paru == 0) paru <- dwarf / min(delta, 0.1)
  par <- min(max(par, parl), paru)
  if (par == 0) par <- gnorm / dxnorm
  repeat {
    iter <- iter + 1L
    if (par == 0) par <- max(dwarf, paru * 0.001)
    qs <- .spmp_qrsolv(n, r, ipvt, sqrt(par) * d, qtb)
    r <- qs$r
    x <- qs$x
    sdiag <- qs$sdiag
    wa2 <- d * x
    dxnorm <- .spmp_enorm(wa2)
    temp <- fp
    fp <- dxnorm - delta
    if (abs(fp) <= 0.1 * delta || (parl == 0 && fp <= temp && temp < 0) || iter == 10L)
      return(list(par = par, x = x))
    wa1 <- d[ipvt] * (wa2[ipvt] / dxnorm)
    for (j in seq_len(n)) {
      wa1[j] <- wa1[j] / sdiag[j]
      temp <- wa1[j]
      if (j < n) {
        i <- (j + 1L):n
        wa1[i] <- wa1[i] - r[i, j] * temp
      }
    }
    temp <- .spmp_enorm(wa1)
    parc <- ((fp / delta) / temp) / temp
    if (fp > 0) parl <- max(parl, par)
    if (fp < 0) paru <- min(paru, par)
    par <- max(parl, par + parc)
  }
}

# HYBRD (jac = NULL, forward differences) and HYBRJ (jac a function giving
# the n by n Jacobian): the Powell hybrid method for fcn(x) = 0
.spmp_hybrd <- function(fcn, x, jac = NULL, xtol = 1.49012e-08, maxfev = 200L * (length(x) + 1L), ml = length(x) - 1L,
                      mu = length(x) - 1L, epsfcn = .spmp_epsmch, diag = NULL, factor = 100) {
  n <- length(x)
  mode <- if (is.null(diag)) 1L else 2L
  if (mode == 1L) diag <- numeric(n)
  out <- function(info) list(x = x, fvec = fvec, info = info, nfev = nfev, njev = njev, fjac = fjac, r = r, qtf = qtf)
  fvec <- NULL; fjac <- NULL; r <- NULL; qtf <- NULL
  nfev <- 0L; njev <- 0L
  if (n <= 0L || xtol < 0 || maxfev <= 0L || (is.null(jac) && (ml < 0L || mu < 0L)) || factor <= 0 ||
      (mode == 2L && any(diag <= 0))) return(out(0L))
  fvec <- fcn(x)
  nfev <- 1L
  fnorm <- .spmp_enorm(fvec)
  msum <- min(ml + mu + 1L, n)
  iter <- 1L
  ncsuc <- 0L; ncfail <- 0L; nslow1 <- 0L; nslow2 <- 0L
  info <- 0L
  repeat {
    jeval <- TRUE
    if (is.null(jac)) {
      fjac <- .spmp_fdjac1(fcn, x, fvec, ml, mu, epsfcn)
      nfev <- nfev + msum
    } else {
      fjac <- matrix(as.numeric(jac(x)), n, n)
      njev <- njev + 1L
    }
    qr <- .spmp_qrfac(fjac, FALSE)
    fjac <- qr$a
    wa1 <- qr$rdiag
    wa2 <- qr$acnorm
    if (iter == 1L) {
      if (mode != 2L) diag <- ifelse(wa2 == 0, 1, wa2)
      xnorm <- .spmp_enorm(diag * x)
      delta <- factor * xnorm
      if (delta == 0) delta <- factor
    }
    qtf <- fvec
    for (j in seq_len(n)) {
      if (fjac[j, j] != 0) {
        rows <- j:n
        temp <- -.spmp_ssum(fjac[rows, j] * qtf[rows]) / fjac[j, j]
        qtf[rows] <- qtf[rows] + fjac[rows, j] * temp
      }
    }
    sing <- FALSE
    r <- numeric((n * (n + 1L)) %/% 2L)
    for (j in seq_len(n)) {
      l <- j
      for (i in seq_len(j - 1L)) {
        r[l] <- fjac[i, j]
        l <- l + n - i
      }
      r[l] <- wa1[j]
      if (wa1[j] == 0) sing <- TRUE
    }
    fjac <- .spmp_qform(fjac)
    if (mode != 2L) diag <- pmax(diag, wa2)
    repeat {
      wa1 <- -.spmp_dogleg(n, r, diag, qtf, delta)
      wa2 <- x + wa1
      pnorm <- .spmp_enorm(diag * wa1)
      if (iter == 1L) delta <- min(delta, pnorm)
      wa4 <- fcn(wa2)
      nfev <- nfev + 1L
      fnorm1 <- .spmp_enorm(wa4)
      actred <- if (fnorm1 < fnorm) 1 - (fnorm1 / fnorm)^2 else -1
      wa3 <- numeric(n)
      l <- 1L
      for (i in seq_len(n)) {
        ssum <- 0
        for (j in i:n) {
          ssum <- ssum + r[l] * wa1[j]
          l <- l + 1L
        }
        wa3[i] <- qtf[i] + ssum
      }
      temp <- .spmp_enorm(wa3)
      prered <- if (temp < fnorm) 1 - (temp / fnorm)^2 else 0
      ratio <- if (prered > 0) actred / prered else 0
      if (ratio < 0.1) {
        ncsuc <- 0L
        ncfail <- ncfail + 1L
        delta <- 0.5 * delta
      } else {
        ncfail <- 0L
        ncsuc <- ncsuc + 1L
        if (ratio >= 0.5 || ncsuc > 1L) delta <- max(delta, pnorm / 0.5)
        if (abs(ratio - 1) <= 0.1) delta <- pnorm / 0.5
      }
      if (ratio >= 1e-4) {
        x <- wa2
        fvec <- wa4
        xnorm <- .spmp_enorm(diag * x)
        fnorm <- fnorm1
        iter <- iter + 1L
      }
      nslow1 <- nslow1 + 1L
      if (actred >= 0.001) nslow1 <- 0L
      if (jeval) nslow2 <- nslow2 + 1L
      if (actred >= 0.1) nslow2 <- 0L
      if (delta <= xtol * xnorm || fnorm == 0) info <- 1L
      if (info != 0L) return(out(info))
      if (nfev >= maxfev) info <- 2L
      if (0.1 * max(0.1 * delta, pnorm) <= .spmp_epsmch * xnorm) info <- 3L
      if (nslow2 == 5L) info <- 4L
      if (nslow1 == 10L) info <- 5L
      if (info != 0L) return(out(info))
      if (ncfail == 2L) break
      ssum <- vapply(seq_len(n), function(j) .spmp_ssum(fjac[, j] * wa4), 0)
      wa2 <- (ssum - wa3) / pnorm
      wa1 <- diag * ((diag * wa1) / pnorm)
      if (ratio >= 1e-4) qtf <- ssum
      up <- .spmp_r1updt(n, n, r, wa1, wa2)
      r <- up$s
      sing <- up$sing
      fjac <- .spmp_r1mpyq(fjac, up$v, up$w)
      qtf <- as.numeric(.spmp_r1mpyq(matrix(qtf, 1L), up$v, up$w))
      jeval <- FALSE
    }
  }
}

# LMDIF (jac = NULL, forward differences) and LMDER (jac a function giving
# the m by n Jacobian): the Levenberg Marquardt method for min sum fcn(x)^2
.spmp_lmdif <- function(fcn, x, jac = NULL, ftol = 1.49012e-08, xtol = 1.49012e-08, gtol = 0, maxfev = NULL,
                      epsfcn = .spmp_epsmch, diag = NULL, factor = 100) {
  n <- length(x)
  mode <- if (is.null(diag)) 1L else 2L
  if (mode == 1L) diag <- numeric(n)
  fvec <- NULL; fjac <- NULL; ipvt <- NULL; qtf <- NULL
  nfev <- 0L; njev <- 0L
  out <- function(info) list(x = x, fvec = fvec, info = info, nfev = nfev, njev = njev, fjac = fjac, ipvt = ipvt,
                             qtf = qtf)
  if (n <= 0L || ftol < 0 || xtol < 0 || gtol < 0 || maxfev <= 0L || factor <= 0 || (mode == 2L && any(diag <= 0)))
    return(out(0L))
  fvec <- as.numeric(fcn(x))
  m <- length(fvec)
  if (m < n) return(out(0L))
  nfev <- 1L
  fnorm <- .spmp_enorm(fvec)
  par <- 0
  iter <- 1L
  delta <- 0
  xnorm <- 0
  info <- 0L
  repeat {
    if (is.null(jac)) {
      fjac <- .spmp_fdjac2(fcn, x, fvec, epsfcn)
      nfev <- nfev + n
    } else {
      fjac <- matrix(as.numeric(jac(x)), m, n)
      njev <- njev + 1L
    }
    qr <- .spmp_qrfac(fjac, TRUE)
    fjac <- qr$a
    ipvt <- qr$ipvt
    wa1 <- qr$rdiag
    wa2 <- qr$acnorm
    if (iter == 1L) {
      if (mode != 2L) diag <- ifelse(wa2 == 0, 1, wa2)
      xnorm <- .spmp_enorm(diag * x)
      delta <- factor * xnorm
      if (delta == 0) delta <- factor
    }
    wa4 <- fvec
    qtf <- numeric(n)
    for (j in seq_len(n)) {
      if (fjac[j, j] != 0) {
        rows <- j:m
        temp <- -.spmp_ssum(fjac[rows, j] * wa4[rows]) / fjac[j, j]
        wa4[rows] <- wa4[rows] + fjac[rows, j] * temp
      }
      fjac[j, j] <- wa1[j]
      qtf[j] <- wa4[j]
    }
    gnorm <- 0
    if (fnorm != 0) {
      for (j in seq_len(n)) {
        l <- ipvt[j]
        if (wa2[l] != 0) {
          ssum <- .spmp_ssum(fjac[seq_len(j), j] * (qtf[seq_len(j)] / fnorm))
          gnorm <- max(gnorm, abs(ssum / wa2[l]))
        }
      }
    }
    if (gnorm <= gtol) info <- 4L
    if (info != 0L) return(out(info))
    if (mode != 2L) diag <- pmax(diag, wa2)
    repeat {
      lp <- .spmp_lmpar(n, fjac, ipvt, diag, qtf, delta, par)
      par <- lp$par
      wa1 <- -lp$x
      wa2 <- x + wa1
      pnorm <- .spmp_enorm(diag * wa1)
      if (iter == 1L) delta <- min(delta, pnorm)
      wa4 <- as.numeric(fcn(wa2))
      nfev <- nfev + 1L
      fnorm1 <- .spmp_enorm(wa4)
      actred <- if (0.1 * fnorm1 < fnorm) 1 - (fnorm1 / fnorm)^2 else -1
      wa3 <- numeric(n)
      for (j in seq_len(n)) {
        i <- seq_len(j)
        wa3[i] <- wa3[i] + fjac[i, j] * wa1[ipvt[j]]
      }
      temp1 <- .spmp_enorm(wa3) / fnorm
      temp2 <- (sqrt(par) * pnorm) / fnorm
      prered <- temp1 * temp1 + temp2 * temp2 / 0.5
      dirder <- -(temp1 * temp1 + temp2 * temp2)
      ratio <- if (prered != 0) actred / prered else 0
      if (ratio <= 0.25) {
        temp <- if (actred < 0) (0.5 * dirder) / (dirder + 0.5 * actred) else 0.5
        if (0.1 * fnorm1 >= fnorm || temp < 0.1) temp <- 0.1
        delta <- temp * min(delta, pnorm / 0.1)
        par <- par / temp
      } else if (par == 0 || ratio >= 0.75) {
        delta <- pnorm / 0.5
        par <- 0.5 * par
      }
      if (ratio >= 1e-4) {
        x <- wa2
        fvec <- wa4
        xnorm <- .spmp_enorm(diag * x)
        fnorm <- fnorm1
        iter <- iter + 1L
      }
      conv <- abs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1
      if (conv) info <- 1L
      if (delta <= xtol * xnorm) info <- 2L
      if (conv && info == 2L) info <- 3L
      if (info != 0L) return(out(info))
      if (nfev >= maxfev) info <- 5L
      if (abs(actred) <= .spmp_epsmch && prered <= .spmp_epsmch && 0.5 * ratio <= 1) info <- 6L
      if (delta <= .spmp_epsmch * xnorm) info <- 7L
      if (gnorm <= .spmp_epsmch) info <- 8L
      if (info != 0L) return(out(info))
      if (ratio >= 1e-4) break
    }
  }
}

# ScipyRoot: scipy.optimize.root with the method of the algorithm
.scipy_run_root <- function(alg, f, grad, x0, kw) {
  method <- alg$method
  fun <- function(x) as.numeric(f(x))
  n <- length(x0)
  jac <- if (!is.null(grad) && method %in% c("hybr", "lm")) grad
  # the keyword arguments of scipy.optimize.root: tol, callback and options
  bad <- setdiff(names(kw), c("tol", "callback", "options"))
  if (length(bad)) stop("root() got an unexpected keyword argument '", bad[1L], "'", call. = FALSE)
  o <- kw$options %||% list()
  if (!is.null(kw$tol)) {
    if (method %in% c("hybr", "lm")) o$xtol <- o$xtol %||% kw$tol
    else if (method == "df-sane") o$ftol <- o$ftol %||% kw$tol
    else {
      o$xtol <- o$xtol %||% kw$tol
      o$xatol <- o$xatol %||% Inf
      o$ftol <- o$ftol %||% Inf
      o$fatol <- o$fatol %||% Inf
    }
  }
  if (method == "hybr") {
    xtol <- o$xtol %||% 1.49012e-08
    maxfev <- o$maxfev %||% 0L
    if (maxfev == 0L) maxfev <- if (is.null(jac)) 200L * (n + 1L) else 100L * (n + 1L)
    band <- o$band
    r <- .spmp_hybrd(fun, x0, jac = jac, xtol = xtol, maxfev = maxfev,
                   ml = if (is.null(band) || band[1L] < 0) n - 1L else as.integer(band[1L]),
                   mu = if (is.null(band) || band[2L] < 0) n - 1L else as.integer(band[2L]),
                   epsfcn = o$eps %||% .spmp_epsmch, diag = o$diag, factor = o$factor %||% 100)
    msg <- switch(as.character(r$info),
                  "0" = "Improper input parameters were entered.",
                  "1" = "The solution converged.",
                  "2" = paste0("The number of calls to function has reached maxfev = ", maxfev, "."),
                  "3" = sprintf("xtol=%f is too small, no further improvement in the approximate\n solution is possible.",
                                xtol),
                  "4" = paste0("The iteration is not making good progress, as measured by the \n improvement from ",
                               "the last five Jacobian evaluations."),
                  "5" = paste0("The iteration is not making good progress, as measured by the \n improvement from ",
                               "the last ten iterations."))
    # root() counts every call of the function: the shape checks of
    # _root_hybr and of the C wrapper and those of MINPACK
    res <- list(x = r$x, success = r$info == 1L, status = r$info, method = "hybr", fun = r$fvec, nfev = r$nfev + 2L,
                fjac = r$fjac, r = r$r, qtf = r$qtf, message = msg)
    if (!is.null(jac)) res$njev <- r$njev
  } else if (method == "lm") {
    xtol <- o$xtol %||% 1.49012e-08
    ftol <- o$ftol %||% 1.49012e-08
    gtol <- o$gtol %||% 0
    maxfev <- o$maxiter %||% 0L
    if (maxfev == 0L) maxfev <- if (is.null(jac)) 200L * (n + 1L) else 100L * (n + 1L)
    r <- .spmp_lmdif(fun, x0, jac = jac, ftol = ftol, xtol = xtol, gtol = gtol, maxfev = maxfev, epsfcn = o$eps %||% 0,
                   diag = o$diag, factor = o$factor %||% 100)
    msg <- .spmp_leastsq_message(r$info, ftol, xtol, gtol, maxfev)
    res <- list(x = r$x, message = msg, status = r$info, success = r$info %in% 1:4, fun = r$fvec, method = "lm",
                fjac = t(r$fjac), ipvt = r$ipvt, qtf = r$qtf, nfev = r$nfev + 2L)
    if (!is.null(jac)) res$njev <- r$njev
  } else if (method == "df-sane") {
    known <- c("ftol", "fatol", "maxfev", "fnorm", "M", "eta_strategy", "sigma_eps", "sigma_0", "line_search", "disp")
    bad <- setdiff(names(o), known)
    if (length(bad)) warning("Unknown solver options: ", paste(bad, collapse = ", "), call. = FALSE)
    res <- do.call(.sp_df_sane, c(list(fun, x0), o[intersect(names(o), setdiff(known, "disp"))]))
  } else {
    res <- .sp_nonlin_solve(fun, x0, method, o)
  }
  res$retcode <- if (isTRUE(res$success)) "Success" else "Failure"
  # the objective of the Julia wrapper: sum(abs2, fun)
  res$objective <- .spmp_ssum(res$fun * res$fun)
  res
}

# a VectorFunction of SciPy: the residuals cached at the last point and the
# forward difference Jacobian taken from the cached residuals
.sp_vector_function <- function(fun, diff_step = NULL, lb = -Inf, ub = Inf) {
  e <- new.env(parent = emptyenv())
  e$fun <- function(x) {
    if (is.null(e$x) || !.sp_same_x(x, e$x)) {
      e$x <- x
      e$f <- as.numeric(fun(x))
    }
    e$f
  }
  e$jac <- function(x) {
    f0 <- e$fun(x)
    matrix(.sp_approx_derivative(fun, x, "2-point", rel_step = diff_step, f0 = f0, lb = lb, ub = ub)$J,
           length(f0), length(x))
  }
  e
}

.sp_lsq_messages <- c(
  "-2" = "Stopped because `callback` function raised `StopIteration` or returned `True`",
  "-1" = "Improper input parameters status returned from `leastsq`",
  "0" = "The maximum number of function evaluations is exceeded.",
  "1" = "`gtol` termination condition is satisfied.",
  "2" = "`ftol` termination condition is satisfied.",
  "3" = "`xtol` termination condition is satisfied.",
  "4" = "Both `ftol` and `xtol` termination conditions are satisfied.")


.spmp_leastsq_message <- function(info, ftol, xtol, gtol, maxfev) {
  switch(as.character(info),
         "0" = "Improper input parameters.",
         "1" = sprintf("Both actual and predicted relative reductions in the sum of squares\n  are at most %f", ftol),
         "2" = sprintf("The relative error between two consecutive iterates is at most %f", xtol),
         "3" = sprintf(paste0("Both actual and predicted relative reductions in the sum of squares\n  are at most %f ",
                              "and the relative error between two consecutive iterates is at \n  most %f"), ftol, xtol),
         "4" = sprintf(paste0("The cosine of the angle between func(x) and any column of the\n  Jacobian is at most ",
                              "%f in absolute value"), gtol),
         "5" = paste0("Number of calls to function has reached maxfev = ", maxfev, "."),
         "6" = sprintf("ftol=%f is too small, no further reduction in the sum of squares\n  is possible.", ftol),
         "7" = sprintf("xtol=%f is too small, no further improvement in the approximate\n  solution is possible.",
                       xtol),
         "8" = sprintf(paste0("gtol=%f is too small, func(x) is orthogonal to the columns of\n  the Jacobian to ",
                              "machine precision."), gtol))
}

# scipy.optimize.root(method = "df-sane"): the derivative free spectral
# method of La Cruz, Martinez and Raydan (2006), _spectral.py, with the
# nonmonotone line searches of _linesearch.py
.sp_df_sane <- function(func, x0, ftol = 1e-8, fatol = 0x1.56e1fc2f8f359p-997, maxfev = 1000L, fnorm = NULL, M = 10L,
                        eta_strategy = NULL, sigma_eps = 1e-10, sigma_0 = 1, line_search = "cruz") {
  if (!(line_search %in% c("cheng", "cruz")))
    stop("Invalid value '", line_search, "' for 'line_search'", call. = FALSE)
  fmerit <- function(F) sqrt(.spmp_ssum(F * F))^2
  nfev <- 1L
  F_k <- as.numeric(func(x0))
  x_k <- as.numeric(x0)
  f_k <- fmerit(F_k)
  f_0 <- f_k
  if (is.null(eta_strategy)) eta_strategy <- function(k, x, F) f_0 / (1 + k)^2
  if (is.null(fnorm)) fnorm <- function(F) f_k^(1 / 2)
  f <- function(x) {
    if (nfev >= maxfev) stop(structure(class = c("sp_no_convergence", "condition"), list(message = "", call = NULL)))
    nfev <<- nfev + 1L
    F <- as.numeric(func(x))
    list(f = fmerit(F), F = F)
  }
  k <- 0L
  sigma_k <- sigma_0
  F_0_norm <- fnorm(F_k)
  prev_fs <- f_k
  Q <- 1
  C <- f_0
  converged <- FALSE
  message <- "too many function evaluations required"
  gamma <- 1e-4; tau_min <- 0.1; tau_max <- 0.5; nu <- 0.85
  repeat {
    F_k_norm <- fnorm(F_k)
    if (F_k_norm < ftol * F_0_norm + fatol) {
      message <- "successful convergence"
      converged <- TRUE
      break
    }
    if (abs(sigma_k) > 1 / sigma_eps) sigma_k <- 1 / sigma_eps * sign(sigma_k)
    else if (abs(sigma_k) < sigma_eps) sigma_k <- sigma_eps
    d <- -sigma_k * F_k
    eta <- eta_strategy(k, x_k, F_k)
    # the nonmonotone line search of Cruz (reference f_bar, the largest of
    # the last M merits) or Cheng (the running average C)
    ref <- if (line_search == "cruz") max(prev_fs) else C
    fk_ls <- if (line_search == "cruz") prev_fs[length(prev_fs)] else f_k
    ls <- tryCatch({
      alpha_p <- 1; alpha_m <- 1
      repeat {
        xp <- x_k + alpha_p * d
        e <- f(xp)
        if (e$f <= ref + eta - gamma * alpha_p^2 * fk_ls) break
        alpha_tp <- alpha_p^2 * fk_ls / (e$f + (2 * alpha_p - 1) * fk_ls)
        xp <- x_k - alpha_m * d
        e <- f(xp)
        if (e$f <= ref + eta - gamma * alpha_m^2 * fk_ls) break
        alpha_tm <- alpha_m^2 * fk_ls / (e$f + (2 * alpha_m - 1) * fk_ls)
        alpha_p <- min(max(alpha_tp, tau_min * alpha_p), tau_max * alpha_p)
        alpha_m <- min(max(alpha_tm, tau_min * alpha_m), tau_max * alpha_m)
      }
      list(xp = xp, fp = e$f, Fp = e$F)
    }, sp_no_convergence = function(cnd) NULL)
    if (is.null(ls)) break
    if (line_search == "cheng") {
      Q_next <- nu * Q + 1
      C <- (nu * Q * (C + eta) + ls$fp) / Q_next
      Q <- Q_next
    }
    s_k <- ls$xp - x_k
    y_k <- ls$Fp - F_k
    sigma_k <- .spmp_ssum(s_k * s_k) / .spmp_ssum(s_k * y_k)
    x_k <- ls$xp
    F_k <- ls$Fp
    f_k <- ls$fp
    if (line_search == "cruz") prev_fs <- utils::tail(c(prev_fs, ls$fp), M)
    k <- k + 1L
  }
  list(x = x_k, success = converged, message = message, fun = F_k, nfev = nfev, nit = k, method = "df-sane")
}
