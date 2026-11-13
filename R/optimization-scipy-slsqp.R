# SLSQP of SciPy 1.18.1 (BSD 3-clause license): the driver
# scipy/optimize/_slsqp_py.py, the C translation scipy/optimize/src/slsqp.c
# of the Fortran code of Dieter Kraft (DFVLR) and the Lawson and Hanson NNLS
# of scipy/optimize/src/nnls.c. The LAPACK and BLAS routines that slsqp.c
# calls (dgerq2, dormr2, dgeqr2, dorm2r, dgelsy with dgeqp3, dlaqp2 and
# dlaic1, dlarfg, dlarfgp, dlarf, dlartgp, dlascl, dtrsm, dtrsv, dtpsv, dtpmv,
# dgemv and dger) follow the reference implementations of LAPACK 3.11 (BSD
# license) in their order of operations. SciPy's wheels link Accelerate,
# whose kernels round differently in the last bits (its ddot, for one, sums
# the products of small vectors in an order that depends on the length).
# The multiply adds that clang fuses in slsqp.c and nnls.c of the arm64
# wheel, and the fused daxpy of Accelerate, are fused here too (.spsq_fma).

.spsq_epsmach <- 2.220446049250313e-16
# DLAMCH('S'), DLAMCH('E') and DLAMCH('P')
.spsq_safmin <- 0x1.0000000000000p-1022
.spsq_lamch_e <- 0x1.0000000000000p-53
.spsq_lamch_p <- 2.220446049250313e-16

.spsq_exit_modes <- c("-1" = "Gradient evaluation required (g & a)",
                      "0" = "Optimization terminated successfully",
                      "1" = "Function evaluation required (f & c)",
                      "2" = "More equality constraints than independent variables",
                      "3" = "More than 3*n iterations in LSQ subproblem",
                      "4" = "Inequality constraints incompatible",
                      "5" = "Singular matrix E in LSQ subproblem",
                      "6" = "Singular matrix C in LSQ subproblem",
                      "7" = "Rank-deficient equality constraint subproblem HFTI",
                      "8" = "Positive directional derivative for linesearch",
                      "9" = "Iteration limit reached")

# nonzero in the sense of LAPACK's scans, where NaN .EQ. 0 is false
.spsq_nz <- function(v) is.na(v) | v != 0

# Fortran SIGN(a, b) as Accelerate evaluates it (the f2c d_sign): a
# negative zero b counts as positive
.spsq_sign <- function(a, b) if (isTRUE(b >= 0)) abs(a) else -abs(a)

# ---------------------------------------------------------------------------
# BLAS (reference order of operations)

# dtpsv('L', 'N', 'U'): solve L x = b with L unit lower triangular, packed
.spsq_dtpsv_lnu <- function(n, ap, x) {
  kk <- 1L
  for (j in seq_len(n)) {
    if (.spsq_nz(x[j]) && j < n) {
      rows <- (j + 1L):n
      x[rows] <- x[rows] - x[j] * ap[kk + seq_along(rows)]
    }
    kk <- kk + n - j + 1L
  }
  x
}

# dtpmv('L', 'T', 'U'): x := L' x
.spsq_dtpmv_ltu <- function(n, ap, x) {
  kk <- 1L
  for (j in seq_len(n)) {
    temp <- x[j]
    k <- kk + 1L
    if (j < n) for (i in (j + 1L):n) {
      temp <- temp + ap[k] * x[i]
      k <- k + 1L
    }
    x[j] <- temp
    kk <- kk + n - j + 1L
  }
  x
}

# dtpmv('L', 'N', 'U'): x := L x
.spsq_dtpmv_lnu <- function(n, ap, x) {
  kk <- (n * (n + 1L)) %/% 2L
  for (j in rev(seq_len(n))) {
    if (.spsq_nz(x[j]) && j < n) {
      temp <- x[j]
      k <- kk
      for (i in n:(j + 1L)) {
        x[i] <- x[i] + temp * ap[k]
        k <- k - 1L
      }
    }
    kk <- kk - (n - j + 1L)
  }
  x
}

# dtrsv('U', 'N', 'N'): solve A x = b, A upper triangular
.spsq_dtrsv_un <- function(A, x) {
  for (j in rev(seq_along(x))) {
    if (.spsq_nz(x[j])) {
      x[j] <- x[j] / A[j, j]
      if (j > 1L) {
        rows <- seq_len(j - 1L)
        x[rows] <- x[rows] - x[j] * A[rows, j]
      }
    }
  }
  x
}

# dtrsv('U', 'T', 'N'): solve A' x = b, A upper triangular
.spsq_dtrsv_ut <- function(A, x) {
  for (j in seq_along(x)) {
    temp <- x[j]
    for (i in seq_len(j - 1L)) temp <- temp - A[i, j] * x[i]
    x[j] <- temp / A[j, j]
  }
  x
}

# y := y + alpha A x (dgemv 'N' with beta = 1)
.spsq_gemv_n <- function(alpha, A, x, y) {
  for (j in seq_along(x)) y <- y + (alpha * x[j]) * A[, j]
  y
}

# y := beta y + alpha A' x (dgemv 'T'); beta is 0 or 1
.spsq_gemv_t <- function(alpha, A, x, beta, y) {
  if (nrow(A) == 0L) return(if (beta == 0) y * 0 else y)
  if (beta == 0) y[] <- 0
  for (j in seq_len(ncol(A))) y[j] <- y[j] + alpha * sum(A[, j] * x)
  y
}

# ---------------------------------------------------------------------------
# LAPACK

.spsq_dlapy2 <- function(x, y) {
  if (is.na(y)) return(y)
  if (is.na(x)) return(x)
  w <- max(abs(x), abs(y))
  z <- min(abs(x), abs(y))
  if (z == 0 || w > 1.7976931348623157e308) return(w)
  q <- z / w
  w * sqrt(1 + q * q)
}

# elementary reflector H with H (alpha, x) = (beta, 0)
.spsq_dlarfg <- function(n, alpha, x) {
  if (n <= 1L) return(list(alpha = alpha, x = x, tau = 0))
  xnorm <- .sptl_nrm2(x)
  if (isTRUE(xnorm == 0)) return(list(alpha = alpha, x = x, tau = 0))
  beta <- -.spsq_sign(.spsq_dlapy2(alpha, xnorm), alpha)
  safmin <- .spsq_safmin / .spsq_lamch_e
  knt <- 0L
  if (isTRUE(abs(beta) < safmin)) {
    rsafmn <- 1 / safmin
    repeat {
      knt <- knt + 1L
      x <- rsafmn * x
      beta <- beta * rsafmn
      alpha <- alpha * rsafmn
      if (!(abs(beta) < safmin && knt < 20L)) break
    }
    xnorm <- .sptl_nrm2(x)
    beta <- -.spsq_sign(.spsq_dlapy2(alpha, xnorm), alpha)
  }
  tau <- (beta - alpha) / beta
  x <- (1 / (alpha - beta)) * x
  for (j in seq_len(knt)) beta <- beta * safmin
  list(alpha = beta, x = x, tau = tau)
}

# dlarfg with a nonnegative beta
.spsq_dlarfgp <- function(n, alpha, x) {
  if (n <= 0L) return(list(alpha = alpha, x = x, tau = 0))
  xnorm <- .sptl_nrm2(x)
  if (isTRUE(xnorm <= .spsq_lamch_p * abs(alpha))) {
    if (isTRUE(alpha >= 0)) return(list(alpha = alpha, x = x, tau = 0))
    x[] <- 0
    return(list(alpha = -alpha, x = x, tau = 2))
  }
  beta <- .spsq_sign(.spsq_dlapy2(alpha, xnorm), alpha)
  smlnum <- .spsq_safmin / .spsq_lamch_e
  knt <- 0L
  if (isTRUE(abs(beta) < smlnum)) {
    bignum <- 1 / smlnum
    repeat {
      knt <- knt + 1L
      x <- bignum * x
      beta <- beta * bignum
      alpha <- alpha * bignum
      if (!(abs(beta) < smlnum && knt < 20L)) break
    }
    xnorm <- .sptl_nrm2(x)
    beta <- .spsq_sign(.spsq_dlapy2(alpha, xnorm), alpha)
  }
  savealpha <- alpha
  alpha <- alpha + beta
  if (isTRUE(beta < 0)) {
    beta <- -beta
    tau <- -alpha / beta
  } else {
    alpha <- xnorm * (xnorm / alpha)
    tau <- alpha / beta
    alpha <- -alpha
  }
  if (isTRUE(abs(tau) <= smlnum)) {
    if (isTRUE(savealpha >= 0)) {
      tau <- 0
    } else {
      tau <- 2
      x[] <- 0
      beta <- -savealpha
    }
  } else {
    x <- (1 / alpha) * x
  }
  for (j in seq_len(knt)) beta <- beta * smlnum
  list(alpha = beta, x = x, tau = tau)
}

# plane rotation with a nonnegative r
.spsq_dlartgp <- function(f, g) {
  safmin <- .spsq_safmin
  safmn2 <- 2^-484
  safmx2 <- 1 / safmn2
  if (isTRUE(g == 0)) return(list(cs = .spsq_sign(1, f), sn = 0, r = abs(f)))
  if (isTRUE(f == 0)) return(list(cs = 0, sn = .spsq_sign(1, g), r = abs(g)))
  f1 <- f
  g1 <- g
  scale <- max(abs(f1), abs(g1))
  if (isTRUE(scale >= safmx2)) {
    count <- 0L
    repeat {
      count <- count + 1L
      f1 <- f1 * safmn2
      g1 <- g1 * safmn2
      scale <- max(abs(f1), abs(g1))
      if (!(scale >= safmx2 && count < 20L)) break
    }
    r <- sqrt(f1 * f1 + g1 * g1)
    cs <- f1 / r
    sn <- g1 / r
    for (i in seq_len(count)) r <- r * safmx2
  } else if (isTRUE(scale <= safmn2)) {
    count <- 0L
    repeat {
      count <- count + 1L
      f1 <- f1 * safmx2
      g1 <- g1 * safmx2
      scale <- max(abs(f1), abs(g1))
      if (!(scale <= safmn2)) break
    }
    r <- sqrt(f1 * f1 + g1 * g1)
    cs <- f1 / r
    sn <- g1 / r
    for (i in seq_len(count)) r <- r * safmn2
  } else {
    r <- sqrt(f1 * f1 + g1 * g1)
    cs <- f1 / r
    sn <- g1 / r
  }
  if (isTRUE(r < 0)) {
    cs <- -cs
    sn <- -sn
    r <- -r
  }
  list(cs = cs, sn = sn, r = r)
}

# C := H C with H = I - tau v v' (dlarf 'Left'); C is a matrix
.spsq_dlarf_left <- function(v, tau, C) {
  if (isTRUE(tau == 0)) return(C)
  lastv <- length(v)
  while (lastv > 0L && isTRUE(v[lastv] == 0)) lastv <- lastv - 1L
  if (lastv == 0L) return(C)
  rows <- seq_len(lastv)
  lastc <- 0L
  for (j in rev(seq_len(ncol(C)))) if (any(.spsq_nz(C[rows, j]))) {
    lastc <- j
    break
  }
  vv <- v[rows]
  for (j in seq_len(lastc)) {
    wj <- sum(C[rows, j] * vv)
    if (.spsq_nz(wj)) C[rows, j] <- C[rows, j] + vv * (-tau * wj)
  }
  C
}

# C := C H with H = I - tau v v' (dlarf 'Right'); C is a matrix
.spsq_dlarf_right <- function(v, tau, C) {
  if (isTRUE(tau == 0)) return(C)
  lastv <- length(v)
  while (lastv > 0L && isTRUE(v[lastv] == 0)) lastv <- lastv - 1L
  if (lastv == 0L || nrow(C) == 0L) return(C)
  cols <- seq_len(lastv)
  lastc <- 0L
  for (j in cols) {
    nzr <- which(.spsq_nz(C[, j]))
    if (length(nzr)) lastc <- max(lastc, nzr[length(nzr)])
  }
  if (lastc == 0L) return(C)
  rows <- seq_len(lastc)
  w <- numeric(lastc)
  for (j in cols) w <- w + v[j] * C[rows, j]
  for (j in cols) if (.spsq_nz(v[j])) C[rows, j] <- C[rows, j] + w * (-tau * v[j])
  C
}

# unblocked QR factorization (dgeqr2)
.spsq_dgeqr2 <- function(A) {
  m <- nrow(A)
  n <- ncol(A)
  k <- min(m, n)
  tau <- numeric(k)
  for (i in seq_len(k)) {
    sub <- if (i < m) (i + 1L):m else integer()
    r <- .spsq_dlarfg(m - i + 1L, A[i, i], A[sub, i])
    A[i, i] <- r$alpha
    A[sub, i] <- r$x
    tau[i] <- r$tau
    if (i < n) {
      aii <- A[i, i]
      A[i, i] <- 1
      A[i:m, (i + 1L):n] <- .spsq_dlarf_left(A[i:m, i], tau[i], A[i:m, (i + 1L):n, drop = FALSE])
      A[i, i] <- aii
    }
  }
  list(A = A, tau = tau)
}

# C := Q' C with Q of k reflectors from dgeqr2 (dorm2r 'L', 'T')
.spsq_dorm2r_lt <- function(A, tau, C, k) {
  m <- nrow(C)
  if (m == 0L || ncol(C) == 0L || k == 0L) return(C)
  for (i in seq_len(k)) {
    v <- A[i:m, i]
    v[1L] <- 1
    C[i:m, ] <- .spsq_dlarf_left(v, tau[i], C[i:m, , drop = FALSE])
  }
  C
}

# unblocked RQ factorization (dgerq2)
.spsq_dgerq2 <- function(A) {
  m <- nrow(A)
  n <- ncol(A)
  k <- min(m, n)
  tau <- numeric(k)
  for (i in rev(seq_len(k))) {
    r0 <- m - k + i
    c0 <- n - k + i
    lead <- seq_len(c0 - 1L)
    r <- .spsq_dlarfg(c0, A[r0, c0], A[r0, lead])
    A[r0, c0] <- r$alpha
    A[r0, lead] <- r$x
    tau[i] <- r$tau
    if (r0 > 1L) {
      v <- A[r0, seq_len(c0)]
      v[c0] <- 1
      rows <- seq_len(r0 - 1L)
      A[rows, seq_len(c0)] <- .spsq_dlarf_right(v, tau[i], A[rows, seq_len(c0), drop = FALSE])
    }
  }
  list(A = A, tau = tau)
}

# C := Q' C (side "L") or C := C Q' (side "R") with Q from dgerq2 (dormr2)
.spsq_dormr2 <- function(side, A, tau, C) {
  k <- nrow(A)
  if (nrow(C) == 0L || ncol(C) == 0L || k == 0L) return(C)
  if (side == "R") {
    nq <- ncol(C)
    for (i in rev(seq_len(k))) {
      l <- nq - k + i
      v <- A[i, seq_len(l)]
      v[l] <- 1
      C[, seq_len(l)] <- .spsq_dlarf_right(v, tau[i], C[, seq_len(l), drop = FALSE])
    }
  } else {
    nq <- nrow(C)
    for (i in seq_len(k)) {
      l <- nq - k + i
      v <- A[i, seq_len(l)]
      v[l] <- 1
      C[seq_len(l), ] <- .spsq_dlarf_left(v, tau[i], C[seq_len(l), , drop = FALSE])
    }
  }
  C
}

# multiply A (type "G": all of it, "U": its upper triangle) by cto / cfrom
# without over or underflow (dlascl)
.spsq_dlascl <- function(type, cfrom, cto, A) {
  smlnum <- .spsq_safmin
  bignum <- 1 / smlnum
  cfromc <- cfrom
  ctoc <- cto
  upper <- if (type == "U") row(A) <= col(A) else matrix(TRUE, nrow(A), ncol(A))
  repeat {
    cfrom1 <- cfromc * smlnum
    if (cfrom1 == cfromc) {
      mul <- ctoc / cfromc
      done <- TRUE
    } else {
      cto1 <- ctoc / bignum
      if (cto1 == ctoc) {
        mul <- ctoc
        done <- TRUE
        cfromc <- 1
      } else if (abs(cfrom1) > abs(ctoc) && ctoc != 0) {
        mul <- smlnum
        done <- FALSE
        cfromc <- cfrom1
      } else if (abs(cto1) > abs(cfromc)) {
        mul <- bignum
        done <- FALSE
        ctoc <- cto1
      } else {
        mul <- ctoc / cfromc
        done <- TRUE
        if (mul == 1) return(A)
      }
    }
    A[upper] <- A[upper] * mul
    if (done) return(A)
  }
}

# incremental condition estimation (dlaic1): job 1 largest, job 2 smallest
# singular value
.spsq_dlaic1 <- function(job, x, sest, w, gamma) {
  eps <- .spsq_lamch_e
  alpha <- sum(x * w)
  absalp <- abs(alpha)
  absgam <- abs(gamma)
  absest <- abs(sest)
  if (isTRUE(job == 1L)) {
    if (isTRUE(sest == 0)) {
      s1 <- max(absgam, absalp)
      if (isTRUE(s1 == 0)) return(list(sestpr = 0, s = 0, c = 1))
      s <- alpha / s1
      c <- gamma / s1
      tmp <- sqrt(s * s + c * c)
      return(list(sestpr = s1 * tmp, s = s / tmp, c = c / tmp))
    }
    if (isTRUE(absgam <= eps * absest)) {
      tmp <- max(absest, absalp)
      s1 <- absest / tmp
      s2 <- absalp / tmp
      return(list(sestpr = tmp * sqrt(s1 * s1 + s2 * s2), s = 1, c = 0))
    }
    if (isTRUE(absalp <= eps * absest)) {
      if (isTRUE(absgam <= absest)) return(list(sestpr = absest, s = 1, c = 0))
      return(list(sestpr = absgam, s = 0, c = 1))
    }
    if (isTRUE(absest <= eps * absalp || absest <= eps * absgam)) {
      s1 <- absgam
      s2 <- absalp
      if (isTRUE(s1 <= s2)) {
        tmp <- s1 / s2
        s <- sqrt(1 + tmp * tmp)
        sestpr <- s2 * s
        c <- (gamma / s2) / s
        s <- .spsq_sign(1, alpha) / s
      } else {
        tmp <- s2 / s1
        c <- sqrt(1 + tmp * tmp)
        sestpr <- s1 * c
        s <- (alpha / s1) / c
        c <- .spsq_sign(1, gamma) / c
      }
      return(list(sestpr = sestpr, s = s, c = c))
    }
    zeta1 <- alpha / absest
    zeta2 <- gamma / absest
    b <- (1 - zeta1 * zeta1 - zeta2 * zeta2) * 0.5
    c <- zeta1 * zeta1
    t <- if (isTRUE(b > 0)) c / (b + sqrt(b * b + c)) else sqrt(b * b + c) - b
    sine <- -zeta1 / t
    cosine <- -zeta2 / (1 + t)
    tmp <- sqrt(sine * sine + cosine * cosine)
    return(list(sestpr = sqrt(t + 1) * absest, s = sine / tmp, c = cosine / tmp))
  }
  if (isTRUE(sest == 0)) {
    if (isTRUE(max(absgam, absalp) == 0)) {
      sine <- 1
      cosine <- 0
    } else {
      sine <- -gamma
      cosine <- alpha
    }
    s1 <- max(abs(sine), abs(cosine))
    s <- sine / s1
    c <- cosine / s1
    tmp <- sqrt(s * s + c * c)
    return(list(sestpr = 0, s = s / tmp, c = c / tmp))
  }
  if (isTRUE(absgam <= eps * absest)) return(list(sestpr = absgam, s = 0, c = 1))
  if (isTRUE(absalp <= eps * absest)) {
    if (isTRUE(absgam <= absest)) return(list(sestpr = absgam, s = 0, c = 1))
    return(list(sestpr = absest, s = 1, c = 0))
  }
  if (isTRUE(absest <= eps * absalp || absest <= eps * absgam)) {
    s1 <- absgam
    s2 <- absalp
    if (isTRUE(s1 <= s2)) {
      tmp <- s1 / s2
      c <- sqrt(1 + tmp * tmp)
      sestpr <- absest * (tmp / c)
      s <- -(gamma / s2) / c
      c <- .spsq_sign(1, alpha) / c
    } else {
      tmp <- s2 / s1
      s <- sqrt(1 + tmp * tmp)
      sestpr <- absest / s
      c <- (alpha / s1) / s
      s <- -.spsq_sign(1, gamma) / s
    }
    return(list(sestpr = sestpr, s = s, c = c))
  }
  zeta1 <- alpha / absest
  zeta2 <- gamma / absest
  norma <- max(1 + zeta1 * zeta1 + abs(zeta1 * zeta2), abs(zeta1 * zeta2) + zeta2 * zeta2)
  test <- 1 + 2 * (zeta1 - zeta2) * (zeta1 + zeta2)
  if (isTRUE(test >= 0)) {
    b <- (zeta1 * zeta1 + zeta2 * zeta2 + 1) * 0.5
    c <- zeta2 * zeta2
    t <- c / (b + sqrt(abs(b * b - c)))
    sine <- zeta1 / (1 - t)
    cosine <- -zeta2 / t
    sestpr <- sqrt(t + 4 * eps * eps * norma) * absest
  } else {
    b <- (zeta2 * zeta2 + zeta1 * zeta1 - 1) * 0.5
    c <- zeta1 * zeta1
    t <- if (isTRUE(b >= 0)) -c / (b + sqrt(b * b + c)) else b - sqrt(b * b + c)
    sine <- -zeta1 / t
    cosine <- -zeta2 / (1 + t)
    sestpr <- sqrt(1 + t + 4 * eps * eps * norma) * absest
  }
  tmp <- sqrt(sine * sine + cosine * cosine)
  list(sestpr = sestpr, s = sine / tmp, c = cosine / tmp)
}

# QR factorization with column pivoting (dgeqp3): columns with jpvt != 0
# are fixed in front, the free ones are pivoted by dlaqp2. dgeqrf and dormqr
# run unblocked as LAPACK does for fewer than 128 (dgeqrf) and 32 (dormqr)
# columns; dgeqp3 switches to the blocked dlaqps beyond 32 free columns.
# ponytail: unblocked for every size, port dlaqps if SLSQP ever needs more
# than 32 variables to agree with LAPACK bit for bit.
.spsq_dgeqp3 <- function(A, jpvt) {
  m <- nrow(A)
  n <- ncol(A)
  minmn <- min(m, n)
  nfxd <- 1L
  for (j in seq_len(n)) {
    if (jpvt[j] != 0L) {
      if (j != nfxd) {
        tmp <- A[, j]
        A[, j] <- A[, nfxd]
        A[, nfxd] <- tmp
        jpvt[j] <- jpvt[nfxd]
        jpvt[nfxd] <- j
      } else {
        jpvt[j] <- j
      }
      nfxd <- nfxd + 1L
    } else {
      jpvt[j] <- j
    }
  }
  nfxd <- nfxd - 1L
  tau <- numeric(minmn)
  if (nfxd > 0L) {
    na <- min(m, nfxd)
    qr <- .spsq_dgeqr2(A[, seq_len(na), drop = FALSE])
    A[, seq_len(na)] <- qr$A
    tau[seq_len(na)] <- qr$tau
    if (na < n) {
      rest <- (na + 1L):n
      A[, rest] <- .spsq_dorm2r_lt(qr$A, qr$tau, A[, rest, drop = FALSE], na)
    }
  }
  if (nfxd < minmn) {
    vn1 <- numeric(n)
    for (j in (nfxd + 1L):n) vn1[j] <- .sptl_nrm2(A[(nfxd + 1L):m, j])
    vn2 <- vn1
    tol3z <- sqrt(.spsq_lamch_e)
    for (i in (nfxd + 1L):minmn) {
      # dlaqp2 with offset nfxd, written in the full column numbering
      offpi <- i
      cand <- abs(vn1[i:n])
      pvt <- i
      dmax <- cand[1L]
      for (q in seq_along(cand)[-1L]) if (isTRUE(cand[q] > dmax)) {
        dmax <- cand[q]
        pvt <- i + q - 1L
      }
      if (pvt != i) {
        tmp <- A[, pvt]
        A[, pvt] <- A[, i]
        A[, i] <- tmp
        itemp <- jpvt[pvt]
        jpvt[pvt] <- jpvt[i]
        jpvt[i] <- itemp
        vn1[pvt] <- vn1[i]
        vn2[pvt] <- vn2[i]
      }
      if (offpi < m) {
        sub <- (offpi + 1L):m
        r <- .spsq_dlarfg(m - offpi + 1L, A[offpi, i], A[sub, i])
        A[offpi, i] <- r$alpha
        A[sub, i] <- r$x
      } else {
        r <- .spsq_dlarfg(1L, A[m, i], numeric())
        A[m, i] <- r$alpha
      }
      tau[i] <- r$tau
      if (i < n) {
        aii <- A[offpi, i]
        A[offpi, i] <- 1
        A[offpi:m, (i + 1L):n] <- .spsq_dlarf_left(A[offpi:m, i], tau[i], A[offpi:m, (i + 1L):n, drop = FALSE])
        A[offpi, i] <- aii
        for (j in (i + 1L):n) {
          if (.spsq_nz(vn1[j])) {
            q <- abs(A[offpi, j]) / vn1[j]
            temp <- max(1 - q * q, 0)
            q <- vn1[j] / vn2[j]
            temp2 <- temp * (q * q)
            if (isTRUE(temp2 <= tol3z)) {
              if (offpi < m) {
                vn1[j] <- .sptl_nrm2(A[(offpi + 1L):m, j])
                vn2[j] <- vn1[j]
              } else {
                vn1[j] <- 0
                vn2[j] <- 0
              }
            } else {
              vn1[j] <- vn1[j] * sqrt(temp)
            }
          }
        }
      }
    }
  }
  list(A = A, jpvt = jpvt, tau = tau)
}

# minimum norm least squares by complete orthogonal factorization (dgelsy)
# for one right hand side b (length nrow(A) >= ncol(A)). A rank deficient
# result is only reported: SLSQP stops with mode 7 and discards the step, so
# the dtzrzf and dormrz stages are not ported.
.spsq_dgelsy <- function(A, b, jpvt, rcond) {
  m <- nrow(A)
  n <- ncol(A)
  mn <- min(m, n)
  if (mn == 0L) return(list(b = b, jpvt = jpvt, rank = 0L))
  smlnum <- .spsq_safmin / .spsq_lamch_p
  bignum <- 1 / smlnum
  anrm <- max(abs(A))
  iascl <- 0L
  if (isTRUE(anrm > 0 && anrm < smlnum)) {
    A <- .spsq_dlascl("G", anrm, smlnum, A)
    iascl <- 1L
  } else if (isTRUE(anrm > bignum)) {
    A <- .spsq_dlascl("G", anrm, bignum, A)
    iascl <- 2L
  } else if (isTRUE(anrm == 0)) {
    return(list(b = b * 0, jpvt = jpvt, rank = 0L))
  }
  bm <- matrix(b, ncol = 1L)
  bnrm <- max(abs(b))
  ibscl <- 0L
  if (isTRUE(bnrm > 0 && bnrm < smlnum)) {
    bm <- .spsq_dlascl("G", bnrm, smlnum, bm)
    ibscl <- 1L
  } else if (isTRUE(bnrm > bignum)) {
    bm <- .spsq_dlascl("G", bnrm, bignum, bm)
    ibscl <- 2L
  }
  qp <- .spsq_dgeqp3(A, jpvt)
  A <- qp$A
  jpvt <- qp$jpvt
  xmin <- 1
  xmax <- 1
  smax <- abs(A[1L, 1L])
  smin <- smax
  if (isTRUE(abs(A[1L, 1L]) == 0)) return(list(b = b * 0, jpvt = jpvt, rank = 0L))
  rank <- 1L
  while (rank < mn) {
    i <- rank + 1L
    lo <- .spsq_dlaic1(2L, xmin, smin, A[seq_len(rank), i], A[i, i])
    hi <- .spsq_dlaic1(1L, xmax, smax, A[seq_len(rank), i], A[i, i])
    if (!isTRUE(hi$sestpr * rcond <= lo$sestpr)) break
    xmin <- c(lo$s * xmin, lo$c)
    xmax <- c(hi$s * xmax, hi$c)
    smin <- lo$sestpr
    smax <- hi$sestpr
    rank <- rank + 1L
  }
  if (rank < n) return(list(b = b, jpvt = jpvt, rank = rank))
  bm <- .spsq_dorm2r_lt(A, qp$tau, bm, mn)
  # dtrsm('L', 'U', 'N', 'N') for one column
  for (k in rev(seq_len(rank))) {
    if (.spsq_nz(bm[k, 1L])) {
      bm[k, 1L] <- bm[k, 1L] / A[k, k]
      if (k > 1L) {
        rows <- seq_len(k - 1L)
        bm[rows, 1L] <- bm[rows, 1L] - bm[k, 1L] * A[rows, k]
      }
    }
  }
  work <- numeric(n)
  work[jpvt] <- bm[seq_len(n), 1L]
  bm[seq_len(n), 1L] <- work
  top <- matrix(bm[seq_len(n), 1L], ncol = 1L)
  if (iascl == 1L) top <- .spsq_dlascl("G", anrm, smlnum, top)
  if (iascl == 2L) top <- .spsq_dlascl("G", anrm, bignum, top)
  if (ibscl == 1L) top <- .spsq_dlascl("G", smlnum, bnrm, top)
  if (ibscl == 2L) top <- .spsq_dlascl("G", bignum, bnrm, top)
  bm[seq_len(n), 1L] <- top
  list(b = bm[, 1L], jpvt = jpvt, rank = rank)
}

# ---------------------------------------------------------------------------
# NNLS of Lawson and Hanson (nnls.c): min |A x - b| subject to x >= 0. The
# column order in `idx` replaces the 0 based `indices` array of the C code.

# distance from u > 0 to the next larger double
.spsq_spacing <- function(u) {
  if (!is.finite(u)) return(NaN)
  e <- floor(log2(u))
  if (2^e > u) e <- e - 1
  if (2^(e + 1) <= u) e <- e + 1
  max(2^(e - 52), 2^-1074)
}

.spsq_nnls <- function(a, b, maxiter) {
  m <- nrow(a)
  n <- ncol(a)
  if (m <= 0L || n <= 0L) return(list(x = numeric(max(n, 0L)), rnorm = 0, info = 2L, idx = seq_len(max(n, 0L))))
  info <- 1L
  idx <- seq_len(n)
  x <- numeric(n)
  w <- numeric(n)
  zz <- numeric(m)
  indz <- 0L
  iteration <- 0L
  izmax <- 0L
  finished <- FALSE
  # back substitution of the triangular system of the columns in set P
  tri_solve <- function(zz, a, idx, indz) {
    jj <- 0L
    for (k in seq_len(indz) - 1L) {
      ip <- indz - k
      if (k != 0L) {
        rows <- seq_len(ip)
        zz[rows] <- .spsq_fma(-a[rows, jj], zz[ip + 1L], zz[rows])
      }
      jj <- idx[ip]
      zz[ip] <- zz[ip] / a[ip, jj]
    }
    zz
  }
  while (indz < min(m, n)) {
    # dual vector of the columns in set Z
    zrows <- (indz + 1L):m
    for (i in (indz + 1L):n) {
      j <- idx[i]
      w[j] <- sum(a[zrows, j] * b[zrows])
    }
    repeat {
      wmax <- 0
      for (k in (indz + 1L):n) {
        j <- idx[k]
        if (isTRUE(w[j] > wmax)) {
          wmax <- w[j]
          izmax <- k
        }
      }
      if (!(wmax > 0)) {
        finished <- TRUE
        break
      }
      iz <- izmax
      j <- idx[iz]
      # dlarfgp scales the subdiagonal of column j in place even when the
      # column is then rejected, as in nnls.c
      sub <- if (indz + 2L <= m) (indz + 2L):m else integer()
      r <- .spsq_dlarfgp(m - indz, a[indz + 1L, j], a[sub, j])
      pivot <- r$alpha
      a[sub, j] <- r$x
      tau <- r$tau
      unorm <- if (indz > 0L) .sptl_nrm2(a[seq_len(indz), j]) else 0
      spacing <- if (isTRUE(unorm > 0)) .spsq_spacing(unorm) else 0
      if (isTRUE(abs(pivot) > 100 * spacing)) {
        zz <- b
        pivot2 <- a[indz + 1L, j]
        a[indz + 1L, j] <- 1
        zz[zrows] <- .spsq_dlarf_left(a[zrows, j], tau, matrix(zz[zrows], ncol = 1L))[, 1L]
        ztest <- zz[indz + 1L] / pivot
        if (isTRUE(ztest > 0)) break
        a[indz + 1L, j] <- pivot2
      }
      w[j] <- 0
    }
    if (finished) break
    # move column j from set Z to set P
    b <- zz
    idx[iz] <- idx[indz + 1L]
    idx[indz + 1L] <- j
    indz <- indz + 1L
    if (indz < n) {
      prow <- indz:m
      v <- a[prow, j]
      for (k in (indz + 1L):n) {
        jj <- idx[k]
        a[prow, jj] <- .spsq_dlarf_left(v, tau, matrix(a[prow, jj], ncol = 1L))[, 1L]
      }
    }
    a[indz, j] <- pivot
    if (indz < m) a[(indz + 1L):m, j] <- 0
    w[j] <- 0
    zz <- tri_solve(zz, a, idx, indz)
    repeat {
      iteration <- iteration + 1L
      if (iteration >= maxiter) {
        info <- 3L
        finished <- TRUE
        break
      }
      # step toward zz as far as the coefficients in set P stay feasible
      alpha <- 2
      jpos <- 0L
      for (ip in seq_len(indz)) {
        k <- idx[ip]
        if (isTRUE(zz[ip] <= 0)) {
          tv <- -x[k] / (zz[ip] - x[k])
          if (isTRUE(alpha > tv)) {
            alpha <- tv
            jpos <- ip
          }
        }
      }
      if (alpha == 2) break
      for (ip in seq_len(indz)) {
        k <- idx[ip]
        x[k] <- .spsq_fma(alpha, zz[ip] - x[k], x[k])
      }
      # move the infeasible coefficients from set P back to set Z
      i <- idx[jpos]
      repeat {
        x[i] <- 0
        if (jpos != indz) {
          jpos <- jpos + 1L
          for (jn in jpos:indz) {
            ii <- idx[jn]
            idx[jn - 1L] <- ii
            g <- .spsq_dlartgp(a[jn - 1L, ii], a[jn, ii])
            a[jn - 1L, ii] <- g$r
            a[jn, ii] <- 0
            ks <- seq_len(n)[-ii]
            tmp <- a[jn - 1L, ks]
            a[jn - 1L, ks] <- .spsq_fma(g$cs, tmp, g$sn * a[jn, ks])
            a[jn, ks] <- .spsq_fma(-g$sn, tmp, g$cs * a[jn, ks])
            tmp <- b[jn - 1L]
            b[jn - 1L] <- .spsq_fma(g$cs, tmp, g$sn * b[jn])
            b[jn] <- .spsq_fma(-g$sn, tmp, g$cs * b[jn])
          }
        }
        indz <- indz - 1L
        idx[indz + 1L] <- i
        # nnls.c would index before its arrays once set P is empty; that
        # cannot happen in exact arithmetic and is stopped here
        if (indz == 0L) break
        allpos <- TRUE
        for (jpos in seq_len(indz)) {
          i <- idx[jpos]
          if (isTRUE(x[i] <= 0)) {
            allpos <- FALSE
            break
          }
        }
        if (allpos) break
      }
      zz <- tri_solve(b, a, idx, indz)
    }
    if (finished) break
    for (k in seq_len(indz)) x[idx[k]] <- zz[k]
  }
  rnorm <- if (indz < m) .sptl_nrm2(b[(indz + 1L):m]) else 0
  list(x = x, rnorm = rnorm, info = info, idx = idx)
}

# ---------------------------------------------------------------------------
# the subproblems of slsqp.c

# least distance programming: min |x| subject to G x >= h (ldp)
.spsq_ldp <- function(G, h, jw) {
  m <- nrow(G)
  n <- ncol(G)
  if (n <= 0L) return(list(x = numeric(0), mode = 2L, mults = numeric(m), jw = jw))
  x <- numeric(n)
  if (m == 0L) return(list(x = x, mode = 1L, mults = numeric(0), jw = jw))
  # the dual problem [G'; h'] y = [0; 1], y >= 0
  r <- .spsq_nnls(rbind(t(G), h), c(numeric(n), 1), 3L * m)
  jw[seq_len(m)] <- r$idx - 1L
  if (r$info != 1L) return(list(x = x, mode = r$info, mults = numeric(m), jw = jw))
  if (isTRUE(r$rnorm <= 0)) return(list(x = x, mode = 4L, mults = numeric(m), jw = jw))
  y <- r$x
  fac <- 1 - sum(h * y)
  if (!isTRUE((1 + fac) - 1 > 0)) return(list(x = x, mode = 4L, mults = numeric(m), jw = jw))
  fac <- 1 / fac
  for (j in seq_len(n)) x[j] <- fac * sum(G[, j] * y)
  list(x = x, mode = 1L, mults = fac * y, jw = jw)
}

# min |A x - b| subject to G x >= h (lsi); A is ma x n with ma >= n
.spsq_lsi <- function(A, b, G, h, jw) {
  n <- ncol(A)
  k <- min(nrow(A), n)
  qr <- .spsq_dgeqr2(A)
  A <- qr$A
  b <- .spsq_dorm2r_lt(A, qr$tau, matrix(b, ncol = 1L), k)[, 1L]
  for (i in seq_len(k)) if (!isTRUE(abs(A[i, i]) >= .spsq_epsmach))
    return(list(x = numeric(n), mode = 5L, mults = numeric(nrow(G)), jw = jw))
  # G := G R^-1 (dtrsm 'R', 'U', 'N', 'N')
  for (j in seq_len(n)) {
    for (kk in seq_len(j - 1L)) if (.spsq_nz(A[kk, j])) G[, j] <- G[, j] - A[kk, j] * G[, kk]
    temp <- 1 / A[j, j]
    G[, j] <- temp * G[, j]
  }
  h <- .spsq_gemv_n(-1, G, b[seq_len(n)], h)
  r <- .spsq_ldp(G, h, jw)
  if (r$mode != 1L) return(r)
  x <- .spsq_dtrsv_un(A[seq_len(n), seq_len(n), drop = FALSE], r$x + b[seq_len(n)])
  list(x = x, mode = 1L, mults = r$mults, jw = r$jw)
}

# min |A x - b| subject to E x = f and G x >= h (lsei)
.spsq_lsei <- function(A, b, E, f, G, h, jw) {
  n <- ncol(A)
  me <- nrow(E)
  mg <- nrow(G)
  x <- numeric(n)
  out <- function(mode, gm = numeric(mg), em = numeric(me)) list(x = x, mode = mode, gmults = gm, emults = em, jw = jw)
  if (me > n) return(out(2L))
  nvars <- n - me
  eq <- nvars + seq_len(me)
  tau <- numeric(0)
  if (me > 0L) {
    rq <- .spsq_dgerq2(E)
    E <- rq$A
    tau <- rq$tau
    A <- .spsq_dormr2("R", E, tau, A)
    G <- .spsq_dormr2("R", E, tau, G)
    for (i in seq_len(me)) if (!isTRUE(abs(E[i, nvars + i]) >= .spsq_epsmach)) return(out(6L))
    x[eq] <- .spsq_dtrsv_un(E[, eq, drop = FALSE], f)
  }
  mode <- 1L
  gmults <- numeric(mg)
  if (me < n) {
    free <- seq_len(nvars)
    wb <- .spsq_gemv_n(-1, A[, eq, drop = FALSE], x[eq], b)
    if (mg == 0L) {
      gl <- .spsq_dgelsy(A[, free, drop = FALSE], wb, jw[free], sqrt(.spsq_epsmach))
      jw[free] <- gl$jpvt
      x[free] <- gl$b[free]
      if (gl$rank < nvars) return(out(7L))
    } else {
      h <- .spsq_gemv_n(-1, G[, eq, drop = FALSE], x[eq], h)
      r <- .spsq_lsi(A[, free, drop = FALSE], wb, G[, free, drop = FALSE], h, jw)
      x[free] <- r$x
      mode <- r$mode
      gmults <- r$mults
      jw <- r$jw
      if (me == 0L || mode != 1L) return(out(mode, gmults))
    }
  }
  if (me == 0L) return(out(mode, gmults))
  # back to the original basis: residual, then the equality multipliers
  b <- .spsq_gemv_n(1, A, x, -b)
  f <- .spsq_gemv_t(1, A[, eq, drop = FALSE], b, 0, numeric(me))
  f <- .spsq_gemv_t(-1, G[, eq, drop = FALSE], gmults, 1, f)
  x <- .spsq_dormr2("L", E, tau, matrix(x, ncol = 1L))[, 1L]
  emults <- .spsq_dtrsv_ut(E[, eq, drop = FALSE], f)
  out(mode, gmults, emults)
}

# the QP subproblem from the LDL' factor Lf of the quasi Newton matrix (lsq);
# `augment` adds the slack variable of an inconsistent linearization
.spsq_lsq <- function(m, meq, n, augment, aug_weight, Lf, gradx, C, d, xl, xu, x, y, jw) {
  orign <- n
  mineq <- m - meq
  ld <- n
  if (augment) {
    ld <- n + 1L
    x[n + 1L] <- 1
    xl[n + 1L] <- 0
    xu[n + 1L] <- 1
  }
  # A = sqrt(D) L' and b = -sqrt(D)^-1 L^-1 gradx
  wA <- matrix(0, ld, ld)
  cursor <- 0L
  for (j in seq_len(n)) {
    diag <- sqrt(Lf[cursor + 1L])
    seg <- Lf[cursor + seq_len(n - j + 1L)]
    wA[j, j:n] <- c(diag, seg[-1L] * diag)
    cursor <- cursor + n - j + 1L
  }
  wb <- numeric(ld)
  wb[seq_len(n)] <- .spsq_dtpsv_lnu(n, Lf, gradx)
  cursor <- 1L
  for (i in seq_len(n)) {
    wb[i] <- wb[i] / -sqrt(Lf[cursor])
    cursor <- cursor + n - i + 1L
  }
  if (augment) {
    wA[ld, ld] <- aug_weight
    n <- n + 1L
  }
  wE <- matrix(0, meq, n)
  wf <- numeric(meq)
  if (meq > 0L) {
    wE[, seq_len(orign)] <- C[seq_len(meq), , drop = FALSE]
    if (augment) wE[, n] <- -d[seq_len(meq)]
    wf <- -d[seq_len(meq)]
  }
  # inequalities, then the finite lower and upper bounds (NaN marks none)
  ineq <- meq + seq_len(mineq)
  low <- which(!is.na(xl[seq_len(n)]))
  upp <- which(!is.na(xu[seq_len(n)]))
  nrows <- mineq + length(low) + length(upp)
  wh <- c(-d[ineq], xl[low], -xu[upp])
  wG <- matrix(0, nrows, n)
  if (mineq > 0L) {
    wG[seq_len(mineq), seq_len(orign)] <- C[ineq, , drop = FALSE]
    if (augment) {
      v <- -d[ineq]
      wG[seq_len(mineq), orign + 1L] <- ifelse(is.na(v), 0, pmax(v, 0))
    }
  }
  wG[cbind(mineq + seq_along(low), low)] <- 1
  wG[cbind(mineq + length(low) + seq_along(upp), upp)] <- -1
  r <- .spsq_lsei(wA, wb, wE, wf, wG, wh, jw)
  xs <- r$x
  if (r$mode == 1L) {
    y[seq_len(meq)] <- r$emults
    y[meq + seq_len(mineq)] <- r$gmults[seq_len(mineq)]
    y[m + seq_len(2L * n)] <- NaN
  }
  # clamp the solution to the finite bounds
  for (i in seq_len(n)) {
    if (!is.na(xl[i]) && isTRUE(xs[i] < xl[i])) xs[i] <- xl[i]
    else if (!is.na(xu[i]) && isTRUE(xs[i] > xu[i])) xs[i] <- xu[i]
  }
  x[seq_len(n)] <- xs
  list(x = x, y = y, xl = xl, xu = xu, mode = r$mode, jw = r$jw)
}

# rank one update of the LDL' factor a (packed) by sigma z z' with the
# composite t method of Fletcher and Powell (ldl_update)
.spsq_ldl_update <- function(n, a, z, sigma, w) {
  if (isTRUE(sigma == 0)) return(list(a = a, z = z, w = w))
  ij <- 1L
  t <- 1 / sigma
  if (isTRUE(sigma <= 0)) {
    w[seq_len(n)] <- z[seq_len(n)]
    for (i in seq_len(n)) {
      v <- w[i]
      t <- t + v * v / a[ij]
      if (i < n) {
        js <- (i + 1L):n
        w[js] <- .spsq_fma(-v, a[ij + seq_along(js)], w[js])
      }
      ij <- ij + n - i + 1L
    }
    if (isTRUE(t >= 0)) t <- .spsq_epsmach / sigma
    for (i in seq_len(n) - 1L) {
      j <- n - i
      ij <- ij - (i + 1L)
      u <- w[j]
      w[j] <- t
      t <- t - u * u / a[ij]
    }
  }
  for (i in seq_len(n)) {
    v <- z[i]
    delta <- v / a[ij]
    tp <- if (isTRUE(sigma < 0)) w[i] else .spsq_fma(delta, v, t)
    alpha <- tp / t
    a[ij] <- alpha * a[ij]
    if (i == n) break
    beta <- delta / tp
    js <- (i + 1L):n
    pos <- ij + seq_along(js)
    if (isTRUE(alpha <= 4)) {
      z[js] <- .spsq_fma(-v, a[pos], z[js])
      a[pos] <- .spsq_fma(beta, z[js], a[pos])
    } else {
      gamma <- t / tp
      u <- a[pos]
      a[pos] <- .spsq_fma(gamma, u, beta * z[js])
      z[js] <- .spsq_fma(-v, u, z[js])
    }
    ij <- ij + n - i + 1L
    t <- tp
  }
  list(a = a, z = z, w = w)
}

# a * b + c rounded once: the multiply adds that clang contracts in
# slsqp.c and nnls.c (the sites of the arm64 wheel) and the daxpy of
# Accelerate; an infinite product or addend passes through as in hardware
.spsq_fma <- function(a, b, c) {
  p <- a * b
  ifelse(is.finite(p) & is.finite(c), .fma(a, b, c), p + c)
}

# fmax of C: a NaN argument is ignored
.spsq_fmax <- function(a, b) ifelse(is.na(a), b, ifelse(is.na(b), a, pmax(a, b)))

# one reverse communication step of __slsqp_body; S holds the state of the
# C struct and of the persistent part of its buffer
.spsq_body <- function(S, funx, gradx, C, d, xl, xu) {
  n <- S$n
  m <- S$m
  meq <- S$meq
  eqmask <- seq_len(m) <= meq
  alfmin <- 0.1
  diagpos <- cumsum(c(1L, n - seq_len(n - 1L) + 1L))
  # the badlin flag is local in slsqp.c, so it is 0 on every new call
  badlin <- 0L
  label <- if (S$mode == -1L) "MODEM1" else if (S$mode == 1L) "MODE1" else "MODE0"
  viol_sum <- function() {
    h <- 0
    for (j in seq_len(m)) h <- h + .spsq_fmax(-d[j], if (eqmask[j]) d[j] else 0)
    h
  }
  repeat {
    switch(label,
      MODE0 = {
        S$acc <- abs(S$acc)
        S$tol <- 10 * S$acc
        S$iter <- 0L
        S$reset <- 0L
        S$s[seq_len(n)] <- 0
        S$mu[] <- 0
        label <- "RESET_BFGS"
      },
      RESET_BFGS = {
        S$reset <- S$reset + 1L
        if (S$reset > 5L) {
          label <- "LABEL255"
        } else {
          S$bfgs[] <- 0
          S$bfgs[diagpos] <- 1
          label <- "ITER_START"
        }
      },
      ITER_START = {
        S$mode <- 9L
        if (S$iter >= S$itermax) return(invisible())
        S$iter <- S$iter + 1L
        S$u[seq_len(n)] <- -S$sol + xl
        S$v[seq_len(n)] <- -S$sol + xu
        S$h4 <- 1
        run_lsq <- function(augment, rho) {
          r <- .spsq_lsq(m, meq, n, augment, rho, S$bfgs, gradx, C, d, S$u, S$v, S$s, S$mult, S$jw)
          S$s <- r$x
          S$mult <- r$y
          S$u <- r$xl
          S$v <- r$xu
          S$jw <- r$jw
          S$mode <- r$mode
        }
        run_lsq(FALSE, 0)
        badlin <- 0L
        if (S$mode == 6L && n == meq) S$mode <- 4L
        if (S$mode == 4L) {
          # inconsistent linearization: the augmented problem of Kraft
          # (1988), section 2.2.3
          badlin <- 1L
          S$s[seq_len(n)] <- 0
          S$h3 <- 0
          rho <- 100
          S$inconsistent <- 0L
          repeat {
            run_lsq(TRUE, rho)
            S$h4 <- 1 - S$s[n + 1L]
            if (S$mode == 4L) {
              rho <- rho * 10
              S$inconsistent <- S$inconsistent + 1L
              if (S$inconsistent > 5L) return(invisible())
              next
            } else if (S$mode != 1L) {
              return(invisible())
            }
            break
          }
        } else if (S$mode != 1L) {
          return(invisible())
        }
        mult <- S$mult[seq_len(m)]
        S$v[seq_len(n)] <- .spsq_gemv_t(-1, C, mult, 1, gradx)
        S$f0 <- funx
        S$x0 <- S$sol
        s <- S$s[seq_len(n)]
        S$gs <- sum(gradx * s)
        S$h1 <- abs(S$gs)
        S$h2 <- 0
        for (j in seq_len(m)) {
          S$h3 <- if (eqmask[j]) d[j] else 0
          S$h2 <- S$h2 + .spsq_fmax(-d[j], S$h3)
          S$h3 <- abs(mult[j])
          S$mu[j] <- .spsq_fmax(S$h3, (S$mu[j] + S$h3) / 2)
          S$h1 <- .spsq_fma(S$h3, abs(d[j]), S$h1)
        }
        S$mode <- 0L
        if (isTRUE(S$h1 < S$acc) && isTRUE(S$h2 < S$acc) && !badlin && isTRUE(funx == funx)) return(invisible())
        S$h1 <- 0
        for (j in seq_len(m)) {
          S$h3 <- if (eqmask[j]) d[j] else 0
          S$h1 <- .spsq_fma(S$mu[j], .spsq_fmax(-d[j], S$h3), S$h1)
        }
        S$t0 <- funx + S$h1
        S$h3 <- .spsq_fma(-S$h1, S$h4, S$gs)
        S$mode <- 8L
        if (isTRUE(S$h3 >= 0)) {
          label <- "RESET_BFGS"
        } else {
          S$line <- 0L
          S$alpha <- 1
          label <- "LINE_SEARCH"
        }
      },
      LINE_SEARCH = {
        S$line <- S$line + 1L
        S$h3 <- S$alpha * S$h3
        S$s[seq_len(n)] <- S$alpha * S$s[seq_len(n)]
        S$sol <- S$x0 + S$s[seq_len(n)]
        S$mode <- 1L
        return(invisible())
      },
      MODE1 = {
        S$t <- funx
        for (j in seq_len(m)) {
          S$h1 <- if (eqmask[j]) d[j] else 0
          S$t <- .spsq_fma(S$mu[j], .spsq_fmax(-d[j], S$h1), S$t)
        }
        S$h1 <- S$t - S$t0
        if (!(isTRUE(S$h1 <= S$h3 / 10) || S$line > 10L)) {
          S$alpha <- .spsq_fmax(S$h3 / (2 * (S$h3 - S$h1)), alfmin)
          label <- "LINE_SEARCH"
        } else {
          S$h3 <- 0
          for (j in seq_len(m)) {
            S$h1 <- if (eqmask[j]) d[j] else 0
            S$h3 <- S$h3 + .spsq_fmax(-d[j], S$h1)
          }
          small_step <- isTRUE(abs(funx - S$f0) < S$acc) || isTRUE(.sptl_nrm2(S$s[seq_len(n)]) < S$acc)
          S$mode <- if (small_step && isTRUE(S$h3 < S$acc) && !badlin && isTRUE(funx == funx)) 0L else -1L
          return(invisible())
        }
      },
      LABEL255 = {
        # relaxed convergence test after a positive directional derivative
        S$h3 <- 0
        for (j in seq_len(m)) {
          S$h1 <- if (eqmask[j]) d[j] else 0
          S$h3 <- S$h3 + .spsq_fmax(-d[j], S$h1)
        }
        small_step <- isTRUE(abs(funx - S$f0) < S$tol) || isTRUE(.sptl_nrm2(S$s[seq_len(n)]) < S$tol)
        S$mode <- if (small_step && isTRUE(S$h3 < S$tol) && !badlin && isTRUE(funx == funx)) 0L else 8L
        return(invisible())
      },
      MODEM1 = {
        # BFGS update of the LDL' factor, damped by Powell's rule
        mult <- S$mult[seq_len(m)]
        u <- .spsq_gemv_t(-1, C, mult, 1, gradx) - S$v[seq_len(n)]
        s <- S$s[seq_len(n)]
        v <- .spsq_dtpmv_ltu(n, S$bfgs, s)
        v <- S$bfgs[diagpos] * v
        v <- .spsq_dtpmv_lnu(n, S$bfgs, v)
        S$h1 <- sum(s * u)
        S$h2 <- sum(s * v)
        S$h3 <- 0.2 * S$h2
        if (isTRUE(S$h1 < S$h3)) {
          S$h4 <- (S$h2 - S$h3) / (S$h2 - S$h1)
          S$h1 <- S$h3
          tmp <- 1 - S$h4
          u <- S$h4 * u
          # the daxpy of Accelerate fuses the multiply and the add
          if (.spsq_nz(tmp)) u <- .spsq_fma(tmp, v, u)
        }
        S$u[seq_len(n)] <- u
        S$v[seq_len(n)] <- v
        if (isTRUE(S$h1 == 0) || isTRUE(S$h2 == 0)) {
          label <- "RESET_BFGS"
        } else {
          r1 <- .spsq_ldl_update(n, S$bfgs, u, 1 / S$h1, v)
          r2 <- .spsq_ldl_update(n, r1$a, r1$w, -1 / S$h2, r1$z)
          S$bfgs <- r2$a
          S$v[seq_len(n)] <- r2$z
          S$u[seq_len(n)] <- r2$w
          label <- "ITER_START"
        }
      })
  }
}

# Python's repr of a double: the shortest of 15 to 17 significant digits
# that reads back to the same number
.spsq_repr <- function(x) {
  if (!is.finite(x)) return(if (is.na(x)) "nan" else if (x > 0) "inf" else "-inf")
  for (dg in 15:17) {
    s <- formatC(x, digits = dg, format = "g")
    if (as.numeric(s) == x) break
  }
  s <- trimws(s)
  if (!grepl("[.e]", s)) s <- paste0(s, ".0")
  s
}

# constraints in the old style of SciPy (standardize_constraints with
# new_constraint_to_old): lists of list(fun, jac) per type, the Jacobian NULL
# when it is to be approximated
.spsq_old_constraints <- function(constraints, x0) {
  if (is.null(constraints)) constraints <- list()
  if (!is.null(constraints$type) || !is.null(constraints$fun)) constraints <- list(constraints)
  eq <- list()
  ineq <- list()
  for (con in constraints) {
    type <- tolower(con$type %||% "nonlinear")
    if (type %in% c("eq", "ineq")) {
      if (is.null(con$fun)) stop("Constraint has no function defined.", call. = FALSE)
      item <- local({
        f <- con$fun
        list(fun = function(x) as.numeric(f(x)), jac = if (is.function(con$jac)) con$jac)
      })
      if (type == "eq") eq <- c(eq, list(item)) else ineq <- c(ineq, list(item))
      next
    }
    if (type != "nonlinear") stop("Unknown constraint type '", con$type, "'.", call. = FALSE)
    # the text of SciPy, missing space included
    if (isTRUE(con$keep_feasible) || !is.null(con$finite_diff_rel_step))
      warning("Constraint options `finite_diff_jac_sparsity`, `finite_diff_rel_step`, `keep_feasible`, and `hess`",
              "are ignored by this method.", call. = FALSE)
    local({
      fun <- con$fun
      jac <- if (is.function(con$jac)) con$jac
      mc <- length(as.numeric(fun(x0)))
      lb <- rep_len(as.numeric(con$lb), mc)
      ub <- rep_len(as.numeric(con$ub), mc)
      i_eq <- lb == ub
      i_bb <- xor(lb != -Inf, i_eq)
      i_ba <- xor(ub != Inf, i_eq)
      if (any(lb == -Inf & ub == Inf))
        warning("At least one constraint is unbounded above and below. Such constraints are ignored.", call. = FALSE)
      jac_rows <- function(x) {
        J <- jac(x)
        if (!is.matrix(J)) J <- matrix(J, nrow = mc)
        J
      }
      if (any(i_eq)) {
        eq <<- c(eq, list(list(fun = function(x) as.numeric(fun(x))[i_eq] - lb[i_eq],
                               jac = if (!is.null(jac)) function(x) jac_rows(x)[i_eq, , drop = FALSE])))
      }
      if (any(i_bb) || any(i_ba)) {
        ineq <<- c(ineq, list(list(
          fun = function(x) {
            y <- as.numeric(fun(x))
            c(y[i_bb] - lb[i_bb], -(y[i_ba] - ub[i_ba]))
          },
          jac = if (!is.null(jac)) function(x) {
            J <- jac_rows(x)
            rbind(J[i_bb, , drop = FALSE], -J[i_ba, , drop = FALSE])
          })))
      }
    })
  }
  list(eq = eq, ineq = ineq)
}

# result of minimize() when the bounds fix every variable
# (_optimize_result_for_equal_bounds)
.sp_equal_bounds_result <- function(fs, lb, constraints) {
  x <- lb
  success <- TRUE
  message <- "All independent variables were fixed by bounds."
  if (length(constraints)) {
    message <- "All independent variables were fixed by bounds at values that satisfy the constraints."
    if (!is.null(constraints$type) || !is.null(constraints$fun)) constraints <- list(constraints)
  }
  maxcv <- 0
  for (con in constraints) {
    type <- tolower(con$type %||% "nonlinear")
    y <- as.numeric(con$fun(x))
    lo <- if (type == "nonlinear") rep_len(con$lb, length(y)) else if (type == "eq") 0 else 0
    up <- if (type == "nonlinear") rep_len(con$ub, length(y)) else if (type == "eq") 0 else Inf
    violation <- pmax(lo - y, 0) + pmax(y - up, 0)
    if (sum(violation) != 0) {
      maxcv <- max(maxcv, max(violation))
      success <- FALSE
      message <- paste0("All independent variables were fixed by bounds, but the independent variables do not ",
                        "satisfy the constraints exactly. (Maximum violation: ", .spsq_repr(maxcv), ").")
    }
  }
  res <- list(x = x, fun = fs(x), nit = 0L, nfev = 1L, njev = 0L, nhev = 0L, success = success, message = message)
  if (success) res$status <- 0L
  res
}

# SciPy's minimize(method = "SLSQP"): the variables fixed by equal bounds are
# removed when derivatives are approximated (_minimize), then
# _minimize_slsqp runs on the rest
.sp_slsqp <- function(fun, x0, jac = NULL, bounds = NULL, constraints = list(), maxiter = 100L, ftol = 1e-6,
                      iprint = 1L, disp = FALSE, eps = 1.4901161193847656e-08, callback = NULL,
                      finite_diff_rel_step = NULL, workers = NULL, ...) {
  extra <- names(list(...))
  if (length(extra)) warning("Unknown solver options: ", paste(extra, collapse = ", "), call. = FALSE)
  x0 <- as.numeric(x0)
  cons <- .spsq_old_constraints(constraints, x0)
  if (!is.null(bounds)) {
    lb <- rep_len(as.numeric(bounds$lb), length(x0))
    ub <- rep_len(as.numeric(bounds$ub), length(x0))
    i_fixed <- lb == ub
    fd_needed <- !is.function(jac) || any(!vapply(c(cons$eq, cons$ineq), function(con) is.function(con$jac), TRUE))
    if (any(i_fixed) && fd_needed) {
      x_fixed <- lb[i_fixed]
      full <- function(x) {
        out <- numeric(length(i_fixed))
        out[i_fixed] <- x_fixed
        out[!i_fixed] <- x
        out
      }
      fun0 <- fun
      jac0 <- jac
      fun_r <- function(x) fun0(full(x))
      jac_r <- if (is.function(jac0)) function(x) as.numeric(jac0(full(x)))[!i_fixed] else jac0
      cb_r <- if (is.function(callback)) function(r) callback(list(x = full(r$x), fun = r$fun))
      reduce <- function(con) {
        f0 <- con$fun
        j0 <- con$jac
        list(fun = function(x) f0(full(x)),
             jac = if (is.function(j0)) function(x) {
               J <- j0(full(x))
               if (!is.matrix(J)) J <- matrix(J, nrow = 1L)
               J[, !i_fixed, drop = FALSE]
             })
      }
      cons_r <- list(eq = lapply(cons$eq, reduce), ineq = lapply(cons$ineq, reduce))
      res <- .spsq_minimize(fun_r, x0[!i_fixed], jac_r, list(lb = lb[!i_fixed], ub = ub[!i_fixed]), cons_r,
                            maxiter, ftol, iprint, disp, eps, cb_r, finite_diff_rel_step)
      res$x <- full(res$x)
      g <- rep(NaN, length(i_fixed))
      g[!i_fixed] <- res$jac
      res$jac <- g
      return(res)
    }
  }
  .spsq_minimize(fun, x0, jac, bounds, cons, maxiter, ftol, iprint, disp, eps, callback, finite_diff_rel_step)
}

# _minimize_slsqp
.spsq_minimize <- function(fun, x0, jac, bounds, cons, maxiter, ftol, iprint, disp, eps, callback,
                           finite_diff_rel_step) {
  if (!isTRUE(disp)) iprint <- 0L
  x <- as.numeric(x0)
  n <- length(x)
  lbv <- rep(-Inf, n)
  ubv <- rep(Inf, n)
  if (!is.null(bounds)) {
    lbv <- rep_len(as.numeric(bounds$lb), n)
    ubv <- rep_len(as.numeric(bounds$ub), n)
  }
  x <- pmin(pmax(x, lbv), ubv)
  clip <- function(x) {
    if (any(x < lbv | x > ubv)) {
      warning("Values in x were outside bounds during a minimize step, clipping to bounds", call. = FALSE)
      x <- pmin(pmax(x, lbv), ubv)
    }
    x
  }
  cjac_factory <- function(cfun) {
    function(x) {
      x <- clip(x)
      r <- if (is.character(jac) && jac %in% c("2-point", "3-point", "cs")) {
        .sp_approx_derivative(cfun, x, method = if (jac == "cs") "2-point" else jac,
                              rel_step = finite_diff_rel_step, lb = lbv, ub = ubv)
      } else {
        .sp_approx_derivative(cfun, x, method = "2-point", abs_step = eps, lb = lbv, ub = ubv)
      }
      r$J
    }
  }
  fill <- function(con) {
    if (!is.function(con$jac)) con$jac <- cjac_factory(con$fun)
    con
  }
  cons$eq <- lapply(cons$eq, fill)
  cons$ineq <- lapply(cons$ineq, fill)
  meq <- sum(vapply(cons$eq, function(con) length(as.numeric(con$fun(x))), 0L))
  mieq <- sum(vapply(cons$ineq, function(con) length(as.numeric(con$fun(x))), 0L))
  m <- meq + mieq
  eval_constraint <- function(x) {
    as.numeric(c(unlist(lapply(cons$eq, function(con) as.numeric(con$fun(x)))),
                 unlist(lapply(cons$ineq, function(con) as.numeric(con$fun(x))))))
  }
  eval_normals <- function(x) {
    C <- matrix(0, m, n)
    row <- 0L
    for (con in c(cons$eq, cons$ineq)) {
      J <- con$jac(x)
      if (!is.matrix(J)) J <- matrix(J, ncol = n)
      C[row + seq_len(nrow(J)), ] <- J
      row <- row + nrow(J)
    }
    C
  }
  if (!is.null(bounds) && any(lbv > ubv))
    stop("SLSQP Error: lb > ub in bounds ", paste(lbv > ubv, collapse = ", "), ".", call. = FALSE)
  xl <- ifelse(is.finite(lbv), lbv, NaN)
  xu <- ifelse(is.finite(ubv), ubv, NaN)
  sf <- .sp_prepare_scalar_function(fun, x, jac, lb = lbv, ub = ubv, epsilon = eps,
                                    finite_diff_rel_step = finite_diff_rel_step)
  if (iprint >= 2L) cat(sprintf("%5s %5s %16s %16s\n", "NIT", "FC", "OBJFUN", "GNORM"))
  fx <- sf$fun(clip(x))
  g <- sf$grad(clip(x))
  C <- eval_normals(x)
  d <- eval_constraint(x)
  S <- new.env(parent = emptyenv())
  S$acc <- ftol
  S$alpha <- 0; S$f0 <- 0; S$gs <- 0; S$h1 <- 0; S$h2 <- 0; S$h3 <- 0; S$h4 <- 0; S$t <- 0; S$t0 <- 0
  S$tol <- 10 * ftol
  S$inconsistent <- 0L; S$reset <- 0L; S$iter <- 0L; S$itermax <- as.integer(maxiter); S$line <- 0L
  S$m <- m; S$meq <- meq; S$mode <- 0L; S$n <- n
  S$bfgs <- numeric((n * (n + 1L)) %/% 2L)
  S$x0 <- numeric(n); S$mu <- numeric(m); S$s <- numeric(n + 1L); S$u <- numeric(n + 1L); S$v <- numeric(n + 1L)
  S$jw <- integer(max(m + 2L * n + 2L, 1L))
  S$mult <- numeric(max(1L, m + 2L * n + 2L))
  S$sol <- x
  iter_prev <- 0L
  repeat {
    .spsq_body(S, fx, g, C, d, xl, xu)
    x <- S$sol
    if (S$mode == 1L) {
      fx <- sf$fun(x)
      d <- eval_constraint(x)
    }
    if (S$mode == -1L) {
      g <- sf$grad(x)
      C <- eval_normals(x)
    }
    if (S$iter > iter_prev) {
      if (is.function(callback) && isTRUE(callback(list(x = x, fun = fx)))) break
      if (iprint >= 2L) cat(sprintf("%5d %5d %16.6E %16.6E\n", S$iter, sf$nfev(), fx, sqrt(sum(g * g))))
    }
    if (abs(S$mode) != 1L) break
    iter_prev <- S$iter
  }
  msg <- .spsq_exit_modes[[as.character(S$mode)]]
  if (iprint >= 1L) {
    cat(msg, "    (Exit mode ", S$mode, ")\n", sep = "")
    cat("            Current function value:", .spsq_repr(fx), "\n")
    cat("            Iterations:", S$iter, "\n")
    cat("            Function evaluations:", sf$nfev(), "\n")
    cat("            Gradient evaluations:", sf$ngev, "\n")
  }
  list(x = x, fun = fx, jac = g, nit = S$iter, nfev = sf$nfev(), njev = sf$ngev, status = S$mode, message = msg,
       success = S$mode == 0L, multipliers = S$mult[seq_len(m)])
}
