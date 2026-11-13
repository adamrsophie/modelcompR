# The wrapper algorithms of NonlinearSolve (src/extension_algs.jl and the
# package extensions in ext/), each implemented natively after the library it
# wraps: MINPACK (hybrd, hybrj, lmdif, lmder), NLsolve 4.5 (trust_region,
# newton, anderson, broyden), LeastSquaresOptim (LevenbergMarquardt, Dogleg),
# FastLevenbergMarquardt (lmsolve!), NLSolvers (LineSearch with Newton),
# SpeedMapping (ACX for maps), FixedPointAcceleration (fixed_point),
# SIAMFANLEquations (nsol, ptcsol, aasol, secant) and PETSc (SNES newtonls
# with the bt line search). Statistics follow what each extension reports.

#' Wrapped nonlinear solvers
#'
#' The algorithms NonlinearSolve takes from other Julia packages, implemented
#' natively after those packages:
#'
#' * `CMINPACK(method)`: MINPACK. `"hybr"` is Powell's hybrid method
#'   (`hybrj` with a Jacobian, `hybrd1` with forward differences), `"lm"` the
#'   Levenberg Marquardt method of More (`lmder`, `lmdif1`); `"auto"` picks
#'   `"lm"` for least squares problems. As in MINPACK.jl the call counts are
#'   recorded only with `store_trace = TRUE`.
#' * `NLsolveJL(method)`: NLsolve with `"trust_region"` (dogleg with
#'   autoscaling), `"newton"` (with `linesearch = "Static"` or
#'   `"BackTracking"`), `"anderson"` and `"broyden"`.
#' * `LeastSquaresOptimJL(alg)`: the `"lm"` and `"dogleg"` methods of
#'   LeastSquaresOptim with a pivoted QR (`"qr"`) or Cholesky solver.
#' * `FastLevenbergMarquardtJL(linsolve)`: FastLevenbergMarquardt. As in the
#'   extension, the residual of the solution is the sum of squares.
#'   `factorupdate = "nielson"` (the spelling NonlinearSolve accepts) and
#'   `"nielsen"` (the spelling of FastLevenbergMarquardt) both select
#'   Nielsen's damping update.
#' * `NLSolversJL(method)`: NLSolvers; `method = "LineSearch"` is Newton's
#'   method with the static unit step (`LineSearch(Newton(), Static(1))`).
#' * `SpeedMappingJL()`: the alternating cyclic extrapolation of SpeedMapping
#'   applied to the map `u + f(u)`.
#' * `FixedPointAccelerationJL(algorithm)`: FixedPointAcceleration on the map
#'   `u + f(u)` with `"Anderson"`, `"MPE"`, `"RRE"`, `"VEA"`, `"SEA"`,
#'   `"Simple"`, `"Aitken"` or `"Newton"`.
#' * `SIAMFANLEquationsJL(method)`: SIAMFANLEquations with `"newton"`,
#'   `"pseudotransient"`, `"secant"` (scalar problems) or `"anderson"`
#'   (`m = 0` is Picard iteration on `u + f(u)`).
#' * `PETScSNES()`: the default PETSc SNES, Newton with the cubic
#'   backtracking line search and the `snes_atol`, `snes_rtol` and
#'   `snes_stol` tests.
#'
#' @param method,alg,algorithm,linsolve The method of the wrapped package.
#' @param autodiff The Jacobian choice (an ADType, a string, or `NULL`;
#'   finite differences are used when the problem has no Jacobian).
#' @param linesearch The NLsolve line search, `"Static"` or `"BackTracking"`.
#' @param beta Relaxation (NLsolve, SIAMFANLEquations Anderson).
#' @param factor,autoscale,m Settings of NLsolve.
#' @param factoraccept,factorreject,factorupdate,minscale,maxscale,minfactor,maxfactor
#'   Settings of FastLevenbergMarquardt.
#' @param sigma_min,stabilize,check_obj,orders Settings of SpeedMapping.
#' @param extrapolation_period,replace_invalids,dampening,condition_number_threshold
#'   Settings of FixedPointAcceleration.
#' @param delta Initial pseudo time step of `"pseudotransient"`.
#' @param petsclib,mpi_comm,... Kept for compatibility (`...` holds SNES
#'   options such as `snes_stol`, `snes_max_funcs`, `snes_linesearch_order`).
#' @return An algorithm object for [solve.NonlinearProblem()].
#' @examples
#' f <- function(u, p) c(u[1]^2 + u[2]^2 - 4, u[1] * u[2] - 1)
#' prob <- NonlinearProblem(f, c(2, 0.5))
#' solve(prob, CMINPACK())$u
#' solve(prob, NLsolveJL())$u
#' solve(prob, SIAMFANLEquationsJL())$u
#' solve(prob, PETScSNES())$u
#' @name nonlinear_wrappers
NULL

.wr_alg <- function(name, run, ...) {
  a <- .nl_alg(name, ...)
  a$run <- function(F, J, x0, opts, lb, ub, ls) run(a, F, J, x0, opts, ls)
  a
}

.wr_out <- function(u, resid, retcode, iters, nf, njacs, nfactors = -1L, nsolve = -1L, nsteps = iters, original = NULL)
  list(u = u, resid = resid, retcode = retcode, iters = nsteps,
       stats = list(nf = nf, njacs = njacs, nfactors = nfactors, nsolve = nsolve), original = original)

.wr_tol <- function(x) x %||% .nl_default_tol()
.wr_jac <- function(F, J) if (is.null(J)) function(x) .jac_central(F, x) else function(x) as.matrix(J(x))
.wr_l2 <- function(v) sqrt(sum(v * v))

# a * b + c with a single rounding (Julia's fma, and muladd on hardware with
# FMA): the exact product by Veltkamp splitting, then an exact sum
.wr_fma <- function(a, b, c) {
  p <- a * b
  t <- 134217729 * a; ah <- t - (t - a); al <- a - ah
  t <- 134217729 * b; bh <- t - (t - b); bl <- b - bh
  e <- ((ah * bh - p) + ah * bl + al * bh) + al * bl
  s <- p + c
  bb <- s - p
  out <- s + (((p - (s - bb)) + (c - bb)) + e)
  bad <- !is.finite(out) | !is.finite(p)
  out[bad] <- (p + c)[bad]
  out
}

# Small BLAS kernels as OpenBLAS evaluates them for short vectors: dot
# products accumulated left to right with fma, axpy with fma, and the column
# oriented upper triangular solve (trsv) with fma
.wr_dot <- function(a, b) {
  s <- 0
  for (i in seq_along(a)) s <- .wr_fma(a[i], b[i], s)
  s
}
# UpperTriangular(view(R, 1:k, 1:k)): entries below the diagonal are ignored
.wr_upper <- function(R, k) {
  Rk <- R[seq_len(k), seq_len(k), drop = FALSE]
  Rk[lower.tri(Rk)] <- 0
  Rk
}
.wr_tmv <- function(A, x) vapply(seq_len(ncol(A)), function(j) .wr_dot(A[, j], x), 0)
.wr_utrisolve <- function(R, b) {
  for (j in rev(seq_along(b))) {
    b[j] <- b[j] / R[j, j]
    if (j > 1L) for (i in seq_len(j - 1L)) b[i] <- .wr_fma(-b[j], R[i, j], b[i])
  }
  b
}

# Julia's sum!(abs2, D', J): column sums of squares accumulated with fma
.wr_colsumsq <- function(J) {
  out <- numeric(ncol(J))
  for (i in seq_len(nrow(J))) out <- .wr_fma(J[i, ], J[i, ], out)
  out
}

# Julia's hypot (the correctly rounded branch for hardware with fma)
.wr_hypot <- function(x, y) {
  ax <- abs(x); ay <- abs(y)
  if (is.infinite(ax) || is.infinite(ay)) return(Inf)
  if (ay > ax) { tt <- ax; ax <- ay; ay <- tt }
  if (ay <= ax * sqrt(.Machine$double.eps / 2)) return(ax)
  scale <- 2^-563
  if (ax > sqrt(.Machine$double.xmax / 2)) {
    ax <- ax * scale; ay <- ay * scale; scale <- 1 / scale
  } else if (ay < sqrt(.Machine$double.xmin)) {
    ax <- ax / scale; ay <- ay / scale
  } else scale <- 1
  h <- sqrt(.wr_fma(ax, ax, ay * ay))
  hs <- h * h; axs <- ax * ax
  h <- h - (.wr_fma(-ay, ay, hs - axs) + .wr_fma(h, h, -hs) - .wr_fma(ax, ax, -axs)) / (2 * h)
  h * scale
}

# Julia's A \ b: LU for a square matrix (an error when exactly singular), the
# pivoted QR least squares solution otherwise.
.wr_ldiv <- function(A, b) {
  if (nrow(A) != ncol(A)) return(.wr_qr_ldiv(A, b))
  if (nrow(A) <= 12L && all(is.finite(A))) return(.wr_lu_solve(A, as.numeric(b)))
  as.numeric(solve(A, b, tol = 0))
}

# Q'b for a LAPACK (geqp3) QR, one reflector at a time as dorm2r and dlarf
# apply them: w = v'b by an fma dot product, then b - (tau w) v with fma
.wr_qty <- function(q, b) {
  QR <- q$qr; tau <- q$qraux; m <- nrow(QR)
  b <- as.numeric(b)
  for (j in seq_len(min(m, ncol(QR)))) {
    if (tau[j] == 0) next
    idx <- j:m
    v <- c(1, QR[idx[-1L], j])
    w <- .wr_dot(v, b[idx])
    b[idx] <- .wr_fma(-tau[j] * w, v, b[idx])
  }
  b
}

# ldiv!(qr(A, ColumnNorm()), b): the pivoted QR least squares solution. A
# numerically rank deficient A (|R[k, k]| below min(m, n) eps |R[1, 1]|) gets
# the minimum norm solution of the leading rows, which is what LAPACK's
# complete orthogonal factorization gives
.wr_qr_ldiv <- function(A, b) {
  A <- as.matrix(A); n <- ncol(A)
  q <- qr(A, LAPACK = TRUE)
  y <- .wr_qty(q, b)
  k <- min(nrow(A), n)
  d <- abs(diag(q$qr)[seq_len(k)])
  rnk <- if (!k || d[1L] == 0) 0L else sum(d > k * .Machine$double.eps * d[1L])
  z <- numeric(n)
  if (rnk == n) {
    z <- .wr_utrisolve(.wr_upper(q$qr, n), y[seq_len(n)])
  } else if (rnk > 0L) {
    Rr <- q$qr[seq_len(rnk), , drop = FALSE]
    Rr[lower.tri(Rr)] <- 0
    z <- as.numeric(.wr_pinv(Rr) %*% y[seq_len(rnk)])
  }
  out <- numeric(n); out[q$pivot] <- z
  out
}

# The LU solve of a small matrix in the operation order of OpenBLAS (getf2
# and getrs): partial pivoting on the first largest entry, multipliers scaled
# by the reciprocal of the pivot, fma updates and fma triangular solves
.wr_lu_solve <- function(A, b) {
  n <- nrow(A)
  for (j in seq_len(n)) {
    p <- j - 1L + which.max(abs(A[j:n, j]))
    if (p != j) { A[c(j, p), ] <- A[c(p, j), ]; b[c(j, p)] <- b[c(p, j)] }
    if (A[j, j] == 0) stop("Lapack routine dgesv: system is exactly singular: U[", j, ",", j, "] = 0", call. = FALSE)
    if (j < n) {
      rows <- (j + 1L):n
      A[rows, j] <- A[rows, j] * (1 / A[j, j])
      for (cc in rows) A[rows, cc] <- .wr_fma(-A[rows, j], A[j, cc], A[rows, cc])
    }
  }
  for (j in seq_len(n - 1L)) {
    rows <- (j + 1L):n
    b[rows] <- .wr_fma(-b[j], A[rows, j], b[rows])
  }
  .wr_utrisolve(A, b)
}

.wr_pinv <- function(A) {
  A <- as.matrix(A)
  s <- svd(A)
  if (!length(s$d) || max(s$d) == 0) return(t(A) * 0)
  tol <- .Machine$double.eps * min(dim(A)) * max(s$d)
  keep <- s$d > tol
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

# ---------------------------------------------------------------------------
# MINPACK

#' @rdname nonlinear_wrappers
#' @export
CMINPACK <- function(method = "auto", autodiff = NULL) {
  method <- match.arg(method, c("auto", "hybr", "lm", "lmdif", "hybrd"))
  .wr_alg("CMINPACK", .mp_run, method = method, autodiff = autodiff)
}

.mp_enorm <- function(v) sqrt(sum(v * v))

# qrfac: Householder QR with optional column pivoting (MINPACK)
.mp_qrfac <- function(a, pivot) {
  m <- nrow(a); n <- ncol(a)
  epsmch <- .Machine$double.eps
  acnorm <- apply(a, 2, .mp_enorm)
  rdiag <- acnorm; wa <- rdiag; ipvt <- seq_len(n)
  for (j in seq_len(min(m, n))) {
    if (pivot) {
      kmax <- j
      for (k in j:n) if (rdiag[k] > rdiag[kmax]) kmax <- k
      if (kmax != j) {
        tmp <- a[, j]; a[, j] <- a[, kmax]; a[, kmax] <- tmp
        rdiag[kmax] <- rdiag[j]; wa[kmax] <- wa[j]
        k <- ipvt[j]; ipvt[j] <- ipvt[kmax]; ipvt[kmax] <- k
      }
    }
    rows <- j:m
    ajnorm <- .mp_enorm(a[rows, j])
    if (ajnorm != 0) {
      if (a[j, j] < 0) ajnorm <- -ajnorm
      a[rows, j] <- a[rows, j] / ajnorm
      a[j, j] <- a[j, j] + 1
      if (j < n) for (k in (j + 1):n) {
        s <- sum(a[rows, j] * a[rows, k])
        temp <- s / a[j, j]
        a[rows, k] <- a[rows, k] - temp * a[rows, j]
        if (pivot && rdiag[k] != 0) {
          temp <- a[j, k] / rdiag[k]
          rdiag[k] <- rdiag[k] * sqrt(max(0, 1 - temp * temp))
          if (0.05 * (rdiag[k] / wa[k])^2 <= epsmch) {
            rdiag[k] <- if (j < m) .mp_enorm(a[(j + 1):m, k]) else 0
            wa[k] <- rdiag[k]
          }
        }
      }
    }
    rdiag[j] <- -ajnorm
  }
  list(a = a, rdiag = rdiag, acnorm = acnorm, ipvt = ipvt)
}

# qform: the orthogonal matrix of a qrfac factorization (square case)
.mp_qform <- function(a) {
  n <- ncol(a)
  if (n >= 2) for (j in 2:n) a[seq_len(j - 1L), j] <- 0
  for (k in n:1) {
    rows <- k:n
    wa <- a[rows, k]
    a[rows, k] <- 0
    a[k, k] <- 1
    if (wa[1L] != 0) for (j in k:n) {
      s <- sum(a[rows, j] * wa)
      temp <- s / wa[1L]
      a[rows, j] <- a[rows, j] - temp * wa
    }
  }
  a
}

# dogleg step of hybrd and hybrj (R upper triangular, full storage)
.mp_dogleg <- function(R, diag, qtb, delta) {
  n <- length(qtb); epsmch <- .Machine$double.eps
  x <- numeric(n)
  for (j in n:1) {
    s <- 0
    if (j < n) for (i in (j + 1):n) s <- s + R[j, i] * x[i]
    temp <- R[j, j]
    if (temp == 0) {
      temp <- max(abs(R[seq_len(j), j]))
      temp <- epsmch * temp
      if (temp == 0) temp <- epsmch
    }
    x[j] <- (qtb[j] - s) / temp
  }
  wa2 <- diag * x
  qnorm <- .mp_enorm(wa2)
  if (qnorm <= delta) return(x)
  wa1 <- numeric(n)
  for (j in seq_len(n)) {
    temp <- qtb[j]
    for (i in j:n) wa1[i] <- wa1[i] + R[j, i] * temp
    wa1[j] <- wa1[j] / diag[j]
  }
  gnorm <- .mp_enorm(wa1)
  sgnorm <- 0
  alpha <- delta / qnorm
  if (gnorm != 0) {
    wa1 <- (wa1 / gnorm) / diag
    for (j in seq_len(n)) {
      s <- 0
      for (i in j:n) s <- s + R[j, i] * wa1[i]
      wa2[j] <- s
    }
    temp <- .mp_enorm(wa2)
    sgnorm <- (gnorm / temp) / temp
    alpha <- 0
    if (sgnorm < delta) {
      bnorm <- .mp_enorm(qtb)
      temp <- (bnorm / gnorm) * (bnorm / qnorm) * (sgnorm / delta)
      temp <- temp - (delta / qnorm) * (sgnorm / delta)^2 +
        sqrt((temp - delta / qnorm)^2 + (1 - (delta / qnorm)^2) * (1 - (sgnorm / delta)^2))
      alpha <- ((delta / qnorm) * (1 - (sgnorm / delta)^2)) / temp
    }
  }
  temp <- (1 - alpha) * min(sgnorm, delta)
  temp * wa1 + alpha * x
}

# r1updt on the upper triangular R (the transpose of MINPACK's lower
# trapezoidal storage); returns R, v and w (the rotations) and sing
.mp_r1updt <- function(R, u, v) {
  n <- length(u); giant <- .Machine$double.xmax
  w <- numeric(n)
  w[n] <- R[n, n]
  if (n >= 2) for (j in (n - 1):1) {
    w[j] <- 0
    if (v[j] != 0) {
      if (abs(v[n]) < abs(v[j])) {
        cotan <- v[n] / v[j]; sn <- 0.5 / sqrt(0.25 + 0.25 * cotan^2); cs <- sn * cotan
        tau <- 1
        if (abs(cs) * giant > 1) tau <- 1 / cs
      } else {
        tn <- v[j] / v[n]; cs <- 0.5 / sqrt(0.25 + 0.25 * tn^2); sn <- cs * tn
        tau <- sn
      }
      v[n] <- sn * v[j] + cs * v[n]
      v[j] <- tau
      for (i in j:n) {
        temp <- cs * R[j, i] - sn * w[i]
        w[i] <- sn * R[j, i] + cs * w[i]
        R[j, i] <- temp
      }
    }
  }
  w <- w + v[n] * u
  sing <- FALSE
  if (n >= 2) for (j in 1:(n - 1)) {
    if (w[j] != 0) {
      if (abs(R[j, j]) < abs(w[j])) {
        cotan <- R[j, j] / w[j]; sn <- 0.5 / sqrt(0.25 + 0.25 * cotan^2); cs <- sn * cotan
        tau <- 1
        if (abs(cs) * giant > 1) tau <- 1 / cs
      } else {
        tn <- w[j] / R[j, j]; cs <- 0.5 / sqrt(0.25 + 0.25 * tn^2); sn <- cs * tn
        tau <- sn
      }
      for (i in j:n) {
        temp <- cs * R[j, i] + sn * w[i]
        w[i] <- -sn * R[j, i] + cs * w[i]
        R[j, i] <- temp
      }
      w[j] <- tau
    }
    if (R[j, j] == 0) sing <- TRUE
  }
  R[n, n] <- w[n]
  if (R[n, n] == 0) sing <- TRUE
  list(R = R, v = v, w = w, sing = sing)
}

# r1mpyq: A times the rotations stored in v and w
.mp_r1mpyq <- function(A, v, w) {
  A <- as.matrix(A)
  n <- ncol(A)
  if (n < 2) return(A)
  cs_sn <- function(z) if (abs(z) > 1) { cs <- 1 / z; c(cs, sqrt(1 - cs * cs)) } else c(sqrt(1 - z * z), z)
  for (j in (n - 1):1) {
    r <- cs_sn(v[j]); cs <- r[1]; sn <- r[2]
    temp <- cs * A[, j] - sn * A[, n]
    A[, n] <- sn * A[, j] + cs * A[, n]
    A[, j] <- temp
  }
  for (j in 1:(n - 1)) {
    r <- cs_sn(w[j]); cs <- r[1]; sn <- r[2]
    temp <- cs * A[, j] + sn * A[, n]
    A[, n] <- -sn * A[, j] + cs * A[, n]
    A[, j] <- temp
  }
  A
}

# forward difference Jacobian of fdjac1 and fdjac2 (epsfcn = 0)
.mp_fdjac <- function(fcn, x, fvec) {
  eps <- sqrt(.Machine$double.eps)
  n <- length(x)
  J <- matrix(0, length(fvec), n)
  for (j in seq_len(n)) {
    temp <- x[j]
    h <- eps * abs(temp)
    if (h == 0) h <- eps
    x[j] <- temp + h
    J[, j] <- (fcn(x, 2L) - fvec) / h
    x[j] <- temp
  }
  J
}

.mp_hybr <- function(fcn, jac, x, xtol, maxfev) {
  n <- length(x); epsmch <- .Machine$double.eps
  diag <- rep(1, n); factor <- 100
  fvec <- fcn(x, 1L); nfev <- 1L
  fnorm <- .mp_enorm(fvec)
  iter <- 1L; ncsuc <- 0L; ncfail <- 0L; nslow1 <- 0L; nslow2 <- 0L
  info <- 0L
  repeat {
    jeval <- TRUE
    if (is.null(jac)) {
      fjac <- .mp_fdjac(fcn, x, fvec)
      nfev <- nfev + n
    } else fjac <- jac(x)
    qf <- .mp_qrfac(fjac, FALSE)
    fjac <- qf$a
    if (iter == 1L) {
      xnorm <- .mp_enorm(diag * x)
      delta <- factor * xnorm
      if (delta == 0) delta <- factor
    }
    qtf <- fvec
    for (j in seq_len(n)) if (fjac[j, j] != 0) {
      rows <- j:n
      s <- sum(fjac[rows, j] * qtf[rows])
      temp <- -s / fjac[j, j]
      qtf[rows] <- qtf[rows] + fjac[rows, j] * temp
    }
    R <- matrix(0, n, n)
    for (j in seq_len(n)) {
      if (j > 1) R[seq_len(j - 1L), j] <- fjac[seq_len(j - 1L), j]
      R[j, j] <- qf$rdiag[j]
    }
    Q <- .mp_qform(fjac)
    repeat {
      p <- -.mp_dogleg(R, diag, qtf, delta)
      wa2 <- x + p
      pnorm <- .mp_enorm(diag * p)
      if (iter == 1L) delta <- min(delta, pnorm)
      wa4 <- fcn(wa2, 1L); nfev <- nfev + 1L
      if (is.null(wa4)) return(list(x = x, fvec = fvec, info = -1L))
      fnorm1 <- .mp_enorm(wa4)
      actred <- -1
      if (fnorm1 < fnorm) actred <- 1 - (fnorm1 / fnorm)^2
      wa3 <- numeric(n)
      for (i in seq_len(n)) {
        s <- 0
        for (j in i:n) s <- s + R[i, j] * p[j]
        wa3[i] <- qtf[i] + s
      }
      temp <- .mp_enorm(wa3)
      prered <- 0
      if (temp < fnorm) prered <- 1 - (temp / fnorm)^2
      ratio <- if (prered > 0) actred / prered else 0
      if (ratio < 0.1) {
        ncsuc <- 0L; ncfail <- ncfail + 1L
        delta <- 0.5 * delta
      } else {
        ncfail <- 0L; ncsuc <- ncsuc + 1L
        if (ratio >= 0.5 || ncsuc > 1L) delta <- max(delta, pnorm / 0.5)
        if (abs(ratio - 1) <= 0.1) delta <- pnorm / 0.5
      }
      if (ratio >= 1e-4) {
        x <- wa2
        xnorm <- .mp_enorm(diag * x)
        fvec <- wa4
        fnorm <- fnorm1
        iter <- iter + 1L
      }
      nslow1 <- nslow1 + 1L
      if (actred >= 0.001) nslow1 <- 0L
      if (jeval) nslow2 <- nslow2 + 1L
      if (actred >= 0.1) nslow2 <- 0L
      if (delta <= xtol * xnorm || fnorm == 0) info <- 1L
      if (info != 0L) return(list(x = x, fvec = fvec, info = info))
      if (nfev >= maxfev) info <- 5L
      if (0.1 * max(0.1 * delta, pnorm) <= epsmch * xnorm) info <- 3L
      if (nslow2 == 5L) info <- 4L
      if (nslow1 == 10L) info <- 5L
      if (info != 0L) return(list(x = x, fvec = fvec, info = info))
      if (ncfail == 2L) break
      qtw <- as.numeric(crossprod(Q, wa4))
      wa2v <- (qtw - wa3) / pnorm
      wa1 <- diag * ((diag * p) / pnorm)
      if (ratio >= 1e-4) qtf <- qtw
      up <- .mp_r1updt(R, wa1, wa2v)
      R <- up$R
      Q <- .mp_r1mpyq(Q, up$v, up$w)
      qtf <- as.numeric(.mp_r1mpyq(matrix(qtf, 1L), up$v, up$w))
      jeval <- FALSE
    }
  }
}

# qrsolv and lmpar of lmder (R holds the qrfac output)
.mp_qrsolv <- function(r, ipvt, diag, qtb) {
  n <- ncol(r)
  x <- numeric(n); wa <- qtb[seq_len(n)]
  for (j in seq_len(n)) {
    if (j < n) r[(j + 1):n, j] <- r[j, (j + 1):n]
    x[j] <- r[j, j]
  }
  sdiag <- numeric(n)
  for (j in seq_len(n)) {
    l <- ipvt[j]
    if (diag[l] != 0) {
      sdiag[j:n] <- 0
      sdiag[j] <- diag[l]
      qtbpj <- 0
      for (k in j:n) {
        if (sdiag[k] != 0) {
          if (abs(r[k, k]) < abs(sdiag[k])) {
            cotan <- r[k, k] / sdiag[k]; sn <- 0.5 / sqrt(0.25 + 0.25 * cotan^2); cs <- sn * cotan
          } else {
            tn <- sdiag[k] / r[k, k]; cs <- 0.5 / sqrt(0.25 + 0.25 * tn^2); sn <- cs * tn
          }
          r[k, k] <- cs * r[k, k] + sn * sdiag[k]
          temp <- cs * wa[k] + sn * qtbpj
          qtbpj <- -sn * wa[k] + cs * qtbpj
          wa[k] <- temp
          if (k < n) for (i in (k + 1):n) {
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
  if (nsing >= 1) for (j in nsing:1) {
    s <- 0
    if (j < nsing) for (i in (j + 1):nsing) s <- s + r[i, j] * wa[i]
    wa[j] <- (wa[j] - s) / sdiag[j]
  }
  x[ipvt] <- wa
  list(x = x, sdiag = sdiag, r = r)
}

.mp_lmpar <- function(r, ipvt, diag, qtb, delta, par) {
  n <- ncol(r); dwarf <- .Machine$double.xmin
  wa1 <- qtb[seq_len(n)]
  nsing <- n
  for (j in seq_len(n)) {
    if (r[j, j] == 0 && nsing == n) nsing <- j - 1L
    if (nsing < n) wa1[j] <- 0
  }
  if (nsing >= 1) for (j in nsing:1) {
    wa1[j] <- wa1[j] / r[j, j]
    temp <- wa1[j]
    if (j > 1) for (i in 1:(j - 1)) wa1[i] <- wa1[i] - r[i, j] * temp
  }
  x <- numeric(n); x[ipvt] <- wa1
  iter <- 0L
  wa2 <- diag * x
  dxnorm <- .mp_enorm(wa2)
  fp <- dxnorm - delta
  if (fp <= 0.1 * delta) return(list(x = x, par = 0))
  parl <- 0
  if (nsing >= n) {
    wa1 <- diag[ipvt] * (wa2[ipvt] / dxnorm)
    for (j in seq_len(n)) {
      s <- 0
      if (j > 1) for (i in 1:(j - 1)) s <- s + r[i, j] * wa1[i]
      wa1[j] <- (wa1[j] - s) / r[j, j]
    }
    temp <- .mp_enorm(wa1)
    parl <- ((fp / delta) / temp) / temp
  }
  for (j in seq_len(n)) {
    s <- 0
    for (i in 1:j) s <- s + r[i, j] * qtb[i]
    wa1[j] <- s / diag[ipvt[j]]
  }
  gnorm <- .mp_enorm(wa1)
  paru <- gnorm / delta
  if (paru == 0) paru <- dwarf / min(delta, 0.1)
  par <- max(par, parl); par <- min(par, paru)
  if (par == 0) par <- gnorm / dxnorm
  repeat {
    iter <- iter + 1L
    if (par == 0) par <- max(dwarf, 0.001 * paru)
    qs <- .mp_qrsolv(r, ipvt, sqrt(par) * diag, qtb)
    x <- qs$x
    wa2 <- diag * x
    dxnorm <- .mp_enorm(wa2)
    temp <- fp
    fp <- dxnorm - delta
    if (abs(fp) <= 0.1 * delta || (parl == 0 && fp <= temp && temp < 0) || iter == 10L) break
    wa1 <- diag[ipvt] * (wa2[ipvt] / dxnorm)
    rr <- qs$r
    for (j in seq_len(n)) {
      wa1[j] <- wa1[j] / qs$sdiag[j]
      temp <- wa1[j]
      if (j < n) for (i in (j + 1):n) wa1[i] <- wa1[i] - rr[i, j] * temp
    }
    temp <- .mp_enorm(wa1)
    parc <- ((fp / delta) / temp) / temp
    if (fp > 0) parl <- max(parl, par)
    if (fp < 0) paru <- min(paru, par)
    par <- max(parl, par + parc)
  }
  list(x = x, par = par)
}

.mp_lm <- function(fcn, jac, x, ftol, xtol, gtol, maxfev) {
  n <- length(x); epsmch <- .Machine$double.eps; factor <- 100
  fvec <- fcn(x, 1L); nfev <- 1L
  m <- length(fvec)
  fnorm <- .mp_enorm(fvec)
  par <- 0; iter <- 1L; info <- 0L
  diag <- numeric(n)
  repeat {
    if (is.null(jac)) {
      fjac <- .mp_fdjac(fcn, x, fvec)
      nfev <- nfev + n
    } else fjac <- jac(x)
    qf <- .mp_qrfac(fjac, TRUE)
    fjac <- qf$a; ipvt <- qf$ipvt
    wa2 <- qf$acnorm
    if (iter == 1L) {
      diag <- ifelse(wa2 == 0, 1, wa2)
      xnorm <- .mp_enorm(diag * x)
      delta <- factor * xnorm
      if (delta == 0) delta <- factor
    }
    wa4 <- fvec
    qtf <- numeric(n)
    for (j in seq_len(n)) {
      if (fjac[j, j] != 0) {
        rows <- j:m
        s <- sum(fjac[rows, j] * wa4[rows])
        temp <- -s / fjac[j, j]
        wa4[rows] <- wa4[rows] + fjac[rows, j] * temp
      }
      fjac[j, j] <- qf$rdiag[j]
      qtf[j] <- wa4[j]
    }
    gnorm <- 0
    if (fnorm != 0) for (j in seq_len(n)) {
      l <- ipvt[j]
      if (wa2[l] != 0) {
        s <- 0
        for (i in 1:j) s <- s + fjac[i, j] * (qtf[i] / fnorm)
        gnorm <- max(gnorm, abs(s / wa2[l]))
      }
    }
    if (gnorm <= gtol) info <- 4L
    if (info != 0L) return(list(x = x, fvec = fvec, info = info))
    diag <- pmax(diag, wa2)
    r <- fjac[seq_len(n), , drop = FALSE]
    repeat {
      lp <- .mp_lmpar(r, ipvt, diag, qtf, delta, par)
      par <- lp$par
      p <- -lp$x
      wa2 <- x + p
      pnorm <- .mp_enorm(diag * p)
      if (iter == 1L) delta <- min(delta, pnorm)
      wa4 <- fcn(wa2, 1L); nfev <- nfev + 1L
      if (is.null(wa4)) return(list(x = x, fvec = fvec, info = -1L))
      fnorm1 <- .mp_enorm(wa4)
      actred <- -1
      if (0.1 * fnorm1 < fnorm) actred <- 1 - (fnorm1 / fnorm)^2
      wa3 <- numeric(n)
      for (j in seq_len(n)) {
        temp <- p[ipvt[j]]
        for (i in 1:j) wa3[i] <- wa3[i] + r[i, j] * temp
      }
      temp1 <- .mp_enorm(wa3) / fnorm
      temp2 <- (sqrt(par) * pnorm) / fnorm
      prered <- temp1^2 + temp2^2 / 0.5
      dirder <- -(temp1^2 + temp2^2)
      ratio <- if (prered != 0) actred / prered else 0
      if (ratio <= 0.25) {
        temp <- if (actred >= 0) 0.5 else 0.5 * dirder / (dirder + 0.5 * actred)
        if (0.1 * fnorm1 >= fnorm || temp < 0.1) temp <- 0.1
        delta <- temp * min(delta, pnorm / 0.1)
        par <- par / temp
      } else if (par == 0 || ratio >= 0.75) {
        delta <- pnorm / 0.5
        par <- 0.5 * par
      }
      if (ratio >= 1e-4) {
        x <- wa2
        xnorm <- .mp_enorm(diag * x)
        fvec <- wa4
        fnorm <- fnorm1
        iter <- iter + 1L
      }
      if (abs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1) info <- 1L
      if (delta <= xtol * xnorm) info <- 2L
      if (abs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1 && info == 2L) info <- 3L
      if (info != 0L) return(list(x = x, fvec = fvec, info = info))
      if (nfev >= maxfev) info <- 5L
      if (abs(actred) <= epsmch && prered <= epsmch && 0.5 * ratio <= 1) info <- 6L
      if (delta <= epsmch * xnorm) info <- 7L
      if (gnorm <= epsmch) info <- 8L
      if (info != 0L) return(list(x = x, fvec = fvec, info = info))
      if (ratio >= 1e-4) break
    }
  }
}

.mp_run <- function(a, F, J, x0, opts, ls) {
  tol <- .wr_tol(opts$abstol)
  maxit <- as.integer(opts$maxiters %||% 1000L)
  method <- if (a$method == "auto") (if (ls) "lm" else "hybr") else a$method
  if (!is.null(J) && method %in% c("lmdif", "hybrd")) stop("unknown method ", method, call. = FALSE)
  tracing <- isTRUE(opts$store_trace)
  cnt <- new.env(parent = emptyenv()); cnt$f <- 0L; cnt$g <- 0L; cnt$prev <- 0L
  # MINPACK.jl counts calls (and stops after `maxiters` residual calls) only when tracing
  fcn <- function(x, iflag) {
    fx <- as.numeric(F(x))
    if (tracing) {
      if (iflag == 2L) { if (cnt$prev == 1L) cnt$g <- cnt$g + 1L } else cnt$f <- cnt$f + 1L
      cnt$prev <- iflag
      if (cnt$f > maxit) return(NULL)
    }
    fx
  }
  jac <- if (is.null(J)) NULL else function(x) {
    if (tracing) { if (cnt$prev == 1L) cnt$g <- cnt$g + 1L; cnt$prev <- 2L }
    as.matrix(J(x))
  }
  n <- length(x0)
  out <- if (method %in% c("hybr", "hybrd")) {
    .mp_hybr(fcn, jac, x0, tol, if (is.null(J) && method == "hybr") 200L * (n + 1L) else .Machine$integer.max)
  } else {
    .mp_lm(fcn, jac, x0, tol, tol, 0, 200L * (n + 1L))
  }
  rc <- if (max(abs(out$fvec)) <= tol) "Success" else "Failure"
  .wr_out(out$x, out$fvec, rc, -1L, cnt$f, cnt$g, cnt$g, cnt$g, original = list(info = out$info))
}

# ---------------------------------------------------------------------------
# NLsolve

#' @rdname nonlinear_wrappers
#' @export
NLsolveJL <- function(method = "trust_region", autodiff = NULL, linesearch = NULL, beta = 1, factor = 1,
                      autoscale = TRUE, m = 10L) {
  method <- match.arg(method, c("trust_region", "newton", "anderson", "broyden"))
  linesearch <- if (is.null(linesearch)) "Static" else match.arg(linesearch, c("Static", "BackTracking"))
  .wr_alg("NLsolveJL", .nls_run, method = method, autodiff = autodiff, linesearch = linesearch, beta = beta,
          factor = factor, autoscale = isTRUE(autoscale), m = as.integer(m))
}

# NLSolversBase caching: a residual or Jacobian is evaluated only at a new point
.nls_df <- function(F, J) {
  o <- new.env(parent = emptyenv())
  o$fc <- 0L; o$jc <- 0L; o$xf <- NULL; o$xdf <- NULL
  same <- function(a, b) !is.null(b) && length(a) == length(b) && all(a == b)
  o$value2 <- function(x) { o$F <- as.numeric(F(x)); o$xf <- x; o$fc <- o$fc + 1L; o$F }
  o$jac2 <- function(x) { o$J <- J(x); o$xdf <- x; o$jc <- o$jc + 1L; o$J }
  o$value <- function(x) { if (!same(x, o$xf)) o$value2(x); o$F }
  o$jacobian <- function(x) { if (!same(x, o$xdf)) o$jac2(x); o$J }
  o$vj2 <- function(x) { o$value2(x); o$jac2(x) }
  o$vj <- function(x) {
    nf <- !same(x, o$xf); nj <- !same(x, o$xdf)
    if (nf) o$value2(x)
    if (nj) o$jac2(x)
  }
  o
}

.nls_assess <- function(x, xold, f, xtol, ftol)
  c(isTRUE(.wr_l2(x - xold) <= xtol), max(abs(f)) <= ftol)

.nls_dogleg <- function(r, d, J, delta) {
  p_i <- tryCatch(.wr_ldiv(J, r), error = function(e) {
    s <- svd(J)
    k <- sum(s$d > .Machine$double.eps)
    dinv <- c(1 / s$d[seq_len(k)], rep(0, length(s$d) - k))
    as.numeric(s$v %*% (dinv * crossprod(s$u, r)))
  })
  p_i <- -p_i
  if (.wr_l2(d * p_i) <= delta) return(p_i)
  g <- as.numeric(crossprod(J, r)) / d^2
  p_c <- -.wr_l2(d * g)^2 / sum(as.numeric(J %*% g)^2) * g
  if (.wr_l2(d * p_c) >= delta) return(-delta / .wr_l2(d * g) * g)
  p_diff <- p_i - p_c
  b <- 2 * sum((d * p_c) * (d * p_diff))
  a <- .wr_l2(d * p_diff)^2
  tau <- (-b + sqrt(b^2 - 4 * a * (.wr_l2(d * p_c)^2 - delta^2))) / (2 * a)
  p_c + tau * p_diff
}

.nls_backtracking <- function(phi, phi0, dphi0) {
  c1 <- 1e-4; rho_hi <- 0.5; rho_lo <- 0.1; iterations <- 1000L
  iterfinitemax <- -log2(.Machine$double.eps)
  if (!(is.finite(phi0) && is.finite(dphi0))) stop("Value and slope at step length = 0 must be finite.", call. = FALSE)
  if (dphi0 > 0) stop("Search direction is not a direction of descent.", call. = FALSE)
  a1 <- a2 <- 1
  phix0 <- phi0
  phix1 <- phi(a2)
  k <- 0
  while (!is.finite(phix1) && k < iterfinitemax) { k <- k + 1; a2 <- a2 / 2; phix1 <- phi(a2) }
  if (!is.finite(phix1)) stop("Backtracking: failed to achieve finite new evaluation point.", call. = FALSE)
  iteration <- 0L
  while (phix1 > phi0 + c1 * a2 * dphi0) {
    iteration <- iteration + 1L
    if (iteration > iterations) stop("Linesearch failed to converge, reached maximum iterations ", iterations, call. = FALSE)
    if (iteration == 1L) {
      atmp <- -(dphi0 * a2^2) / (2 * (phix1 - phi0 - dphi0 * a2))
    } else {
      dv <- 1 / (a1^2 * a2^2 * (a2 - a1))
      aa <- (a1^2 * (phix1 - phi0 - dphi0 * a2) - a2^2 * (phix0 - phi0 - dphi0 * a1)) * dv
      bb <- (-a1^3 * (phix1 - phi0 - dphi0 * a2) + a2^3 * (phix0 - phi0 - dphi0 * a1)) * dv
      atmp <- if (aa == 0) -dphi0 / (2 * bb) else -dphi0 / (bb + sqrt(max(bb^2 - 3 * aa * dphi0, 0)))
    }
    a1 <- a2
    atmp <- if (is.nan(atmp)) a2 * rho_hi else min(atmp, a2 * rho_hi)
    a2 <- if (is.nan(atmp)) a2 * rho_lo else max(atmp, a2 * rho_lo)
    phinew <- phi(a2)
    while (!is.finite(phinew) && iteration < iterations) { iteration <- iteration + 1L; a2 <- a2 / 2; phinew <- phi(a2) }
    if (!is.finite(phinew)) stop("Backtracking: failed to achieve finite new evaluation point.", call. = FALSE)
    phix0 <- phix1; phix1 <- phinew
  }
  a2
}

.nls_run <- function(a, F, J, x0, opts, ls) {
  ftol <- .wr_tol(opts$abstol); xtol <- 0
  iterations <- as.integer(opts$maxiters %||% 1000L)
  df <- .nls_df(F, .wr_jac(F, J))
  x <- as.numeric(x0)
  # NLsolve's IsFiniteException
  fin <- function(v) if (any(!is.finite(v)))
    stop("During the resolution of the non-linear system, the evaluation of the following equation(s) ",
         "resulted in a non-finite number: [", paste(which(!is.finite(v)), collapse = ", "), "]", call. = FALSE)
  out <- function(x, it, conv, gc = df$jc)
    .wr_out(x, as.numeric(F(x)), if (any(conv)) "Success" else "Failure", it, df$fc, gc, gc, gc,
            original = list(x_converged = conv[1], f_converged = conv[2]))
  if (a$method == "trust_region") {
    df$vj2(x); r <- df$F; fin(r)
    conv <- .nls_assess(x, x, r, NaN, ftol)
    if (any(conv)) return(out(x, 0L, conv))
    J0 <- df$J
    d <- if (a$autoscale) { v <- sqrt(colSums(J0^2)); ifelse(v^2 == 0, 1, v) } else rep(1, length(x))
    delta <- a$factor * .wr_l2(d * x)
    if (delta == 0) delta <- a$factor
    it <- 0L; stopped <- any(is.nan(x)) || any(is.nan(r))
    while (!stopped && !any(conv) && it < iterations) {
      it <- it + 1L
      Jc <- df$J
      p <- .nls_dogleg(r, d, Jc, delta)
      xold <- x
      x <- x + p
      fx <- df$value(x)
      r_pred <- as.numeric(Jc %*% p) + r
      rho <- (sum(r^2) - sum(fx^2)) / (sum(r^2) - sum(r_pred^2))
      # comparisons with a NaN ratio are false, as in Julia
      if (isTRUE(rho > 1e-4)) {
        r <- fx
        Jn <- df$jacobian(x)
        if (a$autoscale) { d <- pmax(0.1 * d, sqrt(colSums(Jn^2))); d[d^2 == 0] <- 1 }
        conv <- .nls_assess(x, xold, r, xtol, ftol)
      } else {
        x <- x - p
        conv <- c(FALSE, FALSE)
      }
      if (isTRUE(rho < 0.1)) delta <- delta / 2
      else if (isTRUE(rho >= 0.9)) delta <- 2 * .wr_l2(d * p)
      else if (isTRUE(rho >= 0.5)) delta <- max(delta, 2 * .wr_l2(d * p))
      stopped <- any(is.nan(x)) || any(is.nan(df$F))
    }
    return(out(x, it, conv))
  }
  if (a$method == "newton") {
    n <- length(x)
    df$vj2(x); fin(df$F)
    conv <- .nls_assess(x, x, df$F, NaN, ftol)
    it <- 0L; stopped <- any(is.nan(x)) || any(is.nan(df$F))
    while (!stopped && !any(conv) && it < iterations) {
      it <- it + 1L
      if (it > 1L) df$vj(x)
      Jc <- df$J; Fc <- df$F
      p <- tryCatch(-.wr_ldiv(Jc, Fc), error = function(e) {
        fjac2 <- crossprod(Jc)
        lambda <- 1e6 * sqrt(n * .Machine$double.eps) * sum(abs(fjac2))
        .wr_ldiv(-(fjac2 + lambda * diag(n)), Fc)
      })
      xold <- x
      if (a$linesearch == "Static") {
        x <- x + p
        df$vj(x)
      } else {
        g <- as.numeric(crossprod(Jc, Fc))
        phi <- function(al) { v <- df$value(xold + al * p); sum(v * v) / 2 }
        al <- .nls_backtracking(phi, sum(Fc * Fc) / 2, sum(g * p))
        x <- xold + al * p
      }
      conv <- .nls_assess(x, xold, df$value(x), xtol, ftol)
      stopped <- any(is.nan(x)) || any(is.nan(df$F))
    }
    return(out(x, it, conv))
  }
  if (a$method == "anderson") {
    nn <- length(x); beta <- a$beta
    # AndersonCache keeps at most min(length(x), m) differences
    m <- min(nn, a$m)
    Q <- matrix(0, nn, m); Rm <- matrix(0, m, m); dgs <- vector("list", m)
    m_eff <- 0L; it <- 0L; conv <- c(FALSE, FALSE)
    fxold <- gold <- NULL
    givens <- function(f, g) {
      if (g == 0) return(c(1, 0, f))
      if (f == 0) return(c(0, 1, g))
      r <- sqrt(f * f + g * g); cs <- f / r; sn <- g / r
      if (abs(f) > abs(g) && cs < 0) { cs <- -cs; sn <- -sn; r <- -r }
      c(cs, sn, r)
    }
    qrdelete <- function(k) {
      if (k >= 2) for (i in 2:k) {
        gv <- givens(Rm[i - 1, i], Rm[i, i])
        r1 <- Rm[i - 1, ]; r2 <- Rm[i, ]
        Rm[i - 1, ] <<- gv[1] * r1 + gv[2] * r2
        Rm[i, ] <<- -gv[2] * r1 + gv[1] * r2
        q1 <- Q[, i - 1]; q2 <- Q[, i]
        Q[, i - 1] <<- gv[1] * q1 + gv[2] * q2
        Q[, i] <<- -gv[2] * q1 + gv[1] * q2
      }
      if (k >= 2) for (j in 1:(k - 1)) for (i in 1:(k - 1)) Rm[i, j] <<- Rm[i, j + 1]
    }
    qradd <- function(v, k) {
      if (k >= 2) for (i in 1:(k - 1)) {
        r <- .wr_dot(Q[, i], v)
        Rm[i, k] <<- r
        v <- .wr_fma(-r, Q[, i], v)
      }
      dd <- .wr_l2(v)
      Rm[k, k] <<- dd
      Q[, k] <<- v / dd
    }
    while (it < iterations) {
      it <- it + 1L
      fx <- df$value2(x)
      fin(fx)
      g <- x + beta * fx
      conv <- .nls_assess(g, x, fx, xtol, ftol)
      if (any(is.nan(x)) || any(is.nan(fx))) break
      if (any(conv)) break
      x <- g
      if (m > 0) {
        if (it == 1L) { fxold <- fx; gold <- g }
        else {
          m_eff <- m_eff + 1L
          if (m_eff > m) {
            dgs <- c(dgs[-1L], list(NULL))
            qrdelete(m)
            m_eff <- m
          }
          dgs[[m_eff]] <- g - gold
          qradd(fx - fxold, m_eff)
          fxold <- fx; gold <- g
          Rk <- .wr_upper(Rm, m_eff)
          dk <- svd(Rk)$d
          while (max(dk) / min(dk) > 1e10 && m_eff > 1L) {
            qrdelete(m_eff)
            m_eff <- m_eff - 1L
            Rk <- .wr_upper(Rm, m_eff)
            dk <- svd(Rk)$d
          }
          gam <- .wr_utrisolve(Rk, .wr_tmv(Q[, seq_len(m_eff), drop = FALSE], fx))
          for (i in seq_len(m_eff)) x <- x - gam[i] * dgs[[i]]
        }
      }
    }
    return(out(x, it, conv, 0L))
  }
  # broyden
  nn <- length(x)
  fv <- df$value2(x); fin(fv)
  xold <- numeric(nn)
  conv <- .nls_assess(x, xold, fv, NaN, ftol)
  Jinv <- diag(nn)
  it <- 0L; stopped <- any(is.nan(x)) || any(is.nan(fv))
  resets <- 0L
  while (!stopped && !any(conv) && it < iterations) {
    it <- it + 1L
    xold <- x; fold <- df$F
    p <- -as.numeric(Jinv %*% fold)
    x <- xold + p
    df$value(x)
    if (.wr_l2(df$F) <= 0.9 * .wr_l2(fold) - 0.001 * .wr_l2(p)^2) {
      al <- 1
    } else {
      Fn <- function(z) .wr_l2(df$value(z))
      fxn <- Fn(x)
      if (!is.finite(fxn)) stop("Exiting line search: 2-norm of current residuals is not finite.", call. = FALSE)
      lam <- 1
      # NLsolve's nantrack never succeeds, so a non finite trial is an error
      if (!is.finite(Fn(x - lam * p)))
        stop("Exiting line search: failed to find finite 2-norm of residuals in 5 trials.", call. = FALSE)
      ok <- Fn(x + lam * p) <= fxn - 0.001 * .wr_l2(lam * p)^2 + 0.1 * Fn(x)
      j <- 0L
      while (j < 50L && !ok) {
        j <- j + 1L
        lam <- 0.5 * lam
        Fn(x - lam * p)
        ok <- Fn(x + lam * p) <= fxn - 0.001 * .wr_l2(lam * p)^2 + 0.1 * Fn(x)
      }
      if (j >= 50L && !ok) stop("Exiting line search: failed to satisfy condition (2.4).", call. = FALSE)
      al <- lam
      x <- xold + al * p
    }
    df$value(x)
    dx <- x - xold; dfv <- df$F - fold
    stuck <- all(abs(dx) <= 1e-12) || all(abs(dfv) <= 1e-12)
    if (stuck) {
      Jinv <- diag(nn)
      resets <- resets + 1L
      if (resets > 3L) stuck <- FALSE
    } else {
      Jd <- as.numeric(Jinv %*% dfv)
      Jinv <- Jinv + outer((dx - Jd) / sum(dx * Jd), as.numeric(crossprod(dx, Jinv)))
    }
    if (!stuck) conv <- .nls_assess(x, xold, df$F, xtol, ftol)
    stopped <- any(is.nan(x)) || any(is.nan(df$F))
  }
  out(x, it, conv, 0L)
}

# ---------------------------------------------------------------------------
# LeastSquaresOptim

#' @rdname nonlinear_wrappers
#' @export
LeastSquaresOptimJL <- function(alg = "lm", linsolve = NULL, autodiff = "central") {
  alg <- match.arg(alg, c("lm", "dogleg"))
  if (!is.null(linsolve)) linsolve <- match.arg(linsolve, c("qr", "cholesky", "lsmr"))
  .wr_alg("LeastSquaresOptimJL", .lso_run, alg = alg, linsolve = linsolve, autodiff = autodiff)
}

.lso_run <- function(a, F, J, x0, opts, ls) {
  x_tol <- .wr_tol(opts$reltol); f_tol <- g_tol <- .wr_tol(opts$abstol)
  iterations <- as.integer(opts$maxiters %||% 1000L)
  jac <- .wr_jac(F, J)
  chol_solver <- identical(a$linsolve, "cholesky")
  x <- as.numeric(x0)
  fcur <- as.numeric(F(x)); f_calls <- 1L; g_calls <- 0L
  ssr <- sum(fcur^2); maxabs_gr <- Inf
  assess <- function(dx, ssr, trial, acc) {
    if (acc && abs(trial - ssr) <= f_tol * (abs(ssr) + f_tol)) c(FALSE, TRUE, FALSE)
    else if (max(abs(dx)) <= x_tol) c(TRUE, FALSE, FALSE)
    else if (maxabs_gr <= g_tol) c(FALSE, FALSE, TRUE)
    else c(FALSE, FALSE, FALSE)
  }
  lsq <- function(Jm, y, damp = NULL) {
    if (chol_solver) {
      A <- crossprod(Jm)
      if (!is.null(damp)) diag(A) <- diag(A) + damp
      b <- as.numeric(crossprod(Jm, y))
      if (is.null(damp)) {
        # cholesky!(Symmetric(J'J), RowMaximum()): pivoted, full rank required
        U <- chol(A, pivot = TRUE, tol = 0)
        if (attr(U, "rank") < ncol(A)) stop("RankDeficientException", call. = FALSE)
        piv <- attr(U, "pivot")
        z <- backsolve(U, forwardsolve(t(U), b[piv]))
        out <- numeric(length(z)); out[piv] <- z
        return(out)
      }
      U <- chol(A)
      return(backsolve(U, forwardsolve(t(U), b)))
    }
    A <- if (is.null(damp)) Jm else rbind(Jm, diag(sqrt(damp), length(damp)))
    rhs <- c(y, numeric(nrow(A) - length(y)))
    .wr_qr_ldiv(A, rhs)
  }
  iter <- 0L; conv <- c(FALSE, FALSE, FALSE)
  if (a$alg == "lm") {
    Delta <- 10; decrease <- 2; need_jac <- TRUE
    while (!any(conv) && iter < iterations) {
      iter <- iter + 1L
      if (need_jac) { Jm <- jac(x); g_calls <- g_calls + 1L; need_jac <- FALSE }
      dtd <- colSums(Jm^2)
      mn <- sum(dtd) / length(dtd)
      dtd <- pmin(pmax(dtd, 1e-6 * mn), 1e32 * mn) * (1 / Delta)
      dx <- lsq(Jm, fcur, dtd)
      maxabs_gr <- max(abs(.wr_tmv(Jm, fcur)))
      x <- x - dx
      ftrial <- as.numeric(F(x)); f_calls <- f_calls + 1L
      trial_ssr <- sum(ftrial^2)
      fpred <- .wr_tmv(t(Jm), dx) - fcur
      pr <- abs(ssr - sum(fpred^2))
      rho <- if (pr > 0) (ssr - trial_ssr) / pr else 0
      acc <- rho > 1e-3
      conv <- assess(dx, ssr, trial_ssr, acc)
      if (acc) {
        fcur <- ftrial; ssr <- trial_ssr
        Delta <- min(Delta / max(1 / 3, 1 - (2 * rho - 1)^3), 1e16)
        decrease <- 2; need_jac <- TRUE
      } else {
        x <- x + dx
        Delta <- max(Delta / decrease, 1e-16)
        decrease <- decrease * 2
      }
    }
  } else {
    Delta <- 1; reuse <- FALSE
    while (!any(conv) && iter < iterations) {
      iter <- iter + 1L
      if (!reuse) {
        Jm <- jac(x); g_calls <- g_calls + 1L
        dtd <- pmin(pmax(colSums(Jm^2), 1e-6), 1e32)
        if (iter == 1L) { wx <- sqrt(sum(x * dtd * x)); if (wx > 0) Delta <- Delta * wx }
        dgr <- as.numeric(crossprod(Jm, fcur))
        maxabs_gr <- max(abs(dgr))
        dgr <- dgr / dtd
        wdgr <- sqrt(sum(dgr * dtd * dgr))
        alpha <- wdgr^2 / sum(as.numeric(Jm %*% dgr)^2)
        dgn <- lsq(Jm, fcur)
        wdgn <- sqrt(sum(dgn * dtd * dgn))
      }
      if (wdgn <= Delta) { dx <- dgn; wdx <- wdgn }
      else if (wdgr * alpha >= Delta) { dx <- dgr * (Delta / wdgr); wdx <- Delta }
      else {
        bda <- alpha * sum(dgr * dtd * dgn)
        asq <- (alpha * wdgr)^2
        bmasq <- asq - 2 * bda + wdgn^2
        cc <- bda - asq
        dd <- sqrt(cc^2 + bmasq * (Delta^2 - asq))
        beta <- if (cc <= 0) (dd - cc) / bmasq else (Delta^2 - asq) / (dd + cc)
        dx <- dgn * beta + alpha * (1 - beta) * dgr
        wdx <- sqrt(sum(dx * dtd * dx))
      }
      x <- x - dx
      ftrial <- as.numeric(F(x)); f_calls <- f_calls + 1L
      trial_ssr <- sum(ftrial^2)
      fpred <- .wr_tmv(t(Jm), dx) - fcur
      pr <- abs(ssr - sum(fpred^2))
      rho <- if (pr > 0) (ssr - trial_ssr) / pr else 0
      acc <- rho >= 1e-3
      conv <- assess(dx, ssr, trial_ssr, acc)
      if (acc) { reuse <- FALSE; fcur <- ftrial; ssr <- trial_ssr }
      else { reuse <- TRUE; x <- x + dx }
      if (rho < 0.25) Delta <- max(1e-16, Delta * 0.5)
      else if (rho > 0.75) Delta <- max(Delta, 3 * wdx)
    }
  }
  rc <- if (any(conv)) "Success" else if (iter >= iterations) "MaxIters" else "ConvergenceFailure"
  .wr_out(x, as.numeric(F(x)), rc, iter, f_calls, g_calls, -1L, -1L, original = list(ssr = ssr))
}

# ---------------------------------------------------------------------------
# FastLevenbergMarquardt

#' @rdname nonlinear_wrappers
#' @export
FastLevenbergMarquardtJL <- function(linsolve = "cholesky", factor = 1e-6, factoraccept = 13, factorreject = 3,
                                     factorupdate = "marquardt", minscale = 1e-12, maxscale = 1e16, minfactor = 0x1.fb0f6be506019p-94,
                                     maxfactor = 1e32, autodiff = NULL) {
  linsolve <- match.arg(linsolve, c("cholesky", "qr"))
  factorupdate <- match.arg(factorupdate, c("marquardt", "nielson", "nielsen"))
  # NonlinearSolve accepts only :nielson and passes it on unchanged, so
  # FastLevenbergMarquardt (which tests :nielsen) never uses the Nielsen update
  if (factorupdate == "nielsen") factorupdate <- "nielson"
  .wr_alg("FastLevenbergMarquardtJL", .flm_run, linsolve = linsolve, factor = factor, factoraccept = factoraccept,
          factorreject = factorreject, factorupdate = factorupdate, minscale = minscale, maxscale = maxscale,
          minfactor = minfactor, maxfactor = maxfactor, autodiff = autodiff)
}

# the unrolled solve of FastLevenbergMarquardt for at most 12 unknowns:
# Gaussian elimination of J'J + lambda DtD (lower triangle, no pivoting),
# refused when a diagonal entry is negative
.flm_ldl_solve <- function(JtJ, lam, DtD, g) {
  n <- length(g)
  A <- JtJ
  # every update is a muladd, fused on hardware with fma
  for (k in seq_len(n)) {
    A[k, k] <- .wr_fma(lam, DtD[k], JtJ[k, k])
    if (!(A[k, k] >= 0)) return(NULL)
  }
  x <- g
  for (k in seq_len(n)) {
    if (k > 1L) for (j in 1:(k - 1L)) for (m in k:n) A[m, k] <- .wr_fma(-A[k, j], A[m, j] / A[j, j], A[m, k])
  }
  for (j in seq_len(n)) if (j < n) for (i in (j + 1L):n) x[i] <- .wr_fma(-x[j], A[i, j] / A[j, j], x[i])
  for (j in n:1) {
    x[j] <- x[j] / A[j, j]
    if (j > 1L) for (i in (j - 1L):1) x[i] <- .wr_fma(-x[j], A[j, i], x[i])
  }
  x
}

# qrsolve! of FastLevenbergMarquardt (after MINPACK): solve J p = f,
# sqrt(lambda) D p = 0 from the pivoted QR of J, Givens rotations by hypot
.flm_qrsolve <- function(R, qtf, lam, D, jpvt) {
  n <- ncol(R)
  wf <- qtf[seq_len(n)]
  Jd <- R
  p <- numeric(n)
  for (j in seq_len(n)) {
    p[j] <- R[j, j]
    if (j < n) for (i in (j + 1L):n) Jd[i, j] <- Jd[j, i]
  }
  ld <- numeric(n)
  sl <- sqrt(lam)
  for (j in seq_len(n)) {
    ld[j] <- sl * D[jpvt[j]]
    if (j < n) ld[(j + 1L):n] <- 0
    fpj <- 0
    for (k in j:n) {
      if (ld[k] == 0) next
      rr <- .wr_hypot(Jd[k, k], ld[k])
      cs <- Jd[k, k] / rr; sn <- ld[k] / rr
      Jd[k, k] <- rr
      wfk <- wf[k]
      wf[k] <- cs * wfk + sn * fpj
      fpj <- -sn * wfk + cs * fpj
      if (k < n) for (i in (k + 1L):n) {
        jik <- Jd[i, k]; ldi <- ld[i]
        Jd[i, k] <- cs * jik + sn * ldi
        ld[i] <- -sn * jik + cs * ldi
      }
    }
    ld[j] <- Jd[j, j]
    Jd[j, j] <- p[j]
  }
  nsing <- n
  out <- numeric(n)
  for (i in seq_len(n)) if (ld[i] == 0) {
    wf[i:n] <- 0
    nsing <- i - 1L
    break
  }
  if (nsing >= 1L) for (j in nsing:1) {
    fj <- wf[j]
    if (j < nsing) for (i in (j + 1L):nsing) fj <- fj - wf[i] * Jd[i, j]
    wf[j] <- fj / ld[j]
  }
  out[jpvt[seq_len(nsing)]] <- wf[seq_len(nsing)]
  out
}

.flm_run <- function(a, F, J, x0, opts, ls) {
  ec <- ex <- .wr_tol(opts$reltol); eg <- .wr_tol(opts$abstol)
  maxit <- as.integer(opts$maxiters %||% 1000L)
  jac <- .wr_jac(F, J)
  lam <- a$factor; lrej <- a$factorreject; lacc <- 1 / a$factoraccept
  nu <- lrej
  x <- as.numeric(x0)
  f <- as.numeric(F(x)); nfev <- 1L; njev <- 0L
  Fs <- sum(f^2)
  done <- function(info, iter) {
    .wr_out(x, Fs, if (info == -1L) "MaxIters" else "Success", iter, nfev, njev, -1L, -1L,
            original = list(info = info))
  }
  if (Fs < ec) return(done(1L, 0L))
  Jm <- jac(x); njev <- 1L
  g <- as.numeric(crossprod(Jm, f))
  if (max(abs(g)) < eg) return(done(3L, 0L))
  DtD <- pmax(.wr_colsumsq(Jm), a$minscale)
  # init! and update! of the linear solver: J'J for Cholesky, the pivoted QR
  # of J with Q'f and sqrt(DtD) for QR
  fact <- NULL
  prepare <- function() {
    fact <<- if (a$linsolve == "cholesky") list(JtJ = crossprod(Jm)) else {
      qj <- qr(Jm, LAPACK = TRUE)
      list(R = qj$qr[seq_len(ncol(Jm)), , drop = FALSE], piv = qj$pivot, qtf = .wr_qty(qj, f),
           D = sqrt(DtD))
    }
  }
  solve_step <- function() {
    if (a$linsolve == "cholesky") {
      if (length(g) <= 12L) return(.flm_ldl_solve(fact$JtJ, lam, DtD, g))
      A <- fact$JtJ
      diag(A) <- diag(A) + lam * DtD
      U <- tryCatch(chol(A), error = function(e) NULL)
      if (is.null(U)) return(NULL)
      return(backsolve(U, forwardsolve(t(U), g)))
    }
    .flm_qrsolve(fact$R, fact$qtf, lam, fact$D, fact$piv)
  }
  prepare()
  iter <- 0L; info <- -1L
  while (iter < maxit) {
    iter <- iter + 1L
    accepted <- FALSE
    p <- solve_step()
    if (!is.null(p)) {
      xk <- x - p
      fk <- as.numeric(F(xk)); Fk <- sum(fk^2); nfev <- nfev + 1L
      rho <- 0
      if (Fk < Fs) {
        ac <- Fs - Fk
        pr <- sum(p * (lam * DtD * p + g))
        rho <- ac / pr
        accepted <- TRUE
      }
      if (all(abs(p) < ex * (ex + abs(x)))) {
        if (accepted) { x <- xk; f <- fk; Fs <- Fk }
        info <- 2L
        break
      }
    }
    if (accepted) {
      Fold <- Fs
      x <- xk; f <- fk; Fs <- Fk
      if (abs(ac) < ec * Fold && pr < ec * Fold && rho < 2) { info <- 1L; break }
      Jm <- jac(x); njev <- njev + 1L
      g <- as.numeric(crossprod(Jm, f))
      if (max(abs(g)) < eg) { info <- 3L; break }
      DtD <- pmin(pmax(.wr_colsumsq(Jm), DtD), a$maxscale)
      prepare()
      if (a$factorupdate == "nielson") {
        nu <- 2 * rho - 1
        lam <- lam * max(lacc, 1 - (lrej - 1) * nu^3)
        nu <- lrej
      } else if (rho > 0.75) lam <- lam * lacc
      else if (rho < 1e-3) lam <- lam * lrej
    } else {
      if (a$factorupdate == "nielson") { lam <- lam * nu; nu <- nu * 2 } else lam <- lam * lrej
    }
    lam <- min(max(lam, a$minfactor), a$maxfactor)
  }
  done(info, iter)
}

# ---------------------------------------------------------------------------
# NLSolvers

#' @rdname nonlinear_wrappers
#' @export
NLSolversJL <- function(method = "LineSearch", autodiff = NULL) {
  method <- match.arg(method, "LineSearch")
  .wr_alg("NLSolversJL", .nlsolvers_run, method = method, autodiff = autodiff)
}

.nlsolvers_run <- function(a, F, J, x0, opts, ls) {
  abstol <- .wr_tol(opts$abstol)
  maxiter <- as.integer(opts$maxiters %||% 1000L)
  jac <- .wr_jac(F, J)
  x <- as.numeric(x0)
  Fx <- as.numeric(F(x)); Jx <- jac(x)
  res <- function(z, Fx, it) {
    rc <- if (max(abs(Fx)) <= abstol) "Success" else "MaxIters"
    .wr_out(z, Fx, rc, it, -1L, -1L, -1L, -1L)
  }
  if (max(abs(Fx)) < abstol) return(res(x, Fx, 0L))
  z <- x
  iter <- 1L
  while (iter <= maxiter) {
    x <- z
    d <- .wr_ldiv(Jx, -Fx)
    z <- x + 1 * d
    Fx <- as.numeric(F(z)); Jx <- jac(z)
    rF <- max(abs(Fx))
    if (!is.finite(rF)) break
    if (rF < abstol) break
    iter <- iter + 1L
  }
  res(z, Fx, iter)
}

# ---------------------------------------------------------------------------
# SpeedMapping

#' @rdname nonlinear_wrappers
#' @export
SpeedMappingJL <- function(sigma_min = 0, stabilize = FALSE, check_obj = FALSE, orders = c(3L, 3L, 2L)) {
  orders <- as.integer(orders)
  if (min(orders) < 2L || max(orders) > 3L) stop("orders must be 2 or 3", call. = FALSE)
  if (min(orders) != 2L) stop("orders must contain at least one 2", call. = FALSE)
  .wr_alg("SpeedMappingJL", .sm_run, sigma_min = sigma_min, stabilize = isTRUE(stabilize),
          check_obj = isTRUE(check_obj), orders = orders)
}

# the ACX loop of SpeedMapping for a map (alpha stays 1, pnorm = Inf)
.sm_run <- function(a, F, J, x0, opts, ls) {
  abstol <- .wr_tol(opts$abstol); abstolsq <- abstol^2
  maps_limit <- as.integer(opts$maxiters %||% 1000L)
  m <- function(x) x + as.numeric(F(x))
  orders <- a$orders
  reltol_rr <- 100^2
  lpsq <- function(r) max(abs(r))^2
  x_now <- as.numeric(x0); x_next <- x_best <- x_now
  r1 <- r2 <- r3 <- numeric(length(x_now))
  maps <- 0L; iter <- 0L; i_ord <- 0L
  rr_best <- lp_best <- Inf; al <- 1; sg <- 1
  start_map <- 1L
  while (maps < maps_limit) {
    iter <- iter + 1L
    i_ord <- i_ord %% length(orders) + 1L
    p <- orders[i_ord]
    lp_best <- lpsq(r1)
    i <- start_map
    while (i <= p) {
      x_next <- m(x_now)
      r_now <- x_next - x_now
      rr_now <- .wr_dot(r_now, r_now)
      maps <- maps + 1L
      if (i == 1L) r1 <- r_now else if (i == 2L) r2 <- r_now else r3 <- r_now
      if (rr_now <= rr_best) {
        rr_best <- rr_now
        lp_best <- lpsq(r_now)
        tmp <- x_best; x_best <- x_now; x_now <- tmp
        if (lp_best < abstolsq) break
      }
      tmp <- x_now; x_now <- x_next; x_next <- tmp
      i <- i + 1L
    }
    start_map <- 2L
    if (lp_best < abstolsq) break
    if (p == 3L) r3 <- r3 + (-2 * r2 + r1)
    r2 <- r2 - r1
    sg <- -abs(if (p == 2L) .wr_dot(r2, r1) / .wr_dot(r2, r2) else .wr_dot(r3, r2) / .wr_dot(r3, r3))
    if (!is.finite(sg)) {
      al <- al * 5
      sg <- 1
      x_now <- m(x_best); maps <- maps + 1L
      r1 <- x_now - x_best
      rr_best <- .wr_dot(r1, r1)
    } else {
      # grouped as the broadcast x += a r1 + b r2 (+ c r3); sg^3 as sg * sg * sg
      x_now <- if (p == 2L) x_now + ((-2 * al * (1 + sg)) * r1 + (al * (-1 + sg * sg)) * r2)
               else x_now + ((-3 * al * (1 + sg)) * r1 + (-3 * al * (1 - sg * sg)) * r2 + (-al * (1 + sg * sg * sg)) * r3)
      damp <- 1
      for (attempt in 1:100) {
        x_next <- m(x_now)
        r1 <- x_next - x_now
        rr_now <- .wr_dot(r1, r1)
        maps <- maps + 1L
        if (rr_now < rr_best * reltol_rr) break
        damp <- damp / 2
        x_now <- x_best + damp * (x_now - x_best)
      }
      rr_best <- rr_now
      tmp <- x_best; x_best <- x_now; x_now <- x_next; x_next <- tmp
    }
  }
  conv <- lp_best < abstolsq
  .wr_out(x_best, as.numeric(F(x_best)), if (conv) "Success" else "Failure", maps, maps, 0L, 0L, 0L,
          original = list(maps = maps, iterations = iter))
}

# ---------------------------------------------------------------------------
# FixedPointAcceleration

#' @rdname nonlinear_wrappers
#' @export
FixedPointAccelerationJL <- function(algorithm = "Anderson", m = NULL, condition_number_threshold = NULL,
                                     extrapolation_period = NULL, replace_invalids = "NoAction", dampening = 1) {
  algorithm <- match.arg(algorithm, c("Anderson", "MPE", "RRE", "VEA", "SEA", "Simple", "Aitken", "Newton"))
  replace_invalids <- match.arg(replace_invalids, c("NoAction", "ReplaceInvalids", "ReplaceVector"))
  if (algorithm != "Anderson") {
    if (!is.null(condition_number_threshold))
      stop("`condition_number_threshold` is only valid for Anderson acceleration", call. = FALSE)
    if (!is.null(m)) stop("`m` is only valid for Anderson acceleration", call. = FALSE)
  }
  if (!algorithm %in% c("MPE", "RRE", "VEA", "SEA") && !is.null(extrapolation_period))
    stop("`extrapolation_period` is only valid for MPE, RRE, VEA and SEA", call. = FALSE)
  if (is.null(extrapolation_period)) extrapolation_period <- if (algorithm %in% c("SEA", "VEA")) 6L else 7L
  else if (algorithm %in% c("SEA", "VEA") && extrapolation_period %% 2 != 0)
    stop("`extrapolation_period` must be multiples of 2 for SEA and VEA", call. = FALSE)
  .wr_alg("FixedPointAccelerationJL", .fpa_run, algorithm = algorithm, m = as.integer(m %||% 10L),
          condition_number_threshold = condition_number_threshold %||% 1e3,
          extrapolation_period = as.integer(extrapolation_period), replace_invalids = replace_invalids,
          dampening = dampening)
}

.fpa_cond <- function(A) { d <- svd(A)$d; if (!length(d)) NaN else max(d) / min(d) }

.fpa_together <- function(Inp, Out) {
  k <- ncol(Inp)
  if (k == 1L) return(cbind(Inp, Out))
  dif <- colSums(Inp[, -1L, drop = FALSE] - Out[, -k, drop = FALSE])
  ag <- abs(dif) < 1e-10
  if (all(ag)) return(cbind(Inp, Out[, k]))
  lb <- max(which(!ag))
  cbind(Inp[, (lb + 1L):k, drop = FALSE], Out[, k])
}

.fpa_epsilon <- function(It, alg) {
  if (ncol(It) == 1L) return(It[, 1L])
  if (ncol(It) %% 2 == 0) It <- It[, -1L, drop = FALSE]
  Mat <- It
  prev <- matrix(0, nrow(Mat), ncol(Mat) - 1L)
  for (mc in rev(2:ncol(Mat))) {
    Dm <- Mat[, 2:mc, drop = FALSE] - Mat[, 1:(mc - 1L), drop = FALSE]
    inv <- if (nrow(Dm) == 1L || alg == "SEA") 1 / Dm
           else apply(Dm, 2, function(v) v / sum(v * v))
    inv <- matrix(inv, nrow(Dm))
    newm <- prev + inv
    prev <- Mat[, seq_len(mc - 1L)[-1L], drop = FALSE]
    Mat <- newm
  }
  if (any(is.nan(Mat))) return(.fpa_epsilon(It[, -(1:2), drop = FALSE], alg))
  Mat[, 1L]
}

.fpa_new_input <- function(a, Inp, Out, simple_start) {
  k <- ncol(Out)
  simple <- Out[, k]
  prop <- simple
  alg <- a$algorithm
  if (alg == "Anderson") {
    if (k >= 2L) {
      M <- min(a$m - 1L, k - 1L, nrow(Out))
      ro <- Out[, (k - M):k, drop = FALSE]; ri <- Inp[, (k - M):k, drop = FALSE]
      res <- ro - ri
      dO <- ro[, 2:(M + 1L), drop = FALSE] - ro[, 1:M, drop = FALSE]
      dR <- res[, 2:(M + 1L), drop = FALSE] - res[, 1:M, drop = FALSE]
      lastR <- res[, M + 1L]; lastO <- ro[, M + 1L]
      co <- rep(NaN, ncol(dO))
      while (any(is.nan(co))) {
        if (!length(dR)) break
        if (.fpa_cond(dR) > a$condition_number_threshold) {
          M <- M - 1L
          if (M < 1L) { co <- numeric(); break }
          dO <- dO[, 2:(M + 1L), drop = FALSE]; dR <- dR[, 2:(M + 1L), drop = FALSE]
          co <- rep(NaN, ncol(dO))
          next
        }
        co <- tryCatch(.wr_ldiv(dR, lastR), error = function(e) rep(NaN, ncol(dR)))
        if (any(is.nan(co))) {
          M <- M - 1L
          if (M < 1.5) break
          dO <- dO[, 2:(M + 1L), drop = FALSE]; dR <- dR[, 2:(M + 1L), drop = FALSE]
          co <- rep(NaN, ncol(dO))
        }
      }
      prop <- if (!length(co)) rep(NaN, length(simple)) else lastO - a$dampening * as.numeric(dO %*% co)
    }
  } else if (alg == "Aitken") {
    if ((k + simple_start) %% 3 == 0) {
      x <- Inp[, k - 1L]; fx <- Out[, k - 1L]; ffx <- Out[, k]
      prop <- x - ((fx - x)^2 / (ffx - 2 * fx + x))
    }
  } else if (alg == "Newton") {
    if ((k + simple_start) %% 2 == 1 && k > 1L) {
      xk1 <- Inp[, k - 1L]; gxk1 <- Out[, k - 1L] - xk1
      xk <- Inp[, k]; gxk <- Out[, k] - xk
      prop <- xk - gxk / ((gxk - gxk1) / (xk - xk1))
    }
  } else if (alg %in% c("MPE", "RRE")) {
    It <- .fpa_together(Inp, Out)
    if (ncol(It) %% a$extrapolation_period == 0) {
      tc <- ncol(It)
      if (alg == "MPE") {
        od <- It[, 2:(tc - 1L), drop = FALSE] - It[, 1:(tc - 2L), drop = FALSE]
        ld <- It[, tc] - It[, tc - 1L]
        cv <- c(-as.numeric(.wr_pinv(od) %*% ld), 1)
        prop <- as.numeric(It[, 2:tc, drop = FALSE] %*% cv) / sum(cv)
      } else {
        D <- It[, 2:tc, drop = FALSE] - It[, 1:(tc - 1L), drop = FALSE]
        D2 <- D[, 2:(tc - 1L), drop = FALSE] - D[, 1:(tc - 2L), drop = FALSE]
        prop <- It[, 1L] - as.numeric((D[, 1:(tc - 2L), drop = FALSE] %*% .wr_pinv(D2)) %*% D[, 1L])
      }
    }
  } else if (alg %in% c("VEA", "SEA")) {
    It <- .fpa_together(Inp, Out)
    if (ncol(It) %% a$extrapolation_period == 0) prop <- .fpa_epsilon(It, alg)
  }
  bad <- !is.finite(prop)
  if (any(bad)) {
    if (a$replace_invalids == "ReplaceInvalids") prop[bad] <- simple[bad]
    else if (a$replace_invalids == "ReplaceVector") prop <- simple
  }
  if (a$dampening == 1) return(prop)
  a$dampening * prop + (1 - a$dampening) * simple
}

.fpa_run <- function(a, F, J, x0, opts, ls) {
  tol <- .wr_tol(opts$abstol)
  maxit <- as.integer(opts$maxiters %||% 1000L)
  g <- function(x) x + as.numeric(F(x))
  safe <- function(x) {
    if (any(!is.finite(x))) return(NULL)
    y <- tryCatch(g(x), error = function(e) NULL)
    if (is.null(y) || any(!is.finite(y)) || length(y) != length(x)) return(NULL)
    y
  }
  Inp <- matrix(as.numeric(x0), ncol = 1L)
  y <- safe(Inp[, 1L])
  fail <- function() .wr_out(as.numeric(x0), as.numeric(F(as.numeric(x0))), "Failure", ncol(Inp) - 1L,
                             ncol(Inp) - 1L, 0L, 0L, 0L)
  if (is.null(y)) return(fail())
  Out <- matrix(y, ncol = 1L)
  conv <- max(abs(y - Inp[, 1L]))
  it <- 2L
  while (conv > tol && it <= maxit) {
    xin <- .fpa_new_input(a, Inp, Out, 0L)
    y <- safe(xin)
    if (is.null(y)) return(.wr_out(as.numeric(x0), as.numeric(F(as.numeric(x0))), "Failure", ncol(Out), ncol(Out),
                                   0L, 0L, 0L))
    Inp <- cbind(Inp, xin); Out <- cbind(Out, y)
    conv <- max(abs(y - xin))
    it <- it + 1L
  }
  k <- ncol(Out)
  if (!(conv < tol)) return(.wr_out(as.numeric(x0), as.numeric(F(as.numeric(x0))), "Failure", k, k, 0L, 0L, 0L))
  res <- Out[, k]
  resid <- as.numeric(F(res))
  .wr_out(res, resid, if (max(abs(resid)) <= tol) "Success" else "Failure", k, k, 0L, 0L, 0L)
}

# ---------------------------------------------------------------------------
# SIAMFANLEquations

#' @rdname nonlinear_wrappers
#' @export
SIAMFANLEquationsJL <- function(method = "newton", delta = 1e-3, linsolve = NULL, m = 0L, beta = 1, autodiff = NULL) {
  method <- match.arg(method, c("newton", "pseudotransient", "secant", "anderson"))
  if (!is.null(linsolve)) linsolve <- match.arg(linsolve, c("gmres", "bicgstab"))
  .wr_alg("SIAMFANLEquationsJL", .sf_run, method = method, delta = delta, linsolve = linsolve, m = as.integer(m),
          beta = beta, autodiff = autodiff)
}

.sf_armijo <- function(Fv, x, d, residm, armmax) {
  alpha <- 1e-4; iarm <- -1L; lambda <- 1; lamc <- lambda
  xt <- x + lambda * d; ft <- Fv(xt); residt <- .wr_l2(ft)
  armfail <- residt > (1 - alpha * lambda) * residm
  iarm <- iarm + 1L
  ffc <- residt^2; ff0 <- residm^2; ffm <- ffc
  while (armfail && iarm < armmax) {
    if (iarm == 0L) lambda <- lambda * 0.5
    else {
      lc <- lambda; lm <- lamc
      c2 <- lm * (ffc - ff0) - lc * (ffm - ff0)
      if (c2 >= 0) lambda <- 0.1 * lc
      else {
        c1 <- lc * lc * (ffm - ff0) - lm * lm * (ffc - ff0)
        lambda <- max(0.1 * lc, min(0.5 * lc, -c1 * 0.5 / c2))
      }
    }
    xt <- x + lambda * d; ft <- Fv(xt); residt <- .wr_l2(ft)
    ffm <- ffc; ffc <- residt^2
    iarm <- iarm + 1L
    armfail <- residt > (1 - alpha * lambda) * residm
  }
  list(x = xt, f = ft, resnorm = residt, iarm = iarm, idid = !(iarm >= armmax))
}

.sf_newton_err <- function(resnorm, iline, tol, toosoon, itc, maxit) {
  if (!(resnorm > tol || toosoon)) return(0L)
  code <- 0L
  if (toosoon) code <- -1L
  if (itc >= maxit) code <- 10L
  if (iline) code <- 1L
  code
}

# dirder of SIAMFANLEquations: J v by a forward difference of size 1e-7 / |v|,
# plus v / delta for pseudo transient continuation
.sf_dirder <- function(Fv, xc, FS, v, delta) {
  h <- 1e-7 / .wr_l2(v)
  hinv <- 1 / h
  out <- (-hinv) * FS + hinv * Fv(h * v + xc)
  if (delta != 0) out <- out + (1 / delta) * v
  out
}

# kl_gmres (classical Gram Schmidt twice, Givens rotations, zero initial
# iterate, at most 2 iterations) and kl_bicgstab (at most 5 iterations)
.sf_krylov <- function(atv, b, eta, lsolver) {
  n <- length(b)
  if (lsolver == "gmres") {
    kmax <- 2L
    rho <- .wr_l2(b)
    if (rho == 0) stop("Initial resdiual in kl_gmres is zero. Why?", call. = FALSE)
    errtol <- eta * .wr_l2(b)
    V <- matrix(0, n, kmax + 1L); h <- matrix(0, kmax + 1L, kmax + 1L)
    cs <- sn <- g <- numeric(kmax + 1L)
    g[1L] <- rho
    V[, 1L] <- b * (1 / rho)
    k <- 0L
    while (rho > errtol && k < kmax) {
      k <- k + 1L
      vv <- atv(V[, k])
      Q <- V[, seq_len(k), drop = FALSE]
      rk <- as.numeric(crossprod(Q, vv))
      vv <- vv - as.numeric(Q %*% rk)
      pk <- as.numeric(crossprod(Q, vv))
      vv <- vv - as.numeric(Q %*% pk)
      rk <- rk + pk
      nq <- .wr_l2(vv)
      if (nq > 0) vv <- vv / nq
      h[seq_len(k), k] <- rk; h[k + 1L, k] <- nq
      V[, k + 1L] <- vv
      if (k > 1L) for (i in seq_len(k - 1L)) {
        w1 <- cs[i] * h[i, k] - sn[i] * h[i + 1L, k]
        w2 <- sn[i] * h[i, k] + cs[i] * h[i + 1L, k]
        h[i, k] <- w1; h[i + 1L, k] <- w2
      }
      nu <- .wr_l2(h[k:(k + 1L), k])
      if (nu != 0) {
        cs[k] <- h[k, k] / nu; sn[k] <- -h[k + 1L, k] / nu
        h[k, k] <- cs[k] * h[k, k] - sn[k] * h[k + 1L, k]
        h[k + 1L, k] <- 0
        w1 <- cs[k] * g[k] - sn[k] * g[k + 1L]
        w2 <- sn[k] * g[k] + cs[k] * g[k + 1L]
        g[k] <- w1; g[k + 1L] <- w2
      }
      rho <- abs(g[k + 1L])
    }
    y <- backsolve(h[seq_len(k), seq_len(k), drop = FALSE], g[seq_len(k)])
    return(list(sol = as.numeric(V[, seq_len(k), drop = FALSE] %*% y), lits = k, idid = rho <= errtol))
  }
  lmaxit <- 5L
  r <- b; x <- numeric(n)
  rho <- numeric(lmaxit + 2L); rho[1L] <- 1; rho[2L] <- sum(r * r)
  alpha <- 1; omega <- 1
  r0 <- r; rnorm <- .wr_l2(r0)
  v <- p <- s <- t <- numeric(n)
  tol <- eta * .wr_l2(b)
  k <- 0L
  while (rnorm > tol && k < lmaxit) {
    k <- k + 1L
    beta <- (rho[k + 1L] / rho[k]) * (alpha / omega)
    p <- -omega * v + p
    p <- r + beta * p
    v <- atv(p)
    tau <- sum(r0 * v)
    alpha <- rho[k + 1L] / tau
    s <- -alpha * v + r
    t <- atv(s)
    omega <- sum(t * s) / sum(t * t)
    rho[k + 2L] <- -omega * sum(r0 * t)
    r <- -omega * t + s
    t <- alpha * p + omega * s
    x <- x + t
    rnorm <- .wr_l2(r)
  }
  list(sol = x, lits = k, idid = rnorm <= tol)
}

.sf_run <- function(a, F, J, x0, opts, ls) {
  atol <- .wr_tol(opts$abstol); rtol <- .wr_tol(opts$reltol)
  maxit <- as.integer(opts$maxiters %||% 1000L)
  Fv <- function(x) as.numeric(F(x))
  x <- as.numeric(x0); n <- length(x)
  map_rc <- function(code) switch(as.character(code), "0" = "Success", "10" = "MaxIters", "1" = "Failure",
                                  "-1" = "Default", "-2" = "Unstable")
  jac <- if (!is.null(J)) function(x, fx) as.matrix(J(x)) else function(x, fx) {
    h <- 1e-7 * max(abs(x)) + 1e-8
    JJ <- matrix(0, n, n)
    for (ic in seq_len(n)) { y <- x; y[ic] <- y[ic] + h; JJ[, ic] <- (Fv(y) - fx) / h }
    JJ
  }
  if (a$method == "anderson") {
    m <- a$m; beta <- a$beta
    G <- function(s) { gx <- s + Fv(s); if (beta != 1) gx <- (1 - beta) * s + beta * gx; gx }
    sol <- x
    gx <- G(sol); res <- gx - sol
    resnorm <- .wr_l2(res); up <- 1e4 * resnorm
    tol <- rtol * resnorm + atol
    toosoon <- resnorm <= tol
    k <- 0L
    QP <- matrix(0, n, max(m, 1L)); RP <- matrix(0, max(m, 1L), max(m, 1L)); DG <- matrix(0, n, max(m, 1L))
    point <- function() {
      dg <- -gx; gx <<- G(sol); dg <- dg + gx
      df <- -res; res <<- gx - sol; df <- df + res
      resnorm <<- .wr_l2(res)
      list(dg = dg, df = df)
    }
    if (!toosoon) { sol <- gx; k <- 1L; pt <- point() }
    while (k < maxit && resnorm > tol && !toosoon && resnorm < up) {
      if (m == 0L) sol <- gx
      else {
        kk <- k + 1L
        if (m == 1L) DG[, 1L] <- pt$dg
        else if (kk > m + 1L) { DG[, 1:(m - 1L)] <- DG[, 2:m]; DG[, m] <- pt$dg }
        else DG[, kk - 1L] <- pt$dg
        kq <- k - 1L
        if (kq == 0L) { RP[1, 1] <- .wr_l2(pt$df); QP[, 1L] <- pt$df / .wr_l2(pt$df) }
        else {
          if (kq > m - 1L) {
            qd <- qr(RP[, 2:m, drop = FALSE])
            Rd <- qr.R(qd); Qx <- qr.Q(qd)
            RP[1:(m - 1L), 1:(m - 1L)] <- Rd[1:(m - 1L), , drop = FALSE]
            RP[, m] <- 0
            QP[, 1:(m - 1L)] <- QP %*% Qx
            QP[, m] <- 0
          }
          kq <- min(kq, m - 1L)
          Qk <- QP[, seq_len(kq), drop = FALSE]
          v <- pt$df
          r1 <- as.numeric(crossprod(Qk, v)); v <- v - as.numeric(Qk %*% r1)
          r2 <- as.numeric(crossprod(Qk, v)); v <- v - as.numeric(Qk %*% r2)
          hv <- c(r1 + r2, .wr_l2(v))
          RP[1:(kq + 1L), kq + 1L] <- hv
          QP[, kq + 1L] <- v / hv[kq + 1L]
        }
        mk <- min(m, k)
        theta <- backsolve(RP[1:mk, 1:mk, drop = FALSE], as.numeric(crossprod(QP[, 1:mk, drop = FALSE], res)))
        sol <- gx - as.numeric(DG[, 1:mk, drop = FALSE] %*% theta)
      }
      k <- k + 1L
      pt <- point()
    }
    code <- if (resnorm <= tol) 0L else 10L
    if (!(resnorm < up)) code <- -2L
    if (toosoon) code <- -1L
    return(.wr_out(sol, Fv(sol), map_rc(code), k, -1L, -1L, -1L, -1L, nsteps = -1L))
  }
  fx <- Fv(x); resnorm <- .wr_l2(fx)
  tol <- rtol * resnorm + atol
  toosoon <- resnorm <= tol
  itc <- 0L
  krylov <- !is.null(a$linsolve) && a$method %in% c("newton", "pseudotransient")
  if (krylov) {
    # nsoli and ptcsoli: Jacobian free, J v by the forward difference dirder,
    # fixed forcing term 0.1, GMRES with a basis of 3 vectors (at most 2
    # iterations, as allocated by NonlinearSolve) or BiCGSTAB with 5
    lin <- function(x, fx, delta) .sf_krylov(function(v) .sf_dirder(Fv, x, fx, v, delta), fx, 0.1, a$linsolve)
    if (a$method == "pseudotransient") {
      # ptcsoli; NonlinearSolve does not pass `delta` here (ptcsoli then
      # starts from 1e-6), the port uses the algorithm's `delta`
      delta <- a$delta
      while (resnorm > tol && itc < maxit) {
        residm <- resnorm
        x <- x + (-lin(x, fx, delta)$sol)
        fx <- Fv(x); resnorm <- .wr_l2(fx)
        delta <- delta * (residm / resnorm)
        itc <- itc + 1L
      }
      code <- if (!(resnorm > tol || toosoon)) 0L else if (toosoon) -1L else 10L
      return(.wr_out(x, fx, map_rc(code), itc, -1L, -1L, -1L, -1L, nsteps = -1L))
    }
    ifun <- 1L; ijac <- 0L; iarm_tot <- 0L; armstop <- TRUE; iline <- FALSE
    while (resnorm > tol && itc < maxit && armstop) {
      kout <- lin(x, fx, 0)
      ao <- .sf_armijo(Fv, x, -kout$sol, resnorm, 10L)
      x <- ao$x; fx <- ao$f
      armstop <- ao$idid; iline <- !armstop
      resnorm <- ao$resnorm
      ifun <- ifun + ao$iarm + 1L
      ijac <- ijac + kout$lits * (if (a$linsolve == "gmres") 1L else 2L)
      iarm_tot <- iarm_tot + ao$iarm
      itc <- itc + 1L
    }
    code <- .sf_newton_err(resnorm, iline, tol, toosoon, itc, maxit)
    return(.wr_out(x, fx, map_rc(code), itc, ifun, ijac, 0L, 0L, nsteps = iarm_tot))
  }
  if (a$method == "pseudotransient") {
    dt <- a$delta
    while (resnorm > tol && itc < maxit) {
      if (!(dt > 0)) stop("dt must be > 0 in PTC", call. = FALSE)
      JJ <- jac(x, fx) + diag(1 / dt, n)
      x <- x + (-.wr_ldiv(JJ, fx))
      fx <- Fv(x)
      residm <- resnorm; resnorm <- .wr_l2(fx)
      dt <- dt * (residm / resnorm)
      itc <- itc + 1L
    }
    code <- if (!(resnorm > tol || toosoon)) 0L else if (toosoon) -1L else 10L
    return(.wr_out(x, fx, map_rc(code), itc, -1L, -1L, -1L, -1L, nsteps = -1L))
  }
  ifun <- if (a$method == "secant") 2L else 1L; ijac <- 0L; iarm_tot <- 0L
  armstop <- TRUE; iline <- FALSE
  newiarm <- -1L; residratio <- 1
  if (a$method == "secant") {
    if (n != 1L) stop("the secant method of SIAMFANLEquations is for scalar equations", call. = FALSE)
    xm <- x * 1.0001
    if (xm == 0) xm <- 1e-4
    fm <- Fv(xm)
    while (resnorm > tol && itc < maxit && armstop) {
      dfd <- (fx - fm) / (x - xm)
      d <- -fx / dfd
      ao <- .sf_armijo(Fv, x, d, resnorm, 5L)
      xm <- x; fm <- fx
      x <- ao$x; fx <- ao$f
      armstop <- ao$idid
      iline <- !armstop
      resnorm <- ao$resnorm
      ifun <- ifun + ao$iarm + 1L; iarm_tot <- iarm_tot + ao$iarm
      itc <- itc + 1L
    }
  } else {
    while (resnorm > tol && itc < maxit && armstop) {
      JJ <- jac(x, fx); ijac <- ijac + 1L
      step <- -.wr_ldiv(JJ, fx)
      ao <- .sf_armijo(Fv, x, step, resnorm, 10L)
      x <- ao$x; fx <- ao$f
      armstop <- ao$idid
      iline <- !armstop
      residm <- resnorm; resnorm <- ao$resnorm
      ifun <- ifun + ao$iarm + 1L; iarm_tot <- iarm_tot + ao$iarm
      itc <- itc + 1L
    }
  }
  code <- .sf_newton_err(resnorm, iline, tol, toosoon, itc, maxit)
  .wr_out(x, fx, map_rc(code), itc, ifun, ijac, 0L, 0L, nsteps = iarm_tot)
}

# ---------------------------------------------------------------------------
# PETSc SNES

#' @rdname nonlinear_wrappers
#' @export
PETScSNES <- function(petsclib = NULL, mpi_comm = NULL, autodiff = NULL, ...) {
  .wr_alg("PETScSNES", .snes_run, autodiff = autodiff, snes_options = list(...))
}

.snes_run <- function(a, F, J, x0, opts, ls) {
  o <- a$snes_options
  atol <- .wr_tol(opts$abstol); rtol <- .wr_tol(opts$reltol)
  stol <- o$snes_stol %||% 1e-8
  max_it <- as.integer(opts$maxiters %||% 1000L)
  max_funcs <- as.integer(o$snes_max_funcs %||% 10000L)
  order <- as.integer(o$snes_linesearch_order %||% 3L)
  jacf <- .wr_jac(F, J)
  nf <- 0L; njac <- 0L
  Fv <- function(x) { nf <<- nf + 1L; as.numeric(F(x)) }
  x <- as.numeric(x0)
  fx <- Fv(x); fnorm <- .wr_l2(fx); fnorm0 <- fnorm
  reason <- if (fnorm < atol) "CONVERGED_FNORM_ABS" else "ITERATING"
  it <- 0L
  while (reason == "ITERATING" && it < max_it) {
    it <- it + 1L
    Jx <- jacf(x); njac <- njac + 1L
    y <- tryCatch(.wr_ldiv(Jx, fx), error = function(e) NULL)
    if (is.null(y)) { reason <- "DIVERGED_LINEAR_SOLVE"; break }
    ynorm <- .wr_l2(y)
    # SNESLineSearchApply_BT: cubic (or quadratic) backtracking on ||F||^2 / 2
    alpha <- 1e-4; minlambda <- 1e-12; maxstep <- 1e8
    if (ynorm > maxstep) { y <- y * (maxstep / ynorm); ynorm <- maxstep }
    initslope <- -fnorm * fnorm
    lambda <- 1
    xn <- x - lambda * y; fn <- Fv(xn); gnorm <- .wr_l2(fn)
    if (0.5 * gnorm^2 > 0.5 * fnorm^2 + lambda * alpha * initslope) {
      lambdaprev <- lambda; gnormprev <- gnorm
      lambdatemp <- -initslope / (gnorm^2 - fnorm^2 - 2 * lambda * initslope)
      lambda <- if (lambdatemp > 0.5 * lambdaprev) 0.5 * lambdaprev else max(lambdatemp, 0.1 * lambdaprev)
      repeat {
        if (lambda <= minlambda) break
        xn <- x - lambda * y; fn <- Fv(xn); gnorm <- .wr_l2(fn)
        if (0.5 * gnorm^2 <= 0.5 * fnorm^2 + lambda * alpha * initslope) break
        if (order == 2L) {
          lambdatemp <- -initslope / (gnorm^2 - fnorm^2 - 2 * lambda * initslope)
        } else {
          t1 <- 0.5 * (gnorm^2 - fnorm^2) - lambda * initslope
          t2 <- 0.5 * (gnormprev^2 - fnorm^2) - lambdaprev * initslope
          aa <- (t1 / (lambda * lambda) - t2 / (lambdaprev * lambdaprev)) / (lambda - lambdaprev)
          bb <- (-lambdaprev * t1 / (lambda * lambda) + lambda * t2 / (lambdaprev * lambdaprev)) / (lambda - lambdaprev)
          d <- bb * bb - 3 * aa * initslope
          if (d < 0) d <- 0
          lambdatemp <- if (aa == 0) -initslope / (2 * bb) else (-bb + sqrt(d)) / (3 * aa)
        }
        lambdaprev <- lambda; gnormprev <- gnorm
        lambda <- if (lambdatemp > 0.5 * lambda) 0.5 * lambda else if (lambdatemp <= 0.1 * lambda) 0.1 * lambda else lambdatemp
      }
    }
    xold <- x
    x <- xn; fx <- fn; fnorm <- gnorm
    xnorm <- .wr_l2(x); snorm <- .wr_l2(x - xold)
    if (!is.finite(fnorm)) reason <- "DIVERGED_FNORM_NAN"
    else if (fnorm < atol) reason <- "CONVERGED_FNORM_ABS"
    else if (nf >= max_funcs) reason <- "DIVERGED_FUNCTION_COUNT"
    else if (fnorm <= rtol * fnorm0) reason <- "CONVERGED_FNORM_RELATIVE"
    else if (snorm < stol * xnorm) reason <- "CONVERGED_SNORM_RELATIVE"
  }
  if (reason == "ITERATING") reason <- "DIVERGED_MAX_IT"
  rc <- if (max(abs(fx)) <= atol) "Success" else "Failure"
  .wr_out(x, fx, rc, -1L, nf, njac, -1L, -1L, original = list(reason = reason, iterations = it))
}
