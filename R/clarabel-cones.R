# Cones of the Clarabel interior point solver: zero, nonnegative, second
# order, exponential, power, generalized power and positive semidefinite
# (triangle) cones, with their barriers, scalings, step lengths and
# corrector terms.
#
# Ported from Clarabel.jl 0.11.1 (https://github.com/oxfordcontrol/Clarabel.jl),
# Copyright 2022 University of Oxford Control Group, Apache License 2.0
# (inst/CLARABEL-LICENSE).

.clk_make <- function(spec) {
  K <- new.env(parent = emptyenv())
  switch(spec$type,
    ZeroConeT = {
      K$type <- "Zero"; K$dim <- spec$dim; K$numel <- spec$dim; K$degree <- 0
    },
    NonnegativeConeT = {
      K$type <- "Nonnegative"; K$dim <- spec$dim; K$numel <- spec$dim; K$degree <- spec$dim
      K$w <- numeric(spec$dim); K$lambda <- numeric(spec$dim)
    },
    SecondOrderConeT = {
      if (spec$dim < 2) stop("dimension must be >= 2", call. = FALSE)
      K$type <- "SecondOrder"; K$dim <- spec$dim; K$numel <- spec$dim; K$degree <- 1
      K$w <- numeric(spec$dim); K$lambda <- numeric(spec$dim); K$eta <- 0
      # cones above this size are expanded into sparse KKT columns
      K$expand <- spec$dim > 4
      if (K$expand) {
        K$sp_u <- numeric(spec$dim); K$sp_v <- numeric(spec$dim); K$sp_d <- 0
      }
    },
    ExponentialConeT = {
      K$type <- "Exponential"; K$numel <- 3L; K$degree <- 3
      K$H_dual <- matrix(0, 3, 3); K$Hs <- matrix(0, 3, 3); K$grad <- numeric(3); K$z <- numeric(3)
    },
    PowerConeT = {
      K$type <- "Power"; K$numel <- 3L; K$degree <- 3; K$alpha <- spec$alpha
      K$H_dual <- matrix(0, 3, 3); K$Hs <- matrix(0, 3, 3); K$grad <- numeric(3); K$z <- numeric(3)
    },
    GenPowerConeT = {
      a <- spec$alpha
      K$type <- "GenPower"; K$alpha <- a; K$dim1 <- length(a); K$dim2 <- spec$dim2
      K$numel <- length(a) + spec$dim2; K$degree <- length(a) + 1
      K$grad <- numeric(K$numel); K$z <- numeric(K$numel); K$mu <- 1
      K$p <- numeric(K$numel); K$q <- numeric(K$dim1); K$r <- numeric(K$dim2)
      K$d1 <- numeric(K$dim1); K$d2 <- 0; K$psi <- 1 / .jl_dot(a, a)
    },
    PSDTriangleConeT = {
      n <- spec$dim
      if (n < 0) stop("dimension must be non-negative", call. = FALSE)
      K$type <- "PSDTriangle"; K$n <- n; K$numel <- .cl_triangular_number(n); K$degree <- n
      K$lambda <- numeric(n); K$Lisqrt <- numeric(n)
      K$R <- matrix(0, n, n); K$Rinv <- matrix(0, n, n)
      K$Hs <- matrix(0, K$numel, K$numel)
    },
    stop("Unsupported cone type ", spec$type, call. = FALSE))
  K
}

.clk_is_symmetric <- function(K) K$type %in% c("Zero", "Nonnegative", "SecondOrder", "PSDTriangle")
.clk_allows_primal_dual_scaling <- function(K) K$type != "GenPower"
.clk_sparse_expandable <- function(K) (K$type == "SecondOrder" && K$expand) || K$type == "GenPower"
.clk_Hs_is_diagonal <- function(K) switch(K$type, SecondOrder = K$expand, Exponential = FALSE, Power = FALSE,
                                           PSDTriangle = FALSE, TRUE)

# ---------------------------------------------------------------------------
# Shared helpers

.clk_soc_residual <- function(z) {
  z1norm <- .cl_norm2(z[-1L])
  (z[1L] - z1norm) * (z[1L] + z1norm)
}
.clk_sqrt_soc_residual <- function(z) {
  res <- .clk_soc_residual(z)
  if (res > 0) sqrt(res) else 0
}
.clk_soc_residual_shifted <- function(z, dz, a) {
  x0 <- z[1L] + a * dz[1L]
  x1norm <- .cl_norm_shifted(z[-1L], dz[-1L], a)
  (x0 - x1norm) * (x0 + x1norm)
}
.clk_sumsq <- function(v) .cl_ssum(v * v)
.clk_triangular_index <- function(k) .cl_triangular_number(k)

.clk_svec_to_mat <- function(x, n) {
  M <- matrix(0, n, n)
  isqrt2 <- 1 / sqrt(2)
  idx <- 1L
  for (col in seq_len(n)) for (row in seq_len(col)) {
    if (row == col) M[row, col] <- x[idx] else {
      M[row, col] <- x[idx] * isqrt2
      M[col, row] <- x[idx] * isqrt2
    }
    idx <- idx + 1L
  }
  M
}
.clk_mat_to_svec <- function(M) {
  n <- nrow(M)
  isqrt2 <- 1 / sqrt(2)
  x <- numeric(.cl_triangular_number(n))
  idx <- 1L
  for (col in seq_len(n)) for (row in seq_len(col)) {
    x[idx] <- if (row == col) M[row, col] else (M[row, col] + M[col, row]) * isqrt2
    idx <- idx + 1L
  }
  x
}
.clk_pack_triu <- function(A) A[upper.tri(A, diag = TRUE)]

# Symmetric Kronecker product of A with itself in the svec basis
.clk_skron <- function(A) {
  n <- nrow(A)
  N <- .cl_triangular_number(n)
  out <- matrix(0, N, N)
  sqrt2 <- sqrt(2)
  col <- 1L
  for (l in seq_len(n)) for (k in seq_len(l)) {
    row <- 1L
    kl_eq <- k == l
    for (j in seq_len(n)) {
      Ajl <- A[j, l]; Ajk <- A[j, k]
      for (i in seq_len(j)) {
        if (row > col) break
        ij_eq <- i == j
        out[row, col] <- if (!ij_eq && !kl_eq) A[i, k] * Ajl + A[i, l] * Ajk
          else if (ij_eq && !kl_eq) sqrt2 * Ajl * Ajk
          else if (!ij_eq && kl_eq) sqrt2 * A[i, l] * Ajk
          else Ajl * Ajl
        row <- row + 1L
      }
    }
    col <- col + 1L
  }
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  out
}

.clk_chol_lower <- function(M) tryCatch(t(chol(M)), error = function(e) NULL)
.clk_eigmin <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

# W x for the PSD cone: R' X R (transpose FALSE) or R X R' (TRUE)
.clk_psd_mulW <- function(Rx, x, n, transpose) {
  X <- .clk_svec_to_mat(x, n)
  Y <- if (transpose) Rx %*% (X %*% t(Rx)) else (t(Rx) %*% X) %*% Rx
  .clk_mat_to_svec(Y)
}

# Second order cone helpers
.clk_soc_mulW <- function(K, x) {
  w <- K$w
  zeta <- .jl_dot(w[-1L], x[-1L])
  c0 <- x[1L] + zeta / (1 + w[1L])
  c(K$eta * (w[1L] * x[1L] + zeta), K$eta * (x[-1L] + c0 * w[-1L]))
}
.clk_soc_mulWinv <- function(K, x) {
  w <- K$w
  zeta <- .jl_dot(w[-1L], x[-1L])
  c0 <- -x[1L] + zeta / (1 + w[1L])
  etainv <- 1 / K$eta
  c(etainv * (w[1L] * x[1L] - zeta), etainv * (x[-1L] + c0 * w[-1L]))
}
.clk_soc_circ <- function(y, z) c(.jl_dot(y, z), y[1L] * z[-1L] + z[1L] * y[-1L])

.clk_step_length_soc_component <- function(x, y, amax) {
  if (x[1L] >= 0 && y[1L] < 0) amax <- min(amax, -x[1L] / y[1L])
  a <- .clk_soc_residual(y)
  b <- 2 * (x[1L] * y[1L] - .jl_dot(x[-1L], y[-1L]))
  c0 <- max(0, .clk_soc_residual(x))
  d <- b^2 - 4 * a * c0
  if ((a > 0 && b > 0) || d < 0) return(amax)
  if (a == 0) return(amax)
  if (c0 == 0) return(if (a >= 0) amax else 0)
  t0 <- if (b >= 0) -b - sqrt(d) else -b + sqrt(d)
  r1 <- (2 * c0) / t0
  r2 <- t0 / (2 * a)
  if (r1 < 0) r1 <- .cl_floatmax
  if (r2 < 0) r2 <- .cl_floatmax
  min(amax, r1, r2)
}

# Nonsymmetric cone helpers
.clk_backtrack <- function(dq, q, a_init, a_min, step, feasible) {
  a <- a_init
  repeat {
    if (feasible(q + a * dq)) break
    a <- a * step
    if (a < a_min) {
      a <- 0
      break
    }
  }
  a
}

.clk_newton_onesided <- function(x0, f0, f1) {
  iter <- 0L
  x <- x0
  while (iter < 100L) {
    iter <- iter + 1L
    dfdx <- f1(x)
    dx <- -f0(x) / dfdx
    if (dx < .cl_eps || abs(dx / x) < sqrt(.cl_eps) || abs(dfdx) < .cl_eps) break
    x <- x + dx
  }
  x
}

.clk_dot3 <- function(x, y) x[1L] * y[1L] + x[2L] * y[2L] + x[3L] * y[3L]
# dot of two StaticArrays: the @simd reduction allows contraction, so the
# product sum is a chain of fused multiply adds
.clk_sdot3 <- function(x, y) .jl_dot(x[1:3], y[1:3])

.clk_update_Hs <- function(K, s, z, mu, strategy) {
  if (strategy == "Dual") {
    K$Hs <- mu * K$H_dual
    return(invisible())
  }
  H_dual <- K$H_dual
  st <- K$grad
  zt <- .clk_gradient_primal(K, s)
  dot_sz <- .jl_dot(z, s)
  mu0 <- dot_sz / 3
  mut <- .clk_sdot3(zt, st) / 3
  ds <- s + mu0 * st
  dz <- z + mu0 * zt
  dot_dsz <- .clk_sdot3(ds, dz)
  de1 <- mu0 * mut - 1
  # zt' H_dual zt, column by column as LinearAlgebra.dot(x, A, y)
  de2 <- 0
  for (j in 1:3) if (zt[j] != 0) {
    tmp <- 0
    for (i in 1:3) tmp <- tmp + H_dual[i, j] * zt[i]
    de2 <- de2 + tmp * zt[j]
  }
  de2 <- de2 - 3 * mut * mut
  if (abs(de1) > sqrt(.cl_eps) && abs(de2) > .cl_eps && dot_sz > 0 && dot_dsz > 0) {
    tmp <- vapply(1:3, function(i) mut * st[i] - H_dual[i, 1L] * zt[1L] - H_dual[i, 2L] * zt[2L] -
                    H_dual[i, 3L] * zt[3L], 0)
    Hs <- H_dual
    for (i in 1:3) for (j in 1:3) Hs[i, j] <- Hs[i, j] - (st[i] * st[j] / 3 + tmp[i] * tmp[j] / de2)
    t0 <- mu0 * .cl_norm2(as.vector(Hs))
    axis_z <- c(z[2L] * zt[3L] - z[3L] * zt[2L], z[3L] * zt[1L] - z[1L] * zt[3L], z[1L] * zt[2L] - z[2L] * zt[1L])
    axis_z <- axis_z * (1 / .cl_norm2(axis_z))
    for (i in 1:3) for (j in i:3)
      Hs[i, j] <- s[i] * s[j] / dot_sz + ds[i] * ds[j] / dot_dsz + t0 * axis_z[i] * axis_z[j]
    Hs[2L, 1L] <- Hs[1L, 2L]; Hs[3L, 1L] <- Hs[1L, 3L]; Hs[3L, 2L] <- Hs[2L, 3L]
    K$Hs <- Hs
  } else {
    # the fallback scales by the local mu = <s, z>/3 of this cone, not by the
    # mu of the whole iterate
    K$Hs <- mu0 * K$H_dual
  }
  invisible()
}

.clk_chol3_factor <- function(A) {
  L <- matrix(0, 3, 3)
  t0 <- A[1, 1]
  if (t0 <= 0) return(NULL)
  L[1, 1] <- sqrt(t0)
  L[2, 1] <- A[2, 1] / L[1, 1]
  t0 <- A[2, 2] - L[2, 1] * L[2, 1]
  if (t0 <= 0) return(NULL)
  L[2, 2] <- sqrt(t0)
  L[3, 1] <- A[3, 1] / L[1, 1]
  L[3, 2] <- (A[3, 2] - L[2, 1] * L[3, 1]) / L[2, 2]
  t0 <- A[3, 3] - L[3, 1] * L[3, 1] - L[3, 2] * L[3, 2]
  if (t0 <= 0) return(NULL)
  L[3, 3] <- sqrt(t0)
  L
}
.clk_chol3_solve <- function(L, b) {
  c1 <- b[1] / L[1, 1]
  c2 <- (b[2] * L[1, 1] - b[1] * L[2, 1]) / (L[1, 1] * L[2, 2])
  c3 <- (b[3] * L[1, 1] * L[2, 2] - b[2] * L[1, 1] * L[3, 2] + b[1] * L[2, 1] * L[3, 2] - b[1] * L[2, 2] * L[3, 1]) /
    (L[1, 1] * L[2, 2] * L[3, 3])
  x1 <- (c1 * L[2, 2] * L[3, 3] - c2 * L[2, 1] * L[3, 3] + c3 * L[2, 1] * L[3, 2] - c3 * L[2, 2] * L[3, 1]) /
    (L[1, 1] * L[2, 2] * L[3, 3])
  x2 <- (c2 * L[3, 3] - c3 * L[3, 2]) / (L[2, 2] * L[3, 3])
  x3 <- c3 / L[3, 3]
  c(x1, x2, x3)
}

# ---------------------------------------------------------------------------
# Exponential cone

.clk_wright_omega <- function(z) {
  if (z < 0) stop("argument not in supported range : ", z, call. = FALSE)
  if (z < 1 + pi) {
    zm1 <- z - 1
    p <- zm1
    w <- 1 + 0.5 * p
    p <- p * zm1
    w <- w + (1 / 16.0) * p
    p <- p * zm1
    w <- w - (1 / 192.0) * p
    p <- p * zm1
    w <- w - (1 / 3072.0) * p
    p <- p * zm1
    w <- w + (13 / 61440.0) * p
  } else {
    logz <- .cl_logsafe(z)
    zinv <- 1 / z
    w <- z - logz
    q <- logz * zinv
    w <- w + q
    q <- q * zinv
    w <- w + q * (logz / 2 - 1)
    # Clarabel.jl leaves q unchanged here (the statement `q * zinv` is not an
    # update); kept as upstream
    w <- w + q * (logz * logz / 3. - (3 / 2.) * logz + 1)
  }
  r <- z - w - .cl_logsafe(w)
  for (i in 1:2) {
    wp1 <- w + 1
    t0 <- wp1 * (wp1 + (2. * r) / 3.0)
    w <- w * (1 + (r / wp1) * (t0 - 0.5 * r) / (t0 - r))
    r <- (2 * w * w - 8 * w - 1) / (72.0 * (wp1 * wp1 * wp1 * wp1 * wp1 * wp1)) * r * r * r * r
  }
  w
}

.clk_exp_barrier_dual <- function(z) {
  l <- .cl_logsafe(-z[3L] / z[1L])
  -.cl_logsafe(-z[3L] * z[1L]) - .cl_logsafe(z[2L] - z[1L] - z[1L] * l)
}
.clk_exp_barrier_primal <- function(s) {
  om <- .clk_wright_omega(1 - s[1L] / s[2L] - .cl_logsafe(s[2L] / s[3L]))
  om <- (om - 1) * (om - 1) / om
  -.cl_logsafe(om) - 2 * .cl_logsafe(s[2L]) - .cl_logsafe(s[3L]) - 3
}
.clk_exp_primal_feasible <- function(s) {
  if (s[3L] > 0 && s[2L] > 0) {
    res <- s[2L] * .cl_logsafe(s[3L] / s[2L]) - s[1L]
    if (res > 0) return(TRUE)
  }
  FALSE
}
.clk_exp_dual_feasible <- function(z) {
  if (z[3L] > 0 && z[1L] < 0) {
    res <- z[2L] - z[1L] - z[1L] * .cl_logsafe(-z[3L] / z[1L])
    if (res > 0) return(TRUE)
  }
  FALSE
}
.clk_exp_gradient_primal <- function(s) {
  om <- .clk_wright_omega(1 - s[1L] / s[2L] - .cl_logsafe(s[2L] / s[3L]))
  g1 <- 1 / ((om - 1) * s[2L])
  g2 <- g1 + g1 * .cl_logsafe(om * s[2L] / s[3L]) - 1 / s[2L]
  g3 <- om / ((1 - om) * s[3L])
  c(g1, g2, g3)
}
.clk_exp_update_dual_grad_H <- function(K, z) {
  l <- .cl_logsafe(-z[3L] / z[1L])
  r <- -z[1L] * l - z[1L] + z[2L]
  c2 <- 1 / r
  K$grad <- c(c2 * l - 1 / z[1L], -c2, (c2 * z[1L] - 1) / z[3L])
  H <- matrix(0, 3, 3)
  H[1, 1] <- ((r * r - z[1] * r + l * l * z[1] * z[1]) / (r * z[1] * z[1] * r))
  H[1, 2] <- (-l / (r * r)); H[2, 1] <- H[1, 2]
  H[2, 2] <- (1 / (r * r))
  H[1, 3] <- ((z[2] - z[1]) / (r * r * z[3])); H[3, 1] <- H[1, 3]
  H[2, 3] <- (-z[1] / (r * r * z[3])); H[3, 2] <- H[2, 3]
  H[3, 3] <- ((r * r - z[1] * r + z[1] * z[1]) / (r * r * z[3] * z[3]))
  K$H_dual <- H
}
.clk_exp_higher_correction <- function(K, ds, v) {
  cholH <- .clk_chol3_factor(K$H_dual)
  if (is.null(cholH)) return(numeric(3))
  u <- .clk_chol3_solve(cholH, ds)
  z <- K$z
  eta <- numeric(3)
  eta[2] <- 1
  eta[3] <- -z[1] / z[3]
  eta[1] <- .cl_logsafe(eta[3])
  psi <- z[1] * eta[1] - z[1] + z[2]
  dotpsiu <- .clk_sdot3(eta, u)
  dotpsiv <- .clk_dot3(eta, v)
  coef <- ((u[1] * (v[1] / z[1] - v[3] / z[3]) + u[3] * (z[1] * v[3] / z[3] - v[1]) / z[3]) * psi -
             2 * dotpsiu * dotpsiv) / (psi * psi * psi)
  eta <- eta * coef
  inv_psi2 <- 1 / psi / psi
  # Julia's x += a - b + c adds the whole sum to x
  eta[1] <- eta[1] + ((1 / psi - 2 / z[1]) * u[1] * v[1] / (z[1] * z[1]) - u[3] * v[3] / (z[3] * z[3]) / psi +
    dotpsiu * inv_psi2 * (v[1] / z[1] - v[3] / z[3]) + dotpsiv * inv_psi2 * (u[1] / z[1] - u[3] / z[3]))
  eta[3] <- eta[3] + (2 * (z[1] / psi - 1) * u[3] * v[3] / (z[3] * z[3] * z[3]) -
    (u[3] * v[1] + u[1] * v[3]) / (z[3] * z[3]) / psi +
    dotpsiu * inv_psi2 * (z[1] * v[3] / (z[3] * z[3]) - v[1] / z[3]) +
    dotpsiv * inv_psi2 * (z[1] * u[3] / (z[3] * z[3]) - u[1] / z[3]))
  eta / 2
}

# ---------------------------------------------------------------------------
# Power cone

.clk_pow_primal_feasible <- function(s, a) {
  if (s[1L] > 0 && s[2L] > 0) {
    res <- exp(2 * a * .cl_logsafe(s[1L]) + 2 * (1 - a) * .cl_logsafe(s[2L])) - s[3L] * s[3L]
    if (res > 0) return(TRUE)
  }
  FALSE
}
.clk_pow_dual_feasible <- function(z, a) {
  if (z[1L] > 0 && z[2L] > 0) {
    res <- exp(2 * a * .cl_logsafe(z[1L] / a) + 2 * (1 - a) * .cl_logsafe(z[2L] / (1 - a))) - z[3L] * z[3L]
    if (res > 0) return(TRUE)
  }
  FALSE
}
.clk_newton_powcone <- function(s3, phi, a) {
  x0 <- -1 / s3 + (2 * s3 + sqrt(phi * phi / s3 / s3 + 3 * phi)) / (phi - s3 * s3)
  t0 <- -2 * a * .cl_logsafe(a) - 2 * (1 - a) * .cl_logsafe(1 - a)
  f0 <- function(x) {
    t1 <- x * x; t2 <- 2 * x / s3
    2 * a * .cl_logsafe(2 * a * t1 + (1 + a) * t2) + 2 * (1 - a) * .cl_logsafe(2 * (1 - a) * t1 + (2 - a) * t2) -
      .cl_logsafe(phi) - .cl_logsafe(t1 + t2) - 2 * .cl_logsafe(t2) + t0
  }
  f1 <- function(x) {
    t1 <- x * x; t2 <- x * 2 / s3
    2 * a * a / (a * x + (1 + a) / s3) + 2 * (1 - a) * (1 - a) / ((1 - a) * x + (2 - a) / s3) -
      2 * (x + 1 / s3) / (t1 + t2)
  }
  .clk_newton_onesided(x0, f0, f1)
}
.clk_pow_gradient_primal <- function(K, s) {
  a <- K$alpha
  phi <- (s[1L])^(2 * a) * (s[2L])^(2 - 2 * a)
  g <- numeric(3)
  abs_s <- abs(s[3L])
  if (abs_s > .cl_eps) {
    g[3] <- .clk_newton_powcone(abs_s, phi, a)
    if (s[3L] < 0) g[3] <- -g[3]
    g[1] <- -(a * g[3] * s[3L] + 1 + a) / s[1L]
    g[2] <- -((1 - a) * g[3] * s[3L] + 2 - a) / s[2L]
  } else {
    g[3] <- 0
    g[1] <- -(1 + a) / s[1L]
    g[2] <- -(2 - a) / s[2L]
  }
  g
}
.clk_pow_barrier_dual <- function(K, z) {
  a <- K$alpha
  -.cl_logsafe((z[1] / a)^(2 * a) * (z[2] / (1 - a))^(2 - 2 * a) - z[3] * z[3]) - (1 - a) * .cl_logsafe(z[1]) -
    a * .cl_logsafe(z[2])
}
.clk_pow_barrier_primal <- function(K, s) {
  a <- K$alpha
  g <- .clk_pow_gradient_primal(K, s)
  .cl_logsafe((-g[1] / a)^(2 * a) * (-g[2] / (1 - a))^(2 - 2 * a) - g[3] * g[3]) + (1 - a) * .cl_logsafe(-g[1]) +
    a * .cl_logsafe(-g[2]) - 3
}
.clk_pow_update_dual_grad_H <- function(K, z) {
  a <- K$alpha
  phi <- (z[1] / a)^(2 * a) * (z[2] / (1 - a))^(2 - 2 * a)
  psi <- phi - z[3] * z[3]
  g <- c(2 * a * phi / (z[1] * psi), 2 * (1 - a) * phi / (z[2] * psi), -2 * z[3] / psi)
  H <- matrix(0, 3, 3)
  H[1, 1] <- g[1] * g[1] - 2 * a * (2 * a - 1) * phi / (z[1] * z[1] * psi) + (1 - a) / (z[1] * z[1])
  H[1, 2] <- g[1] * g[2] - 4 * a * (1 - a) * phi / (z[1] * z[2] * psi); H[2, 1] <- H[1, 2]
  H[2, 2] <- g[2] * g[2] - 2 * (1 - a) * (1 - 2 * a) * phi / (z[2] * z[2] * psi) + a / (z[2] * z[2])
  H[1, 3] <- g[1] * g[3]; H[3, 1] <- H[1, 3]
  H[2, 3] <- g[2] * g[3]; H[3, 2] <- H[2, 3]
  H[3, 3] <- g[3] * g[3] + 2 / psi
  K$H_dual <- H
  K$grad <- c(-2 * a * phi / (z[1] * psi) - (1 - a) / z[1], -2 * (1 - a) * phi / (z[2] * psi) - a / z[2],
              2 * z[3] / psi)
}
.clk_pow_higher_correction <- function(K, ds, v) {
  cholH <- .clk_chol3_factor(K$H_dual)
  if (is.null(cholH)) return(numeric(3))
  u <- .clk_chol3_solve(cholH, ds)
  a <- K$alpha
  z <- K$z
  phi <- (z[1] / a)^(2 * a) * (z[2] / (1 - a))^(2 - 2 * a)
  psi <- phi - z[3] * z[3]
  eta <- c(2 * a * phi / z[1], 2 * (1 - a) * phi / z[2], -2 * z[3])
  H11 <- 2 * a * (2 * a - 1) * phi / (z[1] * z[1])
  H12 <- 4 * a * (1 - a) * phi / (z[1] * z[2])
  H22 <- 2 * (1 - a) * (1 - 2 * a) * phi / (z[2] * z[2])
  dotpsiu <- .clk_sdot3(eta, u)
  dotpsiv <- .clk_dot3(eta, v)
  Hv <- c(H11 * v[1] + H12 * v[2], H12 * v[1] + H22 * v[2], -2 * v[3])
  coef <- (.clk_sdot3(u, Hv) * psi - 2 * dotpsiu * dotpsiv) / (psi * psi * psi)
  coef2 <- 4 * a * (2 * a - 1) * (1 - a) * phi * (u[1] / z[1] - u[2] / z[2]) * (v[1] / z[1] - v[2] / z[2]) / psi
  inv_psi2 <- 1 / psi / psi
  eta[1] <- coef * eta[1] - 2 * (1 - a) * u[1] * v[1] / (z[1] * z[1] * z[1]) + coef2 / z[1] +
    Hv[1] * dotpsiu * inv_psi2
  eta[2] <- coef * eta[2] - 2 * a * u[2] * v[2] / (z[2] * z[2] * z[2]) - coef2 / z[2] + Hv[2] * dotpsiu * inv_psi2
  eta[3] <- coef * eta[3] + Hv[3] * dotpsiu * inv_psi2
  Hu <- c(H11 * u[1] + H12 * u[2], H12 * u[1] + H22 * u[2], -2 * u[3])
  (eta + Hu * dotpsiv * inv_psi2) / 2
}

.clk_gradient_primal <- function(K, s) if (K$type == "Exponential") .clk_exp_gradient_primal(s) else
  .clk_pow_gradient_primal(K, s)

# ---------------------------------------------------------------------------
# Generalized power cone

.clk_gp_logprod <- function(x, a, scale) {
  res <- 0
  for (i in seq_along(a)) res <- res + 2 * a[i] * .cl_logsafe(if (scale) x[i] / a[i] else x[i])
  res
}
.clk_gp_primal_feasible <- function(K, s) {
  d1 <- K$dim1
  if (all(s[seq_len(d1)] > 0)) {
    res <- exp(.clk_gp_logprod(s, K$alpha, FALSE)) - .clk_sumsq(s[.cl_rng(d1 + 1L, K$numel)])
    if (res > 0) return(TRUE)
  }
  FALSE
}
.clk_gp_dual_feasible <- function(K, z) {
  d1 <- K$dim1
  if (all(z[seq_len(d1)] > 0)) {
    res <- exp(.clk_gp_logprod(z, K$alpha, TRUE)) - .clk_sumsq(z[.cl_rng(d1 + 1L, K$numel)])
    if (res > 0) return(TRUE)
  }
  FALSE
}
.clk_gp_barrier_dual <- function(K, z) {
  a <- K$alpha
  res <- exp(.clk_gp_logprod(z, a, TRUE)) - .clk_sumsq(z[.cl_rng(K$dim1 + 1L, K$numel)])
  barrier <- -.cl_logsafe(res)
  for (i in seq_len(K$dim1)) barrier <- barrier - (1 - a[i]) * .cl_logsafe(z[i])
  barrier
}
.clk_gp_gradient_primal <- function(K, s) {
  a <- K$alpha
  d1 <- K$dim1
  phi <- 1
  for (i in seq_len(d1)) phi <- phi * s[i]^(2 * a[i])
  p <- s[seq_len(d1)]
  r <- s[.cl_rng(d1 + 1L, K$numel)]
  norm_r <- .cl_norm2(r)
  if (norm_r > .cl_eps) {
    g1 <- .clk_newton_genpowcone(norm_r, p, phi, a, K$psi)
    c(-(1 + a + a * g1 * norm_r) / p, g1 * r / norm_r)
  } else c(-(1 + a) / p, numeric(length(r)))
}
.clk_newton_genpowcone <- function(norm_r, p, phi, a, psi) {
  x0 <- -1 / norm_r + (psi * norm_r + sqrt((phi / norm_r / norm_r + psi * psi - 1) * phi)) / (phi - norm_r * norm_r)
  f0 <- function(x) {
    out <- -.cl_logsafe(2 * x / norm_r + x * x)
    for (i in seq_along(a)) out <- out + 2 * a[i] * (.cl_logsafe(x * norm_r + (1 + a[i]) / a[i]) - .cl_logsafe(p[i]))
    out
  }
  f1 <- function(x) {
    out <- -(2 * x + 2 / norm_r) / (x * x + 2 * x / norm_r)
    for (i in seq_along(a)) out <- out + 2 * a[i] * norm_r / (norm_r * x + (1 + a[i]) / a[i])
    out
  }
  .clk_newton_onesided(x0, f0, f1)
}
.clk_gp_barrier_primal <- function(K, s) -.clk_gp_barrier_dual(K, -.clk_gp_gradient_primal(K, s)) - K$degree
.clk_gp_update_dual_grad_H <- function(K, z) {
  a <- K$alpha
  d1 <- K$dim1
  r2 <- .cl_rng(d1 + 1L, K$numel)
  phi <- 1
  for (i in seq_len(d1)) phi <- phi * (z[i] / a[i])^(2 * a[i])
  norm2w <- .clk_sumsq(z[r2])
  zeta <- phi - norm2w
  if (!(zeta > 0)) stop("AssertionError: zeta > 0", call. = FALSE)
  tau <- 2 * a / z[seq_len(d1)]
  grad <- numeric(K$numel)
  grad[seq_len(d1)] <- -tau * phi / zeta - (1 - a) / z[seq_len(d1)]
  grad[r2] <- 2 * z[r2] / zeta
  K$grad <- grad
  p0 <- sqrt(phi * (phi + norm2w) / 2)
  p1 <- -2 * phi / p0
  q0 <- sqrt(zeta * phi / 2)
  r1 <- 2 * sqrt(zeta / (phi + norm2w))
  K$d1 <- tau * phi / (zeta * z[seq_len(d1)]) + (1 - a) / (z[seq_len(d1)] * z[seq_len(d1)])
  K$d2 <- 2 / zeta
  K$p <- c(p0 * tau / zeta, p1 * z[r2] / zeta)
  K$q <- tau * (q0 / zeta)
  K$r <- r1 * z[r2] / zeta
}

# ---------------------------------------------------------------------------
# Cone operations dispatched on the cone type

.clk_margins <- function(K, z, pd) switch(K$type,
  Zero = c(.cl_floatmax, 0),
  Nonnegative = c(if (length(z)) min(z) else .cl_floatmax, .cl_ssum(ifelse(z > 0, z, 0))),
  SecondOrder = {
    a <- z[1L] - .cl_norm2(z[-1L])
    c(a, max(0, a))
  },
  PSDTriangle = {
    if (!length(z)) c(.cl_floatmax, 0) else {
      e <- eigen(.clk_svec_to_mat(z, K$n), symmetric = TRUE, only.values = TRUE)$values
      c(min(e), .cl_ssum(ifelse(e > 0, e, 0)))
    }
  },
  stop("This function should never be reached.", call. = FALSE))

.clk_scaled_unit_shift <- function(K, z, a, pd) switch(K$type,
  Zero = if (pd == "primal") z * 0 else z,
  Nonnegative = z + a,
  SecondOrder = {
    z[1L] <- z[1L] + a
    z
  },
  PSDTriangle = {
    for (k in seq_len(K$n)) z[.clk_triangular_index(k)] <- z[.clk_triangular_index(k)] + a
    z
  },
  stop("This function should never be reached.", call. = FALSE))

.clk_unit_initialization <- function(K) switch(K$type,
  Zero = list(z = numeric(K$numel), s = numeric(K$numel)),
  Nonnegative = list(z = rep(1, K$numel), s = rep(1, K$numel)),
  SecondOrder = {
    e <- c(1, numeric(K$numel - 1L))
    list(z = e, s = e)
  },
  Exponential = {
    s <- c(-1.051383945322714, 0.556409619469370, 1.258967884768947)
    list(z = s, s = s)
  },
  Power = {
    s <- c(sqrt(1 + K$alpha), sqrt(1 + (1 - K$alpha)), 0)
    list(z = s, s = s)
  },
  GenPower = {
    s <- c(sqrt(1 + K$alpha), numeric(K$dim2))
    list(z = s, s = s)
  },
  PSDTriangle = {
    e <- .clk_scaled_unit_shift(K, numeric(K$numel), 1, "primal")
    list(z = e, s = e)
  })

.clk_set_identity_scaling <- function(K) {
  switch(K$type,
    Nonnegative = K$w <- rep(1, K$numel),
    SecondOrder = {
      K$w <- c(1, numeric(K$dim - 1L))
      K$eta <- 1
      if (K$expand) {
        K$sp_d <- 0.5
        K$sp_u <- c(sqrt(0.5), numeric(K$dim - 1L))
        K$sp_v <- numeric(K$dim)
      }
    },
    PSDTriangle = {
      K$R <- diag(1, K$n); K$Rinv <- diag(1, K$n); K$Hs <- diag(1, K$numel)
    },
    Zero = NULL,
    stop("This function should never be reached.", call. = FALSE))
  invisible()
}

.clk_update_scaling <- function(K, s, z, mu, strategy) {
  switch(K$type,
    Zero = TRUE,
    Nonnegative = {
      K$lambda <- sqrt(s * z)
      K$w <- sqrt(s / z)
      TRUE
    },
    SecondOrder = {
      zscale <- .clk_sqrt_soc_residual(z)
      sscale <- .clk_sqrt_soc_residual(s)
      if (zscale == 0 || sscale == 0) return(FALSE)
      K$eta <- sqrt(sscale / zscale)
      w <- s / sscale
      w[1L] <- w[1L] + z[1L] / zscale
      w[-1L] <- w[-1L] - z[-1L] / zscale
      wscale <- .clk_sqrt_soc_residual(w)
      if (wscale == 0) return(FALSE)
      w <- w / wscale
      w1sq <- .clk_sumsq(w[-1L])
      w[1L] <- sqrt(1 + w1sq)
      K$w <- w
      g <- 0.5 * wscale
      lam <- numeric(K$dim)
      lam[1L] <- g
      lam[-1L] <- ((g + z[1L] / zscale) / sscale) * s[-1L] + ((g + s[1L] / sscale) / zscale) * z[-1L]
      lam[-1L] <- lam[-1L] * (1 / (s[1L] / sscale + z[1L] / zscale + 2 * g))
      K$lambda <- lam * sqrt(sscale * zscale)
      if (K$expand) {
        a2 <- 2 * w[1L]
        wsq <- w[1L] * w[1L] + w1sq
        wsqinv <- 1 / wsq
        K$sp_d <- wsqinv / 2
        u0 <- sqrt(wsq - K$sp_d)
        u1 <- a2 / u0
        v1 <- sqrt(2 * (2 + wsqinv) / (2 * wsq - wsqinv))
        K$sp_u <- c(u0, u1 * w[-1L])
        K$sp_v <- c(0, v1 * w[-1L])
      }
      TRUE
    },
    Exponential = , Power = {
      if (K$type == "Exponential") .clk_exp_update_dual_grad_H(K, z) else .clk_pow_update_dual_grad_H(K, z)
      .clk_update_Hs(K, s, z, mu, strategy)
      K$z <- z
      TRUE
    },
    GenPower = {
      .clk_gp_update_dual_grad_H(K, z)
      K$mu <- mu
      K$z <- z
      TRUE
    },
    PSDTriangle = {
      if (!length(s)) return(TRUE)
      L1 <- .clk_chol_lower(.clk_svec_to_mat(s, K$n))
      L2 <- .clk_chol_lower(.clk_svec_to_mat(z, K$n))
      if (is.null(L1) || is.null(L2)) return(FALSE)
      sv <- svd(t(L2) %*% L1)
      K$lambda <- sv$d
      K$Lisqrt <- 1 / sqrt(sv$d)
      K$R <- (L1 %*% sv$v) %*% diag(K$Lisqrt, K$n)
      K$Rinv <- diag(K$Lisqrt, K$n) %*% (t(sv$u) %*% t(L2))
      K$Hs <- .clk_skron(K$R %*% t(K$R))
      TRUE
    })
}

.clk_get_Hs <- function(K) switch(K$type,
  Zero = numeric(K$numel),
  Nonnegative = K$w^2,
  SecondOrder = {
    if (K$expand) {
      H <- rep(K$eta^2, K$dim)
      H[1L] <- H[1L] * K$sp_d
      H
    } else {
      w <- K$w
      H <- numeric(.cl_triangular_number(K$dim))
      H[1L] <- (sqrt(2) * w[1L] - 1) * (sqrt(2) * w[1L] + 1)
      hidx <- 2L
      for (col in 2:K$dim) {
        for (row in seq_len(col)) {
          H[hidx] <- 2 * w[row] * w[col]
          hidx <- hidx + 1L
        }
        H[hidx - 1L] <- H[hidx - 1L] + 1
      }
      H * K$eta^2
    }
  },
  Exponential = , Power = .clk_pack_triu(K$Hs),
  GenPower = c(K$mu * K$d1, rep(K$mu * K$d2, K$dim2)),
  PSDTriangle = .clk_pack_triu(K$Hs))

.clk_mul_Hs <- function(K, x) switch(K$type,
  Zero = numeric(length(x)),
  Nonnegative = K$w * (K$w * x),
  SecondOrder = {
    c0 <- 2 * .jl_dot(K$w, x)
    y <- x
    y[1L] <- -x[1L]
    y <- y + c0 * K$w
    y * K$eta^2
  },
  Exponential = , Power = {
    H <- K$Hs
    vapply(1:3, function(i) H[i, 1L] * x[1L] + H[i, 2L] * x[2L] + H[i, 3L] * x[3L], 0)
  },
  GenPower = {
    d1 <- K$dim1
    r1 <- seq_len(d1); r2 <- .cl_rng(d1 + 1L, K$numel)
    coef_p <- .jl_dot(K$p, x)
    coef_q <- .jl_dot(K$q, x[r1])
    coef_r <- .jl_dot(K$r, x[r2])
    y <- c(K$d1 * x[r1] - coef_q * K$q, K$d2 * x[r2] - coef_r * K$r)
    y <- y + coef_p * K$p
    y * K$mu
  },
  PSDTriangle = .clk_psd_mulW(K$R, .clk_psd_mulW(K$R, x, K$n, FALSE), K$n, TRUE))

.clk_affine_ds <- function(K, s) switch(K$type,
  Zero = numeric(length(s)),
  Nonnegative = K$lambda^2,
  SecondOrder = .clk_soc_circ(K$lambda, K$lambda),
  PSDTriangle = {
    ds <- numeric(K$numel)
    for (k in seq_len(K$n)) ds[.clk_triangular_index(k)] <- K$lambda[k]^2
    ds
  },
  s)

.clk_combined_ds_shift <- function(K, step_z, step_s, sigmamu) switch(K$type,
  Zero = numeric(length(step_z)),
  Nonnegative = , SecondOrder = , PSDTriangle = {
    if (K$type == "Nonnegative") {
      wz <- step_z * K$w
      ws <- step_s / K$w
      shift <- ws * wz
    } else if (K$type == "SecondOrder") {
      wz <- .clk_soc_mulW(K, step_z)
      ws <- .clk_soc_mulWinv(K, step_s)
      shift <- .clk_soc_circ(ws, wz)
    } else {
      wz <- .clk_psd_mulW(K$R, step_z, K$n, FALSE)
      ws <- .clk_psd_mulW(K$Rinv, step_s, K$n, TRUE)
      Y <- .clk_svec_to_mat(ws, K$n); Z <- .clk_svec_to_mat(wz, K$n)
      shift <- .clk_mat_to_svec(0.5 * (Y %*% Z + Z %*% Y))
    }
    .clk_scaled_unit_shift(K, shift, -sigmamu, "primal")
  },
  Exponential = K$grad * sigmamu - .clk_exp_higher_correction(K, step_s, step_z),
  Power = K$grad * sigmamu - .clk_pow_higher_correction(K, step_s, step_z),
  GenPower = K$grad * sigmamu)

.clk_ds_from_dz_offset <- function(K, ds, z) switch(K$type,
  Zero = numeric(length(ds)),
  Nonnegative = ds / z,
  SecondOrder = {
    resz <- .clk_soc_residual(z)
    lam <- K$lambda; w <- K$w
    l1ds1 <- .jl_dot(lam[-1L], ds[-1L])
    w1ds1 <- .jl_dot(w[-1L], ds[-1L])
    out <- -z
    out[1L] <- z[1L]
    c0 <- lam[1L] * ds[1L] - l1ds1
    out <- out * (c0 / resz)
    out[1L] <- out[1L] + K$eta * w1ds1
    out[-1L] <- out[-1L] + K$eta * (ds[-1L] + w1ds1 / (1 + w[1L]) * w[-1L])
    out * (1 / lam[1L])
  },
  PSDTriangle = {
    Z <- .clk_svec_to_mat(ds, K$n)
    lam <- K$lambda
    X <- Z
    for (i in seq_len(K$n)) for (j in seq_len(K$n)) X[i, j] <- 2 * Z[i, j] / (lam[i] + lam[j])
    .clk_psd_mulW(K$R, .clk_mat_to_svec(X), K$n, TRUE)
  },
  ds)

.clk_step_length <- function(K, dz, ds, z, s, settings, amax) {
  step <- settings$linesearch_backtrack_step
  amin <- settings$min_terminate_step_length
  switch(K$type,
    Zero = c(amax, amax),
    Nonnegative = {
      az <- amax; as <- amax
      for (i in seq_along(ds)) {
        if (dz[i] < 0) az <- min(az, -z[i] / dz[i])
        if (ds[i] < 0) as <- min(as, -s[i] / ds[i])
      }
      c(az, as)
    },
    SecondOrder = c(.clk_step_length_soc_component(z, dz, amax), .clk_step_length_soc_component(s, ds, amax)),
    Exponential = c(.clk_backtrack(dz, z, amax, amin, step, .clk_exp_dual_feasible),
                    .clk_backtrack(ds, s, amax, amin, step, .clk_exp_primal_feasible)),
    Power = c(.clk_backtrack(dz, z, amax, amin, step, function(v) .clk_pow_dual_feasible(v, K$alpha)),
              .clk_backtrack(ds, s, amax, amin, step, function(v) .clk_pow_primal_feasible(v, K$alpha))),
    GenPower = c(.clk_backtrack(dz, z, amax, amin, step, function(v) .clk_gp_dual_feasible(K, v)),
                 .clk_backtrack(ds, s, amax, amin, step, function(v) .clk_gp_primal_feasible(K, v))),
    PSDTriangle = {
      comp <- function(d) {
        if (!length(d)) return(amax)
        D <- .clk_svec_to_mat(d, K$n)
        D <- D * outer(K$Lisqrt, K$Lisqrt)
        g <- .clk_eigmin(D)
        if (g < 0) min(1 / -g, amax) else amax
      }
      c(comp(.clk_psd_mulW(K$R, dz, K$n, FALSE)), comp(.clk_psd_mulW(K$Rinv, ds, K$n, TRUE)))
    })
}

.clk_compute_barrier <- function(K, z, s, dz, ds, a) switch(K$type,
  Zero = 0,
  Nonnegative = {
    barrier <- 0
    for (i in seq_len(K$dim)) barrier <- barrier - .cl_logsafe((s[i] + a * ds[i]) * (z[i] + a * dz[i]))
    barrier
  },
  SecondOrder = {
    res_s <- .clk_soc_residual_shifted(s, ds, a)
    res_z <- .clk_soc_residual_shifted(z, dz, a)
    if (res_s > 0 && res_z > 0) -.cl_logsafe(res_s * res_z) / 2 else Inf
  },
  Exponential = 0 + .clk_exp_barrier_dual(z + a * dz) + .clk_exp_barrier_primal(s + a * ds),
  Power = 0 + .clk_pow_barrier_dual(K, z + a * dz) + .clk_pow_barrier_primal(K, s + a * ds),
  GenPower = 0 + .clk_gp_barrier_primal(K, s + a * ds) + .clk_gp_barrier_dual(K, z + a * dz),
  PSDTriangle = {
    logdet <- function(x, dx) {
      L <- .clk_chol_lower(.clk_svec_to_mat(x + a * dx, K$n))
      if (is.null(L)) Inf else 2 * .cl_ssum(log(diag(L)))
    }
    0 - logdet(z, dz) - logdet(s, ds)
  })
