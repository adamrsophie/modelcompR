# Weak order stochastic Runge-Kutta methods, after StochasticDiffEqWeak
# (perform_step/srk_weak.jl, caches/srk_weak_caches.jl, weak_utils.jl,
# tableaus.jl): the Debrabant-Roessler and Roessler methods DRI1, RI1, RI3,
# RI5, RI6 and RDI1WM to RDI4WM, Roessler's RS1 and RS2, Platen's PL1WM and
# PL1WMA, Komori's NON, NON2 and COM, the SIE/SME methods of Tocino and
# Vigo-Aguiar, and W2Ito1. The constants are generated into
# R/stochastic-weak-constants.R by tools/sde-gen-weak.jl.
#
# The Brownian increments are replaced by the discrete random variables of
# the weak methods, computed from the Gaussian draws: three point variables
# for dW and two point variables for the auxiliary dZ. Diffusions are used
# column by column (the k-th column is the noise of the k-th Wiener process;
# diagonal noise is the diagonal matrix of g).

#' Weak order SDE methods
#'
#' Algorithms for [solve()] on SDE problems with weak (distributional)
#' order 2 or 1, as in StochasticDiffEq. They match the moments of the
#' solution rather than its path, and replace the Brownian increments by
#' discrete random variables.
#'
#' * `DRI1()`, `RI1()`, `RI3()`, `RI5()`, `RI6()`: Debrabant and Roessler's
#'   weak order 2 methods for Ito SDEs (adaptive for scalar, diagonal or one
#'   dimensional noise); `RDI1WM()` (order 1), `RDI2WM()`, `RDI3WM()`,
#'   `RDI4WM()` (order 2) their methods with fewer stages.
#' * `RS1()`, `RS2()`: Roessler's weak order 2 methods for Stratonovich SDEs.
#' * `PL1WM()`, `PL1WMA()`: Platen's weak order 2 method (for general and for
#'   additive noise).
#' * `NON()`, `NON2()`, `COM()`: Komori's methods for Stratonovich SDEs with
#'   non commutative and commutative noise.
#' * `SIEA()`, `SIEB()`, `SMEA()`, `SMEB()`: Tocino and Vigo-Aguiar's weak
#'   order 2 methods.
#' * `W2Ito1()`: Tang and Xiao's efficient weak order 2 method.
#' * `DRI1NM()`: DRI1 for non mixing diagonal and scalar noise.
#' * `IRI1()`: RI1 with drift implicit (theta method) stages, for stiff drift.
#'
#' With several noise channels the stage values use the k-th column of the
#' diffusion for the k-th Wiener process (the diagonal matrix of `g` for
#' diagonal noise). StochasticDiffEq differs here in three places: `RS1()`
#' and `RS2()` with general noise read the first column of the last two
#' diffusion stages, and `RS1()`, `RS2()` and `PL1WM()` with diagonal noise
#' add the k-th diffusion value to every component. For scalar noise and for
#' one Wiener process the results agree with StochasticDiffEq to rounding.
#'
#' `constructDRI1()` and the other `construct*()` functions return the
#' tableaus (`RoesslerRI`, `RoesslerRS`, `KomoriNON` and `KomoriNON2`
#' objects); `checkRIOrder()`, `checkRSOrder()` and `checkNONOrder()` test
#' the order conditions. `Ihat2()` is the approximation of the iterated
#' integral the methods use.
#'
#' @param tableau A weak tableau.
#' @param tol Tolerance of the order conditions.
#' @return An algorithm, a tableau or a logical vector.
#' @examples
#' prob <- SDEProblem(function(u, p, t) 1.01 * u, function(u, p, t) 0.87 * u, 0.5, c(0, 1))
#' sol <- solve(prob, DRI1(), dt = 0.01, adaptive = FALSE, seed = 1)
#' sol$retcode
#' @name weak_sde_algorithms
NULL

# diffusion as an n x m matrix whose columns are the noise channels
.weak_G <- function(gv, type, n) {
  switch(type, diagonal = diag(as.numeric(gv), n), scalar = matrix(as.numeric(gv), n, 1L), as.matrix(gv))
}

.weak_quantile <- -0.9674215661017014

.weak_alg <- function(name, order, adaptive, step, ...) {
  .new_sde_alg(name, order, adaptive, perform = step, needs_dZ = TRUE, noise = c("diagonal", "scalar", "general"),
               weak = TRUE, ...)
}

# Ihat2 for the DRI1 (Roessler) family
.ihat2_ri <- function(dW, dZ, sqdt, k, l) {
  if (k < l) (dW[k] * dW[l] - sqdt * dZ[k]) / 2
  else if (l < k) (dW[k] * dW[l] + sqdt * dZ[l]) / 2
  else (dW[k]^2 - sqdt^2) / 2
}

# perform_step! for DRI1ConstantCache (shared by RI1, RI3, RI5, RI6, RDI2WM,
# RDI3WM and RDI4WM through their constants)
.dri1_step <- function(C) function(f, g, u, p, t, dt, dW, dZ, type, delta) {
  n <- length(u)
  sqdt <- sqrt(abs(dt))
  sq3dt <- sqrt(3 * dt)
  wdW <- calc_threepoint_random(sq3dt, C$NORMAL_ONESIX_QUANTILE %||% .weak_quantile, dW / sqrt(dt))
  m <- length(wdW)
  chi1 <- (wdW^2 - dt) / 2
  gG <- function(x, tt) .weak_G(g(x, p, tt), type, n)
  k1 <- f(u, p, t)
  G1 <- gG(u, t)
  H02 <- u + C$a021 * k1 * dt + C$b021 * as.numeric(G1 %*% wdW)
  k2 <- f(H02, p, t + C$c02 * dt)
  H03 <- u + C$a032 * k2 * dt + C$a031 * k1 * dt + C$b031 * as.numeric(G1 %*% wdW)
  k3 <- f(H03, p, t + C$c03 * dt)
  unew <- u + C$alpha1 * k1 * dt + C$alpha2 * k2 * dt + C$alpha3 * k3 * dt
  wdZ <- calc_twopoint_random(sqdt, dZ)
  G2 <- vector("list", m); G3 <- vector("list", m)
  for (k in seq_len(m)) {
    G2[[k]] <- gG(u + C$a121 * k1 * dt + C$b121 * G1[, k] * sqdt, t + C$c12 * dt)
    G3[[k]] <- gG(u + C$a131 * k1 * dt + C$b131 * G1[, k] * sqdt, t + C$c13 * dt)
  }
  for (k in seq_len(m)) {
    unew <- unew + G1[, k] * wdW[k] * (C$beta11 + (m - 1) * C$beta31) +
      G2[[k]][, k] * wdW[k] * C$beta12 + G3[[k]][, k] * wdW[k] * C$beta13 +
      G2[[k]][, k] * chi1[k] * C$beta22 / sqdt + G3[[k]][, k] * chi1[k] * C$beta23 / sqdt
    if (m > 1L) for (l in seq_len(m)[-k]) {
      ih <- .ihat2_ri(wdW, wdZ, sqdt, k, l)
      H22l <- u + (C$b221 * G1[, l] + C$b222 * G2[[l]][, l] + C$b223 * G3[[l]][, l]) * sqdt
      H23l <- u + (C$b231 * G1[, l] + C$b232 * G2[[l]][, l] + C$b233 * G3[[l]][, l]) * sqdt
      unew <- unew + gG(H22l, t)[, k] * (wdW[k] * C$beta32 + ih * C$beta42 / sqdt) +
        gG(H23l, t)[, k] * (wdW[k] * C$beta33 + ih * C$beta43 / sqdt)
    }
  }
  # error estimate against a lower order scheme; OrdinaryDiffEq only
  # estimates it for scalar, diagonal or one dimensional noise
  if (!(type %in% c("scalar", "diagonal") || m == 1L)) return(list(u = unew, err = 0 * u, nf = 3L))
  if (C$c03 != 0) {
    rat <- C$c02 / C$c03
    d1 <- delta * (rat - 1); d2 <- -delta * rat
    uhat <- u + ((C$alpha1 + d1) * k1 + (C$alpha2 + delta) * k2 + (C$alpha3 + d2) * k3) * dt
  } else uhat <- u + k1 * dt
  uhat <- uhat + as.numeric(G1 %*% wdW)
  list(u = unew, err = unew - uhat, nf = 3L)
}

.ri_alg <- function(name, order = 2) .weak_alg(name, order, TRUE, .dri1_step(.WEAK_CONSTANTS[[name]]))

#' @rdname weak_sde_algorithms
#' @export
DRI1 <- function() .ri_alg("DRI1")
#' @rdname weak_sde_algorithms
#' @export
RI1 <- function() .ri_alg("RI1")
#' @rdname weak_sde_algorithms
#' @export
RI3 <- function() .ri_alg("RI3")
#' @rdname weak_sde_algorithms
#' @export
RI5 <- function() .ri_alg("RI5")
#' @rdname weak_sde_algorithms
#' @export
RI6 <- function() .ri_alg("RI6")
#' @rdname weak_sde_algorithms
#' @export
RDI2WM <- function() .ri_alg("RDI2WM")
#' @rdname weak_sde_algorithms
#' @export
RDI3WM <- function() .ri_alg("RDI3WM")
#' @rdname weak_sde_algorithms
#' @export
RDI4WM <- function() .ri_alg("RDI4WM")

# perform_step! for RDI1WMConstantCache
#' @rdname weak_sde_algorithms
#' @export
RDI1WM <- function() {
  C <- .WEAK_CONSTANTS$RDI1WM
  .weak_alg("RDI1WM", 1, FALSE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    n <- length(u)
    wdW <- calc_threepoint_random(sqrt(3 * dt), C$NORMAL_ONESIX_QUANTILE, dW / sqrt(dt))
    k1 <- f(u, p, t)
    G1 <- .weak_G(g(u, p, t), type, n)
    H02 <- u + C$a021 * k1 * dt + C$b021 * as.numeric(G1 %*% wdW)
    k2 <- f(H02, p, t + C$c02 * dt)
    list(u = u + C$alpha1 * k1 * dt + C$alpha2 * k2 * dt + C$beta11 * as.numeric(G1 %*% wdW), nf = 2L)
  })
}

# ---------------------------------------------------------------------------
# RS1 and RS2 (perform_step! for RSConstantCache), Stratonovich

.ihat2_rs <- function(dW, dZ, k, l) if (k < l) -dW[k] * dZ[l] else if (l < k) dW[k] * dZ[l] else 0

.rs_step <- function(C) function(f, g, u, p, t, dt, dW, dZ, type, delta) {
  n <- length(u)
  sqdt <- sqrt(abs(dt))
  wdW <- calc_threepoint_random(sqrt(3 * dt), C$NORMAL_ONESIX_QUANTILE, dW / sqrt(dt))
  m <- length(wdW)
  gG <- function(x, tt) .weak_G(g(x, p, tt), type, n)
  k1 <- f(u, p, t)
  G1 <- gG(u, t)
  H02 <- u + C$a021 * k1 * dt
  G2 <- lapply(seq_len(m), function(k) gG(u + C$b121 * G1[, k] * wdW[k], t))
  H13 <- lapply(seq_len(m), function(k) {
    h <- u + C$a131 * k1 * dt + C$b131 * G1[, k] * wdW[k] + C$b132 * G2[[k]][, k] * wdW[k]
    for (l in seq_len(m)[-k]) h <- h + C$b331 * G1[, l] * wdW[l] + C$b332 * G2[[l]][, l] * wdW[l]
    h
  })
  G3 <- lapply(seq_len(m), function(k) gG(H13[[k]], t + C$c13 * dt))
  H14 <- lapply(seq_len(m), function(k) {
    h <- u + C$a141 * k1 * dt + C$b141 * G1[, k] * wdW[k] + C$b142 * G2[[k]][, k] * wdW[k] +
      C$b143 * G3[[k]][, k] * wdW[k]
    for (l in seq_len(m)[-k]) h <- h + C$b341 * G1[, l] * wdW[l] + C$b342 * G2[[l]][, l] * wdW[l]
    h
  })
  G4 <- lapply(seq_len(m), function(k) gG(H14[[k]], t + C$c14 * dt))
  k2 <- f(H02, p, t + C$c02 * dt)
  H03 <- u + C$a032 * k2 * dt + C$a031 * k1 * dt + C$b031 * as.numeric(G1 %*% wdW)
  for (k in seq_len(m)) H03 <- H03 + C$b032 * G2[[k]][, k] * wdW[k]
  k3 <- f(H03, p, t + C$c03 * dt)
  unew <- u + (C$alpha1 * k1 + C$alpha2 * k2 + C$alpha3 * k3 + C$alpha4 * k1) * dt
  wdZ <- calc_twopoint_random(sqdt, dZ)
  for (k in seq_len(m)) {
    unew <- unew + (G1[, k] * C$beta11 + G2[[k]][, k] * C$beta12 + G3[[k]][, k] * C$beta13 +
                      G4[[k]][, k] * C$beta14) * wdW[k]
    if (m > 1L) {
      H22 <- u; H23 <- u
      for (l in seq_len(m)[-k]) {
        ih <- .ihat2_rs(wdW, wdZ, k, l)
        H22 <- H22 + C$b221 * G1[, l] * ih / sqdt
        H23 <- H23 + C$b231 * G1[, l] * ih / sqdt
      }
      unew <- unew + gG(H22, t)[, k] * C$beta22 * sqdt + gG(H23, t)[, k] * C$beta23 * sqdt
    }
  }
  list(u = unew, nf = 3L)
}

#' @rdname weak_sde_algorithms
#' @export
RS1 <- function() .weak_alg("RS1", 2, FALSE, .rs_step(.WEAK_CONSTANTS$RS1), interpretation = "Stratonovich")
#' @rdname weak_sde_algorithms
#' @export
RS2 <- function() .weak_alg("RS2", 2, FALSE, .rs_step(.WEAK_CONSTANTS$RS2), interpretation = "Stratonovich")

# ---------------------------------------------------------------------------
# PL1WM and PL1WMA (Platen)

.ihat2_pl <- function(dW, dZ, sqdt, k, l) {
  if (k < l) -dZ[1 + (l - 3) * l / 2 + k]
  else if (l < k) dZ[1 + (k - 3) * k / 2 + l]
  else -sqdt^2
}

#' @rdname weak_sde_algorithms
#' @export
PL1WM <- function() {
  .weak_alg("PL1WM", 2, FALSE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    n <- length(u)
    sqdt <- sqrt(abs(dt))
    wdW <- calc_threepoint_random(sqrt(3 * dt), .weak_quantile, dW / sqdt)
    m <- length(wdW)
    chi1 <- (wdW^2 - dt) / 4
    gG <- function(x) .weak_G(g(x, p, t), type, n)
    k1 <- f(u, p, t)
    G1 <- gG(u)
    Y <- u + k1 * dt + as.numeric(G1 %*% wdW)
    k2 <- f(Y, p, t)
    unew <- u + (k1 + k2) * dt / 2
    # the auxiliary increments need m (m - 1) / 2 two point variables
    need <- m * (m - 1) / 2
    zz <- if (length(dZ) >= need) dZ else c(dZ, stats::rnorm(need - length(dZ), 0, sqdt))
    wdZ <- calc_twopoint_random(sqdt, zz)
    for (k in seq_len(m)) {
      gp <- gG(u + k1 * dt + G1[, k] * sqdt)
      gm <- gG(u + k1 * dt - G1[, k] * sqdt)
      unew <- unew + (gp[, k] + gm[, k] + 2 * G1[, k]) * wdW[k] / 4 + (gp[, k] - gm[, k]) * chi1[k] / sqdt
      for (l in seq_len(m)[-k]) {
        ih <- .ihat2_pl(wdW, wdZ, sqdt, l, k)
        gp <- gG(u + G1[, l] * sqdt)
        gm <- gG(u - G1[, l] * sqdt)
        unew <- unew + (gp[, k] + gm[, k] - 2 * G1[, k]) * wdW[k] / 4 +
          (gp[, k] - gm[, k]) * (wdW[k] * wdW[l] + ih) / (4 * sqdt)
      }
    }
    list(u = unew, nf = 2L)
  })
}

#' @rdname weak_sde_algorithms
#' @export
PL1WMA <- function() {
  .weak_alg("PL1WMA", 2, FALSE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    wdW <- calc_threepoint_random(sqrt(3 * dt), .weak_quantile, dW / sqrt(abs(dt)))
    k1 <- f(u, p, t)
    tmp1 <- as.numeric(.weak_G(g(u, p, t), type, length(u)) %*% wdW)
    k2 <- f(u + k1 * dt + tmp1, p, t)
    list(u = u + (k1 + k2) * dt / 2 + tmp1, nf = 2L)
  })
}

# ---------------------------------------------------------------------------
# Komori's COM (commutative noise; Gaussian increments)

#' @rdname weak_sde_algorithms
#' @export
COM <- function() {
  C <- .WEAK_CONSTANTS$COM
  .weak_alg("COM", 2, FALSE, interpretation = "Stratonovich", function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    n <- length(u); m <- length(dW)
    gG <- function(x) .weak_G(g(x, p, t), type, n)
    Y10 <- f(u, p, t) * dt
    Y1 <- .weak_G(g(u, p, t), type, n) %*% diag(dW, m)
    Y20 <- u + C$a0021 * Y10
    Y2 <- Y1 * 0
    for (j in seq_len(m)) {
      Y20 <- Y20 + C$a0j21 * Y1[, j]
      Y2[, j] <- gG(u + C$aj021 * Y10 + C$ajj21 * Y1[, j])[, j] * dW[j]
    }
    Y20 <- f(Y20, p, t) * dt
    Y30 <- u + C$a0032 * Y20
    Y3 <- Y1 * 0
    for (j in seq_len(m)) {
      tmpu <- u + C$ajj31 * Y1[, j] + C$ajj32 * Y2[, j]
      Y30 <- Y30 + C$a0j31 * Y1[, j] + C$a0j32 * Y2[, j]
      for (l in seq_len(m)[-j]) tmpu <- tmpu + C$ajl31 * Y1[, l] + C$ajl32 * Y2[, l]
      Y3[, j] <- gG(tmpu)[, j] * dW[j]
    }
    Y30 <- f(Y30, p, t) * dt
    Y40 <- u + C$a0043 * Y30
    Y4 <- Y1 * 0
    for (j in seq_len(m)) {
      tmpu <- u + C$ajj41 * Y1[, j] + C$ajj42 * Y2[, j] + C$ajj43 * Y3[, j] + C$aj041 * Y10
      Y40 <- Y40 + C$a0j41 * Y1[, j]
      for (l in seq_len(m)[-j]) tmpu <- tmpu + C$ajl41 * Y1[, l] + C$ajl42 * Y2[, l]
      Y4[, j] <- gG(tmpu)[, j] * dW[j]
    }
    Y40 <- f(Y40, p, t) * dt
    unew <- u + C$c01 * Y10 + C$c02 * Y20 + C$c03 * Y30 + C$c04 * Y40 +
      rowSums(C$cj1 * Y1 + C$cj2 * Y2 + C$cj3 * Y3 + C$cj4 * Y4)
    list(u = unew, nf = 4L)
  })
}

# ---------------------------------------------------------------------------
# SIEA, SIEB, SMEA, SMEB (perform_step! for SIESMEConstantCache); scalar or
# diagonal noise

.siesme_step <- function(C) function(f, g, u, p, t, dt, dW, dZ, type, delta) {
  if (type == "general") stop("the SIE and SME methods take scalar or diagonal noise", call. = FALSE)
  sqdt <- sqrt(abs(dt))
  W2 <- dW^2 / sqdt
  W3 <- C$nu2 * dW^3 / dt
  k0 <- f(u, p, t); g0 <- g(u, p, t)
  k1 <- f(u + C$lambda0 * k0 * dt + C$nu1 * g0 * dW + g0 * W3, p, t + C$mu0 * dt)
  g1 <- g(u + C$lambdabar0 * k0 * dt + C$beta2 * g0 * sqdt + C$beta3 * g0 * W2, p, t + C$mubar0 * dt)
  g2 <- g(u + C$lambdabar0 * k0 * dt + C$delta2 * g0 * sqdt + C$delta3 * g0 * W2, p, t + C$mubar0 * dt)
  unew <- u + (C$alpha1 * k0 + C$alpha2 * k1) * dt + C$gamma1 * g0 * dW +
    (C$lambda1 * dW + C$lambda2 * sqdt + C$lambda3 * W2) * g1 + (C$mu1 * dW + C$mu2 * sqdt + C$mu3 * W2) * g2
  list(u = unew, nf = 2L)
}

#' @rdname weak_sde_algorithms
#' @export
SIEA <- function() .weak_alg("SIEA", 2, FALSE, .siesme_step(.WEAK_CONSTANTS$SIEA))
#' @rdname weak_sde_algorithms
#' @export
SIEB <- function() .weak_alg("SIEB", 2, FALSE, .siesme_step(.WEAK_CONSTANTS$SIEB))
#' @rdname weak_sde_algorithms
#' @export
SMEA <- function() .weak_alg("SMEA", 2, FALSE, .siesme_step(.WEAK_CONSTANTS$SMEA))
#' @rdname weak_sde_algorithms
#' @export
SMEB <- function() .weak_alg("SMEB", 2, FALSE, .siesme_step(.WEAK_CONSTANTS$SMEB))

# ---------------------------------------------------------------------------
# W2Ito1 (Tang and Xiao)

.w2ito_i2 <- function(dW, xi, eta2, k, l) {
  if (k < l) (dW[l] - eta2 * dW[l]) / 2 else if (k > l) (dW[l] + eta2 * dW[l]) / 2 else (dW[k]^2 / xi - xi) / 2
}

#' @rdname weak_sde_algorithms
#' @export
W2Ito1 <- function() {
  C <- .WEAK_CONSTANTS$W2Ito1
  .weak_alg("W2Ito1", 2, TRUE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    n <- length(u)
    sqdt <- sqrt(abs(dt))
    wdW <- calc_threepoint_random(sqrt(3 * dt), C$NORMAL_ONESIX_QUANTILE, dW / sqrt(dt))
    m <- length(wdW)
    zz <- if (length(dZ) >= 2L) dZ else c(dZ, stats::rnorm(2L - length(dZ)))
    wdZ <- calc_twopoint_random(1, zz)
    xi <- sqdt * wdZ[1L]
    chi1 <- (wdW^2 / xi - xi) / 2
    gG <- function(x) .weak_G(g(x, p, t), type, n)
    k1 <- f(u, p, t)
    G1 <- gG(u)
    H02 <- u + C$a021 * k1 * dt + C$b021 * as.numeric(G1 %*% wdW)
    k2 <- f(H02, p, t)
    H03 <- u + C$a032 * k2 * dt + C$a031 * k1 * dt + C$b031 * as.numeric(G1 %*% wdW)
    k3 <- f(H03, p, t)
    G2 <- G1 * 0; G3 <- G1 * 0
    for (k in seq_len(m)) {
      H12 <- u + C$a121 * k1 * dt + C$b121 * G1[, k] * xi
      H13 <- u + C$a131 * k1 * dt + C$b131 * G1[, k] * xi
      for (l in seq_len(m)[-k]) H12 <- H12 + C$b221 * G1[, l] * .w2ito_i2(wdW, xi, wdZ[2L], k, l)
      G2[, k] <- gG(H12)[, k]
      G3[, k] <- gG(H13)[, k]
    }
    unew <- u + C$alpha1 * k1 * dt + C$alpha2 * k2 * dt + C$alpha3 * k3 * dt +
      as.numeric(G1 %*% (wdW * C$beta01 + chi1 * C$beta11)) + as.numeric(G2 %*% (wdW * C$beta02)) +
      as.numeric(G3 %*% (wdW * C$beta03 + chi1 * C$beta13))
    uhat <- u + k1 * dt + as.numeric(G1 %*% wdW)
    list(u = unew, err = unew - uhat, nf = 3L)
  })
}

# ---------------------------------------------------------------------------
# DRI1NM (perform_step! for DRI1NMCache): DRI1 for non mixing diagonal and
# scalar noise, with one diffusion evaluation per stage

#' @rdname weak_sde_algorithms
#' @export
DRI1NM <- function() {
  C <- .WEAK_CONSTANTS$DRI1
  .weak_alg("DRI1NM", 2, TRUE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    if (type == "general") stop("DRI1NM takes non mixing diagonal or scalar noise", call. = FALSE)
    sqdt <- sqrt(abs(dt))
    wdW <- calc_threepoint_random(sqrt(3 * dt), C$NORMAL_ONESIX_QUANTILE, dW / sqrt(dt))
    chi1 <- (wdW^2 - dt) / 2
    k1 <- f(u, p, t)
    g1 <- as.numeric(g(u, p, t))
    H02 <- u + C$a021 * k1 * dt + C$b021 * g1 * wdW
    k2 <- f(H02, p, t + C$c02 * dt)
    H03 <- u + C$a032 * k2 * dt + C$a031 * k1 * dt + C$b031 * g1 * wdW
    k3 <- f(H03, p, t + C$c03 * dt)
    g2 <- as.numeric(g(u + C$a121 * k1 * dt + C$b121 * g1 * sqdt, p, t + C$c12 * dt))
    g3 <- as.numeric(g(u + C$a131 * k1 * dt + C$b131 * g1 * sqdt, p, t + C$c13 * dt))
    unew <- u + (C$alpha1 * k1 + C$alpha2 * k2 + C$alpha3 * k3) * dt +
      (g1 * C$beta11 + g2 * C$beta12 + g3 * C$beta13) * wdW +
      g2 * chi1 * C$beta22 / sqdt + g3 * chi1 * C$beta23 / sqdt
    if (C$c03 != 0) {
      rat <- C$c02 / C$c03
      uhat <- u + ((C$alpha1 + delta * (rat - 1)) * k1 + (C$alpha2 + delta) * k2 + (C$alpha3 - delta * rat) * k3) * dt
    } else uhat <- u + k1 * dt
    uhat <- uhat + g1 * wdW
    list(u = unew, err = unew - uhat, nf = 3L)
  })
}

# ---------------------------------------------------------------------------
# Komori's NON (perform_step! for NONConstantCache): non commutative noise,
# stage values Y[, ja, jb] for every pair of noise channels

.ihat2_non <- function(dW, dZ, k, l) if (k < l) dZ[k] else if (l < k) dW[l] else dW[k]

#' @rdname weak_sde_algorithms
#' @export
NON <- function() {
  C <- .WEAK_CONSTANTS$NON
  .weak_alg("NON", 2, FALSE, interpretation = "Stratonovich", function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    n <- length(u)
    sqdt <- sqrt(abs(dt))
    wdW <- calc_threepoint_random(sqrt(3 * dt), C$NORMAL_ONESIX_QUANTILE, dW / sqdt)
    wdZ <- calc_twopoint_random(sqdt, dZ)
    m <- length(wdW)
    gG <- function(x) .weak_G(g(x, p, t), type, n)
    ih <- function(a, b) .ihat2_non(wdW, wdZ, a, b)
    ms <- seq_len(m)
    G0 <- gG(u)
    Y100 <- f(u, p, t) * dt
    Y1 <- Y2 <- Y3 <- Y4 <- array(0, c(n, m, m))
    for (ja in ms) for (jb in ms) {
      Y1[, ja, jb] <- ih(ja, jb) * G0[, jb]
      if (ja != jb) Y4[, ja, jb] <- Y1[, ja, jb]
    }
    Y200 <- u + C$a0021 * Y100
    for (ja in ms) {
      Y200 <- Y200 + C$a0j21 * Y1[, ja, ja]
      for (jb in ms) {
        x <- if (jb == ja) u + C$aj021 * Y100 + C$ajj21 * Y1[, ja, ja]
             else if (jb < ja) u + C$aljjl21 * Y1[, jb, ja] else u
        Y2[, ja, jb] <- ih(ja, jb) * gG(x)[, jb]
      }
    }
    Y200 <- f(Y200, p, t) * dt
    Y300 <- u + C$a0032 * Y200
    for (ja in ms) {
      Y300 <- Y300 + C$a0j31 * Y1[, ja, ja] + C$a0j32 * Y2[, ja, ja]
      for (jb in ms) {
        if (jb == ja) {
          x <- u + C$ajj31 * Y1[, ja, ja] + C$ajj32 * Y2[, ja, ja]
          for (l in ms[-ja]) x <- x + C$ajl31 * Y1[, l, l] + C$ajl32 * Y2[, l, l]
        } else x <- if (jb < ja) u + C$aljjl31 * Y1[, jb, ja] else u + C$ajljj31 * Y1[, ja, ja]
        Y3[, ja, jb] <- ih(ja, jb) * gG(x)[, jb]
      }
    }
    Y300 <- f(Y300, p, t) * dt
    Y400 <- u + C$a0043 * Y300
    for (ja in ms) {
      Y400 <- Y400 + C$a0j41 * Y1[, ja, ja]
      x <- u + C$ajj41 * Y1[, ja, ja] + C$ajj42 * Y2[, ja, ja] + C$ajj43 * Y3[, ja, ja] + C$aj041 * Y100
      for (l in ms[-ja]) x <- x + C$ajl41 * Y1[, l, l] + C$ajl42 * Y2[, l, l]
      Y4[, ja, ja] <- ih(ja, ja) * gG(x)[, ja]
    }
    Y400 <- f(Y400, p, t) * dt
    unew <- u + C$c01 * Y100 + C$c02 * Y200 + C$c03 * Y300 + C$c04 * Y400
    for (ja in ms) {
      unew <- unew + C$cj1 * Y1[, ja, ja] + C$cj2 * Y2[, ja, ja] + C$cj3 * Y3[, ja, ja] + C$cj4 * Y4[, ja, ja]
      for (jb in ms[-ja]) {
        unew <- unew + if (ja > jb) C$clj2 * Y2[, ja, jb] + C$clj3 * Y3[, ja, jb]
                       else C$cjl2 * Y2[, ja, jb] + C$cjl3 * Y3[, ja, jb]
      }
    }
    list(u = unew, nf = 4L)
  })
}

# ---------------------------------------------------------------------------
# Komori's NON2 (perform_step! for NON2ConstantCache)

#' @rdname weak_sde_algorithms
#' @export
NON2 <- function() {
  C <- .WEAK_CONSTANTS$NON2
  .weak_alg("NON2", 2, FALSE, interpretation = "Stratonovich", function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    n <- length(u)
    sqdt <- sqrt(abs(dt))
    wdW <- calc_threepoint_random(sqrt(3 * dt), C$NORMAL_ONESIX_QUANTILE, dW / sqdt)
    m <- length(wdW)
    ms <- seq_len(m)
    gG <- function(x) .weak_G(g(x, p, t), type, n)
    G0 <- gG(u)
    Y10 <- f(u, p, t) * dt
    Y1 <- sweep(G0, 2L, wdW, `*`)
    Y20 <- u + C$a0021 * Y10
    Y2 <- Y1 * 0
    for (j in ms) {
      Y20 <- Y20 + C$a0j21 * Y1[, j]
      Y2[, j] <- gG(u + C$aj021 * Y10 + C$ajj21 * Y1[, j])[, j] * wdW[j]
    }
    Y20 <- f(Y20, p, t) * dt
    Y30 <- u + C$a0032 * Y20
    Y3 <- Y1 * 0
    for (j in ms) {
      x <- u + C$ajj31 * Y1[, j] + C$ajj32 * Y2[, j]
      Y30 <- Y30 + C$a0j31 * Y1[, j] + C$a0j32 * Y2[, j]
      for (l in ms[-j]) x <- x + C$ajl31 * Y1[, l] + C$ajl32 * Y2[, l]
      Y3[, j] <- gG(x)[, j] * wdW[j]
    }
    Y30 <- f(Y30, p, t) * dt
    Y40 <- u + C$a0043 * Y30
    Y4 <- Y1 * 0
    for (j in ms) {
      x <- u + C$ajj41 * Y1[, j] + C$ajj42 * Y2[, j] + C$ajj43 * Y3[, j] + C$aj041 * Y10
      Y40 <- Y40 + C$a0j41 * Y1[, j]
      for (l in ms[-j]) x <- x + C$ajl41 * Y1[, l] + C$ajl42 * Y2[, l]
      Y4[, j] <- gG(x)[, j] * wdW[j]
    }
    Y40 <- f(Y40, p, t) * dt
    unew <- u + C$c01 * Y10 + C$c02 * Y20 + C$c03 * Y30 + C$c04 * Y40 +
      rowSums(C$cj1 * Y1 + C$cj2 * Y2 + C$cj3 * Y3 + C$cj4 * Y4)
    if (type == "general" && m > 1L) {
      # the non commutative terms Y3^(k(j)j) and Y4^(k(j)j)
      gam <- C$gamma
      wdZ <- calc_twopoint_random(sqdt, dZ)
      ih <- matrix(0, m, m)
      for (j in ms) for (l in seq_len(j - 1L)) {
        ih[j, l] <- -wdZ[j] * wdW[l] / sqdt
        ih[l, j] <- wdW[l] * wdZ[j] / sqdt
      }
      for (j in ms) for (k in ms) {
        tmp <- as.numeric(G0 %*% ih[k, ]) / (4 * gam)
        unew <- unew + gam * sqdt * (gG(u + tmp)[, k] - gG(u - tmp)[, k])
      }
    }
    list(u = unew, nf = 4L)
  })
}

# ---------------------------------------------------------------------------
# IRI1 (perform_step! for IRI1ConstantCache): RI1 with theta method drift
# stages; each stage z = dt f(tmp + theta z, t + theta dt) is solved by
# Newton's method and the drift is evaluated at tmp + z

#' @rdname weak_sde_algorithms
#' @param theta Implicitness of the drift stages of `IRI1()`.
#' @param autodiff,concrete_jac,linsolve,nlsolve,extrapolant,new_jac_conv_bound
#'   Options of the Julia signature of `IRI1()` (the Jacobian is computed by
#'   finite differences or the problem's `jac`).
#' @export
IRI1 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, nlsolve = NULL,
                 extrapolant = "constant", theta = 1, new_jac_conv_bound = 1e-3) {
  C <- .WEAK_CONSTANTS$RI1
  alg <- .weak_alg("IRI1", 2, TRUE, NULL, theta = theta)
  alg$init <- function(ctx, t, u) {
    cache <- new.env(parent = emptyenv())
    cache$nl <- .new_nlsolver()
    cache
  }
  alg$step <- function(ctx, g, u, t, dt, dW, type, delta, cache, dZ = NULL, ...) {
    p <- ctx$p
    n <- length(u)
    out <- list(nf = 0L, njacs = 0L, nsolve = 0L, nw = 0L, cache = cache, k = NULL)
    stage <- function(tmp) {
      nl <- .nlsolve(cache$nl, ctx, t, dt, tmp, .mass_mul(ctx$mass, tmp) / (theta * dt), theta, 1,
                     t + theta * dt, u)
      for (s in c("nf", "njacs", "nsolve", "nw")) out[[s]] <<- out[[s]] + as.integer(nl[[s]] %||% 0L)
      if (!nl$ok) return(NULL)
      tmp + (nl$z - tmp) / theta
    }
    sqdt <- sqrt(abs(dt))
    wdW <- calc_threepoint_random(sqrt(3 * dt), C$NORMAL_ONESIX_QUANTILE, dW / sqrt(dt))
    m <- length(wdW)
    chi1 <- (wdW^2 - dt) / 2
    gG <- function(x, tt) .weak_G(g(x, p, tt), type, n)
    Y1 <- stage(u)
    if (is.null(Y1)) return(c(out, list(fail = TRUE, u = u)))
    k1 <- ctx$f(Y1, p, t)
    G1 <- gG(u, t)
    gw <- as.numeric(G1 %*% wdW)
    Y2 <- stage(u + C$a021 * k1 * dt + C$b021 * gw)
    if (is.null(Y2)) return(c(out, list(fail = TRUE, u = u)))
    k2 <- ctx$f(Y2, p, t + C$c02 * dt)
    Y3 <- stage(u + C$a032 * k2 * dt + C$a031 * k1 * dt + C$b031 * gw)
    if (is.null(Y3)) return(c(out, list(fail = TRUE, u = u)))
    k3 <- ctx$f(Y3, p, t + C$c03 * dt)
    out$nf <- out$nf + 3L
    unew <- u + C$alpha1 * k1 * dt + C$alpha2 * k2 * dt + C$alpha3 * k3 * dt
    wdZ <- calc_twopoint_random(sqdt, dZ %||% numeric(m))
    G2 <- lapply(seq_len(m), function(k) gG(u + C$a121 * k1 * dt + C$b121 * G1[, k] * sqdt, t + C$c12 * dt))
    G3 <- lapply(seq_len(m), function(k) gG(u + C$a131 * k1 * dt + C$b131 * G1[, k] * sqdt, t + C$c13 * dt))
    for (k in seq_len(m)) {
      unew <- unew + G1[, k] * wdW[k] * (C$beta11 + (m - 1) * C$beta31) +
        G2[[k]][, k] * wdW[k] * C$beta12 + G3[[k]][, k] * wdW[k] * C$beta13 +
        G2[[k]][, k] * chi1[k] * C$beta22 / sqdt + G3[[k]][, k] * chi1[k] * C$beta23 / sqdt
      for (l in seq_len(m)[-k]) {
        ih <- .ihat2_ri(wdW, wdZ, sqdt, k, l)
        H22 <- u + (C$b221 * G1[, l] + C$b222 * G2[[l]][, l] + C$b223 * G3[[l]][, l]) * sqdt
        H23 <- u + (C$b231 * G1[, l] + C$b232 * G2[[l]][, l] + C$b233 * G3[[l]][, l]) * sqdt
        unew <- unew + gG(H22, t)[, k] * (wdW[k] * C$beta32 + ih * C$beta42 / sqdt) +
          gG(H23, t)[, k] * (wdW[k] * C$beta33 + ih * C$beta43 / sqdt)
      }
    }
    if (!isTRUE(ctx$adaptive) || !(type %in% c("scalar", "diagonal") || m == 1L)) return(c(out, list(u = unew)))
    if (C$c03 != 0) {
      rat <- C$c02 / C$c03
      uhat <- u + ((C$alpha1 + delta * (rat - 1)) * k1 + (C$alpha2 + delta) * k2 + (C$alpha3 - delta * rat) * k3) * dt
    } else uhat <- u + k1 * dt
    c(out, list(u = unew, err = unew - (uhat + gw)))
  }
  alg
}

# ---------------------------------------------------------------------------
# Iterated integral approximations (weak_utils.jl and the methods' Ihat2)

#' @rdname weak_sde_algorithms
#' @param method A weak algorithm (or its name) whose approximation `Ihat2()`
#'   evaluates.
#' @param dW,dZ The discrete random variables of the step.
#' @param sqdt The square root of the step size.
#' @param k,l Indices of the noise channels.
#' @export
Ihat2 <- function(method, dW, dZ, sqdt, k, l) {
  nm <- if (is.character(method)) method else method$name
  if (nm %in% c("DRI1", "DRI1NM", "RI1", "RI3", "RI5", "RI6", "RDI1WM", "RDI2WM", "RDI3WM", "RDI4WM", "IRI1"))
    return(.ihat2_ri(dW, dZ, sqdt, k, l))
  if (nm %in% c("RS1", "RS2")) return(.ihat2_rs(dW, dZ, k, l))
  if (nm == "PL1WM") return(.ihat2_pl(dW, dZ, sqdt, k, l))
  if (nm == "NON") return(.ihat2_non(dW, dZ, k, l))
  stop("no Ihat2 for ", nm, call. = FALSE)
}

# ---------------------------------------------------------------------------
# Tableaus (tableaus.jl): RoesslerRI, RoesslerRS, KomoriNON, KomoriNON2

#' @rdname weak_sde_algorithms
#' @param c0,c1,c2,A0,A1,A2,B0,B1,B2,B3,alpha,beta1,beta2,beta3,beta4,order,quantile
#'   Coefficients of a Roessler weak tableau.
#' @export
RoesslerRI <- function(c0, c1, c2, A0, A1, A2, B0, B1, B2, alpha, beta1, beta2, beta3, beta4, order = 2,
                       quantile = .weak_quantile)
  structure(list(c0 = c0, c1 = c1, c2 = c2, A0 = as.matrix(A0), A1 = as.matrix(A1), A2 = as.matrix(A2),
                 B0 = as.matrix(B0), B1 = as.matrix(B1), B2 = as.matrix(B2), alpha = alpha, beta1 = beta1,
                 beta2 = beta2, beta3 = beta3, beta4 = beta4, order = order, quantile = quantile),
            class = c("RoesslerRI", "Tableau"))

#' @rdname weak_sde_algorithms
#' @export
RoesslerRS <- function(c0, c1, c2, A0, A1, A2, B0, B1, B2, B3, alpha, beta1, beta2, order = 2,
                       quantile = .weak_quantile)
  structure(list(c0 = c0, c1 = c1, c2 = c2, A0 = as.matrix(A0), A1 = as.matrix(A1), A2 = as.matrix(A2),
                 B0 = as.matrix(B0), B1 = as.matrix(B1), B2 = as.matrix(B2), B3 = as.matrix(B3), alpha = alpha,
                 beta1 = beta1, beta2 = beta2, order = order, quantile = quantile),
            class = c("RoesslerRS", "Tableau"))

#' @rdname weak_sde_algorithms
#' @param cj,cjl,clj,ckj,alpha00,alpha0j,alphaj0,alphajj,alphajl,alphajljj,alphaljjl,alphakjjl
#'   Coefficients of a Komori tableau.
#' @export
KomoriNON <- function(c0, cj, cjl, clj, alpha00, alpha0j, alphaj0, alphajj, alphajl, alphajljj, alphaljjl,
                      order = 2, quantile = .weak_quantile)
  structure(list(c0 = c0, cj = cj, cjl = cjl, clj = clj, alpha00 = as.matrix(alpha00),
                 alpha0j = as.matrix(alpha0j), alphaj0 = as.matrix(alphaj0), alphajj = as.matrix(alphajj),
                 alphajl = as.matrix(alphajl), alphajljj = as.matrix(alphajljj), alphaljjl = as.matrix(alphaljjl),
                 order = order, quantile = quantile), class = c("KomoriNON", "Tableau"))

#' @rdname weak_sde_algorithms
#' @export
KomoriNON2 <- function(c0, cj, ckj, alpha00, alpha0j, alphaj0, alphajj, alphajl, alphakjjl, order = 2,
                       quantile = .weak_quantile)
  structure(list(c0 = c0, cj = cj, ckj = ckj, alpha00 = as.matrix(alpha00), alpha0j = as.matrix(alpha0j),
                 alphaj0 = as.matrix(alphaj0), alphajj = as.matrix(alphajj), alphajl = as.matrix(alphajl),
                 alphakjjl = as.matrix(alphakjjl), order = order, quantile = quantile),
            class = c("KomoriNON2", "Tableau"))

.weak_tab <- function(name) {
  tb <- .SDE_TABLEAUS[[name]]
  ctor <- if (name %in% c("RS1", "RS2")) RoesslerRS else if (name == "NON") KomoriNON else
    if (name == "NON2") KomoriNON2 else RoesslerRI
  do.call(ctor, tb[intersect(names(tb), names(formals(ctor)))])
}

#' @rdname weak_sde_algorithms
#' @export
constructDRI1 <- function() .weak_tab("DRI1")
#' @rdname weak_sde_algorithms
#' @export
constructRI1 <- function() .weak_tab("RI1")
#' @rdname weak_sde_algorithms
#' @export
constructRI3 <- function() .weak_tab("RI3")
#' @rdname weak_sde_algorithms
#' @export
constructRI5 <- function() .weak_tab("RI5")
#' @rdname weak_sde_algorithms
#' @export
constructRI6 <- function() .weak_tab("RI6")
#' @rdname weak_sde_algorithms
#' @export
constructRDI1WM <- function() .weak_tab("RDI1WM")
#' @rdname weak_sde_algorithms
#' @export
constructRDI2WM <- function() .weak_tab("RDI2WM")
#' @rdname weak_sde_algorithms
#' @export
constructRDI3WM <- function() .weak_tab("RDI3WM")
#' @rdname weak_sde_algorithms
#' @export
constructRDI4WM <- function() .weak_tab("RDI4WM")
#' @rdname weak_sde_algorithms
#' @export
constructRS1 <- function() .weak_tab("RS1")
#' @rdname weak_sde_algorithms
#' @export
constructRS2 <- function() .weak_tab("RS2")
#' @rdname weak_sde_algorithms
#' @export
constructNON <- function() .weak_tab("NON")
#' @rdname weak_sde_algorithms
#' @export
constructNON2 <- function() .weak_tab("NON2")

.wdot <- function(a, b) sum(as.numeric(a) * as.numeric(b))

# Order conditions (checkRIOrder, checkRSOrder, checkNONOrder of tableaus.jl)

#' @rdname weak_sde_algorithms
#' @param ps The weak order whose conditions are checked (1 or 2).
#' @export
checkRIOrder <- function(tableau, tol = 1e-6, ps = 2) {
  if (!ps %in% c(1, 2)) stop("only the conditions for first and second order weak convergence are implemented; choose ps = 1 or ps = 2", call. = FALSE)
  with(tableau, {
    e <- rep(1, length(alpha))
    cc <- logical(if (ps == 2) 59L else 9L)
    cc[1] <- abs(.wdot(alpha, e) - 1) < tol
    cc[2] <- abs(.wdot(beta4, e) - 0) < tol
    cc[3] <- abs(.wdot(beta3, e) - 0) < tol
    cc[4] <- abs(.wdot(beta1, e)^2 - 1) < tol
    cc[5] <- abs(.wdot(beta2, e) - 0) < tol
    cc[6] <- abs(.wdot(beta1, (B1 %*% e)) - 0) < tol
    cc[7] <- abs(.wdot(beta4, (A2 %*% e)) - 0) < tol
    cc[8] <- abs(.wdot(beta3, (B2 %*% e)) - 0) < tol
    cc[9] <- abs(.wdot(beta4, (B2 %*% e)^2) - 0) < tol
    if (ps == 2) {
      cc[10] <- abs(.wdot(alpha, (A0 %*% e)) - 0.5) < tol
      cc[11] <- abs(.wdot(alpha, (B0 %*% e)^2) - 0.5) < tol
      cc[12] <- abs(.wdot(beta1, e) %*% .wdot(alpha, (B0 %*% e)) - 0.5) < tol
      cc[13] <- abs(.wdot(beta1, e) %*% .wdot(beta1, (A1 %*% e)) - 0.5) < tol
      cc[14] <- abs(.wdot(beta3, (A2 %*% e)) - 0) < tol
      cc[15] <- abs(.wdot(beta2, (B1 %*% e)) - 1) < tol
      cc[16] <- abs(.wdot(beta4, (B2 %*% e)) - 1) < tol
      cc[17] <- abs(.wdot(beta1, e) %*% .wdot(beta1, (B1 %*% e)^2) - 0.5) < tol
      cc[18] <- abs(.wdot(beta1, e) %*% .wdot(beta3, (B2 %*% e)^2) - 0.5) < tol
      cc[19] <- abs(.wdot(beta1, B1 %*% (B1 %*% e)) - 0) < tol
      cc[20] <- abs(.wdot(beta3, B2 %*% (B1 %*% e)) - 0) < tol
      cc[21] <- abs(.wdot(beta3, B2 %*% (B1 %*% (B1 %*% e))) - 0) < tol
      cc[22] <- abs(.wdot(beta1, A1 %*% (B0 %*% e)) - 0) < tol
      cc[23] <- abs(.wdot(beta3, A2 %*% (B0 %*% e)) - 0) < tol
      cc[24] <- abs(.wdot(beta4, (A2 %*% e)^2) - 0) < tol
      cc[25] <- abs(.wdot(beta4, A2 %*% (A0 %*% e)) - 0) < tol
      cc[26] <- abs(.wdot(alpha, B0 %*% (B1 %*% e)) - 0) < tol
      cc[27] <- abs(.wdot(beta2, A1 %*% e) - 0) < tol
      cc[28] <- abs(.wdot(beta1, (A1 %*% e) * (B1 %*% e)) - 0) < tol
      cc[29] <- abs(.wdot(beta3, (A2 %*% e) * (B2 %*% e)) - 0) < tol
      cc[30] <- abs(.wdot(beta4, A2 %*% (B0 %*% e)) - 0) < tol
      cc[31] <- abs(.wdot(beta2, A1 %*% (B0 %*% e)) - 0) < tol
      cc[32] <- abs(.wdot(beta4, ((B2 %*% e)^2) * (A2 %*% e)) - 0) < tol
      cc[33] <- abs(.wdot(beta4, A2 %*% (B0 %*% e)^2) - 0) < tol
      cc[34] <- abs(.wdot(beta2, A1 %*% (B0 %*% e)^2) - 0) < tol
      cc[35] <- abs(.wdot(beta1, B1 %*% (A1 %*% e)) - 0) < tol
      cc[36] <- abs(.wdot(beta3, B2 %*% (A1 %*% e)) - 0) < tol
      cc[37] <- abs(.wdot(beta2, (B1 %*% e)^2) - 0) < tol
      cc[38] <- abs(.wdot(beta4, B2 %*% (B1 %*% e)) - 0) < tol
      cc[39] <- abs(.wdot(beta2, B1 %*% (B1 %*% e)) - 0) < tol
      cc[40] <- abs(.wdot(beta1, (B1 %*% e)^3) - 0) < tol
      cc[41] <- abs(.wdot(beta3, (B2 %*% e)^3) - 0) < tol
      cc[42] <- abs(.wdot(beta1, B1 %*% ((B1 %*% e)^2)) - 0) < tol
      cc[43] <- abs(.wdot(beta3, B2 %*% ((B1 %*% e)^2)) - 0) < tol
      cc[44] <- abs(.wdot(beta4, (B2 %*% e)^4) - 0) < tol
      cc[45] <- abs(.wdot(beta4, (B2 %*% (B1 %*% e))^2) - 0) < tol
      cc[46] <- abs(.wdot(beta4, (B2 %*% e) * (B2 %*% (B1 %*% e))) - 0) < tol
      cc[47] <- abs(.wdot(alpha, (B0 %*% e) * (B0 %*% (B1 %*% e))) - 0) < tol
      cc[48] <- abs(.wdot(beta1, (A1 %*% (B0 %*% e)) * (B1 %*% e)) - 0) < tol
      cc[49] <- abs(.wdot(beta3, (A2 %*% (B0 %*% e)) * (B2 %*% e)) - 0) < tol
      cc[50] <- abs(.wdot(beta1, A1 %*% (B0 %*% (B1 %*% e))) - 0) < tol
      cc[51] <- abs(.wdot(beta3, A2 %*% (B0 %*% (B1 %*% e))) - 0) < tol
      cc[52] <- abs(.wdot(beta4, (B2 %*% (A1 %*% e)) * (B2 %*% e)) - 0) < tol
      cc[53] <- abs(.wdot(beta1, B1 %*% (A1 %*% (B0 %*% e))) - 0) < tol
      cc[54] <- abs(.wdot(beta3, B2 %*% (A1 %*% (B0 %*% e))) - 0) < tol
      cc[55] <- abs(.wdot(beta1, (B1 %*% e) * (B1 %*% (B1 %*% e))) - 0) < tol
      cc[56] <- abs(.wdot(beta3, (B2 %*% e) * (B2 %*% (B1 %*% e))) - 0) < tol
      cc[57] <- abs(.wdot(beta1, B1 %*% (B1 %*% (B1 %*% e))) - 0) < tol
      cc[58] <- abs(.wdot(beta4, (B2 %*% e) * (B2 %*% ((B1 %*% e)^2))) - 0) < tol
      cc[59] <- abs(.wdot(beta4, (B2 %*% e) * (B2 %*% (B1 %*% (B1 %*% e)))) - 0) < tol
    }
    as.logical(cc)
  })
}

#' @rdname weak_sde_algorithms
#' @export
checkRSOrder <- function(tableau, tol = 1e-6, ps = 2) {
  if (!ps %in% c(1, 2)) stop("only the conditions for first and second order weak convergence are implemented; choose ps = 1 or ps = 2", call. = FALSE)
  with(tableau, {
    e <- rep(1, length(alpha))
    cc <- logical(if (ps == 2) 55L else 6L)
    cc[1] <- abs(.wdot(alpha, e) - 1) < tol
    cc[2] <- abs(.wdot(beta1, e)^2 - 1) < tol
    cc[3] <- abs(.wdot(beta2, e) - 0) < tol
    cc[4] <- abs(.wdot(beta1, (B1 %*% e)) - 0.5) < tol
    cc[5] <- abs(.wdot(beta2, (A2 %*% e)) - 0) < tol
    cc[6] <- abs(.wdot(beta2, (B2 %*% e)^2) - 0) < tol
    if (ps == 2) {
      cc[7] <- abs(.wdot(alpha, (A0 %*% e)) - 0.5) < tol
      cc[8] <- abs(.wdot(alpha, B0 %*% (B1 %*% e)) - 0.25) < tol
      cc[9] <- abs(.wdot(alpha, (B0 %*% e)^2) - 0.5) < tol
      cc[10] <- abs(.wdot(beta1, e) %*% .wdot(alpha, (B0 %*% e)) - 0.5) < tol
      cc[11] <- abs(.wdot(beta1, e) %*% .wdot(beta1, (A1 %*% e)) - 0.5) < tol
      cc[12] <- abs(.wdot(beta1, B1 %*% (A1 %*% e)) - 0.25) < tol
      cc[13] <- abs(.wdot(beta1, (B1 %*% e) * (A1 %*% e)) - 0.25) < tol
      cc[14] <- abs(.wdot(beta1, (B3 %*% e)) - 0.5) < tol
      cc[15] <- abs(.wdot(beta1, e) %*% .wdot(beta1, (B1 %*% e)^2) - 1 / 3) < tol
      cc[16] <- abs(.wdot(beta1, e) %*% .wdot(beta1, (B3 %*% e)^2) - 0.5) < tol
      cc[17] <- abs(.wdot(beta1, B3 %*% (B3 %*% e)) - 0) < tol
      cc[18] <- abs(.wdot(beta2, (B2 %*% e))^2 - 0.25) < tol
      cc[19] <- abs(.wdot(beta1, (B1 %*% e)^3) - 0.25) < tol
      cc[20] <- abs(.wdot(beta1, B1 %*% (B1 %*% e)^2) - 1 / 12) < tol
      cc[21] <- abs(.wdot(beta1, B1 %*% (B3 %*% e)^2) - 0.25) < tol
      cc[22] <- abs(.wdot(beta1, A1 %*% (B0 %*% e)) - 0) < tol
      cc[23] <- abs(.wdot(beta2, (A2 %*% e)^2) - 0) < tol
      cc[24] <- abs(.wdot(beta2, A2 %*% (A0 %*% e)) - 0) < tol
      cc[25] <- abs(.wdot(beta1, B1 %*% B1 %*% (B1 %*% e)) - 1 / 24) < tol
      cc[26] <- abs(.wdot(beta2, A2 %*% (B0 %*% e)) - 0) < tol
      cc[27] <- abs(.wdot(beta2, A2 %*% (B0 %*% e)^2) - 0) < tol
      cc[28] <- abs(.wdot(beta2, (B2 %*% e)^4) - 0) < tol
      cc[29] <- abs(.wdot(beta2, (B2 %*% (B1 %*% e))^2) - 0) < tol
      cc[30] <- abs(.wdot(beta2, (B2 %*% (B3 %*% e))^2) - 0) < tol
      cc[31] <- abs(.wdot(beta1, (B1 %*% e) * (B3 %*% e)^2) - 0.25) < tol
      cc[32] <- abs(.wdot(beta2, (A2 %*% e) * (B2 %*% e)^2) - 0) < tol
      cc[33] <- abs(.wdot(beta1, B1 %*% B3 %*% (B1 %*% e)) - 1 / 8) < tol
      cc[34] <- abs(.wdot(beta1, B3 %*% B3 %*% (B3 %*% e)) - 0) < tol
      cc[35] <- abs(.wdot(beta1, B3 %*% B1 %*% (B3 %*% e)) - 0) < tol
      cc[36] <- abs(.wdot(beta2, A2 %*% B0 %*% (B1 %*% e)) - 0) < tol
      cc[37] <- abs(.wdot(beta1, e) %*% .wdot(beta1, (B3 %*% e) * (B1 %*% e)) - 0.25) < tol
      cc[38] <- abs(.wdot(beta1, e) %*% .wdot(beta1, B1 %*% (B1 %*% e)) - 1 / 6) < tol
      cc[39] <- abs(.wdot(beta1, e) %*% .wdot(beta1, B3 %*% (B1 %*% e)) - 0.25) < tol
      cc[40] <- abs(.wdot(beta1, e) %*% .wdot(beta1, B1 %*% (B3 %*% e)) - 0.25) < tol
      cc[41] <- abs(.wdot(beta1, (B1 %*% e) * (B1 %*% (B1 %*% e))) - 1 / 8) < tol
      cc[42] <- abs(.wdot(beta1, (B1 %*% e) * (B3 %*% (B1 %*% e))) - 1 / 8) < tol
      cc[43] <- abs(.wdot(beta1, (B3 %*% e) * (B1 %*% (B3 %*% e))) - 1 / 4) < tol
      cc[44] <- abs(.wdot(beta1, (B3 %*% e) * (B3 %*% (B3 %*% e))) - 0) < tol
      cc[45] <- abs(.wdot(beta1, B3 %*% ((B3 %*% e) * (B1 %*% e))) - 0) < tol
      cc[46] <- abs(.wdot(beta2, (B2 %*% (A1 %*% e)) * (B2 %*% e)) - 0) < tol
      cc[47] <- abs(.wdot(beta2, (B2 %*% e) * (B2 %*% (B1 %*% e))) - 0) < tol
      cc[48] <- abs(.wdot(beta2, (B2 %*% e) * (B2 %*% (B3 %*% e))) - 0) < tol
      cc[49] <- abs(.wdot(beta2, (B2 %*% e) * (B2 %*% ((B1 %*% e)^2))) - 0) < tol
      cc[50] <- abs(.wdot(beta2, (B2 %*% e) * (B2 %*% ((B3 %*% e)^2))) - 0) < tol
      cc[51] <- abs(.wdot(beta2, (B2 %*% e) * (B2 %*% ((B1 %*% e) * (B3 %*% e)))) - 0) < tol
      cc[52] <- abs(.wdot(beta2, (B2 %*% e) * (B2 %*% B1 %*% (B1 %*% e))) - 0) < tol
      cc[53] <- abs(.wdot(beta2, (B2 %*% e) * (B2 %*% B3 %*% (B1 %*% e))) - 0) < tol
      cc[54] <- abs(.wdot(beta2, (B2 %*% e) * (B2 %*% B3 %*% (B3 %*% e))) - 0) < tol
      cc[55] <- abs(.wdot(beta2, (B2 %*% e) * (B2 %*% B1 %*% (B3 %*% e))) - 0) < tol
    }
    as.logical(cc)
  })
}
#' @rdname weak_sde_algorithms
#' @export
checkNONOrder <- function(tableau, tol = 1e-6) {
  with(tableau, {
    e <- rep(1, length(c0))
    cc <- logical(44L)
    cc[1] <- abs(.wdot(c0, e) - 1) < tol
    cc[2] <- abs(.wdot(cj, e) - 1) < tol
    cc[3] <- abs(.wdot(cj, alphajj %*% e) - 1 / 2) < tol
    cc[4] <- abs(.wdot(cj, alpha0j %*% e) - 1 / 2) < tol
    cc[5] <- abs(.wdot(cj, alphaj0 %*% e) - 1 / 2) < tol
    cc[6] <- abs(.wdot(c0, alpha00 %*% e) - 1 / 2) < tol
    cc[7] <- abs(.wdot(cj, alphajj %*% (alphaj0 %*% e)) - 1 / 4) < tol
    cc[8] <- abs(.wdot(c0, alpha0j %*% (alphajj %*% e)) - 1 / 4) < tol
    cc[9] <- abs(.wdot(cj, alphaj0 %*% (alpha0j %*% e)) - 0) < tol
    cc[10] <- abs(.wdot(c0, (alpha0j %*% e)^2) - 1 / 2) < tol
    cc[11] <- abs(.wdot(cj, (alphaj0 %*% e) * (alphajj %*% e)) - 1 / 4) < tol
    cc[12] <- abs(.wdot(cj, alphajj %*% alphajj %*% (alphajj %*% e)) - 1 / 24) < tol
    cc[13] <- abs(.wdot(cj, alphajj %*% (alphajj %*% e)^2) - 1 / 12) < tol
    cc[14] <- abs(.wdot(cj * (alphajj %*% e), alphajj %*% (alphajj %*% e)) - 1 / 8) < tol
    cc[15] <- abs(.wdot(cj, (alphajj %*% e)^3) - 1 / 4) < tol
    cc[16] <- abs(.wdot(cj, alphajj %*% (alphajj %*% e)) - 1 / 6) < tol
    cc[17] <- abs(.wdot(cj, (alphajj %*% e)^2) - 1 / 3) < tol
    cc[18] <- abs(.wdot(cj, (alphajl %*% e)) - 1 / 2) < tol
    cc[19] <- abs(.wdot(cj * (alphajl %*% e), alphajl %*% (alphajl %*% e)) - 0) < tol
    cc[20] <- abs(.wdot(cj, (alphajl %*% e)^2) - 1 / 2) < tol
    cc[21] <- abs(.wdot(cj, alphajj %*% alphajl %*% (alphajj %*% e)) - 1 / 8) < tol
    cc[22] <- abs(.wdot(cj, alphajl %*% alphajl %*% (alphajl %*% e)) - 0) < tol
    cc[23] <- abs(.wdot(cj, alphajl %*% alphajj %*% (alphajl %*% e)) - 0) < tol
    cc[24] <- abs(.wdot(cj, alphajj %*% (alphajl %*% e)^2) - 1 / 4) < tol
    cc[25] <- abs(.wdot(cj, alphajl %*% ((alphajj %*% e) * (alphajl %*% e))) - 0) < tol
    cc[26] <- abs(.wdot(cj * (alphajj %*% e), alphajl %*% (alphajj %*% e)) - 1 / 8) < tol
    cc[27] <- abs(.wdot(cj * (alphajl %*% e), alphajj %*% (alphajl %*% e)) - 1 / 4) < tol
    cc[28] <- abs(.wdot(cj * (alphajj %*% e), (alphajl %*% e)^2) - 1 / 4) < tol
    cc[29] <- abs(.wdot(cj, alphajj %*% (alphajl %*% e)) - 1 / 4) < tol
    cc[30] <- abs(.wdot(cj, alphajl %*% (alphajj %*% e)) - 1 / 4) < tol
    cc[31] <- abs(.wdot(cj, alphajl %*% (alphajl %*% e)) - 0) < tol
    cc[32] <- abs(.wdot(cj, (alphajj %*% e) * (alphajl %*% e)) - 1 / 4) < tol
    if (inherits(tableau, "KomoriNON")) {
      cc[33] <- abs(.wdot(clj, e) - 0) < tol
      cc[34] <- abs(.wdot(cjl, e) - 0) < tol
      cc[35] <- abs(.wdot(clj, alphaljjl %*% e) - 1 / 2) < tol
      cc[36] <- abs(.wdot(cjl, alphajljj %*% e) + 1 / 2) < tol
      cc[37] <- abs(.wdot(clj * (alphaljjl %*% e), alphaljjl %*% (alphajljj %*% e)) - 0) < tol
      cc[38] <- abs(.wdot(clj, (alphaljjl %*% e)^2) - 0) < tol
    } else if (inherits(tableau, "KomoriNON2")) {
      cc[33] <- abs(ckj[4] + ckj[3] - 0) < tol
      cc[34] <- abs(ckj[2] - 0) < tol
      cc[35] <- abs(alphakjjl[4, 3] - 0) < tol
      cc[36] <- TRUE
      cc[37] <- abs(ckj[3] %*% alphakjjl[3, 2]^2 + ckj[4] %*% alphakjjl[4, 2]^2 - 0) < tol
      cc[38] <- abs(ckj[3] %*% alphakjjl[3, 2] + ckj[4] %*% alphakjjl[4, 2] - 1 / 2) < tol
    } else {
      # tableaus without the extra noise coefficients cannot satisfy conditions 33 to 38
      cc[33:38] <- FALSE
    }
    cc[39] <- abs(.wdot(c0, alpha00 %*% alpha00 %*% (alpha00 %*% e)) - 1 / 24) < tol
    cc[40] <- abs(.wdot(c0, alpha00 %*% (alpha00 %*% e)^2) - 1 / 12) < tol
    cc[41] <- abs(.wdot(c0 * (alpha00 %*% e), alpha00 %*% (alpha00 %*% e)) - 1 / 8) < tol
    cc[42] <- abs(.wdot(c0, (alpha00 %*% e)^3) - 1 / 4) < tol
    cc[43] <- abs(.wdot(c0, alpha00 %*% (alpha00 %*% e)) - 1 / 6) < tol
    cc[44] <- abs(.wdot(c0, (alpha00 %*% e)^2) - 1 / 3) < tol
    as.logical(cc)
  })
}
