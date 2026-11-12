# Pure R SDE integrators, after StochasticDiffEq.jl (perform_step/low_order.jl
# for EM, EulerHeun and RKMil, perform_step/lamba.jl for LambaEM,
# perform_step/sri.jl and caches/rossler_caches.jl for SRIW1 and SOSRI,
# perform_step/sdirk.jl and caches/sdirk_caches.jl for ImplicitEM,
# ImplicitEulerHeun and ImplicitRKMil, default_sde_alg.jl, alg_utils.jl for
# orders, error estimate weights and alg_mass_matrix_compatible) and
# DiffEqNoiseProcess.jl (Wiener increments with rejection sampling with
# memory, RSwM).
#
# The steps run inside the core ODE driver .integrate(): each solve wraps the
# chosen algorithm in a one-off ODE algorithm whose step draws the Wiener
# increments. The driver gives the step only (t, dt), so the noise process
# classifies every call against the previous attempt: a call from the same
# time means the attempt was rejected (its increment goes back on a stack and
# is split with the Brownian bridge when a shorter step is tried), a call
# from the end of the attempt means it was accepted, and a call from inside
# it means an event truncated it (the bridged part up to the event is kept
# and the remainder goes back on the stack, as StochasticDiffEq does when a
# callback changes t).

# ---------------------------------------------------------------------------
# Noise process

.noise_new <- function(t0, m, needs_dZ) {
  ns <- new.env(parent = emptyenv())
  ns$m <- m
  ns$needs_dZ <- needs_dZ
  ns$stack <- list()
  ns$pending <- NULL
  ns$W <- numeric(m)
  ns$log_t <- t0
  ns$log_W <- list(numeric(m))
  ns
}

# Tolerance for comparing times: a few rounding errors of t.
.noise_tol <- function(t) 8 * .Machine$double.eps * max(1, abs(t))

.noise_commit <- function(ns, tend, dW) {
  ns$W <- ns$W + dW
  ns$log_t <- c(ns$log_t, tend)
  ns$log_W[[length(ns$log_W) + 1L]] <- ns$W
}

# Split an increment over [0, h] at s with the Brownian bridge.
.noise_split <- function(ns, e, s) {
  frac <- s / e$h
  sd <- sqrt(s * (e$h - s) / e$h)
  w1 <- stats::rnorm(ns$m, e$dW * frac, sd)
  z1 <- if (ns$needs_dZ) stats::rnorm(ns$m, e$dZ * frac, sd) else NULL
  list(first = list(h = s, dW = w1, dZ = z1),
       rest = list(h = e$h - s, dW = e$dW - w1, dZ = if (ns$needs_dZ) e$dZ - z1 else NULL))
}

# Settle the previous attempt given the time of the next call.
.noise_settle <- function(ns, t) {
  pd <- ns$pending
  if (is.null(pd)) return(invisible())
  ns$pending <- NULL
  tol <- .noise_tol(pd$t + pd$h)
  if (abs(t - pd$t) <= tol) {
    ns$stack <- c(list(pd[c("h", "dW", "dZ")]), ns$stack)
  } else if (abs(t - (pd$t + pd$h)) <= tol) {
    .noise_commit(ns, pd$t + pd$h, pd$dW)
  } else if (t > pd$t && t < pd$t + pd$h) {
    sp <- .noise_split(ns, pd, t - pd$t)
    .noise_commit(ns, t, sp$first$dW)
    ns$stack <- c(list(sp$rest), ns$stack)
  } else {
    .noise_commit(ns, pd$t + pd$h, pd$dW)
  }
  invisible()
}

# Wiener increments over [t, t + h], reusing stored future increments.
.noise_take <- function(ns, t, h) {
  .noise_settle(ns, t)
  need <- h
  dW <- numeric(ns$m)
  dZ <- if (ns$needs_dZ) numeric(ns$m) else NULL
  tol <- 1e-10 * abs(h)
  while (need > tol && length(ns$stack)) {
    e <- ns$stack[[1L]]
    if (e$h <= need + tol) {
      dW <- dW + e$dW
      if (ns$needs_dZ) dZ <- dZ + e$dZ
      need <- need - e$h
      ns$stack <- ns$stack[-1L]
    } else {
      sp <- .noise_split(ns, e, need)
      dW <- dW + sp$first$dW
      if (ns$needs_dZ) dZ <- dZ + sp$first$dZ
      ns$stack[[1L]] <- sp$rest
      need <- 0
    }
  }
  if (need > tol) {
    dW <- dW + stats::rnorm(ns$m, 0, sqrt(need))
    if (ns$needs_dZ) dZ <- dZ + stats::rnorm(ns$m, 0, sqrt(need))
  }
  ns$pending <- list(t = t, h = h, dW = dW, dZ = dZ)
  list(dW = dW, dZ = dZ)
}

# ---------------------------------------------------------------------------
# Algorithms

.new_sde_alg <- function(name, order, adaptive, perform, needs_dZ = FALSE, delta = 1,
                         noise = c("diagonal", "scalar", "general"), interpretation = "Ito", ...) {
  structure(list(name = name, order = order, adaptive = adaptive, perform = perform,
                 needs_dZ = needs_dZ, delta = delta, noise = noise,
                 interpretation = interpretation, ...),
            class = c(name, "SDEAlgorithm"))
}

#' @export
print.SDEAlgorithm <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

# g(u) dW for the three noise shapes.
.gdw <- function(G, dW, type) if (type == "general") as.numeric(G %*% dW) else G * dW

.rms <- function(x) sqrt(mean(x^2))

#' SDE solver algorithms
#'
#' Algorithms for [solve()] on SDE problems, named as in StochasticDiffEq.jl.
#' `EM()` is the Euler-Maruyama method (strong order 0.5) and `EulerHeun()`
#' its Stratonovich counterpart; both are fixed step methods (pass `dt`).
#' `LambaEM()` is Euler-Maruyama with Lamba's error estimate, adaptive for
#' any noise. `RKMil()` is the Runge-Kutta Milstein method (strong order 1)
#' and `SRIW1()` and `SOSRI()` are Rößler's strong order 1.5 SRI methods
#' (SOSRI is stability optimized); these three take diagonal or scalar noise
#' and adapt their step size by default (pass `adaptive = FALSE` and `dt` for
#' fixed steps). Adaptive methods use rejection sampling with memory, so a
#' rejected step does not bias the Brownian path. `LambaEulerHeun()` is the
#' Stratonovich counterpart of `LambaEM()`. `RKMilCommute()` is the
#' Runge-Kutta Milstein method for commutative noise (strong order 1 for any
#' noise shape with commuting diffusion columns). `SRA1()`, `SRA2()`,
#' `SRA3()`, `SOSRA()` and `SOSRA2()` are Rößler's strong order 1.5 methods
#' for additive noise (diffusion independent of the state; SOSRA and SOSRA2
#' are stability optimized). Without `alg_hints`, `SOSRI()` is the default
#' for diagonal noise and `LambaEM()` for other noise; `alg_hints =
#' "additive"` selects `SOSRA()`, `"commutative"` `RKMilCommute()` and
#' `"stratonovich"` `RKMil(interpretation = "Stratonovich")` or
#' `LambaEulerHeun()`, as in StochasticDiffEq.
#'
#' Solver options of `solve()` for SDE problems: `dt`, `adaptive`,
#' `abstol` and `reltol` (default `1e-2`), `delta` (weight of the drift error
#' estimate), `seed` (sets the random number generator for this solve and
#' restores the previous state afterwards), `saveat`, `tstops`, `callback`,
#' `save_everystep`, `dtmax`, `maxiters`. The solution's Brownian path is
#' `sol$W`, a list with the step times `t` and the matrix `W` of cumulative
#' increments. As in StochasticDiffEq, solutions are interpolated linearly
#' between steps and events inside a step are located on that interpolant;
#' the Brownian path then continues from a bridge sample at the event time,
#' so the state at an event carries the interpolation error of its step.
#'
#' @param split For `EM()`: evaluate the diffusion at the Euler predictor
#'   (`TRUE`, the default) rather than at the previous state. For
#'   `LambaEM()`: evaluate it at the previous state (`TRUE`) rather than at
#'   the predictor, as StochasticDiffEq's out of place cache does.
#' @param interpretation `"Ito"` or `"Stratonovich"` for `RKMil()` and
#'   `RKMilCommute()`.
#' @return An algorithm object.
#' @examples
#' prob <- SDEProblem(function(u, p, t) p[1] * u, function(u, p, t) p[2] * u,
#'                    1, c(0, 1), c(1.5, 0.3))
#' sol <- solve(prob, SRIW1(), seed = 2)
#' sol$retcode
#' @name sde_algorithms
NULL

#' @rdname sde_algorithms
#' @export
EM <- function(split = TRUE) {
  .new_sde_alg("EM", 0.5, FALSE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    K <- u + dt * f(u, p, t)
    uc <- if (split) K else u
    list(u = K + .gdw(g(uc, p, t), dW, type), nf = 1L)
  })
}

#' @rdname sde_algorithms
#' @export
EulerHeun <- function() {
  .new_sde_alg("EulerHeun", 0.5, FALSE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    ftmp <- f(u, p, t)
    gtmp <- g(u, p, t)
    tmp <- u + dt * ftmp + .gdw(gtmp, dW, type)
    gtmp2 <- (gtmp + g(tmp, p, t + dt)) / 2
    list(u = u + (dt / 2) * (ftmp + f(tmp, p, t + dt)) + .gdw(gtmp2, dW, type), nf = 2L)
  }, interpretation = "Stratonovich")
}

#' @rdname sde_algorithms
#' @export
LambaEM <- function(split = TRUE) {
  .new_sde_alg("LambaEM", 0.5, TRUE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    # LambaEMConstantCache, the out of place version of lamba.jl
    du1 <- f(u, p, t)
    K <- u + dt * du1
    L <- if (split) g(u, p, t + dt) else g(K, p, t + dt)
    noise <- .gdw(L, dW, type)
    unew <- K + noise
    sq <- sqrt(abs(dt))
    du2 <- f(K, p, t + dt)
    Ed <- dt * (du2 - du1) / 2
    ggprime <- (g(K + noise * sq, p, t) - L) / sq
    En <- if (type == "general") as.numeric(ggprime %*% (dW^2 - dt)) / 2 else ggprime * (dW^2 - dt) / 2
    list(u = unew, err = delta * Ed + En, nf = 2L)
  })
}

#' @rdname sde_algorithms
#' @export
RKMil <- function(interpretation = c("Ito", "Stratonovich")) {
  interpretation <- match.arg(interpretation)
  .new_sde_alg("RKMil", 1, TRUE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    sq <- sqrt(abs(dt))
    du1 <- f(u, p, t)
    K <- u + dt * du1
    L <- g(u, p, t)
    if (interpretation == "Ito") {
      ggprime <- (g(K + L * sq, p, t) - L) / sq
      mil <- ggprime * (dW^2 - abs(dt)) / 2
    } else {
      ggprime <- (g(u + L * sq, p, t) - L) / sq
      mil <- ggprime * dW^2 / 2
    }
    unew <- K + L * dW + mil
    du2 <- f(K, p, t + dt)
    Ed <- dt * (du2 - du1) / 2
    En <- dW^3 * ((du2 - L) / sq)^2 / 6
    list(u = unew, err = delta * Ed + En, nf = 2L)
  }, noise = c("diagonal", "scalar"), interpretation = interpretation)
}

#' @rdname sde_algorithms
#' @export
SRIW1 <- function() {
  .new_sde_alg("SRIW1", 1.5, TRUE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    sq <- sqrt(abs(dt))
    chi1 <- (dW^2 - abs(dt)) / (2 * sq)
    chi2 <- (dW + dZ / sqrt(3)) / 2
    chi3 <- (dW^3 - 3 * dW * dt) / (6 * dt)
    fH01 <- dt * f(u, p, t)
    g1 <- g(u, p, t)
    fH01o4 <- fH01 / 4
    dto4 <- dt / 4
    g1o2 <- g1 / 2
    H0 <- u + 3 * (fH01o4 + chi2 * g1o2)
    H11 <- u + fH01o4 + sq * g1o2
    H12 <- u + fH01 - sq * g1
    g2 <- g(H11, p, t + dto4)
    g3 <- g(H12, p, t + dt)
    H13 <- u + fH01o4 + sq * (-5 * g1 + 3 * g2 + g3 / 2)
    g4 <- g(H13, p, t + dto4)
    fH02 <- dt * f(H0, p, t + 3 * dto4)
    g2o3 <- g2 / 3
    Fg2o3 <- 4 * g2o3
    g3o3 <- g3 / 3
    Tg3o3 <- 2 * g3o3
    mg1 <- -g1
    E2 <- chi2 * (2 * g1 - Fg2o3 - Tg3o3) + chi3 * (2 * mg1 + 5 * g2o3 - Tg3o3 + g4)
    unew <- u + (fH01 + 2 * fH02) / 3 + dW * (mg1 + Fg2o3 + Tg3o3) + chi1 * (mg1 + Fg2o3 - g3o3) + E2
    list(u = unew, err = delta * (fH01 + fH02) + E2, nf = 2L)
  }, needs_dZ = TRUE, delta = 1 / 6, noise = c("diagonal", "scalar"))
}

# SOSRI tableau, SOSRIConstantCache in caches/rossler_caches.jl.
.SOSRI <- list(
  a021 = -0.04199224421316468, a031 = 2.842612915017106, a032 = -2.0527723684000727,
  a041 = 4.338237071435815, a042 = -2.8895936137439793, a043 = 2.3017575594644466,
  a121 = 0.26204282091330466, a131 = 0.20903646383505375, a132 = -0.1502377115150361,
  a141 = 0.05836595312746999, a142 = 0.6149440396332373, a143 = 0.08535117634046772,
  b021 = -0.21641093549612528, b031 = 1.5336352863679572, b032 = 0.26066223492647056,
  b041 = -1.0536037558179159, b042 = 1.7015284721089472, b043 = -0.20725685784180017,
  b121 = -0.5119011827621657, b131 = 2.67767339866713, b132 = -4.9395031322250995,
  b141 = 0.15580956238299215, b142 = 3.2361551006624674, b143 = -1.4223118283355949,
  a1 = 1.140099274172029, a2 = -0.6401334255743456, a3 = 0.4736296532772559, a4 = 0.026404498125060714,
  c02 = -0.04199224421316468, c03 = 0.7898405466170333, c04 = 3.7504010171562823,
  c11 = 0.0, c12 = 0.26204282091330466, c13 = 0.05879875232001766, c14 = 0.758661169101175,
  beta11 = -1.8453464565104432, beta12 = 2.688764531100726, beta13 = -0.2523866501071323,
  beta14 = 0.40896857551684956, beta21 = 0.4969658141589478, beta22 = -0.5771202869753592,
  beta23 = -0.12919702470322217, beta24 = 0.2093514975196336, beta31 = 2.8453464565104425,
  beta32 = -2.688764531100725, beta33 = 0.2523866501071322, beta34 = -0.40896857551684945,
  beta41 = 0.11522663875443433, beta42 = -0.57877086147738, beta43 = 0.2857851028163886,
  beta44 = 0.17775911990655704)

#' @rdname sde_algorithms
#' @export
SOSRI <- function() {
  .new_sde_alg("SOSRI", 1.5, TRUE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    S <- .SOSRI
    sq <- sqrt(abs(dt))
    chi1 <- (dW^2 - abs(dt)) / (2 * sq)
    chi2 <- (dW + dZ / sqrt(3)) / 2
    chi3 <- (dW^3 - 3 * dW * dt) / (6 * dt)
    k1 <- f(u, p, t)
    g1 <- g(u, p, t + S$c11 * dt)
    H01 <- u + dt * S$a021 * k1 + S$b021 * chi2 * g1
    H11 <- u + dt * S$a121 * k1 + sq * S$b121 * g1
    k2 <- f(H01, p, t + S$c02 * dt)
    g2 <- g(H11, p, t + S$c12 * dt)
    H02 <- u + dt * (S$a031 * k1 + S$a032 * k2) + chi2 * (S$b031 * g1 + S$b032 * g2)
    H12 <- u + dt * (S$a131 * k1 + S$a132 * k2) + sq * (S$b131 * g1 + S$b132 * g2)
    k3 <- f(H02, p, t + S$c03 * dt)
    g3 <- g(H12, p, t + S$c13 * dt)
    H03 <- u + dt * (S$a041 * k1 + S$a042 * k2 + S$a043 * k3) +
      chi2 * (S$b041 * g1 + S$b042 * g2 + S$b043 * g3)
    H13 <- u + dt * (S$a141 * k1 + S$a142 * k2 + S$a143 * k3) +
      sq * (S$b141 * g1 + S$b142 * g2 + S$b143 * g3)
    k4 <- f(H03, p, t + S$c04 * dt)
    g4 <- g(H13, p, t + S$c14 * dt)
    E2 <- chi2 * (S$beta31 * g1 + S$beta32 * g2 + S$beta33 * g3 + S$beta34 * g4) +
      chi3 * (S$beta41 * g1 + S$beta42 * g2 + S$beta43 * g3 + S$beta44 * g4)
    unew <- u + dt * (S$a1 * k1 + S$a2 * k2 + S$a3 * k3 + S$a4 * k4) + E2 +
      dW * (S$beta11 * g1 + S$beta12 * g2 + S$beta13 * g3 + S$beta14 * g4) +
      chi1 * (S$beta21 * g1 + S$beta22 * g2 + S$beta23 * g3 + S$beta24 * g4)
    E1 <- dt * (k1 + k2 + k3 + k4)
    list(u = unew, err = delta * E1 + E2, nf = 4L)
  }, needs_dZ = TRUE, noise = c("diagonal", "scalar"))
}

#' @rdname sde_algorithms
#' @export
LambaEulerHeun <- function() {
  # LambaEulerHeunCache of lamba.jl
  .new_sde_alg("LambaEulerHeun", 0.5, TRUE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    du1 <- f(u, p, t)
    K <- u + dt * du1
    L <- g(u, p, t)
    tmp <- K + .gdw(L, dW, type)
    unew <- u + (dt / 2) * (du1 + f(tmp, p, t + dt)) + .gdw((L + g(tmp, p, t + dt)) / 2, dW, type)
    sq <- sqrt(abs(dt))
    if (type == "general") {
      gs <- sqrt(sum(L^2))
      En <- rep((sqrt(sum(g(u + sq * gs, p, t)^2)) - gs) / (2 * sq) * .rms(dW^2), length(u))
    } else {
      En <- (g(u + sq * L, p, t) - L) / (2 * sq) * dW^2
    }
    Ed <- dt * (f(K, p, t + dt) - du1) / 2
    list(u = unew, err = delta * Ed + En, nf = 3L)
  }, interpretation = "Stratonovich")
}

# Strong order 1.5 methods for additive noise (sra.jl, caches/sra_caches.jl):
# SRA1, SRA2 and the three stage tableaus of SRA3, SOSRA and SOSRA2.
.SRA3 <- list(a21 = 1, a31 = 1 / 4, a32 = 1 / 4, b21 = 0, b31 = 1, b32 = 1 / 2, c02 = 1, c03 = 1 / 2,
              c11 = 1, c12 = 0, c13 = 0, a1 = 1 / 6, a2 = 1 / 6, a3 = 2 / 3, beta11 = 1, beta12 = 0,
              beta13 = 0, beta21 = -1, beta22 = 1, beta23 = 0)
.SOSRA <- list(a1 = 0.2889874966892885, a2 = 0.6859880440839937, a3 = 0.025024459226717772,
               c02 = 0.6923962376159507, c03 = 1, c11 = 0, c12 = 0.041248171110700504, c13 = 1,
               beta11 = -16.792534242221663, beta12 = 17.514995785380226, beta13 = 0.27753845684143835,
               beta21 = 0.4237535769069274, beta22 = 0.6010381474428539, beta23 = -1.0247917243497813,
               a21 = 0.6923962376159507, a31 = -3.1609142252828395, a32 = 4.1609142252828395,
               b21 = 1.3371632704399763, b31 = 1.442371048468624, b32 = 1.8632741501139225)
.SOSRA2 <- list(a1 = 0.4999999999999998, a2 = -0.9683897375354181, a3 = 1.4683897375354185,
                c02 = 1, c03 = 1, c11 = 0, c12 = 1, c13 = 1,
                beta11 = 0, beta12 = 0.92438032145683, beta13 = 0.07561967854316998,
                beta21 = 1, beta22 = -0.8169981105823436, beta23 = -0.18300188941765633,
                a21 = 1, a31 = 0.9511849235504364, a32 = 0.04881507644956362,
                b21 = 0.7686101171003622, b31 = 0.43886792994934987, b32 = 0.7490415909204886)

# ThreeStageSRACache step; chi2 = I_(1,0) / h.
.sra3_step <- function(S) function(f, g, u, p, t, dt, dW, dZ, type, delta) {
  chi2 <- (dW + dZ / sqrt(3)) / 2
  g1 <- g(u, p, t + S$c11 * dt)
  k1 <- f(u, p, t)
  H01 <- u + dt * S$a21 * k1 + S$b21 * .gdw(g1, chi2, type)
  g2 <- g(H01, p, t + S$c12 * dt)
  k2 <- f(H01, p, t + S$c02 * dt)
  H02 <- u + dt * (S$a31 * k1 + S$a32 * k2) + .gdw(S$b31 * g1 + S$b32 * g2, chi2, type)
  g3 <- g(H02, p, t + S$c13 * dt)
  k3 <- f(H02, p, t + S$c03 * dt)
  E2 <- .gdw(S$beta21 * g1 + S$beta22 * g2 + S$beta23 * g3, chi2, type)
  unew <- u + dt * (S$a1 * k1 + S$a2 * k2 + S$a3 * k3) + E2 +
    .gdw(S$beta11 * g1 + S$beta12 * g2 + S$beta13 * g3, dW, type)
  list(u = unew, err = delta * dt * (k1 + k2 + k3) + E2, nf = 3L)
}

#' @rdname sde_algorithms
#' @export
SRA1 <- function() {
  .new_sde_alg("SRA1", 1.5, TRUE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    gpdt <- g(u, p, t + dt)
    chi2 <- (dW + dZ / sqrt(3)) / 2
    k1 <- dt * f(u, p, t)
    k2 <- dt * f(u + 3 * k1 / 4 + 3 * .gdw(gpdt, chi2, type) / 2, p, t + 3 * dt / 4)
    E2 <- .gdw(g(u, p, t) - gpdt, chi2, type)
    unew <- u + k1 / 3 + 2 * k2 / 3 + E2 + .gdw(gpdt, dW, type)
    list(u = unew, err = delta * (k1 + k2) + E2, nf = 2L)
  }, needs_dZ = TRUE)
}

#' @rdname sde_algorithms
#' @export
SRA2 <- function() {
  .new_sde_alg("SRA2", 1.5, TRUE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    chi2 <- (dW + dZ / sqrt(3)) / 2
    g1 <- g(u, p, t + dt / 3)
    k1 <- f(u, p, t)
    H01 <- u + dt * 3 / 4 * k1 + 3 / 2 * .gdw(g1, chi2, type)
    g2 <- g(H01, p, t + dt)
    k2 <- f(H01, p, t + 3 / 4 * dt)
    E2 <- .gdw(3 / 2 * g1 - 3 / 2 * g2, chi2, type)
    unew <- u + dt * (k1 / 3 + 2 * k2 / 3) + E2 + .gdw(g2, dW, type)
    list(u = unew, err = delta * dt * (k1 + k2) + E2, nf = 2L)
  }, needs_dZ = TRUE)
}

#' @rdname sde_algorithms
#' @export
SRA3 <- function() .new_sde_alg("SRA3", 1.5, TRUE, .sra3_step(.SRA3), needs_dZ = TRUE)

#' @rdname sde_algorithms
#' @export
SOSRA <- function() .new_sde_alg("SOSRA", 1.5, TRUE, .sra3_step(.SOSRA), needs_dZ = TRUE)

#' @rdname sde_algorithms
#' @export
SOSRA2 <- function() .new_sde_alg("SOSRA2", 1.5, TRUE, .sra3_step(.SOSRA2), needs_dZ = TRUE)

#' @rdname sde_algorithms
#' @export
RKMilCommute <- function(interpretation = c("Ito", "Stratonovich")) {
  interpretation <- match.arg(interpretation)
  # RKMilCommuteCache of low_order.jl with the commutative iterated
  # integrals J = dW dW' / 2 (IICommutative)
  .new_sde_alg("RKMilCommute", 1, TRUE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    sq <- sqrt(abs(dt))
    du1 <- f(u, p, t)
    K <- u + dt * du1
    L <- g(u, p, t)
    ito <- interpretation == "Ito"
    if (type == "general") {
      J <- outer(dW, dW) / 2
      if (ito) J <- J - diag(abs(dt) / 2, length(dW))
      mil <- 0 * u
      gnorm <- 0
      for (j in seq_along(dW)) {
        Dgj <- (g(K + sq * L[, j], p, t) - L) / sq
        gnorm <- gnorm + .rms(Dgj)
        mil <- mil + as.numeric(Dgj %*% J[, j])
      }
      unew <- K + as.numeric(L %*% dW) + mil
    } else {
      J <- dW^2 / 2
      if (ito) J <- J - abs(dt) / 2
      Dgj <- (g((if (ito) K else u) + sq * L, p, t) - L) / sq
      gnorm <- .rms(Dgj)
      unew <- K + L * dW + Dgj * J
    }
    En <- .rms(dW)^3 * gnorm^2 / 6
    Ed <- dt * (f(K, p, t + dt) - du1) / 2
    list(u = unew, err = rep(.rms(delta * Ed) + En, length(u)), nf = 2L)
  }, interpretation = interpretation)
}

# ---------------------------------------------------------------------------
# Implicit methods (perform_step/sdirk.jl, the in place ImplicitEMCache,
# ImplicitEulerHeunCache and ImplicitRKMilCache steps, alg_utils.jl
# alg_mass_matrix_compatible). The drift is implicit with the theta method
# (or the implicit midpoint rule with `symplectic`), solved by the modified
# Newton method of the ODE solvers (.nlsolve). The stage equation
# M z = dt f(tmp + theta z, t + c dt) of StochasticDiffEq is solved for the
# stage state U = tmp + theta z. With a singular mass matrix the algebraic
# equations then hold at the end of each step.

.sde_implicit_alg <- function(name, kind, theta, symplectic, order, interpretation = "Ito",
                              noise = c("diagonal", "scalar", "general"), split = FALSE) {
  if (symplectic) theta <- 1 / 2
  alg <- .new_sde_alg(name, order, TRUE, perform = NULL, noise = noise, interpretation = interpretation,
                      theta = theta, symplectic = symplectic, kind = kind, split = split)
  alg$mass_matrix_ok <- symplectic || theta == 1
  alg$init <- function(ctx, t, u) {
    cache <- new.env(parent = emptyenv())
    cache$nl <- .new_nlsolver()
    cache
  }
  alg$step <- function(ctx, g, u, t, dt, dW, type, delta, cache, ...)
    .sde_implicit_step(alg, ctx, g, u, t, dt, dW, type, delta, cache)
  alg
}

.sde_implicit_step <- function(alg, ctx, g, u, t, dt, dW, type, delta, cache) {
  p <- ctx$p
  theta <- alg$theta
  sq <- sqrt(abs(dt))
  f0 <- ctx$f(u, p, t)
  L <- g(u, p, t)
  gdw <- .gdw(L, dW, type)
  nf <- 1L
  split <- isTRUE(alg$split)
  # the split step methods (implicit_split_step.jl) add the noise after the
  # implicit drift step
  noise_after <- function(ud) {
    if (alg$kind == "EulerHeun") ud + .gdw((g(ud + gdw, p, t) + L) / 2, dW, type) else ud + gdw
  }
  if (split) {
    # no noise inside the implicit stage
  } else if (alg$kind == "EulerHeun") {
    L <- (g(u + gdw, p, t) + L) / 2
    gdw <- .gdw(L, dW, type)
  } else if (alg$kind == "RKMil") {
    if (alg$interpretation == "Ito") {
      ggprime <- (g(u + dt * f0 + sq * L, p, t) - L) / sq
      gdw <- gdw + ggprime * (dW^2 - abs(dt)) / 2
    } else {
      ggprime <- (g(u + sq * L, p, t) - L) / sq
      gdw <- gdw + ggprime * dW^2 / 2
    }
  }
  gstage <- if (split) 0 * u else gdw
  tmp <- if (alg$symplectic) u + gstage / 2 else u + dt * (1 - theta) * f0 + gstage
  # theta = 0 is the explicit limit, u = tmp, with no stage to solve
  if (theta == 0) return(list(u = if (split) noise_after(tmp) else tmp, nf = nf, cache = cache, k = NULL))
  z0 <- if (alg$symplectic) 0 * u else dt * f0
  cstage <- if (alg$symplectic) 1 / 2 else theta
  nl <- .nlsolve(cache$nl, ctx, t, dt, tmp + theta * z0, .mass_mul(ctx$mass, tmp) / (theta * dt),
                 theta, 1, t + cstage * dt, u)
  out <- list(nf = nf + nl$nf, njacs = nl$njacs, nsolve = nl$nsolve, nw = nl$nw, cache = cache, k = NULL)
  if (!nl$ok) return(c(out, list(fail = TRUE, u = u)))
  unew <- if (alg$symplectic) 2 * nl$z - u else nl$z
  if (split) unew <- noise_after(unew)
  if (!isTRUE(ctx$adaptive)) return(c(out, list(u = unew)))
  J <- cache$nl$J %||% .calc_J(ctx, u, t, f0)$J
  Ed <- dt * dt * as.numeric(J %*% f0) / 2
  if (alg$kind == "RKMil" && !split) {
    En <- abs(dW^3) * abs(ggprime)^2 / 6
  } else {
    dw2 <- if (alg$kind == "EM") dW^2 - dt else dW^2
    shift <- if (alg$kind == "EM") u + dt * f0 else u
    if (type == "general") {
      gs <- sqrt(sum(L^2))
      gs2 <- sqrt(sum(g(shift + sq * gs, p, t)^2))
      En <- rep((gs2 - gs) / (2 * sq) * .rms(dw2), length(u))
    } else {
      En <- (g(shift + sq * L, p, t) - L) / (2 * sq) * dw2
    }
  }
  c(out, list(u = unew, err = delta * Ed + En))
}

# SKenCarp (perform_step/kencarp.jl, SKenCarpTableau of tableaus.jl): the
# ESDIRK method of KenCarp3 for the drift with the strong order 1.5 noise
# terms of additive noise. Each stage U = tmp + gamma z with
# M z = dt f(U, t + c dt) is solved by .nlsolve for U.
.SKC <- local({
  gam <- 1767732205903 / 4055673282236
  list(gamma = gam, c2 = 2 * gam, c3 = 3 / 5,
       a31 = 2746238789719 / 10658868560708, a32 = -640167445237 / 6845629431997,
       a41 = 1471266399579 / 7840856788654, a42 = -4482444167858 / 7529755066697,
       a43 = 11266239266428 / 11593286722821,
       btilde = c(681815649026867975666107 / 25159934323302256049469295,
                  18411887981491912264464127 / 167175311446532472108584143,
                  -12719313754959329011138489 / 123410692144842870217698057,
                  -47289384293135913063989 / 1383962894467812063558225),
       nb021 = -12.246764387585055918338744103409192607986567514699471403397969732723452087723101,
       nb043 = -14.432096958608752822047165680776748797565142459789556194474191884258734697161106)
})

.skencarp_step <- function(alg, ctx, g, u, t, dt, dW, type, delta, cache, dZ) {
  S <- .SKC
  gam <- S$gamma
  p <- ctx$p
  chi2 <- (dW + dZ / sqrt(3)) / 2
  out <- list(nf = 1L, njacs = 0L, nsolve = 0L, nw = 0L, cache = cache, k = NULL)
  # one stage: U = tmp + gamma z, M z = dt f(U, tc); extrapolant :min_correct
  stage <- function(tmp, tc) {
    nl <- .nlsolve(cache$nl, ctx, t, dt, tmp, .mass_mul(ctx$mass, tmp) / (gam * dt), gam, 1, tc, u)
    for (s in c("nf", "njacs", "nsolve", "nw")) out[[s]] <<- out[[s]] + as.integer(nl[[s]] %||% 0L)
    if (!nl$ok) return(NULL)
    (nl$z - tmp) / gam
  }
  z1 <- dt * ctx$f(u, p, t)
  g1 <- g(u, p, t)
  z2 <- stage(u + gam * z1 + S$nb021 * .gdw(g1, chi2, type), t + S$c2 * dt)
  if (is.null(z2)) return(c(out, list(fail = TRUE, u = u)))
  z3 <- stage(u + S$a31 * z1 + S$a32 * z2, t + S$c3 * dt)
  if (is.null(z3)) return(c(out, list(fail = TRUE, u = u)))
  # g3 = g1 since the noise stage times c13 = 0
  z4 <- stage(u + S$a41 * z1 + S$a42 * z2 + S$a43 * z3 + S$nb043 * .gdw(g1, chi2, type), t + dt)
  if (is.null(z4)) return(c(out, list(fail = TRUE, u = u)))
  g4 <- g(u, p, t + dt)
  E2 <- .gdw(g1 - g4, chi2, type)
  unew <- u + S$a41 * z1 + S$a42 * z2 + S$a43 * z3 + gam * z4 + .gdw(g4, dW, type) + E2
  if (!isTRUE(ctx$adaptive)) return(c(out, list(u = unew)))
  tmp <- S$btilde[1L] * z1 + S$btilde[2L] * z2 + S$btilde[3L] * z3 + S$btilde[4L] * z4 + E2
  # smooth_est (Shampine): (M - gamma dt J)^-1 tmp; W = J - M / W_gdt
  E1 <- if (!is.null(cache$nl$W)) -.wsolve(cache$nl$W, tmp) / cache$nl$W_gdt else tmp
  out$nsolve <- out$nsolve + 1L
  c(out, list(u = unew, err = delta * E1 + E2))
}

#' Implicit SDE solvers
#'
#' Drift implicit methods of StochasticDiffEq.jl for stiff SDEs and for SDEs
#' with algebraic equations (a singular mass matrix, as `mtkcompile()` leaves
#' when an algebraic equation cannot be solved explicitly). `ImplicitEM()` is
#' the drift implicit Euler-Maruyama method (theta method; `theta = 1/2` is
#' `STrapezoid()` and `symplectic = TRUE` the implicit midpoint rule
#' `SImplicitMidpoint()`), `ImplicitEulerHeun()` its Stratonovich
#' counterpart, and `ImplicitRKMil()` the drift implicit Runge-Kutta Milstein
#' method (strong order 1, diagonal or scalar noise). `ISSEM()` and
#' `ISSEulerHeun()` are the split step variants: the drift step is implicit
#' and the noise is added after it. `SKenCarp()` is the stochastic ESDIRK
#' method of Kennedy and Carpenter for additive noise (strong order 1.5,
#' order 2 for the drift). All adapt their step size (or take `adaptive =
#' FALSE` and `dt`). As in StochasticDiffEq, mass matrices need `theta = 1`
#' or `symplectic = TRUE`; explicit methods do not accept them. The default
#' algorithm follows StochasticDiffEq's `default_algorithm`, with the
#' `alg_hints` option of `solve()` (`"stiff"`, `"stratonovich"`,
#' `"commutative"`, `"additive"`): with a mass matrix or `"stiff"` it is
#' `ImplicitRKMil()` for diagonal noise, `ISSEM()` (or `ImplicitEulerHeun()`
#' for Stratonovich) for other noise and `SKenCarp()` for additive noise.
#' Events inside a step are located on the linear interpolant of the step,
#' as StochasticDiffEq does for every method.
#'
#' @param theta Implicitness of the drift.
#' @param symplectic Use the implicit midpoint rule.
#' @param interpretation `"Ito"` or `"Stratonovich"` for `ImplicitRKMil()`.
#' @return An algorithm object.
#' @examples
#' f <- function(u, p, t) c(-u[1] + u[2], u[2] - 2 * u[1])
#' g <- function(u, p, t) c(0.1 * u[1], 0)
#' prob <- SDEProblem(SDEFunction(f, g, mass_matrix = diag(c(1, 0))), c(1, 2), c(0, 1))
#' sol <- solve(prob, ImplicitEM(), dt = 0.01, adaptive = FALSE, seed = 1)
#' tail(sol$u, 1)
#' @name sde_implicit_algorithms
NULL

#' @rdname sde_implicit_algorithms
#' @export
ImplicitEM <- function(theta = 1, symplectic = FALSE)
  .sde_implicit_alg("ImplicitEM", "EM", theta, symplectic, 0.5)

#' @rdname sde_implicit_algorithms
#' @export
STrapezoid <- function() ImplicitEM(theta = 1 / 2)

#' @rdname sde_implicit_algorithms
#' @export
SImplicitMidpoint <- function() ImplicitEM(theta = 1 / 2, symplectic = TRUE)

#' @rdname sde_implicit_algorithms
#' @export
ImplicitEulerHeun <- function(theta = 1, symplectic = FALSE)
  .sde_implicit_alg("ImplicitEulerHeun", "EulerHeun", theta, symplectic, 0.5, interpretation = "Stratonovich")

#' @rdname sde_implicit_algorithms
#' @export
SKenCarp <- function() {
  alg <- .sde_implicit_alg("SKenCarp", "SKenCarp", 1, FALSE, 1.5)
  alg$order <- 2
  alg$needs_dZ <- TRUE
  alg$step <- function(ctx, g, u, t, dt, dW, type, delta, cache, dZ = NULL, ...)
    .skencarp_step(alg, ctx, g, u, t, dt, dW, type, delta, cache, dZ)
  alg
}

#' @rdname sde_implicit_algorithms
#' @export
ISSEM <- function(theta = 1, symplectic = FALSE)
  .sde_implicit_alg("ISSEM", "EM", theta, symplectic, 0.5, split = TRUE)

#' @rdname sde_implicit_algorithms
#' @export
ISSEulerHeun <- function(theta = 1, symplectic = FALSE)
  .sde_implicit_alg("ISSEulerHeun", "EulerHeun", theta, symplectic, 0.5, interpretation = "Stratonovich",
                    split = TRUE)

#' @rdname sde_implicit_algorithms
#' @export
ImplicitRKMil <- function(theta = 1, symplectic = FALSE, interpretation = c("Ito", "Stratonovich")) {
  interpretation <- match.arg(interpretation)
  .sde_implicit_alg("ImplicitRKMil", "RKMil", theta, symplectic, 1, interpretation = interpretation,
                    noise = c("diagonal", "scalar"))
}

# ---------------------------------------------------------------------------
# Driver

# default_sde_alg.jl: `alg_hints` among "stiff", "stratonovich",
# "commutative" and "additive".
.sde_default_alg <- function(prob, alg_hints = NULL) {
  diag_noise <- (prob$noise_type %||% prob$f$noise_type %||% "diagonal") == "diagonal"
  stiff <- "stiff" %in% alg_hints || !.is_identity_mass(prob$f$mass_matrix)
  strat <- "stratonovich" %in% alg_hints
  alg <- SOSRI()
  if ("commutative" %in% alg_hints) alg <- RKMilCommute()
  if (stiff) alg <- ImplicitRKMil()
  if (strat) alg <- if (stiff) ImplicitRKMil(interpretation = "Stratonovich") else RKMil(interpretation = "Stratonovich")
  if (!diag_noise) {
    alg <- if (strat) { if (stiff) ImplicitEulerHeun() else LambaEulerHeun() }
           else if (stiff) ISSEM() else LambaEM()
  }
  if ("additive" %in% alg_hints) alg <- if (stiff) SKenCarp() else SOSRA()
  alg
}

.sde_linear_interp <- function(theta, dt, y0, y1, k) (1 - theta) * y0 + theta * y1

# Run the core driver with the SDE algorithm `alg` on drift f and diffusion g.
.sde_integrate <- function(prob, alg, opts, f = prob$f$f, g = prob$f$g, u0 = prob$u0,
                           type = prob$noise_type %||% prob$f$noise_type %||% "diagonal",
                           m = NULL) {
  setup <- .sde_setup(prob, alg, opts, f = f, g = g, u0 = u0, type = type, m = m)
  setup$post(.integrate(setup$prob, setup$alg, setup$opts))
}

# The core driver problem, algorithm and options of an SDE solve, and the
# post processing of its solution (the Brownian path); init() uses it too.
.sde_setup <- function(prob, alg, opts, f = prob$f$f, g = prob$f$g, u0 = prob$u0,
                       type = prob$noise_type %||% prob$f$noise_type %||% "diagonal",
                       m = NULL) {
  if (!inherits(alg, "SDEAlgorithm")) stop("not an SDE algorithm: ", class(alg)[1L], call. = FALSE)
  if (!(type %in% alg$noise))
    .stopf("%s() requires diagonal or scalar noise; use EM() or LambaEM() for general noise", alg$name)
  M <- prob$f$mass_matrix
  mm <- !.is_identity_mass(M)
  if (mm) {
    # alg_mass_matrix_compatible (alg_utils.jl)
    if (is.null(alg$mass_matrix_ok)) stop("This solver is not able to use mass matrices.", call. = FALSE)
    if (!alg$mass_matrix_ok) stop("Algorithm must be set as symplectic or theta=1 for mass matrices", call. = FALSE)
  }
  if (prob$tspan[2L] <= prob$tspan[1L]) stop("SDE problems need tspan[2] > tspan[1]", call. = FALSE)
  m <- m %||% switch(type, diagonal = length(u0), scalar = 1L,
                     general = prob$f$nbrownians %||% ncol(as.matrix(g(u0, prob$p, prob$tspan[1L]))))
  ns <- .noise_new(prob$tspan[1L], m, alg$needs_dZ)
  delta <- opts$delta %||% alg$delta
  adaptive <- alg$adaptive && (opts$adaptive %||% TRUE)
  runner <- .new_alg(alg$name, alg$order, adaptive = adaptive,
    step = function(ctx, t, u, dt, cache) {
      inc <- .noise_take(ns, t, dt)
      if (!is.null(alg$step)) return(alg$step(ctx, g, u, t, dt, inc$dW, type, delta, cache, dZ = inc$dZ))
      res <- alg$perform(ctx$f, g, u, ctx$p, t, dt, inc$dW, inc$dZ, type, delta)
      list(u = res$u, err = res$err, k = NULL, cache = cache, nf = res$nf)
    },
    interp = .sde_linear_interp, init = alg$init %||% function(ctx, t, u) list(),
    mass_matrix = mm, on_reject = if (!is.null(alg$step)) .nl_on_reject)
  opts$abstol <- opts$abstol %||% 1e-2
  opts$reltol <- opts$reltol %||% 1e-2
  if (adaptive) {
    opts$qmax <- opts$qmax %||% (9 / 8)
    opts$qmin <- opts$qmin %||% (1 / 5)
  }
  opts$delta <- NULL
  dprob <- prob
  dprob$u0 <- u0
  dprob$f <- list(f = f, mass_matrix = if (mm) M)
  # a passive callback records where each accepted step ends, after events
  t_last <- new.env(parent = emptyenv())
  t_last$t <- prob$tspan[1L]
  rec <- DiscreteCallback(function(u, t, integ) {
    t_last$t <- integ$t
    FALSE
  }, function(integ) NULL, save_positions = c(FALSE, FALSE))
  opts$callback <- CallbackSet(rec, opts$callback)
  post <- function(sol) {
    env <- .sol_env(sol)
    # the last attempt was accepted, truncated or rejected where the driver stopped
    .noise_settle(ns, t_last$t)
    env$W <- list(t = ns$log_t, W = do.call(rbind, ns$log_W))
    env$alg <- alg
    env$prob <- prob
    sol
  }
  list(prob = dprob, alg = runner, opts = opts, post = post)
}

# Run `expr` with the random number generator seeded, restoring the state.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", old, envir = globalenv())
  })
  set.seed(seed)
  expr
}

.as_rode_solution <- function(sol) {
  class(sol) <- c("RODESolution", class(sol))
  sol
}

#' @rdname sde_algorithms
#' @param a An `SDEProblem`.
#' @param b An SDE algorithm; [SOSRI()] or [LambaEM()] by default.
#' @param ... Solver options, see Details.
#' @export
solve.SDEProblem <- function(a, b, ...) {
  opts <- utils::modifyList(a$kwargs %||% list(), list(...))
  alg <- if (missing(b) || is.null(b)) .sde_default_alg(a, opts$alg_hints) else b
  opts$alg_hints <- NULL
  seed <- opts$seed
  opts$seed <- NULL
  .with_seed(seed, .as_rode_solution(.sde_integrate(a, alg, opts)))
}

#' @export
print.RODESolution <- function(x, ...) {
  env <- .sol_env(x)
  cat("RODESolution with retcode ", env$retcode, " (", env$alg$name %||% "", ")\n", sep = "")
  n <- length(env$t)
  cat("Saved time points: ", n, if (n) paste0(" in [", format(env$t[1L]), ", ", format(env$t[n]), "]"), "\n", sep = "")
  if (n) {
    cat("Final state:\n")
    print(env$U[n, ])
  }
  invisible(x)
}
