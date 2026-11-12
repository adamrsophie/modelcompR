# The generic explicit Runge-Kutta stepper for the tableau methods of
# R/solvers-erk-methods.R (the explicit methods of OrdinaryDiffEqLowOrderRK,
# OrdinaryDiffEqHighOrderRK, OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK,
# OrdinaryDiffEqFeagin, OrdinaryDiffEqVerner and OrdinaryDiffEqSIMDRK).
#
# A method is its Butcher tableau (A, b, c), the embedded error weights
# btilde (err = dt sum_i btilde_i k_i) and the controller settings of the
# Julia method. A last stage evaluated at the new state (c = 1, row = b) is
# the first stage of the next step (FSAL); methods without it take a fresh
# first stage and interpolate linearly, as their Julia caches keep no end
# derivative. The SSP methods are run from their full tableaus; the low
# storage methods run the register forms of OrdinaryDiffEqLowStorageRK
# (R/solvers-lowstorage.R), whose rounding differs from the tableau form.

#' Explicit Runge-Kutta methods
#'
#' The explicit Runge-Kutta methods of the OrdinaryDiffEq sublibraries, for
#' [solve()] on ODE problems. All take `stage_limiter(u, integ, p, t)` and
#' `step_limiter(u, integ, p, t)` (functions returning the limited state,
#' applied to every stage and to every step, as the SSP methods use them)
#' and `thread` (ignored; R evaluates the stages in turn). The low storage
#' methods apply only the step limiter, as their out of place forms in
#' OrdinaryDiffEq.
#'
#' * Low order (OrdinaryDiffEqLowOrderRK): `Heun()`, `Ralston()`,
#'   `Ralston4()`, `RKM()`, `OwrenZen3()`, `OwrenZen4()`, `OwrenZen5()`,
#'   `BS5()`, `Alshina2()`, `Alshina3()`, `Alshina6()`, `Stepanov5()`,
#'   `SIR54()`, `MSRK5()`, `MSRK6()`, `PSRK3p5q4()`, `PSRK3p6q5()`,
#'   `PSRK4p7q6()`.
#' * High order: `DP8()`, `TanYam7()`, `TsitPap8()`, `PFRK87()`,
#'   `Feagin10()`, `Feagin12()`, `Feagin14()`, `RKV76IIa()`.
#' * Strong stability preserving (OrdinaryDiffEqSSPRK): `SSPRK22()`,
#'   `SSPRK33()`, `SSPRK53()`, `SSPRK53_2N1()`, `SSPRK53_2N2()`,
#'   `SSPRK53_H()`, `SSPRK63()`, `SSPRK73()`, `SSPRK83()`, `SSPRK43()`,
#'   `SSPRK432()`, `SSPRK932()`, `SSPRK54()`, `SSPRK104()`, `KYKSSPRK42()`.
#' * Low storage (OrdinaryDiffEqLowStorageRK): `ORK256()`,
#'   `CarpenterKennedy2N54()`, `SHLDDRK64()`, `SHLDDRK52()`,
#'   `SHLDDRK_2N()`, `DGLDDRK73_C()`, `DGLDDRK84_C()`, `DGLDDRK84_F()`,
#'   `NDBLSRK124()`, `NDBLSRK134()`, `NDBLSRK144()`, `CFRLDDRK64()`,
#'   `TSLDDRK74()`, `RK46NL()`, the `CKLLSRK*` and `ParsaniKetchesonDeconinck*`
#'   families, the `AlJahdali*` families and the adaptive `RDPK3Sp35()`,
#'   `RDPK3SpFSAL35()`, `RDPK3Sp49()`, `RDPK3SpFSAL49()`, `RDPK3Sp510()`,
#'   `RDPK3SpFSAL510()` (with their PID step size controllers).
#' * `AutoDP5(alg)` switches between `DP5()` and the stiff method `alg`;
#'   `SplitEuler()` is the explicit Euler method for split problems.
#'
#' Each method has the order, embedded error estimate and step size
#' controller of its Julia version; dense output is cubic Hermite
#' interpolation for methods that evaluate the derivative at the new state
#' and linear otherwise.
#'
#' @param stage_limiter,step_limiter Limiters, see Details.
#' @param thread Threading option (ignored).
#' @param ... Controller settings stored in the algorithm (`qmax`, `qmin`,
#'   `gamma`, `beta1`, `beta2`).
#' @return An algorithm object.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
#' solve(prob, OwrenZen5())(1)
#' solve(prob, SSPRK33(), dt = 0.1)(1)
#' solve(prob, RDPK3Sp35())(1)
#' @name explicit_rk_methods
NULL

.erk_alg <- function(name, stage_limiter = NULL, step_limiter = NULL, ..., tableau = NULL) {
  tb <- tableau %||% .ERK_TABLEAUS[[name]]
  s <- tb$s
  A <- matrix(tb$A, s, s, byrow = TRUE)
  b <- tb$b; cc <- tb$c; bt <- tb$btilde
  fsal <- s > 1L && abs(cc[s] - 1) < 1e-14 && max(abs(A[s, ] - b)) < 1e-13 && abs(b[s]) < 1e-14
  adaptive <- isTRUE(tb$adaptive) && length(bt)
  nonlin <- isTRUE(tb$nonlinear_est)
  est <- .ERK_ESTIMATORS[[name]]
  step <- function(ctx, t, u, dt, cache) {
    f <- ctx$f; p <- ctx$p
    n <- length(u)
    K <- matrix(0, n, s)
    K[, 1L] <- if (!is.null(cache$fsal)) cache$fsal else f(u, p, t)
    nf <- if (is.null(cache$fsal)) 1L else 0L
    last <- if (fsal) s - 1L else s
    for (i in seq_len(last)[-1L]) {
      ui <- u + dt * .lincomb(K, A[i, ], which(A[i, seq_len(i - 1L)] != 0))
      ti <- .fma(cc[i], dt, t)
      if (!is.null(stage_limiter)) ui <- stage_limiter(ui, NULL, p, ti)
      K[, i] <- f(ui, p, ti)
      nf <- nf + 1L
    }
    unew <- u + dt * .lincomb(K, b)
    if (!is.null(stage_limiter)) unew <- stage_limiter(unew, NULL, p, t + dt)
    if (fsal) {
      K[, s] <- f(unew, p, t + dt)
      nf <- nf + 1L
    }
    res <- list(u = unew, k = if (fsal) list(K[, 1L], K[, s]) else NULL,
                cache = list(fsal = if (fsal) K[, s]), nf = nf)
    if (adaptive) {
      if (!is.null(est)) {
        res$EEst <- est(K, dt, u, unew, ctx, t)
        if (identical(name, "SSPRK932")) res$nf <- res$nf + 1L
      }
      else res$err <- dt * .lincomb(K, bt)
    }
    res
  }
  ctrl <- list()
  if (identical(tb$controller, "PIDController")) {
    ctrl <- list(controller = "PID", pid_beta = tb$pid_beta, accept_safety = tb$pid_accept_safety %||% 0.81)
  }
  parts <- list(step = step, interp = if (fsal) .hermite_interp else .sde_linear_interp,
                init = function(ctx, t, u) list())
  # the low storage methods run their register forms
  if (is.null(tableau) && !is.null(.LSRK[[name]])) parts <- .lsrk_parts(name)
  alg <- do.call(.new_alg, c(list(name, tb$order, adaptive = adaptive), parts,
                             ctrl,
                             list(adaptive_order = tb$adaptive_order, qmax = tb$qmax, qmin = tb$qmin,
                                  gamma = tb$gamma, beta1 = tb$beta1, beta2 = tb$beta2,
                                  ssp_coefficient = tb$ssp, step_limiter = step_limiter,
                                  stage_limiter = stage_limiter, tableau = tb)))
  utils::modifyList(alg, list(...))
}

# Methods whose Julia error estimate is not a single embedded difference
# (perform_step! of BS5ConstantCache, DP8ConstantCache, SSPRK932ConstantCache).
.erk_nrm <- function(e, u, unew, ctx) .err_norm(e, u, unew, ctx$abstol, ctx$reltol)

.ERK_ESTIMATORS <- list(
  # the larger of the two embedded estimates
  BS5 = function(K, dt, u, unew, ctx, t) {
    bhat <- c(-0.00234375, 0, 0.0103760754048583, -0.016490384615384615, 0.018984375, -0.010526315789473684, 0, 0)
    btil <- c(-0.0019478942125547064, 0, 0x1.2881fd89b5f5dp-7, -0.015478743126634073, 0x1.5bb688355021dp-6,
              -0.013276905111429694, -0.005480370770138146, 0.005912495780636172)
    max(.erk_nrm(dt * .lincomb(K, bhat), u, unew, ctx), .erk_nrm(dt * .lincomb(K, btil), u, unew, ctx))
  },
  # the order 5 and order 3 estimates combined as err5^2 / sqrt(err5^2 + 0.01 err3^2)
  DP8 = function(K, dt, u, unew, ctx, t) {
    z <- numeric(13)
    er <- z; er[c(1, 6:12)] <- c(0.01312004499419488, -1.2251564463762044, -0.4957589496572502, 1.6643771824549864,
                                 -0.35032884874997366, 0.3341791187130175, 0.08192320648511571, -0x1.6e44f50ab6bc2p-6)
    bt <- z; bt[c(1, 6:12)] <- c(-0x1.84b641fbfa1f1p-3, 4.450312892752409, 1.8915178993145003, -5.801203960010585,
                                 -0x1.b0d3a26abb717p-2, -0.1521609496625161, 0x1.9c657697fe72dp-3, 0.022651792198360825)
    e5 <- .erk_nrm(dt * .lincomb(K, er), u, unew, ctx)
    e3 <- .erk_nrm(dt * .lincomb(K, bt), u, unew, ctx)
    if (abs(e5) < 0x1.56e1fc2f8f359p-997 && abs(e3) < 0x1.56e1fc2f8f359p-997) 0 else e5^2 / sqrt(e5^2 + 0.01 * e3^2)
  },
  # the extra stage at the end of the first sweep gives an embedded solution
  SSPRK932 = function(K, dt, u, unew, ctx, t) {
    u6 <- u + dt / 6 * .lincomb(K, rep(1, 6))
    k7 <- ctx$f(u6, ctx$p, t + dt)
    ut <- (u + 6 * u6 + 6 * dt * k7) / 7
    attr(ut, "nf") <- 1L
    .erk_nrm(ut - unew, u, unew, ctx)
  })

#' @rdname explicit_rk_methods
#' @param alg The stiff method of the automatic switch.
#' @export
AutoDP5 <- function(alg, ...) AutoSwitch(DP5(), alg, ...)

#' @rdname explicit_rk_methods
#' @export
SplitEuler <- function() {
  # explicit Euler on f1 + f2 of a split problem (the sum is the problem's f)
  a <- Euler()
  a$name <- "SplitEuler"
  class(a) <- c("SplitEuler", class(a)[-1L])
  a
}
