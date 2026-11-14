# Newmark-beta and generalized alpha methods of OrdinaryDiffEqNewmark
# (algorithms.jl, alg_utils.jl, newmark_nlsolve.jl, newmark_perform_step.jl)
# for second order problems u'' = f1(u', u, p, t) given as
# SecondOrderODEProblem or DynamicalODEProblem (state c(v, u)).
#
# A step solves for the new acceleration a the residual
#   (1 - am) a + am an - f1((1 - af) v1 + af v0, (1 - af) u1 + af u0, t + (1 - af) dt)
# with u1 = u0 + dt v0 + dt^2 / 2 ((1 - 2 beta) an + 2 beta a) and
# v1 = v0 + dt ((1 - gamma) an + gamma a), where an is f1 at the start of
# the step (re-derived each step, as upstream). The nonlinear solve uses the
# `nlsolve` algorithm (NewtonRaphson by default); its f1 calls are not
# counted, as in Julia. The local error is the Zienkiewicz-Xie estimate with
# Julia's coefficients. Dense output is cubic Hermite interpolation from f at
# both ends of a step (OrdinaryDiffEqCore's _ode_addsteps! fallback). Only
# the identity mass matrix: Julia's DynamicalODEFunction rejects others.

#' Newmark-beta and generalized alpha methods
#'
#' Implicit methods of OrdinaryDiffEqNewmark for second order problems
#' `u'' = f1(u', u, p, t)` given as [SecondOrderODEProblem()] or
#' [DynamicalODEProblem()] (state `c(v, u)`). `NewmarkBeta()` is the classical
#' Newmark-beta method, second order for `gamma = 1/2` and first order
#' otherwise. `GeneralizedAlpha()` is the generalized alpha method of Chung
#' and Hulbert, set by the spectral radius at infinity `rho_inf`, by the HHT
#' parameter `alpha_hht` or by its four parameters; as upstream, the step
#' starts from the acceleration of the accepted state, so it is second order
#' only when `gamma (1 - alpha_f) = (1 - alpha_m) / 2`. Both methods are
#' adaptive with the Zienkiewicz-Xie error estimate; parameters for which
#' that estimate vanishes stop unless `adaptive = FALSE`. With
#' `adaptive = FALSE` the starting acceleration is evaluated at the end time
#' `t + dt` of the step, as upstream does with fixed steps, so on problems
#' with explicit time dependence the fixed step methods are first order.
#'
#' @param beta,gamma The Newmark parameters, `0 <= beta <= 1/2` and
#'   `0 <= gamma <= 1`.
#' @param alpha_m,alpha_f The generalized alpha weights of the acceleration
#'   and of the state.
#' @param rho_inf Spectral radius at infinity in `[0, 1]`; sets the four
#'   parameters (Chung and Hulbert).
#' @param alpha_hht HHT alpha in `[-1/3, 0]`; sets `alpha_m = 0`.
#' @param nlsolve Nonlinear solver for the implicit stage, a
#'   [NewtonRaphson()] by default.
#' @return An algorithm object.
#' @examples
#' prob <- SecondOrderODEProblem(function(v, u, p, t) -u, 0, 1, c(0, 1))
#' sol <- solve(prob, NewmarkBeta(), dt = 0.05, adaptive = FALSE)
#' sol$u[[length(sol$t)]] - c(-sin(1), cos(1))
#' solve(prob, GeneralizedAlpha(rho_inf = 0.8), dt = 0.1, adaptive = FALSE)$stats$nf
#' @name newmark_methods
NULL

#' @rdname newmark_methods
#' @export
NewmarkBeta <- function(beta = 0.25, gamma = 0.5, nlsolve = NewtonRaphson()) {
  if (!(beta >= 0 && beta <= 0.5)) stop("Beta outside admissible range [0, 0.5]", call. = FALSE)
  if (!(gamma >= 0 && gamma <= 1)) stop("Gamma outside admissible range [0, 1.0]", call. = FALSE)
  .newmark_alg("NewmarkBeta", 0, 0, beta, gamma, nlsolve,
               order = if (isTRUE(all.equal(gamma, 1 / 2, tolerance = sqrt(.Machine$double.eps)))) 2L else 1L)
}

#' @rdname newmark_methods
#' @export
GeneralizedAlpha <- function(alpha_m = NULL, alpha_f = NULL, beta = NULL, gamma = NULL,
                             rho_inf = NULL, alpha_hht = NULL, nlsolve = NewtonRaphson()) {
  if (!is.null(rho_inf)) {
    if (!(rho_inf >= 0 && rho_inf <= 1)) stop("rho_inf must be in [0, 1]; got ", rho_inf, call. = FALSE)
    r <- rho_inf
    alpha_m <- (2 * r - 1) / (r + 1)
    alpha_f <- r / (r + 1)
    gamma <- 1 / 2 - alpha_m + alpha_f
    beta <- (1 / 2 + alpha_f - alpha_m)^2 / 4
  } else if (!is.null(alpha_hht)) {
    if (!(alpha_hht >= -1 / 3 && alpha_hht <= 0))
      stop("HHT alpha must be in [-1/3, 0]; got ", alpha_hht, call. = FALSE)
    alpha_m <- 0
    alpha_f <- -alpha_hht
    gamma <- (1 - 2 * alpha_hht) / 2
    beta <- (1 - alpha_hht)^2 / 4
  } else if (is.null(alpha_m) || is.null(alpha_f) || is.null(beta) || is.null(gamma)) {
    stop("Must provide either rho_inf, alpha_hht, or all four of alpha_m, alpha_f, beta, gamma", call. = FALSE)
  }
  if (!(alpha_m <= alpha_f))
    stop("Unconditional stability requires alpha_m <= alpha_f; got alpha_m=", alpha_m, ", alpha_f=", alpha_f, call. = FALSE)
  if (!(alpha_f <= 0.5)) stop("Unconditional stability requires alpha_f <= 1/2; got alpha_f=", alpha_f, call. = FALSE)
  beta_min <- (1 / 2 + alpha_f - alpha_m)^2 / 4
  if (!(beta >= beta_min)) stop("Unconditional stability requires beta >= ", beta_min, "; got beta=", beta, call. = FALSE)
  if (!(gamma >= 0 && gamma <= 1)) stop("gamma outside admissible range [0, 1]; got gamma=", gamma, call. = FALSE)
  second <- isTRUE(all.equal(gamma * (1 - alpha_f), (1 - alpha_m) / 2, tolerance = sqrt(.Machine$double.eps)))
  .newmark_alg("GeneralizedAlpha", alpha_m, alpha_f, beta, gamma, nlsolve, order = if (second) 2L else 1L)
}

# zx_band_check: the error estimate is blind when both leading coefficients
# nearly vanish.
.newmark_band_check <- function(beta, gamma, am, af) {
  kappa <- (1 - af) / (1 - am)
  Cp <- beta - 1 / (6 * kappa)
  Cv <- gamma - 1 / (2 * kappa)
  if (abs(Cp) < 1 / 24 && abs(Cv) < 1 / 100)
    stop("The Zienkiewicz-Xie error estimate is blind for these parameters: its position coefficient ",
         "beta - 1/(6 kappa) = ", Cp, " and velocity coefficient gamma - 1/(2 kappa) = ", Cv,
         ", kappa = (1 - alpha_f)/(1 - alpha_m), are both too close to their vanishing loci, so the ",
         "estimated error can be arbitrarily small next to the actual error and the solver would return ",
         "inaccurate results with a Success code. Move beta at least ", 1 / 24, " away from ",
         1 / (6 * kappa), " or gamma at least ", 1 / 100, " away from ", 1 / (2 * kappa),
         ", or pass adaptive = FALSE.", call. = FALSE)
  invisible(NULL)
}

# Hermite interpolation from f at both ends of the step (k holds f, p and the
# start time; the derivatives are evaluated when needed, uncounted as in
# Julia).
.newmark_interp <- function(theta, dt, y0, y1, k)
  .hermite_interp(theta, dt, y0, y1, list(k$f(y0, k$p, k$t), k$f(y1, k$p, k$t + dt)))

.newmark_alg <- function(name, am, af, beta, gamma, nlsolve, order) {
  kappa <- (1 - af) / (1 - am)
  Cp0 <- beta - 1 / (6 * kappa)
  Cv0 <- max(abs(gamma - 1 / (2 * kappa)), 1 / 100)
  .new_alg(name, order, adaptive = TRUE, interp = .newmark_interp, partitioned = TRUE,
           qsteady_max = 6 / 5, newmark = c(alpha_m = am, alpha_f = af, beta = beta, gamma = gamma),
           # is_mass_matrix_alg (the identity is the only mass matrix accepted)
           mass_matrix_ok = TRUE,
           on_event = .part_on_event,
    init = function(ctx, t, z) {
      .dyn_parts(ctx, z)
      if (ctx$adaptive) .newmark_band_check(beta, gamma, am, af)
      .count_nf2(ctx, 0L)
      list(fsal = ctx$f(z, ctx$p, t), nf0 = 1L)
    },
    step = function(ctx, t, z, dt, cache) {
      i <- .dyn_idx(ctx, z)
      f1 <- ctx$f1; p <- ctx$p
      vn <- z[i$iv]; un <- z[i$iu]
      # the acceleration at the start; with fixed steps upstream evaluates
      # it at t + dt
      an <- f1(vn, un, p, if (ctx$adaptive) t else t + dt)
      dtsq2 <- dt^2 / 2
      tf <- t + (1 - af) * dt
      resid <- function(a) {
        u1 <- un + dt * vn + dtsq2 * ((1 - 2 * beta) * an + 2 * beta * a)
        v1 <- vn + dt * ((1 - gamma) * an + gamma * a)
        ((1 - am) * a + am * an) - f1((1 - af) * v1 + af * vn, (1 - af) * u1 + af * un, p, tf)
      }
      # the nonlinear solve steps before its termination test (no check at
      # the start), so a start acceleration that already solves the residual
      # still moves by the Newton correction, as NonlinearSolve
      nl <- .nl_run(nlsolve, resid, NULL, an, list())
      if (!identical(nl$retcode, "Success")) return(list(u = z, fail = TRUE, cache = cache, nf = 1L))
      a <- nl$u
      u1 <- un + dt * vn + dtsq2 * ((1 - 2 * beta) * an + 2 * beta * a)
      v1 <- vn + dt * ((1 - gamma) * an + gamma * a)
      znew <- c(v1, u1)
      EEst <- NULL
      if (ctx$adaptive) {
        da <- a - an
        est <- .err_norm(dt^2 * Cp0 * da, un, u1, .part_tol(ctx$abstol, i$iu), .part_tol(ctx$reltol, i$iu))
        est <- max(est, .err_norm(dt * Cv0 * da, vn, v1, .part_tol(ctx$abstol, i$iv), .part_tol(ctx$reltol, i$iv)))
        EEst <- if (est == 0) 1 else est
      }
      list(u = znew, EEst = EEst, k = list(f = ctx$f, p = p, t = t), cache = cache, nf = 1L)
    })
}
