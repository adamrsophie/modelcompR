# Adjoint based global error estimation and control, after GlobalDiffEq.jl
# v1.7.3: src/adjoint.jl (GlobalAdjoint, adjoint_error_estimate, the
# random projections, scope and control) and
# ext/GlobalDiffEqSciMLSensitivityExt.jl (_adjoint_defect_projection), on
# the adjoint of R/sensitivity-adjoint.R.
#
# For a direction v the error component <v, y(T) - P(T)> is the dual
# weighted residual  int lambda(t)' (f(P(t)) - P'(t)) dt  with the adjoint
# lambda(T) = v. It is the sensitivity of G = <v, u(T)> to a scalar forcing
# q in u' = f(u, p, t) + q r(t), r = f(P) - P', at q = 0: the parameter
# gradient of adjoint_sensitivities on that forced problem.

# _orthogonal_directions: min(samples, n) orthonormal random directions,
# with two Gram-Schmidt passes. rng(n) returns n standard normal draws.
.gda_directions <- function(u0, requested, rng) {
  n <- length(u0)
  dirs <- list()
  for (s in seq_len(min(requested, n))) {
    accepted <- FALSE
    for (attempt in 1:10) {
      d <- as.numeric(rng(n))
      for (pass in 1:2) for (prev in dirs) d <- d - sum(prev * d) * prev
      dn <- sqrt(sum(d * d))
      if (dn > sqrt(.Machine$double.eps)) {
        dirs[[length(dirs) + 1L]] <- d / dn
        accepted <- TRUE
        break
      }
    }
    if (!accepted) stop("failed to generate orthogonal random directions", call. = FALSE)
  }
  dirs
}

# _sphere_projection_expectation: E|<v, x>| for a fixed unit v and x
# uniform on the unit sphere of the dimension, up to the common factor.
.gda_sphere_expectation <- function(dimension) {
  if (!(dimension > 0)) stop("dimension must be positive", call. = FALSE)
  if (dimension == 1) return(1)
  odd <- dimension %% 2 == 1
  e <- if (odd) 1 else 2 / pi
  num <- if (odd) 1 else 2
  while (num <= dimension - 2) {
    e <- e * (num / (num + 1))
    num <- num + 2
  }
  e
}

.gda_validate_problem <- function(prob) .gde_validate_state(prob, "GlobalAdjoint")

.gda_validate_solution <- function(sol) {
  if (!successful_retcode(sol)) stop("the forward solve failed with retcode ", sol$retcode, call. = FALSE)
  .gda_validate_problem(sol$prob)
  if (length(sol$t) < 2L) stop("the forward solution must save at least its two endpoints", call. = FALSE)
  invisible()
}

# _adjoint_defect_projection: the forced problem shares the trajectory of
# `sol` (a new solution object on the same steps) with p = 0; its exact
# Jacobians are the Jacobian of the problem and r(t). The forward events are
# not passed on, as Julia's adjoint_sensitivities gets no callback here.
.gda_projection <- function(sol, sensealg, adjoint_alg, direction, abstol, reltol, terminal_time) {
  env <- .sol_env(sol)
  prob <- env$prob
  realp <- prob$p
  fns <- .sens_fns(prob, sensealg)
  f0 <- fns$f
  jac0 <- fns$jac
  dense <- .gde_sol_dense(sol, right = FALSE)
  resid <- function(t) f0(dense$value(t), realp, t) - dense$deriv(t)
  fprob <- prob
  fprob$f <- ODEFunction(function(u, q, t) f0(u, realp, t) + q[1L] * resid(t),
                         jac_fn = function(u, q, t) jac0(u, realp, t),
                         paramjac = function(u, q, t) matrix(resid(t), ncol = 1L))
  fprob$p <- 0
  fprob$sys <- NULL
  fprob$kwargs$callback <- NULL
  fsol <- .new_ode_solution(t = env$t, u = env$U, prob = fprob, alg = env$alg, retcode = env$retcode,
                            steps = env$steps, stats = env$stats, p = 0, dense = env$dense)
  res <- adjoint_sensitivities(fsol, adjoint_alg, t = terminal_time,
                               dgdu_discrete = function(u, p, t, i) direction,
                               abstol = abstol, reltol = reltol, sensealg = sensealg)
  as.numeric(res$dp)[1L]
}

# _adjoint_error_estimate
.gda_estimate <- function(sol, alg, directions, abstol, reltol, terminal_time = sol$prob$tspan[2L]) {
  proj <- vapply(directions, function(d)
    .gda_projection(sol, alg$sensealg, alg$adjoint_alg, d, abstol, reltol, terminal_time), 0)
  factor <- .gda_sphere_expectation(length(directions)) / .gda_sphere_expectation(length(sol$prob$u0))
  factor * sqrt(sum(proj * proj))
}

# _adjoint_trajectory_errors: the estimate at each saved time, from the
# discrete cost placed there (no error at the initial time)
.gda_trajectory_errors <- function(sol, alg, directions, abstol, reltol) {
  ts <- sol$t
  out <- numeric(length(ts))
  for (j in seq_along(ts)[-1L]) out[j] <- .gda_estimate(sol, alg, directions, abstol, reltol, terminal_time = ts[j])
  out
}

# _solve_on_grid: a fixed step solve stepping on exactly the steps of `grid`
.gda_solve_on_grid <- function(prob, inner_alg, grid, abstol, reltol) {
  dts <- diff(grid)
  integ <- init(prob, inner_alg, adaptive = FALSE, dt = dts[1L], dense = TRUE, save_everystep = TRUE,
                save_start = TRUE, save_end = TRUE, abstol = abstol, reltol = reltol)
  for (dt in dts) {
    set_proposed_dt(integ, dt)
    step_integrator(integ)
  }
  # Julia saves the end point with the last step; the R integrator saves it
  # when it is finalized
  sol <- if (integ$t >= grid[length(grid)]) solve_integrator(integ) else integ$sol
  if (!successful_retcode(sol) && integ$t >= grid[length(grid)]) {
    env <- .sol_env(sol)
    env$retcode <- "Success"
  }
  sol
}

#' Global error estimation and control with adjoints
#'
#' `GlobalAdjoint(alg)` wraps an ODE algorithm with adjoint based global
#' error estimation and control (Cao and Petzold 2004), as GlobalDiffEq.jl.
#' The problem is solved with local tolerances, the 2-norm of the endpoint
#' error is estimated from `samples` orthogonal random projections, and the
#' local tolerances (which default to `gtol`) are tightened until the
#' estimate is at most `gtol`. With more than one state the estimate is
#' probabilistic; two samples are accurate within a factor of ten with
#' probability above 99%. With as many samples as states it is the norm of
#' the full error estimate.
#'
#' For a projection direction `v` the error component `<v, y(T) - P(T)>` is
#' the dual weighted residual, the integral of `lambda(t)' (f(P) - P')`
#' with the adjoint `lambda(T) = v`, `P` the dense numerical solution and
#' `f(P) - P'` its defect. It is computed by [adjoint_sensitivities()] as
#' the sensitivity to a scalar forcing that adds the defect to the
#' dynamics; `sensealg` (by default [InterpolatingAdjoint()]) and
#' `adjoint_alg` select the adjoint method and its solver.
#'
#' `solve(prob, GlobalAdjoint(alg, gtol = ...))` returns a solution with
#' `sol$global_error` (and [has_global_error()] is `TRUE`): with
#' [EndpointError] the estimate in its last element (the others are zero),
#' with [TrajectoryError] the estimate at every saved time. With
#' [ToleranceRefinement] the solution is an adaptive solve at the accepted
#' tolerances; with [StepGridRefinement] a fixed step solve on the accepted
#' step grid. `adjoint_error_estimate(prob, alg)` solves `prob` and returns
#' the estimate. The solver must have a dense output with a first
#' derivative (see [GlobalErrorEstimation()]); forward time ODEs with the
#' standard mass matrix and real vector states are supported.
#'
#' @param alg The ODE algorithm.
#' @param adjoint_alg The algorithm of the adjoint solves.
#' @param sensealg An adjoint method for [adjoint_sensitivities()], or
#'   `NULL` for its default.
#' @param samples The number of random projections (at most the number of
#'   states are used).
#' @param rng A function of `n` returning `n` standard normal draws, which
#'   give the random directions (Julia takes a random number generator).
#' @param scope [EndpointError] or [TrajectoryError].
#' @param control [ToleranceRefinement] or [StepGridRefinement].
#' @param gtol The requested absolute endpoint error, needed by `solve()`.
#' @param maxiters The maximum number of tolerance refinements.
#' @param safety The safety factor of each tolerance reduction.
#' @param adjoint_abstol,adjoint_reltol Tolerances of the adjoint solves.
#' @return An algorithm object.
#' @references Y. Cao and L. Petzold, A posteriori error estimation and
#'   global error control for ordinary differential equations by the adjoint
#'   method, SIAM Journal on Scientific Computing 26, 2004.
#' @examples
#' prob <- ODEProblem(function(u, p, t) p * u, 1, c(0, 1), 1)
#' alg <- GlobalAdjoint(Tsit5(), samples = 1)
#' adjoint_error_estimate(prob, alg, abstol = 1e-4, reltol = 1e-4)
#' sol <- solve(prob, Tsit5(), abstol = 1e-4, reltol = 1e-4)
#' abs(sol$u[[length(sol$t)]] - exp(1))
#' @export
GlobalAdjoint <- function(alg, adjoint_alg = alg, sensealg = NULL, samples = 2, rng = stats::rnorm,
                          scope = EndpointError, control = ToleranceRefinement, gtol = NULL, maxiters = 6,
                          safety = 0.8, adjoint_abstol = if (is.null(gtol)) 1e-10 else min(gtol / 100, 1e-10),
                          adjoint_reltol = if (is.null(gtol)) 1e-8 else min(gtol / 100, 1e-8)) {
  if (!inherits(alg, "ODEAlgorithm")) stop("GlobalAdjoint wraps an ODE algorithm", call. = FALSE)
  if (!(is.numeric(samples) && length(samples) == 1L && is.finite(samples) && samples == round(samples)))
    stop("samples must be an integer", call. = FALSE)
  if (!(samples > 0)) stop("samples must be positive", call. = FALSE)
  if (!inherits(scope, "GlobalAdjointScope")) stop("scope must be EndpointError or TrajectoryError", call. = FALSE)
  if (!inherits(control, "GlobalAdjointControl"))
    stop("control must be ToleranceRefinement or StepGridRefinement", call. = FALSE)
  if (!(is.null(gtol) || .gde_positive_finite(gtol))) stop("gtol must be a positive finite real number", call. = FALSE)
  if (!(is.numeric(maxiters) && length(maxiters) == 1L && maxiters == round(maxiters) && maxiters > 0))
    stop("maxiters must be a positive integer", call. = FALSE)
  if (!(is.numeric(safety) && length(safety) == 1L && is.finite(safety) && safety > 0 && safety < 1))
    stop("safety must be between zero and one", call. = FALSE)
  .gde_validate_tolerances(adjoint_abstol, adjoint_reltol, "adjoint")
  if (!is.function(rng)) stop("rng must be a function of n returning n standard normal draws", call. = FALSE)
  options <- list(maxiters = as.integer(maxiters), safety = safety, adjoint_abstol = adjoint_abstol,
                  adjoint_reltol = adjoint_reltol)
  structure(list(name = "GlobalAdjoint", alg = alg, adjoint_alg = adjoint_alg, sensealg = sensealg,
                 samples = as.integer(samples), rng = rng, scope = scope, control = control, gtol = gtol,
                 options = options, order = alg$order, adaptive = TRUE, has_global_error = TRUE,
                 solve_override = .gda_solve),
            class = c("GlobalAdjoint", "ODEAlgorithm"))
}

#' @rdname GlobalAdjoint
#' @param prob An `ODEProblem`.
#' @param ... Options of the forward solve.
#' @param abstol,reltol Tolerances of the forward solve.
#' @export
adjoint_error_estimate <- function(prob, alg, ..., abstol = 1e-6, reltol = 1e-3,
                                   adjoint_abstol = alg$options$adjoint_abstol,
                                   adjoint_reltol = alg$options$adjoint_reltol) {
  if (!inherits(alg, "GlobalAdjoint")) stop("adjoint_error_estimate(prob, alg) needs a GlobalAdjoint algorithm",
                                           call. = FALSE)
  .gda_validate_problem(prob)
  sopts <- utils::modifyList(list(...), .GDE_DENSE_OPTS)
  sol <- do.call(solve, c(list(prob, alg$alg), list(abstol = abstol, reltol = reltol), sopts))
  .gda_validate_solution(sol)
  directions <- .gda_directions(prob$u0, alg$samples, alg$rng)
  .gda_estimate(sol, alg, directions, adjoint_abstol, adjoint_reltol)
}

# SciMLBase.__solve for GlobalAdjoint
.gda_solve <- function(prob, alg, opts) {
  if (is.null(alg$gtol)) stop("GlobalAdjoint requires a positive `gtol` constructor keyword", call. = FALSE)
  abstol <- opts$abstol %||% alg$gtol
  reltol <- opts$reltol %||% alg$gtol
  opts$abstol <- NULL
  opts$reltol <- NULL
  .gde_validate_tolerances(abstol, reltol, "local")
  .gda_validate_problem(prob)
  o <- alg$options
  directions <- .gda_directions(prob$u0, alg$samples, alg$rng)
  trial_opts <- utils::modifyList(opts, .GDE_DENSE_OPTS)
  solve_with <- function(o2, la, lr) do.call(solve, c(list(prob, alg$alg), o2, list(abstol = la, reltol = lr)))
  trial <- NULL
  estimate <- function(la, lr) {
    trial <<- solve_with(trial_opts, la, lr)
    .gda_validate_solution(trial)
    .gda_estimate(trial, alg, directions, o$adjoint_abstol, o$adjoint_reltol)
  }
  finish <- function(la, lr, est) {
    final <- if (identical(alg$control, StepGridRefinement)) .gda_solve_on_grid(prob, alg$alg, trial$t, la, lr) else
      if (identical(alg$scope, TrajectoryError)) solve_with(trial_opts, la, lr) else solve_with(opts, la, lr)
    ge <- if (identical(alg$scope, TrajectoryError)) {
      .gda_trajectory_errors(final, alg, directions, o$adjoint_abstol, o$adjoint_reltol)
    } else {
      g <- numeric(length(final$t))
      g[length(g)] <- est
      g
    }
    env <- .sol_env(final)
    env$global_error <- ge
    final
  }
  .gde_refine(estimate, finish, alg$gtol, o, abstol, reltol, what = "adjoint global error estimate")
}
