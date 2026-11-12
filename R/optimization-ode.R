# OptimizationODE (lib/OptimizationODE/src/OptimizationODE.jl, v0.1.12): the
# minimizer as the steady state of the gradient flow du/dt = -grad f(u),
# integrated with an ODE solver as the DynamicSS of SteadyStateDiffEq does
# (src/solve.jl): tolerances abstol = 1e-8, reltol = 1e-6 for the
# integration, and the integration stops after the first step whose
# derivative du satisfies the NormTerminationMode of NonlinearSolveBase with
# the maximum norm, max|du| <= 1e-8 or max|du| <= 1e-6 * max|du + u|. The
# constrained variants integrate a mass matrix DAE (DAEMassMatrix) or a fully
# implicit DAE in which the last unknowns follow the equality constraints.

.odeopt_alg <- function(name, alg, dae = FALSE) {
  a <- .opt_alg(name, alg = alg)
  class(a) <- c(name, "OptimizationAlgorithm")
  a$solver <- .odeopt_solve
  a$traits <- list(grad = TRUE, cons = dae, consjac = dae)
  a
}

#' Optimization by integrating the gradient flow
#'
#' The solvers of OptimizationODE find a minimizer as the steady state of
#' the gradient flow `du/dt = -gradient(u)`, integrated from `u0` with an ODE
#' solver until the derivative is small: the largest component of `du` at
#' most `1e-8`, or at most `1e-6` times the largest component of `du + u`,
#' as `DynamicSS` of SteadyStateDiffEq (the integration itself uses
#' `abstol = 1e-8`, `reltol = 1e-6`). The solve options `dt` (the step of
#' [Euler()], the initial step of the adaptive methods) and `maxiters` (the
#' largest number of steps) are passed to the ODE solver.
#'
#' * `ODEOptimizer(alg)` uses the ODE algorithm `alg`;
#'   `ODEGradientDescent()` is [Euler()] (it needs `dt`),
#'   `RKChebyshevDescent()` [ROCK2()], `RKAccelerated()` [Tsit5()] and
#'   `HighOrderDescent()` [Vern7()].
#' * `DAEOptimizer(alg)` handles equality constraints `cons(u) = 0` (the
#'   values of `lcons` and `ucons` are not used, only their number `m`): the
#'   last `m` equations of the gradient flow are replaced by the
#'   constraints, which makes the last `m` unknowns algebraic. With a mass
#'   matrix ODE algorithm the system has the mass matrix
#'   `diag(c(rep(1, n - m), rep(0, m)))` and is integrated to steady state;
#'   `DAEMassMatrix()` is [Rodas5P()] with finite difference Jacobians. With
#'   a fully implicit DAE algorithm (such as [DFBDF()]) the residual
#'   `du + gradient(u)` (with the constraints in the last rows) is
#'   integrated over `[0, 10]`.
#'
#' The gradient of the problem is used, finite differences when it has
#' none. The return code is `"Success"` for `ODEOptimizer` and that of the
#' integration for `DAEOptimizer`; the statistics of `ODEOptimizer` give the
#' accepted steps as iterations and gradient evaluations, the right hand
#' side evaluations as `fevals`. A callback is called after every step.
#'
#' @param alg An ODE algorithm (for `DAEOptimizer`, a mass matrix ODE
#'   algorithm or a DAE algorithm).
#' @return An optimization algorithm. The ODE algorithm is `$alg`.
#' @examples
#' f <- function(u, p) sum((u - p)^2)
#' g <- function(u, p) 2 * (u - p)
#' prob <- OptimizationProblem(f, c(2, 2), c(1, 1), grad = g)
#' solve(prob, RKAccelerated())$u
#' solve(prob, ODEGradientDescent(), dt = 0.01)$u
#' @name ode_optimizers
#' @export
ODEOptimizer <- function(alg) {
  if (!inherits(alg, "ODEAlgorithm")) stop("`alg` must be an ODE algorithm", call. = FALSE)
  .odeopt_alg("ODEOptimizer", alg)
}
#' @rdname ode_optimizers
#' @export
ODEGradientDescent <- function() .odeopt_alg("ODEOptimizer", Euler())
#' @rdname ode_optimizers
#' @export
RKChebyshevDescent <- function() .odeopt_alg("ODEOptimizer", ROCK2())
#' @rdname ode_optimizers
#' @export
RKAccelerated <- function() .odeopt_alg("ODEOptimizer", Tsit5())
#' @rdname ode_optimizers
#' @export
HighOrderDescent <- function() .odeopt_alg("ODEOptimizer", Vern7())
#' @rdname ode_optimizers
#' @export
DAEOptimizer <- function(alg) {
  if (!inherits(alg, "ODEAlgorithm")) stop("`alg` must be an ODE or DAE algorithm", call. = FALSE)
  .odeopt_alg("DAEOptimizer", alg, dae = TRUE)
}
#' @rdname ode_optimizers
#' @export
DAEMassMatrix <- function() .odeopt_alg("DAEOptimizer", Rodas5P(autodiff = AutoFiniteDiff()), dae = TRUE)

#' @export
print.ODEOptimizer <- function(x, ...) {
  cat(x$name, "(", x$alg$name, "())\n", sep = "")
  invisible(x)
}
#' @export
print.DAEOptimizer <- print.ODEOptimizer

.odeopt_solve <- function(prob, alg, opts) {
  fs <- .opt_fns(prob)
  dae <- alg$name == "DAEOptimizer"
  if (!is.null(prob$lb) || !is.null(prob$ub))
    stop(IncompatibleOptimizerError(paste0("The algorithm ", alg$name, " does not support box constraints. ",
                                           "Either remove the `lb` or `ub` bounds passed to ",
                                           "`OptimizationProblem` or use a different algorithm.")))
  if (!dae && is.function(prob$f$cons))
    stop(IncompatibleOptimizerError(paste0("The algorithm ", alg$name, " does not support constraints. Either ",
                                           "remove the `cons` function passed to `OptimizationFunction` or use ",
                                           "a different algorithm.")))
  u0 <- as.numeric(prob$u0)
  p <- prob$p
  n <- length(u0)
  # the objective at every step, for the callback of solve(); the state is
  # unchanged, so the integrator keeps its caches and step size controller
  cb <- DiscreteCallback(function(u, t, integ) TRUE, function(integ) {
    fs$f(integ$u)
    u_modified(integ, FALSE)
  }, save_positions = c(FALSE, FALSE))
  sopts <- list(callback = cb)
  if (!is.null(opts$maxiters)) sopts$maxiters <- .opt_maxiters(opts$maxiters)
  if (!is.null(opts$dt)) sopts$dt <- opts$dt
  if (!dae) {
    rhs <- function(u, p, t) -fs$g(u)
    res <- .odeopt_dynamic_ss(ODEProblem(rhs, u0 = u0, tspan = c(0, .odeopt_tf), p = p), alg$alg, sopts)
    os <- res$sol$stats
    sol <- .opt_solution(prob, alg, res$u, fs$f(res$u), "Success", res$sol)
    return(.opt_set_fields(sol, stats = OptimizationStats(iterations = os$naccept, fevals = os$nf,
                                                          gevals = os$naccept)))
  }
  if (!is.function(prob$f$cons))
    stop("DAEOptimizer requires constraints. Please provide a function with `cons` defined.", call. = FALSE)
  m <- length(prob$ucons)
  if (m > n)
    stop("DAEOptimizer with mass matrix method requires the number of constraints to be less than or equal ",
         "to the number of variables.", call. = FALSE)
  rows <- seq_len(m) + (n - m)
  cons <- function(u) as.numeric(prob$f$cons(u, p))
  if (isTRUE(alg$alg$dae)) {
    # the fully implicit form, integrated over [0, 10] from du = 0. The
    # residual is that of the gradient flow du = -grad f, as in the mass
    # matrix form; OptimizationODE 0.1.12 writes du - grad f
    # (solve_dae_implicit), which climbs the objective and diverges
    resid <- function(du, u, p, t) {
      r <- du + fs$g(u)
      r[rows] <- cons(u)
      r
    }
    dprob <- DAEProblem(resid, du0 = numeric(n), u0 = u0, tspan = c(0, 10), p = p)
    osol <- do.call(solve, c(list(dprob, alg$alg), sopts, list(initializealg = ShampineCollocationInit())))
    u <- utils::tail(osol$u, 1L)[[1L]][seq_len(n)]
    retcode <- osol$retcode
  } else {
    rhs <- function(u, p, t) {
      du <- -fs$g(u)
      du[rows] <- cons(u)
      du
    }
    M <- diag(c(rep(1, n - m), rep(0, m)), n)
    oprob <- ODEProblem(ODEFunction(rhs, mass_matrix = M), u0 = u0, tspan = c(0, .odeopt_tf), p = p)
    sopts$initializealg <- BrownFullBasicInit()
    res <- .odeopt_dynamic_ss(oprob, alg$alg, sopts)
    u <- res$u[seq_len(n)]
    retcode <- res$retcode
    osol <- res$sol
  }
  .opt_solution(prob, alg, u, fs$f(u), retcode, osol)
}

# the end of the integration: Julia integrates over [0, Inf], which the R
# integrator does not take
.odeopt_tf <- 0x1.7e43c8800759cp+996

# DynamicSS of SteadyStateDiffEq with its default termination condition
.odeopt_dynamic_ss <- function(oprob, alg, sopts, abstol = 1e-8, reltol = 1e-6) {
  term <- DiscreteCallback(function(u, t, integ) {
    du <- get_du(integ)
    dn <- max(abs(du))
    isTRUE(dn <= abstol) || isTRUE(dn <= reltol * max(abs(du + integ$u)))
  }, function(integ) {
    terminate(integ)
    u_modified(integ, FALSE)
  }, save_positions = c(TRUE, FALSE))
  sopts$callback <- CallbackSet(term, sopts$callback)
  osol <- do.call(solve, c(list(oprob, alg), sopts,
                           list(abstol = abstol, reltol = reltol, save_everystep = FALSE, save_end = TRUE)))
  # termination_condition_result: only a terminated integration succeeded
  list(u = utils::tail(osol$u, 1L)[[1L]], t = utils::tail(osol$t, 1L),
       retcode = if (identical(osol$retcode, "Terminated")) "Success" else "Failure", sol = osol)
}
