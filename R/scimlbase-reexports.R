# SciMLBase, DiffEqBase and StochasticDiffEqCore names that StochasticDiffEq
# reexports: the solution types of SciMLBase (src/solutions/*.jl), the PDE
# wrappers and discretization hooks (src/problems/pde_problems.jl,
# src/solutions/pde_solutions.jl, src/SciMLBase.jl: discretize and
# symbolic_discretize), ConvexOptimizationProblem
# (src/problems/optimization_problems.jl), HomotopyNonlinearFunction and
# ODEInputFunction (src/scimlfunctions.jl), EnsembleContext and
# EnsembleTestSolution (src/ensemble/), DEVerbosity (DiffEqBase
# src/verbosity.jl), and from StochasticDiffEqCore SDEIntegrator,
# SDEOptions, concrete_prob (src/solve.jl), alg_cache and the integrator
# hooks initialize!, finalize! and addsteps!.

#' SciMLBase solution types
#'
#' Constructors of the solution types of SciMLBase. `DAESolution()` and
#' `RODESolution()` build a solution from saved times and states like
#' [ODESolution()] (interpolated linearly); a `RODESolution` also carries the
#' noise path `W` (a list with the times `t` and the matrix `W`).
#' `NonlinearSolution()`, `SteadyStateSolution()` and
#' `OptimizationSolution()` hold the solution `u` of a nonlinear, steady
#' state or optimization problem with its return code, residual and
#' statistics (and the objective value); `LinearSolution()` the solution of a
#' linear system (as [build_linear_solution()]). `PDETimeSeriesSolution(sol,
#' metadata)` and `PDENoTimeSolution(sol, metadata)` wrap the solution of a
#' discretized PDE with the discretization metadata, time dependent or not;
#' the wrapped solution stays in `original_sol`.
#'
#' @param t,u Saved times and states, or the solution value.
#' @param prob The problem.
#' @param alg The algorithm (any object with a `name`).
#' @param retcode The return code, see [ReturnCode].
#' @param W The noise path of a random or stochastic solution.
#' @param resid The residual.
#' @param stats Solver statistics.
#' @param objective The objective value.
#' @param original The result of the underlying optimizer.
#' @param sol A solution.
#' @param metadata The discretization metadata.
#' @return A solution object of the named class.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
#' sol <- RODESolution(c(0, 0.5, 1), rbind(1, 0.6, 0.37), prob,
#'   W = list(t = c(0, 0.5, 1), W = rbind(0, 0.1, 0.3)))
#' sol(0.25)
#' NonlinearSolution(NonlinearProblem(function(u, p) u^2 - 2, 1), list(name = "manual"), sqrt(2))$u
#' @name scimlbase_solutions
NULL

#' @rdname scimlbase_solutions
#' @export
DAESolution <- function(t, u, prob, alg = list(name = "DAESolution"), retcode = "Success") {
  sol <- ODESolution(t, u, prob, alg, retcode)
  class(sol) <- c("DAESolution", class(sol))
  sol
}

#' @rdname scimlbase_solutions
#' @export
RODESolution <- function(t, u, prob, alg = list(name = "RODESolution"), W = NULL, retcode = "Success") {
  sol <- ODESolution(t, u, prob, alg, retcode)
  if (!is.null(W)) {
    env <- .sol_env(sol)
    env$W <- W
  }
  .as_rode_solution(sol)
}

#' @rdname scimlbase_solutions
#' @export
NonlinearSolution <- function(prob, alg, u, retcode = "Success", resid = NULL, stats = list())
  .nl_new_solution("NonlinearSolution", prob, alg, u, retcode, resid = resid, stats = stats)

#' @rdname scimlbase_solutions
#' @export
SteadyStateSolution <- function(prob, alg, u, retcode = "Success", resid = NULL, stats = list())
  .nl_new_solution(c("SteadyStateSolution", "NonlinearSolution"), prob, alg, u, retcode, resid = resid, stats = stats)

#' @rdname scimlbase_solutions
#' @export
OptimizationSolution <- function(prob, alg, u, objective, retcode = "Success", stats = list(), original = NULL)
  .nl_new_solution(c("OptimizationSolution", "NonlinearSolution"), prob, alg, u, retcode, stats = stats,
                   extra = list(objective = objective, original = original))

#' @rdname scimlbase_solutions
#' @export
LinearSolution <- function(u, resid = NULL, alg = list(name = "LinearSolution"), retcode = "Success", stats = list())
  build_linear_solution(alg, u, resid, retcode, stats)

.pde_solution <- function(cls, sol, metadata) {
  out <- sol
  attr(out, "original_sol") <- sol
  attr(out, "disc_data") <- metadata
  class(out) <- c(cls, "PDESolution", class(sol))
  out
}

#' @rdname scimlbase_solutions
#' @export
PDETimeSeriesSolution <- function(sol, metadata) .pde_solution("PDETimeSeriesSolution", sol, metadata)

#' @rdname scimlbase_solutions
#' @export
PDENoTimeSolution <- function(sol, metadata) .pde_solution("PDENoTimeSolution", sol, metadata)

#' PDE problems and discretization
#'
#' `PDEProblem(prob, extrapolation, space)` wraps a problem obtained from a
#' PDE with the extrapolation and spatial domain data a discretization
#' produced. `discretize(sys, discretizer, ...)` turns a [PDESystem()] into a
#' solver ready problem and `symbolic_discretize(sys, discretizer, ...)` into
#' its discretized system; both are generic functions that dispatch on the
#' class of `discretizer`, for discretization packages to extend (as the
#' methods of MethodOfLines or NeuralPDE extend SciMLBase). modelcompR itself
#' has no PDE discretizer, so the default methods stop with an error.
#'
#' @param prob The problem of the discretized PDE.
#' @param extrapolation,space Extrapolation and spatial domain data.
#' @param sys A PDE system.
#' @param discretizer A discretization algorithm.
#' @param ... Options of the discretizer.
#' @return A `PDEProblem`, or what the method of the discretizer returns.
#' @examples
#' prob <- PDEProblem(ODEProblem(function(u, p, t) -u, c(1, 2), c(0, 1)), space = list(x = c(0, 1)))
#' class(prob)
#' discretize.MyGrid <- function(sys, discretizer, ...) "discretized"
#' discretize(NULL, structure(list(), class = "MyGrid"))
#' @name pde_discretization
NULL

#' @rdname pde_discretization
#' @export
PDEProblem <- function(prob, extrapolation = NULL, space = NULL)
  structure(list(prob = prob, extrapolation = extrapolation, space = space), class = c("PDEProblem", "SciMLProblem"))

#' @rdname pde_discretization
#' @export
discretize <- function(sys, discretizer, ...) UseMethod("discretize", discretizer)

#' @rdname pde_discretization
#' @export
discretize.default <- function(sys, discretizer, ...)
  stop("no discretize() method for a discretizer of class ", class(discretizer)[1L], call. = FALSE)

#' @rdname pde_discretization
#' @export
symbolic_discretize <- function(sys, discretizer, ...) UseMethod("symbolic_discretize", discretizer)

#' @rdname pde_discretization
#' @export
symbolic_discretize.default <- function(sys, discretizer, ...)
  stop("no symbolic_discretize() method for a discretizer of class ", class(discretizer)[1L], call. = FALSE)

#' Further SciMLBase problem and function types
#'
#' `ConvexOptimizationProblem(f, u0, p, constraints, lb, ub, int, sense)` is
#' the experimental convex optimization problem of SciMLBase: an objective
#' `f(u, p)` (or an [OptimizationFunction()]), a start, conic `constraints`,
#' optional bounds (both or neither), a mask `int` of integer variables and
#' the sense `"MinSense"` or `"MaxSense"`. It is solved with
#' [ConvexMOI()] and [ConeConstraint()] constraints. `HomotopyNonlinearFunction(f, polynomialize, unpolynomialize,
#' denominator, polynomial_jac)` is the polynomial form of a nonlinear system
#' for homotopy continuation, with the maps between the original and the
#' polynomial variables. `ODEInputFunction(f, ...)` is the right hand side
#' `f(x, u, p, t)` of a system with states `x` and inputs `u`, with optional
#' derivatives (`jac` with respect to `x`, `controljac` with respect to
#' `u`).
#'
#' @param f The objective, residual or right hand side.
#' @param u0 The start.
#' @param p The parameters.
#' @param constraints Constraints, or `NULL`.
#' @param lb,ub Bounds, or `NULL`.
#' @param int A logical mask of integer variables, or `NULL`.
#' @param sense `"MinSense"` or `"MaxSense"`, or the objects [MinSense] and
#'   [MaxSense].
#' @param polynomialize,unpolynomialize Maps from the original variables to
#'   the polynomial ones and back.
#' @param denominator The denominators of the polynomial form.
#' @param polynomial_jac The Jacobian of the polynomial form.
#' @param mass_matrix,analytic,tgrad,jac,controljac,jvp,vjp,jac_prototype,controljac_prototype,sparsity,paramjac,observed,sys
#'   Optional parts of the function.
#' @param ... Further fields.
#' @return An object of the named class.
#' @examples
#' cp <- ConvexOptimizationProblem(function(u, p) sum(u^2), c(1, 2), lb = c(0, 0), ub = c(3, 3))
#' cp$sense
#' hf <- HomotopyNonlinearFunction(function(u, p) u^2 - 2)
#' hf$f(2, NULL)
#' inf <- ODEInputFunction(function(x, u, p, t) -x + u)
#' inf$f(1, 0.5, NULL, 0)
#' @name scimlbase_more_types
NULL

#' @rdname scimlbase_more_types
#' @export
ConvexOptimizationProblem <- function(f, u0, p = NullParameters(), constraints = NULL, lb = NULL, ub = NULL,
                                      int = NULL, sense = c("MinSense", "MaxSense"), ...) {
  if (xor(is.null(lb), is.null(ub))) stop("If any of `lb` or `ub` is provided, both must be provided.", call. = FALSE)
  sense <- match.arg(sense)
  if (!inherits(f, "OptimizationFunction"))
    f <- structure(list(f = f), class = c("OptimizationFunction", "SciMLFunction"))
  structure(list(f = f, u0 = u0, p = p, constraints = constraints, lb = lb, ub = ub, int = int, sense = sense,
                 kwargs = list(...)),
            class = c("ConvexOptimizationProblem", "SciMLProblem"))
}

#' @rdname scimlbase_more_types
#' @export
HomotopyNonlinearFunction <- function(f, polynomialize = function(u, p) u, unpolynomialize = function(u, p) list(u),
                                      denominator = NULL, polynomial_jac = NULL, ...) {
  if (inherits(f, "HomotopyNonlinearFunction")) return(f)
  structure(list(f = if (is.function(f)) f else f$f, polynomialize = polynomialize, unpolynomialize = unpolynomialize,
                 denominator = denominator, jac = polynomial_jac, ...),
            class = c("HomotopyNonlinearFunction", "NonlinearFunction", "SciMLFunction"))
}

#' @rdname scimlbase_more_types
#' @export
ODEInputFunction <- function(f, mass_matrix = NULL, analytic = NULL, tgrad = NULL, jac = NULL, controljac = NULL,
                             jvp = NULL, vjp = NULL, jac_prototype = NULL, controljac_prototype = NULL,
                             sparsity = jac_prototype, paramjac = NULL, observed = NULL, sys = NULL) {
  if (inherits(f, "ODEInputFunction")) return(f)
  structure(list(f = f, mass_matrix = mass_matrix, analytic = analytic, tgrad = tgrad, jac = jac,
                 controljac = controljac, jvp = jvp, vjp = vjp, jac_prototype = jac_prototype,
                 controljac_prototype = controljac_prototype, sparsity = sparsity, paramjac = paramjac,
                 observed = observed, sys = sys),
            class = c("ODEInputFunction", "SciMLFunction"))
}

#' Ensemble context and test solutions
#'
#' `EnsembleContext(sim_id, repeat_, worker_id, sim_seed, rng, master_rng)`
#' describes one trajectory of an ensemble solve: its index, the number of
#' the repetition (field `repeat`), the worker, the seed of the trajectory
#' and the random number generators. `EnsembleTestSolution()` is an ensemble
#' solution of a problem with an analytic solution together with its error
#' summaries: the errors of each trajectory, the weak errors of the mean and
#' the mean and median errors (as [calculate_ensemble_errors()] computes
#' them); `sim` may be an [EnsembleSolution()] whose trajectories, time and
#' convergence it takes.
#'
#' @param sim_id,repeat_,worker_id The trajectory, repetition and worker.
#' @param sim_seed The seed of the trajectory.
#' @param rng,master_rng Random number generators, or `NULL`.
#' @param sim An ensemble solution or the list of trajectories.
#' @param errors,weak_errors,error_means,error_medians Named error summaries.
#' @param elapsedTime The run time.
#' @param converged Whether the ensemble converged.
#' @return An `EnsembleContext` or an `EnsembleTestSolution`.
#' @examples
#' ctx <- EnsembleContext(3, 1, 1, 42)
#' ctx$sim_id
#' ts <- EnsembleTestSolution(list(1, 2), list(final = c(0.1, 0.2)), list(final = 0.05),
#'                            list(final = 0.15), list(final = 0.15), 0.1, TRUE)
#' ts$error_means$final
#' @name ensemble_context
NULL

#' @rdname ensemble_context
#' @export
EnsembleContext <- function(sim_id, repeat_ = 1L, worker_id = 1L, sim_seed = NULL, rng = NULL, master_rng = NULL)
  structure(list(sim_id = as.integer(sim_id), `repeat` = as.integer(repeat_), worker_id = as.integer(worker_id),
                 sim_seed = sim_seed, rng = rng, master_rng = master_rng), class = "EnsembleContext")

#' @rdname ensemble_context
#' @export
EnsembleTestSolution <- function(sim, errors, weak_errors, error_means, error_medians, elapsedTime = 0,
                                 converged = TRUE) {
  if (inherits(sim, "EnsembleSolution")) {
    elapsedTime <- sim$elapsed
    converged <- sim$converged
    sim <- sim$u
  }
  structure(list(u = sim, errors = errors, weak_errors = weak_errors, error_means = error_means,
                 error_medians = error_medians, elapsed = elapsedTime, converged = converged),
            class = c("EnsembleTestSolution", "EnsembleSolution"))
}

.DE_VERBOSITY <- local({
  toggles <- c("dt_NaN", "init_NaN", "dense_output_saveat", "max_iters", "dt_min_unstable", "instability",
               "symbolic_diagnostic", "newton_convergence", "step_rejected", "step_accepted", "convergence_limit",
               "alg_switch", "stiff_detection", "mismatched_input_output_type", "jacobian_update",
               "w_factorization", "newton_iterations", "rosenbrock_no_differential_states", "shampine_dt",
               "unlimited_dt", "dt_epsilon", "stability_check", "near_singular", "stage_limiter_unused",
               "step_limiter_unused", "sensitivity_vjp_choice", "noise_evaluation", "discontinuity_tracking",
               "delay_evaluation", "constrained_step", "residual_control", "neutral_delay", "state_dependent_delay")
  S <- "Silent"; W <- "WarnLevel"; I <- "InfoLevel"; E <- "ErrorLevel"
  presets <- list(
    None = rep(S, 33),
    Minimal = c(W, W, S, W, W, W, S, W, S, S, S, S, S, S, S, S, S, W, S, W, W, S, W, E, W, S, S, S, S, S, S, S, W),
    Standard = c(W, W, W, W, W, W, W, S, S, S, S, S, S, W, S, S, S, W, S, W, W, S, S, E, W, S, S, S, S, S, S, S, S),
    Detailed = c(W, W, I, W, W, W, W, W, S, S, I, I, S, W, I, I, I, W, I, W, I, I, W, E, W, W, I, I, I, I, I, I, W),
    All = c(W, W, I, W, W, W, W, W, I, I, I, I, I, I, I, I, I, W, I, W, I, I, W, E, W, W, I, I, I, I, I, I, I))
  sub <- c(None = "None", Minimal = "Minimal", Standard = "Minimal", Detailed = "Detailed", All = "All")
  list(toggles = toggles, presets = lapply(presets, stats::setNames, toggles), sub = sub)
})

#' Verbosity of differential equation solvers
#'
#' `DEVerbosity(...)` collects the levels (`"Silent"`, `"DebugLevel"`,
#' `"InfoLevel"`, `"WarnLevel"` or `"ErrorLevel"`) of the diagnostic
#' messages of the differential equation solvers, one per toggle (such as
#' `dt_NaN`, `max_iters`, `alg_switch` or the SDE toggle
#' `noise_evaluation`), with the presets of DiffEqBase: `"Standard"` (the
#' default), `"None"`, `"Minimal"`, `"Detailed"` and `"All"`; `TRUE` and
#' `FALSE` stand for `"Standard"` and `"None"`. The nested
#' `linear_verbosity` and `nonlinear_verbosity` are preset names for the
#' linear and nonlinear solvers. Toggles given in `...` override the preset.
#'
#' @param ... Toggles, such as `alg_switch = "InfoLevel"`.
#' @param preset A preset name or a logical.
#' @param linear_verbosity,nonlinear_verbosity Presets of the inner solvers;
#'   by default those of `preset`.
#' @return An object of class `DEVerbosity`, a named character vector with
#'   the inner presets as attributes.
#' @examples
#' v <- DEVerbosity(alg_switch = "InfoLevel")
#' v[["alg_switch"]]
#' DEVerbosity(preset = "None")[["dt_NaN"]]
#' @export
DEVerbosity <- function(..., preset = "Standard", linear_verbosity = NULL, nonlinear_verbosity = NULL) {
  if (isTRUE(preset)) preset <- "Standard"
  if (identical(preset, FALSE)) preset <- "None"
  V <- .DE_VERBOSITY
  lv <- V$presets[[preset]]
  if (is.null(lv)) stop("unknown verbosity preset ", preset, call. = FALSE)
  over <- unlist(list(...))
  bad <- setdiff(names(over), V$toggles)
  if (length(bad)) stop("unknown verbosity toggle(s): ", paste(bad, collapse = ", "), call. = FALSE)
  badlv <- setdiff(over, c("Silent", "DebugLevel", "InfoLevel", "WarnLevel", "ErrorLevel"))
  if (length(badlv)) stop("unknown verbosity level(s): ", paste(badlv, collapse = ", "), call. = FALSE)
  lv[names(over)] <- over
  structure(lv, linear_verbosity = linear_verbosity %||% V$sub[[preset]],
            nonlinear_verbosity = nonlinear_verbosity %||% V$sub[[preset]], class = "DEVerbosity")
}

#' @export
print.DEVerbosity <- function(x, ...) {
  cat("DEVerbosity (linear ", attr(x, "linear_verbosity"), ", nonlinear ", attr(x, "nonlinear_verbosity"), ")\n",
      sep = "")
  print(stats::setNames(as.character(x), names(x)), quote = FALSE)
  invisible(x)
}

#' SDE integrator pieces
#'
#' `SDEIntegrator` is the class of the integrators [init()] returns for SDE
#' problems (they are also `ODEIntegrator` objects, with the same stepping
#' interface). `SDEOptions(...)` collects solver options of an SDE solve with
#' the defaults of StochasticDiffEq, for `do.call(solve, c(list(prob, alg),
#' opts))`. `concrete_prob(prob)` is the differential equation problem a
#' problem carries: the problem itself, or the wrapped problem of a
#' [JumpProblem()]. `alg_cache(alg, prob, u, t)` builds the cache an
#' algorithm starts a solve with (the state its stepping keeps). The
#' integrator hooks `initialize_callback(cb, u, t, integrator)` and
#' `finalize(cb, u, t, integrator)` run the `initialize` and `finalize`
#' functions of a callback or callback set (`initialize!` and `finalize!`
#' upstream), and `addsteps(integrator)` computes the extra interpolation
#' stages of the current step, of which SDE integrators, interpolated
#' linearly, have none.
#'
#' @param abstol,reltol Tolerances.
#' @param qmax,qmin,gamma Step size controller limits and safety factor.
#' @param delta The weight of the drift error estimate.
#' @param maxiters The largest number of steps.
#' @param dtmax The largest step size, or `NULL`.
#' @param adaptive Adaptive stepping, or `NULL` for the method's default.
#' @param save_everystep Save every step.
#' @param ... Further options.
#' @param prob A problem.
#' @param alg An algorithm.
#' @param u,t The state and time to start from.
#' @param cb A callback or callback set.
#' @param integrator An integrator.
#' @return `SDEIntegrator` is a class name; `SDEOptions()` a list of class
#'   `SDEOptions`; `concrete_prob()` a problem; `alg_cache()` the cache;
#'   the hooks return the integrator invisibly.
#' @examples
#' prob <- SDEProblem(function(u, p, t) -u, function(u, p, t) 0.1 * u, 1, c(0, 1))
#' integ <- init(prob, SOSRI(), seed = 1)
#' inherits(integ, SDEIntegrator)
#' opts <- SDEOptions(abstol = 1e-3)
#' sol <- do.call(solve, c(list(prob, SOSRI()), opts, seed = 1))
#' identical(concrete_prob(prob), prob)
#' names(alg_cache(Tsit5(), ODEProblem(function(u, p, t) -u, 1, c(0, 1))))
#' class(alg_cache(ImplicitEM(), prob))
#' @name sde_integrator_pieces
NULL

#' @rdname sde_integrator_pieces
#' @export
SDEIntegrator <- "SDEIntegrator"

#' @rdname sde_integrator_pieces
#' @export
SDEOptions <- function(abstol = 1e-2, reltol = 1e-2, qmax = 9 / 8, qmin = 1 / 5, gamma = 9 / 10, delta = NULL,
                       maxiters = 1e7, dtmax = NULL, adaptive = NULL, save_everystep = TRUE, ...) {
  o <- list(abstol = abstol, reltol = reltol, qmax = qmax, qmin = qmin, gamma = gamma, delta = delta,
            maxiters = maxiters, dtmax = dtmax, adaptive = adaptive, save_everystep = save_everystep, ...)
  structure(o[!vapply(o, is.null, TRUE)], class = "SDEOptions")
}

#' @rdname sde_integrator_pieces
#' @export
concrete_prob <- function(prob) if (inherits(prob, "JumpProblem")) prob$prob else prob

#' @rdname sde_integrator_pieces
#' @export
alg_cache <- function(alg, prob, u = prob$u0, t = prob$tspan[1L]) {
  prob <- concrete_prob(prob)
  if (is.function(alg$bind)) alg <- alg$bind(prob, list())
  init <- alg$init %||% if (inherits(alg, "SDEAlgorithm")) function(ctx, t, u) list() else .fsal_init
  init(list(f = prob$f$f, p = prob$p, iter = 0L), t, u)
}

.callback_list <- function(cb) {
  s <- .as_callback_set(cb)
  c(s$continuous, s$discrete)
}

#' @rdname sde_integrator_pieces
#' @export
initialize_callback <- function(cb, u, t, integrator) {
  for (x in .callback_list(cb)) if (is.function(x$initialize)) x$initialize(integrator)
  invisible(integrator)
}

#' @rdname sde_integrator_pieces
#' @export
finalize <- function(cb, u, t, integrator) {
  for (x in .callback_list(cb)) if (is.function(x$finalize)) x$finalize(integrator)
  invisible(integrator)
}

#' @rdname sde_integrator_pieces
#' @export
addsteps <- function(integrator, ...) invisible(integrator)
