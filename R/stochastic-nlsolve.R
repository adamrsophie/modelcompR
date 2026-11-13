# The nonlinear solver choices of OrdinaryDiffEqNonlinearSolve (src/type.jl:
# NLFunctional, NLAnderson, NLNewton, NonlinearSolveAlg) that
# StochasticDiffEq reexports for the `nlsolve` argument of its implicit
# methods, and TauLeapingDrift of StochasticDiffEqCore
# (src/integrators/integrator_utils.jl), the drift of a regular jump that
# the implicit tau leaping methods hand to their nonlinear solver.

#' Nonlinear solver choices of implicit methods
#'
#' Option objects for the `nlsolve` argument of implicit methods, with the
#' defaults of OrdinaryDiffEqNonlinearSolve. `NLNewton()` is the quasi
#' Newton iteration with the matrix `W = M / (gamma dt) - J`, kept across
#' steps while the step size changes by less than `new_W_dt_cutoff`;
#' `NLFunctional()` the fixed point iteration `z = g(z)` without Jacobian;
#' `NLAnderson()` the fixed point iteration with Anderson acceleration over
#' the last `max_history` iterates from iteration `aa_start`; and
#' `NonlinearSolveAlg(alg)` hands the stage equations to the nonlinear
#' solver `alg` (for example [NewtonRaphson()]). Each iteration stops when the
#' scaled increment, multiplied by `eta = theta / (1 - theta)`, is below
#' `kappa`, and fails after `max_iter` iterations; a rate `theta` below
#' `fast_convergence_cutoff` counts as fast convergence.
#'
#' @param kappa The relative tolerance on the increment (\ifelse{latex}{\verb{kappa}}{\verb{κ}} upstream).
#' @param max_iter The largest number of iterations per solve.
#' @param fast_convergence_cutoff The convergence rate threshold.
#' @param new_W_dt_cutoff The relative change of `gamma dt` above which `W`
#'   is formed again.
#' @param always_new Form `W` at every solve.
#' @param check_div Stop early when the iteration diverges.
#' @param relax A relaxation in `[0, 1)` damping the Newton update, or
#'   `NULL`.
#' @param max_history,aa_start,droptol The Anderson acceleration memory,
#'   first accelerated iteration and drop tolerance.
#' @param alg A nonlinear solver algorithm.
#' @param predictor A function `predictor(uprev, p, t, dt)` giving the start
#'   of each stage solve, or `NULL`.
#' @return A nonlinear solver choice (class `AbstractNLSolverAlgorithm`).
#' @examples
#' NLNewton(kappa = 1e-3)$kappa
#' NLAnderson()$max_history
#' NonlinearSolveAlg()$alg$name
#' @name nl_solver_choices
NULL

.nlsolver_alg <- function(name, ...) structure(list(name = name, ...), class = c(name, "AbstractNLSolverAlgorithm"))

#' @export
print.AbstractNLSolverAlgorithm <- function(x, ...) {
  cat(x$name, "(kappa = ", format(x$kappa), ", max_iter = ", x$max_iter, ")\n", sep = "")
  invisible(x)
}

.check_predictor <- function(predictor) {
  if (!is.null(predictor) && !is.function(predictor))
    stop("predictor is a function predictor(uprev, p, t, dt) or NULL", call. = FALSE)
  predictor
}

#' @rdname nl_solver_choices
#' @export
NLFunctional <- function(kappa = 1 / 100, max_iter = 10L, fast_convergence_cutoff = 1 / 5, predictor = NULL)
  .nlsolver_alg("NLFunctional", kappa = kappa, fast_convergence_cutoff = fast_convergence_cutoff,
                max_iter = as.integer(max_iter), predictor = .check_predictor(predictor))

#' @rdname nl_solver_choices
#' @export
NLAnderson <- function(kappa = 1 / 100, max_iter = 10L, max_history = 5L, aa_start = 1L, droptol = NULL,
                       fast_convergence_cutoff = 1 / 5, predictor = NULL)
  .nlsolver_alg("NLAnderson", kappa = kappa, fast_convergence_cutoff = fast_convergence_cutoff,
                max_iter = as.integer(max_iter), max_history = as.integer(max_history),
                aa_start = as.integer(aa_start), droptol = droptol, predictor = .check_predictor(predictor))

#' @rdname nl_solver_choices
#' @export
NLNewton <- function(kappa = 1 / 100, max_iter = 10L, fast_convergence_cutoff = 1 / 5, new_W_dt_cutoff = 1 / 5,
                     always_new = FALSE, check_div = TRUE, relax = NULL, predictor = NULL) {
  if (is.numeric(relax) && !(relax >= 0 && relax < 1))
    stop(sprintf("The relaxation parameter must be in [0, 1), got `relax = %s`", format(relax)), call. = FALSE)
  .nlsolver_alg("NLNewton", kappa = kappa, max_iter = as.integer(max_iter),
                fast_convergence_cutoff = fast_convergence_cutoff, new_W_dt_cutoff = new_W_dt_cutoff,
                always_new = always_new, check_div = check_div, relax = relax,
                predictor = .check_predictor(predictor))
}

#' @rdname nl_solver_choices
#' @export
NonlinearSolveAlg <- function(alg = NewtonRaphson(), kappa = 1 / 100, max_iter = 10L,
                              fast_convergence_cutoff = 1 / 5, new_W_dt_cutoff = 1 / 5, always_new = FALSE,
                              check_div = TRUE, predictor = NULL)
  .nlsolver_alg("NonlinearSolveAlg", kappa = kappa, max_iter = as.integer(max_iter),
                fast_convergence_cutoff = fast_convergence_cutoff, new_W_dt_cutoff = new_W_dt_cutoff,
                always_new = always_new, check_div = check_div, alg = alg, predictor = .check_predictor(predictor))

#' Drift of a regular jump
#'
#' `TauLeapingDrift(c, rate)` is the function `drift(u, p, t) = c(u, p, t,
#' rate(u, p, t), NULL)`: the change of the state by the expected jump
#' counts, with `rate` the propensities and `c` the stoichiometry function of
#' a [RegularJump()]. The implicit tau leaping methods solve their stage
#' equations with it.
#'
#' @param c The stoichiometry function `c(u, p, t, counts, mark)`.
#' @param rate The propensity function `rate(u, p, t)`.
#' @param rate_cache Scratch space of the in place form; unused in R.
#' @return A function of `(u, p, t)` of class `TauLeapingDrift`.
#' @examples
#' rate <- function(u, p, t) c(0.5 * u[1], 0.1 * u[1])
#' cfun <- function(u, p, t, counts, mark) c(counts[1] - counts[2])
#' drift <- TauLeapingDrift(cfun, rate)
#' drift(10, NULL, 0)
#' @export
TauLeapingDrift <- function(c, rate, rate_cache = NULL)
  structure(function(u, p, t) c(u, p, t, rate(u, p, t), NULL), c = c, rate = rate,
            class = base::c("TauLeapingDrift", "function"))
