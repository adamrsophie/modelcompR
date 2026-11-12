# The SciMLBase.jl interface functions and types: interpolations
# (src/interpolation.jl), root finding options and algorithm interpretation
# (src/callbacks.jl, src/alg_traits.jl), clocks (src/clock.jl), aliasing
# specifiers (src/problems/*), statistics (DEStats, NLStats), the solution
# builders (src/solutions/*), the function wrappers of
# src/function_wrappers.jl, the traits of src/alg_traits.jl and
# src/scimlfunctions.jl, ensembles (src/ensemble/*), the error types of
# src/utils.jl, initialization data (src/initialization.jl), saved
# subsystems (src/solutions/save_idxs.jl) and parameter time series.

# ---------------------------------------------------------------------------
# Interpolations, root finding, interpretation

#' Interpolations, root finding options and algorithm interpretation
#'
#' `ConstantInterpolation(t, u)`, `LinearInterpolation(t, u)` and
#' `HermiteInterpolation(t, u, du)` return functions of time that interpolate
#' saved states (a matrix with one row per time or a list of states); they
#' are the interpolations solutions use. `BasicInterpolation(t, u, du)` is
#' Hermite interpolation when `du` is given and linear otherwise;
#' `SensitivityInterpolation(t, u)` the linear interpolation of
#' sensitivity solutions. `interp_summary(sol)` describes the interpolation
#' of a solution.
#'
#' `RootfindOpt` lists `NoRootFind`, `LeftRootFind` and `RightRootFind`, the
#' values of the `rootfind` argument of [ContinuousCallback()].
#' `AlgorithmInterpretation` lists `Ito` and `Stratonovich`, and
#' `alg_interpretation(alg)` gives the interpretation of an SDE algorithm.
#'
#' @param t,u,du Times, states and derivatives.
#' @param sol A solution.
#' @param alg An algorithm.
#' @return A function of time, a description, an option or the
#'   interpretation.
#' @examples
#' f <- HermiteInterpolation(c(0, 1), rbind(0, 1), rbind(0, 2))
#' f(0.5)
#' alg_interpretation(EulerHeun())
#' @name interpolations
NULL

.interp_rows <- function(u) if (is.list(u)) do.call(rbind, lapply(u, as.numeric)) else as.matrix(u)

.interp_fun <- function(t, U, step_value, kind) {
  t <- as.numeric(t)
  if (nrow(U) != length(t)) stop("give one state per time", call. = FALSE)
  f <- function(tt) {
    out <- t(vapply(tt, function(x) {
      k <- max(1L, min(length(t) - 1L, findInterval(x, t, rightmost.closed = TRUE)))
      if (length(t) == 1L) return(U[1L, ])
      step_value(k, x)
    }, numeric(ncol(U))))
    if (length(tt) == 1L) out[1L, ] else out
  }
  structure(f, class = c(kind, "DiffEqInterpolation", "function"))
}

#' @rdname interpolations
#' @export
ConstantInterpolation <- function(t, u) {
  U <- .interp_rows(u)
  .interp_fun(t, U, function(k, x) if (x >= t[k + 1L]) U[k + 1L, ] else U[k, ], "ConstantInterpolation")
}

#' @rdname interpolations
#' @export
LinearInterpolation <- function(t, u) {
  U <- .interp_rows(u)
  .interp_fun(t, U, function(k, x) {
    th <- (x - t[k]) / (t[k + 1L] - t[k])
    (1 - th) * U[k, ] + th * U[k + 1L, ]
  }, "LinearInterpolation")
}

#' @rdname interpolations
#' @export
HermiteInterpolation <- function(t, u, du) {
  U <- .interp_rows(u); DU <- .interp_rows(du)
  .interp_fun(t, U, function(k, x) {
    h <- t[k + 1L] - t[k]
    th <- (x - t[k]) / h
    (1 - th) * U[k, ] + th * U[k + 1L, ] +
      th * (th - 1) * ((1 - 2 * th) * (U[k + 1L, ] - U[k, ]) + (th - 1) * h * DU[k, ] + th * h * DU[k + 1L, ])
  }, "HermiteInterpolation")
}

#' @rdname interpolations
#' @export
BasicInterpolation <- function(t, u, du = NULL) if (is.null(du)) LinearInterpolation(t, u) else HermiteInterpolation(t, u, du)

#' @rdname interpolations
#' @export
SensitivityInterpolation <- function(t, u) LinearInterpolation(t, u)

#' @rdname interpolations
#' @export
interp_summary <- function(sol) {
  env <- .sol_env(sol)
  if (!isTRUE(env$dense) || is.null(env$steps)) return("1st order linear")
  alg <- env$alg$name %||% ""
  if (grepl("^Vern|^Tsit5|^DP5|^BS3", alg)) return(paste("specialized order interpolation of", alg))
  if (grepl("^Rodas|^Rosenbrock", alg)) return("specialized 3rd order \"free\" stiffness-aware interpolation")
  "3rd order Hermite"
}

#' @rdname interpolations
#' @export
NoRootFind <- "none"
#' @rdname interpolations
#' @export
LeftRootFind <- "left"
#' @rdname interpolations
#' @export
RightRootFind <- "right"
#' @rdname interpolations
#' @export
RootfindOpt <- list(NoRootFind = "none", LeftRootFind = "left", RightRootFind = "right")

#' @rdname interpolations
#' @export
AlgorithmInterpretation <- list(Ito = "Ito", Stratonovich = "Stratonovich")

#' @rdname interpolations
#' @export
alg_interpretation <- function(alg) alg$interpretation %||% "Ito"

# ---------------------------------------------------------------------------
# Clocks

#' Clocks of SciMLBase
#'
#' `PeriodicClock(dt, phase)` is the periodic clock ([Clock()]),
#' `Continuous()` (or [ContinuousClock()]) the continuous time domain,
#' `EventClock(name)` the clock of an event and `IndexedClock(clock, idx)`
#' the ticks `idx` of a clock. `isclock()`, `iscontinuous()`,
#' `issolverstepclock()`, `iseventclock()` and `is_discrete_time_domain()`
#' classify time domains. `first_clock_tick_time(clock, t0)` is the first tick
#' at or after `t0`; `canonicalize_indexed_clock(ic, sol)` returns the times of
#' the ticks of an indexed clock in a solution (`sol$t` for the continuous
#' clock).
#'
#' @param dt,phase Period and offset.
#' @param name Name of the event.
#' @param clock,ic A clock or an indexed clock.
#' @param idx Tick indices.
#' @param x A time domain.
#' @param t0 A time.
#' @param sol A solution.
#' @return A clock, a flag, a time or times.
#' @examples
#' c1 <- PeriodicClock(0.5)
#' first_clock_tick_time(c1, 0.2)
#' iscontinuous(Continuous())
#' @name scimlbase_clocks
NULL

#' @rdname scimlbase_clocks
#' @export
PeriodicClock <- function(dt, phase = 0) Clock(dt, phase)

#' @rdname scimlbase_clocks
#' @export
Continuous <- function() ContinuousClock()

#' @rdname scimlbase_clocks
#' @export
EventClock <- function(name = NULL) structure(list(type = "event", name = name), class = c("EventClock", "TimeDomain"))

#' @rdname scimlbase_clocks
#' @export
IndexedClock <- function(clock, idx) structure(list(clock = clock, idx = idx), class = "IndexedClock")

#' @rdname scimlbase_clocks
#' @export
isclock <- function(x) inherits(x, "Clock") || (is.list(x) && identical(x$type, "periodic"))

#' @rdname scimlbase_clocks
#' @export
iscontinuous <- function(x) inherits(x, "ContinuousClock")

#' @rdname scimlbase_clocks
#' @export
issolverstepclock <- function(x) inherits(x, "SolverStepClock")

#' @rdname scimlbase_clocks
#' @export
iseventclock <- function(x) inherits(x, "EventClock")

#' @rdname scimlbase_clocks
#' @export
is_discrete_time_domain <- function(x) isclock(x) || issolverstepclock(x) || iseventclock(x)

#' @rdname scimlbase_clocks
#' @export
first_clock_tick_time <- function(clock, t0) {
  dt <- clock$dt %||% stop("only periodic clocks have tick times", call. = FALSE)
  ph <- clock$phase %||% 0
  ph + ceiling((t0 - ph) / dt - 1e-12) * dt
}

#' @rdname scimlbase_clocks
#' @export
canonicalize_indexed_clock <- function(ic, sol) {
  c0 <- ic$clock
  ts <- if (iscontinuous(c0) || issolverstepclock(c0)) sol$t else {
    t0 <- sol$t[1L]; tf <- sol$t[length(sol$t)]
    seq(first_clock_tick_time(c0, t0), tf, by = c0$dt)
  }
  ts[ic$idx]
}

# ---------------------------------------------------------------------------
# Aliasing specifiers

#' Aliasing specifiers
#'
#' SciMLBase's aliasing options say whether a solver may reuse the arrays of
#' a problem (`alias_u0`, `alias_p`, `alias_f`, ...). R values are never
#' shared: a solver works on copies, so every specifier is accepted (as the
#' `alias` option of [solve()]) and has no effect on the results.
#' `ODEAliasSpecifier()`, `SDEAliasSpecifier()`, `DAEAliasSpecifier()`,
#' `DDEAliasSpecifier()`, `SDDEAliasSpecifier()`, `RODEAliasSpecifier()`,
#' `DiscreteAliasSpecifier()`, `ImplicitDiscreteAliasSpecifier()`,
#' `BVPAliasSpecifier()`, `NonlinearAliasSpecifier()`,
#' `SteadyStateAliasSpecifier()`, `OptimizationAliasSpecifier()`,
#' `IntegralAliasSpecifier()`, `AnalyticalAliasSpecifier()` and
#' `LinearAliasSpecifier()` record the flags given.
#'
#' @param ... Flags such as `alias_u0 = TRUE`.
#' @return An aliasing specifier.
#' @examples
#' ODEAliasSpecifier(alias_u0 = TRUE)
#' @name alias_specifiers
NULL

.alias_spec <- function(kind, ...) structure(list(...), class = c(kind, "AliasSpecifier"))

#' @rdname alias_specifiers
#' @export
ODEAliasSpecifier <- function(...) .alias_spec("ODEAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
SDEAliasSpecifier <- function(...) .alias_spec("SDEAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
DAEAliasSpecifier <- function(...) .alias_spec("DAEAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
DDEAliasSpecifier <- function(...) .alias_spec("DDEAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
SDDEAliasSpecifier <- function(...) .alias_spec("SDDEAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
RODEAliasSpecifier <- function(...) .alias_spec("RODEAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
DiscreteAliasSpecifier <- function(...) .alias_spec("DiscreteAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
ImplicitDiscreteAliasSpecifier <- function(...) .alias_spec("ImplicitDiscreteAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
BVPAliasSpecifier <- function(...) .alias_spec("BVPAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
NonlinearAliasSpecifier <- function(...) .alias_spec("NonlinearAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
SteadyStateAliasSpecifier <- function(...) .alias_spec("SteadyStateAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
OptimizationAliasSpecifier <- function(...) .alias_spec("OptimizationAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
IntegralAliasSpecifier <- function(...) .alias_spec("IntegralAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
AnalyticalAliasSpecifier <- function(...) .alias_spec("AnalyticalAliasSpecifier", ...)
#' @rdname alias_specifiers
#' @export
LinearAliasSpecifier <- function(...) .alias_spec("LinearAliasSpecifier", ...)

# ---------------------------------------------------------------------------
# Statistics, errors, ensembles

#' Statistics, error conditions and ensemble defaults
#'
#' `DEStats()` and `NLStats()` are the statistics objects of differential
#' equation and nonlinear solutions (`sol$stats` holds the same counters).
#' `FunctionArgumentsError(nargs)`, `TooFewArgumentsError(nargs)`,
#' `TooManyArgumentsError(nargs)` and `KeywordArgError(kw)` are the condition
#' objects SciMLBase signals for functions with the wrong signature or
#' unknown keywords; `numargs(f)` is the number of arguments of `f`.
#' `DEFAULT_PROB_FUNC`, `DEFAULT_OUTPUT_FUNC` and `DEFAULT_REDUCTION` are
#' the defaults of [EnsembleProblem()]; `generate_sim_seeds(n)` draws `n`
#' seeds and `default_rng_func(seed)` seeds the generator.
#' `WeightedEnsembleProblem(probs, weights)` combines several ensemble
#' problems with weights that sum to one (a mixture of ensembles); `solve()`
#' returns a `WeightedEnsembleSolution` with the solution of each ensemble in
#' `u` and the weights in `weights`.
#' `calculate_ensemble_errors(sim)` summarizes `calculate_solution_errors()`
#' over an ensemble.
#'
#' @param nf,njacs,nfactorizations,nsolve,nsteps,naccept,nreject,nw Counters.
#' @param nargs,kw Number of arguments or keyword name.
#' @param f A function.
#' @param prob,i,repeat_,sol,u,data,I Arguments of the ensemble functions.
#' @param n Number of seeds.
#' @param seed A seed.
#' @param probs A list of ensemble problems.
#' @param weights Their weights.
#' @param sim An ensemble solution.
#' @param ... Options.
#' @return A statistics object, a condition, a count, a value or a problem.
#' @examples
#' numargs(function(u, p, t) u)
#' DEStats(nf = 10L)
#' @name scimlbase_stats
NULL

#' @rdname scimlbase_stats
#' @export
DEStats <- function(nf = 0L, njacs = 0L, nsolve = 0L, nw = 0L, naccept = 0L, nreject = 0L)
  structure(list(nf = as.integer(nf), njacs = as.integer(njacs), nsolve = as.integer(nsolve), nw = as.integer(nw),
                 naccept = as.integer(naccept), nreject = as.integer(nreject)), class = "DEStats")

#' @rdname scimlbase_stats
#' @export
NLStats <- function(nf = 0L, njacs = 0L, nfactorizations = 0L, nsolve = 0L, nsteps = 0L)
  structure(list(nf = as.integer(nf), njacs = as.integer(njacs), nfactors = as.integer(nfactorizations),
                 nsolve = as.integer(nsolve), nsteps = as.integer(nsteps)), class = "NLStats")

.cond <- function(cls, msg) structure(class = c(cls, "error", "condition"), list(message = msg, call = NULL))

#' @rdname scimlbase_stats
#' @export
FunctionArgumentsError <- function(nargs) .cond("FunctionArgumentsError",
  sprintf("the function has %d arguments; problem functions take (u, p, t)", as.integer(nargs)))
#' @rdname scimlbase_stats
#' @export
TooFewArgumentsError <- function(nargs) .cond("TooFewArgumentsError",
  sprintf("the function has %d arguments, too few for (u, p, t)", as.integer(nargs)))
#' @rdname scimlbase_stats
#' @export
TooManyArgumentsError <- function(nargs) .cond("TooManyArgumentsError",
  sprintf("the function has %d arguments, more than (u, p, t)", as.integer(nargs)))
#' @rdname scimlbase_stats
#' @export
KeywordArgError <- function(kw) .cond("KeywordArgError", sprintf("unrecognized keyword argument `%s`", kw))

#' @rdname scimlbase_stats
#' @export
numargs <- function(f) {
  fm <- formals(args(f))
  if ("..." %in% names(fm)) return(Inf)
  length(fm)
}

#' @rdname scimlbase_stats
#' @export
DEFAULT_PROB_FUNC <- function(prob, i, repeat_) prob
#' @rdname scimlbase_stats
#' @export
DEFAULT_OUTPUT_FUNC <- function(sol, i) list(sol, FALSE)
#' @rdname scimlbase_stats
#' @export
DEFAULT_REDUCTION <- function(u, data, I) list(c(u, data), FALSE)

#' @rdname scimlbase_stats
#' @export
generate_sim_seeds <- function(n, seed = NULL) .with_seed(seed, sample.int(.Machine$integer.max, n))

#' @rdname scimlbase_stats
#' @export
default_rng_func <- function(seed) set.seed(seed)

#' @rdname scimlbase_stats
#' @export
WeightedEnsembleProblem <- function(probs, weights) {
  if (length(probs) != length(weights)) stop("give one weight per ensemble problem", call. = FALSE)
  if (abs(sum(weights) - 1) > 1e-8) stop("the weights must sum to one", call. = FALSE)
  structure(list(probs = probs, weights = as.numeric(weights)), class = c("WeightedEnsembleProblem", "SciMLProblem"))
}

#' @export
solve.WeightedEnsembleProblem <- function(a, b = NULL, ...) {
  sims <- lapply(a$probs, function(pr) solve(pr, b, ...))
  structure(list(u = sims, weights = a$weights), class = "WeightedEnsembleSolution")
}

#' @rdname scimlbase_stats
#' @export
WeightedEnsembleSolution <- function(u, weights) structure(list(u = u, weights = weights), class = "WeightedEnsembleSolution")

#' @rdname scimlbase_stats
#' @export
calculate_ensemble_errors <- function(sim, ...) {
  errs <- lapply(.subset2(sim, "u"), function(s) calculate_solution_errors(s, ...)$errors)
  nms <- unique(unlist(lapply(errs, names)))
  stats::setNames(lapply(nms, function(k) {
    v <- vapply(errs, function(e) e[[k]] %||% NA_real_, 0)
    list(mean = mean(v), median = stats::median(v), std = stats::sd(v))
  }), nms)
}

# ---------------------------------------------------------------------------
# Solution builders and solution utilities

#' Solution builders and utilities
#'
#' `build_solution(prob, alg, t, u, ...)` builds the solution of a problem:
#' an `ODESolution` for time dependent problems, a nonlinear solution
#' otherwise (`retcode`, `resid`, `stats`). `build_linear_solution(alg, u,
#' resid)` and `build_eigenvalue_solution(alg, values, vectors)` build linear
#' and eigenvalue solutions. `sensitivity_solution(sol, u, t)` replaces the
#' saved values. `strip_interpolation(sol)` drops the dense output,
#' `solution_new_retcode(sol, retcode)` changes the return code and
#' `wrap_sol(sol)` returns `sol`. `calculate_solution_errors(sol)` compares a
#' solution with the analytic solution of its problem (`analytic = function(u0,
#' p, t)` of [ODEFunction()]): the final error `final`, the
#' errors at the saved times `l2` and `l_inf`, and with `dense_errors` the
#' errors on 100 interpolated points; `has_analytic(f)` says whether a
#' function has one.
#'
#' @param prob,alg,t,u,retcode,resid,stats Fields of the solution.
#' @param values,vectors Eigenvalues and eigenvectors.
#' @param sol A solution.
#' @param timeseries_errors,dense_errors Which errors to compute.
#' @param f A problem function.
#' @param ... Further fields.
#' @return A solution, or the errors.
#' @examples
#' f <- ODEFunction(function(u, p, t) -u, analytic = function(u0, p, t) u0 * exp(-t))
#' sol <- solve(ODEProblem(f, 1, c(0, 1)), Tsit5())
#' calculate_solution_errors(sol)$errors$final
#' @name solution_builders
NULL

#' @rdname solution_builders
#' @export
build_solution <- function(prob, alg, t = NULL, u, retcode = "Success", resid = NULL, stats = list(), ...) {
  if (!is.null(prob$tspan) && !is.null(t)) return(ODESolution(t, u, prob, alg, retcode))
  .nl_new_solution("NonlinearSolution", prob, alg, u, retcode, resid = resid, stats = stats)
}

#' @rdname solution_builders
#' @export
build_linear_solution <- function(alg, u, resid = NULL, retcode = "Success", stats = list())
  .new_solution(c("LinearSolution", "NonlinearSolution"),
                list(t = 0, U = matrix(as.numeric(u), nrow = 1L), prob = NULL, alg = alg, retcode = retcode,
                     resid = resid, stats = stats, p_final = NULL, sys = NULL))

#' @rdname solution_builders
#' @export
build_eigenvalue_solution <- function(alg, values, vectors = NULL, retcode = "Success")
  structure(list(values = values, vectors = vectors, alg = alg, retcode = retcode), class = "EigenvalueSolution")

#' @rdname solution_builders
#' @export
sensitivity_solution <- function(sol, u, t = NULL) {
  env <- .sol_env(sol)
  U <- .interp_rows(u)
  env$U <- U
  if (!is.null(t)) env$t <- as.numeric(t)
  env$steps <- NULL
  env$dense <- FALSE
  sol
}

#' @rdname solution_builders
#' @export
strip_interpolation <- function(sol) {
  env <- .sol_env(sol)
  env$steps <- NULL
  env$dense <- FALSE
  sol
}

#' @rdname solution_builders
#' @export
solution_new_retcode <- function(sol, retcode) {
  env <- .sol_env(sol)
  env$retcode <- retcode
  sol
}

#' @rdname solution_builders
#' @export
wrap_sol <- function(sol) sol

#' @rdname solution_builders
#' @export
has_analytic <- function(f) is.function(f$analytic)

#' @rdname solution_builders
#' @export
calculate_solution_errors <- function(sol, timeseries_errors = TRUE, dense_errors = FALSE) {
  env <- .sol_env(sol)
  an <- env$prob$f$analytic
  if (!is.function(an)) stop("the problem function has no analytic solution", call. = FALSE)
  u0 <- env$prob$u0; p <- env$prob$p
  ex <- function(tt) as.numeric(an(u0, p, tt))
  n <- nrow(env$U)
  errs <- list(final = max(abs(env$U[n, ] - ex(env$t[n]))))
  if (timeseries_errors) {
    d <- vapply(seq_len(n), function(i) sqrt(mean((env$U[i, ] - ex(env$t[i]))^2)), 0)
    errs$l2 <- sqrt(mean(d^2))
    errs$l_inf <- max(d)
  }
  if (dense_errors && isTRUE(env$dense)) {
    ts <- seq(env$t[1L], env$t[n], length.out = 100L)
    d <- vapply(ts, function(tt) sqrt(mean((as.numeric(sol(tt)) - ex(tt))^2)), 0)
    errs$L2 <- sqrt(mean(d^2))
    errs$L_inf <- max(d)
  }
  env$errors <- errs
  list(sol = sol, errors = errs)
}

# ---------------------------------------------------------------------------
# Function wrappers and traits

#' Function wrappers and traits
#'
#' The wrappers fix all but one argument of a problem function, as
#' SciMLBase's differentiation helpers: `UJacobianWrapper(f, t, p)(u)` and
#' `UDerivativeWrapper(f, t, p)(u)` are `f(u, p, t)`, `ParamJacobianWrapper(f,
#' t, u)(p)` is `f(u, p, t)`, `TimeGradientWrapper(f, u, p)(t)` and
#' `TimeDerivativeWrapper(f, u, p)(t)` are `f(u, p, t)`, and
#' `JacobianWrapper(f, p)(u, t)` is `f(u, p, t)`. `unwrapped_f(f)` returns
#' the R function of an `ODEFunction` (or `f` itself).
#'
#' Traits of problem functions: `has_jac`, `has_tgrad`, `has_paramjac`,
#' `has_jvp`, `has_vjp`, `has_vjp_p`, `has_observed`, `has_sys`,
#' `has_symbolic_idxs`, `has_initialization_data`, `has_initializeprob`.
#' Traits of algorithms: `has_init`, `has_step`, `has_reinit`,
#' `has_stats`, `has_lazy_interpolation`, `allows_arbitrary_number_types`,
#' `allows_late_binding_tstops`, `allows_non_wiener_noise`, `allowscomplex`,
#' `requires_additive_noise`, `isdenseplot`, `has_global_error`; and of
#' optimizers `allowsconsjvp`, `allowsconsvjp`, `allowsfg`, `allowsfgh`,
#' `requiresconsjvp`, `requiresconsvjp`, `requireslagh`,
#' `supports_opt_cache_interface`. Of problems: `is_diagonal_noise`,
#' `isconcreteu0`, `problem_type`. Of integrators: `last_step_failed`,
#' `done`, `check_error` (the return code the integrator would stop with).
#'
#' @param f A problem function (or problem).
#' @param t,p,u Fixed arguments.
#' @param alg An algorithm.
#' @param prob A problem.
#' @param integ An integrator.
#' @param ... Unused.
#' @return A function, or a flag.
#' @examples
#' g <- UJacobianWrapper(function(u, p, t) -p * u, 0, 2)
#' g(3)
#' has_jac(ODEFunction(function(u, p, t) -u))
#' @param sol A solution.
#' @name function_wrappers
NULL

#' @rdname function_wrappers
#' @export
UJacobianWrapper <- function(f, t, p) { force(t); force(p); function(u) f(u, p, t) }
#' @rdname function_wrappers
#' @export
UDerivativeWrapper <- function(f, t, p) UJacobianWrapper(f, t, p)
#' @rdname function_wrappers
#' @export
ParamJacobianWrapper <- function(f, t, u) { force(t); force(u); function(p) f(u, p, t) }
#' @rdname function_wrappers
#' @export
TimeGradientWrapper <- function(f, u, p) { force(u); force(p); function(t) f(u, p, t) }
#' @rdname function_wrappers
#' @export
TimeDerivativeWrapper <- function(f, u, p) TimeGradientWrapper(f, u, p)
#' @rdname function_wrappers
#' @export
JacobianWrapper <- function(f, p) { force(p); function(u, t) f(u, p, t) }

#' @rdname function_wrappers
#' @export
unwrapped_f <- function(f) if (is.function(f)) f else f$f

.fn_of <- function(f) if (inherits(f, "SciMLProblem")) f$f else f
#' @rdname function_wrappers
#' @export
has_jac <- function(f) is.function(.fn_of(f)$jac)
#' @rdname function_wrappers
#' @export
has_tgrad <- function(f) is.function(.fn_of(f)$tgrad)
#' @rdname function_wrappers
#' @export
has_paramjac <- function(f) is.function(.fn_of(f)$paramjac)
#' @rdname function_wrappers
#' @export
has_jvp <- function(f) is.function(.fn_of(f)$jvp)
#' @rdname function_wrappers
#' @export
has_vjp <- function(f) is.function(.fn_of(f)$vjp)
#' @rdname function_wrappers
#' @export
has_vjp_p <- function(f) is.function(.fn_of(f)$vjp_p)
#' @rdname function_wrappers
#' @export
has_observed <- function(f) !is.null(.fn_of(f)$observed)
#' @rdname function_wrappers
#' @export
has_sys <- function(f) !is.null(f$sys)
#' @rdname function_wrappers
#' @export
has_symbolic_idxs <- function(f) inherits(f$sys, "System")
#' @rdname function_wrappers
#' @export
has_initialization_data <- function(f) !is.null(f$initialization_data)
#' @rdname function_wrappers
#' @export
has_initializeprob <- function(f) !is.null(f$initialization_data$initializeprob %||% f$initialization_data$prob)

#' @rdname function_wrappers
#' @export
has_init <- function(alg) inherits(alg, c("ODEAlgorithm", "SDEAlgorithm"))
#' @rdname function_wrappers
#' @export
has_step <- function(alg) has_init(alg)
#' @rdname function_wrappers
#' @export
has_reinit <- function(alg) has_init(alg)
#' @rdname function_wrappers
#' @export
has_stats <- function(alg) TRUE
#' @rdname function_wrappers
#' @export
has_lazy_interpolation <- function(alg) grepl("^Vern", alg$name %||% "")
#' @rdname function_wrappers
#' @export
allows_arbitrary_number_types <- function(alg) FALSE
#' @rdname function_wrappers
#' @export
allows_late_binding_tstops <- function(alg) TRUE
#' @rdname function_wrappers
#' @export
allows_non_wiener_noise <- function(alg) inherits(alg, "SDEAlgorithm") && (alg$order %||% 0.5) <= 0.5
#' @rdname function_wrappers
#' @export
allowscomplex <- function(alg) FALSE
#' @rdname function_wrappers
#' @export
requires_additive_noise <- function(alg) (alg$name %||% "") %in% c("SRA1", "SRA2", "SRA3", "SOSRA", "SOSRA2", "SKenCarp")
#' @rdname function_wrappers
#' @export
isdenseplot <- function(sol) isTRUE(.sol_env(sol)$dense)
#' @rdname function_wrappers
#' @export
has_global_error <- function(alg) FALSE
#' @rdname function_wrappers
#' @export
allowsconsjvp <- function(alg) FALSE
#' @rdname function_wrappers
#' @export
allowsconsvjp <- function(alg) FALSE
#' @rdname function_wrappers
#' @export
allowsfg <- function(alg) FALSE
#' @rdname function_wrappers
#' @export
allowsfgh <- function(alg) FALSE
#' @rdname function_wrappers
#' @export
requiresconsjvp <- function(alg) FALSE
#' @rdname function_wrappers
#' @export
requiresconsvjp <- function(alg) FALSE
#' @rdname function_wrappers
#' @export
requireslagh <- function(alg) FALSE
#' @rdname function_wrappers
#' @export
supports_opt_cache_interface <- function(alg) TRUE

#' @rdname function_wrappers
#' @export
is_diagonal_noise <- function(prob) (prob$noise_type %||% prob$f$noise_type %||% "diagonal") == "diagonal"
#' @rdname function_wrappers
#' @export
isconcreteu0 <- function(prob, ...) is.numeric(prob$u0) && all(is.finite(prob$u0))
#' @rdname function_wrappers
#' @export
problem_type <- function(prob) prob$problem_type %||% paste0("Standard", class(prob)[1L])

#' @rdname function_wrappers
#' @export
last_step_failed <- function(integ) {
  I <- .integ_I(integ)
  isFALSE(I$accept_prev)
}
#' @rdname function_wrappers
#' @export
done <- function(integ) isTRUE(.integ_I(integ)$done)
#' @rdname function_wrappers
#' @export
check_error <- function(integ) {
  I <- .integ_I(integ)
  if (isTRUE(integ$terminate)) return("Terminated")
  I$retcode %||% "Default"
}

# ---------------------------------------------------------------------------
# Initialization data, concrete problems, saving

#' Initialization data, concrete problems and saved values
#'
#' SciMLBase's interface between ModelingToolkit and the solvers.
#'
#' * `OverrideInitData(initializeprob, update_initializeprob, initializeprobmap,
#'   initializeprobpmap)` attaches an initialization problem to a problem
#'   built from R functions (`ODEProblem(..., initialization_data = )`):
#'   `get_initial_values(prob)` solves it and maps its solution to the
#'   initial values (`initializeprobmap(sol)`) and parameters
#'   (`initializeprobpmap(prob, sol)`); for problems of a `System` it returns
#'   the values ModelingToolkit's initialization found. `initialize_dae(integ)`
#'   (Julia's `initialize_dae!`) makes the algebraic values of an integrator
#'   consistent again. `is_overdetermined_initialization(prob)` tells whether
#'   the initialization system has more equations than unknowns.
#' * `remake_initialization_data(prob, u0, p)`, `updated_u0_p(prob, u0, p)`
#'   and `late_binding_update_u0_p(prob, u0, p)` give the initialization data
#'   or the numeric initial values and parameters of a problem with new
#'   values; `RemakeInitializationDataContext()` and
#'   `LateBindingUpdateU0PContext()` hold their arguments.
#'   `get_concrete_problem(prob)`, `get_concrete_u0(prob)` and
#'   `get_concrete_p(prob)` return the problem, initial values and parameters
#'   as the solver sees them; `promote_tspan()`, `promote_u0()`,
#'   `prepare_initial_state()` and `unwrap_parameters()` normalize them, and
#'   `get_updated_symbolic_problem(prob, u0, p)` is the problem with new
#'   values.
#' * `solve(prob, alg, save_idxs = )` saves only some components (indices or
#'   variables) in `sol$u`; symbolic indexing still sees all of them.
#'   `get_save_idxs_and_saved_subsystem(prob, save_idxs)` resolves them, and
#'   `SavedSubsystem()`, `SavedSubsystemWithFallback()`,
#'   `get_saved_state_idxs()`, `get_saved_subsystem()` and
#'   `get_saveable_values(integ, idxs)` describe what is saved.
#' * `save_discretes(integ)` (`save_discretes!`) records the current
#'   parameters (the discrete variables) in the solution's parameter history;
#'   `save_discretes_if_enabled()` and `save_final_discretes()` do so when
#'   they changed, and `create_parameter_timeseries_collection(sol)` returns
#'   the history (`t` and a matrix `p`).
#' * `set_ut(integ, u, t)` (`set_ut!`) sets state and time, `postamble(integ)`
#'   (`postamble!`) finishes the integration, `get_sol(integ)` is its
#'   solution so far, and `split_callbacks(cb)` separates continuous and
#'   discrete callbacks.
#' * `getindepsym(sys)` and `getindepsym_defaultt(sys)` give the independent
#'   variable, `get_root_indp(x)` the system behind a problem or solution,
#'   `plottable_indices(x)` the indices that are plotted and
#'   `symbolic_interpolation(sol, idxs, t)` interpolates expressions.
#' * `ImmutableODEProblem()` and `ImmutableNonlinearProblem()` are
#'   [ODEProblem()] and [NonlinearProblem()] (R problems are values);
#'   `StandardODEProblem()`, `StandardDAEProblem()`, `StandardDDEProblem()`,
#'   `StandardNonlinearProblem()`, `StandardBVProblem()` and
#'   `StandardSecondOrderBVProblem()` are the `problem_type` tags.
#'   `DefaultOptimizationCache()` is [OptimizationCache()].
#' * `ODENLStepData(prob)` is the implicit stage system of `ODEProblem(...,
#'   nlstep = TRUE)`; `SymbolicLinearInterface(update_A, update_b, sys)` and
#'   `get_new_A_b(prob, p)` rebuild the matrix and vector of a symbolic
#'   linear problem.
#' * `check_prob_alg_pairing(prob, alg)` signals an error when the algorithm
#'   does not solve that kind of problem, and `detect_cycles(map)` tells
#'   whether a map of variables to expressions refers to itself through a
#'   chain of entries.
#'
#' @param initializeprob,update_initializeprob,initializeprobmap,initializeprobpmap,metadata
#'   Fields of the initialization data.
#' @param prob A problem.
#' @param integ An integrator.
#' @param u0,p,t New values.
#' @param initializealg A DAE initialization algorithm.
#' @param tspan A time span.
#' @param save_idxs Indices or variables.
#' @param x A problem, solution or system.
#' @param idxs Indices.
#' @param sol A solution.
#' @param cb Callbacks.
#' @param sys A system.
#' @param update_A,update_b Functions of the parameters.
#' @param alg An algorithm.
#' @param map A list of formulas or equations `var ~ expr`.
#' @param indexp,saved_subsystem,fallback Fields of a saved subsystem.
#' @param ... Further fields or arguments.
#' @return Values as described above.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, c(1, 2), c(0, 1))
#' sol <- solve(prob, Tsit5(), save_idxs = 2)
#' sol$u[[1]]
#' detect_cycles(list(c("a", "b"), c("b", "a")))
#' @param u A state.
#' @name scimlbase_initialization
NULL

#' @rdname scimlbase_initialization
#' @export
OverrideInitData <- function(initializeprob, update_initializeprob = NULL, initializeprobmap = NULL,
                             initializeprobpmap = NULL, metadata = NULL)
  structure(list(initializeprob = initializeprob, update_initializeprob = update_initializeprob,
                 initializeprobmap = initializeprobmap, initializeprobpmap = initializeprobpmap,
                 metadata = metadata), class = "OverrideInitData")

#' @rdname scimlbase_initialization
#' @export
get_initial_values <- function(prob, ...) {
  d <- prob$initialization_data %||% prob$f$initialization_data
  if (inherits(d, "OverrideInitData")) {
    ip <- d$initializeprob
    if (is.function(d$update_initializeprob)) d$update_initializeprob(ip, prob)
    s <- solve(ip)
    ok <- successful_retcode(s)
    u0 <- if (is.function(d$initializeprobmap)) d$initializeprobmap(s) else prob$u0
    p <- if (is.function(d$initializeprobpmap)) d$initializeprobpmap(prob, s) else prob$p
    return(list(u0 = u0, p = p, success = ok))
  }
  list(u0 = prob$u0, p = prob$p, success = !identical(d$retcode, "InitialFailure"))
}

#' @rdname scimlbase_initialization
#' @export
initialize_dae <- function(integ, initializealg = NULL) {
  I <- .integ_I(integ)
  ini <- .dae_initialize(I$ctx, integ$t, integ$u, integ$du, initializealg, I$dtmax)
  if (isTRUE(ini$ok)) {
    integ$u[] <- ini$u
    if (!is.null(ini$du)) integ$du <- ini$du
    integ$u_modified <- TRUE
  }
  invisible(isTRUE(ini$ok))
}

#' @rdname scimlbase_initialization
#' @export
is_overdetermined_initialization <- function(prob) {
  ip <- (prob$initialization_data %||% list())$initializeprob
  inherits(ip, "NonlinearLeastSquaresProblem")
}

#' @rdname scimlbase_initialization
#' @export
remake_initialization_data <- function(prob, u0 = NULL, p = NULL) remake(prob, u0 = u0, p = p)$initialization_data

#' @rdname scimlbase_initialization
#' @export
updated_u0_p <- function(prob, u0 = NULL, p = NULL) {
  pr <- remake(prob, u0 = u0, p = p)
  list(u0 = pr$u0, p = pr$p)
}

#' @rdname scimlbase_initialization
#' @export
late_binding_update_u0_p <- function(prob, u0 = NULL, p = NULL) updated_u0_p(prob, u0, p)

#' @rdname scimlbase_initialization
#' @export
RemakeInitializationDataContext <- function(...) structure(list(...), class = "RemakeInitializationDataContext")

#' @rdname scimlbase_initialization
#' @export
LateBindingUpdateU0PContext <- function(...) structure(list(...), class = "LateBindingUpdateU0PContext")

#' @rdname scimlbase_initialization
#' @export
get_concrete_problem <- function(prob, ...) prob
#' @rdname scimlbase_initialization
#' @export
get_concrete_u0 <- function(prob, ...) as.numeric(prob$u0)
#' @rdname scimlbase_initialization
#' @export
get_concrete_p <- function(prob, ...) prob$p
#' @rdname scimlbase_initialization
#' @export
promote_tspan <- function(tspan) {
  tspan <- as.numeric(tspan)
  if (length(tspan) == 1L) c(0, tspan) else tspan[1:2]
}
#' @rdname scimlbase_initialization
#' @export
promote_u0 <- function(u0, p = NULL, t = NULL) as.numeric(u0)
#' @rdname scimlbase_initialization
#' @export
prepare_initial_state <- function(u0) as.numeric(u0)
#' @rdname scimlbase_initialization
#' @export
unwrap_parameters <- function(p) p
#' @rdname scimlbase_initialization
#' @export
get_updated_symbolic_problem <- function(prob, u0 = NULL, p = NULL, ...) remake(prob, u0 = u0, p = p)

#' @rdname scimlbase_initialization
#' @export
get_save_idxs_and_saved_subsystem <- function(prob, save_idxs) {
  if (is.null(save_idxs)) return(list(save_idxs = NULL, saved_subsystem = NULL))
  if (is.numeric(save_idxs) && !is_sym(save_idxs)) return(list(save_idxs = as.integer(save_idxs), saved_subsystem = NULL))
  vs <- if (is_sym(save_idxs)) list(save_idxs) else save_idxs
  us <- if (!is.null(prob$sys)) unknowns(prob$sys) else list()
  keys <- .keys(us)
  idx <- vapply(vs, function(v) {
    k <- match(.key(as_sym(v)), keys)
    if (is.na(k)) stop("save_idxs refers to ", format(v), ", which is not an unknown", call. = FALSE)
    k
  }, 0L)
  list(save_idxs = idx, saved_subsystem = SavedSubsystem(idx, NULL, NULL))
}

#' @rdname scimlbase_initialization
#' @export
SavedSubsystem <- function(indexp, t = NULL, saved_subsystem = NULL)
  structure(list(indexp = indexp, t = t, saved_subsystem = saved_subsystem), class = "SavedSubsystem")
#' @rdname scimlbase_initialization
#' @export
SavedSubsystemWithFallback <- function(saved_subsystem, fallback = NULL)
  structure(list(saved_subsystem = saved_subsystem, fallback = fallback), class = "SavedSubsystemWithFallback")
#' @rdname scimlbase_initialization
#' @export
get_saved_state_idxs <- function(x) if (inherits(x, "SavedSubsystem")) x$indexp else .sol_env(x)$save_idxs
#' @rdname scimlbase_initialization
#' @export
get_saved_subsystem <- function(sol) {
  idx <- .sol_env(sol)$save_idxs
  if (is.null(idx)) NULL else SavedSubsystem(idx)
}
#' @rdname scimlbase_initialization
#' @export
get_saveable_values <- function(integ, idxs = NULL) if (is.null(idxs)) integ$u else integ$u[idxs]

#' @rdname scimlbase_initialization
#' @export
save_discretes <- function(integ) {
  integ$param_t <- c(integ$param_t, integ$t)
  integ$param_p[[length(integ$param_p) + 1L]] <- integ$p
  invisible(integ)
}
#' @rdname scimlbase_initialization
#' @export
save_discretes_if_enabled <- function(integ) {
  last <- integ$param_p[[length(integ$param_p)]]
  if (!identical(last, integ$p)) save_discretes(integ)
  invisible(integ)
}
#' @rdname scimlbase_initialization
#' @export
save_final_discretes <- function(integ) save_discretes_if_enabled(integ)
#' @rdname scimlbase_initialization
#' @export
create_parameter_timeseries_collection <- function(sol) {
  env <- .sol_env(sol)
  env$param_history %||% list(t = env$t[1L], p = matrix(as.numeric(env$p_final %||% env$prob$p), nrow = 1L))
}

#' @rdname scimlbase_initialization
#' @export
set_ut <- function(integ, u, t) {
  set_u(integ, u)
  set_t(integ, t)
}
#' @rdname scimlbase_initialization
#' @export
postamble <- function(integ) solve_integrator(integ)
#' @rdname scimlbase_initialization
#' @export
get_sol <- function(integ) integ$sol
#' @rdname scimlbase_initialization
#' @export
split_callbacks <- function(cb) {
  cbs <- .as_callback_set(cb)
  list(continuous = cbs$continuous, discrete = cbs$discrete)
}

#' @rdname scimlbase_initialization
#' @export
getindepsym <- function(x) {
  sys <- get_root_indp(x)
  if (inherits(sys, "System")) get_iv(sys) else NULL
}
#' @rdname scimlbase_initialization
#' @export
getindepsym_defaultt <- function(x) getindepsym(x) %||% Symbol("t")
#' @rdname scimlbase_initialization
#' @export
get_root_indp <- function(x) {
  if (inherits(x, "System")) return(x)
  if (inherits(x, "SciMLSolution")) return(.sol_env(x)$sys %||% .sol_env(x)$prob$sys)
  x$sys
}
#' @rdname scimlbase_initialization
#' @export
plottable_indices <- function(x) {
  n <- if (inherits(x, "SciMLSolution")) ncol(as.matrix(x)) else length(x$u0)
  seq_len(n)
}
#' @rdname scimlbase_initialization
#' @export
symbolic_interpolation <- function(sol, idxs, t) sol(t, idxs = idxs)

#' @rdname scimlbase_initialization
#' @export
ImmutableODEProblem <- function(...) ODEProblem(...)
#' @rdname scimlbase_initialization
#' @export
ImmutableNonlinearProblem <- function(...) NonlinearProblem(...)

.ptype <- function(name) structure(list(name = name), class = c(name, "ProblemType"))
#' @rdname scimlbase_initialization
#' @export
StandardODEProblem <- function() .ptype("StandardODEProblem")
#' @rdname scimlbase_initialization
#' @export
StandardDAEProblem <- function() .ptype("StandardDAEProblem")
#' @rdname scimlbase_initialization
#' @export
StandardDDEProblem <- function() .ptype("StandardDDEProblem")
#' @rdname scimlbase_initialization
#' @export
StandardNonlinearProblem <- function() .ptype("StandardNonlinearProblem")
#' @rdname scimlbase_initialization
#' @export
StandardBVProblem <- function() .ptype("StandardBVProblem")
#' @rdname scimlbase_initialization
#' @export
StandardSecondOrderBVProblem <- function() .ptype("StandardSecondOrderBVProblem")

#' @rdname scimlbase_initialization
#' @export
DefaultOptimizationCache <- function(...) OptimizationCache(...)

#' @rdname scimlbase_initialization
#' @export
ODENLStepData <- function(prob) prob$f$nlstep_data

#' @rdname scimlbase_initialization
#' @export
SymbolicLinearInterface <- function(update_A, update_b, sys = NULL, ...)
  structure(list(update_A = update_A, update_b = update_b, sys = sys, ...), class = "SymbolicLinearInterface")

#' @rdname scimlbase_initialization
#' @export
get_new_A_b <- function(prob, p = prob$p) {
  f <- prob$f
  if (!is.function(f$update_A)) return(list(A = prob$A, b = prob$b))
  list(A = f$update_A(p), b = as.numeric(f$update_b(p)))
}

.prob_alg_kinds <- list(ODEProblem = c("ODEAlgorithm"), SDEProblem = c("SDEAlgorithm"),
                        DAEProblem = c("ODEAlgorithm"), NonlinearProblem = c("NonlinearAlgorithm"),
                        OptimizationProblem = c("OptimizationAlgorithm"))

#' @rdname scimlbase_initialization
#' @export
check_prob_alg_pairing <- function(prob, alg) {
  k <- .prob_alg_kinds[[class(prob)[1L]]]
  if (!is.null(k) && !is.null(alg) && !inherits(alg, k))
    stop(.cond("ProblemSolverPairingError",
               sprintf("%s() does not solve a %s", alg$name %||% class(alg)[1L], class(prob)[1L])))
  invisible(TRUE)
}

#' @rdname scimlbase_initialization
#' @export
detect_cycles <- function(map) {
  # each entry: a variable and the variables its value refers to
  ent <- lapply(map, function(e) {
    if (is.character(e)) return(list(lhs = e[1L], rhs = e[-1L]))
    if (inherits(e, "formula")) e <- as_equation(e)
    list(lhs = .key(as_sym(e$lhs)), rhs = .keys(.collect_vars(as_sym(e$rhs))))
  })
  lhs <- vapply(ent, `[[`, "", "lhs")
  adj <- lapply(ent, function(x) intersect(x$rhs, lhs))
  names(adj) <- lhs
  state <- stats::setNames(integer(length(lhs)), lhs)
  visit <- function(v) {
    if (state[[v]] == 1L) return(TRUE)
    if (state[[v]] == 2L) return(FALSE)
    state[[v]] <<- 1L
    for (w in adj[[v]]) if (visit(w)) return(TRUE)
    state[[v]] <<- 2L
    FALSE
  }
  any(vapply(lhs, visit, TRUE))
}
