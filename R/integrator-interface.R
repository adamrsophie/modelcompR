# The integrator interface of SciMLBase.jl (src/integrator_interface.jl and
# the init/solve!/step! of src/solve.jl), on the integrator of .integrate
# (R/solvers-ode.R): init() builds it, step_integrator() advances it,
# solve_integrator() finishes it, and the mutating functions (set_u!,
# set_t!, set_proposed_dt!, add_saveat!, savevalues!, reinit!, ...) change
# it between steps. The integrator object is the environment `integ` that
# callbacks receive; its element .I is the driver state.

#' The integrator interface
#'
#' `init(prob, alg, ...)` returns an integrator for an ODE, DAE or SDE
#' problem instead of solving it, as SciMLBase's `init`; the options are
#' those of [solve()]. `step_integrator(integ)` (Julia's `step!`) takes one
#' accepted step, or, with `dt`, steps until `dt` has passed (with
#' `stop_at_tdt = TRUE` it stops exactly there). `solve_integrator(integ)`
#' (Julia's `solve!`) integrates to the end of the time span and returns the
#' solution.
#'
#' The integrator is an environment with `u`, `t`, `p`, `uprev`, `tprev`, and
#' read only `dt` (the proposed step), `sol` (the solution so far), `opts`,
#' `alg`, `iter` and `retcode`; `integ$interp(t)` interpolates inside the
#' last step. It is the object callbacks receive, so the functions below also
#' work inside affects:
#' * `set_u(integ, u)`, `set_t(integ, t)` change the state; the next step
#'   restarts the method's caches. `derivative_discontinuity(integ)` marks a
#'   discontinuity without changing the state.
#' * `set_proposed_dt(integ, dt)` and `get_proposed_dt(integ)`;
#'   `set_abstol(integ, x)` and `set_reltol(integ, x)`.
#' * `add_saveat(integ, t)` adds save times; `savevalues(integ)` saves the
#'   current state now. [add_tstop()], [terminate()] and [u_modified()] are
#'   documented with the callbacks.
#' * `reinit(integ, u0, t0, tf)` starts the integration again (by default
#'   from the problem's initial values).
#' * `get_du(integ)` is the derivative of the state at the current time.
#'
#' @param prob An `ODEProblem`, `DAEProblem` or `SDEProblem`.
#' @param alg The algorithm, as for [solve()].
#' @param ... Solver options, as for [solve()].
#' @param integ An integrator returned by `init()` (or given to a callback).
#' @param dt A time span (for `step_integrator()`) or a step size.
#' @param stop_at_tdt Stop exactly after `dt`.
#' @param u,t A state and a time.
#' @param x A tolerance.
#' @param u0,t0,tf Initial values and time span for `reinit()`.
#' @param erase_sol Drop the saved values of the previous integration.
#' @param p New parameters (for an `OptimizationCache`).
#' @param flag Whether there is a discontinuity.
#' @return `init()` returns the integrator; `solve_integrator()` the
#'   solution; the others return the integrator invisibly (or the queried
#'   value).
#' @examples
#' prob <- ODEProblem(function(u, p, t) -p * u, 1, c(0, 1), 2)
#' integ <- init(prob, Tsit5(), reltol = 1e-8)
#' step_integrator(integ)
#' integ$t
#' step_integrator(integ, 0.25, stop_at_tdt = TRUE)
#' set_u(integ, 2 * integ$u)
#' sol <- solve_integrator(integ)
#' tail(sol$u, 1)
#' @name integrator_interface
NULL

#' @rdname integrator_interface
#' @export
init <- function(prob, alg = NULL, ...) UseMethod("init")

#' @export
init.default <- function(prob, alg = NULL, ...) {
  stop("init() is available for ODE, DAE and SDE problems, not for ", class(prob)[1L], call. = FALSE)
}

#' @export
init.ODEProblem <- function(prob, alg = NULL, ...) {
  alg <- alg %||% .default_ode_alg(prob)
  if (!inherits(alg, "ODEAlgorithm")) stop("not an ODE algorithm: ", class(alg)[1L], call. = FALSE)
  if (!.is_identity_mass(prob$f$mass_matrix) && !isTRUE(alg$mass_matrix))
    stop(alg$name, "() does not support mass matrices; use a mass matrix solver such as Rodas5P()", call. = FALSE)
  opts <- utils::modifyList(prob$kwargs %||% list(), list(...))
  .integ_public(.integ_init(prob, alg, opts))
}

#' @export
init.DAEProblem <- function(prob, alg = NULL, ...) {
  alg <- alg %||% DFBDF()
  opts <- utils::modifyList(prob$kwargs %||% list(), list(...))
  .integ_public(.integ_init(prob, alg, opts))
}

#' @export
init.SDEProblem <- function(prob, alg = NULL, ...) {
  opts <- utils::modifyList(prob$kwargs %||% list(), list(...))
  alg <- alg %||% .sde_default_alg(prob, opts$alg_hints)
  opts$alg_hints <- NULL
  # its own generator, seeded with `seed` or, as _resolve_rng of
  # StochasticDiffEq does without one, with a seed drawn from the global stream
  seed <- opts$seed %||% sample.int(.Machine$integer.max, 1L)
  opts$seed <- NULL
  setup <- .sde_setup(prob, alg, opts, rng = .rng_new(seed))
  I <- .integ_init(setup$prob, setup$alg, setup$opts)
  I$reinit_noise <- setup$reinit_noise
  I$post <- function(sol) .as_rode_solution(setup$post(sol))
  integ <- .integ_public(I)
  class(integ) <- c("SDEIntegrator", class(integ))
  integ
}

.integ_public <- function(I) {
  integ <- I$integ
  class(integ) <- c("ODEIntegrator", "DEIntegrator")
  integ
}

# Read only fields of the integrator (active bindings on integ).
.integ_bindings <- function(I) {
  integ <- I$integ
  bind <- function(name, f) {
    if (exists(name, envir = integ, inherits = FALSE)) rm(list = name, envir = integ)
    makeActiveBinding(name, f, integ)
  }
  # dt is the size of the last step (the initial step before the first);
  # get_proposed_dt() is the size of the next
  bind("dt", function() if (!is.null(integ$last_step)) abs(integ$last_step$dt) else I$dt)
  bind("EEst", function() I$EEst %||% 0)
  bind("sol", function() .integ_current_solution(I))
  bind("opts", function(v) {
    if (missing(v))
      return(list(abstol = I$abstol, reltol = I$reltol, dtmax = I$dtmax, dtmin = I$dtmin,
                  maxiters = I$maxiters, adaptive = I$adaptive, saveat = I$saveat,
                  save_everystep = I$save_everystep, verbose = TRUE))
    # integ$opts$abstol <- x and the like
    if (!is.null(v$abstol)) set_abstol(integ, v$abstol)
    if (!is.null(v$reltol)) set_reltol(integ, v$reltol)
    if (!is.null(v$dtmax)) I$dtmax <- v$dtmax
    if (!is.null(v$dtmin)) I$dtmin <- v$dtmin
  })
  bind("alg", function() I$alg)
  bind("iter", function() I$iter)
  bind("retcode", function() I$retcode)
  bind("tdir", function() I$tdir)
  integ$interp <- function(t) {
    st <- integ$last_step
    if (is.null(st)) return(integ$u)
    theta <- if (st$dt == 0) 1 else (t - st$t0) / st$dt
    st$interp(theta, st$dt, st$y0, st$y1, st$k)
  }
  invisible(integ)
}

.integ_I <- function(integ) {
  I <- if (is.environment(integ)) get0(".I", envir = integ, inherits = FALSE)
  if (is.null(I)) stop("not an integrator of init() or of a callback", call. = FALSE)
  I
}

.integ_current_solution <- function(I) {
  if (!is.null(I$final)) return(I$final)
  if (!is.null(I$early)) return(I$early)
  sol <- I$finish(I$retcode, I$ts, I$us, I$steps, I$alg_choice)
  if (!is.null(I$post)) sol <- I$post(sol)
  sol
}

#' @rdname integrator_interface
#' @export
step_integrator <- function(integ, dt = NULL, stop_at_tdt = FALSE) {
  I <- .integ_I(integ)
  if (is.null(dt)) {
    .integ_step(I)
    return(invisible(integ))
  }
  target <- integ$t + I$tdir * abs(dt)
  if (stop_at_tdt) add_tstop(integ, target)
  while (!I$done && (target - integ$t) * I$tdir > 0) .integ_step(I)
  invisible(integ)
}

#' @rdname integrator_interface
#' @export
solve_integrator <- function(integ) {
  if (inherits(integ, "ConvexOptimizationCache")) return(.cvx_solve(integ))
  if (inherits(integ, "OptimizationCache")) return(.opt_solve_cache(integ))
  if (inherits(integ, "LinearCache")) return(.lin_solve_cache(integ))
  I <- .integ_I(integ)
  if (!is.null(I$final)) return(I$final)
  while (!I$done) .integ_step(I)
  sol <- .integ_finalize(I)
  if (!is.null(I$post)) sol <- I$post(sol)
  I$final <- sol
  sol
}

#' @rdname integrator_interface
#' @export
set_u <- function(integ, u) {
  integ$u[] <- as.numeric(u)
  integ$u_modified <- TRUE
  invisible(integ)
}

#' @rdname integrator_interface
#' @export
set_t <- function(integ, t) {
  integ$t <- as.numeric(t)
  integ$u_modified <- TRUE
  invisible(integ)
}

#' @rdname integrator_interface
#' @export
derivative_discontinuity <- function(integ, flag = TRUE) {
  integ$u_modified <- isTRUE(flag)
  integ$u_modified_flag <- isTRUE(flag)
  invisible(integ)
}

#' @rdname integrator_interface
#' @export
set_proposed_dt <- function(integ, dt) {
  I <- .integ_I(integ)
  # inside a callback the driver takes the next step size from dtnew
  I$dt <- I$dtnew <- abs(as.numeric(dt))
  invisible(integ)
}

#' @rdname integrator_interface
#' @export
get_proposed_dt <- function(integ) {
  I <- .integ_I(integ)
  I$dtnew %||% I$dt
}

# has_tstop and first_tstop of SciMLBase
.integ_first_tstop <- function(integ) {
  I <- .integ_I(integ)
  if (is.null(I$tstops) || I$stop_idx > length(I$tstops)) NULL else I$tstops[I$stop_idx]
}

#' @rdname integrator_interface
#' @export
set_abstol <- function(integ, x) {
  I <- .integ_I(integ)
  I$abstol <- x
  I$ctx$abstol <- x
  invisible(integ)
}

#' @rdname integrator_interface
#' @export
set_reltol <- function(integ, x) {
  I <- .integ_I(integ)
  I$reltol <- x
  I$ctx$reltol <- x
  invisible(integ)
}

#' @rdname integrator_interface
#' @export
add_saveat <- function(integ, t) {
  I <- .integ_I(integ)
  sa <- sort(unique(c(I$saveat, as.numeric(t))), decreasing = I$tdir < 0)
  I$saveat <- sa
  # the next save time is the first one after the current time
  I$next_save <- which((sa - integ$t) * I$tdir > 0)[1L] %|% (length(sa) + 1L)
  if (is.null(I$opts$save_everystep)) I$save_everystep <- FALSE
  invisible(integ)
}

`%|%` <- function(a, b) if (is.na(a)) b else a

#' @rdname integrator_interface
#' @export
savevalues <- function(integ) {
  I <- .integ_I(integ)
  if (!(length(I$ts) && I$ts[length(I$ts)] == integ$t)) I$push_save(integ$t, integ$u)
  invisible(integ)
}

#' @rdname integrator_interface
#' @export
reinit <- function(integ, u0 = NULL, t0 = NULL, tf = NULL, erase_sol = TRUE, p = NULL) {
  # the conic cache rebuilds its cone program at the new p (R/optimization-convex.R)
  if (inherits(integ, "ConvexOptimizationCache")) return(.cvx_reinit(integ, u0 = u0, p = p))
  if (inherits(integ, "OptimizationCache")) {
    if (!is.null(u0)) integ$prob$u0[] <- as.numeric(u0)
    if (!is.null(p)) integ$prob$p[] <- as.numeric(p)
    integ$sol <- NULL
    return(invisible(integ))
  }
  I <- .integ_I(integ)
  prob <- I$prob
  if (!is.null(u0)) prob$u0[] <- as.numeric(u0)
  if (!is.null(t0) || !is.null(tf)) prob$tspan <- c(t0 %||% prob$tspan[1L], tf %||% prob$tspan[2L])
  old_ts <- I$ts; old_us <- I$us
  I$prob <- prob
  I$done <- FALSE
  I$early <- NULL
  I$final <- NULL
  # reinit_noise! of an SDE integrator: the noise restarts at the new t0
  if (is.function(I$reinit_noise)) I$reinit_noise(prob$tspan[1L])
  # .INTEG_INIT is the package's own initialization code (R/solvers-ode.R)
  eval(.INTEG_INIT, I)
  # the user's handle stays the same environment
  fresh <- I$integ
  for (nm in setdiff(ls(integ, all.names = TRUE), ls(fresh, all.names = TRUE))) rm(list = nm, envir = integ)
  for (nm in ls(fresh, all.names = TRUE)) assign(nm, get(nm, envir = fresh), envir = integ)
  I$integ <- integ
  integ$.I <- I
  .integ_bindings(I)
  if (!erase_sol) {
    I$ts <- c(old_ts, I$ts)
    I$us <- c(old_us, I$us)
  }
  invisible(integ)
}

#' @rdname integrator_interface
#' @export
get_du <- function(integ) {
  I <- .integ_I(integ)
  if (!is.null(integ$du)) return(integ$du)
  as.numeric(I$ctx$f(integ$u, integ$p, integ$t))
}

#' @export
print.DEIntegrator <- function(x, ...) {
  cat("Integrator at t = ", format(x$t), " (", x$alg$name, ")\n", sep = "")
  print(x$u)
  invisible(x)
}
