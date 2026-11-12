# Algorithm traits and integrator utilities of OrdinaryDiffEqCore and
# StochasticDiffEq (alg_utils.jl, integrator_utils.jl, integrator
# interface): orders, adaptivity, compatibility, caches, stop times, the
# random number generator, resizing the state and stepping internals.

#' Algorithm traits and integrator utilities
#'
#' Traits of ODE and SDE algorithms: `alg_order(alg)` (the strong order for
#' SDE methods), `get_current_alg_order(alg)`, `isadaptive(alg)`,
#' `alg_stability_size(alg)`, `alg_compatible(prob, alg)` (whether the method
#' takes the problem's noise), `alg_mass_matrix_compatible(alg)`,
#' `alg_needs_extra_process(alg)` (methods that use a second Brownian
#' increment), `alg_can_repeat_jac(alg)`, `alg_control_rate(alg)`,
#' `delta_default(alg)` (the weight of the drift error estimate),
#' `is_split_step(alg)`, `isdiscrete(alg)`, `supports_regular_jumps(alg)`,
#' `supports_solve_rng(prob, alg)`, `adaptive_alg(alg)` and `unwrap_alg(x)`
#' (the method of an integrator, or the current one of a composite method).
#' Of functions and problems: `isinplace(f)` (always `FALSE`: R functions
#' return their values), `issplit(prob)`.
#'
#' Integrator utilities: `get_dt(integ)`, `get_tmp_cache(integ)`,
#' `u_cache()`, `du_cache()`, `full_cache()`, `user_cache()`, `rand_cache()`
#' and `ratenoise_cache()` (the vectors of the integrator's state that a
#' callback may use as scratch space), `has_tstop(integ)`,
#' `first_tstop(integ)`, `pop_tstop(integ)` (`pop_tstop!`),
#' `get_rng(integ)`, `set_rng(integ, seed)` (`set_rng!`), `has_rng(integ)`,
#' `auto_dt_reset(integ)` (a fresh initial step size estimate),
#' `change_t_via_interpolation(integ, t)` (move back inside the last step),
#' `reeval_internals_due_to_modification(integ)`, `check_keywords(integ,
#' kwargs, warnlist)`, `addat(integ, idxs, value)` and `deleteat(integ,
#' idxs)` (grow or shrink the state; the right-hand side must follow the new
#' size), and the cache resizing hooks `resize_non_user_cache()`,
#' `addat_non_user_cache()`, `deleteat_non_user_cache()`.
#'
#' @param alg An algorithm.
#' @param prob A problem.
#' @param f A function.
#' @param x An integrator or algorithm.
#' @param integ An integrator.
#' @param seed A seed.
#' @param t A time inside the last step.
#' @param kwargs,warnlist Keyword names and the accepted ones.
#' @param idxs,value Indices and the value of new components.
#' @param ... Unused.
#' @return A trait value, a vector, or the integrator invisibly.
#' @examples
#' alg_order(SRIW1())
#' isadaptive(EM())
#' alg_compatible(SDEProblem(function(u, p, t) u, function(u, p, t) u, 1, c(0, 1)), SOSRA())
#' @name solver_traits
NULL

#' @rdname solver_traits
#' @export
alg_order <- function(alg) {
  if (inherits(alg, "CompositeAlgorithm")) return(max(vapply(alg$algs, alg_order, 0)))
  as.numeric(alg$order %||% NA_real_)
}
#' @rdname solver_traits
#' @export
get_current_alg_order <- function(alg, ...) alg_order(alg)
#' @rdname solver_traits
#' @export
isadaptive <- function(alg) {
  if (inherits(alg, "CompositeAlgorithm")) return(all(vapply(alg$algs, isadaptive, TRUE)))
  isTRUE(alg$adaptive)
}
#' @rdname solver_traits
#' @export
alg_stability_size <- function(alg) alg$stability_size %||% 1
#' @rdname solver_traits
#' @export
alg_compatible <- function(prob, alg) {
  if (!inherits(alg, "SDEAlgorithm")) return(TRUE)
  type <- prob$noise_type %||% prob$f$noise_type %||% "diagonal"
  type %in% alg$noise
}
#' @rdname solver_traits
#' @export
alg_mass_matrix_compatible <- function(alg) isTRUE(alg$mass_matrix_ok %||% alg$mass_matrix)
#' @rdname solver_traits
#' @export
alg_needs_extra_process <- function(alg) isTRUE(alg$needs_dZ)
#' @rdname solver_traits
#' @export
alg_can_repeat_jac <- function(alg) TRUE
#' @rdname solver_traits
#' @export
alg_control_rate <- function(alg) FALSE
#' @rdname solver_traits
#' @export
delta_default <- function(alg) alg$delta %||% 1
#' @rdname solver_traits
#' @export
is_split_step <- function(alg) isTRUE(alg$split)
#' @rdname solver_traits
#' @export
isdiscrete <- function(alg) inherits(alg, c("DiscreteAlgorithm", "FunctionMap", "SSAStepper"))
#' @rdname solver_traits
#' @export
supports_regular_jumps <- function(alg) inherits(alg, "TauLeapingAlgorithm")
#' @rdname solver_traits
#' @export
supports_solve_rng <- function(prob, alg) TRUE
#' @rdname solver_traits
#' @export
adaptive_alg <- function(alg) {
  # DiffEqNoiseProcess: the rejection sampling scheme of an RSWM or a process
  if (inherits(alg, "AbstractNoiseProcess")) alg <- alg$rswm
  if (inherits(alg, "RSWM")) return(alg$adaptivealg)
  alg
}
#' @rdname solver_traits
#' @export
unwrap_alg <- function(x, ...) {
  if (inherits(x, "DEIntegrator")) {
    I <- .integ_I(x)
    if (!is.null(I$algs) && length(I$algs) > 1L) return(I$algs[[I$cur]])
    return(I$alg)
  }
  x
}
#' @rdname solver_traits
#' @export
isinplace <- function(f, ...) FALSE
#' @rdname solver_traits
#' @export
issplit <- function(prob) inherits(prob, c("SplitODEProblem", "SplitSDEProblem"))

#' @rdname solver_traits
#' @export
get_dt <- function(integ) integ$dt
#' @rdname solver_traits
#' @export
get_tmp_cache <- function(integ) list(0 * integ$u, 0 * integ$u)
#' @rdname solver_traits
#' @export
u_cache <- function(integ) list(integ$u, integ$uprev)
#' @rdname solver_traits
#' @export
du_cache <- function(integ) list(get_du(integ))
#' @rdname solver_traits
#' @export
full_cache <- function(integ) c(u_cache(integ), du_cache(integ))
#' @rdname solver_traits
#' @export
user_cache <- function(integ) get_tmp_cache(integ)
#' @rdname solver_traits
#' @export
rand_cache <- function(integ) list()
#' @rdname solver_traits
#' @export
ratenoise_cache <- function(integ) list()

#' @rdname solver_traits
#' @export
has_tstop <- function(integ) !is.null(.integ_first_tstop(integ))
#' @rdname solver_traits
#' @export
first_tstop <- function(integ) .integ_first_tstop(integ)
#' @rdname solver_traits
#' @export
pop_tstop <- function(integ) {
  I <- .integ_I(integ)
  t <- .integ_first_tstop(integ)
  if (!is.null(t)) I$stop_idx <- I$stop_idx + 1L
  t
}
#' @rdname solver_traits
#' @export
get_rng <- function(integ) get0(".Random.seed", envir = globalenv(), inherits = FALSE)
#' @rdname solver_traits
#' @export
set_rng <- function(integ, seed) {
  set.seed(seed)
  invisible(integ)
}
#' @rdname solver_traits
#' @export
has_rng <- function(integ) TRUE

#' @rdname solver_traits
#' @export
auto_dt_reset <- function(integ) {
  I <- .integ_I(integ)
  a <- I$algs[[I$cur]]
  I$dt <- I$start_dt(a, I$caches[[I$cur]])
  I$stats$nf <- I$stats$nf + 2L
  invisible(integ)
}

#' @rdname solver_traits
#' @export
change_t_via_interpolation <- function(integ, t) {
  st <- integ$last_step
  if (is.null(st)) stop("there is no step to interpolate in", call. = FALSE)
  integ$u[] <- integ$interp(t)
  integ$t <- t
  integ$u_modified <- TRUE
  invisible(integ)
}

#' @rdname solver_traits
#' @export
reeval_internals_due_to_modification <- function(integ, ...) u_modified(integ, TRUE)

#' @rdname solver_traits
#' @export
check_keywords <- function(integ, kwargs, warnlist) {
  bad <- setdiff(kwargs, warnlist)
  if (length(bad)) warning("unrecognized keyword arguments: ", paste(bad, collapse = ", "), call. = FALSE)
  invisible(length(bad) == 0L)
}

#' @rdname solver_traits
#' @export
addat <- function(integ, idxs, value = 0) {
  u <- as.numeric(integ$u)
  n <- length(u) + length(idxs)
  out <- numeric(n)
  keep <- setdiff(seq_len(n), idxs)
  out[keep] <- u
  out[idxs] <- value
  .integ_resize(integ, out)
}

#' @rdname solver_traits
#' @export
deleteat <- function(integ, idxs) .integ_resize(integ, as.numeric(integ$u)[-idxs])

.integ_resize <- function(integ, unew) {
  I <- .integ_I(integ)
  integ$u <- unew
  integ$uprev <- unew
  I$ctx$n <- length(unew)
  if (!is.null(I$ctx$mass)) I$ctx$mass <- diag(length(unew))
  integ$last_step <- NULL
  integ$u_modified <- TRUE
  invisible(integ)
}

#' @rdname solver_traits
#' @export
resize_non_user_cache <- function(integ, ...) invisible(integ)
#' @rdname solver_traits
#' @export
addat_non_user_cache <- function(integ, ...) invisible(integ)
#' @rdname solver_traits
#' @export
deleteat_non_user_cache <- function(integ, ...) invisible(integ)
