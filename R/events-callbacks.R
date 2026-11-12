# Symbolic events, after ModelingToolkitBase/src/systems/callbacks.jl
# (SymbolicAffect, AssignmentAffect, AffectSystem, SymbolicContinuousCallback,
# SymbolicDiscreteCallback, to_cb_vector, namespace_callback,
# continuous_events, discrete_events, compile_condition,
# compile_equational_affect, generate_continuous_callbacks,
# generate_discrete_callbacks, process_events) and
# ModelingToolkitBase/src/systems/imperative_affect.jl (ImperativeAffect,
# compile_functional_affect), with ModelingToolkitBase/src/utils.jl
# (modified_unknowns!). Events compile to the numerical callbacks of
# R/callbacks.R through the problem_callbacks hook.

# ---------------------------------------------------------------------------
# Affects

#' Affects of symbolic events
#'
#' An affect says what happens when an event fires. Symbolic affects are
#' equations relating the values after the event to the values before it,
#' written with [Pre()]: `list(v ~ -Pre(v))` reverses `v`, and
#' `list(x ~ Pre(x) + 1)` adds one to `x`. Variables without `Pre` are the
#' values after the event, so the affect equations are solved together with
#' the algebraic and observed equations of the system: unknowns that the
#' equations do not fix keep their values, and the others change to keep the
#' system consistent. Parameters change only when they are discrete
#' parameters: those given in `discrete_parameters` of the callback, and, as
#' an R convenience, variables declared with [Discrete()] or [discretes()]
#' (the `@discretes` of ModelingToolkit) that are the left-hand side of an
#' affect equation. ModelingToolkit requires discrete parameters to be time
#' dependent, `Parameter("c", t)`; modelcompR also accepts time independent
#' ones.
#'
#' `AssignmentAffect(list(x ~ expr))` assigns values computed entirely from
#' the values before the event: each right-hand side is wrapped in `Pre()`,
#' and parameters on the left-hand side become discrete parameters.
#'
#' `ImperativeAffect(f, modified, observed, ctx)` runs an R function. `modified`
#' and `observed` are named lists of variables or expressions (unnamed
#' variables are named after themselves). When the event fires, `f(modified,
#' observed, ctx, integ)` receives their current values as named lists and
#' returns a named list with new values for some of the `modified` entries,
#' which are written back to the unknowns and parameters. `integ` is the
#' integrator, so `terminate(integ)` stops the simulation.
#'
#' @param assignments A list of formulas `lhs ~ rhs` or pairs `list(lhs, rhs)`.
#' @param f An R function `function(modified, observed, ctx, integ)`.
#' @param modified,observed Named lists of variables (modified) or
#'   expressions (observed).
#' @param ctx Any R object passed to `f`.
#' @param skip_checks Skip the checks of the modified variables.
#' @return A `SymbolicAffect` or an `ImperativeAffect`.
#' @seealso [SymbolicContinuousCallback()], [SymbolicDiscreteCallback()].
#' @examples
#' t <- t_nounits
#' x <- Variable("x", t); k <- Parameter("k")
#' AssignmentAffect(list(x ~ 2 * x, k ~ k + 1))
#' ImperativeAffect(function(m, o, ctx, integ) list(x = -m$x), modified = list(x))
#' @export
AssignmentAffect <- function(assignments) {
  eqs <- list()
  dps <- list()
  for (a in .as_assignment_list(assignments)) {
    eqs[[length(eqs) + 1L]] <- Equation(a[[1L]], Pre(a[[2L]]))
    if (.is_param_lhs(a[[1L]])) dps[[length(dps) + 1L]] <- a[[1L]]
  }
  .symbolic_affect(eqs, dps)
}

.is_param_lhs <- function(l) is_sym(l) && .is_declared(l) && isparameter(l) && !isconstant(l)

.as_assignment_list <- function(x) {
  if (inherits(x, c("formula", "Equation"))) x <- list(x)
  out <- list()
  for (el in x) {
    if (inherits(el, c("formula", "Equation"))) {
      for (e in .as_equations(el)) out[[length(out) + 1L]] <- list(e$lhs, e$rhs)
    } else if (is.list(el) && length(el) == 2L && is_sym(el[[1L]])) {
      out[[length(out) + 1L]] <- list(el[[1L]], as_sym(el[[2L]]))
    } else {
      stop("an assignment must be a formula lhs ~ rhs or a pair list(lhs, rhs)", call. = FALSE)
    }
  }
  out
}

#' @rdname AssignmentAffect
#' @export
ImperativeAffect <- function(f, modified = list(), observed = list(), ctx = NULL, skip_checks = FALSE) {
  if (!is.function(f)) stop("the affect of an ImperativeAffect must be an R function", call. = FALSE)
  mod <- .named_syms(modified, "modified")
  obs <- .named_syms(observed, "observed")
  structure(list(f = f, obs = unname(obs), obs_syms = names(obs), modified = unname(mod),
                 mod_syms = names(mod), ctx = ctx, skip_checks = isTRUE(skip_checks)),
            class = c("ImperativeAffect", "mtk_event"))
}

# Named list of expressions; unnamed variables are named after themselves, as
# Julia's `(; x)`.
.named_syms <- function(x, what) {
  if (is.null(x) || !length(x)) return(list())
  if (is_sym(x)) x <- list(x)
  x <- as.list(x)
  nms <- names(x) %||% rep("", length(x))
  for (i in seq_along(x)) {
    x[[i]] <- as_sym(x[[i]])
    if (!nzchar(nms[i])) {
      if (!.is_declared(x[[i]]))
        stop("unnamed ", what, " entries must be variables; name the expression ", format(x[[i]]),
             call. = FALSE)
      nms[i] <- .name(x[[i]])
    }
  }
  names(x) <- nms
  x[!duplicated(nms)]
}

# SymbolicAffect(affect; discrete_parameters): equations, assignment pairs, an
# ImperativeAffect, or the tuple form list(f, modified, observed, ctx).
.symbolic_affect <- function(affect, discrete_parameters = list()) {
  if (is.null(affect)) return(NULL)
  if (inherits(affect, "ImperativeAffect")) return(affect)
  dps <- .flatten_vars(discrete_parameters)
  if (inherits(affect, "SymbolicAffect")) {
    affect$discrete_parameters <- .unique_syms(c(affect$discrete_parameters, dps))
    return(affect)
  }
  if (is.function(affect))
    stop("an affect must be equations, an AssignmentAffect or an ImperativeAffect; ",
         "wrap R functions in ImperativeAffect()", call. = FALSE)
  if (is.list(affect) && length(affect) && is.function(affect[[1L]]) &&
      !inherits(affect, c("formula", "Equation"))) {
    a <- unname(as.list(affect))
    return(ImperativeAffect(a[[1L]], modified = if (length(a) > 1L) a[[2L]] else list(),
                            observed = if (length(a) > 2L) a[[3L]] else list(),
                            ctx = if (length(a) > 3L) a[[4L]] else NULL))
  }
  items <- if (inherits(affect, c("Equation", "formula", "mtk_eqlist"))) list(affect) else as.list(affect)
  eqs <- list()
  for (it in items) {
    if (is.null(it)) next
    if (inherits(it, "SymbolicAffect")) {
      eqs <- c(eqs, it$affect)
      dps <- c(dps, it$discrete_parameters)
    } else if (inherits(it, c("Equation", "formula", "mtk_eqlist"))) {
      eqs <- c(eqs, .as_equations(it))
    } else if (is.list(it) && length(it) == 2L && is_sym(it[[1L]])) {
      # a Pair is an assignment, as in AssignmentAffect
      eqs[[length(eqs) + 1L]] <- Equation(it[[1L]], Pre(as_sym(it[[2L]])))
      if (.is_param_lhs(it[[1L]])) dps[[length(dps) + 1L]] <- it[[1L]]
    } else {
      stop("unrecognized affect format: affects must be equations or pairs list(lhs, rhs)", call. = FALSE)
    }
  }
  # variables declared as discretes on the lhs are discrete parameters
  for (e in eqs) {
    l <- e$lhs
    if (is_sym(l) && identical(.var_kind(l), "discrete")) dps[[length(dps) + 1L]] <- l
  }
  structure(list(affect = eqs, discrete_parameters = .unique_syms(dps)),
            class = c("SymbolicAffect", "mtk_event"))
}

#' @export
format.SymbolicAffect <- function(x, ...) vapply(x$affect, format, "")

#' @export
print.SymbolicAffect <- function(x, ...) {
  cat("Affect equations:\n")
  for (e in x$affect) cat("  ", format(e), "\n", sep = "")
  if (length(x$discrete_parameters))
    cat("Discrete parameters: ", paste(vapply(x$discrete_parameters, format, ""), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @export
print.ImperativeAffect <- function(x, ...) {
  fm <- function(v, n) paste0(vapply(v, format, ""), " => ", n, collapse = ", ")
  cat("ImperativeAffect(observed: [", if (length(x$obs)) fm(x$obs, x$obs_syms), "], modified: [",
      if (length(x$modified)) fm(x$modified, x$mod_syms), "])\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Callbacks

#' Symbolic continuous and discrete events
#'
#' `SymbolicContinuousCallback(conditions, affect)` is an event that fires
#' when one of the condition equations `lhs ~ rhs` becomes satisfied, that is
#' when `lhs - rhs` crosses zero; the time is located by root finding.
#' `affect` runs on upward crossings and `affect_neg` on downward ones (by
#' default the same affect; `affect_neg = NULL` makes downward crossings stop
#' the integrator without changing anything). A callback without an affect
#' only stops the integrator at the crossing, for discontinuities.
#'
#' `SymbolicDiscreteCallback(condition, affect)` is checked after each step.
#' The condition is a number `dt` (periodic event every `dt`), a numeric
#' vector or a list of numbers (events at these preset times; use `list(2)`
#' for a single time), a [Clock()] (periodic with a phase; it also fires at
#' the initial time when the clock ticks there), or a symbolic boolean such
#' as `Eq(t, t1) & (x < 50)` (built with symcompR's `Eq`, `<`, `&`), which
#' fires when it holds after a step; pass `tstops` to `solve()` so that the
#' integrator stops at the times the condition tests.
#'
#' Systems take these through `continuous_events` and `discrete_events`,
#' which also accept a formula `x ~ 0` (a condition only event), a pair
#' `list(conditions, affects)` or `conditions %=>% affects`, and lists of
#' these. A list of single formulas is a list of condition only events, as a
#' vector of equations in ModelingToolkit; write `list(x ~ 0) %=>% affects`
#' or `list(list(x ~ 0), affects)` for an event with an affect.
#'
#' @param conditions,condition Event conditions, see Details.
#' @param affect,affect_neg,initialize,finalize Affects: equations with
#'   [Pre()], an [AssignmentAffect()], an [ImperativeAffect()], or `NULL`.
#'   `initialize` runs when the integration starts and `finalize` when it
#'   ends.
#' @param rootfind `"left"` (default; the affect runs just before the
#'   crossing), `"right"` (just after) or `"none"` (at the end of the step).
#' @param reinitializealg Reinitialization after the affect. Symbolic affects
#'   are consistent by construction; kept for compatibility.
#' @param discrete_parameters Parameters the affect may change.
#' @param zero_crossing_id An identifier of the event.
#' @param initialize_save_discretes Kept for compatibility.
#' @param iv The independent variable; unused, kept for compatibility.
#' @param lhs,rhs Conditions and affects, for `%=>%`.
#' @return A callback object.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t, default = 1); v <- Variable("v", t, default = 0)
#' ball <- mtkcompile(System(list(D(x) ~ v, D(v) ~ -9.8), t, name = "ball",
#'   continuous_events = list(x ~ 0) %=>% list(v ~ -Pre(v))))
#' sol <- solve(ODEProblem(ball, list(), c(0, 2)))
#' min(sol[x]) >= 0
#' @export
SymbolicContinuousCallback <- function(conditions, affect = NULL, affect_neg = affect, initialize = NULL,
                                       finalize = NULL, rootfind = "left", reinitializealg = NULL,
                                       discrete_parameters = list(), zero_crossing_id = NULL,
                                       initialize_save_discretes = TRUE, iv = NULL) {
  if (inherits(conditions, "SymbolicContinuousCallback")) return(conditions)
  neg_same <- missing(affect_neg)
  if (inherits(conditions, "mtk_event_pair") || .looks_like_pair(conditions, "continuous")) {
    if (is.null(affect)) affect <- conditions[[2L]]
    conditions <- conditions[[1L]]
  }
  conds <- .as_equations(conditions)
  if (!all(vapply(conds, inherits, TRUE, "Equation")))
    stop("conditions of a continuous event must be equations", call. = FALSE)
  aff <- .symbolic_affect(affect, discrete_parameters)
  neg <- if (neg_same) aff else .symbolic_affect(affect_neg, discrete_parameters)
  structure(list(conditions = conds, affect = aff, affect_neg = neg, neg_same = neg_same,
                 initialize = .symbolic_affect(initialize, discrete_parameters),
                 finalize = .symbolic_affect(finalize, discrete_parameters),
                 rootfind = .rootfind_mode(rootfind),
                 reinitializealg = reinitializealg %||% .default_reinit(affect, affect_neg, initialize, finalize),
                 zero_crossing_id = zero_crossing_id %||% paste0("zc", .next_id()),
                 initialize_save_discretes = initialize_save_discretes),
            class = c("SymbolicContinuousCallback", "AbstractCallback", "mtk_event"))
}

.default_reinit <- function(...) {
  if (any(vapply(list(...), inherits, TRUE, "ImperativeAffect"))) "CheckInit" else "NoInit"
}

#' @rdname SymbolicContinuousCallback
#' @export
SymbolicDiscreteCallback <- function(condition, affect = NULL, initialize = NULL, finalize = NULL,
                                     reinitializealg = NULL, discrete_parameters = list(),
                                     initialize_save_discretes = TRUE, iv = NULL) {
  if (inherits(condition, "SymbolicDiscreteCallback")) return(condition)
  if (inherits(condition, "mtk_event_pair") || .looks_like_pair(condition, "discrete")) {
    if (is.null(affect)) affect <- condition[[2L]]
    condition <- condition[[1L]]
  }
  cond <- .discrete_condition(condition)
  structure(list(conditions = cond$value, timed = cond$timed,
                 affect = .symbolic_affect(affect, discrete_parameters),
                 initialize = .symbolic_affect(initialize, discrete_parameters),
                 finalize = .symbolic_affect(finalize, discrete_parameters),
                 reinitializealg = reinitializealg %||% .default_reinit(affect, initialize, finalize),
                 initialize_save_discretes = initialize_save_discretes),
            class = c("SymbolicDiscreteCallback", "AbstractCallback", "mtk_event"))
}

# is_timed_condition: a period, preset times, a clock, or a symbolic boolean.
.discrete_condition <- function(c) {
  if (inherits(c, "TimeDomain")) {
    if (!identical(c$type, "periodic"))
      stop("a clock given to SymbolicDiscreteCallback must be a periodic Clock()", call. = FALSE)
    return(list(value = c, timed = "clock"))
  }
  if (is.numeric(c) && !is_sym(c)) {
    if (length(c) == 1L) return(list(value = as.numeric(c), timed = "periodic"))
    return(list(value = sort(as.numeric(c)), timed = "preset"))
  }
  if (is.list(c) && !is_sym(c) && length(c) && all(vapply(c, function(z) is.numeric(z) && !is_sym(z), TRUE)))
    return(list(value = sort(as.numeric(unlist(c))), timed = "preset"))
  if (inherits(c, c("formula", "Equation")))
    stop("the condition of a discrete event is a boolean expression such as Eq(t, 1) or x > 1, ",
         "not an equation", call. = FALSE)
  if (is.list(c) && !is_sym(c))
    stop("vectors of symbolic conditions are not allowed for SymbolicDiscreteCallback; ",
         "combine them with & or |", call. = FALSE)
  list(value = as_sym(c), timed = NULL)
}

#' @rdname SymbolicContinuousCallback
#' @export
`%=>%` <- function(lhs, rhs) structure(list(lhs, rhs), class = "mtk_event_pair")

#' @export
print.SymbolicContinuousCallback <- function(x, ...) {
  cat("SymbolicContinuousCallback:\nConditions:\n")
  for (e in x$conditions) cat("  ", format(e), "\n", sep = "")
  .print_affects(x, c(affect = "Affect", affect_neg = "Negative-edge affect",
                      initialize = "Initialization affect", finalize = "Finalization affect"))
  invisible(x)
}

#' @export
print.SymbolicDiscreteCallback <- function(x, ...) {
  cat("SymbolicDiscreteCallback:\nCondition: ")
  c <- x$conditions
  cat(switch(x$timed %||% "",
             periodic = paste("every", c),
             preset = paste("at t =", paste(c, collapse = ", ")),
             clock = paste0("Clock(", c$dt, ", phase = ", c$phase, ")"),
             format(c)), "\n", sep = "")
  .print_affects(x, c(affect = "Affect", initialize = "Initialization affect",
                      finalize = "Finalization affect"))
  invisible(x)
}

.print_affects <- function(x, labels) {
  for (f in names(labels)) {
    a <- x[[f]]
    if (is.null(a) || (f == "affect_neg" && isTRUE(x$neg_same))) next
    cat(labels[[f]], ":\n", sep = "")
    if (inherits(a, "SymbolicAffect")) for (e in a$affect) cat("  ", format(e), "\n", sep = "")
    else {
      cat("  ")
      print(a)
    }
  }
}

# ---------------------------------------------------------------------------
# Normalizing the events given to System(), as to_cb_vector

.to_cb_vector <- function(ev, type) {
  if (is.null(ev)) return(list())
  if (.is_event_item(ev, type)) return(list(.make_callback(ev, type)))
  if (is.list(ev) && !is_sym(ev)) {
    out <- list()
    for (e in ev) out <- c(out, .to_cb_vector(e, type))
    return(out)
  }
  stop("cannot interpret ", class(ev)[1L], " as a ", type, " event", call. = FALSE)
}

.is_event_item <- function(ev, type) {
  if (inherits(ev, c("AbstractCallback", "mtk_event_pair"))) return(TRUE)
  if (type == "continuous" && inherits(ev, c("formula", "Equation"))) return(TRUE)
  if (type == "discrete" && (is_sym(ev) || inherits(ev, "TimeDomain") || (is.numeric(ev) && !is.list(ev))))
    return(TRUE)
  .looks_like_pair(ev, type)
}

.is_eq_like <- function(z) inherits(z, c("formula", "Equation", "mtk_eqlist"))

.is_condition_like <- function(c, type) {
  if (type == "continuous") {
    if (.is_eq_like(c)) return(TRUE)
    return(is.list(c) && !is_sym(c) && length(c) > 0L && all(vapply(c, .is_eq_like, TRUE)))
  }
  if (is_sym(c) || inherits(c, "TimeDomain") || (is.numeric(c) && !is.list(c))) return(TRUE)
  is.list(c) && length(c) > 0L && all(vapply(c, function(z) is.numeric(z) && !is_sym(z), TRUE))
}

.is_affect_like <- function(a) {
  if (is.null(a) || inherits(a, c("SymbolicAffect", "ImperativeAffect")) || .is_eq_like(a)) return(TRUE)
  if (!is.list(a) || is_sym(a) || !length(a)) return(FALSE)
  if (is.function(a[[1L]])) return(TRUE)
  all(vapply(a, function(z) .is_eq_like(z) || inherits(z, "SymbolicAffect") ||
               (is.list(z) && length(z) == 2L && is_sym(z[[1L]])), TRUE))
}

.looks_like_pair <- function(ev, type) {
  if (!is.list(ev) || is_sym(ev) || inherits(ev, c("mtk_eqlist", "AbstractCallback")) || length(ev) != 2L)
    return(FALSE)
  if (!.is_condition_like(ev[[1L]], type) || !.is_affect_like(ev[[2L]])) return(FALSE)
  # two equations are two condition only events, as a Vector{Equation}
  !(type == "continuous" && .is_eq_like(ev[[1L]]) && .is_eq_like(ev[[2L]]))
}

.make_callback <- function(ev, type) {
  if (type == "continuous") {
    if (inherits(ev, "SymbolicDiscreteCallback"))
      stop("a SymbolicDiscreteCallback was given as a continuous event", call. = FALSE)
    return(SymbolicContinuousCallback(ev))
  }
  if (inherits(ev, "SymbolicContinuousCallback"))
    stop("a SymbolicContinuousCallback was given as a discrete event", call. = FALSE)
  SymbolicDiscreteCallback(ev)
}

# ---------------------------------------------------------------------------
# Namespacing (namespace_callback, namespace_affects)

.namespace_event <- function(x, prefix, scopes = NULL) {
  ns <- function(e) .namespace(e, prefix, scopes)
  if (inherits(x, "SymbolicAffect")) {
    x$affect <- lapply(x$affect, ns)
    x$discrete_parameters <- lapply(x$discrete_parameters, ns)
    return(x)
  }
  if (inherits(x, "ImperativeAffect")) {
    x$obs <- lapply(x$obs, ns)
    x$modified <- lapply(x$modified, ns)
    return(x)
  }
  if (inherits(x, "SymbolicContinuousCallback")) {
    x$conditions <- lapply(x$conditions, ns)
  } else if (inherits(x, "SymbolicDiscreteCallback")) {
    if (is.null(x$timed)) x$conditions <- ns(x$conditions)
  } else {
    return(x)
  }
  for (f in c("affect", "affect_neg", "initialize", "finalize")) {
    if (!is.null(x[[f]])) x[f] <- list(ns(x[[f]]))
  }
  x
}

# Variables referenced by events, for the parameter discovery hook.
.event_vars <- function(events) {
  exprs <- list()
  add_affect <- function(a) {
    if (inherits(a, "SymbolicAffect")) {
      for (e in a$affect) exprs[[length(exprs) + 1L]] <<- e$lhs
      for (e in a$affect) exprs[[length(exprs) + 1L]] <<- e$rhs
      exprs <<- c(exprs, a$discrete_parameters)
    } else if (inherits(a, "ImperativeAffect")) {
      exprs <<- c(exprs, a$modified, a$obs)
    }
  }
  for (cb in events) {
    if (inherits(cb, "SymbolicContinuousCallback")) {
      for (e in cb$conditions) exprs <- c(exprs, list(e$lhs, e$rhs))
    } else if (inherits(cb, "SymbolicDiscreteCallback") && is.null(cb$timed)) {
      exprs <- c(exprs, list(cb$conditions))
    }
    for (f in c("affect", "affect_neg", "initialize", "finalize")) add_affect(cb[[f]])
  }
  Filter(function(v) !.is_derivative(v), .collect_vars(exprs))
}

# Parameters referenced by events become parameters of the system.
.events_parameters <- function(events, unknowns, ps, iv) {
  if (!length(events)) return(list())
  known <- .keys(c(unknowns, ps))
  # parameters in symbolic defaults of these are discovered too (collect_var!)
  vs <- .with_default_vars(.event_vars(events))
  .unique_syms(Filter(function(v) !.is_derivative(v) && isparameter(v) && !(.key(v) %in% known), vs))
}

# ---------------------------------------------------------------------------
# Accessors

#' Events of a system
#'
#' `continuous_events(sys)` and `discrete_events(sys)` return the events of a
#' system and its subsystems, namespaced; the `_toplevel` variants only the
#' events of the system itself (of its hierarchical parent for a compiled
#' system).
#'
#' @param sys A `System`.
#' @return A list of [SymbolicContinuousCallback()] or
#'   [SymbolicDiscreteCallback()] objects.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t)
#' sys <- System(list(D(x) ~ 1), t, name = "s", continuous_events = list(x ~ 1))
#' continuous_events(sys)
#' @export
continuous_events <- function(sys) .collect_hier(sys, function(s) .f(s, "continuous_events"))

#' @rdname continuous_events
#' @export
discrete_events <- function(sys) .collect_hier(sys, function(s) .f(s, "discrete_events"))

#' @rdname continuous_events
#' @export
continuous_events_toplevel <- function(sys) {
  par <- .f(sys, "parent")
  if (!is.null(par)) return(continuous_events_toplevel(par))
  .f(sys, "continuous_events")
}

#' @rdname continuous_events
#' @export
discrete_events_toplevel <- function(sys) {
  par <- .f(sys, "parent")
  if (!is.null(par)) return(discrete_events_toplevel(par))
  .f(sys, "discrete_events")
}

# modified_unknowns! for equations: the variables on the left-hand side,
# among `unknownlist` when given. The exported modified_unknowns() is owned by
# the jumps module (R/jumps-jumps.R), which handles jumps.
.modified_unknowns_eq <- function(munknowns, e, unknownlist = NULL) {
  e <- as_equation(e)
  vs <- Filter(function(v) !.is_derivative(v) && !identical(.var_kind(v), "iv"), .collect_vars(e$lhs))
  if (!is.null(unknownlist)) vs <- Filter(function(v) .contains_key(.flatten_vars(unknownlist), v), vs)
  .unique_syms(c(as.list(munknowns), vs))
}

# ---------------------------------------------------------------------------
# Compilation to numerical callbacks (process_events)

.event_ctx <- function(sys) {
  us <- unknowns(sys)
  ps <- .problem_params(sys)
  obs <- observed(sys)
  list(sys = sys, us = us, ps = ps, iv = .f(sys, "iv"), uk = .keys(us), pk = .keys(ps),
       obs_keys = vapply(obs, function(o) .key(o$lhs), ""), args = .codegen_args(sys),
       obs = .codegen_observed(sys))
}

.events_problem_callbacks <- function(sys) {
  ce <- continuous_events(sys)
  de <- discrete_events(sys)
  clocked <- .f(sys, "metadata")$mtk_clock_partitions
  if (!length(ce) && !length(de) && !length(clocked)) return(NULL)
  ctx <- .event_ctx(sys)
  cbs <- list()
  if (length(ce)) {
    # one vector callback per (rootfind, reinitializealg, initialize_save_discretes)
    modes <- vapply(ce, function(cb) .rootfind_mode(cb$rootfind), "")
    keys <- paste(modes, vapply(ce, function(cb) as.character(cb$reinitializealg)[1L], ""),
                  vapply(ce, function(cb) isTRUE(cb$initialize_save_discretes), TRUE))
    for (k in unique(keys)) {
      sel <- keys == k
      cbs[[length(cbs) + 1L]] <- .generate_vector_callback(ce[sel], ctx, modes[sel][1L])
    }
  }
  for (cb in de) cbs[[length(cbs) + 1L]] <- .generate_discrete_callback(cb, ctx)
  if (length(clocked)) cbs <- c(cbs, .clock_partition_callbacks(clocked, ctx))
  cbs <- Filter(Negate(is.null), cbs)
  if (!length(cbs)) return(NULL)
  do.call(CallbackSet, cbs)
}

.run_all <- function(fs) {
  fs <- Filter(Negate(is.null), fs)
  if (!length(fs)) return(NULL)
  function(integ) for (f in fs) f(integ)
}

# generate_callback for a group of continuous events with the same rootfind:
# one vector callback over all condition equations.
.generate_vector_callback <- function(cbs, ctx, mode) {
  conds <- list()
  eq2ev <- integer()
  for (i in seq_along(cbs)) {
    for (e in cbs[[i]]$conditions) {
      conds[[length(conds) + 1L]] <- e$lhs - e$rhs
      eq2ev <- c(eq2ev, i)
    }
  }
  if (!length(conds)) return(NULL)
  cfn <- .build_closure(conds, ctx$args, observed = ctx$obs, out = "vector")
  affs <- lapply(cbs, function(cb) .compile_affect(cb$affect, ctx))
  negs <- lapply(seq_along(cbs), function(i)
    if (isTRUE(cbs[[i]]$neg_same)) affs[[i]] else .compile_affect(cbs[[i]]$affect_neg, ctx))
  dispatch <- function(fs) {
    force(fs)
    function(integ, idx) {
      f <- fs[[eq2ev[idx]]]
      if (!is.null(f)) f(integ)
      invisible(NULL)
    }
  }
  VectorContinuousCallback(function(u, t, integ) cfn(u, integ$p, t), dispatch(affs), dispatch(negs),
                           initialize = .run_all(lapply(cbs, function(cb) .compile_affect(cb$initialize, ctx))),
                           finalize = .run_all(lapply(cbs, function(cb) .compile_affect(cb$finalize, ctx))),
                           rootfind = mode)
}

.generate_discrete_callback <- function(cb, ctx) {
  aff <- .compile_affect(cb$affect, ctx) %||% function(integ) invisible(NULL)
  init <- .compile_affect(cb$initialize, ctx)
  timed <- cb$timed %||% ""
  out <- if (timed == "preset") {
    PresetTimeCallback(cb$conditions, aff)
  } else if (timed == "periodic") {
    .periodic_event_callback(aff, cb$conditions)
  } else if (timed == "clock") {
    .periodic_event_callback(aff, cb$conditions$dt, cb$conditions$phase)
  } else {
    cfn <- .build_closure(cb$conditions, ctx$args, observed = ctx$obs, out = "scalar")
    DiscreteCallback(function(u, t, integ) isTRUE(as.logical(cfn(u, integ$p, t))), aff)
  }
  if (!is.null(init)) {
    own <- out$initialize
    out$initialize <- if (is.null(own)) init else function(integ) {
      own(integ)
      init(integ)
    }
  }
  out$finalize <- .compile_affect(cb$finalize, ctx)
  out
}

# A periodic event as DiffEqCallbacks' PeriodicCallback: ticks every dt from
# the initial time, or, for a clock with a phase, on the absolute lattice
# phase + k dt, where it also fires at the initial time when the lattice
# contains it. The final time never ticks.
.periodic_event_callback <- function(aff, dt, phase = NULL) {
  dt <- abs(as.numeric(dt))
  times_for <- function(tsp) {
    t0 <- tsp[1L]; tf <- tsp[2L]
    if (tf < t0) return(numeric())
    base <- if (is.null(phase)) t0 else phase
    k <- ceiling((t0 - base) / dt - 1e-9):floor((tf - base) / dt + 1e-9)
    tt <- base + k * dt
    eps <- 1e-12 * max(1, abs(tf))
    tt[tt > t0 + eps & tt < tf - eps]
  }
  st <- new.env(parent = emptyenv())
  cb <- DiscreteCallback(function(u, t, integ) .near_any(st$times, t), aff,
                         initialize = function(integ) {
                           st$times <- times_for(integ$prob$tspan)
                           if (!is.null(phase)) {
                             k <- (integ$t - phase) / dt
                             if (abs(k - round(k)) <= 1e-9) aff(integ)
                           }
                         })
  cb$tstops_fn <- times_for
  cb$periodic <- list(dt = dt, phase = phase %||% 0)
  cb
}

.compile_affect <- function(aff, ctx) {
  if (is.null(aff)) return(NULL)
  if (inherits(aff, "ImperativeAffect")) return(.compile_imperative_affect(aff, ctx))
  .compile_equational_affect(aff, ctx)
}

# Pre(...) nodes in a list of expressions, unique by key.
.find_ops <- function(exprs, op) {
  out <- list()
  seen <- character()
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    if (.is_op(z, op)) {
      k <- .key(z)
      if (!(k %in% seen)) {
        seen <<- c(seen, k)
        out[[length(out) + 1L]] <<- z
      }
      return(invisible())
    }
    for (a in .args(z)) walk(a)
  }
  for (e in exprs) walk(e)
  out
}

.unsolvable_callback <- function(eqs) {
  stop("The callback defined by the following equations:\n\n",
       paste(vapply(eqs, format, ""), collapse = "\n"),
       "\n\nis not solvable. Please check that the algebraic equations and affect equations are ",
       "correct, and that all parameters intended to be changed are passed in as ",
       "`discrete_parameters`.", call. = FALSE)
}

# The AffectSystem of MTK: the affect equations with the algebraic and
# observed equations of the system. Values after the event are found by
# explicit assignment where an equation defines an affected variable, and by
# Newton's method for the remaining equations.
.compile_equational_affect <- function(aff, ctx) {
  sys <- ctx$sys
  user <- aff$affect
  if (!length(user)) return(NULL)
  extra <- c(Filter(function(e) !.is_derivative(e$lhs), equations(sys)), observed(sys))
  all_eqs <- c(user, extra)
  is_user <- c(rep(TRUE, length(user)), rep(FALSE, length(extra)))
  # Pre(v) is the value of v before the event, an input of the affect
  pres <- .find_ops(unlist(lapply(all_eqs, function(e) list(e$lhs, e$rhs)), recursive = FALSE), "Pre")
  pre_vars <- lapply(pres, function(n) default_toterm(.args(n)[[1L]]))
  pre_ph <- lapply(seq_along(pres), function(i) Symbol(paste0(".pre", i)))
  prule <- .rule(pres, pre_ph)
  eqs <- lapply(all_eqs, function(e) Equation(.subst(e$lhs, prule), .subst(e$rhs, prule)))
  disc <- aff$discrete_parameters
  for (d in disc) {
    if (!(.key(d) %in% ctx$pk))
      stop("discrete parameter ", format(d), " of the affect is not a parameter of the system ",
           "with a value (it may be bound to an expression)", call. = FALSE)
  }
  # affected variables: unknowns and observables of the system, and discretes
  dv_keys <- c(ctx$uk, ctx$obs_keys)
  sides <- unlist(lapply(eqs, function(e) list(e$lhs, e$rhs)), recursive = FALSE)
  vars <- Filter(function(v) !.is_derivative(v), .collect_vars(sides))
  A <- .unique_syms(c(Filter(function(v) .key(v) %in% dv_keys, vars), disc))
  akeys <- .keys(A)
  classify <- function(allow) {
    assigned <- character()
    asg <- list()
    res <- list()
    for (i in seq_along(eqs)) {
      e <- eqs[[i]]
      l <- e$lhs
      lk <- if (.is_declared(l)) .key(l) else ""
      if (allow[i] && nzchar(lk) && lk %in% akeys && !(lk %in% assigned) &&
          !(lk %in% .keys(.collect_vars(e$rhs)))) {
        asg[[length(asg) + 1L]] <- e
        assigned <- c(assigned, lk)
      } else {
        res[[length(res) + 1L]] <- e$rhs - e$lhs
      }
    }
    list(asg = asg, res = res, assigned = assigned)
  }
  allow <- rep(TRUE, length(eqs))
  cl <- classify(allow)
  # an algebraic loop among assignments is solved with the residuals
  if (inherits(tryCatch(.toposort_observed(cl$asg), error = function(e) e), "error")) {
    cl <- classify(is_user)
    if (inherits(tryCatch(.toposort_observed(cl$asg), error = function(e) e), "error"))
      cl <- classify(rep(FALSE, length(eqs)))
  }
  asg <- .toposort_observed(cl$asg)
  # unknowns of the Newton solve: affected variables the residuals depend on
  need <- .keys(.collect_vars(cl$res))
  asg_lhs <- vapply(asg, function(e) .key(e$lhs), "")
  repeat {
    hit <- which(asg_lhs %in% need)
    more <- setdiff(.keys(.collect_vars(lapply(asg[hit], `[[`, "rhs"))), need)
    if (!length(more)) break
    need <- c(need, more)
  }
  Z <- Filter(function(v) .key(v) %in% need && !(.key(v) %in% cl$assigned), A)
  changed <- c(cl$assigned, .keys(Z))
  Wu <- Filter(function(v) .key(v) %in% changed && .key(v) %in% ctx$uk, A)
  Wp <- Filter(function(v) .key(v) %in% changed && .key(v) %in% ctx$pk, disc)
  if (!length(Wu) && !length(Wp) && !length(cl$res)) return(NULL)
  args <- list(list(name = "u", vars = ctx$us), list(name = "p", vars = ctx$ps),
               list(name = "pre", vars = pre_ph), list(name = "z", vars = Z),
               list(name = "t", vars = ctx$iv))
  obs_eqs <- c(asg, .param_binding_eqs(sys))
  pre_fn <- if (length(pre_vars)) .build_closure(pre_vars, ctx$args, observed = ctx$obs, out = "vector")
  z0_fn <- if (length(Z)) .build_closure(Z, ctx$args, observed = ctx$obs, out = "vector")
  res_fn <- if (length(cl$res)) .build_closure(cl$res, args, observed = obs_eqs, out = "vector")
  val_fn <- .build_closure(c(Wu, Wp), args, observed = obs_eqs, out = "vector")
  u_idx <- match(.keys(Wu), ctx$uk)
  p_idx <- match(.keys(Wp), ctx$pk)
  nu <- length(u_idx)
  function(integ) {
    u <- integ$u; p <- integ$p; t <- integ$t
    pre <- if (is.null(pre_fn)) numeric() else pre_fn(u, p, t)
    z <- numeric()
    if (!is.null(z0_fn)) {
      z <- .newton_solve(function(zz) res_fn(u, p, pre, zz, t), z0_fn(u, p, t))$x
    }
    vals <- val_fn(u, p, pre, z, t)
    if (!is.null(res_fn)) {
      # the residuals scale with the solved and written values only
      r <- res_fn(u, p, pre, z, t)
      if (any(!is.finite(r)) || max(abs(r)) > 1e-8 * max(1, abs(z), abs(vals))) .unsolvable_callback(user)
    }
    if (nu) integ$u[u_idx] <- vals[seq_len(nu)]
    if (length(p_idx)) integ$p[p_idx] <- vals[nu + seq_along(p_idx)]
    invisible(NULL)
  }
}

# compile_functional_affect: read the modified and observed values, call the
# user function, and write back what it returns.
.compile_imperative_affect <- function(aff, ctx) {
  mod <- aff$modified
  obs <- aff$obs
  for (m in mod) {
    k <- .key(m)
    if (!(k %in% c(ctx$uk, ctx$pk))) {
      what <- if (k %in% ctx$obs_keys) " is an observed variable" else " is not an unknown or parameter"
      stop("Expression ", format(m), " cannot be assigned to: it", what, " of the system; ",
           "only unknowns and parameters may be updated by an affect.", call. = FALSE)
    }
  }
  if (!aff$skip_checks) {
    allk <- c(ctx$uk, ctx$pk, ctx$obs_keys, .keys(.parameters_of(ctx$sys, TRUE)), .key(ctx$iv))
    for (o in obs) {
      bad <- setdiff(.keys(Filter(function(v) !.is_derivative(v), .collect_vars(o))), allk)
      if (length(bad))
        stop("Observed expression ", format(o), " in affect refers to missing variable(s)", call. = FALSE)
    }
  }
  mod_fn <- if (length(mod)) .build_closure(mod, ctx$args, observed = ctx$obs, out = "vector")
  obs_fn <- if (length(obs)) .build_closure(obs, ctx$args, observed = ctx$obs, out = "vector")
  mk <- .keys(mod)
  is_u <- mk %in% ctx$uk
  idx <- ifelse(is_u, match(mk, ctx$uk), match(mk, ctx$pk))
  mnames <- aff$mod_syms
  onames <- aff$obs_syms
  f <- aff$f
  fctx <- aff$ctx
  nf <- if ("..." %in% names(formals(f))) 4L else min(4L, length(formals(f)))
  function(integ) {
    u <- integ$u; p <- integ$p; t <- integ$t
    m <- if (length(mod)) as.list(stats::setNames(mod_fn(u, p, t), mnames)) else list()
    o <- if (length(obs)) as.list(stats::setNames(obs_fn(u, p, t), onames)) else list()
    res <- do.call(f, list(m, o, fctx, integ)[seq_len(nf)])
    if (is.null(res) || !length(res)) return(invisible(NULL))
    res <- as.list(res)
    bad <- setdiff(names(res), mnames)
    if (length(bad) || is.null(names(res)))
      stop("Tried to write back to ", paste(bad, collapse = ", "),
           " from affect; only declared states (", paste(mnames, collapse = ", "),
           ") may be written to.", call. = FALSE)
    for (nm in names(res)) {
      i <- match(nm, mnames)
      v <- as.numeric(res[[nm]])
      if (is_u[i]) integ$u[idx[i]] <- v else integ$p[idx[i]] <- v
    }
    invisible(NULL)
  }
}

# ---------------------------------------------------------------------------
# Hooks

.mtk$hooks$normalize_continuous_events <- function(events) .to_cb_vector(events, "continuous")
.mtk$hooks$normalize_discrete_events <- function(events) .to_cb_vector(events, "discrete")
.mtk$hooks$event_parameters <- function(events, unknowns, ps, iv) .events_parameters(events, unknowns, ps, iv)
.mtk$hooks$namespace_event <- function(x, prefix, scopes) .namespace_event(x, prefix, scopes)
.mtk$hooks$problem_callbacks <- function(sys) .events_problem_callbacks(sys)
