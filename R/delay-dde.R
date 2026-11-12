# Delay differential equations, after ModelingToolkitBase/src/problems/
# ddeproblem.jl and sddeproblem.jl (DDEFunction, DDEProblem, SDDEFunction,
# SDDEProblem), systems/codegen.jl (generate_rhs with wrap_delays,
# generate_history) and DelayDiffEq.jl (MethodOfSteps).
#
# A delayed variable x(t - tau) is written EvalAt(t - tau)(x). In generated
# code every evaluated term x(s) becomes a placeholder filled from the
# history function h(p, s), so the right-hand side is f(u, h, p, t) as in
# Julia. MethodOfSteps(alg) integrates with any ODE algorithm of the core
# driver: a passive callback records each accepted step, the history
# interpolates the recorded steps (and the user history before t0), steps
# whose stages reach into the step itself are repeated by fixed point
# iteration on the step's own interpolant. Discontinuities follow
# DelayDiffEq (solve.jl initialize_tstops_d_discontinuities_propagated,
# integrators/utils.jl handle_discontinuities!, add_next_discontinuities!,
# track.jl track_propagated_discontinuities!; DelayDiffEq in OrdinaryDiffEq.jl
# commit 4460ff8): each has an order, the one at t0 order_discontinuity_t0;
# reaching one propagates it by the constant lags with order + 1 (the same
# order for neutral problems) while that is at most the method's order + 1,
# and after each rejected step the dependent lags are searched for crossings
# of earlier discontinuities, which become new stopping times.

# ---------------------------------------------------------------------------
# Evaluated terms x(s) in generated code

# Evaluated or delayed terms x(s) (s not the independent variable) of
# unknowns or observed variables in `exprs`.
.dde_terms <- function(exprs, iv) {
  out <- list()
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    if (.type(z) == "Function" && identical(.var_kind(z), "unknown") && .is_delay_or_evalat(z, iv)) {
      out[[length(out) + 1L]] <<- z
    }
    for (a in .args(z)) walk(a)
  }
  for (e in exprs) walk(e)
  .unique_syms(out)
}

# Replace evaluated observed variables w(s) by their right-hand side at s.
.dde_expand_observed <- function(exprs, obs, iv) {
  if (!length(obs)) return(exprs)
  onames <- vapply(obs, function(o) .name(o$lhs) %||% "", "")
  for (it in 1:50) {
    terms <- .dde_terms(exprs, iv)
    hit <- Filter(function(d) .name(d) %in% onames, terms)
    if (!length(hit)) break
    to <- lapply(hit, function(d) {
      o <- obs[[match(.name(d), onames)]]
      .subst(o$rhs, .rule(list(iv), list(.args(d)[[1L]])))
    })
    rule <- .rule(hit, to)
    exprs <- lapply(exprs, function(e) .subst(e, rule))
  }
  exprs
}

# Code for expressions holding evaluated terms. `args` are the usual
# argument groups; the result has fn(<args>, vals) with vals the values of the
# terms, tfun(<args>) their times, and idx their state indexes.
.evalat_codegen <- function(exprs, sys, args, observed, out = "vector") {
  iv <- .f(sys, "iv")
  us <- unknowns(sys)
  obs_all <- observed(sys)
  exprs <- .dde_expand_observed(exprs, obs_all, iv)
  observed <- lapply(observed, function(o) Equation(o$lhs, .dde_expand_observed(list(o$rhs), obs_all, iv)[[1L]]))
  terms <- .dde_terms(c(exprs, lapply(observed, `[[`, "rhs")), iv)
  unames <- vapply(us, function(u) .name(u) %||% "", "")
  idx <- match(vapply(terms, .name, ""), unames)
  if (anyNA(idx))
    .stopf("the evaluated variable %s is not an unknown of the system", format(terms[[which(is.na(idx))[1L]]]))
  fn <- .build_closure(exprs, c(args, list(list(name = "vals", vars = terms))), observed = observed, out = out)
  times <- lapply(terms, function(d) .args(d)[[1L]])
  tfun <- if (length(terms)) .build_closure(times, args, observed = observed, out = "vector") else NULL
  list(fn = fn, tfun = tfun, idx = idx, n = length(terms))
}

# Values of the evaluated terms from a state function S(s) (a vector).
.evalat_values <- function(ts, idx, S) {
  vals <- numeric(length(ts))
  ut <- unique(ts)
  for (s in ut) {
    st <- S(s)
    k <- which(ts == s)
    vals[k] <- st[idx[k]]
  }
  vals
}

# f(u, h, p, t) for expressions of a delay system.
.dde_build <- function(exprs, sys, out = "vector") {
  cg <- .evalat_codegen(exprs, sys, .codegen_args(sys), .codegen_observed(sys), out = out)
  fn <- cg$fn; tfun <- cg$tfun; idx <- cg$idx
  if (!cg$n) return(function(u, h, p, t) fn(u, p, t, numeric()))
  function(u, h, p, t) {
    ts <- tfun(u, p, t)
    fn(u, p, t, .evalat_values(ts, idx, function(s) h(p, s)))
  }
}

#' Delay systems
#'
#' `is_dde(sys)` tells whether the equations of a system contain delayed
#' unknowns such as `EvalAt(t - tau)(x)`, written `x(t - τ)` in
#' ModelingToolkit.
#'
#' @param sys A `System`.
#' @return A logical.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)")
#' is_dde(System(list(D(x) ~ -EvalAt(t - 1)(x)), t, name = "dde"))
#' @export
is_dde <- function(sys) isTRUE(.f(sys, "is_dde"))

# ---------------------------------------------------------------------------
# Functions and problems

#' Delay differential equation problems
#'
#' `DDEProblem(sys, op, tspan)` builds a delay differential equation problem
#' from a compiled system whose equations contain delayed unknowns
#' `EvalAt(t - tau)(x)`. The values of unknowns in `op` (or their defaults)
#' form the history: they may be expressions of the independent variable,
#' giving `h(p, t)` for `t <= t0`, and the initial state is `h(p, t0)`.
#' Constant lags are detected from the equations; `constant_lags` adds more.
#' Lags that depend on the state or on time (as in `EvalAt(t - 1 - x)(x)`)
#' become dependent lags `lag(u, p, t)`, whose propagated discontinuities
#' are located during the integration. `DDEFunction(sys)` builds
#' `f(u, h, p, t)`.
#'
#' `SDDEProblem(sys, op, tspan)` does the same for systems with brownians.
#'
#' Without a system, `DDEProblem(f, u0, h, tspan, p)` takes an R function
#' `f(u, h, p, t)` and a history `h(p, t)` returning the state vector;
#' `SDDEProblem(f, g, u0, h, tspan, p)` adds the diffusion `g(u, h, p, t)`.
#'
#' Solve with `solve(prob, MethodOfSteps(Tsit5()))` (the default), any ODE
#' algorithm inside [MethodOfSteps()], or `MethodOfSteps(EM())` with `dt`
#' for SDDEs.
#'
#' @param sys A compiled `System`, or the function `f`.
#' @param ... For systems: `op`, `tspan`, then options stored with the
#'   problem. For functions: `u0`, `h`, `tspan`, `p` (SDDE: `g`, `u0`, `h`,
#'   `tspan`, `p`), then options.
#' @param constant_lags Constant lags: discontinuities propagate from the
#'   initial time by these lags and become stopping times.
#' @param dependent_lags A list of functions `lag(u, p, t)` of state
#'   dependent lags (the delayed time is `t - lag(u, p, t)`).
#' @param neutral Whether the equation is neutral (delays in derivative
#'   terms): discontinuities then keep their order when they propagate.
#' @param order_discontinuity_t0 Order of the discontinuity at the initial
#'   time (0: the history and the solution need not even agree).
#' @return A `DDEProblem`, `SDDEProblem` or `DDEFunction`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("y(t)")
#' sys <- mtkcompile(System(list(D(y) ~ -EvalAt(t - 1)(y)), t, name = "dde"))
#' prob <- DDEProblem(sys, list(y ~ 1), c(0, 2))
#' sol <- solve(prob, MethodOfSteps(Tsit5()), reltol = 1e-8, abstol = 1e-10)
#' tail(sol[y], 1)
#' @export
DDEFunction <- function(sys, ...) {
  if (is.function(sys))
    return(structure(list(f = sys, mass_matrix = NULL, sys = NULL), class = c("DDEFunction", "SciMLFunction")))
  .need_complete(sys, "DDEFunction")
  rhs <- .ode_rhs_exprs(sys, substitute_observed = FALSE)
  M <- calculate_massmatrix(sys)
  if (is.matrix(M) && !isTRUE(all.equal(M, diag(nrow(M)))))
    stop("delay systems with algebraic equations are not supported", call. = FALSE)
  structure(list(f = .dde_build(rhs, sys), mass_matrix = NULL, sys = sys),
            class = c("DDEFunction", "SciMLFunction"))
}

#' @rdname DDEFunction
#' @export
DDEProblem <- function(sys, ..., constant_lags = NULL, dependent_lags = NULL, neutral = FALSE,
                       order_discontinuity_t0 = 0L) {
  lag_opts <- list(constant_lags = constant_lags, dependent_lags = dependent_lags, neutral = neutral,
                   order_discontinuity_t0 = order_discontinuity_t0)
  if (is.function(sys) || inherits(sys, "DDEFunction"))
    return(do.call(.dde_problem_raw, c(list("DDEProblem", sys, NULL), list(...), lag_opts)))
  do.call(.dde_problem_sys, c(list("DDEProblem", sys), list(...), lag_opts))
}

#' @rdname DDEFunction
#' @export
SDDEProblem <- function(sys, ..., constant_lags = NULL, dependent_lags = NULL, neutral = FALSE,
                        order_discontinuity_t0 = 0L) {
  lag_opts <- list(constant_lags = constant_lags, dependent_lags = dependent_lags, neutral = neutral,
                   order_discontinuity_t0 = order_discontinuity_t0)
  if (is.function(sys) || inherits(sys, "DDEFunction")) {
    dots <- list(...)
    return(do.call(.dde_problem_raw, c(list("SDDEProblem", sys, dots[[1L]]), dots[-1L], lag_opts)))
  }
  do.call(.dde_problem_sys, c(list("SDDEProblem", sys), list(...), lag_opts))
}

# The lag fields of a problem: constant lags, dependent lags as a list of
# functions lag(u, p, t), neutral and order_discontinuity_t0.
.dde_lag_fields <- function(constant_lags, dependent_lags, neutral, order_discontinuity_t0) {
  dl <- if (is.function(dependent_lags)) list(dependent_lags) else as.list(dependent_lags)
  if (!all(vapply(dl, is.function, TRUE)))
    stop("dependent_lags must be functions lag(u, p, t)", call. = FALSE)
  o0 <- as.integer(order_discontinuity_t0)
  if (length(o0) != 1L || is.na(o0) || o0 < 0L) stop("order_discontinuity_t0 must be a nonnegative integer", call. = FALSE)
  list(constant_lags = as.numeric(unlist(constant_lags)), dependent_lags = dl, neutral = isTRUE(neutral),
       order_discontinuity_t0 = o0)
}

.dde_problem_raw <- function(type, f, g, u0, h, tspan, p = numeric(), ..., constant_lags = NULL,
                             dependent_lags = NULL, neutral = FALSE, order_discontinuity_t0 = 0L,
                             noise_rate_prototype = NULL, noise = NULL) {
  fn <- if (inherits(f, "DDEFunction")) f else DDEFunction(f)
  ntype <- NULL
  if (type == "SDDEProblem") {
    fn$g <- g
    ntype <- if (!is.null(noise_rate_prototype)) "general" else if (identical(noise, "scalar")) "scalar" else "diagonal"
    fn$noise_type <- ntype
    fn$nbrownians <- switch(ntype, general = ncol(as.matrix(noise_rate_prototype)), scalar = 1L,
                            diagonal = length(unlist(u0)))
  }
  lf <- .dde_lag_fields(constant_lags, dependent_lags, neutral, order_discontinuity_t0)
  .new_problem(type, f = fn, u0 = as.numeric(unlist(u0)), h = h, tspan = .as_tspan(tspan), p = p,
               sys = NULL, constant_lags = lf$constant_lags, dependent_lags = lf$dependent_lags,
               neutral = lf$neutral, order_discontinuity_t0 = lf$order_discontinuity_t0, noise_type = ntype,
               kwargs = list(...))
}

# The history h(p, t) from the values of unknowns in op, after
# generate_history.
.dde_history_fn <- function(sys, op) {
  us <- unknowns(sys)
  allv <- c(us, .parameters_of(sys, TRUE), observables(sys))
  m <- smap_merge(bindings(sys), initial_conditions(sys), .resolve_named(as_symmap(op), allv))
  vals <- lapply(us, function(u) {
    v <- smap_get(m, u)
    if (is.null(v) || (!is_sym(v) && length(v) == 1L && is.na(v)))
      .stopf("no initial value or history for unknown %s", format(u))
    as_sym(v)
  })
  if (length(intersect(.keys(.collect_vars(vals)), .keys(us))))
    stop("history values may depend on parameters and time, not on other unknowns", call. = FALSE)
  .build_closure(vals, list(list(name = "p", vars = .problem_params(sys)),
                            list(name = "t", vars = .f(sys, "iv"))),
                 observed = .param_binding_eqs(sys), out = "vector")
}

# Lags t - s of the delayed terms: constant lags (s is t minus an
# expression of parameters) as numbers, and the others (lags that depend on
# unknowns or on time) as functions lag(u, p, t).
.dde_lags <- function(sys, p) {
  iv <- .f(sys, "iv")
  exprs <- .dde_expand_observed(lapply(equations(sys), `[[`, "rhs"), observed(sys), iv)
  terms <- .dde_terms(exprs, iv)
  lags <- list()
  dep <- list()
  for (d in terms) {
    lag <- symcompR::expand(iv - .args(d)[[1L]])
    vs <- .collect_vars(lag)
    if (!any(vapply(vs, function(v) !isparameter(v), TRUE))) {
      lags[[length(lags) + 1L]] <- lag
    } else if (!length(.dde_terms(list(lag), iv)) && !any(.keys(dep) == .key(lag))) {
      # a lag holding delayed values itself is left to the fixed point iteration
      dep[[length(dep) + 1L]] <- lag
    }
  }
  const <- numeric()
  if (length(lags)) {
    fn <- .build_closure(lags, list(list(name = "p", vars = .problem_params(sys))),
                         observed = .param_binding_eqs(sys), out = "vector")
    out <- fn(p)
    const <- unique(out[out > 0])
  }
  dfun <- lapply(dep, function(l) .build_closure(list(l), .codegen_args(sys), observed = .codegen_observed(sys),
                                                 out = "scalar"))
  list(constant = const, dependent = dfun)
}

.dde_problem_sys <- function(type, sys, op = list(), tspan = NULL, ..., constant_lags = NULL,
                             dependent_lags = NULL, neutral = FALSE, order_discontinuity_t0 = 0L) {
  .need_complete(sys, type)
  if (!is_dde(sys)) .stopf("The system has no delays; use %s", if (type == "DDEProblem") "ODEProblem()" else "SDEProblem()")
  if (length(jumps(sys))) .stopf("A system with jumps cannot be used to construct a %s", type)
  if (type == "DDEProblem" && !is.null(.f(sys, "noise_eqs")))
    stop("A system with noise cannot be used to construct a DDEProblem; use SDDEProblem()", call. = FALSE)
  if (type == "SDDEProblem") .sde_check_noise(sys, type)
  tspan <- .as_tspan(tspan, sys)
  vals <- .process_op(sys, op)
  h <- .dde_history_fn(sys, op)
  u0 <- h(vals$p, tspan[1L])
  names(u0) <- names(vals$u0)
  fn <- DDEFunction(sys)
  ntype <- NULL
  if (type == "SDDEProblem") {
    fn$g <- generate_diffusion_function(sys)
    ne <- .sde_noise_eqs(sys)
    ntype <- .sde_noise_type(ne)
    fn$noise_type <- ntype
    fn$nbrownians <- if (is.matrix(ne)) ncol(ne) else length(ne)
  }
  found <- .dde_lags(sys, vals$p)
  lf <- .dde_lag_fields(sort(unique(c(as.numeric(unlist(constant_lags)), found$constant))),
                        c(found$dependent, if (is.function(dependent_lags)) list(dependent_lags) else dependent_lags),
                        neutral, order_discontinuity_t0)
  prob <- .new_problem(type, f = fn, u0 = u0, h = h, tspan = tspan, p = vals$p, sys = sys,
                       constant_lags = lf$constant_lags, dependent_lags = lf$dependent_lags, neutral = lf$neutral,
                       order_discontinuity_t0 = lf$order_discontinuity_t0, noise_type = ntype,
                       kwargs = .problem_kwargs(sys, list(...)))
  prob$op <- vals$op
  prob
}

# ---------------------------------------------------------------------------
# MethodOfSteps

#' Method of steps for delay equations
#'
#' `MethodOfSteps(alg)` turns an ODE algorithm (or `EM()` for stochastic
#' delay equations) into a delay equation solver, as in DelayDiffEq.jl.
#' Delayed values come from the interpolant of earlier steps; when a lag is
#' shorter than the step, the step is repeated with its own interpolant until
#' the fixed point iteration converges. `constrained = TRUE` instead limits
#' the step to the smallest constant lag.
#'
#' Discontinuities are tracked with their order as in DelayDiffEq: the one
#' at the initial time has order `order_discontinuity_t0` of the problem;
#' when the integrator reaches a discontinuity of order `k` (the lowest
#' order among those within rounding of that time), the constant lags
#' propagate it with order `k + 1` (`k` for neutral problems) as long as
#' that is at most the order of the method plus one, and every such time is
#' a stopping time. After each rejected step the dependent lags are checked
#' at `discontinuity_interp_points` points of the attempted step for times
#' `s` with `s - lag(u(s), p, s)` equal to an earlier discontinuity; a root
#' becomes a new discontinuity and stopping time. The options
#' `discontinuity_interp_points` (10), `discontinuity_abstol` (1e-12) and
#' `d_discontinuities` (extra discontinuity times of order 0) are passed to
#' `solve()`.
#'
#' @param alg An ODE algorithm such as [Tsit5()], or an SDE algorithm.
#' @param constrained Limit steps to the smallest constant lag.
#' @param max_fixedpoint_iters Iterations of the fixed point iteration.
#' @return An algorithm object.
#' @examples
#' MethodOfSteps(Tsit5())
#' @export
MethodOfSteps <- function(alg, constrained = FALSE, max_fixedpoint_iters = 10L) {
  structure(list(alg = alg, constrained = constrained, max_fixedpoint_iters = max_fixedpoint_iters,
                 name = paste0("MethodOfSteps(", alg$name, ")")),
            class = c("MethodOfSteps", "DDEAlgorithm"))
}

#' @export
print.DDEAlgorithm <- function(x, ...) {
  cat(x$name, "\n", sep = "")
  invisible(x)
}

.dde_history_store <- function(prob) {
  hs <- new.env(parent = emptyenv())
  hs$steps <- list(); hs$starts <- numeric(); hs$ends <- numeric()
  hs$current <- NULL; hs$inside <- FALSE
  hs$t0 <- prob$tspan[1L]; hs$u0 <- as.numeric(prob$u0); hs$h <- prob$h
  hs
}

# History h(p, s) during integration.
.dde_H <- function(hs) {
  function(p, s) {
    if (s <= hs$t0) return(as.numeric(hs$h(p, s)))
    n <- length(hs$starts)
    if (n && s <= hs$ends[n] + 1e-12 * max(1, abs(s))) {
      i <- max(1L, findInterval(s, hs$starts))
      return(.step_value(hs$steps[[i]], s))
    }
    if (!is.null(hs$current)) return(.step_value(hs$current, s))
    hs$inside <- TRUE
    if (n) .step_value(hs$steps[[n]], s) else hs$u0
  }
}

# A passive callback records each accepted step and handles the
# discontinuities the step lands on.
.dde_recorder <- function(hs, ds) {
  DiscreteCallback(function(u, t, integ) {
    st <- integ$last_step
    hs$steps[[length(hs$steps) + 1L]] <- st
    hs$starts <- c(hs$starts, st$t0)
    hs$ends <- c(hs$ends, st$tend %||% (st$t0 + st$dt))
    .dde_handle_discontinuities(ds, integ$t)
    FALSE
  }, function(integ) NULL, initialize = function(integ) ds$integ <- integ, save_positions = c(FALSE, FALSE))
}

#' Discontinuities of delay problems
#'
#' `Discontinuity(t, order)` marks a discontinuity of the derivative of order
#' `order` at time `t`, as DelayDiffEq's type; pass a list of them (or plain
#' times, of order 0) as `d_discontinuities` to [solve()] for a DDE problem.
#' A discontinuity of order `k` propagates through each constant lag up to
#' the order of the method.
#'
#' @param t A time.
#' @param order The order of the discontinuity.
#' @return A `Discontinuity`.
#' @examples
#' prob <- DDEProblem(function(u, h, p, t) -h(p, t - 1), 1, function(p, t) 1, c(0, 3),
#'                    constant_lags = 1)
#' sol <- solve(prob, MethodOfSteps(Tsit5()), d_discontinuities = list(Discontinuity(0.5, 1L)))
#' sol$retcode
#' @export
Discontinuity <- function(t, order = 0L)
  structure(list(t = as.numeric(t), order = as.integer(order)), class = "Discontinuity")

# The discontinuity state: pending discontinuities (times and orders), the
# tracked (passed) ones that dependent lags propagate, and the problem's
# lags. `max_order` is alg_maximum_order (the method's order).
.dde_discontinuity_state <- function(prob, max_order, opts) {
  ds <- new.env(parent = emptyenv())
  ds$t0 <- prob$tspan[1L]; ds$tf <- prob$tspan[2L]
  ds$lags <- prob$constant_lags %||% numeric()
  ds$dependent <- prob$dependent_lags %||% list()
  ds$neutral <- isTRUE(prob$neutral)
  ds$max_order <- as.integer(max_order)
  ds$interp_points <- as.integer(opts$discontinuity_interp_points %||% 10L)
  ds$abstol <- opts$discontinuity_abstol %||% 1e-12
  ds$pending_t <- numeric(); ds$pending_o <- integer()
  ds$tracked_t <- numeric(); ds$tracked_o <- integer()
  ds$stops <- numeric()
  o0 <- prob$order_discontinuity_t0 %||% 0L
  # user discontinuities: times (order 0) or Discontinuity(t, order)
  dd <- opts$d_discontinuities
  dd <- if (inherits(dd, "Discontinuity")) list(dd) else as.list(dd)
  extra <- vapply(dd, function(d) if (inherits(d, "Discontinuity")) d$t else as.numeric(d), 0)
  ords <- vapply(dd, function(d) if (inherits(d, "Discontinuity")) as.integer(d$order) else 0L, 0L)
  keep <- extra > ds$t0 & extra < ds$tf
  extra <- extra[keep]
  ds$pending_t <- extra; ds$pending_o <- ords[keep]
  ds$stops <- extra
  # initialize_tstops_d_discontinuities_propagated and the tracked t0
  if (o0 <= ds$max_order) {
    .dde_add_next(ds, o0, ds$t0, check = FALSE)
  }
  ds
}

# add_next_discontinuities!: propagate a discontinuity of order `order` at t
# by the constant lags and track it for the dependent lags.
.dde_add_next <- function(ds, order, t, check = TRUE) {
  next_order <- if (ds$neutral) order else order + 1L
  if (check && next_order > ds$max_order + 1L) return(invisible())
  for (lag in ds$lags) if (lag != 0 && lag < ds$tf - t) {
    ds$pending_t <- c(ds$pending_t, t + lag)
    ds$pending_o <- c(ds$pending_o, next_order)
    ds$stops <- c(ds$stops, t + lag)
  }
  ds$tracked_t <- c(ds$tracked_t, t)
  ds$tracked_o <- c(ds$tracked_o, as.integer(order))
  .dde_flush_stops(ds)
}

# New stopping times go to the integrator (add_tstop, R/jumps-aggregators.R).
.dde_flush_stops <- function(ds) {
  if (length(ds$stops) && !is.null(ds$integ)) {
    ds$integ$tstops_added <- c(ds$integ$tstops_added, ds$stops)
    ds$stops <- numeric()
  }
  invisible()
}

# handle_discontinuities!: the discontinuities at t (within 10 eps(t)) are
# removed and the lowest of their orders propagates. Discontinuities passed
# without landing on them (not expected, as they are stopping times)
# propagate from their own times.
.dde_handle_discontinuities <- function(ds, t) {
  if (!length(ds$pending_t)) return(invisible())
  tol <- 10 * .Machine$double.eps * max(abs(t), .Machine$double.xmin)
  at <- abs(ds$pending_t - t) <= tol
  passed <- !at & ds$pending_t < t
  if (!any(at | passed)) return(invisible())
  pt <- ds$pending_t; po <- ds$pending_o
  ds$pending_t <- pt[!(at | passed)]; ds$pending_o <- po[!(at | passed)]
  for (i in which(passed)) .dde_add_next(ds, po[i], pt[i])
  if (any(at)) {
    .dde_add_next(ds, min(po[at]), t)
    # update_fsal!: the next step starts past the discontinuity
    ds$reset_t <- t
  }
  invisible()
}

# track_propagated_discontinuities!: after a rejected step on [t, t + dt]
# with interpolant `ustep(theta)`, look for a time s with
# T + lag(u(s), p, s) - s = 0 for a tracked discontinuity T; the first one
# found becomes a discontinuity of the next order.
.dde_track <- function(ds, t, dt, ustep, p) {
  thetas <- seq(0, 1, length.out = max(2L, ds$interp_points))
  for (lag in ds$dependent) for (j in seq_along(ds$tracked_t)) {
    Td <- ds$tracked_t[j]
    g <- function(th) {
      s <- t + th * dt
      Td + as.numeric(lag(ustep(th), p, s)) - s
    }
    # discontinuity_interval
    g0 <- g(0)
    prev_sign <- if (!is.finite(g0) || abs(g0) <= ds$abstol) 0 else sign(g0)
    g1 <- g(1)
    if (!is.finite(g1)) next
    iv <- NULL
    if (prev_sign * sign(g1) < 0) {
      iv <- c(0, 1)
    } else {
      prev_th <- 0
      for (th in thetas[-1L]) {
        gs <- g(th)
        if (!is.finite(gs)) break
        if (gs == 0) {
          iv <- c(th, th)
          break
        }
        if (prev_sign * sign(gs) < 0) {
          iv <- c(prev_th, th)
          break
        }
        prev_sign <- sign(gs)
        prev_th <- th
      }
    }
    if (is.null(iv)) next
    # discontinuity_time
    th <- if (iv[1L] == iv[2L]) iv[1L] else
      stats::uniroot(g, iv, tol = 1e-14, maxiter = 200L)$root
    s <- t + th * dt
    if (s <= t || s >= ds$tf) next
    ds$pending_t <- c(ds$pending_t, s)
    ds$pending_o <- c(ds$pending_o, if (ds$neutral) ds$tracked_o[j] else ds$tracked_o[j] + 1L)
    ds$stops <- c(ds$stops, s)
    .dde_flush_stops(ds)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

.dde_setup_opts <- function(prob, mos, opts, hs, max_order) {
  lags <- prob$constant_lags %||% numeric()
  t0 <- prob$tspan[1L]; tf <- prob$tspan[2L]
  if (tf <= t0) stop("delay problems need tspan[2] > tspan[1]", call. = FALSE)
  ds <- .dde_discontinuity_state(prob, max_order, opts)
  opts$tstops <- c(unlist(opts$tstops), ds$stops)
  ds$stops <- numeric()
  opts$d_discontinuities <- NULL
  opts$discontinuity_interp_points <- NULL
  opts$discontinuity_abstol <- NULL
  if (isTRUE(mos$constrained) && length(lags)) opts$dtmax <- min(opts$dtmax %||% Inf, min(lags))
  opts$callback <- CallbackSet(.dde_recorder(hs, ds), opts$callback)
  list(opts = opts, ds = ds)
}

#' @rdname MethodOfSteps
#' @param a A `DDEProblem`.
#' @param b `MethodOfSteps(alg)`, or an ODE algorithm (wrapped
#'   automatically). Defaults to `MethodOfSteps(Tsit5())`.
#' @param ... Solver options, as for ODE problems.
#' @export
solve.DDEProblem <- function(a, b = NULL, ...) {
  mos <- if (is.null(b)) MethodOfSteps(Tsit5()) else if (inherits(b, "MethodOfSteps")) b else MethodOfSteps(b)
  inner <- mos$alg
  if (!inherits(inner, "ODEAlgorithm") || inherits(inner, "CompositeAlgorithm"))
    stop("DDE problems need a single ODE algorithm in MethodOfSteps(), such as Tsit5()", call. = FALSE)
  opts <- utils::modifyList(a$kwargs %||% list(), list(...))
  hs <- .dde_history_store(a)
  H <- .dde_H(hs)
  fdde <- a$f$f
  # alg_maximum_order is alg_order (1 for FBDF, as in Julia)
  setup <- .dde_setup_opts(a, mos, opts, hs, inner$order)
  opts <- setup$opts
  ds <- setup$ds
  runner <- inner
  # handle_callback_modifiers!: a state changed by a callback is an order 0
  # discontinuity at that time
  runner$init <- function(ctx, t, u) {
    if (t > ds$t0) {
      ds$pending_t <- c(ds$pending_t, t)
      ds$pending_o <- c(ds$pending_o, 0L)
      .dde_handle_discontinuities(ds, t)
    }
    inner$init(ctx, t, u)
  }
  runner$step <- function(ctx, t, u, dt, cache) {
    extra_nf <- 0L
    ts <- t
    if (identical(ds$reset_t, t)) {
      # shift_past_discontinuity! and reset_fsal!: the step starts just
      # after the discontinuity (its end stays t + dt), so the first stage
      # sees the solution after it
      ts <- t + max(abs(t) * .Machine$double.eps, .Machine$double.xmin)
      if (!is.null(cache$fsal)) {
        cache$fsal <- ctx$f(u, ctx$p, ts)
        extra_nf <- 1L
      }
    }
    res <- .dde_step(ctx, ts, u, dt - (ts - t), cache)
    res$nf <- as.integer(res$nf %||% 0L) + extra_nf
    # loopfooter!: a step the controller rejects is searched for
    # discontinuities propagated by the dependent lags
    if (length(ds$dependent) && isTRUE(ctx$adaptive) && !isTRUE(res$fail)) {
      EEst <- res$EEst %||% .err_norm(res$err, u, res$u, ctx$abstol, ctx$reltol)
      if (!isTRUE(EEst <= 1)) {
        ustep <- function(th) if (th == 0) u else if (th == 1) res$u else inner$interp(th, dt, u, res$u, res$k)
        .dde_track(ds, t, dt, ustep, ctx$p)
      }
    }
    res
  }
  .dde_step <- function(ctx, t, u, dt, cache) {
    hs$current <- NULL
    hs$inside <- FALSE
    res <- inner$step(ctx, t, u, dt, cache)
    if (!hs$inside) return(res)
    converged <- FALSE
    for (it in seq_len(mos$max_fixedpoint_iters)) {
      hs$current <- list(t0 = t, dt = dt, y0 = u, y1 = res$u, k = res$k, interp = inner$interp)
      hs$inside <- FALSE
      res2 <- inner$step(ctx, t, u, dt, cache)
      d <- .err_norm(res2$u - res$u, u, res2$u, ctx$abstol, ctx$reltol)
      res <- res2
      if (is.finite(d) && d < 0.1) {
        converged <- TRUE
        break
      }
    }
    hs$current <- NULL
    # an unconverged iteration fails the step, so the driver retries it shorter
    if (!converged) res$fail <- TRUE
    res
  }
  dprob <- a
  dprob$f <- list(f = function(u, p, t) fdde(u, H, p, t), mass_matrix = NULL)
  sol <- .integrate(dprob, runner, opts)
  .dde_solution(sol, a, mos)
}

#' @rdname MethodOfSteps
#' @export
solve.SDDEProblem <- function(a, b = NULL, ...) {
  mos <- if (is.null(b)) MethodOfSteps(EM()) else if (inherits(b, "MethodOfSteps")) b else MethodOfSteps(b)
  inner <- mos$alg
  if (!inherits(inner, "SDEAlgorithm")) stop("SDDE problems need an SDE algorithm such as EM() in MethodOfSteps()", call. = FALSE)
  opts <- utils::modifyList(a$kwargs %||% list(), list(...))
  seed <- opts$seed
  opts$seed <- NULL
  lags <- a$constant_lags %||% numeric()
  if (length(lags)) {
    opts$dtmax <- min(opts$dtmax %||% Inf, min(lags))
    if (!is.null(opts$dt) && opts$dt > min(lags) * (1 + 1e-12))
      stop("the step dt must not exceed the smallest lag for stochastic delay equations", call. = FALSE)
  }
  hs <- .dde_history_store(a)
  H <- .dde_H(hs)
  # alg_maximum_order of an SDE algorithm: the integer part of its order
  opts <- .dde_setup_opts(a, mos, opts, hs, floor(inner$order))$opts
  fdde <- a$f$f; gdde <- a$f$g
  .with_seed(seed, {
    sol <- .sde_integrate(a, inner, opts, f = function(u, p, t) fdde(u, H, p, t),
                          g = function(u, p, t) gdde(u, H, p, t))
    .as_rode_solution(.dde_solution(sol, a, mos))
  })
}

.dde_solution <- function(sol, prob, mos) {
  env <- .sol_env(sol)
  env$prob <- prob
  env$alg <- mos
  class(sol) <- c("DDESolution", class(sol))
  sol
}

# Values of expressions with delayed terms along a delay solution.
.dde_series <- function(sol, env, exprs, ts) {
  sys <- env$sys
  iv <- .f(sys, "iv")
  args <- .codegen_args(sys)
  cg <- .evalat_codegen(exprs, sys, args, .codegen_observed(sys), out = "vector")
  prob <- env$prob
  t0 <- prob$tspan[1L]
  p <- env$p_final
  S <- function(s) if (s <= t0) as.numeric(prob$h(p, s)) else as.numeric(sol(s))
  out <- matrix(NA_real_, length(ts), length(exprs))
  for (i in seq_along(ts)) {
    u <- as.numeric(sol(ts[i]))
    vals <- if (cg$n) .evalat_values(cg$tfun(u, p, ts[i]), cg$idx, S) else numeric()
    out[i, ] <- cg$fn(u, p, ts[i], vals)
  }
  out
}

.dde_has_delay <- function(i, sys) {
  if (is.null(sys)) return(FALSE)
  xs <- if (is.list(i) && !is_sym(i)) i else list(i)
  if (!all(vapply(xs, function(z) is_sym(z) || inherits(z, "Equation"), TRUE))) return(FALSE)
  obs <- observed(sys)
  rule <- if (length(obs)) .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs")) else NULL
  ex <- lapply(xs, function(z) if (is.null(rule)) as_sym(z) else .subst_fixpoint(as_sym(z), rule))
  ex <- .dde_expand_observed(ex, obs, .f(sys, "iv"))
  length(.dde_terms(ex, .f(sys, "iv"))) > 0L
}

#' @export
`[.DDESolution` <- function(x, i, j, ...) {
  env <- .sol_env(x)
  if (missing(i) || !.dde_has_delay(i, env$sys)) return(NextMethod())
  rows <- if (missing(j)) seq_along(env$t) else j
  ts <- env$t[rows]
  if (is.list(i) && !is_sym(i)) {
    out <- .dde_series(x, env, lapply(i, as_sym), ts)
    colnames(out) <- vapply(i, format, "")
    return(out)
  }
  .dde_series(x, env, list(as_sym(i)), ts)[, 1L]
}
