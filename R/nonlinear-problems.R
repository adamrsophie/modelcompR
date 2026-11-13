# Nonlinear, interval, linear, steady state and SCC problems built from
# systems, after ModelingToolkitBase/src/problems/nonlinearproblem.jl
# (NonlinearFunction, NonlinearProblem, NonlinearLeastSquaresProblem,
# generate_nonlinear_bounds), intervalnonlinearproblem.jl, linearproblem.jl
# (LinearFunction, LinearProblem, calculate_A_b), odeproblem.jl
# (SteadyStateProblem), systems/system.jl (NonlinearSystem(sys::System)) and
# ModelingToolkit/src/problems/sccnonlinearproblem.jl (SCCNonlinearProblem).
#
# Problems are the S3 lists of R/problems.R: f (a function object), u0 and p
# (named numeric vectors), sys and kwargs. Generated residuals are f(u, p).

# ---------------------------------------------------------------------------
# Shared helpers

# Keys of the declared variables an expression depends on. Applications of
# operators such as Initial(x) are atoms here: Initial(x) is a parameter, not
# an occurrence of x, and x(t) does not count as an occurrence of t.
.nl_atomic_keys <- function(exprs) {
  seen <- new.env(parent = emptyenv())
  out <- character()
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    k <- .key(z)
    if (!is.null(seen[[k]])) return(invisible())
    assign(k, TRUE, envir = seen)
    ty <- .type(z)
    if (ty == "Derivative" || (ty %in% c("Symbol", "Function") && !is.null(.var_kind(z)))) {
      out[length(out) + 1L] <<- k
      return(invisible())
    }
    for (a in .args(z)) walk(a)
  }
  if (is_sym(exprs)) walk(exprs) else for (z in exprs) walk(z)
  out
}

# Operator applications (Initial(x), Pre(x)) inside expressions.
.nl_op_nodes <- function(exprs) {
  seen <- new.env(parent = emptyenv())
  out <- list()
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    k <- .key(z)
    if (!is.null(seen[[k]])) return(invisible())
    assign(k, TRUE, envir = seen)
    if (.type(z) == "Function" && identical(.var_kind(z), "operator")) {
      out[[length(out) + 1L]] <<- z
      return(invisible())
    }
    for (a in .args(z)) walk(a)
  }
  if (is_sym(exprs)) walk(exprs) else for (z in exprs) walk(z)
  out
}

# Substitute variables, leaving operator applications untouched: the rule
# maps every Initial(x) to itself so the top down walker does not rewrite the
# x inside it.
.nl_subst <- function(exprs, from, to) {
  ops <- .nl_op_nodes(exprs)
  rule <- .rule(c(ops, from), c(ops, to))
  if (is_sym(exprs)) .subst(exprs, rule) else lapply(exprs, .subst, rule = rule)
}

.nl_check_system <- function(sys, what, time_dependent = FALSE, cost = FALSE) {
  .need_complete(sys, what)
  iv <- .f(sys, "iv")
  if (time_dependent && is.null(iv))
    .stopf("`%s` requires a time dependent system; `%s` has no independent variable.", what, nameof(sys))
  if (!time_dependent && !is.null(iv))
    .stopf("`%s` requires a time independent system; `%s` depends on %s.", what, nameof(sys), format(iv))
  if (isTRUE(.f(sys, "is_dde"))) .stopf("`%s` does not support systems with delays.", what)
  if (!cost && length(costs(sys)))
    .stopf("`%s` does not support systems with costs; use `OptimizationProblem`.", what)
  if (cost && !length(costs(sys))) .stopf("`%s` requires a system with a cost.", what)
  if (!cost && length(constraints(sys)))
    .stopf("`%s` does not support systems with constraints.", what)
  if (length(jumps(sys))) .stopf("`%s` does not support systems with jumps.", what)
  if (!is.null(.f(sys, "noise_eqs")) || length(.f(sys, "brownians")))
    .stopf("`%s` does not support systems with noise.", what)
  invisible(TRUE)
}

.nl_names <- function(vars) vapply(vars, format, "")

# A deterministic value in [0, 1) from a variable name, the counterpart of
# hash(var) / 2^64 in MissingGuessValue.HashedRandom.
.nl_hashed_random <- function(name) {
  h <- 5381
  for (code in utf8ToInt(name)) h <- (h * 33 + code) %% 4294967296
  h / 4294967296
}

#' Fallback values for missing guesses
#'
#' `MissingGuessValue` holds the constructors of the choices for the
#' `missing_guess_value` argument of the problem constructors, which says
#' what happens to an unknown of a nonlinear solve that has neither a value
#' nor a guess: `MissingGuessValue$Constant(val)` uses `val`,
#' `MissingGuessValue$Random()` a uniform random number,
#' `MissingGuessValue$HashedRandom()` a reproducible pseudo random number
#' computed from the variable name (the default for initialization), and
#' `MissingGuessValue$Error()` stops with an error.
#'
#' @format A list of four constructor functions.
#' @examples
#' MissingGuessValue$Constant(0)
#' MissingGuessValue$HashedRandom()
#' @export
MissingGuessValue <- list(
  Constant = function(val) structure(list(kind = "Constant", val = as.numeric(val)), class = "MissingGuessValue"),
  Random = function(rng = NULL) structure(list(kind = "Random", rng = rng), class = "MissingGuessValue"),
  HashedRandom = function() structure(list(kind = "HashedRandom"), class = "MissingGuessValue"),
  Error = function() structure(list(kind = "Error"), class = "MissingGuessValue"))

#' @export
print.MissingGuessValue <- function(x, ...) {
  cat("MissingGuessValue.", x$kind, if (x$kind == "Constant") paste0("(", x$val, ")") else "()", "\n", sep = "")
  invisible(x)
}

.nl_missing_values <- function(vars, mgv, initialization = FALSE) {
  if (is.null(mgv)) mgv <- MissingGuessValue$Error()
  if (is.numeric(mgv)) mgv <- MissingGuessValue$Constant(mgv)
  nms <- .nl_names(vars)
  switch(mgv$kind,
    Constant = rep(mgv$val, length(vars)),
    Random = if (is.function(mgv$rng)) vapply(vars, function(v) mgv$rng(), 0) else stats::runif(length(vars)),
    HashedRandom = vapply(nms, .nl_hashed_random, 0, USE.NAMES = FALSE),
    Error = if (initialization)
      stop("Cyclic guesses detected in the system. No numeric guess could be found for: ",
           paste(nms, collapse = ", "), ".\nIn order to resolve this, please provide additional numeric ",
           "guesses so that the chain can be resolved to assign numeric values to each variable. ",
           "Alternatively, the `missing_guess_value` argument can be used to set a fallback guess ",
           "for all variables; see `MissingGuessValue`.", call. = FALSE)
    else
      stop("Initial condition underdefined. Some are missing from the variable map.\n",
           "Please provide a default (`u0`), initialization equation, or guess\n",
           "for the following variables:\n", paste(nms, collapse = ", "), call. = FALSE))
}

# Numeric values of `vars` taken from guesses, evaluating guess expressions
# against the already known values.
.nl_eval_guesses <- function(g, vars, known_vars = list(), known_vals = numeric(), defs = list()) {
  m <- symmap()
  # bound parameters and observed variables, which guesses may refer to
  for (e in defs) m[[.key(e$lhs)]] <- list(var = e$lhs, value = e$rhs)
  for (el in unclass(g)) if (!is.null(el$var) && !is.null(el$value)) m[[.key(el$var)]] <- el
  for (i in seq_along(known_vars))
    if (!is.na(known_vals[i])) m[[.key(known_vars[[i]])]] <- list(var = known_vars[[i]], value = known_vals[[i]])
  class(m) <- "mtk_symmap"
  vals <- .evaluate_varmap(m)
  vapply(vars, function(v) {
    z <- vals[[.key(v)]]
    if (is.null(z) || !is.numeric(z) || length(z) != 1L) NA_real_ else z
  }, 0)
}

# u0 and p of a time independent problem: values from the operating point
# and initial conditions first, then guesses, then `missing_guess_value`.
.nl_values <- function(sys, op, guesses = NULL, missing_guess_value = NULL, require_p = TRUE) {
  vals <- .process_op(sys, op, require_p = require_p)
  u0 <- vals$u0
  if (anyNA(u0)) u0 <- .nl_fill_guesses(sys, u0, vals$p, guesses, missing_guess_value)
  vals$u0 <- u0
  vals
}

.nl_fill_guesses <- function(sys, u0, p, guesses = NULL, missing_guess_value = NULL) {
  us <- unknowns(sys)
  miss <- which(is.na(u0))
  if (!length(miss)) return(u0)
  ps <- .problem_params(sys)
  g <- smap_merge(guesses(sys), .resolve_named(as_symmap(guesses), c(us, ps, observables(sys))))
  gv <- .nl_eval_guesses(g, us[miss], c(us[-miss], ps), c(u0[-miss], p), .codegen_observed(sys))
  if (anyNA(gv)) gv[is.na(gv)] <- .nl_missing_values(us[miss][is.na(gv)], missing_guess_value)
  u0[miss] <- gv
  u0
}

# Lower and upper bounds of the unknowns from their `bounds` metadata;
# NULL when no unknown is bounded (generate_nonlinear_bounds).
.nl_bounds <- function(sys, lb = NULL, ub = NULL, vals = NULL) {
  us <- unknowns(sys)
  if (!is.null(lb) || !is.null(ub)) {
    lb <- if (is.null(lb)) rep(-Inf, length(us)) else as.numeric(lb)
    ub <- if (is.null(ub)) rep(Inf, length(us)) else as.numeric(ub)
    if (length(lb) != length(us) || length(ub) != length(us))
      stop("`lb` and `ub` must have one entry per unknown", call. = FALSE)
    return(list(lb = lb, ub = ub))
  }
  bs <- lapply(us, getbounds)
  num <- function(b) {
    if (is.numeric(b)) return(as.numeric(b))
    # symbolic bounds are resolved with the operating point values
    b <- as_sym(b)
    if (!is.null(vals)) {
      vs <- .collect_vars(b)
      known <- vapply(vs, function(v) vals[[.key(v)]] %||% NA_real_, 0)
      if (length(vs) && !anyNA(known)) b <- .subst(b, .rule(vs, lapply(known, as_sym)))
    }
    as.numeric(N(b))
  }
  lo <- vapply(bs, function(b) num(b[[1L]]), 0)
  hi <- vapply(bs, function(b) num(b[[2L]]), 0)
  symbolic <- any(vapply(bs, function(b) is_sym(b[[1L]]) || is_sym(b[[2L]]), TRUE))
  if (all(lo == -Inf) && all(hi == Inf)) return(list(lb = NULL, ub = NULL, symbolic = symbolic))
  list(lb = lo, ub = hi, symbolic = symbolic)
}

# Bounds that are expressions of parameters follow the parameters on remake.
.nl_refresh_bounds <- function(prob) {
  if (!isTRUE(prob$symbolic_bounds)) return(prob)
  sys <- prob$sys
  vars <- c(unknowns(sys), .problem_params(sys))
  b <- .nl_bounds(sys, NULL, NULL, stats::setNames(as.list(c(prob$u0, prob$p)), .keys(vars)))
  prob$lb <- b$lb
  prob$ub <- b$ub
  prob
}

# NonlinearSystem(sys::System): the steady state system of a time dependent
# system. Derivatives of the unknowns become zero and the independent
# variable a parameter with the value Inf.
.nl_from_time_dependent <- function(sys) {
  iv <- .f(sys, "iv")
  us <- unknowns(sys)
  rule <- .rule(lapply(us, .make_derivative, iv = iv, order = 1L), rep(list(.S0()), length(us)))
  eqs <- lapply(equations(sys), function(e) {
    r <- if (.is_derivative(e$lhs)) e$rhs else e$rhs - e$lhs
    Equation(.S0(), .subst(r, rule))
  })
  obs <- observed(sys)
  ps <- .parameters_of(sys, TRUE)
  binds <- bindings(sys)
  ics <- initial_conditions(sys)
  pb <- symmap()
  var_keys <- .keys(c(us, lapply(obs, `[[`, "lhs")))
  for (k in names(binds)) {
    el <- binds[[k]]
    if (is.null(el$var)) next
    # a variable binding is an initialization constraint, which a time
    # independent system cannot enforce; it becomes an initial value
    if (k %in% var_keys) {
      if (is.null(ics[[k]])) ics[[k]] <- el
    } else pb[[k]] <- el
  }
  class(pb) <- "mtk_symmap"
  if (.key(iv) %in% .nl_atomic_keys(c(lapply(eqs, `[[`, "rhs"), lapply(obs, `[[`, "rhs")))) {
    ps <- c(ps, list(iv))
    ics <- smap_set(ics, iv, Inf)
  }
  nsys <- System(eqs, NULL, unknowns = us, ps = ps, observed = obs, bindings = pb,
                 initial_conditions = ics, guesses = guesses(sys), name = nameof(sys),
                 checks = FALSE, discover_from_metadata = FALSE)
  complete(.update_system(nsys, isscheduled = TRUE, metadata = .f(sys, "metadata")))
}

# ---------------------------------------------------------------------------
# NonlinearFunction and NonlinearProblem

#' Nonlinear functions and problems
#'
#' `NonlinearProblem(sys, op)` builds the root finding problem
#' `f(u, p) = 0` from a complete time independent system: the unknowns `u`
#' start from the values in `op` (formulas `x ~ 1`, pairs, or a named list),
#' the initial conditions of the system, and its guesses, in that order.
#' `NonlinearLeastSquaresProblem()` is the same for systems with more
#' equations than unknowns, solved in the least squares sense.
#' `NonlinearFunction(sys)` builds the residual and, with `jac = TRUE`, the
#' symbolic Jacobian. A time dependent system is first turned into its steady
#' state system (derivatives set to zero).
#'
#' `NonlinearProblem(f, u0, p)` with an R function `f(u, p)` builds a problem
#' without a system.
#'
#' Bounds are taken from the `bounds` metadata of the unknowns unless `lb` or
#' `ub` is given; iterates are projected onto the box.
#'
#' @param sys A complete `System`, an R function `f(u, p)`, or a
#'   `NonlinearFunction`.
#' @param op Operating point: values of unknowns and parameters.
#' @param ... Solver options stored with the problem (`abstol`, `maxiters`,
#'   ...).
#' @param jac Generate the symbolic Jacobian. Without it solvers use finite
#'   differences.
#' @param u0,p Initial values and parameters for raw problems.
#' @param lb,ub Lower and upper bounds of the unknowns.
#' @param guesses Further guesses, as a variable map.
#' @param missing_guess_value What to do with unknowns that have no value
#'   and no guess; see [MissingGuessValue].
#' @param check_length Require as many equations as unknowns.
#' @return A `NonlinearFunction`, `NonlinearProblem` or
#'   `NonlinearLeastSquaresProblem`.
#' @examples
#' x <- Variable("x"); y <- Variable("y"); a <- Parameter("a")
#' sys <- complete(System(list(0 ~ x^2 + y^2 - a, 0 ~ x - y), name = "ns"))
#' prob <- NonlinearProblem(sys, list(x ~ 1, y ~ 2, a ~ 2))
#' sol <- solve(prob)
#' sol$u
#' sol[x + y]
#' # without a system
#' sol2 <- solve(NonlinearProblem(function(u, p) u^2 - p, 1, 2))
#' sol2$u
#' @export
NonlinearFunction <- function(sys, jac = FALSE, ...) {
  dots <- list(...)
  if (is.function(sys)) {
    return(structure(list(f = sys, jac = if (is.function(jac)) jac else dots$jac_fn, sys = NULL,
                          observed = NULL, resid_prototype = dots$resid_prototype),
                     class = c("NonlinearFunction", "SciMLFunction")))
  }
  if (inherits(sys, "System") && !is.null(.f(sys, "iv")) && iscomplete(sys)) sys <- .nl_from_time_dependent(sys)
  .nl_check_system(sys, "NonlinearFunction")
  if (.nl_has_call(sys, "limited"))
    stop("the system contains `limited(...)` nodes that were not lowered; build a `NonlinearProblem`, ",
         "which lowers them.", call. = FALSE)
  args <- .codegen_args(sys, with_t = FALSE)
  f <- .build_closure(.ode_rhs_exprs(sys, substitute_observed = FALSE), args,
                      observed = .codegen_observed(sys), out = "vector")
  J <- NULL
  if (isTRUE(jac)) {
    Jm <- .partials(.ode_rhs_exprs(sys), unknowns(sys), .nl_hold(sys))
    J <- .build_closure(Jm, args, out = "matrix")
  }
  structure(list(f = f, jac = J, sys = sys,
                 observed = function(expr) build_explicit_observed_function(sys, expr),
                 resid_prototype = numeric(length(equations(sys)))),
            class = c("NonlinearFunction", "SciMLFunction"))
}

# Everything that must stay fixed when differentiating with respect to the
# unknowns: parameters, including Initial(x) parameters, and operator nodes.
.nl_hold <- function(sys) {
  ps <- .parameters_of(sys, TRUE)
  ops <- .nl_op_nodes(c(lapply(equations(sys), `[[`, "rhs"), lapply(equations(sys), `[[`, "lhs"),
                        lapply(.codegen_observed(sys), `[[`, "rhs")))
  .unique_syms(c(ps, ops, independent_variables(sys)))
}

# Whether a registered function such as homotopy or limited is applied
# anywhere in the equations or observed equations of a system.
.nl_has_call <- function(sys, fname) {
  exprs <- c(lapply(equations(sys), `[[`, "lhs"), lapply(equations(sys), `[[`, "rhs"),
             lapply(observed(sys), `[[`, "rhs"))
  .nl_exprs_have_call(exprs, fname)
}

.nl_exprs_have_call <- function(exprs, fname) {
  seen <- new.env(parent = emptyenv())
  found <- FALSE
  walk <- function(z) {
    if (found || !is_sym(z)) return(invisible())
    k <- .key(z)
    if (!is.null(seen[[k]])) return(invisible())
    assign(k, TRUE, envir = seen)
    if (.type(z) == "Function" && identical(.name(z), fname)) {
      found <<- TRUE
      return(invisible())
    }
    for (a in .args(z)) walk(a)
  }
  for (z in exprs) walk(z)
  found
}

.nl_raw_problem <- function(type, sys, u0, p, op, dots, lb = NULL, ub = NULL) {
  f <- if (inherits(sys, "NonlinearFunction")) sys else NonlinearFunction(sys)
  if (is.null(u0)) u0 <- op
  prob <- .new_problem(type, f = f, u0 = unlist(u0), p = if (is.null(p)) numeric() else p,
                       sys = NULL, kwargs = dots)
  prob$lb <- lb
  prob$ub <- ub
  class(prob) <- c(type, "AbstractNonlinearProblem", "SciMLProblem")
  prob
}

.nl_build_problem <- function(type, sys, op, dots, jac, lb, ub, guesses, missing_guess_value,
                              check_length) {
  if (inherits(sys, "System") && !is.null(.f(sys, "iv")) && iscomplete(sys)) sys <- .nl_from_time_dependent(sys)
  .nl_check_system(sys, type)
  post <- NULL
  if (.nl_has_call(sys, "limited")) {
    low <- .lower_limited(sys)
    sys <- low$sys
    post <- low$postcondition
  }
  neq <- length(equations(sys)); nu <- length(unknowns(sys))
  if (check_length && neq != nu)
    .stopf(paste0("`%s` requires as many equations as unknowns; the system has %d equations and %d ",
                  "unknowns. Use `NonlinearLeastSquaresProblem` for a non square system."), type, neq, nu)
  fn <- NonlinearFunction(sys, jac = jac)
  vals <- .nl_values(sys, op, guesses, missing_guess_value, require_p = FALSE)
  pinit <- .nl_solve_parameters(sys, vals, dots)
  vals$p <- pinit$p
  b <- .nl_bounds(sys, lb, ub, vals$values)
  if (!is.null(post)) dots$postcondition <- post
  prob <- .new_problem(type, f = fn, u0 = vals$u0, p = vals$p, sys = sys, kwargs = dots)
  prob$op <- vals$op
  prob$lb <- b$lb
  prob$ub <- b$ub
  prob$symbolic_bounds <- isTRUE(b$symbolic)
  prob$guesses <- guesses
  prob$missing_guess_value <- missing_guess_value
  prob$initialization_data <- pinit$data
  class(prob) <- c(if (identical(pinit$data$retcode, "InitialFailure")) "InitializationFailure",
                   type, "AbstractNonlinearProblem", "SciMLProblem")
  prob
}

# Parameters without a value: solvable parameters (binding NA) are found by
# the time independent initialization system when the initialization module
# is loaded; anything else is an error.
.nl_solve_parameters <- function(sys, vals, dots = list()) {
  p <- vals$p
  if (!anyNA(p)) return(list(p = p, data = NULL))
  h <- .mtk$hooks$initialize_parameters
  if (!is.null(h)) return(h(sys, vals, dots))
  stop("no value for parameter(s) ", paste(names(p)[is.na(p)], collapse = ", "),
       "; give them in the operating point or as defaults", call. = FALSE)
}

#' @rdname NonlinearFunction
#' @export
NonlinearProblem <- function(sys, op = list(), ..., jac = TRUE, u0 = NULL, p = NULL, lb = NULL,
                             ub = NULL, guesses = list(), missing_guess_value = NULL,
                             check_length = TRUE) {
  dots <- list(...)
  if (inherits(sys, "NonlinearFunction") && !is.null(sys$sys)) sys <- sys$sys
  if (is.function(sys) || inherits(sys, "NonlinearFunction")) {
    if (is.null(p) && length(dots) && (is.null(names(dots)) || !nzchar(names(dots)[1L]))) {
      p <- dots[[1L]]
      dots <- dots[-1L]
    }
    return(.nl_raw_problem("NonlinearProblem", sys, u0, p, op, dots, lb, ub))
  }
  .nl_build_problem("NonlinearProblem", sys, op, dots, jac, lb, ub, guesses, missing_guess_value,
                    check_length)
}

#' @rdname NonlinearFunction
#' @export
NonlinearLeastSquaresProblem <- function(sys, op = list(), ..., jac = TRUE, u0 = NULL, p = NULL,
                                         lb = NULL, ub = NULL, guesses = list(),
                                         missing_guess_value = NULL, check_length = FALSE) {
  dots <- list(...)
  if (inherits(sys, "NonlinearFunction") && !is.null(sys$sys)) sys <- sys$sys
  if (is.function(sys) || inherits(sys, "NonlinearFunction")) {
    if (is.null(p) && length(dots) && (is.null(names(dots)) || !nzchar(names(dots)[1L]))) {
      p <- dots[[1L]]
      dots <- dots[-1L]
    }
    return(.nl_raw_problem("NonlinearLeastSquaresProblem", sys, u0, p, op, dots, lb, ub))
  }
  .nl_build_problem("NonlinearLeastSquaresProblem", sys, op, dots, jac, lb, ub, guesses,
                    missing_guess_value, check_length)
}

# ---------------------------------------------------------------------------
# SteadyStateProblem

#' Steady state problems
#'
#' `SteadyStateProblem(sys, op)` asks for the steady state of a time
#' dependent system: the state where every derivative vanishes.
#' `solve(prob, SSRootfind())` finds it by root finding from the initial
#' values, `solve(prob, DynamicSS())` by integrating the differential
#' equations until the derivatives are small.
#'
#' @param sys A complete time dependent `System`, or an R function
#'   `f(u, p, t)`.
#' @param op Operating point: values of unknowns and parameters.
#' @param ... Passed to [ODEProblem()].
#' @param u0,p Initial values and parameters for raw problems.
#' @return A `SteadyStateProblem`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); k <- Parameter("k"); s <- Parameter("s")
#' sys <- mtkcompile(System(list(D(x) ~ s - k * x), t, name = "decay"))
#' prob <- SteadyStateProblem(sys, list(x ~ 1, k ~ 2, s ~ 3))
#' solve(prob, SSRootfind())$u
#' solve(prob, DynamicSS())$u
#' @export
SteadyStateProblem <- function(sys, op = list(), ..., u0 = NULL, p = NULL) {
  if (inherits(sys, "ODEProblem")) {
    prob <- sys
  } else if (is.function(sys) || inherits(sys, "ODEFunction")) {
    prob <- ODEProblem(sys, op, c(0, Inf), ..., u0 = u0, p = p)
  } else {
    .nl_check_system(sys, "SteadyStateProblem", time_dependent = TRUE)
    prob <- ODEProblem(sys, op, c(0, Inf), ...)
  }
  prob$tspan <- c(prob$tspan[1L], Inf)
  class(prob) <- c(if (inherits(prob, "InitializationFailure")) "InitializationFailure", "SteadyStateProblem",
                   "SciMLProblem")
  prob
}

# ---------------------------------------------------------------------------
# Interval problems

#' Interval nonlinear problems
#'
#' `IntervalNonlinearProblem(sys, uspan, parammap)` builds the problem of
#' finding the root of a single equation in a single unknown inside the
#' bracketing interval `uspan`. It is solved by the bracketing methods
#' [ITP()], [Bisection()], [Brent()] and [Falsi()].
#' `IntervalNonlinearFunction(sys)` builds the scalar residual `f(u, p)`.
#'
#' @param sys A complete `System` with one equation and one unknown, or an R
#'   function `f(u, p)`.
#' @param uspan The interval `c(left, right)`; the residual must change sign
#'   over it.
#' @param parammap Parameter values.
#' @param ... Solver options stored with the problem.
#' @param p Parameters for raw problems.
#' @return An `IntervalNonlinearFunction` or `IntervalNonlinearProblem`.
#' @examples
#' x <- Variable("x"); a <- Parameter("a")
#' sys <- complete(System(list(0 ~ x^2 - a), name = "sq"))
#' prob <- IntervalNonlinearProblem(sys, c(0, 3), list(a ~ 2))
#' solve(prob)$u
#' solve(prob, Brent())$u
#' @export
IntervalNonlinearFunction <- function(sys, ...) {
  if (is.function(sys))
    return(structure(list(f = sys, sys = NULL), class = c("IntervalNonlinearFunction", "SciMLFunction")))
  .nl_check_interval(sys, "IntervalNonlinearFunction")
  f <- .build_closure(.ode_rhs_exprs(sys, substitute_observed = FALSE)[[1L]],
                      .codegen_args(sys, with_t = FALSE), observed = .codegen_observed(sys),
                      out = "scalar")
  structure(list(f = f, sys = sys, observed = function(expr) build_explicit_observed_function(sys, expr)),
            class = c("IntervalNonlinearFunction", "SciMLFunction"))
}

.nl_check_interval <- function(sys, what) {
  .nl_check_system(sys, what)
  if (length(unknowns(sys)) != 1L)
    .stopf("`%s` requires a system with a single unknown. Found %s.", what,
           paste(.nl_names(unknowns(sys)), collapse = ", "))
  if (length(equations(sys)) != 1L)
    .stopf("`%s` requires a system with a single equation. Found %d.", what, length(equations(sys)))
  invisible(TRUE)
}

#' @rdname IntervalNonlinearFunction
#' @export
IntervalNonlinearProblem <- function(sys, uspan, parammap = list(), ..., p = NULL) {
  uspan <- as.numeric(unlist(uspan))
  if (length(uspan) != 2L) stop("`uspan` must be c(left, right)", call. = FALSE)
  if (is.function(sys) || inherits(sys, "IntervalNonlinearFunction")) {
    f <- if (is.function(sys)) IntervalNonlinearFunction(sys) else sys
    if (is.null(p)) p <- if (is.numeric(parammap)) parammap else numeric()
    prob <- .new_problem("IntervalNonlinearProblem", f = f, u0 = uspan[1L], tspan = uspan, p = p,
                         sys = NULL, kwargs = list(...))
    prob$uspan <- uspan
    return(prob)
  }
  .nl_check_interval(sys, "IntervalNonlinearProblem")
  us <- unknowns(sys)
  op <- smap_set(.resolve_named(as_symmap(parammap), .parameters_of(sys, TRUE)), us[[1L]], uspan[1L])
  vals <- .process_op(sys, op)
  prob <- .new_problem("IntervalNonlinearProblem", f = IntervalNonlinearFunction(sys), u0 = vals$u0,
                       tspan = uspan, p = vals$p, sys = sys, kwargs = list(...))
  prob$uspan <- uspan
  prob$op <- vals$op
  prob
}

# ---------------------------------------------------------------------------
# Linear problems

# calculate_A_b: for a system whose equations are affine in the unknowns,
# A u = b with A the Jacobian and b minus the residual at u = 0. Returns NULL
# (or stops, with throw = TRUE) when the system is not affine.
.nl_calculate_A_b <- function(sys, throw = TRUE) {
  us <- unknowns(sys)
  resid <- .ode_rhs_exprs(sys)
  A <- .partials(resid, us, .nl_hold(sys))
  ukeys <- .keys(us)
  for (e in A) {
    if (.is_num(e)) next
    if (any(.nl_atomic_keys(e) %in% ukeys)) {
      if (throw) stop("the system is not affine in its unknowns: the Jacobian entry ", format(e),
                      " depends on them", call. = FALSE)
      return(NULL)
    }
  }
  b <- lapply(.nl_subst(resid, us, rep(list(.S0()), length(us))), function(r) -r)
  list(A = A, b = b)
}

#' Linear problems
#'
#' `LinearProblem(sys, op)` builds `A u = b` from a complete time independent
#' system whose equations are affine in the unknowns: `A` is the Jacobian of
#' the equations and `b` their constant part, both evaluated at the
#' parameter values. `solve()` uses a QR factorization, and least squares for
#' non square systems. `LinearProblem(A, b)` takes a numeric matrix and
#' vector.
#'
#' @param sys A complete `System`, or a numeric matrix `A`.
#' @param op Operating point: parameter values (and optionally initial
#'   values of the unknowns), or the vector `b` for the matrix form.
#' @param ... Options stored with the problem.
#' @param check_length Require a square system.
#' @return A `LinearProblem` with fields `A`, `b`, `u0` and `p`.
#' @examples
#' x <- Variable("x"); y <- Variable("y"); a <- Parameter("a")
#' sys <- complete(System(list(0 ~ 2 * x + y - a, 0 ~ x - y), name = "lin"))
#' prob <- LinearProblem(sys, list(a ~ 3))
#' prob$A; prob$b
#' solve(prob)$u
#' @export
LinearProblem <- function(sys, op = list(), ..., check_length = TRUE) {
  if ((is.matrix(sys) && is.numeric(sys)) || inherits(sys, c("LowRankUpdatedMatrix", "AbstractSciMLOperator")) ||
      is.function(sys)) {
    prob <- .new_problem("LinearProblem", f = NULL, u0 = NULL, p = numeric(), sys = NULL, kwargs = list(...))
    prob$A <- sys
    prob$b <- as.numeric(unlist(op))
    return(prob)
  }
  .nl_check_system(sys, "LinearProblem")
  us <- unknowns(sys)
  if (check_length && length(equations(sys)) != length(us))
    .stopf("`LinearProblem` requires as many equations as unknowns; the system has %d equations and %d unknowns.",
           length(equations(sys)), length(us))
  Ab <- .nl_calculate_A_b(sys)
  pargs <- list(list(name = "p", vars = .problem_params(sys)))
  obs <- .codegen_observed(sys)
  fn <- structure(list(update_A = .build_closure(Ab$A, pargs, observed = obs, out = "matrix"),
                       update_b = .build_closure(Ab$b, pargs, observed = obs, out = "vector"),
                       A = Ab$A, b = Ab$b, sys = sys),
                  class = c("LinearFunction", "SciMLFunction"))
  vals <- .process_op(sys, op)
  prob <- .new_problem("LinearProblem", f = fn, u0 = if (anyNA(vals$u0)) NULL else vals$u0,
                       p = vals$p, sys = sys, kwargs = list(...))
  prob$op <- vals$op
  prob$A <- fn$update_A(vals$p)
  prob$b <- as.numeric(fn$update_b(vals$p))
  prob
}

# ---------------------------------------------------------------------------
# SCC decomposition

# Maximum matching of equations to unknowns by augmenting paths. `adj[[i]]`
# holds the unknown indices incident to equation i. Returns eq_of_var
# (0 when unmatched) and var_of_eq.
.nl_matching <- function(adj, nvar) {
  neq <- length(adj)
  eq_of_var <- integer(nvar)
  var_of_eq <- integer(neq)
  try_eq <- function(i, visited) {
    for (j in adj[[i]]) {
      if (visited[[as.character(j)]] %||% FALSE) next
      visited[[as.character(j)]] <- TRUE
      if (eq_of_var[j] == 0L || try_eq(eq_of_var[j], visited)) {
        eq_of_var[j] <<- i
        var_of_eq[i] <<- j
        return(TRUE)
      }
    }
    FALSE
  }
  for (i in seq_len(neq)) try_eq(i, new.env(parent = emptyenv()))
  list(eq_of_var = eq_of_var, var_of_eq = var_of_eq)
}

# Incidence of equations on unknowns, operator applications being atoms.
.nl_incidence <- function(exprs, vars) {
  vk <- .keys(vars)
  lapply(exprs, function(e) which(vk %in% .nl_atomic_keys(e)))
}

# Strongly connected components of the matched system in topological order
# (Tarjan). Node i is equation i with its matched unknown; i depends on j
# when equation i contains the unknown matched to equation j.
.nl_sccs <- function(adj, var_of_eq) {
  n <- length(adj)
  eq_of_var <- integer(max(c(0L, unlist(adj), var_of_eq)))
  eq_of_var[var_of_eq[var_of_eq > 0L]] <- which(var_of_eq > 0L)
  deps <- lapply(seq_len(n), function(i) setdiff(eq_of_var[adj[[i]]], c(0L, i)))
  index <- integer(n); low <- integer(n); on <- logical(n)
  stack <- integer(); counter <- 0L; comps <- list()
  # iterative Tarjan to stay clear of the R recursion limit on long chains
  for (root in seq_len(n)) {
    if (index[root]) next
    work <- list(c(root, 0L))
    while (length(work)) {
      fr <- work[[length(work)]]
      v <- fr[1L]; pos <- fr[2L]
      if (pos == 0L) {
        counter <- counter + 1L
        index[v] <- counter; low[v] <- counter
        stack <- c(stack, v); on[v] <- TRUE
      }
      advanced <- FALSE
      d <- deps[[v]]
      while (pos < length(d)) {
        pos <- pos + 1L
        w <- d[pos]
        if (!index[w]) {
          work[[length(work)]] <- c(v, pos)
          work[[length(work) + 1L]] <- c(w, 0L)
          advanced <- TRUE
          break
        } else if (on[w]) low[v] <- min(low[v], index[w])
      }
      if (advanced) next
      if (low[v] == index[v]) {
        comp <- integer()
        repeat {
          w <- stack[length(stack)]
          stack <- stack[-length(stack)]
          on[w] <- FALSE
          comp <- c(comp, w)
          if (w == v) break
        }
        comps[[length(comps) + 1L]] <- sort(comp)
      }
      work[[length(work)]] <- NULL
      if (length(work)) {
        parent <- work[[length(work)]][1L]
        low[parent] <- min(low[parent], low[v])
      }
    }
  }
  comps
}

# The blocks of a square system in solving order (SCCDecomposition of
# sccnonlinearproblem.jl): the equation and unknown indices of each strongly
# connected component, and how each block is solved: "linear" when its
# equations are affine in its unknowns (a LinearProblem in ModelingToolkit),
# "homotopy" when they contain homotopy() nodes (a HomotopyProblem), and
# "nonlinear" otherwise. NULL when the system is structurally singular.
.nl_scc_blocks <- function(sys, resid = .ode_rhs_exprs(sys)) {
  us <- unknowns(sys)
  adj <- .nl_incidence(resid, us)
  m <- .nl_matching(adj, length(us))
  if (any(m$var_of_eq == 0L)) return(NULL)
  comps <- .nl_sccs(adj, m$var_of_eq)
  var_sccs <- lapply(comps, function(cp) m$var_of_eq[cp])
  hold <- .nl_hold(sys)
  kinds <- vapply(seq_along(comps), function(k) {
    ei <- comps[[k]]; vi <- var_sccs[[k]]
    if (.nl_exprs_have_call(resid[ei], "homotopy")) return("homotopy")
    vk <- .keys(us[vi])
    for (e in .partials(resid[ei], us[vi], c(hold, us[-vi])))
      if (!.is_num(e) && any(.nl_atomic_keys(e) %in% vk)) return("nonlinear")
    "linear"
  }, "")
  list(eq_sccs = comps, var_sccs = var_sccs, kinds = kinds, matching = m)
}

#' Nonlinear problems solved block by block
#'
#' `SCCNonlinearProblem(sys, op)` decomposes a square nonlinear system into
#' its strongly connected components: equations are matched to unknowns and
#' the blocks are ordered so that each depends only on earlier ones.
#' `solve()` then solves the blocks in sequence, which is faster and more
#' robust than one large Newton iteration. The fields `var_sccs` and
#' `eq_sccs` list the unknown and equation indices of each block, and
#' `scc_kinds` how each is solved: blocks affine in their unknowns by a
#' linear solve, blocks with [homotopy()] nodes by continuation, and the
#' others by the nonlinear algorithm given to `solve()`.
#'
#' @param sys A complete time independent `System` with as many equations
#'   as unknowns.
#' @param op Operating point: values of unknowns and parameters.
#' @param ... Solver options stored with the problem.
#' @param jac Generate the symbolic Jacobian.
#' @param guesses,missing_guess_value As in [NonlinearProblem()].
#' @return An `SCCNonlinearProblem`.
#' @examples
#' x <- Variable("x"); y <- Variable("y"); z <- Variable("z")
#' sys <- complete(System(list(0 ~ x^3 - 8, 0 ~ y - 2 * x, 0 ~ z^2 + z - y), name = "chain"))
#' prob <- SCCNonlinearProblem(sys, list(x ~ 1, y ~ 1, z ~ 1))
#' length(prob$var_sccs)
#' prob$scc_kinds
#' solve(prob)$u
#' @export
SCCNonlinearProblem <- function(sys, op = list(), ..., jac = TRUE, guesses = list(),
                                missing_guess_value = NULL) {
  prob <- .nl_build_problem("NonlinearProblem", sys, op, list(...), jac, NULL, NULL, guesses,
                            missing_guess_value, TRUE)
  .nl_as_scc(prob)
}

# Turn a square nonlinear problem built from a system into its block form.
.nl_as_scc <- function(prob) {
  sys <- prob$sys
  blocks <- .nl_scc_blocks(sys)
  if (is.null(blocks)) {
    m <- .nl_matching(.nl_incidence(.ode_rhs_exprs(sys), unknowns(sys)), length(unknowns(sys)))
    stop("the system is structurally singular: no equation can be matched to ",
         paste(names(prob$u0)[m$eq_of_var == 0L], collapse = ", "), call. = FALSE)
  }
  prob$eq_sccs <- blocks$eq_sccs
  prob$var_sccs <- blocks$var_sccs
  prob$scc_kinds <- blocks$kinds
  if (any(blocks$kinds == "homotopy")) prob$f$homotopy <- .nl_homotopy_residual(sys)
  class(prob) <- c("SCCNonlinearProblem", "AbstractNonlinearProblem", "SciMLProblem")
  prob
}

# ---------------------------------------------------------------------------
# Solutions

# A nonlinear, linear or optimization solution: one saved state. Built on
# .new_solution() so sol[x], sol[expr] and sol$ps work as for ODE solutions.
.nl_new_solution <- function(cls, prob, alg, u, retcode, resid = NULL, stats = list(), extra = list(),
                             t = 0) {
  u <- as.numeric(u)
  U <- matrix(u, nrow = 1L)
  nms <- names(prob$u0)
  if (!is.null(nms) && length(nms) == length(u)) colnames(U) <- nms
  .new_solution(cls, c(list(t = t, U = U, prob = prob, alg = alg, retcode = retcode, resid = resid,
                            stats = stats, p_final = prob$p, sys = prob$sys), extra))
}

.nl_sol_u <- function(x) {
  env <- .sol_env(x)
  u <- env$U[1L, ]
  if (!is.null(colnames(env$U))) names(u) <- colnames(env$U)
  # u has the shape of u0, as in Julia: a matrix for a matrix start (a point
  # of SymmetricPositiveDefinite)
  d <- dim(env$prob$u0)
  if (length(d) == 2L && prod(d) == length(u)) u <- matrix(unname(u), d[1L], d[2L])
  u
}

#' @export
`$.NonlinearSolution` <- function(x, name) {
  if (name == "u") return(.nl_sol_u(x))
  `$.SciMLSolution`(x, name)
}

#' @export
names.NonlinearSolution <- function(x) c("u", "resid", "retcode", "prob", "alg", "stats", "ps", "original")

#' @export
print.NonlinearSolution <- function(x, ...) {
  env <- .sol_env(x)
  cat(class(x)[1L], " with retcode ", env$retcode, if (!is.null(env$alg$name)) paste0(" (", env$alg$name, ")"),
      "\n", sep = "")
  cat("u:\n")
  print(.nl_sol_u(x))
  if (!is.null(env$objective)) cat("Objective: ", paste(format(env$objective), collapse = " "), "\n", sep = "")
  else if (length(env$resid)) cat("Residual norm: ", format(sqrt(sum(env$resid^2))), "\n", sep = "")
  invisible(x)
}

#' Whether a solve succeeded
#'
#' `successful_retcode(sol)` is `TRUE` when the return code of a solution
#' reports success, as `SciMLBase.successful_retcode`.
#'
#' @param sol A solution, or a return code string.
#' @return A logical.
#' @examples
#' sol <- solve(NonlinearProblem(function(u, p) u^2 - 2, 1))
#' successful_retcode(sol)
#' @export
successful_retcode <- function(sol) {
  rc <- if (is.character(sol)) sol else sol$retcode
  rc %in% c("Success", "Terminated", "ExactSolutionLeft", "ExactSolutionRight", "FloatingPointLimit",
            "StalledSuccess")
}
