# Discrete time: clocks and the operators Shift, ShiftIndex, Sample, Hold
# and SampleTime, after ModelingToolkitBase/src/discretedomain.jl and
# SciMLBase's TimeDomain clocks (Clock, SolverStepClock, ContinuousClock);
# shift2term, distribute_shift and simplify_shifts of
# ModelingToolkitBase/src/variables.jl; the shift branch of __mtkcompile in
# ModelingToolkitBase/src/systems/systems.jl; and DiscreteFunction and
# DiscreteProblem of ModelingToolkitBase/src/problems/discreteproblem.jl.
# ImplicitDiscreteProblem and IDSolve are in R/solvers-discrete.R.
#
# Operators are registered symcompR function nodes (see .op_apply):
#   Shift(x, steps, clock)   x shifted by `steps` ticks of clock id `clock`
#   Sample(x, clock)         x sampled by a clock (0: inferred)
#   Hold(x)                  zero order hold of a discrete x
#   SampleTime(clock)        the sample interval of a clock
# Clock ids index the table .mtk$clocks; 0 means inferred.

.mtk$clocks <- list()

# ---------------------------------------------------------------------------
# Clocks

#' Clocks and time domains
#'
#' `Clock(dt, phase)` is a periodic clock ticking at `phase + k dt`;
#' `Clock()` a discrete clock whose period is inferred. `SolverStepClock()`
#' ticks at every step of the integrator, and `ContinuousClock()` is the
#' continuous time domain. `TimeDomain(x)` returns the clock of a variable
#' indexed with a [ShiftIndex()], or `NULL`. These are the clocks of
#' SciMLBase.
#'
#' @param dt Sample interval, or `NULL` to infer it.
#' @param phase Time offset of the ticks.
#' @param x A variable.
#' @return A `TimeDomain` object.
#' @examples
#' Clock(0.1)
#' Clock(0.1, phase = 0.05)
#' SolverStepClock()
#' @export
Clock <- function(dt = NULL, phase = 0) {
  if (!is.null(dt) && (!is.numeric(dt) || length(dt) != 1L || dt <= 0))
    stop("the sample interval of a Clock must be a positive number", call. = FALSE)
  structure(list(type = "periodic", dt = if (is.null(dt)) NULL else as.numeric(dt),
                 phase = as.numeric(phase)), class = c("Clock", "TimeDomain"))
}

#' @rdname Clock
#' @export
SolverStepClock <- function() structure(list(type = "solverstep"), class = c("SolverStepClock", "TimeDomain"))

#' @rdname Clock
#' @export
ContinuousClock <- function() structure(list(type = "continuous"), class = c("ContinuousClock", "TimeDomain"))

#' @rdname Clock
#' @export
TimeDomain <- function(x) {
  m <- attr(as_sym(x), "mtk_meta")$timedomain
  if (is.null(m)) m <- .meta_get(x, "timedomain")
  m
}

.integer_sequence <- function() structure(list(type = "integer"), class = c("IntegerSequence", "TimeDomain"))

#' @export
format.TimeDomain <- function(x, ...) {
  switch(x$type,
         periodic = if (is.null(x$dt)) "Clock()" else
           paste0("Clock(", x$dt, if (x$phase != 0) paste0(", phase = ", x$phase), ")"),
         solverstep = "SolverStepClock()", continuous = "ContinuousClock()",
         integer = "IntegerSequence()", x$type)
}

#' @export
print.TimeDomain <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

# Clock id of a clock object, registering it on first use; 0 for none.
.clock_id <- function(clk) {
  if (is.null(clk)) return(0L)
  for (i in seq_along(.mtk$clocks)) if (identical(unclass(.mtk$clocks[[i]]), unclass(clk))) return(i)
  .mtk$clocks[[length(.mtk$clocks) + 1L]] <- clk
  length(.mtk$clocks)
}

.clock_of <- function(id) if (is.null(id) || id == 0L) NULL else .mtk$clocks[[id]]

.as_clock <- function(c) {
  if (is.null(c) || inherits(c, "TimeDomain")) return(c)
  if (is.numeric(c)) return(Clock(c))
  stop("cannot interpret ", class(c)[1L], " as a clock", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Operator nodes

.op_node <- function(op, args) {
  if (!exists(op, envir = .mtk$vars, inherits = FALSE)) .register_var(op, "operator")
  do.call(Function(op), lapply(args, as_sym))
}

.shift_parts <- function(e) {
  a <- .args(e)
  list(var = a[[1L]], steps = as.integer(as.numeric(a[[2L]])),
       clock = if (length(a) > 2L) as.integer(as.numeric(a[[3L]])) else 0L)
}

.make_shift <- function(x, steps, clock = 0L) {
  if (steps == 0L) return(x)
  args <- list(x, as.integer(steps))
  if (clock > 0L) args[[3L]] <- as.integer(clock)
  .op_node("Shift", args)
}

# Apply a shift to every time dependent variable of e (distribute_shift),
# merging nested shifts (simplify_shifts).
.apply_shift <- function(e, iv, steps, clock = 0L) {
  e <- as_sym(e)
  if (.is_op(e, "Shift")) {
    sp <- .shift_parts(e)
    return(.make_shift(sp$var, sp$steps + steps, if (clock > 0L) clock else sp$clock))
  }
  if (.is_varcall(e)) {
    if (!any(vapply(.args(e), function(a) identical(.key(a), .key(iv)), TRUE))) return(e)
    return(.make_shift(e, steps, clock))
  }
  if (!length(.args(e)) || .is_op(e)) return(e)
  if (.type(e) == "Derivative") return(.make_shift(e, steps, clock))
  .rebuild(e, lapply(.args(e), .apply_shift, iv = iv, steps = steps, clock = clock))
}

#' Discrete time operators
#'
#' `Shift(t, steps)` returns an operator shifting a discrete time variable by
#' `steps` ticks: `Shift(t, -1)(x)` is the previous value of `x`. The
#' operator distributes over expressions and merges with other shifts. Only
#' non-positive shifts are allowed in models.
#'
#' `k <- ShiftIndex(clock)` indexes variables on a clock: `k(x)` is `x` at
#' the current tick and `(k - 1)(x)` (or `k(x, -1)`) its value one tick
#' earlier. `ShiftIndex(t)` counts integer steps of a purely discrete
#' system, `ShiftIndex(dt)` uses `Clock(dt)`, and `ShiftIndex()` infers the
#' clock.
#'
#' `Sample(dt)(x)` (or `Sample(clock)(x)`) samples a continuous variable at
#' the ticks of a clock, `Sample(x)` with an inferred clock. `Hold(xd)` holds
#' a discrete variable constant between ticks (zero order hold), making it
#' available to continuous equations. `SampleTime()` is the sample interval
#' of the clock inferred for its equation.
#'
#' `hasshift()`, `hassample()` and `hashold()` test for the operators.
#'
#' @param t The independent variable.
#' @param steps Number of ticks.
#' @param clock A [Clock()], a sample interval, the independent variable (for
#'   integer steps) or `NULL` to infer it.
#' @param x A variable or expression.
#' @param e An expression or equation.
#' @return An operator (a function), or a `sym` expression.
#' @examples
#' t <- t_nounits
#' x <- Variable("x", t); u <- Variable("u", t)
#' Shift(t, -1)(x)
#' k <- ShiftIndex(Clock(0.1))
#' (k - 1)(x)
#' Hold(u)
#' Sample(0.1)(x)
#' @export
Shift <- function(t, steps = 1L) {
  t <- as_sym(t)
  steps <- as.integer(steps)
  fn <- function(x, clock = 0L) {
    if (inherits(x, "Equation")) return(Equation(fn(x$lhs, clock), fn(x$rhs, clock)))
    if (inherits(x, "formula")) return(fn(as_equation(x), clock))
    if (is.numeric(x) && !is_sym(x)) return(as_sym(x))
    if (steps == 0L) return(as_sym(x))
    .apply_shift(x, t, steps, clock)
  }
  structure(fn, class = c("Shift", "function"), iv = t, steps = steps)
}

#' @export
print.Shift <- function(x, ...) {
  cat("Shift(", format(attr(x, "iv")), ", ", attr(x, "steps"), ")\n", sep = "")
  invisible(x)
}

#' @export
`^.Shift` <- function(e1, e2) Shift(attr(e1, "iv"), attr(e1, "steps") * as.integer(e2))

#' @rdname Shift
#' @export
ShiftIndex <- function(clock = NULL, steps = 0L) {
  clk <- if (is_sym(clock)) .integer_sequence() else .as_clock(clock)
  steps <- as.integer(steps)
  fn <- function(x, shift = 0L) {
    x <- as_sym(x)
    vs <- Filter(function(v) !.is_derivative(v) && !isparameter(v) && !identical(.var_kind(v), "iv"),
                 .collect_vars(x))
    if (length(vs) != 1L)
      stop("Cannot shift a multivariate expression ", format(x), ". Either create a new unknown ",
           "and shift this, or shift the individual variables in the expression.", call. = FALSE)
    v <- vs[[1L]]
    if (.type(v) != "Function") stop("Cannot shift time-independent variable ", format(v), call. = FALSE)
    if (length(.args(v)) != 1L)
      stop("Cannot shift an expression with multiple independent variables ", format(x), call. = FALSE)
    n <- steps + as.integer(shift)
    if (n == 0L) {
      m <- attr(x, "mtk_meta") %||% .meta(x)
      m$timedomain <- clk
      attr(x, "mtk_meta") <- m
      return(x)
    }
    .apply_shift(x, .args(v)[[1L]], n, .clock_id(clk))
  }
  structure(fn, class = c("ShiftIndex", "function"), clock = clk, steps = steps)
}

#' @export
Ops.ShiftIndex <- function(e1, e2) {
  if (!(.Generic %in% c("+", "-")) || !inherits(e1, "ShiftIndex") || !is.numeric(e2))
    stop("a ShiftIndex supports k + i and k - i only", call. = FALSE)
  i <- as.integer(e2)
  ShiftIndex(attr(e1, "clock"), attr(e1, "steps") + if (.Generic == "+") i else -i)
}

#' @export
print.ShiftIndex <- function(x, ...) {
  cat("ShiftIndex(", format(attr(x, "clock") %||% Clock()), ", ", attr(x, "steps"), ")\n", sep = "")
  invisible(x)
}

#' @rdname Shift
#' @export
Sample <- function(clock = NULL) {
  if (is_sym(clock)) return(.op_node("Sample", list(clock)))
  cid <- .clock_id(.as_clock(clock))
  fn <- function(x) {
    x <- as_sym(x)
    .op_node("Sample", if (cid > 0L) list(x, cid) else list(x))
  }
  structure(fn, class = c("Sample", "function"), clock = cid)
}

#' @rdname Shift
#' @export
Hold <- function(x) {
  if (missing(x)) return(structure(function(x) Hold(x), class = c("Hold", "function")))
  if (is.numeric(x) && !is_sym(x)) return(as_sym(x))
  x <- as_sym(x)
  if (.is_num(x)) return(x)
  .op_node("Hold", list(x))
}

#' @rdname Shift
#' @export
SampleTime <- function(clock = NULL) .op_node("SampleTime", list(.clock_id(.as_clock(clock))))

.has_op <- function(e, op) {
  if (inherits(e, "Equation")) return(.has_op(e$lhs, op) || .has_op(e$rhs, op))
  if (is.list(e) && !is_sym(e)) return(any(vapply(e, .has_op, TRUE, op = op)))
  if (!is_sym(e)) return(FALSE)
  if (.is_op(e, op)) return(TRUE)
  for (a in .args(e)) if (.has_op(a, op)) return(TRUE)
  FALSE
}

#' @rdname Shift
#' @export
hasshift <- function(e) .has_op(e, "Shift")
#' @rdname Shift
#' @export
hassample <- function(e) .has_op(e, "Sample")
#' @rdname Shift
#' @export
hashold <- function(e) .has_op(e, "Hold")

# ---------------------------------------------------------------------------
# Lowered names of shifted variables

#' Lowered names of shifted variables
#'
#' `shift2term(Shift(t, -k)(x))` returns the variable \ifelse{latex}{\verb{x_(t-k)(t)}}{\verb{xₜ₋ₖ(t)}} that
#' [mtkcompile()] uses for the value of `x` k ticks earlier; a shift of a
#' lowered variable that cancels returns the original variable.
#'
#' @param var A shifted variable.
#' @return A `sym` expression.
#' @examples
#' t <- t_nounits
#' x <- Variable("x", t)
#' shift2term(Shift(t, -1)(x))
#' @export
shift2term <- function(var) {
  var <- as_sym(var)
  if (!.is_op(var, "Shift")) return(var)
  sp <- .shift_parts(var)
  arg <- sp$var
  root <- .meta_get(arg, "unshifted")
  base <- if (is.null(root)) arg else root
  back <- sp$steps + as.integer(.meta_get(arg, "shift", 0L))
  if (back == 0L) return(base)
  digits <- strsplit(as.character(abs(back)), "")[[1L]]
  sub <- paste(intToUtf8(0x2080 + as.integer(digits), multiple = TRUE), collapse = "")
  newname <- paste0(.name(base), intToUtf8(0x209c), intToUtf8(if (back > 0L) 0x208a else 0x208b), sub)
  meta <- .meta(base)
  meta$kind <- NULL
  meta$default <- NULL
  meta$guess <- NULL
  meta$timedomain <- NULL
  meta$unshifted <- base
  meta$shift <- back
  .make_var(newname, "unknown", .args(base), meta)
}

# Shift nodes of a list of expressions, unique by key.
.shift_nodes <- function(exprs) .find_ops(exprs, "Shift")

.eq_sides <- function(eqs) unlist(lapply(eqs, function(e) list(e$lhs, e$rhs)), recursive = FALSE)

# ---------------------------------------------------------------------------
# The shift branch of __mtkcompile

.compile_shift <- function(sys, eqs, all_dvs, obs0, new_ps, inputs, outputs, binds, ics) {
  if (!any(vapply(eqs, hasshift, TRUE))) return(NULL)
  iv <- .f(sys, "iv")
  low <- .lower_shifts(eqs, all_dvs, iv)
  if (any(vapply(low$eqs, function(e) .has_derivative(e$lhs) || .has_derivative(e$rhs), TRUE)))
    stop("ModelingToolkitBase.jl cannot simplify systems with both `Shift` and `Differential` ",
         "operators.", call. = FALSE)
  # the defining equations of the current values become observed, as far as
  # they can be ordered; the rest are algebraic
  ord <- .partial_toposort(c(obs0, low$defs))
  algeqs <- c(low$alg, ord$rest)
  tr <- .trivial_observed(low$diffvars, algeqs, ord$sorted, all_dvs, iv)
  algeqs <- lapply(tr$algeqs, function(e) Equation(.S0(), e$rhs - e$lhs))
  obseqs <- .toposort_observed(tr$obseqs)
  obs_keys <- vapply(obseqs, function(e) .key(e$lhs), "")
  io_keys <- .keys(inputs)
  alg_vars <- Filter(function(v) !(.key(v) %in% c(.keys(low$diffvars), obs_keys, io_keys)), all_dvs)
  new_eqs <- c(low$diffeqs, algeqs)
  for (e in c(new_eqs, obseqs)) {
    if (hasshift(e$rhs))
      stop("ModelingToolkitBase.jl is unable to simplify such systems. Encountered shift in RHS of ",
           "equation ", .fmt(e), call. = FALSE)
  }
  .update_system(sys, eqs = new_eqs, observed = obseqs, unknowns = c(low$diffvars, alg_vars),
                 ps = new_ps, inputs = inputs, outputs = outputs, isscheduled = TRUE,
                 is_discrete = TRUE, bindings = binds, initial_conditions = ics,
                 schedule = list(var_sccs = list(seq_along(new_eqs)), dummy_sub = low$sub))
}

# Lowest shift per variable, history variables x_{t-i} with their update
# equations Shift(x_{t-i}, 1) ~ x_{t-i+1}, and the equations with every
# shifted variable replaced by its history variable.
.lower_shifts <- function(eqs, all_dvs, iv) {
  lowest <- list()
  vars <- list()
  for (n in .shift_nodes(.eq_sides(eqs))) {
    sp <- .shift_parts(n)
    if (sp$steps > 0L)
      stop("Positive shifts are disallowed in unsimplified equations. Found ", format(n), ".",
           call. = FALSE)
    k <- .key(sp$var)
    lowest[[k]] <- min(lowest[[k]] %||% 0L, sp$steps)
    vars[[k]] <- sp$var
  }
  dv_keys <- .keys(all_dvs)
  is_def <- vapply(eqs, function(e) .key(e$lhs) %in% dv_keys, TRUE)
  diffeqs <- list(); diffvars <- list(); from <- list(); to <- list()
  hist <- list()
  for (k in names(lowest)) {
    var <- vars[[k]]
    order <- -lowest[[k]]
    hv <- list()
    for (i in seq_len(order)) {
      h <- shift2term(.make_shift(var, -i))
      prev <- if (i == 1L) var else hv[[i - 1L]]
      diffeqs[[length(diffeqs) + 1L]] <- Equation(.make_shift(h, 1L), prev)
      diffvars[[length(diffvars) + 1L]] <- h
      hv[[i]] <- h
    }
    hist[[k]] <- list(var = var, vars = hv)
  }
  # substitute every shift node, whatever its clock argument
  for (n in .shift_nodes(.eq_sides(eqs))) {
    sp <- .shift_parts(n)
    from[[length(from) + 1L]] <- n
    to[[length(to) + 1L]] <- shift2term(.make_shift(sp$var, sp$steps))
  }
  rule <- .rule(from, to)
  subeq <- function(e) Equation(.subst(e$lhs, rule), .subst(e$rhs, rule))
  list(diffeqs = diffeqs, diffvars = diffvars, hist = hist, sub = rule,
       defs = lapply(eqs[is_def], subeq), alg = lapply(eqs[!is_def], subeq),
       eqs = lapply(eqs, subeq))
}

# Topological order of explicit equations lhs ~ rhs; equations in cycles are
# returned in `rest` (topsort_equations with check = false).
.partial_toposort <- function(eqs) {
  if (!length(eqs)) return(list(sorted = list(), rest = list()))
  lhs <- vapply(eqs, function(o) .key(o$lhs), "")
  deps <- lapply(eqs, function(o) intersect(.keys(.collect_vars(o$rhs)), lhs))
  done <- character()
  out <- list()
  remaining <- seq_along(eqs)
  repeat {
    ready <- remaining[vapply(remaining, function(i) all(deps[[i]] %in% done) &&
                                !(lhs[i] %in% deps[[i]]), TRUE)]
    if (!length(ready)) break
    for (i in ready) {
      out[[length(out) + 1L]] <- eqs[[i]]
      done <- c(done, lhs[i])
    }
    remaining <- setdiff(remaining, ready)
  }
  list(sorted = out, rest = eqs[remaining])
}

# ---------------------------------------------------------------------------
# Discrete problems

#' Discrete functions and problems
#'
#' `DiscreteProblem(sys, op, tspan)` builds the map `u[k+1] = f(u[k], p,
#' t[k])` of a compiled discrete system, whose unknowns are the past values
#' \ifelse{latex}{\verb{x_(t-1)}}{\verb{xₜ₋₁}}, \ifelse{latex}{\verb{x_(t-2)}}{\verb{xₜ₋₂}}, ... of its variables (see [shift2term()]). Initial
#' conditions give the past: `(k - 1)(x) ~ 1` (or `Shift(t, -1)(x) ~ 1`)
#' sets \ifelse{latex}{\verb{x_(t-1)}}{\verb{xₜ₋₁}} at the initial time. A value (or default) for `x` itself is
#' its value at the initial time, and the past values not given are solved
#' for to match it, from `guesses`, as ModelingToolkit's initialization
#' does. `DiscreteFunction(sys)` builds `f`. Solve with `solve(prob,
#' FunctionMap())`, which steps by `dt` (default 1) and saves every step.
#'
#' `DiscreteProblem(f, u0, tspan, p)` builds a problem from an R function
#' `f(u, p, t)` returning the next state.
#'
#' @param sys A compiled discrete `System`, or an R function `f(u, p, t)`.
#' @param op Values of past unknowns and of parameters.
#' @param tspan Time span `c(t0, t1)`.
#' @param u0,p Initial state and parameters for problems from functions.
#' @param guesses Guesses for past values solved from initial conditions.
#' @param ... Solver options stored with the problem.
#' @return A `DiscreteProblem` or `DiscreteFunction`.
#' @examples
#' t <- t_nounits
#' x <- Variable("x", t)
#' k <- ShiftIndex(t)
#' sys <- mtkcompile(System(list(x ~ 0.5 * (k - 1)(x)), t, name = "halving"))
#' prob <- DiscreteProblem(sys, list((k - 1)(x) ~ 1), c(0, 5))
#' sol <- solve(prob, FunctionMap())
#' sol[x]
#' @export
DiscreteFunction <- function(sys, ...) {
  if (is.function(sys))
    return(structure(list(f = sys, sys = NULL), class = c("DiscreteFunction", "SciMLFunction")))
  .need_complete(sys, "DiscreteFunction")
  if (is.null(.f(sys, "iv"))) stop("a DiscreteFunction needs a time dependent system", call. = FALSE)
  eqs <- equations(sys)
  if (!all(vapply(eqs, function(e) .is_op(e$lhs, "Shift"), TRUE)))
    stop("the system has algebraic equations; it is an implicit discrete system", call. = FALSE)
  f <- .build_closure(lapply(eqs, `[[`, "rhs"), .codegen_args(sys), observed = .codegen_observed(sys),
                      out = "vector")
  structure(list(f = f, sys = sys,
                 observed = function(expr) build_explicit_observed_function(sys, expr)),
            class = c("DiscreteFunction", "SciMLFunction"))
}

#' @rdname DiscreteFunction
#' @export
DiscreteProblem <- function(sys, op = list(), tspan = NULL, ..., u0 = NULL, p = NULL, guesses = list()) {
  dots <- list(...)
  # DiscreteProblem(u0, tspan, p): the identity map, the inner problem of
  # pure jump problems (DISCRETE_OUTOFPLACE_DEFAULT of SciMLBase)
  if (is.numeric(sys)) {
    p <- p %||% if (!missing(tspan) && !is.null(tspan)) tspan
    return(.new_problem("DiscreteProblem", f = DiscreteFunction(function(u, p, t) u), u0 = sys,
                        tspan = .as_tspan(op), p = if (is.null(p)) numeric() else p, sys = NULL, kwargs = dots))
  }
  if (is.function(sys) || inherits(sys, "DiscreteFunction")) {
    f <- if (inherits(sys, "DiscreteFunction")) sys else DiscreteFunction(sys)
    if (is.null(u0)) u0 <- op
    # DiscreteProblem(f, u0, tspan, p)
    if (is.null(p) && length(dots) && (is.null(names(dots)) || !nzchar(names(dots)[1L]))) {
      p <- dots[[1L]]
      dots <- dots[-1L]
    }
    return(.new_problem("DiscreteProblem", f = f, u0 = .as_u0(u0), tspan = .as_tspan(tspan),
                        p = if (is.null(p)) numeric() else p, sys = NULL, kwargs = dots))
  }
  .need_complete(sys, "DiscreteProblem")
  tspan <- .as_tspan(tspan, sys)
  fn <- DiscreteFunction(sys)
  vals <- .process_op(sys, .discrete_op(sys, op))
  u0 <- .solve_discrete_history(sys, vals, tspan[1L], guesses)
  prob <- .new_problem("DiscreteProblem", f = fn, u0 = u0, tspan = tspan, p = vals$p,
                       sys = sys, kwargs = .problem_kwargs(sys, dots))
  prob$op <- vals$op
  prob
}

# Operating point for a discrete system: values given for Shift(x, -i) refer
# to the history variable x_{t-i} (add_toterms!).
.discrete_op <- function(sys, op) {
  m <- as_symmap(op)
  out <- symmap()
  for (k in names(m)) {
    el <- m[[k]]
    if (!is.null(el$var) && .is_op(el$var, "Shift")) {
      out <- smap_set(out, shift2term(el$var), el$value)
    } else {
      out[[k]] <- el
    }
  }
  class(out) <- "mtk_symmap"
  out
}

# Missing history values are found, as by MTK's initialization, from the
# values given for the current variables (op or defaults): x(t0) = value is
# solved for the unknown past values with Newton's method, starting from
# `guesses` (or the given value of the variable).
.solve_discrete_history <- function(sys, vals, t0, guesses = list()) {
  u0 <- vals$u0
  miss <- which(is.na(u0))
  if (!length(miss)) return(u0)
  us <- unknowns(sys)
  roots <- .unique_syms(Filter(Negate(is.null), lapply(us[miss], .meta_get, field = "unshifted")))
  target <- vapply(roots, function(r) {
    v <- vals$values[[.key(r)]]
    if (is.numeric(v) && length(v) == 1L && is.finite(v)) v else NA_real_
  }, 0)
  roots <- roots[!is.na(target)]
  target <- target[!is.na(target)]
  if (length(roots) < length(miss))
    stop("Initial condition for ", paste(names(u0)[miss], collapse = ", "), " not provided.", call. = FALSE)
  g <- smap_merge(guesses(sys), .discrete_op(sys, .resolve_named(as_symmap(guesses), c(us, roots))))
  for (i in miss) {
    gv <- smap_get(g, us[[i]])
    r <- match(.key(.meta_get(us[[i]], "unshifted")), .keys(roots))
    u0[i] <- if (is.numeric(gv)) gv else if (!is.na(r)) target[r] else 1
  }
  obs_fn <- .build_closure(roots, .codegen_args(sys), observed = .codegen_observed(sys), out = "vector")
  res <- function(z) {
    u <- u0
    u[miss] <- z
    obs_fn(u, vals$p, t0) - target
  }
  sol <- .newton_solve(res, u0[miss])
  if (!sol$converged)
    stop("could not find past values ", paste(names(u0)[miss], collapse = ", "),
         " matching the initial conditions; give them with ShiftIndex, as (k - 1)(x) ~ value",
         call. = FALSE)
  u0[miss] <- sol$x
  u0
}

#' @rdname DiscreteFunction
#' @export
FunctionMap <- function() structure(list(name = "FunctionMap"), class = c("FunctionMap", "DiscreteAlgorithm"))

#' @export
solve.DiscreteProblem <- function(a, b, ...) {
  if (!missing(b) && !is.null(b) && !inherits(b, "FunctionMap"))
    stop("a DiscreteProblem is solved with FunctionMap()", call. = FALSE)
  opts <- utils::modifyList(a$kwargs %||% list(), list(...))
  .function_map(a, opts)
}

# FunctionMap: u[k+1] = f(u[k], p, t[k]) with t[k] = t0 + k dt. Discrete
# callbacks run after each step (and may terminate); the saved state is the
# one after the affects, and parameter changes are recorded for sol[expr].
.function_map <- function(prob, opts) {
  t0 <- prob$tspan[1L]; tf <- prob$tspan[2L]
  dt <- abs(opts$dt %||% 1)
  n <- floor((tf - t0) / dt + 1e-9)
  ts <- t0 + (0:n) * dt
  if (abs(ts[length(ts)] - tf) > 1e-9 * max(1, abs(tf))) ts <- c(ts, tf)
  cbs <- .as_callback_set(opts$callback)
  if (length(cbs$continuous)) stop("continuous callbacks need a continuous time problem", call. = FALSE)
  integ <- new.env(parent = emptyenv())
  integ$t <- t0
  integ$u <- .state_vec(prob$u0)
  integ$p <- prob$p
  integ$prob <- prob
  integ$terminate <- FALSE
  integ$u_modified <- FALSE
  integ$sol_events <- list()
  ph_t <- t0
  ph_p <- list()
  .init_callbacks(cbs, integ)
  ph_p[[1L]] <- integ$p
  U <- matrix(NA_real_, length(ts), length(prob$u0))
  U[1L, ] <- integ$u
  f <- .shape_f(prob$f$f, dim(prob$u0))
  last <- 1L
  for (i in seq_len(length(ts) - 1L)) {
    integ$u[] <- .as_num(f(integ$u, integ$p, ts[i]))
    integ$t <- ts[i + 1L]
    for (cb in cbs$discrete) if (isTRUE(cb$condition(integ$u, integ$t, integ))) cb$affect(integ)
    if (!identical(integ$p, ph_p[[length(ph_p)]])) {
      ph_t <- c(ph_t, integ$t)
      ph_p[[length(ph_p) + 1L]] <- integ$p
    }
    U[i + 1L, ] <- integ$u
    last <- i + 1L
    if (integ$terminate) break
  }
  for (cb in cbs$discrete) if (!is.null(cb$finalize)) cb$finalize(integ)
  U <- U[seq_len(last), , drop = FALSE]
  colnames(U) <- names(prob$u0)
  sol <- .new_ode_solution(t = ts[seq_len(last)], u = U, prob = prob, alg = FunctionMap(),
                           retcode = if (integ$terminate) "Terminated" else "Success",
                           stats = list(nf = last - 1L), p = integ$p, events = integ$sol_events)
  if (length(dim(prob$u0)) > 1L) assign("u_dim", dim(prob$u0), envir = .sol_env(sol))
  if (length(ph_p) > 1L) {
    env <- .sol_env(sol)
    env$param_history <- list(t = ph_t, p = do.call(rbind, ph_p))
  }
  sol
}

# ---------------------------------------------------------------------------
# Hooks

.mtk$hooks$compile_shift <- function(...) .compile_shift(...)
.mtk$hooks$shift2term <- function(x) shift2term(x)
