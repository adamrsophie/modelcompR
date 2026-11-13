# Equations, differential operators and the operators Pre, Initial and
# EvalAt, after Symbolics' Equation and ModelingToolkitBase's operators.

#' Symbolic equations
#'
#' `Equation(lhs, rhs)` builds the equation `lhs ~ rhs`. Every function that
#' takes equations also accepts R formulas, so `D(x) ~ -k * x` can be written
#' directly in a `list()`; both sides of a formula are evaluated in the
#' environment where the formula was created.
#'
#' @param lhs,rhs Expressions or numbers.
#' @return An object of class `Equation`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); k <- Parameter("k")
#' Equation(D(x), -k * x)
#' as_equation(D(x) ~ -k * x)
#' @export
Equation <- function(lhs, rhs) {
  structure(list(lhs = .as_eq_side(lhs), rhs = .as_eq_side(rhs)), class = "Equation")
}

.as_eq_side <- function(x) {
  if (is_sym(x)) return(x)
  if (inherits(x, c("Connection", "AnalysisPoint", "mtk_array"))) return(x)
  if (is.list(x) && !is.null(x$connection)) return(x)
  as_sym(x)
}

#' @rdname Equation
#' @param x A formula, an `Equation`, or a list of them.
#' @export
as_equation <- function(x) {
  if (inherits(x, "Equation")) return(x)
  if (inherits(x, "Connection")) return(x)
  if (inherits(x, "formula")) {
    if (length(x) != 3L) stop("an equation formula needs both sides: lhs ~ rhs", call. = FALSE)
    env <- environment(x) %||% parent.frame()
    lhs <- eval(x[[2L]], env)
    rhs <- eval(x[[3L]], env)
    if (inherits(lhs, "mtk_array") || inherits(rhs, "mtk_array") ||
        (is.list(lhs) && !is_sym(lhs)) || (is.list(rhs) && !is_sym(rhs)))
      return(.array_equations(lhs, rhs))
    return(Equation(lhs, rhs))
  }
  stop("cannot interpret ", class(x)[1], " as an equation", call. = FALSE)
}

# Elementwise equations between arrays (lists) of expressions.
.array_equations <- function(lhs, rhs) {
  l <- if (is_sym(lhs) || is.numeric(lhs)) NULL else unclass(lhs)
  r <- if (is_sym(rhs) || is.numeric(rhs)) NULL else unclass(rhs)
  n <- max(length(l), length(r))
  structure(lapply(seq_len(n), function(i) {
    Equation(if (is.null(l)) lhs else l[[i]], if (is.null(r)) rhs[[min(i, length(rhs))]] else r[[i]])
  }), class = "mtk_eqlist")
}

# Normalize equation input into a flat list of Equation (and Connection)
# objects.
.as_equations <- function(eqs) {
  if (is.null(eqs)) return(list())
  if (inherits(eqs, c("Equation", "formula", "Connection"))) eqs <- list(eqs)
  out <- list()
  for (e in eqs) {
    if (inherits(e, "mtk_eqlist")) {
      out <- c(out, unclass(e))
      next
    }
    if (is.list(e) && !inherits(e, c("Equation", "Connection")) && !is_sym(e)) {
      out <- c(out, .as_equations(e))
      next
    }
    eq <- as_equation(e)
    # an equation between arrays is scalarized elementwise
    if (inherits(eq, "Equation") && (inherits(eq$lhs, "mtk_array") || inherits(eq$rhs, "mtk_array")))
      eq <- .array_equations(eq$lhs, eq$rhs)
    if (inherits(eq, "mtk_eqlist")) out <- c(out, unclass(eq)) else out[[length(out) + 1L]] <- eq
  }
  out
}

#' @export
print.Equation <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' @export
format.Equation <- function(x, ...) paste(.fmt_side(x$lhs), "~", .fmt_side(x$rhs))

.fmt_side <- function(s) if (is_sym(s)) format(s) else if (inherits(s, "Connection")) format(s) else "<object>"

#' @export
as.character.Equation <- function(x, ...) format(x)

#' @export
print.mtk_eqlist <- function(x, ...) {
  for (e in x) print(e)
  invisible(x)
}

#' @export
`==.Equation` <- function(e1, e2) {
  inherits(e2, "Equation") && identical(.key(e1$lhs), .key(e2$lhs)) &&
    identical(.key(e1$rhs), .key(e2$rhs))
}

#' Sides of an equation
#'
#' @param eq An equation or formula.
#' @return A `sym` expression.
#' @export
eq_lhs <- function(eq) as_equation(eq)$lhs
#' @rdname eq_lhs
#' @export
eq_rhs <- function(eq) as_equation(eq)$rhs

#' Inequality constraints
#'
#' `Inequality(lhs, rhs, "<=")` states `lhs <= rhs`, used as a constraint of an
#' optimization or optimal control problem. `lhs %<~% rhs` and `lhs %>~% rhs`
#' are shorthands for `<=` and `>=` (Julia's \ifelse{latex}{\verb{lesssim}}{\verb{≲}} and \ifelse{latex}{\verb{gtrsim}}{\verb{≳}}).
#'
#' @param lhs,rhs Expressions.
#' @param relation `"<="` or `">="`.
#' @return An object of class `Inequality`.
#' @export
Inequality <- function(lhs, rhs, relation = c("<=", ">=")) {
  relation <- match.arg(relation)
  structure(list(lhs = as_sym(lhs), rhs = as_sym(rhs), relation = relation),
            class = "Inequality")
}

#' @rdname Inequality
#' @param a,b Expressions.
#' @export
`%<~%` <- function(a, b) Inequality(a, b, "<=")
#' @rdname Inequality
#' @export
`%>~%` <- function(a, b) Inequality(a, b, ">=")

#' @export
format.Inequality <- function(x, ...) paste(format(x$lhs), if (x$relation == "<=") "\u2272" else "\u2273", format(x$rhs))
#' @export
print.Inequality <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

# Residual form lhs - rhs, with <= normalized to lhs - rhs <= 0.
.ineq_residual <- function(x) {
  if (x$relation == "<=") x$lhs - x$rhs else x$rhs - x$lhs
}

# ---------------------------------------------------------------------------
# Differential operators

# An expression of symcompR's differential geometry (diffgeom.R of symcompR):
# scalar and vector fields, forms and the operators on them.
.DIFFGEOM_TYPES <- c("BaseScalarField", "BaseVectorField", "Differential", "TensorProduct", "WedgeProduct",
                     "Commutator", "LieDerivative", "BaseCovarDerivativeOp", "CovarDerivativeOp")
.is_diffgeom <- function(e) {
  if (!is_sym(e)) return(FALSE)
  if (.type(e) %in% .DIFFGEOM_TYPES) return(TRUE)
  for (a in .args(e)) if (.is_diffgeom(a)) return(TRUE)
  FALSE
}

#' Differential operators
#'
#' `Differential(t)` returns the operator `D` with `D(x)` the derivative of
#' `x` with respect to `t`, kept unevaluated. Powers of the operator give
#' higher derivatives: `Differential(t, 2)` or `D(D(x))`.
#' `expand_derivatives()` evaluates derivatives of expressions, leaving
#' derivatives of unknown variables in place.
#'
#' When symcompR exports its own `Differential()` (the exterior derivative of
#' differential geometry), arguments other than a symbol are passed to it.
#'
#' @param iv The independent variable.
#' @param order Order of the derivative.
#' @param ... `form_field` of symcompR's `Differential()`.
#' @return A function of class `Differential`.
#' @examples
#' t <- t_nounits
#' D <- Differential(t)
#' x <- Variable("x", t)
#' D(x^2)
#' expand_derivatives(D(x^2))
#' @export
Differential <- function(iv, order = 1L, ...) {
  # symcompR's Differential() is the exterior derivative of differential
  # geometry, which takes a form field: a call with a differential geometry
  # object (a scalar field, a form, ...) goes to it as written, so attaching
  # both packages is safe.
  if ((missing(iv) || ...length() || .is_diffgeom(iv)) &&
      .symcompr_has("Differential"))
    return(.symcompr_forward("Differential", sys.call(), parent.frame()))
  if (...length()) stop("Differential(iv, order) takes no further arguments", call. = FALSE)
  iv <- as_sym(iv)
  order <- as.integer(order)
  fn <- function(e) {
    if (inherits(e, "mtk_array") || (is.list(e) && !is_sym(e)))
      return(structure(lapply(e, fn), class = class(e)))
    e <- as_sym(e)
    if (.is_num(e)) return(.S0())
    Derivative(e, iv, order)
  }
  structure(fn, class = c("Differential", "function"), iv = iv, order = order)
}

#' @export
print.Differential <- function(x, ...) {
  o <- attr(x, "order")
  cat("Differential(", format(attr(x, "iv")), ")", if (o > 1L) paste0("^", o), "\n", sep = "")
  invisible(x)
}

#' @export
`^.Differential` <- function(e1, e2) Differential(attr(e1, "iv"), attr(e1, "order") * as.integer(e2))

#' @export
Ops.Differential <- function(e1, e2) {
  if (.Generic == "*" && inherits(e1, "Differential") && inherits(e2, "Differential")) {
    if (!identical(.key(attr(e1, "iv")), .key(attr(e2, "iv"))))
      stop("products of differentials in different variables are not supported", call. = FALSE)
    return(Differential(attr(e1, "iv"), attr(e1, "order") + attr(e2, "order")))
  }
  if (.Generic == "^") return(Differential(attr(e1, "iv"), attr(e1, "order") * as.integer(e2)))
  stop("unsupported operation on a Differential", call. = FALSE)
}

#' @rdname Differential
#' @param e An expression, an equation or a list of them.
#' @export
expand_derivatives <- function(e) {
  if (inherits(e, "Equation")) return(Equation(expand_derivatives(e$lhs), expand_derivatives(e$rhs)))
  if (is.list(e) && !is_sym(e)) return(lapply(e, expand_derivatives))
  e <- as_sym(e)
  .expand_derivs(e)
}

# Evaluate Derivative nodes bottom up. A derivative of a declared variable
# stays unevaluated; anything else is differentiated with diff().
.expand_derivs <- function(e) {
  if (!length(.args(e))) return(e)
  ty <- .type(e)
  if (ty == "Derivative") {
    inner <- .expand_derivs(.args(e)[[1L]])
    out <- inner
    for (spec in .args(e)[-1L]) {
      v <- .args(spec)[[1L]]
      n <- as.integer(as.numeric(.args(spec)[[2L]]))
      out <- diff(out, v, n)
    }
    return(out)
  }
  if (!.has_derivative(e)) return(e)
  args <- lapply(.args(e), .expand_derivs)
  .rebuild(e, args)
}

.has_derivative <- function(e) {
  if (.type(e) == "Derivative") return(TRUE)
  for (a in .args(e)) if (.has_derivative(a)) return(TRUE)
  FALSE
}

# Rebuild a node with new arguments using only exported symcompR API: the
# replacement of each old argument by its new value, done in one xreplace.
.rebuild <- function(e, args) {
  old <- .args(e)
  if (identical(.keys(old), .keys(args))) return(e)
  ty <- .type(e)
  if (ty == "Add") return(Reduce(`+`, args))
  if (ty == "Mul") return(Reduce(`*`, args))
  if (ty == "Pow") return(args[[1L]]^args[[2L]])
  if (ty == "Function") return(do.call(Function(.name(e)), args))
  rule <- lapply(seq_along(old), function(i) list(old[[i]], args[[i]]))
  xreplace(e, rule)
}

# ---------------------------------------------------------------------------
# Operators used by events and initialization

#' Operators on variables
#'
#' `Pre(x)` is the value of `x` just before an event, used in event affects:
#' `x ~ Pre(x) + 1`. `Initial(x)` is the initial value of `x`, available as
#' a parameter of every compiled system. `EvalAt(tau)(x)` evaluates a time
#' dependent variable at the time `tau`, for boundary conditions and delays.
#'
#' @param x A variable or expression.
#' @return A `sym` expression.
#' @export
Pre <- function(x) {
  x <- as_sym(x)
  if (.is_num(x)) return(x)
  if (.is_varcall(x) || (.type(x) == "Symbol" && .is_declared(x))) return(.op_apply("Pre", x))
  .map_leaves(x, function(v) if (.is_varcall(v) || (.type(v) == "Symbol" && isparameter(v))) .op_apply("Pre", v) else v)
}

#' @rdname Pre
#' @export
Initial <- function(x) {
  x <- as_sym(x)
  if (.is_num(x)) return(x)
  .op_apply("Initial", x)
}

#' @rdname Pre
#' @export
isinitial <- function(x) is_sym(x) && .type(x) == "Function" && identical(.name(x), "Initial")

#' @rdname Pre
#' @param tau The time at which to evaluate.
#' @export
EvalAt <- function(tau) {
  tau <- as_sym(tau)
  function(x) {
    x <- as_sym(x)
    .map_leaves(x, function(v) {
      if (.is_varcall(v)) return(do.call(Function(.name(v)), list(tau)))
      # D(x) at tau is the lowered derivative variable xˍt at tau
      if (.is_derivative(v) && .is_varcall(.deriv_info(v)$var)) {
        tt <- default_toterm(v)
        if (.is_varcall(tt)) return(do.call(Function(.name(tt)), list(tau)))
      }
      v
    })
  }
}

# Operators are symcompR undefined functions whose names are registered as
# operator kinds, so they are never mistaken for unknowns.
.op_apply <- function(op, x) {
  if (!exists(op, envir = .mtk$vars, inherits = FALSE)) .register_var(op, "operator")
  do.call(Function(op), list(x))
}

.is_op <- function(e, op = NULL) {
  is_sym(e) && .type(e) == "Function" && identical(.var_kind(e), "operator") &&
    (is.null(op) || identical(.name(e), op))
}

# Apply f to every declared variable leaf (and derivative of one).
.map_leaves <- function(e, f) {
  vars <- .collect_vars(e)
  if (!length(vars)) return(e)
  to <- lapply(vars, f)
  .subst(e, .rule(vars, to))
}

#' Equations by kind
#'
#' `is_diff_equation()` tests whether an equation contains a derivative of an
#' unknown and `is_alg_equation()` the converse. `diff_equations()` and
#' `alg_equations()` filter the equations of a system.
#'
#' @param eq An equation.
#' @return A logical, or a list of equations.
#' @export
is_diff_equation <- function(eq) {
  eq <- as_equation(eq)
  .has_derivative(eq$lhs) || .has_derivative(eq$rhs)
}
#' @rdname is_diff_equation
#' @export
is_alg_equation <- function(eq) !is_diff_equation(eq)
