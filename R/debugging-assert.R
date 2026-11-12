# Debugging and assertions, after ModelingToolkitBase/src/debugging.jl
# (LoggedFun, logged_fun, debug_sub, _nan_condition, _debug_assertion,
# ASSERTION_LOG_VARIABLE, get_assertions_expr), debug_system in
# systems/abstractsystem.jl and the assertion term added to generated
# right-hand sides in systems/codegen.jl (generate_rhs).
#
# symcompR has no division or square root node: a / b is a * b^-1 and
# sqrt(a) is a^(1/2). debug_system() therefore maps "sqrt" to powers with
# exponent 1/2, "/" and "inv" to negative powers, and "^" to the remaining
# powers.

#' Parameter controlling assertion logging
#'
#' `ASSERTION_LOG_VARIABLE` is the parameter that [debug_system()] adds to a
#' system. When its value is nonzero (the default set by `debug_system()`),
#' the message of a failing assertion is shown each time the right-hand side
#' is evaluated; set it to 0 with [remake()] to silence the messages.
#'
#' @format A `sym` parameter.
#' @export
ASSERTION_LOG_VARIABLE <- local({
  nm <- "__log_assertions_\u2098\u209C\u2096"
  assign(nm, list(kind = "parameter", meta = list(default = FALSE, tunable = FALSE), ivs = list()),
         envir = .mtk$vars)
  structure(Symbol(nm), mtk_meta = list(kind = "parameter", default = FALSE, tunable = FALSE))
})

# _nan_condition: 0 when the condition holds, NaN otherwise.
.nan_condition <- register_symbolic("_nan_condition", function(cond) if (isTRUE(as.logical(cond))) 0 else NaN)

# _debug_assertion: as _nan_condition, and shows `message` when the
# condition fails and `log` is nonzero. The message is part of the function,
# one registered function per assertion.
.debug_assertion <- function(cond, message) {
  fn <- function(cond, log) {
    if (isTRUE(as.logical(cond))) return(0)
    if (isTRUE(as.logical(log))) base::message(message)
    NaN
  }
  f <- register_symbolic(paste0("_debug_assertion_", .next_id()), fn)
  f(cond, ASSERTION_LOG_VARIABLE)
}

# Assertions as a list of list(cond, msg). Accepted forms: pairs
# list(cond, "msg") (in either order) and named lists list("msg" = cond).
.assertion_pairs <- function(a) {
  out <- list()
  nms <- names(a) %||% rep("", length(a))
  for (i in seq_along(a)) {
    el <- a[[i]]
    if (is_sym(el)) {
      out[[length(out) + 1L]] <- list(cond = el, msg = nms[i])
    } else if (is.list(el) && length(el) == 2L) {
      j <- if (is_sym(el[[1L]])) 1L else 2L
      out[[length(out) + 1L]] <- list(cond = el[[j]], msg = as.character(el[[3L - j]]))
    } else {
      stop("an assertion must be list(condition, \"message\")", call. = FALSE)
    }
  }
  out
}

# get_assertions_expr: the sum of the assertion terms, NaN when one fails.
.assertions_expr <- function(sys) {
  pairs <- .assertion_pairs(assertions(sys))
  if (!length(pairs)) return(.S0())
  logged <- .contains_key(.parameters_of(sys, TRUE), ASSERTION_LOG_VARIABLE)
  terms <- lapply(pairs, function(pr) {
    if (logged) .debug_assertion(pr$cond, paste0("Assertion ", .fmt(pr$cond), " failed:\n", pr$msg))
    else .nan_condition(pr$cond)
  })
  Reduce(`+`, terms)
}

.mtk$hooks$rhs_assertions <- function(sys, rhs) {
  n <- length(rhs)
  rhs[[n]] <- rhs[[n]] + .assertions_expr(sys)
  rhs
}

# ---------------------------------------------------------------------------
# Logged functions

.LOGGED_R_FUNS <- list("^" = `^`, "/" = `/`, sqrt = sqrt, log = log, inv = function(x) 1 / x)

.logged_fun <- function(op, args, error_nonfinite) {
  f <- .LOGGED_R_FUNS[[op]] %||% {
    r <- .R_MATH[op]
    if (is.na(r)) stop("debug_system() cannot wrap the function ", op, call. = FALSE)
    get(r, envir = baseenv())
  }
  labels <- vapply(args, .fmt, "")
  label <- sprintf("%s(%s)", op, paste(labels, collapse = ", "))
  fn <- function(...) {
    vals <- list(...)
    inputs <- paste0("\n  ", labels, " => ", vapply(vals, function(v) format(v), ""), collapse = "")
    val <- tryCatch(suppressWarnings(do.call(f, vals)), error = function(e)
      stop("Function ", label, " errors with input", inputs, "\n", conditionMessage(e), call. = FALSE))
    if (error_nonfinite && any(!is.finite(val)))
      stop("Function ", label, " output non-finite value ", format(val), " with input", inputs, call. = FALSE)
    val
  }
  tag <- switch(op, "^" = "pow", "/" = "div", op)
  do.call(register_symbolic(paste0("logged_", tag, "_", .next_id()), fn), lapply(args, as_sym))
}

# debug_sub: rebuild the expression bottom up, wrapping the selected
# operations.
.debug_sub <- function(ex, functions, error_nonfinite) {
  if (inherits(ex, "Equation"))
    return(Equation(.debug_sub(ex$lhs, functions, error_nonfinite), .debug_sub(ex$rhs, functions, error_nonfinite)))
  if (!is_sym(ex) || !length(.args(ex)) || .type(ex) == "Derivative" || .is_declared(ex)) return(ex)
  args <- lapply(.args(ex), .debug_sub, functions = functions, error_nonfinite = error_nonfinite)
  ty <- .type(ex)
  if (ty == "Pow") {
    b <- args[[1L]]; e <- args[[2L]]
    ev <- if (.is_num(e)) as.numeric(e) else NA
    if (!is.na(ev) && ev == 0.5 && "sqrt" %in% functions) return(.logged_fun("sqrt", list(b), error_nonfinite))
    if (!is.na(ev) && ev < 0 && any(c("/", "inv") %in% functions)) {
      den <- if (ev == -1) b else .debug_sub(b^(-ev), functions, error_nonfinite)
      return(.logged_fun("/", list(1, den), error_nonfinite))
    }
    if ("^" %in% functions) return(.logged_fun("^", list(b, e), error_nonfinite))
  }
  if (ty == "Function" && .name(ex) %in% functions && !.is_registered_fun(ex))
    return(.logged_fun(.name(ex), args, error_nonfinite))
  .rebuild(ex, args)
}

#' Debug numerical errors in a system
#'
#' `debug_system(sys)` wraps the functions listed in `functions` in the
#' equations and observed equations of `sys`, so that an error, or a
#' nonfinite result when `error_nonfinite = TRUE`, stops with the symbolic
#' call and its numeric inputs, for example
#' `Function sqrt(u1(t)) output non-finite value NaN with input u1(t) => -0.1`.
#' `"/"` and `"inv"` wrap divisions (negative powers), `"sqrt"` square roots
#' and `"^"` other powers.
#'
#' It also adds the parameter [ASSERTION_LOG_VARIABLE], with value 1, so the
#' message of every failing assertion of the system is shown when the
#' right-hand side is evaluated. Assertions are given to [System()] as
#' `assertions = list(list(condition, "message"))`; the generated right-hand
#' side becomes `NaN` while one fails, which makes the solver reject the
#' step.
#'
#' @param sys A `System` without subsystems (flatten or compile it first).
#' @param functions Names of the functions to wrap.
#' @param error_nonfinite Also stop on nonfinite results.
#' @return The modified system.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' u <- Variable("u", t)
#' sys <- mtkcompile(System(list(D(u) ~ -sqrt(u)), t, name = "sys"))
#' dsys <- debug_system(sys, functions = "sqrt")
#' prob <- ODEProblem(dsys, list(u ~ 1), c(0, 10))
#' try(solve(prob))
#' @export
debug_system <- function(sys, functions = c("log", "sqrt", "^", "/", "inv", "asin", "acos"),
                         error_nonfinite = TRUE) {
  if (length(.f(sys, "systems")))
    stop("debug_system(sys) only works on systems with no sub-systems! Consider flattening it ",
         "with flatten(sys) or mtkcompile(sys) first.", call. = FALSE)
  functions <- as.character(functions)
  sub <- function(e) if (inherits(e, "Equation")) .debug_sub(e, functions, error_nonfinite) else e
  sys <- .update_system(sys,
    eqs = lapply(.f(sys, "eqs"), sub),
    ps = .unique_syms(c(.f(sys, "ps"), list(ASSERTION_LOG_VARIABLE))),
    initial_conditions = smap_set(.f(sys, "initial_conditions"), ASSERTION_LOG_VARIABLE, 1),
    observed = lapply(.f(sys, "observed"), sub))
  if (iscomplete(sys)) sys <- complete(sys)
  sys
}
