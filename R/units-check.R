# Units and validation, after ModelingToolkitBase/src/systems/unit_check.jl,
# the DynamicQuantities extension ext/MTKDynamicQuantitiesExt.jl (get_unit,
# safe_get_unit, _validate, validate, check_units) and the CheckAll,
# CheckNone, CheckComponents and CheckUnits constants of utils.jl.
#
# Units are given as strings in the `unit` metadata of variables ("m/s^2",
# "kJ", "MW", "ms"). A parsed unit is a scale factor relative to SI and the
# exponents of the seven SI base dimensions. As with DynamicQuantities, two
# units are equivalent only when both agree: kJ and J are not equivalent.

#' System check flags
#'
#' Values for the `checks` argument of [System()]: `CheckNone` disables input
#' validation, `CheckAll` enables all of it (as `TRUE`), and
#' `CheckComponents` and `CheckUnits` enable the structural checks and the
#' unit checks alone. Combine flags with `bitwOr()`.
#'
#' @format Integers.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t, unit = "m")
#' v <- Parameter("v", unit = "s")
#' try(System(list(D(x) ~ v), t, name = "bad"))
#' sys <- System(list(D(x) ~ v), t, name = "unchecked", checks = CheckComponents)
#' @name check_flags
NULL

#' @rdname check_flags
#' @export
CheckNone <- 0L
#' @rdname check_flags
#' @export
CheckAll <- 1L
#' @rdname check_flags
#' @export
CheckComponents <- 2L
#' @rdname check_flags
#' @export
CheckUnits <- 4L

.SI_DIMS <- c("m", "kg", "s", "A", "K", "mol", "cd")

# `.scale` (not `scale`) so that the dimension `s = 1` is not partially
# matched to it.
.unit <- function(.scale = 1, ...) {
  d <- stats::setNames(rep(0, 7L), .SI_DIMS)
  e <- c(...)
  if (length(e)) d[names(e)] <- e
  structure(list(scale = .scale, dims = d), class = "mtk_unit")
}

.unitless <- function() .unit()

.UNIT_TABLE <- local({
  u <- list(
    m = .unit(1, m = 1), g = .unit(1e-3, kg = 1), s = .unit(1, s = 1), A = .unit(1, A = 1),
    K = .unit(1, K = 1), mol = .unit(1, mol = 1), cd = .unit(1, cd = 1),
    Hz = .unit(1, s = -1), N = .unit(1, kg = 1, m = 1, s = -2),
    Pa = .unit(1, kg = 1, m = -1, s = -2), J = .unit(1, kg = 1, m = 2, s = -2),
    W = .unit(1, kg = 1, m = 2, s = -3), C = .unit(1, A = 1, s = 1),
    V = .unit(1, kg = 1, m = 2, s = -3, A = -1), F = .unit(1, kg = -1, m = -2, s = 4, A = 2),
    Ohm = .unit(1, kg = 1, m = 2, s = -3, A = -2), S = .unit(1, kg = -1, m = -2, s = 3, A = 2),
    Wb = .unit(1, kg = 1, m = 2, s = -2, A = -1), T = .unit(1, kg = 1, s = -2, A = -1),
    H = .unit(1, kg = 1, m = 2, s = -2, A = -2), L = .unit(1e-3, m = 3),
    min = .unit(60, s = 1), h = .unit(3600, s = 1), hr = .unit(3600, s = 1),
    d = .unit(86400, s = 1), bar = .unit(1e5, kg = 1, m = -1, s = -2),
    rad = .unit(1), sr = .unit(1), eV = .unit(0x1.7a4da290c1653p-63, kg = 1, m = 2, s = -2),
    M = .unit(1000, mol = 1, m = -3))
  u[["\u03A9"]] <- u$Ohm
  u
})

.UNIT_PREFIX <- c(Y = 1e24, Z = 1e21, E = 1e18, P = 1e15, T = 1e12, G = 1e9, M = 1e6, k = 1e3,
                  h = 1e2, da = 1e1, d = 1e-1, c = 1e-2, m = 1e-3, u = 1e-6, "\u03BC" = 1e-6,
                  "\u00B5" = 1e-6, n = 1e-9, p = 1e-12, f = 1e-15, a = 1e-18, z = 1e-21, y = 0x1.357c299a88ea7p-80)

.unit_mul <- function(a, b) structure(list(scale = a$scale * b$scale, dims = a$dims + b$dims),
                                      class = "mtk_unit")
.unit_pow <- function(a, n) structure(list(scale = a$scale^n, dims = a$dims * n), class = "mtk_unit")
.unit_div <- function(a, b) .unit_mul(a, .unit_pow(b, -1))

.unit_name <- function(nm) {
  hit <- .UNIT_TABLE[[nm]]
  if (!is.null(hit)) return(hit)
  for (p in names(.UNIT_PREFIX)) {
    if (startsWith(nm, p)) {
      base <- .UNIT_TABLE[[substring(nm, nchar(p) + 1L)]]
      if (!is.null(base)) return(.unit_mul(.unit(.UNIT_PREFIX[[p]]), base))
    }
  }
  .validation_error("unknown unit '", nm, "'")
}

# Recursive descent parser: expr = term (("*" | "/") term)*,
# term = atom ("^" number)?, atom = name | number | "(" expr ")".
.parse_unit <- function(s) {
  if (inherits(s, "mtk_unit")) return(s)
  if (is.null(s)) return(.unitless())
  if (!is.character(s) || length(s) != 1L) .validation_error("a unit must be a string")
  toks <- regmatches(s, gregexpr("[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?|[^0-9[:space:]*/^().+-]+|[*/^().+-]", s))[[1L]]
  pos <- 1L
  peek <- function() if (pos <= length(toks)) toks[pos] else ""
  take <- function() { tk <- toks[pos]; pos <<- pos + 1L; tk }
  number <- function() {
    sgn <- 1
    while (peek() %in% c("-", "+")) if (take() == "-") sgn <- -sgn
    if (peek() == "(") {
      take()
      a <- number()
      if (peek() == "/") { take(); a <- a / number() }
      if (take() != ")") .validation_error("cannot parse unit '", s, "'")
      return(sgn * a)
    }
    v <- suppressWarnings(as.numeric(take()))
    if (is.na(v)) .validation_error("cannot parse unit '", s, "'")
    sgn * v
  }
  atom <- function() {
    tk <- peek()
    if (tk == "(") {
      take()
      u <- expr()
      if (take() != ")") .validation_error("cannot parse unit '", s, "'")
      return(u)
    }
    if (grepl("^[0-9.]", tk)) return(.unit(number()))
    take()
    .unit_name(tk)
  }
  term <- function() {
    u <- atom()
    if (peek() == "^") {
      take()
      u <- .unit_pow(u, number())
    }
    u
  }
  expr <- function() {
    u <- term()
    while (peek() %in% c("*", "/", ".")) {
      op <- take()
      u <- if (op == "/") .unit_div(u, term()) else .unit_mul(u, term())
    }
    u
  }
  if (!length(toks)) return(.unitless())
  u <- expr()
  if (pos <= length(toks)) .validation_error("cannot parse unit '", s, "'")
  u
}

#' @export
format.mtk_unit <- function(x, ...) {
  d <- x$dims[abs(x$dims) > 1e-12]
  dims <- if (length(d)) paste0(names(d), ifelse(abs(d - 1) < 1e-12, "", paste0("^", formatC(d, format = "g"))),
                                collapse = " ") else "(dimensionless)"
  if (abs(x$scale - 1) > 1e-12) paste(format(x$scale), dims) else dims
}

#' @export
print.mtk_unit <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

.validation_error <- function(...) {
  msg <- paste0(...)
  stop(structure(class = c("ValidationError", "error", "condition"),
                 list(message = msg, call = NULL)))
}

.unit_equiv <- function(a, b) {
  abs(a$scale - b$scale) <= 1e-9 * max(abs(a$scale), abs(b$scale)) &&
    all(abs(a$dims - b$dims) < 1e-9)
}

.is_unitless <- function(u) .unit_equiv(u, .unitless())

# The unit of a variable: its own metadata, else the registry entry (a leaf
# may carry metadata attached before setmetadata() added the unit).
.unit_meta <- function(v) {
  u <- .meta_get(v, "unit")
  if (!is.null(u)) return(u)
  ent <- .mtk$vars[[.name(v) %||% ""]]
  if (is.null(ent)) return(NULL)
  (ent$forms[[.type(v)]] %||% ent)$meta$unit
}

.var_unit <- function(v) {
  u <- .unit_meta(v)
  if (is.null(u)) .unitless() else .parse_unit(u)
}

.TRANSCENDENTAL <- c("exp", "log", "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh",
                     "tanh", "asinh", "acosh", "atanh", "sec", "csc", "cot", "gamma", "loggamma",
                     "erf", "erfc", "factorial", "beta", "sinc")
.UNIT_PRESERVING <- c("Abs", "floor", "ceiling", "re", "im", "Max", "Min")
.RELATIONALS <- c("Equality", "Unequality", "StrictLessThan", "LessThan", "StrictGreaterThan",
                  "GreaterThan")

#' Units of expressions and unit validation
#'
#' `get_unit(x)` returns the unit of an expression from the `unit` metadata
#' of its variables, and signals a `ValidationError` when the units are
#' inconsistent (terms of a sum with different units, a dimensional argument
#' of `exp`, ...). Units are strings such as `"m/s^2"`, `"kJ"` or `"ms"`,
#' built from SI base and derived units with SI prefixes, products,
#' quotients and powers. As in ModelingToolkit, units must agree exactly:
#' `"kJ"` and `"J"` are different units.
#'
#' `validate(x)` checks an equation, a list of equations, an expression or a
#' system: it warns about each inconsistency and returns `FALSE`, or returns
#' `TRUE`. Terms that are identically zero are exempt, so `0 ~ P` validates.
#' [System()] runs these checks when a variable has a unit, unless the
#' `checks` flags exclude [CheckUnits].
#'
#' @param x An expression, an equation, a list of them or a `System`.
#' @param info A prefix for the warning messages.
#' @return `get_unit()` a unit object; `validate()` a logical.
#' @examples
#' tm <- independent_variable("tm"); tm <- setmetadata(tm, "unit", "ms")
#' tau <- Parameter("tau", unit = "ms")
#' E <- Variable("E", tm, unit = "kJ"); P <- Variable("P", tm, unit = "MW")
#' D <- Differential(tm)
#' get_unit(E / tau)
#' validate(list(D(E) ~ P - E / tau, 0 ~ P))
#' @export
get_unit <- function(x) {
  if (inherits(x, "mtk_unit")) return(x)
  if (is.numeric(x) || is.logical(x)) return(.unitless())
  if (!is_sym(x)) .validation_error("unsupported value of class ", class(x)[1L])
  ty <- .type(x)
  if (ty %in% c("Integer", "Rational", "Float", "Constant", "BooleanAtom", "NumberSymbol")) return(.unitless())
  if (.is_op(x)) return(get_unit(.args(x)[[1L]]))
  if (ty %in% c("Symbol", "Function") && .is_declared(x)) return(.var_unit(x))
  a <- .args(x)
  switch(ty,
    Symbol = , Dummy = .unitless(),
    Derivative = {
      u <- get_unit(a[[1L]])
      for (spec in a[-1L]) {
        sa <- .args(spec)
        u <- .unit_div(u, .unit_pow(get_unit(sa[[1L]]), as.numeric(sa[[2L]])))
      }
      u
    },
    Add = {
      us <- lapply(a, get_unit)
      for (u in us[-1L]) if (!.unit_equiv(u, us[[1L]]))
        .validation_error(", in sum ", format(x), ", units [",
                          paste(vapply(us, format.mtk_unit, ""), collapse = ", "), "] do not match.")
      us[[1L]]
    },
    Mul = Reduce(.unit_mul, lapply(a, get_unit)),
    Pow = {
      b <- get_unit(a[[1L]])
      e <- get_unit(a[[2L]])
      if (!.is_unitless(e)) .validation_error(", in ", format(x), ", the exponent is not dimensionless.")
      if (.is_unitless(b)) return(.unitless())
      if (!.is_num(a[[2L]]))
        .validation_error(", in ", format(x), ", a unitful base has a symbolic exponent.")
      .unit_pow(b, as.numeric(a[[2L]]))
    },
    Piecewise = {
      vals <- list()
      for (pr in a) {
        pa <- .args(pr)
        if (!.is_unitless(get_unit(pa[[2L]])))
          .validation_error(", in ", format(x), ", a condition is not dimensionless.")
        vals[[length(vals) + 1L]] <- pa[[1L]]
      }
      vals <- Filter(function(v) !(.is_num(v) && is.nan(as.numeric(v))), vals)
      us <- lapply(vals, get_unit)
      for (u in us[-1L]) if (!.unit_equiv(u, us[[1L]]))
        .validation_error(", in ", format(x), ", the branches have different units.")
      if (length(us)) us[[1L]] else .unitless()
    },
    Function = {
      nm <- .name(x)
      us <- lapply(a, get_unit)
      if (nm %in% .TRANSCENDENTAL) {
        for (u in us) if (!.is_unitless(u))
          .validation_error(", in ", format(x), ", the argument [", format.mtk_unit(u), "] is not dimensionless.")
        return(.unitless())
      }
      if (nm %in% .UNIT_PRESERVING) {
        for (u in us[-1L]) if (!.unit_equiv(u, us[[1L]]))
          .validation_error(", in ", format(x), ", units do not match.")
        return(us[[1L]])
      }
      if (nm %in% c("sign", "Heaviside", "atan2")) return(.unitless())
      if (nm == "Mod") return(us[[1L]])
      .validation_error("Unable to get unit for operation ", nm, " with arguments ",
                        paste(vapply(a, format, ""), collapse = ", "), ".")
    },
    {
      if (ty %in% .RELATIONALS) {
        us <- lapply(a, get_unit)
        if (!.unit_equiv(us[[1L]], us[[2L]]))
          .validation_error(", in comparison ", format(x), ", units [", format.mtk_unit(us[[1L]]), "] and [",
                            format.mtk_unit(us[[2L]]), "] do not match.")
        return(.unitless())
      }
      if (ty %in% c("And", "Or", "Not")) {
        for (z in a) get_unit(z)
        return(.unitless())
      }
      .validation_error("Unable to get unit for ", format(x), ".")
    })
}

# get_unit() returning NULL and a warning instead of an error.
.safe_get_unit <- function(term, info) {
  tryCatch(get_unit(term), ValidationError = function(e) {
    warning(info, conditionMessage(e), call. = FALSE)
    NULL
  })
}

.validate_terms <- function(terms, labels, info = "") {
  valid <- TRUE
  first <- NULL
  first_label <- NULL
  for (i in seq_along(terms)) {
    u <- .safe_get_unit(terms[[i]], paste0(info, " ", labels[i]))
    if (is.null(u)) {
      valid <- FALSE
      next
    }
    z <- terms[[i]]
    if ((is.numeric(z) && all(z == 0)) || .is_zero(z)) next
    if (is.null(first)) {
      first <- u
      first_label <- labels[i]
    } else if (!.unit_equiv(first, u)) {
      valid <- FALSE
      msg <- sprintf("%s: units [%s] for %s and [%s] for %s do not match.", info, format.mtk_unit(first),
                     first_label, format.mtk_unit(u), labels[i])
      if (all(abs(first$dims - u$dims) < 1e-9))
        msg <- paste(msg, "The units have the same dimensions but different scales.")
      warning(msg, call. = FALSE)
    }
  }
  valid
}

#' @rdname get_unit
#' @export
validate <- function(x, info = "") {
  if (inherits(x, "System")) {
    noise <- .f(x, "noise_eqs")
    eqs <- Filter(function(e) inherits(e, "Equation"), equations(x))
    ok <- validate(eqs, info)
    if (!is.null(noise)) ok <- .validate_noise(eqs, noise, info) && ok
    cons <- constraints(x)
    if (length(cons)) ok <- validate(cons, paste0(info, " in constraints")) && ok
    return(ok)
  }
  if (inherits(x, "formula")) x <- as_equation(x)
  if (inherits(x, c("Equation", "Inequality"))) {
    if (!is_sym(x$lhs) || !is_sym(x$rhs)) return(TRUE)
    return(.validate_terms(list(x$lhs, x$rhs), c("left", "right"), info))
  }
  if (is.list(x) && !is_sym(x)) {
    if (!length(x)) return(TRUE)
    res <- vapply(seq_along(x), function(i) validate(x[[i]], paste0(info, " in eq. #", i)), TRUE)
    return(all(res))
  }
  if (inherits(x, "Connection")) return(TRUE)
  !is.null(.safe_get_unit(x, info))
}

# Each noise term must have the units of its equation (validate(eqs, noise)
# in the extension); a matrix row holds the terms of one equation.
.validate_noise <- function(eqs, noise, info = "") {
  rows <- if (is.matrix(noise)) lapply(seq_len(nrow(noise)), function(i) noise[i, ]) else as.list(noise)
  if (length(rows) != length(eqs)) return(TRUE)
  all(vapply(seq_along(eqs), function(i) {
    terms <- if (is_sym(rows[[i]])) list(rows[[i]]) else as.list(rows[[i]])
    .validate_terms(c(list(eqs[[i]]$lhs, eqs[[i]]$rhs), terms),
                    c("left", "right", paste0("noise #", seq_along(terms))), paste0(info, " in eq. #", i))
  }, TRUE))
}

# check_units: run by System() through the check_units hook when a variable
# of the system has a unit (as __get_unit_type, which only checks systems
# with DynamicQuantities units).
.has_units <- function(vars) any(vapply(vars, function(v) !is.null(.unit_meta(v)), TRUE))

# Only the system's own equations and constraints: subsystems were checked
# when they were built.
.check_units_hook <- function(sys) {
  vars <- c(.f(sys, "unknowns"), .f(sys, "ps"), independent_variables(sys))
  if (!length(vars) || !.has_units(vars)) return(invisible(TRUE))
  eqs <- Filter(function(e) inherits(e, "Equation"), .f(sys, "eqs"))
  ok <- validate(eqs)
  if (!is.null(.f(sys, "noise_eqs"))) ok <- .validate_noise(eqs, .f(sys, "noise_eqs")) && ok
  if (length(.f(sys, "constraints"))) ok <- validate(.f(sys, "constraints"), " in constraints") && ok
  if (!ok) .validation_error("Some equations had invalid units. See warnings for details.")
  invisible(TRUE)
}

.mtk$hooks$check_units <- function(sys) .check_units_hook(sys)
