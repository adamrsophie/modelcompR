# Symbolic variables, parameters and their metadata, after
# ModelingToolkitBase/src/variables.jl and parameters.jl.
#
# A time dependent variable x(t) is the symcompR application of the undefined
# function "x" to the independent variable; a parameter or time independent
# variable is a symcompR Symbol. Metadata lives in an attribute on the
# returned object and in the package registry keyed by name, so it can be
# recovered from bare leaves found inside expressions.

.META_FIELDS <- c("default", "guess", "bounds", "description", "input", "output",
                  "irreducible", "state_priority", "connect", "unit", "disturbance",
                  "tunable", "dist", "nominal", "misc", "scope", "binary", "integer",
                  "variable_type", "dims", "callable")

.register_var <- function(name, kind, meta = list(), ivs = NULL) {
  # One entry per form, so a name can be an independent variable `t` and an
  # unknown `t(u)` at once, as after change_independent_variable(); the top
  # level fields are those of the latest registration.
  form <- if (length(ivs)) "Function" else "Symbol"
  forms <- .mtk$vars[[name]]$forms %||% list()
  forms[[form]] <- list(kind = kind, meta = meta, ivs = ivs)
  assign(name, list(kind = kind, meta = meta, ivs = ivs, forms = forms), envir = .mtk$vars)
  invisible(NULL)
}

.make_var <- function(name, kind, iv = NULL, meta = list()) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name))
    stop("a variable name must be a non-empty string", call. = FALSE)
  ivs <- if (is.null(iv)) list() else if (is_sym(iv)) list(iv) else lapply(iv, as_sym)
  meta <- meta[!vapply(meta, is.null, TRUE)]
  .register_var(name, kind, meta, ivs)
  e <- if (length(ivs)) do.call(Function(name), ivs) else Symbol(name)
  attr(e, "mtk_meta") <- c(list(kind = kind), meta)
  e
}

.meta_list <- function(default = NULL, guess = NULL, bounds = NULL, description = NULL,
                       input = NULL, output = NULL, irreducible = NULL, state_priority = NULL,
                       connect = NULL, unit = NULL, disturbance = NULL, tunable = NULL,
                       dist = NULL, nominal = NULL, misc = NULL, scope = NULL, binary = NULL,
                       integer = NULL, variable_type = NULL) {
  m <- list(default = default, guess = guess, bounds = bounds, description = description,
            input = input, output = output, irreducible = irreducible,
            state_priority = state_priority, connect = connect, unit = unit,
            disturbance = disturbance, tunable = tunable, dist = dist, nominal = nominal,
            misc = misc, scope = scope, binary = binary, integer = integer,
            variable_type = variable_type)
  m[!vapply(m, is.null, TRUE)]
}

#' Symbolic variables and parameters
#'
#' `Variable()` declares one unknown, `Parameter()` one parameter and
#' `Constant()` one constant, each with optional metadata. They correspond to
#' the `@variables`, `@parameters` and `@constants` macros of
#' ModelingToolkit: `Variable("x", t, default = 1)` is `@variables x(t) = 1`
#' and `Parameter("k", default = 0.5)` is `@parameters k = 0.5`.
#'
#' A default that is a number becomes an initial condition of any system the
#' variable is used in. A default that is a symbolic expression, or `NA`
#' (Julia's `missing`), becomes a binding: the variable is then determined by
#' that expression, or solved for during initialization.
#'
#' With `dims`, an array of scalar variables named `x[1]`, `x[2]`, ... is
#' returned as a list.
#'
#' @param name The name, a string.
#' @param iv The independent variable, or a list of independent variables,
#'   for dependent variables such as `x(t)` or `u(t, x)`. `NULL` declares a
#'   time independent variable.
#' @param default Default value: a number, a symbolic expression or `NA`.
#' @param guess Initial guess used by initialization.
#' @param bounds Numeric vector `c(lower, upper)`.
#' @param description A string.
#' @param input,output Mark the variable as a system input or output.
#' @param irreducible Keep the variable as an unknown during `mtkcompile()`.
#' @param state_priority Number; higher values are preferred as states.
#' @param connect Connector type for connector variables: `"Flow"`,
#'   `"Stream"` or `"Equality"`.
#' @param unit A unit, as a string.
#' @param disturbance Mark as a disturbance input.
#' @param tunable For parameters: whether the parameter is tunable (default
#'   `TRUE`).
#' @param dims Array dimensions, for array variables.
#' @param ... Further metadata: `dist`, `nominal`, `misc`, `scope`, `binary`,
#'   `integer`.
#' @return A `sym` expression carrying its metadata, or a list for arrays.
#' @seealso [variables()] declares several variables from strings.
#' @examples
#' t <- t_nounits
#' x <- Variable("x", t, default = 1, description = "position")
#' k <- Parameter("k", default = 0.5)
#' getdefault(x)
#' getdescription(x)
#' @export
Variable <- function(name, iv = NULL, default = NULL, guess = NULL, bounds = NULL,
                     description = NULL, input = NULL, output = NULL, irreducible = NULL,
                     state_priority = NULL, connect = NULL, unit = NULL, disturbance = NULL,
                     tunable = NULL, dims = NULL, ...) {
  # symcompR's Variable() is the codegen AST node of SymPy (a symbol, a type,
  # a value and attributes): hand it its own calls, so attaching both
  # packages is safe
  if ((missing(name) || is_sym(name) || any(names(list(...)) %in% c("symbol", "type", "value", "attrs"))) &&
      "Variable" %in% getNamespaceExports("symcompR")) {
    cl <- sys.call()
    cl[[1L]] <- getExportedValue("symcompR", "Variable")
    return(eval(cl, parent.frame()))
  }
  extra <- list(...)
  meta <- c(.meta_list(default = default, guess = guess, bounds = bounds,
                       description = description, input = input, output = output,
                       irreducible = irreducible, state_priority = state_priority,
                       connect = connect, unit = unit, disturbance = disturbance,
                       tunable = tunable), extra)
  .make_array_or_scalar(name, "unknown", iv, meta, dims)
}

#' @rdname Variable
#' @export
Parameter <- function(name, iv = NULL, default = NULL, guess = NULL, bounds = NULL,
                      description = NULL, tunable = NULL, unit = NULL, dims = NULL, ...) {
  meta <- c(.meta_list(default = default, guess = guess, bounds = bounds,
                       description = description, tunable = tunable, unit = unit), list(...))
  .make_array_or_scalar(name, "parameter", iv, meta, dims)
}

#' @rdname Variable
#' @export
Constant <- function(name, default, description = NULL, unit = NULL) {
  if (missing(default)) stop("a constant needs a value", call. = FALSE)
  meta <- .meta_list(default = default, description = description, unit = unit, tunable = FALSE)
  .make_var(name, "constant", NULL, meta)
}

#' @rdname Variable
#' @export
Brownian <- function(name, description = NULL) {
  .make_var(name, "brownian", NULL, .meta_list(description = description))
}

#' @rdname Variable
#' @param rate The rate of a Poisson counting process.
#' @export
Poissonian <- function(name, rate, description = NULL) {
  .make_var(name, "poissonian", NULL, c(.meta_list(description = description), list(rate = rate)))
}

#' @rdname Variable
#' @export
Discrete <- function(name, iv, default = NULL, description = NULL, ...) {
  meta <- c(.meta_list(default = default, description = description), list(...))
  .make_var(name, "discrete", iv, meta)
}

.make_array_or_scalar <- function(name, kind, iv, meta, dims) {
  if (is.null(dims)) return(.make_var(name, kind, iv, meta))
  dims <- as.integer(dims)
  idx <- as.matrix(expand.grid(lapply(dims, seq_len)))
  default <- meta$default
  out <- vector("list", nrow(idx))
  for (i in seq_len(nrow(idx))) {
    m <- meta
    if (length(default) > 1L) m$default <- default[[i]]
    m$array_parent <- name
    m$array_index <- idx[i, ]
    out[[i]] <- .make_var(paste0(name, "[", paste(idx[i, ], collapse = ","), "]"), kind, iv, m)
  }
  if (length(dims) == 2L) dim(out) <- dims
  .arr_install_ops_hook()
  structure(out, class = c("mtk_array", "list"), mtk_name = name)
}

#' @export
print.mtk_array <- function(x, ...) {
  cat("Array variable ", attr(x, "mtk_name"), " with ", length(x), " elements:\n", sep = "")
  cat(vapply(unclass(x), format, ""), sep = "\n")
  invisible(x)
}

#' @export
format.mtk_array <- function(x, ...) vapply(unclass(x), format, "")

#' The independent variable time and its derivative
#'
#' `t_nounits` is the independent variable `t` and `D_nounits` the
#' differential with respect to it, as in ModelingToolkit. Assign them to
#' short names in your own code: `t <- t_nounits; D <- D_nounits`.
#'
#' @format `t_nounits` is a `sym` symbol; `D_nounits` is a function of
#'   one argument, an expression.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t)
#' D(x)
#' @export
t_nounits <- .make_var("t", "iv", NULL, list(description = "time"))

#' @rdname t_nounits
#' @usage D_nounits
#' @export
D_nounits <- Differential(t_nounits)

#' Declare independent variables
#'
#' `independent_variable("t")` declares an independent variable, the
#' `@independent_variables` macro of ModelingToolkit. Given a system it
#' returns that system's independent variable instead.
#'
#' @param x A name, or a system.
#' @param ... Further names, when declaring.
#' @return A `sym` symbol, or a list of them.
#' @export
independent_variable <- function(x, ...) {
  if (inherits(x, "System")) return(get_iv(x))
  nms <- c(x, unlist(list(...)))
  out <- lapply(nms, function(n) .make_var(n, "iv", NULL, list()))
  if (length(out) == 1L) out[[1L]] else stats::setNames(out, nms)
}

# ---------------------------------------------------------------------------
# Declaring several variables from strings, the R counterpart of the
# @variables and @parameters macros.

.parse_decl <- function(s) {
  ex <- tryCatch(parse(text = s, keep.source = FALSE)[[1L]],
                 error = function(e) stop("cannot parse declaration '", s, "'", call. = FALSE))
  default <- NULL
  if (is.call(ex) && (identical(ex[[1L]], as.name("=")) || identical(ex[[1L]], as.name("<-")))) {
    default <- ex[[3L]]
    ex <- ex[[2L]]
  }
  dims <- NULL
  if (is.call(ex) && identical(ex[[1L]], as.name("["))) {
    dims <- vapply(as.list(ex)[-c(1L, 2L)], function(d) {
      d <- eval(d, baseenv())
      if (length(d) > 1L) length(d) else as.integer(d)
    }, 0L)
    ex <- ex[[2L]]
  }
  ivs <- NULL
  if (is.call(ex)) {
    ivs <- vapply(as.list(ex)[-1L], deparse, "")
    name <- as.character(ex[[1L]])
  } else {
    name <- as.character(ex)
  }
  list(name = name, ivs = ivs, default = default, dims = dims)
}

# Declaration strings are evaluated in the caller's environment, as the
# @variables macro would evaluate its default expressions; they are the
# caller's own code.
.declare <- function(specs, kind, envir, meta, iv_env) {
  specs <- unlist(specs)
  out <- list()
  for (s in specs) {
    d <- .parse_decl(s)
    m <- meta
    if (!is.null(d$default)) {
      dv <- eval(d$default, iv_env)
      m$default <- if (is.numeric(dv) || is_sym(dv) || (length(dv) == 1L && is.na(dv))) dv else as_sym(dv)
    }
    iv <- NULL
    if (length(d$ivs)) {
      iv <- lapply(d$ivs, function(n) {
        v <- get0(n, envir = iv_env)
        if (is.null(v) || !is_sym(v)) v <- .make_var(n, "iv", NULL, list())
        v
      })
    }
    v <- .make_array_or_scalar(d$name, kind, iv, m, d$dims)
    assign(d$name, v, envir = envir)
    out[[d$name]] <- v
  }
  invisible(out)
}

#' Declare several variables or parameters at once
#'
#' `variables("x(t) = 1", "y(t)")` declares `x(t)` with default 1 and `y(t)`
#' and assigns them to `x` and `y` in the calling environment, like the
#' `@variables` macro. `parameters("k = 0.5", "m")` does the same for
#' parameters, and `brownians()`, `constants()` and `discretes()` for the other
#' kinds. The independent variable named inside the parentheses is looked up
#' in the calling environment, and created if missing. Array variables are
#' written `"x(t)[1:3]"`.
#'
#' Given a system, `parameters()` and `brownians()` return its parameters and
#' brownian variables instead, as the accessors of the same name in
#' ModelingToolkit.
#'
#' @param ... Declaration strings, or a system for the accessor forms.
#' @param .meta Metadata applied to every declared variable, as a named list.
#' @param .env Environment to assign into.
#' @return Invisibly, a named list of the declared variables.
#' @examples
#' t <- t_nounits
#' variables("x(t) = 1", "y(t)")
#' parameters("k = 0.5", "m")
#' x; k
#' @export
variables <- function(..., .meta = list(), .env = parent.frame()) {
  .declare(list(...), "unknown", .env, .meta, .env)
}

#' @rdname variables
#' @export
parameters <- function(..., .meta = list(), .env = parent.frame()) {
  dots <- list(...)
  if (length(dots) && inherits(dots[[1L]], "System")) return(do.call(.parameters_of, dots))
  .declare(dots, "parameter", .env, .meta, .env)
}

#' @rdname variables
#' @export
constants <- function(..., .env = parent.frame()) {
  .declare(list(...), "constant", .env, list(tunable = FALSE), .env)
}

#' @rdname variables
#' @export
brownians <- function(..., .env = parent.frame()) {
  dots <- list(...)
  if (length(dots) && inherits(dots[[1L]], "System")) return(get_brownians(dots[[1L]]))
  .declare(dots, "brownian", .env, list(), .env)
}

#' @rdname variables
#' @export
discretes <- function(..., .meta = list(), .env = parent.frame()) {
  .declare(list(...), "discrete", .env, .meta, .env)
}

# ---------------------------------------------------------------------------
# Metadata access

.meta <- function(x) {
  x <- as_sym(x)
  if (.type(x) == "Derivative") return(list())
  m <- attr(x, "mtk_meta")
  if (!is.null(m)) return(m)
  ent <- .mtk$vars[[.name(x) %||% ""]]
  if (is.null(ent)) return(list())
  ent <- ent$forms[[.type(x)]] %||% ent
  c(list(kind = ent$kind), ent$meta)
}

.meta_get <- function(x, field, default = NULL) .meta(x)[[field]] %||% default

#' Set metadata on a variable
#'
#' Returns a copy of the variable with metadata `field` set to `value`. The
#' package registry entry for the name is updated as well.
#'
#' @param x A variable.
#' @param field Metadata field name, such as `"default"`, `"guess"` or
#'   `"description"`.
#' @param value The value.
#' @return The variable with updated metadata.
#' @export
setmetadata <- function(x, field, value) {
  m <- .meta(x)
  m[[field]] <- value
  attr(x, "mtk_meta") <- m
  nm <- .name(x)
  ent <- .mtk$vars[[nm]]
  if (!is.null(ent)) {
    ent$meta[[field]] <- value
    if (!is.null(ent$forms[[.type(x)]])) ent$forms[[.type(x)]]$meta[[field]] <- value
    assign(nm, ent, envir = .mtk$vars)
  }
  x
}

#' @rdname setmetadata
#' @param default Value returned when the field is not set.
#' @export
getmetadata <- function(x, field, default = NULL) .meta_get(x, field, default)

#' Variable metadata accessors
#'
#' Query the metadata attached to a variable when it was declared, as the
#' functions of the same names in ModelingToolkit.
#'
#' @param x A variable.
#' @param value A new value, for the setters.
#' @return The metadata value, or a logical for the `has`/`is` forms.
#' @name metadata
#' @examples
#' x <- Variable("x", t_nounits, default = 1, bounds = c(0, 2), input = TRUE)
#' getdefault(x); getbounds(x); isinput(x)
NULL

#' @rdname metadata
#' @export
getdefault <- function(x) {
  d <- .meta_get(x, "default")
  if (is.null(d)) stop("variable ", format(x), " has no default", call. = FALSE)
  d
}
#' @rdname metadata
#' @export
hasdefault <- function(x) !is.null(.meta_get(x, "default"))
#' @rdname metadata
#' @export
setdefault <- function(x, value) setmetadata(x, "default", value)
#' @rdname metadata
#' @export
getguess <- function(x) .meta_get(x, "guess")
#' @rdname metadata
#' @export
hasguess <- function(x) !is.null(.meta_get(x, "guess"))
#' @rdname metadata
#' @export
setguess <- function(x, value) setmetadata(x, "guess", value)
#' @rdname metadata
#' @export
getbounds <- function(x) .meta_get(x, "bounds", c(-Inf, Inf))
#' @rdname metadata
#' @export
hasbounds <- function(x) !is.null(.meta_get(x, "bounds"))
#' @rdname metadata
#' @export
getnominal <- function(x) .meta_get(x, "nominal", 1)
#' @rdname metadata
#' @export
hasnominal <- function(x) !is.null(.meta_get(x, "nominal"))
#' @rdname metadata
#' @export
setnominal <- function(x, value) setmetadata(x, "nominal", value)
#' @rdname metadata
#' @export
getdescription <- function(x) .meta_get(x, "description", "")
#' @rdname metadata
#' @export
hasdescription <- function(x) nzchar(getdescription(x))
#' @rdname metadata
#' @export
isinput <- function(x) isTRUE(.meta_get(x, "input"))
#' @rdname metadata
#' @export
isoutput <- function(x) isTRUE(.meta_get(x, "output"))
#' @rdname metadata
#' @export
isdisturbance <- function(x) isTRUE(.meta_get(x, "disturbance"))
#' @rdname metadata
#' @export
isirreducible <- function(x) isTRUE(.meta_get(x, "irreducible"))
#' @rdname metadata
#' @export
state_priority <- function(x) .meta_get(x, "state_priority", 0)
#' @rdname metadata
#' @export
istunable <- function(x) {
  if (!isparameter(x)) return(FALSE)
  tu <- .meta_get(x, "tunable")
  if (is.null(tu)) TRUE else isTRUE(tu)
}
#' @rdname metadata
#' @export
getdist <- function(x) .meta_get(x, "dist")
#' @rdname metadata
#' @export
hasdist <- function(x) !is.null(.meta_get(x, "dist"))
#' @rdname metadata
#' @export
getunit <- function(x) .meta_get(x, "unit")
#' @rdname metadata
#' @export
hasunit <- function(x) !is.null(.meta_get(x, "unit"))
#' @rdname metadata
#' @export
getconnect <- function(x) .meta_get(x, "connect")
#' @rdname metadata
#' @export
hasconnect <- function(x) !is.null(.meta_get(x, "connect"))
#' @rdname metadata
#' @export
getmisc <- function(x) .meta_get(x, "misc")
#' @rdname metadata
#' @export
hasmisc <- function(x) !is.null(.meta_get(x, "misc"))
#' @rdname metadata
#' @export
isparameter <- function(x) {
  k <- .var_kind(as_sym(x))
  !is.null(k) && k %in% c("parameter", "constant", "discrete")
}
#' @rdname metadata
#' @export
isvariable <- function(x) identical(.var_kind(as_sym(x)), "unknown")
#' @rdname metadata
#' @export
isconstant <- function(x) {
  # SciMLOperators: an operator whose coefficients do not change
  if (inherits(x, "AbstractSciMLOperator")) return(isTRUE(x$constant))
  identical(.var_kind(as_sym(x)), "constant")
}
#' @rdname metadata
#' @export
isbrownian <- function(x) identical(.var_kind(as_sym(x)), "brownian")

#' Convert between variables and parameters
#'
#' `toparam(x)` re-declares a variable as a parameter, `tovar(x)` a
#' parameter as a variable and `tobrownian(x)` as a brownian variable,
#' keeping the name and metadata.
#'
#' @param x A variable or parameter.
#' @return The converted variable.
#' @export
toparam <- function(x) .rekind(x, "parameter")
#' @rdname toparam
#' @export
tovar <- function(x) .rekind(x, "unknown")
#' @rdname toparam
#' @export
tobrownian <- function(x) .rekind(x, "brownian")

.rekind <- function(x, kind) {
  x <- as_sym(x)
  m <- .meta(x)
  m$kind <- NULL
  ivs <- if (.type(x) == "Function") .args(x) else NULL
  .make_var(.name(x), kind, ivs, m)
}

# ---------------------------------------------------------------------------
# Scopes for hierarchical models, after the SymScope types.

#' Variable scopes
#'
#' By default a variable used in a subsystem is namespaced by that
#' subsystem's name. `ParentScope(x)` makes it refer to the variable of the
#' same name one level up, and `GlobalScope(x)` to the top level variable.
#' `LocalScope(x)` restores the default.
#'
#' @param x A variable or parameter.
#' @return The variable with its scope metadata set.
#' @export
LocalScope <- function(x) .with_scope(x, list(type = "local"))
#' @rdname LocalScope
#' @export
ParentScope <- function(x) {
  sc <- .meta_get(x, "scope")
  depth <- if (!is.null(sc) && identical(sc$type, "parent")) sc$depth + 1L else 1L
  .with_scope(x, list(type = "parent", depth = depth))
}
#' @rdname LocalScope
#' @export
GlobalScope <- function(x) .with_scope(x, list(type = "global"))

.with_scope <- function(x, scope) {
  x <- as_sym(x)
  if (.type(x) == "Derivative") {
    inner <- .with_scope(.args(x)[[1L]], scope)
    di <- .deriv_info(x)
    return(.make_derivative(inner, di$iv, di$order))
  }
  # The scope is kept on this object only, not in the registry, because the
  # same name may be local in the parent. System() records the scope of each
  # leaf it finds; see .collect_scopes().
  m <- .meta(x)
  m$scope <- scope
  attr(x, "mtk_meta") <- m
  x
}
