# Components, connectors and the model DSL, after
# ModelingToolkitBase/src/systems/abstractsystem.jl (@component, @connector,
# component_post_processing) and SciCompDSL/src/model_parsing.jl
# (@mtkmodel, @connector with a body, Model).
#
# The macros become functions: `component(f)` and `connector(f)` wrap a
# constructor function, and `mtkmodel()` builds a constructor from
# declarations, the counterpart of the @mtkmodel blocks.

.wrap_constructor <- function(f, isconnector, type) {
  force(f)
  out <- function(...) {
    res <- f(...)
    if (!inherits(res, "System")) return(res)
    md <- .f(res, "metadata")
    if (!is.null(md$gui_metadata)) return(res)
    if (isconnector) res <- .update_system(res, connector_type = connector_type(res))
    if (!is.null(type)) {
      md$gui_metadata <- type
      res <- .update_system(res, metadata = md)
    }
    res
  }
  structure(out, class = c(if (isconnector) "ConnectorConstructor" else "ComponentConstructor", "function"),
            type = type)
}

#' Components and connectors
#'
#' `component(f)` marks a constructor function `f(name, ...)` returning a
#' `System` as a component, the `@component` macro of ModelingToolkit: the
#' systems it builds record the constructor name in their metadata.
#' `connector(f)` does the same for connectors (`@connector`) and also sets
#' the connector type (see [connector_type()]), which [connect()] needs.
#' `connector(sys)` marks an existing system as a connector, and
#' `connector("Name", variables = ...)` builds a connector constructor from
#' declarations as [mtkmodel()] does.
#'
#' @param x,f A constructor function, a `System`, or a model name.
#' @param type The constructor name recorded in the metadata; by default
#'   the name of the argument.
#' @param ... For the name form, arguments of [mtkmodel()].
#' @return A constructor function, or a `System`.
#' @examples
#' t <- t_nounits
#' Pin <- connector(function(name) {
#'   v <- Variable("v", t); i <- Variable("i", t, connect = Flow)
#'   System(list(), t, list(v, i), list(), name = name)
#' })
#' p <- Pin(name = "p")
#' connector_type(p)
#' Pin2 <- connector("Pin2", variables = list("v(t)", list("i(t)", connect = Flow)))
#' Pin2(name = "q")
#' @export
component <- function(f, type = NULL) {
  if (is.null(type) && is.name(substitute(f))) type <- as.character(substitute(f))
  .wrap_constructor(f, FALSE, type)
}

#' @rdname component
#' @export
connector <- function(x, ..., type = NULL) {
  if (inherits(x, "System")) return(.update_system(x, connector_type = connector_type(x)))
  if (is.character(x)) return(mtkmodel(x, ..., connector = TRUE, .env = parent.frame()))
  if (is.null(type) && is.name(substitute(x))) type <- as.character(substitute(x))
  .wrap_constructor(x, TRUE, type)
}

# Evaluate a DSL block in the instance environment: a quoted expression, or
# a function without arguments whose environment becomes the instance one.
# Blocks are the model author's own code, as the bodies of @mtkmodel are.
.eval_block <- function(x, env) {
  if (is.null(x)) return(NULL)
  if (is.language(x)) return(eval(x, env))
  if (is.function(x)) {
    environment(x) <- env
    return(x())
  }
  x
}

# Declarations of the DSL: strings such as "R = 1.0" or "x(t) = 0", or
# list("i(t)", connect = Flow, ...) with metadata. A factory argument with
# the variable's name overrides the default.
.dsl_declare <- function(specs, kind, env, args) {
  if (is.character(specs)) specs <- as.list(specs)
  out <- list()
  for (s in specs) {
    meta <- list()
    if (is.list(s)) {
      meta <- s[-1L]
      s <- s[[1L]]
    }
    d <- .parse_decl(s)
    default <- if (d$name %in% names(args)) args[[d$name]] else if (!is.null(d$default)) eval(d$default, env)
    if (!is.null(default) && !(is.numeric(default) || is_sym(default) || (length(default) == 1L && is.na(default))))
      default <- as_sym(default)
    meta$default <- default
    iv <- if (length(d$ivs)) lapply(d$ivs, function(n) {
      v <- get0(n, envir = env)
      if (is.null(v) || !is_sym(v)) .make_var(n, "iv", NULL, list()) else v
    })
    meta <- meta[!vapply(meta, is.null, TRUE)]
    v <- .make_array_or_scalar(d$name, kind, iv, meta, d$dims)
    assign(d$name, v, envir = env)
    out[[d$name]] <- v
  }
  out
}

# Build one component or base system from a quoted call, injecting its name
# and the factory arguments `comp__arg`.
.dsl_component <- function(cname, x, env, args) {
  if (inherits(x, "System")) return(.rename_system(x, cname))
  if (is.function(x) && !is.language(x)) {
    environment(x) <- env
    return(x())
  }
  if (!is.call(x)) stop("component ", cname, " must be a quoted call such as quote(Resistor(R = 1))", call. = FALSE)
  cl <- as.list(x)
  if (is.null(names(cl)) || !("name" %in% names(cl))) cl$name <- cname
  pre <- paste0(cname, "__")
  for (n in names(args)) if (startsWith(n, pre)) cl[[substring(n, nchar(pre) + 1L)]] <- args[[n]]
  eval(as.call(cl), env)
}

# Unpack the unknowns, parameters, observed variables and subsystems of a
# base system without namespacing, as @unpack does for @extend.
.dsl_unpack <- function(base, env) {
  for (v in c(.f(base, "unknowns"), .f(base, "ps"), lapply(.f(base, "observed"), `[[`, "lhs")))
    if (!is.null(.name(v))) assign(.name(v), v, envir = env)
  for (s in .f(base, "systems")) assign(nameof(s), s, envir = env)
}

#' Models from declarations
#'
#' `mtkmodel(name, ...)` returns a constructor `f(name, ...)` of systems,
#' the counterpart of ModelingToolkit's `@mtkmodel` blocks. Each call builds
#' the model in a fresh environment, in this order:
#' \enumerate{
#' \item `structural_parameters`: a named list of values;
#' \item `parameters`, `constants` and `variables`: declarations, strings
#'   such as `"R = 1.0"` or `"v(t) = 0"`, or `list("i(t)", connect = Flow)`
#'   to add metadata; each is bound to its name;
#' \item `extend`: a named list of quoted calls building base systems,
#'   whose variables and subsystems become available unnamespaced (as
#'   `@extend` and `@unpack`), and which the model extends;
#' \item `components`: a named list of quoted calls building subsystems,
#'   bound to their names; the name argument is added;
#' \item `equations`, `initial_conditions`, `guesses`, `continuous_events`
#'   and `discrete_events`: quoted expressions (or functions without
#'   arguments) evaluated in that environment, returning lists.
#' }
#' The constructor's arguments override structural parameters and the
#' defaults of declared parameters and variables; `comp__arg = value`
#' passes `arg = value` to the component `comp`; arguments named after the
#' declarations of an extended model are passed to it.
#'
#' @param name The model name.
#' @param structural_parameters Named list of structural parameters.
#' @param parameters,constants,variables Declarations.
#' @param components,extend Named lists of quoted constructor calls.
#' @param equations,initial_conditions,guesses,continuous_events,discrete_events
#'   Quoted expressions or functions without arguments.
#' @param description A string.
#' @param connector Build a connector.
#' @param iv The independent variable.
#' @param metadata A list of system metadata.
#' @param .env The environment in which blocks see other functions and
#'   variables.
#' @return A constructor function of class `Model`.
#' @examples
#' Decay <- mtkmodel("Decay",
#'   parameters = list("k = 1"),
#'   variables = list("x(t) = 1"),
#'   equations = quote(list(D_nounits(x) ~ -k * x)))
#' sys <- Decay(name = "decay", k = 2)
#' equations(sys)
#' @export
mtkmodel <- function(name, structural_parameters = list(), parameters = list(), constants = list(),
                     variables = list(), components = list(), extend = list(), equations = NULL,
                     initial_conditions = NULL, guesses = NULL, continuous_events = NULL,
                     discrete_events = NULL, description = "", connector = FALSE, iv = NULL,
                     metadata = list(), .env = parent.frame()) {
  model_name <- name
  decl_names <- function(specs) vapply(as.list(specs), function(s) .parse_decl(if (is.list(s)) s[[1L]] else s)$name, "")
  structure_info <- list(name = model_name, structural_parameters = names(structural_parameters),
                         parameters = decl_names(parameters), constants = decl_names(constants),
                         variables = decl_names(variables), components = names(components),
                         extend = names(extend), isconnector = connector)
  force(.env)
  # arguments of extended models are arguments of this one too (_arguments)
  inherited <- unlist(lapply(extend, function(cl) {
    ctor <- if (is.call(cl)) tryCatch(eval(cl[[1L]], .env), error = function(e) NULL)
    attr(ctor, "structure")$arguments
  }))
  structure_info$arguments <- unique(c(structure_info$structural_parameters, structure_info$parameters,
                                       structure_info$constants, structure_info$variables, inherited))
  own_args <- c(structure_info$structural_parameters, structure_info$parameters,
                structure_info$constants, structure_info$variables)
  f <- function(name, ...) {
    if (missing(name)) stop("a model needs a name: ", model_name, "(name = \"...\")", call. = FALSE)
    args <- list(...)
    env <- new.env(parent = .env)
    for (n in names(structural_parameters))
      assign(n, if (n %in% names(args)) args[[n]] else .eval_block(structural_parameters[[n]], env), envir = env)
    ps <- .dsl_declare(parameters, "parameter", env, args)
    cs <- .dsl_declare(constants, "constant", env, args)
    vs <- .dsl_declare(variables, "unknown", env, args)
    bases <- list()
    for (bn in names(extend)) {
      cl <- extend[[bn]]
      ctor <- if (is.call(cl)) tryCatch(eval(cl[[1L]], env), error = function(e) NULL)
      fwd <- attr(ctor, "structure")$arguments
      if (length(fwd) && is.call(cl)) {
        cl <- as.list(cl)
        for (n in intersect(fwd, names(args))) if (!(n %in% own_args)) cl[[n]] <- args[[n]]
        cl <- as.call(cl)
      }
      b <- .dsl_component(bn, cl, env, args)
      assign(bn, b, envir = env)
      .dsl_unpack(b, env)
      bases[[bn]] <- b
    }
    comps <- list()
    for (cn in names(components)) {
      s <- .dsl_component(cn, components[[cn]], env, args)
      assign(cn, s, envir = env)
      comps[[cn]] <- s
    }
    iv0 <- iv %||% t_nounits
    sys <- System(.eval_block(equations, env) %||% list(), iv0, unknowns = unname(vs),
                  ps = unname(c(ps, cs)), name = name, systems = unname(comps),
                  initial_conditions = .eval_block(initial_conditions, env),
                  guesses = .eval_block(guesses, env),
                  continuous_events = .eval_block(continuous_events, env) %||% list(),
                  discrete_events = .eval_block(discrete_events, env) %||% list(),
                  description = description,
                  metadata = utils::modifyList(as.list(metadata), list(gui_metadata = model_name)))
    for (b in bases) sys <- extend(sys, b)
    if (connector) sys <- .update_system(sys, connector_type = connector_type(sys))
    sys
  }
  structure(f, class = c("Model", "function"), structure = structure_info)
}

#' @export
print.Model <- function(x, ...) {
  s <- attr(x, "structure")
  cat(if (isTRUE(s$isconnector)) "Connector" else "Model", " ", s$name, "\n", sep = "")
  show <- function(lab, v) if (length(v)) cat("  ", lab, ": ", paste(v, collapse = ", "), "\n", sep = "")
  show("structural parameters", s$structural_parameters)
  show("parameters", c(s$parameters, s$constants))
  show("variables", s$variables)
  show("components", s$components)
  show("extends", s$extend)
  invisible(x)
}
