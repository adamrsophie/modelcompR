# Small public API of ModelingToolkitBase that belongs to no larger module:
# the @named, @nonamespace and @namespace macros (systems/abstractsystem.jl),
# apply_to_variables (abstractsystem.jl), collect_vars!, collect_scoped_vars!,
# collect_var_to_name!, eqtype_supports_collect_vars and
# convert_bindings_for_time_independent_system (utils.jl), Symbolics'
# get_variables, toexpr and substitute (substitute(sys, rules) from
# abstractsystem.jl), the metadata keys of variables.jl and parameters.jl,
# MiscSystemData, SymbolicADDisallowed and check_symbolic_ad_allowed
# (systems/system.jl), Both (systems/problem_utils.jl), get_w
# (systems/connectors.jl) and MTKVerbosity (verbosity.jl).

# ---------------------------------------------------------------------------
# Macros as functions

#' Name a system after the variable it is assigned to
#'
#' `named(sys <- Constructor(...))` (or `named(sys = Constructor(...))`) calls
#' `Constructor(..., name = "sys")` and assigns the result to `sys` in the
#' calling environment, as the `@named` macro of ModelingToolkit. Named
#' arguments whose values are symbolic expressions have their variables
#' wrapped in [ParentScope()], so a symbolic default passed into a component
#' refers to the variables of the enclosing system.
#'
#' @param ... One assignment `sys <- Constructor(...)` or one named argument
#'   `sys = Constructor(...)`.
#' @return Invisibly, the constructed object.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' decay <- function(name, k = 1) {
#'   x <- Variable("x", t)
#'   System(list(D(x) ~ -k * x), t, name = name)
#' }
#' named(sys <- decay(k = 2))
#' nameof(sys)
#' @export
named <- function(...) {
  # named(), nonamespace() and namespace() evaluate the caller's own code in
  # the caller's environment, as the macros they replace expand in place.
  env <- parent.frame()
  dots <- as.list(substitute(list(...)))[-1L]
  if (length(dots) != 1L) stop("named() takes one assignment: named(sys <- Constructor(...))", call. = FALSE)
  nm <- names(dots)
  ex <- dots[[1L]]
  if (!is.null(nm) && nzchar(nm)) {
    name <- nm
    cl <- ex
  } else if (is.call(ex) && (identical(ex[[1L]], as.name("<-")) || identical(ex[[1L]], as.name("=")))) {
    name <- as.character(ex[[2L]])
    cl <- ex[[3L]]
  } else {
    stop("use named(sys <- Constructor(...)) or named(sys = Constructor(...))", call. = FALSE)
  }
  if (!is.call(cl)) stop("named() needs a constructor call", call. = FALSE)
  argn <- names(cl) %||% rep("", length(cl))
  for (i in seq_along(cl)[-1L]) {
    if (!nzchar(argn[i]) || argn[i] == "name") next
    # evaluate once and pass the value, so side effects happen once
    v <- eval(cl[[i]], env)
    if (is_sym(v)) v <- apply_to_variables(ParentScope, v)
    cl[i] <- list(v)
  }
  cl$name <- name
  val <- eval(cl, env)
  assign(name, val, envir = env)
  invisible(val)
}

#' Variable access with or without namespacing
#'
#' `nonamespace(sys$sub$x)` evaluates the `$` accesses in its argument
#' without namespacing, so it returns the variable `x` as it is written in
#' the subsystem, as the `@nonamespace` macro. `namespace()` forces
#' namespacing (the `@namespace` macro), which matters for complete systems.
#'
#' @param expr An expression with `$` accesses on systems.
#' @return The value of `expr`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t)
#' inner <- System(list(D(x) ~ -x), t, name = "inner")
#' outer <- System(list(), t, name = "outer", systems = list(inner))
#' outer$inner$x
#' nonamespace(outer$inner$x)
#' @export
nonamespace <- function(expr) eval(.rewrite_ns(substitute(expr), FALSE), parent.frame())

#' @rdname nonamespace
#' @export
namespace <- function(expr) eval(.rewrite_ns(substitute(expr), TRUE), parent.frame())

.rewrite_ns <- function(ex, ns) {
  if (!is.call(ex)) return(ex)
  if (identical(ex[[1L]], as.name("$")))
    return(as.call(list(.getvar_ns, .rewrite_ns(ex[[2L]], ns), as.character(ex[[3L]]), ns)))
  for (i in seq_along(ex)[-1L]) {
    a <- ex[[i]]
    if (!missing(a) && !is.null(a)) ex[[i]] <- .rewrite_ns(a, ns)
  }
  ex
}

.getvar_ns <- function(obj, name, ns) {
  if (!inherits(obj, "System")) return(obj[[name]])
  # a complete system is flat; namespaced access goes through its hierarchy
  par <- .f(obj, "parent")
  if (ns && iscomplete(obj) && !is.null(par)) obj <- par
  .getproperty(toggle_namespacing(obj, ns), name)
}

# ---------------------------------------------------------------------------
# Walking expressions

#' Apply a function to the variables of an expression
#'
#' `apply_to_variables(f, ex)` replaces every variable in `ex` (anything
#' declared with [Variable()], [Parameter()], [independent_variable()] and
#' the like) by `f(variable)` and rebuilds the expression. Derivatives keep
#' their independent variable.
#'
#' @param f A function of one variable returning an expression.
#' @param ex An expression, an equation, or a list of them.
#' @return The transformed expression.
#' @examples
#' t <- t_nounits
#' x <- Variable("x", t); k <- Parameter("k")
#' apply_to_variables(ParentScope, k * x)
#' @export
apply_to_variables <- function(f, ex) {
  if (inherits(ex, "Equation")) return(Equation(apply_to_variables(f, ex$lhs), apply_to_variables(f, ex$rhs)))
  if (is.list(ex) && !is_sym(ex)) return(lapply(ex, function(e) apply_to_variables(f, e)))
  if (!is_sym(ex)) return(ex)
  if (.is_declared(ex) && !identical(.var_kind(ex), "operator")) return(f(ex))
  a <- .args(ex)
  if (!length(a)) return(ex)
  if (.type(ex) == "Derivative") {
    a[[1L]] <- apply_to_variables(f, a[[1L]])
    return(.rebuild_meta(ex, a))
  }
  .rebuild_meta(ex, lapply(a, function(z) apply_to_variables(f, z)))
}

# .rebuild() keeps `e` when the keys of the new arguments are unchanged, which
# drops metadata such as a scope set on a leaf; rebuild whenever the argument
# objects differ.
.rebuild_meta <- function(e, args) {
  if (identical(.args(e), args)) return(e)
  ty <- .type(e)
  if (ty == "Add") return(Reduce(`+`, args))
  if (ty == "Mul") return(Reduce(`*`, args))
  if (ty == "Pow") return(args[[1L]]^args[[2L]])
  if (ty == "Function") return(do.call(Function(.name(e)), args))
  if (ty == "Derivative") {
    di <- .deriv_info(e)
    return(.make_derivative(args[[1L]], di$iv, di$order))
  }
  .rebuild(e, args)
}

# Symbolic parts of an equation-like object, for collect_vars().
.sym_parts <- function(x) {
  if (is.null(x)) return(list())
  if (inherits(x, "formula")) x <- as_equation(x)
  if (inherits(x, c("Equation", "Inequality"))) return(Filter(is_sym, list(x$lhs, x$rhs)))
  if (is_sym(x)) return(list(x))
  if (is.list(x)) return(do.call(c, c(list(list()), lapply(unname(x), .sym_parts))))
  list()
}

.scope_of <- function(v, scopes = NULL) scopes[[.key(v)]] %||% attr(v, "mtk_meta")$scope

# check_scope_depth: local variables belong to depth 0, ParentScope(...) of
# depth d to depth d, and global variables to depth -1.
.check_scope_depth <- function(scope, depth) {
  if (is.null(scope) || identical(scope$type, "local")) return(depth == 0)
  if (identical(scope$type, "parent")) return(depth == scope$depth)
  if (identical(scope$type, "global")) return(depth == -1)
  FALSE
}

#' Collect the variables of expressions
#'
#' `collect_vars(expr, iv)` returns the unknowns and parameters appearing in
#' `expr` (an expression, equation, inequality, pair `list(a, b)`, or a list
#' of them), as `collect_vars!` in ModelingToolkit. Derivatives contribute
#' their variable; a delayed or evaluated variable `x(t - tau)` contributes
#' `x(t)` and the parameters of its argument. Only variables whose scope
#' matches `depth` are kept (0: local, d: `ParentScope` applied d times, -1:
#' global). `collect_scoped_vars(sys)` collects the variables of the
#' equations, jumps and constraints of `sys` that refer to the enclosing
#' system (depth 1 by default). `eqtype_supports_collect_vars(eq)` tells
#' whether an object is supported.
#'
#' @param expr Expressions or equations.
#' @param iv The independent variable, or `NULL`.
#' @param depth Scope depth to collect.
#' @param scopes Scope map of a system (internal).
#' @return A list with elements `unknowns` and `parameters`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); k <- Parameter("k")
#' collect_vars(D(x) ~ -k * x, t)
#' @export
collect_vars <- function(expr, iv = NULL, depth = 0, scopes = NULL) {
  found <- .collect_vars(.sym_parts(expr))
  extra <- list()
  if (!is.null(iv)) {
    for (v in found) {
      if (.type(v) == "Function" && identical(.var_kind(v), "unknown") && .is_delay_or_evalat(v, iv)) {
        u <- do.call(Function(.name(v)), list(iv))
        attr(u, "mtk_meta") <- attr(v, "mtk_meta")
        extra[[length(extra) + 1L]] <- u
      }
    }
  }
  cls <- .classify_vars(c(found, extra), iv)
  keep <- function(vs) Filter(function(v) .check_scope_depth(.scope_of(v, scopes), depth), vs)
  list(unknowns = keep(cls$unknowns), parameters = keep(cls$ps))
}

#' @rdname collect_vars
#' @param sys A `System`.
#' @export
collect_scoped_vars <- function(sys, iv = get_iv(sys), depth = 1) {
  eqs <- c(Filter(eqtype_supports_collect_vars, equations(sys)),
           Filter(eqtype_supports_collect_vars, jumps(sys)),
           Filter(eqtype_supports_collect_vars, constraints(sys)))
  collect_vars(eqs, iv, depth, scopes = .f(sys, "scopes"))
}

#' @rdname collect_vars
#' @param eq An object.
#' @export
eqtype_supports_collect_vars <- function(eq) {
  inherits(eq, c("Equation", "Inequality", "formula")) ||
    (is.list(eq) && !is_sym(eq) && length(eq) == 2L && all(vapply(eq, is_sym, TRUE)))
}

#' Map variable names to variables
#'
#' `collect_var_to_name(vars)` returns a named list from variable names to
#' variables, as `collect_var_to_name!`. Elements of array variables are
#' named by their array. Two different variables with the same name are an
#' error.
#'
#' @param vars A list of variables.
#' @param out An existing map to add to.
#' @return A named list.
#' @examples
#' t <- t_nounits
#' collect_var_to_name(list(Variable("x", t), Parameter("k")))
#' @export
collect_var_to_name <- function(vars, out = list()) {
  for (x in .flatten_vars(vars)) {
    if (.is_num(x) || is.null(.name(x))) next
    nm <- .meta_get(x, "array_parent") %||% .name(x)
    if (!is.null(.meta_get(x, "array_parent"))) {
      out[[nm]] <- out[[nm]] %||% x
      next
    }
    prev <- out[[nm]]
    if (!is.null(prev) && !identical(.key(prev), .key(x)))
      stop("Found variables with duplicate names! ", format(x), " and ", format(prev),
           " have the same name but are not identical.", call. = FALSE)
    out[[nm]] <- x
  }
  out
}

#' Variables of an expression
#'
#' `get_variables(expr)` returns the variables appearing in an expression or
#' equation, as Symbolics' `get_variables`: applied variables such as `x(t)`
#' count as one variable, derivatives contribute the differentiated variable,
#' and parameters and bare independent variables are included.
#'
#' @param expr An expression, equation or list of them.
#' @return A list of variables, unique, in order of appearance.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); k <- Parameter("k")
#' get_variables(D(x) ~ -k * x + t)
#' @export
get_variables <- function(expr) {
  out <- list()
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    if (.is_declared(z) && !identical(.var_kind(z), "operator")) {
      out[[length(out) + 1L]] <<- z
      return(invisible())
    }
    a <- .args(z)
    if (.type(z) == "Derivative") a <- a[1L]
    for (x in a) walk(x)
  }
  for (p in .sym_parts(expr)) walk(p)
  .unique_syms(out)
}

#' Convert symbolic expressions to R code
#'
#' `toexpr(x)` returns an R language object for an expression or equation,
#' as Symbolics' `toexpr`. Variables become names or calls such as `x(t)`,
#' derivatives `Differential(t)(x(t))`, and equations `lhs ~ rhs`.
#'
#' @param x An expression, equation, or list of them.
#' @return A language object, or a list of them.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); k <- Parameter("k")
#' toexpr(D(x) ~ -k * sin(x))
#' @export
toexpr <- function(x) {
  if (inherits(x, "formula")) x <- as_equation(x)
  if (inherits(x, "Equation")) return(call("~", toexpr(x$lhs), toexpr(x$rhs)))
  if (is.list(x) && !is_sym(x)) return(lapply(x, toexpr))
  if (!is_sym(x)) return(x)
  leaves <- list()
  walk <- function(z) {
    if (.is_declared(z) || .type(z) == "Derivative" ||
        (.type(z) == "Function" && !(.name(z) %in% c(names(.R_MATH), "log", "Max", "Min")))) {
      leaves[[.key(z)]] <<- z
      return(invisible())
    }
    if (.type(z) %in% c("Symbol", "Dummy")) leaves[[.key(z)]] <<- z
    for (a in .args(z)) walk(a)
  }
  walk(x)
  names_map <- list()
  ph <- list()
  for (i in seq_along(leaves)) {
    nm <- paste0(".toexpr", i)
    ph[[i]] <- Symbol(nm)
    names_map[[nm]] <- .leaf_rcall(leaves[[i]])
  }
  low <- if (length(leaves)) .subst(x, .rule(unname(leaves), ph)) else x
  .to_rcall(low, names_map)
}

.leaf_rcall <- function(z) {
  ty <- .type(z)
  if (ty %in% c("Symbol", "Dummy")) return(as.name(.name(z)))
  if (ty == "Derivative") {
    out <- toexpr(.args(z)[[1L]])
    for (spec in .args(z)[-1L]) {
      v <- .args(spec)[[1L]]
      o <- as.integer(as.numeric(.args(spec)[[2L]]))
      op <- if (o == 1L) call("Differential", toexpr(v)) else call("Differential", toexpr(v), o)
      out <- as.call(list(op, out))
    }
    return(out)
  }
  as.call(c(list(as.name(.name(z))), lapply(.args(z), toexpr)))
}

# ---------------------------------------------------------------------------
# Substitution

.as_rules <- function(rules, vars) {
  m <- .resolve_named(as_symmap(rules), vars)
  vals <- lapply(smap_values(m), function(v) if (is_sym(v)) v else as_sym(v))
  .rule(unname(smap_vars(m)), unname(vals))
}

#' Substitute expressions for variables
#'
#' `substitute_vars(x, rules)` replaces subexpressions of `x` as given by
#' `rules` (formulas `x ~ 2`, pairs `list(x, 2)`, or a named list), in one
#' top down pass. `x` may be an expression, an equation, a list, or a
#' `System`, in which case equations, unknowns, parameters, variable maps,
#' observed equations, costs, constraints and subsystems are substituted, as
#' `substitute(sys, rules)` in ModelingToolkit. With a named list of systems,
#' the subsystems of those names are replaced. This is ModelingToolkit's
#' `substitute`, renamed so that `base::substitute` is not masked.
#'
#' @param x An expression, equation, list or `System`.
#' @param rules The substitutions.
#' @return The substituted object.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); k <- Parameter("k")
#' substitute_vars(-k * x, list(k ~ 2))
#' @export
substitute_vars <- function(x, rules) {
  if (inherits(x, "System")) return(.substitute_system(x, rules))
  vars <- get_variables(if (is.list(x) && !is_sym(x)) unlist(lapply(x, .sym_parts), recursive = FALSE) else x)
  .subst(x, .as_rules(rules, vars))
}

.substitute_system <- function(sys, rules) {
  if (is.list(rules) && length(rules) && all(vapply(rules, inherits, TRUE, "System"))) {
    subs <- lapply(.f(sys, "systems"), function(s) {
      s <- .substitute_system(s, rules)
      rules[[nameof(s)]] %||% s
    })
    return(.update_system(sys, systems = subs))
  }
  if (length(.f(sys, "continuous_events")) || length(.f(sys, "discrete_events")))
    warning("substitute_vars() only substitutes in equations; the events of ", nameof(sys),
            " are not updated", call. = FALSE)
  vars <- c(.f(sys, "unknowns"), .f(sys, "ps"), independent_variables(sys))
  rule <- .as_rules(rules, vars)
  sub <- function(e) .subst(e, rule)
  submap <- function(m) {
    out <- symmap()
    for (el in unclass(m)) {
      if (is.null(el$var)) next
      out <- smap_set(out, sub(el$var), if (is_sym(el$value)) sub(el$value) else el$value)
    }
    out
  }
  noise <- .f(sys, "noise_eqs")
  if (!is.null(noise)) {
    noise <- if (is.matrix(noise)) {
      noise[] <- lapply(noise, sub)
      noise
    } else lapply(noise, sub)
  }
  .update_system(sys,
    eqs = lapply(.f(sys, "eqs"), function(e) if (inherits(e, "Equation")) sub(e) else e),
    unknowns = lapply(.f(sys, "unknowns"), sub), ps = lapply(.f(sys, "ps"), sub),
    bindings = submap(.f(sys, "bindings")), initial_conditions = submap(.f(sys, "initial_conditions")),
    guesses = submap(.f(sys, "guesses")), noise_eqs = noise, costs = lapply(.f(sys, "costs"), sub),
    observed = lapply(.f(sys, "observed"), sub),
    initialization_eqs = lapply(.f(sys, "initialization_eqs"), sub),
    constraints = lapply(.f(sys, "constraints"), function(c) {
      if (inherits(c, "Inequality")) {
        c$lhs <- sub(c$lhs)
        c$rhs <- sub(c$rhs)
        c
      } else sub(c)
    }),
    systems = lapply(.f(sys, "systems"), .substitute_system, rules = rules))
}

#' Build a symbolic term
#'
#' `Term(f, ...)` applies `f` to the arguments symbolically: a string names an
#' undefined (or registered) symbolic function, an R function such as `sin`
#' or `` `+` `` is called on the arguments. It plays the role of
#' SymbolicUtils' `Term`.
#'
#' @param f A function name or an R function.
#' @param ... Arguments.
#' @return A `sym` expression.
#' @examples
#' x <- Parameter("x")
#' Term("f", x, 2)
#' Term(sin, x)
#' @export
Term <- function(f, ...) {
  args <- lapply(list(...), as_sym)
  if (is.character(f)) return(do.call(Function(f), args))
  as_sym(do.call(f, args))
}

# ---------------------------------------------------------------------------
# Metadata keys. ModelingToolkit stores metadata under key types; here they
# are the names of the metadata fields read by getmetadata() and
# setmetadata() (variables) or stored in the `metadata` list of a System.

#' Metadata keys
#'
#' Names of the metadata fields of variables (for [getmetadata()] and
#' [setmetadata()]) and of systems (entries of the `metadata` list of
#' [System()]), standing for ModelingToolkit's metadata key types:
#' `VariableUnit` is `"unit"`, `VariableBounds` `"bounds"`, `SymScope`
#' `"scope"`, `MTKVariableTypeCtx` `"kind"`, and so on. `VariableType` maps
#' MTK's variable type enum (`VARIABLE`, `PARAMETER`, `BROWNIAN`,
#' `POISSONIAN`) to the kinds used by modelcompR. `ProblemTypeCtx`,
#' `MiscSystemData`, `SymbolicADDisallowed` and `MTKVerbosityCtx` are system
#' metadata keys.
#'
#' @format Character strings (`VariableType` a named character vector).
#' @examples
#' x <- Variable("x", t_nounits, unit = "m")
#' getmetadata(x, VariableUnit)
#' sys <- System(list(), t_nounits, name = "s",
#'               metadata = stats::setNames(list("my data"), MiscSystemData))
#' get_metadata(sys)[[MiscSystemData]]
#' @name metadata_keys
NULL

#' @rdname metadata_keys
#' @export
VariableBounds <- "bounds"
#' @rdname metadata_keys
#' @export
VariableConnectType <- "connect"
#' @rdname metadata_keys
#' @export
VariableDescription <- "description"
#' @rdname metadata_keys
#' @export
VariableInput <- "input"
#' @rdname metadata_keys
#' @export
VariableIrreducible <- "irreducible"
#' @rdname metadata_keys
#' @export
VariableMisc <- "misc"
#' @rdname metadata_keys
#' @export
VariableOutput <- "output"
#' @rdname metadata_keys
#' @export
VariableStatePriority <- "state_priority"
#' @rdname metadata_keys
#' @export
VariableUnit <- "unit"
#' @rdname metadata_keys
#' @export
MTKVariableTypeCtx <- "kind"
#' @rdname metadata_keys
#' @export
SymScope <- "scope"
#' @rdname metadata_keys
#' @export
VariableType <- c(VARIABLE = "unknown", PARAMETER = "parameter", BROWNIAN = "brownian",
                  POISSONIAN = "poissonian")
#' @rdname metadata_keys
#' @export
ProblemTypeCtx <- "problem_type"
#' @rdname metadata_keys
#' @export
MiscSystemData <- "misc_data"
#' @rdname metadata_keys
#' @export
SymbolicADDisallowed <- "symbolic_ad_disallowed"
#' @rdname metadata_keys
#' @export
MTKVerbosityCtx <- "verbosity"

#' Whether symbolic differentiation is allowed
#'
#' `check_symbolic_ad_allowed(sys)` errors when the system's metadata has an
#' entry [SymbolicADDisallowed], whose value (a string) explains why
#' symbolic Jacobians and similar derivatives are not available.
#'
#' @param sys A `System`.
#' @return Invisibly `TRUE`; errors otherwise.
#' @examples
#' sys <- System(list(), t_nounits, name = "s",
#'               metadata = list(symbolic_ad_disallowed = "it wraps a black box"))
#' try(check_symbolic_ad_allowed(sys))
#' @export
check_symbolic_ad_allowed <- function(sys) {
  reason <- (.f(sys, "metadata") %||% list())[[SymbolicADDisallowed]]
  if (is.null(reason) || identical(reason, FALSE)) return(invisible(TRUE))
  msg <- paste0("System ", nameof(sys), " does not support symbolic automatic differentiation.")
  if (is.character(reason) && nzchar(reason)) msg <- paste(msg, reason)
  stop(msg, call. = FALSE)
}

#' Sentinel for in-place or out-of-place problems
#'
#' In ModelingToolkit, `Both` is the default `iip` type parameter of problem
#' constructors, resolved at construction time. R functions always return
#' their result, so `Both` is kept only as a sentinel value.
#'
#' @format An object of class `Both`.
#' @examples
#' Both
#' @export
Both <- structure(list(), class = "Both")

#' @export
print.Both <- function(x, ...) {
  cat("Both\n")
  invisible(x)
}

#' Move bindings of variables to initial conditions
#'
#' For a time dependent system, `convert_bindings_for_time_independent_system()`
#' moves the bindings of unknowns and observed variables to initial
#' conditions, as needed to turn the system into a time independent one. The
#' system is not modified.
#'
#' @param sys A `System`.
#' @return A list with elements `initial_conditions` and `bindings`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' parameters("a = 2"); x <- Variable("x", t, default = 2 * a)
#' sys <- System(list(D(x) ~ -x), t, name = "s")
#' convert_bindings_for_time_independent_system(sys)$initial_conditions
#' @export
convert_bindings_for_time_independent_system <- function(sys) {
  dvs <- .keys(c(unknowns(sys), observables(sys)))
  ics <- initial_conditions(sys)
  binds <- bindings(sys)
  for (k in intersect(names(binds), dvs)) {
    ics[[k]] <- binds[[k]]
    binds[[k]] <- NULL
  }
  class(ics) <- class(binds) <- "mtk_symmap"
  list(initial_conditions = ics, bindings = binds)
}

#' Angular velocity of a rotation matrix
#'
#' `get_w(R, t)` returns the angular velocity `w` of a 3 by 3 rotation matrix
#' `R` of time dependent expressions, from `R` and its derivative `D(R)`:
#' `w = (R[3, ] . D(R)[2, ], -R[3, ] . D(R)[1, ], R[2, ] . D(R)[1, ])`.
#'
#' @param R A 3 by 3 matrix (a list matrix of expressions or a symbolic
#'   matrix).
#' @param t The independent variable.
#' @return A list of three expressions.
#' @examples
#' t <- t_nounits
#' th <- Variable("th", t)
#' R <- matrix(list(cos(th), sin(th), 0, -sin(th), cos(th), 0, 0, 0, 1), 3, 3)
#' get_w(R, t)[[3]]
#' @export
get_w <- function(R, t) {
  if (inherits(R, "symmatrix")) R <- .symmatrix_to_list_matrix(R)
  R <- matrix(lapply(R, as_sym), 3L, 3L)
  D <- Differential(t)
  DR <- matrix(lapply(R, function(e) expand_derivatives(D(e))), 3L, 3L)
  dot <- function(a, b) Reduce(`+`, lapply(1:3, function(j) a[[j]] * b[[j]]))
  list(dot(R[3L, ], DR[2L, ]), -dot(R[3L, ], DR[1L, ]), dot(R[2L, ], DR[1L, ]))
}

# ---------------------------------------------------------------------------
# Verbosity

.MTK_VERBOSITY_DEFAULTS <- c(
  state_priority_tie = "WarnLevel", underconstrained_variables = "Silent",
  if_lifting_condition_grammar = "WarnLevel", observed_equation_cycle = "Silent",
  singular_initialization = "WarnLevel", overdetermined_initialization = "WarnLevel",
  underdetermined_initialization = "WarnLevel", scc_initialization_unavailable = "WarnLevel",
  cyclic_dependency = "Silent", overdetermined_constraints = "WarnLevel",
  missing_scc_schedule = "WarnLevel", dynamic_opt_time_grid = "WarnLevel",
  empty_operating_point = "WarnLevel", initialization_analysis = "InfoLevel",
  no_unbound_inputs = "WarnLevel", analysis_point_causality = "WarnLevel",
  imperative_affect_specification = "WarnLevel", parameter_not_in_system = "WarnLevel",
  independent_variable_not_parameter = "WarnLevel",
  constraint_variable_without_argument = "WarnLevel", unbalanced_connector = "WarnLevel",
  affect_default_independent_variable = "WarnLevel", substitute_skips_events = "WarnLevel")

.MTK_LEVELS <- c("Silent", "DebugLevel", "InfoLevel", "WarnLevel", "ErrorLevel")

#' Verbosity of ModelingToolkit diagnostics
#'
#' `MTKVerbosity(...)` collects the levels (`"Silent"`, `"DebugLevel"`,
#' `"InfoLevel"`, `"WarnLevel"` or `"ErrorLevel"`) of ModelingToolkit's
#' diagnostic messages, one per toggle, as the `MTKVerbosity` specifier of
#' SciMLLogging. `preset` sets every toggle at once: `"Standard"` (the
#' defaults), `"None"` (all silent), `"Minimal"` (warnings and errors only),
#' `"Detailed"` and `"All"` (everything at least at info level). `TRUE` and
#' `FALSE` stand for `"Standard"` and `"None"`. Individual toggles given in
#' `...` override the preset.
#'
#' @param ... Toggles, such as `state_priority_tie = "Silent"`.
#' @param preset A preset name or a logical.
#' @return An object of class `MTKVerbosity`, a named character vector.
#' @examples
#' v <- MTKVerbosity(state_priority_tie = "Silent")
#' v[["dynamic_opt_time_grid"]]
#' MTKVerbosity(preset = FALSE)[["initialization_analysis"]]
#' @export
MTKVerbosity <- function(..., preset = "Standard") {
  if (isTRUE(preset)) preset <- "Standard"
  if (identical(preset, FALSE)) preset <- "None"
  lv <- .MTK_VERBOSITY_DEFAULTS
  lv[] <- switch(preset,
    Standard = lv,
    None = rep("Silent", length(lv)),
    Minimal = ifelse(lv %in% c("WarnLevel", "ErrorLevel"), lv, "Silent"),
    Detailed = , All = ifelse(lv %in% c("Silent", "DebugLevel"), "InfoLevel", lv),
    stop("unknown verbosity preset ", preset, call. = FALSE))
  over <- unlist(list(...))
  bad <- setdiff(names(over), names(lv))
  if (length(bad)) stop("unknown verbosity toggle(s): ", paste(bad, collapse = ", "), call. = FALSE)
  badlv <- setdiff(over, .MTK_LEVELS)
  if (length(badlv)) stop("unknown verbosity level(s): ", paste(badlv, collapse = ", "), call. = FALSE)
  lv[names(over)] <- over
  structure(lv, class = "MTKVerbosity")
}

#' @export
print.MTKVerbosity <- function(x, ...) {
  cat("MTKVerbosity\n")
  print(unclass(x), quote = FALSE)
  invisible(x)
}

# ---------------------------------------------------------------------------
# Hierarchical accessors of abstractsystem.jl not defined elsewhere

#' State selection hints of a system
#'
#' `state_priorities(sys)` returns the state priorities of `sys` and its
#' subsystems (a variable map, subsystem entries namespaced; the system's own
#' entries win), `irreducibles(sys)` the variables kept as unknowns by
#' [mtkcompile()], `maybe_zeros(sys)` the variables that simplification may
#' constrain to zero, and `analytically_integrated(sys)` the variables of a
#' flattened system that were integrated analytically.
#'
#' @param sys A `System`.
#' @return A variable map (`state_priorities`) or a list of variables.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t, state_priority = 10, irreducible = TRUE)
#' sys <- System(list(D(x) ~ -x), t, name = "s")
#' state_priorities(sys)
#' irreducibles(sys)
#' @export
state_priorities <- function(sys) .hier_map(sys, "state_priorities")

#' @rdname state_priorities
#' @export
irreducibles <- function(sys) .unique_syms(.collect_hier(sys, function(s) .f(s, "irreducibles")))

#' @rdname state_priorities
#' @export
maybe_zeros <- function(sys) .unique_syms(.collect_hier(sys, function(s) .f(s, "maybe_zeros")))

#' @rdname state_priorities
#' @export
analytically_integrated <- function(sys) {
  if (length(.f(sys, "systems")))
    stop("analytically integrated variables are only defined for flattened systems", call. = FALSE)
  .f(sys, "analytically_integrated") %||% list()
}

#' Abstract system type
#'
#' `AbstractSystem` holds the classes of modelcompR's systems, standing for
#' ModelingToolkit's abstract type: `inherits(x, AbstractSystem)` is `TRUE`
#' for a [System()] and a [PDESystem()].
#'
#' @format A character vector of class names.
#' @examples
#' inherits(System(list(), t_nounits, name = "s"), AbstractSystem)
#' @export
AbstractSystem <- c("System", "AbstractSystem")

#' LaTeX for equations and systems
#'
#' `latexify()` returns the LaTeX code of an equation, a list of equations or
#' the equations of a system, as the Latexify recipes of ModelingToolkit:
#' one equation is `lhs = rhs`, several are an `align` environment with one
#' `lhs &= rhs` row each. Other objects are passed to `symcompR::latex()`.
#'
#' @param x An `Equation`, a formula, a list of equations, a `System`, or a
#'   symbolic expression.
#' @param ... Unused.
#' @return A character string.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); k <- Parameter("k")
#' latexify(D(x) ~ -k * x)
#' cat(latexify(System(list(D(x) ~ -k * x), t, name = "s")))
#' @export
latexify <- function(x, ...) UseMethod("latexify")

#' @export
latexify.default <- function(x, ...) {
  if (is.list(x) && !is_sym(x)) {
    rows <- vapply(x, function(e) {
      e <- as_equation(e)
      paste0(symcompR::latex(e$lhs), " &= ", symcompR::latex(e$rhs))
    }, "")
    return(paste0("\\begin{align}\n", paste(rows, collapse = " \\\\\n"), "\n\\end{align}\n"))
  }
  symcompR::latex(x)
}

#' @export
latexify.Equation <- function(x, ...) paste(symcompR::latex(x$lhs), "=", symcompR::latex(x$rhs))

#' @export
latexify.formula <- function(x, ...) latexify(as_equation(x))

#' @export
latexify.System <- function(x, ...) latexify(equations(x))
