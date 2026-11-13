# Partial differential equation systems, after
# ModelingToolkitBase/src/systems/pde/pdesystem.jl (struct PDESystem, its
# getproperty, accessor bridges and show method). As in ModelingToolkit the
# container is symbolic only: discretization is left to other packages.
# Domains pair independent variables with intervals; symcompR's Interval()
# is the interval type (DomainSets' Interval in Julia).

#' Pair variables with their domains
#'
#' `VarDomainPairing(t, Interval(0, 1))` states that `t` ranges over the
#' interval, as \ifelse{latex}{\verb{t in Interval(0, 1)}}{\verb{t ∈ Interval(0, 1)}} in ModelingToolkit. A numeric vector
#' `c(0, 1)` is converted to symcompR's [symcompR::Interval()].
#'
#' @param variables A variable or a list of variables.
#' @param domain An interval, or a numeric vector `c(lower, upper)`.
#' @return An object of class `VarDomainPairing`.
#' @examples
#' x <- Parameter("x")
#' VarDomainPairing(x, c(0, 1))
#' @export
VarDomainPairing <- function(variables, domain) {
  if (is.numeric(domain) || (is.list(domain) && !is_sym(domain))) {
    if (length(domain) != 2L) stop("a domain is an interval c(lower, upper)", call. = FALSE)
    domain <- Interval(domain[[1L]], domain[[2L]])
  }
  structure(list(variables = variables, domain = domain), class = "VarDomainPairing")
}

#' @export
format.VarDomainPairing <- function(x, ...) {
  v <- if (is_sym(x$variables)) format(x$variables) else
    paste0("(", paste(vapply(x$variables, format, ""), collapse = ", "), ")")
  paste(v, "\u2208", format(x$domain))
}

#' @export
print.VarDomainPairing <- function(x, ...) {
  cat(format.VarDomainPairing(x), "\n")
  invisible(x)
}

.as_domain <- function(d) {
  if (inherits(d, "VarDomainPairing")) return(d)
  if (inherits(d, "formula")) {
    # both sides are the caller's own code, evaluated where the formula was
    # written, as as_equation() does
    env <- environment(d) %||% parent.frame()
    return(VarDomainPairing(eval(d[[2L]], env), eval(d[[3L]], env)))
  }
  stop("a domain is VarDomainPairing(var, interval) or a formula var ~ c(lower, upper)", call. = FALSE)
}

#' Systems of partial differential equations
#'
#' `PDESystem(eqs, bcs, domain, ivs, dvs, ps, name = )` collects the
#' equations, boundary and initial conditions, domains, independent and
#' dependent variables and parameters of a PDE problem. ModelingToolkit does
#' not discretize PDE systems itself and neither does modelcompR; the object
#' is the input of discretization packages.
#'
#' `analytic` gives exact solutions `u(t, x) ~ expr`, from which
#' `analytic_func` functions are generated, called with the parameter values
#' first and then the independent variables, as `f(ps, t, x)`.
#'
#' Fields are read with `$`: `pde$eqs`, `pde$bcs`, `pde$domain`, `pde$ivs`,
#' `pde$dvs`, `pde$ps`, ...; `pde$x` returns the variable named `x`.
#'
#' @param eqs Equations (formulas or `Equation` objects).
#' @param bcs Boundary and initial conditions.
#' @param domain Domains: a list of [VarDomainPairing()] objects or formulas
#'   `t ~ c(0, 1)`.
#' @param ivs Independent variables.
#' @param dvs Dependent variables, as applied functions `u(t, x)`.
#' @param ps Parameters.
#' @param initial_conditions Default values (variable map).
#' @param systems Subsystems.
#' @param connector_type Connector type.
#' @param metadata Metadata.
#' @param analytic Analytic solutions, equations `u(t, x) ~ expr`.
#' @param analytic_func Analytic solution functions, a named list.
#' @param checks Check units, as in [System()].
#' @param description A description.
#' @param name The name.
#' @return An object of class `PDESystem`.
#' @examples
#' t <- independent_variable("t"); x <- independent_variable("x")
#' u <- Variable("u", list(t, x))
#' Dt <- Differential(t); Dxx <- Differential(x, 2)
#' eq <- Dt(u) ~ Dxx(u)
#' bcs <- list(Function("u")(0, x) ~ sin(pi * x))
#' pde <- PDESystem(eq, bcs, list(t ~ c(0, 1), x ~ c(0, 1)), list(t, x), list(u),
#'                  name = "heat")
#' pde
#' @export
PDESystem <- function(eqs, bcs, domain, ivs, dvs, ps = list(), initial_conditions = NULL,
                      systems = list(), connector_type = NULL, metadata = NULL, analytic = NULL,
                      analytic_func = NULL, checks = TRUE, description = "", name = NULL) {
  if (is.null(name)) stop("a system needs a name: PDESystem(..., name = \"pde\")", call. = FALSE)
  eqs <- .as_equations(eqs)
  bcs <- .as_equations(bcs)
  if (inherits(domain, c("VarDomainPairing", "formula"))) domain <- list(domain)
  domain <- lapply(domain, .as_domain)
  ivs <- .flatten_vars(ivs)
  dvs <- .flatten_vars(dvs)
  ps <- .flatten_vars(ps)
  if (.checks_on(checks, 4L) && .has_units(c(dvs, ivs, ps))) {
    if (!validate(eqs)) .validation_error("Some equations had invalid units. See warnings for details.")
  }
  if (!is.null(analytic)) {
    analytic <- .as_equations(analytic)
    if (length(analytic) != length(dvs))
      stop("The number of analytic solutions must match the number of dependent variables", call. = FALSE)
    if (is.null(analytic_func)) {
      analytic_func <- lapply(analytic, function(eq) {
        args <- c(list(list(name = "ps", vars = ps)),
                  lapply(.args(eq$lhs), function(a) list(name = .name(a), vars = a)))
        .build_closure(eq$rhs, args, out = "scalar")
      })
      names(analytic_func) <- vapply(analytic, function(eq) format(eq$lhs), "")
    }
  }
  structure(list(eqs = eqs, bcs = bcs, domain = domain, ivs = ivs, dvs = dvs, ps = ps,
                 initial_conditions = as_symmap(initial_conditions), connector_type = connector_type,
                 systems = as.list(systems), analytic = analytic, analytic_func = analytic_func,
                 name = as.character(name), description = description, metadata = metadata,
                 unknowns = dvs, observed = list()),
            class = c("PDESystem", "AbstractSystem"))
}

#' @export
`$.PDESystem` <- function(x, name) {
  if (name == "indvars") return(.subset2(x, "ivs"))
  if (name == "depvars") return(.subset2(x, "dvs"))
  if (name %in% names(unclass(x))) return(.subset2(x, name))
  for (v in c(.subset2(x, "ivs"), .subset2(x, "dvs"), .subset2(x, "ps")))
    if (identical(.name(v), name)) return(v)
  for (s in .subset2(x, "systems")) if (identical(nameof(s), name)) return(s)
  stop("variable ", name, " does not exist in system ", .subset2(x, "name"), call. = FALSE)
}

#' @export
names.PDESystem <- function(x) {
  c(vapply(.subset2(x, "systems"), nameof, ""),
    vapply(c(.subset2(x, "dvs"), .subset2(x, "ps"), .subset2(x, "ivs")), function(v) .name(v), ""))
}

#' @export
print.PDESystem <- function(x, ...) {
  cat("PDESystem ", .subset2(x, "name"), "\n", sep = "")
  show <- function(lbl, xs) cat(lbl, paste(vapply(xs, .fmt, ""), collapse = ", "), "\n")
  show("Equations:", .subset2(x, "eqs"))
  show("Boundary Conditions:", .subset2(x, "bcs"))
  cat("Domain:", paste(vapply(.subset2(x, "domain"), format.VarDomainPairing, ""), collapse = ", "), "\n")
  show("Dependent Variables:", .subset2(x, "dvs"))
  show("Independent Variables:", .subset2(x, "ivs"))
  show("Parameters:", .subset2(x, "ps"))
  ic <- .subset2(x, "initial_conditions")
  cat("Default Parameter Values:", if (length(ic)) paste(vapply(unclass(ic), function(el)
    paste(.fmt(el$var), "=>", .fmt(el$value)), ""), collapse = ", ") else "none", "\n")
  invisible(x)
}
