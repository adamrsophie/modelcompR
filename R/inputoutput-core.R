# Inputs, outputs and control functions, after
# ModelingToolkitBase/src/inputoutput.jl (bound_inputs, unbound_inputs,
# bound_outputs, unbound_outputs, is_bound, IsBoundValidator,
# same_or_inner_namespace, generate_control_function), systems/codegen.jl
# (calculate_control_jacobian, generate_control_jacobian), variables.jl
# (disturbances) and the analysis point method of generate_control_function
# in systems/analysis_points.jl.

#' Bound and unbound inputs and outputs
#'
#' An input or output variable is bound when it appears in an equation
#' together with a variable of another subsystem (or through observed
#' equations with a variable bound further out): it is then internal to
#' `sys`. `unbound_inputs()` are the external inputs a user has to supply.
#'
#' @param sys A `System`.
#' @param u An input or output variable.
#' @return A list of variables, or a logical for `is_bound()`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' u <- Variable("u", t, input = TRUE); x <- Variable("x", t)
#' inner <- System(list(D(x) ~ -x + u), t, name = "inner")
#' r <- Variable("r", t)
#' outer <- System(list(inner$u ~ r, D(r) ~ 0), t, name = "outer", systems = list(inner))
#' is_bound(outer, outer$inner$u)
#' unbound_inputs(inner)
#' @export
is_bound <- function(sys, u) {
  vv <- function(e) {
    vs <- .collect_vars(list(e$lhs, e$rhs))
    .unique_syms(lapply(vs, function(v) if (.is_derivative(v)) .args(v)[[1L]] else v))
  }
  eqs_vars <- lapply(Filter(function(e) inherits(e, "Equation"), equations(sys)), vv)
  obs_vars <- lapply(observed(sys), vv)
  .is_bound_rec(eqs_vars, obs_vars, as_sym(u), character())
}

.is_bound_rec <- function(eqs_vars, obs_vars, u, stack) {
  uk <- .key(u)
  if (uk %in% stack) return(FALSE)
  for (vars in eqs_vars) {
    if (!(uk %in% .keys(vars))) next
    for (v in vars) {
      if (.key(v) == uk) next
      if (!.same_or_inner_namespace(u, v)) return(TRUE)
    }
  }
  for (vars in obs_vars) {
    if (!(uk %in% .keys(vars))) next
    for (v in vars) {
      if (.key(v) == uk) next
      if (!.same_or_inner_namespace(u, v)) return(TRUE)
      if (.is_bound_rec(eqs_vars, obs_vars, v, c(stack, uk)) && !.inner_namespace(u, v)) return(TRUE)
    }
  }
  FALSE
}

.get_namespace <- function(x) {
  parts <- strsplit(.name(x) %||% "", NAMESPACE_SEPARATOR, fixed = TRUE)[[1L]]
  if (length(parts) <= 1L) return("")
  paste(parts[-length(parts)], collapse = NAMESPACE_SEPARATOR)
}

.is_namespaced <- function(x) grepl(NAMESPACE_SEPARATOR, .name(x) %||% "", fixed = TRUE)

.same_or_inner_namespace <- function(u, v) {
  nu <- .get_namespace(u); nv <- .get_namespace(v)
  nu == nv || startsWith(nv, nu) || (.is_namespaced(v) && !.is_namespaced(u))
}

.inner_namespace <- function(u, v) {
  nu <- .get_namespace(u); nv <- .get_namespace(v)
  if (nu == nv) return(FALSE)
  startsWith(nv, nu) || (.is_namespaced(v) && !.is_namespaced(u))
}

#' @rdname is_bound
#' @export
bound_inputs <- function(sys) Filter(function(x) is_bound(sys, x), inputs(sys))
#' @rdname is_bound
#' @export
unbound_inputs <- function(sys) Filter(function(x) !is_bound(sys, x), inputs(sys))
#' @rdname is_bound
#' @export
bound_outputs <- function(sys) Filter(function(x) is_bound(sys, x), outputs(sys))
#' @rdname is_bound
#' @export
unbound_outputs <- function(sys) Filter(function(x) !is_bound(sys, x), outputs(sys))

# disturbances: unknowns and parameters marked as disturbances.
.disturbances <- function(sys) c(Filter(isdisturbance, unknowns(sys)), Filter(isdisturbance, .parameters_of(sys)))

# default_codegen_inputs: declared inputs of compiled systems, external
# inputs otherwise.
.default_codegen_inputs <- function(sys) if (isTRUE(.f(sys, "isscheduled"))) inputs(sys) else unbound_inputs(sys)

#' Control functions of input output systems
#'
#' `generate_control_function(sys, inputs)` compiles `sys` with `inputs` as
#' inputs and returns its dynamics as an R function `f(x, u, p, t)` (or
#' `f(x, u, p, t, w)` with known disturbances), where `x` are the unknowns
#' `dvs`, `u` the inputs and `p` the parameters `ps`. Unknown
#' `disturbance_inputs` keep their dynamics but are set to zero. With
#' `implicit_dae = TRUE` the function is the residual `f(dx, x, u, p, t)`.
#' `inputs` may also be analysis points, whose connections are opened.
#'
#' `calculate_control_jacobian(sys)` is the Jacobian of the right-hand sides
#' of `sys` with respect to its inputs and `generate_control_jacobian(sys)`
#' compiles it to a function of `(u, p, t)`.
#'
#' @param sys A `System`.
#' @param inputs Input variables or analysis points.
#' @param disturbance_inputs Unknown disturbance inputs.
#' @param known_disturbance_inputs Disturbance inputs passed as `w`.
#' @param implicit_dae Generate a residual function.
#' @param simplify Passed to [mtkcompile()].
#' @param system_modifier For analysis point inputs, a function applied to
#'   the opened system.
#' @param sparse Unused.
#' @param ... Unused.
#' @return A list with `f`, `dvs`, `ps` and the compiled system `io_sys`;
#'   for the Jacobian functions a symbolic matrix or an R function.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); u <- Variable("u", t, input = TRUE)
#' k <- Parameter("k")
#' sys <- System(list(D(x) ~ -k * (x + u)), t, name = "s")
#' cf <- generate_control_function(sys, list(u))
#' cf$f(1, 3, 2, 0)
#' @export
generate_control_function <- function(sys, inputs = NULL, disturbance_inputs = NULL,
                                      known_disturbance_inputs = NULL, implicit_dae = FALSE,
                                      simplify = FALSE, system_modifier = identity, ...) {
  if (!is.null(inputs) && .is_ap_like(inputs)) {
    u <- list()
    for (ap in .canonicalize_ap(sys, inputs)) {
      ol <- open_loop(sys, ap)
      sys <- ol$sys
      u <- c(u, ol$du)
    }
    d <- NULL
    if (!is.null(disturbance_inputs) && .is_ap_like(disturbance_inputs)) {
      d <- list()
      for (ap in .canonicalize_ap(sys, disturbance_inputs)) {
        ol <- open_loop(sys, ap)
        sys <- ol$sys
        d <- c(d, ol$du)
      }
    }
    return(generate_control_function(system_modifier(sys), u, d,
                                     known_disturbance_inputs = known_disturbance_inputs,
                                     implicit_dae = implicit_dae, simplify = simplify))
  }
  inputs <- .flatten_vars(if (is.null(inputs)) .default_codegen_inputs(sys) else inputs)
  if (is.null(disturbance_inputs)) disturbance_inputs <- .disturbances(sys)
  disturbance_inputs <- .flatten_vars(disturbance_inputs)
  known <- .flatten_vars(known_disturbance_inputs)
  if (!length(inputs)) message("No unbound inputs were found in system.")
  all_dist <- c(disturbance_inputs, known)
  if (!isTRUE(.f(sys, "isscheduled")))
    sys <- mtkcompile(sys, inputs = inputs, disturbance_inputs = all_dist, simplify = simplify)
  dvs <- unknowns(sys)
  io_keys <- .keys(c(inputs, all_dist))
  ps <- Filter(function(q) !(.key(q) %in% io_keys), .problem_params(sys))
  rhs <- .ode_rhs_exprs(sys, substitute_observed = TRUE)
  if (length(disturbance_inputs))
    rhs <- lapply(rhs, .subst, rule = .rule(disturbance_inputs, rep(list(.S0()), length(disturbance_inputs))))
  iv <- .f(sys, "iv")
  args <- list(list(name = "x", vars = dvs), list(name = "u", vars = inputs),
               list(name = "p", vars = ps), list(name = "t", vars = iv))
  if (length(known)) args[[5L]] <- list(name = "w", vars = known)
  if (implicit_dae) {
    eqs <- equations(sys)
    dx <- lapply(dvs, function(v) .make_derivative(v, iv, 1L))
    rhs <- lapply(seq_along(rhs), function(i) if (.is_derivative(eqs[[i]]$lhs)) rhs[[i]] - eqs[[i]]$lhs else rhs[[i]])
    args <- c(list(list(name = "dx", vars = dx)), args)
  }
  f <- .build_closure(rhs, args, out = "vector")
  list(f = f, dvs = dvs, ps = ps, io_sys = sys)
}

#' @rdname generate_control_function
#' @export
calculate_control_jacobian <- function(sys, sparse = FALSE, simplify = FALSE) {
  rhs <- lapply(full_equations(sys), `[[`, "rhs")
  ctrls <- .default_codegen_inputs(sys)
  others <- Filter(function(q) !(.key(q) %in% .keys(ctrls)),
                   c(unknowns(sys), .parameters_of(sys, TRUE), independent_variables(sys)))
  J <- .partials(rhs, ctrls, others)
  if (simplify) J[] <- lapply(J, symcompR::simplify)
  .as_symmatrix(J)
}

#' @rdname generate_control_function
#' @export
generate_control_jacobian <- function(sys, simplify = FALSE, sparse = FALSE, ...) {
  rhs <- lapply(full_equations(sys), `[[`, "rhs")
  ctrls <- .default_codegen_inputs(sys)
  others <- Filter(function(q) !(.key(q) %in% .keys(ctrls)),
                   c(unknowns(sys), .parameters_of(sys, TRUE), independent_variables(sys)))
  J <- .partials(rhs, ctrls, others)
  if (simplify) J[] <- lapply(J, symcompR::simplify)
  .lin_mat_fn(J, .codegen_args(sys), .codegen_observed(sys))
}
