# Optimization problems built from systems, after
# ModelingToolkitBase/src/problems/optimizationproblem.jl
# (OptimizationFunction, MultiObjectiveOptimizationFunction,
# OptimizationProblem, generate_constraint_fields, system_with_cost_weights,
# constraints_to_penalties), systems/codegen.jl (generate_cost,
# calculate_cost_gradient, generate_cost_gradient, calculate_cost_hessian,
# generate_cost_hessian, constraint_residual, canonical_constraints,
# generate_cons, calculate_constraint_jacobian, generate_constraint_jacobian,
# calculate_constraint_hessian, generate_constraint_hessian and the
# multiobjective variants), variables.jl (tunable_parameters), parameters.jl
# (subset_tunables) and systems/index_cache.jl
# (reorder_dimension_by_tunables).
#
# Generated functions are f(u, p) with u the unknowns (the optimization
# variables) and p the problem parameter vector. Derivatives are symbolic.

.opt_args <- function(sys) .codegen_args(sys, with_t = FALSE)

# Substitute observed variables and bound parameters, for differentiation.
.opt_full <- function(sys, exprs) {
  obs <- .codegen_observed(sys)
  if (!length(obs)) return(exprs)
  rule <- .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs"))
  lapply(exprs, .subst_fixpoint, rule = rule)
}

.opt_hold <- function(sys) .parameters_of(sys, TRUE)

.opt_gradient <- function(sys, expr) {
  us <- unknowns(sys)
  g <- .partials(.opt_full(sys, list(expr)), us, .opt_hold(sys))
  lapply(seq_along(us), function(j) g[[1L, j]])
}

.opt_hessian <- function(sys, expr) {
  us <- unknowns(sys)
  .partials(.opt_gradient(sys, expr), us, .opt_hold(sys))
}

.opt_simplify <- function(x, simplify) {
  if (!isTRUE(simplify)) return(x)
  x[] <- lapply(x, symcompR::simplify)
  x
}

#' Residuals of constraints
#'
#' `constraint_residual(cstr)` is the expression `r` such that an equation
#' constraint reads `0 ~ r` and an inequality constraint `0 >= r`:
#' `lhs - rhs`, with the sides exchanged for `>=`.
#' `canonical_constraints(sys)` returns the residuals of all constraints of a
#' system, one per row of the constraint function, in the order of
#' `constraints(sys)`.
#'
#' @param cstr An `Equation` or `Inequality`.
#' @param sys A `System`.
#' @param scalarize Unused; constraints are always scalar here.
#' @return A `sym` expression, or a list of them.
#' @examples
#' x <- Variable("x"); y <- Variable("y")
#' constraint_residual(Inequality(x^2 + y^2, 1, ">="))
#' sys <- complete(System(list(), costs = list(x^2 + y^2), name = "opt",
#'                        constraints = list(x + y ~ 1, Inequality(x, 0.8))))
#' canonical_constraints(sys)
#' @export
constraint_residual <- function(cstr) {
  if (inherits(cstr, "formula")) cstr <- as_equation(cstr)
  if (inherits(cstr, "Inequality")) return(.ineq_residual(cstr))
  if (inherits(cstr, "Equation")) return(cstr$lhs - cstr$rhs)
  stop("a constraint is an Equation or an Inequality", call. = FALSE)
}

#' @rdname constraint_residual
#' @export
canonical_constraints <- function(sys, scalarize = TRUE) lapply(constraints(sys), constraint_residual)

#' Cost functions and their derivatives
#'
#' For a system with costs (an optimization system): `generate_cost(sys)`
#' returns the R function `f(u, p)` of the consolidated cost;
#' `calculate_cost_gradient(sys)` and `calculate_cost_hessian(sys)` return
#' its symbolic gradient and Hessian with respect to the unknowns, and
#' `generate_cost_gradient(sys)`, `generate_cost_hessian(sys)` the
#' corresponding R functions. `generate_cons(sys)` returns the constraint
#' function (the residuals of [canonical_constraints()]),
#' `calculate_constraint_jacobian(sys)` its Jacobian and
#' `calculate_constraint_hessian(sys)` the list of Hessians of its rows, with
#' `generate_*` counterparts.
#'
#' The multiobjective variants keep the costs separate:
#' `generate_multiobjective_cost(sys)` returns the vector of `costs(sys)`,
#' `calculate_multiobjective_jacobian(sys)` its Jacobian and
#' `calculate_multiobjective_hessian(sys)` one Hessian per cost.
#'
#' @param sys A complete `System` with costs.
#' @param simplify Simplify the expressions.
#' @param sparse Unused; matrices are dense.
#' @param ... Passed to the code generator.
#' @return An R function of `(u, p)`, a list of expressions (gradient), a
#'   symbolic matrix (Hessian, Jacobian), or a list of them.
#' @examples
#' x <- Variable("x"); y <- Variable("y"); a <- Parameter("a")
#' sys <- complete(System(list(), costs = list((a - x)^2 + 100 * (y - x^2)^2), name = "rosenbrock"))
#' f <- generate_cost(sys)
#' f(c(0, 0), 1)
#' calculate_cost_gradient(sys)
#' generate_cost_gradient(sys)(c(0, 0), 1)
#' generate_cost_hessian(sys)(c(1, 1), 1)
#' @export
generate_cost <- function(sys, ...) {
  .build_closure(cost(sys), .opt_args(sys), observed = .codegen_observed(sys), out = "scalar", ...)
}

#' @rdname generate_cost
#' @export
calculate_cost_gradient <- function(sys, simplify = FALSE) {
  g <- .opt_gradient(sys, cost(sys))
  if (isTRUE(simplify)) g <- lapply(g, symcompR::simplify)
  g
}

#' @rdname generate_cost
#' @export
generate_cost_gradient <- function(sys, simplify = FALSE, ...) {
  .build_closure(calculate_cost_gradient(sys, simplify), .opt_args(sys), out = "vector", ...)
}

#' @rdname generate_cost
#' @export
calculate_cost_hessian <- function(sys, sparse = FALSE, simplify = FALSE) {
  .as_symmatrix(.opt_simplify(.opt_hessian(sys, cost(sys)), simplify))
}

#' @rdname generate_cost
#' @export
generate_cost_hessian <- function(sys, simplify = FALSE, sparse = FALSE, ...) {
  .build_closure(.opt_simplify(.opt_hessian(sys, cost(sys)), simplify), .opt_args(sys), out = "matrix", ...)
}

#' @rdname generate_cost
#' @export
generate_cons <- function(sys, ...) {
  .build_closure(canonical_constraints(sys), .opt_args(sys), observed = .codegen_observed(sys),
                 out = "vector", ...)
}

.opt_cons_jacobian <- function(sys) {
  .partials(.opt_full(sys, canonical_constraints(sys)), unknowns(sys), .opt_hold(sys))
}

#' @rdname generate_cost
#' @export
calculate_constraint_jacobian <- function(sys, simplify = FALSE, sparse = FALSE) {
  .as_symmatrix(.opt_simplify(.opt_cons_jacobian(sys), simplify))
}

#' @rdname generate_cost
#' @export
generate_constraint_jacobian <- function(sys, simplify = FALSE, sparse = FALSE, ...) {
  .build_closure(.opt_simplify(.opt_cons_jacobian(sys), simplify), .opt_args(sys), out = "matrix", ...)
}

#' @rdname generate_cost
#' @export
calculate_constraint_hessian <- function(sys, simplify = FALSE, sparse = FALSE) {
  lapply(canonical_constraints(sys), function(r) .as_symmatrix(.opt_simplify(.opt_hessian(sys, r), simplify)))
}

#' @rdname generate_cost
#' @export
generate_constraint_hessian <- function(sys, simplify = FALSE, sparse = FALSE, ...) {
  fns <- lapply(canonical_constraints(sys), function(r)
    .build_closure(.opt_simplify(.opt_hessian(sys, r), simplify), .opt_args(sys), out = "matrix", ...))
  function(u, p) lapply(fns, function(f) f(u, p))
}

#' @rdname generate_cost
#' @export
generate_multiobjective_cost <- function(sys, ...) {
  .build_closure(costs(sys), .opt_args(sys), observed = .codegen_observed(sys), out = "vector", ...)
}

#' @rdname generate_cost
#' @export
calculate_multiobjective_jacobian <- function(sys, simplify = FALSE, sparse = FALSE) {
  J <- .partials(.opt_full(sys, costs(sys)), unknowns(sys), .opt_hold(sys))
  .as_symmatrix(.opt_simplify(J, simplify))
}

#' @rdname generate_cost
#' @export
generate_multiobjective_jacobian <- function(sys, simplify = FALSE, sparse = FALSE, ...) {
  J <- .partials(.opt_full(sys, costs(sys)), unknowns(sys), .opt_hold(sys))
  .build_closure(.opt_simplify(J, simplify), .opt_args(sys), out = "matrix", ...)
}

#' @rdname generate_cost
#' @export
calculate_multiobjective_hessian <- function(sys, simplify = FALSE, sparse = FALSE) {
  lapply(costs(sys), function(obj) .as_symmatrix(.opt_simplify(.opt_hessian(sys, obj), simplify)))
}

#' @rdname generate_cost
#' @export
generate_multiobjective_hessian <- function(sys, simplify = FALSE, sparse = FALSE, ...) {
  lapply(costs(sys), function(obj)
    .build_closure(.opt_simplify(.opt_hessian(sys, obj), simplify), .opt_args(sys), out = "matrix", ...))
}

# ---------------------------------------------------------------------------
# Weights and penalties

# system_with_cost_weights: the consolidated cost becomes sum(weights *
# costs). Symbolic weights must be parameters of the system.
.opt_with_weights <- function(sys, weights) {
  cs <- .f(sys, "costs")
  if (!is.list(weights)) weights <- as.list(weights)
  if (length(weights) != length(cs))
    .stopf("Expected `weights` to have one entry per cost of `sys`, but got %d weights for %d costs.",
           length(weights), length(cs))
  pk <- .keys(.parameters_of(sys, TRUE))
  for (w in weights) {
    if (is.numeric(w)) next
    if (!is_sym(w)) stop("Entries of `weights` must be real numbers or symbolic parameters of `sys`.", call. = FALSE)
    if (!.is_num(w) && !(.key(w) %in% pk))
      stop("Symbolic weight `", format(w), "` is not a parameter of `sys`. Declare it with `Parameter()` in ",
           "the system so that it can be provided and updated through the parameter vector.", call. = FALSE)
  }
  .update_system(sys, consolidate = function(costs) Reduce(`+`, Map(function(w, cc) as_sym(w) * cc, weights, costs)))
}

#' Turn constraints into penalty costs
#'
#' `constraints_to_penalties(sys, weights = 1)` returns a system in which
#' every constraint of `sys`, including those of its subsystems, is removed
#' and added to the costs as a weighted quadratic penalty: an equation
#' `l ~ r` contributes `w * (l - r)^2` and an inequality with residual `res`
#' (`res <= 0`) contributes `w * max(res, 0)^2`. The result can be solved by optimizers
#' that do not handle constraints. `weights` is one number or expression for
#' all constraints, or one per constraint; symbolic weights that are not yet
#' parameters of the system are added to its parameters.
#'
#' @param sys A `System` with constraints.
#' @param weights Penalty weights.
#' @return A `System` without constraints.
#' @examples
#' x <- Variable("x")
#' sys <- complete(System(list(), costs = list((x - 2)^2), name = "opt",
#'                        constraints = list(Inequality(x, 1))))
#' pen <- constraints_to_penalties(sys, weights = 1e3)
#' costs(pen)
#' sol <- solve(OptimizationProblem(pen, list(x ~ 0)))
#' sol$u
#' @export
constraints_to_penalties <- function(sys, weights = 1) {
  n <- length(constraints(sys))
  vec <- is.list(weights) && !is_sym(weights) || (is.numeric(weights) && length(weights) > 1L)
  if (vec && length(weights) != n)
    .stopf(paste0("Expected `weights` to be a scalar or have one entry per constraint of the system ",
                  "(including subsystem constraints). Got %d weights for %d constraints."), length(weights), n)
  out <- .opt_penalties(sys, if (vec) as.list(weights) else weights, vec)
  if (iscomplete(sys)) out <- .update_system(out, index_cache = .build_index_cache(out))
  out
}

.opt_penalties <- function(sys, weights, vec) {
  cstrs <- .f(sys, "constraints")
  own <- if (vec) weights[seq_along(cstrs)] else rep(list(weights), length(cstrs))
  pen <- list(); new_ps <- list()
  known <- .keys(c(.f(sys, "ps"), .f(sys, "unknowns")))
  for (i in seq_along(cstrs)) {
    w <- as_sym(own[[i]])
    res <- constraint_residual(cstrs[[i]])
    pen[[i]] <- if (inherits(cstrs[[i]], "Inequality")) w * Max(res, 0)^2 else w * res^2
    for (v in .collect_vars(w)) if (!(.key(v) %in% c(known, .keys(new_ps)))) new_ps[[length(new_ps) + 1L]] <- v
  }
  subs <- .f(sys, "systems")
  if (length(subs)) {
    offset <- length(cstrs)
    for (j in seq_along(subs)) {
      sw <- weights
      if (vec) {
        nsub <- length(constraints(subs[[j]]))
        sw <- weights[offset + seq_len(nsub)]
        offset <- offset + nsub
      }
      subs[[j]] <- .opt_penalties(subs[[j]], sw, vec)
    }
  }
  .update_system(sys, costs = c(.f(sys, "costs"), pen), constraints = list(), systems = subs,
                 ps = c(.f(sys, "ps"), new_ps))
}

# ---------------------------------------------------------------------------
# Tunable parameters

#' Tunable parameters
#'
#' `tunable_parameters(sys)` returns the parameters of `sys` that are marked
#' tunable (`Parameter("k", tunable = TRUE)`); parameters without the
#' metadata count as tunable when `default = TRUE`.
#' `subset_tunables(sys, new_tunables)` returns the system with only
#' `new_tunables` tunable; the other tunable parameters are marked not
#' tunable. `reorder_dimension_by_tunables(sys, arr, syms)` takes an array
#' whose dimension `dim` is ordered as `tunable_parameters(sys)` and
#' reorders it as `syms`, a permutation of the tunable parameters. R arrays
#' are not modified in place, so there is no separate mutating variant.
#'
#' @param sys A `System`; complete for `subset_tunables()` and
#'   `reorder_dimension_by_tunables()`.
#' @param p The parameters to filter.
#' @param default Whether parameters without `tunable` metadata are tunable.
#' @param new_tunables A list of parameters, a subset of the tunable ones.
#' @param arr A vector, matrix or array.
#' @param syms A list of parameters, a permutation of the tunable ones.
#' @param dim The dimension of `arr` to reorder.
#' @return A list of parameters, a `System`, or the reordered array.
#' @examples
#' x <- Variable("x", t_nounits)
#' a <- Parameter("a"); b <- Parameter("b"); cc <- Parameter("cc", tunable = FALSE)
#' sys <- mtkcompile(System(list(D_nounits(x) ~ -a * x + b + cc), t_nounits, name = "tun"))
#' tunable_parameters(sys)
#' tunable_parameters(subset_tunables(sys, list(b)))
#' reorder_dimension_by_tunables(sys, c(10, 20), list(b, a))
#' @export
tunable_parameters <- function(sys, p = .parameters_of(sys, TRUE), default = TRUE) {
  Filter(function(x) {
    if (!isparameter(x)) return(FALSE)
    tu <- .meta_get(x, "tunable")
    if (is.null(tu)) isTRUE(default) else isTRUE(tu)
  }, .flatten_vars(p))
}

#' @rdname tunable_parameters
#' @export
subset_tunables <- function(sys, new_tunables) {
  if (!iscomplete(sys)) stop("System must be `complete` before changing tunables.", call. = FALSE)
  cur <- tunable_parameters(sys, .parameters_of(sys))
  new_tunables <- .flatten_vars(new_tunables)
  extra <- setdiff(.keys(new_tunables), .keys(cur))
  if (length(extra))
    stop("New tunable parameters must be a subset of the current tunable parameters. Found tunable ",
         "parameters not in the system: ",
         paste(.nl_names(new_tunables[.keys(new_tunables) %in% extra]), collapse = ", "), ".", call. = FALSE)
  drop <- setdiff(.keys(cur), .keys(new_tunables))
  ps <- lapply(.f(sys, "ps"), function(p) {
    if (!(.key(p) %in% drop)) return(p)
    m <- .meta(p)
    m$tunable <- FALSE
    attr(p, "mtk_meta") <- m
    p
  })
  sys <- .update_system(sys, ps = ps)
  .update_system(sys, index_cache = .build_index_cache(sys))
}

#' @rdname tunable_parameters
#' @export
reorder_dimension_by_tunables <- function(sys, arr, syms, dim = 1L) {
  if (!iscomplete(sys))
    stop("A completed system is required. Call `complete` or `mtkcompile` on the system.", call. = FALSE)
  tun <- .keys(tunable_parameters(sys))
  syms <- .flatten_vars(syms)
  idx <- match(.keys(syms), tun)
  if (anyNA(idx))
    stop("`syms` must be a permutation of `tunable_parameters(sys)`. Found ",
         paste(.nl_names(syms[is.na(idx)]), collapse = ", "), " which is not a tunable parameter.", call. = FALSE)
  if (is.null(dim(arr))) return(arr[idx])
  args <- rep(list(quote(expr = )), length(dim(arr)))
  args[[dim]] <- idx
  do.call(`[`, c(list(arr), args, list(drop = FALSE)))
}

# ---------------------------------------------------------------------------
# OptimizationFunction and OptimizationProblem

.opt_check_system <- function(sys, what) {
  .nl_check_system(sys, what, cost = TRUE)
  if (length(Filter(function(e) inherits(e, "Equation"), equations(sys))))
    .stopf("`%s` does not support systems with equations; write them as constraints.", what)
  invisible(TRUE)
}

# mtkcompile() of a system without equations leaves it without unknowns in
# the base pipeline; the optimization variables are those of the model.
.opt_system <- function(sys) {
  if (inherits(sys, "System") && iscomplete(sys) && !length(unknowns(sys)) && length(costs(sys))) {
    par <- .f(sys, "parent")
    if (!is.null(par) && length(unknowns(par))) return(complete(toggle_namespacing(par, TRUE)))
  }
  sys
}

#' Optimization functions and problems
#'
#' `OptimizationProblem(sys, op)` builds the problem of minimizing the cost
#' of a complete system over its unknowns, subject to its constraints
#' (equations `g ~ h` and inequalities, see [Inequality()]) and to the
#' `bounds` metadata of the unknowns. `op` gives the starting values of the
#' unknowns and the parameter values. The gradient of the cost and the
#' Jacobian of the constraints are generated symbolically, and the Hessians
#' with `hess = TRUE` and `cons_h = TRUE`.
#'
#' `OptimizationFunction(sys)` builds the functions without values: fields
#' `f`, `grad`, `hess`, `cons`, `cons_j` and `cons_h`, each a function of
#' `(u, p)`. `MultiObjectiveOptimizationFunction(sys)` keeps the costs
#' separate: `f` returns the vector of costs, `jac` its Jacobian and `hess`
#' one function per cost; `OptimizationProblem(sys, op, multiobjective =
#' TRUE)` builds a problem on it.
#'
#' The problem has `lb` and `ub` (bounds of the unknowns, `NULL` when none is
#' bounded), `lcons` and `ucons` (bounds of the constraint function: 0 and 0
#' for equations, `-Inf` and 0 for inequalities) and `int` (which unknowns
#' are integer).
#'
#' `OptimizationProblem(f, u0, p)` with an R function `f(u, p)` builds a
#' problem without a system; `grad`, `cons`, `lcons`, `ucons` may be given.
#'
#' @param sys A complete `System` with costs, an R function `f(u, p)`, or an
#'   `OptimizationFunction`.
#' @param op Operating point: starting values and parameter values.
#' @param ... Solver options stored with the problem.
#' @param grad,hess Generate the gradient and Hessian of the cost; for raw
#'   functions, R functions of `(u, p)`.
#' @param jac Generate the Jacobian of the cost vector (multiobjective).
#' @param cons_j,cons_h Generate the Jacobian and Hessians of the
#'   constraints.
#' @param cons A constraint function `(u, p)` for raw problems.
#' @param lcons,ucons Bounds of the constraint function for raw problems.
#' @param weights Weights of the costs: the cost becomes
#'   `sum(weights * costs)`. Numbers or parameters of the system.
#' @param lb,ub Bounds of the unknowns; both or neither.
#' @param u0,p Starting values and parameters for raw problems.
#' @param multiobjective Keep the costs as a vector objective.
#' @param guesses,missing_guess_value As in [NonlinearProblem()].
#' @param sense `MinSense` (minimize, the default) or `MaxSense`; see
#'   [optimization_interface].
#' @return An `OptimizationFunction`, `MultiObjectiveOptimizationFunction`
#'   or `OptimizationProblem`.
#' @examples
#' x <- Variable("x"); y <- Variable("y"); a <- Parameter("a"); b <- Parameter("b")
#' sys <- complete(System(list(), costs = list((a - x)^2 + b * (y - x^2)^2), name = "rosenbrock"))
#' prob <- OptimizationProblem(sys, list(x ~ -1, y ~ 2, a ~ 1, b ~ 100))
#' sol <- solve(prob)
#' sol$u; sol$objective
#' # with a constraint: stay inside the unit disk
#' csys <- complete(System(list(), costs = list((a - x)^2 + b * (y - x^2)^2), name = "rosenbrock",
#'                         constraints = list(Inequality(x^2 + y^2, 1))))
#' csol <- solve(OptimizationProblem(csys, list(x ~ 0, y ~ 0, a ~ 1, b ~ 100)))
#' csol$u
#' @export
OptimizationFunction <- function(sys, grad = FALSE, hess = FALSE, cons_j = FALSE, cons_h = FALSE,
                                 weights = NULL, cons = NULL, ...) {
  if (inherits(sys, "AbstractNLPModel")) return(.nlp_optimization_function(sys, ...))
  if (is.function(sys)) {
    return(structure(list(f = sys, grad = if (is.function(grad)) grad, hess = if (is.function(hess)) hess,
                          cons = cons, cons_j = if (is.function(cons_j)) cons_j,
                          cons_h = if (is.function(cons_h)) cons_h, sys = NULL,
                          adtype = if (inherits(grad, "ADType")) grad),
                     class = c("OptimizationFunction", "SciMLFunction")))
  }
  sys <- .opt_system(sys)
  .opt_check_system(sys, "OptimizationFunction")
  if (!is.null(weights)) sys <- .opt_with_weights(sys, weights)
  has_cons <- length(constraints(sys)) > 0L
  structure(list(f = generate_cost(sys),
                 grad = if (isTRUE(grad)) generate_cost_gradient(sys),
                 hess = if (isTRUE(hess)) generate_cost_hessian(sys),
                 cons = if (has_cons) generate_cons(sys),
                 cons_j = if (has_cons && isTRUE(cons_j)) generate_constraint_jacobian(sys),
                 cons_h = if (has_cons && isTRUE(cons_h)) generate_constraint_hessian(sys),
                 cons_expr = if (has_cons) canonical_constraints(sys),
                 expr = cost(sys), sys = sys,
                 observed = function(expr) build_explicit_observed_function(sys, expr)),
            class = c("OptimizationFunction", "SciMLFunction"))
}

#' @rdname OptimizationFunction
#' @export
MultiObjectiveOptimizationFunction <- function(sys, jac = FALSE, hess = FALSE, cons_j = FALSE,
                                               cons_h = FALSE, weights = NULL, cons = NULL, ...) {
  if (is.function(sys)) {
    return(structure(list(f = sys, jac = if (is.function(jac)) jac, hess = if (is.function(hess)) hess,
                          cons = cons, cons_j = if (is.function(cons_j)) cons_j,
                          cons_h = if (is.function(cons_h)) cons_h, sys = NULL),
                     class = c("MultiObjectiveOptimizationFunction", "OptimizationFunction", "SciMLFunction")))
  }
  sys <- .opt_system(sys)
  .opt_check_system(sys, "MultiObjectiveOptimizationFunction")
  if (!is.null(weights))
    stop("`weights` scalarizes the costs of `sys`; a multiobjective objective keeps them separate. ",
         "Use `OptimizationFunction` for a weighted-sum objective.", call. = FALSE)
  has_cons <- length(constraints(sys)) > 0L
  structure(list(f = generate_multiobjective_cost(sys),
                 jac = if (isTRUE(jac)) generate_multiobjective_jacobian(sys),
                 hess = if (isTRUE(hess)) generate_multiobjective_hessian(sys),
                 cons = if (has_cons) generate_cons(sys),
                 cons_j = if (has_cons && isTRUE(cons_j)) generate_constraint_jacobian(sys),
                 cons_h = if (has_cons && isTRUE(cons_h)) generate_constraint_hessian(sys),
                 cons_expr = if (has_cons) canonical_constraints(sys),
                 expr = costs(sys), sys = sys,
                 observed = function(expr) build_explicit_observed_function(sys, expr)),
            class = c("MultiObjectiveOptimizationFunction", "OptimizationFunction", "SciMLFunction"))
}

#' @rdname OptimizationFunction
#' @export
OptimizationProblem <- function(sys, op = list(), ..., lb = NULL, ub = NULL, multiobjective = FALSE,
                                grad = TRUE, hess = FALSE, cons_j = TRUE, cons_h = FALSE, weights = NULL,
                                u0 = NULL, p = NULL, cons = NULL, lcons = NULL, ucons = NULL,
                                guesses = list(), missing_guess_value = NULL, sense = MinSense) {
  if (inherits(sys, "AbstractNLPModel")) return(.nlp_optimization_problem(sys, ...))
  dots <- list(...)
  if (!inherits(sense, "ObjSense")) stop("`sense` must be MinSense or MaxSense", call. = FALSE)
  if (inherits(sys, "OptimizationFunction") && !is.null(sys$sys)) {
    # a function built from a system: the problem is built from that system
    if (inherits(sys, "MultiObjectiveOptimizationFunction")) multiobjective <- TRUE
    sys <- sys$sys
  }
  if (is.function(sys) || inherits(sys, "OptimizationFunction")) {
    f <- if (is.function(sys)) OptimizationFunction(sys, grad = grad, hess = hess, cons = cons, cons_j = cons_j,
                                                    cons_h = cons_h) else sys
    if (is.null(u0)) u0 <- op
    if (is.null(p) && length(dots) && (is.null(names(dots)) || !nzchar(names(dots)[1L]))) {
      p <- dots[[1L]]
      dots <- dots[-1L]
    }
    prob <- .new_problem("OptimizationProblem", f = f, u0 = unlist(u0), p = if (is.null(p)) numeric() else p,
                         sys = NULL, kwargs = dots)
    prob$sense <- sense
    prob$lb <- lb; prob$ub <- ub
    if (!is.null(f$cons)) {
      m <- length(f$cons(prob$u0, prob$p))
      prob$lcons <- if (is.null(lcons)) rep(0, m) else as.numeric(lcons)
      prob$ucons <- if (is.null(ucons)) rep(0, m) else as.numeric(ucons)
    }
    return(prob)
  }
  sys <- .opt_system(sys)
  .opt_check_system(sys, "OptimizationProblem")
  if (xor(is.null(lb), is.null(ub))) stop("Expected both `lb` and `ub` to be supplied", call. = FALSE)
  fn <- if (isTRUE(multiobjective)) MultiObjectiveOptimizationFunction(sys, jac = grad, hess = hess, cons_j = cons_j,
                                                                      cons_h = cons_h, weights = weights)
        else OptimizationFunction(sys, grad = grad, hess = hess, cons_j = cons_j, cons_h = cons_h,
                                  weights = weights)
  sys <- fn$sys
  us <- unknowns(sys)
  vals <- .nl_values(sys, op, guesses, missing_guess_value)
  int <- vapply(us, function(u) isTRUE(.meta_get(u, "integer")) || isTRUE(.meta_get(u, "binary")), TRUE)
  symbolic_bounds <- FALSE
  if (is.null(lb)) {
    b <- .nl_bounds(sys, NULL, NULL, vals$values)
    symbolic_bounds <- isTRUE(b$symbolic)
    lo <- b$lb %||% rep(-Inf, length(us)); hi <- b$ub %||% rep(Inf, length(us))
    bin <- vapply(us, function(u) isTRUE(.meta_get(u, "binary")), TRUE)
    lo[bin] <- 0; hi[bin] <- 1
  } else {
    if (length(lb) != length(us))
      stop("Expected `lb` to be of the same length as the vector of optimization variables", call. = FALSE)
    if (length(ub) != length(us))
      stop("Expected `ub` to be of the same length as the vector of optimization variables", call. = FALSE)
    lo <- as.numeric(lb); hi <- as.numeric(ub)
  }
  if (all(lo == -Inf) && all(hi == Inf)) lo <- hi <- NULL
  cstr <- constraints(sys)
  prob <- .new_problem("OptimizationProblem", f = fn, u0 = vals$u0, p = vals$p, sys = sys, kwargs = dots)
  prob$sense <- sense
  prob$op <- vals$op
  prob$lb <- lo; prob$ub <- hi
  prob$symbolic_bounds <- symbolic_bounds
  prob$int <- int
  if (length(cstr)) {
    prob$lcons <- vapply(cstr, function(cn) if (inherits(cn, "Inequality")) -Inf else 0, 0)
    prob$ucons <- numeric(length(cstr))
  }
  prob$guesses <- guesses
  prob$missing_guess_value <- missing_guess_value
  prob
}
