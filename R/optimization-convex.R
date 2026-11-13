# The conic backend of ConvexOptimization.jl: trace the objective and the
# cone constraints of a ConvexOptimizationProblem symbolically, certify
# convexity, lower norm, exp and log atoms through their epigraphs, extract
# the cone program as an affine function of the parameters (disciplined
# parametrized programming) and solve it with Clarabel.
#
# Ported from ConvexOptimization.jl 0.1.0 (lib/ConvexOptimization of
# Optimization.jl, https://github.com/SciML/Optimization.jl), Copyright (c)
# 2023 Vaibhav Dixit and contributors, MIT License (inst/CONVEX-LICENSES).

#' Cone constraints and the conic backend of a convex problem
#'
#' `ConeConstraint(g, set)` is one convex cone constraint of a
#' [ConvexOptimizationProblem()]: `g(u, p)` returns the affine map whose
#' image must lie in the MathOptInterface vector cone `set` ([MOI]`$Zeros`,
#' `$Nonnegatives`, `$Nonpositives`, `$SecondOrderCone`, ...). The output
#' length of `g` must equal `MOI$dimension(set)`. Each `ConeConstraint` maps
#' to one entry of the solution's `dual`, in the order the constraints are
#' given, already in the user's variables: `>=` is `$Nonnegatives`, `<=`
#' `$Nonpositives` and `==` `$Zeros`. `set` is read once when the cache is
#' built; a set carrying data derived from the parameters
#' (`MOI$PowerCone(p[1])`) is not updated by [reinit()]; rebuild with
#' `solve(remake(prob, p = ...), alg)` instead.
#'
#' `ConvexMOI(optimizer_constructor)` certifies convexity (SymbolicAnalysis
#' rules, see [analyze()]), lowers the objective and each `ConeConstraint` to
#' a cone and solves with the optimizer ([ClarabelOptimizer()]). The
#' objective may be affine or contain atoms lowered through their epigraph:
#' `norm(A %*% u - b, "2")` introduces an epigraph variable `tau` with
#' `(tau, A u - b)` in a second order cone and minimizes `tau`. The atoms are
#' `norm(w, type)` of a vector expression `w` in `u` with the types of
#' [base::norm()] (`"O"` or `"1"` the one norm, `"I"` or `"M"` the infinity
#' norm, `"F"`, `"E"` or `"2"` the two norm, or the order `1`, `2`, `Inf`),
#' lowered to `NormOneCone`, `SecondOrderCone` or `NormInfinityCone`, and
#' `exp(w)` and `log(w)` of scalar affine `w` (`ExponentialCone`); `log` is
#' concave, so it is bounded below and must enter the objective negatively.
#'
#' The functions are traced once, on symbolic vectors `u` and `p` that
#' support arithmetic, indexing, `sum`, `exp`, `log`, `%*%` with numeric
#' matrices and `norm` (the last two inside the traced functions only).
#' Parameters are first class: [reinit()] re-solves the cache returned by
#' [init()] at a new `p` without symbolic work, evaluating the cone data as
#' the affine function of `p` extracted at `init`. The solution has `u`,
#' `objective`, `dual` (one vector per `ConeConstraint`, `NaN` when the
#' solve gives no dual solution), `retcode`, `original` (the solved model,
#' with its MathOptInterface statuses and the Clarabel solver) and `cache`.
#'
#' @param g The map `g(u, p)`.
#' @param set A vector set of [MOI].
#' @param optimizer_constructor [ClarabelOptimizer] or a call of it.
#' @return A `ConeConstraint` or a `ConvexMOI` algorithm.
#' @examples
#' cons <- list(ConeConstraint(function(u, p) u[1] + u[2] - 1, MOI$Zeros(1)),
#'              ConeConstraint(function(u, p) c(u[1], u[2]), MOI$Nonnegatives(2)))
#' prob <- ConvexOptimizationProblem(function(u, p) u[1] + 2 * u[2], c(0.5, 0.5), constraints = cons)
#' sol <- solve(prob, ConvexMOI())
#' sol$u; sol$objective; sol$dual
#' A <- matrix(c(1, 1), 2, 1)
#' sol <- solve(ConvexOptimizationProblem(function(u, p) norm(A %*% u - c(0, 2), "2"), 0), ConvexMOI())
#' sol$u
#' @export
ConeConstraint <- function(g, set) {
  if (!inherits(set, "MOI_AbstractVectorSet"))
    stop("`set` must be a MathOptInterface vector set such as MOI$Zeros(1)", call. = FALSE)
  structure(list(g = g, set = set), class = "ConeConstraint")
}

#' @rdname ConeConstraint
#' @export
ConvexMOI <- function(optimizer_constructor = ClarabelOptimizer)
  structure(list(optimizer_constructor = optimizer_constructor, name = "ConvexMOI"),
            class = c("ConvexMOI", "AbstractConvexOptAlgorithm", "OptimizationAlgorithm"))

#' @export
print.ConeConstraint <- function(x, ...) {
  cat("ConeConstraint in ", .moi_format_set(x$set), "\n", sep = "")
  invisible(x)
}

.CVX_PNOTE <- " (`.p[i]` denotes `p[i]`; `.tau[k]` is the epigraph variable of atom k.)"

#' @export
init.ConvexOptimizationProblem <- function(prob, alg = NULL, ...) {
  if (!inherits(alg, "AbstractConvexOptAlgorithm"))
    stop("a ConvexOptimizationProblem is solved with ConvexMOI()", call. = FALSE)
  tr <- .cvx_trace_problem(prob)
  # with symbolic parameters a product p[1] u[2] has two non constant factors
  # for the DCP rules, so the parametric path is certified by the structural
  # checks of the extraction instead (stronger than DCP)
  analysis <- if (!length(tr$params)) .cvx_certify(prob, tr) else NULL
  dpp <- .cvx_dpp_extract(prob, tr)
  model <- .moi_model_new(.moi_optimizer(alg$optimizer_constructor))
  refs <- .cvx_build_moi(model, dpp, .cvx_theta(prob$p))
  cache <- new.env(parent = emptyenv())
  cache$f <- prob$f; cache$u0 <- prob$u0; cache$p <- .cvx_cachep(prob$p); cache$alg <- alg
  cache$analysis <- analysis; cache$dpp <- dpp; cache$model <- model
  cache$xvars <- refs$x; cache$conrefs <- refs$conrefs; cache$atomrefs <- refs$atomrefs
  class(cache) <- c("ConvexOptimizationCache", "AbstractOptimizationCache")
  cache
}

#' @export
solve.ConvexOptimizationProblem <- function(a, b, ...) {
  if (missing(b)) stop("a ConvexOptimizationProblem is solved with ConvexMOI()", call. = FALSE)
  .cvx_solve(init(a, b, ...))
}

#' @export
print.ConvexOptimizationCache <- function(x, ...) {
  cat("ConvexOptimizationCache: ", x$dpp$n, " variables, ", length(x$dpp$consets), " cone constraints, ",
      length(x$dpp$atomsets), " lowered atoms\n", sep = "")
  invisible(x)
}

.cvx_theta <- function(p) if (inherits(p, "NullParameters")) numeric() else as.numeric(p)
.cvx_cachep <- function(p) if (inherits(p, "NullParameters")) p else as.numeric(p)

.cvx_solve <- function(cache) {
  model <- cache$model
  .moi_optimize(model)
  ret <- .cvx_moi_status_to_retcode(model$termination_status)
  if (model$result_count >= 1L) {
    u <- model$solution$x[cache$xvars]
    objective <- model$objective_value
  } else {
    u <- rep(NaN, length(cache$xvars))
    objective <- NaN
  }
  # a dual is reported only when it is one: on a failed solve the dual status
  # is a certificate or another point, not a dual solution
  dual <- if (model$result_count >= 1L && model$dual_status %in% c("FEASIBLE_POINT", "NEARLY_FEASIBLE_POINT"))
    lapply(cache$conrefs, function(ci) .moi_constraint_dual(model, ci))
  else lapply(cache$dpp$consets, function(s) rep(NaN, s$dimension))
  prob <- list(u0 = cache$u0, p = cache$p, f = cache$f)
  .nl_new_solution(c("OptimizationSolution", "NonlinearSolution"), prob, cache$alg, u, ret,
                   stats = OptimizationStats(),
                   extra = list(objective = objective, original = model, dual = dual, cache = cache))
}

.cvx_moi_status_to_retcode <- function(s) {
  if (s %in% c("OPTIMAL", "LOCALLY_SOLVED", "ALMOST_OPTIMAL", "ALMOST_LOCALLY_SOLVED")) return("Success")
  if (s %in% c("INFEASIBLE", "DUAL_INFEASIBLE", "LOCALLY_INFEASIBLE", "INFEASIBLE_OR_UNBOUNDED")) return("Infeasible")
  if (s == "TIME_LIMIT") return("MaxTime")
  if (s %in% c("ITERATION_LIMIT", "NODE_LIMIT", "SLOW_PROGRESS")) return("MaxIters")
  if (s %in% c("NUMERICAL_ERROR", "INVALID_MODEL", "OTHER_ERROR")) return("Failure")
  "Default"
}

.cvx_is_symbolic_map <- function(v) is.list(v) && (any(vapply(v, function(e) inherits(e, "formula"), TRUE)) ||
                                                      !is.null(names(v)))

.cvx_reinit <- function(cache, u0 = NULL, p = NULL) {
  if (!is.null(u0)) {
    if (.cvx_is_symbolic_map(u0))
      stop("`reinit()` on a `ConvexOptimizationCache` does not support symbolic maps; pass `u0` as a vector of ",
           "values.", call. = FALSE)
    if (length(u0) != length(cache$u0))
      stop("`reinit()` got a `u0` of length ", length(u0), ", but this cache was canonicalized for ",
           length(cache$u0), " optimization variables. The variable count is fixed at `init`; use ",
           "`solve(remake(prob, u0 = ...), alg)` to canonicalize a new problem.", call. = FALSE)
    cache$u0 <- u0
  }
  if (is.null(p)) return(invisible(cache))
  .cvx_validate_p(cache, p)
  theta <- as.numeric(p)
  # build first, assign second: a rejected p leaves the cache as it was
  refs <- .cvx_build_moi(cache$model, cache$dpp, theta)
  cache$xvars <- refs$x; cache$conrefs <- refs$conrefs; cache$atomrefs <- refs$atomrefs
  cache$p <- theta
  invisible(cache)
}

.cvx_validate_p <- function(cache, p) {
  if (.cvx_is_symbolic_map(p))
    stop("`reinit()` on a `ConvexOptimizationCache` does not support symbolic maps (`p = list(sym ~ val, ...)`); ",
         "this problem has no symbolic origin. Pass `p` as a vector of values in the order of the original `p`.",
         call. = FALSE)
  if (inherits(cache$p, "NullParameters"))
    stop("`reinit()` got `p`, but this problem was canonicalized with `NullParameters`: there are no parameters ",
         "in the cached problem data to update. Build the problem with a parameter vector, or use ",
         "`solve(remake(prob, p = ...), alg)`.", call. = FALSE)
  if (!is.numeric(p) || !is.null(dim(p)))
    stop("`reinit()` expects `p` as a numeric vector; got a `", class(p)[1L], "`.", call. = FALSE)
  if (length(p) != cache$dpp$m)
    stop("`reinit()` got a `p` of length ", length(p), ", but this cache was canonicalized for ", cache$dpp$m,
         " parameters. The parameter count is fixed at `init`; use `solve(remake(prob, p = ...), alg)` to ",
         "canonicalize a new problem.", call. = FALSE)
  if (!all(is.finite(p)))
    stop("`reinit()` got a non-finite entry in `p` (", paste(deparse(p), collapse = ""), "); the cone ",
         "program's data would not be well defined.", call. = FALSE)
  invisible()
}

# ---------------------------------------------------------------------------
# Tracing and atom lowering

.cvx_symbolic_vars <- function(prob) {
  n <- length(prob$u0)
  vars <- .symvec(lapply(seq_len(n), function(i) Symbol(paste0("u[", i, "]"))))
  params <- if (inherits(prob$p, "NullParameters")) list() else {
    if (!is.numeric(prob$p) || !is.null(dim(prob$p)))
      stop("This backend traces `p` as a vector of symbolic parameters, so a `ConvexOptimizationProblem`'s `p` ",
           "must be a numeric vector (an `AbstractVector`) or `NullParameters()`; got a `", class(prob$p)[1L],
           "`. Flatten any other container.", call. = FALSE)
    lapply(seq_along(prob$p), function(i) Symbol(paste0(".p[", i, "]")))
  }
  list(vars = vars, cols = .symvec_elems(vars), params = params)
}

.cvx_asvec <- function(v) {
  if (is.null(v)) return(list())
  .symvec_elems(v)
}

.cvx_trace_problem <- function(prob) {
  sv <- .cvx_symbolic_vars(prob)
  psym <- .symvec(sv$params)
  objv <- .cvx_asvec(.cvx_traceable(prob$f$f)(sv$vars, psym))
  if (length(objv) != 1L) stop("the objective must return one value, not ", length(objv), call. = FALSE)
  obj <- objv[[1L]]
  low <- .cvx_epigraph_lower(obj)
  consvals <- if (is.null(prob$constraints)) NULL else
    lapply(prob$constraints, function(con) .cvx_asvec(.cvx_traceable(con$g)(sv$vars, psym)))
  list(vars = sv$vars, cols = sv$cols, params = sv$params, obj = obj, objl = low$objl, taus = low$taus,
       atoms = low$atoms, consvals = consvals, allcols = c(sv$cols, low$taus))
}

.cvx_lowerable <- function(e) is_sym(e) && .type(e) == "Function" && .name(e) %in% c("norm", "exp", "log")

.cvx_collect_atoms <- function(acc, e) {
  if (!is_sym(e) || .cvx_is_num(e) || .type(e) == "Symbol") return(acc)
  if (.cvx_lowerable(e)) {
    if (!(.key(e) %in% .keys(acc))) acc[[length(acc) + 1L]] <- e
    return(acc)
  }
  for (a in .args(e)) acc <- .cvx_collect_atoms(acc, a)
  acc
}

.cvx_has_sqrt <- function(e) {
  if (!is_sym(e) || .cvx_is_num(e) || .type(e) == "Symbol") return(FALSE)
  if (.type(e) == "Pow" && .cvx_is_num(.args(e)[[2L]]) && .cvx_num(.args(e)[[2L]]) == 0.5) return(TRUE)
  any(vapply(.args(e), .cvx_has_sqrt, TRUE))
}

# a sqrt in the traced objective usually is a norm spelled out by hand
.cvx_norm_hint <- function(obj) if (.cvx_has_sqrt(obj))
  paste0(" If this objective computes a norm, write it as `norm(w, type)` with `w` a vector expression in `u` ",
         "(e.g. `A %*% u - b`, `u - c`, `u[1:2] - c`); a norm spelled out with `sqrt` cannot be lowered.") else ""

.cvx_norm_cone <- function(order, dim) {
  pv <- if (identical(.key(order), .key(Symbol(".norm_Inf")))) Inf else if (.cvx_is_num(order)) .cvx_num(order) else
    stop("The order `p` of a `norm(w, p)` atom in the objective must be a constant 1, 2 or Inf; got a ",
         "non-constant expression. Route to a general OptimizationProblem/NLP solver.", call. = FALSE)
  if (pv == 2) return(MOI$SecondOrderCone(dim))
  if (pv == 1) return(MOI$NormOneCone(dim))
  if (is.infinite(pv) && pv > 0) return(MOI$NormInfinityCone(dim))
  stop("This backend lowers `norm(w, p)` in the objective only for p = 1, 2 or Inf; got p = ", format(pv),
       ", which has no corresponding MathOptInterface cone. Reformulate or route to a general ",
       "OptimizationProblem/NLP solver.", call. = FALSE)
}

# convex f(w) <= tau (epigraph, dir +1) or concave f(w) >= tau (hypograph,
# dir -1); ExponentialCone is {(a, b, c): b exp(a / b) <= c, b > 0}, so
# exp(w) <= tau is (w, 1, tau) and log(w) >= tau is (tau, 1, w)
.cvx_atom_lowering <- function(t, tau) {
  args <- .args(t)
  if (.name(t) == "norm") {
    w <- args[-length(args)]
    return(list(rows = c(list(tau), w), set = .cvx_norm_cone(args[[length(args)]], length(w) + 1L), dir = 1L))
  }
  w <- args[[1L]]
  if (.name(t) == "exp") return(list(rows = list(w, as_sym(1), tau), set = MOI$ExponentialCone(), dir = 1L))
  list(rows = list(tau, as_sym(1), w), set = MOI$ExponentialCone(), dir = -1L)
}

.cvx_epigraph_lower <- function(obj) {
  nodes <- .cvx_collect_atoms(list(), obj)
  if (!length(nodes)) return(list(objl = obj, taus = list(), atoms = list()))
  taus <- list(); atoms <- list()
  objl <- obj
  for (k in seq_along(nodes)) {
    tau <- Symbol(paste0(".tau[", k, "]"))
    taus[[k]] <- tau
    atoms[[k]] <- .cvx_atom_lowering(nodes[[k]], tau)
    objl <- subs(objl, nodes[[k]], tau)
  }
  list(objl = objl, taus = taus, atoms = atoms)
}

# ---------------------------------------------------------------------------
# Certification (non parametric path)

.cvx_certify <- function(prob, tr) {
  res <- analyze(tr$objl)
  ok <- if (prob$sense == "MaxSense") res$curvature %in% c("Concave", "Affine") else
    res$curvature %in% c("Convex", "Affine")
  if (!ok)
    stop("Objective is not certified convex for ", prob$sense, ": curvature = ", res$curvature,
         ". Route to a general OptimizationProblem/NLP solver.", .cvx_norm_hint(tr$obj), call. = FALSE)
  cons <- NULL
  if (!is.null(tr$consvals)) {
    cons <- list()
    for (k in seq_along(tr$consvals)) {
      cres <- lapply(tr$consvals[[k]], analyze)
      cv <- vapply(cres, function(r) r$curvature, "")
      if (!all(cv == "Affine"))
        stop("This backend supports affine-in-cone constraints only; got curvatures [", paste(cv, collapse = ", "),
             "] for cone ", .moi_format_set(prob$constraints[[k]]$set), ".", call. = FALSE)
      cons[[k]] <- cres
    }
  }
  list(objective = res, constraints = cons)
}

# ---------------------------------------------------------------------------
# Disciplined parametrized programming: the problem data affine in p

.cvx_fmt <- function(e) paste(format(e), collapse = "")

# an extracted coefficient or constant that still holds an optimization
# variable is the non affine case (the expansion is greedy)
.cvx_check_no_optvars <- function(es, pkeys, tkeys, what) {
  for (e in es) for (v in .cvx_free(e)) {
    if (.key(v) %in% pkeys) next
    if (.key(v) %in% tkeys)
      stop(what, " scales a lowered atom by a non-constant expression: the expansion left the coefficient `",
           .cvx_fmt(e), "`, which contains `", .cvx_fmt(v), "`, the epigraph variable of an atom. An atom may be ",
           "scaled only by a quantity that is affine in `p`. If `", .cvx_fmt(v), "` came from lowering a function ",
           "of `p` alone (as in `exp(p[1]) * u[1]`), pass that quantity as its own parameter instead.", .CVX_PNOTE,
           call. = FALSE)
    stop(what, " is not affine in the optimization variables: the expansion left the coefficient `", .cvx_fmt(e),
         "`, which still contains the optimization variable `", .cvx_fmt(v), "`. Products or nonlinear functions ",
         "of the optimization variables are not disciplined-convex here; route to a general ",
         "OptimizationProblem/NLP solver.", call. = FALSE)
  }
  invisible()
}

.cvx_nonnumeric <- function(...) {
  for (xs in list(...)) for (e in xs) if (!.cvx_is_num(e)) return(e)
  NULL
}

.cvx_floatmat <- function(A) matrix(vapply(A, .cvx_num, 0), nrow(A), ncol(A))

.cvx_theta_affine <- function(es, params, what) {
  le <- .cvx_linear_expansion(es, params)
  allnum <- function(xs) all(vapply(xs, .cvx_is_num, TRUE))
  if (!le$islin || !allnum(le$A) || !allnum(le$b)) {
    bad <- if (le$islin) .cvx_nonnumeric(le$A, le$b) else NULL
    witness <- if (is.null(bad)) "" else paste0(": the expansion left the non-numeric term `", .cvx_fmt(bad), "`")
    stop(what, " is not affine in the parameters `p`", witness, ". Disciplined parametrized programming requires ",
         "the cone program's data to be affine in `p`; rewrite the nonlinear function of `p` as an extra ",
         "parameter, or use `solve(remake(prob, p = ...), alg)` to re-canonicalize at each `p`.", .CVX_PNOTE,
         call. = FALSE)
  }
  list(b0 = vapply(le$b, .cvx_num, 0), B = .cvx_floatmat(le$A))
}

.cvx_dpp_block <- function(rows, cols, params, pkeys, tkeys, what, nonaffine_msg) {
  le <- .cvx_linear_expansion(rows, cols)
  if (!le$islin) stop(nonaffine_msg, call. = FALSE)
  .cvx_check_no_optvars(le$A, pkeys, tkeys, what)
  .cvx_check_no_optvars(le$b, pkeys, tkeys, what)
  bad <- .cvx_nonnumeric(le$A)
  if (!is.null(bad))
    stop(what, " has the parameter-dependent coefficient `", .cvx_fmt(bad), "`. This problem is convex at every ",
         "`p`, but its cone *matrix* moves with `p`, so it cannot be canonicalized once and re-solved from cached ",
         "data. Rewrite so that parameters enter additively (`A %*% u - p` rather than `p * u`), or use ",
         "`solve(remake(prob, p = ...), alg)` to re-canonicalize at each `p`.", .CVX_PNOTE, call. = FALSE)
  th <- .cvx_theta_affine(le$b, params, paste0(what, "'s constant term"))
  list(A = .cvx_floatmat(le$A), b0 = th$b0, B = th$B)
}

.cvx_dpp_extract <- function(prob, tr) {
  n <- length(prob$u0)
  params <- tr$params
  pkeys <- .keys(params); tkeys <- .keys(tr$taus)
  m <- length(params)
  conA <- list(); conb0 <- list(); conB <- list()
  consets <- if (is.null(prob$constraints)) list() else lapply(prob$constraints, function(con) con$set)
  for (k in seq_along(prob$constraints)) {
    con <- prob$constraints[[k]]
    blk <- .cvx_dpp_block(tr$consvals[[k]], tr$cols, params, pkeys, tkeys,
                          paste0("Constraint ", k, " (", .moi_format_set(con$set), ")"),
                          paste0("Constraint ", .moi_format_set(con$set), " is not affine in the variables."))
    if (length(tr$consvals[[k]]) != con$set$dimension)
      stop("Constraint ", k, " returns ", length(tr$consvals[[k]]), " values for a set of dimension ",
           con$set$dimension, ".", call. = FALSE)
    conA[[k]] <- blk$A; conb0[[k]] <- blk$b0; conB[[k]] <- blk$B
  }
  atomA <- list(); atomb0 <- list(); atomB <- list()
  for (j in seq_along(tr$atoms)) {
    at <- tr$atoms[[j]]
    blk <- .cvx_dpp_block(at$rows, tr$allcols, params, pkeys, tkeys,
                          paste0("The argument of atom ", j, " (", .moi_format_set(at$set), ")"),
                          paste0("The argument of a `norm` atom in the objective must be affine in the ",
                                 "optimization variables."))
    atomA[[j]] <- blk$A; atomb0[[j]] <- blk$b0; atomB[[j]] <- blk$B
  }
  lo <- .cvx_linear_expansion(list(tr$objl), tr$allcols)
  if (!lo$islin)
    stop("This backend requires an objective that is affine in the optimization variables and in the epigraph ",
         "variables of its lowered atoms.", call. = FALSE)
  .cvx_check_no_optvars(lo$A, pkeys, tkeys, "The objective")
  .cvx_check_no_optvars(lo$b, pkeys, tkeys, "The objective")
  # the objective coefficients are the one place a parameter may multiply a
  # column: c(theta) is linear data, and a linear objective is convex at every theta
  cth <- .cvx_theta_affine(as.vector(lo$A), params, "An objective coefficient")
  dth <- .cvx_theta_affine(lo$b, params, "The objective's constant term")
  dpp <- list(n = n, m = m, c0 = cth$b0, C = cth$B, d0 = dth$b0[[1L]], dP = as.numeric(dth$B),
              conA = conA, conb0 = conb0, conB = conB, consets = consets,
              atomA = atomA, atomb0 = atomb0, atomB = atomB,
              atomsets = lapply(tr$atoms, function(a) a$set), atomdirs = vapply(tr$atoms, function(a) a$dir, 0L),
              lb = if (is.null(prob$lb)) rep(-Inf, n) else as.numeric(prob$lb),
              ub = if (is.null(prob$ub)) rep(Inf, n) else as.numeric(prob$ub),
              sense = prob$sense)
  .cvx_warn_if_p_unused(dpp)
  dpp
}

# a function that captures a parameter vector from its enclosing scope traces
# to data that does not mention p, so every re-solve would return the first answer
.cvx_warn_if_p_unused <- function(dpp) {
  if (dpp$m == 0) return(invisible())
  used <- any(dpp$C != 0) || any(dpp$dP != 0) || any(vapply(dpp$conB, function(B) any(B != 0), TRUE)) ||
    any(vapply(dpp$atomB, function(B) any(B != 0), TRUE))
  if (used) return(invisible())
  warning("`p` has length ", dpp$m, ", but the traced objective and constraints do not depend on it: the ",
          "canonicalized problem data is constant in `p`, so `reinit(cache, p = ...)` will return the same ",
          "solution at every `p`. Did the objective or a constraint capture a parameter vector from the enclosing ",
          "scope instead of using its own `p` argument?", call. = FALSE)
  invisible()
}

# Numeric assembly shared by the cold build and every reinit, so bounds,
# constraints, atom cones, objective and sense are always refreshed together
.cvx_build_moi <- function(model, dpp, theta) {
  n <- dpp$n; nt <- length(dpp$atomdirs)
  c0 <- dpp$c0 + if (dpp$m) as.numeric(dpp$C %*% theta) else 0
  d <- dpp$d0 + if (dpp$m) .jl_dot(dpp$dP, theta) else 0
  # an epigraph variable bounds its atom on one side only, so substituting it
  # is valid only where the objective pushes it against that bound; checked
  # before the model is touched so a rejected p leaves the cache usable
  sgn <- if (dpp$sense == "MaxSense") -1 else 1
  for (k in seq_along(dpp$atomdirs)) {
    dir <- dpp$atomdirs[k]
    if (!(sgn * dir * c0[n + k] >= 0))
      stop("Lowering an atom through its ", if (dir > 0) "epigraph" else "hypograph", " is valid only when the ",
           "objective is ", if (dir > 0) "nondecreasing" else "nonincreasing", " in it for MinSense (reversed for ",
           "MaxSense); got coefficient ", format(c0[n + k]), " for ", dpp$sense,
           if (dpp$m == 0) "" else paste0(" at p = ", paste(deparse(theta), collapse = "")),
           ". Route to a general OptimizationProblem/NLP solver.", call. = FALSE)
  }
  .moi_empty(model)
  x <- .moi_add_variables(model, n)
  t <- .moi_add_variables(model, nt)
  allx <- c(x, t)
  for (i in seq_len(n)) {
    if (dpp$lb[i] > -Inf) .moi_add_bound(model, x[i], dpp$lb[i], TRUE)
    if (dpp$ub[i] < Inf) .moi_add_bound(model, x[i], dpp$ub[i], FALSE)
  }
  # user constraints over the user columns only and in problem order, so
  # conrefs stays 1:1 with prob$constraints and with sol$dual
  conrefs <- lapply(seq_along(dpp$consets), function(k) {
    b <- dpp$conb0[[k]] + if (dpp$m) as.numeric(dpp$conB[[k]] %*% theta) else 0
    .moi_add_constraint(model, dpp$conA[[k]], x, b, dpp$consets[[k]])
  })
  atomrefs <- lapply(seq_along(dpp$atomsets), function(j) {
    b <- dpp$atomb0[[j]] + if (dpp$m) as.numeric(dpp$atomB[[j]] %*% theta) else 0
    .moi_add_constraint(model, dpp$atomA[[j]], allx, b, dpp$atomsets[[j]])
  })
  coefs <- numeric(length(allx))
  nz <- c0 != 0
  coefs[nz] <- c0[nz]
  .moi_set_objective(model, coefs, d, if (dpp$sense == "MaxSense") "MAX_SENSE" else "MIN_SENSE")
  list(x = x, conrefs = conrefs, atomrefs = atomrefs)
}
