# The conic backend of ConvexOptimization.jl: trace the objective and the
# cone constraints of a ConvexOptimizationProblem symbolically, certify
# convexity, lower norm, exp and log atoms through their epigraphs, extract
# the cone program as an affine function of the parameters (disciplined
# parametrized programming) and solve it with Clarabel.
#
# Ported from ConvexOptimization.jl 0.1.4 (lib/ConvexOptimization of
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
#' sol <- solve(ConvexOptimizationProblem(function(u, p) norm(A %*% u - c(0, 2), "2"), 0),
#'   ConvexMOI())
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
  # for the DCP rules, so the parametric path without nested atoms is
  # certified by the structural checks of the extraction instead
  analysis <- .cvx_certify(prob, tr)
  dpp <- .cvx_dpp_extract(prob, tr)
  opt <- .moi_optimizer(alg$optimizer_constructor)
  # upstream switches Clarabel's chordal decomposition off for the PSD and
  # log determinant cones it builds; the setting survives every rebuild
  psd <- vapply(c(dpp$consets, dpp$atomsets), function(s)
    s$name %in% c("PositiveSemidefiniteConeTriangle", "LogDetConeTriangle"), TRUE)
  if (any(psd)) opt$settings$chordal_decomposition_enable <- FALSE
  model <- .moi_model_new(opt)
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
    stop("`reinit()` expects `p` as an `AbstractVector`; got a `", class(p)[1L], "`.", call. = FALSE)
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
           "must be an `AbstractVector` or `NullParameters`; got a `", class(prob$p)[1L], "`. Wrap a scalar as ",
           "`[p]` and flatten any other container.", call. = FALSE)
    lapply(seq_along(prob$p), function(i) Symbol(paste0(".p[", i, "]")))
  }
  list(vars = vars, cols = .symvec_elems(vars), params = params)
}

.cvx_asvec <- function(v) {
  if (is.null(v)) return(list())
  if (inherits(v, "mtk_symadj")) return(.cvx_adj_row(v))
  .symvec_elems(v)
}

# the scalar objective (`only` of a vector)
.cvx_scalar <- function(v) {
  if (is.matrix(v) || (is.list(v) && !is.null(dim(v))))
    stop("The objective of a ConvexOptimizationProblem must be a scalar; got a (", paste(dim(v), collapse = ", "),
         ") array. Route to a general OptimizationProblem/NLP solver.", call. = FALSE)
  el <- .cvx_asvec(v)
  if (length(el) > 1L) stop("ArgumentError: Collection has multiple elements, must contain exactly 1 element",
                            call. = FALSE)
  if (!length(el)) stop("ArgumentError: Collection is empty, must contain exactly 1 element", call. = FALSE)
  el[[1L]]
}

# Variables of an expression; the order symbol of an infinity norm is a
# constant, not a variable
.cvx_vars <- function(e) {
  vs <- .cvx_free(e)
  vs[.keys(vs) != .key(Symbol(".norm_Inf"))]
}

.cvx_has_optvar <- function(e, optkeys) any(.keys(.cvx_vars(e)) %in% optkeys)

.cvx_is_theta <- function(e, optkeys) is_sym(e) && !.cvx_is_num(e) && length(.cvx_vars(e)) > 0L &&
  !.cvx_has_optvar(e, optkeys)

.CVX_ATOMS <- c("norm", "exp", "log", "abs", "Abs", "abs2", "max", "min", "Max", "Min", "maximum_vec", "minimum_vec",
                "sum_abs", "maximum_abs", "sum_sq", "self_prod", "quad_form", "quad_prod", "eigmax", "eigmin",
                "opnorm", "logdet")

.cvx_is_square <- function(e) .type(e) == "Pow" && .cvx_is_num(.args(e)[[2L]]) && .cvx_num(.args(e)[[2L]]) == 2

.cvx_lowerable <- function(e) is_sym(e) && !.cvx_is_num(e) &&
  ((.type(e) == "Function" && .name(e) %in% .CVX_ATOMS) || .cvx_is_square(e))

.cvx_is_call <- function(e) is_sym(e) && !.cvx_is_num(e) && .type(e) != "Symbol"

.cvx_atom_args <- function(t) {
  nm <- .cvx_fname(t)
  if (nm %in% c("max", "min", "Max", "Min")) .cvx_flatten_fn(.args(t), nm) else .args(t)
}

.cvx_collect_atoms <- function(acc, e) {
  if (!.cvx_is_call(e)) return(acc)
  if (.cvx_lowerable(e)) {
    if (!(.key(e) %in% .keys(acc))) acc[[length(acc) + 1L]] <- e
    return(acc)
  }
  for (a in .args(e)) acc <- .cvx_collect_atoms(acc, a)
  acc
}

# innermost atoms first, each once
.cvx_collect_deep <- function(acc, e) {
  if (!.cvx_is_call(e)) return(acc)
  push <- function(acc) if (.key(e) %in% .keys(acc)) acc else c(acc, list(e))
  if (.cvx_fname(e) %in% c("max", "min", "Max", "Min") && .cvx_lowerable(e)) {
    for (a in .cvx_atom_args(e)) acc <- .cvx_collect_deep(acc, a)
    return(push(acc))
  }
  for (a in .args(e)) acc <- .cvx_collect_deep(acc, a)
  if (.cvx_lowerable(e)) acc <- push(acc)
  acc
}

.cvx_has_inner_atom <- function(t)
  any(vapply(.cvx_atom_args(t), function(a) length(.cvx_collect_atoms(list(), a)) > 0L, TRUE))

# Julia flattens a scalar multiple of u' P u into the product, so the scale
# becomes part of the quadratic form; another factor makes the product an
# ordinary polynomial
.cvx_fold_quad <- function(e) {
  if (!.cvx_is_call(e)) return(e)
  e <- .rebuild(e, lapply(.args(e), .cvx_fold_quad))
  if (.type(e) != "Mul") return(e)
  args <- .args(e)
  isq <- vapply(args, function(a) .type(a) == "Function" && .name(a) %in% c("quad_prod", "self_prod"), TRUE)
  if (!any(isq)) return(e)
  num <- vapply(args, .cvx_is_num, TRUE)
  if (sum(isq) == 1L && all(num | isq)) {
    q <- args[[which(isq)]]
    s <- prod(vapply(args[num], .cvx_num, 0))
    parts <- if (.name(q) == "quad_prod") .cvx_quad_parts(q) else list(v = .args(q), P = diag(1, length(.args(q))))
    return(.cvx_fn("quad_prod", c(parts$v, as.list(as.numeric(s * parts$P)))))
  }
  .rebuild(e, lapply(args, function(a) if (.type(a) == "Function" && .name(a) %in% c("quad_prod", "self_prod")) {
    parts <- if (.name(a) == "quad_prod") .cvx_quad_parts(a) else list(v = .args(a), P = diag(1, length(.args(a))))
    out <- as_sym(0)
    for (i in seq_along(parts$v)) for (j in seq_along(parts$v))
      if (parts$P[i, j] != 0) out <- out + parts$P[i, j] * parts$v[[i]] * parts$v[[j]]
    out
  } else a))
}

# Parameter-only atoms are lifted to `.pdata[k]` columns evaluated at each p
.cvx_lift_param_data <- function(e, optkeys, st) {
  if (!.cvx_is_call(e)) return(e)
  if (.cvx_lowerable(e) && .cvx_is_theta(e, optkeys)) {
    i <- match(.key(e), .keys(st$pdataexprs))
    if (!is.na(i)) return(st$pdatavars[[i]])
    s <- Symbol(paste0(".pdata[", length(st$pdataexprs) + 1L, "]"))
    st$pdatavars[[length(st$pdatavars) + 1L]] <- s
    st$pdataexprs[[length(st$pdataexprs) + 1L]] <- e
    return(s)
  }
  .rebuild(e, lapply(.args(e), .cvx_lift_param_data, optkeys = optkeys, st = st))
}

.cvx_trace_problem <- function(prob) {
  sv <- .cvx_symbolic_vars(prob)
  optkeys <- .keys(sv$cols)
  psym <- .symvec(sv$params)
  .cvx_state$active <- TRUE
  .cvx_state$optkeys <- optkeys
  on.exit({ .cvx_state$active <- FALSE; .cvx_state$optkeys <- character() })
  obj <- .cvx_scalar(.cvx_traceable(prob$f$f)(sv$vars, psym))
  consvals <- if (is.null(prob$constraints)) NULL else
    lapply(prob$constraints, function(con) .cvx_asvec(.cvx_traceable(con$g)(sv$vars, psym)))
  .cvx_state$active <- FALSE
  obj <- .cvx_fold_quad(obj)
  if (!is.null(consvals)) consvals <- lapply(consvals, function(cv) lapply(cv, .cvx_fold_quad))
  st <- new.env(parent = emptyenv())
  st$pdatavars <- list(); st$pdataexprs <- list(); st$taus <- list(); st$atoms <- list()
  if (length(sv$params)) {
    obj <- .cvx_lift_param_data(obj, optkeys, st)
    if (!is.null(consvals)) consvals <- lapply(consvals, function(cv) lapply(cv, .cvx_lift_param_data, optkeys, st))
  }
  low <- .cvx_epigraph_lower(obj, optkeys, st)
  consl <- if (is.null(consvals)) NULL else lapply(seq_along(consvals), function(k) {
    set <- prob$constraints[[k]]$set
    if (!(set$name %in% c("Nonpositives", "Nonnegatives"))) return(consvals[[k]])
    lapply(consvals[[k]], function(gi) .cvx_epigraph_lower(gi, optkeys, st)$e)
  })
  nested <- low$nested || any(vapply(consvals %||% list(), function(cv)
    any(vapply(cv, function(gi) any(vapply(.cvx_collect_deep(list(), gi), .cvx_has_inner_atom, TRUE)), TRUE)), TRUE))
  list(vars = sv$vars, cols = sv$cols, params = sv$params, pdatavars = st$pdatavars, pdataexprs = st$pdataexprs,
       obj = obj, objl = low$e, taus = st$taus, atoms = st$atoms, nested = nested, consvals = consvals,
       consl = consl, allcols = c(sv$cols, st$taus))
}

.cvx_has_sqrt <- function(e) {
  if (!.cvx_is_call(e)) return(FALSE)
  if (.type(e) == "Pow" && .cvx_is_num(.args(e)[[2L]]) && .cvx_num(.args(e)[[2L]]) == 0.5) return(TRUE)
  any(vapply(.args(e), .cvx_has_sqrt, TRUE))
}

# a sqrt in the traced objective usually is a norm spelled out by hand
.cvx_norm_hint <- function(obj) if (.cvx_has_sqrt(obj))
  paste0(" If this objective uses `norm`, spell it `norm(w, type)` with `w` a vector expression built from `u` ",
         "(e.g. `A %*% u - b`, `u - c`, `u[1:2] - c`, or a literal `c(e1, e2, ...)`); `sqrt(sum(w^2))` is a ",
         "different, unsupported expression.") else ""

.cvx_norm_cone <- function(order, dim) {
  pv <- if (identical(.key(order), .key(Symbol(".norm_Inf")))) Inf else if (.cvx_is_num(order)) .cvx_num(order) else
    stop("The order `p` of a `norm(w, p)` atom must be a constant 1, 2 or Inf; got a non-constant expression. ",
         "Route to a general OptimizationProblem/NLP solver.", call. = FALSE)
  if (pv == 2) return(MOI$SecondOrderCone(dim))
  if (pv == 1) return(MOI$NormOneCone(dim))
  if (is.infinite(pv) && pv > 0) return(MOI$NormInfinityCone(dim))
  stop("This backend lowers `norm(w, p)` only for p = 1, 2 or Inf; got p = ", format(pv), ", which has no ",
       "corresponding MathOptInterface cone. Reformulate or route to a general OptimizationProblem/NLP solver.",
       call. = FALSE)
}

.cvx_rsoc <- function(w, tau) list(rows = c(list(tau, as_sym(0.5)), w), set = MOI$RotatedSecondOrderCone(length(w) + 2L),
                                   dir = 1L)

# sqrt of the eigenvalues times the eigenvectors of sym(P): P = L'L
.cvx_psd_factor <- function(P) {
  S <- (P + t(P)) / 2
  E <- eigen(S, symmetric = TRUE)
  # Julia orders the eigenvalues increasingly
  o <- rev(seq_along(E$values))
  lam <- E$values[o]; V <- E$vectors[, o, drop = FALSE]
  lo <- min(lam); hi <- max(lam)
  if (lo < -nrow(P) * .Machine$double.eps * hi)
    stop("`u' * P * u` / `quad_form(u, P)` requires a positive semidefinite `P`, but `sym(P)` has eigenvalue ",
         .jl_show_float(lo), ": the quadratic form is not convex. Route to a general OptimizationProblem/NLP solver.",
         call. = FALSE)
  diag(sqrt(pmax(lam, 0)), length(lam)) %*% t(V)
}

.cvx_quad_lowering <- function(q, tau) {
  parts <- .cvx_quad_parts(q)
  w <- parts$v; P <- parts$P
  if (ncol(P) != length(w))
    stop("`u' * P * u` / `quad_form(u, P)` dimension mismatch: `P` is ", ncol(P), "\u00D7", ncol(P),
         " but `u` has length ", length(w), ".", call. = FALSE)
  L <- .cvx_psd_factor(P)
  Lw <- lapply(seq_len(nrow(L)), function(i) {
    out <- as_sym(0)
    for (j in seq_along(w)) if (L[i, j] != 0) out <- out + L[i, j] * w[[j]]
    out
  })
  .cvx_rsoc(Lw, tau)
}

.cvx_check_symmetric_affine <- function(X, n, m, name) {
  if (m != n) stop("`", name, "` requires a square symmetric matrix; got size (", n, ", ", m, ").", call. = FALSE)
  for (j in seq_len(n)[-1L]) for (i in seq_len(j - 1L)) {
    d <- expand(X[[i, j]] - X[[j, i]])
    if (!(.cvx_is_num(d) && .cvx_num(d) == 0))
      stop("`", name, "` requires a symmetric affine matrix: entries (", i, ", ", j, ") and (", j, ", ", i, ") have ",
           "different coefficients in `u` or `p`.", call. = FALSE)
  }
  n
}

.cvx_matrix_lowering <- function(t, tau) {
  name <- .name(t)
  mp <- .cvx_mat_parts(t)
  m <- mp$nrow; n <- mp$ncol
  X <- matrix(mp$el, m, n)
  rows <- list()
  if (name == "opnorm") {
    for (j in seq_len(m + n)) for (i in seq_len(j)) {
      rows[[length(rows) + 1L]] <- if (j <= m) { if (i == j) tau else as_sym(0) }
        else if (i <= m) X[[i, j - m]] else if (i == j) tau else as_sym(0)
    }
    return(list(rows = rows, set = MOI$PositiveSemidefiniteConeTriangle(m + n), dir = 1L))
  }
  .cvx_check_symmetric_affine(X, m, n, name)
  for (j in seq_len(n)) for (i in seq_len(j)) {
    entry <- X[[i, j]]
    rows[[length(rows) + 1L]] <- switch(name,
      eigmax = (if (i == j) tau else as_sym(0)) - entry,
      eigmin = entry - (if (i == j) tau else as_sym(0)),
      entry)
  }
  if (name == "eigmax") return(list(rows = rows, set = MOI$PositiveSemidefiniteConeTriangle(n), dir = 1L))
  if (name == "eigmin") return(list(rows = rows, set = MOI$PositiveSemidefiniteConeTriangle(n), dir = -1L))
  list(rows = c(list(tau, as_sym(1)), rows), set = MOI$LogDetConeTriangle(n), dir = -1L)
}

# convex f(w) <= tau (epigraph, dir +1) or concave f(w) >= tau (hypograph,
# dir -1); ExponentialCone is {(a, b, c): b exp(a / b) <= c, b > 0}, so
# exp(w) <= tau is (w, 1, tau) and log(w) >= tau is (tau, 1, w)
.cvx_atom_lowering <- function(t, tau) {
  if (.cvx_is_square(t)) return(.cvx_rsoc(list(.args(t)[[1L]]), tau))
  name <- .name(t)
  args <- .args(t)
  switch(name,
    eigmax = , eigmin = , opnorm = , logdet = .cvx_matrix_lowering(t, tau),
    norm = {
      w <- args[-length(args)]
      list(rows = c(list(tau), w), set = .cvx_norm_cone(args[[length(args)]], length(w) + 1L), dir = 1L)
    },
    abs2 = .cvx_rsoc(args[1L], tau),
    sum_sq = , self_prod = .cvx_rsoc(args, tau),
    quad_form = , quad_prod = .cvx_quad_lowering(t, tau),
    max = , Max = , maximum_vec = list(rows = lapply(.cvx_flatten_fn(args, name), function(w) tau - w),
                                   set = MOI$Nonnegatives(length(.cvx_flatten_fn(args, name))), dir = 1L),
    min = , Min = , minimum_vec = list(rows = lapply(.cvx_flatten_fn(args, name), function(w) w - tau),
                                   set = MOI$Nonnegatives(length(.cvx_flatten_fn(args, name))), dir = -1L),
    sum_abs = list(rows = c(list(tau), args), set = MOI$NormOneCone(length(args) + 1L), dir = 1L),
    maximum_abs = list(rows = c(list(tau), args), set = MOI$NormInfinityCone(length(args) + 1L), dir = 1L),
    abs = , Abs = list(rows = list(tau, args[[1L]]), set = MOI$NormOneCone(2L), dir = 1L),
    exp = list(rows = list(args[[1L]], as_sym(1), tau), set = MOI$ExponentialCone(), dir = 1L),
    list(rows = list(tau, as_sym(1), args[[1L]]), set = MOI$ExponentialCone(), dir = -1L))
}

# top down replacement of whole subterms (Symbolics.substitute)
.cvx_replace <- function(e, from, to) {
  if (!length(from) || !is_sym(e)) return(e)
  xreplace(e, lapply(seq_along(from), function(i) list(from[[i]], to[[i]])))
}

.cvx_epigraph_lower <- function(obj, optkeys, st) {
  queue <- .cvx_collect_deep(list(), obj)
  if (!length(queue)) return(list(e = obj, nested = FALSE))
  first <- length(st$atoms) + 1L
  from <- list(); to <- list(); seen <- character()
  decision <- optkeys
  nested <- FALSE
  while (length(queue)) {
    t <- queue[[1L]]; queue <- queue[-1L]
    if (.key(t) %in% seen) next
    if (.cvx_is_theta(t, decision))
      stop("Parameter-only atom `", .cvx_fmt(t), "` could not be lifted to exact parameter data; refusing ",
           "epigraph or hypograph lowering.", call. = FALSE)
    seen <- c(seen, .key(t))
    nested <- nested || .cvx_has_inner_atom(t)
    tsub <- if (length(from)) .rebuild(t, lapply(.args(t), .cvx_replace, from = from, to = to)) else t
    tau <- Symbol(paste0(".tau[", length(st$taus) + 1L, "]"))
    decision <- c(decision, .key(tau))
    low <- .cvx_atom_lowering(tsub, tau)
    st$taus[[length(st$taus) + 1L]] <- tau
    st$atoms[[length(st$atoms) + 1L]] <- low
    from[[length(from) + 1L]] <- t; to[[length(to) + 1L]] <- tau
    for (r in low$rows) {
      found <- .cvx_collect_atoms(list(), r)
      if (!length(found)) next
      nested <- TRUE
      for (f in found) if (!(.key(f) %in% seen) && !(.key(f) %in% .keys(queue))) queue[[length(queue) + 1L]] <- f
    }
  }
  lowered <- .cvx_replace(obj, from, to)
  if (nested) for (j in first:length(st$atoms))
    st$atoms[[j]]$rows <- lapply(st$atoms[[j]]$rows, .cvx_replace, from = from, to = to)
  list(e = lowered, nested = nested)
}

# Numeric value of a parameter-only expression at the values `vals` (named
# by symbol key)
.cvx_eval <- function(e, vals) {
  if (!is_sym(e)) return(as.numeric(e))
  if (.cvx_is_num(e)) return(.cvx_num(e))
  t <- .type(e)
  if (t == "Symbol") {
    if (identical(.key(e), .key(Symbol(".norm_Inf")))) return(Inf)
    v <- vals[[.key(e)]]
    if (is.null(v)) stop("UndefVarError: `", format(e), "` not defined", call. = FALSE)
    return(v)
  }
  a <- lapply(.args(e), .cvx_eval, vals = vals)
  if (t == "Add") return(Reduce(`+`, a))
  if (t == "Mul") return(Reduce(`*`, a))
  if (t == "Pow") return(a[[1L]]^a[[2L]])
  if (t != "Function") return(as.numeric(evalf(xreplace(e, lapply(.cvx_vars(e), function(s)
    list(s, as_sym(vals[[.key(s)]])))))))
  x <- unlist(a)
  last <- function() x[length(x)]
  switch(.name(e),
    exp = exp(x), log = log(x), abs = , Abs = abs(x), abs2 = x^2, log1p = log1p(x), cosh = cosh(x),
    hypot = sqrt(x[1L]^2 + x[2L]^2),
    max = , Max = , maximum_vec = max(x), min = , Min = , minimum_vec = min(x),
    sum_abs = sum(abs(x)), maximum_abs = max(abs(x)), sum_sq = , self_prod = sum(x^2),
    norm = { w <- x[-length(x)]; p <- last(); if (is.infinite(p)) max(abs(w)) else sum(abs(w)^p)^(1 / p) },
    quad_form = , quad_prod = { L <- length(x); n <- as.integer(round((-1 + sqrt(1 + 4 * L)) / 2))
      v <- x[seq_len(n)]; sum(v * (matrix(x[n + seq_len(n * n)], n, n) %*% v)) },
    eigmax = , eigmin = , opnorm = , logdet = {
      m <- last(); M <- matrix(x[-length(x)], m)
      switch(.name(e), eigmax = eigmax(M), eigmin = eigmin(M), opnorm = opnorm(M), logdet = logdet(M))
    },
    rdiv = x[1L] / x[2L],
    as.numeric(evalf(xreplace(e, lapply(.cvx_vars(e), function(s) list(s, as_sym(vals[[.key(s)]])))))))
}

.cvx_extend_theta <- function(dpp, theta) {
  if (!length(dpp$pdataexprs)) return(theta)
  keys <- .keys(dpp$params)[seq_len(dpp$m)]
  vals <- stats::setNames(as.list(theta), keys)
  values <- vapply(dpp$pdataexprs, function(ex) {
    v <- tryCatch(.cvx_eval(ex, vals), error = function(err)
      stop("Could not evaluate parameter-only expression `", .cvx_fmt(ex), "` at p = ", paste(deparse(theta), collapse = ""),
           ": ", conditionMessage(err), call. = FALSE))
    if (length(v) != 1L || !is.finite(v))
      stop("Parameter-only expression `", .cvx_fmt(ex), "` is non-finite at p = ", paste(deparse(theta), collapse = ""),
           ".", call. = FALSE)
    v
  }, 0)
  c(theta, values)
}

# ---------------------------------------------------------------------------
# Certification

.cvx_try_analyze <- function(ex) tryCatch(list(curvature = .sa_curvature(ex), sign = .sa_sign(ex)),
                                          error = function(e) list(curvature = "UnknownCurvature", sign = "AnySign"))

.cvx_curvature_admits <- function(res, sense)
  if (sense == "MaxSense") res$curvature %in% c("Concave", "Affine") else res$curvature %in% c("Convex", "Affine")

.cvx_check_obj <- function(prob, res, tr, suffix = "") {
  if (!.cvx_curvature_admits(res, prob$sense))
    stop("Objective is not certified convex for ", prob$sense, ": curvature = ", res$curvature, ".", suffix,
         " Route to a general OptimizationProblem/NLP solver.", .cvx_norm_hint(tr$obj), call. = FALSE)
  invisible()
}

.cvx_certify <- function(prob, tr) {
  if (tr$nested) return(.cvx_certify_nested(prob, tr))
  if (length(tr$params)) return(NULL)
  res <- .cvx_try_analyze(tr$objl)
  .cvx_check_obj(prob, res, tr)
  list(objective = res, constraints = .cvx_certify_constraints(prob, tr))
}

.cvx_cone_curvature <- function(set) switch(set$name,
  Nonpositives = list(c("Convex", "Affine"), "convex or affine"),
  Nonnegatives = list(c("Concave", "Affine"), "concave or affine"),
  list("Affine", "affine"))

.cvx_certify_constraints <- function(prob, tr) {
  if (is.null(tr$consvals)) return(NULL)
  res <- list()
  for (k in seq_along(tr$consvals)) {
    set <- prob$constraints[[k]]$set
    cres <- lapply(tr$consvals[[k]], .cvx_try_analyze)
    cv <- vapply(cres, function(r) r$curvature, "")
    allowed <- .cvx_cone_curvature(set)
    if (!all(cv %in% allowed[[1L]]))
      stop("Constraint ", k, " in ", .moi_format_set(set), " requires ", allowed[[2L]], " components; got curvatures ",
           "SymbolicAnalysis.Curvature[", paste0("SymbolicAnalysis.", cv, collapse = ", "), "].", if (set$name == "Zeros") " A convex equality is not a convex set." else "",
           " Route to a general OptimizationProblem/NLP solver.", call. = FALSE)
    res[[k]] <- cres
  }
  res
}

.cvx_certify_nested <- function(prob, tr) {
  res <- .cvx_nested_curvature(prob, tr)
  .cvx_check_obj(prob, res, tr, if (!length(tr$params)) "" else paste0(
    " Parameters are kept symbolic in the certificate, so the curvature may also depend on `p`; remaking with a ",
    "numeric `p` does not change that (the nested path never substitutes \u03B8)."))
  list(objective = res, constraints = if (!length(tr$params)) .cvx_certify_constraints(prob, tr) else NULL)
}

# theta subterms replaced by constants of unknown sign; the second form keeps
# the ones of positive sign
.cvx_theta_forms <- function(ex, optkeys) {
  i <- 0L
  subst <- function(keep) {
    walk <- function(t) {
      if (.cvx_is_theta(t, optkeys)) {
        if (!keep(t)) {
          i <<- i + 1L
          return(Symbol(paste0(".thc[", i, "]")))
        }
      }
      if (!.cvx_is_call(t)) return(t)
      .rebuild(t, lapply(.args(t), walk))
    }
    walk(ex)
  }
  list(subst(function(t) FALSE), subst(function(t) .cvx_try_analyze(t)$sign == "Positive"))
}

.cvx_postwalk <- function(e, f) {
  if (!.cvx_is_call(e)) return(f(e))
  f(.rebuild(e, lapply(.args(e), .cvx_postwalk, f = f)))
}

# products drop their theta factors of known sign
.cvx_strip_theta_products <- function(ex, optkeys) .cvx_postwalk(ex, function(t) {
  if (!.cvx_is_call(t) || .type(t) != "Mul") return(t)
  kept <- list(); neg <- FALSE
  for (a in .args(t)) {
    if (.cvx_is_theta(a, optkeys)) {
      s <- .cvx_try_analyze(a)$sign
      if (s == "Positive") next
      if (s == "Negative") { neg <- !neg; next }
    }
    kept[[length(kept) + 1L]] <- a
  }
  r <- if (!length(kept)) as_sym(1) else Reduce(`*`, kept)
  if (neg) -1 * r else r
})

.cvx_nested_curvature <- function(prob, tr) {
  ex <- .cvx_replace(tr$obj, tr$pdatavars, tr$pdataexprs)
  forms <- if (!length(tr$params)) list(ex) else {
    optkeys <- .keys(tr$cols)
    lapply(.cvx_theta_forms(ex, optkeys), .cvx_strip_theta_products, optkeys = optkeys)
  }
  res <- NULL
  for (e in forms) {
    res <- .cvx_analyze_stages(e, prob$sense)
    if (.cvx_curvature_admits(res, prob$sense)) return(res)
  }
  res
}

.cvx_analyze_stages <- function(ex, sense) {
  res <- .cvx_try_analyze(ex)
  if (.cvx_curvature_admits(res, sense)) return(res)
  ex <- tryCatch(.cvx_dcp_reassociate(ex), error = function(e) ex)
  # the third stage of upstream expands array products, which the traced
  # scalar expressions never contain
  .cvx_try_analyze(ex)
}

.cvx_fold <- function(f, es) Reduce(function(a, b) .cvx_fn(f, list(a, b)), es)
.cvx_abs_all <- function(es) lapply(es, function(e) .cvx_fn("abs", list(e)))
.cvx_sum_of <- function(es) Reduce(`+`, es, as_sym(0))

# norm and the reducers rewritten into scalar atoms the DCP rules compose
.cvx_dcp_reassociate <- function(ex) .cvx_postwalk(ex, function(x) {
  if (!.cvx_is_call(x) || .type(x) != "Function") return(x)
  a <- .args(x)
  switch(.name(x),
    norm = {
      p <- a[[length(a)]]; es <- a[-length(a)]
      pv <- if (identical(.key(p), .key(Symbol(".norm_Inf")))) Inf else if (.cvx_is_num(p)) .cvx_num(p) else NA
      if (is.na(pv)) x
      else if (pv == 1) .cvx_sum_of(.cvx_abs_all(es))
      else if (pv == 2) { if (length(es) == 1L) .cvx_fn("abs", es) else .cvx_fold("hypot", es) }
      else if (is.infinite(pv) && pv > 0) .cvx_fold("max", .cvx_abs_all(es))
      else x
    },
    sum_abs = .cvx_sum_of(.cvx_abs_all(a)),
    maximum_abs = .cvx_fold("max", .cvx_abs_all(a)),
    sum_sq = , self_prod = .cvx_sum_of(lapply(a, function(e) .cvx_fn("abs2", list(e)))),
    maximum_vec = .cvx_fold("max", a),
    minimum_vec = .cvx_fold("min", a),
    quad_prod = {
      low <- .cvx_quad_lowering(x, Symbol(".tau[0]"))
      .cvx_sum_of(lapply(low$rows[-(1:2)], function(e) .cvx_fn("abs2", list(e))))
    },
    x)
})

# ---------------------------------------------------------------------------
# Disciplined parametrized programming: the problem data affine in p

.cvx_fmt <- function(e) gsub("\\b(maximum|minimum)_vec\\(", "\\1(", paste(format(e), collapse = ""))

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

.cvx_check_constraint_signs <- function(A, set, n, dirs, k) {
  sgn <- if (set$name == "Nonpositives") 1 else if (set$name == "Nonnegatives") -1 else return(invisible())
  for (j in seq_along(dirs)) for (i in seq_len(nrow(A))) {
    a <- A[i, n + j]
    if (a == 0) next
    if (!(sgn * a * dirs[j] >= 0))
      stop("Component ", i, " of constraint ", k, " in ", .moi_format_set(set), " has coefficient ", .jl_show_float(a),
           " on the ", if (dirs[j] > 0) "epigraph variable of a convex" else "hypograph variable of a concave",
           " atom, which relaxes ", if (set$name == "Nonpositives") "`g(u) <= 0`" else "`g(u) >= 0`",
           " instead of reformulating it: a convex atom may enter `<=` only nonnegatively and `>=` only ",
           "nonpositively (`norm(u) - t <= 0`, `t - norm(u) >= 0`), a concave atom exactly the other way ",
           "(`c - log(u) <= 0`, `log(u) - c >= 0`). Route to a general OptimizationProblem/NLP solver.",
           call. = FALSE)
  }
  invisible()
}

.cvx_dpp_extract <- function(prob, tr) {
  n <- length(prob$u0)
  params <- c(tr$params, tr$pdatavars)
  pkeys <- .keys(params); tkeys <- .keys(tr$taus)
  m <- length(tr$params)
  conA <- list(); conb0 <- list(); conB <- list()
  consets <- if (is.null(prob$constraints)) list() else lapply(prob$constraints, function(con) con$set)
  dirs <- vapply(tr$atoms, function(a) a$dir, 0L)
  for (k in seq_along(prob$constraints)) {
    con <- prob$constraints[[k]]
    blk <- .cvx_dpp_block(tr$consl[[k]], tr$allcols, params, pkeys, tkeys,
                          paste0("Constraint ", k, " (", .moi_format_set(con$set), ")"),
                          paste0("Constraint ", .moi_format_set(con$set), " is not affine in the variables."))
    .cvx_check_constraint_signs(blk$A, con$set, n, dirs, k)
    conA[[k]] <- blk$A; conb0[[k]] <- blk$b0; conB[[k]] <- blk$B
  }
  atomA <- list(); atomb0 <- list(); atomB <- list()
  for (j in seq_along(tr$atoms)) {
    at <- tr$atoms[[j]]
    blk <- .cvx_dpp_block(at$rows, tr$allcols, params, pkeys, tkeys,
                          paste0("The argument of atom ", j, " (", .moi_format_set(at$set), ")"),
                          "The argument of an atom must be affine in the optimization variables.")
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
  for (j in seq_along(tr$atoms)) {
    if (tr$atoms[[j]]$set$name %in% c("PositiveSemidefiniteConeTriangle", "LogDetConeTriangle") &&
        length(params) && any(cth$B[n + j, ] != 0))
      stop("A parameter may not scale a PSD or log-determinant atom's epigraph variable: its sign can change ",
           "across reinit! calls.", call. = FALSE)
  }
  dth <- .cvx_theta_affine(lo$b, params, "The objective's constant term")
  dpp <- list(n = n, m = m, c0 = cth$b0, C = cth$B, d0 = dth$b0[[1L]], dP = as.numeric(dth$B),
              params = params, pdataexprs = tr$pdataexprs,
              conA = conA, conb0 = conb0, conB = conB, consets = consets,
              atomA = atomA, atomb0 = atomb0, atomB = atomB,
              atomsets = lapply(tr$atoms, function(a) a$set), atomdirs = dirs,
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
  thetae <- .cvx_extend_theta(dpp, theta)
  np <- length(thetae)
  c0 <- dpp$c0 + if (np) as.numeric(dpp$C %*% thetae) else 0
  d <- dpp$d0 + if (np) .jl_dot(dpp$dP, thetae) else 0
  # an epigraph variable bounds its atom on one side only, so substituting it
  # is valid only where the objective pushes it against that bound; checked
  # before the model is touched so a rejected p leaves the cache usable
  sgn <- if (dpp$sense == "MaxSense") -1 else 1
  for (k in seq_along(dpp$atomdirs)) {
    dir <- dpp$atomdirs[k]
    if (!(sgn * dir * c0[n + k] >= 0))
      stop("Lowering an atom through its ", if (dir > 0) "epigraph" else "hypograph", " is valid only when the ",
           "objective is ", if (dir > 0) "nondecreasing" else "nonincreasing", " in it for MinSense (reversed for ",
           "MaxSense); got coefficient ", .jl_show_float(c0[n + k]), " for ", dpp$sense,
           if (dpp$m == 0) "" else paste0(" at p = ", .jl_show_vec(theta)),
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
  # user constraints in problem order, so conrefs stays 1:1 with
  # prob$constraints and with sol$dual
  conrefs <- lapply(seq_along(dpp$consets), function(k) {
    b <- dpp$conb0[[k]] + if (np) as.numeric(dpp$conB[[k]] %*% thetae) else 0
    .moi_add_constraint(model, dpp$conA[[k]], allx, b, dpp$consets[[k]])
  })
  atomrefs <- lapply(seq_along(dpp$atomsets), function(j) {
    b <- dpp$atomb0[[j]] + if (np) as.numeric(dpp$atomB[[j]] %*% thetae) else 0
    .moi_add_constraint(model, dpp$atomA[[j]], allx, b, dpp$atomsets[[j]])
  })
  coefs <- numeric(length(allx))
  nz <- c0 != 0
  coefs[nz] <- c0[nz]
  .moi_set_objective(model, coefs, d, if (dpp$sense == "MaxSense") "MAX_SENSE" else "MIN_SENSE")
  list(x = x, conrefs = conrefs, atomrefs = atomrefs)
}
