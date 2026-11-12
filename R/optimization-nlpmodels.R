# OptimizationNLPModels (lib/OptimizationNLPModels/src/OptimizationNLPModels.jl)
# and the parts of NLPModels.jl it reexports (src/nlp/meta.jl, counters.jl,
# tools.jl, utils.jl, api.jl and src/nls/meta.jl, counters.jl, api.jl):
# the model metadata NLPModelMeta and NLSMeta, the evaluation counters
# Counters and NLSCounters, models built from R functions
# (AbstractNLPModel, AbstractNLSModel), the conversion of a model to an
# OptimizationProblem, the adaptor NLPModelsAdaptor that exposes an
# optimization cache as a model, and build_nlpmodel_meta.
#
# Counters are environments, so the evaluation functions of the NLPModels
# module object update them in place as the Julia API does.

.nlp_dim_check <- function(n, ...) {
  vals <- list(...)
  for (nm in names(vals)) {
    if (length(vals[[nm]]) != n)
      .stopf("DimensionError: Input %s should have length %d not %d", nm, as.integer(n), length(vals[[nm]]))
  }
  invisible(TRUE)
}

#' NLPModels metadata
#'
#' `NLPModelMeta(nvar, ...)` describes the problem `optimize obj(x)` subject
#' to `lvar <= x <= uvar` and `lcon <= cons(x) <= ucon` with `nvar`
#' variables and `ncon` general constraints, as in NLPModels.jl. The
#' partitions of the variables (`ifix`, `ilow`, `iupp`, `irng`, `ifree`,
#' `iinf`) and of the constraints (`jfix`, `jlow`, `jupp`, `jrng`, `jfree`,
#' `jinf`) are computed from the bounds, and `nln` is the complement of the
#' linear constraints `lin`. `NLPModelMeta(meta, ...)` copies a metadata
#' object, replacing the fields given.
#'
#' `NLSMeta(nequ, nvar, ...)` describes the residual of a nonlinear least
#' squares model: `nequ` equations, `nnzj` and `nnzh` the sizes of the
#' sparse Jacobian and Hessian of the residual, `lin` the linear residuals.
#'
#' @param nvar Number of variables, or an `NLPModelMeta` to copy.
#' @param x0 Initial guess.
#' @param lvar,uvar Bounds of the variables.
#' @param nlvb,nlvo,nlvc Numbers of nonlinear variables in both, the
#'   objective and the constraints.
#' @param ncon Number of general constraints.
#' @param y0 Initial Lagrange multipliers.
#' @param lcon,ucon Bounds of the constraints.
#' @param nnzo Number of nonzeros of the gradient.
#' @param nnzj,lin_nnzj,nln_nnzj Numbers of entries of the sparse
#'   Jacobian: all, of the linear and of the nonlinear constraints.
#' @param nnzh Number of entries of the sparse Hessian of the Lagrangian.
#' @param lin Indices of the linear constraints (or residuals).
#' @param minimize `TRUE` to minimize, `FALSE` to maximize.
#' @param islp Whether the problem is a linear program.
#' @param name Problem name.
#' @param variable_bounds_analysis,constraint_bounds_analysis Whether to
#'   compute the partitions of the variables and constraints.
#' @param sparse_jacobian,sparse_hessian Whether the Jacobian and Hessian
#'   are sparse.
#' @param grad_available,jac_available,hess_available,jprod_available,jtprod_available,hprod_available
#'   Which derivatives the model provides.
#' @param nequ Number of residual equations.
#' @param jac_residual_available,hess_residual_available,jprod_residual_available,jtprod_residual_available,hprod_residual_available
#'   Which derivatives of the residual the model provides.
#' @return An `NLPModelMeta` or `NLSMeta` object (a list).
#' @examples
#' meta <- NLPModelMeta(3, lvar = c(0, -Inf, 1), uvar = c(0, 2, Inf), ncon = 2,
#'                      lcon = c(0, -Inf), ucon = c(0, 1))
#' meta$ifix; meta$iupp; meta$jfix; meta$nnzh
#' NLSMeta(4, 2)$nnzj
#' @export
NLPModelMeta <- function(nvar, x0 = NULL, lvar = NULL, uvar = NULL, nlvb = NULL, nlvo = NULL, nlvc = NULL,
                         ncon = 0L, y0 = NULL, lcon = NULL, ucon = NULL, nnzo = NULL, nnzj = NULL,
                         lin_nnzj = 0L, nln_nnzj = NULL, nnzh = NULL, lin = integer(), minimize = TRUE,
                         islp = FALSE, name = "Generic", variable_bounds_analysis = TRUE,
                         constraint_bounds_analysis = TRUE, sparse_jacobian = TRUE, sparse_hessian = TRUE,
                         grad_available = TRUE, jac_available = NULL, hess_available = TRUE,
                         jprod_available = NULL, jtprod_available = NULL, hprod_available = TRUE) {
  if (inherits(nvar, "NLPModelMeta")) {
    # NLPModelMeta(meta; kwargs...): the fields of meta, replaced by those given
    m <- nvar
    given <- as.list(match.call())[-1L]
    given <- given[setdiff(names(given), c("", "nvar"))]
    args <- m[c("x0", "lvar", "uvar", "nlvb", "nlvo", "nlvc", "ncon", "y0", "lcon", "ucon", "nnzo", "nnzj",
                "lin_nnzj", "nln_nnzj", "nnzh", "lin", "minimize", "islp", "name", "variable_bounds_analysis",
                "constraint_bounds_analysis", "sparse_jacobian", "sparse_hessian", "grad_available",
                "jac_available", "hess_available", "jprod_available", "jtprod_available", "hprod_available")]
    for (nm in names(given)) args[nm] <- list(get(nm, inherits = FALSE))
    return(do.call(NLPModelMeta, c(list(m$nvar), args)))
  }
  nvar <- as.integer(nvar); ncon <- as.integer(ncon)
  x0 <- as.numeric(x0 %||% numeric(nvar))
  lvar <- as.numeric(lvar %||% rep(-Inf, nvar)); uvar <- as.numeric(uvar %||% rep(Inf, nvar))
  y0 <- as.numeric(y0 %||% numeric(max(ncon, 0L)))
  lcon <- as.numeric(lcon %||% rep(-Inf, max(ncon, 0L))); ucon <- as.numeric(ucon %||% rep(Inf, max(ncon, 0L)))
  nnzj <- as.integer(nnzj %||% (nvar * ncon))
  nnzh <- as.integer(nnzh %||% ((nvar * (nvar + 1L)) %/% 2L))
  lin_nnzj <- as.integer(lin_nnzj)
  nln_nnzj <- as.integer(nln_nnzj %||% (nnzj - lin_nnzj))
  if (nvar < 1L || ncon < 0L || nnzj < 0L || nnzh < 0L) stop("Nonsensical dimensions", call. = FALSE)
  .nlp_dim_check(nvar, x0 = x0, lvar = lvar, uvar = uvar)
  .nlp_dim_check(ncon, y0 = y0, lcon = lcon, ucon = ucon)
  lin <- as.integer(lin)
  if (length(lin) && (any(lin < 1L) || any(lin > ncon)))
    stop("lin elements must be between 1 and ", ncon, call. = FALSE)
  if (nnzj != lin_nnzj + nln_nnzj) stop("AssertionError: nnzj == lin_nnzj + nln_nnzj", call. = FALSE)
  part <- function(lo, hi) list(fix = which(lo == hi), low = which(lo > -Inf & hi == Inf),
                                upp = which(lo == -Inf & hi < Inf), rng = which(lo > -Inf & hi < Inf & lo < hi),
                                free = which(lo == -Inf & hi == Inf), inf = which(lo > hi))
  none <- list(fix = integer(), low = integer(), upp = integer(), rng = integer(), free = integer(),
               inf = integer())
  vp <- if (isTRUE(variable_bounds_analysis)) part(lvar, uvar) else none
  cp <- if (ncon > 0L && isTRUE(constraint_bounds_analysis)) part(lcon, ucon) else none
  nln <- setdiff(seq_len(ncon), lin)
  structure(list(nvar = nvar, x0 = x0, lvar = lvar, uvar = uvar,
                 ifix = vp$fix, ilow = vp$low, iupp = vp$upp, irng = vp$rng, ifree = vp$free, iinf = vp$inf,
                 nlvb = as.integer(nlvb %||% nvar), nlvo = as.integer(nlvo %||% nvar),
                 nlvc = as.integer(nlvc %||% nvar), ncon = ncon, y0 = y0, lcon = lcon, ucon = ucon,
                 jfix = cp$fix, jlow = cp$low, jupp = cp$upp, jrng = cp$rng, jfree = cp$free, jinf = cp$inf,
                 nnzo = as.integer(nnzo %||% nvar), nnzj = nnzj, lin_nnzj = lin_nnzj, nln_nnzj = nln_nnzj,
                 nnzh = nnzh, nlin = length(lin), nnln = length(nln), lin = lin, nln = nln,
                 minimize = isTRUE(minimize), islp = isTRUE(islp), name = as.character(name),
                 variable_bounds_analysis = isTRUE(variable_bounds_analysis),
                 constraint_bounds_analysis = isTRUE(constraint_bounds_analysis),
                 sparse_jacobian = isTRUE(sparse_jacobian), sparse_hessian = isTRUE(sparse_hessian),
                 grad_available = isTRUE(grad_available), jac_available = isTRUE(jac_available %||% (ncon > 0L)),
                 hess_available = isTRUE(hess_available), jprod_available = isTRUE(jprod_available %||% (ncon > 0L)),
                 jtprod_available = isTRUE(jtprod_available %||% (ncon > 0L)),
                 hprod_available = isTRUE(hprod_available)),
            class = c("NLPModelMeta", "AbstractNLPModelMeta"))
}

#' @export
print.NLPModelMeta <- function(x, ...) {
  cat(sprintf("NLPModelMeta: %d variables (%d free, %d bounded, %d fixed), %d constraints (%d equalities), %s\n",
              x$nvar, length(x$ifree), length(x$ilow) + length(x$iupp) + length(x$irng), length(x$ifix),
              x$ncon, length(x$jfix), if (x$minimize) "minimize" else "maximize"))
  invisible(x)
}

#' @rdname NLPModelMeta
#' @export
NLSMeta <- function(nequ, nvar, x0 = NULL, nnzj = NULL, nnzh = NULL, lin = integer(),
                    jac_residual_available = TRUE, hess_residual_available = TRUE, jprod_residual_available = TRUE,
                    jtprod_residual_available = TRUE, hprod_residual_available = TRUE) {
  nequ <- as.integer(nequ); nvar <- as.integer(nvar)
  lin <- as.integer(lin)
  nln <- setdiff(seq_len(nequ), lin)
  structure(list(nequ = nequ, nvar = nvar, x0 = as.numeric(x0 %||% numeric(nvar)),
                 nnzj = max(0L, as.integer(nnzj %||% (nequ * nvar))),
                 nnzh = max(0L, as.integer(nnzh %||% ((nvar * (nvar + 1L)) %/% 2L))),
                 nln = nln, nnln = length(nln), lin = lin, nlin = length(lin),
                 jac_residual_available = isTRUE(jac_residual_available),
                 hess_residual_available = isTRUE(hess_residual_available),
                 jprod_residual_available = isTRUE(jprod_residual_available),
                 jtprod_residual_available = isTRUE(jtprod_residual_available),
                 hprod_residual_available = isTRUE(hprod_residual_available)),
            class = "NLSMeta")
}

#' @export
print.NLSMeta <- function(x, ...) {
  cat(sprintf("NLSMeta: %d residuals, %d variables, nnzj = %d, nnzh = %d\n", x$nequ, x$nvar, x$nnzj, x$nnzh))
  invisible(x)
}

.NLP_COUNTERS <- c("neval_obj", "neval_grad", "neval_cons", "neval_cons_lin", "neval_cons_nln", "neval_jcon",
                   "neval_jgrad", "neval_jac", "neval_jac_lin", "neval_jac_nln", "neval_jprod", "neval_jprod_lin",
                   "neval_jprod_nln", "neval_jtprod", "neval_jtprod_lin", "neval_jtprod_nln", "neval_hess",
                   "neval_hprod", "neval_jhess", "neval_jhprod")
.NLS_COUNTERS <- c("neval_residual", "neval_jac_residual", "neval_jprod_residual", "neval_jtprod_residual",
                   "neval_hess_residual", "neval_jhess_residual", "neval_hprod_residual")

#' Evaluation counters
#'
#' `Counters()` holds the numbers of evaluations of an NLPModels model
#' (`neval_obj`, `neval_grad`, `neval_cons`, `neval_jac`, `neval_jprod`,
#' `neval_jtprod`, `neval_hess`, `neval_hprod` and their linear, nonlinear
#' and per constraint variants), all zero. `NLSCounters()` adds the
#' evaluations of a least squares residual (`neval_residual`,
#' `neval_jac_residual`, ...) and forwards the fields of its `counters`.
#' Both are environments updated in place by the evaluation functions of
#' [NLPModels]; read a field with `c$neval_obj`.
#'
#' @return A `Counters` or `NLSCounters` environment.
#' @examples
#' cnt <- Counters()
#' cnt$neval_obj
#' NLPModels$sum_counters(cnt)
#' nc <- NLSCounters()
#' nc$neval_residual; nc$neval_obj
#' @export
Counters <- function() {
  e <- list2env(stats::setNames(as.list(integer(length(.NLP_COUNTERS))), .NLP_COUNTERS), parent = emptyenv())
  class(e) <- "Counters"
  e
}

#' @rdname Counters
#' @export
NLSCounters <- function() {
  e <- list2env(stats::setNames(as.list(integer(length(.NLS_COUNTERS))), .NLS_COUNTERS), parent = emptyenv())
  assign("counters", Counters(), envir = e)
  class(e) <- "NLSCounters"
  e
}

#' @export
`$.NLSCounters` <- function(x, name) {
  if (name %in% .NLP_COUNTERS) return(get(name, envir = get("counters", envir = x)))
  get(name, envir = x)
}

#' @export
`$<-.NLSCounters` <- function(x, name, value) {
  if (name %in% .NLP_COUNTERS) assign(name, value, envir = get("counters", envir = x)) else assign(name, value, envir = x)
  x
}

#' @export
print.Counters <- function(x, ...) {
  v <- vapply(.NLP_COUNTERS, function(n) get(n, envir = x), 0L)
  cat("Counters:\n")
  print(v[v != 0L])
  invisible(x)
}

#' @export
print.NLSCounters <- function(x, ...) {
  v <- vapply(.NLS_COUNTERS, function(n) get(n, envir = x), 0L)
  cat("NLSCounters:\n")
  print(v[v != 0L])
  print(get("counters", envir = x))
  invisible(x)
}

.nlp_counters <- function(nlp) {
  if (inherits(nlp, c("Counters", "NLSCounters"))) return(nlp)
  nlp$counters
}

.nlp_increment <- function(nlp, s, by = 1L) {
  cnt <- .nlp_counters(nlp)
  if (s %in% .NLP_COUNTERS && inherits(cnt, "NLSCounters")) cnt <- get("counters", envir = cnt)
  if (!exists(s, envir = cnt, inherits = FALSE)) .stopf("unknown counter %s", s)
  assign(s, get(s, envir = cnt) + as.integer(by), envir = cnt)
  invisible(nlp)
}

.nlp_sum_counters <- function(x) {
  cnt <- .nlp_counters(x)
  if (inherits(cnt, "NLSCounters")) {
    return(.nlp_sum_counters(get("counters", envir = cnt)) +
             sum(vapply(.NLS_COUNTERS, function(n) get(n, envir = cnt), 0L)))
  }
  keep <- setdiff(.NLP_COUNTERS, c("neval_cons", "neval_jac", "neval_jprod", "neval_jtprod"))
  sum(vapply(keep, function(n) get(n, envir = cnt), 0L))
}

.nlp_reset <- function(x) {
  cnt <- .nlp_counters(x)
  if (inherits(cnt, "NLSCounters")) {
    for (n in .NLS_COUNTERS) assign(n, 0L, envir = cnt)
    cnt <- get("counters", envir = cnt)
  }
  for (n in .NLP_COUNTERS) assign(n, 0L, envir = cnt)
  invisible(x)
}

#' Models of nonlinear programs
#'
#' `AbstractNLPModel(meta, obj, ...)` builds a model of the NLPModels.jl
#' interface from R functions of `x`: the objective `obj(x)`, its gradient
#' `grad(x)`, the constraints `cons(x)` (the vector bounded by
#' `meta$lcon` and `meta$ucon`), their Jacobian `jac(x)` (a dense matrix)
#' and the Hessian of the Lagrangian `hess(x, y, obj_weight)` (a dense
#' matrix of `obj_weight * H_f + sum(y * H_c)`). Missing derivatives are
#' approximated by central differences. Evaluate the model with the
#' functions of [NLPModels], which count the evaluations in `nlp$counters`
#' as NLPModels.jl does. `OptimizationProblem(nlp)` and
#' `OptimizationFunction(nlp)` convert the model, keeping its bounds,
#' constraint bounds, starting point and sense.
#'
#' `AbstractNLSModel(meta, nls_meta, residual, ...)` is a nonlinear least
#' squares model: the objective is `sum(residual(x)^2) / 2` and its gradient
#' `t(J) %*% residual(x)` with `J = jac_residual(x)`.
#'
#' @param meta An [NLPModelMeta()].
#' @param obj,grad,cons,jac,hess Functions of the model (see Details).
#' @param counters The evaluation counters.
#' @param nls_meta An [NLSMeta()].
#' @param residual,jac_residual The residual of a least squares model and
#'   its Jacobian.
#' @return An `AbstractNLPModel` or `AbstractNLSModel`.
#' @examples
#' meta <- NLPModelMeta(2, x0 = c(-1.2, 1))
#' nlp <- AbstractNLPModel(meta, obj = function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2)
#' NLPModels$obj(nlp, c(1, 1))
#' NLPModels$neval_obj(nlp)
#' sol <- solve(OptimizationProblem(nlp), LBFGSB())
#' round(sol$u, 4)
#' @export
AbstractNLPModel <- function(meta, obj, grad = NULL, cons = NULL, jac = NULL, hess = NULL, counters = Counters()) {
  if (!inherits(meta, "NLPModelMeta")) stop("`meta` must be an NLPModelMeta", call. = FALSE)
  structure(list(meta = meta, counters = counters, obj = obj, grad = grad, cons = cons, jac = jac, hess = hess),
            class = "AbstractNLPModel")
}

#' @rdname AbstractNLPModel
#' @export
AbstractNLSModel <- function(meta, nls_meta, residual, jac_residual = NULL, cons = NULL, jac = NULL,
                             counters = NLSCounters()) {
  if (!inherits(meta, "NLPModelMeta")) stop("`meta` must be an NLPModelMeta", call. = FALSE)
  if (!inherits(nls_meta, "NLSMeta")) stop("`nls_meta` must be an NLSMeta", call. = FALSE)
  structure(list(meta = meta, nls_meta = nls_meta, counters = counters, residual = residual,
                 jac_residual = jac_residual, cons = cons, jac = jac),
            class = c("AbstractNLSModel", "AbstractNLPModel"))
}

#' @export
print.AbstractNLPModel <- function(x, ...) {
  cat(class(x)[1L], "\n")
  print(x$meta)
  invisible(x)
}

# The evaluation API (NLPModels.jl src/nlp/api.jl): counters as the generic
# functions and the concrete models of the JSO ecosystem update them.
.nlp_obj <- function(nlp, x) {
  .nlp_dim_check(nlp$meta$nvar, x = x)
  if (inherits(nlp, "NLPModelsAdaptor")) return(as.numeric(nlp$cache$prob$f$f(x, nlp$cache$prob$p))[1L])
  .nlp_increment(nlp, "neval_obj")
  if (inherits(nlp, "AbstractNLSModel")) {
    Fx <- .nlp_residual(nlp, x)
    return(sum(Fx * Fx) / 2)
  }
  as.numeric(nlp$obj(x))
}

.nlp_grad <- function(nlp, x) {
  .nlp_dim_check(nlp$meta$nvar, x = x)
  if (inherits(nlp, "NLPModelsAdaptor")) {
    g <- nlp$cache$prob$f$grad
    p <- nlp$cache$prob$p
    if (is.function(g)) return(as.numeric(g(x, p)))
    return(.opt_fd_gradient(function(z) as.numeric(nlp$cache$prob$f$f(z, p))[1L], x))
  }
  if (inherits(nlp, "AbstractNLSModel")) {
    # grad!(nls, x, g) counts the gradient once and then calls
    # grad!(nls, x, g, Fx), which counts it again (src/nls/api.jl)
    .nlp_increment(nlp, "neval_grad", 2L)
    return(.nlp_jtprod_residual(nlp, x, .nlp_residual(nlp, x)))
  }
  .nlp_increment(nlp, "neval_grad")
  if (is.function(nlp$grad)) return(as.numeric(nlp$grad(x)))
  .opt_fd_gradient(function(z) as.numeric(nlp$obj(z)), x)
}

.nlp_lin_nln_increment <- function(nlp, what) {
  m <- nlp$meta
  if (m$nlin > 0L) .nlp_increment(nlp, paste0("neval_", what, "_lin"))
  if (m$nnln > 0L) .nlp_increment(nlp, paste0("neval_", what, "_nln"))
}

.nlp_cons <- function(nlp, x) {
  .nlp_dim_check(nlp$meta$nvar, x = x)
  if (inherits(nlp, "NLPModelsAdaptor")) {
    if (!nlp$meta$ncon) return(numeric())
    return(as.numeric(nlp$cache$prob$f$cons(x, nlp$cache$prob$p)))
  }
  .nlp_increment(nlp, "neval_cons")
  .nlp_lin_nln_increment(nlp, "cons")
  if (!nlp$meta$ncon) return(numeric())
  as.numeric(nlp$cons(x))
}

# the dense constraint Jacobian at x
.nlp_jac_dense <- function(nlp, x) {
  m <- nlp$meta
  if (inherits(nlp, "NLPModelsAdaptor")) {
    pr <- nlp$cache$prob
    if (is.function(pr$f$cons_j)) return(matrix(pr$f$cons_j(x, pr$p), m$ncon, m$nvar))
    return(.fd_jacobian(function(z) as.numeric(pr$f$cons(z, pr$p)), x))
  }
  if (is.function(nlp$jac)) return(matrix(nlp$jac(x), m$ncon, m$nvar))
  .fd_jacobian(function(z) as.numeric(nlp$cons(z)), x)
}

# dense structures: the Jacobian column by column, the lower triangle of the
# Hessian column by column (OptimizationNLPModels _enumerate_dense_structure
# and _enumerate_lower_triangle)
.nlp_dense_structure <- function(ncon, nvar) {
  list(rows = rep(seq_len(ncon), times = nvar), cols = rep(seq_len(nvar), each = ncon))
}

.nlp_lower_triangle <- function(n) {
  cols <- rep(seq_len(n), times = n:1)
  rows <- unlist(lapply(seq_len(n), function(j) j:n))
  list(rows = as.integer(rows), cols = as.integer(cols))
}

.nlp_jac_structure <- function(nlp) {
  m <- nlp$meta
  s <- .nlp_dense_structure(m$ncon, m$nvar)
  if (length(s$rows) != m$nnzj)
    .stopf("the model declares nnzj = %d; the dense Jacobian structure has %d entries", m$nnzj, length(s$rows))
  s
}

.nlp_jac_coord <- function(nlp, x) {
  .nlp_dim_check(nlp$meta$nvar, x = x)
  if (!inherits(nlp, "NLPModelsAdaptor")) {
    .nlp_increment(nlp, "neval_jac")
    .nlp_lin_nln_increment(nlp, "jac")
  }
  if (!nlp$meta$ncon) return(numeric())
  s <- .nlp_jac_structure(nlp)
  .nlp_jac_dense(nlp, x)[cbind(s$rows, s$cols)]
}

.nlp_jac <- function(nlp, x) {
  s <- .nlp_jac_structure(nlp)
  vals <- .nlp_jac_coord(nlp, x)
  J <- matrix(0, nlp$meta$ncon, nlp$meta$nvar)
  for (k in seq_along(vals)) J[s$rows[k], s$cols[k]] <- J[s$rows[k], s$cols[k]] + vals[k]
  J
}

.nlp_jprod <- function(nlp, x, v) {
  .nlp_dim_check(nlp$meta$nvar, x = x, v = v)
  if (!inherits(nlp, "NLPModelsAdaptor")) {
    .nlp_increment(nlp, "neval_jprod")
    .nlp_lin_nln_increment(nlp, "jprod")
  }
  if (!nlp$meta$ncon) return(numeric())
  as.numeric(.nlp_jac_dense(nlp, x) %*% v)
}

.nlp_jtprod <- function(nlp, x, v) {
  .nlp_dim_check(nlp$meta$nvar, x = x)
  .nlp_dim_check(nlp$meta$ncon, v = v)
  if (!inherits(nlp, "NLPModelsAdaptor")) {
    .nlp_increment(nlp, "neval_jtprod")
    .nlp_lin_nln_increment(nlp, "jtprod")
  }
  if (!nlp$meta$ncon) return(numeric(nlp$meta$nvar))
  as.numeric(crossprod(.nlp_jac_dense(nlp, x), v))
}

# the dense Hessian of the Lagrangian obj_weight * H_f + sum(y * H_c)
.nlp_hess_dense <- function(nlp, x, y, obj_weight) {
  m <- nlp$meta
  if (inherits(nlp, "NLPModelsAdaptor")) {
    pr <- nlp$cache$prob
    fn <- pr$f
    H <- if (is.function(fn$hess)) fn$hess(x, pr$p)
         else .fd_jacobian(function(z) .nlp_grad(nlp, z), x)
    H <- obj_weight * H
    if (length(y) && is.function(fn$cons_h)) {
      Hc <- fn$cons_h(x, pr$p)
      for (i in seq_along(y)) H <- H + y[i] * Hc[[i]]
    }
    return(H)
  }
  if (is.function(nlp$hess)) return(nlp$hess(x, y, obj_weight))
  if (inherits(nlp, "AbstractNLSModel")) {
    if (length(y) && any(y != 0)) stop("the model has no Hessian of its constraints", call. = FALSE)
    rf <- function(z) as.numeric(nlp$residual(z))
    jr <- function(z) if (is.function(nlp$jac_residual)) matrix(nlp$jac_residual(z), nlp$nls_meta$nequ, m$nvar)
                      else .fd_jacobian(rf, z)
    return(obj_weight * .fd_jacobian(function(z) as.numeric(crossprod(jr(z), rf(z))), x))
  }
  g <- if (is.function(nlp$grad)) function(z) as.numeric(nlp$grad(z))
       else function(z) .opt_fd_gradient(function(w) as.numeric(nlp$obj(w)), z)
  H <- obj_weight * .fd_jacobian(g, x)
  if (length(y) && any(y != 0)) {
    for (i in seq_along(y)) {
      ci <- function(z) as.numeric(nlp$cons(z))[i]
      H <- H + y[i] * .fd_jacobian(function(z) .opt_fd_gradient(ci, z), x)
    }
  }
  (H + t(H)) / 2
}

.nlp_hess_structure <- function(nlp) {
  s <- .nlp_lower_triangle(nlp$meta$nvar)
  if (length(s$rows) != nlp$meta$nnzh)
    .stopf("the model declares nnzh = %d; the dense lower triangle has %d entries", nlp$meta$nnzh, length(s$rows))
  s
}

.nlp_hess_coord <- function(nlp, x, y = numeric(nlp$meta$ncon), obj_weight = 1) {
  .nlp_dim_check(nlp$meta$nvar, x = x)
  .nlp_dim_check(nlp$meta$ncon, y = y)
  if (!inherits(nlp, "NLPModelsAdaptor")) .nlp_increment(nlp, "neval_hess")
  s <- .nlp_hess_structure(nlp)
  .nlp_hess_dense(nlp, x, y, obj_weight)[cbind(s$rows, s$cols)]
}

.nlp_hess <- function(nlp, x, y = numeric(nlp$meta$ncon), obj_weight = 1) {
  s <- .nlp_hess_structure(nlp)
  vals <- .nlp_hess_coord(nlp, x, y, obj_weight)
  n <- nlp$meta$nvar
  H <- matrix(0, n, n)
  for (k in seq_along(vals)) {
    H[s$rows[k], s$cols[k]] <- vals[k]
    H[s$cols[k], s$rows[k]] <- vals[k]
  }
  H
}

.nlp_hprod <- function(nlp, x, v, y = numeric(nlp$meta$ncon), obj_weight = 1) {
  .nlp_dim_check(nlp$meta$nvar, x = x, v = v)
  if (!inherits(nlp, "NLPModelsAdaptor")) .nlp_increment(nlp, "neval_hprod")
  as.numeric(.nlp_hess_dense(nlp, x, y, obj_weight) %*% v)
}

.nlp_residual <- function(nls, x) {
  .nlp_dim_check(nls$meta$nvar, x = x)
  .nlp_increment(nls, "neval_residual")
  as.numeric(nls$residual(x))
}

.nlp_jac_residual_dense <- function(nls, x) {
  if (is.function(nls$jac_residual)) return(matrix(nls$jac_residual(x), nls$nls_meta$nequ, nls$meta$nvar))
  .fd_jacobian(function(z) as.numeric(nls$residual(z)), x)
}

.nlp_jac_residual <- function(nls, x) {
  .nlp_dim_check(nls$meta$nvar, x = x)
  .nlp_increment(nls, "neval_jac_residual")
  .nlp_jac_residual_dense(nls, x)
}

.nlp_jprod_residual <- function(nls, x, v) {
  .nlp_increment(nls, "neval_jprod_residual")
  as.numeric(.nlp_jac_residual_dense(nls, x) %*% v)
}

.nlp_jtprod_residual <- function(nls, x, v) {
  .nlp_increment(nls, "neval_jtprod_residual")
  as.numeric(crossprod(.nlp_jac_residual_dense(nls, x), v))
}

.nlp_meta_of <- function(x) if (inherits(x, "AbstractNLPModelMeta")) x else x$meta

.nlp_has_bounds <- function(x) {
  m <- .nlp_meta_of(x)
  if (m$variable_bounds_analysis) length(m$ifree) < m$nvar else !all(is.infinite(m$lvar)) || !all(is.infinite(m$uvar))
}

.nlp_has_equalities <- function(x) {
  m <- .nlp_meta_of(x)
  if (m$constraint_bounds_analysis) m$ncon > 0L && length(m$jfix) > 0L else m$ncon > 0L && any(m$lcon == m$ucon)
}

.nlp_has_inequalities <- function(x) {
  m <- .nlp_meta_of(x)
  if (m$constraint_bounds_analysis) m$ncon > 0L && m$ncon > length(m$jfix) else m$ncon > 0L && any(m$lcon != m$ucon)
}

.nlp_coo_prod <- function(rows, cols, vals, v, n = max(rows)) {
  Av <- numeric(n)
  for (k in seq_along(rows)) Av[rows[k]] <- Av[rows[k]] + vals[k] * v[cols[k]]
  Av
}

.nlp_coo_sym_prod <- function(rows, cols, vals, v, n = length(v)) {
  Av <- numeric(n)
  for (k in seq_along(rows)) {
    i <- rows[k]; j <- cols[k]
    Av[i] <- Av[i] + vals[k] * v[j]
    if (i != j) Av[j] <- Av[j] + vals[k] * v[i]
  }
  Av
}

.nlp_counter_fns <- function() {
  fns <- lapply(c(.NLP_COUNTERS, .NLS_COUNTERS), function(n) {
    force(n)
    function(nlp) {
      cnt <- .nlp_counters(nlp)
      if (n %in% .NLP_COUNTERS && inherits(cnt, "NLSCounters")) cnt <- get("counters", envir = cnt)
      get(n, envir = cnt)
    }
  })
  stats::setNames(fns, c(.NLP_COUNTERS, .NLS_COUNTERS))
}

.nlp_meta_getters <- function() {
  flds <- names(NLPModelMeta(1L))
  stats::setNames(lapply(flds, function(f) {
    force(f)
    function(nlp) .nlp_meta_of(nlp)[[f]]
  }), paste0("get_", flds))
}

#' The NLPModels module
#'
#' `NLPModels` is the module binding that OptimizationNLPModels reexports:
#' a list holding the evaluation API of NLPModels.jl, reached as
#' `NLPModels$obj(nlp, x)`. It contains `obj`, `grad`, `objgrad`, `objcons`,
#' `cons`, `jac`, `jac_coord`, `jac_structure`, `jprod`, `jtprod`, `hess`
#' (dense), `hess_coord`, `hess_structure`, `hprod`, the least squares
#' functions `residual`, `jac_residual`, `jprod_residual`, `jtprod_residual`,
#' the counter readers `neval_obj`, `neval_grad`, ... and `sum_counters`,
#' `increment!`, `decrement!`, `reset!`, the predicates `has_bounds`,
#' `bound_constrained`, `unconstrained`, `linearly_constrained`,
#' `equality_constrained`, `inequality_constrained`, `has_equalities`,
#' `has_inequalities`, the metadata readers `get_nvar`, `get_ncon`, ...,
#' `nls_meta`, `coo_prod` and `coo_sym_prod`. Evaluations update the counters
#' of the model as in NLPModels.jl; the functions of an [NLPModelsAdaptor()]
#' do not count, as in OptimizationNLPModels.
#'
#' @format A named list of functions.
#' @examples
#' meta <- NLPModelMeta(2, ncon = 1, lcon = 0, ucon = Inf)
#' nlp <- AbstractNLPModel(meta, obj = function(x) sum(x^2), cons = function(x) x[1] + x[2] - 1)
#' NLPModels$cons(nlp, c(1, 2))
#' NLPModels$jac(nlp, c(1, 2))
#' NLPModels$neval_cons(nlp)
#' NLPModels$has_inequalities(nlp)
#' @export
NLPModels <- c(list(
  obj = function(nlp, x) .nlp_obj(nlp, x),
  grad = function(nlp, x) .nlp_grad(nlp, x),
  objgrad = function(nlp, x) list(.nlp_obj(nlp, x), .nlp_grad(nlp, x)),
  objcons = function(nlp, x) list(.nlp_obj(nlp, x), .nlp_cons(nlp, x)),
  cons = function(nlp, x) .nlp_cons(nlp, x),
  jac = function(nlp, x) .nlp_jac(nlp, x),
  jac_coord = function(nlp, x) .nlp_jac_coord(nlp, x),
  jac_structure = function(nlp) .nlp_jac_structure(nlp),
  jprod = function(nlp, x, v) .nlp_jprod(nlp, x, v),
  jtprod = function(nlp, x, v) .nlp_jtprod(nlp, x, v),
  hess = function(nlp, x, y = numeric(nlp$meta$ncon), obj_weight = 1) .nlp_hess(nlp, x, y, obj_weight),
  hess_coord = function(nlp, x, y = numeric(nlp$meta$ncon), obj_weight = 1) .nlp_hess_coord(nlp, x, y, obj_weight),
  hess_structure = function(nlp) .nlp_hess_structure(nlp),
  hprod = function(nlp, x, v, y = numeric(nlp$meta$ncon), obj_weight = 1) .nlp_hprod(nlp, x, v, y, obj_weight),
  residual = function(nls, x) .nlp_residual(nls, x),
  jac_residual = function(nls, x) .nlp_jac_residual(nls, x),
  jprod_residual = function(nls, x, v) .nlp_jprod_residual(nls, x, v),
  jtprod_residual = function(nls, x, v) .nlp_jtprod_residual(nls, x, v),
  nls_meta = function(nls) nls$nls_meta,
  get_counters = function(nlp) .nlp_counters(nlp),
  sum_counters = function(x) .nlp_sum_counters(x),
  `increment!` = function(nlp, s) .nlp_increment(nlp, s, 1L),
  `decrement!` = function(nlp, s) .nlp_increment(nlp, s, -1L),
  `reset!` = function(x) .nlp_reset(x),
  reset_data = function(nlp) nlp,
  has_bounds = function(nlp) .nlp_has_bounds(nlp),
  bound_constrained = function(nlp) .nlp_meta_of(nlp)$ncon == 0L && .nlp_has_bounds(nlp),
  unconstrained = function(nlp) .nlp_meta_of(nlp)$ncon == 0L && !.nlp_has_bounds(nlp),
  linearly_constrained = function(nlp) {
    m <- .nlp_meta_of(nlp)
    m$ncon > 0L && m$nlin == m$ncon
  },
  equality_constrained = function(nlp) {
    m <- .nlp_meta_of(nlp)
    if (m$constraint_bounds_analysis) m$ncon > 0L && length(m$jfix) == m$ncon else m$ncon > 0L && all(m$lcon == m$ucon)
  },
  inequality_constrained = function(nlp) {
    m <- .nlp_meta_of(nlp)
    if (m$constraint_bounds_analysis) m$ncon > 0L && !length(m$jfix) else m$ncon > 0L && all(m$lcon != m$ucon)
  },
  has_equalities = function(nlp) .nlp_has_equalities(nlp),
  has_inequalities = function(nlp) .nlp_has_inequalities(nlp),
  coo_prod = function(rows, cols, vals, v, n = max(rows)) .nlp_coo_prod(rows, cols, vals, v, n),
  coo_sym_prod = function(rows, cols, vals, v, n = length(v)) .nlp_coo_sym_prod(rows, cols, vals, v, n)
), .nlp_counter_fns(), .nlp_meta_getters())

# ---------------------------------------------------------------------------
# Conversion to OptimizationFunction and OptimizationProblem

.nlp_optimization_function <- function(nlp, ...) {
  f <- function(u, p) .nlp_obj(nlp, u)
  grad <- function(u, p) .nlp_grad(nlp, u)
  hess <- function(u, p) .nlp_hess(nlp, u)
  m <- nlp$meta
  if (!NLPModels$unconstrained(nlp) && !NLPModels$bound_constrained(nlp)) {
    fn <- OptimizationFunction(f, grad = grad, hess = hess, cons = function(u, p) .nlp_cons(nlp, u),
                               cons_j = function(u, p) .nlp_jac(nlp, u), ...)
  } else {
    fn <- OptimizationFunction(f, grad = grad, hess = hess, ...)
  }
  fn$hv <- function(u, v, p) .nlp_hprod(nlp, u, v)
  fn$nlpmodel <- nlp
  fn
}

.nlp_optimization_problem <- function(nlp, ..., adtype = NULL) {
  dots <- list(...)
  # a positional AD choice, as OptimizationProblem(nlpmodel, AutoForwardDiff()),
  # is accepted and unused: the derivatives come from the model
  if (length(dots) && (is.null(names(dots)) || !nzchar(names(dots)[1L]))) dots <- dots[-1L]
  f <- .nlp_optimization_function(nlp)
  m <- nlp$meta
  lb <- ub <- NULL
  if (.nlp_has_bounds(nlp)) { lb <- m$lvar; ub <- m$uvar }
  lcons <- ucons <- NULL
  if (.nlp_has_inequalities(nlp) || .nlp_has_equalities(nlp)) { lcons <- m$lcon; ucons <- m$ucon }
  sense <- if (m$minimize) MinSense else MaxSense
  do.call(OptimizationProblem, c(list(f, m$x0, lb = lb, ub = ub, lcons = lcons, ucons = ucons, sense = sense),
                                 dots))
}

# ---------------------------------------------------------------------------
# The adaptor and the metadata of a cache

.nlp_cache_prob <- function(cache) {
  if (inherits(cache, "OptimizationCache")) return(cache$prob)
  if (inherits(cache, "OptimizationProblem")) return(cache)
  stop("expected an OptimizationCache or an OptimizationProblem", call. = FALSE)
}

#' Optimization caches as NLPModels models
#'
#' `build_nlpmodel_meta(cache)` returns the [NLPModelMeta()] of an
#' optimization cache (from [init()]) or problem: its variables and
#' starting point, the bounds of the variables (infinite when absent), the
#' bounds of the constraints, zero multipliers, the sense, and the sizes of
#' the dense Jacobian (`ncon * nvar`) and of the lower triangle of the
#' Hessian (`nvar * (nvar + 1) / 2`).
#'
#' `NLPModelsAdaptor(cache, meta, counters)` exposes the cache as a model of
#' the NLPModels interface: [NLPModels] functions evaluate the objective,
#' gradient, constraints, dense Jacobian (column by column), lower triangle
#' of the Hessian of the Lagrangian and the Jacobian products of the
#' problem, as OptimizationNLPModels does for the solvers that take an
#' NLPModel.
#'
#' @param cache An `OptimizationCache` or `OptimizationProblem`.
#' @param meta An `NLPModelMeta`; `build_nlpmodel_meta(cache)` when `NULL`.
#' @param counters A [Counters()] object.
#' @return An `NLPModelMeta` or an `NLPModelsAdaptor`.
#' @examples
#' prob <- OptimizationProblem(function(u, p) sum((u - p)^2), c(0, 0), c(1, 2),
#'                             lb = c(-5, -5), ub = c(5, 5))
#' cache <- init(prob, LBFGSB())
#' meta <- build_nlpmodel_meta(cache)
#' meta$nnzh; meta$irng
#' nlp <- NLPModelsAdaptor(cache, meta, Counters())
#' NLPModels$grad(nlp, c(0, 0))
#' NLPModels$hess_coord(nlp, c(0, 0))
#' @export
build_nlpmodel_meta <- function(cache) {
  prob <- .nlp_cache_prob(cache)
  u0 <- as.numeric(prob$u0)
  nvar <- length(u0)
  ncon <- length(prob$lcons %||% numeric())
  fn <- prob$f
  nnzj <- if (is.null(fn$cons_jac_prototype)) (if (ncon > 0L) ncon * nvar else 0L) else length(fn$cons_jac_prototype)
  hp <- if (ncon > 0L) fn$lag_hess_prototype else fn$hess_prototype
  nnzh <- if (is.null(hp)) (nvar * (nvar + 1L)) %/% 2L else (nrow(hp) * (nrow(hp) + 1L)) %/% 2L
  NLPModelMeta(nvar, ncon = ncon, nnzj = nnzj, nnzh = nnzh, x0 = u0, y0 = numeric(ncon),
               lvar = prob$lb %||% rep(-Inf, nvar), uvar = prob$ub %||% rep(Inf, nvar),
               lcon = prob$lcons %||% numeric(), ucon = prob$ucons %||% numeric(),
               minimize = !identical(unclass(prob$sense %||% MinSense), "MaxSense"))
}

#' @rdname build_nlpmodel_meta
#' @export
NLPModelsAdaptor <- function(cache, meta = NULL, counters = Counters()) {
  if (inherits(cache, "OptimizationProblem")) cache <- init(cache)
  meta <- meta %||% build_nlpmodel_meta(cache)
  js <- .nlp_dense_structure(meta$ncon, meta$nvar)
  hs <- .nlp_lower_triangle(meta$nvar)
  structure(list(cache = cache, meta = meta, counters = counters, jac_rows = js$rows, jac_cols = js$cols,
                 hess_rows = hs$rows, hess_cols = hs$cols),
            class = c("NLPModelsAdaptor", "AbstractNLPModel"))
}
