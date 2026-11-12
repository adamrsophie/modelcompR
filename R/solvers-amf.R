# Approximate matrix factorization (AMF) for the Rosenbrock methods, after
# OrdinaryDiffEqAMF 2.3.1 (amf.jl, linsolve.jl, operators.jl): the
# iteration matrix W = J - I / (dt gamma) of a Rosenbrock step is the
# operator W_prototype of the ODEFunction, here
#   -(I - gamma J_1)(I - gamma J_2) ... (I - gamma J_k) / gamma
# for a Jacobian split into J_1 + ... + J_k, and SciMLOpFactorization
# solves with it factor by factor (factorize of SciMLOperators: a scaled or
# composed operator factorizes part by part, a tensor product through its
# factors, any other operator through the structured factorize of its
# matrix).
#
# OrdinaryDiffEqAMF needs an in place right-hand side (the out of place
# Rosenbrock caches do not take an operator W), so the steps follow the in
# place perform_step! of Rosenbrock23Cache and RosenbrockCache
# (rosenbrock_perform_step.jl, @muladd: fused multiply adds) with the W
# operator branch of calc_W! and calc_rosenbrock_differentiation!
# (OrdinaryDiffEqDifferentiation, derivative_utils.jl): W is updated at
# (uprev, p, t; gamma = dt gamma) and factorized on every attempt (nw), no
# Jacobian is evaluated (njacs 0), and the time derivative is kept after a
# rejected step of a Rodas method (do_newJW) and fresh on every attempt of
# Rosenbrock23 (a W method with an operator W, _rosenbrock_jac_reuse_decision).

#' Approximate matrix factorization
#'
#' Rosenbrock methods whose linear systems use an approximate matrix
#' factorization (AMF) of the iteration matrix, as OrdinaryDiffEqAMF. For a
#' Jacobian split as `J = J_1 + ... + J_k` (for example into the parts of
#' the directions of a PDE discretization, or into triangular parts) the
#' iteration matrix `W = -(I - gamma J) / gamma` of a step is replaced by
#' `-(I - gamma J_1) ... (I - gamma J_k) / gamma`, whose factors are cheap to
#' solve with.
#'
#' * `build_amf_function(f, jac, split, amf_factors)` returns the
#'   [ODEFunction()] of `f(u, p, t)` with the Jacobian operator `jac` as
#'   `jac_prototype`, the AMF operator as `W_prototype` and the sparsity
#'   pattern (the matrix of `jac` unless `sparsity` is given). Give `split`,
#'   a list of operators `J_1, ..., J_k`, for the AMF product; `amf_factors`,
#'   a list of prebuilt factors (each `I - gamma J_i` with a `gamma` scalar
#'   operator updated by the solver); or neither for the exact `W`.
#' * `AMFOperator(jac)` is the exact operator `-(I - gamma J) / gamma` of an
#'   operator `jac`; `AMFOperator(list(J_1, ..., J_k))` the AMF product of
#'   split operators; `AMFOperator(factors = list(F_1, ..., F_k))` the
#'   operator `-F_1 ... F_k / gamma` of prebuilt factors. The `gamma` of
#'   these operators is a [ScalarOperator()] set by
#'   `update_coefficients(W, u, p, t, gamma = dtgamma)`.
#' * `AMF(alg, ...)` builds the Rosenbrock method `alg(...)`
#'   ([Rosenbrock23()], [Rodas4()] or [Rodas5P()]; `...` such as `autodiff`
#'   go to its constructor) with `SciMLOpFactorization()` as its linear
#'   solver. The problem must come from `build_amf_function()` (an
#'   `ODEFunction` with a `W_prototype` operator). No Jacobian is evaluated;
#'   the time derivative is a forward difference unless `f` has `tgrad`.
#' * `SciMLOpFactorization()` is the linear solver for a [LinearProblem()]
#'   whose matrix is a SciML operator: a product or a multiple of operators is
#'   factorized part by part, a tensor product through its factors, and any
#'   other part through the structured factorization of its matrix (as
#'   `factorize` of Julia's LinearAlgebra: diagonal, bidiagonal or triangular
#'   substitution, LU of a tridiagonal matrix, Bunch-Kaufman of a symmetric
#'   matrix, else LU); the factorization is kept while the operator does not
#'   change.
#'
#' As in OrdinaryDiffEq the steps are those of the in place methods (an
#' AMF solve in Julia needs an in place `f!`): the Rodas methods evaluate
#' `f` at the start of every step, and a fixed step `Rosenbrock23`
#' evaluates the derivative at the new state at the start of the next step.
#'
#' @param alg A Rosenbrock method constructor (`Rosenbrock23`, `Rodas4`,
#'   `Rodas5P`) or a Rosenbrock algorithm object.
#' @param f The right-hand side `f(u, p, t)`.
#' @param jac A square SciML operator: the Jacobian of `f` (for example
#'   `MatrixOperator(J0, update_func)` or a sum of such operators).
#' @param split,amf_factors Lists of square operators of the size of `jac`:
#'   the parts `J_i` of the Jacobian or prebuilt factors `I - gamma J_i`.
#'   When both are given they must have the same length; `amf_factors` then
#'   builds the operator.
#' @param factors Prebuilt factors for `AMFOperator()`.
#' @param jac_cache,w_cache Prototype vectors for [cache_operator()] of the
#'   Jacobian and the AMF operator (zero vectors by default).
#' @param sparsity Sparsity pattern; by default the matrix of `jac`.
#' @param ... Arguments of the Rosenbrock constructor (`autodiff`).
#' @return An algorithm, an operator, a linear solver or an `ODEFunction`.
#' @examples
#' A1 <- rbind(c(-2, 1), c(0, -3))
#' A2 <- rbind(c(0, 0), c(1, 0))
#' f <- function(u, p, t) c(-2 * u[1] + u[2], u[1] - 3 * u[2])
#' J1 <- MatrixOperator(A1)
#' J2 <- MatrixOperator(A2)
#' fun <- build_amf_function(f, jac = J1 + J2, split = list(J1, J2))
#' prob <- ODEProblem(fun, c(1, 1), c(0, 1))
#' sol <- solve(prob, AMF(Rodas5P))
#' sol$u[[length(sol$u)]]
#' W <- update_coefficients(fun$W_prototype, gamma = 0.1)
#' solve(LinearProblem(W, c(1, 2)), SciMLOpFactorization())$u
#' @name amf
NULL

#' @rdname amf
#' @export
AMF <- function(alg, ...) {
  inner <- if (inherits(alg, "ODEAlgorithm")) alg else if (is.function(alg)) alg(...) else
    stop("`alg` must be a Rosenbrock method such as Rodas5P", call. = FALSE)
  if (!inherits(inner, "ODEAlgorithm")) stop("`alg` must be a Rosenbrock method such as Rodas5P", call. = FALSE)
  if (identical(inner$name, "Rosenbrock23")) {
    step <- .amf_ros23_step
    init <- .amf_ros23_init
    on_event <- .amf_ros23_on_event
  } else {
    tab <- if (is.function(inner$step)) environment(inner$step)$tab
    if (!is.list(tab) || is.null(tab$H) || is.null(tab$btilde))
      stop("AMF() wraps a Rosenbrock method: Rosenbrock23, Rodas4 or Rodas5P", call. = FALSE)
    step <- .amf_rodas_step(tab)
    init <- .amf_rodas_init
    on_event <- function(ctx, t, u, cache) list(cache = cache, nf = 0L)
  }
  inner$step <- step
  inner$init <- init
  inner$on_event <- on_event
  inner$linsolve <- SciMLOpFactorization()
  # W_prototype is -(I - gamma J) / gamma: no mass matrix
  inner$mass_matrix <- FALSE
  inner$amf <- TRUE
  inner
}

#' @rdname amf
#' @export
SciMLOpFactorization <- function() .lin_alg("SciMLOpFactorization", "sciml_op")

# solve! of SciMLOpFactorization: factorize when the operator changed
# (isfresh) and apply the factorization.
.amf_lin_solve <- function(A, b, alg, cache) {
  A <- .as_sop(A)
  if (is.null(cache$amf_fact) || !identical(cache$amf_A, A)) {
    cache$amf_fact <- .amf_factorize(A, length(b))
    cache$amf_A <- A
  }
  if (is.null(cache$amf_fact)) stop("SciMLOpFactorization: the operator is singular", call. = FALSE)
  u <- cache$amf_fact(b)
  list(u = u, resid = as.numeric(A$apply(u)) - b, iters = 1L, alg = alg)
}

# factorize of SciMLOperators: a multiple (ScaledOperator) factorizes as the
# multiple of the factorization, a product (ComposedOperator) as the product
# of the factorizations, an InvertibleOperator keeps its own, and any other
# operator is factorized from its matrix. Returns the solve b -> W^-1 b
# (the factor inverses in order, then the division by the scalar), or NULL
# when a factor is singular or not finite.
.amf_factorize <- function(W, n) {
  parts <- .amf_parts(W)
  facs <- lapply(parts$ops, .amf_factor, n = n)
  if (any(vapply(facs, is.null, TRUE)) || !is.finite(parts$lambda) || parts$lambda == 0) return(NULL)
  lambda <- parts$lambda
  function(b) {
    for (fac in facs) b <- fac(b)
    b / lambda
  }
}

.amf_parts <- function(L) {
  if (inherits(L, "ScalarOperator")) return(list(lambda = L$value, ops = list()))
  if (inherits(L, "ComposedOperator")) {
    a <- .amf_parts(L$ops[[1L]])
    b <- .amf_parts(L$ops[[2L]])
    return(list(lambda = a$lambda * b$lambda, ops = c(a$ops, b$ops)))
  }
  list(lambda = 1, ops = list(L))
}

# A part of the operator: a TensorProductOperator factorizes its factors
# (factorize of tensor.jl), any other operator is factorized from its
# matrix by factorize(::AbstractMatrix) of LinearAlgebra: Diagonal, lower
# or upper Bidiagonal and Triangular (substitution), the LU of a tridiagonal
# matrix, Bunch-Kaufman for a symmetric matrix, else LU.
.amf_factor <- function(L, n) {
  if (inherits(L, "TensorProductOperator")) {
    A <- L$A
    B <- L$B
    A$solve_fn <- .amf_factorize(A, A$size[1L])
    B$solve_fn <- .amf_factorize(B, B$size[1L])
    if (is.null(A$solve_fn) || is.null(B$solve_fn)) return(NULL)
    return(TensorProductOperator(A, B)$solve_fn)
  }
  M <- .sop_matrix(L, n)
  if (!all(is.finite(M))) return(NULL)
  # zeros above the first superdiagonal (lower Hessenberg) and below the
  # first subdiagonal (upper Hessenberg)
  lo1 <- all(M[col(M) > row(M) + 1L] == 0)
  up1 <- all(M[row(M) > col(M) + 1L] == 0)
  lo <- lo1 && all(M[col(M) == row(M) + 1L] == 0)
  up <- up1 && all(M[row(M) == col(M) + 1L] == 0)
  d <- diag(M)
  if ((lo || up) && any(d == 0)) return(NULL)
  if (lo && up) return(function(b) b / d)
  if (lo) return(function(b) as.numeric(forwardsolve(M, b)))
  if (up) return(function(b) as.numeric(backsolve(M, b)))
  if (qr(M)$rank < nrow(M)) return(NULL)
  if (!(lo1 && up1) && all(M == t(M))) {
    bk <- Matrix::BunchKaufman(M)
    return(function(b) as.numeric(Matrix::solve(bk, b)))
  }
  function(b) as.numeric(solve(M, b))
}

#' @rdname amf
#' @export
AMFOperator <- function(jac = NULL, factors = NULL) {
  if (!is.null(factors)) {
    if (!is.null(jac)) stop("give `jac` (or a list of split operators) or `factors`, not both", call. = FALSE)
    fs <- .amf_normalize(factors)
    .amf_validate(fs)
    return(.amf_assemble(fs))
  }
  if (inherits(jac, "AbstractSciMLOperator")) {
    n <- .amf_square(jac, "Jacobian operators must be square")
    return(-(IdentityOperator(n) - .amf_gamma() * jac) * .amf_transform())
  }
  if (is.list(jac)) {
    ops <- .amf_normalize(jac)
    .amf_validate(ops)
    return(.amf_assemble(lapply(ops, .amf_factor_operator)))
  }
  stop("`jac` must be an AbstractSciMLOperator or a list of split operators", call. = FALSE)
}

# The gamma of the factors and its inverse, updated through
# update_coefficients(W, u, p, t, gamma = dtgamma).
.amf_gamma <- function() ScalarOperator(1, function(alpha, u, p, t, gamma = 1) gamma)
.amf_transform <- function() ScalarOperator(0, function(alpha, u, p, t, gamma = 1) 1 / gamma)

.amf_square <- function(L, msg) {
  s <- L$size
  if (anyNA(s) || s[1L] != s[2L]) stop(msg, call. = FALSE)
  s[1L]
}

.amf_normalize <- function(ops) {
  ops <- lapply(if (inherits(ops, "AbstractSciMLOperator")) list(ops) else as.list(ops), .as_sop)
  if (!length(ops)) stop("split must contain at least one operator", call. = FALSE)
  ops
}

.amf_validate <- function(ops) {
  n <- .amf_square(ops[[1L]], "split operators must be square")
  for (op in ops[-1L]) {
    s <- op$size
    if (anyNA(s) || s[1L] != n || s[2L] != n)
      stop("all split operators must have the same square size", call. = FALSE)
  }
  n
}

.amf_factor_operator <- function(J) {
  n <- .amf_square(J, "Jacobian operators must be square")
  IdentityOperator(n) - .amf_gamma() * J
}

.amf_assemble <- function(factors) -(Reduce(`*`, factors)) * .amf_transform()

#' @rdname amf
#' @export
build_amf_function <- function(f, jac, split = NULL, amf_factors = NULL, jac_cache = NULL, w_cache = NULL,
                               sparsity = NULL) {
  if (!is.function(f)) stop("`f` must be a function f(u, p, t)", call. = FALSE)
  if (!inherits(jac, "AbstractSciMLOperator")) stop("jac must be an AbstractSciMLOperator", call. = FALSE)
  n <- .amf_square(jac, "jac must be square")
  J_op <- cache_operator(jac, jac_cache %||% numeric(n))
  if (!is.null(split) && !is.null(amf_factors) &&
      length(.amf_normalize(split)) != length(.amf_normalize(amf_factors)))
    stop("split and amf_factors must contain the same number of operators", call. = FALSE)
  W_base <- if (!is.null(amf_factors)) {
    fs <- .amf_normalize(amf_factors)
    if (.amf_validate(fs) != n) stop("amf_factors must have the same size as jac", call. = FALSE)
    AMFOperator(factors = fs)
  } else if (is.null(split)) {
    AMFOperator(J_op)
  } else {
    ops <- .amf_normalize(split)
    if (.amf_validate(ops) != n) stop("split operators must have the same size as jac", call. = FALSE)
    AMFOperator(ops)
  }
  W_op <- cache_operator(W_base, w_cache %||% numeric(n))
  fun <- ODEFunction(f)
  fun$jac_prototype <- J_op
  fun$W_prototype <- W_op
  fun$sparsity <- sparsity %||% concretize(J_op)
  fun
}

# ---------------------------------------------------------------------------
# Steps

# calc_W! with an operator W_prototype: update_coefficients!(W, uprev, p, t;
# gamma = dtgamma), then the factorization of SciMLOpFactorization.
.amf_W <- function(ctx, t, u, dtgamma) {
  W <- ctx$fobj$W_prototype
  if (!inherits(W, "AbstractSciMLOperator"))
    stop("AMF() needs an ODEFunction from build_amf_function() (a W_prototype operator)", call. = FALSE)
  .amf_factorize(update_coefficients(W, u, ctx$p, t, gamma = dtgamma), length(u))
}

.amf_ros23_init <- function(ctx, t, u) {
  cache <- new.env(parent = emptyenv())
  cache$fsal <- ctx$f(u, ctx$p, t)
  cache$stale <- FALSE
  cache$nf0 <- 1L
  cache
}

# reset_fsal! after an event
.amf_ros23_on_event <- function(ctx, t, u, cache) {
  cache$fsal <- ctx$f(u, ctx$p, t)
  cache$stale <- FALSE
  list(cache = cache, nf = 1L)
}

# perform_step! of Rosenbrock23Cache. A fixed step solve does not evaluate
# fsallast in the step; update_fsal! evaluates f at the new state before the
# next step (only_diagonal_mass_matrix and not adaptive), done here at the
# start of that step.
.amf_ros23_step <- function(ctx, t, u, dt, cache) {
  f <- ctx$f; p <- ctx$p
  c32 <- .ROS23$c32; d <- .ROS23$d
  dtg <- dt * d
  neginvdtg <- -(1 / dtg)
  dto2 <- dt / 2
  dto6 <- dt / 6
  nf <- 0L
  if (cache$stale) {
    cache$fsal <- f(u, p, t)
    cache$stale <- FALSE
    nf <- 1L
  }
  fsalfirst <- cache$fsal
  Wsolve <- .amf_W(ctx, t, u, dtg)
  tg <- .calc_tgrad(ctx, u, t, fsalfirst)
  dT <- tg$dT
  out <- list(nf = nf + tg$nf, njacs = 0L, nw = 1L, nsolve = 0L, cache = cache)
  if (is.null(Wsolve)) return(c(out, list(u = u, fail = TRUE)))
  k1 <- Wsolve(fsalfirst + dtg * dT) * neginvdtg
  f1 <- f(.fma(dto2, k1, u), p, t + dto2)
  k2 <- .fma(Wsolve(f1 - k1), neginvdtg, k1)
  unew <- .fma(dt, k2, u)
  out$nf <- out$nf + 1L
  out$nsolve <- 2L
  EEst <- 0
  if (ctx$adaptive) {
    fsallast <- f(unew, p, t + dt)
    out$nf <- out$nf + 1L
    lt <- .fma(-c32, k2 - f1, fsallast)
    lt <- .fma(-2, k1 - fsalfirst, lt)
    lt <- .fma(dt, dT, lt)
    k3 <- Wsolve(lt) * neginvdtg
    out$nsolve <- 3L
    EEst <- .err_norm(dto6 * (k1 - 2 * k2 + k3), u, unew, ctx$abstol, ctx$reltol)
    if (isTRUE(EEst <= 1)) cache$fsal <- fsallast
  } else {
    cache$stale <- TRUE
  }
  c(out, list(u = unew, EEst = EEst, k = list(k1, k2)))
}

.amf_rodas_init <- function(ctx, t, u) {
  cache <- new.env(parent = emptyenv())
  cache$dT <- NULL
  cache$nf0 <- 0L
  cache
}

# sum of w[i] * ks[[i]] over the nonzero weights by fused multiply adds in
# stage order, starting from acc (_madd_terms).
.amf_madd <- function(acc, w, ks) {
  for (i in seq_along(w)) if (w[i] != 0) acc <- .fma(w[i], ks[[i]], acc)
  acc
}

# perform_step! of RosenbrockCache (the Rodas methods): not FSAL, f is
# evaluated at the start of every step; the time derivative is reused after
# a rejected step (do_newJW: new_jac false, calc_tderivative! skipped).
.amf_rodas_step <- function(tab) {
  force(tab)
  function(ctx, t, u, dt, cache) {
    f <- ctx$f; p <- ctx$p
    A <- tab$A; cc <- tab$c
    s <- nrow(A)
    dtC <- tab$C * (1 / dt)
    dtd <- dt * tab$d
    dtgamma <- dt * tab$gamma
    du <- f(u, p, t)
    nf <- 1L
    Wsolve <- .amf_W(ctx, t, u, dtgamma)
    reuse <- !is.null(cache$dT) && ctx$adaptive && ctx$iter > 1L && ctx$EEst_prev > 1 &&
      !isTRUE(ctx$derivative_discontinuity) && !isTRUE(ctx$composite)
    if (!reuse) {
      tg <- .calc_tgrad(ctx, u, t, du)
      cache$dT <- tg$dT
      nf <- nf + tg$nf
    }
    dT <- cache$dT
    out <- list(nf = nf, njacs = 0L, nw = 1L, nsolve = 0L, cache = cache)
    if (is.null(Wsolve)) return(c(out, list(u = u, fail = TRUE)))
    ks <- vector("list", s)
    ks[[1L]] <- -Wsolve(du + dtd[1L] * dT)
    for (stage in 2:s) {
      prev <- seq_len(stage - 1L)
      ui <- .amf_madd(u, A[stage, prev], ks)
      if (!(stage > 2L && cc[stage] == cc[stage - 1L] && all(A[stage, prev] == A[stage - 1L, prev]))) {
        du <- f(ui, p, .fma(cc[stage], dt, t))
        out$nf <- out$nf + 1L
      }
      ks[[stage]] <- -Wsolve(.amf_madd(.fma(dtd[stage], dT, du), dtC[stage, prev], ks))
    }
    out$nsolve <- s
    unew <- .amf_madd(u, tab$b, ks)
    EEst <- if (ctx$adaptive) .err_norm(.amf_madd(numeric(length(u)), tab$btilde, ks), u, unew,
                                        ctx$abstol, ctx$reltol) else 0
    kd <- lapply(seq_len(nrow(tab$H)), function(j) .amf_madd(numeric(length(u)), tab$H[j, ], ks))
    c(out, list(u = unew, EEst = EEst, k = kd))
  }
}
