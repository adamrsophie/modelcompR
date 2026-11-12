# Linear operators, after SciMLOperators.jl (scalar.jl, matrix.jl, basic.jl,
# func.jl, tensor.jl, block.jl, woperator.jl, interface.jl).
#
# An operator is an S3 object of class c("<Kind>", "AbstractSciMLOperator")
# with `size` (rows, columns), `apply(v, u, p, t)` (the action on v after
# updating to (u, p, t)), optional `solve_fn(b)` (the inverse action),
# `as_matrix()` (concretization), `update(u, p, t)` (a copy with updated
# coefficients) and traits.

.sop <- function(kind, size, apply, as_matrix = NULL, solve_fn = NULL, update = NULL, linear = TRUE,
                 constant = TRUE, adjoint = NULL, ...) {
  op <- list(kind = kind, size = as.integer(size), apply = apply, as_matrix = as_matrix, solve_fn = solve_fn,
             update = update, linear = linear, constant = constant, adjoint = adjoint, cached = FALSE, ...)
  structure(op, class = c(kind, "AbstractSciMLOperator"))
}

#' @export
print.AbstractSciMLOperator <- function(x, ...) {
  cat(x$kind, " (", x$size[1L], " x ", x$size[2L], ")\n", sep = "")
  invisible(x)
}

#' SciML operators
#'
#' Linear operators with the interface of SciMLOperators.jl, accepted where
#' the solvers take a matrix (linear problems, mass matrices, Jacobian
#' prototypes, the linear part of split and exponential methods).
#' `MatrixOperator(A, update_func)` wraps a matrix whose entries may depend
#' on `(u, p, t)` through `update_func(A, u, p, t)`; `DiagonalOperator(d)`,
#' `ScalarOperator(alpha)`, `IdentityOperator(n)` and `NullOperator(n)` are
#' the special cases; `AffineOperator(A, B, b)` is `v -> A v + B b`
#' (`AddVector(b)` is `v -> v + b`); `FunctionOperator(op, input, output)`
#' is an operator given by its action `op(v, u, p, t)` (with
#' `op_inverse` and `op_adjoint` when known); `InvertibleOperator(L, F)`
#' attaches a factorization; `TensorProductOperator(A, B)` is the Kronecker
#' product `A x B` acting on vectors of length `ncol(A) ncol(B)`,
#' `TensorSumOperator(A, B)` and `kronsum(A, B)` the Kronecker sum `A x I + I
#' x B`; `BlockDiagonalOperator(...)` stacks operators on the diagonal;
#' `WOperator(mass_matrix, gamma, J)` is the matrix `J - M / gamma` of the
#' implicit solvers, with `jacobian_stale()`, `mark_jacobian_updated()` and
#' `mark_jacobian_current()` tracking changes of `J`; `StaticWOperator(W)`
#' a factorized fixed `W`.
#'
#' Operators combine with `+`, `-`, `*` (composition and scaling) and apply
#' with `L %*% v` (R 4.4 and later), `op_apply(L, v, u, p, t)`, or `L(v, u,
#' p, t)` semantics through `op_apply`; `solve(L, b)` applies the inverse,
#' `as.matrix(L)` or `concretize(L)` gives the matrix and `exp(L)` its
#' exponential. `update_coefficients(L, u, p, t)` returns the operator at
#' `(u, p, t)`. The traits `islinear()`, `isconstant()`, `issquare()`,
#' `has_adjoint()`, `has_mul()`, `has_ldiv()`, `has_exp()`, `has_expmv()`,
#' `iscached()`, `isconvertible()` and `has_concretization()` answer
#' questions about an operator; `cache_operator(L, u)` marks it cached.
#'
#' @param A,B,J Matrices or operators.
#' @param d Diagonal entries.
#' @param alpha Scalar value.
#' @param b Vector.
#' @param n Size.
#' @param update_func Coefficient update `f(A, u, p, t)`.
#' @param op,op_inverse,op_adjoint Actions `f(v, u, p, t)`.
#' @param input,output Prototype input and output vectors.
#' @param islinear Whether the function operator is linear.
#' @param F A factorization (any object with a `solve()` method).
#' @param mass_matrix,gamma Mass matrix and step factor of `WOperator()`.
#' @param W,callinv The matrix of `StaticWOperator()` and whether it is applied
#'   through its factorization.
#' @param L An operator.
#' @param v,u,p,t Vector to act on and the state, parameters and time.
#' @param ... Operators or options.
#' @return An operator, a vector, a matrix or a logical.
#' @examples
#' L <- MatrixOperator(matrix(c(2, 0, 1, 3), 2))
#' op_apply(L, c(1, 1))
#' D <- DiagonalOperator(c(1, 2))
#' as.matrix(L + 2 * D)
#' solve(L, c(3, 3))
#' T <- update_coefficients(MatrixOperator(diag(2), function(A, u, p, t) t * diag(2)), NULL, NULL, 3)
#' as.matrix(T)
#' @name sciml_operators
NULL

.as_sop <- function(x) {
  if (inherits(x, "AbstractSciMLOperator")) return(x)
  if (is.matrix(x)) return(MatrixOperator(x))
  if (is.numeric(x) && length(x) == 1L) return(ScalarOperator(x))
  stop("not an operator: ", class(x)[1L], call. = FALSE)
}

#' @rdname sciml_operators
#' @export
op_apply <- function(L, v, u = NULL, p = NULL, t = NULL) {
  L <- .as_sop(L)
  if (!is.null(L$update) && !(is.null(u) && is.null(p) && is.null(t))) L <- L$update(u, p, t)
  as.numeric(L$apply(as.numeric(v)))
}

#' @rdname sciml_operators
#' @export
MatrixOperator <- function(A, update_func = NULL) {
  A <- as.matrix(A)
  op <- .sop("MatrixOperator", dim(A), function(v) as.numeric(A %*% v), function() A,
             function(b) as.numeric(solve(A, b)),
             update = if (!is.null(update_func)) function(u, p, t) MatrixOperator(update_func(A, u, p, t), update_func),
             constant = is.null(update_func), adjoint = function() MatrixOperator(t(A)), A = A)
  op
}

#' @rdname sciml_operators
#' @export
DiagonalOperator <- function(d, update_func = NULL) {
  d <- as.numeric(d)
  .sop("DiagonalOperator", c(length(d), length(d)), function(v) d * v, function() diag(d, length(d)),
       function(b) b / d,
       update = if (!is.null(update_func)) function(u, p, t) DiagonalOperator(update_func(d, u, p, t), update_func),
       constant = is.null(update_func), adjoint = function() DiagonalOperator(d), diag = d)
}

#' @rdname sciml_operators
#' @export
ScalarOperator <- function(alpha, update_func = NULL) {
  .sop("ScalarOperator", c(NA, NA), function(v) alpha * v, function() alpha, function(b) b / alpha,
       update = if (!is.null(update_func)) function(u, p, t) ScalarOperator(update_func(alpha, u, p, t), update_func),
       constant = is.null(update_func), adjoint = function() ScalarOperator(alpha), value = alpha)
}

#' @rdname sciml_operators
#' @export
IdentityOperator <- function(n)
  .sop("IdentityOperator", c(n, n), function(v) v, function() diag(n), function(b) b,
       adjoint = function() IdentityOperator(n))

#' @rdname sciml_operators
#' @export
NullOperator <- function(n)
  .sop("NullOperator", c(n, n), function(v) 0 * v, function() matrix(0, n, n), NULL,
       adjoint = function() NullOperator(n))

#' @rdname sciml_operators
#' @export
AffineOperator <- function(A, B, b, update_func = NULL) {
  A <- .as_sop(A); B <- .as_sop(B); b <- as.numeric(b)
  .sop("AffineOperator", A$size, function(v) A$apply(v) + B$apply(b), NULL,
       function(y) as.numeric(.sop_solve(A, y - B$apply(b))), linear = FALSE,
       update = function(u, p, t) {
         A2 <- if (!is.null(A$update)) A$update(u, p, t) else A
         B2 <- if (!is.null(B$update)) B$update(u, p, t) else B
         b2 <- if (!is.null(update_func)) update_func(b, u, p, t) else b
         AffineOperator(A2, B2, b2, update_func)
       },
       constant = isTRUE(A$constant) && isTRUE(B$constant) && is.null(update_func), A = A, B = B, b = b)
}

#' @rdname sciml_operators
#' @export
AddVector <- function(b, update_func = NULL) {
  n <- length(b)
  AffineOperator(IdentityOperator(n), IdentityOperator(n), b, update_func)
}

#' @rdname sciml_operators
#' @export
FunctionOperator <- function(op, input, output = input, op_inverse = NULL, op_adjoint = NULL, islinear = TRUE) {
  holder <- new.env(parent = emptyenv())
  holder$u <- NULL; holder$p <- NULL; holder$t <- NULL
  make <- function(u, p, t) {
    .sop("FunctionOperator", c(length(output), length(input)), function(v) as.numeric(op(v, u, p, t)),
         function() {
           n <- length(input)
           do.call(cbind, lapply(seq_len(n), function(j) { e <- numeric(n); e[j] <- 1; as.numeric(op(e, u, p, t)) }))
         },
         if (!is.null(op_inverse)) function(b) as.numeric(op_inverse(b, u, p, t)),
         update = function(u2, p2, t2) make(u2, p2, t2), linear = islinear, constant = FALSE,
         adjoint = if (!is.null(op_adjoint)) function() FunctionOperator(op_adjoint, output, input, NULL, op, islinear))
  }
  make(NULL, NULL, NULL)
}

#' @rdname sciml_operators
#' @export
InvertibleOperator <- function(L, F = NULL) {
  L <- .as_sop(L)
  F <- F %||% (if (!is.null(L$as_matrix)) qr(L$as_matrix()))
  op <- L
  op$kind <- "InvertibleOperator"
  op$solve_fn <- function(b) as.numeric(if (inherits(F, "qr")) qr.coef(F, b) else solve(F, b))
  op$F <- F
  class(op) <- c("InvertibleOperator", "AbstractSciMLOperator")
  op
}

#' @rdname sciml_operators
#' @export
TensorProductOperator <- function(A, B) {
  A <- .as_sop(A); B <- .as_sop(B)
  ma <- A$size; mb <- B$size
  # (A x B) vec(V) = vec(B V A') for V with ncol(B) rows and ncol(A) columns
  app <- function(v) {
    V <- matrix(v, nrow = mb[2L], ncol = ma[2L])
    BV <- apply(V, 2L, B$apply)
    if (is.null(dim(BV))) BV <- matrix(BV, nrow = mb[1L])
    as.numeric(t(apply(t(BV), 2L, A$apply)))
  }
  .sop("TensorProductOperator", ma * mb, app,
       function() kronecker(A$as_matrix(), B$as_matrix()),
       if (!is.null(A$solve_fn) && !is.null(B$solve_fn)) function(b) {
         V <- matrix(b, nrow = mb[1L], ncol = ma[1L])
         BV <- apply(V, 2L, B$solve_fn)
         if (is.null(dim(BV))) BV <- matrix(BV, nrow = mb[2L])
         as.numeric(t(apply(t(BV), 2L, A$solve_fn)))
       },
       update = function(u, p, t) TensorProductOperator(if (!is.null(A$update)) A$update(u, p, t) else A,
                                                        if (!is.null(B$update)) B$update(u, p, t) else B),
       constant = isTRUE(A$constant) && isTRUE(B$constant), A = A, B = B)
}

#' @rdname sciml_operators
#' @export
TensorSumOperator <- function(A, B) {
  A <- .as_sop(A); B <- .as_sop(B)
  TensorProductOperator(A, IdentityOperator(B$size[1L])) + TensorProductOperator(IdentityOperator(A$size[1L]), B)
}

#' @rdname sciml_operators
#' @export
kronsum <- function(A, B) TensorSumOperator(A, B)

#' @rdname sciml_operators
#' @export
BlockDiagonalOperator <- function(...) {
  ops <- lapply(list(...), .as_sop)
  rows <- vapply(ops, function(o) o$size[1L], 0L); cols <- vapply(ops, function(o) o$size[2L], 0L)
  ri <- split(seq_len(sum(rows)), rep(seq_along(ops), rows))
  ci <- split(seq_len(sum(cols)), rep(seq_along(ops), cols))
  .sop("BlockDiagonalOperator", c(sum(rows), sum(cols)),
       function(v) unlist(lapply(seq_along(ops), function(i) ops[[i]]$apply(v[ci[[i]]]))),
       function() {
         M <- matrix(0, sum(rows), sum(cols))
         for (i in seq_along(ops)) M[ri[[i]], ci[[i]]] <- ops[[i]]$as_matrix()
         M
       },
       function(b) unlist(lapply(seq_along(ops), function(i) .sop_solve(ops[[i]], b[ri[[i]]]))),
       constant = all(vapply(ops, function(o) isTRUE(o$constant), TRUE)), ops = ops)
}

# W = J - M / gamma
#' @rdname sciml_operators
#' @export
WOperator <- function(mass_matrix, gamma, J, u = NULL) {
  state <- new.env(parent = emptyenv())
  state$stale <- TRUE
  Jop <- .as_sop(J)
  n <- Jop$size[1L]
  M <- if (is.null(mass_matrix)) IdentityOperator(n) else .as_sop(mass_matrix)
  op <- .sop("WOperator", Jop$size, function(v) Jop$apply(v) - M$apply(v) / gamma,
             function() .sop_matrix(Jop, n) - .sop_matrix(M, n) / gamma,
             function(b) as.numeric(solve(.sop_matrix(Jop, n) - .sop_matrix(M, n) / gamma, b)),
             constant = FALSE, J = Jop, mass_matrix = M, gamma = gamma, state = state)
  op
}

#' @rdname sciml_operators
#' @export
StaticWOperator <- function(W, callinv = TRUE) {
  W <- as.matrix(W)
  F <- qr(W)
  .sop("StaticWOperator", dim(W), function(v) as.numeric(W %*% v), function() W, function(b) qr.coef(F, b), W = W)
}

#' @rdname sciml_operators
#' @export
jacobian_stale <- function(L) isTRUE(L$state$stale)

#' @rdname sciml_operators
#' @export
mark_jacobian_updated <- function(L) {
  L$state$stale <- TRUE
  invisible(L)
}

#' @rdname sciml_operators
#' @export
mark_jacobian_current <- function(L) {
  L$state$stale <- FALSE
  invisible(L)
}

.sop_matrix <- function(L, n = NULL) {
  if (inherits(L, "ScalarOperator")) return(L$value * diag(n %||% 1))
  if (is.null(L$as_matrix)) stop(L$kind, " has no matrix form", call. = FALSE)
  L$as_matrix()
}

.sop_solve <- function(L, b) {
  if (!is.null(L$solve_fn)) return(L$solve_fn(b))
  as.numeric(solve(.sop_matrix(L, length(b)), b))
}

# ---------------------------------------------------------------------------
# Algebra

.sop_compose <- function(A, B) {
  .sop("ComposedOperator", c(if (is.na(A$size[1L])) B$size[1L] else A$size[1L], if (is.na(B$size[2L])) A$size[2L] else B$size[2L]),
       function(v) A$apply(B$apply(v)),
       if (!is.null(A$as_matrix) && !is.null(B$as_matrix)) function() {
         a <- A$as_matrix(); b <- B$as_matrix()
         if (length(a) == 1L || length(b) == 1L) a * b else a %*% b
       },
       if (!is.null(A$solve_fn) && !is.null(B$solve_fn)) function(y) B$solve_fn(A$solve_fn(y)),
       update = function(u, p, t) .sop_compose(if (!is.null(A$update)) A$update(u, p, t) else A,
                                               if (!is.null(B$update)) B$update(u, p, t) else B),
       linear = isTRUE(A$linear) && isTRUE(B$linear), constant = isTRUE(A$constant) && isTRUE(B$constant),
       ops = list(A, B))
}

.sop_add <- function(A, B, sgn = 1) {
  sz <- if (anyNA(A$size)) B$size else A$size
  .sop("AddedOperator", sz, function(v) A$apply(v) + sgn * B$apply(v),
       if (!is.null(A$as_matrix) && !is.null(B$as_matrix)) function() .sop_matrix(A, sz[1L]) + sgn * .sop_matrix(B, sz[1L]),
       NULL,
       update = function(u, p, t) .sop_add(if (!is.null(A$update)) A$update(u, p, t) else A,
                                           if (!is.null(B$update)) B$update(u, p, t) else B, sgn),
       linear = isTRUE(A$linear) && isTRUE(B$linear), constant = isTRUE(A$constant) && isTRUE(B$constant),
       ops = list(A, B))
}

#' @export
Ops.AbstractSciMLOperator <- function(e1, e2) {
  if (.Generic == "-" && missing(e2)) return(.sop_compose(ScalarOperator(-1), .as_sop(e1)))
  a <- .as_sop(e1); b <- .as_sop(e2)
  switch(.Generic,
    "+" = .sop_add(a, b), "-" = .sop_add(a, b, -1), "*" = .sop_compose(a, b),
    "/" = if (inherits(b, "ScalarOperator")) .sop_compose(ScalarOperator(1 / b$value), a)
          else stop("divide an operator by a number only", call. = FALSE),
    stop("operator `", .Generic, "` is not defined for SciML operators", call. = FALSE))
}

#' @export
as.matrix.AbstractSciMLOperator <- function(x, ...) .sop_matrix(x)

#' @export
solve.AbstractSciMLOperator <- function(a, b, ...) .sop_solve(a, as.numeric(b))

#' @export
t.AbstractSciMLOperator <- function(x) {
  if (!is.null(x$adjoint)) return(x$adjoint())
  MatrixOperator(t(.sop_matrix(x)))
}

#' @export
dim.AbstractSciMLOperator <- function(x) x$size

#' @export
exp.AbstractSciMLOperator <- function(x) .expm(.sop_matrix(x))

.sop_matrixops <- function(x, y) {
  if (inherits(x, "AbstractSciMLOperator") && !inherits(y, "AbstractSciMLOperator")) {
    if (is.matrix(y)) return(apply(y, 2L, x$apply))
    return(x$apply(as.numeric(y)))
  }
  as.matrix(.sop_compose(.as_sop(x), .as_sop(y)))
}

#' @rawNamespace if (getRversion() >= "4.4.0") S3method(matrixOps, AbstractSciMLOperator, .sop_matrixops)
NULL

# Matrix exponential by scaling and squaring with a [13/13] Pade approximant
# (Higham 2005).
.expm <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  if (n == 0L) return(A)
  b <- c(64764752532480000, 32382376266240000, 7771770303897600, 1187353796428800, 129060195264000,
         10559470521600, 670442572800, 33522128640, 1323241920, 40840800, 960960, 16380, 182, 1)
  nrm <- max(colSums(abs(A)))
  s <- max(0L, ceiling(log2(nrm / 5.371920351148152)))
  A <- A / 2^s
  I <- diag(n)
  A2 <- A %*% A; A4 <- A2 %*% A2; A6 <- A4 %*% A2
  U <- A %*% (A6 %*% (b[14] * A6 + b[12] * A4 + b[10] * A2) + b[8] * A6 + b[6] * A4 + b[4] * A2 + b[2] * I)
  V <- A6 %*% (b[13] * A6 + b[11] * A4 + b[9] * A2) + b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * I
  E <- solve(V - U, V + U)
  for (i in seq_len(s)) E <- E %*% E
  E
}

# ---------------------------------------------------------------------------
# Interface and traits

#' @rdname sciml_operators
#' @export
update_coefficients <- function(L, u = NULL, p = NULL, t = NULL) {
  L <- .as_sop(L)
  if (is.null(L$update)) L else L$update(u, p, t)
}

#' @rdname sciml_operators
#' @export
cache_operator <- function(L, u = NULL) {
  L <- .as_sop(L)
  L$cached <- TRUE
  L
}

#' @rdname sciml_operators
#' @export
concretize <- function(L) .sop_matrix(.as_sop(L))

#' @rdname sciml_operators
#' @export
islinear <- function(L) isTRUE(.as_sop(L)$linear)

#' @rdname sciml_operators
#' @export
issquare <- function(L) {
  if (is.matrix(L)) return(nrow(L) == ncol(L))
  s <- .as_sop(L)$size
  anyNA(s) || s[1L] == s[2L]
}

#' @rdname sciml_operators
#' @export
has_adjoint <- function(L) !is.null(.as_sop(L)$adjoint) || !is.null(.as_sop(L)$as_matrix)

#' @rdname sciml_operators
#' @export
has_mul <- function(L) TRUE

#' @rdname sciml_operators
#' @export
has_ldiv <- function(L) {
  L <- .as_sop(L)
  !is.null(L$solve_fn)
}

#' @rdname sciml_operators
#' @export
has_exp <- function(L) !is.null(.as_sop(L)$as_matrix)

#' @rdname sciml_operators
#' @export
has_expmv <- function(L) has_exp(L)

#' @rdname sciml_operators
#' @export
iscached <- function(L) isTRUE(.as_sop(L)$cached)

#' @rdname sciml_operators
#' @export
isconvertible <- function(L) !is.null(.as_sop(L)$as_matrix)

#' @rdname sciml_operators
#' @export
has_concretization <- function(L) isconvertible(L)
