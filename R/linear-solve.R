# Linear solvers, after LinearSolve.jl (reexported by NonlinearSolve): the
# factorization algorithms (src/factorization.jl, the extension backends),
# the Krylov methods (src/iterative_wrappers.jl and the KrylovKit,
# IterativeSolvers, Ginkgo, ConjugateGradients wrappers), simple_lu.jl
# (SimpleLUFactorization), simplegmres.jl (SimpleGMRES), the default
# algorithm choice (src/default.jl with OperatorAssumptions and
# OperatorCondition), preconditioners (InvPreconditioner), LinearCache
# (init, solve!), LinearSolveFunction, DirectLdiv!, LowRankUpdatedMatrix
# and the adjoint LinearSolveAdjoint; and the eigenvalue interface of
# SciMLBase (EigenvalueProblem, EigenvalueSolution, EigenvalueTarget) with
# the dense, Arnoldi, Arpack, KrylovKit and Jacobi-Davidson solvers.
#
# R links one BLAS and LAPACK; the backends of LinearSolve that select a
# library or a device (MKL, OpenBLAS, AppleAccelerate, BLIS, CUDA, Metal,
# AMDGPU, sparse direct solvers, PETSc, HYPRE, ...) run the same
# factorization through it. The Krylov methods are implemented here.

.lin_alg <- function(name, method, ...) structure(list(name = name, method = method, ...),
                                                  class = c(name, "LinearAlgorithm"))

#' @export
print.LinearAlgorithm <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

#' Linear solvers
#'
#' Algorithms for `solve()` on a [LinearProblem()] (`A u = b`), with the
#' names of LinearSolve.jl.
#'
#' Direct methods: `LUFactorization()`, `QRFactorization()` (least squares
#' for rectangular `A`), `SVDFactorization()`, `CholeskyFactorization()`
#' (symmetric positive definite `A`), `LDLtFactorization()` and
#' `BunchKaufmanFactorization()` (symmetric `A`), `NormalCholeskyFactorization()`
#' and `NormalBunchKaufmanFactorization()` (the normal equations),
#' `DiagonalFactorization()`, `GenericFactorization(fact)` (any function
#' returning an object for [solve()]), `SimpleLUFactorization()` (LU without
#' pivoting in R code) and `DirectLdiv()` (Julia's `DirectLdiv!`, `solve(A,
#' b)`). The library and device backends `GenericLUFactorization`,
#' `RFLUFactorization`, `FastLUFactorization`, `FastQRFactorization`,
#' `SpecializedLUFactorization`, `SpecializedQRFactorization`,
#' `MKLLUFactorization`, `MKL32MixedLUFactorization`, `OpenBLASLUFactorization`,
#' `OpenBLAS32MixedLUFactorization`, `AppleAccelerateLUFactorization`,
#' `AppleAccelerate32MixedLUFactorization`, `BLISLUFactorization`,
#' `RF32MixedLUFactorization`, `MetalLUFactorization`,
#' `MetalOffload32MixedLUFactorization`, `CudaOffloadFactorization`,
#' `CudaOffloadLUFactorization`, `CudaOffloadQRFactorization`,
#' `CUDAOffload32MixedLUFactorization`, `CUSOLVERRFFactorization`,
#' `AMDGPUOffloadLUFactorization`, `AMDGPUOffloadQRFactorization`,
#' `GESVFactorization`, `ButterflyFactorization`, `LHLFactorization` and the
#' sparse direct solvers `KLUFactorization`, `PureKLUFactorization`,
#' `UMFPACKFactorization`, `PureUMFPACKFactorization`,
#' `SparspakFactorization`, `CHOLMODFactorization`,
#' `SparseColumnPivotedQRFactorization`, `CliqueTreesFactorization`,
#' `ParUFactorization`, `CKTSOFactorization`, `HSLMA57Factorization`,
#' `HSLMA97Factorization`, `MUMPSFactorization`, `MKLPardisoFactorize`,
#' `MKLPardisoIterate`, `PanuaPardisoFactorize`, `PanuaPardisoIterate`,
#' `PardisoJL`, `SuperLUDISTFactorization`, `SupernodalLUFactorization`,
#' `STRUMPACKFactorization`, `SLATEFactorization` and `ElementalJL` select a
#' library or a device in Julia; R has one BLAS and LAPACK, so they compute
#' the same LU, QR or Cholesky factorization through it (dense).
#'
#' Krylov methods (implemented in R; `abstol`, `reltol`, `maxiters`,
#' `restart` and the preconditioners `Pl`, `Pr` are options of `solve()`):
#' `KrylovJL_GMRES()`, `KrylovJL_FGMRES()` (restarted GMRES with right
#' preconditioning), `KrylovJL_CG()` (conjugate gradients), `KrylovJL_MINRES()`
#' and `KrylovJL_MINARES()` (symmetric indefinite), `KrylovJL_BICGSTAB()`,
#' `KrylovJL_LSMR()` (least squares, by CGLS) and `KrylovJL_CRAIGMR()`
#' (minimum norm, by CGNE); `KrylovJL(method)` chooses among them, and the
#' wrappers `IterativeSolversJL(_CG, _GMRES, _BICGSTAB, _MINRES, _IDRS)`,
#' `KrylovKitJL(_CG, _GMRES)`, `GinkgoJL(_CG, _GMRES)`,
#' `ConjugateGradientsJL(_CG, _BICGSTAB)` and `SimpleGMRES()` name the same
#' methods (`_IDRS` runs BiCGSTAB, which is IDR(1)). `AlgebraicMultigridJL()`,
#' `AMGXPreconditioner()`, `HYPREAlgorithm()` and `PETScAlgorithm()` run
#' GMRES (R has no algebraic multigrid). `InvPreconditioner(P)` applies the
#' inverse of `P`.
#'
#' `LinearSolveFunction(f)` solves with `f(A, b, u, p, ...)`. Without an
#' algorithm, `OperatorAssumptions(issquare, condition)` steer the default as
#' LinearSolve's `defaultalg`: square systems use LU unless the condition
#' is `VeryIllConditioned` (QR) or `SuperIllConditioned` (SVD); rectangular
#' ones use QR, or the normal equations by Cholesky when `WellConditioned`
#' (`OperatorCondition` lists the levels). `init(prob, alg)` returns a
#' `LinearCache` that keeps the factorization: change `cache$A` or `cache$b`
#' and call `solve_integrator(cache)` again. `LowRankUpdatedMatrix(A, U, C,
#' V)` represents `A + U C V'` and is solved with the Woodbury identity.
#' `LinearSolveAdjoint()` is the sensitivity method of linear solves: the
#' gradient of a loss `g(u)` is `-lambda u'` for `A` and `lambda` for `b`,
#' `lambda = A' \ g_u` (see `linear_adjoint()`). `LinearVerbosity()` holds
#' verbosity settings.
#'
#' @param method,fact,f,P,A,U,C,V Arguments as described above.
#' @param issquare,condition,nonstructural_zeros Operator assumptions
#'   (`nonstructural_zeros` is a [NonstructuralZeros] value).
#' @param level Verbosity level.
#' @param u,dgdu A solution and the gradient of a loss with respect to it.
#' @param ... Options.
#' @return An algorithm, a cache, an operator or a gradient.
#' @examples
#' A <- matrix(c(4, 1, 1, 3), 2); b <- c(1, 2)
#' prob <- LinearProblem(A, b)
#' solve(prob, KrylovJL_CG())$u
#' solve(prob, CholeskyFactorization())$u
#' cache <- init(prob, LUFactorization())
#' solve_integrator(cache)$u
#' cache$b <- c(2, 4)
#' solve_integrator(cache)$u
#' @name linear_solvers
NULL

#' @rdname linear_solvers
#' @export
LUFactorization <- function() .lin_alg("LUFactorization", "lu")
#' @rdname linear_solvers
#' @export
QRFactorization <- function() .lin_alg("QRFactorization", "qr")
#' @rdname linear_solvers
#' @export
SVDFactorization <- function() .lin_alg("SVDFactorization", "svd")
#' @rdname linear_solvers
#' @export
CholeskyFactorization <- function() .lin_alg("CholeskyFactorization", "chol")
#' @rdname linear_solvers
#' @export
LDLtFactorization <- function() .lin_alg("LDLtFactorization", "ldlt")
#' @rdname linear_solvers
#' @export
BunchKaufmanFactorization <- function() .lin_alg("BunchKaufmanFactorization", "ldlt")
#' @rdname linear_solvers
#' @export
NormalCholeskyFactorization <- function() .lin_alg("NormalCholeskyFactorization", "normal_chol")
#' @rdname linear_solvers
#' @export
NormalBunchKaufmanFactorization <- function() .lin_alg("NormalBunchKaufmanFactorization", "normal_chol")
#' @rdname linear_solvers
#' @export
DiagonalFactorization <- function() .lin_alg("DiagonalFactorization", "diag")
#' @rdname linear_solvers
#' @export
GenericFactorization <- function(fact = qr) .lin_alg("GenericFactorization", "generic", fact = fact)
#' @rdname linear_solvers
#' @export
SimpleLUFactorization <- function() .lin_alg("SimpleLUFactorization", "simple_lu")
#' @rdname linear_solvers
#' @export
DirectLdiv <- function() .lin_alg("DirectLdiv", "lu")


.lin_backend <- function(name, method) function(...) .lin_alg(name, method, ...)

#' @rdname linear_solvers
#' @export
GenericLUFactorization <- .lin_backend("GenericLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
RFLUFactorization <- .lin_backend("RFLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
FastLUFactorization <- .lin_backend("FastLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
FastQRFactorization <- .lin_backend("FastQRFactorization", "qr")
#' @rdname linear_solvers
#' @export
SpecializedLUFactorization <- .lin_backend("SpecializedLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
SpecializedQRFactorization <- .lin_backend("SpecializedQRFactorization", "qr")
#' @rdname linear_solvers
#' @export
MKLLUFactorization <- .lin_backend("MKLLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
MKL32MixedLUFactorization <- .lin_backend("MKL32MixedLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
OpenBLASLUFactorization <- .lin_backend("OpenBLASLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
OpenBLAS32MixedLUFactorization <- .lin_backend("OpenBLAS32MixedLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
AppleAccelerateLUFactorization <- .lin_backend("AppleAccelerateLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
AppleAccelerate32MixedLUFactorization <- .lin_backend("AppleAccelerate32MixedLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
BLISLUFactorization <- .lin_backend("BLISLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
RF32MixedLUFactorization <- .lin_backend("RF32MixedLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
MetalLUFactorization <- .lin_backend("MetalLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
MetalOffload32MixedLUFactorization <- .lin_backend("MetalOffload32MixedLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
CudaOffloadFactorization <- .lin_backend("CudaOffloadFactorization", "lu")
#' @rdname linear_solvers
#' @export
CudaOffloadLUFactorization <- .lin_backend("CudaOffloadLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
CudaOffloadQRFactorization <- .lin_backend("CudaOffloadQRFactorization", "qr")
#' @rdname linear_solvers
#' @export
CUDAOffload32MixedLUFactorization <- .lin_backend("CUDAOffload32MixedLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
CUSOLVERRFFactorization <- .lin_backend("CUSOLVERRFFactorization", "lu")
#' @rdname linear_solvers
#' @export
AMDGPUOffloadLUFactorization <- .lin_backend("AMDGPUOffloadLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
AMDGPUOffloadQRFactorization <- .lin_backend("AMDGPUOffloadQRFactorization", "qr")
#' @rdname linear_solvers
#' @export
GESVFactorization <- .lin_backend("GESVFactorization", "lu")
#' @rdname linear_solvers
#' @export
ButterflyFactorization <- .lin_backend("ButterflyFactorization", "lu")
#' @rdname linear_solvers
#' @export
LHLFactorization <- .lin_backend("LHLFactorization", "lu")
#' @rdname linear_solvers
#' @export
KLUFactorization <- .lin_backend("KLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
PureKLUFactorization <- .lin_backend("PureKLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
UMFPACKFactorization <- .lin_backend("UMFPACKFactorization", "lu")
#' @rdname linear_solvers
#' @export
PureUMFPACKFactorization <- .lin_backend("PureUMFPACKFactorization", "lu")
#' @rdname linear_solvers
#' @export
SparspakFactorization <- .lin_backend("SparspakFactorization", "lu")
#' @rdname linear_solvers
#' @export
CHOLMODFactorization <- .lin_backend("CHOLMODFactorization", "chol")
#' @rdname linear_solvers
#' @export
SparseColumnPivotedQRFactorization <- .lin_backend("SparseColumnPivotedQRFactorization", "qr")
#' @rdname linear_solvers
#' @export
CliqueTreesFactorization <- .lin_backend("CliqueTreesFactorization", "chol")
#' @rdname linear_solvers
#' @export
ParUFactorization <- .lin_backend("ParUFactorization", "lu")
#' @rdname linear_solvers
#' @export
CKTSOFactorization <- .lin_backend("CKTSOFactorization", "lu")
#' @rdname linear_solvers
#' @export
HSLMA57Factorization <- .lin_backend("HSLMA57Factorization", "ldlt")
#' @rdname linear_solvers
#' @export
HSLMA97Factorization <- .lin_backend("HSLMA97Factorization", "ldlt")
#' @rdname linear_solvers
#' @export
MUMPSFactorization <- .lin_backend("MUMPSFactorization", "lu")
#' @rdname linear_solvers
#' @export
MKLPardisoFactorize <- .lin_backend("MKLPardisoFactorize", "lu")
#' @rdname linear_solvers
#' @export
MKLPardisoIterate <- .lin_backend("MKLPardisoIterate", "lu")
#' @rdname linear_solvers
#' @export
PanuaPardisoFactorize <- .lin_backend("PanuaPardisoFactorize", "lu")
#' @rdname linear_solvers
#' @export
PanuaPardisoIterate <- .lin_backend("PanuaPardisoIterate", "lu")
#' @rdname linear_solvers
#' @export
PardisoJL <- .lin_backend("PardisoJL", "lu")
#' @rdname linear_solvers
#' @export
SuperLUDISTFactorization <- .lin_backend("SuperLUDISTFactorization", "lu")
#' @rdname linear_solvers
#' @export
SupernodalLUFactorization <- .lin_backend("SupernodalLUFactorization", "lu")
#' @rdname linear_solvers
#' @export
STRUMPACKFactorization <- .lin_backend("STRUMPACKFactorization", "lu")
#' @rdname linear_solvers
#' @export
SLATEFactorization <- .lin_backend("SLATEFactorization", "lu")
#' @rdname linear_solvers
#' @export
ElementalJL <- .lin_backend("ElementalJL", "lu")

# Krylov methods
.kry <- function(name, method) function(...) .lin_alg(name, "krylov", krylov = method, ...)

#' @rdname linear_solvers
#' @export
KrylovJL <- function(method = c("gmres", "cg", "minres", "bicgstab", "lsmr", "craigmr", "fgmres", "minares"), ...)
  .lin_alg("KrylovJL", "krylov", krylov = match.arg(method), ...)
#' @rdname linear_solvers
#' @export
KrylovJL_GMRES <- .kry("KrylovJL_GMRES", "gmres")
#' @rdname linear_solvers
#' @export
KrylovJL_FGMRES <- .kry("KrylovJL_FGMRES", "fgmres")
#' @rdname linear_solvers
#' @export
KrylovJL_CG <- .kry("KrylovJL_CG", "cg")
#' @rdname linear_solvers
#' @export
KrylovJL_MINRES <- .kry("KrylovJL_MINRES", "minres")
#' @rdname linear_solvers
#' @export
KrylovJL_MINARES <- .kry("KrylovJL_MINARES", "minres")
#' @rdname linear_solvers
#' @export
KrylovJL_BICGSTAB <- .kry("KrylovJL_BICGSTAB", "bicgstab")
#' @rdname linear_solvers
#' @export
KrylovJL_LSMR <- .kry("KrylovJL_LSMR", "lsmr")
#' @rdname linear_solvers
#' @export
KrylovJL_CRAIGMR <- .kry("KrylovJL_CRAIGMR", "craigmr")
#' @rdname linear_solvers
#' @export
IterativeSolversJL <- function(method = "gmres", ...) .lin_alg("IterativeSolversJL", "krylov", krylov = method, ...)
#' @rdname linear_solvers
#' @export
IterativeSolversJL_CG <- .kry("IterativeSolversJL_CG", "cg")
#' @rdname linear_solvers
#' @export
IterativeSolversJL_GMRES <- .kry("IterativeSolversJL_GMRES", "gmres")
#' @rdname linear_solvers
#' @export
IterativeSolversJL_BICGSTAB <- .kry("IterativeSolversJL_BICGSTAB", "bicgstab")
#' @rdname linear_solvers
#' @export
IterativeSolversJL_MINRES <- .kry("IterativeSolversJL_MINRES", "minres")
#' @rdname linear_solvers
#' @export
IterativeSolversJL_IDRS <- .kry("IterativeSolversJL_IDRS", "bicgstab")
#' @rdname linear_solvers
#' @export
KrylovKitJL <- function(method = "gmres", ...) .lin_alg("KrylovKitJL", "krylov", krylov = method, ...)
#' @rdname linear_solvers
#' @export
KrylovKitJL_CG <- .kry("KrylovKitJL_CG", "cg")
#' @rdname linear_solvers
#' @export
KrylovKitJL_GMRES <- .kry("KrylovKitJL_GMRES", "gmres")
#' @rdname linear_solvers
#' @export
GinkgoJL <- function(method = "gmres", ...) .lin_alg("GinkgoJL", "krylov", krylov = method, ...)
#' @rdname linear_solvers
#' @export
GinkgoJL_CG <- .kry("GinkgoJL_CG", "cg")
#' @rdname linear_solvers
#' @export
GinkgoJL_GMRES <- .kry("GinkgoJL_GMRES", "gmres")
#' @rdname linear_solvers
#' @export
ConjugateGradientsJL <- function(method = "cg", ...) .lin_alg("ConjugateGradientsJL", "krylov", krylov = method, ...)
#' @rdname linear_solvers
#' @export
ConjugateGradientsJL_CG <- .kry("ConjugateGradientsJL_CG", "cg")
#' @rdname linear_solvers
#' @export
ConjugateGradientsJL_BICGSTAB <- .kry("ConjugateGradientsJL_BICGSTAB", "bicgstab")
#' @rdname linear_solvers
#' @export
SimpleGMRES <- .kry("SimpleGMRES", "gmres")
#' @rdname linear_solvers
#' @export
AlgebraicMultigridJL <- .kry("AlgebraicMultigridJL", "gmres")
#' @rdname linear_solvers
#' @export
AMGXPreconditioner <- .kry("AMGXPreconditioner", "gmres")
#' @rdname linear_solvers
#' @export
HYPREAlgorithm <- .kry("HYPREAlgorithm", "gmres")
#' @rdname linear_solvers
#' @export
PETScAlgorithm <- .kry("PETScAlgorithm", "gmres")

#' @rdname linear_solvers
#' @export
InvPreconditioner <- function(P) structure(list(P = P), class = "InvPreconditioner")

#' @rdname linear_solvers
#' @export
LinearSolveFunction <- function(f) .lin_alg("LinearSolveFunction", "function", fun = f)

#' @rdname linear_solvers
#' @export
OperatorCondition <- list(IllConditioned = "IllConditioned", VeryIllConditioned = "VeryIllConditioned",
                          SuperIllConditioned = "SuperIllConditioned", WellConditioned = "WellConditioned")

#' @rdname linear_solvers
#' @export
OperatorAssumptions <- function(issquare = NULL, condition = OperatorCondition$IllConditioned,
                                nonstructural_zeros = NonstructuralZeros$Auto)
  structure(list(issquare = issquare, condition = condition, nonstructural_zeros = nonstructural_zeros),
            class = "OperatorAssumptions")

#' @rdname linear_solvers
#' @export
LinearVerbosity <- function(level = c("warn", "none", "info", "all"))
  structure(list(level = match.arg(level)), class = "LinearVerbosity")

#' @rdname linear_solvers
#' @export
LinearSolveAdjoint <- function() structure(list(name = "LinearSolveAdjoint"),
                                           class = c("LinearSolveAdjoint", "AdjointSensitivityAlgorithm", "SensitivityAlgorithm"))

#' @rdname linear_solvers
#' @export
linear_adjoint <- function(A, u, dgdu) {
  lambda <- as.numeric(solve(t(A), as.numeric(dgdu)))
  list(dA = -outer(lambda, as.numeric(u)), db = lambda)
}

#' @rdname linear_solvers
#' @export
LowRankUpdatedMatrix <- function(A, U, C, V) structure(list(A = A, U = as.matrix(U), C = as.matrix(C), V = as.matrix(V)),
                                                       class = "LowRankUpdatedMatrix")

# ---------------------------------------------------------------------------
# Solving

.precond <- function(P) {
  if (is.null(P)) return(function(x) x)
  if (inherits(P, "InvPreconditioner")) return(function(x) as.numeric(P$P %*% x))
  if (is.function(P)) return(P)
  function(x) as.numeric(solve(P, x))
}

# the default algorithm (defaultalg of LinearSolve)
.lin_default <- function(A, assump) {
  sq <- assump$issquare %||% (is.function(A) || nrow(A) == ncol(A))
  cond <- assump$condition %||% OperatorCondition$IllConditioned
  if (cond == "SuperIllConditioned") return(SVDFactorization())
  if (cond == "VeryIllConditioned") return(QRFactorization())
  if (sq) return(LUFactorization())
  if (cond == "WellConditioned") NormalCholeskyFactorization() else QRFactorization()
}

.lin_direct <- function(A, b, alg, cache) {
  m <- alg$method
  fresh <- is.null(cache$fact) || !identical(cache$fact_A, A)
  if (fresh && m %in% c("lu", "qr", "chol", "ldlt", "svd", "normal_chol", "generic")) {
    cache$fact <- switch(m,
      lu = Matrix::lu(A), qr = qr(A),
      chol = chol(A), ldlt = Matrix::BunchKaufman(A), svd = svd(A),
      normal_chol = chol(crossprod(A)), generic = alg$fact(A))
    cache$fact_A <- A
  }
  f <- cache$fact
  switch(m,
    lu = as.numeric(Matrix::solve(f, b)), qr = qr.coef(f, b), ldlt = as.numeric(Matrix::solve(f, b)),
    chol = backsolve(f, forwardsolve(t(f), b)),
    svd = {
      d <- f$d
      tol <- max(dim(A)) * max(d) * .Machine$double.eps
      as.numeric(f$v %*% (ifelse(d > tol, 1 / d, 0) * crossprod(f$u, b)))
    },
    normal_chol = backsolve(f, forwardsolve(t(f), crossprod(A, b))),
    generic = as.numeric(solve(f, b)),
    diag = b / diag(A),
    simple_lu = .simple_lu_solve(A, b))
}

# LU without pivoting (SimpleLUFactorization, simple_lu.jl)
.simple_lu_solve <- function(A, b) {
  n <- nrow(A)
  for (k in seq_len(n - 1L)) {
    for (i in (k + 1L):n) {
      A[i, k] <- A[i, k] / A[k, k]
      A[i, (k + 1L):n] <- A[i, (k + 1L):n] - A[i, k] * A[k, (k + 1L):n]
    }
  }
  y <- b
  for (i in seq_len(n)) if (i > 1L) y[i] <- y[i] - sum(A[i, seq_len(i - 1L)] * y[seq_len(i - 1L)])
  x <- y
  for (i in n:1L) {
    s <- if (i < n) sum(A[i, (i + 1L):n] * x[(i + 1L):n]) else 0
    x[i] <- (y[i] - s) / A[i, i]
  }
  x
}

# Krylov methods: x such that A x = b to |r| <= max(abstol, reltol |b|)
.krylov <- function(A, b, method, opts) {
  mv <- if (is.function(A)) A else function(x) as.numeric(A %*% x)
  n <- length(b)
  nx <- if (is.function(A)) n else ncol(A)
  tol <- max(opts$abstol %||% 0, (opts$reltol %||% sqrt(.Machine$double.eps)) * sqrt(sum(b^2)))
  maxit <- opts$maxiters %||% (10L * n)
  Pl <- .precond(opts$Pl); Pr <- .precond(opts$Pr)
  x <- as.numeric(opts$u0 %||% numeric(nx))
  nrm <- function(v) sqrt(sum(v^2))
  it <- 0L
  if (method %in% c("cg")) {
    r <- b - mv(x); z <- Pl(r); p <- z; rz <- sum(r * z)
    while (nrm(r) > tol && it < maxit) {
      Ap <- mv(p); a <- rz / sum(p * Ap)
      x <- x + a * p; r <- r - a * Ap
      z <- Pl(r); rz2 <- sum(r * z)
      p <- z + (rz2 / rz) * p; rz <- rz2
      it <- it + 1L
    }
  } else if (method == "bicgstab") {
    r <- b - mv(x); rh <- r; rho <- 1; al <- 1; w <- 1; v <- numeric(n); p <- numeric(n)
    while (nrm(r) > tol && it < maxit) {
      rho2 <- sum(rh * r)
      if (rho2 == 0) break
      p <- if (it == 0L) r else r + (rho2 / rho) * (al / w) * (p - w * v)
      ph <- Pl(p); v <- mv(ph); al <- rho2 / sum(rh * v)
      s <- r - al * v
      if (nrm(s) <= tol) { x <- x + al * ph; r <- s; it <- it + 1L; break }
      sh <- Pl(s); tt <- mv(sh); w <- sum(tt * s) / sum(tt * tt)
      x <- x + al * ph + w * sh; r <- s - w * tt
      rho <- rho2
      it <- it + 1L
    }
  } else if (method %in% c("gmres", "fgmres", "minres")) {
    # restarted GMRES with right preconditioning (MINRES, for symmetric A,
    # reaches the same minimal residual iterate)
    restart <- opts$restart %||% min(n, 50L)
    repeat {
      r <- b - mv(x)
      beta <- nrm(r)
      if (beta <= tol || it >= maxit) break
      V <- matrix(0, n, restart + 1L); H <- matrix(0, restart + 1L, restart)
      Z <- matrix(0, n, restart)
      V[, 1L] <- r / beta
      k_used <- 0L
      for (k in seq_len(restart)) {
        Z[, k] <- Pr(V[, k])
        w <- mv(Z[, k])
        for (j in seq_len(k)) {
          H[j, k] <- sum(w * V[, j])
          w <- w - H[j, k] * V[, j]
        }
        H[k + 1L, k] <- nrm(w)
        k_used <- k
        it <- it + 1L
        e1 <- c(beta, numeric(k))
        y <- qr.coef(qr(H[seq_len(k + 1L), seq_len(k), drop = FALSE]), e1)
        res <- nrm(e1 - H[seq_len(k + 1L), seq_len(k), drop = FALSE] %*% y)
        if (H[k + 1L, k] < 1e-14 || res <= tol || it >= maxit) break
        V[, k + 1L] <- w / H[k + 1L, k]
      }
      e1 <- c(beta, numeric(k_used))
      y <- qr.coef(qr(H[seq_len(k_used + 1L), seq_len(k_used), drop = FALSE]), e1)
      y[is.na(y)] <- 0
      x <- x + as.numeric(Z[, seq_len(k_used), drop = FALSE] %*% y)
    }
  } else if (method == "lsmr") {
    # least squares by CGLS
    r <- b - mv(x); s <- as.numeric(crossprod(A, r)); p <- s; g <- sum(s^2)
    while (sqrt(g) > tol * 1e-3 && it < maxit) {
      q <- mv(p); a <- g / sum(q^2)
      x <- x + a * p; r <- r - a * q
      s <- as.numeric(crossprod(A, r)); g2 <- sum(s^2)
      p <- s + (g2 / g) * p; g <- g2
      it <- it + 1L
    }
  } else if (method == "craigmr") {
    # minimum norm solution by CGNE: x = A' y with A A' y = b
    y <- numeric(n)
    AAt <- function(v) mv(as.numeric(crossprod(A, v)))
    r <- b - AAt(y); p <- r; rr <- sum(r^2)
    while (sqrt(rr) > tol && it < maxit) {
      q <- AAt(p); a <- rr / sum(p * q)
      y <- y + a * p; r <- r - a * q
      rr2 <- sum(r^2); p <- r + (rr2 / rr) * p; rr <- rr2
      it <- it + 1L
    }
    x <- as.numeric(crossprod(A, y))
  } else stop("unknown Krylov method ", method, call. = FALSE)
  list(u = x, iters = it, resid = b - mv(x))
}

.lin_solve <- function(prob, alg, opts, cache = new.env(parent = emptyenv())) {
  A <- prob$A; b <- as.numeric(prob$b)
  if (inherits(A, "LowRankUpdatedMatrix")) {
    # Woodbury: (A + U C V')^-1 b
    Ainv_b <- solve(A$A, b); Ainv_U <- solve(A$A, A$U)
    cap <- solve(A$C) + crossprod(A$V, Ainv_U)
    u <- as.numeric(Ainv_b - Ainv_U %*% solve(cap, crossprod(A$V, Ainv_b)))
    full <- A$A + A$U %*% A$C %*% t(A$V)
    return(list(u = u, resid = as.numeric(full %*% u - b), iters = 1L, alg = alg %||% LUFactorization()))
  }
  A <- .lin_prepare_A(A, opts, cache)
  if (!is.matrix(A) && !is.function(A)) A <- matrix(A, nrow = length(b))
  alg <- alg %||% .lin_default(A, opts$assumptions %||% OperatorAssumptions())
  if (identical(alg$method, "function")) {
    u <- as.numeric(alg$fun(A, b, opts$u0 %||% numeric(length(b)), prob$p))
    return(list(u = u, resid = as.numeric(A %*% u - b), iters = 1L, alg = alg))
  }
  if (identical(alg$method, "krylov")) {
    kopts <- utils::modifyList(alg[setdiff(names(alg), c("name", "method", "krylov"))], opts)
    r <- .krylov(A, b, alg$krylov, .lin_warm_start(alg$warm_start, A, b, cache$u_prev, kopts))
    cache$u_prev <- r$u
    return(c(r, list(alg = alg)))
  }
  u <- .lin_direct(A, b, alg, cache)
  list(u = as.numeric(u), resid = as.numeric(A %*% u - b), iters = 1L, alg = alg)
}

#' @export
init.LinearProblem <- function(prob, alg = NULL, ...) {
  cache <- new.env(parent = emptyenv())
  cache$A <- prob$A; cache$b <- prob$b; cache$alg <- alg; cache$opts <- list(...); cache$prob <- prob
  class(cache) <- "LinearCache"
  cache
}

.lin_solve_cache <- function(cache) {
  pr <- cache$prob
  pr$A <- cache$A; pr$b <- cache$b
  r <- .lin_solve(pr, cache$alg, cache$opts, cache)
  .lin_solution(pr, r)
}

.lin_solution <- function(prob, r) {
  ok <- all(is.finite(r$u)) && max(abs(r$resid), 0) <= max(1e-6, 1e-8 * max(1, max(abs(prob$b), 0)))
  if (!is.matrix(prob$A) || nrow(prob$A) > ncol(prob$A)) ok <- all(is.finite(r$u))
  if (!is.null(prob$sys)) prob$u0 <- stats::setNames(r$u, .nl_names(unknowns(prob$sys)))
  .nl_new_solution(c("LinearSolution", "NonlinearSolution"), prob, r$alg, r$u, if (ok) "Success" else "Failure",
                   resid = r$resid, stats = list(iters = r$iters))
}

#' @export
print.LinearCache <- function(x, ...) {
  cat("LinearCache (", x$alg$name %||% "default", ")\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Eigenvalue problems

#' Eigenvalue problems
#'
#' `EigenvalueProblem(A, nev, which, B)` asks for `nev` eigenvalues of `A`
#' (of the generalized problem `A v = lambda B v` with `B`), the ones
#' selected by `which` (an `EigenvalueTarget`: `"LM"` largest magnitude,
#' `"SM"` smallest magnitude, `"LR"`/`"SR"` largest or smallest real part,
#' `"LI"`/`"SI"` imaginary part). `solve(prob, alg)` returns an
#' `EigenvalueSolution` with `values` and `vectors`. `DenseEigen()` uses
#' LAPACK through [eigen()]; `ArnoldiMethodJL()` (also `ArnoldiMethod()`),
#' `ArpackJL()`, `KrylovKitEigen()` and `JacobiDavidsonJL()` run an Arnoldi
#' iteration with explicit restarts in R for large matrices and the dense
#' method for small ones.
#'
#' @param A,B Matrices.
#' @param nev Number of eigenvalues.
#' @param which The `EigenvalueTarget`.
#' @param krylovdim Dimension of the Krylov space of the Arnoldi iteration.
#' @param values,vectors,retcode Fields of an `EigenvalueSolution`.
#' @param ... Options.
#' @return A problem, an algorithm, or the targets.
#' @examples
#' A <- diag(1:6)
#' solve(EigenvalueProblem(A, nev = 2), ArnoldiMethodJL())$values
#' @name eigenvalue_problems
NULL

#' @rdname eigenvalue_problems
#' @export
EigenvalueTarget <- list(LM = "LM", SM = "SM", LR = "LR", SR = "SR", LI = "LI", SI = "SI")

#' @rdname eigenvalue_problems
#' @export
EigenvalueProblem <- function(A, nev = 1L, which = EigenvalueTarget$LM, B = NULL, ...)
  structure(list(A = A, B = B, nev = as.integer(nev), which = which, kwargs = list(...)),
            class = c("EigenvalueProblem", "SciMLProblem"))

.eig_alg <- function(name) structure(list(name = name), class = c(name, "EigenvalueAlgorithm"))
#' @rdname eigenvalue_problems
#' @export
DenseEigen <- function() .eig_alg("DenseEigen")
#' @rdname eigenvalue_problems
#' @export
ArnoldiMethodJL <- function(krylovdim = 30L) structure(list(name = "ArnoldiMethodJL", krylovdim = krylovdim),
                                                       class = c("ArnoldiMethodJL", "EigenvalueAlgorithm"))
#' @rdname eigenvalue_problems
#' @export
ArnoldiMethod <- function(krylovdim = 30L) ArnoldiMethodJL(krylovdim)
#' @rdname eigenvalue_problems
#' @export
ArpackJL <- function(krylovdim = 30L) ArnoldiMethodJL(krylovdim)
#' @rdname eigenvalue_problems
#' @export
KrylovKitEigen <- function(krylovdim = 30L) ArnoldiMethodJL(krylovdim)
#' @rdname eigenvalue_problems
#' @export
JacobiDavidsonJL <- function(krylovdim = 30L) ArnoldiMethodJL(krylovdim)

#' @rdname eigenvalue_problems
#' @export
EigenvalueSolution <- function(values, vectors = NULL, retcode = "Success")
  structure(list(values = values, vectors = vectors, retcode = retcode), class = "EigenvalueSolution")

.eig_order <- function(v, which) {
  key <- switch(which, LM = -Mod(v), SM = Mod(v), LR = -Re(v), SR = Re(v), LI = -Im(v), SI = Im(v), -Mod(v))
  order(key)
}

#' @export
solve.EigenvalueProblem <- function(a, b = NULL, ...) {
  A <- as.matrix(a$A)
  M <- if (is.null(a$B)) A else solve(as.matrix(a$B), A)
  n <- nrow(M)
  if (inherits(b, "ArnoldiMethodJL") && n > 2L * (b$krylovdim %||% 30L)) {
    r <- .arnoldi(M, a$nev, a$which, b$krylovdim)
    return(EigenvalueSolution(r$values, r$vectors))
  }
  e <- eigen(M)
  o <- .eig_order(e$values, a$which)[seq_len(min(a$nev, n))]
  EigenvalueSolution(e$values[o], e$vectors[, o, drop = FALSE])
}

# Arnoldi iteration with explicit restarts on the wanted Ritz vectors
.arnoldi <- function(M, nev, which, m) {
  n <- nrow(M)
  v <- stats::rnorm(n)
  for (restart in seq_len(50L)) {
    V <- matrix(0, n, m + 1L); H <- matrix(0, m + 1L, m)
    V[, 1L] <- v / sqrt(sum(v^2))
    for (k in seq_len(m)) {
      w <- as.numeric(M %*% V[, k])
      for (j in seq_len(k)) { H[j, k] <- sum(w * V[, j]); w <- w - H[j, k] * V[, j] }
      H[k + 1L, k] <- sqrt(sum(w^2))
      if (H[k + 1L, k] < 1e-14) { m <- k; break }
      V[, k + 1L] <- w / H[k + 1L, k]
    }
    e <- eigen(H[seq_len(m), seq_len(m), drop = FALSE])
    o <- .eig_order(e$values, which)[seq_len(min(nev, m))]
    X <- V[, seq_len(m), drop = FALSE] %*% e$vectors[, o, drop = FALSE]
    res <- abs(H[m + 1L, m] * e$vectors[m, o])
    if (all(res < 1e-10 * max(1, Mod(e$values[o])))) return(list(values = e$values[o], vectors = X))
    v <- Re(rowSums(X))
  }
  list(values = e$values[o], vectors = X)
}
