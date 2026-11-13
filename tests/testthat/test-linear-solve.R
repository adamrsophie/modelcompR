test_that("direct factorizations solve square systems", {
  A <- matrix(c(4, 1, 2, 1, 3, 0, 2, 0, 5), 3)
  b <- c(1, 2, 3)
  ref <- solve(A, b)
  prob <- LinearProblem(A, b)
  for (alg in list(LUFactorization(), QRFactorization(), SVDFactorization(), CholeskyFactorization(),
                   LDLtFactorization(), BunchKaufmanFactorization(), NormalCholeskyFactorization(),
                   GenericFactorization(), SimpleLUFactorization(), DirectLdiv(), MKLLUFactorization(),
                   KLUFactorization(), LinearSolveFunction(function(A, b, u, p) solve(A, b)))) {
    sol <- solve(prob, alg)
    expect_equal(as.numeric(sol$u), ref, tolerance = 1e-8, info = alg$name)
    expect_identical(sol$retcode, "Success")
  }
  expect_equal(as.numeric(solve(LinearProblem(diag(c(2, 4)), c(2, 8)), DiagonalFactorization())$u), c(1, 2))
})

test_that("Krylov methods converge", {
  set.seed(3)
  n <- 30
  M <- matrix(stats::rnorm(n * n), n)
  S <- crossprod(M) + n * diag(n)
  b <- stats::rnorm(n)
  ref <- solve(S, b)
  prob <- LinearProblem(S, b)
  for (alg in list(KrylovJL_CG(), KrylovJL_GMRES(), KrylovJL_FGMRES(), KrylovJL_MINRES(), KrylovJL_BICGSTAB(),
                   SimpleGMRES(), IterativeSolversJL_IDRS(), KrylovJL(method = "cg")))
    expect_equal(as.numeric(solve(prob, alg, reltol = 1e-12)$u), ref, tolerance = 1e-6, info = alg$name)
  N <- M + n * diag(n)
  expect_equal(as.numeric(solve(LinearProblem(N, b), KrylovJL_GMRES(), reltol = 1e-12, restart = 5)$u),
               solve(N, b), tolerance = 1e-6)
  expect_equal(as.numeric(solve(prob, KrylovJL_CG(), Pl = InvPreconditioner(diag(1 / diag(S))), reltol = 1e-12)$u),
               ref, tolerance = 1e-6)
})

test_that("rectangular systems give least squares and minimum norm solutions", {
  set.seed(4)
  A <- matrix(stats::rnorm(20), 5, 4)
  b <- stats::rnorm(5)
  ls <- qr.coef(qr(A), b)
  expect_equal(as.numeric(solve(LinearProblem(A, b), QRFactorization())$u), ls, tolerance = 1e-8)
  expect_equal(as.numeric(solve(LinearProblem(A, b), KrylovJL_LSMR(), maxiters = 200)$u), ls, tolerance = 1e-6)
  W <- t(A)
  c4 <- stats::rnorm(4)
  mn <- as.numeric(t(W) %*% solve(W %*% t(W), c4))
  expect_equal(as.numeric(solve(LinearProblem(W, c4), KrylovJL_CRAIGMR(), maxiters = 200)$u), mn, tolerance = 1e-6)
})

test_that("default algorithm follows the operator assumptions", {
  sq <- diag(2); rect <- matrix(1, 3, 2)
  expect_identical(.lin_default(sq, OperatorAssumptions())$name, "LUFactorization")
  expect_identical(.lin_default(rect, OperatorAssumptions())$name, "QRFactorization")
  expect_identical(.lin_default(rect, OperatorAssumptions(condition = OperatorCondition$WellConditioned))$name,
                   "NormalCholeskyFactorization")
  expect_identical(.lin_default(sq, OperatorAssumptions(condition = OperatorCondition$VeryIllConditioned))$name,
                   "QRFactorization")
  expect_identical(.lin_default(sq, OperatorAssumptions(condition = OperatorCondition$SuperIllConditioned))$name,
                   "SVDFactorization")
  expect_setequal(names(OperatorCondition),
                  c("IllConditioned", "VeryIllConditioned", "SuperIllConditioned", "WellConditioned"))
  A <- matrix(c(2, 1, 1, 3), 2)
  sol <- solve(LinearProblem(A, c(1, 2)), assumptions = OperatorAssumptions(condition = "SuperIllConditioned"))
  expect_equal(as.numeric(sol$u), solve(A, c(1, 2)))
})

test_that("LinearCache reuses the factorization and updates b and A", {
  A <- matrix(c(4, 1, 1, 3), 2)
  cache <- init(LinearProblem(A, c(1, 2)), LUFactorization())
  expect_s3_class(cache, "LinearCache")
  expect_equal(as.numeric(solve_integrator(cache)$u), solve(A, c(1, 2)))
  f1 <- cache$fact
  cache$b <- c(2, 4)
  expect_equal(as.numeric(solve_integrator(cache)$u), solve(A, c(2, 4)))
  expect_identical(cache$fact, f1)
  cache$A <- 2 * A
  expect_equal(as.numeric(solve_integrator(cache)$u), solve(2 * A, c(2, 4)))
})

test_that("low rank updated matrices and the linear adjoint", {
  set.seed(5)
  A <- diag(4) * 3
  U <- matrix(stats::rnorm(8), 4); V <- matrix(stats::rnorm(8), 4); C <- diag(2)
  b <- stats::rnorm(4)
  full <- A + U %*% C %*% t(V)
  sol <- solve(LinearProblem(LowRankUpdatedMatrix(A, U, C, V), b), LUFactorization())
  expect_equal(as.numeric(sol$u), solve(full, b), tolerance = 1e-10)
  # gradient of g(u) = sum(u) with respect to b and A, against finite differences
  B <- matrix(c(4, 1, 0, 1, 3, 1, 0, 1, 2), 3); bb <- c(1, 2, 3)
  u <- solve(B, bb)
  g <- linear_adjoint(B, u, rep(1, 3))
  h <- 1e-6
  fd_b <- vapply(1:3, function(i) { e <- bb; e[i] <- e[i] + h; (sum(solve(B, e)) - sum(u)) / h }, 0)
  expect_equal(g$db, fd_b, tolerance = 1e-4)
  E <- B; E[2, 3] <- E[2, 3] + h
  expect_equal(g$dA[2, 3], (sum(solve(E, bb)) - sum(u)) / h, tolerance = 1e-4)
  expect_s3_class(LinearSolveAdjoint(), "AdjointSensitivityAlgorithm")
})

test_that("eigenvalue problems", {
  A <- diag(1:6)
  expect_equal(Re(solve(EigenvalueProblem(A, nev = 2), DenseEigen())$values), c(6, 5))
  expect_equal(Re(solve(EigenvalueProblem(A, nev = 2, which = EigenvalueTarget$SM), DenseEigen())$values), c(1, 2))
  set.seed(6)
  n <- 200
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n)))
  S <- Q %*% diag(c(50, 40, seq(1, 10, length.out = n - 2))) %*% t(Q)
  sol <- solve(EigenvalueProblem(S, nev = 2), ArnoldiMethodJL(krylovdim = 40))
  expect_equal(sort(Re(sol$values), decreasing = TRUE), c(50, 40), tolerance = 1e-6)
  v <- sol$vectors[, 1]
  expect_lt(max(abs(S %*% v - sol$values[1] * v)) / max(abs(v)), 1e-5)
  gen <- solve(EigenvalueProblem(diag(c(2, 6)), nev = 1, B = diag(c(1, 2))), DenseEigen())
  expect_equal(Re(gen$values), 3)
})

test_that("Krylov solves report their convergence, as LinearSolve's KrylovJL wrapper does", {
  # one CG step on an operator stops short of the solution c(1, 0.5)
  op <- function(x) c(x[1], 2 * x[2])
  s <- solve(LinearProblem(op, c(1, 1)), KrylovJL_CG(), maxiters = 1)
  expect_equal(s$retcode, "MaxIters")
  expect_equal(s$u, c(2, 2) / 3)
  s <- solve(LinearProblem(op, c(1, 1)), KrylovJL_CG())
  expect_equal(s$retcode, "Success")
  expect_equal(s$u, c(1, 0.5))
  # rectangular least squares: the iteration limit is reported, the
  # converged least squares solution is a success
  A <- matrix(c(1, 2, 3, 4, 5, 7), 3); b <- c(1, 0, 2)
  expect_equal(solve(LinearProblem(A, b), KrylovJL_LSMR(), maxiters = 1)$retcode, "MaxIters")
  s <- solve(LinearProblem(A, b), KrylovJL_LSMR())
  expect_equal(s$retcode, "Success")
  expect_equal(s$u, qr.coef(qr(A), b), tolerance = 1e-6)
  expect_equal(solve(LinearProblem(diag(3) + 0.1, c(1, 2, 3)), KrylovJL_GMRES(), maxiters = 1)$retcode,
               "MaxIters")
})
