test_that("Krylov warm starts reach the cold start solution", {
  set.seed(1)
  A <- diag(30) * 4 + matrix(runif(900, -0.1, 0.1), 30)
  b <- runif(30)
  ref <- solve(A, b)
  for (mode in c(WarmStart$None, WarmStart$Previous, WarmStart$Hegedus)) {
    cache <- init(LinearProblem(A, b), KrylovJL_GMRES(warm_start = mode))
    s1 <- solve_integrator(cache)
    cache$b <- b * 1.01
    s2 <- solve_integrator(cache)
    expect_equal(unname(s1$u), ref, tolerance = 1e-6)
    expect_equal(unname(s2$u), ref * 1.01, tolerance = 1e-6)
  }
  # a previous solution parallel to the new one: the Hegedus guess is exact
  cache <- init(LinearProblem(A, b), KrylovJL_GMRES(warm_start = WarmStart$Hegedus))
  solve_integrator(cache)
  cache$b <- 2 * b
  expect_lte(solve_integrator(cache)$stats$iters, 1L)
})

test_that("nonstructural zeros are recorded in the assumptions", {
  a <- OperatorAssumptions(TRUE, nonstructural_zeros = NonstructuralZeros$Present)
  expect_identical(a$nonstructural_zeros, "Present")
  expect_identical(OperatorAssumptions()$nonstructural_zeros, "Auto")
})

test_that("PartitionedSolversAlgorithm solves directly or with CG", {
  A <- diag(5) * 3 + 0.2
  b <- 1:5
  expect_equal(unname(solve(LinearProblem(A, b), PartitionedSolversAlgorithm())$u), solve(A, b), tolerance = 1e-12)
  expect_equal(unname(solve(LinearProblem(A, b), PartitionedSolversAlgorithm("cg"))$u), solve(A, b), tolerance = 1e-6)
})

test_that("update_gamma refreshes the shifted operator", {
  J <- matrix(c(-2, 1, 0.5, -3), 2)
  W <- WOperator(NULL, 0.5, J)
  cache <- init(LinearProblem(J, c(1, 2)), LUFactorization())
  cache$A <- W
  u1 <- solve_integrator(cache)$u
  expect_equal(unname(u1), solve(J - diag(2) / 0.5, c(1, 2)), tolerance = 1e-12)
  update_gamma(cache, 0.25)
  u2 <- solve_integrator(cache)$u
  expect_equal(unname(u2), solve(J - diag(2) / 0.25, c(1, 2)), tolerance = 1e-12)
  bad <- init(LinearProblem(J, c(1, 2)))
  expect_error(update_gamma(bad, 0.1), "WOperator")
})

test_that("supernodal panel solves match triangular solves", {
  set.seed(2)
  np <- 4
  W <- matrix(runif(36, 0.5, 1.5), 6) + diag(3, 6)
  B <- matrix(runif(12), 4)
  L <- W[1:np, 1:np]; L[upper.tri(L)] <- 0; diag(L) <- 1
  U <- W[1:np, 1:np]; U[lower.tri(U)] <- 0
  for (algo in c("auto", "kernel", "blas", "triangularsolve")) {
    expect_equal(supernodal_panel_solve(W, B, np, "lower", algo), solve(L, B), tolerance = 1e-12)
    expect_equal(supernodal_panel_solve(W, B, np, "factor_lower", algo), solve(L, B), tolerance = 1e-12)
    expect_equal(supernodal_panel_solve(W, B, np, "upper", algo), solve(U, B), tolerance = 1e-12)
  }
  R <- matrix(runif(12), 3)
  expect_equal(supernodal_panel_solve(W, R, np, "factor_right_upper"), R %*% solve(U), tolerance = 1e-12)
  expect_equal(supernodal_panel_solve_backend("kernel", W, B, np, "upper"), solve(U, B), tolerance = 1e-12)
  expect_error(supernodal_panel_solve(W, B, np, "diagonal"), "unknown supernodal panel operation")
  expect_true(slate_isavailable())
  expect_output(show_algorithm_choices(), "Default Algorithm Choices")
})
