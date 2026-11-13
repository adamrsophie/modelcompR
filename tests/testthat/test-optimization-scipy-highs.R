test_that("the HiGHS simplex solver settles the status of linprog", {
  A <- rbind(c(1, 1, 1), c(2, 0.5, 0), c(0, 1, -1), c(1, -1, 2))
  s <- .hs_lp_solve(c(-1, -2, 0.5), A, c(-Inf, -Inf, -Inf, 0.5), c(4, 5, 1.5, 0.5), c(0, -1, -Inf), c(3, Inf, 2))
  expect_identical(s$status, "Optimal")
  expect_equal(s$col_value, c(1.5, 2, 0.5), tolerance = 1e-14)
  s <- .hs_lp_solve(c(-1, -2, 0.5), A, c(-Inf, -Inf, -Inf, 0.5), c(4, 5, 1.5, 0.5), c(0, -1, -Inf), c(3, Inf, 2),
                    list(iteration_limit = 2))
  expect_identical(s$status, "IterationLimit")
  expect_identical(s$nit, 2L)
  expect_identical(.hs_lp_solve(c(1, 1), rbind(c(1, 1)), -Inf, -1, c(0, 0), c(Inf, Inf))$status, "Infeasible")
  expect_identical(.hs_lp_solve(c(-1, 0), rbind(c(1, -1)), -Inf, 1, c(0, 0), c(Inf, Inf))$status, "Unbounded")
  # scaled LP (entries outside [0.2, 5])
  s <- .hs_lp_solve(c(-1, -1), rbind(c(100, 1), c(0.01, 1)), c(-Inf, -Inf), c(100, 1), c(0, 0), c(Inf, Inf))
  expect_identical(s$status, "Optimal")
  expect_equal(s$col_value, c(1 / (1 - 1e-4) * (1 - 1e-2), 1 - 0.01 * (1 / (1 - 1e-4) * (1 - 1e-2))), tolerance = 1e-12)
  # no rows
  s <- .hs_lp_solve(c(1, -1), matrix(0, 0, 2), numeric(), numeric(), c(0, -Inf), c(1, 3))
  expect_identical(s$status, "Optimal")
  expect_equal(s$col_value, c(0, 3))
})

test_that("the HiGHS LP presolve keeps the optimum of the original LP", {
  set.seed(11)
  for (t in 1:40) {
    m <- sample(2:7, 1); n <- sample(2:7, 1)
    A <- matrix(round(rnorm(m * n), 1) * (runif(m * n) < 0.6), m, n)
    for (i in which(rowSums(A != 0) == 0)) A[i, sample(n, 1)] <- 1
    cost <- round(rnorm(n), 1)
    if (m > 2) A[m, ] <- 2 * A[1, ]
    if (n > 2) {
      A[, n] <- A[, 1]; cost[n] <- cost[1]
    }
    cl <- ifelse(runif(n) < 0.7, round(runif(n, -3, 0), 1), -Inf)
    cu <- ifelse(runif(n) < 0.7 & is.finite(cl), round(cl + runif(n, 0, 4), 1), Inf)
    rl <- ifelse(runif(m) < 0.5, round(runif(m, -5, 0), 1), -Inf)
    ru <- ifelse(runif(m) < 0.6, round(runif(m, 0, 5), 1), Inf)
    ru[is.infinite(rl) & is.infinite(ru)] <- 3
    a <- .hs_highs_run(cost, A, rl, ru, cl, cu, list(presolve = FALSE))
    b <- .hs_highs_run(cost, A, rl, ru, cl, cu)
    if (a$status == "Optimal") {
      expect_identical(b$status, "Optimal")
      expect_equal(b$objective, a$objective, tolerance = 1e-9)
    }
  }
  # reduced to empty: nothing is left for the simplex solver
  A3 <- rbind(c(1, 1, 1, 1), c(1, -2, 0, 1), c(0, 1, 1, -1))
  s <- .hs_highs_run(c(1, -1, 2, 0.5), A3, c(-Inf, -Inf, 1), c(5, 2, 1), c(1, -2, -Inf, 0), c(4, 3, Inf, 2))
  expect_identical(s$status, "Optimal")
  expect_identical(s$nit, 0L)
  expect_equal(s$col_value, c(1, 3, -2, 0))
})
