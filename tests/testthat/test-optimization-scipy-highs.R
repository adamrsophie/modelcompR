# long running: run on CI (NOT_CRAN = true), not on CRAN
skip_on_cran()

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

test_that("postsolve resplits a duplicate column as DuplicateColumn::undoFix does", {
  # C round is half away from zero
  expect_identical(.hs_round(c(-2.5, -0.5, 0.5, 2.5, 2.4)), c(-3, -1, 1, 3, 2))
  r <- list(scale = 2, cl = 0, cu = 1, dl = 0, du = 1)
  expect_false(.hs_dupcol_split_error(r, 0.5, 0.25, 1, 1e-7))
  expect_true(.hs_dupcol_split_error(r, 1.5, 0, 1.5, 1e-7))
  # scale > 0: y at its upper bound leaves x below its lower bound, so x
  # goes to its lower bound and y takes the rest
  expect_identical(.hs_dupcol_undo_fix(r, 1.5, 1e-7), c(0, 0.75))
  # scale < 0: y at its lower bound
  expect_identical(.hs_dupcol_undo_fix(list(scale = -1, cl = 0, cu = 2, dl = 0, du = 1), 0.5, 1e-7), c(0.5, 0))
  # integer x: the first integer value of x that leaves y within its bounds
  expect_identical(.hs_dupcol_undo_fix(list(scale = 1, cl = 0, cu = 3, dl = 0, du = 0.8, ci = TRUE), 2.7, 1e-7),
                   c(2, 0x1.6666666666668p-1))
})

test_that("a rank deficient refactor backtracks as HEkk::getNonsingularInverse does", {
  ek <- .hs_new_ekk(c(-1, -1), matrix(1, 2, 2), c(-Inf, -Inf), c(1, 2), c(0, 0), c(Inf, Inf), list())
  .hs_slack_basis(ek)
  singular <- function(updates) {
    ek$basicIndex <- c(1L, 2L); ek$nonbasicFlag <- c(0L, 0L, 1L, 1L); ek$update_count <- updates
  }
  # a nonsingular refactor keeps the basis for backtracking
  expect_true(.hs_reinvert(ek))
  expect_identical(ek$bt$basicIndex, c(3L, 4L))
  expect_identical(ek$update_limit, 5000L)
  expect_false(ek$backtracking)
  # two equal basic columns: back to the slack basis with half the updates
  singular(6L)
  expect_true(.hs_reinvert(ek))
  expect_identical(ek$basicIndex, c(3L, 4L))
  expect_identical(ek$nonbasicFlag, c(1L, 1L, 0L, 0L))
  expect_identical(ek$update_limit, 3L)
  expect_true(ek$backtracking)
  # no backtracking after a single update, nor without a backtracking basis
  singular(1L)
  expect_false(.hs_reinvert(ek))
  ek$bt <- NULL
  singular(6L)
  expect_false(.hs_reinvert(ek))
})

test_that("lpKktCheck changes the model status as HiGHS does", {
  # min x1 + x2 s.t. x1 + x2 >= 1, x >= 0: optimal x = (1, 0), y = 1
  cost <- c(1, 1); A <- matrix(c(1, 1), 1)
  rl <- 1; ru <- Inf; cl <- c(0, 0); cu <- c(Inf, Inf)
  ok <- list(status = "Optimal", col_value = c(1, 0), row_value = 1, col_dual = c(0, 0), row_dual = 1,
             col_status = c("basic", "lower"))
  chk <- function(s, ...) .hs_lp_kkt_check(s, cost, A, rl, ru, cl, cu, list(), ...)
  r <- chk(ok)
  expect_identical(r$status, "Optimal")
  expect_true(r$primal_feasible)
  expect_equal(r$kkt$primal_dual_objective_error, 0)
  # a primal infeasibility of 1e-4 is 1000 times the tolerance
  bad <- ok; bad$col_value <- c(1, -1e-4)
  r <- chk(bad)
  expect_identical(r$status, "Unknown")
  expect_false(r$primal_feasible)
  expect_identical(r$kkt$num_primal_infeasibilities, 1L)
  # within 100 times the tolerances "Unknown" becomes "Optimal"
  unk <- ok; unk$status <- "Unknown"; unk$col_value <- c(1, -5e-6)
  expect_identical(chk(unk)$status, "Optimal")
  # a primal feasible point turns UnboundedOrInfeasible into Unbounded
  ub <- ok; ub$status <- "UnboundedOrInfeasible"
  expect_identical(chk(ub)$status, "Unbounded")
  ub$col_value <- c(1, -1)
  expect_identical(chk(ub)$status, "UnboundedOrInfeasible")
  # without a basis the residuals count: a row value off by 1e-3
  nb <- ok; nb$row_value <- 1 + 1e-3
  r <- chk(nb, basis_valid = FALSE)
  expect_identical(r$kkt$num_relative_primal_residual_errors, 1L)
  expect_identical(r$status, "Unknown")
  # the dual objective: a dual error moves it from the primal one
  du <- ok; du$row_dual <- 1.5
  r <- chk(du)
  expect_equal(r$kkt$dual_objective, 1.5)
  expect_equal(r$kkt$primal_dual_objective_error, 0.5 / 3.5)
})

test_that("an LP whose matrix has no nonzeros checks its rows", {
  expect_identical(.hs_highs_run_ipm(1, matrix(0, 1, 1), 1, Inf, 0, Inf)$status, "Infeasible")
  expect_identical(.hs_highs_run(1, matrix(0, 1, 1), 1, Inf, 0, Inf)$status, "Infeasible")
  expect_identical(.hs_highs_run(1, matrix(0, 1, 1), -1, Inf, 0, Inf)$status, "Optimal")
  s <- list(status = "Optimal", col_value = 1 + 1e-6, row_value = 1 + 1e-6, col_dual = 0, row_dual = 1,
            col_status = "basic")
  expect_identical(.hs_lp_kkt_check(s, 1, matrix(1), 1, 1, -Inf, Inf)$status, "Optimal")
  expect_identical(.hs_lp_kkt_check(s, 1, matrix(1), 1, 1, -Inf, Inf, list(kkt_tolerance = 1e-10))$status, "Unknown")
})
