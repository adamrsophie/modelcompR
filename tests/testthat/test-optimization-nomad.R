test_that("NOMADOpt reaches the minimizers of NOMAD", {
  # NOMAD_jll has no build for this platform, so the references are the
  # exact minimizers: unconstrained (10/7, -6/7), with x1 >= 1.5 at
  # (1.5, -0.875) and with x1 + x2 >= 1 at (7/4, -3/4)
  quad <- function(x, p) (x[1] - 1) * (x[1] - 1) + 2 * (x[2] + 0.5) * (x[2] + 0.5) + x[1] * x[2]
  s <- solve(OptimizationProblem(quad, c(3, 2), numeric()), NOMADOpt())
  expect_equal(unname(s$u), c(10 / 7, -6 / 7), tolerance = 1e-7)
  expect_equal(s$objective, -11 / 14, tolerance = 1e-12)
  expect_equal(s$retcode, "Default")
  s <- solve(OptimizationProblem(quad, c(3, 2), numeric(), lb = c(1.5, -2), ub = c(4, 4)), NOMADOpt())
  expect_equal(unname(s$u), c(1.5, -0.875), tolerance = 1e-7)
  pc <- OptimizationProblem(quad, c(0, 0), numeric(), cons = function(x, p) x[1] + x[2], lcons = 1, ucons = Inf)
  for (cm in c("ExtremeBarrierMethod", "ProgressiveBarrierMethod")) {
    s <- solve(pc, NOMADOpt(), cons_method = cm)
    expect_equal(unname(s$u), c(1.75, -0.75), tolerance = 1e-6)
    expect_gte(sum(s$u), 1)
  }
  s <- solve(OptimizationProblem(quad, c(3, 2), numeric()), NOMADOpt(), maxiters = 25)
  expect_identical(s$original$nb_evals, 25L)
  expect_error(solve(pc, NOMADOpt(), bogus = 1), "has no field bogus")
})
