nl_poly3 <- function(x, p) {
  a <- x[1] - 1; b <- x[2] - 2; c <- x[3] + 1
  a * a + 2 * b * b + 3 * c * c + a * a * b * b + 0.5 * a * c
}
nl_rosen <- function(x, p) { a <- 1 - x[1]; b <- x[2] - x[1] * x[1]; a * a + 100 * b * b }
nl_roseng <- function(x, p) { a <- 1 - x[1]; b <- x[2] - x[1] * x[1]; c(-2 * a - 400 * x[1] * b, 200 * b) }

test_that("Algorithm, NLopt and Opt follow the NLopt interface", {
  expect_identical(Algorithm(11)$name, "LD_LBFGS")
  expect_identical(Algorithm("NLOPT_LN_COBYLA")$value, 25L)
  expect_identical(NLopt$GN_AGS()$value, 43L)
  expect_length(NLopt, 44L)
  expect_error(Algorithm("LD_FOO"), "unknown NLopt algorithm")
  o <- Opt("LN_NELDERMEAD", 3)
  o$xtol_rel <- 1e-6
  expect_identical(o$algorithm$name, "LN_NELDERMEAD")
  expect_error(solve(OptimizationProblem(nl_rosen, c(-1.2, 1), numeric()), o), "dimension does not match")
  expect_error(solve(OptimizationProblem(nl_rosen, c(-1.2, 1), numeric()), NLopt$LD_LBFGS()), "requires gradients")
  expect_error(solve(OptimizationProblem(nl_rosen, c(-1.2, 1), numeric()), NLopt$GN_CRS2_LM()), "not yet ported")
})

test_that("the NLopt algorithms reach the solutions of NLopt", {
  # julia tools/nlopt-probe.jl (NLopt.jl); the native kernels take other
  # steps, so the points agree and the evaluation counts do not
  o <- Opt("LN_NELDERMEAD", 3)
  o$xtol_rel <- 1e-10
  o$maxeval <- 2000L
  s <- solve(OptimizationProblem(nl_poly3, c(0, 0, 0), numeric()), o)
  expect_equal(s$retcode, "Success")
  expect_equal(unname(s$u), c(1.00000000009593, 2.0000000000848477, -0.9999999999917069), tolerance = 1e-9)
  expect_gt(o$numevals, 0L)
  fr <- OptimizationFunction(nl_rosen, grad = nl_roseng)
  s <- solve(OptimizationProblem(fr, c(-1.2, 1), numeric()), NLopt$LD_LBFGS(), reltol = 1e-14)
  expect_equal(unname(s$u), c(0.9999999999928465, 0.9999999999852197), tolerance = 1e-9)
  s <- solve(OptimizationProblem(fr, c(-1.2, 1), numeric(), lb = c(-2, -2), ub = c(0.8, 2)), Algorithm("LD_LBFGS"),
             reltol = 1e-14)
  expect_equal(unname(s$u), c(0.8, 0.6400000000000001), tolerance = 1e-12)
  expect_equal(s$objective, 0.03999999999999998, tolerance = 1e-12)
  o <- Opt(NLopt$LN_COBYLA(), 3)
  o$xtol_rel <- 1e-10
  o$maxeval <- 5000L
  pc <- OptimizationProblem(nl_poly3, c(0, 0, 0), numeric(), cons = function(x, p) c(sum(x * x), sum(x)),
                            lcons = c(-Inf, 1), ucons = c(4, 1))
  s <- solve(pc, o)
  expect_equal(s$retcode, "Success")
  expect_equal(unname(s$u), c(0.450305841246599, 1.6250586952457444, -1.075364536492343), tolerance = 1e-7)
  expect_equal(s$objective, 0.6635572602138678, tolerance = 1e-7)
})
