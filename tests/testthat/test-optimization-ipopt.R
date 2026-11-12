ip_hs071 <- function() {
  OptimizationProblem(function(x, p) x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3], c(1, 5, 5, 1), numeric(),
                      cons = function(x, p) c(x[1] * x[2] * x[3] * x[4], sum(x * x)), lcons = c(25, 40),
                      ucons = c(Inf, 40), lb = rep(1, 4), ub = rep(5, 4))
}

test_that("IpoptOptimizer and MadNLPOptimizer reach the solution of Ipopt", {
  # julia tools/ipopt-probe.jl (Ipopt.jl, exact derivatives): HS071
  ref_x <- c(0.9999999923235379, 4.742999641809297, 3.8211499817883072, 1.3794082897556985)
  ref_f <- 17.014017145179164
  ref_mult <- c(-0.5522936588816064, 0.16146856313488814)
  for (alg in list(IpoptOptimizer(), MadNLPOptimizer())) {
    s <- solve(ip_hs071(), alg)
    expect_equal(s$retcode, "Success")
    expect_equal(unname(s$u), ref_x, tolerance = 1e-7)
    expect_equal(s$objective, ref_f, tolerance = 1e-8)
    expect_equal(s$original$lambda, ref_mult, tolerance = 1e-5)
  }
  expect_s3_class(IpoptOptimizer(linear_solver = "ma57"), "IpoptOptimizer")
  expect_identical(MadNLPOptimizer()$acceptable_iter, 15L)
  expect_error(MadNLPOptimizer(mu_init = 0.2), "has no field mu_init")
  # additional_options max_iter wins over maxiters, as in OptimizationIpopt.jl
  s <- solve(ip_hs071(), IpoptOptimizer(additional_options = list(max_iter = 2L)), maxiters = 100)
  expect_identical(s$original$iterations, 2L)
})

test_that("the interior point optimizers handle bounds, equalities and the sense", {
  rb <- OptimizationProblem(function(x, p) { a <- 1 - x[1]; b <- x[2] - x[1] * x[1]; a * a + 100 * b * b },
                            c(-1.2, 1), numeric(), lb = c(-2, -2), ub = c(0.8, 2))
  s <- solve(rb, IpoptOptimizer())
  expect_equal(s$retcode, "Success")
  expect_equal(unname(s$u), c(0.8, 0.64), tolerance = 1e-7)
  eqp <- OptimizationProblem(function(x, p) sum((x - c(1, 2, 3))^2), c(0, 0, 0), numeric(),
                             cons = function(x, p) c(x[1] + x[2] + x[3], x[1] - x[2]), lcons = c(3, 0),
                             ucons = c(3, 0))
  s <- solve(eqp, MadNLPOptimizer(), abstol = 1e-10)
  expect_equal(unname(s$u), c(0.5, 0.5, 2), tolerance = 1e-8)
  expect_equal(s$original$lambda, c(2, -1), tolerance = 1e-6)
  mx <- OptimizationProblem(function(x, p) -(x[1] - 1) * (x[1] - 1) - (x[2] + 0.5) * (x[2] + 0.5), c(0, 0),
                            numeric(), sense = MaxSense)
  s <- solve(mx, IpoptOptimizer(), maxiters = 50)
  expect_equal(unname(s$u), c(1, -0.5), tolerance = 1e-8)
  expect_equal(s$objective, 0, tolerance = 1e-12)
})
