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

# HS071 with exact derivatives, as tools/ipopt-options-probe.jl passes them
ip_hs071_exact <- function(u0 = c(1, 5, 5, 1)) {
  gr <- function(x, p) c(x[1] * x[4] + x[4] * (x[1] + x[2] + x[3]), x[1] * x[4], x[1] * x[4] + 1,
                         x[1] * (x[1] + x[2] + x[3]))
  he <- function(x, p) matrix(c(2 * x[4], x[4], x[4], 2 * x[1] + x[2] + x[3], x[4], 0, 0, x[1], x[4], 0, 0, x[1],
                                2 * x[1] + x[2] + x[3], x[1], x[1], 0), 4, 4)
  cj <- function(x, p) rbind(c(x[2] * x[3] * x[4], x[1] * x[3] * x[4], x[1] * x[2] * x[4], x[1] * x[2] * x[3]),
                             2 * x)
  ch <- function(x, p) list(matrix(c(0, x[3] * x[4], x[2] * x[4], x[2] * x[3], x[3] * x[4], 0, x[1] * x[4],
                                     x[1] * x[3], x[2] * x[4], x[1] * x[4], 0, x[1] * x[2], x[2] * x[3],
                                     x[1] * x[3], x[1] * x[2], 0), 4, 4), diag(2, 4))
  fn <- OptimizationFunction(function(x, p) x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3], grad = gr, hess = he,
                             cons = function(x, p) c(x[1] * x[2] * x[3] * x[4], sum(x * x)), cons_j = cj,
                             cons_h = ch)
  OptimizationProblem(fn, u0, numeric(), lcons = c(25, 40), ucons = c(Inf, 40), lb = rep(1, 4), ub = rep(5, 4))
}

test_that("the options of IpoptOptimizer change the iterates as in Ipopt", {
  # julia tools/ipopt-options-probe.jl (Ipopt 3.14.19): iterations and objective
  run <- function(alg, u0 = c(1, 5, 5, 1)) solve(ip_hs071_exact(u0), alg)
  # the far start needs Ipopt's inertia correction (lg(rg) 2.0 in iteration 3)
  s <- run(IpoptOptimizer(), c(5, 1, 1, 5))
  expect_identical(s$original$iterations, 36L)
  expect_equal(s$objective, 27.146427593318933, tolerance = 1e-10)
  s <- run(IpoptOptimizer(warm_start_init_point = "yes"))
  expect_identical(s$original$iterations, 8L)
  expect_equal(s$objective, 17.014017145175, tolerance = 1e-10)
  s <- run(IpoptOptimizer(honor_original_bounds = "yes"))
  expect_equal(s$objective, 17.014017145179164, tolerance = 1e-10)
  # the bound relaxation is capped by constr_viol_tol
  s <- run(IpoptOptimizer(constr_viol_tol = 1e-12, additional_options = list(tol = 1e-3)))
  expect_identical(s$original$iterations, 8L)
  expect_equal(s$objective, 17.014035470693745, tolerance = 1e-12)
  for (cs in list(list("probing", 10L, 17.014017140224137), list("loqo", 9L, 17.014017140320846),
                  list("quality-function", 8L, 17.014017140224134))) {
    s <- run(IpoptOptimizer(mu_strategy = "adaptive", mu_oracle = cs[[1L]]))
    expect_identical(s$original$iterations, cs[[2L]])
    expect_equal(s$objective, cs[[3L]], tolerance = 1e-9)
  }
  # a dual infeasibility tolerance below reach ends with a tiny step
  # (Search_Direction_Becomes_Too_Small, the retcode Infeasible)
  s <- run(IpoptOptimizer(dual_inf_tol = 1e-12, additional_options = list(tol = 1e-3)))
  expect_identical(s$original$iterations, 11L)
  expect_identical(s$retcode, "Infeasible")
  s <- run(IpoptOptimizer(hessian_approximation = "limited-memory", mu_strategy = "adaptive"))
  expect_identical(s$retcode, "Success")
  expect_equal(s$objective, 17.014017140224176, tolerance = 1e-9)
  # an HSL solver is not in Ipopt_jll: Invalid_Option before iteration 0
  s <- run(IpoptOptimizer(linear_solver = "ma27"))
  expect_identical(s$retcode, "Failure")
  expect_identical(s$original$iterations, 0L)
})

test_that("the restoration phase of Ipopt detects a locally infeasible problem", {
  pc <- OptimizationProblem(function(x, p) x[1]^2 + x[2]^2, c(0, 0), numeric(),
                            cons = function(x, p) c(x[1]^2 + x[2]^2, x[2] - x[1]^2), lcons = c(1, 2),
                            ucons = c(1, 2))
  # julia tools/ipopt-resto-probe.jl (Ipopt.jl 1.16.0, Ipopt 3.14.19):
  # LOCALLY_INFEASIBLE after 6 iterations at x = (0, 1)
  s <- solve(pc, IpoptOptimizer())
  expect_identical(s$retcode, "Infeasible")
  expect_identical(s$original$restorations, 1L)
  expect_identical(s$original$iterations, 6L)
  expect_equal(unname(s$u), c(0, 1.000000000055827), tolerance = 1e-8)
  # Waechter and Biegler's example: two restoration phases, the second ends
  # in local infeasibility after 17 iterations at x1 = -1
  wb <- OptimizationProblem(function(x, p) x[1], c(-2, 1, 1), numeric(),
                            cons = function(x, p) c(x[1]^2 - x[2] - 1, x[1] - x[3] - 0.5), lcons = c(0, 0),
                            ucons = c(0, 0), lb = c(-Inf, 0, 0), ub = c(Inf, Inf, Inf))
  s <- solve(wb, IpoptOptimizer())
  expect_identical(s$retcode, "Infeasible")
  expect_identical(s$original$restorations, 2L)
  expect_identical(s$original$iterations, 17L)
  expect_equal(unname(s$u), c(-0.999999995000303, -9.99818181818682e-9, -9.999090909090908e-9),
               tolerance = 1e-8)
})
