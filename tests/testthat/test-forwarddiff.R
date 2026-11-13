test_that("forward mode derivatives match ForwardDiff", {
  # julia tools/forwarddiff-probe.jl (ForwardDiff 1.4.7,
  # DifferentiationInterface 0.7.22, OptimizationBase 5.10.1)
  f <- function(x) x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
  cn <- function(x) c(x[1] * x[2] * x[3] * x[4], sum(x * x))
  ad <- .ad_opt_derivatives(f, cn, 2L)
  x <- c(1, 5, 5, 1)
  expect_identical(ad$grad(x), c(12, 1, 2, 11))
  expect_identical(as.numeric(ad$hess(x)), c(2, 1, 1, 12, 1, 0, 0, 1, 1, 0, 0, 1, 12, 1, 1, 0))
  expect_identical(as.numeric(ad$jac(x)), c(25, 2, 5, 10, 5, 10, 25, 2))
  expect_identical(as.numeric(ad$cons_h(x)[[1L]]), c(0, 5, 5, 25, 5, 0, 1, 5, 5, 1, 0, 5, 25, 5, 5, 0))
  expect_identical(as.numeric(ad$lag_h(x, 1, c(0.5, -0.3))),
                   c(1.4, 3.5, 3.5, 24.5, 3.5, -0.6, 0.5, 3.5, 3.5, 0.5, -0.6, 3.5, 24.5, 3.5, 3.5, -0.6))
  # the matrix form with obj_weight 0 sums the weighted constraint Hessians
  expect_identical(as.numeric(ad$lag_h(x, 0, c(0.7, 0), weighted_zero = TRUE)),
                   c(0, 3.5, 3.5, 17.5, 3.5, 0, 0.7, 3.5, 3.5, 0.7, 0, 3.5, 17.5, 3.5, 3.5, 0))
  x <- c(1.3, 4.1, 3.7, 1.9)
  expect_equal(ad$grad(x), c(19.759999999999998, 2.4699999999999998, 3.4699999999999998, 11.83), tolerance = 1e-15)
  expect_equal(as.numeric(ad$jac(x)), c(28.822999999999997, 2.6, 9.139000000000001, 8.2, 10.126999999999999, 7.4,
                                        19.721, 3.8), tolerance = 1e-15)
  expect_equal(as.numeric(ad$lag_h(x, 0.25, c(1.5, 2.5))),
               c(5.95, 11.02, 12.159999999999998, 25.355, 11.02, 5, 3.7049999999999996, 7.540000000000001,
                 12.159999999999998, 3.7049999999999996, 5, 8.32, 25.355, 7.540000000000001, 8.32, 5),
               tolerance = 1e-15)
  # Rosenbrock with no constraints
  rd <- .ad_opt_derivatives(function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2, function(x) numeric(), 0L)
  expect_equal(rd$grad(c(-1.2, 1)), c(-215.59999999999997, -87.99999999999999), tolerance = 1e-15)
  expect_identical(as.numeric(rd$hess(c(-1.2, 1))), c(1330, 480, 480, 200))
  # 13 unknowns, two chunks of 7
  lf <- function(x) {
    s <- 0
    for (i in 1:13) s <- s + i * x[i]^2 * x[(i %% 13) + 1]
    s + exp(x[1] * x[13])
  }
  ld <- .ad_opt_derivatives(lf, function(x) c(sum(x) * x[2], x[7]^3 - x[12] * x[3]), 2L)
  xl <- 0.1 * (1:13)
  expect_equal(ld$grad(xl), c(23.49047689832201, 0.25000000000000006, 0.8000000000000003, 1.87, 3.6400000000000006,
                              6.290000000000001, 10.000000000000002, 14.950000000000003, 21.32, 29.29,
                              39.040000000000006, 50.75000000000001, 20.773882838332465), tolerance = 1e-15)
  expect_equal(ld$hess(xl), t(ld$hess(xl)))
})

test_that("DAE Jacobians come from dual numbers under AutoForwardDiff", {
  res <- function(du, u, p, t) c(du[1] - p * u[1] * u[2], u[1]^2 + u[2]^2 - 1 + du[2] * u[1])
  ctx <- list(f = res, p = 0.5, autodiff = AutoForwardDiff())
  du <- c(0.3, -0.2); u <- c(0.6, 0.8)
  jj <- .calc_J_dae(ctx, du, u, 0)
  expect_identical(jj$J, rbind(c(-0.5 * 0.8, -0.5 * 0.6), c(2 * 0.6 + -0.2, 2 * 0.8)))
  expect_identical(jj$Jdu, rbind(c(1, 0), c(0, 0.6)))
  expect_identical(jj$nf, 2L)
  # the time derivative of the Rosenbrock methods is ForwardDiff.derivative
  rc <- list(f = function(u, p, t) -1000 * (u - cos(t)), p = NULL, autodiff = AutoForwardDiff())
  expect_identical(.rosen_tgrad(rc, 0.5, 0.25)$dT, -1000 * sin(0.25))
  # finite differences under AutoFiniteDiff
  ctx$autodiff <- AutoFiniteDiff()
  expect_false(identical(.calc_J_dae(ctx, du, u, 0)$J, jj$J))
})
