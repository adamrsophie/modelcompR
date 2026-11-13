# GlobalAdjoint of GlobalDiffEq 1.7.3 (src/adjoint.jl and the
# SciMLSensitivity extension). The Julia values of the same problems are
# printed by the adjoint section of tools/globaldiffeq-probe.jl; here the
# estimates are checked against the known global error, as in GlobalDiffEq's
# test/adjoint_tests.jl. One state with one sample, or as many samples as
# states, makes the estimate independent of the random directions.

gea_lin <- function(...) ODEProblem(function(u, p, t) p * u, 1, c(0, 2), 2, ...)
gea_exact <- function(t) exp(2 * t)

test_that("GlobalAdjoint options and validation", {
  alg <- GlobalAdjoint(Tsit5())
  expect_identical(alg$scope, EndpointError)
  expect_identical(alg$control, ToleranceRefinement)
  expect_identical(GlobalAdjoint(Tsit5(), scope = TrajectoryError)$scope, TrajectoryError)
  expect_identical(GlobalAdjoint(Tsit5(), control = StepGridRefinement)$control, StepGridRefinement)
  expect_output(print(StepGridRefinement), "StepGridRefinement::GlobalAdjointControl = 1", fixed = TRUE)
  expect_identical(GlobalAdjointScope, list(EndpointError = EndpointError, TrajectoryError = TrajectoryError))
  expect_true(has_global_error(alg))
  expect_false(allows_arbitrary_number_types(alg))
  expect_false(allowscomplex(alg))
  expect_identical(alg$options$adjoint_abstol, 1e-10)
  expect_identical(GlobalAdjoint(Tsit5(), gtol = 1e-4)$options$adjoint_reltol, 1e-8)
  expect_equal(GlobalAdjoint(Tsit5(), gtol = 1e-9)$options$adjoint_abstol, 1e-11)
  expect_error(GlobalAdjoint(Tsit5(), samples = 0), "samples must be positive")
  expect_error(GlobalAdjoint(Tsit5(), gtol = -1), "gtol must be a positive finite real number")
  expect_error(GlobalAdjoint(Tsit5(), rng = 1), "rng must be a function")
  expect_error(solve(gea_lin(), alg), "GlobalAdjoint requires a positive `gtol` constructor keyword", fixed = TRUE)
  expect_error(GlobalAdjoint(Tsit5(), adjoint_abstol = 0, adjoint_reltol = 0),
               "adjoint_abstol and adjoint_reltol cannot both be zero")
  expect_error(adjoint_error_estimate(ODEProblem(function(u, p, t) -u, 1, c(1, 0)), alg),
               "GlobalAdjoint requires a forward-time ODEProblem")
})

test_that("the random directions are orthonormal", {
  d <- .gda_directions(c(1, 2, 3), 2, function(n) c(1, 2, 2)[seq_len(n)] + runif(n))
  expect_length(d, 2L)
  expect_equal(sum(d[[1]] * d[[1]]), 1)
  expect_equal(sum(d[[2]] * d[[2]]), 1)
  expect_equal(sum(d[[1]] * d[[2]]), 0, tolerance = 1e-14)
  expect_length(.gda_directions(c(1, 2), 5, stats::rnorm), 2L)
  expect_error(.gda_directions(c(1, 2), 1, function(n) numeric(n)), "failed to generate orthogonal random directions")
  expect_identical(c(.gda_sphere_expectation(1), .gda_sphere_expectation(2), .gda_sphere_expectation(3)),
                   c(1, 2 / pi, 1 / 2))
})

test_that("adjoint_error_estimate tracks the true global error", {
  est <- adjoint_error_estimate(gea_lin(), GlobalAdjoint(Tsit5(), samples = 1), abstol = 1e-4, reltol = 1e-4)
  s <- solve(gea_lin(), Tsit5(), abstol = 1e-4, reltol = 1e-4)
  expect_equal(est / abs(s$u[[length(s$t)]] - gea_exact(2)), 1, tolerance = 0.01)
  # a benign callback is not passed to the adjoint
  cb <- DiscreteCallback(function(u, t, integ) FALSE, function(integ) NULL)
  expect_equal(adjoint_error_estimate(gea_lin(callback = cb), GlobalAdjoint(Tsit5(), samples = 1),
                                      abstol = 1e-4, reltol = 1e-4), est, tolerance = 1e-6)
  # two states, two samples: the norm of the full error estimate
  lv <- function(u, p, t) c(1.5 * u[1] - u[1] * u[2], -3 * u[2] + u[1] * u[2])
  lvjac <- function(u, p, t) matrix(c(1.5 - u[2], u[2], -u[1], -3 + u[1]), 2, 2)
  prob <- ODEProblem(ODEFunction(lv, jac_fn = lvjac), c(1, 1), c(0, 2))
  ref <- solve(prob, Vern9(), abstol = 1e-13, reltol = 1e-13)
  s <- solve(prob, Tsit5(), abstol = 1e-5, reltol = 1e-5)
  est <- adjoint_error_estimate(prob, GlobalAdjoint(Tsit5()), abstol = 1e-5, reltol = 1e-5)
  expect_equal(est / sqrt(sum((s$u[[length(s$t)]] - ref$u[[length(ref$t)]])^2)), 1, tolerance = 0.02)
  # a callable parameter
  fprob <- ODEProblem(function(u, p, t) -u + p(t), 0, c(0, 1), function(t) 0.1 * sin(t))
  est <- adjoint_error_estimate(fprob, GlobalAdjoint(Tsit5(), samples = 1), abstol = 1e-4, reltol = 1e-4)
  s <- solve(fprob, Tsit5(), abstol = 1e-4, reltol = 1e-4)
  expect_equal(est / abs(s$u[[length(s$t)]] - 0.05 * (sin(1) - cos(1) + exp(-1))), 1, tolerance = 0.05)
})

test_that("GlobalAdjoint controls the endpoint error", {
  sol <- solve(gea_lin(), GlobalAdjoint(Tsit5(), gtol = 1e-6, samples = 1), abstol = 1e-3, reltol = 1e-3)
  n <- length(sol$t)
  expect_true(successful_retcode(sol))
  expect_lte(abs(sol$u[[n]] - gea_exact(2)), 1e-6)
  expect_identical(sol$prob$p, 2)
  expect_length(sol$global_error, n)
  expect_gt(sol$global_error[n], 0)
  expect_lte(sol$global_error[n], 1e-6)
  expect_identical(sol$global_error[-n], numeric(n - 1L))
  # the fixed step re-solve on the accepted grid
  grid <- solve(gea_lin(), GlobalAdjoint(Tsit5(), gtol = 1e-5, control = StepGridRefinement, samples = 1),
                abstol = 1e-3, reltol = 1e-3)
  n <- length(grid$t)
  expect_true(successful_retcode(grid))
  expect_identical(grid$t[n], 2)
  expect_lte(abs(grid$u[[n]] - gea_exact(2)), 1e-5)
  expect_lte(grid$global_error[n], 1e-5)
})

test_that("TrajectoryError estimates the error at every saved time", {
  prob <- ODEProblem(function(u, p, t) p * u, 1, c(0, 1), 2)
  sol <- solve(prob, GlobalAdjoint(Tsit5(), gtol = 1e-3, scope = TrajectoryError, samples = 1),
               abstol = 1e-3, reltol = 1e-3)
  n <- length(sol$t)
  expect_length(sol$global_error, n)
  expect_identical(sol$global_error[1], 0)
  expect_lte(sol$global_error[n], 1e-3)
  # the ratio where the error is well above the adjoint noise floor
  true_err <- abs(vapply(sol$u, function(u) u[1], 0) - gea_exact(sol$t))
  ok <- true_err >= 1e-7
  expect_gte(sum(ok), 3L)
  expect_equal(sol$global_error[ok] / true_err[ok], rep(1, sum(ok)), tolerance = 0.05)
})
