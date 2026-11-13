# GlobalAdjoint of GlobalDiffEq 1.7.3 (src/adjoint.jl and the
# SciMLSensitivity extension), against the Julia values of the adjoint
# section of tools/globaldiffeq-probe.jl. One state with one sample, or as
# many samples as states, makes the estimate independent of the random
# directions (Julia draws them from Xoshiro, R from `rng`).
#
# The forward solves, step counts and control decisions match Julia. The
# estimates themselves come from an InterpolatingAdjoint solve whose
# integrand, the defect of the dense output, has a kink at every forward
# node; where an adjoint step straddles one, the embedded error estimate
# misses the kink, so the result moves by about 1e-4 of itself with the
# step placement. Julia's own in place and out of place problems differ by
# that much (ga_lin_est 1.3240616835231346e-4 against ga_lin_est_oop
# 1.3241960614808945e-4), and so does R. The QuadratureAdjoint estimate,
# which has no such steps, matches to 1e-11.

gea_lin <- function(...) ODEProblem(function(u, p, t) p * u, 1, c(0, 2), 2, ...)
gea_lv <- function(u, p, t) c(1.5 * u[1] - u[1] * u[2], -3 * u[2] + u[1] * u[2])

test_that("GlobalAdjoint options and validation", {
  alg <- GlobalAdjoint(Tsit5())
  expect_identical(alg$scope, EndpointError)
  expect_identical(alg$control, ToleranceRefinement)
  expect_identical(GlobalAdjoint(Tsit5(), scope = TrajectoryError)$scope, TrajectoryError)
  expect_identical(GlobalAdjoint(Tsit5(), control = StepGridRefinement)$control, StepGridRefinement)
  expect_output(print(StepGridRefinement), "StepGridRefinement::GlobalAdjointControl = 1", fixed = TRUE)
  expect_identical(GlobalAdjointScope, list(EndpointError = EndpointError, TrajectoryError = TrajectoryError))
  # has_global_error_ga
  expect_true(has_global_error(alg))
  expect_false(allows_arbitrary_number_types(alg))
  expect_false(allowscomplex(alg))
  expect_identical(alg$options$adjoint_abstol, 1e-10)
  expect_identical(GlobalAdjoint(Tsit5(), gtol = 1e-4)$options$adjoint_reltol, 1e-8)
  expect_equal(GlobalAdjoint(Tsit5(), gtol = 1e-9)$options$adjoint_abstol, 1e-11)
  # ga_samples, ga_gtol, ga_nogtol
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

test_that("adjoint_error_estimate matches Julia", {
  # ga_lin_actual: the forward solve
  s <- solve(gea_lin(), Tsit5(), abstol = 1e-4, reltol = 1e-4)
  expect_equal(abs(s$u[[length(s$t)]] - exp(4)), 0.00013244509135290627, tolerance = ref_tol(1e-9))
  # ga_lin_est_quadrature: no adjoint steps over the forward nodes
  expect_equal(adjoint_error_estimate(gea_lin(), GlobalAdjoint(Tsit5(), samples = 1, sensealg = QuadratureAdjoint()),
                                      abstol = 1e-4, reltol = 1e-4), 0.0001324452887558515, tolerance = 1e-9)
  # ga_lin_est (InterpolatingAdjoint, see the note at the top)
  est <- adjoint_error_estimate(gea_lin(), GlobalAdjoint(Tsit5(), samples = 1), abstol = 1e-4, reltol = 1e-4)
  expect_equal(est, 0.00013240616835231346, tolerance = 5e-4)
  # ga_callback: a benign callback is not passed to the adjoint
  cb <- DiscreteCallback(function(u, t, integ) FALSE, function(integ) NULL)
  expect_identical(adjoint_error_estimate(gea_lin(callback = cb), GlobalAdjoint(Tsit5(), samples = 1),
                                          abstol = 1e-4, reltol = 1e-4), est)
  # ga_lv_est: two states, two samples, the norm of the full estimate
  prob <- ODEProblem(gea_lv, c(1, 1), c(0, 2))
  set.seed(5)
  est <- adjoint_error_estimate(prob, GlobalAdjoint(Tsit5()), abstol = 1e-5, reltol = 1e-5)
  expect_equal(est, 0.00011001492535868086, tolerance = 5e-4)
  # ga_forced: a callable parameter
  fprob <- ODEProblem(function(u, p, t) -u + p(t), 0, c(0, 1), function(t) 0.1 * sin(t))
  est <- adjoint_error_estimate(fprob, GlobalAdjoint(Tsit5(), samples = 1), abstol = 1e-4, reltol = 1e-4)
  expect_equal(est, 1.3878314738352899e-8, tolerance = 1e-3)
})

test_that("GlobalAdjoint controls the endpoint error as Julia", {
  # ga_lin_ctl = (25, 147, 54.59815015790747, 1.8049537509688777e-7, 1); the
  # last estimate is near the adjoint's resolution, Julia's and R's differ
  # by 15 percent and both are within a factor 1.5 of the true 1.25e-7
  sol <- solve(gea_lin(), GlobalAdjoint(Tsit5(), gtol = 1e-6, samples = 1), abstol = 1e-3, reltol = 1e-3)
  n <- length(sol$t)
  expect_true(successful_retcode(sol))
  expect_identical(n, 25L)
  expect_equal(sol$stats$nf, 147)
  expect_equal(sol$u[[n]], 54.59815015790747, tolerance = ref_tol(1e-13))
  expect_identical(sol$prob$p, 2)
  expect_length(sol$global_error, n)
  expect_equal(sol$global_error[n], 1.8049537509688777e-7, tolerance = 0.2)
  expect_identical(sum(sol$global_error != 0), 1L)
  # ga_lv_ctl = (31, [6.777188902188865, 2.012753867471386], 5.060473696382924e-7):
  # the refined tolerances follow the estimates, and the adjoint noise of each
  # projection depends on the random direction, so the end state agrees to 1e-8
  set.seed(6)
  sol <- solve(ODEProblem(gea_lv, c(1, 1), c(0, 2)), GlobalAdjoint(Tsit5(), gtol = 1e-6, samples = 2),
               abstol = 1e-3, reltol = 1e-3)
  n <- length(sol$t)
  expect_identical(n, 31L)
  expect_equal(sol$u[[n]], c(6.777188902188865, 2.012753867471386), tolerance = 1e-8)
  expect_equal(sol$global_error[n], 5.060473696382924e-7, tolerance = 0.05)
  # ga_grid = (18, Success, 0.032551292426026744, 54.59814372257128,
  # 6.292812366857436e-6): the fixed step re-solve on the refined grid
  grid <- solve(gea_lin(), GlobalAdjoint(Tsit5(), gtol = 1e-5, control = StepGridRefinement, samples = 1),
                abstol = 1e-3, reltol = 1e-3)
  n <- length(grid$t)
  expect_true(successful_retcode(grid))
  expect_identical(n, 18L)
  expect_identical(grid$t[n], 2)
  expect_equal(grid$t[2], 0.032551292426026744, tolerance = 1e-5)
  expect_equal(grid$u[[n]], 54.59814372257128, tolerance = 1e-11)
  expect_equal(grid$global_error[n], 6.292812366857436e-6, tolerance = 0.02)
  # ga_grid_cheap = (7, Success, 0.0870550563296124, 7.389043134006805, 1.2962803867448026e-5)
  grid <- solve(ODEProblem(function(u, p, t) p * u, 1, c(0, 1), 2),
                GlobalAdjoint(Tsit5(), gtol = 1e-3, control = StepGridRefinement, samples = 1),
                abstol = 1e-3, reltol = 1e-3)
  n <- length(grid$t)
  expect_identical(n, 7L)
  expect_equal(grid$t[2], 0.0870550563296124, tolerance = ref_tol(1e-13))
  expect_equal(grid$u[[n]], 7.389043134006805, tolerance = ref_tol(1e-13))
  expect_equal(grid$global_error[n], 1.2962803867448026e-5, tolerance = 0.01)
})

test_that("TrajectoryError estimates the error at every saved time as Julia", {
  # ga_traj_cheap
  sol <- solve(ODEProblem(function(u, p, t) p * u, 1, c(0, 1), 2),
               GlobalAdjoint(Tsit5(), gtol = 1e-3, scope = TrajectoryError, samples = 1),
               abstol = 1e-3, reltol = 1e-3)
  expect_equal(sol$t, c(0, 0.0870550563296124, 0.2488038859075326, 0.4455147518057155, 0.685192011770674,
                        0.9610195823057235, 1), tolerance = ref_tol(1e-13))
  expect_equal(vapply(sol$u, function(u) u[1], 0),
               c(1, 1.1901866134413397, 1.6447818352387649, 2.4376376334677476, 3.9368600855965195,
                 6.834869724888004, 7.389043134006805), tolerance = ref_tol(1e-13))
  expect_identical(sol$global_error[1], 0)
  expect_equal(sol$global_error[-1],
               c(2.2576753644607202e-10, 3.3595507842574725e-8, 2.8223864568400497e-7, 2.158667720879214e-6,
                 1.1992147262050296e-5, 1.2962803867448026e-5), tolerance = 0.01)
})
