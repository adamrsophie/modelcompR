# Reference values from tools/scimlsensitivity-probe.jl (SciMLSensitivity
# 7.119.14).

dae1_res <- function(du, u, p, t) {
  c(-p[1] * u[1] + p[2] * u[1] * u[2] - du[1], p[3] * u[1] - u[2] * u[3] - du[2], u[3] - u[1]^2 - p[4] * u[2])
}
dae1_prob <- function() {
  DAEProblem(dae1_res, c(-0.5, 0.6, 0), c(1, 0.5, 1.2), c(0, 1.5), c(0.8, 0.6, 1.2, 0.4),
             differential_vars = c(TRUE, TRUE, FALSE))
}

test_that("DAEAdjointProblem differentiates an index-1 DAE", {
  sol <- solve(dae1_prob(), DFBDF(), abstol = 1e-10, reltol = 1e-10)
  ts <- c(0.5, 1, 1.5)
  dg <- function(u, p, t, i) u
  # ForwardDiff through DFBDF (u3 follows from the constraint)
  fd_p <- c(-3.5051415342679215, 2.437199190021155, 2.0929343489243286, 0.5471313405535572)
  fd_u <- c(4.105844498430743, 2.2900303669193605, 0)
  r <- adjoint_sensitivities(sol, FBDF(), t = ts, dgdu_discrete = dg, abstol = 1e-10, reltol = 1e-10)
  expect_equal(unname(r$du0), fd_u, tolerance = 1e-8)
  expect_equal(unname(r$dp), fd_p, tolerance = 1e-8)
  expect_identical(r$du0[[3]], 0)
  r <- adjoint_sensitivities(sol, FBDF(), t = ts, dgdu_discrete = dg, sensealg = QuadratureAdjoint(),
                             abstol = 1e-10, reltol = 1e-10)
  expect_equal(unname(r$dp), fd_p, tolerance = 1e-8)
  # the adjoint problem itself: a singular mass matrix ODE in z = (w, lambda, grad)
  ap <- DAEAdjointProblem(sol, InterpolatingAdjoint(), FBDF(), t = ts, dgdu_discrete = dg)
  expect_s3_class(ap, "ODEProblem")
  expect_identical(diag(ap$f$mass_matrix), c(1, 1, 0, 0, 0, 0, 1, 1, 1, 1))
  expect_s3_class(ap$f$f, "DAEAdjointSensitivityFunction")
  expect_true(inherits(ap$f$f, SensitivityFunction))
  expect_identical(ap$tspan, c(1.5, 0))
  # with a DAE algorithm it is a DAEProblem
  ap <- DAEAdjointProblem(sol, InterpolatingAdjoint(), DFBDF(), t = ts, dgdu_discrete = dg)
  expect_s3_class(ap, "DAEProblem")
  expect_s3_class(ap$f$f, "DAEAdjointResidual")
})

test_that("DAEAdjointProblem integrates a continuous cost", {
  sol <- solve(dae1_prob(), DFBDF(), abstol = 1e-10, reltol = 1e-10)
  r <- adjoint_sensitivities(sol, FBDF(), dgdu_continuous = function(u, p, t) c(u[1], 0, 1),
                             abstol = 1e-10, reltol = 1e-10)
  expect_equal(unname(r$du0), c(2.683085494978922, 0.8695880121398295, 0), tolerance = 1e-8)
  expect_equal(unname(r$dp), c(-1.6536657446011849, 1.120624713925577, 0.48161718547632143, 0.8508994893264923),
               tolerance = 1e-8)
})

test_that("DAEAdjointProblem handles a Hessenberg index-2 DAE", {
  hess2 <- function(du, u, p, t) c(p[1] * u[2] + u[3] - du[1], -p[2] * u[1] - du[2], u[1] - sin(p[3] * t))
  p <- c(1, 2, 3)
  prob <- DAEProblem(hess2, c(3, 0, 0), c(0, 1, 2), c(0, 1), p, differential_vars = c(TRUE, TRUE, FALSE))
  sol <- solve(prob, DImplicitEuler(), adaptive = FALSE, dt = 1e-3, initializealg = NoInit())
  ts <- seq(0.2, 1, by = 0.2)
  dg <- function(u, p, t, i) c(u[1], u[2], 0)
  r <- adjoint_sensitivities(sol, FBDF(), t = ts, dgdu_discrete = dg, abstol = 1e-8, reltol = 1e-8)
  expect_equal(unname(r$du0), c(0, 1.152119134707172, 0), tolerance = 1e-7)
  expect_equal(unname(r$dp), c(0, 0.06142028142309, -0.6751766620394185), tolerance = 1e-4)
  # the analytic gradient of the exact solution, to the O(dt) forward error
  y2 <- 1 + p[2] * (cos(p[3] * ts) - 1) / p[3]
  dy2_dp3 <- p[2] * (-ts * sin(p[3] * ts) * p[3] - (cos(p[3] * ts) - 1)) / p[3]^2
  expect_equal(r$dp[[3]], sum(sin(p[3] * ts) * ts * cos(p[3] * ts) + y2 * dy2_dp3), tolerance = 1e-2)
  expect_error(adjoint_sensitivities(sol, FBDF(), t = ts, dgdu_discrete = dg, sensealg = QuadratureAdjoint()),
               "`QuadratureAdjoint` does not support Hessenberg index-2 DAEs", fixed = TRUE)
  expect_error(adjoint_sensitivities(sol, FBDF(), t = ts, dgdu_discrete = function(u, p, t, i) c(0, 0, 1)),
               "Discrete cost contributions on the algebraic variables of an index-2 DAE")
})

test_that("DAEAdjointProblem rejects the unsupported configurations", {
  f <- function(du, u, p, t) c(du[1] + p[1] * u[1] - u[2], u[2] - p[2] * u[1])
  prob <- DAEProblem(f, c(-0.5, 0), c(1, 0.5), c(0, 1), c(1, 0.5), differential_vars = c(TRUE, FALSE))
  sol <- solve(prob, DFBDF(), abstol = 1e-8, reltol = 1e-8)
  dg <- function(u, p, t, i) c(1, 0)
  expect_error(adjoint_sensitivities(sol, DFBDF(), t = 1, dgdu_discrete = dg, sensealg = BacksolveAdjoint()),
               "`BacksolveAdjoint` is not currently supported for adjoint sensitivity analysis of fully implicit",
               fixed = TRUE)
  expect_error(adjoint_sensitivities(sol, DFBDF(), t = 1, dgdu_discrete = dg,
                                     sensealg = InterpolatingAdjoint(checkpointing = TRUE)),
               "Checkpointing is not currently supported for DAEProblem adjoints. Use a dense forward solution with `InterpolatingAdjoint(checkpointing = false)`.",
               fixed = TRUE)
  expect_error(DAEAdjointProblem(sol, InterpolatingAdjoint(), t = 1),
               "Either `dgdu_discrete`, `dgdu_continuous`, or `g` must be specified.", fixed = TRUE)
  # G = u1(1) = exp(-(p1 - p2)) for the documented example
  r <- adjoint_sensitivities(sol, DFBDF(), t = 1, dgdu_discrete = dg, abstol = 1e-10, reltol = 1e-10)
  expect_equal(unname(r$dp), c(-1, 1) * exp(-0.5), tolerance = 1e-6)
  expect_equal(unname(r$du0), c(exp(-0.5), 0), tolerance = 1e-6)
})
