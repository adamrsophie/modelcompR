# Reference values from tools/scimlsensitivity-probe.jl (SciMLSensitivity
# 7.119.14).

test_that("supports_functor_params, SundialsAdjoint, SensitivityAlg and ReactantDualTag", {
  expect_true(supports_functor_params(GaussAdjoint()))
  expect_true(supports_functor_params(GaussKronrodAdjoint()))
  for (sa in list(QuadratureAdjoint(), InterpolatingAdjoint(), BacksolveAdjoint(), ForwardSensitivity(),
                  ForwardDiffSensitivity(), SteadyStateAdjoint()))
    expect_false(supports_functor_params(sa))
  expect_error(supports_functor_params(1), "sensitivity algorithm")
  sa <- SundialsAdjoint()
  expect_identical(sa$steps, 150)
  expect_identical(sa$interp, "hermite")
  expect_true(sa$quad_error_control)
  expect_error(SundialsAdjoint(interp = "linear"),
               "SundialsAdjoint `interp` must be `:hermite` or `:polynomial`, got `linear`.", fixed = TRUE)
  expect_message(SensitivityAlg(), "The SensitivityAlg choice mechanism was completely overhauled")
  expect_s3_class(ReactantDualTag(), "ReactantDualTag")
  # SundialsAdjoint runs the checkpointed interpolating adjoint
  prob <- ODEProblem(function(u, p, t) c(-p[1] * u[1] + u[2], -p[2] * u[2]), c(1, 0.5), c(0, 1), c(1.3, 0.4))
  sol <- solve(prob, Tsit5(), abstol = 1e-10, reltol = 1e-10)
  dg <- function(u, p, t, i) c(1, 2)
  a <- adjoint_sensitivities(sol, Tsit5(), t = 1, dgdu_discrete = dg, sensealg = SundialsAdjoint())
  b <- adjoint_sensitivities(sol, Tsit5(), t = 1, dgdu_discrete = dg, sensealg = InterpolatingAdjoint(checkpointing = TRUE))
  expect_identical(a, b)
})

lv_sens <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1] * u[2], -p[3] * u[2] + p[4] * u[1] * u[2] + t * p[1])
p_sens <- c(1.5, 1.0, 3.0, 1.0)

test_that("ODEForwardSensitivityFunction is the forward sensitivity right-hand side", {
  z <- c(1.2, 0.7, 0.1, -0.2, 0.3, 0.4, -0.5, 0.6, 0.7, -0.8)
  jl <- c(0.9599999999999999, -0.8099999999999998, 1.52, 0.7300000000000001, -1.08, -0.5100000000000002, -1.12,
          -2.13, 1.5199999999999998, 2.7700000000000005)
  fs <- ODEForwardSensitivityFunction(ODEProblem(lv_sens, c(1, 1), c(0, 2), p_sens))
  expect_s3_class(fs, "ODEForwardSensitivityFunction")
  expect_false(inherits(fs, SensitivityFunction))
  expect_identical(attr(fs, "numindvar"), 2L)
  expect_identical(attr(fs, "numparams"), 4L)
  expect_equal(fs(z, p_sens, 0.3), jl, tolerance = 1e-9)
  # with the symbolic Jacobians of a System (unknowns y, x; parameters a, b, d, c)
  tt <- t_nounits; D <- D_nounits
  variables("x_fs(t)", "y_fs(t)"); parameters("a_fs", "b_fs", "c_fs", "d_fs")
  sys <- mtkcompile(System(list(D(x_fs) ~ a_fs * x_fs - b_fs * x_fs * y_fs,
                                D(y_fs) ~ -c_fs * y_fs + d_fs * x_fs * y_fs + tt * a_fs), tt, name = "lvfs"))
  sprob <- ODEProblem(sys, list(x_fs ~ 1, y_fs ~ 1, a_fs ~ 1.5, b_fs ~ 1, c_fs ~ 3, d_fs ~ 1), c(0, 2))
  ux <- vapply(unknowns(sys), function(v) .key(v) == .key(x_fs), TRUE)
  iu <- if (ux[1L]) c(1, 2) else c(2, 1)
  pk <- vapply(.problem_params(sys), .key, "")
  ip <- match(pk, c(.key(a_fs), .key(b_fs), .key(c_fs), .key(d_fs)))
  zS <- matrix(z[3:10], 2, 4)
  jS <- matrix(jl[3:10], 2, 4)
  out <- ODEForwardSensitivityFunction(sprob)(c(z[iu], zS[iu, ip]), sprob$p, 0.3)
  expect_equal(out, c(jl[iu], jS[iu, ip]), tolerance = 1e-14)
})

test_that("SensitivityFunction tags the adjoint right-hand sides", {
  sol <- solve(ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2), Tsit5())
  ap <- ODEAdjointProblem(sol, InterpolatingAdjoint(), t = 1, dgdu_discrete = function(u, p, t, i) 1)
  expect_true(inherits(ap$f$f, SensitivityFunction))
  expect_s3_class(ap$f$f, "ODEInterpolatingAdjointSensitivityFunction")
  ap <- ODEAdjointProblem(sol, QuadratureAdjoint(), t = 1, dgdu_discrete = function(u, p, t, i) 1)
  expect_s3_class(ap$f$f, "ODEQuadratureAdjointSensitivityFunction")
})

test_that("AdjointSensitivityIntegrand is lambda' f_p + dgdp along the adjoint solution", {
  sol <- solve(ODEProblem(lv_sens, c(1, 1), c(0, 2), p_sens), Tsit5(), abstol = 1e-12, reltol = 1e-12)
  ap <- ODEAdjointProblem(sol, QuadratureAdjoint(), Tsit5(), t = 2, dgdu_discrete = function(u, p, t, i) c(1, 1),
                          dgdu_continuous = function(u, p, t) c(u[1], 0),
                          dgdp_continuous = function(u, p, t) numeric(4))
  adj <- solve(ap, Tsit5(), abstol = 1e-12, reltol = 1e-12)
  ig <- AdjointSensitivityIntegrand(sol, adj)
  expect_s3_class(ig, "AdjointSensitivityIntegrand")
  expect_equal(ig(0), c(4.388385806326229, -4.388385806326229, 3.209238117606435, -3.209238117606435),
               tolerance = 1e-9)
  expect_equal(ig(0.5), c(3.1238049083504036, -2.858171788923313, 2.259427363736097, -3.373393773556924),
               tolerance = 1e-9)
  expect_equal(ig(1.3), c(2.9815255136694954, -6.936503197278934, 2.6574441960460775, -7.2409405730025735),
               tolerance = 1e-9)
  expect_equal(.sens_quadgk(ig, 0, 2, 1e-8, 1e-8),
               c(6.758264897027363, -10.745074463633323, 3.7275618030609388, -7.191688506408433), tolerance = 1e-9)
  ig2 <- AdjointSensitivityIntegrand(sol, adj, dgdp = function(u, p, t) c(u[1] * u[2], 0, 0, t))
  expect_equal(ig2(0.5), c(3.949925808933976, -2.858171788923313, 2.259427363736097, -2.873393773556924),
               tolerance = 1e-9)
})

test_that("StochasticTransformedFunction subtracts the Ito correction", {
  fs <- function(u, p, t) c(p[1] * u[1] - u[2], p[2] * u[2] + u[1] * t)
  gs <- function(u, p, t) c(p[3] * u[1] * u[2], p[4] * sin(u[2]))
  pp <- c(1.01, -0.5, 0.3, 0.7)
  sp <- SDEProblem(fs, gs, c(1, 0.5), c(0, 1), pp)
  ssol <- solve(sp, EM(), dt = 0.1, seed = 1)
  tf <- StochasticTransformedFunction(ssol)
  expect_s3_class(tf, "StochasticTransformedFunction")
  expect_equal(tf(c(1.3, -0.4), pp, 0.25), c(1.69428, 0.7615922422703831), tolerance = 1e-9)
  gm <- function(u, p, t) matrix(c(p[3] * u[1], p[4] * u[1] * u[2], u[2], cos(u[1])), 2, 2)
  mp <- SDEProblem(fs, gm, c(1, 0.5), c(0, 1), pp, noise_rate_prototype = matrix(0, 2, 2))
  tm <- StochasticTransformedFunction(solve(mp, EM(), dt = 0.1, seed = 1))
  expect_equal(tm(c(1.3, -0.4), pp, 0.25), c(1.7518306859107322, 1.25624), tolerance = 1e-9)
  tc <- StochasticTransformedFunction(ssol, corfunc_analytical = function(u, p, t) c(u[1], 2 * u[2]))
  expect_equal(tc(c(1.3, -0.4), pp, 0.25), c(0.41300000000000003, 1.3250000000000002), tolerance = 1e-15)
})

test_that("UnconstrainedOptimizationAdjoint and OptimizationAdjoint differentiate the optimum", {
  fu <- function(x, p) (x[1] - p[1])^2 + p[2] * (x[2] - x[1]^2)^2 + p[3] * x[1]
  gu <- function(x, p) c(2 * (x[1] - p[1]) - 4 * p[2] * x[1] * (x[2] - x[1]^2) + p[3], 2 * p[2] * (x[2] - x[1]^2))
  hu <- function(x, p) rbind(c(2 - 4 * p[2] * (x[2] - 3 * x[1]^2), -4 * p[2] * x[1]), c(-4 * p[2] * x[1], 2 * p[2]))
  pu <- c(1, 3, 0.5)
  probu <- OptimizationProblem(OptimizationFunction(fu, grad = gu, hess = hu), c(0, 0), pu)
  jl <- c(3.999999999999982, 2.960594732333756e-16, -1.999999999999991)
  # at the optimum of Julia's BFGS
  su <- .opt_solution(probu, NULL, c(0.7499999999999978, 0.5624999999999962), 0, "Success", NULL)
  dgdu <- function(u, p) c(1, 2)
  for (sa in list(UnconstrainedOptimizationAdjoint(), OptimizationAdjoint()))
    expect_equal(adjoint_sensitivities(su, sensealg = sa, dgdu = dgdu)$dp, jl, tolerance = 1e-9)
  expect_error(adjoint_sensitivities(su, sensealg = OptimizationAdjoint()), "dgdu must be specified")
  # through solve_gradient; the default is the unconstrained adjoint
  expect_s3_class(.sens_default_alg(probu), "UnconstrainedOptimizationAdjoint")
  # (to the accuracy of the optimum BFGS stops at)
  g <- solve_gradient(function(sol) sum(sol$u * c(1, 2)), probu, BFGS(), g_abstol = 1e-14)
  expect_equal(g$p, jl, tolerance = 1e-6)
  # an equality, an active and an inactive inequality and an active bound
  fc <- function(x, p) (x[1] - p[1])^2 + (x[2] - p[2])^2 + p[4] * (x[3] - p[3])^2 + p[6] * x[1] * x[2]
  cons <- function(x, p) c(x[1] + 2 * x[2] + x[3] - p[5], x[1]^2 + x[2]^2 - 1, x[3] - x[1] - 10)
  pc <- c(2, 1.5, -1, 1, 1.5, 0.3)
  probc <- OptimizationProblem(fc, c(0.1, 0.1, 0.5), pc, cons = cons, lcons = c(0, -Inf, -Inf), ucons = c(0, 0, 0),
                               lb = c(-Inf, -Inf, 0), ub = c(Inf, Inf, Inf))
  expect_s3_class(.sens_default_alg(probc), "OptimizationAdjoint")
  xj <- c(0.9633249580710704, 0.2683375209644633, 3.0837022369139185e-15)
  sc <- .opt_solution(probc, NULL, xj, fc(xj, pc), "Success", NULL)
  expect_equal(adjoint_sensitivities(sc, dgdu = function(u, p) c(1, 2, -1))$dp, c(0, 0, 0, 0, 1, 0),
               tolerance = 1e-9)
  # x1 + 2 x2 = p5 and x1^2 + x2^2 = 1 are active: dx1/dp5 = x2 / (x2 - 2 x1)
  expect_equal(adjoint_sensitivities(sc, dgdu = function(u, p) c(1, 0, 0))$dp,
               c(0, 0, 0, 0, xj[2] / (xj[2] - 2 * xj[1]), 0), tolerance = 1e-9)
})
