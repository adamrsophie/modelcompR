# GlobalErrorEstimation of GlobalDiffEq 1.7.3 (src/companion.jl,
# src/estimation.jl). The Julia values of the same problems are printed by
# the companion section of tools/globaldiffeq-probe.jl; here the estimates
# are checked against the true global error of a tight reference solve, as
# in GlobalDiffEq's test/companion_tests.jl, and the dense derivatives
# against finite differences of the dense output.

gec_lv <- function(u, p, t) c(1.5 * u[1] - u[1] * u[2], -3 * u[2] + u[1] * u[2])
gec_jac <- function(u, p, t) matrix(c(1.5 - u[2], u[2], -u[1], -3 + u[1]), 2, 2)
gec_prob <- function(...) ODEProblem(ODEFunction(gec_lv, jac_fn = gec_jac), c(1, 1), c(0, 2), ...)
gec_ref <- local({
  s <- solve(gec_prob(), Vern9(), abstol = 1e-13, reltol = 1e-13)
  s$u[[length(s$t)]]
})
gec_err <- function(sol) sqrt(sum((sol$u[[length(sol$t)]] - gec_ref)^2))

test_that("the options print as Julia enums", {
  expect_output(print(DefectCorrection), "DefectCorrection::GlobalErrorEquation = 0", fixed = TRUE)
  expect_output(print(SimultaneousMode), "SimultaneousMode::GlobalErrorMode = 1", fixed = TRUE)
  expect_identical(GlobalErrorEquation, list(DefectCorrection = DefectCorrection, ErrorTransport = ErrorTransport))
  expect_identical(GlobalErrorMode, list(InterpolatingMode = InterpolatingMode, SimultaneousMode = SimultaneousMode))
  alg <- GlobalErrorEstimation(Tsit5())
  expect_identical(alg$equation, DefectCorrection)
  expect_identical(alg$mode, InterpolatingMode)
  expect_identical(GlobalErrorEstimation(Tsit5(), equation = ErrorTransport)$equation, ErrorTransport)
  expect_false(has_global_error(alg))
  expect_false(allows_arbitrary_number_types(alg))
  expect_false(allowscomplex(alg))
})

test_that("the dense first derivatives are the derivatives of the dense output", {
  for (alg in list(Tsit5(), DP5(), Vern7(), BS3(), Rosenbrock23(), Rodas5P(), Trapezoid(), FBDF())) {
    sol <- solve(gec_prob(), alg, abstol = 1e-6, reltol = 1e-6)
    expect_true(.gde_has_dense_derivative(sol))
    dn <- .gde_sol_dense(sol, right = FALSE)
    j <- ceiling(length(sol$t) / 2)
    t <- (sol$t[j] + sol$t[j + 1]) / 2
    h <- 1e-3 * (sol$t[j + 1] - sol$t[j])
    fd <- (8 * (dn$value(t + h) - dn$value(t - h)) - (dn$value(t + 2 * h) - dn$value(t - 2 * h))) / (12 * h)
    expect_equal(dn$deriv(t), fd, tolerance = 1e-7, info = alg$name)
    # right continuity at a node takes the next step; at the end point the last
    dr <- .gde_sol_dense(sol, right = TRUE)
    expect_equal(dr$deriv(2), dn$deriv(2), info = alg$name)
    expect_equal(dr$value(sol$t[j]), sol$u[[j]], info = alg$name)
  }
  # the linear dense output of the low storage and SSP methods has none (an
  # R limitation: OrdinaryDiffEq gives them Hermite or quadratic output)
  expect_false(.gde_alg_has_dense_derivative(SSPRK33()))
})

test_that("the companion estimators track the true global error", {
  for (eq in GlobalErrorEquation) for (md in GlobalErrorMode) {
    est <- global_error_estimate(gec_prob(), GlobalErrorEstimation(Tsit5(), equation = eq, mode = md),
                                 abstol = 1e-4, reltol = 1e-4)
    true_err <- gec_err(solve(gec_prob(), Tsit5(), abstol = 1e-4, reltol = 1e-4))
    expect_equal(est / true_err, 1, tolerance = 0.01, info = paste(eq, md))
  }
  # the streaming mode restarts the companion at every step
  a <- global_error_estimate(gec_prob(), GlobalErrorEstimation(Tsit5()), abstol = 1e-6, reltol = 1e-6)
  b <- global_error_estimate(gec_prob(), GlobalErrorEstimation(Tsit5(), mode = SimultaneousMode),
                             abstol = 1e-6, reltol = 1e-6)
  expect_equal(b, a, tolerance = 1e-5)
  # ErrorTransport without a Jacobian: a directional difference
  c1 <- global_error_estimate(ODEProblem(gec_lv, c(1, 1), c(0, 2)), GlobalErrorEstimation(Tsit5(),
                              equation = ErrorTransport), abstol = 1e-6, reltol = 1e-6)
  c2 <- global_error_estimate(gec_prob(), GlobalErrorEstimation(Tsit5(), equation = ErrorTransport),
                              abstol = 1e-6, reltol = 1e-6)
  expect_equal(c1, c2, tolerance = 1e-5)
  # a scalar linear problem: the error of y' = 2 y at t = 1 is known
  lin <- ODEProblem(function(u, p, t) p * u, 1, c(0, 1), 2)
  est <- global_error_estimate(lin, GlobalErrorEstimation(Tsit5()), abstol = 1e-5, reltol = 1e-5)
  s <- solve(lin, Tsit5(), abstol = 1e-5, reltol = 1e-5)
  expect_equal(est, abs(s$u[[length(s$t)]] - exp(2)), tolerance = 1e-3)
})

test_that("other dense solvers and companion solvers", {
  for (alg in list(BS3(), Rodas5P())) {
    est <- global_error_estimate(gec_prob(), GlobalErrorEstimation(alg), abstol = 1e-5, reltol = 1e-5)
    expect_equal(est / gec_err(solve(gec_prob(), alg, abstol = 1e-5, reltol = 1e-5)), 1, tolerance = 0.05,
                 info = alg$name)
  }
  # an implicit Euler companion only samples the defect at the forward
  # nodes, where it vanishes; an explicit companion resolves it
  est <- global_error_estimate(gec_prob(), GlobalErrorEstimation(ImplicitEuler(), companion_alg = Tsit5()),
                               abstol = 1e-3, reltol = 1e-3)
  expect_equal(est / gec_err(solve(gec_prob(), ImplicitEuler(), abstol = 1e-3, reltol = 1e-3)), 1,
               tolerance = 0.01)
})

test_that("global error control meets gtol", {
  for (md in GlobalErrorMode) {
    sol <- solve(gec_prob(), GlobalErrorEstimation(Tsit5(), mode = md, gtol = 1e-7), abstol = 1e-3, reltol = 1e-3)
    expect_true(successful_retcode(sol))
    expect_lte(gec_err(sol), 1e-7)
  }
  expect_error(solve(gec_prob(), GlobalErrorEstimation(Tsit5(), gtol = 1e-9, maxiters = 1), abstol = 1e-3,
                     reltol = 1e-3), "failed to meet gtol=1e-09 after 1 iterations", fixed = TRUE)
})

test_that("GlobalErrorEstimation validates its arguments as Julia", {
  alg <- GlobalErrorEstimation(Tsit5())
  expect_error(GlobalErrorEstimation(Tsit5(), gtol = -1), "gtol must be a positive finite real number")
  expect_error(GlobalErrorEstimation(Tsit5(), maxiters = 0), "maxiters must be a positive integer")
  expect_error(GlobalErrorEstimation(Tsit5(), safety = 2), "safety must be between zero and one")
  expect_error(GlobalErrorEstimation(Tsit5(), companion_abstol = 0, companion_reltol = 0),
               "companion_abstol and companion_reltol cannot both be zero")
  expect_error(solve(gec_prob(), alg), "GlobalErrorEstimation requires a positive `gtol` constructor keyword",
               fixed = TRUE)
  cb <- DiscreteCallback(function(u, t, integ) FALSE, function(integ) NULL)
  expect_error(global_error_estimate(gec_prob(callback = cb), alg),
               "GlobalErrorEstimation does not currently support callbacks")
  expect_error(global_error_estimate(ODEProblem(gec_lv, c(1, 1), c(2, 0)), alg),
               "GlobalErrorEstimation requires a forward-time ODEProblem")
  mprob <- ODEProblem(ODEFunction(gec_lv, mass_matrix = diag(c(1, 2))), c(1, 1), c(0, 2))
  expect_error(global_error_estimate(mprob, GlobalErrorEstimation(Rodas5P())),
               "GlobalErrorEstimation currently requires the standard mass matrix")
  expect_error(solve(gec_prob(), GlobalErrorEstimation(Tsit5(), gtol = 1e-6), abstol = -1),
               "local_abstol must be a nonnegative finite real number")
  # Julia's message for a solution without a derivative interpolation
  expect_error(global_error_estimate(gec_prob(), GlobalErrorEstimation(SSPRK33()), dt = 0.1),
               "requires a solver with dense (continuous) output that provides a first-derivative", fixed = TRUE)
  expect_error(global_error_estimate(gec_prob(), alg, dense = FALSE), NA)
})
