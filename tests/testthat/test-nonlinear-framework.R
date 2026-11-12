nl_prob <- function() NonlinearProblem(function(u, p) c(u[1]^2 + u[2]^2 - 2, u[1] - u[2]), c(2, 0.5))

test_that("composed first order methods converge", {
  prob <- nl_prob()
  algs <- list(GeneralizedFirstOrderAlgorithm(), GeneralizedFirstOrderAlgorithm(descent = DampedNewtonDescent()),
               GeneralizedFirstOrderAlgorithm(descent = SteepestDescent(), linesearch = BackTracking()),
               GeneralizedFirstOrderAlgorithm(descent = GeodesicAcceleration()),
               TrustRegionDogleg(), TrustRegionRobust())
  for (alg in algs) {
    sol <- solve(prob, alg)
    expect_identical(sol$retcode, "Success", info = alg$name)
    expect_equal(as.numeric(sol$u), c(1, 1), tolerance = 1e-8, info = alg$name)
  }
  for (m in unlist(RadiusUpdateSchemes)) for (d in list(Dogleg(), MoreTrustRegionDescent())) {
    alg <- GeneralizedFirstOrderAlgorithm(descent = d, trustregion = GenericTrustRegionScheme(m))
    expect_equal(as.numeric(solve(prob, alg)$u), c(1, 1), tolerance = 1e-8, info = paste(m, d$name))
  }
})

test_that("line searches", {
  prob <- nl_prob()
  for (ls in list(BackTracking(), BackTracking(order = 2), ArmijoLineSearch(), GoldenSection(),
                  LiFukushimaLineSearch(), StrongWolfeLineSearch(), MoreThuenteLineSearch(),
                  HagerZhangLineSearch(), NoLineSearch(), LineSearchesJL("HagerZhang"), ProjectedBackTracking())) {
    sol <- solve(prob, GeneralizedFirstOrderAlgorithm(linesearch = ls))
    expect_equal(as.numeric(sol$u), c(1, 1), tolerance = 1e-8, info = ls$name)
  }
  # on a quadratic merit the Wolfe searches find the exact minimizer along the direction
  phi <- function(a) (a - 0.3)^2
  dphi <- function(a) 2 * (a - 0.3)
  for (ls in list(StrongWolfeLineSearch(), MoreThuenteLineSearch(), HagerZhangLineSearch())) {
    r <- .ls_run(ls, phi, dphi(0), dphi)
    expect_lt(abs(dphi(r$step_size)), 0.9 * abs(dphi(0)))
  }
  expect_s3_class(LineSearchSolution(0.5), "LineSearchSolution")
  expect_error(set_initial_step(list(), -1), "positive")
})

test_that("quasi Newton, spectral and simple methods", {
  prob <- nl_prob()
  for (alg in list(QuasiNewtonAlgorithm(), QuasiNewtonAlgorithm(update_rule = "bad_broyden"),
                   QuasiNewtonAlgorithm(update_rule = "klement", initialization = "true_jacobian"),
                   LimitedMemoryBroyden(), DFSane(), SimpleHalley(), SimpleBroyden(), SimpleKlement(),
                   SimpleDFSane(), SimpleLimitedMemoryBroyden(), SimpleGaussNewton())) {
    sol <- solve(prob, alg)
    expect_identical(sol$retcode, "Success", info = alg$name)
    expect_equal(as.numeric(sol$u), c(1, 1), tolerance = 1e-6, info = alg$name)
  }
  cube <- NonlinearProblem(function(u, p) u^3 - 8, 0.5)
  expect_equal(as.numeric(solve(cube, DFSane())$u), 2, tolerance = 1e-8)
  expect_equal(as.numeric(solve(cube, PseudoTransient())$u), 2, tolerance = 1e-8)
})

test_that("termination modes", {
  prob <- nl_prob()
  for (tm in list(AbsTerminationMode(), RelTerminationMode(), NormTerminationMode(), AbsNormTerminationMode(),
                  RelNormTerminationMode(), AbsNormSafeTerminationMode(), RelNormSafeTerminationMode(),
                  AbsNormSafeBestTerminationMode(), RelNormSafeBestTerminationMode())) {
    sol <- solve(prob, GeneralizedFirstOrderAlgorithm(), termination_condition = tm)
    expect_identical(sol$retcode, "Success", info = class(tm)[1])
  }
  # a loose absolute tolerance stops early
  sol <- solve(prob, GeneralizedFirstOrderAlgorithm(), termination_condition = AbsTerminationMode(), abstol = 0.1)
  expect_lte(max(abs(prob$f$f(sol$u, prob$p))), 0.1)
})

test_that("more bracketing algorithms", {
  prob <- IntervalNonlinearProblem(function(u, p) u^3 - 2, c(0, 2))
  for (alg in list(Alefeld(), ModAB(), Ridder(), Muller()))
    expect_equal(solve(prob, alg, abstol = 1e-12)$u, 2^(1 / 3), tolerance = 1e-9, info = alg$name)
})
