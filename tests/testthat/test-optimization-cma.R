# CMA-ES optimizers: CMAEvolutionStrategyOpt (CMAEvolutionStrategy.jl),
# GCMAESOpt (GCMAES.jl) and PyCMAOpt (pycma). The strategy parameters and
# the bound transformations are deterministic and checked against the
# reference implementations (numbers from tools/cma-probe.jl and
# tools/pycma-probe.py); the runs are random and checked as in the
# Optimization.jl tests, with fixed R seeds.

cma_rosen <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
cma_rgrad <- function(x, p) c(-2 * (p[1] - x[1]) - 4 * p[2] * (x[2] - x[1]^2) * x[1], 2 * p[2] * (x[2] - x[1]^2))

test_that("CMAEvolutionStrategy parameters and bound transformation match CMAEvolutionStrategy.jl", {
  # tools/cma-probe.jl
  ref <- list(
    list(n = 2, lam = 6,
         w = c(0.6370425712412168, 0.28457025743803294, 0.07838717132075033, -0.28638378259655295,
               -0.7649580940851275, -1.1559817781589212), mueff = 2.0286114646100617,
         sigma = c(0.5731731629916698, 1.5731731629916696, 1.2880359845663278, 1.254272742818995),
         cov = c(0.6245545390268264, 0.1548153998964136, 0.08559277942666424, 1.3200984146732668)),
    list(n = 10, lam = 10,
         w = c(0.45627264690340597, 0.2707530970017852, 0.16223111715866978, 0.08523354710016448,
               0.025509591835974777, -0.08001260758087221, -0.22176416099914648, -0.3445549417848209,
               -0.45286408637842174, -0.5497499176973852), mueff = 3.1672992814107035,
         sigma = c(0.3196142529106335, 1.3196142529106334, 1.304253812100349, 3.0847265651690123),
         cov = c(0.29499038303562225, 0.015283824524751714, 0.02355177665041749, 1.2621526878513811)),
    list(n = 2, lam = 20,
         w = c(0.27961472096101825, 0.19718884917099674, 0.14897280508458438, 0.11476297738097527,
               0.08822777350618366, 0.06654693329456288, 0.048216044086412384, 0.03233710559095377,
               0.018330889208150472, 0.005801901716162117, -0.007623169971597677, -0.021881602760196954,
               -0.03499808944350321, -0.047142059308261965, -0.05844782477277602, -0.06902366206845893,
               -0.07895813415027543, -0.08832460357342892, -0.09718452358872108, -0.10558988408103799),
         mueff = 5.938804235601242,
         sigma = c(0.7257469888494528, 2.008816775935778, 2.3435269942493324, 1.254272742818995),
         cov = c(0.5837604822280295, 0.11884385675893783, 0.39721292501393574, 2.2158237094595674)))
  for (r in ref) {
    pr <- modelcompR:::.cmaes_parameters(r$n, r$lam)
    expect_equal(pr$w, r$w, tolerance = 1e-15)
    expect_equal(pr$mueff, r$mueff, tolerance = 1e-15)
    expect_equal(unlist(pr$sigma[c("c", "d", "e", "chin")], use.names = FALSE), r$sigma, tolerance = 1e-15)
    expect_equal(unlist(pr$cov[c("c", "c1", "cmu", "e")], use.names = FALSE), r$cov, tolerance = 1e-15)
  }
  b <- modelcompR:::.cmaes_box(c(-1, -1), c(0.8, 0.8))
  tf <- function(x) modelcompR:::.cmaes_transform(b, x)
  expect_equal(tf(c(-1.03, 0.79)), c(-0.98775, 0.7722222222222223), tolerance = 1e-15)
  expect_equal(tf(c(0.83, -1.2)), c(0.79, -0.9750000000000001), tolerance = 1e-15)
  expect_equal(tf(c(5.0, -7.3)), c(0.7530555555555559, 0.660000000000001), tolerance = 1e-15)
  expect_equal(tf(c(0.76, -0.97)), c(0.7530555555555556, -0.9577499999999999), tolerance = 1e-15)
  expect_equal(modelcompR:::.cmaes_backtransform(b, c(-0.99, 0.79)), c(-1.0367544467966325, 0.83), tolerance = 1e-15)
  expect_equal(modelcompR:::.cmaes_backtransform(b, c(0.8, -1)), c(0.89, -1.1), tolerance = 1e-15)
  b2 <- modelcompR:::.cmaes_box(c(-10, 2), c(Inf, 3))
  expect_equal(modelcompR:::.cmaes_transform(b2, c(-12, 3.4)), c(-9.100000000000001, 2.9499999999999997),
               tolerance = 1e-15)
})

test_that("GCMAES and pycma parameters match GCMAES.jl and pycma", {
  # tools/cma-probe.jl
  g2 <- modelcompR:::.gcmaes_parameters(2)
  expect_equal(c(g2$lam, g2$mu), c(6, 3))
  expect_equal(g2$w, c(0.6370425712412168, 0.28457025743803294, 0.07838717132075033), tolerance = 1e-15)
  expect_equal(unlist(g2[c("mueff", "cc", "cs", "c1", "cmu", "ds", "chin")], use.names = FALSE),
               c(2.0286114646100617, 0.6666666666666666, 0.5731731629916698, 0.1548153998964136,
                 0.057859085071916304, 1.5731731629916696, 1.254272742818995), tolerance = 1e-15)
  g10 <- modelcompR:::.gcmaes_parameters(10)
  expect_equal(c(g10$lam, g10$mu), c(11, 6))
  expect_equal(g10$w, c(0.40240294281871264, 0.25338908403288657, 0.1662215645554205, 0.10437522524706051,
                        0.0564034775763251, 0.017207705769594465), tolerance = 1e-15)
  expect_equal(unlist(g10[c("mueff", "cc", "cs", "c1", "cmu", "ds", "chin")], use.names = FALSE),
               c(3.729458934303068, 0.2857142857142857, 0.34247724070472163, 0.015218446463090412,
                 0.027043953976393947, 1.3424772407047216, 3.0847265651690123), tolerance = 1e-15)
  g1 <- modelcompR:::.gcmaes_parameters(1)
  expect_equal(c(g1$lam, g1$mu), c(4, 2))
  expect_equal(unlist(g1[c("mueff", "cc", "cs", "c1", "cmu", "ds", "chin")], use.names = FALSE),
               c(1.4597898888525862, 0.8, 0.6336855372248923, 0.2963055195692892, 0.027690800889192065,
                 1.6336855372248924, 0.7976190476190477), tolerance = 1e-15)
  # tools/pycma-probe.py
  p2 <- modelcompR:::.pycma_parameters(2)
  expect_equal(c(p2$lam, p2$mu, p2$lam_mirr, p2$maxiter, p2$tolstagnation), c(6, 3, 0, 1620, 146))
  expect_equal(p2$w, c(0.6370425712412167, 0.2845702574380329, 0.07838717132075033, -0.286383782596553,
                       -0.7649580940851276, -1.1559817781589214), tolerance = 1e-15)
  expect_equal(unlist(p2[c("mueff", "cc", "c1", "cmu", "cs", "damps", "chin")], use.names = FALSE),
               c(2.0286114646100626, 0.6245545390268264, 0.1548153998964136, 0.0855927794266643, 0.5, 1.5,
                 1.2533141373155003), tolerance = 1e-15)
  p10 <- modelcompR:::.pycma_parameters(10)
  expect_equal(c(p10$lam, p10$mu, p10$maxiter, p10$tolstagnation), c(10, 5, 7775, 389))
  expect_equal(p10$w, c(0.45627264690340585, 0.2707530970017851, 0.16223111715866975, 0.08523354710016445,
                        0.02550959183597477, -0.08001260758087222, -0.2217641609991465, -0.34455494178482093,
                        -0.45286408637842174, -0.5497499176973852), tolerance = 1e-15)
  expect_equal(unlist(p10[c("mueff", "cc", "c1", "cmu", "cs", "damps", "chin")], use.names = FALSE),
               c(3.1672992814107035, 0.29499038303562225, 0.015283824524751714, 0.02355177665041749,
                 0.3196142529106335, 1.3196142529106334, 3.084327759799864), tolerance = 1e-15)
  p1 <- modelcompR:::.pycma_parameters(1)
  expect_equal(c(p1$lam, p1$mu, p1$lam_mirr, p1$maxiter, p1$tolstagnation), c(4, 2, 1, 1300, 125))
  expect_equal(p1$w, c(0.8041628599327295, 0.19583714006727054, -0.5500162853285434, -1.417877593820792),
               tolerance = 1e-15)
  expect_equal(unlist(p1[c("mueff", "cc", "c1", "cmu", "cs", "damps", "chin")], use.names = FALSE),
               c(1.4597898888525862, 0.6894039888609658, 0.19753701304619278, 0.07549290832280986, 0.5,
                 1.163779122661287, 0.7978845608028655), tolerance = 1e-15)
  p20 <- modelcompR:::.pycma_parameters(2, 20)
  expect_equal(c(p20$maxiter, p20$tolstagnation), c(938, 114))
  expect_equal(unlist(p20[c("mueff", "cc", "c1", "cmu", "damps")], use.names = FALSE),
               c(5.938804235601241, 0.5837604822280295, 0.11884385675893783, 0.39721292501393574,
                 2.0661395741726496), tolerance = 1e-15)
  b <- modelcompR:::.pycma_box(c(-1, -1), c(0.8, 0.8), 2)
  tf <- function(x) modelcompR:::.pycma_transform(b, x)
  expect_equal(tf(c(-1.03, 0.79)), c(-0.998, 0.782), tolerance = 1e-15)
  expect_equal(tf(c(0.83, -1.2)), c(0.798, -0.9000000000000001), tolerance = 1e-15)
  expect_equal(tf(c(5.0, -7.3)), c(0.5000000000000003, 0.3000000000000007), tolerance = 1e-15)
  expect_equal(tf(c(0.76, -0.97)), c(0.7595, -0.968), tolerance = 1e-15)
  b2 <- modelcompR:::.pycma_box(c(-10, 2), c(Inf, 3), 2)
  expect_equal(modelcompR:::.pycma_transform(b2, c(-12, 3.4)), c(-9, 2.8958333333333335), tolerance = 1e-15)
})

test_that("CMAEvolutionStrategyOpt solves the bounded Rosenbrock problem", {
  prob <- OptimizationProblem(cma_rosen, c(0, 0), c(1, 100), lb = c(-1, -1), ub = c(0.8, 0.8))
  l1 <- cma_rosen(c(0, 0), c(1, 100))
  sol <- solve(prob, CMAEvolutionStrategyOpt(), seed = 1)
  expect_lt(10 * sol$objective, l1)
  expect_equal(sol$retcode, "Success")
  expect_true(all(sol$u >= c(-1, -1) & sol$u <= c(0.8, 0.8)))
  expect_equal(solve(prob, CMAEvolutionStrategyOpt(), maxiters = 2, seed = 2)$retcode, "MaxIters")
  expect_equal(solve(prob, CMAEvolutionStrategyOpt(), popsize = 20, maxiters = 5, seed = 3)$original$p$lambda, 20)
  noisy <- solve(prob, CMAEvolutionStrategyOpt(), noise_handling = TRUE, seed = 4)
  expect_lt(10 * noisy$objective, l1)
  expect_gt(noisy$original$p$sigma$sigma, 0)
})

test_that("GCMAESOpt solves Rosenbrock with and without a gradient and needs bounds", {
  lb <- c(-1, -1); ub <- c(1, 1)
  l1 <- cma_rosen(c(0, 0), c(1, 100))
  fg <- OptimizationFunction(cma_rosen, grad = cma_rgrad)
  set.seed(5)
  sol <- solve(OptimizationProblem(fg, c(0, 0), c(1, 100), lb = lb, ub = ub), GCMAESOpt(), maxiters = 1000)
  expect_lt(10 * sol$objective, l1)
  set.seed(6)
  sol <- solve(OptimizationProblem(cma_rosen, c(0, 0), c(1, 100), lb = lb, ub = ub), GCMAESOpt(), maxiters = 1000)
  expect_lt(10 * sol$objective, l1)
  expect_error(solve(OptimizationProblem(cma_rosen, c(0, 0), c(1, 100)), GCMAESOpt()), "requires box constraints")
  # the cache test of OptimizationGCMAES: the minimizer follows the parameter
  objf <- function(u, p) (p[1] - u[1])^2
  set.seed(7)
  expect_equal(solve(OptimizationProblem(objf, 0, 1, lb = -10, ub = 10), GCMAESOpt())$u, 1, tolerance = 1e-3)
  set.seed(8)
  expect_equal(solve(OptimizationProblem(objf, 0, 2, lb = -10, ub = 10), GCMAESOpt())$u, 2, tolerance = 1e-3)
  set.seed(9)
  mx <- solve(OptimizationProblem(function(u, p) -(p[1] - u[1])^2, 0, 2, lb = -10, ub = 10, sense = MaxSense),
              GCMAESOpt())
  expect_equal(mx$u, 2, tolerance = 1e-3)
})

test_that("PyCMAOpt solves the bounded Rosenbrock problem and honors the callback", {
  prob <- OptimizationProblem(cma_rosen, c(0, 0), c(1, 100), lb = c(-1, -1), ub = c(0.8, 0.8))
  l1 <- cma_rosen(c(0, 0), c(1, 100))
  sol <- solve(prob, PyCMAOpt(), seed = 42)
  expect_lt(10 * sol$objective, l1)
  expect_equal(sol$retcode, "Success")
  expect_true("tolfun" %in% sol$original$stop)
  halt <- solve(prob, PyCMAOpt(), seed = 42, maxiters = 25, callback = function(state, l) state$iter > 10)
  expect_equal(halt$retcode, "Terminated")
  expect_equal(solve(prob, PyCMAOpt(), seed = 42, maxiters = 3)$retcode, "MaxIters")
  sph <- solve(OptimizationProblem(function(x, p) sum(x^2), rep(1, 5), NULL), PyCMAOpt(), seed = 1, ftarget = 1e-8)
  expect_lt(sph$objective, 1e-8)
  expect_equal(sph$retcode, "Success")
})
