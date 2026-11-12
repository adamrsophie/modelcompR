# PolyOpt (OptimizationPolyalgorithms) and Sophia (OptimizationSophia). The
# reference rows come from tools/optimizationpolyalgorithms-probe.jl
# (OptimizationPolyalgorithms 0.3.12: problem, maxiters (0: the default),
# retcode, u and objective in hexadecimal, iterations, fevals and gevals) and
# tools/optimizationsophia-probe.jl (OptimizationSophia 1.4.2: problem,
# settings, maxiters, retcode, u and objective in hexadecimal).
poly_ref <- c(
  "rosen\t5\tMaxIters\t-0x1.0a29bfe6015c9p+0 0x1.16fd63ecf958dp+0\t0x1.0ac347d486007p+2\t5\t7\t7",
  "rosen\t20\tMaxIters\t0x1.774812e3e9ddep-1 0x1.09434ed0ea9ecp-1\t0x1.ba5d4a006c148p-4\t20\t31\t31",
  "rosen\t100\tSuccess\t0x1.0000000023cf4p+0 0x1.0000000045c7bp+0\t0x1.9523bc65p-70\t39\t52\t52",
  "rosen\t0\tSuccess\t0x1.000000000040dp+0 0x1.ffffffffffe8ep-1\t0x1.e7b17adp-76\t41\t59\t59"
)
sophia_ref <- c(
  "rosen\tk1\t1\tSuccess\t-0x1.3328b6d86ec18p+0 0x1.fff822bbecaacp-1\t0x1.8333333333332p+4",
  "rosen\tk1\t2\tSuccess\t-0x1.331e3ac2627f7p+0 0x1.fff045ab63793p-1\t0x1.82bb956600ec6p+4",
  "rosen\tk1\t5\tSuccess\t-0x1.32fec81c7d3f9p+0 0x1.ffd8afaef7896p-1\t0x1.81556ef6e21acp+4",
  "rosen\tk1\t50\tSuccess\t-0x1.312831c7a8afdp+0 0x1.fe77beef581f1p-1\t0x1.6cd5bc233072fp+4",
  "rosen\tk1\t0\tSuccess\t-0x1.0c3834a37c7a7p+0 0x1.e2c3c11437132p-1\t0x1.a65915e1e9c36p+2",
  "rosen\tk1eta\t1\tSuccess\t-0x1.31eb851eb851fp+0 0x1.0147ae147ae14p+0\t0x1.8333333333332p+4",
  "rosen\tk1eta\t2\tSuccess\t-0x1.30a3d70a3d70bp+0 0x1.028f5c28f5c28p+0\t0x1.6b689cb252ce5p+4",
  "rosen\tk1eta\t5\tSuccess\t-0x1.2cccccccccccfp+0 0x1.0666666666664p+0\t0x1.29ede1198aed6p+4",
  "rosen\tk1eta\t50\tSuccess\t-0x1.147ae147ae153p+0 0x1.1c28f5c28f5b8p+0\t0x1.1e1f5c28f5c8p+2",
  "rosen\tk1eta\t0\tSuccess\t0x1.e147ae147adfcp-1 0x1.c7ae147ae1479p-1\t0x1.3626b2f22f20fp-7"
)

poly_rosen <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
poly_roseng <- function(x, p) c(-2.0 * (p[1] - x[1]) - 4.0 * p[2] * (x[2] - x[1]^2) * x[1], 2.0 * p[2] * (x[2] - x[1]^2))
poly_prob <- OptimizationProblem(poly_rosen, c(-1.2, 1), c(1, 100), grad = poly_roseng)

test_that("PolyOpt follows Adam and BFGS of Julia", {
  for (row in strsplit(poly_ref, "\t")) {
    n <- as.integer(row[2])
    info <- paste(row[1:2], collapse = " ")
    s <- if (n > 0L) solve(poly_prob, PolyOpt(), maxiters = n) else solve(poly_prob, PolyOpt())
    expect_identical(unname(as.numeric(s$u)), as.numeric(strsplit(row[4], " ")[[1]]), info = info)
    expect_identical(s$objective, as.numeric(row[5]), info = info)
    expect_identical(s$retcode, row[3], info = info)
    expect_identical(c(s$stats$iterations, s$stats$fevals, s$stats$gevals), as.integer(row[6:8]), info = info)
  }
})

test_that("PolyOpt checks determinism and bounds", {
  noisy <- OptimizationProblem(function(x, p) sum((x - 1)^2) + 1e-3 * stats::rnorm(1), c(0, 0), NULL,
                               grad = function(x, p) 2 * (x - 1))
  expect_error(solve(noisy, PolyOpt()), "deterministic loss functions")
  set.seed(1)
  s <- solve(noisy, PolyOpt(), maxiters = 100)
  expect_identical(s$alg$name, "Adam")
  expect_equal(unname(as.numeric(s$u)), c(1, 1), tolerance = 0.05)
  bprob <- OptimizationProblem(poly_rosen, c(0.1, 0.2), c(1, 100), grad = poly_roseng, lb = c(-0.5, -0.5),
                               ub = c(0.5, 0.8))
  expect_error(solve(bprob, PolyOpt()), class = "IncompatibleOptimizerError")
})

test_that("Sophia with k = 1 follows the iterates of Julia", {
  for (row in strsplit(sophia_ref, "\t")) {
    n <- as.integer(row[3])
    info <- paste(row[1:3], collapse = " ")
    alg <- if (row[2] == "k1") Sophia(k = 1) else Sophia(eta = 0.01, k = 1, lambda = 0, rho = 0.5)
    s <- if (n > 0L) solve(poly_prob, alg, maxiters = n) else solve(poly_prob, alg)
    expect_identical(unname(as.numeric(s$u)), as.numeric(strsplit(row[5], " ")[[1]]), info = info)
    expect_identical(s$objective, as.numeric(row[6]), info = info)
    expect_identical(s$retcode, row[4], info = info)
  }
})

test_that("Sophia with Hessian estimates converges", {
  # randomized (Hutchinson's estimator from rnorm()): only convergence is
  # checked, with the Hessian of the problem and with finite differences
  quad <- function(x, p) sum(c(1, 10) * (x - c(1, -2))^2)
  quadg <- function(x, p) 2 * c(1, 10) * (x - c(1, -2))
  for (hess in list(NULL, function(x, p) diag(c(2, 20)))) {
    set.seed(3)
    prob <- OptimizationProblem(quad, c(0, 0), NULL, grad = quadg, hess = hess)
    s <- solve(prob, Sophia(eta = 0.5, rho = 1, lambda = 0), maxiters = 300)
    expect_lt(s$objective, 1e-8)
    expect_identical(s$retcode, "Success")
  }
  bprob <- OptimizationProblem(poly_rosen, c(0.1, 0.2), c(1, 100), grad = poly_roseng, lb = c(-0.5, -0.5),
                               ub = c(0.5, 0.8))
  expect_error(solve(bprob, Sophia()), class = "IncompatibleOptimizerError")
  expect_error(Sophia(betas = 0.9), "two components")
})
