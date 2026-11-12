# SimpleOptimization. The reference rows come from
# tools/simpleoptimization-probe.jl (SimpleOptimization 2.0.2, LineSearch
# 0.1.20), with derivatives that round as ForwardDiff's: problem, method,
# maxiters (0: the default), retcode, u and objective (hexadecimal),
# iterations and fevals.
simple_ref <- c(
  "quad\tSimpleLBFGS\t3\tSuccess\t0x1.3333333333333p-1 -0x1.999999999999ap-1\t-0x1.6666666666666p-1\t3\t0",
  "quad\tSimpleLBFGS\t25\tSuccess\t0x1.3333333333333p-1 -0x1.999999999999ap-1\t-0x1.6666666666666p-1\t3\t0",
  "quad\tSimpleLBFGS\t0\tSuccess\t0x1.3333333333333p-1 -0x1.999999999999ap-1\t-0x1.6666666666666p-1\t3\t0",
  "quad\tSimpleLBFGS3\t3\tSuccess\t0x1.3333333333333p-1 -0x1.999999999999ap-1\t-0x1.6666666666666p-1\t3\t0",
  "quad\tSimpleLBFGS3\t25\tSuccess\t0x1.3333333333333p-1 -0x1.999999999999ap-1\t-0x1.6666666666666p-1\t3\t0",
  "quad\tSimpleLBFGS3\t0\tSuccess\t0x1.3333333333333p-1 -0x1.999999999999ap-1\t-0x1.6666666666666p-1\t3\t0",
  "quad\tSimpleBFGS\t3\tMaxIters\t0x1.2f349615839e1p-1 -0x1.8e7c81fab380ap-1\t-0x1.6632df867836cp-1\t3\t0",
  "quad\tSimpleBFGS\t25\tSuccess\t0x1.3333333333333p-1 -0x1.9999999999999p-1\t-0x1.6666666666666p-1\t25\t0",
  "quad\tSimpleBFGS\t0\tSuccess\t0x1.3333333333333p-1 -0x1.9999999999999p-1\t-0x1.6666666666666p-1\t100\t0",
  "quad\tSimpleNewton\t3\tSuccess\t0x1.3333333333333p-1 -0x1.999999999999ap-1\t-0x1.6666666666666p-1\t3\t0",
  "quad\tSimpleNewton\t25\tSuccess\t0x1.3333333333333p-1 -0x1.999999999999ap-1\t-0x1.6666666666666p-1\t25\t0",
  "quad\tSimpleNewton\t0\tSuccess\t0x1.3333333333333p-1 -0x1.999999999999ap-1\t-0x1.6666666666666p-1\t100\t0",
  "quad\tSimpleGradientDescent\t3\tMaxIters\t0x1.5c1a5870da5dap-1 -0x1.b47cfa26a22b4p-2\t-0x1.0ab190e4f523bp-1\t3\t4",
  "quad\tSimpleGradientDescent\t25\tMaxIters\t0x1.2fb7086fcf517p-1 -0x1.2179ddba1fe5fp-1\t-0x1.4affb9c45b0b5p-1\t25\t26",
  "quad\tSimpleGradientDescent\t0\tMaxIters\t0x1.33332e9c436bap-1 -0x1.9999922c8338bp-1\t-0x1.66666666664c1p-1\t1000\t1001",
  "quad\tSimpleGradientDescent3\t3\tMaxIters\t0x1.655453b26245fp-1 -0x1.9c5b6cc0d11eap-2\t-0x1.f499b063be4dp-2\t3\t4",
  "quad\tSimpleGradientDescent3\t25\tMaxIters\t0x1.5de33514aca54p-1 -0x1.afe3f731702bdp-2\t-0x1.07ab53e50cedp-1\t25\t26",
  "quad\tSimpleGradientDescent3\t0\tMaxIters\t0x1.23345ece00361p-1 -0x1.780ce13e76992p-1\t-0x1.64800a2c1594fp-1\t1000\t1001",
  "quad\tSimpleSOAP\t3\tMaxIters\t0x1.61c295ce92616p-1 -0x1.a2c75a02a476ep-2\t-0x1.ff1ae863a5616p-2\t3\t3",
  "quad\tSimpleSOAP\t25\tMaxIters\t0x1.3fe4abf57fa36p-1 -0x1.e5fb520801953p-2\t-0x1.2b94da4b047p-1\t25\t25",
  "quad\tSimpleSOAP\t0\tMaxIters\t0x1.332e3f8672f58p-1 -0x1.999de65943c0ep-1\t-0x1.666666401e855p-1\t1000\t1000",
  "rosen\tSimpleLBFGS\t3\tMaxIters\t-0x1.021ca77a06897p+0 0x1.01b3b9173f56bp+0\t0x1.02bef67e78ffap+2\t3\t0",
  "rosen\tSimpleLBFGS\t25\tMaxIters\t0x1.c8a41db2c34b7p-1 0x1.959c1c5714a75p-1\t0x1.a16f5163c5d4ap-7\t25\t0",
  "rosen\tSimpleLBFGS\t0\tSuccess\t0x1.fffffffffb49bp-1 0x1.fffffffff6d96p-1\t0x1.dae65d9p-78\t34\t0",
  "rosen\tSimpleLBFGS3\t3\tMaxIters\t-0x1.021ca77a06897p+0 0x1.01b3b9173f56bp+0\t0x1.02bef67e78ffap+2\t3\t0",
  "rosen\tSimpleLBFGS3\t25\tMaxIters\t0x1.b98a0781181eap-1 0x1.7db20bbd54f5dp-1\t0x1.3b971232187e3p-6\t25\t0",
  "rosen\tSimpleLBFGS3\t0\tSuccess\t0x1.0000000003f62p+0 0x1.0000000007eb3p+0\t0x1.f63a4dp-77\t37\t0",
  "rosen\tSimpleBFGS\t3\tMaxIters\t-0x1.12b50fd27493p+0 0x1.c25d600dfc5cap-1\t0x1.760ce75031296p+3\t3\t0",
  "rosen\tSimpleBFGS\t25\tMaxIters\t-0x1.8e0dabb69462dp-1 0x1.352aae5f096f8p-1\t0x1.9465c47ee8afp+1\t25\t0",
  "rosen\tSimpleBFGS\t0\tMaxIters\t-0x1.2de6e3a81b2cdp-2 -0x1.d6425e7bef32p-7\t0x1.59e1f951be9e8p+1\t100\t0",
  "rosen\tSimpleNewton\t3\tMaxIters\t0x1.86e0413f34552p-1 0x1.2a68024c1d94p-1\t0x1.ca7832102d6b5p-5\t3\t0",
  "rosen\tSimpleNewton\t25\tSuccess\t0x1p+0 0x1p+0\t0x0p+0\t25\t0",
  "rosen\tSimpleNewton\t0\tSuccess\t0x1p+0 0x1p+0\t0x0p+0\t100\t0",
  "rosen\tSimpleGradientDescent\t3\tMaxIters\t-0x1.8edba46e18295p+8 0x1.5a8ac20ec5813p+5\t0x1.26794f874ccf2p+41\t3\t4",
  "rosen\tSimpleGradientDescent\t25\tMaxIters\tNaN NaN\tNaN\t25\t26",
  "rosen\tSimpleGradientDescent\t0\tMaxIters\tNaN NaN\tNaN\t1000\t1001",
  "rosen\tSimpleGradientDescent3\t3\tMaxIters\t-0x1.06e1cf0ae142p+0 0x1.0ffb0a2a10a42p+0\t0x1.0754d28ad8f29p+2\t3\t4",
  "rosen\tSimpleGradientDescent3\t25\tMaxIters\t-0x1.0275db8202e1bp+0 0x1.06fd76abb9b19p+0\t0x1.02dffdce24a33p+2\t25\t26",
  "rosen\tSimpleGradientDescent3\t0\tMaxIters\t0x1.4f1df90b0bf7fp-2 0x1.aa09537f3e168p-4\t0x1.d069e9a5afca7p-2\t1000\t1001",
  "rosen\tSimpleSOAP\t3\tMaxIters\t-0x1.30de60375fdd1p+0 0x1.0247db7f45012p+0\t0x1.58dd9cd80b721p+4\t3\t3",
  "rosen\tSimpleSOAP\t25\tMaxIters\t-0x1.1ff249d546bc3p+0 0x1.12d964991a609p+0\t0x1.05da5ad2d5221p+3\t25\t25",
  "rosen\tSimpleSOAP\t0\tMaxIters\t0x1.a78926ffd2dbbp-1 0x1.5dd76a6cf5859p-1\t0x1.eac73370b3e49p-6\t1000\t1000",
  "rosenbox\tSimpleLBFGS\t3\tMaxIters\t0x1.3798994b8843bp-2 0x1.07f177c182a01p-4\t0x1.20660b8ac0d77p-1\t3\t0",
  "rosenbox\tSimpleLBFGS\t0\tFailure\t0x1p-1 0x1.028f5c28f5c29p-2\t0x1.00a3d70a3d70ap-2\t10\t0",
  "rosenbox\tSimpleLBFGS3\t3\tMaxIters\t0x1.3798994b8843bp-2 0x1.07f177c182a01p-4\t0x1.20660b8ac0d77p-1\t3\t0",
  "rosenbox\tSimpleLBFGS3\t0\tSuccess\t0x1p-1 0x1.0000000000001p-2\t0x1p-2\t9\t0",
  "quad\tSimpleLBFGSabstol\t3\tSuccess\t0x1.3333333333333p-1 -0x1.999999999999ap-1\t-0x1.6666666666666p-1\t3\t0",
  "quad\tSimpleLBFGSabstol\t0\tSuccess\t0x1.3333333333333p-1 -0x1.999999999999ap-1\t-0x1.6666666666666p-1\t3\t0",
  "matrix\tSOAPid\t5\tMaxIters\t0x1.ee00beec77a66p-3 0x1.fc7694a9a697bp-3 0x1.fe4c2f419a078p-3 0x1.fedb712fd69bfp-3 0x1.ff192e0902f2ep-3 0x1.ff39460376574p-3\t0x1.b89189249fa98p+1\t5\t6",
  "matrix\tSOAPid\t25\tMaxIters\t0x1.60ffd82792bd2p-3 0x1.59f2c84473012p-1 0x1.01ead7039a2c5p+0 0x1.2522bfa516deap+0 0x1.31d673d02d59bp+0 0x1.3709cc8b92affp+0\t0x1.ab0eb562b8b68p-3\t25\t26",
  "matrix\tSOAPright\t5\tMaxIters\t-0x1.b39082f3cf6cp-4 0x1.722a18122b012p-4 0x1.0b768f7241fbp-2 0x1.42a8d76b1d1dfp-2 0x1.dc9c998d1b00ep-3 0x1.4a2568dd25f0cp-3\t0x1.ea3e51f587422p+1\t5\t6",
  "matrix\tSOAPright\t25\tMaxIters\t0x1.3528b0d3a41c4p-2 0x1.c477d8eade14p-2 0x1.2bf32fda9edfp-1 0x1.1f3d807120cf9p+0 0x1.1ee3186eab2d1p+0 0x1.144f9bf2b7528p+0\t0x1.eecc87b4c9847p-3\t25\t26"
)

simple_quad <- function(x, p) 1.5 * x[1] * x[1] + x[1] * x[2] + x[2] * x[2] - x[1] + x[2]
simple_rosen <- function(x, p) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
# the gradient and Hessian as ForwardDiff evaluates them
simple_rosg <- function(x, p) {
  t2 <- x[2] - x[1] * x[1]
  c(-(2 * (1 - x[1])) + 100 * ((2 * t2) * (-(2 * x[1]))), 100 * (2 * t2))
}
simple_rosh <- function(x, p) {
  t2 <- x[2] - x[1] * x[1]
  h12 <- 100 * (-(4 * x[1]))
  matrix(c(2 + 100 * (8 * (x[1] * x[1]) - 4 * t2), h12, h12, 200), 2, 2)
}
simple_probs <- list(
  quad = OptimizationProblem(simple_quad, c(0.7, -0.4), 1, grad = function(x, p) c(3 * x[1] + x[2] - 1, x[1] + 2 * x[2] + 1),
                             hess = function(x, p) matrix(c(3, 1, 1, 2), 2, 2)),
  rosen = OptimizationProblem(simple_rosen, c(-1.2, 1), 1, grad = simple_rosg, hess = simple_rosh),
  rosenbox = OptimizationProblem(simple_rosen, c(0.1, 0.2), 1, grad = simple_rosg, hess = simple_rosh,
                                 lb = c(-0.5, -0.5), ub = c(0.5, 0.8)),
  matrix = OptimizationProblem(function(x, p) sum((x - p)^2), matrix(0, 3, 2), matrix(1:6, 3, 2) / 4,
                               grad = function(x, p) 2 * (x - p)))
simple_methods <- list(
  SimpleLBFGS = function() SimpleLBFGS(), SimpleLBFGS3 = function() SimpleLBFGS(threshold = 3),
  SimpleBFGS = function() SimpleBFGS(), SimpleNewton = function() SimpleNewton(),
  SimpleGradientDescent = function() SimpleGradientDescent(),
  SimpleGradientDescent3 = function() SimpleGradientDescent(eta = 0.001), SimpleSOAP = function() SimpleSOAP(),
  SimpleLBFGSabstol = function() SimpleLBFGS(), SOAPid = function() SimpleSOAP(eta = 0.05, max_dim = 1),
  SOAPright = function() SimpleSOAP(eta = 0.05, max_dim = 2, freq = 3))

test_that("the SimpleOptimization methods follow the iterates of Julia", {
  for (row in strsplit(simple_ref, "\t")) {
    n <- as.integer(row[3])
    info <- paste(row[1:3], collapse = " ")
    args <- list(simple_probs[[row[1]]], simple_methods[[row[2]]]())
    if (n > 0L) args$maxiters <- n
    if (row[2] == "SimpleLBFGSabstol") args$abstol <- 1e-3
    s <- do.call(solve, args)
    ju <- as.numeric(strsplit(row[5], " ")[[1]])
    jo <- as.numeric(row[6])
    u <- unname(as.numeric(s$u))
    if (all(is.nan(ju))) {
      # the gradient descent with eta = 0.01 diverges on the Rosenbrock
      # function in both
      expect_true(all(is.nan(u)), label = info)
    } else {
      # the matrix rows go through symmetric eigendecompositions (LAPACK
      # rounding); the others agree to the last bit
      tol <- if (row[1] == "matrix") 1e-12 else 0
      expect_lte(max(abs(u - ju) / pmax(1, abs(ju)), abs(s$objective - jo) / max(1, abs(jo))), tol, label = info)
    }
    expect_identical(s$retcode, row[4], info = info)
    expect_identical(as.integer(c(s$stats$iterations, s$stats$fevals)), as.integer(row[7:8]), info = info)
  }
})

test_that("SimpleSOAP with two sided bases converges", {
  # Julia: objective 2.5e-6 after 400 iterations; the null eigendirection of
  # the rank deficient statistics rounds differently, so only convergence is
  # compared
  s <- solve(simple_probs$matrix, SimpleSOAP(eta = 0.05, weight_decay = 0), maxiters = 400)
  expect_lt(s$objective, 1e-5)
  expect_length(s$u, 6L)
})

test_that("SimpleOptimization checks its problems", {
  cprob <- OptimizationProblem(simple_quad, c(0.7, -0.4), 1, cons = function(x, p) x[1] + x[2], lcons = 0,
                               ucons = 0)
  expect_error(solve(cprob, SimpleBFGS()), class = "IncompatibleOptimizerError")
  expect_error(solve(simple_probs$rosenbox, SimpleBFGS()), class = "IncompatibleOptimizerError")
  expect_s3_class(solve(simple_probs$rosenbox, SimpleLBFGS()), "OptimizationSolution")
})
