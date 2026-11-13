nl_poly3 <- function(x, p) {
  a <- x[1] - 1; b <- x[2] - 2; c <- x[3] + 1
  a * a + 2 * b * b + 3 * c * c + a * a * b * b + 0.5 * a * c
}
nl_rosen <- function(x, p) { a <- 1 - x[1]; b <- x[2] - x[1] * x[1]; a * a + 100 * b * b }
nl_roseng <- function(x, p) { a <- 1 - x[1]; b <- x[2] - x[1] * x[1]; c(-2 * a - 400 * x[1] * b, 200 * b) }

test_that("Algorithm, NLopt and Opt follow the NLopt interface", {
  expect_identical(Algorithm(11)$name, "LD_LBFGS")
  expect_identical(Algorithm("NLOPT_LN_COBYLA")$value, 25L)
  expect_identical(NLopt$GN_AGS()$value, 43L)
  expect_length(NLopt, 46L)
  expect_true(all(c("srand", "srand_time") %in% names(NLopt)))
  expect_error(Algorithm("LD_FOO"), "unknown NLopt algorithm")
  o <- Opt("LN_NELDERMEAD", 3)
  o$xtol_rel <- 1e-6
  expect_identical(o$algorithm$name, "LN_NELDERMEAD")
  expect_error(solve(OptimizationProblem(nl_rosen, c(-1.2, 1), numeric()), o), "dimension does not match")
  expect_error(solve(OptimizationProblem(nl_rosen, c(-1.2, 1), numeric()), NLopt$LD_LBFGS()), "requires gradients")
})

test_that("the Mersenne Twister of NLopt", {
  .nlopt_init_genrand(5489)
  expect_equal(c(.nlopt_genrand_int32(), .nlopt_genrand_int32(), .nlopt_genrand_int32()),
               c(3499211612, 581869302, 3890346734))
  NLopt$srand(1234)
  a <- c(.nlopt_urand(0, 1), .nlopt_iurand(10))
  NLopt$srand(1234)
  expect_identical(c(.nlopt_urand(0, 1), .nlopt_iurand(10)), a)
})

test_that("stopval and maxtime stop the NLopt algorithms", {
  # stopval: the first point below it ends the search (NLOPT_STOPVAL_REACHED)
  o <- Opt("LN_NELDERMEAD", 3)
  o$stopval <- Inf
  s <- solve(OptimizationProblem(nl_poly3, c(0, 0, 0), numeric()), o)
  expect_identical(s$original$ret, "STOPVAL_REACHED")
  expect_identical(s$retcode, "Success")
  expect_identical(o$numevals, 1L)
  o <- Opt("LN_SBPLX", 3)
  o$stopval <- 1
  s <- solve(OptimizationProblem(nl_poly3, c(0, 0, 0), numeric()), o)
  expect_identical(s$original$ret, "STOPVAL_REACHED")
  expect_lt(s$objective, 1)
  # for MaxSense the stopval is a lower bound on the maximized objective
  o <- Opt("LN_NELDERMEAD", 3)
  o$stopval <- -1
  s <- solve(OptimizationProblem(function(x, p) -nl_poly3(x, p), c(0, 0, 0), numeric(), sense = MaxSense), o)
  expect_identical(s$original$ret, "STOPVAL_REACHED")
  expect_gt(s$objective, -1)
  # maxtime, from the Opt or from solve()
  slow <- function(x, p) { Sys.sleep(0.01); nl_poly3(x, p) }
  s <- solve(OptimizationProblem(slow, c(0, 0, 0), numeric()), NLopt$LN_NELDERMEAD(), maxtime = 0.05)
  expect_identical(s$original$ret, "MAXTIME_REACHED")
  expect_identical(s$retcode, "MaxTime")
  o <- Opt("LN_COBYLA", 3)
  o$maxtime <- 0.05
  s <- solve(OptimizationProblem(slow, c(0, 0, 0), numeric()), o)
  expect_identical(s$retcode, "MaxTime")
  o <- Opt("LN_COBYLA", 3)
  o$stopval <- 5
  s <- solve(OptimizationProblem(nl_poly3, c(0, 0, 0), numeric()), o)
  expect_identical(s$original$ret, "STOPVAL_REACHED")
  fr <- OptimizationFunction(nl_rosen, grad = nl_roseng)
  s <- solve(OptimizationProblem(fr, c(-1.2, 1), numeric()), NLopt$LD_LBFGS(), maxtime = 1e-9)
  expect_identical(s$retcode, "MaxTime")
  # LD_LBFGS_NOCEDAL has no implementation in NLopt: INVALID_ARGS
  s <- solve(OptimizationProblem(fr, c(-1.2, 1), numeric()), NLopt$LD_LBFGS_NOCEDAL())
  expect_identical(s$original$ret, "INVALID_ARGS")
  expect_identical(s$retcode, "Failure")
})

test_that("the NLopt algorithms reach the solutions of NLopt", {
  # julia tools/nlopt-probe.jl (NLopt.jl); the ported algorithms take the
  # steps of NLopt, the stand in kernels (LD_LBFGS, LN_COBYLA) only reach the
  # same points
  o <- Opt("LN_NELDERMEAD", 3)
  o$xtol_rel <- 1e-10
  o$maxeval <- 2000L
  s <- solve(OptimizationProblem(nl_poly3, c(0, 0, 0), numeric()), o)
  expect_equal(s$retcode, "Success")
  expect_equal(unname(s$u), c(1.00000000009593, 2.0000000000848477, -0.9999999999917069), tolerance = 1e-14)
  expect_identical(o$numevals, 251L)
  fr <- OptimizationFunction(nl_rosen, grad = nl_roseng)
  s <- solve(OptimizationProblem(fr, c(-1.2, 1), numeric()), NLopt$LD_LBFGS(), reltol = 1e-14)
  expect_equal(unname(s$u), c(0.9999999999928465, 0.9999999999852197), tolerance = 1e-9)
  s <- solve(OptimizationProblem(fr, c(-1.2, 1), numeric(), lb = c(-2, -2), ub = c(0.8, 2)), Algorithm("LD_LBFGS"),
             reltol = 1e-14)
  expect_equal(unname(s$u), c(0.8, 0.6400000000000001), tolerance = 1e-12)
  expect_equal(s$objective, 0.03999999999999998, tolerance = 1e-12)
  o <- Opt(NLopt$LN_COBYLA(), 3)
  o$xtol_rel <- 1e-10
  o$maxeval <- 5000L
  pc <- OptimizationProblem(nl_poly3, c(0, 0, 0), numeric(), cons = function(x, p) c(sum(x * x), sum(x)),
                            lcons = c(-Inf, 1), ucons = c(4, 1))
  s <- solve(pc, o)
  expect_equal(s$retcode, "Success")
  expect_equal(unname(s$u), c(0.450305841246599, 1.6250586952457444, -1.075364536492343), tolerance = 1e-7)
  expect_equal(s$objective, 0.6635572602138678, tolerance = 1e-7)
})

test_that("the stochastic global algorithms of NLopt (CRS2, ESCH, ISRES)", {
  run <- function(a, seed = 1234, lb = c(-2, -2), ub = c(2, 2)) {
    NLopt$srand(seed)
    o <- Opt(a, 2)
    o$lower_bounds <- lb; o$upper_bounds <- ub
    o$xtol_rel <- 1e-8; o$maxeval <- 1500L
    s <- solve(OptimizationProblem(nl_rosen, c(-1.2, 1), numeric()), o)
    list(s = s, nevals = o$numevals)
  }
  for (a in c("GN_CRS2_LM", "GN_ESCH", "GN_ISRES")) {
    r1 <- run(a); r2 <- run(a)
    # the same seed gives the same run
    expect_identical(r1$s$u, r2$s$u)
    expect_identical(r1$nevals, r2$nevals)
    expect_true(all(r1$s$u >= -2 & r1$s$u <= 2))
    expect_lt(r1$s$objective, nl_rosen(c(-1.2, 1)))
    # a global algorithm needs a finite domain
    expect_identical(run(a, lb = c(-Inf, -2))$s$original$ret, "INVALID_ARGS")
  }
  expect_identical(run("GN_CRS2_LM")$s$original$ret, "XTOL_REACHED")
  expect_identical(run("GN_ESCH")$nevals, 1500L)
  # a population too small for a simplex
  NLopt$srand(1)
  o <- Opt("GN_CRS2_LM", 2)
  o$lower_bounds <- c(-2, -2); o$upper_bounds <- c(2, 2); o$population <- 2L
  expect_identical(solve(OptimizationProblem(nl_rosen, c(0, 0), numeric()), o)$original$ret, "INVALID_ARGS")
})

test_that("MMA, CCSAQ, AUGLAG and MLSL of NLopt", {
  g <- function(x, p) c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2), 200 * (x[2] - x[1]^2))
  fc <- OptimizationFunction(nl_rosen, grad = g, cons = function(x, p) x[1]^2 + x[2]^2,
                             cons_j = function(x, p) matrix(2 * x, 1))
  pc <- OptimizationProblem(fc, c(0, 0), numeric(), lcons = -Inf, ucons = 0.5)
  # the minimizer of the Rosenbrock function in the disk |x|^2 <= 0.5
  xs <- c(0.6054802, 0.3652311)
  for (a in c("LD_MMA", "LD_CCSAQ")) {
    o <- Opt(a, 2)
    o$xtol_rel <- 1e-8; o$maxeval <- 5000L
    s <- solve(pc, o)
    expect_identical(s$original$ret, "XTOL_REACHED", info = a)
    expect_equal(unname(s$u), xs, tolerance = 1e-6, info = a)
    expect_lte(sum(s$u^2), 0.5 + 1e-6)
  }
  # the parameters of the CCSA family are checked as in NLopt
  o <- Opt("LD_MMA", 2)
  o$params$rho_init <- -1
  expect_identical(solve(pc, o)$original$ret, "INVALID_ARGS")
  # AUGLAG with a local optimizer; AUGLAG_EQ refuses one without
  # inequality constraints
  s <- solve(pc, NLopt$LD_AUGLAG(), local_method = NLopt$LN_NELDERMEAD(), xtol_rel = 1e-8, maxiters = 5000)
  expect_equal(unname(s$u), xs, tolerance = 1e-6)
  expect_lte(sum(s$u^2), 0.5 + 1e-6)
  s <- solve(pc, NLopt$LD_AUGLAG_EQ(), local_method = NLopt$LN_NELDERMEAD(), xtol_rel = 1e-8)
  expect_identical(s$original$ret, "INVALID_ARGS")
  s <- solve(pc, NLopt$LD_AUGLAG_EQ(), local_method = NLopt$LD_MMA(), xtol_rel = 1e-8, maxiters = 5000)
  expect_equal(unname(s$u), xs, tolerance = 1e-6)
  # MLSL: a global search in a box with local searches, deterministic
  # under NLopt$srand
  fr <- OptimizationFunction(nl_rosen, grad = g)
  for (a in c("GN_MLSL", "GD_MLSL_LDS")) {
    run <- function() {
      NLopt$srand(1234)
      o <- Opt(a, 2)
      o$lower_bounds <- c(-2, -2); o$upper_bounds <- c(2, 2)
      o$xtol_rel <- 1e-8; o$maxeval <- 1000L
      list(s = solve(OptimizationProblem(fr, c(-1.2, 1), numeric()), o, local_method = NLopt$LD_MMA()),
           nevals = o$numevals)
    }
    r1 <- run(); r2 <- run()
    expect_identical(r1$s$u, r2$s$u)
    expect_identical(r1$nevals, 1000L)
    expect_lt(r1$s$objective, 0.01)
  }
  o <- Opt("G_MLSL", 2)
  o$lower_bounds <- c(-2, -2); o$upper_bounds <- c(2, 2)
  expect_identical(solve(OptimizationProblem(fr, c(0, 0), numeric()), o)$original$ret, "INVALID_ARGS")
})
