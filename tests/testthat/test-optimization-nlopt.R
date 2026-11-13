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
  # steps of NLopt
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

test_that("the local NLopt algorithms take the steps of NLopt", {
  # julia tools/nlopt-probe.jl (NLopt.jl): the Luksan family, COBYLA, NEWUOA,
  # BOBYQA and SLSQP reproduce NLopt_jll to the last bit on the reference
  # platform (results, return codes and evaluation counts)
  poly3 <- function(x) nl_poly3(x)
  poly3g <- function(x) {
    a <- x[1] - 1; b <- x[2] - 2; c <- x[3] + 1
    c(2 * a + 2 * a * b * b + 0.5 * c, 4 * b + 2 * a * a * b, 6 * c + 0.5 * a)
  }
  rosen <- function(x) nl_rosen(x)
  roseng <- function(x) nl_roseng(x)
  ssum <- function(x) { s <- x[1] + x[2]; for (i in 3:5) s <- s + x[i]; s }
  quart5 <- function(x) {
    s <- ssum(x); v <- 0.1 * s * s * s * s
    for (i in 1:5) { d <- x[i] - 0.1 * i; v <- v + i * d * d }
    v
  }
  quart5g <- function(x) { s <- ssum(x); vapply(1:5, function(i) 2 * i * (x[i] - 0.1 * i) + 0.4 * s * s * s, 0) }
  ballf <- function(x) { s <- 0; for (t in x) s <- s + t * t; s - 4 }
  planef <- function(x) { s <- 0; for (t in x) s <- s + t; s - 1 }
  fns <- list(poly3 = list(poly3, poly3g), rosen = list(rosen, roseng), quart5 = list(quart5, quart5g))
  run <- function(alg, fn, x0, lb = NULL, ub = NULL, cons = "none", xtol = 0, ftol = 0, maxeval = 0) {
    n <- length(x0)
    o <- .nlopt_create(alg, n)
    if (!is.null(lb)) o$lb <- lb
    if (!is.null(ub)) o$ub <- ub
    o$xtol_rel <- xtol; o$ftol_rel <- ftol; o$maxeval <- as.integer(maxeval)
    o$f <- fns[[fn]][[1]]
    o$fg <- function(x) list(val = fns[[fn]][[1]](x), grad = fns[[fn]][[2]](x))
    if (cons != "none") {
      o$fc <- list(list(m = 1L, f = ballf, fg = function(x) list(val = ballf(x), grad = matrix(2 * x, 1)), tol = 1e-8))
      if (cons == "eq") {
        o$h <- list(list(m = 1L, f = planef, fg = function(x) list(val = planef(x), grad = matrix(1, 1, n)), tol = 1e-8))
      }
    }
    .nlopt_optimize(o, x0)
  }
  check <- function(r, x, f, ret, numevals, info) {
    if (ref_platform()) {
      expect_identical(sprintf("%a", r$x), sprintf("%a", as.numeric(x)), info = info)
      expect_identical(sprintf("%a", r$minf), sprintf("%a", as.numeric(f)), info = info)
      expect_identical(r$ret, ret, info = info)
      expect_identical(r$numevals, numevals, info = info)
    } else {
      expect_equal(r$x, as.numeric(x), tolerance = ref_tol(0), info = info)
      expect_equal(r$minf, as.numeric(f), tolerance = ref_tol(0), info = info)
    }
  }
  b5 <- list(rep(-3, 5), rep(3, 5))
  b3 <- list(rep(-5, 3), rep(5, 3))
  check(run("LD_LBFGS", "rosen", c(-1.2, 1), c(-2, -2), c(0.8, 2), ftol = 1e-14, maxeval = 3000),
        c("0x1.999999999999ap-1", "0x1.47ae147ae147cp-1"), "0x1.47ae147ae1478p-5", "SUCCESS", 38L, "LD_LBFGS")
  check(run("LD_VAR1", "rosen", c(-1.2, 1), ftol = 1e-14, maxeval = 3000),
        c("0x1.fffffff60797cp-1", "0x1.ffffffebe101ep-1"), "0x1.9aa565e6608p-60", "SUCCESS", 53L, "LD_VAR1")
  check(run("LD_VAR2", "quart5", rep(1, 5), b5[[1]], b5[[2]], xtol = 1e-9, maxeval = 4000),
        c("-0x1.c6e065b26b565p-4", "0x1.82f63379c076ap-4", "0x1.d652113a2dd51p-3", "0x1.6391f9aa403a1p-2",
          "0x1.d4c6b33d54cbap-2"), "0x1.ac53efc01102ep-3", "SUCCESS", 21L, "LD_VAR2")
  check(run("LD_TNEWTON", "rosen", c(-1.2, 1), ftol = 1e-14, maxeval = 3000),
        c("0x1.fffffffffffa3p-1", "0x1.ffffffffffee7p-1"), "0x1.bce5ap-87", "SUCCESS", 182L, "LD_TNEWTON")
  check(run("LD_TNEWTON_RESTART", "poly3", c(0, 0, 0), xtol = 1e-10, maxeval = 2000),
        c("0x1.0000000000001p+0", "0x1.fffffffffffffp+0", "-0x1.fffffffffffffp-1"), "0x1p-102", "SUCCESS", 19L,
        "LD_TNEWTON_RESTART")
  check(run("LD_TNEWTON_PRECOND", "rosen", c(-1.2, 1), c(-2, -2), c(0.8, 2), ftol = 1e-14, maxeval = 3000),
        c("0x1.999999999999ap-1", "0x1.47ae147ad9ef4p-1"), "0x1.47ae147ae1478p-5", "SUCCESS", 216L,
        "LD_TNEWTON_PRECOND")
  check(run("LD_TNEWTON_PRECOND_RESTART", "quart5", rep(1, 5), b5[[1]], b5[[2]], xtol = 1e-9, maxeval = 4000),
        c("-0x1.c6e065b62c8dcp-4", "0x1.82f6338b5135p-4", "0x1.d652112e70163p-3", "0x1.6391f9a49d79dp-2",
          "0x1.d4c6b33c0302cp-2"), "0x1.ac53efc01102bp-3", "SUCCESS", 29L, "LD_TNEWTON_PRECOND_RESTART")
  check(run("LN_COBYLA", "poly3", c(0, 0, 0), b3[[1]], b3[[2]], cons = "eq", xtol = 1e-10, maxeval = 5000),
        c("0x1.cd1cfa1413d38p-2", "0x1.a003d8a7bff31p+0", "-0x1.134b172cc4e7dp+0"), "0x1.53bdc6f748aa4p-1",
        "XTOL_REACHED", 178L, "LN_COBYLA")
  check(run("LN_NEWUOA", "rosen", c(-1.2, 1), c(-2, -2), c(0.8, 2), ftol = 1e-14, maxeval = 3000),
        c("0x1.000000000000ap+0", "0x1.0000000000014p+0"), "0x1.9p-98", "ROUNDOFF_LIMITED", 466L, "LN_NEWUOA")
  check(run("LN_NEWUOA_BOUND", "rosen", c(-1.2, 1), c(-2, -2), c(0.8, 2), ftol = 1e-14, maxeval = 3000),
        c("0x1.999999999999ap-1", "0x1.47ae1496cb133p-1"), "0x1.47ae147ae151p-5", "FTOL_REACHED", 141L,
        "LN_NEWUOA_BOUND")
  check(run("LN_BOBYQA", "quart5", rep(1, 5), b5[[1]], b5[[2]], xtol = 1e-9, maxeval = 4000),
        c("-0x1.c6e0661df9917p-4", "0x1.82f633ceddd13p-4", "0x1.d652113233c6p-3", "0x1.6391f9a68521ep-2",
          "0x1.d4c6b320e2c11p-2"), "0x1.ac53efc011029p-3", "SUCCESS", 198L, "LN_BOBYQA")
  check(run("LD_SLSQP", "poly3", c(0, 0, 0), b3[[1]], b3[[2]], cons = "eq", xtol = 1e-10, maxeval = 5000),
        c("0x1.cd1cfa12066a8p-2", "0x1.a003d8a80e065p+0", "-0x1.134b172c8fa0fp+0"), "0x1.53bdc6f740f13p-1",
        "XTOL_REACHED", 14L, "LD_SLSQP")
  check(run("LD_SLSQP", "poly3", c(0, 0, 0), b3[[1]], b3[[2]], cons = "ineq", xtol = 1e-10, maxeval = 5000),
        c("0x1.67b94a1321ba4p-1", "0x1.aa65a23bee9e3p+0", "-0x1.b615616f8034bp-1"), "0x1.73bedc8bd3b1bp-2",
        "XTOL_REACHED", 25L, "LD_SLSQP inequality")
})
