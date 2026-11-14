# Nonlinear, interval, linear, SCC, steady state and homotopy problems.
# Reference values are analytic.

test_that("a nonlinear system is solved by every algorithm", {
  x <- Variable("x"); y <- Variable("y"); a <- Parameter("a")
  sys <- complete(System(list(0 ~ x^2 + y^2 - a, 0 ~ x - y), name = "ns"))
  prob <- NonlinearProblem(sys, list(x ~ 1, y ~ 2, a ~ 2))
  expect_s3_class(prob, "NonlinearProblem")
  expect_equal(unname(prob$u0), c(1, 2))
  algs <- list(NewtonRaphson(), NewtonRaphson(linesearch = "BackTracking"), SimpleNewtonRaphson(),
               TrustRegion(), SimpleTrustRegion(nlsolve_update_rule = TRUE), SimpleTrustRegion(),
               LevenbergMarquardt(), GaussNewton(), RobustMultiNewton(), FastShortcutNonlinearPolyalg(u0_len = 2),
               TrustRegion(radius_update_scheme = "Fan"), TrustRegion(radius_update_scheme = "Bastin"),
               Broyden(init_jacobian = "true_jacobian"), Klement(init_jacobian = "true_jacobian"))
  for (alg in algs) {
    sol <- solve(prob, alg)
    expect_equal(sol$retcode, "Success", info = alg$name)
    expect_equal(unname(sol$u), c(1, 1), tolerance = 1e-9, info = alg$name)
    expect_lt(max(abs(sol$resid)), 1e-9)
  }
  sol <- solve(prob)
  expect_true(successful_retcode(sol))
  expect_equal(names(sol$u), c("x", "y"))
  expect_equal(sol[x], 1, tolerance = 1e-10)
  expect_equal(sol[x + 2 * y], 3, tolerance = 1e-10)
  expect_equal(sol$ps[a], 2)
  # a different parameter moves the root to sqrt(a / 2)
  sol2 <- solve(remake(prob, p = list(a ~ 8)))
  expect_equal(unname(sol2$u), c(2, 2), tolerance = 1e-9)
})

test_that("the symbolic Jacobian matches the analytic one", {
  x <- Variable("x"); y <- Variable("y"); a <- Parameter("a")
  sys <- complete(System(list(0 ~ x^2 * y - a, 0 ~ sin(x) + y), name = "jac"))
  fn <- NonlinearFunction(sys, jac = TRUE)
  u <- c(0.7, -1.3)
  expect_equal(fn$f(u, 2), c(u[1]^2 * u[2] - 2, sin(u[1]) + u[2]))
  expect_equal(fn$jac(u, 2), rbind(c(2 * u[1] * u[2], u[1]^2), c(cos(u[1]), 1)))
  expect_null(NonlinearFunction(sys)$jac)
})

test_that("raw nonlinear problems and guesses work", {
  sol <- solve(NonlinearProblem(function(u, p) u^2 - p, 1, 2))
  expect_equal(sol$u, sqrt(2), tolerance = 1e-12)
  # Powell's badly scaled start for the Rosenbrock residuals
  prob <- NonlinearProblem(function(u, p) c(10 * (u[2] - u[1]^2), 1 - u[1]), c(-1.2, 1))
  for (alg in list(NewtonRaphson(), TrustRegion(), LevenbergMarquardt()))
    expect_equal(solve(prob, alg)$u, c(1, 1), tolerance = 1e-8)
  # unknowns without a value start from their guesses
  x <- Variable("x", guess = 3); y <- Variable("y")
  sys <- complete(System(list(0 ~ x^2 - 4, 0 ~ y - x), name = "g", guesses = list(y ~ 2 * x)))
  prob <- NonlinearProblem(sys, list())
  expect_equal(unname(prob$u0), c(3, 6))
  expect_equal(unname(solve(prob)$u), c(2, 2), tolerance = 1e-10)
  z <- Variable("z")
  sys2 <- complete(System(list(0 ~ z^2 - 4), name = "g2"))
  expect_error(NonlinearProblem(sys2, list()), "Initial condition underdefined")
  expect_equal(unname(NonlinearProblem(sys2, list(), missing_guess_value = MissingGuessValue$Constant(1))$u0), 1)
  expect_error(NonlinearProblem(complete(System(list(0 ~ x - 1, 0 ~ x - 2), name = "ov")), list(x ~ 0)),
               "as many equations as unknowns")
})

test_that("bounds from metadata keep the iterates in the box", {
  x <- Variable("x", bounds = c(0, 10))
  sys <- complete(System(list(0 ~ x^2 - 4), name = "b"))
  prob <- NonlinearProblem(sys, list(x ~ 0.1))
  expect_equal(prob$lb, 0); expect_equal(prob$ub, 10)
  # Julia (tools/trustregion-probe.jl, sq_newton): NewtonRaphson through the
  # bounds transform stalls at the start, 0.1; the default algorithm finds 2
  snr <- solve(prob, NewtonRaphson())
  expect_equal(unname(snr$u), 0.1)
  expect_equal(snr$retcode, "Stalled")
  expect_equal(unname(solve(prob)$u), 2, tolerance = 1e-10)
  xm <- Variable("xm", bounds = c(-10, 0))
  probm <- NonlinearProblem(complete(System(list(0 ~ xm^2 - 4), name = "bm")), list(xm ~ -0.1))
  expect_equal(unname(solve(probm, TrustRegion())$u), -2, tolerance = 1e-10)
})

test_that("least squares problems report stationary points", {
  # (x - 1, x + 1, y - 3): the minimizer is x = 0, y = 3 with residual (-1, 1, 0)
  prob <- NonlinearLeastSquaresProblem(function(u, p) c(u[1] - 1, u[1] + 1, u[2] - 3), c(5, 5))
  # as NonlinearSolve: the steps stall at the minimizer (StalledSuccess); Gauss
  # Newton lands there in one step and stalls 32 steps later
  for (alg in list(GaussNewton(), LevenbergMarquardt(), FastShortcutNLLSPolyalg())) {
    sol <- solve(prob, alg)
    expect_equal(sol$retcode, "StalledSuccess", info = alg$name)
    if (inherits(alg, "GaussNewton")) expect_identical(sol$stats$nsteps, 33L)
    expect_equal(sol$u, c(0, 3), tolerance = 1e-7, info = alg$name)
    expect_equal(sol$resid, c(-1, 1, 0), tolerance = 1e-7)
  }
  # TrustRegion (More's subproblem) stops when the rejected steps leave the
  # state unchanged: a stalled success at the stationary point
  for (alg in list(TrustRegion(), TrustRegionDogleg())) {
    sol <- solve(prob, alg)
    expect_true(successful_retcode(sol$retcode), info = alg$name)
    expect_equal(sol$u, c(0, 3), tolerance = 1e-7, info = alg$name)
    expect_equal(sol$resid, c(-1, 1, 0), tolerance = 1e-7)
  }
  x <- Variable("x"); y <- Variable("y")
  sys <- complete(System(list(0 ~ x + y - 3, 0 ~ x - y - 1, 0 ~ 2 * x - 4), name = "ls"))
  sol <- solve(NonlinearLeastSquaresProblem(sys, list(x ~ 0, y ~ 0)))
  expect_equal(unname(sol$u), c(2, 1), tolerance = 1e-9)
})

test_that("interval problems are solved by the bracketing methods", {
  x <- Variable("x"); a <- Parameter("a")
  sys <- complete(System(list(0 ~ x^2 - a), name = "sq"))
  prob <- IntervalNonlinearProblem(sys, c(0, 3), list(a ~ 2))
  for (alg in list(ITP(), Bisection(), Brent(), Falsi())) {
    sol <- solve(prob, alg)
    expect_true(successful_retcode(sol), info = alg$name)
    expect_equal(unname(sol$u), sqrt(2), tolerance = 1e-10, info = alg$name)
  }
  sol <- solve(prob)
  expect_equal(sol$alg$name, "ITP")
  expect_equal(sol[x], sqrt(2), tolerance = 1e-12)
  expect_lte(sol$left, sqrt(2)); expect_gte(sol$right, sqrt(2))
  # a root at the left end, and an interval that does not bracket a root
  expect_equal(solve(IntervalNonlinearProblem(function(u, p) u - p, c(1, 2), p = 1))$retcode, "ExactSolutionLeft")
  expect_warning(bad <- solve(IntervalNonlinearProblem(function(u, p) u^2 + 1, c(0, 1))), "enclosing interval")
  expect_equal(bad$retcode, "InitialFailure")
  # cos(x) = x, the Dottie number
  d <- solve(IntervalNonlinearProblem(function(u, p) cos(u) - u, c(0, 1)), Brent())
  expect_equal(d$u, 0.7390851332151607, tolerance = 1e-12)
  y <- Variable("y")
  expect_error(IntervalNonlinearProblem(complete(System(list(0 ~ x - y, 0 ~ x + y), name = "two")), c(0, 1)),
               "single unknown")
})

test_that("linear problems are assembled and solved", {
  x <- Variable("x"); y <- Variable("y"); a <- Parameter("a")
  sys <- complete(System(list(0 ~ 2 * x + y - a, 0 ~ x - y), name = "lin"))
  prob <- LinearProblem(sys, list(a ~ 3))
  expect_equal(prob$A, rbind(c(2, 1), c(1, -1)))
  expect_equal(prob$b, c(3, 0))
  sol <- solve(prob)
  expect_equal(sol$retcode, "Success")
  expect_equal(unname(sol$u), c(1, 1), tolerance = 1e-12)
  expect_equal(sol[x], 1, tolerance = 1e-12)
  expect_equal(unname(solve(remake(prob, p = list(a ~ 6)))$u), c(2, 2), tolerance = 1e-12)
  expect_equal(solve(LinearProblem(rbind(c(2, 0), c(0, 4)), c(2, 4)))$u, c(1, 1))
  nl <- complete(System(list(0 ~ x^2 + y - a, 0 ~ x - y), name = "nl"))
  expect_error(LinearProblem(nl, list(a ~ 1)), "not affine")
})

test_that("SCC problems solve the blocks in dependency order", {
  x <- Variable("x"); y <- Variable("y"); z <- Variable("z"); w <- Variable("w")
  sys <- complete(System(list(0 ~ z^2 + z - y, 0 ~ x^3 - 8, 0 ~ y - 2 * x + w, 0 ~ w - y + 4), name = "chain"))
  prob <- SCCNonlinearProblem(sys, list(x ~ 1, y ~ 1, z ~ 1, w ~ 1))
  # x alone, then the coupled pair (y, w), then z
  expect_equal(lengths(prob$var_sccs), c(1L, 2L, 1L))
  first <- names(prob$u0)[prob$var_sccs[[1L]]]
  expect_equal(first, "x")
  sol <- solve(prob)
  expect_equal(sol$retcode, "Success")
  expect_equal(sol[x], 2, tolerance = 1e-10)
  expect_equal(sol[y], 4, tolerance = 1e-10)
  expect_equal(sol[w], 0, tolerance = 1e-10)
  expect_equal(sol[z], (-1 + sqrt(17)) / 2, tolerance = 1e-10)
  expect_equal(unname(solve(NonlinearProblem(sys, list(x ~ 1, y ~ 1, z ~ 1, w ~ 1)))$u), unname(sol$u),
               tolerance = 1e-8)
  sing <- complete(System(list(0 ~ x - 1, 0 ~ x - y^2 + y^2 - 1), name = "sing", unknowns = list(x, y)))
  expect_error(SCCNonlinearProblem(sing, list(x ~ 0, y ~ 0)), "structurally singular")
})

test_that("the steady state of a decay equation is s / k", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); k <- Parameter("k"); s <- Parameter("s")
  sys <- mtkcompile(System(list(D(x) ~ s - k * x), t, name = "decay"))
  prob <- SteadyStateProblem(sys, list(x ~ 1, k ~ 2, s ~ 3))
  expect_s3_class(prob, "SteadyStateProblem")
  r <- solve(prob, SSRootfind())
  expect_equal(r$retcode, "Success")
  expect_equal(unname(r$u), 1.5, tolerance = 1e-10)
  expect_equal(r[x], 1.5, tolerance = 1e-10)
  d <- solve(prob, DynamicSS(abstol = 1e-10, reltol = 1e-10), reltol = 1e-8, abstol = 1e-10)
  expect_equal(d$retcode, "Success")
  expect_equal(unname(d$u), 1.5, tolerance = 1e-6)
  expect_equal(unname(solve(prob)$u), 1.5, tolerance = 1e-10)
  # a NonlinearProblem from the time dependent system sets the derivatives to zero
  np <- NonlinearProblem(sys, list(x ~ 1, k ~ 4, s ~ 3))
  expect_equal(solve(np)[x], 0.75, tolerance = 1e-10)
  # two species with a conservation law reach the equilibrium from the start
  a <- Variable("a", t); b <- Variable("b", t)
  rs <- mtkcompile(System(list(D(a) ~ -2 * a + b, D(b) ~ 2 * a - b), t, name = "rev"))
  dd <- solve(SteadyStateProblem(rs, list(a ~ 3, b ~ 0)), DynamicSS(abstol = 1e-10, reltol = 1e-10),
              reltol = 1e-8, abstol = 1e-10)
  expect_equal(c(dd[a], dd[b]), c(1, 2), tolerance = 1e-6)
})

test_that("homotopy continuation rescues a guess outside the basin of Newton", {
  y <- Variable("y")
  sys <- complete(System(list(0 ~ homotopy(atan(y - 3), y)), name = "h"))
  expect_true(homotopy_enabled(sys))
  # numerically the operator is its first argument
  expect_equal(homotopy(2, 5), 2)
  expect_equal(NonlinearFunction(sys)$f(4, numeric()), atan(1))
  plain <- solve(NonlinearProblem(sys, list(y ~ 12)), NewtonRaphson(), maxiters = 50)
  expect_false(isTRUE(all.equal(unname(plain$u), 3, tolerance = 1e-6)) && plain$retcode == "Success")
  prob <- HomotopyProblem(sys, list(y ~ 12))
  expect_s3_class(prob, "HomotopyProblem")
  expect_equal(prob$f$f(4, numeric(), 0), 4)
  expect_equal(prob$f$f(4, numeric(), 1), atan(1))
  expect_equal(prob$f$jac(4, numeric(), 0.5)[1, 1], 0.5 + 0.5 / 2)
  sol <- solve(prob)
  expect_equal(sol$retcode, "Success")
  expect_equal(sol[y], 3, tolerance = 1e-9)
  expect_s3_class(AbstractNonlinearProblem(sys, list(y ~ 12)), "HomotopyProblem")
  stripped <- strip_homotopy(sys)
  expect_false(homotopy_enabled(stripped))
  expect_equal(format(equations(stripped)[[1L]]$rhs), format(atan(y - 3)))
  expect_s3_class(AbstractNonlinearProblem(stripped, list(y ~ 2.5)), "NonlinearProblem")
  expect_error(HomotopyProblem(stripped, list(y ~ 1)), "requires a system containing")
})

test_that("limited quantities are corrected at every iterate", {
  v <- Variable("v")
  sys <- complete(System(list(0 ~ exp(limited(v, Min(limitnew, limitold + 1))) - 5), name = "lim"))
  expect_equal(limited(2, 7), 2)
  prob <- NonlinearProblem(sys, list(v ~ 0))
  # the auxiliary unknown of the limited quantity is added
  expect_length(prob$u0, 2L)
  expect_true(is.function(prob$kwargs$postcondition))
  # the corrector never lets the limited unknown grow by more than one
  expect_equal(prob$kwargs$postcondition(c(0, 9), c(0, 2), prob$p), c(0, 3))
  expect_equal(prob$kwargs$postcondition(c(0, 2.5), c(0, 2), prob$p), c(0, 2.5))
  sol <- solve(prob, NewtonRaphson())
  expect_equal(sol$retcode, "Success")
  expect_equal(sol[v], log(5), tolerance = 1e-9)
  expect_error(NonlinearFunction(sys), "not lowered")
})

test_that("polynomial systems are solved by total degree homotopy continuation", {
  x <- Variable("x"); y <- Variable("y")
  sys <- complete(System(list(0 ~ x^2 + y^2 - 5, 0 ~ x * y - 2), name = "poly"))
  prob <- HomotopyContinuationProblem(sys, list(x ~ 0.9, y ~ 2.2))
  expect_equal(prob$degrees, c(2L, 2L))
  sol <- solve(prob)
  expect_equal(sol$retcode, "Success")
  expect_equal(unname(sol$u), c(1, 2), tolerance = 1e-8)
  roots <- sol$original$real_solutions
  expect_equal(nrow(roots), 4L)
  expect_equal(sol$original$n_paths, 4L)
  ord <- order(roots[, 1])
  expect_equal(unname(roots[ord, ]), rbind(c(-2, -1), c(-1, -2), c(1, 2), c(2, 1)), tolerance = 1e-8)
  # a cubic with roots 1, 2, 3 and a parameter
  a <- Parameter("a")
  cub <- complete(System(list(0 ~ (x - 1) * (x - 2) * (x - a)), name = "cubic"))
  cs <- solve(HomotopyContinuationProblem(cub, list(x ~ 2.2, a ~ 3)))
  expect_equal(sort(cs$original$real_solutions[, 1]), c(1, 2, 3), tolerance = 1e-8)
  expect_equal(cs[x], 2, tolerance = 1e-8)
  # x^2 + 1 has no real root
  none <- solve(HomotopyContinuationProblem(complete(System(list(0 ~ x^2 + 1), name = "none")), list(x ~ 0)))
  expect_equal(none$retcode, "ConvergenceFailure")
  expect_equal(nrow(none$original$solutions), 2L)
  # a rational equation: the pole at x = 1 is not a root
  rat <- complete(System(list(0 ~ (x^2 - 1) / (x - 1)), name = "rat"))
  rs <- solve(HomotopyContinuationProblem(rat, list(x ~ 0)))
  expect_equal(unname(rs$original$real_solutions[, 1]), -1, tolerance = 1e-8)
  expect_error(HomotopyContinuationProblem(complete(System(list(0 ~ sin(x) - x^2), name = "np")), list(x ~ 1)),
               "not a polynomial")
})

test_that("bound constrained least squares stops at the constrained minimizer", {
  # the unconstrained Gauss Newton step leaves the box; the minimizers on
  # the boundary are (0.5, 0) and (0, 0.6)
  p1 <- NonlinearLeastSquaresProblem(function(u, p) c(u[1] - 2 * u[2] - 1, u[1] - u[2]), c(0, 0),
                                     lb = c(0, 0), ub = c(Inf, Inf))
  # an algorithm without native bounds runs on the unbounded problem of
  # bounds_transform.jl; julia tools/trustregion-probe.jl (blsq_*): from the
  # corner u0 = (0, 0) the start is nudged by eps^(3/4), where the map is
  # flat, so only TrustRegion leaves it
  # with the exact Jacobian of the probe (finite differences are useless on
  # the flat map at the corner)
  p1j <- NonlinearLeastSquaresProblem(NonlinearFunction(function(u, p) c(u[1] - 2 * u[2] - 1, u[1] - u[2]),
                                                        jac = function(u, p) matrix(c(1, 1, -2, -1), 2)),
                                      c(0, 0), lb = c(0, 0), ub = c(Inf, Inf))
  sol <- solve(p1j, TrustRegion())
  expect_identical(sol$retcode, "StalledSuccess")
  expect_equal(sol$u, c(0.49999999999999994, 2.567789719906645e-27), tolerance = 1e-9)
  for (alg in list(GaussNewton(), LevenbergMarquardt(), SimpleTrustRegion())) {
    sol <- solve(p1, alg)
    expect_true(all(abs(sol$u) < 1e-10), info = alg$name)
  }
  # with the AutoForwardDiff Jacobian the polyalgorithm reaches its
  # TrustRegion, which leaves the corner as in Julia (StalledSuccess at
  # (0.49999999999999994, 2.567789719906645e-27))
  sol <- solve(p1, FastShortcutNLLSPolyalg())
  expect_identical(sol$retcode, "StalledSuccess")
  expect_equal(sol$u, c(0.5, 0), tolerance = 1e-9)
  # the default of a bounded problem is SobolMultistart(FastShortcutBoundedPolyalg())
  expect_equal(solve(p1)$u, c(0.5, 0), tolerance = 1e-6)
  p2 <- NonlinearLeastSquaresProblem(function(u, p) c(u[1] + 2 * u[2] - 1, u[2] - 1), c(0, 0),
                                     lb = c(0, -Inf), ub = c(Inf, Inf))
  expect_equal(solve(p2)$u, c(0, 0.6), tolerance = 1e-6)
})

test_that("edge cases of the problem constructors", {
  # a residual that is not finite at an end of the interval
  expect_warning(s <- solve(IntervalNonlinearProblem(function(u, p) if (u <= 0) NaN else log(u), c(-1, 2))),
                 "enclosing interval")
  expect_equal(s$retcode, "InitialFailure")
  # guesses may refer to parameters bound to expressions
  x <- Variable("x"); a <- Parameter("a"); b <- Parameter("b", default = 2)
  sys <- complete(System(list(0 ~ x^2 - a), unknowns = list(x), ps = list(a, b), name = "gb",
                         bindings = list(a ~ 2 * b), guesses = list(x ~ a)))
  prob <- NonlinearProblem(sys, list())
  expect_equal(unname(prob$u0), 4)
  expect_equal(unname(solve(prob)$u), 2, tolerance = 1e-10)
  # a function object built from a system carries the system
  expect_s3_class(NonlinearProblem(NonlinearFunction(sys), list())$sys, "System")
  # roots that are complex by a small amount are not reported as real
  u <- Variable("u"); v <- Variable("v")
  hs <- complete(System(list(0 ~ u - 100, 0 ~ v^2 + 1e-4), name = "hc"))
  hr <- solve(HomotopyContinuationProblem(hs, list(u ~ 1, v ~ 1)))
  expect_equal(hr$retcode, "ConvergenceFailure")
  expect_equal(nrow(hr$original$solutions), 2L)
})

test_that("the Lorenz steady state matches NonlinearSolve.jl", {
  # tools/julia/reference.jl: NewtonRaphson from (1, 0, 0) with sigma = 10,
  # rho = 26 and beta = 8 / 3 converges to the origin
  f <- test_path("fixtures", "julia", "nonlinear.csv")
  skip_if_not(file.exists(f), "Julia fixture nonlinear.csv is missing")
  ref <- utils::read.csv(f)
  nx <- Variable("nx"); ny <- Variable("ny"); nz <- Variable("nz")
  sig <- Parameter("sig", default = 10); rho <- Parameter("rho", default = 26)
  bn <- Parameter("bn", default = 8 / 3)
  sys <- mtkcompile(System(list(0 ~ sig * (ny - nx), 0 ~ nx * (rho - nz) - ny, 0 ~ nx * ny - bn * nz), name = "ns"))
  sol <- solve(NonlinearProblem(sys, list(nx ~ 1, ny ~ 0, nz ~ 0)), NewtonRaphson())
  expect_equal(sol$retcode, "Success")
  mine <- c(nx = sol[nx], ny = sol[ny], nz = sol[nz])
  expect_lt(max(abs(mine - ref$value[match(names(mine), ref$name)])), 1e-10)
  # from another start Newton finds the nontrivial root (sqrt(beta (rho - 1)), same, rho - 1)
  far <- solve(NonlinearProblem(sys, list(nx ~ 8, ny ~ 8, nz ~ 24)), NewtonRaphson())
  expect_equal(c(far[nx], far[ny], far[nz]), c(sqrt(8 / 3 * 25), sqrt(8 / 3 * 25), 25), tolerance = 1e-9)
})

test_that("Broyden and Klement solve problems where Newton stalls", {
  # the Jacobian 3 u^2 vanishes at the start u = 0, so Newton's method and
  # the trust region methods cannot move; NonlinearSolve returns Stalled
  # for the default polyalgorithm too (start_index 3 for small problems)
  prob <- NonlinearProblem(function(u, p) u^3 - p, 0, 8)
  expect_false(successful_retcode(solve(prob, NewtonRaphson())))
  expect_false(successful_retcode(solve(prob)))
  for (alg in list(Broyden(), Broyden(update_rule = "bad_broyden"), Klement())) {
    sol <- solve(prob, alg)
    expect_equal(sol$retcode, "Success", info = alg$name)
    expect_equal(sol$u, 2, tolerance = 1e-10)
  }
  # the polyalgorithm without the problem size starts at Broyden
  sol <- solve(prob, FastShortcutNonlinearPolyalg())
  expect_equal(sol$stats$alg, "Broyden")
  expect_equal(sol$u, 2, tolerance = 1e-10)
  # more than 25 unknowns: the default starts at Broyden, as in Julia
  big <- NonlinearProblem(function(u, p) u^3 + u - p, rep(0.5, 30), 2)
  sol <- solve(big)
  expect_equal(sol$stats$alg, "Broyden")
  expect_equal(sol$u, rep(1, 30), tolerance = 1e-9)
  # with a Jacobian function the default leaves out the quasi Newton methods
  nlqn_x <- Variable("nlqn_x")
  sys <- complete(System(list(0 ~ nlqn_x^3 - 8), name = "nlqn"))
  expect_false(successful_retcode(solve(NonlinearProblem(sys, list(nlqn_x ~ 0)))))
  expect_equal(unname(solve(NonlinearProblem(sys, list(nlqn_x ~ 0)), Klement())$u), 2, tolerance = 1e-10)
})

test_that("polyalgorithms have the order of NonlinearSolve", {
  algs <- function(a) vapply(a$algs, function(b) paste(b$name, b$radius_update_scheme %||% b$linesearch %||% ""), "")
  expect_equal(algs(FastShortcutNonlinearPolyalg()),
               c("Broyden ", "Klement ", "NewtonRaphson ", "TrustRegion More", "TrustRegion Fan",
                 "LevenbergMarquardt "))
  expect_equal(FastShortcutNonlinearPolyalg()$start_index, 1L)
  expect_equal(FastShortcutNonlinearPolyalg(u0_len = 10)$start_index, 3L)
  expect_equal(FastShortcutNonlinearPolyalg(u0_len = 30)$start_index, 1L)
  expect_equal(algs(FastShortcutNonlinearPolyalg(must_use_jacobian = TRUE)),
               c("NewtonRaphson ", "TrustRegion More", "TrustRegion Fan", "LevenbergMarquardt "))
  expect_equal(algs(RobustMultiNewton()),
               c("TrustRegion More", "TrustRegion Bastin", "NewtonRaphson ", "NewtonRaphson BackTracking",
                 "TrustRegion NLsolve", "TrustRegion Fan"))
  expect_equal(algs(FastShortcutNLLSPolyalg()),
               c("GaussNewton ", "LevenbergMarquardt ", "TrustRegion More", "GaussNewton BackTracking",
                 "TrustRegion Fan", "LevenbergMarquardt "))
  expect_equal(SimpleTrustRegion(nlsolve_update_rule = TRUE)$radius_update_scheme, "NLsolve")
  expect_error(Broyden(init_jacobian = "true_jacobian", update_rule = "diagonal"), "Unknown combination")
  # Powell's badly scaled Rosenbrock start with every radius update
  prob <- NonlinearProblem(function(u, p) c(10 * (u[2] - u[1]^2), 1 - u[1]), c(-1.2, 1))
  for (sch in c("Simple", "NLsolve", "Fan", "Bastin"))
    expect_equal(solve(prob, TrustRegion(radius_update_scheme = sch))$u, c(1, 1), tolerance = 1e-8, info = sch)
})

test_that("SCC problems solve linear and homotopy blocks", {
  nlscc_x <- Variable("nlscc_x"); nlscc_y <- Variable("nlscc_y"); nlscc_z <- Variable("nlscc_z")
  sys <- complete(System(list(0 ~ nlscc_x^3 - 8, 0 ~ nlscc_y - 2 * nlscc_x,
                              0 ~ homotopy(atan(nlscc_z - 3), nlscc_z) + nlscc_y - 4),
                         name = "nlscc"))
  prob <- SCCNonlinearProblem(sys, list(nlscc_x ~ 1, nlscc_y ~ 1, nlscc_z ~ 12))
  expect_equal(prob$scc_kinds, c("nonlinear", "linear", "homotopy"))
  sol <- solve(prob)
  expect_equal(sol$retcode, "Success")
  expect_equal(unname(sol$u[c("nlscc_x", "nlscc_y", "nlscc_z")]), c(2, 4, 3), tolerance = 1e-9)
})

test_that("Klement returns Unstable when its diagonal Jacobian turns NaN", {
  # NonlinearSolve.jl: retcode Unstable after 166 steps at the best iterate
  prob <- NonlinearProblem(function(u, p) c(10 * (u[2] - u[1]^2), 1 - u[1]), c(-1.2, 1))
  sol <- solve(prob, Klement())
  expect_identical(sol$retcode, "Unstable")
  expect_identical(sol$stats$nsteps, 166L)
  expect_equal(sol$u, c(-0.5014300321370808, 0.6507150160685404), tolerance = 1e-12)
})
