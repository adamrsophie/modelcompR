# Regressions for behaviors found while writing the vignettes.

test_that("mtkcompile tears a system without an independent variable as MTK does", {
  # ModelingToolkit: equations 0 ~ y - z + y^3, observed x ~ a - y, z ~ -b + 2x
  tx <- Variable("tx"); ty <- Variable("ty"); tz <- Variable("tz")
  ta <- Parameter("ta"); tb <- Parameter("tb")
  cs <- mtkcompile(System(list(ty ~ ta - tx, tz ~ ty^3 + ty, 0 ~ tz - 2 * tx + tb), name = "c"))
  expect_length(equations(cs), 1L)
  expect_identical(vapply(unknowns(cs), format, ""), "ty")
  expect_identical(vapply(observed(cs), function(e) format(e$lhs), ""), c("tx", "tz"))
  sol <- solve(NonlinearProblem(cs, list(ta ~ 3, tb ~ 1, ty ~ 1)))
  expect_equal(sol[tx] + sol[ty], 3, tolerance = 1e-10)
})

test_that("generated functions return unnamed vectors", {
  t <- t_nounits; D <- D_nounits
  gx <- Variable("gx", t); gy <- Variable("gy", t); gk <- Parameter("gk")
  sys <- mtkcompile(System(list(D(gx) ~ gy * gx, D(gy) ~ -gk * gx), t, name = "g"))
  u <- c(a = 1, b = 2)
  expect_null(names(generate_rhs(sys)(u, c(k = 3), 0)))
  expect_null(names(ODEFunction(sys)$f(u, c(k = 3), 0)))
  expect_null(names(build_explicit_observed_function(sys, list(gx, gy))(u, c(k = 3), 0)))
  expect_null(names(build_explicit_observed_function(sys, gx)(u, c(k = 3), 0)))
})

test_that("generate_jacobian(sparse = TRUE) returns a sparse matrix", {
  t <- t_nounits; D <- D_nounits
  sx <- Variable("sx", t); sy <- Variable("sy", t); sk <- Parameter("sk")
  sys <- mtkcompile(System(list(D(sx) ~ -sk * sx, D(sy) ~ sx^2 - sy), t, name = "sp"))
  Jd <- generate_jacobian(sys)(c(2, 3), 5, 0)
  Js <- generate_jacobian(sys, sparse = TRUE)(c(2, 3), 5, 0)
  expect_s4_class(Js, "sparseMatrix")
  expect_identical(length(Js@x), 3L)
  expect_equal(as.matrix(Js), Jd)
})

test_that("named() keeps the parent scope of variables inside expression arguments", {
  t <- t_nounits; D <- D_nounits
  Decay <- function(name, k = 1) {
    nx <- Variable("nx", t, default = 1)
    System(list(D(nx) ~ -k * nx), t, name = name)
  }
  nkg <- Parameter("nkg", default = 0.5)
  named(nd <- Decay(k = 2 * nkg))
  top <- flatten(System(list(), t, systems = list(nd), name = "top"))
  expect_identical(vapply(parameters(top), format, ""), "nkg")
  expect_identical(format(equations(top)[[1L]]$rhs), format(-2 * nkg * Variable("nd₊nx", t)))
})

test_that("latexify renders equations and systems", {
  t <- t_nounits; D <- D_nounits
  lx <- Variable("lx", t); lk <- Parameter("lk")
  one <- latexify(D(lx) ~ -lk * lx)
  expect_identical(one, paste(latex(D(lx)), "=", latex(-lk * lx)))
  expect_identical(latexify(Equation(D(lx), -lk * lx)), one)
  al <- latexify(System(list(D(lx) ~ -lk * lx), t, name = "s"))
  expect_match(al, "^\\\\begin\\{align\\}\n")
  expect_match(al, " &= ", fixed = TRUE)
})

test_that("plot methods accept axis labels and fit the ensemble quantile band", {
  sol <- solve(ODEProblem(function(u, p, t) -u, c(1, 2), c(0, 1)), Tsit5())
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(sol, xlab = "time", ylab = "u", col = c(2, 4)))
  expect_no_error(plot(sol, idxs = list(1, 2), ylab = "second"))
  es <- structure(list(t = 1:3, u = matrix(c(0, 0, 0)), qlow = matrix(c(-1, -2, -3)),
                       qhigh = matrix(c(1, 2, 3)), num_monte = 2L), class = "EnsembleSummary")
  plot(es)
  expect_identical(graphics::par("usr")[3:4] <= c(-3, 3), c(TRUE, FALSE))
})

test_that("a time dependent variable cannot take the name of a symcompR function", {
  t <- t_nounits
  expect_error(Variable("gamma", t), "name of a symcompR function")
  expect_error(Variable("yn", t), "name of a symcompR function")
  expect_no_error(Parameter("gamma"))
})

test_that("Ipopt counts objective and gradient evaluations under AD", {
  rb <- function(x, p) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  sol <- solve(OptimizationProblem(OptimizationFunction(rb, AutoForwardDiff()), c(0, 0)), IpoptOptimizer())
  expect_gt(sol$stats$gevals, sol$stats$iterations)
  expect_gte(sol$stats$fevals, sol$stats$gevals)
})

test_that("solve(LinearProblem) defaults to LU and falls back to QR on a singular matrix", {
  s <- solve(LinearProblem(matrix(c(4, 1, 2, 3), 2), c(1, 2)))
  expect_identical(s$alg$name, "LUFactorization")
  expect_equal(s$u, c(-0.1, 0.7))
  s <- solve(LinearProblem(matrix(c(1, 2, 2, 4), 2), c(1, 2)))
  expect_identical(s$alg$name, "QRFactorization")
  expect_equal(s$u, c(0.2, 0.4))
})

test_that("show_trace and store_trace work for the core nonlinear solvers", {
  prob <- NonlinearProblem(function(u, p) c(u[1]^2 + u[2]^2 - 2, u[1] - u[2]), c(2, 0.5))
  out <- capture.output(s <- solve(prob, NewtonRaphson(), show_trace = TRUE))
  expect_true(any(grepl("^Algorithm: NewtonRaphson", out)))
  expect_true(any(grepl("^Final", out)))
  for (alg in list(NewtonRaphson(), SimpleTrustRegion(), Broyden(), Klement(), LevenbergMarquardt())) {
    s <- solve(prob, alg, store_trace = TRUE)
    h <- s$trace$history
    expect_identical(h[[1L]]$iteration, 0L)
    expect_gt(length(h), 1L)
    # Klement stops at MaxIters and returns its best iterate, not the last
    if (s$retcode == "Success") expect_equal(h[[length(h)]]$fnorm, max(abs(s$resid)), tolerance = 1e-12)
  }
})

test_that("stall and root step counts follow NonlinearSolve", {
  # NonlinearSolve.jl: Newton and TrustRegion stall after 33 steps at u = 0
  # (the 32 step stall rule of AbsNormSafeBestTerminationMode); a root start
  # still takes one step; the default polyalgorithm accumulates the steps of
  # NewtonRaphson and the two TrustRegion methods
  cube <- NonlinearProblem(function(u, p) u^3 - p, c(0), 8)
  s <- solve(cube, NewtonRaphson())
  expect_identical(c(s$retcode, s$stats$nsteps), c("Stalled", "33"))
  s <- solve(cube, TrustRegion())
  expect_identical(c(s$retcode, s$stats$nsteps), c("Stalled", "33"))
  root <- NonlinearProblem(function(u, p) u^3 - p, c(2), 8)
  expect_identical(solve(root, NewtonRaphson())$stats$nsteps, 1L)
  expect_identical(solve(root, TrustRegion())$stats$nsteps, 1L)
  # NewtonRaphson then Broyden: 33 + 21 steps, as Julia's NLStats(54, ...)
  s <- solve(cube, NonlinearSolvePolyAlgorithm(list(NewtonRaphson(), Broyden())))
  expect_identical(s$retcode, "Success")
  expect_identical(s$stats$nsteps, 54L)
})

test_that("the quasi Newton, spectral and Levenberg Marquardt solvers follow NonlinearSolve", {
  # julia: circle x^2 + y^2 = 2, x = y with NonlinearSolve defaults
  f <- function(u, p) c(u[1]^2 + u[2]^2 - 2, u[1] - u[2])
  s <- solve(NonlinearProblem(f, c(1, 2)), LimitedMemoryBroyden())
  expect_identical(c(s$retcode, s$stats$nsteps), c("Unstable", "19"))
  s <- solve(NonlinearProblem(f, c(2, 0.5)), LimitedMemoryBroyden())
  expect_identical(s$stats$nsteps, 11L)
  expect_equal(s$u, c(0.999999999999997, 0.9999999999999968), tolerance = 1e-15)
  s <- solve(NonlinearProblem(f, c(2, 0.5)), DFSane())
  expect_identical(c(s$retcode, s$stats$nsteps), c("Success", "58"))
  for (u0 in list(c(1, 2), c(2, 0.5))) {
    s <- solve(NonlinearProblem(f, u0), LevenbergMarquardt())
    expect_identical(c(s$retcode, s$stats$nsteps), c("Success", "10"))
  }
  s <- solve(NonlinearProblem(f, c(2, 0.5)), PseudoTransient(alpha_initial = 10, mass_matrix = diag(c(2, 1))))
  expect_identical(s$stats$nsteps, 6L)
  # at a root every method steps once; the spectral line search fails on
  # eta = 0 / 0 and the switched evolution relaxation on alpha^-1 * 0 / 0
  codes <- vapply(list(Broyden(), Klement(), LimitedMemoryBroyden(), DFSane(), LevenbergMarquardt(), PseudoTransient()),
                  function(a) { s <- solve(NonlinearProblem(f, c(1, 1)), a); paste(s$retcode, s$stats$nsteps) }, "")
  expect_identical(codes, c("Success 1", "Success 1", "Success 1", "InternalLineSearchFailed 1", "Success 1",
                            "Unstable 1"))
  expect_null(AbsNormSafeBestTerminationMode()$protective_threshold)
})

test_that("generated functions keep the c method of symbolic arguments", {
  t <- t_nounits; D <- D_nounits
  xx <- Variable("xx", t); kk <- Parameter("kk")
  s0 <- mtkcompile(System(list(D(xx) ~ -kk * xx), t, name = "s0"))
  s5 <- modelingtoolkitize(ODEProblem(s0, list(xx ~ 1, kk ~ 2), c(0, 1)))
  expect_length(equations(s5), 1L)
})
