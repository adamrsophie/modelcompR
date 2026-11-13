# Interface fixes: remake on plain function problems, preset and periodic
# callbacks, add_tstop, save_idxs in as.data.frame, misuse errors of the
# optimizers, the empty JumpSet and the phase plot.

test_that("remake replaces u0 and p wholesale and merges keywords", {
  prob <- ODEProblem(function(u, p, t) -p$k * u, 1, c(0, 1), list(k = 1))
  pr <- remake(prob, p = list(k = 3))
  expect_identical(pr$p$k, 3)
  expect_equal(tail(solve(pr, Tsit5())$u, 1)[[1]], exp(-3), tolerance = 1e-4)
  prob2 <- ODEProblem(function(u, p, t) -u, c(1, 2), c(0, 1), 2)
  expect_identical(remake(prob2, u0 = c(1, 2, 3))$u0, c(1, 2, 3))
  expect_identical(remake(prob2, u0 = 5)$u0, 5)
  expect_identical(remake(prob2, p = c(2, 3))$p, c(2, 3))
  expect_identical(remake(prob2, u0 = c(4, 5))$u0, c(4, 5))
  prob3 <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
  expect_identical(remake(prob3, p = c(2, 3))$p, c(2, 3))
  # unknown names are solver keywords, as prob.kwargs upstream
  pr3 <- remake(prob3, saveat = c(0, 0.5, 1))
  expect_identical(pr3$kwargs$saveat, c(0, 0.5, 1))
  expect_identical(solve(pr3, Tsit5())$t, c(0, 0.5, 1))
  expect_identical(remake(pr3, kwargs = list())$kwargs, list())
  # a plain function becomes an ODEFunction
  pf <- remake(prob3, f = function(u, p, t) -2 * u)
  expect_s3_class(pf$f, "ODEFunction")
  expect_equal(tail(solve(pf, Tsit5())$u, 1)[[1]], exp(-2), tolerance = 1e-5)
})

test_that("PresetTimeCallback fires at the initial time and takes initialize", {
  hits <- numeric()
  cb <- PresetTimeCallback(c(0, 0.25, 0.5, 1), function(integ) hits <<- c(hits, integ$t))
  sol <- solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), Tsit5(), callback = cb)
  expect_equal(hits, c(0, 0.25, 0.5, 1))
  order <- character()
  cb2 <- PresetTimeCallback(0, function(integ) order <<- c(order, "affect"),
                            initialize = function(integ) order <<- c(order, "init"))
  solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), Tsit5(), callback = cb2)
  expect_identical(order, c("init", "affect"))
  # filter_tstops = FALSE adds every time, and a time before the start is an error
  cb3 <- PresetTimeCallback(c(-1, 0.5), function(integ) NULL, filter_tstops = FALSE)
  expect_error(solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), Tsit5(), callback = cb3),
               "behind the current time")
  hits <- numeric()
  cb4 <- PresetTimeCallback(c(0.5, 2), function(integ) hits <<- c(hits, integ$t), filter_tstops = FALSE)
  sol <- solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), Tsit5(), callback = cb4)
  expect_identical(hits, 0.5)
})

test_that("a symbolic preset event at the initial time fires", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); k <- Parameter("kdose0", default = 1)
  dose <- SymbolicDiscreteCallback(c(0, 2), list(k ~ 2), discrete_parameters = list(k))
  sys <- mtkcompile(System(list(D(x) ~ -k * x), t, name = "dose0", discrete_events = list(dose)))
  sol <- solve(ODEProblem(sys, list(x ~ 1), c(0, 1)), Tsit5(), abstol = 1e-10, reltol = 1e-10)
  expect_equal(tail(sol[x], 1), exp(-2), tolerance = 1e-7)
})

test_that("PeriodicCallback skips the final time unless final_affect", {
  run <- function(...) {
    hits <- numeric()
    cb <- PeriodicCallback(function(integ) hits <<- c(hits, integ$t), 0.25, ...)
    solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), Tsit5(), callback = cb)
    hits
  }
  expect_equal(run(), c(0.25, 0.5, 0.75))
  expect_equal(run(final_affect = TRUE), c(0.25, 0.5, 0.75, 1))
  expect_equal(run(initial_affect = TRUE), c(0, 0.25, 0.5, 0.75))
  expect_equal(run(phase = 0.1), c(0.1, 0.35, 0.6, 0.85))
  calls <- 0
  run(initialize = function(integ) calls <<- calls + 1)
  expect_identical(calls, 1)
  expect_error(PeriodicCallback(function(integ) NULL, 0.25, phase = -1), "non-negative")
})

test_that("add_tstop refuses a time behind the current time", {
  ig <- init(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), Tsit5())
  step_integrator(ig)
  step_integrator(ig)
  expect_error(add_tstop(ig, ig$t - 1e-3), "Tried to add a tstop that is behind the current time")
  expect_null(add_tstop(ig, ig$t + 1e-3))
  igb <- init(ODEProblem(function(u, p, t) -u, 1, c(1, 0)), Tsit5())
  step_integrator(igb)
  expect_error(add_tstop(igb, igb$t + 1e-3), "behind the current time")
})

test_that("as.data.frame keeps only the saved components", {
  prob <- ODEProblem(function(u, p, t) c(-u[1], -2 * u[2]), c(1, 2), c(0, 1))
  df <- as.data.frame(solve(prob, Tsit5(), save_idxs = 2))
  expect_identical(names(df), c("t", "2"))
  expect_identical(names(as.data.frame(solve(prob, Tsit5()))), c("t", "1", "2"))
})

test_that("misused optimizers raise IncompatibleOptimizerError", {
  op <- OptimizationProblem(function(x, p) sum((x - 1)^2), c(0, 0), lb = c(-2, -2), ub = c(2, 2))
  for (alg in list(NSGA2(), NSGA3(), SPEA2(), SMS_EMOA(), WeightedSumScan()))
    expect_error(solve(op, alg), class = "IncompatibleOptimizerError")
  expect_error(solve(OptimizationProblem(function(x, p) sum((x - 1)^2), c(0, 0)), ScipyShgo()),
               "ScipyShgo requires bounds", class = "IncompatibleOptimizerError")
})

test_that("JumpSet() with no jumps is the empty set", {
  expect_s3_class(JumpSet(), "JumpSet")
  expect_length(JumpSet(), 0L)
  expect_length(JumpSet(NULL), 0L)
})

test_that("plot draws a phase plot for a list of two variables", {
  prob <- ODEProblem(function(u, p, t) c(u[2], -u[1]), c(1, 0), c(0, 2 * pi))
  sol <- solve(prob, Tsit5())
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(sol, idxs = list(1, 2)))
  # the axes are the two states, not time
  usr <- graphics::par("usr")
  expect_lt(usr[1L], -0.9)
  expect_invisible(plot(sol, idxs = c(1, 2)))
  expect_gt(graphics::par("usr")[2L], 6)
})
