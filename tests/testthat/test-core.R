test_that("first order ODE matches the analytic solution", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); tau <- Parameter("tau")
  sys <- mtkcompile(System(list(D(x) ~ (1 - x) / tau), t, name = "lag"))
  prob <- ODEProblem(sys, list(x ~ 0, tau ~ 3), c(0, 10))
  sol <- solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-10)
  expect_equal(sol$retcode, "Success")
  expect_equal(tail(sol[x], 1), 1 - exp(-10 / 3), tolerance = 1e-7)
  expect_equal(sol(5)[[1]], 1 - exp(-5 / 3), tolerance = 1e-6)
  expect_equal(tail(sol[(1 - x) / tau], 1), exp(-10 / 3) / 3, tolerance = 1e-6)
})

test_that("second order equations are lowered with MTK names", {
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t); v <- Variable("v", t); k <- Parameter("k", default = 2)
  sys <- mtkcompile(System(list(D(D(y)) ~ -k * y, v ~ 2 * D(y)), t, name = "osc"))
  expect_equal(vapply(unknowns(sys), format, ""), c("y(t)", "yˍt(t)"))
  expect_equal(format(observed(sys)[[1]]$lhs), "v(t)")
  prob <- ODEProblem(sys, list(y ~ 1, D(y) ~ 0), c(0, 5))
  sol <- solve(prob, reltol = 1e-9, abstol = 1e-9)
  expect_equal(tail(sol[y], 1), cos(sqrt(2) * 5), tolerance = 1e-6)
  expect_equal(tail(sol[v], 1), -2 * sqrt(2) * sin(sqrt(2) * 5), tolerance = 1e-6)
})

test_that("numeric defaults are initial conditions and symbolic ones bindings", {
  t <- t_nounits
  a <- Parameter("a", default = 2)
  b <- Parameter("b", default = 2 * a)
  x <- Variable("x", t, default = 1)
  sys <- System(list(D_nounits(x) ~ -a * b * x), t, name = "s")
  expect_true(smap_has(initial_conditions(sys), a))
  expect_true(smap_has(bindings(sys), b))
  prob <- ODEProblem(mtkcompile(sys), list(), c(0, 1))
  expect_equal(unname(prob$p), 2)
  expect_equal(prob$ps[b], 4)
  prob2 <- remake(prob, p = list(a ~ 3))
  expect_equal(prob2$ps[b], 6)
})

test_that("namespacing follows complete()", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); k <- Parameter("k")
  inner <- System(list(D(x) ~ -k * x), t, name = "r1")
  outer <- System(list(), t, name = "top", systems = list(inner))
  expect_equal(format(outer$r1$x), "top₊r1₊x(t)")
  cs <- complete(outer)
  expect_equal(format(cs$r1$x), "r1₊x(t)")
  expect_length(equations(cs), 1L)
})

test_that("flatten delegates to symcompR for non-systems", {
  expect_equal(flatten(list(1, list(2, 3))), symcompR::flatten(list(1, list(2, 3))))
})

# Regressions found in review of the core.

test_that("explicit solvers reject mass matrices", {
  f <- function(u, p, t) 1
  prob <- ODEProblem(f, u0 = 0, tspan = c(0, 1), mass_matrix = matrix(2))
  expect_error(solve(prob, Tsit5()), "mass matrices")
})

test_that("parameter changes reset the cached derivative", {
  f <- function(u, p, t) p[1]
  cb <- PresetTimeCallback(0.5, function(integ) integ$p[1] <- 2)
  prob <- ODEProblem(f, u0 = 0, tspan = c(0, 1), p = c(1))
  sol <- solve(prob, RK4(), dt = 0.5, callback = cb)
  expect_equal(tail(sol$u, 1)[[1]], 1.5, tolerance = 1e-12)
})

test_that("periodic callbacks fire at every period", {
  n <- 0
  cb <- PeriodicCallback(function(integ) n <<- n + 1, 0.1)
  prob <- ODEProblem(function(u, p, t) 0, u0 = 0, tspan = c(0, 1))
  solve(prob, RK4(), dt = 0.3, callback = cb)
  expect_equal(n, 10)
})

test_that("event root finding converges on steep conditions", {
  cb <- ContinuousCallback(function(u, t, integ) exp(500 * (t - 0.99)) - 1,
                           function(integ) terminate(integ))
  prob <- ODEProblem(function(u, p, t) 1, u0 = 0, tspan = c(0, 2))
  sol <- solve(prob, Tsit5(), callback = cb)
  expect_equal(tail(sol$t, 1), 0.99, tolerance = 1e-8)
})

test_that("terminal affects are saved with saveat", {
  cb <- PresetTimeCallback(1, function(integ) integ$u[1] <- 10)
  prob <- ODEProblem(function(u, p, t) 0, u0 = 0, tspan = c(0, 1))
  sol <- solve(prob, Tsit5(), saveat = 0.5, callback = cb)
  expect_equal(tail(sol$u, 1)[[1]], 10)
})

test_that("backward solutions interpolate", {
  prob <- ODEProblem(function(u, p, t) -u, u0 = 1, tspan = c(0, -2))
  sol <- solve(prob, Tsit5(), reltol = 1e-9, abstol = 1e-12)
  expect_equal(sol(-1)[[1]], exp(1), tolerance = 1e-6)
})

test_that("subsystem derivatives are namespaced", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); k <- Parameter("k")
  sub <- System(list(D(x) ~ -k * x), t, name = "r")
  top <- System(list(), t, name = "top", systems = list(sub))
  eq <- equations(top)[[1]]
  expect_equal(format(eq$lhs), "Derivative(r₊x(t), t)")
  expect_equal(vapply(unknowns(top), format, ""), "r₊x(t)")
})

test_that("observed expressions use the problem parameter layout", {
  t <- t_nounits
  x <- Variable("x", t); y <- Variable("y", t)
  a <- Parameter("a", default = 2); b <- Parameter("b", default = 2 * a)
  sys <- mtkcompile(System(list(D_nounits(x) ~ -x, y ~ b * x), t, name = "s"))
  prob <- ODEProblem(sys, list(x ~ 1), c(0, 1))
  expect_equal(prob[y], 4)
  expect_equal(length(parameters(sys)), 2L)
})

test_that("observed expressions use the parameters in force at each time", {
  t <- t_nounits
  x <- Variable("x", t); p <- Parameter("p", default = 1)
  y <- Variable("y", t)
  sys <- mtkcompile(System(list(D_nounits(x) ~ 0, y ~ p * x), t, name = "s"))
  prob <- ODEProblem(sys, list(x ~ 1), c(0, 2))
  pidx <- 1L
  cb <- PresetTimeCallback(1, function(integ) integ$p[pidx] <- 3)
  sol <- solve(prob, Tsit5(), saveat = c(0.5, 1.5), callback = cb)
  expect_equal(sol[p * x], c(1, 3))
})

test_that("variable kinds stay with their objects", {
  t <- t_nounits
  x <- Variable("zz", t)
  Parameter("zz")
  expect_true(isvariable(x))
})

test_that("parameters bound to numbers are fixed and NA marks solvable ones", {
  t <- t_nounits
  x <- Variable("x", t); a <- Parameter("a"); k <- Parameter("k")
  sys <- mtkcompile(System(list(D_nounits(x) ~ -a * k * x), t, name = "s",
                           bindings = list(a ~ 3)))
  prob <- ODEProblem(sys, list(x ~ 1, k ~ 2), c(0, 1))
  expect_equal(prob$ps[a], 3)
  expect_equal(unname(prob$p), 2)
  sol <- solve(prob, reltol = 1e-9, abstol = 1e-12)
  expect_equal(tail(sol[x], 1), exp(-6), tolerance = 1e-7)
  m <- as_symmap(list(k ~ NA))
  expect_true(is.na(smap_get(m, k)))
})
