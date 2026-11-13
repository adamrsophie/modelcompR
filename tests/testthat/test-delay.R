# Delay differential equations with the method of steps.

# y' = -y(t - 1) with history 1: y = 1 - t on [0, 1] and
# y = t^2/2 - 2t + 3/2 on [1, 2], so y(2) = -1/2.
dde_exact <- function(t) ifelse(t <= 1, 1 - t, t^2 / 2 - 2 * t + 1.5)

test_that("y' = -y(t - 1) matches the method of steps polynomial", {
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t); dy <- Variable("dy", t)
  sys <- mtkcompile(System(list(D(y) ~ -EvalAt(t - 1)(y), dy ~ EvalAt(t - 1)(y)), t, name = "dde"))
  expect_true(is_dde(sys))
  expect_length(unknowns(sys), 1L)
  prob <- DDEProblem(sys, list(y ~ 1), c(0, 2))
  expect_equal(prob$constant_lags, 1)
  sol <- solve(prob, MethodOfSteps(Tsit5()), reltol = 1e-8, abstol = 1e-10)
  expect_equal(tail(sol[y], 1), -0.5, tolerance = 1e-9)
  ts <- c(0.3, 0.9, 1.4, 1.8)
  expect_equal(unname(sol(ts)[, 1]), dde_exact(ts), tolerance = 1e-8)
  # delayed observed variables and delayed expressions index the solution
  expect_equal(tail(sol[dy], 1), 0, tolerance = 1e-9)
  expect_equal(tail(sol[EvalAt(t - 0.5)(y)], 1), dde_exact(1.5), tolerance = 1e-8)
  sol2 <- solve(prob, MethodOfSteps(Tsit5(), constrained = TRUE), reltol = 1e-8, abstol = 1e-10)
  expect_equal(tail(sol2[y], 1), -0.5, tolerance = 1e-9)
  expect_s3_class(solve(prob), "DDESolution")
  # any ODE algorithm works inside the method of steps
  sol3 <- solve(prob, MethodOfSteps(RK4()), dt = 0.01)
  expect_equal(tail(sol3[y], 1), -0.5, tolerance = 1e-8)
})

test_that("histories may depend on time and lags may be shorter than steps", {
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t)
  sys <- mtkcompile(System(list(D(y) ~ -EvalAt(t - 1)(y)), t, name = "dde1"))
  # history 1 - t: y' = t - 2 on [0, 1], y = 1 - 2t + t^2/2
  prob <- DDEProblem(sys, list(y ~ 1 - t), c(0, 1))
  expect_equal(prob$h(prob$p, -2), 3)
  sol <- solve(prob, reltol = 1e-9, abstol = 1e-11)
  expect_equal(tail(sol[y], 1), -0.5, tolerance = 1e-9)
  tau <- Parameter("tau", default = 0.01)
  sysb <- mtkcompile(System(list(D(y) ~ -EvalAt(t - tau)(y)), t, name = "dde2"))
  pb <- DDEProblem(sysb, list(y ~ 1), c(0, 0.02))
  expect_equal(pb$constant_lags, 0.01)
  sb <- solve(pb, reltol = 1e-9, abstol = 1e-11)
  expect_equal(tail(sb[y], 1), 0.99 - 0.01 + 0.01^2 / 2, tolerance = 1e-9)
  # a fixed step much longer than a stiff lag: the fixed point iteration
  # fails cleanly instead of accepting an unconverged step
  k <- Parameter("k", default = 50)
  sysc <- mtkcompile(System(list(D(y) ~ -k * EvalAt(t - tau)(y)), t, name = "dde3"))
  pc <- DDEProblem(sysc, list(y ~ 1), c(0, 1))
  sc <- suppressWarnings(solve(pc, MethodOfSteps(RK4()), dt = 0.1, adaptive = FALSE))
  expect_equal(sc$retcode, "ConvergenceFailure")
})

test_that("delay problems from R functions and stochastic delays", {
  rp <- DDEProblem(function(u, h, p, t) -h(p, t - 1), 1, function(p, t) 1, c(0, 2), constant_lags = 1)
  sol <- solve(rp, MethodOfSteps(Tsit5()), reltol = 1e-8, abstol = 1e-10)
  expect_equal(tail(sol[1], 1), -0.5, tolerance = 1e-9)
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t); s <- Parameter("s", default = 0)
  B <- Brownian("ddeB")
  ssys <- mtkcompile(System(list(D(y) ~ -EvalAt(t - 1)(y) + s * B), t, name = "sdde"))
  expect_true(is_dde(ssys))
  sp <- SDDEProblem(ssys, list(y ~ 1), c(0, 2))
  expect_s3_class(sp, "SDDEProblem")
  # without noise Euler-Maruyama is first order accurate
  s0 <- solve(sp, MethodOfSteps(EM()), dt = 0.01, seed = 1)
  expect_lt(abs(tail(s0[y], 1) + 0.5), 0.02)
  sp1 <- SDDEProblem(ssys, list(y ~ 1, s ~ 0.1), c(0, 2))
  s1 <- solve(sp1, MethodOfSteps(EM()), dt = 0.01, seed = 1)
  expect_equal(s1$retcode, "Success")
  expect_error(DDEProblem(ssys, list(y ~ 1), c(0, 2)), "SDDEProblem")
})

# y' = -y(t - 1) with history 1 on [n - 1, n]: sum_k (-1)^k (t - k + 1)^k / k!
dde_exact_sum <- function(t) vapply(t, function(s) {
  k <- 0:(floor(s) + 1)
  sum((-1)^k * (s - k + 1)^k / factorial(k))
}, 0)

dde_stats_ref <- function() {
  f <- test_path("fixtures", "julia", "solver_stats.csv")
  skip_if_not(file.exists(f), "Julia solver statistics not generated")
  utils::read.csv(f)
}

test_that("constant lag discontinuities propagate up to the method's order plus one", {
  prob <- DDEProblem(function(u, h, p, t) -h(p, t - 1), 1, function(p, t) 1, c(0, 8), constant_lags = 1)
  sol <- solve(prob, MethodOfSteps(Tsit5()), reltol = 1e-8, abstol = 1e-10)
  ts <- c(0.5, 1.5, 2.5, 3.7, 5.2, 6.9, 7.9)
  expect_equal(unname(sol(ts)[, 1]), dde_exact_sum(ts), tolerance = 1e-8)
  # DelayDiffEq's error at t = 8 is -7.8e-11
  expect_lt(abs(tail(sol[1], 1) - dde_exact_sum(8)), 2e-10)
  # Tsit5 has order 5: the discontinuities of orders 1 to 6 at t = 1..6
  # are stopping times, t = 7 (order 7) is not (as in DelayDiffEq, whose
  # solution saves exactly 0, 1, ..., 6 and 8 among the integers)
  ints <- sol$t[abs(sol$t - round(sol$t)) < 1e-12]
  expect_equal(ints, c(0:6, 8))
  ref <- dde_stats_ref()
  r <- ref[ref$case == "dde_lag1_Tsit5", ]
  expect_lt(abs(sol$stats$naccept - r$naccept), 0.2 * r$naccept)
  # BS3 (order 3) stops at t = 1..4 only
  sol3 <- solve(prob, MethodOfSteps(BS3()), reltol = 1e-6, abstol = 1e-8)
  expect_true(all(1:4 %in% sol3$t))
  expect_lt(abs(tail(sol3[1], 1) - dde_exact_sum(8)), 1e-6)
  # two lags: every sum of lags up to the order limit is reached
  p2 <- DDEProblem(function(u, h, p, t) -h(p, t - 1) - 0.5 * h(p, t - 1.5), 1, function(p, t) 1, c(0, 6),
                   constant_lags = c(1, 1.5))
  s2 <- solve(p2, MethodOfSteps(Tsit5()), reltol = 1e-8, abstol = 1e-10)
  expect_true(all(seq(1, 6, by = 0.5) %in% s2$t))
  # y = 1 - 1.5 t on [0, 1]; on [1, 1.5] y' = -(1 - 1.5 (t - 1)) - 0.5, so
  # y(1.5) = -0.5 - 0.75 + 1.5 / 8
  expect_equal(s2(1)[1], -0.5, tolerance = 1e-9, ignore_attr = TRUE)
  expect_equal(s2(1.5)[1], -1.0625, tolerance = 1e-9, ignore_attr = TRUE)
  r <- ref[ref$case == "dde_lags2_Tsit5", ]
  expect_lt(abs(s2$stats$naccept - r$naccept), 0.2 * r$naccept)
})

test_that("discontinuity orders: coalescing, neutral problems and the order at t0", {
  st <- function(prob, max_order) .dde_discontinuity_state(prob, max_order, list())
  # lags 1 and 2 reach t = 2 with orders 1 (one lag of 2) and 2 (two of 1):
  # the lowest order propagates
  pr <- list(tspan = c(0, 10), constant_lags = c(1, 2), order_discontinuity_t0 = 0L)
  ds <- st(pr, 1L)
  expect_equal(sort(ds$stops), c(1, 2))
  expect_equal(ds$pending_o, c(1L, 1L))
  .dde_handle_discontinuities(ds, 1)
  expect_equal(sort(ds$pending_t), c(2, 2, 3))
  .dde_handle_discontinuities(ds, 2)
  # order min(1, 2) = 1 at t = 2 propagates with order 2 to 3 and 4
  expect_equal(ds$pending_t[order(ds$pending_t)], c(3, 3, 4))
  expect_equal(ds$tracked_t, c(0, 1, 2))
  .dde_handle_discontinuities(ds, 3)
  # order 2 = max_order + 1 does not propagate further
  expect_equal(ds$pending_t, 4)
  # neutral problems keep the order
  ds <- st(list(tspan = c(0, 10), constant_lags = 1, neutral = TRUE, order_discontinuity_t0 = 0L), 1L)
  for (k in 1:5) .dde_handle_discontinuities(ds, k)
  expect_equal(ds$pending_t, 6)
  expect_equal(ds$pending_o, 0L)
  # a smooth start (order 2 at t0) with a first order method propagates
  # nothing; with order 1 at t0 it propagates once
  expect_length(st(list(tspan = c(0, 10), constant_lags = 1, order_discontinuity_t0 = 2L), 1L)$stops, 0L)
  ds <- st(list(tspan = c(0, 10), constant_lags = 1, order_discontinuity_t0 = 1L), 1L)
  expect_equal(ds$stops, 1)
  .dde_handle_discontinuities(ds, 1)
  expect_length(ds$pending_t, 0L)
  # through the problem constructors
  p <- DDEProblem(function(u, h, p, t) -h(p, t - 1), 1, function(p, t) 1, c(0, 4), constant_lags = 1,
                  neutral = TRUE, order_discontinuity_t0 = 1)
  expect_true(p$neutral)
  expect_equal(p$order_discontinuity_t0, 1L)
  sol <- solve(p, MethodOfSteps(Euler()), dt = 0.25)
  expect_true(all(1:4 %in% sol$t))
  expect_error(DDEProblem(function(u, h, p, t) u, 1, function(p, t) 1, c(0, 1), dependent_lags = list(2)),
               "functions")
})

test_that("state dependent delays: located discontinuities and the analytic solution", {
  # y' = -y(t - 1 - y), history 1: y = 1 - t on [0, 1] (the delayed time
  # stays in the history until t = 1) and y = t - 3 + 2 exp(1 - t) until
  # the next discontinuity at 1 + log(2)
  ex <- function(t) ifelse(t <= 1, 1 - t, t - 3 + 2 * exp(1 - t))
  prob <- DDEProblem(function(u, h, p, t) -h(p, t - 1 - u[1]), 1, function(p, t) 1, c(0, 2),
                     dependent_lags = list(function(u, p, t) 1 + u[1]))
  expect_length(prob$dependent_lags, 1L)
  # an initial step that does not land on t = 1 by chance: the rejected
  # steps around it locate both discontinuities
  sol <- solve(prob, MethodOfSteps(Tsit5()), reltol = 1e-9, abstol = 1e-11, dt = 0.3)
  expect_equal(sol$retcode, "Success")
  ts <- c(0.4, 0.95, 1.05, 1.3, 1.6, 1 + log(2))
  expect_equal(unname(sol(ts)[, 1]), ex(ts), tolerance = 1e-9)
  expect_lt(min(abs(sol$t - 1)), 1e-6)
  expect_lt(min(abs(sol$t - 1 - log(2))), 1e-6)
  # with the default start the first step lands on t = 1 and DelayDiffEq
  # takes the same steps (37 accepted)
  sol2 <- solve(prob, MethodOfSteps(Tsit5()), reltol = 1e-9, abstol = 1e-11)
  ref <- dde_stats_ref()
  r <- ref[ref$case == "dde_sd_Tsit5", ]
  expect_lt(abs(sol2$stats$naccept - r$naccept), 0.2 * r$naccept)
  expect_equal(unname(sol2(ts[1:5])[, 1]), ex(ts[1:5]), tolerance = 1e-9)
  # from a system: the lag 1 + y depends on the state
  t <- t_nounits; D <- D_nounits
  y <- Variable("ddesdy", t)
  sys <- mtkcompile(System(list(D(y) ~ -EvalAt(t - 1 - y)(y)), t, name = "ddesd"))
  sp <- DDEProblem(sys, list(y ~ 1), c(0, 1.6))
  expect_length(sp$constant_lags, 0L)
  expect_length(sp$dependent_lags, 1L)
  expect_equal(sp$dependent_lags[[1]](0.5, sp$p, 0), 1.5)
  ss <- solve(sp, MethodOfSteps(Tsit5()), reltol = 1e-9, abstol = 1e-11, dt = 0.3)
  expect_equal(unname(ss(c(0.5, 1.2, 1.55))[, 1]), ex(c(0.5, 1.2, 1.55)), tolerance = 1e-9)
})

test_that("steps start past a discontinuity and callbacks start new ones", {
  # u' = u(t - 1), history 0, u(0) = 1: u = 1 on [0, 1], u = t on [1, 2];
  # the step from t = 1 takes the derivative from just after t = 0
  pr <- DDEProblem(function(u, h, p, t) h(p, t - 1), 1, function(p, t) 0, c(0, 2), constant_lags = 1)
  sol <- solve(pr, MethodOfSteps(Euler()), dt = 0.5)
  expect_equal(tail(sol[1], 1), 2, tolerance = 1e-12)
  sol <- solve(pr, MethodOfSteps(Tsit5()), reltol = 1e-8, abstol = 1e-10)
  expect_equal(unname(sol(c(0.5, 1.5, 2))[, 1]), c(1, 1.5, 2), tolerance = 1e-9)
  # methods without FSAL start their first stage after the discontinuity:
  # u' = -u(t - 1), history 1, u(0) = 2 gives u(1.5) = 1/8
  pv <- DDEProblem(function(u, h, p, t) -h(p, t - 1), 2, function(p, t) 1, c(0, 1.5), constant_lags = 1)
  for (alg in list(Vern7(), Vern9(), Tsit5())) {
    sv <- solve(pv, MethodOfSteps(alg), dt = 0.5, adaptive = FALSE)
    expect_equal(tail(sv[1], 1), 0.125, tolerance = 1e-12, label = alg$name)
  }
  # u' = -u(t - 1), history 1, and u set to 0 at t = 0.5: the jump
  # propagates to t = 1.5; u(1) = -1/2, u(1.5) = -7/8, u(2) = -3/4
  pc <- DDEProblem(function(u, h, p, t) -h(p, t - 1), 1, function(p, t) 1, c(0, 2), constant_lags = 1)
  cb <- PresetTimeCallback(0.5, function(integ) integ$u[1] <- 0)
  sc <- solve(pc, MethodOfSteps(Tsit5()), reltol = 1e-9, abstol = 1e-11, callback = cb)
  expect_true(1.5 %in% sc$t)
  expect_equal(unname(sc(c(1, 1.5, 2))[, 1]), c(-0.5, -0.875, -0.75), tolerance = 1e-8)
})
