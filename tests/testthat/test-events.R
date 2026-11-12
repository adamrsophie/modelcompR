# Symbolic events (ModelingToolkitBase callbacks.jl, imperative_affect.jl and
# docs/src/basics/Events.md), checked against analytic solutions.

test_that("bouncing ball: continuous event with Pre", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t, default = 1); v <- Variable("v", t, default = 0)
  sys <- System(list(D(x) ~ v, D(v) ~ -9.8), t, name = "ball",
                continuous_events = list(list(x ~ 0), list(v ~ -Pre(v))))
  expect_length(continuous_events(sys), 1L)
  ball <- mtkcompile(sys)
  sol <- solve(ODEProblem(ball, list(), c(0, 5)), Tsit5())
  expect_equal(sol$retcode, "Success")
  expect_true(min(sol[x]) >= 0 && min(sol[x]) <= 1e-10)
  tb <- sqrt(2 / 9.8)
  ev <- vapply(sol$events, `[[`, 0, "t")
  expect_equal(ev, tb * c(1, 3, 5, 7, 9, 11), tolerance = 1e-6)
  # the velocity is reversed: before and after the first bounce
  i <- which(sol$t == ev[1])
  expect_length(i, 2L)
  expect_equal(sol[v][i], c(-1, 1) * 9.8 * tb, tolerance = 1e-6)
})

test_that("events with several conditions and the %=>% form", {
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t, default = 0); vy <- Variable("vy", t, default = 2)
  ev <- list(y ~ -1.5, y ~ 1.5) %=>% list(vy ~ -Pre(vy))
  sys <- mtkcompile(System(list(D(y) ~ vy, D(vy) ~ 0), t, name = "walls", continuous_events = ev))
  sol <- solve(ODEProblem(sys, list(), c(0, 4)))
  expect_true(max(sol[y]) <= 1.5 + 1e-9 && min(sol[y]) >= -1.5 - 1e-9)
  expect_equal(vapply(sol$events, `[[`, 0, "t"), c(0.75, 2.25, 3.75), tolerance = 1e-8)
  # a list of single formulas is a list of condition only events
  cs <- System(list(D(y) ~ vy, D(vy) ~ 0), t, name = "c", continuous_events = list(y ~ 1, vy ~ 0))
  expect_length(continuous_events(cs), 2L)
  expect_null(continuous_events(cs)[[1]]$affect)
})

test_that("preset time dose changes a parameter and sol[expr] follows it", {
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t, default = 1); k <- Parameter("k", default = 1)
  d <- Parameter("d", t, default = 0)
  dose <- SymbolicDiscreteCallback(list(2), list(d ~ 1), discrete_parameters = d)
  sys <- mtkcompile(System(list(D(y) ~ -k * y + d), t, name = "dose", discrete_events = dose))
  expect_true(.contains_key(parameters(sys), d))
  sol <- solve(ODEProblem(sys, list(), c(0, 4)), Tsit5(), reltol = 1e-8, abstol = 1e-10)
  y2 <- exp(-2)
  i <- which(sol$t == 2)
  expect_length(i, 2L)
  # the rows before and after the event pair the state with the parameters
  # in force with it (ModelingToolkit 11.47 gives both rows the parameters
  # after the event, 10 + y2; this is a deliberate difference)
  expect_equal(sol[d][i], c(0, 1))
  expect_equal(sol[10 * d + y][i], c(y2, 10 + y2), tolerance = 1e-7)
  # ModelingToolkit 11.47: 0.8829803555627492
  expect_equal(tail(sol[y], 1), 1 + (y2 - 1) * exp(-2), tolerance = 1e-7)
  expect_equal(sol(1, idxs = d), 0)
  expect_equal(sol(3, idxs = 10 * d), 10)
  expect_equal(sol$ps[d], 1)
  # the problem itself is unchanged
  expect_equal(ODEProblem(sys, list(), c(0, 4))$ps[d], 0)
})

test_that("a parameter changed without being discrete is an unsolvable affect", {
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t, default = 1); d <- Parameter("d", default = 0)
  sys <- mtkcompile(System(list(D(y) ~ d), t, name = "nd",
                           discrete_events = list(list(1), list(d ~ 1))))
  expect_error(solve(ODEProblem(sys, list(), c(0, 2))), "not solvable")
})

test_that("periodic discrete events, numeric and Clock", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t, default = 0); c <- Discrete("c", t, default = 1)
  sys <- mtkcompile(System(list(D(x) ~ c), t, name = "per",
                           discrete_events = list(1.0, list(c ~ Pre(c) + 1))))
  sol <- solve(ODEProblem(sys, list(), c(0, 3.5)), reltol = 1e-9, abstol = 1e-9)
  expect_equal(tail(sol[x], 1), 1 + 2 + 3 + 4 * 0.5, tolerance = 1e-7)
  expect_equal(sol$ps[c], 4)
  # Clock(dt, phase) ticks on the absolute lattice, at t0 too, never at the end
  d <- Discrete("d", t, default = 0)
  for (ph in c(0, 0.05)) {
    ev <- SymbolicDiscreteCallback(Clock(0.1, phase = ph), list(d ~ Pre(d) + 1))
    s2 <- mtkcompile(System(list(D(x) ~ t), t, name = "clk", discrete_events = ev))
    so <- solve(ODEProblem(s2, list(x ~ 0, d ~ 0), c(0, 1)))
    expect_equal(so$ps[d], 10)
    hist <- .sol_env(so)$param_history
    ticks <- seq(ph, 0.95, by = 0.1)
    expect_equal(hist$t[hist$t > 0], ticks[ticks > 0], tolerance = 1e-9)
  }
})

test_that("symbolic discrete conditions fire at tstops", {
  t <- t_nounits; D <- D_nounits
  N <- Variable("N", t, default = 0); M <- Parameter("M", default = 50)
  tinj <- Parameter("tinject", default = 10); alpha <- Parameter("alpha", default = 100)
  ev <- list(Eq(t, tinj), list(N ~ Pre(N) + M))
  sys <- mtkcompile(System(list(D(N) ~ alpha - N), t, name = "cells", discrete_events = ev))
  expect_true(all(c("M", "tinject") %in% vapply(parameters(sys), format, "")))
  sol <- solve(ODEProblem(sys, list(), c(0, 20)), tstops = 10, reltol = 1e-9, abstol = 1e-9)
  n10 <- 100 * (1 - exp(-10))
  expect_equal(sol[N][sol$t == 10], c(n10, n10 + 50), tolerance = 1e-6)
  # the event does not fire when N >= 50 at that time
  ev2 <- list(Eq(t, tinj) & (N < 50), list(N ~ Pre(N) + M))
  s2 <- mtkcompile(System(list(D(N) ~ alpha - N), t, name = "cells2", discrete_events = ev2))
  sol2 <- solve(ODEProblem(s2, list(), c(0, 20)), tstops = 10, reltol = 1e-9, abstol = 1e-9)
  expect_equal(sum(sol2$t == 10), 1L)
})

test_that("a terminating event stops the integration", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t, default = 0)
  stop_at <- ImperativeAffect(function(m, o, ctx, integ) {
    terminate(integ)
    NULL
  })
  sys <- mtkcompile(System(list(D(x) ~ 1), t, name = "term",
                           continuous_events = list(list(x ~ 2), stop_at)))
  sol <- solve(ODEProblem(sys, list(), c(0, 10)))
  expect_equal(sol$retcode, "Terminated")
  expect_equal(tail(sol$t, 1), 2, tolerance = 1e-9)
})

test_that("affect_neg handles downward crossings separately", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t, default = 0)
  up <- Parameter("up", t, default = 0); dn <- Parameter("dn", t, default = 0)
  ev <- SymbolicContinuousCallback(list(x ~ 0.5), list(up ~ Pre(up) + 1),
                                   affect_neg = list(dn ~ Pre(dn) + 1),
                                   discrete_parameters = list(up, dn))
  sys <- mtkcompile(System(list(D(x) ~ cos(t)), t, name = "neg", continuous_events = ev))
  sol <- solve(ODEProblem(sys, list(), c(0, 10)), reltol = 1e-8, abstol = 1e-10)
  expect_equal(c(sol$ps[up], sol$ps[dn]), c(2, 2))
  expect_equal(vapply(sol$events, `[[`, 0, "t"), c(pi / 6, 5 * pi / 6, 13 * pi / 6, 17 * pi / 6),
               tolerance = 1e-6)
  # affect_neg = NULL: downward crossings still stop the integrator but change nothing
  ev2 <- SymbolicContinuousCallback(list(x ~ 0.5), list(up ~ Pre(up) + 1), affect_neg = NULL,
                                    discrete_parameters = up)
  s2 <- mtkcompile(System(list(D(x) ~ cos(t)), t, name = "neg2", continuous_events = ev2))
  sol2 <- solve(ODEProblem(s2, list(x ~ 0, up ~ 0), c(0, 10)))
  expect_equal(sol2$ps[up], 2)
})

test_that("ImperativeAffect: heater with hysteresis", {
  t <- t_nounits; D <- D_nounits
  temp <- Variable("temp", t)
  on_th <- Parameter("furnace_on_threshold", default = 0.5)
  off_th <- Parameter("furnace_off_threshold", default = 0.7)
  power <- Parameter("furnace_power", default = 1)
  leak <- Parameter("leakage", default = 0.1)
  on <- Discrete("furnace_on", t, default = 0)
  off_ev <- SymbolicContinuousCallback(list(temp ~ off_th),
    ImperativeAffect(function(x, o, c, i) list(furnace_on = 0), modified = list(furnace_on = on)))
  on_ev <- SymbolicContinuousCallback(list(temp ~ on_th),
    ImperativeAffect(function(x, o, c, i) list(furnace_on = 1), modified = list(on)))
  sys <- mtkcompile(System(list(D(temp) ~ on * power - temp^2 * leak), t, name = "heater",
                           continuous_events = list(off_ev, on_ev)))
  sol <- solve(ODEProblem(sys, list(temp ~ 0, on ~ 1), c(0, 10)))
  expect_true(max(sol[temp]) <= 0.7 + 1e-8)
  first <- sol$events[[1]]$t
  expect_true(min(sol[temp][sol$t > first]) >= 0.5 - 1e-8)
  expect_true(length(sol$events) >= 3L)
  # writing an undeclared name is an error
  bad <- SymbolicContinuousCallback(list(temp ~ 0.1),
    ImperativeAffect(function(x, o, c, i) list(zzz = 1), modified = list(on)))
  sb <- mtkcompile(System(list(D(temp) ~ on * power - temp^2 * leak), t, name = "h2",
                          continuous_events = bad))
  expect_error(solve(ODEProblem(sb, list(temp ~ 0, on ~ 1), c(0, 1))), "Tried to write back")
})

test_that("affects are solved with observed equations, and events are namespaced", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t, default = 1); y <- Variable("y", t)
  # setting the observed y to 3 moves x to 1.5
  sys <- mtkcompile(System(list(D(x) ~ 0, y ~ 2 * x), t, name = "obs",
                           discrete_events = list(list(1), list(y ~ 3))))
  sol <- solve(ODEProblem(sys, list(), c(0, 2)))
  expect_equal(tail(sol[x], 1), 1.5)
  expect_equal(tail(sol[y], 1), 3)
  # an event of a subsystem
  inner <- System(list(D(x) ~ -1), t, name = "inner",
                  continuous_events = list(list(x ~ 0), list(x ~ Pre(x) + 1)))
  outer <- System(list(), t, name = "outer", systems = list(inner))
  ce <- continuous_events(outer)
  expect_equal(format(ce[[1]]$conditions[[1]]$lhs), "inner\u208Ax(t)")
  expect_length(continuous_events_toplevel(outer), 0L)
  cs <- mtkcompile(outer)
  so <- solve(ODEProblem(cs, list(), c(0, 2.5)))
  expect_equal(tail(so[cs$inner$x], 1), 0.5, tolerance = 1e-8)
  expect_equal(length(so$events), 2L)
})

test_that("AssignmentAffect evaluates every right-hand side before the event", {
  t <- t_nounits; D <- D_nounits
  a <- Variable("a", t, default = 1); b <- Variable("b", t, default = 2)
  sys <- mtkcompile(System(list(D(a) ~ 0, D(b) ~ 0), t, name = "swap",
                           discrete_events = list(list(1), AssignmentAffect(list(a ~ b, b ~ a)))))
  sol <- solve(ODEProblem(sys, list(), c(0, 2)))
  expect_equal(c(tail(sol[a], 1), tail(sol[b], 1)), c(2, 1))
  expect_length(.modified_unknowns_eq(list(), a ~ 2 * b), 1L)
  expect_length(.modified_unknowns_eq(list(), a ~ b, unknownlist = list(b)), 0L)
})

test_that("parameter history: saved rows, interpolation and backward solves", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t, default = 0); q <- Parameter("q", t, default = 1)
  ev <- SymbolicDiscreteCallback(list(2), list(q ~ 5), discrete_parameters = q)
  sys <- mtkcompile(System(list(D(x) ~ q), t, name = "bw", discrete_events = ev))
  sb <- solve(ODEProblem(sys, list(), c(3, 0)))
  expect_equal(tail(sb[x], 1), -11, tolerance = 1e-8)
  i <- which(sb$t == 2)
  expect_equal(sb[q][i], c(1, 5))
  expect_equal(sb[q, i[1L]], 1)
  expect_equal(sb(c(2.5, 1), idxs = q), c(1, 5))
  # an imperative affect may return nothing to write
  keep <- ImperativeAffect(function(m, o, ctx, integ) list(), modified = list(x = x))
  s2 <- mtkcompile(System(list(D(x) ~ 1), t, name = "keep", continuous_events = list(list(x ~ 1), keep)))
  expect_equal(tail(solve(ODEProblem(s2, list(), c(0, 2)))[x], 1), 2, tolerance = 1e-8)
})

test_that("initialize and finalize affects", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t, default = 0)
  a <- Parameter("a", t, default = 0); b <- Parameter("b", t, default = 0)
  ev <- SymbolicContinuousCallback(list(x ~ 100), NULL, initialize = list(a ~ 3),
                                   finalize = list(b ~ Pre(b) + 1), discrete_parameters = list(a, b))
  sys <- mtkcompile(System(list(D(x) ~ a), t, name = "life", continuous_events = ev))
  sol <- solve(ODEProblem(sys, list(), c(0, 1)))
  expect_equal(sol[a][1], 3)
  expect_equal(tail(sol[x], 1), 3, tolerance = 1e-8)
  expect_equal(sol$ps[b], 1)
})

test_that("affects keep algebraic equations satisfied (Events.md)", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t, default = 1); y <- Variable("y", t, guess = 1)
  p <- Parameter("p", t, default = 2)
  eqs <- list(x * y ~ p, D(x) ~ 0)
  # y compensates the change of x
  s1 <- mtkcompile(System(eqs, t, name = "dae1", discrete_events = list(list(1), list(x ~ Pre(x) + 1))))
  sol <- solve(ODEProblem(s1, list(), c(0, 2)))
  expect_equal(tail(sol[y], 1), 1, tolerance = 1e-6)
  expect_equal(tail(sol[x * y], 1), 2, tolerance = 1e-6)
  # with y held, the discrete parameter p changes instead
  ev <- SymbolicDiscreteCallback(list(1), list(x ~ Pre(x) + 1, y ~ Pre(y)), discrete_parameters = p)
  s2 <- mtkcompile(System(eqs, t, name = "dae2", discrete_events = ev))
  sol2 <- solve(ODEProblem(s2, list(), c(0, 2)))
  expect_equal(tail(sol2[y], 1), 2, tolerance = 1e-6)
  expect_equal(sol2$ps[p], 4, tolerance = 1e-6)
})
