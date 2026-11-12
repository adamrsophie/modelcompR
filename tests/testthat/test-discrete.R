# Discrete time systems and clocks (ModelingToolkitBase discretedomain.jl,
# the shift branch of __mtkcompile, discreteproblem.jl, and the
# SampledData tutorial), against analytic recursions.

test_that("Shift, ShiftIndex and shift2term", {
  t <- t_nounits
  x <- Variable("x", t); y <- Variable("y", t)
  expect_true(hasshift(Shift(t, -1)(x)))
  expect_identical(.key(Shift(t, 1)(Shift(t, -1)(x))), .key(x))
  expect_equal(format(shift2term(Shift(t, -2)(x))), "xₜ₋₂(t)")
  expect_identical(.key(shift2term(Shift(t, 1)(shift2term(Shift(t, -1)(x))))), .key(x))
  # shifts distribute over expressions
  e <- Shift(t, -1)(x + 2 * y)
  expect_identical(.key(e), .key(Shift(t, -1)(x) + 2 * Shift(t, -1)(y)))
  k <- ShiftIndex(Clock(0.1))
  expect_identical(.key((k - 1)(x)), .key(k(x, -1)))
  expect_equal(TimeDomain(k(x))$dt, 0.1)
  expect_error(k(x + y), "multivariate")
  expect_true(hassample(Sample(0.1)(x)) && hashold(Hold(x)))
  expect_equal(format(Clock(0.1, phase = 0.05)), "Clock(0.1, phase = 0.05)")
})

test_that("x(k) ~ 0.5 x(k - 1) as a DiscreteProblem", {
  t <- t_nounits
  x <- Variable("x", t)
  k <- ShiftIndex(t)
  sys <- mtkcompile(System(list(x ~ 0.5 * (k - 1)(x)), t, name = "halving"))
  expect_equal(vapply(unknowns(sys), format, ""), "xₜ₋₁(t)")
  expect_equal(format(observed(sys)[[1]]), "x(t) ~ 0.5*xₜ₋₁(t)")
  prob <- DiscreteProblem(sys, list((k - 1)(x) ~ 1), c(0, 5))
  sol <- solve(prob, FunctionMap())
  expect_equal(sol$t, 0:5)
  # ModelingToolkit 11.47 gives the same x and x(k - 1) series
  expect_equal(sol[x], 0.5^(1:6))
  expect_equal(sol[shift2term((k - 1)(x))], 0.5^(0:5))
  # a value of x itself is x at the initial time; the past value follows
  # (ModelingToolkit 11.47: x = 1, 0.5, ... and x(k - 1) = 2, 1, ...)
  sol2 <- solve(DiscreteProblem(sys, list(x ~ 1), c(0, 5)), FunctionMap())
  expect_equal(sol2[x], 0.5^(0:5))
  expect_equal(sol2[shift2term((k - 1)(x))], 2 * 0.5^(0:5))
  expect_error(DiscreteProblem(sys, list(), c(0, 5)), "not provided")
  f <- DiscreteFunction(sys)
  expect_equal(f$f(4, numeric(), 0), 2)
})

test_that("Fibonacci and a forced recursion", {
  t <- t_nounits
  y <- Variable("y", t); u <- Variable("u", t)
  k <- ShiftIndex(t)
  fib <- mtkcompile(System(list(y ~ (k - 1)(y) + (k - 2)(y)), t, name = "fib"))
  expect_length(unknowns(fib), 2L)
  sol <- solve(DiscreteProblem(fib, list((k - 1)(y) ~ 1, (k - 2)(y) ~ 0), c(0, 9)), FunctionMap())
  expect_equal(sol[y], c(1, 2, 3, 5, 8, 13, 21, 34, 55, 89))
  # x ~ x(k - 1) + u with u ~ 1 counts up from 1, as the MTK discrete test
  x <- Variable("x", t)
  de <- mtkcompile(System(list(x ~ (k - 1)(x) + u, u ~ 1), t, name = "de"))
  s2 <- solve(DiscreteProblem(de, list((k - 1)(x) ~ 0), c(0, 10)), FunctionMap())
  expect_equal(s2[x], 1:11)  # as ModelingToolkit 11.47
  # SIR map of the MTK tests, against the direct recursion
  S <- Variable("S", t); I <- Variable("I", t); R <- Variable("R", t)
  beta <- Parameter("beta", default = 0.05); cc <- Parameter("c", default = 10)
  gam <- Parameter("gamma", default = 0.25); dtp <- Parameter("dt", default = 0.1)
  N <- (k - 1)(S) + (k - 1)(I) + (k - 1)(R)
  inf <- (1 - exp(-beta * cc * (k - 1)(I) / N * dtp)) * (k - 1)(S)
  rec <- (1 - exp(-gam * dtp)) * (k - 1)(I)
  sir <- mtkcompile(System(list(S ~ (k - 1)(S) - inf, I ~ (k - 1)(I) + inf - rec,
                                R ~ (k - 1)(R) + rec), t, name = "sir"))
  ss <- solve(DiscreteProblem(sir, list(S ~ 990, I ~ 10, R ~ 0), c(0, 20)), FunctionMap())
  # S, I, R are the values at the initial time; the past is solved for
  u <- c(990, 10, 0); ref <- matrix(NA, 21, 3)
  for (i in 1:21) {
    ref[i, ] <- u
    n <- sum(u)
    a <- (1 - exp(-0.05 * 10 * u[2] / n * 0.1)) * u[1]; b <- (1 - exp(-0.25 * 0.1)) * u[2]
    u <- c(u[1] - a, u[2] + a - b, u[3] + b)
  }
  expect_equal(unname(ss[list(S, I, R)]), ref, tolerance = 1e-10)
})

test_that("positive shifts and mixed systems are rejected", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t)
  k <- ShiftIndex(t)
  expect_error(mtkcompile(System(list(x ~ (k + 1)(x) + (k - 1)(x)), t, name = "pos")), "Positive shifts")
})

test_that("SampledData tutorial: discrete controller and filter on a continuous plant", {
  t <- t_nounits; D <- D_nounits
  dt <- 0.5
  r <- Variable("r", t)
  k <- ShiftIndex(Clock(dt))
  # only x has an initial value: the tutorial also gives y = 0 and u = 0,
  # which contradict y ~ x and u ~ Hold(c.ud) and over determine
  # initialization
  plant <- function(name) {
    x <- Variable("x", t, default = 1); u <- Variable("u", t)
    y <- Variable("y", t)
    System(list(D(x) ~ -x + u, y ~ x), t, name = name)
  }
  filt <- function(name) {
    x <- Variable("x", t, default = 0); u <- Variable("u", t, default = 0)
    y <- Variable("y", t, default = 0)
    a <- 1 / exp(dt)
    System(list(k(x) ~ a * (k - 1)(x) + (1 - a) * k(u), y ~ x), t, name = name)
  }
  controller <- function(kp, name) {
    y <- Variable("y", t); r <- Variable("r", t, default = 0)
    ud <- Variable("ud", t, default = 0); yd <- Variable("yd", t, default = 0)
    kpp <- Parameter("kp", default = kp)
    System(list(yd ~ Sample(y), ud ~ kpp * (r - yd)), t, name = name)
  }
  f <- filt("f"); c <- controller(1, "c"); p <- plant("p")
  conns <- list(r ~ sin(t), f$u ~ r, f$y ~ c$r, Hold(c$ud) ~ p$u, p$y ~ c$y)
  cl <- mtkcompile(System(conns, t, systems = list(f, c, p), name = "cl"))
  # the discrete variables are parameters updated at the ticks
  expect_true(.contains_key(parameters(cl), cl$c$ud))
  sol <- solve(ODEProblem(cl, list(), c(0, 5)), reltol = 1e-10, abstol = 1e-12)
  # exact zero order hold recursion: ticks at 0, 0.5, ..., 4.5
  a <- exp(-dt); xp <- 1; fx <- 0; ud <- numeric()
  for (tk in seq(0, 4.5, by = dt)) {
    fx <- a * fx + (1 - a) * sin(tk)
    u <- fx - xp
    ud <- c(ud, u)
    xp <- u + (xp - u) * exp(-dt)
  }
  expect_equal(tail(sol[cl$p$x], 1), xp, tolerance = 1e-7)
  expect_equal(sol(seq(0.25, 4.75, by = dt), idxs = cl$c$ud), ud, tolerance = 1e-9)
  expect_equal(sol(4.75, idxs = cl$p$u), ud[10], tolerance = 1e-9)
})

test_that("clock inference rejects inconsistent clocks", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t)
  k1 <- ShiftIndex(Clock(0.1)); k2 <- ShiftIndex(Clock(0.2))
  # a discrete equation on a continuous variable
  expect_error(mtkcompile(System(list(D(x) ~ -x, y ~ (k1 - 1)(y) + x), t, name = "bad")),
               "multiple associated clocks|continuous clock")
  # two clocks in one partition
  expect_error(mtkcompile(System(list(D(x) ~ -x + Hold(y), y ~ (k1 - 1)(y) + (k2 - 1)(y)), t,
                                 name = "bad2")), "multiple associated clocks")
})

test_that("FunctionMap runs discrete callbacks; hybrid systems take shifted initial values", {
  p <- DiscreteProblem(function(u, p, t) u + p, 0, c(0, 3), 1)
  s <- solve(p, FunctionMap(), callback = PresetTimeCallback(1, function(i) i$u[] <- 100))
  expect_equal(s[1], c(0, 100, 101, 102))
  t <- t_nounits; D <- D_nounits
  z <- Variable("z", t, default = 0); w <- Variable("w", t, default = 0)
  k <- ShiftIndex(Clock(1))
  hs <- mtkcompile(System(list(D(z) ~ Hold(w), k(w) ~ (k - 1)(w) + (k - 2)(w)), t, name = "hs"))
  sh <- solve(ODEProblem(hs, list((k - 1)(w) ~ 1, (k - 2)(w) ~ 0), c(0, 3.5)))
  expect_equal(sh(c(0.5, 1.5, 2.5, 3.2), idxs = w), c(1, 2, 3, 5))
  expect_equal(tail(sh[z], 1), 1 + 2 + 3 + 5 * 0.5, tolerance = 1e-8)
  # a sampled variable indexed with a zero shift takes the clock of k
  x <- Variable("x", t); y <- Variable("y", t)
  k2 <- ShiftIndex(Clock(0.5))
  sm <- mtkcompile(System(list(D(x) ~ -x + Hold(y), k2(y) ~ -Sample(x)), t, name = "sampled"))
  so <- solve(ODEProblem(sm, list(x ~ 1, y ~ 0), c(0, 1)), reltol = 1e-10, abstol = 1e-12)
  # y = -x(0) on [0, 0.5): x(0.5) = -1 + 2 exp(-0.5)
  x5 <- -1 + 2 * exp(-0.5)
  expect_equal(so(c(0.25, 0.75), idxs = y), c(-1, -x5), tolerance = 1e-8)
})

test_that("clock inference partitions the SampledData tutorial as ModelingToolkit does", {
  # partitions printed by MTKTearing.infer_clocks! (ModelingToolkit 11.47) for
  # docs/src/tutorials/SampledData.md; MTK 11 itself cannot compile the model
  t <- t_nounits; D <- D_nounits
  dt <- 0.5
  k <- ShiftIndex(Clock(dt))
  r <- Variable("s2d_r", t)
  px <- Variable("s2d_px", t); pu <- Variable("s2d_pu", t); py <- Variable("s2d_py", t)
  fx <- Variable("s2d_fx", t); fu <- Variable("s2d_fu", t); fy <- Variable("s2d_fy", t)
  cy <- Variable("s2d_cy", t); cr <- Variable("s2d_cr", t); cud <- Variable("s2d_cud", t)
  cyd <- Variable("s2d_cyd", t); kp <- Parameter("s2d_kp", default = 1)
  a <- exp(-dt)
  eqs <- .as_equations(list(r ~ sin(t), fu ~ r, fy ~ cr, Hold(cud) ~ pu, py ~ cy,
                            k(fx) ~ a * (k - 1)(fx) + (1 - a) * k(fu), fy ~ fx, cyd ~ Sample(cy),
                            cud ~ kp * (cr - cyd), D(px) ~ -px + pu, py ~ px))
  us <- list(r, px, pu, py, fx, fu, fy, cy, cr, cud, cyd)
  inf <- .ci_infer_clocks(eqs, .keys(us))
  expect_equal(inf$eq > 0L, c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
  discrete <- vapply(us, function(v) inf$var(v) > 0L, TRUE)
  expect_equal(vapply(us[discrete], format, ""),
               paste0("s2d_", c("r", "fx", "fu", "fy", "cr", "cud", "cyd"), "(t)"))
  expect_equal(.clock_of(unique(inf$eq[inf$eq > 0L]))$dt, dt)
})

test_that("a two state machine switches on ticks with immediate and delayed transitions", {
  t <- t_nounits; D <- D_nounits
  k <- ShiftIndex(Clock(1))
  x <- GlobalScope(Variable("s2m_x", t, default = 0))
  n <- Variable("s2m_n", t, default = 0)
  up <- System(list(k(x) ~ (k - 1)(x) + 1, k(n) ~ (k - 1)(n) + entry()), t, list(x, n), list(), name = "up")
  down <- System(list(k(x) ~ (k - 1)(x) - 1), t, list(x), list(), name = "down")
  machine <- function(reset, immediate) System(list(
    initial_state(up), transition(up, down, ticksInState() >= 3, reset = reset),
    transition(down, up, ticksInState() >= 2, immediate = immediate, reset = reset)), t,
    systems = list(up, down), name = "sm")
  expect_match(format(equations(machine(FALSE, FALSE))[[2L]]), "up .* down if \\(ticksInState\\(\\) >= 3\\)")
  z <- Variable("s2m_z", t, default = 0)
  # delayed return, no reset: x = 1 2 | 1 0 | 1 2 | 1 0 | 1 2
  cs <- mtkcompile(System(list(D(z) ~ Hold(x)), t, systems = list(machine(FALSE, FALSE)), name = "s2m"))
  expect_equal(vapply(unknowns(cs), format, ""), "s2m_z(t)")
  sol <- solve(ODEProblem(cs, list(), c(0, 10)))
  expect_equal(unname(sol(seq(0.5, 9.5, by = 1), idxs = x)), c(1, 2, 1, 0, 1, 2, 1, 0, 1, 2))
  expect_equal(sol(10, idxs = z), 11, tolerance = 1e-8)
  # entry() counts the entries of state up, at t = 0, 4 and 8
  expect_equal(unname(sol(c(0.5, 3.5, 4.5, 8.5), idxs = cs$sm$up$s2m_n)), c(1, 1, 2, 3))
  # immediate transitions that reset the target state: 1 2 | -1 | 1 2 | -1 ...
  cs <- mtkcompile(System(list(D(z) ~ Hold(x)), t, systems = list(machine(TRUE, TRUE)), name = "s2m2"))
  sol <- solve(ODEProblem(cs, list(), c(0, 6)))
  expect_equal(unname(sol(seq(0.5, 5.5, by = 1), idxs = x)), c(1, 2, -1, 1, 2, -1))
  # a condition on time, and entry(), which holds at the first tick
  b <- System(list(k(x) ~ (k - 1)(x) + 10), t, list(x), list(), name = "b")
  for (cond in list(t >= 2.5, entry())) {
    sm <- System(list(initial_state(up), transition(up, b, cond)), t, systems = list(up, b), name = "smt")
    cs <- mtkcompile(System(list(D(z) ~ Hold(x)), t, systems = list(sm), name = "s2m3"))
    sol <- solve(ODEProblem(cs, list(), c(0, 5)))
    got <- unname(sol(seq(0.5, 4.5, by = 1), idxs = x))
    if (inherits(cond, "sym") && .type(cond) == "GreaterThan") expect_equal(got, c(1, 2, 3, 10, 20))
    else expect_equal(got, c(10, 20, 30, 40, 50))
  }
})
