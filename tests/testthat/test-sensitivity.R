# Sensitivity analysis (R/sensitivity-*.R), after SciMLSensitivity.jl.
# References: analytic derivatives (exponential decay, a bouncing ball, a
# steady state, initial values bound to parameters), central finite
# differences of solves at tight tolerances (Lotka Volterra, events that
# change parameters, a mass matrix DAE, an SDE path with a fixed seed), and
# the nilsas_min example of the Julia test suite.

sensitivity_fd <- function(F, th, h = 1e-6) {
  vapply(seq_along(th), function(j) {
    a <- th; b <- th
    a[j] <- a[j] + h; b[j] <- b[j] - h
    (F(a) - F(b)) / (2 * h)
  }, 0)
}

test_that("forward and adjoint sensitivities of exponential decay match the analytic derivatives", {
  prob <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
  fp <- ODEForwardSensitivityProblem(prob)
  fs <- solve(fp, Tsit5(), reltol = 1e-10, abstol = 1e-12)
  s <- extract_local_sensitivities(fs, length(fs$t))
  expect_equal(s$u[[1]], exp(-2), tolerance = 1e-9)
  expect_equal(s$du[[1]][[1]], -exp(-2), tolerance = 1e-8)
  all_s <- extract_local_sensitivities(fs)
  expect_equal(nrow(all_s$u), length(fs$t))
  expect_equal(all_s$du[[1]][, 1], -fs$t * exp(-2 * fs$t), tolerance = 1e-7)
  fd <- ODEForwardSensitivityProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2,
                                     sensealg = ForwardDiffSensitivity(), du0 = cbind(1, 0), dp = cbind(0, 1))
  m <- extract_local_sensitivities(solve(fd, Tsit5(), reltol = 1e-10, abstol = 1e-12), t = 1, asmatrix = TRUE)
  expect_equal(as.numeric(m$du), c(exp(-2), -exp(-2)), tolerance = 1e-7)

  sol <- solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-12)
  # G = u(1) and G = int_0^1 u^2 dt
  dG2 <- c(2 * (1 - exp(-4)) / 4, (8 * exp(-4) - 2 * (1 - exp(-4))) / 16)
  for (sa in list(InterpolatingAdjoint(), BacksolveAdjoint(), QuadratureAdjoint(abstol = 1e-10, reltol = 1e-10),
                  GaussAdjoint(), GaussKronrodAdjoint(), InterpolatingAdjoint(checkpointing = TRUE),
                  GaussAdjoint(checkpointing = TRUE), InterpolatingAdjoint(autojacvec = FALSE))) {
    r <- adjoint_sensitivities(sol, Tsit5(), t = 1, dgdu_discrete = function(u, p, t, i) 1,
                               abstol = 1e-10, reltol = 1e-10, sensealg = sa)
    expect_equal(c(r$du0, r$dp), c(exp(-2), -exp(-2)), tolerance = 1e-8, label = sa$name)
    r <- adjoint_sensitivities(sol, Tsit5(), g = function(u, p, t) u^2, abstol = 1e-10, reltol = 1e-10,
                               sensealg = sa)
    expect_equal(c(r$du0, r$dp), dG2, tolerance = 1e-7, label = sa$name)
  }
  # the adjoint problem itself
  ap <- ODEAdjointProblem(sol, GaussAdjoint(), Tsit5(), t = 1, dgdu_discrete = function(u, p, t, i) 1)
  expect_s3_class(ap, "ODEProblem")
  r <- ap$finish(solve(ap, Tsit5(), abstol = 1e-10, reltol = 1e-10))
  expect_equal(r$dp, -exp(-2), tolerance = 1e-8)
})

test_that("Lotka Volterra loss gradients agree for forward, adjoint and finite differences", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("sensitivity_lv_x", t); y <- Variable("sensitivity_lv_y", t)
  a <- Parameter("sensitivity_lv_a"); b <- Parameter("sensitivity_lv_b")
  c3 <- Parameter("sensitivity_lv_c"); d <- Parameter("sensitivity_lv_d")
  sys <- mtkcompile(System(list(D(x) ~ a * x - b * x * y, D(y) ~ -c3 * y + d * x * y), t, name = "sensitivity_lv"))
  prob <- ODEProblem(sys, list(x ~ 1, y ~ 1, a ~ 1.5, b ~ 1, c3 ~ 3, d ~ 1), c(0, 3))
  ts <- seq(0, 3, by = 0.5)
  loss <- function(sol) sum((sol[x] - 1)^2) + sum(sol[y])
  F <- function(th) {
    pr <- remake(prob, u0 = th[1:2], p = th[3:6])
    loss(solve(pr, Tsit5(), reltol = 1e-12, abstol = 1e-12, saveat = ts))
  }
  ref <- sensitivity_fd(F, c(prob$u0, prob$p))
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint(), BacksolveAdjoint(),
                  QuadratureAdjoint(abstol = 1e-9, reltol = 1e-9), GaussAdjoint(), GaussKronrodAdjoint(),
                  InterpolatingAdjoint(checkpointing = TRUE), GaussAdjoint(autojacvec = FALSE))) {
    gr <- solve_gradient(loss, prob, Tsit5(), saveat = ts, reltol = 1e-9, abstol = 1e-10, sensealg = sa)
    expect_equal(unname(c(gr$u0, gr$p)), ref, tolerance = 1e-6, label = sa$name)
  }
  expect_equal(names(gr$p), names(prob$p))
  # the symbolic vector Jacobian product equals J' l and fp' l
  fns <- modelcompR:::.sens_fns(prob)
  l <- c(0.3, -0.7)
  v <- fns$vjp(l, prob$u0, prob$p, 0)
  expect_equal(v$u, as.numeric(crossprod(fns$jac(prob$u0, prob$p, 0), l)))
  expect_equal(v$p, as.numeric(crossprod(fns$pjac(prob$u0, prob$p, 0), l)))
})

test_that("gradients through a bouncing ball event match the analytic derivative", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("sensitivity_ball_x", t); v <- Variable("sensitivity_ball_v", t)
  g <- Parameter("sensitivity_ball_g"); e <- Parameter("sensitivity_ball_e")
  sys <- mtkcompile(System(list(D(x) ~ v, D(v) ~ -g), t, name = "sensitivity_ball",
                           continuous_events = list(list(x ~ 0), list(v ~ -e * Pre(v)))))
  prob <- ODEProblem(sys, list(x ~ 1, v ~ 0, g ~ 9.81, e ~ 0.9), c(0, 0.8))
  height <- function(th, T = 0.8) {
    tau <- (th[2] + sqrt(th[2]^2 + 2 * th[3] * th[1])) / th[3]
    vp <- -th[4] * (th[2] - th[3] * tau)
    vp * (T - tau) - th[3] / 2 * (T - tau)^2
  }
  iu <- match(c("sensitivity_ball_x(t)", "sensitivity_ball_v(t)"), names(prob$u0))
  ip <- match(c("sensitivity_ball_g", "sensitivity_ball_e"), names(prob$p))
  ref <- sensitivity_fd(height, c(1, 0, 9.81, 0.9))
  loss <- function(sol) tail(sol[x], 1)
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint(), GaussAdjoint(),
                  QuadratureAdjoint(abstol = 1e-10, reltol = 1e-10))) {
    gr <- solve_gradient(loss, prob, Tsit5(), saveat = 0.8, reltol = 1e-10, abstol = 1e-12, sensealg = sa)
    expect_equal(unname(c(gr$u0[iu], gr$p[ip])), ref, tolerance = 1e-7, label = sa$name)
  }
  # the direct interface needs the events recorded by solve(..., sensealg = )
  sol <- solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-12, sensealg = InterpolatingAdjoint())
  dg <- numeric(2); dg[iu[1]] <- 1
  r <- adjoint_sensitivities(sol, Tsit5(), t = 0.8, dgdu_discrete = function(u, p, t, i) dg,
                             abstol = 1e-10, reltol = 1e-10)
  expect_equal(unname(c(r$du0[iu], r$dp[ip])), ref, tolerance = 1e-7)
  plain <- solve(prob, Tsit5())
  expect_error(adjoint_sensitivities(plain, Tsit5(), t = 0.8, dgdu_discrete = function(u, p, t, i) dg),
               "sensealg")
  # forward sensitivities across the event
  fs <- solve(ODEForwardSensitivityProblem(prob, sensealg = ForwardDiffSensitivity()), Tsit5(),
              reltol = 1e-10, abstol = 1e-12)
  S <- extract_local_sensitivities(fs, length(fs$t), asmatrix = TRUE)$du
  expect_equal(unname(S[iu[1], ip]), ref[3:4], tolerance = 1e-7)
})

test_that("events that change parameters and terminate are differentiated", {
  f <- function(u, p, t) c(u[2], -p[1] + p[3] * u[1])
  cb <- ContinuousCallback(function(u, t, integ) u[1] - 0.2 * t, function(integ) {
    integ$u[2] <- -0.8 * integ$u[2]
    integ$p[1] <- integ$p[1] * integ$p[2]
    integ$p[3] <- 0.1 * integ$u[2]
  })
  prob <- ODEProblem(f, c(1, 0.3), c(0, 1.2), c(9.81, 0.7, 0.5), callback = cb)
  ts <- seq(0, 1.2, by = 0.2)
  loss <- function(sol) {
    U <- as.matrix(sol)
    sum(U[, 1]^2) + U[nrow(U), 2]
  }
  F <- function(th) loss(solve(remake(prob, u0 = th[1:2], p = th[3:5]), Tsit5(), reltol = 1e-12,
                               abstol = 1e-12, saveat = ts))
  ref <- sensitivity_fd(F, c(1, 0.3, 9.81, 0.7, 0.5))
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint(), BacksolveAdjoint(), GaussAdjoint())) {
    gr <- solve_gradient(loss, prob, Tsit5(), saveat = ts, reltol = 1e-10, abstol = 1e-12, sensealg = sa)
    expect_equal(unname(c(gr$u0, gr$p)), ref, tolerance = 1e-6, label = sa$name)
  }
  expect_error(solve_gradient(loss, prob, Tsit5(), saveat = ts, sensealg = QuadratureAdjoint()),
               "depend on the state")
  # an event that terminates: the last saved point moves with the event time
  cbt <- ContinuousCallback(function(u, t, integ) u[1], function(integ) terminate(integ))
  probt <- ODEProblem(function(u, p, t) c(u[2], -p[1]), c(1, 0), c(0, 2), 9.81, callback = cbt)
  losst <- function(sol) {
    U <- as.matrix(sol)
    U[nrow(U), 2] + sum(U[, 1])
  }
  tt <- seq(0, 2, by = 0.1)
  Ft <- function(th) losst(solve(remake(probt, u0 = th[1:2], p = th[3]), Tsit5(), reltol = 1e-12,
                                 abstol = 1e-12, saveat = tt))
  ref <- sensitivity_fd(Ft, c(1, 0, 9.81))
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint(), QuadratureAdjoint(abstol = 1e-10, reltol = 1e-10))) {
    gr <- solve_gradient(losst, probt, Tsit5(), saveat = tt, reltol = 1e-10, abstol = 1e-12, sensealg = sa)
    expect_equal(unname(c(gr$u0, gr$p)), ref, tolerance = 1e-6, label = sa$name)
  }
})

test_that("discrete events that change states and parameters are differentiated", {
  f <- function(u, p, t) c(-p[1] * u[1] + u[2], -p[2] * u[2])
  cb <- CallbackSet(PresetTimeCallback(0.35, function(integ) integ$u[1] <- integ$u[1] * integ$p[3] + integ$u[2]),
                    PresetTimeCallback(0.7, function(integ) integ$p[1] <- integ$p[1] + 0.5 * integ$u[1]))
  prob <- ODEProblem(f, c(1, 0.5), c(0, 1), c(1.2, 0.8, 0.6), callback = cb)
  ts <- seq(0, 1, by = 0.25)
  loss <- function(sol) sum(as.matrix(sol)^2)
  F <- function(th) loss(solve(remake(prob, u0 = th[1:2], p = th[3:5]), Tsit5(), reltol = 1e-12,
                               abstol = 1e-12, saveat = ts))
  ref <- sensitivity_fd(F, c(1, 0.5, 1.2, 0.8, 0.6))
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint(), GaussAdjoint(), BacksolveAdjoint())) {
    gr <- solve_gradient(loss, prob, Tsit5(), saveat = ts, reltol = 1e-10, abstol = 1e-12, sensealg = sa)
    expect_equal(unname(c(gr$u0, gr$p)), ref, tolerance = 1e-7, label = sa$name)
  }
  # without saveat the event time is saved twice (left and right limits),
  # and both samples move with the event: forward and adjoint agree
  fb <- function(u, p, t) c(u[2], -p[1])
  cbb <- ContinuousCallback(function(u, t, integ) u[1], function(integ) integ$u[2] <- -integ$p[2] * integ$u[2])
  probb <- ODEProblem(fb, c(1, 0), c(0, 0.8), c(9.81, 0.9), callback = cbb)
  lb <- function(sol) {
    U <- as.matrix(sol)
    sum(U[, 1]^2 + U[, 2])
  }
  gf <- solve_gradient(lb, probb, RK4(), dt = 0.05, sensealg = ForwardDiffSensitivity())
  expect_equal(sum(duplicated(gf$sol$t)), 1L)
  for (sa in list(InterpolatingAdjoint(), GaussAdjoint(), BacksolveAdjoint())) {
    ga <- solve_gradient(lb, probb, RK4(), dt = 0.05, sensealg = sa)
    expect_equal(c(ga$u0, ga$p), c(gf$u0, gf$p), tolerance = 1e-8, label = sa$name)
  }
})

test_that("solve() callbacks, terminating running costs, unsaved end states and fixed save times", {
  last <- function(sol) tail(as.matrix(sol)[, 1], 1)
  # a callback given to solve_gradient() is differentiated: u(1) = 2 u0
  prob <- ODEProblem(function(u, p, t) 0 * u * p[1], 1, c(0, 1), 1)
  cb <- PresetTimeCallback(0.5, function(integ) integ$u[1] <- 2 * integ$u[1])
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint()))
    expect_equal(unname(solve_gradient(last, prob, Tsit5(), saveat = 1, callback = cb, sensealg = sa)$u0), 2,
                 label = sa$name)
  # G = int_0^tau 1 dt with u' = p terminated at u = 1: G = 1 / p
  probt <- ODEProblem(function(u, p, t) p[1], 0, c(0, 5), 2,
                      callback = ContinuousCallback(function(u, t, integ) u[1] - 1, function(integ) terminate(integ)))
  sol <- solve(probt, Tsit5(), reltol = 1e-10, abstol = 1e-12, sensealg = InterpolatingAdjoint())
  for (sa in list(InterpolatingAdjoint(), GaussAdjoint(), BacksolveAdjoint()))
    expect_equal(adjoint_sensitivities(sol, Tsit5(), g = function(u, p, t) 1, sensealg = sa, abstol = 1e-10,
                                       reltol = 1e-10)$dp, -1 / 4, tolerance = 1e-8, label = sa$name)
  # BacksolveAdjoint starts from the state at the final time, which is not saved
  pe <- ODEProblem(function(u, p, t) p[1] * u, 1, c(0, 1), 1)
  expect_equal(unname(solve_gradient(last, pe, Tsit5(), saveat = c(0, 0.5), reltol = 1e-10, abstol = 1e-12,
                                     sensealg = BacksolveAdjoint())$p), 0.5 * exp(0.5), tolerance = 1e-7)
  # a requested save time that coincides with an event does not move with it
  pf <- ODEProblem(function(u, p, t) 1 + 0 * p[1], 0, c(0, 1), 0.5,
                   callback = ContinuousCallback(function(u, t, integ) u[1] - integ$p[1], function(integ) NULL))
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint()))
    expect_equal(unname(solve_gradient(function(sol) as.matrix(sol)[2, 1], pf, Tsit5(), saveat = c(0, 0.5, 1),
                                       sensealg = sa)$p), 0, tolerance = 1e-8, label = sa$name)
  # decreasing time spans: u(0) = u0 exp(p) for u' = -p u from t = 1 to 0
  pd <- ODEProblem(function(u, p, t) -p[1] * u, 1.5, c(1, 0), 2)
  for (sa in list(ForwardDiffSensitivity(), ForwardSensitivity(), InterpolatingAdjoint(), BacksolveAdjoint(),
                  GaussAdjoint())) {
    gd <- solve_gradient(last, pd, Tsit5(), saveat = c(1, 0), reltol = 1e-10, abstol = 1e-12, sensealg = sa)
    expect_equal(unname(c(gd$u0, gd$p)), c(exp(2), 1.5 * exp(2)), tolerance = 1e-6, label = sa$name)
  }
  pdev <- ODEProblem(function(u, p, t) -p[1] * u, 1.5, c(1, 0), 2,
                     callback = PresetTimeCallback(0.5, function(integ) integ$u <- integ$u + 1))
  expect_error(solve_gradient(last, pdev, Tsit5(), sensealg = InterpolatingAdjoint()), "decreasing")
})

test_that("mass matrix DAE gradients agree with finite differences", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("sensitivity_dae_x", t); y <- Variable("sensitivity_dae_y", t)
  k <- Parameter("sensitivity_dae_k"); a <- Parameter("sensitivity_dae_a")
  sys <- mtkcompile(System(list(D(x) ~ -k * y, 0 ~ y^3 + y - a * x), t, name = "sensitivity_dae"))
  expect_true(any(rowSums(abs(calculate_massmatrix(sys))) == 0))
  prob <- ODEProblem(sys, list(x ~ 1, k ~ 0.7, a ~ 0.5), c(0, 2), jac = TRUE, guesses = list(y ~ 0.5))
  ix <- match("sensitivity_dae_x(t)", names(prob$u0))
  ts <- seq(0, 2, by = 0.5)
  loss <- function(sol) sum(sol[x]^2) + 3 * sum(sol[y])
  F <- function(th) {
    u0 <- prob$u0
    u0[ix] <- th[1]
    loss(solve(remake(prob, u0 = u0, p = th[2:3]), Rodas5P(), reltol = 1e-11, abstol = 1e-12, saveat = ts))
  }
  ref <- sensitivity_fd(F, c(1, prob$p), h = 1e-5)
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint(), GaussAdjoint(),
                  QuadratureAdjoint(abstol = 1e-9, reltol = 1e-9))) {
    gr <- solve_gradient(loss, prob, Rodas5P(), saveat = ts, reltol = 1e-9, abstol = 1e-10, sensealg = sa)
    expect_equal(unname(c(gr$u0[ix], gr$p)), ref, tolerance = 1e-5, label = sa$name)
    # the algebraic initial value is recomputed, so it has no influence
    expect_equal(unname(gr$u0[-ix]), 0)
  }
})

test_that("steady state and nonlinear solutions are differentiated by the implicit function theorem", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("sensitivity_ss_x", t)
  k <- Parameter("sensitivity_ss_k"); s <- Parameter("sensitivity_ss_s")
  sys <- mtkcompile(System(list(D(x) ~ s - k * x^2), t, name = "sensitivity_ss"))
  prob <- SteadyStateProblem(sys, list(x ~ 1, k ~ 2, s ~ 3))
  sol <- solve(prob, SSRootfind())
  # x*^2 = s / k
  r <- adjoint_sensitivities(sol, NULL, sensealg = SteadyStateAdjoint(), g = function(u, p) u[1]^2)
  expect_equal(unname(r$dp[c("sensitivity_ss_k", "sensitivity_ss_s")]), c(-3 / 4, 1 / 2), tolerance = 1e-8)
  r <- adjoint_sensitivities(sol, NULL, dgdu = function(u, p) 2 * u)
  expect_equal(unname(r$dp[c("sensitivity_ss_k", "sensitivity_ss_s")]), c(-3 / 4, 1 / 2), tolerance = 1e-8)
  gr <- solve_gradient(function(sol) sol$u[[1]]^2, prob, SSRootfind())
  expect_equal(unname(gr$p[c("sensitivity_ss_k", "sensitivity_ss_s")]), c(-3 / 4, 1 / 2), tolerance = 1e-8)
  # a nonlinear system: x^2 = a, dx/da = 1 / (2 sqrt(a))
  xn <- Variable("sensitivity_nl_x"); an <- Parameter("sensitivity_nl_a")
  nsys <- complete(System(list(0 ~ xn^2 - an), name = "sensitivity_nl"))
  gr <- solve_gradient(function(sol) sol$u[[1]], NonlinearProblem(nsys, list(xn ~ 1, an ~ 2)), NewtonRaphson())
  expect_equal(unname(gr$p), 1 / (2 * sqrt(2)), tolerance = 1e-8)
  nprob <- NonlinearProblem(function(u, p) c(u[1]^2 - p[1], u[2] - p[2] * u[1]), c(1, 1), p = c(2, 3))
  gr <- solve_gradient(function(sol) sum(sol$u), nprob, NewtonRaphson())
  expect_equal(unname(gr$p), c(4 / (2 * sqrt(2)), sqrt(2)), tolerance = 1e-7)
})

test_that("initial values that depend on parameters through initialization are differentiated", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("sensitivity_init_x", t, default = Symbol("sensitivity_init_a") * 2)
  k <- Parameter("sensitivity_init_k", default = 0.5)
  a <- Parameter("sensitivity_init_a", default = 2)
  sys <- mtkcompile(System(list(D(x) ~ -k * x + 0 * a), t, name = "sensitivity_init"))
  prob <- ODEProblem(sys, list(), c(0, 1))
  loss <- function(sol) tail(sol[x], 1)
  # x(1) = 2 a exp(-k)
  ref <- c(sensitivity_init_k = -4 * exp(-0.5), sensitivity_init_a = 2 * exp(-0.5))
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint())) {
    gr <- solve_gradient(loss, prob, Tsit5(), saveat = 1, reltol = 1e-10, abstol = 1e-12, sensealg = sa)
    expect_equal(gr$p[names(ref)], ref, tolerance = 1e-7, label = sa$name)
  }
  fs <- solve(ODEForwardSensitivityProblem(prob), Tsit5(), reltol = 1e-10, abstol = 1e-12)
  expect_equal(extract_local_sensitivities(fs, 1L)$du$sensitivity_init_a, c(`sensitivity_init_x(t)` = 2))
})

test_that("SDE paths and second order sensitivities", {
  prob <- SDEProblem(function(u, p, t) p[1] * u, function(u, p, t) p[2] * u, 1, c(0, 1), c(1.5, 0.3))
  loss <- function(sol) tail(as.matrix(sol)[, 1], 1)
  gr <- solve_gradient(loss, prob, EM(), dt = 0.02, seed = 3)
  F <- function(th) loss(solve(remake(prob, u0 = th[1], p = th[2:3]), EM(), dt = 0.02, seed = 3))
  expect_equal(c(gr$u0, gr$p), sensitivity_fd(F, c(1, 1.5, 0.3)), tolerance = 1e-7)
  expect_error(solve_gradient(loss, prob, EM(), dt = 0.02, sensealg = GaussAdjoint()), "pathwise")
  prob2 <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
  loss2 <- function(sol) sol$u[[length(sol$t)]]
  H <- second_order_sensitivities(loss2, prob2, Tsit5(), saveat = 1, reltol = 1e-12, abstol = 1e-12)
  expect_equal(as.numeric(H), exp(-2), tolerance = 1e-6)
  Hv <- second_order_sensitivity_product(loss2, 2, prob2, Tsit5(), saveat = 1, reltol = 1e-12, abstol = 1e-12)
  expect_equal(unname(Hv), 2 * exp(-2), tolerance = 1e-6)
})

test_that("shadowing methods give the sensitivity of the Lorenz attractor average", {
  # d<z>/d(rho) is close to 1 (Lea, Allen and Haine 2000; the Julia tests)
  lorenz <- function(u, p, t) c(p[1] * (u[2] - u[1]), u[1] * (p[2] - u[3]) - u[2], u[1] * u[2] - p[3] * u[3])
  prob <- ODEProblem(lorenz, c(-12.808828477348253, -18.197181420445414, 25.88373061170027), c(0, 10),
                     c(10, 28, 8 / 3))
  g <- function(u, p, t) u[3]
  sol <- solve(prob, Tsit5(), saveat = 0.02, reltol = 1e-8, abstol = 1e-8)
  lss <- shadow_forward(ForwardLSSProblem(sol, ForwardLSS(g = g)))
  expect_equal(lss[2], 1, tolerance = 0.1)
  expect_equal(shadow_adjoint(AdjointLSSProblem(sol, AdjointLSS(g = g))), lss, tolerance = 1e-10)
  expect_equal(shadow_forward(ForwardLSSProblem(sol, ForwardLSS(g = g),
                                                dgdu_continuous = function(u, p, t) c(0, 0, 1))),
               lss, tolerance = 1e-10)
  win <- shadow_forward(ForwardLSSProblem(sol, ForwardLSS(LSSregularizer = Cos2Windowing(), g = g)))
  expect_equal(shadow_forward(ForwardLSSProblem(sol, ForwardLSS(LSSregularizer = Cos2Windowing()),
                                                dgdu_continuous = function(u, p, t) c(0, 0, 1))),
               win, tolerance = 1e-10)
  expect_equal(win[2], 1, tolerance = 0.3)
  nilss <- shadow_forward(NILSSProblem(prob, NILSS(10, 51, g = g, seed = 1)), Tsit5(), reltol = 1e-6, abstol = 1e-6)
  expect_equal(nilss[1:2], lss[1:2], tolerance = 0.05)
  dense <- solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
  nilsas <- shadow_adjoint(NILSASProblem(dense, NILSAS(10, 51, M = 2, g = g, seed = 1)), Tsit5(),
                           reltol = 1e-6, abstol = 1e-6)
  expect_equal(nilsas[1:2], lss[1:2], tolerance = 0.05)
  expect_equal(nilsas[3], lss[3], tolerance = 0.2)
  # Julia's default adjoint inside NILSAS: the state is carried backward
  nilsas_b <- shadow_adjoint(NILSASProblem(dense, NILSAS(10, 51, M = 2, g = g, seed = 1,
                                                         adjoint_sensealg = BacksolveAdjoint())),
                             Tsit5(), reltol = 1e-6, abstol = 1e-6)
  expect_equal(nilsas_b, nilsas, tolerance = 1e-3)
  expect_error(NILSAS(10, 51, M = 2, adjoint_sensealg = ForwardSensitivity()), "continuous adjoint")
  # nilsas_min on the example of the Julia test suite
  C <- array(c(1, 0, 0, 1, 4, 0, 0, 1), c(2, 2, 2))
  R <- array(c(Inf, Inf, Inf, Inf, 1, 0, 1, 2), c(2, 2, 2))
  a <- modelcompR:::.nilsas_min(C, R, matrix(c(Inf, Inf, 0, 1), 2), matrix(c(1, 0, 1, 4), 2),
                                matrix(c(1, 1, 3, 1), 2), c(1, 2))
  expect_equal(a, rbind(c(-1, 0), c(-1, -1)), tolerance = 1e-12)
  expect_error(ForwardLSSProblem(sol, ForwardLSS()), "Time dilation needs explicit knowledge of g")
})

test_that("the loss may not interpolate, and the adjoint needs a dense forward solution", {
  prob <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
  expect_error(solve_gradient(function(sol) sol(0.5), prob, Tsit5(), saveat = 0.5, sensealg = GaussAdjoint()),
               "saveat")
  sol <- solve(prob, Tsit5(), saveat = 0.5)
  expect_error(adjoint_sensitivities(sol, Tsit5(), t = 1, dgdu_discrete = function(u, p, t, i) 1), "dense")
  r <- adjoint_sensitivities(sol, Tsit5(), t = 1, dgdu_discrete = function(u, p, t, i) 1,
                             sensealg = InterpolatingAdjoint(checkpointing = TRUE), abstol = 1e-8, reltol = 1e-8)
  expect_equal(r$dp, -exp(-2), tolerance = 1e-4)
  expect_output(print(InterpolatingAdjoint(checkpointing = TRUE)), "checkpointing = TRUE")
  expect_error(solve(prob, Tsit5(), sensealg = "adjoint"), "sensitivity algorithm")
})
