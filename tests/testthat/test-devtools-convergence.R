# Convergence tools of DiffEqDevTools.jl (src/convergence.jl,
# src/test_solution.jl). Reference values from tools/DiffEqDevTools-probe.jl.

dv_ode_prob <- function(analytic = TRUE) {
  f <- function(u, p, t) c(-u[1] * u[1], u[1] * u[2])
  an <- function(u0, p, t) c(1 / (1 + t), 2 * (1 + t))
  if (analytic) ODEProblem(ODEFunction(f, analytic = an), c(1, 2), c(0, 1)) else ODEProblem(f, c(1, 2), c(0, 1))
}

# relative comparison of every entry (expect_equal compares absolute
# differences when the values are below the tolerance)
dv_check_sim <- function(sim, ref, tol) {
  for (k in names(ref)) expect_equal(sim$errors[[k]] / ref[[k]], rep(1, length(ref[[k]])), tolerance = tol, info = k)
}

test_that("test_convergence of ODEs matches DiffEqDevTools", {
  prob <- dv_ode_prob()
  dts <- 1 / 2^(3:6)
  sim <- test_convergence(dts, prob, Euler(), dense_errors = TRUE)
  dv_check_sim(sim, list(
    final = c(0.04773824454478556, 0.02419923396292173, 0.012194281310070193, 0.006122276238434227),
    l2 = c(0.029897093511641503, 0.01485612130512841, 0.007411988281960016, 0.0037027329806756408),
    l_inf = c(0.07228797726857117, 0.03720418642663503, 0.018884429256151414, 0.00951495988860751),
    L2 = c(0.028610962082216473, 0.01455207778753441, 0.007349468014173333, 0.003694472442119961),
    L_inf = c(0.07228797726857117, 0.03720418642663503, 0.018884429256151414, 0.00951495988860751)), 1e-12)
  expect_equal(sim$Oest$final, 0.9876684899043942, tolerance = 1e-12)
  expect_equal(sim$Oest$L2, 0.9843759129041936, tolerance = 1e-12)
  sim <- test_convergence(dts, prob, RK4(), dense_errors = TRUE)
  dv_check_sim(sim, list(
    final = c(1.4397699634405292e-5, 8.515960957522495e-7, 5.17468077276817e-8, 3.1885115658880636e-9),
    l2 = c(1.3615067246749641e-5, 8.023536739292566e-7, 4.86689857368862e-8, 2.996294745592308e-9),
    L2 = c(1.3655078203278797e-5, 8.093100548922644e-7, 4.92097654920547e-8, 3.032608289983433e-9)), 1e-9)
  expect_equal(sim$Oest$l_inf, 4.048608296293494, tolerance = 1e-9)
  sim <- test_convergence(dts, prob, Tsit5(), dense_errors = TRUE)
  dv_check_sim(sim, list(
    final = c(1.0021013646666077e-7, 1.337808253065731e-9, 2.02100558510665e-11, 3.396172232328354e-13),
    L_inf = c(2.0634539044106504e-7, 7.688776948455711e-9, 2.627035255997612e-10, 8.028355757971894e-12)), 1e-12)
  expect_equal(names(test_convergence(dts, prob, Euler())$errors), c("final", "l2", "l_inf"))
  expect_equal(names(test_convergence(dts, prob, Euler(), save_everystep = FALSE)$errors), "final")
  probs <- lapply(c(0.5, 1), function(T) remake(prob, tspan = c(0, T)))
  expect_equal(test_convergence(probs, c(0.5, 1), RK4(), dt = 1 / 8, adaptive = FALSE)$errors$final,
    c(1.0032435632301162e-5, 1.4397699634405292e-5), tolerance = 1e-9)
})

test_that("appxtrue matches DiffEqDevTools", {
  prob <- dv_ode_prob()
  sol <- solve(prob, Euler(), dt = 1 / 16)
  ref <- solve(prob, RK4(), dt = 1 / 256, adaptive = FALSE)
  e <- appxtrue(sol, ref)$errors
  expect_equal(e$final, 0.024199233951308685, tolerance = 1e-12)
  expect_equal(e$l2, 0.014856121296003114, tolerance = 1e-12)
  expect_equal(e$l_inf, 0.03720418640270173, tolerance = 1e-12)
  expect_equal(e$L2, 0.014552077778095205, tolerance = 1e-12)
  expect_equal(appxtrue(sol, TestSolution(ref))$errors$L2, 0.014552077778095205, tolerance = 1e-12)
  expect_equal(appxtrue(sol, TestSolution(ref$t, ref$u, ref))$errors$L2, 0.014552077778095205, tolerance = 1e-12)
  e <- appxtrue(sol, TestSolution(ref$t, ref$u))$errors
  expect_equal(names(e), "final")
  expect_equal(e$final, 0.024199233951308685, tolerance = 1e-12)
  psc <- ODEProblem(function(u, p, t) -u * u, 1, c(0, 1))
  ssc <- solve(psc, Euler(), dt = 1 / 16)
  rsc <- solve(psc, RK4(), dt = 1 / 16, adaptive = FALSE)
  e <- appxtrue(ssc, TestSolution(rsc$t, rsc$u))$errors
  expect_equal(unlist(e[c("final", "l2", "l_inf")]),
    c(final = 0.01119432750443361, l2 = 0.010365268505397391, l_inf = 0.01203355069891865), tolerance = 1e-12)
  e <- appxtrue(ssc, solve(psc, RK4(), dt = 1 / 16, adaptive = FALSE, dense = FALSE))$errors
  expect_equal(e, list(final = 0.01119432750443361), tolerance = 1e-12)
  # Julia throws a DimensionMismatch for vector states here (probe line
  # appxTSsame); R computes the errors at the common times.
  e <- appxtrue(sol, TestSolution(sol$t, solve(prob, RK4(), dt = 1 / 16, adaptive = FALSE)$u))$errors
  expect_equal(names(e), c("final", "l2", "l_inf"))
})

test_that("analyticless_test_convergence of ODEs matches DiffEqDevTools", {
  sim <- analyticless_test_convergence(1 / 2^(3:6), dv_ode_prob(FALSE), Heun(),
    list(alg = RK4(), dt = 1 / 512, adaptive = FALSE))
  dv_check_sim(sim, list(
    final = c(0.004857827674946025, 0.0011532999402482091, 0.00028129159121131897, 6.94813969890351e-5),
    l2 = c(0.004235550139540393, 0.001000799468625751, 0.00024346550985521485, 6.005721150725352e-5),
    l_inf = c(0.008649603534052641, 0.0020510743287491806, 0.0005001192823028688, 0.00012352524501135775),
    L2 = c(0.004210251542038806, 0.0009984093462758906, 0.0002432418316800878, 6.004251719688372e-5)), 1e-9)
  expect_equal(sim$Oest$final, 2.042513580620158, tolerance = 1e-9)
})

# A fixed Brownian path on a grid of step h: u1 = 1 + t^2 + W1/2,
# u2 = 2 + t + W2/4.
dv_grid_path <- function(h) {
  tg <- seq(0, 1, by = h)
  W <- lapply(seq_along(tg), function(k) {
    t <- tg[k]
    c(0.3 * t - 0.7 * t * t + 0.01 * (k %% 3), -0.2 * t + 0.5 * t * t * t - 0.02 * (k %% 5))
  })
  W[[1]] <- c(0, 0)
  list(t = tg, W = W)
}

dv_sde_prob <- function(h) {
  g <- dv_grid_path(h)
  an <- function(u0, p, t, W) c(u0[1] + t * t + 0.5 * W[1], u0[2] + t + 0.25 * W[2])
  SDEProblem(SDEFunction(function(u, p, t) c(2 * t, 1), function(u, p, t) c(0.5, 0.25), analytic = an),
    c(1, 2), c(0, 1), noise = NoiseGrid(g$t, g$W))
}

dv_sde_prob1 <- function(h) {
  g <- dv_grid_path(h)
  an <- function(u0, p, t, W) u0 + t * t + 0.5 * W
  SDEProblem(SDEFunction(function(u, p, t) 2 * t, function(u, p, t) 0.5, analytic = an), 1, c(0, 1),
    noise = NoiseGrid(g$t, vapply(g$W, `[`, 0, 1)))
}

test_that("test_convergence of SDEs matches DiffEqDevTools", {
  dts <- c(1 / 8, 1 / 16, 1 / 32)
  sprob <- dv_sde_prob(1 / 64)
  # Deliberate difference: SciMLBase takes the ith value of the noise grid
  # (sol.W.u[i]) for the ith step when the grid is finer than the steps; R
  # takes W(t[i]). The probe lines "sde at W(t)" recompute the Julia errors
  # at W(t[i]); EM is then exact up to the drift error dt / 2 at t = 1.
  sim <- suppressMessages(test_convergence(dts, sprob, EM(), EnsembleSerial(), trajectories = 2))
  expect_setequal(names(sim$errors), c("final", "l2", "l_inf", "weak_final"))
  expect_equal(sim$errors$final, c(0.0625, 0.031250000000000777, 0.015624999999999778), tolerance = 1e-12)
  expect_equal(sim$errors$l2, c(0.052601588220635827, 0.025911131174651544, 0.012857042746179097),
    tolerance = 1e-12)
  sprob1 <- dv_sde_prob1(1 / 64)
  sim <- suppressMessages(test_convergence(dts, sprob1, EM(), EnsembleSerial(), trajectories = 2,
    dense_errors = TRUE))
  expect_equal(sim$errors$L2, c(0.0712076263362812, 0.03594148368252975, 0.018389235416188355), tolerance = 1e-12)
  expect_equal(sim$errors$L_inf, c(0.125, 0.06250000000000022, 0.035263926951012214), tolerance = 1e-12)
  sim <- suppressMessages(test_convergence(dts, sprob1, EM(), EnsembleSerial(), trajectories = 2,
    weak_dense_errors = TRUE))
  expect_equal(sim$errors$weak_L2, c(0.0712076263362812, 0.03594148368252975, 0.018389235416188355),
    tolerance = 1e-12)
  eprob <- EnsembleProblem(sprob, output_func = function(sol, i) list(sol$u[[length(sol$u)]][1], FALSE))
  sim <- suppressMessages(test_convergence(dts, eprob, EM(), EnsembleSerial(), trajectories = 2, expected_value = 2))
  expect_equal(sim$errors$weak_final, c(0.31499999999999995, 0.25250000000000017, 0.2212499999999995),
    tolerance = 1e-12)
  expect_equal(sim$Oest$weak_final, 0.254837186708461, tolerance = 1e-12)
  eprob2 <- EnsembleProblem(sprob, output_func = function(sol, i) list(sol$u[[length(sol$u)]], FALSE))
  sim <- suppressMessages(test_convergence(dts, eprob2, EM(), EnsembleSerial(), trajectories = 2,
    expected_value = c(2, 3)))
  expect_equal(sim$errors$weak_l2, c(0.11999999999999988, 0.08874999999999933, 0.07312499999999966),
    tolerance = 1e-12)
  # structure only: the noise of the analyticless case is random
  sim <- analyticless_test_convergence(dts, SDEProblem(function(u, p, t) c(2 * t, 1),
    function(u, p, t) c(0.5, 0.25), c(1, 2), c(0, 1)), EM(), 1 / 64, trajectories = 2, verbose = FALSE)
  expect_setequal(names(sim$errors), c("final", "weak_final"))
})

test_that("solution errors on an aligned noise grid match DiffEqDevTools", {
  e <- .dv_solution_errors(solve(dv_sde_prob(1 / 8), EM(), dt = 1 / 8))$errors
  expect_equal(unlist(e[c("final", "l2", "l_inf")]),
    c(final = 0.06250000000000011, l2 = 0.05260158822063593, l_inf = 0.12500000000000022), tolerance = 1e-12)
  e <- .dv_solution_errors(solve(dv_sde_prob1(1 / 8), EM(), dt = 1 / 8), dense_errors = TRUE)$errors
  expect_equal(e$L2, 0.07014600447332302, tolerance = 1e-12)
  expect_equal(e$l2, 0.07438987946398817, tolerance = 1e-12)
})

test_that("test_convergence of DAEs and DDEs matches DiffEqDevTools", {
  fd <- function(du, u, p, t) c(du[1] + u[1] * u[1], u[1] * u[2] - 2)
  an <- function(du0, u0, p, t) c(1 / (1 + t), 2 * (1 + t))
  dprob <- DAEProblem(DAEFunction(fd, analytic = an), c(-1, 2), c(1, 2), c(0, 1), differential_vars = c(TRUE, FALSE))
  sim <- test_convergence(c(1 / 8, 1 / 16, 1 / 32), dprob, DImplicitEuler(), dense_errors = TRUE)
  dv_check_sim(sim, list(
    final = c(0.08850360928444562, 0.04637313806540422, 0.02376828749726939),
    l2 = c(0.07093235107973936, 0.037104784019025845, 0.019021070476860684),
    l_inf = c(0.15663094211439832, 0.08224971321510566, 0.04220646242376658)), 1e-9)
  # the dense error: DImplicitEuler interpolates with the generic Hermite
  # polynomial of OrdinaryDiffEq on its last two derivatives
  expect_equal(sim$errors$L2[1], 0.06993692498469824, tolerance = 1e-10)
  fdd <- function(u, h, p, t) -h(p, t - 1)
  andd <- function(u0, h, p, t) if (t <= 1) 1 - t else -(2 * (t - 1) - (t * t - 1) / 2)
  ddprob <- DDEProblem(DDEFunction(fdd, analytic = andd), 1, function(p, t) 1, c(0, 2), constant_lags = 1)
  e <- .dv_solution_errors(solve(ddprob, MethodOfSteps(Euler()), dt = 1 / 8))$errors
  expect_equal(unlist(e[c("final", "l2", "l_inf")]),
    c(final = 0.06249999999999989, l2 = 0.02706329386826363, l_inf = 0.06249999999999989), tolerance = 1e-12)
})
