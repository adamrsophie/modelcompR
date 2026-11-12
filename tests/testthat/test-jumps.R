# Jump processes: symbolic jumps, poissonians, SSA aggregators and jumps
# coupled to ODEs.

test_that("birth death process mean matches lam/mu (1 - exp(-mu t))", {
  t <- t_nounits
  X <- Variable("X", t)
  lam <- Parameter("lam", default = 10); mu <- Parameter("mu", default = 1)
  birth <- SymbolicMassActionJump(lam, list(), list(X ~ 1))
  death <- SymbolicMassActionJump(mu, list(X ~ 1), list(X ~ -1))
  js <- JumpSystem(list(birth, death), t, name = "bd")
  expect_equal(vapply(unknowns(js), format, ""), "X(t)")
  expect_setequal(vapply(parameters(js), format, ""), c("lam", "mu"))
  cjs <- mtkcompile(js)
  expect_length(unknowns(cjs), 1L)
  jprob <- JumpProblem(cjs, list(X ~ 0), c(0, 3), save_positions = c(FALSE, FALSE))
  expect_s3_class(jprob$prob, "DiscreteProblem")
  ts <- 0:3
  exact <- 10 * (1 - exp(-ts))
  for (agg in list(Direct(), SortingDirect())) {
    jp <- JumpProblem(cjs, list(X ~ 0), c(0, 3), aggregator = agg, save_positions = c(FALSE, FALSE))
    sim <- solve(EnsembleProblem(jp), SSAStepper(), trajectories = 400, saveat = 1, seed = 2024)
    summ <- EnsembleSummary(sim, ts = ts)
    # the count is Poisson with mean exact(t): standard error sqrt(exact / 400)
    expect_true(all(abs(summ$u[-1, 1] - exact[-1]) < 4 * sqrt(exact[-1] / 400)))
    expect_true(all(abs(summ$v[-1, 1] / exact[-1] - 1) < 0.3))
  }
  sol <- solve(JumpProblem(cjs, list(X ~ 0), c(0, 3)), seed = 1)
  # piecewise constant path: the value at any time equals the last saved value
  expect_equal(sol(1.2345)[[1]], sol[X][max(which(sol$t <= 1.2345))])
})

test_that("SIR jump model: rates, affects and the mean field limit", {
  t <- t_nounits; D <- D_nounits
  S <- Variable("S", t); I <- Variable("I", t); R <- Variable("R", t)
  beta <- Parameter("beta"); gamma <- Parameter("gamma")
  j1 <- ConstantRateJump(beta * S * I, list(S ~ Pre(S) - 1, I ~ Pre(I) + 1))
  j2 <- ConstantRateJump(gamma * I, list(I ~ Pre(I) - 1, R ~ Pre(R) + 1))
  sir <- complete(JumpSystem(list(j1, j2), t, list(S, I, R), list(beta, gamma), name = "sir"))
  N <- 500
  op <- list(S ~ N - 50, I ~ 50, R ~ 0, beta ~ 0.25 / N, gamma ~ 0.1)
  jp <- JumpProblem(sir, op, c(0, 30), save_positions = c(FALSE, FALSE))
  js <- jp$jumpset
  u <- c(100, 9, 5); p <- jp$prob$p
  expect_equal(unname(js$rates(u, p, 1)), c(0.25 / N * 900, 0.9))
  expect_equal(js$affect(1L, u, p, 1), c(99, 10, 5))
  expect_equal(js$affect(2L, u, p, 1), c(100, 8, 6))
  sim <- solve(EnsembleProblem(jp), trajectories = 60, saveat = c(10, 20, 30), seed = 7)
  summ <- EnsembleSummary(sim, ts = c(10, 20, 30))
  # mean field ODE
  ode <- mtkcompile(System(list(D(S) ~ -beta * S * I, D(I) ~ beta * S * I - gamma * I, D(R) ~ gamma * I),
                           t, name = "sirode"))
  osol <- solve(ODEProblem(ode, op, c(0, 30)), Tsit5(), reltol = 1e-8, abstol = 1e-8)
  ref <- osol(c(10, 20, 30), idxs = list(S, I, R))
  expect_lt(max(abs(summ$u - ref)) / N, 0.04)
  expect_true(all(abs(rowSums(summ$u) - N) < 1e-9))
})

test_that("variable rate jumps coupled to an ODE integrate their rate", {
  t <- t_nounits
  A <- Variable("A", t); C <- Variable("C", t); B <- Variable("B", t); k <- Parameter("k")
  vrj <- VariableRateJump(k * (sin(t) + 1), list(A ~ Pre(A) + 1, C ~ Pre(C) + 2))
  js <- complete(JumpSystem(list(vrj), t, list(A, C), list(k), name = "vjs", observed = list(B ~ C * A)))
  jp <- JumpProblem(js, list(A ~ 0, C ~ 0, k ~ 1), c(0, 6))
  expect_s3_class(jp$prob, "ODEProblem")
  sol <- solve(jp, Tsit5(), seed = 3)
  expect_equal(sol[B], sol[A] * sol[C])
  expect_equal(sol[C], 2 * sol[A])
  sim <- solve(EnsembleProblem(jp), Tsit5(), trajectories = 80, saveat = c(3, 6), seed = 11)
  summ <- EnsembleSummary(sim, ts = c(3, 6))
  mean_A <- c(3, 6) + 1 - cos(c(3, 6))
  expect_true(all(abs(summ$u[, 1] - mean_A) < 4 * sqrt(mean_A / 80)))
})

test_that("poissonians become jumps during mtkcompile", {
  t <- t_nounits; D <- D_nounits
  X <- Variable("X", t); Y <- Variable("Y", t)
  lam <- Parameter("lam", default = 2); k <- Parameter("k", default = 0.5)
  poissonians(dNj = lam, dMj = k * X)
  sys <- System(list(D(X) ~ dNj - dMj, D(Y) ~ -Y + 3 * dNj), t, name = "pp")
  expect_length(poissonians(sys), 2L)
  cs <- mtkcompile(sys)
  js <- jumps(cs)
  expect_length(js, 2L)
  expect_setequal(vapply(js, function(j) class(j)[1L], ""), c("ConstantRateJump", "VariableRateJump"))
  expect_length(equations(cs), 1L)
  crj <- Filter(function(j) inherits(j, "ConstantRateJump"), js)[[1L]]
  expect_equal(vapply(crj$affect, format, ""), c("X(t) ~ Pre(X(t)) + 1", "Y(t) ~ Pre(Y(t)) + 3"))
  jp <- JumpProblem(cs, list(X ~ 0, Y ~ 0), c(0, 2))
  expect_s3_class(jp$prob, "ODEProblem")
  # drift in unknown order: X changes only by jumps, D(Y) ~ -Y
  us <- names(jp$prob$u0)
  d <- jp$prob$f$f(c(1, 2), jp$prob$p, 0)
  expect_equal(unname(d[us == "X(t)"]), 0)
  expect_equal(unname(d[us == "Y(t)"]), -c(1, 2)[us == "Y(t)"])
  sol <- solve(jp, seed = 2)
  expect_equal(sol$retcode, "Success")
})

test_that("mass action propensities, modified unknowns and checks", {
  t <- t_nounits
  X <- Variable("X", t); Y <- Variable("Y", t); k <- Parameter("k")
  dimer <- SymbolicMassActionJump(k / 2, list(X ~ 2), list(X ~ -2, Y ~ 1))
  expect_equal(format(symcompR::expand(dimer$rate - k / 2 * X * (X - 1))), "0")
  expect_length(modified_unknowns(dimer, list(X, Y)), 2L)
  expect_error(JumpSystem(list(MassActionJump(k, list(X ~ 2), list(X ~ -2))), t, name = "bad"),
               "rescale_rates_on_update")
  expect_error(SymbolicMassActionJump(k, list(X ~ 1), list(X ~ -1), scale_rates = TRUE), "pre-scaled")
  js <- complete(JumpSystem(list(dimer), t, list(X, Y), list(k), name = "dim"))
  jp <- JumpProblem(js, list(X ~ 10, Y ~ 0, k ~ 1), c(0, 50))
  sol <- solve(jp, seed = 1)
  last <- sol$u[[length(sol$t)]]
  expect_equal(unname(last[1] + 2 * last[2]), 10)
})

test_that("jumps with brownians build a jump diffusion", {
  t <- t_nounits; D <- D_nounits
  X <- Variable("X", t); k <- Parameter("k", default = 1)
  B <- Brownian("jmpB")
  jump <- SymbolicMassActionJump(k, list(X ~ 1), list(X ~ -1))
  sys <- mtkcompile(System(list(D(X) ~ -k * X + sqrt(k) * B), t, name = "jd", jumps = list(jump)))
  expect_false(is.null(get_noise_eqs(sys)))
  jp <- JumpProblem(sys, list(X ~ 10), c(0, 1))
  expect_s3_class(jp$prob, "SDEProblem")
  expect_equal(unname(jp$prob$f$f(10, jp$prob$p, 0)), -10)
  expect_equal(unname(jp$prob$f$g(10, jp$prob$p, 0)), 1)
  sol <- solve(jp, SOSRI(), seed = 4)
  expect_equal(sol$retcode, "Success")
  expect_s3_class(sol, "RODESolution")
  expect_equal(ncol(sol$W$W), 1L)
})

test_that("raw jump problems with R functions", {
  rate <- function(u, p, t) p[1]
  affect <- function(integ) integ$u[1] <- integ$u[1] + 1
  prob <- ODEProblem(function(u, p, t) 0 * u, 0, c(0, 10), p = 2)
  jp <- JumpProblem(prob, Direct(), VariableRateJump(rate, affect))
  sim <- solve(EnsembleProblem(jp), Tsit5(), trajectories = 60, saveat = 10, seed = 9)
  m <- mean(vapply(seq_len(60), function(i) tail(sim[[i]][1], 1), 0))
  expect_lt(abs(m - 20), 4 * sqrt(20 / 60))
})

# ---------------------------------------------------------------------------
# Aggregators, dependency graphs, SSAStepper with callbacks, save_positions

jagg_all <- function() list(Direct(), SortingDirect(), NRM(), DirectCR(), RSSA(), RSSACR(), Coevolve(),
                            FRM(), RDirect(), CCNRM(), DirectFW(), FRMFW())

test_that("dependency graphs are built as ModelingToolkit builds them", {
  # reference values printed by JumpProcesses 10 / ModelingToolkit 11 for
  # the same systems (mass action jumps are numbered first)
  t <- t_nounits
  S <- Variable("jdS", t); I <- Variable("jdI", t); R <- Variable("jdR", t)
  beta <- Parameter("jdbeta"); gamma <- Parameter("jdgamma")
  j1 <- ConstantRateJump(beta * S * I, list(S ~ Pre(S) - 1, I ~ Pre(I) + 1))
  maj <- SymbolicMassActionJump(gamma, list(I ~ 1), list(I ~ -1, R ~ 1))
  js <- complete(JumpSystem(list(j1, maj), t, list(S, I, R), list(beta, gamma), name = "jdsir"))
  op <- list(S ~ 990, I ~ 10, R ~ 0, beta ~ 1e-4, gamma ~ 0.01)
  for (agg in list(NRM(), CCNRM(), RDirect(), DirectCR(), SortingDirect(), Coevolve())) {
    jp <- JumpProblem(js, op, c(0, 250), aggregator = agg)
    expect_equal(jp$dependencies$dep_graph, list(1:2, 1:2))
    expect_null(jp$dependencies$vartojumps_map)
  }
  for (agg in list(RSSA(), RSSACR())) {
    jp <- JumpProblem(js, op, c(0, 250), aggregator = agg)
    expect_equal(jp$dependencies$vartojumps_map, list(2L, 1:2, integer()))
    expect_equal(jp$dependencies$jumptovars_map, list(2:3, 1:2))
  }
  jp <- JumpProblem(js, op, c(0, 250))
  expect_s3_class(jp$aggregator, "Direct")
  # mass action jumps first: the recovery is jump 1
  expect_equal(unname(jp$jumpset$rates(c(990, 10, 0), jp$prob$p, 0)), c(0.1, 1e-4 * 9900))
  X <- Variable("jdX", t); lam <- Parameter("jdlam"); mu <- Parameter("jdmu")
  bd <- complete(JumpSystem(list(SymbolicMassActionJump(lam, list(), list(X ~ 1)),
                                 SymbolicMassActionJump(mu, list(X ~ 1), list(X ~ -1))),
                            t, list(X), list(lam, mu), name = "jdbd"))
  bop <- list(X ~ 0, lam ~ 10, mu ~ 1)
  expect_equal(JumpProblem(bd, bop, c(0, 5), aggregator = NRM())$dependencies$dep_graph, list(1:2, 2L))
  d <- JumpProblem(bd, bop, c(0, 5), aggregator = RSSA())$dependencies
  expect_equal(d$vartojumps_map, list(2L))
  expect_equal(d$jumptovars_map, list(1L, 1L))
  # select_aggregator thresholds
  expect_s3_class(.select_aggregator(30, 0, 0, 0, list()), "RSSA")
  expect_s3_class(.select_aggregator(150, 0, 0, 0, list()), "RSSACR")
  expect_s3_class(.select_aggregator(10, 30, 0, 0, list(dep_graph = list())), "SortingDirect")
  expect_s3_class(.select_aggregator(10, 300, 0, 0, list(dep_graph = list())), "DirectCR")
  expect_s3_class(.select_aggregator(10, 30, 0, 0, list()), "Direct")
  expect_s3_class(.select_aggregator(0, 1, 0, 1, list()), "Coevolve")
  expect_s3_class(.select_aggregator(0, 1, 1, 0, list()), "Direct")
})

test_that("priority tables keep each rate in its power of two group", {
  set.seed(4)
  rates <- c(0, 1e-14, 3e-13, 0.7, 1, 2, 5.5, 1e3)
  pt <- .pt_new(-40L, rates)
  consistent <- function() {
    ok <- abs(pt$gsum - sum(rates)) <= 1e-9 * sum(rates)
    for (g in seq_along(pt$gmax)) {
      r <- rates[pt$gpids[[g]][seq_len(pt$gnum[g])]]
      if (g == 1L) ok <- ok && all(r <= .Machine$double.eps)
      if (g > 2L) ok <- ok && all(r < pt$gmax[g] & r >= pt$gmax[g] / 2)
      ok <- ok && abs(pt$gsums[g] - sum(r)) <= 1e-9 * max(1, sum(rates))
    }
    ok
  }
  ok <- consistent()
  for (k in 1:40) {
    i <- sample(length(rates), 1L)
    new <- if (k %% 5 == 0) 0 else stats::rexp(1L, 1 / 10)
    .pt_update(pt, i, rates[i], new)
    rates[i] <- new
    ok <- ok && consistent()
  }
  expect_true(ok)
  draws <- vapply(1:2000, function(i) .pt_sample(pt, rates), 0L)
  freq <- tabulate(draws, length(rates)) / 2000
  expect_lt(max(abs(freq - rates / sum(rates))), 0.04)
})

test_that("every aggregator gives the birth death law of Direct", {
  t <- t_nounits
  X <- Variable("jaX", t)
  lam <- Parameter("jalam", default = 10); mu <- Parameter("jamu", default = 1)
  bd <- complete(JumpSystem(list(SymbolicMassActionJump(lam, list(), list(X ~ 1)),
                                 SymbolicMassActionJump(mu, list(X ~ 1), list(X ~ -1))),
                            t, list(X), list(lam, mu), name = "jabd"))
  ts <- 1:3
  exact <- 10 * (1 - exp(-ts))
  n <- 150
  stats_of <- function(agg, seed) {
    jp <- JumpProblem(bd, list(X ~ 0), c(0, 3), aggregator = agg, save_positions = c(FALSE, FALSE))
    sim <- solve(EnsembleProblem(jp), SSAStepper(), trajectories = n, saveat = 1, seed = seed)
    vals <- t(vapply(seq_len(n), function(i) sim[[i]][X][2:4], numeric(3)))
    list(m = colMeans(vals), v = apply(vals, 2, stats::var))
  }
  ref <- stats_of(Direct(), 11)
  for (k in seq_along(jagg_all())) {
    agg <- jagg_all()[[k]]
    st <- stats_of(agg, 100 + k)
    # X(t) is Poisson with mean exact(t)
    expect_true(all(abs(st$m - exact) < 4 * sqrt(exact / n)), info = agg$name)
    expect_true(all(abs(st$m - ref$m) < 4 * sqrt(2 * exact / n)), info = agg$name)
    expect_true(all(abs(st$v / exact - 1) < 0.4), info = agg$name)
  }
})

test_that("every aggregator gives the SIR law of Direct", {
  t <- t_nounits
  S <- Variable("jsS", t); I <- Variable("jsI", t); R <- Variable("jsR", t)
  beta <- Parameter("jsbeta"); gamma <- Parameter("jsgamma")
  # infection as a constant rate jump, recovery as a mass action jump
  inf <- ConstantRateJump(beta * S * I, list(S ~ Pre(S) - 1, I ~ Pre(I) + 1))
  rec <- SymbolicMassActionJump(gamma, list(I ~ 1), list(I ~ -1, R ~ 1))
  sir <- complete(JumpSystem(list(inf, rec), t, list(S, I, R), list(beta, gamma), name = "jssir"))
  N <- 100
  op <- list(S ~ N - 10, I ~ 10, R ~ 0, beta ~ 0.3 / N, gamma ~ 0.1)
  n <- 60
  run <- function(agg, seed) {
    jp <- JumpProblem(sir, op, c(0, 20), aggregator = agg, save_positions = c(FALSE, FALSE))
    sim <- solve(EnsembleProblem(jp), SSAStepper(), trajectories = n, saveat = c(10, 20), seed = seed)
    t(vapply(seq_len(n), function(i) sim[[i]][R][2:3], numeric(2)))
  }
  ref <- run(Direct(), 21)
  for (k in seq_along(jagg_all())) {
    agg <- jagg_all()[[k]]
    r <- run(agg, 200 + k)
    se <- sqrt((apply(r, 2, stats::var) + apply(ref, 2, stats::var)) / n)
    expect_true(all(abs(colMeans(r) - colMeans(ref)) < 4 * se), info = agg$name)
    # conservation
    expect_true(all(r <= N))
  }
})

test_that("SSAStepper applies callbacks, events and tstops as JumpProcesses", {
  t <- t_nounits
  S <- Variable("jcS", t); I <- Variable("jcI", t); R <- Variable("jcR", t)
  beta <- Parameter("jcbeta"); gamma <- Parameter("jcgamma")
  j1 <- ConstantRateJump(beta * S * I, list(S ~ Pre(S) - 1, I ~ Pre(I) + 1))
  maj <- SymbolicMassActionJump(gamma, list(I ~ 1), list(I ~ -1, R ~ 1))
  js <- complete(JumpSystem(list(j1, maj), t, list(S, I, R), list(beta, gamma), name = "jcsir"))
  op <- list(S ~ 990, I ~ 10, R ~ 0, beta ~ 1e-4, gamma ~ 0.01)
  jp <- JumpProblem(js, op, c(0, 100), aggregator = Direct())
  cb <- DiscreteCallback(function(u, t, integ) t == 50,
                         function(integ) {
                           integ$u[1] <- integ$u[1] + 100
                           reset_aggregated_jumps(integ)
                         })
  sol <- solve(jp, callback = cb, tstops = 50, seed = 3)
  # as in Julia: the state before and after the affect are saved at t = 50
  k <- which(sol$t == 50)
  expect_length(k, 2L)
  expect_equal(unname(sol$u[[k[2]]][1] - sol$u[[k[1]]][1]), 100)
  # only the state after each jump is saved
  expect_equal(anyDuplicated(sol$t[-k]), 0L)
  # the same with an event of the system
  jsev <- complete(JumpSystem(list(j1, maj), t, list(S, I, R), list(beta, gamma), name = "jcsirev",
                              discrete_events = list(SymbolicDiscreteCallback(50, list(S ~ Pre(S) + 100)))))
  sev <- solve(JumpProblem(jsev, op, c(0, 100), aggregator = NRM()), seed = 3)
  k <- which(sev$t == 50)
  expect_length(k, 2L)
  expect_equal(unname(sev$u[[k[2]]][1] - sev$u[[k[1]]][1]), 100)
  # saveat with no jump saves gives the saveat times only
  s2 <- solve(JumpProblem(js, op, c(0, 100), save_positions = c(FALSE, FALSE)), saveat = 10, seed = 1)
  expect_equal(s2$t, seq(0, 100, by = 10))
  # termination and the error for continuous callbacks
  stop_cb <- DiscreteCallback(function(u, t, integ) u[3] >= 5, function(integ) terminate(integ))
  s3 <- solve(jp, callback = stop_cb, seed = 2)
  expect_equal(s3$retcode, "Terminated")
  expect_equal(unname(tail(s3[R], 1)), 5)
  expect_error(solve(jp, callback = ContinuousCallback(function(u, t, integ) t - 1, function(integ) NULL)),
               "continuous callbacks")
  # a pure birth process stopped by setting its rate to zero and resetting
  X <- Variable("jcX", t); lam <- Parameter("jclam", default = 5)
  pb <- complete(JumpSystem(list(SymbolicMassActionJump(lam, list(), list(X ~ 1))), t, list(X), list(lam),
                            name = "jcpb"))
  for (agg in jagg_all()) {
    off <- DiscreteCallback(function(u, t, integ) t == 2, function(integ) {
      integ$p[] <- 0
      reset_aggregated_jumps(integ)
    })
    sb <- solve(JumpProblem(pb, list(X ~ 0), c(0, 6), aggregator = agg), callback = off, tstops = 2, seed = 9)
    expect_equal(sb(2.5)[[1]], sb(6)[[1]], info = agg$name)
    expect_gt(sb(2)[[1]], 0)
  }
  # observed expressions use the parameters in force at each saved time
  rate <- sb[lam * X]
  expect_true(all(rate[sb$t > 2] == 0))
  expect_equal(rate[sb$t < 2], 5 * sb[X][sb$t < 2])
})

test_that("Coevolve thins bounded variable rate jumps in pure jump problems", {
  rate <- function(u, p, t) 1 + sin(t)
  affect <- function(integ) integ$u[1] <- integ$u[1] + 1
  vrj <- VariableRateJump(rate, affect, urate = function(u, p, t) 2, lrate = function(u, p, t) 0,
                          rateinterval = function(u, p, t) 1)
  dprob <- DiscreteProblem(function(u, p, t) u, 0, c(0, 6), numeric())
  jp <- JumpProblem(dprob, Coevolve(), vrj, dependencies = list(dep_graph = list(1L)))
  expect_s3_class(jp$aggregator, "Coevolve")
  sim <- solve(EnsembleProblem(jp), SSAStepper(), trajectories = 150, seed = 13)
  last <- vapply(seq_len(150), function(i) tail(sim[[i]]$u, 1)[[1]], 0)
  m <- 6 + 1 - cos(6)
  expect_lt(abs(mean(last) - m), 4 * sqrt(m / 150))
  # a general variable rate jump needs a continuous problem
  expect_error(JumpProblem(dprob, Direct(), VariableRateJump(rate, affect)), "continuous problem")
  expect_error(JumpProblem(dprob, NRM(), ConstantRateJump(function(u, p, t) 1, affect)), "dependency graph")
})

test_that("jumps coupled to ODEs use the aggregator and save as DiffEqBase", {
  t <- t_nounits; D <- D_nounits
  X <- Variable("joX", t); Y <- Variable("joY", t); k <- Parameter("jok", default = 2)
  crj <- ConstantRateJump(k, list(X ~ Pre(X) + 1))
  sys <- complete(System(list(D(Y) ~ -Y), t, list(X, Y), list(k), jumps = list(crj), name = "jo"))
  jp <- JumpProblem(sys, list(X ~ 0, Y ~ 1), c(0, 5))
  expect_s3_class(jp$prob, "ODEProblem")
  sol <- solve(jp, Tsit5(), saveat = 1, seed = 4)
  nj <- tail(sol[X], 1)
  # saveat points plus the states before and after each jump
  expect_equal(length(sol$t), 6 + 2 * nj)
  sol2 <- solve(JumpProblem(sys, list(X ~ 0, Y ~ 1), c(0, 5), save_positions = c(FALSE, FALSE)),
                Tsit5(), saveat = 1, seed = 4, reltol = 1e-8, abstol = 1e-10)
  expect_equal(sol2$t, 0:5)
  expect_equal(sol2[Y], exp(-(0:5)), tolerance = 1e-5)
  # E[X(t)] = k t
  sim <- solve(EnsembleProblem(JumpProblem(sys, list(X ~ 0, Y ~ 1), c(0, 5), aggregator = NRM(),
                                           save_positions = c(FALSE, FALSE))),
               Tsit5(), trajectories = 100, saveat = 5, seed = 6)
  last <- vapply(seq_len(100), function(i) tail(sim[[i]][X], 1), 0)
  expect_lt(abs(mean(last) - 10), 4 * sqrt(10 / 100))
  # variable rate jumps save after the jump by default (save_positions = (false, true))
  vrj <- VariableRateJump(k * (1 + Y), list(X ~ Pre(X) + 1))
  vsys <- complete(System(list(D(Y) ~ -Y), t, list(X, Y), list(k), jumps = list(vrj), name = "jov"))
  vs <- solve(JumpProblem(vsys, list(X ~ 0, Y ~ 1), c(0, 5)), Tsit5(), saveat = 1, seed = 8)
  expect_equal(length(vs$t), 6 + tail(vs[X], 1))
  # the dense output of the augmented state is projected on the unknowns
  v7 <- solve(JumpProblem(vsys, list(X ~ 0, Y ~ 1), c(0, 5)), Vern7(), seed = 8, reltol = 1e-8, abstol = 1e-10)
  expect_length(v7(2.5), 2L)
  expect_equal(v7(2.5, idxs = Y), exp(-2.5), tolerance = 1e-6)
})

test_that("tau leaping with regular and mass action jumps", {
  # simple_regular_solve.jl docstrings: decay with rate 0.1 u from 20
  rj <- RegularJump(function(u, p, t) 0.1 * u[1], function(u, p, t, counts, mark) -counts[1], 1)
  jp <- JumpProblem(DiscreteProblem(20, c(0, 2)), PureLeaping(), rj)
  s <- solve(jp, SimpleTauLeaping(), dt = 0.1, seed = 1)
  expect_equal(s$t, seq(0, 2, by = 0.1))
  v <- vapply(1:200, function(i) tail(solve(jp, SimpleTauLeaping(), dt = 0.01, seed = i)$u, 1)[[1]], 0)
  expect_lt(abs(mean(v) - 20 * exp(-0.2)), 4 * sqrt(20 * exp(-0.2) * (1 - exp(-0.2)) / 200))
  expect_error(solve(jp, SimpleTauLeaping()), "dt is required")
  expect_error(solve(jp, SSAStepper()), "tau leaping method")
  # isomerization A <-> B with index based mass action jumps (Julia's form)
  maj <- MassActionJump(c(1, 1), list(list(c(1, 1)), list(c(2, 1))),
                        list(list(c(1, -1), c(2, 1)), list(c(1, 1), c(2, -1))))
  jp2 <- JumpProblem(DiscreteProblem(c(150, 50), c(0, 10)), PureLeaping(), maj)
  for (alg in list(SimpleExplicitTauLeaping(), SimpleImplicitTauLeaping(), SimpleTrapezoidalLeaping(),
                   SimpleAdaptiveTauLeaping(), SimpleAdaptiveTauLeaping(eigenvalue_check = TRUE))) {
    r <- t(vapply(1:40, function(i) {
      s <- solve(jp2, alg, seed = i, saveat = c(1, 10))
      expect_equal(s$t, c(1, 10))
      s$u[[2L]]
    }, numeric(2)))
    expect_true(all(rowSums(r) == 200), label = alg$name)
    expect_lt(abs(mean(r[, 1]) - 100), 4 * sqrt(50 / 40), label = alg$name)
  }
  # the same mass action jumps through the SSA aggregators
  jp3 <- JumpProblem(DiscreteProblem(c(150, 50), c(0, 10)), NRM(), maj,
                     dependencies = list(dep_graph = list(1:2, 1:2)))
  r3 <- vapply(1:40, function(i) tail(solve(jp3, SSAStepper(), seed = i)$u, 1)[[1]][1], 0)
  expect_lt(abs(mean(r3) - 100), 4 * sqrt(50 / 40))
  # StochasticDiffEq's tau leaping methods on the SIR model as a RegularJump
  sir_rate <- function(u, p, t) c(0.1 / 1000 * u[1] * u[2], 0.01 * u[2])
  sir_c <- function(u, p, t, counts, mark) c(-counts[1], counts[1] - counts[2], counts[2])
  jsir <- JumpProblem(DiscreteProblem(c(999, 10, 0), c(0, 250)), PureLeaping(), RegularJump(sir_rate, sir_c, 2))
  ref <- mean(vapply(1:20, function(i) tail(solve(jsir, SimpleTauLeaping(), dt = 0.1, seed = i)$u, 1)[[1]][3], 0))
  for (alg in list(TauLeaping(), CaoTauLeaping(), ImplicitTauLeaping(), ThetaTrapezoidalTauLeaping())) {
    fin <- t(vapply(1:20, function(i) tail(solve(jsir, alg, dt = 1, seed = i)$u, 1)[[1]], numeric(3)))
    # the implicit methods move the state off the integers; the total holds
    expect_lt(max(abs(rowSums(fin) - 1009)), 1e-6, label = alg$name)
    expect_lt(abs(mean(fin[, 3]) - ref), 60, label = alg$name)
  }
  expect_error(solve(jsir, ImplicitTauLeaping()), "require a choice of dt")
  expect_error(MassActionJump(c(1, 1), list(list(c(1, 1))), list()), "one reactant")
  expect_error(.leap_hor(list(maj_rs = list(cbind(1L, 4)), n = 1L)), "order 4")
})

test_that("spatial jumps on grids and graphs", {
  # topology.jl docstrings
  grid <- CartesianGrid(c(3, 4))
  expect_equal(num_sites(grid), 12L)
  expect_equal(neighbors(grid, 1), c(2L, 4L))
  expect_equal(outdegree(CartesianGrid(c(2, 2)), 1), 2L)
  expect_equal(neighbors(CartesianGridRej(c(2, 2)), 1), c(2L, 3L))
  expect_equal(num_sites(CartesianGridRej(3, 2)), 8L)
  expect_equal(neighbors(matrix(c(0, 1, 1, 0), 2), 1), 2L)
  # diffusion on a line: the long time distribution is uniform
  u0 <- matrix(c(50, 0, 0, 0, 0), 1)
  for (agg in list(NSM(), DirectCRDirect())) {
    jp <- JumpProblem(DiscreteProblem(u0, c(0, 20)), agg, list(), hopping_constants = 1,
                      spatial_system = CartesianGrid(5))
    fin <- t(vapply(1:20, function(i) tail(solve(jp, SSAStepper(), seed = i)$u, 1)[[1]], numeric(5)))
    expect_true(all(rowSums(fin) == 50), label = agg$name)
    expect_lt(max(abs(colMeans(fin) - 10)), 4 * sqrt(8 / 20), label = agg$name)
  }
  # A + B -> C with hops on a graph against the flattened problem
  adj <- list(c(2, 3), c(1, 4), c(1, 4), c(2, 3))
  smaj <- SpatialMassActionJump(0.01, NULL, list(list(c(1, 1), c(2, 1))), list(list(c(1, -1), c(2, -1), c(3, 1))))
  u0 <- matrix(c(30, 30, 0), 3, 4)
  dp <- DiscreteProblem(u0, c(0, 2))
  jp <- JumpProblem(dp, NSM(), smaj, hopping_constants = c(1, 0.5, 0.1), spatial_system = adj)
  cs <- vapply(1:30, function(i) sum(matrix(tail(solve(jp, SSAStepper(), seed = i)$u, 1)[[1]], 3)[3, ]), 0)
  fl <- flatten_spatial(smaj, dp, adj, c(1, 0.5, 0.1))
  expect_length(fl$jumps, 3 * 8 + 4)
  jf <- JumpProblem(fl$prob, Direct(), fl$jumps)
  cf <- vapply(1:30, function(i) sum(matrix(tail(solve(jf, SSAStepper(), seed = i)$u, 1)[[1]], 3)[3, ]), 0)
  expect_lt(abs(mean(cs) - mean(cf)), 4 * sqrt((stats::var(cs) + stats::var(cf)) / 30))
  # the other forms of hopping constants
  hcl <- matrix(lapply(1:12, function(k) c(1, 2)), 3, 4)
  expect_equal(.sp_hop_weights(hcl, adj, 3)[[2, 3]], c(1, 2))
  w <- .sp_hop_weights(list(c(1, 2, 3), lapply(1:4, function(i) c(0.5, 2))), adj, 3)
  expect_equal(w[[3, 2]], c(1.5, 6))
  expect_equal(.sp_hop_weights(matrix(1:12, 3), adj, 3)[[2, 4]], c(11, 11))
  sm <- SpatialMassActionJump(NULL, matrix(1:4, 1), list(list(c(1, 2))), list(list(c(1, -2))))
  expect_equal(sm$spatial_rates, matrix(1:4 / 2, 1))
  expect_error(JumpProblem(dp, NSM(), smaj, spatial_system = adj), "hopping_constants")
  expect_error(JumpProblem(dp, Direct(), smaj, hopping_constants = 1, spatial_system = adj), "spatial aggregator")
})
