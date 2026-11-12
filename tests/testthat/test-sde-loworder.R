# Reference values: one fixed noise path (NoiseGrid) through StochasticDiffEq 6
# (StochasticDiffEqLowOrder and StochasticDiffEqMilstein), two steps of 0.1.

test_that("low order SDE methods match StochasticDiffEq on a fixed noise path", {
  f <- function(u, p, t) c(-0.5 * u[1] + 0.2 * u[2] * u[1], 0.3 * u[1] - u[2])
  gD <- function(u, p, t) c(0.3 * u[1] + 0.1, 0.2 * u[2]^2)
  gG <- function(u, p, t) matrix(c(0.3 * u[1], 0.05 * u[1] * u[2], 0.1 * u[2], 0.2 * u[2]), 2, 2)
  fs <- function(u, p, t) -0.5 * u + 0.2 * u^2
  gs <- function(u, p, t) 0.3 * u + 0.1 * u^2
  ggprime <- list(diag = function(u, p, t) c(0.3 * (0.3 * u[1] + 0.1), 0.4 * u[2] * 0.2 * u[2]^2),
                  scal = function(u, p, t) (0.3 + 0.2 * u) * gs(u, p, t),
                  gen = function(u, p, t) c(0.01 * u[1], 0.02 * u[2]))
  ts <- c(0, 0.1, 0.2)
  W <- list(c(0, 0), c(0.4, -0.1), c(0.38, 0.35))
  ref <- list(
    SimplifiedEM_diag = c(0.90755691383799364, 0.46399207056780106),
    SimplifiedEM_scal = c(0.73427661445940995),
    SimplifiedEM_gen = c(0.91244150924033396, 0.46318795820064224),
    PCEuler_diag = c(1.05922354983750000, 0.48284629151762459),
    PCEuler_scal = c(0.85149067204505502),
    PCEuler_gen = c(1.05012275594209248, 0.50791515870661819),
    WangLi3SMil_A_diag = c(1.05797464345673031, 0.48053760514443228),
    WangLi3SMil_A_scal = c(0.84879499869771347),
    WangLi3SMil_B_diag = c(1.05693456995189439, 0.48442188327564223),
    WangLi3SMil_B_scal = c(0.85153090804807385),
    WangLi3SMil_C_diag = c(1.05797473763065475, 0.48052271112290934),
    WangLi3SMil_C_scal = c(0.84835052712838388),
    WangLi3SMil_D_diag = c(1.05807475851876420, 0.48015270117314535),
    WangLi3SMil_D_scal = c(0.84819804071569638),
    WangLi3SMil_E_diag = c(1.05702504317356927, 0.48405376216250851),
    WangLi3SMil_E_scal = c(0.85136657003419713),
    WangLi3SMil_F_diag = c(1.05702476589528782, 0.48403790744329428),
    WangLi3SMil_F_scal = c(0.85088802258663176),
    SplitEM_diag = c(1.09025800000000017, 0.47606874999999993)
  )
  for (nm in names(ref)) {
    nt <- sub(".*_", "", nm)
    alg_name <- sub("_(diag|scal|gen)$", "", nm)
    alg <- if (alg_name == "PCEuler") PCEuler(ggprime[[nt]]) else get(alg_name)()
    prob <- switch(nt,
      diag = SDEProblem(f, gD, c(1, 0.5), c(0, 0.2), noise = NoiseGrid(ts, W)),
      scal = SDEProblem(fs, gs, 0.8, c(0, 0.2), noise = NoiseGrid(ts, c(0, 0.4, 0.38))),
      gen = SDEProblem(f, gG, c(1, 0.5), c(0, 0.2), noise = NoiseGrid(ts, W), noise_rate_prototype = matrix(0, 2, 2)))
    if (alg_name == "SplitEM")
      prob <- SplitSDEProblem(function(u, p, t) -0.5 * u, function(u, p, t) 0.2 * u^2, gD, c(1, 0.5), c(0, 0.2),
                              noise = NoiseGrid(ts, W))
    sol <- solve(prob, alg, dt = 0.1, adaptive = FALSE)
    expect_equal(as.numeric(sol$u[[3]]), ref[[nm]], tolerance = 1e-13, info = nm)
  }
})

test_that("Wang and Li's methods converge strongly at order 1 and refuse general noise", {
  a <- 1; b <- 0.5; N <- 200
  prob <- SDEProblem(function(u, p, t) a * u, function(u, p, t) b * u, rep(1, N), c(0, 1))
  dts <- 2^-(3:7)
  for (alg in list(WangLi3SMil_A(), WangLi3SMil_E())) {
    errs <- vapply(seq_along(dts), function(i) {
      sol <- solve(prob, alg, dt = dts[i], adaptive = FALSE, seed = 100 + i)
      W <- sol$W$W
      mean(abs(sol$u[[length(sol$t)]] - exp((a - b^2 / 2) + b * W[nrow(W), ])))
    }, 0)
    expect_gt(unname(coef(lm(log(errs) ~ log(dts)))[2]), 0.8)
  }
  probG <- SDEProblem(function(u, p, t) -u, function(u, p, t) diag(0.1, 2), c(1, 1), c(0, 1),
                      noise_rate_prototype = matrix(0, 2, 2))
  expect_error(solve(probG, WangLi3SMil_A(), dt = 0.1, seed = 1), "diagonal or scalar")
  expect_error(PCEuler(), "ggprime")
})

test_that("SOSRI and SOSRI2 match StochasticDiffEq on a fixed noise path", {
  f <- function(u, p, t) c(-0.5 * u[1] + 0.2 * u[2] * u[1], 0.3 * u[1] - u[2])
  gD <- function(u, p, t) c(0.3 * u[1] + 0.1, 0.2 * u[2]^2)
  ts <- c(0, 0.1, 0.2)
  W <- list(c(0, 0), c(0.4, -0.1), c(0.38, 0.35))
  Z <- list(c(0, 0), c(0.2, -0.3), c(-0.1, -0.5))
  ref <- list(SOSRI = c(1.0589027800400508, 0.4824300514578158), SOSRI2 = c(1.0591196988090228, 0.4825637184059167))
  for (nm in names(ref)) {
    sol <- solve(SDEProblem(f, gD, c(1, 0.5), c(0, 0.2), noise = NoiseGrid(ts, W, Z)), get(nm)(), dt = 0.1,
                 adaptive = FALSE)
    expect_equal(sol$u[[3]], ref[[nm]], tolerance = 1e-13, info = nm)
  }
})
