test_that("Wiener processes simulate, bridge and replay", {
  sol <- solve(NoiseProblem(WienerProcess(0, 0), c(0, 1)), dt = 0.1, seed = 1)
  expect_s3_class(sol, "NoiseProcess")
  expect_equal(sol$t, seq(0, 1, by = 0.1))
  n <- length(sol$t)
  v <- sol(0.55)
  expect_length(sol$t, n + 1L)
  expect_identical(sol(0.55), v)
  set.seed(2)
  w1 <- vapply(1:300, function(i) {
    s <- solve(NoiseProblem(WienerProcess(0, 0), c(0, 1)), dt = 0.25)
    s$W[[length(s$W)]]
  }, 0)
  expect_equal(var(w1), 1, tolerance = 0.25)
  wr <- NoiseWrapper(sol)
  expect_equal(wr(0.5), sol(0.5))
  expect_s3_class(pCN(sol, 0.9), "NoiseWrapper")
  sw <- solve(NoiseProblem(SimpleWienerProcess(0, c(0, 0)), c(0, 1)), dt = 0.5, seed = 3)
  expect_length(sw$W[[3]], 2L)
  expect_error(reject_step(SimpleWienerProcess(0, 0), 0.1), "adaptivity rejections")
})

test_that("special processes have their exact distributions", {
  set.seed(4)
  ou <- vapply(1:300, function(i) {
    s <- solve(NoiseProblem(OrnsteinUhlenbeckProcess(2, 1, 0.5, 0, 3), c(0, 3)), dt = 1)
    s$W[[length(s$W)]]
  }, 0)
  # stationary mean mu, variance sigma^2 / (2 Theta)
  expect_equal(mean(ou), 1, tolerance = 0.1)
  expect_equal(var(ou), 0.5^2 / 4, tolerance = 0.3)
  bb <- solve(NoiseProblem(BrownianBridge(0, 1, 0, 2), c(0, 1)), dt = 0.1, seed = 5)
  expect_equal(bb$W[[length(bb$W)]], 2)
  gbm <- solve(NoiseProblem(GeometricBrownianMotionProcess(0.1, 0.2, 0, 1), c(0, 1)), dt = 0.1, seed = 6)
  expect_true(all(unlist(gbm$W) > 0))
  cp <- solve(NoiseProblem(CompoundPoissonProcess(5, 0, 0), c(0, 2)), dt = 0.1, seed = 7)
  w <- unlist(cp$W)
  expect_true(all(w == round(w)) && all(diff(w) >= 0))
  cw <- solve(NoiseProblem(CorrelatedWienerProcess(matrix(c(1, 0.9, 0.9, 1), 2), 0, c(0, 0)), c(0, 1)),
              dt = 0.01, seed = 8)
  M <- do.call(rbind, cw$W)
  expect_equal(cor(diff(M))[1, 2], 0.9, tolerance = 0.1)
  set.seed(9)
  A <- vapply(1:200, function(i) .bwt_sample(1, 1e-6)$A, 0)
  # the area of BoxWedgeTail is twice the Levy area; its variance is 1 per unit time
  expect_equal(var(A), 1, tolerance = 0.3)
  expect_error(BoxWedgeTail(0, 0), "2d Brownian")
})

test_that("path noises", {
  g <- NoiseGrid(seq(0, 1, 0.5), c(0, 1, 3))
  expect_equal(g(0.75), 2)
  expect_error(g(2), "extrapolate")
  nf <- NoiseFunction(0, function(u, p, t) sin(t))
  expect_equal(nf(1), sin(1))
  vbt <- VirtualBrownianTree(0, 0, tend = 1, seed = 10)
  expect_identical(vbt(0.3), vbt(0.3))
  nt <- NoiseTransport(0, function(u, p, t, rv) rv * t, function() 2)
  expect_equal(nt(3), 6)
})

test_that("SDE and RODE solves take noise processes", {
  f <- function(u, p, t) -u
  g <- function(u, p, t) 0.1
  ou <- OrnsteinUhlenbeckProcess(1, 0, 0.3, 0, 1)
  s <- solve(SDEProblem(f, g, 1, c(0, 1), noise = ou), EM(), dt = 0.01, seed = 2)
  expect_equal(s$retcode, "Success")
  # a grid noise gives a deterministic increment path
  tg <- seq(0, 1, by = 0.01)
  Wg <- cumsum(c(0, rep(0.01, 100)))
  sg <- solve(SDEProblem(function(u, p, t) 0, function(u, p, t) 1, 0, c(0, 1), noise = NoiseGrid(tg, Wg)),
              EM(), dt = 0.1)
  expect_equal(tail(sg$u, 1)[[1]], Wg[101], tolerance = 1e-10)
  # adaptive rejections keep the Brownian path consistent with RSwM
  sa <- solve(SDEProblem(function(u, p, t) 1.01 * u, function(u, p, t) 0.87 * u, 0.5, c(0, 1),
                         noise = WienerProcess(0, 0)), SRIW1(), seed = 3)
  expect_equal(sa$retcode, "Success")
  fr <- function(u, p, t, W) 1.01 * u + 0.87 * u * W
  r <- solve(RODEProblem(fr, 0.5, c(0, 1)), RandomEM(), dt = 0.01, seed = 1)
  expect_s3_class(r, "RODESolution")
  # the RODE with W fixed at 0 is the ODE u' = 1.01 u
  r0 <- solve(RODEProblem(fr, 0.5, c(0, 1), noise = NoiseFunction(0, function(u, p, t) 0)), RandomHeun(), dt = 0.001)
  expect_equal(tail(r0$u, 1)[[1]], 0.5 * exp(1.01), tolerance = 1e-5)
  rt <- solve(RODEProblem(fr, 0.5, c(0, 1), noise = NoiseFunction(0, function(u, p, t) 0)), RandomTamedEM(), dt = 0.001)
  expect_equal(tail(rt$u, 1)[[1]], 0.5 * exp(1.01), tolerance = 1e-2)
  tg <- seq(0, 1, length.out = 1001)
  set.seed(11)
  Wf <- c(0, cumsum(stats::rnorm(1000, 0, sqrt(1 / 1000))))
  r15 <- solve(RODEProblem(function(u, p, t, W) -u + sin(W), 1, c(0, 1), noise = NoiseGrid(tg, Wf)),
               RandomTaylor15(), dt = 0.05)
  ref <- solve(RODEProblem(function(u, p, t, W) -u + sin(W), 1, c(0, 1), noise = NoiseGrid(tg, Wf)),
               RandomHeun(), dt = 0.001)
  expect_equal(tail(r15$u, 1)[[1]], tail(ref$u, 1)[[1]], tolerance = 1e-2)
  expect_error(solve(RODEProblem(fr, 0.5, c(0, 1)), RandomTaylor15(), dt = 0.1), "NoiseGrid")
})

test_that("weak random variables", {
  expect_equal(calc_twopoint_random(0.2, c(-1, 0.5, 0)), c(-0.2, 0.2, -0.2))
  set.seed(12)
  z <- stats::rnorm(6000)
  y <- calc_threepoint_random(sqrt(0.12), stats::qnorm(1 / 6), z)
  expect_setequal(unique(y), c(-sqrt(0.12), 0, sqrt(0.12)))
  expect_equal(mean(y == 0), 2 / 3, tolerance = 0.05)
  expect_equal(mean(y^2), 0.12 / 3 * 3 / 3 * 1, tolerance = 0.1)
})
