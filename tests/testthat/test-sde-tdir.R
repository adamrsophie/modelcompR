# SDE problems backward in time (StochasticDiffEq test/tdir_tests.jl) and
# zero length spans; reference values from tools/review-tdir-probe.jl.

.tdir_z <- function(k) (37 * k) %% 101 / 50 - 1
.tdir_fd <- function(u, p, t) -0.5 * u + 0.1 * t
.tdir_gd <- function(u, p, t) 0.3 * u
.tdir_prob <- function() {
  hb <- 1 / 64
  tb <- 1 - (0:64) * hb
  Wb <- list(c(0, 0)); Zb <- list(c(0, 0))
  for (k in 1:64) {
    Wb[[k + 1]] <- Wb[[k]] + sqrt(hb) * c(.tdir_z(2 * k - 1), .tdir_z(2 * k))
    Zb[[k + 1]] <- Zb[[k]] + sqrt(hb) * c(.tdir_z(7919 + 2 * k - 1), .tdir_z(7919 + 2 * k))
  }
  SDEProblem(.tdir_fd, .tdir_gd, c(1, 0.5), c(1, 0), noise = NoiseGrid(tb, Wb, Zb))
}

test_that("SDE problems integrate backward in time (tdir_tests.jl)", {
  for (prob in list(SDEProblem(function(u, p, t) 1.01 * u, function(u, p, t) 0.87 * u, 0.5, c(1, 0)),
                    SDEProblem(.tdir_fd, .tdir_gd, c(1, 0.5), c(1, 0)))) {
    solEM <- solve(prob, EM(), dt = -1 / 16, tstops = 0.33, seed = 100)
    expect_gt(solEM$t[15], solEM$t[16])
    expect_true(any(abs(solEM$t - 0.33) < 1e-12))
    expect_equal(tail(solEM$t, 1), 0)
    solSRIW1 <- solve(prob, SRIW1(), tstops = 0.33, seed = 100)
    expect_identical(solSRIW1$retcode, "Success")
    expect_gt(solSRIW1$t[15], solSRIW1$t[16])
  }
  # a positive dt on a reversed span steps backward
  s <- solve(.tdir_prob(), EM(), dt = 1 / 64)
  expect_equal(s$t, 1 - (0:64) / 64)
})

test_that("a reversed NoiseGrid retraces the forward path (tdir_tests.jl)", {
  tf <- (0:100) * 0.01
  Wf <- cumsum(c(0, 0.1 * .tdir_z(1:100 + 11)))
  f <- function(u, p, t) u
  sol <- solve(SDEProblem(f, f, 0.5, c(0, 1), noise = NoiseGrid(tf, Wf)), EulerHeun(), dt = 0.01)
  W3 <- NoiseGrid(rev(sol$t), rev(Wf))
  sol2 <- solve(SDEProblem(f, f, tail(sol$u, 1)[[1]], c(1, 0), noise = W3), EulerHeun(), dt = 0.01)
  expect_equal(unlist(sol$u), rev(unlist(sol2$u)), tolerance = 0.1)
  expect_equal(tail(sol2$u, 1)[[1]], 0.500282401688184, tolerance = 1e-13)
})

test_that("zero length spans return the trivial solution", {
  for (alg in list(EM(), SRIW1(), SOSRI(), LambaEM())) {
    s <- solve(SDEProblem(.tdir_fd, .tdir_gd, c(1, 0.5), c(0, 0)), alg, dt = 0.1)
    expect_identical(s$retcode, "Success")
    expect_equal(s$t, 0)
    expect_equal(s$u[[1]], c(1, 0.5))
  }
  expect_identical(solve(SDEProblem(.tdir_fd, .tdir_gd, c(1, 0.5), c(0, 0)))$retcode, "Success")
  ig <- init(SDEProblem(.tdir_fd, .tdir_gd, c(1, 0.5), c(0, 0)), EM(), dt = 0.1)
  expect_equal(solve_integrator(ig)$t, 0)
  d <- solve(DDEProblem(function(u, h, p, t) -h(p, t - 1), 1, function(p, t) 1, c(0, 0), constant_lags = 1),
             MethodOfSteps(Tsit5()))
  expect_identical(d$retcode, "Success")
  expect_equal(d$t, 0)
  expect_equal(d$u[[1]], 1)
})

test_that("backward solves on a reversed NoiseGrid match StochasticDiffEq", {
  # tools/review-tdir-probe.jl, section D; dt = -1/64 on the grid points
  ref <- list(
    EM = c(1.4170096758814468, 0.8184485009933093), EulerHeun = c(1.4423707962196652, 0.8337059595307739),
    LambaEulerHeun = c(1.4423707962196652, 0.8337059595307739), RKMil = c(1.3840249177244819, 0.7977729500732498),
    SimplifiedEM = c(1.753896401225358, 0.6643439388985722), SRIW1 = c(1.3784598973659747, 0.7963010023732773),
    SOSRI = c(1.378464078471383, 0.7962679877436438), SRA1 = c(1.4209767036801484, 0.8207664736494807),
    SOSRA = c(2.3953065289343876, 1.3663431646227349), ImplicitEM = c(1.4239314362762432, 0.823388952287411),
    ImplicitRKMil = c(1.3893548748843436, 0.8027671556149836), SKenCarp = c(1.4209582698168928, 0.8207820279583288),
    SROCK2 = c(1.3779426688409555, 0.7961334067591466), SKSROCK = c(1.4188004108234147, 0.8190137633045653),
    WangLi3SMil_A = c(1.4604351148252095, 0.843869252920433), WangLi3SMil_B = c(1.3747117535631013, 0.7938881321464843),
    COM = c(1.4415306020979513, 0.8327657831501437), SMEB = c(1.4042420927927584, 0.8104323627043933))
  for (nm in names(ref)) {
    s <- solve(.tdir_prob(), get(nm)(), dt = -1 / 64, adaptive = FALSE)
    expect_equal(tail(s$u, 1)[[1]], ref[[nm]], tolerance = 1e-11, info = nm)
  }
  s <- solve(.tdir_prob(), EM(), dt = -1 / 16, tstops = 0.33)
  expect_equal(tail(s$u, 1)[[1]], c(1.4239993087906928, 0.8213308474872024), tolerance = 1e-13)
  # adaptive: steps and final values
  ad <- list(SRIW1 = list(c(1.3791780803018252, 0.7966908967930357), c(33, 0)),
             SOSRI = list(c(1.3782973391321032, 0.7962449085740613), c(145, 1)),
             LambaEM = list(c(1.4216988712792182, 0.8163051987444814), c(18, 0)),
             RKMil = list(c(1.38994550098725, 0.7950800365735086), c(15, 0)))
  for (nm in names(ad)) {
    s <- solve(.tdir_prob(), get(nm)(), dt = -1 / 64)
    expect_equal(tail(s$u, 1)[[1]], ad[[nm]][[1]], tolerance = 1e-12, info = nm)
    expect_equal(c(s$stats$naccept, s$stats$nreject), ad[[nm]][[2]], info = nm)
  }
  # methods whose upstream step takes sqrt(dt) stop (a DomainError upstream)
  for (alg in list(RI1(), DRI1(), PL1WM(), NON(), SROCKC2(), W2Ito1()))
    expect_error(solve(.tdir_prob(), alg, dt = -1 / 64, adaptive = FALSE), "cannot integrate backward")
})

test_that("adaptive RKMil no longer stalls (newest StochasticDiffEqLowOrder estimate)", {
  # tools/review-tdir-probe.jl, section C: 21 accepted steps on the grid
  z <- .tdir_z
  hg <- 1 / 4096
  tg <- (0:4096) * hg
  Wg <- cumsum(c(0, sqrt(hg) * z(1:4096)))
  Zg <- cumsum(c(0, sqrt(hg) * z(1:4096 + 7919)))
  prob <- SDEProblem(function(u, p, t) -u, function(u, p, t) 0.1 * u, 1, c(0, 1), noise = NoiseGrid(tg, Wg, Zg))
  s <- solve(prob, RKMil(), maxiters = 2000)
  expect_identical(s$retcode, "Success")
  expect_equal(c(s$stats$naccept, s$stats$nreject), c(21, 0))
  expect_equal(tail(s$u, 1)[[1]], 0.35825173754918405, tolerance = 1e-3)
  s <- solve(SDEProblem(function(u, p, t) -u, function(u, p, t) 0.1 * u, 1, c(0, 1)), RKMil(), seed = 7, maxiters = 2000)
  expect_identical(s$retcode, "Success")
})
