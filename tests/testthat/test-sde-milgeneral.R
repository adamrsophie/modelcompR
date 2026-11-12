# RKMilGeneral on fixed noise paths. Reference values from
# tools/stochasticdiffeq-probe.jl (StochasticDiffEq). There the second
# process carries the raw Fourier coefficients z(k) per grid interval; here
# it carries sqrt(h) z(k), the increments of a Wiener process, which the
# step divides by sqrt(h) again, so the coefficients agree.

rkmg_f <- function(u, p, t) -u + 0.2 * u^2
rkmg_gdiag <- function(u, p, t) 0.3 * u + 0.1 * u^2
rkmg_gnd <- function(u, p, t) matrix(c(0.2 * u[1], 0.05 * u[1] * u[2], 0.1 * u[2], 0.3 * u[2]), 2, 2)

test_that("RKMilGeneral matches StochasticDiffEq for diagonal noise", {
  gr <- sdeprobe_grid(2, 1 / 64, 80, off = 3)
  mk <- function() SDEProblem(rkmg_f, rkmg_gdiag, c(0.8, 0.4), c(0, 1), noise = NoiseGrid(gr$t, gr$W))
  ito <- solve(mk(), RKMilGeneral(), dt = 1 / 16, adaptive = FALSE)
  expect_equal(unlist(ito$u[c(2, 9, 17)]),
               c(0.77746495488, 0.367671410745, 0.5263152169170552, 0.24834812714898344,
                 0.34652910167782514, 0.14248032419472237), tolerance = 1e-13)
  st <- solve(mk(), RKMilGeneral(interpretation = "Stratonovich"), dt = 1 / 16, adaptive = FALSE)
  expect_equal(unlist(st$u[c(2, 9, 17)]),
               c(0.77962826848, 0.36827015082, 0.5371071619098989, 0.2512746944724939,
                 0.359085943615432, 0.14568809539692457), tolerance = 1e-13)
})

test_that("RKMilGeneral matches StochasticDiffEq for nondiagonal noise", {
  gr <- sdeprobe_grid(2, 1 / 64, 80, off = 13, mz = 15)
  mk <- function() SDEProblem(rkmg_f, rkmg_gnd, c(0.8, 0.4), c(0, 1), noise = NoiseGrid(gr$t, gr$W, gr$Z),
                              noise_rate_prototype = matrix(0, 2, 2))
  # one grid interval per step: the Levy area from the step's coefficients
  s <- solve(mk(), RKMilGeneral(p = 3), dt = 1 / 64, adaptive = FALSE)
  expect_equal(unlist(s$u[c(2, 17, 33, 65)]),
               c(0.7745847563441802, 0.39260255366404, 0.6437721297785541, 0.31227338475944444,
                 0.5269449508110048, 0.2359490817870849, 0.330975447744263, 0.1507155760680965),
               tolerance = 1e-13)
  s <- solve(mk(), RKMilGeneral(p = 3, interpretation = "Stratonovich"), dt = 1 / 64, adaptive = FALSE)
  expect_equal(unlist(s$u[c(2, 17, 33, 65)]),
               c(0.7747613188441802, 0.392778379835915, 0.6460768347671508, 0.31450752163385215,
                 0.530674901915536, 0.23933527055156129, 0.3356318020811785, 0.1550320092725619),
               tolerance = 1e-13)
  # four grid intervals per step: Chen's relation over the grid
  s <- solve(mk(), RKMilGeneral(p = 3), dt = 1 / 16, adaptive = FALSE)
  expect_equal(unlist(s$u[c(2, 5, 9, 17)]),
               c(0.7594362032991216, 0.38834659822900347, 0.644020199931606, 0.31234577962506954,
                 0.5277335533701459, 0.23538905887415304, 0.33115756287320675, 0.15086596250202888),
               tolerance = 1e-13)
})

test_that("adaptive RKMilGeneral reproduces StochasticDiffEq with its defects restored", {
  # the four deliberate differences (see R/stochastic-milgeneral.R) undone:
  # the diagonal error estimate with (du2 - L) / sqrt(dt), raw coefficients,
  # no partial grid intervals and the reset norm of the out of place step
  local_mocked_bindings(.rkmg_diag_en = function(dW, ggprime, du2, L, sqdt) dW^3 * ((du2 - L) / sqdt)^2 / 6,
                        .rkmg_scale = function(dZ, h) dZ,
                        .rkmg_grid_ends = function(first, inner, last) inner,
                        .rkmg_norm_acc = function(acc, x) x)
  gr <- sdeprobe_grid(2, 1 / 1024, 1100, off = 9)
  s <- solve(SDEProblem(rkmg_f, rkmg_gdiag, c(0.8, 0.4), c(0, 1), noise = NoiseGrid(gr$t, gr$W)),
             RKMilGeneral(), dt = 1 / 16, abstol = 1e-3, reltol = 1e-3)
  expect_length(s$t, 847L)
  expect_equal(s$t[1:6], c(0, 0.0625, 0.075, 0.08906249999999999, 0.10488281249999998, 0.11058682714525642),
               tolerance = 1e-14)
  expect_equal(c(s$u[[5]], s$u[[847]]),
               c(0.7249216732520074, 0.36219297567578396, 0.30815890357003783, 0.1476562369361649), tolerance = 1e-13)
  gr <- sdeprobe_grid(2, 1 / 1024, 1100, off = 17, mz = 15, zscale = 1)
  s <- solve(SDEProblem(rkmg_f, rkmg_gnd, c(0.8, 0.4), c(0, 1), noise = NoiseGrid(gr$t, gr$W, gr$Z),
                        noise_rate_prototype = matrix(0, 2, 2)),
             RKMilGeneral(p = 3), dt = 1 / 16, abstol = 1e-3, reltol = 1e-3)
  expect_length(s$t, 23L)
  expect_equal(s$t[1:6], c(0, 0.0625, 0.075, 0.08906249999999999, 0.10488281249999998, 0.12268066406249997),
               tolerance = 1e-14)
  expect_equal(c(s$u[[5]], s$u[[23]]),
               c(0.7303323316226291, 0.366617232730747, 0.3208838384443848, 0.1495330262161068), tolerance = 1e-13)
})

test_that("the corrected diagonal error estimate does not stall", {
  gr <- sdeprobe_grid(2, 1 / 1024, 1100, off = 9)
  s <- solve(SDEProblem(rkmg_f, rkmg_gdiag, c(0.8, 0.4), c(0, 1), noise = NoiseGrid(gr$t, gr$W)),
             RKMilGeneral(), dt = 1 / 16, abstol = 1e-3, reltol = 1e-3)
  expect_lt(length(s$t), 100L)
  expect_equal(utils::tail(s$t, 1), 1)
})

test_that("RKMilGeneral sizes the second process and runs without a grid", {
  expect_equal(RKMilGeneral(p = 3)$z_length(2, 0.1), 15L)
  expect_equal(RKMilGeneral()$z_length(2, NULL), coefficient_length(2, 10, MronRoe()))
  prob <- SDEProblem(rkmg_f, rkmg_gnd, c(0.8, 0.4), c(0, 1), noise_rate_prototype = matrix(0, 2, 2))
  s <- solve(prob, RKMilGeneral(p = 5), dt = 0.1, adaptive = FALSE, seed = 1)
  expect_length(s$t, 11L)
  expect_true(all(is.finite(unlist(s$u))))
  # after the rejected first step the iterated integrals come from the
  # stored noise pieces of the step (Chen's relation)
  s <- solve(prob, RKMilGeneral(), dt = 0.1, seed = 2, abstol = 1e-3, reltol = 1e-3)
  expect_gt(s$stats$nreject, 0)
  expect_equal(utils::tail(s$t, 1), 1)
  expect_true(all(is.finite(unlist(s$u))))
})
