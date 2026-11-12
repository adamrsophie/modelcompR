# Composite SDE algorithms with stiffness switching on a fixed noise path,
# against tools/stochasticdiffeq-probe.jl (StochasticDiffEq). The problem
# turns stiff in the middle of [0, 1].
#
# Two upstream defects are not reproduced:
# * after a switch to the stiff method with dtfac != 1 the step grows, but
#   upstream only shrinks the noise to the step size, so that step uses the
#   increment of the old size; the dtfac = 2 and stiffalgfirst references
#   come from the probe run with the noise following the step both ways;
# * the members' steps read their traits from the composite, so LambaEM()
#   inside a composite runs unsplit upstream; the reference for
#   LambaEM(split = FALSE) here is upstream's AutoSOSRI2(LambaEM()).

sdecomp_run <- function(alg, g = function(u, p, t) 0.1 * u) {
  gr <- sdeprobe_grid(2, 1 / 2048, 2150)
  prob <- SDEProblem(function(u, p, t) -(1 + 8000 * t * t * (1 - t) * (1 - t)) * u, g, c(1, 0.5), c(0, 1),
                     noise = NoiseGrid(gr$t, gr$W, gr$Z))
  sol <- solve(prob, alg, dt = 0.01)
  r <- rle(as.integer(sol$alg_choice))
  list(n = length(sol$t), u = c(sol$u[[100]], sol$u[[300]], sol$u[[length(sol$t)]]), t = sol$t[c(100, 300)],
       choice = rbind(r$values, r$lengths))
}
sdecomp_u300 <- c(0.15445029487845618, 0.07745916423711882)

test_that("AutoAlgSwitch matches StochasticDiffEq", {
  r <- sdecomp_run(AutoAlgSwitch(SOSRI2(), SOSRI(), maxstiffstep = 2, maxnonstiffstep = 1, dtfac = 1))
  expect_equal(r$n, 431L)
  expect_equal(r$choice, rbind(c(1, 2, 1, 2, 1), c(407, 7, 3, 2, 12)), ignore_attr = TRUE)
  expect_equal(r$t, c(0.05304395792194349, 0.0915392254843876), tolerance = 1e-14)
  expect_equal(r$u[3:4], sdecomp_u300, tolerance = 1e-12)
  expect_equal(r$u[5:6], c(1.2951914939328671e-5, 4.898867039117482e-6), tolerance = 1e-10)
  # dtfac = 2: the reference with the noise following the grown step
  r <- sdecomp_run(AutoAlgSwitch(SOSRI2(), SOSRI(), maxstiffstep = 2, maxnonstiffstep = 1))
  expect_equal(r$n, 433L)
  expect_equal(r$choice, rbind(c(1, 2, 1, 2, 1), c(407, 2, 11, 3, 10)), ignore_attr = TRUE)
  expect_equal(r$u[5:6], c(1.4874387222562922e-5, 4.695866477551159e-6), tolerance = 1e-10)
  r <- sdecomp_run(AutoAlgSwitch(SOSRI2(), SOSRI(), maxstiffstep = 2, maxnonstiffstep = 1, stiffalgfirst = TRUE))
  expect_equal(r$n, 441L)
  expect_equal(r$choice, rbind(c(2, 1, 2, 1, 2, 1), c(1, 414, 2, 12, 3, 9)), ignore_attr = TRUE)
  expect_equal(r$t, c(0.05139346224195675, 0.08980637477764176), tolerance = 1e-13)
  expect_equal(r$u[5:6], c(2.7063443908654387e-6, -1.0766188027437436e-7), tolerance = 1e-7)
})

test_that("AutoSOSRI2 and AutoSOSRA2 match StochasticDiffEq", {
  r <- sdecomp_run(AutoSOSRI2(LambaEM(split = FALSE), maxstiffstep = 1, maxnonstiffstep = 1, dtfac = 1))
  expect_equal(r$n, 433L)
  expect_equal(r$choice, rbind(c(1, 2, 1, 2), c(406, 4, 13, 10)), ignore_attr = TRUE)
  expect_equal(r$u[5:6], c(-0.00020130819935125984, -0.0009900470750266485), tolerance = 1e-10)
  r <- sdecomp_run(AutoSOSRA2(SOSRA(), maxstiffstep = 2, maxnonstiffstep = 1, nonstifftol = 0.2, stifftol = 0.2,
                              dtfac = 1), g = function(u, p, t) c(0.1, 0.2))
  expect_equal(r$n, 677L)
  expect_equal(r$choice, rbind(c(1, 2, 1, 2, 1, 2, 1), c(452, 2, 49, 2, 67, 7, 98)), ignore_attr = TRUE)
  expect_equal(r$t, c(0.056576091910853755, 0.10260799825349254), tolerance = 1e-14)
  expect_equal(r$u, c(0.6075051856052301, 0.30501987320895785, 0.07719219716347196, 0.0370544868254499,
                      -0.00010619699301084404, -0.0021875402861042522), tolerance = 1e-12)
})

test_that("StochasticCompositeAlgorithm with a choice function matches StochasticDiffEq", {
  r <- sdecomp_run(StochasticCompositeAlgorithm(list(SOSRI2(), SOSRI()), function(x) if (x$t > 0.5) 2 else 1))
  expect_equal(r$n, 430L)
  expect_equal(r$choice, rbind(c(1, 2), c(415, 15)), ignore_attr = TRUE)
  expect_equal(r$u[5:6], c(-1.2901076611589427e-7, 2.9266750541671574e-7), tolerance = 1e-10)
  expect_error(StochasticCompositeAlgorithm(list(SOSRI2(), Tsit5()), function(x) 1), "not an SDE algorithm")
  expect_error(StochasticCompositeAlgorithm(list(SOSRI2()), 3), "AutoSwitch")
})

test_that("stability region sizes and orders of the switching methods", {
  expect_equal(alg_stability_size(SOSRI2()), 10.6)
  expect_equal(alg_stability_size(SOSRA2()), 5.3)
  expect_equal(SOSRA()$order, 2)
  a <- AutoSOSRI2(SOSRI())
  expect_s3_class(a, "StochasticCompositeAlgorithm")
  expect_equal(a$order, max(SOSRI2()$order, SOSRI()$order))
})
