test_that("isoutofdomain rejects steps, unstable_check ends, internalnorm and force_dtmin are used", {
  pr <- ODEProblem(function(u, p, t) -1, 1, c(0, 2))
  expect_equal(tail(solve(pr, Tsit5())$u, 1)[[1]], -1, tolerance = 1e-12)
  # the state cannot leave u >= 0, so the steps shrink at t = 1 until dt
  # falls below eps(t)
  expect_warning(s <- solve(pr, Tsit5(), isoutofdomain = function(u, p, t) any(u < 0)), "eps\\(t\\)")
  expect_equal(s$retcode, "Unstable")
  expect_true(all(unlist(s$u) >= 0))
  expect_gt(s$stats$nreject, 0L)
  expect_warning(s <- solve(ODEProblem(function(u, p, t) u, 1, c(0, 10)), Tsit5(),
                            unstable_check = function(dt, u, p, t) any(abs(u) > 100)), "unstable_check")
  expect_equal(s$retcode, "Unstable")
  expect_gt(tail(s$t, 1), log(100))
  expect_lt(tail(s$t, 1), 10)
  calls <- 0L
  nrm <- function(u, t) {
    calls <<- calls + 1L
    max(abs(u))
  }
  dp <- ODEProblem(function(u, p, t) -u, c(1, 2), c(0, 1))
  s <- solve(dp, Tsit5(), internalnorm = nrm)
  expect_gt(calls, 0L)
  expect_equal(s$retcode, "Success")
  expect_false(identical(s$t, solve(dp, Tsit5())$t))
  # force_dtmin keeps stepping at dtmin instead of stopping with DtLessThanMin
  bp <- ODEProblem(function(u, p, t) u^2, 1, c(0, 2))
  expect_warning(s <- solve(bp, Tsit5(), dtmin = 1e-10), "dtmin")
  expect_equal(s$retcode, "DtLessThanMin")
  s <- suppressWarnings(solve(bp, Tsit5(), dtmin = 1e-10, force_dtmin = TRUE))
  expect_false(identical(s$retcode, "DtLessThanMin"))
})

test_that("dt = 0 is the automatic step when adaptive and steps between tstops when fixed", {
  pr <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
  expect_identical(solve(pr, Tsit5(), dt = 0)$t, solve(pr, Tsit5())$t)
  ep <- ODEProblem(function(u, p, t) 1, 0, c(0, 1))
  expect_error(solve(ep, Euler(), dt = 0), "InexactError")
  s <- solve(ep, Euler(), dt = 0, saveat = c(0, 1), tstops = 0.5)
  expect_equal(s$t, c(0, 1))
  expect_equal(unlist(s$u), c(0, 1))
})

test_that("infinite time spans", {
  # test/InterfaceI/inf_handling.jl
  pr <- ODEProblem(function(u, p, t) 1.01 * u, 1 / 2, c(0, Inf))
  s <- suppressWarnings(solve(pr, Tsit5(), reltol = 1e-8, abstol = 1e-8, adaptive = FALSE, dt = 0.1))
  expect_equal(s$retcode, "Unstable")
  dp <- ODEProblem(function(u, p, t) -u, 1, c(0, Inf))
  expect_warning(s <- solve(dp, Tsit5(), maxiters = 200), "maximum")
  expect_equal(s$retcode, "MaxIters")
  cb <- DiscreteCallback(function(u, t, integ) t >= 5, function(integ) terminate(integ))
  s <- solve(dp, Tsit5(), callback = cb)
  expect_equal(s$retcode, "Terminated")
  expect_equal(as.numeric(tail(s$u, 1)[[1]]), exp(-tail(s$t, 1)), tolerance = 1e-3)
  expect_error(solve(dp, Tsit5(), saveat = 1), "finite time span")
})

test_that("backward solutions interpolate linearly between saved points", {
  s <- solve(ODEProblem(function(u, p, t) 1, 1, c(1, 0)), Euler(), dt = 0.5, saveat = c(1, 0.5, 0))
  expect_equal(s(0.75), 0.75)
  expect_equal(s(0.25), 0.25)
  expect_equal(as.numeric(s(c(0.1, 0.9))), c(0.1, 0.9))
  expect_equal(s(0.75, deriv = 1), 1)
})

test_that("save_positions[1] is honored when saving every step", {
  cb <- PresetTimeCallback(0.5, function(integ) integ$u[1] <- integ$u[1] + 1, save_positions = c(FALSE, TRUE))
  s <- solve(ODEProblem(function(u, p, t) 0, 1, c(0, 1)), Euler(), dt = 0.25, callback = cb)
  expect_equal(sum(s$t == 0.5), 1L)
  expect_equal(s[1, which(s$t == 0.5)], 2)
  cb2 <- PresetTimeCallback(0.5, function(integ) integ$u[1] <- integ$u[1] + 1)
  s2 <- solve(ODEProblem(function(u, p, t) 0, 1, c(0, 1)), Euler(), dt = 0.25, callback = cb2)
  expect_equal(sum(s2$t == 0.5), 2L)
})

test_that("derivatives of the interpolants", {
  pr <- ODEProblem(function(u, p, t) c(u[2], -u[1]), c(0, 1), c(0, 2))
  tt <- 1.3
  ders <- list(c(cos(tt), -sin(tt)), c(-sin(tt), -cos(tt)), c(-cos(tt), sin(tt)), c(sin(tt), cos(tt)))
  # the highest order of each interpolant and the tolerance of its derivatives
  cases <- list(list(Tsit5(), 4L), list(DP5(), 4L), list(Vern7(), 1L), list(Vern9(), 4L), list(BS3(), 3L),
                list(Rodas5P(), 3L), list(Rosenbrock23(), 1L), list(FBDF(), 1L), list(QNDF(), 1L),
                list(RK4(), 3L), list(SSPRK33(), 2L))
  for (cs in cases) {
    alg <- cs[[1L]]
    s <- solve(pr, alg, reltol = 1e-10, abstol = 1e-12, dt = 0.01)
    expect_equal(unname(s(tt, deriv = 1)), ders[[1L]], tolerance = 1e-3, label = alg$name)
    expect_error(s(tt, deriv = cs[[2L]] + 1L), "Derivative order too high", label = alg$name)
  }
  s <- solve(pr, Vern9(), reltol = 1e-12, abstol = 1e-14)
  for (n in 1:4) expect_equal(unname(s(tt, deriv = n)), ders[[n]], tolerance = 1e-6)
  s <- solve(pr, Tsit5(), reltol = 1e-10, abstol = 1e-12)
  for (n in 1:3) expect_equal(unname(s(tt, deriv = n)), ders[[n]], tolerance = 1e-4)
  expect_equal(unname(s(c(0.5, tt), deriv = 1)), rbind(c(cos(0.5), -sin(0.5)), ders[[1L]]), tolerance = 1e-6)
  expect_equal(s(tt, idxs = 2, deriv = 1), -sin(tt), tolerance = 1e-6)
  # linear interpolation: the slope, nothing higher
  s <- solve(ODEProblem(function(u, p, t) 2 * t, 0, c(0, 1)), Tsit5(), saveat = c(0, 0.5, 1))
  expect_equal(s(0.25, deriv = 1), 0.5)
  expect_error(s(0.25, deriv = 2), "Derivative order too high")
})

test_that("get_du is the derivative for mass matrix problems", {
  mp <- ODEProblem(ODEFunction(function(u, p, t) 2, mass_matrix = matrix(2)), 0, c(0, 1))
  integ <- init(mp, Rodas5P())
  step_integrator(integ)
  expect_equal(get_du(integ), 1, tolerance = 1e-10)
  integ <- init(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), Vern7())
  step_integrator(integ)
  expect_equal(get_du(integ), -integ$u, tolerance = 1e-8)
})

test_that("ImplicitDiscreteProblem takes a positional p and solution builders copy", {
  ip <- ImplicitDiscreteProblem(function(un, u, p, t) un - u - p[1], 0, c(0, 1), 2)
  expect_equal(ip$p, 2)
  expect_length(ip$kwargs, 0L)
  s <- solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), Tsit5())
  s2 <- solution_new_retcode(s, "Failure")
  expect_equal(s$retcode, "Success")
  expect_equal(s2$retcode, "Failure")
  s3 <- sensitivity_solution(s, list(5), 0)
  expect_equal(length(s3$t), 1L)
  expect_gt(length(s$t), 1L)
  expect_true(isTRUE(s$dense))
  expect_false(isTRUE(strip_interpolation(s)$dense))
  expect_true(isTRUE(s$dense))
})

test_that("empty systems, right-hand side lengths, initial values and keywords are checked", {
  e <- ODEProblem(function(u, p, t) numeric(), numeric(), c(0, 1))
  for (alg in list(Tsit5(), Rodas5P(), FBDF(), RadauIIA5())) {
    s <- solve(e, alg)
    expect_equal(s$retcode, "Success", label = alg$name)
    expect_equal(s$t, c(0, 1), label = alg$name)
    expect_length(s$u[[1L]], 0L)
  }
  expect_equal(solve(e, Tsit5(), saveat = 0.25)$t, seq(0, 1, by = 0.25))
  expect_error(solve(ODEProblem(function(u, p, t) -u[1], c(1, 2), c(0, 1)), Tsit5()), "returned 1 values")
  expect_error(solve(ODEProblem(function(u, p, t) c(-u, 1), 1, c(0, 1)), Tsit5()), "returned 2 values")
  expect_error(solve(ODEProblem(function(u, p, t) -u, c(1, NaN), c(0, 1)), Tsit5()),
               "First function call produced NaNs.*not finite: 2")
  expect_error(solve(ODEProblem(function(u, p, t) -u, 1, c(0, NA)), Tsit5()), "without NA")
  pr <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
  expect_error(solve(pr, Tsit5(), relto = 1e-12), class = "CommonKwargError")
  expect_error(init(pr, Tsit5(), relto = 1e-12), class = "CommonKwargError")
  expect_warning(s <- solve(pr, Tsit5(), relto = 1e-12, kwargshandle = "warn"), "relto")
  expect_equal(s$retcode, "Success")
  expect_silent(solve(pr, Tsit5(), relto = 1e-12, kwargshandle = "silent"))
})
