# Dense output of the methods without an interpolant of their own: the free
# SSP interpolant of SSPRK22, SSPRK33, SSPRK43 and SSPRK432 and the generic
# Hermite output of OrdinaryDiffEqCore (R/solvers-erk.R .ERK_DENSE). The
# dense41_* keys of tools/globaldiffeq-probe.jl record the Julia values on the
# same problem (pRRK33 with kappa = 0.5).

dno_lv <- function(u, p, t) c(1.5 * u[1] - u[1] * u[2], -3 * u[2] + u[1] * u[2] + 0.1 * t)
dno_prob <- function() ODEProblem(dno_lv, c(1, 1), c(0, 1))

test_that("the dense derivative is the derivative of the dense value", {
  for (nm in c("SSPRK33", "SSPRK43", "SSPRK53", "SSPRK932", "Alshina3", "PSRK3p5q4", "RKO65", "RKV76IIa",
               "pRRK33", "RK46NL")) {
    alg <- if (nm == "pRRK33") pRRK33(kappa = 0.5) else get(nm)()
    sol <- solve(dno_prob(), alg, dt = 0.1, adaptive = FALSE)
    dn <- .gde_sol_dense(sol, right = FALSE)
    for (t in c(0.03, 0.47)) {
      h <- 1e-4
      fd <- (8 * (dn$value(t + h) - dn$value(t - h)) - (dn$value(t + 2 * h) - dn$value(t - 2 * h))) / (12 * h)
      expect_equal(dn$deriv(t), fd, tolerance = 1e-8, info = nm)
      expect_equal(sol(t), dn$value(t), tolerance = 1e-15, info = nm)
    }
  }
})

test_that("the SSP interpolant is the free second order one", {
  y0 <- c(1, 2); y1 <- c(1.3, 1.7); k <- list(c(2, -4)); dt <- 0.2
  expect_identical(.ssprk_interp(0, dt, y0, y1, k), y0)
  expect_identical(.ssprk_interp(1, dt, y0, y1, k), y1)
  expect_equal(.gde_ssprk_d(0, dt, y0, y1, k), k[[1L]])
  expect_equal(.gde_ssprk_d(1, dt, y0, y1, k), 2 * (y1 - y0) / dt - k[[1L]])
})

test_that("the stage Hermite methods take their first two stages as end slopes", {
  # OrdinaryDiffEq has no interpolant for these and its generic Hermite
  # interpolant reads k[1] and k[2], the first two stages
  sol <- solve(dno_prob(), RKO65(), dt = 0.1, adaptive = FALSE)
  st <- .sol_env(sol)$steps[[1L]]
  expect_equal(st$k[[1L]], dno_lv(c(1, 1), NULL, 0.1 * 2 / 3))
  sol <- solve(dno_prob(), Alshina2(), dt = 0.1, adaptive = FALSE)
  st <- .sol_env(sol)$steps[[1L]]
  expect_equal(st$k[[1L]], dno_lv(c(1, 1), NULL, 0))
})

test_that("DImplicitEuler keeps the derivatives of the last two steps", {
  F <- function(du, u, p, t) c(du[1] - u[2], u[2]^3 + u[2] + u[1]^3 + u[1])
  prob <- DAEProblem(F, c(-1, -1), c(1, -1), c(0, 2), differential_vars = c(TRUE, FALSE))
  sol <- solve(prob, DImplicitEuler(), reltol = 1e-4, abstol = 1e-6, initializealg = NoInit())
  steps <- .sol_env(sol)$steps
  expect_equal(steps[[1L]]$k[[2L]], c(-1, -1))
  expect_true(.gde_has_dense_derivative(sol))
  # the algebraic component is interpolated linearly
  t <- (sol$t[3] + sol$t[4]) / 2
  expect_equal(sol(t)[2], (sol$u[[3]][2] + sol$u[[4]][2]) / 2, tolerance = 1e-14)
})
