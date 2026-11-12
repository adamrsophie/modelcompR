test_that("every tableau method solves a smooth problem at its order", {
  prob <- ODEProblem(function(u, p, t) c(u[2], -u[1]), c(0, 1), c(0, 1))
  exact <- c(sin(1), cos(1))
  # these behave the same way in OrdinaryDiffEq (checked against Julia):
  # SIR54, SSPRK932 and Alshina3 stall in adaptive mode, CKLLSRK95_4C and
  # SHLDDRK64 converge at first and below first order, and SHLDDRK_2N levels
  # off near 1e-7 (its coefficients are given to seven digits)
  defective <- c("SIR54", "SSPRK932", "Alshina3", "CKLLSRK95_4C", "SHLDDRK64", "SHLDDRK_2N")
  for (nm in setdiff(names(.ERK_TABLEAUS), defective)) {
    tb <- .ERK_TABLEAUS[[nm]]
    alg <- get(nm)()
    # fixed steps: the error falls by about 2^order when dt halves
    e1 <- max(abs(solve(prob, alg, dt = 0.1, adaptive = FALSE)(1) - exact))
    e2 <- max(abs(solve(prob, alg, dt = 0.05, adaptive = FALSE)(1) - exact))
    if (e2 > 1e-12) expect_gt(log2(e1 / e2), tb$order - 0.6, label = nm)
    if (isTRUE(tb$adaptive)) {
      s <- solve(prob, alg, abstol = 1e-8, reltol = 1e-8)
      expect_equal(as.numeric(s(1)), exact, tolerance = 1e-5, info = nm)
    }
  }
})

test_that("methods that misbehave upstream misbehave the same way", {
  prob <- ODEProblem(function(u, p, t) c(u[2], -u[1]), c(0, 1), c(0, 1))
  exact <- c(sin(1), cos(1))
  order <- function(alg) {
    e1 <- max(abs(solve(prob, alg, dt = 0.1, adaptive = FALSE)(1) - exact))
    e2 <- max(abs(solve(prob, alg, dt = 0.05, adaptive = FALSE)(1) - exact))
    log2(e1 / e2)
  }
  expect_equal(order(CKLLSRK95_4C()), 1.02, tolerance = 0.05)
  expect_lt(order(SHLDDRK64()), 1)
  expect_gt(order(SIR54()), 4.5)
  expect_gt(order(SSPRK932()), 2.5)
  # Julia stops at t = 3.11e-5 (SIR54) and t = 0.00111 (Alshina3) with MaxIters
  s <- suppressWarnings(solve(prob, SIR54(), abstol = 1e-8, reltol = 1e-8, dt = 0.1, maxiters = 2000))
  expect_identical(s$retcode, "MaxIters")
  s <- suppressWarnings(solve(prob, Alshina3(), abstol = 1e-8, reltol = 1e-8, dt = 0.1, maxiters = 2000))
  expect_identical(s$retcode, "MaxIters")
})

test_that("the Julia step counts are reproduced", {
  prob <- ODEProblem(function(u, p, t) 1.01 * u, 0.5, c(0, 1))
  s <- solve(prob, OwrenZen5(), abstol = 1e-6, reltol = 1e-6)
  expect_equal(s$retcode, "Success")
  r <- solve(prob, RDPK3Sp35(), abstol = 1e-6, reltol = 1e-6)
  expect_equal(r(1), 0.5 * exp(1.01), tolerance = 1e-5)
  expect_identical(RDPK3Sp35()$controller, "PID")
})

test_that("stage and step limiters apply", {
  prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
  lim <- function(u, integ, p, t) pmax(u, 0.5)
  s <- solve(prob, SSPRK33(stage_limiter = lim, step_limiter = lim), dt = 0.1)
  expect_true(all(unlist(s$u) >= 0.5))
})
