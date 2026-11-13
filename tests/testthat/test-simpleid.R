# Reference values from tools/simpleimplicitdiscretesolve-probe.jl
# (SimpleImplicitDiscreteSolve 2.0.4, SimpleNonlinearSolve).

sids_lv <- function(un, u, p, t) {
  h <- 0.1
  c(un[1] - u[1] - h * (p[1] * un[1] - p[2] * un[1] * un[2]),
    un[2] - u[2] - h * (-p[3] * un[2] + p[4] * un[1] * un[2] + t * t))
}
sids_lv_u <- c(1.0708606205651852, 0.8382834403996867, 1.1628383141027294, 0.7090242261365783,
               1.2766923606110165, 0.6082107951691874, 1.4134125787068423, 0.5326942029841446,
               1.5739189318173195, 0.4802120689586657, 1.7586817798161263, 0.4494242214863676,
               1.9672092740867537, 0.4399831679451949, 2.19734831841649, 0.45265105490610164,
               2.4443620886944655, 0.48945507533253535, 2.6998101114602324, 0.5538296686813872,
               2.9504174455276457, 0.6506038088024005)

test_that("SimpleIDSolve matches SimpleImplicitDiscreteSolve", {
  cubic <- function(un, u, p, t) un[1] + un[1] * un[1] * un[1] - u[1]
  sol <- solve(ImplicitDiscreteProblem(cubic, 2, c(0, 3)), SimpleIDSolve())
  expect_identical(sol$retcode, "Default")
  expect_equal(sol$t, c(0, 1, 2, 3))
  expect_equal(unlist(sol$u), c(1, 0.6823278038280193, 0.5318696690666097, 0.44421441630236674),
               tolerance = 1e-14)

  prob <- ImplicitDiscreteProblem(sids_lv, c(1, 1), c(0, 1), p = c(1.5, 1, 3, 1))
  sol <- solve(prob, SimpleIDSolve(), dt = 0.1)
  # the times of the Julia range 0:0.1:1 are k / 10 correctly rounded
  expect_identical(sol$t, (0:10) / 10)
  expect_equal(unlist(sol$u), sids_lv_u, tolerance = 1e-14)

  # the initial solve f(u, u0, p, t0) = 0 replaces u0
  alg2 <- function(un, u, p, t) c(un[1] - 0.5 * u[1], un[2] - un[1] * un[1])
  sol <- solve(ImplicitDiscreteProblem(alg2, c(1, 5), c(0, 2)), SimpleIDSolve())
  expect_equal(unlist(sol$u), c(0.5, 0.25, 0.25, 0.0625, 0.125, 0.015625))
})

test_that("SimpleIDSolve reports failed solves as SimpleImplicitDiscreteSolve", {
  noroot <- function(un, u, p, t) un[1] * un[1] + 1
  sol <- solve(ImplicitDiscreteProblem(noroot, 1.3, c(0, 2)), SimpleIDSolve())
  expect_identical(sol$retcode, "InitialFailure")
  expect_equal(sol$t, 0)
  expect_equal(unlist(sol$u), 1.3)
  # the failed step ends the solution with its last Newton iterate
  shrink <- function(un, u, p, t) un[1] * un[1] - (u[1] - 1)
  sol <- solve(ImplicitDiscreteProblem(shrink, 5, c(0, 5)), SimpleIDSolve())
  expect_identical(sol$retcode, "ConvergenceFailure")
  expect_equal(sol$t, c(0, 1, 2, 3))
  expect_equal(unlist(sol$u), c(2, 1, 2.384185791015625e-7, 10.28111610490308), tolerance = 1e-14)
})

test_that("SimpleIDSolve save options", {
  prob <- ImplicitDiscreteProblem(sids_lv, c(1, 1), c(0, 1), p = c(1.5, 1, 3, 1))
  n <- length(sids_lv_u)
  # Julia returns all 11 times with 2 states here; the saved times are used
  sol <- solve(prob, SimpleIDSolve(), dt = 0.1, save_everystep = FALSE)
  expect_equal(sol$t, c(0, 1))
  expect_equal(unlist(sol$u), sids_lv_u[c(1:2, n - 1, n)], tolerance = 1e-14)
  sol <- solve(prob, SimpleIDSolve(), dt = 0.1, save_everystep = FALSE, save_start = FALSE)
  expect_equal(sol$t, 1)
  expect_equal(unlist(sol$u), sids_lv_u[(n - 1):n], tolerance = 1e-14)
  # Julia raises a BoundsError here
  sol <- solve(prob, SimpleIDSolve(), dt = 0.1, save_start = FALSE)
  expect_identical(sol$t, (1:10) / 10)
  expect_equal(unlist(sol$u), sids_lv_u[-(1:2)], tolerance = 1e-14)
  expect_error(solve(prob, SimpleIDSolve(), adaptive = TRUE), "adaptive")
})

test_that("Julia float ranges", {
  expect_identical(.sids_float_range(0, 0.1, 1), (0:10) / 10)
  expect_identical(.sids_float_range(0, 0.25, 1), c(0, 0.25, 0.5, 0.75, 1))
  expect_identical(.sids_float_range(1, 1, 0), numeric())
  expect_identical(.sids_float_range(0, pi, 7), c(0, pi, 2 * pi))
})
