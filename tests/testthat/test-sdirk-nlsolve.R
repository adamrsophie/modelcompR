# The nlsolve keyword of the SDIRK methods, ImplicitEuler and Trapezoid
# (NLNewton settings and a callable stage predictor). The Julia references
# are the _nlk, _nlpred and _nlnodiv cases of test-sdirk.R.

nls_prob <- ODEProblem(function(u, p, t) -50 * (u - cos(t)), c(0, 1), c(0, 1))

test_that("NLNewton settings reach the Newton iteration", {
  ref <- solve(nls_prob, Kvaerno4())
  tight <- solve(nls_prob, Kvaerno4(nlsolve = NLNewton(kappa = 1e-6)))
  expect_gt(tight$stats$nsolve, ref$stats$nsolve)
  expect_equal(tight(1), ref(1), tolerance = 1e-3)
  expect_identical(Kvaerno4(nlsolve = NLNewton(max_iter = 3L))$nlsolve$max_iter, 3L)
  for (a in list(ImplicitEuler(nlsolve = NLNewton(kappa = 1e-6)), Trapezoid(nlsolve = NLNewton(kappa = 1e-6)),
                 SDIRK22(nlsolve = NLNewton(kappa = 1e-6))))
    expect_equal(solve(nls_prob, a)$retcode, "Success")
})

test_that("a stage predictor seeds every stage solve", {
  seen <- numeric()
  pr <- function(uprev, p, t, dt) {
    seen <<- c(seen, t)
    uprev
  }
  sol <- solve(nls_prob, SFSDIRK4(nlsolve = NLNewton(predictor = pr)), dt = 0.25, adaptive = FALSE)
  # four implicit stages per step, at t + c dt
  expect_length(seen, 16L)
  expect_equal(sol(1), solve(nls_prob, SFSDIRK4(), dt = 0.25, adaptive = FALSE)(1), tolerance = 1e-6)
  ex <- solve(nls_prob, ImplicitEuler(nlsolve = NLNewton(predictor = pr)), dt = 0.1, adaptive = FALSE)
  expect_equal(ex$retcode, "Success")
})

test_that("unsupported nonlinear solvers are refused", {
  expect_error(KenCarp4(nlsolve = NLFunctional()), "NLNewton")
  expect_error(TRBDF2(nlsolve = NLNewton(always_new = TRUE)), "always_new")
  expect_error(ImplicitEuler(nlsolve = NLNewton(relax = 0.5)), "relax")
})
