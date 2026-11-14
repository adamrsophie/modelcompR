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

# References: tools/nlfunctional-probe.jl. The step counts, nsolve and the
# end state match; nf of the Newton runs differs by the f calls of the
# finite difference Jacobians (OrdinaryDiffEq uses ForwardDiff).
nlf_prob <- ODEProblem(function(u, p, t) c(-2 * u[1] + u[2], u[1] - 3 * u[2] + sin(t)), c(1, 0.5), c(0, 2))

test_that("NLFunctional and NLAnderson match OrdinaryDiffEq", {
  # method, solver, nf, naccept, nreject, nsolve, u(2)
  ref <- list(
    list("ImplicitEuler", "fun", 660, 88, 37, 0, c(0.2111945412602614, 0.3918884144162806)),
    list("ImplicitEuler", "and", 527, 88, 38, 146, c(0.21117105133311773, 0.3921375907097098)),
    list("Trapezoid", "fun", 181, 23, 4, 0, c(0.20800168287389625, 0.3915240644060963)),
    list("Trapezoid", "and", 120, 23, 4, 36, c(0.2080042297053185, 0.391527385273913)),
    list("TRBDF2", "fun", 96, 17, 0, 0, c(0.20843057461961587, 0.3917989941055978)),
    list("TRBDF2", "and", 85, 17, 0, 15, c(0.20843107700789817, 0.39179861809637817)),
    list("KenCarp4", "fun", 136, 8, 0, 0, c(0.20889395341340877, 0.3920953256930352)),
    list("KenCarp4", "and", 111, 8, 0, 28, c(0.20889579955938836, 0.3920928592773954)),
    list("SDIRK2", "fun", 502, 78, 0, 0, c(0.20890068319330796, 0.39198403736846704)),
    list("SDIRK2", "and", 470, 78, 0, 77, c(0.20890464050482474, 0.3919856949344779)))
  for (r in ref) {
    nls <- if (r[[2]] == "fun") NLFunctional() else NLAnderson()
    sol <- solve(nlf_prob, get(r[[1]])(nlsolve = nls))
    info <- paste(r[[1]], r[[2]])
    expect_equal(sol$retcode, "Success", info = info)
    expect_equal(unlist(sol$stats[c("nf", "naccept", "nreject", "nsolve", "njacs")]),
                 c(nf = r[[3]], naccept = r[[4]], nreject = r[[5]], nsolve = r[[6]], njacs = 0), info = info)
    expect_equal(sol$u[[length(sol$u)]], r[[7]], tolerance = 1e-13, info = info)
  }
})

test_that("NLNewton relax and always_new match OrdinaryDiffEq", {
  # method, naccept, nreject, njacs, nsolve, nf less the Jacobian calls, u(2)
  relax <- list(
    list("ImplicitEuler", 114, 51, 165, 1026, 1359, c(0.21120475723512905, 0.3919037304155836)),
    list("Trapezoid", 24, 4, 27, 199, 257, c(0.20800163058978657, 0.39152278857950396)),
    list("TRBDF2", 10, 0, 8, 128, 129, c(0.2074435716793809, 0.3911630612841056)),
    list("KenCarp4", 6, 0, 10, 161, 168, c(0.20894150615904092, 0.3921482985433535)),
    list("SDIRK2", 25, 0, 13, 393, 409, c(0.20896999110322007, 0.39143687360437884)))
  # always_new against the in place form of OrdinaryDiffEq (the out of place
  # form overwrites uprev, see PORTING.md); Trapezoid takes a different step
  # sequence in place even with the default Newton settings
  anew <- list(
    list("ImplicitEuler", 88, 38, 247, 247, 623, c(0.2111651341097146, 0.3921287720779027)),
    list("TRBDF2", 8, 3, 43, 54, 89, c(0.2059260916623604, 0.39020979291885294)),
    list("KenCarp4", 6, 0, 42, 48, 87, c(0.20893730488503143, 0.3921442581559428)),
    list("SDIRK2", 26, 0, 61, 87, 151, c(0.20896488947232777, 0.3914506617408477)))
  for (case in list(list(relax, NLNewton(relax = 0.3)), list(anew, NLNewton(always_new = TRUE)))) {
    for (r in case[[1]]) {
      sol <- solve(nlf_prob, get(r[[1]])(nlsolve = case[[2]]))
      st <- sol$stats
      info <- r[[1]]
      expect_equal(c(st$naccept, st$nreject, st$njacs, st$nsolve, st$nf - 2L * st$njacs),
                   unlist(r[2:6]), info = info)
      expect_equal(sol$u[[length(sol$u)]], r[[7]], tolerance = 1e-10, info = info)
    }
  }
})

test_that("a line search relaxation scales the Newton step", {
  seen <- 0L
  ls <- function(phi, dphi, phidphi, a0, phi0, dphi0) {
    seen <<- seen + 1L
    list(a0, phi(a0))
  }
  sol <- solve(nlf_prob, ImplicitEuler(nlsolve = NLNewton(relax = ls)))
  ref <- solve(nlf_prob, ImplicitEuler())
  expect_gt(seen, 0L)
  expect_equal(sol$u[[length(sol$u)]], ref$u[[length(ref$u)]], tolerance = 1e-14)
})

test_that("unsupported nonlinear solvers are refused", {
  expect_error(KenCarp4(nlsolve = NonlinearSolveAlg()), "NLNewton")
  expect_error(ImplicitEuler(nlsolve = NLNewton(relax = "a")), "relax")
})
