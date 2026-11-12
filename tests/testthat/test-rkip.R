# Reference values from tools/rkip-probe.jl (Julia, OrdinaryDiffEqRKIP 2.1.3): final
# states and statistics (nf, nf2, naccept, nreject).
.RKIP_REF <- list(
  `default` = list(u = c(0.026315036106357024, 0.02313801431137406, -0.051789441624872048, -0.0018203776513146144),
    stats = c(122, 64, 7, 0), nt = 8, tend = 2),
  `tight` = list(u = c(0.026315050704621229, 0.023138033522611585, -0.051789430965066642, -0.0018203854884296265),
    stats = c(360, 190, 21, 0), nt = 22, tend = 2),
  `gridb` = list(u = c(0.026315050438202311, 0.023138033147791972, -0.051789431374891011, -0.0018203852719540876),
    stats = c(224, 118, 13, 0), nt = 14, tend = 2),
  `clamplow` = list(u = c(0.026315050716799553, 0.023138033540241576, -0.051789430946544118, -0.0018203854951215684),
    stats = c(479, 253, 28, 0), nt = 29, tend = 2),
  `clamphigh` = list(u = c(-1.9305807313581636e-05, -3.6765897628307673e-05, -4.870645888050543e-05,
    -2.3006046334184908e-06),
    stats = c(3403, 1801, 200, 0), nt = 201, tend = 10),
  `fixed` = list(u = c(0.026315049894676687, 0.023138032014746096, -0.051789431747957658, -0.0018203855398378347),
    stats = c(341, 181, 20, 0), nt = 21, tend = 2),
  `fixedoff` = list(u = c(0.026315042252592458, 0.023138018150651987, -0.051789438899684279, -0.0018203859337150112),
    stats = c(239, 127, 14, 0), nt = 15, tend = 2),
  `reuse1` = list(u = c(0.026315049548181455, 0.023138031984057256, -0.051789432157504485, -0.0018203852472352376),
    stats = c(173, 91, 10, 0), nt = 11, tend = 2),
  `reuse2` = list(u = c(0.0035372146588893145, 0.0006051268741398382, -0.018860924101489405, -0.00074477619423126665),
    stats = c(258, 136, 15, 0), nt = 16, tend = 3)
)

rkip_A4 <- matrix(c(-2, 0.5, 0.1, 0, 0.3, -1.5, 0.4, 0.2, 0, 0.2, -1, 0.5, 0.1, 0, 0.3, -8), 4, byrow = TRUE)
rkip_f2 <- function(u, p, t) c(-u[1] * u[1] * u[1] + 0.2 * u[2] * u[3], 0.5 * u[1] * u[2],
                               -0.2 * u[3] * u[3] + 0.3 * u[4], 0.1 * u[4] * u[1])
rkip_prob <- function(tspan) SplitODEProblem(MatrixOperator(rkip_A4), rkip_f2, c(1, 0.5, -0.5, 0.2), tspan)

rkip_check <- function(sol, nm, tol = 1e-12) {
  ref <- .RKIP_REF[[nm]]
  expect_length(sol$u, ref$nt)
  expect_equal(sol$t[[ref$nt]], ref$tend, tolerance = 1e-14, label = nm)
  expect_equal(as.numeric(sol$u[[ref$nt]]), ref$u, tolerance = tol, label = nm)
  expect_equal(c(sol$stats$nf, sol$stats$nf2, sol$stats$naccept, sol$stats$nreject), ref$stats, label = nm)
}

test_that("RKIP matches OrdinaryDiffEqRKIP", {
  prob <- rkip_prob(c(0, 2))
  rkip_check(solve(prob, RKIP()), "default")
  rkip_check(solve(prob, RKIP(), reltol = 1e-8, abstol = 1e-10), "tight")
  rkip_check(solve(prob, RKIP(1e-2, 10), reltol = 1e-6, abstol = 1e-8), "gridb")
  rkip_check(solve(prob, RKIP(clamp_lower_dt = TRUE), reltol = 1e-9, abstol = 1e-11), "clamplow")
  rkip_check(solve(rkip_prob(c(0, 10)), RKIP(1e-4, 0.05, nb_of_cache_step = 37)), "clamphigh")
  rkip_check(solve(prob, RKIP(), dt = 0.1, adaptive = FALSE), "fixed")
  rkip_check(solve(prob, RKIP(), dt = 0.15, adaptive = FALSE), "fixedoff")
})

test_that("RKIP keeps its cached exponentials across solves", {
  alg <- RKIP()
  rkip_check(solve(rkip_prob(c(0, 2)), alg, reltol = 1e-5, abstol = 1e-7), "reuse1")
  expect_true(any(alg$cache$cached))
  rkip_check(solve(rkip_prob(c(0, 3)), alg, reltol = 1e-6, abstol = 1e-8), "reuse2")
})

test_that("RKIP keeps a step above the grid when not clamping (upstream indexes past it)", {
  sol <- solve(rkip_prob(c(0, 10)), RKIP(1e-4, 0.05, nb_of_cache_step = 37, clamp_higher_dt = FALSE))
  expect_identical(sol$retcode, "Success")
  expect_equal(sol$t[[length(sol$t)]], 10)
  expect_true(any(diff(sol$t) > 0.05 * (1 + 1e-12)))
})

test_that("RKIP grid, tableau and problem checks", {
  g <- .rkip_logrange(1e-3, 1, 100)
  expect_equal(g[c(1, 2, 50, 100)], c(0.001, 0.0010722672220103233, 0.030538555088334154, 1), tolerance = 1e-15)
  tab <- .rkip_verner6()
  expect_equal(rowSums(tab$A), tab$c, tolerance = 1e-14)
  expect_equal(sum(tab$alpha), 1, tolerance = 1e-14)
  expect_error(solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), RKIP()), "semilinear")
  L <- DiagonalOperator(c(-1, -2))
  prob <- SplitODEProblem(L, function(u, p, t) -u * u, c(1, 0.5), c(0, 1))
  sol <- solve(prob, RKIP(), reltol = 1e-10, abstol = 1e-12)
  ref <- solve(prob, Vern9(), reltol = 1e-12, abstol = 1e-14)
  expect_equal(sol$u[[length(sol$u)]], ref$u[[length(ref$u)]], tolerance = 1e-8)
})
