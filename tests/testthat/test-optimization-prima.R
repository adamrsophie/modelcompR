# UOBYQA, NEWUOA, BOBYQA and LINCOA of OptimizationPRIMA. The reference
# values come from PRIMA.jl (libprima 0.7.1, PRIMA_jll) through
# tools/prima-probe.jl; the ports agree with it in every bit on arm64 macOS,
# where the cosines and sines of R are those of libprima (the system libm).
# Elsewhere the last bits may differ, so the exact checks are guarded and a
# tolerance check of the solution runs everywhere.

if (!exists("ref_platform")) {
  ref_platform <- function() identical(Sys.info()[["sysname"]], "Darwin") && R.version$arch %in% c("aarch64", "arm64")
}

pp_rosen <- function(x, p) {
  a <- 1 - x[1]
  b <- x[2] - x[1] * x[1]
  a * a + 100 * b * b
}
pp_poly3 <- function(x, p) {
  a <- x[1] - 1
  b <- x[2] - 2
  cc <- x[3] + 1
  a * a + 2 * b * b + 3 * cc * cc + a * a * b * b + 0.5 * a * cc
}
pp_chained <- function(x, p) {
  s <- 0
  for (i in seq_len(length(x) - 1L)) {
    a <- 1 - x[i]
    b <- x[i + 1] - x[i] * x[i]
    s <- s + (a * a + 4 * (b * b))
  }
  s
}
pp_one <- function(x, p) (x[1] - 2) * (x[1] - 2) + 0.1 * (x[1] * x[1] * x[1] * x[1])
pp_nanrosen <- function(x, p) if (x[1] + x[2] > 1.5) NaN else pp_rosen(x, p)
pp_infrosen <- function(x, p) if (x[1] > 0.9) Inf else pp_rosen(x, p)
pp_status <- c("0" = "SMALL_TR_RADIUS", "1" = "FTARGET_ACHIEVED", "2" = "TRSUBP_FAILED", "3" = "MAXFUN_REACHED",
               "20" = "MAXTR_REACHED", "-1" = "NAN_INF_X", "-2" = "NAN_INF_F", "-3" = "NAN_INF_MODEL",
               "6" = "NO_SPACE_BETWEEN_BOUNDS", "7" = "DAMAGING_ROUNDING")
pp_sweepf <- function(x, w, cc, q) {
  s <- 0
  for (i in seq_along(x)) {
    t <- x[i] - cc[i]
    s <- s + w[i] * (t * t)
  }
  for (i in seq_len(length(x) - 1L)) {
    u <- x[i] * x[i + 1] - 0.5
    s <- s + q[i] * (u * u)
  }
  s
}

# solves fun from x0 with alg and the solve() arguments in ..., and checks
# the result against the reference ref of tools/prima-probe.jl
pp_check <- function(alg, ref, fun, x0, ..., lb = NULL, ub = NULL, cons = NULL, lcons = NULL, ucons = NULL,
                     xsol = NULL, tol = 1e-4) {
  pr <- if (is.null(cons)) OptimizationProblem(fun, x0, numeric(), lb = lb, ub = ub) else
    OptimizationProblem(fun, x0, numeric(), lb = lb, ub = ub, cons = cons, lcons = lcons, ucons = ucons)
  s <- solve(pr, alg, ...)
  if (!is.null(xsol)) expect_equal(unname(s$u), xsol, tolerance = tol)
  if (ref_platform()) {
    expect_identical(unname(s$u), ref$x)
    expect_identical(s$objective, ref$fx)
    expect_identical(s$original$nf, ref$nf)
    expect_identical(s$original$status, pp_status[[as.character(ref$status)]])
  }
  invisible(s)
}

# runs the LINCOA port on fun from x0 with the constraints and options in
# ... (the matrices of PRIMA.lincoa, so that no constraint is linearized
# numerically) and checks the result against ref
pp_check_lc <- function(ref, fun, x0, ..., xsol = NULL, tol = 1e-4) {
  r <- .lc_lincoa(function(x) fun(x, NULL), x0, ...)
  if (!is.null(xsol)) expect_equal(r$x, xsol, tolerance = tol)
  if (ref_platform()) {
    expect_identical(r$x, ref$x)
    expect_identical(r$f, ref$fx)
    expect_identical(r$nf, ref$nf)
    expect_identical(r$info, ref$status)
    expect_identical(r$cstrv, ref$cstrv)
  }
  invisible(r)
}

test_that("UOBYQA matches PRIMA.jl", {
  ref <- list(x = c(0x1.000000000a2ffp+0, 0x1.0000000014018p+0), fx = 0x1.d57b8844p-74, nf = 94L, status = 0L)
  pp_check(UOBYQA(), ref, pp_rosen, c(-1.2, 1), xsol = c(1, 1))
  ref <- list(x = c(0x1.ffffffffffdf3p-1, 0x1.000000000000cp+1, -0x1.ffffffffffff5p-1), fx = 0x1.0f332p-88, nf = 38L, status = 0L)
  pp_check(UOBYQA(), ref, pp_poly3, c(0, 0, 0), xsol = c(1, 2, -1))
  ref <- list(x = c(0x1.0000000000001p+0, 0x1p+1, -0x1p+0), fx = 0x1p-104, nf = 43L, status = 0L)
  pp_check(UOBYQA(), ref, pp_poly3, c(0, 0, 0), rhobeg = 0.5, rhoend = 1e-8, xsol = c(1, 2, -1))
  ref <- list(x = c(0x1.fffffff76d32fp-1, 0x1.fffffff68c66dp-1, 0x1.fffffff4c9755p-1, 0x1.ffffffe80b512p-1, 0x1.ffffffceae651p-1), fx = 0x1.52498a6e0b44cp-56, nf = 171L, status = 0L)
  pp_check(UOBYQA(), ref, pp_chained, rep(-1, 5), xsol = rep(1, 5))
  ref <- list(x = c(0x1.6c5e9749406cep+0), fx = 0x1.7c659563a9f4ep-1, nf = 18L, status = 0L)
  pp_check(UOBYQA(), ref, pp_one, 0, xsol = 1.42331834352847)
  ref <- list(x = c(0x1.6c5e974e9114fp+0), fx = 0x1.7c659563a9f4fp-1, nf = 57L, status = 0L)
  pp_check(UOBYQA(), ref, pp_one, 1000, xsol = 1.42331834352847)
  ref <- list(x = c(0x1.bf3436aa320efp-2, 0x1.5278e253e34d2p-3), fx = 0x1.8741e8d111b77p-2, nf = 30L, status = 3L)
  s <- pp_check(UOBYQA(), ref, pp_rosen, c(-1.2, 1), maxiters = 30)
  expect_identical(s$retcode, "MaxIters")
  expect_identical(s$original$nf, 30L)
  ref <- list(x = c(0x1.f15c011fa7964p-1, 0x1.e31dd9fa78094p-1), fx = 0x1.acc8121b8984bp-11, nf = 63L, status = 1L)
  s <- pp_check(UOBYQA(), ref, pp_rosen, c(-1.2, 1), ftarget = 1e-3)
  expect_identical(s$retcode, "Success")
  expect_lte(s$objective, 1e-3)
  ref <- list(x = c(0x1.000000351bbfap+0, 0x1.0000006a8f2aap+0), fx = 0x1.620754dd17adp-53, nf = 151L, status = 0L)
  pp_check(UOBYQA(), ref, pp_rosen, c(-1.2, 1), scale = c(2, 0.5), xsol = c(1, 1))
  ref <- list(x = c(-0x1.0f5b5568933a7p-1, 0x1.07d3fbe45f20ep-2), fx = 0x1.328d19961261dp+1, nf = 7L, status = 3L)
  expect_warning(s <- pp_check(UOBYQA(), ref, pp_rosen, c(-1.2, 1), maxiters = 3), "Invalid MAXFUN; it should be at least [(]N[+]1[)][(]N[+]2[)]/2 [+] 1; it is set to 7[.]")
  expect_identical(s$original$nf, 7L)
  ref <- list(x = c(0x1.000000000a2ffp+0, 0x1.0000000014018p+0), fx = 0x1.d57b8844p-74, nf = 94L, status = 0L)
  expect_warning(pp_check(UOBYQA(), ref, pp_rosen, c(-1.2, 1), rhoend = 0, xsol = c(1, 1)), "Invalid RHOEND")
  ref <- list(x = c(-0x1.0f5b5568933a7p-1, 0x1.07d3fbe45f20ep-2), fx = 0x1.328d19961261dp+1, nf = 8L, status = 0L)
  s <- pp_check(UOBYQA(), ref, pp_rosen, c(-1.2, 1), rhoend = 1)
  expect_identical(s$retcode, "Success")
  ref <- list(x = c(-0x1.3333333333333p+0, 0x1p+0), fx = 0x1.8333333333332p+4, nf = 1L, status = 1L)
  s <- pp_check(UOBYQA(), ref, pp_rosen, c(-1.2, 1), ftarget = 100)
  expect_identical(s$original$nf, 1L)
  expect_identical(unname(s$u), c(-1.2, 1))
  ref <- list(x = c(0x1.a575de78d2886p-1, 0x1.5a89501270da3p-1), fx = 0x1.00a82aaeeb73p-5, nf = 182L, status = 0L)
  s <- pp_check(UOBYQA(), ref, pp_nanrosen, c(-1.2, 1))
  expect_lte(sum(s$u), 1.5)
  ref <- list(x = c(0x1.cccb7be90926ap-1, 0x1.9e86d0e0b462fp-1), fx = 0x1.482b689ff9a1bp-7, nf = 164L, status = 0L)
  s <- pp_check(UOBYQA(), ref, pp_infrosen, c(-1.2, 1))
  expect_lte(s$u[1], 0.9)
})

test_that("NEWUOA matches PRIMA.jl", {
  ref <- list(x = c(0x1.00000364feea1p+0, 0x1.000006ccf0d8ep+0), fx = 0x1.7126eb710db1cp-45, nf = 168L, status = 0L)
  pp_check(NEWUOA(), ref, pp_rosen, c(-1.2, 1), xsol = c(1, 1))
  ref <- list(x = c(0x1.000000bd12b46p+0, 0x1.00000011234eep+1, -0x1.0000003f59017p+0), fx = 0x1.58e8cfe4710b9p-49, nf = 44L, status = 0L)
  pp_check(NEWUOA(), ref, pp_poly3, c(0, 0, 0), xsol = c(1, 2, -1))
  ref <- list(x = c(0x1.00000126df4a5p+0, 0x1.00000240ae5fep+0), fx = 0x1.964ead2d0aa5cp-48, nf = 223L, status = 0L)
  pp_check(NEWUOA(), ref, pp_rosen, c(-1.2, 1), npt = 4, xsol = c(1, 1))
  ref <- list(x = c(0x1.fffffffc1bf5ep-1, 0x1.fffffff83ba2dp-1), fx = 0x1.e5185f2f4dp-63, nf = 83L, status = 0L)
  pp_check(NEWUOA(), ref, pp_rosen, c(-1.2, 1), npt = 6, xsol = c(1, 1))
  ref <- list(x = c(0x1.00000000c1a75p+0, 0x1.00000000017fep+1, -0x1.0000000046cb8p+0), fx = 0x1.6509382c9ap-65, nf = 71L, status = 0L)
  pp_check(NEWUOA(), ref, pp_poly3, c(0, 0, 0), rhobeg = 0.5, rhoend = 1e-8, xsol = c(1, 2, -1))
  ref <- list(x = c(0x1.fffffb97ed77fp-1, 0x1.fffffaabb226ap-1, 0x1.fffffe1e6db9fp-1, 0x1.0000043bd939ep+0, 0x1.00000b35d5eaap+0), fx = 0x1.2b9abd253aacp-40, nf = 248L, status = 0L)
  pp_check(NEWUOA(), ref, pp_chained, rep(-1, 5), npt = 7, xsol = rep(1, 5))
  ref <- list(x = c(-0x1.b376042c0799bp-1, 0x1.59bb51bc8a358p-1), fx = 0x1.d3f28f496d33cp+1, nf = 30L, status = 3L)
  s <- pp_check(NEWUOA(), ref, pp_rosen, c(-1.2, 1), maxiters = 30)
  expect_identical(s$retcode, "MaxIters")
  expect_identical(s$original$nf, 30L)
  ref <- list(x = c(0x1.f1adaab0943a7p-1, 0x1.e353136615d94p-1), fx = 0x1.bfb0d52df2077p-11, nf = 135L, status = 1L)
  s <- pp_check(NEWUOA(), ref, pp_rosen, c(-1.2, 1), ftarget = 1e-3)
  expect_identical(s$retcode, "Success")
  expect_lte(s$objective, 1e-3)
  ref <- list(x = c(0x1.0000022d58a8p+0, 0x1.00000464406cfp+0), fx = 0x1.384712922151ep-46, nf = 208L, status = 0L)
  pp_check(NEWUOA(), ref, pp_rosen, c(-1.2, 1), scale = c(2, 0.5), xsol = c(1, 1))
  ref <- list(x = c(-0x1.7dd62a87191eap-1, 0x1.beed7e9c28e5p-4), fx = 0x1.708d9f89aa272p+4, nf = 5L, status = 3L)
  expect_warning(expect_warning(s <- pp_check(NEWUOA(), ref, pp_rosen, c(-1.2, 1), maxiters = 3), "it is set to 4[.]"), "Invalid MAXFUN; it should be at least N [+] 3; it is set to 5[.]")
  expect_identical(s$original$nf, 5L)
  ref <- list(x = c(0x1.00000364feea1p+0, 0x1.000006ccf0d8ep+0), fx = 0x1.7126eb710db1cp-45, nf = 168L, status = 0L)
  expect_warning(pp_check(NEWUOA(), ref, pp_rosen, c(-1.2, 1), rhoend = 0, xsol = c(1, 1)), "Invalid RHOEND")
  ref <- list(x = c(-0x1.3333333333333p+0, 0x1p+0), fx = 0x1.8333333333332p+4, nf = 6L, status = 0L)
  s <- pp_check(NEWUOA(), ref, pp_rosen, c(-1.2, 1), rhoend = 1)
  expect_identical(s$retcode, "Success")
  ref <- list(x = c(-0x1.3333333333333p+0, 0x1p+0), fx = 0x1.8333333333332p+4, nf = 1L, status = 1L)
  s <- pp_check(NEWUOA(), ref, pp_rosen, c(-1.2, 1), ftarget = 100)
  expect_identical(s$original$nf, 1L)
  expect_identical(unname(s$u), c(-1.2, 1))
  ref <- list(x = c(0x1.a584226064d33p-1, 0x1.5a7b06e0b1cd5p-1), fx = 0x1.00c625f8ccfe6p-5, nf = 164L, status = 0L)
  s <- pp_check(NEWUOA(), ref, pp_nanrosen, c(-1.2, 1))
  expect_lte(sum(s$u), 1.5)
  ref <- list(x = c(0x1.ccccc4a8088bbp-1, 0x1.9eceec561fdacp-1), fx = 0x1.47c78f6af51d8p-7, nf = 193L, status = 0L)
  s <- pp_check(NEWUOA(), ref, pp_infrosen, c(-1.2, 1))
  expect_lte(s$u[1], 0.9)
})

test_that("the fma emulation of the PRIMA ports is exact at the ends of the range", {
  # a factor above 2^995 and a product below 2^-960 (values of fma in Julia)
  expect_identical(.pw_fma(0x1.dfafcc93d821ep-170, -0x1.ab4c4bb59071dp+1023, 0x1.905498239f5c1p+854),
                   -0x1.9b1f49e71341fp+820)
  expect_identical(.pw_fma(-0x1.3b5fbcc864b96p-638, 0x1.4e524c97a2f83p+1023, 0x1.9bdc4b66206c7p+385),
                   0x1.2fa8a23ed28e4p+357)
  expect_identical(.pw_fma(0x1.f3ece012a188ap-383, 0x1.2ca8843c5f00cp-623, -0x1.25927fbe0d4a3p-1005),
                   -0x1.299b983ce28b1p-1021)
  expect_identical(.pw_fma(c(1, 0x1p1000), c(0x1p-1000, 3), c(-0x1p-1000, -0x1p1001)), c(0, 0x1p1000))
  # subnormal results (R reads subnormal hexadecimal literals as zero, so
  # these values are scaled from normal ones)
  expect_identical(.pw_fma(0x1.22bb9d0f41d2ap-519, 0x1.9e914faa7e29ap-526, 0x1.6a7a35b832348p-1000 * 2^-25),
                   0x1.6a7a532539208p-1000 * 2^-25)
  expect_identical(.pw_fma(-0x1.77887798e002fp-66, -0x1.307f4adf800e6p-957, -0x1.00ff659f1b8c8p-1000 * 2^-25),
                   0x1.7e6cef2666432p-1000 * 2^-23)
  # an infinite c with finite factors whose product overflows, NaN (not NA)
  # with a large factor, and the sign of a zero result
  expect_identical(.pw_fma(c(0x1p600, 5.000000000000016, 1), c(0x1p600, 1.049379054695898e308, 2), c(-Inf, -Inf, Inf)),
                   c(-Inf, -Inf, Inf))
  expect_true(is.nan(.pw_fma(0x1p1000, NaN, 1)))
  expect_identical(1 / .pw_fma(c(0, 2^-600, 1), c(-3, -2^-600, 2^-600), c(-0, 0, -2^-600)), c(-Inf, -Inf, Inf))
  # a product of 2^-1075, which rounds to zero, breaks the tie of a subnormal c
  expect_identical(.pw_fma(2^-500, c(2^-575, -2^-575), 3 * 2^-1074), c(4, 2) * 2^-1074)
  # a tie of a large product broken by the smallest c (fma in Julia)
  expect_identical(.pw_fma((1 + 2^-52) * 2^510, 1.5 * 2^500, c(-2^-1074, 0)),
                   c(0x1.8000000000001p+1010, 0x1.8000000000002p+1010))
})

test_that("UOBYQA and NEWUOA check their options and problems", {
  pr <- OptimizationProblem(pp_rosen, c(-1.2, 1), numeric())
  expect_error(solve(pr, UOBYQA(), npt = 6), "unsupported keyword argument `npt`")
  expect_error(solve(pr, NEWUOA(), npt = 2), "npt")
  expect_error(solve(pr, UOBYQA(), maxtime = 1), "not supported")
  expect_error(solve(pr, UOBYQA(), scale = c(1, -1)), "positive")
  prb <- OptimizationProblem(pp_rosen, c(-1.2, 1), numeric(), lb = c(-2, -2), ub = c(2, 2))
  expect_error(solve(prb, UOBYQA()), "does not support box constraints")
  expect_error(solve(prb, NEWUOA()), "does not support box constraints")
  prc <- OptimizationProblem(pp_rosen, c(-1.2, 1), numeric(), cons = function(u, p) u[1] + u[2], lcons = -Inf,
                             ucons = 1)
  expect_error(solve(prc, UOBYQA()), "does not support constraints")
  # iprint: the messages of libprima, MSG_FEVL (3) for every evaluation and
  # MSG_RHO (2) for every reduction of rho
  out <- capture.output(s <- solve(pr, UOBYQA(), ftarget = 100, iprint = "MSG_FEVL"))
  xl <- "The corresponding X is: -1.200000000000000E+000   1.000000000000000E+000"
  expect_identical(out, c("", "Initialization step with radius =  1.000000000000000E+000",
                          "Function number 1   F =  2.420000000000000E+001", xl, "", "",
                          "Return from UOBYQA because the target function value is achieved.",
                          "Number of function values = 1   Least value of F =  2.420000000000000E+001", xl))
  out <- capture.output(s <- solve(pr, NEWUOA(), iprint = 2))
  expect_identical(sum(startsWith(out, "New RHO = ")), 6L)
  expect_identical(capture.output(s <- solve(pr, UOBYQA(), iprint = 0)), character())
  expect_error(solve(pr, UOBYQA(), iprint = 4), "iprint must be")
  # a callback that returns TRUE stops the solver at the best point so far
  s <- solve(pr, UOBYQA(), callback = function(state, f) state$iter >= 10)
  expect_identical(s$retcode, "Terminated")
  expect_lte(s$objective, pp_rosen(c(-1.2, 1)))
})

test_that("BOBYQA and LINCOA check their options and problems", {
  prc <- OptimizationProblem(pp_rosen, c(-1.2, 1), numeric(), cons = function(u, p) u[1] + u[2], lcons = -Inf,
                             ucons = 1)
  expect_error(solve(prc, BOBYQA()), "does not support constraints")
  prn <- OptimizationProblem(pp_rosen, c(-1.2, 1), numeric(), cons = function(u, p) u[1] * u[1] + u[2],
                             lcons = -Inf, ucons = 1)
  expect_error(solve(prn, LINCOA()), "supports only linear constraints")
  # no space between the bounds of a variable
  pr0 <- OptimizationProblem(pp_rosen, c(-1.2, 1), numeric(), lb = c(-2, 1), ub = c(2, 1))
  expect_warning(s <- solve(pr0, BOBYQA()), "no space between the lower and upper bounds of variable 2")
  expect_identical(s$retcode, "InitialFailure")
  # x0 is moved rhobeg away from the inactive bounds, and rhobeg is revised
  prb <- OptimizationProblem(pp_rosen, c(-1.2, 1), numeric(), lb = c(-2, -2), ub = c(2, 2))
  expect_warning(solve(prb, BOBYQA(), maxiters = 10), "X0 is revised")
  expect_warning(expect_warning(solve(prb, BOBYQA(), maxiters = 10, rhobeg = 3), "Invalid RHOBEG; BOBYQA requires"),
                 "X0 is revised")
  # LINCOA relaxes the constraints that x0 violates and ignores zero rows
  expect_warning(.lc_lincoa(function(x) pp_rosen(x, NULL), c(1, 1), Aineq = matrix(c(1, 1), 1), bineq = 1,
                            maxfun = 6L), "The starting point is infeasible")
  expect_warning(.lc_lincoa(function(x) pp_rosen(x, NULL), c(-1.2, 1), Aineq = rbind(c(0, 0), c(1, 1)),
                            bineq = c(1, 1), maxfun = 6L), "zero gradients")
  # the messages of LINCOA carry the constraint violation and values
  out <- capture.output(r <- .lc_lincoa(function(x) pp_rosen(x, NULL), c(-1.2, 1), Aineq = matrix(c(1, 1), 1),
                                        bineq = 1, ftarget = 100, iprint = 3L))
  xl <- "The corresponding X is: -1.200000000000000E+000   1.000000000000000E+000"
  cl <- "The constraint value is: -1.200000000000000E+000"
  expect_identical(out, c("", "Initialization step with radius =  1.000000000000000E+000",
                          paste0("Function number 1   F =  2.420000000000000E+001   Constraint violation = ",
                                 " 0.000000000000000E+000"), xl, cl, "", "",
                          "Return from LINCOA because the target function value is achieved.",
                          paste0("Number of function values = 1   Least value of F =  2.420000000000000E+001   ",
                                 "Constraint violation =  0.000000000000000E+000"), xl, cl))
  # bounds and linear constraints through solve()
  s <- solve(OptimizationProblem(pp_rosen, c(-1.2, 1), numeric(), lb = c(-2, -2), ub = c(0.5, 2)), LINCOA())
  expect_equal(unname(s$u), c(0.5, 0.25), tolerance = 1e-4)
  s <- solve(OptimizationProblem(pp_rosen, c(0, 0), numeric(), cons = function(u, p) c(u[1] - u[2], u[1] + u[2]),
                                 lcons = c(0, -Inf), ucons = c(0, 1)), LINCOA())
  expect_equal(s$u[1], s$u[2], tolerance = 1e-8)
  expect_lte(s$original$cstrv, 1e-10)
})
