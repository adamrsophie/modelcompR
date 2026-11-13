# UOBYQA, NEWUOA, BOBYQA and LINCOA of OptimizationPRIMA. The reference
# values come from PRIMA.jl (libprima 0.7.1, PRIMA_jll) through
# tools/prima-probe.jl; the ports agree with it in every bit on arm64 macOS,
# where the cosines and sines of R are those of libprima (the system libm).
# Elsewhere the last bits may differ, so the exact checks are guarded and a
# tolerance check of the solution runs everywhere.

# long reference comparisons: run on CI (NOT_CRAN = true), not on CRAN
skip_on_cran()

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

test_that("BOBYQA matches PRIMA.jl", {
  ref <- list(x = c(0x1.0000001fe4fbcp+0, 0x1.000000232960bp+0), fx = 0x1.44197e962a1b1p-48, nf = 164L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_rosen, c(-1.2, 1), lb = c(-2, -2), ub = c(2, 2), xsol = c(1, 1)))
  ref <- list(x = c(0x1p-1, 0x1.00000000002aap-2), fx = 0x1p-2, nf = 91L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_rosen, c(-1.2, 1), lb = c(-2, -2), ub = c(0.5, 2), xsol = c(0.5, 0.25)))
  ref <- list(x = c(0x1.00000005ef306p+0, 0x1.0000000b63ap+0), fx = 0x1.d1ab5a48e32p-59, nf = 138L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_rosen, c(-1.2, 1), xsol = c(1, 1)))
  ref <- list(x = c(0x1.999999999999ap-1, 0x1.47ae147ae14f7p-1), fx = 0x1.47ae147ae1478p-5, nf = 83L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_rosen, c(-0.5, 1.5), lb = c(-2, -2), ub = c(0.8, 4), xsol = c(0.8, 0.64)))
  ref <- list(x = c(0x1.000027f1dc26ep+0, 0x1.000050890fd39p+0), fx = 0x1.9952543dd568ap-38, nf = 57L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_rosen, c(5, -5), lb = c(-2, -2), ub = c(2, 2), xsol = c(1, 1)))
  ref <- list(x = c(0x1p-1, 0x1.fffffffff4e2bp+0, -0x1.eaaaaaaaa321dp-1), fx = 0x1.f555555555555p-3, nf = 42L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_poly3, c(0, 0, 0), lb = c(-1, -1, -3), ub = c(0.5, 3, 3)))
  ref <- list(x = c(0x1.000000001803ep+0, 0x1.ffffffffee46cp+0, -0x1.0000000026858p+0), fx = 0x1.44a25623cp-68, nf = 54L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_poly3, c(0, 0, 0), lb = rep(-3, 3), ub = rep(3, 3), rhobeg = 0.5, rhoend = 1e-8, xsol = c(1, 2, -1)))
  ref <- list(x = c(0x1.ffffee3caa4afp-1, 0x1.ffffdc60de1d3p-1), fx = 0x1.3c703c4326ddep-42, nf = 161L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_rosen, c(-1.2, 1), lb = c(-2, -2), ub = c(2, 2), rhobeg = 3, xsol = c(1, 1)))
  ref <- list(x = c(0x1.999999999999ap-1, 0x1.999999999999ap-1, 0x1.999999999999ap-1, 0x1.6c8b46a5617c2p-1, 0x1.038e309fffd6fp-1), fx = 0x1.b6c3760bf61bp-2, nf = 132L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_chained, rep(-1, 5), lb = rep(-2, 5), ub = rep(0.8, 5)))
  ref <- list(x = c(0x1.fffffb43137cap-1, 0x1.fffff52f1f11dp-1, 0x1.ffffe3157be6ep-1, 0x1.ffffc4e84f0d5p-1, 0x1.ffff70205d283p-1), fx = 0x1.ca4e4897ad6c2p-38, nf = 280L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_chained, rep(-1, 5), lb = rep(-2, 5), ub = rep(2, 5), npt = 7, xsol = rep(1, 5)))
  ref <- list(x = c(0x1.0000000697861p+0, 0x1.0000000e7d59bp+0, 0x1.fffffffd53b7ep-1, 0x1.ffffffd339c5bp-1, 0x1.ffffff8b81933p-1), fx = 0x1.a4006096a4283p-52, nf = 235L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_chained, rep(-1, 5), lb = rep(-2, 5), ub = rep(2, 5), npt = 21, xsol = rep(1, 5)))
  ref <- list(x = c(0x1.6c5e974e966fdp+0), fx = 0x1.7c659563a9f4ep-1, nf = 15L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_one, 0, lb = -5, ub = 5, xsol = 1.42331834352847))
  ref <- list(x = c(0x1.6c5e974ac0719p+0), fx = 0x1.7c659563a9f4ep-1, nf = 60L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_one, 1000, lb = -2000, ub = 2000, xsol = 1.42331834352847))
  ref <- list(x = c(-0x1.2671412c0ce14p-1, 0x1.46d94cc24d5bap-2), fx = 0x1.3f416f65fb7f7p+1, nf = 30L, status = 3L, cstrv = 0x0p+0)
  suppressWarnings(s <- pp_check(BOBYQA(), ref, pp_rosen, c(-1.2, 1), lb = c(-2, -2), ub = c(2, 2), maxiters = 30))
  expect_identical(s$retcode, "MaxIters")
  expect_identical(s$original$nf, 30L)
  ref <- list(x = c(0x1.f6902fa8f37d4p-1, 0x1.ed1dffae46a28p-1), fx = 0x1.71a2d33539bffp-12, nf = 127L, status = 1L, cstrv = 0x0p+0)
  suppressWarnings(s <- pp_check(BOBYQA(), ref, pp_rosen, c(-1.2, 1), lb = c(-2, -2), ub = c(2, 2), ftarget = 1e-3))
  expect_identical(s$retcode, "Success")
  expect_lte(s$objective, 1e-3)
  ref <- list(x = c(0x1.ffffff9b54a2ep-1, 0x1.ffffffe103a9fp-1), fx = 0x1.637cbee77f053p-45, nf = 158L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_rosen, c(-1.2, 1), lb = c(-2, -2), ub = c(2, 2), scale = c(2, 0.5), xsol = c(1, 1)))
  ref <- list(x = c(-0x1p+0, 0x1p+0), fx = 0x1p+2, nf = 5L, status = 3L, cstrv = 0x0p+0)
  suppressWarnings(s <- pp_check(BOBYQA(), ref, pp_rosen, c(-1.2, 1), lb = c(-2, -2), ub = c(2, 2), maxiters = 3))
  expect_identical(s$original$nf, 5L)
  ref <- list(x = c(0x1.0000001fe4fbcp+0, 0x1.000000232960bp+0), fx = 0x1.44197e962a1b1p-48, nf = 164L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(pp_check(BOBYQA(), ref, pp_rosen, c(-1.2, 1), lb = c(-2, -2), ub = c(2, 2), rhoend = 0, xsol = c(1, 1)))
  ref <- list(x = c(0x1.a564c3c8f8d0bp-1, 0x1.5a9b15a8ad9a1p-1), fx = 0x1.00b1570bfcbccp-5, nf = 150L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(s <- pp_check(BOBYQA(), ref, pp_nanrosen, c(-1.2, 1), lb = c(-2, -2), ub = c(2, 2)))
  expect_lte(sum(s$u), 1.5)
  ref <- list(x = c(0x1.cccb87b8e3668p-1, 0x1.9e7ed73cfb3bbp-1), fx = 0x1.485314c17ee22p-7, nf = 171L, status = 0L, cstrv = 0x0p+0)
  suppressWarnings(s <- pp_check(BOBYQA(), ref, pp_infrosen, c(-1.2, 1), lb = c(-2, -2), ub = c(2, 2)))
  expect_lte(s$u[1], 0.9)
})

test_that("LINCOA matches PRIMA.jl", {
  suppressWarnings({
    ref <- list(x = c(0x1.fffff5fdb9a49p-1, 0x1.ffffebb765f31p-1), fx = 0x1.acfa4de2bab4cp-44, nf = 152L, status = 0L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_rosen, c(-1.2, 1), xsol = c(1, 1))
    ref <- list(x = c(0x1.3cd2c785ab46bp-1, 0x1.865a70f4a972ap-2), fx = 0x1.2a340323fee4p-3, nf = 91L, status = 0L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_rosen, c(-1.2, 1), Aineq = matrix(c(1, 1), 1), bineq = 1)
    ref <- list(x = c(0x1.4e82775d88661p-7, 0x1.4e82775d88661p-7), fx = 0x1.fad3c3594d8b6p-1, nf = 26L, status = 0L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_rosen, c(0, 0), Aeq = matrix(c(1, -1), 1), beq = 0, xl = c(-2, -2), xu = c(0.5, 2))
    ref <- list(x = c(0x1p-1, 0x1.0000000000604p-2), fx = 0x1p-2, nf = 89L, status = 0L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_rosen, c(-1.2, 1), xl = c(-2, -2), xu = c(0.5, 2), xsol = c(0.5, 0.25))
    ref <- list(x = c(0x1.723e7af99a1f6p-1, 0x1.dac0ffebe4b3dp+0, -0x1.13e03d68b1c37p+0), fx = 0x1.3210e7619ebabp-3, nf = 40L, status = 0L, cstrv = 0x1p-52)
    pp_check_lc(ref, pp_poly3, c(0, 0, 0), Aineq = rbind(c(1, 1, 1), c(-1, 0, 0)), bineq = c(1.5, 0))
    ref <- list(x = c(0x1.0f7f516c749fcp+0, 0x1.94989894f1ed7p-1, -0x1.71b632fe5a4aap+0), fx = 0x1.c1de1c783c1f2p+1, nf = 46L, status = 0L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_poly3, c(1, 1, 1), Aineq = matrix(1, 1, 3), bineq = 1.5)
    ref <- list(x = c(0x1.3cd2c785ab46bp-1, 0x1.865a70f4a972ap-2), fx = 0x1.2a340323fee4p-3, nf = 91L, status = 0L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_rosen, c(-1.2, 1), Aineq = rbind(c(0, 0), c(1, 1)), bineq = c(1, 1))
    ref <- list(x = c(0x1.c6b99a2ac40c1p-1, 0x1.9e9839617a90ep-1, 0x1.5beea8960c13bp-1, 0x1.dfa0bd3b51642p-2, 0x1.3bbc95003278ap-3), fx = 0x1.d226f25c05aacp-2, nf = 111L, status = 0L, cstrv = 0x1p-49)
    pp_check_lc(ref, pp_chained, rep(-1, 5), Aineq = matrix(1, 1, 5), bineq = 3)
    ref <- list(x = c(0x1.c6b99942dff79p-1, 0x1.9e98397190879p-1, 0x1.5beeaa230189ep-1, 0x1.dfa0c12b4c44cp-2, 0x1.3bbc8a4b9f53cp-3), fx = 0x1.d226f25c06025p-2, nf = 163L, status = 0L, cstrv = 0x1p-50)
    pp_check_lc(ref, pp_chained, rep(-1, 5), Aineq = matrix(1, 1, 5), bineq = 3, npt = 7L)
    ref <- list(x = c(0x1.c6b99a886ad51p-1, 0x1.9e9839b23392fp-1, 0x1.5beea8be9dd7dp-1, 0x1.dfa0bd12f6635p-2, 0x1.3bbc91f52239cp-3), fx = 0x1.d226f25c05aa2p-2, nf = 113L, status = 0L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_chained, rep(-1, 5), Aineq = matrix(1, 1, 5), bineq = 3, npt = 21L)
    ref <- list(x = c(0x1.3cd2c7c0b0c66p-1, 0x1.865a717194388p-2), fx = 0x1.2a34027e6deabp-3, nf = 96L, status = 0L, cstrv = 0x1.e5eb8a8p-27)
    pp_check_lc(ref, pp_rosen, c(-1.2, 1), Aineq = matrix(c(1, 1), 1), bineq = 1, rhobeg = 0.5, rhoend = 1e-8)
    ref <- list(x = c(-0x1.b2e470c62af76p-2, 0x1.141e54021686cp-3), fx = 0x1.1e5d0e501cb1fp+1, nf = 30L, status = 3L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_rosen, c(-1.2, 1), Aineq = matrix(c(1, 1), 1), bineq = 1, maxfun = 30L)
    ref <- list(x = c(0x1.2428bb65ad972p-1, 0x1.4a0caa71a9b88p-2), fx = 0x1.7bcdef1c35071p-3, nf = 68L, status = 1L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_rosen, c(-1.2, 1), Aineq = matrix(c(1, 1), 1), bineq = 1, ftarget = 0.2)
    ref <- list(x = c(0x1.3cd2c785ab46bp-1, 0x1.865a70f4a972ap-2), fx = 0x1.2a340323fee4p-3, nf = 91L, status = 0L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_nanrosen, c(-1.2, 1), Aineq = matrix(c(1, 1), 1), bineq = 1)
    ref <- list(x = c(0x1.3cd2c785ab46bp-1, 0x1.865a70f4a972ap-2), fx = 0x1.2a340323fee4p-3, nf = 91L, status = 0L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_infrosen, c(-1.2, 1), Aineq = matrix(c(1, 1), 1), bineq = 1)
    ref <- list(x = c(0x1p+0), fx = 0x1.199999999999ap+0, nf = 16L, status = 0L, cstrv = 0x0p+0)
    pp_check_lc(ref, pp_one, 0, Aineq = matrix(1, 1, 1), bineq = 1, xsol = 1)
  })
})

test_that("LINCOA solves linearly constrained problems through solve()", {
  ref <- list(x = c(0x1.3cd2c785ab741p-1, 0x1.865a70f4a917ep-2), fx = 0x1.2a340323fee4p-3, nf = 126L, status = 0L, cstrv = 0x0p+0)
  r <- .lc_lincoa(function(x) pp_rosen(c(2, 0.5) * x, NULL), c(-1.2, 1) / c(2, 0.5), Aineq = matrix(c(2, 0.5), 1), bineq = 1)
  if (ref_platform()) {
    expect_identical(r$x * c(2, 0.5), ref$x)
    expect_identical(r$nf, ref$nf)
  }
  prc <- OptimizationProblem(pp_rosen, c(-1.2, 1), numeric(), cons = function(u, p) u[1] + u[2], lcons = -Inf, ucons = 1)
  s <- solve(prc, LINCOA())
  expect_identical(s$retcode, "Success")
  expect_equal(unname(s$u), ref$x, tolerance = 1e-4)
  s <- solve(prc, LINCOA(), scale = c(2, 0.5))
  expect_equal(unname(s$u), ref$x, tolerance = 1e-4)
  expect_lte(s$original$cstrv, 1e-12)
})

test_that("BOBYQA matches PRIMA.jl on a seeded sweep", {
  sweep <- list(
    list(w = c(0x1.7df9b19dace14p-1), c = c(-0x1.1d31659c63073p-1), q = c(), x0 = c(-0x1.7e00e2cee96f8p+0), xl = c(-0x1.58188194fa393p+1), xu = c(-0x1.f153022ef7281p+0), npt = 3L, rhobeg = 0x1.316334c4c876fp-3, x = c(-0x1.f153022ef7281p+0), fx = 0x1.6eb442b4b9b48p+0, nf = 11L, status = 0L),
    list(w = c(0x1.caf827d1da25cp-1, 0x1.832d051100bdcp-1), c = c(0x1.7459d3960fd49p-4, 0x1.9f0f9d7896abcp-1), q = c(0x1.e5b6e7b732ca2p-1), x0 = c(0x1.e3441c5a778e4p-2, 0x1.75cf09608d443p-3), xl = c(-0x1.48193622a2fe2p-2, 0x1.17b9e5dc929p-6), xu = c(0x1.6111c6fe9a208p+0, 0x1.402dd71769dc7p+0), npt = 5L, rhobeg = 0x1p+0, x = c(0x1.38c663b253b77p-2, 0x1.cb636da90d761p-1), fx = 0x1.8694b9e6f86e6p-4, nf = 24L, status = 0L),
    list(w = c(0x1.cd55fdb44f07p-1), c = c(0x1.70f867e9145b6p-2), q = c(), x0 = c(-0x1.efae137110414p+0), xl = c(-0x1.1059880fab70ep+1), xu = c(-0x1.7fbcfd501258p-4), npt = 3L, rhobeg = 0x1p+0, x = c(-0x1.7fbcfd501258p-4), fx = 0x1.7c5e90f7a6557p-3, nf = 14L, status = 0L),
    list(w = c(0x1.97227f3861aa6p+0, 0x1.19ff2ae8a39e4p+0, 0x1.46d6b641c7b54p-1, 0x1.61a288733f6bep+0), c = c(0x1.47df007a718fep+0, -0x1.6b52a0e282a85p-1, 0x1.0866ff41b2feep-2, 0x1.083421674c253p+1), q = c(0x1.baaf85ddc4eep-4, 0x1.0a26370f7b2a8p-3, 0x1.cc325c7ef2d04p-2), x0 = c(0x1.52dd442045726p-2, -0x1.c726beed74504p+0, 0x1.f2cc531939b08p-4, -0x1.527b2480f17a6p+0), xl = c(0x1.2a88a10c0aa4p-6, -0x1.f9b119abd5974p+0, -0x1.c2db7b596d095p-1, -0x1.f9a95cccf0edp+0), xu = c(0x1.70ec65c0c06dep-1, 0x1.d1438eae67536p-1, 0x1.8c6fd6d7bc052p+0, -0x1.e9976426881ap-3), npt = 15L, rhobeg = 0x1p+0, x = c(0x1.70ec65c0c06dep-1, -0x1.4558057c2caa3p-1, 0x1.8f21b7691505dp-4, -0x1.e9976426881ap-3), fx = 0x1.03994b21ae575p+3, nf = 62L, status = 0L),
    list(w = c(0x1.92ad1417cefeap+0, 0x1.23c72cfda491ap+1, 0x1.bc25d54d40978p+0, 0x1.df50cb9f681cep+0), c = c(-0x1.127e08edf2012p+1, 0x1.1923ffb0679f3p+1, 0x1.7bd43509fddbbp-2, -0x1.8610f5048fd32p-3), q = c(0x1.afefb7722c7fap-1, 0x1.7e47040069058p-2, 0x1.0adadfdc839d8p-1), x0 = c(0x1.878d9f6a209bcp+1, -0x1.5249409b8d74cp-1, -0x1.18499b8dc843p-2, 0x1.4fe57c7e7ca4dp+2), xl = c(-0x1.e89425206b12ap-1, -0x1.fd7bd284cf55fp-1, -0x1.0c329ac222facp-2, -0x1.82a5145f3d62bp+0), xu = c(0x1.94fee1478b786p-1, -0x1.341c8d38b39ep-2, 0x1.ca76ff708e8d7p+0, 0x1.4f47689821d96p-1), npt = 6L, rhobeg = 0x1.1c57a3205e059p-3, x = c(-0x1.e89425206b12ap-1, -0x1.341c8d38b39ep-2, 0x1.3c7652f522399p-2, -0x1.2626db86d955fp-3), fx = 0x1.0c77a867c86ebp+4, nf = 53L, status = 0L),
    list(w = c(0x1.c119a61d574bcp+0, 0x1.01bb467a8e374p-1, 0x1.661d1c8e90e8cp-1, 0x1.3c59ccc71f12ep+0, 0x1.0f093bdc80f9cp+0), c = c(-0x1.247fe43514f51p-1, -0x1.1a0eb2b24aecep+0, -0x1.00d6a9d5de582p+0, 0x1.5dda4067413c7p+1, 0x1.6604cf87f3f75p+0), q = c(0x1.49490dc6415ap-5, 0x1.a21c1f27d8feep-1, 0x1.0f6628a48b49p-1, 0x1.fde9140bca04cp-2), x0 = c(-0x1.183137dc3e158p-2, 0x1.aa75db50c409ep-2, -0x1.3d450833d6bap-3, 0x1.ffbb1907676p-7, -0x1.d3a2041676e69p+0), xl = c(-0x1.e4da35258905p-2, -0x1.04339834ada1p-3, -0x1.8ed816170086dp-1, -0x1.6485af272f8dp+0, -0x1.7fb04e6bed4ccp+1), xu = c(0x1.325ca5b5e2d4ap+1, 0x1.0483819983ba2p+1, 0x1.d2e373343733ep+0, 0x1.bbb1784051dc7p+0, -0x1.4a9a69a079ee9p+0), npt = 7L, rhobeg = 0x1p+0, x = c(-0x1.e4da35258905p-2, -0x1.04339834ada1p-3, -0x1.a67b65654c38p-3, 0x1.6fd8402a8e4bap+0, -0x1.4a9a69a079ee9p+0), fx = 0x1.bef2daa4a7fb2p+3, nf = 58L, status = 0L),
    list(w = c(0x1.12384bab3e097p+1, 0x1.5c652e253c398p-1, 0x1.162c3f17702b8p+0, 0x1.b1369c8147504p-1), c = c(-0x1.23f821b69c28cp-1, 0x1.c67fb703f0d0ap-2, -0x1.3afcea25ab31fp-2, 0x1.f7b24f0fd796ap-1), q = c(0x1.ddb591f7a31dap-1, 0x1.4e121afae1d1cp-2, 0x1.31854c6dd19fp-3), x0 = c(-0x1.aff39d0d4b3d2p-3, -0x1.9e0e7d6f3dde1p+0, 0x1.b02fb710244fbp+0, -0x1.e35572e6d624p-2), xl = c(-0x1.924be7de76eb4p-2, -0x1.dfd1ec42e9322p+0, 0x1.61200071b84d6p+0, -0x1.0b94e2a867d06p-1), xu = c(0x1.4373c84651935p+0, -0x1.2dc4f69309356p+0, 0x1.e00699b0e275p+1, -0x1.1f7286a73f58ap-3), npt = 9L, rhobeg = 0x1p+0, x = c(-0x1.924be7de76eb4p-2, -0x1.2dc4f69309356p+0, 0x1.61200071b84d6p+0, -0x1.1f7286a73f58ap-3), fx = 0x1.e45f9de886cabp+2, nf = 25L, status = 0L),
    list(w = c(0x1.09373a67b09d4p+0, 0x1.2e2b56c69b678p-1, 0x1.0bc0d44386796p+0), c = c(0x1.19abe339f7bd6p-1, 0x1.55ee0385cde79p+0, -0x1.e2bd66d6fb4e6p-2), q = c(0x1.fadb1b01d716ap-1, 0x1.6bad4fa14492ap-1), x0 = c(0x1.1d78f61b4caecp-1, -0x1.b256e3184b184p+0, -0x1.5f2519c473f7ap-1), xl = c(-0x1.117861ba7627ap-1, -0x1.de659929c27a8p+0, -0x1.dd9f47e0101e2p+0), xu = c(0x1.1f48404b0aa42p+0, -0x1.4a270322639cbp+0, 0x1.a3273074916d4p-2), npt = 7L, rhobeg = 0x1p+0, x = c(-0x1.a1ce6ee3e92c8p-6, -0x1.4a270322639cbp+0, -0x1.b543befadbb06p-2), fx = 0x1.285fc8111c6b1p+2, nf = 32L, status = 0L),
    list(w = c(0x1.02615ccd71c6fp+1, 0x1.2c55f05c0a3aap+1, 0x1.2fe18f6da258dp+1, 0x1.3e4eff27ea934p+1, 0x1.3093994d35fc7p+1), c = c(-0x1.12479ff13b8e3p+0, -0x1.4bb6bc8cfdff1p-2, 0x1.c1522daa78b98p-2, -0x1.e44cb7ad975bap-2, -0x1.8756440f074d8p+0), q = c(0x1.477a3e4c09a1p-4, 0x1.ccb09a81201aap-1, 0x1.b954b0877d38p-1, 0x1.fa48d12b44d14p-2), x0 = c(-0x1.5335656400d16p-1, -0x1.9fab19ab2f9fp-2, -0x1.26d91ec8318f7p+1, -0x1.dc72abdc0729dp-2, -0x1.b937a9b9afec9p-1), xl = c(-0x1.a971a02ab1438p-1, -0x1.65805e223ec1p-1, -0x1.2e73fe5680a76p+1, -0x1.bc066a4a2aa29p-1, -0x1.6139248cb2c5cp+0), xu = c(-0x1.ea1fee7261838p-2, -0x1.89a7124f4b301p-2, -0x1.00e7958d9d3a6p+1, -0x1.669f6fe4956a8p-2, 0x1.c9dbaabaf5735p+0), npt = 7L, rhobeg = 0x1.011487f75b74cp-4, x = c(-0x1.a971a02ab1438p-1, -0x1.89a7124f4b301p-2, -0x1.00e7958d9d3a6p+1, -0x1.669f6fe4956a8p-2, -0x1.6139248cb2c5cp+0), fx = 0x1.d0b6c5bd1fcbbp+3, nf = 29L, status = 0L),
    list(w = c(0x1.ddf3524f63898p+0, 0x1.9386709837588p+0, 0x1.13f4cf0021f6dp+1, 0x1.5f8751a6d7a28p+0, 0x1.3b9a81b6adf9ep+1), c = c(-0x1.cd6a29ad6b7c8p-1, -0x1.c4077ddab9852p-3, 0x1.25cf5f3bd8264p-5, 0x1.26401a7f46415p-2, -0x1.b1ac968e3a2b6p+0), q = c(0x1.a3e794926c1p-7, 0x1.d41226a2ca38p-6, 0x1.49ac9b0933fdap-1, 0x1.b0c94742904f6p-1), x0 = c(-0x1.a9b6d9a7e9ea8p-2, 0x1.20161bd123776p+0, -0x1.3109c6b72094p-1, 0x1.96dafb43256dap+0, -0x1.ec99ecd00f6a8p-2), xl = c(-0x1.5e5b2b3c80fb2p-1, 0x1.807201521e328p-3, -0x1.7648490c5640cp-1, -0x1.3893e99a3a658p-4, -0x1.733c4760f17cp-1), xu = c(0x1.02412244bf5ccp-1, 0x1.55f7d0862de84p+0, 0x1.c61fd78ac544fp+0, 0x1.28de7ae19d083p+1, 0x1.c6617c3e2d05dp+0), npt = 11L, rhobeg = 0x1p+0, x = c(-0x1.5e5b2b3c80fb2p-1, 0x1.807201521e328p-3, 0x1.7498d438b5eaep-5, 0x1.cee8ca153e33ap-5, -0x1.733c4760f17cp-1), fx = 0x1.94340d0206f6dp+1, nf = 62L, status = 0L),
    list(w = c(0x1.aefbbc2c1e74ep+0, 0x1.2087f5abc4c04p+1, 0x1.9aa4b50647944p+0, 0x1.fac13997e58fp-1, 0x1.7c4477ea5c65cp-1), c = c(-0x1.bb533d0656f04p-1, -0x1.020bcae41ec06p-4, 0x1.1a0929ac36879p-2, -0x1.50bd5d3e48f57p+0, -0x1.e7030db760e5cp+0), q = c(0x1.0fb3681b8ae06p-1, 0x1.d0af21fb891f6p-1, 0x1.a77de463774f6p-1, 0x1.35dbe6ecd0ccap-1), x0 = c(0x1.87bc7ceda51a2p+0, -0x1.9cc1908b1f11ap+0, 0x1.3a5597bbf506p+0, -0x1.4446319622b94p-1, -0x1.beba4952795c3p+0), xl = c(0x1.3b908b65664c8p-3, -0x1.201a0d6521291p+1, -0x1.1942ff6e96c78p-1, -0x1.85581c6b2a231p+0, -0x1.80ae085a7cad1p+1), xu = c(0x1.25b6da5603846p+1, -0x1.cdb1a00f316p-3, 0x1.77e576c8de4b2p+0, -0x1.5cadfa3b7e3c8p-2, -0x1.35bcfca39f344p-1), npt = 21L, rhobeg = 0x1p+0, x = c(0x1.3b908b65664c8p-3, -0x1.cdb1a00f316p-3, 0x1.2fbbf6282d16p-6, -0x1.74644d2b583adp-1, -0x1.89537cffefe5p+0), fx = 0x1.984996c1c352fp+1, nf = 93L, status = 0L),
    list(w = c(0x1.13f410f7f8965p+1), c = c(0x1.257cfea38656bp-1), q = c(), x0 = c(-0x1.2dae5480da8a8p-3), xl = c(-0x1.b198c7aedfcc2p+0), xu = c(0x1.7eef741636a1bp+0), npt = 3L, rhobeg = 0x1.469ce4b53c2bfp-1, x = c(0x1.257cfea38656bp-1), fx = 0x0p+0, nf = 13L, status = 0L),
    list(w = c(0x1.088b778daf25p-1), c = c(0x1.620f206c730bdp+0), q = c(), x0 = c(0x1.9c5f3999d9f68p-2), xl = c(-0x1.24bb18ddc114cp-2), xu = c(0x1.f6efaf59797afp-2), npt = 3L, rhobeg = 0x1p+0, x = c(0x1.f6efaf59797afp-2), fx = 0x1.a4e0c91c2c4cap-2, nf = 10L, status = 0L),
    list(w = c(0x1.70a6839268ac4p-1, 0x1.c69d45cec04c4p+0, 0x1.ac3fd05554962p+0, 0x1.1bc3fb63f05ecp+1), c = c(0x1.cf88a7190d1f3p-4, -0x1.e9c43283528dbp+0, 0x1.1f246c9ad3b8ep-4, 0x1.3518db8d73e56p+0), q = c(0x1.41ebc34d21bf8p-3, 0x1.f6d0aad9d7732p-1, 0x1.d03c5a4f9ca64p-1), x0 = c(0x1.45edd3f21563ap-1, -0x1.e04004c741886p+0, 0x1.16af1563e1b8ep+1, 0x1.5bcc2462bc5e4p-1), xl = c(0x1.2b0f11f349f28p-2, -0x1.205c6a7792d08p+1, 0x1.8ff69a278be6p-5, 0x1.efa80be7e392p-3), xu = c(0x1.03a4fdde88f21p+0, -0x1.3205a24923718p+0, 0x1.2b084babfb9a4p+1, 0x1.68119a46e0851p+0), npt = 15L, rhobeg = 0x1p+0, x = c(0x1.2b0f11f349f28p-2, -0x1.dec183525d064p+0, 0x1.8ff69a278be6p-5, 0x1.37597b8e2a83dp+0), fx = 0x1.6ff9530b46f66p-1, nf = 63L, status = 0L),
    list(w = c(0x1.f9c7f2e5199e4p-1, 0x1.6467e4c0f4b4cp+0, 0x1.9961dd19c58bcp-1), c = c(0x1.baed561e25b03p-2, -0x1.337c63ee8acb5p+0, -0x1.878c259fb524dp-2), q = c(0x1.c3c7850a545a2p-1, 0x1.c923ed8b256a4p-1), x0 = c(0x1.aae3455cabec5p+0, 0x1.74c2844b6fc96p-1, -0x1.6ef97be253c4ep-1), xl = c(0x1.83da136050932p-1, -0x1.758c7279dd54dp+0, -0x1.c6fc767c18f13p+0), xu = c(0x1.b8d4f2a314da2p+0, 0x1.7670bb75c82acp+0, 0x1.b6a38a294feap-4), npt = 7L, rhobeg = 0x1p+0, x = c(0x1.83da136050932p-1, -0x1.79050b1c66628p-1, -0x1.fa1b37ae0c7bep-2), fx = 0x1.6b3efcc6deed1p+0, nf = 46L, status = 0L),
    list(w = c(0x1.548b534756a78p-1, 0x1.15e34da7606f9p+1, 0x1.2eee660d7f508p+0, 0x1.cdae0048a915ep+0), c = c(-0x1.9db2253af1ca1p-3, 0x1.843f4f6a51fc8p+0, -0x1.103a428184043p-4, -0x1.a1861370b50bep+0), q = c(0x1.ceb877c77ea8cp-2, 0x1.2a46d3de607bp-2, 0x1.8049291a4484p-3), x0 = c(-0x1.9681322554a27p-1, -0x1.596677cf620bcp-3, -0x1.20df773773f03p+1, 0x1.fd98795bfba4p-5), xl = c(-0x1.7a07809d5a613p+0, -0x1.5edb88deeee4cp-3, -0x1.33737c5590e8ep+1, -0x1.af1787bda53a4p-1), xu = c(0x1.fafaa8690350cp-1, 0x1.0d3e6351675d9p-2, -0x1.409de0cf2b526p+0, 0x1.e8fa6590d5cabp+0), npt = 9L, rhobeg = 0x1.63bcec9a4bd99p-4, x = c(-0x1.b8d89ed75fee9p-4, 0x1.0d3e6351675d9p-2, -0x1.409de0cf2b526p+0, -0x1.af1787bda53a4p-1), fx = 0x1.a5b305c38975fp+2, nf = 46L, status = 0L),
    list(w = c(0x1.68bddc35e20dp+0, 0x1.45af15f2290ecp+0, 0x1.192d51ccd14d3p+1, 0x1.30bbac1297b8bp+1), c = c(0x1.1cbea31396ed6p+0, 0x1.3853774565b91p-1, 0x1.159bc07870177p-2, 0x1.6ed3ed20a27f4p-3), q = c(0x1.33fe9c6591b9p-1, 0x1.c02eb1de7ed7p-3, 0x1.ca509897ecccap-1), x0 = c(-0x1.8b4f97dab0f07p+1, 0x1.091b89126f2b5p+0, -0x1.887a0eeef5b36p+1, 0x1.27c99f470a6cfp-1), xl = c(-0x1.43a4f64e689p-1, -0x1.8effcb988e5afp+0, -0x1.30e4fce5bda8ap-2, -0x1.9646ab8e2bba2p-1), xu = c(0x1.6d8158f282389p+0, -0x1.84e122381d508p-1, 0x1.f7b0513e62fa6p-1, 0x1.38d174daeaf0bp+0), npt = 6L, rhobeg = 0x1.474b90c732b78p-3, x = c(0x1.865a7f0ec818cp-1, -0x1.84e122381d508p-1, 0x1.08add2868fb7dp-2, 0x1.c6eb2f1005eb6p-3), fx = 0x1.c5cdf72a651e5p+1, nf = 70L, status = 0L),
    list(w = c(0x1.d7bfcbf0a5d78p+0, 0x1.16b1f8df19944p-1), c = c(-0x1.1332ea12861f5p-1, -0x1.12633eb7c1948p-1), q = c(0x1.787cb7a4cced6p-1), x0 = c(-0x1.328f618a33465p+1, 0x1.0491639e050c4p+1), xl = c(0x1.a28668f356542p-1, -0x1.0660283209c1p-3), xu = c(0x1.acaf2c30c9b8dp+1, 0x1.0db9a5afb6138p+0), npt = 4L, rhobeg = 0x1.e40911232545dp-3, x = c(0x1.a28668f356542p-1, 0x1.174db9aa9dep-7), fx = 0x1.dc8ed805a5895p+1, nf = 20L, status = 0L),
    list(w = c(0x1.86f7ba467c8bp+0), c = c(-0x1.7a5bed474c041p-1), q = c(), x0 = c(-0x1.2f5fd7813fb2cp+1), xl = c(-0x1.36c92cb5b151cp+1), xu = c(-0x1.770c1b95476f6p+0), npt = 3L, rhobeg = 0x1.8a7062f02b86ap-3, x = c(-0x1.770c1b95476f6p+0), fx = 0x1.9c31054aa053ap-1, nf = 14L, status = 0L),
    list(w = c(0x1.897049317abcp+0), c = c(0x1.aad9d17877f18p+0), q = c(), x0 = c(0x1.c1381f89fdc2p-2), xl = c(-0x1.e60c32b4cce07p-1), xu = c(0x1.d4c57c787fd63p+0), npt = 3L, rhobeg = 0x1p+0, x = c(0x1.aad9d17877f18p+0), fx = 0x0p+0, nf = 14L, status = 0L),
    list(w = c(0x1.6c85034a8a4c4p+0, 0x1.039076f963e86p+0, 0x1.332bb3d52868bp+1, 0x1.165391e8c6ecfp+1), c = c(0x1.f514d5ba23b32p-2, 0x1.9b3ffe9675e43p-1, 0x1.607080d4e50a5p+0, -0x1.0acff5b8abe51p-1), q = c(0x1.9bb6f3009c9ep-4, 0x1.8fe6a7e1bd6a6p-1, 0x1.824ece0f00cdcp-2), x0 = c(0x1.c722b3df62bafp-1, 0x1.57f999c886aebp+0, -0x1.ed0e7cffe9569p-1, 0x1.1b0b932cec08p-4), xl = c(-0x1.a3261d3446feep-2, 0x1.b73b6596dae4p-2, -0x1.bc05e328bd5c5p+0, -0x1.253ac3a1935a3p+0), xu = c(0x1.3767948affe84p+0, 0x1.58ddfba636754p+1, -0x1.aba956f94006p-1, 0x1.37a79a36abea4p-1), npt = 9L, rhobeg = 0x1p+0, x = c(0x1.fdf60c4440736p-2, 0x1.b73b6596dae4p-2, -0x1.aba956f94006p-1, -0x1.0f18980493012p-1), fx = 0x1.8f02562ee9d37p+3, nf = 37L, status = 0L),
    list(w = c(0x1.3f2aef5343fc2p+0), c = c(-0x1.95eabdcdecc3cp-3), q = c(), x0 = c(-0x1.45cf1812db074p-2), xl = c(-0x1.67e5772f42c94p-2), xu = c(0x1.c462caf2f802cp-4), npt = 3L, rhobeg = 0x1.7a64ee566708p-4, x = c(-0x1.95eabdcdecc3cp-3), fx = 0x0p+0, nf = 13L, status = 0L),
    list(w = c(0x1.54f4cdc39edf8p-1, 0x1.81ccefc06873ep+0, 0x1.d17b0da5e325cp-1, 0x1.26dd946bc4df9p+1), c = c(0x1.30d6fc40ea1fbp-3, -0x1.259fe997cdc6p-1, 0x1.a0a6f3070138dp-2, 0x1.3a606b14eb285p-1), q = c(0x1.c296303bfc7cap-1, 0x1.4b38e763f7668p-3, 0x1.93806b61a7decp-1), x0 = c(0x1.75249bdfedc8bp-1, 0x1.05e6bb7ee2ddfp+0, -0x1.0e3c112b1fba5p-2, -0x1.1d1b7a384e397p+0), xl = c(-0x1.7d421241532dp-4, 0x1.bc8d4dabcf12p-1, -0x1.1c682096ec142p-1, -0x1.576640b241d1p+0), xu = c(0x1.57bd7fccae274p+1, 0x1.5a1f7642ccbbbp+0, 0x1.0780124cb36e7p-2, 0x1.fc34ff9a090dp-2), npt = 6L, rhobeg = 0x1.8c4f648faa3bdp-4, x = c(0x1.72a4e16c00b69p-2, 0x1.bc8d4dabcf12p-1, 0x1.0780124cb36e7p-2, 0x1.fc34ff9a090dp-2), fx = 0x1.aefb8985dba88p+1, nf = 31L, status = 0L),
    list(w = c(0x1.5f3235f4ffe02p+0), c = c(-0x1.4b6bb9148e74p-8), q = c(), x0 = c(0x1.7c57046cbd579p-1), xl = c(-0x1.1a5bf4ab48fa9p-1), xu = c(0x1.9adb8500fd16ap-1), npt = 3L, rhobeg = 0x1p+0, x = c(-0x1.4b6bb9148e8p-8), fx = 0x1.8b187cb39fdc2p-105, nf = 12L, status = 0L)
  )
  for (r in sweep) {
    fun <- function(x) pp_sweepf(x, r$w, r$c, r$q)
    res <- suppressWarnings(.bq_bobyqa(fun, r$x0, xl = r$xl, xu = r$xu, rhobeg = r$rhobeg, npt = r$npt))
    expect_true(is.finite(res$f))
    if (ref_platform()) {
      expect_identical(res$x, r$x)
      expect_identical(res$f, r$fx)
      expect_identical(res$nf, r$nf)
      expect_identical(res$info, r$status)
    }
  }
})

test_that("LINCOA matches PRIMA.jl on a seeded sweep", {
  sweep <- list(
    list(w = c(0x1.7eba1ccb9a742p+0, 0x1.5084e8880df6ep+0, 0x1.2ef9cb45718a4p+1, 0x1.b3ae8ed774c64p+0), c = c(-0x1.5888c65df4e76p+0, 0x1.c8239293cd721p-1, -0x1.7ab36553c1013p-1, -0x1.98db5ecd13fe2p-2), q = c(0x1.223efa0121274p-2, 0x1.bb4287363a1d8p-3, 0x1.eceda924a06cp-6), x0 = c(0x1.303ad6a80097fp+1, 0x1.5274ad8c51bf3p-2, 0x1.79a3f143c210ap-1, -0x1.47c7ca304eb36p-1), Aineq = matrix(c(0x1.01689be9c9713p+0, 0x1.4cc678ef804acp+1, -0x1.8e44e5524463p-1, 0x1.f4c0cd6dc93aap-1), 1, 4), bineq = c(0x1.4ecdc9beae4bcp+1), Aeq = matrix(c(-0x1.43a51ccfc56a8p+0, -0x1.c7e1ad1b49e7p-4, 0x1.66e7f75e5c996p+0, -0x1.3895ffe8b0246p+1), 1, 4), beq = c(-0x1.e5ec9b97aebf1p-2), xl = c(0x1.48159fa1e48aap+0, -0x1.eb659a4f183a4p-1, 0x1.a5627c0b40638p-3, -0x1.c3d25e4742cbap+0), xu = c(0x1.b9a7d06784715p+1, 0x1.afdfd7196a399p+0, 0x1.d4f00164c3afbp+0, 0x1.8861aa262089p-1), npt = 6L, x = c(0x1.46303cbbebb2ep+0, 0x1.6f97a0d0b9c4dp-1, 0x1.a6784d40def63p-3, -0x1.8fc8814604fa2p-2), fx = 0x1.902e829276315p+3, nf = 29L, status = 0L, cstrv = 0x1.b070f9b1c142p-6),
    list(w = c(0x1.74910c8dc325ap+0, 0x1.039c2e91f4c3fp+1, 0x1.b03b943bca4ecp+0), c = c(0x1.069c1f3c67e12p-5, 0x1.65ac113787c3ap+0, -0x1.21c5dc972a559p+1), q = c(0x1.e9b24784cd1c8p-1, 0x1.bb3fe10b498p-8), x0 = c(-0x1.090d7eabe9229p-1, 0x1.0b35bb51ea15fp-3, 0x1.7e749939f192bp-1), Aineq = matrix(c(0x1.fbd398e6c783p-2, -0x1.cb8eac004c74ep-1, -0x1.0a38f13b71f3bp+0), 1, 3), bineq = c(-0x1.4f290e4c8fb7bp+0), Aeq = matrix(c(-0x1.3e9ef74d4feb3p-1, -0x1.4ebe0755861bdp+0, 0x1.903d938438486p-1), 1, 3), beq = c(0x1.a5a70b670e80ep-1), xl = NULL, xu = NULL, npt = 7L, x = c(-0x1.8c802faee396p+0, 0x1.0b2d82c2db778p-2, 0x1.3cdac38008a74p-3), fx = 0x1.0ea4b1cf784bp+4, nf = 34L, status = 0L, cstrv = 0x1.2b13f6cdf348p-3),
    list(w = c(0x1.4dea2b3d4e798p+0, 0x1.4637c2c162c2p-1, 0x1.e58fa6191304cp-1, 0x1.0c40897f52d6cp-1, 0x1.5d60e3c5f5f26p+0), c = c(-0x1.4687fff37b678p-1, -0x1.31a49cdbd3d56p-2, -0x1.b6d4b092b3a44p-1, 0x1.b51e4401784dbp-1, 0x1.c739fac049f8fp-1), q = c(0x1.222195570f9fcp-2, 0x1.af444990d49ap-2, 0x1.2c476db44b602p-1, 0x1.67b33416b513p-3), x0 = c(0x1.274950ddfa97bp-3, -0x1.3e578c0e4b1c8p+0, -0x1.fed1887a4478bp-2, -0x1.1c4b0236e4c67p+0, -0x1.6607255d33c04p+0), Aineq = NULL, bineq = NULL, Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 7L, x = c(-0x1.50cf85b9cc91cp-1, -0x1.c54fd70fed483p-2, -0x1.80b27431403c9p-1, 0x1.3a39395200a6p-2, 0x1.cbc3bee2ceb28p-1), fx = 0x1.0ded407e23743p-1, nf = 99L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.04966f4cbb9acp-1, 0x1.2b0b6087a640ap+0, 0x1.cbb788f08881ap+0, 0x1.c05b2b79ec3e8p+0, 0x1.a631dbba58582p+0), c = c(-0x1.1c62ed77845efp-1, -0x1.da71c88cacceep-1, 0x1.a3dbef34223b1p-4, -0x1.df7b1cabe1befp-2, -0x1.6d018ec81113bp-3), q = c(0x1.edc9be64916cp-4, 0x1.66e618b208afp-1, 0x1.63710fa287b5ap-1, 0x1.3e30419ece34p-2), x0 = c(-0x1.4e63007cf3a3p-2, -0x1.7bc8b474222d5p-1, 0x1.055bda193e614p+1, -0x1.1b11bb8200d1fp+0, 0x1.6d1cf0657dcep+0), Aineq = matrix(c(-0x1.31139ecb39edap+1, -0x1.02b6f2cc89dd8p-1, -0x1.9707dbc0a9416p-2, 0x1.f55b84a3a79b2p-1, -0x1.c1c8d4849ad2fp-2, -0x1.b2996bec6563bp-5, 0x1.daaeca0a059e7p+0, -0x1.7cd898bfbfa5p-2, -0x1.44bf42d1f73adp-5, -0x1.7bccb9f33a14ap-1), 2, 5), bineq = c(-0x1.d82c979076ce8p+0, -0x1.341e8ef29d334p-1), Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 21L, x = c(0x1.7aa01f39da58dp-2, -0x1.bb4203f9238c7p-1, -0x1.9e2fd7f7844e3p-4, -0x1.778547a582866p-1, -0x1.aaef5aea1b9fbp-3), fx = 0x1.00b6e07bbcad8p+0, nf = 78L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.38e3658e17feap+0, 0x1.54e54cab5f90ep+0), c = c(0x1.4325fdaa15ff2p+0, -0x1.10407dde4d62p-2), q = c(0x1.b054427cbbdacp-1), x0 = c(-0x1.2a6f8ab92fe6p+0, 0x1.1638c3be35032p+0), Aineq = NULL, bineq = NULL, Aeq = matrix(c(-0x1.0a18d8cb0b6a6p-1, 0x1.19dca8216fc75p-5), 1, 2), beq = c(0x1.28a0944615c85p-1), xl = c(-0x1.4603ac1c3138cp+1, 0x1.112957a94d1fap-1), xu = c(0x1.a6525adfe1a68p-3, 0x1.39039db8c228p+1), npt = 6L, x = c(-0x1.145467ac0ab06p+0, 0x1.112957a94d1fap-1), fx = 0x1.10fa1b2708d58p+3, nf = 26L, status = 0L, cstrv = 0x1p-53),
    list(w = c(0x1.6fe39d24a43b4p-1, 0x1.1b58be468f5d2p+1, 0x1.423caf8542058p-1, 0x1.377f873a2348cp+1, 0x1.6d38d9e87bee2p+0), c = c(0x1.911a8fa89b541p-4, 0x1.15ffcf0fea4ap-3, 0x1.368edbd4d234ap-2, -0x1.4ac483904defap-3, 0x1.0b800606bb87ap-1), q = c(0x1.978554a12757cp-2, 0x1.318746031c06p-4, 0x1.c568c94be542p-4, 0x1.345329a3222f4p-2), x0 = c(0x1.263f938f52aa4p+0, -0x1.d549e1dd3aba5p-2, -0x1.6138f2fb216a7p-2, -0x1.eb81dffd5773fp-1, -0x1.b39de7bc60554p+0), Aineq = NULL, bineq = NULL, Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 21L, x = c(0x1.1b4371a794f67p-3, 0x1.37cd92483d904p-3, 0x1.3379f8626eaecp-2, -0x1.e729e4b55528fp-4, 0x1.044c00928073bp-1), fx = 0x1.ea8230e256ee6p-3, nf = 70L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.3f3ae87dd662ep+0, 0x1.19f3efc7150fdp+1), c = c(-0x1.84040133597c9p+0, 0x1.9ba1cdd4d849bp-3), q = c(0x1.cdd835a5693c8p-1), x0 = c(0x1.e737ac79d707cp+0, -0x1.c05fd974f4a54p-1), Aineq = matrix(c(-0x1.8ece75467c95bp+0, -0x1.0ab79b8b690d4p+0, -0x1.16b7e7efa1f1cp+0, -0x1.64b768cd58b19p+0, 0x1.84211624cb65bp+0, 0x1.3feff15f4a16p-3), 3, 2), bineq = c(-0x1.a4cd66df8d54ep+0, -0x1.631af3e74ede2p+1, -0x1.33ca344cbf72cp+1), Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 6L, x = c(0x1.739bd63ceb83ep+1, 0x1.fd0134585ad6p-4), fx = 0x1.8613d962a3fb7p+4, nf = 27L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.7fd3da719c05p-1, 0x1.30351fa4b435p+1), c = c(-0x1.1c8516ad7e444p-2, -0x1.6e4bb248c16c9p+1), q = c(0x1.bd02769c898cp-2), x0 = c(0x1.0617e885ca8bdp-1, -0x1.2c93a7fe7d274p+0), Aineq = NULL, bineq = NULL, Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 4L, x = c(-0x1.8ace3e8b1b2e4p-3, -0x1.6e106c1f7941fp+1), fx = 0x1.af19e62ff1136p-8, nf = 35L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.100060efd07f2p+1, 0x1.1ce290d494a36p+0), c = c(0x1.4b67f424ae832p-3, 0x1.a3b9178fab6e4p-1), q = c(0x1.5290f3a34253ep-1), x0 = c(0x1.f34b8fd1ab3fcp-1, -0x1.c4e200f162a48p+0), Aineq = NULL, bineq = NULL, Aeq = matrix(c(0x1.eb36afd1b7f1p-4, 0x1.8041b3e6299f9p-1), 1, 2), beq = c(-0x1.5a83fd722a8dfp+0), xl = NULL, xu = NULL, npt = 6L, x = c(-0x1.968e05ca9808p-6, -0x1.c4e200f162a48p+0), fx = 0x1.eadfdbc6e8619p+2, nf = 33L, status = 0L, cstrv = 0x1.775626f58178p-6),
    list(w = c(0x1.3eaef963f9c6cp+0, 0x1.291e718dad82dp+1), c = c(0x1.eb82d43139aedp-1, 0x1.2b98bba606858p+0), q = c(0x1.f54fb2f73194p-2), x0 = c(-0x1.80ce67681daddp-7, -0x1.f953efc53db8ep-4), Aineq = matrix(c(0x1.04ecdc9d6c201p-1, -0x1.f1f926d4a4e66p+0, 0x1.f97639f37367dp-2, 0x1.e9b366aa916eep-2, -0x1.537f3c04f1ce4p-2, -0x1.5faab2dbf322ap-2), 3, 2), bineq = c(-0x1.f0d24f79bd8eap-3, 0x1.4850f3691561p-5, -0x1.1e29ea43d87bep-6), Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 6L, x = c(-0x1.80ce67681daddp-7, -0x1.1124d92a4cf0fp-3), fx = 0x1.4f7f9c7e9b3c5p+2, nf = 25L, status = 0L, cstrv = 0x1.61f0eeda508p-3),
    list(w = c(0x1.ed83efbd412fcp+0, 0x1.44e35628d00b2p+0, 0x1.0dc3d32b956p-1), c = c(0x1.ebc17812bc631p-8, 0x1.c346b4b103e81p-1, -0x1.9b3eefe3a2dc6p-1), q = c(0x1.9d90fefb88994p-2, 0x1.2fdf19bcf8002p-1), x0 = c(0x1.21f128e917af9p+1, 0x1.7866a6904a669p-1, -0x1.aa1b4f0836d86p-2), Aineq = matrix(c(0x1.1b0861d2c3da9p-1, 0x1.19f0b18a1f0b3p+0, 0x1.0cb991fb593b8p-4), 1, 3), bineq = c(0x1.64e5d88237b7fp+1), Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 7L, x = c(0x1.539c2266f2a99p-4, 0x1.abe1d2a82c2bbp-1, -0x1.7d296d149085dp-3), fx = 0x1.16a45ced34732p-1, nf = 50L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.a7178fefaaee4p-1, 0x1.4d62dcbe79834p+0, 0x1.e74d1c50cdeccp+0, 0x1.e530db1f2bd34p-1), c = c(-0x1.f467b12f03f4ap-4, -0x1.ffd704c5ccb22p-4, -0x1.388b3157e05cbp-2, 0x1.f55bdbfd92065p-2), q = c(0x1.9960e221a7f64p-1, 0x1.6558c32d64aecp-1, 0x1.27eba42108538p-3), x0 = c(-0x1.b5a8ba1f6a12p-1, 0x1.81deb33c18f2p-2, -0x1.10d14162f5071p-1, 0x1.14a56b9419e84p-2), Aineq = matrix(c(-0x1.df0f9355fb427p-3, 0x1.f589a39f38503p-2, 0x1.70631aaf0abafp+0, 0x1.98f4751977fa8p-2), 1, 4), bineq = c(0x1.33e26ed941368p-2), Aeq = NULL, beq = NULL, xl = c(-0x1.ac6e47b0bbb0fp+0, -0x1.65a7073d5ca4p-2, -0x1.f1bd2998d0948p+0, -0x1.66e22ee286766p-1), xu = c(-0x1.a2a55522d54cp-6, 0x1.d34ff33dc932dp+0, 0x1.a2bb4e23460ep-6, 0x1.afd9d648e1b48p-1), npt = 6L, x = c(-0x1.dc66393724df2p-3, -0x1.09e53448ff13cp-2, -0x1.4a165c5418f29p-2, 0x1.d4c438115fe88p-2), fx = 0x1.7c3366b4bd20cp-2, nf = 71L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.5aaa60f182174p+0, 0x1.e5e8aaa416216p+0, 0x1.99d5f287ad3d8p-1, 0x1.0ee83e0234765p+1, 0x1.112704195aff8p+1), c = c(-0x1.7989faae145efp-2, -0x1.0db04febddbebp+0, 0x1.851c0b0447f48p-3, 0x1.2aafa5d8b482fp-1, 0x1.38de9156dac15p-1), q = c(0x1.ff73a8e91162p-3, 0x1.896638c5ac58p-3, 0x1.26e4051f0b7cep-1, 0x1.89de7ebbbf7a8p-1), x0 = c(0x1.bc512ae976e2ep-5, 0x1.d1610ed8abe26p-2, -0x1.fae65640586f6p+0, 0x1.56bc50452165fp+1, -0x1.cc8820fb9ceb1p-3), Aineq = NULL, bineq = NULL, Aeq = NULL, beq = NULL, xl = c(-0x1.27182cf7a8ab7p-1, -0x1.06f0bfa2b7532p+0, -0x1.837ee33b2f552p+1, 0x1.f33e765c1c9ccp+0, -0x1.f1e03ea6ab214p-1), xu = c(0x1.57d41b110e831p-1, 0x1.c6e3427b458c2p+0, -0x1.548c7d4b3a613p+0, 0x1.ad55327b61db4p+1, 0x1.0d6d12df9b6ddp+0), npt = 21L, x = c(-0x1.9054654258625p-2, -0x1.e9807e4978266p-1, -0x1.548c7d4b3a613p+0, 0x1.f33e765c1c9ccp+0, 0x1.9fb9aff76777bp-2), fx = 0x1.738a8510c901p+3, nf = 64L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.e35016d14110cp+0, 0x1.0513ea2adb318p+0, 0x1.d1605230fa0fcp-1, 0x1.3363e2ec48b22p+0), c = c(0x1.db36c0e3fd338p-1, 0x1.c737ad24ddecp-1, 0x1.b41fb30d2ffdep-1, 0x1.8c28c3c39652cp-1), q = c(0x1.b1f0f14293498p-1, 0x1.dfe2b4b7e1f84p-1, 0x1.2b34a2d24dcd8p-3), x0 = c(-0x1.3346090428237p-1, -0x1.15ed7a2f52806p-3, -0x1.5eededbf6ec84p-3, -0x1.3fc5f66be32dbp+0), Aineq = matrix(c(0x1.1f5db2653626dp-1, 0x1.0e2bb5ced8c74p+0, -0x1.d557310f90c7ap-3, 0x1.2066c4a1319a7p+0, -0x1.4259fe574c178p-1, -0x1.bca2f6504333ep-2, 0x1.1afc52388cep+0, 0x1.e82e4fb3d5746p+0), 2, 4), bineq = c(-0x1.3f0b08921784bp+0, -0x1.a04066a16a452p+1), Aeq = NULL, beq = NULL, xl = c(-0x1.2436f274713b9p+0, -0x1.9ee92fc8361b6p-1, -0x1.9f191f7c9be2ap+0, -0x1.d5b1aba27fdf4p+0), xu = c(0x1.b424ec46e63a5p-1, 0x1.9212f28b77dep-1, 0x1.edc1b9862be73p-1, 0x1.be5f5d086f61p-3), npt = 15L, x = c(0x1.aabf19a128ecp-5, -0x1.5c775921b6684p-3, 0x1.e60eb0b285a9bp-1, -0x1.758b140379d5ap+0), fx = 0x1.379936d7ce15fp+3, nf = 50L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.1d99528ba8e9p-1, 0x1.b2c71ce67ff5cp-1, 0x1.f724fbeb8273ep+0, 0x1.d52088f4909f4p+0, 0x1.358187d163e0ep+1), c = c(0x1.cbe61e2c314f8p-1, 0x1.5cb4302bfe59dp+0, -0x1.f29bf5b19687fp-2, -0x1.fbd2ff49c1126p+0, -0x1.a2ab3835db3cap-1), q = c(0x1.43f45324b05fcp-1, 0x1.ebcde3e31fafap-1, 0x1.a8af5fba325fp-3, 0x1.68c83471779ecp-2), x0 = c(-0x1.c6154553abdd6p-1, -0x1.65823be81c1eap-1, 0x1.8d55a4c22cab2p-1, -0x1.1f20fbbc77c0ep+0, 0x1.a47b93b2bf4c1p+0), Aineq = matrix(c(0x1.5eb4e3bc17e74p-1, -0x1.d3584bf983d37p-1, -0x1.877132686c082p-7, 0x1.4af1abbbfb51cp-1, 0x1.e1c3c7f859ec2p-1), 1, 5), bineq = c(0x1.0b8ce5976d323p+0), Aeq = matrix(c(-0x1.424d23116852fp-1, -0x1.c6485a85ecba8p-4, -0x1.5dc612c4e323p+0, 0x1.9deb7acf2f5b5p+0, -0x1.57d4469f262cep-3), 1, 5), beq = c(-0x1.34429fac6d14bp+1), xl = NULL, xu = NULL, npt = 11L, x = c(0x1.927acbce2bffbp-2, 0x1.2336cce31570fp+0, -0x1.1413ee3be81bbp-2, -0x1.8fe41ba57ae9dp+0, -0x1.75f07f1db4c67p-1), fx = 0x1.6522610831088p+0, nf = 67L, status = 0L, cstrv = 0x1p-50),
    list(w = c(0x1.dd9ac060fece2p+0, 0x1.2b399e0dcdf1cp-1, 0x1.859797486af1cp+0, 0x1.ecae8ff8d3b1p+0), c = c(0x1.96019ed922783p-3, 0x1.176fb4d539a89p+1, -0x1.00bf724c74f36p-4, 0x1.7abb8d8497df5p+0), q = c(0x1.e2d27719b384p-5, 0x1.a60badb9502cp-4, 0x1.9436e3439594p-1), x0 = c(-0x1.6130f96ceb2bdp+0, 0x1.29896ca631965p-3, 0x1.a028cc91e9964p+0, 0x1.e7ecd693e538ap-2), Aineq = matrix(c(-0x1.aa5c91504dbe9p-2, 0x1.7102b5b9f2ad1p-1, -0x1.fd146acf8698fp-1, -0x1.4065820c4c446p+0, 0x1.7b60f8bc03578p-3, -0x1.7aa374ff666bfp-1, -0x1.05717f7b8ae77p-2, -0x1.7a0b0f2c697f1p-1), 2, 4), bineq = c(0x1.5d8993a2121c6p-1, -0x1.0478645420558p+1), Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 9L, x = c(0x1.9e242a9421227p-3, 0x1.181fc3c9f2f97p+1, 0x1.492a25117b00ep-3, 0x1.7f1cb52374a49p+0), fx = 0x1.0eb2e8100bd3fp-3, nf = 69L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.b89c208a38e24p-1, 0x1.fe4405f4aa058p-1, 0x1.8a043a4e5168cp+0, 0x1.2ff2e58bdbc69p+1, 0x1.975b504c0489ep+0), c = c(0x1.3d516ba6a980cp-2, 0x1.2766b92ba6734p+0, -0x1.a09a3cb33e429p+0, -0x1.183078e7405fbp-1, 0x1.4483030de8887p-1), q = c(0x1.4b117d06a3c54p-2, 0x1.a4ef27185d078p-2, 0x1.0f56a1ce53d5p-3, 0x1.be590d2e3875p-3), x0 = c(0x1.e80825a999c17p-1, -0x1.de66846d636b6p-1, 0x1.6f368cb07b545p-3, 0x1.91cd8d0f9ac41p+0, 0x1.2748f56fc3642p+0), Aineq = NULL, bineq = NULL, Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 7L, x = c(0x1.788aa144c2c2cp-2, 0x1.e26d9d4a1550dp-2, -0x1.77e41fde10862p+0, -0x1.f36cd2097a665p-2, 0x1.29b61ff3a1fe1p-1), fx = 0x1.470e9fb703458p+0, nf = 88L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.95f3097d908cp+0, 0x1.0f81ba7349db9p+1, 0x1.52c27dd70bc58p-1), c = c(0x1.cd44b1be8d649p-2, -0x1.06c4145db8fp-6, -0x1.d8144e9e8b47dp-3), q = c(0x1.597abe43a165ep-1, 0x1.b762d43c66a5ap-1), x0 = c(0x1.238f9f838c542p-3, -0x1.fc92f9e9f5ce8p-1, 0x1.237c2967417fdp-1), Aineq = matrix(c(0x1.941abe4921c7fp-1, 0x1.c17a30fbfff43p-1, 0x1.2c679992cb467p+0, -0x1.794e21dadc2cfp-3, -0x1.3aba0876d8e25p-1, 0x1.0425bdc81440bp+0, 0x1.d55e59fb9f159p-1, 0x1.a8add02a06452p+0, -0x1.f7797c4028367p-3), 3, 3), bineq = c(0x1.0d3f0e7a70a47p+0, 0x1.ce6be33e17f9ap+0, -0x1.038703ba8d488p+0), Aeq = NULL, beq = NULL, xl = c(-0x1.259ec993de77bp+0, -0x1.25a90b34c65cp+1, -0x1.10f4806eb481dp-1), xu = c(0x1.495c23ff31906p+0, -0x1.2cabc89a8909p-2, 0x1.f342a9547cf23p+0), npt = 10L, x = c(-0x1.cf95d2b715a77p-2, -0x1.3d1ad8dd7604bp-1, -0x1.608fc223d616p-4), fx = 0x1.243ed8c425314p+1, nf = 45L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.a68ce1c6814bp+0, 0x1.78896280b318p+0, 0x1.a4e5434bfd7cap+0, 0x1.a0a250d362d8ap+0, 0x1.61cbd34111362p+0), c = c(0x1.074013c11e8cbp-2, -0x1.50b533bc12b47p-5, -0x1.374ca11ba3571p-3, -0x1.3119b7cb22a35p-2, -0x1.85bb5fe505d05p+0), q = c(0x1.c32db9aa9ccd8p-2, 0x1.aa7eba82e68bcp-1, 0x1.143ebf6d57b4p-4, 0x1.a4b8c2543645ep-1), x0 = c(-0x1.ed64fc5bc3976p-4, -0x1.095ff3ae5e018p-4, 0x1.1f02f2121b837p-4, -0x1.896ff1d3fc039p-1, 0x1.70604f30f2fa1p+0), Aineq = NULL, bineq = NULL, Aeq = matrix(c(0x1.00abdec41326ap+0, -0x1.416326a521213p-2, 0x1.274ce516fadfp+1, 0x1.833fe0a003453p+0, 0x1.4528a59efd174p+0), 1, 5), beq = c(0x1.5f46f4b6457a6p-1), xl = c(-0x1.5cb836ebc35c8p+0, -0x1.410926045fb66p+0, -0x1.81270308539dbp-1, -0x1.5bca9e35a829dp+0, 0x1.314de9e0762acp-1), xu = c(0x1.6f46fb3410ed5p-1, 0x1.0d9662a99e9e4p+0, 0x1.06ef2d75697eep+0, 0x1.955dd2b8b73c4p-3, 0x1.f3eae88e767p+0), npt = 11L, x = c(0x1.1ee3c0f3893ep-2, -0x1.24ccdffeb3f14p-5, -0x1.92695c733d941p-4, -0x1.725281ffd8cb4p-4, 0x1.314de9e0762acp-1), fx = 0x1.b78932dcd1d8ep+2, nf = 48L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.9d5512733bb8p-1, 0x1.086b18139058bp+1), c = c(-0x1.41c1ef003517p-3, -0x1.ccae45db6c48dp+0), q = c(0x1.43fbbfa759ccep-1), x0 = c(-0x1.29c5c3eb5da5ep-1, -0x1.1e9783737ea5dp+0), Aineq = matrix(c(0x1.685e61e83d3d5p-1, -0x1.38c31fe62f9b6p+0), 1, 2), bineq = c(0x1.748c34e9dec97p+0), Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 5L, x = c(-0x1.d0e8187ffc0d9p-2, -0x1.73e4e612bdd8ap+0), fx = 0x1.57cfd24ee371ap-2, nf = 30L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.2027caeb52feep+1, 0x1.df521d53b38bcp-1, 0x1.363b9a2881d65p+1, 0x1.c6399dbb1554p+0, 0x1.d0b11b46cae74p-1), c = c(0x1.6ffd6c2a80586p+0, -0x1.56dcab13ded13p-3, -0x1.2876d9a541f51p+0, 0x1.706345ad4b611p-4, -0x1.1258369ef6aa8p+0), q = c(0x1.8825db9a78558p-3, 0x1.7cf1bb471c176p-1, 0x1.38ee9b27c5198p-1, 0x1.6cf6426dfba4p-1), x0 = c(0x1.fd39ad8791413p-5, 0x1.14d73acf35308p-2, -0x1.dd27b1636c4b1p-6, -0x1.642f1e359922bp+0, -0x1.b2732d46b695ap-1), Aineq = NULL, bineq = NULL, Aeq = NULL, beq = NULL, xl = c(-0x1.0a4724cd1127ep+0, -0x1.e8f551d65acd6p-1, -0x1.667954f80542fp+0, -0x1.023e607f45d7cp+1, -0x1.0dcf66b6f2f22p+1), xu = c(0x1.75b891bc3c9c9p-1, 0x1.8dd669ed8bf4bp+0, 0x1.74f351da53616p+0, -0x1.2e0fc802930cp-4, 0x1.43a375c5d44a8p-3), npt = 7L, x = c(0x1.75b891bc3c9c9p-1, -0x1.035c56fb9aba8p-2, -0x1.2fb65fe294d14p+0, -0x1.5d2413932a8efp-3, -0x1.1cfa5855b3de9p+0), fx = 0x1.7ffba32615661p+0, nf = 81L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.9f57514f2d39p+0, 0x1.303533a83c31p+1, 0x1.dcf5f784f298cp-1, 0x1.5359c9fcc85e2p+0), c = c(0x1.7bd80041201ffp+0, -0x1.deb6e5c31ddbep+0, 0x1.fdbea25939289p-4, 0x1.2b9fe6808b102p+0), q = c(0x1.a6db87dc2b59cp-1, 0x1.1af76b57275bcp-2, 0x1.da9202dae6d3p-1), x0 = c(0x1.3992946033d78p-3, -0x1.f0a1610aa24b4p-1, -0x1.09ed0746d791ap+1, 0x1.a0d112cb59ae9p-1), Aineq = matrix(c(-0x1.0743523351078p-1, -0x1.cbab0704b4c3dp-1, -0x1.da57c0aca98eep-5, -0x1.403233a5e41e9p+0, -0x1.26bb1b7ac1fc6p+0, 0x1.0558638d0e46cp-2, -0x1.b694feec87433p-2, 0x1.48e1b3767ea04p-1, -0x1.336b2c42e15d7p-3, -0x1.c6a870a035a24p+0, 0x1.4f26c537b4e0cp+0, 0x1.7cb2ca400c584p-1), 3, 4), bineq = c(0x1.5bb726b5671edp+0, 0x1.9de7c1e9263ccp-1, 0x1.6ee19ed2950aap+0), Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 6L, x = c(0x1.0783a6cbafd05p+0, -0x1.04925876b02a8p+0, -0x1.191d5a19c958dp-3, 0x1.fb85d57eaf7bdp-2), fx = 0x1.42a70487a9badp+2, nf = 72L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.3230a005884cep+0, 0x1.649933b1f0e22p+0, 0x1.75269761c1adp+0, 0x1.ab4ea9d74929cp+0, 0x1.d236bf40507ap+0), c = c(-0x1.122b8c81102eep-6, -0x1.572f81300d16fp+0, 0x1.ceccc35fa423cp-2, 0x1.855c45637f541p+0, 0x1.0303603fbda5fp-2), q = c(0x1.35dad22c572c8p-1, 0x1.5c3bc563a92d6p-1, 0x1.3d26bb2b9d976p-1, 0x1.80bdfcf95e486p-1), x0 = c(-0x1.33e03f24503c3p+1, 0x1.bde773ee3401ep-1, 0x1.29471bf7ccb44p-1, -0x1.3a2ea117e8671p+0, -0x1.a9151041b7047p+0), Aineq = matrix(c(0x1.57494d95abbe4p-1, -0x1.441ccaacad562p+0, -0x1.6feb897f2f5bap+0, 0x1.dbe830a8089e8p-1, 0x1.09011aa039be9p+0), 1, 5), bineq = c(-0x1.9382655ba559cp+2), Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 11L, x = c(-0x1.4ac7e1ff98737p+0, 0x1.d04c1498ecfa2p-3, 0x1.6caa4b81479f6p+1, 0x1.8d55b3c6d6e46p-4, -0x1.1bb22b423f2d4p+0), fx = 0x1.53421c2597066p+4, nf = 61L, status = 0L, cstrv = 0x0p+0),
    list(w = c(0x1.7ffdb8bf2497ep+0, 0x1.17e4a1f6b152ap+0, 0x1.2c6b0b41e9e4cp+1, 0x1.279a53070fa6bp+1, 0x1.06f9a66ef5811p+1), c = c(-0x1.429d7bd3dc7c2p+0, 0x1.e95f21508f313p-1, -0x1.e48f02debdf7fp-4, -0x1.69c944ea8e7d2p-1, 0x1.7314463c799ap+0), q = c(0x1.7c03c70359e4p-2, 0x1.2d1fec8a7897ap-1, 0x1.a47e54a6731ap-5, 0x1.d5caf816a984p-4), x0 = c(-0x1.6d404f23ce68p-3, -0x1.e7508f38dbb3bp-3, 0x1.20141a6c0a64ap+0, -0x1.25c9cc5e117d1p+0, 0x1.3760d16371e4cp-2), Aineq = matrix(c(0x1.6d1c50e41d3ccp-2, 0x1.70e6e51346189p+0, -0x1.5e87ef713531ep-1, 0x1.4e29597ea477ap-1, 0x1.1937e5710b162p+0, 0x1.e215aed6b1bd3p-1, -0x1.2a84e6b11d8dep-1, -0x1.33c9993cbcdfcp+0, 0x1.264dc404597d7p-2, -0x1.72b34bc8a697ap+1), 2, 5), bineq = c(0x1.4a8aea02bbebfp+1, 0x1.10f415994fbfap+0), Aeq = NULL, beq = NULL, xl = NULL, xu = NULL, npt = 21L, x = c(-0x1.1e45ccb1b0b05p+0, 0x1.0dbd454d4221bp-1, -0x1.c05ef9b578db4p-5, -0x1.39948645edcdap-1, 0x1.672e57dbaed29p+0), fx = 0x1.17f66e9728a37p+0, nf = 66L, status = 0L, cstrv = 0x0p+0)
  )
  for (r in sweep) {
    fun <- function(x) pp_sweepf(x, r$w, r$c, r$q)
    res <- suppressWarnings(.lc_lincoa(fun, r$x0, Aineq = r$Aineq, bineq = r$bineq, Aeq = r$Aeq, beq = r$beq,
                                       xl = r$xl, xu = r$xu, npt = r$npt))
    expect_true(is.finite(res$f))
    if (ref_platform()) {
      expect_identical(res$x, r$x)
      expect_identical(res$f, r$fx)
      expect_identical(res$nf, r$nf)
      expect_identical(res$info, r$status)
      expect_identical(res$cstrv, r$cstrv)
    }
  }
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
