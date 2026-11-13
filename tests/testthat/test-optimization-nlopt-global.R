# The DIRECT family (cdirect.c and the DIRECT 2.0.4 translation), PRAXIS,
# StoGO and AGS of NLopt 2.11.0. The references are the runs of NLopt_jll
# through NLopt.jl in tools/nlopt-probe.jl (the runghex section, srand(1234)):
# the point, the minimum, the return code and the number of evaluations,
# bit for bit on the reference platform.

gl_poly3 <- function(x) {
  a <- x[1] - 1; b <- x[2] - 2; c <- x[3] + 1
  a * a + 2 * b * b + 3 * c * c + a * a * b * b + 0.5 * a * c
}
gl_poly3g <- function(x) {
  a <- x[1] - 1; b <- x[2] - 2; c <- x[3] + 1
  c(2 * a + 2 * a * b * b + 0.5 * c, 4 * b + 2 * a * a * b, 6 * c + 0.5 * a)
}
gl_rosen <- function(x) { a <- 1 - x[1]; b <- x[2] - x[1] * x[1]; a * a + 100 * b * b }
gl_roseng <- function(x) { a <- 1 - x[1]; b <- x[2] - x[1] * x[1]; c(-2 * a - 400 * x[1] * b, 200 * b) }
gl_quart5 <- function(x) {
  s <- sum(x)
  v <- 0.1 * s * s * s * s
  for (i in 1:5) { d <- x[i] - 0.1 * i; v <- v + i * d * d }
  v
}
gl_quart5g <- function(x) { s <- sum(x); vapply(1:5, function(i) 2 * i * (x[i] - 0.1 * i) + 0.4 * s * s * s, 0) }
gl_ball <- function(x) sum(x * x) - 4

# one probe case: the algorithm, the problem, the box, the start and the
# settings of runghex
gl_run <- function(alg, prob, lb, ub, x0, cons = FALSE, xtol_rel = 0, ftol_rel = 0, maxeval = 0,
                   params = list(), upstream_replaceinf = FALSE) {
  f <- switch(prob, poly3 = gl_poly3, rosen = gl_rosen, quart5 = gl_quart5)
  g <- switch(prob, poly3 = gl_poly3g, rosen = gl_roseng, quart5 = gl_quart5g)
  o <- .nlopt_create(alg, length(x0))
  o$lb <- lb; o$ub <- ub
  o$xtol_rel <- xtol_rel; o$ftol_rel <- ftol_rel; o$maxeval <- as.integer(maxeval)
  o$params <- params
  o$upstream_replaceinf <- upstream_replaceinf
  o$f <- f
  o$fg <- function(x) list(val = f(x), grad = g(x))
  if (cons) o$fc <- list(list(m = 1L, f = gl_ball, tol = 1e-8))
  .nlopt_srand(1234)
  .nlopt_optimize(o, x0)
}

gl_ref <- list(
  gn_direct_rosen = list(x = c("0x1.092b59782345ap+0", "0x1.12d88d6720308p+0"), f = "0x1.5d18c0addb99cp-10", ret = "MAXEVAL_REACHED", numevals = 500L),
  gn_direct_l_quart5 = list(x = c("-0x1.c5f77123f642p-4", "0x1.82c03470d764p-4", "0x1.d80dc4e8a2fap-3", "0x1.63562ca5511ap-2", "0x1.d3b5b1090c458p-2"), f = "0x1.ac58109ac3746p-3", ret = "MAXEVAL_REACHED", numevals = 800L),
  gn_direct_l_rand_poly3 = list(x = c("0x1.00000381f891p+0", "0x1.00000381f8914p+1", "-0x1.00000381f8914p+0"), f = "0x1.1afdacafc825dp-41", ret = "MAXEVAL_REACHED", numevals = 600L),
  gn_direct_noscal_rosen = list(x = c("0x1.ff4c33f8fa07ap-1", "0x1.fe91bf3410488p-1"), f = "0x1.21cfe04eb336ep-19", ret = "MAXEVAL_REACHED", numevals = 500L),
  gn_direct_l_noscal_rosen = list(x = c("0x1.ff4c33f8fa07ap-1", "0x1.fe91bf3410488p-1"), f = "0x1.21cfe04eb336ep-19", ret = "MAXEVAL_REACHED", numevals = 500L),
  gn_direct_l_rand_noscal_rosen = list(x = c("0x1.ff4c33f8fa07ap-1", "0x1.fe91bf3410488p-1"), f = "0x1.21cfe04eb336ep-19", ret = "MAXEVAL_REACHED", numevals = 500L),
  gn_direct_magic = list(x = c("0x1.097b425ed0978p+0", "0x1.13506ac1242b6p+0"), f = "0x1.67980e0bf07dp-10", ret = "MAXEVAL_REACHED", numevals = 400L),
  gn_direct_l_fixed = list(x = c("0x1.00000381f891p+0", "0x1p+1", "-0x1.00000381f8914p+0"), f = "0x1.5882bbf010e9dp-43", ret = "MAXEVAL_REACHED", numevals = 400L),
  gn_direct_l_ftol = list(x = c("0x1.00004e9933876p+0", "0x1.00009e1b35b96p+0"), f = "0x1.87459846bd762p-36", ret = "FTOL_REACHED", numevals = 1423L),
  gn_orig_direct_rosen = list(x = c("0x1.092b59782345ap+0", "0x1.12d88d6720308p+0"), f = "0x1.5d18c0addb99cp-10", ret = "MAXEVAL_REACHED", numevals = 515L),
  gn_orig_direct_l_rosen = list(x = c("0x1.001df75680fe8p+0", "0x1.00389a4e1019ep+0"), f = "0x1.f7ee17cd2fd47p-22", ret = "MAXEVAL_REACHED", numevals = 517L),
  gn_orig_direct_xtol = list(x = c("0x1.5555555555554p+0", "0x1.c71c71c71c71cp+0"), f = "0x1.c71c71c71c70ep-4", ret = "XTOL_REACHED", numevals = 121L),
  gn_orig_direct_l_xtol = list(x = c("0x1.5555555555554p+0", "0x1.c71c71c71c71cp+0"), f = "0x1.c71c71c71c70ep-4", ret = "XTOL_REACHED", numevals = 83L),
  gn_orig_direct_magic = list(x = c("0x1.097b425ed0978p+0", "0x1.13506ac1242b6p+0"), f = "0x1.67980e0bf07dp-10", ret = "MAXEVAL_REACHED", numevals = 403L),
  gn_orig_direct_cons = list(x = c("0x1.ba781948b0fep-1", "0x1.7b425ed097b4p+0", "-0x1.f9add3c0ca458p-1"), f = "0x1.1f19a7a1d5aa8p-1", ret = "MAXEVAL_REACHED", numevals = 677L),
  gn_orig_direct_l_cons = list(x = c("0x0p+0", "0x0p+0", "0x0p+0"), f = "0x1.fp+3", ret = "MAXEVAL_REACHED", numevals = 633L),
  ln_praxis_poly3 = list(x = c("0x1p+0", "0x1p+1", "-0x1p+0"), f = "0x0p+0", ret = "XTOL_REACHED", numevals = 52L),
  ln_praxis_poly3_box = list(x = c("0x1.fffffffffffecp-1", "0x1.fffffffffff16p+0", "-0x1.0000000000058p+0"), f = "0x1.0452p-87", ret = "XTOL_REACHED", numevals = 108L),
  ln_praxis_poly3_half = list(x = c("0x1.ffffffffffc3p-1", "0x1.fffffffffff16p+0", "-0x1.ffffffffff8dfp-1"), f = "0x1.402ff6p-83", ret = "XTOL_REACHED", numevals = 393L),
  ln_praxis_poly3_upper = list(x = c("0x1.000000029a4e2p+0", "0x1.8p+0", "-0x1.00000000bc01cp+0"), f = "0x1p-1", ret = "XTOL_REACHED", numevals = 123L),
  ln_praxis_rosen = list(x = c("0x1.0000000000001p+0", "0x1.0000000000002p+0"), f = "0x1p-104", ret = "SUCCESS", numevals = 191L),
  ln_praxis_rosen_box = list(x = c("0x1.999999998d963p-1", "0x1.47ae14e21f54ap-1"), f = "0x1.47ae147b2eacap-5", ret = "FTOL_REACHED", numevals = 792L),
  ln_praxis_rosen_max = list(x = c("-0x1.949134a63ba8ap-1", "0x1.436574a203c3ap-1"), f = "0x1.9ae0d8b94f2bbp+1", ret = "MAXEVAL_REACHED", numevals = 25L),
  ln_praxis_quart5 = list(x = c("-0x1.c6e065ac8e2bcp-4", "0x1.82f6336cd559ap-4", "0x1.d652114d098ap-3", "0x1.6391f9c503648p-2", "0x1.d4c6b303d5a72p-2"), f = "0x1.ac53efc01102bp-3", ret = "XTOL_REACHED", numevals = 221L),
  ln_praxis_t0 = list(x = c("0x1.0000000000148p+0", "0x1.0000000000292p+0"), f = "0x1.a5dp-88", ret = "SUCCESS", numevals = 174L),
  gd_stogo_rosen = list(x = c("0x1.0000023957a69p+0", "0x1.0000049ca8a15p+0"), f = "0x1.e89ad293bc65dp-46", ret = "SUCCESS", numevals = 600L),
  gd_stogo_rand_rosen = list(x = c("0x1.ffffe5c26387dp-1", "0x1.ffffcb81e05bcp-1"), f = "0x1.584af48baae47p-41", ret = "SUCCESS", numevals = 600L),
  gd_stogo_quart5 = list(x = c("-0x1.c6e5c0f659605p-4", "0x1.82ffd379f362ep-4", "0x1.d650c0377f4cfp-3", "0x1.6391eeeaac1cbp-2", "0x1.d4c768c100287p-2"), f = "0x1.ac53efc905e35p-3", ret = "SUCCESS", numevals = 800L),
  gn_ags_rosen = list(x = c("0x1.a9cp-1", "0x1.584p-1"), f = "0x1.09ab105fe4p-4", ret = "MAXEVAL_REACHED", numevals = 500L),
  gn_ags_cons = list(x = c("0x1.c02p-1", "0x1.9cdp+0", "-0x1.716p-1"), f = "0x1.112e2f239452p-1", ret = "MAXEVAL_REACHED", numevals = 132L)
)

gl_check <- function(name, r) {
  e <- gl_ref[[name]]
  if (ref_platform()) {
    # off the reference platform the last bits differ, which can change which
    # stopping test fires first (PRAXIS stops with SUCCESS on Linux)
    expect_identical(r$ret, e$ret, info = name)
    expect_identical(r$x, as.numeric(e$x), info = name)
    expect_identical(r$minf, as.numeric(e$f), info = name)
    expect_identical(r$numevals, e$numevals, info = name)
  } else {
    expect_true(is.finite(r$minf), info = name)
  }
}

B3 <- list(rep(-5, 3), rep(5, 3)); B2 <- list(c(-2, -2), c(2, 2))

test_that("the DIRECT algorithms of cdirect.c follow NLopt", {
  gl_check("gn_direct_rosen", gl_run("GN_DIRECT", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), xtol_rel = 1e-8, maxeval = 500))
  gl_check("gn_direct_l_quart5", gl_run("GN_DIRECT_L", "quart5", rep(-3, 5), rep(2, 5), rep(1, 5), maxeval = 800))
  gl_check("gn_direct_l_rand_poly3",
           gl_run("GN_DIRECT_L_RAND", "poly3", B3[[1]], B3[[2]], c(0, 0, 0), xtol_rel = 1e-8, maxeval = 600))
  # optimize.c runs DIRECT-L for all three NOSCAL variants (its which_alg
  # compares with GN_DIRECT and GN_DIRECT_L_RAND only)
  for (a in c("GN_DIRECT_NOSCAL", "GN_DIRECT_L_NOSCAL", "GN_DIRECT_L_RAND_NOSCAL"))
    gl_check(paste0(tolower(a), "_rosen"), gl_run(a, "rosen", B2[[1]], B2[[2]], c(-1.2, 1), xtol_rel = 1e-8, maxeval = 500))
  gl_check("gn_direct_magic", gl_run("GN_DIRECT", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), maxeval = 400,
                                     params = list(magic_eps = 1e-4)))
  gl_check("gn_direct_l_fixed", gl_run("GN_DIRECT_L", "poly3", c(-5, 2, -5), c(5, 2, 5), c(0, 2, 0), maxeval = 400))
  gl_check("gn_direct_l_ftol", gl_run("GN_DIRECT_L", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), ftol_rel = 1e-3, maxeval = 5000))
  expect_identical(gl_run("GN_DIRECT", "rosen", c(-Inf, -2), B2[[2]], c(-1.2, 1), maxeval = 100)$ret, "INVALID_ARGS")
  # the xtol stop of cdirect.c (the libnlopt 2.11.0 C library, built for
  # arm64 macOS: 29 evaluations)
  r <- gl_run("GN_DIRECT", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), xtol_rel = 0.3, maxeval = 3000)
  expect_identical(r$ret, "XTOL_REACHED")
  if (ref_platform()) {
    expect_identical(r$numevals, 29L)
    expect_identical(r$x, c(0x1.5555555555554p+0, 0x1.c71c71c71c71cp+0))
  }
})

test_that("the DIRECT of Gablonsky (GN_ORIG_DIRECT) follows NLopt", {
  gl_check("gn_orig_direct_rosen",
           gl_run("GN_ORIG_DIRECT", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), xtol_rel = 1e-8, maxeval = 500))
  gl_check("gn_orig_direct_l_rosen",
           gl_run("GN_ORIG_DIRECT_L", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), xtol_rel = 1e-8, maxeval = 500))
  gl_check("gn_orig_direct_xtol",
           gl_run("GN_ORIG_DIRECT", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), xtol_rel = 1e-3, maxeval = 20000))
  gl_check("gn_orig_direct_l_xtol",
           gl_run("GN_ORIG_DIRECT_L", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), xtol_rel = 1e-3, maxeval = 20000))
  gl_check("gn_orig_direct_magic", gl_run("GN_ORIG_DIRECT", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), maxeval = 400,
                                          params = list(magic_eps = 1e-4, magic_eps_abs = 1e-6)))
  # an infeasible center takes the value of the feasible centers around it
  # (dirreplaceinf); DIRsubrout.c reads the side levels transposed, which
  # upstream_replaceinf reproduces
  gl_check("gn_orig_direct_cons", gl_run("GN_ORIG_DIRECT", "poly3", B3[[1]], B3[[2]], c(0, 0, 0), cons = TRUE,
                                         maxeval = 600, upstream_replaceinf = TRUE))
  r <- gl_run("GN_ORIG_DIRECT", "poly3", B3[[1]], B3[[2]], c(0, 0, 0), cons = TRUE, maxeval = 600)
  expect_identical(r$ret, "MAXEVAL_REACHED")
  expect_lte(sum(r$x^2), 4)
  expect_lt(r$minf, 1)
})

test_that("PRAXIS follows NLopt", {
  I3 <- rep(Inf, 3)
  gl_check("ln_praxis_poly3", gl_run("LN_PRAXIS", "poly3", -I3, I3, c(0, 0, 0), xtol_rel = 1e-10, maxeval = 2000))
  gl_check("ln_praxis_poly3_box", gl_run("LN_PRAXIS", "poly3", B3[[1]], B3[[2]], c(0, 0, 0), xtol_rel = 1e-8, maxeval = 2000))
  gl_check("ln_praxis_poly3_half", gl_run("LN_PRAXIS", "poly3", rep(-1, 3), I3, c(0, 0, 0), xtol_rel = 1e-8, maxeval = 2000))
  gl_check("ln_praxis_poly3_upper",
           gl_run("LN_PRAXIS", "poly3", -I3, rep(1.5, 3), c(0, 0, 0), xtol_rel = 1e-8, maxeval = 2000))
  gl_check("ln_praxis_rosen", gl_run("LN_PRAXIS", "rosen", -I3[1:2], I3[1:2], c(-1.2, 1), ftol_rel = 1e-14, maxeval = 3000))
  gl_check("ln_praxis_rosen_box",
           gl_run("LN_PRAXIS", "rosen", c(-2, -2), c(0.8, 2), c(-1.2, 1), ftol_rel = 1e-14, maxeval = 3000))
  gl_check("ln_praxis_rosen_max", gl_run("LN_PRAXIS", "rosen", -I3[1:2], I3[1:2], c(-1.2, 1), maxeval = 25))
  gl_check("ln_praxis_quart5", gl_run("LN_PRAXIS", "quart5", rep(-3, 5), rep(3, 5), rep(1, 5), xtol_rel = 1e-9, maxeval = 4000))
  gl_check("ln_praxis_t0", gl_run("LN_PRAXIS", "rosen", -I3[1:2], I3[1:2], c(-1.2, 1), maxeval = 3000,
                                  params = list(t0_tol = 1e-6)))
})

test_that("StoGO follows NLopt", {
  gl_check("gd_stogo_rosen", gl_run("GD_STOGO", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), xtol_rel = 1e-8, maxeval = 600))
  gl_check("gd_stogo_rand_rosen",
           gl_run("GD_STOGO_RAND", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), xtol_rel = 1e-8, maxeval = 600))
  gl_check("gd_stogo_quart5", gl_run("GD_STOGO", "quart5", rep(-3, 5), rep(2, 5), rep(1, 5), maxeval = 800))
  expect_identical(gl_run("GD_STOGO", "rosen", c(-Inf, -2), B2[[2]], c(-1.2, 1), maxeval = 100)$ret, "INVALID_ARGS")
})

test_that("AGS follows NLopt", {
  gl_check("gn_ags_rosen", gl_run("GN_AGS", "rosen", B2[[1]], B2[[2]], c(-1.2, 1), maxeval = 500))
  # the constraint is evaluated first; maxeval counts its evaluations, the
  # number of evaluations those of the objective
  gl_check("gn_ags_cons", gl_run("GN_AGS", "poly3", B3[[1]], B3[[2]], c(0, 0, 0), cons = TRUE, maxeval = 600))
  # pow(x, 2) of the C library is not always x * x
  expect_identical(.nlopt_libm_pow(3, 2), 9)
  expect_identical(.nlopt_libm_pow(0, 2), 0)
})

test_that("the global NLopt algorithms run through solve()", {
  for (a in c("GN_DIRECT_L", "GN_ORIG_DIRECT", "GN_AGS", "GD_STOGO")) {
    NLopt$srand(1234)
    o <- Opt(a, 2)
    o$lower_bounds <- c(-2, -2); o$upper_bounds <- c(2, 2); o$maxeval <- 300L
    f <- OptimizationFunction(function(x, p) gl_rosen(x), grad = function(x, p) gl_roseng(x))
    s <- solve(OptimizationProblem(f, c(-1.2, 1), numeric()), o)
    expect_true(all(s$u >= -2 & s$u <= 2), info = a)
    expect_lt(s$objective, gl_rosen(c(-1.2, 1)))
  }
  s <- solve(OptimizationProblem(function(x, p) gl_rosen(x), c(-1.2, 1), numeric()), NLopt$LN_PRAXIS(), maxiters = 2000)
  expect_equal(unname(s$u), c(1, 1), tolerance = 1e-6)
})
