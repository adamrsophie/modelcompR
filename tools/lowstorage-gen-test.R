# Writes tests/testthat/test-lowstorage.R from the runs of
# tools/lowstorage-probe.jl:
#
#   julia --project=tools/julia tools/lowstorage-probe.jl - runs.tsv
#   Rscript tools/lowstorage-gen-test.R runs.tsv tests/testthat/test-lowstorage.R
args <- commandArgs(trailingOnly = TRUE)
runs <- read.delim(args[1], header = FALSE, stringsAsFactors = FALSE, quote = "")
names(runs) <- c("alg", "case", "n", "nf", "nacc", "nrej", "u", "t")
hxs <- function(s) strsplit(s, ",", fixed = TRUE)[[1]]
rp <- c("CKLLSRK43_2", "CKLLSRK54_3C", "CKLLSRK95_4C", "CKLLSRK95_4M", "CKLLSRK95_4S", "CKLLSRK54_3C_3R",
        "CKLLSRK54_3M_3R", "CKLLSRK54_3N_3R", "CKLLSRK85_4C_3R", "CKLLSRK85_4M_3R", "CKLLSRK85_4P_3R",
        "CKLLSRK54_3M_4R", "CKLLSRK54_3N_4R", "CKLLSRK65_4M_4R", "CKLLSRK85_4FM_4R", "CKLLSRK75_4M_5R")
lines <- character()
for (i in seq_len(nrow(runs))) {
  r <- runs[i, ]
  # the R+ methods differ on purpose once a step is rejected
  if (r$alg %in% rp && r$nrej > 0) next
  ts <- hxs(r$t)
  pick <- unique(c(length(ts) %/% 2L, length(ts) - 1L))
  pick <- pick[pick >= 1L]
  lines <- c(lines, sprintf('  "%s", "%s", %d, %d, %d, %d, "%s", "%s", "%s",', r$alg, r$case, r$n, r$nf, r$nacc,
                            r$nrej, r$u, paste(pick, collapse = ","), paste(ts[pick], collapse = ",")))
}
lines[length(lines)] <- sub(",$", "", lines[length(lines)])
hdr <- c(
  "# The low storage Runge-Kutta methods (R/solvers-lowstorage.R) against",
  "# OrdinaryDiffEqLowStorageRK: tools/lowstorage-probe.jl runs every method with",
  "# a fixed step, a fixed step with saveat and, for the adaptive ones, with",
  "# tolerances, with step rejections and with saveat. Steps, function",
  "# evaluations, accepted and rejected steps, the final state and two of the",
  "# saved times agree to the bit on the reference platform (helper-platform.R;",
  "# elsewhere the final state is compared to a tolerance). The runs of the R+",
  "# methods with rejected steps differ on purpose (see PORTING.md) and are",
  "# checked against a reference solution.",
  "",
  "lsrk_f <- function(u, p, t) c(u[2], -u[1] + 0.3 * (1 - t / 3) * u[1]^2)",
  "lsrk_g <- function(u, p, t) c(-20 * (u[1] - t * (2 - t)), 0.5 * u[1] * u[2])",
  "lsrk_P1 <- ODEProblem(lsrk_f, c(1, 0.5), c(0, 3))",
  "lsrk_P2 <- ODEProblem(lsrk_g, c(0, 1), c(0, 2))",
  "lsrk_hex <- function(s) as.numeric(strsplit(s, \",\", fixed = TRUE)[[1]])",
  "",
  "lsrk_solve <- function(alg, case) {",
  "  switch(case,",
  "    fixed = solve(lsrk_P1, alg, dt = 0.05, adaptive = FALSE),",
  "    fixed_saveat = solve(lsrk_P1, alg, dt = 0.05, adaptive = FALSE, saveat = 0.37),",
  "    adaptive = solve(lsrk_P1, alg, reltol = 1e-6, abstol = 1e-8),",
  "    reject = solve(lsrk_P2, alg, reltol = 1e-4, abstol = 1e-6, dt = 0.5),",
  "    saveat = solve(lsrk_P1, alg, reltol = 1e-6, abstol = 1e-8, saveat = 0.37))",
  "}",
  "",
  "lsrk_ref <- matrix(c(")
ftr <- c(
  "), ncol = 9, byrow = TRUE)",
  "",
  "test_that(\"the low storage methods reproduce OrdinaryDiffEqLowStorageRK\", {",
  "  for (i in seq_len(nrow(lsrk_ref))) {",
  "    r <- lsrk_ref[i, ]",
  "    lab <- paste(r[1], r[2])",
  "    s <- lsrk_solve(get(r[1])(), r[2])",
  "    if (!ref_platform()) {",
  "      expect_equal(s$u[[length(s$u)]], lsrk_hex(r[7]), tolerance = 1e-6, label = paste(lab, \"u\"))",
  "      next",
  "    }",
  "    expect_identical(length(s$t), as.integer(r[3]), label = paste(lab, \"steps\"))",
  "    expect_identical(s$stats$nf, as.integer(r[4]), label = paste(lab, \"nf\"))",
  "    expect_identical(s$stats$naccept, as.integer(r[5]), label = paste(lab, \"naccept\"))",
  "    expect_identical(s$stats$nreject, as.integer(r[6]), label = paste(lab, \"nreject\"))",
  "    expect_identical(s$u[[length(s$u)]], lsrk_hex(r[7]), label = paste(lab, \"u\"))",
  "    expect_identical(s$t[as.integer(strsplit(r[8], \",\")[[1]])], lsrk_hex(r[9]), label = paste(lab, \"t\"))",
  "  }",
  "})",
  "",
  "test_that(\"the R+ methods start a retried step from the accepted state\", {",
  "  ref <- solve(lsrk_P2, Vern9(), reltol = 1e-12, abstol = 1e-14)",
  "  uref <- ref$u[[length(ref$u)]]",
  sprintf("  for (a in c(%s)) {", paste0("\"", rp, "\"", collapse = ", ")),
  "    s <- lsrk_solve(get(a)(), \"reject\")",
  "    expect_true(s$stats$nreject > 0L, label = a)",
  "    expect_equal(s$u[[length(s$u)]], uref, tolerance = 1e-3, label = a)",
  "  }",
  "})",
  "",
  "test_that(\"the low storage methods interpolate and restart after events\", {",
  "  s <- solve(lsrk_P1, RDPK3Sp35(), reltol = 1e-6, abstol = 1e-8)",
  "  expect_equal(s(1.234), solve(lsrk_P1, Vern9(), reltol = 1e-12, abstol = 1e-14)(1.234), tolerance = 1e-5)",
  "  cb <- DiscreteCallback(function(u, t, integ) t == 1, function(integ) integ$u <- integ$u * 0.5)",
  "  for (a in c(\"SHLDDRK_2N\", \"CarpenterKennedy2N54\", \"RDPK3SpFSAL49\")) {",
  "    s <- solve(lsrk_P1, get(a)(), dt = 0.05, adaptive = FALSE, callback = cb, tstops = 1)",
  "    expect_equal(s(3), solve(lsrk_P1, Vern9(), callback = cb, tstops = 1, reltol = 1e-12, abstol = 1e-14)(3),",
  "                 tolerance = 1e-5, label = a)",
  "  }",
  "})")
writeLines(c(hdr, lines, ftr), args[2])
