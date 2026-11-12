# Writes tests/testthat/test-implicit-julia.R from the output of
# tools/implicit-probe.jl:
#   REF=implicit.tsv OUT=tests/testthat/test-implicit-julia.R Rscript tools/implicit-gen-test.R
ref <- read.delim(Sys.getenv("REF"), header = FALSE, stringsAsFactors = FALSE)
out <- c(
"# Reference values from OrdinaryDiffEq (OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF)",
"# for an out of place nonlinear heat equation with an analytic Jacobian:",
"# steps, function evaluations, Jacobians, factorizations, linear solves,",
"# rejections and the final state, plain, with a tstop and with a discrete",
"# callback that halves the state at the tstop.",
"",
"heat_n <- 8",
"heat_f <- function(u, p, t) {",
"  l <- c(0, u[-heat_n]); r <- c(u[-1], 0)",
"  p[1] * (l - 2 * u + r) * (heat_n + 1)^2 - 0.5 * (u * u * u) + 0.1 * t",
"}",
"heat_jac <- function(u, p, t) {",
"  J <- matrix(0, heat_n, heat_n)",
"  for (i in seq_len(heat_n)) {",
"    J[i, i] <- -2 * p[1] * (heat_n + 1)^2 - 1.5 * u[i]^2",
"    if (i > 1) J[i, i - 1] <- p[1] * (heat_n + 1)^2",
"    if (i < heat_n) J[i, i + 1] <- p[1] * (heat_n + 1)^2",
"  }",
"  J",
"}",
"heat_cases <- list(",
paste0(vapply(seq_len(nrow(ref)), function(k) {
  u <- as.numeric(strsplit(ref[k, 9], ",")[[1]])
  sprintf("  list(alg = \"%s\", case = \"%s\", n = %dL, nf = %dL, njacs = %dL, nw = %dL, nsolve = %dL, nreject = %dL,\n       u = c(%s))",
          ref[k, 1], ref[k, 2], ref[k, 3], ref[k, 4], ref[k, 5], ref[k, 6], ref[k, 7], ref[k, 8],
          paste(sprintf("%.17g", u), collapse = ", "))
}, ""), collapse = ",\n"),
")",
"",
"test_that(\"ImplicitEuler, Trapezoid and FBDF follow OrdinaryDiffEq's out of place path\", {",
"  i <- seq_len(heat_n)",
"  prob <- ODEProblem(ODEFunction(heat_f, jac_fn = heat_jac), i * (heat_n + 1 - i) / heat_n^2, c(0, 0.3), 1)",
"  cb <- DiscreteCallback(function(u, t, integ) t == 0.1, function(integ) integ$u <- integ$u * 0.5)",
"  for (r in heat_cases) {",
"    alg <- get(r$alg)()",
"    sol <- switch(r$case,",
"      plain = solve(prob, alg, reltol = 1e-4, abstol = 1e-6),",
"      tstop = solve(prob, alg, reltol = 1e-4, abstol = 1e-6, tstops = 0.1),",
"      cb = solve(prob, alg, reltol = 1e-4, abstol = 1e-6, tstops = 0.1, callback = cb))",
"    lab <- paste(r$alg, r$case)",
"    st <- sol$stats",
"    expect_equal(length(sol$t), r$n, label = lab)",
"    expect_equal(c(st$nf, st$njacs, st$nw, st$nsolve, st$nreject),",
"                 c(r$nf, r$njacs, r$nw, r$nsolve, r$nreject), label = lab)",
"    expect_equal(sol$u[[length(sol$u)]], r$u, tolerance = 1e-10, ignore_attr = TRUE, label = lab)",
"  }",
"})")
writeLines(out, Sys.getenv("OUT"))
