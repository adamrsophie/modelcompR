# Reference values from OrdinaryDiffEq (OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF)
# for an out of place nonlinear heat equation with an analytic Jacobian:
# steps, function evaluations, Jacobians, factorizations, linear solves,
# rejections and the final state, plain, with a tstop and with a discrete
# callback that halves the state at the tstop.

heat_n <- 8
heat_f <- function(u, p, t) {
  l <- c(0, u[-heat_n]); r <- c(u[-1], 0)
  p[1] * (l - 2 * u + r) * (heat_n + 1)^2 - 0.5 * (u * u * u) + 0.1 * t
}
heat_jac <- function(u, p, t) {
  J <- matrix(0, heat_n, heat_n)
  for (i in seq_len(heat_n)) {
    J[i, i] <- -2 * p[1] * (heat_n + 1)^2 - 1.5 * u[i]^2
    if (i > 1) J[i, i - 1] <- p[1] * (heat_n + 1)^2
    if (i < heat_n) J[i, i + 1] <- p[1] * (heat_n + 1)^2
  }
  J
}
heat_cases <- list(
  list(alg = "ImplicitEuler", case = "plain", n = 448L, nf = 1994L, njacs = 665L, nw = 665L, nsolve = 1326L, nreject = 218L,
       u = c(0.0070821908403659725, 0.013124284747865972, 0.017529460077107317, 0.019851300019821237, 0.019851300019821237, 0.017529460077107317, 0.013124284747865972, 0.0070821908403659725)),
  list(alg = "ImplicitEuler", case = "tstop", n = 448L, nf = 1993L, njacs = 663L, nw = 663L, nsolve = 1326L, nreject = 217L,
       u = c(0.0070821267801815628, 0.01312416437915474, 0.017529297928030475, 0.019851115643890641, 0.019851115643890641, 0.017529297928030478, 0.01312416437915474, 0.0070821267801815628)),
  list(alg = "ImplicitEuler", case = "cb", n = 424L, nf = 1878L, njacs = 626L, nw = 626L, nsolve = 1248L, nreject = 204L,
       u = c(0.0040585782875324399, 0.0074417635377036686, 0.009873437551554164, 0.011145205400949631, 0.011145205400949631, 0.009873437551554164, 0.0074417635377036677, 0.004058578287532439)),
  list(alg = "Trapezoid", case = "plain", n = 47L, nf = 117L, njacs = 45L, nw = 45L, nsolve = 67L, nreject = 1L,
       u = c(0.0069802016519743843, 0.012932608682305596, 0.017271216998855685, 0.019557637851657862, 0.019557637851657862, 0.017271216998855685, 0.012932608682305601, 0.0069802016519743852)),
  list(alg = "Trapezoid", case = "tstop", n = 49L, nf = 109L, njacs = 46L, nw = 46L, nsolve = 58L, nreject = 0L,
       u = c(0.0069812907570279891, 0.012934655486473365, 0.01727397458365721, 0.019560773611721163, 0.019560773611721163, 0.01727397458365721, 0.012934655486473365, 0.0069812907570279883)),
  list(alg = "Trapezoid", case = "cb", n = 55L, nf = 163L, njacs = 57L, nw = 57L, nsolve = 100L, nreject = 6L,
       u = c(0.0040037394416484479, 0.0073387003337072512, 0.0097345810442137273, 0.010987303733258034, 0.010987303733258034, 0.0097345810442137273, 0.0073387003337072512, 0.0040037394416484479)),
  list(alg = "FBDF", case = "plain", n = 53L, nf = 129L, njacs = 23L, nw = 23L, nsolve = 73L, nreject = 1L,
       u = c(0.006994150602061662, 0.012958832353140604, 0.017306560602578562, 0.019597838448755595, 0.019597838448755595, 0.017306560602578545, 0.012958832353140588, 0.0069941506020616507)),
  list(alg = "FBDF", case = "tstop", n = 54L, nf = 137L, njacs = 27L, nw = 27L, nsolve = 79L, nreject = 2L,
       u = c(0.0069939297127621985, 0.012958199715466162, 0.017305381460584777, 0.019596258681694984, 0.019596258681695015, 0.017305381460584926, 0.012958199715466318, 0.0069939297127623269)),
  list(alg = "FBDF", case = "cb", n = 69L, nf = 168L, njacs = 34L, nw = 34L, nsolve = 94L, nreject = 3L,
       u = c(0.0040106448448673764, 0.0073516569615524937, 0.0097520241331879517, 0.011007139237021658, 0.011007139237021639, 0.0097520241331879257, 0.0073516569615525058, 0.0040106448448674146))
)

test_that("ImplicitEuler, Trapezoid and FBDF follow OrdinaryDiffEq's out of place path", {
  i <- seq_len(heat_n)
  prob <- ODEProblem(ODEFunction(heat_f, jac_fn = heat_jac), i * (heat_n + 1 - i) / heat_n^2, c(0, 0.3), 1)
  cb <- DiscreteCallback(function(u, t, integ) t == 0.1, function(integ) integ$u <- integ$u * 0.5)
  for (r in heat_cases) {
    alg <- get(r$alg)()
    sol <- switch(r$case,
      plain = solve(prob, alg, reltol = 1e-4, abstol = 1e-6),
      tstop = solve(prob, alg, reltol = 1e-4, abstol = 1e-6, tstops = 0.1),
      cb = solve(prob, alg, reltol = 1e-4, abstol = 1e-6, tstops = 0.1, callback = cb))
    lab <- paste(r$alg, r$case)
    st <- sol$stats
    expect_equal(length(sol$t), r$n, label = lab)
    expect_equal(c(st$nf, st$njacs, st$nw, st$nsolve, st$nreject),
                 c(r$nf, r$njacs, r$nw, r$nsolve, r$nreject), label = lab)
    expect_equal(sol$u[[length(sol$u)]], r$u, tolerance = 1e-10, ignore_attr = TRUE, label = lab)
  }
})
