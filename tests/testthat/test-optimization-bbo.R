# BBO_borg_moea: the Borg MOEA of BlackBoxOptim.jl through OptimizationBBO,
# ported to R. The epsilon box fitness, the archive and the genetic
# operators (fed with the random numbers BlackBoxOptim drew) match Julia;
# the numbers come from tools/bbo-probe.jl. The runs use R's random
# numbers and are checked as the OptimizationBBO tests check them.

bbo_fun <- function(name) get(name, envir = asNamespace("modelcompR"))

bbo_replay <- function(draws, expr) {
  rng <- bbo_fun(".bbo_rng")
  rng$queue <- draws
  on.exit(rng$queue <- NULL)
  r <- force(expr)
  expect_length(rng$queue, 0L)
  r
}

test_that("BBO epsilon box indexing and comparison match BlackBoxOptim", {
  # tools/bbo-probe.jl, IDX and HAT lines
  sch <- function(eps, minimizing = TRUE)
    bbo_fun(".bbo_scheme")(ParetoFitnessScheme(2, is_minimizing = minimizing), eps)
  idx <- function(u, s) {
    f <- bbo_fun(".bbo_indexed")(u, s)
    c(f$agg, f$index, f$dist)
  }
  expect_equal(idx(c(0.123, 4.56), sch(0.1)), c(4.683, 1, 45, 0.64257295305668849), tolerance = 1e-15)
  expect_equal(idx(c(-0.35, 2), sch(0.1)), c(1.65, -4, 20, 0.50000000000000044), tolerance = 1e-15)
  expect_equal(idx(c(0.3, 0.7), sch(0.1)), c(1, 3, 7, 9.9301366129890925e-16), tolerance = 1e-15)
  expect_equal(idx(c(NaN, 1.5), sch(0.1)), c(NaN, 9223372036854775807, 15, 0))
  expect_equal(idx(c(0.123, -4.56), sch(c(0.1, 0.5), FALSE)), c(-4.437, 2, -9, 0.77929455278476045),
               tolerance = 1e-15)
  expect_equal(idx(c(1e-17, 3), sch(c(0.25, 0.1))), c(3, 0, 30, 4.0000000000000003e-17), tolerance = 1e-15)
  s <- sch(0.1)
  hat <- function(u, v) bbo_fun(".bbo_hat")(bbo_fun(".bbo_indexed")(u, s), bbo_fun(".bbo_indexed")(v, s), TRUE)
  expect_equal(c(hat(c(0.12, 0.33), c(0.25, 0.41)), hat(c(0.12, 0.45), c(0.25, 0.41)),
                 hat(c(0.12, 0.33), c(0.13, 0.34)), hat(c(0.11, 0.31), c(0.19, 0.39)),
                 hat(c(0.15, 0.35), c(0.15, 0.35)), hat(c(0.5, 0.1), c(0.1, 0.5))),
               c(-1, -1, -1, -1, 0, 0))
  expect_error(sch(0), "must be positive")
  expect_error(sch(c(0.1, 0.1, 0.1)), "does not match the specified fitness dimensions (2)", fixed = TRUE)
})

test_that("the BBO epsilon box archive matches BlackBoxOptim", {
  # tools/bbo-probe.jl, ARCH, FRONT and TAGS lines
  s <- bbo_fun(".bbo_scheme")(ParetoFitnessScheme(2), 0.1)
  a <- bbo_fun(".bbo_archive")(s, 2L, 100L)
  state <- list()
  for (k in 1:60) {
    u <- c(abs(sin(1.7 * k)) * 2, abs(cos(0.9 * k)) * (1 + k / 30))
    bbo_fun(".bbo_add_candidate")(a, bbo_fun(".bbo_indexed")(u, s), c(k / 10, -k / 7), k %% 7 + 1L, k)
    if (k == 30) {
      a$n_restarts <- a$n_restarts + 1L
      a$last_restart <- a$num_candidates
    }
    state[[k]] <- c(length(a$FD), a$num_candidates, a$last_progress, bbo_fun(".bbo_noprogress")(a, TRUE),
                    a$best$x, a$best$fit$agg)
  }
  expect_equal(state[[10]], c(1, 10, 2, 8, 0.2, -0.2857142857142857, 0.75343110505962207), tolerance = 1e-15)
  expect_equal(state[[26]], c(4, 26, 26, 0, 2.6, -3.7142857142857144, 0.73295206578158023), tolerance = 1e-15)
  expect_equal(state[[60]], c(4, 60, 37, 23, 2.6, -3.7142857142857144, 0.73295206578158023), tolerance = 1e-15)
  expect_equal(cbind(a$FO, a$tag, a$nfe, a$nrs, a$X),
               rbind(c(0.51108220405366245, 0.24234890100595957, 3, 2, 0, 0.2, -0.2857142857142857),
                     c(0.081386514699729712, 1.664049194845419, 4, 24, 0, 2.4, -3.4285714285714284),
                     c(0.43197453237643429, 0.30097753340514599, 6, 26, 0, 2.6, -3.7142857142857144),
                     c(0.13618838936899627, 0.6882640227204877, 3, 37, 1, 3.7, -5.2857142857142856)),
               tolerance = 1e-15, ignore_attr = TRUE)
  expect_equal(bbo_fun(".bbo_tagcounts")(a, 0.9, 7L), c(0, 0, 1.9, 0.9, 0, 0.9, 0), tolerance = 1e-15)
  expect_error(bbo_fun(".bbo_tagcounts")(a, 1.5, 7L), "should be in (0.0, 1.0] range", fixed = TRUE)
})

test_that("the BBO genetic operators match BlackBoxOptim on its random numbers", {
  # tools/bbo-probe.jl, OP lines: the operator output and the numbers it drew
  pop <- matrix(c(0.1, 0.7, -0.3, 0.45, 0.9, -0.8,
                  0.5, -0.2, 0.6, 0.05, -0.4, 0.3,
                  -0.6, 0.2, 0.8, -0.1, 0.35, 0.0,
                  0.25, 0.9, -0.7, 0.6, -0.15, 0.55), 4, byrow = TRUE)
  lo <- rep(-1, 4); hi <- rep(1, 4)
  d3 <- c(0.33973478919014277, 0.95923289796888112, 0.99997948900346512, 0.092114272085112137, 0.69842599816368145)
  n3 <- c(-0.66376355155512712, 1.3632084599735597, 0.71741058992417861)
  expect_equal(bbo_replay(c(2, d3[2:4]), bbo_fun(".bbo_de")(0.7, 0.6)(pop, c(2, 5, 3))[[1]]),
               c(-0.3, 0.72, 0.8, -0.069999999999999951), tolerance = ref_tol(1e-14, 1e-8))
  expect_equal(bbo_replay(c(0.073366354469292849, 0.34924148955718615, 0, 0.62826474034250168, 0.91492900366283136,
                            0.19280811624587546, 0.77018034788566636, 1),
                          unlist(bbo_fun(".bbo_sbx")(0.2, 3)(pop, c(4, 1)))),
               c(0.43498416743976287, 0.05, -0.1, 0.21246338519109928, 0.1150158325602372, 0.5, -0.6,
                 0.63753661480890067), tolerance = ref_tol(1e-14, 1e-8))
  expect_equal(bbo_replay(c(d3, 1), unlist(bbo_fun(".bbo_sbx")(0.2, 3)(pop, c(4, 1)))),
               c(0.45, 0.05, -0.1, 0.22642156897020563, 0.1, 0.5, -0.6, 0.62357843102979438), tolerance = ref_tol(1e-14, 1e-8))
  expect_equal(bbo_replay(d3[1:2], bbo_fun(".bbo_spx")(3L, 1.1)(pop, c(6, 2, 4))[[1]]),
               c(0.20365227195462099, -0.036330445314106188, 0.13666780372072151, 0.78676607771679885),
               tolerance = ref_tol(1e-14, 1e-8))
  # PCX with two parents is degenerate: the other parent lies on the line of the first
  expect_equal(bbo_replay(n3[1], bbo_fun(".bbo_pcx")(2L, 0.1, 0.5)(pop, c(3, 5))[[1]]),
               c(-0.10087093453346183, 0.43405911211121817, 0.72532660045004826, -0.60873251166116993),
               tolerance = ref_tol(1e-14, 1e-8))
  expect_equal(bbo_replay(n3, bbo_fun(".bbo_pcx")(3L, 0.1, 0.5)(pop, c(1, 4, 6))[[1]]),
               c(0.18631343350166968, 0.5700114587948838, -0.7281702031576992, 0.17628381921797956),
               tolerance = ref_tol(1e-14, 1e-8))
  # UNDX draws along the singular vectors of LAPACK, whose signs may differ
  # between LAPACK builds
  if (ref_platform()) {
    expect_equal(bbo_replay(n3, bbo_fun(".bbo_undx")(3L, 0.5, 0.1)(pop, c(5, 3, 1))[[1]]),
                 c(0.50049482239139476, -0.10152369799888181, 0.24946146563453928, -0.17346921085021264),
                 tolerance = 1e-14)
  }
  clock <- function(pm) {
    mc <- bbo_fun(".bbo_clock")(lo, hi, pm, 1 / 4)
    c(mc$apply(mc$apply(pop[, 2])), mc$nextvar)
  }
  expect_equal(bbo_replay(d3, clock(TRUE)), c(0.7, -0.14244412411892898, 0.16085099985015747, 0.9, 1),
               tolerance = ref_tol(1e-14, 1e-8))
  expect_equal(bbo_replay(d3, clock(FALSE)), c(0.7, 0.91846579593776223, -0.81577145582977573, 0.9, 1),
               tolerance = ref_tol(1e-14, 1e-8))
  expect_equal(bbo_replay(d3[1:2], bbo_fun(".bbo_bound")(c(1.3, -0.2, -1.7, 0.95), pop[, 2], lo, hi)),
               c(0.89807956324295712, -0.2, 0.15107947756265738, 0.95), tolerance = ref_tol(1e-14, 1e-8))
})

test_that("solve with BBO_borg_moea follows the OptimizationBBO tests", {
  # the multiobjective tests of OptimizationBBO (there with the default
  # 10000 steps; tools/bbo-probe.jl, RUN lines, gives the Julia runs)
  u0 <- c(0.25, 0.25)
  mo <- function(f) OptimizationProblem(MultiObjectiveOptimizationFunction(f), u0, lb = c(0, 0), ub = c(2, 2))
  sphere_rastrigin <- function(x, p) c(sum(x^2), sum(x^2 - 10 * cos(2 * pi * x) + 10))
  set.seed(1)
  sol <- solve(mo(sphere_rastrigin), BBO_borg_moea(), NumDimensions = 2,
               FitnessScheme = ParetoFitnessScheme(2, is_minimizing = TRUE), maxiters = 2000)
  expect_lt(abs(sol$objective[1] - 6.9905986e-18), 1e-3)
  expect_lt(abs(sol$objective[2] - 1.7763568e-15), 1e-3)
  expect_equal(sol$retcode, "MaxIters")
  expect_identical(sol$stats$iterations, 2001L)
  expect_identical(sol$original$stop_reason, "Max number of steps (2000) reached")
  expect_gte(sol$original$state$n_restarts, 1L)
  rosen_ackley <- function(x, p)
    c((1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2,
      -20 * exp(-0.2 * sqrt(0.5 * (x[1]^2 + x[2]^2))) - exp(0.5 * (cos(2 * pi * x[1]) + cos(2 * pi * x[2]))) +
        exp(1) + 20)
  set.seed(1)
  sol <- solve(mo(rosen_ackley), BBO_borg_moea(), maxiters = 2000)
  expect_lt(abs(sol$objective[1] - 0.97438), 1e-3)
  expect_lt(abs(sol$objective[2] - 0.04088), 1e-3)
  # the archive: points not dominated in boxes of side 0.1, inside the bounds
  box <- floor(sol$original$frontier$F / 0.1 + 10 * .Machine$double.eps)
  dominated <- outer(seq_len(nrow(box)), seq_len(nrow(box)), Vectorize(function(i, j)
    all(box[i, ] <= box[j, ]) && any(box[i, ] < box[j, ])))
  expect_false(any(dominated))
  expect_true(all(sol$original$frontier$X >= 0 & sol$original$frontier$X <= 2))
  zdt1 <- function(x, p) {
    f1 <- x[1]
    g <- 1 + 9 * sum(x[2:length(x)]) / (length(x) - 1)
    c(f1, g * (1 - sqrt(f1 / g)))
  }
  set.seed(1)
  sol <- solve(mo(zdt1), BBO_borg_moea(), maxiters = 2000)
  expect_identical(sol$objective, zdt1(sol$u, NULL))
  expect_lt(abs(sol$objective[2] - (1 - sqrt(sol$objective[1]))), 1e-3)
})

test_that("BBO_borg_moea runs the callback after every step", {
  sphere_rastrigin <- function(x, p) c(sum(x^2), sum(x^2 - 10 * cos(2 * pi * x) + 10))
  prob <- OptimizationProblem(MultiObjectiveOptimizationFunction(sphere_rastrigin), c(0.25, 0.25),
                              lb = c(0, 0), ub = c(2, 2))
  hist <- list()
  cb <- function(state, objective) {
    hist[[length(hist) + 1L]] <<- list(state = state, objective = objective)
    FALSE
  }
  set.seed(2)
  sol <- solve(prob, BBO_borg_moea(), callback = cb, maxiters = 50)
  expect_length(hist, 51L)
  first <- hist[[1L]]
  expect_identical(first$state$objective, first$objective)
  expect_length(first$objective, 2L)
  expect_identical(first$state$iter, 1L)
  expect_identical(hist[[51L]]$objective, sol$objective)
  # TRUE stops the run (tools/bbo-probe.jl, SHUTDOWN: after 25 steps)
  n <- 0L
  set.seed(2)
  expect_warning(sol <- solve(prob, BBO_borg_moea(), callback = function(state, objective) {
    n <<- n + 1L
    n == 25L
  }), "Unrecognized stop reason: Run explicitly stopped via shutdown method. Defaulting to ReturnCode.Default.",
  fixed = TRUE)
  expect_identical(sol$stats$iterations, 25L)
  expect_equal(sol$retcode, "Default")
  expect_error(solve(prob, BBO_borg_moea(), callback = function(state, objective) 1, maxiters = 5),
               "The callback should return a boolean `halt` for whether to stop the optimization process.",
               fixed = TRUE)
})

test_that("BBO_borg_moea checks its arguments as OptimizationBBO and BlackBoxOptim do", {
  f <- function(x, p) c(sum(x^2), sum((x - 1)^2))
  prob <- OptimizationProblem(MultiObjectiveOptimizationFunction(f), c(0.25, 0.25), lb = c(0, 0), ub = c(2, 2))
  set.seed(3)
  # tools/bbo-probe.jl, SHORT: MaxSteps = 100 stops after 101 steps
  sol <- solve(prob, BBO_borg_moea(), maxiters = 100)
  expect_identical(sol$stats$iterations, 101L)
  expect_gt(sol$stats$fevals, 150L)
  expect_identical(sol$stats$fevals, sol$original$f_calls)
  expect_output(print(sol), paste0("Objective: ", paste(format(sol$objective), collapse = " ")), fixed = TRUE)
  sol <- solve(prob, BBO_borg_moea(), maxtime = 0.2)
  expect_equal(sol$retcode, "MaxTime")
  expect_identical(sol$original$stop_reason, "Max time (0.2 s) reached")
  expect_message(solve(prob, BBO_borg_moea(), maxiters = 5, reltol = 1e-3),
                 "common reltol is currently not used by BBO_borg_moea(:borg_moea)", fixed = TRUE)
  expect_error(solve(prob, BBO_borg_moea(), maxiters = 0),
               "The number of maxiters has to be a non-negative and non-zero number.", fixed = TRUE)
  expect_error(solve(prob, BBO_borg_moea(), maxtime = -1),
               "The maximum time has to be a non-negative and non-zero number.", fixed = TRUE)
  expect_error(solve(prob, BBO_borg_moea(), maxiters = 5, FitnessScheme = ParetoFitnessScheme(3)),
               "does NOT return the expected fitness type", fixed = TRUE)
  expect_error(solve(prob, BBO_borg_moea(), maxiters = 5, PopulationSize = 1),
               "The population size MUST be at least 2", fixed = TRUE)
  far <- OptimizationProblem(MultiObjectiveOptimizationFunction(f), c(3, 0.25), lb = c(0, 0), ub = c(2, 2))
  expect_error(solve(far, BBO_borg_moea(), maxiters = 5),
               "Starting point [3.0, 0.25] is not in the search space/range.", fixed = TRUE)
  expect_error(solve(OptimizationProblem(MultiObjectiveOptimizationFunction(f), c(0.25, 0.25)), BBO_borg_moea()),
               class = "IncompatibleOptimizerError")
  expect_true(requiresbounds(BBO_borg_moea()))
  expect_true(allowscallback(BBO_borg_moea()))
})
