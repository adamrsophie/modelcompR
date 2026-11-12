# OptimizationMetaheuristics: the algorithms of Metaheuristics.jl ported to
# R. The deterministic parts are checked against tools/metaheuristics-probe.jl
# (Metaheuristics.jl 3.6). The algorithms draw R's random numbers, so they are
# checked by convergence with fixed seeds, as the tests of
# OptimizationMetaheuristics do with Julia's.

mh <- function(name) get(name, envir = asNamespace("modelcompR"))

test_that("the operators match Metaheuristics.jl with the same random numbers", {
  # values from tools/metaheuristics-probe.jl
  sbx <- mh(".mh_sbx_core")
  v1 <- c(0.1, 0.5, -0.3, 0.9); v2 <- c(0.4, 0.2, 0.8, 0.9); lb <- rep(-1, 4); ub <- rep(1, 4)
  cc <- sbx(v1, v2, lb, ub, 15, 1,
            c(0.17842277698137288, 0.6312854615169121, 0.9254079413819496, 0.7387916162980288),
            c(0.20047525994825222, 0.010061581092753569, 0.30512028596709184, 0.7225440433669819),
            c(FALSE, TRUE, FALSE, TRUE))
  expect_equal(cc[[1]], c(0.10832786692771235, 0.46750930922061107, -0.2832814536088635, 0.9), tolerance = ref_tol(1e-15, 1e-8))
  expect_equal(cc[[2]], c(0.3916721330722587, 0.2324906907791513, 0.7831647134982135, 0.9), tolerance = ref_tol(1e-15, 1e-8))
  cc <- sbx(v1, v2, lb, ub, 2, 0.9,
            c(0.12396619710885193, 0.689736023123346, 0.3420498402131924, 0.895745573003232),
            c(0.5709455076043957, 0.7993168190264921, 0.5566493338214822, 0.17749190355115108),
            c(TRUE, FALSE, TRUE, TRUE))
  expect_equal(cc[[1]], c(0.407570626450684, 0.14683479177053124, 0.7797915721593236, 0.9), tolerance = ref_tol(1e-15, 1e-8))
  expect_equal(cc[[2]], c(0.0922107430306705, 0.5517176505287152, -0.31263744054308407, 0.9), tolerance = ref_tol(1e-15, 1e-8))
  pm <- mh(".mh_polynomial_core")(v1, lb, ub, 20, 0.7,
                                   c(0.3389033221483959, 0.5756551729622779, 0.9223751031077165, 0.40367251644880686),
                                   c(0.5496834593293771, 0.1597213479747308, 0.048660092210817574))
  expect_equal(pm, c(0.10994253286130609, 0.3942166972752843, -0.3, 0.6899841380097201), tolerance = ref_tol(1e-15, 1e-8))
  vs <- mh(".mh_violations_sum")
  expect_equal(vs(c(1, -2, 0.5), c(1e-5, -0.3), 1e-4), 0.65)
  expect_equal(vs(numeric(0), c(0.2, -0.1)), 0.15000000000000002)
  expect_equal(vs(0.3, numeric(0)), 0.3)
  expect_equal(Metaheuristics$gen_ref_dirs(2, 5), list(c(0, 1), c(0.2, 0.8), c(0.4, 0.6), c(0.6, 0.4), c(0.8, 0.2), c(1, 0)))
  rd <- Metaheuristics$gen_ref_dirs(3, 4)
  expect_length(rd, 15)
  expect_equal(rd[[7]], c(0.25, 0.25, 0.5))
})

test_that("the helpers of ECA, CGSA, SA, RDEx, MOEAD_DE and MCCGA match Metaheuristics.jl", {
  # values from tools/metaheuristics-probe.jl
  expect_equal(mh(".mh_fitness_to_mass")(c(3, -1, 2.5, 0)), c(5, 9, 5.5, 8))
  expect_equal(mh(".mh_fitness_to_mass")(c(0, 0)), c(1, 1))
  U <- list(list(x = c(0.1, 0.2), f = 1), list(x = c(-0.5, 0.3), f = 2), list(x = c(0.7, -0.9), f = 0.5))
  expect_equal(mh(".mh_center")(U, c(0.3, 1.2, 0.5)), c(-0.10999999999999999, -0.015000000000000013), tolerance = ref_tol(1e-15, 1e-8))
  expect_equal(mh(".mh_mass_calculation")(c(3, 1, 2, 5)), c(0.2222222222222222, 0.4444444444444444, 0.3333333333333333, 0))
  expect_equal(mh(".mh_gconstant")(5, 50), 13.53352832366127, tolerance = ref_tol(1e-15, 1e-8))
  # the chaotic maps amplify last bit differences of the elementary functions
  chaos <- vapply(1:10, function(i) mh(".mh_chaos")(i, 7, 50, 3), 0)
  if (ref_platform()) {
    expect_equal(chaos,
                 c(2.4119553229135446, 2.5913697600578995, 2.9999999999999867, 0.3430733698655499, 0.9625317287962338,
                   2.9296874999999636, 1.7728631193203337, 1.934119068028505, 2.6186447346416886, 0.4026100236012814),
                 tolerance = 1e-15)
  } else {
    expect_true(all(is.finite(chaos)))
  }
  expect_equal(mh(".mh_sa_newsol")(c(0.3, -0.7, 0), 10^2.5), c(0.014637368041430806, -0.17505911580742006, 0),
               tolerance = ref_tol(1e-15, 1e-8))
  expect_equal(mh(".mh_rdex_weighted_mean")(c(0.5, 0.9, 0.2, 0.7), c(1, 3, 0.5, 2), 3), 0.8181818181818183,
               tolerance = ref_tol(1e-15, 1e-8))
  expect_equal(mh(".mh_moead_g")(c(0.3, 0.7), c(0.4, 0.6), c(0.1, 0)), 0.42)
  expect_equal(mh(".mh_moead_neighbors")(Metaheuristics$gen_ref_dirs(2, 7), 2L),
               list(c(2L, 3L), c(3L, 1L), c(2L, 4L), c(5L, 3L), c(4L, 6L), c(7L, 5L), c(6L, 8L), c(7L, 6L)))
  expect_equal(mh(".mh_mccga_bits")(0.1),
               c(0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1))
  expect_equal(mh(".mh_mccga_floats")(mh(".mh_mccga_bits")(c(0.1, -2.5, 3.75e-3))),
               c(0.10000000149011612, -2.5, 0.0037499999161809683), tolerance = ref_tol(1e-15, 1e-8))
  h <- mh(".mh_mccga_hj")(function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2, c(-1.2, 1))
  expect_equal(h$iter, 24L)
  expect_equal(h$step, 7.62939453125e-5)
  expect_equal(h$par, c(0.99970703125, 0.9993896484375))
})

mkpop <- function(F, cons = NULL)
  lapply(seq_along(F), function(i) mh(".mh_child")(c(i, 0), list(F[[i]], if (is.null(cons)) 0 else cons[i], 0), mo = TRUE))
F2 <- list(c(0.1, 0.9), c(0.4, 0.5), c(0.2, 0.8), c(0.6, 0.3), c(0.5, 0.6), c(0.9, 0.05), c(0.3, 0.95), c(0.8, 0.2))
F3 <- list(c(0.1, 0.9, 0.5), c(0.4, 0.5, 0.6), c(0.2, 0.8, 0.1), c(0.6, 0.3, 0.4), c(0.5, 0.6, 0.2),
           c(0.9, 0.05, 0.7), c(0.3, 0.7, 0.3), c(0.8, 0.2, 0.9))

test_that("sorting, crowding and the hypervolume contributions match Metaheuristics.jl", {
  # values from tools/metaheuristics-probe.jl
  pop <- mkpop(F2); pop3 <- mkpop(F3)
  expect_equal(Metaheuristics$non_dominated_sort(pop), c(1, 1, 1, 1, 2, 1, 2, 1))
  expect_equal(Metaheuristics$non_dominated_sort(pop3), rep(1, 8))
  popc <- mkpop(F2, ifelse(seq_along(F2) %% 3 == 0, 0.5, -1))
  expect_equal(Metaheuristics$non_dominated_sort(popc), c(1, 1, 3, 1, 2, 3, 2, 1))
  cr <- mh(".mh_crowding")(pop[c(1, 2, 3, 4, 6, 8)])
  expect_equal(vapply(cr$pop, function(s) s$x[1], 0), c(1, 3, 2, 4, 8, 6))
  expect_equal(cr$cd, c(Inf, 0.8455882352941178, 1.088235294117647, 0.8529411764705883, 0.6691176470588236, Inf),
               tolerance = ref_tol(1e-15, 1e-8))
  # three objectives: Metaheuristics.jl adds the distances by position while
  # it re-sorts the front; the scratch variant below reproduces its numbers,
  # the port gives every solution its own sum
  bug <- function(pop) {
    n <- length(pop); Fm <- mh(".mh_fvals")(pop); o <- seq_len(n); cd <- numeric(n)
    for (j in seq_len(ncol(Fm))) {
      o <- o[order(Fm[o, j])]
      cd[c(1L, n)] <- Inf
      for (i in 2:(n - 1L)) cd[i] <- cd[i] + (Fm[o[i + 1L], j] - Fm[o[i - 1L], j]) / (Fm[o[n], j] - Fm[o[1L], j])
    }
    cd
  }
  expect_equal(bug(pop3), c(Inf, 0.7941176470588235, 0.8529411764705883, 0.8529411764705883, 0.7352941176470587,
                            0.8602941176470589, 0.9852941176470591, Inf), tolerance = ref_tol(1e-15, 1e-8))
  cr <- mh(".mh_crowding")(pop3)
  expect_equal(vapply(cr$pop, function(s) s$x[1], 0), c(3, 5, 7, 4, 1, 2, 6, 8))
  Fs <- do.call(rbind, F3)
  own <- numeric(8)
  for (j in 1:3) {
    o <- order(Fs[, j])
    own[o[c(1, 8)]] <- Inf
    own[o[2:7]] <- own[o[2:7]] + (Fs[o[3:8], j] - Fs[o[1:6], j]) / (Fs[o[8], j] - Fs[o[1], j])
  }
  expect_equal(cr$cd, own[c(3, 5, 7, 4, 1, 2, 6, 8)])
  P2 <- rbind(c(0.1, 0.2, 0.4, 0.6, 0.8, 0.9), c(0.9, 0.8, 0.5, 0.3, 0.2, 0.05))
  expect_equal(mh(".mh_hv_contrib_2d")(P2), c(Inf, 0.019999999999999997, 0.06, 0.040000000000000015, 0.009999999999999995, Inf),
               tolerance = ref_tol(1e-15, 1e-8))
  expect_equal(mh(".mh_compute_contribution")(mkpop(F3[c(1, 2, 3, 4, 6, 7, 8)])),
               c(0.00408, 0.015200000000000005, 0.037360000000000004, 0.05320000000000001, 0.004959999999999996,
                 0.016600000000000004, 0.0007999999999999952), tolerance = ref_tol(1e-15, 1e-8))
})

test_that("the SPEA2 and NSGA3 steps match Metaheuristics.jl", {
  # values from tools/metaheuristics-probe.jl
  pop <- mkpop(F2); pop3 <- mkpop(F3)
  expect_equal(mh(".mh_spea2_fitness")(pop),
               c(0.45327725046598033, 0.4380503284503523, 0.4586571987969251, 0.4380503284503523, 1.4317364702520927,
                 0.41832034209348085, 2.4532772504659803, 0.44971979803797996), tolerance = ref_tol(1e-15, 1e-8))
  # three objectives: up to 4 ulp, from the dot products of Distances.jl
  expect_equal(mh(".mh_spea2_fitness")(pop3),
               c(0.41249631709111817, 0.4212005858348719, 0.4212005858348719, 0.4212005858348719, 0.4212005858348719,
                 0.40121280370198137, 0.44544442277200064, 0.3925074168904956), tolerance = ref_tol(1e-15, 1e-8))
  expect_equal(mh(".mh_spea2_truncation")(3, mh(".mh_pairwise_distances")(pop)),
               c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))
  pn <- mh(".mh_fast_nds")(pop3)
  Fn <- mh(".mh_nsga3_normalize")(pn)
  expect_equal(Fn[, 1], c(0, 0.9834710743801653, 0.06611570247933875), tolerance = ref_tol(1e-15, 1e-8))
  expect_equal(Fn[, 6], c(1.0578512396694215, 0, 0.09917355371900811), tolerance = ref_tol(1e-15, 1e-8))
  expect_equal(Fn[3, ], c(0.06611570247933875, 0.08264462809917343, 0, 0.04958677685950407, 0.016528925619834687,
                          0.09917355371900811, 0.033057851239669374, 0.1322314049586775), tolerance = ref_tol(1e-15, 1e-8))
  a <- mh(".mh_nsga3_associate")(Fn, Metaheuristics$gen_ref_dirs(3, 3), 4L)
  expect_equal(a$nich, c(4L, 7L, 4L, 9L, 7L, 10L, 7L, 10L))
  expect_equal(a$freq, c(0L, 0L, 0L, 2L, 0L, 0L, 1L, 0L, 0L, 0L))
  expect_equal(a$dist, c(0.00437128611433644, 0.021706167611501938, 0.01748514445734581, 0.0038248753500443886,
                         0.035803565330237014, 0.009835393757256988, 0.011051157707806835, 0.04760603783894539),
               tolerance = ref_tol(1e-15, 1e-8))
  expect_equal(mh(".mh_generational_distance")(rbind(c(0.1, 0.9), c(0.5, 0.5), c(0.95, 0.1)),
                                               rbind(c(0, 1), c(0.25, 0.75), c(0.5, 0.5), c(0.75, 0.25), c(1, 0)),
                                               inverted = TRUE),
               0.14307135789365263, tolerance = ref_tol(1e-15, 1e-8))
})

rosen <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
rprob <- OptimizationProblem(rosen, c(0, 0), c(1, 100), lb = c(-1, -1), ub = c(1, 1))
l1 <- rosen(c(0, 0), c(1, 100))

test_that("the single objective algorithms minimize the Rosenbrock function", {
  # seeded convergence tests, as the OptimizationMetaheuristics tests
  cases <- list(ECA = list(ECA(), 60), DE = list(DE(), 60), PSO = list(PSO(), 60), ABC = list(ABC(N = 20, options = Metaheuristics$Options(f_calls_limit = 3000)), 60),
                CGSA = list(CGSA(N = 20), 40), SA = list(SA(N = 100), 20), WOA = list(WOA(N = 20), 60),
                SHADE = list(SHADE(N = 20), 60), epsilonDE = list(epsilonDE(), 60), CSO = list(CSO(N = 20), 80),
                GA = list(GA(N = 30), 60), RDEx = list(RDEx(), 30), MCCGA = list(MCCGA(), 2000))
  for (nm in names(cases)) {
    set.seed(1)
    sol <- solve(rprob, cases[[nm]][[1]], maxiters = cases[[nm]][[2]])
    expect_true(10 * sol$objective < l1, info = nm)
    expect_equal(sol$retcode, "Default", info = nm)
    expect_s3_class(sol$original, "MetaheuristicsState")
  }
  for (alg in list(ECA(), DE(), PSO(), ABC(N = 20, options = Metaheuristics$Options(f_calls_limit = 3000)), CGSA(N = 20),
                   WOA(N = 20))) {
    set.seed(2)
    sol <- solve(rprob, alg, maxiters = 30, use_initial = TRUE)
    expect_true(10 * sol$objective < l1, info = alg$name)
  }
  set.seed(3)
  sol <- solve(rprob, ECA(resize_population = TRUE, options = Metaheuristics$Options(parallel_evaluation = TRUE)), maxiters = 60)
  expect_true(sol$objective < 1e-3)
  set.seed(3)
  expect_true(solve(rprob, ECA(adaptive = TRUE), maxiters = 60)$objective < 1e-3)
  for (s in c("rand2", "best1", "randToBest1", "best2")) {
    set.seed(4)
    expect_true(solve(rprob, DE(strategy = s), maxiters = 60)$objective < 1e-2, info = s)
  }
  set.seed(5)
  expect_true(solve(rprob, DE(F_min = 0.5, F_max = 0.9, CR_min = 0.3, CR_max = 0.9), maxiters = 60)$objective < 1e-2)
})

test_that("the solve arguments are mapped as in OptimizationMetaheuristics", {
  set.seed(1)
  sol <- solve(rprob, DE(), maxiters = 7)
  expect_equal(sol$original$iteration, 7L)
  set.seed(1)
  sol <- solve(rprob, DE(), f_calls_limit = 100, maxiters = 1000)
  expect_true(sol$original$f_calls <= 120)
  expect_equal(sol$original$termination_status_code, "EVALUATIONS_LIMIT")
  sol <- solve(rprob, DE(options = Metaheuristics$Options(seed = 11)), maxiters = 5)
  expect_identical(solve(rprob, DE(options = Metaheuristics$Options(seed = 11)), maxiters = 5)$u, sol$u)
  expect_error(solve(rprob, ECA(), f_optimum = 0),
               "Set f_optimum by directly passing it to Information Structure which is passed to Metaheuristics.Algorithm{ECA} algorithms when calling solve().",
               fixed = TRUE)
  expect_error(solve(rprob, epsilonDE(), x_optimum = c(1, 1)),
               "which is passed to Metaheuristics.Algorithm{\u03B5DE} algorithms", fixed = TRUE)
  expect_error(solve(rprob, ECA(), foo = 1),
               "foo keyword is not a valid option for Metaheuristics.AbstractAlgorithm algorithm.", fixed = TRUE)
  set.seed(1)
  sol <- solve(rprob, ECA(information = Metaheuristics$Information(f_optimum = 0)), abstol = 1e-6, maxiters = 500)
  expect_equal(sol$original$termination_status_code, "ACCURACY_LIMIT")
  expect_true(sol$objective < 1e-6)
  expect_equal(Metaheuristics$termination_status_message(sol), "The desired accuracy was obtained.")
  expect_equal(Metaheuristics$minimizer(sol), sol$u)
  expect_equal(Metaheuristics$minimum(sol), sol$objective)
  expect_equal(nrow(Metaheuristics$positions(sol)), length(sol$original$population))
  set.seed(1)
  sol <- solve(rprob, PSO(termination = Metaheuristics$RobustConvergence(ftol = 1e-3, period = 5)), maxiters = 100)
  expect_true(sol$original$iteration <= 100)
})

test_that("Metaheuristics$optimize, the GA operators and BRKGA work", {
  sphere <- function(x) sum((x - 0.3)^2)
  set.seed(1)
  res <- Metaheuristics$optimize(sphere, rbind(c(-1, -1, -1), c(1, 1, 1)), ECA(options = Metaheuristics$Options(iterations = 60)))
  expect_s3_class(res, "MetaheuristicsState")
  expect_true(Metaheuristics$minimum(res) < 1e-4)
  expect_output(print(res), "Minimum")
  expect_output(print(ECA()), "ECA\\(eta_max=2")
  M <- Metaheuristics
  set.seed(2)
  ga <- GA(N = 20, selection = M$RouletteWheelSelection(N = 20), crossover = M$BinomialCrossover(p = 0.5),
           mutation = M$PolynomialMutation(eta = 20), environmental_selection = M$ElitistReplacement(),
           options = M$Options(iterations = 40))
  expect_true(M$minimum(M$optimize(sphere, rbind(rep(-1, 3), rep(1, 3)), ga)) < 0.05)
  expect_error(M$optimize(sphere, rbind(rep(-1, 3), rep(1, 3)),
                          GA(N = 10, environmental_selection = M$RankAndCrowding(), options = M$Options(iterations = 3))),
               "RankAndCrowding needs a multiobjective population")
  set.seed(3)
  res <- M$optimize(sphere, rbind(rep(0, 3), rep(1, 3)), BRKGA(options = M$Options(iterations = 20)))
  expect_true(M$minimum(res) < 0.05)
  expect_equal(res$parameters$N, 100L)
  # a permutation problem: the identity minimizes the displacement
  disp <- function(x) sum(abs(x - seq_along(x)))
  set.seed(4)
  ga <- GA(N = 30, initializer = M$RandomPermutation(N = 30), crossover = M$OrderCrossover(), mutation = M$SlightMutation(),
           options = M$Options(iterations = 60))
  res <- M$optimize(disp, rbind(rep(1, 6), rep(6, 6)), ga)
  expect_true(all(sort(M$minimizer(res)) == 1:6))
  expect_true(M$minimum(res) <= 4)
  set.seed(5)
  ga <- GA(N = 30, p_mutation = 0.05, p_crossover = 0.5, initializer = M$RandomBinary(N = 30),
           options = M$Options(iterations = 30))
  ones <- function(x) -sum(x)
  expect_true(M$minimum(M$optimize(ones, rbind(rep(0, 8), rep(1, 8)), ga)) <= -7)
  set.seed(6)
  expect_true(M$minimum(M$optimize(sphere, rbind(rep(-1, 2), rep(1, 2)), GA(N = 20, environmental_selection = M$GenerationalReplacement(),
                                                                         options = M$Options(iterations = 20)))) < 0.2)
})

test_that("the local search methods search the given neighborhoods", {
  target <- c(0.3, 0.7, 0.1)
  f <- function(x, p) sum((x - target)^2)
  step <- function(x, i) {
    n <- 2 * length(x)
    if (is.null(i)) i <- sample.int(n, 1)
    if (i > n) return(NULL)
    k <- (i + 1) %/% 2
    x[k] <- x[k] + c(-0.1, 0.1)[(i - 1) %% 2 + 1]
    x
  }
  prob <- OptimizationProblem(f, c(0, 0, 0), lb = rep(-1, 3), ub = rep(1, 3))
  sol <- solve(prob, VND(initial = c(0, 0, 0), neighborhood = list(step)))
  expect_equal(sol$u, target, tolerance = 1e-12)
  set.seed(1)
  sol <- solve(prob, VNS(initial = c(0, 0, 0), neighborhood_shaking = list(step), neighborhood_local = list(step)), maxiters = 5)
  expect_true(sol$objective < 1e-20)
  set.seed(2)
  sol <- solve(prob, GRASP(constructor = function() round(runif(3), 1),
                           local_search = Metaheuristics$BestImproveSearch(neighborhood = step)), maxiters = 5)
  expect_true(sol$objective < 0.05)
  expect_warning(solve(prob, GRASP(constructor = function() runif(3)), maxiters = 2), "BoxConstrainedSpace")
  # two opt and the greedy randomized constructor on a permutation
  disp <- function(x, p) sum(abs(x - seq_along(x)))
  pprob <- OptimizationProblem(disp, as.numeric(5:1), lb = rep(1, 5), ub = rep(5, 5))
  sol <- solve(pprob, VND(initial = c(2, 1, 3, 5, 4), neighborhood = list(Metaheuristics$TwoOptNeighborhood())))
  expect_equal(sol$u, as.numeric(1:5))
  sol <- solve(pprob, VND(initial = as.numeric(5:1), neighborhood = list(Metaheuristics$TwoOptNeighborhood(k = 5))))
  expect_equal(sol$objective, 0)
  set.seed(3)
  cons <- Metaheuristics$GreedyRandomizedConstructor(candidates = as.numeric(1:5), alpha = 0,
                                                     compute_cost = function(cand, constructor, instance) cand)
  sol <- solve(pprob, GRASP(constructor = cons, local_search = Metaheuristics$FirstImproveSearch(
    neighborhood = Metaheuristics$TwoOptNeighborhood())), maxiters = 3)
  expect_equal(sol$u, as.numeric(1:5))
})

test_that("the constraints of the problem are handled as in OptimizationMetaheuristics", {
  r <- sqrt(0.5)
  phi <- seq(0, 2 * pi, length.out = 20000)
  f_ref <- min(vapply(phi, function(a) rosen(r * c(cos(a), sin(a)), c(1, 100)), 0))
  optf <- OptimizationFunction(rosen, cons = function(x, p) x[1]^2 + x[2]^2)
  prob <- OptimizationProblem(optf, c(0, 0), c(1, 100), lb = c(-1, -1), ub = c(1, 1), lcons = -Inf, ucons = 0.5)
  for (alg in list(DE(), PSO(), SHADE(N = 40), epsilonDE())) {
    set.seed(1)
    sol <- solve(prob, alg, maxiters = 150)
    expect_true(sum(sol$u^2) <= 0.5, info = alg$name)
    expect_true(sol$objective <= f_ref + 1e-2, info = alg$name)
  }
  set.seed(1)
  sol <- solve(prob, ECA(), maxiters = 100)
  expect_true(sum(sol$u^2) <= 0.5)
  set.seed(1)
  prob2 <- OptimizationProblem(OptimizationFunction(rosen, cons = function(x, p) c(x[1]^2 + x[2]^2, x[1])), c(0, 0), c(1, 100),
                               lb = c(-2, -2), ub = c(2, 2), lcons = c(2.5, -Inf), ucons = c(Inf, 1))
  sol <- solve(prob2, DE(), maxiters = 150)
  expect_true(sum(sol$u^2) >= 2.5)
  expect_true(sol$u[1] <= 1)
  set.seed(1)
  prob3 <- OptimizationProblem(OptimizationFunction(rosen, cons = function(x, p) x[1] + x[2]), c(0, 0), c(1, 100),
                               lb = c(-2, -2), ub = c(2, 2), lcons = 1, ucons = 1)
  sol <- solve(prob3, DE(options = Metaheuristics$Options(h_tol = 1e-4)), maxiters = 300)
  expect_true(abs(sum(sol$u) - 1) <= 1e-4 + 1e-12)
  line <- seq(-2, 2, length.out = 40001)
  expect_true(sol$objective <= min(vapply(line, function(x) rosen(c(x, 1 - x), c(1, 100)), 0)) + 1e-2)
  expect_error(solve(prob, ABC()), class = "IncompatibleOptimizerError")
  expect_error(solve(prob, SA()), "SA has no method for constrained solutions")
  for (alg in list(CSO(N = 20), GA(N = 20))) {
    set.seed(1)
    sol <- solve(prob, alg, maxiters = 60)
    expect_true(sum(sol$u^2) <= 0.5 + 1e-12, info = alg$name)
  }
})

test_that("the multiobjective algorithms find Pareto sets", {
  M <- Metaheuristics
  sphere2 <- function(x, p) list(c(sum(x^2), sum((x - 2)^2)), 0, 0)
  prob <- OptimizationProblem(MultiObjectiveOptimizationFunction(sphere2), c(0, 0, 0), lb = c(0, 0, 0), ub = c(1, 1, 1))
  algs <- list(M$NSGA2(N = 20), NSGA3(N = 20, partitions = 19L), SPEA2(N = 20), CCMO(M$NSGA2(N = 20, p_m = 0.001)),
               CCMO(NSGA2(population = 20)), MOEAD_DE(M$gen_ref_dirs(2, 19)), SMS_EMOA(N = 20))
  for (alg in algs) {
    set.seed(42)
    sol <- if (alg$name == "CCMO") solve(prob, alg, maxiters = 15) else solve(prob, alg, maxiters = 15, use_initial = TRUE)
    expect_true(NROW(sol$u) > 0, info = alg$name)
    # on this problem the Pareto set is the segment x1 = x2 = x3 in [0, 1]
    # (a random point has a mean spread of 0.5)
    expect_true(mean(apply(sol$u, 1, function(x) max(x) - min(x))) < 0.3, info = alg$name)
    expect_equal(nrow(M$pareto_front(sol)), NROW(sol$u), info = alg$name)
  }
  three <- function(x, p) c(sum(x^2), sum((x - 1)^2), sum((x - c(1, 0, 0))^2))
  prob3 <- OptimizationProblem(MultiObjectiveOptimizationFunction(three), c(0, 0, 0), lb = c(0, 0, 0), ub = c(1, 1, 1))
  for (alg in list(NSGA3(N = 15, partitions = 4L), SMS_EMOA(N = 10), MOEAD_DE(M$gen_ref_dirs(3, 4)), SPEA2(N = 15))) {
    set.seed(7)
    sol <- solve(prob3, alg, maxiters = 5)
    expect_equal(ncol(sol$objective), 3L, info = alg$name)
  }
  set.seed(8)
  sol <- solve(prob, SMS_EMOA(N = 10, options = M$Options(parallel_evaluation = TRUE)), maxiters = 5)
  expect_true(NROW(sol$u) > 0)
  # constraints of the problem are appended to those of the objective
  mof <- MultiObjectiveOptimizationFunction(function(x, p) list(c(sum(x^2), sum((x - 2)^2)), 0, 0),
                                            cons = function(x, p) x[1])
  cprob <- OptimizationProblem(mof, c(0, 0), lb = c(-1, -1), ub = c(3, 3), lcons = -Inf, ucons = 0.5)
  set.seed(1)
  sol <- solve(cprob, M$NSGA2(N = 30), maxiters = 30)
  expect_true(all(M$positions(sol$original)[, 1] <= 0.5))
  expect_error(MOEAD_DE(M$gen_ref_dirs(2, 1)), "at least 3 weight vectors")
  expect_output(print(sol$original), "Nondominated solutions")
})
