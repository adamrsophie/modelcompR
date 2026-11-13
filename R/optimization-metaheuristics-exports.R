# The algorithm constructors of OptimizationMetaheuristics (the exports of
# Metaheuristics.jl that Optimization.jl re-exports) and of the operators
# and helpers in the `Metaheuristics` list.

.mh_bits <- function(...) list(...)

#' Single objective metaheuristics of Metaheuristics.jl
#'
#' The population based optimizers of Metaheuristics.jl that
#' OptimizationMetaheuristics makes available to [solve()], ported to R.
#' They need box bounds (`lb` and `ub`) and use no derivatives.
#'
#' * `ECA()`: evolutionary centers algorithm.
#' * `DE()`: differential evolution (strategies `"rand1"`, `"rand2"`,
#'   `"best1"`, `"randToBest1"` and `"best2"`).
#' * `epsilonDE()`: the epsilon constrained differential evolution (named
#'   with a Greek epsilon in Metaheuristics.jl).
#' * `SHADE()`: success history based adaptive differential evolution.
#' * `RDEx()`: the reconstructed differential evolution with elite based
#'   mutation and adaptive memories.
#' * `PSO()`: particle swarm optimization.
#' * `ABC()`: artificial bee colony.
#' * `SA()`: simulated annealing.
#' * `WOA()`: whale optimization algorithm.
#' * `CGSA()`: chaotic gravitational search algorithm.
#' * `CSO()`: competitive swarm optimizer.
#' * `GA()`: a genetic algorithm assembled from operators of
#'   [Metaheuristics] (initializer, selection, crossover, mutation and
#'   environmental selection); `BRKGA()` is the biased random key genetic
#'   algorithm built from them.
#' * `MCCGA()`: machine coded compact genetic algorithm, with a final
#'   Hooke and Jeeves search.
#'
#' ECA, DE, epsilonDE, SHADE, PSO, SA, CSO, GA and BRKGA accept the
#' constraints of the problem (`cons` with `lcons` and `ucons`), which become
#' the inequality and equality constraints of Metaheuristics.jl and are
#' handled by its feasibility rules. The `solve` arguments map as in
#' OptimizationMetaheuristics: `maxiters` to the number of iterations,
#' `maxtime` to the time limit, `abstol` to `f_tol`; every field of
#' `Metaheuristics$Options()` can be given too, and `use_initial = TRUE`
#' seeds the population with `u0`. The random numbers are R's: call
#' [set.seed()] (or give `seed`) for reproducible runs.
#'
#' @param N Population size (`0` lets the algorithm choose from the
#'   dimension).
#' @param eta_max,K,N_init,p_exploit,p_bin,p_cr,epsilon,adaptive,resize_population
#'   Parameters of ECA: the step scale, the size of the neighborhood, the
#'   initial population, the exploitation share of the budget, the
#'   crossover probability, the crossover probabilities per coordinate, an
#'   unused tolerance, the adaptation of `p_cr` and the shrinking of the
#'   population.
#' @param F,CR,CR_min,CR_max,F_min,F_max,strategy Parameters of DE: the step
#'   size, the crossover rate, the ranges they are drawn from in each
#'   iteration when the minimum is below the maximum, and the mutation
#'   strategy.
#' @param cp Exponent of the epsilon level control of `epsilonDE()`.
#' @param population_size,memory_size Parameters of RDEx: the population
#'   per dimension and the size of the success memories.
#' @param C1,C2,omega,v,flock Parameters of PSO: the cognitive and social
#'   factors, the inertia, and the initial velocities and swarm (empty to
#'   start from rest).
#' @param Ne,No,limit Parameters of ABC: employed and onlooker bees and the
#'   trials before a bee becomes a scout.
#' @param x_initial,tol_fun Parameters of SA: the starting point (random
#'   when empty) and the tolerance of the acceptance rule; `N` is the number
#'   of moves per iteration.
#' @param chValueInitial,chaosIndex,ElitistCheck,Rpower,Rnorm,wMax,wMin,X,V,fitness
#'   Parameters of CGSA: the chaotic map value and index (1 to 10), the
#'   elitism, the distance power and norm, the range of the chaotic term and
#'   the initial positions, velocities and fitness (empty for random).
#' @param phi Social factor of CSO.
#' @param p_mutation,p_crossover Probabilities of the bit operators of
#'   Metaheuristics.GA: when given (and `mutation` or `crossover` is `NULL`)
#'   `GA()` uses `BitFlipMutation(p = p_mutation)` or
#'   `UniformCrossover(p = p_crossover)`.
#' @param initializer,selection,crossover,mutation,environmental_selection
#'   Operators of GA, made with the functions of [Metaheuristics]. `NULL`
#'   gives `RandomInBounds(N)`, `TournamentSelection(N = N)`, `SBX()`,
#'   `PolynomialMutation()` and `ElitistReplacement()`, the operators
#'   Metaheuristics.jl chooses for a box search space.
#' @param num_elites,num_mutants,num_offsprings,bias Parameters of BRKGA.
#' @param maxsamples,use_local_search Parameters of MCCGA: the samples of
#'   the initial probabilities and the final local search; `mutation` is the
#'   step of the probabilities.
#' @param information `Metaheuristics$Information()`, the known optimum.
#' @param options `Metaheuristics$Options()`.
#' @param termination A termination criterion of [Metaheuristics] (or a
#'   `Termination()` of several).
#' @return An algorithm object for [solve()].
#' @examples
#' rosen <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
#' prob <- OptimizationProblem(rosen, c(0, 0), c(1, 100), lb = c(-1, -1), ub = c(1, 1))
#' set.seed(1)
#' sol <- solve(prob, DE(), maxiters = 200)
#' round(sol$u, 2)
#' set.seed(1)
#' solve(prob, ECA(), maxiters = 100)$objective < 0.1
#' @name metaheuristics_single
NULL

#' @rdname metaheuristics_single
#' @export
ECA <- function(eta_max = 2, K = 7L, N = 0L, N_init = N, p_exploit = 0.95, p_bin = 0.02, p_cr = numeric(0),
                epsilon = 0, adaptive = FALSE, resize_population = FALSE, information = NULL, options = NULL,
                termination = NULL)
  .mh_algorithm("ECA", list(eta_max = eta_max, K = as.integer(K), N = as.integer(N), N_init = as.integer(N),
                            p_exploit = p_exploit, p_bin = p_bin, p_cr = p_cr, epsilon = epsilon, adaptive = adaptive,
                            resize_population = resize_population),
                information, options, termination, cons = TRUE)

.mh_de_parameters <- function(N, F, CR, CR_min, CR_max, F_min, F_max, strategy)
  list(N = as.integer(N), F = as.numeric(F), CR = as.numeric(CR), CR_min = as.numeric(CR_min), CR_max = as.numeric(CR_max),
       F_min = as.numeric(F_min), F_max = as.numeric(F_max), strategy = as.character(strategy))

#' @rdname metaheuristics_single
#' @export
DE <- function(N = 0L, F = 0.7, CR = 0.5, CR_min = CR, CR_max = CR, F_min = F, F_max = F, strategy = "rand1",
               information = NULL, options = NULL, termination = NULL)
  .mh_algorithm("DE", .mh_de_parameters(N, F, CR, CR_min, CR_max, F_min, F_max, strategy), information, options,
                termination, cons = TRUE)

#' @rdname metaheuristics_single
#' @export
epsilonDE <- function(cp = 5, N = 0L, F = 0.7, CR = 0.5, CR_min = CR, CR_max = CR, F_min = F, F_max = F,
                      strategy = "rand1", information = NULL, options = NULL, termination = NULL) {
  a <- .mh_algorithm("epsilonDE", list(de = .mh_de_parameters(N, F, CR, CR_min, CR_max, F_min, F_max, strategy), N = 0L,
                                       eps = 0, eps0 = 0, Tc = 0, cp = cp),
                     information, options, termination, cons = TRUE)
  a$jl_name <- "\u03B5DE"
  a
}

#' @rdname metaheuristics_single
#' @export
SHADE <- function(N = 100L, information = NULL, options = NULL, termination = NULL)
  .mh_algorithm("SHADE", list(N = as.integer(N), archive = list(), MF = numeric(0), F = numeric(0), MCR = numeric(0),
                              CR = numeric(0), k = 1L),
                information, options, termination, cons = TRUE)

#' @rdname metaheuristics_single
#' @export
RDEx <- function(population_size = 18L, memory_size = 5L, information = NULL, options = NULL)
  .mh_algorithm("RDEx", list(N = 0L, population_size = as.integer(population_size), memory_size = as.integer(memory_size)),
                information, options)

#' @rdname metaheuristics_single
#' @export
PSO <- function(N = 0L, C1 = 2, C2 = 2, omega = 0.8, v = numeric(0), flock = list(), information = NULL, options = NULL,
                termination = NULL)
  .mh_algorithm("PSO", list(N = as.integer(N), C1 = C1, C2 = C2, omega = omega, v = v, flock = flock),
                information, options, termination, cons = TRUE)

#' @rdname metaheuristics_single
#' @export
ABC <- function(N = 50L, Ne = (N + 1L) %/% 2L, No = (N + 1L) %/% 2L, limit = 10L, information = NULL, options = NULL,
                termination = NULL)
  .mh_algorithm("ABC", list(N = as.integer(N), Ne = as.integer(Ne), No = as.integer(No), limit = as.integer(limit)),
                information, options, termination)

#' @rdname metaheuristics_single
#' @export
SA <- function(x_initial = numeric(0), N = 500L, tol_fun = 1e-4, information = NULL, options = NULL, termination = NULL)
  .mh_algorithm("SA", list(N = as.integer(N), x_initial = as.numeric(x_initial), tol_fun = tol_fun, x = numeric(0), fx = Inf),
                information, options, termination, cons = TRUE)

#' @rdname metaheuristics_single
#' @export
WOA <- function(N = 30L, information = NULL, options = NULL, termination = NULL)
  .mh_algorithm("WOA", list(N = as.integer(N)), information, options, termination)

#' @rdname metaheuristics_single
#' @export
CGSA <- function(N = 30L, chValueInitial = 20, chaosIndex = 8, ElitistCheck = 1L, Rpower = 1L, Rnorm = 2L,
                 wMax = chValueInitial, wMin = 1e-10, X = matrix(0, 0, 0), V = matrix(0, 0, 0), fitness = numeric(0),
                 information = NULL, options = NULL, termination = NULL)
  .mh_algorithm("CGSA", list(N = as.integer(N), chValueInitial = chValueInitial, chaosIndex = chaosIndex,
                             ElitistCheck = ElitistCheck, Rpower = Rpower, Rnorm = Rnorm, wMax = wMax, wMin = wMin,
                             X = X, V = V, fitness = fitness),
                information, options, termination)

#' @rdname metaheuristics_single
#' @export
CSO <- function(N = 400L, phi = 0.1, information = NULL, options = NULL, termination = NULL)
  .mh_algorithm("CSO", list(N = as.integer(N), phi = phi, V_loser = matrix(0, 0, 0)), information, options, termination,
                cons = TRUE)

#' @rdname metaheuristics_single
#' @export
GA <- function(N = 100L, p_mutation = NULL, p_crossover = NULL, initializer = NULL, selection = NULL, crossover = NULL,
               mutation = NULL, environmental_selection = NULL, information = NULL, options = NULL,
               termination = NULL) {
  N <- as.integer(N)
  crossover <- crossover %mh_or% (if (is.null(p_crossover)) .mh_SBX() else .mh_UniformCrossover(p = p_crossover))
  mutation <- mutation %mh_or% (if (is.null(p_mutation)) .mh_PolynomialMutation() else .mh_BitFlipMutation(p = p_mutation))
  .mh_algorithm("GA", list(N = N, initializer = initializer %mh_or% .mh_RandomInBounds(N),
                           selection = selection %mh_or% .mh_TournamentSelection(N = N), crossover = crossover,
                           mutation = mutation, environmental_selection = environmental_selection %mh_or% .mh_op("ElitistReplacement")),
                information, options, termination, cons = TRUE)
}

#' @rdname metaheuristics_single
#' @export
BRKGA <- function(num_elites = 20L, num_mutants = 10L, num_offsprings = 70L, N = num_elites + num_mutants + num_offsprings,
                  bias = 0.7, information = NULL, options = NULL, termination = NULL) {
  a <- GA(N = N, initializer = .mh_RandomInBounds(N), selection = .mh_op("BiasedSelection", num_elites = as.integer(num_elites),
                                                                          num_offsprings = as.integer(num_offsprings)),
          crossover = .mh_op("BinomialCrossover", p = bias, n_offsprings = 1L),
          mutation = .mh_op("InsertRandomMutation", num_mutants = as.integer(num_mutants)),
          environmental_selection = .mh_op("ElitistReplacement"), information = information, options = options,
          termination = termination)
  class(a)[1L] <- "BRKGA"
  a$name <- "BRKGA"
  a$impl <- "GA"
  a$jl_name <- "GA"
  a
}

#' @rdname metaheuristics_single
#' @export
MCCGA <- function(N = 100L, maxsamples = 10000L, mutation = 1 / N, use_local_search = TRUE, information = NULL,
                  options = NULL, termination = NULL)
  .mh_algorithm("MCCGA", list(N = as.integer(N), maxsamples = as.integer(maxsamples), mutation = mutation,
                              probvector = numeric(0), use_local_search = use_local_search),
                information, options, termination)

#' Multiobjective metaheuristics of Metaheuristics.jl
#'
#' The multiobjective optimizers of Metaheuristics.jl that
#' OptimizationMetaheuristics makes available, ported to R, for problems
#' built from a [MultiObjectiveOptimizationFunction()] with box bounds. The
#' objective returns the vector of the costs, or a list `(F, g, h)` with the
#' inequality (`g <= 0`) and equality (`h = 0`) constraints as in
#' Metaheuristics.jl; constraints given to the problem are appended.
#'
#' * `NSGA3()`: reference point based nondominated sorting.
#' * `SPEA2()`: strength Pareto evolutionary algorithm.
#' * `SMS_EMOA()`: hypervolume selection (exact contributions for two and
#'   three objectives, Monte Carlo with `n_samples` points beyond).
#' * `MOEAD_DE()`: decomposition with differential evolution for the
#'   given weight vectors (see `Metaheuristics$gen_ref_dirs()`).
#' * `CCMO()`: coevolutionary constrained optimization around an NSGA2
#'   (`Metaheuristics$NSGA2()` or [NSGA2()]).
#'
#' The solution's `u` is the set of the nondominated solutions of the final
#' population and `objective` their costs; `original` is the final state.
#'
#' @param N Population size.
#' @param eta_cr,p_cr,eta_m,p_m Distribution index and probability of the
#'   simulated binary crossover and of the polynomial mutation (`p_m < 0`
#'   means `1 / n`).
#' @param partitions,reference_points Partitions of the simplex for the
#'   reference points of NSGA3, or the points themselves (a list).
#' @param n_samples Monte Carlo samples of SMS_EMOA beyond three objectives.
#' @param weights The weight vectors of MOEAD_DE (a list); the population
#'   has one solution per weight.
#' @param F,CR Step size and crossover rate of MOEAD_DE.
#' @param lambda Ignored, as in Metaheuristics.jl (the weights are used).
#' @param eta Distribution index of the MOEAD_DE mutation.
#' @param T,delta,n_r Neighborhood size, probability of mating in the
#'   neighborhood and maximal number of replaced solutions.
#' @param z,B Initial ideal point and neighborhoods (computed when empty).
#' @param s1,s2 Penalty factors of the constraint violation.
#' @param base_algorithm The NSGA2 algorithm wrapped by CCMO.
#' @param information,options,termination As for [metaheuristics_single].
#' @return An algorithm object for [solve()].
#' @examples
#' f <- function(x, p) c(sum(x^2), sum((x - 2)^2))
#' prob <- OptimizationProblem(MultiObjectiveOptimizationFunction(f), c(0, 0),
#'   lb = c(0, 0), ub = c(2, 2))
#' set.seed(1)
#' sol <- solve(prob, SPEA2(N = 20), maxiters = 20)
#' dim(sol$u)
#' @name metaheuristics_multi
NULL

.mh_nsga_parameters <- function(N, eta_cr, p_cr, eta_m, p_m)
  list(N = as.integer(N), eta_cr = as.numeric(eta_cr), p_cr = as.numeric(p_cr), eta_m = as.numeric(eta_m),
       p_m = as.numeric(p_m))

.mh_NSGA2 <- function(N = 100L, eta_cr = 20, p_cr = 0.9, eta_m = 20, p_m = -1, information = NULL, options = NULL,
                      termination = NULL)
  .mh_algorithm("NSGA2", .mh_nsga_parameters(N, eta_cr, p_cr, eta_m, p_m), information, options, termination, cons = TRUE)

#' @rdname metaheuristics_multi
#' @export
NSGA3 <- function(N = 100L, eta_cr = 20, p_cr = 0.9, eta_m = 20, p_m = -1, partitions = 12L, reference_points = list(),
                  information = NULL, options = NULL, termination = NULL)
  .mh_algorithm("NSGA3", c(.mh_nsga_parameters(N, eta_cr, p_cr, eta_m, p_m),
                           list(partitions = as.integer(partitions), reference_points = reference_points)),
                information, options, termination, cons = TRUE)

#' @rdname metaheuristics_multi
#' @export
SPEA2 <- function(N = 100L, eta_cr = 20, p_cr = 0.9, eta_m = 20, p_m = -1, information = NULL, options = NULL,
                  termination = NULL)
  .mh_algorithm("SPEA2", c(.mh_nsga_parameters(N, eta_cr, p_cr, eta_m, p_m), list(fitness = numeric(0))),
                information, options, termination, cons = TRUE)

#' @rdname metaheuristics_multi
#' @export
SMS_EMOA <- function(N = 100L, eta_cr = 20, p_cr = 0.9, eta_m = 20, p_m = -1, n_samples = 10000L, information = NULL,
                     options = NULL, termination = NULL)
  .mh_algorithm("SMS_EMOA", c(.mh_nsga_parameters(N, eta_cr, p_cr, eta_m, p_m), list(n_samples = as.integer(n_samples))),
                information, options, termination, cons = TRUE)

#' @rdname metaheuristics_multi
#' @export
MOEAD_DE <- function(weights, F = 0.5, CR = 1, lambda = list(), eta = 20, p_m = -1,
                     T = min(max(round(0.2 * length(weights)), 2), length(weights) - 1), delta = 0.9,
                     n_r = max(1, round(0.02 * length(weights))), z = numeric(0), B = list(), s1 = 0.01, s2 = 20,
                     information = NULL, options = NULL, termination = NULL) {
  if (!length(weights)) stop("Provide weights points", call. = FALSE)
  weights <- lapply(weights, as.numeric)
  N <- length(weights)
  if (N < 3L) stop("MOEAD_DE needs at least 3 weight vectors (given ", N, ").", call. = FALSE)
  if (!(2 <= T && T <= N - 1))
    stop("MOEAD_DE needs 2 \u2264 T \u2264 N - 1 = ", N - 1, " (given T = ", T, ").", call. = FALSE)
  if (!(n_r >= 1)) stop("MOEAD_DE needs n_r \u2265 1 (given n_r = ", n_r, ").", call. = FALSE)
  M <- length(weights[[1L]])
  if (!length(z)) z <- rep(Inf, M)
  if (!length(B)) B <- .mh_moead_neighbors(weights, as.integer(T))
  .mh_algorithm("MOEAD_DE", list(nobjectives = M, N = N, F = F, CR = CR, lambda = weights, eta = eta, p_m = p_m,
                                 T = as.integer(T), delta = delta, n_r = n_r, z = as.numeric(z), B = B, s1 = s1, s2 = s2,
                                 tau = 0),
                information, options, termination)
}

#' @rdname metaheuristics_multi
#' @export
CCMO <- function(base_algorithm, information = NULL, options = NULL, termination = NULL) {
  base <- if (inherits(base_algorithm, "MetaheuristicsAlgorithm") && base_algorithm$name == "NSGA2") {
    base_algorithm$parameters
  } else if (inherits(base_algorithm, "NSGA2")) {
    .mh_nsga_parameters(base_algorithm$population, base_algorithm$eta_c, base_algorithm$crossover_rate,
                        base_algorithm$eta_m, base_algorithm$mutation_rate %mh_or% -1)
  } else {
    stop("CCMO is implemented for an NSGA2 base algorithm", call. = FALSE)
  }
  .mh_algorithm("CCMO", list(base = base, phelper = list(), fitness = numeric(0), fitness_helper = numeric(0)),
                information, options, termination, cons = TRUE)
}

#' Combinatorial metaheuristics of Metaheuristics.jl
#'
#' The local search methods of Metaheuristics.jl that OptimizationMetaheuristics
#' exports, ported to R. They search the neighborhoods that the user
#' defines; in R a neighborhood is a function `function(x, i)` that returns
#' the `i`-th neighbor of `x` (`NULL` after the last one) or a random
#' neighbor when `i` is `NULL`, or `Metaheuristics$TwoOptNeighborhood()`.
#'
#' * `GRASP()`: greedy randomized adaptive search: each iteration builds a
#'   point with `constructor` (a function of no arguments, or
#'   `Metaheuristics$GreedyRandomizedConstructor()`) and improves it with
#'   `local_search`.
#' * `VND()`: variable neighborhood descent over the list `neighborhood`.
#' * `VNS()`: variable neighborhood search: shaking in
#'   `neighborhood_shaking`, then a VND over `neighborhood_local`.
#'
#' @param initial The starting point (random in the bounds when `NULL`).
#' @param constructor The randomized greedy constructor of GRASP.
#' @param local_search `Metaheuristics$BestImproveSearch()` or
#'   `Metaheuristics$FirstImproveSearch()` (for GRASP with a
#'   `neighborhood`), or for GRASP a function of a point that returns the
#'   improved point.
#' @param neighborhood A list of neighborhoods of VND.
#' @param neighborhood_shaking,neighborhood_local Lists of neighborhoods of
#'   VNS.
#' @param neighborhood_change `Metaheuristics$SequentialChange()` or
#'   `Metaheuristics$CyclicChange()`.
#' @param options,information As for [metaheuristics_single].
#' @return An algorithm object for [solve()].
#' @examples
#' f <- function(x, p) sum((x - c(0.3, 0.7))^2)
#' step <- function(x, i) {
#'   if (is.null(i)) i <- sample.int(4, 1)
#'   if (i > 4) return(NULL)
#'   x[(i + 1) %/% 2] <- x[(i + 1) %/% 2] + c(-0.1, 0.1)[(i - 1) %% 2 + 1]
#'   x
#' }
#' prob <- OptimizationProblem(f, c(0, 0), lb = c(-1, -1), ub = c(1, 1))
#' sol <- solve(prob, VND(initial = c(0, 0), neighborhood = list(step)))
#' sol$u
#' @name metaheuristics_combinatorial
NULL

#' @rdname metaheuristics_combinatorial
#' @export
GRASP <- function(initial = NULL, constructor = NULL, local_search = .mh_op("BestImproveSearch"), options = NULL,
                  information = NULL) {
  if (is.null(constructor)) stop("Provide a constructor.", call. = FALSE)
  .mh_algorithm("GRASP", list(initial = initial, constructor = constructor, local_search = local_search),
                information, options)
}

#' @rdname metaheuristics_combinatorial
#' @export
VND <- function(initial = NULL, neighborhood = NULL, local_search = .mh_op("BestImproveSearch"), options = NULL,
                information = NULL)
  .mh_algorithm("VND", list(initial = initial, neighborhood = neighborhood, local_search = local_search),
                information, options)

#' @rdname metaheuristics_combinatorial
#' @export
VNS <- function(initial = NULL, neighborhood_shaking = NULL, neighborhood_local = NULL,
                local_search = .mh_op("FirstImproveSearch"), neighborhood_change = .mh_op("SequentialChange"),
                options = NULL, information = NULL) {
  if (is.null(neighborhood_shaking) && is.null(neighborhood_local))
    stop("Provide neighborhood_shaking and neighborhood_local.", call. = FALSE)
  .mh_algorithm("VNS", list(initial = initial, neighborhood_shaking = neighborhood_shaking,
                            neighborhood_local = neighborhood_local, local_search = local_search,
                            neighborhood_change = neighborhood_change),
                information, options)
}

# ---------------------------------------------------------------------------
# Operators and helpers of the Metaheuristics list

.mh_SBX <- function(eta = 15, p = 0.9, bounds = NULL)
  .mh_op("SBX", eta = eta, p = p, lb = if (is.null(bounds)) NULL else bounds[1L, ], ub = if (is.null(bounds)) NULL else bounds[2L, ])

.mh_PolynomialMutation <- function(eta = 15, bounds = NULL, p = NULL)
  .mh_op("PolynomialMutation", eta = eta, p = p, lb = if (is.null(bounds)) NULL else bounds[1L, ],
         ub = if (is.null(bounds)) NULL else bounds[2L, ])

.mh_UniformCrossover <- function(p = 0.5) .mh_op("UniformCrossover", p = p)

.mh_BinomialCrossover <- function(p = 0.5, n_offsprings = 2L) {
  if (!n_offsprings %in% c(1L, 2L)) stop("n_offsprings must be 1 or 2, got ", n_offsprings, call. = FALSE)
  .mh_op("BinomialCrossover", p = p, n_offsprings = as.integer(n_offsprings))
}

.mh_BitFlipMutation <- function(p = 1e-2) .mh_op("BitFlipMutation", p = p)

.mh_TournamentSelection <- function(K = 2L, N = 0L) .mh_op("TournamentSelection", K = as.integer(K), N = as.integer(N))

.mh_RandomInBounds <- function(N = 0L) .mh_op("RandomInBounds", N = as.integer(N))

.mh_GreedyRandomizedConstructor <- function(candidates, instance = NULL, alpha = 0.6,
                                            compute_cost = function(candidates, constructor, instance) numeric(length(candidates)))
  .mh_op("GreedyRandomizedConstructor", candidates = candidates, instance = instance, alpha = alpha,
         compute_cost = compute_cost)

.mh_state_of <- function(s) {
  if (inherits(s, "OptimizationSolution") || inherits(s, "MultiObjectiveOptimizationSolution")) s <- s$original
  s
}

.mh_bounds <- function(bounds) {
  if (is.list(bounds) && !is.null(bounds$lb)) return(list(lb = as.numeric(bounds$lb), ub = as.numeric(bounds$ub)))
  bounds <- as.matrix(bounds)
  if (nrow(bounds) == 2L) return(list(lb = bounds[1L, ], ub = bounds[2L, ]))
  if (ncol(bounds) == 2L) return(list(lb = bounds[, 1L], ub = bounds[, 2L]))
  stop("bounds must be a matrix with the rows lb and ub", call. = FALSE)
}

.mh_mo_names <- c("NSGA2", "NSGA3", "SPEA2", "SMS_EMOA", "MOEAD_DE", "CCMO")

.mh_api_optimize <- function(f, bounds, method) {
  b <- .mh_bounds(bounds)
  .mh_optimize(f, b$lb, b$ub, method, method$name %in% .mh_mo_names)
}

.mh_api_pop <- function(s) {
  s <- .mh_state_of(s)
  if (inherits(s, "MetaheuristicsState")) s$population else s
}

#' The Metaheuristics module
#'
#' `Metaheuristics` is the module binding that OptimizationMetaheuristics
#' exports: a list with the parts of Metaheuristics.jl used with its
#' algorithms, reached as `Metaheuristics$Options(...)`. It contains
#'
#' * the algorithms (also exported on their own, see
#'   [metaheuristics_single], [metaheuristics_multi] and
#'   [metaheuristics_combinatorial]) and `NSGA2(N, eta_cr, p_cr, eta_m,
#'   p_m)`, the NSGA2 of Metaheuristics.jl;
#' * `Options(x_tol, f_tol, g_tol, h_tol, f_tol_rel, f_tol_abs,
#'   f_calls_limit, iterations, time_limit, store_convergence,
#'   parallel_evaluation, seed, ...)` and `Information(f_optimum,
#'   x_optimum)`;
#' * the termination criteria `CheckConvergence()`, `RobustConvergence()`,
#'   `AbsoluteFunctionConvergence()`, `RelativeFunctionConvergence()`,
#'   `SmallStandardDeviation()`, `RelativeParameterConvergence()` and
#'   `Termination(checkall, checkany)`;
#' * `optimize(f, bounds, method)`, which minimizes `f(x)` in the box
#'   `bounds` (a matrix with the rows `lb` and `ub`) and returns the final
#'   state, and `boxconstraints(lb, ub)`;
#' * the readers of a state (or of a solution of [solve()]) `minimizer`,
#'   `minimum`, `positions`, `fvals`, `nfes`, `pareto_front`,
#'   `get_non_dominated_solutions`, `termination_status_message`, and
#'   `non_dominated_sort` (the ranks of a population), `gen_ref_dirs(m, H)`
#'   and `gen_weights(m, H)` (the points of the simplex);
#' * the operators of `GA()`: `RandomInBounds`, `RandomPermutation`,
#'   `RandomBinary`, `TournamentSelection`, `RouletteWheelSelection`,
#'   `BiasedSelection`, `SBX`, `UniformCrossover`, `BinomialCrossover`,
#'   `OrderCrossover`, `PolynomialMutation`, `BitFlipMutation`,
#'   `SlightMutation`, `InsertRandomMutation`, `GenerationalReplacement`,
#'   `ElitistReplacement` and `RankAndCrowding`;
#' * the parts of the local search methods: `BestImproveSearch`,
#'   `FirstImproveSearch`, `TwoOptNeighborhood(k)`,
#'   `GreedyRandomizedConstructor(candidates, instance, alpha,
#'   compute_cost)`, `SequentialChange` and `CyclicChange`.
#'
#' @format A named list of functions.
#' @examples
#' set.seed(1)
#' res <- Metaheuristics$optimize(function(x) sum(x^2), rbind(c(-1, -1), c(1, 1)),
#'                               ECA(options = Metaheuristics$Options(iterations = 50)))
#' Metaheuristics$minimum(res) < 1e-3
#' Metaheuristics$gen_ref_dirs(2, 4)
#' @export
Metaheuristics <- list(
  ABC = ABC, BRKGA = BRKGA, CCMO = CCMO, CGSA = CGSA, CSO = CSO, DE = DE, ECA = ECA, GA = GA, GRASP = GRASP,
  MCCGA = MCCGA, MOEAD_DE = MOEAD_DE, NSGA2 = .mh_NSGA2, NSGA3 = NSGA3, PSO = PSO, RDEx = RDEx, SA = SA,
  SHADE = SHADE, SMS_EMOA = SMS_EMOA, SPEA2 = SPEA2, VND = VND, VNS = VNS, WOA = WOA, epsilonDE = epsilonDE,
  Options = .mh_Options, Information = .mh_Information,
  CheckConvergence = .mh_CheckConvergence, RobustConvergence = .mh_RobustConvergence,
  AbsoluteFunctionConvergence = function(ftol = 0) .mh_criterion("AbsoluteFunctionConvergence", ftol = ftol),
  RelativeFunctionConvergence = function(ftol = .Machine$double.eps) .mh_criterion("RelativeFunctionConvergence", ftol = ftol),
  SmallStandardDeviation = function(ftol = 1e-6) .mh_criterion("SmallStandardDeviation", ftol = ftol),
  RelativeParameterConvergence = function(xtol = 1e-8) .mh_criterion("RelativeParameterConvergence", xtol = xtol),
  Termination = .mh_Termination,
  optimize = .mh_api_optimize,
  boxconstraints = function(lb, ub) list(lb = as.numeric(lb), ub = as.numeric(ub)),
  minimizer = function(s) .mh_state_of(s)$best_sol$x,
  minimum = function(s) .mh_state_of(s)$best_sol$f,
  positions = function(s) .mh_positions(.mh_api_pop(s)),
  fvals = function(s) .mh_fvals(.mh_api_pop(s)),
  nfes = function(s) .mh_state_of(s)$f_calls,
  get_non_dominated_solutions = function(s) {
    pop <- .mh_api_pop(s)
    pop[.mh_get_non_dominated_perm(pop)]
  },
  pareto_front = function(s) {
    pop <- .mh_api_pop(s)
    .mh_fvals(pop[.mh_get_non_dominated_perm(pop)])
  },
  non_dominated_sort = function(s) .mh_non_dominated_sort(.mh_api_pop(s)),
  termination_status_message = function(s) .mh_termination_message(.mh_state_of(s)$termination_status_code),
  gen_ref_dirs = function(m, H) .mh_gen_weights(as.integer(m), as.integer(H)),
  gen_weights = function(m, H) .mh_gen_weights(as.integer(m), as.integer(H)),
  RandomInBounds = .mh_RandomInBounds,
  RandomPermutation = function(N = 0L) .mh_op("RandomPermutation", N = as.integer(N)),
  RandomBinary = function(N = 0L) .mh_op("RandomBinary", N = as.integer(N)),
  TournamentSelection = .mh_TournamentSelection,
  RouletteWheelSelection = function(N = 0L) .mh_op("RouletteWheelSelection", N = as.integer(N)),
  BiasedSelection = function(num_elites, num_offsprings)
    .mh_op("BiasedSelection", num_elites = as.integer(num_elites), num_offsprings = as.integer(num_offsprings)),
  SBX = .mh_SBX, UniformCrossover = .mh_UniformCrossover, BinomialCrossover = .mh_BinomialCrossover,
  OrderCrossover = function() .mh_op("OrderCrossover"),
  PolynomialMutation = .mh_PolynomialMutation, BitFlipMutation = .mh_BitFlipMutation,
  SlightMutation = function() .mh_op("SlightMutation"),
  InsertRandomMutation = function(num_mutants) .mh_op("InsertRandomMutation", num_mutants = as.integer(num_mutants)),
  GenerationalReplacement = function() .mh_op("GenerationalReplacement"),
  ElitistReplacement = function() .mh_op("ElitistReplacement"),
  RankAndCrowding = function() .mh_op("RankAndCrowding"),
  BestImproveSearch = function(neighborhood = NULL) .mh_op("BestImproveSearch", neighborhood = neighborhood),
  FirstImproveSearch = function(neighborhood = NULL) .mh_op("FirstImproveSearch", neighborhood = neighborhood),
  TwoOptNeighborhood = function(k = 2L) .mh_op("TwoOptNeighborhood", k = as.integer(k)),
  GreedyRandomizedConstructor = .mh_GreedyRandomizedConstructor,
  SequentialChange = function() .mh_op("SequentialChange"),
  CyclicChange = function() .mh_op("CyclicChange"))
