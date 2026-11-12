# Levy areas and iterated integrals against StochasticDiffEqLevyArea and the
# evaluators of StochasticDiffEqCore. Reference numbers from
# tools/levyarea-probe.jl, with the fixed random numbers
# z(k) = mod(37 k, 101) / 50 - 1 drawn in the order of generate_coefficients.

zgen <- function() {
  k <- 0L
  function(n) {
    v <- ((37L * (k + seq_len(n))) %% 101L) / 50 - 1
    k <<- k + n
    v
  }
}
cf <- function(m, n, alg) generate_coefficients(m, n, alg, rng = zgen())
W <- c(0.3, -0.7, 1.1)

test_that("the four Levy area algorithms match StochasticDiffEqLevyArea", {
  ref <- list(
    Fourier = c(0.0, 0.0721894068981071, 0.10258505331300247, -0.0721894068981071, 0.0, 0.16220187996758675,
                -0.10258505331300247, -0.16220187996758675, 0.0),
    Milstein = c(0.0, 0.07774831530062978, 0.09664277191720234, -0.07774831530062978, 0.0, 0.1556845390818705,
                 -0.09664277191720234, -0.1556845390818705, 0.0),
    Wiktorsson = c(0.0, 0.03652388482828246, 0.18205520315913828, -0.03652388482828246, 0.0, 0.033426332548990995,
                   -0.18205520315913828, -0.033426332548990995, 0.0),
    MronRoe = c(0.0, -0.0016970708928252332, -0.08600474949443807, 0.0016970708928252332, 0.0, 0.02862853807047939,
                0.08600474949443807, -0.02862853807047939, 0.0))
  lens <- c(Fourier = 30, Milstein = 33, Wiktorsson = 33, MronRoe = 36)
  for (alg in list(Fourier(), Milstein(), Wiktorsson(), MronRoe())) {
    A <- levyarea(W, 5, alg, cf(3, 5, alg))
    expect_lt(max(abs(as.numeric(A) - ref[[alg$name]])), 1e-15)
    expect_equal(coefficient_length(3, 5, alg), lens[[alg$name]])
    # the random form draws as generate_coefficients
    expect_equal(levyarea(W, 5, alg, rng = zgen()), A)
  }
  expect_lt(max(abs(trigamma(2:9) - c(0.6449340668482265, 0.39493406684822646, 0.2838229557371153,
                                       0.22132295573711533, 0.18132295573711535, 0.15354517795933756,
                                       0.13313701469403144, 0.11751201469403143))), 1e-15)
  expect_s3_class(MronRoe(), AbstractIteratedIntegralAlgorithm)
  expect_s3_class(MaxL2(), AbstractErrorNorm)
})

test_that("iterated integrals, terms and the optimal algorithm match", {
  h <- 0.04
  dW <- W * sqrt(h)
  cM <- cf(3, 5, MronRoe())
  expect_lt(max(abs(as.numeric(iterated_integrals(dW, h, cM)) -
    c(-0.0182, -0.004267882835713009, 0.003159810020222478, -0.0041321171642869905, -0.010200000000000002,
      -0.014254858477180825, 0.010040189979777524, -0.016545141522819176, 0.004200000000000006))), 1e-16)
  expect_lt(max(abs(as.numeric(iterated_integrals(dW, h, cM, ito_correction = FALSE)) -
    c(0.0018, -0.004267882835713009, 0.003159810020222478, -0.0041321171642869905, 0.009799999999999998,
      -0.014254858477180825, 0.010040189979777524, -0.016545141522819176, 0.024200000000000006))), 1e-16)
  expect_lt(max(abs(as.numeric(iterated_integrals(dW, h, cf(3, 5, Fourier()), alg = Fourier())) -
    c(-0.0182, -0.0013124237240757157, 0.0107034021325201, -0.007087576275924283, -0.010200000000000002,
      -0.008911924801296534, 0.002496597867479902, -0.021888075198703467, 0.004200000000000006))), 1e-16)
  expect_equal(c(iterated_integrals(0.3, 0.04), iterated_integrals(0.3, 0.04, ito_correction = FALSE),
                 iterated_integrals(0.3, 0.04, 0.01, q_12 = 0.5)), c(0.024999999999999998, 0.045, 0.04))
  cases <- list(list(2, 0.1, c(2, 1, 1, 1), c(4, 2, 2, 1), "Wiktorsson"),
                list(3, 0.01, c(16, 6, 4, 2), c(92, 31, 9, 4), "MronRoe"),
                list(5, 0.001, c(152, 51, 15, 7), c(3040, 1014, 65, 30), "MronRoe"),
                list(10, 0.5, c(1, 1, 1, 1), c(28, 10, 9, 4), "Fourier"),
                list(4, 1e-4, c(1520, 507, 42, 19), c(18238, 6080, 143, 64), "MronRoe"))
  algs <- list(Fourier(), Milstein(), Wiktorsson(), MronRoe())
  for (cs in cases) {
    e <- cs[[2]]^(3 / 2)
    expect_equal(vapply(algs, function(a) terms_needed(cs[[1]], cs[[2]], e, a, MaxL2()), 0L), cs[[3]])
    expect_equal(vapply(algs, function(a) terms_needed(cs[[1]], cs[[2]], e, a, FrobeniusL2()), 0L), cs[[4]])
    expect_equal(optimal_algorithm(cs[[1]], cs[[2]])$name, cs[[5]])
  }
  q <- c(1, 0.5, 0.25)
  expect_equal(vapply(algs, function(a) terms_needed(3, 0.01, 0.001, a, FrobeniusL2(), q_12 = q), 0L), c(10, 4, 3, 2))
  expect_equal(vapply(algs, function(a) terms_needed(3, 0.01, 0.001, a, MaxL2(), q_12 = q), 0L), c(4, 2, 2, 1))
  expect_equal(optimal_algorithm(3, 0.01, 0.001, q_12 = q)$name, "MronRoe")
  expect_error(terms_needed(2, 0.01, 0.001, MronRoe(), MaxL2(), q_12 = q), "Length of q_12")
})

test_that("the Q-Wiener iterated integrals use the corrected Levy area call", {
  # upstream passes the algorithm as a keyword levyarea does not take (a
  # MethodError); the reference is the intended formula assembled in
  # tools/levyarea-probe.jl from the same coefficients
  h <- 0.04
  dW <- W * sqrt(h)
  q <- c(1, 0.5, 0.25)
  I <- iterated_integrals(dW, h, 0.01, alg = MronRoe(), q_12 = q, rng = zgen())
  expect_lt(max(abs(as.numeric(I) -
    c(-0.0182, -0.005630649862205798, 0.007422868809449023, -0.002769350137794202, 0.004799999999999998,
      -0.011313250250542955, 0.005777131190550979, -0.019486749749457046, 0.022950000000000005))), 1e-16)
})

test_that("path reconstruction and sub-interval integrals match", {
  h <- 0.04
  dW <- W * sqrt(h)
  pts <- c(0, 0.01, 0.025, 0.04)
  expect_lt(max(abs(unlist(reconstruct_path(dW, h, cf(3, 5, MronRoe()), pts)) -
    c(0.0, 0.0, 0.0, -0.008157664300549601, 0.05426465644347476, -0.03142733490568743, 0.02010734315274284,
      -0.09538414100069605, 0.0901481213819783, 0.06, -0.13999999999999999, 0.22000000000000003))), 1e-15)
  expect_lt(max(abs(unlist(reconstruct_path(dW, h, cf(3, 5, Fourier()), pts)) -
    c(0.0, 0.0, 0.0, -0.005579078661121178, -0.05254980659160564, -0.002386265861507224, 0.0037462106105000358,
      -0.16331576972782558, 0.042821016564656034, 0.05999999999999999, -0.13999999999999999,
      0.22000000000000003))), 1e-15)
  expect_lt(max(abs(as.numeric(iterated_integrals_subinterval(dW, h, cf(3, 5, MronRoe()), 0.01, 0.03)) -
    c(-0.009408645237493934, -0.004033665876869242, 0.0033979822856464352, -0.000425036646560765,
      -0.0017712140869020152, -0.006144202269194235, 0.0023055615337005185, -0.015022126096683387,
      0.0032389590289477396))), 1e-14)
  expect_lt(max(abs(as.numeric(iterated_integrals_subinterval(dW, h, cf(3, 5, MronRoe()), 0, 0.04,
                                                              n_quadrature = 16, ito_correction = FALSE)) -
    c(0.001101657511355336, -0.004868786863777209, 0.00783729963504859, -0.0012136563330785953,
      0.0057093451990012254, -0.013824470238381181, 0.004227956976430769, -0.012678134884464565,
      0.021204054673387263))), 1e-14)
})

test_that("the iterated integral evaluators match StochasticDiffEqCore", {
  d2 <- c(0.1, -0.2)
  expect_equal(get_iterated_I(0.01, d2, NULL, JDiagonal_oop()), c(0.005000000000000001, 0.020000000000000004))
  expect_equal(as.numeric(get_iterated_I(0.01, d2, NULL, JCommute_oop())),
               c(0.005000000000000001, -0.010000000000000002, -0.010000000000000002, 0.020000000000000004))
  ji <- JDiagonal_iip(d2)
  get_iterated_I(0.01, d2, NULL, ji)
  expect_equal(ji$J, c(0.005000000000000001, 0.020000000000000004))
  ref <- list(Fourier = c(0.005000000000000001, -0.01255015159518915, -0.0074498484048108535, 0.020000000000000004),
              MronRoe = c(0.005000000000000001, -0.004895249775839353, -0.015104750224160652, 0.020000000000000004))
  for (alg in list(Fourier(), MronRoe())) {
    J <- get_iterated_I(0.01, d2, NULL, alg, rng = zgen())
    expect_lt(max(abs(as.numeric(J) - ref[[alg$name]])), 1e-16)
    it <- IteratedIntegralAlgorithm_iip(d2, alg)
    get_iterated_I(0.01, d2, NULL, it, rng = zgen())
    expect_equal(it$J, J)
  }
  diagp <- SDEProblem(function(u, p, t) u, function(u, p, t) u, c(1, 1), c(0, 1))
  genp <- SDEProblem(function(u, p, t) u, function(u, p, t) matrix(u, 2, 3), c(1, 1), c(0, 1),
                     noise_rate_prototype = matrix(0, 2, 3))
  expect_s3_class(get_Jalg(c(0.1, 0.2), 0.01, diagp, list(ii_approx = IILevyArea())), "JDiagonal_oop")
  expect_s3_class(get_Jalg(c(0.1, 0.2, 0.3), 0.01, genp, list(ii_approx = IILevyArea())), "MronRoe")
  expect_s3_class(get_Jalg(c(0.1, 0.2, 0.3), 0.01, genp, list(ii_approx = IICommutative())), "JCommute_oop")
  expect_s3_class(get_Jalg(c(0.1, 0.2, 0.3), 0.01, genp, list(ii_approx = Fourier())), "Fourier")
  expect_s3_class(IICommutative(), IteratedIntegralApprox)
})
