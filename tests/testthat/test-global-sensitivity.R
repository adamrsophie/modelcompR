# QuasiMonteCarlo.jl and GlobalSensitivity.jl against the Julia references of
# tools/globalsensitivity-probe.jl (raw Float64, read in the probe's order).

ref_con <- function() {
  f <- test_path("fixtures", "julia", "globalsensitivity.bin")
  skip_if_not(file.exists(f))
  file(f, "rb")
}

ishi_gsa <- function(X) sin(X[1]) + 7 * sin(X[2])^2 + 0.1 * X[3]^4 * sin(X[1])

test_that("QuasiMonteCarlo samplers, ranges and the analyses match Julia", {
  con <- ref_con()
  on.exit(close(con))
  R <- function(n) readBin(con, "double", n)
  same <- function(got) {
    got <- as.numeric(unlist(got))
    expect_identical(got, R(length(got)))
  }
  near <- function(got, tol = 1e-13) {
    got <- as.numeric(unlist(got))
    expect_equal(got, R(length(got)), tolerance = tol)
  }
  same(qmc_sample(37, 5, SobolSample()))
  same(qmc_sample(50, 4, HaltonSample()))
  same(qmc_sample(50, 5, FaureSample()))
  same(qmc_sample(18, 3, FaureSample()))
  same(qmc_sample(40, 6, KroneckerSample()))
  same(qmc_sample(33, 7, LatticeRuleSample()))
  same(qmc_sample(5, 300, LatticeRuleSample()))
  same(qmc_sample(30, 1, VanDerCorputSample(3)))
  same(qmc_sample(125, 1, VanDerCorputSample(5)))
  same(qmc_sample(10, c(0, -1, 2), c(1, 3, 2.5), HaltonSample()))
  same(.qmc_unif2bits(t(qmc_sample(27, 3, FaureSample())), 3, 20))
  expect_warning(dm <- generate_design_matrices(16, c(0, 0), c(1, 2), SobolSample(), 3), "independent")
  same(do.call(cbind, dm))
  for (cs in list(c(0, 10, 7), c(-pi, pi, 1000), c(0.1, 0.7, 11), c(-1.2345, 3.3, 2048),
                  c(2.5e-3, 17.25, 33), c(-0.37, 1e5 + 0.1, 101))) {
    r <- .jl_linrange(cs[1], cs[2], cs[3])
    same(as.numeric(r))
    same(attr(r, "step"))
  }

  lb <- rep(-pi, 3); ub <- rep(pi, 3)
  AB <- suppressWarnings(generate_design_matrices(512, lb, ub, SobolSample(), 2))
  for (est in c("Jansen1999", "Homma1996", "Sobol2007", "Janon2014")) {
    r <- gsa(ishi_gsa, Sobol(order = 0:2), AB[[1]], AB[[2]], Ei_estimator = est)
    near(r$S1); near(r$ST); near(r$S2)
  }
  r <- gsa(ishi_gsa, Sobol(order = 0:2, nboot = 4), AB[[1]], AB[[2]])
  near(r$S1); near(r$ST); near(r$S2); near(r$S1_Conf_Int); near(r$ST_Conf_Int); near(r$S2_Conf_Int)
  r <- gsa(function(X) c(ishi_gsa(X), X[1] * X[2] + X[3]), Sobol(order = 0:2), AB[[1]], AB[[2]])
  expect_equal(dim(r$S1), c(2L, 3L))
  near(r$S1); near(r$ST); near(r$S2)

  X <- matrix(R(1800), 3); Y <- R(600)
  expect_identical(X, qmc_sample(600, lb, ub, SobolSample()))
  r <- gsa(X, matrix(Y, 1), RegressionGSA(rank = TRUE))
  near(r$pearson); near(r$standard_regression); near(r$partial_correlation)
  r <- gsa(X, Y, EASI()); near(r$S1); near(r$S1_Corr)
  r <- gsa(X, Y, EASI(dct_method = TRUE)); near(r$S1); near(r$S1_Corr)
  grid <- as.numeric(.jl_linrange(min(Y), max(Y), 2048))
  near(.kde_pdf(Y, grid)); near(.kde_pdf(Y[1:137], grid))
  cuts <- as.numeric(.jl_linrange(0, 600, 21))
  near(vapply(1:3, function(i) .gsa_calc_delta(X[i, ], Y, grid, cuts), 0))
  nb <- round(sqrt(600))
  same(vapply(1:3, function(i) .gsa_entropy(X[i, ], nb), 0))
  same(.gsa_entropy(Y, nb))
  same(vapply(1:3, function(i) .gsa_entropy(cbind(X[i, ], Y), nb), 0))
  Xd <- rbind(X, qmc_sample(600, c(0, 0), c(1, 1), SobolSample()))
  r <- .gsa_rsa_analysis(Xd, Y, RSA(n_dummy_parameters = 2))
  same(r$S); same(r$Sd[1]); same(Xd[4:5, ])
  s <- lapply(1:3, function(i) R(301)); Yb <- R(301)
  near(.gsa_rbdfast_analysis(Yb, s, 6))
  us <- R(3); k <- 0
  rng <- function(n) { k <<- k + n; us[(k - n + 1):k] }
  r <- gsa(ishi_gsa, eFAST(), rep(list(c(-pi, pi)), 3), samples = 400, rng = rng)
  near(r$S1, 1e-12); near(r$ST, 1e-12)
  ff <- function(x) 2 * x[1] - 3 * x[2] + 0.5 * x[3] * x[1] + x[4]
  r <- gsa(ff, FractionalFactorial(), num_params = 3, p_range = list(c(0, 1), c(-1, 2), c(1, 3)))
  same(r[[1]]); same(r[[2]])
  expect_length(readBin(con, "double", 10), 0L)
})

test_that("the random GSA methods recover known sensitivities", {
  set.seed(1)
  rg <- rep(list(c(-pi, pi)), 3)
  # Ishigami: S1 = (0.3139, 0.4424, 0), ST = (0.5576, 0.4424, 0.2437)
  r <- suppressWarnings(gsa(ishi_gsa, Sobol(), rg, samples = 20000))
  expect_equal(r$S1, c(0.3139, 0.4424, 0), tolerance = 0.02)
  expect_equal(r$ST, c(0.5576, 0.4424, 0.2437), tolerance = 0.02)
  expect_equal(gsa(ishi_gsa, eFAST(), rg, samples = 2000)$S1[1, ], c(0.3139, 0.4424, 0), tolerance = 0.02)
  expect_equal(gsa(ishi_gsa, EASI(), rg, samples = 4000)$S1, c(0.3139, 0.4424, 0), tolerance = 0.03)
  m <- gsa(function(x) 2 * x[1] - 3 * x[2] + 0.5 * x[3], Morris(num_trajectory = 10), rep(list(c(0, 1)), 3))
  expect_equal(as.numeric(m$means), c(2, -3, 0.5))
  expect_equal(as.numeric(m$means_star), c(2, 3, 0.5))
  dg <- gsa(function(x) x[1] + 2 * x[2]^2, DGSM(crossed = TRUE),
            rep(list(univariate_distribution("unif", 0, 1)), 2), samples = 200)
  expect_equal(dg$a[1], 1)
  expect_equal(dg$crossed, matrix(0, 2, 2))
  Di <- SklarDist(IndependentCopula(3), rep(list(univariate_distribution("unif", -pi, pi)), 3))
  sh <- suppressMessages(gsa(ishi_gsa, Shapley(n_var = 2000, n_outer = 300, n_inner = 3), Di))
  expect_equal(sum(sh$shapley_effects), 1, tolerance = 0.1)
  D <- SklarDist(GaussianCopula(matrix(c(1, 0.5, 0.5, 1), 2)), rep(list(univariate_distribution("norm")), 2))
  sh <- suppressMessages(gsa(function(x) x[1] + x[2], Shapley(n_var = 2000, n_outer = 200), D))
  expect_equal(sh$shapley_effects, c(0.5, 0.5), tolerance = 0.1)
  rs <- gsa(ishi_gsa, RSA(), rg, samples = 1000)
  expect_true(all(rs$S[1:2] > rs$Sd[1]))
  mi <- gsa(ishi_gsa, MutualInformation(n_bootstraps = 50), rg, samples = 1000)
  expect_true(mi$S[2] > mi$S[3])
  dm <- gsa(ishi_gsa, DeltaMoment(nboot = 10), rg, samples = 600)
  expect_true(dm$deltas[2] > dm$deltas[3])
  set.seed(5)
  bat <- gsa(function(X) apply(X, 2, ishi_gsa), Sobol(), rg, samples = 512, batch = TRUE)
  set.seed(5)
  ser <- gsa(ishi_gsa, Sobol(), rg, samples = 512)
  expect_equal(bat$S1, ser$S1)
})

test_that("the estimators are unbiased on known indices", {
  set.seed(3)
  bat <- function(X) sin(X[1, ]) + 7 * sin(X[2, ])^2 + 0.1 * X[3, ]^4 * sin(X[1, ])
  rg <- rep(list(c(-pi, pi)), 3)
  # independent scrambled replicates: their mean is within a few standard
  # errors of the analytic indices
  R <- replicate(40, {
    r <- gsa(bat, Sobol(order = 0:2), rg, samples = 256, batch = TRUE)
    c(r$S1, r$ST, r$S2[1, 3])
  })
  z <- (rowMeans(R) - c(0.3139, 0.4424, 0, 0.5576, 0.4424, 0.2437, 0.2437)) / (apply(R, 1, stats::sd) / sqrt(40))
  expect_true(all(abs(z) < 4))
  r <- gsa(bat, Sobol(nboot = 10), rg, samples = 256, batch = TRUE)
  expect_true(all(abs(r$ST - c(0.5576, 0.4424, 0.2437)) < 4 * r$ST_Conf_Int / 1.96 + 0.01))
  # RBD-FAST with the Tissot and Prieur correction: a null input near 0
  lin <- function(x) x[1] + 2 * x[2] + 0 * x[3]
  B <- replicate(30, gsa(lin, RBDFAST(), num_params = 3, samples = 200))
  expect_equal(rowMeans(B), c(0.2, 0.8, 0), tolerance = 0.02)
  E <- replicate(30, gsa(function(x) x[1] + 2 * x[2], EASI(), rep(list(c(0, 1)), 3), samples = 1000)$S1_Corr)
  expect_equal(rowMeans(E), c(0.2, 0.8, 0), tolerance = 0.02)
  # DGSM tau and sigma on uniform inputs: equal to the total variance part
  # for linear terms, an upper bound otherwise, NA for other laws
  f <- function(x) x[1] + 2 * x[2]^2 + x[3]
  d <- list(univariate_distribution("unif", 0, 1), univariate_distribution("unif", 0, 1),
            univariate_distribution("unif", -1, 3))
  dg <- gsa(f, DGSM(), d, samples = 4000)
  expect_equal(dg$tao[c(1, 3)], c(1, 16) / 12, tolerance = 0.03)
  expect_equal(dg$sigma[c(1, 3)], c(1, 16) / 12, tolerance = 0.03)
  expect_gt(dg$sigma[2], 4 * (1 / 5 - 1 / 9))
  expect_true(is.na(gsa(f, DGSM(), list(univariate_distribution("norm"), d[[2]], d[[3]]), samples = 20)$tao[1]))
})

test_that("randomized QMC keeps points in the unit cube", {
  set.seed(2)
  x <- qmc_sample(16, 2, SobolSample())
  for (R in list(OwenScramble(2), MatousekScramble(2), DigitalShift(2), QMCShift())) {
    y <- randomize(x, R)
    expect_equal(dim(y), dim(x))
    expect_true(all(y >= 0 & y < 1))
  }
  expect_equal(lengths(generate_design_matrices(16, c(0, 0), c(1, 1), SobolSample(R = OwenScramble(2)), 3)), rep(32L, 3))
  expect_equal(dim(qmc_sample(10, 3, LatinHypercubeSample())), c(3L, 10L))
  expect_equal(dim(qmc_sample(10, 3, RandomizedHaltonSample())), c(3L, 10L))
})
