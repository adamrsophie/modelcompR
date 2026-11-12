# Stochastic differential equations: lowering of brownians, transformations,
# SDE solvers (strong orders against exact solutions) and ensembles.

test_that("brownians are lowered to diagonal, scalar and general noise", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t); z <- Variable("z", t)
  a <- Parameter("a"); b <- Parameter("b")
  Bx <- Brownian("stBx"); By <- Brownian("stBy"); Bz <- Brownian("stBz")
  s1 <- mtkcompile(System(list(D(x) ~ -x + a * x * Bx, D(y) ~ -y + b * By, D(z) ~ x - z + Bz), t, name = "d"))
  ne <- get_noise_eqs(s1)
  expect_false(is.matrix(ne))
  expect_length(ne, 3L)
  # equations, unknowns and noise rows share one order
  us <- vapply(unknowns(s1), format, "")
  drift <- stats::setNames(vapply(equations(s1), function(e) format(e$rhs), ""), us)
  noise <- stats::setNames(vapply(ne, format, ""), us)
  expect_equal(unname(drift[c("x(t)", "y(t)", "z(t)")]), c("-x(t)", "-y(t)", "x(t) - z(t)"))
  expect_equal(unname(noise[c("x(t)", "y(t)", "z(t)")]), c("a*x(t)", "b", "1"))
  s2 <- mtkcompile(System(list(D(x) ~ -x + a * x * Bx, D(y) ~ -y + b * Bx), t, name = "sc"))
  expect_equal(dim(get_noise_eqs(s2)), c(2L, 1L))
  # two equations share one brownian: the noise stays a matrix (correlated)
  s4 <- mtkcompile(System(list(D(x) ~ -x + a * Bx, D(y) ~ -y + b * Bx, D(z) ~ -z + By), t, name = "g2"))
  expect_equal(dim(get_noise_eqs(s4)), c(3L, 2L))
  prob <- SDEProblem(s4, list(x ~ 1, y ~ 2, z ~ 3, a ~ 2, b ~ 3), c(0, 1))
  expect_equal(prob$noise_type, "general")
  G <- prob$f$g(prob$u0, prob$p, 0)
  rownames(G) <- names(prob$u0)
  expect_equal(unname(G[c("x(t)", "y(t)", "z(t)"), ]), matrix(c(2, 3, 0, 0, 0, 1), 3, 2))
  expect_error(mtkcompile(System(list(D(x) ~ -x + Bx^2), t, name = "nl")), "linearly")
})

test_that("noise equations, Ito-Stratonovich and Girsanov transforms", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); a <- Parameter("a"); b <- Parameter("b")
  sde <- SDESystem(list(D(x) ~ a * x), list(b * x), t, name = "gbm")
  nb <- noise_to_brownians(sde, names = "w")
  expect_equal(format(equations(nb)[[1]]$rhs), format(a * x + b * x * Brownian("w_1")))
  cs <- mtkcompile(sde)
  expect_equal(format(equations(cs)[[1]]$rhs), "a*x(t)")
  expect_equal(format(get_noise_eqs(cs)[[1, 1]]), "b*x(t)")
  st <- stochastic_integral_transform(sde, -1/2)
  expect_equal(format(symcompR::simplify(equations(st)[[1]]$rhs - (a * x - b^2 * x / 2))), "0")
  back <- stochastic_integral_transform(st, 1/2)
  expect_equal(format(symcompR::simplify(equations(back)[[1]]$rhs - a * x)), "0")
  gm <- Girsanov_transform(sde, x, theta0 = 0.1)
  expect_length(unknowns(gm), 2L)
  expect_equal(format(symcompR::simplify(equations(gm)[[1]]$rhs - (a + b^2) * x)), "0")
  expect_equal(dim(get_noise_eqs(gm)), c(2L, 1L))
  expect_equal(format(observed(gm)[[1]]$lhs), "weight(t)")
})

test_that("geometric Brownian motion has ensemble mean exp(mu t)", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t)
  mu <- Parameter("mu", default = 0.5); sigma <- Parameter("sigma", default = 0.3)
  B <- Brownian("stB")
  sys <- mtkcompile(System(list(D(x) ~ mu * x + sigma * x * B), t, name = "gbm"))
  prob <- SDEProblem(sys, list(x ~ 1), c(0, 1))
  expect_s3_class(prob, "SDEProblem")
  sim <- solve(EnsembleProblem(prob), EM(), trajectories = 200, dt = 0.02, seed = 42)
  expect_length(sim, 200L)
  summ <- EnsembleSummary(sim, ts = c(0, 0.5, 1))
  se <- sqrt(summ$v[, 1] / 200)
  expect_equal(unname(summ$u[1, 1]), 1)
  expect_lt(abs(summ$u[3, 1] - exp(0.5)), 4 * se[3])
  expect_lt(abs(summ$u[2, 1] - exp(0.25)), 4 * se[2])
  # exact variance exp(2 mu t) (exp(sigma^2 t) - 1) within a wide band
  expect_lt(abs(summ$v[3, 1] / (exp(1) * (exp(0.09) - 1)) - 1), 0.35)
})

test_that("strong convergence: EM order 1/2, RKMil 1, SRIW1 and SOSRI 3/2", {
  a <- 1; b <- 0.5; N <- 200
  prob <- SDEProblem(function(u, p, t) a * u, function(u, p, t) b * u, rep(1, N), c(0, 1))
  dts <- 2^-(3:7)
  slope <- function(alg) {
    errs <- vapply(seq_along(dts), function(i) {
      sol <- solve(prob, alg, dt = dts[i], adaptive = FALSE, seed = 100 + i)
      W <- sol$W$W
      exact <- exp((a - b^2 / 2) + b * W[nrow(W), ])
      mean(abs(sol$u[[length(sol$t)]] - exact))
    }, 0)
    unname(coef(lm(log(errs) ~ log(dts)))[2])
  }
  s_em <- slope(EM())
  expect_gt(s_em, 0.3); expect_lt(s_em, 0.8)
  expect_gt(slope(RKMil()), 0.8)
  s_sri <- slope(SRIW1())
  expect_gt(s_sri, 1.25); expect_lt(s_sri, 1.8)
  expect_gt(slope(SOSRI()), 1.25)
})

test_that("adaptive SRIW1 with rejections follows the exact path", {
  # dX = X/2 dt + sqrt(1 + X^2) dW, X(0) = 0 has the solution sinh(W)
  prob <- SDEProblem(function(u, p, t) u / 2, function(u, p, t) sqrt(1 + u^2), 0, c(0, 1))
  n <- 150
  set.seed(77)
  res <- t(vapply(seq_len(n), function(i) {
    sol <- solve(prob, SRIW1(), abstol = 1e-3, reltol = 1e-3)
    W <- sol$W$W
    c(W[nrow(W), 1], sol$u[[length(sol$t)]], sol$stats$nreject, all.equal(sol$W$t, sol$t) == TRUE)
  }, numeric(4)))
  expect_gt(sum(res[, 3]), 20)
  expect_true(all(res[, 4] == 1))
  expect_lt(max(abs(res[, 2] - sinh(res[, 1]))), 0.02)
  # rejection sampling with memory keeps W(1) standard normal
  expect_lt(abs(mean(res[, 1])), 4 / sqrt(n))
  expect_lt(abs(var(res[, 1]) - 1), 4 * sqrt(2 / n))
})

test_that("raw SDE problems, noise shapes, seeds and solver checks", {
  f <- function(u, p, t) -u
  g <- function(u, p, t) matrix(c(0.1, 0, 0.05, 0.2), 2, 2)
  prob <- SDEProblem(f, g, c(1, 1), c(0, 1), noise_rate_prototype = matrix(0, 2, 2))
  s1 <- solve(prob, LambaEM(), seed = 3)
  s2 <- solve(prob, LambaEM(), seed = 3)
  expect_equal(s1$u, s2$u)
  expect_equal(ncol(s1$W$W), 2L)
  expect_error(solve(prob, SRIW1()), "diagonal or scalar")
  expect_error(solve(prob, EM()), "dt")
  set.seed(1); r1 <- runif(1)
  set.seed(1); invisible(solve(prob, EM(), dt = 0.1, seed = 99)); r2 <- runif(1)
  expect_equal(r1, r2)
  sp <- SDEProblem(function(u, p, t) 0 * u, function(u, p, t) c(1, 2), c(0, 0), c(0, 1), noise = "scalar")
  s3 <- solve(sp, EM(), dt = 0.1, seed = 5)
  expect_equal(2 * s3[1], s3[2])
  expect_s3_class(s3, "RODESolution")
})

test_that("symbolic indexing and observed variables of SDE solutions", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); w <- Variable("w", t); a <- Parameter("a", default = 0.1)
  B <- Brownian("stB2")
  sys <- mtkcompile(System(list(D(x) ~ -w + a * B, w ~ 2 * x), t, name = "o"))
  expect_length(observed(sys), 1L)
  sol <- solve(SDEProblem(sys, list(x ~ 1), c(0, 1)), SOSRI(), seed = 1)
  expect_equal(sol[w], 2 * sol[x])
  expect_equal(sol$retcode, "Success")
})

test_that("noise equations of subsystems are dropped when flattening, as in ModelingToolkit", {
  # ModelingToolkit's flatten() builds a new System without noise_eqs, and
  # mtkcompile() calls noise_to_brownians on the unflattened system
  t <- t_nounits; D <- D_nounits
  x <- Variable("nsx", t); y <- Variable("nsy", t); a <- Parameter("nsa"); b <- Parameter("nsb")
  sub <- System(list(D(x) ~ -a * x), t, list(x), list(a, b), noise_eqs = list(b * x), name = "nssub")
  root <- System(list(D(y) ~ -y), t, list(y), list(), systems = list(sub), name = "nsroot")
  expect_null(get_noise_eqs(flatten(root)))
  expect_length(equations(flatten(root)), 2L)
  cs <- mtkcompile(root)
  expect_null(get_noise_eqs(cs))
  expect_length(equations(cs), 2L)
  root2 <- System(list(D(y) ~ -y), t, list(y), list(b), noise_eqs = list(b * y), systems = list(sub),
                  name = "nsroot2")
  expect_length(get_noise_eqs(root2), 1L)
  expect_null(get_noise_eqs(flatten(root2)))
  expect_null(get_noise_eqs(complete(root2)))
  expect_error(mtkcompile(root2), "The system must be flattened.")
  expect_error(SDEProblem(complete(root2), list(y ~ 1, x ~ 1, a ~ 1, b ~ 0.1), c(0, 1)), "without noise")
  # a system without subsystems keeps its noise
  expect_length(get_noise_eqs(flatten(sub)), 1L)
})

test_that("implicit methods solve SDEs with a singular mass matrix", {
  # x' = -x + z, 0 = z - 2 x with noise s x dW: x is a geometric Brownian
  # motion with drift 1, x(t) = exp((1 - s^2 / 2) t + s W(t)), and z = 2 x
  s <- 0.5; N <- 50
  odd <- seq(1, 2 * N, by = 2)
  f <- function(u, p, t) {
    out <- numeric(2 * N)
    out[odd] <- -u[odd] + u[odd + 1]
    out[odd + 1] <- u[odd + 1] - 2 * u[odd]
    out
  }
  g <- function(u, p, t) {
    out <- numeric(2 * N)
    out[odd] <- s * u[odd]
    out
  }
  M <- diag(rep(c(1, 0), N))
  prob <- SDEProblem(SDEFunction(f, g, mass_matrix = M), rep(c(1, 2), N), c(0, 1))
  dts <- 2^-(3:6)
  slope <- function(alg) {
    errs <- vapply(seq_along(dts), function(i) {
      sol <- solve(prob, alg, dt = dts[i], adaptive = FALSE, seed = 300 + i)
      W <- sol$W$W
      uend <- sol$u[[length(sol$t)]]
      # the algebraic equation holds at the end of every step
      expect_lt(max(abs(uend[odd + 1] - 2 * uend[odd])), 1e-6)
      exact <- exp((1 - s^2 / 2) + s * W[nrow(W), odd])
      mean(abs(uend[odd] - exact))
    }, 0)
    unname(coef(lm(log(errs) ~ log(dts)))[2])
  }
  s_em <- slope(ImplicitEM())
  expect_gt(s_em, 0.25); expect_lt(s_em, 0.85)
  expect_gt(slope(ImplicitRKMil()), 0.75)
  expect_gt(slope(SImplicitMidpoint()), 0.25)
  # StochasticDiffEq's compatibility rules
  expect_error(solve(prob, STrapezoid(), dt = 0.1, adaptive = FALSE), "symplectic or theta=1")
  expect_error(solve(prob, EM(), dt = 0.1), "not able to use mass matrices")
  sol <- solve(prob, seed = 5)
  expect_s3_class(sol$alg, "ImplicitRKMil")
  expect_equal(sol$retcode, "Success")
})

test_that("SDE systems with algebraic equations get a mass matrix", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("mmsx", t); z <- Variable("mmsz", t); B <- Brownian("mmsB")
  sys <- mtkcompile(System(list(D(x) ~ -x + z + 0.1 * x * B, 0 ~ z^3 + z - x), t, name = "mms"))
  prob <- SDEProblem(sys, list(x ~ 1), c(0, 1), guesses = list(z ~ 0.5))
  expect_equal(unname(prob$f$mass_matrix), diag(c(1, 0)))
  z0 <- prob[z]
  expect_equal(z0^3 + z0, 1, tolerance = 1e-8)
  sol <- solve(prob, seed = 1, abstol = 1e-4, reltol = 1e-4)
  expect_equal(sol$retcode, "Success")
  # StochasticDiffEq chooses ISSEM for a mass matrix with non diagonal noise
  expect_s3_class(sol$alg, "ISSEM")
  # ISSEM adds the noise after its implicit drift step, so x moves by
  # g dW without z following: the constraint holds to the size of the noise
  zx <- sol[z]; xx <- sol[x]
  expect_lt(max(abs(zx^3 + zx - xx)), 0.05)
  # a fully implicit method solves the constraint at the end of each step;
  # the Newton iterations stop at the tolerance of the step
  sol1 <- solve(prob, ImplicitEM(), seed = 1, abstol = 1e-4, reltol = 1e-4)
  zx <- sol1[z]; xx <- sol1[x]
  expect_lt(max(abs(zx^3 + zx - xx)), 1e-3)
  sol2 <- solve(prob, ImplicitRKMil(), seed = 1)
  expect_equal(sol2$retcode, "Success")
})

test_that("split step, additive noise and commutative noise methods converge", {
  # geometric Brownian motion against its exact solution on the same path
  a <- 1.01; b <- 0.87
  gbm <- SDEProblem(function(u, p, t) a * u, function(u, p, t) b * u, 0.5, c(0, 1))
  err <- function(alg, dt) mean(vapply(1:6, function(i) {
    s <- solve(gbm, alg, dt = dt, adaptive = FALSE, seed = i)
    W <- tail(s$W$W, 1)[1]
    strat <- alg$interpretation == "Stratonovich"
    abs(tail(s$u, 1)[[1]] - 0.5 * exp(if (strat) a + b * W else a - b^2 / 2 + b * W))
  }, 0))
  for (alg in list(ISSEM(), ISSEulerHeun(), LambaEulerHeun(), RKMilCommute())) {
    expect_lt(err(alg, 1 / 128), err(alg, 1 / 16), label = alg$name)
    expect_equal(solve(gbm, alg, seed = 1)$retcode, "Success", label = alg$name)
  }
  # additive noise: u' = cos(t) + 0.3 W', exact 1 + sin(t) + 0.3 W(t)
  add <- SDEProblem(function(u, p, t) cos(t) + 0 * u, function(u, p, t) 0.3 + 0 * u, 1, c(0, 1))
  for (alg in list(SRA1(), SRA2(), SRA3(), SOSRA(), SOSRA2(), SKenCarp())) {
    s <- solve(add, alg, dt = 1 / 16, adaptive = FALSE, seed = 2)
    expect_equal(tail(s$u, 1)[[1]], 1 + sin(1) + 0.3 * tail(s$W$W, 1)[1], tolerance = 1e-3, label = alg$name)
  }
  # Ornstein-Uhlenbeck moments with the adaptive additive noise methods
  ou <- SDEProblem(function(u, p, t) -2 * u, function(u, p, t) 0.5 + 0 * u, 1, c(0, 1))
  v <- vapply(1:150, function(i) tail(solve(ou, SOSRA(), seed = i)$u, 1)[[1]], 0)
  expect_lt(abs(mean(v) - exp(-2)), 4 * sqrt(0.0614 / 150))
  # the default algorithm follows StochasticDiffEq's alg_hints
  expect_s3_class(.sde_default_alg(ou, "additive"), "SOSRA")
  expect_s3_class(.sde_default_alg(ou, c("additive", "stiff")), "SKenCarp")
  expect_s3_class(.sde_default_alg(ou, "commutative"), "RKMilCommute")
  expect_equal(.sde_default_alg(ou, "stratonovich")$interpretation, "Stratonovich")
  gen <- ou
  gen$noise_type <- "general"
  expect_s3_class(.sde_default_alg(gen, "stiff"), "ISSEM")
  expect_s3_class(.sde_default_alg(gen, "stratonovich"), "LambaEulerHeun")
  expect_s3_class(.sde_default_alg(gen, c("stratonovich", "stiff")), "ImplicitEulerHeun")
})
