# Reference values from tools/manopt-probe.jl (Manopt 0.6.12, Manifolds 0.11.38).
# The R iterates equal the Julia ones bit for bit except where the exponential
# map calls cos and sin (a few ulps from the system math library); the bound
# 3e-16 also absorbs the last digit that R's parser can lose on 17 digit
# literals.

mo_f <- function(M, p) 2 * p[1]^2 + p[2]^2 + 3 * p[3]^2 + p[1] * p[2] + 0.2 * p[1] * p[3] + 0.6 * p[2] * p[3] +
  p[1] * p[2] * p[3]
mo_egrad <- function(p) c(4 * p[1] + p[2] + 0.2 * p[3] + p[2] * p[3],
                          2 * p[2] + p[1] + 0.6 * p[3] + p[1] * p[3],
                          6 * p[3] + 0.2 * p[1] + 0.6 * p[2] + p[1] * p[2])
mo_ehess <- function(p, v) c(4 * v[1] + (1 + p[3]) * v[2] + (0.2 + p[2]) * v[3],
                             (1 + p[3]) * v[1] + 2 * v[2] + (0.6 + p[1]) * v[3],
                             (0.2 + p[2]) * v[1] + (0.6 + p[1]) * v[2] + 6 * v[3])
mo_grad <- function(M, p) { G <- mo_egrad(p); G - sum(p * G) * p }
mo_grad_jl <- function(M, p) { G <- mo_egrad(p); G - modelcompR:::.jl_dot(p, G) * p }
mo_hess_jl <- function(M, p, X) {
  dot <- modelcompR:::.jl_dot
  H <- mo_ehess(p, X)
  Y <- H - dot(p, H) * p
  G <- mo_egrad(p)
  Y + (-dot(p, G - (G - dot(p, G) * p))) * X
}
mo_p0 <- c(0.6, -0.48, 0.64)
mo_dev <- function(a, b) max(abs(a - b))

test_that("Manopt solvers on the sphere follow Manopt.jl iterate by iterate", {
  M <- Sphere(2)
  S <- Manopt$StopAfterIteration
  expect_lt(mo_dev(Manopt$gradient_descent(M, mo_f, mo_grad_jl, mo_p0, stopping_criterion = S(25)),
                   c(0.39124950596404234, -0.90330559854597, 0.17596255203423183)), 3e-16)
  expect_lt(mo_dev(Manopt$gradient_descent(M, mo_f, mo_grad_jl, mo_p0, stopping_criterion = S(25),
                                           retraction_method = Manopt$ProjectionRetraction()),
                   c(0.3912478498277492, -0.9033046610653426, 0.17597104677414904)), 3e-16)
  expect_lt(mo_dev(Manopt$conjugate_gradient_descent(M, mo_f, mo_grad_jl, mo_p0, stopping_criterion = S(25)),
                   c(0.6713152328342332, -0.7124620055445724, 0.2042883961954886)), 3e-16)
  expect_lt(mo_dev(Manopt$quasi_Newton(M, mo_f, mo_grad_jl, mo_p0, stopping_criterion = S(12)),
                   c(0.39112903367280033, -0.9032372863251186, 0.17657996378473167)), 3e-16)
  expect_lt(mo_dev(Manopt$trust_regions(M, mo_f, mo_grad_jl, mo_hess_jl, mo_p0, stopping_criterion = S(10)),
                   c(0.39112903356505907, -0.9032372863824442, 0.17657996373015128)), 3e-16)
  expect_lt(mo_dev(Manopt$adaptive_regularization_with_cubics(M, mo_f, mo_grad_jl, mo_hess_jl, mo_p0,
                                                              stopping_criterion = S(8)),
                   c(0.39112903356505907, -0.9032372863824442, 0.1765799637301513)), 1e-15)
  expect_lt(mo_dev(Manopt$Frank_Wolfe_method(M, mo_f, mo_grad_jl, mo_p0, stopping_criterion = S(6)),
                   c(-0.5885778106963631, 0.39117979814477166, -0.7074987818217725)), 3e-16)
  simplex <- list(c(0.6, -0.48, 0.64), c(0, 0.6, 0.8), c(-0.48, 0.64, 0.6))
  expect_lt(mo_dev(Manopt$NelderMead(M, mo_f, simplex, stopping_criterion = S(30)),
                   c(0.3911386562926317, -0.9032507871367662, 0.17648956651996975)), 3e-16)

  # the default stopping criteria: iterations until the gradient norm is small
  st <- Manopt$gradient_descent(M, mo_f, mo_grad_jl, mo_p0, return_state = TRUE)
  expect_equal(st$iterations, 70L)
  expect_true(st$converged)
  expect_lt(mo_dev(Manopt$get_solver_result(st), c(0.39112903377215263, -0.9032372865000773, 0.1765799626697204)),
            1e-15)
  st <- Manopt$quasi_Newton(M, mo_f, mo_grad_jl, mo_p0, return_state = TRUE)
  expect_equal(st$iterations, 7L)
  expect_lt(mo_dev(st$p, c(0.3911290505705838, -0.9032372702955959, 0.1765800083495141)), 3e-16)
  st <- Manopt$trust_regions(M, mo_f, mo_grad_jl, mo_hess_jl, mo_p0, return_state = TRUE)
  expect_equal(st$iterations, 5L)
  expect_equal(st$stopped_by$type, "StopWhenGradientNormLess")
})

test_that("Manopt solvers on Euclidean(2) follow Manopt.jl on the Rosenbrock function", {
  R2 <- Euclidean(2)
  pp <- c(1, 100)
  rf <- function(M, x) (pp[1] - x[1])^2 + pp[2] * (x[2] - x[1]^2)^2
  rg <- function(M, x) c(-2.0 * (pp[1] - x[1]) - 4.0 * pp[2] * (x[2] - x[1]^2) * x[1], 2.0 * pp[2] * (x[2] - x[1]^2))
  rh <- function(M, x, v) c((2 - 4 * pp[2] * (x[2] - 3 * x[1]^2)) * v[1] + (-4 * pp[2] * x[1]) * v[2],
                            (-4 * pp[2] * x[1]) * v[1] + (2 * pp[2]) * v[2])
  S <- Manopt$StopAfterIteration
  expect_lt(mo_dev(Manopt$gradient_descent(R2, rf, rg, c(0, 0), stopping_criterion = S(50)),
                   c(0.47122696121008584, 0.22722515495628925)), 3e-16)
  expect_lt(mo_dev(Manopt$quasi_Newton(R2, rf, rg, c(0, 0), stopping_criterion = S(20)),
                   c(0.9914139549210506, 0.982768046285956)), 3e-16)
  # with its defaults on a two dimensional problem the cubic regularization
  # stops after one step: all max(2, 1) - 1 Lanczos vectors are used
  st <- Manopt$adaptive_regularization_with_cubics(R2, rf, rg, rh, c(0, 0), return_state = TRUE)
  expect_lt(mo_dev(st$p, c(0.1546307035639236, 0)), 3e-16)
  expect_equal(st$stopped_by$type, "StopWhenAllLanczosVectorsUsed")
  st <- Manopt$trust_regions(R2, rf, rg, rh, c(0, 0), return_state = TRUE)
  expect_equal(st$iterations, 18L)
  expect_equal(rf(R2, st$p), 3.639224896962477e-27, tolerance = 1e-12)
  st <- Manopt$quasi_Newton(R2, rf, rg, c(0, 0), return_state = TRUE,
                            stopping_criterion = Manopt$StopWhenAny(S(30), Manopt$StopWhenGradientNormLess(1e-6)))
  expect_equal(st$iterations, 26L)
  expect_lt(mo_dev(st$p, c(0.9999999903261467, 0.9999999796894145)), 3e-16)
})

test_that("the OptimizationManopt optimizers solve the Rosenbrock problem", {
  set.seed(42)
  rosen <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
  rgrad <- function(x, p) c(-2 * (p[1] - x[1]) - 4 * p[2] * (x[2] - x[1]^2) * x[1], 2 * p[2] * (x[2] - x[1]^2))
  rhess <- function(x, p) matrix(c(2 - 4 * p[2] * (x[2] - 3 * x[1]^2), -4 * p[2] * x[1], -4 * p[2] * x[1],
                                   2 * p[2]), 2)
  R2 <- Euclidean(2)
  expect_error(solve(OptimizationProblem(rosen, c(0, 0), c(1, 100)), GradientDescentOptimizer()),
               "Manifold not specified")
  cons <- OptimizationFunction(rosen, grad = rgrad, cons = function(x, p) c(x[1]^2 + x[2]^2, x[1] * x[2]))
  expect_error(solve(OptimizationProblem(cons, c(0, 0), c(1, 100), lcons = c(0, 0), ucons = c(1, 1)),
                     GradientDescentOptimizer()), class = "IncompatibleOptimizerError")
  expect_true(requiresgradient(QuasiNewtonOptimizer()))
  expect_true(requireshessian(TrustRegionsOptimizer()))
  expect_false(requiresgradient(NelderMeadOptimizer()))
  expect_false(allowsbounds(CMAESOptimizer()))

  of <- OptimizationFunction(rosen, grad = rgrad, hess = rhess)
  prob <- OptimizationProblem(of, c(0, 0), c(1, 100), manifold = R2)
  seen <- 0L
  sol <- solve(prob, QuasiNewtonOptimizer(), maxiters = 30, callback = function(state, l) {
    seen <<- seen + 1L
    FALSE
  })
  expect_lt(sol$objective, 1e-12)
  expect_equal(sol$retcode, "Success")
  expect_gt(seen, 0L)
  sol <- solve(prob, TrustRegionsOptimizer())
  expect_lt(sol$objective, 0.1)
  expect_equal(sol$retcode, "Success")
  # maxiters replaces the Manopt defaults, so the cubic regularization runs
  sol <- solve(prob, AdaptiveRegularizationCubicOptimizer(), maxiters = 100)
  expect_lt(sol$objective, 0.1)
  expect_equal(sol$retcode, "Success")
  # a user criterion is added and suppresses the default fallback
  sol <- solve(prob, GradientDescentOptimizer(), maxiters = 2000,
               stopping_criterion = Manopt$StopWhenGradientNormLess(1))
  expect_equal(sol$retcode, "Success")
  expect_lt(sol$original$iterations, 2000L)
  sol <- solve(prob, ConvexBundleOptimizer())
  expect_lt(sol$objective, 0.1)
  expect_lt(solve(prob, ConjugateGradientDescentOptimizer(), maxiters = 50)$objective, 0.5)

  # derivative free methods on the cost alone
  dfp <- OptimizationProblem(rosen, c(0, 0), c(1, 100), manifold = R2)
  expect_lt(solve(dfp, NelderMeadOptimizer())$objective, 0.7)
  expect_lt(solve(dfp, CMAESOptimizer())$objective, 0.1)
  # no convergence criterion is added for the swarm: maxiters ends the run
  sol <- solve(dfp, ParticleSwarmOptimizer(), population_size = 20, maxiters = 100)
  expect_lt(sol$objective, 0.1)
  expect_equal(sol$retcode, "MaxIters")
})

test_that("problems on the sphere go through the Riemannian gradient and Hessian", {
  M <- Sphere(2)
  f <- function(x, p) mo_f(M, x)
  prob <- OptimizationProblem(OptimizationFunction(f, grad = function(x, p) mo_egrad(x)), mo_p0, manifold = M)
  sol <- solve(prob, GradientDescentOptimizer(), stopping_criterion = Manopt$StopAfterIteration(25))
  expect_lt(mo_dev(sol$u, c(0.39124950596404234, -0.90330559854597, 0.17596255203423183)), 3e-16)
  expect_equal(sol$retcode, "MaxIters")
  # the Hessian from Hessian vector products of the Euclidean Hessian
  H <- function(x, p) matrix(c(4, 1 + x[3], 0.2 + x[2], 1 + x[3], 2, 0.6 + x[1], 0.2 + x[2], 0.6 + x[1], 6), 3)
  hp <- OptimizationProblem(OptimizationFunction(f, grad = function(x, p) mo_egrad(x), hess = H), mo_p0, manifold = M)
  sol <- solve(hp, TrustRegionsOptimizer(), maxiters = 50)
  expect_equal(sol$u, c(0.3911290335650558, -0.9032372863824442, 0.17657996373015825), tolerance = 1e-10)
  expect_equal(sol$retcode, "Success")
  # maximizing the negated cost lands on the same point; the objective keeps its sign
  mx <- OptimizationProblem(OptimizationFunction(function(x, p) -f(x, p), grad = function(x, p) -mo_egrad(x)),
                            mo_p0, manifold = M, sense = MaxSense)
  sol <- solve(mx, QuasiNewtonOptimizer(), maxiters = 100)
  expect_equal(sol$u, c(0.3911290335650558, -0.9032372863824442, 0.17657996373015825), tolerance = 1e-6)
  expect_lt(sol$objective, 0)
  expect_error(solve(OptimizationProblem(f, c(1, 0), manifold = M), GradientDescentOptimizer()), "have length 3")
  expect_output(print(M), "Sphere\\(2\\)")
  expect_output(print(Euclidean(3)), "Euclidean\\(3\\)")
})

test_that("SymmetricPositiveDefinite carries the affine invariant geometry", {
  M <- SymmetricPositiveDefinite(2)
  A <- matrix(c(2, 0.5, 0.5, 1), 2); B <- matrix(c(1, -0.3, -0.3, 3), 2)
  X <- matrix(c(0.1, 0.2, 0.2, -0.3), 2)
  mf <- function(name, ...) get(name, envir = asNamespace("modelcompR"))(M, ...)
  q <- mf(".mf_exp", as.numeric(A), as.numeric(X))
  expect_equal(mf(".mf_log", as.numeric(A), q), as.numeric(X), tolerance = 1e-12)
  Ai <- solve(A)
  expect_equal(mf(".mf_inner", as.numeric(A), X, X), sum(diag(Ai %*% X %*% Ai %*% X)))
  expect_equal(mf(".mf_distance", as.numeric(A), q), sqrt(sum(diag(Ai %*% X %*% Ai %*% X))))
  Y <- mf(".mf_transport", as.numeric(A), as.numeric(X), as.numeric(B))
  expect_equal(mf(".mf_norm", as.numeric(B), Y), mf(".mf_norm", as.numeric(A), as.numeric(X)))
  expect_equal(mf(".mf_get_vector", as.numeric(A), mf(".mf_get_coordinates", as.numeric(A), as.numeric(X))),
               as.numeric(X))
  G <- matrix(c(1, 3, 2, 4), 2)
  expect_equal(mf(".mf_rgrad", as.numeric(A), G), as.numeric(A %*% ((G + t(G)) / 2) %*% A))
  # the Karcher mean of two points is their geodesic midpoint
  mid <- matrix(mf(".mf_exp", as.numeric(A), 0.5 * mf(".mf_log", as.numeric(A), as.numeric(B))), 2)
  f <- function(x, p) distance(M, x, A)^2 + distance(M, x, B)^2
  for (alg in list(GradientDescentOptimizer(), QuasiNewtonOptimizer(), TrustRegionsOptimizer())) {
    sol <- solve(OptimizationProblem(f, diag(2), manifold = M), alg, maxiters = 200)
    expect_identical(sol$retcode, "Success")
    expect_true(is.matrix(sol$u))
    expect_equal(sol$u, mid, tolerance = 1e-6)
  }
  kg <- function(M, p) -2 * (matrix(mf(".mf_log", as.numeric(p), as.numeric(A)), 2) +
                               matrix(mf(".mf_log", as.numeric(p), as.numeric(B)), 2))
  p <- Manopt$gradient_descent(M, function(M, p) f(p, NULL), kg, diag(2))
  expect_equal(p, mid, tolerance = 1e-6)
})

test_that("SymmetricPositiveDefinite matches Manifolds.jl and Manopt.jl", {
  # values of tools/manopt-probe.jl; the maps go through eigen decompositions
  # whose LAPACK results differ from Julia's OpenBLAS by a few ulps
  M <- SymmetricPositiveDefinite(2)
  A <- matrix(c(2, 0.5, 0.5, 1), 2); B <- matrix(c(1, -0.3, -0.3, 3), 2)
  X <- c(0.1, 0.2, 0.2, -0.3)
  mf <- function(name, ...) get(name, envir = asNamespace("modelcompR"))(M, ...)
  a <- as.numeric(A); b <- as.numeric(B)
  expect_equal(mf(".mf_exp", a, X), c(2.1182432448980046, 0.6747567208287562, 0.6747567208287562, 0.7687298032404923),
               tolerance = 1e-14)
  expect_equal(mf(".mf_log", a, b), c(-1.543705951730838, -0.7418335371745315, -0.7418335371745314, 0.8459063388701346),
               tolerance = 1e-14)
  expect_equal(mf(".mf_inner", a, X, c(0.3, -0.1, -0.1, 0.2)), -0.1657142857142857, tolerance = 1e-14)
  expect_equal(mf(".mf_distance", a, b), 1.5426995089660311, tolerance = 1e-14)
  expect_equal(mf(".mf_transport", a, X, b),
               c(-0.03392772177023339, 0.40074482220507573, 0.40074482220507573, -1.3026637280123452), tolerance = 1e-14)
  expect_equal(mf(".mf_get_coordinates", a, X), c(8.252620095397387e-18, 0.2828427124746189, -0.4000000000000001),
               tolerance = 1e-14)
  expect_equal(mf(".mf_rgrad", a, c(1, 3, 2, 4)), c(10, 8.625, 8.625, 6.75))
  kf <- function(M, p) distance(M, p, A)^2 + distance(M, p, B)^2
  kg <- function(M, p) -2 * (matrix(mf(".mf_log", as.numeric(p), a), 2) + matrix(mf(".mf_log", as.numeric(p), b), 2))
  s15 <- Manopt$StopAfterIteration(15)
  expect_equal(Manopt$gradient_descent(M, kf, kg, diag(2), stopping_criterion = s15),
               matrix(c(1.3754378460243553, 0.13981645384132613, 0.13981645384132624, 1.657483693475966), 2),
               tolerance = 1e-12)
  expect_equal(Manopt$quasi_Newton(M, kf, kg, diag(2), stopping_criterion = s15),
               matrix(c(1.3775242983950455, 0.132693266240773, 0.1326932662407731, 1.6509794599157068), 2),
               tolerance = 1e-12)
  expect_equal(Manopt$trust_regions(M, kf, kg, diag(2), stopping_criterion = Manopt$StopAfterIteration(8)),
               matrix(c(1.3775242983950484, 0.1326932662407666, 0.1326932662407666, 1.6509794599156857), 2),
               tolerance = 1e-12)
  # the Armijo tests of the conjugate gradient and of the default gradient
  # descent come down to cost differences of a few ulps, decided by the
  # rounding of the eigen decompositions: Julia stops after 63 iterations, R
  # after 67
  expect_equal(Manopt$conjugate_gradient_descent(M, kf, kg, diag(2), stopping_criterion = s15),
               matrix(c(1.095091799168431, 0.2706275432363794, 0.2706275432363794, 1.459217149903989), 2),
               tolerance = 1e-7)
  st <- Manopt$gradient_descent(M, kf, kg, diag(2), return_state = TRUE)
  expect_true(st$converged)
  expect_equal(st$p, matrix(c(1.3775242993164076, 0.13269326414674465, 0.13269326414674454, 1.6509794584079749), 2),
               tolerance = 1e-8)
})

test_that("Lorentz is R^n with the Minkowski metric", {
  L <- Lorentz(3)
  ns <- asNamespace("modelcompR")
  expect_equal(get(".mf_inner", envir = ns)(L, NULL, c(1, 2, 3), c(1, 1, 1)), 0)
  expect_error(get(".mf_norm", envir = ns)(L, NULL, c(0, 0, 1)), "DomainError")
  # Manifolds.jl cannot change the representer of a Euclidean gradient to the
  # Minkowski metric, so the Riemannian gradient and Hessian fail as upstream
  msg <- "local_metric(::AbstractManifold, ::Any, ::AbstractBasis) is ambiguous"
  expect_error(get(".mf_rgrad", envir = ns)(L, c(0, 0, 1), c(1, 2, 3)), msg, fixed = TRUE)
  expect_error(get(".mf_rhess", envir = ns)(L, c(0, 0, 1), c(1, 2, 3), c(1, 0, 0), c(0, 1, 0)), msg, fixed = TRUE)
  g <- function(x, p) (x[1] - 1)^2 + (x[2] + 2)^2
  expect_error(solve(OptimizationProblem(g, c(0, 0, 1), manifold = L), GradientDescentOptimizer(), maxiters = 100),
               msg, fixed = TRUE)
  # with the gradient in the Minkowski representation the solver runs
  lg <- function(M, x) c(2 * (x[1] - 1), 2 * (x[2] + 2), -0)
  p <- Manopt$gradient_descent(L, function(M, x) g(x, NULL), lg, c(0, 0, 1),
                               stopping_criterion = Manopt$StopAfterIteration(20))
  expect_equal(p, c(1, -2, 1), tolerance = 1e-2)
})
