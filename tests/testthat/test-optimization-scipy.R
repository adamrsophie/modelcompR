# The native ports of the SciPy optimizers (OptimizationSciPy). Reference
# values come from SciPy 1.18.1 through tools/scipy-probe.py (section named
# in each comment); the iteration and evaluation counts must agree exactly.

sp_rosen <- function(x) { a <- 1 - x[1]; b <- x[2] - x[1] * x[1]; a * a + 100 * b * b }
sp_rosen_grad <- function(x) { b <- x[2] - x[1] * x[1]; c(-2 * (1 - x[1]) - 400 * x[1] * b, 200 * b) }
sp_rosen_hess <- function(x) matrix(c(2 - 400 * x[2] + 1200 * x[1] * x[1], -400 * x[1], -400 * x[1], 200), 2, 2)
sp_quartic <- function(x) x * x * x * x - 3 * x * x * x + 2
sp_poly3 <- function(x) {
  a <- x[1] - 1; b <- x[2] - 2; c <- x[3] + 1
  a * a + 2 * b * b + 3 * c * c + a * a * b * b + 0.5 * a * c
}
sp_poly3_grad <- function(x) {
  a <- x[1] - 1; b <- x[2] - 2; c <- x[3] + 1
  c(2 * a + 2 * a * b * b + 0.5 * c, 4 * b + 2 * a * a * b, 6 * c + 0.5 * a)
}
sp_poly3_hess <- function(x) {
  a <- x[1] - 1; b <- x[2] - 2
  matrix(c(2 + 2 * b * b, 4 * a * b, 0.5, 4 * a * b, 4 + 2 * a * a, 0, 0.5, 0, 6), 3, 3)
}

sp_check <- function(res, ref, tol = 1e-10, label = "") {
  expect_equal(res$x, ref$x, tolerance = tol, info = label)
  expect_lte(abs(res$fun - ref$fun), tol * max(1, abs(ref$fun)), label = paste(label, "fun"))
  for (f in intersect(c("nit", "nfev", "njev", "nhev", "status"), names(ref)))
    expect_identical(as.integer(res[[f]]), ref[[f]], info = paste(label, f))
}

test_that("Nelder-Mead, Powell, CG, BFGS and Newton-CG follow SciPy step by step", {
  x0 <- c(-1.2, 1)
  # python3 tools/scipy-probe.py minimize
  sp_check(.sp_minimize(sp_rosen, x0, "Nelder-Mead"), list(x = c(1.0000220217835696, 1.0000422197517715),
    fun = 8.177661197416674e-10, nit = 85L, nfev = 159L, status = 0L), label = "NM")
  sp_check(.sp_minimize(sp_rosen, x0, "Powell"), list(x = c(1.0000000000000888, 1.0000000000001876),
    fun = 1.787262988391355e-26, nit = 23L, nfev = 607L, status = 0L), label = "Powell")
  sp_check(.sp_minimize(sp_rosen, x0, "CG", jac = sp_rosen_grad), list(x = c(1.0000000054854676, 0.9999999972757344),
    fun = 1.878594302092354e-14, nit = 36L, nfev = 78L, njev = 77L, status = 0L), label = "CG")
  sp_check(.sp_minimize(sp_rosen, x0, "BFGS", jac = sp_rosen_grad), list(x = c(0.9999999710013947, 0.999999946119085),
    fun = 2.5353073555215442e-15, nit = 32L, nfev = 39L, njev = 39L, status = 0L), label = "BFGS")
  sp_check(.sp_minimize(sp_poly3, c(0, 0, 0), "Nelder-Mead"), list(x = c(1.000044414318398, 2.00002316080351,
    -1.0000074537872385), fun = 3.0466267110382427e-09, nit = 158L, nfev = 288L, status = 0L), label = "NM3")
  sp_check(.sp_minimize(sp_poly3, c(0, 0, 0), "Powell"), list(x = c(1.0, 2.0, -0.9999999999999999),
    fun = 3.697785493223493e-32, nit = 3L, nfev = 107L, status = 0L), label = "Powell3")
  sp_check(.sp_minimize(sp_poly3, c(0, 0, 0), "BFGS"), list(x = c(0.9999998324943379, 2.0000000184440454,
    -0.9999998303576163), fun = 1.0086609759947917e-13, nit = 8L, nfev = 44L, njev = 11L, status = 0L), label = "BFGSfd")
  sp_check(.sp_minimize(sp_poly3, c(0, 0, 0), "CG"), list(x = c(0.999997719749209, 2.000001502156387,
    -0.9999990015337875), fun = 1.1564918937770767e-11, nit = 7L, nfev = 56L, njev = 14L, status = 0L), label = "CGfd")
  sp_check(.sp_minimize(sp_rosen, x0, "Newton-CG", jac = sp_rosen_grad, hess = sp_rosen_hess),
    list(x = c(0.9999826029170281, 0.9999651365366965), fun = 3.0314291218378575e-10, nit = 83L, nfev = 105L,
         njev = 105L, nhev = 83L, status = 0L), label = "NewtonCG")
  b <- list(lb = c(-2, 0.5), ub = c(0.8, 2))
  sp_check(.sp_minimize(sp_rosen, x0, "Nelder-Mead", bounds = b), list(x = c(-0.6984289444237956, 0.5),
    fun = 2.8995375835529753, nit = 36L, nfev = 66L, status = 0L), label = "NMb")
  sp_check(.sp_minimize(sp_rosen, x0, "Powell", bounds = b), list(x = c(-0.6857643396781444, 0.5000234086049975),
    fun = 2.9303116992098355, nit = 4L, nfev = 152L, status = 0L), label = "Powellb")
  sp_check(.sp_minimize(sp_rosen, x0, "BFGS", jac = sp_rosen_grad, options = list(maxiter = 5)),
    list(x = c(-0.6382321094586818, 0.3779285684797238), fun = 2.770309001588565, nit = 5L, nfev = 8L, njev = 8L,
         status = 1L), label = "BFGSmax")
})

test_that("the scalar minimizers and brute follow SciPy", {
  # python3 tools/scipy-probe.py scalar
  sp_check(.sp_minimize_scalar(sp_quartic, method = "brent"), list(x = 2.2500000105090607, fun = -6.54296875,
    nit = 10L, nfev = 14L), label = "brent")
  sp_check(.sp_minimize_scalar(sp_quartic, method = "golden"), list(x = 2.250000012987402,
    fun = -6.5429687500000036, nit = 38L, nfev = 44L), label = "golden")
  sp_check(.sp_minimize_scalar(sp_quartic, bounds = c(-1, 4), method = "bounded"), list(x = 2.250000082565353,
    fun = -6.5429687499999325, nit = 13L, nfev = 13L, status = 0L), label = "bounded")
  sp_check(.sp_brute(sp_rosen, c(-2, -1), c(2, 3)), list(x = c(1.0000194819756114, 1.000038536536097),
    fun = 3.978482019500388e-10), label = "brute")
})

test_that("the trust region methods follow SciPy", {
  x0 <- c(-1.2, 1)
  tro <- list(gtol = 1e-10, maxiter = 50000)
  # python3 tools/scipy-probe.py trust
  sp_check(.sp_minimize(sp_rosen, x0, "dogleg", jac = sp_rosen_grad, hess = sp_rosen_hess, options = tro),
    list(x = c(1, 1), fun = 0, nit = 25L, nfev = 26L, njev = 23L, nhev = 22L, status = 0L), label = "dogleg")
  sp_check(.sp_minimize(sp_rosen, x0, "trust-ncg", jac = sp_rosen_grad, hess = sp_rosen_hess, options = tro),
    list(x = c(0.9999999999999261, 0.9999999999997595), fun = 8.62806211982294e-25, nit = 30L, nfev = 31L,
         njev = 28L, nhev = 27L, status = 0L), label = "ncg")
  sp_check(.sp_minimize(sp_rosen, x0, "trust-exact", jac = sp_rosen_grad, hess = sp_rosen_hess, options = tro),
    list(x = c(1, 1), fun = 0, nit = 26L, nfev = 27L, njev = 24L, nhev = 27L, status = 0L), label = "exact")
  sp_check(.sp_minimize(sp_rosen, x0, "trust-ncg", jac = sp_rosen_grad, hess = .sp_hessian_update("BFGS"),
                        options = tro),
    list(x = c(0.9999999999999986, 0.9999999999999981), fun = 1.0192329414488354e-28, nit = 50L, nfev = 51L,
         njev = 53L, nhev = 122L, status = 0L), label = "ncg BFGS")
  sp_check(.sp_minimize(sp_poly3, c(3, -1, 0.5), "trust-exact", jac = sp_poly3_grad, hess = sp_poly3_hess,
                        options = tro),
    list(x = c(1, 2, -1), fun = 0, nit = 7L, nfev = 8L, njev = 8L, nhev = 8L, status = 0L), label = "exact3")
  # an indefinite Hessian at x0: the Cholesky factorization fails at once
  sp_check(.sp_minimize(sp_poly3, c(3, -1, 0.5), "dogleg", jac = sp_poly3_grad, hess = sp_poly3_hess, options = tro),
    list(x = c(3, -1, 0.5), fun = 66.25, nit = 0L, nfev = 1L, njev = 1L, nhev = 1L, status = 3L), label = "dogleg3")
})

test_that("L-BFGS-B follows SciPy", {
  x0 <- c(-1.2, 1)
  L <- "L-BFGS-B"
  # python3 tools/scipy-probe.py lbfgsb
  sp_check(.sp_minimize(sp_rosen, x0, L, jac = sp_rosen_grad), list(x = c(0.9999989453492616, 0.9999980209056968),
    fun = 2.807650019071083e-12, nit = 36L, nfev = 44L, njev = 44L, status = 0L), label = "plain")
  sp_check(.sp_minimize(sp_rosen, x0, L, jac = sp_rosen_grad, options = list(maxcor = 2)),
    list(x = c(1.0000001654861872, 1.0000003405967195), fun = 3.664842744018677e-14, nit = 38L, nfev = 51L,
         njev = 51L, status = 0L), label = "m2")
  sp_check(.sp_minimize(sp_rosen, x0, L, jac = sp_rosen_grad, bounds = list(lb = c(-2, 1.5), ub = c(0.8, 3))),
    list(x = c(-1.2210262420670972, 1.5), fun = 4.941229317989185, nit = 4L, nfev = 6L, njev = 6L, status = 0L),
    label = "box")
  sp_check(.sp_minimize(sp_rosen, c(2, 3), L, jac = sp_rosen_grad, bounds = list(lb = c(1.2, -Inf), ub = c(Inf, Inf))),
    list(x = c(1.2, 1.44), fun = 0.03999999999999998, nit = 17L, nfev = 21L, njev = 21L, status = 0L), label = "lower")
  sp_check(.sp_minimize(sp_poly3, c(0, 0, 0), L, bounds = list(lb = c(-1, -1, -0.5), ub = c(0.5, 1.5, 1))),
    list(x = c(0.5, 1.5, -0.5), fun = 1.4375, nit = 1L, nfev = 8L, njev = 2L, status = 0L), label = "fd box")
  sp_check(.sp_minimize(sp_poly3, c(0, 0, 0), L, jac = sp_poly3_grad,
                        bounds = list(lb = rep(-Inf, 3), ub = c(Inf, 1, Inf)), options = list(maxcor = 1)),
    list(x = c(1.0000031960435822, 1.0, -1.0000013198340452), fun = 2.000000000023546, nit = 8L, nfev = 9L,
         njev = 9L, status = 0L), label = "upper m1")
  sp_check(.sp_minimize(sp_rosen, x0, L, jac = sp_rosen_grad, options = list(maxiter = 7)),
    list(x = c(-0.9493671678671748, 0.8669746872519423), fun = 3.917841468298796, nit = 7L, nfev = 9L, njev = 9L,
         status = 1L), label = "maxiter")
  # SciPy's compiled lbfgsb.c fuses stp * d + t into one multiply add, so the
  # finite difference run drifts from the second iterate on in the last bit;
  # the counts still agree and the minimizer to 1e-7
  sp_check(.sp_minimize(sp_rosen, x0, L), list(x = c(0.9999959932564302, 0.9999921141462681),
    fun = 1.7682612930421444e-11, nit = 36L, nfev = 132L, njev = 44L, status = 0L), tol = 1e-6, label = "fd")
})

test_that("the Scipy optimizers solve problems through solve()", {
  rosen <- function(u, p) { a <- p[1] - u[1]; b <- u[2] - u[1] * u[1]; a * a + p[2] * b * b }
  rg <- function(u, p) { b <- u[2] - u[1] * u[1]; c(-2 * (p[1] - u[1]) - 4 * p[2] * u[1] * b, 2 * p[2] * b) }
  rh <- function(u, p) matrix(c(2 - 4 * p[2] * u[2] + 12 * p[2] * u[1] * u[1], -4 * p[2] * u[1], -4 * p[2] * u[1],
                                2 * p[2]), 2, 2)
  prob <- OptimizationProblem(rosen, c(-1.2, 1), c(1, 100), grad = rg, hess = rh)
  for (alg in list(ScipyNelderMead(), ScipyPowell(), ScipyCG(), ScipyBFGS(), ScipyNewtonCG(), ScipyLBFGSB(),
                   ScipyDogleg(), ScipyTrustNCG(), ScipyTrustExact())) {
    sol <- solve(prob, alg)
    expect_equal(sol$retcode, "Success", info = alg$method)
    expect_equal(unname(sol$u), c(1, 1), tolerance = 1e-4, info = alg$method)
  }
  expect_equal(solve(prob, ScipyBFGS())$stats$iterations, 32L)
  bx <- OptimizationProblem(rosen, c(-1.2, 1), c(1, 100), lb = c(-2, -1), ub = c(2, 3))
  expect_equal(unname(solve(bx, ScipyBrute())$u), c(1, 1), tolerance = 1e-4)
  q <- OptimizationProblem(function(u, p) u * u * u * u - 3 * u * u * u + 2, 0.5, numeric(), lb = -1, ub = 4)
  for (alg in list(ScipyBrent(), ScipyGolden(), ScipyBounded()))
    expect_equal(unname(solve(q, alg)$u), 2.25, tolerance = 1e-6, info = alg$method)
  expect_error(solve(OptimizationProblem(rosen, c(-1.2, 1), c(1, 100)), ScipyBrute()), "requires bounds")
  expect_error(solve(prob, ScipyBrent()), "exactly 1 variable")
  expect_error(ScipyMinimize("nope"), "Invalid method")
})

test_that("the scalar root finders follow SciPy", {
  cub <- function(x) x * x * x - 2 * x - 5
  cub_p <- function(x) 3 * x * x - 2
  cub_pp <- function(x) 6 * x
  rc <- function(r, root, iterations, calls, converged = TRUE, tol = 0) {
    expect_equal(r$root, root, tolerance = tol)
    expect_equal(r$iterations, iterations)
    expect_equal(r$function_calls, calls)
    expect_identical(r$converged, converged)
  }
  # python3 tools/scipy-probe.py zeros
  rc(.sp_root_scalar(cub, "brentq", bracket = c(2, 3)), 2.094551481542327, 7, 8)
  rc(.sp_root_scalar(function(x) cos(x) - x, "brentq", bracket = c(0, 1)), 0.7390851332151559, 7, 8)
  rc(.sp_root_scalar(cub, "brentq", bracket = c(-1, 4)), 2.0945514815423163, 11, 12)
  rc(.sp_root_scalar(cub, "brenth", bracket = c(2, 3)), 2.0945514815423265, 7, 8)
  rc(.sp_root_scalar(cub, "brenth", bracket = c(-1, 4)), 2.0945514815423265, 11, 12)
  rc(.sp_root_scalar(cub, "bisect", bracket = c(2, 3)), 2.0945514815412025, 39, 41)
  rc(.sp_root_scalar(cub, "bisect", bracket = c(-1, 4)), 2.09455148154143, 42, 44)
  rc(.sp_root_scalar(cub, "ridder", bracket = c(2, 3)), 2.094551481543327, 6, 14)
  rc(.sp_root_scalar(function(x) cos(x) - x, "ridder", bracket = c(0, 1)), 0.7390851332151607, 5, 12)
  # one unit in the last place: the compiled ridder.c contracts to a fused multiply add
  rc(.sp_root_scalar(cub, "ridder", bracket = c(-1, 4)), 2.0945514815414166, 6, 14, tol = 1e-14)
  rc(.sp_root_scalar(cub, "toms748", bracket = c(2, 3)), 2.094551481542327, 5, 11)
  rc(.sp_root_scalar(function(x) cos(x) - x, "toms748", bracket = c(0, 1)), 0.7390851332151607, 4, 8)
  rc(.sp_root_scalar(cub, "toms748", bracket = c(-1, 4)), 2.0945514815425765, 7, 16)
  rc(.sp_root_scalar(cub, "newton", x0 = 1.5, fprime = cub_p), 2.0945514815423265, 6, 12)
  rc(.sp_root_scalar(cub, "newton", x0 = 1.5), 2.0945514815423265, 6, 12)
  rc(.sp_root_scalar(cub, "secant", x0 = 1.5), 2.094551481542326, 8, 9)
  rc(.sp_root_scalar(cub, "secant", x0 = 1.5, x1 = 2.5), 2.094551481542326, 7, 8)
  rc(.sp_root_scalar(cub, "halley", x0 = 1.5, fprime = cub_p, fprime2 = cub_pp), 2.0945514815423265, 4, 12)
  rc(.sp_root_scalar(cub, "brentq", bracket = c(2, 3), xtol = 1e-3, maxiter = 3), 2.094528891117347, 3, 5, FALSE)
  expect_error(.sp_root_scalar(cub, "brentq", bracket = c(3, 4)), "different signs")
  # through solve(): the bracket comes from the bounds, x0 from u0
  pr <- OptimizationProblem(function(u, p) u * u * u - 2 * u - 5, 1.5, numeric(), lb = 2, ub = 3,
                            grad = function(u, p) 3 * u * u - 2, hess = function(u, p) matrix(6 * u, 1, 1))
  for (m in c("brentq", "brenth", "bisect", "ridder", "toms748", "newton", "secant", "halley")) {
    sol <- solve(pr, ScipyRootScalar(m))
    expect_equal(sol$retcode, "Success", info = m)
    expect_equal(unname(sol$u), 2.0945514815423265, tolerance = 1e-9, info = m)
  }
})

test_that("differential evolution, basin hopping and dual annealing converge", {
  # SciPy (tools/scipy-probe.py is not used: the random streams differ) finds
  # the same optima: Rosenbrock (1, 1); x^4 - 4 x^2 + x per coordinate has its
  # global minimum -10.888384133221 at x = -1.4729976 and a local one at 1.3469974
  rosen <- function(u, p) { a <- 1 - u[1]; b <- u[2] - u[1] * u[1]; a * a + 100 * b * b }
  w <- function(u, p) sum(u * u * u * u - 4 * u * u + u)
  pb <- OptimizationProblem(rosen, c(-1.2, 1), numeric(), lb = c(-2, -1), ub = c(2, 3))
  pw <- OptimizationProblem(w, c(1.5, 1.5), numeric(), lb = c(-3, -3), ub = c(3, 3))
  de <- solve(pw, ScipyDifferentialEvolution(), seed = 1)
  expect_equal(de$retcode, "Success")
  expect_equal(unname(de$u), rep(-1.4729976, 2), tolerance = 1e-6)
  expect_equal(de$objective, -10.888384133221, tolerance = 1e-10)
  expect_equal(unname(solve(pb, ScipyDifferentialEvolution(), seed = 2)$u), c(1, 1), tolerance = 1e-6)
  da <- solve(pw, ScipyDualAnnealing(), seed = 3)
  expect_equal(unname(da$u), rep(-1.4729976, 2), tolerance = 1e-6)
  expect_equal(da$original$nit, 1000L)
  # basin hopping with the default step 0.5 stays in the basin of x0, as in SciPy
  bh <- solve(pw, ScipyBasinhopping(), seed = 4, maxiters = 20)
  expect_equal(unname(bh$u), rep(1.3469974, 2), tolerance = 1e-6)
  expect_equal(unname(solve(pb, ScipyBasinhopping(), seed = 5, maxiters = 10)$u), c(1, 1), tolerance = 1e-4)
  expect_error(solve(OptimizationProblem(rosen, c(0, 0), numeric()), ScipyDualAnnealing()), "requires bounds")
})

test_that("COBYLA follows SciPy and PRIMA", {
  nlc <- .sp_nonlinear_constraint
  # python3 tools/scipy-probe.py cobyla: with active bounds the runs agree exactly
  r <- .sp_minimize(sp_poly3, c(0, 0, 0), method = "COBYLA", bounds = list(lb = c(-1, -1, -0.5), ub = c(0.5, 1.5, 1)))
  expect_equal(r$x, c(0.5, 1.5, -0.5))
  expect_equal(c(r$nfev, r$status), c(14L, 0L))
  r <- .sp_minimize(sp_poly3, c(0, 0, 0), method = "COBYLA", bounds = list(lb = c(-1, 0.25, -0.5), ub = c(0.5, 0.25, 1)))
  expect_equal(r$x, c(0.5, 0.25, -0.5))
  expect_equal(c(r$nfev, r$status), c(8L, 0L))
  # elsewhere the evaluated points agree with SciPy to rounding for the first
  # 30 to 60 evaluations; then last bit differences of numpy's BLAS products
  # change the path (SciPy started 1e-15 away diverges at the same place) and
  # the results agree within rhoend = 1e-4
  r <- .sp_minimize(sp_poly3, c(0, 0, 0), method = "COBYLA")
  expect_equal(r$x, c(1.0000222657589217, 1.9999751185226917, -0.9999636052688081), tolerance = 2e-4)
  expect_lt(r$fun, 1e-7)
  expect_equal(r$status, 0L)
  r <- .sp_minimize(sp_poly3, c(0, 0, 0), method = "COBYLA",
                    constraints = list(nlc(function(x) c(x[1] + x[2] * x[2], x[1] * x[3]), c(1, -0.5), c(1, 0.5))))
  expect_equal(r$x, c(0.2592873533505713, 0.8606466521080087, -0.9383308308111745), tolerance = 1e-3)
  expect_equal(r$fun, 3.8457006607956044, tolerance = 1e-7)
  expect_lt(r$maxcv, sqrt(.Machine$double.eps))
  expect_true(r$success)
  r <- .sp_minimize(sp_rosen, c(-1.2, 1), method = "COBYLA",
                    constraints = list(nlc(function(x) x[1] * x[1] + x[2] * x[2], -Inf, 1.5)))
  expect_equal(c(r$nfev, r$status), c(1000L, 3L))
  expect_lt(r$fun, 0.03)
  # PRIMA.jl (libprima 0.7.1, julia tools/prima-probe.jl) reaches the same
  # points. Its evaluated points agree with the port to rounding for the
  # first 8 to 28 evaluations and then split, earlier than libprima started
  # 1e-15 away splits from itself (28 and 111 evaluations): libprima 0.7.1
  # differs in some step from the newer PRIMA that pyprima translates, whose
  # sources this port follows, so the counts differ
  calc <- function(x) list(f = sp_poly3(x), c = sum(x * x) - 4)
  r <- .pr_cobyla(calc, 1L, c(0, 0, 0), f0 = sp_poly3(c(0, 0, 0)), nlconstr0 = -4, rhobeg = 1, rhoend = 1e-6,
                  maxfun = 1500L, fortran = TRUE)
  expect_equal(r$x, c(0.702592382241251, 1.6656076887922722, -0.8556371738634042), tolerance = 1e-5)
  expect_equal(r$info, 0L)
  calc <- function(x) list(f = sp_poly3(x), c = x[1] * x[3] - 0.5)
  r <- .pr_cobyla(calc, 1L, c(0, 0, 0), Aeq = matrix(1, 1, 3), beq = 1, xl = c(-1, -1, -2), xu = c(2, 1.5, 2),
                  f0 = sp_poly3(c(0, 0, 0)), nlconstr0 = -0.5, rhobeg = 1, rhoend = 1e-6, maxfun = 1500L,
                  fortran = TRUE)
  expect_equal(r$x, c(0.6333331479738197, 1.5, -1.1333331479738198), tolerance = 1e-7)
  expect_equal(r$f, 0.7458333333334622, tolerance = 1e-10)
  # through solve(): ScipyCOBYLA and COBYLA with the constraint bounds
  pr <- OptimizationProblem(function(u, p) sp_poly3(u), c(0, 0, 0), numeric(),
                            cons = function(u, p) c(u[1] + u[2] + u[3], sum(u * u)), lcons = c(1, -Inf),
                            ucons = c(1, 4))
  s1 <- solve(pr, ScipyCOBYLA(), maxiters = 3000)
  s2 <- solve(pr, COBYLA())
  expect_equal(s1$retcode, "Success")
  expect_equal(s2$retcode, "Success")
  expect_equal(unname(s1$u), unname(s2$u), tolerance = 1e-3)
  expect_equal(sum(s2$u), 1, tolerance = 1e-8)
  expect_lte(sum(s2$u * s2$u), 4 + 1e-8)
  expect_error(solve(pr, NEWUOA()), "not yet ported")
})

test_that("COBYLA evaluates the same points as SciPy and libprima at the start", {
  nlc <- .sp_nonlinear_constraint
  pts <- list()
  rec <- function(g) function(x) { pts[[length(pts) + 1L]] <<- x; g(x) }
  # the 25th evaluated point of SciPy (trace of the objective in SciPy 1.18.1)
  .sp_minimize(rec(sp_poly3), c(0, 0, 0), method = "COBYLA",
               constraints = list(nlc(function(x) c(x[1] + x[2] * x[2], x[1] * x[3]), c(1, -0.5), c(1, 0.5))))
  expect_equal(pts[[25]], c(0.2646842756093336, 0.8575881244263656, -0.919174120610123), tolerance = 1e-12)
  pts <- list()
  .sp_minimize(rec(sp_rosen), c(-1.2, 1), method = "COBYLA")
  expect_equal(pts[[50]], c(-0.9600662134855054, 0.9842520051894493), tolerance = 1e-12)
  # the 25th evaluated point of PRIMA.jl (libprima 0.7.1)
  pts <- list()
  f <- rec(sp_rosen)
  calc <- function(x) list(f = f(x), c = -x[1])
  .pr_cobyla(calc, 1L, c(-1.2, 1), f0 = f(c(-1.2, 1)), nlconstr0 = 1.2, rhobeg = 0.5, rhoend = 1e-8, maxfun = 400L,
             fortran = TRUE)
  expect_equal(pts[[25]], c(0.3547375871480715, 0.1258050185441), tolerance = 1e-12)
})

test_that("ScipyRoot and ScipyLeastSquaresLM follow MINPACK of SciPy", {
  # python3 tools/scipy-probe.py root (SciPy 1.18.1)
  sys2 <- function(x, p) { d <- x[1] - x[2]; c(x[1] + 0.5 * d * d * d - 1, -0.5 * d * d * d + x[2]) }
  sys2j <- function(x, p) {
    d <- x[1] - x[2]
    matrix(c(1 + 1.5 * d * d, -1.5 * d * d, -1.5 * d * d, 1 + 1.5 * d * d), 2, 2, byrow = TRUE)
  }
  sys3 <- function(x, p) c(x[1] * x[1] + x[2] * x[2] - 4, x[1] * x[2] - 1, x[3] * x[3] * x[3] + x[1] - 2)
  chk <- function(s, x, nfev, status, njev = NULL, tol = 1e-15) {
    expect_equal(unname(s$u), x, tolerance = tol)
    expect_identical(s$original$nfev, nfev)
    expect_identical(s$original$status, status)
    if (!is.null(njev)) expect_identical(s$original$njev, njev)
  }
  s <- solve(OptimizationProblem(sys2, c(0, 0), numeric()), ScipyRoot())
  chk(s, c(0.8411639019140096, 0.1588360980859903), 14L, 1L)
  expect_equal(s$retcode, "Success")
  expect_identical(s$original$fun, c(-2^-53, 0))
  pj <- OptimizationProblem(OptimizationFunction(sys2, grad = sys2j), c(0, 0), numeric())
  chk(solve(pj, ScipyRoot("hybr")), c(0.8411639019140096, 0.1588360980859903), 12L, 1L, njev = 1L)
  chk(solve(pj, ScipyRoot("lm")), c(0.8411639019140097, 0.15883609808599033), 9L, 2L, njev = 6L)
  p3 <- OptimizationProblem(sys3, c(1, 1, 1), numeric())
  chk(solve(p3, ScipyRoot("hybr")), c(0.517638090204806, 1.9318516525792895, 1.1402097464162466), 25L, 1L)
  # the Jacobian at (1, 1, 1) is singular: SciPy's FMA contracted MINPACK
  # takes another path (nfev 34) than the C code without contraction, which
  # the port follows to 1 ulp in all 25 points (tools/minpack-nofma-check.c
  # gives the point and nfev 25 + 2)
  s <- solve(p3, ScipyRoot("lm"))
  expect_equal(unname(s$u), c(0.5176380902050415, 1.9318516525781366, 1.1402097464162713), tolerance = 1e-12)
  expect_identical(s$original$nfev, 27L)
  p3f <- OptimizationProblem(sys3, c(3, -2, 5), numeric())
  chk(solve(p3f, ScipyRoot("hybr")), c(1.9318516525781135, 0.5176380902054615, 0.4084621098126088), 32L, 1L)
  chk(solve(p3f, ScipyRoot("lm")), c(1.9318516525781366, 0.5176380902050415, 0.40846210981256387), 47L, 2L)
  s <- solve(OptimizationProblem(sys2, c(0, 0), numeric()), ScipyRoot("lm"))
  expect_equal(unname(s$u), c(0.7252038049458914, 1.1102230246251565e-16), tolerance = 1e-9)
  expect_identical(s$original$nfev, 21L)
  expect_error(solve(p3, ScipyRoot(), xtol = 1e-3), "unexpected keyword argument 'xtol'")
  expect_error(solve(p3, ScipyRoot("krylov")), "not yet ported")
  res3 <- function(x, p) c(10 * (x[2] - x[1] * x[1]), 1 - x[1], 0.5 * (x[1] * x[2] - 1))
  s <- solve(OptimizationProblem(res3, c(-1.2, 1), numeric()), ScipyLeastSquaresLM())
  chk(s, c(1, 1), 23L, 3L)
  expect_equal(s$objective, 0)
  s <- solve(p3f, ScipyLeastSquaresLM())
  chk(s, c(1.9318516525781366, 0.5176380902050416, 0.40846210981256365), 12L, 3L)
  expect_equal(s$original$optimality, 4.2895723756260546e-16, tolerance = 1e-6)
})

test_that("ScipyLeastSquares trf and dogbox follow SciPy", {
  # python3 tools/scipy-probe.py lsq (SciPy 1.18.1)
  res3 <- function(x, p) c(10 * (x[2] - x[1] * x[1]), 1 - x[1], 0.5 * (x[1] * x[2] - 1))
  sys3 <- function(x, p) c(x[1] * x[1] + x[2] * x[2] - 4, x[1] * x[2] - 1, x[3] * x[3] * x[3] + x[1] - 2)
  tt <- seq(0, 2, length.out = 9)
  yy <- 1 + 0.5 * tt - 0.3 * tt * tt + c(0.1, -0.2, 0.05, 3.0, -0.1, 0.0, 0.2, -2.5, 0.1)
  fit <- function(x, p) x[1] + x[2] * tt + x[3] * tt * tt * x[1] - yy
  chk <- function(s, x, nfev, njev, status, tol = 1e-10) {
    expect_equal(unname(s$u), x, tolerance = tol)
    expect_identical(s$original$nfev, nfev)
    expect_identical(s$original$njev, njev)
    expect_identical(s$original$status, status)
  }
  pr <- OptimizationProblem(res3, c(-1.2, 1), numeric())
  prb <- OptimizationProblem(res3, c(-1.2, 1), numeric(), lb = c(-2, -2), ub = c(0.8, 2))
  ps <- OptimizationProblem(sys3, c(3, -2, 5), numeric())
  s <- solve(pr, ScipyLeastSquares())
  chk(s, c(1.0000000000000209, 0.9999999999997212), 28L, 20L, 1L, tol = 1e-15)
  expect_equal(s$retcode, "Success")
  chk(solve(prb, ScipyLeastSquaresTRF()), c(0.7999999999999999, 0.6409744407720932), 26L, 21L, 2L)
  chk(solve(ps, ScipyLeastSquaresTRF()), c(1.9318516525781366, 0.5176380902050416, 0.40846210981613446), 11L, 11L, 1L)
  chk(solve(pr, ScipyLeastSquaresDogbox()), c(0.999999999999921, 0.9999999999872115), 20L, 16L, 1L)
  s <- solve(prb, ScipyLeastSquaresDogbox())
  chk(s, c(0.8, 0.640974440890955), 18L, 14L, 1L, tol = 1e-15)
  expect_identical(s$original$active_mask, c(1L, 0L))
  # robust losses: the forward difference Jacobians amplify last bit
  # differences, and SciPy itself started 1e-15 away moves by 1e-7 and may
  # change its counts (trf huber: 5, 5, 2 from 1 and 7, 4, 3 from 1 - 1e-15,
  # the counts of the port), so only the points are compared
  pf <- OptimizationProblem(fit, c(1, 0, 0), numeric())
  s <- solve(pf, ScipyLeastSquares(loss = "cauchy"))
  expect_equal(unname(s$u), c(1.0122144703015732, 0.616196867559221, -0.37557678302508185), tolerance = 1e-7)
  expect_equal(s$objective, 2.191113801044832, tolerance = 1e-12)
  s <- solve(pf, ScipyLeastSquares(method = "dogbox", loss = "arctan"))
  expect_equal(unname(s$u), c(0.9988789180083442, 0.4218010377324867, -0.2295639729307764), tolerance = 1e-7)
  s <- solve(OptimizationProblem(fit, c(1, 0.7, 0), numeric(), lb = c(0, 0.6, -1), ub = c(2, 2, 1)),
             ScipyLeastSquares(loss = "soft_l1"))
  expect_equal(unname(s$u), c(1.0217155317180437, 1.1621250193083985, -0.7612480455712755), tolerance = 1e-6)
  expect_error(solve(pr, ScipyLeastSquares(), tr_solver = "lsmr"), "not yet ported")
  # df-sane
  for (ls in c("cruz", "cheng")) {
    s <- solve(ps, ScipyRoot("df-sane"), options = list(line_search = ls))
    expect_equal(unname(s$u), c(1.9318516509800148, 0.5176381034564093, 0.4084640918648213), tolerance = 1e-15)
    expect_identical(s$original$nfev, 31L)
    expect_identical(s$original$nit, 26L)
  }
})

test_that("ScipyLinprog and ScipyMilp reach the optimum of HiGHS", {
  # python3 tools/scipy-probe.py lp (SciPy 1.18.1, HiGHS); iteration and
  # node counts differ from HiGHS
  cost <- function(x, p) c(-1, -2, 0.5)
  A1 <- rbind(c(1, 1, 1), c(2, 0.5, 0), c(0, 1, -1))
  pr <- OptimizationProblem(cost, c(0, 0, 0), numeric(), lb = c(0, -1, -Inf), ub = c(3, Inf, 2))
  for (m in c("highs", "highs-ds", "highs-ipm", "simplex")) {
    s <- solve(pr, ScipyLinprog(m), A_ub = A1, b_ub = c(4, 5, 1.5), A_eq = rbind(c(1, -1, 2)), b_eq = 0.5)
    expect_equal(unname(s$u), c(1.5, 2, 0.5), tolerance = 1e-14)
    expect_equal(s$objective, -5.25, tolerance = 1e-14)
    expect_equal(s$retcode, "Success")
  }
  o <- s$original
  expect_equal(o$slack, c(0, 1, 0), tolerance = 1e-14)
  expect_equal(o$ineqlin$marginals, c(-0.5, 0, -2), tolerance = 1e-14)
  expect_equal(o$eqlin$marginals, -0.5, tolerance = 1e-14)
  expect_equal(o$lower$marginals, c(0, 0, 0))
  expect_equal(o$upper$marginals, c(0, 0, 0))
  expect_identical(o$message, "Optimization terminated successfully.")
  # default bounds x >= 0
  pf <- OptimizationProblem(function(x, p) c(1, 1), c(0, 0), numeric())
  s <- solve(pf, ScipyLinprog(), A_ub = rbind(c(-1, -2), c(-3, -1)), b_ub = c(-4, -6))
  expect_equal(unname(s$u), c(1.6, 1.2), tolerance = 1e-14)
  expect_equal(s$original$ineqlin$marginals, c(-0.4, -0.2), tolerance = 1e-14)
  expect_identical(s$original$message, "Optimization terminated successfully. (HiGHS Status 7: Optimal)")
  # the maximum with MaxSense
  s <- solve(OptimizationProblem(function(x, p) c(1, 1), c(0, 0), numeric(), sense = MaxSense), ScipyLinprog(),
             A_ub = rbind(c(1, 2), c(3, 1)), b_ub = c(4, 6))
  expect_equal(unname(s$u), c(1.6, 1.2), tolerance = 1e-14)
  expect_equal(s$objective, 2.8, tolerance = 1e-14)
  s <- solve(pf, ScipyLinprog(), A_ub = rbind(c(1, 1)), b_ub = -1)
  expect_identical(s$original$status, 2L)
  expect_true(all(is.nan(s$u)))
  s <- solve(OptimizationProblem(function(x, p) c(-1, 0), c(0, 0), numeric()), ScipyLinprog(),
             A_ub = rbind(c(1, -1)), b_ub = 1)
  expect_identical(s$original$status, 3L)
  expect_error(solve(pf, ScipyLinprog(), A_ub = rbind(c(1, 1))), "must be provided together")
  # milp
  pm <- OptimizationProblem(function(x, p) c(-5, -4, -3, -1), numeric(4), numeric(), lb = numeric(4),
                            ub = c(10, 10, 10, 2.5))
  s <- solve(pm, ScipyMilp(), integrality = c(1, 1, 1, 0), A = rbind(c(2, 3, 1, 1), c(4, 1, 2, 0.5), c(3, 4, 2, 1)),
             lb_con = rep(-Inf, 3), ub_con = c(5, 11, 8))
  expect_equal(unname(s$u), c(2, 0, 1, 0), tolerance = 1e-12)
  expect_equal(s$objective, -13)
  expect_equal(s$original$mip_gap, 0)
  s <- solve(OptimizationProblem(function(x, p) c(-1, -1), c(0, 0), numeric(), lb = c(0, 0), ub = c(Inf, Inf)),
             ScipyMilp(), integrality = c(1, 1), A = rbind(c(2, 2), c(1, -1)), lb_con = c(-Inf, 0), ub_con = c(7, 0))
  expect_equal(unname(s$u), c(1, 1), tolerance = 1e-12)
  s <- solve(OptimizationProblem(function(x, p) 1, 0.5, numeric(), lb = 0.2, ub = 0.8), ScipyMilp(), integrality = 1)
  expect_identical(s$original$status, 2L)
  expect_equal(s$retcode, "Failure")
  expect_error(solve(pm, ScipyMilp(), A = diag(4)), "must all be provided")
})

test_that("ScipyRoot follows the quasi Newton methods of SciPy", {
  # python3 tools/scipy-probe.py nonlin (SciPy 1.18.1); the iterates agree to
  # rounding (BLAS dot products) and the iteration counts exactly
  fq <- function(x, p) c(x[1] + 0.1 * x[1] * x[1] * x[1] - 0.2 * x[2] - 1,
                         x[2] + 0.1 * x[2] * x[2] * x[2] + 0.2 * x[1] - 0.5, x[3] + 0.05 * x[1] * x[3] - 0.3)
  pr <- OptimizationProblem(fq, c(0, 0, 0), numeric())
  chk <- function(s, x, nit, status, tol = 1e-13) {
    expect_equal(unname(s$u), x, tolerance = tol)
    expect_identical(s$original$nit, nit)
    expect_identical(s$original$status, status)
  }
  chk(solve(pr, ScipyRoot("broyden1")), c(0.9695267682599797, 0.30330450866644454, 0.28613046567461), 12L, 1L)
  chk(solve(pr, ScipyRoot("anderson")), c(0.9695278760331499, 0.3033051791379116, 0.2861295472241734), 9L, 1L)
  chk(solve(pr, ScipyRoot("diagbroyden")), c(0.9695261675941134, 0.30330430003049885, 0.28612952972334454), 13L, 1L)
  # the mixing methods diverge on this system, along SciPy's path
  s <- solve(pr, ScipyRoot("excitingmixing"))
  chk(s, c(-707.0615197960445, -91.84936439513639, -2.9715431643343093), 400L, 2L)
  expect_equal(s$retcode, "Failure")
  chk(solve(pr, ScipyRoot("linearmixing"), options = list(nit = 4)),
      c(-3.378163311653144, -2.0445152309148957, -0.9251138496924428), 5L, 2L)
  chk(solve(pr, ScipyRoot("broyden1"), options = list(line_search = "wolfe")),
      c(0.9695269927537431, 0.30330982586672206, 0.28612627583123224), 12L, 1L, tol = 1e-11)
  chk(solve(pr, ScipyRoot("broyden1"), options = list(jac_options = list(reduction_method = "svd", max_rank = 2))),
      c(0.9695261068487376, 0.30330930164218173, 0.28612863832394325), 85L, 1L)
  chk(solve(pr, ScipyRoot("broyden2"), options = list(jac_options = list(reduction_method = "simple", max_rank = 2))),
      c(0.969527297407631, 0.30330196780840407, 0.2861292977795033), 27L, 1L)
  expect_error(solve(pr, ScipyRoot("krylov")), "not yet ported")
  expect_error(solve(pr, ScipyRoot("anderson"), options = list(jac_options = list(bogus = 1))), "unexpected keyword")
})

test_that("ScipyLinprog marginals of active bounds and ScipyMilp semicontinuous variables follow HiGHS", {
  # python3 tools/scipy-probe.py lp2 (SciPy 1.18.1)
  pr <- OptimizationProblem(function(x, p) c(1, -1, 2, 0.5), numeric(4), numeric(), lb = c(1, -2, -Inf, 0),
                            ub = c(4, 3, Inf, 2))
  s <- solve(pr, ScipyLinprog(), A_ub = rbind(c(1, 1, 1, 1), c(1, -2, 0, 1)), b_ub = c(5, 2),
             A_eq = rbind(c(0, 1, 1, -1)), b_eq = 1)
  expect_equal(unname(s$u), c(1, 3, -2, 0), tolerance = 1e-14)
  expect_equal(s$objective, -6, tolerance = 1e-14)
  o <- s$original
  expect_equal(o$ineqlin$marginals, c(0, 0), tolerance = 1e-14)
  expect_equal(o$eqlin$marginals, 2, tolerance = 1e-14)
  expect_equal(o$lower$marginals, c(1, 0, 0, 2.5), tolerance = 1e-14)
  expect_equal(o$upper$marginals, c(0, -3, 0, 0), tolerance = 1e-14)
  s <- solve(OptimizationProblem(function(x, p) c(1, 1), c(0, 0), numeric()), ScipyLinprog("simplex"),
             A_ub = rbind(c(-1, -2), c(-3, -1)), b_ub = c(-4, -6))
  expect_equal(unname(s$u), c(1.6, 1.2), tolerance = 1e-14)
  expect_identical(s$original$message, "Optimization terminated successfully.")
  s <- solve(OptimizationProblem(function(x, p) c(-1, -2), c(0, 0), numeric(), lb = c(2, 1.5), ub = c(5, 4)),
             ScipyMilp(), integrality = c(2, 3), A = rbind(c(1, 1), c(1, 3)), lb_con = c(-Inf, -Inf),
             ub_con = c(5.5, 10))
  expect_equal(unname(s$u), c(3.5, 2), tolerance = 1e-12)
  expect_equal(s$objective, -7.5, tolerance = 1e-12)
})
