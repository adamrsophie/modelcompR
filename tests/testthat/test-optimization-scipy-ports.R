# The native ports of SLSQP, TNC, COBYQA, trust-constr, trust-krylov, shgo,
# DIRECT, root(method = "krylov") and least_squares(tr_solver = "lsmr") of
# SciPy 1.18.1. Reference values come from tools/scipy-probe.py (the section
# is named in each comment). SciPy's compiled kernels link Accelerate and
# fuse multiply adds, so on the reference platform (ref_platform(), arm64
# macOS) the points agree to rounding and the counts exactly; elsewhere the
# last bits differ and only the converged results are checked.

spp_fq <- function(x, p) c(x[1] + 0.1 * x[1] * x[1] * x[1] - 0.2 * x[2] - 1,
                           x[2] + 0.1 * x[2] * x[2] * x[2] + 0.2 * x[1] - 0.5, x[3] + 0.05 * x[1] * x[3] - 0.3)
spp_res3 <- function(x, p) c(10 * (x[2] - x[1] * x[1]), 1 - x[1], 0.5 * (x[1] * x[2] - 1))
spp_sys3 <- function(x, p) c(x[1] * x[1] + x[2] * x[2] - 4, x[1] * x[2] - 1, x[3] * x[3] * x[3] + x[1] - 2)

test_that("root(method = 'krylov') follows the KrylovJacobian of SciPy", {
  # python3 tools/scipy-probe.py krylov; the Jacobian vector products are
  # forward differences, which amplify rounding, so x agrees to 1e-10
  chk <- function(s, x, nit, status) {
    expect_equal(unname(s$u), x, tolerance = 1e-6)
    if (ref_platform()) {
      expect_equal(unname(s$u), x, tolerance = 1e-10)
      expect_identical(s$original$nit, nit)
      expect_identical(s$original$status, status)
    }
  }
  pr <- OptimizationProblem(spp_fq, c(0, 0, 0), numeric())
  chk(solve(pr, ScipyRoot("krylov")), c(0.9695276947461309, 0.3033042689763239, 0.28612948013707656), 4L, 1L)
  chk(solve(OptimizationProblem(spp_fq, c(2, -1, 3), numeric()), ScipyRoot("krylov"),
            options = list(line_search = NULL)),
      c(0.9695270210708323, 0.3033043909157576, 0.28612948656265985), 5L, 1L)
  sys2 <- function(x, p) { d <- x[1] - x[2]; c(x[1] + 0.5 * d * d * d - 1, -0.5 * d * d * d + x[2]) }
  chk(solve(OptimizationProblem(sys2, c(0.3, 0.1), numeric()), ScipyRoot("krylov")),
      c(0.8411644184209905, 0.15883558158391173), 5L, 1L)
  chk(solve(pr, ScipyRoot("krylov"), options = list(jac_options = list(outer_k = 2, inner_inner_m = 2))),
      c(0.9695277744535771, 0.30330425477072726, 0.2861294752185819), 4L, 1L)
  # a singular Jacobian at the start gives a zero Newton step, as in SciPy
  expect_error(solve(OptimizationProblem(spp_sys3, c(1, 1, 1), numeric()), ScipyRoot("krylov")),
               "Jacobian inversion yielded zero vector")
  expect_error(solve(pr, ScipyRoot("krylov"), options = list(jac_options = list(bogus = 1))), "Unknown parameter")
})

test_that("least_squares with tr_solver = 'lsmr' follows SciPy", {
  # python3 tools/scipy-probe.py krylov
  tt <- seq(0, 2, length.out = 9)
  yy <- 1 + 0.5 * tt - 0.3 * tt * tt + c(0.1, -0.2, 0.05, 3.0, -0.1, 0.0, 0.2, -2.5, 0.1)
  fit <- function(x, p) x[1] + x[2] * tt + x[3] * tt * tt * x[1] - yy
  chk <- function(s, x, nfev, njev, status, tol = 1e-10) {
    expect_equal(unname(s$u), x, tolerance = 1e-5)
    if (ref_platform()) {
      expect_equal(unname(s$u), x, tolerance = tol)
      expect_identical(s$original$nfev, nfev)
      expect_identical(s$original$njev, njev)
      expect_identical(s$original$status, status)
    }
  }
  pr <- OptimizationProblem(spp_res3, c(-1.2, 1), numeric())
  prb <- OptimizationProblem(spp_res3, c(-1.2, 1), numeric(), lb = c(-2, -2), ub = c(0.8, 2))
  ps <- OptimizationProblem(spp_sys3, c(3, -2, 5), numeric())
  pf <- OptimizationProblem(fit, c(1, 0.7, 0), numeric(), lb = c(0, 0.6, -1), ub = c(2, 2, 1))
  L <- "lsmr"
  chk(solve(pr, ScipyLeastSquaresTRF(), tr_solver = L), c(1.0000000000000213, 0.9999999999997223), 28L, 20L, 1L)
  chk(solve(prb, ScipyLeastSquaresTRF(), tr_solver = L), c(0.7999999999999999, 0.6409744407630403), 26L, 21L, 2L)
  chk(solve(ps, ScipyLeastSquaresTRF(), tr_solver = L), c(1.9318516525781366, 0.5176380902050415, 0.4084621098137116),
      11L, 11L, 1L)
  chk(solve(pr, ScipyLeastSquaresDogbox(), tr_solver = L), c(0.999999999999921, 0.9999999999872115), 20L, 16L, 1L)
  chk(solve(prb, ScipyLeastSquaresDogbox(), tr_solver = L), c(0.8, 0.640974440890955), 18L, 14L, 1L)
  chk(solve(ps, ScipyLeastSquaresDogbox(), tr_solver = L),
      c(1.9318516525781366, 0.5176380902050416, 0.4084621098161344), 11L, 11L, 1L)
  chk(solve(pr, ScipyLeastSquaresTRF(), tr_solver = L, tr_options = list(regularize = FALSE)),
      c(1.0000000000000213, 0.9999999999997223), 28L, 20L, 1L)
  # stops by ftol on flat minima: the points agree to 1e-6 as with the exact
  # solver (see test-optimization-scipy.R)
  chk(solve(pf, ScipyLeastSquares(loss = "soft_l1"), tr_solver = L),
      c(1.0217180353360171, 1.1621145550054288, -0.7612408765607316), 9L, 9L, 2L, tol = 1e-6)
  chk(solve(pf, ScipyLeastSquares(method = "dogbox", loss = "soft_l1"), tr_solver = L),
      c(1.021721757030659, 1.1620953419690758, -0.7612279002206009), 8L, 8L, 2L, tol = 1e-6)
  chk(solve(OptimizationProblem(fit, c(1, 0, 0), numeric()), ScipyLeastSquares(), tr_solver = L,
            tr_options = list(damp = 0.1)),
      c(1.0627581167304667, 2.086758822812002, -1.3299199645359276), 10L, 9L, 2L, tol = 1e-6)
  expect_error(solve(pr, ScipyLeastSquares(), tr_solver = "qr"), "must be None, 'exact' or 'lsmr'")
})

spp_rosen <- function(x) { a <- 1 - x[1]; b <- x[2] - x[1] * x[1]; a * a + 100 * b * b }
spp_rosen_grad <- function(x) { b <- x[2] - x[1] * x[1]; c(-2 * (1 - x[1]) - 400 * x[1] * b, 200 * b) }
spp_rosen_hess <- function(x) matrix(c(2 - 400 * x[2] + 1200 * x[1] * x[1], -400 * x[1], -400 * x[1], 200), 2, 2)
spp_poly3 <- function(x) {
  a <- x[1] - 1; b <- x[2] - 2; c <- x[3] + 1
  a * a + 2 * b * b + 3 * c * c + a * a * b * b + 0.5 * a * c
}
spp_poly3_grad <- function(x) {
  a <- x[1] - 1; b <- x[2] - 2; c <- x[3] + 1
  c(2 * a + 2 * a * b * b + 0.5 * c, 4 * b + 2 * a * a * b, 6 * c + 0.5 * a)
}
spp_poly3_hess <- function(x) {
  a <- x[1] - 1; b <- x[2] - 2
  matrix(c(2 + 2 * b * b, 4 * a * b, 0.5, 4 * a * b, 4 + 2 * a * a, 0, 0.5, 0, 6), 3, 3)
}
# the trajectory (points to `tol`, the counts exactly) on the reference
# platform; elsewhere the converged point to 1e-5 (`converged = FALSE` for runs
# that stop before convergence, which depend on the trajectory)
spp_check <- function(res, ref, tol = 1e-10, label = "", converged = TRUE) {
  if (converged) expect_equal(res$x, ref$x, tolerance = 1e-5, info = label)
  if (ref_platform()) {
    expect_equal(res$x, ref$x, tolerance = tol, info = label)
    expect_lte(abs(res$fun - ref$fun), tol * max(1, abs(ref$fun)), label = paste(label, "fun"))
    for (f in intersect(c("nit", "nfev", "njev", "nhev", "status"), names(ref)))
      expect_identical(as.integer(res[[f]]), ref[[f]], info = paste(label, f))
  }
}

test_that("trust-krylov follows the trlib solver of SciPy", {
  # python3 tools/scipy-probe.py trlib
  chain <- function(x) {
    f <- 0
    for (i in seq_len(length(x) - 1L)) { a <- 1 - x[i]; b <- x[i + 1] - x[i] * x[i]; f <- f + a * a + 100 * b * b }
    f
  }
  chain_grad <- function(x) {
    g <- numeric(length(x))
    for (i in seq_len(length(x) - 1L)) {
      b <- x[i + 1] - x[i] * x[i]
      g[i] <- g[i] + (-2 * (1 - x[i]) - 400 * x[i] * b)
      g[i + 1] <- g[i + 1] + 200 * b
    }
    g
  }
  chain_hess <- function(x) {
    n <- length(x)
    H <- matrix(0, n, n)
    for (i in seq_len(n - 1L)) {
      H[i, i] <- H[i, i] + (2 - 400 * x[i + 1] + 1200 * x[i] * x[i])
      H[i, i + 1] <- H[i, i + 1] + (-400 * x[i])
      H[i + 1, i] <- H[i + 1, i] + (-400 * x[i])
      H[i + 1, i + 1] <- H[i + 1, i + 1] + 200
    }
    H
  }
  # a double well in each coordinate: indefinite Hessians at the start
  well <- function(x) {
    f <- 0.1 * x[1]
    for (i in seq_along(x)) { u <- x[i] * x[i] - 1; f <- f + u * u }
    for (i in seq_len(length(x) - 1L)) f <- f + 0.5 * x[i] * x[i + 1]
    f
  }
  well_grad <- function(x) {
    n <- length(x)
    g <- 4 * x * (x * x - 1)
    g[1] <- g[1] + 0.1
    g[-n] <- g[-n] + 0.5 * x[-1]
    g[-1] <- g[-1] + 0.5 * x[-n]
    g
  }
  well_hess <- function(x) {
    n <- length(x)
    H <- diag(12 * x * x - 4, n)
    for (i in seq_len(n - 1L)) { H[i, i + 1] <- 0.5; H[i + 1, i] <- 0.5 }
    H
  }
  tro <- list(gtol = 1e-10, maxiter = 50000)
  trx <- list(gtol = 1e-10, maxiter = 50000, inexact = FALSE)
  c0 <- c(-1.2, 1, -1.2, 1, -1.2, 1)
  w0 <- c(0.1, -0.1, 0.1, -0.1, 0.1)
  km <- "trust-krylov"
  spp_check(.sp_minimize(spp_poly3, c(3, -1, 0.5), km, jac = spp_poly3_grad, hess = spp_poly3_hess, options = trx),
            list(x = c(0.9999999999931417, 1.9999999999993214, -0.9999999999994285), fun = 4.697715501359636e-23,
                 nit = 6L, nfev = 7L, njev = 7L, nhev = 6L, status = 0L), label = "poly3 exact")
  spp_check(.sp_minimize(chain, c0, km, jac = chain_grad, hess = chain_hess, options = tro),
            list(x = c(0.9999999999999987, 0.9999999999999972, 0.9999999999999936, 0.9999999999999843,
                       0.9999999999999574, 0.9999999999998689), fun = 2.260304662802299e-25,
                 nit = 57L, nfev = 58L, njev = 58L, nhev = 50L, status = 0L), label = "chain")
  spp_check(.sp_minimize(chain, c0, km, jac = chain_grad, hess = chain_hess, options = trx),
            list(x = c(-0.9865749795427969, 0.9833982287783201, 0.9721066699302146, 0.9474374365520487,
                       0.8986511841975399, 0.8075739503763605), fun = 3.973940500930295,
                 nit = 20L, nfev = 22L, njev = 22L, nhev = 19L, status = 2L), label = "chain exact",
            converged = FALSE)
  spp_check(.sp_minimize(well, w0, km, jac = well_grad, hess = well_hess, options = trx),
            list(x = c(1.052916435817234, -1.1150519417519669, 1.1177849094642143, -1.1155379560170513,
                       1.0635373750188475), fun = -2.1114330390557465,
                 nit = 8L, nfev = 10L, njev = 10L, nhev = 8L, status = 2L), label = "well exact", converged = FALSE)
  spp_check(.sp_minimize(spp_rosen, c(-1.2, 1), km, jac = spp_rosen_grad, hess = spp_rosen_hess),
            list(x = c(0.9999303920730152, 0.9998605099754639), fun = 4.853048482460746e-09,
                 nit = 36L, nfev = 37L, njev = 37L, nhev = 32L, status = 0L), label = "rosen")
  rosen <- function(u, p) { a <- p[1] - u[1]; b <- u[2] - u[1] * u[1]; a * a + p[2] * b * b }
  rg <- function(u, p) { b <- u[2] - u[1] * u[1]; c(-2 * (p[1] - u[1]) - 4 * p[2] * u[1] * b, 2 * p[2] * b) }
  rh <- function(u, p) matrix(c(2 - 4 * p[2] * u[2] + 12 * p[2] * u[1] * u[1], -4 * p[2] * u[1], -4 * p[2] * u[1],
                                2 * p[2]), 2, 2)
  # with gtol = 1e-10 the inexact steps stagnate near the minimum (SciPy runs
  # to maxiter as in the trust section of the probe), so cap the iterations
  sol <- solve(OptimizationProblem(rosen, c(-1.2, 1), c(1, 100), grad = rg, hess = rh), ScipyTrustKrylov(),
               maxiters = 60L)
  expect_equal(unname(sol$u), c(1, 1), tolerance = 1e-6)
})

test_that("SLSQP follows slsqp.c of SciPy", {
  # python3 tools/scipy-probe.py slsqp; forward difference derivatives (the
  # fd and mixed cases) amplify rounding, so x agrees to 1e-7 there
  sm <- "SLSQP"
  rg <- spp_rosen_grad
  pg <- spp_poly3_grad
  nlc <- .sp_nonlinear_constraint
  ceq_lin <- list(type = "eq", fun = function(x) x[1] + x[2] + x[3] - 1, jac = function(x) matrix(c(1, 1, 1), 1))
  q2 <- function(x) { a <- x[1] - 2; b <- x[2] - 1; a * a + b * b }
  q2g <- function(x) c(2 * (x[1] - 2), 2 * (x[2] - 1))
  spp_check(.sp_minimize(spp_rosen, c(-1.2, 1), sm, jac = rg),
            list(x = c(0.9999219216952796, 0.999851010199817), fun = 1.1223802786679747e-08,
                 nit = 34L, nfev = 47L, njev = 34L, status = 0L), label = "rosen")
  spp_check(.sp_minimize(spp_rosen, c(-1.2, 1), sm),
            list(x = c(0.9999176162751068, 0.9998424625361283), fun = 1.2004538270608679e-08,
                 nit = 34L, nfev = 115L, njev = 34L, status = 0L), tol = 1e-7, label = "rosen fd")
  spp_check(.sp_minimize(spp_poly3, c(3, -1, 0.5), sm),
            list(x = c(0.9997359702214009, 1.9999099710361217, -1.0001251180542656), fun = 1.4940318181901482e-07,
                 nit = 10L, nfev = 43L, njev = 10L, status = 0L), tol = 1e-9, label = "poly3 fd")
  spp_check(.sp_minimize(spp_rosen, c(-1.2, 1), sm, jac = rg, options = list(maxiter = 5)),
            list(x = c(-0.4852625613017814, 0.20704089293198674), fun = 2.286881754483268,
                 nit = 5L, nfev = 12L, njev = 6L, status = 9L), label = "maxiter", converged = FALSE)
  spp_check(.sp_minimize(spp_poly3, c(0.5, 1, 0), sm, jac = pg,
                         bounds = list(lb = c(0, -Inf, -0.5), ub = c(0.8, 1.5, Inf))),
            list(x = c(0.7999999999999927, 1.4999999999999805, -0.49999999999997335), fun = 1.250000000000119,
                 nit = 2L, nfev = 2L, njev = 2L, status = 0L), label = "bounds")
  spp_check(.sp_minimize(spp_poly3, c(3, -1, 0.5), sm, jac = pg, constraints = list(ceq_lin)),
            list(x = c(0.4442497577819467, 1.7229323824717433, -1.1671821402536902), fun = 0.6164065848350512,
                 nit = 8L, nfev = 10L, njev = 8L, status = 0L), label = "eq")
  sph <- function(x) x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 2
  spp_check(.sp_minimize(spp_poly3, c(1, 0.5, 0), sm, jac = pg,
                         bounds = list(lb = c(-Inf, -1, -Inf), ub = c(Inf, 0.9, Inf)),
                         constraints = list(list(type = "eq", fun = sph),
                                            list(type = "ineq", fun = function(x) x[1] - x[3] - 0.2))),
            list(x = c(0.7453760613974816, 0.8999999999999995, -0.7965015012589937), fun = 2.6616088260868693,
                 nit = 7L, nfev = 8L, njev = 7L, status = 0L), tol = 1e-7, label = "mixed")
  # dependent equality constraints: dgelsy finds rank 1 < 2
  dep <- .sp_minimize(spp_poly3, c(3, -1, 0.5), sm, jac = pg,
                      constraints = list(ceq_lin, list(type = "eq", fun = function(x) 2 * (x[1] + x[2] + x[3] - 1),
                                                       jac = function(x) matrix(c(2, 2, 2), 1))))
  spp_check(dep, list(x = c(3, -1, 0.5), fun = 66.25, nit = 1L, nfev = 1L, njev = 1L, status = 6L),
            label = "dependent")
  expect_identical(dep$message, "Singular matrix C in LSQ subproblem")
  # linearizations inconsistent at the start: the augmented problem
  c_two <- function(x) c(x[1] * x[1] + x[2] * x[2] - 2, x[1] - x[2])
  c_two_jac <- function(x) matrix(c(2 * x[1], 1, 2 * x[2], -1), 2)
  spp_check(.sp_minimize(q2, c(0, 0), sm, jac = q2g, constraints = list(list(type = "eq", fun = c_two, jac = c_two_jac))),
            list(x = c(1.0000000004573606, 1.0000000004573606), fun = 0.9999999990852788,
                 nit = 7L, nfev = 9L, njev = 7L, status = 0L), label = "square")
  spp_check(.sp_minimize(function(x) (x[1] - 0.5) * (x[1] - 0.5) + (x[2] - 1) * (x[2] - 1), c(0.5, 0.5), sm,
                         jac = function(x) c(2 * (x[1] - 0.5), 2 * (x[2] - 1)),
                         bounds = list(lb = c(0, 0), ub = c(1, 2.5)),
                         constraints = list(list(type = "ineq", fun = function(x) x[1] * x[1] + x[2] * x[2] - 4,
                                                 jac = function(x) matrix(c(2 * x[1], 2 * x[2]), 1)))),
            list(x = c(0.8944272167217433, 1.7888543692486438), fun = 0.7778640451735367,
                 nit = 7L, nfev = 7L, njev = 7L, status = 0L), label = "inconsistent")
  # NonlinearConstraint, with a variable fixed by its bounds (removed because
  # the constraint Jacobian is a finite difference)
  two <- nlc(function(x) c(x[1] * x[2] + x[3], x[1] - x[2]), c(-Inf, -0.25), c(1, -0.25))
  spp_check(.sp_minimize(spp_poly3, c(0.5, 1, -0.5), sm, jac = pg,
                         bounds = list(lb = c(0, -Inf, -0.5), ub = c(2, Inf, -0.5)), constraints = list(two)),
            list(x = c(1.1061072254371762, 1.3561072254371762, -0.5), fun = 1.6216492113001604,
                 nit = 5L, nfev = 5L, njev = 5L, status = 0L), label = "fixed")
  ref_two <- list(x = c(0.8910057222281111, 1.9749982766962093, -1.0097347668292682), fun = 0.013952164787935163,
                  nit = 7L, nfev = 9L, njev = 7L, status = 0L)
  spp_check(.sp_minimize(spp_poly3, c(0.5, 1, 0), sm, jac = pg,
                         constraints = list(nlc(function(x) x[1] * x[2] + x[3], -0.5, 0.75))),
            ref_two, label = "two sided")
  # every variable fixed by its bounds, at a point that violates the constraint
  fx <- .sp_minimize(spp_poly3, c(0, 0, 0), sm, bounds = list(lb = c(1, 2, 0.5), ub = c(1, 2, 0.5)),
                     constraints = list(two))
  expect_identical(fx[c("x", "fun", "nfev", "njev", "success")],
                   list(x = c(1, 2, 0.5), fun = 6.75, nfev = 1L, njev = 0L, success = FALSE))
  expect_identical(fx$message, paste("All independent variables were fixed by bounds, but the independent variables",
                                     "do not satisfy the constraints exactly. (Maximum violation: 1.5)."))
  # through solve(): the constraint of the problem is the two sided
  # NonlinearConstraint above (its finite_diff_rel_step draws the SciPy warning)
  pr <- OptimizationProblem(function(u, p) spp_poly3(u), c(0.5, 1, 0), numeric(),
                            grad = function(u, p) spp_poly3_grad(u),
                            cons = function(u, p) u[1] * u[2] + u[3], lcons = -0.5, ucons = 0.75)
  expect_warning(s <- solve(pr, ScipySLSQP()), "are ignored by this method")
  expect_identical(s$retcode, "Success")
  spp_check(s$original, ref_two, label = "solve")
  expect_warning(.sp_minimize(spp_rosen, c(-1.2, 1), sm, jac = rg, options = list(bogus = 1)),
                 "Unknown solver options: bogus")
})

test_that("TNC follows tnc.c of SciPy", {
  # python3 tools/scipy-probe.py tnc
  tm <- "TNC"
  rg <- spp_rosen_grad
  pg <- spp_poly3_grad
  chain <- function(x) {
    f <- 0
    for (i in seq_len(length(x) - 1L)) {
      a <- 1 - x[i]
      b <- x[i + 1] - x[i] * x[i]
      f <- f + (a * a + 100 * b * b)
    }
    f
  }
  chain_grad <- function(x) {
    g <- numeric(length(x))
    for (i in seq_len(length(x) - 1L)) {
      b <- x[i + 1] - x[i] * x[i]
      g[i] <- g[i] + (-2 * (1 - x[i]) - 400 * x[i] * b)
      g[i + 1] <- g[i + 1] + 200 * b
    }
    g
  }
  c0 <- c(-1.2, 1, -1.2, 1, -1.2, 1)
  b3 <- list(lb = c(0, -Inf, -0.5), ub = c(0.8, 1.5, Inf))
  r <- .sp_minimize(spp_rosen, c(-1.2, 1), tm, jac = rg)
  spp_check(r, list(x = c(0.9999942161808862, 0.9999884466230691), fun = 3.347280669520409e-11,
                    nit = 23L, nfev = 80L, status = 1L), label = "rosen")
  expect_identical(r[c("message", "success")], list(message = "Converged (|f_n-f_(n-1)| ~= 0)", success = TRUE))
  spp_check(.sp_minimize(spp_rosen, c(-1.2, 1), tm),
            list(x = c(0.9999969695596286, 0.9999939300050315), fun = 9.191892504372929e-12,
                 nit = 25L, nfev = 297L, status = 1L), label = "rosen fd")
  spp_check(.sp_minimize(spp_poly3, c(3, -1, 0.5), tm),
            list(x = c(0.9999999943497029, 1.9999999997870983, -1.0000000069259247), fun = 1.9548857784227748e-16,
                 nit = 9L, nfev = 100L, status = 0L), label = "poly3 fd")
  spp_check(.sp_minimize(spp_poly3, c(0.5, 1, 0), tm, jac = pg, bounds = b3),
            list(x = c(0.8, 1.5, -0.5), fun = 1.25, nit = 3L, nfev = 7L, status = 0L), label = "bounds")
  spp_check(.sp_minimize(spp_poly3, c(0.5, 1, 0), tm, bounds = b3),
            list(x = c(0.8, 1.5, -0.5), fun = 1.25, nit = 3L, nfev = 28L, status = 0L), label = "bounds fd")
  # a variable fixed by its bounds is removed; its gradient entry is NaN
  fx <- .sp_minimize(spp_poly3, c(0.5, 1, 0), tm, jac = pg, bounds = list(lb = c(0, 1.25, -Inf), ub = c(2, 1.25, Inf)))
  spp_check(fx, list(x = c(1.0000000014680737, 1.25, -0.9999999995281191), fun = 1.125,
                     nit = 2L, nfev = 5L, status = 0L), label = "fixed")
  expect_true(is.nan(fx$jac[2]) && !anyNA(fx$jac[-2]))
  spp_check(.sp_minimize(chain, c0, tm, jac = chain_grad,
                         bounds = list(lb = c(-1, -1, -1, -Inf, -Inf, -Inf), ub = c(0.9, 0.9, 0.9, Inf, Inf, Inf))),
            list(x = c(-0.9432393198066283, 0.8999999999999999, 0.8285012150056523, 0.6965427561359374,
                       0.49027000909588925, 0.24036494379935613), fun = 4.22519729078256,
                 nit = 20L, nfev = 89L, status = 1L), label = "chain")
  spp_check(.sp_minimize(chain, c0, tm, jac = chain_grad, options = list(maxfun = 2000)),
            list(x = c(1.0000004379094263, 1.0000010829097568, 1.0000010377997601, 1.0000000689432178,
                       0.9999997909914137, 0.9999995433209297), fun = 5.488721759953213e-10,
                 nit = 23L, nfev = 109L, status = 1L), label = "chain free")
  c9 <- rep(c(-1.2, 1), length.out = 9L)
  spp_check(.sp_minimize(chain, c9, tm, jac = chain_grad, bounds = list(lb = c(-2, -2, 0.5, rep(-2, 6)), ub = rep(1.5, 9)),
                         options = list(maxfun = 3000)),
            list(x = c(-0.9931933065964428, 0.9964670225577223, 0.9979634020046241, 0.9984340585588997,
                       0.9981159925498893, 0.9968502151083372, 0.9940029702070174, 0.9881591602229495,
                       0.9764577574371524), fun = 3.9864428389452544, nit = 24L, nfev = 126L, status = 1L),
            label = "chain nine")
  mf <- .sp_minimize(spp_rosen, c(-1.2, 1), tm, jac = rg, options = list(maxfun = 20))
  spp_check(mf, list(x = c(-0.7892632179336032, 0.58676874644523), fun = 3.332272976044362,
                     nit = 4L, nfev = 21L, status = 3L), label = "maxfun", converged = FALSE)
  expect_identical(mf[c("message", "success")],
                   list(message = "Max. number of function evaluations reached", success = FALSE))
  spp_check(.sp_minimize(spp_rosen, c(-1.2, 1), tm, jac = rg, options = list(scale = c(2, 0.5), offset = c(0.1, -0.2))),
            list(x = c(0.9999774665833416, 0.9999548200011878), fun = 5.090470270819159e-10,
                 nit = 20L, nfev = 74L, status = 1L), label = "scale")
  spp_check(.sp_minimize(spp_rosen, c(-1.2, 1), tm, jac = rg, options = list(maxCGit = 0, maxfun = 300)),
            list(x = c(0.8021234662693005, 0.6428871667983225), fun = 0.0391816336017036,
                 nit = 38L, nfev = 301L, status = 3L), label = "steepest", converged = FALSE)
  spp_check(.sp_minimize(spp_rosen, c(-1.2, 1), tm, jac = rg,
                         options = list(eta = 0.5, stepmx = 1, accuracy = 1e-10, ftol = 1e-12, xtol = 1e-12,
                                        gtol = 1e-12, rescale = 0, minfev = -1)),
            list(x = c(-0.9452072363660899, 0.906068480169706), fun = 3.799837896762872,
                 nit = 49L, nfev = 101L, status = 3L), label = "options", converged = FALSE)
  spp_check(.sp_minimize(spp_poly3, c(3, 5, -4), tm, jac = pg, bounds = list(lb = c(0, 0, -2), ub = c(2, 2.5, 0))),
            list(x = c(1.000002421762131, 1.9999993154154208, -1.0000010636658057), fun = 8.908425965949463e-12,
                 nit = 6L, nfev = 14L, status = 1L), label = "outside")
  expect_error(.sp_minimize(spp_rosen, c(-1.2, 1), tm, jac = rg, options = list(scale = c(1, 1, 1))),
               "tnc: vector sizes must be equal")
  # through solve(), with the bounds of the problem
  pr <- OptimizationProblem(function(u, p) spp_poly3(u), c(0.5, 1, 0), numeric(),
                            grad = function(u, p) spp_poly3_grad(u), lb = b3$lb, ub = b3$ub)
  s <- solve(pr, ScipyTNC())
  expect_identical(s$retcode, "Success")
  spp_check(s$original, list(x = c(0.8, 1.5, -0.5), fun = 1.25, nit = 3L, nfev = 7L, status = 0L), label = "solve")
  expect_warning(.sp_minimize(spp_rosen, c(-1.2, 1), tm, jac = rg, options = list(maxiter = 5)),
                 "Unknown solver options: maxiter")
})

test_that("DIRECT follows DIRect.c of SciPy", {
  # python3 tools/scipy-probe.py direct: the drf_ cases (SciPy with the log
  # of its C code on, which drops the rectangles off the convex hull as the
  # port does; the dr_ cases of plain direct() keep them, see the porting
  # notes)
  dr <- function(f, lb, ub, ...) {
    a <- utils::modifyList(list(eps = 1e-4, maxiter = 1000L, locally_biased = TRUE, vol_tol = 1e-16, len_tol = 1e-6),
                           list(...))
    do.call(.sp_direct, c(list(f, lb, ub), a))
  }
  l2 <- c(-2, -1)
  u2 <- c(2, 3)
  l3 <- c(-1, 0, -2)
  u3 <- c(2.5, 3, 1)
  vol <- "The volume of the hyperrectangle containing the lowest function value found is below vol_tol="
  len <- "The side length measure of the hyperrectangle containing the lowest function value found is below len_tol="
  r <- dr(spp_rosen, l2, u2)
  spp_check(r, list(x = c(1.0000056450292685, 1.000012544509488), fun = 1.892230794801337e-10,
                    nit = 47L, nfev = 713L, status = 5L), label = "rosen")
  expect_identical(r[c("message", "success")], list(message = paste0(len, "1e-06"), success = TRUE))
  r <- dr(spp_poly3, l3, u3)
  spp_check(r, list(x = c(0.9999929437134136, 1.9999915324560957, -0.9999971774853652), fun = 2.0713131044621974e-10,
                    nit = 33L, nfev = 551L, status = 4L), label = "poly3")
  expect_identical(r$message, paste0(vol, "1e-16"))
  spp_check(dr(spp_rosen, l2, u2, locally_biased = FALSE),
            list(x = c(1.0000508052634243, 1.0001028649777997), fun = 2.7378925852873623e-09,
                 nit = 42L, nfev = 793L, status = 5L), label = "original")
  r <- dr(spp_poly3, l3, u3, locally_biased = FALSE, maxfun = 600)
  spp_check(r, list(x = c(0.9948559670781898, 2.002057613168724, -0.9979423868312757), fun = 4.23378315512262e-05,
                    nit = 16L, nfev = 683L, status = 1L), label = "original maxfun", converged = FALSE)
  expect_identical(r[c("message", "success")],
                   list(message = "Number of function evaluations done is larger than maxfun=600", success = FALSE))
  spp_check(dr(spp_poly3, l3, u3, maxfun = 200),
            list(x = c(0.9996570644718796, 1.9993141289437586, -0.9993141289437586), fun = 2.35209558427236e-06,
                 nit = 19L, nfev = 217L, status = 1L), label = "maxfun", converged = FALSE)
  seen <- list()
  r <- dr(spp_rosen, l2, u2, maxiter = 12, callback = function(x) seen[[length(seen) + 1L]] <<- x)
  spp_check(r, list(x = c(1.333333333333333, 1.7736625514403292), fun = 0.11280461989195328,
                    nit = 12L, nfev = 73L, status = 2L), label = "maxiter", converged = FALSE)
  expect_identical(r$message, "Number of iterations is larger than maxiter=12")
  # the callback sees the last evaluated point after each iteration
  expect_length(seen, 10L)
  expect_equal(seen[[10L]], c(1.333333333333333, 1.7736625514403292), tolerance = 1e-15)
  r <- dr(spp_poly3, l3, u3, vol_tol = 1e-8)
  spp_check(r, list(x = c(0.9948559670781898, 1.9979423868312756, -0.9979423868312757), fun = 4.2337831551229856e-05,
                    nit = 16L, nfev = 167L, status = 4L), tol = 1e-10, label = "vol_tol", converged = FALSE)
  expect_identical(r$message, paste0(vol, "1e-08"))
  r <- dr(spp_rosen, l2, u2, f_min = 0, f_min_rtol = 1e-3)
  spp_check(r, list(x = c(0.9876543209876534, 0.9780521262002744), fun = 0.0008237793536014661,
                    nit = 25L, nfev = 237L, status = 3L), label = "f_min", converged = FALSE)
  expect_identical(r$message, paste("The best function value found is within a relative error=0.001",
                                    "of the (known) global optimum f_min"))
  # a negative eps: the Jones update of epsilon
  r <- dr(spp_poly3, l3, u3, eps = -1e-3, len_tol = 1e-4)
  spp_check(r, list(x = c(1.0000127013158568, 1.999923792104862, -0.999923792104862), fun = 2.968351010579187e-08,
                    nit = 26L, nfev = 367L, status = 5L), label = "eps")
  expect_identical(r$message, paste0(len, "0.0001"))
  qrt <- function(x) x[1] * x[1] * x[1] * x[1] - 3 * x[1] * x[1] * x[1] + 2
  spp_check(dr(qrt, -1, 4), list(x = 2.249970363596335, fun = -6.542968741107202, nit = 97L, nfev = 1003L, status = 1L),
            label = "one variable", converged = FALSE)
  expect_error(.sp_direct(spp_rosen, c(-2, -Inf), u2), "Bounds must not be inf.")
  expect_error(.sp_direct(spp_rosen, c(-2, 3), u2), "Bounds are not consistent min < max")
  expect_error(.sp_direct(spp_rosen, l2, u2, maxfun = 2.5), "maxfun must be of type int.")
  expect_error(.sp_direct(spp_rosen, l2, u2, vol_tol = 2), "vol_tol must be between 0 and 1.")
  # through solve(): the options of OptimizationSciPy.jl
  pr <- OptimizationProblem(function(u, p) spp_rosen(u), c(0, 0), numeric(), lb = l2, ub = u2)
  s <- solve(pr, ScipyDirect())
  expect_identical(s$retcode, "Success")
  spp_check(s$original, list(x = c(1.0000056450292685, 1.000012544509488), fun = 1.892230794801337e-10,
                             nit = 47L, nfev = 713L, status = 5L), label = "solve")
  s <- solve(pr, ScipyDirect(), maxiters = 12)
  expect_identical(s$retcode, "Failure")
  expect_identical(s$original$nfev, 73L)
  expect_error(solve(OptimizationProblem(function(u, p) spp_rosen(u), c(0, 0), numeric()), ScipyDirect()),
               "requires bounds")
})

test_that("shgo follows the simplicial SHGO of SciPy", {
  # python3 tools/scipy-probe.py shgo. The complex (its vertices, their
  # order, the minimizer pool) and every count are those of SciPy. The local
  # minima come from SLSQP with forward differences, whose Accelerate kernels
  # round differently, so on the reference platform the points agree to 1e-8
  # where a local minimization ran and exactly where none did.
  rosen <- spp_rosen
  poly3 <- spp_poly3
  wells <- function(x) { a <- x[1] * x[1] - 1; b <- x[2] * x[2] - 1; a * a + b * b + 0.1 * x[1] + 0.2 * x[2] }
  sym <- function(x) { a <- x[1] * x[1] - 1; b <- x[2] * x[2] - 1; a * a + b * b }
  wells3 <- function(x) {
    a <- x[1] * x[1] - 1; b <- x[2] * x[2] - 1; c3 <- x[3] - 0.5
    a * a + b * b + c3 * c3 * c3 * c3 + 0.1 * x[1] + 0.3 * x[2] * x[3]
  }
  quart1 <- function(x) x[1] * x[1] * x[1] * x[1] - 3 * x[1] * x[1] * x[1] + 2 * x[1] * x[1] + 0.5 * x[1]
  ref <- list()
  ref$sh_rosen <- list(x = c(2.0, 3.0), fun = 101.0, nfev = 5L, nlfev = 0L, nljev = 0L, nit = 1L, success = FALSE)
  ref$sh_rosen$message <- "Failed to find a feasible minimizer point. Lowest sampling point = 101.0"
  ref$sh_poly3 <- list(x = c(0.9999999942760373, 1.9999999924983323, -1.0000000070095656), fun = 3.127770626978238e-16,
                       nfev = 36L, nlfev = 27L, nljev = 6L, nit = 1L, success = TRUE, funl = c(3.127770626978238e-16))
  ref$sh_wells <- list(x = c(-1.012273143356282, -1.0241203228752345), fun = -0.3030578109829779, nfev = 56L,
                       nlfev = 51L, nljev = 15L, nit = 1L, success = TRUE, funl = c(-0.3030578109829779))
  ref$sh_sym4 <- list(x = c(1.0, -1.0), fun = 0.0, nfev = 30L, nlfev = 15L, nljev = 5L, nit = 4L, success = TRUE,
                      funl = c(0.0, 0.0, 0.0, 0.0, 2.0))
  ref$sh_quart5 <- list(x = c(-0.10101176075176221), fun = -0.026903037009023013, nfev = 46L, nlfev = 37L, nljev = 17L,
                        nit = 5L, success = TRUE,
                        funl = c(-0.026903037009023013, 0.18030783319220933, 0.18030783319220944))
  ref$sh_wells3 <- list(x = c(-1.0122729833506914, -1.0330777411691763, 0.9263160692978034), fun = -0.3501504319573555,
                        nfev = 508L, nlfev = 388L, nljev = 89L, nit = 3L, success = TRUE)
  ref$sh_rosen_none <- list(x = c(0.9999999999655058, 0.999999999991398), fun = 3.658411273465293e-19, nfev = 321L,
                            nlfev = 280L, nljev = 82L, nit = 3L, success = TRUE)
  ref$sh_cons <- list(x = c(-1.012273143356282, -1.0241203228752345), fun = -0.3030578109829779, nfev = 76L,
                      nlfev = 69L, nljev = 20L, nit = 2L, success = TRUE,
                      funl = c(-0.3030578109829779, -0.1030734488208584))
  ref$sh_cons_r <- list(x = c(0.9999999999655058, 0.999999999991398), fun = 3.658411273465293e-19, nfev = 169L,
                        nlfev = 164L, nljev = 46L, nit = 3L, success = TRUE,
                        funl = c(3.658411273465293e-19, 2.0116040095907022e-11, 2.4357664392236632e-11))
  ref$sh_eq <- list(x = c(-1.0182475186842106, -1.0182475186842106), fun = -0.30276165106629627, nfev = 42L,
                    nlfev = 38L, nljev = 11L, nit = 1L, success = TRUE, funl = c(-0.30276165106629627))
  ref$sh_infeasible <- list(x = NULL, fun = NULL, nfev = 0L, nlfev = 0L, nljev = 0L, nit = 1L, success = FALSE)
  ref$sh_infeasible$message <- "Failed to find a feasible minimizer point. Lowest sampling point = None"
  ref$sh_maxfev <- list(x = c(-1.0122731344944333, -1.0241203027655008), fun = -0.3030578109829807, nfev = 157L,
                        nlfev = 127L, nljev = 37L, nit = 8L, success = TRUE)
  ref$sh_maxiter <- list(x = c(-1.0122731344944333, -1.0241203027655008), fun = -0.3030578109829807, nfev = 142L,
                         nlfev = 127L, nljev = 37L, nit = 3L, success = TRUE)
  ref$sh_maxev <- list(x = c(-1.0122731344944333, -1.0241203027655008), fun = -0.3030578109829807, nfev = 157L,
                       nlfev = 127L, nljev = 37L, nit = 6L, success = TRUE)
  ref$sh_fmin <- list(x = c(-1.012273143356282, -1.0241203228752345), fun = -0.3030578109829779, nfev = 56L,
                      nlfev = 51L, nljev = 15L, nit = 1L, success = TRUE, funl = c(-0.3030578109829779))
  ref$sh_minhgrd <- list(x = c(0.9872574590865341, -1.0241203283987728), fun = -0.1030734488208584, nfev = 28L,
                         nlfev = 18L, nljev = 5L, nit = 1L, success = TRUE, funl = c(-0.1030734488208584))
  ref$sh_local <- list(x = c(-1.012273143356282, -1.0241203228752345), fun = -0.3030578109829779, nfev = 79L,
                       nlfev = 69L, nljev = 20L, nit = 2L, success = TRUE,
                       funl = c(-0.3030578109829779, -0.1030734488208584))
  ref$sh_every <- list(x = c(-1.0122731344944333, -1.0241203027655008), fun = -0.3030578109829807, nfev = 91L,
                       nlfev = 76L, nljev = 22L, nit = 3L, success = TRUE)
  ref$sh_lbfgsb <- list(x = c(-1.012273126942616, -1.0241203915591932), fun = -0.3030578109829449, nfev = 58L,
                        nlfev = 48L, nljev = 16L, nit = 2L, success = TRUE,
                        funl = c(-0.3030578109829449, -0.10307344881943697))
  ref$sh_symmetry <- list(x = c(0.0, 0.0), fun = 2.0, nfev = 8L, nlfev = 3L, nljev = 1L, nit = 1L, success = TRUE,
                          funl = c(2.0))
  ref$sh_two <- list(x = c(-1.012273143356282, -1.0241203228752345), fun = -0.3030578109829779, nfev = 79L,
                     nlfev = 69L, nljev = 20L, nit = 2L, success = TRUE,
                     funl = c(-0.3030578109829779, -0.1030734488208584))
  ok <- "Optimization terminated successfully."
  chk <- function(res, r, tol = 1e-8, label = "") {
    if (is.null(r$x)) expect_null(res$x, label = label) else expect_equal(res$x, r$x, tolerance = 1e-6, info = label)
    expect_identical(res$success, r$success, info = label)
    expect_identical(res$message, r$message %||% ok, info = label)
    if (ref_platform()) {
      if (!is.null(r$x)) {
        expect_equal(res$x, r$x, tolerance = tol, info = label)
        expect_lte(abs(res$fun - r$fun), tol * max(1, abs(r$fun)), label = paste(label, "fun"))
      }
      for (f in c("nfev", "nlfev", "nljev", "nit")) expect_identical(as.integer(res[[f]]), r[[f]], info = paste(label, f))
      if (!is.null(r$funl)) expect_equal(res$funl, r$funl, tolerance = tol, info = paste(label, "funl"))
    }
  }
  sh <- function(f, lb, ub, ...) {
    a <- utils::modifyList(list(n = 100L, iters = 1L, sampling_method = "simplicial"), list(...), keep.null = TRUE)
    do.call(.sp_shgo, c(list(f, lb, ub), a))
  }
  w2 <- c(-2, -2)
  cons_w <- list(list(type = "ineq", fun = function(x) 1.5 - x[1] - x[2]))
  # no local minimization ran: exact
  chk(sh(rosen, c(-2, -1), c(2, 3)), ref$sh_rosen, 0, "rosen")
  chk(sh(sym, w2, -w2, iters = 4L), ref$sh_sym4, 1e-12, "sym4")
  chk(sh(quart1, -1, 3, iters = 5L), ref$sh_quart5, 1e-12, "quart5")
  chk(sh(poly3, c(-1, 0, -2), c(2.5, 3, 1)), ref$sh_poly3, 1e-12, "poly3")
  chk(sh(wells, w2, -w2), ref$sh_wells, label = "wells")
  chk(sh(wells3, c(-2, -2, -1), c(2, 2, 2), n = 40L, iters = 3L), ref$sh_wells3, label = "wells3")
  chk(sh(rosen, c(-2, -1), c(2, 3), n = NULL, iters = 3L), ref$sh_rosen_none, label = "n NULL")
  chk(sh(wells, w2, -w2, constraints = cons_w, iters = 2L), ref$sh_cons, label = "constraint")
  chk(sh(rosen, c(-2, -1), c(2, 3), iters = 3L,
         constraints = list(list(type = "ineq", fun = function(x) 2 - x[1] * x[1] - x[2] * x[2]))),
      ref$sh_cons_r, label = "disk")
  # an equality constraint is a feasibility test x0 - x1 >= 0 of the sampling
  chk(sh(wells, w2, -w2, constraints = list(list(type = "eq", fun = function(x) x[1] - x[2]))), ref$sh_eq,
      1e-12, "eq")
  chk(sh(wells, w2, -w2, constraints = list(list(type = "ineq", fun = function(x) x[1] + x[2] - 10))),
      ref$sh_infeasible, label = "infeasible")
  chk(sh(wells, w2, -w2, options = list(maxfev = 40)), ref$sh_maxfev, label = "maxfev")
  chk(sh(wells, w2, -w2, options = list(maxiter = 3)), ref$sh_maxiter, label = "maxiter")
  chk(sh(wells, w2, -w2, options = list(maxev = 30)), ref$sh_maxev, label = "maxev")
  # each stopping test warns (Python shows the repeated text once)
  msgs <- character()
  fm <- withCallingHandlers(sh(wells, w2, -w2, options = list(f_min = -0.3, f_tol = 1e-3)),
                           warning = function(w) {
                             msgs <<- c(msgs, conditionMessage(w))
                             invokeRestart("muffleWarning")
                           })
  expect_match(msgs, "A much lower value than expected f\\* = -0.3 was found f_lowest = -0.3030578109",
               all = TRUE)
  chk(fm, ref$sh_fmin, label = "f_min")
  chk(sh(wells, w2, -w2, n = 10L, options = list(minhgrd = 0)), ref$sh_minhgrd, 1e-12, "minhgrd")
  chk(sh(wells, w2, -w2, iters = 2L, options = list(local_iter = 1)), ref$sh_local, label = "local_iter")
  chk(sh(wells, w2, -w2, iters = 3L, options = list(minimize_every_iter = FALSE)), ref$sh_every, label = "every")
  chk(sh(wells, w2, -w2, iters = 2L, minimizer_kwargs = list(method = "L-BFGS-B")), ref$sh_lbfgsb, label = "L-BFGS-B")
  expect_warning(sy <- sh(sym, w2, -w2, options = list(symmetry = TRUE)), "Variable 1 was specified as symmetric")
  chk(sy, ref$sh_symmetry, 1e-12, "symmetry")
  # the callback goes to the local minimizer (SLSQP calls it with x and fun)
  seen <- list()
  chk(sh(wells, w2, -w2, iters = 2L, callback = function(r) seen[[length(seen) + 1L]] <<- r$x), ref$sh_two,
      label = "callback")
  expect_gt(length(seen), 0L)
  if (ref_platform()) {
    expect_length(seen, 20L)
    expect_equal(seen[[20L]], c(0.9872574590865341, -1.0241203283987728), tolerance = 1e-8)
  }
  # the parallel evaluation of the sampling points
  chk(sh(wells, w2, -w2, iters = 2L, workers = 2L), ref$sh_two, label = "workers")
  expect_error(sh(wells, w2, -w2, sampling_method = "sobol"), "Delaunay triangulation of Qhull")
  expect_error(sh(wells, w2, -w2, sampling_method = "grid"), "Unknown sampling_method")
  # through solve(): lcons <= cons(u) <= ucons becomes 1.5 - u1 - u2 >= 0
  pr <- OptimizationProblem(function(u, p) wells(u), c(0, 0), numeric(), lb = w2, ub = -w2,
                            cons = function(u, p) u[1] + u[2], lcons = -Inf, ucons = 1.5)
  s <- solve(pr, ScipyShgo(), iters = 2L)
  expect_identical(s$retcode, "Success")
  chk(s$original, ref$sh_cons, label = "solve")
})

test_that("trust-constr follows minimize_trustregion_constr of SciPy", {
  # python3 tools/scipy-probe.py trustconstr: the options of OptimizationSciPy.jl,
  # the BFGS or SR1 updates, NonlinearConstraint objects (finite_diff_rel_step
  # 1e-8) and LinearConstraint. Equality constraints alone run the SQP method,
  # inequalities and bounds the interior point method. In the cases with another
  # path a last bit difference (dense LU and R's LAPACK against SuperLU and
  # Accelerate) changes a branch, so they are compared at the solution only.
  ref <- list()
  ref$tc_none <- list(x = c(0.9999999999999999, 0.9999999999999996), fun = 4.942706609275402e-30, nit = 63L,
                      nfev = 63L, njev = 63L, nhev = 0L, cg_niter = 114L, status = 1L, constr_nfev = integer(),
                      constr_njev = integer(), nwarn = 0L)
  ref$tc_none$message <- "`gtol` termination condition is satisfied."
  ref$tc_hess <- list(x = c(0.9999999999999901, 0.999999999999968), fun = 1.50120353523075e-26, nit = 44L, nfev = 44L,
                      njev = 30L, nhev = 30L, cg_niter = 79L, status = 1L, constr_nfev = integer(),
                      constr_njev = integer(), nwarn = 0L)
  ref$tc_hess$message <- "`gtol` termination condition is satisfied."
  ref$tc_sr1 <- list(x = c(0.999999999999994, 0.9999999999999787), fun = 8.733133955055788e-27, nit = 81L, nfev = 81L,
                     njev = 81L, nhev = 0L, cg_niter = 141L, status = 1L, constr_nfev = integer(),
                     constr_njev = integer(), nwarn = 0L)
  ref$tc_sr1$message <- "`gtol` termination condition is satisfied."
  ref$tc_fd <- list(x = c(0.9999968597880621, 0.9999937087091308), fun = 9.872761610800063e-12, nit = 68L,
                    nfev = 204L, njev = 68L, nhev = 0L, cg_niter = 116L, status = 2L, constr_nfev = integer(),
                    constr_njev = integer(), nwarn = 0L)
  ref$tc_fd$message <- "`xtol` termination condition is satisfied."
  ref$tc_maxiter <- list(x = c(-1.0737422042011235, 1.0702674978984783), fun = 4.983588709046221, nit = 5L, nfev = 5L,
                         njev = 5L, nhev = 0L, cg_niter = 5L, status = 0L, constr_nfev = integer(),
                         constr_njev = integer(), nwarn = 0L)
  ref$tc_maxiter$message <- "The maximum number of function evaluations is exceeded."
  ref$tc_bounds <- list(x = c(0.894755849106932, 0.799999912941749), fun = 0.011110919407113487, nit = 58L,
                        nfev = 50L, njev = 50L, nhev = 0L, cg_niter = 80L, status = 1L, constr_nfev = c(0L),
                        constr_njev = c(0L), nwarn = 0L)
  ref$tc_bounds$message <- "`gtol` termination condition is satisfied."
  ref$tc_ineq <- list(x = c(0.9072339053843999, 0.8227553560949962), fun = 0.008615660900015936, nit = 63L,
                      nfev = 69L, njev = 69L, nhev = 0L, cg_niter = 91L, status = 1L, constr_nfev = c(69L),
                      constr_njev = c(69L), nwarn = 0L)
  ref$tc_ineq$message <- "`gtol` termination condition is satisfied."
  ref$tc_ineq_fd <- list(x = c(0.907233958306085, 0.8227554523057607), fun = 0.008615651069508626, nit = 69L,
                         nfev = 71L, njev = 71L, nhev = 0L, cg_niter = 100L, status = 1L, constr_nfev = c(213L),
                         constr_njev = c(0L), nwarn = 2L)
  ref$tc_ineq_fd$message <- "`gtol` termination condition is satisfied."
  ref$tc_greater <- list(x = c(-0.6030197635496816, 0.36927925220338564), fun = 2.572860564792901, nit = 96L,
                         nfev = 75L, njev = 75L, nhev = 0L, cg_niter = 124L, status = 2L, constr_nfev = c(75L),
                         constr_njev = c(75L), nwarn = 0L)
  ref$tc_greater$message <- "`xtol` termination condition is satisfied."
  ref$tc_eq <- list(x = c(-0.9048735015447967, 0.8253508019030809), fun = 3.6328395291597175, nit = 14L, nfev = 15L,
                    njev = 15L, nhev = 0L, cg_niter = 13L, status = 1L, constr_nfev = c(15L), constr_njev = c(15L),
                    nwarn = 0L)
  ref$tc_eq$message <- "`gtol` termination condition is satisfied."
  ref$tc_both <- list(x = c(0.8947546748933857, 0.7999978261756493), fun = 0.011111164861026009, nit = 65L,
                      nfev = 62L, njev = 62L, nhev = 0L, cg_niter = 102L, status = 1L, constr_nfev = c(62L, 0L),
                      constr_njev = c(62L, 0L), nwarn = 0L)
  ref$tc_both$message <- "`gtol` termination condition is satisfied."
  ref$tc_mixed <- list(x = c(0.8660161260737592, 0.8660161260737592), fun = 1.364298716658835, nit = 73L, nfev = 70L,
                       njev = 70L, nhev = 0L, cg_niter = 68L, status = 1L, constr_nfev = c(70L), constr_njev = c(70L),
                       nwarn = 0L)
  ref$tc_mixed$message <- "`gtol` termination condition is satisfied."
  ref$tc_mixed2 <- list(x = c(-0.6030197631206492, 0.3692792514531478), fun = 2.5728605631545, nit = 79L, nfev = 52L,
                        njev = 52L, nhev = 0L, cg_niter = 100L, status = 2L, constr_nfev = c(52L),
                        constr_njev = c(52L), nwarn = 0L)
  ref$tc_mixed2$message <- "`xtol` termination condition is satisfied."
  ref$tc_keep <- list(x = c(0.9072339053844001, 0.8227553560949962), fun = 0.008615660900015915, nit = 46L,
                      nfev = 35L, njev = 35L, nhev = 0L, cg_niter = 55L, status = 1L, constr_nfev = c(36L),
                      constr_njev = c(36L), nwarn = 0L)
  ref$tc_keep$message <- "`gtol` termination condition is satisfied."
  ref$tc_svd <- list(x = c(0.9072339053843999, 0.8227553560949962), fun = 0.008615660900015936, nit = 63L, nfev = 69L,
                     njev = 69L, nhev = 0L, cg_niter = 91L, status = 1L, constr_nfev = c(69L), constr_njev = c(69L),
                     nwarn = 0L)
  ref$tc_svd$message <- "`gtol` termination condition is satisfied."
  ref$tc_singular <- list(x = c(-0.6, 1.5999999999999999), fun = 156.31999999999994, nit = 13L, nfev = 5L, njev = 5L,
                          nhev = 0L, cg_niter = 0L, status = 2L, constr_nfev = c(5L), constr_njev = c(5L),
                          nwarn = 11L)
  ref$tc_singular$message <- "`xtol` termination condition is satisfied."
  ref$tc_linear <- list(x = c(0.6187955520987578, 0.3812042976339819), fun = 0.1456070692284466, nit = 62L,
                        nfev = 59L, njev = 59L, nhev = 0L, cg_niter = 93L, status = 1L, constr_nfev = c(0L),
                        constr_njev = c(0L), nwarn = 0L)
  ref$tc_linear$message <- "`gtol` termination condition is satisfied."
  ref$tc_three <- list(x = c(-0.9999999995392004, 0.6666666664850571, 1.333333333054143), fun = 7.666666667076266,
                       nit = 29L, nfev = 17L, njev = 17L, nhev = 0L, cg_niter = 31L, status = 1L,
                       constr_nfev = c(17L, 0L), constr_njev = c(17L, 0L), nwarn = 5L)
  ref$tc_three$message <- "`gtol` termination condition is satisfied."
  ref$tc_three_eq <- list(x = c(0.3562069497862754, 0.5614713584897776, 0.08232169172394717),
                          fun = 11.013145192878136, nit = 10L, nfev = 10L, njev = 10L, nhev = 0L, cg_niter = 9L,
                          status = 1L, constr_nfev = c(10L), constr_njev = c(10L), nwarn = 0L)
  ref$tc_three_eq$message <- "`gtol` termination condition is satisfied."
  ref$tc_hess_fd <- list(x = c(0.9999999999999972, 0.9999999999999917), fun = 7.780756975324433e-28, nit = 44L,
                         nfev = 44L, njev = 170L, nhev = 0L, cg_niter = 79L, status = 1L, constr_nfev = integer(),
                         constr_njev = integer(), nwarn = 0L)
  ref$tc_hess_fd$message <- "`gtol` termination condition is satisfied."
  ref$tc_hess_fd3 <- list(x = c(0.9999999999999977, 0.9999999999999915), fun = 1.5153648210746314e-27, nit = 44L,
                          nfev = 44L, njev = 310L, nhev = 0L, cg_niter = 79L, status = 1L, constr_nfev = integer(),
                          constr_njev = integer(), nwarn = 0L)
  ref$tc_hess_fd3$message <- "`gtol` termination condition is satisfied."
  ref$tc_hessp <- list(x = c(0.9999999999999901, 0.999999999999968), fun = 1.50120353523075e-26, nit = 44L,
                       nfev = 44L, njev = 30L, nhev = 30L, cg_niter = 79L, status = 1L, constr_nfev = integer(),
                       constr_njev = integer(), nwarn = 0L)
  ref$tc_hessp$message <- "`gtol` termination condition is satisfied."
  ref$tc_con_hess <- list(x = c(0.9072339053843147, 0.8227553560950847), fun = 0.008615660900016277, nit = 62L,
                          nfev = 61L, njev = 44L, nhev = 44L, cg_niter = 84L, status = 1L, constr_nfev = c(61L),
                          constr_njev = c(44L), nwarn = 0L)
  ref$tc_con_hess$message <- "`gtol` termination condition is satisfied."
  ref$tc_con_hess_fd <- list(x = c(0.9072339053844, 0.8227553560949961), fun = 0.008615660900015929, nit = 71L,
                             nfev = 71L, njev = 71L, nhev = 0L, cg_niter = 108L, status = 1L, constr_nfev = c(71L),
                             constr_njev = c(331L), nwarn = 0L)
  ref$tc_con_hess_fd$message <- "`gtol` termination condition is satisfied."
  ref$tc_stop <- list(x = c(-1.2, 1.0), fun = 24.199999999999996, nit = 3L, nfev = 3L, njev = 3L, nhev = 0L,
                      cg_niter = 2L, status = 3L, constr_nfev = integer(), constr_njev = integer(), nwarn = 0L)
  ref$tc_stop$message <- "`callback` raised `StopIteration`."
  tcx0 <- c(-1.2, 1)
  tcbox <- list(lb = c(-2, -2), ub = c(2, 0.8))
  circ <- function(x) x[1] * x[1] + x[2] * x[2]
  circ_jac <- function(x) matrix(c(2 * x[1], 2 * x[2]), 1)
  circ_hess <- function(x, v) matrix(c(2 * v[1], 0, 0, 2 * v[1]), 2, 2)
  two <- function(x) c(x[1] * x[1] + x[2] * x[2], x[1] - x[2])
  two_jac <- function(x) matrix(c(2 * x[1], 1, 2 * x[2], -1), 2)
  dup <- function(x) c(x[1] + x[2], x[1] + x[2])
  dup_jac <- function(x) matrix(1, 2, 2)
  f3 <- function(x) {
    a <- x[1] - 1; b <- x[2] - 2; cc <- x[3] - 3
    a * a + b * b + cc * cc + x[1] * x[2] * x[3]
  }
  f3_grad <- function(x) c(2 * (x[1] - 1) + x[2] * x[3], 2 * (x[2] - 2) + x[1] * x[3], 2 * (x[3] - 3) + x[1] * x[2])
  g3 <- function(x) c(x[1] + x[2] + x[3], x[1] * x[2])
  g3_jac <- function(x) matrix(c(1, x[2], 1, x[1], 1, 0), 2)
  rhessp <- function(x, p) {
    h11 <- 2 - 400 * x[2] + 1200 * x[1] * x[1]
    h12 <- -400 * x[1]
    c(h11 * p[1] + h12 * p[2], h12 * p[1] + 200 * p[2])
  }
  rg <- spp_rosen_grad
  nlc <- function(fun, lb, ub, jac = "2-point", keep_feasible = FALSE, hess = NULL)
    .sp_nonlinear_constraint(fun, lb, ub, jac = jac, keep_feasible = keep_feasible, finite_diff_rel_step = 1e-8,
                             hess = hess)
  tc <- function(fun, x0, jac = NULL, hess = "BFGS", bounds = NULL, constraints = list(), options = list()) {
    o <- utils::modifyList(list(initial_tr_radius = 1, verbose = 0, finite_diff_rel_step = 1e-8, gtol = 1e-10,
                                maxiter = 50000L), options)
    if (is.character(hess) && hess %in% c("BFGS", "SR1")) hess <- .sp_hessian_update(hess)
    nwarn <- 0L
    r <- withCallingHandlers(
      .sp_minimize(fun, x0, method = "trust-constr", jac = jac, hess = hess, bounds = bounds,
                   constraints = constraints, options = o),
      warning = function(w) {
        nwarn <<- nwarn + 1L
        invokeRestart("muffleWarning")
      })
    r$nwarn <- nwarn
    r
  }
  ints <- c("nit", "nfev", "njev", "nhev", "cg_niter", "status", "nwarn", "constr_nfev", "constr_njev")
  chk <- function(r, rf, label, tol = 1e-12) {
    expect_equal(r$x, rf$x, tolerance = 1e-6, info = label)
    if (ref_platform()) {
      expect_equal(r$x, rf$x, tolerance = tol, info = label)
      expect_lte(abs(r$fun - rf$fun), tol * max(1, abs(rf$fun)), label = paste(label, "fun"))
      for (f in ints) expect_identical(as.integer(r[[f]]), as.integer(rf[[f]]), info = paste(label, f))
      expect_identical(r$message, rf$message, info = label)
    }
  }
  # the solution of the cases with another path
  near <- function(r, rf, label, tol) {
    expect_equal(r$x, rf$x, tolerance = 1e-6, info = label)
    if (ref_platform()) expect_equal(r$x, rf$x, tolerance = tol, info = label)
  }
  chk(tc(spp_rosen, tcx0, rg), ref$tc_none, "none")
  chk(tc(spp_rosen, tcx0, rg, spp_rosen_hess), ref$tc_hess, "hess")
  chk(tc(spp_rosen, tcx0, rg, "SR1"), ref$tc_sr1, "SR1")
  chk(tc(spp_rosen, tcx0), ref$tc_fd, "finite difference gradient")
  chk(tc(spp_rosen, tcx0, rg, options = list(maxiter = 5L)), ref$tc_maxiter, "maxiter")
  chk(tc(spp_rosen, tcx0, rg, constraints = list(nlc(circ, -Inf, 1.5))), ref$tc_ineq_fd, "constraint differences")
  chk(tc(spp_rosen, tcx0, rg, constraints = list(nlc(circ, 1.5, 1.5, circ_jac))), ref$tc_eq, "equality")
  chk(tc(spp_rosen, tcx0, rg, bounds = tcbox, constraints = list(nlc(circ, -Inf, 1.5, circ_jac))), ref$tc_both,
      "bounds and constraint")
  chk(tc(spp_rosen, tcx0, rg, constraints = list(nlc(two, c(0.5, 0), c(1.5, 0), two_jac))), ref$tc_mixed,
      "interval and equality")
  chk(tc(spp_rosen, c(0.5, 0.5), rg, constraints = list(nlc(circ, -Inf, 1.5, circ_jac, keep_feasible = TRUE))),
      ref$tc_keep, "keep_feasible")
  chk(tc(spp_rosen, tcx0, rg, constraints = list(nlc(dup, 1, 1, dup_jac))), ref$tc_singular, "singular Jacobian")
  chk(tc(f3, c(0.5, 0.5, 0.5), f3_grad, constraints = list(nlc(g3, c(1, 0.2), c(1, 0.2), g3_jac))), ref$tc_three_eq,
      "three variables, two equalities")
  chk(tc(spp_rosen, tcx0, rg, "2-point"), ref$tc_hess_fd, "2-point Hessian")
  chk(tc(spp_rosen, tcx0, rg, "3-point"), ref$tc_hess_fd3, "3-point Hessian")
  chk(tc(spp_rosen, tcx0, rg, NULL, options = list(hessp = rhessp)), ref$tc_hessp, "hessp")
  chk(tc(spp_rosen, tcx0, rg, spp_rosen_hess, constraints = list(nlc(circ, -Inf, 1.5, circ_jac, hess = circ_hess))),
      ref$tc_con_hess, "constraint Hessian")
  chk(tc(spp_rosen, tcx0, rg, constraints = list(nlc(circ, -Inf, 1.5, circ_jac, hess = "2-point"))),
      ref$tc_con_hess_fd, "2-point constraint Hessian")
  chk(tc(spp_rosen, tcx0, rg, options = list(callback = function(intermediate_result) intermediate_result$nit == 3L)),
      ref$tc_stop, "callback")
  near(tc(spp_rosen, tcx0, rg, bounds = tcbox), ref$tc_bounds, "bounds", 1e-12)
  near(tc(spp_rosen, tcx0, rg, constraints = list(nlc(circ, -Inf, 1.5, circ_jac))), ref$tc_ineq, "inequality", 1e-12)
  near(tc(spp_rosen, tcx0, rg, constraints = list(nlc(circ, -Inf, 1.5, circ_jac)),
          options = list(factorization_method = "SVDFactorization")), ref$tc_svd, "SVD", 1e-12)
  near(tc(spp_rosen, tcx0, rg, constraints = list(nlc(circ, 0.5, Inf, circ_jac))), ref$tc_greater, "greater", 1e-8)
  near(tc(spp_rosen, tcx0, rg, constraints = list(nlc(two, c(0.5, -Inf), c(1.5, 0.3), two_jac))), ref$tc_mixed2,
       "interval and upper", 1e-10)
  near(tc(spp_rosen, tcx0, rg, constraints = list(list(type = "linear", A = matrix(c(1, 1), 1), lb = -Inf, ub = 1))),
       ref$tc_linear, "linear", 1e-12)
  near(tc(f3, c(0.5, 0.5, 0.5), f3_grad, bounds = list(lb = rep(-1, 3), ub = rep(2, 3)),
          constraints = list(nlc(g3, c(1, -Inf), c(1, 0.5), g3_jac))), ref$tc_three, "three variables", 1e-12)
  # the report of verbose = 3 without its last line (the run time)
  report <- function(...) {
    out <- utils::capture.output(r <- tc(spp_rosen, tcx0, rg, options = list(verbose = 3L, maxiter = 3L), ...))
    out[-length(out)]
  }
  sqp3 <- c("| niter |f evals|CG iter|  obj func   |tr radius |   opt    |  c viol  | penalty  |CG stop|",
            "|-------|-------|-------|-------------|----------|----------|----------|----------|-------|",
            "|   1   |   1   |   0   | +2.4200e+01 | 1.00e+00 | 2.16e+02 | 0.00e+00 | 1.00e+00 |   0   |",
            "|   2   |   2   |   1   | +2.4200e+01 | 5.00e-01 | 2.16e+02 | 0.00e+00 | 1.00e+00 |   2   |",
            "|   3   |   3   |   2   | +2.4200e+01 | 2.50e-01 | 2.16e+02 | 0.00e+00 | 1.00e+00 |   2   |",
            "",
            "The maximum number of function evaluations is exceeded.")
  ip3 <- c("| niter |f evals|CG iter|  obj func   |tr radius |   opt    |  c viol  | penalty  |barrier param|CG stop|",
           "|-------|-------|-------|-------------|----------|----------|----------|----------|-------------|-------|",
           "|   1   |   1   |   0   | +2.4200e+01 | 1.00e+00 | 1.51e+02 | 9.40e-01 | 1.00e+00 |  1.00e-01   |   0   |",
           "|   2   |   2   |   1   | +2.4200e+01 | 5.00e-01 | 1.51e+02 | 9.40e-01 | 1.00e+00 |  1.00e-01   |   2   |",
           "|   3   |   3   |   2   | +2.5258e+01 | 5.00e-01 | 1.32e+02 | 0.00e+00 | 1.00e+00 |  1.00e-01   |   2   |",
           "",
           "The maximum number of function evaluations is exceeded.")
  expect_identical(report(), sqp3)
  expect_identical(report(constraints = list(nlc(circ, -Inf, 1.5, circ_jac))), ip3)
  # the errors of SciPy: the ftol of the wrapper's abstol is not an option of
  # trust-constr (a TypeError in SciPy)
  ftol <- "_minimize_trustregion_constr() got an unexpected keyword argument 'ftol'"
  expect_error(tc(spp_rosen, tcx0, rg, options = list(ftol = 1e-8)), ftol, fixed = TRUE)
  expect_error(tc(spp_rosen, tcx0, rg, constraints = list(nlc(circ, 1.5, 1.5, circ_jac)),
                  options = list(factorization_method = "LU")), "Method not allowed for dense array.", fixed = TRUE)
  expect_error(tc(spp_rosen, tcx0, NULL, "2-point"),
               paste("Whenever the gradient is estimated via finite-differences, we require the Hessian to be",
                     "estimated using one of the quasi-Newton strategies."), fixed = TRUE)
  # through solve(): the wrapper's options and the BFGS update without a Hessian
  pr <- OptimizationProblem(function(u, p) spp_rosen(u), tcx0, numeric(), grad = function(u, p) rg(u))
  s <- solve(pr, ScipyTrustConstr())
  expect_identical(s$retcode, "Success")
  chk(c(s$original, nwarn = 0L), ref$tc_none, "solve")
  pc <- OptimizationProblem(function(u, p) spp_rosen(u), tcx0, numeric(), grad = function(u, p) rg(u),
                            cons = function(u, p) circ(u), cons_j = function(u, p) circ_jac(u), lcons = -Inf,
                            ucons = 1.5)
  near(solve(pc, ScipyTrustConstr())$original, ref$tc_ineq, "solve with a constraint", 1e-12)
  expect_error(solve(pr, ScipyTrustConstr(), abstol = 1e-8), ftol, fixed = TRUE)
})

test_that("COBYQA follows the cobyqa package vendored by SciPy", {
  # python3 tools/scipy-probe.py cobyqa: minimize(method = "COBYQA") with the
  # defaults of _minimize_cobyqa. The interpolation system is solved by an
  # eigendecomposition (Accelerate in SciPy, R's LAPACK here), so the points
  # agree to rounding; in the "same" runs the iteration and evaluation
  # counts agree too, in the "path" runs a last bit difference moves a
  # rounding sensitive test (a geometry step, a step on the boundary) and
  # only the solution is compared.
  ref <- list()
  ref$cq_quad <- list(x = c(1.4285713381318552, -0.8571425770641983), fun = -0.7857142857141459, maxcv = 0.0,
    nfev = 31L, nit = 43L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_quad$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_rosen <- list(x = c(1.0000000130122655, 1.0000000295072946), fun = 1.3822831508130037e-15, maxcv = 0.0,
    nfev = 129L, nit = 119L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_rosen$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_poly3 <- list(x = c(0.6005315218149659, -0.5000004757775581, 0.8662979328342875), fun = -1.8095668506407185,
    maxcv = 0.0, nfev = 56L, nit = 58L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_poly3$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_bounds <- list(x = c(0.8, 0.639999999999971), fun = 0.03999999999999998, maxcv = 0.0, nfev = 118L,
    nit = 118L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_bounds$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_scale <- list(x = c(1.428571434867686, -0.8571429450586968), fun = -0.7857142857142708, maxcv = 0.0,
    nfev = 31L, nit = 42L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_scale$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_fixed <- list(x = c(0.6005314023337625, 0.3, 0.8662979985901115), fun = -1.1695668506410366, maxcv = 0.0,
    nfev = 38L, nit = 54L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_fixed$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_linear <- list(x = c(1.1499999999996677, -0.9499999999996677), fun = -0.665, maxcv = 5.551115123125783e-17,
    nfev = 27L, nit = 35L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_linear$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_lineq <- list(x = c(0.6100994253145352, -0.48153121320867026, 0.8714317878941352),
    fun = -1.8089555333374205, maxcv = 2.220446049250313e-16, nfev = 47L, nit = 52L, status = 0L, success = TRUE,
    nwarn = 0L)
  ref$cq_lineq$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_ineq <- list(x = c(0.9072338717684225, 0.8227555546820372), fun = 0.00861565063402164,
    maxcv = 8.440692589317678e-10, nfev = 122L, nit = 114L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_ineq$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_greater <- list(x = c(1.4285713579339068, -0.8571428164665101), fun = -0.7857142857142803, maxcv = 0.0,
    nfev = 31L, nit = 41L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_greater$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_eq <- list(x = c(0.8303742521805225, -0.5572060672202017), fun = -0.42737160880150354,
    maxcv = 3.136602089170992e-11, nfev = 30L, nit = 41L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_eq$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_mixed <- list(x = c(-0.6030200344023003, 0.36927880807571317), fun = 2.572860562779956,
    maxcv = 1.562905360685818e-11, nfev = 92L, nit = 111L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_mixed$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_maxfev <- list(x = c(-0.6693947020320618, 0.414944311764094), fun = 2.896737477651598, maxcv = 0.0,
    nfev = 20L, nit = 12L, status = 5L, success = FALSE, nwarn = 0L)
  ref$cq_maxfev$message <- 'The maximum number of function evaluations has been exceeded'
  ref$cq_maxiter <- list(x = c(-0.8943715453758639, 0.7940942605348666), fun = 3.5920147485205844, maxcv = 0.0,
    nfev = 12L, nit = 5L, status = 6L, success = FALSE, nwarn = 0L)
  ref$cq_maxiter$message <- 'The maximum number of iterations has been exceeded'
  ref$cq_target <- list(x = c(0.5, -0.5), fun = 0.0, maxcv = 0.0, nfev = 5L, nit = 0L, status = 1L, success = TRUE,
    nwarn = 0L)
  ref$cq_target$message <- 'The target objective function value has been reached'
  ref$cq_radius <- list(x = c(1.4285714227539912, -0.857142855638113), fun = -0.7857142857142858, maxcv = 0.0,
    nfev = 43L, nit = 51L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_radius$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_smallbox <- list(x = c(0.4, -0.5999999999999921), fun = 0.13999999999999999, maxcv = 0.0, nfev = 26L,
    nit = 27L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_smallbox$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_initfev <- list(x = c(0.5, 0.5), fun = 2.5, maxcv = 0.0, nfev = 3L, nit = 0L, status = 6L, success = FALSE,
    nwarn = 0L)
  ref$cq_initfev$message <- 'The maximum number of iterations has been exceeded'
  ref$cq_scale_nl <- list(x = c(1.1303748773566087, -0.6995394003463574), fun = -0.6941122106457108, maxcv = 0.0,
    nfev = 153L, nit = 59L, status = 0L, success = TRUE, nwarn = 0L)
  ref$cq_scale_nl$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_fixed_nl <- list(x = c(0.5, 0.3, 0.5), fun = -0.735, maxcv = 0.0, nfev = 56L, nit = 38L, status = 0L,
    success = TRUE, nwarn = 0L)
  ref$cq_fixed_nl$message <- 'The lower bound for the trust-region radius has been reached'
  ref$cq_cb_x <- list(x = c(1.3828125000001013, -0.9238281250000439), fun = -0.7716751098632835, maxcv = 0.0,
    nfev = 9L, nit = 9L, status = 3L, success = TRUE, nwarn = 0L)
  ref$cq_cb_x$message <- 'The callback requested to stop the optimization procedure'
  ref$cq_cb_x$seen <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.7500000000000008, -0.6250000000000051,
    1.2499999999999962, -0.8750000000000275, 1.3828125000001013, -0.9238281250000439, 1.3828125000001013,
    -0.9238281250000439)
  ref$cq_cb_ir <- list(x = c(1.2499999999999962, -0.8750000000000275), fun = -0.7499999999999917, maxcv = 0.0,
    nfev = 7L, nit = 7L, status = 3L, success = TRUE, nwarn = 0L)
  ref$cq_cb_ir$message <- 'The callback requested to stop the optimization procedure'
  ref$cq_report <- c('Starting the optimization procedure.', 'Initial trust-region radius: 1.0.',
    'Final trust-region radius: 1e-06.', 'Maximum number of function evaluations: 1000.',
    'Maximum number of iterations: 12.', '', 'fun([ 5.000e-01  5.000e-01]) = 2.5',
    'fun([ 1.500e+00  5.000e-01]) = 3.0', 'fun([ 5.000e-01  1.500e+00]) = 9.0', 'fun([-5.000e-01  5.000e-01]) = 4.0',
    'fun([ 5.000e-01 -5.000e-01]) = 0.0', '', 'New trust-region radius: 0.1.', 'Number of function evaluations: 5.',
    'Number of iterations: 5.', 'Least value of fun: 0.0.', 'Maximum constraint violation: 0.0.',
    'Corresponding point: [ 5.000e-01 -5.000e-01].')
  ref$cq_report7 <- c('Starting the optimization procedure.', 'Initial trust-region radius: 1.0.',
    'Final trust-region radius: 1e-06.', 'Maximum number of function evaluations: 3500.',
    'Maximum number of iterations: 1.', '', 'fun([-2.000e-05  9.998e-02 ...  5.000e-01  6.000e-01]) = 0.2800000028',
    'fun([ 1.000e+00  9.998e-02 ...  5.000e-01  6.000e-01]) = 0.6799600028000001',
    'fun([-2.000e-05  1.100e+00 ...  5.000e-01  6.000e-01]) = 0.8799600028',
    'fun([-2.000e-05  9.998e-02 ...  5.000e-01  6.000e-01]) = 1.0799600028')
  ref$cq_err_radius <- 'The initial trust-region radius must be positive.'
  ref$cq_warn_radius <- c()
  ref$cq_err_final <- 'The final trust-region radius must be nonnegative.'
  ref$cq_warn_final <- c()
  ref$cq_err_order <- 'The initial trust-region radius must be greater than or equal to the final trust-region radius.'
  ref$cq_warn_order <- c()
  ref$cq_err_maxfev <- 'The maximum number of function evaluations must be positive.'
  ref$cq_warn_maxfev <- c()
  ref$cq_err_maxiter <- 'The maximum number of iterations must be positive.'
  ref$cq_warn_maxiter <- c()
  ref$cq_err_allfixed <- '`lb` and `ub` must be non-empty.'
  ref$cq_warn_allfixed <- c()
  ref$cq_err_unknown <- NULL
  ref$cq_warn_unknown <- c('Unknown solver options: foo, bar')
  quad <- function(x) (x[1] - 1) * (x[1] - 1) + 2 * (x[2] + 0.5) * (x[2] + 0.5) + x[1] * x[2]
  rosen <- function(x) (1 - x[1]) * (1 - x[1]) + 100 * (x[2] - x[1] * x[1]) * (x[2] - x[1] * x[1])
  poly3 <- function(x) x[1] * x[1] * x[1] * x[1] + x[2] * x[2] + x[3] * x[3] * x[3] * x[3] - x[1] * x[3] + x[2] -
    2 * x[3]
  circ <- function(x) x[1] * x[1] + x[2] * x[2]
  two <- function(x) c(x[1] * x[1] + x[2] * x[2], x[1] - x[2])
  nlc <- .sp_nonlinear_constraint
  lc <- function(A, lb, ub) list(type = "linear", A = matrix(A, nrow = 1L), lb = lb, ub = ub)
  bx <- function(lb, ub) list(lb = lb, ub = ub)
  cq <- function(fun, x0, bounds = NULL, constraints = list(), ...)
    .sp_minimize(fun, x0, method = "COBYQA", bounds = bounds, constraints = constraints, options = list(...))
  x2 <- c(0.5, 0.5); x3 <- c(0.5, 0.5, 0.5); xr <- c(-1.2, 1)
  runs <- list(
    cq_quad = list(cq(quad, x2), "same", 1e-10),
    cq_rosen = list(cq(rosen, xr), "path", 1e-6),
    cq_poly3 = list(cq(poly3, x3), "same", 1e-10),
    cq_bounds = list(cq(rosen, xr, bx(c(-2, -0.5), c(0.8, 2))), "path", 1e-6),
    cq_scale = list(cq(quad, x2, bx(c(-2, -3), c(2, 1)), scale = TRUE), "same", 1e-12),
    cq_fixed = list(cq(poly3, x3, bx(c(-2, 0.3, -2), c(2, 0.3, 2))), "same", 1e-10),
    cq_linear = list(cq(quad, x2, constraints = list(lc(c(1, 1), -Inf, 0.2))), "path", 1e-6),
    cq_lineq = list(cq(poly3, x3, constraints = list(lc(c(1, 1, 1), 1, 1))), "path", 1e-6),
    cq_ineq = list(cq(rosen, xr, constraints = list(nlc(circ, -Inf, 1.5))), "path", 1e-6),
    cq_greater = list(cq(quad, x2, constraints = list(nlc(circ, 2, Inf))), "same", 1e-7),
    cq_eq = list(cq(quad, x2, constraints = list(nlc(circ, 1, 1))), "same", 1e-7),
    cq_mixed = list(cq(rosen, xr, bx(c(-2, -2), c(2, 0.9)), list(nlc(two, c(0.5, -Inf), c(1.5, 0.3)))), "same",
                    1e-5),
    cq_maxfev = list(cq(rosen, xr, maxfev = 20L), "same", 1e-11),
    cq_maxiter = list(cq(rosen, xr, maxiter = 5L), "same", 1e-11),
    cq_target = list(cq(quad, x2, f_target = 0.5), "same", 0),
    cq_radius = list(cq(quad, x2, initial_tr_radius = 0.1, final_tr_radius = 1e-8), "path", 1e-7),
    cq_smallbox = list(cq(quad, x2, bx(c(0, -0.6), c(0.4, 0.6))), "same", 1e-12),
    cq_initfev = list(cq(quad, x2, maxfev = 3L), "same", 0),
    # SciPy's wrong runs, reproduced by .upstream_bugs = TRUE
    cq_scale_nl = list(cq(quad, x2, bx(c(-2, -3), c(2, 1)), list(nlc(circ, -Inf, 1)), scale = TRUE,
                          .upstream_bugs = TRUE), "same", 1e-12),
    cq_fixed_nl = list(cq(poly3, x3, bx(c(-2, 0.3, -2), c(2, 0.3, 2)), list(nlc(circ, -Inf, 0.5)),
                          .upstream_bugs = TRUE), "same", 1e-12))
  for (nm in names(runs)) {
    r <- runs[[nm]][[1]]
    e <- ref[[nm]]
    expect_equal(r$x, e$x, tolerance = 1e-4, info = nm)
    expect_equal(r$fun, e$fun, tolerance = 1e-6, info = nm)
    expect_identical(r$status, e$status, info = nm)
    expect_identical(r$success, e$success, info = nm)
    expect_identical(r$message, e$message, info = nm)
    if (ref_platform()) {
      expect_lte(max(abs(r$x - e$x)), runs[[nm]][[3]])
      expect_lte(abs(r$fun - e$fun), 1e-9)
      expect_lte(abs(r$maxcv - e$maxcv), 1e-8)
      if (runs[[nm]][[2]] == "same") {
        expect_identical(r$nfev, e$nfev, info = nm)
        expect_identical(r$nit, e$nit, info = nm)
      }
    }
  }
  # the defects of SciPy fixed: the constraint holds at the point of the user
  r <- cq(quad, x2, bx(c(-2, -3), c(2, 1)), list(nlc(circ, -Inf, 1)), scale = TRUE)
  expect_equal(r$x, ref$cq_eq$x, tolerance = 1e-6)
  expect_lte(circ(r$x), 1 + 1e-8)
  r <- cq(poly3, x3, bx(c(-2, 0.3, -2), c(2, 0.3, 2)), list(nlc(circ, -Inf, 0.5)))
  expect_equal(r$x, ref$cq_fixed$x, tolerance = 1e-6)
  # callbacks: the best point after every evaluation; StopIteration stops
  seen <- numeric()
  cb_x <- function(x) {
    seen <<- c(seen, x)
    if (length(seen) == 18L) stop(structure(class = c("StopIteration", "error", "condition"),
                                            list(message = "", call = NULL)))
  }
  r <- cq(quad, x2, callback = cb_x)
  expect_identical(r[c("status", "message")], ref$cq_cb_x[c("status", "message")])
  expect_equal(seen, ref$cq_cb_x$seen, tolerance = 1e-8)
  if (ref_platform()) {
    expect_lte(max(abs(seen - ref$cq_cb_x$seen)), 1e-12)
    expect_identical(r[c("nfev", "nit")], ref$cq_cb_x[c("nfev", "nit")])
  }
  cb_ir <- function(intermediate_result) {
    if (intermediate_result$fun < -0.7) stop(structure(class = c("StopIteration", "error", "condition"),
                                                      list(message = "", call = NULL)))
  }
  r <- cq(quad, x2, callback = cb_ir)
  expect_identical(r$status, 3L)
  expect_lt(r$fun, -0.7)
  if (ref_platform()) expect_equal(r[c("x", "nfev", "nit")], ref$cq_cb_ir[c("x", "nfev", "nit")], tolerance = 1e-10)
  # the report of disp = TRUE (its numbers print with all digits)
  sq7 <- function(x) {
    d <- x - 0.3
    s <- d[1] * d[1]
    for (i in 2:7) s <- s + d[i] * d[i]
    s
  }
  out <- capture.output(cq(quad, x2, constraints = list(nlc(circ, -Inf, 1)), disp = TRUE, maxiter = 12L))
  expect_identical(out[1:6], ref$cq_report[1:6])
  out7 <- capture.output(cq(sq7, (0:6) / 10 - 2e-5, disp = TRUE, maxiter = 1L))
  expect_identical(out7[1:6], ref$cq_report7[1:6])
  if (ref_platform()) {
    expect_identical(out[1:18], ref$cq_report)
    expect_identical(out7[1:10], ref$cq_report7)
  }
  # the errors and warnings of SciPy
  expect_error(cq(quad, x2, initial_tr_radius = 0), ref$cq_err_radius, fixed = TRUE)
  expect_error(cq(quad, x2, final_tr_radius = -1), ref$cq_err_final, fixed = TRUE)
  expect_error(cq(quad, x2, initial_tr_radius = 1e-3, final_tr_radius = 1e-2), ref$cq_err_order, fixed = TRUE)
  expect_error(cq(quad, x2, maxfev = 0L), ref$cq_err_maxfev, fixed = TRUE)
  expect_error(cq(quad, x2, maxiter = 0L), ref$cq_err_maxiter, fixed = TRUE)
  expect_error(cq(quad, x2, bx(c(1, 2), c(1, 2))), ref$cq_err_allfixed, fixed = TRUE)
  expect_warning(cq(quad, x2, maxiter = 2L, foo = 1, bar = 2), ref$cq_warn_unknown, fixed = TRUE)
  # through solve(): the abstol of the wrapper is the feasibility tolerance
  pr <- OptimizationProblem(function(u, p) quad(u), x2, numeric())
  s <- solve(pr, ScipyCOBYQA())
  expect_identical(s$retcode, "Success")
  expect_equal(s$u, ref$cq_quad$x, tolerance = 1e-8)
  pc <- OptimizationProblem(function(u, p) quad(u), x2, numeric(), cons = function(u, p) circ(u), lcons = 1,
                            ucons = 1)
  s <- solve(pc, ScipyCOBYQA(), abstol = 1e-8)
  expect_equal(s$u, ref$cq_eq$x, tolerance = 1e-6)
})
