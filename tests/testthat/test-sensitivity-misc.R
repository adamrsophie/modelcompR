# Reference values from tools/scimlsensitivity-probe.jl (SciMLSensitivity
# 7.119.14).

test_that("supports_functor_params, SundialsAdjoint, SensitivityAlg and ReactantDualTag", {
  expect_true(supports_functor_params(GaussAdjoint()))
  expect_true(supports_functor_params(GaussKronrodAdjoint()))
  for (sa in list(QuadratureAdjoint(), InterpolatingAdjoint(), BacksolveAdjoint(), ForwardSensitivity(),
                  ForwardDiffSensitivity(), SteadyStateAdjoint()))
    expect_false(supports_functor_params(sa))
  expect_error(supports_functor_params(1), "sensitivity algorithm")
  sa <- SundialsAdjoint()
  expect_identical(sa$steps, 150)
  expect_identical(sa$interp, "hermite")
  expect_true(sa$quad_error_control)
  expect_error(SundialsAdjoint(interp = "linear"),
               "SundialsAdjoint `interp` must be `:hermite` or `:polynomial`, got `linear`.", fixed = TRUE)
  expect_message(SensitivityAlg(), "The SensitivityAlg choice mechanism was completely overhauled")
  expect_s3_class(ReactantDualTag(), "ReactantDualTag")
  # SundialsAdjoint runs the checkpointed interpolating adjoint
  prob <- ODEProblem(function(u, p, t) c(-p[1] * u[1] + u[2], -p[2] * u[2]), c(1, 0.5), c(0, 1), c(1.3, 0.4))
  sol <- solve(prob, Tsit5(), abstol = 1e-10, reltol = 1e-10)
  dg <- function(u, p, t, i) c(1, 2)
  a <- adjoint_sensitivities(sol, Tsit5(), t = 1, dgdu_discrete = dg, sensealg = SundialsAdjoint())
  b <- adjoint_sensitivities(sol, Tsit5(), t = 1, dgdu_discrete = dg, sensealg = InterpolatingAdjoint(checkpointing = TRUE))
  expect_identical(a, b)
})

lv_sens <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1] * u[2], -p[3] * u[2] + p[4] * u[1] * u[2] + t * p[1])
p_sens <- c(1.5, 1.0, 3.0, 1.0)

test_that("ODEForwardSensitivityFunction is the forward sensitivity right-hand side", {
  z <- c(1.2, 0.7, 0.1, -0.2, 0.3, 0.4, -0.5, 0.6, 0.7, -0.8)
  jl <- c(0.9599999999999999, -0.8099999999999998, 1.52, 0.7300000000000001, -1.08, -0.5100000000000002, -1.12,
          -2.13, 1.5199999999999998, 2.7700000000000005)
  fs <- ODEForwardSensitivityFunction(ODEProblem(lv_sens, c(1, 1), c(0, 2), p_sens))
  expect_s3_class(fs, "ODEForwardSensitivityFunction")
  expect_false(inherits(fs, SensitivityFunction))
  expect_identical(attr(fs, "numindvar"), 2L)
  expect_identical(attr(fs, "numparams"), 4L)
  expect_equal(fs(z, p_sens, 0.3), jl, tolerance = 1e-9)
  # with the symbolic Jacobians of a System (unknowns y, x; parameters a, b, d, c)
  tt <- t_nounits; D <- D_nounits
  variables("x_fs(t)", "y_fs(t)"); parameters("a_fs", "b_fs", "c_fs", "d_fs")
  sys <- mtkcompile(System(list(D(x_fs) ~ a_fs * x_fs - b_fs * x_fs * y_fs,
                                D(y_fs) ~ -c_fs * y_fs + d_fs * x_fs * y_fs + tt * a_fs), tt, name = "lvfs"))
  sprob <- ODEProblem(sys, list(x_fs ~ 1, y_fs ~ 1, a_fs ~ 1.5, b_fs ~ 1, c_fs ~ 3, d_fs ~ 1), c(0, 2))
  ux <- vapply(unknowns(sys), function(v) .key(v) == .key(x_fs), TRUE)
  iu <- if (ux[1L]) c(1, 2) else c(2, 1)
  pk <- vapply(.problem_params(sys), .key, "")
  ip <- match(pk, c(.key(a_fs), .key(b_fs), .key(c_fs), .key(d_fs)))
  zS <- matrix(z[3:10], 2, 4)
  jS <- matrix(jl[3:10], 2, 4)
  out <- ODEForwardSensitivityFunction(sprob)(c(z[iu], zS[iu, ip]), sprob$p, 0.3)
  expect_equal(out, c(jl[iu], jS[iu, ip]), tolerance = 1e-14)
})

test_that("SensitivityFunction tags the adjoint right-hand sides", {
  sol <- solve(ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2), Tsit5())
  ap <- ODEAdjointProblem(sol, InterpolatingAdjoint(), t = 1, dgdu_discrete = function(u, p, t, i) 1)
  expect_true(inherits(ap$f$f, SensitivityFunction))
  expect_s3_class(ap$f$f, "ODEInterpolatingAdjointSensitivityFunction")
  ap <- ODEAdjointProblem(sol, QuadratureAdjoint(), t = 1, dgdu_discrete = function(u, p, t, i) 1)
  expect_s3_class(ap$f$f, "ODEQuadratureAdjointSensitivityFunction")
})

test_that("AdjointSensitivityIntegrand is lambda' f_p + dgdp along the adjoint solution", {
  sol <- solve(ODEProblem(lv_sens, c(1, 1), c(0, 2), p_sens), Tsit5(), abstol = 1e-12, reltol = 1e-12)
  ap <- ODEAdjointProblem(sol, QuadratureAdjoint(), Tsit5(), t = 2, dgdu_discrete = function(u, p, t, i) c(1, 1),
                          dgdu_continuous = function(u, p, t) c(u[1], 0),
                          dgdp_continuous = function(u, p, t) numeric(4))
  adj <- solve(ap, Tsit5(), abstol = 1e-12, reltol = 1e-12)
  ig <- AdjointSensitivityIntegrand(sol, adj)
  expect_s3_class(ig, "AdjointSensitivityIntegrand")
  expect_equal(ig(0), c(4.388385806326229, -4.388385806326229, 3.209238117606435, -3.209238117606435),
               tolerance = 1e-9)
  expect_equal(ig(0.5), c(3.1238049083504036, -2.858171788923313, 2.259427363736097, -3.373393773556924),
               tolerance = 1e-9)
  expect_equal(ig(1.3), c(2.9815255136694954, -6.936503197278934, 2.6574441960460775, -7.2409405730025735),
               tolerance = 1e-9)
  expect_equal(.sens_quadgk(ig, 0, 2, 1e-8, 1e-8),
               c(6.758264897027363, -10.745074463633323, 3.7275618030609388, -7.191688506408433), tolerance = 1e-9)
  ig2 <- AdjointSensitivityIntegrand(sol, adj, dgdp = function(u, p, t) c(u[1] * u[2], 0, 0, t))
  expect_equal(ig2(0.5), c(3.949925808933976, -2.858171788923313, 2.259427363736097, -2.873393773556924),
               tolerance = 1e-9)
})

test_that("StochasticTransformedFunction subtracts the Ito correction", {
  fs <- function(u, p, t) c(p[1] * u[1] - u[2], p[2] * u[2] + u[1] * t)
  gs <- function(u, p, t) c(p[3] * u[1] * u[2], p[4] * sin(u[2]))
  pp <- c(1.01, -0.5, 0.3, 0.7)
  sp <- SDEProblem(fs, gs, c(1, 0.5), c(0, 1), pp)
  ssol <- solve(sp, EM(), dt = 0.1, seed = 1)
  tf <- StochasticTransformedFunction(ssol)
  expect_s3_class(tf, "StochasticTransformedFunction")
  expect_equal(tf(c(1.3, -0.4), pp, 0.25), c(1.69428, 0.7615922422703831), tolerance = 1e-9)
  gm <- function(u, p, t) matrix(c(p[3] * u[1], p[4] * u[1] * u[2], u[2], cos(u[1])), 2, 2)
  mp <- SDEProblem(fs, gm, c(1, 0.5), c(0, 1), pp, noise_rate_prototype = matrix(0, 2, 2))
  tm <- StochasticTransformedFunction(solve(mp, EM(), dt = 0.1, seed = 1))
  expect_equal(tm(c(1.3, -0.4), pp, 0.25), c(1.7518306859107322, 1.25624), tolerance = 1e-9)
  tc <- StochasticTransformedFunction(ssol, corfunc_analytical = function(u, p, t) c(u[1], 2 * u[2]))
  expect_equal(tc(c(1.3, -0.4), pp, 0.25), c(0.41300000000000003, 1.3250000000000002), tolerance = 1e-15)
})

test_that("UnconstrainedOptimizationAdjoint and OptimizationAdjoint differentiate the optimum", {
  fu <- function(x, p) (x[1] - p[1])^2 + p[2] * (x[2] - x[1]^2)^2 + p[3] * x[1]
  gu <- function(x, p) c(2 * (x[1] - p[1]) - 4 * p[2] * x[1] * (x[2] - x[1]^2) + p[3], 2 * p[2] * (x[2] - x[1]^2))
  hu <- function(x, p) rbind(c(2 - 4 * p[2] * (x[2] - 3 * x[1]^2), -4 * p[2] * x[1]), c(-4 * p[2] * x[1], 2 * p[2]))
  pu <- c(1, 3, 0.5)
  probu <- OptimizationProblem(OptimizationFunction(fu, grad = gu, hess = hu), c(0, 0), pu)
  jl <- c(3.999999999999982, 2.960594732333756e-16, -1.999999999999991)
  # at the optimum of Julia's BFGS
  su <- .opt_solution(probu, NULL, c(0.7499999999999978, 0.5624999999999962), 0, "Success", NULL)
  dgdu <- function(u, p) c(1, 2)
  for (sa in list(UnconstrainedOptimizationAdjoint(), OptimizationAdjoint()))
    expect_equal(adjoint_sensitivities(su, sensealg = sa, dgdu = dgdu)$dp, jl, tolerance = 1e-9)
  expect_error(adjoint_sensitivities(su, sensealg = OptimizationAdjoint()), "dgdu must be specified")
  # through solve_gradient; the default is the unconstrained adjoint
  expect_s3_class(.sens_default_alg(probu), "UnconstrainedOptimizationAdjoint")
  # (to the accuracy of the optimum BFGS stops at)
  g <- solve_gradient(function(sol) sum(sol$u * c(1, 2)), probu, BFGS(), g_abstol = 1e-14)
  expect_equal(g$p, jl, tolerance = 1e-6)
  # an equality, an active and an inactive inequality and an active bound
  fc <- function(x, p) (x[1] - p[1])^2 + (x[2] - p[2])^2 + p[4] * (x[3] - p[3])^2 + p[6] * x[1] * x[2]
  cons <- function(x, p) c(x[1] + 2 * x[2] + x[3] - p[5], x[1]^2 + x[2]^2 - 1, x[3] - x[1] - 10)
  pc <- c(2, 1.5, -1, 1, 1.5, 0.3)
  probc <- OptimizationProblem(fc, c(0.1, 0.1, 0.5), pc, cons = cons, lcons = c(0, -Inf, -Inf), ucons = c(0, 0, 0),
                               lb = c(-Inf, -Inf, 0), ub = c(Inf, Inf, Inf))
  expect_s3_class(.sens_default_alg(probc), "OptimizationAdjoint")
  xj <- c(0.9633249580710704, 0.2683375209644633, 3.0837022369139185e-15)
  sc <- .opt_solution(probc, NULL, xj, fc(xj, pc), "Success", NULL)
  expect_equal(adjoint_sensitivities(sc, dgdu = function(u, p) c(1, 2, -1))$dp, c(0, 0, 0, 0, 1, 0),
               tolerance = 1e-9)
  # x1 + 2 x2 = p5 and x1^2 + x2^2 = 1 are active: dx1/dp5 = x2 / (x2 - 2 x1)
  expect_equal(adjoint_sensitivities(sc, dgdu = function(u, p) c(1, 0, 0))$dp,
               c(0, 0, 0, 0, xj[2] / (xj[2] - 2 * xj[1]), 0), tolerance = 1e-9)
})

test_that("the sensitivity fixes of the final review hold", {
  # a mass matrix whose differential equation is not in the row of its
  # variable: lambda is indexed by the equations
  f <- function(u, p, t) c(u[2] - p[1] * u[1], 0)
  prob <- ODEProblem(ODEFunction(f, mass_matrix = matrix(c(0, 1, 0, 0), 2)), c(1, 2), c(0, 1), 2)
  sol <- solve(prob, Rodas5P(), abstol = 1e-10, reltol = 1e-10, dense = TRUE)
  r <- adjoint_sensitivities(sol, Rodas5P(), t = 1, dgdu_discrete = function(u, p, t, i) c(0, 1),
                             sensealg = InterpolatingAdjoint(), abstol = 1e-10, reltol = 1e-10)
  expect_equal(c(r$du0, r$dp), c(2, 0, 1), tolerance = 1e-6)
  # an affect applied at the initial time by the initialization
  mk <- function(cb) ODEProblem(function(u, p, t) 0 * u, 1, c(0, 1), 2, callback = cb)
  last <- function(sol) tail(as.matrix(sol)[, 1], 1)
  for (cb in list(PresetTimeCallback(0, function(i) i$u[1] <- i$p[1] * i$u[1]),
                  PeriodicCallback(function(i) i$u[1] <- i$p[1] * i$u[1], 5, initial_affect = TRUE)))
    for (sa in list(ForwardSensitivity(), ForwardDiffSensitivity(), InterpolatingAdjoint(), QuadratureAdjoint())) {
      g <- solve_gradient(last, mk(cb), Tsit5(), sensealg = sa)
      expect_equal(unname(c(g$u0, g$p)), c(2, 1), tolerance = 1e-8, label = sa$name)
    }
  # saveat stored in the problem marks fixed sample times
  cb <- ContinuousCallback(function(u, t, integ) u[1] - integ$p[1], function(integ) NULL)
  pp <- ODEProblem(function(u, p, t) 1, 0, c(0, 1), 0.5, callback = cb, saveat = c(0, 0.5, 1))
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint())) {
    g <- solve_gradient(function(sol) as.matrix(sol)[2, 1], pp, Tsit5(), sensealg = sa)
    expect_equal(unname(c(g$u0, g$p)), c(1, 0), tolerance = 1e-8)
  }
  # a running cost given with dgdu_continuous still differentiates g in p
  s0 <- solve(ODEProblem(function(u, p, t) 0 * u, 1, c(0, 1), 3), Tsit5(), dense = TRUE)
  r <- adjoint_sensitivities(s0, Tsit5(), g = function(u, p, t) p[1], dgdu_continuous = function(u, p, t) 0,
                             sensealg = InterpolatingAdjoint())
  expect_equal(r$dp, 1, tolerance = 1e-8)
})

test_that("losses of the event times keep the parameters of each saved row", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("mvr_x", t); b <- Parameter("mvr_b"); a <- Parameter("mvr_a")
  sys <- mtkcompile(System(list(D(x) ~ a), t, name = "mvr",
    continuous_events = list(SymbolicContinuousCallback(list(x ~ 0.47), list(b ~ Pre(b) + 1),
                                                        discrete_parameters = list(b)))))
  prob <- ODEProblem(sys, list(x ~ 0, b ~ 1, a ~ 1), c(0, 1))
  loss <- function(sol) {
    i <- which(duplicated(sol$t))
    sum(sol[b * x][c(i - 1, i)]) + sum(sol[x][i]^2) + tail(sol[x], 1)^2
  }
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint(), BacksolveAdjoint())) {
    g <- solve_gradient(loss, prob, RK4(), dt = 0.05, adaptive = FALSE, sensealg = sa)
    expect_equal(unname(c(g$u0, g$p)), c(2, 2, 0.94), tolerance = 1e-6, label = sa$name)
  }
})

test_that("DAEProblem gradients keep the callbacks and the interpolation in the DAE variables", {
  fr <- function(du, u, p, t) c(du[1] - p[1], u[2] - u[1])
  cb <- ContinuousCallback(function(u, t, integ) sum(u) - 1, function(integ) terminate(integ))
  prob <- DAEProblem(fr, du0 = c(0.25, 0), u0 = c(0, 0), tspan = c(0, 5), p = 0.25,
                     differential_vars = c(TRUE, FALSE), callback = cb)
  loss <- function(sol) tail(sol$t, 1) + tail(as.matrix(sol)[, 1], 1)
  g <- solve_gradient(loss, prob, DFBDF(), reltol = 1e-8, abstol = 1e-10, sensealg = InterpolatingAdjoint())
  expect_equal(c(unname(g$p), g$value), c(-8, 2.5), tolerance = 1e-6)
  fr2 <- function(du, u, p, t) c(du[1] - p[1] * u[1], u[2] - u[1])
  prob2 <- DAEProblem(fr2, du0 = c(1, 0), u0 = c(1, 1), tspan = c(0, 1), p = 1, differential_vars = c(TRUE, FALSE))
  g <- solve_gradient(function(sol) sol(0.5)[1], prob2, DFBDF(), saveat = c(0, 1), dense = TRUE,
                      reltol = 1e-8, abstol = 1e-10, sensealg = ForwardDiffSensitivity())
  expect_equal(c(unname(g$p), unname(g$value)), c(0.5, 1) * exp(0.5), tolerance = 1e-6)
  s <- solve_gradient(function(sol) sol(0.5)[1], prob2, DFBDF(), reltol = 1e-8, abstol = 1e-10)$sol
  expect_equal(unname(s(0.25)), rep(exp(0.25), 2), tolerance = 1e-4)
})

test_that("events saving one side, integrator side effects and late parameter changes", {
  f <- function(u, p, t) c(-p[1] * u[1] + u[2], -p[2] * u[2])
  u0 <- c(1, 0.5); p <- c(1.2, 0.8, 0.6)
  aff <- function(integ) integ$u[1] <- integ$u[1] * integ$p[3] + integ$u[2]
  fdg <- function(mk, loss, dots, alg) {
    F <- function(th) loss(do.call(solve, c(list(mk(th[1:2], th[3:5]), alg), dots)))
    th <- c(u0, p)
    vapply(seq_along(th), function(j) {
      a <- th; b <- th; h <- 1e-6 * max(1, abs(th[j]))
      a[j] <- a[j] + h; b[j] <- b[j] - h
      (F(a) - F(b)) / (2 * h)
    }, 0)
  }
  chk <- function(mk, loss, dots, sas = list(ForwardDiffSensitivity(), InterpolatingAdjoint(), GaussAdjoint()),
                  tol = 1e-5, alg = Tsit5()) {
    ref <- fdg(mk, loss, c(dots, list(reltol = 1e-10, abstol = 1e-12)), alg)
    for (sa in sas) {
      g <- do.call(solve_gradient, c(list(loss, mk(u0, p), alg), dots,
                                     list(reltol = 1e-10, abstol = 1e-12, sensealg = sa)))
      expect_equal(unname(c(g$u0, g$p)), ref, tolerance = tol, label = sa$name)
    }
  }
  # a continuous event saving only its left limit, read by the loss
  cb1 <- ContinuousCallback(function(u, t, integ) u[1] - 0.8, function(integ) {
    integ$u[1] <- 0.4 * integ$u[1] + integ$p[3]
    integ$u[2] <- 2 * integ$u[2]
  }, save_positions = c(TRUE, FALSE))
  chk(function(u0, p) ODEProblem(f, u0, c(0, 1), p, callback = cb1),
      function(sol) {
        k <- which(abs(sol$t / 0.05 - round(sol$t / 0.05)) > 1e-9)[1L]
        U <- as.matrix(sol)
        sum(U[k, ]^2) + 3 * U[k, 2] + sol$t[k]
      }, list(dt = 0.05, adaptive = FALSE), alg = RK4())
  # an affect that adds stops (IterativeCallback), also at the initial time
  for (ia in c(FALSE, TRUE)) {
    cb2 <- IterativeCallback(function(integ) if (integ$t + 0.4 < 1) integ$t + 0.4 else NULL, aff, initial_affect = ia)
    chk(function(u0, p) ODEProblem(f, u0, c(0, 1), p, callback = cb2), function(sol) sum(as.matrix(sol)^2),
        list(saveat = 0.25))
  }
  # parameters changed after the last saved time, and the initial ones
  cb3 <- PresetTimeCallback(0.7, function(integ) integ$p[1] <- integ$p[1] + 0.5 * integ$u[1] * integ$p[3])
  mk3 <- function(u0, p) ODEProblem(f, u0, c(0, 1), p, callback = cb3)
  chk(mk3, function(sol) sum(as.matrix(sol)^2) + .sol_env(sol)$p_final[1]^2, list(saveat = c(0, 0.25, 0.5)))
  chk(mk3, function(sol) sum(as.matrix(sol)^2) + sol$prob$p[1]^2, list(saveat = c(0, 0.25, 0.5, 1)))
  # a step saveat that does not divide the span, with an interpolating loss
  pe <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
  g <- solve_gradient(function(sol) sum(as.matrix(sol)) + sol(0.33), pe, Tsit5(), saveat = 0.3,
                      reltol = 1e-10, abstol = 1e-12)
  ts <- c(0, 0.3, 0.6, 0.9, 1)
  expect_equal(unname(g$value), sum(exp(-2 * ts)) + exp(-0.66), tolerance = 1e-7)
  # checkpointing when the start and end times are not saved
  l2 <- function(sol) sum(as.matrix(sol)^2)
  pr <- ODEProblem(f, u0, c(0, 1), p[1:2])
  a <- solve_gradient(l2, pr, Tsit5(), saveat = c(0.2, 0.5, 0.9), reltol = 1e-10, abstol = 1e-12,
                      sensealg = InterpolatingAdjoint(checkpointing = TRUE))
  b <- solve_gradient(l2, pr, Tsit5(), saveat = c(0.2, 0.5, 0.9), reltol = 1e-10, abstol = 1e-12,
                      sensealg = ForwardDiffSensitivity())
  expect_equal(c(a$u0, a$p), c(b$u0, b$p), tolerance = 1e-8)
})

test_that("an affect at the initial time of a DAE reinitializes the algebraic variables", {
  yc <- function(x, a) uniroot(function(y) y + y^3 - a * x, c(-5, 5), tol = 1e-14)$root
  fa <- function(u, p, t) c(-p[1] * u[1] + u[2], u[2] + u[2]^3 - p[2] * u[1])
  cb <- PresetTimeCallback(c(0, 0.6), function(integ) integ$u[1] <- integ$u[1] * integ$p[3])
  pr <- ODEProblem(ODEFunction(fa, mass_matrix = diag(c(1, 0))), c(1, yc(1, 0.5)), c(0, 1.5), c(0.7, 0.5, 1.3),
                   callback = cb)
  s <- solve(pr, Rodas5P(), reltol = 1e-8, abstol = 1e-10, saveat = c(0, 1.5))
  expect_identical(format(s$retcode), "Success")
  U <- as.matrix(s)
  expect_equal(U[1, 2], yc(1.3, 0.5), tolerance = 1e-8)
})
