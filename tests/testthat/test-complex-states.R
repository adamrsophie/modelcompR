# Complex and array states (OrdinaryDiffEq test/InterfaceI/complex_tests.jl,
# StochasticDiffEq test/complex_tests.jl, SciMLBase ode_problems.jl), the
# allowscomplex trait, interpolant derivatives and the keyword checks of SDE
# and DDE solves. Reference values: tools/complex-probe.jl.

cx_H <- function(t) -1i * matrix(c(t, 1, 1, -t), 2)
cx_fun <- function(u, p, t) as.vector(cx_H(t) %*% u)
cx_T <- 0.1

test_that("explicit methods integrate the complex Landau-Zener problem", {
  for (alg in list(Midpoint(), RK4(), DP5(), Tsit5(), Vern7())) {
    prob <- ODEProblem(cx_fun, c(1 + 0i, 0), c(-cx_T, cx_T))
    sol <- solve(prob, alg)
    expect_equal(sol$retcode, "Success", label = alg$name)
    expect_true(is.complex(sol$u[[length(sol$t)]]), label = alg$name)
    expect_equal(sqrt(sum(Mod(sol(cx_T))^2)), 1, tolerance = 1e-2, label = alg$name)
  }
})

test_that("implicit methods with finite differences integrate complex states", {
  for (alg in list(ImplicitEuler(autodiff = AutoFiniteDiff()), Trapezoid(autodiff = AutoFiniteDiff()),
                   Kvaerno3(autodiff = AutoFiniteDiff()), Rosenbrock23(autodiff = AutoFiniteDiff()))) {
    prob <- ODEProblem(cx_fun, c(1 + 0i, 0), c(-cx_T, cx_T))
    sol <- solve(prob, alg)
    expect_equal(sol$retcode, "Success", label = alg$name)
    expect_equal(sqrt(sum(Mod(sol(cx_T))^2)), 1, tolerance = 1e-2, label = alg$name)
  }
})

test_that("a scalar complex ODE keeps complex values, dense output and derivatives", {
  lam <- -1 + 2i
  prob <- ODEProblem(function(u, p, t) lam * u, 1 + 1i, c(0, 1))
  sol <- solve(prob, Tsit5(), abstol = 1e-10, reltol = 1e-10)
  expect_true(is.complex(prob$u0))
  expect_equal(sol(1), (1 + 1i) * exp(lam), tolerance = 1e-8)
  expect_equal(sol(0.37), (1 + 1i) * exp(lam * 0.37), tolerance = 1e-8)
  expect_equal(sol(0.37, deriv = 1), lam * (1 + 1i) * exp(lam * 0.37), tolerance = 1e-7)
  integ <- init(prob, Tsit5())
  step_integrator(integ)
  expect_true(is.complex(integ$u))
  expect_true(is.complex(get_du(integ)))
})

test_that("complex SDE states stay complex with complex noise", {
  sde_cf <- function(u, p, t) 1.01 * u
  for (alg in list(SRIW1(), SRIW2(), SOSRI(), SOSRI2(), SRA2(), SRA3(), SOSRA(), SOSRA2(), SKenCarp(),
                   ImplicitEulerHeun())) {
    sol <- solve(SDEProblem(sde_cf, sde_cf, 1 + 1i, c(0, 1)), alg, seed = 3)
    expect_true(is.complex(sol$u[[length(sol$t)]]), label = alg$name)
    expect_true(is.double(sol$t), label = alg$name)
  }
  sp <- SDEProblem(function(u, p, t) { stopifnot(is.matrix(u)); 1.01 * u },
                   function(u, p, t) { stopifnot(is.matrix(u)); 1.01 * u }, matrix(1 + 1i, 2, 4), c(0, 1))
  for (alg in list(SRIW1(), ImplicitEM(), ISSEM())) {
    sol <- solve(sp, alg, adaptive = FALSE, dt = 0.01, seed = 3)
    expect_true(is.double(sol$t))
    uend <- sol$u[[length(sol$t)]]
    expect_true(is.complex(uend), label = alg$name)
    expect_equal(dim(uend), c(2L, 4L))
  }
})

test_that("FunctionMap and NewtonRaphson take complex states", {
  sol <- solve(DiscreteProblem(function(u, p, t) 0.5i * u, 1 + 0i, c(0, 3)), FunctionMap())
  expect_equal(unlist(sol$u), c(1, 0.5i, -0.25, -0.125i))
  s1 <- solve(NonlinearProblem(function(u, p) u^2 + 1, 1 + 1i), NewtonRaphson())
  expect_equal(s1$retcode, "Success")
  expect_equal(s1$u, 1i, tolerance = 1e-10)
  s2 <- solve(NonlinearProblem(function(u, p) u^2 + 1, c(0.5 - 2i, 1 + 0.5i)), NewtonRaphson())
  expect_equal(s2$u, c(-1i, 1i), tolerance = 1e-10)
})

test_that("allowscomplex reports the upstream traits and solve checks them", {
  expect_true(allowscomplex(Tsit5()))
  expect_true(allowscomplex(Rodas5P()))
  expect_true(allowscomplex(SRIW1()))
  expect_true(allowscomplex(MethodOfSteps(Tsit5())))
  expect_true(allowscomplex(GlobalRichardson(Tsit5())))
  expect_true(allowscomplex(FunctionMap()))
  expect_true(allowscomplex(MIRK4()))
  expect_false(allowscomplex(GlobalErrorEstimation(Tsit5())))
  expect_false(allowscomplex(NewtonRaphson()))
  expect_error(solve(ODEProblem(function(u, p, t) -u, 1 + 1i, c(0, 1)), GlobalErrorEstimation(Tsit5())),
               class = "ComplexSupportError")
  expect_error(solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1 + 1i)), Tsit5()), class = "ComplexTspanError")
})

test_that("a matrix state is passed to f with its shape and restored in the solution", {
  cx_A <- matrix(c(-1, 0, 0.5, -2), 2)
  u0 <- matrix(c(1, 3, 2, 4), 2)
  prob <- ODEProblem(function(u, p, t) { stopifnot(is.matrix(u)); cx_A %*% u }, u0, c(0, 1))
  expect_equal(dim(prob$u0), c(2L, 2L))
  sol <- solve(prob, Tsit5(), abstol = 1e-10, reltol = 1e-10)
  uend <- sol$u[[length(sol$t)]]
  expect_equal(dim(uend), c(2L, 2L))
  ex <- function(t) {
    e <- eigen(cx_A)
    Re(e$vectors %*% diag(exp(e$values * t)) %*% solve(e$vectors) %*% u0)
  }
  expect_equal(uend, ex(1), tolerance = 1e-8)
  expect_equal(sol(0.5), ex(0.5), tolerance = 1e-8)
  expect_equal(dim(as.matrix(sol)), c(length(sol$t), 4L))
})

test_that("RK4 is adaptive with the residual error estimate, fixed with adaptive = FALSE", {
  prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
  sol <- solve(prob, RK4())
  expect_equal(sol$retcode, "Success")
  expect_gt(sol$stats$nreject + sol$stats$naccept, 0)
  expect_equal(sol(1), exp(-1), tolerance = 1e-4)
  fixed <- solve(prob, RK4(), dt = 0.25, adaptive = FALSE)
  expect_equal(fixed$t, seq(0, 1, by = 0.25))
  expect_equal(fixed$stats$nf, 4L * 4L + 1L)
})

test_that("interpolant derivatives of Radau, Nordsieck BDF, MIRK and the generic fallback", {
  rf <- function(u, p, t) c(-0.5 * u[1] + sin(t), u[1] - 2 * u[2])
  prob <- ODEProblem(rf, c(1, 0.5), c(0, 2))
  ref <- solve(prob, Vern9(), abstol = 1e-13, reltol = 1e-13)
  d1 <- rf(ref(0.73), NULL, 0.73)
  for (alg in list(RadauIIA5(), NordsieckBDF())) {
    sol <- solve(prob, alg, abstol = 1e-8, reltol = 1e-8)
    expect_equal(sol(0.73, deriv = 1), d1, tolerance = 1e-4, label = alg$name)
    expect_length(sol(0.73, deriv = 3), 2L)
    expect_error(sol(0.73, deriv = 4), "Derivative order too high")
  }
  t <- t_nounits; D <- D_nounits
  cxy <- Variable("cxy", t)
  bsys <- mtkcompile(System(list(D(D(cxy)) ~ -cxy), t, name = "cxbvp",
                            constraints = list(EvalAt(0)(cxy) ~ 0, EvalAt(pi / 2)(cxy) ~ 1)))
  bs <- solve(BVProblem(bsys, list(), c(0, pi / 2)), MIRK4(), dt = 0.05)
  h <- 1e-6
  expect_equal(bs(0.4, deriv = 1), (bs(0.4 + h) - bs(0.4 - h)) / (2 * h), tolerance = 1e-6)
})

test_that("linear interpolation extrapolates outside the saved times, save_idxs selects components", {
  rf <- function(u, p, t) c(-0.5 * u[1] + sin(t), u[1] - 2 * u[2])
  sol <- solve(ODEProblem(rf, c(1, 0.5), c(0, 2)), Tsit5(), saveat = 0.5)
  U <- as.matrix(sol)
  slope1 <- (U[2, ] - U[1, ]) / 0.5
  expect_equal(sol(-0.25), U[1, ] - 0.25 * slope1)
  expect_equal(sol(-0.25, deriv = 1), slope1)
  slope_end <- (U[5, ] - U[4, ]) / 0.5
  expect_equal(sol(2.4), U[5, ] + 0.4 * slope_end)
  si <- solve(ODEProblem(rf, c(1, 0.5), c(0, 2)), Tsit5(), save_idxs = 2)
  expect_length(si(0.3), 1L)
  expect_equal(ncol(si(c(0.3, 0.6))), 1L)
  sl <- solve(ODEProblem(rf, c(1, 0.5), c(0, 2)), Tsit5(), save_idxs = 2, saveat = 0.5)
  expect_length(sl(0.3), 1L)
})

test_that("linear interpolation helpers handle decreasing times", {
  expect_equal(modelcompR:::.approx_hold(c(1, 0.5, 0), c(10, 5, 0), 0.25), 2.5)
  expect_equal(modelcompR:::.approx_hold(c(0, 0.5, 1), c(0, 5, 10), 2), 10)
})

test_that("get_du has no derivative for a FunctionMap stepper", {
  integ <- init(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), Tsit5())
  integ$.I$alg <- FunctionMap()
  integ$.I$algs <- list(FunctionMap())
  expect_error(get_du(integ), "Derivatives are not defined for this stepper.")
})

test_that("SDE and DDE solves check their keywords", {
  sp <- SDEProblem(function(u, p, t) -u, function(u, p, t) 0.1, 1, c(0, 1))
  expect_error(solve(sp, SRIW1(), bogus_kw = 1), class = "CommonKwargError")
  expect_error(init(sp, SRIW1(), bogus_kw = 1), class = "CommonKwargError")
  expect_warning(solve(sp, SRIW1(), bogus_kw = 1, kwargshandle = "warn", seed = 1), "Unrecognized keyword")
  dp <- DDEProblem(function(u, h, p, t) -h(p, t - 1), 1, function(p, t) 1, c(0, 2), constant_lags = 1)
  expect_error(solve(dp, MethodOfSteps(Tsit5()), bogus_kw = 1), class = "CommonKwargError")
  expect_equal(solve(dp, MethodOfSteps(Tsit5()), discontinuity_interp_points = 5)$retcode, "Success")
})

test_that("complex Landau-Zener steps, stats and dense output match OrdinaryDiffEq", {
  # tools/complex-probe.jl (OrdinaryDiffEq 7.8.1): nf, naccept, nreject, then
  # Re and Im of sol(0.1) and of sol(0, Val{1})
  ref <- list(
    RK4 = list(c(39, 6, 0), c(0.9800667478447638, 0.0013278977808525714, -4.783384817798678e-9, -0.19866404926000028),
               c(-0.09983567117864944, 0.004991226679411247, -0.0001637021634408667, -0.9949916372776025)),
    DP5 = list(c(21, 3, 0), c(0.9800667532024202, 0.001327992601022011, -2.3465359863528157e-11, -0.19866401724611554),
               c(-0.09983274078328715, 0.00499152584493696, -0.00016641609750538765, -0.9949916617851498)),
    Tsit5 = list(c(27, 4, 0), c(0.980066755082459, 0.001327992453512979, 8.672073129045721e-11, -0.19866401774861042),
                 c(-0.09983283646166007, 0.004991609932318649, -0.0001665202935498709, -0.9949915486324429)),
    Vern7 = list(c(42, 4, 0), c(0.9800667551110478, 0.001327992420508986, 3.582316734918045e-13, -0.1986640177847921),
                 c(-0.09983283389313885, 0.004991650021317035, -0.00016649952229670054, -0.9949916889021381)),
    ImplicitEuler = list(c(147, 20, 4),
                         c(0.9776848864061006, 0.0013271783035259908, -0.0024104969757055717, -0.198135343205285),
                         c(-0.10317834275112926, 0.0041881335996154755, -0.003673886235039028, -0.9935043704533458)),
    Trapezoid = list(c(73, 12, 0), c(0.9800711722194663, 0.00135979072023623, 2.596198219341463e-6, -0.1986301273529706),
                     c(-0.09979295057018431, 0.00530196452672628, -0.0001440441913630186, -0.9946819935767525)),
    Kvaerno3 = list(c(71, 8, 0), c(0.9800679641937237, 0.0011461037100289448, 0.004946387690579716, -0.19866756387645532),
                    c(-0.09990380751494794, 0.0036333460853231357, 0.026931857344668908, -0.995004257807481)),
    Rosenbrock23 = list(c(102, 11, 0),
                        c(0.980068681092442, 0.0013050128719378114, -1.429282363101924e-6, -0.19865686211432437),
                        c(-0.09982638917299473, 0.004860299100356238, -0.00017373999828990713, -0.9949457835803587)))
  ad <- c("ImplicitEuler", "Trapezoid", "Kvaerno3", "Rosenbrock23")
  for (n in names(ref)) {
    alg <- if (n %in% ad) get(n)(autodiff = AutoFiniteDiff()) else get(n)()
    sol <- solve(ODEProblem(cx_fun, c(1 + 0i, 0), c(-cx_T, cx_T)), alg)
    expect_equal(c(sol$stats$nf, sol$stats$naccept, sol$stats$nreject), ref[[n]][[1]], label = n)
    z <- sol(cx_T)
    expect_equal(c(Re(z), Im(z)), ref[[n]][[2]], tolerance = 1e-10, label = n)
    dz <- sol(0, deriv = 1)
    expect_equal(c(Re(dz), Im(dz)), ref[[n]][[3]], tolerance = 1e-10, label = n)
  }
})

test_that("interpolant derivatives and linear extrapolation match OrdinaryDiffEq", {
  rf <- function(u, p, t) c(-0.5 * u[1] + sin(t), u[1] - 2 * u[2])
  prob <- ODEProblem(rf, c(1, 0.5), c(0, 2))
  ref <- list(
    RadauIIA5 = c(0.9201623202516104, 0.46529542642995053, 0.2067866808274573, -0.010437885810577375,
                  0.6417427683598894, 0.2274428114864923, -0.981905473688492, 0.217384447333271),
    NordsieckBDF = c(0.9201623243662688, 0.46529534489118696, 0.20678847778018597, -0.010428415873450822,
                     0.641779843965809, 0.22764855576723642, -0.9877885070844494, 0.18677358359523932))
  for (n in names(ref)) {
    sol <- solve(prob, get(n)(), abstol = 1e-8, reltol = 1e-8)
    expect_equal(unlist(lapply(0:3, function(d) sol(0.73, deriv = d))), ref[[n]], tolerance = 1e-10, label = n)
    expect_error(sol(0.73, deriv = 4), "Derivative order too high")
  }
  sol <- solve(prob, Tsit5(), saveat = 0.5)
  expect_equal(sol(-0.25), c(1.0542271956037026, 0.5134177210319026), tolerance = 1e-12)
  expect_equal(sol(-0.25, deriv = 1), c(-0.21690878241481037, -0.05367088412761045), tolerance = 1e-12)
  expect_equal(sol(2.4), c(1.491739599741374, 0.6654722890906228), tolerance = 1e-12)
  expect_equal(sol(2.4, deriv = 1), c(0.33230040683844697, 0.15762931683670756), tolerance = 1e-12)
})
