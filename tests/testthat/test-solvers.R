# ODE and DAE solvers (R/solvers-*.R). References: SUNDIALS CVODE example
# cvRoberts_dns.out (examples/cvode/serial, Robertson kinetics), the Julia
# series in fixtures/julia (tools/julia/reference.jl: Robertson with Rodas5P,
# an index 1 DAE, Lotka Volterra), analytic solutions, and the asymptotic
# Van der Pol relaxation period (3 - 2 log 2) mu.

rober <- function(u, p, t) c(-0.04 * u[1] + 1e4 * u[2] * u[3],
                             0.04 * u[1] - 1e4 * u[2] * u[3] - 3e7 * u[2]^2,
                             3e7 * u[2]^2)

test_that("Robertson kinetics match the CVODE reference", {
  prob <- ODEProblem(rober, c(1, 0, 0), c(0, 4e4))
  # cvRoberts_dns.out: t = 0.4 and t = 4e4
  ref04 <- c(9.851641e-01, 3.386242e-05, 1.480205e-02)
  ref4e4 <- c(3.898902e-02, 1.622006e-07, 9.610108e-01)
  for (alg in list(Rodas5P(), Rosenbrock23(), FBDF())) {
    sol <- solve(prob, alg, reltol = 1e-5, abstol = 1e-10)
    expect_equal(sol$retcode, "Success")
    expect_equal(sol(0.4), ref04, tolerance = 2e-3, ignore_attr = TRUE)
    expect_equal(sol$u[[length(sol$t)]], ref4e4, tolerance = 2e-3, ignore_attr = TRUE)
    # mass conservation
    expect_lt(max(abs(rowSums(sol[]) - 1)), 1e-6)
    expect_gt(sol$stats$njacs, 0L)
  }
  # the CVODE root: y3 = 0.01 at t = 0.26391
  t_ev <- NULL
  cb <- ContinuousCallback(function(u, t, integ) u[3] - 0.01, function(integ) {
    t_ev <<- integ$t
    terminate(integ)
  })
  sol <- solve(prob, Rodas5P(), reltol = 1e-6, abstol = 1e-10, callback = cb)
  expect_equal(sol$retcode, "Terminated")
  # CVODE ran with rtol 1e-4 and printed five digits
  expect_equal(t_ev, 0.26391, tolerance = 1e-3)
})

test_that("Robertson as a mass matrix DAE and as a DAEProblem", {
  rober_mm <- function(u, p, t) c(-0.04 * u[1] + 1e4 * u[2] * u[3],
                                  0.04 * u[1] - 1e4 * u[2] * u[3] - 3e7 * u[2]^2,
                                  u[1] + u[2] + u[3] - 1)
  prob <- ODEProblem(rober_mm, c(1, 0, 0), c(0, 40), mass_matrix = diag(c(1, 1, 0)))
  sol <- solve(prob, Rodas5P(), reltol = 1e-6, abstol = 1e-10)
  expect_equal(sol$retcode, "Success")
  expect_equal(sol$u[[length(sol$t)]], c(7.158017e-01, 9.185037e-06, 2.841892e-01), tolerance = 2e-3,
               ignore_attr = TRUE)
  dae <- function(du, u, p, t) c(-0.04 * u[1] + 1e4 * u[2] * u[3] - du[1],
                                 0.04 * u[1] - 3e7 * u[2]^2 - 1e4 * u[2] * u[3] - du[2],
                                 u[1] + u[2] + u[3] - 1)
  dprob <- DAEProblem(dae, c(-0.04, 0.04, 0), c(1, 0, 0), c(0, 40), differential_vars = c(TRUE, TRUE, FALSE))
  dsol <- solve(dprob, DFBDF(), reltol = 1e-6, abstol = 1e-10)
  expect_s3_class(dsol, "DAESolution")
  expect_equal(dsol$retcode, "Success")
  expect_equal(dsol$u[[length(dsol$t)]], c(7.158017e-01, 9.185037e-06, 2.841892e-01), tolerance = 2e-3,
               ignore_attr = TRUE)
  expect_equal(dsol(4), c(9.055097e-01, 2.240338e-05, 9.446793e-02), tolerance = 2e-3, ignore_attr = TRUE)
})

test_that("Van der Pol with mu = 1000 relaxes with the asymptotic period", {
  mu <- 1000
  # MATLAB's vdp1000: x'' = mu (1 - x^2) x' - x
  vdp <- function(u, p, t) c(u[2], p * (1 - u[1]^2) * u[2] - u[1])
  vdp_jac <- function(u, p, t) matrix(c(0, 1, -2 * p * u[1] * u[2] - 1, p * (1 - u[1]^2)), 2, byrow = TRUE)
  prob <- ODEProblem(ODEFunction(vdp, jac_fn = vdp_jac), u0 = c(2, 0), tspan = c(0, 1000), p = mu)
  half <- (3 / 2 - log(2)) * mu
  t_cross <- list()
  for (alg in list(Rodas5P(), FBDF())) {
    tc <- NULL
    cb <- ContinuousCallback(function(u, t, integ) u[1], affect = NULL,
                             affect_neg = function(integ) if (is.null(tc)) tc <<- integ$t)
    sol <- solve(prob, alg, reltol = 1e-5, abstol = 1e-7, callback = cb)
    expect_equal(sol$retcode, "Success")
    expect_equal(tc, half, tolerance = 0.01)
    t_cross[[alg$name]] <- tc
    # the slow manifold: x stays within [-2, 2] up to the tolerance
    expect_lt(max(abs(sol[1])), 2.01)
  }
  expect_equal(t_cross$Rodas5P, t_cross$FBDF, tolerance = 1e-3)
  # the default algorithm detects the stiffness and switches
  sol <- solve(ODEProblem(vdp, c(2, 0), c(0, 20), p = mu), reltol = 1e-4)
  expect_equal(sol$retcode, "Success")
  expect_true(3L %in% sol$alg_choice)
})

test_that("explicit methods converge with their order on y' = y", {
  prob <- ODEProblem(function(u, p, t) u, 1, c(0, 2))
  err <- function(alg, dt) abs(tail(solve(prob, alg, dt = dt, adaptive = FALSE)[1], 1) - exp(2))
  # steps in the asymptotic range (Tsit5's optimized error constants make
  # its observed order irregular for larger steps)
  for (case in list(list(Euler(), 1, 0.01), list(BS3(), 3, 1 / 16), list(RK4(), 4, 0.1),
                    list(Tsit5(), 5, 1 / 16), list(DP5(), 5, 1 / 8), list(Vern7(), 7, 1 / 4),
                    list(Vern9(), 9, 1 / 2))) {
    e1 <- err(case[[1]], case[[3]])
    e2 <- err(case[[1]], case[[3]] / 2)
    expect_equal(log2(e1 / e2), case[[2]], tolerance = 0.1, label = case[[1]]$name)
  }
  # adaptive solutions and their dense output
  for (alg in list(BS3(), DP5(), Tsit5(), Vern7(), Vern9())) {
    sol <- solve(prob, alg, reltol = 1e-9, abstol = 1e-12)
    tol <- if (alg$name == "BS3") 1e-6 else 1e-7
    expect_equal(tail(sol[1], 1), exp(2), tolerance = tol, label = alg$name)
    tq <- c(0.3, 1.1, 1.77)
    expect_equal(sol(tq)[, 1], exp(tq), tolerance = tol, ignore_attr = TRUE, label = alg$name)
  }
})

test_that("implicit and Rosenbrock methods are accurate on a stiff linear problem", {
  f <- function(u, p, t) -50 * (u - cos(t))
  exact <- function(t) {
    c1 <- 2500 / 2501; c2 <- 50 / 2501
    c1 * cos(t) + c2 * sin(t) - c1 * exp(-50 * t)
  }
  prob <- ODEProblem(f, 0, c(0, 2))
  for (alg in list(Rosenbrock23(), Rodas4(), Rodas5P(), ImplicitEuler(), Trapezoid(), FBDF())) {
    rt <- switch(alg$name, ImplicitEuler = 1e-4, Trapezoid = 1e-5, 1e-6)
    tol <- switch(alg$name, ImplicitEuler = 1e-2, Trapezoid = , Rosenbrock23 = 1e-3, 1e-5)
    sol <- solve(prob, alg, reltol = rt, abstol = rt / 100)
    expect_equal(sol$retcode, "Success", label = alg$name)
    expect_equal(tail(sol[1], 1), exact(2), tolerance = tol, label = alg$name)
    expect_equal(sol(c(0.5, 1.5))[, 1], exact(c(0.5, 1.5)), tolerance = 10 * tol, ignore_attr = TRUE,
                 label = alg$name)
  }
  # Jacobians are reused across steps by the Newton based methods
  sol <- solve(prob, FBDF(), reltol = 1e-6, abstol = 1e-8)
  expect_lt(sol$stats$njacs, sol$stats$naccept)
  expect_gt(sol$stats$nsolve, 0L)
})

test_that("an index 1 mass matrix DAE from mtkcompile matches its analytic solution", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t)
  # 0 = y^3 + y + x^3 + x forces y = -x, so x = exp(-t)
  sys <- mtkcompile(System(list(D(x) ~ y, 0 ~ y^3 + y + x^3 + x), t, name = "dae1"))
  M <- calculate_massmatrix(sys)
  expect_true(any(rowSums(abs(M)) == 0))
  prob <- ODEProblem(sys, list(x ~ 1, y ~ -1), c(0, 2), jac = TRUE)
  for (alg in list(Rodas5P(), Rodas4(), Rosenbrock23(), FBDF())) {
    sol <- solve(prob, alg, reltol = 1e-7, abstol = 1e-9)
    tol <- if (alg$name == "Rosenbrock23") 1e-3 else 1e-5
    expect_equal(sol$retcode, "Success", label = alg$name)
    expect_equal(tail(sol[x], 1), exp(-2), tolerance = tol, label = alg$name)
    expect_equal(tail(sol[y], 1), -exp(-2), tolerance = tol, label = alg$name)
    expect_equal(sol(1.3, idxs = x), exp(-1.3), tolerance = 10 * tol, label = alg$name)
  }
  # inconsistent algebraic initial values are repaired (BrownFullBasicInit)
  f <- function(u, p, t) c(u[2], u[2]^3 + u[2] + u[1]^3 + u[1])
  prob2 <- ODEProblem(f, c(1, 0.3), c(0, 1), mass_matrix = diag(c(1, 0)))
  sol <- solve(prob2, Rodas5P(), reltol = 1e-7, abstol = 1e-9)
  expect_equal(sol$u[[1]], c(1, -1), tolerance = 1e-8)
  expect_equal(sol$u[[length(sol$t)]], c(exp(-1), -exp(-1)), tolerance = 1e-5)
  sol <- solve(prob2, Rodas5P(), initializealg = ShampineCollocationInit())
  expect_equal(sol$retcode, "Success")
  expect_warning(sol <- solve(prob2, Rodas5P(), initializealg = CheckInit()), "consistent")
  expect_equal(sol$retcode, "InitialFailure")
  # the default algorithm uses Rodas5P for mass matrix problems
  sol <- solve(prob)
  expect_true(all(sol$alg_choice == 4L))
})

test_that("DAEProblem from a system is solved by DFBDF and DImplicitEuler", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t)
  sys <- mtkcompile(System(list(D(x) ~ y, 0 ~ y^3 + y + x^3 + x), t, name = "dae2"))
  prob <- DAEProblem(sys, list(x ~ 1, y ~ -1, D(x) ~ -1), c(0, 2))
  expect_equal(prob$differential_vars, vapply(unknowns(sys), function(u) .key(u) == .key(x), TRUE))
  sol <- solve(prob, DFBDF(), reltol = 1e-7, abstol = 1e-9)
  expect_s3_class(sol, "DAESolution")
  expect_equal(sol$retcode, "Success")
  expect_equal(tail(sol[x], 1), exp(-2), tolerance = 1e-5)
  expect_equal(tail(sol[y], 1), -exp(-2), tolerance = 1e-5)
  expect_equal(sol(1, idxs = x), exp(-1), tolerance = 1e-4)
  sol <- solve(prob, DImplicitEuler(), reltol = 1e-4, abstol = 1e-6)
  expect_equal(sol$retcode, "Success")
  expect_equal(tail(sol[x], 1), exp(-2), tolerance = 2e-2)
  # analytic Jacobian dF/du + gamma dF/d(du)
  J <- generate_dae_jacobian(sys)
  u <- c(0.5, -0.5)[match(names(prob$u0), c("x(t)", "y(t)"))]
  Jn <- J(c(0, 0), u, prob$p, 2, 0)
  expect_equal(dim(Jn), c(2L, 2L))
  # default solver and inconsistent du0: BrownFullBasicInit
  prob2 <- DAEProblem(sys, list(x ~ 1, y ~ -1), c(0, 2), jac = TRUE)
  sol <- solve(prob2, reltol = 1e-7, abstol = 1e-9)
  expect_equal(sol$alg$name, "DFBDF")
  expect_equal(tail(sol[x], 1), exp(-2), tolerance = 1e-5)
})

test_that("a bouncing ball bounces at the analytic times with every method", {
  g <- 9.81
  ball <- function(u, p, t) c(u[2], -g)
  prob <- ODEProblem(ball, c(1, 0), c(0, 2))
  t1 <- sqrt(2 / g)
  bounce <- function(hits_env) ContinuousCallback(function(u, t, integ) u[1], function(integ) {
    hits_env$t <- c(hits_env$t, integ$t)
    integ$u[2] <- -integ$u[2]
  })
  for (alg in list(Tsit5(), BS3(), DP5(), Vern7(), Vern9(), Rosenbrock23(), Rodas4(), Rodas5P(), FBDF(),
                   Trapezoid(), ImplicitEuler())) {
    h <- new.env()
    rt <- if (alg$name == "ImplicitEuler") 1e-4 else 1e-8
    sol <- solve(prob, alg, callback = bounce(h), reltol = rt, abstol = rt / 100)
    expect_equal(sol$retcode, "Success", label = alg$name)
    expect_equal(h$t[1:2], c(t1, 3 * t1), tolerance = if (rt > 1e-6) 1e-2 else 1e-6, label = alg$name)
  }
  # the same ball as a mass matrix DAE (0 = a + g) and as a fully implicit DAE:
  # the algebraic variable is reinitialized after each bounce
  mm <- ODEProblem(function(u, p, t) c(u[2], u[3], u[3] + g), c(1, 0, 0), c(0, 2),
                   mass_matrix = diag(c(1, 1, 0)))
  h <- new.env()
  sol <- solve(mm, Rodas5P(), callback = bounce(h), reltol = 1e-8, abstol = 1e-10)
  expect_equal(sol$retcode, "Success")
  expect_equal(h$t[1:2], c(t1, 3 * t1), tolerance = 1e-6)
  expect_equal(sol$u[[length(sol$t)]][3], -g, tolerance = 1e-8)
  dae <- DAEProblem(function(du, u, p, t) c(du[1] - u[2], du[2] - u[3], u[3] + g), c(0, -g, 0), c(1, 0, -g),
                    c(0, 2), differential_vars = c(TRUE, TRUE, FALSE))
  h <- new.env()
  sol <- solve(dae, DFBDF(), callback = bounce(h), reltol = 1e-8, abstol = 1e-10)
  expect_equal(sol$retcode, "Success")
  expect_equal(h$t[1:2], c(t1, 3 * t1), tolerance = 1e-5)
})

test_that("the default algorithm follows OrdinaryDiffEqDefault", {
  prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
  sol <- solve(prob)
  expect_true(all(sol$alg_choice == 1L))
  expect_equal(tail(sol[1], 1), exp(-1), tolerance = 1e-4)
  sol <- solve(prob, reltol = 1e-8)
  expect_true(all(sol$alg_choice == 2L))
  # stiff: switches from Tsit5 to Rosenbrock23 (Rodas5P at tight tolerances)
  rp <- ODEProblem(rober, c(1, 0, 0), c(0, 10))
  sol <- solve(rp)
  expect_equal(sol$alg_choice[1], 1L)
  expect_true(3L %in% sol$alg_choice)
  expect_equal(sol(10), solve(rp, Rodas5P(), reltol = 1e-8, abstol = 1e-12)(10), tolerance = 1e-2,
               ignore_attr = TRUE)
  # choices of stiffchoice / nonstiffchoice
  expect_equal(.stiffchoice(1e-3, 3, NULL), 3L)
  expect_equal(.stiffchoice(1e-8, 3, NULL), 4L)
  expect_equal(.stiffchoice(1e-3, 3, diag(c(1, 0, 1))), 4L)
  expect_equal(.stiffchoice(1e-3, 60, NULL), 5L)
  expect_equal(.nonstiffchoice(1e-7), 2L)
  big <- ODEProblem(function(u, p, t) c(-u[1:59], sum(u[1:59]) - u[60]), rep(1, 60), c(0, 0.1),
                    mass_matrix = diag(c(rep(1, 59), 0)))
  sol <- solve(big)
  expect_true(all(sol$alg_choice == 5L))
  expect_equal(sol$u[[length(sol$t)]][60], 59 * exp(-0.1), tolerance = 1e-3)
  expect_s3_class(AutoTsit5(Rosenbrock23()), "CompositeAlgorithm")
  sol <- solve(rp, AutoTsit5(Rosenbrock23()))
  expect_true(all(c(1L, 2L) %in% sol$alg_choice))
})

test_that("retcodes and statistics follow SciMLBase", {
  prob <- ODEProblem(function(u, p, t) u^2, 1, c(0, 2))
  expect_warning(sol <- solve(prob, Tsit5()), "dtmin")
  expect_equal(sol$retcode, "DtLessThanMin")
  expect_warning(sol <- solve(ODEProblem(function(u, p, t) -u, 1, c(0, 100)), Tsit5(), maxiters = 3), "maximum")
  expect_equal(sol$retcode, "MaxIters")
  expect_warning(sol <- solve(ODEProblem(function(u, p, t) NaN * u, 1, c(0, 1)), Euler(), dt = 0.1), "non finite")
  expect_equal(sol$retcode, "Unstable")
  sol <- solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), Rodas5P())
  expect_named(sol$stats, c("nf", "njacs", "nsolve", "nw", "naccept", "nreject"), ignore.order = TRUE)
  expect_equal(sol$stats$naccept, length(sol$t) - 1L)
  cb <- DiscreteCallback(function(u, t, integ) t >= 0.5, function(integ) terminate(integ))
  sol <- solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), FBDF(), callback = cb)
  expect_equal(sol$retcode, "Terminated")
  # a domain error in f rejects steps instead of failing
  for (alg in list(Rosenbrock23(), Rodas4(), Rodas5P(), Trapezoid(), FBDF())) {
    sol <- suppressWarnings(solve(ODEProblem(function(u, p, t) -sqrt(u), 1, c(0, 3)), alg))
    expect_true(sol$retcode %in% c("Success", "DtLessThanMin", "Unstable"), label = alg$name)
  }
  # a right-hand side that turns NaN (as a failed assertion does) rejects the
  # steps of every adaptive path until dt falls below its minimum
  nanf <- ODEProblem(function(u, p, t) if (t > 0.5) c(NaN, NaN) else c(-u[1], -1000 * (u[2] - u[1])),
                     c(1, 0), c(0, 1))
  for (alg in list(DefaultODEAlgorithm(), AutoTsit5(Rosenbrock23()), Tsit5(), Rosenbrock23(), Rodas5P(),
                   FBDF(), ImplicitEuler(), Trapezoid())) {
    expect_warning(sol <- solve(nanf, alg), "dtmin", label = alg$name)
    expect_equal(sol$retcode, "DtLessThanMin", label = alg$name)
    expect_equal(tail(sol$t, 1), 0.5, tolerance = 1e-6, label = alg$name)
  }
  nand <- DAEProblem(function(du, u, p, t) if (t > 0.5) NaN else du + u, -1, 1, c(0, 1), differential_vars = TRUE)
  expect_warning(sol <- solve(nand, DFBDF()), "dtmin")
  expect_equal(sol$retcode, "DtLessThanMin")
  # a constant state does not count as stiff for the default algorithm
  sol <- solve(ODEProblem(function(u, p, t) c(-u[1], 0), c(1, 1), c(0, 5)))
  expect_true(all(sol$alg_choice == 1L))
  bad <- ImplicitDiscreteProblem(function(un, u, p, t) un^2 + 1, 1, c(0, 2))
  expect_equal(suppressWarnings(solve(bad, IDSolve()))$retcode, "InitialFailure")
})

test_that("implicit discrete problems, semilinear problems and generate_W", {
  # implicit map: u_next + u_next^3 = u
  prob <- ImplicitDiscreteProblem(function(un, u, p, t) un + un^3 - u, 2, c(0, 3))
  sol <- solve(prob, IDSolve())
  us <- unlist(sol$u)
  expect_equal(length(us), 4L)
  expect_equal(us[-1] + us[-1]^3, us[-4], tolerance = 1e-10)
  expect_equal(us[1] + us[1]^3, 2, tolerance = 1e-10)
  # from a compiled discrete system with an algebraic equation
  t <- t_nounits
  x <- Variable("x", t)
  k <- ShiftIndex(t)
  sys <- mtkcompile(System(list(0 ~ x^3 + x - (k - 1)(x)), t, name = "imp"))
  dprob <- ImplicitDiscreteProblem(sys, list((k - 1)(x) ~ 2, x ~ 1), c(0, 3))
  dsol <- solve(dprob, IDSolve())
  xs <- dsol[x]
  expect_equal(xs[1], 1, tolerance = 1e-10)
  expect_equal(xs[-1]^3 + xs[-1], xs[-length(xs)], tolerance = 1e-8)
  # semilinear split: A x + x' B x + C
  D <- D_nounits
  y <- Variable("y", t); z <- Variable("z", t); a <- Parameter("a")
  s2 <- mtkcompile(System(list(D(y) ~ -a * y + y * z, D(z) ~ sin(y) - 2 * z), t, name = "semi"))
  sq <- calculate_semiquadratic_form(s2)
  expect_false(is.null(sq$A)); expect_false(is.null(sq$B)); expect_false(is.null(sq$C))
  sp <- SemilinearODEProblem(s2, list(y ~ 1, z ~ 0.5, a ~ 3), c(0, 1))
  expect_s3_class(sp, "ODEProblem")
  u <- sp$u0
  full <- ODEProblem(s2, list(y ~ 1, z ~ 0.5, a ~ 3), c(0, 1))
  expect_equal(sp$f$f1(u, sp$p, 0) + sp$f$f2(u, sp$p, 0), full$f$f(u, full$p, 0), ignore_attr = TRUE)
  iy <- match("y(t)", names(u))
  expect_equal(sp$f$f1(u, sp$p, 0)[iy], -3 * u[["y(t)"]], ignore_attr = TRUE)
  s1 <- solve(sp, Rodas5P(), reltol = 1e-8); s0 <- solve(full, Rodas5P(), reltol = 1e-8)
  expect_equal(tail(s1[y], 1), tail(s0[y], 1), tolerance = 1e-7)
  # W = gamma M + J
  W <- generate_W(s2)
  J <- generate_jacobian(s2)
  expect_equal(W(u, sp$p, -2, 0), -2 * diag(2) + J(u, sp$p, 0))
})

test_that("compiled systems match the Julia references (tools/julia/reference.jl)", {
  fx <- function(name) {
    f <- test_path("fixtures", "julia", paste0(name, ".csv"))
    skip_if_not(file.exists(f), "Julia fixtures not generated")
    utils::read.csv(f, check.names = FALSE, encoding = "UTF-8")
  }
  t <- t_nounits; D <- D_nounits
  # Robertson, Rodas5P at tolerance 1e-10 in Julia
  ref <- fx("robertson")
  y1 <- Variable("y1", t); y2 <- Variable("y2", t); y3 <- Variable("y3", t)
  k1 <- Parameter("k1", default = 0.04); k2 <- Parameter("k2", default = 3e7)
  k3 <- Parameter("k3", default = 1e4)
  rober_sys <- mtkcompile(System(list(D(y1) ~ -k1 * y1 + k3 * y2 * y3,
                                      D(y2) ~ k1 * y1 - k2 * y2^2 - k3 * y2 * y3,
                                      D(y3) ~ k2 * y2^2), t, name = "rober"))
  prob <- ODEProblem(rober_sys, list(y1 ~ 1, y2 ~ 0, y3 ~ 0), c(0, 1e5), jac = TRUE)
  for (alg in list(Rodas5P(), FBDF())) {
    sol <- solve(prob, alg, reltol = 1e-7, abstol = 1e-12, saveat = ref$t)
    expect_equal(sol$t, ref$t)
    expect_equal(sol[y1], ref[["y1(t)"]], tolerance = 1e-4, label = alg$name)
    expect_equal(sol[y2], ref[["y2(t)"]], tolerance = 1e-4, label = alg$name)
    expect_equal(sol[y3], ref[["y3(t)"]], tolerance = 1e-4, label = alg$name)
  }
  # index 1 DAE kept with a singular mass matrix; b starts from its guess
  ref <- fx("dae1")
  a <- Variable("a", t); b <- Variable("b", t, guess = 0.5); cc <- Parameter("c", default = 2)
  dae1 <- mtkcompile(System(list(D(a) ~ -a + b, 0 ~ b^3 + b - cc * a), t, name = "dae1"))
  sol <- solve(ODEProblem(dae1, list(a ~ 1), c(0, 3)), Rodas5P(), reltol = 1e-8, abstol = 1e-10,
               saveat = ref$t)
  expect_equal(sol[a], ref[["a(t)"]], tolerance = 1e-8)
  expect_equal(sol[b], ref[["b(t)"]], tolerance = 1e-8)
  # Lotka Volterra with the explicit methods
  ref <- fx("lotka_volterra")
  prey <- Variable("prey", t); pred <- Variable("pred", t)
  al <- Parameter("alpha", default = 1.5); be <- Parameter("beta", default = 1)
  ga <- Parameter("gamma", default = 3); de <- Parameter("delta", default = 1)
  lv <- mtkcompile(System(list(D(prey) ~ al * prey - be * prey * pred,
                               D(pred) ~ -ga * pred + de * prey * pred), t, name = "lv"))
  lprob <- ODEProblem(lv, list(prey ~ 1, pred ~ 1), c(0, 10))
  for (alg in list(Tsit5(), DP5(), Vern7(), Vern9())) {
    sol <- solve(lprob, alg, reltol = 1e-9, abstol = 1e-11)
    expect_equal(sol(ref$t, idxs = prey), ref[["prey(t)"]], tolerance = 1e-6, label = alg$name)
    expect_equal(sol(ref$t, idxs = pred), ref[["pred(t)"]], tolerance = 1e-6, label = alg$name)
  }
})

solvers_stats_ref <- function() {
  f <- test_path("fixtures", "julia", "solver_stats.csv")
  skip_if_not(file.exists(f), "Julia solver statistics not generated")
  utils::read.csv(f)
}

# naccept within 15% of OrdinaryDiffEq's, nreject within 50% or 10 (the
# R DAE problems use finite difference Jacobians, Julia ForwardDiff). The
# measured accepted counts are within 4%.
expect_steps_like_julia <- function(sol, case) {
  ref <- solvers_stats_ref()
  r <- ref[ref$case == case, ]
  expect_equal(nrow(r), 1L, label = case)
  expect_lt(abs(sol$stats$naccept - r$naccept), 0.15 * r$naccept, label = case)
  expect_lte(abs(sol$stats$nreject - r$nreject), max(10, 0.5 * r$nreject), label = case)
}

rober_jac <- function(u, p, t) matrix(c(-0.04, 1e4 * u[3], 1e4 * u[2],
                                        0.04, -1e4 * u[3] - 6e7 * u[2], -1e4 * u[2],
                                        0, 6e7 * u[2], 0), 3, byrow = TRUE)

test_that("FBDF and DFBDF take OrdinaryDiffEq's steps on Robertson and Van der Pol", {
  # Robertson on [0, 1e5]: OrdinaryDiffEq takes 385 steps (18 rejected)
  prob <- ODEProblem(ODEFunction(rober, jac_fn = rober_jac), c(1, 0, 0), c(0, 1e5))
  sol <- solve(prob, FBDF(), reltol = 1e-6, abstol = 1e-10)
  expect_steps_like_julia(sol, "rober_FBDF_1e-6")
  expect_equal(sol$u[[length(sol$t)]], c(1.786592e-02, 7.274751e-08, 9.821340e-01), tolerance = 1e-5,
               ignore_attr = TRUE)
  # out of place, as R functions are, every new W comes with a fresh
  # Jacobian (calc_W of OrdinaryDiffEqDifferentiation)
  expect_equal(sol$stats$njacs, sol$stats$nw)
  expect_gt(sol$stats$nw, 50L)
  # the time filters of DeCaria et al.
  solf <- solve(prob, FBDF(time_filter = TRUE), reltol = 1e-6, abstol = 1e-10)
  expect_steps_like_julia(solf, "rober_FBDF_tf_1e-6")
  expect_equal(solf$u[[length(solf$t)]], sol$u[[length(sol$t)]], tolerance = 1e-4)
  # a smaller largest growth factor
  solq <- solve(prob, FBDF(qmax = 2), reltol = 1e-6, abstol = 1e-10)
  expect_steps_like_julia(solq, "rober_FBDF_qmax2_1e-6")
  # Van der Pol, mu = 1000, on [0, 1000]
  vdp <- function(u, p, t) c(u[2], p * (1 - u[1]^2) * u[2] - u[1])
  vdp_jac <- function(u, p, t) matrix(c(0, 1, -2 * p * u[1] * u[2] - 1, p * (1 - u[1]^2)), 2, byrow = TRUE)
  vp <- ODEProblem(ODEFunction(vdp, jac_fn = vdp_jac), c(2, 0), c(0, 1000), p = 1000)
  sv <- solve(vp, FBDF(), reltol = 1e-5, abstol = 1e-7)
  expect_steps_like_julia(sv, "vdp_FBDF_1e-5")
  # the fully implicit forms with DFBDF
  dae <- function(du, u, p, t) c(-0.04 * u[1] + 1e4 * u[2] * u[3] - du[1],
                                 0.04 * u[1] - 3e7 * u[2]^2 - 1e4 * u[2] * u[3] - du[2],
                                 u[1] + u[2] + u[3] - 1)
  dp <- DAEProblem(dae, c(-0.04, 0.04, 0), c(1, 0, 0), c(0, 1e5), differential_vars = c(TRUE, TRUE, FALSE))
  sd <- solve(dp, DFBDF(), reltol = 1e-6, abstol = 1e-10)
  expect_steps_like_julia(sd, "roberdae_DFBDF_1e-6")
  expect_equal(sd$u[[length(sd$t)]], c(1.786592e-02, 7.274751e-08, 9.821340e-01), tolerance = 1e-4,
               ignore_attr = TRUE)
  vdae <- DAEProblem(function(du, u, p, t) c(u[2] - du[1], p * (1 - u[1]^2) * u[2] - u[1] - du[2]),
                     c(0, -2), c(2, 0), c(0, 1000), 1000, differential_vars = c(TRUE, TRUE))
  sdv <- solve(vdae, DFBDF(), reltol = 1e-5, abstol = 1e-7)
  expect_steps_like_julia(sdv, "vdpdae_DFBDF_1e-5")
  expect_equal(sdv$u[[length(sdv$t)]], sv$u[[length(sv$t)]], tolerance = 1e-3)
  expect_error(FBDF(max_order = 6), "between 1 and 5")
  expect_error(solve(ODEProblem(rober, c(1, 0, 0), c(0, 1), mass_matrix = diag(3)), FBDF(time_filter = TRUE)),
               NA)
  expect_error(solve(ODEProblem(function(u, p, t) c(u[2], u[1] - 1), c(1, 1), c(0, 1),
                                mass_matrix = diag(c(1, 0))), FBDF(time_filter = TRUE)), "identity mass matrix")
})

test_that("stability limit detection matches OrdinaryDiffEqBDF's _stald_detect", {
  # random and recorded scaled derivative data with Julia's verdicts (flags
  # 4 to 6 order a reduction); tools/solvers-julia-references.jl
  f <- test_path("fixtures", "julia", "stald_detect.csv")
  skip_if_not(file.exists(f), "STALD fixtures not generated")
  d <- utils::read.csv(f)
  flags <- vapply(seq_len(nrow(d)), function(i) {
    st <- .stald_new()
    st$ssdat <- matrix(as.numeric(d[i, 2:16]), 5L, 3L)
    .stald_detect(st, d$q[i])
  }, 0L)
  expect_equal(flags, d$flag)
  expect_true(any(d$flag > 3L))
  # degenerate (zero) data: comparisons with NaN are false, as in Julia
  expect_equal(.stald_detect(.stald_new(), 3L), -4L)
  # in the controller: a detected violation at order 5 lowers the order
  viol <- which(d$flag > 3L & d$q == 5L)[1L]
  cache <- .bdf_init(list(f = function(u, p, t) -u, p = NULL), 0, 1, 5L, FALSE,
                     list(stald = list(enabled = TRUE)))
  cache$order <- 5L
  cache$qwait <- 2L
  cache$terkm2 <- 1e-3; cache$terkm1 <- 1e-4; cache$terk <- 1e-5; cache$terkp1 <- 0
  cache$step_info <- list(t = 1, dt = 0.1, u = 1, uprev = 1, ctx = list(abstol = 1e-6, reltol = 1e-3))
  st <- cache$stald
  st$last_order <- 5L
  st$nscon <- 20L
  ss <- matrix(as.numeric(d[viol, 2:16]), 5L, 3L)
  # the controller adds the step's norms as the newest row
  st$ssdat <- rbind(ss[2:5, ], 0)
  cache$terkm2 <- sqrt(ss[1, 1]); cache$terkm1 <- sqrt(ss[1, 2]); cache$terk <- sqrt(ss[1, 3])
  co <- .ctrl_opts(FBDF(), list())
  q <- .bdf_stepsize(cache, 0.5, co, FALSE)
  expect_equal(st$nor, 1L)
  expect_equal(cache$order, 4L)
  expect_true(is.finite(q))
  # disabled detection keeps the order
  cache$stald$enabled <- FALSE
  cache$order <- 5L
  .bdf_stepsize(cache, 0.5, co, FALSE)
  expect_equal(cache$order, 5L)
  # a restart clears the data
  .bdf_restart(cache)
  expect_equal(sum(cache$stald$ssdat), 0)
})

test_that("step limiters change the accepted state", {
  # project the harmonic oscillator onto the unit circle after every step
  osc <- ODEProblem(function(u, p, t) c(u[2], -u[1]), c(1, 0), c(0, 20))
  proj <- function(u, integ, p, t) u / sqrt(sum(u^2))
  for (alg in list(Tsit5(), Rodas5P(), FBDF())) {
    sol <- solve(osc, alg, reltol = 1e-3, abstol = 1e-6, step_limiter = proj)
    expect_equal(rowSums(sol[]^2), rep(1, length(sol$t)), tolerance = 1e-12, label = alg$name)
  }
  # the algorithm field of FBDF (deprecated in Julia in favor of the option)
  sol <- solve(osc, FBDF(step_limiter = proj), reltol = 1e-3, abstol = 1e-6)
  expect_equal(rowSums(sol[]^2), rep(1, length(sol$t)), tolerance = 1e-12)
  free <- solve(osc, FBDF(), reltol = 1e-3, abstol = 1e-6)
  expect_gt(max(abs(rowSums(free[]^2) - 1)), 1e-6)
  expect_error(FBDF(step_limiter = 1), "function")
})

test_that("interpolation only Verner stages count as function evaluations", {
  # 1 for the derivative at the start, 2 for the initial step size estimate
  prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
  for (case in list(list(Vern7(), 10L, 6L), list(Vern9(), 16L, 10L))) {
    alg <- case[[1]]
    sol <- solve(prob, alg, reltol = 1e-8)
    n <- sol$stats$naccept + sol$stats$nreject
    expect_equal(sol$stats$nf, 3L + n * case[[2]], label = alg$name)
    # the extra stages of a step are computed once, on the first interpolation
    sol(0.05)
    expect_equal(sol$stats$nf, 3L + n * case[[2]] + case[[3]], label = alg$name)
    sol(0.051)
    expect_equal(sol$stats$nf, 3L + n * case[[2]] + case[[3]], label = alg$name)
    # saveat interpolates during the solve, once per step it falls in
    sv <- solve(prob, alg, reltol = 1e-8, saveat = c(0.05, 0.5))
    nsteps <- length(unique(findInterval(c(0.05, 0.5), sol$t, left.open = TRUE)))
    expect_equal(sv$stats$nf, 3L + n * case[[2]] + nsteps * case[[3]], label = alg$name)
  }
})

test_that("remake() rebuilds DAE problems", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("rmdx", t); y <- Variable("rmdy", t); k <- Parameter("rmdk", default = 1)
  sys <- mtkcompile(System(list(D(x) ~ y, 0 ~ y^3 + y + k * (x^3 + x)), t, name = "rmdae"))
  prob <- DAEProblem(sys, list(x ~ 1, y ~ -1, D(x) ~ -1), c(0, 1))
  ix <- match("rmdx(t)", names(prob$u0))
  # a new initial value and derivative
  p2 <- remake(prob, u0 = list(x ~ 0.5, D(x) ~ -0.5))
  expect_s3_class(p2, "DAEProblem")
  expect_equal(unname(p2$u0[ix]), 0.5)
  expect_equal(unname(p2$du0[ix]), -0.5)
  expect_named(p2$du0, names(prob$u0))
  sol <- solve(p2, DFBDF(), reltol = 1e-8, abstol = 1e-10)
  expect_equal(sol$retcode, "Success")
  expect_equal(tail(sol[x], 1), 0.5 * exp(-1), tolerance = 1e-5)
  # parameters only: values and derivatives are kept; equal to building anew
  p3 <- remake(prob, p = list(k ~ 3))
  expect_equal(p3$p, DAEProblem(sys, list(x ~ 1, y ~ -1, D(x) ~ -1, k ~ 3), c(0, 1))$p)
  expect_equal(p3$du0, prob$du0)
  expect_equal(p3$u0, prob$u0)
  s3 <- solve(p3, DFBDF(), reltol = 1e-8, abstol = 1e-10)
  expect_equal(s3$retcode, "Success")
  expect_equal(tail(s3[y]^3 + s3[y] + 3 * (s3[x]^3 + s3[x]), 1), 0, tolerance = 1e-6)
  # an explicit du0 wins; tspan changes
  p4 <- remake(prob, tspan = c(0, 2), du0 = c(-2, 0))
  expect_equal(unname(p4$du0), c(-2, 0))
  expect_named(p4$du0, names(prob$u0))
  expect_equal(p4$tspan, c(0, 2))
  # a parameter remake keeps a state set numerically before
  p5 <- remake(remake(prob, u0 = c(0.8, -0.8)[order(names(prob$u0) != "rmdx(t)")]), p = list(k ~ 2))
  expect_equal(unname(p5$u0[ix]), 0.8)
  # problems without a system
  fp <- DAEProblem(function(du, u, p, t) du + p * u, -1, 1, c(0, 1), 1, differential_vars = TRUE)
  r <- remake(fp, p = 2, u0 = 3, du0 = -6)
  sr <- solve(r, DFBDF(), reltol = 1e-8, abstol = 1e-10)
  expect_equal(tail(sr[1], 1), 3 * exp(-2), tolerance = 1e-6)
})

test_that("Vern6 and Vern8 converge and interpolate", {
  p <- ODEProblem(function(u, p, t) c(u[2], -u[1]), c(0, 1), c(0, 20))
  for (alg in list(Vern6(), Vern8())) {
    e <- sapply(c(0.8, 0.4), function(dt) abs(tail(solve(p, alg, dt = dt, adaptive = FALSE)$u, 1)[[1]][1] - sin(20)))
    expect_gt(log2(e[1] / e[2]), alg$order - 0.5, label = alg$name)
    s <- solve(p, alg, reltol = 1e-10, abstol = 1e-12)
    tt <- seq(0.05, 19.95, by = 0.7)
    expect_lt(max(abs(vapply(tt, function(x) s(x)[1] - sin(x), 0))), 1e-8, label = alg$name)
  }
  for (alg in list(AutoVern6(Rodas5P()), AutoVern8(Rodas5P()), AutoVern9(Rodas5P())))
    expect_equal(tail(solve(p, alg, reltol = 1e-8, abstol = 1e-10)$u, 1)[[1]][1], sin(20), tolerance = 1e-5)
})
