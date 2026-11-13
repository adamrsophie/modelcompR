# PDESystem, debugging and assertions, units, dynamic optimization and the
# small API of the misc module. Reference values are analytic.

test_that("PDESystem collects the problem and builds analytic functions", {
  t <- independent_variable("t"); x <- independent_variable("x")
  u <- Variable("u", list(t, x)); a <- Parameter("a")
  Dt <- Differential(t); Dxx <- Differential(x, 2)
  eq <- Dt(u) ~ a * Dxx(u)
  bcs <- list(Function("u")(0, x) ~ sin(pi * x))
  pde <- PDESystem(eq, bcs, list(t ~ c(0, 1), VarDomainPairing(x, Interval(0, 1))), list(t, x), list(u),
                   list(a), analytic = list(u ~ exp(-a * pi^2 * t) * sin(pi * x)), name = "heat")
  expect_s3_class(pde, "PDESystem")
  expect_length(pde$eqs, 1L)
  expect_length(pde$bcs, 1L)
  expect_identical(format(pde$domain[[1L]]$domain), "Interval(0, 1)")
  expect_identical(format(pde$x), "x")
  expect_identical(format(pde$u), "u(t, x)")
  expect_identical(pde$ivs, pde$indvars)
  expect_length(equations(pde), 1L)
  f <- pde$analytic_func[[1L]]
  expect_equal(f(0.5, 0.1, 0.25), exp(-0.5 * pi^2 * 0.1) * sin(pi / 4))
  expect_output(print(pde), "Boundary Conditions")
  expect_error(PDESystem(eq, bcs, list(t ~ c(0, 1)), list(t, x), list(u), name = "p", analytic = list()),
               "number of analytic")
})

test_that("debug_system reports the failing call and assertions reject steps", {
  t <- t_nounits; D <- D_nounits
  u1 <- Variable("u1", t); u2 <- Variable("u2", t); u3 <- Variable("u3", t)
  eqs <- list(D(u1) ~ -sqrt(u1), D(u2) ~ -sqrt(u2), D(u3) ~ -sqrt(u3))
  ics <- list(u1 ~ 1, u2 ~ 2, u3 ~ 3)
  sys <- mtkcompile(System(eqs, t, initial_conditions = ics, name = "sys"))
  dsys <- debug_system(sys, functions = "sqrt")
  expect_error(solve(ODEProblem(dsys, list(), c(0, 10))),
               "Function sqrt\\(u1\\(t\\)\\) output non-finite value NaN with input\n  u1\\(t\\) => -")
  expect_error(debug_system(System(list(), t, name = "o", systems = list(System(eqs, t, name = "i")))),
               "no sub-systems")
  # u1 + u2 = (1 - t/2)^2 + (sqrt(2) - t/2)^2 reaches 2 at the root below
  s <- (1 + sqrt(2) - sqrt((1 + sqrt(2))^2 - 2)) / 2
  asys <- mtkcompile(System(eqs, t, initial_conditions = ics, name = "sys",
                            assertions = list(list(u1 + u2 >= 2, "Oh no!"))))
  prob <- ODEProblem(asys, list(), c(0, 10))
  expect_warning(sol <- solve(prob, Tsit5()), "eps\\(t\\)")
  expect_identical(sol$retcode, "Unstable")
  expect_equal(tail(sol$t, 1), 2 * s, tolerance = 1e-3)
  dprob <- ODEProblem(debug_system(asys, functions = character()), list(), c(0, 10))
  msgs <- character()
  withCallingHandlers(suppressWarnings(solve(dprob, Tsit5(), dtmin = 0.1)), message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  expect_true(length(msgs) > 0L && all(grepl("Assertion u1(t) + u2(t) >= 2 failed:\nOh no!", msgs, fixed = TRUE)))
  quiet <- remake(dprob, p = list(ASSERTION_LOG_VARIABLE ~ 0))
  expect_silent(suppressWarnings(solve(quiet, Tsit5(), dtmin = 0.1)))
  # divisions and powers are wrapped too
  x <- Variable("x", t)
  ds <- debug_system(mtkcompile(System(list(D(x) ~ 1 / (x - 1) + x^2), t, name = "d")))
  dp <- ODEProblem(ds, list(x ~ 1), c(0, 1), check_initialization = FALSE)
  expect_error(dp$f$f(1, dp$p, 0), "Function /\\(1, x\\(t\\) - 1\\) output non-finite value Inf")
})

test_that("units are parsed, propagated and validated", {
  tm <- setmetadata(independent_variable("tm"), "unit", "ms")
  tau <- Parameter("tau", unit = "ms")
  E <- Variable("E", tm, unit = "kJ"); P <- Variable("P", tm, unit = "MW")
  D <- Differential(tm)
  ue <- get_unit(E / tau)
  expect_equal(ue$scale, 1e6)
  expect_equal(unname(ue$dims), c(2, 1, -3, 0, 0, 0, 0))
  expect_true(.unit_equiv(get_unit(D(E)), .parse_unit("MW")))
  expect_true(validate(list(D(E) ~ P - E / tau, 0 ~ P)))
  E2 <- Variable("E2", tm, unit = "J")
  expect_warning(ok <- validate(D(E2) ~ P - E2 / tau), "do not match")
  expect_false(ok)
  expect_error(suppressWarnings(System(list(D(E2) ~ P - E2 / tau), tm, name = "bad")),
               class = "ValidationError")
  expect_s3_class(System(list(D(E2) ~ P - E2 / tau), tm, name = "unchecked", checks = CheckComponents), "System")
  expect_s3_class(System(list(D(E2) ~ P - E2 / tau), tm, name = "unchecked2", checks = CheckNone), "System")
  expect_error(suppressWarnings(System(list(D(E2) ~ P - E2 / tau), tm, name = "b2", checks = CheckUnits)),
               class = "ValidationError")
  expect_identical(c(CheckNone, CheckAll, CheckComponents, CheckUnits), c(0L, 1L, 2L, 4L))
  # parser: products, quotients, powers and prefixes
  expect_true(.unit_equiv(.parse_unit("kg*m^2/s^2"), .parse_unit("J")))
  expect_true(.unit_equiv(.parse_unit("1/s"), .parse_unit("Hz")))
  expect_true(.unit_equiv(.parse_unit("N*m"), .parse_unit("J")))
  expect_equal(.parse_unit("mmol/L")$scale, 1)
  expect_equal(.parse_unit("min")$scale, 60)
  expect_equal(.parse_unit("cm^(1/2)")$dims[["m"]], 0.5)
  # function rules: sqrt halves, exp needs dimensionless, abs keeps units
  L <- Parameter("L", unit = "m")
  expect_equal(get_unit(sqrt(L))$dims[["m"]], 0.5)
  expect_true(.unit_equiv(get_unit(abs(L)), .parse_unit("m")))
  expect_error(get_unit(exp(L)), class = "ValidationError")
  expect_error(get_unit(L + tau), class = "ValidationError")
  expect_true(.is_unitless(get_unit(L > 2 * L)))
  # noise terms carry the units of their equation
  ts <- setmetadata(independent_variable("ts"), "unit", "s")
  xs <- Variable("xs", ts, unit = "m")
  sg <- Parameter("sg", unit = "m/s"); sk <- Parameter("sk", unit = "kg")
  Ds <- Differential(ts)
  expect_s3_class(System(list(Ds(xs) ~ 0), ts, noise_eqs = list(sg), name = "n1"), "System")
  expect_error(suppressWarnings(System(list(Ds(xs) ~ 0), ts, noise_eqs = list(sk), name = "n2")),
               class = "ValidationError")
})

test_that("dynamic optimization respects bounds on observed variables", {
  t <- t_nounits; D <- D_nounits
  z <- Variable("z", t); power <- Variable("power", t, bounds = c(0, 1))
  w <- Variable("w", t, input = TRUE, bounds = c(-1, 1))
  lim <- mtkcompile(System(list(D(z) ~ w, power ~ 2 * z), t, costs = list(-EvalAt(1)(z)), name = "limited"),
                    inputs = list(w))
  prob <- DynamicOptProblem(lim, list(z ~ 0, w ~ 0), c(0, 1), dt = 0.1)
  expect_true(inherits(prob, AbstractDynamicOptProblem))
  sol <- solve(prob)
  expect_s3_class(sol, "DynamicOptSolution")
  expect_identical(sol$sol$retcode, "Default")
  expect_equal(tail(sol$sol[z], 1), 0.5, tolerance = 1e-5)
  expect_true(max(sol$sol[power]) <= 1 + 1e-5)
  expect_equal(sol$model$objective, -0.5, tolerance = 1e-5)
  expect_true(all(abs(sol$input_sol[, 1]) <= 1 + 1e-9))
  expect_error(DynamicOptProblem(lim, list(z ~ 0), c(0, 1)), "dt or the number")
  # every backend builds and solves the native problem with its own scheme
  for (be in list(list(JuMPDynamicOptProblem, JuMPCollocation, "RadauIIA5"),
                  list(CasADiDynamicOptProblem, CasADiCollocation, "RadauIIA5"),
                  list(InfiniteOptDynamicOptProblem, InfiniteOptCollocation, "Backward"),
                  list(PyomoDynamicOptProblem, PyomoCollocation, "LagrangeRadau(5)"))) {
    bp <- be[[1L]](lim, list(z ~ 0, w ~ 0), c(0, 1), dt = 0.1)
    expect_true(inherits(bp, AbstractDynamicOptProblem))
    bs <- solve(bp, be[[2L]]("ipopt"))
    expect_identical(bs$sol$retcode, "Default")
    expect_identical(bs$model$method, be[[3L]])
    expect_equal(tail(bs$sol[z], 1), 0.5, tolerance = 1e-5)
  }
  expect_true(inherits(JuMPCollocation("ipopt"), AbstractCollocation))
  expect_error(solve(JuMPDynamicOptProblem(lim, list(z ~ 0, w ~ 0), c(0, 1), dt = 0.1), CasADiCollocation()),
               "JuMPCollocation")
  expect_error(PyomoDynamicOptProblem(lim, list(z ~ 0, w ~ 0), c(0, 1), dt = 0.1,
                                      initial_trajectory = list(z ~ t)), "not supported by the Pyomo")
})

test_that("free final time: the double integrator reaches rest in time 2", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); v <- Variable("v", t); u <- Variable("u", t, bounds = c(-1, 1), input = TRUE)
  tf <- Parameter("tf")
  block <- System(list(D(x) ~ v, D(v) ~ u), t, costs = list(tf), name = "block",
                  constraints = list(EvalAt(tf)(v) ~ 0, EvalAt(tf)(x) ~ 0))
  block <- mtkcompile(block, inputs = list(u))
  prob <- suppressWarnings(DynamicOptProblem(block, list(x ~ 1, v ~ 0, u ~ 0, tf ~ 1), list(0, tf), steps = 11))
  sol <- solve(prob)
  expect_identical(sol$sol$retcode, "Default")
  expect_equal(sol$model$final_time, 2, tolerance = 2e-2)
  expect_equal(tail(sol$sol$t, 1), sol$model$final_time)
  expect_equal(tail(sol$sol[x], 1), 0, tolerance = 1e-5)
  expect_error(DynamicOptProblem(block, list(x ~ 1), list(0, tf), dt = 0.1), "number of steps")
  # the final time may appear in the dynamics; check the exact Jacobian
  y <- Variable("y", t)
  b2 <- mtkcompile(System(list(D(x) ~ v, D(v) ~ u, D(y) ~ x / tf + t), t, costs = list(tf), name = "b2",
                          constraints = list(EvalAt(tf)(x) ~ 0)), inputs = list(u))
  p2 <- suppressWarnings(DynamicOptProblem(b2, list(x ~ 1, v ~ 0, y ~ 0, u ~ 0, tf ~ 1.5), list(0, tf), steps = 4))
  ctx <- p2$ctx
  z0 <- ctx$z0 + seq_len(ctx$nz) / (10 * ctx$nz)
  ev <- .dopt_eval(ctx, z0)
  fd <- vapply(seq_len(ctx$nz), function(i) {
    e <- numeric(ctx$nz); e[i] <- 1e-6
    (.dopt_eval(ctx, z0 + e)$cE - .dopt_eval(ctx, z0 - e)$cE) / 2e-6
  }, numeric(length(ev$cE)))
  expect_equal(ev$JE, fd, tolerance = 1e-6)
})

test_that("algebraic states without initial values stay free", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); z <- Variable("z", t, guess = 1)
  dae <- mtkcompile(System(list(D(x) ~ z, 0 ~ z^2 - x - 1), t, name = "dae"))
  prob <- DynamicOptProblem(dae, list(x ~ 0), c(0, 1), dt = 0.1, guesses = list(z ~ 1))
  sol <- solve(prob)
  expect_identical(sol$sol$retcode, "Default")
  # x' = sqrt(x + 1), x(0) = 0 gives x(1) = (3/2)^2 - 1
  expect_equal(tail(sol$sol[x], 1), 1.25, tolerance = 5e-3)
})

test_that("parameter estimation recovers the Lotka Volterra parameters", {
  t <- t_nounits; D <- D_nounits
  al <- Parameter("al", default = 1.5); be <- Parameter("be", default = 1, tunable = FALSE)
  ga <- Parameter("ga", default = 3, tunable = FALSE); de <- Parameter("de", default = 1)
  xp <- Variable("xp", t); yp <- Variable("yp", t)
  eqs <- list(D(xp) ~ al * xp - be * xp * yp, D(yp) ~ -ga * yp + de * xp * yp)
  # data from the trapezoidal recursion itself, so the estimate is exact
  h <- 0.1; X <- matrix(c(4, 2), 1)
  rhs <- function(s) c(2.5 * s[1] - s[1] * s[2], -3 * s[2] + 1.8 * s[1] * s[2])
  for (i in 1:10) {
    s0 <- X[i, ]
    s1 <- .newton_solve(function(s) s - s0 - h / 2 * (rhs(s0) + rhs(s)), s0)$x
    X <- rbind(X, s1)
  }
  ts <- seq(0, 1, by = h)
  costs <- lapply(seq_along(ts), function(i) (EvalAt(ts[i])(xp) - X[i, 1])^2)
  spe <- mtkcompile(System(eqs, t, costs = costs, name = "spe"))
  prob <- DynamicOptProblem(spe, list(xp ~ 4, yp ~ 2), c(0, 1), dt = h, tune_parameters = TRUE)
  sol <- solve(prob)
  expect_equal(unname(sol$model$parameters), c(2.5, 1.8), tolerance = 1e-4)
  expect_equal(sol$sol$ps[al], 2.5, tolerance = 1e-4)
})

test_that("named, nonamespace and namespace mirror the macros", {
  t <- t_nounits; D <- D_nounits
  decay <- function(name, k = 1) {
    x <- Variable("x", t)
    System(list(D(x) ~ -k * x), t, name = name)
  }
  named(sys <- decay(k = 2))
  expect_identical(nameof(sys), "sys")
  named(other = decay())
  expect_identical(nameof(other), "other")
  kp <- Parameter("kp")
  named(comp <- decay(k = kp))
  expect_identical(attr(parameters(comp)[[1L]], "mtk_meta")$scope$type, "parent")
  expect_error(named(decay()), "named\\(sys <- ")
  calls <- 0L
  bump <- function() { calls <<- calls + 1L; 3 }
  named(once <- decay(k = bump()))
  expect_identical(calls, 1L)
  outer <- System(list(), t, name = "outer", systems = list(sys))
  expect_identical(format(outer$sys$x), "outer₊sys₊x(t)")
  expect_identical(format(nonamespace(outer$sys$x)), "x(t)")
  co <- complete(outer)
  expect_identical(format(namespace(co$sys$x)), "outer₊sys₊x(t)")
})

test_that("variable collection, substitution and conversion helpers", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t); k <- Parameter("k"); tau <- Parameter("tau")
  cv <- collect_vars(list(D(x) ~ -k * x + EvalAt(t - tau)(y)), t)
  expect_setequal(vapply(cv$unknowns, format, ""), c("x(t)", "y(t)"))
  expect_setequal(vapply(cv$parameters, format, ""), c("k", "tau"))
  expect_true(eqtype_supports_collect_vars(D(x) ~ x))
  expect_true(eqtype_supports_collect_vars(list(x, k)))
  expect_false(eqtype_supports_collect_vars(1))
  pk <- ParentScope(Parameter("pk"))
  inner <- System(list(D(x) ~ -pk * x), t, name = "inner")
  sc <- collect_scoped_vars(inner)
  expect_identical(vapply(sc$parameters, format, ""), "pk")
  expect_length(collect_vars(list(D(x) ~ -pk * x), t)$parameters, 0L)
  m <- collect_var_to_name(list(x, k))
  expect_identical(names(m), c("x", "k"))
  expect_error(collect_var_to_name(list(x, Variable("x"))), "duplicate names")
  gv <- get_variables(D(x) ~ -k * x + t)
  expect_identical(vapply(gv, format, ""), c("x(t)", "k", "t"))
  expect_identical(deparse(toexpr(D(x) ~ -k * sin(x))), "Differential(t)(x(t)) ~ -(sin(x(t)) * k)")
  sp <- independent_variable("sp"); fld <- Variable("fld", list(t, sp))
  expect_identical(deparse(toexpr(Differential(sp)(D(fld)))),
                   "Differential(sp)(Differential(t)(fld(t, sp)))")
  expect_identical(format(substitute_vars(-k * x, list(k ~ 2))), "-2*x(t)")
  expect_identical(format(apply_to_variables(function(v) if (identical(.name(v), "k")) tau else v, k * x)),
                   "tau*x(t)")
  sys <- System(list(D(x) ~ -k * x), t, name = "s")
  s2 <- substitute_vars(sys, list(k ~ 3))
  expect_identical(format(equations(s2)[[1L]]$rhs), "-3*x(t)")
  rep <- System(list(D(y) ~ -y), t, name = "s")
  top <- System(list(), t, name = "top", systems = list(sys))
  expect_identical(format(equations(substitute_vars(top, list(s = rep)))[[1L]]$lhs), "Derivative(s\u208ay(t), t)")
  expect_identical(format(Term("f", x, 2)), "f(x(t), 2)")
  expect_identical(format(Term(sin, x)), "sin(x(t))")
  a <- Parameter("a", default = 2)
  xb <- Variable("xb", t, default = 2 * a)
  cb <- convert_bindings_for_time_independent_system(System(list(D(xb) ~ -xb), t, name = "b"))
  expect_true(smap_has(cb$initial_conditions, xb))
  expect_false(smap_has(cb$bindings, xb))
  th <- Variable("th", t)
  R <- matrix(list(cos(th), sin(th), 0, -sin(th), cos(th), 0, 0, 0, 1), 3, 3)
  w <- get_w(R, t)
  expect_true(.is_zero(symcompR::simplify(w[[3L]] + D(th))))
})

test_that("metadata keys, verbosity and sentinels", {
  x <- Variable("x", t_nounits, unit = "m", bounds = c(0, 1), description = "pos")
  expect_identical(getmetadata(x, VariableUnit), "m")
  expect_identical(getmetadata(x, VariableBounds), c(0, 1))
  expect_identical(getmetadata(x, VariableDescription), "pos")
  expect_identical(getmetadata(x, MTKVariableTypeCtx), VariableType[["VARIABLE"]])
  expect_identical(SymScope, "scope")
  v <- MTKVerbosity(state_priority_tie = "Silent")
  expect_identical(v[["state_priority_tie"]], "Silent")
  expect_identical(v[["dynamic_opt_time_grid"]], "WarnLevel")
  expect_identical(unique(unclass(MTKVerbosity(preset = FALSE))), "Silent")
  expect_error(MTKVerbosity(nope = "Silent"), "unknown verbosity toggle")
  sys <- System(list(), t_nounits, name = "s", metadata = stats::setNames(list("black box"), SymbolicADDisallowed))
  expect_error(check_symbolic_ad_allowed(sys), "does not support symbolic automatic differentiation. black box")
  expect_true(check_symbolic_ad_allowed(System(list(), t_nounits, name = "ok")))
  expect_s3_class(Both, "Both")
  x2 <- Variable("x2", t_nounits, state_priority = 10, irreducible = TRUE)
  sp <- System(list(D_nounits(x2) ~ -x2), t_nounits, name = "sp")
  expect_equal(smap_get(state_priorities(sp), x2), 10)
  expect_identical(vapply(irreducibles(sp), format, ""), "x2(t)")
  expect_length(maybe_zeros(sp), 0L)
  expect_length(analytically_integrated(sp), 0L)
  expect_true(inherits(sp, AbstractSystem))
  expect_identical(c(ProblemTypeCtx, MiscSystemData), c("problem_type", "misc_data"))
})

test_that("collocation schemes, integral costs and the collocation polynomial", {
  t <- t_nounits; D <- D_nounits
  # min int_0^1 x^2 + u^2 with x' = u, x(0) = 1: u = -tanh(1 - t) x, the
  # cost tanh(1) and x(t) = cosh(1 - t) / cosh(1)
  x <- Variable("dq_x", t); u <- Variable("dq_u", t, input = TRUE)
  sys <- mtkcompile(System(list(D(x) ~ u), t, costs = list(Integral(x^2 + u^2, list(t, 0, 1))), name = "dqlqr"),
                    inputs = list(u))
  prob <- DynamicOptProblem(sys, list(x ~ 1, u ~ 0), c(0, 1), dt = 0.1)
  tq <- c(0.05, 0.55, 0.73)
  errs <- list(trapezoidal = 3e-3, hermite_simpson = 1e-6, RadauIIA5 = 1e-8, "LagrangeLegendre(2)" = 1e-6)
  for (sch in names(errs)) {
    s <- solve(prob, sch)
    expect_identical(s$sol$retcode, "Default")
    expect_equal(s$model$method, sch)
    expect_lt(abs(s$model$objective - tanh(1)), errs[[sch]])
    # between the nodes the states follow the collocation polynomial
    expect_lt(max(abs(s$sol(tq) - cosh(1 - tq) / cosh(1))), 10 * errs[[sch]])
    expect_equal(s$sol(0.6), s$sol$u[[7L]], tolerance = 1e-12)
  }
  s <- solve(prob, "hermite_simpson")
  expect_lt(max(abs(s$input_sol(tq) + sinh(1 - tq) / cosh(1))), 1e-3)
  # an explicit tableau given as a list is used as is; JuMP holds the input
  # over each step, which costs accuracy here
  rk4 <- list(A = rbind(c(0, 0, 0, 0), c(0.5, 0, 0, 0), c(0, 0.5, 0, 0), c(0, 0, 1, 0)), b = c(1, 2, 2, 1) / 6,
              c = c(0, 0.5, 0.5, 1))
  expect_lt(abs(solve(prob, rk4)$model$objective - tanh(1)), 1e-5)
  js <- solve(prob, JuMPCollocation("ipopt", rk4))
  expect_equal(js$input_sol(0.03), js$input_sol(0.07))
  # with an explicit tableau JuMP has no stage variables (one equality per
  # step), fixes x(0) and integrates the cost by the trapezoidal rule on the
  # nodes; the input at t = 0 carries half the weight of the others in the
  # cost, so the optimum pushes it hard. julia tools/ipopt-resto-probe.jl
  # (case lqr rk4_dt0.1, Ipopt 3.14.19): 21 free variables, 10 equalities,
  # one iteration; the states agree to 1e-10, the second derivatives being
  # differences of exact Jacobians here
  jp <- JuMPDynamicOptProblem(sys, list(x ~ 1, u ~ 0), c(0, 1), dt = 0.1)
  js <- solve(jp, JuMPCollocation("ipopt", rk4))
  expect_identical(js$model$iterations, 1L)
  expect_equal(js$model$objective, 0.71499799333216096, tolerance = 1e-12)
  expect_equal(unname(unlist(js$sol$u)),
               c(1, 0.8670004013335678, 0.8091706060136874, 0.759432516753944, 0.7172887526617397,
                 0.6823178760961531, 0.6541701782915279, 0.632564182269818, 0.6172838280708064,
                 0.6081763121525028, 0.605150559355724), tolerance = 1e-9)
  expect_identical(length(js$model$z), 22L)
  # implicit tableaus: JuMP adds the stage derivatives K (free at the last
  # node too), CasADi only on the steps; InfiniteOpt adds the derivative
  # variables at every support and integrates the cost by the trapezoidal
  # rule over all supports. julia tools/ipopt-resto-probe.jl (Ipopt
  # 3.14.19): the free variables (x(0) is fixed), equalities, iterations,
  # objective and the states at the nodes
  tabs <- list(list("backward_euler", 32L, 20L, 0.71499799333216107,
                    c(1, 0.867000401360674, 0.8091706060401394, 0.7594325167798005, 0.7172887526870726,
                      0.6823178761210441, 0.6541701783160673, 0.6325641822941027, 0.6172838280949389,
                      0.6081763121765895, 0.605150559379875)),
               list(list(A = rbind(c(0.25, -0.25), c(0.25, 5 / 12)), b = c(0.25, 0.75), c = c(0, 2 / 3)), 43L, 30L,
                    0.71499799333216107,
                    c(1, 0.8670004014195365, 0.8091706060964615, 0.7594325168338132, 0.7172887527390414,
                      0.6823178761712617, 0.6541701783648484, 0.6325641823417787, 0.6172838281418538,
                      0.6081763122230969, 0.6051505594263347)),
               list("RadauIIA5", 54L, 40L, 0.71499799333216085,
                    c(1, 0.8670004015310333, 0.8091706062025418, 0.7594325169349724, 0.717288752835844,
                      0.6823178762643262, 0.6541701784548348, 0.6325641824293795, 0.617283828227785,
                      0.608176312308092, 0.6051505595111382)))
  for (tb in tabs) {
    js <- solve(jp, JuMPCollocation("ipopt", tb[[1L]]))
    expect_identical(js$sol$retcode, "Default")
    expect_identical(js$model$iterations, 3L)
    expect_identical(length(js$model$z) - 1L, tb[[2L]])
    expect_identical(length(js$model$multipliers$equality), tb[[3L]])
    expect_equal(js$model$objective, tb[[4L]], tolerance = 1e-12)
    expect_equal(unname(unlist(js$sol$u)), tb[[5L]], tolerance = 1e-9)
  }
  cs <- solve(CasADiDynamicOptProblem(sys, list(x ~ 1, u ~ 0), c(0, 1), dt = 0.1), CasADiCollocation("ipopt"))
  expect_equal(cs$model$objective, 0.71499799333216085, tolerance = 1e-12)
  ip <- InfiniteOptDynamicOptProblem(sys, list(x ~ 1, u ~ 0), c(0, 1), dt = 0.1)
  is <- solve(ip, InfiniteOptCollocation("ipopt", "Backward"))
  expect_identical(is$model$iterations, 1L)
  expect_identical(length(is$model$z) - 1L, 32L)
  expect_identical(length(is$model$multipliers$equality), 21L)
  expect_equal(is$model$objective, 0.76226593057200254, tolerance = 1e-9)
  expect_equal(unname(unlist(is$sol$u)),
               c(1, 0.9287734069427998, 0.8668345479550275, 0.8135640344468056, 0.7684291612830516,
                 0.7309785797321282, 0.700837783978526, 0.6777053660647092, 0.6613500018115395,
                 0.6516081375764851, 0.645156571857906), tolerance = 1e-9)
  is <- solve(ip, InfiniteOptCollocation("ipopt", "OrthogonalCollocation(3)"))
  expect_identical(length(is$model$z) - 1L, 62L)
  expect_identical(length(is$model$multipliers$equality), 41L)
  expect_equal(is$model$objective, 0.6191490749537345, tolerance = 1e-9)
  expect_equal(unname(unlist(is$sol$u)),
               c(1, 0.8834502291613356, 0.7845971052785251, 0.7014604728475077, 0.6323749979483105,
                 0.5759568094598302, 0.5310757783373129, 0.4968328796707264, 0.4725421840568683,
                 0.45771711754918076, 0.45206000639967636), tolerance = 1e-9)
  expect_error(solve(prob, "nonsense"), "unknown collocation method")
  # Radau IIA with three stages is the default tableau of ModelingToolkit
  sc <- .dopt_scheme("RadauIIA5")
  expect_equal(sc$c, c(2 / 5 - sqrt(6) / 10, 2 / 5 + sqrt(6) / 10, 1), tolerance = 1e-12)
  expect_equal(sc$b, c(4 / 9 - sqrt(6) / 36, 4 / 9 + sqrt(6) / 36, 1 / 9), tolerance = 1e-12)
  expect_equal(sc$A[1L, ], c(11 / 45 - 7 * sqrt(6) / 360, 37 / 225 - 169 * sqrt(6) / 1800, -2 / 225 + sqrt(6) / 75),
               tolerance = 1e-12)
  expect_equal(.dopt_scheme("hermite_simpson")$b, c(1, 4, 1) / 6, tolerance = 1e-12)
})

test_that("free final time with an integral cost reaches tf = 10", {
  # the free final time test of ModelingToolkit: x' = -2 + u / 2, u in [0, 1],
  # max int_0^tf (x - u) dt + x(tf): u = 1 until tf - 1, tf = 10, cost -92.75
  t <- t_nounits; D <- D_nounits
  xf <- Variable("dq_xf", t); uf <- Variable("dq_uf", t, input = TRUE, bounds = c(0, 1)); tf <- Parameter("dq_tf")
  fsys <- mtkcompile(System(list(D(xf) ~ -2 + 0.5 * uf), t, name = "dqfft",
                            costs = list(-Integral(xf - uf, list(t, 0, tf)), -EvalAt(tf)(xf))), inputs = list(uf))
  fp <- suppressWarnings(JuMPDynamicOptProblem(fsys, list(xf ~ 17.5, uf ~ 0, tf ~ 8), list(0, tf), steps = 51))
  fs <- solve(fp, JuMPCollocation("ipopt", "hermite_simpson"))
  expect_identical(fs$sol$retcode, "Default")
  expect_equal(fs$model$final_time, 10, tolerance = 1e-3)
  # the tolerance of ModelingToolkit's test; the input is held over each step
  expect_lt(abs(fs$model$objective + 92.75), 0.25)
  hs <- solve(fp, "hermite_simpson")
  expect_equal(hs$model$final_time, 10, tolerance = 1e-3)
  expect_equal(hs$model$objective, -92.75, tolerance = 1e-5)
  expect_equal(tail(fs$sol$t, 1), fs$model$final_time)
  # exact first and second derivatives of the transcription
  ctx <- .dopt_grid(.subset2(DynamicOptProblem(fsys, list(xf ~ 17.5, uf ~ 0.5, tf ~ 8), list(0, tf), steps = 4),
                             "ctx"), .dopt_scheme("hermite_simpson"))
  z <- ctx$z0 + seq_len(ctx$nz) / (7 * ctx$nz)
  e1 <- .dopt_eval(ctx, z)
  nz <- ctx$nz
  fd <- function(g, h = 1e-6) sapply(seq_len(nz), function(i) {
    e <- numeric(nz); e[i] <- h
    (g(z + e) - g(z - e)) / (2 * h)
  })
  expect_equal(e1$JE, fd(function(zz) .dopt_eval(ctx, zz, 0L)$cE), tolerance = 1e-6)
  expect_equal(e1$gf, fd(function(zz) .dopt_eval(ctx, zz, 0L)$f), tolerance = 1e-6)
  lam <- seq_along(e1$cE) / 3
  e2 <- .dopt_eval(ctx, z, 2L, dense = FALSE, lamE = lam, obj = 0.5)
  H <- matrix(0, nz, nz)
  for (q in seq_along(e2$H$i)) H[e2$H$i[q], e2$H$j[q]] <- H[e2$H$i[q], e2$H$j[q]] + e2$H$v[q]
  gL <- function(zz) {
    e <- .dopt_eval(ctx, zz)
    0.5 * e$gf + as.numeric(crossprod(e$JE, lam))
  }
  expect_equal(H, fd(gL, 1e-5), tolerance = 1e-4)
})

test_that("the rocket of the dynamic optimization tutorial", {
  t <- t_nounits; D <- D_nounits
  h_c <- Parameter("rk_h_c"); m0 <- Parameter("rk_m0"); h0 <- Parameter("rk_h0"); g0 <- Parameter("rk_g0")
  D_c <- Parameter("rk_D_c"); cc <- Parameter("rk_c"); Tm <- Parameter("rk_Tm"); m_c <- Parameter("rk_m_c")
  h <- Variable("rk_h", t); v <- Variable("rk_v", t); m <- Variable("rk_m", t, bounds = list(m_c, 1))
  Th <- Variable("rk_T", t, input = TRUE, bounds = list(0, Tm))
  drag <- D_c * v^2 * exp(-h_c * (h - h0) / h0)
  eqs <- list(D(h) ~ v, D(v) ~ (Th - drag) / m - g0 * (h0 / h), D(m) ~ -Th / cc)
  te <- 0.2
  rocket <- System(eqs, t, list(h, v, m, Th), list(h_c, m0, h0, g0, D_c, cc, Tm, m_c), name = "rocket",
                   costs = list(-EvalAt(te)(h)), constraints = list(EvalAt(te)(Th) ~ 0, EvalAt(te)(m) ~ m_c))
  rocket <- mtkcompile(rocket, inputs = list(Th))
  op <- list(h ~ 1, m ~ 1, v ~ 0, g0 ~ 1, m0 ~ 1, h_c ~ 500, cc ~ 0.5, D_c ~ 0.5 * 620, Tm ~ 3.5, Th ~ 0, h0 ~ 1,
             m_c ~ 0.6)
  # the tutorial uses dt = 0.001; ModelingToolkit's test asks for h(te) > 1.012
  jprob <- suppressWarnings(JuMPDynamicOptProblem(rocket, op, c(0, te), dt = 0.005))
  jsol <- solve(jprob, JuMPCollocation("ipopt"))
  expect_identical(jsol$sol$retcode, "Default")
  expect_identical(jsol$model$method, "RadauIIA5")
  expect_gt(tail(jsol$sol[h], 1), 1.012)
  expect_equal(tail(jsol$sol[m], 1), 0.6, tolerance = 1e-6)
  expect_lt(abs(tail(jsol$input_sol$u, 1)[[1L]]), 1e-4)
  expect_true(all(unlist(jsol$input_sol$u) >= -1e-8 & unlist(jsol$input_sol$u) <= 3.5 + 1e-8))
  # integrating the dynamics with the optimal thrust reproduces the trajectory
  f <- function(y, p, tt) {
    thr <- jsol$input_sol(tt)[[1L]]
    hh <- y[3L]; vv <- y[1L]; mm <- y[2L]
    c((thr - 310 * vv^2 * exp(-500 * (hh - 1))) / mm - 1 / hh, -thr / 0.5, vv)
  }
  osol <- solve(ODEProblem(f, c(0, 1, 1), c(0, te)), Tsit5(), abstol = 1e-10, reltol = 1e-10,
                tstops = jsol$sol$t)
  yend <- tail(osol$u, 1)[[1L]]
  expect_equal(yend[3L], tail(jsol$sol[h], 1), tolerance = 2e-3)
  expect_equal(yend[2L], 0.6, tolerance = 2e-3)
})

test_that("the cart pole swings up", {
  # the cart pole test of ModelingToolkit: reach the upright rest position
  # at tf = 5 with the least integral of u^2; the swing up needs a starting
  # trajectory (initial_trajectory of the tutorial) for the local solver
  t <- t_nounits; D <- D_nounits
  g <- Parameter("cp_g"); l <- Parameter("cp_l"); mp <- Parameter("cp_mp"); mk <- Parameter("cp_mk")
  x <- Variable("cp_x", t); th <- Variable("cp_th", t); u <- Variable("cp_u", t, input = TRUE, bounds = c(-10, 10))
  s <- sin(th); c <- cos(th)
  # qdd = -H \ (C qd + G - B u) for H = [mk + mp, mp l c; mp l c, mp l^2]
  v1 <- -mp * l * s * D(th)^2 - u; v2 <- mp * g * l * s
  det <- (mk + mp) * mp * l^2 - mp^2 * l^2 * c^2
  tf <- 5
  cons <- list(EvalAt(tf)(th) ~ pi, EvalAt(tf)(x) ~ 0, EvalAt(tf)(D(th)) ~ 0, EvalAt(tf)(D(x)) ~ 0)
  cartpole <- System(list(D(D(x)) ~ -(mp * l^2 * v1 - mp * l * c * v2) / det,
                          D(D(th)) ~ -(-mp * l * c * v1 + (mk + mp) * v2) / det),
                     t, list(x, th, u), list(g, l, mp, mk), name = "cartpole",
                     costs = list(Integral(u^2, list(t, 0, tf))), constraints = cons)
  cartpole <- mtkcompile(cartpole, inputs = list(u))
  op <- list(D(x) ~ 0, D(th) ~ 0, th ~ 0, x ~ 0, mk ~ 1, mp ~ 0.2, l ~ 0.5, g ~ 9.81, u ~ 0)
  prob <- suppressWarnings(JuMPDynamicOptProblem(cartpole, op, c(0, tf), dt = 0.25,
                                                 initial_trajectory = list(th ~ pi * t / 5)))
  sol <- solve(prob, JuMPCollocation("ipopt", "hermite_simpson"))
  expect_identical(sol$sol$retcode, "Default")
  # the states are x, theta and their lowered derivatives
  final <- tail(sol$sol$u, 1)[[1L]]
  expect_equal(unname(final[c("cp_th(t)", "cp_x(t)")]), c(pi, 0), tolerance = 1e-7)
  expect_equal(max(abs(final[setdiff(names(final), c("cp_th(t)", "cp_x(t)"))])), 0, tolerance = 1e-7)
  expect_equal(sol$sol(tf)[["cp_th(t)"]], pi, tolerance = 1e-7)
  expect_true(all(abs(unlist(sol$input_sol$u)) <= 10 + 1e-8))
  # from rest with the RK4 tableau of ModelingToolkit's test, which needs no
  # starting trajectory: julia tools/ipopt-resto-probe.jl (Ipopt 3.14.19,
  # Linux) reports 101 variables (the 4 fixed initial states removed), 84
  # equality constraints and Optimal Solution Found after 60 iterations
  skip_on_cran()
  prob <- suppressWarnings(JuMPDynamicOptProblem(cartpole, op, c(0, tf), dt = 0.25))
  sol <- solve(prob, JuMPCollocation("ipopt", "RK4"))
  expect_identical(sol$sol$retcode, "Default")
  expect_identical(length(sol$model$z), 105L)
  expect_identical(length(sol$model$multipliers$equality), 84L)
  expect_identical(sol$model$iterations, 60L)
  expect_equal(sol$model$objective, 27.398119, tolerance = 1e-7)
  final <- tail(sol$sol$u, 1)[[1L]]
  expect_equal(unname(final[c("cp_th(t)", "cp_x(t)")]), c(pi, 0), tolerance = 1e-8)
})

test_that("collocation derivatives are exact for every kind of term", {
  # path constraints with point values and integrals, integral costs,
  # inputs held over a step (JuMP) and Gauss points with tied inputs
  t <- t_nounits; D <- D_nounits
  x <- Variable("dd_x", t); y <- Variable("dd_y", t); u <- Variable("dd_u", t, input = TRUE)
  sys <- mtkcompile(System(list(D(x) ~ u - x * y, D(y) ~ x * u), t, name = "dddiff",
                           costs = list(Integral(x^2 * u + sin(y), list(t, 0, 1)), EvalAt(1)(y)^2),
                           constraints = list(Inequality(x * EvalAt(1)(x) + Integral(u^2, list(t, 0, 1)) * y, 2),
                                              EvalAt(0.5)(x) * EvalAt(1)(y) ~ 0.3)),
                    inputs = list(u))
  prob <- suppressWarnings(DynamicOptProblem(sys, list(x ~ 1, y ~ 0, u ~ 0.2), c(0, 1), dt = 0.5))
  for (sch in list(.dopt_scheme("hermite_simpson"), JuMPCollocation()$scheme, .dopt_scheme("LagrangeLegendre(2)"),
                   CasADiCollocation()$scheme, InfiniteOptCollocation("ipopt", "OrthogonalCollocation(3)")$scheme,
                   InfiniteOptCollocation("ipopt", "Backward")$scheme)) {
    ctx <- .dopt_grid(.subset2(prob, "ctx"), sch)
    nz <- ctx$nz
    z <- ctx$z0 + sin(seq_len(nz)) / 5
    e1 <- .dopt_eval(ctx, z)
    fd <- function(g, h = 1e-6) sapply(seq_len(nz), function(i) {
      e <- numeric(nz); e[i] <- h
      (g(z + e) - g(z - e)) / (2 * h)
    })
    expect_equal(e1$JE, fd(function(zz) .dopt_eval(ctx, zz, 0L)$cE), tolerance = 1e-6, info = sch$name)
    expect_equal(e1$JI, fd(function(zz) .dopt_eval(ctx, zz, 0L)$cI), tolerance = 1e-6, info = sch$name)
    expect_equal(e1$gf, fd(function(zz) .dopt_eval(ctx, zz, 0L)$f), tolerance = 1e-6, info = sch$name)
    lE <- cos(seq_along(e1$cE)); lI <- seq_along(e1$cI) / 4
    e2 <- .dopt_eval(ctx, z, 2L, dense = FALSE, lamE = lE, lamI = lI, obj = 0.7)
    H <- matrix(0, nz, nz)
    for (q in seq_along(e2$H$i)) H[e2$H$i[q], e2$H$j[q]] <- H[e2$H$i[q], e2$H$j[q]] + e2$H$v[q]
    gL <- function(zz) {
      e <- .dopt_eval(ctx, zz)
      0.7 * e$gf + as.numeric(crossprod(e$JE, lE)) + as.numeric(crossprod(e$JI, lI))
    }
    expect_equal(H, fd(gL, 1e-5), tolerance = 1e-4, info = sch$name)
  }
})
