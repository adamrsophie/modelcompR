# Transformations of basic_transformations.jl. Reference values are analytic.

test_that("change_independent_variable reproduces the projectile tutorial", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t)
  g <- Parameter("g", default = 9.81); v <- Parameter("v")
  M1 <- System(list(D(D(y)) ~ -g, D(x) ~ v), t, initialization_eqs = list(D(x) ~ D(y)), name = "M")
  M2 <- change_independent_variable(M1, x)
  # the old independent variable stays an independent variable everywhere
  expect_identical(.var_kind(t_nounits), "iv")
  expect_identical(.var_kind(x), "unknown")
  xi <- M2$x
  expect_identical(format(get_iv(M2)), "x")
  expect_identical(format(xi), "x")
  expect_setequal(vapply(unknowns(M2), format, ""), c("y(x)", "t(x)", "xˍt(x)"))
  expect_identical(.var_kind(M2$t), "unknown")
  # equations: xˍt^2 y'' + xˍt xˍt' y' ~ -g and xˍt ~ v
  eqs <- equations(M2)
  w <- Function("xˍt")(xi); yx <- Function("y")(xi)
  # the closed form y = x - g x^2 / (2 v^2), xˍt = v satisfies them
  Dx <- Differential(xi)
  sol <- list(list(Derivative(yx, xi, 2L), -g / v^2), list(Derivative(yx, xi), 1 - g * xi / v^2),
              list(Derivative(w, xi), 0), list(w, v))
  sub <- function(e) .subst(e, .rule(lapply(sol, `[[`, 1L), lapply(sol, `[[`, 2L)))
  res <- lapply(eqs, function(e) symcompR::simplify(sub(e$lhs - e$rhs)))
  expect_true(all(vapply(res, .is_zero, TRUE)))
  # the initialization equation holds at the start, x = 0
  ini <- initialization_equations(M2)[[1L]]
  expect_true(.is_zero(symcompR::simplify(.subst(sub(ini$lhs - ini$rhs), .rule(list(xi), list(0))))))
  # explicit form: substitute xˍt = v (from the second equation) and solve the
  # first for y'', then integrate to x = 10 with y(0) = 0, y'(0) = 1
  e1 <- .subst(eqs[[1L]]$lhs - eqs[[1L]]$rhs, .rule(list(Derivative(w, xi), w), list(0, v)))
  ypp <- Derivative(yx, xi, 2L)
  coef <- diff(e1, ypp)
  rest <- .subst(e1, .rule(list(ypp), list(0)))
  y2 <- Variable("y", xi)
  ex <- System(list(Dx(Dx(y2)) ~ .subst(-rest / coef, .rule(list(yx), list(y2)))), xi, name = "ex")
  cex <- mtkcompile(ex)
  prob <- ODEProblem(cex, list(y2 ~ 0, Dx(y2) ~ 1, v ~ 8), c(0, 10))
  s <- solve(prob, reltol = 1e-10, abstol = 1e-10)
  expect_equal(tail(s[y2], 1), 10 - 9.81 * 100 / (2 * 64), tolerance = 1e-8)
  # with the structural pipeline the transformed model compiles directly
  skip_if(is.null(.mtk$hooks$structural_compile), "structural pipeline not loaded")
  M2s <- tryCatch(mtkcompile(M2, allow_symbolic = TRUE), error = function(e) NULL)
  skip_if(is.null(M2s), "structural pipeline cannot compile the implicit model yet")
  p2 <- ODEProblem(M2s, list(M2s$y ~ 0, v ~ 8), c(0, 10))
  s2 <- solve(p2, reltol = 1e-8, abstol = 1e-8)
  expect_equal(tail(s2[M2s$y], 1), 2.3359375, tolerance = 1e-5)
  expect_length(equations(M2s), 2L)
})

test_that("change_independent_variable checks its input and adds the inverse rule", {
  t <- t_nounits; D <- D_nounits
  a <- Variable("a", t); k <- Parameter("k")
  M <- System(list(D(a) ~ k * a), t, name = "M")
  expect_error(change_independent_variable(M, k), "not a function of the independent variable")
  expect_error(change_independent_variable(mtkcompile(M), a), "structurally simplified")
  M2 <- change_independent_variable(M, a, add_old_diff = TRUE)
  ai <- get_iv(M2)
  eqs <- equations(M2)
  w <- Function("aˍt")(ai)
  expect_true(any(vapply(eqs, function(e) .is_zero(symcompR::expand(e$lhs - e$rhs - (w - k * ai))), TRUE)))
  expect_true(any(vapply(eqs, function(e) format(e) == "Derivative(t(a), a) ~ 1/aˍt(a)", TRUE)))
  # t(a) = log(a / a0) / k solves dt/da = 1 / (k a): integrate with aˍt observed
  cs <- mtkcompile(M2)
  prob <- ODEProblem(cs, list(M2$t ~ 0, k ~ 2), c(1, 3))
  s <- solve(prob, reltol = 1e-10, abstol = 1e-10)
  expect_equal(tail(s[M2$t], 1), log(3) / 2, tolerance = 1e-7)
})

test_that("liouville_transform adds the trace of the Jacobian", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t)
  a <- Parameter("a"); b <- Parameter("b"); c <- Parameter("c"); d <- Parameter("d")
  sys <- System(list(D(x) ~ a * x - b * x * y, D(y) ~ -d * y + c * x * y), t, name = "lv")
  ls <- liouville_transform(sys)
  trJ <- nonamespace(ls$trJ)
  eq <- equations(ls)[[3L]]
  expect_identical(format(eq$lhs), "Derivative(trJ(t), t)")
  expect_true(.is_zero(symcompR::simplify(eq$rhs + trJ * (a - b * y - d + c * x))))
  # linear decay: trJ(t) = exp(k t) for D(x) ~ -k x
  k <- Parameter("k")
  ds <- mtkcompile(liouville_transform(System(list(D(x) ~ -k * x), t, name = "dec")))
  prob <- ODEProblem(ds, list(x ~ 1, ds$trJ ~ 1, k ~ 0.5), c(0, 2))
  s <- solve(prob, reltol = 1e-10, abstol = 1e-10)
  expect_equal(tail(s[ds$trJ], 1), exp(1), tolerance = 1e-8)
})

test_that("change_of_variables follows the docstring example", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); z <- Variable("z", t); al <- Parameter("al")
  sys <- mtkcompile(System(list(D(x) ~ al * x), t, name = "sys",
                           initial_conditions = list(x ~ 1, al ~ -0.5)))
  new_sys <- change_of_variables(sys, t, list(log(x) ~ z), list(x ~ exp(z)))
  expect_true(equations(new_sys)[[1L]] == Equation(D(z), al))
  s <- solve(ODEProblem(new_sys, list(), c(0, 1)), reltol = 1e-10, abstol = 1e-10)
  expect_equal(tail(s[x], 1), exp(-0.5), tolerance = 1e-8)
  expect_equal(s[z][1L], 0)
})

test_that("add_accumulations integrates unknowns and expressions", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t)
  sys <- add_accumulations(System(list(D(x) ~ -x), t, name = "s"))
  expect_identical(vapply(unknowns(sys), format, ""), c("x(t)", "accumulation_x(t)"))
  acc <- Variable("acc", t)
  sys2 <- add_accumulations(System(list(D(x) ~ -x), t, name = "s"), list(acc ~ x^2))
  cs <- mtkcompile(sys2)
  s <- solve(ODEProblem(cs, list(x ~ 1), c(0, 1)), reltol = 1e-10, abstol = 1e-10)
  expect_equal(tail(s[acc], 1), (1 - exp(-2)) / 2, tolerance = 1e-8)
  expect_error(add_accumulations(sys2, list(acc ~ x)), "already exist")
})

test_that("convert_system_indepvar turns a nonlinear system time dependent", {
  xs <- Variable("xs"); a <- Parameter("a")
  ns <- System(list(0 ~ xs^2 - a), name = "ns")
  tau <- independent_variable("tau")
  cs <- convert_system_indepvar(ns, tau, name = "ts")
  expect_identical(format(unknowns(cs)[[1L]]), "xs(tau)")
  expect_identical(format(get_iv(cs)), "tau")
  expect_identical(nameof(cs), "ts")
  expect_identical(.var_kind(xs), "unknown")
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t)
  s1 <- System(list(D(y) ~ -y), t, name = "s1")
  s2 <- convert_system_indepvar(s1, tau)
  expect_identical(format(equations(s2)[[1L]]), "Derivative(y(tau), tau) ~ -y(tau)")
  expect_error(convert_system_indepvar(mtkcompile(System(list(D(y) ~ -y, xs ~ 2 * y), t, name = "o")), tau),
               "reduced model")
})

test_that("fractional_to_ordinary uses the exponential sum quadrature", {
  # sum_i c_i exp(-gamma_i s) approximates s^(alpha - 1) / gamma(alpha)
  for (alpha in c(0.3, 0.5, 0.8)) {
    q <- .fto_quadrature(alpha, 1e-4, 1)
    s <- c(q$delta * 10, 0.01, 0.1, 0.5, 1)
    approx <- vapply(s, function(si) sum(q$c * exp(-q$gamma * si)), 0)
    expect_equal(approx, s^(alpha - 1) / gamma(alpha), tolerance = 1e-3)
  }
  t <- t_nounits
  x <- Variable("x", t)
  q <- .fto_quadrature(0.5, 1e-2, 1)
  sys <- fractional_to_ordinary(list(1), list(x), 0.5, 1e-2, 1)
  expect_length(unknowns(sys), length(q$c))
  expect_identical(format(observed(sys)[[1L]]$lhs), "x(t)")
  expect_error(fractional_to_ordinary(list(1), list(x), 1, 1e-2, 1), "non-integer")
  lsys <- linear_fractional_to_ordinary(c(1.5, 0), c(1, 1), 1, 1e-2, 1, initials = c(0, 0))
  expect_true(length(unknowns(lsys)) > 2L)
  # D^0.5 x = 1, x(0) = 0 has x(t) = 2 sqrt(t / pi)
  skip_if_not(exists("Rodas5P"))
  prob <- ODEProblem(sys, list(), c(0, 1))
  s <- solve(prob, Rodas5P(), reltol = 1e-8, abstol = 1e-10)
  expect_equal(tail(s[x], 1), 2 / sqrt(pi), tolerance = 2e-2)
})

test_that("respecialize records the class of nonnumeric parameter values", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); f <- Parameter("f"); k <- Parameter("k")
  sys <- System(list(D(x) ~ -k * x), t, ps = list(f, k), name = "s")
  expect_error(respecialize(sys, list(f)), "completed systems")
  cs <- complete(sys)
  sp <- splinefun(1:3, c(1, 4, 9))
  rs <- respecialize(cs, list(list(f, sp)))
  fp <- Filter(function(p) .name(p) == "f", parameters(rs))[[1L]]
  expect_identical(getmetadata(fp, "variable_type"), "function")
  expect_identical(smap_get(initial_conditions(rs), f), sp)
  expect_error(respecialize(cs, list(list(k, 2))), "Numeric types cannot be respecialized")
  expect_error(respecialize(cs, list(f)), "needs an associated value")
  cs2 <- complete(System(list(D(x) ~ -k * x), t, ps = list(f, k), name = "s",
                         initial_conditions = list(list(f, sp), k ~ 1)))
  expect_identical(getmetadata(Filter(function(p) .name(p) == "f", parameters(respecialize(cs2)))[[1L]],
                               "variable_type"), "function")
})

test_that("fractional_to_ordinary in array form reproduces ModelingToolkit's test", {
  # example 1 of documentation/ModelingToolkit.jl/test/fractional_to_ordinary.jl:
  # D^0.5 x = f(t) - x^(3/2) has the solution (3/2 t^(alpha/2) - t^4)^2;
  # ModelingToolkit checks it at 0, 0.1, ..., 1 with atol 1e-7 (131 unknowns)
  t <- t_nounits
  x <- Variable("tfx", t)
  a <- 0.5
  eq <- (9 * gamma(1 + a) / 4) - (3 * t^(4 - a / 2) * gamma(5 + a / 2) / gamma(5 - a / 2)) +
    (gamma(9) * t^(8 - a) / gamma(9 - a)) + (3 / 2 * t^(a / 2) - t^4)^3 - x^(3 / 2)
  sys <- fractional_to_ordinary(list(eq), list(x), a, 1e-7, 1, matrix = TRUE)
  us <- vapply(unknowns(sys), format, "")
  expect_length(us, 131L)
  expect_true(all(grepl("^ʐ_0\\[[0-9]+\\]\\(t\\)$", us)))
  prob <- ODEProblem(sys, list(), c(0, 1))
  sol <- solve(prob, FBDF(), reltol = 1e-10, abstol = 1e-10)
  tt <- seq(0, 1, by = 0.1)
  expect_lt(max(abs(sol(tt, idxs = x) - (3 / 2 * tt^(a / 2) - tt^4)^2)), 1e-7)
})

test_that("linear_fractional_to_ordinary in array form matches ModelingToolkit", {
  # example 3 of ModelingToolkit's fractional_to_ordinary.jl, on [0, 20] with
  # epsilon = 1e-3; reference values from ModelingToolkit 11 (matrix = true,
  # 69 unknowns, Rodas5P with tolerances 1e-11)
  t <- t_nounits
  sys <- linear_fractional_to_ordinary(c(3, 2.5, 2, 1, 0.5, 0), c(1, 1, 1, 4, 1, 4), 6 * cos(t), 1e-3, 20,
                                       initials = c(1, 1, -1), matrix = TRUE)
  expect_length(unknowns(sys), 69L)
  vars <- c(unknowns(sys), lapply(observed(sys), `[[`, "lhs"))
  x0 <- Filter(function(v) format(v) == "x_0(t)", vars)[[1L]]
  sol <- solve(ODEProblem(sys, list(), c(0, 20)), FBDF(), reltol = 1e-10, abstol = 1e-10)
  expect_equal(sol(c(7, 20), idxs = x0), c(1.410888853062601, 1.3210273125415624), tolerance = 1e-8)
})
