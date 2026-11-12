# Compatibility with symcompR: modelcompR expressions are symcompR
# expressions, symcompR objects work in models, and attaching both packages
# changes nothing in either.

test_that("model expressions work with symcompR functions", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t); k <- Parameter("k")
  sys <- mtkcompile(System(list(D(x) ~ -k * x + y, y ~ 2 * x), t, name = "s"))
  # full_equations() substitutes the observed y ~ 2x
  rhs <- full_equations(sys)[[1]]$rhs
  expect_true(is_sym(rhs))
  expect_equal(format(symcompR::simplify(rhs - (2 - k) * x)), "0")
  expect_match(symcompR::latex(rhs), "x")
  expect_equal(format(diff(rhs, x)), format(2 - k))
  J <- calculate_jacobian(sys)
  expect_s3_class(J, "symmatrix")
  expect_equal(format(symcompR::sym_det(J)), format(2 - k))
  f <- symcompR::lambdify(list(k), symcompR::subs(rhs, x, 1))
  expect_equal(f(3), -1)
})

test_that("plain symcompR symbols and functions take part in models", {
  # time independent: undeclared symbols are unknowns, as in ModelingToolkit
  a <- symcompR::Symbol("aplain"); b <- symcompR::Symbol("bplain")
  ns <- System(list(0 ~ a^2 - 2, 0 ~ b - a), name = "plain")
  expect_setequal(vapply(unknowns(ns), format, ""), c("aplain", "bplain"))
  # undeclared applications f(t) are unknowns of a time dependent system
  t <- t_nounits; D <- D_nounits
  z <- symcompR::Function("zfun")(t)
  k <- Parameter("kz")
  sys <- System(list(D(z) ~ -k * z), t, name = "plainfun")
  expect_equal(vapply(unknowns(sys), format, ""), "zfun(t)")
  prob <- ODEProblem(mtkcompile(sys), list(z ~ 1, k ~ 2), c(0, 1))
  sol <- solve(prob, Tsit5(), reltol = 1e-9, abstol = 1e-12)
  expect_equal(tail(sol[z], 1), exp(-2), tolerance = 1e-7)
  # mathematical functions are never unknowns
  e <- sin(z) + symcompR::erf(z)
  found <- Filter(function(v) !identical(.var_kind(v), "iv"), .collect_vars(e))
  expect_equal(vapply(found, format, ""), "zfun(t)")
  # an undeclared symbol without time dependence is rejected, as check_variables does
  q <- symcompR::Symbol("qundeclared")
  expect_error(System(list(D(z) ~ -q * z), t, name = "bad"), "not a function of the independent variable")
})

test_that("solve() keeps working for matrices and symcompR equations", {
  expect_equal(solve(matrix(c(2, 0, 0, 4), 2), c(2, 4)), c(1, 1))
  x <- symcompR::Symbol("x")
  r <- solve(x^2 - 4, x)
  expect_setequal(vapply(r, format, ""), c("-2", "2"))
})

test_that("attaching modelcompR masks nothing from symcompR or base R", {
  ns <- getNamespaceExports("modelcompR")
  shared <- intersect(ns, getNamespaceExports("symcompR"))
  expect_true(all(c("flatten", "compose") %in% shared))
  expect_true(all(shared %in% c("flatten", "compose", "Differential")))
  x <- symcompR::Symbol("x")
  expect_equal(format(compose(x^2, x + 1)), format(symcompR::compose(x^2, x + 1)))
})

test_that("undeclared functions named like R or symcompR functions are unknowns", {
  t <- t_nounits; D <- D_nounits
  S <- symcompR::Function("S")(t); I <- symcompR::Function("I")(t); R <- symcompR::Function("R")(t)
  b <- Parameter("bsir"); g <- Parameter("gsir")
  sys <- System(list(D(S) ~ -b * S * I, D(I) ~ b * S * I - g * I, D(R) ~ g * I), t, name = "sir")
  expect_setequal(vapply(unknowns(sys), format, ""), c("S(t)", "I(t)", "R(t)"))
  prob <- ODEProblem(mtkcompile(sys), list(S ~ 0.99, I ~ 0.01, R ~ 0, b ~ 0.5, g ~ 0.25), c(0, 50))
  sol <- solve(prob)
  N <- sol[S] + sol[I] + sol[R]
  expect_equal(N, rep(1, length(N)), tolerance = 1e-6)
  # mathematical heads stay functions
  expect_null(.var_kind(symcompR::erf(S)))
  expect_null(.var_kind(sin(S)))
})

test_that("Differential delegates differential geometry to symcompR", {
  skip_if_not("Differential" %in% getNamespaceExports("symcompR"))
  rn <- getExportedValue("symcompR", "diffgeom_rn")()
  x <- getExportedValue("symcompR", "base_scalar")(rn$R2_r, 0L)
  expect_equal(format(Differential(x)), format(getExportedValue("symcompR", "Differential")(x)))
  expect_s3_class(Differential(t_nounits), "Differential")
})
