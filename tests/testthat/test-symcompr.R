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
  expect_true(all(shared %in% c("flatten", "compose", "Differential", "extend", "System", "Variable", "WienerProcess", "Stream", "distance", "hadamard_product", "IdentityOperator")))
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

# every name both packages export: a call that is not for modelcompR reaches
# symcompR exactly as written, named arguments and all
test_that("shared names forward symcompR calls with their original arguments", {
  sc <- function(name) getExportedValue("symcompR", name)
  has <- function(name) name %in% getNamespaceExports("symcompR")
  x <- symcompR::Symbol("x"); y <- symcompR::Symbol("y")
  nested <- list(1, list(2, list(3)))
  expect_identical(flatten(x = nested), symcompR::flatten(x = nested))
  expect_identical(flatten(nested), symcompR::flatten(nested))
  expect_equal(format(compose(f = x^2, g = x + 1)), format(symcompR::compose(f = x^2, g = x + 1)))
  expect_equal(format(compose(g = x + 1, f = x^2)), "x^2 + 2*x + 1")
  expect_equal(format(compose(x^2, y + 1, x = x)), format(symcompR::compose(x^2, y + 1, x = x)))
  k <- symcompR::Symbol("k")
  A <- symcompR::MatrixSymbol("A", k, k); B <- symcompR::MatrixSymbol("B", k, k)
  expect_equal(format(hadamard_product(A, B)), format(symcompR::hadamard_product(A, B)))
  # numeric matrices stay modelcompR's g-convex atom
  expect_equal(hadamard_product(matrix(1:4, 2), diag(2)), diag(c(1, 4)))
  if (length(formals(symcompR::hadamard_product)) == 1L)
    expect_equal(format(hadamard_product(B, A, A)), format(symcompR::hadamard_product(B, A, A)))
  p <- symcompR::Point(c(0, 0)); q <- symcompR::Point(c(3, 4))
  expect_equal(format(distance(p, q)), "5")
  expect_equal(format(distance(p = p, q = q)), "5")
  expect_equal(format(distance(p, q = q)), "5")
  expect_equal(format(distance(q = q, p = p)), "5")
  expect_s3_class(WienerProcess(0, 0), "WienerProcess")
  expect_s3_class(Differential(t_nounits), "Differential")
  if (has("WienerProcess")) {
    expect_equal(format(WienerProcess(sym = "W")), format(sc("WienerProcess")(sym = "W")))
    expect_equal(format(WienerProcess("W")), format(sc("WienerProcess")("W")))
  }
  if (has("Differential")) {
    rn <- getExportedValue("symcompR", "diffgeom_rn")()
    f <- getExportedValue("symcompR", "base_scalar")(rn$R2_r, 0L)
    expect_equal(format(Differential(form_field = f)), format(sc("Differential")(form_field = f)))
    expect_equal(format(Differential(f^2)), format(sc("Differential")(f^2)))
  }
  # a dependent variable is a modelcompR operator, whatever symcompR exports
  xt <- Variable("xdg", t_nounits)
  Dx <- Differential(xt)
  expect_s3_class(Dx, "Differential")
  expect_equal(format(Dx(xt^2)), format(Derivative(xt^2, xt, 1L)))
  if (has("extend")) {
    ds <- sc("DimensionSystem")(list(symcompR::Dimension("length")))
    e1 <- extend(x = ds, new_base_dims = list(symcompR::Dimension("time")))
    e2 <- sc("extend")(x = ds, new_base_dims = list(symcompR::Dimension("time")))
    expect_s3_class(e1, "DimensionSystem")
    expect_equal(e1, e2)
  }
  if (has("Variable")) {
    expect_equal(Variable(symbol = x), sc("Variable")(symbol = x))
    expect_equal(Variable(x), sc("Variable")(x))
    real <- getExportedValue("symcompR", "real")
    expect_equal(Variable("x", real), sc("Variable")("x", real))
    expect_equal(getdefault(Variable("vsc", t_nounits, default = 2)), 2)
  }
  if (has("Stream")) expect_equal(format(Stream("stdout")), format(sc("Stream")("stdout")))
  expect_identical(Stream, "Stream")
  if (has("System")) {
    expect_s3_class(System(), "physics_obj")
    expect_equal(class(System(frame = NULL)), class(sc("System")(frame = NULL)))
  }
})

test_that("IdentityOperator without a numeric size is symcompR's quantum operator", {
  skip_if_not("IdentityOperator" %in% getNamespaceExports("symcompR"))
  expect_identical(format(IdentityOperator()), format(symcompR::IdentityOperator()))
  N <- symcompR::Symbol("N")
  expect_identical(format(IdentityOperator(N)), format(symcompR::IdentityOperator(N)))
  # a numeric size is the SciMLOperators identity
  expect_equal(as.matrix(IdentityOperator(2)), diag(2))
})
