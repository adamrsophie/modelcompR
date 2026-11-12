# Cross-module behavior found in the integration review.

test_that("Pre() works on plain symcompR function variables in events", {
  t <- t_nounits; D <- D_nounits
  z <- symcompR::Function("zev")(t)
  expect_true(.is_op(Pre(z), "Pre"))
  ev <- list(list(z ~ 1), list(z ~ Pre(z) - 1))
  sys <- mtkcompile(System(list(D(z) ~ 1), t, name = "saw", continuous_events = ev))
  sol <- solve(ODEProblem(sys, list(z ~ 0), c(0, 2.5)), Tsit5(), reltol = 1e-9, abstol = 1e-12)
  expect_equal(sol(2.4)[[1]], 0.4, tolerance = 1e-6)
})

test_that("bindings on retained unknowns are enforced at initialization", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t)
  sys <- mtkcompile(System(list(D(x) ~ -x, D(y) ~ -y), t, name = "b2",
                           bindings = list(x ~ y + 1)))
  prob <- ODEProblem(sys, list(y ~ 2), c(0, 1))
  expect_equal(prob[x], 3)
})

test_that("the index cache follows the problem parameter layout", {
  t <- t_nounits
  x <- Variable("x", t); k <- Parameter("kk"); q <- Parameter("qq", default = 2 * k)
  sys <- mtkcompile(System(list(D_nounits(x) ~ -q * x), t, name = "ic"))
  ic <- get_index_cache(sys)
  expect_equal(names(ic$param_index), .keys(.problem_params(sys)))
  prob <- ODEProblem(sys, list(x ~ 1, k ~ 3), c(0, 1))
  expect_equal(unname(prob$p[ic$param_index[[.key(k)]]]), 3)
})

test_that("linearization at a solution uses parameters changed by events", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); k <- Parameter("klin", default = 1)
  dose <- SymbolicDiscreteCallback(c(0.5), list(k ~ 2), discrete_parameters = list(k))
  sys <- System(list(D(x) ~ -k * x), t, name = "lin_ev", discrete_events = list(dose))
  csys <- mtkcompile(sys)
  sol <- solve(ODEProblem(csys, list(x ~ 1), c(0, 1)), Tsit5())
  op <- LinearizationOpPoint(sol, 0.9)
  lin <- linearize(sys, list(), list(x), op = op)
  expect_equal(c(lin$A), -2, tolerance = 1e-8)
})
