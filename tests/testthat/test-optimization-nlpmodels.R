# Reference values from tools/nlpmodels-probe.jl (NLPModels 0.21).

test_that("NLPModelMeta partitions the variables and constraints as NLPModels.jl", {
  m1 <- NLPModelMeta(4, lvar = c(0, -Inf, 1, -2), uvar = c(0, 2, Inf, 3), ncon = 3,
                     lcon = c(0, -Inf, 1), ucon = c(0, 1, Inf), lin = 2)
  expect_equal(list(m1$ifix, m1$ilow, m1$iupp, m1$irng, m1$ifree, m1$iinf),
               list(1L, 3L, 2L, 4L, integer(), integer()))
  expect_equal(list(m1$jfix, m1$jlow, m1$jupp, m1$jrng, m1$jfree, m1$jinf),
               list(1L, 3L, 2L, integer(), integer(), integer()))
  expect_equal(c(m1$nnzj, m1$nnzh, m1$nlin, m1$nnln), c(12L, 10L, 1L, 2L))
  expect_equal(m1$nln, c(1L, 3L))
  m2 <- NLPModelMeta(m1, minimize = FALSE, name = "copy")
  expect_false(m2$minimize)
  expect_equal(m2$name, "copy")
  expect_equal(m2$irng, 4L)
  expect_equal(c(NLPModels$has_bounds(m1), NLPModels$bound_constrained(m1), NLPModels$unconstrained(m1),
                 NLPModels$linearly_constrained(m1), NLPModels$equality_constrained(m1),
                 NLPModels$inequality_constrained(m1), NLPModels$has_equalities(m1),
                 NLPModels$has_inequalities(m1)),
               c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))
  m3 <- NLPModelMeta(2)
  expect_equal(list(NLPModels$has_bounds(m3), NLPModels$unconstrained(m3), m3$nnzj, m3$nnzh, m3$jac_available),
               list(FALSE, TRUE, 0L, 3L, FALSE))
  ns <- NLSMeta(5, 3, lin = c(1, 4))
  expect_equal(list(ns$nnzj, ns$nnzh, ns$nln, ns$nnln, ns$nlin), list(15L, 6L, c(2L, 3L, 5L), 3L, 2L))
  expect_error(NLPModelMeta(0), "Nonsensical")
  expect_error(NLPModelMeta(2, x0 = 1), "DimensionError")
})

test_that("evaluations update the counters as NLPModels.jl", {
  meta <- NLPModelMeta(2, ncon = 1, lcon = -Inf, ucon = 1, x0 = c(0.5, 0.5))
  toy <- AbstractNLPModel(meta,
    obj = function(x) (1 - x[1])^2 + 10 * (x[2] - x[1]^2)^2,
    grad = function(x) c(-2 * (1 - x[1]) - 40 * x[1] * (x[2] - x[1]^2), 20 * (x[2] - x[1]^2)),
    cons = function(x) x[1]^2 + x[2]^2,
    jac = function(x) matrix(c(2 * x[1], 2 * x[2]), 1))
  x <- c(0.3, -0.4)
  expect_equal(NLPModels$obj(toy, x), 2.8909999999999996, tolerance = 1e-15)
  expect_equal(NLPModels$grad(toy, x), c(4.48, -9.8), tolerance = 1e-15)
  expect_equal(NLPModels$cons(toy, x), 0.25, tolerance = 1e-15)
  expect_equal(NLPModels$jac_coord(toy, x), c(0.6, -0.8), tolerance = 1e-15)
  expect_equal(NLPModels$jprod(toy, x, c(1, 2)), -1, tolerance = 1e-15)
  expect_equal(NLPModels$jtprod(toy, x, 3), c(1.7999999999999998, -2.4000000000000004), tolerance = 1e-15)
  cnt <- toy$counters
  expect_equal(c(cnt$neval_obj, cnt$neval_grad, cnt$neval_cons, cnt$neval_cons_nln, cnt$neval_jac,
                 cnt$neval_jac_nln, cnt$neval_jprod, cnt$neval_jprod_nln, cnt$neval_jtprod, cnt$neval_jtprod_nln),
               rep(1L, 10))
  expect_equal(NLPModels$sum_counters(toy), 6L)
  NLPModels$`reset!`(toy)
  expect_equal(NLPModels$neval_obj(toy), 0L)

  # least squares: grad! counts the gradient twice (src/nls/api.jl)
  nls <- AbstractNLSModel(NLPModelMeta(2, x0 = c(-1.2, 1)), NLSMeta(2, 2),
    residual = function(x) c(1 - x[1], 10 * (x[2] - x[1]^2)),
    jac_residual = function(x) matrix(c(-1, -20 * x[1], 0, 10), 2))
  expect_equal(NLPModels$obj(nls, x), 12.250000000000002, tolerance = 1e-15)
  expect_equal(NLPModels$grad(nls, x), c(28.700000000000003, -49), tolerance = 1e-15)
  nc <- nls$counters
  expect_equal(c(nc$neval_obj, nc$neval_grad, nc$neval_residual, nc$neval_jtprod_residual), c(1L, 2L, 2L, 1L))
  expect_equal(NLPModels$sum_counters(nls), 6L)
})

test_that("models convert to OptimizationProblem and back through the adaptor", {
  # Hock and Schittkowski problem 13: box and inequality constraints
  meta <- NLPModelMeta(2, x0 = c(-2, -2), lvar = c(0, 0), uvar = c(Inf, Inf), ncon = 1, lcon = 0, ucon = Inf,
                       name = "HS13")
  hs13 <- AbstractNLPModel(meta, obj = function(x) (x[1] - 2)^2 + x[2]^2,
                           cons = function(x) (1 - x[1])^3 - x[2])
  prob <- OptimizationProblem(hs13)
  expect_equal(prob$lb, c(0, 0))
  expect_equal(prob$lcons, 0)
  expect_equal(prob$ucons, Inf)
  expect_equal(as.numeric(prob$u0), c(-2, -2))
  # an unconstrained model maximized: the sense carries over
  mx <- AbstractNLPModel(NLPModelMeta(1, x0 = 0, minimize = FALSE), obj = function(x) -(x - 3)^2)
  sol <- solve(OptimizationProblem(mx), LBFGSB())
  expect_equal(as.numeric(sol$u), 3, tolerance = 1e-6)
  expect_equal(sol$objective, 0, tolerance = 1e-10)

  # build_nlpmodel_meta of a problem with bounds and two constraints, maximized
  cp <- OptimizationProblem(function(u, p) sum(u), c(1, 2, 3), lb = c(-1, -1, -1), ub = c(1, 1, 1),
                            cons = function(u, p) c(u[1] - u[2], u[3]), lcons = c(0, -Inf), ucons = c(0, 1),
                            sense = MaxSense)
  m4 <- build_nlpmodel_meta(init(cp, AugLag()))
  expect_equal(list(m4$irng, m4$jfix, m4$jupp, m4$nnzj, m4$nnzh, m4$nln, m4$minimize),
               list(1:3, 1L, 2L, 6L, 6L, c(1L, 2L), FALSE))
  ad <- NLPModelsAdaptor(init(cp, AugLag()), m4, Counters())
  expect_equal(NLPModels$jac_structure(ad), list(rows = c(1L, 2L, 1L, 2L, 1L, 2L), cols = c(1L, 1L, 2L, 2L, 3L, 3L)))
  expect_equal(NLPModels$jac_coord(ad, c(0.1, 0.2, 0.3)), c(1, 0, -1, 0, 0, 1), tolerance = 1e-6)
  expect_equal(NLPModels$hess_structure(ad), list(rows = c(1L, 2L, 3L, 2L, 3L, 3L), cols = c(1L, 1L, 1L, 2L, 2L, 3L)))
  expect_equal(NLPModels$cons(ad, c(0.1, 0.2, 0.3)), c(-0.1, 0.3))
  expect_equal(NLPModels$neval_cons(ad), 0L)
})
