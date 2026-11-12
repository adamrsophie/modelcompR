# Optimization problems. Reference values are analytic, except the
# constrained Rosenbrock optimum, which is the solution reported by the
# Optimization.jl documentation for the same problem.

rosenbrock_system <- function(constraints = list(), name = "rosenbrock") {
  x <- Variable("x"); y <- Variable("y"); a <- Parameter("a"); b <- Parameter("b")
  list(x = x, y = y, a = a, b = b,
       sys = complete(System(list(), costs = list((a - x)^2 + b * (y - x^2)^2), unknowns = list(x, y),
                             ps = list(a, b), name = name, constraints = constraints)))
}

test_that("the Rosenbrock function is minimized at (a, a^2)", {
  r <- rosenbrock_system()
  prob <- OptimizationProblem(r$sys, list(r$x ~ -1, r$y ~ 2, r$a ~ 1, r$b ~ 100))
  expect_s3_class(prob, "OptimizationProblem")
  expect_null(prob$lb); expect_null(prob$lcons)
  expect_equal(prob$f$f(c(-1, 2), c(1, 100)), 4 + 100)
  for (alg in list(LBFGSB(), LBFGS(), BFGS(), Nlminb(), AugLag())) {
    sol <- solve(prob, alg)
    expect_equal(sol$retcode, "Success", info = alg$name)
    expect_equal(unname(sol$u), c(1, 1), tolerance = 1e-6, info = alg$name)
    expect_lt(sol$objective, 1e-10)
  }
  nm <- solve(prob, NelderMead())
  expect_equal(unname(nm$u), c(1, 1), tolerance = 1e-4)
  sol <- solve(prob)
  expect_equal(sol$alg$name, "LBFGSB")
  expect_equal(names(sol$u), c("x", "y"))
  expect_equal(sol[r$x], 1, tolerance = 1e-6)
  expect_equal(sol[r$x + r$y], 2, tolerance = 1e-6)
  expect_equal(sol$ps[r$b], 100)
  sol2 <- solve(remake(prob, p = list(r$a ~ 2)))
  expect_equal(unname(sol2$u), c(2, 4), tolerance = 1e-5)
})

test_that("symbolic gradients and Hessians of the cost are exact", {
  r <- rosenbrock_system()
  u <- c(-0.3, 0.8); p <- c(1.5, 10)
  g <- c(-2 * (p[1] - u[1]) - 4 * p[2] * u[1] * (u[2] - u[1]^2), 2 * p[2] * (u[2] - u[1]^2))
  H <- rbind(c(2 - 4 * p[2] * (u[2] - u[1]^2) + 8 * p[2] * u[1]^2, -4 * p[2] * u[1]), c(-4 * p[2] * u[1], 2 * p[2]))
  expect_equal(generate_cost(r$sys)(u, p), (p[1] - u[1])^2 + p[2] * (u[2] - u[1]^2)^2)
  expect_equal(generate_cost_gradient(r$sys)(u, p), g)
  expect_equal(generate_cost_hessian(r$sys)(u, p), H)
  expect_length(calculate_cost_gradient(r$sys), 2L)
  expect_equal(dim(calculate_cost_hessian(r$sys)), c(2L, 2L))
  fn <- OptimizationFunction(r$sys, grad = TRUE, hess = TRUE)
  expect_equal(fn$grad(u, p), g)
  expect_equal(fn$hess(u, p), H)
  expect_null(fn$cons)
  expect_null(OptimizationFunction(r$sys)$grad)
})

test_that("the Rosenbrock function constrained to the unit disk", {
  x <- Variable("x"); y <- Variable("y")
  r <- rosenbrock_system(list(Inequality(x^2 + y^2, 1)))
  prob <- OptimizationProblem(r$sys, list(r$x ~ 0, r$y ~ 0, r$a ~ 1, r$b ~ 100))
  expect_equal(prob$lcons, -Inf); expect_equal(prob$ucons, 0)
  expect_equal(prob$f$cons(c(1, 2), c(1, 100)), 4)
  expect_equal(prob$f$cons_j(c(1, 2), c(1, 100)), rbind(c(2, 4)))
  sol <- solve(prob)
  expect_equal(sol$alg$name, "AugLag")
  expect_equal(sol$retcode, "Success")
  expect_equal(unname(sol$u), c(0.7864151541684254, 0.6176983125233897), tolerance = 1e-6)
  expect_equal(sol$objective, 0.045674808764, tolerance = 1e-6)
  expect_equal(sum(sol$u^2), 1, tolerance = 1e-7)
  expect_error(solve(prob, LBFGSB()), "does not handle constraints")
  # the same constraint written with >= has the same residual
  expect_equal(format(constraint_residual(Inequality(1, x^2 + y^2, ">="))),
               format(constraint_residual(Inequality(x^2 + y^2, 1))))
  expect_equal(format(canonical_constraints(r$sys)[[1L]]), format(x^2 + y^2 - 1))
  expect_length(calculate_constraint_hessian(r$sys), 1L)
  expect_equal(generate_constraint_hessian(r$sys)(c(1, 2), c(1, 100))[[1L]], diag(2, 2))
  expect_equal(dim(calculate_constraint_jacobian(r$sys)), c(1L, 2L))
})

test_that("equality constraints, mixed constraints and bounds", {
  x <- Variable("x"); y <- Variable("y")
  # the point of the line x + y = 1 closest to the origin
  es <- complete(System(list(), costs = list(x^2 + y^2), name = "eq", constraints = list(x + y ~ 1)))
  ep <- OptimizationProblem(es, list(x ~ 0, y ~ 0))
  expect_equal(ep$lcons, 0)
  sol <- solve(ep)
  expect_equal(sol$retcode, "Success")
  expect_equal(unname(sol$u), c(0.5, 0.5), tolerance = 1e-7)
  # adding x <= 0.2 moves it to (0.2, 0.8)
  ms <- complete(System(list(), costs = list(x^2 + y^2), name = "mix",
                        constraints = list(x + y ~ 1, Inequality(x, 0.2))))
  sol <- solve(OptimizationProblem(ms, list(x ~ 0, y ~ 0)))
  expect_equal(unname(sol$u), c(0.2, 0.8), tolerance = 1e-6)
  # bounds from metadata: the Rosenbrock minimum with x <= 0.5 is (0.5, 0.25)
  xb <- Variable("xb", bounds = c(-1, 0.5)); yb <- Variable("yb")
  bs <- complete(System(list(), costs = list((1 - xb)^2 + 100 * (yb - xb^2)^2), name = "rb"))
  bp <- OptimizationProblem(bs, list(xb ~ 0, yb ~ 0))
  expect_equal(bp$lb, c(-1, -Inf)); expect_equal(bp$ub, c(0.5, Inf))
  for (alg in list(LBFGSB(), Nlminb(), AugLag())) {
    sol <- solve(bp, alg)
    expect_equal(unname(sol$u), c(0.5, 0.25), tolerance = 1e-6, info = alg$name)
    expect_equal(sol$objective, 0.25, tolerance = 1e-8)
  }
  # as in OptimizationOptimJL, a method without bounds runs inside Fminbox
  fb <- solve(bp, BFGS())
  expect_s3_class(fb$alg, "Fminbox")
  expect_equal(unname(fb$u), c(0.5, 0.25), tolerance = 1e-6)
  expect_error(OptimizationProblem(bs, list(xb ~ 0, yb ~ 0), lb = c(0, 0)), "both `lb` and `ub`")
  bp2 <- OptimizationProblem(bs, list(xb ~ 0, yb ~ 0), lb = c(-2, -2), ub = c(0, 2))
  expect_equal(unname(solve(bp2)$u), c(0, 0), tolerance = 1e-6)
  # an infeasible pair of constraints is reported
  inf <- complete(System(list(), costs = list(x^2), name = "inf", constraints = list(x ~ 1, x ~ 2)))
  expect_equal(solve(OptimizationProblem(inf, list(x ~ 0)))$retcode, "Infeasible")
})

test_that("constraints become penalties", {
  x <- Variable("x")
  sys <- complete(System(list(), costs = list((x - 2)^2), name = "opt", constraints = list(Inequality(x, 1))))
  pen <- constraints_to_penalties(sys, weights = 1e3)
  expect_length(constraints(pen), 0L)
  expect_length(costs(pen), 2L)
  # minimum of (x - 2)^2 + w * (x - 1)^2 for x > 1 is (2 + w) / (1 + w)
  sol <- solve(OptimizationProblem(pen, list(x ~ 0)))
  expect_equal(unname(sol$u), (2 + 1e3) / (1 + 1e3), tolerance = 1e-6)
  w <- Parameter("w")
  pen2 <- constraints_to_penalties(sys, weights = w)
  expect_true(any(vapply(parameters(pen2), function(p) identical(format(p), "w"), TRUE)))
  sol2 <- solve(OptimizationProblem(pen2, list(x ~ 0, w ~ 3)))
  expect_equal(unname(sol2$u), 5 / 4, tolerance = 1e-6)
  expect_error(constraints_to_penalties(sys, weights = c(1, 2)), "one entry per constraint")
})

test_that("multiobjective functions keep the costs separate", {
  x <- Variable("x"); y <- Variable("y"); w1 <- Parameter("w1")
  sys <- complete(System(list(), costs = list((x - 1)^2, (y + 2)^2), ps = list(w1), name = "mo"))
  prob <- OptimizationProblem(sys, list(x ~ 0, y ~ 0, w1 ~ 2), multiobjective = TRUE)
  expect_s3_class(prob$f, "MultiObjectiveOptimizationFunction")
  expect_equal(prob$f$f(c(0, 0), 2), c(1, 4))
  expect_equal(prob$f$jac(c(0, 0), 2), rbind(c(-2, 0), c(0, 4)))
  expect_equal(generate_multiobjective_cost(sys)(c(2, 0), 2), c(1, 4))
  H <- generate_multiobjective_hessian(sys)
  expect_length(H, 2L)
  expect_equal(H[[2L]](c(0, 0), 2), rbind(c(0, 0), c(0, 2)))
  expect_length(calculate_multiobjective_hessian(sys), 2L)
  expect_equal(dim(calculate_multiobjective_jacobian(sys)), c(2L, 2L))
  expect_error(solve(prob, LBFGSB()), "multiobjective")
  # the costs are independent, so the Pareto front is the single point (1, -2)
  front <- solve(prob, WeightedSumScan(n = 3))
  expect_equal(unname(front$u[1L, ]), c(1, -2), tolerance = 1e-5)
  # a weighted sum with a symbolic weight is a scalar problem
  wp <- OptimizationProblem(sys, list(x ~ 0, y ~ 0, w1 ~ 2), weights = list(w1, 1))
  expect_equal(wp$f$f(c(0, 0), 2), 2 * 1 + 4)
  expect_equal(unname(solve(wp)$u), c(1, -2), tolerance = 1e-6)
  expect_error(OptimizationProblem(sys, list(x ~ 0, y ~ 0, w1 ~ 2), weights = list(1)), "one entry per cost")
  expect_error(MultiObjectiveOptimizationFunction(sys, weights = list(1, 1)), "scalarizes")
})

test_that("tunable parameters are listed, subset and reordered", {
  t <- t_nounits
  x <- Variable("x", t)
  a <- Parameter("a"); b <- Parameter("b"); cc <- Parameter("cc", tunable = FALSE)
  sys <- mtkcompile(System(list(D_nounits(x) ~ -a * x + b + cc), t, name = "tun"))
  expect_equal(vapply(tunable_parameters(sys), format, ""), c("a", "b"))
  sub <- subset_tunables(sys, list(b))
  expect_equal(vapply(tunable_parameters(sub), format, ""), "b")
  expect_error(subset_tunables(sys, list(cc)), "subset of the current tunable")
  expect_equal(reorder_dimension_by_tunables(sys, c(10, 20), list(b, a)), c(20, 10))
  M <- rbind(c(1, 2), c(3, 4))
  expect_equal(reorder_dimension_by_tunables(sys, M, list(b, a), dim = 2), M[, c(2, 1)])
  expect_error(reorder_dimension_by_tunables(sys, c(1, 2), list(cc)), "not a tunable parameter")
})

test_that("raw problems and compiled optimization systems work", {
  sol <- solve(OptimizationProblem(function(u, p) sum((u - p)^2), c(0, 0), c(1, 2)))
  expect_equal(sol$u, c(1, 2), tolerance = 1e-7)
  # projection of (1, 2) on the line u1 + u2 = 1
  con <- solve(OptimizationProblem(function(u, p) sum((u - p)^2), c(0, 0), c(1, 2),
                                   cons = function(u, p) u[1] + u[2], lcons = 1, ucons = 1))
  expect_equal(con$u, c(0, 1), tolerance = 1e-6)
  x <- Variable("x"); y <- Variable("y")
  ms <- mtkcompile(System(list(), costs = list((1 - x)^2 + 100 * (y - x^2)^2), name = "rosenbrock"))
  sol <- solve(OptimizationProblem(ms, list(x ~ -1, y ~ 2)))
  expect_equal(unname(sol$u), c(1, 1), tolerance = 1e-6)
  eqs <- complete(System(list(0 ~ x - y), costs = list(x^2), name = "witheqs"))
  expect_error(OptimizationProblem(eqs, list(x ~ 1, y ~ 1)), "does not support systems with equations")
  xi <- Variable("xi", integer = TRUE)
  ip <- OptimizationProblem(complete(System(list(), costs = list((xi - 1.4)^2), name = "int")), list(xi ~ 0))
  expect_true(ip$int)
  expect_equal(unname(solve(ip)$u), 1)
  expect_error(solve(ip, LBFGSB()), "BranchAndBound")
})

test_that("symbolic bounds follow the parameters and function objects carry their system", {
  am <- Parameter("am"); x <- Variable("x", bounds = list(0, am))
  sys <- complete(System(list(), costs = list((x - 5)^2), unknowns = list(x), ps = list(am), name = "sb"))
  prob <- OptimizationProblem(sys, list(x ~ 1, am ~ 10))
  expect_equal(prob$ub, 10)
  expect_equal(unname(solve(prob)$u), 5, tolerance = 1e-7)
  prob2 <- remake(prob, p = list(am ~ 1))
  expect_equal(prob2$ub, 1)
  expect_equal(unname(solve(prob2)$u), 1, tolerance = 1e-7)
  fn <- OptimizationFunction(sys, grad = TRUE)
  prob3 <- OptimizationProblem(fn, list(x ~ 1, am ~ 2))
  expect_s3_class(prob3$sys, "System")
  expect_equal(unname(solve(prob3)$u), 2, tolerance = 1e-7)
})

test_that("multiobjective problems return a Pareto front", {
  # f1 = x^2, f2 = (x - 2)^2: the Pareto set is [0, 2] and the front
  # f2 = (sqrt(f1) - 2)^2
  x <- Variable("opmo_x", bounds = c(-1, 3))
  sys <- complete(System(list(), costs = list(x^2, (x - 2)^2), name = "opmo"))
  prob <- OptimizationProblem(sys, list(x ~ 1), multiobjective = TRUE)
  ws <- solve(prob, WeightedSumScan(n = 5))
  expect_s3_class(ws, "MultiObjectiveOptimizationSolution")
  expect_equal(unname(ws$u[, 1L]), c(0, 0.5, 1, 1.5, 2), tolerance = 1e-5)
  expect_equal(unname(ws$objective[, 2L]), (sqrt(ws$objective[, 1L]) - 2)^2, tolerance = 1e-8)
  expect_equal(ws[x], unname(ws$u[, 1L]))
  ns <- solve(prob, NSGA2(population = 40L, generations = 60L, seed = 1L))
  expect_identical(ns$retcode, "Success")
  expect_true(all(ns$u >= -1e-3 & ns$u <= 2 + 1e-3))
  expect_gt(nrow(ns$u), 20L)
  expect_equal(unname(ns$objective[, 2L]), (sqrt(ns$objective[, 1L]) - 2)^2, tolerance = 1e-6)
  # a seed makes the run reproducible
  expect_equal(solve(prob, NSGA2(population = 40L, generations = 60L, seed = 1L))$u, ns$u)
  # ZDT2: a nonconvex front f2 = 1 - f1^2 that a weighted sum cannot trace
  zs <- lapply(1:3, function(i) Variable(paste0("opz", i), bounds = c(0, 1)))
  g <- 1 + 9 * (zs[[2L]] + zs[[3L]]) / 2
  zsys <- complete(System(list(), costs = list(zs[[1L]], g * (1 - (zs[[1L]] / g)^2)), name = "opzdt2"))
  zp <- OptimizationProblem(zsys, lapply(zs, function(z) z ~ 0.5), multiobjective = TRUE)
  zf <- solve(zp, NSGA2(population = 60L, generations = 150L, seed = 2L))$objective
  expect_lt(stats::median(abs(zf[, 2L] - (1 - zf[, 1L]^2))), 1e-3)
  expect_gt(diff(range(zf[, 1L])), 0.9)
  # constraints restrict the front: x >= 0.5
  csys <- complete(System(list(), costs = list(x^2, (x - 2)^2), name = "opmoc",
                          constraints = list(Inequality(0.5, x))))
  cp <- OptimizationProblem(csys, list(x ~ 1), multiobjective = TRUE)
  expect_gt(min(solve(cp, NSGA2(population = 30L, generations = 40L, seed = 3L))$u), 0.5 - 1e-6)
  expect_equal(min(solve(cp, WeightedSumScan(n = 5))$u), 0.5, tolerance = 1e-5)
  expect_error(solve(OptimizationProblem(complete(System(list(), costs = list(x^2, x), name = "opmof")),
                                         list(x ~ 1), multiobjective = TRUE, lb = -Inf, ub = Inf),
                     NSGA2()), "finite bounds")
})

test_that("branch and bound solves integer and binary problems", {
  xi <- Variable("opbb_x", integer = TRUE, bounds = c(-10, 10)); yi <- Variable("opbb_y")
  isys <- complete(System(list(), costs = list((xi - 2.6)^2 + (yi - 1.3)^2), name = "opbbq"))
  s <- solve(OptimizationProblem(isys, list(xi ~ 0, yi ~ 0)))
  expect_identical(s$retcode, "Success")
  expect_equal(unname(s$u[c("opbb_x", "opbb_y")]), c(3, 1.3), tolerance = 1e-6)
  # knapsack: values 10, 13, 7, 8, weights 3, 4, 2, 3, capacity 7; the best
  # choice (by enumeration) is the first two items, value 23
  bs <- lapply(1:4, function(i) Variable(paste0("opbb_b", i), binary = TRUE))
  v <- c(10, 13, 7, 8); w <- c(3, 4, 2, 3)
  ksys <- complete(System(list(), costs = list(-Reduce(`+`, Map(`*`, v, bs))), name = "opknap",
                          constraints = list(Inequality(Reduce(`+`, Map(`*`, w, bs)), 7))))
  ks <- solve(OptimizationProblem(ksys, lapply(bs, function(b) b ~ 0)), BranchAndBound())
  expect_identical(ks$retcode, "Success")
  expect_equal(unname(ks$u[paste0("opbb_b", 1:4)]), c(1, 1, 0, 0))
  expect_equal(ks$objective, -23)
  best <- max(apply(as.matrix(expand.grid(rep(list(0:1), 4))), 1L, function(z) if (sum(w * z) <= 7) sum(v * z) else 0))
  expect_equal(-ks$objective, best)
  # an infeasible integer problem
  esys <- complete(System(list(), costs = list(xi^2), name = "opbbe",
                          constraints = list(xi ~ 0.5)))
  expect_identical(solve(OptimizationProblem(esys, list(xi ~ 0)))$retcode, "Infeasible")
})
