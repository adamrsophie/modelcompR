# modelingtoolkitize of numerical problems (modelingtoolkitize/*.jl).

test_that("the Lotka Volterra ODEProblem becomes the symbolic system", {
  lotka <- function(u, p, t) {
    x <- u[1]; y <- u[2]
    c(p[1] * x - p[2] * x * y, -p[3] * y + p[4] * x * y)
  }
  prob <- ODEProblem(lotka, c(1, 1), c(0, 10), c(1.5, 1, 3, 1))
  sys <- modelingtoolkitize(prob, name = "lv")
  expect_identical(nameof(sys), "lv")
  expect_identical(vapply(unknowns(sys), format, ""), c("x₁(t)", "x₂(t)"))
  expect_identical(vapply(parameters(sys), format, ""), paste0("α", c("₁", "₂", "₃", "₄")))
  u <- unknowns(sys); p <- parameters(sys)
  eqs <- equations(sys)
  expect_true(.is_zero(symcompR::expand(eqs[[1L]]$rhs - (p[[1L]] * u[[1L]] - p[[2L]] * u[[1L]] * u[[2L]]))))
  expect_true(.is_zero(symcompR::expand(eqs[[2L]]$rhs - (-p[[3L]] * u[[2L]] + p[[4L]] * u[[1L]] * u[[2L]]))))
  expect_equal(smap_get(initial_conditions(sys), p[[3L]]), 3)
  # same trajectory as the numerical problem, and the symbolic Jacobian
  cs <- mtkcompile(sys)
  prob2 <- ODEProblem(cs, list(), c(0, 10), jac = TRUE)
  s1 <- solve(prob, reltol = 1e-9, abstol = 1e-9)
  s2 <- solve(prob2, reltol = 1e-9, abstol = 1e-9)
  expect_equal(unname(s2(10, idxs = u)), unname(s1(10)), tolerance = 1e-6)
  # the compiled system may order its unknowns differently
  ord <- match(.keys(u), .keys(unknowns(cs)))
  uval <- numeric(2); uval[ord] <- c(2, 3)
  J <- prob2$f$jac(uval, prob2$p, 0)[ord, ord]
  expect_equal(J, matrix(c(1.5 - 3, 3, -2, -3 + 2), 2, 2))
})

test_that("names come from u0, p, the arguments or the system", {
  pend <- function(u, p, t) {
    du <- numeric(2)
    du[1] <- u[2]
    du[2] <- -p["g"] * sin(u[1])
    du
  }
  prob <- ODEProblem(pend, c(th = 0.1, w = 0), c(0, 1), c(g = 9.81))
  sys <- modelingtoolkitize(prob)
  expect_identical(vapply(unknowns(sys), format, ""), c("th(t)", "w(t)"))
  expect_identical(format(equations(sys)[[2L]]), "Derivative(w(t), t) ~ -g*sin(th(t))")
  sys2 <- modelingtoolkitize(prob, u_names = c("a", "b"), p_names = "k")
  expect_identical(vapply(unknowns(sys2), format, ""), c("a(t)", "b(t)"))
  expect_error(modelingtoolkitize(prob, u_names = "a"), "does not match")
  # list parameters, scalar parameter, numeric first in c(), sums and time
  rob <- function(u, p, t) list(-p$k1 * u[1], p$k1 * u[1] - p$k2 * u[2]^2, p$k2 * u[2]^2)
  s3 <- modelingtoolkitize(ODEProblem(rob, c(1, 0, 0), c(0, 1), list(k1 = 0.04, k2 = 3e7)))
  expect_identical(vapply(parameters(s3), format, ""), c("k1", "k2"))
  f5 <- function(u, p, t) {
    a <- c(v = u[1], w = u[2])
    c(a["w"], -a["v"])
  }
  s5n <- modelingtoolkitize(ODEProblem(f5, c(1, 2), c(0, 1)))
  expect_identical(format(equations(s5n)[[2L]]$rhs), "-x\u2081(t)")
  f4 <- function(u, p, t) c(0, -p * u[1] + sum(u) * t)
  s4 <- modelingtoolkitize(ODEProblem(f4, c(1, 2), c(0, 1), 2))
  expect_identical(vapply(parameters(s4), format, ""), "α")
  expect_true(.is_zero(equations(s4)[[1L]]$rhs))
  # problems built from a system keep its names
  t <- t_nounits; D <- D_nounits
  xx <- Variable("xx", t, bounds = c(0, 2)); kk <- Parameter("kk")
  s0 <- mtkcompile(System(list(D(xx) ~ -kk * xx), t, name = "s0"))
  p0 <- ODEProblem(s0, list(xx ~ 1, kk ~ 2), c(0, 1))
  s5 <- modelingtoolkitize(p0)
  expect_identical(format(equations(s5)[[1L]]), "Derivative(xx(t), t) ~ -kk*xx(t)")
  expect_identical(getbounds(unknowns(s5)[[1L]]), c(0, 2))
})

test_that("mass matrices give algebraic equations", {
  f <- function(u, p, t) c(-u[1] + u[2], u[1] + u[2] - 1)
  prob <- ODEProblem(f, c(0.5, 0.5), c(0, 1), mass_matrix = diag(c(1, 0)))
  sys <- modelingtoolkitize(prob)
  expect_true(.is_zero(equations(sys)[[2L]]$lhs))
  expect_error(modelingtoolkitize(ODEProblem(f, c(0.5, 0.5), c(0, 1), mass_matrix = matrix(c(1, 1, 0, 1), 2))),
               "Non-permutation")
})

test_that("tracing errors are explained", {
  f <- function(u, p, t) if (u[1] > 0) c(u[1]) else c(-u[1])
  expect_error(modelingtoolkitize(ODEProblem(f, 1, c(0, 1))), "could not trace")
})

test_that("nonlinear problems are traced", {
  skip_if_not(exists("NonlinearProblem"))
  nf <- function(u, p) c(u[1]^2 - p[1], u[2] - u[1], 0)
  np <- NonlinearProblem(nf, c(1, 1, 0), c(2))
  ns <- modelingtoolkitize(np)
  expect_length(equations(ns), 2L)
  expect_null(get_iv(ns))
  sol <- solve(NonlinearProblem(mtkcompile(ns), list()))
  expect_equal(sol[unknowns(ns)[[2L]]], sqrt(2), tolerance = 1e-8)
})

test_that("optimization problems are traced", {
  skip_if_not(exists("OptimizationProblem"))
  rosen <- function(u, p) (p[1] - u[1])^2 + p[2] * (u[2] - u[1]^2)^2
  cons <- function(u, p) c(u[1]^2 + u[2]^2)
  op <- tryCatch(OptimizationProblem(rosen, c(0, 0), c(1, 100), lb = c(-1, -1), ub = c(1, 1),
                                     cons = cons, lcons = -Inf, ucons = 1),
                 error = function(e) NULL)
  skip_if(is.null(op), "raw OptimizationProblem() constructor has a different signature")
  os <- modelingtoolkitize(op)
  expect_length(costs(os), 1L)
  expect_length(constraints(os), 1L)
  expect_identical(getbounds(unknowns(os)[[1L]]), c(-1, 1))
  expect_true(inherits(constraints(os)[[1L]], "Inequality"))
})

test_that("SDE problems keep their noise", {
  skip_if_not(exists("SDEProblem"))
  f <- function(u, p, t) -p[1] * u
  g <- function(u, p, t) 0.1 * u
  sp <- SDEProblem(f, g, c(1, 2), c(0, 1), c(0.5))
  ss <- modelingtoolkitize(sp)
  ne <- get_noise_eqs(ss)
  expect_length(ne, 2L)
  expect_identical(format(ne[[2L]]), "0.1*x₂(t)")
  expect_identical(format(equations(ss)[[1L]]), "Derivative(x₁(t), t) ~ -α*x₁(t)")
})
