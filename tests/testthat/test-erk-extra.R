# Reference values from OrdinaryDiffEq (Julia), u' = (u2, -u1), u(0) = (0, 1).
osc <- function(tf = 1) ODEProblem(function(u, p, t) c(u[2], -u[1]), c(0, 1), c(0, tf))

test_that("parametric methods reproduce OrdinaryDiffEq", {
  prob <- osc()
  ref <- list(Anas5 = c(0.8414709847129438, 0.5403023059702371),
              pRRK22 = c(0.8429951874439476, 0.5424217179809776),
              pRRK33 = c(0.8415126846366356, 0.5402831687445417),
              pRRK54 = c(0.8414708097663632, 0.5403017688798579))
  algs <- list(Anas5 = Anas5(w = 1.3), pRRK22 = pRRK22(kappa = 0.5), pRRK33 = pRRK33(kappa = 0.5),
               pRRK54 = pRRK54(kappa = 0.5))
  for (nm in names(algs))
    expect_equal(as.numeric(solve(prob, algs[[nm]], dt = 0.1, adaptive = FALSE)(1)), ref[[nm]], tolerance = 1e-13,
                 info = nm)
  # FRK65 with omega != 0 diverges with fixed steps in OrdinaryDiffEq too
  s <- solve(prob, FRK65(omega = 1.3), dt = 0.1, adaptive = FALSE)
  expect_equal(s(1)[1], -1.0236958555327328e19, tolerance = 1e-10)
  a <- solve(prob, FRK65(omega = 1.3), abstol = 1e-8, reltol = 1e-8)
  expect_length(a$t, 389L)
  expect_equal(a(1)[1], 0.8414704548928075, tolerance = 1e-10)
})

test_that("SSP multistep methods reproduce OrdinaryDiffEq", {
  prob <- ODEProblem(function(u, p, t) c(u[2], -u[1] + 0.1 * t), c(0, 1), c(0, 1))
  expect_equal(as.numeric(solve(prob, SSPRKMSVS32(), dt = 0.05)(1)), c(0.8578828731328382, 0.5857687510422109),
               tolerance = 1e-13)
  expect_equal(as.numeric(solve(prob, SSPRKMSVS43(), dt = 0.05)(1)), c(0.8574292822032692, 0.5863272018962368),
               tolerance = 1e-13)
})

test_that("SIMD tableau methods reproduce OrdinaryDiffEq", {
  prob <- osc()
  s <- solve(prob, MER5v2(), abstol = 1e-8, reltol = 1e-8)
  expect_length(s$t, 11L)
  expect_equal(s(1)[1], 0.8414709850198119, tolerance = 1e-12)
  expect_equal(solve(prob, MER5v2(), dt = 0.1, adaptive = FALSE)(1)[1], 0.841470984832923, tolerance = 1e-13)
  s <- solve(prob, MER6v2(), abstol = 1e-8, reltol = 1e-8)
  expect_length(s$t, 9L)
  expect_equal(s(1)[1], 0.8414709876346385, tolerance = 1e-12)
  expect_equal(solve(prob, RK6v4(), abstol = 1e-8, reltol = 1e-8)(1)[1], sin(1), tolerance = 1e-7)
})

test_that("ExplicitRK runs a tableau", {
  prob <- osc()
  expect_equal(solve(prob, ExplicitRK(), abstol = 1e-8, reltol = 1e-8)(1)[1], sin(1), tolerance = 1e-7)
  tab <- ExplicitRKTableau(matrix(c(0, 0.5, 0, 0), 2), c(0, 0.5), c(0, 1), 2)
  expect_equal(solve(prob, ExplicitRK(tab), dt = 0.01)(1)[1], sin(1), tolerance = 1e-4)
})

test_that("sensitivity backends select exact vector Jacobian products", {
  expect_s3_class(ZygoteVJP(), "VJPChoice")
  expect_s3_class(ZygoteAdjoint(), "ForwardSensitivityAlgorithm")
  t <- t_nounits; D <- D_nounits
  x <- Variable("backend_x", t); k <- Parameter("backend_k")
  sys <- mtkcompile(System(list(D(x) ~ -k * x), t, name = "backend"))
  prob <- ODEProblem(sys, list(x ~ 1, k ~ 2), c(0, 1))
  loss <- function(sol) sum(sol[x])
  g1 <- solve_gradient(loss, prob, Tsit5(), saveat = 0.5, sensealg = InterpolatingAdjoint(autojacvec = ZygoteVJP()))
  g2 <- solve_gradient(loss, prob, Tsit5(), saveat = 0.5, sensealg = InterpolatingAdjoint())
  expect_equal(g1$p, g2$p)
  # the discretized solve's derivatives agree with the adjoint at tight tolerances
  g4 <- solve_gradient(loss, prob, Tsit5(), saveat = 0.5, reltol = 1e-10, abstol = 1e-10,
                       sensealg = InterpolatingAdjoint())
  g3 <- solve_gradient(loss, prob, Tsit5(), saveat = 0.5, reltol = 1e-10, abstol = 1e-10,
                       sensealg = ZygoteAdjoint())
  expect_equal(unname(g3$p), unname(g4$p), tolerance = 1e-6)
})
