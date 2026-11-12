# Boundary value problems: collocation and shooting against y = sin(t).

test_that("y'' = -y, y(0) = 0, y(pi/2) = 1 gives sin(t) for every solver", {
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t)
  sys <- mtkcompile(System(list(D(D(y)) ~ -y), t, name = "bvp",
                           constraints = list(EvalAt(0)(y) ~ 0, EvalAt(pi / 2)(y) ~ 1)))
  bvp <- BVProblem(sys, list(), c(0, pi / 2), guesses = list(y ~ 0, D(y) ~ 0.5))
  expect_s3_class(bvp, "BVProblem")
  expect_length(bvp$u0_idxs, 0L)
  ts <- c(0.2, 0.7, 1.3)
  for (alg in list(MIRK4(), Shooting(Tsit5()), MultipleShooting(4, Tsit5()))) {
    sol <- solve(bvp, alg, dt = 0.05)
    expect_s3_class(sol, "BVPSolution")
    expect_equal(sol$retcode, "Success")
    expect_equal(unname(sol(ts)[, 1]), sin(ts), tolerance = 1e-6)
    expect_equal(unname(sol[default_toterm(D(y))]), cos(sol$t), tolerance = 1e-6)
  }
})

test_that("MIRK4 converges with fourth order", {
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t)
  sys <- mtkcompile(System(list(D(D(y)) ~ -y), t, name = "bvp4",
                           constraints = list(EvalAt(0)(y) ~ 0, EvalAt(pi / 2)(y) ~ 1)))
  bvp <- BVProblem(sys, list(), c(0, pi / 2))
  err <- vapply(c(pi / 20, pi / 40), function(dt) {
    sol <- solve(bvp, MIRK4(), dt = dt)
    max(abs(sol[y] - sin(sol$t)))
  }, 0)
  expect_gt(log2(err[1] / err[2]), 3.5)
})

test_that("fixed initial values from op and interior constraints", {
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t)
  # y(0) = 0 from op, y(0.5) = 1 inside the interval: y = sin(t) / sin(0.5)
  sys <- mtkcompile(System(list(D(D(y)) ~ -y), t, name = "bvpi",
                           constraints = list(EvalAt(0.5)(y) ~ 1)))
  bvp <- BVProblem(sys, list(y ~ 0), c(0, 1), guesses = list(D(y) ~ 1))
  expect_equal(bvp$u0_idxs, 1L)
  for (alg in list(MIRK4(), Shooting(Tsit5()))) {
    sol <- solve(bvp, alg, dt = 0.02)
    expect_equal(sol(1)[[1]], sin(1) / sin(0.5), tolerance = 1e-6)
    expect_equal(sol[y][1], 0)
  }
})

test_that("consistent overdetermined conditions are solved in least squares", {
  t <- t_nounits; D <- D_nounits
  y <- Variable("y", t)
  sys <- mtkcompile(System(list(D(D(y)) ~ -y), t, name = "bvpo",
                           constraints = list(EvalAt(pi / 2)(y) ~ 1)))
  expect_warning(bvp <- BVProblem(sys, list(y ~ 0, D(y) ~ 1), c(0, pi / 2)), "overdetermined")
  for (alg in list(MIRK4(), Shooting(Tsit5()))) {
    sol <- solve(bvp, alg, dt = 0.05)
    expect_equal(sol(1)[[1]], sin(1), tolerance = 1e-6)
  }
})

test_that("boundary value problems from R functions", {
  f <- function(u, p, t) c(u[2], -u[1])
  bc <- function(sol, p, t) c(sol(t[1])[1], sol(t[2])[1] - 1)
  bvp <- BVProblem(f, bc, c(0, 0), c(0, pi / 2))
  sol <- solve(bvp, MIRK4(), dt = 0.05)
  expect_equal(sol(pi / 4)[[1]], sin(pi / 4), tolerance = 1e-6)
  expect_error(solve(bvp, MIRK4()), "dt")
})

test_that("MIRK5 and MIRK6 converge with their orders on fixed meshes", {
  f <- function(u, p, t) c(u[2], -u[1])
  bc <- function(sol, p, t) c(sol(0)[1], sol(pi / 2)[1] - 1)
  bvp <- BVProblem(f, bc, c(0, 0), c(0, pi / 2))
  for (case in list(list(MIRK4(), 4), list(MIRK5(), 5), list(MIRK6(), 6))) {
    err <- vapply(c(pi / 8, pi / 16), function(dt) {
      sol <- solve(bvp, case[[1]], dt = dt, adaptive = FALSE)
      max(abs(sol[1] - sin(sol$t)))
    }, 0)
    expect_equal(log2(err[1] / err[2]), case[[2]], tolerance = 0.05, label = case[[1]]$name)
    # the continuous extension between the nodes
    sol <- solve(bvp, case[[1]], dt = pi / 16, adaptive = FALSE)
    tt <- seq(0.01, pi / 2 - 0.01, length.out = 37)
    expect_lt(max(abs(sol(tt)[, 1] - sin(tt))), 1e-5, label = case[[1]]$name)
  }
  # the weights of the continuous extensions: w(0) = 0, w'(0) = e1, and at
  # tau = 1 the discrete weights for MIRK4 and MIRK6 (MIRK5's extension is
  # only accurate there to its order, as in BoundaryValueDiffEq)
  for (nm in c("MIRK4", "MIRK5", "MIRK6")) {
    tab <- .MIRK_TABS[[nm]]
    expect_equal(.mirk_weights(nm, 0)$w, rep(0, tab$s_star), label = nm)
    expect_equal(.mirk_weights(nm, 0)$wp, c(1, rep(0, tab$s_star - 1L)), label = nm)
    expect_equal(sum(.mirk_weights(nm, 0.37)$w), 0.37, label = nm)
  }
  expect_equal(.mirk_weights("MIRK4", 1)$w, c(.MIRK_TABS$MIRK4$b, 0))
  expect_equal(.mirk_weights("MIRK6", 1)$w, c(.MIRK_TABS$MIRK6$b, 0, 0, 0, 0), tolerance = 1e-12)
})

test_that("defect control refines the mesh of y'' = 100 y with boundary layers", {
  # y(0) = y(1) = 1: y = cosh(10 (t - 1/2)) / cosh(5)
  bvp <- BVProblem(function(u, p, t) c(u[2], 100 * u[1]),
                   function(sol, p, t) c(sol(0)[1] - 1, sol(1)[1] - 1), c(1, 0), c(0, 1))
  ex <- function(t) cosh(10 * (t - 0.5)) / cosh(5)
  tt <- seq(0, 1, length.out = 101)
  fixed <- solve(bvp, MIRK4(), dt = 0.25, adaptive = FALSE)
  expect_gt(max(abs(fixed(tt)[, 1] - ex(tt))), 1e-2)
  for (alg in list(MIRK4(), MIRK5(), MIRK6())) {
    sol <- solve(bvp, alg, dt = 0.25)
    expect_equal(sol$retcode, "Success", label = alg$name)
    expect_gt(length(sol$t), 5L)
    expect_lt(max(abs(sol(tt)[, 1] - ex(tt))), 1e-6)
    expect_equal(sol[2], 10 * sinh(10 * (sol$t - 0.5)) / cosh(5), tolerance = 1e-6, label = alg$name)
    # the final defect is below abstol
    tab <- .MIRK_TABS[[alg$name]]
    mesh <- sol$t; h <- diff(mesh); Ym <- t(as.matrix(sol))
    K <- .mirk_stages(tab, bvp$f$f, bvp$p, mesh, h, Ym, extra = TRUE)
    expect_lt(max(.mirk_defect(tab, alg$name, bvp$f$f, bvp$p, mesh, h, Ym, K)), 1e-6)
  }
  # the mesh is denser in the boundary layers
  sol <- solve(bvp, MIRK4(), dt = 0.25)
  h <- diff(sol$t)
  expect_lt(min(h[c(1, length(h))]), min(h[sol$t[-1] > 0.4 & sol$t[-1] < 0.6]))
  # a cap on the number of intervals fails
  expect_warning(capped <- solve(bvp, MIRK4(max_num_subintervals = 8L), dt = 0.25), "converge")
  expect_equal(capped$retcode, "Failure")
  # a looser tolerance needs fewer points
  expect_lt(length(solve(bvp, MIRK4(), dt = 0.25, abstol = 1e-3)$t), length(sol$t))
})

test_that("mesh selection follows BoundaryValueDiffEq", {
  # an even defect halves the mesh
  m <- c(0, 0.5, 1)
  expect_equal(.mirk_mesh_selector(m, diff(m), c(1e-3, 1e-3), 4L, 1e-6, 3000L), c(0, 0.25, 0.5, 0.75, 1))
  # an uneven defect redistributes: more points where the defect is large
  m <- seq(0, 1, length.out = 5)
  nm <- .mirk_mesh_selector(m, diff(m), c(1, 1e-6, 1e-6, 1e-6), 4L, 1e-6, 3000L)
  expect_equal(range(nm), c(0, 1))
  expect_true(all(diff(nm) > 0))
  expect_gt(sum(nm < 0.25), sum(nm > 0.75))
  expect_null(.mirk_mesh_selector(m, diff(m), c(1, 1e-6, 1e-6, 1e-6), 4L, 1e-6, 3L))
  expect_equal(.mirk_half_mesh(c(0, 1, 3)), c(0, 0.5, 1, 2, 3))
})
