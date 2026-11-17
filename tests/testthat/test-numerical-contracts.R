# Mathematical invariants and representation equivalences: each case
# checks a result against its definition (an analytic value, the same
# quantity in another representation, or a finite difference) rather than
# against a reference run.

test_that("mass matrices keep their meaning in every representation", {
  f <- function(u, p, t) -u
  for (M in list(matrix(2), Matrix::Diagonal(1, 2), Matrix::Matrix(matrix(2), sparse = TRUE), 2,
                 MatrixOperator(matrix(2)))) {
    s <- solve(ODEProblem(ODEFunction(f, mass_matrix = M), 1, c(0, 1)), Rodas5P(), reltol = 1e-10, abstol = 1e-12)
    expect_equal(unname(tail(as.matrix(s), 1)[1L]), exp(-0.5), tolerance = 1e-7)
  }
  fa <- function(u, p, t) c(-u[1] + u[2], u[1] - 2 * u[2])
  sd <- solve(ODEProblem(ODEFunction(fa, mass_matrix = diag(c(1, 0))), c(1, 0.5), c(0, 1)), Rodas5P())
  ss <- solve(ODEProblem(ODEFunction(fa, mass_matrix = Matrix::Diagonal(2, c(1, 0))), c(1, 0.5), c(0, 1)), Rodas5P())
  expect_equal(as.matrix(ss), as.matrix(sd))
  expect_error(ODEProblem(ODEFunction(f, mass_matrix = diag(2)), 1, c(0, 1)), "1 states")
  expect_error(ODEProblem(ODEFunction(f, mass_matrix = "M"), 1, c(0, 1)), "unsupported mass matrix")
  # the adjoint accepts a perfectly conditioned mass matrix of any scale
  M <- diag(1e-16, 3)
  pr <- ODEProblem(ODEFunction(function(u, p, t) -1e-16 * p[1] * u, mass_matrix = M), rep(1, 3), c(0, 1), 1)
  s <- solve(pr, Rodas5P(), reltol = 1e-9, abstol = 1e-11)
  g <- adjoint_sensitivities(s, t = 1, dgdu_discrete = function(u, p, t, i) rep(1, 3))
  expect_equal(as.numeric(g$dp), -3 * exp(-1), tolerance = 1e-5)
})

test_that("observed chains and parameter maps resolve at any length, cycles are errors", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("nc_obs_x", t)
  z <- lapply(seq_len(120L), function(i) Variable(paste0("nc_obs_z", i), t))
  obs <- c(list(Equation(z[[1L]], x)), lapply(2:120, function(i) Equation(z[[i]], z[[i - 1L]] + 1)))
  s <- complete(System(list(Equation(D(x), z[[120L]])), t, unknowns = list(x), observed = obs, name = "nc_chain"))
  expect_equal(generate_rhs(s)(2, numeric(), 0), 121)
  expect_equal(generate_jacobian(s)(2, numeric(), 0), matrix(1))
  a <- Variable("nc_cyc_a", t); b <- Variable("nc_cyc_b", t)
  expect_error(.subst_fixpoint(a + 1, .rule(list(a, b), list(b + 1, a))), "cyclic")
  p <- lapply(1:61, function(i) Parameter(paste0("nc_chain_p", i)))
  m <- c(list(list(p[[1L]], 1)), lapply(2:61, function(i) list(p[[i]], p[[i - 1L]] + 1)))
  expect_equal(unname(varmap_to_vars(m, p)), 1:61)
})

test_that("compiled functions do not share state and caches follow registrations", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("nc_br_x", t); z <- Variable("nc_br_z", t, guess = 1); q <- Parameter("nc_br_q")
  cs <- mtkcompile(System(list(Equation(D(x), z), Equation(0, z^3 - z - q)), t, name = "nc_branches"),
                   inline_nonlinear_sccs = TRUE)
  f1 <- generate_rhs(cs); f2 <- generate_rhs(cs)
  a <- f1(0, 0, 0)
  f2(0, -1, 0)
  expect_identical(f1(0, 0, 0), a)
  xr <- Variable("nc_reg_x", t)
  h <- register_symbolic("nc_registered_fn", function(x) x^2, derivatives = list(function(x) 2 * x))
  cr <- mtkcompile(System(list(Equation(D(xr), h(xr))), t, name = "nc_registered"))
  expect_equal(.sens_symbolic(cr)$jac(2, numeric(), 0), matrix(4))
  register_symbolic("nc_registered_fn", function(x) x^3, derivatives = list(function(x) 3 * x^2))
  expect_equal(.sens_symbolic(cr)$jac(2, numeric(), 0), matrix(12))
})

test_that("graph searches run on long chains without recursion", {
  n <- 5000L
  sc <- .tarjan_scc(n, function(v) if (v < n) v + 1L else integer())
  expect_length(sc, n)
  expect_identical(sc[[1L]], n)
  # a cycle is one component, and the order matches the recursive search
  expect_identical(.tarjan_scc(4L, function(v) c(2L, 3L, 1L, 3L)[v]), list(c(3L, 2L, 1L), 4L))
  st <- new.env(parent = emptyenv())
  st$m <- integer(n); st$dcolor <- logical(n)
  adj <- lapply(seq_len(n), function(i) if (i < n) c(i, i + 1L) else i)
  for (e in seq_len(n)) { st$dcolor[] <- FALSE; .augment(st, adj, e, rep(TRUE, n)) }
  expect_true(all(st$m > 0L))
})

test_that("dual number reductions keep base R semantics", {
  x <- c(-100, 2, 100)
  r <- .ad_jacobian(function(z) mean(z, trim = 1 / 3), x)
  expect_equal(r$value, 2); expect_equal(r$jacobian, matrix(c(0, 1, 0), 1))
  r <- .ad_jacobian(function(z) median(z), x)
  expect_equal(r$value, 2); expect_equal(r$jacobian, matrix(c(0, 1, 0), 1))
  expect_true(is.na(.ad_jacobian(function(z) max(c(z, NA)), x)$value))
  expect_equal(.ad_jacobian(function(z) sum(c(z, NA), na.rm = TRUE), x)$value, 2)
  # a function that drops the dual numbers is not differentiated as a constant
  expect_null(.ad_try_jacobian(function(z) mean(c(0, z)), 2))
  expect_equal(.ad_try_jacobian(function(z) 3, 2)$jacobian, matrix(0))
})

test_that("event values, step size requests and gradients after failures", {
  cb <- PresetTimeCallback(1, function(i) i$u[1] <- 10)
  s <- solve(ODEProblem(function(u, p, t) 0 * u, 1, c(0, 1), callback = cb), Tsit5())
  expect_equal(c(tail(as.matrix(s), 1)[1L], s(1), s(1 - 1e-9)), c(10, 10, 1))
  cc <- ContinuousCallback(function(u, t, integ) u[1] - 0.5, function(integ) set_proposed_dt(integ, 0.001))
  i <- init(ODEProblem(function(u, p, t) 1, 0, c(0, 1)), Tsit5(), callback = cc, dt = 0.75)
  step_integrator(i); te <- i$t
  step_integrator(i)
  expect_equal(i$t - te, 0.001)
  p <- ODEProblem(function(u, p, t) p[1] * u, 1, c(0, 1), 1)
  l <- function(sol) tail(as.matrix(sol)[, 1], 1)
  expect_error(solve_gradient(l, p, Tsit5(), dt = 0.01, maxiters = 1), "MaxIters")
  expect_equal(unname(solve_gradient(l, p, Tsit5(), dt = 0.01, maxiters = 1, allow_failure = TRUE)$p),
               exp(0.01) * 0.01, tolerance = 1e-6)
  # controls in a differentiated affect
  cb2 <- PresetTimeCallback(0.5, function(i) { i$u <- 2 * i$u; set_proposed_dt(i, 0.01) })
  p2 <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 1, callback = cb2)
  expect_equal(unname(solve_gradient(l, p2, Tsit5(), reltol = 1e-10, abstol = 1e-12)$p), -2 * exp(-1), tolerance = 1e-6)
})

test_that("tolerances and right-hand side shapes are checked", {
  p <- ODEProblem(function(u, p, t) -u, c(1, 2), c(0, 1))
  expect_identical(format(solve(p, reltol = c(1e-6, 1e-8))$retcode), "Success")
  expect_error(solve(p, Tsit5(), reltol = c(1e-6, 1e-8, 1e-4)), "one per state")
  bad <- ODEProblem(function(u, p, t) if (t == 0) c(1, 2) else 1, c(0, 0), c(0, 1))
  expect_error(solve(bad, Tsit5()), "returned 1 values")
  calls <- 0L
  of <- ODEFunction(function(u, p, t) -p * u, vjp = function(v, u, p, t) { calls <<- calls + 1L; -p * v })
  r <- .sens_fns(ODEProblem(of, rep(1, 4), c(0, 1), 1))$vjp(rep(1, 4), rep(1, 4), 1, 0)
  expect_equal(r$u, rep(-1, 4)); expect_equal(calls, 1L)
})

test_that("optimizers respect integer domains and the objective sense", {
  x <- Variable("nc_int_x", integer = TRUE, bounds = c(0, 1))
  pr <- OptimizationProblem(complete(System(list(), costs = list((x - 0.4)^2), name = "nc_int")), list(x ~ 0))
  expect_error(solve(pr, CMAEvolutionStrategyOpt()), "integer")
  expect_equal(unname(solve(pr)$u), 0)
  fn <- MultiObjectiveOptimizationFunction(function(u, p) c(u[1], 2 * u[1]))
  for (sense in list(MaxSense, MinSense)) {
    sol <- solve(OptimizationProblem(fn, 0.5, lb = 0, ub = 1, sense = sense), WeightedSumScan(n = 3))
    best <- if (identical(sense, MaxSense)) 1 else 0
    expect_true(all(abs(sol$u - best) < 1e-6))
    expect_equal(unname(sol$objective[1L, ]), c(best, 2 * best), tolerance = 1e-6)
  }
})

test_that("linear solves agree across operator representations", {
  A <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(1, 2, 3))
  expect_equal(solve(LinearProblem(A, c(1, 2, 3)))$u, c(1, 1, 1))
  U <- cbind(c(1, 0, 0), c(0, 1, 0)); b <- c(1, 2, 3)
  for (C in list(diag(c(1, 0)), diag(2), matrix(0, 2, 2))) {
    ref <- base::solve(diag(3) + U %*% C %*% t(U), b)
    L <- LowRankUpdatedMatrix(diag(3), U, C, U)
    expect_equal(solve(LinearProblem(L, b), LUFactorization())$u, ref)
    expect_equal(solve(LinearProblem(L, b))$u, ref)
  }
  sing <- solve(LinearProblem(LowRankUpdatedMatrix(diag(3), U, diag(c(-1, 0)), U), b))
  expect_false(identical(format(sing$retcode), "Success"))
  used <- 0L
  s <- solve(LinearProblem(diag(c(1, 100)), c(1, 1)), KrylovJL_GMRES(),
             Pl = function(v) { used <<- used + 1L; v / c(1, 100) }, maxiters = 1, reltol = 1e-10)
  expect_gt(used, 0L); expect_equal(s$u, c(1, 0.01))
  set.seed(1)
  A <- diag(1:100)
  e <- solve(EigenvalueProblem(A, nev = 1, which = EigenvalueTarget$SM), ArnoldiMethodJL(krylovdim = 2))
  expect_equal(e$values, 1, tolerance = 1e-8)
  e2 <- solve(EigenvalueProblem(A, nev = 1, which = EigenvalueTarget$SR), ArnoldiMethodJL(krylovdim = 2), maxiters = 2)
  v <- e2$vectors[, 1L]
  res <- sqrt(sum(Mod(A %*% v - e2$values * v)^2)) / sqrt(sum(Mod(v)^2))
  expect_true(identical(e2$retcode, "Success") == (res <= 1e-10 * max(1, Mod(e2$values))))
  expect_equal(e2$residuals, res)
})

test_that("ensembles account for failing trajectories", {
  ep <- EnsembleProblem(ODEProblem(function(u, p, t) -u, 1, c(0, 1)),
                        prob_func = function(prob, i, repeat_) { if (i == 2L) stop("bad trajectory"); prob })
  expect_error(solve(ep, Tsit5(), trajectories = 3), "trajectory 2 failed: bad trajectory")
  s <- solve(ep, Tsit5(), trajectories = 3, on_error = "collect")
  expect_equal(s$retcodes, c("Success", "Error", "Success"))
  expect_equal(s$failures, 1L)
  expect_s3_class(s$u[[2L]], "error")
  em <- EnsembleProblem(ODEProblem(function(u, p, t) -u, 1, c(0, 1)))
  s <- solve(em, Tsit5(), trajectories = 2, maxiters = 2)
  expect_equal(s$failures, 2L)
})

test_that("sparsity patterns and loss cotangents avoid dense work", {
  expect_equal(max(column_coloring(Matrix::Diagonal(3000))), 1L)
  A <- matrix(c(1, 1, 0, 0, 1, 0, 0, 0, 1), 3)
  expect_equal(column_coloring(Matrix::Matrix(A, sparse = TRUE)), column_coloring(A))
  sol <- solve(ODEProblem(function(u, p, t) -u, rep(1, 10), c(0, 1)), Tsit5(), saveat = seq(0, 1, length.out = 1001),
               dense = FALSE)
  calls <- 0L
  ct <- .sens_loss_cotangent(function(s) { calls <<- calls + 1L; sum(tail(as.matrix(s), 1)^2) }, sol)
  expect_lt(calls, 200L)
  expect_equal(ct$dU[1001, ], 2 * unname(tail(as.matrix(sol), 1)[1, ]), tolerance = 1e-8)
  expect_equal(sum(abs(ct$dU[-1001, ])), 0)
  # a loss reading every row still gets every entry
  ct2 <- .sens_loss_cotangent(function(s) sum(as.matrix(s)), sol)
  expect_equal(ct2$dU, matrix(1, 1001, 10), tolerance = 1e-6)
})

test_that("SDE forward sensitivities include an event after the last saved time", {
  f <- function(u, p, t) -p[1] * u; g <- function(u, p, t) 0.1 * u
  cb <- PresetTimeCallback(0.8, function(i) i$p[1] <- i$p[1] * i$p[2])
  loss <- function(sol) sum(as.matrix(sol)) + .sol_env(sol)$p_final[1]^2
  mk <- function(p) SDEProblem(f, g, 1, c(0, 1), p, callback = cb)
  G <- solve_gradient(loss, mk(c(1, 2)), EM(), dt = 0.01, seed = 3, saveat = c(0, 0.5),
                      sensealg = ForwardDiffSensitivity())
  fd <- vapply(1:2, function(j) {
    a <- c(1, 2); b <- a; a[j] <- a[j] + 1e-6; b[j] <- b[j] - 1e-6
    (loss(solve(mk(a), EM(), dt = 0.01, seed = 3, saveat = c(0, 0.5))) -
       loss(solve(mk(b), EM(), dt = 0.01, seed = 3, saveat = c(0, 0.5)))) / 2e-6
  }, 0)
  expect_equal(unname(G$p), fd, tolerance = 1e-6)
})

test_that("inline nonlinear solves follow the branch of each integrator", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("nc_ws_x", t); y <- Variable("nc_ws_y", t); z <- Variable("nc_ws_z", t, guess = 1)
  cs <- mtkcompile(System(list(D(x) ~ 1, D(y) ~ z, 0 ~ z^3 - z - x), t, name = "nc_fold"),
                   inline_nonlinear_sccs = TRUE)
  # z follows the lower root of z^3 - z = x from x = -1 until the fold at
  # x = 2 / (3 sqrt(3)), then the upper root: the integral of z is -3/4
  pr <- ODEProblem(cs, list(x ~ -1, y ~ 0), c(0, 2))
  s <- solve(pr, Rodas5P(), reltol = 1e-10, abstol = 1e-12, tstops = 1 + 2 / (3 * sqrt(3)))
  expect_equal(s(2, idxs = y), -0.75, tolerance = 1e-5)
  p2 <- remake(pr, u0 = list(x ~ 0.5, y ~ 0))
  a <- init(pr, Tsit5(), reltol = 1e-8); b <- init(p2, Tsit5(), reltol = 1e-8)
  for (k in 1:400) { if (a$t < 2) step_integrator(a); if (b$t < 2) step_integrator(b) }
  expect_equal(a$u, tail(solve(pr, Tsit5(), reltol = 1e-8)$u, 1)[[1L]])
  expect_equal(b$u, tail(solve(p2, Tsit5(), reltol = 1e-8)$u, 1)[[1L]])
})

test_that("cyclic parameter definitions are reported", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("nc_cy_x", t); a <- Parameter("nc_cy_a"); b <- Parameter("nc_cy_b")
  s <- mtkcompile(System(list(D(x) ~ -a * b * x), t, name = "nc_cy"))
  expect_error(ODEProblem(s, list(x ~ 1, a ~ b + 1, b ~ a * 2), c(0, 1)), "cyclic definitions")
  expect_error(ODEProblem(s, list(x ~ 1, a ~ 1), c(0, 1)), "no value for parameter")
})

test_that("native arithmetic keeps the tolerances of the reference arithmetic", {
  old <- arithmetic_mode("native")
  on.exit(arithmetic_mode(old))
  expect_identical(arithmetic_mode(), "native")
  stf <- ODEProblem(function(u, p, t) -c(1, 10, 100) * u, c(1, 1, 1), c(0, 3))
  for (alg in list(Rodas5P(), FBDF(), Tsit5(), KenCarp4())) {
    s <- solve(stf, alg, reltol = 1e-8, abstol = 1e-10)
    expect_equal(tail(s$u, 1)[[1L]], exp(-c(1, 10, 100) * 3), tolerance = 1e-6, info = class(alg)[1L])
  }
  # extreme scales: the norms and step control do not overflow or underflow
  for (sc in c(1e-200, 1e200)) {
    s <- solve(ODEProblem(function(u, p, t) -u, sc, c(0, 1)), Tsit5(), reltol = 1e-10, abstol = sc * 1e-12)
    expect_equal(tail(s$u, 1)[[1L]] / sc, exp(-1), tolerance = 1e-8)
  }
  cb <- ContinuousCallback(function(u, t, i) u[1] - 0.5, function(i) i$u[1] <- 1)
  s <- solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1), callback = cb), Tsit5(), reltol = 1e-10, abstol = 1e-12)
  expect_equal(tail(s$u, 1)[[1L]], exp(-(1 - log(2))), tolerance = 1e-8)
  # optimizers: the same optimum, and the KKT conditions of a bounded problem
  f <- function(x, p) sum((x - 1:4)^2) + x[1] * x[2]
  of <- OptimizationFunction(f, AutoForwardDiff())
  r <- solve(OptimizationProblem(of, rep(0, 4)), Newton())
  expect_equal(r$u, base::solve(matrix(c(2, 1, 0, 0, 1, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2), 4), 2 * (1:4)), tolerance = 1e-8)
  rb <- solve(OptimizationProblem(of, rep(0, 4), lb = rep(0, 4), ub = rep(2, 4)), Fminbox())
  g <- c(2 * (rb$u[1] - 1) + rb$u[2], 2 * (rb$u[2] - 2) + rb$u[1], 2 * (rb$u[3:4] - 3:4))
  free <- rb$u > 1e-6 & rb$u < 2 - 1e-6
  expect_lt(max(abs(g[free])), 1e-5)
  expect_true(all(g[rb$u >= 2 - 1e-6] <= 1e-6))
})

test_that("sparse Jacobians stay sparse through the implicit solvers", {
  n <- 400
  L <- Matrix::bandSparse(n, k = -1:1, diagonals = list(rep(1, n - 1), rep(-2, n), rep(1, n - 1))) * 100
  f <- function(u, p, t) as.numeric(L %*% u) - u^3
  u0 <- sin(pi * seq_len(n) / (n + 1))
  ref <- solve(ODEProblem(f, u0, c(0, 0.2)), Rodas5P(), reltol = 1e-9, abstol = 1e-11)
  ua <- tail(ref$u, 1)[[1L]]
  sj <- ODEFunction(f, jac = function(u, p, t) L - Matrix::Diagonal(n, 3 * u^2))
  for (alg in list(Rodas5P(), FBDF(), KenCarp4())) {
    s <- solve(ODEProblem(sj, u0, c(0, 0.2)), alg, reltol = 1e-9, abstol = 1e-11)
    expect_equal(tail(s$u, 1)[[1L]], ua, tolerance = 1e-6, info = class(alg)[1L])
  }
  s <- solve(ODEProblem(ODEFunction(f, sparsity = L), u0, c(0, 0.2)),
             Rodas5P(autodiff = AutoSparse(AutoFiniteDiff())), reltol = 1e-9, abstol = 1e-11)
  expect_equal(tail(s$u, 1)[[1L]], ua, tolerance = 1e-6)
  t <- t_nounits; D <- D_nounits
  xs <- lapply(1:30, function(i) Variable(paste0("nc_sp_x", i), t))
  sys <- complete(System(lapply(1:30, function(i) D(xs[[i]]) ~ -i * xs[[i]]), t, name = "nc_sp"))
  M <- generate_jacobian(sys, sparse = TRUE)(rep(1, 30), numeric(), 0)
  expect_s4_class(M, "dgCMatrix")
  expect_equal(Matrix::diag(M), -(1:30))
})

test_that("inline solves have implicit derivatives and generated code shares subexpressions", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("nc_im_x", t); z <- Variable("nc_im_z", t, guess = 1); w <- Variable("nc_im_w", t, guess = 1)
  q <- Parameter("nc_im_q")
  cs <- mtkcompile(System(list(D(x) ~ z * w - x, 0 ~ z^3 + z - q * x - w, 0 ~ w^3 - z), t, name = "nc_im"),
                   inline_nonlinear_sccs = TRUE)
  J <- generate_jacobian(cs); f <- generate_rhs(cs)
  expect_equal(J(0.7, 1.3, 0)[1, 1], (f(0.7 + 1e-6, 1.3, 0) - f(0.7 - 1e-6, 1.3, 0)) / 2e-6, tolerance = 1e-7)
  y <- Variable("nc_il_y", t); a <- Variable("nc_il_a", t); b <- Variable("nc_il_b", t); c3 <- Variable("nc_il_c", t)
  k <- Parameter("nc_il_k")
  cl <- mtkcompile(System(list(D(y) ~ a + b + c3, 0 ~ k * y * a + y^2 * b + c3 - 1, 0 ~ a - y * b + y^3 * c3 - 2,
                               0 ~ y^2 * a + b - y * c3), t, name = "nc_il"),
                   inline_linear_sccs = TRUE, analytical_linear_scc_limit = 0)
  J3 <- generate_jacobian(cl); f3 <- generate_rhs(cl)
  expect_equal(J3(0.4, 2, 0)[1, 1], (f3(0.4 + 1e-6, 2, 0) - f3(0.4 - 1e-6, 2, 0)) / 2e-6, tolerance = 1e-7)
  expect_false(any(grepl("fd_partial", deparse(body(J3)))))
  u <- Variable("nc_cs_x", t); v <- Variable("nc_cs_y", t); oa <- Variable("nc_cs_a", t); ob <- Variable("nc_cs_b", t)
  sc <- mtkcompile(System(list(D(u) ~ oa + sin(u * v)^2, D(v) ~ ob - sin(u * v)^2, oa ~ exp(sin(u * v)),
                               ob ~ oa * cos(u)), t, name = "nc_cs"))
  expect_equal(sum(grepl("sin(", deparse(body(generate_rhs(sc))), fixed = TRUE)), 1L)
})

test_that("loss derivatives, reduced storage, ensembles and caches", {
  prob <- ODEProblem(function(u, p, t) -p[1] * u, c(1, 2), c(0, 1), 2)
  loss <- function(sol) sum(tail(as.matrix(sol), 1)^2) + 3 * sol$prob$p[1]
  lg <- function(sol) {
    M <- as.matrix(sol); G <- M * 0; G[nrow(M), ] <- 2 * M[nrow(M), ]
    list(u = G, p = 3)
  }
  g <- solve_gradient(loss, prob, Tsit5(), saveat = 0.1, reltol = 1e-10, abstol = 1e-12, loss_grad = lg)
  expect_equal(unname(g$p), -10 * exp(-4) + 3, tolerance = 1e-7)
  expect_error(solve_gradient(loss, prob, Tsit5(), save_idxs = 1, save_idxs_storage = "reduced"), "every state")
  big <- ODEProblem(function(u, p, t) -u, as.numeric(1:200), c(0, 1))
  s <- solve(big, Tsit5(), saveat = 0.01, save_idxs = c(1, 5), save_idxs_storage = "reduced")
  expect_equal(dim(.sol_env(s)$U), c(101L, 2L))
  expect_equal(tail(s$u, 1)[[1L]], c(1, 5) * exp(-1), tolerance = 1e-6)
  ep <- EnsembleProblem(SDEProblem(function(u, p, t) -u, function(u, p, t) 0.5, 1, c(0, 1)),
                        output_func = function(sol, i) list(tail(sol$u, 1)[[1L]], FALSE))
  a <- solve(ep, EM(), dt = 0.01, trajectories = 4, seed = 4, rng_streams = "trajectory")
  b <- solve(ep, EM(), EnsembleThreads(2), dt = 0.01, trajectories = 4, seed = 4, batch_size = 2,
             rng_streams = "trajectory")
  expect_identical(unlist(a$u), unlist(b$u))
  expect_true(a$converged); expect_false(b$converged)
  f <- tempfile(fileext = ".rds")
  full <- solve(ep, EM(), dt = 0.01, trajectories = 4, seed = 4, batch_size = 2, rng_streams = "trajectory", checkpoint = f)
  ck <- readRDS(f); ck$done <- 2L; ck$u <- ck$u[1:2]; ck$retcodes[3:4] <- ""; saveRDS(ck, f)
  res <- solve(ep, EM(), dt = 0.01, trajectories = 4, seed = 4, batch_size = 2, rng_streams = "trajectory", checkpoint = f)
  expect_identical(unlist(res$u), unlist(full$u))
  flaky <- EnsembleProblem(ODEProblem(function(u, p, t) -u, 1, c(0, 1)),
                           prob_func = function(prob, i, r) { if (i == 2L && r == 1L) stop("flaky"); prob })
  expect_equal(solve(flaky, Tsit5(), trajectories = 3, retries = 1)$failures, 0L)
  expect_true(all(c("sensitivity_systems", "registered_functions") %in% names(cache_info())))
  expect_equal(.simplex_compositions(3, 2), unname(as.matrix(subset(expand.grid(0:3, 0:3), Var1 + Var2 <= 3))),
               ignore_attr = TRUE)
})
