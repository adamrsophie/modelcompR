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
  # a declared terminal loss is differentiated in its row alone
  calls <- 0L
  old <- .mtk$sens_loss_rows
  .mtk$sens_loss_rows <- "last"
  ct <- .sens_loss_cotangent(function(s) { calls <<- calls + 1L; sum(tail(as.matrix(s), 1)^2) }, sol)
  .mtk$sens_loss_rows <- old
  expect_lt(calls, 30L)
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

test_that("loss gradients keep every coordinate and the public solution view", {
  # a weighted cancellation that a probe of the rows would miss (65 rows)
  prob <- ODEProblem(function(u, p, t) p[1] + 0 * u, 0, c(0, 1), p = 0)
  h0 <- .Machine$double.eps^(1 / 3)
  hm <- h0 * (1 + 0.5 * sin(1:65 * 0.7548776662))
  loss <- function(sol) { U <- as.matrix(sol); (hm[2] * U[1, 1] - hm[1] * U[2, 1]) / h0 }
  for (ns in c(64, 65)) {
    tt <- seq(0, 1, length.out = ns)
    for (sa in list(ForwardSensitivity(), InterpolatingAdjoint())) {
      g <- solve_gradient(loss, prob, Tsit5(), saveat = tt, sensealg = sa, reltol = 1e-10, abstol = 1e-12)
      expect_equal(unname(g$p), -hm[1] / h0 / (ns - 1), tolerance = 1e-8)
      expect_equal(unname(g$u0), (hm[2] - hm[1]) / h0, tolerance = 1e-8)
    }
  }
  # loss_rows declares a terminal loss; a wrong declaration is rejected
  p2 <- ODEProblem(function(u, p, t) -p[1] * u, c(1, 2), c(0, 1), p = 1)
  tl <- function(sol) sum(tail(as.matrix(sol), 1))
  a <- solve_gradient(tl, p2, Tsit5(), saveat = 0.01, reltol = 1e-10, abstol = 1e-12)
  b <- solve_gradient(tl, p2, Tsit5(), saveat = 0.01, reltol = 1e-10, abstol = 1e-12, loss_rows = "last")
  expect_equal(b$p, a$p, tolerance = 1e-8)
  expect_error(solve_gradient(tl, p2, Tsit5(), saveat = 0.5, loss_rows = 9), "loss_rows")
  # save_idxs: the loss and loss_grad see the saved component, the gradient
  # is that of the visible objective
  lg <- function(sol) { U <- as.matrix(sol); U[] <- 0; U[nrow(U), ] <- 1; U }
  for (sa in list(ForwardSensitivity(), InterpolatingAdjoint())) for (use_lg in c(FALSE, TRUE)) {
    g <- solve_gradient(tl, p2, Tsit5(), save_idxs = 2, saveat = c(0, 1), sensealg = sa, reltol = 1e-10,
                        abstol = 1e-12, loss_grad = if (use_lg) lg)
    expect_equal(dim(as.matrix(g$sol)), c(2L, 1L))
    expect_equal(g$value, 2 * exp(-1), tolerance = 1e-8)
    expect_equal(unname(g$p), -2 * exp(-1), tolerance = 1e-7)
    expect_equal(unname(g$u0), c(0, exp(-1)), tolerance = 1e-7)
  }
  # reordered selection of two outputs
  g <- solve_gradient(function(sol) { U <- as.matrix(sol); U[2, 1] + 10 * U[2, 2] }, p2, Tsit5(),
                      save_idxs = c(2, 1), saveat = c(0, 1), reltol = 1e-10, abstol = 1e-12)
  expect_equal(unname(g$u0), c(10, 1) * exp(-1), tolerance = 1e-7)
  # reduced storage given through the problem is rejected too
  pr <- ODEProblem(function(u, p, t) -p[1] * u, c(1, 2), c(0, 1), p = 1, save_idxs = 2, save_idxs_storage = "reduced")
  expect_error(solve_gradient(tl, pr, Tsit5(), saveat = c(0, 1)), "every state")
  # loss_grad still differences the samples read by interpolation
  p3 <- ODEProblem(function(u, p, t) p[1] + 0 * u, 0, c(0, 1), p = 2)
  g <- solve_gradient(function(sol) sum(sol(0.5)), p3, Tsit5(), saveat = c(0, 1), sensealg = ForwardSensitivity(),
                      loss_grad = function(sol) 0 * as.matrix(sol))
  expect_equal(c(unname(g$p), unname(g$u0)), c(0.5, 1), tolerance = 1e-8)
})

test_that("every numerical stage of a gradient is checked", {
  prob <- ODEProblem(function(u, p, t) p[1] * u, 1, c(0, 1), p = 20)
  loss <- function(sol) tail(as.matrix(sol)[, 1], 1)
  e <- tryCatch(suppressWarnings(solve_gradient(loss, prob, Tsit5(), saveat = c(0, 1), sensealg = InterpolatingAdjoint(),
                                                dt = 0.1, maxiters = 5, reltol = 1e-9, abstol = 1e-11)),
                error = function(e) e)
  # the forward solve already stops at maxiters, before the backward one
  expect_s3_class(e, "sens_stage_failure")
  expect_equal(e$stage, "forward")
  expect_equal(e$retcode, "MaxIters")
  expect_lt(e$t, 1)
  g <- solve_gradient(loss, prob, Tsit5(), saveat = c(0, 1), sensealg = InterpolatingAdjoint(), reltol = 1e-9,
                      abstol = 1e-11)
  expect_equal(unname(g$u0) / exp(20), 1, tolerance = 1e-6)
  # a failure while replaying a checkpoint
  bad <- new.env(); bad$on <- FALSE
  pc <- ODEProblem(function(u, p, t) if (bad$on && t < 0.5) NaN * u else -p[1] * u, 1, c(0, 1), p = 1)
  sa <- InterpolatingAdjoint(checkpointing = TRUE)
  sol <- solve(pc, Tsit5(), saveat = c(0, 0.5, 1), sensealg = sa)
  bad$on <- TRUE
  e2 <- tryCatch(suppressWarnings(adjoint_sensitivities(sol, Tsit5(), t = 1, sensealg = sa,
                                                        dgdu_discrete = function(u, p, t, i) 1)), error = function(e) e)
  expect_s3_class(e2, "sens_stage_failure")
  expect_equal(e2$stage, "checkpoint replay")
})

test_that("reduced storage maps original states to stored columns everywhere", {
  prob <- ODEProblem(function(u, p, t) c(0, 0, 0), c(10, 20, 30), c(0, 1))
  s <- solve(prob, Tsit5(), saveat = c(0, 1), save_idxs = c(3, 2), save_idxs_storage = "reduced")
  expect_equal(s[2], c(20, 20))
  expect_equal(s(0.5, idxs = 2), 20)
  expect_equal(s(0.5, idxs = 3), 30)
  expect_equal(unname(s(0.5, idxs = list(2, 3))), c(20, 30))
  expect_error(s[1], "not saved")
  expect_error(s(0.5, idxs = 1), "not saved")
})

test_that("sparse prototypes, masks and Matrix storage give the same Jacobian", {
  n <- 3L
  A <- matrix(c(-100, 0, 1, 0, -2, 0, 0, 3, -1), 3)
  f <- function(u, p, t) as.numeric(A %*% u)
  alg <- Rodas5P(autodiff = AutoSparse(AutoFiniteDiff()))
  protos <- list(Matrix::sparseMatrix(i = c(1, 2, 3, 2, 3), j = c(1, 2, 1, 3, 3), x = 0, dims = c(3, 3)),
                 Matrix::Matrix(A, sparse = TRUE), A * 0)
  ref <- solve(ODEProblem(f, rep(1, n), c(0, 0.1)), Rodas5P(), dt = 0.02, adaptive = FALSE)
  for (pr in protos) {
    sp <- ODEFunction(f, jac_prototype = pr)
    ctx <- list(f = f, p = numeric(), autodiff = alg$autodiff, sparsity = sp$sparsity)
    J <- .calc_J(ctx, rep(1, n), 0, f(rep(1, n), numeric(), 0), sparse_ok = TRUE)$J
    expect_equal(as.matrix(J), A, tolerance = 1e-6, ignore_attr = TRUE)
    s <- solve(ODEProblem(sp, rep(1, n), c(0, 0.1)), alg, dt = 0.02, adaptive = FALSE)
    expect_equal(tail(as.matrix(s), 1), tail(as.matrix(ref), 1), tolerance = 1e-6, ignore_attr = TRUE)
  }
  # mirrored and implicit unit diagonal entries take part in the coloring
  valid <- function(M, cl) { D <- as.matrix(M) != 0
    all(vapply(seq_len(nrow(D)), function(i) !anyDuplicated(cl[D[i, ]]), TRUE)) }
  S <- Matrix::bandSparse(5, k = c(0, 1), diagonals = list(rep(2, 5), rep(1, 4)), symmetric = TRUE)
  Tu <- methods::as(Matrix::bandSparse(4, k = 1, diagonals = list(rep(1, 3))), "dtCMatrix"); Tu@diag <- "U"
  for (M in list(S, Tu, Matrix::t(S), methods::as(S, "generalMatrix"), as.matrix(S) != 0))
    expect_true(valid(M, column_coloring(M)))
  fs <- function(u) as.numeric(S %*% u)
  expect_equal(as.matrix(.jac_colored(fs, rep(1, 5), fs(rep(1, 5)), S)$J), as.matrix(S), tolerance = 1e-6,
               ignore_attr = TRUE)
})

test_that("native arithmetic keeps complex states in the implicit methods", {
  p <- ODEProblem(ODEFunction(function(u, p, t) 1i * u, jac_fn = function(u, p, t) matrix(1i, 1, 1),
                              tgrad_fn = function(u, p, t) 0i * u), 1 + 0i, c(0, 1))
  old <- arithmetic_mode(); on.exit(arithmetic_mode(old))
  for (alg in list(Rodas5P(), Rosenbrock23(), RadauIIA5(), ImplicitEuler(), FBDF(), KenCarp4())) {
    arithmetic_mode("reference")
    r <- tail(as.matrix(solve(p, alg, reltol = 1e-8, abstol = 1e-10)), 1)
    arithmetic_mode("native")
    nv <- tail(as.matrix(solve(p, alg, reltol = 1e-8, abstol = 1e-10)), 1)
    expect_true(is.complex(nv))
    expect_equal(nv, r, tolerance = 1e-6, ignore_attr = TRUE)
  }
})

test_that("inline nonlinear branches commit only accepted states and replay with the solution", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("bc_x", t); y <- Variable("bc_y", t); z <- Variable("bc_z", t, guess = 1)
  cs <- mtkcompile(System(list(D(x) ~ -100 * x, D(y) ~ z, 0 ~ z^3 - z - x), t, name = "bc_rej"),
                   inline_nonlinear_sccs = TRUE)
  p <- ODEProblem(cs, list(x ~ 0.1, y ~ 0), c(0, 0.1))
  a <- solve(p, BS3(), dt = 0.08, reltol = 1e-6, abstol = 1e-8, maxiters = 3000)
  b <- solve(p, BS3(), dt = 1e-5, reltol = 1e-6, abstol = 1e-8, maxiters = 3000)
  expect_gt(a$stats$nreject, 0L)
  expect_gt(a(0.1, idxs = y), 0.1)
  expect_equal(a(0.1, idxs = y), b(0.1, idxs = y), tolerance = 1e-5)
  # the branch of z^3 - z = t - 1 that starts at the only root at t = 0
  w <- Variable("bc_w", t); v <- Variable("bc_v", t, guess = 1); q <- Parameter("bc_q")
  cs2 <- mtkcompile(System(list(D(w) ~ v, 0 ~ v^3 - v - (t - 1 + q)), t, name = "bc_live"), inline_nonlinear_sccs = TRUE)
  p2 <- ODEProblem(cs2, list(w ~ 0, q ~ 0), c(0, 1))
  i <- init(p2, Tsit5(), dt = 0.01, adaptive = FALSE)
  s <- solve_integrator(i)
  expect_equal(s(1, idxs = v), -1, tolerance = 1e-8)
  expect_equal(unname(s[v][length(s$t)]), -1, tolerance = 1e-8)
  expect_equal(get_du(i), -1, tolerance = 1e-8)
  # reinit starts over from the guesses, as a fresh solve
  reinit(i, u0 = 0, t0 = 0.9, tf = 1)
  r <- solve_integrator(i)
  f <- solve(remake(p2, u0 = 0, tspan = c(0.9, 1)), Tsit5(), dt = 0.01, adaptive = FALSE)
  expect_equal(r(1, idxs = w), f(1, idxs = w), tolerance = 1e-10)
  # adjoints replay the forward branch: d w(1) / d q = -1 - z(0)
  truth <- -1 - uniroot(function(z) z^3 - z + 1, c(-2, -1), tol = 1e-13)$root
  loss <- function(s) tail(as.matrix(s), 1)[1]
  for (sa in list(ForwardSensitivity(), InterpolatingAdjoint(), InterpolatingAdjoint(checkpointing = TRUE))) {
    g <- solve_gradient(loss, p2, Tsit5(), sensealg = sa, dt = 0.01, adaptive = FALSE, saveat = 0.1,
                        reltol = 1e-8, abstol = 1e-10)
    expect_equal(unname(g$p), truth, tolerance = 1e-5)
  }
})

test_that("a problem keeps the registrations it was built with", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("rs_x", t)
  h <- register_symbolic("rs_h", function(x) x, derivatives = list(function(x) 1))
  cs <- mtkcompile(System(list(D(x) ~ h(x)), t, name = "rs_sys"))
  p <- ODEProblem(cs, list(x ~ 1), c(0, 1))
  loss <- function(s) tail(as.matrix(s), 1)[1]
  h <- register_symbolic("rs_h", function(x) x^2, derivatives = list(function(x) 2 * x))
  g <- solve_gradient(loss, p, Tsit5(), sensealg = ForwardSensitivity(), saveat = 1, reltol = 1e-9, abstol = 1e-11)
  expect_equal(c(g$value, unname(g$u0)), rep(exp(1), 2), tolerance = 1e-6)
  p2 <- ODEProblem(mtkcompile(System(list(D(x) ~ h(x)), t, name = "rs_sys2")), list(x ~ 0.5), c(0, 1))
  g2 <- solve_gradient(loss, p2, Tsit5(), sensealg = ForwardSensitivity(), saveat = 1, reltol = 1e-10, abstol = 1e-12)
  expect_equal(c(g2$value, unname(g2$u0)), c(1, 4), tolerance = 1e-6)
})

test_that("changing and nonfinite mass matrices", {
  calls <- 0L
  M <- MatrixOperator(matrix(1), update_func = function(A, u, p, t) { calls <<- calls + 1L; matrix(1 + t) })
  prob <- ODEProblem(ODEFunction(function(u, p, t) -u, mass_matrix = M), 1, c(0, 1))
  for (alg in list(Rodas5P(), Tsit5(), FBDF())) {
    s <- solve(prob, alg, reltol = 1e-10, abstol = 1e-12)
    expect_equal(tail(as.matrix(s)[, 1], 1), 0.5, tolerance = 1e-7)
  }
  expect_gt(calls, 0L)
  Ms <- MatrixOperator(matrix(c(1, 0, 0, 0), 2), update_func = function(A, u, p, t) A)
  expect_error(ODEProblem(ODEFunction(function(u, p, t) -u, mass_matrix = Ms), c(1, 1), c(0, 1)), "nonsingular")
  for (v in c(NA_real_, NaN, Inf))
    expect_error(ODEProblem(ODEFunction(function(u, p, t) -u, mass_matrix = matrix(v)), 1, c(0, 1)), "non-finite")
})

test_that("Krylov solves, low rank updates and eigenpairs keep their contracts", {
  for (alg in list(KrylovJL_GMRES(), KrylovJL_CG(), KrylovJL_BICGSTAB(), KrylovJL_MINRES(), KrylovJL_FGMRES()))
    for (sc in c(1, 1e-200, 1e200, 1e-300, 1e300)) {
      A <- matrix(c(4, 1, 0, 1, 3, 1, 0, 1, 2), 3)
      s <- solve(LinearProblem(A, as.numeric(A %*% c(1, -2, 3)) * sc), alg, abstol = 0, reltol = 1e-10)
      expect_equal(s$retcode, "Success")
      expect_equal(s$u / sc, c(1, -2, 3), tolerance = 1e-8)
    }
  # a scaled left preconditioner does not stop the inner iteration early
  A <- diag(c(1, 2, 4)); b <- rep(1, 3)
  for (sc in c(1, 1e-12, 1e12)) {
    z <- solve(LinearProblem(A, b), KrylovJL_GMRES(), Pl = function(x) sc * x, maxiters = 3, restart = 3, reltol = 1e-10)
    expect_equal(z$retcode, "Success")
    expect_equal(z$u, c(1, 0.5, 0.25), tolerance = 1e-8)
  }
  # an ill conditioned base with a well conditioned update
  A <- matrix(c(1, 1, 1, 1 + 1e-8), 2); U <- matrix(c(0, 1), 2, 1)
  M <- A + U %*% t(U)
  s <- solve(LinearProblem(LowRankUpdatedMatrix(A, U, matrix(1), U), as.numeric(M %*% c(1, 1))), LUFactorization())
  expect_equal(s$u, c(1, 1), tolerance = 1e-8)
  s <- solve(LinearProblem(LowRankUpdatedMatrix(diag(3), matrix(1, 3, 1), matrix(0), matrix(1, 3, 1)), c(1, 2, 3)))
  expect_equal(s$u, c(1, 2, 3))
  # eigenpairs: repeated eigenvalues, the complete count, impossible requests
  set.seed(1)
  e <- solve(EigenvalueProblem(diag(10), nev = 2), ArnoldiMethodJL(krylovdim = 2))
  expect_equal(e$retcode, "Success"); expect_length(e$values, 2L); expect_equal(qr(e$vectors)$rank, 2L)
  set.seed(2)
  e <- solve(EigenvalueProblem(diag(c(5, 5, 5, 4, 3, 2, 1, 1, 1, 1, 0.5, 0.2)), nev = 4), ArnoldiMethodJL(krylovdim = 5))
  expect_equal(e$values, c(5, 5, 5, 4), tolerance = 1e-8)
  set.seed(3)
  B <- matrix(rnorm(400), 20)
  e <- solve(EigenvalueProblem(B, nev = 3), ArnoldiMethodJL(krylovdim = 6))
  expect_equal(e$retcode, "Success")
  expect_equal(e$values, eigen(B)$values[1:3], tolerance = 1e-8)
  e <- solve(EigenvalueProblem(diag(3), nev = 5))
  expect_equal(c(e$retcode, e$nconverged), c("Failure", "3"))
})

test_that("regression sensitivity is invariant to shifts, units, order and ties", {
  set.seed(5)
  X <- matrix(runif(60), 3); Y <- rbind(X[1, ] + 2 * X[2, ]^2, X[3, ] - X[1, ])
  a <- gsa(X, Y, RegressionGSA(rank = TRUE))
  b <- gsa(1e-6 * X + 1e3, 1e6 * Y - 3, RegressionGSA(rank = TRUE))
  for (k in names(a)) expect_equal(b[[k]], a[[k]], tolerance = 1e-6, info = k)
  b2 <- gsa(X + 1e8, Y, RegressionGSA())
  for (k in c("pearson", "standard_regression", "partial_correlation")) expect_equal(b2[[k]], a[[k]], tolerance = 1e-6, info = k)
  ord <- sample(20)
  d <- gsa(X[, ord], Y[, ord], RegressionGSA(rank = TRUE))
  for (k in names(a)) expect_equal(d[[k]], a[[k]], tolerance = 1e-12, info = k)
  Xt <- matrix(c(1, 1, 2, 2, 3), 1); Yt <- c(2, 1, 3, 5, 4)
  expect_equal(c(gsa(Xt, Yt, RegressionGSA(rank = TRUE))$pearson_rank),
               stats::cor(c(Xt), Yt, method = "spearman"))
  cst <- gsa(rbind(rep(1, 5), 1:5), (1:5)^2, RegressionGSA(rank = TRUE))
  expect_true(all(is.na(c(cst$pearson[1], cst$standard_regression[1], cst$partial_correlation[1],
                          cst$pearson_rank[1]))))
  expect_true(all(is.finite(c(cst$pearson[2], cst$standard_regression[2]))))
  # collinear inputs: undetermined coefficients are NA
  col <- gsa(rbind(1:6, 2 * (1:6), c(3, 1, 4, 1, 5, 9)), c(1, 4, 2, 8, 5, 7), RegressionGSA())
  expect_true(anyNA(col$standard_regression[1:2]))
})

test_that("Shapley handles one input and keeps its design apart from the model", {
  D <- SklarDist(IndependentCopula(1), list(univariate_distribution("unif", 0, 1)))
  m <- Shapley(n_var = 10, n_outer = 3, n_inner = 2)
  set.seed(1); b <- suppressMessages(gsa(function(X) X[1, ], m, D, batch = TRUE))
  set.seed(1); a <- suppressMessages(gsa(function(x) x[1], m, D))
  expect_equal(a$shapley_effects, 1); expect_equal(b$shapley_effects, 1)
  D3 <- SklarDist(IndependentCopula(3), rep(list(univariate_distribution("unif", 0, 1)), 3))
  m3 <- Shapley(n_perms = 4, n_var = 50, n_outer = 5, n_inner = 3)
  set.seed(3); b3 <- suppressMessages(gsa(function(X) X[1, ] + 2 * X[2, ] * X[3, ], m3, D3, batch = TRUE))
  set.seed(3); a3 <- suppressMessages(gsa(function(x) { stats::runif(1); x[1] + 2 * x[2] * x[3] }, m3, D3))
  expect_equal(a3$shapley_effects, b3$shapley_effects)
})

test_that("the end state obeys the event side policy whether or not it is saved", {
  prob <- ODEProblem(function(u, p, t) 0 * u, 1, c(0, 1))
  for (sp in list(c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, FALSE))) for (se in c(TRUE, FALSE)) {
    cb <- PresetTimeCallback(1, function(i) i$u[] <- 10, save_positions = sp)
    s <- solve(prob, Tsit5(), callback = cb, dense = TRUE, saveat = c(0, 0.5), save_end = se)
    expect_equal(s(1), 10)
    expect_equal(s(1 - 1e-9), 1)
  }
})

test_that("ensemble checkpoints are versioned, identified and resume exactly", {
  base <- ODEProblem(function(u, p, t) -u, 1, c(0, 0.1))
  flag <- new.env(); flag$abort <- FALSE
  red <- function(u, batch, I) { if (flag$abort && max(I) > 2) stop("interrupted"); list(c(u, batch), FALSE) }
  ep <- EnsembleProblem(base, output_func = function(sol, i) list(stats::runif(1), FALSE), reduction = red)
  for (bs in 1:3) for (streams in c("trajectory", "session")) for (sd in list(NULL, 7L)) {
    f <- tempfile(fileext = ".rds")
    set.seed(42); full <- solve(ep, Tsit5(), trajectories = 6, batch_size = bs, rng_streams = streams, seed = sd)
    set.seed(42); flag$abort <- TRUE
    expect_error(solve(ep, Tsit5(), trajectories = 6, batch_size = bs, rng_streams = streams, seed = sd, checkpoint = f),
                 "interrupted")
    flag$abort <- FALSE
    set.seed(999)
    res <- solve(ep, Tsit5(), trajectories = 6, batch_size = bs, rng_streams = streams, seed = sd, checkpoint = f)
    expect_identical(res$u, full$u)
    unlink(f)
  }
  f <- tempfile(fileext = ".rds")
  ep2 <- function(rate) EnsembleProblem(ODEProblem(function(u, p, t) -p * u, 1, c(0, 1), rate))
  solve(ep2(1), Tsit5(), trajectories = 4, batch_size = 2, checkpoint = f)
  expect_error(solve(ep2(2), Tsit5(), trajectories = 4, batch_size = 2, checkpoint = f), "different run")
  expect_error(solve(ep2(1), Tsit5(), trajectories = 4, batch_size = 1, checkpoint = f), "batch_size")
  expect_error(solve(ep2(1), Tsit5(), trajectories = 4, batch_size = 2, checkpoint = f, run_id = "a"), "run_id")
  writeBin(readBin(f, "raw", 20), f)
  expect_error(solve(ep2(1), Tsit5(), trajectories = 4, batch_size = 2, checkpoint = f), "truncated")
  saveRDS(list(a = 1), f)
  expect_error(solve(ep2(1), Tsit5(), trajectories = 4, batch_size = 2, checkpoint = f), "not an ensemble checkpoint")
})

test_that("dual means of large finite values stay finite", {
  r <- .ad_try_jacobian(function(x) mean(x), rep(1e308, 2))
  expect_equal(r$value, 1e308)
  expect_equal(c(r$jacobian), c(0.5, 0.5))
})

# ---------------------------------------------------------------------------
# Third numerical contract review

test_that("a singular Rosenbrock stage on a fixed step fails in both arithmetic modes", {
  J <- matrix(1, 2, 2) + 4 * diag(2)
  pr <- ODEProblem(ODEFunction(function(u, p, t) drop(J %*% u) + c(1, 0), jac_fn = function(u, p, t) J,
                               tgrad_fn = function(u, p, t) c(0, 0)), c(0, 0), c(0, 1))
  old <- arithmetic_mode()
  on.exit(arithmetic_mode(old), add = TRUE)
  for (mode in c("reference", "native")) {
    arithmetic_mode(mode)
    s <- suppressWarnings(solve(pr, Rodas4(), dt = 1, adaptive = FALSE))
    expect_identical(s$retcode, "ConvergenceFailure")
  }
})

test_that("inline roots take the same branch at any scale of the residual", {
  t <- t_nounits
  w <- Variable("nc_inl_w", t)
  z <- Variable("nc_inl_z", t, guess = 1)
  q <- Parameter("nc_inl_q")
  sys <- mtkcompile(System(list(D_nounits(w) ~ z, Equation(0, q * (z^2 - 4))), t, name = "nc_inline"),
                    inline_nonlinear_sccs = TRUE)
  for (a in c(1, 1e-16)) {
    s <- solve(ODEProblem(sys, list(w ~ 0, q ~ a), c(0, 1)), Tsit5())
    expect_equal(unname(c(s(1, idxs = w), s(1, idxs = z))), c(2, 2), tolerance = 1e-6)
  }
})

test_that("checkpoints compare exact descriptions, record the backend and survive a failed replace", {
  span <- c(0, 1)
  ep <- function(k) EnsembleProblem(ODEProblem(function(u, p, t) p[1], 0, span, k),
                                    output_func = function(s, i) list(s(1), FALSE))
  f <- tempfile(fileext = ".rds")
  solve(ep(c(1L, 3L, 2L)), Tsit5(), trajectories = 1L, checkpoint = f)
  expect_error(solve(ep(c(2L, 1L, 3L)), Tsit5(), trajectories = 1L, checkpoint = f), "different run")
  ck <- readRDS(f)
  expect_identical(ck$version, 2L)
  p <- ODEProblem(function(u, p, t) 0, 0, c(0, 0.1))
  e2 <- EnsembleProblem(p, output_func = function(s, i) list(stats::runif(1), FALSE))
  f2 <- tempfile(fileext = ".rds")
  solve(e2, Tsit5(), EnsembleSerial(), trajectories = 2, seed = 42, checkpoint = f2)
  expect_error(solve(e2, Tsit5(), EnsembleThreads(cores = 1), trajectories = 2, seed = 42, checkpoint = f2),
               "ensemblealg")
  # an interrupted replacement leaves the old checkpoint as .bak, which is read
  file.rename(f2, paste0(f2, ".bak"))
  expect_silent(suppressMessages(solve(e2, Tsit5(), EnsembleSerial(), trajectories = 2, seed = 42, checkpoint = f2)))
  old <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
  RNGkind("Mersenne-Twister", "Box-Muller", "Rejection")
  flag <- new.env(); flag$stop <- FALSE
  e3 <- EnsembleProblem(p, output_func = function(s, i) list(stats::rnorm(1), FALSE),
                        reduction = function(u, batch, I) { if (flag$stop && I == 2) stop("interrupted"); list(c(u, batch), FALSE) })
  run <- function(ck = NULL) solve(e3, Tsit5(), trajectories = 2, batch_size = 1, seed = 42, checkpoint = ck)
  full <- run()
  f3 <- tempfile(fileext = ".rds")
  flag$stop <- TRUE; expect_error(run(f3), "interrupted"); flag$stop <- FALSE
  expect_identical(run(f3)$u, full$u)
})

test_that("Shapley designs do not share the Box-Muller cache with the model", {
  old <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
  RNGkind("Mersenne-Twister", "Box-Muller", "Rejection")
  D <- SklarDist(GaussianCopula(matrix(c(1, .4, .4, 1), 2)), rep(list(univariate_distribution("unif", 0, 1)), 2))
  m <- Shapley(n_var = 11, n_outer = 3, n_inner = 3)
  set.seed(1); a <- suppressMessages(gsa(function(X) X[1, ] + 2 * X[2, ], m, D, batch = TRUE))
  set.seed(1); b <- suppressMessages(gsa(function(x) { stats::rnorm(1); x[1] + 2 * x[2] }, m, D))
  expect_identical(a$shapley_effects, b$shapley_effects)
})

test_that("adjoints fail on unfinished quadrature, forward or backward solves", {
  span <- c(0, 1)
  last1 <- function(s) tail(as.matrix(s)[, 1], 1)
  h <- 2^-60; k <- function(t) (t + h)^(-0.99)
  pr <- ODEProblem(ODEFunction(function(u, p, t) p[1] * k(t), jac_fn = function(u, p, t) matrix(0),
                               paramjac = function(u, p, t) matrix(k(t))), 0, span, 0)
  g <- solve_gradient(last1, pr, Tsit5(), saveat = span, sensealg = QuadratureAdjoint(), abstol = 1e-8, reltol = 1e-8)
  expect_equal(unname(g$p), ((1 + h)^0.01 - h^0.01) / 0.01, tolerance = 1e-6)
  pr <- ODEProblem(function(u, p, t) p[1] * u, 1, span, 1)
  s <- suppressWarnings(solve(pr, Tsit5(), dt = 0.1, adaptive = FALSE, maxiters = 1, dense = TRUE))
  expect_error(adjoint_sensitivities(s, Tsit5(), t = 1, dgdu_discrete = function(u, p, t, i) 1),
               class = "sens_stage_failure")
  s <- solve(ODEProblem(function(u, p, t) u, 1, span), Tsit5(), reltol = 1e-10, abstol = 1e-12)
  cb <- PresetTimeCallback(.5, function(i) terminate(i))
  expect_error(adjoint_sensitivities(s, Tsit5(), t = 1, callback = cb, dgdu_discrete = function(u, p, t, i) 1),
               "before t = 0")
  # the checkpoint replay takes the fixed step of the forward solve
  pr <- ODEProblem(function(u, p, t) -u, 1, span)
  g <- solve_gradient(last1, pr, Euler(), dt = .1, saveat = span, sensealg = InterpolatingAdjoint(checkpointing = TRUE))
  expect_equal(unname(g$u0), 0.9^10, tolerance = 1e-10)
})

test_that("the loss cotangent covers every state, repeated columns and interpolated reads", {
  span <- c(0, 1)
  flat <- ODEProblem(function(u, p, t) 0 * u, 2, span)
  two <- ODEProblem(function(u, p, t) 0 * u, c(2, 3), span)
  G <- function(loss, pr = flat, ...) solve_gradient(loss, pr, Tsit5(), saveat = span, ...)
  expect_equal(unname(G(function(s) tail(s[1], 1), two, save_idxs = 2)$u0), c(1, 0))
  expect_equal(unname(G(function(s) sum(tail(as.matrix(s), 1)), two, save_idxs = c(1, 1),
                        loss_grad = function(s) rbind(c(0, 0), c(1, 1)))$u0), c(2, 0))
  expect_equal(unname(G(function(s) as.numeric(s(1)), loss_grad = function(s) 0 * as.matrix(s))$u0), 1)
  expect_equal(unname(G(function(s) s(0), loss_rows = "last")$u0), 1)
  g <- G(function(s) s(.1) + s(.2) + s(.3) + s(.4))
  expect_equal(unname(c(g$value, g$u0)), c(8, 4))
  run <- function(T0) {
    pr <- ODEProblem(function(u, p, t) 1, 0, T0 + c(0, 2))
    solve_gradient(function(s) as.numeric(s(T0 + 0.5)), pr, Tsit5(), dt = 0.25, adaptive = FALSE,
                   saveat = T0 + c(0, .5, 1, 1.5, 2))$value
  }
  expect_equal(unname(run(1e12)), 0.5)
  p <- ODEProblem(function(u, p, t) p, c(1, 2), span, c(3, 4))
  for (sa in list(ForwardSensitivity(), InterpolatingAdjoint()))
    expect_equal(unname(solve_gradient(function(s) tail(s[1], 1), p, Tsit5(), saveat = span, save_idxs = 2,
                                       sensealg = sa)$p), c(1, 0), tolerance = 1e-6)
})

test_that("optimization gradients follow the objective and the sense", {
  pr <- OptimizationProblem(function(x, p) (x[1] - p[1])^2 + p[1], 0, 2)
  g <- solve_gradient(function(s) s$objective, pr, BFGS())
  expect_equal(unname(c(g$value, g$p)), c(2, 1), tolerance = 1e-6)
  p <- OptimizationProblem(function(x, p) -(x[1] - p[1])^2, .5, -1, lb = 0, ub = Inf, sense = MaxSense)
  g <- solve_gradient(function(s) s$u[1], p, LBFGSB())
  expect_equal(unname(g$p), 0, tolerance = 1e-6)
})

test_that("sol(t) gives the right side of an event in dense and linear output", {
  flat <- ODEProblem(function(u, p, t) 0 * u, 2, c(0, 1))
  cb <- PresetTimeCallback(0.5, function(i) i$u <- 12)
  for (h in c(.5, .25)) for (d in c(TRUE, FALSE))
    expect_equal(unname(solve(flat, RK4(), dt = h, dense = d, callback = cb)(.5)), 12)
  cb <- PresetTimeCallback(1, function(i) i$u[] <- 10, save_positions = c(TRUE, FALSE))
  for (d in c(TRUE, FALSE)) expect_equal(unname(solve(flat, Tsit5(), callback = cb, dense = d)(1)), 10)
})

test_that("Krylov least squares and rescaling report their real status", {
  s <- solve(LinearProblem(matrix(1e-12), 1), KrylovJL_LSMR())
  expect_equal(s$u, 1e12, tolerance = 1e-10)
  expect_identical(s$retcode, "Success")
  z <- solve(LinearProblem(matrix(0.5), 1e308), KrylovJL_GMRES())
  expect_identical(z$retcode, "Failure")
  A <- matrix(c(1, 1, 1, 1, 2, 3), 3)
  expect_equal(solve(LinearProblem(A, c(1, 2, 2)), KrylovJL_LSMR())$u, qr.solve(A, c(1, 2, 2)), tolerance = 1e-8)
})

test_that("changing mass matrices survive remake and split problems", {
  f <- function(u, p, t) -u
  M <- MatrixOperator(matrix(1), function(A, u, p, t) matrix(1 + t))
  p <- remake(ODEProblem(f, 1, c(0, 1)), f = ODEFunction(f, mass_matrix = M))
  q <- SplitODEProblem(f, function(u, p, t) 0 * u, u0 = 1, tspan = c(0, 1), mass_matrix = M)
  expect_equal(solve(p, Tsit5(), reltol = 1e-10, abstol = 1e-12)(1), 0.5, tolerance = 1e-7)
  expect_equal(solve(q, KenCarp4(), reltol = 1e-10, abstol = 1e-12)(1), 0.5, tolerance = 1e-6)
  g <- ODEFunction(function(u, p, t) 1, mass_matrix = M, tgrad = function(u, p, t) 0)
  expect_equal(solve(ODEProblem(g, 0, c(0, 1)), Rosenbrock23(), reltol = 1e-10, abstol = 1e-12)(1), log(2),
               tolerance = 1e-6)
  M0 <- MatrixOperator(matrix(0), function(A, u, p, t) matrix(2))
  expect_equal(solve(ODEProblem(ODEFunction(function(u, p, t) 1, mass_matrix = M0), 0, c(0, 1)), Tsit5())(1), 0.5)
})

test_that("mass matrices of pattern and scalar operator type", {
  P <- Matrix::sparseMatrix(i = c(1, 2, 1), j = c(1, 2, 2), dims = c(2, 2))
  run <- function(M) solve(ODEProblem(function(u, p, t) -u, c(1, 2), c(0, 1), mass_matrix = M), Rodas5P())(1)
  expect_equal(run(P), run(as.matrix(P) + 0))
  expect_equal(run(ScalarOperator(2)), run(2))
  expect_true(.is_identity_mass(.as_mass(Matrix::Diagonal(500))))
  expect_false(.is_identity_mass(matrix(c(1, 0, 0, 0), 2)))
})

test_that("tau leaping reaches the end time or says it did not", {
  j <- RegularJump(function(u, p, t) 100, function(u, p, t, k, mark) k[1], 1)
  pr <- JumpProblem(DiscreteProblem(0, c(0, .5)), PureLeaping(), j)
  v <- sapply(c(1, .5), function(h) solve(pr, SimpleTauLeaping(), dt = h, seed = 17)(.5))
  expect_identical(v[1], v[2])
  expect_gt(v[1], 0)
  j <- MassActionJump(100, list(c(1, 1)), list(c(1, -1)))
  pr <- JumpProblem(DiscreteProblem(1, c(0, 2)), PureLeaping(), j)
  s <- suppressWarnings(solve(pr, SimpleAdaptiveTauLeaping(), dtmin = 1, seed = 17))
  expect_identical(s$retcode, "DtLessThanMin")
  expect_true(max(s$t) < 2)
})

test_that("a Newton step that stops moving is not a root", {
  pr <- NonlinearProblem(function(x, p) (x - 1e20) + 1, 1e20, jac = function(x, p) matrix(1))
  expect_identical(solve(pr, NLsolveJL(method = "newton"), abstol = 1e-8)$retcode, "Stalled")
  expect_identical(solve(NonlinearProblem(function(x, p) x^2 - 2, 1), NLsolveJL())$retcode, "Success")
})

test_that("regression sensitivity is shift and scale safe and leaves collinear pairs undetermined", {
  x <- 1:6
  expect_true(all(is.na(gsa(rbind(x, x), x, RegressionGSA())$standard_regression)))
  set.seed(1); z <- stats::runif(6)
  r <- gsa(rbind(x, z), x + z, RegressionGSA())$standard_regression
  expect_false(anyNA(r))
  X <- matrix(1:4, 1)
  for (Y in list(X, 1e-200 * X, 1e200 * X, 1e14 + X)) expect_equal(c(gsa(Y, 1:4, RegressionGSA())$pearson), 1)
  x <- 1e15 + (0:7)
  r <- gsa(matrix(x, 1), x, RegressionGSA())
  expect_equal(c(r$pearson, r$standard_regression, r$partial_correlation), c(1, 1, 1))
})

test_that("save_idxs is checked against the states", {
  flat <- ODEProblem(function(u, p, t) 0 * u, 2, c(0, 1))
  for (st in c("full", "reduced"))
    expect_error(solve(flat, Tsit5(), save_idxs = 2, save_idxs_storage = st), "between 1 and 1")
  p <- ODEProblem(function(u, p, t) 0 * u, c(10, 20, 30), c(0, 1))
  s <- solve(p, Tsit5(), save_idxs = c(3, 2), save_idxs_storage = "reduced")
  expect_equal(unname(s(1)), c(30, 20))
  expect_identical(get_saved_state_idxs(s), c(3L, 2L))
  expect_identical(get_saved_subsystem(s)$indexp, c(3L, 2L))
})

test_that("terminal parameter resets are differentiated on their right side", {
  p <- ODEProblem(function(u, p, t) p + 0 * u, 0, c(0, 1), 3)
  for (sp in list(c(TRUE, FALSE), c(TRUE, TRUE))) {
    cb <- PresetTimeCallback(1, function(i) i$p[] <- 2 * i$p, save_positions = sp)
    for (sa in list(ForwardSensitivity(), InterpolatingAdjoint())) {
      g <- solve_gradient(function(s) s$p_final[1], p, Tsit5(), callback = cb, save_end = FALSE, sensealg = sa)
      expect_equal(unname(c(g$value, g$p)), c(6, 2), tolerance = 1e-6)
    }
  }
})

test_that("sensitivities call the right hand side with the array shape of u0", {
  f <- function(u, p, t) if (is.matrix(u)) c(p, p) else c(-p, -p)
  p <- ODEProblem(f, matrix(c(1, 2), 2, 1), c(0, 1), 3)
  for (sa in list(ForwardSensitivity(), InterpolatingAdjoint())) {
    g <- solve_gradient(function(s) tail(s[1], 1), p, Tsit5(), saveat = c(0, 1), sensealg = sa)
    expect_equal(unname(c(g$value, g$p)), c(4, 1), tolerance = 1e-6)
  }
})

test_that("a later registration does not change the derivatives of an existing problem", {
  t <- t_nounits
  x <- Variable("nc_rev_sin", t)
  sys <- mtkcompile(System(list(D_nounits(x) ~ sin(x)), t, name = "nc_rev_sin"))
  p <- ODEProblem(sys, list(x ~ 0), c(0, 1))
  old <- get0("cos", envir = .mtk$funs, inherits = FALSE)
  on.exit(if (is.null(old)) rm("cos", envir = .mtk$funs) else assign("cos", old, envir = .mtk$funs), add = TRUE)
  register_symbolic("cos", function(x) 0 * x)
  g <- solve_gradient(function(s) tail(s[x], 1), p, Tsit5(), saveat = c(0, 1), reltol = 1e-10, abstol = 1e-12)
  expect_equal(unname(g$u0), exp(1), tolerance = 1e-6)
})

test_that("forward sensitivities take per state tolerances", {
  two <- ODEProblem(function(u, p, t) -u, c(2, 3), c(0, 1))
  g <- solve_gradient(function(s) sum(tail(as.matrix(s), 1)), two, Tsit5(), saveat = c(0, 1),
                      abstol = c(1e-6, 1e-7), sensealg = ForwardSensitivity())
  expect_equal(unname(g$u0), rep(exp(-1), 2), tolerance = 1e-5)
})

test_that("complex values reach linear solves, sparse Jacobians and gradients", {
  s <- solve(LinearProblem(matrix(1), 1i))
  expect_equal(s$u, 1i)
  expect_identical(s$retcode, "Success")
  A <- matrix(c(1, 1), 2, 1)
  expect_equal(solve(LinearProblem(Matrix::Matrix(A, sparse = TRUE), c(1, 3)))$u, 2)
  t <- t_nounits
  x <- Variable("nc_rev_cx", t)
  sys <- mtkcompile(System(list(D_nounits(x) ~ x^2), t, name = "nc_rev_cx"))
  expect_equal(as.matrix(generate_jacobian(sys, sparse = TRUE)(1 + 2i, numeric(), 0)),
               generate_jacobian(sys)(1 + 2i, numeric(), 0))
  pr <- ODEProblem(function(u, p, t) -p[1] * u, 1i, c(0, 1), 2)
  last1 <- function(s) tail(as.matrix(s)[, 1], 1)
  for (sa in list(ForwardSensitivity(), InterpolatingAdjoint())) {
    g <- solve_gradient(function(s) Im(last1(s)), pr, Tsit5(), saveat = c(0, 1), sensealg = sa,
                        reltol = 1e-10, abstol = 1e-12)
    expect_equal(unname(c(g$value, g$p)), c(exp(-2), -exp(-2)), tolerance = 1e-6)
    expect_equal(unname(g$u0), exp(-2) * 1i, tolerance = 1e-6)
  }
})

test_that("dual medians and trimmed means stay finite", {
  expect_equal(.ad_try_jacobian(function(x) mean(x, trim = .5), rep(1e308, 2))$value, 1e308)
  expect_equal(c(.ad_try_jacobian(function(x) stats::median(x), c(1, 2, 3, 4))$jacobian), c(0, .5, .5, 0))
})

test_that("Krylov methods apply operators matrix free and MINRES and MINARES are their own methods", {
  n <- 64L
  calls <- 0L
  L <- FunctionOperator(function(v, u, p, t) { calls <<- calls + 1L; 2 * v }, numeric(n))
  s <- solve(LinearProblem(L, rep(1, n)), KrylovJL_GMRES())
  expect_equal(s$u, rep(0.5, n))
  expect_lt(calls, 10L)
  for (a in list(KrylovJL_MINRES(), KrylovJL_MINARES(), KrylovJL("minares"))) {
    s <- solve(LinearProblem(diag(c(-1, 1)), c(1, 1)), a, restart = 1, maxiters = 2)
    expect_equal(s$u, c(-1, 1))
    expect_identical(s$retcode, "Success")
  }
  set.seed(3); m <- 40
  Q <- qr.Q(qr(matrix(stats::rnorm(m * m), m)))
  A <- Q %*% diag(c(-3:-1, seq(1, 10, length.out = m - 3))) %*% t(Q); A <- (A + t(A)) / 2
  b <- stats::rnorm(m)
  r_of <- function(a) { s <- suppressWarnings(solve(LinearProblem(A, b), a, maxiters = 6)); b - A %*% s$u }
  rm_ <- r_of(KrylovJL_MINRES()); ra <- r_of(KrylovJL_MINARES())
  # each minimizes its own norm over the Krylov space
  expect_lt(sqrt(sum((A %*% ra)^2)), sqrt(sum((A %*% rm_)^2)) * (1 + 1e-8))
})

test_that("sparse partials and the observed sort are linear in size", {
  x <- lapply(seq_len(2000L), function(i) Symbol(paste0("nc_ps", i)))
  J <- .partials_sparse(lapply(x, function(v) -v), x)
  expect_identical(length(J$x), 2000L)
  expect_identical(J$i, seq_len(2000L))
  E <- Equation
  xs <- lapply(1:6, function(i) Symbol(paste0("nc_to", i))); a <- Symbol("nc_toa")
  obs <- list(E(xs[[3]], xs[[1]] + xs[[2]]), E(xs[[1]], a), E(xs[[2]], 2 * a), E(xs[[4]], xs[[3]]),
              E(xs[[5]], a), E(xs[[6]], xs[[4]]))
  expect_identical(vapply(.toposort_observed(obs), function(o) .key(o$lhs), ""),
                   vapply(xs[c(1, 2, 5, 3, 4, 6)], .key, ""))
  expect_error(.toposort_observed(list(E(xs[[1]], xs[[2]]), E(xs[[2]], xs[[1]]))), "cyclic")
})

test_that("neutral DDEs query derivative history and DDE gradients match their analytic values", {
  pr <- DDEProblem(function(u, h, p, t) h(p, t - 1, deriv = 1L), 1,
                   function(p, t, deriv = 0L) if (deriv == 1L) 1 else 1 + t, c(0, .5), constant_lags = 1, neutral = TRUE)
  expect_equal(solve(pr, MethodOfSteps(Tsit5()))(.5), 1.5)
  pr2 <- DDEProblem(function(u, h, p, t) -h(p, t - 0.5, deriv = 1L) - u, 1,
                    function(p, t, deriv = 0L) if (deriv == 1L) 0 else 1, c(0, 2), constant_lags = .5, neutral = TRUE)
  s <- solve(pr2, MethodOfSteps(Tsit5()), reltol = 1e-10, abstol = 1e-12)
  expect_equal(s(1), exp(-1) + 0.5 * exp(-0.5), tolerance = 1e-7)
  a <- MethodOfSteps(Tsit5())
  p1 <- DDEProblem(function(u, h, p, t) p[1] * h(p, t - 1), 1, function(p, t) 1, c(0, .5), 2, constant_lags = 1)
  g <- solve_gradient(function(s) tail(s[1], 1), p1, a)
  expect_equal(unname(c(g$value, g$p, g$u0)), c(2, 0.5, 1), tolerance = 1e-6)
  # u(2) = u0 - p - p (u0 - p / 2) for u' = -p u(t - 1), h = 1
  p2 <- DDEProblem(function(u, h, p, t) -p[1] * h(p, t - 1), 1, function(p, t) 1, c(0, 2), 0.5, constant_lags = 1)
  g <- solve_gradient(function(s) tail(s[1], 1), p2, a, reltol = 1e-10, abstol = 1e-12)
  expect_equal(unname(c(g$value, g$p, g$u0)), c(0.125, -1.5, 0.5), tolerance = 1e-6)
  p3 <- DDEProblem(function(u, h, p, t) h(p, t - 1), 0, function(p, t) p[1] * t, c(0, 1), 3, constant_lags = 1)
  g <- solve_gradient(function(s) s(1)[1], p3, a)
  expect_equal(unname(c(g$value, g$p)), c(-1.5, -0.5), tolerance = 1e-6)
})
