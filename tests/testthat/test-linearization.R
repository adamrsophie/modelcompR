# Linearization, analysis points and input output functions. Reference
# values are analytic or from ModelingToolkit's tests
# (lib/ModelingToolkitBase/test/analysis_points.jl) and documentation
# (docs/src/tutorials/linear_analysis.md, the linearize docstring).

.blocks <- function() {
  t <- t_nounits; D <- D_nounits
  RealInput <- connector("RealInput", variables = list(list("u(t)", input = TRUE)))
  RealOutput <- connector("RealOutput", variables = list(list("u(t)", output = TRUE)))
  SISO <- mtkmodel("SISO", variables = list("u(t)", "y(t)"),
    components = list(input = quote(RealInput()), output = quote(RealOutput())),
    equations = quote(list(u ~ input$u, y ~ output$u)))
  FirstOrder <- mtkmodel("FirstOrder", extend = list(siso = quote(SISO())),
    parameters = list("T = 1", "k = 1"), variables = list("x(t) = 0"),
    equations = quote(list(D(x) ~ (k * u - x) / T, y ~ x)))
  Gain <- mtkmodel("Gain", extend = list(siso = quote(SISO())), parameters = list("k = 1"),
    equations = quote(list(y ~ k * u)))
  list(FirstOrder = FirstOrder, Gain = Gain)
}

.feedback <- function() {
  b <- .blocks()
  P <- b$FirstOrder(name = "P", k = 1, T = 1)
  C <- b$Gain(name = "C", k = -1)
  System(list(connect(P$output, "plant_output", C$input), connect(C$output, "plant_input", P$input)),
         t_nounits, systems = list(P, C), name = "hej")
}

test_that("a mass spring damper linearizes to the analytic state space model", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); F <- Variable("F", t, input = TRUE); y <- Variable("y", t, output = TRUE)
  m <- Parameter("m", default = 2); k <- Parameter("k", default = 4); c <- Parameter("c", default = 1)
  msd <- System(list(D(D(x)) ~ (F - k * x - c * D(x)) / m, y ~ x), t, name = "msd")
  lin <- linearize(msd, list(F), list(y), op = list(x ~ 0.3, D(x) ~ 0, F ~ 0))
  states <- lin$states
  xi <- match("x(t)", states); vi <- match("x\u02CDt(t)", states)
  expect_false(anyNA(c(xi, vi)))
  lin <- reorder_unknowns(lin, states, states[c(xi, vi)])
  expect_equal(lin$A, rbind(c(0, 1), c(-2, -0.5)), tolerance = 1e-12)
  expect_equal(lin$B, matrix(c(0, 0.5), 2), tolerance = 1e-12)
  expect_equal(lin$C, matrix(c(1, 0), 1), tolerance = 1e-12)
  expect_equal(lin$D, matrix(0, 1, 1))
  sym <- linearize_symbolic(msd, list(F), list(y))
  expect_equal(format(sym$B[[2L]]) %in% c("1/m", "m^(-1)", "m**(-1)"), TRUE)
})

test_that("the linearize docstring example reorders to the documented matrices", {
  t <- t_nounits; D <- D_nounits
  plant <- function(name) {
    x <- Variable("x", t, default = 1); u <- Variable("u", t, default = 0)
    y <- Variable("y", t, default = 0)
    System(list(D(x) ~ -x + u, y ~ x), t, name = name)
  }
  ref_filt <- function(name) {
    x <- Variable("x", t, default = 0); y <- Variable("y", t, default = 0)
    u <- Variable("u", t, default = 0, input = TRUE)
    System(list(D(x) ~ -2 * x + u, y ~ x), t, name = name)
  }
  controller <- function(kp, name) {
    y <- Variable("y", t, default = 0); r <- Variable("r", t, default = 0)
    u <- Variable("u", t, default = 0); k <- Parameter("kp", default = kp)
    System(list(u ~ k * (r - y)), t, name = name)
  }
  f <- ref_filt("f"); c <- controller(1, "c"); p <- plant("p")
  cl <- System(list(f$y ~ c$r, c$u ~ p$u, p$y ~ c$y), t, systems = list(f, c, p), name = "cl")
  n <- toggle_namespacing(cl, FALSE)
  lsys0 <- linearize(cl, list(n$f$u), list(n$p$x), verbose = FALSE)
  # the structural pipeline may keep aliases of f.x and p.x as the states
  state_of <- function(v) {
    if (v %in% lsys0$states) return(v)
    for (o in observed(lsys0$sys)) if (format(o$lhs) == v && format(o$rhs) %in% lsys0$states) return(format(o$rhs))
    v
  }
  lsys <- reorder_unknowns(lsys0, lsys0$states, c(state_of("f\u208Ax(t)"), state_of("p\u208Ax(t)")))
  expect_equal(lsys$A, rbind(c(-2, 0), c(1, -2)), tolerance = 1e-12)
  expect_equal(lsys$B, matrix(c(1, 0), 2), tolerance = 1e-12)
  expect_equal(lsys$C, matrix(c(0, 1), 1), tolerance = 1e-12)
  expect_equal(lsys$D[1, 1], 0)
  expect_error(linearize(cl, list(cl$f$u), list(n$p$x), verbose = FALSE), "namespaced with the name")
})

test_that("sensitivity functions match the linear analysis tutorial", {
  sys <- .feedback()
  expect_s3_class(sys$plant_input, "AnalysisPoint")
  expect_identical(sys$plant_input$name, "hej\u208Aplant_input")
  S <- get_sensitivity(sys, sys$plant_input, verbose = FALSE)
  expect_equal(c(S$A, S$B * S$C, S$D), c(-2, -1, 1), tolerance = 1e-12)
  T <- get_comp_sensitivity(sys, sys$plant_input, verbose = FALSE)
  expect_equal(c(T$A, T$B * T$C, T$D), c(-2, 1, 0), tolerance = 1e-12)
  L <- get_looptransfer(sys, sys$plant_output, verbose = FALSE)
  expect_equal(c(L$A, L$B * L$C, L$D), c(-1, -1, 0), tolerance = 1e-12)
  PS <- linearize(sys, sys$plant_input, sys$plant_output, verbose = FALSE)
  expect_equal(c(PS$A, PS$B * PS$C, PS$D), c(-2, 1, 0), tolerance = 1e-12)
  # names of analysis points, the old style
  S2 <- get_sensitivity(sys, "plant_input", verbose = FALSE)
  expect_equal(c(S2$A, S2$B * S2$C, S2$D), c(-2, -1, 1), tolerance = 1e-12)
  expect_error(linearize(sys, sys$plant_input, sys$plant_input, verbose = FALSE), "both as an input")
})

test_that("analysis points work through the hierarchy", {
  sys <- .feedback()
  nested <- System(list(), t_nounits, systems = list(sys), name = "nested_sys")
  ap <- nested$hej$plant_input
  expect_identical(ap$name, "nested_sys\u208Ahej\u208Aplant_input")
  S <- get_sensitivity(nested, ap, verbose = FALSE)
  expect_equal(c(S$A, S$B * S$C, S$D), c(-2, -1, 1), tolerance = 1e-12)
  L <- get_looptransfer(nested, "hej\u208Aplant_input", verbose = FALSE)
  expect_equal(c(L$A, L$B * L$C, L$D), c(-1, -1, 0), tolerance = 1e-12)
  expect_identical(complete(nested)$hej$plant_input$name, "hej\u208Aplant_input")
})

test_that("open_loop breaks the connection of an analysis point", {
  sys <- .feedback()
  ol <- open_loop(sys, sys$plant_input)
  expect_identical(format(ol$u), "C\u208Aoutput\u208Au(t)")
  expect_identical(format(ol$du[[1L]]), "P\u208Ainput\u208Au(t)")
  ex <- expand_connections(ol$sys)
  expect_false(any(vapply(equations(ex), function(e) format(e) == "P\u208Ainput\u208Au(t) ~ C\u208Aoutput\u208Au(t)", TRUE)))
  m <- linearize(ol$sys, ol$du, list(ol$u), verbose = FALSE)
  expect_equal(c(m$A, m$B * m$C, m$D), c(-1, -1, 0), tolerance = 1e-12)
  # the same connection opened as a loop opening of linearize
  P <- linearize(sys, sys$plant_input, sys$plant_output, loop_openings = list(sys$plant_output),
                 verbose = FALSE)
  expect_equal(c(P$A, P$B * P$C, P$D), c(-1, 1, 0), tolerance = 1e-12)
})

test_that("loop openings hold opened signals at their drivers", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); u <- Variable("u", t, input = TRUE)
  y <- Variable("y", t, output = TRUE); e <- Variable("e", t, output = TRUE)
  # a plant nonlinear in its input, closed by e = 1 - y
  sys <- System(list(D(x) ~ -x + u^2, y ~ x, e ~ 1 - y, connect(e, "loop", u)), t, name = "fb")
  lin <- linearize(sys, list(), list(y), op = list(x ~ 0.25), loop_openings = list(sys$loop), verbose = FALSE)
  expect_equal(lin$A, matrix(-1))
  expect_equal(unname(lin$extras$p[["u(t)"]]), 0.75, tolerance = 1e-10)
})

test_that("linearization problems relinearize at new operating points", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); u <- Variable("u", t, input = TRUE)
  a <- Parameter("a", default = 1)
  sys <- System(list(D(x) ~ -a * x^2 + u), t, name = "s")
  prob <- LinearizationProblem(sys, list(u), list(x), op = list(x ~ 1, u ~ 0))
  expect_equal(solve(prob)$A, matrix(-2))
  prob[x] <- 3
  prob[a] <- 2
  expect_equal(prob[x], 3)
  expect_equal(solve(prob)$A, matrix(-12))
  b <- Parameter("b")
  sysb <- System(list(D(x) ~ -b * x + u), t, name = "sb")
  expect_equal(linearize(sysb, list(u), list(x), op = list(x ~ 1), p = list(b ~ 2))$A, matrix(-2))
  lf <- linearization_function(sys, list(u), list(x), op = list(x ~ 1))
  expect_equal(linearize(lf$sys, lf$lin_fun, op = list(x ~ 2))$A, matrix(-4))
  sol <- solve(ODEProblem(mtkcompile(sys, inputs = list(u)), list(x ~ 1, u ~ 0), c(0, 2)),
               Tsit5(), reltol = 1e-10, abstol = 1e-12)
  res <- linearize(sys, list(u), list(x), op = LinearizationOpPoint(sol, c(0, 1)))
  expect_length(res, 2L)
  expect_length(linearize(lf$sys, lf$lin_fun, op = LinearizationOpPoint(sol, c(0, 1, 2))), 3L)
  expect_equal(res[[2L]]$A[1, 1], -2 / 2, tolerance = 1e-6)
})

test_that("algebraic unknowns are eliminated, or kept in the ModelingToolkit form", {
  # x' = -x + z, 0 = z - 2x + u, y = z
  b <- list(f_x = matrix(-1), f_z = matrix(1), g_x = matrix(-2), g_z = matrix(1),
            f_u = matrix(0), g_u = matrix(1), h_x = matrix(0), h_z = matrix(1), h_u = matrix(0))
  el <- .lin_matrices(b, list(Variable("u", t_nounits)))
  expect_equal(unlist(el), c(A = 1, B = -1, C = 2, D = -1))
  expect_error(.lin_matrices(b, list(Variable("u", t_nounits)), eliminate = FALSE), "Input derivatives")
  full <- .lin_matrices(b, list(Variable("u", t_nounits)), allow_input_derivatives = TRUE, eliminate = FALSE)
  expect_equal(full$A, rbind(c(-1, 1), c(-2, 2)))
  expect_equal(full$B, rbind(c(0, 0), c(0, -1)))
  expect_equal(full$C, matrix(c(0, 1), 1))
  expect_equal(full$D, matrix(0, 1, 2))
  sing <- utils::modifyList(b, list(g_z = matrix(0)))
  expect_error(.lin_matrices(sing, list()), "index > 1")
  T <- matrix(c(2, 1, 0, 1), 2)
  lsys <- list(A = matrix(c(-1, 0, 1, -2), 2), B = matrix(c(1, 0), 2), C = matrix(c(0, 1), 1), D = matrix(0))
  tr <- similarity_transform(lsys, T)
  expect_equal(tr$A, solve(T) %*% lsys$A %*% T)
  expect_equal(tr$C, lsys$C %*% T)
})

test_that("differential algebraic systems linearize numerically and symbolically", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); z <- Variable("z", t); u <- Variable("u", t, input = TRUE)
  y <- Variable("y", t, output = TRUE)
  sys <- System(list(D(x) ~ -x + z, 0 ~ z - 2 * x + u, y ~ z), t, name = "dae")
  lin <- linearize(sys, list(u), list(y), op = list(x ~ 1, u ~ 0))
  expect_equal(unlist(lin[c("A", "B", "C", "D")]), c(A = 1, B = -1, C = 2, D = -1), tolerance = 1e-12)
  # the guess of z selects the branch z = 1 of 0 = z (z - 1)
  w <- Variable("w", t); zg <- Variable("zg", t, guess = 1)
  br <- System(list(D(w) ~ zg * w, 0 ~ zg * (zg - 1)), t, name = "br")
  # ModelingToolkit keeps the algebraic unknown in the state: A = [f_x f_z; G f_x G f_z]
  expect_equal(linearize(br, list(), list(w), op = list(w ~ 1))$A, matrix(c(1, 0, 1, 0), 2), tolerance = 1e-10)
  expect_equal(linearize(br, list(), list(w), op = list(w ~ 1), eliminate_algebraic = TRUE)$A,
               matrix(1), tolerance = 1e-10)
  sym <- linearize_symbolic(sys, list(u), list(y))
  num <- function(M) as.numeric(vapply(M, function(e) as.numeric(e), 0))
  expect_equal(c(num(sym$A), num(sym$B), num(sym$C), num(sym$D)), c(1, -1, 2, -1))
})

test_that("control functions and bound inputs", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); u <- Variable("u", t, input = TRUE); k <- Parameter("k")
  sys <- System(list(D(x) ~ -k * (x + u)), t, name = "s")
  cf <- generate_control_function(sys, list(u))
  expect_equal(cf$f(1, 3, 2, 0), -8)
  expect_identical(vapply(cf$ps, format, ""), "k")
  expect_equal(format(calculate_control_jacobian(cf$io_sys)[1, 1]), "-k")
  expect_equal(format(calculate_control_jacobian(sys)[1, 1]), "-k")
  expect_equal(generate_control_jacobian(cf$io_sys)(1, 2, 0), matrix(-2))
  expect_length(unbound_inputs(sys), 1L)
  r <- Variable("r", t)
  outer <- System(list(sys$u ~ r, D(r) ~ 0), t, systems = list(sys), name = "o")
  expect_true(is_bound(outer, toggle_namespacing(outer, FALSE)$s$u))
  expect_length(unbound_inputs(outer), 0L)
  d <- Variable("d", t, disturbance = TRUE); w <- Variable("w", t)
  sys2 <- System(list(D(x) ~ -k * (x + u) + d + w), t, name = "s2")
  cf3 <- generate_control_function(sys2, list(u), list(d), known_disturbance_inputs = list(w))
  expect_equal(cf3$f(1, 3, 2, 0, 5), -8 + 5)
  cfi <- generate_control_function(sys, list(u), implicit_dae = TRUE)
  expect_equal(cfi$f(0.5, 1, 3, 2, 0), -8 - 0.5)
  fb <- .feedback()
  cf2 <- generate_control_function(fb, fb$plant_input)
  expect_identical(vapply(inputs(cf2$io_sys), format, ""), "P\u208Ainput\u208Au(t)")
  n <- length(cf2$dvs)
  expect_length(cf2$f(rep(0, n), 1, rep(1, length(cf2$ps)), 0), n)
})

test_that("analysis points between components of a chain", {
  t <- t_nounits; D <- D_nounits
  blk <- component(function(name) {
    u <- Variable("u", t, input = TRUE); y <- Variable("y", t, output = TRUE)
    x <- Variable("x", t)
    System(list(D(x) ~ -x + u, y ~ x), t, name = name)
  }, type = "Lag")
  a <- blk(name = "a"); b <- blk(name = "b"); c <- blk(name = "c")
  expect_identical(get_metadata(a)$gui_metadata, "Lag")
  sys <- System(list(a$u ~ 0, connect(a$y, "ab", b$u), connect(b$y, "bc", c$u)), t,
                systems = list(a, b, c), name = "chain")
  lin <- linearize(sys, sys$ab, sys$bc, op = list(), verbose = FALSE)
  expect_length(lin$states, 3L)
  # only b lies between the points: the transfer function is 1/(s + 1)
  expect_equal(c(-lin$C %*% solve(lin$A, lin$B) + lin$D), 1, tolerance = 1e-12)
  expect_equal(c(lin$C %*% solve(diag(3) - lin$A, lin$B) + lin$D), 1 / 2, tolerance = 1e-12)
})
