# Initialization of systems, after docs/src/tutorials/initialization.md of
# ModelingToolkit. Reference values are analytic.

# A semi explicit index one DAE: x' = -x + y on the unit circle. It is
# compiled with the ModelingToolkitBase pipeline so the tests of the
# initialization system do not depend on state selection; the examples
# further down use the default pipeline.
circle_dae <- function() {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t)
  list(x = x, y = y,
       sys = mtkcompile(System(list(D(x) ~ -x + y, 0 ~ x^2 + y^2 - 1), t, name = "dae"), structural = FALSE))
}

test_that("the initialization system has the documented unknowns and equations", {
  m <- circle_dae()
  isys <- generate_initializesystem(m$sys, op = list(m$x ~ 0.6), guesses = list(m$y ~ 1))
  expect_true(modelcompR:::.f(isys, "is_initializesystem"))
  expect_null(get_iv(isys))
  expect_equal(vapply(unknowns(isys), format, ""), c("x(t)", "y(t)", "x\u02CDt(t)"))
  eqs <- vapply(equations(isys), format, "")
  expect_length(eqs, 3L)
  expect_true("x(t) ~ Initial(x(t))" %in% eqs)
  expect_true(any(grepl("^x\u02CDt\\(t\\) ~ ", eqs)))
  # the constant of the condition is the value of the Initial parameter
  ix <- Initial(m$x)
  expect_true(isinitial(ix))
  expect_equal(smap_get(initial_conditions(isys), ix), 0.6)
  expect_equal(smap_get(guesses(isys), m$y), 1)
  # mtkcompile() accepts it when it is not required to be balanced
  expect_s3_class(mtkcompile(isys, fully_determined = FALSE), "System")
})

test_that("guesses choose the branch of the algebraic unknown", {
  m <- circle_dae()
  prob <- ODEProblem(m$sys, list(m$x ~ 0.6), c(0, 1), guesses = list(m$y ~ 1))
  expect_equal(unname(prob$u0), c(0.6, 0.8), tolerance = 1e-10)
  expect_equal(prob$initialization_data$retcode, "Success")
  prob <- ODEProblem(m$sys, list(m$x ~ 0.6), c(0, 1), guesses = list(m$y ~ -1))
  expect_equal(unname(prob$u0), c(0.6, -0.8), tolerance = 1e-10)
  # the algebraic constraint holds at the initial point
  expect_equal(prob[m$x^2 + m$y^2], 1, tolerance = 1e-10)
  # remake() runs initialization again
  p2 <- remake(prob, u0 = list(m$x ~ 0.8))
  expect_equal(unname(p2$u0), c(0.8, -0.6), tolerance = 1e-10)
  # numeric values: x is the condition, y the new starting value
  p3 <- remake(prob, u0 = c(0, 0.5))
  expect_equal(unname(p3$u0), c(0, 1), tolerance = 1e-10)
  # the initialization problem can be inspected and solved directly
  iprob <- InitializationProblem(m$sys, 0, list(m$x ~ 0.6), guesses = list(m$y ~ 1))
  # a square system is solved block by block (use_scc = TRUE, as ModelingToolkit)
  expect_s3_class(iprob, "SCCNonlinearProblem")
  expect_s3_class(InitializationProblem(m$sys, 0, list(m$x ~ 0.6), guesses = list(m$y ~ 1), use_scc = FALSE),
                  "NonlinearProblem")
  sol <- solve(iprob)
  expect_equal(sol$retcode, "Success")
  expect_equal(sol[m$y], 0.8, tolerance = 1e-10)
  expect_equal(sol[m$x], 0.6, tolerance = 1e-10)
})

test_that("consistent values pass and inconsistent values fail", {
  m <- circle_dae()
  prob <- ODEProblem(m$sys, list(m$x ~ 0.6, m$y ~ 0.8), c(0, 1))
  expect_equal(unname(prob$u0), c(0.6, 0.8))
  expect_warning(expect_warning(bad <- ODEProblem(m$sys, list(m$x ~ 0.6, m$y ~ 0.5), c(0, 1)),
                                "overdetermined"), "Initialization failed")
  expect_equal(bad$initialization_data$retcode, "InitialFailure")
  # with initialization switched off the values are taken as they are
  raw <- ODEProblem(m$sys, list(m$x ~ 0.6, m$y ~ 0.5), c(0, 1), check_initialization = FALSE)
  expect_null(raw$initialization_data)
  expect_equal(unname(raw$u0), c(0.6, 0.5))
  # an inconsistent least squares problem: nonzero residual at its minimum
  iprob <- suppressWarnings(InitializationProblem(m$sys, 0, list(m$x ~ 0.6, m$y ~ 0.5)))
  expect_s3_class(iprob, "NonlinearLeastSquaresProblem")
  expect_gt(max(abs(solve(iprob)$resid)), 0.1)
})

test_that("a failed initialization is the return code of solve", {
  m <- circle_dae()
  bad <- suppressWarnings(ODEProblem(m$sys, list(m$x ~ 0.6, m$y ~ 0.5), c(0, 1)))
  expect_s3_class(bad, "InitializationFailure")
  expect_s3_class(bad, "ODEProblem")
  sol <- solve(bad)
  expect_equal(sol$retcode, "InitialFailure")
  expect_false(successful_retcode(sol))
  expect_length(sol$t, 1L)
  # fixing the conditions clears the mark
  fixed <- suppressWarnings(remake(bad, u0 = list(m$y ~ 0.8)))
  expect_false(inherits(fixed, "InitializationFailure"))
  skip_if(!exists("Rodas5P", mode = "function"), "no mass matrix solver available")
  good <- ODEProblem(m$sys, list(m$x ~ 0.6), c(0, 1), guesses = list(m$y ~ 1))
  sol <- solve(good, Rodas5P(), reltol = 1e-8, abstol = 1e-8)
  expect_equal(sol$retcode, "Success")
  expect_equal(max(abs(sol[m$x^2 + m$y^2] - 1)), 0, tolerance = 1e-6)
})

test_that("underdetermined and overdetermined systems warn or stop", {
  m <- circle_dae()
  expect_warning(prob <- ODEProblem(m$sys, list(), c(0, 1), guesses = list(m$x ~ 0.5, m$y ~ 1)),
                 "Initialization system is underdetermined")
  # least squares from the guess (0.5, 1) lands on the circle
  expect_equal(sum(prob$u0^2), 1, tolerance = 1e-8)
  expect_equal(prob$initialization_data$retcode, "Success")
  expect_error(ODEProblem(m$sys, list(), c(0, 1), guesses = list(m$x ~ 0.5, m$y ~ 1), fully_determined = TRUE),
               "Initialization system is underdetermined")
  expect_error(ODEProblem(m$sys, list(m$x ~ 0.6, m$y ~ 0.8, D_nounits(m$x) ~ 0.2), c(0, 1), fully_determined = TRUE),
               "Initialization system is overdetermined")
  # consistent redundant conditions: overdetermined but solvable
  expect_warning(ok <- ODEProblem(m$sys, list(m$x ~ 0.6, m$y ~ 0.8, D_nounits(m$x) ~ 0.2), c(0, 1)),
                 "Initialization system is overdetermined")
  expect_equal(ok$initialization_data$retcode, "Success")
  expect_silent(ODEProblem(m$sys, list(), c(0, 1), guesses = list(m$x ~ 0.5, m$y ~ 1),
                           warn_initialize_determined = FALSE))
  # the rank analysis names the free unknowns
  iprob <- suppressWarnings(InitializationProblem(m$sys, 0, list(), guesses = list(m$x ~ 0.5, m$y ~ 1)))
  res <- analyze_initialization_jacobian(iprob, verbose = FALSE)
  expect_equal(res$rank, 1L); expect_equal(res$nullity, 1L); expect_equal(res$redundancy, 0L)
  expect_equal(sum(res$underdetermined_unknowns$weight), 1, tolerance = 1e-8)
  expect_message(analyze_initialization_jacobian(iprob), "Underdetermined unknowns")
  full <- analyze_initialization_jacobian(ODEProblem(m$sys, list(m$x ~ 0.6), c(0, 1), guesses = list(m$y ~ 1)),
                                          verbose = FALSE)
  expect_equal(full$nullity, 0L)
})

test_that("missing guesses follow MissingGuessValue", {
  m <- circle_dae()
  # the default is a reproducible pseudo random start in [0, 1)
  prob <- ODEProblem(m$sys, list(m$x ~ 0.6), c(0, 1))
  expect_equal(prob[m$x], 0.6)
  expect_equal(abs(prob[m$y]), 0.8, tolerance = 1e-9)
  expect_identical(ODEProblem(m$sys, list(m$x ~ 0.6), c(0, 1))$u0, prob$u0)
  expect_error(ODEProblem(m$sys, list(m$x ~ 0.6), c(0, 1), missing_guess_value = MissingGuessValue$Error()),
               "Cyclic guesses detected")
  up <- ODEProblem(m$sys, list(m$x ~ 0.6), c(0, 1), missing_guess_value = MissingGuessValue$Constant(0.5))
  expect_equal(up[m$y], 0.8, tolerance = 1e-9)
  down <- ODEProblem(m$sys, list(m$x ~ 0.6), c(0, 1), missing_guess_value = MissingGuessValue$Constant(-0.5))
  expect_equal(down[m$y], -0.8, tolerance = 1e-9)
})

test_that("derivatives and observed variables can be initial conditions", {
  # the Lotka Volterra example of the tutorial
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t); z <- Variable("z", t)
  al <- Parameter("al", default = 1.5); be <- Parameter("be", default = 1)
  ga <- Parameter("ga", default = 3); de <- Parameter("de", default = 1)
  sys <- mtkcompile(System(list(D(x) ~ al * x - be * x * y, D(y) ~ -ga * y + de * x * y, z ~ x + y), t, name = "lv"))
  # start at the steady state by setting both derivatives to zero
  prob <- ODEProblem(sys, list(D(x) ~ 0, D(y) ~ 0), c(0, 10), guesses = list(x ~ 2.5, y ~ 2))
  expect_equal(c(prob[x], prob[y]), c(3, 1.5), tolerance = 1e-8)
  # a condition on the observed z ~ x + y, with D(x) = 0: y = al / be, x = z - y
  prob2 <- ODEProblem(sys, list(D(x) ~ 0, z ~ 2), c(0, 10), guesses = list(x ~ 1, y ~ 1))
  expect_equal(c(prob2[x], prob2[y]), c(0.5, 1.5), tolerance = 1e-8)
  expect_equal(prob2[al * x - be * x * y], 0, tolerance = 1e-8)
  # only the constant changes: the compiled initialization system is kept
  prob3 <- remake(prob2, u0 = list(z ~ 3))
  expect_equal(c(prob3[x], prob3[y]), c(1.5, 1.5), tolerance = 1e-8)
  expect_identical(prob3$initialization_data$signature, prob2$initialization_data$signature)
  # initial conditions may be expressions of parameters
  prob4 <- ODEProblem(sys, list(x ~ 2 * al, y ~ 1), c(0, 10))
  expect_equal(c(prob4[x], prob4[y]), c(3, 1))
})

test_that("initialization equations of the problem and of the system are used", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t)
  eqs <- list(D(x) ~ -x * y, D(y) ~ x - y)
  sys <- mtkcompile(System(eqs, t, name = "ie"))
  prob <- ODEProblem(sys, list(x ~ 1), c(0, 1), guesses = list(y ~ 1), initialization_eqs = list(y ~ 2 * x))
  expect_equal(c(prob[x], prob[y]), c(1, 2), tolerance = 1e-9)
  sys2 <- mtkcompile(System(eqs, t, name = "ie", initialization_eqs = list(x + y ~ 5)))
  prob2 <- ODEProblem(sys2, list(x ~ 1), c(0, 1), guesses = list(y ~ 1))
  expect_equal(c(prob2[x], prob2[y]), c(1, 4), tolerance = 1e-9)
  # a nonlinear initialization equation: x^2 + y^2 = 25 with x = 3
  prob3 <- ODEProblem(sys, list(x ~ 3), c(0, 1), guesses = list(y ~ 1), initialization_eqs = list(x^2 + y^2 ~ 25))
  expect_equal(c(prob3[x], prob3[y]), c(3, 4), tolerance = 1e-8)
})

test_that("parameters bound to NA are solved for", {
  # the parameter initialization example of the tutorial
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); y <- Variable("y", t); total <- Parameter("total", default = NA)
  sys <- mtkcompile(System(list(D(x) ~ -x, total ~ x + y), t, name = "paraminit"))
  prob <- ODEProblem(sys, list(x ~ 1, y ~ 2), c(0, 1))
  expect_equal(prob$ps[total], 3, tolerance = 1e-10)
  # solve for x given total and y: the value of x has to be removed
  prob2 <- remake(prob, u0 = list(y ~ 1, list(x, NULL)), p = list(total ~ 4))
  expect_equal(prob2[x], 3, tolerance = 1e-10)
  expect_equal(prob2$ps[total], 4)
  # keeping x makes the system overdetermined and unsatisfiable
  expect_warning(expect_warning(prob3 <- remake(prob, u0 = list(y ~ 1), p = list(total ~ 4)), "overdetermined"),
                 "Initialization failed")
  expect_equal(prob3$initialization_data$retcode, "InitialFailure")
})

test_that("Initial parameters hold the initial value of a variable", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); z <- Variable("z", t); r <- Parameter("r")
  sys <- mtkcompile(System(list(D(x) ~ -r * x, 0 ~ z^3 + z - x), t, name = "ib",
                           bindings = list(r ~ Initial(z)), guesses = list(z ~ 1)))
  # Initial(z) is hidden from parameters() unless asked for
  expect_length(parameters(sys), 1L)
  expect_true(any(vapply(parameters(sys, initial_parameters = TRUE), isinitial, TRUE)))
  prob <- ODEProblem(sys, list(x ~ 2), c(0, 1))
  # z solves z^3 + z = 2, and r is bound to that initial value
  expect_equal(c(prob[x], prob[z]), c(2, 1), tolerance = 1e-9)
  expect_equal(prob$ps[Initial(z)], 1, tolerance = 1e-9)
  expect_equal(prob$ps[r], 1, tolerance = 1e-9)
  prob2 <- remake(prob, u0 = list(x ~ 10))
  expect_equal(prob2$ps[r], 2, tolerance = 1e-9)
  # only the Initial parameters a model refers to are added on completion
  y <- Variable("y", t); k <- Parameter("k")
  plain <- mtkcompile(System(list(D(y) ~ -k * y), t, name = "plain"))
  expect_length(parameters(plain, initial_parameters = TRUE), 1L)
  all_init <- add_initialization_parameters(System(list(D(y) ~ -k * y), t, name = "plain"))
  inits <- vapply(Filter(isinitial, parameters(all_init, initial_parameters = TRUE)), format, "")
  expect_setequal(inits, c("Initial(y(t))", "Initial(y\u02CDt(t))"))
})

test_that("ordinary differential equations are left alone", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); k <- Parameter("k", default = 2)
  sys <- mtkcompile(System(list(D(x) ~ -k * x), t, name = "ode"))
  prob <- ODEProblem(sys, list(x ~ 1), c(0, 1))
  expect_null(prob$initialization_data)
  expect_equal(unname(prob$p), 2)
  # an unknown without a value is an underdetermined initialization
  expect_warning(p2 <- ODEProblem(sys, list(), c(0, 1), guesses = list(x ~ 0.25)), "underdetermined")
  expect_equal(unname(p2$u0), 0.25)
})

test_that("time independent initialization solves for parameters", {
  x <- Variable("x"); a <- Parameter("a", default = NA); b <- Parameter("b")
  sys <- complete(System(list(0 ~ x^2 - a), unknowns = list(x), ps = list(a, b), name = "pinit",
                         initialization_eqs = list(a ~ 2 * b)))
  isys <- generate_initializesystem(sys, op = list(b ~ 8))
  expect_equal(vapply(unknowns(isys), format, ""), "a")
  prob <- NonlinearProblem(sys, list(x ~ 1, b ~ 8))
  expect_equal(prob$ps[a], 16, tolerance = 1e-10)
  expect_equal(solve(prob)[x], 4, tolerance = 1e-9)
})

test_that("the pendulum of the tutorial initializes", {
  skip_if(is.null(modelcompR:::.mtk$hooks$structural_compile), "needs the structural compiler (index reduction)")
  # the setup of tools/julia/reference.jl, whose solution starts at (1, 0)
  ref <- utils::read.csv(test_path("fixtures", "julia", "pendulum.csv"), check.names = FALSE)
  t <- t_nounits; D <- D_nounits
  g <- Parameter("g", default = 9.81); L <- Parameter("L", default = 1)
  px <- Variable("px", t); py <- Variable("py", t); lam <- Variable("lam", t)
  pend <- mtkcompile(System(list(D(D(px)) ~ lam * px, D(D(py)) ~ lam * py - g, px^2 + py^2 ~ L^2), t, name = "pend"))
  # px = 1 and a derivative condition determine the state; lambda and py are solved for
  prob <- suppressWarnings(ODEProblem(pend, list(px ~ 1, D(px) ~ 0), c(0, 2), guesses = list(lam ~ 0, py ~ 0)))
  expect_equal(prob[px], ref[["px(t)"]][1L], tolerance = 1e-8)
  expect_equal(prob[px^2 + py^2], 1, tolerance = 1e-8)
  # py = 0 is a double root of the constraint, so it is found less accurately
  expect_lt(abs(prob[py] - ref[["py(t)"]][1L]), 1e-3)
  # at rest in the far right position the tension vanishes
  expect_lt(abs(prob[lam]), 1e-3)
})

test_that("the guesses example matches the Julia reference", {
  # tools/julia/reference.jl: D(a) ~ -a + b, 0 ~ b^3 + b - c a with c = 2,
  # a(0) = 1 and the guess b = 0.5; b(0) = 1 and the solution is constant
  ref <- utils::read.csv(test_path("fixtures", "julia", "dae1.csv"), check.names = FALSE)
  t <- t_nounits; D <- D_nounits
  a <- Variable("a", t); b <- Variable("b", t); cc <- Parameter("c", default = 2)
  sys <- mtkcompile(System(list(D(a) ~ -a + b, 0 ~ b^3 + b - cc * a), t, name = "dae1"))
  prob <- ODEProblem(sys, list(a ~ 1), c(0, 3), guesses = list(b ~ 0.5))
  expect_equal(prob[a], ref[["a(t)"]][1L], tolerance = 1e-10)
  expect_equal(prob[b], ref[["b(t)"]][1L], tolerance = 1e-10)
  # a different condition: b(0) = 2 gives a(0) = (8 + 2) / 2
  prob2 <- ODEProblem(sys, list(b ~ 2), c(0, 3), guesses = list(a ~ 1))
  expect_equal(prob2[a], 5, tolerance = 1e-10)
  skip_if(!exists("Rodas5P", mode = "function"), "no mass matrix solver available")
  sol <- solve(prob, Rodas5P(), abstol = 1e-10, reltol = 1e-10, saveat = ref$t)
  expect_equal(sol[a], ref[["a(t)"]], tolerance = 1e-8)
  expect_equal(sol[b], ref[["b(t)"]], tolerance = 1e-8)
})

test_that("changed conditions are seen even when the numbers stay the same", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t)
  sys <- mtkcompile(System(list(D(x) ~ -x), t, name = "c1"))
  # x = 1 with D(x) = -1 is redundant but consistent
  prob <- suppressWarnings(ODEProblem(sys, list(x ~ 1, D(x) ~ -1), c(0, 1)))
  expect_equal(prob$initialization_data$retcode, "Success")
  # D(x) = 0 contradicts it, although u0 and p do not change
  bad <- suppressWarnings(remake(prob, u0 = list(D(x) ~ 0)))
  expect_s3_class(bad, "InitializationFailure")
  good <- suppressWarnings(remake(bad, u0 = list(D(x) ~ -1)))
  expect_false(inherits(good, "InitializationFailure"))
  # numeric values given to remake() are conditions for a fully specified DAE
  m <- circle_dae()
  full <- ODEProblem(m$sys, list(m$x ~ 0.6, m$y ~ 0.8), c(0, 1))
  expect_equal(unname(remake(full, u0 = c(0.8, 0.6))$u0), c(0.8, 0.6))
  off <- suppressWarnings(remake(full, u0 = c(0.8, 0.5)))
  expect_equal(off$initialization_data$retcode, "InitialFailure")
})

test_that("Initial(x) of a solved variable is that variable during initialization", {
  t <- t_nounits; D <- D_nounits
  x <- Variable("x", t); z <- Variable("z", t); r <- Parameter("r")
  sys <- mtkcompile(System(list(D(x) ~ -r * x, 0 ~ z^3 + z - x), t, name = "ic",
                           bindings = list(r ~ Initial(z)), guesses = list(z ~ 1.3)))
  # x = 2 gives z = 1 and r = 1, so D(x) = -2 is consistent and -3 is not
  ok <- suppressWarnings(ODEProblem(sys, list(x ~ 2, D(x) ~ -2), c(0, 1)))
  expect_equal(ok$initialization_data$retcode, "Success")
  expect_equal(c(ok[x], ok[z]), c(2, 1), tolerance = 1e-8)
  expect_equal(ok$ps[r], 1, tolerance = 1e-8)
  bad <- suppressWarnings(ODEProblem(sys, list(x ~ 2, D(x) ~ -3), c(0, 1)))
  expect_equal(bad$initialization_data$retcode, "InitialFailure")
  # the derivative condition alone determines x through r = z: z * (z^3 + z) = 2
  der <- ODEProblem(sys, list(D(x) ~ -2), c(0, 1), guesses = list(x ~ 1.5))
  expect_equal(c(der[x], der[z]), c(2, 1), tolerance = 1e-8)
})

test_that("initialization uses the problem types ModelingToolkit chooses", {
  t <- t_nounits; D <- D_nounits
  # an affine initialization system is a LinearProblem
  lx <- Variable("in2_lx", t); ly <- Variable("in2_ly", t); k <- Parameter("in2_k", default = 2)
  lsys <- mtkcompile(System(list(D(lx) ~ -k * lx + ly, D(ly) ~ -ly), t, name = "in2lin",
                            initialization_eqs = list(2 * lx + ly ~ 3, lx - ly ~ 0)))
  lp <- ODEProblem(lsys, list(), c(0, 1))
  expect_s3_class(lp$initialization_data$initializeprob, "LinearProblem")
  expect_equal(unname(lp$u0[c("in2_lx(t)", "in2_ly(t)")]), c(1, 1), tolerance = 1e-12)
  ip <- InitializationProblem(lsys, 0, list())
  expect_s3_class(ip, "LinearProblem")
  sol <- solve(ip)
  expect_equal(c(sol[lx], sol[ly]), c(1, 1), tolerance = 1e-12)
  # overdetermined and affine: least squares, which fails for inconsistent values
  lp2 <- suppressWarnings(ODEProblem(lsys, list(ly ~ 0.5), c(0, 1), warn_initialize_determined = FALSE))
  expect_s3_class(lp2$initialization_data$initializeprob, "LinearProblem")
  expect_s3_class(lp2, "InitializationFailure")
  # a square nonlinear system: SCC blocks by default, a plain problem without
  cx <- Variable("in2_cx", t); cy <- Variable("in2_cy", t)
  csys <- mtkcompile(System(list(D(cx) ~ -cx + cy, 0 ~ cx^2 + cy^2 - 1), t, name = "in2c"))
  cp <- ODEProblem(csys, list(cx ~ 0.6), c(0, 1), guesses = list(cy ~ 1))
  expect_s3_class(cp$initialization_data$initializeprob, "SCCNonlinearProblem")
  expect_equal(unname(cp$u0[c("in2_cx(t)", "in2_cy(t)")]), c(0.6, 0.8), tolerance = 1e-10)
  np <- ODEProblem(csys, list(cx ~ 0.6), c(0, 1), guesses = list(cy ~ 1), use_scc = FALSE)
  expect_s3_class(np$initialization_data$initializeprob, "NonlinearProblem")
  expect_equal(unname(np$u0), unname(cp$u0), tolerance = 1e-10)
  # polynomial homotopy continuation when requested: the real root nearest
  # to the guess
  hp <- ODEProblem(csys, list(cx ~ 0.6), c(0, 1), guesses = list(cy ~ -1), use_homotopy_continuation = TRUE)
  expect_s3_class(hp$initialization_data$initializeprob, "HomotopyContinuationProblem")
  expect_equal(hp[cy], -0.8, tolerance = 1e-10)
  # homotopy() nodes: a HomotopyProblem without SCC, a homotopy block with
  hx <- Variable("in2_hx", t); hy <- Variable("in2_hy", t)
  hsys <- mtkcompile(System(list(D(hx) ~ -hx, 0 ~ homotopy(atan(hy - 3), hy) + hx - 1), t, name = "in2h"))
  hpp <- ODEProblem(hsys, list(hx ~ 1), c(0, 1), guesses = list(hy ~ 12), use_scc = FALSE)
  expect_s3_class(hpp$initialization_data$initializeprob, "HomotopyProblem")
  expect_equal(hpp[hy], 3, tolerance = 1e-9)
  hps <- ODEProblem(hsys, list(hx ~ 1), c(0, 1), guesses = list(hy ~ 12))
  expect_equal(hps$initialization_data$initializeprob$scc_kinds, "homotopy")
  expect_equal(hps[hy], 3, tolerance = 1e-9)
  # the underdetermined pendulum of the structural tests is least squares
  expect_s3_class(InitializationProblem(csys, 0, list(), guesses = list(cx ~ 0.6, cy ~ 0.8),
                                        warn_initialize_determined = FALSE), "NonlinearLeastSquaresProblem")
})

test_that("Initial parameters can be set and remake can skip initialization", {
  t <- t_nounits; D <- D_nounits
  cx <- Variable("in3_cx", t); cy <- Variable("in3_cy", t); a <- Parameter("in3_a", default = 1)
  csys <- mtkcompile(System(list(D(cx) ~ -a * cx + cy, 0 ~ cx^2 + cy^2 - 1), t, name = "in3c"))
  prob <- ODEProblem(csys, list(cx ~ 0.6), c(0, 1), guesses = list(cy ~ 1))
  expect_equal(prob$ps[Initial(cx)], 0.6)
  # the constant of the condition x ~ Initial(x) changes, and initialization runs again
  prob$ps[Initial(cx)] <- 0.8
  expect_equal(prob[cx], 0.8, tolerance = 1e-12)
  expect_equal(prob[cy], 0.6, tolerance = 1e-10)
  expect_equal(prob$ps[Initial(cx)], 0.8)
  # Initial of a variable that is solved for has no effect (initialization tutorial)
  before <- prob$u0
  prob$ps[Initial(cy)] <- 0.3
  expect_equal(prob$u0, before)
  # plain parameters
  prob$ps[a] <- 3
  expect_equal(prob$ps[a], 3)
  expect_equal(unname(prob$f$f(prob$u0, prob$p, 0)[1L]), -3 * 0.8 + 0.6, tolerance = 1e-10)
  expect_error(prob$ps[Variable("in3_nope", t)] <- 1, "not a parameter")
  # remake with check_initialization = FALSE keeps the values as given
  p2 <- remake(prob, u0 = c(0.1, 0.2), check_initialization = FALSE)
  expect_equal(unname(p2$u0), c(0.1, 0.2))
  p3 <- remake(prob, u0 = list(cx ~ 0.6), check_initialization = FALSE)
  expect_equal(p3[cy], prob[cy])
  p4 <- remake(prob, u0 = list(cx ~ 0.6))
  expect_equal(p4[cy], 0.8, tolerance = 1e-10)
  # other fields are written as before
  prob$tspan <- c(0, 2)
  expect_equal(prob$tspan, c(0, 2))
  # the form of the initialization tutorial: Initial(D(y)) of a second order
  # equation, and Initial(y) after a numeric remake
  y <- Variable("in4_y", t); w <- Parameter("in4_w", default = 2)
  osys <- mtkcompile(System(list(D(D(y)) ~ -w * y), t, name = "in4"))
  op <- ODEProblem(osys, list(y ~ 1, D(y) ~ 0), c(0, 1))
  expect_equal(op$ps[Initial(D(y))], 0)
  op$ps[Initial(D(y))] <- 0.5
  expect_equal(unname(op$u0), c(1, 0.5))
  expect_equal(op$ps[Initial(D(y))], 0.5)
  o2 <- remake(op, u0 = c(2, 0.25))
  expect_equal(o2$ps[Initial(y)], 2)
})
