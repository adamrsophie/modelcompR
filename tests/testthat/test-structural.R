# Structural simplification: the mtkcompile pipeline of ModelingToolkit
# v11.47 (alias elimination, index reduction, state selection, tearing).
# Expected values are analytic or are the expectations of the
# ModelingToolkit.jl test suite, cited by file and line at each test.

t <- t_nounits; D <- D_nounits
fmt <- function(xs) vapply(xs, format, "")
lhs <- function(eqs) fmt(lapply(eqs, `[[`, "lhs"))
is_diff <- function(sys) vapply(equations(sys), function(e) .is_derivative(e$lhs), TRUE)

# Models compiled once and shared by the tests below.
rc_eqs <- local({
  variables("vs(t)", "is(t)", "vr(t)", "ir(t)", "vc(t)", "ic(t)", "va(t)", "vb(t)", "vg(t)")
  parameters("R = 2", "C = 1", "V = 1")
  list(vs ~ V, vs ~ va - vg, vr ~ va - vb, vr ~ R * ir, vc ~ vb - vg, D(vc) ~ ic / C,
       ir ~ ic, 0 ~ is + ir, vg ~ 0, 0 ~ -ir - is)
})
rc <- mtkcompile(System(rc_eqs[1:9], t, name = "rc"))
rc_redundant <- System(rc_eqs, t, name = "rc")

loop <- local({
  variables("h(t)", "u3(t)"); parameters("p = 0.3")
  u1 <- Variable("u1", t, guess = 1); u2 <- Variable("u2", t, guess = 0.5)
  System(list(D(h) ~ u1 - h, 0 ~ u1 - sin(u2) - h, 0 ~ u2 - cos(u1) + p, u3 ~ u1 + u2), t, name = "loop")
})
loop_c <- mtkcompile(loop)

pend_vars <- local({
  list(x = Variable("x", t, state_priority = 10, guess = 1), y = Variable("y", t, guess = -1),
       T = Variable("T", t, guess = 1), g = Parameter("g", default = 9.81), L = Parameter("L", default = 1))
})
pend <- with(pend_vars, System(list(D(D(x)) ~ T * x, D(D(y)) ~ T * y - g, 0 ~ x^2 + y^2 - L^2), t, name = "pend"))
pend_c <- mtkcompile(pend)

# Structures and trajectories written by tools/julia/reference.jl with
# ModelingToolkit 11.46.1.
julia_fixture <- function(name) {
  f <- test_path("fixtures", "julia", name)
  if (!file.exists(f)) return(NULL)
  utils::read.csv(f, check.names = FALSE, encoding = "UTF-8", stringsAsFactors = FALSE)
}
julia_names <- function(st, kind) st$name[st$kind == kind]

test_that("unknowns and observed variables match ModelingToolkit for the reference models", {
  # Cartesian pendulum (index 3): the same states, in the same order
  st <- julia_fixture("pendulum_structure.csv")
  skip_if(is.null(st))
  px <- Variable("px", t); py <- Variable("py", t); lam <- Variable("λ", t)
  g <- Parameter("g", default = 9.81); L <- Parameter("L", default = 1)
  jp <- mtkcompile(System(list(D(D(px)) ~ lam * px, D(D(py)) ~ lam * py - g, px^2 + py^2 ~ L^2), t, name = "pend"))
  expect_equal(fmt(unknowns(jp)), julia_names(st, "unknown"))
  expect_equal(fmt(observables(jp)), julia_names(st, "observed"))
  expect_setequal(fmt(parameters(jp)), julia_names(st, "parameter"))
  # index 1 system kept as a mass matrix problem
  st <- julia_fixture("dae1_structure.csv")
  a <- Variable("a", t); b <- Variable("b", t); cc <- Parameter("c", default = 2)
  dae1 <- mtkcompile(System(list(D(a) ~ -a + b, 0 ~ b^3 + b - cc * a), t, name = "dae1"))
  expect_equal(fmt(unknowns(dae1)), julia_names(st, "unknown"))
  expect_length(observed(dae1), 0L)
  # second order oscillator with an observed variable
  st <- julia_fixture("osc_structure.csv")
  y <- Variable("y", t); v <- Variable("v", t); k <- Parameter("k", default = 2)
  osc <- mtkcompile(System(list(D(D(y)) ~ -k * y, v ~ 2 * D(y)), t, name = "osc"))
  expect_equal(fmt(unknowns(osc)), julia_names(st, "unknown"))
  expect_equal(fmt(observables(osc)), julia_names(st, "observed"))
  skip_if_not(exists("Rodas5P"))
  # Julia starts the pendulum horizontally; with py as the state the reduced
  # system is singular at the bottom (px = 0), so compare up to t = 0.5
  ref <- julia_fixture("pendulum.csv")
  ref <- ref[ref$t <= 0.5, ]
  sol <- solve(ODEProblem(jp, list(px ~ 1, D(px) ~ 0), c(0, 0.5), guesses = list(lam ~ 0, py ~ 0)),
               Rodas5P(), abstol = 1e-10, reltol = 1e-10)
  got <- t(vapply(ref$t, function(ti) as.numeric(sol(ti, idxs = list(px, py))), numeric(2)))
  expect_equal(got[, 1L], ref[["px(t)"]], tolerance = 1e-5)
  expect_equal(got[, 2L], ref[["py(t)"]], tolerance = 1e-5)
  ref <- julia_fixture("dae1.csv")
  sol <- solve(ODEProblem(dae1, list(a ~ 1), c(0, 3), guesses = list(b ~ 0.5)), Rodas5P(), abstol = 1e-10, reltol = 1e-10)
  got <- t(vapply(ref$t, function(ti) as.numeric(sol(ti, idxs = list(a, b))), numeric(2)))
  expect_equal(got[, 1L], ref[["a(t)"]], tolerance = 1e-7)
  expect_equal(got[, 2L], ref[["b(t)"]], tolerance = 1e-7)
  # Robertson kinetics: three states, compared on the first hundred time units
  ref <- julia_fixture("robertson.csv")
  ref <- ref[ref$t <= 100, ]
  y1 <- Variable("y1", t); y2 <- Variable("y2", t); y3 <- Variable("y3", t)
  k1 <- Parameter("k1", default = 0.04); k2 <- Parameter("k2", default = 3e7); k3 <- Parameter("k3", default = 1e4)
  rob <- mtkcompile(System(list(D(y1) ~ -k1 * y1 + k3 * y2 * y3, D(y2) ~ k1 * y1 - k2 * y2^2 - k3 * y2 * y3,
                                D(y3) ~ k2 * y2^2), t, name = "rober"))
  expect_setequal(fmt(unknowns(rob)), c("y1(t)", "y2(t)", "y3(t)"))
  sol <- solve(ODEProblem(rob, list(y1 ~ 1, y2 ~ 0, y3 ~ 0), c(0, 100)), Rodas5P(), abstol = 1e-10, reltol = 1e-8)
  got <- t(vapply(ref$t, function(ti) as.numeric(sol(ti, idxs = list(y1, y2, y3))), numeric(3)))
  expect_equal(got, unname(as.matrix(ref[, -1L])), tolerance = 1e-6)
})

test_that("state order and alias targets are those of ModelingToolkit", {
  # unknowns and observed variables printed by ModelingToolkit 11.46.1
  # (tools/julia) for the same models: states are not kept in declaration
  # order, and an alias group without priorities keeps its first member
  x <- Variable("x", t); y <- Variable("y", t); z <- Variable("z", t)
  al <- Parameter("al"); be <- Parameter("be"); ga <- Parameter("ga"); de <- Parameter("de")
  lv <- mtkcompile(System(list(D(x) ~ al * x - be * x * y, D(y) ~ -ga * y + de * x * y, z ~ x + y), t, name = "lv"))
  expect_equal(fmt(unknowns(lv)), c("y(t)", "x(t)"))
  expect_equal(fmt(observables(lv)), "z(t)")
  lin <- mtkcompile(System(list(D(x) ~ -x + y, 0 ~ y - 2 * x), t, name = "lin"))
  expect_equal(fmt(unknowns(lin)), "x(t)")
  expect_equal(fmt(observables(lin)), "y(t)")
  # the feedback loop of the linearize docstring, compiled with the input f.u
  plant <- System(list(D(x) ~ -x + Variable("u", t), y ~ x), t, name = "p")
  filt <- System(list(D(x) ~ -2 * x + Variable("u", t, input = TRUE), y ~ x), t, name = "f")
  ctrl <- System(list(Variable("u", t) ~ Parameter("kp", default = 1) * (Variable("r", t) - y)), t, name = "c")
  cl <- System(list(filt$y ~ ctrl$r, ctrl$u ~ plant$u, plant$y ~ ctrl$y), t, systems = list(filt, ctrl, plant), name = "cl")
  n <- toggle_namespacing(cl, FALSE)
  cs <- mtkcompile(cl, inputs = list(n$f$u), outputs = list(n$p$x))
  expect_equal(fmt(unknowns(cs)), c("c₊y(t)", "c₊r(t)"))
  expect_equal(fmt(observables(cs)), c("c₊u(t)", "f₊y(t)", "f₊x(t)", "p₊y(t)", "p₊x(t)", "p₊u(t)"))
})

test_that("the component RC circuit has the unknown and observed set of ModelingToolkit", {
  st <- julia_fixture("rc_structure.csv")
  skip_if(is.null(st))
  skip_if_not(exists("connector") && exists("connect"))
  Pin <- connector(function(name) {
    System(list(), t, list(Variable("v", t), Variable("i", t, connect = Flow)), list(), name = name)
  })
  OnePort <- function(name) {
    p <- Pin(name = "p"); n <- Pin(name = "n")
    v <- Variable("v", t); i <- Variable("i", t)
    System(list(v ~ p$v - n$v, 0 ~ p$i + n$i, i ~ p$i), t, list(v, i), list(), systems = list(p, n), name = name)
  }
  part <- function(name, par, eq) {
    op <- OnePort("oneport"); u <- toggle_namespacing(op, FALSE)
    q <- Parameter(par, default = 1)
    extend(System(list(eq(u, q)), t, list(), list(q), name = name), op)
  }
  resistor <- part("resistor", "R", function(u, R) u$v ~ u$i * R)
  capacitor <- part("capacitor", "C", function(u, C) D(u$v) ~ u$i / C)
  source <- part("source", "V", function(u, V) V ~ u$v)
  gp <- Pin(name = "g")
  ground <- System(list(gp$v ~ 0), t, list(), list(), systems = list(gp), name = "ground")
  model <- System(list(connect(source$p, resistor$p), connect(resistor$n, capacitor$p),
                       connect(capacitor$n, source$n), connect(capacitor$n, ground$g)),
                  t, systems = list(resistor, capacitor, source, ground), name = "rc_model")
  cs <- mtkcompile(model)
  expect_equal(fmt(unknowns(cs)), julia_names(st, "unknown"))
  expect_setequal(fmt(observables(cs)), julia_names(st, "observed"))
  expect_setequal(fmt(parameters(cs)), julia_names(st, "parameter"))
})

test_that("a flat RC circuit with aliases reduces to the capacitor voltage", {
  # one storage element: ModelingToolkit reduces the RC circuit of
  # docs/src/tutorials/acausal_components.md to the single unknown capacitor.v
  expect_equal(fmt(unknowns(rc)), "vc(t)")
  expect_equal(fmt(lapply(equations(rc), `[[`, "lhs")), "Derivative(vc(t), t)")
  expect_setequal(lhs(observed(rc)), c("vs(t)", "is(t)", "vr(t)", "ir(t)", "ic(t)", "va(t)", "vb(t)", "vg(t)"))
  full <- equations(tearing_substitution(rc))[[1L]]
  vc <- rc$vc
  expect_true(.is_zero(symcompR::simplify(full$rhs - (rc$V - vc) / (rc$R * rc$C))))
  sol <- solve(ODEProblem(rc, list(vc ~ 0), c(0, 4)), Tsit5(), reltol = 1e-9, abstol = 1e-11)
  expect_equal(tail(sol[vc], 1), 1 - exp(-2), tolerance = 1e-7)
  expect_equal(tail(sol[rc$ir], 1), exp(-2) / 2, tolerance = 1e-7)
})

test_that("alias elimination removes redundant integer linear equations", {
  # ten equations for nine unknowns: only the Bareiss pass can balance this
  red <- alias_elimination(rc_redundant)
  expect_length(equations(red), 9L)
  cs <- mtkcompile(rc_redundant)
  expect_equal(fmt(unknowns(cs)), "vc(t)")
  expect_length(observed(cs), 8L)
})

test_that("an algebraic loop is torn into one residual", {
  expect_equal(fmt(unknowns(loop_c)), c("h(t)", "u1(t)"))
  expect_equal(is_diff(loop_c), c(TRUE, FALSE))
  expect_setequal(lhs(observed(loop_c)), c("u2(t)", "u3(t)"))
  expect_equal(diag(calculate_massmatrix(loop_c)), c(1, 0))
  modia <- mtkcompile(loop, tearing_alg = "modia")
  expect_length(unknowns(modia), 2L)
  expect_equal(sum(is_diff(modia)), 1L)
  skip_if_not(exists("Rodas5P"))
  n <- toggle_namespacing(loop, FALSE)
  prob <- ODEProblem(loop_c, list(n$h ~ 0.2), c(0, 2))
  sol <- solve(prob, Rodas5P(), reltol = 1e-8, abstol = 1e-10)
  expect_equal(sol$retcode, "Success")
  expect_lt(max(abs(sol[n$u1 - sin(n$u2) - n$h])), 1e-6)
  expect_lt(max(abs(sol[n$u2 - cos(n$u1) + n$p])), 1e-8)
})

test_that("the Cartesian pendulum is reduced from index 3 by dummy derivatives", {
  x <- pend_vars$x; y <- pend_vars$y; L <- pend_vars$L
  # two differential equations for the prioritized coordinate, three
  # algebraic unknowns: the other coordinate, its velocity and the tension
  # (the order ModelingToolkit 11.46.1 gives with the same priority)
  expect_equal(fmt(unknowns(pend_c)), c("x(t)", "y(t)", "xˍt(t)", "yˍt(t)", "T(t)"))
  expect_equal(lhs(equations(pend_c))[is_diff(pend_c)], c("Derivative(x(t), t)", "Derivative(xˍt(t), t)"))
  alg <- lapply(equations(pend_c)[!is_diff(pend_c)], `[[`, "rhs")
  expect_length(alg, 3L)
  # the position constraint and its first derivative survive as residuals
  expect_true(any(vapply(alg, function(e) .is_zero(symcompR::expand(e - (x^2 + y^2 - L^2))), TRUE)))
  dx <- unknowns(pend_c)[[match("xˍt(t)", fmt(unknowns(pend_c)))]]
  dy <- unknowns(pend_c)[[match("yˍt(t)", fmt(unknowns(pend_c)))]]
  expect_true(any(vapply(alg, function(e) .is_zero(symcompR::expand(e - (2 * x * dx + 2 * y * dy))), TRUE)))
  expect_equal(lhs(observed(pend_c)), "yˍtt(t)")
  expect_length(.f(pend_c, "schedule")$var_sccs, 5L)
  skip_if_not(exists("Rodas5P"))
  th0 <- 0.5
  prob <- ODEProblem(pend_c, list(x ~ sin(th0), D(x) ~ 0, y ~ -cos(th0), D(y) ~ 0, pend_vars$T ~ -9.81 * cos(th0)),
                     c(0, 1))
  sol <- solve(prob, Rodas5P(), reltol = 1e-8, abstol = 1e-8)
  expect_equal(sol$retcode, "Success")
  expect_lt(max(abs(sol[x]^2 + sol[y]^2 - 1)), 1e-6)
  # the angle obeys theta'' = -(g / L) sin(theta); classical Runge-Kutta
  f <- function(u) c(u[2], -9.81 * sin(u[1]))
  u <- c(th0, 0); hh <- 1e-3
  for (i in seq_len(1000)) {
    k1 <- f(u); k2 <- f(u + hh / 2 * k1); k3 <- f(u + hh / 2 * k2); k4 <- f(u + hh * k3)
    u <- u + hh / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  }
  expect_equal(tail(sol[x], 1), sin(u[1]), tolerance = 1e-5)
  expect_equal(tail(sol[y], 1), -cos(u[1]), tolerance = 1e-5)
})

test_that("index one form without dummy derivatives keeps four states", {
  cs <- mtkcompile(pend, dummy_derivative = FALSE)
  expect_equal(sum(is_diff(cs)), 4L)
  expect_length(unknowns(cs), 5L)
})

test_that("structurally singular and unbalanced systems give clear errors", {
  variables("x(t)", "u1(t)", "u2(t)", "a(t)", "b(t)"); parameters("y1", "y2")
  # Example 5 of Pantelides' paper, test/reduction.jl:226-236
  p5 <- System(list(0 ~ x + sin(u1 + u2), D(x) ~ x + y1, cos(x) ~ sin(y2)), t, list(x, u1, u2), list(y1, y2), name = "p5")
  err <- tryCatch(mtkcompile(p5), error = identity)
  expect_s3_class(err, "InvalidSystemException")
  expect_match(conditionMessage(err), "structurally singular")
  expect_match(conditionMessage(err), "u2(t)", fixed = TRUE)
  err <- tryCatch(mtkcompile(System(list(D(a) ~ b, 0 ~ a + 1, 0 ~ a - 2), t, name = "over")), error = identity)
  expect_s3_class(err, "ExtraEquationsSystemException")
  expect_match(conditionMessage(err), "2 highest order derivative variables and 3 equations")
  expect_match(conditionMessage(err), "0 ~ a(t) - 2", fixed = TRUE)
  err <- tryCatch(mtkcompile(System(list(D(a) ~ a + b), t, name = "under")), error = identity)
  expect_s3_class(err, "ExtraVariablesSystemException")
  expect_match(conditionMessage(err), "More variables than equations")
  expect_match(conditionMessage(err), "b(t)", fixed = TRUE)
  # without the requirement, the extra variables stay unknowns
  cs <- mtkcompile(System(list(D(a) ~ a + b), t, name = "under"), fully_determined = FALSE)
  expect_setequal(fmt(unknowns(cs)), c("a(t)", "b(t)"))
  expect_error(mtkcompile(System(list(D(a) ~ -a), t, name = "io"), inputs = list(b)), "do not occur")
  # an equation without unknowns that cannot hold
  err <- tryCatch(mtkcompile(System(list(D(a) ~ -a, Equation(0, 1)), t, name = "never")), error = identity)
  expect_s3_class(err, "InvalidSystemException")
  expect_match(conditionMessage(err), "contains no unknowns")
  # conservative solving keeps 0 ~ 2 a + b intact
  cs <- mtkcompile(System(list(D(a) ~ -a, 0 ~ 2 * a + b), t, name = "cons"), conservative = TRUE)
  ob <- observed(cs)[[1L]]
  expect_equal(format(ob), "b(t) ~ -2*a(t)")
})

test_that("alias targets follow state priorities and irreducibility", {
  variables("x(t)", "y(t)", "z(t)"); parameters("p", "h")
  # test/structural_transformation/tearing.jl:171-182
  xp <- Variable("x", t, state_priority = 1)
  cs <- mtkcompile(System(list(D(xp) ~ z * h, 0 ~ xp - y, 0 ~ sin(z) + y - p * t), t, name = "dae"))
  expect_setequal(fmt(unknowns(cs)), c("x(t)", "z(t)"))
  expect_setequal(fmt(equations(cs)), c("Derivative(x(t), t) ~ h*z(t)", "0 ~ -p*t + x(t) + sin(z(t))"))
  # test/reduction.jl:384-388: a differential variable is swapped for an
  # alias of higher priority
  variables("x(t)")
  cs <- mtkcompile(System(list(D(x) ~ 2 * x, y ~ x), t, name = "swap", state_priorities = list(y ~ 10)))
  expect_equal(fmt(unknowns(cs)), "y(t)")
  expect_equal(fmt(equations(cs)), "Derivative(y(t), t) ~ 2*y(t)")
  # test/reduction.jl:420-439: irreducible variables survive alias groups
  e <- Variable("e", t, irreducible = TRUE); cc <- Variable("c", t, irreducible = TRUE)
  d <- Variable("d", t, irreducible = TRUE)
  cs <- mtkcompile(System(list(D(x) ~ x, D(cc) ~ -cc, e ~ x, cc ~ d, y ~ cc), t, name = "irr"))
  expect_setequal(fmt(unknowns(cs)), c("e(t)", "c(t)", "d(t)"))
  expect_setequal(lhs(observed(cs)), c("x(t)", "y(t)"))
  # contradictory aliases force both variables to zero
  cs <- mtkcompile(System(list(x ~ y, x ~ -y, D(z) ~ x + y - z), t, name = "conflict"))
  expect_equal(fmt(unknowns(cs)), "z(t)")
  expect_setequal(fmt(observed(cs)), c("x(t) ~ 0", "y(t) ~ 0"))
})

test_that("reductions of the ModelingToolkit test suite are reproduced", {
  variables("x(t)", "y(t)", "z(t)", "a(t)", "u(t)"); parameters("sigma", "rho", "beta")
  # test/reduction.jl:38-59
  cs <- mtkcompile(System(list(D(x) ~ sigma * (y - x), D(y) ~ x * (rho - z) - y + beta, 0 ~ z - x + y,
                               0 ~ a + z, u ~ z + a), t, name = "lorenz1"))
  expect_setequal(fmt(unknowns(cs)), c("x(t)", "y(t)"))
  expect_length(parameters(cs), 3L)
  # test/reduction.jl:377-382
  so <- System(list(D(D(x)) ~ -x), t, list(x), list(), name = "so")
  ae <- alias_elimination(so)
  expect_length(equations(ae), 1L); expect_length(unknowns(ae), 1L)
  cs <- mtkcompile(so)
  expect_equal(fmt(unknowns(cs)), c("x(t)", "xˍt(t)"))
  expect_equal(fmt(equations(cs)), c("Derivative(x(t), t) ~ xˍt(t)", "Derivative(xˍt(t), t) ~ -x(t)"))
  # issue 963, test/reduction.jl:244-272; the unknown and the order of the
  # observed variables are those of ModelingToolkit 11.46.1
  variables("v47(t)", "v57(t)", "v66(t)", "v25(t)", "i74(t)", "i75(t)", "i64(t)", "i71(t)", "v1(t)", "v2(t)")
  cs <- mtkcompile(System(list(v47 ~ v1, v47 ~ sin(10 * t), v57 ~ v1 - v2, v57 ~ 10 * i64, v66 ~ v2,
                               v66 ~ 5 * i74, v25 ~ v2, i75 ~ 0.005 * D(v25), 0 ~ i74 + i75 - i64,
                               0 ~ i64 + i71), t, name = "i963"))
  expect_length(equations(cs), 1L)
  expect_equal(fmt(unknowns(cs)), "v25(t)")
  expect_equal(fmt(observables(cs)), paste0(c("v1", "i74", "v2", "v66", "v57", "v47", "i64", "i75", "i71"), "(t)"))
  vv <- unknowns(cs)[[1L]]
  eq <- equations(tearing_substitution(cs))[[1L]]
  expect_equal(.key(eq$lhs), .key(D(vv)))
  expect_equal(as.numeric(diff(eq$rhs, vv)), -60)
})

test_that("zero variables are eliminated and integrated on request", {
  variables("x(t)", "y(t)")
  # test/reduction.jl:367-375
  sys <- System(list(x ~ 0, D(x) ~ x + y), t, list(x, y), list(), name = "zero")
  cs <- mtkcompile(sys, eliminate_mm_zeros = TRUE)
  expect_length(equations(cs), 0L)
  expect_setequal(fmt(observed(cs)), c("x(t) ~ 0", "xˍt(t) ~ 0", "y(t) ~ 0"))
  # D(x) ~ 0 integrates to the constant x#0, a parameter solved at initialization
  cs <- mtkcompile(System(list(D(x) ~ 0, D(y) ~ x), t, name = "int"), eliminate_mm_zeros = TRUE)
  expect_equal(fmt(unknowns(cs)), "y(t)")
  expect_true("x(t) ~ x#0" %in% fmt(observed(cs)))
  expect_true("x#0" %in% fmt(parameters(cs)))
  expect_length(.f(cs, "analytically_integrated"), 1L)
  # the default keeps such states (see R/structural-compile.R)
  expect_setequal(fmt(unknowns(mtkcompile(System(list(D(x) ~ 0, D(y) ~ x), t, name = "int")))), c("x(t)", "y(t)"))
  # the constants of a third order integration are the derivatives at t = 0
  cs <- mtkcompile(System(list(D(D(D(x))) ~ 0, D(y) ~ x), t, name = "int3"), eliminate_mm_zeros = TRUE)
  xo <- Filter(function(e) format(e$lhs) == "x(t)", observed(cs))[[1L]]
  expect_true(.is_zero(symcompR::expand(xo$rhs - (cs[["xˍtt#0"]] * t^2 / 2 + cs[["xˍt#0"]] * t + cs[["x#0"]]))))
  # test/reduction.jl:898-905: an irreducible zero derivative is retained
  yi <- Variable("y", t, irreducible = TRUE)
  cs <- mtkcompile(System(list(D(x) ~ x - x * yi, D(yi) ~ 0), t, name = "i4764"), eliminate_mm_zeros = TRUE)
  expect_true("y(t)" %in% fmt(unknowns(cs)))
  expect_length(.f(cs, "analytically_integrated"), 0L)
})

test_that("TearingState records incidence and derivative chains", {
  variables("x(t)", "y(t)"); parameters("k")
  ts <- TearingState(System(list(D(D(x)) ~ -k * y, 0 ~ x - y), t, name = "ts"))
  expect_equal(fmt(ts$fullvars), c("Derivative(x(t), (t, 2))", "Derivative(x(t), t)", "x(t)", "y(t)"))
  expect_equal(ts$structure$var_to_diff, c(NA, 1L, 2L, NA))
  expect_equal(ts$structure$diff_to_var, c(2L, 3L, NA, NA))
  # the alias equation sorts before the differential equation
  expect_equal(ts$structure$graph$fadjlist, list(c(3L, 4L), c(1L, 4L)))
  expect_equal(highest_order_variable_mask(ts)(1:4), c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(lowest_order_variable_mask(ts)(1:4), c(FALSE, FALSE, TRUE, TRUE))
  find_solvables(ts)
  # y has the parameter coefficient -k, which may be divided by
  expect_equal(ts$structure$solvable_graph$fadjlist, list(c(3L, 4L), c(1L, 4L)))
  expect_error(TearingState(System(list(D(x) ~ y), t, list(x), list(), name = "bad")), "neither an unknown nor a parameter")
})

test_that("index lowering, tearing and the incidence utilities work on their own", {
  variables("x(t)", "y(t)", "z(t)"); parameters("k")
  idx2 <- System(list(D(x) ~ y, 0 ~ x - sin(t)), t, name = "idx")
  low <- dae_index_lowering(idx2)
  expect_setequal(fmt(equations(low)), c("Derivative(x(t), t) ~ y(t)", "0 ~ -cos(t) + Derivative(x(t), t)"))
  dd <- dummy_derivative(idx2)
  expect_length(unknowns(dd), 0L)
  expect_setequal(fmt(observed(dd)), c("x(t) ~ sin(t)", "xˍt(t) ~ cos(t)", "y(t) ~ xˍt(t)"))
  torn <- tearing(System(list(D(x) ~ -k * y, 0 ~ x - y + z, 0 ~ z - sin(x)), t, name = "torn"))
  expect_equal(fmt(unknowns(torn)), "x(t)")
  expect_equal(fmt(observed(torn)), c("z(t) ~ sin(x(t))", "y(t) ~ x(t) + z(t)"))
  expect_equal(sorted_incidence_matrix(loop_c), rbind(c(1L, 1L, 1L), c(1L, 1L, 0L)))
  expect_error(sorted_incidence_matrix(loop), "simplified")
  boi <- but_ordered_incidence(TearingState(loop))
  expect_equal(dim(boi[[1L]]), c(4L, 4L))
  expect_equal(boi[[2L]], c(1L, 3L, 4L, 5L))
  mp <- map_variables_to_equations(loop_c)
  expect_setequal(names(mp), c("h(t)", "u1(t)", "u2(t)", "u3(t)"))
  expect_equal(format(mp[["u1(t)"]]), "0 ~ -h(t) + u1(t) - sin(u2(t))")
  expect_equal(format(mp[["h(t)"]]$lhs), "Derivative(h(t), t)")
  expect_error(map_variables_to_equations(loop), "simplified")
  io <- mtkcompile(System(list(D(x) ~ -x + z, y ~ 2 * x), t, name = "io"), inputs = list(z), outputs = list(y))
  expect_true("z(t)" %in% fmt(parameters(io)))
  expect_equal(fmt(unknowns(io)), "x(t)")
})

test_that("graphs, matchings and the Bareiss factorization are exact", {
  g <- BipartiteGraph(fadjlist = list(c(1L, 2L), 2L, c(2L, 3L)), ndsts = 3L)
  expect_equal(g$badjlist, list(1L, c(1L, 2L, 3L), 3L))
  m <- .maximal_matching(g)
  expect_equal(m, c(1L, 2L, 3L))
  # x2 is needed by the equations of x1 and x3, so it comes first
  expect_equal(.find_var_sccs(g, m), list(2L, 3L, 1L))
  loop_g <- BipartiteGraph(fadjlist = list(c(1L, 2L), c(1L, 2L)), ndsts = 2L)
  expect_equal(.find_var_sccs(loop_g, .maximal_matching(loop_g)), list(c(1L, 2L)))
  # rank 2 integer matrix: the third column is the sum of the first two
  b <- .bareiss_rank(rbind(c(1, 2, 3), c(2, 1, 3), c(3, 3, 6)))
  expect_equal(b$rank, 2L)
  # x + y = 0 and x + y + z = 0 make z zero and leave one free variable
  s <- .ss_structure(BipartiteGraph(fadjlist = list(c(1L, 2L), c(1L, 2L, 3L)), ndsts = 3L), 3L)
  r <- .do_bareiss(list(nzrows = 1:2, row_cols = list(c(1L, 2L), c(1L, 2L, 3L)), row_vals = list(c(1, 1), c(1, 1, 1))),
                   rep(TRUE, 3), rep(TRUE, 3), s$canonical_ranks)
  expect_equal(r$rank3, 2L)
  expect_equal(r$M$row_cols[[2L]], 3L)
  expect_error(.bareiss_guard(2^60), "overflow")
  variables("x(t)", "y(t)"); parameters("k")
  dep <- System(list(D(x) ~ -k * y, D(y) ~ x), t, name = "dep")
  expect_equal(lapply(equation_dependencies(dep), fmt), list("y(t)", "x(t)"))
  eg <- asgraph(dep); vg <- variable_dependencies(dep)
  expect_equal(eg$fadjlist, list(2L, 1L))
  expect_equal(vg$fadjlist, list(1L, 2L))
  expect_equal(eqeq_dependencies(eg, vg)$fadjlist, list(2L, 1L))
  expect_equal(varvar_dependencies(eg, vg)$fadjlist, list(2L, 1L))
  expect_equal(asdigraph(eg, dep)$fadjlist, list(4L, 3L, integer(), integer()))
})

test_that("components are substituted and subsystems isolated", {
  x <- Variable("x", t); k <- Parameter("k")
  decay <- System(list(D(x) ~ -k * x), t, name = "part")
  growth <- System(list(D(x) ~ k * x), t, name = "part")
  model <- System(list(), t, systems = list(decay), name = "model")
  expect_equal(fmt(equations(substitute_component(model, model$part, growth))),
               "Derivative(part₊x(t), t) ~ part₊k*part₊x(t)")
  expect_error(substitute_component(model, model$part, System(list(D(x) ~ x), t, name = "other")), "same name")
  y <- Variable("y", t)
  expect_error(substitute_component(model, model$part, System(list(D(y) ~ k * y), t, name = "part")),
               "does not contain unknown")
  expect_error(substitute_component(complete(model), model$part, growth), "completed")
  skip_if_not(exists("connect") && exists("AnalysisPoint"))
  # test "basic plant isolation" and "internal connections are preserved" of
  # lib/ModelingToolkitBase/test/analysis_points.jl:1201-1260
  block <- function(name, eq) {
    u <- Variable("u", t, input = TRUE); y <- Variable("y", t, output = TRUE)
    System(list(eq(u, y)), t, name = name)
  }
  P <- block("P", function(u, y) D(y) ~ u - y); C <- block("C", function(u, y) y ~ -u)
  cl <- System(list(connect(C$y, "u", P$u), connect(P$y, "y", C$u)), t, systems = list(P, C), name = "cl")
  iso <- isolate_subsystem(cl, "u", "y")
  expect_equal(vapply(get_systems(iso$sys), nameof, ""), "P")
  expect_length(get_eqs(iso$sys), 0L)
  expect_equal(fmt(iso$input_vars), "P₊u(t)")
  expect_equal(fmt(iso$output_vars), "P₊y(t)")
  expect_equal(fmt(isolate_subsystem(cl, cl$u, cl$y)$input_vars), "P₊u(t)")
  P1 <- block("P1", function(u, y) D(y) ~ u - y); P2 <- block("P2", function(u, y) D(y) ~ u - 2 * y)
  cl2 <- System(list(connect(C$y, "u", P1$u), connect(P1$y, P2$u), connect(P2$y, "y", C$u)), t,
                systems = list(P1, P2, C), name = "cl")
  iso <- isolate_subsystem(cl2, "u", "y")
  expect_setequal(vapply(get_systems(iso$sys), nameof, ""), c("P1", "P2"))
  expect_length(get_eqs(iso$sys), 1L)
  expect_s3_class(get_eqs(iso$sys)[[1L]], "Connection")
  expect_equal(fmt(iso$input_vars), "P1₊u(t)")
  expect_equal(fmt(iso$output_vars), "P2₊y(t)")
})

test_that("array variables scalarize through the structural pipeline", {
  # elementwise arithmetic, %*% and array equations, as in MTK's TearingState
  x <- Variable("s2a_x", t, dims = 3); y <- Variable("s2a_y", t, dims = 3)
  k <- Parameter("s2a_k")
  expect_s3_class(-x, "mtk_array")
  expect_equal(fmt(k * x), paste0("s2a_k*s2a_x[", 1:3, "](t)"))
  expect_equal(fmt(x * k), fmt(k * x))
  expect_equal(format(sum(x)), "s2a_x[1](t) + s2a_x[2](t) + s2a_x[3](t)")
  # vectors follow the orientation rules of base R's %*%
  expect_equal(format(x %*% x), "s2a_x[1](t)^2 + s2a_x[2](t)^2 + s2a_x[3](t)^2")
  expect_equal(dim(x %*% diag(3)), c(1L, 3L))
  expect_length(.as_equations(list(D(x) ~ -x)), 3L)
  variables("s2a_z(t)[1:3]")
  cs <- mtkcompile(System(list(D(s2a_z) ~ -k * s2a_z), t, name = "s2arr"))
  expect_length(unknowns(cs), 3L)
  sol <- solve(ODEProblem(cs, list(s2a_z ~ c(1, 2, 3), k ~ 0.5), c(0, 1)), abstol = 1e-10, reltol = 1e-10)
  expect_equal(unname(sol(1, idxs = s2a_z)), c(1, 2, 3) * exp(-0.5), tolerance = 1e-7)
  # D(x) ~ x, p y ~ x with a 3 x 3 parameter matrix: y = p^-1 x0 e^t
  p <- Parameter("s2a_p", dims = c(3, 3))
  pval <- matrix(c(2, 0.5, 1, 1, 3, 0.2, 0.3, 1, 4), 3, 3)
  cs <- mtkcompile(System(list(D(x) ~ x, p %*% y ~ x), t, name = "s2arr2"))
  expect_length(equations(cs), 5L)
  prob <- ODEProblem(cs, list(x ~ c(1, 2, 3), p ~ pval), c(0, 1), guesses = list(y ~ c(1, 1, 1)))
  sol <- solve(prob, Rodas5P(), abstol = 1e-10, reltol = 1e-10)
  expect_equal(unname(sol(1, idxs = y)), solve(pval, c(1, 2, 3)) * exp(1), tolerance = 1e-7)
})

test_that("linear algebraic blocks are solved inline (inline_linear_sccs)", {
  # test/structural_transformation/tearing.jl:241-275; the observed solve
  # and its 2 x 2 matrix are those printed by ModelingToolkit 11.47
  x <- Variable("s2i_x", t, dims = 3); y <- Variable("s2i_y", t, dims = 3)
  p <- Parameter("s2i_p", dims = c(3, 3))
  sys <- System(list(D(x) ~ x, p %*% y ~ x), t, name = "s2inl")
  s1 <- mtkcompile(sys, inline_linear_sccs = TRUE, analytical_linear_scc_limit = 1)
  s2 <- mtkcompile(sys, inline_linear_sccs = TRUE, analytical_linear_scc_limit = 5)
  expect_length(equations(s1), 3L)
  expect_length(equations(s2), 3L)
  uses_solve <- function(s) any(grepl("inline_linsolve", fmt(lapply(observed(s), `[[`, "rhs"))))
  expect_true(uses_solve(s1))
  expect_false(uses_solve(s2))
  pval <- matrix(c(2, 0.5, 1, 1, 3, 0.2, 0.3, 1, 4), 3, 3)
  ref <- solve(pval, c(1, 2, 3)) * exp(1)
  # y(1) of ModelingToolkit 11.47 with analytical_linear_scc_limit = 1, Tsit5
  expect_equal(ref, c(0.5278945842666167, 1.107077504919326, 1.851383850034579), tolerance = 1e-9)
  for (s in list(s1, s2)) {
    prob <- ODEProblem(s, list(x ~ c(1, 2, 3), p ~ pval), c(0, 1))
    for (alg in list(Tsit5(), Rodas5P())) {
      sol <- solve(prob, alg, abstol = 1e-10, reltol = 1e-10)
      expect_equal(unname(sol(1, idxs = y)), ref, tolerance = 1e-8)
    }
  }
  # a scalar loop: ModelingToolkit gives u1 ~ (-z - sin(t)/b)/(-a - 1/b)
  u1 <- Variable("s2i_u1", t); u2 <- Variable("s2i_u2", t); z <- Variable("s2i_z", t)
  a <- Parameter("s2i_a", default = 2); b <- Parameter("s2i_b", default = 3)
  sl <- System(list(D(z) ~ -z + u1, a * u1 + u2 ~ z, u1 - b * u2 ~ sin(t)), t, name = "s2inls")
  q0 <- mtkcompile(sl)
  q1 <- mtkcompile(sl, inline_linear_sccs = TRUE)
  expect_length(equations(q0), 2L)
  expect_length(equations(q1), 1L)
  u1o <- Filter(function(e) .key(e$lhs) == .key(u1), observed(q1))[[1L]]
  expect_true(.is_zero(symcompR::simplify(u1o$rhs - (-z - sin(t) / b) / (-a - 1 / b))))
  v0 <- solve(ODEProblem(q0, list(z ~ 1), c(0, 2)), Rodas5P(), abstol = 1e-10, reltol = 1e-10)
  v1 <- solve(ODEProblem(q1, list(z ~ 1), c(0, 2)), abstol = 1e-10, reltol = 1e-10)
  expect_equal(v1(2, idxs = u1), v0(2, idxs = u1), tolerance = 1e-7)
})

test_that("nonlinear algebraic blocks become nested solves (inline_nonlinear_sccs)", {
  xn <- Variable("s2n_x", t); yn <- Variable("s2n_y", t, guess = 0.5)
  nl <- System(list(D(xn) ~ -xn + yn, 0 ~ yn^3 + yn - xn), t, name = "s2nl")
  n1 <- mtkcompile(nl, inline_nonlinear_sccs = TRUE)
  expect_length(equations(n1), 1L)
  expect_equal(fmt(observables(n1)), "s2n_y(t)")
  # ModelingToolkit 11.47 solving the DAE with Rodas5P, tolerances 1e-12
  ref <- c(0.7560471618879038, 0.5704323907945071, 0.6078644422996997, 0.4901252219074364)
  for (alg in list(Tsit5(), Rodas5P())) {
    so <- solve(ODEProblem(n1, list(xn ~ 1), c(0, 2)), alg, abstol = 1e-12, reltol = 1e-12)
    expect_equal(c(so(1, idxs = xn), so(1, idxs = yn), so(2, idxs = xn), so(2, idxs = yn)), ref,
                 tolerance = 1e-8)
  }
  # the incidence graph follows the inputs of the nested solve: D(z) ~ y
  # depends on x through y
  xg <- Variable("s2n_xg", t); zg <- Variable("s2n_zg", t); yg <- Variable("s2n_yg", t)
  cg <- mtkcompile(System(list(D(xg) ~ -xg, D(zg) ~ yg, 0 ~ yg^3 + yg - xg), t, name = "s2ng"),
                   inline_nonlinear_sccs = TRUE)
  ts <- .f(cg, "tearing_state")
  row <- which(vapply(equations(cg), function(e) .key(e$lhs) == .key(D(zg)), TRUE))
  expect_true(match(.key(xg), .keys(ts$fullvars)) %in% ts$structure$graph$fadjlist[[row]])
  # a block without inputs is solved at the first call
  xq <- Variable("s2n_xq", t); yq <- Variable("s2n_yq", t)
  cq <- mtkcompile(System(list(D(xq) ~ yq, 0 ~ yq^3 + yq - 2), t, name = "s2nq"), inline_nonlinear_sccs = TRUE)
  sq <- solve(ODEProblem(cq, list(xq ~ 0), c(0, 1)))
  expect_equal(c(sq(1, idxs = xq), sq(1, idxs = yq)), c(1, 1), tolerance = 1e-8)
  # a torn block whose solved equation is folded into the residual (the loop
  # model above); ModelingToolkit 11.47 values at t = 2 with Rodas5P
  l1 <- mtkcompile(loop, inline_nonlinear_sccs = TRUE)
  expect_length(equations(l1), 1L)
  n <- toggle_namespacing(loop, FALSE)
  s1 <- solve(ODEProblem(l1, list(n$h ~ 0.5), c(0, 2)), Tsit5(), abstol = 1e-12, reltol = 1e-12)
  s0 <- solve(ODEProblem(loop_c, list(n$h ~ 0.5), c(0, 2)), Rodas5P(), abstol = 1e-12, reltol = 1e-12)
  vs <- list(n$h, n$u1, n$u2, n$u3)
  expect_equal(vapply(vs, function(v) s1(2, idxs = v), 0), vapply(vs, function(v) s0(2, idxs = v), 0),
               tolerance = 1e-8)
  expect_equal(vapply(vs, function(v) s1(2, idxs = v), 0),
               c(0.9679947818471124, 1.1111407845780275, 0.14363942837727733, 1.2547802129553047),
               tolerance = 1e-8)
})

test_that("IfLifting turns conditions into parameters switched by events", {
  # test/if_lifting.jl:6-52
  x <- Variable("s2l_x", t); y <- Variable("s2l_y", t)
  sys <- System(list(D(x) ~ abs(y), y ~ sin(t)), t, name = "s2abs")
  ss2 <- mtkcompile(sys, additional_passes = list(IfLifting))
  expect_length(equations(ss2), 1L)
  expect_length(parameters(ss2), 1L)
  expect_equal(.type(equations(ss2)[[1L]]$rhs), "Piecewise")
  cv <- parameters(ss2)[[1L]]
  sol <- solve(ODEProblem(ss2, list(x ~ 0), c(0, 5)), Tsit5())
  expect_equal(sum(abs(sol$t - pi) < 1e-8), 2L)
  vals <- sol[cv]
  expect_equal(vals[c(1L, length(vals))], c(0, 1))
  # x = 1 - cos(t) on [0, pi), 3 + cos(t) on [pi, 2 pi)
  expect_equal(sol(pi + 1, idxs = x), 3 + cos(pi + 1), tolerance = 1e-3)
  # test/if_lifting.jl:145-201, nested conditions
  xr <- Variable("s2l_r", t); yr <- Variable("s2l_ry", t)
  st <- Parameter("s2l_st", default = 1); du <- Parameter("s2l_du", default = 1)
  h <- Parameter("s2l_h", default = 1)
  ramp <- System(list(yr ~ IfElse(st < t, IfElse(t < st + du, (t - st) * h / du, h), 0), D(xr) ~ yr), t,
                 name = "s2ramp")
  s2 <- mtkcompile(ramp, additional_passes = list(IfLifting))
  so <- solve(ODEProblem(s2, list(xr ~ 1), c(0, 3)))
  expect_equal(so(0.99, idxs = xr), 1)
  expect_equal(so(1.01, idxs = xr), 1.00005, tolerance = 1e-6)
  expect_equal(so(2, idxs = xr), 1.5, tolerance = 1e-6)
  expect_equal(so(3, idxs = xr), 2.5, tolerance = 1e-6)
  # no_if_lift and time independent conditions are left alone
  pz <- Parameter("s2l_p", default = 1); w <- Variable("s2l_w", t); z <- Variable("s2l_z", t)
  big <- System(list(D(w) ~ IfElse(pz < 3, 1, 2), D(z) ~ no_if_lift(IfElse(w < z, w, z))), t, name = "s2big")
  sb <- mtkcompile(big, additional_passes = list(IfLifting))
  expect_length(parameters(sb), 1L)
  # mod stops the integrator at every multiple of the period
  xm <- Variable("s2l_m", t)
  sm <- mtkcompile(System(list(D(xm) ~ t %% 1), t, name = "s2mod"), additional_passes = list(IfLifting))
  expect_length(parameters(sm), 1L)
  expect_false(any(grepl("ifl_root", fmt(equations(sm)))))
  som <- solve(ODEProblem(sm, list(xm ~ 0), c(0, 3.5)), Tsit5(), reltol = 1e-9, abstol = 1e-9)
  expect_equal(sum(abs(som$t - round(som$t)) < 1e-9 & som$t > 0), 6L)
  expect_equal(som(3.5, idxs = xm), 1.625, tolerance = 1e-7)
  sf <- mtkcompile(System(list(D(xm) ~ t %/% 1), t, name = "s2div"), additional_passes = list(IfLifting))
  sof <- solve(ODEProblem(sf, list(xm ~ 0), c(0, 3.5)), Tsit5(), reltol = 1e-9, abstol = 1e-9)
  expect_equal(sof(3.5, idxs = xm), 4.5, tolerance = 1e-7)
})

test_that("events update the constants of analytically integrated variables", {
  # AffectSystem of ModelingToolkitBase/src/systems/callbacks.jl: x ~ x#0
  x <- Variable("s2z_x", t); y <- Variable("s2z_y", t)
  ev <- SymbolicDiscreteCallback(list(1.0), list(x ~ Pre(x) + 1))
  cs <- mtkcompile(System(list(D(x) ~ 0, D(y) ~ x), t, discrete_events = list(ev), name = "s2mmz"),
                   eliminate_mm_zeros = TRUE)
  expect_equal(fmt(unknowns(cs)), "s2z_y(t)")
  sol <- solve(ODEProblem(cs, list(x ~ 2, y ~ 0), c(0, 2)))
  expect_equal(sol(c(0.5, 1.5), idxs = x), c(2, 3))
  expect_equal(sol(2, idxs = y), 5, tolerance = 1e-8)
  # after D(D(x)) ~ 0 the event moves x#0 and keeps the slope
  cs <- mtkcompile(System(list(D(D(x)) ~ 0, D(y) ~ x), t, discrete_events = list(ev), name = "s2mmz2"),
                   eliminate_mm_zeros = TRUE)
  sol <- solve(ODEProblem(cs, list(x ~ 2, D(x) ~ 0.5, y ~ 0), c(0, 2)))
  expect_equal(sol(c(0.5, 1.5, 2), idxs = x), c(2.25, 3.75, 4))
})

test_that("Initial(x) cannot collide with a parameter named like it", {
  # with Initial(x) named by its printed form, q and Initial(x) were the same
  # symbol in the initialization system and y came out as 9
  x <- Variable("s2c_x", t); y <- Variable("s2c_y", t)
  q <- Parameter("Initial(s2c_x(t))", default = 3)
  cs <- mtkcompile(System(list(D(x) ~ -q * x, D(y) ~ -y), t, name = "s2col",
                          initialization_eqs = list(y ~ q * x)))
  prob <- ODEProblem(cs, list(x ~ 2), c(0, 1), guesses = list(y ~ 1))
  expect_equal(unname(prob$u0[c("s2c_x(t)", "s2c_y(t)")]), c(2, 6))
})

test_that("the implicit stage system reproduces an implicit Euler step (nlstep)", {
  # src/systems/solver_nlprob.jl: 0 = outer + g1 f(g2 z + inner, p, c) - g3 M z
  x <- Variable("s2s_x", t); y <- Variable("s2s_y", t); a <- Parameter("s2s_a", default = 2)
  sys <- mtkcompile(System(list(D(x) ~ -a * x^3 + y, D(y) ~ -y), t, name = "s2s"))
  prob <- ODEProblem(sys, list(x ~ 1, y ~ 0.5), c(0, 1), nlstep = TRUE)
  nd <- prob$f$nlstep_data
  expect_s3_class(nd, "ODENLStepData")
  h <- 0.1
  np <- nd$set_gamma_c(nd$nlprob, c(h, 1, 1, h))
  np <- nd$set_outer_tmp(np, prob$u0)
  np <- nd$set_inner_tmp(np, c(0, 0))
  z <- nd$nlprobmap(solve(np)$u, np$p)
  ref <- .newton_solve(function(zz) prob$u0 + h * prob$f$f(zz, prob$p, h) - zz, prob$u0)$x
  expect_equal(unname(z), unname(ref), tolerance = 1e-9)
})
