# Connectors, components and connection expansion. The electrical library
# follows docs/src/tutorials/acausal_components.md of ModelingToolkit.

.rc_library <- function() {
  t <- t_nounits; D <- D_nounits
  Pin <- connector(function(name) {
    v <- Variable("v", t); i <- Variable("i", t, connect = Flow)
    System(list(), t, list(v, i), list(), name = name)
  })
  Ground <- function(name) {
    g <- Pin(name = "g")
    System(list(g$v ~ 0), t, list(), list(), systems = list(g), name = name)
  }
  OnePort <- function(name) {
    p <- Pin(name = "p"); n <- Pin(name = "n")
    v <- Variable("v", t); i <- Variable("i", t)
    System(list(v ~ p$v - n$v, 0 ~ p$i + n$i, i ~ p$i), t, list(v, i), list(),
           systems = list(p, n), name = name)
  }
  Resistor <- function(name, R = 1) {
    op <- OnePort("oneport"); u <- toggle_namespacing(op, FALSE)
    Rp <- Parameter("R", default = R)
    extend(System(list(u$v ~ u$i * Rp), t, list(), list(Rp), name = name), op)
  }
  Capacitor <- function(name, C = 1) {
    op <- OnePort("oneport"); u <- toggle_namespacing(op, FALSE)
    Cp <- Parameter("C", default = C)
    extend(System(list(D(u$v) ~ u$i / Cp), t, list(), list(Cp), name = name), op)
  }
  ConstantVoltage <- function(name, V = 1) {
    op <- OnePort("oneport"); u <- toggle_namespacing(op, FALSE)
    Vp <- Parameter("V", default = V)
    extend(System(list(Vp ~ u$v), t, list(), list(Vp), name = name), op)
  }
  list(Pin = Pin, Ground = Ground, OnePort = OnePort, Resistor = Resistor,
       Capacitor = Capacitor, ConstantVoltage = ConstantVoltage)
}

.rc_model <- function(lib, R = 2, C = 1, V = 1) {
  resistor <- lib$Resistor("resistor", R = R)
  capacitor <- lib$Capacitor("capacitor", C = C)
  source <- lib$ConstantVoltage("source", V = V)
  ground <- lib$Ground("ground")
  System(list(connect(source$p, resistor$p), connect(resistor$n, capacitor$p),
              connect(capacitor$n, source$n), connect(capacitor$n, ground$g)),
         t_nounits, systems = list(resistor, capacitor, source, ground), name = "rc_model")
}

# Simulate a compiled index 1 DAE through its state space reduction: the
# algebraic unknowns are solved by the linearization function at each call.
.reduced_rhs <- function(lin_fun) {
  st <- attr(lin_fun, "state")
  d <- st$diff_idxs
  u0 <- st$u0
  function(x, p, t) {
    u <- u0
    u[d] <- x
    b <- lin_fun(u, st$p0, t)
    st$F(b$x, b$p, t)[d]
  }
}

test_that("the RC circuit expands, compiles and matches the analytic response", {
  lib <- .rc_library()
  rc <- .rc_model(lib)
  ex <- expand_connections(rc)
  eqs <- equations(ex)
  expect_length(eqs, 20L)
  expect_length(unknowns(ex), 20L)
  fmt <- vapply(eqs, format, "")
  expect_true("0 ~ capacitor\u208An\u208Ai(t) + ground\u208Ag\u208Ai(t) + source\u208An\u208Ai(t)" %in% fmt)
  expect_true("source\u208Ap\u208Av(t) ~ resistor\u208Ap\u208Av(t)" %in% fmt)
  expect_false(any(vapply(eqs, inherits, TRUE, "Connection")))

  cs <- mtkcompile(rc)
  expect_true(iscomplete(cs))
  rn <- toggle_namespacing(rc, FALSE)
  vc <- rn$capacitor$v
  expect_true(.contains_key(unknowns(cs), vc))
  lf <- linearization_function(rc, list(), list(vc), op = list(vc ~ 0))
  sol <- solve(ODEProblem(.reduced_rhs(lf$lin_fun), 0, c(0, 10)), Tsit5(), reltol = 1e-9, abstol = 1e-11)
  tt <- c(1, 2, 5, 10)
  vals <- vapply(tt, function(ti) sol(ti)[[1L]], 0)
  expect_equal(vals, 1 - exp(-tt / 2), tolerance = 1e-7)
  if (all(vapply(equations(cs), function(e) .is_derivative(e$lhs), TRUE))) {
    prob <- ODEProblem(cs, list(cs$capacitor$v ~ 0), c(0, 10))
    s2 <- solve(prob, Tsit5(), reltol = 1e-9, abstol = 1e-11)
    expect_equal(tail(s2[cs$capacitor$v], 1), 1 - exp(-5), tolerance = 1e-6)
  }
  lin <- linearize(rc, list(), list(vc, rn$resistor$i), op = list(vc ~ 0.5))
  expect_equal(lin$A, matrix(-0.5), tolerance = 1e-12)
  expect_equal(lin$C[, 1], c(1, -0.5), tolerance = 1e-12)
  ri <- build_explicit_observed_function(lin$sys, rn$resistor$i)
  expect_equal(unname(ri(lin$extras$x, lin$extras$p, 0)), 0.25, tolerance = 1e-10)
})

test_that("the RC circuit matches ModelingToolkit (tools/julia/reference.jl, section 7)", {
  fx <- test_path("fixtures", "julia", "rc.csv")
  skip_if_not(file.exists(fx), "Julia RC fixture not generated")
  ref <- utils::read.csv(fx, check.names = FALSE, encoding = "UTF-8")
  lib <- .rc_library()
  rc <- .rc_model(lib, R = 1, C = 1, V = 1)
  rn <- toggle_namespacing(rc, FALSE)
  vc <- rn$capacitor$v
  lf <- linearization_function(rc, list(), list(vc), op = list(vc ~ 0))
  sol <- solve(ODEProblem(.reduced_rhs(lf$lin_fun), 0, c(0, 10)), Tsit5(), reltol = 1e-10, abstol = 1e-12)
  st <- attr(lf$lin_fun, "state")
  ri <- build_explicit_observed_function(lf$sys, rn$resistor$i)
  vals <- t(vapply(ref$t, function(ti) {
    u <- st$u0
    u[st$diff_idxs] <- sol(ti)[[1L]]
    b <- lf$lin_fun(u, st$p0, ti)
    c(b$x[[st$diff_idxs]], unname(ri(b$x, b$p, ti)))
  }, c(0, 0)))
  expect_equal(vals[, 1], ref[["capacitor\u208Av(t)"]], tolerance = 1e-7)
  expect_equal(vals[, 2], ref[["resistor\u208Ai(t)"]], tolerance = 1e-7)
  fs <- test_path("fixtures", "julia", "rc_structure.csv")
  skip_if_not(file.exists(fs), "Julia RC structure fixture not generated")
  jl <- utils::read.csv(fs, encoding = "UTF-8")
  # the expanded model has the variables of ModelingToolkit's compiled system
  expect_setequal(vapply(unknowns(expand_connections(rc)), format, ""), jl$name[jl$kind != "parameter"])
  cs <- mtkcompile(rc)
  if (all(vapply(equations(cs), function(e) .is_derivative(e$lhs), TRUE)))
    expect_setequal(vapply(unknowns(cs), format, ""), jl$name[jl$kind == "unknown"])
})

test_that("mtkmodel builds the same circuit from declarations", {
  t <- t_nounits; D <- D_nounits
  Pin <- connector("Pin", variables = list("v(t)", list("i(t)", connect = Flow)))
  Ground <- mtkmodel("Ground", components = list(g = quote(Pin())), equations = quote(list(g$v ~ 0)))
  OnePort <- mtkmodel("OnePort", components = list(p = quote(Pin()), n = quote(Pin())),
                      variables = list("v(t)", "i(t)"),
                      equations = quote(list(v ~ p$v - n$v, 0 ~ p$i + n$i, i ~ p$i)))
  Resistor <- mtkmodel("Resistor", extend = list(oneport = quote(OnePort())),
                       parameters = list("R = 1"), equations = quote(list(v ~ i * R)))
  Capacitor <- mtkmodel("Capacitor", extend = list(oneport = quote(OnePort())),
                        parameters = list("C = 1"), equations = quote(list(D(v) ~ i / C)))
  ConstantVoltage <- mtkmodel("ConstantVoltage", extend = list(oneport = quote(OnePort())),
                              parameters = list("V = 1"), equations = quote(list(V ~ v)))
  RC <- mtkmodel("RC", structural_parameters = list(R0 = 2),
    components = list(resistor = quote(Resistor(R = R0)), capacitor = quote(Capacitor()),
                      source = quote(ConstantVoltage()), ground = quote(Ground())),
    equations = quote(list(connect(source$p, resistor$p), connect(resistor$n, capacitor$p),
                           connect(capacitor$n, source$n), connect(capacitor$n, ground$g))))
  expect_s3_class(RC, "Model")
  expect_identical(connector_type(Pin(name = "p")), "RegularConnector")
  rc <- RC(name = "rc", R0 = 4, capacitor__C = 0.5)
  expect_equal(smap_get(initial_conditions(rc), toggle_namespacing(rc, FALSE)$resistor$R), 4)
  expect_equal(smap_get(initial_conditions(rc), toggle_namespacing(rc, FALSE)$capacitor$C), 0.5)
  r1 <- Resistor(name = "r1", R = 3); r2 <- Resistor(name = "r2")
  expect_equal(smap_get(initial_conditions(r1), toggle_namespacing(r1, FALSE)$R), 3)
  expect_equal(smap_get(initial_conditions(r2), toggle_namespacing(r2, FALSE)$R), 1)
  expect_length(equations(expand_connections(rc)), 20L)
  # arguments reach declarations inherited through several extends
  Base <- mtkmodel("Base", parameters = list("k = 1"), variables = list("x(t) = 1"),
                   equations = quote(list(D(x) ~ -k * x)))
  Mid <- mtkmodel("Mid", extend = list(base = quote(Base())))
  Top <- mtkmodel("Top", extend = list(mid = quote(Mid())))
  top <- Top(name = "top", k = 3)
  expect_equal(smap_get(initial_conditions(top), toggle_namespacing(top, FALSE)$k), 3)
  vc <- toggle_namespacing(rc, FALSE)$capacitor$v
  lin <- linearize(rc, list(), list(vc), op = list(vc ~ 0))
  expect_equal(lin$A, matrix(-1 / (4 * 0.5)), tolerance = 1e-12)
})

test_that("causal variable connections root at the output", {
  t <- t_nounits
  y <- Variable("y", t, output = TRUE); u <- Variable("u", t, input = TRUE)
  src <- System(list(y ~ sin(t)), t, name = "src")
  k1 <- System(list(), t, list(u), list(), name = "k1")
  k2 <- System(list(), t, list(u), list(), name = "k2")
  sys <- System(list(connect(src$y, k1$u, k2$u)), t, systems = list(src, k1, k2), name = "s")
  fmt <- vapply(equations(expand_connections(sys)), format, "")
  expect_true(all(c("k1\u208Au(t) ~ src\u208Ay(t)", "k2\u208Au(t) ~ src\u208Ay(t)") %in% fmt))
  # top level variables are outside connectors: two outer inputs have no source
  u2 <- Variable("u2", t, input = TRUE)
  bad <- System(list(connect(y, u, u2)), t, name = "b")
  expect_error(expand_connections(bad), "two outer input")
  z <- Variable("z", t)
  expect_error(connect(y, z), "non-causal")
})

test_that("outer connectors enter flow balances with a negative sign", {
  lib <- .rc_library()
  t <- t_nounits
  TwoResistors <- function(name) {
    p <- lib$Pin(name = "p"); n <- lib$Pin(name = "n")
    r1 <- lib$Resistor("r1"); r2 <- lib$Resistor("r2")
    System(list(connect(p, r1$p), connect(r1$n, r2$p), connect(r2$n, n)), t,
           systems = list(p, n, r1, r2), name = name)
  }
  sys <- System(list(), t, systems = list(TwoResistors("tr")), name = "top")
  fmt <- vapply(equations(expand_connections(sys)), format, "")
  expect_true("0 ~ -tr\u208Ap\u208Ai(t) + tr\u208Ar1\u208Ap\u208Ai(t)" %in% fmt)
  # the unconnected pins of the top level are inside connectors with zero flow
  expect_true("0 ~ tr\u208Ap\u208Ai(t)" %in% fmt)
  expect_true("0 ~ tr\u208An\u208Ai(t)" %in% fmt)
})

test_that("stream connections substitute instream", {
  t <- t_nounits
  Port <- connector(function(name) {
    m <- Variable("m", t, connect = Flow); h <- Variable("h", t, connect = Stream)
    P <- Variable("P", t)
    System(list(), t, list(P, m, h), list(), name = name)
  })
  expect_identical(connector_type(Port(name = "x")), "StreamConnector")
  Src <- function(name, h0) {
    port <- Port(name = "port"); hin <- Variable("hin", t)
    System(list(port$h ~ h0, hin ~ instream(port$h), port$m ~ -1, port$P ~ 1), t, list(hin), list(),
           systems = list(port), name = name)
  }
  a <- Src("a", 2); b <- Src("b", 3)
  sys <- System(list(connect(a$port, b$port)), t, systems = list(a, b), name = "s")
  fmt <- vapply(equations(expand_connections(sys)), format, "")
  expect_true("a\u208Ahin(t) ~ b\u208Aport\u208Ah(t)" %in% fmt)
  expect_true("b\u208Ahin(t) ~ a\u208Aport\u208Ah(t)" %in% fmt)
  expect_true("a\u208Aport\u208AP(t) ~ b\u208Aport\u208AP(t)" %in% fmt)
  expect_true("0 ~ a\u208Aport\u208Am(t) + b\u208Aport\u208Am(t)" %in% fmt)
  # three inner ports: the mixing rule of the Modelica specification
  c <- Src("c", 4)
  sys3 <- System(list(connect(a$port, b$port, c$port)), t, systems = list(a, b, c), name = "s3")
  eqs <- equations(expand_connections(sys3))
  rhs <- Filter(function(e) format(e$lhs) == "a\u208Ahin(t)", eqs)[[1L]]$rhs
  n <- toggle_namespacing(sys3, FALSE)
  f <- build_function(rhs, list(n$b$port$m, n$c$port$m, n$b$port$h, n$c$port$h))
  expect_equal(f(c(-1, -3, 10, 20)), (1 * 10 + 3 * 20) / 4)
  # an outside connector joined to two inside ones gets the mixing equation;
  # as in the Julia source, positivemax receives `tol` as its second argument,
  # which makes the mixture an average for small tolerances
  J <- function(name) {
    o <- Port(name = "o"); s1 <- Src("s1", 2); s2 <- Src("s2", 3)
    System(list(connect(o, s1$port, s2$port)), t, systems = list(o, s1, s2), name = name)
  }
  top <- System(list(), t, systems = list(J("j")), name = "top")
  eqs <- equations(expand_connections(top))
  oe <- Filter(function(e) format(e$lhs) == "j\u208Ao\u208Ah(t)", eqs)[[1L]]
  n <- toggle_namespacing(top, FALSE)$j
  g <- build_function(oe$rhs, list(n$s1$port$m, n$s2$port$m, n$s1$port$h, n$s2$port$h))
  expect_equal(g(c(-1, -3, 2, 3)), 2.5, tolerance = 1e-6)
})

test_that("domain connections bind connector parameters to the domain source", {
  t <- t_nounits
  Fluid <- connector(function(name, rho = 1000) {
    dm <- Variable("dm", t, connect = Flow); r <- Parameter("rho", default = rho)
    System(list(), t, list(dm), list(r), name = name)
  })
  Port <- connector(function(name) {
    p <- Variable("p", t); dm <- Variable("dm", t, connect = Flow)
    r <- Parameter("rho", default = 1)
    System(list(), t, list(p, dm), list(r), name = name)
  })
  Tank <- function(name) {
    port <- Port(name = "port")
    System(list(port$p ~ 1, port$dm ~ port$rho / 1000), t, systems = list(port), name = name)
  }
  fluid <- Fluid(name = "fluid", rho = 876)
  expect_identical(connector_type(fluid), "DomainConnector")
  tank <- Tank("tank")
  sys <- System(list(connect(fluid, tank$port)), t, systems = list(fluid, tank), name = "s")
  ex <- expand_connections(sys)
  b <- smap_get(bindings(ex), toggle_namespacing(sys, FALSE)$tank$port$rho)
  expect_identical(format(b), "fluid\u208Arho")
  fmt <- vapply(equations(ex), format, "")
  expect_true("0 ~ fluid\u208Adm(t) + tank\u208Aport\u208Adm(t)" %in% fmt)
  prob <- ODEProblem(mtkcompile(sys), list(), c(0, 1))
  expect_equal(prob$ps[toggle_namespacing(sys, FALSE)$tank$port$rho], 876)
  expect_equal(prob[toggle_namespacing(sys, FALSE)$fluid$dm], -0.876)
  f2 <- Fluid(name = "f2")
  bad <- System(list(connect(fluid, f2, tank$port)), t, systems = list(fluid, f2, tank), name = "s2")
  expect_error(expand_connections(bad), "multiple source domain connectors")
})

test_that("unbalanced connectors warn once per type", {
  t <- t_nounits
  Bad <- connector(function(name) {
    System(list(), t, list(Variable("i", t, connect = Flow)), list(), name = name)
  })
  Bad2 <- connector(function(name) {
    System(list(), t, list(Variable("a", t), Variable("b", t), Variable("i", t, connect = Flow)),
           list(), name = name)
  }, type = "Bad2")
  sys <- System(list(), t, systems = list(Bad2(name = "x"), Bad2(name = "y")), name = "s")
  expect_warning(expand_connections(sys), "Bad2 contains 1 flow variables, yet 2 regular")
  expect_identical(connector_type(Bad(name = "z")), "DomainConnector")
})

test_that("connections are namespaced with their system", {
  lib <- .rc_library()
  a <- lib$Pin(name = "a"); b <- lib$Pin(name = "b")
  inner <- System(list(connect(a, b)), t_nounits, systems = list(a, b), name = "inner")
  outer <- System(list(), t_nounits, systems = list(inner), name = "outer")
  conn <- Filter(function(e) inherits(e, "Connection"), equations(outer))[[1L]]
  expect_identical(format(conn), "connect(inner.a, inner.b)")
  # inside inner, a and b are outside connectors; from outer they are
  # unconnected inside connectors, whose flows are zero
  fmt <- vapply(equations(expand_connections(outer)), format, "")
  expect_setequal(fmt, c("inner\u208Aa\u208Av(t) ~ inner\u208Ab\u208Av(t)", "0 ~ -inner\u208Aa\u208Ai(t) - inner\u208Ab\u208Ai(t)",
                         "0 ~ inner\u208Aa\u208Ai(t)", "0 ~ inner\u208Ab\u208Ai(t)"))
})
