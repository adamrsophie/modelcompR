# Continuation algorithms. Reference values printed by
# tools/nonlinearsolve-homotopy-probe.jl (NonlinearSolve 4.32.0,
# NonlinearSolveBase 2.55.1); LAM is the sequence of distinct lambda values at
# which the residual was evaluated.

.hx_problems <- list(
  atan = list(function(u, l) c(l * atan(u[1] - 3) + (1 - l) * u[1]),
              function(u, l) matrix(l / (1 + (u[1] - 3)^2) + (1 - l), 1, 1), 12),
  scurve = list(function(u, l) c(u[1]^3 - 3 * u[1] - (8 * l - 4)), function(u, l) matrix(3 * u[1]^2 - 3, 1, 1), -2.2),
  circ = list(function(u, l) c(l * (u[1]^2 + u[2]^2 - 4) + (1 - l) * (u[1] - 2), l * (u[1] * u[2] - 1) + (1 - l) * (u[2] - 0.5)),
              function(u, l) rbind(c(l * 2 * u[1] + (1 - l), l * 2 * u[2]), c(l * u[2], l * u[1] + (1 - l))), c(2, 0.5)))

.hx_alg <- function(tag) switch(tag,
  SimpleSweep = SimpleHomotopySweep(),
  SimpleSweepConst = SimpleHomotopySweep(predictor = "constant"),
  SimpleSweepFixed = SimpleHomotopySweep(adaptive = FALSE, nsteps = 4),
  KantNR = KantorovichHomotopy(inner = NewtonRaphson()),
  KantNRSecant = KantorovichHomotopy(inner = NewtonRaphson(), predictor = "secant"),
  Kant = KantorovichHomotopy(),
  ArcLength = ArcLengthContinuation(),
  ArcLengthTangent = ArcLengthContinuation(predictor = "tangent"),
  ArcLengthNR = ArcLengthContinuation(inner = NewtonRaphson()),
  Poly = HomotopyPolyAlgorithm(),
  PolyNoHandoff = HomotopyPolyAlgorithm(warm_handoff = FALSE),
  FSHP = FastShortcutHomotopyPolyalg())

.hx_run <- function(pn, tag) {
  pr <- .hx_problems[[pn]]
  lams <- numeric()
  at <- function(l) list(F = function(x) {
    if (!length(lams) || lams[length(lams)] != l) lams <<- c(lams, l)
    pr[[1]](x, l)
  }, J = function(x) pr[[2]](x, l))
  r <- .nl_sweep(at, pr[[3]], .hx_alg(tag), list(), 0, 1)
  r$lams <- lams
  r
}

.hx_ref <- list(
  list("atan", "SimpleSweep", "Success", c(3.0),
       c(0.0, 0.1, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0)),
  list("atan", "SimpleSweepConst", "Success", c(3.0),
       c(0.0, 0.1, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0)),
  list("atan", "SimpleSweepFixed", "Success", c(3.0),
       c(0.0, 0.25, 0.5, 0.75, 1.0)),
  list("atan", "KantNR", "Success", c(3.0),
       c(0.0, 0.1, 0.40943899315001697, 0.718877986300034, 1.0)),
  list("atan", "KantNRSecant", "Success", c(3.0),
       c(0.0, 0.1, 0.2714546714127428, 0.5654217149052589, 0.8593887583977751, 1.0)),
  list("atan", "Kant", "Success", c(3.0),
       c(0.0, 0.1, 0.40943899315001697, 0.718877986300034, 1.0)),
  list("atan", "ArcLength", "Success", c(3.0),
       NULL),
  list("atan", "ArcLengthTangent", "Success", c(3.0000000000000004),
       NULL),
  list("atan", "ArcLengthNR", "Success", c(3.0),
       NULL),
  list("atan", "Poly", "Success", c(3.0),
       c(0.0, 0.1, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0)),
  list("atan", "PolyNoHandoff", "Success", c(3.0),
       c(0.0, 0.1, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0)),
  list("atan", "FSHP", "Success", c(3.0),
       c(0.0, 0.1, 0.2, 0.35000000000000003, 0.5, 0.65, 0.8, 1.0)),
  list("scurve", "SimpleSweep", "Success", c(2.195823345445647),
       c(0.0, 0.1, 0.2, 0.4, 0.6000000000000001, 0.8, 0.7000000000000001, 0.8, 0.9, 1.0)),
  list("scurve", "SimpleSweepConst", "Success", c(2.195823345445647),
       c(0.0, 0.1, 0.2, 0.4, 0.6000000000000001, 0.8, 0.7000000000000001, 0.8, 0.9, 1.0)),
  list("scurve", "SimpleSweepFixed", "MaxIters", c(-1.7320508075688774),
       c(0.0, 0.25, 0.5, 0.75)),
  list("scurve", "KantNR", "Success", c(2.195823345445647),
       c(0.0, 0.1, 0.40943899315001697, 1.0)),
  list("scurve", "KantNRSecant", "ConvergenceFailure", c(-1.0001705616187537),
       c(0.0, 0.1, 0.2714546714127428, 0.5654217149052589, 0.8593887583977751, 0.6242151236037622, 0.7250191693001319, 0.8258232149965016, 0.7451799784394059, 0.7653407875786798, 0.7492121402672607, 0.7532443020951155, 0.7500185726328317, 0.7493734267403749, 0.7496499599328861, 0.7501240890094533, 0.7497447857481995, 0.7499073690382596, 0.7500699523283197, 0.7499398856962717, 0.7499956370254206, 0.7500513883545695, 0.7500067872912504, 0.7499978670785865, 0.7500016511513756, 0.7499986238931443, 0.7499999214870576, 0.750001219080971, 0.7500001810058403, 0.7499999733908141, 0.7500000518875187, 0.7499999890901551, 0.7500000140696756)),
  list("scurve", "Kant", "Success", c(2.195823345445647),
       c(0.0, 0.1, 0.40943899315001697, 1.0)),
  list("scurve", "ArcLength", "Success", c(2.195823345445647),
       NULL),
  list("scurve", "ArcLengthTangent", "Success", c(2.195823345445649),
       NULL),
  list("scurve", "ArcLengthNR", "Success", c(2.195823345445647),
       NULL),
  list("scurve", "Poly", "Success", c(2.195823345445647),
       c(0.0, 0.1, 0.2, 0.4, 0.6000000000000001, 0.8, 0.7000000000000001, 0.8, 0.8500000000000001, 0.9000000000000001, 0.9750000000000001, 1.0)),
  list("scurve", "PolyNoHandoff", "Success", c(2.195823345445647),
       c(0.0, 0.1, 0.2, 0.4, 0.6000000000000001, 0.8, 0.7000000000000001, 0.8, 0.8500000000000001, 0.9000000000000001, 0.9750000000000001, 1.0)),
  list("scurve", "FSHP", "Success", c(2.1958233454456413),
       NULL),
  list("circ", "SimpleSweep", "Success", c(1.9318516525781366, 0.5176380902050416),
       c(0.0, 0.1, 0.2, 0.30000000000000004, 0.4, 0.6000000000000001, 0.8, 1.0)),
  list("circ", "SimpleSweepConst", "Success", c(1.9318516525781366, 0.5176380902050416),
       c(0.0, 0.1, 0.2, 0.30000000000000004, 0.4, 0.6000000000000001, 0.8, 1.0)),
  list("circ", "SimpleSweepFixed", "Success", c(1.9318516525781366, 0.5176380902050416),
       c(0.0, 0.25, 0.5, 0.75, 1.0)),
  list("circ", "KantNR", "Success", c(1.9318516525781366, 0.5176380902050415),
       c(0.0, 0.1, 0.40943899315001697, 0.718877986300034, 1.0)),
  list("circ", "KantNRSecant", "Success", c(1.9318516525781364, 0.5176380902050416),
       c(0.0, 0.1, 0.2714546714127428, 0.44290934282548566, 0.6143640142382285, 0.7858186856509713, 1.0)),
  list("circ", "Kant", "Success", c(1.9318516525781366, 0.5176380902050415),
       c(0.0, 0.1, 0.40943899315001697, 0.718877986300034, 1.0)),
  list("circ", "ArcLength", "Success", c(1.9318516525781588, 0.5176380902050195),
       NULL),
  list("circ", "ArcLengthTangent", "Success", c(1.9318516525781597, 0.5176380902050188),
       NULL),
  list("circ", "ArcLengthNR", "Success", c(1.9318516525781588, 0.5176380902050195),
       NULL),
  list("circ", "Poly", "Success", c(1.9318516525781366, 0.5176380902050415),
       c(0.0, 0.1, 0.2, 0.30000000000000004, 0.4, 0.6000000000000001, 0.8, 1.0)),
  list("circ", "PolyNoHandoff", "Success", c(1.9318516525781366, 0.5176380902050415),
       c(0.0, 0.1, 0.2, 0.30000000000000004, 0.4, 0.6000000000000001, 0.8, 1.0)),
  list("circ", "FSHP", "Success", c(1.9318516525781415, 0.5176380902050367),
       c(0.0, 0.1, 0.15000000000000002, 0.17500000000000002, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.8999999999999999, 1.0)))

test_that("continuation algorithms follow the NonlinearSolve path", {
  for (rf in .hx_ref) {
    r <- .hx_run(rf[[1]], rf[[2]])
    lab <- paste(rf[[1]], rf[[2]])
    expect_identical(r$retcode, rf[[3]], label = lab)
    expect_lte(max(abs(r$u - rf[[4]])), 1e-12 * max(1, abs(rf[[4]])), label = lab)
    # the arclength path differs in its dH/dlambda evaluations (central
    # differences in R, ForwardDiff in Julia), so only the sweeps compare paths
    if (!is.null(rf[[5]])) expect_equal(r$lams, rf[[5]], tolerance = 1e-12, label = lab)
  }
})

test_that("continuation algorithms solve HomotopyProblem and validate settings", {
  y <- Variable("y")
  sys <- complete(System(list(0 ~ homotopy(atan(y - 3), y)), name = "hx"))
  prob <- HomotopyProblem(sys, list(y ~ 12))
  for (alg in list(SimpleHomotopySweep(), KantorovichHomotopy(inner = NewtonRaphson()), ArcLengthContinuation(),
                   HomotopyPolyAlgorithm(), FastShortcutHomotopyPolyalg())) {
    sol <- solve(prob, alg)
    expect_identical(sol$retcode, "Success")
    expect_equal(unname(sol$u), 3, tolerance = 1e-10)
  }
  expect_error(SimpleHomotopySweep(adaptive = FALSE), "explicit `nsteps`")
  expect_error(KantorovichHomotopy(qmin = 2), "qmin")
  expect_error(KantorovichHomotopy(theta_min = 0.9), "theta_min")
  expect_error(ArcLengthContinuation(theta = 1), "theta")
  expect_error(HomotopyPolyAlgorithm(list()), "at least one")
  expect_s3_class(HomotopyPolyAlgorithm()$stages[[2]], "ArcLengthContinuation")
})
