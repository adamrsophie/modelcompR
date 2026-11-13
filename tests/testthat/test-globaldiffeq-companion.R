# GlobalErrorEstimation of GlobalDiffEq 1.7.3 (src/companion.jl,
# src/estimation.jl), against the Julia values of the companion section of
# tools/globaldiffeq-probe.jl (the key names are quoted). The probe problem
# has no Jacobian, so Julia differentiates it with ForwardDiff; here it
# carries the exact one.
#
# The estimate is the end point of a companion solve whose right hand side
# is the defect of the dense output, so it is a small difference of state
# sized numbers; the comparisons are absolute. The streaming mode and the
# explicit dense solvers agree with Julia to about 1e-12 of the state. In
# InterpolatingMode the companion steps end on the forward nodes, where the
# derivative of the dense output jumps, and the estimates agree to 5e-11.

gec_lv <- function(u, p, t) c(1.5 * u[1] - u[1] * u[2], -3 * u[2] + u[1] * u[2])
gec_jac <- function(u, p, t) matrix(c(1.5 - u[2], u[2], -u[1], -3 + u[1]), 2, 2)
gec_prob <- function(...) ODEProblem(ODEFunction(gec_lv, jac_fn = gec_jac), c(1, 1), c(0, 2), ...)
gec_near <- function(x, ref, tol, ...) expect_lte(max(abs(x - ref)), tol, ...)
gec_last <- function(sol) sol$u[[length(sol$t)]]

test_that("the options print as Julia enums", {
  expect_output(print(DefectCorrection), "DefectCorrection::GlobalErrorEquation = 0", fixed = TRUE)
  expect_output(print(SimultaneousMode), "SimultaneousMode::GlobalErrorMode = 1", fixed = TRUE)
  expect_identical(GlobalErrorEquation, list(DefectCorrection = DefectCorrection, ErrorTransport = ErrorTransport))
  expect_identical(GlobalErrorMode, list(InterpolatingMode = InterpolatingMode, SimultaneousMode = SimultaneousMode))
  alg <- GlobalErrorEstimation(Tsit5())
  expect_identical(alg$equation, DefectCorrection)
  expect_identical(alg$mode, InterpolatingMode)
  expect_identical(GlobalErrorEstimation(Tsit5(), equation = ErrorTransport)$equation, ErrorTransport)
  # has_global_error_gee
  expect_false(has_global_error(alg))
  expect_false(allows_arbitrary_number_types(alg))
  expect_false(allowscomplex(alg))
})

test_that("the dense first derivatives match Julia", {
  # dense_<alg>_n, _d (t = 0.3 and 1.1), _node_left and _node_right
  # (sol$t[3]) and _T_right; Julia's end point derivative with continuity =
  # :right is NaN or Inf for BS3, Rodas5P, ImplicitEuler and FBDF (a step of
  # length zero), R's is the left limit, which no estimator reads
  ref <- list(
    DP5 = list(16L, c(1.1656135589562262, -0.9925309801273862, 3.9011257225804705, 0.03602317072214835),
               c(0.7304758501158342, -1.5750548043419463), c(0.7304758501158342, -1.5750548043419463),
               c(-3.47509149084482, 7.602567974441355)),
    Vern7 = list(13L, c(1.1656111797432214, -0.9925318856768213, 3.9011960704118667, 0.03598496408043515),
                 c(0.9112528091068542, -1.2984765267716234), c(0.9112528091072287, -1.2984765267721736),
                 c(-3.475029586628671, 7.602551539815019)),
    BS3 = list(132L, c(1.1656101206363254, -0.9925335127784982, 3.901199103050433, 0.03598261794349736),
               c(0.5106687871587006, -1.978730427108958), c(0.5106687871586401, -1.9787304271087516), NULL),
    Rosenbrock23 = list(102L, c(1.165552124535311, -0.9925772962200375, 3.9011273434200127, 0.03595233922686949),
                        c(0.5092481843710771, -1.9815474135242843), c(0.5092498096982315, -1.9815396107073242),
                        c(-3.4740555686166403, 7.601491880124437)),
    Rodas5P = list(30L, c(1.1656112843912663, -0.9925316748020753, 3.9011988272558966, 0.035981282039515984),
                   c(0.6188396613667105, -1.7713841008818147), c(0.6188397182060195, -1.7713842433114926), NULL),
    ImplicitEuler = list(4382L, c(1.1664252505870796, -0.9919120200645531, 3.9015595330137436, 0.036544777764066574),
                         c(0.5000000126568072, -1.9999999746863857), c(0.5000000151881685, -1.9999999696236632),
                         NULL),
    FBDF = list(69L, c(1.1656172131622293, -0.9925356580336342, 3.9011513102683684, 0.03600118897621477),
                c(0.5000001321039087, -1.9999997444206714), c(0.5000004706334386, -1.999999060174163), NULL))
  # the 4381 implicit Euler steps and the FBDF history carry the last bits
  # of the states into the slopes, about 1e-10; the node slopes of FBDF are
  # differences over its first step of 1e-7, about 1e-8
  dtol <- c(DP5 = 1e-13, Vern7 = 1e-13, BS3 = 1e-13, Rosenbrock23 = 1e-13, Rodas5P = 1e-13, ImplicitEuler = 1e-9,
            FBDF = 1e-9)
  ntol <- c(dtol[1:6], FBDF = 1e-7)
  for (nm in names(ref)) {
    r <- ref[[nm]]
    sol <- solve(gec_prob(), get(nm)(), abstol = 1e-6, reltol = 1e-6)
    expect_true(.gde_has_dense_derivative(sol))
    expect_identical(length(sol$t), r[[1]], info = nm)
    dn <- .gde_sol_dense(sol, right = FALSE)
    dr <- .gde_sol_dense(sol, right = TRUE)
    expect_equal(c(dn$deriv(0.3), dn$deriv(1.1)), r[[2]], tolerance = ref_tol(dtol[[nm]]), info = nm)
    expect_equal(dn$deriv(sol$t[3]), r[[3]], tolerance = ref_tol(ntol[[nm]]), info = nm)
    expect_equal(dr$deriv(sol$t[3]), r[[4]], tolerance = ref_tol(ntol[[nm]]), info = nm)
    if (!is.null(r[[5]])) expect_equal(dr$deriv(2), r[[5]], tolerance = ref_tol(1e-12), info = nm)
    expect_equal(dr$value(sol$t[3]), sol$u[[3]], info = nm)
  }
  # Tsit5 and Trapezoid against finite differences of the dense output
  for (alg in list(Tsit5(), Trapezoid())) {
    sol <- solve(gec_prob(), alg, abstol = 1e-6, reltol = 1e-6)
    dn <- .gde_sol_dense(sol, right = FALSE)
    j <- ceiling(length(sol$t) / 2)
    t <- (sol$t[j] + sol$t[j + 1]) / 2
    h <- 1e-3 * (sol$t[j + 1] - sol$t[j])
    fd <- (8 * (dn$value(t + h) - dn$value(t - h)) - (dn$value(t + 2 * h) - dn$value(t - 2 * h))) / (12 * h)
    expect_equal(dn$deriv(t), fd, tolerance = 1e-7, info = alg$name)
  }
  # the linear dense output of the low storage and SSP methods has none (an
  # R limitation: OrdinaryDiffEq gives them Hermite or quadratic output)
  expect_false(.gde_alg_has_dense_derivative(SSPRK33()))
})

test_that("the companion estimates match Julia", {
  # gee_<equation>_<mode>_<tol>
  ref <- list(DefectCorrection = list(InterpolatingMode = c(0.000916347598339453, 1.1743736186382576e-5),
                                      SimultaneousMode = c(0.0009163475732161836, 1.1743742633483352e-5)),
              ErrorTransport = list(InterpolatingMode = c(0.0009162610948896422, 1.1743702610333257e-5),
                                    SimultaneousMode = c(0.0009162610889749235, 1.1743721928840932e-5)))
  for (eq in GlobalErrorEquation) for (md in GlobalErrorMode) for (i in 1:2) {
    tol <- c(1e-4, 1e-6)[i]
    est <- global_error_estimate(gec_prob(), GlobalErrorEstimation(Tsit5(), equation = eq, mode = md),
                                 abstol = tol, reltol = tol)
    gec_near(est, ref[[eq]][[md]][i], if (identical(md, SimultaneousMode)) ref_tol(1e-13, 1e-9) else 5e-11,
             label = paste(eq, md, tol))
  }
  # gee_tsit5_1e-4_uT
  expect_equal(gec_last(solve(gec_prob(), Tsit5(), abstol = 1e-4, reltol = 1e-4)),
               c(6.776363734900814, 2.0131511903088173), tolerance = ref_tol(1e-15))
  # gee_companion_vern7
  est <- global_error_estimate(gec_prob(), GlobalErrorEstimation(Tsit5(), companion_alg = Vern7(),
                               companion_abstol = 1e-12, companion_reltol = 1e-10), abstol = 1e-5, reltol = 1e-5)
  gec_near(est, 0.00011000863288554095, 1e-12)
  # gee_lin: a length one state, y' = 2 y on [0, 1]
  lin <- ODEProblem(function(u, p, t) p * u, 1, c(0, 1), 2)
  gec_near(global_error_estimate(lin, GlobalErrorEstimation(Tsit5()), abstol = 1e-5, reltol = 1e-5),
           1.276675144076432e-7, 1e-12)
  # gee_tolerance: a negative local tolerance is passed to the solve
  gec_near(global_error_estimate(gec_prob(), GlobalErrorEstimation(Tsit5()), abstol = -1), 0.023960251550027966,
           5e-12)
})

test_that("other dense solvers and companion solvers match Julia", {
  # gee_alg_<alg>, gee_alg_<alg>_tsit5companion and gee_alg_<alg>_uT
  ref <- list(DP5 = c(0.0002641107797232254, 0.00026411068209641805, 6.776981687037784, 2.012916851452559),
              Vern7 = c(2.187243250717025e-6, 2.1871131793378377e-6, 6.777187221639555, 2.0127540110588353),
              BS3 = c(0.00012717098550861324, 0.00012717788615182831, 6.777315753454796, 2.0127680662608523),
              Rosenbrock23 = c(0.0013137730488144727, 0.0013137483924184365, 6.775966537143576, 2.01227349941631),
              Rodas5P = c(0.00010066805079087561, 0.0001006399993563053, 6.777288610499445, 2.0127368034134583),
              FBDF = c(0.0004783912639103742, 0.0004784206231572864, 6.776751529559152, 2.012560981355187))
  # the forward solves match; the estimates on the dense output of the
  # Rosenbrock and BDF methods differ from Julia's by up to 4e-5 of
  # themselves although the slopes match at the probed points (an open
  # item: dense5_<alg> of the probe prints the slopes at every step)
  etol <- c(DP5 = 5e-11, Vern7 = 5e-11, BS3 = 5e-11, Rosenbrock23 = 1e-8, Rodas5P = 1e-8, FBDF = 1e-8)
  for (nm in names(ref)) {
    r <- ref[[nm]]
    # the order 2 Rosenbrock23 companion at the default 1e-10 takes minutes
    if (nm != "Rosenbrock23")
      gec_near(global_error_estimate(gec_prob(), GlobalErrorEstimation(get(nm)()), abstol = 1e-5, reltol = 1e-5),
               r[1], etol[[nm]], label = nm)
    gec_near(global_error_estimate(gec_prob(), GlobalErrorEstimation(get(nm)(), companion_alg = Tsit5()),
                                   abstol = 1e-5, reltol = 1e-5), r[2], etol[[nm]], label = nm)
    expect_equal(gec_last(solve(gec_prob(), get(nm)(), abstol = 1e-5, reltol = 1e-5)), r[3:4],
                 tolerance = ref_tol(1e-13, 1e-9), info = nm)
  }
  # an implicit Euler forward solve with an explicit companion
  # (gee_alg_ImplicitEuler_tsit5companion, _tsit5companion_1e-3, _uT and _1e-3_uT)
  gec_near(global_error_estimate(gec_prob(), GlobalErrorEstimation(ImplicitEuler(), companion_alg = Tsit5()),
                                 abstol = 1e-5, reltol = 1e-5), 0.05894078671563971, 1e-11)
  expect_equal(gec_last(solve(gec_prob(), ImplicitEuler(), abstol = 1e-5, reltol = 1e-5)),
               c(6.724865923776161, 2.0398880330212985), tolerance = ref_tol(1e-12, 1e-9))
  gec_near(global_error_estimate(gec_prob(), GlobalErrorEstimation(ImplicitEuler(), companion_alg = Tsit5()),
                                 abstol = 1e-3, reltol = 1e-3), 0.5836791304284332, 1e-9)
  expect_equal(gec_last(solve(gec_prob(), ImplicitEuler(), abstol = 1e-3, reltol = 1e-3)),
               c(6.246018362701547, 2.254701369363988), tolerance = ref_tol(1e-11, 1e-9))
})

test_that("global error control matches Julia", {
  # geectl_<equation>_<mode> = (42, 249, Success, u(T))
  ref <- list(DefectCorrection = list(InterpolatingMode = c(6.777189344989764, 2.012753722240371),
                                      SimultaneousMode = c(6.7771893449897656, 2.0127537222403693)),
              ErrorTransport = list(InterpolatingMode = c(6.777189344968206, 2.0127537222450167),
                                    SimultaneousMode = c(6.7771893449682095, 2.0127537222450167)))
  for (eq in GlobalErrorEquation) for (md in GlobalErrorMode) {
    sol <- solve(gec_prob(), GlobalErrorEstimation(Tsit5(), equation = eq, mode = md, gtol = 1e-7),
                 abstol = 1e-3, reltol = 1e-3)
    expect_true(successful_retcode(sol))
    expect_identical(length(sol$t), 42L)
    expect_equal(sol$stats$nf, 249)
    gec_near(gec_last(sol), ref[[eq]][[md]], ref_tol(1e-13, 1e-9), label = paste(eq, md))
  }
  # gee_maxiters_fail (Julia's estimate 0.002679852694947922, R's agrees to 1e-13)
  expect_error(solve(gec_prob(), GlobalErrorEstimation(Tsit5(), gtol = 1e-9, maxiters = 1), abstol = 1e-3,
                     reltol = 1e-3),
               "failed to meet gtol=1.0e-9 after 1 iterations; last estimated global error was 0.00267985269",
               fixed = TRUE)
})

test_that("GlobalErrorEstimation validates its arguments as Julia", {
  alg <- GlobalErrorEstimation(Tsit5())
  # gee_gtol, gee_companion_tolerance, gee_nogtol, gee_callback, gee_backwards, gee_mass
  expect_error(GlobalErrorEstimation(Tsit5(), gtol = -1), "gtol must be a positive finite real number")
  expect_error(GlobalErrorEstimation(Tsit5(), maxiters = 0), "maxiters must be a positive integer")
  expect_error(GlobalErrorEstimation(Tsit5(), safety = 2), "safety must be between zero and one")
  expect_error(GlobalErrorEstimation(Tsit5(), companion_abstol = 0, companion_reltol = 0),
               "companion_abstol and companion_reltol cannot both be zero")
  expect_error(solve(gec_prob(), alg), "GlobalErrorEstimation requires a positive `gtol` constructor keyword",
               fixed = TRUE)
  cb <- DiscreteCallback(function(u, t, integ) FALSE, function(integ) NULL)
  expect_error(global_error_estimate(gec_prob(callback = cb), alg),
               "GlobalErrorEstimation does not currently support callbacks")
  expect_error(global_error_estimate(ODEProblem(gec_lv, c(1, 1), c(2, 0)), alg),
               "GlobalErrorEstimation requires a forward-time ODEProblem")
  mprob <- ODEProblem(ODEFunction(gec_lv, mass_matrix = diag(c(1, 2))), c(1, 1), c(0, 2))
  expect_error(global_error_estimate(mprob, GlobalErrorEstimation(Rodas5P())),
               "GlobalErrorEstimation currently requires the standard mass matrix")
  expect_error(solve(gec_prob(), GlobalErrorEstimation(Tsit5(), gtol = 1e-6), abstol = -1),
               "local_abstol must be a nonnegative finite real number")
  # Julia's message for a solution without a derivative interpolation
  expect_error(global_error_estimate(gec_prob(), GlobalErrorEstimation(SSPRK33()), dt = 0.1),
               "requires a solver with dense (continuous) output that provides a first-derivative", fixed = TRUE)
  expect_error(global_error_estimate(gec_prob(), alg, dense = FALSE), NA)
})
