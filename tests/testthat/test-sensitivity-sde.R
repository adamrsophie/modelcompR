# Reference values from tools/scimlsensitivity-probe.jl (SciMLSensitivity
# 7.119.14): the forward solves and the adjoints run on one fixed noise path.

sde_tg <- seq(0, 1, by = 0.05)
sde_w1 <- function(t) 0.3 * sin(3.1 * t) + 0.25 * sin(7.3 * t) - 0.1 * sin(13 * t)
sde_w2 <- function(t) 0.4 * sin(5.2 * t) - 0.2 * sin(11.1 * t)
sde_w3 <- function(t) 0.35 * sin(2.3 * t) + 0.15 * sin(17 * t)
sde_f <- function(u, p, t) c(p[1] * u[1] - u[1] * u[2], -p[2] * u[2] + u[1] * u[2] + 0.1 * t)
sde_gdiag <- function(u, p, t) c(p[3] * u[1], p[4] * u[2]^2)
sde_p <- c(1.2, 0.9, 0.3, 0.25)
sde_u0 <- c(1, 0.8)
sde_dg <- function(u, p, t, i) u
sde_grid2 <- function() NoiseGrid(sde_tg, cbind(sde_w1(sde_tg), sde_w2(sde_tg)))

sde_check <- function(sol, sensealg, du0, dp, alg = EulerHeun(), tol = 1e-8, ...) {
  r <- adjoint_sensitivities(sol, alg, t = c(0.5, 1), dgdu_discrete = sde_dg, sensealg = sensealg, ...)
  expect_equal(unname(r$du0), du0, tolerance = tol)
  expect_equal(unname(r$dp), dp, tolerance = tol)
}

test_that("SDEAdjointProblem matches Julia for diagonal noise", {
  sol <- solve(SDEProblem(sde_f, sde_gdiag, sde_u0, c(0, 1), sde_p, noise = sde_grid2()), EulerHeun(), dt = 0.05)
  expect_equal(sol$u[[21]], c(1.3646737278407568, 1.1273419105106566), tolerance = 1e-14)
  sde_check(sol, BacksolveAdjoint(), c(3.9754009961255847, -0.8296080105425372),
            c(3.0873998756660734, -0.4910055594652466, 0.7119511251266634, -0.14061091021802846))
  sde_check(sol, InterpolatingAdjoint(), c(3.9738277598131146, -0.8293275698406211),
            c(3.087144217004288, -0.491384288040898, 0.7115648193667982, -0.14019538680328994))
  # an Ito forward solve. Julia's adjoints are biased here: its
  # BacksolveAdjoint differentiates the transformed drift f - c for lambda
  # and the gradient (du0 = (3.6237449643094606, -0.8929885564022412), dp =
  # (2.954889559093849, -0.37388803081650507, -1.1470841971095112,
  # -0.49575315139355697)) and its InterpolatingAdjoint leaves out the Ito
  # correction (du0 = (3.9685419462739335, -0.8958133522845659)). Both
  # adjoints here approach the pathwise derivative of the forward scheme
  prob <- SDEProblem(sde_f, sde_gdiag, sde_u0, c(0, 1), sde_p, noise = sde_grid2())
  sol <- solve(prob, EM(), dt = 0.05)
  fw <- solve_gradient(function(s) sum(as.matrix(s)^2) / 2, prob, EM(), dt = 0.05, saveat = c(0.5, 1),
                       sensealg = ForwardDiffSensitivity())
  ref <- c(fw$u0, fw$p)
  for (sa in list(BacksolveAdjoint(), InterpolatingAdjoint())) {
    r <- adjoint_sensitivities(sol, EM(), t = c(0.5, 1), dgdu_discrete = sde_dg, sensealg = sa)
    expect_lt(max(abs(c(r$du0, r$dp) - ref)) / max(abs(ref)), 0.05)
  }
})

test_that("SDEAdjointProblem matches Julia for scalar, general and mixing noise", {
  sol <- solve(SDEProblem(sde_f, sde_gdiag, sde_u0, c(0, 1), sde_p, noise = NoiseGrid(sde_tg, sde_w1(sde_tg))),
               EulerHeun(), dt = 0.05)
  sde_check(sol, BacksolveAdjoint(), c(3.9304329241588385, -0.7821902391304539),
            c(2.9992018678915366, -0.5821637737813367, 0.7033535558004187, -0.2887577295603916))
  sde_check(sol, InterpolatingAdjoint(), c(3.9292575015371294, -0.7818523098752052),
            c(2.9989857147175463, -0.582515658229029, 0.7026601985731322, -0.28811285395328684))
  ggen <- function(u, p, t) rbind(c(p[3] * u[1], 0.2 * u[2], p[4] * u[1] * u[2]),
                                  c(0.1, p[3] * u[2]^2, p[4] * cos(u[1])))
  noise3 <- NoiseGrid(sde_tg, cbind(sde_w1(sde_tg), sde_w2(sde_tg), sde_w3(sde_tg)))
  sol <- solve(SDEProblem(sde_f, ggen, sde_u0, c(0, 1), sde_p, noise = noise3, noise_rate_prototype = matrix(0, 2, 3)),
               EulerHeun(), dt = 0.05)
  sde_check(sol, BacksolveAdjoint(), c(4.236479918122186, -0.20902919608540593),
            c(3.1416374093137747, -0.663498381188943, 0.6433542933234166, 1.2716388202065443))
  sde_check(sol, InterpolatingAdjoint(), c(4.233996552862944, -0.20973176816197123),
            c(3.1424942089905947, -0.6635199620566226, 0.6416833456601871, 1.272104733718106))
  gmix <- function(u, p, t) c(p[3] * u[1] * u[2], p[4] * u[1])
  sol <- solve(SDEProblem(sde_f, gmix, sde_u0, c(0, 1), sde_p, noise = sde_grid2()), EulerHeun(), dt = 0.05)
  sde_check(sol, BacksolveAdjoint(noisemixing = TRUE), c(3.6091032050839607, -0.2661054072677112),
            c(2.793609536947338, -0.5213581012052876, 0.2607686249119221, -0.2868971104739691))
  sde_check(sol, InterpolatingAdjoint(noisemixing = TRUE), c(3.6084612640672358, -0.2666027154957468),
            c(2.793914476751979, -0.5212719043054403, 0.2624683543310595, -0.2877672307564562))
  # without noisemixing Julia uses only the diagonal of the noise Jacobian
  # (du0 = (3.5643, -0.4712), dp = (2.7952, -0.4321, 0.2064, -0.0347), a
  # wrong gradient for this noise); here the off diagonal entries at u0
  # switch the full Jacobian on
  sde_check(sol, BacksolveAdjoint(), c(3.6091032050839607, -0.2661054072677112),
            c(2.793609536947338, -0.5213581012052876, 0.2607686249119221, -0.2868971104739691))
})

test_that("costs at every saved time, solve_gradient and continuous costs match Julia", {
  prob <- SDEProblem(sde_f, sde_gdiag, sde_u0, c(0, 1), sde_p, noise = sde_grid2())
  sol <- solve(prob, EulerHeun(), dt = 0.05)
  # t0 included; at a checkpoint time the backsolved state is reset before the cost
  back_du0 <- c(72.3664595968196, 1.1503691155715159)
  back_dp <- c(38.939225741444744, -6.622173159715885, 15.566450444224984, -3.933164700851843)
  interp_du0 <- c(72.34504748307933, 1.1525402045924895)
  interp_dp <- c(38.93881076983239, -6.626444310742121, 15.562355358880955, -3.9289801611911863)
  r <- adjoint_sensitivities(sol, EulerHeun(), t = sol$t, dgdu_discrete = function(u, p, t, i) 2 * u,
                             sensealg = BacksolveAdjoint())
  expect_equal(unname(r$du0), back_du0, tolerance = 1e-8)
  expect_equal(unname(r$dp), back_dp, tolerance = 1e-8)
  # the adjoint of the solve: the cotangent of the saved states (Julia: the
  # same discrete cost with the cotangent fixed at the saved states)
  loss <- function(s) sum(as.matrix(s)^2)
  g <- solve_gradient(loss, prob, EulerHeun(), dt = 0.05, sensealg = BacksolveAdjoint())
  expect_equal(unname(g$u0), back_du0, tolerance = 1e-8)
  expect_equal(unname(g$p), back_dp, tolerance = 1e-8)
  g <- solve_gradient(loss, prob, EulerHeun(), dt = 0.05, sensealg = InterpolatingAdjoint())
  expect_equal(unname(g$u0), interp_du0, tolerance = 1e-8)
  expect_equal(unname(g$p), interp_dp, tolerance = 1e-8)
  # with saveat the adjoint runs on the grid of a solve that saves every step:
  # sum(u(0.5)^2 + u(1)^2) / 2 is the cost of diag_strat_interp
  g <- solve_gradient(function(s) sum(as.matrix(s)^2) / 2, prob, EulerHeun(), dt = 0.05, saveat = c(0.5, 1),
                      sensealg = InterpolatingAdjoint())
  expect_equal(unname(g$u0), c(3.9738277598131146, -0.8293275698406211), tolerance = 1e-8)
  expect_equal(unname(g$p), c(3.087144217004288, -0.491384288040898, 0.7115648193667982, -0.14019538680328994),
               tolerance = 1e-8)
  # Julia's adjoint diffusion also subtracts the gradient of a continuous cost
  # (Back du0 = (1.7488039427632456, -0.09842622852296067), Interp du0 =
  # (1.7486563482957083, -0.09726997281047389)); R keeps it in the drift only
  r <- adjoint_sensitivities(sol, EulerHeun(), dgdu_continuous = function(u, p, t) u, sensealg = BacksolveAdjoint())
  expect_equal(unname(r$du0), c(1.7308737239350551, 0.036850919204772001), tolerance = 1e-8)
  expect_equal(unname(r$dp), c(0.91653214046085629, -0.16112787650984703, 0.37871008233323139, -0.09527275315558649),
               tolerance = 1e-8)
  r <- adjoint_sensitivities(sol, EulerHeun(), dgdu_continuous = function(u, p, t) u,
                             dgdp_continuous = function(u, p, t) c(u[1] * u[2], 0, 0, 0), sensealg = InterpolatingAdjoint())
  expect_equal(unname(r$du0), c(1.7306865427947586, 0.037601753726542689), tolerance = 1e-8)
  expect_equal(unname(r$dp), c(2.0883122749910155, -0.16154104458073865, 0.37862039402266023, -0.095234128915909377),
               tolerance = 1e-8)
  # with zero diffusion the adjoint has no noise term, so the gradient of a
  # continuous cost does not depend on the noise path; Julia's does, as its
  # accumulate_cost! also runs in the adjoint diffusion (contzero_a_* and
  # contzero_b_* of the probe)
  gzero <- function(u, p, t) c(0, 0)
  grid2x <- NoiseGrid(sde_tg, cbind(2 * sde_w1(sde_tg), -sde_w2(sde_tg)))
  for (sa in list(BacksolveAdjoint(), InterpolatingAdjoint())) {
    ra <- adjoint_sensitivities(solve(SDEProblem(sde_f, gzero, sde_u0, c(0, 1), sde_p, noise = sde_grid2()),
                                      EulerHeun(), dt = 0.05), EulerHeun(), dgdu_continuous = function(u, p, t) u,
                                sensealg = sa)
    rb <- adjoint_sensitivities(solve(SDEProblem(sde_f, gzero, sde_u0, c(0, 1), sde_p, noise = grid2x),
                                      EulerHeun(), dt = 0.05), EulerHeun(), dgdu_continuous = function(u, p, t) u,
                                sensealg = sa)
    expect_equal(rb$du0, ra$du0, tolerance = 1e-12)
    expect_equal(rb$dp, ra$dp, tolerance = 1e-12)
  }
  # InterpolatingAdjoint solved with EM: Julia samples the cost at t = 0.5 with
  # the stale state u(0.55) (du0 = (3.9647706834542507, -0.8995838768260501));
  # ReverseLossCallback reads the buffer of the last drift evaluation, which
  # upstream already replaces by the forward solution for DAEs. The cost with
  # the cotangent fixed at the forward states is the same cost (Julia's
  # em_interpfix of the probe, which differs from its em_interp)
  sol <- solve(prob, EM(), dt = 0.05)
  r <- adjoint_sensitivities(sol, EM(), t = c(0.5, 1), dgdu_discrete = sde_dg, sensealg = InterpolatingAdjoint())
  fw <- solve_gradient(function(s) sum(as.matrix(s)^2) / 2, prob, EM(), dt = 0.05, saveat = c(0.5, 1),
                       sensealg = ForwardDiffSensitivity())
  expect_lt(max(abs(c(r$du0, r$dp) - c(fw$u0, fw$p))) / max(abs(c(fw$u0, fw$p))), 0.05)
  rfix <- adjoint_sensitivities(sol, EM(), t = c(0.5, 1), dgdu_discrete = function(u, p, t, i) sol(t),
                                sensealg = InterpolatingAdjoint())
  expect_equal(rfix$du0, r$du0, tolerance = 1e-12)
  expect_equal(rfix$dp, r$dp, tolerance = 1e-12)
})

test_that("RODEAdjointProblem matches Julia", {
  frode <- function(u, p, t, W) c(p[1] * u[1] * sin(W[1]) - p[2] * u[2], p[3] * u[2] * W[2] + p[4] * u[1])
  rprob <- RODEProblem(frode, sde_u0, c(0, 1), sde_p, noise = sde_grid2())
  sol <- solve(rprob, RandomEM(), dt = 0.05)
  # solve_gradient: the cotangent is fixed at the forward states (rode_backfix
  # of the probe), not the backsolved ones of the cost function u
  g <- solve_gradient(function(s) sum(as.matrix(s)^2) / 2, rprob, RandomEM(), dt = 0.05, saveat = c(0.5, 1),
                      sensealg = BacksolveAdjoint())
  ufix <- list(sol(0.5), sol(1))
  r <- adjoint_sensitivities(sol, RandomEM(), t = c(0.5, 1), dgdu_discrete = function(u, p, t, i) ufix[[i]],
                             sensealg = BacksolveAdjoint())
  expect_equal(unname(c(g$u0, g$p)), unname(c(r$du0, r$dp)), tolerance = 1e-9)
  expect_equal(sol$u[[21]], c(0.3561407427805689, 0.9873657791924934), tolerance = 1e-14)
  sde_check(sol, BacksolveAdjoint(), c(1.7229698376181326, 1.1110203516733446),
            c(0.20387964283431442, -0.8534151734082994, 0.09463000292955186, 0.862264400467529), alg = RandomEM())
  # Julia's InterpolatingAdjoint samples the cost at t = 0.5 with the state of
  # its last drift evaluation, u(0.55) for RandomEM (du0 = (1.6752235059118823,
  # 1.1432754857128713)); R samples u(0.5), the same scheme otherwise
  sde_check(sol, InterpolatingAdjoint(), c(1.71484463694468, 1.111062280276795),
            c(0.20291043396132113, -0.8502265178080268, 0.09432078427923676, 0.86151471863776), alg = RandomEM())
  # the same cost with the cotangent fixed at the forward states (Julia's
  # rode_interpfix of the probe differs from its rode_interp)
  r <- adjoint_sensitivities(sol, RandomEM(), t = c(0.5, 1), dgdu_discrete = function(u, p, t, i) sol(t),
                             sensealg = InterpolatingAdjoint())
  expect_equal(unname(r$du0), c(1.71484463694468, 1.111062280276795), tolerance = 1e-8)
  expect_equal(unname(r$dp), c(0.20291043396132113, -0.8502265178080268, 0.09432078427923676, 0.86151471863776),
               tolerance = 1e-8)
  ap <- RODEAdjointProblem(sol, InterpolatingAdjoint(), t = c(0.5, 1), dgdu_discrete = sde_dg)
  expect_s3_class(ap, "RODEProblem")
  expect_identical(ap$tspan, c(0, 1))
})

test_that("SDEAdjointProblem builds the reversed time problem and rejects the unsupported cases", {
  sol <- solve(SDEProblem(sde_f, sde_gdiag, sde_u0, c(0, 1), sde_p, noise = sde_grid2()), EulerHeun(), dt = 0.05)
  ap <- SDEAdjointProblem(sol, BacksolveAdjoint(), t = c(0.5, 1), dgdu_discrete = sde_dg)
  expect_s3_class(ap, "SDEProblem")
  expect_s3_class(ap$alg, "EulerHeun")
  expect_equal(ap$dt, 0.05)
  expect_equal(ap$time_map(0.25), 0.75)
  expect_true(inherits(ap$f$f, SensitivityFunction))
  # z = (lambda, grad, y), the cost at T already applied
  expect_equal(ap$u0, c(sol$u[[21]], numeric(4), sol$u[[21]]))
  expect_error(SDEAdjointProblem(sol, QuadratureAdjoint(), t = 1, dgdu_discrete = sde_dg),
               "BacksolveAdjoint() and InterpolatingAdjoint()", fixed = TRUE)
  expect_error(SDEAdjointProblem(sol, InterpolatingAdjoint(), t = 1),
               "Either `dgdu_discrete`, `dgdu_continuous`, or `g` must be specified.", fixed = TRUE)
  expect_error(SDEAdjointProblem(sol, InterpolatingAdjoint(), t = 1, dgdu_continuous = function(u, p, t) u),
               "no specified `dgdu_discrete` or `dgdp_discrete`", fixed = TRUE)
  sol0 <- solve(SDEProblem(function(u, p, t) -u, function(u, p, t) 0.1 * u, 1, c(0, 1), NULL,
                           noise = NoiseGrid(sde_tg, sde_w1(sde_tg))), EulerHeun(), dt = 0.05)
  expect_error(SDEAdjointProblem(sol0, InterpolatingAdjoint(), t = 1, dgdu_discrete = function(u, p, t, i) 1),
               "Your model does not have parameters", fixed = TRUE)
})

test_that("pathwise SDE sensitivities differentiate discrete events", {
  f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1] * u[2], -p[3] * u[2] + p[4] * u[1] * u[2])
  cb <- CallbackSet(PresetTimeCallback(0.4, function(integ) integ$u[1] <- integ$u[1] * integ$p[2] + 0.1 * integ$u[2]),
                    PresetTimeCallback(0.7, function(integ) integ$p[3] <- integ$p[3] + 0.5 * integ$u[1]))
  prob <- SDEProblem(f, function(u, p, t) 0.1 * u, c(1, 1), c(0, 1), c(1.5, 1, 3, 1), callback = cb)
  ts <- seq(0, 1, by = 0.1)
  loss <- function(sol) sum((as.matrix(sol) - 1)^2)
  F <- function(th) loss(solve(remake(prob, u0 = th[1:2], p = th[3:6]), EM(), dt = 0.01, saveat = ts, seed = 3))
  h <- 1e-6
  th <- c(1, 1, 1.5, 1, 3, 1)
  ref <- vapply(seq_along(th), function(j) {
    a <- th; b <- th
    a[j] <- a[j] + h; b[j] <- b[j] - h
    (F(a) - F(b)) / (2 * h)
  }, 0)
  gr <- solve_gradient(loss, prob, EM(), dt = 0.01, saveat = ts, seed = 3, sensealg = ForwardDiffSensitivity())
  expect_equal(unname(c(gr$u0, gr$p)), ref, tolerance = 1e-7)
  pc <- remake(prob, callback = ContinuousCallback(function(u, t, integ) u[1] - 2, function(integ) NULL))
  expect_error(solve_gradient(loss, pc, EM(), dt = 0.01, saveat = ts, seed = 3, sensealg = ForwardDiffSensitivity()),
               "not differentiable")
  # the adjoints pull back through the same events and converge to the
  # pathwise derivative as dt falls (each cost applied once, also where the
  # accumulated step times fall within rounding of the save times)
  for (dt in c(1e-3, 2.5e-4)) {
    fw <- solve_gradient(loss, prob, EM(), dt = dt, saveat = ts, seed = 3, sensealg = ForwardDiffSensitivity())
    ga <- solve_gradient(loss, prob, EM(), dt = dt, saveat = ts, seed = 3, sensealg = InterpolatingAdjoint())
    expect_lt(max(abs(c(ga$u0, ga$p) - c(fw$u0, fw$p))) / max(abs(c(fw$u0, fw$p))), 2 * dt)
  }
  expect_error(solve_gradient(loss, pc, EM(), dt = 0.01, saveat = ts, seed = 3, sensealg = InterpolatingAdjoint()),
               "not differentiable")
})

test_that("a nonsingular mass matrix multiplies the drift in the SDE sensitivities", {
  # the solvers read M du = f dt + M g dW: the sensitivities use M^-1 f
  M <- matrix(c(2, 0.5, 0, 1), 2)
  f <- function(u, p, t) c(-p[1] * u[1] + u[2], -p[2] * u[2])
  g <- function(u, p, t) c(0.3 * u[1], p[3] + 0.1 * u[2])
  mk <- function(u0, p) SDEProblem(SDEFunction(f, g, mass_matrix = M), u0, c(0, 1), p, seed = 7)
  ts <- seq(0, 1, by = 0.25)
  loss <- function(sol) sum(as.matrix(sol)^2)
  dt <- 1 / 256
  F <- function(th) loss(solve(mk(th[1:2], th[3:5]), ImplicitEM(), dt = dt, adaptive = FALSE, saveat = ts))
  th <- c(1, 0.5, 0.7, 0.4, 0.2)
  ref <- vapply(seq_along(th), function(j) {
    a <- th; b <- th
    a[j] <- a[j] + 1e-6; b[j] <- b[j] - 1e-6
    (F(a) - F(b)) / 2e-6
  }, 0)
  prob <- mk(c(1, 0.5), c(0.7, 0.4, 0.2))
  fw <- solve_gradient(loss, prob, ImplicitEM(), dt = dt, adaptive = FALSE, saveat = ts,
                       sensealg = ForwardDiffSensitivity())
  expect_equal(unname(c(fw$u0, fw$p)), ref, tolerance = 1e-7)
  ga <- solve_gradient(loss, prob, ImplicitEM(), dt = dt, adaptive = FALSE, saveat = ts,
                       sensealg = InterpolatingAdjoint())
  expect_lt(max(abs(c(ga$u0, ga$p) - ref)) / max(abs(ref)), 0.01)
})

test_that("the Ito adjoints converge to the pathwise derivative for noise nonlinear in the state", {
  # the backward Ito solve evaluates the adjoint diffusion with the forward
  # state at the end of each step; the correction sum_jl q_jl D_{g_l} B_j
  # removes the bias that does not vanish with dt (about 10% here before)
  f <- function(u, p, t) c(-p[1] * u[1] + u[2], -p[2] * u[2])
  gs <- list(diagonal = function(u, p, t) c(0.3 * u[1] * u[2], p[3] * sin(u[1])),
             general = function(u, p, t) matrix(c(0.3 * u[1], 0.1 * u[2]^2, p[3] * u[2], 0.2 * sin(u[1]),
                                                  0.1, 0.1 * u[1]), 2),
             scalar = function(u, p, t) c(0.3 * u[1]^2, p[3] * u[2]))
  args <- list(diagonal = list(), general = list(noise_rate_prototype = matrix(0, 2, 3)),
               scalar = list(noise = "scalar"))
  ts <- seq(0, 1, by = 0.25)
  loss <- function(sol) sum(as.matrix(sol)^2)
  for (nm in names(gs)) {
    prob <- do.call(SDEProblem, c(list(f, gs[[nm]], c(1, 0.5), c(0, 1), c(0.7, 0.4, 0.2)), args[[nm]]))
    for (sa in list(InterpolatingAdjoint(noisemixing = TRUE), BacksolveAdjoint(noisemixing = TRUE))) {
      fw <- solve_gradient(loss, prob, EM(), dt = 1 / 512, saveat = ts, seed = 3, sensealg = ForwardDiffSensitivity())
      ga <- solve_gradient(loss, prob, EM(), dt = 1 / 512, saveat = ts, seed = 3, sensealg = sa)
      expect_lt(max(abs(c(ga$u0, ga$p) - c(fw$u0, fw$p))) / max(abs(c(fw$u0, fw$p))), 2e-3,
                label = paste(nm, sa$name))
    }
  }
})

test_that("SDE adjoints keep the forward steps, terminal events and the final parameters", {
  # saveat keeps fewer times than the forward steps
  prob6 <- SDEProblem(function(u, p, t) p[1] * u, function(u, p, t) 0 * u, 1, c(0, 1), 1)
  sol6 <- solve(prob6, EM(), dt = 0.001, saveat = c(0, 1), dense = TRUE)
  r <- adjoint_sensitivities(sol6, EM(), t = 1, dgdu_discrete = function(u, p, t, i) 1, sensealg = InterpolatingAdjoint())
  expect_equal(r$du0, 1.001^1000, tolerance = 1e-6)
  # an event at the final time
  prob7 <- SDEProblem(function(u, p, t) 0 * u, function(u, p, t) 0 * u, 1, c(0, 1), 2,
                      callback = PresetTimeCallback(1, function(i) i$u[1] <- i$p[1] * i$u[1]))
  for (sa in list(InterpolatingAdjoint(), BacksolveAdjoint(), ForwardDiffSensitivity())) {
    g <- solve_gradient(function(sol) tail(as.matrix(sol)[, 1], 1), prob7, EM(), dt = 0.01, sensealg = sa)
    expect_equal(unname(c(g$u0, g$p)), c(2, 1), label = sa$name)
  }
  # a loss of the parameters after an event changed them
  prob8 <- SDEProblem(function(u, p, t) 0 * u, function(u, p, t) 0 * u, 1, c(0, 1), 1,
                      callback = PresetTimeCallback(0.5, function(i) i$p[1] <- 2 * i$p[1]))
  loss8 <- function(sol) .sol_env(sol)$p_final[1]^2 + tail(as.matrix(sol)[, 1], 1)
  for (sa in list(ForwardDiffSensitivity(), InterpolatingAdjoint())) {
    g <- solve_gradient(loss8, prob8, EM(), dt = 0.01, sensealg = sa)
    expect_equal(unname(c(g$u0, g$p, g$value)), c(1, 8, 5), label = sa$name)
  }
  # mixed interpretations need the quadratic variation of a Wiener path
  tg <- seq(0, 1, by = 0.001)
  prob5 <- SDEProblem(function(u, p, t) 0 * u, function(u, p, t) p[1] * u, 1, c(0, 1), 1, noise = NoiseGrid(tg, 0 * tg))
  sol5 <- solve(prob5, EulerHeun(), dt = 0.001)
  expect_error(adjoint_sensitivities(sol5, EM(), t = 1, dgdu_discrete = function(u, p, t, i) 1,
                                     sensealg = InterpolatingAdjoint()), "same interpretation")
  r <- adjoint_sensitivities(sol5, EulerHeun(), t = 1, dgdu_discrete = function(u, p, t, i) 1, sensealg = InterpolatingAdjoint())
  expect_equal(c(r$du0, r$dp), c(1, 0), tolerance = 1e-10)
  # a well conditioned small mass matrix
  pm <- SDEProblem(SDEFunction(function(u, p, t) -p[1] * u, function(u, p, t) 0.1 * u, mass_matrix = 1e-8 * diag(2)),
                   c(1, 1), c(0, 0.1), 1, seed = 3)
  expect_no_error(solve_gradient(function(sol) sum(as.matrix(sol)^2), pm, ImplicitEM(), dt = 0.01,
                                 adaptive = FALSE, sensealg = InterpolatingAdjoint()))
})
