# Work precision tools of DiffEqDevTools.jl (src/benchmark.jl,
# src/autoplot.jl). Reference values from tools/DiffEqDevTools-probe.jl.
# Only errors, names, tags and counts are compared; times are not.

dv_wp_prob <- function(analytic = TRUE) {
  f <- function(u, p, t) c(-u[1] * u[1], u[1] * u[2])
  an <- function(u0, p, t) c(1 / (1 + t), 2 * (1 + t))
  if (analytic) ODEProblem(ODEFunction(f, analytic = an), c(1, 2), c(0, 1)) else ODEProblem(f, c(1, 2), c(0, 1))
}

dv_tols <- c(1e-3, 1e-4, 1e-5)

test_that("WorkPrecisionSet and WorkPrecision of ODEs match DiffEqDevTools", {
  prob <- dv_wp_prob()
  setups <- list(list(alg = Euler(), dts = c(1 / 8, 1 / 16, 1 / 32)),
    list(alg = RK4(), dts = c(1 / 4, 1 / 8, 1 / 16), adaptive = FALSE))
  wps <- WorkPrecisionSet(prob, dv_tols, dv_tols, setups, numruns = 2, error_estimates = c("final", "l2", "L2"))
  expect_equal(wps$names, c("Euler", "RK4"))
  expect_equal(wps[[1]]$errors$final, c(0.04773824454478556, 0.02419923396292173, 0.012194281310070193),
    tolerance = 1e-12)
  expect_equal(wps[[1]]$errors$L2, c(0.028610962082216473, 0.01455207778753441, 0.007349468014173333),
    tolerance = 1e-12)
  expect_equal(wps[[2]]$errors$l2, c(0.00024379759703634882, 1.3615067246749641e-5, 8.023536739292566e-7),
    tolerance = 1e-9)
  expect_equal(wps[[1]]$tags, c("order_1", "fixed_step", "explicit"))
  # RK4 is adaptive (Shampine's residual estimate), as Julia tags it
  expect_equal(wps[[2]]$tags, c("order_4", "adaptive", "explicit"))
  expect_equal(available_errors(wps), c("final", "l2", "L2"))
  wp <- WorkPrecision(prob, RK4(), dv_tols, dv_tols, c(1 / 4, 1 / 8, 1 / 16), numruns = 2, error_estimate = "l_inf",
    name = "R", adaptive = FALSE)
  expect_equal(wp$errors$l_inf, c(0.0005015858798054218, 2.8078247553686708e-5, 1.657186966319557e-6),
    tolerance = 1e-9)
  expect_equal(names(wp$errors), c("final", "l2", "l_inf"))
  wp <- WorkPrecision(prob, Tsit5(), dv_tols, dv_tols, numruns = 2)
  expect_equal(wp$errors$final, c(1.9972291112613227e-6, 4.0023480279982593e-7, 7.08805942273294e-8),
    tolerance = 1e-8)
  expect_equal(vapply(wp$stats, function(s) as.numeric(s$naccept), 0), c(5, 6, 8))
  ref <- solve(prob, RK4(), dt = 1 / 256, adaptive = FALSE)
  wp <- WorkPrecision(dv_wp_prob(FALSE), Heun(), dv_tols, dv_tols, c(1 / 8, 1 / 16, 1 / 32), numruns = 2,
    appxsol = ref, error_estimate = "L2")
  expect_equal(wp$errors$L2, c(0.0004924589093903681, 6.096040639626874e-5, 6.60937193000892e-6), tolerance = 1e-9)
  sh <- Shootout(prob, list(list(alg = Euler()), list(alg = RK4())), dt = 1 / 16, numruns = 2, adaptive = FALSE)
  expect_equal(sh$names, c("Euler", "RK4"))
  expect_equal(as.numeric(sh$errors) / c(0.02419923396292173, 8.515960957522495e-7), c(1, 1), tolerance = 1e-9)
})

test_that("WorkPrecisionSet with several problems matches DiffEqDevTools", {
  prob <- dv_wp_prob()
  probs <- list(prob, remake(prob, tspan = c(0, 0.5)))
  setups <- list(list(alg = RK4(), dts = c(1 / 4, 1 / 8, 1 / 16), adaptive = FALSE),
    list(alg = RK4(), dts = c(1 / 4, 1 / 8, 1 / 16), adaptive = FALSE, prob_choice = 2))
  refs <- lapply(probs, function(p) solve(p, Tsit5(), abstol = 1e-12, reltol = 1e-12))
  w <- WorkPrecisionSet(probs, dv_tols, dv_tols, setups, numruns = 2, appxsol = refs)
  expect_equal(w[[1]]$errors$final, c(0.00025580159361310306, 1.4397699622636928e-5, 8.515960839838854e-7),
    tolerance = 1e-8)
  expect_equal(w[[2]]$errors$final, c(0.00017578294121217697, 1.003243561631395e-5, 5.960429766416198e-7),
    tolerance = 1e-8)
  # Julia fails without appxsol (getindex(::Nothing, ::Int64), probe line
  # "pc noappx"); R uses the analytic solution of each problem.
  w <- WorkPrecisionSet(probs, dv_tols, dv_tols, setups, numruns = 2)
  expect_equal(w[[1]]$errors$final[3], 8.515960957522495e-7, tolerance = 1e-9)
})

test_that("WorkPrecisionSet of DDEs and nonlinear problems matches DiffEqDevTools", {
  fdd <- function(u, h, p, t) -h(p, t - 1)
  andd <- function(u0, h, p, t) if (t <= 1) 1 - t else -(2 * (t - 1) - (t * t - 1) / 2)
  ddprob <- DDEProblem(DDEFunction(fdd, analytic = andd), 1, function(p, t) 1, c(0, 2), constant_lags = 1)
  dw <- WorkPrecisionSet(ddprob, c(1e-3, 1e-4), c(1e-3, 1e-4),
    list(list(alg = MethodOfSteps(Euler()), dts = c(1 / 8, 1 / 16))), numruns = 2, error_estimate = "l2")
  expect_equal(dw$names, "Euler")
  expect_equal(dw[[1]]$errors$l2, c(0.02706329386826363, 0.013150397055158887), tolerance = 1e-12)
  expect_equal(dw[[1]]$tags, c("order_1", "fixed_step"))
  fnl <- function(u, p) c(u[1] * u[1] * u[1] + u[1] - p[1], u[2] * u[2] - p[2])
  nprob <- NonlinearProblem(fnl, c(1, 1), c(2, 3))
  nwp <- WorkPrecisionSet(nprob, c(1e-6, 1e-10), c(1e-6, 1e-10), list(list(alg = NewtonRaphson())), numruns = 2,
    error_estimate = "l2")
  # The residual after the last Newton step differs in the last digits from
  # Julia (u[2] differs by 2e-13), so the comparison is loose.
  expect_equal(nwp[[1]]$errors$l2[1] / 8.47267411785424e-9, 1, tolerance = 1e-3)
  expect_lt(nwp[[1]]$errors$l2[2], 1e-14)
  expect_error(WorkPrecision(nprob, NewtonRaphson(), 1e-6, 1e-6, numruns = 2, error_estimate = "final"),
    "Unsupported norm used: final.", fixed = TRUE)
})

dv_wp_sde16 <- function() {
  tg <- seq(0, 1, by = 1 / 16)
  W <- lapply(seq_along(tg), function(k) {
    t <- tg[k]
    c(0.3 * t - 0.7 * t * t + 0.01 * (k %% 3), -0.2 * t + 0.5 * t * t * t - 0.02 * (k %% 5))
  })
  W[[1]] <- c(0, 0)
  an <- function(u0, p, t, W) c(u0[1] + t * t + 0.5 * W[1], u0[2] + t + 0.25 * W[2])
  SDEProblem(SDEFunction(function(u, p, t) c(2 * t, 1), function(u, p, t) c(0.5, 0.25), analytic = an),
    c(1, 2), c(0, 1), noise = NoiseGrid(tg, W))
}

test_that("WorkPrecisionSet of SDEs and ensembles matches DiffEqDevTools", {
  sprob <- dv_wp_sde16()
  setups <- list(list(alg = EM(), dts = rep(1 / 16, 3)), list(alg = EulerHeun(), dts = rep(1 / 16, 3)))
  ref <- list(final = c(0.03125000000000078, 7.771561172376096e-16),
    l2 = c(0.025911131174651544, 5.596643484763003e-16), weak_final = c(0.06250000000000022, 1.3506446028928519e-15))
  for (ee in names(ref)) {
    rw <- suppressMessages(WorkPrecisionSet(sprob, dv_tols, dv_tols, setups, numruns = 2, numruns_error = 2,
      error_estimate = ee))
    expect_equal(rw$names, c("EM", "EulerHeun"))
    expect_equal(rw[[1]]$errors[[ee]], rep(ref[[ee]][1], 3), tolerance = 1e-12)
    expect_equal(rw[[2]]$errors[[ee]], rep(ref[[ee]][2], 3), tolerance = 1e-6)
    # Julia passes the number of setups (2) as N; R the number of tolerances.
    expect_equal(rw[[1]]$N, 3)
    expect_equal(rw[[1]]$tags, c("order_1_2", "fixed_step", "sde"))
  }
  expect_equal(names(rw[[1]]$errors), "weak_final")
  eset <- function(p, ee, ...) suppressMessages(WorkPrecisionSet(p, dv_tols, dv_tols, setups, numruns = 2,
    trajectories = 2, ensemblealg = EnsembleSerial(), error_estimate = ee, ...))
  eprob <- EnsembleProblem(sprob, output_func = function(sol, i) list(sol$u[[length(sol$u)]][1], FALSE))
  ew <- eset(eprob, "weak_final", expected_value = 2)
  expect_equal(ew[[1]]$errors$weak_final, rep(0.25250000000000017, 3), tolerance = 1e-12)
  expect_equal(ew[[2]]$errors$weak_final, rep(0.19000000000000017, 3), tolerance = 1e-12)
  eprob2 <- EnsembleProblem(sprob, output_func = function(sol, i) list(sol$u[[length(sol$u)]], FALSE))
  ew <- eset(eprob2, "weak_l2", expected_value = c(2, 3))
  expect_equal(ew[[1]]$errors$weak_l2, rep(0.09375000000000078, 3), tolerance = 1e-12)
  expect_equal(ew[[2]]$errors$weak_l2, rep(0.06250000000000078, 3), tolerance = 1e-12)
  ew <- eset(eprob, "weak_final", appxsol_setup = list(alg = EM(), dt = 1 / 16))
  expect_equal(ew[[2]]$errors$weak_final, rep(0.0625, 3), tolerance = 1e-12)
  expect_error(eset(eprob, "final", expected_value = 2),
    "use RODEProblem instead of EnsembleProblem for strong errors.", fixed = TRUE)
  expect_error(eset(eprob, "weak_l_inf", expected_value = 2), "is not implemented yet.", fixed = TRUE)
})

dv_raw_set <- function() {
  prob <- dv_wp_prob()
  mk <- function(name, errs, times, tags) WorkPrecision(prob, dv_tols, dv_tols, list(final = errs), times, NULL,
    NULL, name, "final", 3, tags)
  WorkPrecisionSet(list(mk("A", c(1e-2, 1e-4, 1e-6), c(1e-3, 2e-3, 5e-3), c("rk", "explicit")),
    mk("B", c(1e-3, 1e-5, 1e-7), c(2e-3, 4e-3, 9e-3), c("rk", "explicit")),
    mk("C", c(1e-2, NaN, 1e-6), c(3e-3, 5e-3, 7e-3), c("rosenbrock", "implicit")),
    mk("D", c(2e-2, 3e-4, 1e-5), c(1e-3, 2e-3, 3e-3), c("rosenbrock", "implicit", "reference"))),
    4, dv_tols, dv_tols, prob, rep(list(list(alg = Euler())), 4), c("A", "B", "C", "D"), "final", NULL)
}

test_that("tag selection and autoplot match DiffEqDevTools", {
  set <- dv_raw_set()
  expect_equal(vapply(set$wps, wp_area, 0),
    c(-10.698970004336019, -9.540607512240769, -9.355561410532161, -9.05385249113878), tolerance = 1e-14)
  expect_equal(best_by_tag(set, "rk")$names, "A")
  expect_equal(best_by_tag(set, "rk", n = 2)$names, c("A", "B"))
  expect_equal(best_by_tag(set, "rosenbrock")$names, "D")
  expect_equal(best_of_families(set, c("rk", "rosenbrock"))$names, c("A", "D"))
  expect_error(best_of_families(set, "sdirk"), "No entries found for any of the family tags")
  expect_error(best_by_tag(set, "rk", metric = "time"), "Unknown metric: time. Supported metrics:")
  expect_equal(filter_by_tags(set, "rk", "explicit")$names, c("A", "B"))
  expect_equal(filter_by_tags(set, "implicit")$names, c("C", "D"))
  expect_equal(exclude_by_tags(set, "reference")$names, c("A", "B", "C"))
  expect_equal(exclude_by_tags(set, "rk", "reference")$names, "C")
  expect_equal(sort(unique_tags(set)), c("explicit", "implicit", "reference", "rk", "rosenbrock"))
  expect_equal(get_tags(set)[[4]], c("rosenbrock", "implicit", "reference"))
  ap <- autoplot(set, reference_tags = "reference", best_n = 1)
  expect_setequal(names(ap), c("all", "best_of_families", "family_rk", "family_rosenbrock"))
  expect_equal(ap$best_of_families$names, c("A", "D"))
  expect_equal(ap$family_rk$names, c("A", "B", "D"))
  expect_equal(ap$family_rosenbrock$names, c("C", "D"))
  expect_equal(ap$all$names, c("A", "B", "C", "D"))
  ap <- autoplot(set)
  expect_equal(ap$best_of_families$names, c("A", "B", "C", "D"))
  expect_equal(ap$family_rk$names, c("A", "B"))
  m <- merge_wp_sets(set, best_by_tag(set, "rk"))
  expect_equal(m$names, c("A", "B", "C", "D", "A"))
  expect_equal(m$N, 5)
  expect_equal(suppressWarnings(filter_by_tags(set, "nothing"))$N, 0)
})

test_that("auto_tags, tag_kind and with_autodiff_variants match DiffEqDevTools", {
  expect_equal(auto_tags(Euler()), c("order_1", "fixed_step", "explicit"))
  expect_equal(auto_tags(Tsit5()), c("order_5", "adaptive", "explicit"))
  expect_equal(auto_tags(Rodas5P()), c("order_5", "adaptive", "implicit", "rosenbrock"))
  expect_equal(auto_tags(DFBDF()), c("order_1", "adaptive", "dae", "implicit"))
  expect_equal(auto_tags(EM()), c("order_1_2", "fixed_step", "sde"))
  expect_equal(auto_tags(SRIW1()), c("order_3_2", "adaptive", "sde"))
  expect_equal(auto_tags(RKMil()), c("order_1_1", "adaptive", "sde"))
  expect_equal(auto_tags(SRA1()), c("order_2_1", "adaptive", "sde"))
  expect_equal(auto_tags(SKenCarp()), c("order_2_1", "adaptive", "sde"))
  expect_equal(auto_tags(AutoTsit5(Rosenbrock23())), c("order_5", "adaptive", "implicit", "composite"))
  expect_equal(auto_tags(MethodOfSteps(Tsit5())), c("order_5", "adaptive"))
  expect_equal(auto_tags(NewtonRaphson()), character())
  expect_equal(auto_tags(DImplicitEuler()), c("order_1", "adaptive", "dae", "implicit"))
  expect_equal(auto_tags(Rodas4()), c("order_4", "adaptive", "implicit", "rosenbrock"))
  expect_equal(auto_tags(DRI1()), c("order_1_1", "adaptive", "sde"))
  expect_equal(auto_tags(RS1()), c("order_1_1", "fixed_step", "sde"))
  expect_equal(auto_tags(SOSRA()), c("order_2_1", "adaptive", "sde"))
  expect_equal(auto_tags(SplitEM()), c("order_1_2", "fixed_step", "sde", "split"))
  expect_equal(auto_tags(RandomEM()), c("order_1_2", "fixed_step", "rode"))
  expect_equal(auto_tags(RandomTaylor15()), c("order_3_2", "fixed_step", "rode"))
  # boundary value algorithms: no order, adaptive by the SciMLBase default
  expect_equal(auto_tags(MIRK4()), "adaptive")
  expect_equal(auto_tags(Shooting(Tsit5())), "adaptive")
  expect_equal(auto_tags(MultipleShooting(4, Tsit5())), "adaptive")
  expect_equal(tag_kind(c("rk", "sdirk", "adaptive", "order_5", "reference", "sundials", "autodiff_x", "dense",
    "stiff", "foo")), c("family", "family", "trait", "trait", "role", "provider", "variant", "variant", "domain",
    "unknown"))
  v <- with_autodiff_variants(list(list(alg = Rodas5P()), list(alg = Tsit5(), tags = "nonstiff")),
    ad_backends = list(AutoFiniteDiff(), AutoForwardDiff()))
  expect_length(v, 4)
  expect_null(v[[1]]$name)
  expect_equal(v[[2]]$name, "Rodas5P (FiniteDiff)")
  expect_equal(v[[3]]$name, "Rodas5P (ForwardDiff)")
  expect_equal(lapply(v, `[[`, "tags"), list("autodiff_default", "autodiff_finitediff", "autodiff_forwarddiff",
    c("nonstiff", "autodiff_default")))
})

test_that("autoplot answers the ggplot2 generic", {
  skip_if_not_installed("ggplot2")
  set <- dv_raw_set()
  expect_setequal(names(ggplot2::autoplot(set)), c("all", "best_of_families", "family_rk", "family_rosenbrock"))
  # other objects go to ggplot2 (no recursion through an R default method)
  expect_error(autoplot(data.frame(x = 1)), "not supported by autoplot")
})
