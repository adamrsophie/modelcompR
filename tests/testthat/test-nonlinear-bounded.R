# Reference values from tools/nonlinearsolve-bounded-probe.jl (NonlinearSolve
# 4.32.0, NonlinearSolveFirstOrder 2.10.0), analytic Jacobians. Columns: case,
# algorithm, maxiters, retcode, u, nf, njacs, nfactors, nsolve, nsteps. The
# powell_bd runs of the active set methods with maxiters = 1000 are left out:
# they stall at a minimum on a bound with a singular Jacobian, where the
# unpivoted QR of LAPACK and the one of R differ in the last bits and the
# trajectories part at the rounding level (their final points agree to 2e-9).
bounded_reference <- "
rosen_in BGN_ls 1000 Success 1.0;1.0 33 12 10 10 10
rosen_in BGN_ls 3 MaxIters -0.69189453125;0.05187187194824233 15 5 3 3 3
rosen_in BGN_tr 1000 Success 1.0;1.0 15 16 14 14 14
rosen_in BGN_tr 3 MaxIters -0.41644037913809306;0.11463546336574905 4 5 3 3 3
rosen_in BGN_trls 1000 Success 1.0;1.0 29 15 13 13 13
rosen_in BGN_trls 3 MaxIters -0.3351045374243133;0.12849717757968815 12 5 3 3 3
rosen_in BLM 1000 Success 0.9999999999999999;0.9999999999999999 33 14 15 15 12
rosen_in BLM 3 MaxIters -0.6321813582269329;0.009765625 13 5 6 6 3
rosen_in TRR 1000 Success 1.0;1.0 28 29 27 27 27
rosen_in TRR 3 MaxIters -0.9282397435917771;0.8574119336842506 4 5 3 3 3
rosen_in Dogbox 1000 Success 1.0;1.0 17 18 16 16 16
rosen_in Dogbox 3 MaxIters -0.34762948628790835;0.06507862784260471 4 5 3 3 3
rosen_in BTR 1000 Success 1.0;1.0 26 21 20 26 26
rosen_in BTR 3 MaxIters -0.9065030721167389;0.7423965191945423 3 3 2 3 3
rosen_bd BGN_ls 1000 Stalled 0.5;0.25 28 12 9 9 10
rosen_bd BGN_ls 3 MaxIters -0.69189453125;0.05187187194824233 15 5 3 3 3
rosen_bd BGN_tr 1000 Stalled 0.5;0.25 20 22 19 19 20
rosen_bd BGN_tr 3 MaxIters -0.41644037913809306;0.11463546336574905 4 5 3 3 3
rosen_bd BGN_trls 1000 Stalled 0.5;0.25 30 14 11 11 12
rosen_bd BGN_trls 3 MaxIters -0.3351045374243133;0.12849717757968815 12 5 3 3 3
rosen_bd BLM 1000 Stalled 0.5;0.25 21 12 18 18 10
rosen_bd BLM 3 MaxIters -0.5006886556828071;-0.15988048137274236 12 5 6 6 3
rosen_bd TRR 1000 Stalled 0.49999999999999994;0.25 22 23 21 21 21
rosen_bd TRR 3 MaxIters -0.9749984838000925;0.9501778071598939 4 5 3 3 3
rosen_bd Dogbox 1000 Stalled 0.5;0.25 11 13 10 10 11
rosen_bd Dogbox 3 MaxIters -0.34762948628790835;0.06507862784260471 4 5 3 3 3
rosen_bd BTR 1000 Stalled 0.5;0.24870073444322444 55 18 17 54 54
rosen_bd BTR 3 MaxIters -0.9065030721167389;0.7423965191945423 3 3 2 3 3
powell BGN_ls 1000 Success 2.838316418035712e-7;-2.838316418035712e-8;4.541306270347587e-8;4.5413062703475875e-8 24 25 23 23 23
powell BGN_ls 3 MaxIters 0.29761904761904795;-0.02976190476190481;0.047619047619047845;0.04761904761904785 4 5 3 3 3
powell BGN_tr 1000 Success 3.1529734223291857e-7;-3.152973422329187e-8;6.703375908662907e-8;6.703375908662905e-8 24 25 23 23 23
powell BGN_tr 3 MaxIters 0.3306132259014981;-0.03306132259014979;0.07028999093155684;0.07028999093155681 4 5 3 3 3
powell BGN_trls 1000 Success 3.1529734223291857e-7;-3.152973422329187e-8;6.703375908662907e-8;6.703375908662905e-8 24 25 23 23 23
powell BGN_trls 3 MaxIters 0.3306132259014981;-0.03306132259014979;0.07028999093155684;0.07028999093155681 4 5 3 3 3
powell BLM 1000 Success 2.846532153234431e-7;-2.8465321532344316e-8;4.617450669552875e-8;4.6174506695528766e-8 24 25 23 23 23
powell BLM 3 MaxIters 0.2984179958926056;-0.02984148279383022;0.04834741888551129;0.048360546111865604 4 5 3 3 3
powell TRR 1000 Success 3.2212574734029783e-7;-3.221257473402979e-8;1.4168690063680744e-7;1.4168690063680747e-7 25 26 24 24 24
powell TRR 3 MaxIters 0.7153128736639435;-0.06808764915054359;0.24490005503191487;0.3809795503100979 4 5 3 3 3
powell Dogbox 1000 Success 2.974414806415404e-7;-2.9744148064154046e-8;6.444461612574505e-8;6.444461612574503e-8 24 25 23 23 23
powell Dogbox 3 MaxIters 0.31188999800815687;-0.031188999800815683;0.06757507779617244;0.06757507779617243 4 5 3 3 3
powell BTR 1000 Success 2.838316417875744e-7;-2.838316417875744e-8;4.541306268601189e-8;4.541306268601189e-8 23 24 23 23 23
powell BTR 3 MaxIters 0.2976190476190476;-0.029761904761904753;0.047619047619047616;0.047619047619047616 3 4 3 3 3
powell_bd BGN_ls 3 MaxIters 0.5;-0.02976190476190481;0.047619047619047845;0.04761904761904785 4 5 3 3 3
powell_bd BGN_tr 3 MaxIters 0.5469005833336418;-0.06272979384217449;0.13379427918686262;0.25490585033246915 4 5 3 3 3
powell_bd BGN_trls 3 MaxIters 0.5469005833336418;-0.06272979384217449;0.13379427918686262;0.25490585033246915 4 5 3 3 3
powell_bd TRR 1000 Stalled 0.5000000000000001;-0.048452014344084075;0.18885708230235912;0.25077650853899597 27 28 26 26 26
powell_bd TRR 3 MaxIters 1.4421209957177434;-0.12948269375289698;0.3630983331452792;0.9514472415505227 4 5 3 3 3
powell_bd Dogbox 3 MaxIters 0.5;-0.05658345434918112;0.12310137347868795;0.22536573240358215 4 5 3 3 3
powell_bd BTR 1000 Stalled 0.5;-0.048486568857933274;0.18886201648319417;0.2507808961170696 132 101 100 131 131
powell_bd BTR 3 MaxIters 0.5;-0.05816516576964494;0.09252080772976616;0.20392959557125737 3 4 3 3 3
mm_ls BGN_ls 1000 Success 2.5;1.117491238129395 13 14 11 11 12
mm_ls BGN_ls 3 MaxIters 2.5;1.0874564356332974 5 5 3 3 3
mm_ls BGN_tr 1000 Success 2.5;1.1174912572634745 12 14 11 11 12
mm_ls BGN_tr 3 MaxIters 2.5;1.1229822931551527 4 5 3 3 3
mm_ls BGN_trls 1000 Success 2.5;1.1174912572634745 12 14 11 11 12
mm_ls BGN_trls 3 MaxIters 2.5;1.1229822931551527 4 5 3 3 3
mm_ls BLM 1000 Success 2.5;1.117491258298655 12 14 22 22 12
mm_ls BLM 3 MaxIters 2.5;1.1247606870552225 4 5 6 6 3
mm_ls TRR 1000 Success 2.4999999999999996;1.1174912475621503 13 15 12 12 13
mm_ls TRR 3 MaxIters 2.3939390204474074;1.031842389607402 4 5 3 3 3
mm_ls Dogbox 1000 Success 2.5;1.1174912512338095 11 13 10 10 11
mm_ls Dogbox 3 MaxIters 2.5;1.1183224508004646 4 5 3 3 3
mm_ls BTR 1000 Success 2.5;1.1174912631931133 11 13 11 11 12
mm_ls BTR 3 MaxIters 2.5;1.1329423691606066 3 4 3 3 3
mm_ls_free BGN_ls 1000 Success 3.2429548933214907;2.3184966040784016 9 11 8 8 9
mm_ls_free BGN_ls 3 MaxIters 3.2253913131346406;2.2157337156866115 4 5 3 3 3
mm_ls_free BGN_tr 1000 Success 3.2429549116947456;2.318496652872286 8 10 7 7 8
mm_ls_free BGN_tr 3 MaxIters 3.2445213393440073;2.3249852804247366 4 5 3 3 3
mm_ls_free BGN_trls 1000 Success 3.2429549116947456;2.318496652872286 8 10 7 7 8
mm_ls_free BGN_trls 3 MaxIters 3.2445213393440073;2.3249852804247366 4 5 3 3 3
mm_ls_free BLM 1000 Success 3.2429548875168632;2.3184965886634354 9 11 8 8 9
mm_ls_free BLM 3 MaxIters 3.235297645255168;2.2622518375452545 4 5 3 3 3
mm_ls_free TRR 1000 Success 3.2429548748276407;2.318496554964197 9 11 8 8 9
mm_ls_free TRR 3 MaxIters 3.2385373976590213;2.309370183556025 4 5 3 3 3
mm_ls_free Dogbox 1000 Success 3.2429548681153424;2.3184965371383006 9 11 8 8 9
mm_ls_free Dogbox 3 MaxIters 3.2450080764995626;2.3269978273521557 4 5 3 3 3
mm_ls_free BTR 1000 Success 3.2429548755563595;2.3184965568994604 8 10 8 8 9
mm_ls_free BTR 3 MaxIters 3.2418913176607385;2.3072679486316012 3 4 3 3 3
circ BGN_ls 1000 Success 1.9318516525781366;0.5176380902050415 6 7 5 5 5
circ BGN_ls 3 MaxIters 1.9319035496002646;0.5175908220067006 4 5 3 3 3
circ BGN_tr 1000 Success 1.9318516525781366;0.5176380902050415 6 7 5 5 5
circ BGN_tr 3 MaxIters 1.9319035496002646;0.5175908220067006 4 5 3 3 3
circ BGN_trls 1000 Success 1.9318516525781366;0.5176380902050415 6 7 5 5 5
circ BGN_trls 3 MaxIters 1.9319035496002646;0.5175908220067006 4 5 3 3 3
circ BLM 1000 Success 1.9318516525781366;0.5176380902050415 6 7 5 5 5
circ BLM 3 MaxIters 1.9319038395005337;0.5175905220986053 4 5 3 3 3
circ TRR 1000 Success 1.931851652578137;0.517638090205041 7 8 6 6 6
circ TRR 3 MaxIters 1.9416317306441053;0.5064899315322111 4 5 3 3 3
circ Dogbox 1000 Success 1.9318516525781366;0.5176380902050415 6 7 5 5 5
circ Dogbox 3 MaxIters 1.9319035496002646;0.5175908220067006 4 5 3 3 3
circ BTR 1000 Success 1.9318516525781366;0.5176380902050415 5 6 5 5 5
circ BTR 3 MaxIters 1.9319035496002646;0.5175908220067005 3 4 3 3 3
"

bounded_cases <- local({
  rosen <- function(u, p) c(10 * (u[2] - u[1] * u[1]), 1 - u[1])
  rosenj <- function(u, p) matrix(c(-20 * u[1], -1, 10, 0), 2)
  s5 <- sqrt(5); s10 <- sqrt(10)
  powell <- function(u, p) c(u[1] + 10 * u[2], s5 * (u[3] - u[4]), (u[2] - 2 * u[3]) * (u[2] - 2 * u[3]),
                             s10 * (u[1] - u[4]) * (u[1] - u[4]))
  powellj <- function(u, p) rbind(c(1, 10, 0, 0), c(0, 0, s5, -s5), c(0, 2 * (u[2] - 2 * u[3]), -4 * (u[2] - 2 * u[3]), 0),
                                  c(2 * s10 * (u[1] - u[4]), 0, 0, -2 * s10 * (u[1] - u[4])))
  tt <- c(1, 2, 4, 8, 16); yy <- c(0.9, 1.5, 2.1, 2.6, 2.75)
  mm <- function(u, p) u[1] * tt / (u[2] + tt) - yy
  mmj <- function(u, p) cbind(tt / (u[2] + tt), -u[1] * tt / ((u[2] + tt) * (u[2] + tt)))
  circ <- function(u, p) c(u[1] * u[1] + u[2] * u[2] - 4, u[1] * u[2] - 1)
  circj <- function(u, p) rbind(c(2 * u[1], 2 * u[2]), c(u[2], u[1]))
  list(
    rosen_in = list(NonlinearProblem, rosen, rosenj, c(-1.2, 1), c(-2, -2), c(2, 2)),
    rosen_bd = list(NonlinearProblem, rosen, rosenj, c(-1.2, 1), c(-2, -2), c(0.5, 2)),
    powell = list(NonlinearProblem, powell, powellj, c(3, -1, 0, 1), rep(-5, 4), rep(5, 4)),
    powell_bd = list(NonlinearProblem, powell, powellj, c(3, -1, 0, 1), c(0.5, -5, -5, -5), rep(5, 4)),
    mm_ls = list(NonlinearLeastSquaresProblem, mm, mmj, c(1, 1), c(0, 0), c(2.5, 10)),
    mm_ls_free = list(NonlinearLeastSquaresProblem, mm, mmj, c(1, 1), c(0, 0), c(10, 10)),
    circ = list(NonlinearProblem, circ, circj, c(3, 0.5), c(0, 0), c(3, 3)))
})

bounded_algs <- list(
  BGN_ls = function() BoundedGaussNewton(),
  BGN_tr = function() BoundedGaussNewton(globalization = "trustregion"),
  BGN_trls = function() BoundedGaussNewton(globalization = "trustregion_linesearch"),
  BLM = function() BoundedLevenbergMarquardt(),
  TRR = function() TrustRegionReflective(),
  Dogbox = function() Dogbox(),
  BTR = function() BoundedTrustRegion())

test_that("the bounded solvers match NonlinearSolve iterate for iterate", {
  for (line in strsplit(trimws(bounded_reference), "\n")[[1L]]) {
    f <- strsplit(line, " ")[[1L]]
    cs <- bounded_cases[[f[1L]]]
    prob <- cs[[1L]](NonlinearFunction(cs[[2L]], jac = cs[[3L]]), u0 = cs[[4L]], lb = cs[[5L]], ub = cs[[6L]])
    sol <- solve(prob, bounded_algs[[f[2L]]](), maxiters = as.integer(f[3L]))
    label <- paste(f[1:3], collapse = " ")
    expect_identical(sol$retcode, f[4L], label = label)
    expect_lt(max(abs(unname(sol$u) - as.numeric(strsplit(f[5L], ";")[[1L]]))), 1e-12, label = label)
    st <- sol$stats
    expect_identical(c(st$nf, st$njacs, st$nfactors, st$nsolve, st$nsteps), as.integer(f[6:10]), label = label)
  }
})

test_that("bounded solvers keep the iterates in the box and validate input", {
  f <- function(u, p) c(u[1] * u[1] - 2, u[2] - 1)
  prob <- NonlinearProblem(f, u0 = c(1, 0.5), lb = c(0, 0), ub = c(1.2, 2))
  for (alg in list(TrustRegionReflective(), Dogbox(), BoundedLevenbergMarquardt(), BoundedGaussNewton())) {
    sol <- solve(prob, alg)
    expect_true(all(sol$u >= c(0, 0) & sol$u <= c(1.2, 2)))
    expect_identical(sol$retcode, "Stalled")
  }
  expect_error(solve(NonlinearProblem(f, u0 = c(3, 0.5), lb = c(0, 0), ub = c(1.2, 2)), Dogbox()), "lower|upper")
  expect_error(BoundedGaussNewton(globalization = "none"), "globalization")
  expect_error(TrustRegionReflective(initial_trust_radius = 0), "positive")
})

test_that("bound aware finite differences never leave the box", {
  seen <- numeric()
  f <- function(u, p) { seen <<- c(seen, u[1]); c(sqrt(u[1]) - 0.5, u[2] - 1) }
  prob <- NonlinearProblem(f, u0 = c(0, 0), lb = c(0, 0), ub = c(4, 4))
  sol <- solve(prob, TrustRegionReflective(autodiff = AutoFiniteDiff()))
  expect_true(all(seen >= 0))
  expect_equal(unname(sol$u), c(0.25, 1), tolerance = 1e-8)
})
