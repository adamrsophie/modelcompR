# BAOAB and the IIF methods on fixed noise paths. Reference values from
# tools/stochasticdiffeq-probe.jl (StochasticDiffEq); the IIF1Mil values come
# from the in place step with f1(u, p, t), since upstream calls the old
# f1(t, u, du) there (a MethodError) and its out of place step stops.

test_that("BAOAB matches StochasticDiffEq", {
  gr <- sdeprobe_grid(2, 1 / 64, 80, off = 11)
  mk <- function() DynamicalSDEProblem(function(v, u, p, t) -u - 0.1 * u^3, function(v, u, p, t) v,
                                       function(v, u, p, t) 0.3 + 0.1 * u, c(0.5, -0.2), c(1, 0.3), c(0, 1),
                                       noise = NoiseGrid(gr$t, gr$W))
  sol <- solve(mk(), BAOAB(gamma = 0.7), dt = 1 / 16)
  expect_length(sol$t, 17L)
  expect_equal(unlist(sol$u[c(2, 9, 17)]),
               c(0.4101534733546766, -0.16478201273874563, 1.0284786904936123, 0.288589129099692,
                 -0.1251661519200029, -0.2629833934751265, 1.0843909298206758, 0.20335161819839473,
                 -0.5119625252985771, -0.21638654050637549, 0.9200846654139105, 0.06664879009158717),
               tolerance = 1e-13)
  sol <- solve(mk(), BAOAB(gamma = 0.7, scale_noise = FALSE), dt = 1 / 16)
  expect_equal(unlist(sol$u[c(2, 9, 17)]),
               c(0.4101534733546766, -0.1708824123500849, 1.0284786904936123, 0.28839830060283866,
                 -0.12861229837333593, -0.26027784735662846, 1.0838292511601388, 0.20099220118344768,
                 -0.5196588340380462, -0.21725937756675814, 0.9158669950852719, 0.06673552873092928),
               tolerance = 1e-13)
  bad <- DynamicalSDEProblem(function(v, u, p, t) -u, function(v, u, p, t) 2 * v, function(v, u, p, t) 1,
                             0.5, 1, c(0, 1))
  expect_error(solve(bad, BAOAB(), dt = 0.1, seed = 1), "Algorithm BAOAB is not applicable if f2(p, q, t) != p",
               fixed = TRUE)
  expect_error(solve(SDEProblem(function(u, p, t) -u, function(u, p, t) 1, 1, c(0, 1)), BAOAB(), dt = 0.1),
               "DynamicalSDEProblem")
})

test_that("IIF1M and IIF2M match StochasticDiffEq for scalar noise", {
  h <- 1 / 64
  W <- cumsum(c(0, sqrt(h) * sdeprobe_z((1:80) + 23)))
  mk <- function() SplitSDEProblem(function(u, p, t) -1.3, function(u, p, t) 0.4 * u + 0.1 * t,
                                   function(u, p, t) 0.3 * u, 0.8, c(0, 1), noise = NoiseGrid((0:80) * h, W))
  s1 <- solve(mk(), IIF1M(), dt = 1 / 16)
  expect_equal(unlist(s1$u)[c(2, 8, 12, 17)],
               c(0.7790103596834107, 0.5834325879724453, 0.47422308845092875, 0.38103282366441765),
               tolerance = 1e-13)
  s2 <- solve(mk(), IIF2M(), dt = 1 / 16)
  expect_equal(unlist(s2$u)[c(2, 8, 12, 17)],
               c(0.778288020161726, 0.5801338542127465, 0.46971556197090725, 0.37500131110074875),
               tolerance = 1e-13)
  # IIF1Mil: upstream's out of place step stops ("Milstein correction does
  # not work."); here the scalar step is the in place formula
  s3 <- solve(mk(), IIF1Mil(), dt = 1 / 16)
  dt <- 1 / 16
  dW <- W[5] - W[1]
  M <- exp(-1.3 * dt)
  r2 <- M * 0.3 * 0.8
  Dg <- (M * 0.3 * (0.8 + sqrt(dt) * r2) - r2) / sqrt(dt)
  tmp <- M * 0.8 + r2 * dW + Dg * (dW^2 - dt) / 2
  expect_equal(s3$u[[2]], (tmp + dt * 0.1 * dt) / (1 - 0.4 * dt), tolerance = 1e-12)
  expect_error(solve(SDEProblem(function(u, p, t) -u, function(u, p, t) 1, 1, c(0, 1)), IIF1M(), dt = 0.1),
               "SplitSDEProblem")
})

test_that("IIF methods match StochasticDiffEq for nondiagonal noise", {
  gr <- sdeprobe_grid(2, 1 / 64, 80, off = 5)
  A <- matrix(c(-2, 1, 1, -2), 2, 2)
  mk <- function() SplitSDEProblem(function(u, p, t) A, function(u, p, t) 1.01 * u + 0.05 * u^2,
                                   function(u, p, t) cbind(0.87 * u, 0.2 * u^2), c(0.3, 0.6), c(0, 1),
                                   noise = NoiseGrid(gr$t, gr$W), noise_rate_prototype = matrix(0, 2, 2))
  pick <- function(sol) unlist(sol$u[c(2, 9, 17)])
  expect_equal(pick(solve(mk(), IIF1M(), dt = 1 / 16)),
               c(0.3316481897568543, 0.6061561720617155, 0.4860887905909153, 0.6227375648772285,
                 0.5596400540220217, 0.6135418838296022), tolerance = 1e-12)
  expect_equal(pick(solve(mk(), IIF2M(), dt = 1 / 16)),
               c(0.33054580442814363, 0.6041560253314623, 0.4753683612212472, 0.6088608594347621,
                 0.5373785560090096, 0.5890021333947493), tolerance = 1e-12)
  expect_equal(pick(solve(mk(), IIF1Mil(), dt = 1 / 16)),
               c(0.32431802059995457, 0.5932370312763131, 0.41055146451228086, 0.5274278131307271,
                 0.3964847169338756, 0.4354347626854535), tolerance = 1e-12)
})

test_that("NLSOLVEJL_SETUP solves the IIF stage equation", {
  s <- NLSOLVEJL_SETUP()
  expect_s3_class(s, "NLSOLVEJL_SETUP")
  expect_equal(s(IIFNLSolveFunc(function(u) u^2 - 2), 1), sqrt(2), tolerance = 1e-12)
})
