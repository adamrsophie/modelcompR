# SpeedMappingOpt (OptimizationSpeedMapping 0.2.9 with SpeedMapping 0.4.1).
# The reference rows come from tools/optimizationspeedmapping-probe.jl:
# problem, maxiters (0: the default), retcode, u and objective (hexadecimal),
# maps, f_calls, iterations, status, residual norm and last learning rate
# (hexadecimal).
speed_ref <- c(
  "rosen\t5\tMaxIters\t-0x1.075aeb0b91a88p+0 0x1.10f878916d7f1p+0\t0x1.07d166533e4b8p+2\t6\t5\t2\tmax_eval\t0x1.330eba1df2001p+1\t0x1.27cdd3527c4d5p-11",
  "rosen\t10\tMaxIters\t0x1.0c8f0c18f8ab2p+2 -0x1.c0246ed0c3ddp+2\t0x1.d944898941788p+15\t12\t5\t4\tmax_eval\t0x1.c65dcb00c6d91p+0\t0x1.c024a02734a34p-11",
  "rosen\t20\tMaxIters\t-0x1.47b80507962f8p-5 -0x1.258941f769b35p+2\t0x1.0744d0835fb8fp+11\t20\t5\t7\tmax_eval\t0x1.33185a9b71c4dp+10\t0x1.6683b3529082ap-13",
  "rosen\t50\tMaxIters\t0x1.c0d3ac9d27a2fp-1 0x1.892d20aa737c1p-1\t0x1.f3c7ebe4d7a93p-7\t50\t5\t18\tmax_eval\t0x1.f580b17ff5337p-4\t0x1.d37e685e024ccp-11",
  "rosen\t0\tSuccess\t0x1.0000000006034p+0 0x1.fffffffff9eddp-1\t0x1.62fd0d1c44p-66\t85\t5\t32\tfirst_order\t0x1.a5a33fc91837cp-28\t0x1.65aaf52fb4d18p-10",
  "rosenbox\t5\tMaxIters\t0x1.5893197d3627fp-3 0x1.a6b40d7643f32p-6\t0x1.628743aaaf6a9p-1\t6\t3\t2\tmax_eval\t0x1.998bc9c1f0b8p+0\t0x1.295449bfb51ccp-8",
  "rosenbox\t10\tMaxIters\t0x1.b1a8d9d2e8fcap-2 0x1.6b9abee658b16p-3\t0x1.54ab1565ad219p-2\t11\t3\t4\tmax_eval\t0x1.00e17b630a43ep+3\t0x1.295449bfb51ccp-8",
  "rosenbox\t20\tMaxIters\t0x1.fcb83f4cbc2b1p-2 0x1.f98e839b9423cp-3\t0x1.034a802dd0625p-2\t22\t3\t8\tmax_eval\t0x1.12aab1a2ab5fap+1\t0x1.295449bfb51ccp-8",
  "rosenbox\t50\tSuccess\t0x1p-1 0x1p-2\t0x1p-2\t36\t3\t14\tfirst_order\t0x0p+0\t0x1.94684cf60947ep-9",
  "rosenbox\t0\tSuccess\t0x1p-1 0x1p-2\t0x1p-2\t36\t3\t14\tfirst_order\t0x0p+0\t0x1.94684cf60947ep-9",
  "rosenlow\t5\tMaxIters\t-0x1.38952ce9c5c1ap+0 0x1.8p+0\t0x1.3c3d19e4b3336p+2\t6\t5\t2\tmax_eval\t0x1.b5615bcfe7p-6\t0x1.27cdd3527c4d5p-11",
  "rosenlow\t10\tSuccess\t-0x1.38952d0168644p+0 0x1.8p+0\t0x1.3c3d19e4b3322p+2\t10\t5\t4\tfirst_order\t0x1.8p-49\t0x1.2ac3156f786cdp-11",
  "rosenlow\t20\tSuccess\t-0x1.38952d0168644p+0 0x1.8p+0\t0x1.3c3d19e4b3322p+2\t10\t5\t4\tfirst_order\t0x1.8p-49\t0x1.2ac3156f786cdp-11",
  "rosenlow\t50\tSuccess\t-0x1.38952d0168644p+0 0x1.8p+0\t0x1.3c3d19e4b3322p+2\t10\t5\t4\tfirst_order\t0x1.8p-49\t0x1.2ac3156f786cdp-11",
  "rosenlow\t0\tSuccess\t-0x1.38952d0168644p+0 0x1.8p+0\t0x1.3c3d19e4b3322p+2\t10\t5\t4\tfirst_order\t0x1.8p-49\t0x1.2ac3156f786cdp-11",
  "quad3\t5\tMaxIters\t0x1.1259719c40656p+1 0x1.ffffffffff735p-1 -0x1.0743d824d6bfcp-1\t0x1.274d77331368ep+2\t6\t3\t2\tmax_eval\t0x1.451ed821e9301p+5\t0x1.b993fe00d5376p-8",
  "quad3\t10\tMaxIters\t0x1.bp-46 0x1.ffffffffdfaa1p-1 -0x1.fffffffffd54cp-2\t0x1.4c5070f6dp-69\t11\t3\t4\tmax_eval\t0x1.07cf3430903d4p+2\t0x1.b993fe00d5376p-8",
  "quad3\t20\tSuccess\t0x1.aa2dac86bd306p-46 0x1.ffffffffe4059p-1 -0x1.000000000076cp-1\t0x1.ea95b53090296p-70\t12\t3\t5\tfirst_order\t0x1.1b9a1d13ae995p-32\t0x1.b993fe00d5376p-8",
  "quad3\t50\tSuccess\t0x1.aa2dac86bd306p-46 0x1.ffffffffe4059p-1 -0x1.000000000076cp-1\t0x1.ea95b53090296p-70\t12\t3\t5\tfirst_order\t0x1.1b9a1d13ae995p-32\t0x1.b993fe00d5376p-8",
  "quad3\t0\tSuccess\t0x1.aa2dac86bd306p-46 0x1.ffffffffe4059p-1 -0x1.000000000076cp-1\t0x1.ea95b53090296p-70\t12\t3\t5\tfirst_order\t0x1.1b9a1d13ae995p-32\t0x1.b993fe00d5376p-8",
  "quad3opt\t5\tSuccess\t0x0p+0 0x1p+0 -0x1p-1\t0x0p+0\t1\t0\t0\tfirst_order\t0x0p+0\t0x1p+0",
  "quad3opt\t10\tSuccess\t0x0p+0 0x1p+0 -0x1p-1\t0x0p+0\t1\t0\t0\tfirst_order\t0x0p+0\t0x1p+0",
  "quad3opt\t20\tSuccess\t0x0p+0 0x1p+0 -0x1p-1\t0x0p+0\t1\t0\t0\tfirst_order\t0x0p+0\t0x1p+0",
  "quad3opt\t50\tSuccess\t0x0p+0 0x1p+0 -0x1p-1\t0x0p+0\t1\t0\t0\tfirst_order\t0x0p+0\t0x1p+0",
  "quad3opt\t0\tSuccess\t0x0p+0 0x1p+0 -0x1p-1\t0x0p+0\t1\t0\t0\tfirst_order\t0x0p+0\t0x1p+0"
)

speed_rosen <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
speed_roseng <- function(x, p) c(-2.0 * (p[1] - x[1]) - 4.0 * p[2] * (x[2] - x[1]^2) * x[1], 2.0 * p[2] * (x[2] - x[1]^2))
speed_quad3 <- function(x, p) x[1]^2 + 10.0 * (x[2] - 1.0)^2 + 100.0 * (x[3] + 0.5)^2
speed_quad3g <- function(x, p) c(2.0 * x[1], 20.0 * (x[2] - 1.0), 200.0 * (x[3] + 0.5))
speed_probs <- list(
  rosen = OptimizationProblem(speed_rosen, c(-1.2, 1), c(1, 100), grad = speed_roseng),
  rosenbox = OptimizationProblem(speed_rosen, c(0.1, 0.2), c(1, 100), grad = speed_roseng, lb = c(-0.5, -0.5),
                                 ub = c(0.5, 0.8)),
  rosenlow = OptimizationProblem(speed_rosen, c(-1.2, 1), c(1, 100), grad = speed_roseng, lb = c(-2, 1.5),
                                 ub = c(2, 3)),
  quad3 = OptimizationProblem(speed_quad3, c(3, -2, 1), NULL, grad = speed_quad3g),
  quad3opt = OptimizationProblem(speed_quad3, c(0, 1, -0.5), NULL, grad = speed_quad3g))

test_that("SpeedMappingOpt follows the iterates of SpeedMapping.jl", {
  for (row in strsplit(speed_ref, "\t")) {
    n <- as.integer(row[2])
    info <- paste(row[1:2], collapse = " ")
    s <- suppressMessages(if (n > 0L) solve(speed_probs[[row[1]]], SpeedMappingOpt(), maxiters = n)
                          else solve(speed_probs[[row[1]]], SpeedMappingOpt()))
    expect_identical(unname(as.numeric(s$u)), as.numeric(strsplit(row[4], " ")[[1]]), info = info)
    expect_identical(s$objective, as.numeric(row[5]), info = info)
    expect_identical(s$retcode, row[3], info = info)
    o <- s$original
    expect_identical(c(o$maps, o$f_calls, o$iterations), as.integer(row[6:8]), info = info)
    expect_identical(o$status, row[9], info = info)
    expect_identical(c(o$residual_norm, o$last_learning_rate), as.numeric(row[10:11]), info = info)
  }
})

test_that("SpeedMappingOpt reports its options and checks its problems", {
  expect_message(solve(speed_probs$quad3, SpeedMappingOpt(), maxiters = 5), "maximum gradient calls")
  expect_message(solve(speed_probs$quad3, SpeedMappingOpt(), abstol = 1e-3), "abstol is currently not used")
  cprob <- OptimizationProblem(speed_rosen, c(0.1, 0.2), c(1, 100), grad = speed_roseng,
                               cons = function(x, p) x[1] + x[2], lcons = 0, ucons = 1)
  expect_error(solve(cprob, SpeedMappingOpt()), class = "IncompatibleOptimizerError")
  # without a gradient, central differences take its place
  s <- solve(OptimizationProblem(speed_quad3, c(3, -2, 1), NULL), SpeedMappingOpt())
  expect_identical(s$retcode, "Success")
  expect_equal(unname(as.numeric(s$u)), c(0, 1, -0.5), tolerance = 1e-6)
})
