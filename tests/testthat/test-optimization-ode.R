# OptimizationODE (0.1.12): gradient flows integrated to a steady state. The
# reference rows come from tools/optimizationode-probe.jl: problem, method,
# maxiters (0: 1e6 steps), retcode, u and objective (hexadecimal), and, for
# the ODE methods, iterations, fevals, gevals and time.
ode_ref <- c(
  "quad\tODEGradientDescent\t1\tSuccess\t0x1.ff7ced916872bp+0 0x1.ff7ced916872bp+0\t0x1.fdf43c7d5ed07p+0\t1\t2\t1\t0x0p+0",
  "quad\tODEGradientDescent\t5\tSuccess\t0x1.fd7341967a58bp+0 0x1.fd7341967a58bp+0\t0x1.f5da070ed7a5fp+0\t5\t6\t5\t0x0p+0",
  "quad\tODEGradientDescent\t50\tSuccess\t0x1.e79d7cca7cc63p+0 0x1.e79d7cca7cc63p+0\t0x1.a31b302f74bfp+0\t50\t51\t50\t0x0p+0",
  "quad\tRKChebyshevDescent\t1\tSuccess\t0x1.ff7d0f1d2aep+0 0x1.ff7d0f1d2aep+0\t0x1.fdf4c267bd952p+0\t1\t6\t1\t0x0p+0",
  "quad\tRKChebyshevDescent\t5\tSuccess\t0x1.fe723921cd0ddp+0 0x1.fe723921cd0ddp+0\t0x1.f9cdb8ac03843p+0\t5\t18\t5\t0x0p+0",
  "quad\tRKChebyshevDescent\t50\tSuccess\t0x1.e5fced8a81e64p+0 0x1.e5fced8a81e64p+0\t0x1.9d3cf5bca2fedp+0\t50\t155\t50\t0x0p+0",
  "quad\tRKAccelerated\t1\tSuccess\t0x1.ff7d0f199e6bbp+0 0x1.ff7d0f199e6bbp+0\t0x1.fdf4c25993068p+0\t1\t7\t1\t0x0p+0",
  "quad\tRKAccelerated\t5\tSuccess\t0x1.bbcd05ea9cbd1p+0 0x1.bbcd05ea9cbd1p+0\t0x1.138a55ae75bc8p+0\t5\t31\t5\t0x0p+0",
  "quad\tRKAccelerated\t50\tSuccess\t0x1.0000064b89aedp+0 0x1.0000064b89aedp+0\t0x1.3d05f141d39p-42\t35\t211\t35\t0x0p+0",
  "quad\tRKAccelerated\t0\tSuccess\t0x1.0000064b89aedp+0 0x1.0000064b89aedp+0\t0x1.3d05f141d39p-42\t35\t211\t35\t0x0p+0",
  "quad\tHighOrderDescent\t1\tSuccess\t0x1.ff7d0f199e6bbp+0 0x1.ff7d0f199e6bbp+0\t0x1.fdf4c25993068p+0\t1\t10\t1\t0x0p+0",
  "quad\tHighOrderDescent\t5\tSuccess\t0x1.7f29d0f243587p+0 0x1.7f29d0f243587p+0\t0x1.f9542127b7af2p-2\t5\t50\t5\t0x0p+0",
  "quad\tHighOrderDescent\t50\tSuccess\t0x1.000001586afd8p+0 0x1.000001586afd8p+0\t0x1.cf5fb5fee9064p-47\t22\t220\t22\t0x0p+0",
  "quad\tHighOrderDescent\t0\tSuccess\t0x1.000001586afd8p+0 0x1.000001586afd8p+0\t0x1.cf5fb5fee9064p-47\t22\t220\t22\t0x0p+0",
  "rosen\tODEGradientDescent\t1\tSuccess\t0x1.595810624dd2fp+0 0x1.0666666666666p+1\t0x1.5af0fdd9a59bcp+2\t1\t2\t1\t0x0p+0",
  "rosen\tODEGradientDescent\t5\tSuccess\t0x1.63f799fcfac28p+0 0x1.03be1ecdeb83ap+1\t0x1.11c7532cfb582p+0\t5\t6\t5\t0x0p+0",
  "rosen\tODEGradientDescent\t50\tSuccess\t0x1.6a75669333542p+0 0x1.00c2ae5caef7p+1\t0x1.6283ec7f2ea6fp-3\t50\t51\t50\t0x0p+0",
  "rosen\tRKChebyshevDescent\t1\tSuccess\t0x1.8p+0 0x1p+1\t0x1.ap+2\t0\t6\t0\t0x0p+0",
  "rosen\tRKChebyshevDescent\t5\tSuccess\t0x1.8p+0 0x1p+1\t0x1.ap+2\t0\t26\t0\t0x0p+0",
  "rosen\tRKChebyshevDescent\t50\tSuccess\t0x1.79d801fe8618dp+0 0x1.0106b8519d74bp+1\t0x1.90bc34abbf85ap+1\t45\t165\t45\t0x0p+0",
  "rosen\tRKAccelerated\t1\tSuccess\t0x1.8p+0 0x1p+1\t0x1.ap+2\t0\t7\t0\t0x0p+0",
  "rosen\tRKAccelerated\t5\tSuccess\t0x1.72455f7d93cep+0 0x1.024df76a1d19bp+1\t0x1.7e45616cf9396p-1\t3\t31\t3\t0x0p+0",
  "rosen\tRKAccelerated\t50\tSuccess\t0x1.6ab9ec9d48389p+0 0x1.0123fde21e8fcp+1\t0x1.644d7ea2f8ef6p-3\t48\t301\t48\t0x0p+0",
  "rosen\tHighOrderDescent\t1\tSuccess\t0x1.8p+0 0x1p+1\t0x1.ap+2\t0\t10\t0\t0x0p+0",
  "rosen\tHighOrderDescent\t5\tSuccess\t0x1.713c8b90fb739p+0 0x1.027ac6dd81c38p+1\t0x1.2257cdaffff36p-1\t2\t50\t2\t0x0p+0",
  "rosen\tHighOrderDescent\t50\tSuccess\t0x1.69f70a3f9b57ep+0 0x1.000ff6b23c4f1p+1\t0x1.5f3d660fad8abp-3\t47\t500\t47\t0x0p+0",
  "rosen\tODEGradientDescent2\t1\tSuccess\t0x1.32b020c49ba5ep+0 0x1.0cccccccccccdp+1\t0x1.61e0697773f24p+5\t1\t2\t1\t0x0p+0",
  "rosen\tODEGradientDescent2\t5\tSuccess\t-0x1.52d1b15974421p+0 0x1.7716924b13c3bp+0\t0x1.b36a3ff02c155p+3\t5\t6\t5\t0x0p+0",
  "rosen\tODEGradientDescent2\t50\tSuccess\t-0x1.d29af7677b3d3p-1 0x1.ad866ea5eedb8p-1\t0x1.d48291bcc5402p+1\t50\t51\t50\t0x0p+0",
  "quad\tRKAcceleratednodt\t0\tSuccess\t0x1.0000046249ba8p+0 0x1.0000046249ba8p+0\t0x1.3380a5c05e4dep-43\t34\t207\t34\t0x0p+0",
  "quad\tHighOrderDescentnodt\t0\tSuccess\t0x1.000001c4b8483p+0 0x1.000001c4b8483p+0\t0x1.904da1c8acad8p-46\t20\t202\t20\t0x0p+0",
  "cons\tDAEMassMatrix\t5\tFailure\t0x1.794efaf25b39bp-1 0x1.0d620a1b498cdp-2\t0x1.397bc196fbea7p-1",
  "cons\tDAEMassMatrix\t50\tSuccess\t0x1.0b3d3dc5a90dbp-21 0x1.ffffef4c2c23bp-1\t0x1.ffffde98595e6p-1",
  "cons\tDAEMassMatrix\t0\tSuccess\t0x1.0b3d3dc5a90dbp-21 0x1.ffffef4c2c23bp-1\t0x1.ffffde98595e6p-1"
)

odeopt_quad <- function(x, p) (x[1] - p[1])^2 + (x[2] - p[2])^2
odeopt_rosen <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
odeopt_probs <- list(
  quad = OptimizationProblem(odeopt_quad, c(2, 2), c(1, 1),
                             grad = function(x, p) c(2.0 * (x[1] - p[1]), 2.0 * (x[2] - p[2]))),
  rosen = OptimizationProblem(odeopt_rosen, c(1.5, 2), c(1, 100),
                              grad = function(x, p) c(-2.0 * (p[1] - x[1]) - 4.0 * p[2] * (x[2] - x[1]^2) * x[1],
                                                      2.0 * p[2] * (x[2] - x[1]^2))),
  cons = OptimizationProblem(function(x, p) x[1]^2 + x[2]^2, c(1, 1), 1, grad = function(x, p) c(2.0 * x[1], 2.0 * x[2]),
                             cons = function(x, p) x[1] + x[2] - 1, lcons = -10, ucons = 10))
odeopt_methods <- list(
  ODEGradientDescent = ODEGradientDescent, ODEGradientDescent2 = ODEGradientDescent,
  RKChebyshevDescent = RKChebyshevDescent, RKAccelerated = RKAccelerated, HighOrderDescent = HighOrderDescent,
  RKChebyshevDescentnodt = RKChebyshevDescent, RKAcceleratednodt = RKAccelerated,
  HighOrderDescentnodt = HighOrderDescent, DAEMassMatrix = DAEMassMatrix)

test_that("the ODE optimizers follow the integrations of OptimizationODE", {
  for (row in strsplit(ode_ref, "\t")) {
    n <- as.integer(row[3])
    info <- paste(row[1:3], collapse = " ")
    args <- list(odeopt_probs[[row[1]]], odeopt_methods[[row[2]]]())
    if (!endsWith(row[2], "nodt")) {
      args$dt <- if (row[2] == "ODEGradientDescent2") 0.002 else if (row[1] == "cons") 0.01 else 0.001
      args$maxiters <- if (n > 0L) n else 1e6
    }
    s <- suppressWarnings(do.call(solve, args))
    ju <- as.numeric(strsplit(row[5], " ")[[1]])
    jo <- as.numeric(row[6])
    rel <- max(abs(unname(as.numeric(s$u)) - ju) / pmax(1, abs(ju)), abs(s$objective - jo) / max(1, abs(jo)))
    # Euler and ROCK2 agree to the last bits. The first steps of Tsit5 and
    # Vern7 from dt = 0.001 have error estimates at the level of rounding.
    # OptimizationODE builds an in place problem, whose Tsit5 fuses the
    # multiply adds of its broadcasts; the R Tsit5 follows the out of place
    # method, whose step sizes it matches here to the last bits. So the next
    # step sizes and the iterates after 5 steps differ by about 1e-4 (the
    # converged results by less than 1e-9). The in place and out of place
    # Vern7 of Julia differ by 10% after 2 steps; the R Vern7 follows the
    # out of place sums, and its iterates differ from the in place ones by
    # 1e-2 after 5 steps and 6e-9 at convergence. The Rodas5P iterates
    # differ by 1.3e-6.
    hod <- row[2] == "HighOrderDescent"
    tol <- if (row[1] == "quad" && n == 5L && hod) 2e-2
           else if (row[1] == "quad" && n == 5L && row[2] == "RKAccelerated") 2e-4
           else if (row[1] == "cons") 2e-6 else if (hod) 1e-8 else 1e-9
    expect_lt(rel, ref_tol(tol, max(tol, 1e-4)), label = info)
    expect_identical(s$retcode, row[4], info = info)
    if (length(row) > 6L) {
      expect_identical(c(s$stats$iterations, s$stats$fevals, s$stats$gevals),
                       as.integer(row[c(7L, 8L, 9L)]), info = info)
    }
  }
})

test_that("DAEOptimizer with an implicit DAE solver follows the gradient flow", {
  # the consimpl rows of tools/optimizationode-probe.jl: the DAE of
  # solve_dae_implicit with the residual du + grad f (OptimizationODE writes
  # du - grad f, which diverges; see R/optimization-ode.R). DFBDF agrees to
  # the last bits for 5 steps; the converged points agree to 1e-8 (the R and
  # Julia DImplicitEuler choose different steps)
  ref <- c("DFBDF\t5\tMaxIters\t0x1.553691b7f9bf1p-2 0x1.5564b72403207p-1",
           "DFBDF\t0\tSuccess\t-0x1.001df886f9202p-26 0x1.00000040077e2p+0",
           "DImplicitEuler\t0\tSuccess\t0x1.f7e0731e42c94p-25 0x1.fffffe081f8cep-1")
  for (row in strsplit(ref, "\t")) {
    n <- as.integer(row[2])
    alg <- DAEOptimizer(get(row[1])())
    s <- suppressWarnings(if (n == 0L) solve(odeopt_probs$cons, alg, dt = 0.01)
                          else solve(odeopt_probs$cons, alg, dt = 0.01, maxiters = n))
    ju <- as.numeric(strsplit(row[4], " ")[[1]])
    expect_lt(max(abs(as.numeric(s$u) - ju)), if (n == 5L) 1e-15 else 5e-8, label = paste(row[1:2], collapse = " "))
    expect_identical(s$retcode, row[3])
  }
})

test_that("the ODE optimizers check their problems", {
  expect_s3_class(ODEOptimizer(Tsit5()), "ODEOptimizer")
  expect_s3_class(DAEOptimizer(Rodas5P()), "DAEOptimizer")
  bprob <- OptimizationProblem(odeopt_quad, c(2, 2), c(1, 1), lb = c(0, 0), ub = c(3, 3))
  expect_error(solve(bprob, RKAccelerated()), class = "IncompatibleOptimizerError")
  expect_error(solve(odeopt_probs$cons, RKAccelerated()), class = "IncompatibleOptimizerError")
  expect_error(solve(odeopt_probs$quad, DAEMassMatrix()), "requires constraints")
  # a callback that stops the integration
  calls <- 0L
  cb <- function(state, l) {
    calls <<- calls + 1L
    calls >= 3L
  }
  s <- solve(odeopt_probs$quad, ODEGradientDescent(), dt = 0.001, maxiters = 100, callback = cb)
  expect_identical(calls, 3L)
  expect_identical(s$retcode, "Terminated")
})
