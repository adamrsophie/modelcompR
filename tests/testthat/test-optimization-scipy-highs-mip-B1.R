# HiGHS 1.12.0 MIP driver (group B1): the shared MIP utilities, the
# HighsObjectiveFunction, and the HighsMipSolver and HighsMipSolverData
# pieces that run without the other MIP classes. Expected values are worked
# by hand from the C++ (mip/HighsMipSolverData.cpp, HighsObjectiveFunction.cpp,
# util/HighsIntegers.h, pdqsort.h, libstdc++ stl_algo.h).

test_that("pdqsort, partial_sort and partition follow the C++ algorithms", {
  set.seed(1)
  for (n in c(0, 1, 5, 23, 24, 50, 129, 300, 1000)) {
    v <- sample.int(1000L, n, replace = TRUE)
    expect_identical(.hmip_pdqsort(v, function(a, b) a < b), sort(v))
    expect_identical(.hmip_pdqsort(v, function(a, b) a < b, branchless = TRUE), sort(v))
  }
  v <- sample.int(100L, 40L)
  expect_identical(.hmip_std_partial_sort(v, 7L, function(a, b) a < b)[1:7], sort(v)[1:7])
  # libstdc++ bidirectional std::partition
  p <- .hmip_std_partition(1:10, function(i) i %% 2L == 0L)
  expect_identical(p$v, c(10L, 2L, 8L, 4L, 6L, 5L, 7L, 3L, 9L, 1L))
  expect_identical(p$k, 5L)
})

test_that(".hmip_round is std::round", {
  x <- c(0.49999999999999994, -0.49999999999999994, 2.5, -2.5, 2^52 - 0.5, 0.5, -0.5, 3, NaN, Inf)
  expect_identical(.hmip_round(x), c(0, -0, 3, -3, 2^52, 1, -1, 3, NaN, Inf))
  expect_identical(.hmip_fractionality(0.49999999999999994), 0.49999999999999994)
})

test_that("integralScale and nextafter match HighsIntegers", {
  expect_identical(.hmip_integralScale(c(0.5, 1.5, 3), 1e-9, 1e-9), 2)
  expect_identical(.hmip_integralScale(c(1, 2, 4), 1e-9, 1e-9), 1)
  expect_identical(.hmip_integralScale(c(1 / 3, 1), 1e-9, 1e-9), 3)
  expect_identical(.hmip_nextafter_down(1), 1 - 2^-53)
  expect_identical(.hmip_nextafter_down(3), 3 - 2^-51)
})

mip_lp <- function() {
  # min -x0 - 2 x1 + 0.5 x2 + 3 x3 with x0, x1 binary, x2 integer in [0, 5],
  # x3 continuous; rows x0 + x1 + x2 <= 4 and x1 - x3 >= -1
  A <- rbind(c(1, 1, 1, 0), c(0, 1, 0, -1))
  .hmip_HighsLp_from_dense(c(-1, -2, 0.5, 3), A, c(-Inf, -1), c(4, Inf),
                           c(0, 0, 0, 0), c(1, 1, 5, 10), integrality = c(1L, 1L, 1L, 0L))
}

test_that("HighsObjectiveFunction partitions the objective columns as C++", {
  ms <- .hmip_new_HighsMipSolver(list(), mip_lp())
  of <- .hmip_new_HighsObjectiveFunction(ms)
  # std::partition integral first, then binaries first among those
  expect_identical(of$numIntegral, 3L)
  expect_identical(of$numBinary, 2L)
  expect_identical(sort(of$objectiveNonzeros[1:2]), c(0L, 1L))
  expect_identical(of$objectiveNonzeros[4], 3L)
  expect_identical(of$objectiveVals, ms$model_$col_cost_[of$objectiveNonzeros + 1L])
  .hmip_HighsObjectiveFunction_checkIntegrality(of, 1e-9)
  expect_false(.hmip_HighsObjectiveFunction_isIntegral(of))  # continuous column in the objective
  lp <- mip_lp(); lp$col_cost_[4] <- 0
  of2 <- .hmip_new_HighsObjectiveFunction(.hmip_new_HighsMipSolver(list(), lp))
  .hmip_HighsObjectiveFunction_checkIntegrality(of2, 1e-9)
  expect_identical(.hmip_HighsObjectiveFunction_integralScale(of2), 2)
  lp$col_cost_[] <- 0
  of3 <- .hmip_new_HighsObjectiveFunction(.hmip_new_HighsMipSolver(list(), lp))
  expect_true(.hmip_HighsObjectiveFunction_isZero(of3))
})

test_that("solutionFeasible measures violations and the objective in HighsCDouble", {
  ms <- .hmip_new_HighsMipSolver(list(), mip_lp())
  r <- .hmip_HighsMipSolver_solutionFeasible(ms, ms$model_, c(1, 1, 2, 0))
  expect_true(r$feasible)
  expect_identical(.cd_val(r$obj), -2)
  r <- .hmip_HighsMipSolver_solutionFeasible(ms, ms$model_, c(1, 1, 2.5, 0))
  expect_false(r$feasible)
  expect_identical(r$integrality_violation, 0.5)
  expect_identical(r$row_violation, 0.5)
  r <- .hmip_HighsMipSolver_solutionFeasible(ms, ms$model_, c(1, 0, 0, 11))
  expect_identical(r$bound_violation, 1)
  # an initial solution makes the constructor record it
  ms2 <- .hmip_new_HighsMipSolver(list(), mip_lp(), list(value_valid = TRUE, col_value = c(1, 1, 2, 0)))
  expect_identical(ms2$solution_objective_, -2)
})

limits_md <- function(offset = 0, cost = c(-1, -2, 0.5, 0)) {
  lp <- mip_lp(); lp$offset_ <- offset; lp$col_cost_ <- cost
  ms <- .hmip_new_HighsMipSolver(list(), lp)
  md <- new.env(parent = emptyenv())
  md$mipsolver <- ms; md$epsilon <- 1e-9; md$feastol <- 1e-6
  md$objectiveFunction <- .hmip_new_HighsObjectiveFunction(ms)
  .hmip_HighsObjectiveFunction_checkIntegrality(md$objectiveFunction, 1e-9)
  md
}

test_that("limitsToGap and computeNewUpperLimit follow HighsMipSolverData", {
  md <- limits_md(offset = 1)
  g <- .hmip_HighsMipSolverData_limitsToGap(md, -5, -3)
  expect_identical(c(g$lb, g$ub, g$gap), c(-4, -2, 1))
  g <- .hmip_HighsMipSolverData_limitsToGap(md, -5, Inf)
  expect_identical(c(g$lb, g$ub, g$gap), c(-4, Inf, Inf))
  g <- .hmip_HighsMipSolverData_limitsToGap(md, -1, -1)
  expect_identical(c(g$lb, g$ub, g$gap), c(0, 0, 0))
  # integral objective with scale 2: floor(2 * ub - 0.5) / 2 + feastol
  expect_identical(.hmip_HighsMipSolverData_computeNewUpperLimit(md, -3, 0, 0), -3.5 + 1e-6)
  # relative gap: ub - ceil(1e-4 * |ub + offset| * 2 - eps) / 2
  expect_identical(.hmip_HighsMipSolverData_computeNewUpperLimit(md, -3, 1e-6, 1e-4),
                   min(-3.5, -3 - 0.5, -3 - 0.5) + 1e-6)
  md2 <- limits_md(cost = c(-1, -2, 0.5, 0.3))
  expect_identical(.hmip_HighsMipSolverData_computeNewUpperLimit(md2, -3, 0, 0), -3 - 1e-6)
  expect_identical(.hmip_HighsMipSolverData_computeNewUpperLimit(md2, -3, 1e-6, 1e-4),
                   min(-3 - 1e-6, -3 - 1e-4 * 3, -3 - 1e-6))
})

test_that("checkLimits sets the model status of the first limit hit", {
  md <- limits_md()
  ms <- md$mipsolver
  md$num_nodes <- 0; md$num_leaves <- 0; md$numImprovingSols <- 0L
  expect_false(.hmip_HighsMipSolverData_checkLimits(md))
  ms$options_mip_$mip_max_nodes <- 3
  md$num_nodes <- 3
  expect_true(.hmip_HighsMipSolverData_checkLimits(md))
  expect_identical(ms$modelstatus_, "SolutionLimit")
  ms$modelstatus_ <- "Notset"
  ms$options_mip_$mip_max_nodes <- .hmip_kHighsIInf
  ms$solution_objective_ <- -2; ms$options_mip_$objective_target <- -1
  expect_true(.hmip_HighsMipSolverData_checkLimits(md))
  expect_identical(ms$modelstatus_, "ObjectiveTarget")
})

test_that("withoutSemiVariables adds a binary and two rows per semi variable", {
  lp <- .hmip_HighsLp_from_dense(c(1, 2), matrix(c(1, 1), 1), 2, Inf, c(0, 1), c(5, 4),
                                 integrality = c(0L, 2L))
  out <- .hmip_withoutSemiVariables(lp)
  expect_identical(out$num_col_, 3L); expect_identical(out$num_row_, 3L)
  expect_identical(out$integrality_, c(0L, 0L, 1L))
  expect_identical(out$col_lower_, c(0, 0, 0)); expect_identical(out$col_upper_, c(5, 4, 1))
  expect_identical(out$row_lower_, c(2, 0, -Inf)); expect_identical(out$row_upper_, c(Inf, Inf, 0))
  A <- .hmip_HighsLp_dense(out)
  expect_identical(A, rbind(c(1, 1, 0), c(0, 1, -1), c(0, 1, -4)))
})

test_that("row values are HighsCDouble sums in column order", {
  lp <- .hmip_HighsLp_from_dense(c(0, 0, 0), matrix(c(1e16, 1, -1e16), 1), -Inf, Inf, c(0, 0, 0), c(1, 1, 1))
  expect_identical(.hmip_row_values_quad(lp, c(1, 1, 1)), 1)
  tr <- .hmip_sparse_transpose(2L, 3L, c(0L, 1L, 3L, 4L), c(1L, 0L, 1L, 0L), c(1, 2, 3, 4))
  expect_identical(tr$ARstart_, c(0L, 2L, 4L))
  expect_identical(tr$ARindex_, c(1L, 2L, 0L, 1L))
  expect_identical(tr$ARvalue_, c(2, 4, 1, 3))
})

test_that("the display line formats numbers as convertToPrintString", {
  expect_identical(.hmip_print_int(19), "19")
  expect_identical(.hmip_print_int(1234567), "1234k")
  expect_identical(.hmip_print_double(-12.5), "-12.5")
  expect_identical(.hmip_print_double(Inf), "inf")
})

# linprog(method = "highs-ipm"): case = list(lp, variant, status, nit, x)
# from tools/highs-probe.py (scratchpad highs-probe-linux.txt); x NULL when
# SciPy returns none
.b1_ipm_refs <- list(
  t_lp1_ipm = list("t_lp1", "", 0, 6, c(0x1.8p+0, 0x1p+1, 0x1p-1)),
  t_lp1_ipm_nopresolve = list("t_lp1", "nopresolve", 0, 6, c(0x1.8p+0, 0x1p+1, 0x1p-1)),
  t_lp1_ipm_nocrossover = list("t_lp1", "nocrossover", 0, 6, c(0x1.8000007420e08p+0, 0x1.ffffff8a8250cp+0, 0x1.fffffe2ca0eadp-2)),
  t_lp2_ipm = list("t_lp2", "", 0, 4, c(0x1.999999999999ap+0, 0x1.3333333333333p+0)),
  t_lp2_ipm_nocrossover = list("t_lp2", "nocrossover", 0, 4, c(0x1.9999999d225b1p+0, 0x1.33333334604f4p+0)),
  t_lpmax_ipm = list("t_lpmax", "", 0, 5, c(0x1.999999999999ap+0, 0x1.3333333333333p+0)),
  t_lpmax_ipm_nopresolve = list("t_lpmax", "nopresolve", 0, 5, c(0x1.999999999999ap+0, 0x1.3333333333333p+0)),
  t_lpmax_ipm_nocrossover = list("t_lpmax", "nocrossover", 0, 4, c(0x1.999999203f36cp+0, 0x1.33333366bf18fp+0)),
  t_lpinf_ipm = list("t_lpinf", "", 2, 0, NULL),
  t_lpinf_ipm_nopresolve = list("t_lpinf", "nopresolve", 2, 5, NULL),
  t_lpinf_ipm_nocrossover = list("t_lpinf", "nocrossover", 2, 0, NULL),
  t_lpunb_ipm = list("t_lpunb", "", 3, 1, NULL),
  t_lpunb_ipm_nopresolve = list("t_lpunb", "nopresolve", 3, 2, NULL),
  t_lpunb_ipm_nocrossover = list("t_lpunb", "nocrossover", 3, 1, NULL),
  t_lp3_ipm = list("t_lp3", "", 0, 0, c(0x1p+0, 0x1.8p+1, -0x1p+1, 0x0p+0)),
  t_lp3_ipm_nopresolve = list("t_lp3", "nopresolve", 0, 4, c(0x1p+0, 0x1.8p+1, -0x1p+1, 0x0p+0)),
  t_lp3_ipm_nocrossover = list("t_lp3", "nocrossover", 0, 0, c(0x1p+0, 0x1.8p+1, -0x1p+1, 0x0p+0)),
  lp_dense_ipm_nopresolve = list("lp_dense", "nopresolve", 0, 7, c(0x0p+0, 0x1.b20de3bae3bfbp+0, 0x0p+0, 0x1.a50de51db82ap-1, 0x1.4f8dd4d589216p+0, 0x0p+0, 0x0p+0, 0x1.ab3baa1231347p+0, 0x0p+0, 0x1.7c04c5457dd6p+0, 0x1.5a8abbaa9f563p+0, 0x1.0df15c16f1b56p+0)),
  lp_dense_ipm_nocrossover = list("lp_dense", "nocrossover", 0, 7, c(0x1.4182374bp-36, 0x1.b20de3baad749p+0, 0x1.a79ae28dp-40, 0x1.a50de51dc1498p-1, 0x1.4f8dd4d590fc5p+0, 0x1.4aa828ecp-40, 0x1.bbac5547p-40, 0x1.ab3baa120be9ep+0, 0x1.84f121a1p-36, 0x1.7c04c545990d3p+0, 0x1.5a8abbab14ef9p+0, 0x1.0df15c16f28e7p+0)),
  lp_klee_minty_ipm = list("lp_klee_minty", "", 0, 0, c(0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.86ap+11)),
  lp_klee_minty_ipm_nopresolve = list("lp_klee_minty", "nopresolve", 0, 6, c(0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.86ap+11)),
  lp_klee_minty_ipm_nocrossover = list("lp_klee_minty", "nocrossover", 0, 0, c(0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.86ap+11)),
  lp_free_fixed_ipm = list("lp_free_fixed", "", 0, 0, c(-0x0p+0, -0x1p+0, 0x0p+0, 0x0p+0)),
  lp_free_fixed_ipm_nopresolve = list("lp_free_fixed", "nopresolve", 0, 4, c(0x0p+0, -0x1p+0, 0x0p+0, 0x0p+0)),
  lp_free_fixed_ipm_nocrossover = list("lp_free_fixed", "nocrossover", 0, 0, c(0x0p+0, -0x1p+0, 0x0p+0, 0x0p+0)),
  lp_transport_ipm = list("lp_transport", "", 0, 6, c(0x0p+0, 0x1.4p+4, 0x0p+0, 0x0p+0, 0x1.4p+3, 0x0p+0, 0x1.ep+3, 0x1.4p+2, 0x0p+0, 0x1.4p+2, 0x0p+0, 0x1.4p+4)),
  lp_transport_ipm_nopresolve = list("lp_transport", "nopresolve", 0, 6, c(0x0p+0, 0x1.4p+4, 0x0p+0, 0x0p+0, 0x1.4p+3, 0x0p+0, 0x1.ep+3, 0x1.4p+2, 0x0p+0, 0x1.4p+2, 0x0p+0, 0x1.4p+4)),
  lp_transport_ipm_nocrossover = list("lp_transport", "nocrossover", 0, 6, c(0x1.00cb0a74p-38, 0x1.3ffffffff9d77p+4, 0x1.3ae9c2f5p-35, 0x1.a8cdbda8p-35, 0x1.3ffffffff84bfp+3, 0x1.9b270815p-37, 0x1.dffffffff204fp+3, 0x1.400000002822fp+2, 0x1.b3df0bf1p-35, 0x1.40000000153dfp+2, 0x1.2abb02c4p-34, 0x1.3ffffffff2949p+4)),
  lp_sparse_ipm = list("lp_sparse", "", 0, 8, c(0x1.8p+1, 0x1.8p+1, -0x1p+0, 0x1.9fc93124a93cp-4, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.8p+1, 0x1.8p+1, -0x1.70b25a57bf038p+2, -0x1.57dca04ccc8eap+1, -0x1.6cbd455c34049p+4, 0x1.3f835ac95bcb5p+1, 0x0p+0, 0x1.a9d466dea2781p+3, 0x0p+0, 0x1p+1, -0x1.382edcbb28311p+4, -0x1.6f10519255de3p+3)),
  lp_sparse_ipm_nopresolve = list("lp_sparse", "nopresolve", 0, 8, c(0x1.8p+1, 0x1.8p+1, -0x1p+0, 0x1.9fc93124a9443p-4, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.8p+1, 0x1.8p+1, -0x1.70b25a57bf036p+2, -0x1.57dca04ccc8eap+1, -0x1.6cbd455c3404bp+4, 0x1.3f835ac95bcb9p+1, 0x0p+0, 0x1.a9d466dea277ep+3, 0x0p+0, 0x1p+1, -0x1.382edcbb2831p+4, -0x1.6f10519255de8p+3)),
  lp_sparse_ipm_nocrossover = list("lp_sparse", "nocrossover", 0, 8, c(0x1.7ffffffffffdp+1, 0x1.7fffffffffeecp+1, -0x1.ffffffffff69fp-1, 0x1.9fc93124a550dp-4, -0x1.fffffffffff03p-1, -0x1.fffffffffd1ddp-1, -0x1.ffffffffffe5ep-1, -0x1.ffffffffffecap-1, 0x1.7fffffffffe43p+1, 0x1.7ffffffffff43p+1, -0x1.70b25a57bf441p+2, -0x1.57dca04ccd871p+1, -0x1.6cbd455c34114p+4, 0x1.3f835ac95b8d8p+1, 0x1.136e8b9cp-45, 0x1.a9d466dea2a52p+3, 0x1.51bde455p-46, 0x1.fffffffffff04p+0, -0x1.382edcbb282cap+4, -0x1.6f10519255c5dp+3)),
  lp_infeas_eq_ipm = list("lp_infeas_eq", "", 2, 0, NULL),
  lp_infeas_eq_ipm_nopresolve = list("lp_infeas_eq", "nopresolve", 2, 5, NULL),
  lp_infeas_eq_ipm_nocrossover = list("lp_infeas_eq", "nocrossover", 2, 0, NULL),
  lp_unb_free_ipm = list("lp_unb_free", "", 3, 0, NULL),
  lp_unb_free_ipm_nocrossover = list("lp_unb_free", "nocrossover", 3, 0, NULL),
  lp_dense_ipm_maxiter0 = list("lp_dense", "maxiter0", 1, 0, NULL),
  lp_dense_ipm_maxiter1 = list("lp_dense", "maxiter1", 1, 1, NULL),
  lp_dense_ipm_maxiter2 = list("lp_dense", "maxiter2", 1, 2, NULL),
  lp_dense_ipm_maxiter3 = list("lp_dense", "maxiter3", 1, 3, NULL),
  lp_dense_ipm_maxiter5 = list("lp_dense", "maxiter5", 1, 5, NULL)
)
# IPX without crossover: presolve sparsifies with the equations
# (lp_presolve_requires_basis_postsolve false) and postsolve recomputes the
# row activities; the reference slacks and marginals
.b1_ipm_ineqlin <- list(
  t_lp1_ipm_nocrossover = list(
    m = c(-0x1.00000000000b1p-1, -0x1.2cf67aee41000p-46, -0x1.ffffffffffea9p+0),
    r = c(0x1.d8d2500000000p-26, 0x1.fffffea4fa2d0p-1, 0x1.4bd3e00000000p-33)),
  lp_dense_ipm_nocrossover = list(
    m = c(-0x1.5d0a321000000p-38, -0x1.2b2794772b68ep-1, -0x1.b8f4fc4700000p-38, -0x1.1eda73f0cadbep-1,
          -0x1.f46259a8758f2p-2, -0x1.651754ad00000p-39, -0x1.fc00cce186773p-4, -0x1.aa96370800000p-37),
    r = c(0x1.9b0769597adedp+0, 0x1.0e20000000000p-38, 0x1.173623c050a12p+0, 0x1.8210000000000p-40,
          0x1.74c0000000000p-39, 0x1.1f15f86497b8ep+1, 0x1.1b4ad80000000p-33, 0x1.41e4d317d470cp-1)))

test_that("highs-ipm runs IPX, the simplex clean up and postsolve as Highs::run", {
  ex <- parse(test_path("test-optimization-scipy-highs-probe.R"))
  # evaluates only the .hs_probe_lps data assignment of the probe test file
  for (e in ex) if (is.call(e) && identical(e[[2]], as.name(".hs_probe_lps"))) eval(e)
  for (case in names(.b1_ipm_refs)) {
    ref <- .b1_ipm_refs[[case]]
    p <- .hs_probe_lps[[ref[[1]]]]
    kw <- list(A_ub = p$A_ub, b_ub = p$b_ub, A_eq = p$A_eq, b_eq = p$b_eq)
    opts <- list()
    if (ref[[2]] == "nopresolve") kw$presolve <- FALSE
    if (ref[[2]] == "nocrossover") kw$run_crossover <- "off"
    if (startsWith(ref[[2]], "maxiter")) opts$maxiters <- as.integer(sub("maxiter", "", ref[[2]]))
    r <- .scipy_run_linprog(list(method = "highs-ipm"), function(u) p$c, numeric(length(p$c)),
                            p$lb, p$ub, opts, kw)
    expect_identical(r$status, as.integer(ref[[3]]), info = case)
    expect_identical(as.integer(r$nit), as.integer(ref[[4]]), info = case)
    if (!is.null(ref[[5]]))
      expect_identical(sprintf("%a", r$x), sprintf("%a", ref[[5]]), info = case)
    if (!is.null(.b1_ipm_ineqlin[[case]])) {
      expect_identical(sprintf("%a", r$ineqlin$marginals), sprintf("%a", .b1_ipm_ineqlin[[case]]$m), info = case)
      expect_identical(sprintf("%a", r$ineqlin$residual), sprintf("%a", .b1_ipm_ineqlin[[case]]$r), info = case)
    }
  }
})

test_that("ported MILP solver matches the probe (needs every MIP group)", {
  needed <- c(".hmip_new_HighsDomain", ".hmip_new_HighsLpRelaxation", ".hmip_new_HighsCutPool",
              ".hmip_new_HighsPrimalHeuristics", ".hmip_new_HPresolve", ".hmip_new_HighsSymmetries")
  skip_if_not(all(vapply(needed, exists, TRUE)), "MIP groups not integrated")
  lp <- .hmip_HighsLp_from_dense(c(-1, -2), rbind(c(1, 1), c(1, -1)), c(-Inf, -Inf), c(3.5, 1.5),
                                 c(0, 0), c(10, 10), integrality = c(1L, 1L))
  r <- .hmip_callSolveMip(lp)
  expect_identical(r$model_status, "Optimal")
  expect_identical(r$objective, -6)
})
