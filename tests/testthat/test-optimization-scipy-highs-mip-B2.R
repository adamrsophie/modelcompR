# HiGHS 1.12.0 MIP port, group B2: HighsPseudocost, HighsNodeQueue,
# HighsHashTable, the Highs lpsolver adapter, HighsLpRelaxation and the
# objective bound bailout of the dual simplex method.

ns <- asNamespace("modelcompR")
hfn <- function(name) get(name, envir = ns)

# milp_cover30 of tools/highs-probe.py (30 covering rows sum(A x) >= 1 over
# 20 binaries); the root LP takes 16 simplex iterations in the HiGHS log
cover30_cost <- c(0x1.075810fc58d95p+1, 0x1.28b02f5d3216cp+1, 0x1.123060e03dd86p+1, 0x1.aef3c6e563e05p+0,
                  0x1.041685f24d234p+1, 0x1.3f9ea1fb1f7afp+1, 0x1.19a77d12bbdc6p+0, 0x1.8bcb50976db4cp+0,
                  0x1.6619e293b37b9p+0, 0x1.a369ce4e27c16p+0, 0x1.c2112a8ab35cp+0, 0x1.7b0117a2d755cp+0,
                  0x1.655c610e22526p+1, 0x1.6603f1e614451p+1, 0x1.d1e2cdcbd4bfdp+0, 0x1.7b319dbbc324ep+0,
                  0x1.2fca2355d2b3ap+1, 0x1.b58272e460e78p+0, 0x1.b57b53f02d2a5p+0, 0x1.1aa4c9f3bc74ap+1)
cover30_bits <- paste0(
  "0101000010001010100001000010101101000000000100010000000001001001000000101000001011010000000100000000",
  "0001101000000010000001000000100000100000000100000010000010000001000000001001000000000000000010010001",
  "0000000001010011000000000001100000000000010000000000001000000000000000010101000001000110010100010010",
  "1001000110100010100010000000010100000010011110000000000000010000000011000000100000010000000001001000",
  "0000000100000001000010001000010100101011000100000100011100000010100000000010000000001000011000010100",
  "0000001000000000001000000100001000001000000001000001100001001100000100000000000010000000001001000000")
cover30_A <- matrix(as.numeric(strsplit(cover30_bits, "")[[1L]]), 30L, 20L, byrow = TRUE)
cover30_obj <- 11.61475723

cover30_highs <- function() {
  h <- hfn(".hmip_new_Highs")()
  hfn(".hmip_Highs_setOptionValue")(h, "primal_feasibility_tolerance", 1e-6)
  hfn(".hmip_Highs_setOptionValue")(h, "dual_feasibility_tolerance", 1e-7)
  hfn(".hmip_Highs_passModel")(h, list(num_col_ = 20L, num_row_ = 30L, col_cost_ = cover30_cost,
                                        col_lower_ = numeric(20), col_upper_ = rep(1, 20),
                                        row_lower_ = rep(1, 30), row_upper_ = rep(Inf, 30), A = cover30_A))
  h
}

# stand-ins for the members of other groups this file needs, only where
# those groups have not landed yet: the namespace is locked, but its lookup
# chain ends in the global environment, where the stand-ins go for the test
# names in `always` replace the real function too (they record what the
# relaxation passes on instead of running the whole solver)
with_stubs <- function(stubs, code, always = character()) {
  added <- character()
  for (nm in intersect(always, names(stubs))) if (exists(nm, envir = ns))
    do.call(testthat::local_mocked_bindings, c(stubs[nm], list(.package = "modelcompR")))
  for (nm in names(stubs)) if (!exists(nm, envir = ns)) {
    assign(nm, stubs[[nm]], envir = globalenv())
    added <- c(added, nm)
  }
  on.exit(if (length(added)) rm(list = added, envir = globalenv()), add = TRUE)
  force(code)
}

test_that("HighsPseudocost observations and pseudocosts follow the C++ updates", {
  pc <- hfn(".hmip_new_HighsPseudocost")()
  ms <- list(model_ = list(num_col_ = 2L), options_mip_ = list(mip_pscost_minreliable = 8L), pscostinit = NULL)
  hfn(".hmip_HighsPseudocost_init")(pc, ms)
  hfn(".hmip_HighsPseudocost_addObservation")(pc, 0L, 0.5, 1)
  expect_identical(pc$pseudocostup[1L], 2)
  expect_identical(pc$cost_total, 2)
  hfn(".hmip_HighsPseudocost_addObservation")(pc, 0L, -0.25, 1)
  expect_identical(pc$pseudocostdown[1L], 4)
  expect_identical(pc$cost_total, 3)
  # two argument overload: (ceil(frac) - frac) times the column cost, or the
  # average for a column without samples
  expect_identical(hfn(".hmip_HighsPseudocost_getPseudocostUp")(pc, 0L, 0.25), 0.75 * 2)
  expect_identical(hfn(".hmip_HighsPseudocost_getPseudocostUp")(pc, 1L, 0.25), 0.75 * 3)
  # offset overload with nsamples 1 < minreliable 8: weight 0.9 + 0.1 / 8
  w <- 0.9 + 0.1 * 1 / 8
  expect_identical(hfn(".hmip_HighsPseudocost_getPseudocostDown")(pc, 0L, 0.25, 0.5),
                   0.25 * (0.5 + (w * 4 + (1 - w) * 3)))
  expect_false(hfn(".hmip_HighsPseudocost_isReliable")(pc, 0L))
  hfn(".hmip_HighsPseudocost_setMinReliable")(pc, 1L)
  expect_true(hfn(".hmip_HighsPseudocost_isReliable")(pc, 0L))
})

test_that("HighsNodeQueue orders nodes, reuses free slots and bounds", {
  nq <- hfn(".hmip_new_HighsNodeQueue")()
  hfn(".hmip_HighsNodeQueue_setNumCol")(nq, 2L)
  dc <- function(v, c, t) list(boundval = v, column = c, boundtype = t)
  emplace <- hfn(".hmip_HighsNodeQueue_emplaceNode")
  expect_identical(emplace(nq, dc(1, 0L, 0L), 0L, 5, 7, 2L), 0)
  emplace(nq, dc(0, 0L, 1L), 0L, 3, 9, 2L)
  emplace(nq, dc(c(1, 1), c(0L, 1L), c(0L, 0L)), c(0L, 1L), 4, 4, 3L)
  expect_identical(hfn(".hmip_HighsNodeQueue_numNodesUp")(nq, 0L), 2L)
  expect_identical(hfn(".hmip_HighsNodeQueue_numNodesDown")(nq, 0L), 1L)
  expect_identical(hfn(".hmip_HighsNodeQueue_getBestLowerBound")(nq), 3)
  # best estimate node: 0.5 lb + 0.5 estimate is 6, 6 and 4
  best <- hfn(".hmip_HighsNodeQueue_popBestNode")(nq)
  expect_identical(best$lower_bound, 4)
  # the freed slot 2 is reused
  emplace(nq, dc(0, 1L, 1L), 0L, 10, 10, 4L)
  expect_length(nq$nodes, 3L)
  expect_identical(nq$nodes[[3L]]$lower_bound, 10)
  # bounding at 6 prunes the node of lower bound 10 at depth 4
  expect_identical(hfn(".hmip_HighsNodeQueue_performBounding")(nq, 6), 2^(1 - 4))
  expect_identical(hfn(".hmip_HighsNodeQueue_numNodes")(nq), 2L)
  b <- hfn(".hmip_HighsNodeQueue_popBestBoundNode")(nq)
  expect_identical(b$lower_bound, 3)
})

test_that("HighsHashTable keeps the HiGHS slot order", {
  skip_if_not(exists(".hmip_hash", envir = ns), "group B1 hash helpers not landed")
  ht <- hfn(".hmip_new_HighsHashTable")(0)
  for (k in c(5, 17, 3, 200, 64)) {
    s <- hfn(".hmip_HighsHashTable_ref")(ht, k, 0)
    hfn(".hmip_HighsHashTable_set")(ht, s, k * 2)
  }
  expect_identical(hfn(".hmip_HighsHashTable_size")(ht), 5)
  expect_setequal(hfn(".hmip_HighsHashTable_keys")(ht), c(5, 17, 3, 200, 64))
  expect_identical(hfn(".hmip_HighsHashTable_find")(ht, 200), 400)
  expect_null(hfn(".hmip_HighsHashTable_find")(ht, 4))
})

test_that("the Highs adapter solves the cover30 root LP and warm starts", {
  h <- cover30_highs()
  expect_identical(hfn(".hmip_Highs_run")(h), "Ok")
  expect_identical(h$model_status, "Optimal")
  expect_identical(h$info$simplex_iteration_count, 16L)
  expect_equal(h$info$objective_function_value, cover30_obj, tolerance = 1e-9)
  expect_true(h$ek_valid)
  hfn(".hmip_Highs_setOptionValue")(h, "presolve", "off")
  # a bound change keeps the HEkk; the resolve takes a few dual iterations
  hfn(".hmip_Highs_changeColsBounds")(h, 0L, 0, 0)
  hfn(".hmip_Highs_run")(h)
  expect_identical(h$model_status, "Optimal")
  expect_lt(h$info$simplex_iteration_count, 16L)
  expect_gt(h$info$objective_function_value, cover30_obj)
  obj_down <- h$info$objective_function_value
  # an added redundant row is basic in the warm start: no iterations
  hfn(".hmip_Highs_addRows")(h, 1L, 1, Inf, 20L, 0L, 0:19, rep(1, 20))
  expect_null(h$ek)
  hfn(".hmip_Highs_run")(h)
  expect_identical(h$info$simplex_iteration_count, 0L)
  expect_equal(h$info$objective_function_value, obj_down, tolerance = 1e-12)
  expect_identical(hfn(".hmip_Highs_getLp")(h)$num_row_, 31L)
  r <- hfn(".hmip_Highs_getBasisInverseRowSparse")(h, 0L)
  expect_true(r$ok)
  expect_identical(r$count, length(r$index))
  mask <- hfn(".hmip_Highs_deleteRowsMask")(h, c(integer(30), 1L))
  expect_identical(mask, c(0:29, -1L))
  expect_identical(h$lp$num_row_, 30L)
})

test_that("changeColsCost sets masked costs, keeps the HEkk and invalidates", {
  h <- cover30_highs()
  hfn(".hmip_Highs_run")(h)
  hfn(".hmip_Highs_setOptionValue")(h, "presolve", "off")
  lp0 <- hfn(".hmip_Highs_getLp")(h)
  newcost <- lp0$col_cost_ * 2
  mask <- c(1L, integer(29))
  hfn(".hmip_Highs_changeColsCost")(h, mask, newcost)
  expect_identical(h$model_status, "Notset")
  expect_false(h$solution$value_valid)
  expect_true(h$ek_valid)
  lp1 <- hfn(".hmip_Highs_getLp")(h)
  expect_identical(lp1$col_cost_, c(newcost[1L], lp0$col_cost_[-1L]))
  hfn(".hmip_Highs_run")(h)
  warm <- h$info$objective_function_value
  # a cold solve of the same costs reaches the same optimum
  h2 <- cover30_highs()
  hfn(".hmip_Highs_changeColsCost")(h2, mask, newcost)
  hfn(".hmip_Highs_run")(h2)
  expect_identical(h$model_status, "Optimal")
  expect_equal(warm, h2$info$objective_function_value, tolerance = 1e-9)
  # getLp follows added and deleted rows (a_matrix_ and row bounds)
  hfn(".hmip_Highs_addRows")(h, 1L, 1, 5, 2L, 0L, c(2L, 7L), c(1, 3))
  lp2 <- hfn(".hmip_Highs_getLp")(h)
  expect_identical(lp2$num_row_, 31L)
  expect_identical(lp2$row_lower_[31L], 1)
  expect_identical(lp2$row_upper_[31L], 5)
  in_row <- which(lp2$a_matrix_$index_ == 30L)
  expect_identical(lp2$a_matrix_$value_[in_row], c(1, 3))
  hfn(".hmip_Highs_deleteRowsMask")(h, c(integer(30), 1L))
  lp3 <- hfn(".hmip_Highs_getLp")(h)
  expect_identical(lp3$num_row_, 30L)
  expect_identical(lp3$a_matrix_, lp0$a_matrix_)
})

test_that("the objective bound stops the dual simplex method", {
  h <- cover30_highs()
  hfn(".hmip_Highs_run")(h)
  hfn(".hmip_Highs_setOptionValue")(h, "presolve", "off")
  hfn(".hmip_Highs_setOptionValue")(h, "objective_bound", cover30_obj + 0.01)
  hfn(".hmip_Highs_changeColsBounds")(h, c(0L, 1L, 3L), c(0, 0, 0), c(0, 0, 0))
  hfn(".hmip_Highs_run")(h)
  expect_identical(h$model_status, "ObjectiveBound")
  expect_true(h$solution$dual_valid)
})

cover30_mipsolver <- function() {
  nz <- which(cover30_A != 0, arr.ind = TRUE)
  nz <- nz[order(nz[, 2L], nz[, 1L]), , drop = FALSE]
  am <- list(format_ = "colwise", start_ = c(0L, cumsum(tabulate(nz[, 2L], 20L))), index_ = nz[, 1L] - 1L,
             value_ = cover30_A[nz])
  dom <- new.env(parent = emptyenv())
  dom$col_lower_ <- numeric(20); dom$col_upper_ <- rep(1, 20); dom$changed <- integer()
  md <- new.env(parent = emptyenv())
  md$domain <- dom; md$integral_cols <- 0:19; md$feastol <- 1e-6; md$epsilon <- 1e-9
  md$uplocks <- numeric(20); md$downlocks <- colSums(cover30_A)
  md$incumbents <- list()
  md$cliquetable <- hfn(".hmip_new_HighsCliqueTable")(20L)
  ms <- new.env(parent = emptyenv())
  ms$model_ <- list(num_col_ = 20L, num_row_ = 30L, col_cost_ = cover30_cost, col_lower_ = numeric(20),
                    col_upper_ = rep(1, 20), row_lower_ = rep(1, 30), row_upper_ = rep(Inf, 30),
                    integrality_ = rep(1L, 20), a_matrix_ = am)
  ms$options_mip_ <- list(random_seed = 0L, mip_feasibility_tolerance = 1e-6, small_matrix_value = 1e-9,
                          mip_lp_age_limit = 10L)
  ms$mipdata_ <- md
  ms$submip <- FALSE
  ms
}

lprel_stubs <- list(
  .hmip_hash = function(x, hi = 0) hfn(".hs_pair_hash")(0L, x, hi),
  .hmip_HighsCliqueTable_getSubstitution = function(ct, col) NULL,
  .hmip_HighsCliqueTable_getSubstitutions = function(ct) list(),
  .hmip_HighsMipSolverData_addIncumbent = function(md, sol, obj, src, ...) {
    md$incumbents[[length(md$incumbents) + 1L]] <- list(sol = sol, obj = obj, src = src)
    TRUE
  },
  .hmip_HighsMipSolverData_checkSolution = function(md, sol) FALSE,
  .hmip_HighsDomain_getChangedCols = function(dom) dom$changed,
  .hmip_HighsDomain_clearChangedCols = function(dom, start = 0L) dom$changed <- dom$changed[seq_len(start)],
  .hmip_HighsDomain_removeContinuousChangedCols = function(dom) invisible(),
  .hmip_HighsCutSet_numCuts = function(cs) length(cs$cutindices),
  .hmip_HighsCutSet_clear = function(cs) cs$cutindices <- integer(),
  .hmip_HighsCutPool_lpCutRemoved = function(cp, cut) invisible())

lprel_always <- c(".hmip_HighsMipSolverData_addIncumbent", ".hmip_HighsMipSolverData_checkSolution",
                  ".hmip_HighsCutSet_numCuts", ".hmip_HighsCutSet_clear", ".hmip_HighsCutPool_lpCutRemoved")

test_that("HighsLpRelaxation solves the cover30 root and rounds it", {
  with_stubs(lprel_stubs, always = lprel_always, code = {
    ms <- cover30_mipsolver()
    lp <- hfn(".hmip_new_HighsLpRelaxation")(ms)
    hfn(".hmip_HighsLpRelaxation_setAdjustSymmetricBranchingCol")(lp, FALSE)
    hfn(".hmip_HighsLpRelaxation_loadModel")(lp)
    expect_identical(lp$lpsolver$options$dual_feasibility_tolerance, 1e-6 * 0.1)
    st <- hfn(".hmip_HighsLpRelaxation_resolveLp")(lp)
    expect_identical(st, "Optimal")
    expect_identical(hfn(".hmip_HighsLpRelaxation_getNumLpIterations")(lp), 16)
    expect_equal(lp$objective, cover30_obj, tolerance = 1e-9)
    fi <- hfn(".hmip_HighsLpRelaxation_getFractionalIntegers")(lp)
    expect_setequal(fi$first, c(0L, 1L, 3L, 4L, 7L, 8L, 10L, 11L, 12L, 14L, 15L, 16L))
    x <- lp$lpsolver$solution$col_value
    expect_equal(fi$second, x[fi$first + 1L])
    # every fractional column has no up locks, so the rounded point (all
    # fractional columns up) goes to addIncumbent with source SolveLp
    inc <- ms$mipdata_$incumbents
    expect_length(inc, 1L)
    expect_identical(inc[[1L]]$src, 9L)
    expect_identical(inc[[1L]]$sol[fi$first + 1L], rep(1, length(fi$first)))
    expect_false(hfn(".hmip_HighsLpRelaxation_integerFeasible")(lp))
  })
})

test_that("HighsLpRelaxation adds and removes cuts and plays in the playground", {
  with_stubs(lprel_stubs, always = lprel_always, code = {
    ms <- cover30_mipsolver()
    lp <- hfn(".hmip_new_HighsLpRelaxation")(ms)
    hfn(".hmip_HighsLpRelaxation_setAdjustSymmetricBranchingCol")(lp, FALSE)
    hfn(".hmip_HighsLpRelaxation_loadModel")(lp)
    hfn(".hmip_HighsLpRelaxation_resolveLp")(lp)
    obj0 <- lp$objective
    cs <- new.env(parent = emptyenv())
    cs$cutindices <- 0L; cs$lower_ <- -Inf; cs$upper_ <- 100
    cs$ARstart_ <- c(0L, 20L); cs$ARindex_ <- 0:19; cs$ARvalue_ <- rep(1, 20)
    hfn(".hmip_HighsLpRelaxation_addCuts")(lp, cs)
    expect_identical(hfn(".hmip_HighsLpRelaxation_numRows")(lp), 31L)
    expect_identical(lp$lprows$origin[31L], "cut")
    cp <- hfn(".hmip_HighsLpRelaxation_getCutPool")(lp)
    expect_identical(cp$num_cut, 1L)
    expect_identical(cp$cut_matrix$index_, 0:19)
    hfn(".hmip_HighsLpRelaxation_resolveLp")(lp)
    expect_identical(lp$objective, obj0)
    # the cut is basic, so removeObsoleteRows drops it and resolves
    hfn(".hmip_HighsLpRelaxation_removeObsoleteRows")(lp)
    expect_identical(hfn(".hmip_HighsLpRelaxation_numRows")(lp), 30L)
    expect_length(lp$lprows$index, 30L)
    hfn(".hmip_HighsLpRelaxation_resolveLp")(lp)
    expect_identical(lp$objective, obj0)
    # playground: a probing bound change, then the destructor restores the
    # iterate and resolves
    dom <- new.env(parent = emptyenv())
    dom$col_lower_ <- numeric(20); dom$col_upper_ <- rep(1, 20); dom$col_upper_[1L] <- 0
    dom$changed <- 0L
    pg <- hfn(".hmip_HighsLpRelaxation_playground")(lp)
    st <- hfn(".hmip_Playground_solveLp")(pg, dom)
    expect_identical(st, "Optimal")
    expect_gt(lp$lpsolver$info$objective_function_value, obj0)
    dom$col_upper_[1L] <- 1; dom$changed <- 0L
    hfn(".hmip_HighsLpRelaxation_flushDomain")(lp, dom)
    hfn(".hmip_Playground_destroy")(pg)
    expect_equal(lp$lpsolver$info$objective_function_value, obj0, tolerance = 1e-12)
  })
})

test_that("HighsSearch flips a branching decision to the other child as backtrack does", {
  nd <- hfn(".hmip_b2_nodedata")(1, 2)
  nd$branching_point <- 2.4
  nd$branchingdecision <- list(boundval = 3, column = 4L, boundtype = 0L)
  expect_false(hfn(".hmip_b2_flip")(nd))
  expect_identical(nd$branchingdecision, list(boundval = 2, column = 4L, boundtype = 1L))
  # a fallback branch has branching_point equal to the bound value
  nd$branching_point <- 2
  expect_true(hfn(".hmip_b2_flip")(nd))
  expect_identical(nd$branchingdecision$boundval, 3)
  expect_identical(nd$branchingdecision$boundtype, 0L)
  nd$branching_point <- 0.3
  expect_identical(hfn(".hmip_b2_choose")(nd, TRUE, 5, 6), 6)
  expect_identical(nd$branchingdecision$boundval, 1)
  expect_identical(nd$other_child_lb, 5)
})
