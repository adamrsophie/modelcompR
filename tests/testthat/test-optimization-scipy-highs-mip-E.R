# Group E of the HiGHS MIP port: the MIP rules of HPresolve, the MIP
# postsolve and the HPresolve and HighsPostsolveStack classes. Expected
# values are the presolved models HiGHS 1.12.0 gives on Linux x86-64
# (tools/highs-probe.py raw_presolve, CI run highsE2u), compared bit for bit.

hx_same <- function(a, b) {
  length(a) == length(b) && all(sprintf("%a", as.numeric(a)) == sprintf("%a", as.numeric(b)))
}
expect_model <- function(lp, out) {
  expect_true(hx_same(lp$cost, out$cost))
  expect_true(hx_same(lp$col_lower, out$col_lower))
  expect_true(hx_same(lp$col_upper, out$col_upper))
  expect_true(hx_same(lp$row_lower, out$row_lower))
  expect_true(hx_same(lp$row_upper, out$row_upper))
  expect_true(hx_same(lp$offset, out$offset))
  expect_identical(as.integer(lp$index), out$index)
  expect_true(hx_same(lp$value, out$value))
}

test_that("integer bounds are rounded and shifted to zero (mpre_round)", {
  r <- .hs_mip_presolve(c(1, -1, 2), rbind(c(1, 1, 1), c(2, -1, 1)), c(-Inf, 1), c(5.5, Inf),
                        c(0.2, -1.5, 0), c(3.7, 4.2, 10), c(1L, 1L, 0L))
  expect_identical(r$status, "Reduced")
  expect_model(r$lp, list(cost = c(1, -1, 2), col_lower = c(0, 0, 0), col_upper = c(2, 5, 5.5),
                          row_lower = c(-Inf, -2), row_upper = c(5.5, Inf), offset = 2,
                          index = c(0L, 1L, 0L, 1L, 0L, 1L), value = c(1, 2, 1, -1, 1, 1)))
  types <- vapply(r$postsolve_stack$reductions, `[[`, "", "type")
  expect_identical(types, c("lineartransform", "lineartransform"))
  # the probe's optimal x = (2, 3, 0) maps to the presolved space and back
  red <- .hmip_HighsPostsolveStack_getReducedPrimalSolution(r$postsolve_stack, c(2, 3, 0))
  expect_identical(red, c(1, 4, 0))
  expect_identical(.hmip_HighsPostsolveStack_undoPrimal(r$postsolve_stack, list(), red), c(2, 3, 0))
})

test_that("strengthenInequalities without other reductions leaves the model not reduced", {
  r <- .hs_mip_presolve(c(-5, -3, -2, -1), rbind(c(10, 3, 2, 1), c(1, 1, 1, 1)), c(-Inf, -Inf), c(12, 3),
                        rep(0, 4), rep(1, 4), rep(1L, 4))
  expect_identical(r$status, "NotReduced")
  expect_true("Strengthened 1 coefficients" %in% r$pr$log)
  # 10 x1 + 3 x2 + 2 x3 + x4 <= 12 becomes 4 x1 + 3 x2 + x3 + x4 <= 5
  expect_identical(r$lp$row_upper, c(5, 3))
  expect_identical(r$lp$value, c(4, 1, 3, 1, 1, 1, 1, 1))
})

test_that("dominated columns, coefficient strengthening and scaleMIP (mpre_domcol)", {
  r <- .hs_mip_presolve(c(-2, -1, -3, -1), rbind(c(1, 1, 1, 1), c(2, 2, 1, 0), c(1, 1, 3, 2)), rep(-Inf, 3),
                        c(2.5, 3, 4), rep(0, 4), c(1, 1, 1, 4), c(1L, 1L, 1L, 0L))
  expect_identical(r$status, "Reduced")
  expect_identical(r$pr$log, c("3 rows, 4 cols, 11 nonzeros", "Strengthened 1 coefficients",
                               "Fixed 1 dominated columns and strengthened 0 bounds",
                               "2 rows, 3 cols, 6 nonzeros"))
  expect_model(r$lp, list(cost = c(-2, -3, -1), col_lower = c(0, 0, 0), col_upper = c(1, 1, 2),
                          row_lower = c(-Inf, -Inf), row_upper = c(2.5, 2), offset = 0,
                          index = c(0L, 1L, 0L, 1L, 0L, 1L), value = c(1, 0.5, 1, 1.5, 1, 1)))
  red <- .hmip_HighsPostsolveStack_getReducedPrimalSolution(r$postsolve_stack, c(1, 0, 1, 0))
  expect_identical(.hmip_HighsPostsolveStack_undoPrimal(r$postsolve_stack, list(), red), c(1, 0, 1, 0))
})

test_that("presolve statuses match HiGHS", {
  st <- function(...) .hs_mip_presolve(...)$status
  expect_identical(st(c(-1, -1), rbind(c(1, 1)), -Inf, 1.5, c(0, 0), c(Inf, Inf), c(1L, 1L)), "ReducedToEmpty")
  expect_identical(st(c(-1, -2, 0.5, 0.25), rbind(c(1, 2, 1, 0), c(3, 1, 0, -1)), c(9, 2), c(9, 2),
                      rep(0, 4), c(6, 6, Inf, Inf), c(1L, 1L, 0L, 0L)), "NotReduced")
})

hx_fixtures <- list(
  milp_genint = list(
    c = c(-3, -8, -10, -6, -5, -3, -8, -6, -9, -1),
    A = matrix(c(5, -1, 7, 7, -1, 5, 0, -3, 8, 4, -2, -3, 0, 5, 0, 7, 0, 4, 2, 2, 6, -2, 6, 3, -3, 1, 4, -1, 8, 6, -2, -1, -2, 3, 5, 0, 7, 1, -1, 8, 0, 8, 8, 6, 4, 6, -2, 0, 8, 5, 1, -3, 3, 6, 5, -1, 3, 5, 7, 2), 6, 10),
    lbc = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), ubc = c(16, 29, 18, 39, 16, 24),
    lb = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ub = c(8, 8, 8, 8, 8, 8, 8, 8, 8, 8),
    int = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    out = list(cost = c(-3, -8, -10, -6, -5, -3, -8, -6, -9), col_lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0), col_upper = c(5, 4, 8, 6, 4, 7, 4, 5, 5), row_lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), row_upper = c(16, 29, 18, 39, 16, 24),
               offset = c(0), index = c(0L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 3L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 0L, 1L, 2L, 3L, 5L, 0L, 1L, 2L, 3L, 4L, 0L, 1L, 2L, 3L, 4L, 5L), value = c(5, -1, 7, 7, -1, 5, -3, 8, 4, -2, -3, 5, 7, 4, 2, 2, 6, -2, 6, 3, -3, 1, 4, -1, 8, 6, -2, -1, -2, 3, 5, 7, 1, -1, 8, 8, 8, 6, 4, 6, -2, 8, 5, 1, -3, 3, 6))),
  milp_int8 = list(
    c = c(-6, -6, -1, -8, -7, -3, -8, -4),
    A = matrix(c(-1, 2, 5, 5, 3, -3, 5, 5, 0, 3, 4, 0, 2, 1, -1, 0, 5, 2, 3, 0, 0, 5, 0, -2, 2, -1, 4, 3, -1, 2, 1, 5, -3, -1, 0, 0, 2, 4, -3, 2, 2, 0, 3, 0, 0, 0, -3, 5), 6, 8),
    lbc = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), ubc = c(14, 12, 16, 12, 13, 10),
    lb = c(0, 0, 0, 0, 0, 0, 0, 0), ub = c(5, 5, 5, 5, 5, 5, 5, 5),
    int = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    out = list(cost = c(-6, -1, -8, -7, -3, -8, -4), col_lower = c(0, 0, 0, 0, 0, 0, 0), col_upper = c(3, 5, 3, 5, 3, 4, 5), row_lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), row_upper = c(14, 12, 16, 12, 13, 10),
               offset = c(0), index = c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 4L, 5L, 0L, 3L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L, 4L, 0L, 4L, 5L), value = c(-1, 2, 5, 5, 3, -3, 2, 1, -1, 5, 2, 3, 5, -2, 2, -1, 4, 3, -1, 2, 1, 5, -3, -1, 2, 4, -3, 2, 2, 3, -3, 5))),
  milp_knapsack = list(
    c = c(-4, -2, -1, -2, -10, -7, -8, -3),
    A = matrix(c(12, 2, 1, 1, 4, 7, 9, 3), 1, 8),
    lbc = c(-Inf), ubc = c(15),
    lb = c(0, 0, 0, 0, 0, 0, 0, 0), ub = c(1, 1, 1, 1, 1, 1, 1, 1),
    int = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    out = list(cost = c(-1, -2, -10, -8), col_lower = c(1, 0, 0, 0), col_upper = c(13, 1, 1, 1), row_lower = c(-Inf), row_upper = c(15),
               offset = c(0), index = c(0L, 0L, 0L, 0L), value = c(1, 1, 4, 9))),
  milp_cover30 = list(
    c = c(0x1.075810fc58d95p+1, 0x1.28b02f5d3216cp+1, 0x1.123060e03dd86p+1, 0x1.aef3c6e563e05p+0, 0x1.041685f24d234p+1, 0x1.3f9ea1fb1f7afp+1, 0x1.19a77d12bbdc6p+0, 0x1.8bcb50976db4cp+0, 0x1.6619e293b37b9p+0, 0x1.a369ce4e27c16p+0, 0x1.c2112a8ab35cp+0, 0x1.7b0117a2d755cp+0, 0x1.655c610e22526p+1, 0x1.6603f1e614451p+1, 0x1.d1e2cdcbd4bfdp+0, 0x1.7b319dbbc324ep+0, 0x1.2fca2355d2b3ap+1, 0x1.b58272e460e78p+0, 0x1.b57b53f02d2a5p+0, 0x1.1aa4c9f3bc74ap+1),
    A = matrix(c(0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 30, 20),
    lbc = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), ubc = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
    lb = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ub = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    int = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    out = list(cost = c(0x1.075810fc58d95p+1, 0x1.28b02f5d3216cp+1, 0x1.aef3c6e563e05p+0, 0x1.041685f24d234p+1, 0x1.3f9ea1fb1f7afp+1, 0x1.19a77d12bbdc6p+0, 0x1.8bcb50976db4cp+0, 0x1.6619e293b37b9p+0, 0x1.a369ce4e27c16p+0, 0x1.c2112a8ab35cp+0, 0x1.7b0117a2d755cp+0, 0x1.655c610e22526p+1, 0x1.6603f1e614451p+1, 0x1.d1e2cdcbd4bfdp+0, 0x1.7b319dbbc324ep+0, 0x1.2fca2355d2b3ap+1, 0x1.b58272e460e78p+0, 0x1.b57b53f02d2a5p+0, 0x1.1aa4c9f3bc74ap+1), col_lower = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), col_upper = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), row_lower = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), row_upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
               offset = c(0), index = c(3L, 4L, 15L, 16L, 21L, 28L, 29L, 0L, 1L, 4L, 6L, 12L, 14L, 17L, 28L, 0L, 2L, 3L, 4L, 5L, 7L, 8L, 15L, 17L, 19L, 22L, 5L, 17L, 21L, 23L, 24L, 14L, 26L, 27L, 1L, 5L, 14L, 25L, 2L, 11L, 15L, 20L, 28L, 0L, 1L, 6L, 11L, 15L, 18L, 10L, 14L, 16L, 18L, 21L, 22L, 24L, 1L, 3L, 7L, 15L, 24L, 26L, 29L, 1L, 4L, 10L, 13L, 14L, 16L, 21L, 27L, 0L, 3L, 8L, 9L, 27L, 1L, 13L, 19L, 22L, 29L, 0L, 5L, 6L, 10L, 12L, 15L, 21L, 22L, 23L, 8L, 9L, 10L, 13L, 14L, 20L, 22L, 24L, 0L, 7L, 15L, 18L, 19L, 21L, 26L, 2L, 24L, 27L, 3L, 14L, 16L, 21L, 25L, 9L, 17L, 21L), value = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
)

test_that("presolved models match HiGHS bit for bit (probing proper needs a mipsolver)", {
  for (nm in names(hx_fixtures)) {
    f <- hx_fixtures[[nm]]
    r <- .hs_mip_presolve(f$c, f$A, f$lbc, f$ubc, f$lb, f$ub, f$int)
    expect_identical(r$status, "Reduced", info = nm)
    expect_model(r$lp, f$out)
  }
})

test_that("DuplicateColumn::okMerge accepts and refuses merges as HiGHS does", {
  expect_true(.hs_dupcol_ok_merge(1, 0, 1, 0, 1, TRUE, TRUE, 1e-6))
  expect_false(.hs_dupcol_ok_merge(0.5, 0, 1, 0, 1, TRUE, TRUE, 1e-6))
  expect_false(.hs_dupcol_ok_merge(3, 0, 1, 0, 1, TRUE, TRUE, 1e-6))
  expect_false(.hs_dupcol_ok_merge(1, 0, 1, 2, 2, TRUE, FALSE, 1e-6))
  expect_false(.hs_dupcol_ok_merge(3, 0, 1, 0, 5, FALSE, TRUE, 1e-6))
  expect_true(.hs_dupcol_ok_merge(0.7, 0, 1, 0, 5, FALSE, FALSE, 1e-6))
})

test_that("postsolve stack index maps follow HighsPostsolveStack", {
  s <- .hmip_new_HighsPostsolveStack()
  .hmip_HighsPostsolveStack_initializeIndexMaps(s, 3L, 4L)
  expect_identical(.hmip_HighsPostsolveStack_getOrigNumCol(s), 4L)
  expect_identical(.hmip_HighsPostsolveStack_getOrigNumRow(s), 3L)
  .hmip_HighsPostsolveStack_compressIndexMaps(s, c(0L, -1L, 1L), c(-1L, 1L, 0L, -1L))
  expect_identical(.hmip_HighsPostsolveStack_getOrigColIndex(s, 0L), 2L)
  expect_identical(.hmip_HighsPostsolveStack_getOrigColIndex(s, 1L), 1L)
  expect_identical(.hmip_HighsPostsolveStack_getOrigRowIndex(s, 1L), 2L)
  .hmip_HighsPostsolveStack_appendCutsToModel(s, 2L)
  expect_identical(s$origRowIndex, c(0L, 2L, 3L, 4L))
  expect_identical(.hmip_HighsPostsolveStack_getOrigNumRow(s), 5L)
  .hmip_HighsPostsolveStack_removeCutsFromModel(s, 2L)
  expect_identical(s$origRowIndex, c(0L, 2L))
  expect_identical(.hmip_HighsPostsolveStack_getOrigNumRow(s), 3L)
})

test_that("the dominated column signature bit is the top five hash bits", {
  for (r in 0:40) expect_identical(.hs_row_sig_bit(r) - 1, .hs_u64_top(.hmip_hash(r), 5L))
})

test_that("liftingForProbing lifts the best clique and drops rows implied by both values", {
  # group A's clique table and domain are replaced by stubs: one clique of all candidates, nothing fixed
  stubs <- list(.hmip_HighsCliqueTable_computeMaximalCliques = function(ct, vars, feastol) list(vars),
                .hmip_HighsDomain_isFixed = function(domain, col) FALSE)
  ns <- asNamespace("modelcompR")
  tenv <- environment()
  missing <- character()
  for (nm in names(stubs)) {
    if (exists(nm, envir = ns, inherits = FALSE)) {
      do.call(local_mocked_bindings, c(stubs[nm], list(.package = "modelcompR", .env = tenv)))
    } else {
      assign(nm, stubs[[nm]], envir = globalenv())
      missing <- c(missing, nm)
    }
  }
  withr::defer(rm(list = missing, envir = globalenv()))
  r <- .hs_mip_presolve(c(-2, -1, -3, -1), rbind(c(1, 1, 1, 1), c(2, 2, 1, 0), c(1, 1, 3, 2)), rep(-Inf, 3),
                        c(2.5, 3, 4), rep(0, 4), c(1, 1, 1, 4), c(1L, 1L, 1L, 0L))
  pr <- r$pr
  pr$mipsolver <- list(mipdata_ = list(cliquetable = NULL, domain = NULL),
                       options_mip_ = list(mip_lifting_for_probing = 1L))
  row1 <- pr$ridx[1L]
  row2 <- pr$ridx[2L]
  # (col, val) repeats: insert_or_get keeps the first coefficient
  pr$liftingOpportunities[[as.character(row1)]] <- list(col = c(0L, 1L, 0L), val = c(1L, 0L, 1L), coef = c(0.5, 0.25, 9))
  pr$liftingOpportunities[[as.character(row2)]] <- list(col = c(2L, 2L), val = c(1L, 0L), coef = c(1, 1))
  .hs_lifting_for_probing(pr)
  expect_identical(tail(pr$log, 1L),
                   "Lifting for probing removed 1 and modified 1 row(s), added 0 new and modified 2 existing nonzero(s)")
  expect_true(pr$rowDeleted[row2])
  expect_identical(pr$ru[row1], 2.75)
  vals <- vapply(pr$cidx, function(j) pr$Avalue[.hs_find_nonzero(pr, row1, j)], 0)
  expect_identical(vals, c(1.5, 1.25, 1))
})
