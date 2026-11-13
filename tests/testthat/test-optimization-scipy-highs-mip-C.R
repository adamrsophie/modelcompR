# Group C of the HiGHS 1.12.0 MIP port: cut pool, cut generation and the
# separators (R/optimization-scipy-highs-mip-cuts.R). Hand computed values.

hmipC <- function(name) get(name, envir = asNamespace("modelcompR"))

test_that("std::frexp matches the C library", {
  frexp <- hmipC(".hmip_frexp")
  expect_identical(frexp(8), c(0.5, 4))
  expect_identical(frexp(-3), c(-0.75, 2))
  expect_identical(frexp(0.1), c(0.8, -3))
  expect_identical(frexp(1), c(0.5, 1))
  expect_identical(frexp(0), c(0, 0))
  for (k in c(-1021, -3, 0, 1, 52, 1023)) {
    expect_identical(frexp(2^k), c(0.5, k + 1))
    expect_identical(frexp(2^k * (1 + 2^-52)), c(0.5 + 2^-53, k + 1))
    expect_identical(frexp(-2^k * (1 - 2^-53)), c(-(1 - 2^-53), k))
  }
  # subnormals and the extremes
  expect_identical(frexp(2^-1074), c(0.5, -1073))
  expect_identical(frexp(3 * 2^-1074), c(0.75, -1072))
  expect_identical(frexp(2^-1022 - 2^-1074), c(1 - 2^-52, -1022))
  expect_identical(frexp(.Machine$double.xmax), c(1 - 2^-53, 1024))
  expect_identical(frexp(Inf), c(Inf, 0))
  expect_true(is.nan(frexp(NaN)[1L]))
})

test_that("the path delta rounds log2 as a correctly rounded libm does", {
  f <- hmipC(".hmip_exp2_ceil_log2")
  expect_identical(f(4), 4)
  expect_identical(f(4 - 2^-50), 4)
  # log2(4 (1 + 2^-52)) = 2 + 1.44 2^-52 rounds above 2, log2 of 16 (1 + 2^-52)
  # = 4 + 1.44 2^-52 rounds to 4
  expect_identical(f(4 * (1 + 2^-52)), 8)
  expect_identical(f(16 * (1 + 2^-52)), 16)
  expect_identical(f(16 * (1 + 2^-50)), 32)
  expect_identical(f(1.5), 2)
  expect_identical(f(5), 8)
})

test_that("HighsDynamicRowMatrix links columns by sign and reuses freed space", {
  m <- hmipC(".hmip_new_HighsDynamicRowMatrix")(3L)
  addRow <- hmipC(".hmip_HighsDynamicRowMatrix_addRow")
  expect_identical(addRow(m, c(0L, 2L), c(1, -2), 2L), 0L)
  expect_identical(addRow(m, c(1L, 2L), c(3, 4), 2L), 1L)
  pos <- integer(); neg <- integer()
  hmipC(".hmip_HighsDynamicRowMatrix_forEachPositiveColumnEntry")(m, 2L, function(r, v) {
    pos <<- c(pos, r); TRUE
  })
  hmipC(".hmip_HighsDynamicRowMatrix_forEachNegativeColumnEntry")(m, 2L, function(r, v) {
    neg <<- c(neg, r); TRUE
  })
  expect_identical(pos, 1L)
  expect_identical(neg, 0L)
  hmipC(".hmip_HighsDynamicRowMatrix_removeRow")(m, 0L)
  expect_identical(hmipC(".hmip_HighsDynamicRowMatrix_getNumDelRows")(m), 1L)
  # the deleted row index and its two nonzero slots come back
  expect_identical(addRow(m, c(0L, 1L), c(5, 6), 2L), 0L)
  expect_identical(hmipC(".hmip_HighsDynamicRowMatrix_getRowStart")(m, 0L), 0L)
  expect_identical(hmipC(".hmip_HighsDynamicRowMatrix_nonzeroCapacity")(m), 4L)
})

test_that("HighsCutPool adds, rejects duplicates and separates violated cuts", {
  skip_if_not(exists(".hmip_numNonzero", envir = asNamespace("modelcompR")))
  mip <- new.env()
  mip$model_ <- list(num_col_ = 2L, num_row_ = 1L,
                     a_matrix_ = list(start_ = c(0L, 1L, 2L), index_ = c(0L, 0L), value_ = c(1, 1)))
  mip$mipdata_ <- new.env()
  cp <- hmipC(".hmip_new_HighsCutPool")(2L, 10L, 100L)
  addCut <- hmipC(".hmip_HighsCutPool_addCut")
  r <- addCut(cp, mip, c(1L, 0L), c(1, 1), 2L, 1, TRUE, FALSE, FALSE)
  expect_identical(r$index, 0L)
  expect_identical(r$Rindex, c(0L, 1L))
  expect_identical(addCut(cp, mip, c(0L, 1L), c(1, 1), 2L, 1, TRUE, FALSE, FALSE)$index, -1L)
  expect_identical(addCut(cp, mip, c(0L, 1L), c(1, -1), 2L, 0, TRUE, FALSE, FALSE)$index, 1L)
  expect_identical(hmipC(".hmip_HighsCutPool_getNumCuts")(cp), 2L)
  dom <- new.env(); dom$col_lower_ <- c(0, 0); dom$col_upper_ <- c(1, 1)
  cs <- hmipC(".hmip_new_HighsCutSet")()
  hmipC(".hmip_HighsCutPool_separate")(cp, c(0.8, 0.8), dom, cs, 1e-6)
  # only x0 + x1 <= 1 is violated
  expect_identical(cs$cutindices, 0L)
  expect_identical(cs$upper_, 1)
})

test_that("HighsSparseVectorSum accumulates in double double", {
  vs <- hmipC(".hmip_new_HighsSparseVectorSum")(4L)
  add <- hmipC(".hmip_HighsSparseVectorSum_add")
  add(vs, 2L, 1); add(vs, 0L, 3); add(vs, 2L, 1e-20); add(vs, 2L, -1)
  expect_identical(hmipC(".hmip_HighsSparseVectorSum_getNonzeros")(vs), c(2L, 0L))
  expect_identical(vs$hi[3L] + vs$lo[3L], 1e-20)
  hmipC(".hmip_HighsSparseVectorSum_cleanup")(vs, function(i, v) abs(v) <= 1e-14)
  expect_identical(hmipC(".hmip_HighsSparseVectorSum_getNonzeros")(vs), 0L)
})

test_that("lifted knapsack cover of 3x0 + 3x1 + 3x2 <= 5 is x0 + x1 + x2 <= 1", {
  cg <- new.env()
  cg$feastol <- 1e-6; cg$epsilon <- 1e-9
  cg$randgen <- modelcompR:::.hs_rng(0)
  cg$rowlen <- 3L; cg$vals <- c(3, 3, 3); cg$inds <- 0:2
  cg$upper <- c(1, 1, 1); cg$solval <- c(1, 1, 0.5); cg$isintegral <- c(TRUE, TRUE, TRUE)
  cg$complementation <- integer(); cg$rhs <- c(5, 0)
  expect_true(hmipC(".hmip_HighsCutGeneration_determineCover")(cg, TRUE))
  expect_identical(cg$cover, c(0L, 1L))
  expect_identical(cg$lambda, c(1, 0))
  hmipC(".hmip_HighsCutGeneration_separateLiftedKnapsackCover")(cg)
  expect_identical(cg$vals, c(1, 1, 1))
  expect_identical(cg$rhs, c(1, 0))
})

test_that("libstdc++ heap emulation pops by the first member", {
  push <- hmipC(".hmip_gfk_heap_push"); pop <- hmipC(".hmip_gfk_heap_pop")
  h <- list(first = integer(), second = integer())
  for (q in seq_along(c(5L, 2L, 7L, 2L, 1L, 9L, 3L))) h <- push(h, c(5L, 2L, 7L, 2L, 1L, 9L, 3L)[q], q)
  out <- integer()
  while (length(h$first)) { out <- c(out, h$first[1L]); h <- pop(h) }
  expect_identical(out, c(1L, 2L, 2L, 3L, 5L, 7L, 9L))
})

gfk_solutions <- function(k, Aval, Aindex, Astart, numRow, rhsrow) {
  g <- hmipC(".hmip_new_HighsGFkSolve")()
  hmipC(".hmip_HighsGFkSolve_fromCSC")(g, k, Aval, Aindex, Astart, numRow)
  hmipC(".hmip_HighsGFkSolve_setRhs")(g, rhsrow, 1)
  sols <- list()
  hmipC(".hmip_HighsGFkSolve_solve")(g, function(w) sols[[length(sols) + 1L]] <<- w)
  sols
}

test_that("HighsGFkSolve finds the mod 2 and mod 3 solutions", {
  # GF(2): columns are integer rows over variables 0..2 and the rhs row 3
  sols <- gfk_solutions(2, c(1, 1, 1, 1, 1, 1, 1, 1, 1), c(0L, 1L, 3L, 1L, 2L, 3L, 0L, 2L, 3L),
                        c(0L, 3L, 6L, 9L), 4L, 3L)
  expect_gt(length(sols), 0L)
  for (w in sols) {
    o <- order(w$index)
    expect_identical(w$index[o], 0:2)
    expect_identical(w$weight[o], c(1, 1, 1))
  }
  # GF(3): w0 + 2 w1 = 0 and w0 + w1 = 1 give w = (2, 2)
  sols <- gfk_solutions(3, c(1, 1, 2, 1), c(0L, 1L, 0L, 1L), c(0L, 2L, 4L), 2L, 1L)
  expect_gt(length(sols), 0L)
  w <- sols[[1L]]
  o <- order(w$index)
  expect_identical(w$index[o], 0:1)
  expect_identical(w$weight[o], c(2, 2))
})

test_that("an inconsistent system reports no solution", {
  # w0 = 0 and w0 = 1 mod 2
  sols <- gfk_solutions(2, c(1, 1), c(0L, 1L), c(0L, 2L), 2L, 1L)
  expect_length(sols, 0L)
})
