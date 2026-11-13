# Group A of the HiGHS MIP port: HighsHashTree, HighsDomain, HighsConflictPool,
# HighsCliqueTable and HighsImplications on small hand checked models.

hmipA_stub <- function(cost, lower, upper, integrality, Astart, Aindex, Avalue, rlo, rup) {
  ncol <- length(cost); nrow <- length(rlo)
  model <- new.env(parent = emptyenv())
  model$num_col_ <- ncol; model$num_row_ <- nrow
  model$col_cost_ <- cost; model$col_lower_ <- lower; model$col_upper_ <- upper
  model$row_lower_ <- rlo; model$row_upper_ <- rup
  model$integrality_ <- integrality
  model$a_matrix_ <- list(start_ = Astart, index_ = Aindex, value_ = Avalue)
  # row-wise copy
  rs <- integer(nrow + 1L); ri <- integer(); rv <- numeric()
  for (r in seq_len(nrow) - 1L) {
    for (c in seq_len(ncol) - 1L) {
      s <- Astart[c + 1L]; e <- Astart[c + 2L]
      if (e > s) for (k in (s + 1L):e) if (Aindex[k] == r) { ri <- c(ri, c); rv <- c(rv, Avalue[k]) }
    }
    rs[r + 2L] <- length(ri)
  }
  ms <- new.env(parent = emptyenv())
  ms$model_ <- model
  md <- new.env(parent = emptyenv())
  md$feastol <- 1e-6; md$epsilon <- 1e-9; md$upper_limit <- Inf
  md$ARstart_ <- rs; md$ARindex_ <- ri; md$ARvalue_ <- rv
  md$integral_cols <- which(integrality != 0L) - 1L
  md$pruned_treeweight <- c(0, 0)
  ms$mipdata_ <- md
  md$cliquetable <- .hmip_new_HighsCliqueTable(ncol)
  md$domain <- .hmip_new_HighsDomain(ms)
  .hmip_HighsDomain_computeRowActivities(md$domain)
  ms
}

test_that("hash tree inserts, finds and erases like a set", {
  keys <- c(5L, 900L, 3L, 77L, 1L, 64L, 4096L, 12L, 300L, 2L, 31L, 63L)
  set.seed(1)
  big <- unique(c(keys, sample(0:5000, 200)))
  tree <- NULL
  for (k in big) {
    r <- .hmip_hht_insert(tree, k, k * 2L)
    expect_true(r$inserted)
    tree <- r$node
  }
  expect_false(.hmip_hht_insert(tree, big[3], 0L)$inserted)
  e <- .hmip_hht_entries(tree)
  expect_setequal(e$k, big)
  expect_equal(unlist(e$v), e$k * 2L)
  for (k in big) expect_equal(.hmip_hht_find(tree, k)$value, k * 2L)
  expect_null(.hmip_hht_find(tree, 99999L))
  half <- big[seq(1, length(big), 2)]
  for (k in half) tree <- .hmip_hht_erase(tree, k)
  expect_setequal(.hmip_hht_keys(tree), setdiff(big, half))
  t2 <- NULL
  for (k in c(7000L, 7001L, big[2])) t2 <- .hmip_hht_insert(t2, k)$node
  expect_equal(.hmip_hht_find_common(tree, t2)$key, big[2])
  for (k in setdiff(big, half)) tree <- .hmip_hht_erase(tree, k)
  expect_true(.hmip_hht_empty(tree))
})

test_that("domain propagates a packing row, backtracks and keeps the stack", {
  # x0 + x1 + x2 <= 1, binaries; row 1: x0 + y3 >= 1.5 with y3 in [0, 1]
  ms <- hmipA_stub(c(1, 1, 1, 0), c(0, 0, 0, 0), c(1, 1, 1, 1), c(1L, 1L, 1L, 0L),
                   c(0L, 2L, 3L, 4L, 5L), c(0L, 1L, 0L, 0L, 1L), c(1, 1, 1, 1, 1), c(-Inf, 1.5), c(1, Inf))
  dom <- ms$mipdata_$domain
  expect_equal(.hmip_HighsDomain_getMinActivity(dom, 0L), 0)
  expect_equal(.hmip_HighsDomain_getMaxActivity(dom, 1L), 2)
  .hmip_HighsDomain_changeBound(dom, .hmip_domchg(1, 0L, 0L))
  expect_true(.hmip_HighsDomain_propagate(dom))
  expect_false(dom$infeasible_)
  expect_equal(dom$col_upper_[2:3], c(0, 0))
  expect_equal(dom$col_lower_[4], 0.5)
  expect_equal(dom$domchgreason_type[1], -1L)
  expect_true(all(dom$domchgreason_type[-1] %in% c(-3L, -4L)))
  r <- .hmip_HighsDomain_getColUpperPos(dom, 1L, 0L)
  expect_equal(r$bound, 1); expect_equal(r$pos, -1L)
  b <- .hmip_HighsDomain_backtrack(dom)
  expect_equal(b$column, 0L)
  expect_equal(dom$col_upper_, c(1, 1, 1, 1))
  expect_equal(dom$col_lower_, c(0, 0, 0, 0))
  expect_length(dom$domchgstack_column, 0L)
  expect_equal(.hmip_HighsDomain_getMinActivity(dom, 0L), 0)
})

test_that("conflicting bounds mark the domain infeasible and backtrack clears it", {
  ms <- hmipA_stub(c(1, 1), c(0, 0), c(1, 1), c(1L, 1L), c(0L, 1L, 2L), c(0L, 0L), c(1, 1), -Inf, 1)
  dom <- ms$mipdata_$domain
  .hmip_HighsDomain_changeBound(dom, .hmip_domchg(1, 0L, 0L))
  .hmip_HighsDomain_changeBound(dom, .hmip_domchg(1, 1L, 0L))
  expect_true(dom$infeasible_)
  expect_equal(dom$infeasible_reason$type, -3L)
  .hmip_HighsDomain_backtrack(dom)
  expect_false(dom$infeasible_)
  .hmip_HighsDomain_backtrackToGlobal(dom)
  expect_equal(dom$col_lower_, c(0, 0))
})

test_that("clique table stores cliques and fixes through implications", {
  ms <- hmipA_stub(c(1, 1, 1, 1), rep(0, 4), rep(1, 4), rep(1L, 4),
                   c(0L, 0L, 0L, 0L, 0L), integer(), numeric(), numeric(), numeric())
  md <- ms$mipdata_; ct <- md$cliquetable; dom <- md$domain
  .hmip_HighsCliqueTable_doAddClique(ct, c(.hmip_cv(0L, 1L), .hmip_cv(1L, 1L), .hmip_cv(2L, 1L)))
  .hmip_HighsCliqueTable_doAddClique(ct, c(.hmip_cv(2L, 0L), .hmip_cv(3L, 1L)))
  expect_equal(.hmip_HighsCliqueTable_numCliques(ct), 2L)
  expect_equal(ct$numEntries, 5L)
  expect_equal(.hmip_HighsCliqueTable_findCommonClique(ct, .hmip_cv(0L, 1L), .hmip_cv(2L, 1L)),
               c(.hmip_cv(0L, 1L), .hmip_cv(1L, 1L), .hmip_cv(2L, 1L)))
  expect_null(.hmip_HighsCliqueTable_findCommonClique(ct, .hmip_cv(0L, 1L), .hmip_cv(3L, 1L)))
  expect_equal(.hmip_HighsCliqueTable_findCommonCliqueId(ct, .hmip_cv(3L, 1L), .hmip_cv(2L, 0L)), 1L)
  expect_equal(.hmip_HighsCliqueTable_getNumImplications(ct, 2L), 3L)
  # fixing x0 = 1 fixes x1 = x2 = 0 through the clique table, which fixes x3 = 0
  .hmip_HighsDomain_changeBound(dom, .hmip_domchg(1, 0L, 0L))
  expect_equal(dom$col_upper_, c(1, 0, 0, 0))
  expect_true(all(dom$domchgreason_type[-1] == -5L))
  .hmip_HighsDomain_backtrack(dom)
  expect_equal(dom$col_upper_, rep(1, 4))
  .hmip_HighsCliqueTable_removeClique(ct, 0L)
  expect_equal(.hmip_HighsCliqueTable_numCliques(ct), 1L)
  expect_equal(ct$free_first, 3L); expect_equal(ct$free_second, 0L)
  .hmip_HighsCliqueTable_doAddClique(ct, c(.hmip_cv(0L, 0L), .hmip_cv(1L, 0L)))
  expect_equal(ct$c_start[1], 0L)
  expect_equal(ct$free_first, 1L); expect_equal(ct$free_second, 2L)
})

test_that("conflict pool adds, propagates and ages out a conflict", {
  ms <- hmipA_stub(c(1, 1, 1), rep(0, 3), rep(1, 3), rep(1L, 3),
                   c(0L, 0L, 0L, 0L), integer(), numeric(), numeric(), numeric())
  dom <- ms$mipdata_$domain
  pool <- .hmip_new_HighsConflictPool(3L, 100L)
  .hmip_HighsDomain_addConflictPool(dom, pool)
  # conflict: x0 >= 1 and x1 >= 1 cannot both hold
  fr <- list(pos = c(0L, 1L), col = c(0L, 1L), bt = c(0L, 0L), bv = c(1, 1))
  .hmip_HighsConflictPool_addConflictCut(pool, dom, fr)
  expect_equal(.hmip_HighsConflictPool_getNumConflicts(pool), 1L)
  # upstream linkWatchedLiteral only moves the head of a nonempty list, so
  # watches never link and a bound change does not wake the conflict
  expect_equal(dom$conflictPoolPropagation[[1]]$colLowerWatched_, c(-1L, -1L, -1L))
  .hmip_HighsDomain_changeBound(dom, .hmip_domchg(1, 0L, 0L))
  .hmip_HighsDomain_propagate(dom)
  expect_equal(dom$col_upper_[2], 1)
  # a conflict added while x0 >= 1 holds is propagated on addition
  .hmip_HighsConflictPool_addConflictCut(pool, dom, fr)
  .hmip_HighsDomain_propagate(dom)
  expect_equal(dom$col_upper_[2], 0)
  expect_equal(dom$domchgreason_type[2], 0L)
  .hmip_HighsDomain_backtrack(dom)
  for (i in 1:4) .hmip_HighsConflictPool_performAging(pool)
  expect_equal(.hmip_HighsConflictPool_getNumConflicts(pool), 0L)
  expect_equal(pool$range_first, c(-1L, -1L))
  expect_equal(pool$free_first, c(2L, 2L)); expect_equal(pool$free_second, c(0L, 2L))
})

test_that("implications store variable bounds and keep the stronger one", {
  ms <- hmipA_stub(c(1, 1), c(0, 0), c(1, 10), c(1L, 0L),
                   c(0L, 0L, 0L), integer(), numeric(), numeric(), numeric())
  im <- .hmip_new_HighsImplications(ms)
  ms$mipdata_$implications <- im
  .hmip_HighsImplications_addVUB(im, 1L, 0L, 5, 3)
  expect_equal(im$numVarBounds, 1)
  .hmip_HighsImplications_addVUB(im, 1L, 0L, 8, 1)
  expect_equal(.hmip_hht_find(im$vubs[[2]], 0L)$value, c(8, 1))
  .hmip_HighsImplications_addVUB(im, 1L, 0L, 1, 9)
  expect_equal(.hmip_hht_find(im$vubs[[2]], 0L)$value, c(8, 1))
  .hmip_HighsImplications_addVUB(im, 1L, 0L, 0, 10)
  expect_equal(im$numVarBounds, 1)
  r <- .hmip_HighsImplications_cleanupVub(im, 1L, 0L, list(coef = 20, constant = 1), 10)
  expect_equal(r$vub, list(coef = 9, constant = 1))
  expect_false(r$redundant)
})

test_that("clique algorithms find, partition and merge a triangle", {
  ms <- hmipA_stub(c(1, 1, 1, 1), rep(0, 4), rep(1, 4), rep(1L, 4),
                   c(0L, 0L, 0L, 0L, 0L), integer(), numeric(), numeric(), numeric())
  md <- ms$mipdata_; ct <- md$cliquetable
  v <- c(.hmip_cv(0L, 1L), .hmip_cv(1L, 1L), .hmip_cv(2L, 1L), .hmip_cv(3L, 1L))
  .hmip_HighsCliqueTable_doAddClique(ct, v[c(1, 2)], 2L, FALSE, 0L)
  .hmip_HighsCliqueTable_doAddClique(ct, v[c(2, 3)], 2L, FALSE, 0L)
  .hmip_HighsCliqueTable_doAddClique(ct, v[c(1, 3)], 2L, FALSE, 0L)
  expect_equal(.hmip_HighsCliqueTable_cliquesSize(ct), 3L)
  expect_equal(.hmip_HighsCliqueTable_computeMaximalCliques(ct, v, 1e-6), list(v[1:3]))
  p <- .hmip_HighsCliqueTable_cliquePartition(ct, v)
  expect_equal(p$clqVars, c(1L, 5L, 3L, 7L))
  expect_equal(p$partitionStart, c(0L, 3L, 4L))
  .hmip_HighsCliqueTable_runCliqueMerging(ct, md$domain)
  expect_equal(.hmip_HighsCliqueTable_numCliques(ct), 1L)
  k <- which(ct$c_start != -1L)
  expect_setequal(ct$cliqueentries[(ct$c_start[k] + 1L):ct$c_end[k]], v[1:3])
})

test_that("cleanupVarbounds tightens a VUB and drops a redundant VLB", {
  ms <- hmipA_stub(c(1, 1), c(0, 0), c(1, 10), c(1L, 0L),
                   c(0L, 0L, 0L), integer(), numeric(), numeric(), numeric())
  im <- .hmip_new_HighsImplications(ms)
  ms$mipdata_$implications <- im
  .hmip_HighsImplications_addVUB(im, 1L, 0L, 20, 1)
  .hmip_HighsImplications_addVLB(im, 1L, 0L, -5, 0)
  r <- .hmip_HighsImplications_cleanupVlb(im, 1L, 0L, list(coef = -5, constant = 0), 0)
  expect_true(r$redundant)
  .hmip_HighsImplications_cleanupVarbounds(im, 1L)
  expect_equal(.hmip_hht_find(im$vubs[[2]], 0L)$value, c(9, 1))
  expect_null(.hmip_hht_find(im$vlbs[[2]], 0L))
})

test_that("destroy deregisters a domain copy from its pools and is idempotent", {
  ms <- hmipA_stub(c(1, 1), c(0, 0), c(1, 1), c(1L, 1L),
                   c(0L, 0L, 0L), integer(), numeric(), numeric(), numeric())
  dom <- ms$mipdata_$domain
  pool <- .hmip_new_HighsConflictPool(3L, 100L)
  .hmip_HighsDomain_addConflictPool(dom, pool)
  expect_length(pool$propagationDomains, 1L)
  cp <- .hmip_HighsDomain_copy(dom)
  expect_length(pool$propagationDomains, 2L)
  .hmip_HighsDomain_destroy(cp)
  expect_length(pool$propagationDomains, 1L)
  expect_identical(pool$propagationDomains[[1]], dom$conflictPoolPropagation[[1]])
  .hmip_HighsDomain_destroy(cp)
  expect_length(pool$propagationDomains, 1L)
})
