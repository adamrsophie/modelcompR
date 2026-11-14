# HiGHS 1.12.0 symmetry detection and handling (presolve/HighsSymmetry.cpp),
# R/optimization-scipy-highs-mip-symmetry.R. Expected values are computed by
# hand from the C++ logic.

ns <- asNamespace("modelcompR")
has_sort <- function() {
  exists(".hmip_pdqsort", envir = ns) && exists(".hmip_std_partition", envir = ns)
}

symF_binpack <- function(size, nb, cap) {
  ni <- length(size)
  nc <- ni * nb + nb
  nr <- ni + nb
  A <- matrix(0, nr, nc)
  for (i in seq_len(ni)) A[i, (i - 1) * nb + seq_len(nb)] <- 1
  for (b in seq_len(nb)) {
    A[ni + b, seq(b, ni * nb, nb)] <- size
    A[ni + b, ni * nb + b] <- -cap
  }
  start <- 0L
  index <- integer(0)
  value <- numeric(0)
  for (j in seq_len(nc)) {
    nz <- which(A[, j] != 0)
    index <- c(index, nz - 1L)
    value <- c(value, A[nz, j])
    start <- c(start, length(index))
  }
  list(num_col_ = nc, num_row_ = nr, col_cost_ = c(rep(0, ni * nb), rep(1, nb)), col_lower_ = rep(0, nc),
       col_upper_ = rep(1, nc), row_lower_ = c(rep(1, ni), rep(-Inf, nb)), row_upper_ = c(rep(1, ni), rep(0, nb)),
       integrality_ = rep(1L, nc), a_matrix_ = list(start_ = start, index_ = index, value_ = value), A = A)
}

symF_detect <- function(m) {
  sd <- .hmip_new_HighsSymmetryDetection()
  .hmip_HighsSymmetryDetection_loadModelAsGraph(sd, m, 1e-9)
  ok <- .hmip_HighsSymmetryDetection_initializeDetection(sd)
  sym <- .hmip_new_HighsSymmetries()
  if (ok) .hmip_HighsSymmetryDetection_run(sd, sym)
  list(ok = ok, sym = sym)
}

test_that("M31 arithmetic follows HighsHashHelpers", {
  expect_identical(.hmip_mulM31(2147483646, 2147483646), 1)
  expect_identical(.hmip_mulM31(2^30, 2), 1)
  # 123456789 * 987654321 = 121932631112635269 (exact, needs more than 53 bits)
  expect_identical(.hmip_mulM31(123456789, 987654321), 2137109934)
  # modexp_M31 squares before testing the exponent bit, as written
  expect_identical(.hmip_modexpM31(3, 2), 9)
  expect_identical(.hmip_modexpM31(3, 3), 27)
  expect_identical(.hmip_modexpM31(3, 4), 81)
  expect_identical(.hmip_modexpM31(3, 6), 243)
  expect_identical(.hmip_sparse_combine32(2147483640, 10), 3)
})

test_that("HighsDisjointSets merges by size", {
  ds <- .hmip_new_HighsDisjointSets(5L)
  .hmip_HighsDisjointSets_merge(ds, 0L, 1L)
  expect_identical(.hmip_HighsDisjointSets_getSet(ds, 0L), 1L)
  .hmip_HighsDisjointSets_merge(ds, 2L, 0L)
  expect_identical(.hmip_HighsDisjointSets_getSet(ds, 2L), 1L)
  expect_identical(.hmip_HighsDisjointSets_getSetSize(ds, 1L), 3L)
})

test_that("HighsMatrixColoring numbers values in insertion order", {
  cm <- .hmip_new_HighsMatrixColoring(1e-9)
  expect_identical(.hmip_HighsMatrixColoring_color(cm, 0), 1)
  expect_identical(.hmip_HighsMatrixColoring_color(cm, 1 + 1e-12), 2)
  expect_identical(.hmip_HighsMatrixColoring_color(cm, -Inf), 3)
  expect_identical(.hmip_HighsMatrixColoring_color(cm, 6), 5)
  expect_identical(.hmip_HighsMatrixColoring_color(cm, -10), 6)
  expect_identical(.hmip_HighsMatrixColoring_color(cm, 6), 5)
})

test_that("symmetry detection finds the bin and item symmetries", {
  skip_if_not(has_sort())
  m <- symF_binpack(c(6, 5, 4, 4, 3, 3), 4, 10)
  r <- symF_detect(m)
  s <- r$sym
  expect_true(r$ok)
  # 3 bin transpositions and the swaps of the two pairs of equal items: one
  # component, so no full orbitope
  expect_identical(s$numGenerators, 5L)
  expect_identical(s$numPerms, 5L)
  expect_length(s$orbitopes, 0L)
  expect_identical(s$permutationColumns, 0:27)
  for (p in seq_len(s$numPerms) - 1L) {
    perm <- s$permutations[p * 28L + seq_len(28L)]
    key <- function(M) sort(apply(cbind(M, m$row_lower_, m$row_upper_), 1L, paste, collapse = ","))
    expect_identical(key(m$A[, perm + 1L]), key(m$A))
  }
})

test_that("a fixed orbitope block next to a generator block compresses the permutations", {
  skip_if_not(has_sort())
  m1 <- symF_binpack(c(7, 6, 5), 3, 10)
  m2 <- symF_binpack(c(6, 5, 4, 4, 3, 3), 4, 10)
  A <- rbind(cbind(m1$A, matrix(0, nrow(m1$A), ncol(m2$A))), cbind(matrix(0, nrow(m2$A), ncol(m1$A)), m2$A))
  nc <- ncol(A)
  start <- 0L
  index <- integer(0)
  value <- numeric(0)
  for (j in seq_len(nc)) {
    nz <- which(A[, j] != 0)
    index <- c(index, nz - 1L)
    value <- c(value, A[nz, j])
    start <- c(start, length(index))
  }
  m <- list(num_col_ = nc, num_row_ = nrow(A), col_cost_ = c(m1$col_cost_, m2$col_cost_), col_lower_ = rep(0, nc),
            col_upper_ = rep(1, nc), row_lower_ = c(m1$row_lower_, m2$row_lower_),
            row_upper_ = c(m1$row_upper_, m2$row_upper_), integrality_ = rep(1L, nc),
            a_matrix_ = list(start_ = start, index_ = index, value_ = value))
  s <- symF_detect(m)$sym
  expect_length(s$orbitopes, 1L)
  expect_identical(s$numGenerators, 7L)
  expect_identical(s$numPerms, 5L)
  expect_identical(s$permutationColumns, 12:39)
  expect_length(s$columnPosition, nc)
  expect_false(anyNA(s$permutations))
  expect_length(s$permutations, 5L * 28L)
  key <- function(M) sort(apply(cbind(M, m$row_lower_, m$row_upper_), 1L, paste, collapse = ","))
  for (p in 0:4) {
    full <- 0:(nc - 1L)
    full[s$permutationColumns + 1L] <- s$permutations[p * 28L + seq_len(28L)]
    expect_identical(key(A[, full + 1L]), key(A))
  }
})

test_that("milp_binpackmix matches the Linux log (Found 12 generator(s))", {
  skip_if_not(has_sort())
  s <- symF_detect(symF_binpack(c(4, 4, 4, 4, 4, 6, 6, 7, 7), 7, 10))$sym
  expect_identical(s$numGenerators, 12L)
  expect_length(s$orbitopes, 0L)
})

test_that("identical bins with distinct items form a full orbitope", {
  skip_if_not(has_sort())
  r <- symF_detect(symF_binpack(c(7, 6, 5), 3, 10))
  s <- r$sym
  expect_identical(s$numGenerators, 2L)
  expect_identical(s$numPerms, 0L)
  expect_length(s$orbitopes, 1L)
  om <- s$orbitopes[[1L]]
  expect_identical(c(om$numRows, om$rowLength), c(4L, 3L))
  # matrix order as this port produces it (not yet checked against Linux)
  expect_identical(om$matrix, c(1L, 4L, 7L, 10L, 2L, 5L, 8L, 11L, 0L, 3L, 6L, 9L))
  expect_identical(.hmip_HighsSymmetries_columnToOrbitope_size(s), 12L)
})

# a minimal HighsDomain stand-in with the group A member and method names
symF_mock_domain <- function(n) {
  d <- new.env()
  d$col_lower_ <- rep(0, n)
  d$col_upper_ <- rep(1, n)
  d$branchPos_ <- integer(0)
  d$domchgstack_column <- integer(0)
  d$domchgstack_boundtype <- integer(0)
  d$infeasible_ <- FALSE
  d
}
symF_with_mocks <- function(code) {
  mocks <- list(
    .hmip_HighsDomain_changeBound3 = function(d, bt, col, val, reason) {
      if (bt == 0L) d$col_lower_[col + 1L] <- val else d$col_upper_[col + 1L] <- val
      d$domchgstack_column <- c(d$domchgstack_column, col)
      d$domchgstack_boundtype <- c(d$domchgstack_boundtype, bt)
      if (d$col_lower_[col + 1L] > d$col_upper_[col + 1L]) d$infeasible_ <- TRUE
    },
    .hmip_HighsDomain_propagate = function(d) TRUE,
    .hmip_HighsDomain_markInfeasible = function(d, reason) d$infeasible_ <- TRUE,
    .hmip_HighsDomain_isGlobalBinary = function(d, col) TRUE,
    .hmip_HighsDomain_variableType = function(d, col) 1L)
  # mock inside the namespace when group A's functions exist, else provide
  # them through the global environment the namespace falls back to
  added <- character(0)
  inns <- names(mocks)[vapply(names(mocks), exists, NA, envir = ns, inherits = FALSE)]
  if (length(inns)) do.call(testthat::local_mocked_bindings, c(mocks[inns], list(.package = "modelcompR", .env = environment())))
  for (nm in setdiff(names(mocks), inns)) {
    assign(nm, mocks[[nm]], envir = globalenv())
    added <- c(added, nm)
  }
  on.exit(rm(list = added, envir = globalenv()), add = TRUE)
  code
}
symF_branch <- function(d, bt, col, val) {
  d$branchPos_ <- c(d$branchPos_, length(d$domchgstack_column))
  .hmip_HighsDomain_changeBound3(d, bt, col, val, NULL)
}

test_that("full and packing orbitope fixing", {
  skip_if_not(has_sort())
  symF_with_mocks({
    s <- symF_detect(symF_binpack(c(7, 6, 5), 3, 10))$sym
    om <- s$orbitopes[[1L]]
    om$rowIsSetPacking <- rep(0L, 4L)
    for (j in 0:2) for (i in 0:3) {
      om$columnToRow <- .hmip_sym_insert(om$columnToRow, .hmip_HighsOrbitopeMatrix_entry(om, i, j), i)
    }
    # x(0, 1) = 1 forces x(0, 0) = 1 in a lexicographically sorted orbitope
    d <- symF_mock_domain(12L)
    symF_branch(d, 0L, .hmip_HighsOrbitopeMatrix_entry(om, 0L, 1L), 1)
    expect_identical(.hmip_HighsSymmetries_propagateOrbitopes(s, d), 1L)
    expect_identical(d$col_lower_[2L], 1)
    # packing rows: the first branched row may hold its 1 only in column 0
    om$rowIsSetPacking <- rep(1L, 4L)
    om$numSetPackingRows <- 4L
    d <- symF_mock_domain(12L)
    symF_branch(d, 0L, .hmip_HighsOrbitopeMatrix_entry(om, 1L, 1L), 1)
    expect_identical(.hmip_HighsSymmetries_propagateOrbitopes(s, d), 0L)
    expect_true(d$infeasible_)
    expect_identical(.hmip_HighsSymmetries_getBranchingColumn(s, rep(0, 12), rep(1, 12), 0L), 1L)
  })
})

test_that("stabilizer orbits and orbital fixing", {
  skip_if_not(has_sort())
  symF_with_mocks({
    s <- symF_detect(symF_binpack(c(6, 5, 4, 4, 3, 3), 4, 10))$sym
    d <- symF_mock_domain(28L)
    st <- .hmip_HighsSymmetries_computeStabilizerOrbits(s, d)
    expect_identical(st$orbitStarts, c(0L, 4L, 8L, 16L, 24L, 28L))
    expect_identical(sort(st$orbitCols[25:28]), 24:27)
    symF_branch(d, 1L, 24L, 0)
    expect_identical(.hmip_StabilizerOrbits_orbitalFixing(st, d), 3L)
    expect_identical(d$col_upper_[25:28], c(0, 0, 0, 0))
  })
})
