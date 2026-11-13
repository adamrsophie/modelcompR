# Symmetry detection and symmetry handling of the HiGHS 1.12.0 MIP solver
# (presolve/HighsSymmetry.cpp and .h, the M31 hashing of util/HighsHash.h
# and util/HighsDisjointSets.h), as SciPy 1.18.1 vendors them. MIT license
# notice in R/optimization-scipy-highs-util.R.
#
# Conventions (scratchpad API-F.md): every index value stored in a vector is
# 0-based as in C++ (vertices, columns, cells, positions), and an R vector
# indexed by such a value is read as v[i + 1L]. The C++ pair vector Gedge is
# held as two vectors Gfirst and Gsecond. Hash tables keyed by small
# nonnegative integers are vectors with NA for an absent key; HiGHS never
# iterates over any of the hash tables used here, so their internal order
# does not matter. HighsSymmetryDetection::run executes in a separate task
# in HiGHS; the task only reads the presolved model and writes its own
# HighsSymmetries, so the order in which it runs does not change results
# and the port runs it inline.

# ---------------------------------------------------------------------------
# HighsHashHelpers modulo the Mersenne prime M31 = 2^31 - 1

.hmip_M31 <- 2147483647

# multiply_modM31(a, b) for a, b < 2^31: the 62 bit product split so that
# every partial result is exact in a double
.hmip_mulM31 <- function(a, b) {
  bl <- b %% 65536
  bh <- (b - bl) / 65536
  p1 <- a * bl
  p2 <- a * bh
  q <- p2 %% 32768
  t <- q * 65536 + p1
  lo <- t %% 2147483648
  r <- (p2 - q) / 32768 + (t - lo) / 2147483648 + lo
  if (r >= .hmip_M31) r <- r - .hmip_M31
  r
}

# modexp_M31 as written (it squares first, like modexp_M61)
.hmip_modexpM31 <- function(a, e) {
  res <- a
  while (e != 1) {
    res <- .hmip_mulM31(res, res)
    if (e %% 2 == 1) res <- .hmip_mulM31(res, a)
    e <- e %/% 2
  }
  res
}

# the multiplier of sparse_combine32(hash, index, value):
# ((pair_hash<0>(value, 0) >> 33) | 1) * modexp_M31(c[index & 63] & M31,
# (index >> 6) + 1) modulo M31
.hmip_combine32_term <- function(index, value) {
  v <- .hs_u64_shift(.hs_pair_hash(0L, value, 0), 33L)
  v <- v[1L] + 65536 * v[2L]
  if (v %% 2 == 0) v <- v + 1
  a <- .hs_lo32(.hs_hash_cl()[[bitwAnd(index, 63L) + 1L]]) %% 2147483648
  .hmip_mulM31(v, .hmip_modexpM31(a, index %/% 64 + 1))
}

.hmip_sparse_combine32 <- function(hash, term) {
  r <- hash + term
  r <- r %/% 2147483648 + r %% 2147483648
  if (r >= .hmip_M31) r <- r - .hmip_M31
  r
}

# ---------------------------------------------------------------------------
# HighsDisjointSets<false> (sizes decide the representative), 0-based items

.hmip_new_HighsDisjointSets <- function(numItems) {
  ds <- new.env(parent = emptyenv())
  ds$sizes <- rep(1L, numItems)
  ds$sets <- seq_len(numItems) - 1L
  ds
}

.hmip_HighsDisjointSets_getSet <- function(ds, item) {
  repr <- ds$sets[item + 1L]
  if (repr != ds$sets[repr + 1L]) {
    stack <- integer(0)
    repeat {
      stack <- c(stack, item)
      item <- repr
      repr <- ds$sets[repr + 1L]
      if (repr == ds$sets[repr + 1L]) break
    }
    ds$sets[stack + 1L] <- repr
    ds$sets[item + 1L] <- repr
  }
  repr
}

.hmip_HighsDisjointSets_getSetSize <- function(ds, set) ds$sizes[set + 1L]

.hmip_HighsDisjointSets_merge <- function(ds, item1, item2) {
  r1 <- .hmip_HighsDisjointSets_getSet(ds, item1)
  r2 <- .hmip_HighsDisjointSets_getSet(ds, item2)
  if (r1 == r2) return(invisible(NULL))
  if (ds$sizes[r1 + 1L] > ds$sizes[r2 + 1L]) {
    ds$sets[r2 + 1L] <- r1
    ds$sizes[r1 + 1L] <- ds$sizes[r1 + 1L] + ds$sizes[r2 + 1L]
  } else {
    ds$sets[r1 + 1L] <- r2
    ds$sizes[r2 + 1L] <- ds$sizes[r2 + 1L] + ds$sizes[r1 + 1L]
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# HighsMatrixColoring: std::map<double, u32> as sorted keys with colors

.hmip_new_HighsMatrixColoring <- function(tolerance) {
  cm <- new.env(parent = emptyenv())
  cm$keys <- c(-Inf, 0, 1, Inf)
  cm$colors <- c(3, 1, 2, 4)
  cm$tolerance <- tolerance
  cm
}

.hmip_HighsMatrixColoring_color <- function(cm, value) {
  keys <- cm$keys
  # lower_bound(value - tolerance)
  pos <- sum(keys < value - cm$tolerance) + 1L
  if (pos > length(keys) || keys[pos] > value + cm$tolerance) {
    col <- length(keys) + 1
    cm$keys <- append(keys, value, pos - 1L)
    cm$colors <- append(cm$colors, col, pos - 1L)
    return(col)
  }
  cm$colors[pos]
}

# ---------------------------------------------------------------------------
# small hash table helpers (integer keys, NA for absent)

.hmip_sym_find <- function(tab, key) if (key + 1L > length(tab)) NA_integer_ else tab[key + 1L]

.hmip_sym_insert <- function(tab, key, value) {
  if (key + 1L > length(tab)) tab[(length(tab) + 1L):(key + 1L)] <- NA_integer_
  if (is.na(tab[key + 1L])) tab[key + 1L] <- value
  tab
}

# ---------------------------------------------------------------------------
# HighsSymmetries

.hmip_new_HighsSymmetries <- function() {
  sym <- new.env(parent = emptyenv())
  .hmip_HighsSymmetries_clear(sym)
  sym
}

.hmip_HighsSymmetries_clear <- function(sym) {
  sym$permutationColumns <- integer(0)
  sym$permutations <- integer(0)
  sym$orbitPartition <- integer(0)
  sym$orbitSize <- integer(0)
  sym$columnPosition <- integer(0)
  sym$linkCompressionStack <- integer(0)
  sym$columnToOrbitope <- integer(0)
  sym$orbitopes <- list()
  sym$numPerms <- 0L
  sym$numGenerators <- 0L
  invisible(sym)
}

.hmip_HighsSymmetries_columnToOrbitope_size <- function(sym) sum(!is.na(sym$columnToOrbitope))

.hmip_HighsSymmetries_mergeOrbits <- function(sym, v1, v2) {
  if (v1 == v2) return(invisible(NULL))
  o1 <- .hmip_HighsSymmetries_getOrbit(sym, v1)
  o2 <- .hmip_HighsSymmetries_getOrbit(sym, v2)
  if (o1 == o2) return(invisible(NULL))
  if (sym$orbitSize[o2 + 1L] < sym$orbitSize[o1 + 1L]) {
    sym$orbitPartition[o2 + 1L] <- o1
    sym$orbitSize[o1 + 1L] <- sym$orbitSize[o1 + 1L] + sym$orbitSize[o2 + 1L]
  } else {
    sym$orbitPartition[o1 + 1L] <- o2
    sym$orbitSize[o2 + 1L] <- sym$orbitSize[o2 + 1L] + sym$orbitSize[o1 + 1L]
  }
  invisible(NULL)
}

.hmip_HighsSymmetries_getOrbit <- function(sym, col) {
  i <- sym$columnPosition[col + 1L]
  if (i == -1L) return(-1L)
  op <- sym$orbitPartition
  orbit <- op[i + 1L]
  if (orbit != op[orbit + 1L]) {
    stack <- integer(0)
    repeat {
      stack <- c(stack, i)
      i <- orbit
      orbit <- op[orbit + 1L]
      if (orbit == op[orbit + 1L]) break
    }
    op[stack + 1L] <- orbit
    sym$orbitPartition <- op
  }
  orbit
}

.hmip_HighsSymmetries_getBranchingColumn <- function(sym, colLower, colUpper, col) {
  if (.hmip_HighsSymmetries_columnToOrbitope_size(sym) == 0L) return(col)
  o <- .hmip_sym_find(sym$columnToOrbitope, col)
  if (is.na(o) || sym$orbitopes[[o + 1L]]$numSetPackingRows == 0L) return(col)
  .hmip_HighsOrbitopeMatrix_getBranchingColumn(sym$orbitopes[[o + 1L]], colLower, colUpper, col)
}

.hmip_HighsSymmetries_propagateOrbitopes <- function(sym, domain) {
  if (.hmip_HighsSymmetries_columnToOrbitope_size(sym) == 0L || length(domain$branchPos_) == 0L) return(0L)
  # std::set: ascending unique orbitope indices
  orbs <- integer(0)
  for (pos in domain$branchPos_) {
    o <- .hmip_sym_find(sym$columnToOrbitope, domain$domchgstack_column[pos + 1L])
    if (!is.na(o)) orbs <- c(orbs, o)
  }
  numFixed <- 0L
  for (i in sort(unique(orbs))) {
    numFixed <- numFixed + .hmip_HighsOrbitopeMatrix_orbitalFixing(sym$orbitopes[[i + 1L]], domain)
    if (domain$infeasible_) break
  }
  numFixed
}

.hmip_HighsSymmetries_computeStabilizerOrbits <- function(sym, localdom) {
  st <- new.env(parent = emptyenv())
  st$orbitCols <- integer(0)
  st$orbitStarts <- integer(0)
  st$stabilizedCols <- integer(0)
  pc <- sym$permutationColumns
  cp <- sym$columnPosition
  stab <- integer(0)
  for (i in localdom$branchPos_) {
    col <- localdom$domchgstack_column[i + 1L]
    if (cp[col + 1L] == -1L) next
    if (localdom$domchgstack_boundtype[i + 1L] == 0L || !.hmip_HighsDomain_isGlobalBinary(localdom, col)) {
      stab <- c(stab, cp[col + 1L])
    }
  }
  permLength <- length(pc)
  sym$orbitPartition <- seq_len(permLength) - 1L
  sym$orbitSize <- rep(1L, permLength)
  for (i in seq_len(sym$numPerms) - 1L) {
    perm <- sym$permutations[i * permLength + seq_len(permLength)]
    if (any(pc[stab + 1L] != perm[stab + 1L])) next
    for (j in seq_len(permLength)) .hmip_HighsSymmetries_mergeOrbits(sym, pc[j], perm[j])
  }
  stab <- integer(0)
  orbitCols <- integer(0)
  for (i in seq_len(permLength)) {
    col <- pc[i]
    if (.hmip_HighsDomain_variableType(localdom, col) == 0L) next
    orbit <- .hmip_HighsSymmetries_getOrbit(sym, col)
    if (sym$orbitSize[orbit + 1L] == 1L) {
      stab <- c(stab, col)
    } else if (.hmip_HighsDomain_isGlobalBinary(localdom, col)) {
      orbitCols <- c(orbitCols, col)
    }
  }
  st$symmetries <- sym
  st$stabilizedCols <- sort(stab)
  if (length(orbitCols)) {
    orbitCols <- .hmip_pdqsort(orbitCols, function(c1, c2) {
      .hmip_HighsSymmetries_getOrbit(sym, c1) < .hmip_HighsSymmetries_getOrbit(sym, c2)
    })
    n <- length(orbitCols)
    starts <- 0L
    if (n > 1L) for (i in 1:(n - 1L)) {
      if (.hmip_HighsSymmetries_getOrbit(sym, orbitCols[i + 1L]) != .hmip_HighsSymmetries_getOrbit(sym, orbitCols[i])) {
        starts <- c(starts, i)
      }
    }
    st$orbitCols <- orbitCols
    st$orbitStarts <- c(starts, n)
  }
  st
}

# ---------------------------------------------------------------------------
# StabilizerOrbits

.hmip_StabilizerOrbits_isStabilized <- function(st, col) {
  st$symmetries$columnPosition[col + 1L] == -1L || col %in% st$stabilizedCols
}

.hmip_StabilizerOrbits_orbitalFixing <- function(st, domain) {
  numFixed <- .hmip_HighsSymmetries_propagateOrbitopes(st$symmetries, domain)
  if (domain$infeasible_ || length(st$orbitCols) == 0L) return(numFixed)
  oc <- st$orbitCols
  os <- st$orbitStarts
  numOrbits <- length(os) - 1L
  unspecified <- list(type = -2L, index = 0L)
  i <- 0L
  while (i < numOrbits) {
    range <- seq.int(os[i + 1L], length.out = os[i + 2L] - os[i + 1L])
    fixcol <- -1L
    for (j in range) {
      if (domain$col_lower_[oc[j + 1L] + 1L] == domain$col_upper_[oc[j + 1L] + 1L]) {
        fixcol <- oc[j + 1L]
        break
      }
    }
    if (fixcol != -1L) {
      oldNumFixed <- numFixed
      oldSize <- length(domain$domchgstack_column)
      if (domain$col_lower_[fixcol + 1L] == 1) {
        for (j in range) {
          if (domain$col_lower_[oc[j + 1L] + 1L] == 1) next
          numFixed <- numFixed + 1L
          .hmip_HighsDomain_changeBound3(domain, 0L, oc[j + 1L], 1, unspecified)
          if (domain$infeasible_) return(numFixed)
        }
      } else {
        for (j in range) {
          if (domain$col_upper_[oc[j + 1L] + 1L] == 0) next
          numFixed <- numFixed + 1L
          .hmip_HighsDomain_changeBound3(domain, 1L, oc[j + 1L], 0, unspecified)
          if (domain$infeasible_) return(numFixed)
        }
      }
      newFixed <- numFixed - oldNumFixed
      if (newFixed != 0L) {
        .hmip_HighsDomain_propagate(domain)
        if (domain$infeasible_) return(numFixed)
        if (length(domain$domchgstack_column) - oldSize > newFixed) i <- -1L
      }
    }
    i <- i + 1L
  }
  numFixed
}

# ---------------------------------------------------------------------------
# HighsOrbitopeMatrix: matrix holds column ids, entry(i, j) at
# matrix[i + j * numRows + 1]; columnToRow is a NA-filled vector

.hmip_new_HighsOrbitopeMatrix <- function() {
  om <- new.env(parent = emptyenv())
  om$rowLength <- 0L
  om$numRows <- 0L
  om$numSetPackingRows <- 0L
  om$columnToRow <- integer(0)
  om$rowIsSetPacking <- integer(0)
  om$matrix <- integer(0)
  om
}

.hmip_HighsOrbitopeMatrix_entry <- function(om, i, j) om$matrix[i + j * om$numRows + 1L]

# HighsCliqueTable::findCommonClique({c0, val}, {c1, val}) is expected to
# return the CliqueVar indices (2 * col + val) of the common clique, or NULL
.hmip_HighsOrbitopeMatrix_packingPass <- function(om, cliquetable, val, mark) {
  nr <- om$numRows
  rl <- om$rowLength
  for (j in seq_len(rl - 1L)) {
    for (j0 in 0:(j - 1L)) {
      for (i in 0:(nr - 1L)) {
        if (om$rowIsSetPacking[i + 1L] != -1L) next
        x0 <- .hmip_HighsOrbitopeMatrix_entry(om, i, j0)
        x1 <- .hmip_HighsOrbitopeMatrix_entry(om, i, j)
        clq <- .hmip_HighsCliqueTable_findCommonClique(cliquetable, 2L * x0 + val, 2L * x1 + val)
        if (is.null(clq)) {
          om$rowIsSetPacking[i + 1L] <- 0L
          next
        }
        overlap <- 0L
        for (cv in clq) {
          if (cv %% 2L != val) next
          r <- .hmip_sym_find(om$columnToRow, cv %/% 2L)
          if (!is.na(r) && r == i) overlap <- overlap + 1L
        }
        if (overlap == rl) {
          om$rowIsSetPacking[i + 1L] <- mark
          om$numSetPackingRows <- om$numSetPackingRows + 1L
          if (om$numSetPackingRows == nr) break
        }
      }
      if (om$numSetPackingRows == nr) break
    }
    if (om$numSetPackingRows == nr) break
  }
}

.hmip_HighsOrbitopeMatrix_determineOrbitopeType <- function(om, cliquetable) {
  for (j in seq_len(om$rowLength) - 1L) {
    for (i in seq_len(om$numRows) - 1L) {
      om$columnToRow <- .hmip_sym_insert(om$columnToRow, .hmip_HighsOrbitopeMatrix_entry(om, i, j), i)
    }
  }
  om$rowIsSetPacking <- rep(-1L, om$numRows)
  om$numSetPackingRows <- 0L
  .hmip_HighsOrbitopeMatrix_packingPass(om, cliquetable, 1L, 1L)
  om$rowIsSetPacking[om$rowIsSetPacking == 0L] <- -1L
  .hmip_HighsOrbitopeMatrix_packingPass(om, cliquetable, 0L, 2L)
  invisible(om)
}

.hmip_HighsOrbitopeMatrix_getBranchingColumn <- function(om, colLower, colUpper, col) {
  i <- .hmip_sym_find(om$columnToRow, col)
  if (!is.na(i) && om$rowIsSetPacking[i + 1L] != 0L) {
    for (j in seq_len(om$rowLength) - 1L) {
      b <- .hmip_HighsOrbitopeMatrix_entry(om, i, j)
      if (b == col) break
      if (colLower[b + 1L] != colUpper[b + 1L]) return(b)
    }
  }
  col
}

.hmip_HighsOrbitopeMatrix_orbitalFixing <- function(om, domain) {
  rowUsed <- logical(om$numRows)
  rows <- integer(0)
  isPacking <- TRUE
  for (pos in domain$branchPos_) {
    i <- .hmip_sym_find(om$columnToRow, domain$domchgstack_column[pos + 1L])
    if (!is.na(i) && !rowUsed[i + 1L]) {
      rowUsed[i + 1L] <- TRUE
      isPacking <- isPacking && om$rowIsSetPacking[i + 1L] != 0L
      rows <- c(rows, i)
    }
  }
  if (length(rows) == 0L) return(0L)
  if (isPacking) return(.hmip_HighsOrbitopeMatrix_orbitalFixingForPackingOrbitope(om, rows, domain))
  .hmip_HighsOrbitopeMatrix_orbitalFixingForFullOrbitope(om, rows, domain)
}

.hmip_HighsOrbitopeMatrix_orbitalFixingForPackingOrbitope <- function(om, rows, domain) {
  nd <- length(rows)
  rl <- om$rowLength
  E <- function(r, j) .hmip_HighsOrbitopeMatrix_entry(om, r, j)
  neg <- om$rowIsSetPacking[rows + 1L] == 2L
  unspecified <- list(type = -2L, index = 0L)
  notZeroFixed <- function(k, j) {
    col <- E(rows[k + 1L], j)
    if (neg[k + 1L]) domain$col_lower_[col + 1L] < 0.5 else domain$col_upper_[col + 1L] > 0.5
  }
  # fix the entry (k, j) to zero in the (possibly negated) row
  fixZero <- function(k, j) {
    col <- E(rows[k + 1L], j)
    if (neg[k + 1L]) {
      if (domain$col_lower_[col + 1L] > 0.5) return(FALSE)
      .hmip_HighsDomain_changeBound3(domain, 0L, col, 1, unspecified)
    } else {
      if (domain$col_upper_[col + 1L] < 0.5) return(FALSE)
      .hmip_HighsDomain_changeBound3(domain, 1L, col, 0, unspecified)
    }
    TRUE
  }
  firstOneInRow <- rep(-1L, nd)
  for (j in seq_len(rl) - 1L) {
    for (i in seq_len(nd) - 1L) {
      if (firstOneInRow[i + 1L] != -1L) next
      colrj <- E(rows[i + 1L], j)
      # C++ tests rowIsSetPacking == 1 here (a row left at -1 takes the
      # upper bound branch), but == 2 for the negation everywhere else
      if (om$rowIsSetPacking[rows[i + 1L] + 1L] == 1L) {
        if (domain$col_lower_[colrj + 1L] > 0.5) firstOneInRow[i + 1L] <- j
      } else {
        if (domain$col_upper_[colrj + 1L] < 0.5) firstOneInRow[i + 1L] <- j
      }
    }
  }
  j <- 0L
  numFixed <- 0L
  for (i in seq_len(nd) - 1L) {
    if (firstOneInRow[i + 1L] > j) {
      .hmip_HighsDomain_markInfeasible(domain, unspecified)
      return(numFixed)
    }
    col_ij <- E(rows[i + 1L], j)
    if (notZeroFixed(i, j)) {
      j0 <- j
      k <- i + 1L
      while (k < nd) {
        if (firstOneInRow[k + 1L] > j0) {
          if (neg[i + 1L]) {
            .hmip_HighsDomain_changeBound3(domain, 1L, col_ij, 0, unspecified)
          } else {
            .hmip_HighsDomain_changeBound3(domain, 0L, col_ij, 1, unspecified)
          }
          numFixed <- numFixed + 1L
          if (domain$infeasible_) return(numFixed)
          break
        }
        if (notZeroFixed(k, j0)) {
          j0 <- j0 + 1L
          if (j0 == rl) break
        }
        k <- k + 1L
      }
      j <- j + 1L
      if (j == rl) break
      for (k in 0:i) {
        if (!fixZero(k, j)) next
        numFixed <- numFixed + 1L
        if (domain$infeasible_) return(numFixed)
      }
    }
  }
  repeat {
    j <- j + 1L
    if (j >= rl) break
    for (k in seq_len(nd) - 1L) {
      if (!fixZero(k, j)) next
      numFixed <- numFixed + 1L
      if (domain$infeasible_) return(numFixed)
    }
  }
  if (!domain$infeasible_ && numFixed) .hmip_HighsDomain_propagate(domain)
  numFixed
}

.hmip_HighsOrbitopeMatrix_orbitalFixingForFullOrbitope <- function(om, rows, domain) {
  nd <- length(rows)
  rl <- om$rowLength
  nr <- om$numRows
  unspecified <- list(type = -2L, index = 0L)
  # Mminimal and Mmaximal as nd x rl matrices with entries -1, 0, 1
  Mmin <- matrix(-1L, nd, rl)
  for (j in seq_len(rl)) {
    for (i in seq_len(nd)) {
      colij <- om$matrix[rows[i] + (j - 1L) * nr + 1L]
      if (domain$col_lower_[colij + 1L] == 1) {
        Mmin[i, j] <- 1L
      } else if (domain$col_upper_[colij + 1L] == 0) {
        Mmin[i, j] <- 0L
      }
    }
  }
  Mmax <- Mmin
  Mmin[Mmin[, rl] == -1L, rl] <- 0L
  Mmax[Mmax[, 1L] == -1L, 1L] <- 1L
  # returns the 0-based row or nd
  i_fixed <- function(c0, c1) {
    w <- which(c0 != -1L & c1 != -1L & c0 != c1)
    if (length(w)) w[1L] - 1L else nd
  }
  i_discr <- function(c0, c1, i_f) {
    i <- i_f
    while (i >= 0L) {
      # c0[i] read past the end at i == nd is never 0 or 1 in C++ either
      # because i_f < nd whenever this is called
      if (c0[i + 1L] != 0L && c1[i + 1L] != 1L) return(i)
      i <- i - 1L
    }
    -1L
  }
  # colj0[k] += (isFree & colj1[k]) + isFree on rows k in ks
  addFree <- function(c0, c1, ks) {
    for (k in ks) {
      isFree <- as.integer(c0[k] == -1L)
      c0[k] <- c0[k] + bitwAnd(isFree, c1[k]) + isFree
    }
    c0
  }
  if (rl >= 2L) for (j in (rl - 1L):1L) {
    c0 <- Mmin[, j]
    c1 <- Mmin[, j + 1L]
    i_f <- i_fixed(c0, c1)
    if (i_f == nd) {
      c0 <- addFree(c0, c1, seq_len(nd))
    } else {
      i_d <- i_discr(c0, c1, i_f)
      if (i_d == -1L) {
        .hmip_HighsDomain_markInfeasible(domain, unspecified)
        return(0L)
      }
      c0 <- addFree(c0, c1, seq_len(i_d))
      c0[i_d + 1L] <- 1L
      ks <- seq.int(i_d + 2L, length.out = nd - i_d - 1L)
      c0[ks] <- c0[ks] + as.integer(c0[ks] == -1L)
    }
    Mmin[, j] <- c0
  }
  if (rl >= 2L) for (j in 2:rl) {
    c0 <- Mmax[, j]
    c1 <- Mmax[, j - 1L]
    i_f <- i_fixed(c1, c0)
    if (i_f == nd) {
      c0 <- addFree(c0, c1, seq_len(nd))
    } else {
      i_d <- i_discr(c1, c0, i_f)
      if (i_d == -1L) {
        .hmip_HighsDomain_markInfeasible(domain, unspecified)
        return(0L)
      }
      c0 <- addFree(c0, c1, seq_len(i_d))
      c0[i_d + 1L] <- 0L
      ks <- seq.int(i_d + 2L, length.out = nd - i_d - 1L)
      c0[ks] <- c0[ks] + 2L * as.integer(c0[ks] == -1L)
    }
    Mmax[, j] <- c0
  }
  numFixed <- 0L
  for (j in seq_len(rl)) {
    for (i in seq_len(nd)) {
      if (Mmin[i, j] != Mmax[i, j]) break
      colrj <- om$matrix[rows[i] + (j - 1L) * nr + 1L]
      if (domain$col_lower_[colrj + 1L] == domain$col_upper_[colrj + 1L]) next
      numFixed <- numFixed + 1L
      if (Mmin[i, j] == 1L) {
        .hmip_HighsDomain_changeBound3(domain, 0L, colrj, 1, unspecified)
      } else {
        .hmip_HighsDomain_changeBound3(domain, 1L, colrj, 0, unspecified)
      }
      if (domain$infeasible_) break
    }
    if (domain$infeasible_) break
  }
  if (!domain$infeasible_) .hmip_HighsDomain_propagate(domain)
  numFixed
}

# ---------------------------------------------------------------------------
# HighsSymmetryDetection

.hmip_new_HighsSymmetryDetection <- function() {
  sd <- new.env(parent = emptyenv())
  sd$firstLeavePartition <- integer(0)
  sd$bestLeavePartition <- integer(0)
  sd$firstLeaveCertificate <- numeric(0)
  sd$bestLeaveCertificate <- numeric(0)
  sd$firstLeaveGraph <- character(0)
  sd$bestLeaveGraph <- character(0)
  sd$firstLeavePrefixLen <- 0L
  sd$bestLeavePrefixLen <- 0L
  sd$firstPathDepth <- 0L
  sd$bestPathDepth <- 0L
  sd$linkCompressionStack <- integer(0)
  sd$refinementQueue <- integer(0)
  sd$cellCreationStack <- integer(0)
  sd$currNodeCertificate <- numeric(0)
  sd$numAutomorphisms <- 0L
  sd$termcache <- new.env(parent = emptyenv())
  sd
}

.hmip_sd_cellSize <- function(sd, cell) sd$currentPartitionLinks[cell + 1L] - cell

# vertexHash[vertex] via operator[], then sparse_combine32 with (cell, color)
.hmip_sd_combine <- function(sd, vertex, cell, color) {
  key <- paste(cell, color)
  term <- sd$termcache[[key]]
  if (is.null(term)) {
    term <- .hmip_combine32_term(cell, color)
    assign(key, term, envir = sd$termcache)
  }
  sd$vhset[vertex + 1L] <- TRUE
  sd$vh[vertex + 1L] <- .hmip_sparse_combine32(sd$vh[vertex + 1L], term)
}

.hmip_sd_clearHash <- function(sd) {
  sd$vh[] <- 0
  sd$vhset[] <- FALSE
}

.hmip_HighsSymmetryDetection_loadModelAsGraph <- function(sd, model, epsilon) {
  sd$model <- model
  numCol <- as.integer(model$num_col_)
  numRow <- as.integer(model$num_row_)
  nv <- numRow + numCol
  sd$numCol <- numCol
  sd$numRow <- numRow
  sd$numVertices <- nv
  sd$cellInRefinementQueue <- logical(nv)
  vtc <- integer(nv)
  coloring <- .hmip_new_HighsMatrixColoring(epsilon)
  am <- model$a_matrix_
  idx <- as.integer(am$index_)
  val <- am$value_
  numNz <- length(idx)
  Gfirst <- integer(2L * numNz)
  Gsecond <- numeric(2L * numNz)
  Gfirst[seq_len(numNz)] <- idx + numCol
  Gstart <- integer(nv + 1L)
  Gstart[seq_len(numCol + 1L)] <- as.integer(am$start_)
  rowSizes <- integer(numRow)
  for (j in seq_len(Gstart[numCol + 1L])) {
    Gsecond[j] <- .hmip_HighsMatrixColoring_color(coloring, val[j])
    rowSizes[idx[j] + 1L] <- rowSizes[idx[j] + 1L] + 1L
  }
  offset <- numNz
  for (i in seq_len(numRow)) {
    Gstart[numCol + i] <- offset
    offset <- offset + rowSizes[i]
  }
  Gstart[nv + 1L] <- offset
  sd$Gend <- Gstart[-1L]
  for (i in seq_len(numCol) - 1L) {
    for (j in seq.int(Gstart[i + 1L], length.out = Gstart[i + 2L] - Gstart[i + 1L])) {
      row <- idx[j + 1L]
      ARpos <- Gstart[numCol + row + 2L] - rowSizes[row + 1L]
      rowSizes[row + 1L] <- rowSizes[row + 1L] - 1L
      Gfirst[ARpos + 1L] <- i
      Gsecond[ARpos + 1L] <- Gsecond[j + 1L]
    }
  }
  # column and row classes (HighsHashTable<MatrixColumn/MatrixRow, HighsInt>)
  columnSet <- new.env(parent = emptyenv())
  nset <- 0L
  indexOffset <- numCol + 1L
  for (i in seq_len(numCol)) {
    key <- paste(.hmip_HighsMatrixColoring_color(coloring, model$col_cost_[i]),
                 .hmip_HighsMatrixColoring_color(coloring, model$col_lower_[i]),
                 .hmip_HighsMatrixColoring_color(coloring, model$col_upper_[i]),
                 model$integrality_[i], Gstart[i + 1L] - Gstart[i])
    cell <- columnSet[[key]]
    if (is.null(cell)) {
      nset <- nset + 1L
      cell <- nset
      if (model$col_lower_[i] != 0 || model$col_upper_[i] != 1 || model$integrality_[i] == 0L) {
        cell <- cell + indexOffset
      }
      assign(key, cell, envir = columnSet)
    }
    vtc[i] <- cell
  }
  rowSet <- new.env(parent = emptyenv())
  nset <- 0L
  indexOffset <- 2L * numCol + 1L
  for (i in seq_len(numRow)) {
    key <- paste(.hmip_HighsMatrixColoring_color(coloring, model$row_lower_[i]),
                 .hmip_HighsMatrixColoring_color(coloring, model$row_upper_[i]),
                 Gstart[numCol + i + 1L] - Gstart[numCol + i])
    cell <- rowSet[[key]]
    if (is.null(cell)) {
      nset <- nset + 1L
      cell <- nset
      assign(key, cell, envir = rowSet)
    }
    vtc[numCol + i] <- indexOffset + cell
  }
  cp <- .hmip_pdqsort(seq_len(nv) - 1L, function(v1, v2) vtc[v1 + 1L] < vtc[v2 + 1L])
  links <- integer(nv)
  cellStart <- 0L
  cellNumber <- 0L
  for (i in seq_len(nv) - 1L) {
    vertex <- cp[i + 1L]
    if (cellNumber != vtc[vertex + 1L]) {
      cellNumber <- vtc[vertex + 1L]
      links[cellStart + 1L] <- i
      cellStart <- i
    }
    vtc[vertex + 1L] <- cellStart
    links[i + 1L] <- cellStart
  }
  if (nv > 0L) links[cellStart + 1L] <- nv
  sd$Gstart <- Gstart
  sd$Gfirst <- Gfirst
  sd$Gsecond <- Gsecond
  sd$vertexToCell <- vtc
  sd$currentPartition <- cp
  sd$currentPartitionLinks <- links
  sd$vh <- numeric(nv)
  sd$vhset <- logical(nv)
  invisible(sd)
}

.hmip_HighsSymmetryDetection_markCellForRefinement <- function(sd, cell) {
  if (.hmip_sd_cellSize(sd, cell) == 1L || sd$cellInRefinementQueue[cell + 1L]) return(invisible(NULL))
  sd$cellInRefinementQueue[cell + 1L] <- TRUE
  # a min heap of distinct cells: only the pop order matters
  sd$refinementQueue <- c(sd$refinementQueue, cell)
  invisible(NULL)
}

.hmip_HighsSymmetryDetection_initializeHashValues <- function(sd) {
  for (i in seq_len(sd$numVertices) - 1L) {
    cell <- sd$vertexToCell[i + 1L]
    for (j in seq.int(sd$Gstart[i + 1L], length.out = sd$Gend[i + 1L] - sd$Gstart[i + 1L])) {
      .hmip_sd_combine(sd, sd$Gfirst[j + 1L], cell, sd$Gsecond[j + 1L])
    }
    .hmip_HighsSymmetryDetection_markCellForRefinement(sd, cell)
  }
}

.hmip_HighsSymmetryDetection_updateCellMembership <- function(sd, i, cell, markForRefinement = TRUE) {
  vertex <- sd$currentPartition[i + 1L]
  if (sd$vertexToCell[vertex + 1L] != cell) {
    sd$vertexToCell[vertex + 1L] <- cell
    if (i != cell) sd$currentPartitionLinks[i + 1L] <- cell
    if (markForRefinement) {
      for (j in seq.int(sd$Gstart[vertex + 1L], length.out = sd$Gend[vertex + 1L] - sd$Gstart[vertex + 1L])) {
        dest <- sd$Gfirst[j + 1L]
        dc <- sd$vertexToCell[dest + 1L]
        if (.hmip_sd_cellSize(sd, dc) == 1L) next
        .hmip_sd_combine(sd, dest, cell, sd$Gsecond[j + 1L])
        .hmip_HighsSymmetryDetection_markCellForRefinement(sd, dc)
      }
    }
    return(TRUE)
  }
  FALSE
}

.hmip_HighsSymmetryDetection_getVertexHash <- function(sd, v) if (sd$vhset[v + 1L]) sd$vh[v + 1L] else 0

.hmip_HighsSymmetryDetection_splitCell <- function(sd, cell, splitPoint) {
  cp <- sd$currentPartition
  links <- sd$currentPartitionLinks
  hSplit <- .hmip_HighsSymmetryDetection_getVertexHash(sd, cp[splitPoint + 1L])
  hCell <- .hmip_HighsSymmetryDetection_getVertexHash(sd, cp[cell + 1L])
  s <- .hs_u64_add(.hs_u64_add(.hs_pair_hash(0L, hSplit, hCell),
                               .hs_pair_hash(1L, cell, links[cell + 1L] - splitPoint)),
                   .hs_pair_hash(2L, splitPoint, splitPoint - cell))
  cv <- .hs_hi32(s)
  if (length(sd$firstLeaveCertificate)) {
    n <- length(sd$currNodeCertificate)
    if (sd$firstLeavePrefixLen == n && cv == sd$firstLeaveCertificate[n + 1L]) {
      sd$firstLeavePrefixLen <- sd$firstLeavePrefixLen + 1L
    }
    if (sd$bestLeavePrefixLen == n && cv == sd$bestLeaveCertificate[n + 1L]) {
      sd$bestLeavePrefixLen <- sd$bestLeavePrefixLen + 1L
    }
    blp <- sd$bestLeavePrefixLen
    if (sd$firstLeavePrefixLen <= n && blp <= n) {
      diffVal <- if (blp == n) cv else sd$currNodeCertificate[blp + 1L]
      if (diffVal > sd$bestLeaveCertificate[blp + 1L]) return(FALSE)
    }
  }
  links[splitPoint + 1L] <- links[cell + 1L]
  links[cell + 1L] <- splitPoint
  sd$currentPartitionLinks <- links
  sd$cellCreationStack <- c(sd$cellCreationStack, splitPoint)
  sd$currNodeCertificate <- c(sd$currNodeCertificate, cv)
  TRUE
}

.hmip_sd_abortRefinement <- function(sd) {
  sd$cellInRefinementQueue[sd$refinementQueue + 1L] <- FALSE
  sd$refinementQueue <- integer(0)
  .hmip_sd_clearHash(sd)
}

.hmip_HighsSymmetryDetection_partitionRefinement <- function(sd) {
  upd <- .hmip_HighsSymmetryDetection_updateCellMembership
  while (length(sd$refinementQueue)) {
    k <- which.min(sd$refinementQueue)
    cellStart <- sd$refinementQueue[k]
    sd$refinementQueue <- sd$refinementQueue[-k]
    firstCellStart <- cellStart
    sd$cellInRefinementQueue[cellStart + 1L] <- FALSE
    if (.hmip_sd_cellSize(sd, cellStart) == 1L) next
    cellEnd <- sd$currentPartitionLinks[cellStart + 1L]
    rng <- (cellStart + 1L):cellEnd
    vhset <- sd$vhset
    part <- .hmip_std_partition(sd$currentPartition[rng], function(v) !vhset[v + 1L])
    sd$currentPartition[rng] <- part[[1L]]
    refineStart <- cellStart + part[[2L]]
    if (refineStart == cellEnd) next
    rng2 <- (refineStart + 1L):cellEnd
    vh <- sd$vh
    sd$currentPartition[rng2] <- .hmip_pdqsort(sd$currentPartition[rng2],
                                               function(v1, v2) vh[v1 + 1L] < vh[v2 + 1L])
    if (refineStart != cellStart) {
      if (!.hmip_HighsSymmetryDetection_splitCell(sd, cellStart, refineStart)) {
        .hmip_sd_abortRefinement(sd)
        return(FALSE)
      }
      cellStart <- refineStart
      upd(sd, cellStart, cellStart)
    }
    prune <- FALSE
    lastHash <- sd$vh[sd$currentPartition[cellStart + 1L] + 1L]
    i <- cellStart + 1L
    while (i < cellEnd) {
      vertex <- sd$currentPartition[i + 1L]
      hash <- sd$vh[vertex + 1L]
      if (hash != lastHash) {
        if (!.hmip_HighsSymmetryDetection_splitCell(sd, cellStart, i)) {
          prune <- TRUE
          break
        }
        cellStart <- i
        lastHash <- hash
      }
      upd(sd, i, cellStart)
      i <- i + 1L
    }
    if (prune) {
      .hmip_sd_abortRefinement(sd)
      sd$currentPartitionLinks[firstCellStart + 1L] <- cellEnd
      i <- i - 1L
      while (i >= refineStart) {
        upd(sd, i, firstCellStart, FALSE)
        i <- i - 1L
      }
      return(FALSE)
    }
  }
  .hmip_sd_clearHash(sd)
  TRUE
}

.hmip_HighsSymmetryDetection_removeFixPoints <- function(sd) {
  nv <- sd$numVertices
  Gend <- integer(nv)
  for (i in seq_len(nv)) {
    rng <- seq.int(sd$Gstart[i] + 1L, length.out = sd$Gstart[i + 1L] - sd$Gstart[i])
    if (length(rng)) {
      ord <- .hmip_std_partition(seq_along(rng), function(e) {
        .hmip_sd_cellSize(sd, sd$vertexToCell[sd$Gfirst[rng[e]] + 1L]) > 1L
      })
      sd$Gfirst[rng] <- sd$Gfirst[rng][ord[[1L]]]
      sd$Gsecond[rng] <- sd$Gsecond[rng][ord[[1L]]]
      Gend[i] <- sd$Gstart[i] + ord[[2L]]
    } else {
      Gend[i] <- sd$Gstart[i]
    }
  }
  sd$Gend <- Gend
  unitCellIndex <- nv
  keep <- logical(length(sd$currentPartition))
  for (k in seq_along(sd$currentPartition)) {
    v <- sd$currentPartition[k]
    if (.hmip_sd_cellSize(sd, sd$vertexToCell[v + 1L]) == 1L) {
      unitCellIndex <- unitCellIndex - 1L
      sd$vertexToCell[v + 1L] <- unitCellIndex
    } else {
      keep[k] <- TRUE
    }
  }
  sd$currentPartition <- sd$currentPartition[keep]
  for (i in seq_len(nv)) {
    if (Gend[i] == sd$Gstart[i + 1L]) next
    rng <- (Gend[i] + 1L):sd$Gstart[i + 1L]
    sd$Gfirst[rng] <- sd$vertexToCell[sd$Gfirst[rng] + 1L]
  }
  if (length(sd$currentPartition) < nv) {
    nv <- length(sd$currentPartition)
    sd$numVertices <- nv
    if (nv == 0L) {
      sd$numActiveCols <- 0L
      return(invisible(NULL))
    }
    sd$currentPartitionLinks <- sd$currentPartitionLinks[seq_len(nv)]
    sd$cellInRefinementQueue <- logical(nv)
    sd$refinementQueue <- integer(0)
    cellStart <- 0L
    cellNumber <- 0L
    for (i in seq_len(nv) - 1L) {
      vertex <- sd$currentPartition[i + 1L]
      if (cellNumber != sd$vertexToCell[vertex + 1L]) {
        cellNumber <- sd$vertexToCell[vertex + 1L]
        sd$currentPartitionLinks[cellStart + 1L] <- i
        cellStart <- i
      }
      .hmip_HighsSymmetryDetection_updateCellMembership(sd, i, cellStart, FALSE)
    }
    sd$currentPartitionLinks[cellStart + 1L] <- nv
    # std::partition_point on a range partitioned by v < numCol
    lo <- 0L
    len <- nv
    while (len > 0L) {
      half <- len %/% 2L
      if (sd$currentPartition[lo + half + 1L] < sd$numCol) {
        lo <- lo + half + 1L
        len <- len - half - 1L
      } else {
        len <- half
      }
    }
    sd$numActiveCols <- lo
  } else {
    sd$numActiveCols <- sd$numCol
  }
  invisible(NULL)
}

.hmip_HighsSymmetryDetection_initializeDetection <- function(sd) {
  .hmip_HighsSymmetryDetection_initializeHashValues(sd)
  .hmip_HighsSymmetryDetection_partitionRefinement(sd)
  .hmip_HighsSymmetryDetection_removeFixPoints(sd)
  sd$numActiveCols != 0L
}

.hmip_HighsSymmetryDetection_initializeGroundSet <- function(sd) {
  nv <- sd$numVertices
  sd$vertexGroundSet <- sort(sd$currentPartition)
  sd$vertexPosition <- rep(-1L, length(sd$vertexToCell))
  sd$vertexPosition[sd$vertexGroundSet + 1L] <- seq_len(nv) - 1L
  sd$orbitPartition <- seq_len(nv) - 1L
  sd$orbitSize <- rep(1L, nv)
  sd$automorphisms <- integer(nv * 64L)
  sd$numAutomorphisms <- 0L
}

.hmip_HighsSymmetryDetection_getOrbit <- function(sd, vertex) {
  i <- sd$vertexPosition[vertex + 1L]
  op <- sd$orbitPartition
  orbit <- op[i + 1L]
  if (orbit != op[orbit + 1L]) {
    stack <- integer(0)
    repeat {
      stack <- c(stack, i)
      i <- orbit
      orbit <- op[orbit + 1L]
      if (orbit == op[orbit + 1L]) break
    }
    op[stack + 1L] <- orbit
    sd$orbitPartition <- op
  }
  orbit
}

.hmip_HighsSymmetryDetection_mergeOrbits <- function(sd, v1, v2) {
  if (v1 == v2) return(FALSE)
  o1 <- .hmip_HighsSymmetryDetection_getOrbit(sd, v1)
  o2 <- .hmip_HighsSymmetryDetection_getOrbit(sd, v2)
  if (o1 == o2) return(FALSE)
  if (o1 < o2) {
    sd$orbitPartition[o2 + 1L] <- o1
    sd$orbitSize[o1 + 1L] <- sd$orbitSize[o1 + 1L] + sd$orbitSize[o2 + 1L]
  } else {
    sd$orbitPartition[o1 + 1L] <- o2
    sd$orbitSize[o2 + 1L] <- sd$orbitSize[o2 + 1L] + sd$orbitSize[o1 + 1L]
  }
  TRUE
}

# node stack as four parallel vectors
.hmip_HighsSymmetryDetection_createNode <- function(sd) {
  sd$ns_stackStart <- c(sd$ns_stackStart, length(sd$cellCreationStack))
  sd$ns_certificateEnd <- c(sd$ns_certificateEnd, length(sd$currNodeCertificate))
  sd$ns_targetCell <- c(sd$ns_targetCell, -1L)
  sd$ns_lastDistiguished <- c(sd$ns_lastDistiguished, -1L)
}

.hmip_sd_resizeNodes <- function(sd, n) {
  k <- seq_len(n)
  sd$ns_stackStart <- sd$ns_stackStart[k]
  sd$ns_certificateEnd <- sd$ns_certificateEnd[k]
  sd$ns_targetCell <- sd$ns_targetCell[k]
  sd$ns_lastDistiguished <- sd$ns_lastDistiguished[k]
}

.hmip_HighsSymmetryDetection_selectTargetCell <- function(sd) {
  i <- 0L
  n <- length(sd$ns_targetCell)
  if (n > 1L) i <- sd$ns_targetCell[n - 1L]
  while (i < sd$numVertices) {
    if (.hmip_sd_cellSize(sd, i) > 1L) return(i)
    i <- i + 1L
  }
  -1L
}

.hmip_HighsSymmetryDetection_checkStoredAutomorphism <- function(sd, vertex) {
  nv <- sd$numVertices
  numCheck <- min(sd$numAutomorphisms, 64L)
  n <- length(sd$ns_targetCell)
  vp <- sd$vertexPosition
  for (i in seq_len(numCheck) - 1L) {
    aut <- sd$automorphisms[i * nv + seq_len(nv)]
    useful <- TRUE
    j <- n - 2L
    while (j >= sd$firstPathDepth) {
      fixPos <- vp[sd$ns_lastDistiguished[j + 1L] + 1L]
      if (aut[fixPos + 1L] != sd$vertexGroundSet[fixPos + 1L]) {
        useful <- FALSE
        break
      }
      j <- j - 1L
    }
    if (!useful) next
    if (aut[vp[vertex + 1L] + 1L] < vertex) return(FALSE)
  }
  TRUE
}

# sets sd$distinguishPos (0-based position in currentPartition)
.hmip_HighsSymmetryDetection_determineNextToDistinguish <- function(sd) {
  n <- length(sd$ns_targetCell)
  tc <- sd$ns_targetCell[n]
  last <- sd$ns_lastDistiguished[n]
  pos <- tc:(sd$currentPartitionLinks[tc + 1L] - 1L)
  vals <- sd$currentPartition[pos + 1L]
  if (last == -1L) {
    sd$distinguishPos <- pos[which.min(vals)]
    return(TRUE)
  }
  ok <- logical(length(vals))
  for (k in seq_along(vals)) {
    v <- vals[k]
    if (v <= last) next
    ok[k] <- if (n > sd$firstPathDepth) {
      .hmip_HighsSymmetryDetection_checkStoredAutomorphism(sd, v)
    } else {
      sd$vertexGroundSet[.hmip_HighsSymmetryDetection_getOrbit(sd, v) + 1L] == v
    }
  }
  if (!any(ok)) return(FALSE)
  sd$distinguishPos <- pos[ok][which.min(vals[ok])]
  TRUE
}

.hmip_HighsSymmetryDetection_distinguishVertex <- function(sd, targetCell) {
  newCell <- sd$currentPartitionLinks[targetCell + 1L] - 1L
  p <- sd$distinguishPos
  tmp <- sd$currentPartition[p + 1L]
  sd$currentPartition[p + 1L] <- sd$currentPartition[newCell + 1L]
  sd$currentPartition[newCell + 1L] <- tmp
  sd$ns_lastDistiguished[length(sd$ns_lastDistiguished)] <- sd$currentPartition[newCell + 1L]
  if (!.hmip_HighsSymmetryDetection_splitCell(sd, targetCell, newCell)) return(FALSE)
  .hmip_HighsSymmetryDetection_updateCellMembership(sd, newCell, newCell)
  TRUE
}

.hmip_HighsSymmetryDetection_getCellStart <- function(sd, pos) {
  links <- sd$currentPartitionLinks
  startPos <- links[pos + 1L]
  if (startPos > pos) return(pos)
  if (links[startPos + 1L] < startPos) {
    stack <- integer(0)
    repeat {
      stack <- c(stack, pos)
      pos <- startPos
      startPos <- links[startPos + 1L]
      if (!(links[startPos + 1L] < startPos)) break
    }
    links[stack + 1L] <- startPos
    sd$currentPartitionLinks <- links
  }
  startPos
}

.hmip_HighsSymmetryDetection_backtrack <- function(sd, newEnd, end) {
  stackPos <- end - 1L
  while (stackPos >= newEnd) {
    cell <- sd$cellCreationStack[stackPos + 1L]
    newStart <- .hmip_HighsSymmetryDetection_getCellStart(sd, cell - 1L)
    currEnd <- sd$currentPartitionLinks[cell + 1L]
    sd$currentPartitionLinks[cell + 1L] <- newStart
    sd$currentPartitionLinks[newStart + 1L] <- currEnd
    stackPos <- stackPos - 1L
  }
}

.hmip_HighsSymmetryDetection_cleanupBacktrack <- function(sd, pos) {
  stackPos <- length(sd$cellCreationStack) - 1L
  while (stackPos >= pos) {
    cell <- sd$cellCreationStack[stackPos + 1L]
    cs <- .hmip_HighsSymmetryDetection_getCellStart(sd, cell)
    ce <- sd$currentPartitionLinks[cs + 1L]
    v <- cell
    while (v < ce && sd$vertexToCell[sd$currentPartition[v + 1L] + 1L] == cell) {
      .hmip_HighsSymmetryDetection_updateCellMembership(sd, v, cs, FALSE)
      v <- v + 1L
    }
    stackPos <- stackPos - 1L
  }
  sd$cellCreationStack <- sd$cellCreationStack[seq_len(pos)]
}

.hmip_HighsSymmetryDetection_switchToNextNode <- function(sd, backtrackDepth) {
  stackEnd <- length(sd$cellCreationStack)
  .hmip_sd_resizeNodes(sd, backtrackDepth)
  if (backtrackDepth == 0L) return(invisible(NULL))
  repeat {
    n <- length(sd$ns_targetCell)
    .hmip_HighsSymmetryDetection_backtrack(sd, sd$ns_stackStart[n], stackEnd)
    stackEnd <- sd$ns_stackStart[n]
    sd$firstPathDepth <- min(n, sd$firstPathDepth)
    sd$bestPathDepth <- min(n, sd$bestPathDepth)
    sd$firstLeavePrefixLen <- min(sd$ns_certificateEnd[n], sd$firstLeavePrefixLen)
    sd$bestLeavePrefixLen <- min(sd$ns_certificateEnd[n], sd$bestLeavePrefixLen)
    sd$currNodeCertificate <- sd$currNodeCertificate[seq_len(sd$ns_certificateEnd[n])]
    if (!.hmip_HighsSymmetryDetection_determineNextToDistinguish(sd)) {
      .hmip_sd_resizeNodes(sd, n - 1L)
      if (n == 1L) break
      next
    }
    .hmip_HighsSymmetryDetection_cleanupBacktrack(sd, stackEnd)
    if (!.hmip_HighsSymmetryDetection_distinguishVertex(sd, sd$ns_targetCell[n])) {
      .hmip_sd_resizeNodes(sd, n - 1L)
      if (n == 1L) break
      next
    }
    if (!.hmip_HighsSymmetryDetection_partitionRefinement(sd)) {
      stackEnd <- length(sd$cellCreationStack)
      next
    }
    .hmip_HighsSymmetryDetection_createNode(sd)
    break
  }
  invisible(NULL)
}

# the HighsHashTable<tuple<HighsInt, HighsInt, HighsUInt>> of the current
# graph: one key per edge of every column, in column order, with the cell
# of the column (compareCurrentGraph reports the first missing one)
.hmip_sd_graphKeys <- function(sd) {
  gs <- sd$Gstart
  ge <- sd$Gend
  keys <- character(0)
  cells <- integer(0)
  for (i in seq_len(sd$numCol)) {
    colCell <- sd$vertexToCell[i]
    r1 <- seq.int(gs[i] + 1L, length.out = ge[i] - gs[i])
    r2 <- seq.int(ge[i] + 1L, length.out = gs[i + 1L] - ge[i])
    if (length(r1) + length(r2) == 0L) next
    a <- c(sd$vertexToCell[sd$Gfirst[r1] + 1L], sd$Gfirst[r2])
    keys <- c(keys, paste(a, colCell, c(sd$Gsecond[r1], sd$Gsecond[r2])))
    cells <- c(cells, rep(colCell, length(a)))
  }
  list(keys = keys, cells = cells)
}

.hmip_HighsSymmetryDetection_dumpCurrentGraph <- function(sd) unique(.hmip_sd_graphKeys(sd)$keys)

# returns list(equal, wrongCell)
.hmip_HighsSymmetryDetection_compareCurrentGraph <- function(sd, otherGraph) {
  g <- .hmip_sd_graphKeys(sd)
  miss <- which(!(g$keys %in% otherGraph))
  if (length(miss)) return(list(equal = FALSE, wrongCell = g$cells[miss[1L]]))
  list(equal = TRUE, wrongCell = -1L)
}

.hmip_HighsSymmetryDetection_isFromBinaryColumn <- function(sd, pos) {
  if (pos >= sd$numActiveCols) return(FALSE)
  col <- sd$currentPartition[pos + 1L]
  m <- sd$model
  !(m$col_lower_[col + 1L] != 0 || m$col_upper_[col + 1L] != 1 || m$integrality_[col + 1L] == 0L)
}

# store the automorphism mapping the current leaf to leafPartition and
# merge orbits; returns TRUE when it acts on an active column
.hmip_sd_storeAutomorphism <- function(sd, sym, leafPartition) {
  nv <- sd$numVertices
  k <- sd$numAutomorphisms %% 64L
  sd$numAutomorphisms <- sd$numAutomorphisms + 1L
  perm <- integer(nv)
  perm[sd$vertexPosition[sd$currentPartition + 1L] + 1L] <- leafPartition
  sd$automorphisms[k * nv + seq_len(nv)] <- perm
  report <- FALSE
  for (i in seq_len(nv)) {
    if (.hmip_HighsSymmetryDetection_mergeOrbits(sd, perm[i], sd$vertexGroundSet[i]) && i <= sd$numActiveCols) {
      report <- TRUE
    }
  }
  if (report) {
    sym$permutations <- c(sym$permutations, perm[seq_len(sd$numActiveCols)])
    sym$numPerms <- sym$numPerms + 1L
  }
  report
}

.hmip_HighsSymmetryDetection_run <- function(sd, sym) {
  .hmip_HighsSymmetryDetection_initializeGroundSet(sd)
  sd$currNodeCertificate <- numeric(0)
  sd$cellCreationStack <- integer(0)
  sd$ns_stackStart <- integer(0)
  sd$ns_certificateEnd <- integer(0)
  sd$ns_targetCell <- integer(0)
  sd$ns_lastDistiguished <- integer(0)
  .hmip_HighsSymmetryDetection_createNode(sd)
  maxPerms <- 64000000L %/% sd$numActiveCols
  while (length(sd$ns_targetCell)) {
    targetCell <- .hmip_HighsSymmetryDetection_selectTargetCell(sd)
    if (targetCell == -1L) {
      nc <- length(sd$currNodeCertificate)
      if (length(sd$firstLeavePartition) == 0L) {
        sd$firstLeavePartition <- sd$currentPartition
        sd$firstLeaveCertificate <- sd$currNodeCertificate
        sd$bestLeaveCertificate <- sd$currNodeCertificate
        sd$firstLeaveGraph <- .hmip_HighsSymmetryDetection_dumpCurrentGraph(sd)
        sd$firstPathDepth <- length(sd$ns_targetCell)
        sd$bestPathDepth <- length(sd$ns_targetCell)
        sd$firstLeavePrefixLen <- nc
        sd$bestLeavePrefixLen <- nc
        backtrackDepth <- sd$firstPathDepth - 1L
        while (backtrackDepth > 0L &&
               !.hmip_HighsSymmetryDetection_isFromBinaryColumn(sd, sd$ns_targetCell[backtrackDepth])) {
          backtrackDepth <- backtrackDepth - 1L
        }
        .hmip_HighsSymmetryDetection_switchToNextNode(sd, backtrackDepth)
      } else {
        wrongCell <- -1L
        backtrackDepth <- length(sd$ns_targetCell) - 1L
        if (sd$firstLeavePrefixLen == nc || sd$bestLeavePrefixLen == nc) {
          handled <- FALSE
          if (sd$firstLeavePrefixLen == nc) {
            cmp <- .hmip_HighsSymmetryDetection_compareCurrentGraph(sd, sd$firstLeaveGraph)
            wrongCell <- cmp$wrongCell
            if (cmp$equal) {
              handled <- TRUE
              if (.hmip_sd_storeAutomorphism(sd, sym, sd$firstLeavePartition) && sym$numPerms == maxPerms) break
              backtrackDepth <- min(backtrackDepth, sd$firstPathDepth)
            }
          }
          if (!handled && length(sd$bestLeavePartition) && sd$bestLeavePrefixLen == nc) {
            cmp <- .hmip_HighsSymmetryDetection_compareCurrentGraph(sd, sd$bestLeaveGraph)
            if (!cmp$equal) wrongCell <- cmp$wrongCell
            if (cmp$equal) {
              handled <- TRUE
              if (.hmip_sd_storeAutomorphism(sd, sym, sd$bestLeavePartition) && sym$numPerms == maxPerms) break
              backtrackDepth <- min(backtrackDepth, sd$bestPathDepth)
            }
          }
          if (!handled) {
            blp <- sd$bestLeavePrefixLen
            if (blp < nc && sd$currNodeCertificate[blp + 1L] > sd$bestLeaveCertificate[blp + 1L]) {
              pbd <- sd$firstPathDepth - 1L
              while (sd$ns_certificateEnd[pbd + 1L] <= blp) pbd <- pbd + 1L
              backtrackDepth <- min(pbd, backtrackDepth)
            } else {
              pbd <- backtrackDepth
              while (pbd >= 0L) {
                if (sd$ns_targetCell[pbd + 1L] == wrongCell) {
                  backtrackDepth <- pbd
                  break
                }
                pbd <- pbd - 1L
              }
            }
          }
        } else {
          sd$bestLeaveCertificate <- sd$currNodeCertificate
          sd$bestLeaveGraph <- .hmip_HighsSymmetryDetection_dumpCurrentGraph(sd)
          sd$bestLeavePartition <- sd$currentPartition
          sd$bestPathDepth <- length(sd$ns_targetCell)
          sd$bestLeavePrefixLen <- nc
        }
        .hmip_HighsSymmetryDetection_switchToNextNode(sd, backtrackDepth)
      }
    } else {
      sd$ns_targetCell[length(sd$ns_targetCell)] <- targetCell
      .hmip_HighsSymmetryDetection_determineNextToDistinguish(sd)
      if (!.hmip_HighsSymmetryDetection_distinguishVertex(sd, targetCell)) {
        .hmip_HighsSymmetryDetection_switchToNextNode(sd, length(sd$ns_targetCell) - 1L)
        next
      }
      if (!.hmip_HighsSymmetryDetection_partitionRefinement(sd)) {
        .hmip_HighsSymmetryDetection_switchToNextNode(sd, length(sd$ns_targetCell))
        next
      }
      .hmip_HighsSymmetryDetection_createNode(sd)
    }
  }
  .hmip_HighsSymmetryDetection_finishRun(sd, sym)
}

.hmip_HighsSymmetryDetection_finishRun <- function(sd, sym) {
  sym$numGenerators <- sym$numPerms
  if (sym$numPerms == 0L) return(invisible(sym))
  nac <- sd$numActiveCols
  sd$vertexPosition <- sd$vertexPosition[seq_len(sd$numCol)]
  cd <- .hmip_HighsSymmetryDetection_computeComponentData(sd, sym)
  for (i in seq_len(length(cd$componentStarts) - 1L) - 1L) {
    if (cd$componentStarts[i + 2L] - cd$componentStarts[i + 1L] == 1L) next
    .hmip_HighsSymmetryDetection_isFullOrbitope(sd, cd, i, sym)
  }
  vgs <- sd$vertexGroundSet
  deleted <- logical(sym$numPerms)
  for (p in seq_len(sym$numPerms) - 1L) {
    perm <- sym$permutations[p * nac + seq_len(nac)]
    for (i in seq_len(nac)) {
      if (perm[i] != vgs[i] && !is.na(.hmip_sym_find(sym$columnToOrbitope, perm[i]))) {
        deleted[p + 1L] <- TRUE
        break
      }
    }
  }
  numFixed <- 0L
  for (i in seq_len(nac)) {
    v <- vgs[i]
    if (sd$orbitSize[.hmip_HighsSymmetryDetection_getOrbit(sd, v) + 1L] == 1L ||
        !is.na(.hmip_sym_find(sym$columnToOrbitope, v))) {
      sd$vertexPosition[v + 1L] <- -1L
      vgs[i] <- -1L
      numFixed <- numFixed + 1L
    }
  }
  # only the active columns (the first nac ground set entries) are compressed
  vgs <- vgs[seq_len(nac)]
  if (numFixed != 0L) {
    keepCol <- vgs != -1L
    out <- integer(0)
    np <- sym$numPerms
    for (p in seq_len(np) - 1L) {
      if (!deleted[p + 1L]) {
        out <- c(out, sym$permutations[p * nac + seq_len(nac)][keepCol])
      } else {
        sym$numPerms <- sym$numPerms - 1L
      }
    }
    sym$permutations <- out
    vgs <- vgs[keepCol]
    sd$vertexPosition[vgs + 1L] <- seq_along(vgs) - 1L
    nac <- nac - numFixed
    sd$numActiveCols <- nac
  }
  sd$vertexGroundSet <- vgs[seq_len(nac)]
  sym$permutationColumns <- sd$vertexGroundSet
  sym$columnPosition <- sd$vertexPosition
  sym$permutations <- sym$permutations[seq_len(sym$numPerms * nac)]
  invisible(sym)
}

.hmip_HighsSymmetryDetection_computeComponentData <- function(sd, sym) {
  nac <- sd$numActiveCols
  vgs <- sd$vertexGroundSet
  vp <- sd$vertexPosition
  cd <- list()
  ds <- .hmip_new_HighsDisjointSets(nac)
  firstUnfixed <- rep(-1L, sym$numPerms)
  numUnfixed <- integer(sym$numPerms)
  for (i in seq_len(sym$numPerms)) {
    perm <- sym$permutations[(i - 1L) * nac + seq_len(nac)]
    for (j in seq_len(nac)) {
      if (perm[j] != vgs[j]) {
        pos <- vp[perm[j] + 1L]
        numUnfixed[i] <- numUnfixed[i] + 1L
        if (firstUnfixed[i] != -1L) {
          .hmip_HighsDisjointSets_merge(ds, firstUnfixed[i], pos)
        } else {
          firstUnfixed[i] <- pos
        }
      }
    }
  }
  getSet <- function(x) .hmip_HighsDisjointSets_getSet(ds, x)
  setSize <- function(s) .hmip_HighsDisjointSets_getSetSize(ds, s)
  compKey <- function(u) {
    s <- getSet(vp[u + 1L])
    c(as.numeric(setSize(s) == 1L), s)
  }
  componentSets <- .hmip_pdqsort(vgs[seq_len(nac)], function(u, v) {
    a <- compKey(u)
    b <- compKey(v)
    a[1L] < b[1L] || (a[1L] == b[1L] && a[2L] < b[2L])
  })
  componentStarts <- integer(0)
  componentNumber <- integer(0)
  componentNumOrbits <- integer(0)
  currentComponent <- -1L
  orbitsSeen <- integer(0)
  for (i in seq_len(nac)) {
    comp <- getSet(vp[componentSets[i] + 1L])
    if (setSize(comp) == 1L) break
    if (comp != currentComponent) {
      currentComponent <- comp
      componentStarts <- c(componentStarts, i - 1L)
      componentNumber <- c(componentNumber, comp)
      componentNumOrbits <- c(componentNumOrbits, 0L)
      orbitsSeen <- integer(0)
    }
    o <- .hmip_HighsSymmetryDetection_getOrbit(sd, componentSets[i])
    if (!(o %in% orbitsSeen)) {
      orbitsSeen <- c(orbitsSeen, o)
      componentNumOrbits[length(componentNumOrbits)] <- componentNumOrbits[length(componentNumOrbits)] + 1L
    }
  }
  permComponents <- which(firstUnfixed != -1L) - 1L
  permComponents <- .hmip_pdqsort(permComponents, function(i, j) {
    si <- getSet(firstUnfixed[i + 1L])
    sj <- getSet(firstUnfixed[j + 1L])
    si < sj || (si == sj && numUnfixed[i + 1L] < numUnfixed[j + 1L])
  })
  permComponentStarts <- integer(0)
  currentComponent <- -1L
  for (i in seq_along(permComponents)) {
    comp <- getSet(firstUnfixed[permComponents[i] + 1L])
    if (comp != currentComponent) {
      currentComponent <- comp
      permComponentStarts <- c(permComponentStarts, i - 1L)
    }
  }
  list(components = ds, componentStarts = c(componentStarts, nac), componentSets = componentSets,
       componentNumOrbits = componentNumOrbits, componentNumber = componentNumber,
       permComponentStarts = c(permComponentStarts, length(permComponents)),
       permComponents = permComponents, firstUnfixed = firstUnfixed, numUnfixed = numUnfixed)
}

.hmip_HighsSymmetryDetection_isFullOrbitope <- function(sd, cd, component, sym) {
  nac <- sd$numActiveCols
  vgs <- sd$vertexGroundSet
  vp <- sd$vertexPosition
  m <- sd$model
  cs <- cd$componentStarts
  pcs <- cd$permComponentStarts
  componentSize <- cs[component + 2L] - cs[component + 1L]
  if (componentSize == 1L) return(FALSE)
  for (i in seq.int(cs[component + 1L] + 1L, length.out = componentSize)) {
    col <- cd$componentSets[i]
    if (m$integrality_[col + 1L] == 0L || m$col_lower_[col + 1L] != 0 || m$col_upper_[col + 1L] != 1) return(FALSE)
  }
  p0 <- cd$permComponents[pcs[component + 1L] + 1L]
  if (cd$numUnfixed[p0 + 1L] %% 2L == 1L) return(FALSE)
  others <- seq.int(pcs[component + 1L] + 2L, length.out = pcs[component + 2L] - pcs[component + 1L] - 1L)
  for (k in others) {
    if (cd$numUnfixed[cd$permComponents[k] + 1L] != cd$numUnfixed[p0 + 1L]) return(FALSE)
  }
  numRows <- cd$numUnfixed[p0 + 1L] %/% 2L
  orbitSize <- componentSize %/% numRows
  if (orbitSize * numRows != componentSize) return(FALSE)
  if (pcs[component + 2L] - pcs[component + 1L] != orbitSize - 1L) return(FALSE)
  getPerm <- function(p) sym$permutations[p * nac + seq_len(nac)]
  mat <- rep(-1L, componentSize)
  perm <- getPerm(p0)
  colSet <- integer(0)
  mrow <- 0L
  for (j in seq_len(nac) - 1L) {
    jImagePos <- vp[perm[j + 1L] + 1L]
    if (jImagePos <= j) next
    if (mrow == numRows) return(FALSE)
    if (perm[jImagePos + 1L] != vgs[j + 1L]) return(FALSE)
    mat[mrow + 1L] <- vgs[j + 1L]
    mat[numRows + mrow + 1L] <- perm[j + 1L]
    if (vgs[j + 1L] %in% colSet) return(FALSE)
    colSet <- c(colSet, vgs[j + 1L])
    if (perm[j + 1L] %in% colSet) return(FALSE)
    colSet <- c(colSet, perm[j + 1L])
    mrow <- mrow + 1L
  }
  numColsAdded <- 2L
  triedLeftExtension <- FALSE
  while (numColsAdded < orbitSize) {
    if (length(colSet) != numColsAdded * numRows) return(FALSE)
    thisOff <- numColsAdded * numRows
    prevOff <- (numColsAdded - 1L) * numRows
    foundCand <- FALSE
    repeat {
      movePos <- vp[mat[prevOff + 1L] + 1L]
      perm <- NULL
      for (k in others) {
        perm <- getPerm(cd$permComponents[k])
        if (perm[movePos + 1L] != vgs[movePos + 1L] && !(perm[movePos + 1L] %in% colSet)) {
          foundCand <- TRUE
          break
        }
      }
      if (!foundCand && !triedLeftExtension) {
        prevOff <- 0L
        triedLeftExtension <- TRUE
        next
      }
      break
    }
    if (!foundCand) return(FALSE)
    for (j in seq_len(numRows) - 1L) {
      prev <- mat[prevOff + j + 1L]
      nextVertex <- perm[vp[prev + 1L] + 1L]
      mat[thisOff + j + 1L] <- nextVertex
      if (perm[vp[nextVertex + 1L] + 1L] != prev) return(FALSE)
      if (nextVertex %in% colSet) return(FALSE)
      colSet <- c(colSet, nextVertex)
    }
    numColsAdded <- numColsAdded + 1L
  }
  if (length(colSet) != componentSize) return(FALSE)
  idx <- length(sym$orbitopes)
  for (col in mat) sym$columnToOrbitope <- .hmip_sym_insert(sym$columnToOrbitope, col, idx)
  om <- .hmip_new_HighsOrbitopeMatrix()
  om$matrix <- mat
  om$numRows <- numRows
  om$rowLength <- orbitSize
  sym$orbitopes[[idx + 1L]] <- om
  TRUE
}

# HighsMipSolverData::startSymmetryDetection and finishSymmetryDetection
# (HighsMipSolverData.cpp:488-549) without the timing; called by group B1.
# Returns the detection data (sd with member symmetries) or NULL.
.hmip_symmetry_start <- function(mipdata) {
  sd <- .hmip_new_HighsSymmetryDetection()
  .hmip_HighsSymmetryDetection_loadModelAsGraph(sd, mipdata$presolvedModel,
                                                 mipdata$mipsolver$options_mip_$small_matrix_value)
  mipdata$detectSymmetries <- .hmip_HighsSymmetryDetection_initializeDetection(sd)
  if (!mipdata$detectSymmetries) return(NULL)
  sd$symmetries <- .hmip_new_HighsSymmetries()
  .hmip_HighsSymmetryDetection_run(sd, sd$symmetries)
  sd
}

# returns the log lines HiGHS prints (without timing)
.hmip_symmetry_finish <- function(mipdata, sd) {
  sym <- sd$symmetries
  mipdata$symmetries <- sym
  nc2o <- .hmip_HighsSymmetries_columnToOrbitope_size(sym)
  msg <- if (sym$numGenerators == 0L) {
    mipdata$detectSymmetries <- FALSE
    "No symmetry present"
  } else if (length(sym$orbitopes) == 0L) {
    sprintf("Found %d generator(s)", sym$numGenerators)
  } else if (sym$numPerms != 0L) {
    sprintf("Found %d generator(s) and %d full orbitope(s) acting on %d columns",
            sym$numPerms, length(sym$orbitopes), nc2o)
  } else {
    sprintf("Found %d full orbitope(s) acting on %d columns", length(sym$orbitopes), nc2o)
  }
  for (om in sym$orbitopes) .hmip_HighsOrbitopeMatrix_determineOrbitopeType(om, mipdata$cliquetable)
  if (sym$numPerms != 0L) mipdata$globalOrbits <- .hmip_HighsSymmetries_computeStabilizerOrbits(sym, mipdata$domain)
  msg
}
