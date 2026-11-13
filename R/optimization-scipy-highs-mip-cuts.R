# The cut separation of the HiGHS 1.12.0 MIP solver as SciPy 1.18.1 vendors
# it: mip/HighsDynamicRowMatrix, HighsCutPool, util/HighsSparseVectorSum,
# HighsLpAggregator, HighsTransformedLp, HighsCutGeneration, HighsSeparator,
# HighsTableauSeparator, HighsPathSeparator, HighsModkSeparator,
# HighsGFkSolve, HighsSeparation, and the separation entry points
# HighsCliqueTable::separateCliques and
# HighsImplications::separateImpliedBounds. MIT license notice in
# R/optimization-scipy-highs-util.R.
#
# Conventions (scratchpad API-B1.md): every index value stored in a vector is
# 0-based as in C++ (-1 keeps its C++ meaning) and read as v[i + 1L]; the
# HighsGFkSolve internals alone use 1-based positions with 0 for -1 so that
# the splay tree helpers of R/optimization-scipy-highs-presolve.R apply.
# HighsCDouble values are c(hi, lo). Sums are always loops in C++ order (R's
# sum() accumulates in long double and is never used on C++ sums).

# ---------------------------------------------------------------------------
# HighsCDouble operations the presolve port does not provide

# CD - CD, CD * CD, CD / CD, floor and ceil are .cd_sub_cd, .cd_mul_cd,
# .cd_div_cd, .cd_floor and .cd_ceil of R/optimization-scipy-highs-mip-util.R
# double - CD (friend operator-: -b + a)
# std::frexp: c(mantissa, exponent), mantissa in [0.5, 1) in magnitude. The
# log2 guess is corrected by exact power of two comparisons, and the mantissa
# is formed as (x / 2^(e - 1)) / 2, exact for subnormals and near DBL_MAX.
# Zero, infinities and NaN come back unchanged with exponent 0 as in glibc.
.hmip_frexp <- function(x) {
  if (is.na(x) || x == 0 || !is.finite(x)) return(c(x, 0))
  a <- abs(x)
  e <- floor(log2(a)) + 1
  while (a >= 2^e) e <- e + 1
  while (a < 2^(e - 1)) e <- e - 1
  c((x / 2^(e - 1)) / 2, e)
}
# std::exp2(std::ceil(std::log2(y))) for y >= 1 with a correctly rounded
# log2, as glibc gives near powers of two: with 2^k <= y < 2^(k + 1) the
# rounded log2 is k exactly when y = 2^k or when log2(y / 2^k) is at most half
# an ulp of k (ties go to the even k); otherwise the ceiling is k + 1.
.hmip_exp2_ceil_log2 <- function(y) {
  k <- .hmip_frexp(y)[2L] - 1
  m <- y / 2^k
  if (m == 1) return(2^k)
  t <- log1p(m - 1) / log(2)
  if (k >= 1 && t <= 2^(.hmip_frexp(k)[2L] - 54)) 2^k else 2^(k + 1)
}
# HighsIntegers::isIntegral
.hmip_cuts_isIntegral <- function(x, eps) {
  y <- abs(x - trunc(x))
  min(y, 1 - y) <= eps
}
# HighsRandom::integer() of the 32 bit build: pair_hash<0>(lo, hi) >> 33
.hmip_rand_int31 <- function(r) {
  .hs_rng_advance(r)
  .hs_u64_top(.hs_pair_hash(0L, .hs_lo32(r$state), .hs_hi32(r$state)), 31L)
}
# HighsRandom::integer(min, sup)
.hmip_rand_range <- function(r, lo, sup) lo + .hs_rng_integer(r, sup - lo)
# HighsHashHelpers::hash of a nonnegative int64 value (.hmip_hash of the
# shared helpers takes the low and high 32 bit words)
.hmip_cuts_hash_i64 <- function(v) .hmip_hash(v %% 4294967296, v %/% 4294967296)
# a > b for two hash values
.hmip_cuts_hash_gt <- function(a, b) .hmip_u64_cmp(a, b) > 0L

# HighsHashHelpers::double_hash_code
.hmip_double_hash_code <- function(val) {
  fr <- .hmip_frexp(val * 0x1.3c6ef372fe94fp-1) # 0.61803398874989484 (golden ratio minus 1)
  (fr[2L] %% 65536) * 65536 + (trunc(fr[1L] * 32768) %% 65536)
}

# the LP relaxation and solver data behind a HighsLpRelaxation
.hmip_cuts_mip <- function(lp) .hmip_HighsLpRelaxation_getMipSolver(lp)
.hmip_cuts_sol <- function(lp) .hmip_HighsLpRelaxation_getSolution(lp)
.hmip_cuts_lp <- function(lp) .hmip_HighsLpRelaxation_getLp(lp)
# fractional integers as a vector of 0-based columns (C++ vector of pairs
# held either as list(first, second) or as a list of pairs)
.hmip_cuts_frac_cols <- function(fi) {
  if (is.null(fi) || length(fi) == 0L) return(integer())
  if (!is.null(fi$first)) return(fi$first)
  vapply(fi, function(p) as.integer(p[[1L]]), 0L)
}

# ---------------------------------------------------------------------------
# HighsDynamicRowMatrix (mip/HighsDynamicRowMatrix.cpp)

.hmip_new_HighsDynamicRowMatrix <- function(ncols) {
  m <- new.env(parent = emptyenv())
  m$ARrange_first <- integer(); m$ARrange_second <- integer()
  m$ARindex_ <- integer(); m$ARvalue_ <- numeric(); m$ARrowindex_ <- integer()
  m$AnextPos_ <- integer(); m$AprevPos_ <- integer()
  m$AnextNeg_ <- integer(); m$AprevNeg_ <- integer()
  m$AheadPos_ <- rep(-1L, ncols); m$AheadNeg_ <- rep(-1L, ncols)
  m$colsLinked <- logical()
  # std::set<std::pair<len, start>>
  m$freespaces_len <- integer(); m$freespaces_start <- integer()
  m$deletedrows_ <- integer()
  m
}
.hmip_HighsDynamicRowMatrix_columnsLinked <- function(m, rowindex) m$colsLinked[rowindex + 1L]
.hmip_HighsDynamicRowMatrix_getNumRows <- function(m) length(m$ARrange_first)
.hmip_HighsDynamicRowMatrix_getNumDelRows <- function(m) length(m$deletedrows_)
.hmip_HighsDynamicRowMatrix_getRowStart <- function(m, row) m$ARrange_first[row + 1L]
.hmip_HighsDynamicRowMatrix_getRowEnd <- function(m, row) m$ARrange_second[row + 1L]
.hmip_HighsDynamicRowMatrix_getARindex <- function(m) m$ARindex_
.hmip_HighsDynamicRowMatrix_getARvalue <- function(m) m$ARvalue_
.hmip_HighsDynamicRowMatrix_nonzeroCapacity <- function(m) length(m$ARvalue_)

.hmip_HighsDynamicRowMatrix_forEachPositiveColumnEntry <- function(m, col, f) {
  iter <- m$AheadPos_[col + 1L]
  while (iter != -1L) {
    if (!f(m$ARrowindex_[iter + 1L], m$ARvalue_[iter + 1L])) break
    iter <- m$AnextPos_[iter + 1L]
  }
  invisible()
}
.hmip_HighsDynamicRowMatrix_forEachNegativeColumnEntry <- function(m, col, f) {
  iter <- m$AheadNeg_[col + 1L]
  while (iter != -1L) {
    if (!f(m$ARrowindex_[iter + 1L], m$ARvalue_[iter + 1L])) break
    iter <- m$AnextNeg_[iter + 1L]
  }
  invisible()
}

.hmip_HighsDynamicRowMatrix_addRow <- function(m, Rindex, Rvalue, Rlen, linkCols = TRUE) {
  fs <- which(m$freespaces_len >= Rlen)
  if (length(fs) == 0L) {
    start <- length(m$ARindex_)
    end <- start + Rlen
    grow <- end - start
    if (grow > 0L) {
      m$ARindex_ <- c(m$ARindex_, integer(grow))
      m$ARvalue_ <- c(m$ARvalue_, numeric(grow))
      m$ARrowindex_ <- c(m$ARrowindex_, integer(grow))
      m$AprevPos_ <- c(m$AprevPos_, rep(-1L, grow))
      m$AnextPos_ <- c(m$AnextPos_, rep(-1L, grow))
      m$AprevNeg_ <- c(m$AprevNeg_, rep(-1L, grow))
      m$AnextNeg_ <- c(m$AnextNeg_, rep(-1L, grow))
    }
  } else {
    # lower_bound(make_pair(Rlen, -1)): the smallest (len, start) with len >= Rlen
    k <- fs[order(m$freespaces_len[fs], m$freespaces_start[fs])[1L]]
    flen <- m$freespaces_len[k]; fstart <- m$freespaces_start[k]
    m$freespaces_len <- m$freespaces_len[-k]; m$freespaces_start <- m$freespaces_start[-k]
    start <- fstart
    end <- start + Rlen
    if (flen > Rlen) {
      m$freespaces_len <- c(m$freespaces_len, flen - Rlen)
      m$freespaces_start <- c(m$freespaces_start, end)
    }
  }
  start <- as.integer(start); end <- as.integer(end)

  nd <- length(m$deletedrows_)
  if (nd == 0L) {
    rowindex <- length(m$ARrange_first)
    m$ARrange_first <- c(m$ARrange_first, start)
    m$ARrange_second <- c(m$ARrange_second, end)
    m$colsLinked <- c(m$colsLinked, linkCols)
  } else {
    rowindex <- m$deletedrows_[nd]
    m$deletedrows_ <- m$deletedrows_[-nd]
    m$ARrange_first[rowindex + 1L] <- start
    m$ARrange_second[rowindex + 1L] <- end
    m$colsLinked[rowindex + 1L] <- linkCols
  }
  rowindex <- as.integer(rowindex)

  if (Rlen > 0L) {
    pos <- (start + 1L):end
    m$ARindex_[pos] <- as.integer(Rindex[seq_len(Rlen)])
    m$ARvalue_[pos] <- Rvalue[seq_len(Rlen)]
    m$ARrowindex_[pos] <- rowindex
  }

  if (!linkCols) return(rowindex)

  if (Rlen > 0L) for (i in start:(end - 1L)) {
    col <- m$ARindex_[i + 1L]
    if (m$ARvalue_[i + 1L] > 0) {
      m$AprevPos_[i + 1L] <- -1L
      head <- m$AheadPos_[col + 1L]
      m$AheadPos_[col + 1L] <- i
      m$AnextPos_[i + 1L] <- head
      if (head != -1L) m$AprevPos_[head + 1L] <- i
    } else {
      m$AprevNeg_[i + 1L] <- -1L
      head <- m$AheadNeg_[col + 1L]
      m$AheadNeg_[col + 1L] <- i
      m$AnextNeg_[i + 1L] <- head
      if (head != -1L) m$AprevNeg_[head + 1L] <- i
    }
  }
  rowindex
}

.hmip_drm_unlink_range <- function(m, start, end) {
  if (end <= start) return(invisible())
  for (i in start:(end - 1L)) {
    col <- m$ARindex_[i + 1L]
    if (m$ARvalue_[i + 1L] > 0) {
      prev <- m$AprevPos_[i + 1L]; nxt <- m$AnextPos_[i + 1L]
      if (nxt != -1L) m$AprevPos_[nxt + 1L] <- prev
      if (prev != -1L) m$AnextPos_[prev + 1L] <- nxt else m$AheadPos_[col + 1L] <- nxt
    } else {
      prev <- m$AprevNeg_[i + 1L]; nxt <- m$AnextNeg_[i + 1L]
      if (nxt != -1L) m$AprevNeg_[nxt + 1L] <- prev
      if (prev != -1L) m$AnextNeg_[prev + 1L] <- nxt else m$AheadNeg_[col + 1L] <- nxt
    }
  }
  invisible()
}

.hmip_HighsDynamicRowMatrix_unlinkColumns <- function(m, rowindex) {
  if (!m$colsLinked[rowindex + 1L]) return(invisible())
  m$colsLinked[rowindex + 1L] <- FALSE
  .hmip_drm_unlink_range(m, m$ARrange_first[rowindex + 1L], m$ARrange_second[rowindex + 1L])
}

.hmip_HighsDynamicRowMatrix_removeRow <- function(m, rowindex) {
  start <- m$ARrange_first[rowindex + 1L]
  end <- m$ARrange_second[rowindex + 1L]
  if (m$colsLinked[rowindex + 1L]) .hmip_drm_unlink_range(m, start, end)
  m$deletedrows_ <- c(m$deletedrows_, as.integer(rowindex))
  m$freespaces_len <- c(m$freespaces_len, end - start)
  m$freespaces_start <- c(m$freespaces_start, start)
  m$ARrange_first[rowindex + 1L] <- -1L
  m$ARrange_second[rowindex + 1L] <- -1L
  invisible()
}

# ---------------------------------------------------------------------------
# HighsCutSet (mip/HighsCutPool.h)

.hmip_new_HighsCutSet <- function() {
  cs <- new.env(parent = emptyenv())
  cs$cutindices <- integer(); cs$ARstart_ <- integer(); cs$ARindex_ <- integer()
  cs$ARvalue_ <- numeric(); cs$lower_ <- numeric(); cs$upper_ <- numeric()
  cs
}
.hmip_HighsCutSet_numCuts <- function(cs) length(cs$cutindices)
.hmip_HighsCutSet_resize <- function(cs, nnz) {
  ncuts <- length(cs$cutindices)
  cs$lower_ <- c(cs$lower_, rep(-Inf, max(0L, ncuts - length(cs$lower_))))[seq_len(ncuts)]
  cs$upper_ <- c(cs$upper_, numeric(max(0L, ncuts - length(cs$upper_))))[seq_len(ncuts)]
  cs$ARstart_ <- c(cs$ARstart_, integer(max(0L, ncuts + 1L - length(cs$ARstart_))))[seq_len(ncuts + 1L)]
  cs$ARindex_ <- c(cs$ARindex_, integer(max(0L, nnz - length(cs$ARindex_))))[seq_len(nnz)]
  cs$ARvalue_ <- c(cs$ARvalue_, numeric(max(0L, nnz - length(cs$ARvalue_))))[seq_len(nnz)]
  invisible()
}
# clear() keeps lower_ as the C++ does
.hmip_HighsCutSet_clear <- function(cs) {
  cs$cutindices <- integer(); cs$upper_ <- numeric(); cs$ARstart_ <- integer()
  cs$ARindex_ <- integer(); cs$ARvalue_ <- numeric()
  invisible()
}
.hmip_HighsCutSet_empty <- function(cs) length(cs$cutindices) == 0L

# ---------------------------------------------------------------------------
# HighsCutPool (mip/HighsCutPool.h/.cpp)
#
# propRows (std::set<std::pair<age, row>>) holds each row at most once, so it
# is kept as propAge[row + 1] (NA when absent); its descending iteration is
# the order by (age, row) decreasing.

.hmip_new_HighsCutPool <- function(ncols, agelim, softlimit) {
  cp <- new.env(parent = emptyenv())
  cp$matrix_ <- .hmip_new_HighsDynamicRowMatrix(ncols)
  cp$rhs_ <- numeric(); cp$ages_ <- integer(); cp$rownormalization_ <- numeric()
  cp$maxabscoef_ <- numeric(); cp$rowintegral <- logical()
  cp$hashToCutMap <- new.env(parent = emptyenv())
  cp$propagationDomains <- list()
  cp$propAge <- integer()
  cp$agelim_ <- as.integer(agelim); cp$softlimit_ <- as.integer(softlimit)
  cp$numLpCuts <- 0L; cp$numPropNzs <- 0L; cp$numPropRows <- 0L
  cp$ageDistribution <- integer(agelim + 1L)
  cp$minScoreFactor <- 0.9
  cp$bestObservedScore <- 0
  cp$minDensityLim <- 0.1 * ncols
  cp
}

.hmip_cp_prop_erase <- function(cp, age, row) {
  if (row + 1L <= length(cp$propAge) && !is.na(cp$propAge[row + 1L]) && cp$propAge[row + 1L] == age)
    cp$propAge[row + 1L] <- NA_integer_
  invisible()
}
.hmip_cp_prop_emplace <- function(cp, age, row) {
  if (row + 1L > length(cp$propAge)) cp$propAge[row + 1L] <- NA_integer_
  if (is.na(cp$propAge[row + 1L])) cp$propAge[row + 1L] <- as.integer(age)
  invisible()
}
.hmip_cp_prop_size <- function(cp) sum(!is.na(cp$propAge))
.hmip_cp_dist_add <- function(cp, age, d) cp$ageDistribution[age + 1L] <- cp$ageDistribution[age + 1L] + d

.hmip_HighsCutPool_getMatrix <- function(cp) cp$matrix_
.hmip_HighsCutPool_getRhs <- function(cp) cp$rhs_
.hmip_HighsCutPool_cutIsIntegral <- function(cp, cut) cp$rowintegral[cut + 1L]
.hmip_HighsCutPool_getNumCuts <- function(cp)
  .hmip_HighsDynamicRowMatrix_getNumRows(cp$matrix_) - .hmip_HighsDynamicRowMatrix_getNumDelRows(cp$matrix_)
.hmip_HighsCutPool_getNumAvailableCuts <- function(cp) .hmip_HighsCutPool_getNumCuts(cp) - cp$numLpCuts
.hmip_HighsCutPool_getMaxAbsCutCoef <- function(cp, cut) cp$maxabscoef_[cut + 1L]
.hmip_HighsCutPool_getRowNormalization <- function(cp, cut) cp$rownormalization_[cut + 1L]
.hmip_HighsCutPool_getRowLength <- function(cp, row)
  cp$matrix_$ARrange_second[row + 1L] - cp$matrix_$ARrange_first[row + 1L]
.hmip_HighsCutPool_getCut <- function(cp, cut) {
  start <- cp$matrix_$ARrange_first[cut + 1L]
  end <- cp$matrix_$ARrange_second[cut + 1L]
  idx <- if (end > start) (start + 1L):end else integer()
  list(len = end - start, inds = cp$matrix_$ARindex_[idx], vals = cp$matrix_$ARvalue_[idx])
}
.hmip_HighsCutPool_addPropagationDomain <- function(cp, domain) {
  cp$propagationDomains[[length(cp$propagationDomains) + 1L]] <- domain
  invisible()
}
.hmip_HighsCutPool_removePropagationDomain <- function(cp, domain) {
  k <- length(cp$propagationDomains)
  while (k >= 1L) {
    if (identical(cp$propagationDomains[[k]], domain)) {
      cp$propagationDomains[[k]] <- NULL
      return(invisible())
    }
    k <- k - 1L
  }
  invisible()
}
.hmip_HighsCutPool_setAgeLimit <- function(cp, agelim) {
  cp$agelim_ <- as.integer(agelim)
  n <- agelim + 1L
  cp$ageDistribution <- c(cp$ageDistribution, integer(max(0L, n - length(cp$ageDistribution))))[seq_len(n)]
  invisible()
}
.hmip_HighsCutPool_resetAge <- function(cp, cut) {
  age <- cp$ages_[cut + 1L]
  if (age > 0L) {
    if (cp$matrix_$colsLinked[cut + 1L]) {
      .hmip_cp_prop_erase(cp, age, cut)
      .hmip_cp_prop_emplace(cp, 0L, cut)
    }
    .hmip_cp_dist_add(cp, age, -1L)
    .hmip_cp_dist_add(cp, 0L, 1L)
    cp$ages_[cut + 1L] <- 0L
  }
  invisible()
}

# compute_cut_hash: the map key. Two cuts share it exactly when the 64 bit
# hash of the C++ would be equal by construction (same support, same
# double_hash_code of every scaled value); the map is only searched by key.
.hmip_cut_hash_key <- function(Rindex, Rvalue, maxabscoef, Rlen) {
  scale <- 1 / maxabscoef
  codes <- vapply(seq_len(Rlen), function(i) .hmip_double_hash_code(scale * Rvalue[i]), 0)
  paste(c(Rindex[seq_len(Rlen)], "|", codes), collapse = ",")
}

.hmip_HighsCutPool_isDuplicate <- function(cp, hash, norm, Rindex, Rvalue, Rlen, rhs) {
  rows <- cp$hashToCutMap[[hash]]
  m <- cp$matrix_
  for (rowindex in rows) {
    start <- m$ARrange_first[rowindex + 1L]
    end <- m$ARrange_second[rowindex + 1L]
    if (end - start != Rlen) next
    idx <- if (Rlen > 0L) (start + 1L):end else integer()
    if (all(Rindex[seq_len(Rlen)] == m$ARindex_[idx])) {
      dotprod <- 0
      for (i in seq_len(Rlen)) dotprod <- dotprod + Rvalue[i] * m$ARvalue_[start + i]
      parallelism <- dotprod * cp$rownormalization_[rowindex + 1L] * norm
      if (parallelism >= 1 - 1e-6) return(TRUE)
    }
  }
  FALSE
}

.hmip_HighsCutPool_getParallelism <- function(cp, row1, row2) {
  m <- cp$matrix_
  i1 <- m$ARrange_first[row1 + 1L]; end1 <- m$ARrange_second[row1 + 1L]
  i2 <- m$ARrange_first[row2 + 1L]; end2 <- m$ARrange_second[row2 + 1L]
  dotprod <- 0
  while (i1 != end1 && i2 != end2) {
    col1 <- m$ARindex_[i1 + 1L]; col2 <- m$ARindex_[i2 + 1L]
    if (col1 < col2) i1 <- i1 + 1L
    else if (col2 < col1) i2 <- i2 + 1L
    else {
      dotprod <- dotprod + m$ARvalue_[i1 + 1L] * m$ARvalue_[i2 + 1L]
      i1 <- i1 + 1L; i2 <- i2 + 1L
    }
  }
  dotprod * cp$rownormalization_[row1 + 1L] * cp$rownormalization_[row2 + 1L]
}

.hmip_HighsCutPool_lpCutRemoved <- function(cp, cut) {
  if (cp$matrix_$colsLinked[cut + 1L]) {
    .hmip_cp_prop_erase(cp, -1L, cut)
    .hmip_cp_prop_emplace(cp, 1L, cut)
  }
  cp$ages_[cut + 1L] <- 1L
  cp$numLpCuts <- cp$numLpCuts - 1L
  .hmip_cp_dist_add(cp, 1L, 1L)
  invisible()
}

.hmip_cp_cut_deleted <- function(cp, i, onlyProp = FALSE) {
  for (pd in cp$propagationDomains) .hmip_CutpoolPropagation_cutDeleted(pd, i, onlyProp)
}

.hmip_HighsCutPool_performAging <- function(cp) {
  cutIndexEnd <- .hmip_HighsDynamicRowMatrix_getNumRows(cp$matrix_)
  agelim <- cp$agelim_
  numAvailableCuts <- .hmip_HighsCutPool_getNumAvailableCuts(cp)
  while (agelim > 5L && numAvailableCuts > cp$softlimit_) {
    numAvailableCuts <- numAvailableCuts - cp$ageDistribution[agelim + 1L]
    agelim <- agelim - 1L
  }
  for (i in seq_len(cutIndexEnd) - 1L) {
    if (cp$ages_[i + 1L] < 0L) next
    isPropagated <- cp$matrix_$colsLinked[i + 1L]
    if (isPropagated) .hmip_cp_prop_erase(cp, cp$ages_[i + 1L], i)
    .hmip_cp_dist_add(cp, cp$ages_[i + 1L], -1L)
    cp$ages_[i + 1L] <- cp$ages_[i + 1L] + 1L
    if (cp$ages_[i + 1L] > agelim) {
      .hmip_cp_cut_deleted(cp, i)
      if (isPropagated) {
        cp$numPropRows <- cp$numPropRows - 1L
        cp$numPropNzs <- cp$numPropNzs - .hmip_HighsCutPool_getRowLength(cp, i)
      }
      .hmip_HighsDynamicRowMatrix_removeRow(cp$matrix_, i)
      cp$ages_[i + 1L] <- -1L
      cp$rhs_[i + 1L] <- Inf
    } else {
      if (isPropagated) .hmip_cp_prop_emplace(cp, cp$ages_[i + 1L], i)
      .hmip_cp_dist_add(cp, cp$ages_[i + 1L], 1L)
    }
  }
  invisible()
}

.hmip_cp_map_remove <- function(cp, h, i) {
  rows <- cp$hashToCutMap[[h]]
  k <- match(i, rows)
  if (!is.na(k)) {
    rows <- rows[-k]
    if (length(rows)) assign(h, rows, envir = cp$hashToCutMap) else rm(list = h, envir = cp$hashToCutMap)
  }
}

# HighsCutPool::separate(sol, domain, cutset, feastol); sol is the column
# value vector, cutset an env from .hmip_new_HighsCutSet
.hmip_HighsCutPool_separate <- function(cp, sol, domain, cutset, feastol) {
  m <- cp$matrix_
  nrows <- .hmip_HighsDynamicRowMatrix_getNumRows(m)
  eff_score <- numeric(); eff_idx <- integer()

  agelim <- cp$agelim_
  numCuts <- .hmip_HighsCutPool_getNumCuts(cp) - cp$numLpCuts
  while (agelim > 1L && numCuts > cp$softlimit_) {
    numCuts <- numCuts - cp$ageDistribution[agelim + 1L]
    agelim <- agelim - 1L
  }

  for (i in seq_len(nrows) - 1L) {
    if (cp$ages_[i + 1L] < 0L) next
    start <- m$ARrange_first[i + 1L]
    end <- m$ARrange_second[i + 1L]
    viol <- -cp$rhs_[i + 1L]
    if (end > start) for (j in (start + 1L):end) viol <- viol + m$ARvalue_[j] * sol[m$ARindex_[j] + 1L]

    .hmip_cp_dist_add(cp, cp$ages_[i + 1L], -1L)
    isPropagated <- m$colsLinked[i + 1L]
    if (isPropagated) .hmip_cp_prop_erase(cp, cp$ages_[i + 1L], i)
    if (viol <= feastol) {
      cp$ages_[i + 1L] <- cp$ages_[i + 1L] + 1L
      if (cp$ages_[i + 1L] >= agelim) {
        idx <- if (end > start) (start + 1L):end else integer()
        h <- .hmip_cut_hash_key(m$ARindex_[idx], m$ARvalue_[idx], cp$maxabscoef_[i + 1L], end - start)
        .hmip_cp_cut_deleted(cp, i)
        if (isPropagated) {
          cp$numPropRows <- cp$numPropRows - 1L
          cp$numPropNzs <- cp$numPropNzs - (end - start)
        }
        .hmip_HighsDynamicRowMatrix_removeRow(m, i)
        cp$ages_[i + 1L] <- -1L
        cp$rhs_[i + 1L] <- 0
        .hmip_cp_map_remove(cp, h, i)
      } else {
        if (isPropagated) .hmip_cp_prop_emplace(cp, cp$ages_[i + 1L], i)
        .hmip_cp_dist_add(cp, cp$ages_[i + 1L], 1L)
      }
      next
    }

    rownorm <- c(0, 0)
    numActiveNzs <- 0L
    if (end > start) for (j in (start + 1L):end) {
      col <- m$ARindex_[j]
      solval <- sol[col + 1L]
      a <- m$ARvalue_[j]
      if (a > 0) {
        if (solval > domain$col_lower_[col + 1L] + feastol) {
          rownorm <- .cd_add(rownorm, a * a); numActiveNzs <- numActiveNzs + 1L
        }
      } else {
        if (solval < domain$col_upper_[col + 1L] - feastol) {
          rownorm <- .cd_add(rownorm, a * a); numActiveNzs <- numActiveNzs + 1L
        }
      }
    }
    cp$ages_[i + 1L] <- 0L
    .hmip_cp_dist_add(cp, 0L, 1L)
    if (isPropagated) .hmip_cp_prop_emplace(cp, 0L, i)
    score <- viol / (numActiveNzs * sqrt(rownorm[1L] + rownorm[2L]))
    eff_score <- c(eff_score, score); eff_idx <- c(eff_idx, i)
  }
  ne <- length(eff_idx)
  if (ne == 0L) return(invisible())

  # the comparator is a strict total order (score, hash, index), so any
  # correct sort gives the pdqsort result
  hh <- numeric(ne); hl <- numeric(ne)
  for (q in seq_len(ne)) {
    h <- .hmip_hash(ne, eff_idx[q])
    hh[q] <- .hs_hi32(h); hl[q] <- .hs_lo32(h)
  }
  o <- order(eff_score, hh, hl, eff_idx, decreasing = TRUE, method = "radix")
  eff_score <- eff_score[o]; eff_idx <- eff_idx[o]

  cp$bestObservedScore <- max(eff_score[1L], cp$bestObservedScore)
  minScore <- cp$minScoreFactor * cp$bestObservedScore
  numefficacious <- ne
  for (q in seq_len(ne)) if (minScore > eff_score[q]) { numefficacious <- q - 1L; break }

  lowerThreshold <- ne %/% 20L
  upperThreshold <- ne - 1L
  if (numefficacious <= lowerThreshold) {
    numefficacious <- max(ne %/% 2L, 1L)
    cp$minScoreFactor <- eff_score[numefficacious] / cp$bestObservedScore
  } else if (numefficacious > upperThreshold) {
    cp$minScoreFactor <- eff_score[upperThreshold + 1L] / cp$bestObservedScore
  }
  eff_idx <- eff_idx[seq_len(numefficacious)]

  selectednnz <- 0L
  for (cut in eff_idx) {
    discard <- FALSE
    for (k in cutset$cutindices) {
      if (.hmip_HighsCutPool_getParallelism(cp, k, cut) > 0.1) { discard <- TRUE; break }
    }
    if (discard) next
    .hmip_cp_dist_add(cp, cp$ages_[cut + 1L], -1L)
    cp$numLpCuts <- cp$numLpCuts + 1L
    if (m$colsLinked[cut + 1L]) {
      .hmip_cp_prop_erase(cp, cp$ages_[cut + 1L], cut)
      .hmip_cp_prop_emplace(cp, -1L, cut)
    }
    cp$ages_[cut + 1L] <- -1L
    cutset$cutindices <- c(cutset$cutindices, cut)
    selectednnz <- selectednnz + m$ARrange_second[cut + 1L] - m$ARrange_first[cut + 1L]
  }
  .hmip_cp_fill_cutset(cp, cutset, selectednnz)
}

.hmip_cp_fill_cutset <- function(cp, cutset, nnz) {
  m <- cp$matrix_
  .hmip_HighsCutSet_resize(cutset, nnz)
  offset <- 0L
  for (i in seq_along(cutset$cutindices)) {
    cutset$ARstart_[i] <- offset
    cut <- cutset$cutindices[i]
    start <- m$ARrange_first[cut + 1L]; end <- m$ARrange_second[cut + 1L]
    cutset$upper_[i] <- cp$rhs_[cut + 1L]
    len <- end - start
    if (len > 0L) {
      cutset$ARvalue_[offset + seq_len(len)] <- m$ARvalue_[start + seq_len(len)]
      cutset$ARindex_[offset + seq_len(len)] <- m$ARindex_[start + seq_len(len)]
      offset <- offset + len
    }
  }
  cutset$ARstart_[length(cutset$cutindices) + 1L] <- offset
  invisible()
}

.hmip_HighsCutPool_separateLpCutsAfterRestart <- function(cp, cutset) {
  m <- cp$matrix_
  numcuts <- .hmip_HighsDynamicRowMatrix_getNumRows(m)
  cutset$cutindices <- seq_len(numcuts) - 1L
  for (i in seq_len(numcuts) - 1L) {
    .hmip_cp_dist_add(cp, cp$ages_[i + 1L], -1L)
    cp$numLpCuts <- cp$numLpCuts + 1L
    if (m$colsLinked[i + 1L]) {
      .hmip_cp_prop_erase(cp, cp$ages_[i + 1L], i)
      .hmip_cp_prop_emplace(cp, -1L, i)
    }
    cp$ages_[i + 1L] <- -1L
  }
  .hmip_cp_fill_cutset(cp, cutset, .hmip_HighsDynamicRowMatrix_nonzeroCapacity(m))
}

# HighsCutPool::addCut. Rindex and Rvalue are sorted by index in place in the
# C++; the sorted arrays are returned with the new row index (-1 when the cut
# is a duplicate): list(index, Rindex, Rvalue).
.hmip_HighsCutPool_addCut <- function(cp, mipsolver, Rindex, Rvalue, Rlen, rhs, integral = FALSE,
                                      propagate = TRUE, extractCliques = TRUE, isConflict = FALSE) {
  Rlen <- as.integer(Rlen)
  Rindex <- as.integer(Rindex[seq_len(Rlen)]); Rvalue <- Rvalue[seq_len(Rlen)]
  norm <- 0
  maxabscoef <- 0
  for (i in seq_len(Rlen)) {
    norm <- norm + Rvalue[i] * Rvalue[i]
    maxabscoef <- max(maxabscoef, abs(Rvalue[i]))
  }
  # pdqsort_branchless by index: indices are distinct, the order is unique
  o <- order(Rindex, method = "radix")
  Rindex <- Rindex[o]; Rvalue <- Rvalue[o]
  h <- .hmip_cut_hash_key(Rindex, Rvalue, maxabscoef, Rlen)
  normalization <- 1 / sqrt(norm)

  if (.hmip_HighsCutPool_isDuplicate(cp, h, normalization, Rindex, Rvalue, Rlen, rhs))
    return(list(index = -1L, Rindex = Rindex, Rvalue = Rvalue))

  numNonzero <- .hmip_numNonzero(mipsolver)
  if (propagate) {
    newPropNzs <- cp$numPropNzs + Rlen
    avgModelNzs <- numNonzero / .hmip_numRow(mipsolver)
    newAvgPropNzs <- newPropNzs / (cp$numPropRows + 1)
    alpha <- 2
    if (isConflict) {
      if (newAvgPropNzs > max(alpha * avgModelNzs, cp$minDensityLim)) {
        propagate <- FALSE
      } else {
        cp$numPropRows <- cp$numPropRows + 1L
        cp$numPropNzs <- newPropNzs
      }
    } else {
      if (Rlen >= max(alpha * avgModelNzs, cp$minDensityLim)) {
        propagate <- FALSE
      } else {
        cp$numPropRows <- cp$numPropRows + 1L
        cp$numPropNzs <- newPropNzs
      }
    }
  }

  propRowExcessNzs <- cp$numPropNzs - 2 * numNonzero
  if (propRowExcessNzs > 0) {
    pr <- which(!is.na(cp$propAge)) - 1L
    pr <- pr[order(cp$propAge[pr + 1L], pr, decreasing = TRUE, method = "radix")]
    nit <- 0L
    while (propRowExcessNzs > 0 && nit < length(pr)) {
      len <- .hmip_HighsCutPool_getRowLength(cp, pr[nit + 1L])
      propRowExcessNzs <- propRowExcessNzs - len
      cp$numPropNzs <- cp$numPropNzs - len
      cp$numPropRows <- cp$numPropRows - 1L
      nit <- nit + 1L
    }
    for (row in pr[seq_len(nit)]) {
      .hmip_HighsDynamicRowMatrix_unlinkColumns(cp$matrix_, row)
      .hmip_cp_cut_deleted(cp, row, TRUE)
    }
    cp$propAge[pr[seq_len(nit)] + 1L] <- NA_integer_
  }

  rowindex <- .hmip_HighsDynamicRowMatrix_addRow(cp$matrix_, Rindex, Rvalue, Rlen, propagate)
  assign(h, c(cp$hashToCutMap[[h]], rowindex), envir = cp$hashToCutMap)

  if (rowindex == length(cp$rhs_)) {
    cp$rhs_[rowindex + 1L] <- 0; cp$ages_[rowindex + 1L] <- 0L
    cp$rownormalization_[rowindex + 1L] <- 0; cp$maxabscoef_[rowindex + 1L] <- 0
    cp$rowintegral[rowindex + 1L] <- FALSE
  }
  cp$rhs_[rowindex + 1L] <- rhs
  cp$ages_[rowindex + 1L] <- max(0L, cp$agelim_ - 5L)
  .hmip_cp_dist_add(cp, cp$ages_[rowindex + 1L], 1L)
  cp$rowintegral[rowindex + 1L] <- integral
  if (propagate) .hmip_cp_prop_emplace(cp, cp$ages_[rowindex + 1L], rowindex)
  cp$rownormalization_[rowindex + 1L] <- normalization
  cp$maxabscoef_[rowindex + 1L] <- maxabscoef

  for (pd in cp$propagationDomains) .hmip_CutpoolPropagation_cutAdded(pd, rowindex, propagate)

  if (extractCliques && identical(cp, mipsolver$mipdata_$cutpool)) {
    if (Rlen <= 100L)
      .hmip_HighsCliqueTable_extractCliquesFromCut(mipsolver$mipdata_$cliquetable, mipsolver,
                                                   Rindex, Rvalue, Rlen, rhs)
  }
  list(index = rowindex, Rindex = Rindex, Rvalue = Rvalue)
}

# ---------------------------------------------------------------------------
# HighsSparseVectorSum (util/HighsSparseVectorSum.h): values as hi and lo
# vectors, the nonzero index list in a buffer of length dimension

.hmip_new_HighsSparseVectorSum <- function(dimension = 0L) {
  vs <- new.env(parent = emptyenv())
  .hmip_HighsSparseVectorSum_setDimension(vs, dimension)
  vs
}
.hmip_HighsSparseVectorSum_setDimension <- function(vs, dimension) {
  vs$hi <- numeric(dimension); vs$lo <- numeric(dimension)
  vs$nz <- integer(dimension); vs$nnz <- 0L
  invisible()
}
# add(index, value) for a double value or a HighsCDouble c(hi, lo)
.hmip_HighsSparseVectorSum_add <- function(vs, index, value) {
  i <- index + 1L
  h <- vs$hi[i]; l <- vs$lo[i]
  if (h + l != 0) {
    r <- .cd_two_sum(value[1L], h)
    h <- r[1L]
    l <- if (length(value) == 1L) l + r[2L] else l + r[2L] + value[2L]
  } else {
    h <- value[1L]
    l <- if (length(value) == 1L) 0 else value[2L]
    vs$nnz <- vs$nnz + 1L
    vs$nz[vs$nnz] <- as.integer(index)
  }
  if (h + l == 0) { h <- .Machine$double.xmin; l <- 0 }
  vs$hi[i] <- h; vs$lo[i] <- l
  invisible()
}
.hmip_HighsSparseVectorSum_getNonzeros <- function(vs) vs$nz[seq_len(vs$nnz)]
.hmip_HighsSparseVectorSum_getValue <- function(vs, index) vs$hi[index + 1L] + vs$lo[index + 1L]
.hmip_HighsSparseVectorSum_isEmpty <- function(vs) vs$nnz == 0L
.hmip_HighsSparseVectorSum_clear <- function(vs) {
  k <- vs$nz[seq_len(vs$nnz)] + 1L
  vs$hi[k] <- 0; vs$lo[k] <- 0
  vs$nnz <- 0L
  invisible()
}
# cleanup(isZero): isZero(index, value) returns TRUE to drop the entry
.hmip_HighsSparseVectorSum_cleanup <- function(vs, isZero) {
  numNz <- vs$nnz
  nz <- vs$nz
  i <- numNz
  while (i >= 1L) {
    pos <- nz[i]
    val <- vs$hi[pos + 1L] + vs$lo[pos + 1L]
    if (isZero(pos, val)) {
      vs$hi[pos + 1L] <- 0; vs$lo[pos + 1L] <- 0
      t <- nz[numNz]; nz[numNz] <- nz[i]; nz[i] <- t
      numNz <- numNz - 1L
    }
    i <- i - 1L
  }
  vs$nz <- nz
  vs$nnz <- numNz
  invisible()
}

# ---------------------------------------------------------------------------
# HighsLpAggregator (mip/HighsLpAggregator.cpp)

.hmip_new_HighsLpAggregator <- function(lprelaxation) {
  agg <- new.env(parent = emptyenv())
  agg$lprelaxation <- lprelaxation
  lp <- .hmip_cuts_lp(lprelaxation)
  agg$vectorsum <- .hmip_new_HighsSparseVectorSum(lp$num_row_ + lp$num_col_)
  agg
}
.hmip_HighsLpAggregator_addRow <- function(agg, row, weight) {
  r <- .hmip_HighsLpRelaxation_getRow(agg$lprelaxation, row)
  for (i in seq_len(r$len)) .hmip_HighsSparseVectorSum_add(agg$vectorsum, r$inds[i], weight * r$vals[i])
  .hmip_HighsSparseVectorSum_add(agg$vectorsum, .hmip_cuts_lp(agg$lprelaxation)$num_col_ + row, -weight)
  invisible()
}
# returns list(inds, vals)
.hmip_HighsLpAggregator_getCurrentAggregation <- function(agg, negate) {
  droptol <- .hmip_cuts_mip(agg$lprelaxation)$options_mip_$small_matrix_value
  numCol <- .hmip_HighsLpRelaxation_numCols(agg$lprelaxation)
  vs <- agg$vectorsum
  .hmip_HighsSparseVectorSum_cleanup(vs, function(col, val) col < numCol && abs(val) <= droptol)
  inds <- .hmip_HighsSparseVectorSum_getNonzeros(vs)
  vals <- vs$hi[inds + 1L] + vs$lo[inds + 1L]
  if (negate) vals <- -vals
  list(inds = inds, vals = vals)
}
.hmip_HighsLpAggregator_clear <- function(agg) .hmip_HighsSparseVectorSum_clear(agg$vectorsum)
.hmip_HighsLpAggregator_isEmpty <- function(agg) agg$vectorsum$nnz == 0L

# ---------------------------------------------------------------------------
# HighsTransformedLp (mip/HighsTransformedLp.cpp)
#
# BoundType codes: 0 kSimpleUb, 1 kSimpleLb, 2 kVariableUb, 3 kVariableLb
# (enum order; value initialization gives kSimpleUb). bestVub/bestVlb are
# held as column vectors (-1 for none) and the VarBound coef and constant.

.hmip_bt_SimpleUb <- 0L
.hmip_bt_SimpleLb <- 1L
.hmip_bt_VariableUb <- 2L
.hmip_bt_VariableLb <- 3L

# VarBound::minValue / maxValue
.hmip_vb_minValue <- function(coef, constant) .cd_val(.cd_add(c(constant, 0), min(coef, 0)))
.hmip_vb_maxValue <- function(coef, constant) .cd_val(.cd_add(c(constant, 0), max(coef, 0)))

.hmip_new_HighsTransformedLp <- function(lprelaxation, implications) {
  tl <- new.env(parent = emptyenv())
  tl$lprelaxation <- lprelaxation
  mipsolver <- implications$mipsolver
  mipdata <- mipsolver$mipdata_
  lpSolution <- .hmip_cuts_sol(lprelaxation)
  n <- .hmip_HighsLpRelaxation_numCols(lprelaxation) + .hmip_HighsLpRelaxation_numRows(lprelaxation)
  tl$boundDist <- numeric(n); tl$simpleLbDist <- numeric(n); tl$simpleUbDist <- numeric(n)
  tl$lbDist <- numeric(n); tl$ubDist <- numeric(n)
  tl$bestVlb_col <- rep(-1L, n); tl$bestVlb_coef <- numeric(n); tl$bestVlb_constant <- numeric(n)
  tl$bestVub_col <- rep(-1L, n); tl$bestVub_coef <- numeric(n); tl$bestVub_constant <- numeric(n)
  tl$boundTypes <- integer(n)
  tl$vectorsum <- .hmip_new_HighsSparseVectorSum(n)
  feastol <- mipdata$feastol
  dom <- mipdata$domain
  setVub <- function(col, r) {
    tl$bestVub_col[col + 1L] <- r$first
    tl$bestVub_coef[col + 1L] <- r$second$coef
    tl$bestVub_constant[col + 1L] <- r$second$constant
  }
  setVlb <- function(col, r) {
    tl$bestVlb_col[col + 1L] <- r$first
    tl$bestVlb_coef[col + 1L] <- r$second$coef
    tl$bestVlb_constant[col + 1L] <- r$second$constant
  }

  for (col in mipdata$continuous_cols) {
    j <- col + 1L
    .hmip_HighsImplications_cleanupVarbounds(mipdata$implications, col)
    if (.hmip_HighsDomain_infeasible(dom)) return(tl)
    if (.hmip_HighsDomain_isFixed(dom, col)) next
    x <- lpSolution$col_value[j]

    bestub <- dom$col_upper_[j]
    tl$simpleUbDist[j] <- bestub - x
    if (tl$simpleUbDist[j] <= feastol) tl$simpleUbDist[j] <- 0
    r <- .hmip_HighsImplications_getBestVub(implications, col, lpSolution, bestub)
    setVub(col, r); bestub <- r$bestUb

    bestlb <- dom$col_lower_[j]
    tl$simpleLbDist[j] <- x - bestlb
    if (tl$simpleLbDist[j] <= feastol) tl$simpleLbDist[j] <- 0
    r <- .hmip_HighsImplications_getBestVlb(implications, col, lpSolution, bestlb)
    setVlb(col, r); bestlb <- r$bestLb

    tl$lbDist[j] <- x - bestlb
    if (tl$lbDist[j] <= feastol) tl$lbDist[j] <- 0
    tl$ubDist[j] <- bestub - x
    if (tl$ubDist[j] <= feastol) tl$ubDist[j] <- 0
    tl$boundDist[j] <- min(tl$lbDist[j], tl$ubDist[j])
  }

  for (col in mipdata$integral_cols) {
    j <- col + 1L
    bestub <- dom$col_upper_[j]
    bestlb <- dom$col_lower_[j]
    .hmip_HighsImplications_cleanupVarbounds(mipdata$implications, col)
    if (.hmip_HighsDomain_infeasible(dom)) return(tl)
    x <- lpSolution$col_value[j]
    tl$simpleUbDist[j] <- bestub - x
    if (tl$simpleUbDist[j] <= feastol) tl$simpleUbDist[j] <- 0
    tl$simpleLbDist[j] <- x - bestlb
    if (tl$simpleLbDist[j] <= feastol) tl$simpleLbDist[j] <- 0
    simpleBndDist <- min(tl$simpleLbDist[j], tl$simpleUbDist[j])
    if (simpleBndDist > 0 && abs(.hmip_nearestInteger(x) - x) < feastol) {
      r <- .hmip_HighsImplications_getBestVub(mipdata$implications, col, lpSolution, bestub)
      setVub(col, r); bestub <- r$bestUb
      r <- .hmip_HighsImplications_getBestVlb(mipdata$implications, col, lpSolution, bestlb)
      setVlb(col, r); bestlb <- r$bestLb
      tl$lbDist[j] <- x - bestlb
      if (tl$lbDist[j] <= feastol) tl$lbDist[j] <- 0
      tl$ubDist[j] <- bestub - x
      if (tl$ubDist[j] <= feastol) tl$ubDist[j] <- 0
      tl$boundDist[j] <- min(tl$lbDist[j], tl$ubDist[j])
      if (tl$boundDist[j] > simpleBndDist + feastol) {
        tl$lbDist[j] <- tl$simpleLbDist[j]
        tl$ubDist[j] <- tl$simpleUbDist[j]
        tl$boundDist[j] <- simpleBndDist
        tl$bestVub_col[j] <- -1L
        tl$bestVlb_col[j] <- -1L
      }
    } else {
      tl$lbDist[j] <- tl$simpleLbDist[j]
      tl$ubDist[j] <- tl$simpleUbDist[j]
      tl$boundDist[j] <- simpleBndDist
    }
  }

  numLpRow <- .hmip_HighsLpRelaxation_numRows(lprelaxation)
  indexOffset <- .hmip_numCol(mipsolver)
  for (row in seq_len(numLpRow) - 1L) {
    s <- indexOffset + row + 1L
    bestub <- .hmip_HighsLpRelaxation_slackUpper(lprelaxation, row)
    bestlb <- .hmip_HighsLpRelaxation_slackLower(lprelaxation, row)
    if (bestlb == bestub) next
    rv <- lpSolution$row_value[row + 1L]
    tl$lbDist[s] <- rv - bestlb
    if (tl$lbDist[s] <= feastol) tl$lbDist[s] <- 0
    tl$simpleLbDist[s] <- tl$lbDist[s]
    tl$ubDist[s] <- bestub - rv
    if (tl$ubDist[s] <= feastol) tl$ubDist[s] <- 0
    tl$simpleUbDist[s] <- tl$ubDist[s]
    tl$boundDist[s] <- min(tl$lbDist[s], tl$ubDist[s])
  }
  tl
}

.hmip_HighsTransformedLp_boundDistance <- function(tl, col) tl$boundDist[col + 1L]

# transform(vals, upper, solval, inds, rhs, integersPositive, preferVbds):
# returns list(ok, vals, upper, solval, inds, rhs, integersPositive)
.hmip_HighsTransformedLp_transform <- function(tl, vals, upper, solval, inds, rhs, integersPositive,
                                               preferVbds = FALSE) {
  lpr <- tl$lprelaxation
  mip <- .hmip_cuts_mip(lpr)
  mipdata <- mip$mipdata_
  feastol <- mipdata$feastol
  smallval <- mip$options_mip_$small_matrix_value
  dom <- mipdata$domain
  vs <- tl$vectorsum
  out <- function(ok) list(ok = ok, vals = vals, upper = upper, solval = solval, inds = inds,
                           rhs = rhs, integersPositive = integersPositive)

  tmpRhs <- c(rhs, 0)
  slackOffset <- .hmip_HighsLpRelaxation_numCols(lpr)
  numNz <- length(inds)
  getLb <- function(col) if (col < slackOffset) dom$col_lower_[col + 1L] else
    .hmip_HighsLpRelaxation_slackLower(lpr, col - slackOffset)
  getUb <- function(col) if (col < slackOffset) dom$col_upper_[col + 1L] else
    .hmip_HighsLpRelaxation_slackUpper(lpr, col - slackOffset)
  isColIntegral <- function(col) .hmip_HighsLpRelaxation_isColIntegral(lpr, col)

  i <- 0L
  while (i < numNz) {
    p <- i + 1L
    col <- inds[p]
    j <- col + 1L
    lb <- getLb(col)
    ub <- getUb(col)
    remove <- FALSE

    if (ub - lb < smallval) {
      tmpRhs <- .cd_add(tmpRhs, -(min(lb, ub) * vals[p]))
      remove <- TRUE
    } else {
      if (lb == -Inf && ub == Inf) {
        .hmip_HighsSparseVectorSum_clear(vs)
        return(out(FALSE))
      }
      if (tl$bestVub_col[j] != -1L &&
          .hmip_vb_maxValue(tl$bestVub_coef[j], tl$bestVub_constant[j]) > ub + feastol) {
        r <- .hmip_HighsImplications_cleanupVub(mipdata$implications, col, tl$bestVub_col[j],
                                                list(coef = tl$bestVub_coef[j], constant = tl$bestVub_constant[j]),
                                                ub, FALSE)
        tl$bestVub_coef[j] <- r$vub$coef; tl$bestVub_constant[j] <- r$vub$constant
      }
      if (tl$bestVlb_col[j] != -1L &&
          .hmip_vb_minValue(tl$bestVlb_coef[j], tl$bestVlb_constant[j]) < lb - feastol) {
        r <- .hmip_HighsImplications_cleanupVlb(mipdata$implications, col, tl$bestVlb_col[j],
                                                list(coef = tl$bestVlb_coef[j], constant = tl$bestVlb_constant[j]),
                                                lb, FALSE)
        tl$bestVlb_coef[j] <- r$vlb$coef; tl$bestVlb_constant[j] <- r$vlb$constant
      }

      oldBoundType <- tl$boundTypes[j]
      skipSwitch <- FALSE
      if (isColIntegral(col)) {
        if (ub - lb <= 1.5 || tl$boundDist[j] != 0 || tl$simpleLbDist[j] == 0 || tl$simpleUbDist[j] == 0) {
          if (tl$simpleLbDist[j] < tl$simpleUbDist[j] - feastol) tl$boundTypes[j] <- .hmip_bt_SimpleLb
          else if (tl$simpleUbDist[j] < tl$simpleLbDist[j] - feastol) tl$boundTypes[j] <- .hmip_bt_SimpleUb
          else if (vals[p] > 0) tl$boundTypes[j] <- .hmip_bt_SimpleLb
          else tl$boundTypes[j] <- .hmip_bt_SimpleUb
          i <- i + 1L
          next
        }
        if (tl$bestVlb_col[j] == -1L || tl$ubDist[j] < tl$lbDist[j] - feastol)
          tl$boundTypes[j] <- .hmip_bt_VariableUb
        else if (tl$bestVub_col[j] == -1L || tl$lbDist[j] < tl$ubDist[j] - feastol)
          tl$boundTypes[j] <- .hmip_bt_VariableLb
        else if (vals[p] > 0)
          tl$boundTypes[j] <- .hmip_bt_VariableUb
        else
          tl$boundTypes[j] <- .hmip_bt_VariableLb
      } else {
        if (tl$lbDist[j] < tl$ubDist[j] - feastol) {
          if (tl$bestVlb_col[j] == -1L) tl$boundTypes[j] <- .hmip_bt_SimpleLb
          else if (preferVbds || vals[p] > 0 || tl$simpleLbDist[j] > tl$lbDist[j] + feastol)
            tl$boundTypes[j] <- .hmip_bt_VariableLb
          else tl$boundTypes[j] <- .hmip_bt_SimpleLb
        } else if (tl$ubDist[j] < tl$lbDist[j] - feastol) {
          if (tl$bestVub_col[j] == -1L) tl$boundTypes[j] <- .hmip_bt_SimpleUb
          else if (preferVbds || vals[p] < 0 || tl$simpleUbDist[j] > tl$ubDist[j] + feastol)
            tl$boundTypes[j] <- .hmip_bt_VariableUb
          else tl$boundTypes[j] <- .hmip_bt_SimpleUb
        } else if (vals[p] > 0) {
          if (tl$bestVlb_col[j] != -1L) tl$boundTypes[j] <- .hmip_bt_VariableLb
          else if (preferVbds && tl$bestVub_col[j] != -1L) tl$boundTypes[j] <- .hmip_bt_VariableUb
          else tl$boundTypes[j] <- .hmip_bt_SimpleLb
        } else {
          if (tl$bestVub_col[j] != -1L) tl$boundTypes[j] <- .hmip_bt_VariableUb
          else if (preferVbds && tl$bestVlb_col[j] != -1L) tl$boundTypes[j] <- .hmip_bt_VariableLb
          else tl$boundTypes[j] <- .hmip_bt_SimpleUb
        }
      }

      bt <- tl$boundTypes[j]
      if (bt == .hmip_bt_SimpleLb) {
        if (vals[p] > 0) {
          tmpRhs <- .cd_add(tmpRhs, -(lb * vals[p]))
          tl$boundTypes[j] <- oldBoundType
          remove <- TRUE
        }
      } else if (bt == .hmip_bt_SimpleUb) {
        if (vals[p] < 0) {
          tmpRhs <- .cd_add(tmpRhs, -(ub * vals[p]))
          tl$boundTypes[j] <- oldBoundType
          remove <- TRUE
        }
      } else if (bt == .hmip_bt_VariableLb) {
        tmpRhs <- .cd_add(tmpRhs, -(tl$bestVlb_constant[j] * vals[p]))
        .hmip_HighsSparseVectorSum_add(vs, tl$bestVlb_col[j], vals[p] * tl$bestVlb_coef[j])
        if (vals[p] > 0) {
          tl$boundTypes[j] <- oldBoundType
          remove <- TRUE
        }
      } else {
        tmpRhs <- .cd_add(tmpRhs, -(tl$bestVub_constant[j] * vals[p]))
        .hmip_HighsSparseVectorSum_add(vs, tl$bestVub_col[j], vals[p] * tl$bestVub_coef[j])
        vals[p] <- -vals[p]
        if (vals[p] > 0) {
          tl$boundTypes[j] <- oldBoundType
          remove <- TRUE
        }
      }
    }
    if (remove) {
      inds[p] <- inds[numNz]; vals[p] <- vals[numNz]
      inds[numNz] <- 0L; vals[numNz] <- 0
      numNz <- numNz - 1L
      next
    }
    i <- i + 1L
  }

  if (vs$nnz != 0L) {
    for (q in seq_len(numNz)) if (vals[q] != 0) .hmip_HighsSparseVectorSum_add(vs, inds[q], vals[q])
    .hmip_HighsSparseVectorSum_cleanup(vs, function(col, val) abs(val) <= smallval)
    inds <- .hmip_HighsSparseVectorSum_getNonzeros(vs)
    numNz <- length(inds)
    vals <- vs$hi[inds + 1L] + vs$lo[inds + 1L]
    .hmip_HighsSparseVectorSum_clear(vs)
  } else {
    vals <- vals[seq_len(numNz)]
    inds <- inds[seq_len(numNz)]
  }

  for (q in seq_len(numNz)) {
    col <- inds[q]
    j <- col + 1L
    if (!isColIntegral(col)) next
    if (tl$boundTypes[j] == .hmip_bt_VariableLb || tl$boundTypes[j] == .hmip_bt_VariableUb) next
    lb <- getLb(col); ub <- getUb(col)
    if (integersPositive) {
      tl$boundTypes[j] <- if ((lb != -Inf && vals[q] > 0) || ub == Inf) .hmip_bt_SimpleLb else .hmip_bt_SimpleUb
    } else {
      tl$boundTypes[j] <- if (tl$lbDist[j] < tl$ubDist[j]) .hmip_bt_SimpleLb else .hmip_bt_SimpleUb
    }
  }

  upper <- numeric(numNz)
  solval <- numeric(numNz)
  for (q in seq_len(numNz)) {
    col <- inds[q]
    j <- col + 1L
    lb <- getLb(col); ub <- getUb(col)
    upper[q] <- ub - lb
    bt <- tl$boundTypes[j]
    if (bt == .hmip_bt_SimpleLb) {
      tmpRhs <- .cd_add(tmpRhs, -(lb * vals[q]))
      solval[q] <- tl$lbDist[j]
    } else if (bt == .hmip_bt_SimpleUb) {
      tmpRhs <- .cd_add(tmpRhs, -(ub * vals[q]))
      vals[q] <- -vals[q]
      solval[q] <- tl$ubDist[j]
    } else if (bt == .hmip_bt_VariableLb) {
      solval[q] <- tl$lbDist[j]
    } else {
      solval[q] <- tl$ubDist[j]
    }
    if (isColIntegral(col)) integersPositive <- integersPositive && vals[q] > 0
  }

  rhs <- tmpRhs[1L] + tmpRhs[2L]
  if (numNz == 0L && rhs >= -feastol) return(out(FALSE))
  out(TRUE)
}

# untransform(vals, inds, rhs, integral): returns list(ok, vals, inds, rhs)
.hmip_HighsTransformedLp_untransform <- function(tl, vals, inds, rhs, integral = FALSE) {
  lpr <- tl$lprelaxation
  mip <- .hmip_cuts_mip(lpr)
  mipdata <- mip$mipdata_
  dom <- mipdata$domain
  vs <- tl$vectorsum
  tmpRhs <- c(rhs, 0)
  slackOffset <- .hmip_numCol(mip)
  add <- function(i, v) .hmip_HighsSparseVectorSum_add(vs, i, v)

  for (q in seq_along(inds)) {
    if (vals[q] == 0) next
    col <- inds[q]
    j <- col + 1L
    bt <- tl$boundTypes[j]
    if (bt == .hmip_bt_VariableLb) {
      tmpRhs <- .cd_add(tmpRhs, tl$bestVlb_constant[j] * vals[q])
      add(tl$bestVlb_col[j], -vals[q] * tl$bestVlb_coef[j])
      add(col, vals[q])
    } else if (bt == .hmip_bt_VariableUb) {
      tmpRhs <- .cd_add(tmpRhs, -(tl$bestVub_constant[j] * vals[q]))
      add(tl$bestVub_col[j], vals[q] * tl$bestVub_coef[j])
      add(col, -vals[q])
    } else if (bt == .hmip_bt_SimpleLb) {
      if (col < slackOffset) {
        tmpRhs <- .cd_add(tmpRhs, vals[q] * dom$col_lower_[j])
        add(col, vals[q])
      } else {
        row <- col - slackOffset
        tmpRhs <- .cd_add(tmpRhs, vals[q] * .hmip_HighsLpRelaxation_slackLower(lpr, row))
        r <- .hmip_HighsLpRelaxation_getRow(lpr, row)
        for (k in seq_len(r$len)) add(r$inds[k], vals[q] * r$vals[k])
      }
    } else {
      if (col < slackOffset) {
        tmpRhs <- .cd_add(tmpRhs, -(vals[q] * dom$col_upper_[j]))
        add(col, -vals[q])
      } else {
        row <- col - slackOffset
        tmpRhs <- .cd_add(tmpRhs, -(vals[q] * .hmip_HighsLpRelaxation_slackUpper(lpr, row)))
        vals[q] <- -vals[q]
        r <- .hmip_HighsLpRelaxation_getRow(lpr, row)
        for (k in seq_len(r$len)) add(r$inds[k], vals[q] * r$vals[k])
      }
    }
  }

  if (integral) {
    .hmip_HighsSparseVectorSum_cleanup(vs, function(col, val) abs(val) < 0.5)
    rhs <- .hmip_round(tmpRhs[1L] + tmpRhs[2L])
  } else {
    abort <- FALSE
    smallval <- mip$options_mip_$small_matrix_value
    feastol <- mipdata$feastol
    isZero <- function(col, val) {
      absval <- abs(val)
      if (absval <= smallval) return(TRUE)
      if (absval <= feastol) {
        if (val > 0) {
          if (dom$col_lower_[col + 1L] == -Inf) abort <<- TRUE
          else tmpRhs <<- .cd_add(tmpRhs, -(val * dom$col_lower_[col + 1L]))
        } else {
          if (dom$col_upper_[col + 1L] == Inf) abort <<- TRUE
          else tmpRhs <<- .cd_add(tmpRhs, -(val * dom$col_upper_[col + 1L]))
        }
        return(TRUE)
      }
      FALSE
    }
    .hmip_HighsSparseVectorSum_cleanup(vs, isZero)
    if (abort) {
      .hmip_HighsSparseVectorSum_clear(vs)
      return(list(ok = FALSE, vals = vals, inds = inds, rhs = rhs))
    }
    rhs <- tmpRhs[1L] + tmpRhs[2L]
  }

  inds <- .hmip_HighsSparseVectorSum_getNonzeros(vs)
  vals <- vs$hi[inds + 1L] + vs$lo[inds + 1L]
  if (integral) vals <- vapply(vals, .hmip_round, 0)
  .hmip_HighsSparseVectorSum_clear(vs)
  list(ok = TRUE, vals = vals, inds = inds, rhs = rhs)
}

# ---------------------------------------------------------------------------
# HighsCutGeneration (mip/HighsCutGeneration.cpp)
#
# The C++ works on raw pointers vals/inds into caller vectors; here cg$vals
# and cg$inds are the working arrays (the content of the vector the pointer
# designates) and the public methods return the caller vectors in a list.
# complementation of length 0 is the empty vector of the C++.

.hmip_new_HighsCutGeneration <- function(lpRelaxation, cutpool) {
  cg <- new.env(parent = emptyenv())
  mip <- .hmip_cuts_mip(lpRelaxation)
  cg$lpRelaxation <- lpRelaxation
  cg$cutpool <- cutpool
  seed <- (mip$options_mip_$random_seed + .hmip_HighsLpRelaxation_getNumLpIterations(lpRelaxation) +
             .hmip_HighsCutPool_getNumCuts(cutpool)) %% 4294967296
  cg$randgen <- .hs_rng(seed)
  cg$feastol <- mip$mipdata_$feastol
  cg$epsilon <- mip$mipdata_$epsilon
  cg$cover <- integer(); cg$coverweight <- c(0, 0); cg$lambda <- c(0, 0)
  cg$upper <- numeric(); cg$solval <- numeric(); cg$complementation <- integer()
  cg$isintegral <- logical(); cg$vals <- numeric(); cg$inds <- integer()
  cg$rhs <- c(0, 0); cg$integralSupport <- FALSE; cg$integralCoefficients <- FALSE
  cg$rowlen <- 0L; cg$initialScale <- 1
  cg$integerinds <- integer(); cg$deltas <- numeric()
  cg
}

.hmip_HighsCutGeneration_determineCover <- function(cg, lpSol = TRUE) {
  feastol <- cg$feastol
  if (.cd_val(cg$rhs) <= 10 * feastol) return(FALSE)
  rowlen <- cg$rowlen
  vals <- cg$vals; upper <- cg$upper; solval <- cg$solval; inds <- cg$inds
  cover <- integer()
  for (j in seq_len(rowlen) - 1L) {
    if (!cg$isintegral[j + 1L]) next
    if (lpSol && solval[j + 1L] <= feastol) next
    cover <- c(cover, j)
  }
  maxCoverSize <- length(cover)
  coversize <- 0L
  r <- .hmip_rand_int31(cg$randgen)
  coverweight <- c(0, 0)
  hcache <- new.env(parent = emptyenv())
  hashOf <- function(j) {
    key <- as.character(j)
    h <- hcache[[key]]
    if (is.null(h)) { h <- .hmip_hash(inds[j + 1L], r); assign(key, h, envir = hcache) }
    h
  }
  if (lpSol) {
    p <- .hmip_std_partition(cover, function(j) solval[j + 1L] >= upper[j + 1L] - feastol)
    cover <- p$v
    coversize <- p$k
    for (i in seq_len(coversize)) {
      j <- cover[i]
      coverweight <- .cd_add(coverweight, vals[j + 1L] * upper[j + 1L])
    }
    if (maxCoverSize > coversize) {
      less <- function(i, j) {
        if (upper[i + 1L] < 1.5 && upper[j + 1L] > 1.5) return(TRUE)
        if (upper[i + 1L] > 1.5 && upper[j + 1L] < 1.5) return(FALSE)
        contributionA <- solval[i + 1L] * vals[i + 1L]
        contributionB <- solval[j + 1L] * vals[j + 1L]
        if (contributionA > contributionB + feastol) return(TRUE)
        if (contributionA < contributionB - feastol) return(FALSE)
        if (abs(vals[i + 1L] - vals[j + 1L]) <= feastol)
          return(.hmip_u64_cmp(hashOf(i), hashOf(j)) > 0L)
        vals[i + 1L] > vals[j + 1L]
      }
      rng <- (coversize + 1L):maxCoverSize
      cover[rng] <- .hmip_pdqsort(cover[rng], less)
    }
  } else {
    nodequeue <- .hmip_cuts_mip(cg$lpRelaxation)$mipdata_$nodequeue
    compl <- cg$complementation
    numNodes <- function(i) {
      if (length(compl) && compl[i + 1L]) .hmip_HighsNodeQueue_numNodesDown(nodequeue, inds[i + 1L])
      else .hmip_HighsNodeQueue_numNodesUp(nodequeue, inds[i + 1L])
    }
    if (maxCoverSize > coversize) {
      less <- function(i, j) {
        if (solval[i + 1L] > feastol && solval[j + 1L] <= feastol) return(TRUE)
        if (solval[i + 1L] <= feastol && solval[j + 1L] > feastol) return(FALSE)
        numNodesA <- numNodes(i)
        numNodesB <- numNodes(j)
        if (numNodesA > numNodesB) return(TRUE)
        if (numNodesA < numNodesB) return(FALSE)
        .hmip_u64_cmp(hashOf(i), hashOf(j)) > 0L
      }
      rng <- (coversize + 1L):maxCoverSize
      cover[rng] <- .hmip_pdqsort(cover[rng], less)
    }
  }

  minlambda <- max(10 * feastol, feastol * abs(.cd_val(cg$rhs)))
  while (coversize != maxCoverSize) {
    lambda <- .cd_val(.cd_sub_cd(coverweight, cg$rhs))
    if (lambda > minlambda) break
    j <- cover[coversize + 1L]
    coverweight <- .cd_add(coverweight, vals[j + 1L] * upper[j + 1L])
    coversize <- coversize + 1L
  }
  cg$cover <- cover
  if (coversize == 0L) return(FALSE)
  coverweight <- .cd_two_sum(coverweight[1L], coverweight[2L])
  cg$coverweight <- coverweight
  cg$lambda <- .cd_sub_cd(coverweight, cg$rhs)
  if (.cd_val(cg$lambda) <= minlambda) return(FALSE)
  cg$cover <- cover[seq_len(coversize)]
  TRUE
}

.hmip_cg_sort_cover_desc <- function(cg) {
  vals <- cg$vals
  cg$cover <- .hmip_pdqsort(cg$cover, function(a, b) vals[a + 1L] > vals[b + 1L], branchless = TRUE)
}

.hmip_HighsCutGeneration_separateLiftedKnapsackCover <- function(cg) {
  feastol <- cg$feastol; epsilon <- cg$epsilon
  rowlen <- cg$rowlen
  coversize <- length(cg$cover)
  S <- numeric(coversize)
  coverflag <- integer(rowlen)
  .hmip_cg_sort_cover_desc(cg)
  cover <- cg$cover; vals <- cg$vals

  abartmp <- c(vals[cover[1L] + 1L], 0)
  sigma <- cg$lambda
  if (coversize >= 2L) for (i in 1:(coversize - 1L)) {
    delta <- .cd_add(abartmp, -vals[cover[i + 1L] + 1L])
    kdelta <- .cd_mul(delta, i)
    if (.cd_val(kdelta) < .cd_val(sigma)) {
      abartmp <- c(vals[cover[i + 1L] + 1L], 0)
      sigma <- .cd_sub_cd(sigma, kdelta)
    } else {
      abartmp <- .cd_sub_cd(abartmp, .cd_mul(sigma, 1 / i))
      sigma <- c(0, 0)
      break
    }
  }
  if (.cd_val(sigma) > 0) abartmp <- .cd_div(cg$rhs, coversize)
  abar <- .cd_val(abartmp)

  sum <- c(0, 0)
  cplussize <- 0L
  for (i in seq_len(coversize)) {
    v <- vals[cover[i] + 1L]
    sum <- .cd_add(sum, min(abar, v))
    S[i] <- .cd_val(sum)
    if (v > abar + feastol) {
      cplussize <- cplussize + 1L
      coverflag[cover[i] + 1L] <- 1L
    } else coverflag[cover[i] + 1L] <- -1L
  }
  halfintegral <- FALSE
  g <- function(z) {
    hfrac <- z / abar
    coef <- 0
    h <- floor(hfrac + 0.5)
    if (h != 0 && abs(hfrac - h) * max(1, abar) <= epsilon && h <= cplussize - 1L) {
      halfintegral <<- TRUE
      coef <- 0.5
    }
    h <- max(h - 1, 0)
    while (h < coversize) {
      if (z <= S[h + 1L] + feastol) break
      h <- h + 1
    }
    coef + h
  }

  cg$rhs <- c(coversize - 1, 0)
  for (i in seq_len(rowlen)) {
    if (vals[i] == 0) next
    if (coverflag[i] == -1L) vals[i] <- 1 else vals[i] <- g(vals[i])
  }
  if (halfintegral) {
    cg$rhs <- .cd_mul(cg$rhs, 2)
    vals[seq_len(rowlen)] <- vals[seq_len(rowlen)] * 2
  }
  cg$vals <- vals
  cg$integralSupport <- TRUE
  cg$integralCoefficients <- TRUE
  invisible()
}

.hmip_HighsCutGeneration_separateLiftedMixedBinaryCover <- function(cg) {
  cg$integralSupport <- FALSE
  cg$integralCoefficients <- FALSE
  epsilon <- cg$epsilon
  rowlen <- cg$rowlen
  coversize <- length(cg$cover)
  S <- numeric(coversize)
  coverflag <- logical(rowlen)
  if (coversize == 0L) return(FALSE)
  coverflag[cg$cover + 1L] <- TRUE
  .hmip_cg_sort_cover_desc(cg)
  cover <- cg$cover; vals <- cg$vals; lambda <- cg$lambda

  sum <- c(0, 0)
  p <- coversize
  for (i in seq_len(coversize)) {
    v <- vals[cover[i] + 1L]
    if (.cd_val(.hmip_cd_rsub(v, lambda)) <= epsilon) { p <- i - 1L; break }
    sum <- .cd_add(sum, v)
    S[i] <- .cd_val(sum)
  }
  if (p == 0L) return(FALSE)
  phi <- function(a) {
    for (i in seq_len(p) - 1L) {
      if (a <= .cd_val(.hmip_cd_rsub(S[i + 1L], lambda))) return(.cd_val(.cd_mul(lambda, i)))
      if (a <= S[i + 1L])
        return(.cd_val(.cd_add_cd(.cd_mul(lambda, i + 1), .cd_add(c(a, 0), -S[i + 1L]))))
    }
    .cd_val(.cd_add_cd(.cd_mul(lambda, p), .cd_add(c(a, 0), -S[p])))
  }

  rhs <- -lambda
  integralSupport <- TRUE
  for (i in seq_len(rowlen)) {
    if (!cg$isintegral[i]) {
      if (vals[i] < 0) integralSupport <- FALSE else vals[i] <- 0
      next
    }
    if (coverflag[i]) {
      vals[i] <- min(vals[i], .cd_val(lambda))
      rhs <- .cd_add(rhs, vals[i])
    } else {
      vals[i] <- phi(vals[i])
    }
  }
  cg$vals <- vals
  cg$rhs <- rhs
  cg$integralCoefficients <- FALSE
  cg$integralSupport <- integralSupport
  TRUE
}

.hmip_HighsCutGeneration_separateLiftedMixedIntegerCover <- function(cg) {
  cg$integralSupport <- FALSE
  cg$integralCoefficients <- FALSE
  feastol <- cg$feastol; epsilon <- cg$epsilon
  rowlen <- cg$rowlen
  coversize <- length(cg$cover)
  coverflag <- logical(rowlen)
  coverflag[cg$cover + 1L] <- TRUE
  .hmip_cg_sort_cover_desc(cg)
  cover <- cg$cover; vals <- cg$vals; upper <- cg$upper; solval <- cg$solval
  lambda <- cg$lambda

  a <- vector("list", coversize)
  u <- vector("list", coversize + 1L)
  m <- vector("list", coversize + 1L)
  usum <- c(0, 0); msum <- c(0, 0)
  for (c0 in seq_len(coversize)) {
    i <- cover[c0]
    u[[c0]] <- usum
    m[[c0]] <- msum
    a[[c0]] <- c(vals[i + 1L], 0)
    ub <- upper[i + 1L]
    usum <- .cd_add(usum, ub)
    msum <- .cd_add_cd(msum, .cd_mul(a[[c0]], ub))
  }
  u[[coversize + 1L]] <- usum
  m[[coversize + 1L]] <- msum

  lpos <- -1L
  bestlCplusend <- -1L
  bestlVal <- 0
  bestlAtUpper <- TRUE
  coverVals <- vals[cover + 1L]
  for (i in seq_len(coversize) - 1L) {
    j <- cover[i + 1L]
    ub <- upper[j + 1L]
    atUpper <- solval[j + 1L] >= ub - feastol
    if (atUpper && !bestlAtUpper) next
    mju <- ub * vals[j + 1L]
    mu <- .hmip_cd_rsub(mju, lambda)
    if (.cd_val(mu) <= 10 * feastol) next
    if (abs(vals[j + 1L]) < 1000 * feastol) next
    mudival <- .cd_val(.cd_div(mu, vals[j + 1L]))
    if (.hmip_cuts_isIntegral(mudival, feastol)) next
    eta <- ceiling(mudival)
    ulminusetaplusone <- .cd_add(.cd_add(c(ub, 0), -eta), 1)
    cplusthreshold <- .cd_val(.cd_mul(ulminusetaplusone, vals[j + 1L]))
    # std::upper_bound with comp(t, i) = t > vals[i] on the sorted cover
    w <- which(cplusthreshold > coverVals)
    cplusend <- if (length(w)) w[1L] - 1L else coversize
    mcplus <- m[[cplusend + 1L]]
    if (i < cplusend) mcplus <- .cd_add(mcplus, -mju)
    jlVal <- .cd_val(.cd_add(mcplus, eta * vals[j + 1L]))
    if (jlVal > bestlVal || (!atUpper && bestlAtUpper)) {
      lpos <- i
      bestlCplusend <- cplusend
      bestlVal <- jlVal
      bestlAtUpper <- atUpper
    }
  }
  if (lpos == -1L) return(FALSE)

  l <- cover[lpos + 1L]
  al <- c(vals[l + 1L], 0)
  upperl <- upper[l + 1L]
  mlu <- .cd_mul(al, upperl)
  mu <- .cd_sub_cd(mlu, lambda)

  a <- a[seq_len(bestlCplusend)]
  cover <- cover[seq_len(bestlCplusend)]
  u <- u[seq_len(bestlCplusend + 1L)]
  m <- m[seq_len(bestlCplusend + 1L)]
  if (lpos < bestlCplusend) {
    a <- a[-(lpos + 1L)]
    cover <- cover[-(lpos + 1L)]
    u <- u[-(lpos + 2L)]
    m <- m[-(lpos + 2L)]
    if (lpos + 1L < bestlCplusend) for (i in (lpos + 1L):(bestlCplusend - 1L)) {
      u[[i + 1L]] <- .cd_add(u[[i + 1L]], -upperl)
      m[[i + 1L]] <- .cd_sub_cd(m[[i + 1L]], mlu)
    }
  }
  cg$cover <- cover
  cplussize <- length(a)

  mudival <- .cd_val(.cd_div_cd(mu, al))
  eta <- ceiling(mudival)
  r <- .cd_sub_cd(mu, .cd_mul(al, floor(mudival)))
  if (.cd_val(r) < 0) r <- c(0, 0)
  ulminusetaplusone <- .cd_add(.cd_add(c(upperl, 0), -eta), 1)
  cplusthreshold <- .cd_mul_cd(ulminusetaplusone, al)
  kmin <- floor(eta - upperl - 0.5)
  alminusr <- .cd_sub_cd(al, r)
  alv <- .cd_val(al)

  phi_l <- function(av) {
    k <- min(trunc(av / alv), -1)
    while (k >= kmin) {
      kal <- .cd_mul(al, k)
      if (av >= .cd_val(.cd_add_cd(kal, r))) return(.cd_val(.hmip_cd_rsub(av, .cd_mul(r, k + 1))))
      if (av >= .cd_val(kal)) return(.cd_val(.cd_mul(alminusr, k)))
      k <- k - 1
    }
    .cd_val(.cd_mul(alminusr, kmin))
  }

  kmax <- floor(upperl - eta + 0.5)
  gamma_l <- function(z) {
    for (i in seq_len(cplussize)) {
      upperi <- trunc(upper[cover[i] + 1L])
      h <- 0
      while (h <= upperi) {
        mih <- .cd_add_cd(m[[i]], .cd_mul(a[[i]], h))
        uih <- .cd_add(u[[i]], h)
        mihplusdeltai <- .cd_sub_cd(.cd_add_cd(mih, a[[i]]), cplusthreshold)
        if (z <= .cd_val(mihplusdeltai))
          return(.cd_val(.cd_mul_cd(.cd_mul_cd(uih, ulminusetaplusone), alminusr)))
        k <- trunc(.cd_val(.cd_div_cd(.hmip_cd_rsub(z, mihplusdeltai), al))) - 1
        while (k <= kmax) {
          if (z <= .cd_val(.cd_add_cd(.cd_add_cd(mihplusdeltai, .cd_mul(al, k)), r)))
            return(.cd_val(.cd_mul_cd(.cd_add(.cd_mul_cd(uih, ulminusetaplusone), k), alminusr)))
          if (z <= .cd_val(.cd_add_cd(mihplusdeltai, .cd_mul(al, k + 1)))) {
            t <- .cd_add(.cd_mul_cd(.cd_mul_cd(uih, ulminusetaplusone), alminusr), z)
            t <- .cd_sub_cd(t, mih)
            t <- .cd_sub_cd(t, a[[i]])
            t <- .cd_add_cd(t, cplusthreshold)
            t <- .cd_sub_cd(t, .cd_mul(r, k + 1))
            return(.cd_val(t))
          }
          k <- k + 1
        }
        h <- h + 1
      }
    }
    mc <- m[[cplussize + 1L]]; uc <- u[[cplussize + 1L]]
    p <- trunc(.cd_val(.cd_div_cd(.hmip_cd_rsub(z, mc), al))) - 1
    repeat {
      if (z <= .cd_val(.cd_add_cd(.cd_add_cd(mc, .cd_mul(al, p)), r)))
        return(.cd_val(.cd_mul_cd(.cd_add(.cd_mul_cd(uc, ulminusetaplusone), p), alminusr)))
      if (z <= .cd_val(.cd_add_cd(mc, .cd_mul(al, p + 1)))) {
        t <- .cd_add(.cd_mul_cd(.cd_mul_cd(uc, ulminusetaplusone), alminusr), z)
        t <- .cd_sub_cd(t, mc)
        t <- .cd_sub_cd(t, .cd_mul(r, p + 1))
        return(.cd_val(t))
      }
      p <- p + 1
    }
  }

  rhs <- .cd_sub_cd(.cd_mul_cd(.cd_add(c(upperl, 0), -eta), r), lambda)
  integralSupport <- TRUE
  for (i in seq_len(rowlen)) {
    if (vals[i] == 0) next
    if (!cg$isintegral[i]) {
      if (vals[i] < 0) integralSupport <- FALSE else vals[i] <- 0
      next
    }
    if (coverflag[i]) {
      vals[i] <- -phi_l(-vals[i])
      rhs <- .cd_add(rhs, vals[i] * upper[i])
    } else {
      vals[i] <- gamma_l(vals[i])
    }
  }
  cg$vals <- vals
  cg$rhs <- rhs
  cg$integralSupport <- integralSupport
  cg$integralCoefficients <- FALSE
  TRUE
}

# updateViolationAndNorm: returns c(violation, norm)
.hmip_cg_viol_norm <- function(cg, index, aj, violation, norm) {
  sv <- cg$solval[index + 1L]
  violation <- violation + aj * sv
  if (aj > 0 && sv <= cg$feastol) return(c(violation, norm))
  if (aj < 0 && sv >= cg$upper[index + 1L] - cg$feastol) return(c(violation, norm))
  c(violation, norm + aj * aj)
}

.hmip_HighsCutGeneration_flipComplementation <- function(cg, index) {
  p <- index + 1L
  cg$complementation[p] <- 1L - cg$complementation[p]
  cg$solval[p] <- cg$upper[p] - cg$solval[p]
  cg$rhs <- .cd_add(cg$rhs, -(cg$upper[p] * cg$vals[p]))
  cg$vals[p] <- -cg$vals[p]
  invisible()
}

.hmip_HighsCutGeneration_removeComplementation <- function(cg) {
  if (length(cg$complementation) == 0L) return(invisible())
  for (i in seq_len(cg$rowlen) - 1L)
    if (cg$complementation[i + 1L]) .hmip_HighsCutGeneration_flipComplementation(cg, i)
  invisible()
}

# efficacy of the cmir cut with scale delta (the loop body the C++ repeats)
.hmip_cg_cmir_efficacy <- function(cg, delta, continuouscontribution, continuoussqrnorm) {
  scale <- 1 / delta
  scalrhs <- .cd_val(cg$rhs) * scale
  downrhs <- floor(scalrhs)
  f0 <- scalrhs - downrhs
  if (f0 < 0.005 || f0 > 0.995) return(NA_real_)
  oneoveroneminusf0 <- 1 / (1 - f0)
  if (oneoveroneminusf0 > 1e6) return(NA_real_)
  contscale <- scale * oneoveroneminusf0
  vn <- c(contscale * continuouscontribution - downrhs, contscale * contscale * continuoussqrnorm)
  for (j in cg$integerinds) {
    scalaj <- cg$vals[j + 1L] * scale
    downaj <- floor(scalaj + 1e-14)
    fj <- scalaj - downaj
    aj <- downaj + max(0, (fj - f0) * oneoveroneminusf0)
    vn <- .hmip_cg_viol_norm(cg, j, aj, vn[1L], vn[2L])
  }
  vn[1L] / sqrt(vn[2L])
}

.hmip_HighsCutGeneration_cmirCutGenerationHeuristic <- function(cg, minEfficacy, onlyInitialCMIRScale = FALSE) {
  cg$integralSupport <- FALSE
  cg$integralCoefficients <- FALSE
  feastol <- cg$feastol
  rowlen <- cg$rowlen
  continuouscontribution <- 0
  continuoussqrnorm <- 0
  deltas <- numeric()
  integerinds <- integer()
  maxabsdelta <- 0
  if (length(cg$complementation) < rowlen)
    cg$complementation <- c(cg$complementation, integer(rowlen - length(cg$complementation)))
  cg$complementation <- cg$complementation[seq_len(rowlen)]

  for (i in seq_len(rowlen) - 1L) {
    p <- i + 1L
    if (cg$isintegral[p]) {
      integerinds <- c(integerinds, i)
      if (cg$upper[p] < 2 * cg$solval[p]) .hmip_HighsCutGeneration_flipComplementation(cg, i)
      if (onlyInitialCMIRScale) next
      if (cg$solval[p] > feastol) {
        delta <- abs(cg$vals[p])
        if (delta <= 1e-4 || delta == maxabsdelta) next
        maxabsdelta <- max(maxabsdelta, delta)
        deltas <- c(deltas, delta)
      }
    } else {
      vn <- .hmip_cg_viol_norm(cg, i, cg$vals[p], continuouscontribution, continuoussqrnorm)
      continuouscontribution <- vn[1L]; continuoussqrnorm <- vn[2L]
    }
  }
  cg$integerinds <- integerinds

  if (continuoussqrnorm == 0 && length(deltas) > 1L) {
    intScale <- .hmip_integralScale(deltas, feastol, 1e-14)
    if (intScale != 0 && intScale <= 1e4) {
      scalrhs <- .cd_val(cg$rhs) * intScale
      downrhs <- floor(scalrhs)
      f0 <- scalrhs - downrhs
      if (f0 >= 0.005 && f0 <= 0.995) deltas <- c(deltas, 1 / intScale)
    }
  }
  deltas <- c(deltas, min(1, cg$initialScale))
  if (!onlyInitialCMIRScale) deltas <- c(deltas, maxabsdelta + min(1, cg$initialScale))

  deltas <- sort(deltas, method = "radix")
  curdelta <- deltas[1L]
  if (length(deltas) >= 2L) for (i in 2:length(deltas)) {
    if (deltas[i] - curdelta <= 10 * feastol) deltas[i] <- 0 else curdelta <- deltas[i]
  }
  deltas <- deltas[deltas != 0]
  cg$deltas <- deltas
  bestdelta <- -1
  bestefficacy <- minEfficacy

  for (delta in deltas) {
    efficacy <- .hmip_cg_cmir_efficacy(cg, delta, continuouscontribution, continuoussqrnorm)
    if (is.na(efficacy)) next
    if (efficacy > bestefficacy) { bestdelta <- delta; bestefficacy <- efficacy }
  }
  if (bestdelta == -1) return(FALSE)

  if (!onlyInitialCMIRScale) for (k in 1:3) {
    delta <- bestdelta * 2^k
    efficacy <- .hmip_cg_cmir_efficacy(cg, delta, continuouscontribution, continuoussqrnorm)
    if (is.na(efficacy)) next
    if (efficacy > bestefficacy) { bestdelta <- delta; bestefficacy <- efficacy }
  }

  for (k in integerinds) {
    if (cg$upper[k + 1L] == Inf) next
    if (cg$solval[k + 1L] <= feastol) next
    .hmip_HighsCutGeneration_flipComplementation(cg, k)
    efficacy <- .hmip_cg_cmir_efficacy(cg, bestdelta, continuouscontribution, continuoussqrnorm)
    if (is.na(efficacy)) {
      .hmip_HighsCutGeneration_flipComplementation(cg, k)
      next
    }
    if (efficacy > bestefficacy) bestefficacy <- efficacy
    else .hmip_HighsCutGeneration_flipComplementation(cg, k)
  }

  scale <- .cd_div_cd(c(1, 0), c(bestdelta, 0))
  scalrhs <- .cd_mul_cd(cg$rhs, scale)
  downrhs <- floor(.cd_val(scalrhs))
  f0 <- .cd_add(scalrhs, -downrhs)
  oneoveroneminusf0 <- .cd_div_cd(c(1, 0), .hmip_cd_rsub(1, f0))

  cg$rhs <- c(downrhs * bestdelta, 0)
  integralSupport <- TRUE
  vals <- cg$vals
  for (j in seq_len(rowlen)) {
    if (vals[j] == 0) next
    if (!cg$isintegral[j]) {
      if (vals[j] > 0) vals[j] <- 0
      else {
        vals[j] <- .cd_val(.cd_mul(oneoveroneminusf0, vals[j]))
        integralSupport <- FALSE
      }
    } else {
      scalaj <- .cd_mul(scale, vals[j])
      downaj <- floor(.cd_val(.cd_add(scalaj, 1e-14)))
      fj <- .cd_add(scalaj, -downaj)
      aj <- c(downaj, 0)
      if (.cd_val(fj) > .cd_val(f0)) aj <- .cd_add_cd(aj, .cd_mul_cd(.cd_sub_cd(fj, f0), oneoveroneminusf0))
      vals[j] <- .cd_val(.cd_mul(aj, bestdelta))
    }
  }
  cg$vals <- vals
  cg$integralSupport <- integralSupport
  cg$integralCoefficients <- FALSE
  TRUE
}

# scale(val): returns the scaling factor
.hmip_HighsCutGeneration_scale <- function(cg, val) {
  expshift <- min(10L, -.hmip_frexp_exp(val))
  f <- 2^expshift
  cg$rhs <- cg$rhs * f
  rl <- seq_len(cg$rowlen)
  cg$vals[rl] <- cg$vals[rl] * f
  f
}

.hmip_HighsCutGeneration_postprocessCut <- function(cg) {
  feastol <- cg$feastol; epsilon <- cg$epsilon
  rhsv <- .cd_val(cg$rhs)
  if (rhsv < 0 && rhsv > -epsilon) cg$rhs <- c(0, 0)

  removeZeros <- function() {
    i <- cg$rowlen
    while (i >= 1L) {
      if (cg$vals[i] == 0) {
        rl <- cg$rowlen
        cg$inds[i] <- cg$inds[rl]; cg$vals[i] <- cg$vals[rl]
        cg$rowlen <- rl - 1L
      }
      i <- i - 1L
    }
  }

  if (cg$integralSupport && cg$integralCoefficients) {
    removeZeros()
    return(TRUE)
  }

  globaldomain <- .hmip_cuts_mip(cg$lpRelaxation)$mipdata_$domain
  rl <- seq_len(cg$rowlen)
  maxAbsValue <- 0
  for (i in rl) maxAbsValue <- max(abs(cg$vals[i]), maxAbsValue)
  minCoefficientValue <- 100 * feastol * max(maxAbsValue, 1e-3)

  cg$integralSupport <- TRUE
  i <- cg$rowlen
  while (i >= 1L) {
    v <- cg$vals[i]
    if (v != 0) {
      if (abs(v) <= minCoefficientValue) {
        if (v < 0) {
          ub <- globaldomain$col_upper_[cg$inds[i] + 1L]
          if (ub == Inf) return(FALSE)
          cg$rhs <- .cd_add(cg$rhs, -(ub * v))
        } else {
          lb <- globaldomain$col_lower_[cg$inds[i] + 1L]
          if (lb == -Inf) return(FALSE)
          cg$rhs <- .cd_add(cg$rhs, -(lb * v))
        }
        cg$vals[i] <- 0
      } else if (cg$integralSupport && !.hmip_HighsLpRelaxation_isColIntegral(cg$lpRelaxation, cg$inds[i])) {
        cg$integralSupport <- FALSE
      }
    }
    i <- i - 1L
  }
  removeZeros()
  if (cg$rowlen == 0L) return(FALSE)

  if (cg$integralSupport) {
    rl <- seq_len(cg$rowlen)
    intscale <- .hmip_integralScale(cg$vals[rl], feastol, epsilon)
    scaleSmallestValToOne <- TRUE
    if (intscale != 0 && intscale * max(1, maxAbsValue) <= 2^52) {
      rhs <- .cd_two_sum(cg$rhs[1L], cg$rhs[2L])
      rhs <- .cd_mul(rhs, intscale)
      maxAbsValue <- .hmip_nearestInteger(maxAbsValue * intscale)
      for (i in rl) {
        scaleval <- .cd_mul(c(cg$vals[i], 0), intscale)
        intval <- .hmip_nearestInteger(.cd_val(scaleval))
        delta <- .cd_val(.cd_add(scaleval, -intval))
        cg$vals[i] <- intval
        if (delta < 0) {
          ub <- globaldomain$col_upper_[cg$inds[i] + 1L]
          if (ub == Inf) { cg$rhs <- rhs; return(FALSE) }
          rhs <- .cd_add(rhs, -(delta * ub))
        } else {
          lb <- globaldomain$col_lower_[cg$inds[i] + 1L]
          if (lb == -Inf) { cg$rhs <- rhs; return(FALSE) }
          rhs <- .cd_add(rhs, -(delta * lb))
        }
      }
      cg$rhs <- .cd_floor(.cd_add(rhs, feastol))
      if (intscale * maxAbsValue * feastol < 0.5) {
        scaleSmallestValToOne <- FALSE
        cg$integralCoefficients <- TRUE
      }
    }
    if (scaleSmallestValToOne) {
      minAbsValue <- Inf
      for (i in rl) minAbsValue <- min(abs(cg$vals[i]), minAbsValue)
      .hmip_HighsCutGeneration_scale(cg, minAbsValue - epsilon)
    }
  } else {
    .hmip_HighsCutGeneration_scale(cg, maxAbsValue - epsilon)
  }
  TRUE
}

# returns list(ok, hasUnboundedInts, hasGeneralInts, hasContinuous)
.hmip_HighsCutGeneration_preprocessBaseInequality <- function(cg) {
  feastol <- cg$feastol
  hasUnboundedInts <- FALSE; hasContinuous <- FALSE; hasGeneralInts <- FALSE
  res <- function(ok) list(ok = ok, hasUnboundedInts = hasUnboundedInts,
                           hasGeneralInts = hasGeneralInts, hasContinuous = hasContinuous)
  numZeros <- 0L
  maxact <- -feastol
  maxAbsVal <- 0
  rowlen <- cg$rowlen
  for (i in seq_len(rowlen)) maxAbsVal <- max(abs(cg$vals[i]), maxAbsVal)
  cg$initialScale <- .hmip_HighsCutGeneration_scale(cg, maxAbsVal)

  isintegral <- logical(rowlen)
  for (i in seq_len(rowlen)) {
    isintegral[i] <- .hmip_HighsLpRelaxation_isColIntegral(cg$lpRelaxation, cg$inds[i]) &&
      abs(cg$vals[i]) > 10 * feastol
    if (!isintegral[i]) {
      if (cg$upper[i] < 2 * cg$solval[i]) {
        if (length(cg$complementation) == 0L) cg$complementation <- integer(rowlen)
        .hmip_HighsCutGeneration_flipComplementation(cg, i - 1L)
      }
      v <- cg$vals[i]
      if (v > 0 || abs(v) * cg$upper[i] <= 10 * feastol) {
        if (v < 0) {
          if (cg$upper[i] == Inf) { cg$isintegral <- isintegral; return(res(FALSE)) }
          cg$rhs <- .cd_add(cg$rhs, -(v * cg$upper[i]))
        }
        numZeros <- numZeros + 1L
        cg$vals[i] <- 0
        next
      }
      hasContinuous <- TRUE
      if (v > 0) {
        if (cg$upper[i] == Inf) maxact <- Inf else maxact <- maxact + v * cg$upper[i]
      }
    } else {
      if (cg$upper[i] == Inf) {
        hasUnboundedInts <- TRUE
        hasGeneralInts <- TRUE
      } else if (cg$upper[i] != 1) {
        hasGeneralInts <- TRUE
      }
      if (cg$vals[i] > 0) maxact <- maxact + cg$vals[i] * cg$upper[i]
    }
  }
  cg$isintegral <- isintegral

  maxLen <- trunc(100 + 0.15 * .hmip_HighsLpRelaxation_numCols(cg$lpRelaxation))
  if (rowlen - numZeros > maxLen) {
    numCancel <- rowlen - numZeros - maxLen
    cancelNzs <- integer()
    for (i in seq_len(rowlen) - 1L) {
      p <- i + 1L
      cancelSlack <- if (cg$vals[p] > 0) cg$solval[p] else cg$upper[p] - cg$solval[p]
      if (cancelSlack <= feastol) cancelNzs <- c(cancelNzs, i)
    }
    if (length(cancelNzs) < numCancel) return(res(FALSE))
    if (length(cancelNzs) > numCancel) {
      vals <- cg$vals
      cancelNzs <- .hmip_std_partial_sort(cancelNzs, numCancel,
                                          function(a, b) abs(vals[a + 1L]) < abs(vals[b + 1L]))
    }
    for (i in seq_len(numCancel)) {
      p <- cancelNzs[i] + 1L
      if (cg$vals[p] < 0) cg$rhs <- .cd_add(cg$rhs, -(cg$vals[p] * cg$upper[p]))
      else maxact <- maxact - cg$vals[p] * cg$upper[p]
      cg$vals[p] <- 0
    }
    numZeros <- numZeros + numCancel
  }

  if (numZeros != 0L) {
    withCompl <- length(cg$complementation) != 0L
    i <- rowlen
    while (i >= 1L) {
      if (cg$vals[i] == 0) {
        rl <- cg$rowlen
        cg$inds[i] <- cg$inds[rl]; cg$vals[i] <- cg$vals[rl]
        cg$upper[i] <- cg$upper[rl]; cg$solval[i] <- cg$solval[rl]
        cg$isintegral[i] <- cg$isintegral[rl]
        if (withCompl) cg$complementation[i] <- cg$complementation[rl]
        cg$rowlen <- rl - 1L
        numZeros <- numZeros - 1L
        if (numZeros == 0L) break
      }
      i <- i - 1L
    }
  }
  res(maxact > .cd_val(cg$rhs))
}

# tryGenerateCut: on entry cg$vals/cg$inds are the caller vectors vals_ and
# inds_; on a TRUE return they hold the cut (the C++ vals_ and inds_).
.hmip_HighsCutGeneration_tryGenerateCut <- function(cg, hasUnboundedInts, hasGeneralInts, hasContinuous,
                                                    minEfficacy, onlyInitialCMIRScale = FALSE,
                                                    allowRejectCut = TRUE, lpSol = TRUE) {
  if (hasUnboundedInts)
    return(.hmip_HighsCutGeneration_cmirCutGenerationHeuristic(cg, minEfficacy, onlyInitialCMIRScale))

  rl <- seq_len(cg$rowlen)
  tmpVals <- cg$vals[rl]
  tmpInds <- cg$inds[rl]
  tmpComplementation <- cg$complementation
  tmpSolval <- cg$solval
  tmpRhs <- cg$rhs

  success <- FALSE
  saveIntegalSupport <- FALSE
  saveIntegralCoefficients <- FALSE
  if (.hmip_HighsCutGeneration_determineCover(cg, lpSol)) {
    if (!hasContinuous && !hasGeneralInts) {
      .hmip_HighsCutGeneration_separateLiftedKnapsackCover(cg)
      success <- TRUE
    } else if (hasGeneralInts) {
      success <- .hmip_HighsCutGeneration_separateLiftedMixedIntegerCover(cg)
    } else {
      success <- .hmip_HighsCutGeneration_separateLiftedMixedBinaryCover(cg)
    }
  }

  minMirEfficacy <- minEfficacy
  if (success) {
    saveIntegalSupport <- cg$integralSupport
    saveIntegralCoefficients <- cg$integralCoefficients
    vn <- c(-.cd_val(cg$rhs), 0)
    for (i in rl) vn <- .hmip_cg_viol_norm(cg, i - 1L, cg$vals[i], vn[1L], vn[2L])
    efficacy <- vn[1L] / sqrt(vn[2L])
    if (allowRejectCut && efficacy <= minEfficacy) {
      success <- FALSE
      cg$rhs <- tmpRhs
    } else {
      minMirEfficacy <- minMirEfficacy + efficacy
      t <- tmpRhs; tmpRhs <- cg$rhs; cg$rhs <- t
    }
  }

  # the lifted cut stays in the caller vectors; cmir works on the copies
  liftedVals <- cg$vals
  liftedInds <- cg$inds
  cg$vals <- tmpVals
  cg$inds <- tmpInds

  if (.hmip_HighsCutGeneration_cmirCutGenerationHeuristic(cg, minMirEfficacy, onlyInitialCMIRScale)) {
    TRUE
  } else if (success) {
    cg$rhs <- tmpRhs
    cg$complementation <- tmpComplementation
    cg$solval <- tmpSolval
    cg$vals <- liftedVals
    cg$inds <- liftedInds
    cg$integralSupport <- saveIntegalSupport
    cg$integralCoefficients <- saveIntegralCoefficients
    TRUE
  } else {
    FALSE
  }
}

.hmip_cg_remove_zeros <- function(cg) {
  i <- cg$rowlen
  while (i >= 1L) {
    if (cg$vals[i] == 0) {
      rl <- cg$rowlen
      cg$inds[i] <- cg$inds[rl]; cg$vals[i] <- cg$vals[rl]
      cg$rowlen <- rl - 1L
    }
    i <- i - 1L
  }
}

# violation check, coefficient tightening and addCut shared by generateCut
# and finalizeAndAddCut; returns list(success, inds, vals, rhs)
.hmip_cg_add_violated <- function(cg, inds_, vals_, rhs_) {
  mip <- .hmip_cuts_mip(cg$lpRelaxation)
  sol <- .hmip_cuts_sol(cg$lpRelaxation)$col_value
  violation <- c(-rhs_, 0)
  for (i in seq_along(inds_)) violation <- .cd_add(violation, sol[inds_[i] + 1L] * vals_[i])
  if (.cd_val(violation) <= 10 * cg$feastol) return(list(success = FALSE, inds = inds_, vals = vals_, rhs = rhs_))
  tc <- .hmip_HighsDomain_tightenCoefficients(mip$mipdata_$domain, inds_, vals_, length(inds_), rhs_)
  vals_ <- tc$vals; rhs_ <- tc$rhs
  r <- .hmip_HighsCutPool_addCut(cg$cutpool, mip, inds_, vals_, length(inds_), rhs_,
                                 cg$integralSupport && cg$integralCoefficients)
  list(success = r$index != -1L, inds = r$Rindex, vals = r$Rvalue, rhs = rhs_)
}

# generateCut(transLp, inds_, vals_, rhs_, onlyInitialCMIRScale): returns
# list(success, inds, vals, rhs)
.hmip_HighsCutGeneration_generateCut <- function(cg, transLp, inds_, vals_, rhs_, onlyInitialCMIRScale = FALSE) {
  fail <- function() list(success = FALSE, inds = cg$inds, vals = cg$vals, rhs = rhs_)
  tr <- .hmip_HighsTransformedLp_transform(transLp, vals_, cg$upper, cg$solval, inds_, rhs_, TRUE)
  cg$upper <- tr$upper; cg$solval <- tr$solval
  if (!tr$ok) return(list(success = FALSE, inds = tr$inds, vals = tr$vals, rhs = tr$rhs))
  intsPositive <- tr$integersPositive
  inds_ <- tr$inds; vals_ <- tr$vals; rhs_ <- tr$rhs

  cg$rowlen <- length(inds_)
  cg$inds <- inds_
  cg$vals <- vals_
  cg$rhs <- c(rhs_, 0)
  cg$complementation <- integer()
  pp <- .hmip_HighsCutGeneration_preprocessBaseInequality(cg)
  if (!pp$ok) return(fail())

  if (!pp$hasUnboundedInts && !intsPositive) {
    cg$complementation <- integer(cg$rowlen)
    for (i in seq_len(cg$rowlen)) {
      if (cg$vals[i] > 0 || !cg$isintegral[i]) next
      .hmip_HighsCutGeneration_flipComplementation(cg, i - 1L)
    }
  }

  if (!.hmip_HighsCutGeneration_tryGenerateCut(cg, pp$hasUnboundedInts, pp$hasGeneralInts, pp$hasContinuous,
                                               10 * cg$feastol, onlyInitialCMIRScale))
    return(fail())

  .hmip_HighsCutGeneration_removeComplementation(cg)
  .hmip_cg_remove_zeros(cg)

  rhs_ <- .cd_val(cg$rhs)
  rl <- seq_len(cg$rowlen)
  ut <- .hmip_HighsTransformedLp_untransform(transLp, cg$vals[rl], cg$inds[rl], rhs_)
  if (!ut$ok) return(list(success = FALSE, inds = ut$inds, vals = ut$vals, rhs = ut$rhs))
  inds_ <- ut$inds; vals_ <- ut$vals; rhs_ <- ut$rhs

  cg$rowlen <- length(inds_)
  cg$inds <- inds_
  cg$vals <- vals_
  cg$rhs <- c(rhs_, 0)
  if (!.hmip_HighsCutGeneration_postprocessCut(cg)) return(fail())
  rhs_ <- .cd_val(cg$rhs)
  rl <- seq_len(cg$rowlen)
  .hmip_cg_add_violated(cg, cg$inds[rl], cg$vals[rl], rhs_)
}

# generateConflict(localdomain, proofinds, proofvals, proofrhs): returns
# list(success, proofinds, proofvals, proofrhs)
.hmip_HighsCutGeneration_generateConflict <- function(cg, localdomain, proofinds, proofvals, proofrhs) {
  cg$inds <- proofinds
  cg$vals <- proofvals
  cg$rhs <- c(proofrhs, 0)
  rowlen <- length(proofinds)
  cg$rowlen <- rowlen
  fail <- function() list(success = FALSE, proofinds = cg$inds, proofvals = cg$vals, proofrhs = proofrhs)
  cg$complementation <- integer(rowlen)
  cg$upper <- numeric(rowlen)
  cg$solval <- numeric(rowlen)
  mip <- .hmip_cuts_mip(cg$lpRelaxation)
  globaldomain <- mip$mipdata_$domain
  activity <- 0
  for (i in seq_len(rowlen)) {
    col <- cg$inds[i] + 1L
    gl <- globaldomain$col_lower_[col]; gu <- globaldomain$col_upper_[col]
    cg$upper[i] <- gu - gl
    cg$solval[i] <- if (cg$vals[i] < 0) min(gu, localdomain$col_upper_[col]) else max(gl, localdomain$col_lower_[col])
    if (cg$vals[i] < 0 && gu != Inf) {
      cg$rhs <- .cd_add(cg$rhs, -(gu * cg$vals[i]))
      cg$vals[i] <- -cg$vals[i]
      cg$complementation[i] <- 1L
      cg$solval[i] <- gu - cg$solval[i]
    } else {
      cg$rhs <- .cd_add(cg$rhs, -(gl * cg$vals[i]))
      cg$complementation[i] <- 0L
      cg$solval[i] <- cg$solval[i] - gl
    }
    activity <- activity + cg$solval[i] * cg$vals[i]
  }
  if (activity > .cd_val(cg$rhs)) {
    solScale <- .cd_val(cg$rhs) / activity
    cg$solval <- cg$solval * solScale
  }

  pp <- .hmip_HighsCutGeneration_preprocessBaseInequality(cg)
  if (!pp$ok) return(fail())
  if (!.hmip_HighsCutGeneration_tryGenerateCut(cg, pp$hasUnboundedInts, pp$hasGeneralInts, pp$hasContinuous,
                                               cg$feastol, FALSE, FALSE, FALSE))
    return(fail())

  if (length(cg$complementation) != 0L) {
    for (i in seq_len(cg$rowlen)) {
      col <- cg$inds[i] + 1L
      if (cg$complementation[i]) {
        cg$rhs <- .cd_add(cg$rhs, -(globaldomain$col_upper_[col] * cg$vals[i]))
        cg$vals[i] <- -cg$vals[i]
      } else {
        cg$rhs <- .cd_add(cg$rhs, globaldomain$col_lower_[col] * cg$vals[i])
      }
    }
  }
  if (!.hmip_HighsCutGeneration_postprocessCut(cg)) return(fail())

  rl <- seq_len(cg$rowlen)
  proofvals <- cg$vals[rl]
  proofinds <- cg$inds[rl]
  proofrhs <- .cd_val(cg$rhs)
  cutintegral <- cg$integralSupport && cg$integralCoefficients
  tc <- .hmip_HighsDomain_tightenCoefficients(globaldomain, proofinds, proofvals, cg$rowlen, proofrhs)
  proofvals <- tc$vals; proofrhs <- tc$rhs
  r <- .hmip_HighsCutPool_addCut(cg$cutpool, mip, proofinds, proofvals, cg$rowlen, proofrhs,
                                 cutintegral, TRUE, TRUE, TRUE)
  list(success = r$index != -1L, proofinds = r$Rindex, proofvals = r$Rvalue, proofrhs = proofrhs)
}

# finalizeAndAddCut(inds_, vals_, rhs_): returns list(success, inds, vals, rhs)
.hmip_HighsCutGeneration_finalizeAndAddCut <- function(cg, inds_, vals_, rhs_) {
  cg$complementation <- integer()
  cg$rowlen <- length(inds_)
  cg$inds <- inds_
  cg$vals <- vals_
  cg$rhs <- c(rhs_, 0)
  cg$integralSupport <- TRUE
  cg$integralCoefficients <- FALSE
  i <- cg$rowlen
  while (i >= 1L) {
    if (cg$vals[i] == 0) {
      rl <- cg$rowlen
      cg$inds[i] <- cg$inds[rl]; cg$vals[i] <- cg$vals[rl]
      cg$rowlen <- rl - 1L
    } else {
      cg$integralSupport <- cg$integralSupport &&
        .hmip_HighsLpRelaxation_isColIntegral(cg$lpRelaxation, cg$inds[i])
    }
    i <- i - 1L
  }
  rl <- seq_len(cg$rowlen)
  if (!.hmip_HighsCutGeneration_postprocessCut(cg))
    return(list(success = FALSE, inds = cg$inds[rl], vals = cg$vals[rl], rhs = rhs_))
  rhs_ <- .cd_val(cg$rhs)
  rl <- seq_len(cg$rowlen)
  .hmip_cg_add_violated(cg, cg$inds[rl], cg$vals[rl], rhs_)
}

# ---------------------------------------------------------------------------
# HighsSeparator (mip/HighsSeparator.cpp) and the three separators
# HighsSeparation constructs: kind is "tableau", "path" or "modk". The timer
# clock index has no effect on results and is not kept.

.hmip_new_HighsSeparator <- function(mipsolver, kind) {
  s <- new.env(parent = emptyenv())
  s$kind <- kind
  s$numCutsFound <- 0L
  s$numCalls <- 0L
  if (kind == "tableau") s$numTries <- 0
  if (kind == "path") s$randgen <- .hs_rng(mipsolver$options_mip_$random_seed)
  s
}
.hmip_new_HighsTableauSeparator <- function(mipsolver) .hmip_new_HighsSeparator(mipsolver, "tableau")
.hmip_new_HighsPathSeparator <- function(mipsolver) .hmip_new_HighsSeparator(mipsolver, "path")
.hmip_new_HighsModkSeparator <- function(mipsolver) .hmip_new_HighsSeparator(mipsolver, "modk")
.hmip_HighsSeparator_getNumCutsFound <- function(s) s$numCutsFound
.hmip_HighsSeparator_getNumCalls <- function(s) s$numCalls

.hmip_HighsSeparator_run <- function(s, lpRelaxation, lpAggregator, transLp, cutpool) {
  s$numCalls <- s$numCalls + 1L
  currNumCuts <- .hmip_HighsCutPool_getNumCuts(cutpool)
  f <- switch(s$kind,
              tableau = .hmip_HighsTableauSeparator_separateLpSolution,
              path = .hmip_HighsPathSeparator_separateLpSolution,
              modk = .hmip_HighsModkSeparator_separateLpSolution)
  f(s, lpRelaxation, lpAggregator, transLp, cutpool)
  s$numCutsFound <- s$numCutsFound + .hmip_HighsCutPool_getNumCuts(cutpool) - currNumCuts
  invisible()
}

# ---------------------------------------------------------------------------
# HighsTableauSeparator (mip/HighsTableauSeparator.cpp)

.hmip_HighsTableauSeparator_separateLpSolution <- function(s, lpRelaxation, lpAggregator, transLp, cutpool) {
  lpSolver <- .hmip_HighsLpRelaxation_getLpSolver(lpRelaxation)
  if (!.hmip_Highs_hasInvert(lpSolver)) return(invisible())
  mip <- .hmip_cuts_mip(lpRelaxation)
  mipdata <- mip$mipdata_
  if (.hmip_HighsCutPool_getNumAvailableCuts(cutpool) > mip$options_mip_$mip_pool_soft_limit) return(invisible())

  basisinds <- .hmip_Highs_getBasicVariablesArray(lpSolver)
  numRow <- .hmip_HighsLpRelaxation_numRows(lpRelaxation)
  numCol <- .hmip_HighsLpRelaxation_numCols(lpRelaxation)
  cutGen <- .hmip_new_HighsCutGeneration(lpRelaxation, cutpool)
  lpSolution <- .hmip_cuts_sol(lpRelaxation)

  fb_index <- integer(); fb_frac <- numeric()
  for (i in seq_len(numRow) - 1L) {
    b <- basisinds[i + 1L]
    if (b >= numCol) {
      row <- b - numCol
      if (!.hmip_HighsLpRelaxation_isRowIntegral(lpRelaxation, row)) next
      my_fractionality <- .hmip_fractionality(lpSolution$row_value[row + 1L])
    } else {
      if (.hmip_variableType(mip, b) == 0L) next
      my_fractionality <- .hmip_fractionality(lpSolution$col_value[b + 1L])
    }
    if (my_fractionality < 1000 * mipdata$feastol) next
    fb_index <- c(fb_index, i); fb_frac <- c(fb_frac, my_fractionality)
  }
  if (length(fb_index) == 0L) return(invisible())

  maxTries <- 5000 + s$numCalls * 50 +
    (mipdata$total_lp_iterations - mipdata$heuristic_lp_iterations) %/% 10
  if (s$numTries >= maxTries) return(invisible())
  maxTries <- maxTries - s$numTries
  maxTries <- min(maxTries, 200 + trunc(0.1 * min(numRow, length(mipdata$integral_cols))))

  if (length(fb_index) > maxTries) {
    edgeWt <- .hmip_Highs_getDualEdgeWeights(lpSolver)
    key <- if (!is.null(edgeWt)) fb_frac * (1 - fb_frac) / edgeWt[fb_index + 1L] else fb_frac
    hh <- numeric(length(fb_index)); hl <- numeric(length(fb_index))
    for (q in seq_along(fb_index)) {
      h <- .hmip_cuts_hash_i64(s$numTries + fb_index[q])
      hh[q] <- .hs_hi32(h); hl[q] <- .hs_lo32(h)
    }
    # a strict total order: any correct sort gives the pdqsort result
    o <- order(key, hh, hl, decreasing = TRUE, method = "radix")[seq_len(maxTries)]
    fb_index <- fb_index[o]; fb_frac <- fb_frac[o]
  }

  s$numTries <- s$numTries + length(fb_index)
  nf <- length(fb_index)
  score <- rep(-1, nf)
  row_ep <- vector("list", nf)
  for (q in seq_len(nf)) {
    be <- .hmip_Highs_getBasisInverseRowSparse(lpSolver, fb_index[q])
    if (identical(be$ok, FALSE)) next
    if (be$count == 1L) next
    norm2 <- 0
    minWeight <- Inf; maxWeight <- 0
    rows <- integer(); weights <- numeric()
    for (j in seq_len(be$count)) {
      row <- be$index[j]
      weight <- be$array[row + 1L]
      scaledWeight <- .hmip_HighsLpRelaxation_getMaxAbsRowVal(lpRelaxation, row) * abs(weight)
      if (scaledWeight <= mipdata$feastol) next
      minWeight <- min(minWeight, scaledWeight)
      maxWeight <- max(maxWeight, scaledWeight)
      norm2 <- norm2 + scaledWeight * scaledWeight
      rows <- c(rows, row); weights <- c(weights, weight)
    }
    row_ep[[q]] <- list(rows = rows, weights = weights)
    if (length(rows) <= 1L) next
    if (maxWeight / minWeight <= 1e4) score[q] <- fb_frac[q] * (1 - fb_frac[q]) / norm2
  }
  keep <- which(!(score <= mipdata$feastol))
  if (length(keep) == 0L) return(invisible())
  score <- score[keep]; row_ep <- row_ep[keep]

  ord <- .hmip_pdqsort(seq_along(keep), function(a, b) score[a] > score[b], branchless = TRUE)
  bestScore <- -1
  numCuts <- .hmip_HighsCutPool_getNumCuts(cutpool)
  bestScoreFac <- c(0.0025, 0.01)
  mnumCol <- .hmip_numCol(mip)

  for (q in ord) {
    found <- .hmip_HighsCutPool_getNumCuts(cutpool) - numCuts
    if (found >= 1000) break
    if (score[q] < bestScoreFac[(found >= 50) + 1L] * bestScore) break
    rep_q <- row_ep[[q]]
    for (k in seq_along(rep_q$rows)) .hmip_HighsLpAggregator_addRow(lpAggregator, rep_q$rows[k], rep_q$weights[k])
    ag <- .hmip_HighsLpAggregator_getCurrentAggregation(lpAggregator, FALSE)
    nb <- length(ag$inds); ne <- length(rep_q$rows)
    # size_t arithmetic: a smaller aggregation wraps around and is skipped
    if (nb < ne || 10 * (nb - ne) > 10000 + mnumCol) {
      .hmip_HighsLpAggregator_clear(lpAggregator)
      next
    }
    if (nb > ne) {
      maxAbsVal <- 0; minAbsVal <- Inf
      for (i in seq_len(nb)) if (ag$inds[i] < mnumCol) {
        maxAbsVal <- max(abs(ag$vals[i]), maxAbsVal)
        minAbsVal <- min(abs(ag$vals[i]), minAbsVal)
      }
      if (maxAbsVal / minAbsVal > 1e6) {
        .hmip_HighsLpAggregator_clear(lpAggregator)
        next
      }
    }
    .hmip_HighsCutGeneration_generateCut(cutGen, transLp, ag$inds, ag$vals, 0)
    if (.hmip_HighsDomain_infeasible(mipdata$domain)) break
    ag <- .hmip_HighsLpAggregator_getCurrentAggregation(lpAggregator, TRUE)
    .hmip_HighsCutGeneration_generateCut(cutGen, transLp, ag$inds, ag$vals, 0)
    if (.hmip_HighsDomain_infeasible(mipdata$domain)) break
    .hmip_HighsLpAggregator_clear(lpAggregator)
    if (bestScore == -1 && .hmip_HighsCutPool_getNumCuts(cutpool) != numCuts) bestScore <- score[q]
  }
  invisible()
}

# ---------------------------------------------------------------------------
# HighsPathSeparator (mip/HighsPathSeparator.cpp)
# RowType: -2 kUnusuable, -1 kGeq, 0 kEq, 1 kLeq

.hmip_HighsPathSeparator_separateLpSolution <- function(s, lpRelaxation, lpAggregator, transLp, cutpool) {
  mip <- .hmip_cuts_mip(lpRelaxation)
  mipdata <- mip$mipdata_
  feastol <- mipdata$feastol
  lp <- .hmip_cuts_lp(lpRelaxation)
  lpSolution <- .hmip_cuts_sol(lpRelaxation)
  nrow <- lp$num_row_; ncol <- lp$num_col_
  Astart <- lp$a_matrix_$start_; Aindex <- lp$a_matrix_$index_; Avalue <- lp$a_matrix_$value_
  bdist <- transLp$boundDist

  rowtype <- integer(nrow)
  for (i in seq_len(nrow)) {
    if (lp$row_lower_[i] == lp$row_upper_[i]) { rowtype[i] <- 0L; next }
    lowerslack <- Inf; upperslack <- Inf
    if (lp$row_lower_[i] != -Inf) lowerslack <- lpSolution$row_value[i] - lp$row_lower_[i]
    if (lp$row_upper_[i] != Inf) upperslack <- lp$row_upper_[i] - lpSolution$row_value[i]
    if (lowerslack > feastol && upperslack > feastol) rowtype[i] <- -2L
    else if (lowerslack < upperslack) rowtype[i] <- -1L
    else rowtype[i] <- 1L
  }

  numContinuous <- integer(nrow)
  for (col in mipdata$continuous_cols) {
    if (bdist[col + 1L] == 0) next
    for (k in seq_len(Astart[col + 2L] - Astart[col + 1L])) {
      r <- Aindex[Astart[col + 1L] + k] + 1L
      numContinuous[r] <- numContinuous[r] + 1L
    }
  }

  colSubst_row <- rep(-1L, ncol); colSubst_val <- numeric(ncol)
  for (i in seq_len(nrow) - 1L) {
    if (rowtype[i + 1L] != 0L || numContinuous[i + 1L] != 1L) next
    rr <- .hmip_HighsLpRelaxation_getRow(lpRelaxation, i)
    col <- -1L; val <- 0
    for (j in seq_len(rr$len)) {
      if (.hmip_variableType(mip, rr$inds[j]) != 0L) next
      if (bdist[rr$inds[j] + 1L] == 0) next
      col <- rr$inds[j]; val <- rr$vals[j]
      break
    }
    if (colSubst_row[col + 1L] != -1L) next
    colSubst_row[col + 1L] <- i
    colSubst_val[col + 1L] <- val
    rowtype[i + 1L] <- -2L
  }

  inArc_row <- integer(); inArc_val <- numeric()
  outArc_row <- integer(); outArc_val <- numeric()
  colIn_first <- integer(ncol); colIn_second <- integer(ncol)
  colOut_first <- integer(ncol); colOut_second <- integer(ncol)
  for (col in mipdata$continuous_cols) {
    j <- col + 1L
    if (bdist[j] == 0) next
    if (colSubst_row[j] != -1L) next
    colIn_first[j] <- length(inArc_row)
    colOut_first[j] <- length(outArc_row)
    for (k in seq_len(Astart[j + 1L] - Astart[j])) {
      p <- Astart[j] + k
      r <- Aindex[p]; v <- Avalue[p]
      t <- rowtype[r + 1L]
      if (t == -2L) next
      if (t == 1L) {
        if (v < 0) { inArc_row <- c(inArc_row, r); inArc_val <- c(inArc_val, v) }
        else { outArc_row <- c(outArc_row, r); outArc_val <- c(outArc_val, v) }
      } else if (t == -1L) {
        if (v > 0) { inArc_row <- c(inArc_row, r); inArc_val <- c(inArc_val, v) }
        else { outArc_row <- c(outArc_row, r); outArc_val <- c(outArc_val, v) }
      } else {
        inArc_row <- c(inArc_row, r); inArc_val <- c(inArc_val, v)
        outArc_row <- c(outArc_row, r); outArc_val <- c(outArc_val, v)
      }
    }
    colIn_second[j] <- length(inArc_row)
    colOut_second[j] <- length(outArc_row)
  }
  arcs <- list(
    "in" = list(first = colIn_first, second = colIn_second, row = inArc_row, val = inArc_val),
    out = list(first = colOut_first, second = colOut_second, row = outArc_row, val = outArc_val))

  cutGen <- .hmip_new_HighsCutGeneration(lpRelaxation, cutpool)
  maxPathLen <- 6L
  randgen <- s$randgen
  maxWeight <- 1 / feastol
  minWeight <- feastol

  for (i in seq_len(nrow) - 1L) {
    t <- rowtype[i + 1L]
    if (t == -2L) next
    if (t == 0L) {
      scales <- if (lpSolution$row_dual[i + 1L] <= mipdata$epsilon) c(1, -1) else c(-1, 1)
    } else if (t == 1L) scales <- c(1, -1) else scales <- c(-1, 1)

    success <- FALSE
    for (sidx in 1:2) {
      .hmip_HighsLpAggregator_addRow(lpAggregator, i, scales[sidx])
      currentPath <- i
      tryNegatedScale <- FALSE
      isRowInCurrentPath <- function(row) any(currentPath == row)
      checkWeight <- function(w) { w <- abs(w); w <= maxWeight && w >= minWeight }
      skipCol <- function(col, a, o) {
        currPathLen <- length(currentPath)
        j <- col + 1L
        if (currPathLen == 1L && !tryNegatedScale) {
          if (a$second[j] - a$first[j] <= currPathLen) {
            k <- a$first[j]
            while (k < a$second[j]) {
              if (a$row[k + 1L] != i) { tryNegatedScale <<- TRUE; break }
              k <- k + 1L
            }
          } else tryNegatedScale <<- TRUE
        }
        if (o$first[j] == o$second[j]) return(TRUE)
        if (o$second[j] - o$first[j] <= currPathLen) {
          k <- o$first[j]
          while (k < o$second[j]) {
            if (!isRowInCurrentPath(o$row[k + 1L])) return(FALSE)
            k <- k + 1L
          }
          return(TRUE)
        }
        FALSE
      }
      # returns c(row, weight) or NULL
      findRow <- function(bestArcCol, val, a) {
        j <- bestArcCol + 1L
        arcRow <- .hmip_rand_range(randgen, a$first[j], a$second[j])
        try1 <- function(k) {
          r <- a$row[k + 1L]
          w <- -val / a$val[k + 1L]
          if (!isRowInCurrentPath(r) && checkWeight(w)) c(r, w) else NULL
        }
        res <- try1(arcRow)
        if (!is.null(res)) return(res)
        k <- arcRow + 1L
        while (k < a$second[j]) { res <- try1(k); if (!is.null(res)) return(res); k <- k + 1L }
        k <- a$first[j]
        while (k < arcRow) { res <- try1(k); if (!is.null(res)) return(res); k <- k + 1L }
        NULL
      }

      aggregatedPath <- list()
      while (length(currentPath) != maxPathLen) {
        ag <- .hmip_HighsLpAggregator_getCurrentAggregation(lpAggregator, FALSE)
        baseRowInds <- ag$inds; baseRowVals <- ag$vals
        addedSubstitutionRows <- FALSE
        bestOutArcCol <- -1L; outArcColVal <- 0; outArcColBoundDist <- 0
        bestInArcCol <- -1L; inArcColVal <- 0; inArcColBoundDist <- 0
        for (j in seq_along(baseRowInds)) {
          col <- baseRowInds[j]
          if (col >= ncol || bdist[col + 1L] == 0 || .hmip_HighsLpRelaxation_isColIntegral(lpRelaxation, col)) next
          if (colSubst_row[col + 1L] != -1L) {
            addedSubstitutionRows <- TRUE
            .hmip_HighsLpAggregator_addRow(lpAggregator, colSubst_row[col + 1L],
                                           -baseRowVals[j] / colSubst_val[col + 1L])
            next
          }
          if (addedSubstitutionRows) next
          if (baseRowVals[j] < 0) {
            if (skipCol(col, arcs$out, arcs[["in"]])) next
            if (bestOutArcCol == -1L || bdist[col + 1L] > outArcColBoundDist) {
              bestOutArcCol <- col; outArcColVal <- baseRowVals[j]; outArcColBoundDist <- bdist[col + 1L]
            }
          } else {
            if (skipCol(col, arcs[["in"]], arcs$out)) next
            if (bestInArcCol == -1L || bdist[col + 1L] > inArcColBoundDist) {
              bestInArcCol <- col; inArcColVal <- baseRowVals[j]; inArcColBoundDist <- bdist[col + 1L]
            }
          }
        }
        if (addedSubstitutionRows) next

        success <- .hmip_HighsCutGeneration_generateCut(cutGen, transLp, baseRowInds, baseRowVals, 0)$success
        ag <- .hmip_HighsLpAggregator_getCurrentAggregation(lpAggregator, TRUE)
        if (length(aggregatedPath) || bestOutArcCol != -1L || bestInArcCol != -1L)
          aggregatedPath[[length(aggregatedPath) + 1L]] <- ag
        success <- .hmip_HighsCutGeneration_generateCut(cutGen, transLp, ag$inds, ag$vals, 0)$success || success
        if (success || (bestOutArcCol == -1L && bestInArcCol == -1L)) break

        rw <- NULL
        if (bestInArcCol == -1L ||
            (bestOutArcCol != -1L && outArcColBoundDist >= inArcColBoundDist - feastol)) {
          rw <- findRow(bestOutArcCol, outArcColVal, arcs[["in"]])
          if (is.null(rw)) {
            if (bestInArcCol == -1L) break
            rw <- findRow(bestInArcCol, inArcColVal, arcs$out)
            if (is.null(rw)) break
          }
        } else {
          rw <- findRow(bestInArcCol, inArcColVal, arcs$out)
          if (is.null(rw)) break
        }
        currentPath <- c(currentPath, as.integer(rw[1L]))
        .hmip_HighsLpAggregator_addRow(lpAggregator, as.integer(rw[1L]), rw[2L])
      }

      pathLen <- length(aggregatedPath)
      if (pathLen > 1L) {
        success <- .hmip_path_mixing(lpRelaxation, transLp, cutGen, aggregatedPath, pathLen, nrow, ncol,
                                     feastol, mipdata$epsilon) || success
      }
      .hmip_HighsLpAggregator_clear(lpAggregator)
      if (!tryNegatedScale) break
    }
  }
  invisible()
}

# the path mixing cut of HighsPathSeparator::separateLpSolution; returns the
# success flag of finalizeAndAddCut (FALSE when no cut is tried)
.hmip_path_mixing <- function(lpRelaxation, transLp, cutGen, aggregatedPath, pathLen, nrow, ncol, feastol, epsilon) {
  success <- FALSE
  indexPos <- integer(ncol + nrow)
  inds <- integer(); solval <- numeric(); upper <- numeric(); isIntegral <- logical()
  rhs <- numeric(pathLen)
  delta <- 1
  for (k in seq_len(pathLen)) {
    ap <- aggregatedPath[[k]]
    tr <- .hmip_HighsTransformedLp_transform(transLp, ap$vals, numeric(), numeric(), ap$inds, 0, FALSE)
    aggregatedPath[[k]] <- list(inds = tr$inds, vals = tr$vals)
    rhs[k] <- tr$rhs
    if (!tr$ok) { pathLen <- k - 1L; break }
    if (k == 1L) {
      if (rhs[k] > 1e-14) { pathLen <- k - 1L; break }
      if (rhs[k] >= -feastol) rhs[k] <- 0
    } else if (rhs[k] >= rhs[k - 1L] - feastol) {
      pathLen <- k - 1L
      break
    }
    delta <- max(abs(rhs[k]), delta)
    ai <- tr$inds; av <- tr$vals
    for (j in seq_along(ai)) {
      index <- ai[j]
      pos <- indexPos[index + 1L]
      if (pos == 0L) {
        inds <- c(inds, index); solval <- c(solval, tr$solval[j]); upper <- c(upper, tr$upper[j])
        isIntegral <- c(isIntegral, .hmip_HighsLpRelaxation_isColIntegral(lpRelaxation, index))
        if (isIntegral[length(isIntegral)]) delta <- max(abs(av[j]), delta)
        indexPos[index + 1L] <- length(inds)
      } else if (isIntegral[pos]) {
        delta <- max(abs(av[j]), delta)
      }
    }
  }
  if (pathLen <= 1L) return(FALSE)

  delta <- .hmip_exp2_ceil_log2(delta + 1)
  numInds <- length(inds)
  cutRhs <- c(0, 0)
  cutVals <- numeric(numInds)
  maxFrac <- numeric(numInds)
  downSum <- matrix(0, 2L, numInds)
  fSum <- matrix(0, 2L, numInds)
  fLast <- 0
  scale <- -1 / delta
  for (k in seq_len(pathLen)) {
    f <- rhs[k] * scale
    fDiff <- .cd_add(c(f, 0), -fLast)
    ap <- aggregatedPath[[k]]
    cutRhs <- .cd_add_cd(cutRhs, fDiff)
    for (j in seq_along(ap$inds)) {
      i <- indexPos[ap$inds[j] + 1L]
      gj <- ap$vals[j] * scale
      if (!isIntegral[i]) {
        cutVals[i] <- max(cutVals[i], gj)
      } else {
        gjdown <- floor(gj)
        hj <- gj - gjdown
        maxFrac[i] <- max(maxFrac[i], hj)
        downSum[, i] <- .cd_add_cd(downSum[, i], .cd_mul(fDiff, gjdown))
        fSum[, i] <- .cd_add_cd(fSum[, i], fDiff)
        if (.cd_val(fSum[, i]) < maxFrac[i]) cutVals[i] <- .cd_val(.cd_add_cd(downSum[, i], fSum[, i]))
        else cutVals[i] <- .cd_val(.cd_add(downSum[, i], maxFrac[i]))
      }
    }
    if (k > 1L) {
      viol <- .cd_val(cutRhs)
      for (j in seq_len(numInds)) viol <- viol - solval[j] * cutVals[j]
      viol <- viol * delta
      if (viol > 10 * feastol) {
        scale <- -delta
        crhs <- .cd_val(cutRhs) * scale
        cutVals <- cutVals * scale
        j <- numInds
        while (j >= 1L) {
          if (abs(cutVals[j]) <= epsilon) {
            t <- cutVals[j]; cutVals[j] <- cutVals[numInds]; cutVals[numInds] <- t
            t <- inds[j]; inds[j] <- inds[numInds]; inds[numInds] <- t
            numInds <- numInds - 1L
          }
          j <- j - 1L
        }
        cutVals <- cutVals[seq_len(numInds)]
        inds <- inds[seq_len(numInds)]
        ut <- .hmip_HighsTransformedLp_untransform(transLp, cutVals, inds, crhs)
        if (ut$ok)
          success <- .hmip_HighsCutGeneration_finalizeAndAddCut(cutGen, ut$inds, ut$vals, ut$rhs)$success || success
        break
      }
    }
    fLast <- f
  }
  success
}

# ---------------------------------------------------------------------------
# HighsGFkSolve (mip/HighsGFkSolve.h/.cpp). Internally 1-based: rows,
# columns and nonzero positions are C++ value + 1 and 0 stands for -1, so
# the row splay trees use .hs_splay, .hs_splay_link, .hs_splay_unlink_at and
# .hs_find_nonzero (fields ARleft, ARright, Acol, rowroot).

# HighsGFk<k>::inverse(a) = a^(k - 2) mod k
.hmip_gfk_inverse <- function(a, k) {
  if (k == 2) return(1)
  r <- 1
  for (e in seq_len(k - 2)) r <- r * a
  r %% k
}

.hmip_new_HighsGFkSolve <- function() new.env(parent = emptyenv())

# fromCSC<k>(Aval, Aindex, Astart, numRow): Aindex and Astart 0-based
.hmip_HighsGFkSolve_fromCSC <- function(g, k, Aval, Aindex, Astart, numRow) {
  numCol <- length(Astart) - 1L
  g$k <- k; g$numCol <- numCol; g$numRow <- numRow
  g$colhead <- integer(numCol); g$colsize <- integer(numCol)
  g$rhs <- numeric(numRow); g$rowroot <- integer(numRow); g$rowsize <- integer(numRow)
  Avalue <- numeric(); Acol <- integer(); Arow <- integer()
  for (i in seq_len(numCol)) {
    for (j in seq_len(Astart[i + 1L] - Astart[i])) {
      p <- Astart[i] + j
      val <- Aval[p] %% k
      if (val == 0) next
      Avalue <- c(Avalue, val); Acol <- c(Acol, i); Arow <- c(Arow, Aindex[p] + 1L)
    }
  }
  nnz <- length(Avalue)
  g$Avalue <- Avalue; g$Acol <- Acol; g$Arow <- Arow
  g$Anext <- integer(nnz); g$Aprev <- integer(nnz)
  g$ARleft <- integer(nnz); g$ARright <- integer(nnz)
  g$freeslots <- integer()
  for (pos in seq_len(nnz)) .hmip_gfk_link(g, pos)
  invisible()
}
.hmip_HighsGFkSolve_setRhs <- function(g, row, val) {
  g$rhs[row + 1L] <- abs(val) %% g$k
  invisible()
}

.hmip_gfk_link <- function(g, pos) {
  col <- g$Acol[pos]
  g$Anext[pos] <- g$colhead[col]
  g$Aprev[pos] <- 0L
  g$colhead[col] <- pos
  if (g$Anext[pos] != 0L) g$Aprev[g$Anext[pos]] <- pos
  g$colsize[col] <- g$colsize[col] + 1L
  .hs_splay_link(g, pos, g$Arow[pos])
  g$rowsize[g$Arow[pos]] <- g$rowsize[g$Arow[pos]] + 1L
}
.hmip_gfk_unlink <- function(g, pos) {
  nxt <- g$Anext[pos]; prev <- g$Aprev[pos]
  if (nxt != 0L) g$Aprev[nxt] <- prev
  if (prev != 0L) g$Anext[prev] <- nxt else g$colhead[g$Acol[pos]] <- nxt
  g$colsize[g$Acol[pos]] <- g$colsize[g$Acol[pos]] - 1L
  row <- g$Arow[pos]
  .hs_splay_unlink_at(g, pos, g$rowroot[row], function(v) g$rowroot[row] <- v)
  g$rowsize[row] <- g$rowsize[row] - 1L
  g$Avalue[pos] <- 0
  g$freeslots <- c(g$freeslots, pos)
}
.hmip_gfk_addNonzero <- function(g, row, col, val) {
  if (length(g$freeslots) == 0L) {
    pos <- length(g$Avalue) + 1L
    g$Avalue[pos] <- val; g$Arow[pos] <- row; g$Acol[pos] <- col
    g$Anext[pos] <- 0L; g$Aprev[pos] <- 0L; g$ARleft[pos] <- 0L; g$ARright[pos] <- 0L
  } else {
    # priority_queue with std::greater: the smallest free position
    q <- which.min(g$freeslots)
    pos <- g$freeslots[q]
    g$freeslots <- g$freeslots[-q]
    g$Avalue[pos] <- val; g$Arow[pos] <- row; g$Acol[pos] <- col
    g$Aprev[pos] <- 0L
  }
  .hmip_gfk_link(g, pos)
}
.hmip_gfk_rowpositions <- function(g, pos) {
  out <- integer()
  if (pos == 0L) return(out)
  stack <- pos
  repeat {
    n <- length(stack)
    pos <- stack[n]; stack <- stack[-n]
    out <- c(out, pos)
    if (g$ARleft[pos] != 0L) stack <- c(stack, g$ARleft[pos])
    if (g$ARright[pos] != 0L) stack <- c(stack, g$ARright[pos])
    if (length(stack) == 0L) break
  }
  out
}

# std::priority_queue<pair<int,int>> with comp(a, b) = a.first > b.first
# (libstdc++ push_heap / pop_heap); h is list(first, second)
.hmip_gfk_heap_push <- function(h, a, b) {
  hole <- length(h$first)
  h$first[hole + 1L] <- a; h$second[hole + 1L] <- b
  parent <- (hole - 1L) %/% 2L
  while (hole > 0L && h$first[parent + 1L] > a) {
    h$first[hole + 1L] <- h$first[parent + 1L]; h$second[hole + 1L] <- h$second[parent + 1L]
    hole <- parent
    parent <- (hole - 1L) %/% 2L
  }
  h$first[hole + 1L] <- a; h$second[hole + 1L] <- b
  h
}
.hmip_gfk_heap_pop <- function(h) {
  n <- length(h$first)
  if (n > 1L) {
    vf <- h$first[n]; vs <- h$second[n]
    h$first[n] <- h$first[1L]; h$second[n] <- h$second[1L]
    len <- n - 1L
    hole <- 0L; second <- 0L
    while (second < (len - 1L) %/% 2L) {
      second <- 2L * (second + 1L)
      if (h$first[second + 1L] > h$first[second]) second <- second - 1L
      h$first[hole + 1L] <- h$first[second + 1L]; h$second[hole + 1L] <- h$second[second + 1L]
      hole <- second
    }
    if (len %% 2L == 0L && second == (len - 2L) %/% 2L) {
      second <- 2L * (second + 1L)
      h$first[hole + 1L] <- h$first[second]; h$second[hole + 1L] <- h$second[second]
      hole <- second - 1L
    }
    parent <- (hole - 1L) %/% 2L
    while (hole > 0L && h$first[parent + 1L] > vf) {
      h$first[hole + 1L] <- h$first[parent + 1L]; h$second[hole + 1L] <- h$second[parent + 1L]
      hole <- parent
      parent <- (hole - 1L) %/% 2L
    }
    h$first[hole + 1L] <- vf; h$second[hole + 1L] <- vs
  }
  h$first <- h$first[-n]; h$second <- h$second[-n]
  h
}

# solve<k>(reportSolution): reportSolution(weights) gets list(index (0-based
# system columns), weight)
.hmip_HighsGFkSolve_solve <- function(g, reportSolution) {
  k <- g$k
  numCol <- g$numCol; numRow <- g$numRow
  pq <- list(first = integer(), second = integer())
  for (i in seq_len(numCol)) pq <- .hmip_gfk_heap_push(pq, g$colsize[i], i)
  maxPivot <- min(numRow, numCol)
  factorColPerm <- integer(); factorRowPerm <- integer()
  colBasisStatus <- integer(numCol)
  rowUsed <- logical(numRow)
  numPivot <- 0L

  while (length(pq$first)) {
    oldColSize <- pq$first[1L]; pivotCol <- pq$second[1L]
    pq <- .hmip_gfk_heap_pop(pq)
    if (g$colsize[pivotCol] == 0L) next
    if (g$colsize[pivotCol] != oldColSize) {
      pq <- .hmip_gfk_heap_push(pq, g$colsize[pivotCol], pivotCol)
      next
    }
    pivot <- 0L; pivotRow <- 0L; pivotRowLen <- .Machine$integer.max
    coliter <- g$colhead[pivotCol]
    while (coliter != 0L) {
      row <- g$Arow[coliter]
      if (!rowUsed[row] && g$rowsize[row] < pivotRowLen) {
        pivotRowLen <- g$rowsize[row]; pivotRow <- row; pivot <- coliter
      }
      coliter <- g$Anext[coliter]
    }
    pivotInverse <- .hmip_gfk_inverse(g$Avalue[pivot], k)
    rowpositions <- .hmip_gfk_rowpositions(g, g$rowroot[pivotRow])
    rowposColsizes <- g$colsize[g$Acol[rowpositions]]

    coliter <- g$colhead[pivotCol]
    while (coliter != 0L) {
      nxt <- g$Anext[coliter]
      if (coliter != pivot) {
        row <- g$Arow[coliter]
        if (!rowUsed[row]) {
          pivotRowScale <- pivotInverse * (k - g$Avalue[coliter])
          g$rhs[row] <- (g$rhs[row] + pivotRowScale * g$rhs[pivotRow]) %% k
          for (prp in rowpositions) {
            nonzeroPos <- .hs_find_nonzero(g, g$Arow[coliter], g$Acol[prp])
            if (nonzeroPos == 0L) {
              val <- (pivotRowScale * g$Avalue[prp]) %% k
              if (val != 0) .hmip_gfk_addNonzero(g, row, g$Acol[prp], val)
            } else {
              g$Avalue[nonzeroPos] <- (g$Avalue[nonzeroPos] + pivotRowScale * g$Avalue[prp]) %% k
              if (g$Avalue[nonzeroPos] == 0) .hmip_gfk_unlink(g, nonzeroPos)
            }
          }
        }
      }
      coliter <- nxt
    }

    numPivot <- numPivot + 1L
    factorColPerm <- c(factorColPerm, pivotCol)
    factorRowPerm <- c(factorRowPerm, pivotRow)
    colBasisStatus[pivotCol] <- 1L
    rowUsed[pivotRow] <- TRUE
    if (numPivot == maxPivot) break

    for (q in seq_along(rowpositions)) {
      col <- g$Acol[rowpositions[q]]
      oldsize <- rowposColsizes[q]
      g$colsize[col] <- g$colsize[col] - 1L
      if (g$colsize[col] == 0L) next
      if (g$colsize[col] < oldsize) pq <- .hmip_gfk_heap_push(pq, g$colsize[col], col)
    }
  }

  for (i in seq_len(numRow)) if (!rowUsed[i] && g$rhs[i] != 0) return(invisible())

  numFactorRows <- length(factorRowPerm)
  bs_index <- integer(); bs_col <- integer()
  i <- numFactorRows
  while (i >= 1L) {
    row <- factorRowPerm[i]
    stack <- g$rowroot[row]
    while (length(stack)) {
      n <- length(stack)
      rowpos <- stack[n]; stack <- stack[-n]
      if (g$ARleft[rowpos] != 0L) stack <- c(stack, g$ARleft[rowpos])
      if (g$ARright[rowpos] != 0L) stack <- c(stack, g$ARright[rowpos])
      col <- g$Acol[rowpos]
      if (colBasisStatus[col] != 0L) next
      colBasisStatus[col] <- -1L
      bs_index <- c(bs_index, i); bs_col <- c(bs_col, col)
    }
    i <- i - 1L
  }

  basisSwapPos <- 0L
  repeat {
    performedBasisSwap <- FALSE
    sol_index <- integer(); sol_weight <- numeric()
    i <- numFactorRows
    while (i >= 1L) {
      row <- factorRowPerm[i]
      solval <- 0
      for (q in seq_along(sol_index)) {
        pos <- .hs_find_nonzero(g, row, sol_index[q])
        if (pos != 0L) solval <- solval + g$Avalue[pos] * sol_weight[q]
      }
      solval <- g$rhs[row] + k - (solval %% k)
      col <- factorColPerm[i]
      pos <- .hs_find_nonzero(g, row, col)
      colValInverse <- .hmip_gfk_inverse(g$Avalue[pos], k)
      solval <- (solval * colValInverse) %% k
      if (solval != 0) { sol_index <- c(sol_index, col); sol_weight <- c(sol_weight, solval) }
      i <- i - 1L
    }
    reportSolution(list(index = sol_index - 1L, weight = sol_weight))

    if (basisSwapPos < length(bs_index)) {
      basisIndex <- bs_index[basisSwapPos + 1L]
      enteringCol <- bs_col[basisSwapPos + 1L]
      leavingCol <- factorColPerm[basisIndex]
      factorColPerm[basisIndex] <- enteringCol
      colBasisStatus[enteringCol] <- 1L
      colBasisStatus[leavingCol] <- 0L
      performedBasisSwap <- TRUE
      basisSwapPos <- basisSwapPos + 1L
    }
    if (!performedBasisSwap) break
  }
  invisible()
}

# ---------------------------------------------------------------------------
# HighsModkSeparator (mip/HighsModkSeparator.cpp)

.hmip_HighsModkSeparator_separateLpSolution <- function(s, lpRelaxation, lpAggregator, transLp, cutpool) {
  mip <- .hmip_cuts_mip(lpRelaxation)
  mipdata <- mip$mipdata_
  feastol <- mipdata$feastol; epsilon <- mipdata$epsilon
  lp <- .hmip_cuts_lp(lpRelaxation)
  nrow <- lp$num_row_; ncol <- lp$num_col_
  Astart <- lp$a_matrix_$start_; Aindex <- lp$a_matrix_$index_
  skipRow <- logical(nrow)
  for (col in mipdata$continuous_cols) {
    if (transLp$boundDist[col + 1L] == 0) next
    for (q in seq_len(Astart[col + 2L] - Astart[col + 1L])) skipRow[Aindex[Astart[col + 1L] + q] + 1L] <- TRUE
  }
  cutGen <- .hmip_new_HighsCutGeneration(lpRelaxation, cutpool)
  is_row <- integer(); is_scale <- numeric()
  intSystemValue <- numeric(); intSystemIndex <- integer(); intSystemStart <- 0L
  lpSolution <- .hmip_cuts_sol(lpRelaxation)
  numNonzeroRhs <- 0L
  maxIntRowLen <- trunc(1000 + 0.1 * ncol)
  isCont <- function(col) .hmip_variableType(mip, col) == 0L

  for (row in seq_len(nrow) - 1L) {
    p <- row + 1L
    if (skipRow[p]) next
    if (lp$row_upper_[p] - lpSolution$row_value[p] <= feastol) leqRow <- TRUE
    else if (lpSolution$row_value[p] - lp$row_lower_[p] <= feastol) leqRow <- FALSE
    else next
    rr <- .hmip_HighsLpRelaxation_getRow(lpRelaxation, row)
    inds <- rr$inds[seq_len(rr$len)]
    if (leqRow) { rhs <- lp$row_upper_[p]; vals <- rr$vals[seq_len(rr$len)] }
    else { rhs <- -lp$row_lower_[p]; vals <- -rr$vals[seq_len(rr$len)] }

    tr <- .hmip_HighsTransformedLp_transform(transLp, vals, numeric(), numeric(), inds, rhs, FALSE, TRUE)
    if (!tr$ok) next
    inds <- tr$inds; vals <- tr$vals; rhs <- tr$rhs; solval <- tr$solval
    rowlen <- length(inds)
    if (rowlen > maxIntRowLen) {
      intRowLen <- 0L
      for (i in seq_len(rowlen)) {
        if (solval[i] <= feastol) next
        if (isCont(inds[i])) next
        intRowLen <- intRowLen + 1L
      }
      if (intRowLen > maxIntRowLen || (intRowLen == 0L && abs(rhs) <= epsilon)) next
    }

    if (!.hmip_HighsLpRelaxation_isRowIntegral(lpRelaxation, row)) {
      scaleVals <- numeric()
      for (i in seq_len(rowlen)) {
        if (isCont(inds[i])) next
        if (solval[i] > feastol) scaleVals <- c(scaleVals, vals[i])
      }
      if (abs(rhs) > epsilon) scaleVals <- c(scaleVals, -rhs)
      if (length(scaleVals) == 0L) next
      intscale <- .hmip_integralScale(scaleVals, feastol, epsilon)
      if (intscale == 0 || intscale > 1e6) next
      intrhs <- .hmip_nearestInteger(intscale * rhs)
      for (i in seq_len(rowlen)) {
        if (isCont(inds[i])) next
        if (solval[i] > feastol) {
          intSystemIndex <- c(intSystemIndex, inds[i])
          intSystemValue <- c(intSystemValue, .hmip_nearestInteger(intscale * vals[i]))
        }
      }
    } else {
      intscale <- 1
      intrhs <- .hmip_nearestInteger(rhs)
      for (i in seq_len(rowlen)) {
        if (solval[i] > feastol) {
          intSystemIndex <- c(intSystemIndex, inds[i])
          intSystemValue <- c(intSystemValue, .hmip_nearestInteger(vals[i]))
        }
      }
    }
    numNonzeroRhs <- numNonzeroRhs + (intrhs != 0)
    intSystemIndex <- c(intSystemIndex, ncol)
    intSystemValue <- c(intSystemValue, intrhs)
    intSystemStart <- c(intSystemStart, length(intSystemValue))
    is_row <- c(is_row, row); is_scale <- c(is_scale, intscale)
  }
  if (length(is_row) == 0L || numNonzeroRhs == 0L) return(invisible())

  usedWeights <- new.env(parent = emptyenv())
  k <- 2
  foundCut <- function(w) {
    if (length(w$index) == 0L) return()
    o <- order(w$index, method = "radix")
    w$index <- w$index[o]; w$weight <- w$weight[o]
    key <- paste(c(w$index, "|", w$weight), collapse = ",")
    if (!is.null(usedWeights[[key]])) return()
    assign(key, TRUE, envir = usedWeights)
    for (q in seq_along(w$index)) {
      weight <- is_scale[w$index[q] + 1L] * (((w$weight[q] * (k - 1)) %% k) / k)
      .hmip_HighsLpAggregator_addRow(lpAggregator, is_row[w$index[q] + 1L], weight)
    }
    ag <- .hmip_HighsLpAggregator_getCurrentAggregation(lpAggregator, FALSE)
    .hmip_HighsCutGeneration_generateCut(cutGen, transLp, ag$inds, ag$vals, 0, TRUE)
    if (k != 2) {
      .hmip_HighsLpAggregator_clear(lpAggregator)
      for (q in seq_along(w$index)) {
        weight <- is_scale[w$index[q] + 1L] * (w$weight[q] / k)
        .hmip_HighsLpAggregator_addRow(lpAggregator, is_row[w$index[q] + 1L], weight)
      }
    }
    ag <- .hmip_HighsLpAggregator_getCurrentAggregation(lpAggregator, TRUE)
    .hmip_HighsCutGeneration_generateCut(cutGen, transLp, ag$inds, ag$vals, 0, TRUE)
    .hmip_HighsLpAggregator_clear(lpAggregator)
  }

  for (kk in c(2, 3, 5, 7)) {
    rm(list = ls(usedWeights), envir = usedWeights)
    k <- kk
    numCuts <- .hmip_HighsCutPool_getNumCuts(cutpool)
    g <- .hmip_new_HighsGFkSolve()
    .hmip_HighsGFkSolve_fromCSC(g, kk, intSystemValue, intSystemIndex, intSystemStart, ncol + 1L)
    .hmip_HighsGFkSolve_setRhs(g, ncol, 1)
    .hmip_HighsGFkSolve_solve(g, foundCut)
    if (.hmip_HighsCutPool_getNumCuts(cutpool) != numCuts) break
  }
  invisible()
}

# ---------------------------------------------------------------------------
# HighsSeparation (mip/HighsSeparation.cpp)

.hmip_new_HighsSeparation <- function(mipsolver) {
  sep <- new.env(parent = emptyenv())
  sep$separators <- list(.hmip_new_HighsTableauSeparator(mipsolver),
                         .hmip_new_HighsPathSeparator(mipsolver),
                         .hmip_new_HighsModkSeparator(mipsolver))
  sep$cutset <- .hmip_new_HighsCutSet()
  sep$lp <- NULL
  sep
}
.hmip_HighsSeparation_setLpRelaxation <- function(sep, lp) {
  sep$lp <- lp
  invisible()
}

# separationRound(propdomain, status): returns list(ncuts, status)
.hmip_HighsSeparation_separationRound <- function(sep, propdomain, status) {
  lp <- sep$lp
  mipdata <- .hmip_cuts_mip(lp)$mipdata_
  isGlobal <- identical(propdomain, mipdata$domain)
  redcost <- function(st) {
    if (isGlobal && .hmip_HighsLpRelaxation_unscaledDualFeasible(st)) {
      .hmip_HighsRedcostFixing_addRootRedcost(mipdata$redcostfixing, mipdata$mipsolver,
                                              .hmip_cuts_sol(lp)$col_dual, .hmip_HighsLpRelaxation_getObjective(lp))
      if (mipdata$upper_limit != Inf)
        .hmip_HighsRedcostFixing_propagateRootRedcost(mipdata$redcostfixing, mipdata$mipsolver)
    }
  }
  propagateAndResolve <- function() {
    if (.hmip_HighsDomain_infeasible(propdomain) || .hmip_HighsDomain_infeasible(mipdata$domain)) {
      status <<- "Infeasible"
      .hmip_HighsDomain_clearChangedCols(propdomain)
      return(-1L)
    }
    .hmip_HighsDomain_propagate(propdomain)
    if (.hmip_HighsDomain_infeasible(propdomain)) {
      status <<- "Infeasible"
      .hmip_HighsDomain_clearChangedCols(propdomain)
      return(-1L)
    }
    .hmip_HighsCliqueTable_cleanupFixed(mipdata$cliquetable, mipdata$domain)
    if (.hmip_HighsDomain_infeasible(mipdata$domain)) {
      status <<- "Infeasible"
      .hmip_HighsDomain_clearChangedCols(propdomain)
      return(-1L)
    }
    numBoundChgs <- length(.hmip_HighsDomain_getChangedCols(propdomain))
    while (length(.hmip_HighsDomain_getChangedCols(propdomain))) {
      .hmip_HighsLpRelaxation_setObjectiveLimit(lp, mipdata$upper_limit)
      status <<- .hmip_HighsLpRelaxation_resolveLp(lp, propdomain)
      if (!.hmip_HighsLpRelaxation_scaledOptimal(status)) return(-1L)
      redcost(status)
    }
    numBoundChgs
  }
  done <- function() list(ncuts = 0L, status = status)

  .hmip_HighsImplications_separateImpliedBounds(mipdata$implications, lp, .hmip_cuts_sol(lp)$col_value,
                                                mipdata$cutpool, mipdata$feastol)
  ncuts <- 0L
  nb <- propagateAndResolve()
  if (nb == -1L) return(done())
  ncuts <- ncuts + nb

  .hmip_HighsCliqueTable_separateCliques(mipdata$cliquetable, mipdata$mipsolver, .hmip_cuts_sol(lp)$col_value,
                                         mipdata$cutpool, mipdata$feastol)
  nb <- propagateAndResolve()
  if (nb == -1L) return(done())
  ncuts <- ncuts + nb

  if (!isGlobal) .hmip_HighsLpRelaxation_computeBasicDegenerateDuals(lp, mipdata$feastol, propdomain)

  transLp <- .hmip_new_HighsTransformedLp(lp, mipdata$implications)
  if (.hmip_HighsDomain_infeasible(mipdata$domain)) {
    status <- "Infeasible"
    return(done())
  }
  lpAggregator <- .hmip_new_HighsLpAggregator(lp)
  for (separator in sep$separators) {
    .hmip_HighsSeparator_run(separator, lp, lpAggregator, transLp, mipdata$cutpool)
    if (.hmip_HighsDomain_infeasible(mipdata$domain)) {
      status <- "Infeasible"
      return(done())
    }
  }

  nb <- propagateAndResolve()
  if (nb == -1L) return(done())
  ncuts <- ncuts + nb

  .hmip_HighsCutPool_separate(mipdata$cutpool, .hmip_cuts_sol(lp)$col_value, propdomain, sep$cutset, mipdata$feastol)
  nc <- .hmip_HighsCutSet_numCuts(sep$cutset)
  if (nc > 0L) {
    ncuts <- ncuts + nc
    .hmip_HighsLpRelaxation_addCuts(lp, sep$cutset)
    status <- .hmip_HighsLpRelaxation_resolveLp(lp, propdomain)
    .hmip_HighsLpRelaxation_performAging(lp, TRUE)
    redcost(status)
  }
  list(ncuts = ncuts, status = status)
}

.hmip_HighsSeparation_separate <- function(sep, propdomain) {
  lp <- sep$lp
  status <- .hmip_HighsLpRelaxation_getStatus(lp)
  mipdata <- .hmip_cuts_mip(lp)$mipdata_
  nfrac <- function() length(.hmip_cuts_frac_cols(.hmip_HighsLpRelaxation_getFractionalIntegers(lp)))
  if (.hmip_HighsLpRelaxation_scaledOptimal(status) && nfrac() > 0L) {
    firstobj <- mipdata$rootlpsolobj
    while (.hmip_HighsLpRelaxation_getObjective(lp) < mipdata$optimality_limit) {
      lastobj <- .hmip_HighsLpRelaxation_getObjective(lp)
      it0 <- .hmip_HighsLpRelaxation_getNumLpIterations(lp)
      r <- .hmip_HighsSeparation_separationRound(sep, propdomain, status)
      status <- r$status
      nlpiters <- .hmip_HighsLpRelaxation_getNumLpIterations(lp) - it0
      mipdata$sepa_lp_iterations <- mipdata$sepa_lp_iterations + nlpiters
      mipdata$total_lp_iterations <- mipdata$total_lp_iterations + nlpiters
      if (r$ncuts == 0L || !.hmip_HighsLpRelaxation_scaledOptimal(status) || nfrac() == 0L) break
      if ((.hmip_HighsLpRelaxation_getObjective(lp) - firstobj) <= max(lastobj - firstobj, mipdata$feastol) * 1.01)
        break
    }
  } else {
    .hmip_HighsLpRelaxation_performAging(lp, TRUE)
    .hmip_HighsCutPool_performAging(mipdata$cutpool)
  }
  invisible()
}

# ---------------------------------------------------------------------------
# HighsCliqueTable::separateCliques(mipsolver, sol, cutpool, feastol)
# (HighsCliqueTable.cpp:1687). CliqueVar values are 2 * col + val. The
# BronKerboschData env carries the C++ member names; bronKerboschRecurse,
# runCliqueSubsumption and doAddClique are group A's.


.hmip_HighsCliqueTable_separateCliques <- function(ct, mipsolver, sol, cutpool, feastol) {
  data <- .hmip_new_BronKerboschData(sol)
  data$feastol <- feastol
  mipdata <- mipsolver$mipdata_
  data$maxNeighbourhoodQueries <- 1000000 + 100 * .hmip_numNonzero(mipsolver) + mipdata$total_lp_iterations * 1000
  if (ct$numNeighbourhoodQueries > data$maxNeighbourhoodQueries) return(invisible())
  data$maxNeighbourhoodQueries <- data$maxNeighbourhoodQueries - ct$numNeighbourhoodQueries
  globaldom <- mipdata$domain

  P <- integer()
  for (i in mipdata$integral_cols) {
    if (ct$colsubstituted[i + 1L] || ct$colDeleted[i + 1L]) next
    v0 <- 2L * i; v1 <- 2L * i + 1L
    if (ct$numcliquesvar[v0 + 1L] != 0L && .hmip_cv_weight(v0, sol) > feastol) P <- c(P, v0)
    if (ct$numcliquesvar[v1 + 1L] != 0L && .hmip_cv_weight(v1, sol) > feastol) P <- c(P, v1)
  }
  data$P <- P
  .hmip_HighsCliqueTable_bronKerboschRecurse(ct, data, length(P), integer(), 0L)

  runcliquesubsumption <- FALSE
  for (clique in data$cliques) {
    rhs <- 1
    runcliquesubsumption <- .hmip_HighsCliqueTable_cliquesSize(ct) > 2L
    inds <- clique %/% 2L
    vals <- numeric(length(clique))
    for (q in seq_along(clique)) {
      if (clique[q] %% 2L == 0L) { vals[q] <- -1; rhs <- rhs - 1 } else vals[q] <- 1
    }
    rhs <- floor(rhs + 0.5)
    .hmip_HighsCutPool_addCut(cutpool, mipsolver, inds, vals, length(inds), rhs, TRUE, FALSE, FALSE)
  }
  ct$numNeighbourhoodQueries <- ct$numNeighbourhoodQueries + data$numNeighbourhoodQueries

  if (runcliquesubsumption) {
    nc <- .hmip_HighsCliqueTable_cliquesSize(ct)
    if (length(ct$cliquehits) < nc) ct$cliquehits <- c(ct$cliquehits, integer(nc - length(ct$cliquehits)))
    for (clique in data$cliques) {
      rs <- .hmip_HighsCliqueTable_runCliqueSubsumption(ct, globaldom, clique)
      clique <- rs$clique
      if (length(clique) == 0L) next
      if (rs$nremoved != 0L) .hmip_HighsCliqueTable_doAddClique(ct, clique, length(clique), FALSE, -1L)
    }
  }
  invisible()
}

# ---------------------------------------------------------------------------
# HighsImplications::separateImpliedBounds(lpRelaxation, sol, cutpool,
# feastol) (HighsImplications.cpp:541)

.hmip_HighsImplications_separateImpliedBounds <- function(impl, lpRelaxation, sol, cutpool, feastol) {
  mipsolver <- impl$mipsolver
  mipdata <- mipsolver$mipdata_
  globaldomain <- mipdata$domain
  ct <- mipdata$cliquetable
  numboundchgs <- 0L
  fraccols <- .hmip_cuts_frac_cols(.hmip_HighsLpRelaxation_getFractionalIntegers(lpRelaxation))

  if (!.hmip_HighsCliqueTable_isFull(ct)) {
    oldNumQueries <- ct$numNeighbourhoodQueries
    oldNumEntries <- .hmip_HighsCliqueTable_getNumEntries(ct)
    for (col in fraccols) {
      if (globaldomain$col_lower_[col + 1L] != 0 || globaldomain$col_upper_[col + 1L] != 1 ||
          (.hmip_HighsImplications_implicationsCached(impl, col, 0L) &&
           .hmip_HighsImplications_implicationsCached(impl, col, 1L)))
        next
      pr <- .hmip_HighsImplications_runProbing(impl, col, numboundchgs)
      numboundchgs <- pr$numReductions
      if (isTRUE(pr$result) && .hmip_HighsDomain_infeasible(globaldomain)) return(invisible())
      if (.hmip_HighsCliqueTable_isFull(ct)) break
    }
    numNewEntries <- .hmip_HighsCliqueTable_getNumEntries(ct) - oldNumEntries
    impl$nextCleanupCall <- impl$nextCleanupCall - max(0L, numNewEntries)
    if (impl$nextCleanupCall < 0) {
      .hmip_HighsCliqueTable_runCliqueMerging(ct, globaldomain)
      impl$nextCleanupCall <- min(mipdata$numCliqueEntriesAfterFirstPresolve, .hmip_HighsCliqueTable_getNumEntries(ct))
    }
    ct$numNeighbourhoodQueries <- oldNumQueries
  }

  addImplicCuts <- function(col, val) {
    gi <- .hmip_HighsImplications_getImplications(impl, col, val)
    if (.hmip_HighsDomain_infeasible(globaldomain)) return(FALSE)
    if (isTRUE(gi$infeasible)) {
      if (val == 1L) .hmip_HighsCutPool_addCut(cutpool, mipsolver, col, 1, 1L, 0, FALSE, TRUE, FALSE)
      else .hmip_HighsCutPool_addCut(cutpool, mipsolver, col, -1, 1L, -1, FALSE, TRUE, FALSE)
      return(NA)
    }
    im <- gi$implics
    for (q in seq_along(im$column)) {
      icol <- im$column[q]; bv <- im$boundval[q]
      lb <- globaldomain$col_lower_[icol + 1L]; ub <- globaldomain$col_upper_[icol + 1L]
      if (im$boundtype[q] == .hmip_kUpper) {
        if (bv + feastol >= ub) next
        vals <- if (val == 1L) c(1, ub - bv) else c(1, bv - ub)
        rhs <- if (val == 1L) ub else bv
      } else {
        if (bv - feastol <= lb) next
        vals <- if (val == 1L) c(-1, bv - lb) else c(-1, lb - bv)
        rhs <- if (val == 1L) -lb else -bv
      }
      inds <- c(icol, col)
      viol <- sol[inds[1L] + 1L] * vals[1L] + sol[inds[2L] + 1L] * vals[2L] - rhs
      if (viol > feastol)
        .hmip_HighsCutPool_addCut(cutpool, mipsolver, inds, vals, 2L, rhs,
                                  .hmip_variableType(mipsolver, icol) != 0L, FALSE, FALSE, FALSE)
    }
    TRUE
  }

  for (col in fraccols) {
    if (globaldomain$col_lower_[col + 1L] != 0 || globaldomain$col_upper_[col + 1L] != 1) next
    if (.hmip_HighsImplications_implicationsCached(impl, col, 1L)) {
      r <- addImplicCuts(col, 1L)
      if (identical(r, FALSE)) return(invisible())
      if (is.na(r)) next
    }
    if (.hmip_HighsImplications_implicationsCached(impl, col, 0L)) {
      r <- addImplicCuts(col, 0L)
      if (identical(r, FALSE)) return(invisible())
      if (is.na(r)) next
    }
  }
  invisible()
}
