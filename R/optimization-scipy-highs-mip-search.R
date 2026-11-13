# The branch and bound search of the HiGHS 1.12.0 MIP solver as SciPy
# 1.18.1 vendors it: mip/HighsPseudocost, mip/HighsNodeQueue,
# mip/HighsRedcostFixing, mip/HighsLpRelaxation and mip/HighsSearch, with
# the members of class Highs that the LP relaxation uses (over the ported
# HEkk simplex solver) and the HighsHashTable of util/HighsHash.h that
# decides the order of the fractional integers. HiGHS is MIT licensed; the
# notice is in R/optimization-scipy-highs-util.R.
#
# Conventions (API-B2.md): column, row, node and stack indices stored in
# vectors are 0-based as in C++ and read as v[i + 1L]; a HighsDomainChange
# is list(boundval, column, boundtype) with boundtype 0 (kLower) or 1
# (kUpper), a vector of them is the same list with vector fields; a
# HighsCDouble is c(hi, lo) (the .cd_* helpers).

.hmip_b2_lower <- 0L
.hmip_b2_upper <- 1L
.hmip_b2_unspecified <- function() list(type = -2L, index = 0L)
.hmip_b2_iinf <- 2147483647

# ---------------------------------------------------------------------------
# HighsPseudocost (mip/HighsPseudocost.h, HighsPseudocost.cpp)

.hmip_b2_minThreshold <- 1e-6

.hmip_new_HighsPseudocost <- function(mipsolver = NULL) {
  pc <- new.env(parent = emptyenv())
  pc$pseudocostup <- pc$pseudocostdown <- numeric()
  pc$nsamplesup <- pc$nsamplesdown <- integer()
  pc$inferencesup <- pc$inferencesdown <- numeric()
  pc$ninferencesup <- pc$ninferencesdown <- integer()
  pc$ncutoffsup <- pc$ncutoffsdown <- integer()
  pc$conflictscoreup <- pc$conflictscoredown <- numeric()
  pc$conflict_weight <- 0; pc$conflict_avg_score <- 0
  pc$cost_total <- 0; pc$inferences_total <- 0
  pc$nsamplestotal <- 0; pc$ninferencestotal <- 0; pc$ncutoffstotal <- 0
  pc$minreliable <- 0L; pc$degeneracyFactor <- 1
  if (!is.null(mipsolver)) .hmip_HighsPseudocost_init(pc, mipsolver)
  pc
}

# HighsPseudocost::HighsPseudocost(const HighsMipSolver&), assigned in place
.hmip_HighsPseudocost_init <- function(pc, mipsolver) {
  n <- mipsolver$model_$num_col_
  pc$pseudocostup <- pc$pseudocostdown <- numeric(n)
  pc$nsamplesup <- pc$nsamplesdown <- integer(n)
  pc$inferencesup <- pc$inferencesdown <- numeric(n)
  pc$ninferencesup <- pc$ninferencesdown <- integer(n)
  pc$ncutoffsup <- pc$ncutoffsdown <- integer(n)
  pc$conflictscoreup <- pc$conflictscoredown <- numeric(n)
  pc$conflict_weight <- 1; pc$conflict_avg_score <- 0
  pc$cost_total <- 0; pc$inferences_total <- 0
  pc$nsamplestotal <- 0; pc$ninferencestotal <- 0; pc$ncutoffstotal <- 0
  pc$minreliable <- as.integer(mipsolver$options_mip_$mip_pscost_minreliable)
  pc$degeneracyFactor <- 1
  init <- mipsolver$pscostinit
  if (!is.null(init)) {
    pc$cost_total <- init$cost_total
    pc$inferences_total <- init$inferences_total
    pc$nsamplestotal <- init$nsamplestotal
    pc$ninferencestotal <- init$ninferencestotal
    pc$conflict_avg_score <- init$conflict_avg_score * n
    for (i in seq_len(n)) {
      o <- .hmip_HighsPostsolveStack_getOrigColIndex(mipsolver$mipdata_$postSolveStack, i - 1L) + 1L
      pc$pseudocostup[i] <- init$pseudocostup[o]; pc$nsamplesup[i] <- init$nsamplesup[o]
      pc$pseudocostdown[i] <- init$pseudocostdown[o]; pc$nsamplesdown[i] <- init$nsamplesdown[o]
      pc$inferencesup[i] <- init$inferencesup[o]; pc$ninferencesup[i] <- init$ninferencesup[o]
      pc$inferencesdown[i] <- init$inferencesdown[o]; pc$ninferencesdown[i] <- init$ninferencesdown[o]
      pc$conflictscoreup[i] <- init$conflictscoreup[o]; pc$conflictscoredown[i] <- init$conflictscoredown[o]
    }
  }
  invisible(pc)
}

# HighsPseudocostInitialization (both constructors)
.hmip_new_HighsPseudocostInitialization <- function(pscost, maxCount, postsolveStack = NULL) {
  ncol <- length(pscost$pseudocostup)
  out <- list(cost_total = pscost$cost_total, inferences_total = pscost$inferences_total,
              conflict_avg_score = pscost$conflict_avg_score / (ncol * pscost$conflict_weight),
              nsamplestotal = min(1, pscost$nsamplestotal), ninferencestotal = min(1, pscost$ninferencestotal))
  if (is.null(postsolveStack)) {
    out$pseudocostup <- pscost$pseudocostup; out$pseudocostdown <- pscost$pseudocostdown
    out$nsamplesup <- pmin(pscost$nsamplesup, maxCount); out$nsamplesdown <- pmin(pscost$nsamplesdown, maxCount)
    out$inferencesup <- pscost$inferencesup; out$inferencesdown <- pscost$inferencesdown
    out$ninferencesup <- pmin(pscost$ninferencesup, 1L); out$ninferencesdown <- pmin(pscost$ninferencesdown, 1L)
    out$conflictscoreup <- pscost$conflictscoreup / pscost$conflict_weight
    out$conflictscoredown <- pscost$conflictscoredown / pscost$conflict_weight
    return(out)
  }
  no <- .hmip_HighsPostsolveStack_getOrigNumCol(postsolveStack)
  out$pseudocostup <- out$pseudocostdown <- numeric(no)
  out$nsamplesup <- out$nsamplesdown <- integer(no)
  out$inferencesup <- out$inferencesdown <- numeric(no)
  out$ninferencesup <- out$ninferencesdown <- integer(no)
  out$conflictscoreup <- out$conflictscoredown <- numeric(no)
  for (i in seq_len(ncol)) {
    o <- .hmip_HighsPostsolveStack_getOrigColIndex(postsolveStack, i - 1L) + 1L
    out$pseudocostup[o] <- pscost$pseudocostup[i]; out$pseudocostdown[o] <- pscost$pseudocostdown[i]
    out$nsamplesup[o] <- min(maxCount, pscost$nsamplesup[i])
    out$nsamplesdown[o] <- min(maxCount, pscost$nsamplesdown[i])
    out$inferencesup[o] <- pscost$inferencesup[i]; out$inferencesdown[o] <- pscost$inferencesdown[i]
    out$ninferencesup[o] <- 1L; out$ninferencesdown[o] <- 1L
    out$conflictscoreup[o] <- pscost$conflictscoreup[i] / pscost$conflict_weight
    out$conflictscoredown[o] <- pscost$conflictscoredown[i] / pscost$conflict_weight
  }
  out
}

.hmip_HighsPseudocost_subtractBase <- function(pc, base) {
  pc$pseudocostup <- pc$pseudocostup - base$pseudocostup
  pc$pseudocostdown <- pc$pseudocostdown - base$pseudocostdown
  pc$nsamplesup <- pc$nsamplesup - base$nsamplesup
  pc$nsamplesdown <- pc$nsamplesdown - base$nsamplesdown
  invisible(pc)
}
.hmip_HighsPseudocost_increaseConflictWeight <- function(pc) {
  pc$conflict_weight <- pc$conflict_weight * 1.02
  if (pc$conflict_weight > 1000) {
    scale <- 1 / pc$conflict_weight
    pc$conflict_weight <- 1
    pc$conflict_avg_score <- pc$conflict_avg_score * scale
    pc$conflictscoreup <- pc$conflictscoreup * scale
    pc$conflictscoredown <- pc$conflictscoredown * scale
  }
  invisible(pc)
}
.hmip_HighsPseudocost_setDegeneracyFactor <- function(pc, f) pc$degeneracyFactor <- f
.hmip_HighsPseudocost_increaseConflictScoreUp <- function(pc, col) {
  pc$conflictscoreup[col + 1L] <- pc$conflictscoreup[col + 1L] + pc$conflict_weight
  pc$conflict_avg_score <- pc$conflict_avg_score + pc$conflict_weight
}
.hmip_HighsPseudocost_increaseConflictScoreDown <- function(pc, col) {
  pc$conflictscoredown[col + 1L] <- pc$conflictscoredown[col + 1L] + pc$conflict_weight
  pc$conflict_avg_score <- pc$conflict_avg_score + pc$conflict_weight
}
.hmip_HighsPseudocost_setMinReliable <- function(pc, m) pc$minreliable <- as.integer(m)
.hmip_HighsPseudocost_getMinReliable <- function(pc) pc$minreliable
.hmip_HighsPseudocost_getNumObservations <- function(pc, col) pc$nsamplesup[col + 1L] + pc$nsamplesdown[col + 1L]
.hmip_HighsPseudocost_getNumObservationsUp <- function(pc, col) pc$nsamplesup[col + 1L]
.hmip_HighsPseudocost_getNumObservationsDown <- function(pc, col) pc$nsamplesdown[col + 1L]
.hmip_HighsPseudocost_addCutoffObservation <- function(pc, col, upbranch) {
  pc$ncutoffstotal <- pc$ncutoffstotal + 1
  if (upbranch) pc$ncutoffsup[col + 1L] <- pc$ncutoffsup[col + 1L] + 1L
  else pc$ncutoffsdown[col + 1L] <- pc$ncutoffsdown[col + 1L] + 1L
}
.hmip_HighsPseudocost_addObservation <- function(pc, col, delta, objdelta) {
  j <- col + 1L
  if (delta > 0) {
    unit_gain <- objdelta / delta
    d <- unit_gain - pc$pseudocostup[j]
    pc$nsamplesup[j] <- pc$nsamplesup[j] + 1L
    pc$pseudocostup[j] <- pc$pseudocostup[j] + d / pc$nsamplesup[j]
  } else {
    unit_gain <- -objdelta / delta
    d <- unit_gain - pc$pseudocostdown[j]
    pc$nsamplesdown[j] <- pc$nsamplesdown[j] + 1L
    pc$pseudocostdown[j] <- pc$pseudocostdown[j] + d / pc$nsamplesdown[j]
  }
  d <- unit_gain - pc$cost_total
  pc$nsamplestotal <- pc$nsamplestotal + 1
  pc$cost_total <- pc$cost_total + d / pc$nsamplestotal
  invisible(pc)
}
.hmip_HighsPseudocost_addInferenceObservation <- function(pc, col, ninferences, upbranch) {
  j <- col + 1L
  d <- ninferences - pc$inferences_total
  pc$ninferencestotal <- pc$ninferencestotal + 1
  pc$inferences_total <- pc$inferences_total + d / pc$ninferencestotal
  if (upbranch) {
    d <- ninferences - pc$inferencesup[j]
    pc$ninferencesup[j] <- pc$ninferencesup[j] + 1L
    pc$inferencesup[j] <- pc$inferencesup[j] + d / pc$ninferencesup[j]
  } else {
    d <- ninferences - pc$inferencesdown[j]
    pc$ninferencesdown[j] <- pc$ninferencesdown[j] + 1L
    pc$inferencesdown[j] <- pc$inferencesdown[j] + d / pc$ninferencesdown[j]
  }
  invisible(pc)
}
.hmip_HighsPseudocost_isReliable <- function(pc, col)
  min(pc$nsamplesup[col + 1L], pc$nsamplesdown[col + 1L]) >= pc$minreliable
.hmip_HighsPseudocost_isReliableUp <- function(pc, col) pc$nsamplesup[col + 1L] >= pc$minreliable
.hmip_HighsPseudocost_isReliableDown <- function(pc, col) pc$nsamplesdown[col + 1L] >= pc$minreliable
.hmip_HighsPseudocost_getAvgPseudocost <- function(pc) pc$cost_total

.hmip_b2_pscost <- function(pc, ns, cost, mult, offset) {
  if (is.null(offset)) {
    if (ns == 0L) return(mult * pc$cost_total)
    return(mult * cost)
  }
  if (ns == 0L || ns < pc$minreliable) {
    weightPs <- if (ns == 0L) 0 else 0.9 + 0.1 * ns / pc$minreliable
    cst <- weightPs * cost
    cst <- cst + (1 - weightPs) * pc$cost_total
  } else cst <- cost
  mult * (offset + cst)
}
.hmip_HighsPseudocost_getPseudocostUp <- function(pc, col, frac, offset = NULL)
  .hmip_b2_pscost(pc, pc$nsamplesup[col + 1L], pc$pseudocostup[col + 1L], ceiling(frac) - frac, offset)
.hmip_HighsPseudocost_getPseudocostDown <- function(pc, col, frac, offset = NULL)
  .hmip_b2_pscost(pc, pc$nsamplesdown[col + 1L], pc$pseudocostdown[col + 1L], frac - floor(frac), offset)
.hmip_HighsPseudocost_getConflictScoreUp <- function(pc, col) pc$conflictscoreup[col + 1L] / pc$conflict_weight
.hmip_HighsPseudocost_getConflictScoreDown <- function(pc, col) pc$conflictscoredown[col + 1L] / pc$conflict_weight

.hmip_b2_mapScore <- function(score) 1 - 1 / (1 + score)

.hmip_HighsPseudocost_getScore <- function(pc, col, upcost, downcost) {
  j <- col + 1L; mt <- .hmip_b2_minThreshold
  costScore <- max(upcost, mt) * max(downcost, mt) / max(mt, pc$cost_total * pc$cost_total)
  inferenceScore <- max(pc$inferencesup[j], mt) * max(pc$inferencesdown[j], mt) /
    max(mt, pc$inferences_total * pc$inferences_total)
  cutOffScoreUp <- pc$ncutoffsup[j] / max(1, pc$ncutoffsup[j] + pc$nsamplesup[j])
  cutOffScoreDown <- pc$ncutoffsdown[j] / max(1, pc$ncutoffsdown[j] + pc$nsamplesdown[j])
  avgCutoffs <- pc$ncutoffstotal / max(1, pc$ncutoffstotal + pc$nsamplestotal)
  cutoffScore <- max(cutOffScoreUp, mt) * max(cutOffScoreDown, mt) / max(mt, avgCutoffs * avgCutoffs)
  conflictScoreUp <- pc$conflictscoreup[j] / pc$conflict_weight
  conflictScoreDown <- pc$conflictscoredown[j] / pc$conflict_weight
  conflictScoreAvg <- pc$conflict_avg_score / (pc$conflict_weight * length(pc$conflictscoreup))
  conflictScore <- max(conflictScoreUp, mt) * max(conflictScoreDown, mt) /
    max(mt, conflictScoreAvg * conflictScoreAvg)
  m <- .hmip_b2_mapScore
  m(costScore) / pc$degeneracyFactor +
    pc$degeneracyFactor * (1e-2 * m(conflictScore) + 1e-4 * (m(cutoffScore) + m(inferenceScore)))
}
# getScore(col, frac)
.hmip_HighsPseudocost_getScore2 <- function(pc, col, frac)
  .hmip_HighsPseudocost_getScore(pc, col, .hmip_HighsPseudocost_getPseudocostUp(pc, col, frac),
                                 .hmip_HighsPseudocost_getPseudocostDown(pc, col, frac))

.hmip_b2_score_dir <- function(pc, col, pcost, inf, ncut, nsamp, cscore, size) {
  mt <- .hmip_b2_minThreshold
  costScore <- pcost / max(mt, pc$cost_total)
  inferenceScore <- inf / max(mt, pc$inferences_total)
  cutOffScore <- ncut / max(1, ncut + nsamp)
  avgCutoffs <- pc$ncutoffstotal / max(1, pc$ncutoffstotal + pc$nsamplestotal)
  cutoffScore <- cutOffScore / max(mt, avgCutoffs)
  conflictScoreDir <- cscore / pc$conflict_weight
  conflictScoreAvg <- pc$conflict_avg_score / (pc$conflict_weight * size)
  conflictScore <- conflictScoreDir / max(mt, conflictScoreAvg)
  m <- .hmip_b2_mapScore
  m(costScore) + (1e-2 * m(conflictScore) + 1e-4 * (m(cutoffScore) + m(inferenceScore)))
}
.hmip_HighsPseudocost_getScoreUp <- function(pc, col, frac) {
  j <- col + 1L
  .hmip_b2_score_dir(pc, col, .hmip_HighsPseudocost_getPseudocostUp(pc, col, frac), pc$inferencesup[j],
                     pc$ncutoffsup[j], pc$nsamplesup[j], pc$conflictscoreup[j], length(pc$conflictscoreup))
}
.hmip_HighsPseudocost_getScoreDown <- function(pc, col, frac) {
  j <- col + 1L
  .hmip_b2_score_dir(pc, col, .hmip_HighsPseudocost_getPseudocostDown(pc, col, frac), pc$inferencesdown[j],
                     pc$ncutoffsdown[j], pc$nsamplesdown[j], pc$conflictscoredown[j], length(pc$conflictscoredown))
}
.hmip_HighsPseudocost_getAvgInferencesUp <- function(pc, col) pc$inferencesup[col + 1L]
.hmip_HighsPseudocost_getAvgInferencesDown <- function(pc, col) pc$inferencesdown[col + 1L]

# ---------------------------------------------------------------------------
# HighsNodeQueue (mip/HighsNodeQueue.h, HighsNodeQueue.cpp). The three red
# black trees order unique keys that end in the node id, so their minimum,
# last element and predecessor are found from the keys directly. The per
# column std::set<pair<double, int64>> are kept as unsorted pairs and
# ordered by (value, node) when traversed.

.hmip_new_HighsNodeQueue <- function() {
  nq <- new.env(parent = emptyenv())
  nq$nodes <- list()
  nq$freeslots <- numeric()
  nq$colLowerNodes <- list(); nq$colUpperNodes <- list()
  nq$lowerSet <- numeric(); nq$estimSet <- numeric(); nq$suboptSet <- numeric()
  nq$numSuboptimal <- 0
  nq$optimality_limit <- Inf
  nq$numCol <- 0L
  nq
}

.hmip_b2_nodeset <- function() list(first = numeric(), second = numeric())

.hmip_HighsNodeQueue_setNumCol <- function(nq, numCol) {
  if (nq$numCol == numCol) return(invisible())
  nq$numCol <- as.integer(numCol)
  nq$colLowerNodes <- rep(list(.hmip_b2_nodeset()), numCol)
  nq$colUpperNodes <- rep(list(.hmip_b2_nodeset()), numCol)
  invisible()
}
.hmip_HighsNodeQueue_setOptimalityLimit <- function(nq, x) nq$optimality_limit <- x

.hmip_b2_node <- function(nq, id) nq$nodes[[id + 1L]]

# key orders of NodeLowerRbTree, NodeHybridEstimRbTree and SuboptimalNodeRbTree
.hmip_b2_tree_order <- function(nq, ids, tree) {
  if (!length(ids)) return(ids)
  nd <- nq$nodes[ids + 1L]
  lb <- vapply(nd, `[[`, 0, "lower_bound")
  est <- vapply(nd, `[[`, 0, "estimate")
  sz <- vapply(nd, function(x) length(x$domchgstack$column), 0)
  o <- switch(tree,
    lower = order(lb, sz, est, ids),
    estim = order(0.5 * lb + 0.5 * est, -sz, ids),
    subopt = order(lb, ids))
  ids[o]
}
.hmip_b2_tree_min <- function(nq, tree) {
  ids <- nq[[paste0(tree, "Set")]]
  if (!length(ids)) return(-1)
  .hmip_b2_tree_order(nq, ids, tree)[1L]
}

.hmip_b2_set_order <- function(s) order(s$first, s$second)

.hmip_HighsNodeQueue_link_domchgs <- function(nq, node) {
  nd <- nq$nodes[[node + 1L]]
  st <- nd$domchgstack
  for (i in seq_along(st$column)) {
    col <- st$column[i] + 1L
    if (st$boundtype[i] == .hmip_b2_lower) {
      s <- nq$colLowerNodes[[col]]
      nq$colLowerNodes[[col]] <- list(first = c(s$first, st$boundval[i]), second = c(s$second, node))
    } else {
      s <- nq$colUpperNodes[[col]]
      nq$colUpperNodes[[col]] <- list(first = c(s$first, st$boundval[i]), second = c(s$second, node))
    }
  }
}
.hmip_HighsNodeQueue_unlink_domchgs <- function(nq, node) {
  nd <- nq$nodes[[node + 1L]]
  st <- nd$domchgstack
  for (i in seq_along(st$column)) {
    col <- st$column[i] + 1L
    key <- if (st$boundtype[i] == .hmip_b2_lower) "colLowerNodes" else "colUpperNodes"
    s <- nq[[key]][[col]]
    k <- which(s$first == st$boundval[i] & s$second == node)[1L]
    if (!is.na(k)) {
      lst <- nq[[key]]
      lst[[col]] <- list(first = s$first[-k], second = s$second[-k])
      nq[[key]] <- lst
    }
  }
}

.hmip_HighsNodeQueue_link <- function(nq, node) {
  nd <- nq$nodes[[node + 1L]]
  if (nd$lower_bound > nq$optimality_limit) {
    nq$nodes[[node + 1L]]$estimate <- Inf
    nq$suboptSet <- c(nq$suboptSet, node); nq$numSuboptimal <- nq$numSuboptimal + 1
    .hmip_HighsNodeQueue_link_domchgs(nq, node)
    return(2^(1 - nd$depth))
  }
  nq$estimSet <- c(nq$estimSet, node)
  nq$lowerSet <- c(nq$lowerSet, node)
  .hmip_HighsNodeQueue_link_domchgs(nq, node)
  0
}
.hmip_HighsNodeQueue_unlink <- function(nq, node) {
  if (nq$nodes[[node + 1L]]$estimate == Inf) {
    nq$suboptSet <- nq$suboptSet[nq$suboptSet != node]; nq$numSuboptimal <- nq$numSuboptimal - 1
  } else {
    nq$estimSet <- nq$estimSet[nq$estimSet != node]
    nq$lowerSet <- nq$lowerSet[nq$lowerSet != node]
  }
  .hmip_HighsNodeQueue_unlink_domchgs(nq, node)
  nq$freeslots <- c(nq$freeslots, node)
}

.hmip_HighsNodeQueue_checkGlobalBounds <- function(nq, col, lb, ub, feastol, treeweight) {
  lo <- nq$colLowerNodes[[col + 1L]]; up <- nq$colUpperNodes[[col + 1L]]
  # lower_bound(make_pair(ub + feastol, -1)) .. end and begin .. upper_bound(make_pair(lb - feastol, kHighsIInf))
  delnodes <- sort(unique(c(lo$second[lo$first >= ub + feastol], up$second[up$first <= lb - feastol])))
  for (d in delnodes) {
    if (nq$nodes[[d + 1L]]$estimate != Inf) treeweight <- .cd_add(treeweight, 2^(1 - nq$nodes[[d + 1L]]$depth))
    .hmip_HighsNodeQueue_unlink(nq, d)
  }
  treeweight
}

.hmip_HighsNodeQueue_pruneInfeasibleNodes <- function(nq, globaldomain, feastol) {
  treeweight <- c(0, 0)
  repeat {
    if (.hmip_HighsDomain_infeasible(globaldomain)) break
    numchgs <- length(.hmip_HighsDomain_getDomainChangeStack(globaldomain)$column)
    for (i in seq_len(nq$numCol) - 1L)
      treeweight <- .hmip_HighsNodeQueue_checkGlobalBounds(nq, i, globaldomain$col_lower_[i + 1L],
                                                           globaldomain$col_upper_[i + 1L], feastol, treeweight)
    numopennodes <- .hmip_HighsNodeQueue_numNodes(nq)
    if (numopennodes == 0) break
    for (i in seq_len(nq$numCol) - 1L) {
      lo <- nq$colLowerNodes[[i + 1L]]
      if (length(lo$first) == numopennodes) {
        globallb <- lo$first[.hmip_b2_set_order(lo)[1L]]
        if (globallb > globaldomain$col_lower_[i + 1L]) {
          .hmip_HighsDomain_changeBound3(globaldomain, .hmip_b2_lower, i, globallb, .hmip_b2_unspecified())
          if (.hmip_HighsDomain_infeasible(globaldomain)) break
        }
      }
      up <- nq$colUpperNodes[[i + 1L]]
      if (length(up$first) == numopennodes) {
        o <- .hmip_b2_set_order(up)
        globalub <- up$first[o[length(o)]]
        if (globalub < globaldomain$col_upper_[i + 1L]) {
          .hmip_HighsDomain_changeBound3(globaldomain, .hmip_b2_upper, i, globalub, .hmip_b2_unspecified())
          if (.hmip_HighsDomain_infeasible(globaldomain)) break
        }
      }
    }
    .hmip_HighsDomain_propagate(globaldomain)
    if (numchgs == length(.hmip_HighsDomain_getDomainChangeStack(globaldomain)$column)) break
  }
  .cd_val(treeweight)
}

.hmip_HighsNodeQueue_pruneNode <- function(nq, nodeId) {
  nd <- nq$nodes[[nodeId + 1L]]
  tw <- if (nd$estimate != Inf) 2^(1 - nd$depth) else 0
  .hmip_HighsNodeQueue_unlink(nq, nodeId)
  tw
}

.hmip_HighsNodeQueue_performBounding <- function(nq, upper_limit) {
  if (!length(nq$lowerSet)) return(0)
  treeweight <- c(0, 0)
  ord <- .hmip_b2_tree_order(nq, nq$lowerSet, "lower")
  k <- length(ord)
  while (k >= 1L) {
    node <- ord[k]
    if (nq$nodes[[node + 1L]]$lower_bound < upper_limit) break
    treeweight <- .cd_add(treeweight, .hmip_HighsNodeQueue_pruneNode(nq, node))
    k <- k - 1L
  }
  if (nq$optimality_limit < upper_limit) {
    while (k >= 1L) {
      node <- ord[k]
      if (nq$nodes[[node + 1L]]$lower_bound < nq$optimality_limit) break
      nq$estimSet <- nq$estimSet[nq$estimSet != node]
      nq$lowerSet <- nq$lowerSet[nq$lowerSet != node]
      treeweight <- .cd_add(treeweight, 2^(1 - nq$nodes[[node + 1L]]$depth))
      nq$nodes[[node + 1L]]$estimate <- Inf
      nq$suboptSet <- c(nq$suboptSet, node); nq$numSuboptimal <- nq$numSuboptimal + 1
      k <- k - 1L
    }
  }
  if (nq$numSuboptimal) {
    ord <- .hmip_b2_tree_order(nq, nq$suboptSet, "subopt")
    for (node in rev(ord)) {
      if (nq$nodes[[node + 1L]]$lower_bound < upper_limit) break
      .hmip_HighsNodeQueue_unlink(nq, node)
    }
  }
  .cd_val(treeweight)
}

.hmip_HighsNodeQueue_emplaceNode <- function(nq, domchgs, branchings, lower_bound, estimate, depth) {
  node <- list(domchgstack = domchgs, branchings = as.integer(branchings), lower_bound = lower_bound,
               estimate = estimate, depth = depth)
  if (!length(nq$freeslots)) {
    pos <- length(nq$nodes)
  } else {
    pos <- min(nq$freeslots)
    nq$freeslots <- nq$freeslots[-which(nq$freeslots == pos)[1L]]
  }
  nq$nodes[[pos + 1L]] <- node
  .hmip_HighsNodeQueue_link(nq, pos)
}

.hmip_HighsNodeQueue_popBestNode <- function(nq) {
  best <- .hmip_b2_tree_min(nq, "estim")
  .hmip_HighsNodeQueue_unlink(nq, best)
  nq$nodes[[best + 1L]]
}
.hmip_HighsNodeQueue_popBestBoundNode <- function(nq) {
  best <- .hmip_b2_tree_min(nq, "lower")
  .hmip_HighsNodeQueue_unlink(nq, best)
  nq$nodes[[best + 1L]]
}

.hmip_HighsNodeQueue_numNodesUp <- function(nq, col, val = NULL) {
  s <- nq$colLowerNodes[[col + 1L]]
  if (is.null(val)) return(length(s$first))
  sum(s$first > val)
}
.hmip_HighsNodeQueue_numNodesDown <- function(nq, col, val = NULL) {
  s <- nq$colUpperNodes[[col + 1L]]
  if (is.null(val)) return(length(s$first))
  sum(s$first < val)
}
.hmip_HighsNodeQueue_getUpNodes <- function(nq, col) {
  s <- nq$colLowerNodes[[col + 1L]]; o <- .hmip_b2_set_order(s)
  list(first = s$first[o], second = s$second[o])
}
.hmip_HighsNodeQueue_getDownNodes <- function(nq, col) {
  s <- nq$colUpperNodes[[col + 1L]]; o <- .hmip_b2_set_order(s)
  list(first = s$first[o], second = s$second[o])
}

.hmip_HighsNodeQueue_getBestLowerBound <- function(nq) {
  lm <- .hmip_b2_tree_min(nq, "lower")
  lb <- if (lm == -1) Inf else nq$nodes[[lm + 1L]]$lower_bound
  sm <- .hmip_b2_tree_min(nq, "subopt")
  if (sm == -1) return(lb)
  min(nq$nodes[[sm + 1L]]$lower_bound, lb)
}
.hmip_HighsNodeQueue_getBestBoundDomchgStackSize <- function(nq) {
  lm <- .hmip_b2_tree_min(nq, "lower")
  sz <- if (lm == -1) .hmip_b2_iinf else length(nq$nodes[[lm + 1L]]$domchgstack$column)
  sm <- .hmip_b2_tree_min(nq, "subopt")
  if (sm == -1) return(sz)
  min(length(nq$nodes[[sm + 1L]]$domchgstack$column), sz)
}

# HighsNodeQueue::clear keeps numCol and resets everything else, including
# optimality_limit
.hmip_HighsNodeQueue_clear <- function(nq) {
  n <- nq$numCol
  fresh <- .hmip_new_HighsNodeQueue()
  for (k in ls(fresh, all.names = TRUE)) assign(k, get(k, fresh), envir = nq)
  .hmip_HighsNodeQueue_setNumCol(nq, n)
  invisible(nq)
}
.hmip_HighsNodeQueue_numNodes <- function(nq) length(nq$nodes) - length(nq$freeslots)
.hmip_HighsNodeQueue_numActiveNodes <- function(nq) length(nq$nodes) - length(nq$freeslots) - nq$numSuboptimal
.hmip_HighsNodeQueue_empty <- function(nq) .hmip_HighsNodeQueue_numActiveNodes(nq) == 0

# ---------------------------------------------------------------------------
# HighsRedcostFixing (mip/HighsRedcostFixing.cpp). A std::multimap is the
# pair of vectors list(first, second) in map order; emplace_hint at
# lower_bound inserts before the equal keys.

.hmip_new_HighsRedcostFixing <- function() {
  rf <- new.env(parent = emptyenv())
  rf$lurkingColUpper <- list(); rf$lurkingColLower <- list()
  rf
}

.hmip_HighsRedcostFixing_getLurkingBounds <- function(rf, mipsolver) {
  out <- list()
  if (!length(rf$lurkingColLower)) return(out)
  md <- mipsolver$mipdata_
  for (col in md$integral_cols) {
    lo <- rf$lurkingColLower[[col + 1L]]
    for (k in seq_along(lo$first)) if (lo$second[k] > md$domain$col_lower_[col + 1L])
      out[[length(out) + 1L]] <- list(first = lo$first[k], second = list(boundval = lo$second[k], column = col,
                                                                          boundtype = .hmip_b2_lower))
    up <- rf$lurkingColUpper[[col + 1L]]
    for (k in seq_along(up$first)) if (up$second[k] < md$domain$col_upper_[col + 1L])
      out[[length(out) + 1L]] <- list(first = up$first[k], second = list(boundval = up$second[k], column = col,
                                                                          boundtype = .hmip_b2_upper))
  }
  out
}

.hmip_HighsRedcostFixing_propagateRootRedcost <- function(rf, mipsolver) {
  if (!length(rf$lurkingColLower)) return(invisible())
  md <- mipsolver$mipdata_
  dom <- md$domain
  for (col in md$integral_cols) {
    j <- col + 1L
    for (side in c("lurkingColLower", "lurkingColUpper")) {
      s <- rf[[side]][[j]]
      keep <- s$first > md$lower_bound
      lst <- rf[[side]]; lst[[j]] <- list(first = s$first[keep], second = s$second[keep]); rf[[side]] <- lst
    }
    lo <- rf$lurkingColLower[[j]]
    for (k in which(lo$first >= md$upper_limit)) {
      if (lo$second[k] > dom$col_lower_[j]) {
        .hmip_HighsDomain_changeBound3(dom, .hmip_b2_lower, col, lo$second[k], .hmip_b2_unspecified())
        if (.hmip_HighsDomain_infeasible(dom)) return(invisible())
      }
    }
    up <- rf$lurkingColUpper[[j]]
    for (k in which(up$first >= md$upper_limit)) {
      if (up$second[k] < dom$col_upper_[j]) {
        .hmip_HighsDomain_changeBound3(dom, .hmip_b2_upper, col, up$second[k], .hmip_b2_unspecified())
        if (.hmip_HighsDomain_infeasible(dom)) return(invisible())
      }
    }
  }
  .hmip_HighsDomain_propagate(dom)
  invisible()
}

.hmip_HighsRedcostFixing_propagateRedCost <- function(mipsolver, localdomain, lp) {
  md <- mipsolver$mipdata_
  lpredcost <- .hmip_HighsLpRelaxation_getSolution(lp)$col_dual
  lpobjective <- .hmip_HighsLpRelaxation_getObjective(lp)
  gapcd <- .cd_add(c(md$upper_limit, 0), -lpobjective)
  gap <- .cd_val(gapcd)
  tolerance <- max(10 * md$feastol, md$epsilon * gap)
  bc <- list(boundval = numeric(), column = integer(), boundtype = integer())
  push <- function(v, c, t) { bc$boundval <<- c(bc$boundval, v); bc$column <<- c(bc$column, c); bc$boundtype <<- c(bc$boundtype, t) }
  for (col in md$integral_cols) {
    j <- col + 1L
    lo <- localdomain$col_lower_[j]; up <- localdomain$col_upper_[j]
    if (up == lo) next
    rc <- lpredcost[j]
    if (abs(rc) <= tolerance) next
    maxIncrease <- rc * (up - lo)
    if (maxIncrease > gap) {
      newub <- .cd_val(.hmip_b2_cd_floor(.cd_add(.cd_add(.cd_div(gapcd, rc), lo), md$feastol)))
      if (newub >= up) next
      if (.hmip_HighsDomain_isBinary(md$domain, col)) push(newub, col, .hmip_b2_upper)
      else {
        .hmip_HighsDomain_changeBound3(localdomain, .hmip_b2_upper, col, newub, .hmip_b2_unspecified())
        if (.hmip_HighsDomain_infeasible(localdomain)) return(invisible())
      }
    } else if (maxIncrease < -gap) {
      newlb <- .cd_val(.hmip_b2_cd_ceil(.cd_add(.cd_add(.cd_div(gapcd, rc), up), -md$feastol)))
      if (newlb <= lo) next
      if (.hmip_HighsDomain_isBinary(md$domain, col)) push(newlb, col, .hmip_b2_lower)
      else {
        .hmip_HighsDomain_changeBound3(localdomain, .hmip_b2_lower, col, newlb, .hmip_b2_unspecified())
        if (.hmip_HighsDomain_infeasible(localdomain)) return(invisible())
      }
    }
  }
  nb <- length(bc$column)
  if (!nb) return(invisible())
  apply_all <- function(b) {
    for (k in seq_along(b$column)) {
      .hmip_HighsDomain_changeBound(localdomain, list(boundval = b$boundval[k], column = b$column[k],
                                                      boundtype = b$boundtype[k]), .hmip_b2_unspecified())
      if (.hmip_HighsDomain_infeasible(localdomain)) break
    }
    if (!.hmip_HighsDomain_infeasible(localdomain)) .hmip_HighsDomain_propagate(localdomain)
  }
  proof <- if (nb <= 100) .hmip_HighsLpRelaxation_computeDualProof(lp, md$domain, md$upper_limit, FALSE) else list(ok = FALSE)
  if (proof$ok) {
    old <- .hmip_HighsConflictPool_getNumConflicts(md$conflictPool)
    for (k in seq_len(nb)) {
      dc <- list(boundval = bc$boundval[k], column = bc$column[k], boundtype = bc$boundtype[k])
      if (.hmip_HighsDomain_isActive(localdomain, dc)) next
      .hmip_HighsDomain_conflictAnalyzeReconvergence(localdomain, dc, proof$inds, proof$vals, length(proof$inds),
                                                     proof$rhs, md$conflictPool)
    }
    if (.hmip_HighsConflictPool_getNumConflicts(md$conflictPool) != old) {
      .hmip_HighsDomain_propagate(localdomain)
      if (.hmip_HighsDomain_infeasible(localdomain)) return(invisible())
      act <- vapply(seq_len(nb), function(k) .hmip_HighsDomain_isActive(localdomain, list(
        boundval = bc$boundval[k], column = bc$column[k], boundtype = bc$boundtype[k])), TRUE)
      bc <- list(boundval = bc$boundval[!act], column = bc$column[!act], boundtype = bc$boundtype[!act])
    }
    if (length(bc$column)) apply_all(bc)
  } else apply_all(bc)
  invisible()
}

# floor and ceil of a HighsCDouble (HighsCDouble.h)
.hmip_b2_cd_floor <- function(x) {
  v <- .cd_val(x)
  if (abs(v) < 1) return(if (v >= 0) c(0, 0) else c(-1, 0))
  f <- floor(v)
  .cd_two_sum(f, floor(.cd_val(.cd_add(x, -f))))
}
.hmip_b2_cd_ceil <- function(x) {
  v <- .cd_val(x)
  if (abs(v) < 1) return(if (v <= 0) c(0, 0) else c(1, 0))
  f <- ceiling(v)
  .cd_two_sum(f, ceiling(.cd_val(.cd_add(x, -f))))
}

.hmip_HighsRedcostFixing_addRootRedcost <- function(rf, mipsolver, lpredcost, lpobjective) {
  md <- mipsolver$mipdata_
  n <- mipsolver$model_$num_col_
  empty <- list(first = numeric(), second = numeric())
  if (length(rf$lurkingColLower) < n) {
    rf$lurkingColLower <- c(rf$lurkingColLower, rep(list(empty), n - length(rf$lurkingColLower)))
    rf$lurkingColUpper <- c(rf$lurkingColUpper, rep(list(empty), n - length(rf$lurkingColUpper)))
  }
  .hmip_HighsLpRelaxation_computeBasicDegenerateDuals(md$lp, md$feastol)
  feastol <- md$feastol
  dom <- md$domain
  numLarge <- 0L
  for (col in md$integral_cols)
    if ((dom$col_upper_[col + 1L] - dom$col_lower_[col + 1L]) >= 512 && abs(lpredcost[col + 1L]) > feastol)
      numLarge <- numLarge + 1L
  maxNumStepsExp <- 10L
  q <- numLarge %/% 10L
  expshift <- if (q == 0L) 0L else as.integer(floor(log2(q))) + 1L
  if (expshift > 5L) {
    expshift <- min(expshift, maxNumStepsExp)
    maxNumStepsExp <- maxNumStepsExp - expshift + 5L
  }
  maxNumSteps <- 2^maxNumStepsExp
  for (col in md$integral_cols) {
    j <- col + 1L
    rc <- lpredcost[j]
    if (rc > feastol) {
      lb <- trunc(dom$col_lower_[j])
      s <- rf$lurkingColUpper[[j]]
      if (rc == Inf) {
        rf$lurkingColLower[[j]] <- empty
        rf$lurkingColUpper[[j]] <- list(first = -Inf, second = lb)
        next
      }
      maxub <- if (dom$col_upper_[j] == Inf) lb + maxNumSteps else trunc(floor(dom$col_upper_[j] - 0.5))
      step <- 1
      if (maxub - lb > maxNumSteps) step <- (maxub - lb + maxNumSteps - 1) %/% 2^maxNumStepsExp
      lurkub <- lb
      while (lurkub <= maxub) {
        fracbound <- (lurkub - lb + 1) - 10 * feastol
        req <- fracbound * rc + lpobjective
        if (req >= md$lower_bound) {
          pos <- sum(s$first < req) + 1L
          later <- if (pos <= length(s$first)) s$second[pos:length(s$first)] else numeric()
          if (!any(later < lurkub + step)) {
            s$first <- append(s$first, req, pos - 1L); s$second <- append(s$second, lurkub, pos - 1L)
            if (pos > 1L) {
              del <- which(s$second[seq_len(pos - 1L)] >= lurkub)
              if (length(del)) { s$first <- s$first[-del]; s$second <- s$second[-del] }
            }
          }
        }
        lurkub <- lurkub + step
      }
      rf$lurkingColUpper[[j]] <- s
    } else if (rc < -feastol) {
      ub <- trunc(dom$col_upper_[j])
      s <- rf$lurkingColLower[[j]]
      if (rc == -Inf) {
        rf$lurkingColUpper[[j]] <- empty
        rf$lurkingColLower[[j]] <- list(first = -Inf, second = ub)
        next
      }
      minlb <- if (dom$col_lower_[j] == -Inf) ub - maxNumSteps else trunc(dom$col_lower_[j] + 1.5)
      step <- 1
      if (ub - minlb > maxNumSteps) step <- (ub - minlb + maxNumSteps - 1) %/% 2^maxNumStepsExp
      lurklb <- minlb
      while (lurklb <= ub) {
        fracbound <- (lurklb - ub - 1) + 10 * feastol
        req <- fracbound * rc + lpobjective - feastol
        if (req >= md$lower_bound) {
          pos <- sum(s$first < req) + 1L
          later <- if (pos <= length(s$first)) s$second[pos:length(s$first)] else numeric()
          if (!any(later > lurklb - step)) {
            s$first <- append(s$first, req, pos - 1L); s$second <- append(s$second, lurklb, pos - 1L)
            if (pos > 1L) {
              del <- which(s$second[seq_len(pos - 1L)] <= lurklb)
              if (length(del)) { s$first <- s$first[-del]; s$second <- s$second[-del] }
            }
          }
        }
        lurklb <- lurklb + step
      }
      rf$lurkingColLower[[j]] <- s
    }
  }
  invisible()
}

# ---------------------------------------------------------------------------
# HighsHashTable<HighsInt, V> (util/HighsHash.h): Robin Hood hashing with
# 7 bits of the hash in the metadata byte; iteration runs over the slots.

.hmip_new_HighsHashTable <- function(minCapacity = NULL) {
  ht <- new.env(parent = emptyenv())
  cap <- if (is.null(minCapacity)) 128 else 2^ceiling(log2(max(128, 8 * minCapacity / 7)))
  .hmip_b2_ht_empty(ht, cap)
  ht
}
.hmip_b2_ht_empty <- function(ht, cap) {
  ht$cap <- cap
  ht$bits <- as.integer(round(log2(cap)))
  ht$meta <- integer(cap)
  ht$keys <- numeric(cap)
  ht$vals <- vector("list", cap)
  ht$numElements <- 0
}
# start slot (hash >> numHashShift) and metadata byte of a key
.hmip_b2_ht_slot <- function(ht, key) {
  h <- .hmip_hash(key)
  start <- .hs_u64_top(h, ht$bits)
  list(start = start, meta = bitwOr(as.integer(start %% 128), 128L))
}
.hmip_b2_ht_dist <- function(ht, pos) (pos - ht$meta[pos + 1L]) %% 128
.hmip_b2_ht_findpos <- function(ht, key) {
  s <- .hmip_b2_ht_slot(ht, key)
  mask <- ht$cap
  pos <- s$start
  maxPos <- (s$start + 127) %% mask
  repeat {
    m <- ht$meta[pos + 1L]
    if (m < 128L) return(list(found = FALSE, pos = pos, start = s$start, maxPos = maxPos, meta = s$meta))
    if (m == s$meta && ht$keys[pos + 1L] == key) return(list(found = TRUE, pos = pos))
    if (((pos - s$start) %% mask) > .hmip_b2_ht_dist(ht, pos))
      return(list(found = FALSE, pos = pos, start = s$start, maxPos = maxPos, meta = s$meta))
    pos <- (pos + 1) %% mask
    if (pos == maxPos) break
  }
  list(found = FALSE, pos = pos, start = s$start, maxPos = maxPos, meta = s$meta)
}
.hmip_b2_ht_grow <- function(ht) {
  ok <- which(ht$meta >= 128L)
  keys <- ht$keys[ok]; vals <- ht$vals[ok]
  .hmip_b2_ht_empty(ht, 2 * ht$cap)
  for (i in seq_along(keys)) .hmip_b2_ht_insert(ht, keys[i], vals[[i]])
}
# insert(entry): returns the slot of the new key, or -1 when present
.hmip_b2_ht_insert <- function(ht, key, value) {
  f <- .hmip_b2_ht_findpos(ht, key)
  if (f$found) return(-1)
  if (ht$numElements == (ht$cap * 7) %/% 8 || f$pos == f$maxPos) {
    .hmip_b2_ht_grow(ht)
    return(.hmip_b2_ht_insert(ht, key, value))
  }
  ht$numElements <- ht$numElements + 1
  pos <- f$pos; start <- f$start; maxPos <- f$maxPos; meta <- f$meta
  mask <- ht$cap
  ekey <- key; eval <- value
  first <- -1
  repeat {
    if (ht$meta[pos + 1L] < 128L) {
      ht$meta[pos + 1L] <- meta; ht$keys[pos + 1L] <- ekey; ht$vals[pos + 1L] <- list(eval)
      return(if (first < 0) pos else first)
    }
    cd <- (pos - start) %% mask
    dco <- .hmip_b2_ht_dist(ht, pos)
    if (cd > dco) {
      if (first < 0) first <- pos
      tk <- ht$keys[pos + 1L]; tv <- ht$vals[[pos + 1L]]; tm <- ht$meta[pos + 1L]
      ht$keys[pos + 1L] <- ekey; ht$vals[pos + 1L] <- list(eval); ht$meta[pos + 1L] <- meta
      ekey <- tk; eval <- tv; meta <- tm
      start <- (pos - dco) %% mask
      maxPos <- (start + 127) %% mask
    }
    pos <- (pos + 1) %% mask
    if (pos == maxPos) break
  }
  .hmip_b2_ht_grow(ht)
  .hmip_b2_ht_insert(ht, ekey, eval)
  .hmip_b2_ht_findpos(ht, key)$pos
}
.hmip_HighsHashTable_find <- function(ht, key) {
  f <- .hmip_b2_ht_findpos(ht, key)
  if (f$found) ht$vals[[f$pos + 1L]] else NULL
}
# operator[]: the slot of key, inserted with value `default` when absent
.hmip_HighsHashTable_ref <- function(ht, key, default) {
  f <- .hmip_b2_ht_findpos(ht, key)
  if (f$found) return(f$pos)
  .hmip_b2_ht_insert(ht, key, default)
}
.hmip_HighsHashTable_get <- function(ht, slot) ht$vals[[slot + 1L]]
.hmip_HighsHashTable_set <- function(ht, slot, value) ht$vals[slot + 1L] <- list(value)
.hmip_HighsHashTable_keys <- function(ht) ht$keys[ht$meta >= 128L]
.hmip_HighsHashTable_values <- function(ht) ht$vals[ht$meta >= 128L]
.hmip_HighsHashTable_size <- function(ht) ht$numElements
.hmip_HighsHashTable_clear <- function(ht) if (ht$numElements) .hmip_b2_ht_empty(ht, 128)

# ---------------------------------------------------------------------------
# HighsSparseVectorSum (util/HighsSparseVectorSum.h): HighsCDouble values
# with the nonzero indices in the order of their first addition.

.hmip_b2_svs <- function(n) {
  s <- new.env(parent = emptyenv())
  s$hi <- numeric(n); s$lo <- numeric(n); s$nz <- integer()
  s
}
.hmip_b2_svs_add <- function(s, i, v) {
  j <- i + 1L
  cur <- c(s$hi[j], s$lo[j])
  if (.cd_val(cur) != 0) cur <- .cd_add(cur, v)
  else { cur <- c(v, 0); s$nz <- c(s$nz, i) }
  if (.cd_val(cur) == 0) cur <- c(0x1.0000000000000p-1022, 0)  # DBL_MIN, 2.2250738585072014e-308
  s$hi[j] <- cur[1L]; s$lo[j] <- cur[2L]
}
.hmip_b2_svs_get <- function(s, i) s$hi[i + 1L] + s$lo[i + 1L]
.hmip_b2_svs_cleanup <- function(s, is_zero) {
  nz <- s$nz; num <- length(nz)
  for (i in rev(seq_along(nz))) {
    pos <- nz[i]
    if (is_zero(pos, .hmip_b2_svs_get(s, pos))) {
      s$hi[pos + 1L] <- 0; s$lo[pos + 1L] <- 0
      tmp <- nz[num]; nz[num] <- nz[i]; nz[i] <- tmp
      num <- num - 1L
    }
  }
  s$nz <- nz[seq_len(num)]
}

# ---------------------------------------------------------------------------
# The members of class Highs that HighsLpRelaxation calls, over the ported
# simplex solver. The LP is held densely (A) with the bounds and costs of
# HighsLp; the HEkk instance (ek) is kept across runs while only column
# bounds change and cleared when rows are added or deleted, as
# Highs::addRows and Highs::deleteRows do with kExtendInvertWhenAddingRows
# false. The scale factors are found once (considerScaling) and extended
# for added rows by HighsSparseMatrix::considerRowScaling.

.hmip_new_Highs <- function() {
  h <- new.env(parent = emptyenv())
  h$options <- list(presolve = "choose", solver = "choose", simplex_strategy = 1L,
                    simplex_iteration_limit = Inf, objective_bound = Inf, time_limit = Inf,
                    primal_feasibility_tolerance = 1e-7, dual_feasibility_tolerance = 1e-7,
                    random_seed = 0L, simplex_initial_condition_check = TRUE, output_flag = TRUE,
                    parallel = "choose", small_matrix_value = 1e-9, allowed_matrix_scale_factor = 20L)
  h$lp <- list(num_col_ = 0L, num_row_ = 0L, col_cost_ = numeric(), col_lower_ = numeric(),
               col_upper_ = numeric(), row_lower_ = numeric(), row_upper_ = numeric(),
               A = matrix(0, 0, 0), offset_ = 0)
  h$scale <- NULL; h$scale_tried <- FALSE
  .hmip_b2_highs_invalidate(h)
  h$basis <- list(valid = FALSE, col_status = character(), row_status = character())
  h$ek <- NULL; h$ek_valid <- FALSE; h$iterate <- NULL
  h
}
.hmip_b2_highs_invalidate <- function(h) {
  h$model_status <- "Notset"
  h$solution <- list(value_valid = FALSE, dual_valid = FALSE, col_value = numeric(), col_dual = numeric(),
                     row_value = numeric(), row_dual = numeric())
  h$info <- list(valid = FALSE, simplex_iteration_count = 0, objective_function_value = 0,
                 basis_validity = 0L, primal_solution_status = 0L, dual_solution_status = 0L,
                 num_primal_infeasibilities = -1L, max_primal_infeasibility = Inf,
                 sum_primal_infeasibilities = Inf, num_dual_infeasibilities = -1L,
                 max_dual_infeasibility = Inf, sum_dual_infeasibilities = Inf)
}
.hmip_b2_highs_clear_ekk <- function(h) { h$ek <- NULL; h$ek_valid <- FALSE; h$iterate <- NULL }

.hmip_Highs_setOptionValue <- function(h, name, value) {
  o <- h$options; o[[name]] <- value; h$options <- o
  invisible("Ok")
}
.hmip_Highs_getOptionValue <- function(h, name) h$options[[name]]
.hmip_Highs_getOptions <- function(h) h$options
.hmip_Highs_passOptions <- function(h, options) { h$options <- options; invisible("Ok") }
.hmip_Highs_getRunTime <- function(h) 0

# passModel from a HighsLp (B1 layout: column-wise a_matrix_ with 0-based
# starts and indices)
.hmip_Highs_passModel <- function(h, lp) {
  n <- lp$num_col_; m <- lp$num_row_
  if (!is.null(lp$A)) A <- lp$A else {
    A <- matrix(0, m, n)
    am <- lp$a_matrix_
    for (j in seq_len(n)) {
      k <- if (am$start_[j + 1L] > am$start_[j]) (am$start_[j] + 1L):am$start_[j + 1L] else integer()
      A[am$index_[k] + 1L, j] <- am$value_[k]
    }
  }
  h$lp_cache <- NULL
  h$lp <- list(num_col_ = n, num_row_ = m, col_cost_ = as.numeric(lp$col_cost_),
               col_lower_ = as.numeric(lp$col_lower_), col_upper_ = as.numeric(lp$col_upper_),
               row_lower_ = as.numeric(lp$row_lower_), row_upper_ = as.numeric(lp$row_upper_),
               A = A, offset_ = lp$offset_ %||% 0)
  h$scale <- NULL; h$scale_tried <- FALSE
  h$basis <- list(valid = FALSE, col_status = character(), row_status = character())
  .hmip_b2_highs_clear_ekk(h)
  .hmip_b2_highs_invalidate(h)
  invisible("Ok")
}
.hmip_Highs_clearSolver <- function(h) {
  .hmip_b2_highs_clear_ekk(h)
  h$basis <- list(valid = FALSE, col_status = character(), row_status = character())
  .hmip_b2_highs_invalidate(h)
  invisible("Ok")
}

# getLp: the HighsLp with the column-wise matrix (entries of each column in
# row order, which is the order addRows and deleteRows keep)
.hmip_Highs_getLp <- function(h) {
  if (!is.null(h$lp_cache)) return(h$lp_cache)
  lp <- h$lp
  A <- lp$A
  nz <- if (length(A)) which(A != 0, arr.ind = TRUE) else matrix(0L, 0, 2)
  nz <- nz[order(nz[, 2L], nz[, 1L]), , drop = FALSE]
  cnt <- tabulate(nz[, 2L], lp$num_col_)
  lp$a_matrix_ <- list(format_ = "colwise", start_ = c(0L, cumsum(cnt)), index_ = nz[, 1L] - 1L,
                       value_ = A[nz])
  lp$A <- NULL
  h$lp_cache <- lp
  lp
}
.hmip_Highs_getNumRow <- function(h) h$lp$num_row_
.hmip_Highs_getNumCol <- function(h) h$lp$num_col_
.hmip_Highs_getNumNz <- function(h) sum(h$lp$A != 0)
.hmip_Highs_getModelStatus <- function(h) h$model_status
.hmip_Highs_getInfo <- function(h) h$info
.hmip_Highs_getSolution <- function(h) h$solution
.hmip_Highs_getBasis <- function(h) h$basis
.hmip_Highs_hasInvert <- function(h) !is.null(h$ek) && h$ek_valid

# accommodateAlienBasis (lp_data/HighsSolution.cpp): factor the basic
# columns of the scaled LP and keep at most num_row of them basic in their
# order, completing the basis with the logicals of the rows without a pivot.
# A singular basis that results is repaired at the next invert.
.hmip_accommodate_alien_basis <- function(h, basis) {
  lp <- h$lp; m <- lp$num_row_; n <- lp$num_col_
  cs <- basis$col_status; rs <- basis$row_status
  basic_index <- c(which(cs == "basic"), n + which(rs == "basic"))
  if (!h$scale_tried) { h$scale <- .hs_scale(lp$A, lp$col_cost_); h$scale_tried <- TRUE }
  A <- if (is.null(h$scale)) lp$A else lp$A * outer(h$scale$row, h$scale$col)
  cols <- lapply(basic_index, function(v) {
    if (v > n) return(list(index = v - n, value = 1, logical = TRUE))
    nz <- which(A[, v] != 0)
    list(index = nz, value = A[nz, v])
  })
  f <- .hs_factor_build(m, cols)
  cs[cs == "basic"] <- "nonbasic"; rs[rs == "basic"] <- "nonbasic"
  use <- basic_index[seq_len(min(m, length(basic_index)))]
  cs[use[use <= n]] <- "basic"; rs[use[use > n] - n] <- "basic"
  num_missing <- m - length(use)
  rs[f$row_with_no_pivot[f$rank_deficiency + seq_len(num_missing)]] <- "basic"
  list(col_status = cs, row_status = rs)
}

.hmip_Highs_setBasis <- function(h, basis, origin = "") {
  if (isTRUE(basis$alien)) {
    if (h$lp$num_row_ == 0L) basis$col_status[basis$col_status == "basic"] <- "nonbasic"
    else basis <- .hmip_accommodate_alien_basis(h, basis)
  }
  h$basis <- list(valid = TRUE, col_status = basis$col_status, row_status = basis$row_status)
  .hmip_b2_highs_clear_ekk(h)
  .hmip_b2_highs_invalidate(h)
  invisible("Ok")
}

# changeColsBounds (set form, cols 0-based); the HEkk keeps its invert
.hmip_Highs_changeColsBounds <- function(h, cols, lower, upper) {
  if (!length(cols)) return(invisible("Ok"))
  lp <- h$lp
  lp$col_lower_[cols + 1L] <- lower; lp$col_upper_[cols + 1L] <- upper
  h$lp <- lp
  if (!is.null(h$lp_cache)) {
    h$lp_cache$col_lower_[cols + 1L] <- lower; h$lp_cache$col_upper_[cols + 1L] <- upper
  }
  if (!is.null(h$ek) && h$ek_valid) {
    ek <- h$ek
    cs <- if (is.null(h$scale)) rep(1, lp$num_col_) else h$scale$col
    ek$lp$col_lower[cols + 1L] <- lower / cs[cols + 1L]
    ek$lp$col_upper[cols + 1L] <- upper / cs[cols + 1L]
  }
  .hmip_b2_highs_invalidate(h)
  invisible("Ok")
}
.hmip_Highs_changeColsBoundsRange <- function(h, from, to, lower, upper) {
  k <- seq_len(to - from + 1L)
  .hmip_Highs_changeColsBounds(h, seq.int(from, to), lower[k], upper[k])
}

# changeColsCost (mask form): set the cost of every column whose mask entry is
# nonzero; cost has one entry per column. The HEkk keeps its invert
# (LpAction::kNewCosts) and the model status and solution are invalidated.
.hmip_Highs_changeColsCost <- function(h, mask, cost) {
  cols <- which(mask != 0)
  if (!length(cols)) return(invisible("Ok"))
  lp <- h$lp
  lp$col_cost_[cols] <- cost[cols]
  h$lp <- lp
  if (!is.null(h$lp_cache)) h$lp_cache$col_cost_[cols] <- cost[cols]
  if (!is.null(h$ek) && h$ek_valid) {
    ek <- h$ek
    cs <- if (is.null(h$scale)) rep(1, lp$num_col_) else h$scale$col
    ek$lp$col_cost[cols] <- cost[cols] * cs[cols]
  }
  .hmip_b2_highs_invalidate(h)
  invisible("Ok")
}

# addRows from a row-wise matrix (0-based starts and indices)
.hmip_Highs_addRows <- function(h, num, lower, upper, nnz, starts, index, value) {
  lp <- h$lp
  n <- lp$num_col_
  R <- matrix(0, num, n)
  ends <- c(starts[-1L], nnz)[seq_len(num)]
  for (i in seq_len(num)) {
    k <- if (ends[i] > starts[i]) (starts[i] + 1L):ends[i] else integer()
    R[i, index[k] + 1L] <- value[k]
  }
  lp$A <- rbind(lp$A, R)
  lp$row_lower_ <- c(lp$row_lower_, lower[seq_len(num)])
  lp$row_upper_ <- c(lp$row_upper_, upper[seq_len(num)])
  h$lp_cache <- NULL
  lp$num_row_ <- lp$num_row_ + num
  h$lp <- lp
  if (!is.null(h$scale)) {
    # applyColScale, then considerRowScaling of the new rows
    mx <- 2^h$options$allowed_matrix_scale_factor
    rs <- numeric(num)
    for (i in seq_len(num)) {
      rmax <- max(abs(R[i, ] * h$scale$col), 0)
      if (rmax) {
        v <- 2^floor(log(1 / rmax) / log(2) + 0.5)
        rs[i] <- min(max(1 / mx, v), mx)
      } else rs[i] <- 1
    }
    h$scale$row <- c(h$scale$row, rs)
  }
  if (length(h$basis$col_status))
    h$basis$row_status <- c(h$basis$row_status, rep("basic", num))
  .hmip_b2_highs_invalidate(h)
  .hmip_b2_highs_clear_ekk(h)
  invisible("Ok")
}

# deleteRows by mask (Highs::deleteRows(HighsInt* mask)): returns the mask
# with the new index of each kept row and -1 for deleted rows
.hmip_Highs_deleteRowsMask <- function(h, mask) {
  del <- mask != 0
  newmask <- rep(-1L, length(mask))
  newmask[!del] <- seq_len(sum(!del)) - 1L
  .hmip_b2_highs_delete_rows(h, which(del))
  newmask
}
.hmip_Highs_deleteRowsRange <- function(h, from, to) {
  if (to >= from) .hmip_b2_highs_delete_rows(h, (from:to) + 1L)
  invisible("Ok")
}
.hmip_b2_highs_delete_rows <- function(h, rows) {
  if (!length(rows)) return(invisible())
  lp <- h$lp
  lp$A <- lp$A[-rows, , drop = FALSE]
  lp$row_lower_ <- lp$row_lower_[-rows]; lp$row_upper_ <- lp$row_upper_[-rows]
  h$lp_cache <- NULL
  lp$num_row_ <- lp$num_row_ - length(rows)
  h$lp <- lp
  if (!is.null(h$scale)) h$scale$row <- h$scale$row[-rows]
  if (length(h$basis$row_status)) {
    deleted_nonbasic <- any(h$basis$row_status[rows] != "basic")
    h$basis$row_status <- h$basis$row_status[-rows]
    if (deleted_nonbasic) h$basis$valid <- FALSE
  }
  .hmip_b2_highs_invalidate(h)
  .hmip_b2_highs_clear_ekk(h)
}

# the HighsBasis of the HEkk (HEkk::getHighsBasis)
.hmip_b2_ek_basis <- function(ek) {
  s <- .hs_get_solution(ek)
  list(valid = TRUE, col_status = s$col_status, row_status = s$row_status)
}

# getKktFailures: the primal and dual infeasibilities of the unscaled
# solution with respect to the LP bounds
.hmip_b2_kkt <- function(lp, sol, ptol, dtol) {
  val <- c(sol$col_value, sol$row_value)
  dual <- c(sol$col_dual, sol$row_dual)
  lo <- c(lp$col_lower_, lp$row_lower_); up <- c(lp$col_upper_, lp$row_upper_)
  pinf <- pmax(lo - val, val - up, 0)
  dinf <- numeric(length(val))
  for (k in seq_along(val)) {
    if (lo[k] == up[k]) next
    atl <- lo[k] > -Inf && val[k] <= lo[k] + ptol
    atu <- up[k] < Inf && val[k] >= up[k] - ptol
    dinf[k] <- if (atl && !atu) max(-dual[k], 0) else if (atu && !atl) max(dual[k], 0)
               else if (!atl && !atu) abs(dual[k]) else 0
  }
  list(num_primal_infeasibilities = sum(pinf > ptol), max_primal_infeasibility = max(c(pinf, 0)),
       sum_primal_infeasibilities = sum(pinf[pinf > 0]),
       num_dual_infeasibilities = sum(dinf > dtol), max_dual_infeasibility = max(c(dinf, 0)),
       sum_dual_infeasibilities = sum(dinf[dinf > 0]))
}

# Highs::run for the LP: presolve only without a valid basis, otherwise the
# simplex solver from the basis or the retained HEkk
.hmip_Highs_run <- function(h) {
  lp <- h$lp; o <- h$options
  .hmip_b2_highs_invalidate(h)
  opts <- list(ptol = o$primal_feasibility_tolerance, dtol = o$dual_feasibility_tolerance,
               iteration_limit = o$simplex_iteration_limit, objective_bound = o$objective_bound,
               strategy = if (identical(as.integer(o$simplex_strategy), 4L)) "primal" else "dual")
  use_presolve <- FALSE
  if (lp$num_row_ == 0L) {
    sol <- .hs_unconstrained(lp$col_cost_, lp$col_lower_, lp$col_upper_, opts$ptol, opts$dtol)
    sol$row_status <- character()
    .hmip_b2_highs_clear_ekk(h)
  } else {
    use_presolve <- !identical(o$presolve, "off") && !isTRUE(h$basis$valid)
    if (!h$scale_tried) { h$scale <- .hs_scale(lp$A, lp$col_cost_); h$scale_tried <- TRUE }
    if (use_presolve) {
      sol <- .hs_highs_run(lp$col_cost_, lp$A, lp$row_lower_, lp$row_upper_, lp$col_lower_, lp$col_upper_, opts)
    } else {
      opts$fixed_scale <- h$scale
      if (!is.null(h$ek) && h$ek_valid) {
        ek <- h$ek
        ek$iteration_count <- 0L
        ek$iteration_limit <- opts$iteration_limit
        ek$objective_bound <- opts$objective_bound
        ek$dual_ray <- NULL
        opts$ek <- ek
      } else if (isTRUE(h$basis$valid)) {
        opts$basis <- list(col_status = h$basis$col_status, row_status = h$basis$row_status)
      }
      sol <- .hs_lp_solve(lp$col_cost_, lp$A, lp$row_lower_, lp$row_upper_, lp$col_lower_, lp$col_upper_, opts)
    }
    if (!is.null(sol$ek)) {
      h$ek <- sol$ek
      h$ek_valid <- !isTRUE(sol$ek$moved_unscaled) && identical(sol$scale, h$scale)
    } else .hmip_b2_highs_clear_ekk(h)
  }
  st <- sol$status
  if (st == "UnboundedOrInfeasible") st <- "Infeasible"
  h$model_status <- st
  if (st == "SolveError") return("Error")
  if (!is.null(sol$col_status))
    h$basis <- list(valid = TRUE, col_status = sol$col_status, row_status = sol$row_status)
  h$solution <- list(value_valid = TRUE, dual_valid = TRUE, col_value = sol$col_value, col_dual = sol$col_dual,
                     row_value = sol$row_value, row_dual = sol$row_dual)
  k <- .hmip_b2_kkt(lp, h$solution, opts$ptol, opts$dtol)
  h$info <- c(list(valid = TRUE, simplex_iteration_count = sol$nit, objective_function_value = sol$objective,
                   basis_validity = 1L,
                   primal_solution_status = if (k$num_primal_infeasibilities == 0L) 2L else 1L,
                   dual_solution_status = if (k$num_dual_infeasibilities == 0L) 2L else 1L), k)
  if (st %in% c("Optimal", "Infeasible", "Unbounded", "ObjectiveBound")) "Ok" else "Warning"
}

# btran of rhs (indexed by pivot row) in the space of the unscaled LP,
# HSimplexNla::btran: applyBasisMatrixColScale, the factor BTRAN, then
# applyBasisMatrixRowScale. Returns list(count, index (0-based), array).
.hmip_b2_btran <- function(h, rhs) {
  ek <- h$ek
  sc <- if (isTRUE(ek$moved_unscaled)) NULL else h$scale
  if (!is.null(sc)) rhs <- rhs * c(sc$col, 1 / sc$row)[ek$basicIndex]
  y <- .hs_btran(ek, rhs, "row_ep", pack = TRUE)
  idx <- ek$last_index
  if (is.null(idx)) idx <- which(y != 0)
  if (!is.null(sc)) y <- y * sc$row
  idx <- idx[y[idx] != 0]
  list(count = length(idx), index = idx - 1L, array = y)
}

.hmip_Highs_getBasicVariablesArray <- function(h) h$ek$basicIndex - 1L
# returns list(ok, count, index (0-based), array); ok is FALSE (C++ kError)
# without an invertible representation
.hmip_Highs_getBasisInverseRowSparse <- function(h, row) {
  if (!.hmip_Highs_hasInvert(h)) return(list(ok = FALSE, count = 0L, index = integer(), array = numeric()))
  e <- numeric(h$ek$m); e[row + 1L] <- 1
  c(list(ok = TRUE), .hmip_b2_btran(h, e))
}
# Highs::getDualEdgeWeights: the dual steepest edge weights of the basic
# rows (scaled space, as HEkk holds them), NULL without such weights
.hmip_Highs_getDualEdgeWeights <- function(h) {
  ek <- h$ek
  if (is.null(ek) || !isTRUE(ek$has_dse_weights)) return(NULL)
  ek$dual_edge_weight
}
.hmip_Highs_getDualRaySparse <- function(h) {
  ek <- h$ek
  if (is.null(ek) || is.null(ek$dual_ray)) return(list(has = FALSE, row_ep = NULL))
  e <- numeric(ek$m); e[ek$dual_ray$row] <- ek$dual_ray$sign
  list(has = TRUE, row_ep = .hmip_b2_btran(h, e))
}

# HEkk::putIterate and HEkk::getIterate: the basis, the invert and the dual
# edge weights
.hmip_Highs_putIterate <- function(h) {
  ek <- h$ek
  h$iterate <- list(basicIndex = ek$basicIndex, nonbasicFlag = ek$nonbasicFlag, nonbasicMove = ek$nonbasicMove,
                    fa = list2env(as.list(ek$fa, all.names = TRUE), parent = emptyenv()),
                    update_count = ek$update_count, has_dse_weights = isTRUE(ek$has_dse_weights),
                    dual_edge_weight = ek$dual_edge_weight)
  invisible("Ok")
}
.hmip_Highs_getIterate <- function(h) {
  it <- h$iterate
  if (is.null(it) || is.null(h$ek)) return(invisible("Error"))
  ek <- h$ek
  ek$basicIndex <- it$basicIndex; ek$nonbasicFlag <- it$nonbasicFlag; ek$nonbasicMove <- it$nonbasicMove
  ek$fa <- list2env(as.list(it$fa, all.names = TRUE), parent = emptyenv())
  ek$update_count <- it$update_count
  if (it$has_dse_weights) ek$dual_edge_weight <- it$dual_edge_weight else ek$has_dse_weights <- FALSE
  h$ek_valid <- TRUE
  h$basis <- .hmip_b2_ek_basis(ek)
  .hmip_b2_highs_invalidate(h)
  invisible("Ok")
}

# HEkkDual::bailoutOnDualObjective and reachedExactObjectiveBound. The port
# checks after each rebuild of phase 2, where the dual objective is fresh
# and update_count is zero; HiGHS also checks the updated dual objective
# after each iteration.
.hs_bailout_on_dual_objective <- function(ek, ds) {
  if (ek$solve_bailout) return(TRUE)
  bound <- ek$objective_bound %||% Inf
  if (ds$solve_phase != 2L) return(FALSE)
  if (is.null(ek$dual_objective_value) || !(ek$dual_objective_value > bound)) return(FALSE)
  # computeExactDualObjectiveValue
  n <- ek$n; m <- ek$m
  cost <- ek$lp$col_cost
  dual_col <- numeric(m)
  bc <- ek$basicIndex <= n
  dual_col[bc] <- cost[ek$basicIndex[bc]]
  dual_row <- numeric(n)
  if (any(dual_col != 0)) {
    dual_col <- .hs_btran(ek, dual_col, "dual_col")
    dual_row <- as.numeric(crossprod(ek$A, dual_col))
  }
  small <- 1e-9
  obj <- c(0, 0)
  for (j in seq_len(n)) {
    if (ek$nonbasicFlag[j] == 0L) next
    d <- cost[j] - dual_row[j]
    av <- if (d > small) ek$lp$col_lower[j] else if (d < -small) ek$lp$col_upper[j] else ek$workValue[j]
    if (is.infinite(av)) return(FALSE)
    obj <- .cd_add(obj, av * d)
  }
  for (i in seq_len(m)) {
    if (ek$nonbasicFlag[n + i] == 0L) next
    d <- dual_col[i]
    av <- if (d > small) ek$lp$row_lower[i] else if (d < -small) ek$lp$row_upper[i] else -ek$workValue[n + i]
    if (is.infinite(av)) return(FALSE)
    obj <- .cd_add(obj, av * d)
  }
  if (!(.cd_val(obj) > bound)) return(FALSE)
  if (isTRUE(ek$costs_perturbed) || isTRUE(ek$costs_shifted)) .hs_init_cost(ek, "dual", FALSE)
  ek$workDual <- c(ek$workCost[seq_len(n)] - dual_row, -dual_col)
  ds$force_phase2 <- FALSE
  .hs_correct_dual(ek, ds)
  ek$model_status <- "ObjectiveBound"
  ek$solve_bailout <- TRUE
  TRUE
}

# ---------------------------------------------------------------------------
# HighsLpRelaxation (mip/HighsLpRelaxation.h/.cpp)

.hmip_b2_src_SolveLp <- 9L
.hmip_b2_src_Unbounded <- 11L

# std::round (half away from zero) and fractionality() of HighsIntegers.h
.hmip_b2_round <- function(x) {
  r <- trunc(x)
  ifelse(abs(x - r) >= 0.5, r + sign(x), r)
}
.hmip_b2_frac <- function(x) abs(x - .hmip_b2_round(x))

.hmip_new_HighsLpRelaxation <- function(mipsolver) {
  lp <- new.env(parent = emptyenv())
  lp$mipsolver <- mipsolver
  lp$lpsolver <- .hmip_new_Highs()
  o <- mipsolver$options_mip_
  .hmip_Highs_setOptionValue(lp$lpsolver, "output_flag", FALSE)
  .hmip_Highs_setOptionValue(lp$lpsolver, "random_seed", o$random_seed)
  .hmip_Highs_setOptionValue(lp$lpsolver, "primal_feasibility_tolerance", o$mip_feasibility_tolerance)
  .hmip_Highs_setOptionValue(lp$lpsolver, "dual_feasibility_tolerance", o$mip_feasibility_tolerance * 0.1)
  lp$lprows <- list(origin = character(), index = integer(), age = integer())
  lp$fractionalints <- list(first = integer(), second = numeric())
  lp$dualproofvals <- numeric(); lp$dualproofinds <- integer(); lp$dualproofrhs <- Inf
  lp$hasdualproof <- FALSE
  lp$status <- "NotSet"
  lp$numlpiters <- 0; lp$avgSolveIters <- 0; lp$numSolved <- 0; lp$epochs <- 0
  lp$maxNumFractional <- 0L; lp$lastAgeCall <- 0
  lp$objective <- -Inf
  lp$basischeckpoint <- NULL
  lp$currentbasisstored <- FALSE
  lp$adjustSymBranchingCol <- TRUE
  lp$solved_first_lp <- TRUE
  lp
}

# the copy constructor (LP based heuristics take a local LP relaxation)
.hmip_HighsLpRelaxation_copy <- function(other) {
  lp <- new.env(parent = emptyenv())
  lp$mipsolver <- other$mipsolver
  lp$lprows <- other$lprows
  lp$fractionalints <- other$fractionalints
  lp$objective <- other$objective
  lp$basischeckpoint <- other$basischeckpoint
  lp$currentbasisstored <- other$currentbasisstored
  lp$adjustSymBranchingCol <- other$adjustSymBranchingCol
  lp$lpsolver <- .hmip_new_Highs()
  .hmip_Highs_setOptionValue(lp$lpsolver, "output_flag", FALSE)
  .hmip_Highs_passOptions(lp$lpsolver, .hmip_Highs_getOptions(other$lpsolver))
  .hmip_Highs_passModel(lp$lpsolver, .hmip_Highs_getLp(other$lpsolver))
  .hmip_Highs_setBasis(lp$lpsolver, .hmip_Highs_getBasis(other$lpsolver))
  lp$dualproofvals <- numeric(); lp$dualproofinds <- integer(); lp$dualproofrhs <- Inf
  lp$hasdualproof <- FALSE
  lp$status <- "NotSet"
  lp$numlpiters <- 0; lp$avgSolveIters <- 0; lp$numSolved <- 0; lp$epochs <- 0
  lp$maxNumFractional <- 0L; lp$lastAgeCall <- 0
  lp$objective <- -Inf
  lp$solved_first_lp <- TRUE
  lp
}

.hmip_HighsLpRelaxation_getCutPool <- function(lp) {
  l <- .hmip_Highs_getLp(lp$lpsolver)
  num_col <- l$num_col_
  num_model_row <- lp$mipsolver$model_$num_row_
  num_cut <- l$num_row_ - num_model_row
  cut_row_index <- rep(-1L, l$num_row_)
  cut_lower <- numeric(num_cut); cut_upper <- numeric(num_cut)
  cut_num <- 0L
  for (iRow in seq_len(l$num_row_)) {
    if (lp$lprows$origin[iRow] != "cut") next
    cut_row_index[iRow] <- cut_num
    cut_lower[cut_num + 1L] <- l$row_lower_[iRow]
    cut_upper[cut_num + 1L] <- l$row_upper_[iRow]
    cut_num <- cut_num + 1L
  }
  am <- l$a_matrix_
  len <- integer(num_cut)
  for (iCol in seq_len(num_col)) for (iEl in seq_len(am$start_[iCol + 1L] - am$start_[iCol]) + am$start_[iCol]) {
    iCut <- cut_row_index[am$index_[iEl] + 1L]
    if (iCut >= 0L) len[iCut + 1L] <- len[iCut + 1L] + 1L
  }
  start <- c(0L, cumsum(len))
  pos <- start[seq_len(num_cut)]
  index <- integer(start[num_cut + 1L]); value <- numeric(start[num_cut + 1L])
  for (iCol in seq_len(num_col)) for (iEl in seq_len(am$start_[iCol + 1L] - am$start_[iCol]) + am$start_[iCol]) {
    iCut <- cut_row_index[am$index_[iEl] + 1L]
    if (iCut >= 0L) {
      index[pos[iCut + 1L] + 1L] <- iCol - 1L
      value[pos[iCut + 1L] + 1L] <- am$value_[iEl]
      pos[iCut + 1L] <- pos[iCut + 1L] + 1L
    }
  }
  list(num_col = num_col, num_cut = num_cut, cut_lower = cut_lower, cut_upper = cut_upper,
       cut_matrix = list(num_col_ = num_col, num_row_ = num_cut, format_ = "rowwise", start_ = start,
                         index_ = index, value_ = value))
}

# LpRow::get, getRowLen, isIntegral, getMaxAbsVal
.hmip_b2_lprow_get <- function(lp, row) {
  md <- lp$mipsolver$mipdata_
  if (lp$lprows$origin[row + 1L] == "cut")
    return(.hmip_HighsCutPool_getCut(md$cutpool, lp$lprows$index[row + 1L]))
  r <- .hmip_HighsMipSolverData_getRow(md, lp$lprows$index[row + 1L])
  list(len = r$rowlen, inds = r$rowinds, vals = r$rowvals)
}
.hmip_HighsLpRelaxation_getRow <- function(lp, row) .hmip_b2_lprow_get(lp, row)
.hmip_HighsLpRelaxation_getRowLen <- function(lp, row) {
  md <- lp$mipsolver$mipdata_; idx <- lp$lprows$index[row + 1L]
  if (lp$lprows$origin[row + 1L] == "cut") .hmip_HighsCutPool_getRowLength(md$cutpool, idx)
  else md$ARstart_[idx + 2L] - md$ARstart_[idx + 1L]
}
.hmip_HighsLpRelaxation_isRowIntegral <- function(lp, row) {
  md <- lp$mipsolver$mipdata_; idx <- lp$lprows$index[row + 1L]
  if (lp$lprows$origin[row + 1L] == "cut") .hmip_HighsCutPool_cutIsIntegral(md$cutpool, idx)
  else md$rowintegral[idx + 1L] != 0
}
.hmip_HighsLpRelaxation_getMaxAbsRowVal <- function(lp, row) {
  md <- lp$mipsolver$mipdata_; idx <- lp$lprows$index[row + 1L]
  if (lp$lprows$origin[row + 1L] == "cut") .hmip_HighsCutPool_getMaxAbsCutCoef(md$cutpool, idx)
  else md$maxAbsRowCoef[idx + 1L]
}

.hmip_HighsLpRelaxation_rowLower <- function(lp, row) .hmip_Highs_getLp(lp$lpsolver)$row_lower_[row + 1L]
.hmip_HighsLpRelaxation_rowUpper <- function(lp, row) .hmip_Highs_getLp(lp$lpsolver)$row_upper_[row + 1L]
.hmip_HighsLpRelaxation_slackLower <- function(lp, row) {
  md <- lp$mipsolver$mipdata_; idx <- lp$lprows$index[row + 1L]
  if (lp$lprows$origin[row + 1L] == "cut") return(.hmip_HighsDomain_getMinCutActivity(md$domain, md$cutpool, idx))
  rowlower <- .hmip_HighsLpRelaxation_rowLower(lp, row)
  if (rowlower != -Inf) return(rowlower)
  .hmip_HighsDomain_getMinActivity(md$domain, idx)
}
.hmip_HighsLpRelaxation_slackUpper <- function(lp, row) {
  rowupper <- .hmip_HighsLpRelaxation_rowUpper(lp, row)
  if (lp$lprows$origin[row + 1L] == "cut") return(rowupper)
  if (rowupper != Inf) return(rowupper)
  md <- lp$mipsolver$mipdata_
  .hmip_HighsDomain_getMaxActivity(md$domain, lp$lprows$index[row + 1L])
}
.hmip_HighsLpRelaxation_colLower <- function(lp, col) {
  l <- .hmip_Highs_getLp(lp$lpsolver)
  if (col < l$num_col_) l$col_lower_[col + 1L] else .hmip_HighsLpRelaxation_slackLower(lp, col - l$num_col_)
}
.hmip_HighsLpRelaxation_colUpper <- function(lp, col) {
  l <- .hmip_Highs_getLp(lp$lpsolver)
  if (col < l$num_col_) l$col_upper_[col + 1L] else .hmip_HighsLpRelaxation_slackUpper(lp, col - l$num_col_)
}
.hmip_HighsLpRelaxation_isColIntegral <- function(lp, col) {
  n <- .hmip_Highs_getNumCol(lp$lpsolver)
  if (col < n) lp$mipsolver$model_$integrality_[col + 1L] != 0L
  else .hmip_HighsLpRelaxation_isRowIntegral(lp, col - n)
}
.hmip_HighsLpRelaxation_solutionValue <- function(lp, col) {
  n <- .hmip_Highs_getNumCol(lp$lpsolver)
  s <- lp$lpsolver$solution
  if (col < n) s$col_value[col + 1L] else s$row_value[col - n + 1L]
}

.hmip_HighsLpRelaxation_setAdjustSymmetricBranchingCol <- function(lp, x) lp$adjustSymBranchingCol <- x
.hmip_HighsLpRelaxation_getAvgSolveIters <- function(lp) lp$avgSolveIters
.hmip_HighsLpRelaxation_getLp <- function(lp) .hmip_Highs_getLp(lp$lpsolver)
.hmip_HighsLpRelaxation_getSolution <- function(lp) .hmip_Highs_getSolution(lp$lpsolver)
.hmip_HighsLpRelaxation_getStatus <- function(lp) lp$status
.hmip_HighsLpRelaxation_getSolverInfo <- function(lp) .hmip_Highs_getInfo(lp$lpsolver)
.hmip_HighsLpRelaxation_getNumLpIterations <- function(lp) lp$numlpiters
.hmip_HighsLpRelaxation_integerFeasible <- function(lp)
  lp$status %in% c("Optimal", "UnscaledPrimalFeasible") && length(lp$fractionalints$first) == 0L
.hmip_HighsLpRelaxation_scaledOptimal <- function(status)
  status %in% c("Optimal", "UnscaledDualFeasible", "UnscaledPrimalFeasible", "UnscaledInfeasible")
.hmip_HighsLpRelaxation_unscaledPrimalFeasible <- function(status) status %in% c("Optimal", "UnscaledPrimalFeasible")
.hmip_HighsLpRelaxation_unscaledDualFeasible <- function(status) status %in% c("Optimal", "UnscaledDualFeasible")
.hmip_HighsLpRelaxation_storeBasis <- function(lp) {
  b <- .hmip_Highs_getBasis(lp$lpsolver)
  if (!lp$currentbasisstored && isTRUE(b$valid)) {
    lp$basischeckpoint <- b
    lp$currentbasisstored <- TRUE
  }
}
.hmip_HighsLpRelaxation_getStoredBasis <- function(lp) lp$basischeckpoint
.hmip_HighsLpRelaxation_setStoredBasis <- function(lp, basis) {
  lp$basischeckpoint <- basis
  lp$currentbasisstored <- FALSE
}
.hmip_HighsLpRelaxation_getMipSolver <- function(lp) lp$mipsolver
.hmip_HighsLpRelaxation_getNumModelRows <- function(lp) lp$mipsolver$model_$num_row_
.hmip_HighsLpRelaxation_numRows <- function(lp) .hmip_Highs_getNumRow(lp$lpsolver)
.hmip_HighsLpRelaxation_numCols <- function(lp) .hmip_Highs_getNumCol(lp$lpsolver)
.hmip_HighsLpRelaxation_numNonzeros <- function(lp) .hmip_Highs_getNumNz(lp$lpsolver)
.hmip_HighsLpRelaxation_getDualProof <- function(lp)
  list(inds = lp$dualproofinds, vals = lp$dualproofvals, rhs = lp$dualproofrhs, len = length(lp$dualproofinds))
.hmip_HighsLpRelaxation_getLpSolver <- function(lp) lp$lpsolver
.hmip_HighsLpRelaxation_getFractionalIntegers <- function(lp) lp$fractionalints
.hmip_HighsLpRelaxation_getObjective <- function(lp) lp$objective
.hmip_HighsLpRelaxation_setIterationLimit <- function(lp, limit = .hmip_b2_iinf)
  .hmip_Highs_setOptionValue(lp$lpsolver, "simplex_iteration_limit", limit)
.hmip_HighsLpRelaxation_setSolvedFirstLp <- function(lp, x) lp$solved_first_lp <- x

.hmip_HighsLpRelaxation_loadModel <- function(lp) {
  ms <- lp$mipsolver
  m <- ms$model_
  dom <- ms$mipdata_$domain
  lpmodel <- list(num_col_ = m$num_col_, num_row_ = m$num_row_, col_cost_ = m$col_cost_,
                  col_lower_ = dom$col_lower_, col_upper_ = dom$col_upper_, row_lower_ = m$row_lower_,
                  row_upper_ = m$row_upper_, a_matrix_ = m$a_matrix_, offset_ = 0)
  lp$lprows <- list(origin = rep("model", m$num_row_), index = seq_len(m$num_row_) - 1L,
                    age = integer(m$num_row_))
  .hmip_Highs_clearSolver(lp$lpsolver)
  .hmip_Highs_passModel(lp$lpsolver, lpmodel)
  invisible()
}

.hmip_HighsLpRelaxation_resetToGlobalDomain <- function(lp) {
  ms <- lp$mipsolver
  .hmip_Highs_changeColsBoundsRange(lp$lpsolver, 0L, ms$model_$num_col_ - 1L,
                                    ms$mipdata_$domain$col_lower_, ms$mipdata_$domain$col_upper_)
}

.hmip_HighsLpRelaxation_computeBasicDegenerateDuals <- function(lp, threshold, localdom = NULL) {
  h <- lp$lpsolver
  if (!.hmip_Highs_hasInvert(h)) return(invisible())
  md <- lp$mipsolver$mipdata_
  feastol <- md$feastol
  l <- .hmip_Highs_getLp(h)
  basis <- .hmip_Highs_getBasis(h)
  sol <- h$solution
  k <- 0L
  for (col in md$integral_cols) {
    j <- col + 1L
    if (basis$col_status[j] != "basic") next
    lb <- l$col_lower_[j]; ub <- l$col_upper_[j]
    if (ub - lb < feastol) next
    if (sol$col_value[j] - lb < ub - sol$col_value[j]) {
      if (sol$col_value[j] > lb + feastol) next
      sol$col_dual[j] <- 1
    } else {
      if (sol$col_value[j] < ub - feastol) next
      sol$col_dual[j] <- -1
    }
    k <- k + 1L
  }
  h$solution <- sol
  if (k == 0L) return(invisible())
  num_col <- l$num_col_
  basicIndex <- .hmip_Highs_getBasicVariablesArray(h)
  row <- 0L
  while (k > 0L) {
    var <- basicIndex[row + 1L]
    row <- row + 1L
    if (var >= num_col) next
    if (abs(h$solution$col_dual[var + 1L]) != 1) next
    k <- k - 1L
    row_ep <- .hmip_Highs_getBasisInverseRowSparse(h, row - 1L)
    sign <- h$solution$col_dual[var + 1L]
    h$solution$col_dual[var + 1L] <- 0
    sol <- h$solution
    degenerateColDual <- Inf
    for (iRow in row_ep$index) {
      lb <- l$row_lower_[iRow + 1L]; ub <- l$row_upper_[iRow + 1L]
      if (lb == ub) next
      dual <- sol$row_dual[iRow + 1L]
      val <- -sign * row_ep$array[iRow + 1L]
      if (val > 0) {
        if (sol$row_value[iRow + 1L] - lb > feastol) {
          degenerateColDual <- min(degenerateColDual, -dual / val)
          if (degenerateColDual < threshold) break
        }
      } else {
        if (ub - sol$row_value[iRow + 1L] > feastol) {
          degenerateColDual <- min(degenerateColDual, -dual / val)
          if (degenerateColDual < threshold) break
        }
      }
    }
    if (degenerateColDual < threshold) next
    row_ap <- .hmip_b2_svs(num_col)
    for (iRow in row_ep$index) {
      r <- .hmip_HighsLpRelaxation_getRow(lp, iRow)
      for (q in seq_len(r$len)) .hmip_b2_svs_add(row_ap, r$inds[q], row_ep$array[iRow + 1L] * r$vals[q])
    }
    eps <- md$epsilon
    .hmip_b2_svs_cleanup(row_ap, function(i, v) abs(v) <= eps)
    for (iCol in row_ap$nz) {
      if (iCol == var) next
      lb <- l$col_lower_[iCol + 1L]; ub <- l$col_upper_[iCol + 1L]
      if (lb == ub) next
      dual <- sol$col_dual[iCol + 1L]
      val <- sign * .hmip_b2_svs_get(row_ap, iCol)
      if (val > md$epsilon) {
        if (sol$col_value[iCol + 1L] - lb > feastol) {
          degenerateColDual <- min(degenerateColDual, -dual / val)
          if (degenerateColDual < threshold) break
        }
      } else if (val < -md$epsilon) {
        if (ub - sol$col_value[iCol + 1L] > feastol) {
          degenerateColDual <- min(degenerateColDual, -dual / val)
          if (degenerateColDual < threshold) break
        }
      }
    }
    if (degenerateColDual < threshold) next
    if (degenerateColDual == Inf && !is.null(localdom)) {
      rhs <- c(0, 0)
      for (iRow in row_ep$index) {
        lb <- l$row_lower_[iRow + 1L]; ub <- l$row_upper_[iRow + 1L]
        val <- sign * row_ep$array[iRow + 1L]
        if (ub == lb || val > md$epsilon) rhs <- .cd_add(rhs, val * ub)
        else if (val < -md$epsilon) rhs <- .cd_add(rhs, val * lb)
        else rhs <- .cd_add(rhs, val * sol$row_value[iRow + 1L])
      }
      vals <- vapply(row_ap$nz, function(i) sign * .hmip_b2_svs_get(row_ap, i), 0)
      domchg <- if (sign == 1) list(boundval = l$col_lower_[var + 1L], column = var, boundtype = .hmip_b2_upper)
                else list(boundval = l$col_upper_[var + 1L], column = var, boundtype = .hmip_b2_lower)
      lp$dualproofvals <- vals
      .hmip_HighsDomain_conflictAnalyzeReconvergence(localdom, domchg, row_ap$nz, vals, length(row_ap$nz),
                                                     .cd_val(rhs), md$conflictPool)
      next
    }
    h$solution$col_dual[var + 1L] <- sign * degenerateColDual
  }
  invisible()
}

.hmip_HighsLpRelaxation_computeBestEstimate <- function(lp, ps) {
  estimate <- c(lp$objective, 0)
  fi <- lp$fractionalints
  if (length(fi$first)) {
    md <- lp$mipsolver$mipdata_
    increase <- c(0, 0)
    offset <- md$feastol * max(abs(lp$objective), 1) / length(md$integral_cols)
    for (q in seq_along(fi$first))
      increase <- .cd_add(increase, min(.hmip_HighsPseudocost_getPseudocostUp(ps, fi$first[q], fi$second[q], offset),
                                        .hmip_HighsPseudocost_getPseudocostDown(ps, fi$first[q], fi$second[q], offset)))
    estimate <- .cd_add(estimate, .cd_val(increase))
  }
  .cd_val(estimate)
}

.hmip_HighsLpRelaxation_computeLPDegneracy <- function(lp, localdomain) {
  h <- lp$lpsolver
  sol <- h$solution; basis <- h$basis
  if (!isTRUE(sol$dual_valid) || !isTRUE(basis$valid)) return(1)
  dualFeasTol <- h$info$max_dual_infeasibility
  l <- .hmip_Highs_getLp(h)
  nr <- l$num_row_; nc <- l$num_col_
  numFixedRows <- 0L; numInequalities <- 0L; numBasicEqualities <- 0L
  for (i in seq_len(nr)) {
    if (l$row_lower_[i] != l$row_upper_[i]) {
      numInequalities <- numInequalities + 1L
      if (basis$row_status[i] != "basic" && abs(sol$row_dual[i]) > dualFeasTol) numFixedRows <- numFixedRows + 1L
    } else numBasicEqualities <- numBasicEqualities + (basis$row_status[i] == "basic")
  }
  numAlreadyFixedCols <- 0L; numFixedCols <- 0L
  for (i in seq_len(nc)) {
    if (basis$col_status[i] != "basic") {
      if (abs(sol$col_dual[i]) > dualFeasTol) numFixedCols <- numFixedCols + 1L
      else if (localdomain$col_lower_[i] == localdomain$col_upper_[i]) numAlreadyFixedCols <- numAlreadyFixedCols + 1L
    }
  }
  base <- nc - numAlreadyFixedCols + numInequalities + numBasicEqualities - nr
  degenerateColumnShare <- if (base > 0) 1 - (numFixedCols + numFixedRows) / base else 0
  varConsRatio <- if (nr > 0) (nc + numInequalities + numBasicEqualities - numFixedCols - numFixedRows -
                                 numAlreadyFixedCols) / nr else 1
  fac1 <- if (degenerateColumnShare < 0.8) 1 else 10^(10 * (degenerateColumnShare - 0.7))
  fac2 <- if (varConsRatio < 2) 1 else 10 * varConsRatio
  fac1 * fac2
}

.hmip_HighsLpRelaxation_addCuts <- function(lp, cutset) {
  numcuts <- .hmip_HighsCutSet_numCuts(cutset)
  if (numcuts > 0L) {
    lp$status <- "NotSet"
    lp$currentbasisstored <- FALSE
    lp$basischeckpoint <- NULL
    lp$lprows$origin <- c(lp$lprows$origin, rep("cut", numcuts))
    lp$lprows$index <- c(lp$lprows$index, as.integer(cutset$cutindices[seq_len(numcuts)]))
    lp$lprows$age <- c(lp$lprows$age, integer(numcuts))
    .hmip_Highs_addRows(lp$lpsolver, numcuts, cutset$lower_, cutset$upper_, length(cutset$ARvalue_),
                        cutset$ARstart_, cutset$ARindex_, cutset$ARvalue_)
    .hmip_HighsCutSet_clear(cutset)
  }
  invisible()
}

.hmip_HighsLpRelaxation_removeObsoleteRows <- function(lp, notifyPool = TRUE) {
  nlprows <- .hmip_HighsLpRelaxation_numRows(lp)
  nummodelrows <- .hmip_HighsLpRelaxation_getNumModelRows(lp)
  deletemask <- integer(nlprows)
  ndelcuts <- 0L
  rs <- .hmip_Highs_getBasis(lp$lpsolver)$row_status
  for (i in seq_len(nlprows - nummodelrows) + nummodelrows) {
    if (rs[i] == "basic") {
      ndelcuts <- ndelcuts + 1L
      deletemask[i] <- 1L
      if (notifyPool) .hmip_HighsCutPool_lpCutRemoved(lp$mipsolver$mipdata_$cutpool, lp$lprows$index[i])
    }
  }
  .hmip_HighsLpRelaxation_removeCuts(lp, ndelcuts, deletemask)
}

# removeCuts(ndelcuts, deletemask)
.hmip_HighsLpRelaxation_removeCuts <- function(lp, ndelcuts, deletemask) {
  if (ndelcuts > 0L) {
    basis <- .hmip_Highs_getBasis(lp$lpsolver)
    nlprows <- .hmip_Highs_getNumRow(lp$lpsolver)
    deletemask <- .hmip_Highs_deleteRowsMask(lp$lpsolver, deletemask)
    for (i in seq_len(nlprows - lp$mipsolver$model_$num_row_) + lp$mipsolver$model_$num_row_) {
      d <- deletemask[i]
      if (d >= 0L) {
        lp$lprows$origin[d + 1L] <- lp$lprows$origin[i]
        lp$lprows$index[d + 1L] <- lp$lprows$index[i]
        lp$lprows$age[d + 1L] <- lp$lprows$age[i]
        basis$row_status[d + 1L] <- basis$row_status[i]
      }
    }
    keep <- seq_len(length(basis$row_status) - ndelcuts)
    basis$row_status <- basis$row_status[keep]
    keep <- seq_len(length(lp$lprows$index) - ndelcuts)
    lp$lprows <- list(origin = lp$lprows$origin[keep], index = lp$lprows$index[keep], age = lp$lprows$age[keep])
    .hmip_Highs_setBasis(lp$lpsolver, basis)
    .hmip_Highs_run(lp$lpsolver)
  }
  invisible(deletemask)
}

# the argumentless removeCuts()
.hmip_HighsLpRelaxation_removeAllCuts <- function(lp) {
  nlprows <- .hmip_Highs_getNumRow(lp$lpsolver)
  modelrows <- lp$mipsolver$model_$num_row_
  .hmip_Highs_deleteRowsRange(lp$lpsolver, modelrows, nlprows - 1L)
  for (i in seq_len(nlprows - modelrows) + modelrows)
    if (lp$lprows$origin[i] == "cut") .hmip_HighsCutPool_lpCutRemoved(lp$mipsolver$mipdata_$cutpool, lp$lprows$index[i])
  keep <- seq_len(modelrows)
  lp$lprows <- list(origin = lp$lprows$origin[keep], index = lp$lprows$index[keep], age = lp$lprows$age[keep])
  invisible()
}

.hmip_b2_aging_skip <- function(lp) {
  info <- .hmip_Highs_getInfo(lp$lpsolver)
  info$basis_validity == 0L || info$max_dual_infeasibility > lp$mipsolver$mipdata_$feastol ||
    !isTRUE(lp$lpsolver$solution$dual_valid)
}

.hmip_HighsLpRelaxation_performAging <- function(lp, deleteRows = FALSE) {
  if (.hmip_b2_aging_skip(lp)) return(invisible())
  if (deleteRows) {
    agelimit <- lp$mipsolver$options_mip_$mip_lp_age_limit
    lp$epochs <- lp$epochs + 1
    if (lp$epochs %% max(agelimit %/% 2L, 2L) != 0) agelimit <- .hmip_b2_iinf
    else if (lp$epochs < agelimit) agelimit <- lp$epochs
  } else {
    if (lp$lastAgeCall == lp$numlpiters) return(invisible())
    agelimit <- .hmip_b2_iinf
  }
  lp$lastAgeCall <- lp$numlpiters
  nlprows <- .hmip_HighsLpRelaxation_numRows(lp)
  nummodelrows <- .hmip_HighsLpRelaxation_getNumModelRows(lp)
  deletemask <- integer(nlprows)
  ndelcuts <- 0L
  rs <- .hmip_Highs_getBasis(lp$lpsolver)$row_status
  rd <- lp$lpsolver$solution$row_dual
  dtol <- lp$lpsolver$options$dual_feasibility_tolerance
  for (i in seq_len(nlprows - nummodelrows) + nummodelrows) {
    if (rs[i] == "basic") {
      lp$lprows$age[i] <- lp$lprows$age[i] + as.integer(deleteRows || lp$lprows$age[i] != 0L)
      if (lp$lprows$age[i] > agelimit) {
        ndelcuts <- ndelcuts + 1L
        deletemask[i] <- 1L
        .hmip_HighsCutPool_lpCutRemoved(lp$mipsolver$mipdata_$cutpool, lp$lprows$index[i])
      }
    } else if (abs(rd[i]) > dtol) {
      lp$lprows$age[i] <- 0L
    }
  }
  .hmip_HighsLpRelaxation_removeCuts(lp, ndelcuts, deletemask)
  invisible()
}

.hmip_HighsLpRelaxation_resetAges <- function(lp) {
  if (.hmip_b2_aging_skip(lp)) return(invisible())
  nlprows <- .hmip_HighsLpRelaxation_numRows(lp)
  nummodelrows <- .hmip_HighsLpRelaxation_getNumModelRows(lp)
  rs <- .hmip_Highs_getBasis(lp$lpsolver)$row_status
  rd <- lp$lpsolver$solution$row_dual
  dtol <- lp$lpsolver$options$dual_feasibility_tolerance
  for (i in seq_len(nlprows - nummodelrows) + nummodelrows)
    if (rs[i] != "basic" && abs(rd[i]) > dtol) lp$lprows$age[i] <- 0L
  invisible()
}

.hmip_HighsLpRelaxation_flushDomain <- function(lp, domain, continuous = FALSE) {
  if (length(.hmip_HighsDomain_getChangedCols(domain))) {
    if (identical(domain, lp$mipsolver$mipdata_$domain)) continuous <- TRUE
    lp$currentbasisstored <- FALSE
    if (!continuous) .hmip_HighsDomain_removeContinuousChangedCols(domain)
    chg <- .hmip_HighsDomain_getChangedCols(domain)
    if (length(chg) == 0L) return(invisible())
    .hmip_Highs_changeColsBounds(lp$lpsolver, chg, domain$col_lower_[chg + 1L], domain$col_upper_[chg + 1L])
    .hmip_HighsDomain_clearChangedCols(domain)
  }
  invisible()
}

# computeDualProof (the active #else branch); returns list(ok, inds, vals, rhs)
.hmip_HighsLpRelaxation_computeDualProof <- function(lp, globaldomain, upperbound, extractCliques = TRUE) {
  ms <- lp$mipsolver; md <- ms$mipdata_
  row_dual <- lp$lpsolver$solution$row_dual
  l <- .hmip_Highs_getLp(lp$lpsolver)
  upper <- c(upperbound, 0)
  for (i in seq_len(l$num_row_)) {
    if (row_dual[i] > 0) {
      if (l$row_lower_[i] != -Inf) upper <- .cd_add(upper, -(row_dual[i] * l$row_lower_[i])) else row_dual[i] <- 0
    } else if (row_dual[i] < 0) {
      if (l$row_upper_[i] != Inf) upper <- .cd_add(upper, -(row_dual[i] * l$row_upper_[i])) else row_dual[i] <- 0
    }
  }
  inds <- integer(); vals <- numeric()
  am <- l$a_matrix_
  colval <- lp$lpsolver$solution$col_value
  small <- ms$options_mip_$small_matrix_value
  for (i in seq_len(l$num_col_)) {
    sum <- c(l$col_cost_[i], 0)
    for (j in seq_len(am$start_[i + 1L] - am$start_[i]) + am$start_[i]) {
      r <- am$index_[j] + 1L
      if (row_dual[r] == 0) next
      sum <- .cd_add(sum, -(am$value_[j] * row_dual[r]))
    }
    val <- .cd_val(sum)
    if (abs(val) <= small) next
    removeValue <- abs(val) <= md$feastol
    if (!removeValue && (globaldomain$col_lower_[i] == globaldomain$col_upper_[i] ||
                         ms$model_$integrality_[i] == 0L)) {
      removeValue <- if (val > 0) colval[i] - globaldomain$col_lower_[i] <= md$feastol
                     else globaldomain$col_upper_[i] - colval[i] <= md$feastol
    }
    if (removeValue) {
      if (val < 0) {
        if (globaldomain$col_upper_[i] == Inf) return(list(ok = FALSE, inds = inds, vals = vals, rhs = NA_real_))
        upper <- .cd_add(upper, -(val * globaldomain$col_upper_[i]))
      } else {
        if (globaldomain$col_lower_[i] == -Inf) return(list(ok = FALSE, inds = inds, vals = vals, rhs = NA_real_))
        upper <- .cd_add(upper, -(val * globaldomain$col_lower_[i]))
      }
      next
    }
    vals <- c(vals, val); inds <- c(inds, i - 1L)
  }
  rhs <- .cd_val(upper)
  t <- .hmip_HighsDomain_tightenCoefficients(globaldomain, inds, vals, length(inds), rhs)
  vals <- t$vals; rhs <- t$rhs
  if (extractCliques)
    .hmip_HighsCliqueTable_extractCliquesFromCut(md$cliquetable, ms, inds, vals, length(inds), rhs)
  list(ok = TRUE, inds = inds, vals = vals, rhs = rhs)
}

.hmip_HighsLpRelaxation_storeDualInfProof <- function(lp) {
  lp$hasdualproof <- FALSE
  h <- lp$lpsolver
  if (h$info$basis_validity == 0L) return(invisible())
  ray <- .hmip_Highs_getDualRaySparse(h)
  lp$hasdualproof <- ray$has
  if (!lp$hasdualproof) return(invisible())
  ms <- lp$mipsolver; md <- ms$mipdata_
  row_ep <- ray$row_ep
  lp$dualproofinds <- integer(); lp$dualproofvals <- numeric(); lp$dualproofrhs <- Inf
  l <- .hmip_Highs_getLp(h)
  upper <- c(0, 0)
  row_ap <- .hmip_b2_svs(l$num_col_)
  for (iRow in row_ep$index) {
    weight <- -row_ep$array[iRow + 1L]
    if (abs(weight) * .hmip_HighsLpRelaxation_getMaxAbsRowVal(lp, iRow) <= md$epsilon) next
    if (weight > 0) {
      if (l$row_upper_[iRow + 1L] == Inf) next
      upper <- .cd_add(upper, weight * l$row_upper_[iRow + 1L])
    } else {
      if (l$row_lower_[iRow + 1L] == -Inf) next
      upper <- .cd_add(upper, weight * l$row_lower_[iRow + 1L])
    }
    r <- .hmip_HighsLpRelaxation_getRow(lp, iRow)
    for (q in seq_len(r$len)) .hmip_b2_svs_add(row_ap, r$inds[q], weight * r$vals[q])
  }
  gd <- md$domain
  inds <- integer(); vals <- numeric()
  for (i in row_ap$nz) {
    val <- .hmip_b2_svs_get(row_ap, i)
    if (abs(val) <= md$epsilon) next
    removeValue <- abs(val) <= md$feastol
    j <- i + 1L
    if (!removeValue && (gd$col_lower_[j] == gd$col_upper_[j] || ms$model_$integrality_[j] == 0L)) {
      removeValue <- if (val > 0) l$col_lower_[j] - gd$col_lower_[j] <= md$feastol
                     else gd$col_upper_[j] - l$col_upper_[j] <= md$feastol
    }
    if (removeValue) {
      if (val < 0) {
        if (gd$col_upper_[j] == Inf) { lp$hasdualproof <- FALSE; return(invisible()) }
        upper <- .cd_add(upper, -(val * gd$col_upper_[j]))
      } else {
        if (gd$col_lower_[j] == -Inf) { lp$hasdualproof <- FALSE; return(invisible()) }
        upper <- .cd_add(upper, -(val * gd$col_lower_[j]))
      }
      next
    }
    vals <- c(vals, val); inds <- c(inds, i)
  }
  rhs <- .cd_val(upper)
  t <- .hmip_HighsDomain_tightenCoefficients(gd, inds, vals, length(inds), rhs)
  lp$dualproofinds <- inds; lp$dualproofvals <- t$vals; lp$dualproofrhs <- t$rhs
  .hmip_HighsCliqueTable_extractCliquesFromCut(md$cliquetable, ms, inds, t$vals, length(inds), t$rhs)
  invisible()
}

.hmip_HighsLpRelaxation_storeDualUBProof <- function(lp) {
  lp$dualproofinds <- integer(); lp$dualproofvals <- numeric()
  if (isTRUE(lp$lpsolver$solution$dual_valid)) {
    md <- lp$mipsolver$mipdata_
    r <- .hmip_HighsLpRelaxation_computeDualProof(lp, md$domain, md$upper_limit)
    lp$hasdualproof <- r$ok
    if (r$ok) { lp$dualproofinds <- r$inds; lp$dualproofvals <- r$vals; lp$dualproofrhs <- r$rhs }
  } else lp$hasdualproof <- FALSE
  if (!lp$hasdualproof) lp$dualproofrhs <- Inf
  invisible()
}

.hmip_HighsLpRelaxation_checkDualProof <- function(lp) {
  if (!lp$hasdualproof) return(TRUE)
  if (lp$dualproofrhs == Inf) return(FALSE)
  l <- .hmip_Highs_getLp(lp$lpsolver)
  viol <- c(-lp$dualproofrhs, 0)
  for (i in seq_along(lp$dualproofinds)) {
    col <- lp$dualproofinds[i] + 1L; v <- lp$dualproofvals[i]
    if (v > 0) {
      if (l$col_lower_[col] == -Inf) return(FALSE)
      viol <- .cd_add(viol, v * l$col_lower_[col])
    } else {
      if (l$col_upper_[col] == Inf) return(FALSE)
      viol <- .cd_add(viol, v * l$col_upper_[col])
    }
  }
  .cd_val(viol) > lp$mipsolver$mipdata_$feastol
}

.hmip_HighsLpRelaxation_computeDualInfProof <- function(lp, globaldomain) {
  if (!lp$hasdualproof) return(list(ok = FALSE))
  list(ok = TRUE, inds = lp$dualproofinds, vals = lp$dualproofvals, rhs = lp$dualproofrhs)
}

.hmip_HighsLpRelaxation_recoverBasis <- function(lp) {
  if (!is.null(lp$basischeckpoint)) {
    .hmip_Highs_setBasis(lp$lpsolver, lp$basischeckpoint, "HighsLpRelaxation::recoverBasis")
    lp$currentbasisstored <- TRUE
  }
  invisible()
}

.hmip_HighsLpRelaxation_setObjectiveLimit <- function(lp, objlim = Inf) {
  md <- lp$mipsolver$mipdata_
  of <- md$objectiveFunction
  offset <- if (.hmip_HighsObjectiveFunction_isIntegral(of)) 0.5 / .hmip_HighsObjectiveFunction_integralScale(of)
            else max(1000 * md$feastol, abs(objlim) * 1e-14)
  .hmip_Highs_setOptionValue(lp$lpsolver, "objective_bound", objlim + offset)
}

# run: the simplex solver always (mip_lp_solver "choose" never selects an
# interior point solver; the root analytic center is computed by B1).
# The time limit is not modeled.
.hmip_HighsLpRelaxation_run <- function(lp, resolve_on_error = TRUE) {
  h <- lp$lpsolver
  callstatus <- .hmip_Highs_run(h)
  lp$solved_first_lp <- TRUE
  info <- h$info
  itercount <- max(0, info$simplex_iteration_count)
  lp$numlpiters <- lp$numlpiters + itercount
  if (callstatus == "Error") {
    .hmip_Highs_clearSolver(h)
    if (resolve_on_error) {
      .hmip_Highs_setOptionValue(h, "simplex_strategy", 1L)
      .hmip_Highs_setOptionValue(h, "presolve", "on")
      retval <- .hmip_HighsLpRelaxation_run(lp, FALSE)
      .hmip_Highs_setOptionValue(h, "presolve", "off")
      return(retval)
    }
    .hmip_HighsLpRelaxation_recoverBasis(lp)
    return("Error")
  }
  ms <- lp$mipsolver; md <- ms$mipdata_
  model_status <- h$model_status
  solved <- function() {
    lp$numSolved <- lp$numSolved + 1
    lp$avgSolveIters <- lp$avgSolveIters + (itercount - lp$avgSolveIters) / lp$numSolved
  }
  switch(model_status,
    ObjectiveBound = {
      solved()
      .hmip_HighsLpRelaxation_storeDualUBProof(lp)
      "Infeasible"
    },
    Infeasible = {
      solved()
      .hmip_HighsLpRelaxation_storeDualInfProof(lp)
      "Infeasible"
    },
    Unbounded = {
      if (info$primal_solution_status == 2L)
        .hmip_HighsMipSolverData_trySolution(md, h$solution$col_value, .hmip_b2_src_Unbounded)
      "Unbounded"
    },
    Unknown = ,
    Optimal = {
      if (model_status == "Unknown" && info$basis_validity == 0L) return("Error")
      solved()
      feastol <- md$feastol
      if (info$max_primal_infeasibility <= feastol && info$max_dual_infeasibility <= feastol) return("Optimal")
      if (info$max_primal_infeasibility <= feastol) return("UnscaledPrimalFeasible")
      if (info$max_dual_infeasibility <= feastol) return("UnscaledDualFeasible")
      if (model_status == "Optimal") return("UnscaledInfeasible")
      "Error"
    },
    IterationLimit = {
      if (!isTRUE(ms$submip) && resolve_on_error)
        stop("not ported: HighsLpRelaxation::run IPM restart after the simplex iteration limit")
      "Error"
    },
    "Error")
}

.hmip_HighsLpRelaxation_resolveLp <- function(lp, domain = NULL) {
  lp$fractionalints <- list(first = integer(), second = numeric())
  ms <- lp$mipsolver; md <- ms$mipdata_
  h <- lp$lpsolver
  feastol <- md$feastol
  ncol <- ms$model_$num_col_
  cost <- ms$model_$col_cost_
  repeat {
    solveagain <- FALSE
    if (!is.null(domain)) .hmip_HighsLpRelaxation_flushDomain(lp, domain)
    lp$status <- .hmip_HighsLpRelaxation_run(lp)
    if (.hmip_HighsLpRelaxation_scaledOptimal(lp$status)) {
      fracints <- .hmip_new_HighsHashTable(lp$maxNumFractional)
      sol <- h$solution
      l <- .hmip_Highs_getLp(h)
      roundable <- TRUE
      for (i in md$integral_cols) {
        j <- i + 1L
        val <- max(min(sol$col_value[j], l$col_upper_[j]), l$col_lower_[j])
        if (.hmip_b2_frac(val) > feastol) {
          col <- i
          roundable <- roundable && (md$uplocks[col + 1L] == 0 || md$downlocks[col + 1L] == 0)
          subst <- .hmip_HighsCliqueTable_getSubstitution(md$cliquetable, col)
          while (!is.null(subst)) {
            rc <- .hmip_cv_col(subst$replace)
            if (l$col_lower_[rc + 1L] == l$col_upper_[rc + 1L]) {
              if (!is.null(domain)) {
                w <- if (.hmip_cv_val(subst$replace) != 0) l$col_lower_[rc + 1L] else 1 - l$col_lower_[rc + 1L]
                .hmip_HighsDomain_fixCol(domain, col, w)
                if (.hmip_HighsDomain_infeasible(domain)) {
                  lp$objective <- Inf
                  lp$status <- "Infeasible"
                  return(lp$status)
                }
              } else break
            }
            if (!is.null(domain)) {
              replace_val <- if (.hmip_cv_val(subst$replace) == 0) 1 - val else val
              if (replace_val < domain$col_lower_[rc + 1L] - feastol ||
                  replace_val > domain$col_upper_[rc + 1L] + feastol) break
            }
            col <- rc
            if (.hmip_cv_val(subst$replace) == 0) val <- 1 - val
            subst <- .hmip_HighsCliqueTable_getSubstitution(md$cliquetable, col)
          }
          if (lp$adjustSymBranchingCol)
            col <- .hmip_HighsSymmetries_getBranchingColumn(md$symmetries, l$col_lower_, l$col_upper_, col)
          slot <- .hmip_HighsHashTable_ref(fracints, col, c(0, 0))
          pr <- .hmip_HighsHashTable_get(fracints, slot)
          .hmip_HighsHashTable_set(fracints, slot, c(pr[1L] + val, pr[2L] + 1))
        } else {
          if (h$basis$col_status[j] == "basic") next
          glb <- md$domain$col_lower_[j]; gub <- md$domain$col_upper_[j]
          if (min(gub - sol$col_value[j], sol$col_value[j] - glb) <= feastol) next
          ms_am <- ms$model_$a_matrix_
          colStart <- l$a_matrix_$start_[j] + (ms_am$start_[j + 1L] - ms_am$start_[j])
          colEnd <- l$a_matrix_$start_[j + 1L]
          if (colStart == colEnd) next
          for (q in (colStart + 1L):colEnd) {
            row <- l$a_matrix_$index_[q]
            if (lp$lprows$age[row + 1L] == 0L) next
            if (sol$row_value[row + 1L] < l$row_upper_[row + 1L] - feastol) next
            lp$lprows$age[row + 1L] <- 0L
          }
        }
      }
      lp$maxNumFractional <- max(.hmip_HighsHashTable_size(fracints), lp$maxNumFractional)
      if (!is.null(domain) && length(.hmip_HighsDomain_getChangedCols(domain))) {
        solveagain <- TRUE
        next
      }
      keys <- .hmip_HighsHashTable_keys(fracints)
      pv <- .hmip_HighsHashTable_values(fracints)
      lp$fractionalints <- list(first = as.integer(keys),
                                second = vapply(pv, function(p) p[1L] / p[2L], 0))
      fi <- lp$fractionalints
      if (roundable && length(fi$first)) {
        roundsol <- sol$col_value
        for (q in seq_along(fi$first)) {
          col <- fi$first[q] + 1L
          if (md$uplocks[col] == 0 && (cost[col] < 0 || md$downlocks[col] != 0)) {
            roundsol[col] <- min(ceiling(fi$second[q] - feastol),
                                 if (l$col_upper_[col] == Inf) Inf else floor(l$col_upper_[col] + feastol))
          } else {
            roundsol[col] <- max(floor(fi$second[q] + feastol),
                                 if (l$col_lower_[col] == -Inf) -Inf else ceiling(l$col_lower_[col] - feastol))
          }
        }
        cs <- .hmip_HighsCliqueTable_getSubstitutions(md$cliquetable)
        for (k in rev(seq_along(cs$substcol))) {
          rc <- .hmip_cv_col(cs$replace[k])
          roundsol[cs$substcol[k] + 1L] <- if (.hmip_cv_val(cs$replace[k]) == 0L) 1 - roundsol[rc + 1L]
                                           else roundsol[rc + 1L]
        }
        objsum <- c(0, 0)
        for (i in seq_len(ncol)) objsum <- .cd_add(objsum, roundsol[i] * cost[i])
        .hmip_HighsMipSolverData_addIncumbent(md, roundsol, .cd_val(objsum), .hmip_b2_src_SolveLp)
      }
      objsum <- c(0, 0)
      for (i in seq_len(ncol)) objsum <- .cd_add(objsum, sol$col_value[i] * cost[i])
      if (!length(fi$first) && !.hmip_HighsLpRelaxation_unscaledPrimalFeasible(lp$status)) {
        fixSol <- sol$col_value
        for (i in seq_len(ncol)) {
          if (fixSol[i] < l$col_lower_[i]) fixSol[i] <- l$col_lower_[i]
          else if (fixSol[i] > l$col_upper_[i]) fixSol[i] <- l$col_upper_[i]
          else if (ms$model_$integrality_[i] != 0L) fixSol[i] <- .hmip_b2_round(fixSol[i])
        }
        if (.hmip_HighsMipSolverData_checkSolution(md, fixSol)) {
          h$solution$col_value <- fixSol
          lp$status <- if (.hmip_HighsLpRelaxation_unscaledDualFeasible(lp$status)) "Optimal" else "UnscaledPrimalFeasible"
        }
      }
      lp$objective <- .cd_val(objsum)
    } else if (lp$status == "Infeasible") {
      lp$objective <- Inf
    }
    if (!solveagain) break
  }
  lp$status
}

# Playground: putIterate on the first solve, getIterate afterwards, and the
# destructor restores the iterate and resolves
.hmip_HighsLpRelaxation_playground <- function(lp) {
  pg <- new.env(parent = emptyenv())
  pg$lp <- lp; pg$iterateStored <- FALSE
  pg
}
.hmip_Playground_solveLp <- function(pg, localdom) {
  lp <- pg$lp
  if (pg$iterateStored) {
    .hmip_HighsLpRelaxation_flushDomain(lp, localdom)
    .hmip_Highs_getIterate(lp$lpsolver)
  } else {
    .hmip_Highs_putIterate(lp$lpsolver)
    .hmip_HighsLpRelaxation_flushDomain(lp, localdom)
    pg$iterateStored <- TRUE
  }
  .hmip_HighsLpRelaxation_run(lp, FALSE)
}
.hmip_Playground_destroy <- function(pg) {
  if (pg$iterateStored) {
    .hmip_Highs_getIterate(pg$lp$lpsolver)
    .hmip_HighsLpRelaxation_run(pg$lp)
    pg$iterateStored <- FALSE
  }
  invisible()
}

# ---------------------------------------------------------------------------
# HighsSearch (mip/HighsSearch.h/.cpp)

.hmip_b2_branching <- function() list(type = -1L, index = 0L)
.hmip_b2_src_Branching <- 0L
.hmip_b2_src_Heuristic <- 3L
.hmip_b2_src_EvaluateNode <- 10L
.hmip_b2_kUpReliable <- 1L
.hmip_b2_kDownReliable <- 2L

# HighsRandom::bit
.hmip_b2_rng_bit <- function(r) {
  .hs_rng_advance(r)
  .hs_u64_top(.hs_pair_hash(0L, .hs_lo32(r$state), .hs_hi32(r$state)), 1L) == 1
}

# NodeData(parentlb, parentestimate, parentBasis, stabilizerOrbits)
.hmip_b2_nodedata <- function(parentlb = -Inf, parentestimate = -Inf, parentBasis = NULL, stabilizerOrbits = NULL) {
  nd <- new.env(parent = emptyenv())
  nd$lower_bound <- parentlb
  nd$estimate <- parentestimate
  nd$branching_point <- 0
  nd$lp_objective <- -Inf
  nd$other_child_lb <- parentlb
  nd$nodeBasis <- parentBasis
  nd$stabilizerOrbits <- stabilizerOrbits
  nd$branchingdecision <- list(boundval = 0, column = -1L, boundtype = .hmip_b2_lower)
  nd$domgchgStackPos <- -1L
  nd$skipDepthCount <- 0L
  nd$opensubtrees <- 2L
  nd
}
.hmip_b2_top <- function(s) s$nodestack[[length(s$nodestack)]]
# nodes are environments, so `.hmip_b2_top(s)$field <- v` writes into the
# top node in place; the replacement function only hands the stack back
`.hmip_b2_top<-` <- function(s, value) s
.hmip_b2_push <- function(s, nd) s$nodestack[[length(s$nodestack) + 1L]] <- nd
.hmip_b2_pop <- function(s) s$nodestack[[length(s$nodestack)]] <- NULL
.hmip_b2_nchg <- function(dom) length(.hmip_HighsDomain_getChangedCols(dom))
.hmip_b2_stacksize <- function(dom) .hmip_HighsDomain_getNumDomainChanges(dom)
.hmip_b2_chgbound <- function(dom, boundtype, col, val, reason = .hmip_b2_branching())
  .hmip_HighsDomain_changeBound3(dom, boundtype, col, val, reason)

.hmip_new_HighsSearch <- function(mipsolver, pseudocost) {
  s <- new.env(parent = emptyenv())
  s$mipsolver <- mipsolver
  s$lp <- NULL
  s$localdom <- .hmip_HighsDomain_copy(mipsolver$mipdata_$domain)
  s$pseudocost <- pseudocost
  s$random <- .hs_rng(0)
  s$nnodes <- 0
  s$treeweight <- c(0, 0)
  s$depthoffset <- 0L
  s$lpiterations <- 0; s$heurlpiterations <- 0; s$sblpiterations <- 0
  s$upper_limit <- Inf
  s$inheuristic <- FALSE
  s$inbranching <- FALSE
  s$countTreeWeight <- TRUE
  s$childselrule <- if (isTRUE(mipsolver$submip)) "HybridInferenceCost" else "RootSol"
  s$inds <- integer(); s$vals <- numeric()
  s$subrootsol <- numeric()
  s$nodestack <- list()
  s$reliableatnode <- new.env(parent = emptyenv())
  .hmip_HighsDomain_setDomainChangeStack(s$localdom, list(boundval = numeric(), column = integer(),
                                                          boundtype = integer()))
  s
}

.hmip_b2_reliable_flags <- function(s, col) s$reliableatnode[[as.character(col)]] %||% 0L
.hmip_b2_mark_reliable <- function(s, col, flag) {
  k <- as.character(col)
  s$reliableatnode[[k]] <- bitwOr(s$reliableatnode[[k]] %||% 0L, flag)
}

.hmip_HighsSearch_checkSol <- function(s, sol) {
  ms <- s$mipsolver
  objval <- c(0, 0)
  integerfeasible <- TRUE
  cost <- ms$model_$col_cost_; integ <- ms$model_$integrality_
  for (i in seq_len(ms$model_$num_col_)) {
    objval <- .cd_add(objval, sol[i] * cost[i])
    if (!integerfeasible || integ[i] != 1L) next
    if (.hmip_b2_frac(sol[i]) > ms$mipdata_$feastol) integerfeasible <- FALSE
  }
  list(objval = .cd_val(objval), integerfeasible = integerfeasible)
}

.hmip_b2_orbitsValidInChildNode <- function(s, branchChg) {
  cur <- .hmip_b2_top(s)
  so <- cur$stabilizerOrbits
  if (is.null(so) || length(so$orbitCols) == 0L || .hmip_StabilizerOrbits_isStabilized(so, branchChg$column))
    return(TRUE)
  branchChg$boundtype == .hmip_b2_upper && .hmip_HighsDomain_isGlobalBinary(s$localdom, branchChg$column)
}

.hmip_HighsSearch_getCutoffBound <- function(s) min(s$mipsolver$mipdata_$upper_limit, s$upper_limit)
.hmip_HighsSearch_setLpRelaxation <- function(s, lp) s$lp <- lp

.hmip_HighsSearch_setRINSNeighbourhood <- function(s, basesol, relaxsol) {
  ms <- s$mipsolver; dom <- s$localdom
  for (i in seq_len(ms$model_$num_col_)) {
    if (ms$model_$integrality_[i] != 1L) next
    if (dom$col_lower_[i] == dom$col_upper_[i]) next
    intval <- floor(basesol[i] + 0.5)
    if (abs(relaxsol[i] - intval) < ms$mipdata_$feastol) {
      if (dom$col_lower_[i] < intval)
        .hmip_b2_chgbound(dom, .hmip_b2_lower, i - 1L, min(intval, dom$col_upper_[i]), .hmip_b2_unspecified())
      if (dom$col_upper_[i] > intval)
        .hmip_b2_chgbound(dom, .hmip_b2_upper, i - 1L, max(intval, dom$col_lower_[i]), .hmip_b2_unspecified())
    }
  }
  invisible()
}

.hmip_HighsSearch_setRENSNeighbourhood <- function(s, lpsol) {
  ms <- s$mipsolver; dom <- s$localdom; feastol <- ms$mipdata_$feastol
  for (i in seq_len(ms$model_$num_col_)) {
    if (ms$model_$integrality_[i] != 1L) next
    if (dom$col_lower_[i] == dom$col_upper_[i]) next
    downval <- floor(lpsol[i] + feastol)
    upval <- ceiling(lpsol[i] - feastol)
    if (dom$col_lower_[i] < downval) {
      .hmip_b2_chgbound(dom, .hmip_b2_lower, i - 1L, min(downval, dom$col_upper_[i]), .hmip_b2_unspecified())
      if (.hmip_HighsDomain_infeasible(dom)) return(invisible())
    }
    if (dom$col_upper_[i] > upval) {
      .hmip_b2_chgbound(dom, .hmip_b2_upper, i - 1L, max(upval, dom$col_lower_[i]), .hmip_b2_unspecified())
      if (.hmip_HighsDomain_infeasible(dom)) return(invisible())
    }
  }
  invisible()
}

.hmip_HighsSearch_createNewNode <- function(s) {
  nd <- .hmip_b2_nodedata()
  nd$domgchgStackPos <- .hmip_b2_stacksize(s$localdom)
  .hmip_b2_push(s, nd)
}
.hmip_HighsSearch_cutoffNode <- function(s) .hmip_b2_top(s)$opensubtrees <- 0L
.hmip_HighsSearch_setMinReliable <- function(s, minreliable) .hmip_HighsPseudocost_setMinReliable(s$pseudocost, minreliable)
.hmip_HighsSearch_setHeuristic <- function(s, inheuristic) {
  s$inheuristic <- inheuristic
  if (inheuristic) s$childselrule <- "HybridInferenceCost"
}

# branchDownwards / branchUpwards
.hmip_b2_branch_dir <- function(s, col, newbound, branchpoint, boundtype) {
  cur <- .hmip_b2_top(s)
  cur$opensubtrees <- 1L
  cur$branching_point <- branchpoint
  cur$branchingdecision <- list(boundval = newbound, column = col, boundtype = boundtype)
  domchgPos <- .hmip_b2_stacksize(s$localdom)
  pass <- .hmip_b2_orbitsValidInChildNode(s, cur$branchingdecision)
  .hmip_HighsDomain_changeBound(s$localdom, cur$branchingdecision)
  nd <- .hmip_b2_nodedata(cur$lower_bound, cur$estimate, cur$nodeBasis, if (pass) cur$stabilizerOrbits)
  nd$domgchgStackPos <- domchgPos
  .hmip_b2_push(s, nd)
}
.hmip_HighsSearch_branchDownwards <- function(s, col, newub, branchpoint)
  .hmip_b2_branch_dir(s, col, newub, branchpoint, .hmip_b2_upper)
.hmip_HighsSearch_branchUpwards <- function(s, col, newlb, branchpoint)
  .hmip_b2_branch_dir(s, col, newlb, branchpoint, .hmip_b2_lower)

.hmip_HighsSearch_addBoundExceedingConflict <- function(s) {
  md <- s$mipsolver$mipdata_
  if (md$upper_limit != Inf) {
    r <- .hmip_HighsLpRelaxation_computeDualProof(s$lp, md$domain, md$upper_limit)
    if (r$ok) {
      s$inds <- r$inds; s$vals <- r$vals
      if (.hmip_HighsDomain_infeasible(md$domain)) return(invisible())
      .hmip_HighsDomain_conflictAnalysisProof(s$localdom, r$inds, r$vals, length(r$inds), r$rhs, md$conflictPool)
      cg <- .hmip_new_HighsCutGeneration(s$lp, md$cutpool)
      .hmip_HighsCutGeneration_generateConflict(cg, s$localdom, r$inds, r$vals, r$rhs)
    }
  }
  invisible()
}

.hmip_HighsSearch_addInfeasibleConflict <- function(s) {
  md <- s$mipsolver$mipdata_
  if (.hmip_Highs_getModelStatus(s$lp$lpsolver) == "ObjectiveBound") .hmip_HighsLpRelaxation_performAging(s$lp)
  r <- .hmip_HighsLpRelaxation_computeDualInfProof(s$lp, md$domain)
  if (r$ok) {
    s$inds <- r$inds; s$vals <- r$vals
    if (.hmip_HighsDomain_infeasible(md$domain)) return(invisible())
    .hmip_HighsDomain_conflictAnalysisProof(s$localdom, r$inds, r$vals, length(r$inds), r$rhs, md$conflictPool)
    cg <- .hmip_new_HighsCutGeneration(s$lp, md$cutpool)
    .hmip_HighsCutGeneration_generateConflict(cg, s$localdom, r$inds, r$vals, r$rhs)
  }
  invisible()
}

# selectBranchingCandidate: returns list(candidate, downNodeLb, upNodeLb);
# candidate is the 0-based position in the fractional integers or -1
.hmip_HighsSearch_selectBranchingCandidate <- function(s, maxSbIters) {
  ms <- s$mipsolver; md <- ms$mipdata_
  lp <- s$lp; dom <- s$localdom; pc <- s$pseudocost
  feastol <- md$feastol
  fi <- lp$fractionalints
  fcol <- fi$first; fval <- fi$second
  numfrac <- length(fcol)
  curlb <- .hmip_HighsSearch_getCurrentLowerBound(s)
  upscore <- rep(Inf, numfrac); downscore <- rep(Inf, numfrac)
  upbound <- rep(curlb, numfrac); downbound <- rep(curlb, numfrac)
  upscorereliable <- logical(numfrac); downscorereliable <- logical(numfrac)
  for (k in seq_len(numfrac)) {
    col <- fcol[k]; fracval <- fval[k]
    if (.hmip_HighsPseudocost_isReliable(pc, col)) {
      upscore[k] <- .hmip_HighsPseudocost_getPseudocostUp(pc, col, fracval)
      downscore[k] <- .hmip_HighsPseudocost_getPseudocostDown(pc, col, fracval)
      upscorereliable[k] <- TRUE; downscorereliable[k] <- TRUE
    } else {
      flags <- .hmip_b2_reliable_flags(s, col)
      if (bitwAnd(flags, .hmip_b2_kUpReliable)) {
        upscore[k] <- .hmip_HighsPseudocost_getPseudocostUp(pc, col, fracval)
        upscorereliable[k] <- TRUE
      }
      if (bitwAnd(flags, .hmip_b2_kDownReliable)) {
        downscore[k] <- .hmip_HighsPseudocost_getPseudocostDown(pc, col, fracval)
        downscorereliable[k] <- TRUE
      }
    }
  }
  minScore <- feastol
  selectBestScore <- function(finalSelection) {
    best <- -1L; bestscore <- -1; bestnodes <- -1; bestnumnodes <- 0
    oldminscore <- minScore
    for (k in seq_len(numfrac)) {
      if (upscore[k] <= oldminscore) upscorereliable[k] <<- TRUE
      if (downscore[k] <= oldminscore) downscorereliable[k] <<- TRUE
      sc <- 1e-3 * min(if (upscorereliable[k]) upscore[k] else 0, if (downscorereliable[k]) downscore[k] else 0)
      minScore <<- max(sc, minScore)
      if (upscore[k] <= oldminscore || downscore[k] <= oldminscore) {
        score <- .hmip_HighsPseudocost_getScore(pc, fcol[k], min(upscore[k], oldminscore), min(downscore[k], oldminscore))
      } else if (upscore[k] == Inf || downscore[k] == Inf) {
        score <- if (finalSelection) .hmip_HighsPseudocost_getScore2(pc, fcol[k], fval[k]) else Inf
      } else {
        score <- .hmip_HighsPseudocost_getScore(pc, fcol[k], upscore[k], downscore[k])
      }
      upnodes <- .hmip_HighsNodeQueue_numNodesUp(md$nodequeue, fcol[k])
      downnodes <- .hmip_HighsNodeQueue_numNodesDown(md$nodequeue, fcol[k])
      nodes <- 0
      numnodes <- upnodes + downnodes
      if (upnodes != 0 || downnodes != 0) nodes <- (downnodes / numnodes) * (upnodes / numnodes)
      if (score > bestscore ||
          (score > bestscore - feastol && (nodes > bestnodes || (nodes == bestnodes && numnodes > bestnumnodes)))) {
        bestscore <- score; best <- k; bestnodes <- nodes; bestnumnodes <- numnodes
      }
    }
    best
  }

  pg <- .hmip_HighsLpRelaxation_playground(lp)
  on.exit(.hmip_Playground_destroy(pg), add = TRUE)

  repeat {
    mustStop <- .hmip_HighsSearch_getStrongBranchingLpIterations(s) >= maxSbIters ||
      .hmip_HighsMipSolverData_checkLimits(md)
    candidate <- selectBestScore(mustStop)
    if ((upscorereliable[candidate] && downscorereliable[candidate]) || mustStop)
      return(list(candidate = candidate - 1L, downNodeLb = downbound[candidate], upNodeLb = upbound[candidate]))
    .hmip_HighsLpRelaxation_setObjectiveLimit(lp, md$upper_limit)
    col <- fcol[candidate]; fracval <- fval[candidate]
    upval <- ceiling(fracval); downval <- floor(fracval)

    # true when the probed solution survives the propagation of fixing the
    # other candidate k (or no fixing was needed)
    probe <- function(k, boundtype, bound, sol, numChangedCols, domchgStackSize) {
      .hmip_b2_chgbound(dom, boundtype, fcol[k], bound)
      if (!.hmip_HighsDomain_infeasible(dom)) .hmip_HighsDomain_propagate(dom)
      if (.hmip_HighsDomain_infeasible(dom)) {
        .hmip_HighsDomain_conflictAnalysis(dom, md$conflictPool)
        .hmip_HighsDomain_backtrack(dom)
        .hmip_HighsDomain_clearChangedCols(dom, numChangedCols)
        return(FALSE)
      }
      newStackSize <- .hmip_b2_stacksize(dom)
      solutionValid <- TRUE
      j <- domchgStackSize + 1L
      while (j < newStackSize) {
        c1 <- dom$domchgstack_column[j + 1L] + 1L; bv <- dom$domchgstack_boundval[j + 1L]
        if (dom$domchgstack_boundtype[j + 1L] == .hmip_b2_lower) {
          if (bv > sol[c1] + feastol) { solutionValid <- FALSE; break }
        } else {
          if (bv < sol[c1] - feastol) { solutionValid <- FALSE; break }
        }
        j <- j + 1L
      }
      .hmip_HighsDomain_backtrack(dom)
      .hmip_HighsDomain_clearChangedCols(dom, numChangedCols)
      solutionValid
    }

    analyzeSolution <- function(objdelta, sol) {
      numChangedCols <- .hmip_b2_nchg(dom)
      domchgStackSize <- .hmip_b2_stacksize(dom)
      for (k in seq_len(numfrac)) {
        if (fcol[k] == col) next
        c1 <- fcol[k] + 1L
        otherfracval <- fval[k]
        otherdownval <- floor(otherfracval); otherupval <- ceiling(otherfracval)
        if (sol[c1] <= otherdownval + feastol) {
          if (dom$col_upper_[c1] > otherdownval + feastol &&
              !probe(k, .hmip_b2_upper, otherdownval, sol, numChangedCols, domchgStackSize)) next
          if (objdelta <= feastol) {
            .hmip_HighsPseudocost_addObservation(pc, fcol[k], otherdownval - otherfracval, objdelta)
            .hmip_b2_mark_reliable(s, fcol[k], .hmip_b2_kDownReliable)
          }
          downscore[k] <<- min(downscore[k], objdelta)
        } else if (sol[c1] >= otherupval - feastol) {
          if (dom$col_lower_[c1] < otherupval - feastol &&
              !probe(k, .hmip_b2_lower, otherupval, sol, numChangedCols, domchgStackSize)) next
          if (objdelta <= feastol) {
            .hmip_HighsPseudocost_addObservation(pc, fcol[k], otherupval - otherfracval, objdelta)
            .hmip_b2_mark_reliable(s, fcol[k], .hmip_b2_kUpReliable)
          }
          upscore[k] <<- min(upscore[k], objdelta)
        }
      }
    }

    # the other child after the probed one is cut off
    cutoffChild <- function(upbranch, opensub, other_lb = NULL) {
      if (upbranch) .hmip_HighsSearch_branchDownwards(s, col, downval, fracval)
      else .hmip_HighsSearch_branchUpwards(s, col, upval, fracval)
      parent <- s$nodestack[[length(s$nodestack) - 1L]]
      parent$opensubtrees <- opensub
      if (!is.null(other_lb)) parent$other_child_lb <- other_lb
      parent$skipDepthCount <- 1L
      s$depthoffset <- s$depthoffset - 1L
      TRUE
    }

    strongBranch <- function(upbranch) {
      inferences <- -.hmip_b2_stacksize(dom) - 1
      domchg <- list(boundval = if (upbranch) upval else downval, column = col,
                     boundtype = if (upbranch) .hmip_b2_lower else .hmip_b2_upper)
      top <- .hmip_b2_top(s)
      orbitalFixing <- !is.null(top$stabilizerOrbits) && .hmip_b2_orbitsValidInChildNode(s, domchg)
      .hmip_HighsDomain_changeBound(dom, domchg)
      .hmip_HighsDomain_propagate(dom)
      if (!.hmip_HighsDomain_infeasible(dom)) {
        if (orbitalFixing) .hmip_StabilizerOrbits_orbitalFixing(top$stabilizerOrbits, dom)
        else .hmip_HighsSymmetries_propagateOrbitopes(md$symmetries, dom)
      }
      inferences <- inferences + .hmip_b2_stacksize(dom)
      if (.hmip_HighsDomain_infeasible(dom)) {
        .hmip_HighsDomain_conflictAnalysis(dom, md$conflictPool)
        .hmip_HighsPseudocost_addCutoffObservation(pc, col, upbranch)
        .hmip_HighsDomain_backtrack(dom)
        .hmip_HighsDomain_clearChangedCols(dom)
        return(cutoffChild(upbranch, 0L))
      }
      .hmip_HighsPseudocost_addInferenceObservation(pc, col, inferences, upbranch)
      numiters <- lp$numlpiters
      status <- .hmip_Playground_solveLp(pg, dom)
      numiters <- lp$numlpiters - numiters
      s$lpiterations <- s$lpiterations + numiters
      s$sblpiterations <- s$sblpiterations + numiters
      if (.hmip_HighsLpRelaxation_scaledOptimal(status)) {
        .hmip_HighsLpRelaxation_performAging(lp)
        delta <- if (upbranch) upval - fracval else downval - fracval
        sol <- lp$lpsolver$solution$col_value
        cs <- .hmip_HighsSearch_checkSol(s, sol)
        solobj <- cs$objval
        objdelta <- max(solobj - lp$objective, 0)
        if (objdelta <= md$epsilon) objdelta <- 0
        if (upbranch) {
          upscore[candidate] <<- objdelta; upscorereliable[candidate] <<- TRUE
          .hmip_b2_mark_reliable(s, col, .hmip_b2_kUpReliable)
        } else {
          downscore[candidate] <<- objdelta; downscorereliable[candidate] <<- TRUE
          .hmip_b2_mark_reliable(s, col, .hmip_b2_kDownReliable)
        }
        .hmip_HighsPseudocost_addObservation(pc, col, delta, objdelta)
        analyzeSolution(objdelta, sol)
        if (.hmip_HighsLpRelaxation_unscaledPrimalFeasible(status) && cs$integerfeasible) {
          cutoffbnd <- .hmip_HighsSearch_getCutoffBound(s)
          .hmip_HighsMipSolverData_addIncumbent(md, lp$lpsolver$solution$col_value, solobj,
                                                if (s$inheuristic) .hmip_b2_src_Heuristic else .hmip_b2_src_Branching)
          if (md$upper_limit < cutoffbnd) .hmip_HighsLpRelaxation_setObjectiveLimit(lp, md$upper_limit)
        }
        if (.hmip_HighsLpRelaxation_unscaledDualFeasible(status)) {
          if (upbranch) upbound[candidate] <<- solobj else downbound[candidate] <<- solobj
          if (solobj > md$optimality_limit) {
            .hmip_HighsSearch_addBoundExceedingConflict(s)
            pruned <- solobj > .hmip_HighsSearch_getCutoffBound(s)
            .hmip_HighsDomain_backtrack(dom)
            .hmip_HighsLpRelaxation_flushDomain(lp, dom)
            return(cutoffChild(upbranch, if (pruned) 0L else 1L, solobj))
          }
        } else if (solobj > .hmip_HighsSearch_getCutoffBound(s)) {
          .hmip_HighsSearch_addBoundExceedingConflict(s)
          .hmip_HighsDomain_propagate(dom)
          if (.hmip_HighsDomain_infeasible(dom)) {
            .hmip_HighsDomain_backtrack(dom)
            .hmip_HighsLpRelaxation_flushDomain(lp, dom)
            return(cutoffChild(upbranch, 0L))
          }
        }
      } else if (status == "Infeasible") {
        .hmip_HighsSearch_addInfeasibleConflict(s)
        .hmip_HighsPseudocost_addCutoffObservation(pc, col, upbranch)
        .hmip_HighsDomain_backtrack(dom)
        .hmip_HighsLpRelaxation_flushDomain(lp, dom)
        return(cutoffChild(upbranch, 0L))
      } else {
        downscore[candidate] <<- 0; upscore[candidate] <<- 0
        downscorereliable[candidate] <<- TRUE; upscorereliable[candidate] <<- TRUE
        .hmip_b2_mark_reliable(s, col, .hmip_b2_kUpReliable)
        .hmip_b2_mark_reliable(s, col, .hmip_b2_kDownReliable)
      }
      .hmip_HighsDomain_backtrack(dom)
      .hmip_HighsLpRelaxation_flushDomain(lp, dom)
      FALSE
    }

    dinf <- .hmip_HighsPseudocost_getAvgInferencesDown(pc, col)
    uinf <- .hmip_HighsPseudocost_getAvgInferencesUp(pc, col)
    evalDown <- !downscorereliable[candidate] &&
      (upscorereliable[candidate] || downscore[candidate] > upscore[candidate] ||
         (downscore[candidate] == upscore[candidate] && dinf >= uinf))
    if (strongBranch(!evalDown))
      return(list(candidate = -1L, downNodeLb = downbound[candidate], upNodeLb = upbound[candidate]))
  }
}

.hmip_HighsSearch_getParentNodeData <- function(s) {
  n <- length(s$nodestack)
  if (n <= 1L) NULL else s$nodestack[[n - 1L]]
}

# push the current node (top of the stack) to the queue or prune it
.hmip_b2_node_to_queue <- function(s, nodequeue) {
  dom <- s$localdom; md <- s$mipsolver$mipdata_
  top <- .hmip_b2_top(s)
  oldchangedcols <- .hmip_b2_nchg(dom)
  prune <- top$lower_bound > .hmip_HighsSearch_getCutoffBound(s)
  if (!prune) {
    .hmip_HighsDomain_propagate(dom)
    .hmip_HighsDomain_clearChangedCols(dom, oldchangedcols)
    prune <- .hmip_HighsDomain_infeasible(dom)
    if (prune) .hmip_HighsDomain_conflictAnalysis(dom, md$conflictPool)
  }
  if (!prune) {
    red <- .hmip_HighsDomain_getReducedDomainChangeStack(dom)
    tw <- .hmip_HighsNodeQueue_emplaceNode(nodequeue, red$stack, red$branchingPositions,
                                           max(top$lower_bound, .hmip_HighsDomain_getObjectiveLowerBound(dom)),
                                           top$estimate, .hmip_HighsSearch_getCurrentDepth(s))
    if (s$countTreeWeight) s$treeweight <- .cd_add(s$treeweight, tw)
  } else if (s$countTreeWeight) {
    s$treeweight <- .cd_add(s$treeweight, 2^(1 - .hmip_HighsSearch_getCurrentDepth(s)))
  }
  top$opensubtrees <- 0L
}

.hmip_HighsSearch_currentNodeToQueue <- function(s, nodequeue) .hmip_b2_node_to_queue(s, nodequeue)

.hmip_HighsSearch_openNodesToQueue <- function(s, nodequeue) {
  if (!length(s$nodestack)) return(invisible())
  basis <- NULL
  for (nd in s$nodestack) if (!is.null(nd$nodeBasis)) {
    basis <- nd$nodeBasis
    nd$nodeBasis <- NULL
    break
  }
  if (.hmip_b2_top(s)$opensubtrees == 0L) .hmip_HighsSearch_backtrack(s, FALSE)
  while (length(s$nodestack)) {
    .hmip_b2_node_to_queue(s, nodequeue)
    .hmip_HighsSearch_backtrack(s, FALSE)
  }
  .hmip_HighsLpRelaxation_flushDomain(s$lp, s$localdom)
  if (!is.null(basis)) {
    if (length(basis$row_status) == .hmip_HighsLpRelaxation_numRows(s$lp))
      .hmip_HighsLpRelaxation_setStoredBasis(s$lp, basis)
    .hmip_HighsLpRelaxation_recoverBasis(s$lp)
  }
  invisible()
}

.hmip_HighsSearch_flushStatistics <- function(s) {
  md <- s$mipsolver$mipdata_
  md$num_nodes <- md$num_nodes + s$nnodes; s$nnodes <- 0
  md$pruned_treeweight <- .cd_iadd(md$pruned_treeweight, s$treeweight); s$treeweight <- c(0, 0)
  md$total_lp_iterations <- md$total_lp_iterations + s$lpiterations; s$lpiterations <- 0
  md$heuristic_lp_iterations <- md$heuristic_lp_iterations + s$heurlpiterations; s$heurlpiterations <- 0
  md$sb_lp_iterations <- md$sb_lp_iterations + s$sblpiterations; s$sblpiterations <- 0
  invisible()
}
.hmip_HighsSearch_getHeuristicLpIterations <- function(s) s$heurlpiterations + s$mipsolver$mipdata_$heuristic_lp_iterations
.hmip_HighsSearch_getTotalLpIterations <- function(s) s$lpiterations + s$mipsolver$mipdata_$total_lp_iterations
.hmip_HighsSearch_getLocalLpIterations <- function(s) s$lpiterations
.hmip_HighsSearch_getLocalNodes <- function(s) s$nnodes
.hmip_HighsSearch_getStrongBranchingLpIterations <- function(s) s$sblpiterations + s$mipsolver$mipdata_$sb_lp_iterations
.hmip_HighsSearch_hasNode <- function(s) length(s$nodestack) > 0L
.hmip_HighsSearch_currentNodePruned <- function(s) .hmip_b2_top(s)$opensubtrees == 0L
.hmip_HighsSearch_getCurrentEstimate <- function(s) .hmip_b2_top(s)$estimate
.hmip_HighsSearch_getCurrentLowerBound <- function(s) .hmip_b2_top(s)$lower_bound
.hmip_HighsSearch_getCurrentDepth <- function(s) length(s$nodestack) + s$depthoffset
.hmip_HighsSearch_getLocalDomain <- function(s) s$localdom
.hmip_HighsSearch_getPseudoCost <- function(s) s$pseudocost

.hmip_HighsSearch_resetLocalDomain <- function(s) {
  .hmip_HighsLpRelaxation_resetToGlobalDomain(s$lp)
  .hmip_HighsDomain_assign(s$localdom, s$mipsolver$mipdata_$domain)
  invisible()
}

.hmip_HighsSearch_installNode <- function(s, node) {
  md <- s$mipsolver$mipdata_
  dom <- s$localdom
  .hmip_HighsDomain_setDomainChangeStack2(dom, node$domchgstack, node$branchings)
  globalSymmetriesValid <- TRUE
  if (!is.null(md$globalOrbits)) {
    for (i in .hmip_HighsDomain_getBranchingPositions(dom)) {
      col <- dom$domchgstack_column[i + 1L]
      if (md$symmetries$columnPosition[col + 1L] == -1L) next
      if (!.hmip_HighsDomain_isBinary(md$domain, col) ||
          (dom$domchgstack_boundtype[i + 1L] == .hmip_b2_lower && dom$domchgstack_boundval[i + 1L] == 1)) {
        globalSymmetriesValid <- FALSE
        break
      }
    }
  }
  .hmip_b2_push(s, .hmip_b2_nodedata(node$lower_bound, node$estimate, NULL,
                                     if (globalSymmetriesValid) md$globalOrbits))
  s$subrootsol <- numeric()
  s$depthoffset <- node$depth - 1L
  invisible()
}

# the pseudocost cutoff observation for the branching that led here
.hmip_b2_parent_cutoff <- function(s, parent) {
  if (!is.null(parent) && parent$lp_objective != -Inf &&
      parent$branching_point != parent$branchingdecision$boundval)
    .hmip_HighsPseudocost_addCutoffObservation(s$pseudocost, parent$branchingdecision$column,
                                               parent$branchingdecision$boundtype == .hmip_b2_lower)
}
.hmip_b2_domain_infeasible <- function(s, parent) {
  .hmip_HighsDomain_clearChangedCols(s$localdom)
  .hmip_b2_parent_cutoff(s, parent)
  .hmip_HighsDomain_conflictAnalysis(s$localdom, s$mipsolver$mipdata_$conflictPool)
  "DomainInfeasible"
}

.hmip_HighsSearch_evaluateNode <- function(s) {
  ms <- s$mipsolver; md <- ms$mipdata_
  lp <- s$lp; dom <- s$localdom; pc <- s$pseudocost
  currnode <- .hmip_b2_top(s)
  parent <- .hmip_HighsSearch_getParentNodeData(s)
  if (!s$inheuristic && currnode$lower_bound > md$optimality_limit) return("SubOptimal")
  .hmip_HighsDomain_propagate(dom)
  if (!s$inheuristic && !.hmip_HighsDomain_infeasible(dom)) {
    if (md$symmetries$numPerms > 0 && is.null(currnode$stabilizerOrbits) &&
        (is.null(parent) || is.null(parent$stabilizerOrbits) || length(parent$stabilizerOrbits$orbitCols) > 0L))
      currnode$stabilizerOrbits <- .hmip_HighsSymmetries_computeStabilizerOrbits(md$symmetries, dom)
    if (!is.null(currnode$stabilizerOrbits)) .hmip_StabilizerOrbits_orbitalFixing(currnode$stabilizerOrbits, dom)
    else .hmip_HighsSymmetries_propagateOrbitopes(md$symmetries, dom)
  }
  if (!is.null(parent)) {
    inferences <- .hmip_b2_stacksize(dom) - (currnode$domgchgStackPos + 1L)
    .hmip_HighsPseudocost_addInferenceObservation(pc, parent$branchingdecision$column, inferences,
                                                  parent$branchingdecision$boundtype == .hmip_b2_lower)
  }
  result <- "Open"
  if (.hmip_HighsDomain_infeasible(dom)) {
    result <- .hmip_b2_domain_infeasible(s, parent)
  } else {
    .hmip_HighsLpRelaxation_flushDomain(lp, dom)
    .hmip_HighsLpRelaxation_setObjectiveLimit(lp, md$upper_limit)
    oldnumiters <- lp$numlpiters
    status <- .hmip_HighsLpRelaxation_resolveLp(lp, dom)
    s$lpiterations <- s$lpiterations + lp$numlpiters - oldnumiters
    currnode$lower_bound <- max(.hmip_HighsDomain_getObjectiveLowerBound(dom), currnode$lower_bound)
    if (.hmip_HighsDomain_infeasible(dom)) {
      result <- .hmip_b2_domain_infeasible(s, parent)
    } else if (.hmip_HighsLpRelaxation_scaledOptimal(status)) {
      .hmip_HighsLpRelaxation_storeBasis(lp)
      .hmip_HighsLpRelaxation_performAging(lp)
      currnode$nodeBasis <- lp$basischeckpoint
      currnode$estimate <- .hmip_HighsLpRelaxation_computeBestEstimate(lp, pc)
      currnode$lp_objective <- lp$objective
      if (!is.null(parent) && parent$lp_objective != -Inf &&
          parent$branching_point != parent$branchingdecision$boundval) {
        delta <- parent$branchingdecision$boundval - parent$branching_point
        objdelta <- max(0, currnode$lp_objective - parent$lp_objective)
        .hmip_HighsPseudocost_addObservation(pc, parent$branchingdecision$column, delta, objdelta)
      }
      if (.hmip_HighsLpRelaxation_unscaledPrimalFeasible(status) && length(lp$fractionalints$first) == 0L) {
        cutoffbnd <- .hmip_HighsSearch_getCutoffBound(s)
        .hmip_HighsMipSolverData_addIncumbent(md, lp$lpsolver$solution$col_value, lp$objective,
                                              if (s$inheuristic) .hmip_b2_src_Heuristic else .hmip_b2_src_EvaluateNode)
        if (md$upper_limit < cutoffbnd) .hmip_HighsLpRelaxation_setObjectiveLimit(lp, md$upper_limit)
        if (.hmip_HighsLpRelaxation_unscaledDualFeasible(status)) {
          .hmip_HighsSearch_addBoundExceedingConflict(s)
          result <- "BoundExceeding"
        }
      }
      if (result == "Open") {
        if (.hmip_HighsLpRelaxation_unscaledDualFeasible(status)) {
          currnode$lower_bound <- max(currnode$lp_objective, currnode$lower_bound)
          if (currnode$lower_bound > .hmip_HighsSearch_getCutoffBound(s)) {
            result <- "BoundExceeding"
            .hmip_HighsSearch_addBoundExceedingConflict(s)
          } else if (md$upper_limit != Inf) {
            if (!s$inheuristic) {
              gap <- md$upper_limit - lp$objective
              .hmip_HighsLpRelaxation_computeBasicDegenerateDuals(lp, gap + max(10 * md$feastol, md$epsilon * gap), dom)
            }
            .hmip_HighsRedcostFixing_propagateRedCost(ms, dom, lp)
            .hmip_HighsDomain_propagate(dom)
            if (.hmip_HighsDomain_infeasible(dom)) result <- .hmip_b2_domain_infeasible(s, parent)
            else if (.hmip_b2_nchg(dom)) return(.hmip_HighsSearch_evaluateNode(s))
          } else if (!s$inheuristic) {
            .hmip_HighsLpRelaxation_computeBasicDegenerateDuals(lp, Inf, dom)
            .hmip_HighsDomain_propagate(dom)
            if (.hmip_HighsDomain_infeasible(dom)) result <- .hmip_b2_domain_infeasible(s, parent)
            else if (.hmip_b2_nchg(dom)) return(.hmip_HighsSearch_evaluateNode(s))
          }
        } else if (lp$objective > .hmip_HighsSearch_getCutoffBound(s)) {
          .hmip_HighsSearch_addBoundExceedingConflict(s)
          .hmip_HighsDomain_propagate(dom)
          if (.hmip_HighsDomain_infeasible(dom)) result <- "BoundExceeding"
        }
      }
    } else if (status == "Infeasible") {
      result <- if (.hmip_Highs_getModelStatus(lp$lpsolver) == "ObjectiveBound") "BoundExceeding" else "LpInfeasible"
      .hmip_HighsSearch_addInfeasibleConflict(s)
      .hmip_b2_parent_cutoff(s, parent)
    }
  }
  if (result != "Open") {
    s$treeweight <- .cd_add(s$treeweight, 2^(1 - .hmip_HighsSearch_getCurrentDepth(s)))
    currnode$opensubtrees <- 0L
  } else if (!s$inheuristic) {
    if (currnode$lower_bound > md$optimality_limit) {
      result <- "SubOptimal"
      .hmip_HighsSearch_addBoundExceedingConflict(s)
    }
  }
  result
}

# set the branching decision of currnode to the up (TRUE) or down child
.hmip_b2_choose <- function(currnode, up, downNodeLb, upNodeLb) {
  bp <- currnode$branching_point
  if (up) {
    currnode$branchingdecision$boundtype <- .hmip_b2_lower
    currnode$branchingdecision$boundval <- ceiling(bp)
    currnode$other_child_lb <- downNodeLb
    upNodeLb
  } else {
    currnode$branchingdecision$boundtype <- .hmip_b2_upper
    currnode$branchingdecision$boundval <- floor(bp)
    currnode$other_child_lb <- upNodeLb
    downNodeLb
  }
}

.hmip_HighsSearch_branch <- function(s) {
  ms <- s$mipsolver; md <- ms$mipdata_
  dom <- s$localdom; pc <- s$pseudocost
  .hmip_b2_top(s)$branchingdecision$column <- -1L
  s$inbranching <- TRUE
  minrel <- .hmip_HighsPseudocost_getMinReliable(pc)
  childLb <- .hmip_HighsSearch_getCurrentLowerBound(s)
  result <- "Open"
  while (.hmip_b2_top(s)$opensubtrees == 2L && .hmip_HighsLpRelaxation_scaledOptimal(s$lp$status) &&
         length(s$lp$fractionalints$first) > 0L) {
    sbmaxiters <- 0
    if (minrel > 0) {
      sbiters <- .hmip_HighsSearch_getStrongBranchingLpIterations(s)
      sbmaxiters <- 100000 + (.hmip_HighsSearch_getTotalLpIterations(s) - .hmip_HighsSearch_getHeuristicLpIterations(s) -
                                .hmip_HighsSearch_getStrongBranchingLpIterations(s)) %/% 2
      if (sbiters > sbmaxiters) {
        .hmip_HighsPseudocost_setMinReliable(pc, 0L)
      } else if (sbiters > sbmaxiters %/% 2) {
        reductionratio <- (sbiters - sbmaxiters %/% 2) / (sbmaxiters - sbmaxiters %/% 2)
        minrelreduced <- trunc(minrel - reductionratio * (minrel - 1))
        .hmip_HighsPseudocost_setMinReliable(pc, min(minrel, minrelreduced))
      }
    }
    degeneracyFac <- .hmip_HighsLpRelaxation_computeLPDegneracy(s$lp, dom)
    .hmip_HighsPseudocost_setDegeneracyFactor(pc, degeneracyFac)
    if (degeneracyFac >= 10) .hmip_HighsPseudocost_setMinReliable(pc, 0L)
    sel <- .hmip_HighsSearch_selectBranchingCandidate(s, sbmaxiters)
    downNodeLb <- sel$downNodeLb; upNodeLb <- sel$upNodeLb
    currnode <- .hmip_b2_top(s)
    childLb <- currnode$lower_bound
    if (sel$candidate != -1L) {
      col <- s$lp$fractionalints$first[sel$candidate + 1L]
      currnode$branchingdecision$column <- col
      currnode$branching_point <- s$lp$fractionalints$second[sel$candidate + 1L]
      bp <- currnode$branching_point
      up <- switch(s$childselrule,
        Up = TRUE,
        Down = FALSE,
        RootSol = {
          downPrio <- .hmip_HighsPseudocost_getAvgInferencesDown(pc, col) + md$epsilon
          upPrio <- .hmip_HighsPseudocost_getAvgInferencesUp(pc, col) + md$epsilon
          downVal <- floor(bp); upVal <- ceiling(bp)
          rootsol <- NULL
          if (length(s$subrootsol)) {
            rootsol <- s$subrootsol[col + 1L]
          } else {
            if (currnode$lp_objective != -Inf) s$subrootsol <- s$lp$lpsolver$solution$col_value
            if (length(md$rootlpsol)) rootsol <- md$rootlpsol[col + 1L]
          }
          if (!is.null(rootsol)) {
            if (rootsol < downVal) rootsol <- downVal else if (rootsol > upVal) rootsol <- upVal
            upPrio <- upPrio * (1 + (bp - rootsol))
            downPrio <- downPrio * (1 + (rootsol - bp))
          }
          upPrio + md$epsilon >= downPrio
        },
        Obj = ms$model_$col_cost_[col + 1L] >= 0,
        Random = .hmip_b2_rng_bit(s$random),
        BestCost = !(.hmip_HighsPseudocost_getPseudocostUp(pc, col, bp, md$feastol) >
                       .hmip_HighsPseudocost_getPseudocostDown(pc, col, bp, md$feastol)),
        WorstCost = .hmip_HighsPseudocost_getPseudocostUp(pc, col, bp) >=
          .hmip_HighsPseudocost_getPseudocostDown(pc, col, bp),
        Disjunction = {
          nu <- .hmip_HighsNodeQueue_numNodesUp(md$nodequeue, col)
          nd <- .hmip_HighsNodeQueue_numNodesDown(md$nodequeue, col)
          if (nu > nd) TRUE else if (nd > nu) FALSE else ms$model_$col_cost_[col + 1L] >= 0
        },
        HybridInferenceCost = {
          upScore <- (1 + .hmip_HighsPseudocost_getAvgInferencesUp(pc, col)) /
            .hmip_HighsPseudocost_getPseudocostUp(pc, col, bp, md$feastol)
          downScore <- (1 + .hmip_HighsPseudocost_getAvgInferencesDown(pc, col)) /
            .hmip_HighsPseudocost_getPseudocostDown(pc, col, bp, md$feastol)
          upScore >= downScore
        })
      childLb <- .hmip_b2_choose(currnode, up, downNodeLb, upNodeLb)
      result <- "Branched"
      break
    }
    result <- .hmip_HighsSearch_evaluateNode(s)
    if (result == "SubOptimal") break
  }
  s$inbranching <- FALSE
  currnode <- .hmip_b2_top(s)
  .hmip_HighsPseudocost_setMinReliable(pc, minrel)
  .hmip_HighsPseudocost_setDegeneracyFactor(pc, 1)
  if (currnode$opensubtrees != 2L || result == "SubOptimal") return(result)

  if (currnode$branchingdecision$column == -1L) {
    bestscore <- -1
    .hmip_HighsPseudocost_setDegeneracyFactor(pc, 1e6)
    for (i in md$integral_cols) {
      j <- i + 1L
      lo <- dom$col_lower_[j]; hi <- dom$col_upper_[j]
      if (hi - lo < 0.5) next
      fracval <- if (lo != -Inf) lo + 0.5 else if (hi != Inf) hi - 0.5 else 0.5
      score <- .hmip_HighsPseudocost_getScore2(pc, i, fracval)
      if (score > bestscore) {
        bestscore <- score
        cost <- if (.hmip_HighsLpRelaxation_unscaledDualFeasible(s$lp$status)) s$lp$lpsolver$solution$col_dual[j]
                else ms$model_$col_cost_[j]
        iu <- .hmip_HighsPseudocost_getAvgInferencesUp(pc, i)
        id <- .hmip_HighsPseudocost_getAvgInferencesDown(pc, i)
        if (abs(cost) > md$feastol && .hmip_HighsSearch_getCutoffBound(s) < Inf) {
          branchUp <- cost > 0
        } else if (iu > id + md$feastol) {
          branchUp <- TRUE
        } else if (iu < id - md$feastol) {
          branchUp <- FALSE
        } else {
          nchg <- .hmip_HighsDomain_getNumDomainChanges(dom)
          lpos <- .hmip_HighsDomain_getColLowerPos(dom, i, nchg)$pos
          upos <- .hmip_HighsDomain_getColUpperPos(dom, i, nchg)$pos
          branchUp <- lpos <= upos
        }
        if (branchUp) {
          v <- ceiling(fracval)
          currnode$branching_point <- v
          currnode$branchingdecision <- list(boundval = v, column = i, boundtype = .hmip_b2_lower)
        } else {
          v <- floor(fracval)
          currnode$branching_point <- v
          currnode$branchingdecision <- list(boundval = v, column = i, boundtype = .hmip_b2_upper)
        }
      }
    }
    .hmip_HighsPseudocost_setDegeneracyFactor(pc, 1)
  }

  if (currnode$branchingdecision$column == -1L) {
    if (s$lp$status == "Optimal") {
      currnode$opensubtrees <- 0L
      return("LpInfeasible")
    }
    .hmip_HighsLpRelaxation_setIterationLimit(s$lp)
    lpCopy <- .hmip_new_HighsLpRelaxation(ms)
    .hmip_HighsLpRelaxation_loadModel(lpCopy)
    .hmip_Highs_changeColsBoundsRange(lpCopy$lpsolver, 0L, ms$model_$num_col_ - 1L, dom$col_lower_, dom$col_upper_)
    tmpLp <- s$lp
    s$lp <- lpCopy
    on.exit(s$lp <- tmpLp, add = TRUE)
    h <- lpCopy$lpsolver
    .hmip_Highs_setOptionValue(h, "presolve", "on")
    result <- .hmip_HighsSearch_evaluateNode(s)
    if (result == "Open") {
      .hmip_Highs_clearSolver(h)
      .hmip_Highs_setOptionValue(h, "simplex_strategy", 4L)
      result <- .hmip_HighsSearch_evaluateNode(s)
      .hmip_Highs_setOptionValue(h, "simplex_strategy", 1L)
      if (result == "Open") {
        # ponytail: the adapter has no IPM run; "ipm" falls back to the simplex solver
        .hmip_Highs_clearSolver(h)
        .hmip_Highs_setOptionValue(h, "solver", "ipm")
        result <- .hmip_HighsSearch_evaluateNode(s)
        if (result == "Open") {
          currnode$opensubtrees <- 0L
          result <- "LpInfeasible"
        }
      }
    }
    return(result)
  }

  domchgPos <- .hmip_b2_stacksize(dom)
  pass <- .hmip_b2_orbitsValidInChildNode(s, currnode$branchingdecision)
  .hmip_HighsDomain_changeBound(dom, currnode$branchingdecision)
  currnode$opensubtrees <- 1L
  nd <- .hmip_b2_nodedata(max(childLb, currnode$lower_bound), currnode$estimate, currnode$nodeBasis,
                          if (pass) currnode$stabilizerOrbits)
  nd$domgchgStackPos <- domchgPos
  .hmip_b2_push(s, nd)
  "Branched"
}

# the pop loop shared by backtrack and backtrackPlunge; returns FALSE when
# the stack empties (after restoring the global domain)
.hmip_b2_pop_closed <- function(s, recoverBasis) {
  md <- s$mipsolver$mipdata_; dom <- s$localdom
  while (.hmip_b2_top(s)$opensubtrees == 0L) {
    s$countTreeWeight <- TRUE
    top <- .hmip_b2_top(s)
    s$depthoffset <- s$depthoffset + top$skipDepthCount
    if (length(s$nodestack) == 1L) {
      if (recoverBasis && !is.null(top$nodeBasis)) .hmip_HighsLpRelaxation_setStoredBasis(s$lp, top$nodeBasis)
      .hmip_b2_pop(s)
      .hmip_HighsDomain_backtrackToGlobal(dom)
      .hmip_HighsLpRelaxation_flushDomain(s$lp, dom)
      if (recoverBasis) .hmip_HighsLpRelaxation_recoverBasis(s$lp)
      return(FALSE)
    }
    .hmip_b2_pop(s)
    .hmip_HighsDomain_backtrack(dom)
    top <- .hmip_b2_top(s)
    if (top$opensubtrees != 0L) {
      s$countTreeWeight <- top$skipDepthCount == 0L
      oldNumDomchgs <- .hmip_HighsDomain_getNumDomainChanges(dom)
      oldNumChangedCols <- .hmip_b2_nchg(dom)
      .hmip_HighsDomain_propagate(dom)
      if (!.hmip_HighsDomain_infeasible(dom) && oldNumDomchgs != .hmip_HighsDomain_getNumDomainChanges(dom)) {
        if (!is.null(top$stabilizerOrbits)) .hmip_StabilizerOrbits_orbitalFixing(top$stabilizerOrbits, dom)
        else .hmip_HighsSymmetries_propagateOrbitopes(md$symmetries, dom)
      }
      if (.hmip_HighsDomain_infeasible(dom)) {
        .hmip_HighsDomain_clearChangedCols(dom, oldNumChangedCols)
        if (s$countTreeWeight) s$treeweight <- .cd_add(s$treeweight, 2^(-.hmip_HighsSearch_getCurrentDepth(s)))
        top$opensubtrees <- 0L
      }
    }
  }
  TRUE
}

# flip the branching decision of currnode to the other child
.hmip_b2_flip <- function(currnode) {
  bd <- currnode$branchingdecision
  fallbackbranch <- bd$boundval == currnode$branching_point
  if (bd$boundtype == .hmip_b2_lower) {
    bd$boundtype <- .hmip_b2_upper; bd$boundval <- floor(bd$boundval - 0.5)
  } else {
    bd$boundtype <- .hmip_b2_lower; bd$boundval <- ceiling(bd$boundval + 0.5)
  }
  currnode$branchingdecision <- bd
  fallbackbranch
}

# change to the other child and propagate; returns list(prune, nodelb, pass)
.hmip_b2_open_other <- function(s, currnode, numChangedCols) {
  md <- s$mipsolver$mipdata_; dom <- s$localdom
  pass <- .hmip_b2_orbitsValidInChildNode(s, currnode$branchingdecision)
  .hmip_HighsDomain_changeBound(dom, currnode$branchingdecision)
  nodelb <- max(currnode$lower_bound, currnode$other_child_lb)
  prune <- nodelb > .hmip_HighsSearch_getCutoffBound(s) || .hmip_HighsDomain_infeasible(dom)
  if (!prune) {
    .hmip_HighsDomain_propagate(dom)
    prune <- .hmip_HighsDomain_infeasible(dom)
    if (prune) .hmip_HighsDomain_conflictAnalysis(dom, md$conflictPool)
  }
  if (!prune) {
    .hmip_HighsSymmetries_propagateOrbitopes(md$symmetries, dom)
    prune <- .hmip_HighsDomain_infeasible(dom)
  }
  if (!prune && pass && !is.null(currnode$stabilizerOrbits)) {
    .hmip_StabilizerOrbits_orbitalFixing(currnode$stabilizerOrbits, dom)
    prune <- .hmip_HighsDomain_infeasible(dom)
  }
  if (prune) {
    .hmip_HighsDomain_backtrack(dom)
    .hmip_HighsDomain_clearChangedCols(dom, numChangedCols)
    if (s$countTreeWeight) s$treeweight <- .cd_add(s$treeweight, 2^(-.hmip_HighsSearch_getCurrentDepth(s)))
  }
  list(prune = prune, nodelb = nodelb, pass = pass)
}

.hmip_HighsSearch_backtrack <- function(s, recoverBasis = TRUE) {
  if (!length(s$nodestack)) return(FALSE)
  dom <- s$localdom
  repeat {
    if (!.hmip_b2_pop_closed(s, recoverBasis)) return(FALSE)
    currnode <- .hmip_b2_top(s)
    currnode$opensubtrees <- 0L
    fallbackbranch <- .hmip_b2_flip(currnode)
    domchgPos <- .hmip_b2_stacksize(dom)
    if (fallbackbranch) currnode$branching_point <- currnode$branchingdecision$boundval
    numChangedCols <- .hmip_b2_nchg(dom)
    o <- .hmip_b2_open_other(s, currnode, numChangedCols)
    if (o$prune) next
    nd <- .hmip_b2_nodedata(o$nodelb, currnode$estimate, currnode$nodeBasis, if (o$pass) currnode$stabilizerOrbits)
    .hmip_b2_push(s, nd)
    .hmip_HighsLpRelaxation_flushDomain(s$lp, dom)
    nd$domgchgStackPos <- domchgPos
    break
  }
  top <- .hmip_b2_top(s)
  if (recoverBasis && !is.null(top$nodeBasis)) {
    .hmip_HighsLpRelaxation_setStoredBasis(s$lp, top$nodeBasis)
    .hmip_HighsLpRelaxation_recoverBasis(s$lp)
  }
  TRUE
}

.hmip_HighsSearch_backtrackPlunge <- function(s, nodequeue) {
  if (!length(s$nodestack)) return(FALSE)
  md <- s$mipsolver$mipdata_; dom <- s$localdom; pc <- s$pseudocost
  repeat {
    if (!.hmip_b2_pop_closed(s, TRUE)) return(FALSE)
    currnode <- .hmip_b2_top(s)
    currnode$opensubtrees <- 0L
    bp0 <- currnode$branching_point
    fallbackbranch <- .hmip_b2_flip(currnode)
    bd <- currnode$branchingdecision
    nodeScore <- if (bd$boundtype == .hmip_b2_upper)
      .hmip_HighsPseudocost_getScoreDown(pc, bd$column, if (fallbackbranch) 0.5 else bp0)
    else .hmip_HighsPseudocost_getScoreUp(pc, bd$column, if (fallbackbranch) 0.5 else bp0)
    if (fallbackbranch) currnode$branching_point <- bd$boundval
    domchgPos <- .hmip_b2_stacksize(dom)
    numChangedCols <- .hmip_b2_nchg(dom)
    o <- .hmip_b2_open_other(s, currnode, numChangedCols)
    if (o$prune) next
    nodelb <- max(o$nodelb, .hmip_HighsDomain_getObjectiveLowerBound(dom))
    nodeToQueue <- nodelb > md$optimality_limit
    if (!nodeToQueue) {
      for (i in rev(seq_len(length(s$nodestack) - 1L))) {
        an <- s$nodestack[[i]]
        if (an$opensubtrees == 0L) next
        fb <- an$branchingdecision$boundval == an$branching_point
        branchpoint <- if (fb) 0.5 else an$branching_point
        acol <- an$branchingdecision$column
        if (an$branchingdecision$boundtype == .hmip_b2_lower) {
          inactive <- .hmip_HighsPseudocost_getScoreDown(pc, acol, branchpoint)
          active <- .hmip_HighsPseudocost_getScoreUp(pc, acol, branchpoint)
        } else {
          active <- .hmip_HighsPseudocost_getScoreDown(pc, acol, branchpoint)
          inactive <- .hmip_HighsPseudocost_getScoreUp(pc, acol, branchpoint)
        }
        nodeToQueue <- inactive - active > nodeScore + md$feastol
        break
      }
    }
    if (nodeToQueue) {
      red <- .hmip_HighsDomain_getReducedDomainChangeStack(dom)
      tw <- .hmip_HighsNodeQueue_emplaceNode(nodequeue, red$stack, red$branchingPositions, nodelb,
                                             .hmip_b2_top(s)$estimate, .hmip_HighsSearch_getCurrentDepth(s) + 1L)
      if (s$countTreeWeight) s$treeweight <- .cd_add(s$treeweight, tw)
      .hmip_HighsDomain_backtrack(dom)
      .hmip_HighsDomain_clearChangedCols(dom, numChangedCols)
      next
    }
    nd <- .hmip_b2_nodedata(nodelb, currnode$estimate, currnode$nodeBasis, if (o$pass) currnode$stabilizerOrbits)
    .hmip_b2_push(s, nd)
    .hmip_HighsLpRelaxation_flushDomain(s$lp, dom)
    nd$domgchgStackPos <- domchgPos
    break
  }
  top <- .hmip_b2_top(s)
  if (!is.null(top$nodeBasis)) {
    .hmip_HighsLpRelaxation_setStoredBasis(s$lp, top$nodeBasis)
    .hmip_HighsLpRelaxation_recoverBasis(s$lp)
  }
  TRUE
}

.hmip_HighsSearch_backtrackUntilDepth <- function(s, targetDepth) {
  if (!length(s$nodestack)) return(FALSE)
  dom <- s$localdom
  if (.hmip_HighsSearch_getCurrentDepth(s) >= targetDepth) .hmip_b2_top(s)$opensubtrees <- 0L
  while (.hmip_b2_top(s)$opensubtrees == 0L) {
    s$depthoffset <- s$depthoffset + .hmip_b2_top(s)$skipDepthCount
    .hmip_b2_pop(s)
    .hmip_HighsDomain_backtrack(dom)
    if (!length(s$nodestack)) {
      .hmip_HighsLpRelaxation_flushDomain(s$lp, dom)
      return(FALSE)
    }
    if (.hmip_HighsSearch_getCurrentDepth(s) >= targetDepth) .hmip_b2_top(s)$opensubtrees <- 0L
  }
  currnode <- .hmip_b2_top(s)
  currnode$opensubtrees <- 0L
  fallbackbranch <- .hmip_b2_flip(currnode)
  if (fallbackbranch) currnode$branching_point <- currnode$branchingdecision$boundval
  domchgPos <- .hmip_b2_stacksize(dom)
  pass <- .hmip_b2_orbitsValidInChildNode(s, currnode$branchingdecision)
  .hmip_HighsDomain_changeBound(dom, currnode$branchingdecision)
  nd <- .hmip_b2_nodedata(currnode$lower_bound, currnode$estimate, currnode$nodeBasis,
                          if (pass) currnode$stabilizerOrbits)
  .hmip_b2_push(s, nd)
  .hmip_HighsLpRelaxation_flushDomain(s$lp, dom)
  nd$domgchgStackPos <- domchgPos
  if (!is.null(nd$nodeBasis) && length(nd$nodeBasis$row_status) == .hmip_HighsLpRelaxation_numRows(s$lp))
    .hmip_HighsLpRelaxation_setStoredBasis(s$lp, nd$nodeBasis)
  .hmip_HighsLpRelaxation_recoverBasis(s$lp)
  TRUE
}

.hmip_HighsSearch_dive <- function(s) {
  s$reliableatnode <- new.env(parent = emptyenv())
  repeat {
    s$nnodes <- s$nnodes + 1
    result <- .hmip_HighsSearch_evaluateNode(s)
    if (.hmip_HighsMipSolverData_checkLimits(s$mipsolver$mipdata_, s$nnodes)) return(result)
    if (result != "Open") return(result)
    result <- .hmip_HighsSearch_branch(s)
    if (result != "Branched") return(result)
  }
}

.hmip_HighsSearch_solveDepthFirst <- function(s, maxbacktracks = 1) {
  repeat {
    if (maxbacktracks == 0) break
    result <- .hmip_HighsSearch_dive(s)
    if (result == "Open") break
    maxbacktracks <- maxbacktracks - 1
    if (!.hmip_HighsSearch_backtrack(s)) break
  }
  invisible()
}
