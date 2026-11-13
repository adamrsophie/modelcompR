# The domain of the HiGHS 1.12.0 MIP solver (mip/HighsDomain.cpp,
# mip/HighsDomainChange.h) with its propagators for model rows, cut pools,
# conflict pools and the objective, the conflict analysis of
# HighsDomain::ConflictSet, and the conflict pool (mip/HighsConflictPool.cpp),
# as SciPy 1.18.1 vendors them (release build: asserts and the debug solution
# hooks are no-ops). MIT license notice in R/optimization-scipy-highs-util.R.
# Every stored index is 0-based exactly as in C++ (columns, rows, stack
# positions, cuts, conflicts) and -1 keeps its C++ meaning; vectors are read
# as v[i + 1L]. HighsCDouble values are c(hi, lo).

.hmip_kLower <- 0L
.hmip_kUpper <- 1L
.hmip_Reason_kBranching <- -1L
.hmip_Reason_kUnknown <- -2L
.hmip_Reason_kModelRowUpper <- -3L
.hmip_Reason_kModelRowLower <- -4L
.hmip_Reason_kCliqueTable <- -5L
.hmip_Reason_kConflictingBounds <- -6L
.hmip_Reason_kObjective <- -7L

.hmip_reason <- function(type, index = 0L) list(type = as.integer(type), index = as.integer(index))
.hmip_reason_branching <- function() .hmip_reason(-1L, 0L)
.hmip_reason_unspecified <- function() .hmip_reason(-2L, 0L)
.hmip_reason_cliqueTable <- function(col, val) .hmip_reason(-5L, 2L * col + val)
.hmip_domchg <- function(boundval, column, boundtype) {
  list(boundval = boundval, column = as.integer(column), boundtype = as.integer(boundtype))
}
.hmip_domchg_eq <- function(a, b) {
  a$boundtype == b$boundtype && a$column == b$column && a$boundval == b$boundval
}

# HighsCDouble helpers on top of .cd_* (presolve) and the B1 util file
# a += v for a double v (two_sum(v, hi) order of operator+=(double))
.hmip_cd_iaddd <- function(a, v) .cd_iadd(a, c(v, 0))
# a -= b for HighsCDouble b
.hmip_cd_isub <- function(a, b) {
  r <- .cd_two_sum(-b[1L], a[1L])
  c(r[1L], (a[2L] + r[2L]) - b[2L])
}
# double v - HighsCDouble b (friend operator-(double, CD) is -b + v)
.hmip_cd_rsub <- function(v, b) .cd_add(-b, v)
.hmip_cd_renorm <- function(a) .cd_two_sum(a[1L], a[2L])

.hmip_activityContributionMin <- function(coef, lb, ub) {
  if (coef < 0) { if (ub == Inf) -Inf else coef * ub } else { if (lb == -Inf) -Inf else coef * lb }
}
.hmip_activityContributionMax <- function(coef, lb, ub) {
  if (coef < 0) { if (lb == -Inf) Inf else coef * lb } else { if (ub == Inf) Inf else coef * ub }
}
# computeDelta: list(d = delta (CD), dinf = change of numinfs)
.hmip_computeDelta <- function(val, oldbound, newbound, inf) {
  if (oldbound == inf) return(list(d = .cd_mul(c(newbound, 0), val), dinf = -1L))
  if (newbound == inf) return(list(d = .cd_mul(c(-oldbound, 0), val), dinf = 1L))
  list(d = .cd_mul(.cd_sub_cd(c(newbound, 0), c(oldbound, 0)), val), dinf = 0L)
}
.hmip_boundRange <- function(upper_bound, lower_bound, tolerance, var_type) {
  range <- upper_bound - lower_bound
  range - (if (var_type == 0L) max(0.3 * range, 1000.0 * tolerance) else tolerance)
}

.hmip_dom_vtype <- function(dom, col) {
  it <- dom$mipsolver$model_$integrality_
  if (length(it)) it[col + 1L] else 0L
}
.hmip_HighsDomain_variableType <- function(dom, col) .hmip_dom_vtype(dom, col)
.hmip_HighsDomain_feastol <- function(dom) dom$mipsolver$mipdata_$feastol

# ---------------------------------------------------------------------------
# HighsDomain
# ---------------------------------------------------------------------------

.hmip_new_HighsDomain <- function(mipsolver) {
  dom <- new.env(parent = emptyenv())
  ncol <- length(mipsolver$model_$col_lower_)
  dom$mipsolver <- mipsolver
  dom$col_lower_ <- mipsolver$model_$col_lower_
  dom$col_upper_ <- mipsolver$model_$col_upper_
  dom$colLowerPos_ <- rep(-1L, ncol)
  dom$colUpperPos_ <- rep(-1L, ncol)
  dom$changedcolsflags_ <- logical(ncol)
  dom$changedcols_ <- integer()
  dom$domchgstack_boundval <- numeric()
  dom$domchgstack_column <- integer()
  dom$domchgstack_boundtype <- integer()
  dom$domchgreason_type <- integer()
  dom$domchgreason_index <- integer()
  dom$prevboundval_val <- numeric()
  dom$prevboundval_pos <- integer()
  dom$branchPos_ <- integer()
  dom$activitymin_hi <- numeric(); dom$activitymin_lo <- numeric()
  dom$activitymax_hi <- numeric(); dom$activitymax_lo <- numeric()
  dom$activitymininf_ <- integer(); dom$activitymaxinf_ <- integer()
  dom$capacityThreshold_ <- numeric()
  dom$propagateflags_ <- logical()
  dom$propagateinds_ <- integer()
  dom$objProp_ <- NULL
  dom$cutpoolpropagation <- list()
  dom$conflictPoolPropagation <- list()
  dom$infeasible_ <- FALSE
  dom$infeasible_reason <- .hmip_reason_unspecified()
  dom$infeasible_pos <- 0L
  # ponytail: insertion ordered instead of HighsHashTable slot order; only
  # iterated when mip_lifting_for_probing != -1 (default -1)
  dom$redundantRows_ <- integer()
  dom$recordRedundantRows_ <- FALSE
  dom
}

.hmip_dom_copy_fields <- c("mipsolver", "col_lower_", "col_upper_", "colLowerPos_", "colUpperPos_",
  "changedcolsflags_", "changedcols_", "domchgstack_boundval", "domchgstack_column", "domchgstack_boundtype",
  "domchgreason_type", "domchgreason_index", "prevboundval_val", "prevboundval_pos", "branchPos_",
  "activitymin_hi", "activitymin_lo", "activitymax_hi", "activitymax_lo", "activitymininf_", "activitymaxinf_",
  "capacityThreshold_", "propagateflags_", "propagateinds_", "infeasible_", "infeasible_reason")

.hmip_copy_env <- function(src) {
  dst <- new.env(parent = emptyenv())
  for (nm in ls(src, all.names = TRUE)) assign(nm, get(nm, envir = src), envir = dst)
  dst
}

# HighsDomain(const HighsDomain& other)
.hmip_HighsDomain_copy <- function(other) {
  dom <- .hmip_new_HighsDomain(other$mipsolver)
  for (nm in .hmip_dom_copy_fields) assign(nm, get(nm, envir = other), envir = dom)
  dom$infeasible_pos <- other$infeasible_pos
  dom$objProp_ <- if (is.null(other$objProp_)) NULL else .hmip_copy_env(other$objProp_)
  if (!is.null(dom$objProp_)) dom$objProp_$domain <- dom
  dom$cutpoolpropagation <- lapply(other$cutpoolpropagation, function(p) .hmip_CutpoolPropagation_copy(p, dom))
  dom$conflictPoolPropagation <- lapply(other$conflictPoolPropagation, function(p) .hmip_ConflictPoolPropagation_copy(p, dom))
  dom
}

# HighsDomain::operator=; infeasible_pos and the redundant row state stay
.hmip_HighsDomain_assign <- function(dst, src) {
  for (nm in .hmip_dom_copy_fields) assign(nm, get(nm, envir = src), envir = dst)
  dst$objProp_ <- if (is.null(src$objProp_)) NULL else .hmip_copy_env(src$objProp_)
  if (!is.null(dst$objProp_)) dst$objProp_$domain <- dst
  assign_props <- function(dprops, sprops, copyfun, fields) {
    nd <- length(dprops); ns <- length(sprops)
    out <- vector("list", ns)
    for (i in seq_len(ns)) {
      if (i <= nd) {
        d <- dprops[[i]]
        for (nm in fields) assign(nm, get(nm, envir = sprops[[i]]), envir = d)
        d$domain <- dst
        out[[i]] <- d
      } else out[[i]] <- copyfun(sprops[[i]], dst)
    }
    if (nd > ns) for (i in (ns + 1L):nd) .hmip_prop_destroy(dprops[[i]])
    out
  }
  dst$cutpoolpropagation <- assign_props(dst$cutpoolpropagation, src$cutpoolpropagation,
    .hmip_CutpoolPropagation_copy, c("cutpoolindex", "cutpool", "activitycuts_hi", "activitycuts_lo",
      "activitycutsinf_", "propagatecutflags_", "propagatecutinds_", "capacityThreshold_"))
  dst$conflictPoolPropagation <- assign_props(dst$conflictPoolPropagation, src$conflictPoolPropagation,
    .hmip_ConflictPoolPropagation_copy, c("conflictpoolindex", "conflictpool_", "colLowerWatched_",
      "colUpperWatched_", "conflictFlag_", "propagateConflictInds_", "wl_boundval", "wl_column",
      "wl_boundtype", "wl_prev", "wl_next"))
  invisible(dst)
}

.hmip_prop_destroy <- function(p) {
  if (!is.null(p$cutpool)) .hmip_HighsCutPool_removePropagationDomain(p$cutpool, p)
  else .hmip_HighsConflictPool_removePropagationDomain(p$conflictpool_, p)
}

# the destructor: deregisters the propagators of a domain copy from their
# pools; call it where the C++ object leaves scope
.hmip_HighsDomain_destroy <- function(dom) {
  for (p in dom$cutpoolpropagation) .hmip_prop_destroy(p)
  for (p in dom$conflictPoolPropagation) .hmip_prop_destroy(p)
  dom$cutpoolpropagation <- list()
  dom$conflictPoolPropagation <- list()
  invisible(NULL)
}

.hmip_HighsDomain_addCutpool <- function(dom, cutpool) {
  idx <- length(dom$cutpoolpropagation)
  dom$cutpoolpropagation[[idx + 1L]] <- .hmip_new_CutpoolPropagation(idx, dom, cutpool)
  invisible(NULL)
}
.hmip_HighsDomain_addConflictPool <- function(dom, conflictPool) {
  idx <- length(dom$conflictPoolPropagation)
  dom$conflictPoolPropagation[[idx + 1L]] <- .hmip_new_ConflictPoolPropagation(idx, dom, conflictPool)
  invisible(NULL)
}

# simple accessors
.hmip_HighsDomain_infeasible <- function(dom) dom$infeasible_
.hmip_HighsDomain_markInfeasible <- function(dom, reason = .hmip_reason_unspecified()) {
  dom$infeasible_ <- TRUE
  dom$infeasible_pos <- length(dom$domchgstack_column)
  dom$infeasible_reason <- reason
  invisible(NULL)
}
.hmip_HighsDomain_isFixed <- function(dom, col) dom$col_lower_[col + 1L] == dom$col_upper_[col + 1L]
.hmip_HighsDomain_isBinary <- function(dom, col) {
  .hmip_dom_vtype(dom, col) != 0L && dom$col_lower_[col + 1L] == 0 && dom$col_upper_[col + 1L] == 1
}
.hmip_HighsDomain_isGlobalBinary <- function(dom, col) {
  m <- dom$mipsolver$model_
  .hmip_dom_vtype(dom, col) != 0L && m$col_lower_[col + 1L] == 0 && m$col_upper_[col + 1L] == 1
}
.hmip_HighsDomain_isActive <- function(dom, domchg) {
  if (domchg$boundtype == 0L) domchg$boundval <= dom$col_lower_[domchg$column + 1L]
  else domchg$boundval >= dom$col_upper_[domchg$column + 1L]
}
.hmip_HighsDomain_getChangedCols <- function(dom) dom$changedcols_
.hmip_HighsDomain_clearChangedCols <- function(dom, start = 0L) {
  n <- length(dom$changedcols_)
  if (start < n) {
    dom$changedcolsflags_[dom$changedcols_[(start + 1L):n] + 1L] <- FALSE
    dom$changedcols_ <- dom$changedcols_[seq_len(start)]
  }
  invisible(NULL)
}
.hmip_HighsDomain_removeContinuousChangedCols <- function(dom) {
  for (i in dom$changedcols_) dom$changedcolsflags_[i + 1L] <- .hmip_dom_vtype(dom, i) != 0L
  dom$changedcols_ <- dom$changedcols_[dom$changedcolsflags_[dom$changedcols_ + 1L]]
  invisible(NULL)
}
.hmip_HighsDomain_getBranchingPositions <- function(dom) dom$branchPos_
.hmip_HighsDomain_getPreviousBounds <- function(dom) list(first = dom$prevboundval_val, second = dom$prevboundval_pos)
.hmip_HighsDomain_getDomainChangeStack <- function(dom) {
  list(boundval = dom$domchgstack_boundval, column = dom$domchgstack_column, boundtype = dom$domchgstack_boundtype)
}
.hmip_HighsDomain_getDomainChangeReason <- function(dom) list(type = dom$domchgreason_type, index = dom$domchgreason_index)
.hmip_HighsDomain_getNumDomainChanges <- function(dom) length(dom$domchgstack_column)
.hmip_HighsDomain_getBranchDepth <- function(dom) length(dom$branchPos_)
.hmip_HighsDomain_colBoundsAreGlobal <- function(dom, col) {
  dom$colLowerPos_[col + 1L] == -1L && dom$colUpperPos_[col + 1L] == -1L
}
.hmip_dom_stackchg <- function(dom, k) {
  .hmip_domchg(dom$domchgstack_boundval[k + 1L], dom$domchgstack_column[k + 1L], dom$domchgstack_boundtype[k + 1L])
}
.hmip_HighsDomain_getMinActivity <- function(dom, row) {
  if (dom$activitymininf_[row + 1L] == 0L) dom$activitymin_hi[row + 1L] + dom$activitymin_lo[row + 1L] else -Inf
}
.hmip_HighsDomain_getMaxActivity <- function(dom, row) {
  if (dom$activitymaxinf_[row + 1L] == 0L) dom$activitymax_hi[row + 1L] + dom$activitymax_lo[row + 1L] else Inf
}
.hmip_dom_amin <- function(dom, r) c(dom$activitymin_hi[r + 1L], dom$activitymin_lo[r + 1L])
.hmip_dom_amax <- function(dom, r) c(dom$activitymax_hi[r + 1L], dom$activitymax_lo[r + 1L])
.hmip_dom_set_amin <- function(dom, r, v) { dom$activitymin_hi[r + 1L] <- v[1L]; dom$activitymin_lo[r + 1L] <- v[2L] }
.hmip_dom_set_amax <- function(dom, r, v) { dom$activitymax_hi[r + 1L] <- v[1L]; dom$activitymax_lo[r + 1L] <- v[2L] }

.hmip_HighsDomain_getReducedDomainChangeStack <- function(dom) {
  n <- length(dom$domchgstack_column)
  keep <- integer(); bp <- integer()
  for (i in seq_len(n) - 1L) {
    col <- dom$domchgstack_column[i + 1L]
    if ((dom$domchgstack_boundtype[i + 1L] == 0L && dom$colLowerPos_[col + 1L] != i) ||
        (dom$domchgstack_boundtype[i + 1L] == 1L && dom$colUpperPos_[col + 1L] != i)) next
    if (dom$domchgreason_type[i + 1L] == -1L) {
      bp <- c(bp, length(keep))
    } else {
      k <- i
      while (dom$prevboundval_pos[k + 1L] != -1L) {
        k <- dom$prevboundval_pos[k + 1L]
        if (dom$domchgreason_type[k + 1L] == -1L) { bp <- c(bp, length(keep)); break }
      }
    }
    keep <- c(keep, i)
  }
  list(stack = list(boundval = dom$domchgstack_boundval[keep + 1L], column = dom$domchgstack_column[keep + 1L],
                    boundtype = dom$domchgstack_boundtype[keep + 1L]),
       branchingPositions = bp)
}

.hmip_HighsDomain_isFixing <- function(dom, domchg) {
  other <- if (domchg$boundtype == 1L) dom$col_lower_[domchg$column + 1L] else dom$col_upper_[domchg$column + 1L]
  abs(domchg$boundval - other) <= dom$mipsolver$mipdata_$epsilon
}

.hmip_HighsDomain_flip <- function(dom, domchg) {
  ft <- dom$mipsolver$mipdata_$feastol
  if (domchg$boundtype == 0L) {
    v <- domchg$boundval - ft
    if (.hmip_dom_vtype(dom, domchg$column) != 0L) v <- floor(v)
    .hmip_domchg(v, domchg$column, 1L)
  } else {
    v <- domchg$boundval + ft
    if (.hmip_dom_vtype(dom, domchg$column) != 0L) v <- ceiling(v)
    .hmip_domchg(v, domchg$column, 0L)
  }
}

# computeMinActivity / computeMaxActivity over entries start..end-1 of the
# 0-based arrays ARindex / ARvalue; list(ninf, activity)
.hmip_HighsDomain_computeMinActivity <- function(dom, start, end, ARindex, ARvalue) {
  act <- c(0, 0); ninf <- 0L
  j <- start
  while (j != end) {
    col <- ARindex[j + 1L]; val <- ARvalue[j + 1L]
    if (dom$infeasible_) {
      lb <- .hmip_HighsDomain_getColLowerPos(dom, col, dom$infeasible_pos - 1L)$bound
      ub <- .hmip_HighsDomain_getColUpperPos(dom, col, dom$infeasible_pos - 1L)$bound
    } else {
      lb <- dom$col_lower_[col + 1L]; ub <- dom$col_upper_[col + 1L]
    }
    cm <- .hmip_activityContributionMin(val, lb, ub)
    if (cm == -Inf) ninf <- ninf + 1L else act <- .hmip_cd_iaddd(act, cm)
    j <- j + 1L
  }
  list(ninf = ninf, activity = .hmip_cd_renorm(act))
}
.hmip_HighsDomain_computeMaxActivity <- function(dom, start, end, ARindex, ARvalue) {
  act <- c(0, 0); ninf <- 0L
  j <- start
  while (j != end) {
    col <- ARindex[j + 1L]; val <- ARvalue[j + 1L]
    if (dom$infeasible_) {
      lb <- .hmip_HighsDomain_getColLowerPos(dom, col, dom$infeasible_pos - 1L)$bound
      ub <- .hmip_HighsDomain_getColUpperPos(dom, col, dom$infeasible_pos - 1L)$bound
    } else {
      lb <- dom$col_lower_[col + 1L]; ub <- dom$col_upper_[col + 1L]
    }
    cm <- .hmip_activityContributionMax(val, lb, ub)
    if (cm == Inf) ninf <- ninf + 1L else act <- .hmip_cd_iaddd(act, cm)
    j <- j + 1L
  }
  list(ninf = ninf, activity = .hmip_cd_renorm(act))
}

# adjustedUb / adjustedLb with boundVal a HighsCDouble; list(bound, accept)
.hmip_HighsDomain_adjustedUb <- function(dom, col, boundVal) {
  md <- dom$mipsolver$mipdata_
  lo <- dom$col_lower_[col + 1L]; up <- dom$col_upper_[col + 1L]
  if (.hmip_dom_vtype(dom, col) != 0L) {
    bound <- .cd_val(.cd_floor(.cd_add(boundVal, md$feastol)))
    accept <- bound < up && up - bound > 1000.0 * md$feastol * abs(bound)
  } else {
    bv <- .cd_val(boundVal)
    bound <- if (abs(bv - lo) <= md$epsilon) lo else bv
    if (up == Inf) accept <- TRUE
    else if (bound + 1000.0 * md$feastol < up) {
      rel <- up - bound
      if (lo != -Inf) rel <- rel / (up - lo) else rel <- rel / max(abs(up), abs(bound))
      accept <- rel >= 0.3
    } else accept <- FALSE
  }
  list(bound = bound, accept = accept)
}
.hmip_HighsDomain_adjustedLb <- function(dom, col, boundVal) {
  md <- dom$mipsolver$mipdata_
  lo <- dom$col_lower_[col + 1L]; up <- dom$col_upper_[col + 1L]
  if (.hmip_dom_vtype(dom, col) != 0L) {
    bound <- .cd_val(.cd_ceil(.cd_sub(boundVal, md$feastol)))
    accept <- bound > lo && bound - lo > 1000.0 * md$feastol * abs(bound)
  } else {
    bv <- .cd_val(boundVal)
    bound <- if (abs(up - bv) <= md$epsilon) up else bv
    if (lo == -Inf) accept <- TRUE
    else if (bound - 1000.0 * md$feastol > lo) {
      rel <- bound - lo
      if (up != Inf) rel <- rel / (up - lo) else rel <- rel / max(abs(lo), abs(bound))
      accept <- rel >= 0.3
    } else accept <- FALSE
  }
  list(bound = bound, accept = accept)
}

.hmip_kHighsTiny <- 1e-14

# propagateRowUpper / propagateRowLower: the bound changes as an SoA list
.hmip_HighsDomain_propagateRowUpper <- function(dom, Rindex, Rvalue, Rlen, Rupper, minactivity, ninfmin) {
  out <- list(boundval = numeric(), column = integer(), boundtype = integer())
  if (ninfmin > 1L) return(out)
  ft <- dom$mipsolver$mipdata_$feastol
  for (i in seq_len(Rlen)) {
    col <- Rindex[i]; v <- Rvalue[i]
    ac <- .hmip_activityContributionMin(v, dom$col_lower_[col + 1L], dom$col_upper_[col + 1L])
    if (ninfmin == 1L) {
      if (ac != -Inf) next
      minresact <- minactivity
    } else minresact <- .cd_sub(minactivity, ac)
    boundVal <- .cd_div(.hmip_cd_rsub(Rupper, minresact), v)
    if (abs(.cd_val(boundVal) * .hmip_kHighsTiny) > ft) next
    if (v > 0) {
      a <- .hmip_HighsDomain_adjustedUb(dom, col, boundVal)
      if (a$accept) { out$boundval <- c(out$boundval, a$bound); out$column <- c(out$column, col); out$boundtype <- c(out$boundtype, 1L) }
    } else {
      a <- .hmip_HighsDomain_adjustedLb(dom, col, boundVal)
      if (a$accept) { out$boundval <- c(out$boundval, a$bound); out$column <- c(out$column, col); out$boundtype <- c(out$boundtype, 0L) }
    }
  }
  out
}
.hmip_HighsDomain_propagateRowLower <- function(dom, Rindex, Rvalue, Rlen, Rlower, maxactivity, ninfmax) {
  out <- list(boundval = numeric(), column = integer(), boundtype = integer())
  if (ninfmax > 1L) return(out)
  ft <- dom$mipsolver$mipdata_$feastol
  for (i in seq_len(Rlen)) {
    col <- Rindex[i]; v <- Rvalue[i]
    ac <- .hmip_activityContributionMax(v, dom$col_lower_[col + 1L], dom$col_upper_[col + 1L])
    if (ninfmax == 1L) {
      if (ac != Inf) next
      maxresact <- maxactivity
    } else maxresact <- .cd_sub(maxactivity, ac)
    boundVal <- .cd_div(.hmip_cd_rsub(Rlower, maxresact), v)
    if (abs(.cd_val(boundVal) * .hmip_kHighsTiny) > ft) next
    if (v < 0) {
      a <- .hmip_HighsDomain_adjustedUb(dom, col, boundVal)
      if (a$accept) { out$boundval <- c(out$boundval, a$bound); out$column <- c(out$column, col); out$boundtype <- c(out$boundtype, 1L) }
    } else {
      a <- .hmip_HighsDomain_adjustedLb(dom, col, boundVal)
      if (a$accept) { out$boundval <- c(out$boundval, a$bound); out$column <- c(out$column, col); out$boundtype <- c(out$boundtype, 0L) }
    }
  }
  out
}

# updateThresholdLbChange / UbChange return the new threshold
.hmip_dom_thrLb <- function(dom, col, newbound, val, threshold) {
  up <- dom$col_upper_[col + 1L]
  if (newbound != up) {
    ft <- dom$mipsolver$mipdata_$feastol
    tn <- abs(val) * .hmip_boundRange(up, newbound, ft, .hmip_dom_vtype(dom, col))
    threshold <- max(threshold, tn, ft)
  }
  threshold
}
.hmip_dom_thrUb <- function(dom, col, newbound, val, threshold) {
  lo <- dom$col_lower_[col + 1L]
  if (newbound != lo) {
    ft <- dom$mipsolver$mipdata_$feastol
    tn <- abs(val) * .hmip_boundRange(newbound, lo, ft, .hmip_dom_vtype(dom, col))
    threshold <- max(threshold, tn, ft)
  }
  threshold
}

.hmip_dom_setinf <- function(dom, type, index) {
  dom$infeasible_ <- TRUE
  dom$infeasible_pos <- length(dom$domchgstack_column)
  dom$infeasible_reason <- .hmip_reason(type, index)
}

# updateActivityLbChange (lb = TRUE) and updateActivityUbChange (lb = FALSE)
.hmip_dom_updateActivity <- function(dom, col, oldbound, newbound, lb) {
  ms <- dom$mipsolver; mip <- ms$model_; md <- ms$mipdata_
  A <- mip$a_matrix_
  start <- A$start_[col + 1L]; end <- A$start_[col + 2L]
  op <- dom$objProp_
  if (!is.null(op)) {
    if (lb) .hmip_ObjProp_updateActivityLbChange(op, col, oldbound, newbound)
    else .hmip_ObjProp_updateActivityUbChange(op, col, oldbound, newbound)
    if (dom$infeasible_) return(invisible(NULL))
  }
  inf <- if (lb) -Inf else Inf
  ft <- md$feastol
  i <- start
  while (i != end) {
    v <- A$value_[i + 1L]; r <- A$index_[i + 1L]
    rl <- mip$row_lower_[r + 1L]; ru <- mip$row_upper_[r + 1L]
    # for a lower bound change a positive coefficient moves the min activity,
    # for an upper bound change a positive coefficient moves the max activity
    usemin <- (v > 0) == lb
    if (usemin) {
      cd <- .hmip_computeDelta(v, oldbound, newbound, inf)
      dom$activitymininf_[r + 1L] <- dom$activitymininf_[r + 1L] + cd$dinf
      .hmip_dom_set_amin(dom, r, .cd_iadd(.hmip_dom_amin(dom, r), cd$d))
      if (dom$recordRedundantRows_ && rl != -Inf && ru == Inf)
        .hmip_dom_updateRedundantRows(dom, r, 1L, dom$activitymininf_[r + 1L], .hmip_dom_amin(dom, r), rl)
      if (.cd_val(cd$d) <= 0) {
        dom$capacityThreshold_[r + 1L] <- if (lb) .hmip_dom_thrLb(dom, col, newbound, v, dom$capacityThreshold_[r + 1L])
          else .hmip_dom_thrUb(dom, col, newbound, v, dom$capacityThreshold_[r + 1L])
        i <- i + 1L
        next
      }
      if (ru != Inf && dom$activitymininf_[r + 1L] == 0L && .cd_val(.cd_sub(.hmip_dom_amin(dom, r), ru)) > ft) {
        .hmip_dom_setinf(dom, -3L, r)
        end <- i + 1L
        break
      }
      if (dom$activitymininf_[r + 1L] <= 1L && !dom$propagateflags_[r + 1L] && ru != Inf)
        .hmip_HighsDomain_markPropagate(dom, r)
    } else {
      cd <- .hmip_computeDelta(v, oldbound, newbound, inf)
      dom$activitymaxinf_[r + 1L] <- dom$activitymaxinf_[r + 1L] + cd$dinf
      .hmip_dom_set_amax(dom, r, .cd_iadd(.hmip_dom_amax(dom, r), cd$d))
      if (dom$recordRedundantRows_ && rl == -Inf && ru != Inf)
        .hmip_dom_updateRedundantRows(dom, r, -1L, dom$activitymaxinf_[r + 1L], .hmip_dom_amax(dom, r), ru)
      if (.cd_val(cd$d) >= 0) {
        dom$capacityThreshold_[r + 1L] <- if (lb) .hmip_dom_thrLb(dom, col, newbound, v, dom$capacityThreshold_[r + 1L])
          else .hmip_dom_thrUb(dom, col, newbound, v, dom$capacityThreshold_[r + 1L])
        i <- i + 1L
        next
      }
      if (rl != -Inf && dom$activitymaxinf_[r + 1L] == 0L && .cd_val(.hmip_cd_rsub(rl, .hmip_dom_amax(dom, r))) > ft) {
        .hmip_dom_setinf(dom, -4L, r)
        end <- i + 1L
        break
      }
      if (dom$activitymaxinf_[r + 1L] <= 1L && !dom$propagateflags_[r + 1L] && rl != -Inf)
        .hmip_HighsDomain_markPropagate(dom, r)
    }
    i <- i + 1L
  }
  if (!dom$infeasible_) {
    for (cp in dom$cutpoolpropagation) {
      if (lb) .hmip_CutpoolPropagation_updateActivityLbChange(cp, col, oldbound, newbound)
      else .hmip_CutpoolPropagation_updateActivityUbChange(cp, col, oldbound, newbound)
    }
  }
  if (dom$infeasible_) {
    t <- oldbound; oldbound <- newbound; newbound <- t
    i <- start
    while (i != end) {
      v <- A$value_[i + 1L]; r <- A$index_[i + 1L]
      cd <- .hmip_computeDelta(v, oldbound, newbound, inf)
      if ((v > 0) == lb) {
        dom$activitymininf_[r + 1L] <- dom$activitymininf_[r + 1L] + cd$dinf
        .hmip_dom_set_amin(dom, r, .cd_iadd(.hmip_dom_amin(dom, r), cd$d))
      } else {
        dom$activitymaxinf_[r + 1L] <- dom$activitymaxinf_[r + 1L] + cd$dinf
        .hmip_dom_set_amax(dom, r, .cd_iadd(.hmip_dom_amax(dom, r), cd$d))
      }
      i <- i + 1L
    }
    if (!is.null(op)) {
      if (lb) .hmip_ObjProp_updateActivityLbChange(op, col, oldbound, newbound)
      else .hmip_ObjProp_updateActivityUbChange(op, col, oldbound, newbound)
    }
    return(invisible(NULL))
  }
  for (cp in dom$conflictPoolPropagation) {
    if (lb) .hmip_ConflictPoolPropagation_updateActivityLbChange(cp, col, oldbound, newbound)
    else .hmip_ConflictPoolPropagation_updateActivityUbChange(cp, col, oldbound, newbound)
  }
  invisible(NULL)
}

.hmip_dom_recomputeCapacityThreshold <- function(dom, row) {
  md <- dom$mipsolver$mipdata_
  ft <- md$feastol
  start <- md$ARstart_[row + 1L]; end <- md$ARstart_[row + 2L]
  thr <- -ft
  i <- start
  while (i < end) {
    col <- md$ARindex_[i + 1L]
    if (dom$col_upper_[col + 1L] != dom$col_lower_[col + 1L]) {
      t <- abs(md$ARvalue_[i + 1L]) * .hmip_boundRange(dom$col_upper_[col + 1L], dom$col_lower_[col + 1L], ft, .hmip_dom_vtype(dom, col))
      thr <- max(thr, t, ft)
    }
    i <- i + 1L
  }
  dom$capacityThreshold_[row + 1L] <- thr
}

.hmip_dom_updateRedundantRows <- function(dom, row, direction, numinf, activity, bound) {
  if (numinf != 0L || .cd_val(.cd_mul(activity, direction)) <= direction * bound + dom$mipsolver$mipdata_$feastol) return()
  if (!(row %in% dom$redundantRows_)) dom$redundantRows_ <- c(dom$redundantRows_, row)
}
.hmip_HighsDomain_clearRedundantRows <- function(dom) dom$redundantRows_ <- integer()
.hmip_HighsDomain_getRedundantRows <- function(dom) dom$redundantRows_
.hmip_HighsDomain_setRecordRedundantRows <- function(dom, val) dom$recordRedundantRows_ <- val
.hmip_HighsDomain_getRedundantRowValue <- function(dom, row) {
  m <- dom$mipsolver$model_
  if (m$row_lower_[row + 1L] != -Inf) .cd_val(.cd_sub(.hmip_dom_amin(dom, row), m$row_lower_[row + 1L]))
  else .cd_val(.cd_sub(.hmip_dom_amax(dom, row), m$row_upper_[row + 1L]))
}

.hmip_HighsDomain_markPropagateCut <- function(dom, reason) {
  t <- reason$type
  if (t < 0L) return(invisible(NULL))
  ncp <- length(dom$cutpoolpropagation)
  if (t < ncp) .hmip_CutpoolPropagation_markPropagateCut(dom$cutpoolpropagation[[t + 1L]], reason$index)
  else .hmip_ConflictPoolPropagation_markPropagateConflict(dom$conflictPoolPropagation[[t - ncp + 1L]], reason$index)
  invisible(NULL)
}

.hmip_HighsDomain_markPropagate <- function(dom, row) {
  if (dom$propagateflags_[row + 1L]) return(invisible(NULL))
  m <- dom$mipsolver$model_; ft <- dom$mipsolver$mipdata_$feastol
  rl <- m$row_lower_[row + 1L]; ru <- m$row_upper_[row + 1L]
  amin <- .cd_val(.hmip_dom_amin(dom, row)); amax <- .cd_val(.hmip_dom_amax(dom, row))
  ct <- dom$capacityThreshold_[row + 1L]
  proplower <- rl != -Inf && (dom$activitymininf_[row + 1L] != 0L || amin < rl - ft) &&
    (dom$activitymaxinf_[row + 1L] == 1L || (amax - rl) <= ct)
  propupper <- ru != Inf && (dom$activitymaxinf_[row + 1L] != 0L || amax > ru + ft) &&
    (dom$activitymininf_[row + 1L] == 1L || (ru - amin) <= ct)
  if (proplower || propupper) {
    dom$propagateinds_ <- c(dom$propagateinds_, row)
    dom$propagateflags_[row + 1L] <- TRUE
  }
  invisible(NULL)
}

.hmip_HighsDomain_computeRowActivities <- function(dom) {
  ms <- dom$mipsolver; md <- ms$mipdata_; m <- ms$model_
  nrow <- length(m$row_lower_)
  dom$activitymin_hi <- numeric(nrow); dom$activitymin_lo <- numeric(nrow)
  dom$activitymax_hi <- numeric(nrow); dom$activitymax_lo <- numeric(nrow)
  dom$activitymininf_ <- integer(nrow); dom$activitymaxinf_ <- integer(nrow)
  dom$capacityThreshold_ <- numeric(nrow)
  if (length(dom$propagateflags_) < nrow) dom$propagateflags_ <- c(dom$propagateflags_, logical(nrow - length(dom$propagateflags_)))
  for (i in seq_len(nrow) - 1L) {
    s <- md$ARstart_[i + 1L]; e <- md$ARstart_[i + 2L]
    a <- .hmip_HighsDomain_computeMinActivity(dom, s, e, md$ARindex_, md$ARvalue_)
    dom$activitymininf_[i + 1L] <- a$ninf; .hmip_dom_set_amin(dom, i, a$activity)
    a <- .hmip_HighsDomain_computeMaxActivity(dom, s, e, md$ARindex_, md$ARvalue_)
    dom$activitymaxinf_[i + 1L] <- a$ninf; .hmip_dom_set_amax(dom, i, a$activity)
    .hmip_dom_recomputeCapacityThreshold(dom, i)
    if ((dom$activitymininf_[i + 1L] <= 1L && m$row_upper_[i + 1L] != Inf) ||
        (dom$activitymaxinf_[i + 1L] <= 1L && m$row_lower_[i + 1L] != -Inf))
      .hmip_HighsDomain_markPropagate(dom, i)
  }
  invisible(NULL)
}

.hmip_dom_doChangeBound <- function(dom, boundval, col, boundtype) {
  if (boundtype == 0L) {
    oldbound <- dom$col_lower_[col + 1L]
    dom$col_lower_[col + 1L] <- boundval
    if (oldbound != boundval) {
      if (!dom$infeasible_) .hmip_dom_updateActivity(dom, col, oldbound, boundval, TRUE)
      if (!dom$changedcolsflags_[col + 1L]) {
        dom$changedcolsflags_[col + 1L] <- TRUE
        dom$changedcols_ <- c(dom$changedcols_, col)
      }
    }
  } else {
    oldbound <- dom$col_upper_[col + 1L]
    dom$col_upper_[col + 1L] <- boundval
    if (oldbound != boundval) {
      if (!dom$infeasible_) .hmip_dom_updateActivity(dom, col, oldbound, boundval, FALSE)
      if (!dom$changedcolsflags_[col + 1L]) {
        dom$changedcolsflags_[col + 1L] <- TRUE
        dom$changedcols_ <- c(dom$changedcols_, col)
      }
    }
  }
  oldbound
}

.hmip_HighsDomain_changeBound <- function(dom, boundchg, reason = .hmip_reason_branching()) {
  bv <- boundchg$boundval; col <- boundchg$column; bt <- boundchg$boundtype
  ft <- dom$mipsolver$mipdata_$feastol
  n <- length(dom$domchgstack_column)
  if (bt == 0L) {
    if (bv <= dom$col_lower_[col + 1L]) {
      if (reason$type != -1L) return(invisible(NULL))
      bv <- dom$col_lower_[col + 1L]
    }
    if (bv > dom$col_upper_[col + 1L]) {
      if (bv - dom$col_upper_[col + 1L] > ft) {
        if (!dom$infeasible_) {
          dom$infeasible_pos <- n
          dom$infeasible_ <- TRUE
          dom$infeasible_reason <- .hmip_reason(-6L, n)
        }
      } else {
        bv <- dom$col_upper_[col + 1L]
        if (bv == dom$col_lower_[col + 1L]) return(invisible(NULL))
      }
    }
    prevPos <- dom$colLowerPos_[col + 1L]
    dom$colLowerPos_[col + 1L] <- n
  } else {
    if (bv >= dom$col_upper_[col + 1L]) {
      if (reason$type != -1L) return(invisible(NULL))
      bv <- dom$col_upper_[col + 1L]
    }
    if (bv < dom$col_lower_[col + 1L]) {
      if (dom$col_lower_[col + 1L] - bv > ft) {
        if (!dom$infeasible_) {
          dom$infeasible_pos <- n
          dom$infeasible_ <- TRUE
          dom$infeasible_reason <- .hmip_reason(-6L, n)
        }
      } else {
        bv <- dom$col_lower_[col + 1L]
        if (bv == dom$col_upper_[col + 1L]) return(invisible(NULL))
      }
    }
    prevPos <- dom$colUpperPos_[col + 1L]
    dom$colUpperPos_[col + 1L] <- n
  }
  if (reason$type == -1L) dom$branchPos_ <- c(dom$branchPos_, n)
  binary <- .hmip_HighsDomain_isBinary(dom, col)
  oldbound <- .hmip_dom_doChangeBound(dom, bv, col, bt)
  dom$prevboundval_val[n + 1L] <- oldbound
  dom$prevboundval_pos[n + 1L] <- prevPos
  dom$domchgstack_boundval[n + 1L] <- bv
  dom$domchgstack_column[n + 1L] <- col
  dom$domchgstack_boundtype[n + 1L] <- bt
  dom$domchgreason_type[n + 1L] <- reason$type
  dom$domchgreason_index[n + 1L] <- reason$index
  if (binary && !dom$infeasible_ && .hmip_HighsDomain_isFixed(dom, col))
    .hmip_HighsCliqueTable_addImplications(dom$mipsolver$mipdata_$cliquetable, dom, col,
                                          as.integer(dom$col_lower_[col + 1L] > 0.5))
  invisible(NULL)
}
.hmip_HighsDomain_changeBound3 <- function(dom, boundtype, col, boundval, reason = .hmip_reason_branching()) {
  .hmip_HighsDomain_changeBound(dom, .hmip_domchg(boundval, col, boundtype), reason)
}

.hmip_HighsDomain_fixCol <- function(dom, col, val, reason = .hmip_reason_unspecified()) {
  if (dom$col_lower_[col + 1L] < val) {
    .hmip_HighsDomain_changeBound(dom, .hmip_domchg(val, col, 0L), reason)
    if (!dom$infeasible_) .hmip_HighsDomain_propagate(dom)
  }
  if (!dom$infeasible_ && dom$col_upper_[col + 1L] > val)
    .hmip_HighsDomain_changeBound(dom, .hmip_domchg(val, col, 1L), reason)
  invisible(NULL)
}

.hmip_dom_reset_stack <- function(dom) {
  n <- length(dom$domchgstack_column)
  for (k in seq_len(n)) {
    if (dom$domchgstack_boundtype[k] == 0L) dom$colLowerPos_[dom$domchgstack_column[k] + 1L] <- -1L
    else dom$colUpperPos_[dom$domchgstack_column[k] + 1L] <- -1L
  }
  dom$prevboundval_val <- numeric(); dom$prevboundval_pos <- integer()
  dom$domchgstack_boundval <- numeric(); dom$domchgstack_column <- integer(); dom$domchgstack_boundtype <- integer()
  dom$domchgreason_type <- integer(); dom$domchgreason_index <- integer()
  dom$branchPos_ <- integer()
}

.hmip_dom_redundant <- function(dom, bv, col, bt) {
  (bt == 0L && bv <= dom$col_lower_[col + 1L]) || (bt == 1L && bv >= dom$col_upper_[col + 1L])
}

.hmip_HighsDomain_setDomainChangeStack <- function(dom, stack) {
  dom$infeasible_ <- FALSE
  .hmip_dom_reset_stack(dom)
  for (k in seq_along(stack$column)) {
    bv <- stack$boundval[k]; col <- stack$column[k]; bt <- stack$boundtype[k]
    if (.hmip_dom_redundant(dom, bv, col, bt)) next
    .hmip_HighsDomain_changeBound(dom, .hmip_domchg(bv, col, bt), .hmip_reason_unspecified())
    if (dom$infeasible_) break
  }
  invisible(NULL)
}

.hmip_HighsDomain_setDomainChangeStack2 <- function(dom, stack, branchingPositions) {
  dom$infeasible_ <- FALSE
  .hmip_dom_reset_stack(dom)
  stacksize <- length(stack$column)
  k <- 0L
  for (branchPos in branchingPositions) {
    while (k < branchPos) {
      bv <- stack$boundval[k + 1L]; col <- stack$column[k + 1L]; bt <- stack$boundtype[k + 1L]
      if (!.hmip_dom_redundant(dom, bv, col, bt)) {
        .hmip_HighsDomain_changeBound(dom, .hmip_domchg(bv, col, bt), .hmip_reason_unspecified())
        if (!dom$infeasible_) .hmip_HighsDomain_propagate(dom)
        if (dom$infeasible_) return(invisible(NULL))
      }
      k <- k + 1L
    }
    if (k == stacksize) return(invisible(NULL))
    bv <- stack$boundval[k + 1L]; col <- stack$column[k + 1L]; bt <- stack$boundtype[k + 1L]
    if (bt == 0L) {
      if (bv <= dom$col_lower_[col + 1L]) {
        if (bv < dom$col_lower_[col + 1L]) next
        if (dom$colLowerPos_[col + 1L] != -1L) next
      }
    } else {
      if (bv >= dom$col_upper_[col + 1L]) {
        if (bv > dom$col_upper_[col + 1L]) next
        if (dom$colUpperPos_[col + 1L] != -1L) next
      }
    }
    .hmip_HighsDomain_changeBound(dom, .hmip_domchg(bv, col, bt), .hmip_reason_branching())
    if (!dom$infeasible_) .hmip_HighsDomain_propagate(dom)
    if (dom$infeasible_) return(invisible(NULL))
  }
  while (k < stacksize) {
    bv <- stack$boundval[k + 1L]; col <- stack$column[k + 1L]; bt <- stack$boundtype[k + 1L]
    if (!.hmip_dom_redundant(dom, bv, col, bt)) {
      .hmip_HighsDomain_changeBound(dom, .hmip_domchg(bv, col, bt), .hmip_reason_unspecified())
      if (!dom$infeasible_) .hmip_HighsDomain_propagate(dom)
      if (dom$infeasible_) break
    }
    k <- k + 1L
  }
  invisible(NULL)
}

# one step of the backtrack loops: undo stack entry k
.hmip_dom_undo <- function(dom, k, always) {
  prevbound <- dom$prevboundval_val[k + 1L]
  prevpos <- dom$prevboundval_pos[k + 1L]
  col <- dom$domchgstack_column[k + 1L]; bt <- dom$domchgstack_boundtype[k + 1L]
  if (bt == 0L) dom$colLowerPos_[col + 1L] <- prevpos else dom$colUpperPos_[col + 1L] <- prevpos
  if (always || prevbound != dom$domchgstack_boundval[k + 1L]) .hmip_dom_doChangeBound(dom, prevbound, col, bt)
  if (dom$infeasible_ && dom$infeasible_pos == k) {
    dom$infeasible_ <- FALSE
    dom$infeasible_reason <- .hmip_reason_unspecified()
  }
}

.hmip_HighsDomain_backtrackToGlobal <- function(dom) {
  n <- length(dom$domchgstack_column)
  k <- n - 1L
  old_infeasible <- dom$infeasible_
  old_reason <- dom$infeasible_reason
  if (dom$infeasible_ && dom$infeasible_pos == n) {
    dom$infeasible_ <- FALSE
    dom$infeasible_reason <- .hmip_reason_unspecified()
  }
  while (k >= 0L) {
    .hmip_dom_undo(dom, k, FALSE)
    k <- k - 1L
  }
  if (old_infeasible) {
    .hmip_HighsDomain_markPropagateCut(dom, old_reason)
    dom$infeasible_reason <- .hmip_reason_unspecified()
    dom$infeasible_ <- FALSE
  }
  numreason <- length(dom$domchgreason_type)
  i <- k + 1L
  while (i < numreason) {
    .hmip_HighsDomain_markPropagateCut(dom, .hmip_reason(dom$domchgreason_type[i + 1L], dom$domchgreason_index[i + 1L]))
    i <- i + 1L
  }
  dom$domchgstack_boundval <- numeric(); dom$domchgstack_column <- integer(); dom$domchgstack_boundtype <- integer()
  dom$prevboundval_val <- numeric(); dom$prevboundval_pos <- integer()
  dom$domchgreason_type <- integer(); dom$domchgreason_index <- integer()
  dom$branchPos_ <- integer()
  invisible(NULL)
}

.hmip_HighsDomain_backtrack <- function(dom) {
  n <- length(dom$domchgstack_column)
  k <- n - 1L
  old_infeasible <- dom$infeasible_
  old_reason <- dom$infeasible_reason
  if (dom$infeasible_ && dom$infeasible_pos == n) {
    dom$infeasible_ <- FALSE
    dom$infeasible_reason <- .hmip_reason_unspecified()
  }
  while (k >= 0L) {
    .hmip_dom_undo(dom, k, TRUE)
    if (dom$domchgreason_type[k + 1L] == -1L) {
      dom$branchPos_ <- dom$branchPos_[-length(dom$branchPos_)]
      break
    }
    k <- k - 1L
  }
  if (old_infeasible) {
    .hmip_HighsDomain_markPropagateCut(dom, old_reason)
    dom$infeasible_reason <- .hmip_reason_unspecified()
    dom$infeasible_ <- FALSE
  }
  numreason <- length(dom$domchgreason_type)
  i <- k + 1L
  while (i < numreason) {
    .hmip_HighsDomain_markPropagateCut(dom, .hmip_reason(dom$domchgreason_type[i + 1L], dom$domchgreason_index[i + 1L]))
    i <- i + 1L
  }
  if (k < 0L) {
    dom$domchgstack_boundval <- numeric(); dom$domchgstack_column <- integer(); dom$domchgstack_boundtype <- integer()
    dom$prevboundval_val <- numeric(); dom$prevboundval_pos <- integer()
    dom$domchgreason_type <- integer(); dom$domchgreason_index <- integer()
    dom$branchPos_ <- integer()
    return(.hmip_domchg(0, -1L, 0L))
  }
  res <- .hmip_dom_stackchg(dom, k)
  keep <- seq_len(k)
  dom$domchgstack_boundval <- dom$domchgstack_boundval[keep]
  dom$domchgstack_column <- dom$domchgstack_column[keep]
  dom$domchgstack_boundtype <- dom$domchgstack_boundtype[keep]
  dom$domchgreason_type <- dom$domchgreason_type[keep]
  dom$domchgreason_index <- dom$domchgreason_index[keep]
  dom$prevboundval_val <- dom$prevboundval_val[keep]
  dom$prevboundval_pos <- dom$prevboundval_pos[keep]
  res
}

.hmip_dom_havePropagationRows <- function(dom) {
  if (length(dom$propagateinds_)) return(TRUE)
  if (!is.null(dom$objProp_) && .hmip_ObjProp_shouldBePropagated(dom$objProp_)) return(TRUE)
  for (cp in dom$cutpoolpropagation) if (length(cp$propagatecutinds_)) return(TRUE)
  for (cp in dom$conflictPoolPropagation) if (length(cp$propagateConflictInds_)) return(TRUE)
  FALSE
}

.hmip_dom_apply_changes <- function(dom, chg, reason) {
  for (j in seq_along(chg$column)) {
    if (dom$infeasible_) break
    .hmip_HighsDomain_changeBound(dom, .hmip_domchg(chg$boundval[j], chg$column[j], chg$boundtype[j]), reason)
  }
}

.hmip_HighsDomain_propagate <- function(dom) {
  if (!.hmip_dom_havePropagationRows(dom)) return(FALSE)
  ms <- dom$mipsolver; md <- ms$mipdata_; m <- ms$model_
  ft <- md$feastol
  while (.hmip_dom_havePropagationRows(dom)) {
    if (!is.null(dom$objProp_)) .hmip_ObjProp_propagate(dom$objProp_)

    for (cpp in dom$conflictPoolPropagation) {
      while (length(cpp$propagateConflictInds_)) {
        inds <- cpp$propagateConflictInds_
        cpp$propagateConflictInds_ <- integer()
        for (conflict in inds) .hmip_ConflictPoolPropagation_propagateConflict(cpp, conflict)
      }
    }

    if (length(dom$propagateinds_)) {
      inds <- dom$propagateinds_
      dom$propagateinds_ <- integer()
      dom$propagateflags_[inds + 1L] <- FALSE
      if (!dom$infeasible_) {
        numproprows <- length(inds)
        upchg <- vector("list", numproprows); lochg <- vector("list", numproprows)
        for (k in seq_len(numproprows)) {
          i <- inds[k]
          s <- md$ARstart_[i + 1L]; e <- md$ARstart_[i + 2L]
          idx <- if (e > s) (s + 1L):e else integer()
          Rindex <- md$ARindex_[idx]; Rvalue <- md$ARvalue_[idx]
          recomp <- FALSE
          ru <- m$row_upper_[i + 1L]; rl <- m$row_lower_[i + 1L]
          if (ru != Inf && (dom$activitymaxinf_[i + 1L] != 0L || .cd_val(.hmip_dom_amax(dom, i)) > ru + ft)) {
            .hmip_dom_set_amin(dom, i, .hmip_cd_renorm(.hmip_dom_amin(dom, i)))
            upchg[[k]] <- .hmip_HighsDomain_propagateRowUpper(dom, Rindex, Rvalue, length(idx), ru,
                                                             .hmip_dom_amin(dom, i), dom$activitymininf_[i + 1L])
            recomp <- TRUE
          }
          if (rl != -Inf && (dom$activitymininf_[i + 1L] != 0L || .cd_val(.hmip_dom_amin(dom, i)) < rl - ft)) {
            .hmip_dom_set_amax(dom, i, .hmip_cd_renorm(.hmip_dom_amax(dom, i)))
            lochg[[k]] <- .hmip_HighsDomain_propagateRowLower(dom, Rindex, Rvalue, length(idx), rl,
                                                             .hmip_dom_amax(dom, i), dom$activitymaxinf_[i + 1L])
            recomp <- TRUE
          }
          if (recomp) .hmip_dom_recomputeCapacityThreshold(dom, i)
        }
        for (k in seq_len(numproprows)) {
          i <- inds[k]
          if (length(upchg[[k]]$column)) {
            .hmip_dom_apply_changes(dom, upchg[[k]], .hmip_reason(-3L, i))
            if (dom$infeasible_) break
          }
          if (length(lochg[[k]]$column)) {
            .hmip_dom_apply_changes(dom, lochg[[k]], .hmip_reason(-4L, i))
            if (dom$infeasible_) break
          }
        }
      }
    }

    numpools <- length(dom$cutpoolpropagation)
    for (cutpool in seq_len(numpools) - 1L) {
      cpp <- dom$cutpoolpropagation[[cutpool + 1L]]
      if (!length(cpp$propagatecutinds_)) next
      inds <- cpp$propagatecutinds_
      cpp$propagatecutinds_ <- integer()
      cpp$propagatecutflags_[inds + 1L] <- bitwAnd(cpp$propagatecutflags_[inds + 1L], 2L)
      if (!dom$infeasible_) {
        numproprows <- length(inds)
        chg <- vector("list", numproprows)
        rhs <- .hmip_HighsCutPool_getRhs(cpp$cutpool)
        for (k in seq_len(numproprows)) {
          # the C++ tests the flag at the list position k, not at the cut
          if (bitwAnd(cpp$propagatecutflags_[k], 2L) != 0L) next
          i <- inds[k]
          cut <- .hmip_HighsCutPool_getCut(cpp$cutpool, i)
          a <- .hmip_cd_renorm(c(cpp$activitycuts_hi[i + 1L], cpp$activitycuts_lo[i + 1L]))
          cpp$activitycuts_hi[i + 1L] <- a[1L]; cpp$activitycuts_lo[i + 1L] <- a[2L]
          chg[[k]] <- .hmip_HighsDomain_propagateRowUpper(dom, cut$inds, cut$vals, cut$len, rhs[i + 1L], a,
                                                         cpp$activitycutsinf_[i + 1L])
          .hmip_CutpoolPropagation_recomputeCapacityThreshold(cpp, i)
        }
        for (k in seq_len(numproprows)) {
          i <- inds[k]
          if (length(chg[[k]]$column)) {
            .hmip_HighsCutPool_resetAge(cpp$cutpool, i)
            .hmip_dom_apply_changes(dom, chg[[k]], .hmip_reason(cutpool, i))
          }
          if (dom$infeasible_) break
        }
      }
    }
  }
  TRUE
}

# getColLowerPos / getColUpperPos: list(bound, pos)
.hmip_HighsDomain_getColLowerPos <- function(dom, col, stackpos) {
  lb <- dom$col_lower_[col + 1L]
  pos <- dom$colLowerPos_[col + 1L]
  while (pos > stackpos || (pos != -1L && dom$prevboundval_val[pos + 1L] == lb)) {
    lb <- dom$prevboundval_val[pos + 1L]
    pos <- dom$prevboundval_pos[pos + 1L]
  }
  list(bound = lb, pos = pos)
}
.hmip_HighsDomain_getColUpperPos <- function(dom, col, stackpos) {
  ub <- dom$col_upper_[col + 1L]
  pos <- dom$colUpperPos_[col + 1L]
  while (pos > stackpos || (pos != -1L && dom$prevboundval_val[pos + 1L] == ub)) {
    ub <- dom$prevboundval_val[pos + 1L]
    pos <- dom$prevboundval_pos[pos + 1L]
  }
  list(bound = ub, pos = pos)
}

.hmip_HighsDomain_tightenCoefficients <- function(dom, inds, vals, len, rhs) {
  maxactivity <- c(0, 0)
  for (i in seq_len(len)) {
    if (vals[i] > 0) {
      if (dom$col_upper_[inds[i] + 1L] == Inf) return(list(vals = vals, rhs = rhs))
      maxactivity <- .hmip_cd_iaddd(maxactivity, dom$col_upper_[inds[i] + 1L] * vals[i])
    } else {
      if (dom$col_lower_[inds[i] + 1L] == -Inf) return(list(vals = vals, rhs = rhs))
      maxactivity <- .hmip_cd_iaddd(maxactivity, dom$col_lower_[inds[i] + 1L] * vals[i])
    }
  }
  maxabscoef <- .cd_sub(maxactivity, rhs)
  if (.cd_val(maxabscoef) > dom$mipsolver$mipdata_$feastol) {
    upper <- c(rhs, 0)
    tightened <- 0L
    mac <- .cd_val(maxabscoef)
    for (i in seq_len(len)) {
      if (.hmip_dom_vtype(dom, inds[i]) == 0L) next
      if (vals[i] > mac) {
        delta <- .hmip_cd_rsub(vals[i], maxabscoef)
        upper <- .hmip_cd_isub(upper, .cd_mul(delta, dom$col_upper_[inds[i] + 1L]))
        vals[i] <- mac
        tightened <- tightened + 1L
      } else if (vals[i] < -mac) {
        delta <- .hmip_cd_rsub(-vals[i], maxabscoef)
        upper <- .cd_iadd(upper, .cd_mul(delta, dom$col_lower_[inds[i] + 1L]))
        vals[i] <- -mac
        tightened <- tightened + 1L
      }
    }
    if (tightened != 0L) rhs <- .cd_val(upper)
  }
  list(vals = vals, rhs = rhs)
}

.hmip_HighsDomain_getMinCutActivity <- function(dom, cutpool, cut) {
  for (cp in dom$cutpoolpropagation) {
    if (identical(cp$cutpool, cutpool)) {
      return(if (cut < length(cp$propagatecutflags_) && bitwAnd(cp$propagatecutflags_[cut + 1L], 2L) == 0L &&
                 cp$activitycutsinf_[cut + 1L] == 0L) cp$activitycuts_hi[cut + 1L] + cp$activitycuts_lo[cut + 1L] else -Inf)
    }
  }
  -Inf
}

.hmip_HighsDomain_setupObjectivePropagation <- function(dom) {
  dom$objProp_ <- .hmip_new_ObjectivePropagation(dom)
  invisible(NULL)
}
.hmip_HighsDomain_getObjectiveLowerBound <- function(dom) {
  op <- dom$objProp_
  if (!is.null(op) && op$numInfObjLower == 0L) return(.cd_val(op$objectiveLower))
  -Inf
}
.hmip_HighsDomain_getCutoffConstraint <- function(dom) {
  .hmip_ObjProp_getPropagationConstraint(dom$objProp_, length(dom$domchgstack_column))
}

# ---------------------------------------------------------------------------
# CutpoolPropagation (cut pool from group C)
# ---------------------------------------------------------------------------

.hmip_new_CutpoolPropagation <- function(cutpoolindex, domain, cutpool) {
  p <- new.env(parent = emptyenv())
  p$cutpoolindex <- cutpoolindex
  p$domain <- domain
  p$cutpool <- cutpool
  p$activitycuts_hi <- numeric(); p$activitycuts_lo <- numeric()
  p$activitycutsinf_ <- integer()
  p$propagatecutflags_ <- integer()
  p$propagatecutinds_ <- integer()
  p$capacityThreshold_ <- numeric()
  .hmip_HighsCutPool_addPropagationDomain(cutpool, p)
  p
}
.hmip_CutpoolPropagation_copy <- function(other, domain) {
  p <- .hmip_copy_env(other)
  p$domain <- domain
  .hmip_HighsCutPool_addPropagationDomain(p$cutpool, p)
  p
}

.hmip_CutpoolPropagation_recomputeCapacityThreshold <- function(p, cut) {
  dom <- p$domain
  mat <- .hmip_HighsCutPool_getMatrix(p$cutpool)
  start <- .hmip_HighsDynamicRowMatrix_getRowStart(mat, cut)
  end <- .hmip_HighsDynamicRowMatrix_getRowEnd(mat, cut)
  arindex <- .hmip_HighsDynamicRowMatrix_getARindex(mat)
  arvalue <- .hmip_HighsDynamicRowMatrix_getARvalue(mat)
  ft <- dom$mipsolver$mipdata_$feastol
  thr <- -ft
  i <- start
  while (i < end) {
    col <- arindex[i + 1L]
    if (dom$col_upper_[col + 1L] != dom$col_lower_[col + 1L]) {
      t <- abs(arvalue[i + 1L]) * .hmip_boundRange(dom$col_upper_[col + 1L], dom$col_lower_[col + 1L], ft, .hmip_dom_vtype(dom, col))
      thr <- max(thr, t, ft)
    }
    i <- i + 1L
  }
  p$capacityThreshold_[cut + 1L] <- thr
}

.hmip_CutpoolPropagation_cutAdded <- function(p, cut, propagate) {
  dom <- p$domain
  if (!propagate && !identical(dom, dom$mipsolver$mipdata_$domain)) return(invisible(NULL))
  mat <- .hmip_HighsCutPool_getMatrix(p$cutpool)
  start <- .hmip_HighsDynamicRowMatrix_getRowStart(mat, cut)
  end <- .hmip_HighsDynamicRowMatrix_getRowEnd(mat, cut)
  if (length(p$activitycutsinf_) <= cut) {
    n0 <- length(p$activitycutsinf_); add <- cut + 1L - n0
    p$activitycuts_hi <- c(p$activitycuts_hi, numeric(add)); p$activitycuts_lo <- c(p$activitycuts_lo, numeric(add))
    p$activitycutsinf_ <- c(p$activitycutsinf_, integer(add))
    p$propagatecutflags_ <- c(p$propagatecutflags_, rep(2L, add))
    p$capacityThreshold_ <- c(p$capacityThreshold_, numeric(add))
  }
  p$propagatecutflags_[cut + 1L] <- bitwAnd(p$propagatecutflags_[cut + 1L], bitwNot(2L))
  a <- .hmip_HighsDomain_computeMinActivity(dom, start, end, .hmip_HighsDynamicRowMatrix_getARindex(mat),
                                           .hmip_HighsDynamicRowMatrix_getARvalue(mat))
  p$activitycutsinf_[cut + 1L] <- a$ninf
  p$activitycuts_hi[cut + 1L] <- a$activity[1L]; p$activitycuts_lo[cut + 1L] <- a$activity[2L]
  if (propagate) {
    .hmip_CutpoolPropagation_recomputeCapacityThreshold(p, cut)
    .hmip_CutpoolPropagation_markPropagateCut(p, cut)
  }
  invisible(NULL)
}

.hmip_CutpoolPropagation_cutDeleted <- function(p, cut, deletedOnlyForPropagation = FALSE) {
  if (deletedOnlyForPropagation && identical(p$domain, p$domain$mipsolver$mipdata_$domain)) return(invisible(NULL))
  if (cut < length(p$propagatecutflags_)) p$propagatecutflags_[cut + 1L] <- bitwOr(p$propagatecutflags_[cut + 1L], 2L)
  invisible(NULL)
}

.hmip_CutpoolPropagation_markPropagateCut <- function(p, cut) {
  if (p$propagatecutflags_[cut + 1L] == 0L &&
      (p$activitycutsinf_[cut + 1L] == 1L ||
       (.hmip_HighsCutPool_getRhs(p$cutpool)[cut + 1L] - (p$activitycuts_hi[cut + 1L] + p$activitycuts_lo[cut + 1L]) <=
        p$capacityThreshold_[cut + 1L]))) {
    p$propagatecutinds_ <- c(p$propagatecutinds_, cut)
    p$propagatecutflags_[cut + 1L] <- bitwOr(p$propagatecutflags_[cut + 1L], 1L)
  }
  invisible(NULL)
}

# updateActivityLbChange (lb = TRUE) / UbChange (lb = FALSE)
.hmip_cutprop_update <- function(p, col, oldbound, newbound, lb) {
  dom <- p$domain
  mat <- .hmip_HighsCutPool_getMatrix(p$cutpool)
  fpos <- .hmip_HighsDynamicRowMatrix_forEachPositiveColumnEntry
  fneg <- .hmip_HighsDynamicRowMatrix_forEachNegativeColumnEntry
  thrside <- if (lb) fneg else fpos
  actside <- if (lb) fpos else fneg
  inf <- if (lb) -Inf else Inf
  if ((lb && newbound < oldbound) || (!lb && newbound > oldbound)) {
    thrside(mat, col, function(row, val) {
      p$capacityThreshold_[row + 1L] <- if (lb) .hmip_dom_thrLb(dom, col, newbound, val, p$capacityThreshold_[row + 1L])
        else .hmip_dom_thrUb(dom, col, newbound, val, p$capacityThreshold_[row + 1L])
      TRUE
    })
  }
  ft <- dom$mipsolver$mipdata_$feastol
  rhs <- .hmip_HighsCutPool_getRhs(p$cutpool)
  actside(mat, col, function(row, val) {
    cd <- .hmip_computeDelta(val, oldbound, newbound, inf)
    p$activitycutsinf_[row + 1L] <- p$activitycutsinf_[row + 1L] + cd$dinf
    a <- .cd_iadd(c(p$activitycuts_hi[row + 1L], p$activitycuts_lo[row + 1L]), cd$d)
    p$activitycuts_hi[row + 1L] <- a[1L]; p$activitycuts_lo[row + 1L] <- a[2L]
    if (.cd_val(cd$d) <= 0) {
      p$capacityThreshold_[row + 1L] <- if (lb) .hmip_dom_thrLb(dom, col, newbound, val, p$capacityThreshold_[row + 1L])
        else .hmip_dom_thrUb(dom, col, newbound, val, p$capacityThreshold_[row + 1L])
      return(TRUE)
    }
    if (p$activitycutsinf_[row + 1L] == 0L && .cd_val(.cd_sub(a, rhs[row + 1L])) > ft) {
      .hmip_dom_setinf(dom, p$cutpoolindex, row)
      return(FALSE)
    }
    .hmip_CutpoolPropagation_markPropagateCut(p, row)
    TRUE
  })
  if (dom$infeasible_) {
    t <- oldbound; oldbound <- newbound; newbound <- t
    actside(mat, col, function(row, val) {
      cd <- .hmip_computeDelta(val, oldbound, newbound, inf)
      p$activitycutsinf_[row + 1L] <- p$activitycutsinf_[row + 1L] + cd$dinf
      a <- .cd_iadd(c(p$activitycuts_hi[row + 1L], p$activitycuts_lo[row + 1L]), cd$d)
      p$activitycuts_hi[row + 1L] <- a[1L]; p$activitycuts_lo[row + 1L] <- a[2L]
      dom$infeasible_reason$index != row
    })
  }
  invisible(NULL)
}
.hmip_CutpoolPropagation_updateActivityLbChange <- function(p, col, oldbound, newbound) .hmip_cutprop_update(p, col, oldbound, newbound, TRUE)
.hmip_CutpoolPropagation_updateActivityUbChange <- function(p, col, oldbound, newbound) .hmip_cutprop_update(p, col, oldbound, newbound, FALSE)

# ---------------------------------------------------------------------------
# ConflictPoolPropagation
# ---------------------------------------------------------------------------

.hmip_new_ConflictPoolPropagation <- function(conflictpoolindex, domain, conflictpool) {
  p <- new.env(parent = emptyenv())
  ncol <- length(domain$mipsolver$model_$col_lower_)
  p$conflictpoolindex <- conflictpoolindex
  p$domain <- domain
  p$conflictpool_ <- conflictpool
  p$cutpool <- NULL
  p$colLowerWatched_ <- rep(-1L, ncol)
  p$colUpperWatched_ <- rep(-1L, ncol)
  p$conflictFlag_ <- integer()
  p$propagateConflictInds_ <- integer()
  p$wl_boundval <- numeric(); p$wl_column <- integer(); p$wl_boundtype <- integer()
  p$wl_prev <- integer(); p$wl_next <- integer()
  .hmip_HighsConflictPool_addPropagationDomain(conflictpool, p)
  p
}
.hmip_ConflictPoolPropagation_copy <- function(other, domain) {
  p <- .hmip_copy_env(other)
  p$domain <- domain
  .hmip_HighsConflictPool_addPropagationDomain(p$conflictpool_, p)
  p
}

.hmip_cpp_link <- function(p, linkPos) {
  col <- p$wl_column[linkPos + 1L]
  lower <- p$wl_boundtype[linkPos + 1L] == 0L
  head <- if (lower) p$colLowerWatched_[col + 1L] else p$colUpperWatched_[col + 1L]
  p$wl_prev[linkPos + 1L] <- -1L
  p$wl_next[linkPos + 1L] <- head
  # as in the C++ the head only moves when the list was not empty
  if (head != -1L) {
    p$wl_prev[head + 1L] <- linkPos
    if (lower) p$colLowerWatched_[col + 1L] <- linkPos else p$colUpperWatched_[col + 1L] <- linkPos
  }
}
.hmip_cpp_unlink <- function(p, linkPos) {
  col <- p$wl_column[linkPos + 1L]
  if (col == -1L) return()
  lower <- p$wl_boundtype[linkPos + 1L] == 0L
  p$wl_column[linkPos + 1L] <- -1L
  prev <- p$wl_prev[linkPos + 1L]; nxt <- p$wl_next[linkPos + 1L]
  if (prev != -1L) p$wl_next[prev + 1L] <- nxt
  else if (lower) p$colLowerWatched_[col + 1L] <- nxt else p$colUpperWatched_[col + 1L] <- nxt
  if (nxt != -1L) p$wl_prev[nxt + 1L] <- prev
}
.hmip_cpp_setwatch <- function(p, pos, dc) {
  p$wl_boundval[pos + 1L] <- dc$boundval; p$wl_column[pos + 1L] <- dc$column; p$wl_boundtype[pos + 1L] <- dc$boundtype
}
.hmip_cpp_watch <- function(p, pos) .hmip_domchg(p$wl_boundval[pos + 1L], p$wl_column[pos + 1L], p$wl_boundtype[pos + 1L])

.hmip_ConflictPoolPropagation_conflictDeleted <- function(p, conflict) {
  p$conflictFlag_[conflict + 1L] <- bitwOr(p$conflictFlag_[conflict + 1L], 8L)
  .hmip_cpp_unlink(p, 2L * conflict)
  .hmip_cpp_unlink(p, 2L * conflict + 1L)
  invisible(NULL)
}

.hmip_ConflictPoolPropagation_conflictAdded <- function(p, conflict) {
  pool <- p$conflictpool_
  dom <- p$domain
  start <- pool$range_first[conflict + 1L]; end <- pool$range_second[conflict + 1L]
  if (length(p$conflictFlag_) <= conflict) {
    nw <- 2L * conflict + 2L; add <- nw - length(p$wl_column)
    p$wl_boundval <- c(p$wl_boundval, numeric(add)); p$wl_column <- c(p$wl_column, rep(-1L, add))
    p$wl_boundtype <- c(p$wl_boundtype, integer(add)); p$wl_prev <- c(p$wl_prev, rep(-1L, add))
    p$wl_next <- c(p$wl_next, rep(-1L, add))
    p$conflictFlag_ <- c(p$conflictFlag_, integer(conflict + 1L - length(p$conflictFlag_)))
  }
  entry <- function(i) .hmip_domchg(pool$entry_boundval[i + 1L], pool$entry_column[i + 1L], pool$entry_boundtype[i + 1L])
  entpos <- function(i) {
    col <- pool$entry_column[i + 1L]
    if (pool$entry_boundtype[i + 1L] == 0L) dom$colLowerPos_[col + 1L] else dom$colUpperPos_[col + 1L]
  }
  numWatched <- 0L
  i <- start
  while (i != end) {
    e <- entry(i)
    if (!.hmip_HighsDomain_isActive(dom, e)) {
      wp <- 2L * conflict + numWatched
      .hmip_cpp_setwatch(p, wp, e)
      .hmip_cpp_link(p, wp)
      numWatched <- numWatched + 1L
      if (numWatched == 2L) break
    }
    i <- i + 1L
  }
  if (numWatched == 0L) {
    la_pos <- c(0L, 0L); la_i <- c(0L, 0L)
    numActive <- 0L
    i <- start
    while (i != end) {
      pos <- entpos(i)
      if (numActive == 0L) {
        la_pos[1L] <- pos; la_i[1L] <- i; numActive <- 1L
      } else if (numActive == 1L) {
        la_pos[2L] <- pos; la_i[2L] <- i; numActive <- 2L
        if (la_pos[1L] < la_pos[2L]) { la_pos <- rev(la_pos); la_i <- rev(la_i) }
      } else if (pos > la_pos[2L]) {
        la_pos[2L] <- pos; la_i[2L] <- i
        if (la_pos[1L] < la_pos[2L]) { la_pos <- rev(la_pos); la_i <- rev(la_i) }
      }
      i <- i + 1L
    }
    for (j in seq_len(numActive)) {
      wp <- 2L * conflict + j - 1L
      .hmip_cpp_setwatch(p, wp, entry(la_i[j]))
      .hmip_cpp_link(p, wp)
    }
  } else if (numWatched == 1L) {
    latestActive <- -1L; latestPos <- -1L
    i <- start
    while (i != end) {
      pos <- entpos(i)
      if (pos > latestPos) { latestActive <- i; latestPos <- pos }
      i <- i + 1L
    }
    if (latestActive != -1L) {
      wp <- 2L * conflict + 1L
      .hmip_cpp_setwatch(p, wp, entry(latestActive))
      .hmip_cpp_link(p, wp)
    }
  }
  p$conflictFlag_[conflict + 1L] <- bitwOr(numWatched, bitwAnd(p$conflictFlag_[conflict + 1L], 4L))
  .hmip_ConflictPoolPropagation_markPropagateConflict(p, conflict)
  invisible(NULL)
}

.hmip_ConflictPoolPropagation_markPropagateConflict <- function(p, conflict) {
  if (p$conflictFlag_[conflict + 1L] < 2L) {
    p$propagateConflictInds_ <- c(p$propagateConflictInds_, conflict)
    p$conflictFlag_[conflict + 1L] <- bitwOr(p$conflictFlag_[conflict + 1L], 4L)
  }
  invisible(NULL)
}

.hmip_cpp_update <- function(p, col, oldbound, newbound, lb) {
  i <- if (lb) p$colLowerWatched_[col + 1L] else p$colUpperWatched_[col + 1L]
  while (i != -1L) {
    conflict <- i %/% 2L
    bv <- p$wl_boundval[i + 1L]
    delta <- if (lb) (bv > newbound) - (bv > oldbound) else (bv < newbound) - (bv < oldbound)
    if (delta != 0L) {
      # uint8_t arithmetic
      p$conflictFlag_[conflict + 1L] <- (p$conflictFlag_[conflict + 1L] + delta) %% 256L
      .hmip_ConflictPoolPropagation_markPropagateConflict(p, conflict)
    }
    i <- p$wl_next[i + 1L]
  }
  invisible(NULL)
}
.hmip_ConflictPoolPropagation_updateActivityLbChange <- function(p, col, oldbound, newbound) .hmip_cpp_update(p, col, oldbound, newbound, TRUE)
.hmip_ConflictPoolPropagation_updateActivityUbChange <- function(p, col, oldbound, newbound) .hmip_cpp_update(p, col, oldbound, newbound, FALSE)

.hmip_ConflictPoolPropagation_propagateConflict <- function(p, conflict) {
  p$conflictFlag_[conflict + 1L] <- bitwAnd(p$conflictFlag_[conflict + 1L], 11L)
  if (p$conflictFlag_[conflict + 1L] >= 2L) return(invisible(NULL))
  dom <- p$domain
  if (dom$infeasible_) return(invisible(NULL))
  pool <- p$conflictpool_
  start <- pool$range_first[conflict + 1L]
  if (start == -1L) {
    .hmip_cpp_unlink(p, 2L * conflict)
    .hmip_cpp_unlink(p, 2L * conflict + 1L)
    return(invisible(NULL))
  }
  end <- pool$range_second[conflict + 1L]
  entry <- function(i) .hmip_domchg(pool$entry_boundval[i + 1L], pool$entry_column[i + 1L], pool$entry_boundtype[i + 1L])
  inactive <- integer()
  i <- start
  while (i != end) {
    if (!.hmip_HighsDomain_isActive(dom, entry(i))) {
      inactive <- c(inactive, i)
      if (length(inactive) == 2L) break
    }
    i <- i + 1L
  }
  numInactive <- length(inactive)
  p$conflictFlag_[conflict + 1L] <- numInactive
  reason <- .hmip_reason(length(dom$cutpoolpropagation) + p$conflictpoolindex, conflict)
  if (numInactive == 0L) {
    dom$infeasible_ <- TRUE
    dom$infeasible_reason <- reason
    dom$infeasible_pos <- length(dom$domchgstack_column)
    .hmip_HighsConflictPool_resetAge(pool, conflict)
  } else if (numInactive == 1L) {
    dc <- .hmip_HighsDomain_flip(dom, entry(inactive[1L]))
    if (!.hmip_HighsDomain_isActive(dom, dc)) {
      .hmip_HighsDomain_changeBound(dom, .hmip_HighsDomain_flip(dom, entry(inactive[1L])), reason)
      .hmip_HighsConflictPool_resetAge(pool, conflict)
    }
  } else {
    w0 <- 2L * conflict
    if (!.hmip_domchg_eq(.hmip_cpp_watch(p, w0), entry(inactive[1L]))) {
      .hmip_cpp_unlink(p, w0)
      .hmip_cpp_setwatch(p, w0, entry(inactive[1L]))
      .hmip_cpp_link(p, w0)
    }
    if (!.hmip_domchg_eq(.hmip_cpp_watch(p, w0 + 1L), entry(inactive[2L]))) {
      .hmip_cpp_unlink(p, w0 + 1L)
      .hmip_cpp_setwatch(p, w0 + 1L, entry(inactive[2L]))
      .hmip_cpp_link(p, w0 + 1L)
    }
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# HighsConflictPool
# ---------------------------------------------------------------------------

.hmip_new_HighsConflictPool <- function(agelim, softlimit) {
  cp <- new.env(parent = emptyenv())
  cp$agelim_ <- as.integer(agelim)
  cp$softlimit_ <- as.integer(softlimit)
  cp$ageDistribution_ <- integer(agelim + 1L)
  cp$ages_ <- integer()
  cp$modification_ <- integer()
  cp$entry_boundval <- numeric(); cp$entry_column <- integer(); cp$entry_boundtype <- integer()
  cp$range_first <- integer(); cp$range_second <- integer()
  # std::set<std::pair<HighsInt, HighsInt>> kept sorted
  cp$free_first <- integer(); cp$free_second <- integer()
  cp$deletedConflicts_ <- integer()
  cp$propagationDomains <- list()
  cp
}

.hmip_HighsConflictPool_addPropagationDomain <- function(pool, prop) {
  pool$propagationDomains[[length(pool$propagationDomains) + 1L]] <- prop
  invisible(NULL)
}
.hmip_HighsConflictPool_removePropagationDomain <- function(pool, prop) {
  for (k in rev(seq_along(pool$propagationDomains))) {
    if (identical(pool$propagationDomains[[k]], prop)) {
      pool$propagationDomains[[k]] <- NULL
      return(invisible(NULL))
    }
  }
  invisible(NULL)
}
.hmip_HighsConflictPool_resetAge <- function(pool, conflict) {
  a <- pool$ages_[conflict + 1L]
  if (a > 0L) {
    pool$ageDistribution_[a + 1L] <- pool$ageDistribution_[a + 1L] - 1L
    pool$ageDistribution_[1L] <- pool$ageDistribution_[1L] + 1L
    pool$ages_[conflict + 1L] <- 0L
  }
  invisible(NULL)
}
.hmip_HighsConflictPool_setAgeLimit <- function(pool, agelim) {
  pool$agelim_ <- as.integer(agelim)
  n <- agelim + 1L
  ad <- pool$ageDistribution_
  pool$ageDistribution_ <- if (length(ad) >= n) ad[seq_len(n)] else c(ad, integer(n - length(ad)))
  invisible(NULL)
}
.hmip_HighsConflictPool_getModificationCount <- function(pool, cut) pool$modification_[cut + 1L]
.hmip_HighsConflictPool_getConflictEntryVector <- function(pool) {
  list(boundval = pool$entry_boundval, column = pool$entry_column, boundtype = pool$entry_boundtype)
}
.hmip_HighsConflictPool_getConflictRanges <- function(pool) list(first = pool$range_first, second = pool$range_second)
.hmip_HighsConflictPool_getNumConflicts <- function(pool) length(pool$range_first) - length(pool$deletedConflicts_)

.hmip_cpool_free_insert <- function(pool, a, b) {
  f <- pool$free_first; s <- pool$free_second
  if (any(f == a & s == b)) return()
  k <- sum(f < a | (f == a & s < b))
  pool$free_first <- append(f, a, after = k)
  pool$free_second <- append(s, b, after = k)
}

# the common head of addConflictCut and addReconvergenceCut
.hmip_cpool_alloc <- function(pool, conflictLen) {
  f <- pool$free_first; s <- pool$free_second
  it <- if (length(f)) which(f > conflictLen | (f == conflictLen & s >= -1L))[1L] else NA
  if (is.na(it)) {
    start <- length(pool$entry_column)
    end <- start + conflictLen
    add <- conflictLen
    pool$entry_boundval <- c(pool$entry_boundval, numeric(add))
    pool$entry_column <- c(pool$entry_column, integer(add))
    pool$entry_boundtype <- c(pool$entry_boundtype, integer(add))
  } else {
    fs1 <- f[it]; fs2 <- s[it]
    pool$free_first <- f[-it]; pool$free_second <- s[-it]
    start <- fs2
    end <- start + conflictLen
    if (fs1 > conflictLen) .hmip_cpool_free_insert(pool, fs1 - conflictLen, end)
  }
  if (!length(pool$deletedConflicts_)) {
    conflictIndex <- length(pool$range_first)
    pool$range_first <- c(pool$range_first, start)
    pool$range_second <- c(pool$range_second, end)
    pool$ages_ <- c(pool$ages_, 0L)
    pool$modification_ <- c(pool$modification_, 0L)
  } else {
    nd <- length(pool$deletedConflicts_)
    conflictIndex <- pool$deletedConflicts_[nd]
    pool$deletedConflicts_ <- pool$deletedConflicts_[-nd]
    pool$range_first[conflictIndex + 1L] <- start
    pool$range_second[conflictIndex + 1L] <- end
  }
  pool$modification_[conflictIndex + 1L] <- pool$modification_[conflictIndex + 1L] + 1L
  pool$ages_[conflictIndex + 1L] <- 0L
  pool$ageDistribution_[1L] <- pool$ageDistribution_[1L] + 1L
  list(index = conflictIndex, start = start, end = end)
}

.hmip_cpool_fill <- function(pool, domain, i, frontier) {
  feastol <- domain$mipsolver$mipdata_$feastol
  for (k in seq_along(frontier$pos)) {
    bv <- frontier$bv[k]; col <- frontier$col[k]; bt <- frontier$bt[k]
    if (.hmip_dom_vtype(domain, col) == 0L) bv <- if (bt == 0L) bv + feastol else bv - feastol
    pool$entry_boundval[i + 1L] <- bv; pool$entry_column[i + 1L] <- col; pool$entry_boundtype[i + 1L] <- bt
    i <- i + 1L
  }
}

.hmip_HighsConflictPool_addConflictCut <- function(pool, domain, frontier) {
  a <- .hmip_cpool_alloc(pool, length(frontier$pos))
  .hmip_cpool_fill(pool, domain, a$start, frontier)
  for (prop in pool$propagationDomains) .hmip_ConflictPoolPropagation_conflictAdded(prop, a$index)
  invisible(NULL)
}

.hmip_HighsConflictPool_addReconvergenceCut <- function(pool, domain, frontier, domchg) {
  a <- .hmip_cpool_alloc(pool, length(frontier$pos) + 1L)
  f <- .hmip_HighsDomain_flip(domain, domchg)
  pool$entry_boundval[a$start + 1L] <- f$boundval
  pool$entry_column[a$start + 1L] <- f$column
  pool$entry_boundtype[a$start + 1L] <- f$boundtype
  .hmip_cpool_fill(pool, domain, a$start + 1L, frontier)
  for (prop in pool$propagationDomains) .hmip_ConflictPoolPropagation_conflictAdded(prop, a$index)
  invisible(NULL)
}

.hmip_HighsConflictPool_removeConflict <- function(pool, conflict) {
  for (prop in pool$propagationDomains) .hmip_ConflictPoolPropagation_conflictDeleted(prop, conflict)
  a <- pool$ages_[conflict + 1L]
  if (a >= 0L) {
    pool$ageDistribution_[a + 1L] <- pool$ageDistribution_[a + 1L] - 1L
    pool$ages_[conflict + 1L] <- -1L
  }
  start <- pool$range_first[conflict + 1L]; end <- pool$range_second[conflict + 1L]
  pool$deletedConflicts_ <- c(pool$deletedConflicts_, conflict)
  .hmip_cpool_free_insert(pool, end - start, start)
  pool$range_first[conflict + 1L] <- -1L
  pool$range_second[conflict + 1L] <- -1L
  pool$modification_[conflict + 1L] <- pool$modification_[conflict + 1L] + 1L
  invisible(NULL)
}

.hmip_HighsConflictPool_performAging <- function(pool) {
  conflictMaxIndex <- length(pool$range_first)
  agelim <- pool$agelim_
  numActive <- .hmip_HighsConflictPool_getNumConflicts(pool)
  while (agelim > 5L && numActive > pool$softlimit_) {
    numActive <- numActive - pool$ageDistribution_[agelim + 1L]
    agelim <- agelim - 1L
  }
  for (i in seq_len(conflictMaxIndex) - 1L) {
    a <- pool$ages_[i + 1L]
    if (a < 0L) next
    pool$ageDistribution_[a + 1L] <- pool$ageDistribution_[a + 1L] - 1L
    a <- a + 1L
    pool$ages_[i + 1L] <- a
    if (a > agelim) {
      pool$ages_[i + 1L] <- -1L
      .hmip_HighsConflictPool_removeConflict(pool, i)
    } else pool$ageDistribution_[a + 1L] <- pool$ageDistribution_[a + 1L] + 1L
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# ObjectivePropagation (HighsObjectiveFunction getters from group B1)
# ---------------------------------------------------------------------------

# ObjectiveContributionTree of a partition: the linked positions sorted by
# the key (-contribution, col)
.hmip_ot_less <- function(op, a, b) {
  ca <- -op$contrib[a + 1L]; cb <- -op$contrib[b + 1L]
  ca < cb || (ca == cb && op$ccol[a + 1L] < op$ccol[b + 1L])
}
.hmip_ot_link <- function(op, part, j) {
  v <- op$tree[[part + 1L]]
  k <- 0L
  while (k < length(v) && .hmip_ot_less(op, v[k + 1L], j)) k <- k + 1L
  op$tree[[part + 1L]] <- append(v, j, after = k)
}
.hmip_ot_unlink <- function(op, part, j) {
  v <- op$tree[[part + 1L]]
  op$tree[part + 1L] <- list(v[v != j])
}
.hmip_ot_first <- function(op, part) { v <- op$tree[[part + 1L]]; if (length(v)) v[1L] else -1L }
.hmip_ot_last <- function(op, part) { v <- op$tree[[part + 1L]]; if (length(v)) v[length(v)] else -1L }
.hmip_ot_successor <- function(op, part, j) {
  v <- op$tree[[part + 1L]]
  k <- match(j, v)
  if (k < length(v)) v[k + 1L] else -1L
}

.hmip_new_ObjectivePropagation <- function(domain) {
  op <- new.env(parent = emptyenv())
  op$domain <- domain
  ms <- domain$mipsolver
  op$objFunc <- ms$mipdata_$objectiveFunction
  op$cost <- ms$model_$col_cost_
  of <- op$objFunc
  objNonzeros <- .hmip_HighsObjectiveFunction_getObjectiveNonzeros(of)
  partitionStarts <- .hmip_HighsObjectiveFunction_getCliquePartitionStarts(of)
  numPartitions <- .hmip_HighsObjectiveFunction_getNumCliquePartitions(of)
  op$propagationConsBuffer <- numeric()
  op$pcd_multiplier <- numeric(numPartitions)
  op$pcd_rhs <- integer(numPartitions)
  if (numPartitions != 0L) op$propagationConsBuffer <- .hmip_HighsObjectiveFunction_getObjectiveValuesPacked(of)
  op$isPropagated <- FALSE
  op$capacityThreshold <- Inf
  op$objectiveLower <- c(0, 0)
  op$numInfObjLower <- 0L
  ncontrib <- partitionStarts[numPartitions + 1L]
  op$contrib <- numeric(ncontrib); op$ccol <- integer(ncontrib); op$cpart <- integer(ncontrib)
  op$tree <- vector("list", numPartitions)
  for (i in seq_len(numPartitions) - 1L) op$tree[i + 1L] <- list(integer())
  cost <- op$cost
  for (i in seq_len(numPartitions) - 1L) {
    op$pcd_rhs[i + 1L] <- 1L
    j <- partitionStarts[i + 1L]
    while (j < partitionStarts[i + 2L]) {
      col <- objNonzeros[j + 1L]
      op$ccol[j + 1L] <- col; op$cpart[j + 1L] <- i
      if (cost[col + 1L] > 0) {
        op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, cost[col + 1L])
        op$contrib[j + 1L] <- cost[col + 1L]
        op$pcd_rhs[i + 1L] <- op$pcd_rhs[i + 1L] - 1L
        if (domain$col_lower_[col + 1L] == 0) .hmip_ot_link(op, i, j)
      } else {
        op$contrib[j + 1L] <- -cost[col + 1L]
        if (domain$col_upper_[col + 1L] == 1) .hmip_ot_link(op, i, j)
      }
      j <- j + 1L
    }
    worst <- .hmip_ot_first(op, i)
    if (worst != -1L) op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, -op$contrib[worst + 1L])
  }
  i <- partitionStarts[numPartitions + 1L]
  while (i < length(objNonzeros)) {
    col <- objNonzeros[i + 1L]
    if (cost[col + 1L] > 0) {
      if (domain$col_lower_[col + 1L] == -Inf) op$numInfObjLower <- op$numInfObjLower + 1L
      else op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, domain$col_lower_[col + 1L] * cost[col + 1L])
    } else {
      if (domain$col_upper_[col + 1L] == Inf) op$numInfObjLower <- op$numInfObjLower + 1L
      else op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, domain$col_upper_[col + 1L] * cost[col + 1L])
    }
    i <- i + 1L
  }
  .hmip_ObjProp_recomputeCapacityThreshold(op)
  op
}

# getPropagationConstraint: list(vals, inds, len, rhs)
.hmip_ObjProp_getPropagationConstraint <- function(op, domchgStackSize, domchgCol = -1L) {
  of <- op$objFunc; dom <- op$domain
  numPartitions <- .hmip_HighsObjectiveFunction_getNumCliquePartitions(of)
  inds <- .hmip_HighsObjectiveFunction_getObjectiveNonzeros(of)
  len <- length(inds)
  if (numPartitions == 0L) {
    return(list(vals = .hmip_HighsObjectiveFunction_getObjectiveValuesPacked(of), inds = inds, len = len,
                rhs = dom$mipsolver$mipdata_$upper_limit))
  }
  partitionStarts <- .hmip_HighsObjectiveFunction_getCliquePartitionStarts(of)
  cost <- op$cost
  tmpRhs <- c(dom$mipsolver$mipdata_$upper_limit, 0)
  packed <- .hmip_HighsObjectiveFunction_getObjectiveValuesPacked(of)
  for (i in seq_len(numPartitions) - 1L) {
    start <- partitionStarts[i + 1L]; end <- partitionStarts[i + 2L]
    largest <- 0
    j <- start
    while (j < end) {
      c0 <- inds[j + 1L]
      if (c0 != domchgCol) {
        if (cost[c0 + 1L] > 0) {
          lb <- .hmip_HighsDomain_getColLowerPos(dom, c0, domchgStackSize)$bound
          if (lb < 1) largest <- max(largest, cost[c0 + 1L])
        } else {
          ub <- .hmip_HighsDomain_getColUpperPos(dom, c0, domchgStackSize)$bound
          if (ub > 0) largest <- max(largest, -cost[c0 + 1L])
        }
      }
      j <- j + 1L
    }
    tmpRhs <- .hmip_cd_iaddd(tmpRhs, largest * op$pcd_rhs[i + 1L])
    if (op$pcd_multiplier[i + 1L] != largest) {
      op$pcd_multiplier[i + 1L] <- largest
      idx <- if (end > start) (start + 1L):end else integer()
      op$propagationConsBuffer[idx] <- packed[idx] - ifelse(packed[idx] < 0 | (packed[idx] == 0 & 1 / packed[idx] < 0), -largest, largest)
    }
  }
  list(vals = op$propagationConsBuffer, inds = inds, len = len, rhs = .cd_val(tmpRhs))
}

.hmip_ObjProp_recomputeCapacityThreshold <- function(op) {
  of <- op$objFunc; dom <- op$domain
  partitionStarts <- .hmip_HighsObjectiveFunction_getCliquePartitionStarts(of)
  numPartitions <- .hmip_HighsObjectiveFunction_getNumCliquePartitions(of)
  ft <- dom$mipsolver$mipdata_$feastol
  thr <- -ft
  for (i in seq_len(numPartitions) - 1L) {
    worst <- .hmip_ot_first(op, i)
    if (worst == -1L) next
    if (.hmip_HighsDomain_isFixed(dom, op$ccol[worst + 1L])) next
    contribution <- op$contrib[worst + 1L]
    best <- .hmip_ot_last(op, i)
    if (best != worst) contribution <- contribution - op$contrib[best + 1L]
    thr <- max(thr, contribution * (1.0 - ft))
  }
  objNonzeros <- .hmip_HighsObjectiveFunction_getObjectiveNonzeros(of)
  i <- partitionStarts[numPartitions + 1L]
  while (i < length(objNonzeros)) {
    col <- objNonzeros[i + 1L]
    thr <- max(thr, abs(op$cost[col + 1L]) *
      .hmip_boundRange(dom$col_upper_[col + 1L], dom$col_lower_[col + 1L], ft, .hmip_dom_vtype(dom, col)))
    i <- i + 1L
  }
  op$capacityThreshold <- thr
}

# a clique partition column changes its contribution (relax: the bound moved
# away from the cost side and the column is linked again; tighten: unlinked)
.hmip_ObjProp_partition_relax <- function(op, partitionPos) {
  ft <- op$domain$mipsolver$mipdata_$feastol
  partition <- op$cpart[partitionPos + 1L]
  currFirst <- .hmip_ot_first(op, partition)
  .hmip_ot_link(op, partition, partitionPos)
  oldContribution <- if (currFirst != -1L) op$contrib[currFirst + 1L] else 0
  pc <- op$contrib[partitionPos + 1L]
  if (partitionPos == .hmip_ot_first(op, partition) && pc != oldContribution) {
    op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, oldContribution)
    op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, -pc)
    best <- .hmip_ot_last(op, partition)
    delta <- pc
    if (best != partitionPos) delta <- delta - op$contrib[best + 1L]
    op$capacityThreshold <- max(delta * (1.0 - ft), op$capacityThreshold)
  } else {
    op$capacityThreshold <- max((oldContribution - pc) * (1.0 - ft), op$capacityThreshold)
  }
}
.hmip_ObjProp_partition_tighten <- function(op, partitionPos) {
  partition <- op$cpart[partitionPos + 1L]
  wasFirst <- .hmip_ot_first(op, partition) == partitionPos
  if (wasFirst) op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, op$contrib[partitionPos + 1L])
  .hmip_ot_unlink(op, partition, partitionPos)
  if (wasFirst) {
    nw <- .hmip_ot_first(op, partition)
    if (nw != -1L) op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, -op$contrib[nw + 1L])
  }
}
.hmip_ObjProp_cutoff <- function(op) {
  dom <- op$domain
  if (op$numInfObjLower == 0L && .cd_val(op$objectiveLower) > dom$mipsolver$mipdata_$upper_limit) {
    .hmip_dom_setinf(dom, -7L, 0L)
    return(TRUE)
  }
  FALSE
}

.hmip_ObjProp_updateActivityLbChange <- function(op, col, oldbound, newbound) {
  dom <- op$domain; cost <- op$cost[col + 1L]
  ft <- dom$mipsolver$mipdata_$feastol
  if (cost <= 0) {
    if (cost != 0 && newbound < oldbound) {
      op$capacityThreshold <- max(op$capacityThreshold,
        -cost * .hmip_boundRange(dom$col_upper_[col + 1L], newbound, ft, .hmip_dom_vtype(dom, col)))
      op$isPropagated <- FALSE
    }
    return(invisible(NULL))
  }
  op$isPropagated <- FALSE
  partitionPos <- .hmip_HighsObjectiveFunction_getColCliquePartition(op$objFunc, col)
  if (partitionPos == -1L) {
    if (oldbound == -Inf) op$numInfObjLower <- op$numInfObjLower - 1L
    else op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, -(oldbound * cost))
    if (newbound == -Inf) op$numInfObjLower <- op$numInfObjLower + 1L
    else op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, newbound * cost)
    if (newbound < oldbound) {
      op$capacityThreshold <- max(op$capacityThreshold,
        cost * .hmip_boundRange(dom$col_upper_[col + 1L], dom$col_lower_[col + 1L], ft, .hmip_dom_vtype(dom, col)))
    } else if (.hmip_ObjProp_cutoff(op)) {
      .hmip_ObjProp_updateActivityLbChange(op, col, newbound, oldbound)
    }
  } else if (newbound == 0) {
    .hmip_ObjProp_partition_relax(op, partitionPos)
  } else {
    .hmip_ObjProp_partition_tighten(op, partitionPos)
    if (.hmip_ObjProp_cutoff(op)) .hmip_ObjProp_updateActivityLbChange(op, col, newbound, oldbound)
  }
  invisible(NULL)
}

.hmip_ObjProp_updateActivityUbChange <- function(op, col, oldbound, newbound) {
  dom <- op$domain; cost <- op$cost[col + 1L]
  ft <- dom$mipsolver$mipdata_$feastol
  if (cost >= 0) {
    if (cost != 0 && newbound > oldbound) {
      op$capacityThreshold <- max(op$capacityThreshold,
        cost * .hmip_boundRange(newbound, dom$col_lower_[col + 1L], ft, .hmip_dom_vtype(dom, col)))
      op$isPropagated <- FALSE
    }
    return(invisible(NULL))
  }
  op$isPropagated <- FALSE
  partitionPos <- .hmip_HighsObjectiveFunction_getColCliquePartition(op$objFunc, col)
  if (partitionPos == -1L) {
    if (oldbound == Inf) op$numInfObjLower <- op$numInfObjLower - 1L
    else op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, -(oldbound * cost))
    if (newbound == Inf) op$numInfObjLower <- op$numInfObjLower + 1L
    else op$objectiveLower <- .hmip_cd_iaddd(op$objectiveLower, newbound * cost)
    if (newbound > oldbound) {
      op$capacityThreshold <- max(op$capacityThreshold,
        -cost * .hmip_boundRange(dom$col_upper_[col + 1L], dom$col_lower_[col + 1L], ft, .hmip_dom_vtype(dom, col)))
    } else if (.hmip_ObjProp_cutoff(op)) {
      .hmip_ObjProp_updateActivityUbChange(op, col, newbound, oldbound)
    }
  } else if (newbound == 1) {
    .hmip_ObjProp_partition_relax(op, partitionPos)
  } else {
    .hmip_ObjProp_partition_tighten(op, partitionPos)
    if (.hmip_ObjProp_cutoff(op)) .hmip_ObjProp_updateActivityUbChange(op, col, newbound, oldbound)
  }
  invisible(NULL)
}

.hmip_ObjProp_shouldBePropagated <- function(op) {
  if (op$isPropagated) return(FALSE)
  if (op$numInfObjLower > 1L) return(FALSE)
  dom <- op$domain
  if (dom$infeasible_) return(FALSE)
  upperLimit <- dom$mipsolver$mipdata_$upper_limit
  if (upperLimit == Inf) return(FALSE)
  if (upperLimit - .cd_val(op$objectiveLower) > op$capacityThreshold) return(FALSE)
  TRUE
}

.hmip_ObjProp_propagate <- function(op) {
  if (!.hmip_ObjProp_shouldBePropagated(op)) return(invisible(NULL))
  dom <- op$domain; of <- op$objFunc; cost <- op$cost
  ft <- dom$mipsolver$mipdata_$feastol
  upperLimit <- dom$mipsolver$mipdata_$upper_limit
  if (op$numInfObjLower == 0L && .cd_val(op$objectiveLower) > upperLimit) {
    .hmip_dom_setinf(dom, -7L, 0L)
    return(invisible(NULL))
  }
  objNonzeros <- .hmip_HighsObjectiveFunction_getObjectiveNonzeros(of)
  robj <- .hmip_reason(-7L, 0L)
  capacity <- .hmip_cd_rsub(upperLimit, op$objectiveLower)
  if (op$numInfObjLower == 1L) {
    numCol <- length(objNonzeros)
    i <- .hmip_HighsObjectiveFunction_getNumBinariesInObjective(of)
    while (i < numCol) {
      col <- objNonzeros[i + 1L]
      i <- i + 1L
      if (cost[col + 1L] > 0) {
        if (dom$col_lower_[col + 1L] > -Inf) next
        boundVal <- .cd_div(capacity, cost[col + 1L])
        if (abs(.cd_val(boundVal) * .hmip_kHighsTiny) > ft) next
        a <- .hmip_HighsDomain_adjustedUb(dom, col, boundVal)
        if (a$accept) .hmip_HighsDomain_changeBound(dom, .hmip_domchg(a$bound, col, 1L), robj)
        break
      } else {
        if (dom$col_upper_[col + 1L] < Inf) next
        boundVal <- .cd_div(capacity, cost[col + 1L])
        if (abs(.cd_val(boundVal) * .hmip_kHighsTiny) > ft) next
        a <- .hmip_HighsDomain_adjustedLb(dom, col, boundVal)
        if (a$accept) .hmip_HighsDomain_changeBound(dom, .hmip_domchg(a$bound, col, 0L), robj)
        break
      }
    }
  } else {
    numPartitions <- .hmip_HighsObjectiveFunction_getNumCliquePartitions(of)
    partitionStarts <- .hmip_HighsObjectiveFunction_getCliquePartitionStarts(of)
    currLb <- .cd_val(op$objectiveLower)
    repeat {
      for (i in seq_len(numPartitions) - 1L) {
        worst <- .hmip_ot_first(op, i)
        if (worst == -1L) next
        contribution <- op$contrib[worst + 1L]
        secondWorst <- .hmip_ot_successor(op, i, worst)
        if (secondWorst != -1L) contribution <- contribution - op$contrib[secondWorst + 1L]
        if (contribution > .cd_val(capacity)) {
          col <- op$ccol[worst + 1L]
          if (cost[col + 1L] > 0) {
            if (dom$col_upper_[col + 1L] > 0) {
              .hmip_HighsDomain_changeBound(dom, .hmip_domchg(0, col, 1L), robj)
              if (dom$infeasible_) break
            }
          } else {
            if (dom$col_lower_[col + 1L] < 1) {
              .hmip_HighsDomain_changeBound(dom, .hmip_domchg(1, col, 0L), robj)
              if (dom$infeasible_) break
            }
          }
        } else if (secondWorst != -1L) {
          best <- .hmip_ot_last(op, i)
          while (best != .hmip_ot_first(op, i)) {
            if (op$contrib[.hmip_ot_first(op, i) + 1L] - op$contrib[best + 1L] > .cd_val(capacity)) {
              col <- op$ccol[best + 1L]
              if (cost[col + 1L] > 0) {
                .hmip_HighsDomain_changeBound(dom, .hmip_domchg(1, col, 0L), robj)
                if (dom$infeasible_) break
              } else {
                .hmip_HighsDomain_changeBound(dom, .hmip_domchg(0, col, 1L), robj)
                if (dom$infeasible_) break
              }
            } else break
            best <- .hmip_ot_last(op, i)
          }
          if (dom$infeasible_) break
        }
      }
      if (dom$infeasible_) break
      i <- partitionStarts[numPartitions + 1L]
      while (i < length(objNonzeros)) {
        col <- objNonzeros[i + 1L]
        i <- i + 1L
        if (cost[col + 1L] > 0) {
          boundVal <- .cd_div(.hmip_cd_iaddd(capacity, dom$col_lower_[col + 1L] * cost[col + 1L]), cost[col + 1L])
          if (abs(.cd_val(boundVal) * .hmip_kHighsTiny) > ft) next
          a <- .hmip_HighsDomain_adjustedUb(dom, col, boundVal)
          if (a$accept) {
            .hmip_HighsDomain_changeBound(dom, .hmip_domchg(a$bound, col, 1L), robj)
            if (dom$infeasible_) break
          }
        } else {
          boundVal <- .cd_div(.hmip_cd_iaddd(capacity, dom$col_upper_[col + 1L] * cost[col + 1L]), cost[col + 1L])
          if (abs(.cd_val(boundVal) * .hmip_kHighsTiny) > ft) next
          a <- .hmip_HighsDomain_adjustedLb(dom, col, boundVal)
          if (a$accept) {
            .hmip_HighsDomain_changeBound(dom, .hmip_domchg(a$bound, col, 0L), robj)
            if (dom$infeasible_) break
          }
        }
      }
      if (dom$infeasible_) break
      newLb <- .cd_val(op$objectiveLower)
      if (newLb == currLb) break
      currLb <- newLb
      capacity <- .hmip_cd_rsub(upperLimit, op$objectiveLower)
    }
  }
  .hmip_ObjProp_recomputeCapacityThreshold(op)
  op$isPropagated <- TRUE
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# ConflictSet (conflict analysis)
# ---------------------------------------------------------------------------

# std::set<LocalDomChg> ordered by pos: an env with parallel vectors
.hmip_new_frontier <- function() {
  f <- new.env(parent = emptyenv())
  f$pos <- integer(); f$col <- integer(); f$bt <- integer(); f$bv <- numeric()
  f
}
.hmip_fr_find <- function(f, pos) { k <- match(pos, f$pos); if (is.na(k)) 0L else k }
# insert(i): returns TRUE when inserted
.hmip_fr_insert <- function(f, pos, col, bt, bv) {
  if (!is.na(match(pos, f$pos))) return(FALSE)
  k <- sum(f$pos < pos)
  f$pos <- append(f$pos, pos, after = k); f$col <- append(f$col, col, after = k)
  f$bt <- append(f$bt, bt, after = k); f$bv <- append(f$bv, bv, after = k)
  TRUE
}
.hmip_fr_erase <- function(f, pos) {
  k <- match(pos, f$pos)
  if (is.na(k)) return()
  f$pos <- f$pos[-k]; f$col <- f$col[-k]; f$bt <- f$bt[-k]; f$bv <- f$bv[-k]
}
.hmip_fr_clear <- function(f) { f$pos <- integer(); f$col <- integer(); f$bt <- integer(); f$bv <- numeric() }

.hmip_new_ConflictSet <- function(localdom) {
  cs <- new.env(parent = emptyenv())
  cs$localdom <- localdom
  cs$globaldom <- localdom$mipsolver$mipdata_$domain
  cs$reasonSideFrontier <- .hmip_new_frontier()
  cs$reconvergenceFrontier <- .hmip_new_frontier()
  cs$resolveQueue <- integer()
  # resolvedDomainChanges as parallel vectors
  cs$rd_pos <- integer(); cs$rd_col <- integer(); cs$rd_bt <- integer(); cs$rd_bv <- numeric()
  cs$rb <- NULL
  cs
}
.hmip_cs_rd_clear <- function(cs) { cs$rd_pos <- integer(); cs$rd_col <- integer(); cs$rd_bt <- integer(); cs$rd_bv <- numeric() }
.hmip_cs_rd_push <- function(cs, pos, col, bt, bv) {
  cs$rd_pos <- c(cs$rd_pos, pos); cs$rd_col <- c(cs$rd_col, col); cs$rd_bt <- c(cs$rd_bt, bt); cs$rd_bv <- c(cs$rd_bv, bv)
}
.hmip_cs_rd_push_stack <- function(cs, pos) {
  d <- cs$localdom
  .hmip_cs_rd_push(cs, pos, d$domchgstack_column[pos + 1L], d$domchgstack_boundtype[pos + 1L], d$domchgstack_boundval[pos + 1L])
}

.hmip_computePrio <- function(val, bound, globalbound, numnodes) abs(val * (bound - globalbound) * (1 + numnodes))

# sort the resolve buffer by prio descending then boundPos ascending
.hmip_cs_sort_rb <- function(rb) {
  if (!length(rb$boundPos)) return(rb)
  o <- order(-rb$prio, rb$boundPos)
  lapply(rb, function(x) x[o])
}
.hmip_cs_empty_rb <- function() list(delta = numeric(), baseBound = numeric(), prio = numeric(), boundPos = integer(), valuePos = integer())

# the candidate loop of explainBoundChangeLeq (leq = TRUE) / Geq
.hmip_cs_explainBoundChange_lin <- function(cs, frontier, dpos, dcol, dbt, dbv, inds, vals, len, rhs, act, leq) {
  if ((leq && act == -Inf) || (!leq && act == Inf)) return(FALSE)
  ld <- cs$localdom; gd <- cs$globaldom
  nq <- ld$mipsolver$mipdata_$nodequeue
  domchgVal <- 0
  rb <- .hmip_cs_empty_rb()
  M <- c(act, 0)
  for (i in seq_len(len) - 1L) {
    col <- inds[i + 1L]; v <- vals[i + 1L]
    if (col == dcol) { domchgVal <- v; next }
    uselower <- (v > 0) == leq
    if (uselower) {
      r <- .hmip_HighsDomain_getColLowerPos(ld, col, dpos)
      glb <- gd$col_lower_[col + 1L]
      if (glb >= r$bound || r$pos == -1L) next
      k <- .hmip_fr_find(frontier, r$pos)
      if (k) {
        base <- frontier$bv[k]
        if (base != glb) M <- .hmip_cd_iaddd(M, v * (base - glb))
        if (base >= r$bound) next
      } else base <- glb
      delta <- v * (r$bound - base)
      prio <- .hmip_computePrio(v, r$bound, glb, .hmip_HighsNodeQueue_numNodesUp(nq, col))
    } else {
      r <- .hmip_HighsDomain_getColUpperPos(ld, col, dpos)
      gub <- gd$col_upper_[col + 1L]
      if (gub <= r$bound || r$pos == -1L) next
      k <- .hmip_fr_find(frontier, r$pos)
      if (k) {
        base <- frontier$bv[k]
        if (base != gub) M <- .hmip_cd_iaddd(M, v * (base - gub))
        if (base <= r$bound) next
      } else base <- gub
      delta <- v * (r$bound - base)
      prio <- .hmip_computePrio(v, r$bound, gub, .hmip_HighsNodeQueue_numNodesDown(nq, col))
    }
    rb$delta <- c(rb$delta, delta); rb$baseBound <- c(rb$baseBound, base); rb$prio <- c(rb$prio, prio)
    rb$boundPos <- c(rb$boundPos, r$pos); rb$valuePos <- c(rb$valuePos, i)
  }
  if (domchgVal == 0) return(FALSE)
  cs$rb <- .hmip_cs_sort_rb(rb)
  md <- ld$mipsolver$mipdata_
  b0 <- dbv
  if (.hmip_dom_vtype(ld, dcol) != 0L) {
    if (dbt == 0L) b0 <- b0 - (1.0 - 10 * md$feastol) else b0 <- b0 + (1.0 - 10 * md$feastol)
  } else {
    if (dbt == 0L) b0 <- b0 - md$epsilon else b0 <- b0 + md$epsilon
  }
  b0 <- b0 * domchgVal
  Mlim <- rhs - b0
  if (leq) {
    M <- .hmip_cd_iaddd(M, -(domchgVal * (if (domchgVal < 0) gd$col_upper_[dcol + 1L] else gd$col_lower_[dcol + 1L])))
    .hmip_cs_resolveLinear(cs, M, Mlim, vals, TRUE)
  } else {
    M <- .hmip_cd_iaddd(M, -(domchgVal * (if (domchgVal < 0) gd$col_lower_[dcol + 1L] else gd$col_upper_[dcol + 1L])))
    .hmip_cs_resolveLinear(cs, M, Mlim, vals, FALSE)
  }
}

# resolveLinearLeq (leq = TRUE, Mlim = Mlower) / resolveLinearGeq (Mupper)
.hmip_cs_resolveLinear <- function(cs, M, Mlim, vals, leq) {
  .hmip_cs_rd_clear(cs)
  ld <- cs$localdom; md <- ld$mipsolver$mipdata_
  ft <- md$feastol; eps <- md$epsilon
  rb <- cs$rb
  # covered is M - Mlim; Leq needs covered >= 0, Geq needs covered <= 0
  sgn <- if (leq) 1 else -1
  covered <- .cd_val(.cd_sub(M, Mlim))
  if (sgn * covered < 0) {
    for (k in seq_along(rb$boundPos)) {
      .hmip_cs_rd_push_stack(cs, rb$boundPos[k])
      M <- .hmip_cd_iaddd(M, rb$delta[k])
      covered <- .cd_val(.cd_sub(M, Mlim))
      if (sgn * covered >= 0) break
    }
    if (sgn * covered < 0) return(FALSE)
    if (sgn * covered > ft) {
      k <- length(cs$rd_pos) - 1L
      while (k >= 0L) {
        i <- rb$valuePos[k + 1L]
        col <- cs$rd_col[k + 1L]
        delta <- rb$delta[k + 1L]
        relax <- .cd_val(.cd_add(.cd_div(.hmip_cd_rsub(Mlim, .cd_sub(M, delta)), vals[i + 1L]), rb$baseBound[k + 1L]))
        if (cs$rd_bt[k + 1L] == 0L) {
          lb <- cs$rd_bv[k + 1L]; glb <- rb$baseBound[k + 1L]
          if (.hmip_dom_vtype(ld, col) != 0L) relax <- ceiling(relax)
          if (relax - lb >= -ft) { k <- k - 1L; next }
          cs$rd_bv[k + 1L] <- relax
          if (relax - glb <= eps) {
            last <- length(cs$rd_pos)
            for (nm in c("rd_pos", "rd_col", "rd_bt", "rd_bv")) {
              x <- get(nm, envir = cs); t <- x[last]; x[last] <- x[k + 1L]; x[k + 1L] <- t
              assign(nm, x[-last], envir = cs)
            }
            M <- .hmip_cd_iaddd(M, -delta)
          } else {
            while (relax <= ld$prevboundval_val[cs$rd_pos[k + 1L] + 1L]) cs$rd_pos[k + 1L] <- ld$prevboundval_pos[cs$rd_pos[k + 1L] + 1L]
            M <- .hmip_cd_iaddd(M, vals[i + 1L] * (relax - lb))
          }
        } else {
          ub <- cs$rd_bv[k + 1L]; gub <- rb$baseBound[k + 1L]
          if (.hmip_dom_vtype(ld, col) != 0L) relax <- floor(relax)
          if (relax - ub <= ft) { k <- k - 1L; next }
          cs$rd_bv[k + 1L] <- relax
          if (relax - gub >= -eps) {
            last <- length(cs$rd_pos)
            for (nm in c("rd_pos", "rd_col", "rd_bt", "rd_bv")) {
              x <- get(nm, envir = cs); t <- x[last]; x[last] <- x[k + 1L]; x[k + 1L] <- t
              assign(nm, x[-last], envir = cs)
            }
            M <- .hmip_cd_iaddd(M, -delta)
          } else {
            while (relax >= ld$prevboundval_val[cs$rd_pos[k + 1L] + 1L]) cs$rd_pos[k + 1L] <- ld$prevboundval_pos[cs$rd_pos[k + 1L] + 1L]
            M <- .hmip_cd_iaddd(M, vals[i + 1L] * (relax - ub))
          }
        }
        covered <- .cd_val(.cd_sub(M, Mlim))
        if (sgn * covered <= ft) break
        k <- k - 1L
      }
    }
  }
  TRUE
}

# explainInfeasibilityLeq (leq = TRUE) / Geq
.hmip_cs_explainInfeasibility_lin <- function(cs, inds, vals, len, rhs, act, leq) {
  if ((leq && act == -Inf) || (!leq && act == Inf)) return(FALSE)
  ld <- cs$localdom; gd <- cs$globaldom
  infeasible_pos <- if (ld$infeasible_) ld$infeasible_pos else .Machine$integer.max
  nq <- ld$mipsolver$mipdata_$nodequeue
  rb <- .hmip_cs_empty_rb()
  for (i in seq_len(len) - 1L) {
    col <- inds[i + 1L]; v <- vals[i + 1L]
    if ((v > 0) == leq) {
      r <- .hmip_HighsDomain_getColLowerPos(ld, col, infeasible_pos)
      base <- gd$col_lower_[col + 1L]
      if (base >= r$bound || r$pos == -1L) next
      prio <- .hmip_computePrio(v, r$bound, base, .hmip_HighsNodeQueue_numNodesUp(nq, col))
    } else {
      r <- .hmip_HighsDomain_getColUpperPos(ld, col, infeasible_pos)
      base <- gd$col_upper_[col + 1L]
      if (base <= r$bound || r$pos == -1L) next
      prio <- .hmip_computePrio(v, r$bound, base, .hmip_HighsNodeQueue_numNodesDown(nq, col))
    }
    rb$delta <- c(rb$delta, v * (r$bound - base)); rb$baseBound <- c(rb$baseBound, base); rb$prio <- c(rb$prio, prio)
    rb$boundPos <- c(rb$boundPos, r$pos); rb$valuePos <- c(rb$valuePos, i)
  }
  cs$rb <- .hmip_cs_sort_rb(rb)
  ft <- ld$mipsolver$mipdata_$feastol
  if (leq) .hmip_cs_resolveLinear(cs, c(act, 0), rhs + max(10.0, abs(rhs)) * ft, vals, TRUE)
  else .hmip_cs_resolveLinear(cs, c(act, 0), rhs - max(10.0, abs(rhs)) * ft, vals, FALSE)
}

.hmip_cs_getRow <- function(md, row) {
  s <- md$ARstart_[row + 1L]; e <- md$ARstart_[row + 2L]
  idx <- if (e > s) (s + 1L):e else integer()
  list(len = length(idx), inds = md$ARindex_[idx], vals = md$ARvalue_[idx])
}

.hmip_cs_conflict_entries <- function(pool, conflict) {
  s <- pool$range_first[conflict + 1L]; e <- pool$range_second[conflict + 1L]
  idx <- if (e > s) (s + 1L):e else integer()
  list(boundval = pool$entry_boundval[idx], column = pool$entry_column[idx], boundtype = pool$entry_boundtype[idx])
}

.hmip_cs_explainInfeasibility <- function(cs) {
  ld <- cs$localdom; gd <- cs$globaldom; ms <- ld$mipsolver; md <- ms$mipdata_
  t <- ld$infeasible_reason$type
  if (t == -2L || t == -1L) return(FALSE)
  if (t == -6L) {
    .hmip_cs_rd_clear(cs)
    p <- ld$infeasible_reason$index
    col <- ld$domchgstack_column[p + 1L]
    .hmip_cs_rd_push_stack(cs, p)
    other <- if (ld$domchgstack_boundtype[p + 1L] == 0L) .hmip_HighsDomain_getColUpperPos(ld, col, p)$pos
      else .hmip_HighsDomain_getColLowerPos(ld, col, p)$pos
    if (other != -1L) .hmip_cs_rd_push_stack(cs, other)
    return(TRUE)
  }
  if (t == -5L) return(FALSE)
  if (t == -4L) {
    r <- ld$infeasible_reason$index; row <- .hmip_cs_getRow(md, r)
    return(.hmip_cs_explainInfeasibility_lin(cs, row$inds, row$vals, row$len, ms$model_$row_lower_[r + 1L],
                                             .hmip_HighsDomain_getMaxActivity(gd, r), FALSE))
  }
  if (t == -3L) {
    r <- ld$infeasible_reason$index; row <- .hmip_cs_getRow(md, r)
    return(.hmip_cs_explainInfeasibility_lin(cs, row$inds, row$vals, row$len, ms$model_$row_upper_[r + 1L],
                                             .hmip_HighsDomain_getMinActivity(gd, r), TRUE))
  }
  if (t == -7L) {
    pc <- .hmip_ObjProp_getPropagationConstraint(ld$objProp_, ld$infeasible_pos)
    a <- .hmip_HighsDomain_computeMinActivity(gd, 0L, pc$len, pc$inds, pc$vals)
    if (a$ninf > 0L) return(FALSE)
    return(.hmip_cs_explainInfeasibility_lin(cs, pc$inds, pc$vals, pc$len, pc$rhs, .cd_val(a$activity), TRUE))
  }
  ncp <- length(ld$cutpoolpropagation)
  if (t < ncp) {
    cp <- ld$cutpoolpropagation[[t + 1L]]; cut <- ld$infeasible_reason$index
    cc <- .hmip_HighsCutPool_getCut(cp$cutpool, cut)
    minAct <- .hmip_HighsDomain_getMinCutActivity(gd, cp$cutpool, cut)
    return(.hmip_cs_explainInfeasibility_lin(cs, cc$inds, cc$vals, cc$len, .hmip_HighsCutPool_getRhs(cp$cutpool)[cut + 1L], minAct, TRUE))
  }
  cpp <- ld$conflictPoolPropagation[[t - ncp + 1L]]
  conflict <- ld$infeasible_reason$index
  if (bitwAnd(cpp$conflictFlag_[conflict + 1L], 8L) != 0L) return(FALSE)
  .hmip_cs_explainInfeasibilityConflict(cs, .hmip_cs_conflict_entries(cpp$conflictpool_, conflict))
}

.hmip_cs_explainInfeasibilityConflict <- function(cs, conflict) {
  .hmip_cs_rd_clear(cs)
  ld <- cs$localdom; gd <- cs$globaldom
  for (i in seq_along(conflict$column)) {
    e <- .hmip_domchg(conflict$boundval[i], conflict$column[i], conflict$boundtype[i])
    if (.hmip_HighsDomain_isActive(gd, e)) next
    if (e$boundtype == 0L) {
      r <- .hmip_HighsDomain_getColLowerPos(ld, e$column, ld$infeasible_pos)
      if (r$pos == -1L || r$bound < e$boundval) return(FALSE)
      pos <- r$pos
      while (ld$prevboundval_val[pos + 1L] >= e$boundval) pos <- ld$prevboundval_pos[pos + 1L]
    } else {
      r <- .hmip_HighsDomain_getColUpperPos(ld, e$column, ld$infeasible_pos)
      if (r$pos == -1L || r$bound > e$boundval) return(FALSE)
      pos <- r$pos
      while (ld$prevboundval_val[pos + 1L] <= e$boundval) pos <- ld$prevboundval_pos[pos + 1L]
    }
    .hmip_cs_rd_push(cs, pos, e$column, e$boundtype, e$boundval)
  }
  TRUE
}

.hmip_cs_explainBoundChange <- function(cs, frontier, dpos, dcol, dbt, dbv) {
  ld <- cs$localdom; gd <- cs$globaldom; ms <- ld$mipsolver; md <- ms$mipdata_
  t <- ld$domchgreason_type[dpos + 1L]; idx <- ld$domchgreason_index[dpos + 1L]
  if (t == -2L || t == -1L || t == -6L) return(FALSE)
  if (t == -5L) {
    col <- idx %/% 2L; val <- idx %% 2L
    .hmip_cs_rd_clear(cs)
    bp <- if (val) .hmip_HighsDomain_getColLowerPos(ld, col, dpos)$pos else .hmip_HighsDomain_getColUpperPos(ld, col, dpos)$pos
    if (bp != -1L) .hmip_cs_rd_push_stack(cs, bp)
    return(TRUE)
  }
  if (t == -4L) {
    row <- .hmip_cs_getRow(md, idx)
    return(.hmip_cs_explainBoundChange_lin(cs, frontier, dpos, dcol, dbt, dbv, row$inds, row$vals, row$len,
                                           ms$model_$row_lower_[idx + 1L], .hmip_HighsDomain_getMaxActivity(gd, idx), FALSE))
  }
  if (t == -3L) {
    row <- .hmip_cs_getRow(md, idx)
    return(.hmip_cs_explainBoundChange_lin(cs, frontier, dpos, dcol, dbt, dbv, row$inds, row$vals, row$len,
                                           ms$model_$row_upper_[idx + 1L], .hmip_HighsDomain_getMinActivity(gd, idx), TRUE))
  }
  if (t == -7L) {
    pc <- .hmip_ObjProp_getPropagationConstraint(ld$objProp_, dpos, dcol)
    a <- .hmip_HighsDomain_computeMinActivity(gd, 0L, pc$len, pc$inds, pc$vals)
    if (a$ninf == 1L) return(FALSE)
    return(.hmip_cs_explainBoundChange_lin(cs, frontier, dpos, dcol, dbt, dbv, pc$inds, pc$vals, pc$len, pc$rhs,
                                           .cd_val(a$activity), TRUE))
  }
  ncp <- length(ld$cutpoolpropagation)
  if (t < ncp) {
    cp <- ld$cutpoolpropagation[[t + 1L]]
    cc <- .hmip_HighsCutPool_getCut(cp$cutpool, idx)
    minAct <- .hmip_HighsDomain_getMinCutActivity(gd, cp$cutpool, idx)
    return(.hmip_cs_explainBoundChange_lin(cs, frontier, dpos, dcol, dbt, dbv, cc$inds, cc$vals, cc$len,
                                           .hmip_HighsCutPool_getRhs(cp$cutpool)[idx + 1L], minAct, TRUE))
  }
  cpp <- ld$conflictPoolPropagation[[t - ncp + 1L]]
  if (bitwAnd(cpp$conflictFlag_[idx + 1L], 8L) != 0L) return(FALSE)
  .hmip_cs_explainBoundChangeConflict(cs, dpos, .hmip_domchg(dbv, dcol, dbt), .hmip_cs_conflict_entries(cpp$conflictpool_, idx))
}

.hmip_cs_explainBoundChangeConflict <- function(cs, lpos, ldomchg, conflict) {
  .hmip_cs_rd_clear(cs)
  ld <- cs$localdom; gd <- cs$globaldom
  domchg <- .hmip_HighsDomain_flip(ld, ldomchg)
  found <- FALSE
  for (i in seq_along(conflict$column)) {
    e <- .hmip_domchg(conflict$boundval[i], conflict$column[i], conflict$boundtype[i])
    if (!found && e$column == domchg$column && e$boundtype == domchg$boundtype) {
      if (e$boundtype == 0L) {
        if (e$boundval <= domchg$boundval) { found <- TRUE; next }
      } else {
        if (e$boundval >= domchg$boundval) { found <- TRUE; next }
      }
    }
    if (.hmip_HighsDomain_isActive(gd, e)) next
    if (e$boundtype == 0L) {
      r <- .hmip_HighsDomain_getColLowerPos(ld, e$column, lpos - 1L)
      if (r$pos == -1L || r$bound < e$boundval) return(FALSE)
      pos <- r$pos
      while (ld$prevboundval_val[pos + 1L] >= e$boundval) pos <- ld$prevboundval_pos[pos + 1L]
    } else {
      r <- .hmip_HighsDomain_getColUpperPos(ld, e$column, lpos - 1L)
      if (r$pos == -1L || r$bound > e$boundval) return(FALSE)
      pos <- r$pos
      while (ld$prevboundval_val[pos + 1L] <= e$boundval) pos <- ld$prevboundval_pos[pos + 1L]
    }
    .hmip_cs_rd_push_stack(cs, pos)
  }
  found
}

.hmip_cs_resolvable <- function(cs, pos) {
  t <- cs$localdom$domchgreason_type[pos + 1L]
  !(t == -1L || t == -2L)
}

# resolveQueue holds distinct positions; pop_heap returns the largest
.hmip_cs_popQueue <- function(cs) {
  q <- cs$resolveQueue
  k <- which.max(q)
  cs$resolveQueue <- q[-k]
  q[k]
}

.hmip_cs_resolveDepth <- function(cs, frontier, depthLevel, stopSize, minResolve = 0L, increaseConflictScore = FALSE) {
  ld <- cs$localdom
  cs$resolveQueue <- integer()
  bp <- ld$branchPos_
  startPos <- if (depthLevel == 0L) 0L else bp[depthLevel] + 1L
  while (depthLevel < length(bp)) {
    b <- bp[depthLevel + 1L]
    if (ld$domchgstack_boundval[b + 1L] != ld$prevboundval_val[b + 1L]) break
    depthLevel <- depthLevel + 1L
  }
  endPos <- if (depthLevel == length(bp)) Inf else bp[depthLevel + 1L]
  sel <- frontier$pos[frontier$pos >= startPos & frontier$pos <= endPos]
  if (!length(sel)) return(-1L)
  for (p in sel) if (.hmip_cs_resolvable(cs, p)) cs$resolveQueue <- c(cs$resolveQueue, p)
  numResolved <- 0L
  pc <- ld$mipsolver$mipdata_$pseudocost
  while (length(cs$resolveQueue) > stopSize || (length(cs$resolveQueue) > 0L && numResolved < minResolve)) {
    pos <- .hmip_cs_popQueue(cs)
    k <- .hmip_fr_find(frontier, pos)
    if (!.hmip_cs_explainBoundChange(cs, frontier, pos, frontier$col[k], frontier$bt[k], frontier$bv[k])) next
    numResolved <- numResolved + 1L
    .hmip_fr_erase(frontier, pos)
    for (j in seq_along(cs$rd_pos)) {
      ip <- cs$rd_pos[j]
      if (.hmip_fr_insert(frontier, ip, cs$rd_col[j], cs$rd_bt[j], cs$rd_bv[j])) {
        if (increaseConflictScore) {
          scol <- ld$domchgstack_column[ip + 1L]
          if (ld$domchgstack_boundtype[ip + 1L] == 0L) .hmip_HighsPseudocost_increaseConflictScoreUp(pc, scol)
          else .hmip_HighsPseudocost_increaseConflictScoreDown(pc, scol)
        }
        if (ip >= startPos && .hmip_cs_resolvable(cs, ip)) cs$resolveQueue <- c(cs$resolveQueue, ip)
      } else {
        kk <- .hmip_fr_find(frontier, ip)
        if (cs$rd_bt[j] == 0L) frontier$bv[kk] <- max(frontier$bv[kk], cs$rd_bv[j])
        else frontier$bv[kk] <- min(frontier$bv[kk], cs$rd_bv[j])
      }
    }
  }
  numResolved
}

.hmip_cs_computeCuts <- function(cs, depthLevel, conflictPool) {
  ld <- cs$localdom
  numResolved <- .hmip_cs_resolveDepth(cs, cs$reasonSideFrontier, depthLevel, 1L,
                                       if (depthLevel == length(ld$branchPos_)) 1L else 0L, TRUE)
  if (numResolved == -1L) return(-1L)
  numConflicts <- 0L
  if (numResolved > 0L) {
    .hmip_HighsConflictPool_addConflictCut(conflictPool, ld, cs$reasonSideFrontier)
    numConflicts <- numConflicts + 1L
  }
  if (length(cs$resolveQueue) == 1L) {
    f <- cs$reasonSideFrontier
    up <- .hmip_cs_popQueue(cs)
    k <- .hmip_fr_find(f, up)
    uip <- list(pos = up, col = f$col[k], bt = f$bt[k], bv = f$bv[k])
    cs$resolveQueue <- integer()
    rf <- cs$reconvergenceFrontier
    .hmip_fr_clear(rf)
    .hmip_fr_insert(rf, uip$pos, uip$col, uip$bt, uip$bv)
    nr <- .hmip_cs_resolveDepth(cs, rf, depthLevel, 0L)
    if (nr > 0L && !.hmip_fr_find(rf, uip$pos)) {
      .hmip_HighsConflictPool_addReconvergenceCut(conflictPool, ld, rf, .hmip_domchg(uip$bv, uip$col, uip$bt))
      numConflicts <- numConflicts + 1L
    }
  }
  numConflicts
}

# the depth loop shared by both conflictAnalysis entry points
.hmip_cs_depthLoop <- function(cs, conflictPool) {
  ld <- cs$localdom
  numConflicts <- 0L
  lastDepth <- length(ld$branchPos_)
  currDepth <- lastDepth
  while (currDepth >= 0L) {
    if (currDepth > 0L) {
      b <- ld$branchPos_[currDepth]
      if (ld$domchgstack_boundval[b + 1L] == ld$prevboundval_val[b + 1L]) {
        lastDepth <- lastDepth - 1L
        currDepth <- currDepth - 1L
        next
      }
    }
    nn <- .hmip_cs_computeCuts(cs, currDepth, conflictPool)
    if (nn == -1L) {
      lastDepth <- lastDepth - 1L
      currDepth <- currDepth - 1L
      next
    }
    numConflicts <- numConflicts + nn
    if (numConflicts == 0L) break
    if (lastDepth - currDepth >= 4L && nn == 0L) break
    currDepth <- currDepth - 1L
  }
  if (currDepth == lastDepth) .hmip_HighsConflictPool_addConflictCut(conflictPool, ld, cs$reasonSideFrontier)
}

.hmip_cs_after_explain <- function(cs, conflictPool) {
  ld <- cs$localdom; md <- ld$mipsolver$mipdata_
  pc <- md$pseudocost
  .hmip_HighsPseudocost_increaseConflictWeight(pc)
  for (j in seq_along(cs$rd_pos)) {
    if (cs$rd_bt[j] == 0L) .hmip_HighsPseudocost_increaseConflictScoreUp(pc, cs$rd_col[j])
    else .hmip_HighsPseudocost_increaseConflictScoreDown(pc, cs$rd_col[j])
  }
  if (10 * length(cs$rd_pos) > 1000 + 3 * length(md$integral_cols)) return(invisible(NULL))
  for (j in seq_along(cs$rd_pos)) .hmip_fr_insert(cs$reasonSideFrontier, cs$rd_pos[j], cs$rd_col[j], cs$rd_bt[j], cs$rd_bv[j])
  .hmip_cs_depthLoop(cs, conflictPool)
  invisible(NULL)
}

.hmip_HighsDomain_conflictAnalysis <- function(dom, conflictPool) {
  gd <- dom$mipsolver$mipdata_$domain
  if (identical(gd, dom)) return(invisible(NULL))
  if (gd$infeasible_ || !dom$infeasible_) return(invisible(NULL))
  .hmip_HighsDomain_propagate(gd)
  if (gd$infeasible_) return(invisible(NULL))
  cs <- .hmip_new_ConflictSet(dom)
  if (!.hmip_cs_explainInfeasibility(cs)) return(invisible(NULL))
  .hmip_cs_after_explain(cs, conflictPool)
}

.hmip_HighsDomain_conflictAnalysisProof <- function(dom, proofinds, proofvals, prooflen, proofrhs, conflictPool) {
  gd <- dom$mipsolver$mipdata_$domain
  if (identical(gd, dom)) return(invisible(NULL))
  if (gd$infeasible_) return(invisible(NULL))
  .hmip_HighsDomain_propagate(gd)
  if (gd$infeasible_) return(invisible(NULL))
  cs <- .hmip_new_ConflictSet(dom)
  a <- .hmip_HighsDomain_computeMinActivity(gd, 0L, prooflen, proofinds, proofvals)
  if (a$ninf != 0L) return(invisible(NULL))
  if (!.hmip_cs_explainInfeasibility_lin(cs, proofinds, proofvals, prooflen, proofrhs, .cd_val(a$activity), TRUE))
    return(invisible(NULL))
  .hmip_cs_after_explain(cs, conflictPool)
}

.hmip_HighsDomain_conflictAnalyzeReconvergence <- function(dom, domchg, proofinds, proofvals, prooflen, proofrhs, conflictPool) {
  gd <- dom$mipsolver$mipdata_$domain
  if (identical(gd, dom)) return(invisible(NULL))
  if (gd$infeasible_) return(invisible(NULL))
  .hmip_HighsDomain_propagate(gd)
  if (gd$infeasible_) return(invisible(NULL))
  cs <- .hmip_new_ConflictSet(dom)
  a <- .hmip_HighsDomain_computeMinActivity(gd, 0L, prooflen, proofinds, proofvals)
  if (a$ninf != 0L) return(invisible(NULL))
  if (!.hmip_cs_explainBoundChange_lin(cs, cs$reconvergenceFrontier, length(dom$domchgstack_column), domchg$column,
                                       domchg$boundtype, domchg$boundval, proofinds, proofvals, prooflen, proofrhs,
                                       .cd_val(a$activity), TRUE)) return(invisible(NULL))
  if (10 * length(cs$rd_pos) > 1000 + 3 * length(dom$mipsolver$mipdata_$integral_cols)) return(invisible(NULL))
  for (j in seq_along(cs$rd_pos)) .hmip_fr_insert(cs$reconvergenceFrontier, cs$rd_pos[j], cs$rd_col[j], cs$rd_bt[j], cs$rd_bv[j])
  depth <- length(dom$branchPos_)
  while (depth > 0L) {
    b <- dom$branchPos_[depth]
    if (dom$domchgstack_boundval[b + 1L] != dom$prevboundval_val[b + 1L]) break
    depth <- depth - 1L
  }
  .hmip_cs_resolveDepth(cs, cs$reconvergenceFrontier, depth, 0L)
  .hmip_HighsConflictPool_addReconvergenceCut(conflictPool, dom, cs$reconvergenceFrontier, domchg)
  invisible(NULL)
}
