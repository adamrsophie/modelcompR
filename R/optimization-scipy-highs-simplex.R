# The simplex solver of HiGHS 1.12.0 (simplex/HEkk*.cpp, HApp.h and the
# scaling of lp_data/HighsLpUtils.cpp) as SciPy 1.18.1 vendors it behind
# linprog(method = "highs-ds") and the node LPs of milp. HiGHS is MIT
# licensed; the notice is in R/optimization-scipy-highs-util.R.
#
# The dual simplex method runs in the scaled space of HiGHS with the slack
# basis, the random vectors and cost perturbation of HEkk, dual steepest edge
# pricing (HEkkDualRHS), the bound flipping ratio test (HEkkDualRow), dual
# phase 1 on the artificial box of HEkk::initialiseBound and the cleanup by
# the primal simplex method (HEkkPrimal: Devex pricing, bound perturbation
# and shifting, the two pass Harris ratio test). The basis matrix is
# factored and updated by the port of HFactor
# (R/optimization-scipy-highs-factor.R) with the expected densities HEkk
# tracks.
#
# State lives in an environment `ek`; variables are 1..n (columns) and
# n + 1..n + m (the logicals of the rows, valued -A x as in HiGHS).

.hs_inf <- Inf

.hs_tight <- function(x) {
  x[abs(x) < 1e-14] <- 0
  x
}

.hs_colj <- function(ek, v) {
  if (v <= ek$n) return(ek$A[, v])
  e <- numeric(ek$m)
  e[v - ek$n] <- 1
  e
}

# HEkk::computeFactor: HFactor::build orders the basic variables by pivot
# row. A rank deficient basis is repaired as
# HEkk::initialiseSimplexLpBasisAndFactor does with buildMarkSingC and
# HEkk::handleRankDeficiency; with repair = FALSE (the refactor of a
# rebuild) the basis is left alone and the rank deficiency returned for
# HEkk::getNonsingularInverse to backtrack.
.hs_invert <- function(ek, repair = TRUE) {
  n <- ek$n
  cols <- lapply(ek$basicIndex, function(v) {
    if (v > n) return(list(index = v - n, value = 1, logical = TRUE))
    col <- ek$A[, v]; nz <- which(col != 0)
    list(index = nz, value = col[nz])
  })
  # computeFactor clears the bad basis changes
  ek$bad <- list()
  f <- .hs_factor_build(ek$m, cols)
  if (!repair && f$rank_deficiency > 0L) {
    ek$update_count <- 0L
    return(f$rank_deficiency)
  }
  repaired <- NULL
  if (f$rank_deficiency > 0L) {
    k <- seq_len(f$rank_deficiency)
    repaired <- list(row = f$row_with_no_pivot[k], vout = ek$basicIndex[f$col_with_no_pivot[k]])
    ek$basicIndex[f$col_with_no_pivot[k]] <- n + repaired$row
  }
  bi <- integer(ek$m)
  bi[f$permute] <- ek$basicIndex
  ek$basicIndex <- bi
  if (!is.null(repaired)) {
    ek$nonbasicFlag[n + repaired$row] <- 0L
    ek$nonbasicFlag[repaired$vout] <- 1L
    b <- .hs_lp_bounds(ek)
    lo <- b$lower; up <- b$upper
    ek$nonbasicMove <- ifelse(ek$nonbasicFlag == 0L | lo == up, 0L,
                              ifelse(is.finite(lo), ifelse(is.finite(up), ifelse(abs(lo) < abs(up), 1L, -1L), 1L),
                                     ifelse(is.finite(up), -1L, 0L)))
    # addBadBasisChange(row_out, variable_in, variable_out) as handleRankDeficiency calls it
    for (k in seq_along(repaired$row))
      .hs_add_bad_basis_change(ek, repaired$row[k], n + repaired$row[k], repaired$vout[k], "singular", TRUE)
  }
  ek$fa <- f
  ek$update_count <- 0L
  0L
}

# HEkk::initialiseControl: the running averages of the result densities
.hs_init_densities <- function(ek) {
  ek$dens <- list(col_aq = 0, row_ep = 0, row_ap = 0, row_DSE = 0, col_basic_feasibility_change = 0,
                  row_basic_feasibility_change = 0, col_BFRT = 0, primal_col = 0, dual_col = 1)
  ek$allow_dse_to_devex <- !identical(ek$dual_edge_weight_strategy, "steepest")
  ek$control_iteration_count0 <- ek$iteration_count %||% 0L
  ek$costly_dse_frequency <- 0
  ek$num_costly_dse_iteration <- 0L
  ek$costly_dse_measure <- 0
  ek$average_log_low_dse_weight_error <- 0
  ek$average_log_high_dse_weight_error <- 0
}

# HEkk::assessDSEWeightError
.hs_assess_dse_weight_error <- function(ek, computed, updated) {
  if (updated < computed) {
    ek$average_log_low_dse_weight_error <- 0.99 * ek$average_log_low_dse_weight_error + 0.01 * log(computed / updated)
  } else {
    ek$average_log_high_dse_weight_error <- 0.99 * ek$average_log_high_dse_weight_error + 0.01 * log(updated / computed)
  }
}

# HEkk::switchToDevex (dual_steepest_edge_weight_log_error_threshold 10)
.hs_switch_to_devex <- function(ek) {
  d <- ek$dens
  den <- max(max(d$row_ep, d$col_aq), d$row_ap)
  if (den > 0) {
    ek$costly_dse_measure <- d$row_DSE / den
    ek$costly_dse_measure <- ek$costly_dse_measure * ek$costly_dse_measure
  } else ek$costly_dse_measure <- 0
  costly <- ek$costly_dse_measure > 1000 && d$row_DSE > 0.01
  ek$costly_dse_frequency <- (1 - 0.05) * ek$costly_dse_frequency
  switch_to_devex <- FALSE
  if (costly) {
    ek$num_costly_dse_iteration <- ek$num_costly_dse_iteration + 1L
    ek$costly_dse_frequency <- ek$costly_dse_frequency + 0.05 * 1.0
    local_iteration_count <- ek$iteration_count - ek$control_iteration_count0
    switch_to_devex <- ek$allow_dse_to_devex &&
      ek$num_costly_dse_iteration > local_iteration_count * 0.05 &&
      local_iteration_count > 0.1 * (ek$n + ek$m)
  }
  if (!switch_to_devex) {
    switch_to_devex <- ek$allow_dse_to_devex &&
      ek$average_log_low_dse_weight_error + ek$average_log_high_dse_weight_error > 10
  }
  switch_to_devex
}

.hs_density_update <- function(ek, dens, count) {
  if (is.null(ek$dens)) .hs_init_densities(ek)
  ek$dens[[dens]] <- (1 - 0.05) * ek$dens[[dens]] + 0.05 * (count / ek$m)
}

# FTRAN and BTRAN of r (a dense vector or an HVector list) with the
# expected density dens; the index list of the result is left in
# ek$last_index and pack keeps the vectors of the FT update
.hs_ftran <- function(ek, r, dens = "col_aq", pack = FALSE) {
  if (is.null(ek$dens)) .hs_init_densities(ek)
  hv <- if (is.list(r)) r else .hs_hv(r)
  hv <- .hs_factor_ftran(ek$fa, hv, ek$dens[[dens]], pack)
  .hs_density_update(ek, dens, hv$count)
  ek$last_index <- hv$index
  ek$last_tick <- hv$tick %||% 0
  hv$array
}
.hs_btran <- function(ek, r, dens = "row_ep", pack = FALSE) {
  if (is.null(ek$dens)) .hs_init_densities(ek)
  hv <- if (is.list(r)) r else .hs_hv(r)
  hv <- .hs_factor_btran(ek$fa, hv, ek$dens[[dens]], pack)
  .hs_density_update(ek, dens, hv$count)
  ek$last_index <- hv$index
  ek$last_tick <- hv$tick %||% 0
  hv$array
}

# a sum of doubles in the given order, without the extended precision
# accumulator of sum()
.hs_dsum <- function(x) {
  s <- 0
  for (v in x) s <- s + v
  s
}

# the row of the tableau over the nonbasic columns
# HighsSparseMatrix::priceByColumn: each column sums its nonzeros times y
# in ascending row order; sums of magnitude at most kHighsTiny are zero
.hs_price_col <- function(ek, y) {
  A <- ek$A
  ap <- numeric(ek$n)
  for (i in seq_len(ek$m)) {
    nz <- A[i, ] != 0
    ap[nz] <- ap[nz] + y[i] * A[i, nz]
  }
  ap[abs(ap) <= 1e-14] <- 0
  ap
}

# HEkk::tableauRowPrice with the default kSimplexPriceStrategyRowSwitchColSwitch:
# HighsSparseMatrix::priceByColumn (rows ascending in each column) when more
# than 0.75 of row_ep is nonzero, else priceByRowWithSwitch, which adds the
# rows of the nonbasic columns in the order of the row_ep index list, keeps
# each partial sum below kHighsTiny at kHighsZero and, after a switch to
# priceByRowDenseResult, zeroes the final sums below kHighsTiny. index is the
# HVector index list of row_ep; dens the expected density of the result.
.hs_price <- function(ek, row_ep, index = which(row_ep != 0), dens = ek$dens$row_ap %||% 0) {
  n <- ek$n; A <- ek$A
  nonbasic <- ek$nonbasicFlag[seq_len(n)] != 0L
  if (length(index) / ek$m > 0.75) {
    ap <- .hs_price_col(ek, row_ep)
    ap[!nonbasic] <- 0
    return(ap)
  }
  ap <- numeric(n)
  cnt <- 0L
  sparse <- dens <= 0.1
  for (i in index) {
    cols <- which(A[i, ] != 0 & nonbasic)
    if (sparse && (cnt + length(cols) >= n || cnt / n > 0.1)) sparse <- FALSE
    mult <- row_ep[i]
    if (sparse) {
      if (!mult) next
      cnt <- cnt + sum(ap[cols] == 0)
    }
    v <- ap[cols] + mult * A[i, cols]
    v[abs(v) < 1e-14] <- 1e-50
    ap[cols] <- v
  }
  if (!sparse) ap[abs(ap) < 1e-14] <- 0
  ap
}

# HEkk::updateFactor: the FT update with the packed spike of the last
# col_aq FTRAN and the packed row eta of the last row_ep BTRAN
.hs_update_inverse <- function(ek, col_aq, r) {
  .hs_factor_update_ft(ek$fa, col_aq[r], r)
  ek$update_count <- ek$update_count + 1L
}

# the rebuild reason HEkk::updateFactor hints after an update: the update
# limit (6), or the synthetic clock (10) once the ticks of the solves since
# the rebuild reach those of the last build and 50 updates are done
# (kSyntheticTickReinversionMinUpdateCount)
.hs_update_hint <- function(ek, reason) {
  if (ek$update_count >= ek$update_limit) reason <- 6L
  if ((ek$total_tick %||% 0) >= (ek$fa$build_tick %||% Inf) && ek$update_count >= 50L) reason <- 10L
  reason
}

# HEkk::resetSyntheticClock
.hs_reset_clock <- function(ek) ek$total_tick <- 0

.hs_nearest_pow2_scale <- function(v) {
  if (v <= 0) return(1)
  e <- ceiling(log2(v))
  if (2^e < v) e <- e + 1
  if (2^(e - 1) >= v) e <- e - 1
  2^(-e)
}

# HEkk::initialiseSimplexLpRandomVectors
.hs_random_vectors <- function(ek) {
  n <- ek$n; tot <- n + ek$m
  if (n) ek$numColPermutation <- .hs_rng_shuffle(ek$rng, seq_len(n) - 1L)
  ek$numTotPermutation <- .hs_rng_shuffle(ek$rng, seq_len(tot) - 1L)
  ek$numTotRandomValue <- vapply(seq_len(tot), function(i) .hs_rng_fraction(ek$rng), 0)
}

.hs_basis_key <- function(ek) paste(sort(ek$basicIndex), collapse = ",")

# HEkk::setBasis for the slack basis with HEkk::setNonbasicMove
.hs_slack_basis <- function(ek) {
  n <- ek$n; m <- ek$m
  ek$basicIndex <- n + seq_len(m)
  ek$nonbasicFlag <- c(rep(1L, n), rep(0L, m))
  mv <- integer(n + m)
  for (j in seq_len(n)) {
    lo <- ek$lp$col_lower[j]; up <- ek$lp$col_upper[j]
    mv[j] <- if (lo == up) 0L else if (is.finite(lo)) {
      if (is.finite(up)) { if (abs(lo) < abs(up)) 1L else -1L } else 1L
    } else if (is.finite(up)) -1L else 0L
  }
  ek$nonbasicMove <- mv
  .hs_invert(ek)
}

.hs_lp_bounds <- function(ek) {
  list(lower = c(ek$lp$col_lower, -ek$lp$row_upper), upper = c(ek$lp$col_upper, -ek$lp$row_lower))
}

# HEkk::initialiseCost
.hs_init_cost <- function(ek, algorithm, perturb = FALSE) {
  n <- ek$n; tot <- n + ek$m
  ek$workCost <- c(ek$lp$col_cost, numeric(ek$m))
  ek$workShift <- numeric(tot)
  ek$costs_shifted <- FALSE
  ek$costs_perturbed <- FALSE
  if (algorithm == "primal") return(invisible())
  if (!perturb || ek$dual_cost_perturbation_multiplier == 0) return(invisible())
  max_abs_cost <- if (n) max(abs(ek$workCost[seq_len(n)])) else 0
  if (max_abs_cost > 100) max_abs_cost <- sqrt(sqrt(max_abs_cost))
  boxed_rate <- sum((ek$workUpper - ek$workLower) < 1e30) / tot
  if (boxed_rate < 0.01) max_abs_cost <- min(max_abs_cost, 1)
  base <- ek$dual_cost_perturbation_multiplier * 5e-7 * max_abs_cost
  rv <- ek$numTotRandomValue
  for (i in seq_len(n)) {
    lo <- ek$lp$col_lower[i]; up <- ek$lp$col_upper[i]
    xpert <- (1 + rv[i]) * (abs(ek$workCost[i]) + 1) * base
    if (lo == -Inf && up == Inf) {
    } else if (up == Inf) {
      ek$workCost[i] <- ek$workCost[i] + xpert
    } else if (lo == -Inf) {
      ek$workCost[i] <- ek$workCost[i] - xpert
    } else if (lo != up) {
      ek$workCost[i] <- ek$workCost[i] + (if (ek$workCost[i] >= 0) xpert else -xpert)
    }
  }
  rb <- ek$dual_cost_perturbation_multiplier * 1e-12
  for (i in (n + 1L):tot) ek$workCost[i] <- ek$workCost[i] + (0.5 - rv[i]) * rb
  ek$costs_perturbed <- TRUE
}

# HEkk::initialiseBound; phase is 1, 2 or 0 (unknown)
.hs_init_bound <- function(ek, algorithm, phase, perturb = FALSE) {
  b <- .hs_lp_bounds(ek)
  ek$workLower <- b$lower
  ek$workUpper <- b$upper
  tot <- ek$n + ek$m
  ek$workLowerShift <- numeric(tot)
  ek$workUpperShift <- numeric(tot)
  ek$bounds_perturbed <- FALSE
  if (algorithm == "primal") {
    if (!perturb || ek$primal_bound_perturbation_multiplier == 0) return(invisible())
    base <- ek$primal_bound_perturbation_multiplier * 5e-7
    for (v in seq_len(tot)) {
      lo <- ek$workLower[v]; up <- ek$workUpper[v]
      if (ek$nonbasicFlag[v] == 1L && lo == up) next
      r <- ek$numTotRandomValue[v]
      if (lo > -Inf) {
        lo <- if (lo < -1) lo - r * base * (-lo) else if (lo < 1) lo - r * base else lo - r * base * lo
        ek$workLower[v] <- lo
      }
      if (up < Inf) {
        up <- if (up < -1) up + r * base * (-up) else if (up < 1) up + r * base else up + r * base * up
        ek$workUpper[v] <- up
      }
      if (ek$nonbasicFlag[v] == 0L) next
      if (ek$nonbasicMove[v] > 0) ek$workValue[v] <- lo
      else if (ek$nonbasicMove[v] < 0) ek$workValue[v] <- up
    }
    ek$baseLower <- ek$workLower[ek$basicIndex]
    ek$baseUpper <- ek$workUpper[ek$basicIndex]
    ek$bounds_perturbed <- TRUE
    return(invisible())
  }
  if (phase == 2L) return(invisible())
  lo <- ek$workLower; up <- ek$workUpper
  free <- lo == -Inf & up == Inf
  upper <- !free & lo == -Inf
  lower <- !free & up == Inf
  boxed <- !free & !upper & !lower
  ek$workLower[free] <- -1000; ek$workUpper[free] <- 1000
  ek$workLower[upper] <- -1; ek$workUpper[upper] <- 0
  ek$workLower[lower] <- 0; ek$workUpper[lower] <- 1
  ek$workLower[boxed] <- 0; ek$workUpper[boxed] <- 0
}

# HEkk::initialiseNonbasicValueAndMove
.hs_init_nonbasic <- function(ek) {
  tot <- ek$n + ek$m
  if (is.null(ek$workValue)) ek$workValue <- numeric(tot)
  for (v in seq_len(tot)) {
    if (ek$nonbasicFlag[v] == 0L) { ek$nonbasicMove[v] <- 0L; next }
    lo <- ek$workLower[v]; up <- ek$workUpper[v]
    om <- ek$nonbasicMove[v]
    if (lo == up) { val <- lo; mv <- 0L }
    else if (lo > -Inf) {
      if (up < Inf) {
        if (om == 1L) { val <- lo; mv <- 1L }
        else if (om == -1L) { val <- up; mv <- -1L }
        else { val <- lo; mv <- 1L }
      } else { val <- lo; mv <- 1L }
    } else if (up < Inf) { val <- up; mv <- -1L }
    else { val <- 0; mv <- 0L }
    ek$nonbasicMove[v] <- mv
    ek$workValue[v] <- val
  }
}

# HEkk::computePrimal
.hs_compute_primal <- function(ek) {
  tot <- ek$n + ek$m
  hv <- list(array = numeric(ek$m), index = integer(), count = 0L)
  for (v in seq_len(tot)) if (ek$nonbasicFlag[v] == 1L && ek$workValue[v] != 0) hv <- .hs_hv_collect(ek, hv, v, ek$workValue[v])
  pc <- if (hv$count) .hs_ftran(ek, hv, "primal_col") else hv$array
  ek$baseValue <- -pc
  ek$baseLower <- ek$workLower[ek$basicIndex]
  ek$baseUpper <- ek$workUpper[ek$basicIndex]
}

# HEkk::computeDual
.hs_compute_dual <- function(ek) {
  tot <- ek$n + ek$m
  cb <- ek$workCost[ek$basicIndex] + ek$workShift[ek$basicIndex]
  ek$workDual <- ek$workCost + ek$workShift
  if (any(cb != 0)) {
    y <- .hs_btran(ek, cb, "dual_col")
    row <- .hs_price_col(ek, y)
    ek$workDual <- ek$workDual - c(row, y)
  }
}

.hs_primal_infeas <- function(value, lower, upper, tol) {
  ifelse(value < lower - tol, lower - value, ifelse(value > upper + tol, value - upper, 0))
}

# HEkk::computeSimplexPrimalInfeasible
.hs_simplex_primal_infeasible <- function(ek) {
  tol <- ek$ptol
  nb <- ek$nonbasicFlag == 1L
  a <- .hs_primal_infeas(ek$workValue[nb], ek$workLower[nb], ek$workUpper[nb], tol)
  b <- .hs_primal_infeas(ek$baseValue, ek$baseLower, ek$baseUpper, tol)
  inf <- c(a, b)
  ek$num_primal_infeasibilities <- sum(inf > tol)
  ek$max_primal_infeasibility <- if (length(inf)) max(0, inf) else 0
  ek$sum_primal_infeasibilities <- .hs_dsum(inf)
}

# HEkk::computeSimplexDualInfeasible (also
# computeDualInfeasibilitiesWithFixedVariableFlips, which is the same)
.hs_simplex_dual_infeasible <- function(ek) {
  nb <- which(ek$nonbasicFlag == 1L)
  d <- ek$workDual[nb]
  free <- ek$workLower[nb] == -Inf & ek$workUpper[nb] == Inf
  inf <- ifelse(free, abs(d), -ek$nonbasicMove[nb] * d)
  inf <- inf[inf > 0]
  ek$num_dual_infeasibilities <- sum(inf >= ek$dtol)
  ek$max_dual_infeasibility <- if (length(inf)) max(inf) else 0
  ek$sum_dual_infeasibilities <- .hs_dsum(inf)
}

# HEkk::computeSimplexLpDualInfeasible: dual infeasibilities with respect
# to the bounds of the LP
.hs_lp_dual_infeasible_count <- function(ek) {
  n <- ek$n; cnt <- 0L
  dinf <- function(dual, lo, up) {
    if (up == Inf) { if (lo == -Inf) abs(dual) else -dual }
    else if (lo == -Inf) dual else 0
  }
  for (j in seq_len(n)) if (ek$nonbasicFlag[j] == 1L) {
    if (dinf(ek$workDual[j], ek$lp$col_lower[j], ek$lp$col_upper[j]) >= ek$dtol) cnt <- cnt + 1L
  }
  for (i in seq_len(ek$m)) if (ek$nonbasicFlag[n + i] == 1L) {
    if (dinf(-ek$workDual[n + i], ek$lp$row_lower[i], ek$lp$row_upper[i]) >= ek$dtol) cnt <- cnt + 1L
  }
  cnt
}

.hs_compute_dual_objective <- function(ek) {
  nb <- ek$nonbasicFlag == 1L
  ek$dual_objective_value <- .hs_dsum(ek$workValue[nb] * ek$workDual[nb])
}

.hs_compute_primal_objective <- function(ek) {
  n <- ek$n
  v <- 0
  for (i in seq_len(ek$m)) {
    j <- ek$basicIndex[i]
    if (j <= n) v <- v + ek$baseValue[i] * ek$lp$col_cost[j]
  }
  for (j in seq_len(n)) if (ek$nonbasicFlag[j] == 1L) v <- v + ek$workValue[j] * ek$lp$col_cost[j]
  ek$primal_objective_value <- v
}

.hs_flip_bound <- function(ek, v) {
  mv <- -ek$nonbasicMove[v]
  ek$nonbasicMove[v] <- mv
  ek$workValue[v] <- if (mv == 1L) ek$workLower[v] else ek$workUpper[v]
}

# HEkk::updatePivots, with the update of the inverse
.hs_update_pivots <- function(ek, vin, row_out, move_out, col_aq) {
  vout <- ek$basicIndex[row_out]
  ek$basicIndex[row_out] <- vin
  ek$visited[[.hs_basis_key(ek)]] <- TRUE
  ek$nonbasicFlag[vin] <- 0L
  ek$nonbasicMove[vin] <- 0L
  ek$baseLower[row_out] <- ek$workLower[vin]
  ek$baseUpper[row_out] <- ek$workUpper[vin]
  ek$nonbasicFlag[vout] <- 1L
  if (ek$workLower[vout] == ek$workUpper[vout]) {
    ek$workValue[vout] <- ek$workLower[vout]; ek$nonbasicMove[vout] <- 0L
  } else if (move_out == -1L) {
    ek$workValue[vout] <- ek$workLower[vout]; ek$nonbasicMove[vout] <- 1L
  } else {
    ek$workValue[vout] <- ek$workUpper[vout]; ek$nonbasicMove[vout] <- -1L
  }
  .hs_update_inverse(ek, col_aq, row_out)
  ek$has_fresh_rebuild <- FALSE
}

# HEkk::isBadBasisChange, with the visited bases kept by their sorted basic
# sets instead of a hash
.hs_is_bad_basis_change <- function(ek, vin, row_out, reason) {
  if (reason != 0L) return(FALSE)
  if (vin < 1L || row_out < 1L) return(FALSE)
  vout <- ek$basicIndex[row_out]
  bi <- ek$basicIndex
  bi[row_out] <- vin
  key <- paste(sort(bi), collapse = ",")
  if (!is.null(ek$visited[[key]])) {
    if (ek$iteration_count == ek$previous_cycling + 1L) {
      .hs_add_bad_basis_change(ek, row_out, vout, vin, "cycling", TRUE)
      return(TRUE)
    }
    ek$previous_cycling <- ek$iteration_count
    return(FALSE)
  }
  for (k in seq_along(ek$bad)) {
    b <- ek$bad[[k]]
    if (b$vout == vout && b$vin == vin && b$row == row_out) {
      ek$bad[[k]]$taboo <- TRUE
      return(TRUE)
    }
  }
  FALSE
}

.hs_add_bad_basis_change <- function(ek, row_out, vout, vin, reason, taboo) {
  for (k in seq_along(ek$bad)) {
    b <- ek$bad[[k]]
    if (b$row == row_out && b$vout == vout && b$vin == vin && b$reason == reason) {
      ek$bad[[k]]$taboo <- taboo
      return(invisible())
    }
  }
  ek$bad[[length(ek$bad) + 1L]] <- list(row = row_out, vout = vout, vin = vin, reason = reason, taboo = taboo)
}

.hs_taboo <- function(ek) any(vapply(ek$bad, function(b) b$taboo, TRUE))
.hs_clear_taboo <- function(ek) for (k in seq_along(ek$bad)) ek$bad[[k]]$taboo <- FALSE

.hs_update_bad_basis_change <- function(ek, col_aq, theta_primal) {
  if (!length(ek$bad)) return(invisible())
  keep <- vapply(ek$bad, function(b) abs(col_aq[b$row] * theta_primal) < ek$ptol, TRUE)
  ek$bad <- ek$bad[keep]
}

.hs_bailout <- function(ek) {
  if (ek$solve_bailout) return(TRUE)
  if (ek$iteration_count >= ek$iteration_limit) {
    ek$solve_bailout <- TRUE
    ek$model_status <- "IterationLimit"
  }
  ek$solve_bailout
}

# HEkk::rebuildRefactor: with no_unnecessary_rebuild_refactor the solution
# error of the dense inverse never asks for a refactorization
.hs_rebuild_refactor <- function(ek, reason) {
  if (ek$update_count == 0L) return(FALSE)
  !(reason %in% c(0L, 1L, 2L, 3L, 4L, 5L))
}

# rebuild reasons: 0 no, 1 possibly optimal, 2 possibly phase 1 feasible,
# 3 possibly primal unbounded, 4 possibly dual unbounded, 5 primal
# infeasible in primal simplex, 6 update limit, 7 possibly singular basis,
# 8 choose column fail, 9 excessive primal value

# HEkk::putBacktrackingBasis: the basis (with basic_index in its order
# before the refactor), the cost shift and perturbation flags, the cost
# shifts and the dual edge weights scattered over the variables
.hs_put_backtracking <- function(ek, basic_index = ek$basicIndex, scattered = NULL) {
  if (is.null(scattered)) {
    scattered <- numeric(ek$n + ek$m)
    if (!is.null(ek$dual_edge_weight)) scattered[ek$basicIndex] <- ek$dual_edge_weight
  }
  ek$bt <- list(basicIndex = basic_index, nonbasicFlag = ek$nonbasicFlag, nonbasicMove = ek$nonbasicMove,
                costs_shifted = ek$costs_shifted, costs_perturbed = ek$costs_perturbed,
                workShift = ek$workShift, weight = scattered)
}

# HEkk::getNonsingularInverse: refactor and keep the basis as the
# backtracking basis, or after a rank deficient refactor go back to the
# last nonsingular basis with half the update limit; FALSE when that is not
# possible. The dual edge weights follow their basic variables through the
# permutation of the refactor.
.hs_reinvert <- function(ek) {
  before <- ek$basicIndex
  updates <- ek$update_count
  scattered <- numeric(ek$n + ek$m)
  if (!is.null(ek$dual_edge_weight)) scattered[ek$basicIndex] <- ek$dual_edge_weight
  if (.hs_invert(ek, repair = FALSE) > 0L) {
    if (is.null(ek$bt)) return(FALSE)
    deficient <- .hs_basis_key(ek)
    bt <- ek$bt
    ek$basicIndex <- bt$basicIndex; ek$nonbasicFlag <- bt$nonbasicFlag; ek$nonbasicMove <- bt$nonbasicMove
    ek$costs_shifted <- bt$costs_shifted; ek$costs_perturbed <- bt$costs_perturbed
    ek$workShift <- bt$workShift
    scattered <- bt$weight
    ek$backtracking <- TRUE
    ek$visited <- new.env(parent = emptyenv())
    ek$visited[[.hs_basis_key(ek)]] <- TRUE
    ek$visited[[deficient]] <- TRUE
    ek$has_fresh_rebuild <- FALSE
    if (.hs_invert(ek, repair = FALSE) > 0L) return(FALSE)
    if (updates <= 1L) return(FALSE)
    ek$update_limit <- updates %/% 2L
  } else {
    .hs_put_backtracking(ek, before, scattered)
    ek$backtracking <- FALSE
    ek$update_limit <- 5000L
  }
  if (!is.null(ek$dual_edge_weight)) ek$dual_edge_weight <- scattered[ek$basicIndex]
  TRUE
}

# HEkkDual::computeDualSteepestEdgeWeights for a nonlogical basis
.hs_compute_dse_weights <- function(ek) {
  ek$dual_edge_weight <- vapply(seq_len(ek$m), function(i) {
    e <- numeric(ek$m); e[i] <- 1
    y <- .hs_btran(ek, e)
    .hs_dsum(y[ek$last_index]^2)
  }, 0)
}
