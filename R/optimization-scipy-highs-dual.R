# The dual simplex method of HiGHS 1.12.0 (simplex/HEkkDual.cpp,
# HEkkDualRow.cpp, HEkkDualRHS.cpp), serial strategy, as SciPy 1.18.1
# vendors it. MIT license notice in R/optimization-scipy-highs-util.R.
# Phases: 1, 2, 0 unknown, -1 exit, -2 optimal, -3 optimal cleanup,
# -4 primal infeasible cleanup, -5 error, -6 taboo basis.

# HEkkDualRHS::createArrayOfPrimalInfeasibilities and createInfeasList
.hs_rhs_create <- function(ek) {
  inf <- .hs_primal_infeas(ek$baseValue, ek$baseLower, ek$baseUpper, ek$ptol)
  ek$work_infeasibility <- if (ek$store_squared) inf * inf else abs(inf)
  m <- ek$m
  ek$workMark <- ek$work_infeasibility != 0
  ek$workIndex <- which(ek$workMark)
  # ponytail: the cutoff of HiGHS for more than max(0.01 m, 500) infeasible
  # rows is left out; it only acts on LPs with over 2500 rows
  ek$workCountDense <- length(ek$workIndex) > 0.2 * m
}

# HEkkDualRHS::updatePrimal
.hs_rhs_update_primal <- function(ek, column, theta) {
  rows <- which(column != 0)
  if (!length(rows)) return(TRUE)
  ek$baseValue[rows] <- ek$baseValue[rows] - theta * column[rows]
  inf <- .hs_primal_infeas(ek$baseValue[rows], ek$baseLower[rows], ek$baseUpper[rows], ek$ptol)
  ek$work_infeasibility[rows] <- if (ek$store_squared) inf * inf else abs(inf)
  !any(abs(ek$baseValue[rows]) >= 1e25)
}

# HEkkDualRHS::updateInfeasList
.hs_rhs_update_list <- function(ek, column) {
  if (ek$workCountDense) return(invisible())
  for (i in which(column != 0)) {
    if (!ek$workMark[i] && ek$work_infeasibility[i] != 0) {
      ek$workIndex <- c(ek$workIndex, i)
      ek$workMark[i] <- TRUE
    }
  }
}

.hs_rhs_update_pivots <- function(ek, r, value) {
  ek$baseValue[r] <- value
  inf <- .hs_primal_infeas(value, ek$baseLower[r], ek$baseUpper[r], ek$ptol)
  ek$work_infeasibility[r] <- if (ek$store_squared) inf * inf else abs(inf)
}

# HEkkDualRHS::chooseNormal
.hs_rhs_choose <- function(ek) {
  rows <- if (ek$workCountDense) seq_len(ek$m) else ek$workIndex
  cnt <- length(rows)
  if (cnt == 0L) return(-1L)
  start <- .hs_rng_integer(ek$rng, cnt)
  ord <- c(if (start < cnt) (start + 1L):cnt, if (start > 0L) seq_len(start))
  best_merit <- 0; best <- -1L
  w <- ek$dual_edge_weight
  for (k in ord) {
    i <- rows[k]
    inf <- ek$work_infeasibility[i]
    if (inf > 0x1.dee7a4ad4b81fp-167 && best_merit * w[i] < inf) {
      best_merit <- inf / w[i]
      best <- i
    }
  }
  best
}

.hs_dual_init_devex <- function(ek) {
  ek$devex_index <- 1L - ek$nonbasicFlag * ek$nonbasicFlag
  ek$dual_edge_weight <- rep(1, ek$m)
  ek$num_devex_iterations <- 0L
  ek$new_devex_framework <- FALSE
}

# HEkkDual::correctDualInfeasibilities
.hs_correct_dual <- function(ek, ds) {
  free_cnt <- 0L
  tot <- ek$n + ek$m
  for (v in seq_len(tot)) {
    if (ek$nonbasicFlag[v] == 0L) next
    lo <- ek$workLower[v]; up <- ek$workUpper[v]
    d <- ek$workDual[v]; mv <- ek$nonbasicMove[v]
    fixed <- lo == up
    boxed <- lo > -Inf && up < Inf
    if (lo == -Inf && up == Inf) {
      if (abs(d) >= ek$dtol) free_cnt <- free_cnt + 1L
      next
    }
    if (-mv * d < ek$dtol) next
    if (fixed || (boxed && !ds$force_phase2)) {
      .hs_flip_bound(ek, v)
      next
    }
    ek$costs_shifted <- TRUE
    nd <- (1 + .hs_rng_fraction(ek$rng)) * ek$dtol
    if (mv != 1L) nd <- -nd
    shift <- nd - d
    ek$workDual[v] <- nd
    ek$workCost[v] <- ek$workCost[v] + shift
  }
  ds$force_phase2 <- FALSE
  ds$dualInfeasCount <- free_cnt
}

# HEkkDual::rebuild
.hs_dual_rebuild <- function(ek, ds) {
  .hs_clear_taboo(ek)
  if (.hs_rebuild_refactor(ek, ds$rebuild_reason) && !.hs_reinvert(ek)) {
    ds$solve_phase <- -5L
    return(invisible())
  }
  ds$rebuild_reason <- 0L
  .hs_compute_dual(ek)
  if (isTRUE(ek$backtracking)) {
    ds$solve_phase <- 0L
    return(invisible())
  }
  .hs_correct_dual(ek, ds)
  .hs_compute_primal(ek)
  .hs_rhs_create(ek)
  .hs_compute_dual_objective(ek)
  ek$has_fresh_rebuild <- TRUE
  .hs_reset_clock(ek)
}

# HEkkDual::cleanup
.hs_dual_cleanup <- function(ek, ds) {
  if (ds$solve_phase == 1L) ek$dual_phase1_cleanup_level <- ek$dual_phase1_cleanup_level + 1L
  .hs_init_cost(ek, "dual", FALSE)
  ek$allow_cost_perturbation <- FALSE
  .hs_init_bound(ek, "dual", ds$solve_phase)
  .hs_compute_dual(ek)
  .hs_simplex_dual_infeasible(ek)
  ds$dualInfeasCount <- ek$num_dual_infeasibilities
  .hs_compute_dual_objective(ek)
}

# HEkkDual::chooseRow
.hs_dual_choose_row <- function(ek, ds) {
  if (ds$rebuild_reason != 0L) return(invisible())
  saved <- list()
  for (b in ek$bad) if (b$taboo) {
    saved[[length(saved) + 1L]] <- c(b$row, ek$work_infeasibility[b$row])
    ek$work_infeasibility[b$row] <- 0
  }
  repeat {
    r <- .hs_rhs_choose(ek)
    ds$row_out <- r
    if (r < 0L) {
      ds$rebuild_reason <- 1L
      return(invisible())
    }
    e <- numeric(ek$m); e[r] <- 1
    ds$row_ep <- .hs_btran(ek, e, "row_ep", pack = TRUE)
    ds$row_ep_index <- ek$last_index
    ds$row_ep_tick <- ek$last_tick
    if (ds$edge_weight_mode == "dse") {
      updated <- ek$dual_edge_weight[r]
      computed <- .hs_dsum(ds$row_ep[ds$row_ep_index]^2)
      ek$dual_edge_weight[r] <- computed
      .hs_assess_dse_weight_error(ek, computed, updated)
      if (updated >= 0.25 * computed) break
    } else break
  }
  for (s in rev(saved)) ek$work_infeasibility[s[1L]] <- s[2L]
  r <- ds$row_out
  ds$variable_out <- ek$basicIndex[r]
  ds$delta_primal <- if (ek$baseValue[r] < ek$baseLower[r]) ek$baseValue[r] - ek$baseLower[r] else ek$baseValue[r] - ek$baseUpper[r]
  ds$move_out <- if (ds$delta_primal < 0) -1L else 1L
}

# HEkkDualRow::chooseFinal with chooseFinalWorkGroupQuad and
# chooseFinalLargeAlpha. Returns the pivot, alpha, theta and the flips, or
# NULL when the group search fails.
.hs_choose_final <- function(ek, wi, wa, work_theta, work_delta) {
  move <- ek$nonbasicMove; dual <- ek$workDual
  range <- ek$workUpper - ek$workLower
  full <- length(wi)
  cnt <- 0L
  total_change <- 0
  total_delta <- abs(work_delta)
  select_theta <- 10 * work_theta + 1e-7
  repeat {
    i <- cnt + 1L
    while (i <= full) {
      j <- wi[i]
      if (wa[i] * select_theta >= move[j] * dual[j]) {
        cnt <- cnt + 1L
        t <- wi[cnt]; wi[cnt] <- wi[i]; wi[i] <- t
        t <- wa[cnt]; wa[cnt] <- wa[i]; wa[i] <- t
        total_change <- total_change + range[wi[cnt]] * wa[cnt]
      }
      i <- i + 1L
    }
    select_theta <- select_theta * 10
    if (total_change >= total_delta || cnt == full) break
  }
  # chooseFinalWorkGroupQuad over the cnt preselected candidates
  full <- cnt
  cnt <- 0L
  total_change <- 1e-12
  select_theta <- work_theta
  groups <- 0L
  prev_cnt <- cnt; prev_remain <- 0x1.249ad2594c37dp+332; prev_select <- select_theta
  while (select_theta < 1e18) {
    remain <- 0x1.249ad2594c37dp+332
    i <- cnt + 1L
    while (i <= full) {
      j <- wi[i]
      d <- move[j] * dual[j]
      if (d <= select_theta * wa[i]) {
        cnt <- cnt + 1L
        t <- wi[cnt]; wi[cnt] <- wi[i]; wi[i] <- t
        t <- wa[cnt]; wa[cnt] <- wa[i]; wa[i] <- t
        total_change <- total_change + wa[cnt] * range[wi[cnt]]
      } else if (d + ek$dtol < remain * wa[i]) {
        remain <- (d + ek$dtol) / wa[i]
      }
      i <- i + 1L
    }
    groups <- c(groups, cnt)
    select_theta <- remain
    if (cnt == prev_cnt && prev_select == select_theta && prev_remain == remain) return(NULL)
    prev_cnt <- cnt; prev_remain <- remain; prev_select <- select_theta
    if (total_change >= total_delta || cnt == full) break
  }
  if (length(groups) <= 1L) return(NULL)
  # chooseFinalLargeAlpha
  final_compare <- if (cnt) min(0.1 * max(0, wa[seq_len(cnt)]), 1) else 0
  ng <- length(groups) - 1L
  bi <- -1L; bg <- -1L
  perm <- ek$numTotPermutation
  for (g in ng:1L) {
    dmax <- 0; imax <- -1L
    if (groups[g + 1L] > groups[g]) for (i in (groups[g] + 1L):groups[g + 1L]) {
      if (dmax < wa[i]) { dmax <- wa[i]; imax <- i }
      else if (dmax == wa[i] && imax > 0L && perm[wi[i]] < perm[wi[imax]]) imax <- i
    }
    if (imax > 0L && wa[imax] > final_compare) { bi <- imax; bg <- g; break }
  }
  if (bi < 0L) return(NULL)
  move_out <- if (work_delta < 0) -1L else 1L
  piv <- wi[bi]
  alpha <- wa[bi] * move_out * move[piv]
  theta <- if (dual[piv] * move[piv] > 0) dual[piv] / alpha else 0
  nf <- groups[bg]
  flips <- if (nf > 0L && theta != 0) wi[seq_len(nf)] else integer()
  flips <- sort(flips)
  list(pivot = piv, alpha = alpha, theta = theta, flips = flips, flip_change = move[flips] * range[flips])
}

# HEkkDual::chooseColumn
.hs_dual_choose_column <- function(ek, ds) {
  if (ds$rebuild_reason != 0L) return(invisible())
  n <- ek$n
  row_ap <- .hs_price(ek, ds$row_ep, ds$row_ep_index)
  ds$row_ap <- row_ap
  # ponytail: the row_ap count is that of the nonzeros over the nonbasic
  # columns, which priceByRowWithSwitch reports unless it switches to dense
  .hs_density_update(ek, "row_ap", sum(row_ap != 0) * ek$m / max(1L, ek$n))
  # createFreemove
  Ta <- if (ek$update_count < 10L) 1e-9 else if (ek$update_count < 20L) 3e-8 else 1e-6
  for (v in ds$free_list) {
    alpha <- if (v <= n) row_ap[v] else ds$row_ep[v - n]
    if (abs(alpha) > Ta) ek$nonbasicMove[v] <- if (alpha * ds$move_out > 0) 1L else -1L
  }
  pidx <- c(which(row_ap != 0), n + which(ds$row_ep != 0))
  pval <- c(row_ap[row_ap != 0], ds$row_ep[ds$row_ep != 0])
  row_ep_scale <- .hs_nearest_pow2_scale(if (length(pval)) max(abs(pval)) else 0)
  pass <- 0L
  repeat {
    # choosePossible
    alpha <- pval * ds$move_out * ek$nonbasicMove[pidx]
    ok <- alpha > Ta
    wi <- pidx[ok]; wa <- alpha[ok]
    relax <- ek$workDual[wi] * ek$nonbasicMove[wi] + ek$dtol
    work_theta <- Inf
    for (k in seq_along(wi)) if (work_theta * wa[k] > relax[k]) work_theta <- relax[k] / wa[k]
    ds$variable_in <- -1L
    if (work_theta <= 0 || !length(wi)) {
      ds$rebuild_reason <- 4L
      return(invisible())
    }
    cf <- .hs_choose_final(ek, wi, wa, work_theta, ds$delta_primal)
    if (is.null(cf)) {
      ds$rebuild_reason <- 8L
      return(invisible())
    }
    piv <- cf$pivot
    if (abs(row_ep_scale * cf$alpha) <= 1e-9) {
      if (pass > 0L) {
        k <- which(pidx == piv)[1L]
        last <- length(pidx)
        pidx[k] <- pidx[last]; pval[k] <- pval[last]
        pidx <- pidx[-last]; pval <- pval[-last]
      }
      piv <- -1L
    }
    if (piv >= 1L || !length(pidx)) break
    pass <- pass + 1L
  }
  for (v in ds$free_list) ek$nonbasicMove[v] <- 0L
  if (piv < 1L) {
    ds$rebuild_reason <- 4L
    return(invisible())
  }
  ds$variable_in <- piv
  ds$alpha_row <- cf$alpha
  ds$theta_dual <- cf$theta
  ds$flips <- cf$flips
  ds$flip_change <- cf$flip_change
  ds$pack_index <- pidx
  ds$pack_value <- pval
  if (ds$edge_weight_mode == "devex" && !ek$new_devex_framework) {
    nb <- ek$nonbasicFlag[pidx] == 1L
    pv <- ek$devex_index[pidx[nb]] * pval[nb]
    ds$computed_edge_weight <- max(1, .hs_dsum(pv * pv))
  }
}

# HEkk::proofOfPrimalInfeasibility
.hs_proof_of_infeasibility <- function(ek, row_ep, move_out) {
  n <- ek$n
  y <- row_ep
  proof_lower <- 0
  for (i in which(y != 0)) {
    if (abs(y[i] * max(abs(ek$A[i, ]))) <= 1e-9) { y[i] <- 0; next }
    y[i] <- y[i] * move_out
    if (y[i] > 0) {
      bound <- ek$lp$row_lower[i]
      if (bound == -Inf) { y[i] <- 0; next }
    } else {
      bound <- ek$lp$row_upper[i]
      if (bound == Inf) { y[i] <- 0; next }
    }
    proof_lower <- proof_lower + y[i] * bound
  }
  pv <- as.numeric(crossprod(ek$A, y))
  implied_upper <- 0
  sum_inf <- 0
  for (j in which(abs(pv) > 1e-14)) {
    v <- pv[j]
    if (v > 0) {
      if (ek$lp$col_upper[j] == Inf) {
        sum_inf <- sum_inf + v
        if (sum_inf > 1e-9) break
        next
      }
      implied_upper <- implied_upper + v * ek$lp$col_upper[j]
    } else {
      if (ek$lp$col_lower[j] == -Inf) {
        sum_inf <- sum_inf - v
        if (sum_inf > 1e-9) break
        next
      }
      implied_upper <- implied_upper + v * ek$lp$col_lower[j]
    }
  }
  !(sum_inf > 1e-9) && (proof_lower - implied_upper) > ek$ptol
}

# HEkkDual::iterate
.hs_dual_iterate <- function(ek, ds) {
  .hs_dual_choose_row(ek, ds)
  .hs_dual_choose_column(ek, ds)
  if (.hs_is_bad_basis_change(ek, ds$variable_in, ds$row_out, ds$rebuild_reason)) return(invisible())
  if (ds$rebuild_reason != 0L) return(invisible())
  n <- ek$n
  vin <- ds$variable_in; r <- ds$row_out
  # updateFtranBFRT
  # updateFtran
  hv <- .hs_hv_collect(ek, list(array = numeric(ek$m), index = integer(), count = 0L), vin, 1)
  col_aq <- .hs_ftran(ek, hv, "col_aq", pack = TRUE)
  col_aq_tick <- ek$last_tick
  alpha_col <- col_aq[r]
  hv <- list(array = numeric(ek$m), index = integer(), count = 0L)
  for (k in seq_along(ds$flips)) {
    v <- ds$flips[k]
    .hs_flip_bound(ek, v)
    hv <- .hs_hv_collect(ek, hv, v, ds$flip_change[k])
  }
  col_bfrt <- if (hv$count) .hs_ftran(ek, hv, "col_BFRT") else hv$array
  # the DSE FTRAN works in place on row_ep, whose clock it carries on
  dse_tick <- ds$row_ep_tick
  dse <- if (ds$edge_weight_mode == "dse") {
    out <- .hs_ftran(ek, list(array = ds$row_ep, index = ds$row_ep_index, count = length(ds$row_ep_index),
                              tick = ds$row_ep_tick), "row_DSE")
    dse_tick <- ek$last_tick
    out
  }
  # updateVerify
  trouble <- abs(abs(alpha_col) - abs(ds$alpha_row)) / min(abs(alpha_col), abs(ds$alpha_row))
  if (trouble > 1e-7 && ek$update_count > 0L) {
    ds$rebuild_reason <- 7L
    return(invisible())
  }
  # updateDual
  theta_dual <- ds$theta_dual
  if (theta_dual == 0) {
    ek$costs_shifted <- TRUE
    if (ek$workDual[vin] != 0) ek$workShift[vin] <- -ek$workDual[vin]
  } else {
    ek$workDual[ds$pack_index] <- ek$workDual[ds$pack_index] - theta_dual * ds$pack_value
  }
  vout <- ds$variable_out
  ek$workDual[vin] <- 0
  ek$workDual[vout] <- -theta_dual
  if (ek$workShift[vout] != 0) {
    ek$workDual[vout] <- ek$workDual[vout] - ek$workShift[vout]
    ek$workShift[vout] <- 0
  }
  # updatePrimal
  if (ds$edge_weight_mode == "devex") {
    updated <- ek$dual_edge_weight[r]
    ek$dual_edge_weight[r] <- ds$computed_edge_weight
    ratio <- max(updated / ds$computed_edge_weight, ds$computed_edge_weight / updated)
    ek$new_devex_framework <- !(ratio <= 9) || !(ek$num_devex_iterations <= max(25L, as.integer(ek$m / 1e-2)))
  }
  .hs_rhs_update_primal(ek, col_bfrt, 1)
  .hs_rhs_update_list(ek, col_bfrt)
  bound <- if (ds$delta_primal < 0) ek$baseLower[r] else ek$baseUpper[r]
  theta_primal <- (ek$baseValue[r] - bound) / alpha_col
  if (!.hs_rhs_update_primal(ek, col_aq, theta_primal)) {
    ds$rebuild_reason <- 9L
    return(invisible())
  }
  .hs_update_bad_basis_change(ek, col_aq, theta_primal)
  w <- ek$dual_edge_weight
  if (ds$edge_weight_mode == "dse") {
    npw <- w[r] / (alpha_col * alpha_col)
    kai <- -2 / alpha_col
    for (i in which(col_aq != 0)) w[i] <- max(1e-4, w[i] + col_aq[i] * (npw * col_aq[i] + kai * dse[i]))
    w[r] <- npw
  } else if (ds$edge_weight_mode == "devex") {
    npw <- max(1, w[r] / (alpha_col * alpha_col))
    for (i in which(col_aq != 0)) w[i] <- max(w[i], npw * col_aq[i] * col_aq[i])
    w[r] <- npw
    ek$num_devex_iterations <- ek$num_devex_iterations + 1L
  }
  ek$dual_edge_weight <- w
  .hs_rhs_update_list(ek, col_aq)
  ek$total_tick <- (ek$total_tick %||% 0) + col_aq_tick + dse_tick
  # updatePivots
  value_in <- ek$workValue[vin] + theta_primal
  .hs_update_pivots(ek, vin, r, ds$move_out, col_aq)
  ek$iteration_count <- ek$iteration_count + 1L
  ds$rebuild_reason <- .hs_update_hint(ek, ds$rebuild_reason)
  ds$free_list <- setdiff(ds$free_list, vin)
  .hs_rhs_update_pivots(ek, r, value_in)
  if (isTRUE(ek$new_devex_framework)) .hs_dual_init_devex(ek)
  # iterationAnalysis
  if (ds$edge_weight_mode == "dse" && .hs_switch_to_devex(ek)) {
    ds$edge_weight_mode <- "devex"
    .hs_dual_init_devex(ek)
  }
}

.hs_dual_assess_unbounded <- function(ek, ds) {
  if (ds$solve_phase != 2L || !ek$has_fresh_rebuild) return(invisible())
  if (.hs_proof_of_infeasibility(ek, ds$row_ep, ds$move_out)) {
    ds$solve_phase <- -1L
    ek$dual_ray <- list(row = ds$row_out, sign = ds$move_out)
    ek$model_status <- "Infeasible"
  } else {
    .hs_add_bad_basis_change(ek, ds$row_out, ds$variable_out, ds$variable_in, "proof", TRUE)
    ds$rebuild_reason <- 0L
  }
}

# the rebuild and iterate loop shared by the two phases
.hs_dual_loop <- function(ek, ds, phase2) {
  repeat {
    .hs_dual_rebuild(ek, ds)
    if (ds$solve_phase %in% c(0L, -5L)) return(invisible())
    if (.hs_bailout(ek)) break
    # HEkkDual::bailoutOnDualObjective (the MIP LP relaxation sets a bound)
    if (phase2 && is.finite(ek$objective_bound %||% Inf) &&
        exists(".hs_bailout_on_dual_objective") && .hs_bailout_on_dual_objective(ek, ds)) break
    if (phase2 && ds$dualInfeasCount > 0L) break
    repeat {
      .hs_dual_iterate(ek, ds)
      if (.hs_bailout(ek)) break
      if (phase2 && ds$rebuild_reason == 4L) .hs_dual_assess_unbounded(ek, ds)
      if (ds$rebuild_reason != 0L) break
    }
    if (ek$solve_bailout) break
    finished <- ek$has_fresh_rebuild && !.hs_rebuild_refactor(ek, ds$rebuild_reason)
    if (finished && .hs_taboo(ek)) {
      ds$solve_phase <- -6L
      return(invisible())
    }
    if (finished) break
  }
}

.hs_dual_assess_phase1 <- function(ek, ds) {
  if (ek$costs_perturbed) .hs_dual_cleanup(ek, ds)
  if (ds$dualInfeasCount == 0L) {
    if (ek$dual_objective_value == 0) {
      ds$solve_phase <- 2L
    } else if (.hs_lp_dual_infeasible_count(ek) == 0L) {
      ds$solve_phase <- 2L
    } else {
      ek$model_status <- "UnboundedOrInfeasible"
      ds$solve_phase <- -1L
    }
  }
  if (ds$dualInfeasCount == 0L && ds$solve_phase == 2L) {
    # exitPhase1ResetDuals
    if (!ek$costs_perturbed) {
      .hs_init_cost(ek, "dual", TRUE)
      .hs_compute_dual(ek)
    }
    lpb <- c(ek$lp$col_lower, ek$lp$row_lower); lpu <- c(ek$lp$col_upper, ek$lp$row_upper)
    for (v in which(ek$nonbasicFlag == 1L)) if (lpb[v] == -Inf && lpu[v] == Inf) {
      shift <- -ek$workDual[v]
      ek$workDual[v] <- 0
      ek$workCost[v] <- ek$workCost[v] + shift
      ek$costs_shifted <- TRUE
    }
  }
}

.hs_dual_phase1 <- function(ek, ds) {
  ds$rebuild_reason <- 0L
  if (.hs_bailout(ek)) return(invisible())
  .hs_init_bound(ek, "dual", 1L)
  .hs_init_nonbasic(ek)
  if (is.null(ek$bt)) .hs_put_backtracking(ek)
  .hs_dual_loop(ek, ds, FALSE)
  if (ds$solve_phase %in% c(0L, -5L, -6L)) return(invisible())
  if (ek$solve_bailout) return(invisible())
  if (ds$row_out < 0L) {
    if (ek$dual_objective_value == 0) ds$solve_phase <- 2L else .hs_dual_assess_phase1(ek, ds)
  } else if (ds$rebuild_reason %in% c(8L, 9L)) {
    ds$solve_phase <- -5L
  } else if (ds$variable_in < 0L) {
    if (ek$costs_perturbed) {
      .hs_dual_cleanup(ek, ds)
      if (ds$dualInfeasCount == 0L) ds$solve_phase <- 2L
    } else {
      ds$solve_phase <- -5L
    }
  }
  if (ds$solve_phase %in% c(2L, -1L, -5L)) {
    .hs_init_bound(ek, "dual", 2L)
    .hs_init_nonbasic(ek)
    if (ds$solve_phase == 2L && ek$dual_phase1_cleanup_level < 2L) ek$allow_cost_perturbation <- TRUE
  }
}

.hs_dual_phase2 <- function(ek, ds) {
  ds$rebuild_reason <- 0L
  ds$solve_phase <- 2L
  ek$solve_bailout <- FALSE
  if (.hs_bailout(ek)) return(invisible())
  tot <- ek$n + ek$m
  ds$free_list <- which(ek$nonbasicFlag == 1L & ek$workLower == -Inf & ek$workUpper == Inf)
  if (is.null(ek$bt)) .hs_put_backtracking(ek)
  .hs_dual_loop(ek, ds, TRUE)
  if (ds$solve_phase %in% c(0L, -5L, -6L)) return(invisible())
  if (ek$solve_bailout) return(invisible())
  if (ds$dualInfeasCount > 0L) {
    ds$solve_phase <- 1L
  } else if (ds$row_out < 0L) {
    .hs_dual_cleanup(ek, ds)
    if (ds$dualInfeasCount > 0L) {
      ds$solve_phase <- -3L
    } else {
      ds$solve_phase <- -2L
      ek$model_status <- "Optimal"
    }
  } else if (ds$rebuild_reason %in% c(8L, 9L)) {
    ds$solve_phase <- -5L
  }
}

# HEkkDual::solve
.hs_dual_solve <- function(ek, force_phase2 = FALSE) {
  ds <- new.env(parent = emptyenv())
  ds$rebuild_reason <- 0L
  ds$row_out <- -1L
  ds$variable_in <- -1L
  ds$free_list <- integer()
  ds$dualInfeasCount <- 0L
  ek$model_status <- "Notset"
  ek$solve_bailout <- FALSE
  ek$exit_algorithm <- "dual"
  ek$bt <- NULL; ek$backtracking <- FALSE
  ds$edge_weight_mode <- if (identical(ek$dual_edge_weight_strategy, "steepest")) "dse" else ek$dual_edge_weight_strategy
  # possiblyUseLiDualSteepestEdge
  nz <- ek$A != 0
  cnt <- colSums(nz)
  ek$store_squared <- !(ek$n > 0 && all(cnt <= 24) && all(abs(ek$A[nz]) == 1) && sum(cnt) / ek$n <= 6)
  .hs_init_cost(ek, "dual", FALSE)
  .hs_compute_dual(ek)
  .hs_simplex_dual_infeasible(ek)
  dual_feasible <- ek$num_dual_infeasibilities == 0L
  ds$force_phase2 <- force_phase2 || ek$max_dual_infeasibility^2 < ek$dtol
  near_optimal <- (dual_feasible || ds$force_phase2) && ek$num_primal_infeasibilities < 1000L &&
    ek$max_primal_infeasibility < 1e-3
  perturb <- !near_optimal
  .hs_init_cost(ek, "dual", perturb)
  if (!isTRUE(ek$has_dse_weights)) {
    ek$dual_edge_weight <- rep(1, ek$m)
    if (ds$edge_weight_mode == "dse") {
      if (all(ek$basicIndex > ek$n)) {
        ek$has_dse_weights <- TRUE
      } else if (near_optimal) {
        ds$edge_weight_mode <- "devex"
      } else {
        .hs_compute_dse_weights(ek)
        ek$has_dse_weights <- TRUE
      }
    }
    if (ds$edge_weight_mode == "devex") .hs_dual_init_devex(ek)
  }
  if (perturb) {
    .hs_compute_dual(ek)
    .hs_simplex_dual_infeasible(ek)
    ds$dualInfeasCount <- ek$num_dual_infeasibilities
  }
  ds$solve_phase <- if (ds$force_phase2) 2L else if (ds$dualInfeasCount > 0L) 1L else 2L
  while (ds$solve_phase > 0L || ds$solve_phase == 0L) {
    if (ds$solve_phase == 0L) {
      .hs_init_bound(ek, "dual", 0L)
      .hs_init_nonbasic(ek)
      .hs_simplex_dual_infeasible(ek)
      ds$dualInfeasCount <- ek$num_dual_infeasibilities
      ds$solve_phase <- if (ds$dualInfeasCount > 0L) 1L else 2L
      if (isTRUE(ek$backtracking)) {
        .hs_init_bound(ek, "dual", ds$solve_phase)
        .hs_init_nonbasic(ek)
        ek$backtracking <- FALSE
      }
    }
    if (ds$solve_phase == 1L) .hs_dual_phase1(ek, ds) else .hs_dual_phase2(ek, ds)
    if (ek$solve_bailout) return(.hs_return_from_solve(ek))
    if (ds$solve_phase == -6L) { ek$model_status <- "Unknown"; return(.hs_return_from_solve(ek)) }
    if (ds$solve_phase == -5L) { ek$model_status <- "SolveError"; return(.hs_return_from_solve(ek)) }
    if (ds$solve_phase %in% c(-1L, -3L, -4L)) break
  }
  if (ds$solve_phase %in% c(-3L, -4L)) {
    ek$dual_cleanup_level <- ek$dual_cleanup_level + 1L
    if (ek$dual_cleanup_level > 1L) {
      ek$model_status <- if (ds$solve_phase == -3L) "Optimal" else "Infeasible"
    } else {
      save <- ek$primal_bound_perturbation_multiplier
      ek$primal_bound_perturbation_multiplier <- 0
      .hs_primal_solve(ek, TRUE)
      ek$primal_bound_perturbation_multiplier <- save
      return(invisible())
    }
  }
  .hs_return_from_solve(ek)
}

# HEkk::returnFromSolve
.hs_return_from_solve <- function(ek) {
  st <- ek$model_status
  if (st == "Infeasible") {
    if (ek$exit_algorithm == "primal") {
      .hs_init_cost(ek, "dual", FALSE)
      .hs_compute_dual(ek)
    }
  } else if (st == "UnboundedOrInfeasible") {
    .hs_init_bound(ek, "dual", 2L)
    .hs_compute_primal(ek)
  } else if (st %in% c("IterationLimit", "Unknown")) {
    .hs_init_bound(ek, "dual", 2L)
    .hs_init_nonbasic(ek)
    .hs_compute_primal(ek)
    .hs_init_cost(ek, "dual", FALSE)
    .hs_compute_dual(ek)
  }
  .hs_simplex_primal_infeasible(ek)
  .hs_simplex_dual_infeasible(ek)
  .hs_compute_primal_objective(ek)
  invisible()
}
