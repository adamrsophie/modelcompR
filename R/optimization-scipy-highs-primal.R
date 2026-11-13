# The primal simplex method of HiGHS 1.12.0 (simplex/HEkkPrimal.cpp) with
# Devex pricing, as SciPy 1.18.1 vendors it; HiGHS calls it to clean up
# after the dual simplex method and to settle an LP that dual phase 1 finds
# unbounded or infeasible. MIT license notice in
# R/optimization-scipy-highs-util.R. Phases as in the dual file.

.hs_primal_init_devex <- function(ek, ps) {
  tot <- ek$n + ek$m
  ps$edge_weight <- rep(1, tot)
  ps$devex_index <- ek$nonbasicFlag * ek$nonbasicFlag
  ps$num_bad_devex <- 0L
}

.hs_shift_bound <- function(ek, lower, value, random_value, bound) {
  feas <- (1 + random_value) * ek$ptol
  if (lower) {
    shift <- bound - value + feas
    list(bound = bound - shift, shift = shift)
  } else {
    shift <- value - bound + feas
    list(bound = bound + shift, shift = shift)
  }
}

# HEkkPrimal::correctPrimal (primal correction strategy always)
.hs_primal_correct <- function(ek) {
  for (i in seq_len(ek$m)) {
    v <- ek$baseValue[i]
    viol <- if (v < ek$baseLower[i] - ek$ptol) -1L else if (v > ek$baseUpper[i] + ek$ptol) 1L else 0L
    if (!viol) next
    if (!ek$allow_bound_perturbation) next
    j <- ek$basicIndex[i]
    if (viol > 0) {
      s <- .hs_shift_bound(ek, FALSE, v, ek$numTotRandomValue[j], ek$workUpper[j])
      ek$workUpper[j] <- s$bound; ek$baseUpper[i] <- s$bound
      ek$workUpperShift[j] <- ek$workUpperShift[j] + s$shift
    } else {
      s <- .hs_shift_bound(ek, TRUE, v, ek$numTotRandomValue[j], ek$workLower[j])
      ek$workLower[j] <- s$bound; ek$baseLower[i] <- s$bound
      ek$workLowerShift[j] <- ek$workLowerShift[j] + s$shift
    }
    ek$bounds_perturbed <- TRUE
  }
}

# HEkkPrimal::phase1ComputeDual
.hs_primal_phase1_dual <- function(ek) {
  n <- ek$n; tot <- n + ek$m
  buf <- numeric(ek$m)
  base <- 5e-7
  for (i in seq_len(ek$m)) {
    v <- ek$baseValue[i]
    viol <- if (v < ek$baseLower[i] - ek$ptol) -1 else if (v > ek$baseUpper[i] + ek$ptol) 1 else 0
    if (!viol) next
    buf[i] <- viol * (1 + base * ek$numTotRandomValue[i])
  }
  ek$workCost <- numeric(tot)
  ek$workDual <- numeric(tot)
  ek$workCost[ek$basicIndex] <- buf
  y <- .hs_btran(ek, buf, "dual_col")
  row <- .hs_price_col(ek, y)
  ek$workDual <- -ek$nonbasicFlag * c(row, y)
}

.hs_primal_rebuild <- function(ek, ps) {
  .hs_clear_taboo(ek)
  if (.hs_rebuild_refactor(ek, ps$rebuild_reason) && !.hs_reinvert(ek)) {
    ps$solve_phase <- -5L
    return(invisible())
  }
  ps$rebuild_reason <- 0L
  if (isTRUE(ek$backtracking)) {
    ps$solve_phase <- 0L
    return(invisible())
  }
  .hs_compute_primal(ek)
  if (ps$solve_phase == 2L) .hs_primal_correct(ek)
  inf <- .hs_primal_infeas(ek$baseValue, ek$baseLower, ek$baseUpper, ek$ptol)
  ek$num_primal_infeasibilities <- sum(inf > ek$ptol)
  if (ek$num_primal_infeasibilities > 0L) {
    ps$solve_phase <- 1L
    .hs_primal_phase1_dual(ek)
  } else {
    if (ps$solve_phase == 1L) {
      .hs_init_cost(ek, "primal")
      ps$solve_phase <- 2L
    }
    .hs_compute_dual(ek)
  }
  .hs_simplex_dual_infeasible(ek)
  .hs_compute_primal_objective(ek)
  ps$num_flip_since_rebuild <- 0L
  ek$has_fresh_rebuild <- TRUE
  .hs_reset_clock(ek)
}

.hs_primal_chuzc <- function(ek, ps) {
  dual <- ek$workDual
  for (b in ek$bad) if (b$taboo) dual[b$vin] <- 0
  best <- 0; vin <- -1L
  w <- ps$edge_weight
  for (v in ps$free_set) {
    di <- abs(dual[v])
    if (di > ek$dtol && di * di > best * w[v]) { vin <- v; best <- di * di / w[v] }
  }
  for (v in seq_len(ek$n + ek$m)) {
    di <- -ek$nonbasicMove[v] * dual[v]
    if (di > ek$dtol && di * di > best * w[v]) { vin <- v; best <- di * di / w[v] }
  }
  ps$variable_in <- vin
}

.hs_primal_use_variable_in <- function(ek, ps) {
  vin <- ps$variable_in
  updated <- ek$workDual[vin]
  ps$move_in <- if (updated > 0) -1L else 1L
  hv <- .hs_hv_collect(ek, list(array = numeric(ek$m), index = integer(), count = 0L), vin, 1)
  ps$col_aq <- .hs_ftran(ek, hv, "col_aq", pack = TRUE)
  ps$col_aq_tick <- ek$last_tick
  computed <- ek$workCost[vin]
  for (i in ek$last_index) computed <- computed - ps$col_aq[i] * ek$workCost[ek$basicIndex[i]]
  ek$workDual[vin] <- computed
  ps$theta_dual <- computed
  small <- abs(computed) <= ek$dtol
  sign_error <- updated * computed <= 0
  if (small) ek$num_dual_infeasibilities <- ek$num_dual_infeasibilities - 1L
  if (small || sign_error) {
    if (!small && ek$update_count > 0L) ps$rebuild_reason <- 7L
    return(FALSE)
  }
  TRUE
}

.hs_primal_phase1_choose_row <- function(ek, ps) {
  tol <- if (ek$update_count < 10L) 1e-9 else if (ek$update_count < 20L) 1e-8 else 1e-7
  pt <- ek$ptol
  col <- ps$col_aq
  rt <- numeric(); ri <- integer(); tt <- numeric(); ti <- integer()
  for (i in which(col != 0)) {
    a <- col[i] * ps$move_in
    v <- ek$baseValue[i]; lo <- ek$baseLower[i]; up <- ek$baseUpper[i]
    if (a > tol) {
      if (v > up + pt) {
        th <- (v - up - pt) / a
        rt <- c(rt, th); ri <- c(ri, i); tt <- c(tt, th); ti <- c(ti, i)
      }
      if (v > lo - pt && lo > -Inf) {
        rt <- c(rt, (v - lo + pt) / a); ri <- c(ri, i - ek$m)
        tt <- c(tt, (v - lo) / a); ti <- c(ti, i - ek$m)
      }
    }
    if (a < -tol) {
      if (v < lo - pt) {
        th <- (v - lo + pt) / a
        rt <- c(rt, th); ri <- c(ri, i - ek$m); tt <- c(tt, th); ti <- c(ti, i - ek$m)
      }
      if (v < up + pt && up < Inf) {
        rt <- c(rt, (v - up - pt) / a); ri <- c(ri, i)
        tt <- c(tt, (v - up) / a); ti <- c(ti, i)
      }
    }
  }
  ps$row_out <- -1L
  ps$move_out <- 0L
  if (!length(rt)) return(invisible())
  o <- order(rt, ri); rt <- rt[o]; ri <- ri[o]
  max_theta <- rt[1L]
  gradient <- abs(ps$theta_dual)
  for (k in seq_along(rt)) {
    row <- if (ri[k] > 0L) ri[k] else ri[k] + ek$m
    gradient <- gradient - abs(col[row])
    if (gradient <= 0) break
    max_theta <- rt[k]
  }
  o <- order(tt, ti); tt <- tt[o]; ti <- ti[o]
  max_alpha <- 0
  last <- length(tt)
  for (k in seq_along(tt)) {
    row <- if (ti[k] > 0L) ti[k] else ti[k] + ek$m
    if (tt[k] > max_theta) { last <- k - 1L; break }
    max_alpha <- max(max_alpha, abs(col[row]))
  }
  if (last >= 1L) for (k in last:1L) {
    row <- if (ti[k] > 0L) ti[k] else ti[k] + ek$m
    if (abs(col[row]) > max_alpha * 0.1) {
      ps$row_out <- row
      ps$move_out <- if (ti[k] > 0L) 1L else -1L
      break
    }
  }
}

.hs_primal_choose_row <- function(ek, ps) {
  tol <- if (ek$update_count < 10L) 1e-9 else if (ek$update_count < 20L) 1e-8 else 1e-7
  col <- ps$col_aq
  relax <- 0x1.249ad2594c37dp+332
  nz <- which(col != 0)
  for (i in nz) {
    a <- col[i] * ps$move_in
    if (a > tol) {
      sp <- ek$baseValue[i] - ek$baseLower[i] + ek$ptol
      if (sp < relax * a) relax <- sp / a
    } else if (a < -tol) {
      sp <- ek$baseValue[i] - ek$baseUpper[i] - ek$ptol
      if (sp > relax * a) relax <- sp / a
    }
  }
  best <- 0
  ps$row_out <- -1L
  for (i in nz) {
    a <- col[i] * ps$move_in
    if (a > tol) {
      if (ek$baseValue[i] - ek$baseLower[i] < relax * a && best < a) { best <- a; ps$row_out <- i }
    } else if (a < -tol) {
      if (ek$baseValue[i] - ek$baseUpper[i] > relax * a && best < -a) { best <- -a; ps$row_out <- i }
    }
  }
}

.hs_primal_bound_swap <- function(ek, ps) {
  vin <- ps$variable_in
  if (ps$row_out < 0L) {
    theta <- ps$move_in * Inf
    ps$move_out <- 0L
  } else {
    r <- ps$row_out
    ps$alpha_col <- ps$col_aq[r]
    if (ps$solve_phase == 2L) ps$move_out <- if (ps$alpha_col * ps$move_in > 0) -1L else 1L
    theta <- if (ps$move_out == 1L) (ek$baseValue[r] - ek$baseUpper[r]) / ps$alpha_col
             else (ek$baseValue[r] - ek$baseLower[r]) / ps$alpha_col
  }
  flipped <- FALSE
  lo <- ek$workLower[vin]; up <- ek$workUpper[vin]
  ps$value_in <- ek$workValue[vin] + theta
  if (ps$move_in > 0) {
    if (ps$value_in > up + ek$ptol) {
      flipped <- TRUE; ps$row_out <- -1L; ps$value_in <- up; theta <- up - lo
    }
  } else if (ps$value_in < lo - ek$ptol) {
    flipped <- TRUE; ps$row_out <- -1L; ps$value_in <- lo; theta <- lo - up
  }
  ps$theta_primal <- theta
  if (ps$solve_phase == 2L && !(ps$row_out >= 1L || flipped)) ps$rebuild_reason <- 3L
}

.hs_primal_update <- function(ek, ps) {
  n <- ek$n
  vin <- ps$variable_in
  flipped <- ps$row_out < 1L
  base <- 5e-7
  if (flipped) {
    ps$variable_out <- vin
    ek$workValue[vin] <- ps$value_in
    ek$nonbasicMove[vin] <- -ps$move_in
  } else if (ek$bounds_perturbed) {
    # adjustPerturbedEquationOut
    vout <- ps$variable_out
    b <- .hs_lp_bounds(ek)
    if (b$lower[vout] >= b$upper[vout]) {
      fixed <- b$lower[vout]
      ps$theta_primal <- (ek$baseValue[ps$row_out] - fixed) / ps$alpha_col
      ek$workLower[vout] <- fixed; ek$workUpper[vout] <- fixed
      ps$value_in <- ek$workValue[vin] + ps$theta_primal
    }
  }
  col <- ps$col_aq
  if (ps$solve_phase == 1L) {
    # phase1UpdatePrimal
    change <- numeric(ek$m)
    for (i in which(col != 0)) {
      ek$baseValue[i] <- ek$baseValue[i] - ps$theta_primal * col[i]
      j <- ek$basicIndex[i]
      was <- ek$workCost[j]
      v <- ek$baseValue[i]
      viol <- if (v < ek$baseLower[i] - ek$ptol) -1 else if (v > ek$baseUpper[i] + ek$ptol) 1 else 0
      cost <- viol
      if (viol) cost <- cost * (1 + base * ek$numTotRandomValue[i])
      ek$workCost[j] <- cost
      if (was) { if (!cost) ek$num_primal_infeasibilities <- ek$num_primal_infeasibilities - 1L }
      else if (cost) ek$num_primal_infeasibilities <- ek$num_primal_infeasibilities + 1L
      dc <- cost - was
      if (dc) {
        change[i] <- dc
        if (j > n) ek$workDual[j] <- ek$workDual[j] + dc
      }
    }
    # basicFeasibilityChangeUpdateDual
    if (any(change != 0)) {
      y <- .hs_btran(ek, change, "col_basic_feasibility_change")
      pr <- .hs_price(ek, y, ek$last_index, ek$dens$row_basic_feasibility_change)
      ek$workDual[seq_len(n)] <- ek$workDual[seq_len(n)] - pr
      ek$workDual[n + seq_len(ek$m)] <- ek$workDual[n + seq_len(ek$m)] - y
    }
  } else {
    # phase2UpdatePrimal with bound shifting
    for (i in seq_len(ek$m)) {
      if (col[i] != 0) ek$baseValue[i] <- ek$baseValue[i] - ps$theta_primal * col[i]
      v <- ek$baseValue[i]
      viol <- if (v < ek$baseLower[i] - ek$ptol) -1L else if (v > ek$baseUpper[i] + ek$ptol) 1L else 0L
      if (!viol) next
      j <- ek$basicIndex[i]
      if (viol > 0) {
        s <- .hs_shift_bound(ek, FALSE, v, ek$numTotRandomValue[j], ek$workUpper[j])
        ek$workUpper[j] <- s$bound; ek$baseUpper[i] <- s$bound
        ek$workUpperShift[j] <- ek$workUpperShift[j] + s$shift
      } else {
        s <- .hs_shift_bound(ek, TRUE, v, ek$numTotRandomValue[j], ek$workLower[j])
        ek$workLower[j] <- s$bound; ek$baseLower[i] <- s$bound
        ek$workLowerShift[j] <- ek$workLowerShift[j] + s$shift
      }
    }
  }
  if (flipped) {
    ps$num_flip_since_rebuild <- ps$num_flip_since_rebuild + 1L
    ek$total_tick <- (ek$total_tick %||% 0) + ps$col_aq_tick
    return(invisible())
  }
  r <- ps$row_out
  ek$baseValue[r] <- ps$value_in
  # considerInfeasibleValueIn
  lo <- ek$workLower[vin]; up <- ek$workUpper[vin]
  viol <- if (ps$value_in < lo - ek$ptol) -1 else if (ps$value_in > up + ek$ptol) 1 else 0
  if (viol) {
    if (ps$solve_phase == 1L) {
      ek$num_primal_infeasibilities <- ek$num_primal_infeasibilities + 1L
      cost <- viol * (1 + base * ek$numTotRandomValue[r])
      ek$workCost[vin] <- cost
      ek$workDual[vin] <- ek$workDual[vin] + cost
    } else {
      if (viol > 0) {
        s <- .hs_shift_bound(ek, FALSE, ps$value_in, ek$numTotRandomValue[vin], up)
        ek$workUpper[vin] <- s$bound; ek$workUpperShift[vin] <- ek$workUpperShift[vin] + s$shift
      } else {
        s <- .hs_shift_bound(ek, TRUE, ps$value_in, ek$numTotRandomValue[vin], lo)
        ek$workLower[vin] <- s$bound; ek$workLowerShift[vin] <- ek$workLowerShift[vin] + s$shift
      }
      ek$bounds_perturbed <- TRUE
    }
  }
  # updateDual
  vout <- ps$variable_out
  theta_dual <- ek$workDual[vin] / ps$alpha_col
  ek$workDual[seq_len(n)] <- ek$workDual[seq_len(n)] - theta_dual * ps$row_ap
  ek$workDual[n + seq_len(ek$m)] <- ek$workDual[n + seq_len(ek$m)] - theta_dual * ps$row_ep
  ek$workDual[vin] <- 0
  ek$workDual[vout] <- -theta_dual
  # updateDevex
  w <- ps$edge_weight; di <- ps$devex_index
  pw <- 0
  for (i in which(col != 0)) {
    a <- di[ek$basicIndex[i]] * col[i]
    pw <- pw + a * a
  }
  pw <- pw + di[vin]
  if (w[vin] > 3 * pw) ps$num_bad_devex <- ps$num_bad_devex + 1L
  pw <- pw / (col[r] * col[r])
  for (j in which(ps$row_ap != 0)) {
    dv <- pw * ps$row_ap[j]^2 + di[j]
    if (w[j] < dv) w[j] <- dv
  }
  for (i in which(ps$row_ep != 0)) {
    j <- n + i
    dv <- pw * ps$row_ep[i]^2 + di[j]
    if (w[j] < dv) w[j] <- dv
  }
  w[vout] <- max(1, pw)
  w[vin] <- 1
  ps$edge_weight <- w
  if (ek$nonbasicMove[vin] == 0L) ps$free_set <- setdiff(ps$free_set, vin)
  if (isTRUE(ek$has_dse_weights)) {
    dse <- .hs_ftran(ek, list(array = ps$row_ep, index = ps$row_ep_index, count = length(ps$row_ep_index)), "row_DSE")
    dw <- ek$dual_edge_weight
    dw[r] <- .hs_dsum(ps$row_ep[ps$row_ep_index]^2)
    npw <- dw[r] / (ps$alpha_col^2)
    kai <- -2 / ps$alpha_col
    for (i in which(col != 0)) dw[i] <- max(1e-4, dw[i] + col[i] * (npw * col[i] + kai * dse[i]))
    dw[r] <- npw
    ek$dual_edge_weight <- dw
  }
  .hs_update_pivots(ek, vin, r, ps$move_out, col)
  ps$rebuild_reason <- .hs_update_hint(ek, ps$rebuild_reason)
  ek$iteration_count <- ek$iteration_count + 1L
  if (ps$num_bad_devex > 3L) .hs_primal_init_devex(ek, ps)
  ek$total_tick <- (ek$total_tick %||% 0) + ps$col_aq_tick + ps$row_ep_tick
}

.hs_primal_iterate <- function(ek, ps) {
  ps$row_out <- -2L
  .hs_primal_chuzc(ek, ps)
  if (ps$variable_in < 0L) {
    ps$rebuild_reason <- 1L
    return(invisible())
  }
  if (!.hs_primal_use_variable_in(ek, ps)) return(invisible())
  if (ps$solve_phase == 1L) {
    .hs_primal_phase1_choose_row(ek, ps)
    if (ps$row_out < 1L) {
      ps$solve_phase <- -5L
      return(invisible())
    }
  } else {
    .hs_primal_choose_row(ek, ps)
  }
  .hs_primal_bound_swap(ek, ps)
  if (ps$rebuild_reason == 3L) return(invisible())
  if (ps$row_out >= 1L) {
    r <- ps$row_out
    ps$alpha_col <- ps$col_aq[r]
    ps$variable_out <- ek$basicIndex[r]
    e <- numeric(ek$m); e[r] <- 1
    ps$row_ep <- .hs_btran(ek, e, "row_ep", pack = TRUE)
    ps$row_ep_index <- ek$last_index
    ps$row_ep_tick <- ek$last_tick
    ps$row_ap <- .hs_price(ek, ps$row_ep, ps$row_ep_index)
    vin <- ps$variable_in
    alpha_row <- if (vin <= ek$n) ps$row_ap[vin] else ps$row_ep[vin - ek$n]
    trouble <- abs(abs(ps$alpha_col) - abs(alpha_row)) / min(abs(ps$alpha_col), abs(alpha_row))
    if (trouble > 1e-7 && ek$update_count > 0L) {
      ps$rebuild_reason <- 7L
      return(invisible())
    }
  }
  if (.hs_is_bad_basis_change(ek, ps$variable_in, ps$row_out, ps$rebuild_reason)) return(invisible())
  .hs_primal_update(ek, ps)
  if (ek$num_primal_infeasibilities == 0L && ps$solve_phase == 1L) ps$rebuild_reason <- 2L
}

.hs_primal_cleanup <- function(ek, ps) {
  if (!ek$bounds_perturbed) return(invisible())
  .hs_init_bound(ek, "primal", ps$solve_phase, FALSE)
  .hs_init_nonbasic(ek)
  ek$allow_bound_perturbation <- FALSE
  .hs_compute_primal(ek)
  .hs_simplex_primal_infeasible(ek)
  .hs_compute_primal_objective(ek)
  .hs_simplex_dual_infeasible(ek)
}

.hs_primal_loop <- function(ek, ps, phase) {
  if (is.null(ek$bt)) .hs_put_backtracking(ek)
  repeat {
    .hs_primal_rebuild(ek, ps)
    if (ps$solve_phase %in% c(0L, -5L)) return(invisible())
    if (.hs_bailout(ek)) return(invisible())
    if (ps$solve_phase != phase) break
    repeat {
      .hs_primal_iterate(ek, ps)
      if (.hs_bailout(ek)) return(invisible())
      if (ps$solve_phase == -5L) return(invisible())
      if (ps$rebuild_reason != 0L) break
    }
    finished <- ek$has_fresh_rebuild && ps$num_flip_since_rebuild == 0L && !.hs_rebuild_refactor(ek, ps$rebuild_reason)
    if (finished && .hs_taboo(ek)) {
      ps$solve_phase <- -6L
      return(invisible())
    }
    if (finished) break
  }
}

# HEkkPrimal::solve
.hs_primal_solve <- function(ek, force_phase2 = FALSE) {
  ps <- new.env(parent = emptyenv())
  ps$rebuild_reason <- 0L
  ek$model_status <- "Notset"
  ek$solve_bailout <- FALSE
  ek$exit_algorithm <- "primal"
  ek$bt <- NULL; ek$backtracking <- FALSE
  if (!isTRUE(ek$has_dse_weights)) ek$dual_edge_weight <- rep(1, ek$m)
  .hs_primal_init_devex(ek, ps)
  ps$free_set <- which(ek$nonbasicFlag == 1L & ek$workLower == -Inf & ek$workUpper == Inf)
  feasible <- ek$num_primal_infeasibilities == 0L
  force <- force_phase2 || ek$max_primal_infeasibility^2 < ek$ptol
  near_optimal <- ek$num_dual_infeasibilities < 1000L && ek$max_dual_infeasibility < 1e-3 && (feasible || force)
  if (!near_optimal && ek$primal_bound_perturbation_multiplier != 0) {
    .hs_init_bound(ek, "primal", 0L, TRUE)
    .hs_init_nonbasic(ek)
    .hs_compute_primal(ek)
    .hs_simplex_primal_infeasible(ek)
  }
  if (.hs_bailout(ek)) return(.hs_return_from_solve(ek))
  ps$solve_phase <- if (ek$num_primal_infeasibilities > 0L) 1L else 2L
  if (force) ps$solve_phase <- 2L
  repeat {
    if (ps$solve_phase == 0L) {
      .hs_simplex_primal_infeasible(ek)
      ps$solve_phase <- if (ek$num_primal_infeasibilities > 0L) 1L else 2L
      if (isTRUE(ek$backtracking)) {
        .hs_init_cost(ek, "primal")
        .hs_init_nonbasic(ek)
        ek$backtracking <- FALSE
      }
    }
    if (ps$solve_phase == 1L) {
      .hs_primal_loop(ek, ps, 1L)
      if (!ek$solve_bailout && ps$solve_phase == 1L && ps$variable_in < 0L) {
        if (ek$bounds_perturbed) {
          .hs_primal_cleanup(ek, ps)
        } else {
          ek$model_status <- "Infeasible"
          ps$solve_phase <- -1L
        }
      }
    } else if (ps$solve_phase == 2L) {
      .hs_primal_loop(ek, ps, 2L)
      if (!ek$solve_bailout && ps$solve_phase == 2L) {
        if (ps$variable_in < 0L) {
          .hs_primal_cleanup(ek, ps)
          if (ek$num_primal_infeasibilities > 0L) {
            ps$solve_phase <- -3L
          } else {
            ps$solve_phase <- -2L
            ek$model_status <- "Optimal"
          }
        } else if (ek$bounds_perturbed) {
          .hs_primal_cleanup(ek, ps)
          if (ek$num_primal_infeasibilities > 0L) ps$solve_phase <- 1L
        } else {
          ps$solve_phase <- -1L
          ek$model_status <- "Unbounded"
        }
      }
    }
    if (ek$solve_bailout) return(.hs_return_from_solve(ek))
    if (ps$solve_phase == -6L) { ek$model_status <- "Unknown"; return(.hs_return_from_solve(ek)) }
    if (ps$solve_phase == -5L) { ek$model_status <- "SolveError"; return(.hs_return_from_solve(ek)) }
    if (ps$solve_phase %in% c(-1L, -2L, -3L)) break
  }
  if (ps$solve_phase == -3L) {
    save <- ek$dual_cost_perturbation_multiplier
    ek$dual_cost_perturbation_multiplier <- 0
    .hs_dual_solve(ek, TRUE)
    ek$dual_cost_perturbation_multiplier <- save
    return(invisible())
  }
  .hs_return_from_solve(ek)
}
