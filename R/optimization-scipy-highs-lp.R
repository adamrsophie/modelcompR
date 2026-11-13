# The LP driver of HiGHS 1.12.0 for linprog: solveLpSimplex (simplex/HApp.h),
# the equilibration scaling of scaleLp (lp_data/HighsLpUtils.cpp),
# solveUnconstrainedLp (lp_data/HighsSolve.cpp) and the solution and basis
# that HEkk::getSolution and HEkk::getHighsBasis hand to SciPy's
# _highs_wrapper. MIT license notice in R/optimization-scipy-highs-util.R.

# equilibrationScaleMatrix with the choice of scaleLp (strategy 2)
.hs_scale <- function(A, cost) {
  m <- nrow(A); n <- ncol(A)
  nz <- A != 0
  av <- abs(A[nz])
  if (!length(av) || (min(av) >= 0.2 && max(av) <= 5)) return(NULL)
  big <- 0x1.4e718d7d7625ap+664
  nzc <- cost[cost != 0]
  include_cost <- length(nzc) && min(abs(nzc)) < 0.1
  lo_s <- 2^-20; hi_s <- 2^20
  cs <- rep(1, n); rs <- rep(1, m)
  for (pass in 1:6) {
    rmin <- rep(big, m); rmax <- rep(1 / big, m)
    for (j in seq_len(n)) {
      cmin <- big; cmax <- 1 / big
      if (include_cost && cost[j] != 0) { cmin <- min(cmin, abs(cost[j])); cmax <- max(cmax, abs(cost[j])) }
      rows <- which(nz[, j])
      for (i in rows) { v <- abs(A[i, j]) * rs[i]; cmin <- min(cmin, v); cmax <- max(cmax, v) }
      cs[j] <- min(max(lo_s, 1 / sqrt(cmin * cmax)), hi_s)
      for (i in rows) { v <- abs(A[i, j]) * cs[j]; rmin[i] <- min(rmin[i], v); rmax[i] <- max(rmax[i], v) }
    }
    rs <- pmin(pmax(lo_s, 1 / sqrt(rmin * rmax)), hi_s)
  }
  cs <- 2^floor(log(cs) / log(2) + 0.5)
  rs <- 2^floor(log(rs) / log(2) + 0.5)
  S <- A * outer(rs, cs)
  equil <- function(M, margin) {
    vapply(seq_len(dim(M)[margin]), function(k) {
      v <- abs(if (margin == 1L) M[k, ] else M[, k])
      v <- v[v != 0]
      if (!length(v)) 1 else 1 / sqrt(min(v) * max(v))
    }, 0)
  }
  oc <- equil(A, 2L); orr <- equil(A, 1L); nc <- equil(S, 2L); nr <- equil(S, 1L)
  gm <- function(x) { g <- exp(sum(log(x)) / length(x)); max(g, 1 / g) }
  mean_imp <- sqrt((gm(oc) * gm(orr)) / (gm(nc) * gm(nr)))
  ext_imp <- (max(oc) / min(oc) + max(orr) / min(orr)) / (max(nc) / min(nc) + max(nr) / min(nr))
  sv <- abs(S[nz])
  ratio_imp <- (max(av) / min(av)) / (max(sv) / min(sv))
  if (ext_imp * mean_imp * ratio_imp < 1) return(NULL)
  list(col = cs, row = rs)
}

# solveUnconstrainedLp
.hs_unconstrained <- function(cost, col_lower, col_upper, ptol, dtol) {
  n <- length(cost)
  x <- numeric(n); stat <- character(n)
  npi <- 0L; ndi <- 0L
  for (j in seq_len(n)) {
    d <- cost[j]; lo <- col_lower[j]; up <- col_upper[j]
    pinf <- 0; dinf <- 0
    if (lo > up) {
      if (lo == Inf) {
        if (up == -Inf) { x[j] <- 0; stat[j] <- "zero"; pinf <- Inf; dinf <- abs(d) }
        else { x[j] <- up; stat[j] <- "upper"; pinf <- lo - up; dinf <- max(d, 0) }
      } else { x[j] <- lo; stat[j] <- "lower"; pinf <- lo - up; dinf <- max(-d, 0) }
    } else if (lo == -Inf && up == Inf) {
      x[j] <- 0; stat[j] <- "zero"; dinf <- abs(d)
    } else if (d >= dtol) {
      if (lo > -Inf) { x[j] <- lo; stat[j] <- "lower" } else { x[j] <- up; stat[j] <- "upper"; dinf <- d }
    } else if (d <= -dtol) {
      if (up < Inf) { x[j] <- up; stat[j] <- "upper" } else { x[j] <- lo; stat[j] <- "lower"; dinf <- -d }
    } else {
      if (lo == -Inf) { x[j] <- up; stat[j] <- "upper" } else { x[j] <- lo; stat[j] <- "lower" }
      dinf <- abs(d)
    }
    if (pinf > ptol) npi <- npi + 1L
    if (dinf > dtol) ndi <- ndi + 1L
  }
  status <- if (npi) "Infeasible" else if (ndi) "Unbounded" else "Optimal"
  list(status = status, col_value = x, col_dual = cost, row_value = numeric(), row_dual = numeric(),
       col_status = stat, objective = sum(x * cost), nit = 0L, primal_feasible = npi == 0L)
}

.hs_new_ekk <- function(cost, A, row_lower, row_upper, col_lower, col_upper, opts) {
  ek <- new.env(parent = emptyenv())
  ek$n <- length(cost); ek$m <- nrow(A)
  ek$A <- A
  ek$lp <- list(col_cost = cost, col_lower = col_lower, col_upper = col_upper,
                row_lower = row_lower, row_upper = row_upper)
  ek$ptol <- opts$ptol %||% 1e-7
  ek$dtol <- opts$dtol %||% 1e-7
  ek$iteration_limit <- opts$iteration_limit %||% Inf
  ek$dual_edge_weight_strategy <- opts$edge_weight %||% "dse"
  ek$start_basis <- opts$basis
  ek$dual_cost_perturbation_multiplier <- opts$perturb %||% 1
  ek$primal_bound_perturbation_multiplier <- 1
  ek$update_limit <- 5000L
  ek$iteration_count <- as.integer(opts$iteration_start %||% 0L)
  ek$previous_cycling <- -2147483647L
  ek$bad <- list()
  ek$visited <- new.env(parent = emptyenv())
  ek$dual_phase1_cleanup_level <- 0L
  ek$dual_cleanup_level <- 0L
  ek$rng <- .hs_rng(0)
  # HEkk::initialiseEkk
  .hs_random_vectors(ek)
  ek
}

# HEkk::solve: initialiseForSolve, then the dual simplex method (or the
# primal one when strategy is "primal")
.hs_ekk_solve <- function(ek, force_phase2 = FALSE, strategy = "dual") {
  ek$dual_phase1_cleanup_level <- 0L
  ek$dual_cleanup_level <- 0L
  ek$previous_cycling <- -2147483647L
  .hs_init_densities(ek)
  if (is.null(ek$basicIndex)) {
    if (is.null(ek$start_basis)) .hs_slack_basis(ek) else .hs_set_basis(ek, ek$start_basis)
  }
  .hs_random_vectors(ek)
  .hs_init_cost(ek, "primal")
  .hs_init_bound(ek, "primal", 0L)
  .hs_init_nonbasic(ek)
  .hs_compute_primal(ek)
  .hs_compute_dual(ek)
  .hs_simplex_primal_infeasible(ek)
  .hs_simplex_dual_infeasible(ek)
  .hs_compute_dual_objective(ek)
  .hs_compute_primal_objective(ek)
  ek$visited <- new.env(parent = emptyenv())
  ek$visited[[.hs_basis_key(ek)]] <- TRUE
  ek$model_status <- "Notset"
  ek$solve_bailout <- FALSE
  if (ek$num_primal_infeasibilities == 0L && ek$num_dual_infeasibilities == 0L) {
    ek$model_status <- "Optimal"
    return(invisible())
  }
  ek$allow_cost_perturbation <- TRUE
  ek$allow_bound_perturbation <- TRUE
  if (strategy == "primal") {
    .hs_primal_solve(ek, force_phase2)
  } else {
    .hs_dual_solve(ek, force_phase2)
    if (ek$model_status == "UnboundedOrInfeasible") .hs_primal_solve(ek)
  }
}

# HEkk::getSolution and HEkk::getHighsBasis in the space of the ekk LP
.hs_get_solution <- function(ek) {
  n <- ek$n; m <- ek$m
  wv <- ek$workValue; wd <- ek$workDual
  wv[ek$basicIndex] <- ek$baseValue
  wd[ek$basicIndex] <- 0
  stat <- character(n)
  for (j in seq_len(n)) {
    stat[j] <- if (ek$nonbasicFlag[j] == 0L) "basic"
    else if (ek$nonbasicMove[j] == 1L) "lower"
    else if (ek$nonbasicMove[j] == -1L) "upper"
    else if (ek$lp$col_lower[j] == ek$lp$col_upper[j]) { if (wd[j] >= 0) "lower" else "upper" }
    else "zero"
  }
  rstat <- character(m)
  for (i in seq_len(m)) {
    k <- n + i
    rstat[i] <- if (ek$nonbasicFlag[k] == 0L) "basic"
    else if (ek$nonbasicMove[k] == 1L) "upper"
    else if (ek$nonbasicMove[k] == -1L) "lower"
    else if (ek$lp$row_lower[i] == ek$lp$row_upper[i]) { if (wd[k] >= 0) "lower" else "upper" }
    else "zero"
  }
  list(col_value = wv[seq_len(n)], col_dual = wd[seq_len(n)], row_value = -wv[n + seq_len(m)],
       row_dual = -wd[n + seq_len(m)], col_status = stat, row_status = rstat)
}

# HEkk::setBasis from a HighsBasis (statuses "basic", "lower", "upper", "zero")
.hs_set_basis <- function(ek, basis) {
  n <- ek$n; m <- ek$m
  st <- c(basis$col_status, basis$row_status)
  b <- .hs_lp_bounds(ek)
  ek$basicIndex <- which(st == "basic")
  ek$nonbasicFlag <- as.integer(st != "basic")
  mv <- integer(n + m)
  for (k in which(st != "basic")) {
    up_dir <- if (k <= n) 1L else -1L
    mv[k] <- if (b$lower[k] == b$upper[k]) 0L
             else if (st[k] == "lower") up_dir else if (st[k] == "upper") -up_dir else 0L
  }
  ek$nonbasicMove <- mv
  .hs_invert(ek)
}

# solveLpSimplex for an LP with rows; returns the model status, the
# unscaled solution and basis, the objective, the iteration count and the
# primal solution status
.hs_lp_solve <- function(cost, A, row_lower, row_upper, col_lower, col_upper, opts = list()) {
  ptol <- opts$ptol %||% 1e-7; dtol <- opts$dtol %||% 1e-7
  if (nrow(A) == 0L) return(.hs_unconstrained(cost, col_lower, col_upper, ptol, dtol))
  sc <- if (isFALSE(opts$scale)) NULL else .hs_scale(A, cost)
  if (is.null(sc)) {
    ek <- .hs_new_ekk(cost, A, row_lower, row_upper, col_lower, col_upper, opts)
    .hs_ekk_solve(ek, FALSE, opts$strategy %||% "dual")
    sol <- .hs_get_solution(ek)
  } else {
    As <- A * outer(sc$row, sc$col)
    ek <- .hs_new_ekk(cost * sc$col, As, row_lower * sc$row, row_upper * sc$row,
                      col_lower / sc$col, col_upper / sc$col, opts)
    .hs_ekk_solve(ek, FALSE, opts$strategy %||% "dual")
    sol <- .hs_get_solution(ek)
    status <- ek$model_status
    # getUnscaledInfeasibilities
    nb <- which(ek$nonbasicFlag == 1L & ek$workLower != ek$workUpper)
    mu <- c(1 / sc$col, sc$row)
    ud <- ek$workDual[nb] * mu[nb]
    free <- ek$workLower[nb] == -Inf & ek$workUpper[nb] == Inf
    dinf <- ifelse(free, abs(ud), -ek$nonbasicMove[nb] * ud)
    ndi <- sum(dinf >= dtol)
    bmu <- c(sc$col, 1 / sc$row)[ek$basicIndex]
    npi <- sum(.hs_primal_infeas(ek$baseValue * bmu, ek$baseLower * bmu, ek$baseUpper * bmu, ptol) > 0)
    refine <- (status == "Optimal" && (npi + ndi) > 0) ||
      status %in% c("Infeasible", "UnboundedOrInfeasible", "Unbounded", "Unknown")
    if (!refine) {
      sol$col_value <- sol$col_value * sc$col
      sol$col_dual <- sol$col_dual / sc$col
      sol$row_value <- sol$row_value / sc$row
      sol$row_dual <- sol$row_dual * sc$row
      sol$status <- status
      sol$objective <- ek$primal_objective_value
      sol$nit <- ek$iteration_count
      sol$primal_feasible <- npi == 0L
      return(sol)
    }
    # moveLp: the unscaled LP into the ekk with the basis kept
    ek$A <- A
    ek$lp <- list(col_cost = cost, col_lower = col_lower, col_upper = col_upper,
                  row_lower = row_lower, row_upper = row_upper)
    .hs_invert(ek)
    solve_unscaled <- TRUE
    if (status == "Infeasible" && !is.null(ek$dual_ray)) {
      e <- numeric(ek$m); e[ek$dual_ray$row] <- 1
      if (.hs_proof_of_infeasibility(ek, .hs_btran(ek, e, "row_ep", pack = TRUE), ek$dual_ray$sign)) solve_unscaled <- FALSE
    }
    if (solve_unscaled) {
      strategy <- if (npi == 0L) "primal" else "dual"
      if (strategy == "dual" && !isTRUE(ek$has_dse_weights)) ek$dual_edge_weight_strategy <- "devex"
      .hs_ekk_solve(ek, TRUE, strategy)
      sol <- .hs_get_solution(ek)
    } else {
      sol$col_value <- sol$col_value * sc$col
      sol$col_dual <- sol$col_dual / sc$col
      sol$row_value <- sol$row_value / sc$row
      sol$row_dual <- sol$row_dual * sc$row
    }
  }
  sol$status <- ek$model_status
  sol$objective <- ek$primal_objective_value
  sol$nit <- ek$iteration_count
  sol$primal_feasible <- ek$num_primal_infeasibilities == 0L
  sol
}
