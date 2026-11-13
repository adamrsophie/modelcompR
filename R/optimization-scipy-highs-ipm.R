# Highs::run for an LP with solver = "ipm" (HiGHS 1.12.0 lp_data/Highs.cpp
# 1466-1800, as SciPy 1.18.1 linprog(method = "highs-ipm") runs it):
# presolve, solveLp with IPX (lp_data/HighsSolve.cpp), the simplex clean up
# of an unwelcome IPX status (kUnknown, or kUnboundedOrInfeasible since
# allow_unbounded_or_infeasible is false) when run_crossover is not "off",
# postsolve, the simplex hot start of the original LP from the postsolved
# basis, and the objective value that lpKktCheck stores
# (HighsLp::objectiveValue). SciPy reports nit as simplex_nit or ipm_nit.
# MIT license notice in R/optimization-scipy-highs-util.R.

.hs_ipx_status <- c(kOptimal = "Optimal", kUnknown = "Unknown", kInfeasible = "Infeasible",
                    kUnboundedOrInfeasible = "UnboundedOrInfeasible",
                    kIterationLimit = "IterationLimit", kTimeLimit = "TimeLimit",
                    kSolveError = "SolveError")

# solveLp with solver = "ipm" on one LP; info accumulates the IPM and
# crossover iteration counts as HighsInfo does
.hs_ipm_solve_lp <- function(cost, A, row_lower, row_upper, col_lower, col_upper, opts, info) {
  if (nrow(A) == 0L || all(A == 0))
    return(.hs_unconstrained(cost, col_lower, col_upper, opts$ptol %||% 1e-7, opts$dtol %||% 1e-7))
  crossover <- opts$run_crossover %||% "on"
  ipx_opts <- list(run_crossover = crossover, primal_feasibility_tolerance = opts$ptol,
                   dual_feasibility_tolerance = opts$dtol)
  if (is.finite(opts$iteration_limit %||% Inf)) ipx_opts$ipm_iteration_limit <- opts$iteration_limit
  r <- .highs_ipx_solve(list(c = cost, A = A, row_lower = row_lower, row_upper = row_upper,
                             col_lower = col_lower, col_upper = col_upper, offset = 0, sense = 1L),
                        ipx_opts)
  info$ipm_nit <- info$ipm_nit + r$ipm_iterations
  info$crossover_nit <- info$crossover_nit + r$crossover_iterations
  st <- unname(.hs_ipx_status[r$status])
  basis <- !is.null(r$col_basis)
  sol <- list(status = st, col_value = r$x, row_value = r$row_activity, col_dual = r$col_dual,
              row_dual = r$row_dual, objective = r$objective, nit = 0L, basis_valid = basis,
              primal_feasible = st == "Optimal")
  if (basis) {
    sol$col_status <- tolower(sub("^k", "", r$col_basis))
    sol$row_status <- tolower(sub("^k", "", r$row_basis))
  }
  if (st %in% c("Unknown", "UnboundedOrInfeasible") && crossover != "off") {
    o <- opts
    if (basis) o$basis <- list(col_status = sol$col_status, row_status = sol$row_status)
    sol <- .hs_lp_solve(cost, A, row_lower, row_upper, col_lower, col_upper, o)
    sol$basis_valid <- TRUE
  }
  sol
}

# calculateRowValuesQuad on a dense matrix: in column order, each rounded
# product added to a HighsCDouble row sum
.hs_row_values_quad <- function(A, x) {
  acc <- matrix(0, nrow(A), 2L)
  for (j in seq_len(ncol(A))) for (i in which(A[, j] != 0))
    acc[i, ] <- .cd_add(acc[i, ], x[j] * A[i, j])
  acc[, 1L] + acc[, 2L]
}

.hs_highs_run_ipm <- function(cost, A, row_lower, row_upper, col_lower, col_upper, opts = list()) {
  info <- new.env(parent = emptyenv())
  info$ipm_nit <- 0L
  info$crossover_nit <- 0L
  done <- function(sol) {
    force(sol)  # the solve must run before the counts are read
    sol$ipm_nit <- info$ipm_nit
    sol$crossover_nit <- info$crossover_nit
    sol$simplex_nit <- sol$nit %||% 0L
    sol$nit <- if (sol$simplex_nit != 0L) sol$simplex_nit else sol$ipm_nit
    if (!is.null(sol$col_value) && length(sol$col_value) == length(cost)) {
      v <- 0
      for (j in seq_along(cost)) v <- v + cost[j] * sol$col_value[j]
      sol$objective <- v
    }
    sol
  }
  solve_orig <- function()
    .hs_ipm_solve_lp(cost, A, row_lower, row_upper, col_lower, col_upper, opts, info)
  if (isFALSE(opts$presolve) || nrow(A) == 0L || all(A == 0)) return(done(solve_orig()))
  # lp_no_solution_basis (IPX without crossover) lets presolve use the rules
  # without basis postsolve (lp_presolve_requires_basis_postsolve = false)
  pre <- .hs_presolve(cost, A, row_lower, row_upper, col_lower, col_upper,
                      opts$ptol %||% 1e-7, opts$dtol %||% 1e-7,
                      basis_postsolve = !identical(opts$run_crossover, "off"))
  switch(pre$status,
    notreduced = return(done(solve_orig())),
    infeasible = return(done(list(status = "Infeasible", nit = 0L, primal_feasible = FALSE))),
    unboundedorinfeasible = {
      o <- opts
      o$strategy <- "primal"
      return(done(.hs_lp_solve(cost, A, row_lower, row_upper, col_lower, col_upper, o)))
    })
  if (pre$status == "empty") {
    red <- list(status = "Optimal", nit = 0L, col_value = numeric(), col_dual = numeric(),
                row_value = numeric(), row_dual = numeric(), col_status = character(),
                row_status = character(), basis_valid = TRUE)
  } else {
    lp <- pre$lp
    red <- .hs_ipm_solve_lp(lp$cost, lp$A, lp$row_lower, lp$row_upper, lp$col_lower, lp$col_upper,
                            opts, info)
    if (!red$status %in% c("Optimal", "Unknown"))
      return(done(list(status = red$status, nit = red$nit, primal_feasible = FALSE)))
    if (is.null(red$col_status)) red$col_status <- character(length(red$col_value))
    if (is.null(red$row_status)) red$row_status <- character(length(red$row_value))
  }
  post <- .hs_postsolve(pre, red, opts$dtol %||% 1e-7)
  # lp_no_solution_basis (IPX without crossover) invalidates the basis
  # before postsolve, so there is no simplex hot start
  if (!isTRUE(red$basis_valid) || identical(opts$run_crossover, "off")) {
    # IPX without crossover: the postsolved primal dual point and no basis
    # runPostsolve recomputes the row activities (calculateRowValuesQuad)
    return(done(list(status = red$status, nit = red$nit, col_value = post$col_value,
                     row_value = .hs_row_values_quad(A, post$col_value), col_dual = post$col_dual,
                     row_dual = post$row_dual, col_status = rep("", length(cost)),
                     primal_feasible = red$status == "Optimal")))
  }
  o <- opts
  o$basis <- list(col_status = post$col_status, row_status = post$row_status)
  o$iteration_start <- red$nit
  done(.hs_lp_solve(cost, A, row_lower, row_upper, col_lower, col_upper, o))
}
