# scipy.optimize.linprog and scipy.optimize.milp as OptimizationSciPy.jl
# calls them. The HiGHS methods run the native port of the HiGHS LP
# presolve, postsolve and simplex solver (R/optimization-scipy-highs-*.R);
# "highs-ipm" runs IPX with crossover and the simplex clean up of HiGHS. milp branches and bounds
# depth first over node LPs of that solver, without the presolve, cuts and
# heuristics of the HiGHS MIP solver, so node counts differ. The legacy
# methods use the dense two phase simplex method with Bland's rule below.

# min c'z subject to A z = b, z >= 0. Returns z, the duals y (the derivative
# of the optimum with respect to b), the status of linprog (0 optimal,
# 1 iteration limit, 2 infeasible, 3 unbounded) and the pivot count.
.sp_lp_simplex <- function(A, b, cost, maxiter = Inf, tol = 1e-9) {
  m <- nrow(A); n <- ncol(A)
  neg <- b < 0
  A[neg, ] <- -A[neg, ]
  b[neg] <- -b[neg]
  Tb <- cbind(A, diag(1, m), b)
  rhs <- n + m + 1L
  basis <- n + seq_len(m)
  nit <- 0L
  run <- function(cst, allowed) {
    repeat {
      cb <- cst[basis]
      red <- cst - colSums(Tb[, -rhs, drop = FALSE] * cb)
      cand <- which(red < -tol & allowed)
      if (!length(cand)) return(0L)
      if (nit >= maxiter) return(1L)
      k <- cand[1L]
      col <- Tb[, k]
      pos <- which(col > tol)
      if (!length(pos)) return(3L)
      ratio <- Tb[pos, rhs] / col[pos]
      ties <- pos[ratio <= min(ratio) + tol]
      r <- ties[which.min(basis[ties])]
      Tb[r, ] <<- Tb[r, ] / Tb[r, k]
      for (i in seq_len(nrow(Tb))) if (i != r && Tb[i, k] != 0) Tb[i, ] <<- Tb[i, ] - Tb[i, k] * Tb[r, ]
      basis[r] <<- k
      nit <<- nit + 1L
    }
  }
  st <- run(c(numeric(n), rep(1, m)), rep(TRUE, n + m))
  if (st == 1L) return(list(status = 1L, nit = nit))
  if (sum(Tb[basis > n, rhs]) > 1e-7 * max(1, sum(b))) return(list(status = 2L, nit = nit))
  # pivot the artificial variables left at zero out of the basis
  keep <- rep(TRUE, m)
  for (r in which(basis > n)) {
    k <- which(abs(Tb[r, seq_len(n)]) > tol)
    if (!length(k)) { keep[r] <- FALSE; next }
    k <- k[1L]
    Tb[r, ] <- Tb[r, ] / Tb[r, k]
    for (i in seq_len(m)) if (i != r && Tb[i, k] != 0) Tb[i, ] <- Tb[i, ] - Tb[i, k] * Tb[r, ]
    basis[r] <- k
  }
  cst <- c(cost, numeric(m))
  allowed <- c(rep(TRUE, n), rep(FALSE, m))
  Tb <- Tb[keep, , drop = FALSE]
  basis <- basis[keep]
  st <- run(cst, allowed)
  z <- numeric(n)
  z[basis] <- Tb[, rhs]
  y <- unname(colSums(Tb[, n + seq_len(m), drop = FALSE] * cst[basis]))
  y[neg] <- -y[neg]
  list(z = z, y = y, status = st, nit = nit)
}

# linprog's problem min c'x, A_ub x <= b_ub, A_eq x = b_eq, lb <= x <= ub in
# the standard form of .sp_lp_simplex, with the marginals of HiGHS
.sp_linprog <- function(cvec, A_ub, b_ub, A_eq, b_eq, lb, ub, maxiter = Inf) {
  n <- length(cvec)
  A_ub <- matrix(as.numeric(A_ub %||% numeric()), ncol = n)
  A_eq <- matrix(as.numeric(A_eq %||% numeric()), ncol = n)
  b_ub <- as.numeric(b_ub %||% numeric()); b_eq <- as.numeric(b_eq %||% numeric())
  mu <- nrow(A_ub); me <- nrow(A_eq)
  # x = off + M z with z >= 0
  off <- numeric(n)
  M <- matrix(0, n, 0)
  up <- integer()
  for (j in seq_len(n)) {
    e <- replace(numeric(n), j, 1)
    if (is.finite(lb[j])) {
      off[j] <- lb[j]; M <- cbind(M, e)
      if (is.finite(ub[j])) up <- c(up, j)
    } else if (is.finite(ub[j])) {
      off[j] <- ub[j]; M <- cbind(M, -e)
    } else {
      M <- cbind(M, e, -e)
    }
  }
  nz <- ncol(M)
  mup <- length(up)
  U <- matrix(0, mup, nz)
  for (i in seq_len(mup)) U[i, which(M[up[i], ] != 0)] <- 1
  As <- rbind(cbind(A_ub %*% M, diag(1, mu), matrix(0, mu, mup)),
              cbind(A_eq %*% M, matrix(0, me, mu + mup)),
              cbind(U, matrix(0, mup, mu), diag(1, mup)))
  bs <- c(b_ub - as.numeric(A_ub %*% off), b_eq - as.numeric(A_eq %*% off), ub[up] - lb[up])
  cs <- c(as.numeric(crossprod(M, cvec)), numeric(mu + mup))
  dimnames(As) <- NULL
  r <- .sp_lp_simplex(As, bs, cs, maxiter)
  if (r$status != 0L) return(list(status = r$status, nit = r$nit))
  x <- off + as.numeric(M %*% r$z[seq_len(nz)])
  m_ub <- r$y[seq_len(mu)]
  m_eq <- r$y[mu + seq_len(me)]
  m_up <- numeric(n)
  m_up[up] <- r$y[mu + me + seq_len(mup)]
  red <- cvec - as.numeric(crossprod(A_ub, m_ub)) - as.numeric(crossprod(A_eq, m_eq)) - m_up
  m_lo <- ifelse(is.finite(lb), red, 0)
  only_ub <- !is.finite(lb) & is.finite(ub)
  m_up[only_ub] <- red[only_ub]
  slack <- b_ub - as.numeric(A_ub %*% x)
  con <- b_eq - as.numeric(A_eq %*% x)
  list(x = x, fun = sum(cvec * x), slack = slack, con = con, status = 0L, nit = r$nit,
       ineqlin = list(residual = slack, marginals = m_ub), eqlin = list(residual = con, marginals = m_eq),
       lower = list(residual = x - lb, marginals = m_lo), upper = list(residual = ub - x, marginals = m_up))
}

.sp_highs_messages <- c(
  "0" = "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
  "1" = "Iteration limit reached. (HiGHS Status 14: model_status is Iteration limit reached; primal_status is Feasible)",
  "2" = "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
  "3" = "The problem is unbounded. (HiGHS Status 10: model_status is Unbounded; primal_status is Feasible)",
  "4" = "The problem is unbounded or infeasible. (HiGHS Status 9: model_status is Primal infeasible or unbounded; primal_status is None)")

.sp_legacy_lp_messages <- c(
  "0" = "Optimization terminated successfully.", "1" = "Iteration limit reached.",
  "2" = "Optimization failed. Unable to find a feasible starting point.",
  "3" = "Optimization failed. The problem appears to be unbounded.")

.sp_pad_bounds <- function(v, n, fill) {
  if (is.null(v)) return(rep(fill, n))
  v <- as.numeric(v)
  if (length(v) != n) v <- c(v, rep(fill, n))[seq_len(n)]
  v
}

# the result of linprog from a solve of .hs_lp_solve, as _linprog_highs
# assembles it from _highs_wrapper
.sp_highs_result <- function(sol, n, nub, b_ub, b_eq, lo, hi) {
  code <- c(Optimal = 7L, Infeasible = 8L, UnboundedOrInfeasible = 9L, Unbounded = 10L,
            IterationLimit = 14L, SolveError = 4L, Unknown = 15L)[[sol$status]]
  text <- c(Optimal = "Optimal", Infeasible = "Infeasible", UnboundedOrInfeasible = "Primal infeasible or unbounded",
            Unbounded = "Unbounded", IterationLimit = "Iteration limit reached", SolveError = "Solve error",
            Unknown = "Unknown")[[sol$status]]
  status <- c(Optimal = 0L, Infeasible = 2L, UnboundedOrInfeasible = 4L, Unbounded = 3L,
              IterationLimit = 1L, SolveError = 4L, Unknown = 4L)[[sol$status]]
  lead <- c(Optimal = "Optimization terminated successfully. ", Infeasible = "The problem is infeasible. ",
            UnboundedOrInfeasible = "The problem is unbounded or infeasible. ", Unbounded = "The problem is unbounded. ",
            IterationLimit = "Iteration limit reached. ", SolveError = "", Unknown = "The HiGHS status code was not recognized. ")[[sol$status]]
  if (status == 0L) {
    msg <- paste0(lead, "(HiGHS Status 7: Optimal)")
  } else {
    pstat <- if (sol$status %in% c("Infeasible", "UnboundedOrInfeasible")) "None"
             else if (isTRUE(sol$primal_feasible)) "Feasible" else "Infeasible"
    msg <- paste0(lead, "(HiGHS Status ", code, ": model_status is ", text, "; primal_status is ", pstat, ")")
  }
  res <- list(status = status, message = msg, success = status == 0L, nit = sol$nit)
  if (status != 0L) return(c(list(x = rep(NaN, n), fun = NaN), res))
  x <- sol$col_value
  rhs <- c(b_ub, b_eq)
  slack_all <- rhs - sol$row_value
  ml <- ifelse(sol$col_status == "lower", sol$col_dual, 0)
  mu <- ifelse(sol$col_status == "upper", sol$col_dual, 0)
  iu <- seq_len(nub); ie <- nub + seq_along(b_eq)
  c(list(x = x, fun = sol$objective, slack = slack_all[iu], con = slack_all[ie]), res,
    list(ineqlin = list(residual = slack_all[iu], marginals = sol$row_dual[iu]),
         eqlin = list(residual = slack_all[ie], marginals = sol$row_dual[ie]),
         lower = list(residual = x - lo, marginals = ml), upper = list(residual = hi - x, marginals = mu)))
}

# ScipyLinprog: the cost vector is the objective at u0; bounds default to
# linprog's x >= 0 when the problem has none
.scipy_run_linprog <- function(alg, f, x0, lb, ub, opts, kw) {
  cvec <- as.numeric(f(x0))
  n <- length(x0)
  if (is.null(lb) && is.null(ub)) {
    lo <- rep(0, n); hi <- rep(Inf, n)
  } else {
    lo <- .sp_pad_bounds(lb, n, -Inf); hi <- .sp_pad_bounds(ub, n, Inf)
  }
  if (is.null(kw$A_ub) != is.null(kw$b_ub)) stop("Both A_ub and b_ub must be provided together", call. = FALSE)
  if (is.null(kw$A_eq) != is.null(kw$b_eq)) stop("Both A_eq and b_eq must be provided together", call. = FALSE)
  maxiter <- if (is.null(opts$maxiters)) Inf else as.integer(opts$maxiters)
  if (startsWith(alg$method, "highs")) {
    A_ub <- matrix(as.numeric(kw$A_ub %||% numeric()), ncol = n)
    A_eq <- matrix(as.numeric(kw$A_eq %||% numeric()), ncol = n)
    b_ub <- as.numeric(kw$b_ub %||% numeric()); b_eq <- as.numeric(kw$b_eq %||% numeric())
    run <- if (alg$method == "highs-ipm") .hs_highs_run_ipm else .hs_highs_run
    sol <- run(cvec, rbind(A_ub, A_eq), c(rep(-Inf, length(b_ub)), b_eq), c(b_ub, b_eq), lo, hi,
               list(iteration_limit = maxiter, presolve = kw$presolve, ptol = kw$primal_feasibility_tolerance,
                    dtol = kw$dual_feasibility_tolerance, run_crossover = kw$run_crossover,
                    edge_weight = switch(kw$simplex_dual_edge_weight_strategy %||% "",
                                         dantzig = "dantzig", devex = "devex", steepest = "steepest", "dse")))
    return(.sp_highs_result(sol, n, length(b_ub), b_ub, b_eq, lo, hi))
  }
  res <- .sp_linprog(cvec, kw$A_ub, kw$b_ub, kw$A_eq, kw$b_eq, lo, hi, maxiter)
  res$message <- .sp_legacy_lp_messages[[as.character(res$status)]]
  res$success <- res$status == 0L
  if (!res$success) {
    res$x <- rep(NaN, n)
    res$fun <- NaN
  }
  res
}

# ScipyMilp: bounds default to free variables; integrality 0 continuous,
# 1 integer, 2 semicontinuous, 3 semi-integer
.scipy_run_milp <- function(alg, f, x0, lb, ub, a, kw) {
  cvec <- as.numeric(f(x0))
  n <- length(cvec)
  lo <- .sp_pad_bounds(lb, n, -Inf); hi <- .sp_pad_bounds(ub, n, Inf)
  typ <- rep_len(as.integer(kw$integrality %||% 0L), n)
  given <- c(!is.null(kw$A), !is.null(kw$lb_con), !is.null(kw$ub_con))
  if (any(given) && !all(given))
    stop("A, lb_con, and ub_con must all be provided for linear constraints", call. = FALSE)
  A <- matrix(0, 0, n); lc <- uc <- numeric()
  if (all(given)) {
    A <- matrix(as.numeric(kw$A), ncol = n)
    lc <- rep_len(as.numeric(kw$lb_con), nrow(A)); uc <- rep_len(as.numeric(kw$ub_con), nrow(A))
  }
  # internal switch: the HiGHS MIP solver port (R/optimization-scipy-highs-mip-*.R)
  # replaces the depth first search once it matches every probe MILP
  if (isTRUE(getOption("modelcompR.highs_mip", FALSE)))
    return(.sp_milp_highs(cvec, A, lc, uc, lo, hi, typ, kw))
  best <- NULL
  nodes <- 0L
  status <- 2L
  stack <- list(list(lo = lo, hi = hi, typ = typ))
  while (length(stack)) {
    nd <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL
    nodes <- nodes + 1L
    semi <- nd$typ >= 2L
    rlo <- ifelse(semi, pmin(0, nd$lo), nd$lo)
    sol <- .hs_lp_solve(cvec, A, lc, uc, rlo, nd$hi)
    if (sol$status %in% c("Unbounded", "UnboundedOrInfeasible")) {
      status <- if (any(typ != 0L)) 4L else 3L; best <- NULL; break
    }
    if (sol$status != "Optimal") next
    r <- list(x = sol$col_value, fun = sol$objective)
    if (!is.null(best) && r$fun >= best$fun - 1e-9 * max(1, abs(best$fun))) next
    x <- r$x
    frac <- nd$typ %in% c(1L, 3L) & abs(x - round(x)) > 1e-6
    gap <- semi & abs(x) > 1e-6 & x < nd$lo - 1e-6
    j <- which(frac | gap)[1L]
    if (is.na(j)) { best <- r; status <- 0L; next }
    up <- nd; down <- nd
    if (gap[j]) {
      down$lo[j] <- 0; down$hi[j] <- 0; down$typ[j] <- 0L
      up$typ[j] <- nd$typ[j] - 2L
    } else {
      down$hi[j] <- floor(x[j]); up$lo[j] <- ceiling(x[j])
      if (semi[j]) up$typ[j] <- 1L
    }
    stack[[length(stack) + 1L]] <- up
    stack[[length(stack) + 1L]] <- down
  }
  res <- list(status = status, success = status == 0L, message = .sp_highs_messages[[as.character(status)]])
  if (status == 0L) {
    res <- c(list(x = best$x, fun = best$fun), res,
             list(mip_node_count = nodes, mip_dual_bound = best$fun, mip_gap = 0))
  } else {
    res$x <- rep(NaN, n)
    res$fun <- NaN
  }
  res$retcode <- if (res$success) "Success" else "Failure"
  res
}

# scipy.optimize.milp through Highs::callSolveMip (the HiGHS MIP solver port)
# with the result checks of _highs_wrapper and SciPy's status and message
# conversion (_highs_to_scipy_status_message)
.sp_milp_highs <- function(cvec, A, lc, uc, lo, hi, typ, kw) {
  n <- length(cvec)
  lp <- .hmip_HighsLp_from_dense(cvec, A, lc, uc, lo, hi, integrality = typ)
  opts <- list()
  if (!is.null(kw$node_limit)) opts$mip_max_nodes <- kw$node_limit
  if (!is.null(kw$mip_rel_gap)) opts$mip_rel_gap <- kw$mip_rel_gap
  if (!is.null(kw$time_limit)) opts$time_limit <- kw$time_limit
  if (isFALSE(kw$presolve)) opts$presolve <- "off"
  r <- .hmip_callSolveMip(lp, opts)
  st <- r$model_status
  code <- c(Notset = 0L, SolveError = 4L, Optimal = 7L, Infeasible = 8L, UnboundedOrInfeasible = 9L,
            Unbounded = 10L, ObjectiveTarget = 12L, TimeLimit = 13L, IterationLimit = 14L,
            Unknown = 15L, SolutionLimit = 16L, Interrupt = 17L, MemoryLimit = 18L)[[st]]
  text <- c(Notset = "Not Set", SolveError = "Solve error", Optimal = "Optimal", Infeasible = "Infeasible",
            UnboundedOrInfeasible = "Primal infeasible or unbounded", Unbounded = "Unbounded",
            ObjectiveTarget = "Objective target reached", TimeLimit = "Time limit reached",
            IterationLimit = "Iteration limit reached", Unknown = "Unknown",
            SolutionLimit = "Solution limit reached", Interrupt = "Interrupted by user",
            MemoryLimit = "Memory limit reached")[[st]]
  lead <- c(Optimal = "Optimization terminated successfully. ", TimeLimit = "Time limit reached. ",
            IterationLimit = "Iteration limit reached. ", Infeasible = "The problem is infeasible. ",
            Unbounded = "The problem is unbounded. ",
            UnboundedOrInfeasible = "The problem is unbounded or infeasible. ")
  status <- c(Optimal = 0L, TimeLimit = 1L, IterationLimit = 1L, Infeasible = 2L, Unbounded = 3L)[st]
  if (is.na(status)) status <- 4L
  ld <- if (st %in% names(lead)) lead[[st]]
        else if (st %in% c("Notset", "SolveError", "ObjectiveTarget")) ""
        else "The HiGHS status code was not recognized. "
  have_sol <- r$objective != Inf
  fail <- !(st %in% c("Optimal", "TimeLimit", "IterationLimit", "SolutionLimit")) ||
    (st != "Optimal" && !have_sol)
  hmsg <- if (fail) paste0("model_status is ", text, "; primal_status is ",
                           if (have_sol) "Feasible" else "None") else text
  res <- list(status = unname(status), success = status == 0L,
              message = paste0(ld, "(HiGHS Status ", code, ": ", hmsg, ")"))
  if (!fail) {
    res <- c(list(x = r$col_value, fun = r$objective), res,
             list(mip_node_count = r$mip_node_count, mip_dual_bound = r$mip_dual_bound, mip_gap = r$mip_gap))
  } else {
    res$x <- rep(NaN, n)
    res$fun <- NaN
  }
  res$retcode <- if (res$success) "Success" else "Failure"
  res
}
