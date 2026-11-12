# scipy.optimize.linprog and scipy.optimize.milp as OptimizationSciPy.jl
# calls them. SciPy solves both with the HiGHS library (C++); here a dense
# two phase simplex method with Bland's rule solves the linear programs and
# a depth first branch and bound solves the mixed integer programs. At a
# unique optimum the solutions, the objective and the marginals agree with
# HiGHS; iteration and node counts differ.

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
  res <- .sp_linprog(cvec, kw$A_ub, kw$b_ub, kw$A_eq, kw$b_eq, lo, hi, maxiter)
  msgs <- if (startsWith(alg$method, "highs")) .sp_highs_messages else .sp_legacy_lp_messages
  res$message <- msgs[[as.character(res$status)]]
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
  A_ub <- A_eq <- matrix(0, 0, n)
  b_ub <- b_eq <- numeric()
  if (all(given)) {
    A <- matrix(as.numeric(kw$A), ncol = n)
    lc <- rep_len(as.numeric(kw$lb_con), nrow(A)); uc <- rep_len(as.numeric(kw$ub_con), nrow(A))
    for (i in seq_len(nrow(A))) {
      if (lc[i] == uc[i]) { A_eq <- rbind(A_eq, A[i, ]); b_eq <- c(b_eq, uc[i]); next }
      if (is.finite(uc[i])) { A_ub <- rbind(A_ub, A[i, ]); b_ub <- c(b_ub, uc[i]) }
      if (is.finite(lc[i])) { A_ub <- rbind(A_ub, -A[i, ]); b_ub <- c(b_ub, -lc[i]) }
    }
  }
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
    r <- .sp_linprog(cvec, A_ub, b_ub, A_eq, b_eq, rlo, nd$hi)
    if (r$status == 3L) { status <- if (any(typ != 0L)) 4L else 3L; best <- NULL; break }
    if (r$status != 0L) next
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
