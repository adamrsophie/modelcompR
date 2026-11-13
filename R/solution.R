# Solutions, after SciMLBase's ODESolution and SymbolicIndexingInterface.
#
# A solution is a closure of class c("ODESolution", "SciMLSolution",
# "function"): calling it interpolates, sol(2.5) or sol(c(1, 2), idxs = x).
# Its data live in the closure environment and are read with `$`:
#   sol$t        saved times
#   sol$u        list of saved state vectors (as Julia's sol.u)
#   sol$retcode  "Success", "Terminated", "MaxIters", ...
#   sol$prob, sol$alg, sol$stats, sol$ps
# `sol[x]` returns the time series of a variable or expression.

.new_solution <- function(cls, data) {
  env <- list2env(data, parent = emptyenv())
  env$obs_cache <- new.env(parent = emptyenv())
  fn <- function(t, idxs = NULL, deriv = 0L) .sol_interp(env, t, idxs, deriv)
  environment(fn) <- list2env(list(env = env), parent = asNamespace("modelcompR"))
  class(fn) <- c(cls, "SciMLSolution", "function")
  fn
}

.sol_env <- function(sol) get("env", envir = environment(unclass(sol)))

# A copy of a solution with some fields replaced: a new data environment, so
# the original keeps its values (the SciMLBase solution builders return new
# solutions).
.sol_copy <- function(sol, ...) {
  data <- as.list(.sol_env(sol), all.names = TRUE)
  data$obs_cache <- NULL
  new <- list(...)
  for (nm in names(new)) data[nm] <- list(new[[nm]])
  .new_solution(setdiff(class(sol), c("SciMLSolution", "function")), data)
}

.new_ode_solution <- function(t, u, prob, alg, retcode, steps = NULL, stats = list(), p = NULL,
                              events = list(), dense = FALSE) {
  .new_solution("ODESolution", list(t = t, U = u, prob = prob, alg = alg, retcode = retcode,
                                    steps = steps, stats = stats, p_final = p %||% prob$p,
                                    events = events, dense = dense, sys = prob$sys))
}

#' @export
`$.SciMLSolution` <- function(x, name) {
  env <- .sol_env(x)
  switch(name,
    t = env$t,
    u = {
      # save_idxs: the saved components (symbolic indexing keeps them all)
      U <- if (is.null(env$save_idxs)) env$U else env$U[, env$save_idxs, drop = FALSE]
      d <- if (is.null(env$save_idxs)) env$u_dim
      lapply(seq_len(nrow(U)), function(i) if (is.null(d)) U[i, ] else array(U[i, ], d))
    },
    retcode = env$retcode,
    prob = env$prob,
    alg = env$alg,
    stats = env$stats,
    sys = env$sys,
    ps = .sol_param_accessor(env),
    dense = env$dense,
    events = env$events,
    objective = env$objective,
    resid = env$resid,
    original = env$original,
    get0(name, envir = env, inherits = FALSE))
}

#' @export
names.SciMLSolution <- function(x) c("t", "u", "retcode", "prob", "alg", "stats", "ps")

.sol_param_accessor <- function(env) {
  prob <- env$prob
  prob$p <- env$p_final
  .param_accessor(prob)
}

#' @export
length.SciMLSolution <- function(x) length(.sol_env(x)$t)

#' @export
print.ODESolution <- function(x, ...) {
  env <- .sol_env(x)
  cat("ODESolution with retcode ", env$retcode, " (", env$alg$name %||% "", ")\n", sep = "")
  n <- length(env$t)
  cat("Saved time points: ", n, if (n) paste0(" in [", format(env$t[1L]), ", ", format(env$t[n]), "]"), "\n", sep = "")
  if (n) {
    cat("Final state:\n")
    print(env$U[n, ])
  }
  invisible(x)
}

#' @export
as.matrix.SciMLSolution <- function(x, ...) {
  env <- .sol_env(x)
  if (is.null(env$save_idxs)) env$U else env$U[, env$save_idxs, drop = FALSE]
}

#' @export
as.data.frame.SciMLSolution <- function(x, row.names = NULL, optional = FALSE, ...) {
  env <- .sol_env(x)
  # only the saved components (save_idxs), named by their state index when
  # the states have no names
  U <- env$U
  if (is.null(colnames(U))) colnames(U) <- seq_len(ncol(U))
  if (!is.null(env$save_idxs)) U <- U[, env$save_idxs, drop = FALSE]
  df <- data.frame(t = env$t, U, check.names = FALSE)
  sys <- env$sys
  if (!is.null(sys) && length(observed(sys))) {
    for (o in observed(sys)) df[[format(o$lhs)]] <- .sol_series(env, o$lhs)
  }
  df
}

# Values of a variable or expression at the saved times.
.sol_series <- function(env, idx, rows = NULL) {
  U <- env$U
  if (is.null(rows)) rows <- seq_len(nrow(U))
  if (is.numeric(idx) && !is_sym(idx)) return(U[rows, idx])
  sys <- env$sys
  if (is.null(sys)) stop("symbolic indexing needs a problem built from a system", call. = FALSE)
  v <- as_sym(idx)
  ui <- match(.key(v), .keys(unknowns(sys)))
  if (!is.na(ui)) return(unname(U[rows, ui]))
  fn <- .sol_observed_fn(env, v)
  # parameters of the saved rows, so a row saved before an event keeps them
  p <- if (!is.null(env$param_history)) .sol_params_at(env, env$t, saved = TRUE)[rows, , drop = FALSE]
  .eval_observed(fn, env, U[rows, , drop = FALSE], env$t[rows], p = p)
}

.sol_observed_fn <- function(env, v) {
  k <- .key(v)
  fn <- env$obs_cache[[k]]
  if (is.null(fn)) {
    sys <- env$sys
    fn <- .build_closure(list(v), .codegen_args(sys, with_t = !is.null(.f(sys, "iv"))),
                         observed = .codegen_observed(sys), out = "scalar")
    assign(k, fn, envir = env$obs_cache)
  }
  fn
}

.eval_observed <- function(fn, env, U, ts, p = NULL) {
  p <- p %||% .sol_params_at(env, ts)
  has_t <- length(formals(fn)) >= 3L
  vapply(seq_len(nrow(U)), function(i) {
    pp <- if (is.matrix(p)) p[i, ] else p
    if (has_t) fn(U[i, ], pp, ts[i]) else fn(U[i, ], pp)
  }, 0)
}

# Parameters can change at events; the solution records the history when
# the events module provides it, otherwise the final values are used.
.sol_params_at <- function(env, ts, saved = FALSE) {
  h <- env$param_history
  if (is.null(h)) return(env$p_final)
  # the history is in integration order, decreasing for backward solves
  back <- length(h$t) > 1L && h$t[2L] < h$t[1L]
  idx <- if (back) findInterval(-ts, -h$t) else findInterval(ts, h$t, left.open = FALSE)
  idx[idx < 1L] <- 1L
  # an event time is saved twice, before and after the affect: the first
  # saved copy takes the parameters in force before the event
  n <- length(ts)
  if (saved && n > 1L) {
    first <- which(ts[-n] == ts[-1L] & idx[-n] > 1L)
    first <- first[h$t[idx[first]] == ts[first]]
    idx[first] <- idx[first] - 1L
  }
  h$p[idx, , drop = FALSE]
}

#' Index a solution
#'
#' `sol[x]` is the time series of a variable or expression `x`;
#' `sol[list(x, y)]` a matrix with one column each; `sol[x, j]` the value at
#' the `j`th saved time. Numeric indices select state columns.
#'
#' @param x A solution.
#' @param i A variable, expression, list of them, or column index.
#' @param j Time indices.
#' @param ... Unused.
#' @return A numeric vector or matrix.
#' @name solution-indexing
#' @export
`[.SciMLSolution` <- function(x, i, j, ...) {
  env <- .sol_env(x)
  rows <- if (missing(j)) NULL else j
  if (missing(i)) {
    U <- env$U
    return(if (is.null(rows)) U else U[rows, , drop = FALSE])
  }
  if (is.list(i) && !is_sym(i)) {
    cols <- lapply(i, function(z) .sol_series(env, z, rows))
    out <- do.call(cbind, cols)
    colnames(out) <- vapply(i, function(z) if (is_sym(z)) format(z) else as.character(z), "")
    return(out)
  }
  .sol_series(env, i, rows)
}

#' @export
`[[.SciMLSolution` <- function(x, i) .sol_series(.sol_env(x), i)

# Interpolation at arbitrary times.
.sol_interp <- function(env, tq, idxs = NULL, deriv = 0L) {
  deriv <- as.integer(deriv)
  if (length(deriv) != 1L || is.na(deriv) || deriv < 0L) stop("deriv must be 0, 1, 2, ...", call. = FALSE)
  # the copies solve_gradient() differentiates a loss on (R/sensitivity-interface.R)
  if (isTRUE(env$sens_frozen))
    stop("the loss of solve_gradient() may use the saved values (sol$u, sol[x], as.matrix(sol)) but not ",
         "interpolation sol(t); give the times it needs as `saveat`", call. = FALSE)
  U <- .interp_states(env, tq, deriv)
  if (deriv > 0L && !is.null(idxs)) {
    # error_if_observed_derivative of SciMLBase
    ids <- if (is.list(idxs) && !is_sym(idxs)) idxs else list(idxs)
    observed_idx <- vapply(ids, function(z) {
      if (is.numeric(z) && !is_sym(z)) return(FALSE)
      is.na(match(.key(as_sym(z)), .keys(unknowns(env$sys))))
    }, TRUE)
    if (any(observed_idx))
      stop("Cannot interpolate derivatives of observed variables. A possible solution could be\n",
           "interpolating the symbolic expression that evaluates to the derivative of the\n",
           "observed variable or using DataInterpolations.jl.", call. = FALSE)
  }
  if (is.null(idxs)) {
    # save_idxs: the interpolation gives the saved components only
    if (!is.null(env$save_idxs)) U <- U[, env$save_idxs, drop = FALSE]
    if (length(tq) == 1L) {
      v <- U[1L, ]
      # an array state keeps its shape
      if (is.null(env$save_idxs) && !is.null(env$u_dim)) v <- array(v, env$u_dim)
      return(v)
    }
    return(U)
  }
  if (is.list(idxs) && !is_sym(idxs)) {
    out <- vapply(idxs, function(z) .series_from(env, z, U, tq), numeric(length(tq)))
    if (length(tq) == 1L) out <- matrix(out, nrow = 1L)
    colnames(out) <- vapply(idxs, function(z) if (is_sym(z)) format(z) else as.character(z), "")
    return(if (length(tq) == 1L) out[1L, ] else out)
  }
  .series_from(env, idxs, U, tq)
}

.series_from <- function(env, idx, U, tq) {
  if (is.numeric(idx) && !is_sym(idx)) return(U[, idx])
  v <- as_sym(idx)
  ui <- match(.key(v), .keys(unknowns(env$sys)))
  if (!is.na(ui)) return(unname(U[, ui]))
  .eval_observed(.sol_observed_fn(env, v), env, U, tq)
}

.interp_states <- function(env, tq, deriv = 0L) {
  steps <- env$steps
  n <- ncol(env$U)
  out <- matrix(NA_real_, length(tq), n)
  colnames(out) <- colnames(env$U)
  if (is.null(steps) || !length(steps)) {
    # linear interpolation between saved points; a solve that stopped at its
    # first point holds that state (its slope is 0 / 0)
    if (length(env$t) == 1L) {
      if (deriv > 1L) .deriv_order_error()
      out[] <- rep(if (deriv == 1L) NaN else env$U[1L, ], each = length(tq))
      return(out)
    }
    out[] <- .linear_states(env$t, env$U, tq, deriv)
    return(out)
  }
  starts <- vapply(steps, `[[`, 0, "t0")
  ends <- vapply(steps, function(s) s$tend %||% (s$t0 + s$dt), 0)
  fwd <- ends[length(ends)] >= starts[1L]
  for (q in seq_along(tq)) {
    tt <- tq[q]
    # starts increase forward in time and decrease for backward solves
    i <- if (fwd) findInterval(tt, starts, rightmost.closed = TRUE) else
      findInterval(-tt, -starts, rightmost.closed = TRUE)
    i <- max(1L, min(length(steps), i))
    # after an event, prefer the step that starts at the event time
    dti <- ends[i] - starts[i]
    out[q, ] <- if (deriv == 0L) .step_value(steps[[i]], tt, dti) else .step_deriv(steps[[i]], tt, deriv, dti)
  }
  out
}

#' Plot a solution
#'
#' Plots the unknowns (or the variables in `idxs`) against time with base
#' graphics.
#'
#' @param x A solution.
#' @param idxs Variables or expressions to plot against time (a vector of
#'   state indices plots each of them); a list of two plots the first against
#'   the second (a phase plot).
#' @param n Number of interpolation points for dense solutions.
#' @param ... Passed to [graphics::matplot()].
#' @return Invisibly `x`.
#' @export
plot.SciMLSolution <- function(x, idxs = NULL, n = 500L, ...) {
  env <- .sol_env(x)
  ts <- env$t
  if (isTRUE(env$dense) && length(ts) > 1L) ts <- sort(unique(c(seq(ts[1L], ts[length(ts)], length.out = n), env$t)))
  if (is.null(idxs)) {
    # the saved components (save_idxs), labeled by name or by state index
    cols <- env$save_idxs %||% seq_len(ncol(env$U))
    vals <- if (isTRUE(env$dense)) x(ts) else env$U[, cols, drop = FALSE]
    if (!is.matrix(vals)) vals <- matrix(vals, ncol = length(cols))
    labs <- colnames(env$U)[cols] %||% as.character(cols)
  } else {
    idl <- if (is.list(idxs) && !is_sym(idxs)) idxs else if (is.atomic(idxs)) as.list(idxs) else list(idxs)
    vals <- if (isTRUE(env$dense)) x(ts, idxs = idl) else x[idl]
    if (!is.matrix(vals)) vals <- matrix(vals, ncol = length(idl))
    labs <- vapply(idl, function(z) if (is_sym(z)) format(z) else as.character(z), "")
    # idxs = (x, y) in Plots.jl: one variable against the other
    if (is.list(idxs) && !is_sym(idxs) && length(idxs) == 2L) {
      graphics::plot(vals[, 1L], vals[, 2L], type = "l", xlab = labs[1L], ylab = labs[2L], ...)
      return(invisible(x))
    }
  }
  graphics::matplot(ts, vals, type = "l", lty = 1, xlab = "t", ylab = "", ...)
  graphics::legend("topright", legend = labs, col = seq_along(labs), lty = 1, bty = "n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Derivatives of the interpolants, sol(t, deriv = n): the Val{n} methods of
# _ode_interpolant (OrdinaryDiffEqCore dense/generic_dense.jl and the
# interpolants.jl of the solver libraries). The first derivatives are those
# of R/globaldiffeq-companion.R (.gde_dinterp_of); the higher ones follow
# the same polynomials. Each interpolant supports the orders its Julia
# cache does (Hermite 3, Tsit5, DP5 and Vern9 4, the Rodas H matrix 3,
# SSPRK 2, the others 1).

.DERIVATIVE_ORDER_NOT_POSSIBLE_MESSAGE <- paste0(
  "Derivative order too high for interpolation order. An interpolation derivative is\n",
  "only accurate to a certain derivative. For example, a second order interpolation\n",
  "is a quadratic polynomial, and thus third derivatives cannot be computed (will be\n",
  "incorrectly zero). Thus use a solver with a higher order interpolation or compute\n",
  "the higher order derivative through other means.\n\n",
  "You can find the list of available ODE/DAE solvers with their documented interpolations at:\n\n",
  "* https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/\n",
  "* https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/\n")

.deriv_order_error <- function() stop(.DERIVATIVE_ORDER_NOT_POSSIBLE_MESSAGE, call. = FALSE)

# The n-th derivative at theta of the polynomial sum_i c[i] theta^(i - 1).
.dpoly <- function(c, theta, n) {
  p <- seq_along(c) - 1L
  keep <- p >= n & c != 0
  if (!any(keep)) return(0)
  sum(c[keep] * choose(p[keep], n) * factorial(n) * theta^(p[keep] - n))
}

# hermite_interpolant at Val{1}, Val{2} and Val{3}
.hermite_dn <- function(theta, dt, y0, y1, k, n) {
  k1 <- k[[1L]]; k2 <- k[[2L]]
  switch(n,
    .gde_hermite_d(theta, dt, y0, y1, k),
    (-4 * dt * k1 - 2 * dt * k2 - 6 * y0 + theta * (6 * dt * k1 + 6 * dt * k2 + 12 * y0 - 12 * y1) + 6 * y1) /
      (dt * dt),
    (6 * dt * k1 + 6 * dt * k2 + 12 * y0 - 12 * y1) / (dt * dt * dt))
}

# with differential_vars: the algebraic components (dv FALSE, or a missing
# slope when dv is NULL) take the derivative of the linear interpolant,
# (y1 - y0) / dt and then 0
.hermite_dv_dn <- function(theta, dt, y0, y1, k, n, dv) {
  lin <- if (n == 1L) (y1 - y0) / dt else 0 * y0
  if (identical(dv, FALSE)) return(lin)
  h <- .hermite_dn(theta, dt, y0, y1, k, n)
  if (is.null(dv)) return(ifelse(is.na(h), lin, h))
  ifelse(dv, h, lin)
}

# y0 + dt * sum_j K_j b_j(theta) with b_j(theta) = sum_p C[j, p + 1] theta^p:
# the n-th time derivative dt^(1 - n) sum_j K_j b_j^(n)(theta)
.rk_poly_dn <- function(theta, dt, K, C, n) {
  acc <- 0
  for (j in seq_along(K)) {
    w <- .dpoly(C[j, ], theta, n)
    if (w != 0) acc <- acc + K[[j]] * w
  }
  acc * dt^(1 - n)
}

.tsit5_coefs <- function() {
  r <- .TSIT5$r
  C <- matrix(0, 7L, 5L)
  C[1L, 2:5] <- r[[1L]]
  for (i in 2:7) C[i, 3:5] <- r[[i]]
  C
}

# dp5pre: k1 theta + k2 b20 + k3 theta b20 + k4 b20^2 with b20 = theta (1 - theta)
.DP5_COEFS <- rbind(c(0, 1, 0, 0, 0), c(0, 1, -1, 0, 0), c(0, 0, 1, -1, 0), c(0, 0, 1, -2, 1))

# The derivative of order n of an interpolant: list(max, f) with f(theta,
# dt, y0, y1, k, n) and the highest order max (NA when f checks it), or NULL
# when the interpolant has no derivative in the port.
.interp_deriv_spec <- function(interp) {
  if (!is.function(interp)) return(NULL)
  spec <- function(max, f) list(max = max, f = f)
  same_body <- function(g) identical(body(interp), body(g))
  if (identical(interp, .hermite_interp) || identical(interp, .hermite_interp_fused))
    return(spec(3L, .hermite_dn))
  if (identical(interp, .hermite_dae_interp))
    return(spec(3L, function(theta, dt, y0, y1, k, n) .hermite_dv_dn(theta, dt, y0, y1, k, n, NULL)))
  if (identical(interp, .dabdf2_interp))
    return(spec(3L, function(theta, dt, y0, y1, k, n) {
      if (length(k) < 2L) return(if (n == 1L) (y1 - y0) / dt else 0 * y0)
      dv <- k[[3L]]
      if (is.null(dv) || all(dv)) .hermite_dn(theta, dt, y0, y1, k, n) else .hermite_dv_dn(theta, dt, y0, y1, k, n, dv)
    }))
  if (identical(interp, .lsrk_interp) || identical(interp, .lsrk_interp_fused))
    return(spec(3L, function(theta, dt, y0, y1, k, n) .hermite_dn(theta, dt, y0, y1, .lsrk_ends(k, y0, y1), n)))
  if (identical(interp, .newmark_interp))
    return(spec(3L, function(theta, dt, y0, y1, k, n)
      .hermite_dn(theta, dt, y0, y1, list(k$f(y0, k$p, k$t), k$f(y1, k$p, k$t + dt)), n)))
  # DPRKN6 has no Val{n} interpolant: the generic Hermite fallback
  if (same_body(.dprkn6_interp(0, 0)))
    return(spec(3L, function(theta, dt, y0, y1, k, n)
      .hermite_dn(theta, dt, y0, y1, list(c(k[[1L]], k[[2L]]), c(k[[3L]], k[[4L]])), n)))
  if (identical(interp, .tsit5_interp))
    return(spec(4L, function(theta, dt, y0, y1, k, n)
      if (n == 1L) .gde_tsit5_d(theta, dt, y0, y1, k) else .rk_poly_dn(theta, dt, k, .tsit5_coefs(), n)))
  if (identical(interp, .dp5_interp))
    return(spec(4L, function(theta, dt, y0, y1, k, n)
      if (n == 1L) .gde_dp5_d(theta, dt, y0, y1, k) else .rk_poly_dn(theta, dt, k, .DP5_COEFS, n)))
  if (identical(interp, .verner_interp))
    return(spec(NA_integer_, function(theta, dt, y0, y1, k, n) {
      # Vern9 has Val{1} to Val{4}, Vern6, Vern7 and Vern8 only Val{1}
      if (n > (if (identical(k$tab$name, "Vern9")) 4L else 1L)) .deriv_order_error()
      if (n == 1L) return(.gde_verner_d(theta, dt, y0, y1, k))
      if (!k$lazy$done) .verner_interp(theta, dt, y0, y1, k)
      K <- lapply(seq_len(ncol(k$lazy$K)), function(j) k$lazy$K[, j])
      .rk_poly_dn(theta, dt, K, cbind(0, k$tab$R)[seq_along(K), , drop = FALSE], n)
    }))
  if (identical(interp, .ssprk_interp))
    return(spec(2L, function(theta, dt, y0, y1, k, n)
      if (n == 1L) .gde_ssprk_d(theta, dt, y0, y1, k) else (2 * y1 - 2 * y0) / (dt * dt) - 2 * k[[1L]] / dt))
  if (same_body(.rosen_interp(2L))) {
    order <- environment(interp)$order
    if (order < 0L)
      return(spec(3L, function(theta, dt, y0, y1, k, n) {
        m <- attr(k, "dvars")
        if (is.null(m)) .hermite_dn(theta, dt, y0, y1, k, n) else .hermite_dv_dn(theta, dt, y0, y1, k, n, m)
      }))
    # (1 - theta) y0 + theta y1 + theta (1 - theta) sum_j k_j theta^(j - 1)
    return(spec(3L, function(theta, dt, y0, y1, k, n) {
      if (n == 1L) return(.gde_rodas_d(order, theta, dt, y0, y1, k))
      acc <- 0 * y0
      for (j in seq_len(order)) {
        cj <- numeric(j + 2L)
        cj[j + 1L] <- 1
        cj[j + 2L] <- -1
        acc <- acc + k[[j]] * .dpoly(cj, theta, n)
      }
      acc / dt^n
    }))
  }
  if (same_body(.rosen_interp23(0))) {
    d <- environment(interp)$d
    return(spec(1L, function(theta, dt, y0, y1, k, n) .gde_rosenbrock23_d(d)(theta, dt, y0, y1, k)))
  }
  if (identical(interp, .qndf_interp))
    return(spec(1L, function(theta, dt, y0, y1, k, n) {
      # the QNDF Val{1}: phi'_(j+1) = (phi'_j (s + j) + phi_j) / (j + 1)
      s <- theta - 1
      phi <- s
      dphi <- 1
      out <- dphi * k[[1L]]
      if (length(k) >= 2L) for (j in 2:length(k)) {
        dphi_new <- (dphi * (s + j - 1) + phi) / j
        phi <- phi * ((s + j - 1) / j)
        dphi <- dphi_new
        out <- out + dphi * k[[j]]
      }
      out * (1 / dt)
    }))
  if (identical(interp, .bdf_interp))
    return(spec(1L, function(theta, dt, y0, y1, k, n) .gde_bdf_d(theta, dt, y0, y1, k)))
  if (identical(interp, .taylor2_interp))
    return(spec(1L, function(theta, dt, y0, y1, k, n) k[[1L]] + k[[2L]] * (theta * dt)))
  if (identical(interp, .taylor_interp))
    return(spec(1L, function(theta, dt, y0, y1, k, n) {
      P <- length(k)
      h <- theta * dt
      res <- P * k[[P]]
      if (P > 1L) for (i in (P - 1L):1L) res <- i * k[[i]] + h * res
      res
    }))
  if (identical(interp, .sde_linear_interp))
    return(spec(1L, function(theta, dt, y0, y1, k, n) (y1 - y0) / dt))
  # RadauIIA3, 5, 9 and AdaptiveRadau have only Val{0}: the generic
  # fallback of _ode_interpolant, the Hermite interpolant of k[1] = f(y0)
  # and k[2] = f(y1)
  if (identical(interp, .firk_radau_interp))
    return(spec(3L, function(theta, dt, y0, y1, k, n) .hermite_dn(theta, dt, y0, y1, k[1:2], n)))
  if (identical(interp, .nordbdf_interp))
    return(spec(3L, .nordbdf_dinterp))
  # the MIRK interpolant has Val{0} and Val{1} (BoundaryValueDiffEqMIRK
  # interpolation.jl): sum_stages with the derivative weights w'
  if (identical(interp, .mirk_interp))
    return(spec(1L, function(theta, dt, y0, y1, k, n) .mirk_eval(k$name, k$K, y0, dt, theta, deriv = TRUE)))
  if (identical(interp, .dopt_interp))
    return(spec(1L, .dopt_dinterp))
  NULL
}

# d^m/dt^m of the Nordsieck polynomial sum_j zn[j] s^j (bdf_interpolants.jl
# of OrdinaryDiffEqBDF): Horner over j = n, ..., m + 1 with the falling
# factorial (j - 1) ... (j - m), divided by dt^m. Without columns (before
# the first step) the interpolant is linear.
.nordbdf_dinterp <- function(theta, dt, y0, y1, k, n) {
  kk <- if (is.environment(k)) k$k else NULL
  if (is.null(kk)) return(if (n == 1L) (y1 - y0) / dt else 0 * y0)
  s <- theta - 1
  nk <- length(kk)
  out <- 0 * kk[[1L]]
  if (nk > n) for (j in nk:(n + 1L)) {
    cj <- 1
    for (r in 0:(n - 1L)) cj <- cj * (j - 1L - r)
    out <- .fma(out, s, cj * kk[[j]])
  }
  out / switch(n, dt, dt * dt, dt * dt * dt)
}

# The time derivative of the collocation interpolant of a dynamic
# optimization solution (.dopt_interp): of the Lagrange polynomial through
# the stage values, of x_k + h sum_l beta_l(theta) f_l for the differential
# states, or of the cubic Hermite form.
.dopt_dinterp <- function(theta, dt, y0, y1, k, n) {
  th <- k$th
  if (length(unique(round(th, 12))) == length(th)) {
    dl <- vapply(seq_along(th), function(i) {
      others <- th[-i]
      sum(vapply(seq_along(others), function(m)
        prod((theta - others[-m]) / (th[i] - others[-m])) / (th[i] - others[m]), 0))
    }, 0)
    out <- .as_num(crossprod(k$V, dl)) / dt
  } else {
    out <- (y1 - y0) / dt
  }
  if (length(k$diff)) {
    j <- k$diff
    if (!is.null(k$K)) {
      q <- ncol(k$beta)
      bt <- .as_num(k$beta %*% (seq_len(q) * theta^(seq_len(q) - 1L)))
      out[j] <- .as_num(crossprod(k$K, bt))
    } else {
      d00 <- 6 * theta^2 - 6 * theta; d10 <- 3 * theta^2 - 4 * theta + 1
      d01 <- -6 * theta^2 + 6 * theta; d11 <- 3 * theta^2 - 2 * theta
      out[j] <- (d00 * y0[j] + d01 * y1[j]) / dt + d10 * k$f0 + d11 * k$f1
    }
  }
  names(out) <- names(y0)
  out
}

# The generic _ode_interpolant of OrdinaryDiffEqCore for a cache without
# its own Val{n} method: the linear interpolant (Val{1} only) without two
# dense output vectors, else the Hermite interpolant of k[1] and k[2] (orders up to 3).
# NULL when the step's k is not a list of state sized vectors.
.generic_deriv_spec <- function(step) {
  k <- step$k
  n <- length(step$y0)
  if (!length(k) || (is.list(k) && length(k) < 2L))
    return(list(max = 1L, f = function(theta, dt, y0, y1, k, n) (y1 - y0) / dt))
  ok <- is.list(k) && !is.environment(k) && all(vapply(k[1:2], function(v) (is.numeric(v) || is.complex(v)) &&
                                                         length(v) == n, TRUE))
  if (!ok) return(NULL)
  list(max = 3L, f = function(theta, dt, y0, y1, k, n) .hermite_dn(theta, dt, y0, y1, k[1:2], n))
}

# The n-th derivative of the interpolant of a dense step at time t.
.step_deriv <- function(step, t, n, dt = step$dt) {
  # dt as in .step_value: sol(t) passes the difference of the saved times
  sp <- .interp_deriv_spec(step$interp) %||% .generic_deriv_spec(step)
  if (is.null(sp)) stop("the interpolant of this method has no derivatives in the port", call. = FALSE)
  if (!is.na(sp$max) && n > sp$max) .deriv_order_error()
  theta <- if (dt == 0) 1 else (t - step$t0) / dt
  v <- sp$f(theta, dt, step$y0, step$y1, step$k, n)
  if (is.complex(v)) as.vector(v) else .as_num(v)
}

# Linear interpolation between the saved points (linear_interpolant), the
# value or (deriv = 1) the slope, with OrdinaryDiffEq's left continuity
# ts[i-] < t <= ts[i+] in the direction of integration, so a backward
# solution and a time saved twice at an event are handled. Outside the
# saved times the end intervals extrapolate (ode_interpolation does not
# clamp theta).
.linear_states <- function(ts, U, tq, deriv) {
  n <- length(ts)
  if (deriv > 1L) .deriv_order_error()
  tdir <- if (ts[n] < ts[1L]) -1 else 1
  i <- findInterval(tdir * tq, tdir * ts, left.open = TRUE)
  lo <- pmin(pmax(i, 1L), n - 1L)
  hi <- lo + 1L
  dt <- ts[hi] - ts[lo]
  if (deriv == 1L) return((U[hi, , drop = FALSE] - U[lo, , drop = FALSE]) / dt)
  theta <- ifelse(dt == 0, 1, (tq - ts[lo]) / dt)
  .fma(theta, U[hi, , drop = FALSE], (1 - theta) * U[lo, , drop = FALSE])
}
