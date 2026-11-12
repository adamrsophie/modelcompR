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
      lapply(seq_len(nrow(U)), function(i) U[i, ])
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
  df <- data.frame(t = env$t, env$U, check.names = FALSE)
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
  if (deriv != 0L) stop("derivatives of the interpolant are not supported", call. = FALSE)
  # the copies solve_gradient() differentiates a loss on (R/sensitivity-interface.R)
  if (isTRUE(env$sens_frozen))
    stop("the loss of solve_gradient() may use the saved values (sol$u, sol[x], as.matrix(sol)) but not ",
         "interpolation sol(t); give the times it needs as `saveat`", call. = FALSE)
  U <- .interp_states(env, tq)
  if (is.null(idxs)) {
    if (length(tq) == 1L) return(U[1L, ])
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

.interp_states <- function(env, tq) {
  steps <- env$steps
  n <- ncol(env$U)
  out <- matrix(NA_real_, length(tq), n)
  colnames(out) <- colnames(env$U)
  if (is.null(steps) || !length(steps)) {
    # linear interpolation between saved points
    for (j in seq_len(n)) out[, j] <- stats::approx(env$t, env$U[, j], xout = tq, rule = 2, ties = "ordered")$y
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
    out[q, ] <- .step_value(steps[[i]], tt)
  }
  out
}

#' Plot a solution
#'
#' Plots the unknowns (or the variables in `idxs`) against time with base
#' graphics.
#'
#' @param x A solution.
#' @param idxs Variables or expressions to plot; a list of two plots one
#'   against the other.
#' @param n Number of interpolation points for dense solutions.
#' @param ... Passed to [graphics::matplot()].
#' @return Invisibly `x`.
#' @export
plot.SciMLSolution <- function(x, idxs = NULL, n = 500L, ...) {
  env <- .sol_env(x)
  ts <- env$t
  if (isTRUE(env$dense) && length(ts) > 1L) ts <- sort(unique(c(seq(ts[1L], ts[length(ts)], length.out = n), env$t)))
  if (is.null(idxs)) {
    vals <- if (isTRUE(env$dense)) x(ts) else env$U
    if (!is.matrix(vals)) vals <- matrix(vals, ncol = ncol(env$U))
    labs <- colnames(env$U)
  } else {
    idl <- if (is.list(idxs) && !is_sym(idxs)) idxs else list(idxs)
    vals <- if (isTRUE(env$dense)) x(ts, idxs = idl) else x[idl]
    if (!is.matrix(vals)) vals <- matrix(vals, ncol = length(idl))
    labs <- vapply(idl, function(z) if (is_sym(z)) format(z) else as.character(z), "")
  }
  graphics::matplot(ts, vals, type = "l", lty = 1, xlab = "t", ylab = "", ...)
  graphics::legend("topright", legend = labs, col = seq_along(labs), lty = 1, bty = "n")
  invisible(x)
}
