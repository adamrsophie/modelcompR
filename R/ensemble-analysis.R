# The EnsembleAnalysis module of SciMLBase.jl
# (src/ensemble/ensemble_analysis.jl): statistics of ensemble solutions at a
# saved step index, at a time point, along all saved steps and along given
# times, and the componentwise statistics they are built on.
#
# Observations are numeric vectors (one per trajectory); a "DiffEqArray"
# result is a list with the times `t` and a matrix `u` with one row per time.

#' Ensemble analysis
#'
#' Statistics over the trajectories of an [EnsembleProblem()] solution, as
#' SciMLBase's `EnsembleAnalysis`. `get_timestep(sim, i)` is the list of the
#' `i`th saved states and `get_timepoint(sim, t)` of the states at time `t`;
#' `componentwise_vectors_timestep()` and `componentwise_vectors_timepoint()`
#' regroup them by component.
#'
#' * At a saved step `i`: `timestep_mean`, `timestep_median`,
#'   `timestep_quantile(sim, q, i)`, `timestep_meanvar`, and between steps
#'   `i` and `j` `timestep_meancov`, `timestep_meancor`,
#'   `timestep_weighted_meancov(sim, W, i, j)` (`W` one weight per
#'   trajectory).
#' * At a time `t` (solutions are interpolated): the `timepoint_` versions.
#' * Along all saved steps (the steps of the first trajectory):
#'   `timeseries_steps_mean`, `_median`, `_quantile`, `_meanvar`,
#'   `_meancov`, `_meancor`, `_weighted_meancov`.
#' * Along times `ts`: the `timeseries_point_` versions.
#' * On lists of observations: `componentwise_mean(A)`,
#'   `componentwise_meanvar(A, bessel)`, `componentwise_meancov(A, B,
#'   bessel)`, `componentwise_meancor(A, B, bessel)` and
#'   `componentwise_weighted_meancov(A, B, W, weight_type)` with
#'   `weight_type` one of `"reliability"`, `"frequency"`, `"population"`.
#'
#' Mean and variance are returned as `list(mean, var)`, covariances as
#' `list(meanx, meany, cov)`. The series functions return `list(t, u)` with
#' one row of `u` per time (two of them for mean and variance); the
#' covariance series return matrices of such lists, entry `(i, j)` for steps
#' (or times) `i` and `j`.
#'
#' @param sim An `EnsembleSolution` of solutions.
#' @param i,j Saved step indices.
#' @param t,t1,t2 Times.
#' @param ts,ts1,ts2 Vectors of times.
#' @param q A quantile level.
#' @param W Weights, one per trajectory.
#' @param A,B Lists of observations (numeric vectors of the same length).
#' @param bessel Use the sample (`n - 1`) denominator.
#' @param weight_type How weights enter the covariance.
#' @return See Details.
#' @examples
#' prob <- SDEProblem(function(u, p, t) u, function(u, p, t) 0.2 * u, 1, c(0, 1))
#' sim <- solve(EnsembleProblem(prob), EM(), trajectories = 20, dt = 0.1, seed = 1)
#' timepoint_meanvar(sim, 0.5)
#' componentwise_meanvar(list(c(1, 2), c(3, 4)))
#' @name ensemble_analysis
NULL

.ea_sols <- function(sim) .subset2(sim, "u")
.ea_rows <- function(A) do.call(rbind, lapply(A, as.numeric))
.ea_shape <- function(x, like) if (length(like) == 1L) x[[1L]] else x

#' @rdname ensemble_analysis
#' @export
get_timestep <- function(sim, i) lapply(.ea_sols(sim), function(s) s$u[[i]])

#' @rdname ensemble_analysis
#' @export
get_timepoint <- function(sim, t) lapply(.ea_sols(sim), function(s) as.numeric(s(t)))

#' @rdname ensemble_analysis
#' @export
componentwise_vectors_timestep <- function(sim, i) {
  X <- .ea_rows(get_timestep(sim, i))
  if (ncol(X) == 1L) X[, 1L] else lapply(seq_len(ncol(X)), function(k) X[, k])
}

#' @rdname ensemble_analysis
#' @export
componentwise_vectors_timepoint <- function(sim, t) {
  X <- .ea_rows(get_timepoint(sim, t))
  if (ncol(X) == 1L) X[, 1L] else lapply(seq_len(ncol(X)), function(k) X[, k])
}

#' @rdname ensemble_analysis
#' @export
componentwise_mean <- function(A) colMeans(.ea_rows(A))

#' @rdname ensemble_analysis
#' @export
componentwise_meanvar <- function(A, bessel = TRUE) {
  X <- .ea_rows(A)
  n <- nrow(X)
  if (n < 2L) return(NaN)
  m <- colMeans(X)
  v <- colSums(sweep(X, 2L, m)^2) / (if (bessel) n - 1 else n)
  list(mean = m, var = v)
}

#' @rdname ensemble_analysis
#' @export
componentwise_meancov <- function(A, B, bessel = TRUE) {
  X <- .ea_rows(A); Y <- .ea_rows(B)
  n <- nrow(X)
  if (n < 2L) return(NaN)
  mx <- colMeans(X); my <- colMeans(Y)
  C <- colSums(sweep(X, 2L, mx) * sweep(Y, 2L, my)) / (if (bessel) n - 1 else n)
  list(meanx = mx, meany = my, cov = C)
}

#' @rdname ensemble_analysis
#' @export
componentwise_meancor <- function(A, B, bessel = TRUE) {
  mc <- componentwise_meancov(A, B, bessel)
  if (!is.list(mc)) return(mc)
  vx <- componentwise_meanvar(A, bessel)$var
  vy <- componentwise_meanvar(B, bessel)$var
  list(meanx = mc$meanx, meany = mc$meany, cor = mc$cov / sqrt(vx * vy))
}

#' @rdname ensemble_analysis
#' @export
componentwise_weighted_meancov <- function(A, B, W, weight_type = c("reliability", "frequency", "population")) {
  weight_type <- match.arg(weight_type)
  X <- .ea_rows(A); Y <- .ea_rows(B)
  if (nrow(X) < 2L) return(NaN)
  w <- as.numeric(W)
  ws <- sum(w); ws2 <- sum(w^2)
  mx <- colSums(w * X) / ws; my <- colSums(w * Y) / ws
  C <- colSums(w * sweep(X, 2L, mx) * sweep(Y, 2L, my))
  C <- C / switch(weight_type, population = ws, reliability = ws - ws2 / ws, frequency = ws - 1)
  list(meanx = mx, meany = my, cov = C)
}

.ea_median <- function(A) apply(.ea_rows(A), 2L, stats::median)
.ea_quantile <- function(A, q) apply(.ea_rows(A), 2L, stats::quantile, probs = q, names = FALSE, type = 7)

#' @rdname ensemble_analysis
#' @export
timestep_mean <- function(sim, i) if (missing(i) || is.null(i)) timeseries_steps_mean(sim) else componentwise_mean(get_timestep(sim, i))

#' @rdname ensemble_analysis
#' @export
timestep_median <- function(sim, i) .ea_median(get_timestep(sim, i))

#' @rdname ensemble_analysis
#' @export
timestep_quantile <- function(sim, q, i) .ea_quantile(get_timestep(sim, i), q)

#' @rdname ensemble_analysis
#' @export
timestep_meanvar <- function(sim, i) componentwise_meanvar(get_timestep(sim, i))

#' @rdname ensemble_analysis
#' @export
timestep_meancov <- function(sim, i, j) componentwise_meancov(get_timestep(sim, i), get_timestep(sim, j))

#' @rdname ensemble_analysis
#' @export
timestep_meancor <- function(sim, i, j) componentwise_meancor(get_timestep(sim, i), get_timestep(sim, j))

#' @rdname ensemble_analysis
#' @export
timestep_weighted_meancov <- function(sim, W, i, j)
  componentwise_weighted_meancov(get_timestep(sim, i), get_timestep(sim, j), W)

#' @rdname ensemble_analysis
#' @export
timepoint_mean <- function(sim, t) componentwise_mean(get_timepoint(sim, t))

#' @rdname ensemble_analysis
#' @export
timepoint_median <- function(sim, t) .ea_median(get_timepoint(sim, t))

#' @rdname ensemble_analysis
#' @export
timepoint_quantile <- function(sim, q, t) .ea_quantile(get_timepoint(sim, t), q)

#' @rdname ensemble_analysis
#' @export
timepoint_meanvar <- function(sim, t) componentwise_meanvar(get_timepoint(sim, t))

#' @rdname ensemble_analysis
#' @export
timepoint_meancov <- function(sim, t1, t2) componentwise_meancov(get_timepoint(sim, t1), get_timepoint(sim, t2))

#' @rdname ensemble_analysis
#' @export
timepoint_meancor <- function(sim, t1, t2) componentwise_meancor(get_timepoint(sim, t1), get_timepoint(sim, t2))

#' @rdname ensemble_analysis
#' @export
timepoint_weighted_meancov <- function(sim, W, t1, t2)
  componentwise_weighted_meancov(get_timepoint(sim, t1), get_timepoint(sim, t2), W)

.ea_steps <- function(sim) seq_along(.ea_sols(sim)[[1L]]$t)
.ea_series <- function(t, rows) list(t = t, u = do.call(rbind, rows))
.ea_series2 <- function(t, pairs) list(.ea_series(t, lapply(pairs, `[[`, 1L)), .ea_series(t, lapply(pairs, `[[`, 2L)))
.ea_grid <- function(n1, n2, f) {
  out <- matrix(list(), n1, n2)
  for (a in seq_len(n1)) for (b in seq_len(n2)) out[[a, b]] <- f(a, b)
  out
}

#' @rdname ensemble_analysis
#' @export
timeseries_steps_mean <- function(sim) .ea_series(.ea_sols(sim)[[1L]]$t, lapply(.ea_steps(sim), function(i) timestep_mean(sim, i)))

#' @rdname ensemble_analysis
#' @export
timeseries_steps_median <- function(sim) .ea_series(.ea_sols(sim)[[1L]]$t, lapply(.ea_steps(sim), function(i) timestep_median(sim, i)))

#' @rdname ensemble_analysis
#' @export
timeseries_steps_quantile <- function(sim, q)
  .ea_series(.ea_sols(sim)[[1L]]$t, lapply(.ea_steps(sim), function(i) timestep_quantile(sim, q, i)))

#' @rdname ensemble_analysis
#' @export
timeseries_steps_meanvar <- function(sim)
  .ea_series2(.ea_sols(sim)[[1L]]$t, lapply(.ea_steps(sim), function(i) timestep_meanvar(sim, i)))

#' @rdname ensemble_analysis
#' @export
timeseries_steps_meancov <- function(sim) {
  n <- length(.ea_steps(sim))
  .ea_grid(n, n, function(i, j) timestep_meancov(sim, i, j))
}

#' @rdname ensemble_analysis
#' @export
timeseries_steps_meancor <- function(sim) {
  n <- length(.ea_steps(sim))
  .ea_grid(n, n, function(i, j) timestep_meancor(sim, i, j))
}

#' @rdname ensemble_analysis
#' @export
timeseries_steps_weighted_meancov <- function(sim, W) {
  n <- length(.ea_steps(sim))
  .ea_grid(n, n, function(i, j) timestep_weighted_meancov(sim, W, i, j))
}

#' @rdname ensemble_analysis
#' @export
timeseries_point_mean <- function(sim, ts) .ea_series(ts, lapply(ts, function(t) timepoint_mean(sim, t)))

#' @rdname ensemble_analysis
#' @export
timeseries_point_median <- function(sim, ts) .ea_series(ts, lapply(ts, function(t) timepoint_median(sim, t)))

#' @rdname ensemble_analysis
#' @export
timeseries_point_quantile <- function(sim, q, ts) .ea_series(ts, lapply(ts, function(t) timepoint_quantile(sim, q, t)))

#' @rdname ensemble_analysis
#' @export
timeseries_point_meanvar <- function(sim, ts) .ea_series2(ts, lapply(ts, function(t) timepoint_meanvar(sim, t)))

#' @rdname ensemble_analysis
#' @export
timeseries_point_meancov <- function(sim, ts1, ts2 = ts1)
  .ea_grid(length(ts1), length(ts2), function(a, b) timepoint_meancov(sim, ts1[a], ts2[b]))

#' @rdname ensemble_analysis
#' @export
timeseries_point_meancor <- function(sim, ts1, ts2 = ts1)
  .ea_grid(length(ts1), length(ts2), function(a, b) timepoint_meancor(sim, ts1[a], ts2[b]))

#' @rdname ensemble_analysis
#' @export
timeseries_point_weighted_meancov <- function(sim, W, ts1, ts2 = ts1)
  .ea_grid(length(ts1), length(ts2), function(a, b) timepoint_weighted_meancov(sim, W, ts1[a], ts2[b]))
