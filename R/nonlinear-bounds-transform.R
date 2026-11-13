# The bounds transform of NonlinearSolveBase (src/bounds_transform.jl): an
# algorithm without native bounds (allowsbounds false) solves a bounded
# NonlinearProblem or NonlinearLeastSquaresProblem in unbounded coordinates
# t, with u = lb + (ub - lb) logistic(t) for a two sided bound, u = lb +
# exp(t) or u = ub - exp(t) for a one sided bound, and u = t otherwise. The
# start is nudged into the interior first, the Jacobian is scaled by the
# derivative of the map, and the solution is mapped back.

.nlbt_logistic <- function(x) {
  # LogExpFunctions.logistic for Float64
  e <- exp(x)
  ifelse(x < -744.4400719213812, 0, ifelse(x > 36.7368005696771, 1, e / (1 + e)))
}

.nlbt_to_unbounded <- function(u, lb, ub) {
  hl <- is.finite(lb); hu <- is.finite(ub)
  out <- u
  b <- hl & hu
  if (any(b)) { q <- (u[b] - lb[b]) / (ub[b] - lb[b]); out[b] <- log(q / (1 - q)) }
  l <- hl & !hu
  if (any(l)) out[l] <- log(u[l] - lb[l])
  r <- !hl & hu
  if (any(r)) out[r] <- log(ub[r] - u[r])
  out
}

.nlbt_from_unbounded <- function(t, lb, ub) {
  hl <- is.finite(lb); hu <- is.finite(ub)
  out <- t
  b <- hl & hu
  if (any(b)) out[b] <- lb[b] + (ub[b] - lb[b]) * .nlbt_logistic(t[b])
  l <- hl & !hu
  if (any(l)) out[l] <- lb[l] + exp(t[l])
  r <- !hl & hu
  if (any(r)) out[r] <- ub[r] - exp(t[r])
  out
}

.nlbt_derivative <- function(t, lb, ub) {
  hl <- is.finite(lb); hu <- is.finite(ub)
  out <- rep(1, length(t))
  b <- hl & hu
  if (any(b)) out[b] <- (ub[b] - lb[b]) * .nlbt_logistic(t[b]) * .nlbt_logistic(-t[b])
  l <- hl & !hu
  if (any(l)) out[l] <- exp(t[l])
  r <- !hl & hu
  if (any(r)) out[r] <- -exp(t[r])
  out
}

# _clamp_to_bounds: into the strict interior, by eps^(3/4) of the width (or
# of max(|bound|, 1) for a one sided bound)
.nlbt_clamp <- function(u, lb, ub) {
  ef <- .Machine$double.eps^(3 / 4)
  hl <- is.finite(lb); hu <- is.finite(ub)
  out <- u
  b <- hl & hu
  if (any(b)) {
    margin <- (ub[b] - lb[b]) * ef
    out[b] <- pmin(pmax(u[b], lb[b] + margin), ub[b] - margin)
  }
  l <- hl & !hu
  if (any(l)) out[l] <- pmax(u[l], lb[l] + ef * pmax(abs(lb[l]), 1))
  r <- !hl & hu
  if (any(r)) out[r] <- pmin(u[r], ub[r] - ef * pmax(abs(ub[r]), 1))
  out
}

# needs_bounds_transform: bounds, an explicit algorithm, and no native
# support for bounds
.nlbt_needed <- function(alg, lb, ub) (!is.null(lb) || !is.null(ub)) && !isTRUE(alg$allows_bounds)

# Run `alg` on the transformed problem and map the result back.
.nlbt_run <- function(alg, F, J, u0, opts, lb, ub, ls) {
  n <- length(u0)
  lb <- if (is.null(lb)) rep(-Inf, n) else rep_len(as.numeric(lb), n)
  ub <- if (is.null(ub)) rep(Inf, n) else rep_len(as.numeric(ub), n)
  from <- function(t) .nlbt_from_unbounded(t, lb, ub)
  Ft <- function(t) F(from(t))
  Jt <- if (!is.null(J)) function(t) {
    Jx <- J(from(t))
    Jx <- if (is.matrix(Jx)) Jx else matrix(as.numeric(Jx), ncol = n)
    sweep(Jx, 2L, .nlbt_derivative(t, lb, ub), `*`)
  }
  if (is.function(opts$postcondition)) {
    # a corrector written in the original coordinates
    post <- opts$postcondition
    opts$postcondition <- function(t, tprev)
      .nlbt_to_unbounded(.nlbt_clamp(post(from(t), from(tprev)), lb, ub), lb, ub)
  }
  t0 <- .nlbt_to_unbounded(.nlbt_clamp(as.numeric(u0), lb, ub), lb, ub)
  res <- .nl_run(alg, Ft, Jt, t0, opts, NULL, NULL, ls)
  res$u <- from(res$u)
  res
}
