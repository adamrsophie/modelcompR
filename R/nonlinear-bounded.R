# The bound constrained solvers of NonlinearSolveFirstOrder (v2.10.0):
# native_bounded.jl (the shared NativeBoundedAlgorithm driver: strictly
# feasible starts, the trust radius update, the projected Armijo line search
# through LineSearch.ProjectedBackTracking, the scaled and augmented least
# squares solve and its two dimensional trust region subproblem),
# box_constraints.jl and bounded_jacobian.jl (projected gradients, active
# sets and the bound aware finite difference Jacobian),
# trust_region_reflective.jl (TrustRegionReflective, Coleman and Li),
# bounded_levenberg_marquardt.jl (BoundedLevenbergMarquardt, Kanzow,
# Yamashita and Fukushima), dogbox.jl (Dogbox, Voglis and Lagaris),
# bounded_gauss_newton.jl (BoundedGaussNewton) and the BoundedTrustRegion
# of trust_region.jl (BoundedTrustRegionScheme driven by the
# GeneralizedFirstOrderAlgorithm loop of solve.jl with the Dogleg descent of
# NonlinearSolveBase/src/descent/dogleg.jl).
#
# The kernels count the statistics of SciMLBase.NLStats (nf, njacs,
# nfactors, nsolve, nsteps) as the Julia caches do, so they can be compared
# with Julia run for run.

#' Bound constrained nonlinear solvers
#'
#' Solvers for [NonlinearProblem()]s and [NonlinearLeastSquaresProblem()]s
#' with bounds `lb <= u <= ub`, as in NonlinearSolveFirstOrder. The initial
#' guess must satisfy the bounds.
#'
#' * `TrustRegionReflective()` is the interior reflective trust region
#'   method of Coleman and Li: the iterates stay strictly inside the box (a
#'   start on a bound is moved inside), the model uses distance to bound
#'   scaling, and each step is the best of a scaled trust region step, its
#'   reflection at the first bound it hits and a scaled gradient step.
#' * `BoundedLevenbergMarquardt(damping, max_backtracks)` solves the
#'   regularized Gauss Newton model over the box with an active set method
#'   and globalizes with a projected Armijo search, falling back to the
#'   projected gradient.
#' * `Dogbox()` follows the dogleg path from a constrained Cauchy point to
#'   the Gauss Newton step inside the intersection of the box and an infinity
#'   norm trust region.
#' * `BoundedGaussNewton(globalization)` is the active set Gauss Newton
#'   method with a projected line search (`"linesearch"`), a spherical trust
#'   region with a projected Cauchy safeguard (`"trustregion"`) or both
#'   (`"trustregion_linesearch"`).
#' * `BoundedTrustRegion()` compares a projected dogleg step with a feasible
#'   projected Cauchy step and keeps the one with the larger predicted
#'   reduction.
#'
#' A stationary point that is not a root ends a `NonlinearProblem` with
#' `"Stalled"`; for least squares problems the projected gradient test with
#' `gtol` (default `sqrt(eps)`) ends the solve with `"Success"`. The default
#' termination of the first four methods is `AbsNormTerminationMode` with
#' the 2-norm. With `autodiff = AutoFiniteDiff()` and no Jacobian the
#' Jacobian is the bound aware finite difference of NonlinearSolve, which
#' never steps outside the box; otherwise the problem's Jacobian or central
#' differences are used. `sol$stats` holds `nf`, `njacs`, `nfactors`,
#' `nsolve` and `nsteps`.
#'
#' @param autodiff An ADType (see [AutoFiniteDiff()]), or `NULL`.
#' @param jvp_autodiff,vjp_autodiff,linsolve,concrete_jac Options of the
#'   Julia signatures; the linear least squares problems are solved by QR.
#' @param gtol Tolerance of the projected gradient test.
#' @param jacobian_reuse The Jacobian reuse policy of `BoundedTrustRegion()`:
#'   `NULL` (reuse for problems with at least 16 unknowns), a logical, or a
#'   [JacobianReuse()] object.
#' @param initial_trust_radius,max_trust_radius Trust radius settings.
#' @param damping Initial relative damping of `BoundedLevenbergMarquardt()`.
#' @param max_backtracks Largest number of trial points of each line search.
#' @param globalization `"linesearch"`, `"trustregion"` or
#'   `"trustregion_linesearch"`.
#' @param step_threshold,shrink_threshold,expand_threshold,shrink_factor,expand_factor,max_shrink_times
#'   Settings of the trust region of `BoundedTrustRegion()` (a radius of 0
#'   selects the default).
#' @return An algorithm object.
#' @examples
#' f <- function(u, p) c(10 * (u[2] - u[1] * u[1]), 1 - u[1])
#' prob <- NonlinearProblem(f, u0 = c(-1.2, 1), lb = c(-2, -2), ub = c(2, 2))
#' solve(prob, TrustRegionReflective())$u
#' solve(prob, BoundedGaussNewton())$u
#' sol <- solve(prob, Dogbox())
#' sol$stats$nf
#' @name bounded_nonlinear
NULL

# ---------------------------------------------------------------------------
# Floating point helpers (Julia eps(x), nextfloat, prevfloat, isapprox, clamp)

.nb_tiny <- 4.9406564584124654e-324

# eps(x): the distance from |x| to the next larger double
.nb_ulp <- function(x) {
  x <- abs(x)
  if (is.na(x) || is.infinite(x)) return(NaN)
  if (x < 0x1.0000000000000p-1022) return(.nb_tiny)
  e <- floor(log2(x))
  if (2^e > x) e <- e - 1 else if (2^(e + 1) <= x) e <- e + 1
  2^(e - 52)
}

# the largest double below a positive finite x
.nb_below <- function(x) {
  if (x <= 0x1.0000000000000p-1022) return(x - .nb_tiny)
  u <- .nb_ulp(x)
  e <- log2(x)
  if (e == floor(e)) x - u / 2 else x - u
}

.nb_nextfloat <- function(x) {
  if (is.infinite(x)) return(if (x > 0) x else -.Machine$double.xmax)
  if (x == 0) return(.nb_tiny)
  if (x > 0) x + .nb_ulp(x) else -.nb_below(-x)
}

.nb_prevfloat <- function(x) -.nb_nextfloat(-x)

.nb_isapprox <- function(x, y, rtol) x == y || (is.finite(x) && is.finite(y) && abs(x - y) <= rtol * max(abs(x), abs(y)))

# Julia's clamp, elementwise: hi wins when the interval is empty
.nb_clamp <- function(x, lo, hi) ifelse(x > hi, hi, ifelse(x < lo, lo, x))

.nb_norm <- function(v) sqrt(sum(v * v))

.nb_default_tol <- 3.0e-13

# ---------------------------------------------------------------------------
# box_constraints.jl

.nb_projected_gradient <- function(x, g, lb, ub) x - .nb_clamp(x - g, lb, ub)

.nb_active <- function(x, g, lb, ub) lb == ub | (x <= lb & g > 0) | (x >= ub & g < 0)

.nb_model_value <- function(Jstep, step, g) sum(step * g) + sum(Jstep * Jstep) / 2

.nb_model <- function(Jm, g, s) .nb_model_value(as.numeric(Jm %*% s), s, g)

.nb_cauchy_length <- function(slope, curvature, limit)
  if (curvature > 0) .nb_clamp(-slope / curvature, 0, limit) else limit

.nb_update_radius <- function(radius, ratio, stepnorm, max_radius, shrink_threshold = 0.25,
                              expand_threshold = 0.75, shrink_radius = radius / 4, expand_factor = 2) {
  shrunk <- ratio < shrink_threshold
  if (shrunk) radius <- shrink_radius
  else if (ratio > expand_threshold && stepnorm >= 0.95 * radius) radius <- radius * expand_factor
  list(radius = min(radius, max_radius), shrunk = shrunk)
}

.nb_limit <- function(x, d, lb, ub) {
  a <- Inf
  for (i in seq_along(x)) {
    if (d[i] > 0) a <- min(a, (ub[i] - x[i]) / d[i])
    else if (d[i] < 0) a <- min(a, (lb[i] - x[i]) / d[i])
  }
  max(0, a)
}

.nb_interior <- function(x, lb, ub, initial = FALSE) {
  vapply(seq_along(x), function(i) {
    xi <- x[i]; lo <- lb[i]; hi <- ub[i]
    if (lo == hi) return(lo)
    left <- if (is.finite(lo)) .nb_nextfloat(lo) else lo
    right <- if (is.finite(hi)) .nb_prevfloat(hi) else hi
    if (!(left <= right)) stop("A nonfixed interval has no floating-point interior.", call. = FALSE)
    if (initial) {
      margin <- sqrt(.Machine$double.eps) * max(1, abs(xi))
      margin <- min(margin, hi / 2 - lo / 2)
      xi <- .nb_clamp(xi, max(left, lo + margin), min(right, hi - margin))
    }
    .nb_clamp(xi, left, right)
  }, 0)
}

.nb_bounds <- function(b, fill, n) {
  if (is.null(b)) return(rep(fill, n))
  b <- as.numeric(b)
  if (length(b) == 1L) b <- rep(b, n)
  if (length(b) != n) stop("Bounds must have the same length as u0.", call. = FALSE)
  b
}

.nb_validate <- function(u, lb, ub) {
  if (any(!(u >= lb))) stop("The initial guess must satisfy the lower bounds.", call. = FALSE)
  if (any(!(u <= ub))) stop("The initial guess must satisfy the upper bounds.", call. = FALSE)
  if (any(!(lb <= ub))) stop("Each lower bound must be less than or equal to its upper bound.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Statistics and Jacobians

.nb_stats <- function() {
  st <- new.env(parent = emptyenv())
  st$nf <- 0L; st$njacs <- 0L; st$nfactors <- 0L; st$nsolve <- 0L; st$nsteps <- 0L
  st
}

.nb_stats_list <- function(st) list(nf = st$nf, njacs = st$njacs, nfactors = st$nfactors, nsolve = st$nsolve,
                                    nsteps = st$nsteps)

.nb_dense_ad <- function(ad) if (inherits(ad, "AutoSparse")) ad$dense_ad else ad

# One column of the bound aware difference Jacobian (_box_difference_column).
.nb_diff_column <- function(F, x, lb, ub, f0, i, st) {
  if (lb[i] == ub[i]) return(f0 * 0)
  h <- 6.055454452393343e-06 * max(1, abs(x[i]))  # cbrt(eps())
  plus <- min(x[i] + h, ub[i]); minus <- max(x[i] - h, lb[i])
  if (plus > x[i] && minus < x[i] && min(plus - x[i], x[i] - minus) >= h / 2) {
    xp <- x; xm <- x
    xp[i] <- plus; xm[i] <- minus
    fp <- F(xp); fm <- F(xm)
    st$nf <- st$nf + 2L
    hp <- plus - x[i]; hm <- x[i] - minus
    return((hm / hp * (fp - f0) + hp / hm * (f0 - fm)) / (hp + hm))
  }
  xp <- x
  xp[i] <- if (plus - x[i] >= x[i] - minus) plus else minus
  if (xp[i] == x[i]) return(f0 * 0)
  fp <- F(xp)
  st$nf <- st$nf + 1L
  (fp - f0) / (xp[i] - x[i])
}

# The Jacobian evaluator of a bounded solve. An analytic Jacobian, or the
# automatic differentiation of Julia (central differences here), is
# evaluated once when the cache is built (construct_jacobian_cache) and once
# per call; the bound aware finite differences (BoxDifferenceCache) count
# their residual evaluations.
.nb_jacobian <- function(autodiff, F, J, lb, ub, st, nres) {
  dense <- .nb_dense_ad(autodiff)
  if (is.null(J) && inherits(dense, "NoAutoDiff")) stop(NoAutoDiffSelectedError())
  if (is.null(J) && identical(class(dense)[1L], "AutoFiniteDiff")) {
    return(function(x) {
      f0 <- F(x)
      st$nf <- st$nf + 1L
      Jm <- matrix(0, length(f0), length(x))
      for (i in seq_along(x)) Jm[, i] <- .nb_diff_column(F, x, lb, ub, f0, i, st)
      st$njacs <- st$njacs + 1L
      Jm
    })
  }
  st$njacs <- st$njacs + 1L
  function(x) {
    st$njacs <- st$njacs + 1L
    if (!is.null(J)) {
      Jm <- J(x)
      if (!is.matrix(Jm)) Jm <- matrix(Jm, nrow = nres)
      return(Jm)
    }
    .jac_central(F, x)
  }
}

# ---------------------------------------------------------------------------
# Linear least squares by unpivoted QR, as the default LinearSolve algorithm
# for a tall matrix; a rank deficient factor (smallest |R_ii| at most
# min(m, n) eps times the largest) falls back to the minimum norm solution.

.nb_qr_lsq <- function(A, b) {
  n <- ncol(A)
  q <- qr(A, tol = 0)
  d <- abs(diag(qr.R(q))[seq_len(min(nrow(A), n))])
  if (q$rank == n && length(d) == n && min(d) > min(dim(A)) * .Machine$double.eps * max(d))
    return(as.numeric(qr.coef(q, b)))
  s <- svd(A)
  keep <- s$d > min(dim(A)) * .Machine$double.eps * max(s$d, 0)
  if (!any(keep)) return(numeric(n))
  as.numeric(s$v[, keep, drop = FALSE] %*% (crossprod(s$u[, keep, drop = FALSE], b) / s$d[keep]))
}

# _box_linear_solve!: the scaled and augmented system [J D; sqrt(C)] p = -[b; r].
.nb_linear_solve <- function(cc, Jm, b, scale, diagonal, lower_rhs) {
  A <- rbind(Jm * rep(scale, each = nrow(Jm)), diag(sqrt(diagonal), length(scale)))
  cc$st$nsolve <- cc$st$nsolve + 1L
  cc$st$nfactors <- cc$st$nfactors + 1L
  p <- .nb_qr_lsq(A, -c(b, lower_rhs))
  if (!all(is.finite(p))) {
    cc$retcode <- "InternalLinearSolveFailed"; cc$force_stop <- TRUE
    return(scale * 0)
  }
  p * (scale != 0)
}

# _box_lsq: the least squares step, restricted to the trust radius by the
# two dimensional subspace of the scaled gradient and the step.
.nb_lsq <- function(cc, Jm, b, scale, diagonal, radius = Inf, lower_rhs = scale * 0) {
  p <- .nb_linear_solve(cc, Jm, b, scale, diagonal, lower_rhs)
  if (cc$force_stop || isTRUE(.nb_norm(p) <= radius)) return(p)
  g <- scale * as.numeric(crossprod(Jm, b)) + sqrt(diagonal) * lower_rhs
  q1 <- -g / .nb_norm(g)
  q2 <- p - sum(q1 * p) * q1
  Q <- if (.nb_norm(q2) > sqrt(.Machine$double.eps) * .nb_norm(p)) cbind(q1, q2 / .nb_norm(q2)) else matrix(q1, ncol = 1L)
  AQ <- apply(Q, 2L, function(q) as.numeric(Jm %*% (scale * q)))
  if (!is.matrix(AQ)) AQ <- matrix(AQ, ncol = ncol(Q))
  H <- crossprod(AQ) + crossprod(Q, diagonal * Q)
  eig <- eigen((H + t(H)) / 2, symmetric = TRUE)
  values <- pmax(eig$values, 0)
  rhs <- as.numeric(crossprod(eig$vectors, crossprod(Q, g)))
  step <- function(lambda) -as.numeric(eig$vectors %*% (rhs / (values + lambda)))
  lo <- 0
  hi <- max(1, .nb_norm(rhs) / radius)
  for (k in 1:60) {
    mid <- lo / 2 + hi / 2
    if (.nb_norm(step(mid)) > radius) lo <- mid else hi <- mid
    if (hi - lo <= .Machine$double.eps * max(1, hi)) break
  }
  as.numeric(Q %*% step(hi))
}

# ---------------------------------------------------------------------------
# The NativeBoundedAlgorithm driver

.nb_check <- function(cc, fu, u, uprev) {
  if (.term_check(cc$tc, fu, u, uprev)) {
    cc$retcode <- cc$tc$retcode
    if (isTRUE(cc$tc$mode$best) && !is.null(cc$tc$u) && !identical(cc$tc$u, cc$u)) {
      cc$u <- cc$tc$u; cc$fu <- cc$F(cc$u); cc$st$nf <- cc$st$nf + 1L
    }
    cc$force_stop <- TRUE
  }
  invisible(NULL)
}

.nb_trial <- function(cc, x) {
  u <- .nb_clamp(x, cc$lb, cc$ub)
  fu <- cc$F(u)
  cc$st$nf <- cc$st$nf + 1L
  list(u = u, fu = fu)
}

.nb_accept <- function(cc, u, fu) {
  uprev <- cc$u
  cc$du <- u - uprev
  cc$u <- u; cc$fu <- fu
  .nb_check(cc, fu, u, uprev)
}

# LineSearch.ProjectedBackTracking along P(u + alpha d), with phi0 given.
.nb_linesearch <- function(cc, g, direction) {
  u <- cc$u
  cost <- sum(cc$fu * cc$fu) / 2
  alpha <- 1
  for (k in seq_len(cc$alg$options$max_backtracks)) {
    trial <- .nb_clamp(u + alpha * direction, cc$lb, cc$ub)
    slope <- sum(g * (trial - u))
    if (is.finite(slope) && slope < 0) {
      ft <- cc$F(trial)
      cc$st$nf <- cc$st$nf + 1L
      value <- .nb_norm(ft)^2 / 2
      if (is.finite(value) && value <= cost + 1e-4 * slope) {
        .nb_accept(cc, trial, ft)
        return(TRUE)
      }
    }
    alpha <- alpha * 0.5
  }
  FALSE
}

.nb_trust_update <- function(cc, u, fu, predicted, stepnorm, correction = 0) {
  actual <- (sum(cc$fu * cc$fu) - sum(fu * fu)) / 2
  ratio <- if (predicted > 0 && all(is.finite(fu))) (actual - correction) / predicted else -Inf
  r <- .nb_update_radius(cc$radius, ratio, stepnorm, cc$alg$max_trust_radius,
                         shrink_radius = min(cc$radius / 4, max(stepnorm / 4, .nb_ulp(cc$radius))))
  cc$radius <- r$radius
  accepted <- ratio > 1e-4 && actual > 0
  if (accepted) .nb_accept(cc, u, fu)
  if (!accepted && (stepnorm == 0 || cc$radius <= .Machine$double.eps)) {
    cc$retcode <- "Stalled"; cc$force_stop <- TRUE
  }
  accepted
}

.nb_step <- function(cc) {
  .nb_check(cc, cc$fu, cc$u, cc$u)
  if (cc$force_stop) return(invisible(NULL))
  if (!all(is.finite(cc$fu))) { cc$retcode <- "Unstable"; cc$force_stop <- TRUE; return(invisible(NULL)) }
  Jm <- cc$jac(cc$u)
  if (!all(is.finite(Jm))) { cc$retcode <- "Unstable"; cc$force_stop <- TRUE; return(invisible(NULL)) }
  x <- cc$u; f <- cc$fu
  g <- as.numeric(crossprod(Jm, f))
  if (!all(is.finite(g))) { cc$retcode <- "Unstable"; cc$force_stop <- TRUE; return(invisible(NULL)) }
  if (max(abs(.nb_projected_gradient(x, g, cc$lb, cc$ub))) <= cc$gtol) {
    cc$retcode <- if (cc$ls) "Success" else "Stalled"
    cc$force_stop <- TRUE
    return(invisible(NULL))
  }
  accepted <- switch(cc$alg$method,
                     reflective = .nb_step_reflective(cc, Jm, x, f, g),
                     lm = .nb_step_lm(cc, Jm, x, f, g),
                     dogbox = .nb_step_dogbox(cc, Jm, x, f, g),
                     gauss_newton = .nb_step_gauss_newton(cc, Jm, x, f, g))
  .trace_update(cc$trace, cc$st$nsteps + 1L, cc$u, cc$fu, cc$du, isTRUE(accepted))
  invisible(NULL)
}

.nb_run <- function(alg, F0, J, x0, opts, lb, ub, ls) {
  st <- .nb_stats()
  n <- length(x0)
  x0 <- as.numeric(x0)
  lb <- .nb_bounds(lb, -Inf, n); ub <- .nb_bounds(ub, Inf, n)
  .nb_validate(x0, lb, ub)
  if (!all(is.finite(x0))) stop("u0 must be finite.", call. = FALSE)
  F <- function(x) as.numeric(F0(x))
  cc <- new.env(parent = emptyenv())
  cc$alg <- alg; cc$F <- F; cc$st <- st; cc$lb <- lb; cc$ub <- ub; cc$ls <- ls
  cc$u <- if (alg$method == "reflective") .nb_interior(x0, lb, ub, initial = TRUE) else x0
  cc$fu <- F(cc$u)
  st$nf <- 1L
  cc$jac <- .nb_jacobian(alg$autodiff, F, J, lb, ub, st, length(cc$fu))
  mode <- opts$termination_condition %||% AbsNormTerminationMode(.nb_norm)
  cc$tc <- .term_init(mode, cc$fu, cc$u, opts$abstol %||% .nb_default_tol, opts$reltol %||% .nb_default_tol, ls)
  cc$jac(cc$u)
  cc$du <- cc$u * 0
  cc$trace <- .trace_init(opts, ls, alg$name, cc$u, cc$fu, cc$du)
  cc$radius <- alg$initial_trust_radius
  cc$damping <- alg$options$damping %||% 1
  cc$gtol <- alg$gtol %||% (if (ls) sqrt(.Machine$double.eps) else 0)
  cc$retcode <- "Default"; cc$force_stop <- FALSE
  maxiters <- as.integer(opts$maxiters %||% 1000L)
  nsteps <- 0L
  while (!cc$force_stop && nsteps < maxiters) {
    .nb_step(cc)
    nsteps <- nsteps + 1L
    st$nsteps <- nsteps
  }
  if (cc$retcode == "Default") cc$retcode <- if (nsteps >= maxiters) "MaxIters" else "Success"
  list(u = cc$u, resid = cc$fu, retcode = cc$retcode, iters = nsteps, stats = .nb_stats_list(st),
       trace = .trace_last(cc$trace, cc$fu))
}

.nb_alg <- function(name, method, autodiff, gtol, initial_trust_radius = 1, max_trust_radius = Inf,
                    options = list()) {
  if (!(initial_trust_radius > 0 && is.finite(initial_trust_radius)))
    stop("initial_trust_radius must be positive and finite.", call. = FALSE)
  if (!(max_trust_radius >= initial_trust_radius))
    stop("max_trust_radius must be at least initial_trust_radius.", call. = FALSE)
  if (!is.null(gtol) && !(is.finite(gtol) && gtol >= 0)) stop("gtol must be finite and nonnegative.", call. = FALSE)
  a <- .nl_alg(name, method = method, autodiff = autodiff, gtol = gtol,
               initial_trust_radius = initial_trust_radius, max_trust_radius = max_trust_radius,
               options = options, allows_bounds = TRUE)
  a$run <- function(F, J, x0, opts, lb, ub, ls) .nb_run(a, F, J, x0, opts, lb, ub, ls)
  a
}

# ---------------------------------------------------------------------------
# The four step rules

.nb_reflective_scaling <- function(x, g, lb, ub) {
  v <- rep(1, length(x)); cc <- numeric(length(x))
  for (i in seq_along(x)) {
    if (lb[i] == ub[i]) v[i] <- 0
    else if (g[i] < 0 && is.finite(ub[i])) { v[i] <- ub[i] - x[i]; cc[i] <- -g[i] }
    else if (g[i] > 0 && is.finite(lb[i])) { v[i] <- x[i] - lb[i]; cc[i] <- g[i] }
  }
  list(d = sqrt(v), c = cc)
}

.nb_step_reflective <- function(cc, Jm, x, f, g) {
  lb <- cc$lb; ub <- cc$ub
  sc <- .nb_reflective_scaling(x, g, lb, ub)
  d <- sc$d; cv <- sc$c
  h <- .nb_lsq(cc, Jm, f, d, cv + (d == 0), cc$radius)
  if (cc$force_stop) return(FALSE)
  gh <- d * g
  product <- function(q) as.numeric(Jm %*% (d * q))
  model <- function(q) { pq <- product(q); sum(gh * q) + (sum(pq * pq) + sum(cv * (q * q))) / 2 }
  pg <- max(abs(.nb_projected_gradient(x, g, lb, ub)))
  theta <- min(1 - 2^-53, max(0.995, 1 - pg))
  hit <- .nb_limit(x, d * h, lb, ub)
  best <- if (hit >= 1) h else theta * hit * h
  bestvalue <- model(best)
  if (hit < 1) {
    at_hit <- hit * h
    direction <- h
    point <- x + d * at_hit
    for (i in seq_along(x)) {
      t <- if (h[i] > 0 && d[i] > 0) (ub[i] - x[i]) / (d[i] * h[i])
           else if (h[i] < 0 && d[i] > 0) (lb[i] - x[i]) / (d[i] * h[i]) else Inf
      if (.nb_isapprox(t, hit, 8 * .Machine$double.eps)) direction[i] <- -direction[i]
    }
    aa <- sum(direction * direction); ab <- sum(at_hit * direction)
    sphere <- (sqrt(max(0, ab * ab + aa * (cc$radius * cc$radius - sum(at_hit * at_hit)))) - ab) / aa
    limit <- min(sphere, .nb_limit(.nb_clamp(point, lb, ub), d * direction, lb, ub))
    pa <- product(at_hit); pd <- product(direction)
    slope <- sum(gh * direction) + sum(pa * pd) + sum((cv * at_hit) * direction)
    curvature <- sum(pd * pd) + sum(cv * (direction * direction))
    t <- .nb_cauchy_length(slope, curvature, limit)
    reflected <- theta * (at_hit + t * direction)
    value <- model(reflected)
    if (value < bestvalue) { best <- reflected; bestvalue <- value }
  }
  direction <- -gh
  limit <- min(cc$radius / .nb_norm(direction), theta * .nb_limit(x, d * direction, lb, ub))
  pd <- product(direction)
  curvature <- sum(pd * pd) + sum(cv * (direction * direction))
  t <- .nb_cauchy_length(-sum(gh * gh), curvature, limit)
  cauchy <- t * direction
  if (model(cauchy) < bestvalue) best <- cauchy
  trial <- .nb_interior(x + d * best, lb, ub)
  # the model of the representable step after the strictly feasible rounding
  s <- trial - x
  actual_h <- ifelse(d == 0, s * 0, s / d)
  tr <- .nb_trial(cc, trial)
  correction <- sum(cv * (actual_h * actual_h)) / 2
  .nb_trust_update(cc, tr$u, tr$fu, -model(actual_h), .nb_norm(actual_h), correction)
}

.nb_constrained_lsq <- function(cc, Jm, f, lo, hi, damping) {
  n <- length(lo)
  s <- .nb_clamp(numeric(n), lo, hi)
  active <- lo == hi
  for (k in seq_len(10L * n + 100L)) {
    free <- !active
    candidate <- s
    if (any(free)) {
      scale <- as.numeric(free)
      diagonal <- damping * scale + !free
      candidate <- candidate + .nb_lsq(cc, Jm, f + as.numeric(Jm %*% s), scale, diagonal,
                                       lower_rhs = sqrt(damping) * s * free)
      if (cc$force_stop) return(s)
    }
    direction <- candidate - s
    alpha <- min(1, .nb_limit(s, direction, lo, hi))
    s <- .nb_clamp(s + alpha * direction, lo, hi)
    if (alpha < 1) {
      for (i in seq_len(n))
        if ((direction[i] < 0 && s[i] <= lo[i]) || (direction[i] > 0 && s[i] >= hi[i])) active[i] <- TRUE
      next
    }
    gradient <- as.numeric(crossprod(Jm, as.numeric(Jm %*% s) + f)) + damping * s
    violation <- vapply(seq_len(n), function(i) {
      if (!active[i] || lo[i] == hi[i]) 0
      else if (s[i] <= lo[i]) max(-gradient[i], 0) else max(gradient[i], 0)
    }, 0)
    index <- which.max(violation)
    if (violation[index] <= sqrt(.Machine$double.eps) * max(1, max(abs(gradient)))) break
    active[index] <- FALSE
  }
  s
}

.nb_step_lm <- function(cc, Jm, x, f, g) {
  lambda <- cc$damping * max(.nb_norm(f), .Machine$double.eps)
  s <- .nb_constrained_lsq(cc, Jm, f, cc$lb - x, cc$ub - x, lambda)
  if (cc$force_stop) return(FALSE)
  cost <- sum(f * f) / 2
  accepted <- .nb_linesearch(cc, g, s)
  if (!accepted) accepted <- .nb_linesearch(cc, g, -g)
  if (accepted) {
    step <- cc$u - x
    predicted <- -.nb_model(Jm, g, step)
    ratio <- if (predicted > 0) (cost - sum(cc$fu * cc$fu) / 2) / predicted else 0
    if (ratio > 0.75) cc$damping <- max(cc$damping / 2, .Machine$double.eps)
    else if (ratio < 0.25) cc$damping <- cc$damping * 2
  } else {
    cc$retcode <- "Stalled"; cc$force_stop <- TRUE
  }
  accepted
}

.nb_dogbox_segment <- function(x, direction, alpha, lo, hi) {
  trial <- .nb_clamp(x + alpha * direction, lo, hi)
  for (i in seq_along(x)) {
    if (direction[i] == 0) next
    bound <- if (direction[i] < 0) lo[i] else hi[i]
    hit <- (bound - x[i]) / direction[i]
    near <- is.finite(bound) && abs(trial[i] - bound) <= .nb_ulp(max(abs(x[i]), abs(bound)))
    if (alpha >= hit || .nb_isapprox(alpha, hit, 8 * .Machine$double.eps) || near) trial[i] <- bound
  }
  trial
}

.nb_step_dogbox <- function(cc, Jm, x, f, g) {
  free <- !.nb_active(x, g, cc$lb, cc$ub)
  gn <- .nb_lsq(cc, Jm, f, as.numeric(free), as.numeric(!free))
  if (cc$force_stop) return(FALSE)
  lo <- pmax(cc$lb, x - cc$radius)
  hi <- pmin(cc$ub, x + cc$radius)
  if (all(lo <= x + gn & x + gn <= hi)) {
    trial <- x + gn
  } else {
    direction <- -g * free
    Jd <- as.numeric(Jm %*% direction)
    alpha <- .nb_cauchy_length(sum(g * direction), sum(Jd * Jd), .nb_limit(x, direction, lo, hi))
    cauchy <- .nb_dogbox_segment(x, direction, alpha, lo, hi)
    toward <- x + gn - cauchy
    beta <- min(1, .nb_limit(cauchy, toward, lo, hi))
    trial <- .nb_dogbox_segment(cauchy, toward, beta, lo, hi)
  }
  tr <- .nb_trial(cc, trial)
  step <- tr$u - x
  .nb_trust_update(cc, tr$u, tr$fu, -.nb_model(Jm, g, step), max(abs(step)))
}

.nb_step_gauss_newton <- function(cc, Jm, x, f, g) {
  free <- !.nb_active(x, g, cc$lb, cc$ub)
  kind <- cc$alg$options$globalization
  radius <- if (kind == "linesearch") Inf else cc$radius
  step <- .nb_lsq(cc, Jm, f, as.numeric(free), as.numeric(!free), radius)
  if (cc$force_stop) return(FALSE)
  if (kind == "linesearch") {
    accepted <- .nb_linesearch(cc, g, step)
  } else {
    step <- .nb_clamp(x + step, cc$lb, cc$ub) - x
    direction <- -g * free
    Jd <- as.numeric(Jm %*% direction)
    alpha <- .nb_cauchy_length(sum(g * direction), sum(Jd * Jd), radius / .nb_norm(direction))
    cauchy <- .nb_clamp(x + alpha * direction, cc$lb, cc$ub) - x
    if (.nb_model(Jm, g, cauchy) < .nb_model(Jm, g, step)) step <- cauchy
    tr <- .nb_trial(cc, x + step)
    step <- tr$u - x
    accepted <- .nb_trust_update(cc, tr$u, tr$fu, -.nb_model(Jm, g, step), .nb_norm(step))
    if (kind == "trustregion") return(accepted)
  }
  if (!accepted && !cc$force_stop) {
    accepted <- .nb_linesearch(cc, g, -g)
    if (!accepted) { cc$retcode <- "Stalled"; cc$force_stop <- TRUE }
  }
  accepted
}

# ---------------------------------------------------------------------------
# Constructors

#' @rdname bounded_nonlinear
#' @export
TrustRegionReflective <- function(autodiff = NULL, jvp_autodiff = NULL, vjp_autodiff = NULL, linsolve = NULL,
                                  concrete_jac = NULL, gtol = NULL, initial_trust_radius = 1,
                                  max_trust_radius = Inf)
  .nb_alg("TrustRegionReflective", "reflective", autodiff, gtol, initial_trust_radius, max_trust_radius)

#' @rdname bounded_nonlinear
#' @export
BoundedLevenbergMarquardt <- function(autodiff = NULL, jvp_autodiff = NULL, vjp_autodiff = NULL, linsolve = NULL,
                                      concrete_jac = NULL, gtol = NULL, damping = 1e-3, max_backtracks = 40L) {
  if (!(is.finite(damping) && damping > 0)) stop("damping must be positive and finite.", call. = FALSE)
  if (!(max_backtracks > 0)) stop("max_backtracks must be positive.", call. = FALSE)
  .nb_alg("BoundedLevenbergMarquardt", "lm", autodiff, gtol,
          options = list(damping = damping, max_backtracks = as.integer(max_backtracks)))
}

#' @rdname bounded_nonlinear
#' @export
Dogbox <- function(autodiff = NULL, jvp_autodiff = NULL, vjp_autodiff = NULL, linsolve = NULL, concrete_jac = NULL,
                   gtol = NULL, initial_trust_radius = 1, max_trust_radius = Inf)
  .nb_alg("Dogbox", "dogbox", autodiff, gtol, initial_trust_radius, max_trust_radius)

#' @rdname bounded_nonlinear
#' @export
BoundedGaussNewton <- function(globalization = "linesearch", autodiff = NULL, jvp_autodiff = NULL,
                               vjp_autodiff = NULL, linsolve = NULL, concrete_jac = NULL, gtol = NULL,
                               initial_trust_radius = 1, max_trust_radius = Inf, max_backtracks = 40L) {
  globalization <- sub("^:", "", as.character(globalization))
  if (!globalization %in% c("linesearch", "trustregion", "trustregion_linesearch"))
    stop("globalization must be :linesearch, :trustregion, or :trustregion_linesearch.", call. = FALSE)
  if (!(max_backtracks > 0)) stop("max_backtracks must be positive.", call. = FALSE)
  .nb_alg("BoundedGaussNewton", "gauss_newton", autodiff, gtol, initial_trust_radius, max_trust_radius,
          options = list(globalization = globalization, max_backtracks = as.integer(max_backtracks)))
}

# ---------------------------------------------------------------------------
# BoundedTrustRegion: GeneralizedFirstOrderAlgorithm with the Dogleg descent
# and the BoundedTrustRegionScheme

#' @rdname bounded_nonlinear
#' @export
BoundedTrustRegion <- function(concrete_jac = NULL, linsolve = NULL, max_trust_radius = 0, initial_trust_radius = 0,
                               step_threshold = 1e-4, shrink_threshold = 1 / 4, expand_threshold = 3 / 4,
                               shrink_factor = 1 / 4, expand_factor = 2, max_shrink_times = 32L, gtol = NULL,
                               autodiff = NULL, vjp_autodiff = NULL, jvp_autodiff = NULL, jacobian_reuse = NULL) {
  a <- .nl_alg("BoundedTrustRegion", autodiff = autodiff, gtol = gtol, max_trust_radius = max_trust_radius,
               jacobian_reuse = .jr_normalize(jacobian_reuse),
               initial_trust_radius = initial_trust_radius, step_threshold = step_threshold,
               shrink_threshold = shrink_threshold, expand_threshold = expand_threshold,
               shrink_factor = shrink_factor, expand_factor = expand_factor,
               max_shrink_times = as.integer(max_shrink_times), allows_bounds = TRUE)
  a$run <- function(F, J, x0, opts, lb, ub, ls) .nb_btr_run(a, F, J, x0, opts, lb, ub, ls)
  a
}

# A linear solve of the Newton descent: LU for a square Jacobian, QR for a
# rectangular one (the minimum norm solution when singular).
.nb_newton_solve <- function(Jm, fu) {
  if (nrow(Jm) == ncol(Jm)) {
    x <- tryCatch(solve(Jm, fu), error = function(e) NULL)
    if (!is.null(x) && all(is.finite(x))) return(as.numeric(x))
  }
  .nb_qr_lsq(Jm, fu)
}

# The Dogleg descent (NewtonDescent and SteepestDescent with the L2 norm).
.nb_dogleg <- function(Jm, fu, radius, newton) {
  if (.nb_norm(newton) <= radius) return(newton)
  cauchy <- -as.numeric(crossprod(Jm, fu))
  l_grad <- .nb_norm(cauchy)
  Jc <- as.numeric(Jm %*% cauchy)
  d_cauchy <- l_grad^3 / sum(Jc * Jc)
  if (d_cauchy >= radius) return((radius / l_grad) * cauchy)
  d1 <- (d_cauchy / l_grad) * cauchy
  d2 <- newton - d1
  a <- sum(d2 * d2); b <- 2 * sum(d1 * d2); c <- d_cauchy^2 - radius^2
  tau <- (-b + sqrt(max(0, b^2 - 4 * a * c))) / (2 * a)
  d1 + tau * d2
}

.nb_btr_run <- function(alg, F0, J, x0, opts, lb, ub, ls) {
  st <- .nb_stats()
  n <- length(x0)
  u <- as.numeric(x0)
  lb <- .nb_bounds(lb, -Inf, n); ub <- .nb_bounds(ub, Inf, n)
  .nb_validate(u, lb, ub)
  F <- function(x) as.numeric(F0(x))
  fu <- F(u)
  tc <- .nlf_term(utils::modifyList(opts, list(abstol = opts$abstol %||% .nb_default_tol,
                                               reltol = opts$reltol %||% .nb_default_tol)), fu, u, ls)
  jac <- .nb_jacobian(alg$autodiff, F, J, lb, ub, st, length(fu))
  max_radius <- if (alg$max_trust_radius == 0) Inf else alg$max_trust_radius
  radius <- min(if (alg$initial_trust_radius == 0) max(.nb_norm(u), 1) else alg$initial_trust_radius, max_radius)
  gtol <- if (is.null(alg$gtol)) (if (ls) sqrt(.Machine$double.eps)) else alg$gtol
  shrink_counter <- 0L
  make_new <- TRUE
  Jm <- NULL
  newton <- NULL
  # JacobianReuseCache (jacobian_reuse.jl): age of the Jacobian and the
  # residual norm when it was computed
  policy <- .jr_resolve(alg$jacobian_reuse %||% .jr_normalize(NULL), n)
  reuses <- policy$max_age > 1L
  jr_age <- 0L
  jr_rnorm <- .nb_norm(fu)
  stale <- function() reuses && jr_age > 0L
  trace <- .trace_init(opts, ls, alg$name, u, fu, u * 0, jac0 = function() if (is.null(J)) .jac_central(F, u) else J(u))
  retcode <- "Default"; force_stop <- FALSE
  u_cache <- u
  maxiters <- as.integer(opts$maxiters %||% 1000L)
  nsteps <- 0L
  check <- function(fu, u, uprev) {
    if (.term_check(tc, fu, u, uprev)) {
      retcode <<- tc$retcode
      force_stop <<- TRUE
      rollback()
    }
  }
  rollback <- function() {
    if (isTRUE(tc$mode$best) && !is.null(tc$u) && !identical(tc$u, u)) {
      u <<- tc$u; fu <<- F(u); st$nf <- st$nf + 1L
    }
  }
  while (!force_stop && nsteps < maxiters) {
    nsteps <- nsteps + 1L
    st$nsteps <- nsteps
    policy_driven <- TRUE
    repeat {
      new_jacobian <- if (policy_driven) make_new else TRUE
      if (new_jacobian) {
        Jm <- jac(u)
        jr_age <- 0L
        if (reuses) jr_rnorm <- .nb_norm(fu)
      }
      gradient <- as.numeric(crossprod(Jm, fu))
      if (!is.null(gtol) && max(abs(.nb_projected_gradient(u, gradient, lb, ub))) <= gtol) {
        retcode <- if (ls) "Success" else "Stalled"; force_stop <- TRUE
        break
      }
      # NewtonDescent: the factorization is reused with a stale Jacobian
      st$nsolve <- st$nsolve + 1L
      if (new_jacobian) st$nfactors <- st$nfactors + 1L
      newton <- -.nb_newton_solve(Jm, fu)
      if (all(is.finite(newton))) break
      if (new_jacobian) { retcode <- "InternalLinearSolveFailed"; force_stop <- TRUE; break }
      # a failed solve with a stale Jacobian is retried with a new one
      make_new <- TRUE; policy_driven <- FALSE
    }
    if (force_stop) break
    du <- .nb_dogleg(Jm, fu, radius, newton)
    # BoundedTrustRegionScheme
    dogleg_step <- .nb_clamp(u + du, lb, ub) - u
    Jd <- as.numeric(Jm %*% dogleg_step)
    dogleg_model <- .nb_model_value(Jd, dogleg_step, gradient)
    cdir <- ifelse(.nb_active(u, gradient, lb, ub), 0, -gradient)
    dnorm <- .nb_norm(cdir)
    if (dnorm == 0) {
      cauchy_step <- u * 0; Jc <- Jd * 0; cauchy_model <- 0
    } else {
      Jc <- as.numeric(Jm %*% cdir)
      alpha <- .nb_cauchy_length(sum(cdir * gradient), sum(Jc * Jc), radius / dnorm)
      cauchy_step <- .nb_clamp(u + alpha * cdir, lb, ub) - u
      Jc <- as.numeric(Jm %*% cauchy_step)
      cauchy_model <- .nb_model_value(Jc, cauchy_step, gradient)
    }
    if (cauchy_model < dogleg_model) {
      dogleg_step <- cauchy_step
      predicted <- -cauchy_model
    } else predicted <- -dogleg_model
    u_new <- u + dogleg_step
    fu_new <- F(u_new)
    st$nf <- st$nf + 1L
    actual <- (.nb_norm(fu)^2 - .nb_norm(fu_new)^2) / 2
    rho <- if (predicted > 0 && all(is.finite(fu_new))) actual / predicted else -1
    accepted <- rho > alg$step_threshold && actual > 0
    r <- .nb_update_radius(radius, rho, .nb_norm(dogleg_step), max_radius, alg$shrink_threshold,
                           alg$expand_threshold, radius * alg$shrink_factor, alg$expand_factor)
    radius <- r$radius
    shrink_counter <- if (r$shrunk) shrink_counter + 1L else 0L
    if (accepted) {
      u <- u_new; fu <- fu_new
      if (!is.null(opts$postcondition)) { u <- opts$postcondition(u, u_cache); fu <- F(u); st$nf <- st$nf + 1L }
    }
    if (shrink_counter > alg$max_shrink_times) { retcode <- "ShrinkThresholdExceeded"; force_stop <- TRUE }
    if (accepted) {
      # prepare_next_jacobian!
      if (!reuses) make_new <- TRUE
      else {
        rn <- .nb_norm(fu)
        jr_age <- jr_age + 1L
        improved <- is.finite(rn) && is.finite(jr_rnorm) && rn < policy$max_residual_ratio * jr_rnorm
        jr_rnorm <- rn
        make_new <- !(improved && jr_age < policy$max_age)
      }
    } else make_new <- stale()
    check(fu, u, u_cache)
    .trace_update(trace, nsteps, u, fu, du, accepted)
    u_cache <- u
  }
  if (retcode == "Default") retcode <- if (nsteps >= maxiters) "MaxIters" else "Success"
  rollback()
  list(u = u, resid = fu, retcode = retcode, iters = nsteps, stats = .nb_stats_list(st),
       trace = .trace_last(trace, fu))
}
