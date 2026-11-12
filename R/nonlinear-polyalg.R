# Polyalgorithms, multistart and solver options of NonlinearSolveBase
# (v2.55.1: polyalg.jl and the generated __solve of solve.jl, tracing.jl,
# verbosity.jl, conditioning.jl) and NonlinearSolveFirstOrder (v2.10.0:
# poly_algs.jl FastShortcutBoundedPolyalg, multistart.jl SobolMultistart,
# jacobian_reuse.jl, eisenstat_walker.jl), with the Sobol sequence of
# Sobol.jl (src/Sobol.jl).

# ---------------------------------------------------------------------------
# Statistics shared by the stages of a solve (SciMLBase.NLStats, merged by
# adding the fields)

.nl_stats_fields <- c("nf", "njacs", "nfactors", "nsolve", "nsteps")

.nl_stats_add <- function(a, b) {
  out <- a
  for (k in .nl_stats_fields) out[[k]] <- as.integer((a[[k]] %||% 0L) + (b[[k]] %||% 0L))
  out
}

.nl_stats_of <- function(res) {
  st <- res$stats %||% list()
  if (is.null(st$nsteps)) st$nsteps <- as.integer(res$iters %||% 0L)
  st
}

.nl_resid_norm <- function(r, ls) {
  v <- if (ls) sqrt(sum(r * r)) else max(abs(r), 0)
  if (is.na(v)) Inf else v
}

# ---------------------------------------------------------------------------
# NonlinearSolvePolyAlgorithm

#' Polyalgorithms and multistart for nonlinear problems
#'
#' `NonlinearSolvePolyAlgorithm(algs)` tries the algorithms in `algs` in
#' turn, each from the initial guess, and returns the first successful
#' solution (`"StalledSuccess"` does not count); when none succeeds it
#' returns the one with the smallest residual (infinity norm for
#' `NonlinearProblem`, 2-norm for least squares). The statistics in
#' `sol$stats` add up the work of every stage tried. When every stage handles
#' bounds natively an initial guess outside the box is projected onto it.
#'
#' `FastShortcutBoundedPolyalg()` is the polyalgorithm of NonlinearSolve for
#' problems with bounds: [BoundedTrustRegion()], then
#' [BoundedGaussNewton()] with its projected line search (only the first
#' with `must_support_postcondition = TRUE`).
#'
#' `SobolMultistart(alg)` runs `alg` from `nstarts` start points: the initial
#' guess, then points of a Sobol sequence over the box (components with an
#' infinite bound are sampled `search_scale * max(1, |u0_i|)` around
#' `u0_i`). With `early_exit` the first successful solve is returned;
#' otherwise every start runs and the feasible solution with the smallest
#' residual 2-norm is returned. With `restoration`, a start that stalls on an
#' active bound is continued without the bounds that are not fixed: a root
#' found inside the box replaces the stalled result, and a root outside it
#' shows that the stalled point is a genuine constrained solution (a warning
#' says so). When no start succeeds the retcode is `"Stalled"`. This is the
#' default of NonlinearSolve for problems with `lb` or `ub`, wrapping
#' `FastShortcutBoundedPolyalg()`. The Sobol points (Joe and Kuo direction
#' numbers, as Sobol.jl) are available for up to 1024 unknowns.
#'
#' @param algs A list of algorithms.
#' @param start_index The stage to start at.
#' @param store_original Keep the solution of the winning stage in
#'   `sol$original`.
#' @param concrete_jac,linsolve,autodiff,jvp_autodiff,vjp_autodiff,gtol
#'   Passed to the stages.
#' @param must_support_postcondition Only use stages that apply the
#'   `postcondition` corrector.
#' @param alg The local algorithm of `SobolMultistart()`.
#' @param nstarts The number of start points, the initial guess included.
#' @param search_scale Half width factor of the sampling window of unbounded
#'   components.
#' @param early_exit Stop at the first successful start.
#' @param restoration Continue stalled starts without the bounds.
#' @param restoration_alg The algorithm of the restoration probes (default
#'   `alg`).
#' @return An algorithm object.
#' @examples
#' f <- function(u, p) c(10 * (u[2] - u[1] * u[1]), 1 - u[1])
#' prob <- NonlinearProblem(f, u0 = c(-1.2, 1), lb = c(-2, -2), ub = c(2, 2))
#' sol <- solve(prob, NonlinearSolvePolyAlgorithm(list(Dogbox(), TrustRegionReflective())))
#' sol$u; sol$stats$nf
#' solve(prob, FastShortcutBoundedPolyalg())$retcode
#' solve(prob, SobolMultistart(nstarts = 4))$u
#' @name nonlinear_polyalgorithms
NULL

#' @rdname nonlinear_polyalgorithms
#' @export
NonlinearSolvePolyAlgorithm <- function(algs, start_index = 1L, store_original = FALSE) {
  if (inherits(algs, "NonlinearAlgorithm")) algs <- list(algs)
  if (!(start_index >= 1 && start_index <= length(algs))) stop("start_index must be in 1:length(algs)", call. = FALSE)
  a <- .nl_alg("NonlinearSolvePolyAlgorithm", stages = algs, start_index = as.integer(start_index),
               store_original = isTRUE(store_original),
               allows_bounds = all(vapply(algs, function(z) isTRUE(z$allows_bounds), TRUE)),
               supports_postcondition = all(vapply(algs, .nl_supports_post, TRUE)))
  a$run <- function(F, J, x0, opts, lb, ub, ls) .nl_poly_run(a, F, J, x0, opts, lb, ub, ls)
  a
}

.nl_supports_post <- function(alg) {
  if (!is.null(alg$supports_postcondition)) return(isTRUE(alg$supports_postcondition))
  class(alg)[1L] %in% c("BoundedTrustRegion", "GeneralizedFirstOrderAlgorithm", "NewtonRaphson", "TrustRegion",
                        "GaussNewton", "LevenbergMarquardt")
}

.nl_box <- function(x, lb, ub) {
  n <- length(x)
  if (is.null(lb) && is.null(ub)) return(x)
  .nb_clamp(x, .nb_bounds(lb, -Inf, n), .nb_bounds(ub, Inf, n))
}

.nl_poly_run <- function(alg, F, J, x0, opts, lb, ub, ls) {
  if (isTRUE(alg$allows_bounds)) x0 <- .nl_box(as.numeric(x0), lb, ub)
  stats <- list()
  results <- list()
  for (i in seq.int(alg$start_index, length(alg$stages))) {
    res <- .nl_run(alg$stages[[i]], F, J, x0, opts, lb, ub, ls)
    stats <- .nl_stats_add(stats, .nl_stats_of(res))
    res$alg <- res$alg %||% alg$stages[[i]]$name
    results[[length(results) + 1L]] <- res
    if (successful_retcode(res$retcode) && res$retcode != "StalledSuccess") {
      res$stats <- stats; res$iters <- stats$nsteps
      if (alg$store_original) res$original <- res
      return(res)
    }
  }
  norms <- vapply(results, function(r) .nl_resid_norm(r$resid, ls), 0)
  best <- results[[which.min(norms)]]
  best$stats <- stats; best$iters <- stats$nsteps
  best
}

#' @rdname nonlinear_polyalgorithms
#' @export
FastShortcutBoundedPolyalg <- function(concrete_jac = NULL, linsolve = NULL, autodiff = NULL, jvp_autodiff = NULL,
                                       vjp_autodiff = NULL, gtol = NULL, must_support_postcondition = FALSE) {
  first <- BoundedTrustRegion(concrete_jac = concrete_jac, linsolve = linsolve, autodiff = autodiff,
                              jvp_autodiff = jvp_autodiff, vjp_autodiff = vjp_autodiff, gtol = gtol)
  algs <- if (isTRUE(must_support_postcondition)) list(first)
          else list(first, BoundedGaussNewton(concrete_jac = concrete_jac, linsolve = linsolve, autodiff = autodiff,
                                              jvp_autodiff = jvp_autodiff, vjp_autodiff = vjp_autodiff, gtol = gtol))
  NonlinearSolvePolyAlgorithm(algs)
}

# ---------------------------------------------------------------------------
# The Sobol sequence of Sobol.jl, with unsigned 32 bit arithmetic on doubles

.u32 <- function(x) x %% 4294967296

.xor32 <- function(a, b) {
  ah <- a %/% 65536; al <- a %% 65536
  bh <- b %/% 65536; bl <- b %% 65536
  bitwXor(as.integer(ah), as.integer(bh)) * 65536 + bitwXor(as.integer(al), as.integer(bl))
}

.shl32 <- function(a, k) .u32(a * 2^k)

.trailing_zeros <- function(n) {
  k <- 0L
  while (n %% 2 == 0) { n <- n %/% 2; k <- k + 1L }
  k
}

.sobol_new <- function(N) {
  if (N < 0 || N > length(.sobol_a) + 1L)
    stop("Sobol direction numbers are available for up to ", length(.sobol_a) + 1L, " dimensions", call. = FALSE)
  m <- matrix(1, N, 32L)
  offset <- 0L
  for (i in seq_len(N)[-1L]) {
    a <- .sobol_a[i - 1L]
    d <- as.integer(floor(log2(a)))
    m[i, seq_len(d)] <- .sobol_minit_flat[offset + seq_len(d)]
    offset <- offset + d
    for (j in seq.int(d + 1L, 32L)) {
      ac <- a
      m[i, j] <- m[i, j - d]
      for (k in 0:(d - 1L)) {
        if (ac %% 2 == 1) m[i, j] <- .xor32(m[i, j], .shl32(m[i, j - d + k], d - k))
        ac <- ac %/% 2
      }
    }
  }
  s <- new.env(parent = emptyenv())
  s$m <- m; s$x <- numeric(N); s$b <- integer(N); s$n <- 0
  s
}

.sobol_next <- function(s) {
  s$n <- s$n + 1
  c <- .trailing_zeros(s$n)
  N <- length(s$x)
  out <- numeric(N)
  for (i in seq_len(N)) {
    b <- s$b[i]
    if (b >= c) {
      s$x[i] <- .xor32(s$x[i], .shl32(s$m[i, c + 1L], b - c))
      out[i] <- s$x[i] * 2^(-b - 1L)
    } else {
      s$x[i] <- .xor32(.shl32(s$x[i], c - b), s$m[i, c + 1L])
      s$b[i] <- c
      out[i] <- s$x[i] * 2^(-c - 1L)
    }
  }
  out
}

# ---------------------------------------------------------------------------
# SobolMultistart

#' @rdname nonlinear_polyalgorithms
#' @export
SobolMultistart <- function(alg = FastShortcutBoundedPolyalg(), nstarts = 16L, search_scale = 10, early_exit = TRUE,
                            restoration = TRUE, restoration_alg = NULL) {
  if (!(nstarts >= 1)) stop("`nstarts` must be at least 1.", call. = FALSE)
  if (!(search_scale > 0)) stop("`search_scale` must be positive.", call. = FALSE)
  a <- .nl_alg("SobolMultistart", alg = alg, restoration_alg = restoration_alg %||% alg, nstarts = as.integer(nstarts),
               search_scale = search_scale, early_exit = isTRUE(early_exit), restoration = isTRUE(restoration),
               allows_bounds = TRUE, supports_postcondition = .nl_supports_post(alg))
  a$run <- function(F, J, x0, opts, lb, ub, ls) .nl_multistart_run(a, F, J, x0, opts, lb, ub, ls)
  a
}

.ms_tol <- 2^-39  # eps()^(3/4)

.ms_feasible <- function(x, lb, ub) {
  for (i in seq_along(x)) {
    slack <- .ms_tol * max(1, abs(lb[i]), abs(ub[i]), abs(x[i]))
    if (x[i] < lb[i] - slack || x[i] > ub[i] + slack) return(FALSE)
  }
  TRUE
}

.ms_on_active_bound <- function(x, lb, ub) {
  for (i in seq_along(x)) {
    if (lb[i] == ub[i]) next
    slack <- .ms_tol * max(1, abs(lb[i]), abs(ub[i]), abs(x[i]))
    if (is.finite(lb[i]) && x[i] - lb[i] <= slack) return(TRUE)
    if (is.finite(ub[i]) && ub[i] - x[i] <= slack) return(TRUE)
  }
  FALSE
}

.nl_multistart_run <- function(alg, F, J, x0, opts, lb0, ub0, ls) {
  x0 <- as.numeric(x0)
  n <- length(x0)
  lb <- .nb_bounds(lb0, -Inf, n); ub <- .nb_bounds(ub0, Inf, n)
  halfwidth <- alg$search_scale * pmax(1, abs(x0))
  lo <- ifelse(is.finite(lb), lb, x0 - halfwidth)
  hi <- ifelse(is.finite(ub), ub, x0 + halfwidth)
  lo2 <- pmin(lo, hi); hi2 <- pmax(lo, hi)
  starts <- list(.nb_clamp(x0, lb, ub))
  if (alg$nstarts > 1L) {
    sobol <- .sobol_new(n)
    for (i in 2:alg$nstarts) starts[[i]] <- lo2 + (hi2 - lo2) * .sobol_next(sobol)
  }
  stats <- list(); best <- NULL; best_norm <- Inf; first_error <- NULL; warned <- FALSE
  finish <- function(sol, retcode = sol$retcode) {
    sol$retcode <- retcode; sol$stats <- stats; sol$iters <- stats$nsteps %||% 0L
    sol$original <- sol
    sol
  }
  for (x in starts) {
    sol <- tryCatch(.nl_run(alg$alg, F, J, x, opts, lb0, ub0, ls), error = function(e) e)
    if (inherits(sol, "error")) { if (is.null(first_error)) first_error <- sol; next }
    stats <- .nl_stats_add(stats, .nl_stats_of(sol))
    if (alg$restoration && identical(sol$retcode, "Stalled") && !is.null(sol$resid) &&
        sqrt(sum(sol$resid * sol$resid)) != 0 && .ms_on_active_bound(sol$u, lb, ub)) {
      fixed <- lb == ub
      plb <- if (any(fixed)) ifelse(fixed, lb, -Inf) else NULL
      pub <- if (any(fixed)) ifelse(fixed, ub, Inf) else NULL
      probe <- tryCatch(.nl_run(alg$restoration_alg, F, J, sol$u, opts, plb, pub, ls), error = function(e) NULL)
      if (!is.null(probe)) {
        stats <- .nl_stats_add(stats, .nl_stats_of(probe))
        if (successful_retcode(probe$retcode)) {
          if (.ms_feasible(probe$u, lb, ub)) sol <- probe
          else if (!warned) {
            warning("SobolMultistart restoration converged to a root outside the bounds at u = ",
                    paste(format(probe$u), collapse = ", "), "; the stalled point is a genuine constrained solution.",
                    call. = FALSE)
            warned <- TRUE
          }
        }
      }
    }
    if (alg$early_exit && successful_retcode(sol$retcode)) return(finish(sol))
    rn <- if (.ms_feasible(sol$u, lb, ub)) sqrt(sum(sol$resid * sol$resid)) else Inf
    if (isTRUE(rn < best_norm)) { best <- sol; best_norm <- rn }
  }
  if (is.null(best)) {
    if (!is.null(first_error)) stop(first_error)
    u <- .nb_clamp(x0, lb, ub)
    return(list(u = u, resid = F(u), retcode = "Stalled", iters = stats$nsteps %||% 0L, stats = stats))
  }
  if (successful_retcode(best$retcode)) finish(best) else finish(best, "Stalled")
}

# ---------------------------------------------------------------------------
# JacobianReuse and the Eisenstat Walker forcing term

#' Jacobian reuse and forcing terms
#'
#' `JacobianReuse(max_age, max_residual_ratio)` is the policy of the first
#' order solvers for reusing a Jacobian over accepted steps (a modified
#' Newton method): the Jacobian is refreshed when it has served `max_age`
#' accepted steps, when the residual norm fails to drop below
#' `max_residual_ratio` times the previous one, or when a step fails with the
#' stale Jacobian. `max_age = 0` (or 1) recomputes it at every step. The
#' default of the solvers (`jacobian_reuse = NULL`) reuses Jacobians only for
#' problems with at least 16 unknowns; [BoundedTrustRegion()] applies it.
#'
#' `EisenstatWalkerForcing2(eta_0, eta_max, gamma, alpha, safeguard,
#' safeguard_threshold)` is choice 2 of Eisenstat and Walker (1996) for the
#' relative tolerance of the inner linear solves of a Newton Krylov method,
#' `eta_k = gamma * (|r_k| / |r_(k-1)|)^alpha`, safeguarded by `gamma *
#' eta_(k-1)^alpha` when that exceeds `safeguard_threshold` and clamped to
#' `[0, eta_max]`. The linear systems of the R solvers are solved directly,
#' to full accuracy, so the forcing term records the tolerance sequence
#' without loosening any solve.
#'
#' @param max_age The number of accepted steps a Jacobian may serve.
#' @param max_residual_ratio The required residual reduction factor.
#' @param eta_0,eta_max,gamma,alpha,safeguard,safeguard_threshold Settings of
#'   the forcing term.
#' @return A `JacobianReuse` or `EisenstatWalkerForcing2` object.
#' @examples
#' JacobianReuse(max_age = 5)
#' EisenstatWalkerForcing2(eta_0 = 0.1)
#' @name jacobian_reuse
NULL

#' @rdname jacobian_reuse
#' @export
JacobianReuse <- function(max_age = 10L, max_residual_ratio = 0.1) {
  if (!(max_age >= -1)) stop("`max_age` must be nonnegative or `JACOBIAN_REUSE_AUTO`, got ", max_age, ".", call. = FALSE)
  if (!(max_residual_ratio >= 0)) stop("`max_residual_ratio` must be nonnegative, got ", max_residual_ratio, ".",
                                       call. = FALSE)
  structure(list(max_age = as.integer(max_age), max_residual_ratio = max_residual_ratio), class = "JacobianReuse")
}

.jr_normalize <- function(reuse) {
  if (is.null(reuse)) return(JacobianReuse(-1L, 0.1))
  if (inherits(reuse, "JacobianReuse")) return(reuse)
  if (is.logical(reuse) && length(reuse) == 1L) return(JacobianReuse(if (reuse) 10L else 0L, 0.1))
  stop("`jacobian_reuse` must be `NULL`, a logical, or a `JacobianReuse`.", call. = FALSE)
}

.jr_resolve <- function(policy, n) {
  if (policy$max_age != -1L) return(policy)
  JacobianReuse(if (n >= 16L) 10L else 0L, policy$max_residual_ratio)
}

#' @export
print.JacobianReuse <- function(x, ...) {
  cat("JacobianReuse(max_age = ", x$max_age, ", max_residual_ratio = ", x$max_residual_ratio, ")\n", sep = "")
  invisible(x)
}

#' @rdname jacobian_reuse
#' @export
EisenstatWalkerForcing2 <- function(eta_0 = 0.5, eta_max = 0.9, gamma = 0.9, alpha = 2, safeguard = TRUE,
                                    safeguard_threshold = 0.1)
  structure(list(eta_0 = eta_0, eta_max = eta_max, gamma = gamma, alpha = alpha, safeguard = isTRUE(safeguard),
                 safeguard_threshold = safeguard_threshold), class = "EisenstatWalkerForcing2")

# The forcing cache: init with the first residual, pre_step before the
# descent (iteration 0 gives eta_0) and post_step after an accepted step.
.ew_init <- function(p, fu, internalnorm = function(v) sqrt(sum(v * v))) {
  cc <- new.env(parent = emptyenv())
  cc$p <- p; cc$eta <- p$eta_0; cc$rnorm <- cc$rnorm_prev <- internalnorm(fu); cc$norm <- internalnorm
  cc
}

.ew_pre <- function(cc, fu, iter) {
  if (iter == 0L) {
    cc$eta <- cc$p$eta_0
    cc$rnorm <- cc$rnorm_prev <- cc$norm(fu)
    return(cc$eta)
  }
  prev <- cc$eta
  cc$eta <- cc$p$gamma * (cc$rnorm / cc$rnorm_prev)^cc$p$alpha
  if (cc$p$safeguard) {
    sg <- cc$p$gamma * prev^cc$p$alpha
    if (sg > cc$p$safeguard_threshold && sg > cc$eta) cc$eta <- sg
  }
  cc$eta <- min(max(cc$eta, 0), cc$p$eta_max)
  cc$eta
}

.ew_post <- function(cc, fu) {
  cc$rnorm_prev <- cc$rnorm
  cc$rnorm <- cc$norm(fu)
}

# ---------------------------------------------------------------------------
# Tracing

#' Solver traces
#'
#' Trace levels for the `trace_level` option of `solve()` with the bounded
#' solvers ([TrustRegionReflective()] and the others, and
#' [BoundedTrustRegion()]). `show_trace = TRUE` prints one line per step,
#' `store_trace = TRUE` keeps the entries in `sol$trace`. `TraceMinimal()`
#' records the iteration, the residual norm (infinity norm, or 2-norm for
#' least squares) and the 2-norm of the step;
#' `TraceWithJacobianConditionNumber()` adds the condition number of the
#' Jacobian (-1 when the solver keeps no Jacobian, as in Julia) and
#' `TraceAll()` also stores `u`, `fu`, the step and the Jacobian. A single
#' `freq` sets both frequencies.
#'
#' @param freq Print and store frequency.
#' @param print_frequency,store_frequency Print every so many steps, store
#'   every so many steps.
#' @return A `NonlinearSolveTracing` object.
#' @examples
#' f <- function(u, p) c(u[1] * u[1] - 2, u[2] - 1)
#' prob <- NonlinearProblem(f, u0 = c(1, 0), lb = c(0, 0), ub = c(2, 2))
#' sol <- solve(prob, Dogbox(), store_trace = TRUE, trace_level = TraceMinimal())
#' sol$trace
#' @name nonlinear_trace
NULL

.trace_level <- function(mode, freq, print_frequency, store_frequency)
  structure(list(trace_mode = mode, print_frequency = as.integer(freq %||% print_frequency),
                 store_frequency = as.integer(freq %||% store_frequency)), class = "NonlinearSolveTracing")

#' @rdname nonlinear_trace
#' @export
TraceMinimal <- function(freq = NULL, print_frequency = 1L, store_frequency = 1L)
  .trace_level("minimal", freq, print_frequency, store_frequency)

#' @rdname nonlinear_trace
#' @export
TraceWithJacobianConditionNumber <- function(freq = NULL, print_frequency = 1L, store_frequency = 1L)
  .trace_level("condition_number", freq, print_frequency, store_frequency)

#' @rdname nonlinear_trace
#' @export
TraceAll <- function(freq = NULL, print_frequency = 1L, store_frequency = 1L)
  .trace_level("all", freq, print_frequency, store_frequency)

.trace_entry <- function(ls, iteration, fu, du, J, u, mode) {
  e <- list(iteration = iteration, fnorm = if (ls) sqrt(sum(fu * fu)) else max(abs(fu), 0),
            stepnorm = sqrt(sum(du * du)), condJ = NULL, norm_type = if (ls) "L2" else "Inf")
  if (mode != "minimal") e$condJ <- if (is.matrix(J)) kappa(J, exact = TRUE) else -1
  if (mode == "all") e$storage <- list(u = u, fu = fu, du = du, J = J)
  structure(e, class = "NonlinearSolveTraceEntry")
}

# `jac0` computes the Jacobian at the initial point for the first entry of
# the solvers that keep one.
.trace_init <- function(opts, ls, alg_name, u, fu, du, jac0 = NULL) {
  show <- isTRUE(opts$show_trace); store <- isTRUE(opts$store_trace)
  if (!show && !store) return(NULL)
  level <- opts$trace_level %||% TraceMinimal()
  tr <- new.env(parent = emptyenv())
  tr$show <- show; tr$store <- store; tr$level <- level; tr$ls <- ls; tr$history <- list()
  if (show) cat("\nAlgorithm: ", alg_name, "\n\n", sep = "")
  J <- if (!is.null(jac0) && level$trace_mode != "minimal") jac0()
  e <- .trace_entry(ls, 0L, fu, du, J, u, level$trace_mode)
  if (show) print(e, header = TRUE)
  if (store) tr$history <- list(e)
  tr
}

.trace_update <- function(tr, iter, u, fu, du, alpha = TRUE) {
  if (is.null(tr)) return(invisible(NULL))
  show_now <- tr$show && ((iter - 1L) %% tr$level$print_frequency + 1L) == 1L
  store_now <- tr$store && ((iter - 1L) %% tr$level$store_frequency + 1L) == 1L
  if (show_now || store_now) {
    e <- .trace_entry(tr$ls, iter, fu, if (isTRUE(alpha)) du else rep(0, length(du)), NULL, u, tr$level$trace_mode)
    if (show_now) print(e)
    if (store_now) tr$history[[length(tr$history) + 1L]] <- e
  }
  invisible(NULL)
}

.trace_last <- function(tr, fu) {
  if (is.null(tr)) return(NULL)
  if (tr$show) cat(sprintf("%-8s\t%-20.8e\n%-28s\n", "Final", if (tr$ls) sqrt(sum(fu * fu)) else max(abs(fu), 0),
                           "----------------------"))
  structure(list(history = if (tr$store) tr$history, trace_level = tr$level), class = "NonlinearSolveTrace")
}

#' @export
print.NonlinearSolveTraceEntry <- function(x, header = x$iteration == 0L, ...) {
  nm <- if (x$norm_type == "L2") "f(u) 2-norm" else "f(u) inf-norm"
  if (header) {
    if (is.null(x$condJ)) {
      cat(sprintf("%-8s\t%-20s\t%-20s\n", "----", "-------------", "-----------"))
      cat(sprintf("%-8s\t%-20s\t%-20s\n", "Iter", nm, "Step 2-norm"))
      cat(sprintf("%-8s\t%-20s\t%-20s\n", "----", "-------------", "-----------"))
    } else {
      cat(sprintf("%-8s\t%-20s\t%-20s\t%-20s\n", "----", "-------------", "-----------", "-------"))
      cat(sprintf("%-8s\t%-20s\t%-20s\t%-20s\n", "Iter", nm, "Step 2-norm", "cond(J)"))
      cat(sprintf("%-8s\t%-20s\t%-20s\t%-20s\n", "----", "-------------", "-----------", "-------"))
    }
  }
  if (is.null(x$condJ)) cat(sprintf("%-8d\t%-20.8e\t%-20.8e\n", x$iteration, x$fnorm, x$stepnorm))
  else cat(sprintf("%-8d\t%-20.8e\t%-20.8e\t%-20.8e\n", x$iteration, x$fnorm, x$stepnorm, x$condJ))
  invisible(x)
}

#' @export
print.NonlinearSolveTrace <- function(x, ...) {
  if (is.null(x$history)) { cat("Tracing Disabled\n"); return(invisible(x)) }
  for (e in x$history) print(e)
  invisible(x)
}

# ---------------------------------------------------------------------------
# NonlinearVerbosity

.NL_VERBOSITY_PRESETS <- list(
  None = c(non_enclosing_interval = "Silent", alias_u0_immutable = "Silent", linsolve_failed_noncurrent = "Silent",
           termination_condition = "Silent", threshold_state = "Silent", forcing = "Silent",
           sensitivity_vjp_choice = "Silent", unsupported_postcondition = "Silent"),
  Minimal = c(non_enclosing_interval = "WarnLevel", alias_u0_immutable = "Silent",
              linsolve_failed_noncurrent = "WarnLevel", termination_condition = "Silent", threshold_state = "Silent",
              forcing = "Silent", sensitivity_vjp_choice = "Silent", unsupported_postcondition = "ErrorLevel"),
  Standard = c(non_enclosing_interval = "WarnLevel", alias_u0_immutable = "WarnLevel",
               linsolve_failed_noncurrent = "WarnLevel", termination_condition = "WarnLevel",
               threshold_state = "WarnLevel", forcing = "InfoLevel", sensitivity_vjp_choice = "WarnLevel",
               unsupported_postcondition = "ErrorLevel"),
  Detailed = c(non_enclosing_interval = "WarnLevel", alias_u0_immutable = "WarnLevel",
               linsolve_failed_noncurrent = "WarnLevel", termination_condition = "WarnLevel",
               threshold_state = "WarnLevel", forcing = "InfoLevel", sensitivity_vjp_choice = "WarnLevel",
               unsupported_postcondition = "ErrorLevel"),
  All = c(non_enclosing_interval = "WarnLevel", alias_u0_immutable = "WarnLevel",
          linsolve_failed_noncurrent = "WarnLevel", termination_condition = "WarnLevel", threshold_state = "InfoLevel",
          forcing = "InfoLevel", sensitivity_vjp_choice = "WarnLevel", unsupported_postcondition = "ErrorLevel"))

.NL_VERBOSITY_LINEAR <- c(None = "None", Minimal = "None", Standard = "None", Detailed = "Detailed", All = "Detailed")

.NL_VERBOSITY_GROUPS <- list(
  error_control = c("non_enclosing_interval", "alias_u0_immutable", "linsolve_failed_noncurrent",
                    "termination_condition"),
  numerical = c("threshold_state", "forcing"),
  sensitivity = "sensitivity_vjp_choice")

#' Verbosity of the nonlinear solvers
#'
#' `NonlinearVerbosity()` holds the message levels (`"Silent"`,
#' `"DebugLevel"`, `"InfoLevel"`, `"WarnLevel"`, `"ErrorLevel"`) of the
#' diagnostics of NonlinearSolve, as its SciMLLogging specifier. `preset` is
#' one of `"None"`, `"Minimal"`, `"Standard"` (the default), `"Detailed"` and
#' `"All"` (`TRUE` and `FALSE` stand for `"Standard"` and `"None"`); a group
#' (`error_control`, `numerical`, `sensitivity`) sets all its toggles, and
#' individual toggles in `...` override both. `linear_verbosity` is the
#' verbosity of the linear solves.
#'
#' @param preset A preset name or a logical.
#' @param error_control,numerical,sensitivity Levels of the groups.
#' @param linear_verbosity Verbosity of the linear solver.
#' @param ... Individual toggles, such as `forcing = "Silent"`.
#' @return A `NonlinearVerbosity` object.
#' @examples
#' v <- NonlinearVerbosity(numerical = "Silent", termination_condition = "ErrorLevel")
#' v$toggles
#' NonlinearVerbosity("None")$toggles[["forcing"]]
#' @export
NonlinearVerbosity <- function(preset = "Standard", error_control = NULL, numerical = NULL, sensitivity = NULL,
                               linear_verbosity = NULL, ...) {
  if (isTRUE(preset)) preset <- "Standard"
  if (identical(preset, FALSE)) preset <- "None"
  lv <- .NL_VERBOSITY_PRESETS[[preset]]
  if (is.null(lv)) stop("unknown verbosity preset ", preset, call. = FALSE)
  groups <- list(error_control = error_control, numerical = numerical, sensitivity = sensitivity)
  for (g in names(groups)) if (!is.null(groups[[g]])) lv[.NL_VERBOSITY_GROUPS[[g]]] <- groups[[g]]
  over <- unlist(list(...))
  bad <- setdiff(names(over), names(lv))
  if (length(bad)) stop("unknown verbosity toggle(s): ", paste(bad, collapse = ", "), call. = FALSE)
  lv[names(over)] <- over
  badlv <- setdiff(lv, c("Silent", "DebugLevel", "InfoLevel", "WarnLevel", "ErrorLevel"))
  if (length(badlv)) stop("unknown verbosity level(s): ", paste(badlv, collapse = ", "), call. = FALSE)
  structure(list(toggles = lv, linear_verbosity = linear_verbosity %||% .NL_VERBOSITY_LINEAR[[preset]]),
            class = "NonlinearVerbosity")
}

#' @export
print.NonlinearVerbosity <- function(x, ...) {
  cat("NonlinearVerbosity\n")
  print(x$toggles, quote = FALSE)
  cat("linear_verbosity: ", format(x$linear_verbosity), "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Postcondition specifiers

#' Coordinates of a postcondition corrector
#'
#' `PostconditionSpecifier(corrector, space)` declares in which coordinates
#' the iterate corrector `corrector(u, u_prev, p, cache)` of the
#' `postcondition` solver option is written when the problem has bounds:
#' `PostconditionSpace$Original` (the default) or
#' `PostconditionSpace$Transformed`, the unconstrained variable of solvers
#' that reparameterize the bounds. The R solvers handle bounds natively, in
#' the original coordinates, so both spaces refer to the original variable
#' (as for the native bounded solvers of NonlinearSolve). The result is a
#' function that can be passed as `postcondition`.
#'
#' @param corrector A function `(u, u_prev, p, cache)`.
#' @param space A value of `PostconditionSpace`.
#' @return `PostconditionSpace` is a list of the two spaces;
#'   `PostconditionSpecifier()` returns a corrector function.
#' @examples
#' H <- function(u, u_prev, p, cache) pmax(u, 0)
#' post <- PostconditionSpecifier(H, space = PostconditionSpace$Transformed)
#' attr(post, "space")
#' prob <- NonlinearProblem(function(u, p) u * u - 2, u0 = 1)
#' solve(prob, BoundedTrustRegion(), postcondition = post)$u
#' @name postcondition_specifier
NULL

#' @rdname postcondition_specifier
#' @export
PostconditionSpace <- list(Original = "Original", Transformed = "Transformed")

#' @rdname postcondition_specifier
#' @export
PostconditionSpecifier <- function(corrector, space = PostconditionSpace$Original) {
  space <- match.arg(space, c("Original", "Transformed"))
  force(corrector)
  structure(function(u, u_prev, p, cache = NULL) corrector(u, u_prev, p, cache),
            space = space, class = c("PostconditionSpecifier", "function"))
}
