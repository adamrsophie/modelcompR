# Continuation algorithms for HomotopyProblem, after NonlinearSolveBase
# (v2.55.1: homotopy_sweep.jl, the driver shared by HomotopySweep and
# KantorovichHomotopy; kantorovich_homotopy.jl; arclength.jl; homotopy_polyalg.jl),
# SimpleNonlinearSolve (simple_homotopy_sweep.jl) and NonlinearSolve
# (poly_algs.jl FastShortcutHomotopyPolyalg).
#
# An algorithm here carries a `sweep(at, u0, opts, l0, l1)` function, which
# solve.HomotopyProblem and the SCC solver reach through the HomotopySweep
# class. `at(lambda)` returns the residual F(u) and Jacobian J(u) at lambda;
# a sweep returns list(u, retcode, nsteps, iters, lambda_last), lambda_last
# being the last accepted parameter of a failed continuation.

#' Continuation algorithms
#'
#' Algorithms for [HomotopyProblem()], as in NonlinearSolve:
#'
#' * `SimpleHomotopySweep()` follows the solution from the start of
#'   `lambdaspan` to its end in adaptive steps (natural parameter
#'   continuation), with a secant predictor whose trust is monitored, the
#'   corrector `inner` (default [SimpleNewtonRaphson()]) capped at
#'   `tracking_maxiters` iterations at interior points, step growth after
#'   `expand_threshold` cheap or well predicted steps and bisection after a
#'   failure.
#' * `KantorovichHomotopy()` is the same continuation with the step size of
#'   the Newton Kantorovich controller of Deuflhard: the next increment is
#'   scaled by `gamma * (g(theta_bar) / g(Theta))^(1 / p)`, `g(x) = sqrt(1 +
#'   4x) - 1`, from the first contraction ratio `Theta` of the corrector
#'   (`theta_min` when it is not measured), clamped to `[qmin, qmax]`. With
#'   `strict` a corrector contracting worse than `theta_reject` is rejected.
#'   The contraction is measured for `inner = NewtonRaphson()`; other inner
#'   algorithms use `theta_min`.
#' * `ArcLengthContinuation()` is pseudo arclength continuation: the
#'   unknowns and `lambda` are corrected together on the augmented system
#'   `H(u, lambda) = 0`, `<tau, x - x_k>_theta = ds` (the theta weighted inner
#'   product), so it can pass folds where `lambda` turns back. The predictor
#'   is the secant (or the tangent of the bordered system with `predictor =
#'   "tangent"`); steps turning more than `max_angle` are rejected. When a
#'   step crosses the end of the span the end point is solved at fixed
#'   `lambda`.
#' * `HomotopyPolyAlgorithm(algs)` tries the stages in turn; with
#'   `warm_handoff` a later stage starts from the last point a sweep
#'   reached (backed off by a twentieth of the span), with its step caps
#'   rescaled to the shorter span. `HomotopyPolyAlgorithm()` is
#'   `HomotopySweep()` then `ArcLengthContinuation()`, and
#'   `FastShortcutHomotopyPolyalg()` uses [FastShortcutNonlinearPolyalg()] as
#'   corrector. Inside a polyalgorithm `HomotopySweep()` runs the controller of
#'   NonlinearSolve (secant trust, corrector effort bands and the
#'   `tracking_maxiters` cap) with the settings it has.
#'
#' @param inner The corrector algorithm (`NULL`: the default polyalgorithm).
#' @param nsteps Optional number of equal divisions setting the first step.
#' @param adaptive Adapt the step size (`FALSE` needs `nsteps`).
#' @param initial_step_factor First increment as a fraction of the span.
#' @param min_dlambda,min_ds Smallest increment (default `sqrt(eps)`).
#' @param max_step_factor Largest increment as a fraction of the span.
#' @param expand_factor,expand_threshold,expand_quality Step growth: factor,
#'   accepted steps in a row, largest relative prediction error allowing it.
#' @param predictor `"secant"` or `"constant"` (`"secant"` or `"tangent"`
#'   for `ArcLengthContinuation()`).
#' @param tracking_maxiters Iteration cap of the interior correctors (`NULL`
#'   for none).
#' @param tracking_abstol Optional loose tolerance of the interior correctors.
#' @param maxsteps Largest number of continuation attempts.
#' @param qmin,qmax,theta_min,theta_reject,theta_bar,gamma,strict,predictor_order
#'   Settings of the Kantorovich controller.
#' @param max_angle Largest turn of the path in one step.
#' @param theta Weight of the unknowns in the arclength inner product.
#' @param autodiff,linsolve,concrete_jac,vjp_autodiff,jvp_autodiff Options of
#'   the correctors.
#' @param algs A list of continuation algorithms.
#' @param warm_handoff Start later stages where the sweep stopped.
#' @param store_original Kept for compatibility.
#' @return An algorithm object.
#' @examples
#' y <- Variable("y")
#' sys <- complete(System(list(0 ~ homotopy(atan(y - 3), y)), name = "h"))
#' prob <- HomotopyProblem(sys, list(y ~ 12))
#' solve(prob, SimpleHomotopySweep())$u
#' solve(prob, KantorovichHomotopy(inner = NewtonRaphson()))$u
#' solve(prob, ArcLengthContinuation())$u
#' solve(prob, HomotopyPolyAlgorithm())$retcode
#' @name continuation_algorithms
NULL

.hs_check <- function(name, ok, msg) if (!isTRUE(ok)) stop(name, " ", msg, call. = FALSE)

.hs_alg <- function(name, ...) {
  a <- .nl_alg(name, ...)
  class(a) <- c(name, "HomotopySweep", "NonlinearAlgorithm")
  a$sweep <- function(at, u0, opts, l0, l1) .hs_run(a, at, u0, opts, l0, l1)
  a
}

.hs_validate_common <- function(name, nsteps, adaptive, initial_step_factor, min_dlambda, max_step_factor,
                                tracking_maxiters, tracking_abstol, maxsteps) {
  .hs_check(name, is.null(nsteps) || nsteps >= 1, paste0("`nsteps` must be >= 1, got ", nsteps))
  .hs_check(name, adaptive || !is.null(nsteps),
            "with `adaptive = FALSE` takes fixed-size lambda steps, so an explicit `nsteps` is required.")
  .hs_check(name, initial_step_factor > 0 && initial_step_factor <= 1,
            paste0("`initial_step_factor` must be in (0, 1], got ", initial_step_factor))
  .hs_check(name, is.null(min_dlambda) || min_dlambda > 0, paste0("`min_dlambda` must be positive, got ", min_dlambda))
  .hs_check(name, max_step_factor > 0 && max_step_factor <= 1,
            paste0("`max_step_factor` must be in (0, 1], got ", max_step_factor))
  .hs_check(name, is.null(tracking_maxiters) || tracking_maxiters >= 1,
            paste0("`tracking_maxiters` must be >= 1 (NULL disables the cap), got ", tracking_maxiters))
  .hs_check(name, is.null(tracking_abstol) || tracking_abstol > 0,
            paste0("`tracking_abstol` must be positive (NULL disables the loose tracking tolerance), got ",
                   tracking_abstol))
  .hs_check(name, maxsteps >= 1, paste0("`maxsteps` must be >= 1, got ", maxsteps))
}

#' @rdname continuation_algorithms
#' @export
SimpleHomotopySweep <- function(inner = SimpleNewtonRaphson(), nsteps = NULL, adaptive = TRUE, initial_step_factor = 0.1,
                                min_dlambda = NULL, max_step_factor = 1, expand_factor = 2, expand_threshold = 2L,
                                expand_quality = 0.25, predictor = c("secant", "constant"), tracking_maxiters = 10L,
                                tracking_abstol = NULL, maxsteps = 10000L) {
  predictor <- match.arg(predictor)
  .hs_validate_common("SimpleHomotopySweep", nsteps, adaptive, initial_step_factor, min_dlambda, max_step_factor,
                      tracking_maxiters, tracking_abstol, maxsteps)
  .hs_check("SimpleHomotopySweep", expand_factor >= 1,
            paste0("`expand_factor` must be >= 1 (1 disables expansion), got ", expand_factor))
  .hs_check("SimpleHomotopySweep", expand_threshold >= 1, paste0("`expand_threshold` must be >= 1, got ", expand_threshold))
  .hs_check("SimpleHomotopySweep", expand_quality > 0,
            paste0("`expand_quality` must be positive (Inf disables the gate), got ", expand_quality))
  .hs_alg("SimpleHomotopySweep", inner = inner, nsteps = nsteps, adaptive = isTRUE(adaptive),
          initial_step_factor = initial_step_factor, min_dlambda = min_dlambda, max_step_factor = max_step_factor,
          expand_factor = expand_factor, expand_threshold = as.integer(expand_threshold), expand_quality = expand_quality,
          predictor = predictor, tracking_maxiters = tracking_maxiters, tracking_abstol = tracking_abstol,
          maxsteps = as.integer(maxsteps), controller = "effort", default_inner = "simple")
}

#' @rdname continuation_algorithms
#' @export
KantorovichHomotopy <- function(inner = NULL, nsteps = NULL, initial_step_factor = 0.1, min_dlambda = NULL,
                                max_step_factor = 1, qmin = 1 / 5, qmax = 5, theta_min = 1 / 8, theta_reject = 0.95,
                                theta_bar = 0.5, gamma = 0.95, strict = TRUE, predictor = c("constant", "secant"),
                                predictor_order = NULL, expand_quality = 0.25, tracking_maxiters = 10L,
                                tracking_abstol = NULL, maxsteps = 10000L, store_original = FALSE) {
  predictor <- match.arg(predictor)
  nm <- "KantorovichHomotopy"
  .hs_validate_common(nm, nsteps, TRUE, initial_step_factor, min_dlambda, max_step_factor, tracking_maxiters,
                      tracking_abstol, maxsteps)
  .hs_check(nm, qmin > 0 && qmin < 1, paste0("`qmin` must be in (0, 1), got ", qmin))
  .hs_check(nm, qmax >= 1, paste0("`qmax` must be >= 1, got ", qmax))
  if (!(0 < theta_min && theta_min <= theta_bar && theta_bar < theta_reject && theta_reject < 1))
    stop(sprintf("KantorovichHomotopy requires 0 < theta_min <= theta_bar < theta_reject < 1; got theta_min = %s, theta_bar = %s, theta_reject = %s",
                 theta_min, theta_bar, theta_reject), call. = FALSE)
  .hs_check(nm, gamma > 0 && gamma < 1, paste0("`gamma` must be in (0, 1), got ", gamma))
  predictor_order <- predictor_order %||% (if (predictor == "secant") 2L else 1L)
  .hs_check(nm, predictor_order >= 1, paste0("`predictor_order` must be >= 1, got ", predictor_order))
  .hs_check(nm, expand_quality > 0, paste0("`expand_quality` must be positive, got ", expand_quality))
  .hs_alg(nm, inner = inner, nsteps = nsteps, adaptive = TRUE, initial_step_factor = initial_step_factor,
          min_dlambda = min_dlambda, max_step_factor = max_step_factor, qmin = qmin, qmax = qmax,
          theta_min = theta_min, theta_reject = theta_reject, theta_bar = theta_bar, gamma = gamma,
          strict = isTRUE(strict), predictor = predictor, predictor_order = as.integer(predictor_order),
          expand_quality = expand_quality, tracking_maxiters = tracking_maxiters, tracking_abstol = tracking_abstol,
          maxsteps = as.integer(maxsteps), expand_factor = 1, controller = "kantorovich", default_inner = "poly")
}

# Settings of the NonlinearSolve HomotopySweep controller for an R
# HomotopySweep() object (which has fewer fields).
.hs_from_sweep <- function(alg) {
  if (!is.null(alg$controller)) return(alg)
  alg$nsteps <- alg$nsteps; alg$adaptive <- alg$adaptive %||% TRUE
  alg$expand_quality <- alg$expand_quality %||% 0.25
  alg$predictor <- alg$predictor %||% "secant"
  if (!"tracking_maxiters" %in% names(alg)) alg$tracking_maxiters <- 10L
  alg$controller <- "effort"; alg$default_inner <- "poly"
  alg
}

.hs_effort_growth <- function(nit, budget, expand_factor) {
  if (nit < 0) return(expand_factor)
  if (nit <= max(3L, budget %/% 20L)) return(expand_factor)
  if (nit <= max(6L, budget %/% 4L)) return(1 + (expand_factor - 1) / 2)
  1
}

.hs_effort_shrink <- function(nit, budget) nit >= 0 && nit <= budget && 4 * nit >= 3 * budget

.hs_kfactor <- function(alg, Theta, has) {
  if (has && (!is.finite(Theta) || Theta < 0)) return(alg$qmin)
  eff <- if (has) max(Theta, alg$theta_min) else alg$theta_min
  g <- function(x) sqrt(1 + 4 * x) - 1
  q <- alg$gamma * (g(alg$theta_bar) / g(eff))^(1 / alg$predictor_order)
  min(max(q, alg$qmin), alg$qmax)
}

.hs_l2 <- function(v) sqrt(sum(v * v))

# Newton corrector recording the contraction ratios of the residual norm
# (the measured corrector of KantorovichHomotopy).
.hs_newton <- function(F, J, x, opts) {
  abstol <- opts$abstol %||% .nl_default_tol()
  maxiters <- as.integer(opts$maxiters %||% 1000L)
  fx <- F(x)
  prev <- .hs_l2(fx)
  thetas <- numeric()
  if (max(abs(fx)) <= abstol) return(list(u = x, resid = fx, retcode = "Success", iters = 0L, thetas = thetas))
  for (k in seq_len(maxiters)) {
    Jx <- if (is.null(J)) .jac_central(F, x) else J(x)
    dx <- tryCatch(solve(Jx, fx), error = function(e) .nl_lstsq(Jx, fx))
    x <- x - as.numeric(dx)
    fx <- F(x)
    rn <- .hs_l2(fx)
    thetas <- c(thetas, if (prev == 0) 0 else rn / prev)
    prev <- rn
    if (!all(is.finite(fx))) return(list(u = x, resid = fx, retcode = "Unstable", iters = k, thetas = thetas))
    if (max(abs(fx)) <= abstol) return(list(u = x, resid = fx, retcode = "Success", iters = k, thetas = thetas))
  }
  list(u = x, resid = fx, retcode = "MaxIters", iters = maxiters, thetas = thetas)
}

.hs_inner <- function(alg, u0) {
  alg$inner %||% (if (identical(alg$default_inner, "simple")) SimpleNewtonRaphson() else .nl_default_alg(u0, TRUE))
}

.hs_solve_at <- function(alg, at, lam, guess, opts) {
  fn <- at(lam)
  inner <- .hs_inner(alg, guess)
  if (identical(alg$controller, "kantorovich") && inherits(inner, "NewtonRaphson"))
    return(.hs_newton(fn$F, fn$J, guess, opts))
  res <- .nl_run(inner, fn$F, fn$J, guess, opts, ls = FALSE)
  res
}

# _homotopy_sweep_solve! (and the SimpleHomotopySweep solve)
.hs_run <- function(alg, at, u0, opts, l0, l1) {
  alg <- .hs_from_sweep(alg)
  lam <- as.numeric(l0); lend <- as.numeric(l1)
  span <- lend - lam
  dl <- if (is.null(alg$nsteps)) alg$initial_step_factor * span else span / alg$nsteps
  min_dl <- alg$min_dlambda %||% sqrt(.Machine$double.eps)
  max_dl <- alg$max_step_factor * span
  expand_factor <- if (alg$controller == "kantorovich") 1 else alg$expand_factor
  if (abs(dl) > abs(max_dl)) dl <- max_dl
  # _tracking_budget and _tracking_tolerance
  user_iters <- !is.null(opts$maxiters)
  budget <- if (user_iters) as.integer(opts$maxiters) else if (is.null(alg$tracking_maxiters)) 1000L
            else as.integer(alg$tracking_maxiters)
  cap_active <- !user_iters && !is.null(alg$tracking_maxiters)
  tol_active <- !is.null(alg$tracking_abstol) && is.null(opts$abstol) && is.null(opts$reltol)
  interior_opts <- opts
  if (cap_active) interior_opts$maxiters <- budget
  if (tol_active) interior_opts$abstol <- alg$tracking_abstol
  ok <- function(r) successful_retcode(r$retcode)
  out <- function(u, retcode, lambda_last, nsteps)
    list(u = u, retcode = retcode, nsteps = nsteps, iters = nsteps, lambda_last = lambda_last)
  u <- as.numeric(u0)
  simple <- inherits(alg, "SimpleHomotopySweep")
  # the start point: Julia's sweep cache carries the tracking cap, the
  # simple sweep solves with the user options
  first_opts <- if (simple || tol_active) opts else interior_opts
  last <- .hs_solve_at(alg, at, lam, u, first_opts)
  if (!simple && cap_active && !tol_active && !ok(last)) last <- .hs_solve_at(alg, at, lam, u, opts)
  if (!ok(last)) return(out(u, last$retcode, NULL, 0L))
  u <- last$u
  if (lam == lend) return(out(u, "Success", NULL, 0L))
  u_prev <- u; l_prev <- lam
  streak <- 0L; trust <- 2L; disp_prev <- 0; attempts <- 0L
  repeat {
    attempts <- attempts + 1L
    if (attempts > alg$maxsteps) return(out(u, "MaxIters", lam, attempts - 1L))
    next_l <- if (abs(lend - lam) <= abs(dl)) lend else lam + dl
    if (next_l == lam && next_l != lend) return(out(u, "Stalled", lam, attempts))
    used_secant <- alg$predictor == "secant" && trust >= 2L && l_prev != lam
    guess <- if (used_secant) u + ((next_l - lam) / (lam - l_prev)) * (u - u_prev) else u
    interior <- next_l != lend
    sopts <- if (interior || !simple) interior_opts else opts
    if (!interior && simple) sopts <- opts
    last <- .hs_solve_at(alg, at, next_l, guess, sopts)
    thetas <- last$thetas
    has_T <- length(thetas) > 0L
    first_T <- if (has_T) thetas[1L] else 0
    rejected <- has_T && any(!(thetas <= alg$theta_reject %||% Inf))
    rejected_T <- if (rejected) thetas[which(!(thetas <= alg$theta_reject))[1L]] else 0
    if (!simple && next_l == lend) {
      if (cap_active && !ok(last)) {
        last <- .hs_solve_at(alg, at, next_l, guess, opts)
        has_T <- FALSE; first_T <- 0; rejected <- FALSE; rejected_T <- 0
      } else if (tol_active && ok(last)) {
        last <- .hs_solve_at(alg, at, next_l, last$u, opts)
      }
    }
    strict_rejection <- alg$controller == "kantorovich" && isTRUE(alg$strict) && rejected
    if (ok(last) && !strict_rejection) {
      th <- NULL
      if (l_prev != lam) {
        virtual <- u + ((next_l - lam) / (lam - l_prev)) * (u - u_prev)
        correction <- .hs_l2(last$u - virtual)
        disp <- .hs_l2(last$u - u)
        scale <- max(disp, disp_prev, sqrt(.Machine$double.eps) * (1 + .hs_l2(last$u)))
        th <- correction / scale
        trust <- if (th < 1 / 2) trust + 1L else 0L
        disp_prev <- disp
      } else disp_prev <- .hs_l2(last$u - u)
      u_prev <- u; l_prev <- lam
      u <- last$u; lam <- next_l
      if (lam == lend) break
      nit <- if (simple && inherits(.hs_inner(alg, u), "SimpleNewtonRaphson")) -1L else as.integer(last$iters %||% -1L)
      if (alg$controller == "kantorovich") {
        q <- .hs_kfactor(alg, first_T, has_T)
        cheap <- nit >= 0 && nit <= 2
        if (q > 1 && !(is.null(th) || th <= alg$expand_quality || cheap)) q <- 1
        proposed <- q * dl
        if (abs(proposed) > abs(max_dl)) proposed <- max_dl
        else if (abs(proposed) < min_dl) proposed <- sign(dl) * min_dl
        dl <- proposed; streak <- 0L
      } else if (alg$adaptive) {
        if (.hs_effort_shrink(nit, budget)) {
          if (abs(dl) / 2 >= min_dl) dl <- dl / 2
          streak <- 0L
        } else {
          streak <- streak + 1L
          cheap <- nit >= 0 && nit <= 2
          if (streak >= alg$expand_threshold && (is.null(th) || th <= alg$expand_quality || cheap)) {
            g <- .hs_effort_growth(nit, budget, expand_factor)
            if (g > 1) {
              grown <- g * dl
              dl <- if (abs(grown) > abs(max_dl)) max_dl else grown
              streak <- 0L
            }
          }
        }
      }
    } else {
      attempted <- next_l - lam
      if (alg$controller == "kantorovich") {
        q <- if (ok(last) && rejected && is.finite(rejected_T) && rejected_T >= 0) .hs_kfactor(alg, rejected_T, TRUE)
             else alg$qmin
        proposed <- q * attempted
        can_retry <- abs(proposed) >= min_dl && proposed != 0
      } else if (alg$adaptive) {
        proposed <- dl / 2
        can_retry <- abs(proposed) >= min_dl
      } else {
        proposed <- dl; can_retry <- FALSE
      }
      if (can_retry) {
        dl <- proposed; streak <- 0L; trust <- 0L
      } else {
        return(out(u, if (strict_rejection) "ConvergenceFailure" else last$retcode, lam, attempts))
      }
    }
  }
  out(u, "Success", NULL, attempts)
}

# ---------------------------------------------------------------------------
# ArcLengthContinuation

#' @rdname continuation_algorithms
#' @export
ArcLengthContinuation <- function(inner = NULL, initial_step_factor = 0.1, adaptive = TRUE, min_ds = NULL,
                                  max_step_factor = 1, expand_factor = 2, expand_threshold = 2L, max_angle = pi / 6,
                                  predictor = c("secant", "tangent"), autodiff = NULL, linsolve = NULL,
                                  tracking_maxiters = 10L, maxsteps = 10000L, theta = 0.5, store_original = FALSE) {
  predictor <- match.arg(predictor)
  nm <- "ArcLengthContinuation"
  .hs_check(nm, initial_step_factor > 0 && initial_step_factor <= 1,
            paste0("`initial_step_factor` must be in (0, 1], got ", initial_step_factor))
  .hs_check(nm, is.null(min_ds) || min_ds > 0, paste0("`min_ds` must be positive, got ", min_ds))
  .hs_check(nm, max_step_factor > 0 && max_step_factor <= 1,
            paste0("`max_step_factor` must be in (0, 1], got ", max_step_factor))
  .hs_check(nm, expand_factor >= 1, paste0("`expand_factor` must be >= 1, got ", expand_factor))
  .hs_check(nm, expand_threshold >= 1, paste0("`expand_threshold` must be >= 1, got ", expand_threshold))
  .hs_check(nm, max_angle > 0 && max_angle <= pi, paste0("`max_angle` must be in (0, pi], got ", max_angle))
  .hs_check(nm, is.null(tracking_maxiters) || tracking_maxiters >= 1,
            paste0("`tracking_maxiters` must be >= 1, got ", tracking_maxiters))
  .hs_check(nm, maxsteps >= 1, paste0("`maxsteps` must be >= 1, got ", maxsteps))
  .hs_check(nm, theta > 0 && theta < 1, paste0("`theta` must be in (0, 1), got ", theta))
  a <- .nl_alg(nm, inner = inner, initial_step_factor = initial_step_factor, adaptive = isTRUE(adaptive),
               min_ds = min_ds, max_step_factor = max_step_factor, expand_factor = expand_factor,
               expand_threshold = as.integer(expand_threshold), max_angle = max_angle, predictor = predictor,
               autodiff = autodiff, linsolve = linsolve, tracking_maxiters = tracking_maxiters,
               maxsteps = as.integer(maxsteps), theta = theta)
  class(a) <- c(nm, "HomotopySweep", "NonlinearAlgorithm")
  a$sweep <- function(at, u0, opts, l0, l1) .alc_run(a, at, u0, opts, l0, l1)
  a
}

.alc_run <- function(alg, at, u0, opts, l0, l1) {
  inner <- alg$inner %||% .nl_default_alg(u0, TRUE)
  lam <- as.numeric(l0); lend <- as.numeric(l1)
  span <- lend - lam
  out <- function(u, retcode, nsteps) list(u = u, retcode = retcode, nsteps = nsteps, iters = nsteps, lambda_last = NULL)
  fixed <- function(u, l) { fn <- at(l); .nl_run(inner, fn$F, fn$J, u, opts, ls = FALSE) }
  start <- fixed(as.numeric(u0), lam)
  if (!successful_retcode(start$retcode)) return(out(as.numeric(u0), start$retcode, 0L))
  u <- start$u
  if (span == 0) return(out(u, "Success", 0L))
  n <- length(u)
  xcur <- c(u, lam); xprev <- xcur; have_prev <- FALSE
  min_ds <- alg$min_ds %||% sqrt(.Machine$double.eps)
  max_ds <- alg$max_step_factor * abs(span)
  ds <- min(alg$initial_step_factor * abs(span), max_ds)
  cos_reject <- cos(alg$max_angle); cos_grow <- cos(alg$max_angle / 3)
  streak <- 0L
  user_iters <- !is.null(opts$maxiters)
  budget <- if (user_iters) as.integer(opts$maxiters) else if (is.null(alg$tracking_maxiters)) 1000L
            else as.integer(alg$tracking_maxiters)
  copts <- opts
  if (!user_iters && !is.null(alg$tracking_maxiters)) copts$maxiters <- budget
  wu <- alg$theta / n; wl <- 1 - alg$theta
  tdot <- function(x, y) wu * sum(x[seq_len(n)] * y[seq_len(n)]) + wl * x[n + 1L] * y[n + 1L]
  tnorm <- function(x) sqrt(tdot(x, x))
  tau_seed <- c(rep(0, n), sign(span) / sqrt(wl))
  tau <- tau_seed
  # the path Jacobian [dH/du, dH/dlambda] (dH/dlambda by central differences)
  pjac <- function(x) {
    l <- x[n + 1L]; uu <- x[seq_len(n)]
    fn <- at(l)
    Ju <- if (is.null(fn$J)) .jac_central(fn$F, uu) else fn$J(uu)
    h <- (.Machine$double.eps)^(1 / 3) * max(1, abs(l))
    dl <- (at(l + h)$F(uu) - at(l - h)$F(uu)) / (2 * h)
    cbind(matrix(Ju, n, n), dl)
  }
  tangent <- function(x, tprev) {
    Jp <- pjac(x)
    B <- rbind(Jp, c(wu * tprev[seq_len(n)], wl * tprev[n + 1L]))
    t <- tryCatch(as.numeric(solve(B, c(rep(0, n), 1))), error = function(e) NULL)
    solved <- FALSE
    if (!is.null(t) && all(is.finite(t))) {
      r <- as.numeric(Jp %*% t)
      resid <- sqrt(sum(r * r) + (tdot(t, tprev) - 1)^2)
      nt <- tnorm(t)
      if (is.finite(nt) && nt > 0 && resid <= sqrt(.Machine$double.eps) * (1 + nt)) { t <- t / nt; solved <- TRUE }
    }
    if (!solved) {
      sv <- svd(Jp, nv = n + 1L)
      t <- sv$v[, n + 1L]
      t <- t / tnorm(t)
    }
    if (tdot(t, tprev) < 0) t <- -t
    t
  }
  steps <- 0L
  for (k in seq_len(alg$maxsteps)) {
    steps <- k
    if (alg$predictor == "tangent") tau <- tangent(xcur, tau)
    else if (have_prev) {
      d <- xcur - xprev
      dn <- tnorm(d)
      tau <- if (dn > 0) d / dn else tau_seed
    } else tau <- tau_seed
    xpred <- xcur + ds * tau
    x0 <- xcur; tau_k <- tau; ds_k <- ds
    Faug <- function(x) c(at(x[n + 1L])$F(x[seq_len(n)]), tdot(tau_k, x - x0) - ds_k)
    Jaug <- function(x) rbind(pjac(x), c(wu * tau_k[seq_len(n)], wl * tau_k[n + 1L]))
    last <- .nl_run(inner, Faug, Jaug, xpred, copts, ls = FALSE)
    if (successful_retcode(last$retcode)) {
      xnew <- last$u
      chord <- xnew - xcur
      nchord <- tnorm(chord)
      trusted <- alg$predictor == "tangent" || have_prev
      cosang <- if (trusted && nchord > 0) min(max(tdot(tau, chord) / nchord, -1), 1) else 1
      if (trusted && cosang < cos_reject && alg$adaptive && ds / 2 >= min_ds) {
        ds <- ds / 2; streak <- 0L
        next
      }
      lold <- lam; lnew <- xnew[n + 1L]
      xprev <- xcur; xcur <- xnew
      u <- xcur[seq_len(n)]; lam <- lnew
      have_prev <- TRUE
      if ((lold - lend) * (lnew - lend) <= 0) {
        den <- lnew - lold
        frac <- if (den == 0) 1 else (lend - lold) / den
        uland <- xprev[seq_len(n)] + frac * (xcur[seq_len(n)] - xprev[seq_len(n)])
        fin <- fixed(uland, lend)
        if (successful_retcode(fin$retcode)) return(out(fin$u, "Success", k))
      }
      if (alg$adaptive) {
        nit <- as.integer(last$iters %||% -1L)
        if (.hs_effort_shrink(nit, budget)) {
          if (ds / 2 >= min_ds) ds <- ds / 2
          streak <- 0L
        } else {
          streak <- streak + 1L
          if (streak >= alg$expand_threshold && cosang >= cos_grow) {
            g <- .hs_effort_growth(nit, budget, alg$expand_factor)
            if (g > 1) { ds <- min(g * ds, max_ds); streak <- 0L }
          }
        }
      }
    } else if (alg$adaptive && ds / 2 >= min_ds) {
      ds <- ds / 2; streak <- 0L
    } else return(out(u, last$retcode, k))
  }
  out(u, "MaxIters", steps)
}

# ---------------------------------------------------------------------------
# HomotopyPolyAlgorithm

#' @rdname continuation_algorithms
#' @export
HomotopyPolyAlgorithm <- function(algs = NULL, inner = NULL, warm_handoff = TRUE, store_original = FALSE) {
  if (is.null(algs)) {
    algs <- list(HomotopySweep(inner = inner), ArcLengthContinuation(inner = inner))
  }
  if (inherits(algs, "NonlinearAlgorithm")) algs <- list(algs)
  if (!length(algs)) stop("HomotopyPolyAlgorithm requires at least one algorithm", call. = FALSE)
  a <- .nl_alg("HomotopyPolyAlgorithm", stages = algs, warm_handoff = isTRUE(warm_handoff))
  class(a) <- c("HomotopyPolyAlgorithm", "HomotopySweep", "NonlinearAlgorithm")
  a$sweep <- function(at, u0, opts, l0, l1) .hpa_run(a, at, u0, opts, l0, l1)
  a
}

.hpa_rescale <- function(stage, scale) {
  if (inherits(stage, "ArcLengthContinuation")) {
    stage$max_step_factor <- min(stage$max_step_factor * scale, 1)
    return(ArcLengthContinuation(stage$inner, stage$initial_step_factor, stage$adaptive, stage$min_ds,
                                 stage$max_step_factor, stage$expand_factor, stage$expand_threshold, stage$max_angle,
                                 stage$predictor, stage$autodiff, stage$linsolve, stage$tracking_maxiters,
                                 stage$maxsteps, stage$theta))
  }
  if (inherits(stage, c("KantorovichHomotopy", "SimpleHomotopySweep")) || class(stage)[1L] == "HomotopySweep") {
    if (!is.null(stage$nsteps)) stage$nsteps <- max(1L, as.integer(ceiling(stage$nsteps / scale)))
    stage$max_step_factor <- min(stage$max_step_factor * scale, 1)
    if (!is.null(stage$sweep)) {
      s <- stage
      stage$sweep <- function(at, u0, opts, l0, l1) .hs_run(s, at, u0, opts, l0, l1)
    }
  }
  stage
}

.hpa_stage <- function(stage, at, u0, opts, l0, l1) {
  if (!is.null(stage$sweep)) return(stage$sweep(at, u0, opts, l0, l1))
  if (inherits(stage, "HomotopySweep")) return(.hs_run(stage, at, u0, opts, l0, l1))
  # a plain nonlinear algorithm solves at the end of the span
  fn <- at(l1)
  r <- .nl_run(stage, fn$F, fn$J, u0, opts, ls = FALSE)
  list(u = r$u, retcode = r$retcode, nsteps = r$iters, iters = r$iters, lambda_last = NULL)
}

.hpa_run <- function(alg, at, u0, opts, l0, l1) {
  handoff <- NULL
  nst <- length(alg$stages)
  for (i in seq_len(nst)) {
    stage <- alg$stages[[i]]
    if (alg$warm_handoff && !is.null(handoff)) {
      scale <- abs(l1 - l0) / abs(l1 - handoff$lambda)
      hs <- .hpa_stage(.hpa_rescale(stage, scale), at, handoff$u, opts, handoff$lambda, l1)
      if (successful_retcode(hs$retcode)) return(hs)
    }
    sweepish <- inherits(stage, c("KantorovichHomotopy", "SimpleHomotopySweep")) || class(stage)[1L] == "HomotopySweep"
    sol <- .hpa_stage(stage, at, u0, opts, l0, l1)
    if (sweepish && !is.null(sol$lambda_last)) {
      backoff <- (l1 - l0) / 20
      if (abs(sol$lambda_last - l0) > abs(backoff)) handoff <- list(u = sol$u, lambda = sol$lambda_last - backoff)
    }
    if (successful_retcode(sol$retcode) || i == nst) return(sol)
  }
}

#' @rdname continuation_algorithms
#' @export
FastShortcutHomotopyPolyalg <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, vjp_autodiff = NULL,
                                        jvp_autodiff = NULL, warm_handoff = TRUE, store_original = FALSE) {
  inner <- FastShortcutNonlinearPolyalg()
  HomotopyPolyAlgorithm(inner = inner, warm_handoff = warm_handoff)
}
