# Continuous adjoint sensitivity analysis of ODE solutions, after
# SciMLSensitivity.jl v7.119 src/sensitivity_interface.jl
# (adjoint_sensitivities, _adjoint_sensitivities), src/interpolating_adjoint.jl
# (ODEInterpolatingAdjointSensitivityFunction, ODEAdjointProblem, checkpoint
# re-solves), src/backsolve_adjoint.jl (ODEBacksolveSensitivityFunction,
# backsolve_checkpoint_callbacks), src/quadrature_adjoint.jl
# (AdjointSensitivityIntegrand), src/gauss_adjoint.jl (GaussIntegrand and
# the per step integration), src/adjoint_common.jl (generate_callbacks,
# ReverseLossCallback with the semi-explicit DAE jumps, separate_nonunique)
# and src/callback_tracking.jl (track_callbacks, setup_reverse_callbacks).
#
# Conventions as in Julia: with lambda(T) = 0 and mu(T) = 0 the adjoint
#   M' dlambda/dt = -J' lambda - g_u,   dmu/dt = -f_p' lambda - g_p
# is integrated from T back to t0; a discrete cost sum_i g_i(u(t_i)) adds
# the jumps lambda += g_u(t_i). Then dG/du0 = M' lambda(t0) and
# dG/dp = mu(t0). InterpolatingAdjoint carries (lambda, mu), BacksolveAdjoint
# (lambda, mu, u), QuadratureAdjoint and the Gauss methods lambda alone.
#
# Events recorded during the forward solve (see .sens_solve_hook) are
# pulled back at their times with the transposed event linearization of
# .sens_event_jacobian(); the forward state is read from the segment
# between events that the backward solve is in.

#' Adjoint sensitivities of an ODE solution
#'
#' `adjoint_sensitivities(sol, alg, ...)` returns the gradient of a cost
#' functional of the solution with respect to the initial values and the
#' parameters, by solving the adjoint equation backward in time:
#' `G = int g(u, p, t) dt` (continuous cost `g`, or its gradients
#' `dgdu_continuous` and `dgdp_continuous`) plus `sum_i g_i(u(t_i))` at the
#' times `t` (discrete cost, given by its gradients `dgdu_discrete(u, p, t,
#' i)` and optionally `dgdp_discrete(u, p, t, i)`). Gradients that are not
#' given are taken by finite differences of `g`.
#'
#' The forward solution must be dense (solve without `saveat`, with `dense =
#' TRUE`, or with `sensealg = ...`) unless the method uses checkpointing.
#' Problems with events must be solved with `solve(prob, alg, sensealg =
#' ...)`, which records the events so they are differentiated: the adjoint
#' jumps at each event by the transposed linearization of the event,
#' including the dependence of the event time on the state and parameters,
#' and changes of parameters by the event. Mass matrix (semi-explicit index
#' 1) DAEs are supported by the interpolating, quadrature and Gauss methods.
#' For nonlinear and steady state solutions, see [SteadyStateAdjoint()].
#'
#' `ODEAdjointProblem()` returns the backward problem without solving it;
#' its element `finish` turns the backward solution into the gradients.
#'
#' The result follows Julia: `dp` is the gradient with respect to the
#' problem parameters with the initial values held fixed. When the initial
#' values depend on the parameters through initialization,
#' [solve_gradient()] adds that dependence.
#'
#' @param sol A solution of an `ODEProblem` (or a nonlinear or steady state
#'   solution, see [SteadyStateAdjoint()]).
#' @param alg ODE algorithm of the backward solve; the forward algorithm when
#'   `NULL`.
#' @param t Times of the discrete cost.
#' @param dgdu_discrete,dgdp_discrete Functions `(u, p, t, i)` returning the
#'   gradient of the `i`th discrete cost term.
#' @param dgdu_continuous,dgdp_continuous Functions `(u, p, t)` returning
#'   the gradients of the continuous cost.
#' @param g The continuous cost `g(u, p, t)`.
#' @param abstol,reltol Tolerances of the backward solve.
#' @param checkpoints Saved times used by checkpointing methods (default
#'   `sol$t`).
#' @param callback Further callbacks for the backward solve.
#' @param no_start Leave out the discrete cost at the initial time.
#' @param sensealg An adjoint method, [InterpolatingAdjoint()] by default
#'   (or the `sensealg` the solution was computed with).
#' @param dgdu,dgdp For steady state and nonlinear solutions: functions
#'   `(u, p)` returning the gradients of the cost `g(u, p)` (or give `g`).
#' @param ... Options for the backward [solve()].
#' @return A list with `du0` and `dp`.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
#' sol <- solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-12)
#' # G = u(1): dG/du0 = exp(-2), dG/dp = -exp(-2)
#' res <- adjoint_sensitivities(sol, Tsit5(), t = 1, dgdu_discrete = function(u, p, t, i) 1,
#'                              abstol = 1e-10, reltol = 1e-10)
#' c(res$du0, res$dp) / exp(-2)
#' @export
adjoint_sensitivities <- function(sol, alg = NULL, t = NULL, dgdu_discrete = NULL, dgdp_discrete = NULL,
                                  dgdu_continuous = NULL, dgdp_continuous = NULL, g = NULL,
                                  abstol = 1e-6, reltol = 1e-3, checkpoints = NULL, callback = NULL,
                                  no_start = FALSE, sensealg = NULL, dgdu = NULL, dgdp = NULL, ...) {
  env <- .sol_env(sol)
  prob <- env$prob
  sensealg <- sensealg %||% if (.is_sensealg(env$sensealg, "AdjointSensitivityAlgorithm")) env$sensealg
  if (inherits(prob, c("AbstractNonlinearProblem", "SteadyStateProblem")) ||
      inherits(sensealg, "SteadyStateAdjoint"))
    return(.sens_steady_adjoint(sol, sensealg %||% SteadyStateAdjoint(), dgdu = dgdu, dgdp = dgdp, g = g))
  sensealg <- sensealg %||% InterpolatingAdjoint()
  if (!inherits(prob, "ODEProblem"))
    stop("Continuous adjoint sensitivities are only supported for ODE problems (and mass matrix DAEs).",
         call. = FALSE)
  if (!.is_sensealg(sensealg, "AdjointSensitivityAlgorithm"))
    stop("sensealg must be an adjoint method such as InterpolatingAdjoint()", call. = FALSE)
  if (inherits(sensealg, "QuadratureAdjoint") && missing(abstol)) abstol <- sensealg$abstol
  if (inherits(sensealg, "QuadratureAdjoint") && missing(reltol)) reltol <- sensealg$reltol
  ap <- ODEAdjointProblem(sol, sensealg, alg, t, dgdu_discrete, dgdp_discrete, dgdu_continuous,
                          dgdp_continuous, g, checkpoints = checkpoints, callback = callback,
                          no_start = no_start, abstol = abstol, reltol = reltol)
  adj <- solve(ap, ap$alg, abstol = abstol, reltol = reltol, ...)
  ap$finish(adj)
}

#' @rdname adjoint_sensitivities
#' @export
ODEAdjointProblem <- function(sol, sensealg = InterpolatingAdjoint(), alg = NULL, t = NULL,
                              dgdu_discrete = NULL, dgdp_discrete = NULL, dgdu_continuous = NULL,
                              dgdp_continuous = NULL, g = NULL, checkpoints = NULL, callback = NULL,
                              no_start = FALSE, abstol = 1e-6, reltol = 1e-3) {
  env <- .sol_env(sol)
  prob <- env$prob
  if (is.null(dgdu_discrete) && is.null(dgdp_discrete) && is.null(dgdu_continuous) && is.null(g) &&
      is.null(dgdp_continuous))
    stop("Either `dgdu_discrete`, `dgdu_continuous`, or `g` must be specified.", call. = FALSE)
  if (!is.null(t) && is.null(dgdu_discrete) && is.null(dgdp_discrete))
    stop("It looks like you're using the direct `adjoint_sensitivities` interface with a discrete cost ",
         "function but no specified `dgdu_discrete` or `dgdp_discrete`.", call. = FALSE)
  method <- class(sensealg)[1L]
  quad <- method %in% c("QuadratureAdjoint", "GaussAdjoint", "GaussKronrodAdjoint")
  backsolve <- method == "BacksolveAdjoint"
  ckpt <- isTRUE(sensealg$checkpointing) && !backsolve
  alg <- alg %||% env$alg %||% Tsit5()
  fns <- .sens_fns(prob, sensealg)
  n <- fns$n; np <- fns$np
  ms <- .sens_mass_structure(fns$mass, n)
  events <- .sens_events_of(sol)
  if (length(events) && !is.null(fns$mass))
    stop("events in mass matrix problems are not supported by the sensitivity methods", call. = FALSE)
  if (prob$tspan[2L] < prob$tspan[1L] && (length(events) || ckpt || method == "QuadratureAdjoint"))
    stop("a decreasing time span is supported without events, checkpointing and QuadratureAdjoint",
         call. = FALSE)
  if (length(events) && is.null(g) && (!is.null(dgdu_continuous) || !is.null(dgdp_continuous)))
    stop("a continuous cost with events needs `g`: the cost jumps at the event times", call. = FALSE)
  if (length(events) && ckpt)
    stop("checkpointing together with events is not supported; solve densely instead", call. = FALSE)
  ev_t <- vapply(events, `[[`, 0, "t")
  t0 <- prob$tspan[1L]
  tend <- if (identical(env$retcode, "Terminated")) env$t[length(env$t)] else prob$tspan[2L]
  tol <- function(x) 1e-12 * max(1, abs(x))
  # strictly between t0 and tend, for either direction of the time span
  dir <- if (tend < t0) -1 else 1
  inside <- function(x) dir * (x - t0) > tol(t0) & dir * (tend - x) > tol(tend)

  # forward trajectory
  tr <- if (backsolve) NULL else if (ckpt) {
    cps <- sort(unique(checkpoints %||% env$t))
    if (abs(cps[1L] - t0) > tol(t0)) cps <- c(t0, cps)
    if (abs(cps[length(cps)] - tend) > tol(tend)) cps <- c(cps, tend)
    states <- lapply(cps, function(tc) {
      k <- which(abs(env$t - tc) <= tol(tc))
      if (!length(k)) stop("checkpoint ", tc, " is not a saved time of the solution", call. = FALSE)
      env$U[k[length(k)], ]
    })
    .sens_trajectory(sol, checkpoint = list(times = cps, states = states, prob = prob, alg = env$alg,
                                            opts = list(abstol = abstol, reltol = reltol)))
  } else .sens_trajectory(sol, ev_t)

  # costs
  cont <- !is.null(g) || !is.null(dgdu_continuous) || !is.null(dgdp_continuous)
  gu_c <- dgdu_continuous %||% if (!is.null(g)) function(u, p, tt) .sens_fd_grad(function(z) g(z, p, tt), u)
  gp_c <- dgdp_continuous %||% if (!is.null(g) && is.null(dgdu_continuous) && np)
    function(u, p, tt) .sens_fd_grad(function(q) {
      qq <- p
      qq[] <- q
      g(u, qq, tt)
    }, p)
  ts <- if (is.null(t)) numeric() else as.numeric(t)
  ord <- order(ts)
  # samples at an event time move with it, unless marked as fixed times
  # (requested saveat times, see solve_gradient())
  fixed <- attr(t, "fixed") %||% rep(FALSE, length(ts))
  cost <- list(t = ts[ord], i = ord, fixed = fixed[ord])
  pfinal <- env$p_final %||% prob$p

  # adjoint state layout
  il <- seq_len(n)
  imu <- if (quad) integer() else n + seq_len(np)
  iy <- if (backsolve) n + np + seq_len(n) else integer()
  acc <- new.env(parent = emptyenv())
  acc$mu <- numeric(np)       # parameter gradient of the Gauss methods
  acc$records <- list()       # event records of QuadratureAdjoint
  acc$seg <- length(events)
  state_y <- function(z, tt) if (backsolve) z[iy] else tr$at(tt)

  rhs <- function(z, p, tt) {
    l <- z[il]
    y <- state_y(z, tt)
    v <- fns$vjp(l, y, p, tt)
    dl <- -v$u
    dmu <- -v$p
    if (cont) {
      if (!is.null(gu_c)) dl <- dl - as.numeric(gu_c(y, p, tt))
      if (!is.null(gp_c)) dmu <- dmu - as.numeric(gp_c(y, p, tt))
    }
    c(dl, if (!quad) dmu, if (backsolve) fns$f(y, p, tt))
  }
  jac <- if (backsolve) NULL else function(z, p, tt) {
    y <- tr$at(tt)
    J <- fns$jac(y, p, tt)
    if (quad) return(-t(J))
    rbind(cbind(-t(J), matrix(0, n, np)), cbind(-t(fns$pjac(y, p, tt)), matrix(0, np, np)))
  }
  Madj <- if (is.null(fns$mass)) NULL else
    .sens_blockdiag(c(list(t(fns$mass)), if (!quad && np) list(diag(1, np)), if (backsolve) list(fns$mass)))

  # one discrete cost sample: the jump of lambda (with the semi-explicit DAE
  # correction of ReverseLossCallback) and of the parameter gradient
  sample_jump <- function(z, p, tt, i, y, moving = NULL) {
    gu <- if (!is.null(dgdu_discrete)) as.numeric(dgdu_discrete(y, p, tt, i)) else numeric(n)
    gp <- if (!is.null(dgdp_discrete) && np) as.numeric(dgdp_discrete(y, p, tt, i)) else numeric(np)
    if (ms$dae) {
      J <- fns$jac(y, p, tt)
      ha <- J[ms$alg_eq, ms$alg, drop = FALSE]
      hd <- J[ms$alg_eq, ms$diff, drop = FALSE]
      eta <- -solve(t(ha), gu[ms$alg])
      z[ms$diff] <- z[ms$diff] + solve(t(ms$Mdd), gu[ms$diff] + as.numeric(crossprod(hd, eta)))
      if (np) gp <- gp + as.numeric(crossprod(fns$pjac(y, p, tt)[ms$alg_eq, , drop = FALSE], eta))
    } else {
      dl <- if (is.null(fns$mass)) gu else solve(t(fns$mass), gu)
      z[il] <- z[il] + dl
    }
    add_mu(z, gp)
  }
  add_mu <- function(z, v) {
    if (!np || !any(v != 0)) return(z)
    if (quad) {
      if (method == "QuadratureAdjoint") acc$records[[length(acc$records) + 1L]] <- list(seg = acc$seg, add = v)
      else acc$mu <- acc$mu + v
      return(z)
    }
    z[imu] <- z[imu] + v
    z
  }
  get_mu <- function(z) if (!np) numeric() else if (quad) acc$mu else z[imu]

  # all discrete work at one time point of the backward solve: samples
  # after the events at this time (right limits), event pullbacks (latest
  # first), samples before them (left limits) and checkpoint resets
  ck_t <- if (backsolve && isTRUE(sensealg$checkpointing)) {
    tt_ <- env$t[!duplicated(env$t) & !(env$t %in% env$t[duplicated(env$t)])]
    tt_[inside(tt_) & !vapply(tt_, function(x) any(abs(ev_t - x) <= tol(x)), TRUE)]
  } else numeric()
  point <- function(tt, z, p, skip_cost = FALSE) {
    eidx <- which(abs(ev_t - tt) <= tol(tt))
    cidx <- if (skip_cost) integer() else which(abs(cost$t - tt) <= tol(tt))
    if (backsolve && length(ck_t) && any(abs(ck_t - tt) <= tol(tt))) {
      k <- which(abs(env$t - tt) <= tol(tt))
      z[iy] <- env$U[k[1L], ]
    }
    if (!length(eidx)) {
      for (k in rev(cidx)) z <- sample_jump(z, p, tt, cost$i[k], state_y(z, tt))
      return(list(z = z, p = p))
    }
    Ds <- lapply(events[eidx], .sens_event_jacobian, fns = fns, prob = prob)
    # two samples at an event time are its left and right limits; one is
    # the value after the event, as the forward solver saves it
    nleft <- if (length(cidx) >= 2L) 1L else 0L
    left <- cidx[seq_len(nleft)]
    right <- setdiff(cidx, left)
    last_ev <- events[[eidx[length(eidx)]]]
    lr <- numeric(n)
    for (k in rev(right)) {
      y <- if (backsolve) z[iy] else as.numeric(last_ev$u_right)
      zb <- z[il]
      z <- sample_jump(z, p, tt, cost$i[k], y)
      if (!cost$fixed[k]) lr <- lr + z[il] - zb
    }
    for (j in rev(seq_along(eidx))) {
      ev <- events[[eidx[j]]]
      D <- Ds[[j]]
      tx <- attr(D, "tau_x")
      w <- as.numeric(crossprod(D, c(z[il], get_mu(z), if (quad && !np) numeric())))
      if (!is.null(tx)) {
        # samples saved at the event time move with it (Lu_right of
        # callback_tracking.jl), and a running cost jumps there
        fp <- fns$f(as.numeric(ev$u_right), .sens_pvec(p, ev$p_right), tt)
        w <- w + tx * sum(fp * lr)
        lr <- numeric(n)
        # (after a terminating event nothing is integrated)
        if (cont && !is.null(g))
          w <- w + tx * (g(ev$u, .sens_pvec(p, ev$p), tt) -
                           if (isTRUE(ev$terminate)) 0 else g(as.numeric(ev$u_right), .sens_pvec(p, ev$p_right), tt))
      }
      if (np && method == "QuadratureAdjoint") {
        if (any(abs(D[n + seq_len(np), il]) > 1e-10))
          stop("QuadratureAdjoint does not support events whose parameter changes depend on the state; ",
               "use InterpolatingAdjoint() or GaussAdjoint()", call. = FALSE)
        acc$records[[length(acc$records) + 1L]] <-
          list(seg = acc$seg, A = t(D[n + seq_len(np), n + seq_len(np), drop = FALSE]), add = w[n + seq_len(np)])
      }
      z[il] <- w[il]
      if (np) {
        if (!quad) z[imu] <- w[n + seq_len(np)] else if (method != "QuadratureAdjoint") acc$mu <- w[n + seq_len(np)]
      }
      if (backsolve) z[iy] <- ev$u
      p <- .sens_pvec(p, ev$p)
      acc$seg <- acc$seg - 1L
      if (!is.null(tr)) tr$seg <- acc$seg
    }
    tx1 <- attr(Ds[[1L]], "tau_x")
    for (k in rev(left)) {
      y <- if (backsolve) z[iy] else events[[eidx[1L]]]$u
      zb <- z[il]
      z <- sample_jump(z, p, tt, cost$i[k], y)
      if (!is.null(tx1)) {
        s <- sum(fns$f(y, p, tt) * (z[il] - zb))
        z[il] <- z[il] + tx1[il] * s
        if (np) z <- add_mu(z, tx1[n + seq_len(np)] * s)
      }
    }
    list(z = z, p = p)
  }

  # terminal conditions
  z0 <- numeric(n + length(imu) + length(iy))
  if (backsolve) {
    # the state at the final time, which need not be a saved time
    z0[iy] <- if (abs(env$t[length(env$t)] - tend) <= tol(tend)) env$U[nrow(env$U), ] else
      if (length(env$steps)) .interp_states(env, tend)[1L, ] else
        stop("BacksolveAdjoint needs the state at the final time: save it or solve densely", call. = FALSE)
  }
  pT <- pfinal
  r <- point(tend, z0, pT)
  z0 <- r$z
  pT <- r$p
  # interior time points
  inner <- sort(unique(c(cost$t, ev_t, ck_t)))
  inner <- inner[inside(inner)]
  cb_rev <- if (length(inner)) {
    cb <- DiscreteCallback(function(u, tt, integ) .near_any(inner, tt), function(integ) {
      r <- point(integ$t, integ$u, integ$p)
      integ$u <- r$z
      integ$p <- r$p
    }, save_positions = c(FALSE, FALSE))
    cb$tstops <- inner
    cb
  }
  gauss_cb <- if (quad && method != "QuadratureAdjoint" && np) {
    # integrate mu over each accepted step; the condition never fires, so the
    # driver keeps its step size history
    DiscreteCallback(function(u, tt, integ) {
      st <- integ$last_step
      a <- st$t0; b <- st$tend %||% (st$t0 + st$dt)
      integrand <- function(s) mu_integrand(.step_value(st, s)[il], integ$p, s)
      v <- if (method == "GaussAdjoint") .sens_gauss7(integrand, b, a) else
        .sens_quadgk(integrand, b, a, abstol, reltol)
      acc$mu <- acc$mu + v
      FALSE
    }, function(integ) NULL)
  }
  mu_integrand <- function(l, p, s) {
    y <- tr$at(s)
    v <- fns$vjp(l, y, p, s)$p
    if (cont && !is.null(gp_c)) v <- v + as.numeric(gp_c(y, p, s))
    v
  }
  fn <- ODEFunction(rhs, mass_matrix = Madj, jac_fn = jac)
  aprob <- .new_problem("ODEProblem", f = fn, u0 = z0, tspan = c(tend, t0), p = pT, sys = NULL,
                        kwargs = list(callback = CallbackSet(gauss_cb, cb_rev, callback),
                                      dense = method == "QuadratureAdjoint",
                                      save_everystep = method == "QuadratureAdjoint"))
  aprob$alg <- alg
  aprob$finish <- function(adj) {
    aenv <- .sol_env(adj)
    zf <- aenv$U[nrow(aenv$U), ]
    if (!no_start) zf <- point(t0, zf, aenv$p_final)$z
    lam <- zf[il]
    du0 <- if (is.null(fns$mass)) lam else as.numeric(crossprod(fns$mass, lam))
    dp <- if (!np) numeric() else if (!quad) zf[imu] else if (method != "QuadratureAdjoint") acc$mu else
      .sens_quadrature_dp(adj, events, acc, tr, mu_integrand, cost, prob, pfinal, t0, tend, abstol, reltol, il)
    names(du0) <- names(prob$u0)
    if (np) names(dp) <- names(prob$p)
    list(du0 = du0, dp = dp)
  }
  aprob
}

.sens_pvec <- function(template, values) {
  out <- template
  out[] <- as.numeric(values)
  out
}

# Events recorded by the forward solve (.sens_solve_hook); a problem with
# callbacks solved without recording cannot be differentiated through them.
.sens_events_of <- function(sol) {
  env <- .sol_env(sol)
  if (!is.null(env$sens_events)) return(env$sens_events$events)
  cbs <- .as_callback_set(env$prob$kwargs$callback)
  if (length(cbs$continuous) + length(cbs$discrete))
    stop("the problem has events: solve it with `solve(prob, alg, sensealg = ...)` so that they are ",
         "recorded for the adjoint", call. = FALSE)
  list()
}

# QuadratureAdjoint: the parameter gradient by adaptive Gauss-Kronrod
# quadrature of lambda' f_p + g_p over the pieces between jumps, chained
# backward through the recorded events.
.sens_quadrature_dp <- function(adj, events, acc, tr, mu_integrand, cost, prob, pfinal, t0, tend, abstol,
                                reltol, il) {
  aenv <- .sol_env(adj)
  ev_t <- vapply(events, `[[`, 0, "t")
  nseg <- length(events)
  bounds <- c(t0, ev_t, tend)
  pseg <- c(lapply(events, `[[`, "p"), list(pfinal))
  mu <- 0 * as.numeric(pfinal)
  recs <- acc$records
  for (s in nseg:0) {
    a <- bounds[s + 1L]; b <- bounds[s + 2L]
    cuts <- sort(unique(c(a, cost$t[cost$t > a & cost$t < b], b)))
    tr$seg <- s
    ps <- .sens_pvec(pfinal, pseg[[s + 1L]])
    lam_at <- function(x) .interp_states(aenv, x)[1L, il]
    for (k in seq_len(length(cuts) - 1L)) {
      if (cuts[k + 1L] > cuts[k])
        mu <- mu + .sens_quadgk(function(x) mu_integrand(lam_at(x), ps, x), cuts[k], cuts[k + 1L],
                                abstol, reltol)
    }
    for (rec in recs) if (identical(rec$seg, s) && is.null(rec$A)) mu <- mu + rec$add
    if (s > 0L) {
      rec <- Filter(function(r) identical(r$seg, s) && !is.null(r$A), recs)
      for (r in rec) mu <- as.numeric(r$A %*% mu) + r$add
    }
  }
  mu
}
