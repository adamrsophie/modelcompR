# The high level sensitivity interface, after SciMLSensitivity.jl v7.119
# src/concrete_solve.jl (_concrete_solve_adjoint, _concrete_solve_forward,
# automatic_sensealg_choice, the initialization pullback), src/
# steadystate_adjoint.jl (SteadyStateAdjointProblem), src/second_order.jl
# (_second_order_sensitivities, _second_order_sensitivity_product),
# src/sensitivity_interface.jl (second_order_sensitivities,
# second_order_sensitivity_product) and src/callback_tracking.jl
# (track_callbacks, TrackedAffect).
#
# In Julia the gradient of a loss of a solution is taken by an automatic
# differentiation package over `solve(prob, alg; sensealg)`; the sensealg
# chooses the rule used for the solve. In R, `solve_gradient(loss, prob,
# alg, ..., sensealg)` plays that role: the loss is an R function of the
# solution, its derivative with respect to the saved values (Zygote's
# cotangent of the solution) is taken by central differences, and the rule
# of the sensealg pulls it back to the initial values and parameters.

# solve(prob, alg, sensealg = ...) on ODE problems (hook of solve.ODEProblem):
# records the events for the adjoint (track_callbacks), keeps a dense forward
# solution unless the method checkpoints, and stores the sensealg.
.sens_solve_hook <- function(prob, alg, opts) {
  sa <- opts$sensealg
  if (!.is_sensealg(sa)) stop("`sensealg` must be a sensitivity algorithm such as GaussAdjoint()", call. = FALSE)
  opts$sensealg <- NULL
  rec <- new.env(parent = emptyenv())
  rec$events <- list()
  cbs <- .as_callback_set(opts$callback)
  if (length(cbs$continuous) + length(cbs$discrete))
    opts$callback <- .sens_track_callbacks(cbs, rec, length(prob$u0))
  if (is.null(opts$dense) && !(isTRUE(sa$checkpointing) && !inherits(sa, "BacksolveAdjoint")))
    opts$dense <- TRUE
  sol <- .integrate(prob, alg, opts)
  env <- .sol_env(sol)
  env$sensealg <- sa
  env$sens_events <- rec
  sol
}
.mtk$hooks$sensitivity_solve <- .sens_solve_hook

# track_callbacks: each affect records the state and parameters before and
# after it, so the adjoint can differentiate it.
.sens_track_callbacks <- function(cbs, rec, n) {
  wrap <- function(aff, cb, continuous) {
    if (is.null(aff)) return(NULL)
    force(aff)
    function(integ, ...) {
      ev <- .sens_event_record(integ, aff, list(...), cb, continuous, n)
      aff(integ, ...)
      ev$u_right <- as.numeric(integ$u)
      ev$p_right <- integ$p
      ev$terminate <- isTRUE(integ$terminate)
      rec$events[[length(rec$events) + 1L]] <- ev
      invisible(NULL)
    }
  }
  structure(list(
    continuous = lapply(cbs$continuous, function(cb) {
      cb2 <- cb
      cb2$affect <- wrap(cb$affect, cb, TRUE)
      cb2$affect_neg <- wrap(cb$affect_neg, cb, TRUE)
      cb2
    }),
    discrete = lapply(cbs$discrete, function(cb) {
      cb2 <- cb
      cb2$affect <- wrap(cb$affect, cb, FALSE)
      cb2
    })), class = "CallbackSet")
}

# The automatic choice of automatic_sensealg_choice: forward for small
# problems, the Gauss adjoint otherwise.
.sens_default_alg <- function(prob) {
  if (inherits(prob, c("AbstractNonlinearProblem", "SteadyStateProblem"))) return(SteadyStateAdjoint())
  # the KKT adjoint when constraints or bounds enter the optimality conditions
  if (inherits(prob, "OptimizationProblem"))
    return(if (is.null(prob$f$cons) && is.null(prob$lb) && is.null(prob$ub)) UnconstrainedOptimizationAdjoint()
           else OptimizationAdjoint())
  if (length(prob$u0) + length(prob$p) <= 100L) ForwardDiffSensitivity() else GaussAdjoint()
}

#' Gradient of a loss of a solution
#'
#' `solve_gradient(loss, prob, alg, ..., sensealg)` solves `prob` with
#' `solve(prob, alg, ...)` and returns the gradient of `loss(sol)` with
#' respect to the initial values and the parameters of the problem. It is
#' the R counterpart of differentiating a solve in Julia,
#' `Zygote.gradient((u0, p) -> loss(solve(remake(prob; u0, p), alg; sensealg,
#' kwargs...)), u0, p)`, the "differentiating ODE solutions" of
#' SciMLSensitivity: `sensealg` selects the rule by which the solve is
#' differentiated, as the `sensealg` keyword of Julia's `solve`.
#'
#' The loss receives the solution and may use its saved values (`sol$u`,
#' `sol[x]`, `sol[expr]`, `as.matrix(sol)`) and parameters (`sol$ps[k]`);
#' give the times it needs as `saveat`. It may also interpolate values
#' `sol(t)`: the solve is repeated with those times saved (added to `saveat`,
#' or as stops without `saveat`), hidden from the saved values the loss
#' reads, so they are differentiated like the others; the returned `sol`
#' includes them.
#' Its derivatives with respect to the saved values and the parameters are
#' taken by central differences, then pulled back by
#' * forward sensitivities ([ForwardSensitivity()],
#'   [ForwardDiffSensitivity()], the default when `length(u0) +
#'   length(p) <= 100`),
#' * an adjoint ([InterpolatingAdjoint()], [BacksolveAdjoint()],
#'   [QuadratureAdjoint()], [GaussAdjoint()] (the default for larger
#'   problems), [GaussKronrodAdjoint()]) through [adjoint_sensitivities()],
#' * [SteadyStateAdjoint()] for `NonlinearProblem` and `SteadyStateProblem`
#'   (the implicit function theorem),
#' * pathwise forward sensitivities for `SDEProblem` (with a fixed step
#'   method, a `dt` and a `seed`, the derivative of the sampled path), or
#'   for `SDEProblem` and `RODEProblem` the adjoints [BacksolveAdjoint()]
#'   and [InterpolatingAdjoint()] of [SDEAdjointProblem()].
#'
#' Events are differentiated, including the dependence of event times on
#' the parameters. When the problem was built from a System whose initial
#' values depend on the parameters through initialization, the gradient
#' with respect to `p` includes that dependence.
#'
#' @param loss A function of the solution returning a number.
#' @param prob An `ODEProblem`, `NonlinearProblem`, `SteadyStateProblem`,
#'   `OptimizationProblem`, `SDEProblem` or `RODEProblem`.
#' @param alg The algorithm for [solve()].
#' @param ... Options for [solve()] (`saveat`, `abstol`, `reltol`, `dt`,
#'   `seed`, ...); the adjoint is solved with the same tolerances.
#' @param sensealg A sensitivity algorithm; chosen automatically when
#'   `NULL` (or taken from the problem's `sensealg` option).
#' @return A list with the gradients `u0` and `p`, the loss value `value`
#'   and the solution `sol`.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
#' loss <- function(sol) sum(as.matrix(sol))
#' g1 <- solve_gradient(loss, prob, Tsit5(), saveat = 0.25, reltol = 1e-10, abstol = 1e-12)
#' g2 <- solve_gradient(loss, prob, Tsit5(), saveat = 0.25, reltol = 1e-10, abstol = 1e-12,
#'                      sensealg = GaussAdjoint())
#' c(g1$p, g2$p)
#' @export
solve_gradient <- function(loss, prob, alg = NULL, ..., sensealg = NULL) {
  # a loss that interpolates sol(t): solve again saving those times (as
  # `saveat` when it is given, else as stops), hidden from the saved values
  # the loss sees
  old <- list(.mtk$sens_interp_t, .mtk$sens_hide_t)
  on.exit({
    .mtk$sens_interp_t <- old[[1L]]
    .mtk$sens_hide_t <- old[[2L]]
  })
  dots <- list(...)
  base_saveat <- dots$saveat
  tq <- numeric()
  for (attempt in 1:4) {
    r <- tryCatch(.solve_gradient(loss, prob, alg, dots, sensealg),
                  sens_needs_times = function(e) e)
    if (!inherits(r, "sens_needs_times")) return(r)
    if (attempt == 4L) stop("the loss interpolates the solution at changing times", call. = FALSE)
    tq <- sort(unique(c(tq, as.numeric(r$times))))
    if (!is.null(base_saveat)) {
      sv <- if (length(base_saveat) == 1L) seq(prob$tspan[1L], prob$tspan[2L], by = base_saveat *
                                                 sign(prob$tspan[2L] - prob$tspan[1L])) else base_saveat
      dots$saveat <- sort(unique(c(sv, tq)), decreasing = prob$tspan[2L] < prob$tspan[1L])
      .mtk$sens_hide_t <- tq[!vapply(tq, function(x) any(abs(sv - x) <= 1e-12 * max(1, abs(x))), TRUE)]
    } else {
      dots$tstops <- sort(unique(c(dots$tstops, tq)))
      .mtk$sens_hide_t <- tq
    }
    .mtk$sens_interp_t <- tq
  }
}

.solve_gradient <- function(loss, prob, alg, dots, sensealg) {
  sensealg <- sensealg %||% prob$kwargs$sensealg %||% .sens_default_alg(prob)
  if (!.is_sensealg(sensealg)) stop("`sensealg` must be a sensitivity algorithm", call. = FALSE)
  if (inherits(prob, c("AbstractNonlinearProblem", "SteadyStateProblem")))
    return(.sens_nl_gradient(loss, prob, alg, dots, sensealg))
  if (inherits(prob, "OptimizationProblem")) return(.sens_opt_gradient(loss, prob, alg, dots, sensealg))
  # the SDE and RODE adjoints (R/sensitivity-sde.R)
  if (inherits(prob, c("SDEProblem", "RODEProblem")) && inherits(sensealg, c("BacksolveAdjoint", "InterpolatingAdjoint")))
    return(.sens_stochastic_gradient(loss, prob, alg, dots, sensealg))
  if (inherits(prob, "SDEProblem")) return(.sens_sde_gradient(loss, prob, alg, dots, sensealg))
  if (!inherits(prob, "ODEProblem"))
    stop("solve_gradient() supports ODE, nonlinear, steady state, optimization and SDE problems, and RODE problems ",
         "with BacksolveAdjoint() or InterpolatingAdjoint()", call. = FALSE)
  if (.is_sensealg(sensealg, "ShadowingSensitivityAlgorithm"))
    stop("shadowing methods differentiate long time averages: use ForwardLSSProblem(), NILSSProblem(), ",
         "AdjointLSSProblem() or NILSASProblem()", call. = FALSE)
  n <- length(prob$u0); np <- length(prob$p)
  if (!is.null(dots$callback)) {
    # callbacks given to solve() join those of the problem, so the
    # sensitivity methods see and differentiate them
    prob$kwargs$callback <- CallbackSet(prob$kwargs$callback, dots$callback)
    dots$callback <- NULL
  }
  ini <- .sens_init_jacobian(prob)
  if (.is_sensealg(sensealg, "ForwardSensitivityAlgorithm")) {
    # theta = (u0, p): the initial values directly and through initialization
    S0 <- cbind(diag(1, n), ini$du0)
    Q0 <- cbind(matrix(0, np, n), ini$dp)
    fprob <- .sens_forward_problem(prob, S0, Q0, sensealg)
    fsol <- do.call(solve, c(list(fprob, alg), dots))
    fenv <- .sol_env(fsol)
    m <- n + np
    U <- fenv$U
    sol <- .new_ode_solution(t = fenv$t, u = U[, seq_len(n), drop = FALSE], prob = prob, alg = fenv$alg,
                             retcode = fenv$retcode, stats = fenv$stats, p = fenv$p_final)
    senv <- .sol_env(sol)
    if (!is.null(fenv$param_history)) senv$param_history <- fenv$param_history
    ct <- .sens_loss_cotangent(loss, sol)
    grad <- numeric(m)
    Ss <- lapply(seq_len(nrow(U)), function(r) matrix(U[r, n + seq_len(n * m)], n, m))
    # samples saved at an event time move with it: the left limit (the
    # first of two samples) changes by f- dtau, the right limit by f+ dtau
    fixed <- .sens_fixed_times(fenv$t, dots$saveat, prob$tspan)
    for (ev in fprob$sensitivity$events$events) {
      rows <- which(abs(fenv$t - ev$t) <= 1e-12 * max(1, abs(ev$t)) & !fixed)
      if (!length(rows)) next
      r_right <- rows[length(rows)]
      Ss[[r_right]] <- Ss[[r_right]] + outer(ev$fp, ev$dtau)
      if (length(rows) >= 2L) Ss[[rows[1L]]] <- Ss[[rows[1L]]] + outer(ev$fm, ev$dtau)
      grad <- grad + sum(ct$dT[rows]) * ev$dtau
    }
    for (r in seq_len(nrow(U))) grad <- grad + as.numeric(ct$dU[r, ] %*% Ss[[r]])
    if (np) {
      Qend <- if (fprob$sensitivity$nq) matrix(U[nrow(U), n + n * m + seq_len(np * m)], np, m) else Q0
      grad <- grad + as.numeric(ct$dp %*% Qend)
      # parameters in force before a later change: Q0 for the first, else the
      # Q after the event that set them
      if (!is.null(ct$dp_hist)) {
        evs <- fprob$sensitivity$events$events
        for (k in seq_len(nrow(ct$dp_hist))) {
          Qk <- if (k == 1L) Q0 else {
            e <- which(vapply(evs, function(ev) abs(ev$t - ct$t_hist[k]) <= 1e-12 * max(1, abs(ev$t)), TRUE))
            evs[[e[length(e)]]]$Q_right
          }
          grad <- grad + as.numeric(ct$dp_hist[k, ] %*% Qk)
        }
      }
    }
    return(list(u0 = stats::setNames(grad[seq_len(n)], names(prob$u0)),
                p = stats::setNames(grad[n + seq_len(np)], names(prob$p)), value = ct$value, sol = sol))
  }
  if (!.is_sensealg(sensealg, "AdjointSensitivityAlgorithm"))
    stop("`sensealg` must be a forward or adjoint sensitivity algorithm", call. = FALSE)
  sol <- do.call(solve, c(list(prob, alg), dots, list(sensealg = sensealg)))
  env <- .sol_env(sol)
  ct <- .sens_loss_cotangent(loss, sol)
  nt <- length(env$t)
  last <- nt
  # the backward solve takes the step options of the forward one
  step_opts <- dots[intersect(names(dots), c("dt", "dtmax", "dtmin", "maxiters", "adaptive"))]
  tsave <- env$t
  fixed <- .sens_fixed_times(env$t, dots$saveat, prob$tspan)
  # the parameters in force before a later change enter as parameter costs
  # at a time inside their segment
  hp <- NULL
  if (!is.null(ct$dp_hist)) {
    tb <- c(ct$t_hist, if (identical(env$retcode, "Terminated")) env$t[nt] else prob$tspan[2L])
    hp <- vapply(seq_len(nrow(ct$dp_hist)), function(k) (tb[k] + tb[k + 1L]) / 2, 0)
    tsave <- c(tsave, hp)
    fixed <- c(fixed, rep(TRUE, length(hp)))
  }
  attr(tsave, "fixed") <- fixed
  res <- do.call(adjoint_sensitivities, c(list(sol, alg, t = tsave,
                               dgdu_discrete = function(u, p, tt, i) if (i <= nt) ct$dU[i, ] else numeric(n),
                               dgdp_discrete = if (np && (any(ct$dp != 0) || !is.null(hp))) function(u, p, tt, i)
                                 if (i == last) ct$dp else if (i > nt) ct$dp_hist[i - nt, ] else numeric(np),
                               abstol = dots$abstol %||% 1e-6, reltol = dots$reltol %||% 1e-3,
                               sensealg = sensealg,
                               dgdt_discrete = if (any(ct$dT != 0)) function(tt)
                                 sum(ct$dT[abs(env$t - tt) <= 1e-12 * max(1, abs(tt)) & !fixed[seq_len(nt)]]),
                               no_start = !isTRUE(abs(env$t[1L] - prob$tspan[1L]) <= 1e-14 * max(1, abs(prob$tspan[1L])))),
                    step_opts))
  dp <- res$dp
  if (np && ini$active) dp <- as.numeric(crossprod(ini$dp, dp) + crossprod(ini$du0, res$du0))
  list(u0 = stats::setNames(res$du0, names(prob$u0)), p = stats::setNames(dp, names(prob$p)),
       value = ct$value, sol = sol)
}

# Saved times that were requested with `saveat`: they stay fixed when an
# event at the same time moves (the samples of the event itself move).
.sens_fixed_times <- function(ts, saveat, tspan) {
  if (is.null(saveat)) return(rep(FALSE, length(ts)))
  sa <- as.numeric(saveat)
  if (length(sa) == 1L) sa <- c(seq(tspan[1L], tspan[2L], by = sign(tspan[2L] - tspan[1L]) * abs(sa)), tspan[2L])
  tol <- 1e-12 * max(1, abs(tspan))
  vapply(ts, function(x) any(abs(sa - x) <= tol), TRUE)
}

# The derivative of the loss with respect to the saved values (dU, one row
# per saved time) and to the parameters, by central differences of the loss
# on copies of the solution that hold the perturbed values.
.sens_loss_cotangent <- function(loss, sol) {
  env <- .sol_env(sol)
  data <- mget(setdiff(ls(env, all.names = TRUE), c("obs_cache", "steps")), envir = env)
  data$steps <- NULL
  data$dense <- FALSE
  data$sens_frozen <- TRUE
  cls <- setdiff(class(sol), c("SciMLSolution", "function"))
  U0 <- env$U
  p0 <- env$p_final
  # rows saved only for the interpolation of the loss (solve_gradient())
  hide <- .mtk$sens_hide_t
  hidden <- if (length(hide)) vapply(env$t, function(x) any(abs(hide - x) <= 1e-12 * max(1, abs(x))), TRUE) else
    logical(length(env$t))
  itimes <- .mtk$sens_interp_t
  irows <- if (length(itimes)) vapply(itimes, function(x) {
    i <- which(abs(env$t - x) <= 1e-12 * max(1, abs(x)))
    if (length(i)) i[length(i)] else NA_integer_
  }, 0L) else integer()
  ev <- function(U, p) {
    d <- data
    if (length(irows) && !anyNA(irows)) d$sens_interp <- list(t = itimes, U = U[irows, , drop = FALSE])
    d$U <- U[!hidden, , drop = FALSE]
    d$t <- env$t[!hidden]
    d$p_final <- p
    if (!identical(p, p0) && !is.null(d$param_history)) {
      d$param_history$p[nrow(d$param_history$p), ] <- p
    }
    if (!identical(p, p0) && !is.null(d$prob)) d$prob$p <- p
    v <- loss(.new_solution(cls, d))
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) stop("the loss must return a finite number", call. = FALSE)
    v
  }
  value <- ev(U0, p0)
  dU <- U0
  for (k in seq_along(U0)) {
    h <- .sens_h(U0[k])
    Up <- U0; Um <- U0
    Up[k] <- U0[k] + h
    Um[k] <- U0[k] - h
    dU[k] <- (ev(Up, p0) - ev(Um, p0)) / (2 * h)
  }
  dp <- if (length(p0)) .sens_fd_grad(function(q) ev(U0, .sens_pvec(p0, q)), as.numeric(p0)) else numeric()
  # the cotangent of the saved times (dgdt_discrete): it only matters for
  # samples that move with an event
  tvis <- env$t[!hidden]
  dT <- numeric(length(env$t))
  ev_t <- function(tv) {
    d <- data
    if (length(irows) && !anyNA(irows)) d$sens_interp <- list(t = itimes, U = U0[irows, , drop = FALSE])
    d$U <- U0[!hidden, , drop = FALSE]
    d$t <- tv
    loss(.new_solution(cls, d))
  }
  vis <- which(!hidden)
  for (j in seq_along(tvis)) {
    hh <- 1e-6 * max(1, abs(tvis[j]))
    tp <- tvis; tm <- tvis
    tp[j] <- tp[j] + hh; tm[j] <- tm[j] - hh
    dT[vis[j]] <- (ev_t(tp) - ev_t(tm)) / (2 * hh)
  }
  # the cotangents of the parameters in force before the last change (rows
  # of the parameter history, in integration order); the sensitivity
  # methods pull each back through its own segment
  h <- data$param_history
  dp_hist <- NULL
  if (length(p0) && !is.null(h) && nrow(h$p) > 1L) {
    K <- nrow(h$p) - 1L
    dp_hist <- t(vapply(seq_len(K), function(k) .sens_fd_grad(function(q) {
      d <- data
      if (length(irows) && !anyNA(irows)) d$sens_interp <- list(t = itimes, U = U0[irows, , drop = FALSE])
      d$U <- U0[!hidden, , drop = FALSE]
      d$t <- env$t[!hidden]
      d$param_history$p[k, ] <- d$param_history$p[k, ] + q
      loss(.new_solution(cls, d))
    }, numeric(length(p0))), numeric(length(p0))))
    if (length(p0) == 1L) dp_hist <- matrix(dp_hist, ncol = 1L)
    if (all(abs(dp_hist) <= 1e-12 * max(1, abs(value)))) dp_hist <- NULL
  }
  list(value = value, dU = unname(dU), dp = dp, dp_hist = dp_hist, t_hist = h$t, dT = dT)
}

.sens_no_hist <- function(ct) {
  if (!is.null(ct$dp_hist))
    stop("the loss depends on parameter values in force before an event changed them; for SDEs this ",
         "dependence is not differentiated", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Steady states and nonlinear problems: the implicit function theorem.
# With F(u*, p) = 0, du*/dp = -F_u^-1 F_p; the adjoint lambda = F_u'^-1 g_u
# gives dG/dp = g_p - lambda' F_p (SteadyStateAdjointProblem).

.sens_nl_fns <- function(prob, sensealg) {
  ss <- inherits(prob, "SteadyStateProblem")
  f <- prob$f$f
  F <- if (ss) function(u, p) as.numeric(f(u, p, Inf)) else function(u, p) as.numeric(f(u, p))
  sys <- prob$sys
  if (!is.null(sys) && !identical(sensealg$autojacvec, FALSE)) {
    d <- .sens_symbolic(sys)
    np <- length(prob$p)
    return(list(F = F, jac = function(u, p) d$jac(u, p, Inf),
                pjac = function(u, p) if (np) d$pjac(u, p, Inf) else matrix(0, length(u), 0L)))
  }
  uj <- prob$f$jac
  list(F = F,
       jac = if (is.function(uj)) function(u, p) as.matrix(if (ss) uj(u, p, Inf) else uj(u, p)) else
         function(u, p) .sens_fd_jac(function(z) F(z, p), u),
       pjac = function(u, p) .sens_fd_jac(function(q) F(u, .sens_pvec(p, q)), as.numeric(p)))
}

.sens_steady_adjoint <- function(sol, sensealg, dgdu = NULL, dgdp = NULL, g = NULL, ...) {
  env <- .sol_env(sol)
  prob <- env$prob
  if (is.null(dgdu) && is.null(dgdp) && is.null(g)) stop("Either `dgdu`, `dgdp`, or `g` must be specified.",
                                                         call. = FALSE)
  if (!length(prob$p))
    stop("Your model does not have parameters, and thus it is impossible to calculate the derivative of ",
         "the solution with respect to the parameters.", call. = FALSE)
  y <- as.numeric(env$U[1L, ])
  p <- env$p_final %||% prob$p
  nf <- .sens_nl_fns(prob, sensealg)
  gu <- if (!is.null(dgdu)) as.numeric(dgdu(y, p)) else if (!is.null(g))
    .sens_fd_grad(function(z) g(z, p), y) else numeric(length(y))
  gp <- if (!is.null(dgdp)) as.numeric(dgdp(y, p)) else if (!is.null(g))
    .sens_fd_grad(function(q) g(y, .sens_pvec(p, q)), as.numeric(p)) else numeric(length(p))
  dp <- .sens_ift_pullback(nf, y, p, gu, gp)
  list(du0 = stats::setNames(numeric(length(y)), names(prob$u0)), dp = stats::setNames(dp, names(prob$p)))
}

.sens_ift_pullback <- function(nf, y, p, gu, gp) {
  J <- nf$jac(y, p)
  if (nrow(J) != ncol(J)) stop("SteadyStateAdjoint needs as many equations as unknowns", call. = FALSE)
  lam <- if (length(y)) solve(t(J), gu) else numeric()
  gp - as.numeric(crossprod(nf$pjac(y, p), lam))
}

.sens_nl_gradient <- function(loss, prob, alg, dots, sensealg) {
  if (!inherits(sensealg, "SteadyStateAdjoint"))
    stop("nonlinear and steady state problems are differentiated with SteadyStateAdjoint()", call. = FALSE)
  sol <- do.call(solve, c(list(prob, alg), dots))
  ct <- .sens_loss_cotangent(loss, sol)
  env <- .sol_env(sol)
  y <- as.numeric(env$U[1L, ])
  dp <- .sens_ift_pullback(.sens_nl_fns(prob, sensealg), y, prob$p, as.numeric(ct$dU[1L, ]), ct$dp)
  list(u0 = stats::setNames(numeric(length(y)), names(prob$u0)), p = stats::setNames(dp, names(prob$p)),
       value = ct$value, sol = sol)
}

# ---------------------------------------------------------------------------
# SDEs: pathwise forward sensitivities. The augmented SDE
#   dS = (f_u S + f_p Q) dt + (g_u S + g_p Q) dW
# driven by the same Wiener increments is, for fixed step methods such as
# EM(), the derivative of the computed path (what ForwardDiffSensitivity
# gives in Julia).

.sens_sde_gradient <- function(loss, prob, alg, dots, sensealg) {
  if (!.is_sensealg(sensealg, "ForwardSensitivityAlgorithm"))
    stop("SDE problems are differentiated pathwise with ForwardDiffSensitivity() or with the adjoints ",
         "BacksolveAdjoint() and InterpolatingAdjoint()", call. = FALSE)
  n <- length(prob$u0); np <- length(prob$p); m <- n + np
  f <- prob$f$f; gfun <- prob$f$g
  type <- prob$noise_type %||% prob$f$noise_type %||% "diagonal"
  G <- function(u, p, t) {
    v <- gfun(u, p, t)
    if (type == "diagonal") diag(as.numeric(v), n) else matrix(v, nrow = n)
  }
  nw <- ncol(G(prob$u0, prob$p, prob$tspan[1L]))
  fp <- prob
  fp$kwargs <- list()
  fns <- .sens_fns(fp, ForwardSensitivity(autojacvec = FALSE))
  ini <- .sens_init_jacobian(prob)
  S0 <- cbind(diag(1, n), ini$du0)
  Q0 <- cbind(matrix(0, np, n), ini$dp)
  # discrete events map (u, S, Q) through their linearization, so the
  # parameter sensitivities Q are carried as states when there are events
  cbs <- .as_callback_set(dots$callback %||% prob$kwargs$callback)
  if (length(cbs$continuous))
    stop("continuous callbacks are not supported by the SDE sensitivities: the time at which a Brownian ",
         "path crosses a level is not differentiable; use a DiscreteCallback", call. = FALSE)
  has_cb <- length(cbs$discrete) > 0L
  nq <- if (has_cb) np * m else 0L
  Qof <- function(z) if (has_cb) matrix(z[n + n * m + seq_len(nq)], np, m) else Q0
  drift <- function(z, p, t) {
    u <- z[seq_len(n)]
    S <- matrix(z[n + seq_len(n * m)], n, m)
    dS <- fns$jac(u, p, t) %*% S + fns$pjac(u, p, t) %*% Qof(z)
    c(f(u, p, t), as.numeric(dS), numeric(nq))
  }
  diffusion <- function(z, p, t) {
    u <- z[seq_len(n)]
    S <- matrix(z[n + seq_len(n * m)], n, m)
    Gu <- .sens_fd_jac(function(x) as.numeric(G(x, p, t)), u)            # (n nw) x n
    Gp <- if (np) .sens_fd_jac(function(q) as.numeric(G(u, .sens_pvec(p, q), t)), as.numeric(p)) else
      matrix(0, n * nw, 0L)
    dG <- Gu %*% S + Gp %*% Qof(z)                                          # (n nw) x m
    rbind(G(u, p, t), do.call(rbind, lapply(seq_len(m), function(j) matrix(dG[, j], n, nw))),
          matrix(0, nq, nw))
  }
  z0 <- c(as.numeric(prob$u0), as.numeric(S0), if (nq) as.numeric(Q0))
  # a mass matrix multiplies the drift of u and of each sensitivity column
  M <- prob$f$mass_matrix
  Maug <- if (is.matrix(M) && !.is_identity_mass(M))
    .sens_blockdiag(c(list(M), rep(list(M), m), if (nq) list(diag(1, nq))))
  aprob <- SDEProblem(SDEFunction(drift, diffusion, mass_matrix = Maug), z0, prob$tspan, prob$p,
                      noise_rate_prototype = matrix(0, n + n * m + nq, nw))
  aprob$kwargs <- prob$kwargs
  # the same noise path (a NoiseGrid or other given process)
  aprob$noise_process <- prob$noise_process
  if (has_cb) {
    rec <- new.env(parent = emptyenv())
    cb <- .sens_forward_callbacks(cbs, fp, fns, n, np, m, rec)
    aprob$kwargs$callback <- cb
    dots$callback <- cb
  }
  asol <- do.call(solve, c(list(aprob, alg), dots))
  aenv <- .sol_env(asol)
  U <- aenv$U
  sol <- .new_ode_solution(t = aenv$t, u = U[, seq_len(n), drop = FALSE], prob = prob, alg = aenv$alg,
                           retcode = aenv$retcode, stats = aenv$stats, p = prob$p)
  ct <- .sens_loss_cotangent(loss, sol)
  .sens_no_hist(ct)
  grad <- numeric(m)
  for (r in seq_len(nrow(U))) grad <- grad + as.numeric(ct$dU[r, ] %*% matrix(U[r, n + seq_len(n * m)], n, m))
  if (np) grad <- grad + as.numeric(ct$dp %*% Qof(U[nrow(U), ]))
  list(u0 = stats::setNames(grad[seq_len(n)], names(prob$u0)),
       p = stats::setNames(grad[n + seq_len(np)], names(prob$p)), value = ct$value, sol = sol)
}

# ---------------------------------------------------------------------------
# Second order

#' Second order sensitivities
#'
#' `second_order_sensitivities(loss, prob, alg, ...)` returns the Hessian of
#' `loss(solve(prob, alg, ...))` with respect to the parameters, and
#' `second_order_sensitivity_product(loss, v, prob, alg, ...)` its product
#' with `v`. Julia's `ForwardDiffOverAdjoint` runs forward mode automatic
#' differentiation over an adjoint gradient; here the adjoint gradient of
#' [solve_gradient()] is differentiated by central differences, so the
#' solves should use tight tolerances.
#'
#' @param loss A function of the solution, as for [solve_gradient()].
#' @param v A vector of the length of the parameters.
#' @param prob,alg The problem and algorithm.
#' @param ... Options for [solve()].
#' @param sensealg [ForwardDiffOverAdjoint()] with a first order method.
#' @return A matrix (the Hessian) or a vector (its product with `v`).
#' @examples
#' prob <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
#' loss <- function(sol) sol$u[[length(sol$t)]]
#' H <- second_order_sensitivities(loss, prob, Tsit5(), saveat = 1, reltol = 1e-12, abstol = 1e-12)
#' H - exp(-2)   # d2 u(1) / dp2 = exp(-2)
#' @export
second_order_sensitivities <- function(loss, prob, alg, ..., sensealg = ForwardDiffOverAdjoint()) {
  np <- length(prob$p)
  H <- sapply(seq_len(np), function(j) {
    e <- numeric(np)
    e[j] <- 1
    .sens_grad_direction(loss, prob, alg, e, sensealg, list(...))
  })
  H <- matrix(H, np, np, dimnames = list(names(prob$p), names(prob$p)))
  (H + t(H)) / 2
}

#' @rdname second_order_sensitivities
#' @export
second_order_sensitivity_product <- function(loss, v, prob, alg, ..., sensealg = ForwardDiffOverAdjoint()) {
  stats::setNames(.sens_grad_direction(loss, prob, alg, as.numeric(v), sensealg, list(...)), names(prob$p))
}

.sens_grad_direction <- function(loss, prob, alg, v, sensealg, dots) {
  p0 <- prob$p
  h <- 1e-4 * max(1, sqrt(sum(p0^2))) / max(1, sqrt(sum(v^2)))
  gr <- function(s) {
    pr <- suppressWarnings(remake(prob, p = as.numeric(p0) + s * v))
    do.call(solve_gradient, c(list(loss, pr, alg), dots, list(sensealg = sensealg$adjalg)))$p
  }
  (gr(h) - gr(-h)) / (2 * h)
}
