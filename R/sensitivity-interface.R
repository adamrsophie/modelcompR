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
      .cb_rebind_init(cb2)
    }),
    discrete = lapply(cbs$discrete, function(cb) {
      cb2 <- cb
      cb2$affect <- wrap(cb$affect, cb, FALSE)
      .cb_rebind_init(cb2)
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
#' @param prob An `ODEProblem`, `DAEProblem`, `NonlinearProblem`,
#'   `SteadyStateProblem`, `OptimizationProblem`, `SDEProblem` or
#'   `RODEProblem`. A `DAEProblem` `F(du, u, p, t) = 0` of index 1 is solved
#'   and differentiated as the mass matrix problem `u_d' = w, 0 = F(w, u, p,
#'   t)` (with [FBDF()] for [DFBDF()], [Rodas5P()] for other DAE methods), so
#'   its events and every sensitivity method apply.
#' @param alg The algorithm for [solve()].
#' @param ... Options for [solve()] (`saveat`, `abstol`, `reltol`, `dt`,
#'   `seed`, ...); the adjoint is solved with the same tolerances.
#' @param sensealg A sensitivity algorithm; chosen automatically when
#'   `NULL` (or taken from the problem's `sensealg` option).
#' @param allow_failure By default a solve that does not succeed (a retcode
#'   other than `Success` or `Terminated`, such as `MaxIters`) is an error;
#'   `TRUE` differentiates the trajectory that was computed.
#' @param loss_grad Optional derivative of the loss: a function of the
#'   solution returning the derivative with respect to the saved values (a
#'   matrix shaped as `as.matrix(sol)`, or a list with that matrix as `u`
#'   and the derivative with respect to the final parameters as `p`). It
#'   replaces the finite differences of the loss in the saved values. With
#'   `save_idxs` the matrix has the saved components only, as the loss sees
#'   them. Values the loss reads by interpolation, `sol(t)`, are not saved
#'   values: they are still differenced, so `loss_grad` gives the
#'   derivative in the saved values alone.
#' @details The `value` is the loss of the solution the gradient
#'   differentiates: the forward solve itself for the adjoint methods, and
#'   the solve of the state with its sensitivities for the forward methods,
#'   whose step sizes also control the error of the sensitivities (as in
#'   SciMLSensitivity), so it can differ from `loss(solve(...))` within the
#'   tolerances. A failed numerical stage (the forward solve, the backward
#'   adjoint solve or a checkpoint replay) is an error naming the stage, its
#'   retcode and the time where it stopped, unless `allow_failure = TRUE`.
#' @param loss_rows Optional declaration of the saved rows the loss reads
#'   (indices of the rows of `as.matrix(sol)`, or `"last"` for a terminal
#'   loss): the loss is differentiated in those rows only, and treated as
#'   independent of the others. Without it every saved value is
#'   differentiated.
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
solve_gradient <- function(loss, prob, alg = NULL, ..., sensealg = NULL, allow_failure = FALSE,
                           loss_grad = NULL, loss_rows = NULL) {
  # a loss that interpolates sol(t): solve again saving those times (as
  # `saveat` when it is given, else as stops), hidden from the saved values
  # the loss sees
  old <- list(.mtk$sens_interp_t, .mtk$sens_hide_t, .mtk$sens_loss_grad, .mtk$sens_loss_rows)
  on.exit({
    .mtk$sens_interp_t <- old[[1L]]
    .mtk$sens_hide_t <- old[[2L]]
    .mtk$sens_loss_grad <- old[[3L]]
    .mtk$sens_loss_rows <- old[[4L]]
  })
  .mtk$sens_loss_rows <- loss_rows
  old_af <- .mtk$sens_allow_failure
  on.exit(.mtk$sens_allow_failure <- old_af, add = TRUE)
  .mtk$sens_allow_failure <- allow_failure
  if (!is.null(loss_grad) && !is.function(loss_grad)) stop("`loss_grad` must be a function of the solution", call. = FALSE)
  .mtk$sens_loss_grad <- loss_grad
  dots <- list(...)
  # the options of the solve are those given here over those of the problem
  if (identical(dots$save_idxs_storage %||% prob$kwargs$save_idxs_storage, "reduced"))
    stop("solve_gradient() needs every state of the solution; use save_idxs_storage = \"full\"", call. = FALSE)
  # save_idxs: the loss sees the saved components, as in solve(); the
  # sensitivity solves keep every state and the cotangents are scattered
  # back to the selected states
  si <- dots$save_idxs %||% prob$kwargs$save_idxs
  old_si <- .mtk$sens_save_idxs
  on.exit(.mtk$sens_save_idxs <- old_si, add = TRUE)
  .mtk$sens_save_idxs <- if (!is.null(si)) get_save_idxs_and_saved_subsystem(prob, si)$save_idxs
  base_saveat <- dots$saveat %||% prob$kwargs$saveat
  tq <- numeric()
  for (attempt in 1:4) {
    r <- tryCatch(.solve_gradient(loss, prob, alg, dots, sensealg),
                  sens_needs_times = function(e) e)
    if (!inherits(r, "sens_needs_times")) {
      # the derivative of a solve that stopped early (MaxIters, DtLessThanMin,
      # Unstable, ...) is not the derivative of the loss over the time span
      rc <- if (!is.null(r$sol)) r$sol$retcode
      if (!is.null(rc) && !successful_retcode(rc) && !isTRUE(allow_failure))
        stop("the solve ended with retcode ", rc, " before the end of the time span, so its gradient is not ",
             "the gradient of the loss; pass allow_failure = TRUE to differentiate the partial trajectory",
             call. = FALSE)
      # what ran: the method, the arithmetic, how the loss was differentiated
      # and the view it saw
      r$execution <- list(sensealg = class(sensealg %||% prob$kwargs$sensealg %||% .sens_default_alg(prob))[1L],
                          arithmetic = if (.arith$native) "native" else "reference",
                          loss_derivative = if (!is.null(loss_grad)) "loss_grad" else if (!is.null(loss_rows))
                            "finite differences in loss_rows" else "finite differences in every saved value",
                          save_idxs = .mtk$sens_save_idxs, forward_retcode = rc,
                          registrations = prob$registry$version)
      return(r)
    }
    if (attempt == 4L) stop("the loss interpolates the solution at changing times", call. = FALSE)
    tq <- sort(unique(c(tq, as.numeric(r$times))))
    if (!is.null(base_saveat)) {
      # a step saveat also saves the final time (.sens_fixed_times)
      sv <- if (length(base_saveat) == 1L) c(seq(prob$tspan[1L], prob$tspan[2L], by = base_saveat *
                                                   sign(prob$tspan[2L] - prob$tspan[1L])), prob$tspan[2L]) else base_saveat
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
  if (inherits(prob, "DAEProblem")) return(.sens_dae_gradient(loss, prob, alg, dots, sensealg))
  if (!inherits(prob, "ODEProblem"))
    stop("solve_gradient() supports ODE, DAE, nonlinear, steady state, optimization and SDE problems, and RODE problems ",
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
    fixed <- .sens_fixed_times(fenv$t, dots$saveat %||% prob$kwargs$saveat, prob$tspan)
    for (ev in fprob$sensitivity$events$events) {
      rows <- which(abs(fenv$t - ev$t) <= 1e-12 * max(1, abs(ev$t)) & !fixed)
      if (!length(rows)) next
      if (length(rows) == 1L && .sens_left_only(ev)) {
        # the event saves only its left limit
        Ss[[rows]] <- Ss[[rows]] + outer(ev$fm, ev$dtau)
      } else {
        r_right <- rows[length(rows)]
        Ss[[r_right]] <- Ss[[r_right]] + outer(ev$fp, ev$dtau)
        if (length(rows) >= 2L) Ss[[rows[1L]]] <- Ss[[rows[1L]]] + outer(ev$fm, ev$dtau)
      }
      if (any(ev$dtau != 0)) grad <- grad + sum(ct$dT(rows)) * ev$dtau
    }
    for (r in seq_len(nrow(U))) grad <- grad + as.numeric(ct$dU[r, ] %*% Ss[[r]])
    if (np) {
      Qend <- if (fprob$sensitivity$nq) matrix(U[nrow(U), n + n * m + seq_len(np * m)], np, m) else Q0
      # an event after the last saved time changes the final parameters
      evs_all <- fprob$sensitivity$events$events
      if (length(evs_all)) {
        tl <- fenv$t[nrow(U)]
        dir <- sign(prob$tspan[2L] - prob$tspan[1L])
        late <- which(vapply(evs_all, function(ev) dir * (ev$t - tl) > 1e-12 * max(1, abs(tl)), TRUE))
        if (length(late)) Qend <- evs_all[[late[length(late)]]]$Q_right
      }
      grad <- grad + as.numeric(ct$dp %*% Qend)
      if (!is.null(ct$dp_init)) grad[n + seq_len(np)] <- grad[n + seq_len(np)] + ct$dp_init
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
  fixed <- .sens_fixed_times(env$t, dots$saveat %||% prob$kwargs$saveat, prob$tspan)
  # the parameters in force before a later change enter as parameter costs
  # at a time inside their segment
  hp <- NULL
  if (!is.null(ct$dp_hist)) {
    tb <- c(ct$t_hist, if (identical(env$retcode, "Terminated")) env$t[nt] else prob$tspan[2L])
    hp <- vapply(seq_len(nrow(ct$dp_hist)), function(k) (tb[k] + tb[k + 1L]) / 2, 0)
    tsave <- c(tsave, hp)
    fixed <- c(fixed, rep(TRUE, length(hp)))
  }
  # the final parameters enter at the end time (after any event there)
  tfin <- if (identical(env$retcode, "Terminated")) env$t[nt] else prob$tspan[2L]
  if (!isTRUE(abs(env$t[nt] - tfin) <= 1e-12 * max(1, abs(tfin)))) {
    tsave <- c(tsave, tfin)
    fixed <- c(fixed, TRUE)
    last <- length(tsave)
  }
  attr(tsave, "fixed") <- fixed
  res <- do.call(adjoint_sensitivities, c(list(sol, alg, t = tsave,
                               dgdu_discrete = function(u, p, tt, i) if (i <= nt) ct$dU[i, ] else numeric(n),
                               dgdp_discrete = if (np && (any(ct$dp != 0) || !is.null(hp))) function(u, p, tt, i)
                                 if (i == last) ct$dp else if (i > nt && i - nt <= length(hp)) ct$dp_hist[i - nt, ]
                                 else numeric(np),
                               abstol = dots$abstol %||% 1e-6, reltol = dots$reltol %||% 1e-3,
                               sensealg = sensealg,
                               dgdt_discrete = function(tt)
                                 sum(ct$dT(which(abs(env$t - tt) <= 1e-12 * max(1, abs(tt)) & !fixed[seq_len(nt)]))),
                               no_start = !isTRUE(abs(env$t[1L] - prob$tspan[1L]) <= 1e-14 * max(1, abs(prob$tspan[1L])))),
                    step_opts))
  dp <- res$dp
  if (np && ini$active) dp <- as.numeric(crossprod(ini$dp, dp) + crossprod(ini$du0, res$du0))
  if (!is.null(ct$dp_init)) dp <- dp + ct$dp_init
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
  if (!is.null(.mtk$sens_save_idxs) && is.null(env$save_idxs)) env$save_idxs <- .mtk$sens_save_idxs
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
  # sol$prob$p holds the initial parameters: when events change them it gets
  # its own cotangent (dp_init), which belongs to p directly
  pinit <- data$prob$p
  hh <- data$param_history
  split <- length(p0) && !is.null(pinit) &&
    ((!is.null(hh) && nrow(hh$p) > 1L) || !identical(as.numeric(pinit), as.numeric(p0)))
  ev <- function(U, p, pi = NULL) {
    d <- data
    if (!is.null(pi)) d$prob$p <- pi
    if (length(irows) && !anyNA(irows)) d$sens_interp <- list(t = itimes, U = U[irows, , drop = FALSE])
    d$U <- U[!hidden, , drop = FALSE]
    d$t <- env$t[!hidden]
    d$p_final <- p
    if (!identical(p, p0) && !is.null(d$param_history)) {
      d$param_history$p[nrow(d$param_history$p), ] <- p
    }
    if (!split && !identical(p, p0) && !is.null(d$prob)) d$prob$p <- p
    v <- loss(.new_solution(cls, d))
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) stop("the loss must return a finite number", call. = FALSE)
    v
  }
  value <- ev(U0, p0)
  dU <- U0
  dU[] <- 0
  lg <- .mtk$sens_loss_grad
  if (!is.null(lg)) {
    d <- data
    if (length(irows) && !anyNA(irows)) d$sens_interp <- list(t = itimes, U = U0[irows, , drop = FALSE])
    d$U <- U0[!hidden, , drop = FALSE]; d$t <- env$t[!hidden]
    g <- lg(.new_solution(cls, d))
    gu <- if (is.list(g)) g$u else g
    gu <- as.matrix(gu)
    # the cotangent is in the view of the loss, the saved components
    cols <- env$save_idxs %||% seq_len(ncol(d$U))
    if (!identical(dim(gu), c(nrow(d$U), length(cols)))) stop("`loss_grad` must return a ", nrow(d$U), " by ",
                                                             length(cols), " matrix, as as.matrix(sol)", call. = FALSE)
    dU[!hidden, cols] <- gu
    gp <- if (is.list(g) && !is.null(g$p)) as.numeric(g$p) else numeric(length(p0))
    if (length(gp) != length(p0)) stop("the `p` of `loss_grad` must have ", length(p0), " values", call. = FALSE)
  }
  # Every saved value is differentiated unless the caller declares which
  # saved rows the loss reads (loss_rows of solve_gradient(), rows of the
  # visible solution); a black box probe cannot prove that a loss ignores
  # a value.
  # loss_grad covers the saved values: the rows saved only for sol(t) are
  # still differentiated
  rows <- if (!is.null(lg)) which(hidden) else seq_len(nrow(U0))
  lr <- .mtk$sens_loss_rows
  if (is.null(lg) && !is.null(lr)) {
    vis <- which(!hidden)
    lr <- if (identical(lr, "last")) length(vis) else as.integer(lr)
    if (anyNA(lr) || any(lr < 1L | lr > length(vis)))
      stop("`loss_rows` must index the ", length(vis), " saved rows", call. = FALSE)
    # the rows saved for interpolation only are read through sol(t)
    rows <- c(vis[lr], which(hidden))
  }
  scols <- env$save_idxs %||% seq_len(ncol(U0))
  for (k in which(row(U0) %in% rows & col(U0) %in% scols)) {
    h <- .sens_h(U0[k])
    Up <- U0; Um <- U0
    Up[k] <- U0[k] + h
    Um[k] <- U0[k] - h
    dU[k] <- (ev(Up, p0) - ev(Um, p0)) / (2 * h)
  }
  dp <- if (!is.null(lg)) gp else if (length(p0)) .sens_fd_grad(function(q) ev(U0, .sens_pvec(p0, q)), as.numeric(p0)) else numeric()
  dp_init <- if (split) .sens_fd_grad(function(q) ev(U0, p0, .sens_pvec(pinit, q)), as.numeric(pinit))
  if (!is.null(dp_init) && all(dp_init == 0)) dp_init <- NULL
  # the cotangent of the saved times (dgdt_discrete): it only matters for
  # samples that move with an event, so dT(rows) differentiates the loss in
  # the times of those rows only, when a sensitivity method asks
  tvis <- env$t[!hidden]
  pidx <- if (!is.null(data$param_history)) .sol_params_idx(list(param_history = data$param_history), tvis, TRUE)
  ev_t <- function(tv) {
    d <- data
    d$sens_param_idx <- pidx
    if (length(irows) && !anyNA(irows)) d$sens_interp <- list(t = itimes, U = U0[irows, , drop = FALSE])
    d$U <- U0[!hidden, , drop = FALSE]
    d$t <- tv
    loss(.new_solution(cls, d))
  }
  vis <- which(!hidden)
  dT <- function(rows) {
    vapply(rows, function(r) {
      j <- match(r, vis)
      if (is.na(j)) return(0)
      hh <- 1e-6 * max(1, abs(tvis[j]))
      tp <- tvis; tm <- tvis
      tp[j] <- tp[j] + hh; tm[j] <- tm[j] - hh
      (ev_t(tp) - ev_t(tm)) / (2 * hh)
    }, 0)
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
  list(value = value, dU = unname(dU), dp = dp, dp_init = dp_init, dp_hist = dp_hist, t_hist = h$t, dT = dT)
}

# A fully implicit DAE F(du, u, p, t) = 0 of index 1 is the semi-explicit
# DAE u_d' = w, 0 = F(w, u, p, t) in (u, w), with w the derivatives of the
# differential variables: a mass matrix ODEProblem, whose sensitivities
# (with events, reinitialization of the algebraic variables and every
# method) apply. The differential equation of variable d sits in row d, so
# the mass matrix is diagonal.
.sens_dae_gradient <- function(loss, prob, alg, dots, sensealg) {
  F <- prob$f$f %||% prob$f
  u0 <- as.numeric(prob$u0)
  du0 <- as.numeric(prob$du0 %||% numeric(length(u0)))
  n <- length(u0)
  dv <- prob$differential_vars
  if (is.null(dv)) {
    Jd <- .sens_fd_jac(function(x) as.numeric(F(x, u0, prob$p, prob$tspan[1L])), du0)
    dv <- colSums(abs(Jd)) > 0
  }
  d <- which(as.logical(dv))
  nd <- length(d)
  rest <- setdiff(seq_len(n + nd), d)
  f_aug <- function(z, p, t) {
    du <- numeric(n)
    du[d] <- z[n + seq_len(nd)]
    out <- numeric(n + nd)
    out[d] <- z[n + seq_len(nd)]
    out[rest] <- as.numeric(F(du, z[seq_len(n)], p, t))
    out
  }
  M <- matrix(0, n + nd, n + nd)
  M[cbind(d, d)] <- 1
  aprob <- ODEProblem(ODEFunction(f_aug, mass_matrix = M), c(u0, du0[d]), prob$tspan, prob$p)
  # the callbacks see the DAE variables u of the state (u, w); after an
  # affect w is kept as the start of the algebraic re-solve
  lift <- function(cbs) {
    if (is.null(cbs)) return(NULL)
    cs <- .as_callback_set(cbs)
    wc <- function(cond) if (is.function(cond)) function(u, t, integ, ...) cond(u[seq_len(n)], t, integ, ...)
    wa <- function(aff) if (is.function(aff)) function(integ, ...) {
      z <- as.numeric(integ$u)
      integ$u <- z[seq_len(n)]
      on.exit(integ$u <- c(as.numeric(integ$u), z[n + seq_len(nd)]))
      aff(integ, ...)
    }
    lc <- function(cb, cont) {
      cb$condition <- wc(cb$condition)
      cb$affect <- wa(cb$affect)
      if (cont) cb$affect_neg <- wa(cb$affect_neg)
      .cb_rebind_init(cb)
    }
    structure(list(continuous = lapply(cs$continuous, lc, TRUE), discrete = lapply(cs$discrete, lc, FALSE)),
              class = "CallbackSet")
  }
  aprob$kwargs <- prob$kwargs
  aprob$kwargs$callback <- lift(prob$kwargs$callback)
  dots$callback <- lift(dots$callback)
  # a mass matrix method for the DAE method
  oalg <- if (inherits(alg, "DFBDF")) FBDF() else if (inherits(alg, "DImplicitEuler")) ImplicitEuler() else
    if (inherits(alg, "DABDF2")) ABDF2() else if (is.null(alg) || isTRUE(alg$dae) || !isTRUE(alg$mass_matrix))
      Rodas5P() else alg
  # the solution in the DAE variables, with its interpolation (the steps
  # interpolate (u, w); .interp_states keeps the first n values) and the
  # samples solve_gradient() saved for an interpolating loss
  restrict <- function(s) {
    e <- .sol_env(s)
    data <- mget(setdiff(ls(e, all.names = TRUE), "obs_cache"), envir = e)
    data$U <- e$U[, seq_len(n), drop = FALSE]
    if (!is.null(data$sens_interp)) data$sens_interp$U <- data$sens_interp$U[, seq_len(n), drop = FALSE]
    data$prob <- prob
    data$prob$p <- e$prob$p
    data$alg <- alg %||% e$alg
    .new_solution(setdiff(class(s), c("SciMLSolution", "function")), data)
  }
  g <- .solve_gradient(function(s) loss(restrict(s)), aprob, oalg, dots, sensealg)
  list(u0 = stats::setNames(g$u0[seq_len(n)], names(prob$u0)), p = g$p, value = g$value, sol = restrict(g$sol))
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
                           retcode = aenv$retcode, stats = aenv$stats, p = aenv$p_final %||% prob$p)
  # the parameters in force at each saved time, so the loss sees them and
  # .sens_no_hist() can refuse a loss that reads earlier ones
  if (!is.null(aenv$param_history)) {
    senv <- .sol_env(sol)
    senv$param_history <- aenv$param_history
  }
  ct <- .sens_loss_cotangent(loss, sol)
  .sens_no_hist(ct)
  grad <- numeric(m)
  for (r in seq_len(nrow(U))) grad <- grad + as.numeric(ct$dU[r, ] %*% matrix(U[r, n + seq_len(n * m)], n, m))
  if (np) {
    Qend <- Qof(U[nrow(U), ])
    # an event after the last saved time changes the final parameters
    if (has_cb && length(rec$events)) {
      tl <- aenv$t[nrow(U)]
      dir <- sign(prob$tspan[2L] - prob$tspan[1L])
      late <- which(vapply(rec$events, function(ev) dir * (ev$t - tl) > 1e-12 * max(1, abs(tl)), TRUE))
      if (length(late)) Qend <- rec$events[[late[length(late)]]]$Q_right
    }
    grad <- grad + as.numeric(ct$dp %*% Qend)
  }
  if (!is.null(ct$dp_init)) grad[n + seq_len(np)] <- grad[n + seq_len(np)] + ct$dp_init
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
