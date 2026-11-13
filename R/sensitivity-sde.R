# Adjoint sensitivities of stochastic and random differential equations,
# after SciMLSensitivity.jl v7.119 src/backsolve_adjoint.jl and
# src/interpolating_adjoint.jl (SDEAdjointProblem, RODEAdjointProblem),
# src/derivative_wrappers.jl (vecjacobian!, jacNoise!: the diffusion of the
# adjoint SDE), src/sde_tools.jl (StochasticTransformedFunction) and
# src/sensitivity_interface.jl (the reverse solve on the forward grid).
#
# The adjoint of du = f dt + g dW is the backward SDE
#   dlambda = -(f_u' lambda + g_u) dt - (sigma_u' lambda) o dW,
#   dgrad   = -(f_p' lambda + g_p) dt - (sigma_p' lambda) o dW
# driven by the forward noise path, with the forward state either read from
# the saved solution (InterpolatingAdjoint) or solved backward with it
# (BacksolveAdjoint, whose drift is the transformed drift f - sigma_u' sigma
# of StochasticTransformedFunction when the forward solve was Ito; Julia's
# InterpolatingAdjoint keeps f for an Ito solve as well). The diffusion has one
# column per Wiener process: for diagonal noise column i holds
# (sigma_u' lambda)_i on the diagonal (lambda_i grad sigma_i with
# `noisemixing`), lambda_i d sigma_i / dp and sigma_i(y); for scalar noise it
# is a vector; for a noise matrix column j holds the products with column j.
#
# Julia integrates the adjoint backward in time with the reversed noise
# reverse(sol.W). R's SDE and RODE solvers integrate forward in time, so the
# adjoint problem is posed in the reversed time tau = t0 + T - t, which runs
# from t0 to T: its drift is the negated backward drift, its diffusion is
# unchanged and its noise is the forward path read backward, the same steps
# the backward solve takes.
#
# Julia's adjoint diffusion also subtracts the continuous cost gradient (its
# noise term function calls accumulate_cost! like the drift), which does not
# belong to the adjoint SDE; here the cost enters the drift only. Julia's
# InterpolatingAdjoint samples a discrete cost with the state its last drift
# evaluation interpolated (the stale buffer S.y of ReverseLossCallback), the
# state one step later for one stage methods such as EM and RandomEM; here
# the cost reads the forward state at its own time. The discrete cost of
# BacksolveAdjoint reads the reset state at a checkpoint time and the
# backsolved state at a cost time above it, in the order of Julia's
# callbacks.

# Central difference Jacobian with step `rel * max(1, |x|)`; the nested
# differences of the transformed drift use the fourth root of the machine
# epsilon, which balances the truncation and rounding errors of both levels.
.sens_fd_jac_step <- function(F, x, rel) {
  x <- as.numeric(x)
  cols <- lapply(seq_along(x), function(j) {
    h <- rel * max(1, abs(x[j]))
    xp <- x; xm <- x
    xp[j] <- x[j] + h
    xm[j] <- x[j] - h
    (as.numeric(F(xp)) - as.numeric(F(xm))) / (2 * h)
  })
  if (!length(cols)) return(matrix(0, length(F(x)), 0L))
  do.call(cbind, cols)
}

# The derivative bundle of a stochastic problem: the drift `a` the adjoint
# differentiates, its Jacobians, and the noise: G(u, p, t) (a vector for
# diagonal and scalar noise, a matrix otherwise) with the Jacobians of
# vec(G) by u and p. For an RODE the functions take the noise value W.
.sens_sde_fns <- function(prob, transform, corfunc_analytical, diffusion_jac, diffusion_paramjac) {
  rode <- inherits(prob, "RODEProblem")
  n <- length(prob$u0)
  np <- length(prob$p)
  f <- prob$f$f
  out <- list(n = n, np = np, rode = rode)
  if (rode) {
    out$a <- function(u, p, t, W) as.numeric(f(u, p, t, W))
    out$Ja <- function(u, p, t, W) .sens_fd_jac(function(z) f(z, p, t, W), u)
    out$Pa <- function(u, p, t, W) if (!np) matrix(0, n, 0L) else
      .sens_fd_jac(function(q) f(u, .sens_pvec(p, q), t, W), as.numeric(p))
    return(out)
  }
  g <- prob$f$g
  if (inherits(g, "SDEFunction")) g <- g$g
  out$type <- prob$noise_type %||% prob$f$noise_type %||% "diagonal"
  out$G <- function(u, p, t) {
    v <- g(u, p, t)
    if (out$type == "general") as.matrix(v) else as.numeric(v)
  }
  out$m <- switch(out$type, scalar = 1L, diagonal = n, general = ncol(out$G(prob$u0, prob$p, prob$tspan[1L])))
  out$Gu <- if (is.function(diffusion_jac)) function(u, p, t) as.matrix(diffusion_jac(u, p, t)) else
    function(u, p, t) .sens_fd_jac(function(z) as.numeric(g(z, p, t)), u)
  out$Gp <- if (is.function(diffusion_paramjac)) function(u, p, t) as.matrix(diffusion_paramjac(u, p, t)) else
    function(u, p, t) if (!np) matrix(0, length(as.numeric(g(u, p, t))), 0L) else
      .sens_fd_jac(function(q) as.numeric(g(u, .sens_pvec(p, q), t)), as.numeric(p))
  ujac <- prob$f$jac
  if (!transform) {
    out$a <- function(u, p, t, W = NULL) as.numeric(f(u, p, t))
    out$Ja <- if (is.function(ujac)) function(u, p, t, W = NULL) as.matrix(ujac(u, p, t)) else
      function(u, p, t, W = NULL) .sens_fd_jac(function(z) f(z, p, t), u)
    out$Pa <- function(u, p, t, W = NULL) if (!np) matrix(0, n, 0L) else
      .sens_fd_jac(function(q) f(u, .sens_pvec(p, q), t), as.numeric(p))
    return(out)
  }
  # the Ito drift transformed for the Stratonovich backward solve
  rel <- .Machine$double.eps^(1 / 4)
  cor <- if (is.function(corfunc_analytical)) function(u, p, t) as.numeric(corfunc_analytical(u, p, t)) else
    function(u, p, t) as.numeric(crossprod(.sens_fd_jac_step(function(z) as.numeric(g(z, p, t)), u, rel),
                                           as.numeric(g(u, p, t))))
  h <- if (is.function(corfunc_analytical)) .sens_fd_jac else function(F, x) .sens_fd_jac_step(F, x, rel)
  out$a <- function(u, p, t, W = NULL) as.numeric(f(u, p, t)) - cor(u, p, t)
  out$Ja <- function(u, p, t, W = NULL) h(function(z) out$a(z, p, t), u)
  out$Pa <- function(u, p, t, W = NULL) if (!np) matrix(0, n, 0L) else
    h(function(q) out$a(u, .sens_pvec(p, q), t), as.numeric(p))
  out
}

# The diffusion of the adjoint for the state rows (lambda, grad[, y]): a
# vector for scalar noise, a matrix with a column per Wiener process
# otherwise (vecjacobian! and jacNoise! of SciMLSensitivity).
.sens_sde_diffusion <- function(fns, l, y, p, t, backsolve, mixing) {
  n <- fns$n
  np <- fns$np
  G <- fns$G(y, p, t)
  Gu <- fns$Gu(y, p, t)
  Gp <- fns$Gp(y, p, t)
  if (fns$type == "scalar")
    return(c(-as.numeric(crossprod(Gu, l)), -as.numeric(crossprod(Gp, l)), if (backsolve) G))
  m <- fns$m
  rows <- n + np + if (backsolve) n else 0L
  out <- matrix(0, rows, m)
  if (fns$type == "diagonal") {
    if (mixing) out[seq_len(n), ] <- -t(Gu * l) else diag(out)[seq_len(n)] <- -as.numeric(crossprod(Gu, l))
    if (np) out[n + seq_len(np), ] <- -t(Gp * l)
    if (backsolve) out[cbind(n + np + seq_len(n), seq_len(n))] <- G
    return(out)
  }
  for (j in seq_len(m)) {
    rj <- (j - 1L) * n + seq_len(n)
    out[seq_len(n), j] <- -as.numeric(crossprod(Gu[rj, , drop = FALSE], l))
    if (np) out[n + seq_len(np), j] <- -as.numeric(crossprod(Gp[rj, , drop = FALSE], l))
    if (backsolve) out[n + np + seq_len(n), j] <- G[, j]
  }
  out
}

#' Adjoint problems of SDEs and RODEs
#'
#' `SDEAdjointProblem(sol, sensealg, alg, ...)` builds the adjoint of the
#' solution of an [SDEProblem()] `du = f dt + g dW`, the backward SDE
#' `dlambda = -(f_u' lambda + g_u) dt - (g_u' lambda) o dW` with the
#' parameter gradient `dgrad = -(f_p' lambda + g_p) dt - (g_p' lambda) o dW`,
#' driven by the noise path of the forward solution.
#' `RODEAdjointProblem()` builds the adjoint of an [RODEProblem()] `du = f(u,
#' p, t, W) dt`, the random ODE `dlambda = -(f_u' lambda + g_u) dt` with the
#' forward noise values. [InterpolatingAdjoint()] reads the forward state
#' from the saved solution; [BacksolveAdjoint()] solves it backward with
#' the adjoint (with `checkpointing`, resetting it to the saved values at the
#' saved times) and, for an Ito forward solve, uses the transformed drift
#' `f - sigma_u' sigma` of [StochasticTransformedFunction()] (with
#' `corfunc_analytical` as the correction when given), as Julia does.
#'
#' Julia integrates the adjoint backward in time with the reversed noise.
#' R's SDE and RODE solvers integrate forward in time, so the returned
#' problem is posed in the reversed time `tau = t0 + T - t`, which runs from
#' `t0` to `T`, with the negated drift and the forward noise read backward:
#' the steps are those of the backward solve. Solve it with fixed steps on
#' the forward grid (`dt` and `tstops` of the problem, which
#' [adjoint_sensitivities()] passes). The derivatives are central
#' differences (Julia: automatic differentiation); `diffusion_jac(u, p, t)`
#' and `diffusion_paramjac(u, p, t)` give the Jacobians of `vec(g)`. The
#' continuous cost enters the drift only (Julia's diffusion also subtracts
#' its gradient), and [InterpolatingAdjoint()] samples a discrete cost with
#' the forward state at its time (Julia: the state of the last drift
#' evaluation, one step later for [EM()] and [RandomEM()]). Events of the
#' forward solve and mass matrices are not supported.
#'
#' @inheritParams adjoint_sensitivities
#' @param sol The solution of an `SDEProblem` or an `RODEProblem`.
#' @param sensealg [InterpolatingAdjoint()] or [BacksolveAdjoint()];
#'   `noisemixing = TRUE` uses the full Jacobian of diagonal noise.
#' @param alg The algorithm of the adjoint solve; the forward algorithm when
#'   `NULL`.
#' @param checkpoints Unused (the saved times are the checkpoints).
#' @param callback Must be `NULL`.
#' @param corfunc_analytical The Ito correction `sigma_u' sigma` as a
#'   function `(u, p, t)`.
#' @param diffusion_jac,diffusion_paramjac Jacobians of `vec(g(u, p, t))` by
#'   `u` and by `p`.
#' @return An `SDEProblem` or `RODEProblem` in the reversed time, with
#'   elements `alg`, `dt`, `tstops` and `finish` (a function of its solution
#'   returning `list(du0, dp)`).
#' @examples
#' f <- function(u, p, t) p[1] * u
#' g <- function(u, p, t) p[2] * u
#' W <- c(0, 0.3, 0.1, -0.2, 0.25)
#' prob <- SDEProblem(f, g, 1, c(0, 1), c(1, 0.5), noise = NoiseGrid(seq(0, 1, 0.25), W))
#' sol <- solve(prob, EulerHeun(), dt = 0.25)
#' # d u(1) / d(u0, p) along the path
#' adjoint_sensitivities(sol, EulerHeun(), t = 1, dgdu_discrete = function(u, p, t, i) 1,
#'                       sensealg = BacksolveAdjoint())
#' @export
SDEAdjointProblem <- function(sol, sensealg = InterpolatingAdjoint(), alg = NULL, t = NULL, dgdu_discrete = NULL,
                              dgdp_discrete = NULL, dgdu_continuous = NULL, dgdp_continuous = NULL, g = NULL,
                              checkpoints = NULL, callback = NULL, corfunc_analytical = NULL, diffusion_jac = NULL,
                              diffusion_paramjac = NULL, no_start = FALSE) {
  .sens_stochastic_adjoint(sol, "SDEProblem", sensealg, alg, t, dgdu_discrete, dgdp_discrete, dgdu_continuous,
                           dgdp_continuous, g, callback, corfunc_analytical, diffusion_jac, diffusion_paramjac,
                           no_start)
}

#' @rdname SDEAdjointProblem
#' @export
RODEAdjointProblem <- function(sol, sensealg = InterpolatingAdjoint(), alg = NULL, t = NULL, dgdu_discrete = NULL,
                               dgdp_discrete = NULL, dgdu_continuous = NULL, dgdp_continuous = NULL, g = NULL,
                               checkpoints = NULL, callback = NULL, no_start = FALSE) {
  .sens_stochastic_adjoint(sol, "RODEProblem", sensealg, alg, t, dgdu_discrete, dgdp_discrete, dgdu_continuous,
                           dgdp_continuous, g, callback, NULL, NULL, NULL, no_start)
}

.sens_stochastic_adjoint <- function(sol, kind, sensealg, alg, t, dgdu_discrete, dgdp_discrete, dgdu_continuous,
                                     dgdp_continuous, g, callback, corfunc_analytical, diffusion_jac,
                                     diffusion_paramjac, no_start) {
  env <- .sol_env(sol)
  prob <- env$prob
  if (!inherits(prob, kind)) stop(sprintf("%s needs the solution of an %s",
                                          if (kind == "SDEProblem") "SDEAdjointProblem()" else "RODEAdjointProblem()",
                                          kind), call. = FALSE)
  if (is.null(dgdu_discrete) && is.null(dgdu_continuous) && is.null(g))
    stop("Either `dgdu_discrete`, `dgdu_continuous`, or `g` must be specified.", call. = FALSE)
  if (!is.null(t) && is.null(dgdu_discrete) && is.null(dgdp_discrete))
    stop("It looks like you're using the direct `adjoint_sensitivities` interface with a discrete cost ",
         "function but no specified `dgdu_discrete` or `dgdp_discrete`.", call. = FALSE)
  method <- class(sensealg)[1L]
  if (!method %in% c("BacksolveAdjoint", "InterpolatingAdjoint"))
    stop("the adjoints of ", kind, "s are BacksolveAdjoint() and InterpolatingAdjoint()", call. = FALSE)
  if (!length(prob$p))
    stop("Your model does not have parameters, and thus it is impossible to calculate the derivative of the ",
         "solution with respect to the parameters. Your model must have parameters to use parameter ",
         "sensitivity calculations!", call. = FALSE)
  cbs <- .as_callback_set(callback)
  fcb <- .as_callback_set(prob$kwargs$callback)
  if (length(cbs$continuous) + length(cbs$discrete) + length(fcb$continuous) + length(fcb$discrete))
    stop("events are not supported by the adjoints of ", kind, "s", call. = FALSE)
  M <- prob$f$mass_matrix
  if (is.matrix(M) && !.is_identity_mass(M))
    stop("mass matrices are not supported by the adjoints of ", kind, "s", call. = FALSE)
  if (is.null(env$W) || !length(env$W$t))
    stop("the adjoint needs the noise path of the forward solution", call. = FALSE)
  backsolve <- method == "BacksolveAdjoint"
  rode <- kind == "RODEProblem"
  fwd_alg <- env$alg
  transform <- !rode && backsolve && !identical(fwd_alg$interpretation, "Stratonovich")
  fns <- .sens_sde_fns(prob, transform, corfunc_analytical, diffusion_jac, diffusion_paramjac)
  n <- fns$n
  np <- fns$np
  p <- prob$p
  t0 <- prob$tspan[1L]
  tend <- if (identical(env$retcode, "Terminated")) env$t[length(env$t)] else prob$tspan[2L]
  # tau = s - t runs from t0 to tend
  s <- t0 + tend
  tol <- function(x) 1e-12 * max(1, abs(x))
  inside <- function(x) x - t0 > tol(t0) & tend - x > tol(tend)
  tr <- if (!backsolve) .sens_trajectory(sol)

  cont <- !is.null(g) || !is.null(dgdu_continuous) || !is.null(dgdp_continuous)
  gu_c <- dgdu_continuous %||% if (!is.null(g)) function(u, p, tt) .sens_fd_grad(function(z) g(z, p, tt), u)
  gp_c <- dgdp_continuous %||% if (!is.null(g) && is.null(dgdu_continuous))
    function(u, p, tt) .sens_fd_grad(function(q) g(u, .sens_pvec(p, q), tt), as.numeric(p))
  ts <- if (is.null(t)) numeric() else as.numeric(t)
  ord <- order(ts)
  cost <- list(t = ts[ord], i = ord)

  # z = (lambda, grad[, y])
  il <- seq_len(n)
  ig <- n + seq_len(np)
  iy <- if (backsolve) n + np + seq_len(n) else integer()
  state_y <- function(z, tt) if (backsolve) z[iy] else tr$at(tt)
  ftimes <- env$t
  ck_t <- if (backsolve && isTRUE(sensealg$checkpointing)) unique(ftimes[inside(ftimes)]) else numeric()

  # the drift in tau: the negated backward drift
  drift_real <- function(z, p, tt, W = NULL) {
    l <- z[il]
    y <- state_y(z, tt)
    dl <- -as.numeric(crossprod(fns$Ja(y, p, tt, W), l))
    dg <- -as.numeric(crossprod(fns$Pa(y, p, tt, W), l))
    if (cont) {
      if (!is.null(gu_c)) dl <- dl - as.numeric(gu_c(y, p, tt))
      if (!is.null(gp_c)) dg <- dg - as.numeric(gp_c(y, p, tt))
    }
    c(dl, dg, if (backsolve) fns$a(y, p, tt, W))
  }
  drift <- if (rode) function(z, p, tau, W) -drift_real(z, p, s - tau, W) else
    function(z, p, tau) -drift_real(z, p, s - tau)
  drift <- structure(drift, class = c(.sens_function_class(method), SensitivityFunction, "function"))

  # the discrete costs and the checkpoint reset at one time. Julia's
  # backsolve_checkpoint_callbacks come before the ReverseLossCallback at the
  # same time (the cost reads the reset state), and a cost time above the
  # checkpoint time (0.5 against the saved 0.49999999999999994) is reached
  # first (the cost reads the backsolved state).
  ck_all <- if (backsolve && isTRUE(sensealg$checkpointing)) unique(ftimes) else numeric()
  point <- function(tt, z, p) {
    ks <- rev(which(abs(cost$t - tt) <= tol(tt)))
    ck <- ck_all[abs(ck_all - tt) <= tol(tt)]
    apply_cost <- function(z, ks) {
      for (k in ks) {
        y <- state_y(z, tt)
        i <- cost$i[k]
        if (!is.null(dgdu_discrete)) z[il] <- z[il] + as.numeric(dgdu_discrete(y, p, tt, i))
        if (!is.null(dgdp_discrete)) z[ig] <- z[ig] + as.numeric(dgdp_discrete(y, p, tt, i))
      }
      z
    }
    if (!length(ck)) return(apply_cost(z, ks))
    z <- apply_cost(z, ks[cost$t[ks] > ck[1L]])
    k <- which(ftimes == ck[1L])
    z[iy] <- env$U[k[length(k)], ]
    apply_cost(z, ks[cost$t[ks] <= ck[1L]])
  }
  z0 <- numeric(n + np + length(iy))
  if (backsolve) z0[iy] <- if (abs(ftimes[length(ftimes)] - tend) <= tol(tend)) env$U[nrow(env$U), ] else
    .interp_states(env, tend)[1L, ]
  z0 <- point(tend, z0, p)
  inner <- sort(unique(c(cost$t, ck_t)))
  inner <- inner[inside(inner)]
  cb_rev <- if (length(inner)) {
    tau_inner <- sort(s - inner)
    cb <- DiscreteCallback(function(u, tau, integ) .near_any(tau_inner, tau), function(integ) {
      integ$u <- point(s - integ$t, integ$u, integ$p)
    }, save_positions = c(FALSE, FALSE))
    cb$tstops <- tau_inner
    cb
  }
  # the forward noise read backward: W~(tau) = W(s - tau)
  wt <- env$W$t
  W <- as.matrix(env$W$W)
  keep <- !duplicated(wt, fromLast = TRUE)
  wt <- wt[keep]
  W <- W[keep, , drop = FALSE]
  noise <- NoiseGrid(s - rev(wt), W[rev(seq_along(wt)), , drop = FALSE])
  tstops <- sort(unique(s - ftimes[inside(ftimes)]))
  fdt <- abs(diff(ftimes))
  dt <- if (length(fdt)) fdt[fdt > 0][1L] else abs(tend - t0)
  if (rode) {
    aprob <- RODEProblem(RODEFunction(drift), z0, c(t0, tend), p, noise = noise, callback = CallbackSet(cb_rev))
  } else {
    mixing <- isTRUE(sensealg$noisemixing)
    diffusion <- function(z, p, tau) {
      tt <- s - tau
      .sens_sde_diffusion(fns, z[il], state_y(z, tt), p, tt, backsolve, mixing)
    }
    diffusion <- structure(diffusion, class = c(.sens_function_class(method), SensitivityFunction, "function"))
    nrp <- if (fns$type == "scalar") NULL else matrix(0, length(z0), fns$m)
    aprob <- SDEProblem(drift, diffusion, z0, c(t0, tend), p, noise_rate_prototype = nrp, noise = noise,
                        callback = CallbackSet(cb_rev))
  }
  aprob$alg <- alg %||% fwd_alg
  aprob$dt <- dt
  aprob$tstops <- tstops
  aprob$time_map <- function(tau) s - tau
  aprob$finish <- function(adj) {
    aenv <- .sol_env(adj)
    zf <- aenv$U[nrow(aenv$U), ]
    if (!no_start) zf <- point(t0, zf, p)
    du0 <- zf[il]
    dp <- zf[ig]
    names(du0) <- names(prob$u0)
    names(dp) <- names(prob$p)
    list(du0 = du0, dp = dp)
  }
  aprob
}

# solve_gradient() with BacksolveAdjoint or InterpolatingAdjoint for an SDE
# or RODE (concrete_solve.jl): the loss cotangent of every saved state is a
# discrete cost of the adjoint, solved on the grid of a forward solve that
# saves every step (the same path) with the forward algorithm.
.sens_stochastic_gradient <- function(loss, prob, alg, dots, sensealg) {
  sol <- do.call(solve, c(list(prob, alg), dots))
  env <- .sol_env(sol)
  ct <- .sens_loss_cotangent(loss, sol)
  full <- dots[setdiff(names(dots), c("saveat", "save_everystep", "save_start", "save_end", "save_idxs", "dense"))]
  fsol <- if (length(full) == length(dots)) sol else do.call(solve, c(list(prob, alg), full))
  np <- length(prob$p)
  last <- length(env$t)
  at_start <- isTRUE(abs(env$t[1L] - prob$tspan[1L]) <= 1e-14 * max(1, abs(prob$tspan[1L])))
  dgdp <- if (np && any(ct$dp != 0)) function(u, p, tt, i) if (i == last) ct$dp else numeric(np)
  res <- do.call(adjoint_sensitivities,
                 c(list(fsol, env$alg, t = env$t, dgdu_discrete = function(u, p, tt, i) ct$dU[i, ],
                        dgdp_discrete = dgdp, sensealg = sensealg, no_start = !at_start),
                   dots[intersect(names(dots), "dt")]))
  list(u0 = stats::setNames(res$du0, names(prob$u0)), p = stats::setNames(res$dp, names(prob$p)),
       value = ct$value, sol = sol)
}
