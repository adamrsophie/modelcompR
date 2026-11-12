# Tau leaping, after JumpProcesses.jl: jumps.jl (RegularJump),
# aggregators/aggregators.jl (PureLeaping) and simple_regular_solve.jl
# (SimpleTauLeaping, SimpleExplicitTauLeaping, SimpleImplicitTauLeaping,
# SimpleTrapezoidalLeaping, SimpleAdaptiveTauLeaping, _process_saveat,
# compute_hor, precompute_reaction_conditions, compute_gi, compute_tau,
# implicit_equation!, solve_implicit, compute_tau_implicit,
# compute_drift_jacobian, is_stiff and the stepping loops).
#
# A leaping step fires each jump a Poisson number of times with mean
# rate * tau. Mass action jumps carry their stoichiometry (maj_rs, maj_ns
# and maj_ratefun of the numeric jump set); a RegularJump gives the rates of
# all its jumps at once and applies the counts with its function c.

#' Regular jumps and tau leaping
#'
#' `RegularJump(rate, c, numjumps)` encodes `numjumps` jump processes for
#' tau leaping, as in JumpProcesses.jl: `rate(u, p, t)` returns the vector of
#' their rates and `c(u, p, t, counts, mark)` returns the change of the state
#' when jump `i` fires `counts[i]` times. `JumpProblem(prob, PureLeaping(),
#' jumps)`, with `prob` a [DiscreteProblem()] and `jumps` a `RegularJump` or
#' mass action jumps, builds a pure leaping problem.
#'
#' Solvers:
#' * `SimpleTauLeaping()`: fixed step tau leaping for a `RegularJump`;
#'   `solve()` needs `dt`.
#' * `SimpleExplicitTauLeaping(epsilon)`: adaptive explicit tau leaping for
#'   mass action jumps (or a `RegularJump` together with mass action jumps
#'   giving the stoichiometry), with the step size of Cao, Gillespie and
#'   Petzold (2006); a step that would make a population negative is halved.
#' * `SimpleImplicitTauLeaping(epsilon)` and `SimpleTrapezoidalLeaping(epsilon)`:
#'   the implicit (Rathinam et al. 2003) and trapezoidal tau leaping methods
#'   for stiff mass action systems.
#' * `SimpleAdaptiveTauLeaping()`: switches between the explicit step and
#'   `implicit_alg` according to a stiffness test (Cao et al. 2007).
#'
#' The tau leaping methods of StochasticDiffEq.jl take a `RegularJump`:
#' * `TauLeaping()`: adaptive tau leaping; a step is accepted when the rates
#'   change by less than `50 reltol` relative to their value (StochasticDiffEq's
#'   error estimate), otherwise it is retried with a shorter step whose counts
#'   are a binomial thinning of the rejected ones (the Poisson bridge).
#' * `CaoTauLeaping()`: each step size is chosen a posteriori by the formula of
#'   Cao, Gillespie and Petzold (2006), with the net changes of the jumps read
#'   from `c`; steps are never rejected.
#' * `ImplicitTauLeaping()` and `ThetaTrapezoidalTauLeaping(theta)`: the
#'   implicit and theta trapezoidal methods (Rathinam et al. 2003; Hu, Li and
#'   Min 2011) with step `dt`: `X1 = X0 + nu k + theta dt (drift(X1) -
#'   drift(X0))`, `drift(u) = c(u, p, t, rate(u, p, t), NULL)`.
#'
#' `solve()` takes `saveat`, `save_start`, `save_end`, `dtmin` (default
#' `1e-10`), `seed`, and for the StochasticDiffEq methods `dt`, `adaptive`,
#' `reltol` (default `1e-2`), `dtmax` and `gamma` (default `0.9`). Solutions
#' are piecewise constant.
#'
#' @param rate Function `(u, p, t)` returning the rates of all jumps.
#' @param c Function `(u, p, t, counts, mark)` returning the state change.
#' @param numjumps Number of jumps.
#' @param mark_dist Not supported (as in JumpProcesses).
#' @param epsilon Error control parameter of the step size selection.
#' @param implicit_alg Algorithm of the implicit steps of
#'   `SimpleAdaptiveTauLeaping()`.
#' @param eigenvalue_check Judge stiffness from the eigenvalues of the drift
#'   Jacobian instead of the spread of the rates.
#' @param stiffness_ratio_threshold Eigenvalue ratio above which the system
#'   counts as stiff.
#' @param implicit_epsilon_factor Factor relaxing `epsilon` on implicit steps.
#' @param theta Implicitness of `ThetaTrapezoidalTauLeaping()`.
#' @return A `RegularJump`, an aggregator or an algorithm object.
#' @examples
#' rj <- RegularJump(function(u, p, t) 0.1 * u[1],
#'                   function(u, p, t, counts, mark) -counts[1], 1)
#' jprob <- JumpProblem(DiscreteProblem(20, c(0, 2)), PureLeaping(), rj)
#' sol <- solve(jprob, SimpleTauLeaping(), dt = 0.1, seed = 1)
#' tail(sol$u, 1)
#' maj <- MassActionJump(c(1, 1), list(list(c(1, 1)), list(c(2, 1))),
#'                       list(list(c(1, -1), c(2, 1)), list(c(1, 1), c(2, -1))))
#' jp2 <- JumpProblem(DiscreteProblem(c(100, 100), c(0, 10)), PureLeaping(), maj)
#' tail(solve(jp2, SimpleImplicitTauLeaping(), seed = 1)$u, 1)
#' @name tau_leaping
NULL

#' @rdname tau_leaping
#' @export
RegularJump <- function(rate, c, numjumps, mark_dist = NULL) {
  if (!is.function(rate) || !is.function(c)) stop("RegularJump needs R functions `rate` and `c`", call. = FALSE)
  # the argument `c` (named as in JumpProcesses) masks base::c here
  structure(list(rate = rate, c = c, numjumps = as.integer(numjumps), mark_dist = mark_dist),
            class = base::c("RegularJump", "AbstractJump"))
}

#' @rdname tau_leaping
#' @export
PureLeaping <- function() .new_aggregator("PureLeaping")

.new_leaping_alg <- function(name, ...) structure(list(name = name, ...), class = c(name, "TauLeapingAlgorithm"))

#' @export
print.TauLeapingAlgorithm <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

#' @rdname tau_leaping
#' @export
SimpleTauLeaping <- function() .new_leaping_alg("SimpleTauLeaping")

#' @rdname tau_leaping
#' @export
SimpleExplicitTauLeaping <- function(epsilon = 0.05) .new_leaping_alg("SimpleExplicitTauLeaping", epsilon = epsilon)

#' @rdname tau_leaping
#' @export
SimpleImplicitTauLeaping <- function(epsilon = 0.05) .new_leaping_alg("SimpleImplicitTauLeaping", epsilon = epsilon)

#' @rdname tau_leaping
#' @export
SimpleTrapezoidalLeaping <- function(epsilon = 0.05) .new_leaping_alg("SimpleTrapezoidalLeaping", epsilon = epsilon)

#' @rdname tau_leaping
#' @export
SimpleAdaptiveTauLeaping <- function(epsilon = 0.05, implicit_alg = SimpleImplicitTauLeaping(),
                                     eigenvalue_check = FALSE, stiffness_ratio_threshold = 1e4,
                                     implicit_epsilon_factor = 10) {
  .new_leaping_alg("SimpleAdaptiveTauLeaping", epsilon = epsilon, implicit_alg = implicit_alg,
                   eigenvalue_check = eigenvalue_check, stiffness_ratio_threshold = stiffness_ratio_threshold,
                   implicit_epsilon_factor = implicit_epsilon_factor)
}

#' @rdname tau_leaping
#' @export
TauLeaping <- function() .new_leaping_alg("TauLeaping")

#' @rdname tau_leaping
#' @export
CaoTauLeaping <- function(epsilon = 0.03) .new_leaping_alg("CaoTauLeaping", epsilon = epsilon)

#' @rdname tau_leaping
#' @export
ImplicitTauLeaping <- function() .new_leaping_alg("ImplicitTauLeaping", theta = 1)

#' @rdname tau_leaping
#' @export
ThetaTrapezoidalTauLeaping <- function(theta = 0.5) .new_leaping_alg("ThetaTrapezoidalTauLeaping", theta = theta)

# _process_saveat: interior save times and the endpoint flags.
.leap_saveat <- function(saveat, tspan, save_start, save_end) {
  t0 <- tspan[1L]; tf <- tspan[2L]
  if (is.null(saveat) || !length(saveat))
    return(list(times = numeric(), start = save_start %||% TRUE, end = save_end %||% TRUE))
  if (length(saveat) == 1L) {
    v <- seq(t0 + saveat, tf, by = saveat)
    if (length(v) && v[length(v)] == tf) v <- v[-length(v)]
    return(list(times = v, start = save_start %||% TRUE, end = save_end %||% TRUE))
  }
  v <- sort(as.numeric(saveat))
  list(times = v[v > t0 & v < tf], start = save_start %||% (t0 %in% v), end = save_end %||% (tf %in% v))
}

# The jumps of a pure leaping problem: the regular jump and the stoichiometry
# of the mass action jumps (net change matrix nu, reactant stoichiometry).
.leap_parts <- function(jp, alg) {
  js <- jp$jumpset
  rj <- jp$regular_jump
  if (!inherits(jp$aggregator, "PureLeaping"))
    warning("When using ", alg$name, ", please pass PureLeaping() as the aggregator to the JumpProblem, i.e. ",
            "call JumpProblem(<DiscreteProblem>, PureLeaping(), ...).", call. = FALSE)
  if (!inherits(jp$prob, "DiscreteProblem"))
    stop(alg$name, " needs a JumpProblem built from a DiscreteProblem", call. = FALSE)
  cbs <- .as_callback_set(jp$kwargs$callback)
  if (js$ncrj || length(js$vrj) || js$nbvrj || length(cbs$continuous) + length(cbs$discrete))
    stop(alg$name, " can only be used with PureLeaping JumpProblems with only RegularJumps or ",
         "MassActionJumps.", call. = FALSE)
  rj_algs <- c("SimpleTauLeaping", "TauLeaping", "CaoTauLeaping", "ImplicitTauLeaping", "ThetaTrapezoidalTauLeaping")
  if (alg$name %in% rj_algs && is.null(rj))
    stop(alg$name, " can only be used with PureLeaping JumpProblems with only RegularJumps.", call. = FALSE)
  if (!(alg$name %in% rj_algs) && !js$nmaj)
    stop(alg$name, " can only be used with PureLeaping JumpProblem with a MassActionJump.", call. = FALSE)
  n <- js$n
  nu <- matrix(0, n, js$nmaj)
  for (j in seq_len(js$nmaj)) {
    ns <- js$maj_ns[[j]]
    nu[ns[, 1L], j] <- ns[, 2L]
  }
  list(js = js, rj = rj, nu = nu, n = n, nj = if (!is.null(rj)) rj$numjumps else js$nmaj)
}

# compute_hor and precompute_reaction_conditions
.leap_hor <- function(js) {
  hor <- vapply(js$maj_rs, function(rs) sum(rs[, 2L]), 0)
  if (any(hor > 3)) .stopf("Reaction %d has order %d, which is not supported (maximum order is 3).",
                           which(hor > 3)[1L], hor[hor > 3][1L])
  max_hor <- numeric(js$n); max_stoich <- numeric(js$n)
  for (j in seq_along(js$maj_rs)) {
    rs <- js$maj_rs[[j]]
    for (r in seq_len(nrow(rs))) {
      i <- rs[r, 1L]; s <- rs[r, 2L]
      if (s <= 0) next
      if (hor[j] > max_hor[i]) {
        max_hor[i] <- hor[j]
        max_stoich[i] <- s
      } else if (hor[j] == max_hor[i]) max_stoich[i] <- max(max_stoich[i], s)
    }
  }
  list(hor = hor, max_hor = max_hor, max_stoich = max_stoich)
}

# compute_gi, equation (27) of Cao et al. (2006)
.leap_gi <- function(x, mh, ms) {
  if (mh <= 1) return(1)
  if (mh == 2) return(if (ms == 1) 2 else if (x > 1) 2 + 1 / (x - 1) else 2)
  if (ms == 1) return(3)
  if (ms == 2) return(if (x > 1) 1.5 * (2 + 1 / (x - 1)) else 3)
  if (x > 2) 3 + 1 / (x - 1) + 2 / (x - 2) else 3
}

# compute_tau, equation (20) of Cao et al. (2006)
.leap_tau <- function(u, rates, nu, hc, epsilon, dtmin) {
  if (all(rates <= 0)) return(dtmin)
  mu <- as.numeric(nu %*% rates)
  sigma2 <- as.numeric(nu^2 %*% rates)
  tau <- Inf
  for (i in seq_along(u)) {
    bound <- max(epsilon * u[i] / .leap_gi(u[i], hc$max_hor[i], hc$max_stoich[i]), 1)
    if (abs(mu[i]) > 0) tau <- min(tau, bound / abs(mu[i]))
    if (sigma2[i] > 0) tau <- min(tau, bound^2 / sigma2[i])
  }
  max(tau, dtmin)
}

.leap_counts <- function(lambda) {
  out <- numeric(length(lambda))
  pos <- lambda > 0
  out[pos] <- stats::rpois(sum(pos), lambda[pos])
  out
}

# solve_implicit: SimpleNewtonRaphson with a finite difference Jacobian on
# u_new - u - tau * nu a(u_new) (implicit) or the trapezoidal average.
.leap_implicit <- function(u, rate, nu, p, t, tau, trapezoidal) {
  a0 <- if (trapezoidal) rate(u, p, t) else NULL
  resid <- function(x) {
    a1 <- rate(x, p, t + tau)
    x - u - tau * as.numeric(nu %*% (if (trapezoidal) (a1 + a0) / 2 else a1))
  }
  x <- as.numeric(u)
  for (it in seq_len(100L)) {
    r <- resid(x)
    if (!all(is.finite(r))) return(list(u = x, ok = FALSE))
    if (max(abs(r)) <= 1e-6) return(list(u = x, ok = TRUE))
    J <- matrix(0, length(x), length(x))
    for (i in seq_along(x)) {
      h <- sqrt(.Machine$double.eps) * max(abs(x[i]), 1)
      xp <- x
      xp[i] <- xp[i] + h
      J[, i] <- (resid(xp) - r) / h
    }
    dx <- tryCatch(solve(J, -r), error = function(e) NULL)
    if (is.null(dx)) return(list(u = x, ok = FALSE))
    x <- x + dx
    if (max(abs(dx)) <= 1e-6 + 1e-6 * max(abs(x))) return(list(u = x, ok = all(is.finite(resid(x)))))
  }
  list(u = x, ok = FALSE)
}

# compute_drift_jacobian and is_stiff
.leap_is_stiff <- function(rates, u, alg, rate, nu, p, t) {
  pos <- rates[rates > 0]
  if (length(pos) <= 1L) return(FALSE)
  if (isTRUE(alg$eigenvalue_check)) {
    a0 <- rate(u, p, t)
    da <- matrix(0, length(a0), length(u))
    for (i in seq_along(u)) {
      h <- sqrt(.Machine$double.eps) * max(abs(u[i]), 1)
      up <- as.numeric(u)
      up[i] <- up[i] + h
      da[, i] <- (rate(up, p, t) - a0) / h
    }
    mags <- abs(Re(eigen(nu %*% da, only.values = TRUE)$values))
    mags <- mags[mags > 1e-10]
    if (length(mags) <= 1L) return(FALSE)
    return(max(mags) / min(mags) > alg$stiffness_ratio_threshold)
  }
  max(pos) / min(pos) > alg$epsilon * sum(u)
}

.leap_solve <- function(jp, alg, opts) {
  parts <- .leap_parts(jp, alg)
  prob <- jp$prob
  js <- parts$js
  p <- prob$p
  t0 <- prob$tspan[1L]; tf <- prob$tspan[2L]
  sv <- .leap_saveat(opts$saveat, prob$tspan, opts$save_start, opts$save_end)
  dtmin <- opts$dtmin %||% 1e-10
  ts <- numeric(); us <- list()
  save <- function(tt, uu) {
    ts[length(ts) + 1L] <<- tt
    us[[length(us) + 1L]] <<- uu
  }
  u <- as.numeric(prob$u0)
  if (sv$start) save(t0, u)
  rj <- parts$rj
  rate <- if (!is.null(rj)) function(x, pp, tt) as.numeric(rj$rate(x, pp, tt)) else {
    function(x, pp, tt) {
      k <- as.numeric(js$maj_ratefun(x, pp, tt))
      vapply(seq_len(js$nmaj), function(j) .evalrxrate(x, js$maj_rs[[j]], k[j]), 0)
    }
  }
  apply_counts <- function(x, tt, counts) {
    if (!is.null(rj)) as.numeric(rj$c(x, p, tt, counts, NULL)) else as.numeric(parts$nu %*% counts)
  }
  if (alg$name %in% c("TauLeaping", "CaoTauLeaping", "ImplicitTauLeaping", "ThetaTrapezoidalTauLeaping")) {
    r <- .sdeleap_loop(alg, u, t0, tf, p, rate, rj, sv, opts, dtmin, save)
    u <- r$u
  } else if (alg$name == "SimpleTauLeaping") {
    dt <- opts$dt %||% stop("dt is required for SimpleTauLeaping.", call. = FALSE)
    if (!is.null(rj$mark_dist)) stop("Mark distributions are currently not supported in SimpleTauLeaping", call. = FALSE)
    nsteps <- as.integer((tf - t0) / dt) + 1L
    si <- 1L
    for (i in seq_len(nsteps - 1L) + 1L) {
      tprev <- t0 + (i - 2L) * dt
      tnew <- tprev + dt
      counts <- .leap_counts(rate(u, p, tprev) * dt)
      u <- u + apply_counts(u, tprev, counts)
      if (!length(sv$times)) save(tnew, u) else {
        while (si <= length(sv$times) && tnew >= sv$times[si]) {
          save(sv$times[si], u)
          si <- si + 1L
        }
      }
    }
  } else {
    hc <- .leap_hor(js)
    nu <- parts$nu
    implicit_alg <- function(a) a$name %in% c("SimpleImplicitTauLeaping", "SimpleTrapezoidalLeaping")
    t <- t0
    si <- 1L
    tau_cap <- Inf
    while (t < tf) {
      rates <- rate(u, p, t)
      if (all(rates <= 0)) break
      use_implicit <- implicit_alg(alg) ||
        (alg$name == "SimpleAdaptiveTauLeaping" && .leap_is_stiff(rates, u, alg, rate, nu, p, t))
      eps <- alg$epsilon
      tau <- if (alg$name == "SimpleAdaptiveTauLeaping" && use_implicit) {
        # compute_tau_implicit: the step from a predicted state with a
        # relaxed epsilon
        te <- .leap_tau(u, rates, nu, hc, eps, dtmin)
        up <- pmax(u + as.numeric(nu %*% rates) * te, 0)
        max(.leap_tau(up, rate(up, p, t + te), nu, hc, eps * alg$implicit_epsilon_factor, dtmin), dtmin)
      } else .leap_tau(u, rates, nu, hc, eps, dtmin)
      tau <- min(tau, tau_cap, tf - t)
      if (si <= length(sv$times) && t + tau > sv$times[si]) tau <- sv$times[si] - t
      if (use_implicit) {
        ialg <- if (alg$name == "SimpleAdaptiveTauLeaping") alg$implicit_alg else alg
        sol_i <- .leap_implicit(u, rate, nu, p, t, tau, ialg$name == "SimpleTrapezoidalLeaping")
        if (!sol_i$ok) {
          if (tau <= dtmin)
            .stopf("%s failed to converge at t = %g with the smallest permitted step dtmin = %g.",
                   alg$name, t, dtmin)
          tau_cap <- tau / 2
          next
        }
        rates <- rate(sol_i$u, p, t + tau)
      }
      counts <- .leap_counts(pmax(rates * tau, 0))
      unew <- u + (if (use_implicit) as.numeric(nu %*% counts) else apply_counts(u, t, counts))
      if (any(unew < 0)) {
        # a negative population rejects the step. The implicit loops halve
        # tau; the explicit loop of JumpProcesses halves a tau that the next
        # pass recomputes, so it resamples the counts with the same step.
        if (alg$name != "SimpleExplicitTauLeaping" && tau <= dtmin) break
        if (use_implicit) tau_cap <- tau / 2
        next
      }
      tnew <- t + tau
      if (!length(sv$times) || (si <= length(sv$times) && tnew >= sv$times[si])) {
        save(tnew, unew)
        if (length(sv$times) && tnew >= sv$times[si]) si <- si + 1L
      }
      u <- unew
      t <- tnew
      tau_cap <- Inf
    }
  }
  if (sv$end && (!length(ts) || ts[length(ts)] != tf)) save(tf, u)
  U <- if (length(us)) do.call(rbind, us) else matrix(numeric(), 0L, length(u))
  colnames(U) <- names(prob$u0)
  const <- function(theta, dt, y0, y1, k) y0
  steps <- if (length(ts)) lapply(seq_along(ts), function(i) {
    te <- if (i < length(ts)) ts[i + 1L] else ts[i]
    list(t0 = ts[i], dt = max(te - ts[i], 0), tend = te, y0 = U[i, ], y1 = U[i, ], k = NULL, interp = const)
  })
  .new_ode_solution(t = ts, u = U, prob = prob, alg = alg, retcode = "Success", steps = steps, p = p,
                    dense = TRUE)
}

# The tau leaping methods of StochasticDiffEq.jl (perform_step/tau_leaping.jl,
# the TauLeaping and CaoTauLeaping step controllers of
# integrators/stepsize_controllers.jl, TauLeapingDrift).
.sdeleap_loop <- function(alg, u, t0, tf, p, rate, rj, sv, opts, dtmin, save) {
  reltol <- opts$reltol %||% 1e-2
  gamma <- opts$gamma %||% 0.9
  qmax <- opts$qmax %||% 10
  dtmax <- opts$dtmax %||% (tf - t0)
  name <- alg$name
  fixed <- name %in% c("ImplicitTauLeaping", "ThetaTrapezoidalTauLeaping")
  adaptive <- !fixed && (opts$adaptive %||% TRUE)
  if ((fixed || !adaptive) && is.null(opts$dt)) stop("Fixed timestep methods require a choice of dt", call. = FALSE)
  cfun <- function(x, tt, k) as.numeric(rj$c(x, p, tt, k, NULL))
  drift <- function(x, tt) cfun(x, tt, rate(x, p, tt))
  # rate_constants of solve.jl: the rates at u = 1
  rate_constants <- rate(rep(1, length(u)), p, t0)
  # net change of each jump, for Cao's step size
  nu_of <- function(x, tt) vapply(seq_len(rj$numjumps), function(j) {
    e <- numeric(rj$numjumps)
    e[j] <- 1
    cfun(x, tt, e)
  }, numeric(length(x)))
  hc <- list(max_hor = numeric(length(u)), max_stoich = numeric(length(u)))
  dt <- opts$dt %||% (dtmax / 1e6)
  t <- t0
  si <- 1L
  currate <- rate(u, p, t)
  while (t < tf) {
    if (name == "CaoTauLeaping") dt <- min(.leap_tau(u, currate, matrix(nu_of(u, t), length(u)), hc,
                                                     alg$epsilon, dtmin), dtmax)
    h <- min(dt, tf - t)
    if (si <= length(sv$times) && t + h > sv$times[si]) h <- sv$times[si] - t
    k <- .leap_counts(currate * h)
    repeat {
      if (fixed) {
        # X1 = tmp + theta z, z = dt drift(X1) (Newton on the stage equation)
        th <- alg$theta
        tmp <- u + cfun(u, t, k) - th * h * drift(u, t)
        x <- tmp + th * h * drift(u, t)
        for (it in seq_len(50L)) {
          r <- x - tmp - th * h * drift(x, t + h)
          if (max(abs(r)) <= 1e-8 * max(1, max(abs(x)))) break
          J <- matrix(0, length(x), length(x))
          for (i in seq_along(x)) {
            hh <- sqrt(.Machine$double.eps) * max(abs(x[i]), 1)
            xp <- x
            xp[i] <- xp[i] + hh
            J[, i] <- (xp - tmp - th * h * drift(xp, t + h) - r) / hh
          }
          x <- x - as.numeric(solve(J, r))
        }
        unew <- x
      } else {
        unew <- u + cfun(u, t, k)
      }
      newrate <- rate(unew, p, t + h)
      if (name != "TauLeaping" || !adaptive) break
      EEst <- max(abs(newrate - currate) / pmax(50 * reltol * currate, rate_constants / h))
      if (!is.finite(EEst) || EEst <= 1) {
        dtnew <- if (is.finite(EEst) && EEst > 0) h * min(gamma / EEst, qmax) else h * qmax
        break
      }
      # rejection: a shorter step with the Poisson bridge of the counts
      hnew <- max(gamma * h / EEst, dtmin)
      k <- stats::rbinom(length(k), k, hnew / h)
      h <- hnew
    }
    t <- t + h
    u <- unew
    currate <- newrate
    if (name == "TauLeaping" && adaptive) dt <- min(dtnew, dtmax)
    if (!length(sv$times) || (si <= length(sv$times) && t >= sv$times[si])) {
      save(t, u)
      if (length(sv$times) && t >= sv$times[si]) si <- si + 1L
    }
  }
  list(u = u)
}
