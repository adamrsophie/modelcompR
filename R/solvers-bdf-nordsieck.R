# NordsieckBDF and DNordsieckBDF of OrdinaryDiffEqBDF: variable order (1 to
# 5) BDF on a propagated Nordsieck history array zn[j] = h^j / j! y^(j)(t_n),
# following SUNDIALS CVODE (nordsieck_utils.jl, the NordsieckBDF steps of
# bdf_perform_step.jl and dae_perform_step.jl, the NordsieckBDFAlgs hooks of
# controllers.jl, nordsieck_stald! and the Nordsieck interpolant). The array
# is the list cache$zn with Julia's indices (zn[[1]] = y, zn[[2]] = h y', ...);
# l and tau have max_order + 3 entries and tq 6, also 1 based as in Julia.
# The step leaves the array predicted; the accept hook commits it with the
# rank one update and stores the dense output columns, the reject and Newton
# failure hooks undo the prediction. Every operation is plain arithmetic
# (the Julia functions are not under @muladd), except the interpolant.

# The norm of the step: calculate_residuals with uprev and u
.nordbdf_wrms <- function(cache, x) .err_norm(x, cache$uprev, cache$unew, cache$abstol, cache$reltol)

# cvSetBDF and cvSetTqBDF: the correction vector l and the test quantities tq
.nordbdf_set_coeffs <- function(cache, dt) {
  l <- cache$l
  tau <- cache$tau
  tq <- cache$tq
  q <- cache$order
  l[1L] <- 1
  l[2L] <- 1
  xi_inv <- 1
  xistar_inv <- 1
  if (q >= 2L) for (i in 2:q) l[i + 1L] <- 0
  alpha0 <- -1
  alpha0_hat <- -1
  hsum <- dt
  if (q > 1L) {
    if (q >= 3L) for (j in 2:(q - 1L)) {
      hsum <- hsum + tau[j - 1L]
      xi_inv <- dt / hsum
      alpha0 <- alpha0 - 1 / j
      for (i in j:1) l[i + 1L] <- l[i + 1L] + l[i] * xi_inv
    }
    alpha0 <- alpha0 - 1 / q
    xistar_inv <- -l[2L] - alpha0
    hsum <- hsum + tau[q - 1L]
    xi_inv <- dt / hsum
    alpha0_hat <- -l[2L] - xi_inv
    for (i in q:1) l[i + 1L] <- l[i + 1L] + l[i] * xistar_inv
  }
  A1 <- 1 - alpha0_hat + alpha0
  A2 <- 1 + q * A1
  tq[2L] <- abs(A1 / (alpha0 * A2))
  tq[5L] <- abs(A2 * xistar_inv / (l[q + 1L] * xi_inv))
  if (cache$qwait == 1L) {
    if (q > 1L) {
      C <- xistar_inv / l[q + 1L]
      A3 <- alpha0 + 1 / q
      A4 <- alpha0_hat + xi_inv
      tq[1L] <- abs(C * (1 - A4 + A3) / A3)
    } else {
      tq[1L] <- 1
    }
    hsum <- hsum + tau[q]
    xi_inv <- dt / hsum
    A5 <- alpha0 - 1 / (q + 1)
    A6 <- alpha0_hat - xi_inv
    tq[3L] <- abs((1 - A6 + A5) / A2 / (xi_inv * (q + 2) * A5))
  }
  cache$l <- l
  cache$tq <- tq
  invisible()
}

# cvPredict (a Pascal triangle shift) and cvRestore, both idempotent
.nordbdf_predict <- function(cache) {
  if (cache$predicted) return(invisible())
  q <- cache$order
  zn <- cache$zn
  for (k in seq_len(q)) for (j in q:k) zn[[j]] <- zn[[j]] + zn[[j + 1L]]
  cache$zn <- zn
  cache$predicted <- TRUE
  invisible()
}

.nordbdf_restore <- function(cache) {
  if (!cache$predicted) return(invisible())
  q <- cache$order
  zn <- cache$zn
  for (k in seq_len(q)) for (j in q:k) zn[[j]] <- zn[[j]] - zn[[j + 1L]]
  cache$zn <- zn
  cache$predicted <- FALSE
  invisible()
}

# cvRescale: zn[j] = eta^j zn[j]
.nordbdf_rescale <- function(cache, eta) {
  zn <- cache$zn
  factor <- eta
  for (j in seq_len(cache$order)) {
    zn[[j + 1L]] <- zn[[j + 1L]] * factor
    factor <- factor * eta
  }
  cache$zn <- zn
  cache$hscale <- cache$hscale * eta
  invisible()
}

# cvCompleteStep: the rank one update and the step size history
.nordbdf_complete <- function(cache, dt, acor) {
  q <- cache$order
  tau <- cache$tau
  if (q >= 2L) for (i in q:2) tau[i] <- tau[i - 1L]
  if (q == 1L && cache$nst > 0L) tau[2L] <- tau[1L]
  tau[1L] <- dt
  cache$tau <- tau
  zn <- cache$zn
  l <- cache$l
  for (j in 0:q) zn[[j + 1L]] <- zn[[j + 1L]] + l[j + 1L] * acor
  cache$nst <- cache$nst + 1L
  cache$predicted <- FALSE
  cache$qwait <- cache$qwait - 1L
  mo <- cache$max_order
  if (cache$qwait == 1L && q != mo) {
    zn[[mo + 1L]] <- acor
    cache$saved_tq5 <- cache$tq[5L]
    cache$index_acor <- mo
  }
  cache$zn <- zn
  invisible()
}

# cvIncreaseBDF: a new highest column from the saved correction
.nordbdf_increase <- function(cache) {
  mo <- cache$max_order
  l <- cache$l
  l[seq_len(mo + 1L)] <- 0
  l[3L] <- 1
  alpha1 <- 1
  prod <- 1
  xiold <- 1
  alpha0 <- -1
  hsum <- cache$hscale
  q <- cache$order
  if (q > 1L) for (j in 1:(q - 1L)) {
    hsum <- hsum + cache$tau[j + 1L]
    xi <- hsum / cache$hscale
    prod <- prod * xi
    alpha0 <- alpha0 - 1 / (j + 1)
    alpha1 <- alpha1 + 1 / xi
    for (i in (j + 2L):2) l[i + 1L] <- l[i + 1L] * xiold + l[i]
    xiold <- xi
  }
  cache$l <- l
  A1 <- (-alpha0 - alpha1) / prod
  zn <- cache$zn
  zn[[q + 2L]] <- A1 * zn[[cache$index_acor + 1L]]
  dst <- zn[[q + 2L]]
  if (q >= 2L) for (j in 2:q) zn[[j + 1L]] <- zn[[j + 1L]] + l[j + 1L] * dst
  cache$zn <- zn
  invisible()
}

# cvDecreaseBDF: drop the highest column
.nordbdf_decrease <- function(cache) {
  mo <- cache$max_order
  l <- cache$l
  l[seq_len(mo + 1L)] <- 0
  l[3L] <- 1
  hsum <- 0
  q <- cache$order
  if (q >= 3L) for (j in 1:(q - 2L)) {
    hsum <- hsum + cache$tau[j]
    xi <- hsum / cache$hscale
    for (i in (j + 2L):2) l[i + 1L] <- l[i + 1L] * xi + l[i]
  }
  cache$l <- l
  zn <- cache$zn
  zq <- zn[[q + 1L]]
  if (q >= 3L) for (j in 2:(q - 1L)) zn[[j + 1L]] <- zn[[j + 1L]] - l[j + 1L] * zq
  cache$zn <- zn
  invisible()
}

.nordbdf_adjust_order <- function(cache, deltaq) {
  if (cache$order == 2L && deltaq != 1L) return(invisible())
  if (deltaq == 1L) .nordbdf_increase(cache) else if (deltaq == -1L) .nordbdf_decrease(cache)
  invisible()
}

# cvAdjustParams: a pending order change and the rescaling to the step size
.nordbdf_prepare <- function(cache, dt) {
  if (cache$qprime != cache$order) {
    .nordbdf_adjust_order(cache, cache$qprime - cache$order)
    cache$order <- cache$qprime
    cache$qwait <- cache$order + 1L
  }
  if (dt != cache$hscale && cache$hscale != 0) .nordbdf_rescale(cache, dt / cache$hscale)
  cache$hscale <- dt
  invisible()
}

# The candidate step ratios of orders q - 1 and q + 1 (BIAS1 = 6, BIAS3 =
# 10, ADDON = 1 / 1000000)
.nordbdf_etaqm1 <- function(cache) {
  if (cache$order <= 1L) return(0)
  ddn <- .nordbdf_wrms(cache, cache$zn[[cache$order + 1L]]) * cache$tq[1L]
  1 / ((6 * ddn)^(1 / cache$order) + 1e-6)
}

.nordbdf_etaqp1 <- function(cache, dt) {
  if (cache$order == cache$max_order) return(0)
  if (cache$saved_tq5 == 0) return(0)
  L <- cache$order + 1L
  # (dt / tau[2])^L: Julia's power by squaring for an integer exponent
  cquot <- (cache$tq[5L] / cache$saved_tq5) * .simple_powi(dt / cache$tau[2L], L)
  dup <- .nordbdf_wrms(cache, cache$acor - cquot * cache$zn[[cache$max_order + 1L]]) * cache$tq[3L]
  1 / ((10 * dup)^(1 / (L + 1)) + 1e-6)
}

# cvChooseEta (ETA_MAX_FX = 3 / 2)
.nordbdf_choose_eta <- function(cache) {
  etam <- max(cache$etaqm1, max(cache$etaq, cache$etaqp1))
  if (etam < 1.5) {
    cache$eta <- 1
    cache$qprime <- cache$order
    return(invisible())
  }
  if (etam == cache$etaq) {
    cache$eta <- cache$etaq
    cache$qprime <- cache$order
  } else if (etam == cache$etaqm1) {
    cache$eta <- cache$etaqm1
    cache$qprime <- cache$order - 1L
  } else {
    cache$eta <- cache$etaqp1
    cache$qprime <- cache$order + 1L
    zn <- cache$zn
    zn[[cache$max_order + 1L]] <- cache$acor
    cache$zn <- zn
  }
  invisible()
}

.nordbdf_set_eta <- function(cache) {
  cache$eta <- if (cache$eta < 1.5) 1 else min(cache$eta, cache$etamax)
  invisible()
}

# The dense output columns of the committed array (unused columns zero)
.nordbdf_store_k <- function(cache) {
  n <- length(cache$zn[[1L]])
  cache$kenv$k <- lapply(0:cache$max_order, function(j) if (j <= cache$order) cache$zn[[j + 1L]] else numeric(n))
  invisible()
}

# The cold start (order 1): zn = (uprev, dt du) with du = f(uprev) for an
# ODE and integrator.du for a DAE
.nordbdf_start <- function(cache, u, dt, du) {
  n <- length(u)
  zn <- rep(list(numeric(n)), cache$max_order + 1L)
  zn[[1L]] <- u
  zn[[2L]] <- dt * du
  cache$zn <- zn
  cache$order <- 1L
  cache$qprime <- 1L
  cache$qwait <- 2L
  cache$nef <- 0L
  cache$ncf <- 0L
  cache$hscale <- dt
  cache$eta <- 1
  cache$etamax <- 10000
  cache$saved_tq5 <- 0
  cache$index_acor <- cache$max_order
  cache$predicted <- FALSE
  .stald_reset(cache$stald)
  cache$tau[] <- 0
  invisible()
}

# cvBDFStab through nordsieck_stald!: whether STALD lowers the order
.nordbdf_stald <- function(cache, dsm) {
  if (!cache$stald$enabled) return(FALSE)
  q <- cache$order
  if (q < 3L) return(FALSE)
  eps <- .Machine$double.eps
  fact <- 1
  for (i in seq_len(q - 1L)) fact <- fact * i
  acnrm <- dsm / max(cache$tq[2L], eps)
  sq <- fact * q * (q + 1) * acnrm / max(cache$tq[5L], eps)
  sqm1 <- fact * q * .nordbdf_wrms(cache, cache$zn[[q + 1L]])
  sqm2 <- fact * .nordbdf_wrms(cache, cache$zn[[q]])
  .stald_collect(cache$stald, q, sqm2, sqm1, sq)
  .stald_check(cache$stald, q)
}

.nordbdf_init <- function(ctx, t, u, max_order, stald, dae) {
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  # NLNewton(max_iter = 3): CVODE's MAXCOR
  cache$nl$maxiters <- 3L
  mo <- max_order
  n <- length(u)
  cache$max_order <- mo
  cache$zn <- rep(list(numeric(n)), mo + 1L)
  cache$acor <- numeric(n)
  cache$l <- numeric(mo + 3L)
  cache$tau <- numeric(mo + 3L)
  cache$tq <- numeric(6L)
  cache$order <- 1L
  cache$qprime <- 1L
  cache$qwait <- 2L
  cache$nst <- 0L
  cache$nef <- 0L
  cache$ncf <- 0L
  cache$index_acor <- mo
  cache$hscale <- 0
  cache$eta <- 1
  cache$etamax <- 10000
  cache$etaq <- 1
  cache$etaqm1 <- 0
  cache$etaqp1 <- 0
  cache$saved_tq5 <- 0
  cache$predicted <- FALSE
  cache$disc <- FALSE
  cache$stald <- .stald_new(enabled = stald)
  if (dae) {
    cache$fsal <- ctx$du
    cache$nf0 <- 0L
  } else {
    cache$fsal <- ctx$f(u, ctx$p, t)
  }
  cache
}

.nordbdf_step <- function(ctx, t, u, dt, cache, dae) {
  # nordsieck_needs_start: the first step, after an event and after a
  # restart requested by the reject hook
  disc <- isTRUE(ctx$derivative_discontinuity) || isTRUE(cache$disc)
  cache$disc <- FALSE
  if (disc) ctx$derivative_discontinuity <- TRUE
  if (cache$nst == 0L || disc) .nordbdf_start(cache, u, dt, if (dae) ctx$du else cache$fsal)
  .nordbdf_prepare(cache, dt)
  .nordbdf_predict(cache)
  .nordbdf_set_coeffs(cache, dt)
  zn <- cache$zn
  ypred <- zn[[1L]]
  l1 <- cache$l[2L]
  # has_special_newton_error: the Newton increment norm times tq[2]
  ctx$newton_scale <- cache$tq[2L]
  if (dae) {
    # f((zn[2] + l1 z) / dt, ypred + z) = 0: tmp = zn[2] / l1, u0 = ypred
    nl <- .nlsolve(cache$nl, ctx, t, dt, numeric(length(u)), zn[[2L]] / l1, 1 / l1, 1, t + dt, u, u0 = ypred)
  } else {
    tmp <- .mass_mul(ctx$mass, (l1 * ypred - zn[[2L]]) / dt)
    nl <- .nlsolve(cache$nl, ctx, t, dt, ypred, tmp, 1 / l1, 1, t + dt, u)
  }
  out <- .bdfk_counts(nl, cache)
  if (!nl$ok) return(.bdfk_fail(out, u))
  if (dae) {
    acor <- nl$z
    unew <- ypred + nl$z
  } else {
    unew <- nl$z
    acor <- unew - ypred
  }
  cache$acor <- acor
  cache$tdir <- if (dt < 0) -1 else 1
  cache$uprev <- u
  cache$unew <- unew
  cache$abstol <- ctx$abstol
  cache$reltol <- ctx$reltol
  EEst <- 1
  if (ctx$adaptive) EEst <- cache$tq[2L] * .nordbdf_wrms(cache, acor)
  cache$fsallast <- (zn[[2L]] + l1 * acor) / dt
  # the dense output columns, filled when the step is committed
  cache$kenv <- new.env(parent = emptyenv())
  if (!ctx$adaptive) {
    # _nordsieck_finish_fixed!
    .nordbdf_complete(cache, dt, acor)
    .nordbdf_store_k(cache)
    if (cache$qwait <= 0L) {
      cache$qwait <- 2L
      if (cache$order < cache$max_order) cache$qprime <- cache$order + 1L
    }
    if (!dae) cache$fsal <- cache$fsallast
  }
  res <- c(out, list(u = unew, EEst = EEst, k = cache$kenv))
  # DNordsieckBDF never updates integrator.du
  if (dae) res$du <- ctx$du
  res
}

# step_accept_controller! of the Nordsieck methods; `dt` is the magnitude of
# integrator.dt, the hooks return magnitudes
.nordbdf_accept <- function(cache, q, dt, co) {
  dt <- cache$tdir * dt
  dsm <- cache$EEst_cur
  stald_reduce <- .nordbdf_stald(cache, dsm)
  .nordbdf_complete(cache, dt, cache$acor)
  cache$nef <- 0L
  cache$ncf <- 0L
  .nordbdf_store_k(cache)
  if (!isTRUE(cache$dae)) cache$fsal <- cache$fsallast
  if (cache$etamax == 1) {
    # a failure earlier in this step forbids growth
    cache$qwait <- max(cache$qwait, 2L)
    cache$qprime <- cache$order
    cache$eta <- 1
  } else {
    cache$etaq <- 1 / ((6 * dsm)^(1 / (cache$order + 1)) + 1e-6)
    if (cache$qwait != 0L) {
      cache$eta <- cache$etaq
      cache$qprime <- cache$order
      .nordbdf_set_eta(cache)
    } else {
      cache$qwait <- 2L
      cache$etaqm1 <- .nordbdf_etaqm1(cache)
      cache$etaqp1 <- .nordbdf_etaqp1(cache, dt)
      .nordbdf_choose_eta(cache)
      .nordbdf_set_eta(cache)
    }
  }
  if (stald_reduce && cache$qprime > 1L) {
    cache$qprime <- min(cache$qprime, cache$order - 1L)
    cache$eta <- min(cache$eta, 1)
  }
  cache$etamax <- 10
  abs(dt) * min(cache$eta, if (isTRUE(cache$first)) co$qmax_first else co$qmax)
}

# step_reject_controller! (MXNEF1 = 2, SMALL_NEF = 2, ETA_MIN_EF = 1 / 10,
# ETA_MAX_EF = 1 / 5, LONG_WAIT = 10)
.nordbdf_reject <- function(cache, EEst, dt, co) {
  .nordbdf_restore(cache)
  cache$nef <- cache$nef + 1L
  cache$etamax <- 1
  if (cache$nef <= 2L) {
    eta <- 1 / ((6 * EEst)^(1 / (cache$order + 1)) + 1e-6)
    eta <- max(0.1, eta)
    if (cache$nef >= 2L) eta <- min(eta, 0.2)
    cache$eta <- eta
  } else if (cache$order > 1L) {
    cache$eta <- 0.1
    .nordbdf_adjust_order(cache, -1L)
    cache$order <- cache$order - 1L
    cache$qprime <- cache$order
    cache$qwait <- cache$order + 1L
  } else {
    # restart the history: derivative_discontinuity!(integrator, true)
    cache$eta <- 0.1
    cache$qwait <- 10L
    cache$disc <- TRUE
  }
  .nordbdf_rescale(cache, cache$eta)
  abs(cache$hscale)
}

# post_newton_controller! (ETA_CF = 1 / 4)
.nordbdf_newton_fail <- function(cache, dt, co) {
  .nordbdf_restore(cache)
  cache$etamax <- 1
  cache$ncf <- cache$ncf + 1L
  if (cache$ncf >= 3L && cache$order > 1L) {
    .nordbdf_adjust_order(cache, -1L)
    cache$order <- cache$order - 1L
    cache$qprime <- cache$order
    cache$qwait <- cache$order + 1L
    cache$ncf <- 0L
  }
  cache$eta <- 0.25
  .nordbdf_rescale(cache, cache$eta)
  abs(cache$hscale)
}

# The Nordsieck polynomial y(t_n + s h) = sum_j zn[j] s^j, s = theta - 1
.nordbdf_interp <- function(theta, dt, y0, y1, k) {
  kk <- if (is.environment(k)) k$k else NULL
  if (is.null(kk)) return((1 - theta) * y0 + theta * y1)
  s <- theta - 1
  n <- length(kk)
  out <- kk[[n]]
  for (j in (n - 1L):1) out <- .fma(out, s, kk[[j]])
  out
}

.nordbdf_alg <- function(name, dae, max_order, stald, qmax, qsteady_min, qsteady_max, step_limiter = NULL) {
  max_order <- as.integer(max_order)
  if (length(max_order) != 1L || is.na(max_order) || max_order < 1L || max_order > 5L)
    stop("max_order must be between 1 and 5", call. = FALSE)
  if (!is.null(step_limiter) && !is.function(step_limiter))
    stop("step_limiter must be a function(u, integ, p, t)", call. = FALSE)
  stald <- isTRUE(stald)
  .new_alg(name, 1L, step = function(ctx, t, u, dt, cache) .nordbdf_step(ctx, t, u, dt, cache, dae),
           interp = .nordbdf_interp,
           init = function(ctx, t, u) {
             cache <- .nordbdf_init(ctx, t, u, max_order, stald, dae)
             cache$dae <- dae
             cache
           },
           on_event = if (!dae) .implicit_on_event, stiff = TRUE, dae = dae, mass_matrix = !dae,
           controller = "BDF", qmax = qmax, qsteady_min = qsteady_min, qsteady_max = qsteady_max,
           gamma = 1, bdf_stepsize = .bdfk_stepsize, bdf_accept = .nordbdf_accept,
           bdf_reject = .nordbdf_reject, newton_fail = .nordbdf_newton_fail, on_reject = .nl_on_reject,
           max_order = max_order, stald = stald, step_limiter = step_limiter)
}

#' Nordsieck BDF methods
#'
#' `NordsieckBDF()` is the variable order (1 to `max_order`) BDF method of
#' SUNDIALS CVODE as implemented in OrdinaryDiffEqBDF: the history is the
#' Nordsieck array `zn[j] = h^j / j! y^(j)`, predicted by a Pascal triangle
#' shift, rescaled by `eta^j` when the step size changes and updated by a
#' rank one correction, with CVODE's step size and order selection. The
#' Newton iteration stops after 3 iterations (CVODE's `MAXCOR`) and its
#' increment norm is scaled into the units of the local error test. The
#' dense output is the Nordsieck polynomial. `DNordsieckBDF()` is the same
#' method for a fully implicit [DAEProblem()] (the corrector solves
#' `F((zn[2] + l1 z) / h, ypred + z, p, t) = 0`, IDA's leading coefficient).
#'
#' @param max_order Maximum order (1 to 5).
#' @param stald Stability limit detection for orders 3 to 5 (CVODE's
#'   algorithm, off by default as in CVODE).
#' @param qmax Largest step growth factor.
#' @param qsteady_min,qsteady_max Settings of the controller kept for
#'   compatibility; the method keeps the step size unless it may grow by
#'   1.5 (CVODE's rule).
#' @param step_limiter A function `limiter(u, integ, p, t)` returning the
#'   state after each accepted step.
#' @param autodiff As for [ImplicitEuler()].
#' @return An algorithm object.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -50 * (u - t * t), 0, c(0, 1))
#' sol <- solve(prob, NordsieckBDF(), reltol = 1e-6)
#' sol(1)
#' f <- function(du, u, p, t) c(du[1] + u[1], u[2] - u[1])
#' dprob <- DAEProblem(f, c(-1, 0), c(1, 1), c(0, 1), differential_vars = c(TRUE, FALSE))
#' solve(dprob, DNordsieckBDF())$retcode
#' @name nordsieck_bdf_methods
NULL

#' @rdname nordsieck_bdf_methods
#' @export
NordsieckBDF <- function(max_order = 5L, stald = FALSE, qmax = 10, qsteady_min = 1, qsteady_max = 1,
                         step_limiter = NULL, autodiff = NULL) .with_autodiff(autodiff, {
  .nordbdf_alg("NordsieckBDF", FALSE, max_order, stald, qmax, qsteady_min, qsteady_max, step_limiter)
})

#' @rdname nordsieck_bdf_methods
#' @export
DNordsieckBDF <- function(max_order = 5L, stald = FALSE, qmax = 10, qsteady_min = 1, qsteady_max = 1,
                          autodiff = NULL) .with_autodiff(autodiff, {
  .nordbdf_alg("DNordsieckBDF", TRUE, max_order, stald, qmax, qsteady_min, qsteady_max)
})
