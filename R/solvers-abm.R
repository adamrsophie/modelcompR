# Adams-Bashforth and Adams-Moulton methods of
# OrdinaryDiffEqAdamsBashforthMoulton 2.2.1: adams_bashforth_moulton_perform_step.jl
# (the out of place ConstantCache steps), adams_utils.jl (the divided
# differences and the coefficients of Hairer, Norsett and Wanner, Solving
# ODEs I, III.5), adams_bashforth_moulton_caches.jl and alg_utils.jl; the
# starting steps are the BS3 and RK4 steps of OrdinaryDiffEqLowOrderRK.
#
# Deliberate differences from the installed Julia code (tools/abm-probe.jl
# runs Julia with the same corrections for the reference values):
# * VCAB3, VCAB4, VCAB5, VCABM3, VCABM4, VCABM5: a rejected starting step
#   leaves the step counter, the step size history and the divided
#   differences as they were (.abm_startup_undo). Julia keeps them, so the
#   retried step builds its history on the rejected one.
# * VCABM: at order 1 the corrector starts from uprev (.vcabm_base); Julia
#   adds the Euler predictor term twice, so u' = 1 gives u(0.1) = 0.2.
# * VCABM: with adaptive = FALSE the derivative at the end of each step is
#   evaluated (.vcabm_fixed_fsal); Julia never sets it, so every later step
#   starts from a zero derivative.
# * VCABM: a derivative discontinuity (a callback that changes u) restarts
#   at order 1 (.vcabm_restart), as the fixed order methods restart; Julia
#   keeps the divided differences from before the discontinuity.
# The installed ABM32, ABM43 and ABM54 predict with the starting method on
# every step (their step counter stops at 2, 3 and 4); this is reproduced.

# fma(a, b, c): Julia's muladd of two Float64 scalars, which LLVM fuses on
# processors with a fused multiply add, correctly rounded. The product is
# split exactly (Dekker and Veltkamp) into p + e, and p + c into s + u
# (Knuth's TwoSum), so the exact result is s + (u + e). Rounding u + e to t
# and then s + t to r rounds twice only when t is inexact and s + t falls on
# a rounding tie (r + 2d, d the error of s + t, is then the neighbor of r).
# A zero factor gives s, except for a zero result, whose sign the sum may
# lose. The ties, the zero results of a zero factor, a tiny product (whose
# low parts underflow) and a huge result (where the split may overflow,
# which a NaN shows) go to the exact .pw_fma() (rounding to odd,
# R/optimization-prima-powell.R).
.fma <- function(a, b, c) {
  p <- a * b
  ta <- a * 134217729; ah <- ta - (ta - a); al <- a - ah
  tb <- b * 134217729; bh <- tb - (tb - b); bl <- b - bh
  e <- ((ah * bh - p) + ah * bl + al * bh) + al * bl
  s <- p + c
  z <- s - p
  u <- (p - (s - z)) + (c - z)
  t <- u + e
  r <- s + t
  z <- t - u
  w <- r - s
  d <- (s - (r - w)) + (t - w)
  ok <- (abs(p) >= 0x1p-960 | (a == 0 | b == 0) & r != 0) & abs(r) <= 0x1p1000 &
    ((u - (t - z)) + (e - z) == 0 | d == 0 | (r + 2 * d) - r != 2 * d)
  if (!anyNA(ok) && all(ok)) return(r)
  bad <- which(!ok | is.na(ok))
  m <- length(r)
  a <- rep_len(a, m)[bad]
  b <- rep_len(b, m)[bad]
  # a zero factor gives the exact zero product plus c
  zero <- which((a == 0 | b == 0) & rep_len(p, m)[bad] == 0)
  if (length(zero)) {
    r[bad[zero]] <- s[bad[zero]]
    bad <- bad[-zero]; a <- a[-zero]; b <- b[-zero]
  }
  if (length(bad)) r[bad] <- .pw_fma(a, b, rep_len(c, m)[bad])
  r
}

# ---------------------------------------------------------------------------
# Fixed step Adams-Bashforth and Adams-Bashforth-Moulton methods

# Ralston's second order step (the start of AB3 and ABM32).
.abm_ralston <- function(ctx, t, u, dt, k1) {
  ralk2 <- ctx$f(u + 2 * dt / 3 * k1, ctx$p, t + 2 * dt / 3)
  u + (dt / 4) * (k1 + 3 * ralk2)
}

# The classical Runge-Kutta step (the start of the other methods).
.abm_rk4 <- function(ctx, t, u, dt, k1) {
  f <- ctx$f; p <- ctx$p
  halfdt <- dt / 2
  k2 <- f(u + halfdt * k1, p, t + halfdt)
  k3 <- f(u + halfdt * k2, p, t + halfdt)
  k4 <- f(u + dt * k3, p, t + dt)
  u + (dt / 6) * (2 * (k2 + k3) + (k1 + k4))
}

# The steps of AB3ConstantCache, AB4ConstantCache and AB5ConstantCache on the
# history s = list(step, k2, k3, k4, k5): the new state, history and count
# of evaluations (without the derivative at the end of the step).
.ab3_core <- function(ctx, t, u, dt, k1, s) {
  if (s$step <= 2L) {
    cnt <- s$step
    s$step <- cnt + 1L
    unew <- .abm_ralston(ctx, t, u, dt, k1)
    if (cnt == 1L) s$k3 <- k1 else s$k2 <- k1
    return(list(u = unew, s = s, nf = 1L))
  }
  unew <- u + (dt / 12) * (23 * k1 - 16 * s$k2 + 5 * s$k3)
  s$k3 <- s$k2; s$k2 <- k1
  list(u = unew, s = s, nf = 0L)
}

.ab4_core <- function(ctx, t, u, dt, k1, s) {
  if (s$step <= 3L) {
    cnt <- s$step
    s$step <- cnt + 1L
    unew <- .abm_rk4(ctx, t, u, dt, k1)
    if (cnt == 1L) s$k4 <- k1 else if (cnt == 2L) s$k3 <- k1 else s$k2 <- k1
    return(list(u = unew, s = s, nf = 3L))
  }
  unew <- u + (dt / 24) * (55 * k1 - 59 * s$k2 + 37 * s$k3 - 9 * s$k4)
  s$k4 <- s$k3; s$k3 <- s$k2; s$k2 <- k1
  list(u = unew, s = s, nf = 0L)
}

.ab5_core <- function(ctx, t, u, dt, k1, s) {
  if (s$step <= 4L) {
    cnt <- s$step
    s$step <- cnt + 1L
    unew <- .abm_rk4(ctx, t, u, dt, k1)
    s[[c("k5", "k4", "k3", "k2")[cnt]]] <- k1
    return(list(u = unew, s = s, nf = 3L))
  }
  unew <- u + (dt / 720) * (1901 * k1 - 2774 * s$k2 + 2616 * s$k3 - 1274 * s$k4 + 251 * s$k5)
  s$k5 <- s$k4; s$k4 <- s$k3; s$k3 <- s$k2; s$k2 <- k1
  list(u = unew, s = s, nf = 0L)
}

# The predictor-corrector steps of ABM32ConstantCache, ABM43ConstantCache and
# ABM54ConstantCache. The predictor runs on a copy of the history whose step
# counter is the corrector's, which stays at 2, 3 and 4 after the start.
.abm32_core <- function(ctx, t, u, dt, k1, s) {
  if (s$step == 1L) {
    s$step <- 2L
    unew <- .abm_ralston(ctx, t, u, dt, k1)
    s$k2 <- k1
    return(list(u = unew, s = s, nf = 1L))
  }
  pr <- .ab3_core(ctx, t, u, dt, k1, s)
  k <- ctx$f(pr$u, ctx$p, t + dt)
  unew <- u + (dt / 12) * (5 * k + 8 * k1 - s$k2)
  s$k3 <- s$k2; s$k2 <- k1
  list(u = unew, s = s, nf = pr$nf + 1L)
}

.abm43_core <- function(ctx, t, u, dt, k1, s) {
  if (s$step <= 2L) {
    cnt <- s$step
    s$step <- cnt + 1L
    unew <- .abm_rk4(ctx, t, u, dt, k1)
    if (cnt == 1L) s$k3 <- k1 else s$k2 <- k1
    return(list(u = unew, s = s, nf = 3L))
  }
  pr <- .ab4_core(ctx, t, u, dt, k1, s)
  k <- ctx$f(pr$u, ctx$p, t + dt)
  unew <- u + (dt / 24) * (9 * k + 19 * k1 - 5 * s$k2 + s$k3)
  s$k4 <- s$k3; s$k3 <- s$k2; s$k2 <- k1
  list(u = unew, s = s, nf = pr$nf + 1L)
}

.abm54_core <- function(ctx, t, u, dt, k1, s) {
  if (s$step <= 3L) {
    cnt <- s$step
    s$step <- cnt + 1L
    unew <- .abm_rk4(ctx, t, u, dt, k1)
    s[[c("k4", "k3", "k2")[cnt]]] <- k1
    return(list(u = unew, s = s, nf = 3L))
  }
  pr <- .ab5_core(ctx, t, u, dt, k1, s)
  k <- ctx$f(pr$u, ctx$p, t + dt)
  unew <- u + (dt / 720) * (251 * k + 646 * k1 - 264 * s$k2 + 106 * s$k3 - 19 * s$k4)
  s$k5 <- s$k4; s$k4 <- s$k3; s$k3 <- s$k2; s$k2 <- k1
  list(u = unew, s = s, nf = pr$nf + 1L)
}

# The derivative at the start of the step is fsal; a discontinuity restarts
# the history (cache.step = 1). Each step ends with the derivative at its
# end (Hermite dense output).
.ab_alg <- function(name, order, core) {
  init <- function(ctx, t, u) {
    z <- 0 * u
    list(fsal = ctx$f(u, ctx$p, t), step = 1L, k2 = z, k3 = z, k4 = z, k5 = z)
  }
  step <- function(ctx, t, u, dt, cache) {
    if (isTRUE(ctx$derivative_discontinuity)) cache$step <- 1L
    k1 <- cache$fsal
    r <- core(ctx, t, u, dt, k1, cache)
    fnew <- ctx$f(r$u, ctx$p, t + dt)
    s <- r$s
    s$fsal <- fnew
    list(u = r$u, err = NULL, k = list(k1, fnew), cache = s, nf = r$nf + 1L)
  }
  .new_alg(name, order, adaptive = FALSE, step = step, init = init, on_event = .abm_on_event)
}

# After a callback changed u: the derivative is evaluated again and the next
# step sees the discontinuity.
.abm_on_event <- function(ctx, t, u, cache) {
  cache$fsal <- ctx$f(u, ctx$p, t)
  list(cache = cache, nf = 1L)
}

#' Adams-Bashforth and Adams-Moulton methods
#'
#' Explicit linear multistep methods of OrdinaryDiffEq.jl
#' (OrdinaryDiffEqAdamsBashforthMoulton), for nonstiff problems whose right
#' hand side is expensive: after a few starting steps each step costs one or
#' two evaluations of `f`.
#'
#' * `AB3()`, `AB4()`, `AB5()`: the fixed step Adams-Bashforth methods of
#'   orders 3, 4 and 5 (pass `dt`). AB3 starts with Ralston's second order
#'   method, AB4 and AB5 with the classical Runge-Kutta method.
#' * `ABM32()`, `ABM43()`, `ABM54()`: fixed step predictor-corrector pairs
#'   (Adams-Bashforth predictor, Adams-Moulton corrector) of orders 3, 4
#'   and 5. As in the Julia version the predictor is the starting method
#'   (Ralston's method for ABM32, the classical Runge-Kutta method for the
#'   others) on every step.
#' * `VCAB3()`, `VCAB4()`, `VCAB5()`: adaptive variable step Adams-Bashforth
#'   methods of orders 3, 4 and 5 in divided difference form, started with
#'   the Bogacki-Shampine method (VCAB3) or the classical Runge-Kutta method
#'   with Shampine's residual error estimate.
#' * `VCABM3()`, `VCABM4()`, `VCABM5()`: adaptive variable step
#'   Adams-Bashforth-Moulton predictor-corrector methods of orders 3, 4
#'   and 5.
#' * `VCABM()`: the adaptive variable step, variable order (1 to 12)
#'   Adams-Bashforth-Moulton method, with the I step size controller.
#'
#' A callback that changes the state restarts the methods (the variable
#' order method at order 1). The dense output is the cubic Hermite
#' interpolant of each step.
#'
#' @param thread Threading option of the Julia methods (ignored).
#' @return An ODE algorithm.
#' @examples
#' prob <- ODEProblem(function(u, p, t) c(u[2], -u[1]), c(0, 1), c(0, 1))
#' solve(prob, AB4(), dt = 0.01)$u[[101]]
#' sol <- solve(prob, VCABM(), reltol = 1e-8, abstol = 1e-8)
#' sol(1)
#' @name adams_methods
NULL

#' @rdname adams_methods
#' @export
AB3 <- function(thread = NULL) .ab_alg("AB3", 3L, .ab3_core)

#' @rdname adams_methods
#' @export
AB4 <- function(thread = NULL) .ab_alg("AB4", 4L, .ab4_core)

#' @rdname adams_methods
#' @export
AB5 <- function(thread = NULL) .ab_alg("AB5", 5L, .ab5_core)

#' @rdname adams_methods
#' @export
ABM32 <- function(thread = NULL) .ab_alg("ABM32", 3L, .abm32_core)

#' @rdname adams_methods
#' @export
ABM43 <- function(thread = NULL) .ab_alg("ABM43", 4L, .abm43_core)

#' @rdname adams_methods
#' @export
ABM54 <- function(thread = NULL) .ab_alg("ABM54", 5L, .abm54_core)

# ---------------------------------------------------------------------------
# Divided differences (adams_utils.jl). The cache is an environment holding
# dts (the last step sizes, newest first), phi_n, phistar_n, phistar_nm1,
# phi_np1 (lists of vectors), beta, c, g, xi and xi0.

# phi_and_phistar!: the divided differences phi_n[1..k] of the derivatives
# and their scaled versions phistar_n from the previous phistar.
.phi_and_phistar <- function(cache, du, k) {
  dts <- cache$dts
  xi <- dts[1L]; xi0 <- 0
  beta <- cache$beta; phi <- cache$phi_n; phis <- cache$phistar_n; old <- cache$phistar_nm1
  beta[1L] <- 1
  phi[[1L]] <- du; phis[[1L]] <- du
  for (i in seq_len(k)[-1L]) {
    xi0 <- xi0 + dts[i]
    beta[i] <- beta[i - 1L] * xi / xi0
    xi <- xi + dts[i]
    phi[[i]] <- phi[[i - 1L]] - old[[i - 1L]]
    phis[[i]] <- beta[i] * phi[[i]]
  }
  cache$beta <- beta; cache$phi_n <- phi; cache$phistar_n <- phis
  cache$xi <- xi; cache$xi0 <- xi0
  invisible(cache)
}

# expand_phi_and_phistar!: one more divided difference (VCABM order choice).
.expand_phi <- function(cache, i) {
  xi0 <- cache$xi0 + cache$dts[i]
  cache$beta[i] <- cache$beta[i - 1L] * cache$xi / xi0
  cache$phi_n[[i]] <- cache$phi_n[[i - 1L]] - cache$phistar_nm1[[i - 1L]]
  cache$phistar_n[[i]] <- cache$beta[i] * cache$phi_n[[i]]
  invisible(cache)
}

# phi_np1!: the divided differences at the new point from its derivative.
.phi_np1 <- function(cache, du, k) {
  phi <- cache$phi_np1; phis <- cache$phistar_n
  phi[[1L]] <- du
  for (i in seq_len(k)[-1L]) phi[[i]] <- phi[[i - 1L]] - phis[[i - 1L]]
  cache$phi_np1 <- phi
  invisible(cache)
}

# g_coefs!: the integration coefficients g[1..k] (scaled by dt) from the
# recurrence for c[i, q] (formulae 5.9 and 5.10).
.g_coefs <- function(cache, k) {
  dts <- cache$dts; cc <- cache$c; g <- cache$g
  dt <- dts[1L]; xi <- dt
  for (i in seq_len(k)) {
    if (i > 2L) xi <- xi + dts[i - 1L]
    for (q in seq_len(k - (i - 1L))) {
      cc[i, q] <- if (i > 2L) .fma(-dt / xi, cc[i - 1L, q + 1L], cc[i - 1L, q])
                  else if (i == 1L) 1 / q else 1 / (q * (q + 1))
    }
    g[i] <- cc[i, 1L] * dt
  }
  cache$c <- cc; cache$g <- g
  invisible(cache)
}

# Coefficients of the implicit Adams methods (gamma star).
.ABM_GAMMASTAR <- c(1, -1 / 2, -1 / 12, -1 / 24, -19 / 720, -3 / 160, -863 / 60480, -275 / 24192,
                    -33953 / 3628800, -0.00789255, -0.00678585, -0.00592406, -0.00523669, -0.0046775,
                    -0.00421495, -0.0038269)

.abm_cache <- function(ctx, t, u, n_dts, n_c, n_phi, order) {
  cache <- new.env(parent = emptyenv())
  z <- 0 * u
  cache$fsal <- ctx$f(u, ctx$p, t)
  cache$fsallast <- z
  cache$step <- 1L
  cache$order <- order
  cache$dts <- numeric(n_dts)
  cache$c <- matrix(0, n_c, n_c)
  cache$g <- numeric(n_c)
  cache$beta <- numeric(n_phi)
  cache$xi <- 0; cache$xi0 <- 0
  cache$phi_n <- cache$phistar_n <- cache$phistar_nm1 <- rep(list(z), n_phi)
  cache$phi_np1 <- rep(list(z), n_phi + 1L)
  cache$undo <- NULL
  cache
}

# A rejected step restores what it changed (the derivative at its start and,
# for the starting steps, the history).
.abm_on_reject <- function(cache) {
  for (nm in names(cache$undo)) cache[[nm]] <- cache$undo[[nm]]
  cache$undo <- NULL
  invisible(cache)
}

.abm_swap <- function(cache) {
  tmp <- cache$phistar_nm1
  cache$phistar_nm1 <- cache$phistar_n
  cache$phistar_n <- tmp
  invisible(cache)
}

# What a rejected starting step restores (see the header).
.abm_startup_undo <- function(cache, k) {
  list(step = k, dts = cache$dts, phistar_n = cache$phistar_n, phistar_nm1 = cache$phistar_nm1)
}

# The classical Runge-Kutta step of RK4ConstantCache with Shampine's
# residual error estimate (the start of VCAB4, VCAB5, VCABM4 and VCABM5).
.rk4_shampine <- function(ctx, t, u, dt, k1) {
  f <- ctx$f; p <- ctx$p
  unew <- .abm_rk4(ctx, t, u, dt, k1)
  k5 <- f(unew, p, t + dt)
  if (!isTRUE(ctx$adaptive)) return(list(u = unew, fsal = k5, EEst = NULL, nf = 4L))
  est <- function(s) {
    p1 <- (1 - s) * u + s * unew +
      s * (s - 1) * ((1 - 2 * s) * (unew - u) + (s - 1) * dt * k1 + s * dt * k5)
    pp1 <- k1 + s * (-4 * dt * k1 - 2 * dt * k5 - 6 * u +
                       s * (3 * dt * k1 + 3 * dt * k5 + 6 * u - 6 * unew) + 6 * unew) / dt
    .err_norm(dt * (f(p1, p, .fma(s, dt, t)) - pp1), u, unew, ctx$abstol, ctx$reltol)
  }
  e1 <- est(1 / 2 - sqrt(3) / 6)
  e2 <- est(1 / 2 + sqrt(3) / 6)
  list(u = unew, fsal = k5, EEst = 2.1342 * max(e1, e2), nf = 6L)
}

.abm_start_step <- function(ctx, t, u, dt, k1, bs3) {
  if (!bs3) return(.rk4_shampine(ctx, t, u, dt, k1))
  r <- .bs3_step(ctx, t, u, dt, list(fsal = k1))
  EEst <- if (isTRUE(ctx$adaptive)) .err_norm(r$err, u, r$u, ctx$abstol, ctx$reltol) else NULL
  list(u = r$u, fsal = r$cache$fsal, EEst = EEst, nf = 3L)
}

# ---------------------------------------------------------------------------
# VCAB3, VCAB4, VCAB5 and VCABM3, VCABM4, VCABM5: n - 1 starting steps
# (Bogacki-Shampine for n = 3, else Runge-Kutta), then the Adams method of
# order n in divided difference form.

.vc_step <- function(n, moulton) function(ctx, t, u, dt, cache) {
  k1 <- cache$fsal
  if (isTRUE(ctx$derivative_discontinuity)) cache$step <- 1L
  k <- cache$step
  cache$undo <- c(list(fsal = k1), if (k < n) .abm_startup_undo(cache, k))
  dts <- cache$dts
  tmp <- dts[n]
  if (k > 1L) dts[2:k] <- dts[seq_len(k - 1L)]
  dts[1L] <- dt
  cache$dts <- dts
  if (k < n) cache$step <- k + 1L
  .phi_and_phistar(cache, k1, k)
  if (k < n) {
    r <- .abm_start_step(ctx, t, u, dt, k1, n == 3L)
    .abm_swap(cache)
    cache$fsal <- r$fsal
    return(list(u = r$u, EEst = r$EEst %||% 0, k = list(k1, r$fsal), cache = cache, nf = r$nf))
  }
  nf <- 0L
  if (moulton) {
    .g_coefs(cache, k + 1L)
    g <- cache$g; phis <- cache$phistar_n
    unew <- u
    for (i in seq_len(k - 1L)) unew <- unew + g[i] * phis[[i]]
    du <- ctx$f(unew, ctx$p, t + dt)
    nf <- 1L
    .phi_np1(cache, du, k + 1L)
    phi1 <- cache$phi_np1
    unew <- unew + g[k] * phi1[[k]]
    utilde <- if (isTRUE(ctx$adaptive)) (g[k + 1L] - g[k]) * phi1[[k + 1L]]
  } else {
    .g_coefs(cache, k)
    g <- cache$g; phis <- cache$phistar_n
    unew <- u
    for (i in seq_len(k)) unew <- unew + g[i] * phis[[i]]
    utilde <- if (isTRUE(ctx$adaptive)) g[k] * phis[[k]]
  }
  EEst <- 0
  if (isTRUE(ctx$adaptive)) {
    EEst <- .err_norm(utilde, u, unew, ctx$abstol, ctx$reltol)
    if (EEst > 1) {
      cache$dts <- c(cache$dts[-1L], tmp)
      return(list(u = unew, EEst = EEst, k = list(k1, k1), cache = cache, nf = nf))
    }
  }
  fnew <- ctx$f(unew, ctx$p, t + dt)
  .abm_swap(cache)
  cache$fsal <- fnew
  list(u = unew, EEst = EEst, k = list(k1, fnew), cache = cache, nf = nf + 1L)
}

.vc_alg <- function(name, n, moulton, thread) {
  # VCAB3's initialize! does not count its evaluation of f (as in Julia)
  nf0 <- if (n == 3L && !moulton) 0L else 1L
  init <- function(ctx, t, u) {
    cache <- .abm_cache(ctx, t, u, n, n + moulton, n, n)
    cache$nf0 <- nf0
    cache
  }
  .new_alg(name, n, step = .vc_step(n, moulton), init = init, on_reject = .abm_on_reject,
           on_event = .abm_on_event)
}

#' @rdname adams_methods
#' @export
VCAB3 <- function(thread = NULL) .vc_alg("VCAB3", 3L, FALSE)

#' @rdname adams_methods
#' @export
VCAB4 <- function(thread = NULL) .vc_alg("VCAB4", 4L, FALSE)

#' @rdname adams_methods
#' @export
VCAB5 <- function(thread = NULL) .vc_alg("VCAB5", 5L, FALSE)

#' @rdname adams_methods
#' @export
VCABM3 <- function(thread = NULL) .vc_alg("VCABM3", 3L, TRUE)

#' @rdname adams_methods
#' @export
VCABM4 <- function(thread = NULL) .vc_alg("VCABM4", 4L, TRUE)

#' @rdname adams_methods
#' @export
VCABM5 <- function(thread = NULL) .vc_alg("VCABM5", 5L, TRUE)

# ---------------------------------------------------------------------------
# VCABM: variable order 1 to 12 (VCABMConstantCache). The order rises to 3
# over the first steps, then follows the error estimates of the neighboring
# orders; the I controller uses the current order.

# The corrections described in the header.
.vcabm_restart <- function(ctx, cache) {
  if (isTRUE(ctx$derivative_discontinuity)) cache$order <- 1L
  invisible(cache)
}

.vcabm_base <- function(k, uprev, upred) if (k == 1L) uprev else upred

.vcabm_fixed_fsal <- function(ctx, cache, u, t) list(fsal = ctx$f(u, ctx$p, t), nf = 1L)

.vcabm_step <- function(ctx, t, u, dt, cache) {
  f <- ctx$f; p <- ctx$p
  k1 <- cache$fsal
  cache$undo <- list(fsal = k1)
  .vcabm_restart(ctx, cache)
  step <- ctx$iter
  k <- order <- cache$order
  dts <- cache$dts
  tmp <- dts[13L]
  cache$dts <- c(dt, dts[1:12])
  .phi_and_phistar(cache, k1, k)
  .g_coefs(cache, k + 1L)
  g <- cache$g; phis <- cache$phistar_n
  upred <- g[1L] * phis[[1L]] + u
  for (i in seq_len(max(k - 2L, 0L)) + 1L) upred <- g[i] * phis[[i]] + upred
  du <- f(upred, p, t + dt)
  nf <- 1L
  .phi_np1(cache, du, k + 1L)
  phi1 <- cache$phi_np1
  unew <- g[k] * phi1[[k]] + .vcabm_base(k, u, upred)
  EEst <- 0
  if (isTRUE(ctx$adaptive)) {
    utilde <- (g[k + 1L] - g[k]) * phi1[[k + 1L]]
    EEst <- .err_norm(utilde, u, unew, ctx$abstol, ctx$reltol)
    if (EEst > 1) {
      cache$dts <- c(cache$dts[-1L], tmp)
      return(list(u = unew, EEst = EEst, k = list(k1, k1), cache = cache, nf = nf))
    }
    fnew <- f(unew, p, t + dt)
    nf <- nf + 1L
    if (step <= 4L || order < 3L) {
      cache$order <- min(order + 1L, 3L)
    } else {
      utm2 <- (g[k - 1L] - g[k - 2L]) * phi1[[k - 1L]]
      utm1 <- (g[k] - g[k - 1L]) * phi1[[k]]
      .expand_phi(cache, k + 1L)
      .phi_np1(cache, fnew, k + 2L)
      utp1 <- dt * .ABM_GAMMASTAR[k + 2L] * cache$phi_np1[[k + 2L]]
      errm2 <- .err_norm(utm2, u, unew, ctx$abstol, ctx$reltol)
      errm1 <- .err_norm(utm1, u, unew, ctx$abstol, ctx$reltol)
      errp1 <- .err_norm(utp1, u, unew, ctx$abstol, ctx$reltol)
      if (max(errm2, errm1) <= EEst) {
        cache$order <- order - 1L
      } else if (errp1 < EEst) {
        cache$order <- min(order + 1L, 12L)
        # keeps the step size in the next step
        EEst <- 1
      }
    }
  } else {
    r <- .vcabm_fixed_fsal(ctx, cache, unew, t + dt)
    fnew <- r$fsal
    nf <- nf + r$nf
  }
  .abm_swap(cache)
  cache$fsal <- fnew
  list(u = unew, EEst = EEst, k = list(k1, fnew), cache = cache, nf = nf)
}

#' @rdname adams_methods
#' @export
VCABM <- function(thread = NULL) {
  init <- function(ctx, t, u) .abm_cache(ctx, t, u, 13L, 13L, 13L, 1L)
  .new_alg("VCABM", 1L, step = .vcabm_step, init = init, controller = "I",
           on_reject = .abm_on_reject, on_event = .abm_on_event)
}
