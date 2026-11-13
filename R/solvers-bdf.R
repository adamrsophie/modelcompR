# Backward differentiation and numerical differentiation formulas of
# OrdinaryDiffEqBDF (the out of place ConstantCache steps):
#   bdf_perform_step.jl   ABDF2, SBDF (IMEXEuler, IMEXEulerARK, SBDF2 to
#                         SBDF4), QNDF1, QNDF2, QNDF, MEBDF2
#   dae_perform_step.jl   DABDF2
#   bdf_caches.jl, dae_caches.jl, bdf_utils.jl (R!, U!, update_D!, gamma_k,
#                         _sbdf_dt_changed, error_constant), controllers.jl
#                         (the QNDF order and step size controller),
#                         bdf_interpolants.jl (QNDF), alg_utils.jl, algorithms.jl
# NordsieckBDF and DNordsieckBDF are in R/solvers-bdf-nordsieck.R. The
# reference values of tests/testthat/test-bdf.R come from tools/bdf-probe.jl.
#
# Deliberate differences from the installed Julia code:
# * SBDF: after an event the derivatives f1(u) and f2(u) of the new state
#   are evaluated (.sbdf_on_event); OrdinaryDiffEqBDF keeps du1 and du2 of
#   the state before the event and only evaluates fsalfirst = f(u).
# * QNDF: a step retried after a Newton failure starts from the backward
#   differences of the last accepted step (.qndf_newton_fail); Julia keeps
#   the differences the failed attempt rescaled and rescales them again.
# * DABDF2: the dense output is the cubic Hermite interpolant from the
#   derivatives at the start and the end of the step (.dabdf2_interp);
#   OrdinaryDiffEqBDF stores k = [du(t1), du(t0)], so its interpolant uses
#   them in the reversed order (and du(t0) of the rejected attempt after a
#   rejection).

# ---------------------------------------------------------------------------
# Shared pieces

# Julia's default kappa coefficients are Rationals: 1 - kappa and the
# products with other Rationals are exact and rounded once when they meet a
# Float64. A numeric kappa that is the double nearest to num / den with a
# denominator up to 10^6 is taken as that fraction; any other value is used
# as a Float64, as Julia does with a Float64 kappa.
.bdfk_rational <- function(x) {
  if (!is.finite(x)) return(NULL)
  if (x == round(x)) return(c(x, 1))
  h0 <- 0; h1 <- 1; k0 <- 1; k1 <- 0
  r <- x
  for (i in seq_len(40L)) {
    a <- floor(r)
    h <- a * h1 + h0; k <- a * k1 + k0
    if (k > 1e6) return(NULL)
    if (h / k == x) return(c(h, k))
    h0 <- h1; h1 <- h; k0 <- k1; k1 <- k
    if (r == a) return(NULL)
    r <- 1 / (r - a)
  }
  NULL
}

# The Float64 values of 1 - kappa and of kappa * m / d (m / d a Rational
# factor) as Julia rounds them.
.bdfk_one_minus <- function(kappa) {
  r <- .bdfk_rational(kappa)
  if (is.null(r)) 1 - kappa else (r[2L] - r[1L]) / r[2L]
}

.bdfk_times <- function(kappa, m, d) {
  r <- .bdfk_rational(kappa)
  if (is.null(r)) kappa * (m / d) else (r[1L] * m) / (r[2L] * d)
}

# The statistics of a Newton solve and a failed step (a Newton failure is
# not counted as a rejection: force_stepfail).
.bdfk_counts <- function(nl, cache) {
  list(nf = nl$nf, njacs = nl$njacs, nsolve = nl$nsolve, nw = nl$nw, cache = cache)
}

.bdfk_fail <- function(out, u) c(out, list(fail = TRUE, u = u, nreject = 0L))

# The step size controller hooks of the "BDF" controller type for methods
# whose controller acts after the step: the error estimate and whether this
# is the first step are kept for the accept hook.
.bdfk_stepsize <- function(cache, EEst, co, first) {
  cache$EEst_cur <- EEst
  cache$first <- first
  1
}

# ---------------------------------------------------------------------------
# ABDF2: the adaptive second order BDF with a first step of implicit Euler
# (the ESDIRK :ie_dd2 tableau of OrdinaryDiffEqSDIRK).

.abdf2_init <- function(ctx, t, u) {
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  cache$fsal <- ctx$f(u, ctx$p, t)
  # fsalfirstprev starts as the rate prototype (u0)
  cache$fsalprev <- u
  cache$dtprev <- 1
  cache
}

# The local error estimate of ABDF2 and DABDF2 from the derivatives at the
# last three points.
.abdf2_est <- function(ctx, cache, t, dt, u, unew, fnew) {
  dtp <- cache$dtprev
  r <- dt / dtp
  tmp <- r * cache$fsalprev + (-(1 + r) * cache$fsal + fnew)
  est <- (dtp + dt) / 6 * tmp
  .err_norm(est, u, unew, ctx$abstol, ctx$reltol)
}

.abdf2_step <- function(ctx, t, u, dt, cache, extrapolant = "linear") {
  fprev <- cache$fsal
  if (ctx$iter == 1L && !isTRUE(ctx$derivative_discontinuity)) {
    # the implicit Euler start: z = dt f(u) as the guess, DIRK form
    cache$dtprev <- dt
    nl <- .nlsolve(cache$nl, ctx, t, dt, dt * fprev, u, 1, 1, t + dt, u, method = "dirk")
    out <- .bdfk_counts(nl, cache)
    if (!nl$ok) return(.bdfk_fail(out, u))
    unew <- u + nl$z
    fnew <- ctx$f(unew, ctx$p, t + dt)
    out$nf <- out$nf + 1L
    EEst <- 1
    if (ctx$adaptive) EEst <- EEst + .alg_penalty(ctx, fnew)
    cache$fsalprev <- fprev
    if (!ctx$adaptive || EEst <= 1) cache$fsal <- fnew
    return(c(out, list(u = unew, EEst = EEst, k = list(.mm_deriv(ctx, fprev), .mm_deriv(ctx, fnew)))))
  }
  rho <- dt / cache$dtprev
  b0 <- 2 / 3
  b1 <- -(rho - 1) / 3
  ah <- rho * rho / 3
  a1 <- 1 + ah
  a2 <- -ah
  z0 <- if (identical(extrapolant, "linear")) .fma(dt, fprev, u) else u
  if (is.null(ctx$mass)) {
    tmp <- .fma(dt * b1, fprev, .fma(a2, ctx$uprev2, a1 * u)) / (dt * b0)
  } else {
    mt <- .mass_mul(ctx$mass, .fma(a2, ctx$uprev2, a1 * u))
    tmp <- .fma(dt * b1, fprev, mt) / (dt * b0)
  }
  nl <- .nlsolve(cache$nl, ctx, t, dt, z0, tmp, b0, 1, t + dt, u)
  out <- .bdfk_counts(nl, cache)
  if (!nl$ok) return(.bdfk_fail(out, u))
  unew <- nl$z
  fnew <- ctx$f(unew, ctx$p, t + dt)
  out$nf <- out$nf + 1L
  EEst <- 1
  if (ctx$adaptive) EEst <- .abdf2_est(ctx, cache, t, dt, u, unew, fnew)
  if (EEst < 1) {
    cache$fsalprev <- fprev
    cache$dtprev <- dt
  }
  if (!ctx$adaptive || EEst <= 1) cache$fsal <- fnew
  c(out, list(u = unew, EEst = EEst, k = list(.mm_deriv(ctx, fprev), .mm_deriv(ctx, fnew))))
}

# ---------------------------------------------------------------------------
# QNDF1 and QNDF2: the first and second order numerical differentiation
# formulas with the quasi constant step size implementation.

.qndf12_init <- function(ctx, t, u) {
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  cache$fsal <- ctx$f(u, ctx$p, t)
  n <- length(u)
  cache$uprev2 <- u
  cache$D1 <- numeric(n)
  cache$D2 <- numeric(n)
  cache$dtprev <- 0
  cache
}

.qndf1_step <- function(ctx, t, u, dt, cache, kappa) {
  fprev <- cache$fsal
  D <- cache$D1
  if (ctx$iter > 1L && cache$dtprev != 0) {
    rho <- dt / cache$dtprev
    D <- u - cache$uprev2
    # R * U of the 1 by 1 matrices is (-rho) (-1) = rho
    if (rho != 1) D <- D * rho
    cache$D1 <- D
    a0 <- .bdfk_one_minus(kappa)
    a1 <- .bdfk_one_minus(2 * kappa)
    if (!is.null(r <- .bdfk_rational(kappa))) a1 <- (r[2L] - 2 * r[1L]) / r[2L]
    a2 <- kappa
    ecoef <- kappa + 1 / 2
  } else {
    a0 <- 1; a1 <- 1; a2 <- 0; ecoef <- 0 + 1 / 2
  }
  uprev2 <- u - D
  z0 <- u + D
  tmp <- a1 * u + a2 * uprev2
  tmp <- .mass_mul(ctx$mass, tmp) / dt
  nl <- .nlsolve(cache$nl, ctx, t, dt, z0, tmp, 1, a0, t + dt, u)
  out <- .bdfk_counts(nl, cache)
  if (!nl$ok) return(.bdfk_fail(out, u))
  unew <- nl$z
  EEst <- 1
  if (ctx$adaptive && ctx$success_iter > 0L) {
    d21 <- unew - u
    d22 <- d21 - D
    EEst <- .err_norm(ecoef * d22, u, unew, ctx$abstol, ctx$reltol)
  }
  if (ctx$adaptive && EEst > 1) return(c(out, list(u = unew, EEst = EEst, k = list(fprev, fprev))))
  cache$dtprev <- dt
  cache$uprev2 <- u
  fnew <- ctx$f(unew, ctx$p, t + dt)
  out$nf <- out$nf + 1L
  cache$fsal <- fnew
  c(out, list(u = unew, EEst = EEst, k = list(.mm_deriv(ctx, fprev), .mm_deriv(ctx, fnew))))
}

.qndf2_step <- function(ctx, t, u, dt, cache, kappa) {
  fprev <- cache$fsal
  startup <- ctx$iter <= 2L || cache$dtprev == 0
  if (startup) {
    g2 <- 1; b0 <- 1; ecoef <- 0 + 1 / 3
  } else {
    g2 <- 1.5
    # beta0 = inv((1 - kappa) * 3 // 2) and kappa * 3 // 2, exact for a
    # Rational kappa
    r <- .bdfk_rational(kappa)
    b0 <- if (is.null(r)) 1 / ((1 - kappa) * 1.5) else (2 * r[2L]) / (3 * (r[2L] - r[1L]))
    ecoef <- .bdfk_times(kappa, 3, 2) + 1 / 3
  }
  D1 <- cache$D1; D2 <- cache$D2
  if (!startup && dt != cache$dtprev) {
    rho <- dt / cache$dtprev
    # R!(2, rho) and R * U (matmul2x2, without fused operations)
    R11 <- -1 * rho; R12 <- -2 * rho
    R21 <- R11 * (1 - 1 * rho) / 2
    R22 <- R12 * (1 - 2 * rho) / 2
    U11 <- -1; U12 <- -2; U21 <- -1 * (1 - 1) / 2; U22 <- -2 * (1 - 2) / 2
    RU11 <- R11 * U11 + R12 * U21
    RU12 <- R11 * U12 + R12 * U22
    RU21 <- R21 * U11 + R22 * U21
    RU22 <- R21 * U12 + R22 * U22
    d1 <- D1 * RU11 + D2 * RU21
    d2 <- D1 * RU12 + D2 * RU22
  } else {
    d1 <- D1; d2 <- D2
  }
  u0 <- u + d1 + d2
  phi <- (1 * d1 + g2 * d2) * b0
  tmp <- .mass_mul(ctx$mass, u0 - phi) / (dt * b0)
  nl <- .nlsolve(cache$nl, ctx, t, dt, u0, tmp, b0, 1, t + dt, u)
  out <- .bdfk_counts(nl, cache)
  if (!nl$ok) return(.bdfk_fail(out, u))
  unew <- nl$z
  EEst <- 1
  if (ctx$adaptive) {
    if (ctx$success_iter == 1L) {
      utilde <- (unew - u) - ((u - cache$uprev2) * dt / cache$dtprev)
      EEst <- .err_norm(utilde, u, unew, ctx$abstol, ctx$reltol)
    } else if (ctx$success_iter > 1L) {
      e1 <- unew - u
      e2 <- e1 - d1
      e3 <- e2 - d2
      EEst <- .err_norm(ecoef * e3, u, unew, ctx$abstol, ctx$reltol)
    }
  }
  if (ctx$adaptive && EEst > 1) return(c(out, list(u = unew, EEst = EEst, k = list(fprev, fprev))))
  delta <- unew - u
  cache$D1 <- delta
  cache$D2 <- if (ctx$success_iter == 0L) numeric(length(u)) else delta - d1
  cache$uprev2 <- u
  cache$dtprev <- dt
  fnew <- ctx$f(unew, ctx$p, t + dt)
  out$nf <- out$nf + 1L
  cache$fsal <- fnew
  c(out, list(u = unew, EEst = EEst, k = list(.mm_deriv(ctx, fprev), .mm_deriv(ctx, fnew))))
}

# ---------------------------------------------------------------------------
# QNDF: variable order (1 to max_order) NDF on backward differences, with
# the order and step size controller of Alberdi Celaya et al.

.qndf_init <- function(ctx, t, u, max_order, kappa) {
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  cache$fsal <- ctx$f(u, ctx$p, t)
  mo <- max_order
  n <- length(u)
  cache$max_order <- mo
  cache$D <- rep(list(numeric(n)), mo + 2L)
  cache$prevD <- cache$D
  cache$order <- 1L
  cache$prevorder <- 1L
  cache$dtprev <- 1
  cache$nconsteps <- 0L
  cache$consfailcnt <- 0L
  cache$EEst1 <- 1
  cache$EEst2 <- 1
  # gamma_k = sum(1 / j), summed in Float64 from the left
  gam <- numeric(mo)
  s <- 0
  for (k in seq_len(mo)) {
    s <- if (k == 1L) 1 else s + 1 / k
    gam[k] <- s
  }
  cache$gam <- gam
  # 1 - kappa_k and kappa_k as Julia rounds them
  cache$omk <- vapply(seq_len(mo), function(k) .bdfk_one_minus(kappa[k]), 0)
  cache$kf <- vapply(seq_len(mo), function(k) {
    r <- .bdfk_rational(kappa[k])
    if (is.null(r)) kappa[k] else r[1L] / r[2L]
  }, 0)
  U <- matrix(0, mo, mo)
  for (r in seq_len(mo)) {
    U[1L, r] <- -r
    if (mo >= 2L) for (j in 2:mo) U[j, r] <- U[j - 1L, r] * ((j - 1) - r) / j
  }
  cache$U <- U
  cache
}

# error_constant(integrator, k) = kappa_k gamma_k + 1 / (k + 1)
.qndf_ec <- function(cache, k) cache$kf[k] * cache$gam[k] + 1 / (k + 1)

# calc_R! and mul!(RU, R, U): matrices up to 3 by 3 are multiplied without
# fused operations (matmul2x2!, matmul3x3!), larger ones by BLAS as a fused
# multiply add chain over all max_order columns.
.qndf_RU <- function(cache, rho, k) {
  mo <- cache$max_order
  U <- cache$U
  R <- matrix(0, mo, mo)
  for (r in seq_len(k)) {
    R[1L, r] <- -r * rho
    if (k >= 2L) for (j in 2:k) R[j, r] <- R[j - 1L, r] * ((j - 1) - r * rho) / j
  }
  RU <- matrix(0, mo, mo)
  for (i in seq_len(mo)) for (j in seq_len(mo)) {
    if (mo <= 3L) {
      s <- R[i, 1L] * U[1L, j]
      if (mo >= 2L) for (l in 2:mo) s <- s + R[i, l] * U[l, j]
    } else {
      s <- 0
      for (l in seq_len(mo)) s <- .fma(R[i, l], U[l, j], s)
    }
    RU[i, j] <- s
  }
  RU
}

.qndf_step <- function(ctx, t, u, dt, cache) {
  disc <- isTRUE(ctx$derivative_discontinuity) || isTRUE(cache$disc)
  cache$disc <- FALSE
  if (disc) ctx$derivative_discontinuity <- TRUE
  dtprev <- cache$dtprev
  k <- cache$order
  n <- length(u)
  if (disc) {
    # the order 1 restart of the step; cache$order is kept
    dtprev <- 1
    k <- 1L
    cache$nconsteps <- 0L
    cache$consfailcnt <- 0L
    cache$D <- rep(list(numeric(n)), cache$max_order + 2L)
    cache$prevD <- cache$D
  }
  if (cache$consfailcnt > 0L || isTRUE(cache$restore_D)) cache$D <- cache$prevD
  cache$restore_D <- FALSE
  D <- cache$D
  if (dt != dtprev || cache$prevorder != k) {
    rho <- dt / dtprev
    cache$nconsteps <- 0L
    RU <- .qndf_RU(cache, rho, k)
    Dn <- vector("list", k)
    for (j in seq_len(k)) {
      s <- D[[1L]] * RU[1L, j]
      if (k >= 2L) for (i in 2:k) s <- s + D[[i]] * RU[i, j]
      Dn[[j]] <- s
    }
    D[seq_len(k)] <- Dn
    cache$D <- D
  }
  gam <- cache$gam
  b0 <- 1 / (cache$omk[k] * gam[k])
  u0 <- D[[1L]]
  if (k >= 2L) for (i in 2:k) u0 <- u0 + D[[i]]
  u0 <- u0 + u
  phi <- numeric(n)
  for (i in seq_len(k)) phi <- phi + gam[i] * D[[i]]
  tmp <- .mass_mul(ctx$mass, (u0 / b0 - phi) / dt)
  ctx$newton_scale <- .qndf_ec(cache, 1L)
  nl <- .nlsolve(cache$nl, ctx, t, dt, u0, tmp, b0, 1, t + dt, u)
  out <- .bdfk_counts(nl, cache)
  if (!nl$ok) return(.bdfk_fail(out, u))
  unew <- nl$z
  dd <- unew - u0
  # update_D!
  D[[k + 2L]] <- dd - D[[k + 1L]]
  D[[k + 1L]] <- dd
  for (i in k:1) D[[i]] <- D[[i]] + D[[i + 1L]]
  cache$D <- D
  EEst <- 1
  if (ctx$adaptive) {
    nrm <- function(x) .err_norm(x, u, unew, ctx$abstol, ctx$reltol)
    e <- if (cache$consfailcnt > 1L && !is.null(ctx$mass)) .mass_mul(ctx$mass, dd) else dd
    EEst <- .qndf_ec(cache, k) * nrm(e)
    if (k > 1L) cache$EEst1 <- .qndf_ec(cache, k - 1L) * nrm(D[[k]])
    if (k < cache$max_order) cache$EEst2 <- .qndf_ec(cache, k + 1L) * nrm(D[[k + 2L]])
  }
  if (!ctx$adaptive || EEst <= 1) {
    cache$prevD <- D
    cache$dtprev <- dt
    cache$prevorder <- k
    # opts.dense, put in ctx by solve.ODEProblem (no isTRUE: a driver
    # without that field fails here instead of skipping the f call)
    if (ctx$dense) {
      cache$fsal <- ctx$f(unew, ctx$p, t + dt)
      out$nf <- out$nf + 1L
    }
  }
  kk <- c(D[seq_len(k)], rep(list(numeric(n)), cache$max_order - k))
  c(out, list(u = unew, EEst = EEst, k = kk))
}

# Newton backward difference interpolation through the stored differences
.qndf_interp <- function(theta, dt, y0, y1, k) {
  s <- theta - 1
  phi <- s
  out <- .fma(phi, k[[1L]], y1)
  if (length(k) >= 2L) for (j in 2:length(k)) {
    phi <- phi * ((s + j - 1) / j)
    out <- .fma(phi, k[[j]], out)
  }
  out
}

# step_accept_controller! of QNDF
.qndf_accept <- function(cache, q, dt, co) {
  cache$consfailcnt <- 0L
  cache$nconsteps <- cache$nconsteps + 1L
  est <- cache$EEst_cur
  if (est == 0) return(dt * (if (isTRUE(cache$first)) co$qmax_first else co$qmax))
  h <- abs(dt)
  k <- cache$order
  prefer_const <- cache$nconsteps < k + 2L
  z <- co$gamma * est^(1 / (k + 1))
  hk <- if (z <= 0.1) 10 * h else (1 / z) * h
  hkm1 <- 0
  hkp1 <- 0
  hn <- hk
  kn <- k
  if (k > 1L) {
    zk1 <- 1.3 * cache$EEst1^(1 / k)
    if (zk1 <= 0.1) hkm1 <- 10 * h else if (zk1 <= 1.3) hkm1 <- (1 / zk1) * h
    if (hkm1 > hk) {
      hn <- hkm1
      kn <- k - 1L
    }
  }
  if (k < cache$max_order) {
    zkp1 <- 1.4 * cache$EEst2^(1 / (k + 2))
    if (zkp1 <= 0.1) hkp1 <- 10 * h else if (0.1 < zkp1 && zkp1 <= 1.4) hkp1 <- (1 / zkp1) * h
    if (hkp1 > hn) {
      hn <- hkp1
      kn <- k + 1L
    }
  }
  cache$order <- kn
  q <- h / hn
  if (prefer_const && 0.6 < q && q < 1.2) h
  else if (q <= co$qsteady_max && q >= co$qsteady_min) h
  else h / q
}

# bdf_step_reject_controller! with EEst1
.qndf_reject <- function(cache, EEst, dt, co) {
  k <- cache$order
  h <- abs(dt)
  cache$consfailcnt <- cache$consfailcnt + 1L
  cache$nconsteps <- 0L
  if (cache$consfailcnt > 1L) h <- h / 2
  z <- co$gamma * EEst^(1 / (k + 1))
  hk <- if (z <= 10) (1 / z) * h else 0.1 * h
  hn <- hk
  kn <- k
  if (k > 1L) {
    zk1 <- 1.3 * cache$EEst1^(1 / k)
    hk1 <- if (zk1 <= 10) (1 / zk1) * h else 0.1 * h
    if (cache$consfailcnt > 2L || abs(hk1) > abs(hk)) {
      hn <- min(abs(h), abs(hk1))
      kn <- k - 1L
    }
  }
  # derivative_discontinuity!(integrator, true): the next attempt restarts
  if (kn == 1L && cache$consfailcnt > 3L) cache$disc <- TRUE
  cache$order <- kn
  hn
}

# post_newton_controller!: dt / failfactor. The retry restarts from the
# differences of the last accepted step (see the header).
.qndf_newton_fail <- function(cache, dt, co) {
  cache$restore_D <- TRUE
  dt / co$failfactor
}

# ---------------------------------------------------------------------------
# MEBDF2: the modified extended BDF of order 2 (fixed step), three implicit
# Euler type solves per step.

.mebdf2_step <- function(ctx, t, u, dt, cache, extrapolant = "constant") {
  fprev <- cache$fsal
  nlsv <- cache$nl
  z <- if (identical(extrapolant, "linear")) dt * fprev else numeric(length(u))
  cnt <- list(nf = 0L, njacs = 0L, nsolve = 0L, nw = 0L)
  add <- function(x) for (s in names(cnt)) cnt[[s]] <<- cnt[[s]] + as.integer(x[[s]] %||% 0L)
  r <- .nlsolve(nlsv, ctx, t, dt, z, u, 1, 1, t + dt, u, method = "dirk", eta_first = FALSE)
  add(r)
  if (!r$ok) return(.bdfk_fail(c(cnt, list(cache = cache)), u))
  z1 <- u + r$z
  r <- .nlsolve(nlsv, ctx, t, dt, r$z, z1, 1, 1, t + 2 * dt, u, method = "dirk", firststage = FALSE, eta_first = FALSE)
  add(r)
  if (!r$ok) return(.bdfk_fail(c(cnt, list(cache = cache)), u))
  z2 <- z1 + r$z
  tmp2 <- (0.5 * u + z1) - 0.5 * z2
  r <- .nlsolve(nlsv, ctx, t, dt, r$z, tmp2, 1, 1, t + dt, u, method = "dirk", firststage = FALSE, eta_first = FALSE)
  add(r)
  if (!r$ok) return(.bdfk_fail(c(cnt, list(cache = cache)), u))
  unew <- tmp2 + r$z
  fnew <- ctx$f(unew, ctx$p, t + dt)
  cnt$nf <- cnt$nf + 1L
  cache$fsal <- fnew
  c(cnt, list(cache = cache, u = unew, EEst = 1, k = list(.mm_deriv(ctx, fprev), .mm_deriv(ctx, fnew))))
}

# ---------------------------------------------------------------------------
# SBDF: the implicit explicit BDF methods of orders 1 to 4 for
# SplitODEProblem (f1 implicit, f2 explicit).

.sbdf_init <- function(ctx, t, u) {
  if (is.null(ctx$f1) || is.null(ctx$f2))
    stop("SBDF methods integrate a SplitODEProblem (f1 implicit, f2 explicit)", call. = FALSE)
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  cache$du1 <- ctx$f1(u, ctx$p, t)
  cache$du2 <- ctx$f2(u, ctx$p, t)
  cache$fsal <- cache$du1 + cache$du2
  # the histories start as the rate and state prototypes (u0)
  cache$uprev2 <- u; cache$uprev3 <- u; cache$uprev4 <- u
  cache$k1 <- u; cache$k2 <- u; cache$k3 <- u
  cache$dtprev <- 0
  cache$nf0 <- 1L
  cache$nf20 <- 1L
  cache
}

# _sbdf_dt_changed: a change of dt within 100 eps(t) is roundoff
.sbdf_dt_changed <- function(dt, dtprev, t) {
  if (dt == dtprev) return(FALSE)
  if (!is.finite(t)) return(TRUE)
  tdt <- t + dt
  if (!is.finite(tdt)) return(TRUE)
  m <- max(abs(t), abs(tdt))
  ulp <- if (m == 0) 2^-1074 else 2^(floor(log2(m)) - 52)
  # log2 may round up just below a power of two
  if (2^(floor(log2(m))) > m) ulp <- ulp / 2
  abs(dt - dtprev) > 100 * ulp
}

.sbdf_step <- function(ctx, t, u, dt, cache, order, ark, extrapolant = "linear") {
  f1 <- ctx$f1; f2 <- ctx$f2; p <- ctx$p
  nf2 <- 0L
  cnt <- min(order, ctx$iter)
  if (ctx$iter == 1L && !isTRUE(ctx$derivative_discontinuity)) cnt <- 1L
  if (.sbdf_dt_changed(dt, cache$dtprev, t)) cnt <- 1L
  # inv(gamma_k[cnt]) of the Rationals 1, 3 // 2, 11 // 6, 25 // 12
  g <- c(1, 2 / 3, 6 / 11, 12 / 25)[cnt]
  du1 <- cache$du1; du2 <- cache$du2
  if (ark) {
    du2 <- f2(u + dt * du1, p, t)
    nf2 <- nf2 + 1L
  }
  up2 <- cache$uprev2; up3 <- cache$uprev3; up4 <- cache$uprev4
  k1 <- cache$k1; k2 <- cache$k2; k3 <- cache$k3
  tmp <- switch(cnt,
    u + dt * du2,
    g * ((2 * u - 0.5 * up2) + dt * (2 * du2 - k1)),
    g * (((3 * u - 1.5 * up2) + (1 / 3) * up3) + dt * (3 * (du2 - k1) + k2)),
    g * ((((4 * u - 3 * up2) + (4 / 3) * up3) - 0.25 * up4) + dt * (((4 * du2 - 6 * k1) + 4 * k2) - k3)))
  # the implicit part: Newton on f1 with its Jacobian; SBDF is not an
  # adaptive method
  ctx1 <- ctx
  ctx1$f <- f1
  ctx1$jac <- ctx$jac1
  ctx1$tgrad <- NULL
  z0 <- if (identical(extrapolant, "linear")) dt * du1 else numeric(length(u))
  nl <- .nlsolve(cache$nl, ctx1, t, dt, z0, tmp, g, 1, t + dt, u, method = "dirk", eta_first = FALSE)
  out <- .bdfk_counts(nl, cache)
  out$nf2 <- nf2
  if (!nl$ok) return(.bdfk_fail(out, u))
  unew <- tmp + g * nl$z
  if (order == 4L) {
    cache$uprev4 <- up3
    cache$k3 <- k2
  }
  if (order >= 3L) {
    cache$uprev3 <- up2
    cache$k2 <- k1
  }
  cache$uprev2 <- u
  cache$k1 <- du2
  cache$dtprev <- dt
  fprev <- cache$fsal
  cache$du1 <- f1(unew, p, t + dt)
  cache$du2 <- f2(unew, p, t + dt)
  out$nf <- out$nf + 1L
  out$nf2 <- out$nf2 + 1L
  cache$fsal <- cache$du1 + cache$du2
  c(out, list(u = unew, EEst = 1, k = list(fprev, cache$fsal)))
}

# After an event the derivatives belong to the new state (see the header).
.sbdf_on_event <- function(ctx, t, u, cache) {
  cache$du1 <- ctx$f1(u, ctx$p, t)
  cache$du2 <- ctx$f2(u, ctx$p, t)
  cache$fsal <- cache$du1 + cache$du2
  list(cache = cache, nf = 1L, nf2 = 1L)
}

# ---------------------------------------------------------------------------
# DABDF2: ABDF2 for fully implicit DAEs, started with DImplicitEuler.
# Neither step calls markfirststage!, so after the first solve the
# nonlinear solver is never in its first stage: an error rejection alone
# does not form a new W (firststage = FALSE).

.dabdf2_init <- function(ctx, t, u) {
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  cache$dtprev <- 1
  # fsalfirst and fsalfirstprev start as du0 (the rate prototype)
  cache$fsal <- ctx$du
  cache$fsalprev <- ctx$du
  cache$nf0 <- 0L
  cache
}

.dabdf2_step <- function(ctx, t, u, dt, cache) {
  n <- length(u)
  if (ctx$iter == 1L && !isTRUE(ctx$derivative_discontinuity)) {
    cache$dtprev <- dt
    nl <- .nlsolve(cache$nl, ctx, t, dt, numeric(n), numeric(n), 1, 1, t + dt, u, u0 = u, firststage = FALSE)
    out <- .bdfk_counts(nl, cache)
    if (!nl$ok) return(.bdfk_fail(out, u))
    unew <- u + nl$z
    du <- (unew - u) / dt
    cache$fsal <- du
    cache$fsalprev <- du
    return(c(out, list(u = unew, EEst = 1, du = du, k = list(ctx$du, du, ctx$differential_vars))))
  }
  rho <- dt / cache$dtprev
  den <- .fma(2, rho, 1)
  c1 <- rho * rho / den
  g <- (1 + rho) / den
  tmp <- -c1 * u + c1 * ctx$uprev2
  nl <- .nlsolve(cache$nl, ctx, t, dt, numeric(n), tmp, g, 1, t + dt, u, u0 = u, firststage = FALSE)
  out <- .bdfk_counts(nl, cache)
  if (!nl$ok) return(.bdfk_fail(out, u))
  z <- nl$z
  unew <- u + z
  du <- (1 * z + tmp) * (1 / (g * dt))
  EEst <- 1
  if (ctx$adaptive) EEst <- .abdf2_est(ctx, cache, t, dt, u, unew, du)
  if (EEst < 1) {
    cache$fsalprev <- cache$fsal
    cache$dtprev <- dt
  }
  if (!ctx$adaptive || EEst <= 1) cache$fsal <- du
  c(out, list(u = unew, EEst = EEst, du = du, k = list(ctx$du, du, ctx$differential_vars)))
}

# hermite_interpolant (generic_dense.jl, out of place, @muladd) with
# k = list(du(t0), du(t1), differential_vars): the algebraic components are
# interpolated linearly. With some algebraic components @muladd sums
# theta * y1 + ((1 - theta) * y0 + dv .* h).
.dabdf2_interp <- function(theta, dt, y0, y1, k) {
  if (length(k) < 2L) return((1 - theta) * y0 + theta * y1)
  h <- theta * (theta - 1) *
    ((1 - 2 * theta) * (y1 - y0) + (theta - 1) * dt * k[[1L]] + theta * dt * k[[2L]])
  dv <- k[[3L]]
  if (is.null(dv) || all(dv)) return((1 - theta) * y0 + theta * y1 + h)
  theta * y1 + ((1 - theta) * y0 + ifelse(dv, h, 0))
}

# ---------------------------------------------------------------------------
# Constructors

.bdfk_check_limiter <- function(step_limiter) {
  if (!is.null(step_limiter) && !is.function(step_limiter))
    stop("step_limiter must be a function(u, integ, p, t)", call. = FALSE)
}

.bdfk_check_extrapolant <- function(extrapolant) {
  if (!(is.character(extrapolant) && length(extrapolant) == 1L && extrapolant %in% c("linear", "constant")))
    stop("extrapolant must be \"linear\" or \"constant\"", call. = FALSE)
  extrapolant
}

#' BDF and NDF methods of OrdinaryDiffEqBDF
#'
#' Multistep methods for stiff problems, as in OrdinaryDiffEqBDF.
#' `ABDF2()` is the adaptive second order backward differentiation formula
#' (started with one implicit Euler step); `QNDF1()` and `QNDF2()` are the
#' first and second order numerical differentiation formulas (NDF) with a
#' quasi constant step size, and `QNDF()` the variable order (1 to 5) NDF of
#' MATLAB's ode15s with the order and step size controller of Alberdi
#' Celaya, Anza Aguirrezabala and Chatzipantelidis. `QBDF1()`, `QBDF2()` and
#' `QBDF()` are the same methods with `kappa = 0`, i.e. the backward
#' differentiation formulas. `MEBDF2()` is the fixed step modified extended
#' BDF of order 2. Each step is solved by a modified Newton iteration
#' (see [ImplicitEuler()]); problems with a mass matrix are supported.
#'
#' @param extrapolant The initial guess of the Newton iteration: `"linear"`
#'   (the explicit Euler predictor) or `"constant"` (the previous state).
#' @param kappa The NDF coefficient: one number for `QNDF1()` and `QNDF2()`,
#'   one per order for `QNDF()`. Julia's defaults are the Rationals
#'   `-37/200`, `-1/9`, `-823/10000`, `-83/2000` and `0`; a value that is
#'   the double nearest to a fraction with denominator up to `10^6` is
#'   treated as that fraction (exact arithmetic, as Julia's Rationals).
#' @param max_order Maximum order of `QNDF()` (1 to 5).
#' @param qmax,qsteady_min,qsteady_max Largest step growth factor and the
#'   range of step ratios that keep the step size (`QNDF()`).
#' @param step_limiter A function `limiter(u, integ, p, t)` returning the
#'   state after each accepted step.
#' @param autodiff As for [ImplicitEuler()].
#' @return An algorithm object.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -50 * (u - t * t), 1, c(0, 1))
#' sol <- solve(prob, QNDF(), reltol = 1e-6)
#' sol(1)
#' solve(prob, ABDF2())$retcode
#' sol2 <- solve(prob, MEBDF2(), dt = 1 / 50)
#' sol2(1)
#' @name bdf_methods
NULL

#' @rdname bdf_methods
#' @export
ABDF2 <- function(extrapolant = "linear", step_limiter = NULL, autodiff = NULL) .with_autodiff(autodiff, {
  .bdfk_check_limiter(step_limiter)
  ex <- .bdfk_check_extrapolant(extrapolant)
  .new_alg("ABDF2", 2L, step = function(ctx, t, u, dt, cache) .abdf2_step(ctx, t, u, dt, cache, ex),
           interp = .hermite_dae_interp, init = .abdf2_init, on_event = .implicit_on_event,
           stiff = TRUE, mass_matrix = TRUE, controller = "I", qsteady_max = 6 / 5,
           on_reject = .nl_on_reject, step_limiter = step_limiter)
})

.qndf12_alg <- function(name, order, kappa, step_limiter, step) {
  .bdfk_check_limiter(step_limiter)
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa))
    stop("kappa must be a single number", call. = FALSE)
  .new_alg(name, order, step = function(ctx, t, u, dt, cache) step(ctx, t, u, dt, cache, kappa),
           interp = .hermite_dae_interp, init = .qndf12_init, on_event = .implicit_on_event,
           stiff = TRUE, mass_matrix = TRUE, controller = "I", qsteady_min = 9 / 10,
           qsteady_max = 2, on_reject = .nl_on_reject, step_limiter = step_limiter, kappa = kappa)
}

#' @rdname bdf_methods
#' @export
QNDF1 <- function(kappa = -37 / 200, step_limiter = NULL, autodiff = NULL) .with_autodiff(autodiff, {
  .qndf12_alg("QNDF1", 1L, kappa, step_limiter, .qndf1_step)
})

#' @rdname bdf_methods
#' @export
QNDF2 <- function(kappa = -1 / 9, step_limiter = NULL, autodiff = NULL) .with_autodiff(autodiff, {
  .qndf12_alg("QNDF2", 2L, kappa, step_limiter, .qndf2_step)
})

#' @rdname bdf_methods
#' @export
QBDF1 <- function(step_limiter = NULL, autodiff = NULL) .with_autodiff(autodiff, {
  .qndf12_alg("QBDF1", 1L, 0, step_limiter, .qndf1_step)
})

#' @rdname bdf_methods
#' @export
QBDF2 <- function(step_limiter = NULL, autodiff = NULL) .with_autodiff(autodiff, {
  .qndf12_alg("QBDF2", 2L, 0, step_limiter, .qndf2_step)
})

.qndf_alg <- function(name, max_order, kappa, qmax, qsteady_min, qsteady_max, step_limiter) {
  .bdfk_check_limiter(step_limiter)
  max_order <- as.integer(max_order)
  if (length(max_order) != 1L || is.na(max_order) || max_order < 1L || max_order > 5L)
    stop("max_order must be between 1 and 5", call. = FALSE)
  if (!is.numeric(kappa) || length(kappa) < max_order || any(!is.finite(kappa)))
    stop("kappa must give one coefficient for each order up to max_order", call. = FALSE)
  .new_alg(name, 1L, step = .qndf_step, interp = .qndf_interp,
           init = function(ctx, t, u) .qndf_init(ctx, t, u, max_order, kappa),
           on_event = .implicit_on_event, stiff = TRUE, mass_matrix = TRUE, controller = "BDF",
           qmax = qmax, qsteady_min = qsteady_min, qsteady_max = qsteady_max, gamma = 12 / 10,
           bdf_stepsize = .bdfk_stepsize, bdf_accept = .qndf_accept, bdf_reject = .qndf_reject,
           newton_fail = .qndf_newton_fail, on_reject = .nl_on_reject, max_order = max_order,
           step_limiter = step_limiter, kappa = kappa)
}

#' @rdname bdf_methods
#' @export
QNDF <- function(max_order = 5L, kappa = c(-37 / 200, -1 / 9, -823 / 10000, -83 / 2000, 0),
                 qmax = 10, qsteady_min = 1, qsteady_max = 2, step_limiter = NULL,
                 autodiff = NULL) .with_autodiff(autodiff, {
  .qndf_alg("QNDF", max_order, kappa, qmax, qsteady_min, qsteady_max, step_limiter)
})

#' @rdname bdf_methods
#' @export
QBDF <- function(max_order = 5L, qmax = 10, qsteady_min = 1, qsteady_max = 2, step_limiter = NULL,
                 autodiff = NULL) .with_autodiff(autodiff, {
  .qndf_alg("QBDF", max_order, rep(0, 5L), qmax, qsteady_min, qsteady_max, step_limiter)
})

#' @rdname bdf_methods
#' @export
MEBDF2 <- function(extrapolant = "constant", autodiff = NULL) .with_autodiff(autodiff, {
  ex <- .bdfk_check_extrapolant(extrapolant)
  .new_alg("MEBDF2", 2L, adaptive = FALSE,
           step = function(ctx, t, u, dt, cache) .mebdf2_step(ctx, t, u, dt, cache, ex),
           interp = .hermite_dae_interp, init = .implicit_init, on_event = .implicit_on_event,
           stiff = TRUE, mass_matrix = TRUE)
})

#' Implicit explicit BDF methods
#'
#' The IMEX multistep methods of Ascher, Ruuth and Wetton (OrdinaryDiffEqBDF)
#' for a [SplitODEProblem()] `u' = f1(u, p, t) + f2(u, p, t)`: `f1` (stiff)
#' is treated implicitly by a BDF formula and `f2` explicitly by
#' extrapolation. `SBDF(order)` has orders 1 to 4 and fixed steps;
#' `IMEXEuler()` is `SBDF(1)` and `SBDF2()`, `SBDF3()`, `SBDF4()` the
#' higher orders. With `ark = TRUE` the explicit term is evaluated at the
#' explicit Euler predictor (an additive Runge-Kutta form);
#' `IMEXEulerARK()` is `SBDF(1, ark = TRUE)`. The step restarts at order 1
#' when the step size changes. The Jacobian of `f1` is the `jac` of the
#' [SplitFunction()] or finite differences.
#'
#' @param order The order (1 to 4).
#' @param ark Evaluate `f2` at the explicit Euler predictor.
#' @param extrapolant The initial guess of the Newton iteration: `"linear"`
#'   or `"constant"`.
#' @param autodiff As for [ImplicitEuler()].
#' @param ... Arguments of `SBDF()` other than `order`.
#' @return An algorithm object.
#' @examples
#' f1 <- function(u, p, t) -100 * u
#' f2 <- function(u, p, t) 1 + 0 * u
#' prob <- SplitODEProblem(f1, f2, 1, c(0, 1))
#' sol <- solve(prob, SBDF2(), dt = 1 / 20)
#' sol(1)
#' solve(prob, IMEXEuler(), dt = 1 / 20)$retcode
#' @name sbdf_methods
NULL

#' @rdname sbdf_methods
#' @export
SBDF <- function(order, ark = FALSE, extrapolant = "linear", autodiff = NULL) .with_autodiff(autodiff, {
  order <- as.integer(order)
  if (length(order) != 1L || is.na(order) || order < 1L || order > 4L)
    stop("order must be between 1 and 4", call. = FALSE)
  ark <- isTRUE(ark)
  ex <- .bdfk_check_extrapolant(extrapolant)
  name <- if (order == 1L) if (ark) "IMEXEulerARK" else "IMEXEuler" else paste0("SBDF", order)
  .new_alg(name, order, adaptive = FALSE,
           step = function(ctx, t, u, dt, cache) .sbdf_step(ctx, t, u, dt, cache, order, ark, ex),
           init = .sbdf_init, on_event = .sbdf_on_event, stiff = TRUE, ark = ark, sbdf_order = order)
})

#' @rdname sbdf_methods
#' @export
SBDF2 <- function(...) SBDF(2L, ...)

#' @rdname sbdf_methods
#' @export
SBDF3 <- function(...) SBDF(3L, ...)

#' @rdname sbdf_methods
#' @export
SBDF4 <- function(...) SBDF(4L, ...)

#' @rdname sbdf_methods
#' @export
IMEXEuler <- function(...) SBDF(1L, ...)

#' @rdname sbdf_methods
#' @export
IMEXEulerARK <- function(...) SBDF(1L, ark = TRUE, ...)

#' Second order BDF for fully implicit DAEs
#'
#' `DABDF2()` is the adaptive second order backward differentiation formula
#' of OrdinaryDiffEqBDF for a fully implicit [DAEProblem()]
#' `F(du, u, p, t) = 0`: the first step is a [DImplicitEuler()] step, the
#' later ones solve `F((z + tmp) / (gamma dt), u + z, p, t) = 0` by a
#' modified Newton iteration, and the local error is estimated from the
#' derivatives of the last three steps as in [ABDF2()]. The dense output is
#' the cubic Hermite interpolant of the differential variables (linear for
#' the algebraic ones).
#'
#' @param autodiff As for [ImplicitEuler()].
#' @return An algorithm object.
#' @examples
#' f <- function(du, u, p, t) c(du[1] + u[1], u[2] - u[1])
#' dprob <- DAEProblem(f, c(-1, 0), c(1, 1), c(0, 1), differential_vars = c(TRUE, FALSE))
#' sol <- solve(dprob, DABDF2())
#' sol(0.5)
#' @export
DABDF2 <- function(autodiff = NULL) .with_autodiff(autodiff, {
  .new_alg("DABDF2", 2L, step = .dabdf2_step, interp = .dabdf2_interp, init = .dabdf2_init,
           stiff = TRUE, dae = TRUE, controller = "I", on_reject = .nl_on_reject)
})
