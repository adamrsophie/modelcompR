# Nordsieck form Adams and BDF methods of OrdinaryDiffEqNordsieck 2.3.2:
# nordsieck_perform_step.jl (the out of place ConstantCache steps),
# nordsieck_utils.jl, nordsieck_caches.jl, controllers.jl and alg_utils.jl.
# The fixed leading coefficient BDF formulas of JVODE_BDF follow cvSetBDF,
# cvSetTqBDF, cvIncreaseBDF and cvDecreaseBDF of SUNDIALS CVODE, as the Julia
# code does.
#
# Deliberate differences from the installed Julia code (tools/nordsieck-probe.jl
# runs Julia with the same corrections for the reference values):
# * AN5: a step rejected by its error estimate rewinds the Nordsieck vector
#   and restores the step size history (.an5_reject). Julia restores the
#   history but corrects the predicted vector and keeps it, so the retried
#   step starts from the vector of the rejected step.
# * AN5: a failed functional iteration also restores the step size history
#   (.an5_fail); Julia only rewinds the vector.
# * AN5: the starting step clears the corrector coefficients, as in a fresh
#   cache (.an5_start_coeffs). After a restart Julia applies the coefficients
#   of the last step before the callback to the new Nordsieck vector.
# * AN5: the dense output of the starting step uses the derivative at its end
#   (.an5_start_k2); Julia passes the second Tsit5 stage to the Hermite
#   interpolant.
# * JVODE: an Adams order decrease measures the past steps in units of the
#   step the Nordsieck vector is scaled to (.jvode_decrease_scale), as CVODE
#   does with hscale; Julia divides by the new step.
# * JVODE: a restart after a callback also resets L and n_wait as in a fresh
#   cache (.jvode_restart); Julia keeps them, so the first order change
#   after the restart can use the saved error vector of a step before it.

# The product a * b and its rounding error (Dekker and Veltkamp), as a
# double double list(hi, lo).
.dd_two_prod <- function(a, b) {
  p <- a * b
  ta <- a * 134217729; ah <- ta - (ta - a); al <- a - ah
  tb <- b * 134217729; bh <- tb - (tb - b); bl <- b - bh
  list(p, ((ah * bh - p) + ah * bl + al * bh) + al * bl)
}

# Double double arithmetic (about 106 bits), which stands in for the 256 bit
# BigFloat arithmetic of Julia below.
.dd_norm <- function(h, l) {
  s <- h + l
  list(s, l - (s - h))
}

.dd_add <- function(a, b) {
  s <- a[[1L]] + b[[1L]]
  bb <- s - a[[1L]]
  .dd_norm(s, ((a[[1L]] - (s - bb)) + (b[[1L]] - bb)) + (a[[2L]] + b[[2L]]))
}

.dd_mul <- function(a, d) {
  tp <- .dd_two_prod(a[[1L]], d)
  .dd_norm(tp[[1L]], tp[[2L]] + a[[2L]] * d)
}

.dd_div <- function(a, d) {
  q1 <- a[[1L]] / d
  tp <- .dd_two_prod(q1, d)
  .dd_norm(q1, (((a[[1L]] - tp[[1L]]) - tp[[2L]]) + a[[2L]]) / d)
}

# Julia's pow_body(x::Float64, n::Integer) (base/special/pow.jl), the
# compensated power by squaring of x^n for an integer n >= 2. The last
# muladd is not fused: the compiled code shares x * y with the branch for a
# zero error, so it multiplies and adds separately.
.ns_pow_int <- function(x, n) {
  if (n == 3L) return(x * x * x)
  y <- 1; xnlo <- -0; ynlo <- 0
  while (n > 1L) {
    if (n %% 2L == 1L) {
      err <- .fma(y, xnlo, x * ynlo)
      tm <- .dd_two_prod(x, y)
      y <- tm[[1L]]; ynlo <- tm[[2L]] + err
    }
    err <- x * 2 * xnlo
    tm <- .dd_two_prod(x, x)
    x <- tm[[1L]]; xnlo <- tm[[2L]] + err
    n <- n %/% 2L
  }
  err <- .fma(y, xnlo, x * ynlo)
  if (is.finite(x) && is.finite(err) && err != 0) x * y + err else x * y
}

# The default internal norm of OrdinaryDiffEq (a fused sum of squares) and
# calculate_residuals (err / (abstol + max(|u0|, |u1|) * reltol), fused).
.ns_norm <- function(x) {
  s <- 0
  for (v in x) s <- .fma(v, v, s)
  sqrt(s / max(length(x), 1L))
}

.ns_err <- function(err, u0, u1, ctx) {
  .ns_norm(err / .fma(pmax(abs(u0), abs(u1)), ctx$reltol, ctx$abstol))
}

# perform_predict!: the Pascal triangle operator (rewind = TRUE undoes it).
.ns_predict <- function(z, order, rewind = FALSE) {
  for (i in seq_len(order)) for (j in order:i)
    z[[j]] <- if (rewind) z[[j]] - z[[j + 1L]] else z[[j]] + z[[j + 1L]]
  z
}

# nordsieck_rescale!: the vector scaled from the step dts[2] to dts[1].
.ns_rescale <- function(z, dts, order, rewind = FALSE) {
  eta <- if (rewind) dts[2L] / dts[1L] else dts[1L] / dts[2L]
  factor <- eta
  for (i in 2:(order + 1L)) {
    z[[i]] <- z[[i]] * factor
    factor <- factor * eta
  }
  z
}

.ns_rewind <- function(cache) {
  cache$z <- .ns_rescale(.ns_predict(cache$z, cache$order, TRUE), cache$dts, cache$order, TRUE)
  invisible(cache)
}

.ns_restore_dts <- function(cache, oldest) {
  cache$dts <- c(cache$dts[-1L], oldest)
  invisible(cache)
}

# update_nordsieck_vector!: the corrector z[i] + l[i] * Delta.
.ns_correct <- function(cache) {
  for (i in seq_len(cache$order + 1L)) cache$z[[i]] <- .fma(cache$l[i], cache$Delta, cache$z[[i]])
  invisible(cache)
}

# The integral from -1 to 0 of the polynomial with coefficients a (degree
# deg), each power raised by k.
.ns_int <- function(a, deg, k) {
  s <- 0; sg <- 1
  for (i in 0:deg) {
    s <- s + sg * (a[i + 1L] / (i + k))
    sg <- -sg
  }
  s
}

# calc_coeff!: the Adams coefficients l and the error constants for the step
# sizes in dts (Byrne and Hindmarsh). vode adds the constants of the order
# selection.
.ns_coeff <- function(cache, vode) {
  m <- cache$m; l <- cache$l; dts <- cache$dts; order <- cache$order
  varorder <- vode && cache$n_wait == 1L
  dt <- dtsum <- dts[1L]
  if (order == 1L) {
    l[1L] <- l[2L] <- 1
    cache$l <- l
    cache$c_LTEm1 <- cache$c_D <- 1
    cache$c_LTE <- 0.5
    cache$c_LTEp1 <- 1 / 12
    cache$c_conv <- 0.1 / cache$c_LTE
    return(invisible(cache))
  }
  m[1L] <- 1
  m[2:(order + 1L)] <- 0
  for (j in seq_len(order - 1L)) {
    if (varorder && j == order - 1L) {
      Mm1 <- .ns_int(m, order - 2L, 2L)
      cache$c_LTEm1 <- order * Mm1 / m[order - 1L]
    }
    xi_inv <- dt / dtsum
    for (i in j:1) m[i + 1L] <- .fma(m[i], xi_inv, m[i + 1L])
    dtsum <- dtsum + dts[j + 1L]
  }
  xi_inv <- dt / dtsum
  M0 <- .ns_int(m, order - 1L, 1L)
  M1 <- .ns_int(m, order - 1L, 2L)
  M0_inv <- 1 / M0
  l[1L] <- 1
  for (i in seq_len(order)) l[i + 1L] <- M0_inv * m[i] / i
  cache$c_LTE <- M1 * M0_inv * xi_inv
  if (vode) cache$c_D <- (1 / xi_inv) / l[order + 1L]
  if (varorder) {
    for (i in order:1) m[i + 1L] <- .fma(xi_inv, m[i], m[i + 1L])
    M2 <- .ns_int(m, order, 2L)
    cache$c_LTEp1 <- M2 * M0_inv / (order + 1L)
  }
  cache$c_conv <- 0.1 / cache$c_LTE
  cache$m <- m; cache$l <- l
  invisible(cache)
}

# calc_coeff_bdf!: the fixed leading coefficient BDF coefficients.
.ns_coeff_bdf <- function(cache) {
  l <- cache$l; dts <- cache$dts; q <- cache$order
  dt <- dts[1L]
  l[seq_len(q + 1L)] <- 0
  l[1L] <- l[2L] <- 1
  xi_inv <- xistar_inv <- 1
  a0 <- a0_hat <- -1
  dtsum <- dt
  if (q > 1L) {
    if (q > 2L) for (j in 2:(q - 1L)) {
      dtsum <- dtsum + dts[j]
      xi_inv <- dt / dtsum
      a0 <- a0 - 1 / j
      for (i in j:1) l[i + 1L] <- .fma(l[i], xi_inv, l[i + 1L])
    }
    a0 <- a0 - 1 / q
    xistar_inv <- -l[2L] - a0
    dtsum <- dtsum + dts[q]
    xi_inv <- dt / dtsum
    a0_hat <- -l[2L] - xi_inv
    for (i in q:1) l[i + 1L] <- .fma(l[i], xistar_inv, l[i + 1L])
  }
  A1 <- 1 - a0_hat + a0
  A2 <- 1 + q * A1
  cache$c_LTE <- abs(A1 / (a0 * A2))
  cache$c_D <- abs(A2 * xistar_inv / (l[q + 1L] * xi_inv))
  if (cache$n_wait == 1L) {
    if (q > 1L) {
      A3 <- a0 + 1 / q
      A4 <- a0_hat + xi_inv
      cache$c_LTEm1 <- abs(xistar_inv / l[q + 1L] * (1 - A4 + A3) / A3)
    } else cache$c_LTEm1 <- 1
    dtsum <- dtsum + dts[q + 1L]
    xi_inv <- dt / dtsum
    A5 <- a0 - 1 / (q + 1L)
    A6 <- a0_hat - xi_inv
    cache$c_LTEp1 <- abs((1 - A6 + A5) / A2 / (xi_inv * (q + 2L) * A5))
  }
  cache$c_conv <- 0.1 / cache$c_LTE
  cache$l <- l
  invisible(cache)
}

# nlsolve_functional!: at most three functional iterations of the corrector.
# The state, Delta (u - z[1]) and the number of evaluations of f; ok is
# FALSE when the iteration does not converge.
.ns_functional <- function(ctx, t, dt, cache) {
  z <- cache$z
  rate <- ctx$f(z[[1L]], ctx$p, dt + t)
  nf <- 1L
  inv_l2 <- 1 / cache$l[2L]
  Delta <- 0 * z[[1L]]
  conv_rate <- 1
  d_prev <- 0
  for (k in 1:3) {
    rate <- inv_l2 * .fma(dt, rate, -z[[2L]])
    u <- rate + z[[1L]]
    d <- .ns_norm(rate - Delta)
    Delta <- rate
    conv_rate <- max(0.1 * conv_rate, d / d_prev)
    if (d * min(1, conv_rate) / cache$c_conv <= 1) {
      cache$Delta <- Delta
      return(list(ok = TRUE, u = u, nf = nf))
    }
    if (k == 3L || (k >= 2L && d > 2 * d_prev)) return(list(ok = FALSE, nf = nf))
    d_prev <- d
    nf <- nf + 1L
    rate <- ctx$f(u, ctx$p, dt + t)
  }
}

# After a callback changed u the derivative is evaluated again and the next
# step restarts.
.ns_on_event <- function(ctx, t, u, cache) {
  cache$fsal <- ctx$f(u, ctx$p, t)
  list(cache = cache, nf = 1L)
}

# ---------------------------------------------------------------------------
# AN5: the fixed order 5 Adams method in Nordsieck form, started with Tsit5

.an5_init <- function(ctx, t, u) {
  cache <- new.env(parent = emptyenv())
  cache$fsal <- ctx$f(u, ctx$p, t)
  cache$z <- rep(list(0 * u), 6L)
  cache$l <- cache$m <- numeric(6L)
  cache$c_LTE <- cache$c_conv <- 0
  cache$dts <- numeric(6L)
  cache$Delta <- u
  cache$order <- 1L
  cache
}

# The Tsit5 tableau of OrdinaryDiffEqTsit5 (Tsit5ConstantCacheActual) as
# hexadecimal literals: the Nordsieck vector is built from the Tsit5 stages,
# and decimal literals can be parsed one unit in the last place away from
# Julia's values.
.AN5_TSIT5 <- list(
  c = c(0x1.49ba5e353f7cfp-3, 0x1.4ed916872b021p-2, 0x1.ccccccccccccdp-1, 0x1.f5c5e85edbf1dp-1),
  a = list(0x1.49ba5e353f7cfp-3,
           c(-0x1.15e4e4fe8123ep-7, 0x1.57883daf1f0b3p-2),
           c(0x1.72d5e94fe3235p+1, -0x1.9701343d4287bp+2, 0x1.172fd92eea8fap+2),
           c(0x1.54daf827440c7p+2, -0x1.77f6daaa052d3p+3, 0x1.dfb6eaa5e5f6ap+2, -0x1.7adc1b510a9a4p-4),
           c(0x1.772216029aa0ep+2, -0x1.9d7894a598ba8p+3, 0x1.05198ab50d56dp+3, -0x1.253648f20ef47p-4,
             -0x1.cf28fdbb9db09p-6),
           c(0x1.8b1a71ed4a6eep-4, 0x1.47ae147ae147bp-7, 0x1.eb683148e0688p-2, 0x1.6106b4b6adddbp+0,
             -0x1.a520ff752ea91p+1, 0x1.29901d4cd6786p+1)),
  btilde = c(-0x1.d29e8505f7d7fp-10, -0x1.ac0bfa6dad476p-11, 0x1.023d98a8af021p-7, -0x1.285e3e99716ep-3,
             0x1.2a2ab7e86b834p-1, -0x1.d51379bf6508ap-2, 0x1.f07c1f07c1f08p-7))

# The weights of the stages in the derivatives of order 2, 3 and 4 of the
# Tsit5 interpolant at theta = 0 (2 r_i2, 6 r_i3 and 24 r_i4). AN5 calls the
# interpolant without the end values, so Julia takes Tsit5Interp(Any): most
# coefficients are 256 bit BigFloat values (given here as double doubles
# hi + lo), the others (big = FALSE) Float64 values whose products with the
# stages are rounded to Float64.
.AN5_TAYLOR <- list(
  list(hi = c(-0x1.61c11fea8621ap+2, 0x1.0db8bac710cb2p-2, 0x1.f713f27487049p+2, -0x1.8d278b4d83996p+4,
              0x1.2c1312e99aa1bp+6, -0x1.be582bf341b49p+5, 0x1.8p+1),
       lo = c(-0x1.ed9eeeca1c7f9p-53, 0x1.f2ef67c8869e8p-56, -0x1.5f9d7ca6bfa24p-53, -0x1.6bff2d7dac927p-50,
              -0x1.3502cbec53efp-48, -0x1.b95d59d247653p-49, 0),
       big = c(rep(TRUE, 6L), FALSE)),
  list(hi = c(0x1.17ac2a8e082cep+4, -0x1.572474538ef35p+0, -0x1.1d2b6c99bb2fdp+5, 0x1.6c0eea5c634cap+7,
              -0x1.0889662230592p+9, 0x1.868d2674d97ep+8, -0x1.8p+4),
       lo = c(-0x1.4bfe7f57e0592p-50, 0x1.afad0fd1e4b23p-55, -0x1.8e3307db4b282p-50, -0x1.410e8e9fdb391p-48,
              0x1.aa65214749541p-45, 0x1.ffdd838a719b2p-47, 0),
       big = c(rep(TRUE, 6L), FALSE)),
  list(hi = c(-0x1.9462cfca4f542p+4, 0x1.386c226809d49p+1, 0x1.de334ff5cc4c1p+5, -0x1.8d278b4d83996p+8,
              0x1.1c46f21bae89p+10, -0x1.a272a9340d995p+9, 0x1.ep+5),
       lo = c(0, 0, -0x1.14f1ee49c1263p-53, -0x1.8c6e8474f873cp-46, 0x1.735bc58d03825p-44,
              0x1.4bfab7a9a25c3p-46, 0),
       big = c(FALSE, FALSE, rep(TRUE, 4L), FALSE)))

# sum_i a[i] * ks[[i]], left to right as in Julia's @muladd.
.an5_lc <- function(a, ks) {
  s <- a[1L] * ks[[1L]]
  for (i in seq_along(a)[-1L]) s <- s + a[i] * ks[[i]]
  s
}

# The Tsit5 step of tsit_perform_step.jl (stage times fused) from the
# derivative k1 at the start.
.an5_tsit5 <- function(ctx, t, u, dt, k1) {
  f <- ctx$f; p <- ctx$p; A <- .AN5_TSIT5$a; cc <- .AN5_TSIT5$c
  k <- list(k1)
  k[[2L]] <- f(u + (dt * A[[1L]]) * k1, p, .fma(cc[1L], dt, t))
  for (s in 2:4) k[[s + 1L]] <- f(u + dt * .an5_lc(A[[s]], k), p, .fma(cc[s], dt, t))
  k[[6L]] <- f(u + dt * .an5_lc(A[[5L]], k), p, t + dt)
  unew <- u + dt * .an5_lc(A[[6L]], k)
  k[[7L]] <- f(unew, p, t + dt)
  list(u = unew, err = dt * .an5_lc(.AN5_TSIT5$btilde, k), k = k)
}

# The Taylor coefficients z[3], z[4], z[5] of the Tsit5 dense output at the
# start of its step: sum_i k_i w_i, times dt^(1 - d) and dt^d / d!, in the
# mixed BigFloat and Float64 arithmetic of Julia, rounded to Float64.
.an5_tsit5_taylor <- function(k, dt) {
  idt <- 1 / dt
  taylor <- function(w, s1, s2, fac) {
    s <- NULL
    for (i in 1:7) {
      if (w$big[i]) {
        tp <- .dd_two_prod(k[[i]], w$hi[i])
        tm <- .dd_norm(tp[[1L]], tp[[2L]] + k[[i]] * w$lo[i])
      } else tm <- k[[i]] * w$hi[i]
      s <- if (is.null(s)) tm
        else if (!is.list(s) && !is.list(tm)) tm + s
        else .dd_add(if (is.list(s)) s else list(s, 0 * s), if (is.list(tm)) tm else list(tm, 0 * tm))
    }
    .dd_div(.dd_mul(.dd_mul(s, s1), s2), fac)[[1L]]
  }
  list(taylor(.AN5_TAYLOR[[1L]], idt, dt * dt, 2),
       taylor(.AN5_TAYLOR[[2L]], idt * idt, dt * dt * dt, 6),
       taylor(.AN5_TAYLOR[[3L]], idt * idt * idt, .ns_pow_int(dt, 4L), 24))
}

.an5_start_coeffs <- function(cache) {
  cache$l[] <- 0
  invisible(cache)
}

.an5_start_k2 <- function(k) k[[7L]]

.an5_reject <- function(cache, oldest) {
  .ns_rewind(cache)
  .ns_restore_dts(cache, oldest)
  TRUE
}

.an5_fail <- function(cache, oldest) {
  .ns_rewind(cache)
  .ns_restore_dts(cache, oldest)
}

.an5_step <- function(ctx, t, u, dt, cache) {
  adaptive <- isTRUE(ctx$adaptive)
  if (isTRUE(ctx$derivative_discontinuity)) cache$order <- 1L
  if (cache$order == 1L) {
    # build the Nordsieck vector from one Tsit5 step
    r <- .an5_tsit5(ctx, t, u, dt, cache$fsal)
    k <- r$k
    cache$order <- 4L
    tay <- .an5_tsit5_taylor(k, dt)
    cache$z <- c(list(u, k[[1L]] * dt), tay, list(0 * u))
    cache$dts[] <- dt
    cache$z <- .ns_predict(cache$z, cache$order)
    cache$Delta <- r$u - u
    .an5_start_coeffs(cache)
    .ns_correct(cache)
    EEst <- if (adaptive) .ns_err(r$err, u, r$u, ctx) else NULL
    if (adaptive && EEst >= 1) cache$order <- 1L
    return(list(u = r$u, EEst = EEst, k = list(k[[1L]], .an5_start_k2(k)), cache = cache, nf = 6L))
  }
  oldest <- cache$dts[6L]
  cache$dts <- c(dt, cache$dts[1:5])
  if (dt != cache$dts[2L]) cache$z <- .ns_rescale(cache$z, cache$dts, cache$order)
  k1 <- cache$z[[2L]] / dt
  cache$z <- .ns_predict(cache$z, cache$order)
  .ns_coeff(cache, FALSE)
  nl <- .ns_functional(ctx, t, dt, cache)
  if (!nl$ok) {
    .an5_fail(cache, oldest)
    return(list(fail = TRUE, nreject = 0L, nf = nl$nf, cache = cache))
  }
  EEst <- NULL
  if (adaptive) {
    EEst <- .ns_err(cache$Delta, u, nl$u, ctx) * cache$c_LTE
    if (EEst > 1 && .an5_reject(cache, oldest)) return(list(u = nl$u, EEst = EEst, cache = cache, nf = nl$nf))
  }
  cache$order <- 5L
  .ns_correct(cache)
  list(u = nl$u, EEst = EEst, k = list(k1, cache$z[[2L]] / dt), cache = cache, nf = nl$nf)
}

# ---------------------------------------------------------------------------
# JVODE: variable order Adams (orders 1 to 12) or BDF (1 to 5) methods in
# Nordsieck form with CVODE's order and step size selection

.jvode_init <- function(ctx, t, u) {
  cache <- new.env(parent = emptyenv())
  cache$fsal <- ctx$f(u, ctx$p, t)
  cache$z <- rep(list(0 * u), 13L)
  cache$l <- cache$m <- numeric(13L)
  cache$c_LTEp1 <- cache$c_LTE <- cache$c_LTEm1 <- cache$c_conv <- cache$c_D <- cache$prev_D <- 0
  cache$dts <- numeric(13L)
  cache$Delta <- u
  cache$L <- 2L
  cache$order <- cache$nextorder <- 1L
  cache$n_wait <- 2L
  cache$eta <- cache$etaold <- cache$etaq <- cache$eta_p1 <- cache$eta_m1 <- 0
  cache
}

.jvode_restart <- function(cache) {
  cache$L <- 2L
  cache$n_wait <- 2L
  invisible(cache)
}

.jvode_decrease_scale <- function(dts) dts[2L]

# nordsieck_adjust_order!: the Adams vector at order + dorder.
.jvode_adjust_adams <- function(cache, dorder) {
  order <- cache$order; dts <- cache$dts
  if (dorder == 1L) {
    cache$z[[order + 2L]] <- 0 * cache$z[[order + 2L]]
    return(invisible(cache))
  }
  l <- numeric(length(cache$l))
  l[2L] <- 1
  h <- .jvode_decrease_scale(dts)
  hsum <- 0
  for (j in seq_len(order - 2L)) {
    hsum <- hsum + dts[j + 1L]
    xi <- hsum / h
    for (i in (j + 1L):1) l[i + 1L] <- l[i + 1L] * xi + l[i]
  }
  if (order >= 3L) {
    for (j in 2:(order - 1L)) l[j + 1L] <- order * l[j] / j
    for (j in 3:order) cache$z[[j]] <- .fma(-l[j], cache$z[[order + 1L]], cache$z[[j]])
  }
  cache$l <- l
  invisible(cache)
}

# nordsieck_adjust_order_bdf! (cvIncreaseBDF and cvDecreaseBDF).
.jvode_adjust_bdf <- function(cache, dorder) {
  z <- cache$z; dts <- cache$dts; q <- cache$order
  scale <- dts[2L]
  l <- numeric(length(cache$l))
  l[3L] <- 1
  if (dorder == 1L) {
    a0 <- -1
    a1 <- prd <- xiold <- 1
    dtsum <- scale
    for (j in seq_len(q - 1L)) {
      dtsum <- dtsum + dts[j + 2L]
      xi <- dtsum / scale
      prd <- prd * xi
      a0 <- a0 - 1 / (j + 1L)
      a1 <- a1 + 1 / xi
      for (i in (j + 2L):2) l[i + 1L] <- .fma(l[i + 1L], xiold, l[i])
      xiold <- xi
    }
    A1 <- (-a0 - a1) / prd
    z[[q + 2L]] <- A1 * z[[13L]]
    for (j in seq_len(q)[-1L]) z[[j + 1L]] <- .fma(l[j + 1L], z[[q + 2L]], z[[j + 1L]])
  } else if (q > 2L) {
    dtsum <- 0
    for (j in 1:(q - 2L)) {
      dtsum <- dtsum + dts[j + 1L]
      xi <- dtsum / scale
      for (i in (j + 2L):2) l[i + 1L] <- .fma(l[i + 1L], xi, l[i])
    }
    for (j in 2:(q - 1L)) z[[j + 1L]] <- .fma(-l[j + 1L], z[[q + 1L]], z[[j + 1L]])
  }
  cache$l <- l; cache$z <- z
  invisible(cache)
}

.jvode_eta <- function(bias, err, expo, addon) 1 / ((bias * err)^(1 / expo) + addon)

# The step of JVODE(); par holds the constructor arguments.
.jvode_step_fn <- function(par) {
  bdf <- par$algorithm == "BDF"
  qmax_order <- if (bdf) 5L else 12L
  set_eta <- function(ctx, cache) {
    if (cache$eta < par$qsteady_max) cache$eta <- 1
    else if (ctx$iter == 1L || isTRUE(ctx$derivative_discontinuity)) cache$eta <- min(1e5, cache$eta)
    else cache$eta <- min(par$qmax, max(par$qmin, cache$eta))
  }
  choose_eta <- function(cache) {
    eta <- max(cache$etaq, cache$eta_m1, cache$eta_p1)
    order <- cache$order
    if (eta < par$qsteady_max) {
      cache$eta <- 1
      cache$nextorder <- order
      if (bdf) return(invisible())
    }
    if (eta == cache$etaq) {
      cache$eta <- cache$etaq
      cache$nextorder <- order
    } else if (eta == cache$eta_m1) {
      cache$eta <- cache$eta_m1
      cache$nextorder <- order - 1L
    } else {
      cache$eta <- cache$eta_p1
      cache$nextorder <- order + 1L
      if (bdf) cache$z[[13L]] <- cache$Delta
    }
  }
  # nordsieck_prepare_next!: the step size ratio eta and the next order
  prepare_next <- function(ctx, cache, EEst, u, unew) {
    order <- cache$order
    if (EEst > 1) {
      .ns_rewind(cache)
      cache$n_wait <- max(2L, cache$n_wait)
      cache$nextorder <- order
      cache$eta <- .jvode_eta(par$bias2, EEst, cache$L, par$addon)
      return(invisible())
    }
    cache$etaq <- .jvode_eta(par$bias2, EEst, order + 1L, par$addon)
    if (cache$n_wait != 0L) {
      cache$eta <- cache$etaq
      cache$nextorder <- order
      set_eta(ctx, cache)
      return(invisible())
    }
    cache$n_wait <- 2L
    # stepsize_eta+1!: order + 1 from the change of the saved error vector
    cache$eta_p1 <- 0
    if (order != qmax_order && cache$prev_D != 0) {
      L <- order + 1L
      cquot <- (cache$c_D / cache$prev_D) * .ns_pow_int(cache$dts[1L] / cache$dts[2L], L)
      dup <- .ns_err(.fma(-cquot, cache$z[[13L]], cache$Delta), u, unew, ctx) * cache$c_LTEp1
      cache$eta_p1 <- .jvode_eta(par$bias3, dup, L + 1L, par$addon)
    }
    # stepsize_eta-1!: order - 1 from the last column of the vector
    cache$eta_m1 <- 0
    if (order > 1L) {
      approx <- .ns_err(cache$z[[order + 1L]], u, unew, ctx) * cache$c_LTEm1
      cache$eta_m1 <- .jvode_eta(par$bias1, approx, order, par$addon)
    }
    choose_eta(cache)
    set_eta(ctx, cache)
    cache$Delta <- cache$c_LTE * cache$Delta
    invisible()
  }
  function(ctx, t, u, dt, cache) {
    nf <- 0L
    if (isTRUE(ctx$derivative_discontinuity) || ctx$iter == 1L) {
      cache$order <- cache$nextorder <- 1L
      if (isTRUE(ctx$derivative_discontinuity)) .jvode_restart(cache)
      cache$z[[1L]] <- u
      cache$z[[2L]] <- ctx$f(u, ctx$p, t) * dt
      nf <- 1L
      cache$dts[1L] <- dt
    }
    oldest <- cache$dts[13L]
    cache$dts <- c(dt, cache$dts[1:12])
    if (dt != cache$dts[2L] || (bdf && cache$nextorder != cache$order)) {
      # nordsieck_adjust!
      if (cache$nextorder != cache$order) {
        dorder <- cache$nextorder - cache$order
        if (bdf) .jvode_adjust_bdf(cache, dorder) else .jvode_adjust_adams(cache, dorder)
        cache$order <- cache$nextorder
        cache$L <- cache$order + 1L
        cache$n_wait <- cache$L
      }
      cache$z <- .ns_rescale(cache$z, cache$dts, cache$order)
    }
    k1 <- cache$z[[2L]] / dt
    cache$z <- .ns_predict(cache$z, cache$order)
    if (bdf) .ns_coeff_bdf(cache) else .ns_coeff(cache, TRUE)
    nl <- .ns_functional(ctx, t, dt, cache)
    nf <- nf + nl$nf
    if (!nl$ok) {
      .ns_rewind(cache)
      .ns_restore_dts(cache, oldest)
      return(list(fail = TRUE, nreject = 0L, nf = nf, cache = cache))
    }
    # without adaptivity the error estimate keeps its initial value 1
    EEst <- if (isTRUE(ctx$adaptive)) .ns_err(cache$Delta, u, nl$u, ctx) * cache$c_LTE else 1
    # nordsieck_complete_step!
    if (EEst > 1) {
      prepare_next(ctx, cache, EEst, u, nl$u)
      .ns_restore_dts(cache, oldest)
      return(list(u = nl$u, EEst = EEst, cache = cache, nf = nf))
    }
    .ns_correct(cache)
    cache$n_wait <- cache$n_wait - 1L
    if (cache$n_wait == 1L && cache$order != qmax_order) {
      cache$z[[13L]] <- cache$Delta
      cache$prev_D <- cache$c_D
    }
    prepare_next(ctx, cache, EEst, u, nl$u)
    list(u = nl$u, EEst = EEst, k = list(k1, cache$z[[2L]] / dt), cache = cache, nf = nf)
  }
}

#' Nordsieck form Adams and BDF methods
#'
#' Multistep methods of OrdinaryDiffEq.jl (OrdinaryDiffEqNordsieck) that
#' store the solution history as a Nordsieck vector (the scaled derivatives
#' of the interpolating polynomial) and solve the corrector by functional
#' iteration, so they suit nonstiff problems. The Julia documentation marks
#' them as experimental and recommends `VCABM()` (Adams) and `FBDF()` (BDF).
#'
#' * `AN5()`: the adaptive fixed order 5 Adams method. The first step (and
#'   the first step after a callback that changes `u`) is a Tsit5 step that
#'   builds the Nordsieck vector; the PI controller sets the step size.
#' * `JVODE()`: the variable step, variable order method of CVODE
#'   (Brown, Byrne and Hindmarsh), Adams of orders 1 to 12 or BDF of orders
#'   1 to 5, which chooses the order and step size from the error estimates
#'   at orders q - 1, q and q + 1. `JVODE_Adams()` and `JVODE_BDF()` are
#'   `JVODE("Adams")` and `JVODE("BDF")`. An Adams order change is applied
#'   only together with a step size change, as in the Julia version; with a
#'   fixed `dt` the Adams method therefore stays at order 1.
#'
#' The dense output is the cubic Hermite interpolant of each step.
#'
#' @param algorithm `"Adams"` or `"BDF"`.
#' @param bias1,bias2,bias3 Safety factors of CVODE for the error estimates
#'   at orders q - 1, q and q + 1.
#' @param addon Constant added to the step size ratio denominators.
#' @param qmin,qmax Smallest and largest step size ratio.
#' @param qsteady_min,qsteady_max The step size is kept when the ratio of the
#'   old to the new step lies in this interval (a ratio of new to old below
#'   `qsteady_max` keeps it too).
#' @param ... Arguments of `JVODE()`.
#' @return An ODE algorithm.
#' @examples
#' prob <- ODEProblem(function(u, p, t) c(u[2], -u[1]), c(0, 1), c(0, 1))
#' sol <- solve(prob, JVODE_Adams(), reltol = 1e-8, abstol = 1e-8)
#' sol(1)
#' solve(prob, AN5(), reltol = 1e-6, abstol = 1e-8)$stats$nf
#' @name nordsieck_methods
NULL

#' @rdname nordsieck_methods
#' @export
AN5 <- function() {
  .new_alg("AN5", 5L, step = .an5_step, init = .an5_init, adaptive_order = 5L, qsteady_max = 1.5,
           on_event = .ns_on_event)
}

#' @rdname nordsieck_methods
#' @export
JVODE <- function(algorithm = "Adams", bias1 = 6, bias2 = 6, bias3 = 10, addon = 1e-6, qmin = 0.2,
                  qmax = 10, qsteady_min = 1, qsteady_max = 1.5) {
  algorithm <- match.arg(algorithm, c("Adams", "BDF"))
  par <- list(algorithm = algorithm, bias1 = bias1, bias2 = bias2, bias3 = bias3, addon = addon, qmin = qmin,
              qmax = qmax, qsteady_min = qsteady_min, qsteady_max = qsteady_max)
  # JVODEController: eta from the step, dt / q with q = 1 / eta kept at 1
  # inside [qsteady_min, qsteady_max]
  .new_alg("JVODE", 1L, step = .jvode_step_fn(par), init = .jvode_init, controller = "BDF",
           algorithm = algorithm, qmin = qmin, qmax = qmax, qsteady_min = qsteady_min,
           qsteady_max = qsteady_max, on_event = .ns_on_event,
           bdf_stepsize = function(cache, EEst, co, first) {
             if (EEst == 0) return(par$qmax)
             cache$etaold <- cache$eta
             cache$eta
           },
           bdf_accept = function(cache, eta, dt, co) {
             q <- 1 / eta
             if (q <= par$qsteady_max && q >= par$qsteady_min) q <- 1
             dt / q
           },
           bdf_reject = function(cache, EEst, dt, co) dt * cache$eta,
           newton_fail = function(cache, dtcur, co) dtcur / co$failfactor)
}

#' @rdname nordsieck_methods
#' @export
JVODE_Adams <- function(...) JVODE("Adams", ...)

#' @rdname nordsieck_methods
#' @export
JVODE_BDF <- function(...) JVODE("BDF", ...)
