# Stabilized explicit Runge-Kutta methods of OrdinaryDiffEqStabilizedRK
# (rkc_perform_step.jl, rkc_utils.jl, rkc_caches.jl, alg_utils.jl): the
# Chebyshev methods ROCK2, ROCK4, RKC, TSRKC2, TSRKC3 and RKMC2, the
# extrapolated methods ESERK4, ESERK5 and SERK2, and the Legendre and
# Gegenbauer super time stepping methods RKL1, RKL2, RKG1 and RKG2. The
# coefficient tables of ROCK2, ROCK4, ESERK4, ESERK5 and SERK2 are
# .STAB_TABLES in R/sysdata.rda (tools/stabilized-dump.jl and
# tools/stabilized-gen.R).

#' Stabilized explicit Runge-Kutta methods
#'
#' Explicit methods whose number of stages grows with the stiffness, for
#' mildly stiff problems whose Jacobian has eigenvalues near the negative
#' real axis (diffusion, parabolic PDEs after a space discretization). Each
#' step estimates the spectral radius of the Jacobian by a nonlinear power
#' iteration (or takes it from `eigen_est`) and picks the number of stages
#' that keeps the step stable.
#'
#' * `ROCK2()`, `ROCK4()`: the orthogonal Runge-Kutta Chebyshev methods of
#'   Abdulle (order 2) and Abdulle (order 4).
#' * `RKC()`: the Runge-Kutta Chebyshev method of Sommeijer, Shampine and
#'   Verwer (order 2).
#' * `TSRKC2()`, `TSRKC3()`: two step Runge-Kutta Chebyshev methods of
#'   orders 2 and 3; `TSRKC3()` starts with RKC steps.
#' * `RKMC2()`: the monotone second order Chebyshev method.
#' * `ESERK4()`, `ESERK5()`, `SERK2()`: the extrapolated stabilized methods
#'   of Martin-Vaquero (orders 4, 5 and 2).
#' * `RKL1()`, `RKL2()`: Runge-Kutta Legendre super time stepping (orders 1
#'   and 2); the stage counts are odd.
#' * `RKG1()`, `RKG2()`: Runge-Kutta Gegenbauer super time stepping (orders
#'   1 and 2).
#'
#' ROCK2, ROCK4, ESERK4, ESERK5 use the PI step size controller; the others
#' use the predictive (Gustafsson) controller, as in OrdinaryDiffEq. The
#' methods with a stability bound cap the proposed step at the largest
#' stable step for their maximal number of stages.
#'
#' @param min_stages,max_stages Bounds on the number of stages.
#' @param eigen_est An upper bound on the spectral radius of the Jacobian:
#'   a number, or a function `(u, p, t)` returning one. `NULL` (the default)
#'   estimates it by power iteration.
#' @param eigen_est_interval With the power iteration, the estimate is
#'   recomputed every `eigen_est_interval` steps, and always on the first
#'   step, after a rejected step and after a discontinuity.
#' @return An ODE algorithm.
#' @examples
#' n <- 20
#' heat <- function(u, p, t) {
#'   (c(0, u[-n]) - 2 * u + c(u[-1], 0)) * (n + 1)^2
#' }
#' prob <- ODEProblem(heat, sin(pi * seq_len(n) / (n + 1)), c(0, 0.1))
#' sol <- solve(prob, ROCK2(), reltol = 1e-4)
#' sol(0.1)[1:3]
#' @name stabilized_rk_methods
NULL

# ---------------------------------------------------------------------------
# Spectral radius estimation (maxeig!, maybe_maxeig!)

# eps(x) of Julia: the spacing of doubles at |x|.
.eps64 <- function(x) {
  x <- abs(x)
  if (is.na(x) || is.infinite(x)) return(NaN)
  hi <- readBin(writeBin(x, raw(), size = 8L, endian = "little")[5:8], "integer", size = 4L, endian = "little")
  e <- bitwAnd(bitwShiftR(hi, 20L), 0x7FFL)
  if (e == 0L) 2^-1074 else 2^(e - 1075)
}

.rms_norm <- function(x) sqrt(sum(x^2) / length(x))

# The nonlinear power iteration of maxeig!: sets st$eig and returns the
# number of function evaluations. `rkc` selects the RKCAlgs variant (start
# from the derivative, convergence to 1 percent, safety factor 1.2 after
# convergence); RKMC2 uses the generic one.
.stab_maxeig <- function(ctx, st, t, u, fsal, rkc, maxiter) {
  f <- ctx$f; p <- ctx$p
  nf <- 0L
  if (ctx$iter == 1L || st$disc) {
    if (rkc) z <- fsal else {
      z <- f(fsal, p, t)
      nf <- 1L
    }
  } else {
    z <- st$zprev
  }
  safe <- if (rkc) 1 else 1.2
  u_norm <- .rms_norm(u)
  z_norm <- .rms_norm(z)
  pert <- .eps64(u_norm)
  sqrt_pert <- sqrt(pert)
  if (u_norm != 0 && z_norm != 0) {
    dz_u <- u_norm * sqrt_pert
    z <- u + (dz_u / z_norm) * z
  } else if (u_norm != 0) {
    dz_u <- u_norm * sqrt_pert
    z <- u + u * dz_u
  } else if (z_norm != 0) {
    dz_u <- pert
    z <- z * (dz_u / z_norm)
  } else {
    dz_u <- pert
    z <- rep(dz_u, length(z))
  }
  eig <- 0
  for (it in seq_len(maxiter)) {
    fz <- f(z, p, t)
    nf <- nf + 1L
    tmp <- fz - fsal
    Delta <- .rms_norm(tmp)
    eig_prev <- eig
    eig <- Delta / dz_u * safe
    if (rkc) {
      if (it >= 2L && abs(eig_prev - eig) < max(eig, 1 / abs(ctx$dtmax)) * 0.01) {
        st$eig <- eig * 1.2
        st$zprev <- z - u
        return(nf)
      }
    } else if (it >= 2L && abs(eig_prev - eig) < eig * 0.05) {
      st$eig <- eig
      st$zprev <- z - u
      return(nf)
    }
    if (Delta != 0) {
      z <- u + (dz_u / Delta) * tmp
    } else if (length(z) != 1L) {
      ind <- 1L + it %% length(z)
      z[ind] <- -z[ind]
    } else {
      z <- -z
    }
  }
  st$eig <- eig
  nf
}

# The estimate of a step: the user's bound, the power iteration every step,
# or (interval > 0, maybe_maxeig!) the power iteration when the estimate is
# old, on the first step, after a rejection or a discontinuity.
.stab_eig <- function(ctx, st, t, u, fsal, a) {
  nf <- 0L
  e <- a$eigen_est
  if (!is.null(e)) {
    st$eig <- as.numeric(if (is.function(e)) e(u, ctx$p, t) else e)
  } else if (is.null(a$eigen_est_interval)) {
    nf <- .stab_maxeig(ctx, st, t, u, fsal, !isTRUE(a$generic_maxeig), a$maxeig_iters %||% 50L)
  } else if (st$eig_age >= a$eigen_est_interval || ctx$iter <= 1L || st$disc || ctx$EEst_prev > 1) {
    nf <- .stab_maxeig(ctx, st, t, u, fsal, TRUE, 50L)
    st$eig_age <- 1L
  } else {
    st$eig_age <- st$eig_age + 1L
  }
  st$disc <- FALSE
  nf
}

# The cache: the FSAL derivative (replaced on acceptance) and the state that
# persists through rejected steps (st). A new cache after a change of the
# state by a callback marks a discontinuity.
.stab_init <- function(zero_z = FALSE) function(ctx, t, u) {
  st <- new.env(parent = emptyenv())
  st$zprev <- if (zero_z) u * 0 else u
  st$eig <- 1
  st$eig_age <- 0L
  st$disc <- TRUE
  st$mdeg <- 0L
  st$deg_index <- 1L
  st$start <- 1L
  st$internal_deg <- 0L
  st$w0 <- 0
  st$w1 <- 0
  list(fsal = ctx$f(u, ctx$p, t), st = st)
}

.stab_result <- function(u, unew, fsal, fnew, st, nf, err = NULL, EEst = NULL) {
  list(u = unew, err = err, EEst = EEst, k = list(fsal, fnew), cache = list(fsal = fnew, st = st),
       nf = nf, eigen_est = st$eig)
}

# choosedeg!: the first tabulated degree at least mdeg and the start of its
# recurrence coefficients.
.stab_choosedeg <- function(st, ms, mdeg) {
  start <- 1L
  st$mdeg <- mdeg
  for (i in seq_along(ms)) {
    if (ms[i] >= mdeg) {
      st$deg_index <- i
      st$mdeg <- ms[i]
      st$start <- start
      break
    }
    start <- start + ms[i] * 2L - 1L
  }
  invisible(st)
}

# choosedeg_SERK!: the degree, the start of its coefficients and the
# internal degree of the extrapolated methods.
.stab_choosedeg_serk <- function(st, ms, mdeg, name) {
  start <- 1L
  st$mdeg <- mdeg
  for (i in seq_along(ms)) {
    if (ms[i] < mdeg) {
      start <- start + ms[i] + 1L
    } else {
      st$start <- start
      st$mdeg <- ms[i]
      break
    }
  }
  m <- st$mdeg
  if (name == "ESERK5") {
    br <- c(20, 50, 100, 500, 1000, 2000); id <- c(2L, 5L, 10L, 50L, 100L, 200L)
    j <- which(m <= br)
    if (length(j)) st$internal_deg <- id[j[1L]]
  } else if (name == "ESERK4") {
    br <- c(20, 100, 500, 1000, 4000); id <- c(2L, 10L, 25L, 100L, 200L)
    j <- which(m <= br)
    if (length(j)) st$internal_deg <- id[j[1L]]
  } else {
    st$internal_deg <- m %/% 10L
  }
  invisible(st)
}

# The cap of the proposed step (dtnew_modification).
.stab_cap <- function(bound) function(dtnew, eig) min(abs(dtnew), abs(bound(eig)))

.stab_alg <- function(name, order, step, controller = "PI", adaptive_order = NULL, gamma = NULL, qmax = NULL,
                      dtnew_mod = NULL, zero_z = FALSE, ...) {
  a <- .new_alg(name, order, step = NULL, init = .stab_init(zero_z), controller = controller,
                adaptive_order = adaptive_order, gamma = gamma, qmax = qmax, dtnew_mod = dtnew_mod, ...)
  # the step reads the settings of the algorithm
  a$step <- function(ctx, t, u, dt, cache) step(ctx, t, u, dt, cache, a)
  a
}

.rkl_clamp_odd <- function(min_stages, max_stages) {
  lo <- max(3L, as.integer(min_stages))
  if (lo %% 2L == 0L) lo <- lo + 1L
  hi <- max(as.integer(max_stages), lo)
  if (hi %% 2L == 0L) hi <- hi - 1L
  c(lo, max(hi, lo))
}

# ---------------------------------------------------------------------------
# ROCK2 and ROCK4

.rock2_step <- function(ctx, t, u, dt, cache, a) {
  st <- cache$st; fsal <- cache$fsal
  f <- ctx$f; p <- ctx$p
  TB <- .STAB_TABLES$ROCK2
  nf <- .stab_eig(ctx, st, t, u, fsal, a)
  mdeg <- floor(sqrt((1.5 + abs(dt) * st$eig) / 0.811)) + 1
  mdeg <- min(max(mdeg, a$min_stage), a$max_stage)
  .stab_choosedeg(st, TB$ms, max(mdeg, 3) - 2)
  md <- st$mdeg; s0 <- st$start; recf <- TB$recf
  tim1 <- t + dt * recf[s0]
  tim2 <- t + dt * recf[s0]
  tim3 <- t
  uim2 <- u
  uim1 <- u + (dt * recf[s0]) * fsal
  unew <- uim1
  for (i in seq_len(md)[-1L]) {
    mu <- recf[s0 + (i - 2L) * 2L + 1L]; kappa <- recf[s0 + (i - 2L) * 2L + 2L]
    nu <- -1 - kappa
    unew <- f(uim1, p, tim1)
    tim1 <- dt * mu - nu * tim2 - kappa * tim3
    unew <- (dt * mu) * unew - nu * uim1 - kappa * uim2
    if (i < md) {
      uim2 <- uim1
      uim1 <- unew
    }
    tim3 <- tim2
    tim2 <- tim1
  }
  nf <- nf + md - 1L
  # two stage finishing procedure
  dt1 <- dt * TB$fp1[st$deg_index]
  dt2 <- dt * TB$fp2[st$deg_index]
  uim2 <- f(unew, p, tim1)
  uim1 <- unew + dt1 * uim2
  tim1 <- tim1 + dt1
  unew <- f(uim1, p, tim1)
  tmp <- dt2 * (unew - uim2)
  unew <- uim1 + dt1 * unew + tmp
  fnew <- f(unew, p, t + dt)
  .stab_result(u, unew, fsal, fnew, st, nf + 3L, err = tmp)
}

.rock4_step <- function(ctx, t, u, dt, cache, a) {
  st <- cache$st; fsal <- cache$fsal
  f <- ctx$f; p <- ctx$p
  TB <- .STAB_TABLES$ROCK4
  nf <- .stab_eig(ctx, st, t, u, fsal, a)
  mdeg <- floor(sqrt((3 + abs(dt) * st$eig) / 0.353)) + 1
  mdeg <- min(max(mdeg, a$min_stage), a$max_stage)
  .stab_choosedeg(st, TB$ms, max(mdeg, 5) - 4)
  md <- st$mdeg; s0 <- st$start; recf <- TB$recf; di <- st$deg_index
  tim1 <- t + dt * recf[s0]
  tim2 <- t + dt * recf[s0]
  tim3 <- t
  uim2 <- u
  uim1 <- u + (dt * recf[s0]) * fsal
  unew <- uim1
  for (i in seq_len(md)[-1L]) {
    mu <- recf[s0 + (i - 2L) * 2L + 1L]; kappa <- recf[s0 + (i - 2L) * 2L + 2L]
    nu <- -1 - kappa
    unew <- f(uim1, p, tim1)
    tim1 <- dt * mu - nu * tim2 - kappa * tim3
    unew <- (dt * mu) * unew - nu * uim1 - kappa * uim2
    if (i < md) {
      uim2 <- uim1
      uim1 <- unew
    }
    tim3 <- tim2
    tim2 <- tim1
  }
  nf <- nf + md - 1L
  fpa <- TB$fpa[di, ]; fpb <- TB$fpb[di, ]; fpbe <- TB$fpbe[di, ]
  a21 <- dt * fpa[1L]; a31 <- dt * fpa[2L]; a32 <- dt * fpa[3L]
  a41 <- dt * fpa[4L]; a42 <- dt * fpa[5L]; a43 <- dt * fpa[6L]
  B1 <- dt * fpb[1L]; B2 <- dt * fpb[2L]; B3 <- dt * fpb[3L]; B4 <- dt * fpb[4L]
  Bh1 <- dt * (fpbe[1L] - fpb[1L]); Bh2 <- dt * (fpbe[2L] - fpb[2L])
  Bh3 <- dt * (fpbe[3L] - fpb[3L]); Bh4 <- dt * (fpbe[4L] - fpb[4L]); Bh5 <- dt * fpbe[5L]
  # four stage finishing procedure
  uim1 <- f(unew, p, tim1)
  uim2 <- unew + a31 * uim1
  uim3 <- unew + a41 * uim1
  unew <- unew + B1 * uim1
  tmp <- Bh1 * uim1
  uim1 <- unew + (a21 - B1) * uim1
  tim2 <- tim1 + a21
  uim1 <- f(uim1, p, tim2)
  uim2 <- uim2 + a32 * uim1
  uim3 <- uim3 + a42 * uim1
  unew <- unew + B2 * uim1
  tmp <- tmp + Bh2 * uim1
  tim2 <- tim1 + (a31 + a32)
  uim2 <- f(uim2, p, tim2)
  uim3 <- uim3 + a43 * uim2
  unew <- unew + B3 * uim2
  tmp <- tmp + Bh3 * uim2
  tim2 <- tim1 + (a41 + a42 + a43)
  uim3 <- f(uim3, p, tim2)
  unew <- unew + B4 * uim3
  tmp <- tmp + Bh4 * uim3
  fnew <- f(unew, p, t + dt)
  tmp <- tmp + Bh5 * fnew
  .stab_result(u, unew, fsal, fnew, st, nf + 5L, err = tmp)
}

#' @rdname stabilized_rk_methods
#' @export
ROCK2 <- function(min_stages = 0, max_stages = 200, eigen_est = NULL, eigen_est_interval = 25) {
  max_stage <- if (max_stages < 1 || max_stages > 200) 200L else as.integer(max_stages)
  .stab_alg("ROCK2", 2L, .rock2_step, eigen_est = eigen_est, eigen_est_interval = eigen_est_interval,
            min_stage = if (min_stages > max_stage) max_stage else as.integer(min_stages), max_stage = max_stage,
            dtnew_mod = .stab_cap(function(eig) ((min(max_stages, 200)^2) * 0.811 - 1.5) / eig))
}

#' @rdname stabilized_rk_methods
#' @export
ROCK4 <- function(min_stages = 0, max_stages = 152, eigen_est = NULL, eigen_est_interval = 25) {
  max_stage <- if (max_stages < 1 || max_stages > 152) 152L else as.integer(max_stages)
  .stab_alg("ROCK4", 4L, .rock4_step, eigen_est = eigen_est, eigen_est_interval = eigen_est_interval,
            min_stage = if (min_stages > max_stage) max_stage else as.integer(min_stages), max_stage = max_stage,
            dtnew_mod = .stab_cap(function(eig) ((min(max_stages, 152)^2) * 0.353 - 3) / eig))
}

# ---------------------------------------------------------------------------
# RKC, TSRKC2, TSRKC3 and RKMC2

# The Runge-Kutta Chebyshev stages of RKC (and of the RKC steps of TSRKC3).
.rkc_core <- function(f, p, t, u, dt, fsal, mdeg) {
  w0m1 <- (2 / 13) / (mdeg^2)
  w0 <- 1 + w0m1
  w0sqm1 <- w0m1 * (w0m1 + 2)
  temp <- sqrt(w0sqm1)
  arg <- mdeg * log(w0 + temp)
  w1 <- (sinh(arg) * w0sqm1) / (cosh(arg) * mdeg * temp - w0 * sinh(arg))
  b1 <- 1 / ((2 * w0)^2)
  b2 <- b1
  gprev2 <- u
  mus <- w1 * b1
  gprev <- u + dt * mus * fsal
  unew <- gprev
  th2 <- 0; th1 <- mus
  z1 <- w0; z2 <- 1; dz1 <- 1; dz2 <- 0; d2z1 <- 0; d2z2 <- 0
  for (iter in seq_len(mdeg)[-1L]) {
    z <- 2 * w0 * z1 - z2
    dz <- 2 * w0 * dz1 - dz2 + 2 * z1
    d2z <- 2 * w0 * d2z1 - d2z2 + 4 * dz1
    b <- d2z / (dz^2)
    nus <- 1 - z1 * b1
    mu <- 2 * w0 * b / b1
    nu <- -b / b2
    mus <- mu * w1 / w0
    unew <- f(gprev, p, t + dt * th1)
    unew <- mu * gprev + nu * gprev2 + (1 - mu - nu) * u + dt * mus * (unew - nus * fsal)
    th <- mu * th1 + nu * th2 + mus * (1 - nus)
    if (iter < mdeg) {
      gprev2 <- gprev; gprev <- unew
      th2 <- th1; th1 <- th
      b2 <- b1; b1 <- b
      z2 <- z1; z1 <- z
      dz2 <- dz1; dz1 <- dz
      d2z2 <- d2z1; d2z1 <- d2z
    }
  }
  unew
}

.rkc_step <- function(ctx, t, u, dt, cache, a) {
  st <- cache$st; fsal <- cache$fsal
  f <- ctx$f; p <- ctx$p
  nf <- .stab_eig(ctx, st, t, u, fsal, a)
  mdeg <- floor(sqrt(1.54 * abs(dt) * st$eig + 1)) + 1
  unew <- .rkc_core(f, p, t, u, dt, fsal, mdeg)
  fnew <- f(unew, p, t + dt)
  tmp <- (4 * (u - unew) + 2 * dt * (fsal + fnew)) / 5
  .stab_result(u, unew, fsal, fnew, st, nf + mdeg, err = tmp)
}

#' @rdname stabilized_rk_methods
#' @export
RKC <- function(eigen_est = NULL, eigen_est_interval = 25) {
  .stab_alg("RKC", 2L, .rkc_step, controller = "Predictive", adaptive_order = 2, gamma = 8 / 10,
            eigen_est = eigen_est, eigen_est_interval = eigen_est_interval)
}

.tsrkc2_step <- function(ctx, t, u, dt, cache, a) {
  st <- cache$st; fsal <- cache$fsal
  f <- ctx$f; p <- ctx$p
  q <- (t - ctx$tprev) / dt
  onemq <- 1 - q
  nf <- .stab_eig(ctx, st, t, u, fsal, a)
  mdeg <- floor(sqrt(1 + 0.759782816506459 * abs(dt) * st$eig *
                       (onemq + sqrt(1 + q * (q - 0.598626091572911))))) + 1
  tsw0 <- 1.1
  acoshtsw0 <- acosh(tsw0)
  sinhacoshtsw0 <- sinh(acoshtsw0)
  w0m1 <- 2 * (sinh(acoshtsw0 / (2 * mdeg))^2)
  w0 <- 1 + w0m1
  w0sqm1 <- w0m1 * (w0m1 + 2)
  dtsw0 <- mdeg * sinhacoshtsw0 / sqrt(w0sqm1)
  d2tsw0 <- ((mdeg^2) * tsw0 - w0 * dtsw0) / w0sqm1
  w1 <- (onemq * dtsw0 + sqrt((onemq * dtsw0)^2 + 4 * q * tsw0 * d2tsw0)) / (2 * d2tsw0)
  gprev2 <- u
  mus <- w1 / w0
  gprev <- u + dt * mus * fsal
  unew <- gprev
  th2 <- 0; th1 <- mus
  z1 <- w0; z2 <- 1
  z <- 2 * w0 * z1 - z2
  for (iter in seq_len(mdeg)[-1L]) {
    mu <- 2 * w0 * z1 / z
    nu <- -z2 / z
    mus <- mu * w1 / w0
    unew <- f(gprev, p, t + dt * th1)
    unew <- mu * gprev + nu * gprev2 + dt * mus * unew
    if (iter < mdeg) {
      gprev2 <- gprev; gprev <- unew
      th <- mu * th1 + nu * th2 + mus
      th2 <- th1; th1 <- th
      z2 <- z1; z1 <- z
      z <- 2 * w0 * z1 - z2
    }
  }
  g <- (1 + q) * tsw0 / (q * tsw0 + w1 * dtsw0)
  unew <- (1 - g) * ctx$uprev2 + g * unew
  fnew <- f(unew, p, t + dt)
  tmp <- (u - unew + dt * fnew) / 3
  .stab_result(u, unew, fsal, fnew, st, nf + mdeg, err = tmp)
}

.tsrkc3_step <- function(ctx, t, u, dt, cache, a) {
  st <- cache$st; fsal <- cache$fsal
  f <- ctx$f; p <- ctx$p
  uprev2 <- ctx$uprev2
  q <- (t - ctx$tprev) / dt
  onemq <- 1 - q
  onepq <- 1 + q
  onepq2 <- onepq * onepq
  # the first and possibly the second steps are RKC steps
  rkcstep <- q < 0.49
  nf <- .stab_eig(ctx, st, t, u, fsal, a)
  dtsw0 <- 0; d2tsw0 <- 0
  if (rkcstep) {
    mdeg <- floor(sqrt(1.54 * abs(dt) * st$eig + 1)) + 1
    w0m1 <- (2 / 13) / (mdeg^2)
    w0 <- 1 + w0m1
    w0sqm1 <- w0m1 * (w0m1 + 2)
    temp <- sqrt(w0sqm1)
    arg <- mdeg * log(w0 + temp)
    w1 <- (sinh(arg) * w0sqm1) / (cosh(arg) * mdeg * temp - w0 * sinh(arg))
    b1 <- 1 / ((2 * w0)^2)
    b <- b1
    b2 <- b1
  } else {
    mdeg <- floor(sqrt(4 + 1.267029788142009 * abs(dt) * st$eig *
                         (onemq + sqrt(1 + q * (0.44256220745562963 + q))))) + 1
    mdeg2 <- mdeg^2
    tsw0 <- 1.25
    acoshtsw0 <- acosh(tsw0)
    sinhacoshtsw0 <- sinh(acoshtsw0)
    acoshtsw0dm <- acoshtsw0 / mdeg
    w0m1 <- 2 * (sinh(acoshtsw0dm / 2)^2)
    w0 <- 1 + w0m1
    w0sq <- w0^2
    w0sqm1 <- w0m1 * (w0m1 + 2)
    dtsw0 <- mdeg * sinhacoshtsw0 / sqrt(w0sqm1)
    d2tsw0 <- (mdeg2 * tsw0 - w0 * dtsw0) / w0sqm1
    d3tsw0 <- ((1 + 2 * w0sq + mdeg2 * w0sqm1) * dtsw0 - 3 * mdeg2 * w0 * tsw0) / (w0sqm1^2)
    w1 <- (onemq * d2tsw0 + sqrt((onemq * d2tsw0)^2 + 4 * q * dtsw0 * d3tsw0)) / (2 * d3tsw0)
    b1 <- sinh((mdeg - 2) * acoshtsw0dm) / (4 * sinh((mdeg - 1) * acoshtsw0dm))
    b <- 15 / ((8 * w0)^2)
    b2 <- b1
  }
  gprev2 <- u
  mus <- w1 * b1
  gprev <- u + dt * mus * fsal
  unew <- gprev
  th2 <- 0; th1 <- mus
  z1 <- w0; z2 <- 1; dz1 <- 1; dz2 <- 0; d2z1 <- 0; d2z2 <- 0
  z <- 2 * w0 * z1 - z2
  dz <- 2 * w0 * dz1 - dz2 + 2 * z1
  d2z <- 2 * w0 * d2z1 - d2z2 + 4 * dz1
  for (iter in seq_len(mdeg)[-1L]) {
    nus <- 1 - z1 * b1
    mu <- 2 * w0 * b / b1
    nu <- -b / b2
    mus <- mu * w1 / w0
    unew <- f(gprev, p, t + dt * th1)
    unew <- mu * gprev + nu * gprev2 + (1 - mu - nu) * u + dt * mus * (unew - nus * fsal)
    if (iter < mdeg) {
      gprev2 <- gprev; gprev <- unew
      th <- mu * th1 + nu * th2 + mus * (1 - nus)
      th2 <- th1; th1 <- th
      b2 <- b1; b1 <- b
      z2 <- z1; z1 <- z
      dz2 <- dz1; dz1 <- dz
      d2z2 <- d2z1; d2z1 <- d2z
      z <- 2 * w0 * z1 - z2
      dz <- 2 * w0 * dz1 - dz2 + 2 * z1
      d2z <- 2 * w0 * d2z1 - d2z2 + 4 * dz1
      b <- d2z / (dz^2)
    }
  }
  if (!rkcstep) {
    aa <- onepq / (q * dtsw0 + w1 * d2tsw0)
    g <- (aa * dtsw0 - 1) / q
    bb <- aa / (b * w1)
    unew <- (1 - g - bb) * u + g * uprev2 + bb * unew
  }
  fnew <- f(unew, p, t + dt)
  tmp <- if (rkcstep) (3 / 5) * (2 * (u - unew) + dt * (fsal + fnew)) else
    (3 / 5) * (u / q - uprev2 / (q * onepq2) - unew * (2 + q) / onepq2 + dt * fnew / onepq)
  .stab_result(u, unew, fsal, fnew, st, nf + mdeg, err = tmp)
}

#' @rdname stabilized_rk_methods
#' @export
TSRKC2 <- function(eigen_est = NULL) {
  .stab_alg("TSRKC2", 2L, .tsrkc2_step, controller = "Predictive", adaptive_order = 1, gamma = 8 / 10,
            eigen_est = eigen_est)
}

#' @rdname stabilized_rk_methods
#' @export
TSRKC3 <- function(eigen_est = NULL) {
  .stab_alg("TSRKC3", 3L, .tsrkc3_step, controller = "Predictive", adaptive_order = 2, gamma = 8 / 10,
            qmax = 2, eigen_est = eigen_est)
}

.rkmc2_step <- function(ctx, t, u, dt, cache, a) {
  st <- cache$st; fsal <- cache$fsal
  f <- ctx$f; p <- ctx$p
  nf <- .stab_eig(ctx, st, t, u, fsal, a)
  mdeg <- max(3, ceiling(-0.8306782178712795 + 1.8547887825836553 * (abs(dt) * st$eig)^0.533871357807877))
  mdeg <- min(max(mdeg, a$min_stages), a$max_stages)
  if (mdeg != st$mdeg) {
    # w0 = cosh(alpha) shifts the Chebyshev polynomial to be monotone: bisection on alpha
    st$mdeg <- mdeg
    s <- mdeg
    sign_s <- if (s %% 2 == 0) 1 else -1
    alo <- 1e-10; ahi <- 2
    for (k in 1:50) {
      al <- (alo + ahi) / 2
      Ts <- cosh(s * al)
      Tsm1 <- cosh((s - 1) * al)
      Tsm2 <- cosh((s - 2) * al)
      dTsm1 <- (s - 1) * sinh((s - 1) * al) / sinh(al)
      fval <- 1 + sign_s / (s * (s - 2)) + cosh(al) + Ts / (2 * s) - Tsm2 / (2 * (s - 2)) - (1 + Tsm1)^2 / dTsm1
      if (fval > 0) alo <- al else ahi <- al
    }
    al <- (alo + ahi) / 2
    Tsm1 <- cosh((s - 1) * al)
    dTsm1 <- (s - 1) * sinh((s - 1) * al) / sinh(al)
    st$w0 <- cosh(al)
    st$w1 <- (1 + Tsm1) / dTsm1
  }
  w0 <- st$w0; w1 <- st$w1
  Tjm2 <- 1; Tjm1 <- w0
  bjm2 <- 1 / (1 + Tjm2)
  bjm1 <- 1 / (1 + Tjm1)
  gprev2 <- u
  mus <- bjm1 * w1
  gprev <- u + dt * mus * fsal
  unew <- gprev
  th2 <- 0; th1 <- mus
  bs <- bjm1
  for (j in seq_len(mdeg)[-1L]) {
    Tj <- 2 * w0 * Tjm1 - Tjm2
    bj <- 1 / (1 + Tj)
    bs <- bj
    mu <- 2 * w0 * bj / bjm1
    nu <- -bj / bjm2
    mut <- 2 * w1 * bj / bjm1
    unew <- f(gprev, p, t + dt * th1)
    unew <- (1 - mu - nu) * u + mu * gprev + nu * gprev2 + dt * mut * (unew - bjm1 * fsal)
    th <- mu * th1 + nu * th2 + mut * (1 - bjm1)
    if (j < mdeg) {
      gprev2 <- gprev; gprev <- unew
      th2 <- th1; th1 <- th
      Tjm2 <- Tjm1; Tjm1 <- Tj
      bjm2 <- bjm1; bjm1 <- bj
    }
  }
  gs <- bjm1 / (2 * mdeg * w1)
  ds <- -bjm1 / (2 * (mdeg - 2) * w1)
  unew <- (1 - gs / bs - ds / bjm2) * u + (gs / bs) * unew + (ds / bjm2) * gprev2 + dt * bjm1 * fsal
  fnew <- f(unew, p, t + dt)
  tmp <- (u - unew + dt * fnew) / 10
  .stab_result(u, unew, fsal, fnew, st, nf + mdeg, err = tmp)
}

#' @rdname stabilized_rk_methods
#' @export
RKMC2 <- function(min_stages = 3, max_stages = 1000, eigen_est = NULL) {
  min_stages <- max(3L, as.integer(min_stages))
  .stab_alg("RKMC2", 2L, .rkmc2_step, controller = "Predictive", adaptive_order = 2, gamma = 8 / 10,
            eigen_est = eigen_est, min_stages = min_stages, max_stages = as.integer(max_stages),
            generic_maxeig = TRUE, zero_z = TRUE,
            dtnew_mod = .stab_cap(function(eig) 0.31 * (min(max_stages, 1000) + 0.83)^1.87 / eig))
}

# ---------------------------------------------------------------------------
# ESERK4, ESERK5 and SERK2

.eserk_step <- function(name, nsub, scale, cap, degree, alpha) function(ctx, t, u, dt, cache, a) {
  st <- cache$st; fsal <- cache$fsal
  f <- ctx$f; p <- ctx$p
  TB <- .STAB_TABLES[[name]]
  nf <- .stab_eig(ctx, st, t, u, fsal, a)
  mdeg <- degree(abs(dt), st$eig)
  if (mdeg > cap) mdeg <- cap
  .stab_choosedeg_serk(st, TB$ms, mdeg, name)
  mdeg <- st$mdeg; start <- st$start; id <- st$internal_deg
  Bi <- TB$Bi
  al <- alpha(mdeg)
  unew <- u * 0
  tmp <- u * 0
  for (i in seq_len(nsub)) {
    hi <- dt / i
    ti <- t
    Si <- u * 0
    uim1 <- u
    uim2 <- u * 0
    for (j in seq_len(i)) {
      r <- ti
      Si <- Bi[start] * uim1
      for (s in seq_len(mdeg)) {
        k <- f(uim1, p, r)
        ui <- if (s %% id == 1L) uim1 + al * hi * k else 2 * uim1 - uim2 + 2 * al * hi * k
        q <- s %/% id
        r <- ti + al * (s^2 + q * id^2) * hi
        Si <- Si + Bi[start + s] * ui
        if (s < mdeg) {
          uim2 <- uim1
          uim1 <- ui
        }
      }
      if (j < i) {
        ti <- ti + hi
        uim1 <- Si
      }
    }
    unew <- unew + TB$Cu[i] * Si
    tmp <- tmp + TB$Ce[i] * Si
  }
  nf <- nf + mdeg * ((nsub * (nsub + 1L)) %/% 2L)
  unew <- unew / scale
  tmp <- tmp / scale
  fnew <- f(unew, p, t + dt)
  .stab_result(u, unew, fsal, fnew, st, nf + 1L, err = tmp)
}

#' @rdname stabilized_rk_methods
#' @export
ESERK4 <- function(eigen_est = NULL) {
  .stab_alg("ESERK4", 4L, .eserk_step("ESERK4", 4L, 6, 4000, function(h, e) floor(sqrt(h * e)) + 1,
                                                 function(m) 2 / m^2), eigen_est = eigen_est,
            maxeig_iters = 100L, dtnew_mod = .stab_cap(function(eig) 0.98 * 4000 * 4000 / eig))
}

#' @rdname stabilized_rk_methods
#' @export
ESERK5 <- function(eigen_est = NULL) {
  .stab_alg("ESERK5", 5L, .eserk_step("ESERK5", 5L, 24, 2000, function(h, e) floor(sqrt(h * e / 0.98)) + 1,
                                                 function(m) 100 / (49 * m^2)), eigen_est = eigen_est,
            maxeig_iters = 100L, dtnew_mod = .stab_cap(function(eig) 0.98 * 2000 * 2000 / eig))
}

.serk2_step <- function(ctx, t, u, dt, cache, a) {
  st <- cache$st; fsal <- cache$fsal
  f <- ctx$f; p <- ctx$p
  TB <- .STAB_TABLES$SERK2
  nf <- .stab_eig(ctx, st, t, u, fsal, a)
  mdeg <- floor(sqrt(abs(dt) * st$eig / 0.8)) + 1
  if (mdeg > 250) mdeg <- 250
  .stab_choosedeg_serk(st, TB$ms, mdeg, "SERK2")
  mdeg <- st$mdeg; start <- st$start; id <- st$internal_deg
  Bi <- TB$Bi
  al <- (5 / 2) / mdeg^2
  uim1 <- u
  uim2 <- u
  Si <- Bi[start] * u
  # the stage produced at index j of block i sits at
  # t + (j^2 + (i - 1) internal_deg^2) alpha dt, so f is evaluated one stage back
  for (i in 1:10) {
    k <- f(uim1, p, t + ((i - 1) * id^2) * al * dt)
    unew <- uim1 + al * dt * k
    Si <- Si + Bi[start + (i - 1) * id + 1] * unew
    uim2 <- uim1
    uim1 <- unew
    for (j in seq_len(id)[-1L]) {
      k <- f(uim1, p, t + ((j - 1)^2 + (i - 1) * id^2) * al * dt)
      unew <- 2 * uim1 - uim2 + 2 * al * dt * k
      Si <- Si + Bi[start + j + (i - 1) * id] * unew
      if (j * i < mdeg) {
        uim2 <- uim1
        uim1 <- unew
      }
    }
  }
  nf <- nf + 10L * id
  unew <- Si
  fnew <- f(unew, p, t + dt)
  tmp <- unew - u - dt * fnew
  .stab_result(u, unew, fsal, fnew, st, nf + 1L, err = tmp)
}

#' @rdname stabilized_rk_methods
#' @export
SERK2 <- function(eigen_est = NULL) {
  .stab_alg("SERK2", 2L, .serk2_step, controller = "Predictive", adaptive_order = 1, eigen_est = eigen_est,
            maxeig_iters = 100L, dtnew_mod = .stab_cap(function(eig) 0.8 * 250 * 250 / (eig + 1)))
}

# ---------------------------------------------------------------------------
# RKL1, RKL2, RKG1 and RKG2

# The stages shared by the super time stepping methods: coefficient
# functions mu(j), nu(j), mut(j), gamt(j) (gamt NULL for the first order
# methods), the first stage factor mu1 and the number of stages s.
.sts_stages <- function(f, p, t, u, dt, fsal, s, mu1, cf) {
  uim2 <- u
  uim1 <- u + (dt * mu1) * fsal
  cjm2 <- 0; cjm1 <- mu1
  unew <- uim1
  for (j in seq_len(s)[-1L]) {
    co <- cf(j)
    fy <- f(uim1, p, t + dt * cjm1)
    if (is.null(co$gt)) {
      unew <- co$mu * uim1 + co$nu * uim2 + co$mt * dt * fy
      cj <- co$mu * cjm1 + co$nu * cjm2 + co$mt
    } else {
      unew <- co$mu * uim1 + co$nu * uim2 + (1 - co$mu - co$nu) * u + dt * co$mt * fy + dt * co$gt * fsal
      cj <- co$mu * cjm1 + co$nu * cjm2 + co$mt + co$gt
    }
    uim2 <- uim1; uim1 <- unew
    cjm2 <- cjm1; cjm1 <- cj
  }
  unew
}

.sts_step <- function(stages, coefs) function(ctx, t, u, dt, cache, a) {
  st <- cache$st; fsal <- cache$fsal
  f <- ctx$f; p <- ctx$p
  nf <- .stab_eig(ctx, st, t, u, fsal, a)
  s <- stages(abs(dt), st$eig, a)
  st$mdeg <- s
  cf <- coefs(s)
  unew <- .sts_stages(f, p, t, u, dt, fsal, s, cf$mu1, cf$fn)
  tmp <- unew - (u + dt * fsal)
  fnew <- f(unew, p, t + dt)
  EEst <- if (isTRUE(ctx$adaptive)) .err_norm(tmp, u, unew, ctx$abstol, ctx$reltol) / s
  .stab_result(u, unew, fsal, fnew, st, nf + s, EEst = EEst)
}

.rkl_odd_stages <- function(x, a) {
  s <- max(x, a$min_stage)
  if (s %% 2 == 0) s <- s + 1
  min(s, a$max_stage)
}

.rkl1_coefs <- function(s) {
  w1 <- 2 / (s^2 + s)
  list(mu1 = w1, fn = function(j) {
    mu <- (2 * j - 1) / j
    list(mu = mu, nu = -(j - 1) / j, mt = mu * w1)
  })
}

.rkl2_b <- function(j) if (j <= 2) 1 / 3 else (j * j + j - 2) / (2 * j * (j + 1))

.rkl2_coefs <- function(s) {
  w1 <- 4 / (s^2 + s - 2)
  list(mu1 = (1 / 3) * w1, fn = function(j) {
    bj <- .rkl2_b(j); bjm1 <- .rkl2_b(j - 1); bjm2 <- .rkl2_b(j - 2)
    mu <- (2 * j - 1) / j * bj / bjm1
    mt <- mu * w1
    list(mu = mu, nu = -(j - 1) / j * bj / bjm2, mt = mt, gt = -(1 - bjm1) * mt)
  })
}

.rkg1_b <- function(j) 2 / ((j + 1) * (j + 2))

.rkg1_coefs <- function(s) {
  w1 <- 4 / (s * (s + 3))
  list(mu1 = (3 * .rkg1_b(1) / .rkg1_b(0)) * w1, fn = function(j) {
    bj <- .rkg1_b(j)
    mu <- (2 * j + 1) / j * bj / .rkg1_b(j - 1)
    list(mu = mu, nu = -(j + 1) / j * bj / .rkg1_b(j - 2), mt = mu * w1)
  })
}

.rkg2_coefs <- function(s) {
  w1 <- 6 / ((s + 4) * (s - 1))
  list(mu1 = w1, fn = function(j) {
    bj <- (4 * (j - 1) * (j + 4)) / (3 * j * (j + 1) * (j + 2) * (j + 3))
    bjm1 <- if (j == 2) 1 / 3 else (4 * (j - 2) * (j + 3)) / (3 * (j - 1) * j * (j + 1) * (j + 2))
    bjm2 <- if (j == 2) 1 else if (j == 3) 1 / 3 else (4 * (j - 3) * (j + 2)) / (3 * (j - 2) * (j - 1) * j * (j + 1))
    ajm1 <- 1 - ((j * (j + 1)) %/% 2) * bjm1
    mu <- (2 * j + 1) / j * bj / bjm1
    mt <- mu * w1
    list(mu = mu, nu = -(j + 1) / j * bj / bjm2, mt = mt, gt = -mt * ajm1)
  })
}

#' @rdname stabilized_rk_methods
#' @export
RKL1 <- function(min_stages = 3, max_stages = 200, eigen_est = NULL) {
  ms <- .rkl_clamp_odd(min_stages, max_stages)
  .stab_alg("RKL1", 1L, .sts_step(function(h, e, a) .rkl_odd_stages(ceiling((sqrt(1 + 4 * h * e) - 1) / 2), a),
                                  .rkl1_coefs),
            controller = "Predictive", adaptive_order = 1, gamma = 8 / 10, eigen_est = eigen_est,
            min_stage = ms[1L], max_stage = ms[2L], zero_z = TRUE,
            dtnew_mod = .stab_cap(function(eig) (ms[2L]^2 + ms[2L]) / eig))
}

#' @rdname stabilized_rk_methods
#' @export
RKL2 <- function(min_stages = 3, max_stages = 200, eigen_est = NULL) {
  ms <- .rkl_clamp_odd(min_stages, max_stages)
  .stab_alg("RKL2", 2L, .sts_step(function(h, e, a) .rkl_odd_stages(ceiling((sqrt(9 + 8 * h * e) - 1) / 2), a),
                                  .rkl2_coefs),
            controller = "Predictive", adaptive_order = 2, gamma = 8 / 10, eigen_est = eigen_est,
            min_stage = ms[1L], max_stage = ms[2L], zero_z = TRUE,
            dtnew_mod = .stab_cap(function(eig) (ms[2L]^2 + ms[2L] - 2) / (2 * eig)))
}

#' @rdname stabilized_rk_methods
#' @export
RKG1 <- function(min_stages = 2, max_stages = 200, eigen_est = NULL) {
  lo <- max(2L, as.integer(min_stages)); hi <- max(as.integer(max_stages), lo)
  .stab_alg("RKG1", 1L, .sts_step(function(h, e, a) min(max(ceiling((sqrt(9 + 16 * h * e) - 3) / 2), a$min_stage),
                                                     a$max_stage), .rkg1_coefs),
            controller = "Predictive", adaptive_order = 1, gamma = 8 / 10, eigen_est = eigen_est,
            min_stage = lo, max_stage = hi, zero_z = TRUE,
            dtnew_mod = .stab_cap(function(eig) hi * (hi + 3) / (2 * eig)))
}

#' @rdname stabilized_rk_methods
#' @export
RKG2 <- function(min_stages = 3, max_stages = 200, eigen_est = NULL) {
  lo <- max(3L, as.integer(min_stages)); hi <- max(as.integer(max_stages), lo)
  .stab_alg("RKG2", 2L, .sts_step(function(h, e, a) min(max(ceiling((sqrt(25 + 24 * h * e) - 3) / 2), a$min_stage),
                                                     a$max_stage), .rkg2_coefs),
            controller = "Predictive", adaptive_order = 2, gamma = 8 / 10, eigen_est = eigen_est,
            min_stage = lo, max_stage = hi, zero_z = TRUE,
            dtnew_mod = .stab_cap(function(eig) (hi + 4) * (hi - 1) / (3 * eig)))
}
