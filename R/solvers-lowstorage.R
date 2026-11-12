# The low storage Runge-Kutta methods of OrdinaryDiffEqLowStorageRK in the
# out of place form of low_storage_rk_perform_step.jl (the ConstantCache
# methods): the 2N, 2C, 3S, 3S* and 3S* FSAL families, the 2R+ to 5R+
# families and the two SHLDDRK schemes, with the coefficients of their
# caches (.LSRK in solvers-lowstorage-tableaus.R, from
# tools/lowstorage-probe.jl). The multiply adds of the vector arithmetic are
# unfused, as for Julia arrays under @muladd; the stage times t + c dt are
# fused. Like the out of place methods of OrdinaryDiffEq, they do not apply a
# stage limiter (the step limiter is applied by the driver to every method).
#
# The dense output is the cubic Hermite interpolant. The 3S* FSAL and SHLDDRK
# methods keep both end derivatives of the step; for the others OrdinaryDiffEq
# keeps one and its generic _ode_addsteps! evaluates both ends again on first
# use (not counted in the statistics), which is done here too.
#
# The R+ forms of Julia start a step from integrator.u, which after a rejected
# step still holds the rejected result; the port starts every step from the
# accepted state (see PORTING.md).

# the k of a step: the derivatives at both ends, evaluated on demand
.lsrk_k <- function(ctx, t, dt, f0 = NULL, f1 = NULL) {
  k <- new.env(parent = emptyenv())
  k$f0 <- f0
  k$f1 <- f1
  k$f <- ctx$f
  k$p <- ctx$p
  k$t0 <- t
  k$t1 <- t + dt
  k
}

.lsrk_ends <- function(k, y0, y1) {
  if (is.null(k$f1)) {
    k$f0 <- k$f(y0, k$p, k$t0)
    k$f1 <- k$f(y1, k$p, k$t1)
  }
  list(k$f0, k$f1)
}

.lsrk_interp <- function(theta, dt, y0, y1, k) .hermite_interp(theta, dt, y0, y1, .lsrk_ends(k, y0, y1))

.lsrk_interp_fused <- function(theta, dt, y0, y1, k) {
  .hermite_interp_fused(theta, dt, y0, y1, .lsrk_ends(k, y0, y1))
}

# a zero vector shaped like x (zero(x) of Julia)
.lsrk_zero <- function(x) {
  x[] <- 0
  x
}

# perform_step! of the ConstantCache of family `fam` with coefficients `co`
.lsrk_step <- function(fam, co) {
  switch(fam,
    LowStorageRK2NConstantCache = function(ctx, t, u, dt, cache) {
      f <- ctx$f; p <- ctx$p
      tmp <- dt * cache$fsal
      u <- u + co$B1 * tmp
      for (i in seq_along(co$A2end)) {
        k <- f(u, p, .fma(co$c2end[i], dt, t))
        tmp <- co$A2end[i] * tmp + dt * k
        u <- u + co$B2end[i] * tmp
      }
      f1 <- f(u, p, t + dt)
      list(u = u, k = .lsrk_k(ctx, t, dt), cache = list(fsal = f1), nf = length(co$A2end) + 1L)
    },
    LowStorageRK2CConstantCache = function(ctx, t, u, dt, cache) {
      f <- ctx$f; p <- ctx$p
      k <- f0 <- f(u, p, t)
      u <- u + (co$B1 * dt) * k
      for (i in seq_along(co$A2end)) {
        tmp <- u + (co$A2end[i] * dt) * k
        k <- f(tmp, p, .fma(co$c2end[i], dt, t))
        u <- u + (co$B2end[i] * dt) * k
      }
      # the method is FSAL in OrdinaryDiffEq although every step evaluates
      # its first derivative afresh: the stored value is unused, but a step
      # limiter evaluates it again (and counts it) as there
      list(u = u, k = .lsrk_k(ctx, t, dt), cache = list(fsal = f0), nf = length(co$A2end) + 1L)
    },
    LowStorageRK3SConstantCache = function(ctx, t, u, dt, cache) {
      f <- ctx$f; p <- ctx$p
      uprev <- u
      tmp <- u
      u <- tmp + (co$beta1 * dt) * cache$fsal
      for (i in seq_along(co$gamma12end)) {
        k <- f(u, p, .fma(co$c2end[i], dt, t))
        tmp <- tmp + co$delta2end[i] * u
        u <- co$gamma12end[i] * u + co$gamma22end[i] * tmp + co$gamma32end[i] * uprev + (co$beta2end[i] * dt) * k
      }
      f1 <- f(u, p, t + dt)
      list(u = u, k = .lsrk_k(ctx, t, dt), cache = list(fsal = f1), nf = length(co$gamma12end) + 1L)
    },
    LowStorageRK3SpConstantCache = ,
    LowStorageRK3SpFSALConstantCache = function(ctx, t, u, dt, cache) {
      f <- ctx$f; p <- ctx$p
      fsal <- fam == "LowStorageRK3SpFSALConstantCache"
      uprev <- u
      nf <- length(co$gamma12end)
      if (fsal) f0 <- cache$fsal else {
        f0 <- f(uprev, p, t)
        nf <- nf + 1L
      }
      tmp <- uprev
      u <- tmp + (co$beta1 * dt) * f0
      utilde <- (co$bhat1 * dt) * f0
      for (i in seq_along(co$gamma12end)) {
        k <- f(u, p, .fma(co$c2end[i], dt, t))
        tmp <- tmp + co$delta2end[i] * u
        u <- co$gamma12end[i] * u + co$gamma22end[i] * tmp + co$gamma32end[i] * uprev + (co$beta2end[i] * dt) * k
        utilde <- utilde + (co$bhat2end[i] * dt) * k
      }
      if (!fsal) return(list(u = u, err = utilde, k = .lsrk_k(ctx, t, dt), cache = list(), nf = nf))
      f1 <- f(u, p, t + dt)
      utilde <- utilde + (co$bhatfsal * dt) * f1
      list(u = u, err = utilde, k = .lsrk_k(ctx, t, dt, f0, f1), cache = list(fsal = f1), nf = nf + 1L)
    },
    LowStorageRK2RPConstantCache = ,
    LowStorageRK3RPConstantCache = ,
    LowStorageRK4RPConstantCache = ,
    LowStorageRK5RPConstantCache = function(ctx, t, u, dt, cache) {
      f <- ctx$f; p <- ctx$p
      # r registers: Ai1 to Ai(r - 1) weight the derivatives of the stages
      # i - 1, i - 2, ..., i - r + 1, added to the state of stage i - r + 1
      r <- switch(fam, LowStorageRK2RPConstantCache = 2L, LowStorageRK3RPConstantCache = 3L,
                  LowStorageRK4RPConstantCache = 4L, 5L)
      A <- if (r == 2L) list(co$Ai) else co[paste0("Ai", seq_len(r - 1L))]
      uprev <- u
      k <- cache$fsal
      # f_{i-2}, f_{i-3}, ... (zero at the start) and u_{i-1}, u_{i-2}, ...
      # (uprev at the start)
      fs <- rep(list(.lsrk_zero(k)), r - 2L)
      us <- rep(list(uprev), r - 1L)
      tmp <- .lsrk_zero(uprev)
      for (i in seq_along(A[[1L]])) {
        tmp <- tmp + ((co$Bi[i] - co$Bhi[i]) * dt) * k
        if (r == 2L) {
          gprev <- u + (A[[1L]][i] * dt) * k
        } else {
          s <- A[[1L]][i] * k
          for (j in seq_len(r - 2L)) s <- s + A[[j + 1L]][i] * fs[[j]]
          gprev <- us[[r - 1L]] + s * dt
        }
        u <- u + (co$Bi[i] * dt) * k
        if (r > 2L) {
          fs <- c(list(k), fs)[seq_len(r - 2L)]
          us <- c(list(u), us)[seq_len(r - 1L)]
        }
        k <- f(gprev, p, .fma(co$Ci[i], dt, t))
      }
      tmp <- tmp + ((co$Bl - co$Bhl) * dt) * k
      u <- u + (co$Bl * dt) * k
      f1 <- f(u, p, t + dt)
      list(u = u, err = tmp, k = .lsrk_k(ctx, t, dt), cache = list(fsal = f1), nf = length(A[[1L]]) + 1L)
    },
    SHLDDRK52ConstantCache = ,
    SHLDDRK_2NConstantCache = function(ctx, t, u, dt, cache) {
      f <- ctx$f; p <- ctx$p
      if (fam == "SHLDDRK52ConstantCache") {
        sfx <- ""
        last <- 5L
        nstep <- NULL
      } else {
        # the two schemes alternate; a derivative discontinuity restarts them
        nstep <- if (isTRUE(ctx$derivative_discontinuity)) 1L else cache$step %||% 1L
        sfx <- if (nstep %% 2L == 1L) "1" else "2"
        last <- if (sfx == "1") 5L else 6L
        nstep <- nstep + 1L
      }
      cf <- function(nm, j) co[[paste0(nm, j, sfx)]]
      f0 <- cache$fsal
      tmp <- dt * f0
      u <- u + cf("beta", 1L) * tmp
      for (j in 2:last) {
        tmp <- cf("alpha", j) * tmp + dt * f(u, p, .fma(cf("c", j), dt, t))
        u <- u + cf("beta", j) * tmp
      }
      f1 <- f(u, p, t + dt)
      list(u = u, k = .lsrk_k(ctx, t, dt, f0, f1), cache = list(fsal = f1, step = nstep), nf = last)
    },
    stop("unknown low storage family ", fam, call. = FALSE))
}

# The step, interpolant and initialization of a low storage method of .LSRK.
# initialize! evaluates the derivative at the initial state (counted). The 3S*
# methods evaluate it afresh in every step and have no FSAL derivative to
# evaluate again after an event.
.lsrk_parts <- function(name) {
  co <- .LSRK[[name]]
  parts <- list(step = .lsrk_step(co$family, co), interp = .lsrk_interp, init = .fsal_init)
  if (co$family == "LowStorageRK3SpConstantCache") {
    parts$init <- function(ctx, t, u) {
      ctx$f(u, ctx$p, t)
      list()
    }
    parts$on_event <- function(ctx, t, u, cache) list(cache = cache, nf = 0L)
  }
  parts
}
