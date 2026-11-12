# Runge-Kutta-Nystrom methods of OrdinaryDiffEqRKN (algorithms.jl,
# alg_utils.jl, rkn_caches.jl, rkn_perform_step.jl, interpolants.jl; the
# tableaus of rkn_tableaus.jl are in R/solvers-rkn-tableaus.R) for second
# order problems u'' = f1(u', u, p, t), given as SecondOrderODEProblem or
# DynamicalODEProblem (state c(v, u), velocities first).
#
# The velocity independent methods (Nystrom4VelocityIndependent,
# Nystrom5VelocityIndependent, DPRKN*, ERKN*, IRKN*) hold the velocity fixed
# at its start value in every stage; at initialization they call f1 once with
# a sentinel velocity that stops on any use, as Julia's
# RKNFirstDerivativeSentinel probe. The velocity dependent methods (Nystrom4,
# FineRKN4, FineRKN5, RKN4) also carry the velocity through the stages (abar).
# Each step keeps c(f1, f2) at its end for the next step and for cubic Hermite
# dense output; DPRKN6 has its own sixth order interpolant from its stages.
# sol$stats$nf counts the f1 evaluations and sol$stats$nf2 the f2 ones, with
# Julia's conventions (tab.nf_per_step for the velocity dependent methods).

#' Runge-Kutta-Nystrom methods
#'
#' Runge-Kutta-Nystrom methods of OrdinaryDiffEqRKN for second order problems
#' `u'' = f1(u', u, p, t)` given as [SecondOrderODEProblem()] or
#' [DynamicalODEProblem()] (state `c(v, u)`).
#'
#' * Velocity dependent (`f1` may use `v`): `Nystrom4()` (fourth order,
#'   fixed step), `RKN4()` (fourth order, fixed step, Nystrom form of the
#'   classical method), `FineRKN4()` and `FineRKN5()` (adaptive embedded
#'   pairs of order 4 and 5 by Fine).
#' * Velocity independent (`f1` must not use `v`; the methods stop at
#'   initialization if it does): `Nystrom4VelocityIndependent()` and
#'   `Nystrom5VelocityIndependent()` (fixed step), `IRKN3()` and `IRKN4()`
#'   (fixed step two step improved methods, bootstrapped by
#'   Nystrom4VelocityIndependent), and the adaptive embedded pairs `DPRKN4()`,
#'   `DPRKN5()`, `DPRKN6()` (with a free sixth order interpolant),
#'   `DPRKN6FM()`, `DPRKN8()`, `DPRKN12()` (Dormand, El-Mikkawy and Prince),
#'   `ERKN4()`, `ERKN5()` and `ERKN7()` (Demba, Senu and Ismail; ERKN5 controls
#'   the error of the positions only).
#'
#' `sol$stats$nf` counts evaluations of `f1` and `sol$stats$nf2` those of
#' `f2`, with Julia's conventions. Upstream IRKN3 and IRKN4 build some stages
#' from the wrong accelerations or times and fall to second order on problems
#' with explicit time dependence; here they keep orders 3 and 4 with the same
#' numbers of evaluations.
#'
#' @return An algorithm object.
#' @examples
#' prob <- SecondOrderODEProblem(function(v, u, p, t) -u, 0, 1, c(0, 1))
#' sol <- solve(prob, DPRKN6(), reltol = 1e-8, abstol = 1e-8)
#' sol$u[[length(sol$t)]] - c(-sin(1), cos(1))
#' solve(prob, Nystrom4(), dt = 0.1)$stats$nf
#' @name rkn_methods
NULL

# The velocity sentinel of the independence probe: any arithmetic, math or
# indexing on it stops with an RKNVelocityDependenceError.
.rkn_velocity_error <- function(...) {
  msg <- paste0("This RKN/Nystrom method requires a second order ODE of the form u'' = f(u, p, t) ",
                "that does not depend on the first derivative u' (the velocity). The provided ",
                "acceleration function `f1` accesses u', so this method would silently integrate ",
                "the problem at order 1. Use a method that supports velocity dependence instead, ",
                "for example Nystrom4, FineRKN4, FineRKN5 or RKN4, or reformulate the problem as ",
                "a first order ODE.")
  stop(structure(class = c("RKNVelocityDependenceError", "error", "condition"),
                 list(message = msg, call = NULL)))
}

#' @export
Ops.RKNFirstDerivativeSentinel <- function(e1, e2) .rkn_velocity_error()
#' @export
Math.RKNFirstDerivativeSentinel <- function(x, ...) .rkn_velocity_error()
#' @export
Summary.RKNFirstDerivativeSentinel <- function(..., na.rm = FALSE) .rkn_velocity_error()
#' @export
`[.RKNFirstDerivativeSentinel` <- function(x, ...) .rkn_velocity_error()
#' @export
`[[.RKNFirstDerivativeSentinel` <- function(x, ...) .rkn_velocity_error()
#' @export
as.double.RKNFirstDerivativeSentinel <- function(x, ...) .rkn_velocity_error()

# _probe_rkn_velocity_independence: f1 with the sentinel as the velocity.
.rkn_probe <- function(ctx, t, u) {
  ctx$f1(structure(list(), class = "RKNFirstDerivativeSentinel"), u, ctx$p, t)
  invisible(NULL)
}

# initialize!: f1 and f2 at the start (nf 1, nf2 1), after the probe for
# the velocity independent methods.
.rkn_init <- function(probe) function(ctx, t, z) {
  d <- .dyn_parts(ctx, z)
  if (probe) .rkn_probe(ctx, t, d$u)
  .count_nf2(ctx, 1L)
  list(fsal = c(ctx$f1(d$v, d$u, ctx$p, t), ctx$f2(d$v, d$u, ctx$p, t)), nf0 = 1L)
}

.part_tol <- function(tol, idx) if (length(tol) > 1L) tol[idx] else tol

# The error norm of an ArrayPartition((v, u)): the sum over the velocities
# plus the sum over the positions (DiffEqBase UNITLESS_ABS2).
.rkn_eest <- function(ctx, i, ev, eu, zprev, z) {
  sv <- ev / (.part_tol(ctx$abstol, i$iv) + pmax(abs(zprev[i$iv]), abs(z[i$iv])) * .part_tol(ctx$reltol, i$iv))
  su <- eu / (.part_tol(ctx$abstol, i$iu) + pmax(abs(zprev[i$iu]), abs(z[i$iu])) * .part_tol(ctx$reltol, i$iu))
  sqrt((sum(sv * sv) + sum(su * su)) / (length(sv) + length(su)))
}

# perform_step! of NystromVIConstantCache and NystromVDConstantCache.
.rkn_tab_step <- function(tab, dense6 = FALSE) {
  a <- tab$a; b <- tab$b; bp <- tab$bp; cc <- tab$c
  btilde <- tab$btilde; bptilde <- tab$bptilde
  vd <- identical(tab$type, "VD")
  abar <- tab$abar
  s <- length(b)
  function(ctx, t, z, dt, cache) {
    i <- .dyn_idx(ctx, z)
    f1 <- ctx$f1; p <- ctx$p
    duprev <- z[i$iv]; uprev <- z[i$iu]
    dtsq <- dt^2
    ks <- vector("list", s)
    ks[[1L]] <- cache$fsal[i$iv]
    for (st in 2:s) {
      acc <- a[st, 1L] * ks[[1L]]
      if (st > 2L) for (j in 2:(st - 1L)) acc <- acc + a[st, j] * ks[[j]]
      ku <- uprev + dt * (cc[st - 1L] * duprev + dt * acc)
      kdu <- duprev
      if (vd) {
        if (st == 2L) kdu <- duprev + dt * abar[2L, 1L] * ks[[1L]] else {
          accb <- abar[st, 1L] * ks[[1L]]
          for (j in 2:(st - 1L)) accb <- accb + abar[st, j] * ks[[j]]
          kdu <- duprev + dt * accb
        }
      }
      ks[[st]] <- f1(kdu, ku, p, t + dt * cc[st - 1L])
    }
    u <- uprev + dt * duprev
    du <- duprev
    for (st in seq_len(s)) {
      u <- u + dtsq * b[st] * ks[[st]]
      du <- du + dt * bp[st] * ks[[st]]
    }
    fsal <- c(f1(du, u, p, t + dt), ctx$f2(du, u, p, t + dt))
    .count_nf2(ctx, 1L)
    znew <- c(du, u)
    EEst <- NULL
    if (ctx$adaptive && length(btilde)) {
      uhat <- 0 * uprev
      duhat <- 0 * duprev
      for (st in seq_len(s)) {
        if (!vd || btilde[st] != 0) uhat <- uhat + dtsq * btilde[st] * ks[[st]]
        if (length(bptilde) && (!vd || bptilde[st] != 0)) duhat <- duhat + dt * bptilde[st] * ks[[st]]
      }
      EEst <- if (isTRUE(tab$pos_only_error)) {
        .err_norm(uhat, uprev, u, .part_tol(ctx$abstol, i$iu), .part_tol(ctx$reltol, i$iu))
      } else .rkn_eest(ctx, i, duhat, uhat, z, znew)
    }
    k <- if (dense6) ks else list(cache$fsal, fsal)
    list(u = znew, EEst = EEst, k = k, cache = list(fsal = fsal),
         nf = if (vd) tab$nf_per_step else s)
  }
}

# The free sixth order interpolant of DPRKN6 (interpolants.jl, dprkn6pre0):
# k holds the six stage accelerations.
.dprkn6_interp <- function(R, Rp) function(theta, dt, y0, y1, k) {
  poly <- function(M, r) M[r, 1L] + theta * (M[r, 2L] + theta * (M[r, 3L] + theta * (M[r, 4L] + theta * M[r, 5L])))
  nv <- length(k[[1L]])
  duprev <- y0[seq_len(nv)]; uprev <- y0[-seq_len(nv)]
  kv <- function(M) poly(M, 1L) * k[[1L]] + poly(M, 3L) * k[[3L]] + poly(M, 4L) * k[[4L]] +
    poly(M, 5L) * k[[5L]] + poly(M, 6L) * k[[6L]]
  c(duprev + dt * theta * kv(Rp), uprev + dt * theta * (duprev + dt * theta * kv(R)))
}

.rkn_alg <- function(name, order, adaptive, step, init, interp = .hermite_interp, ...)
  .new_alg(name, order, adaptive = adaptive, step = step, init = init, interp = interp,
           partitioned = TRUE, on_event = .part_on_event, ...)

.rkn_tab_alg <- function(name, order, adaptive) {
  tab <- .RKN_TABLEAUS[[name]]
  dense6 <- !is.null(tab$R)
  .rkn_alg(name, order, adaptive, step = .rkn_tab_step(tab, dense6),
           init = .rkn_init(identical(tab$type, "VI")),
           interp = if (dense6) .dprkn6_interp(tab$R, tab$Rp) else .hermite_interp)
}

#' @rdname rkn_methods
#' @export
Nystrom4 <- function() .rkn_tab_alg("Nystrom4", 4L, FALSE)
#' @rdname rkn_methods
#' @export
FineRKN4 <- function() .rkn_tab_alg("FineRKN4", 4L, TRUE)
#' @rdname rkn_methods
#' @export
FineRKN5 <- function() .rkn_tab_alg("FineRKN5", 5L, TRUE)
#' @rdname rkn_methods
#' @export
RKN4 <- function() .rkn_tab_alg("RKN4", 4L, FALSE)
#' @rdname rkn_methods
#' @export
Nystrom4VelocityIndependent <- function() .rkn_tab_alg("Nystrom4VelocityIndependent", 4L, FALSE)
#' @rdname rkn_methods
#' @export
Nystrom5VelocityIndependent <- function() .rkn_tab_alg("Nystrom5VelocityIndependent", 5L, FALSE)
#' @rdname rkn_methods
#' @export
DPRKN4 <- function() .rkn_tab_alg("DPRKN4", 4L, TRUE)
#' @rdname rkn_methods
#' @export
DPRKN5 <- function() .rkn_tab_alg("DPRKN5", 5L, TRUE)
#' @rdname rkn_methods
#' @export
DPRKN6 <- function() .rkn_tab_alg("DPRKN6", 6L, TRUE)
#' @rdname rkn_methods
#' @export
DPRKN6FM <- function() .rkn_tab_alg("DPRKN6FM", 6L, TRUE)
#' @rdname rkn_methods
#' @export
DPRKN8 <- function() .rkn_tab_alg("DPRKN8", 8L, TRUE)
#' @rdname rkn_methods
#' @export
DPRKN12 <- function() .rkn_tab_alg("DPRKN12", 12L, TRUE)
#' @rdname rkn_methods
#' @export
ERKN4 <- function() .rkn_tab_alg("ERKN4", 4L, TRUE)
#' @rdname rkn_methods
#' @export
ERKN5 <- function() .rkn_tab_alg("ERKN5", 5L, TRUE)
#' @rdname rkn_methods
#' @export
ERKN7 <- function() .rkn_tab_alg("ERKN7", 7L, TRUE)

# perform_step! of Nystrom4VelocityIndependentConstantCache, the first step
# of IRKN3 and IRKN4: the velocity, position, end derivatives and counts.
.nystrom4vi_step <- function(ctx, t, duprev, uprev, k1, dt) {
  f1 <- ctx$f1; p <- ctx$p
  halfdt <- dt / 2
  dtsq <- dt^2
  ku <- uprev + halfdt * duprev + (dtsq / 8) * k1
  k2 <- f1(duprev, ku, p, t + halfdt)
  ku <- uprev + dt * duprev + (dtsq / 2) * k2
  k3 <- f1(duprev, ku, p, t + dt)
  u <- uprev + (dtsq / 6) * (k1 + 2 * k2) + dt * duprev
  du <- duprev + (dt / 6) * (k1 + k3 + 4 * k2)
  .count_nf2(ctx, 1L)
  list(du = du, u = u, fsal = c(f1(du, u, p, t + dt), ctx$f2(du, u, p, t + dt)))
}

# The two step improved RKN methods (perform_step! of IRKN3Cache and
# IRKN4ConstantCache). The cache carries the accelerations at the last two
# points (k1c, k1c2) and the stages of the last step (k2o, k3o). The first
# step, a step after a derivative discontinuity and a step after a restart
# of the cache (an event) take a Nystrom4VelocityIndependent step and
# compute the stages of that step from its start point.
.irkn_alg <- function(name, order, four) {
  tab <- .RKN_TABLEAUS[[name]]
  .rkn_alg(name, order, FALSE, init = .rkn_init(TRUE),
    step = function(ctx, t, z, dt, cache) {
      i <- .dyn_idx(ctx, z)
      f1 <- ctx$f1; p <- ctx$p
      duprev <- z[i$iv]; uprev <- z[i$iu]
      c1 <- tab$c1; a21 <- tab$a21
      if (ctx$iter < 2L || isTRUE(ctx$derivative_discontinuity) || is.null(cache$k1c)) {
        st <- .nystrom4vi_step(ctx, t, duprev, uprev, cache$fsal[i$iv], dt)
        # the acceleration at the start (upstream evaluates it at t + c1 dt)
        k1c2 <- f1(duprev, uprev, p, t)
        kdu <- uprev + dt * (c1 * duprev + dt * a21 * k1c2)
        k2o <- f1(duprev, kdu, p, t + c1 * dt)
        k3o <- NULL
        if (four) {
          # the third stage from the second (upstream uses the start
          # acceleration and the node c1 here)
          kdu <- uprev + dt * (tab$c2 * duprev + dt * tab$a32 * k2o)
          k3o <- f1(duprev, kdu, p, t + tab$c2 * dt)
        }
        return(list(u = c(st$du, st$u), err = NULL, k = list(cache$fsal, st$fsal),
                    cache = list(fsal = st$fsal, k1c = st$fsal[i$iv], k1c2 = k1c2, k2o = k2o, k3o = k3o),
                    nf = if (four) 6L else 5L))
      }
      uprev2 <- ctx$uprev2
      duprev2 <- uprev2[i$iv]
      k1c <- cache$k1c; k1c2 <- cache$k1c2
      # the second stage from the current point (upstream IRKN3 uses the
      # acceleration of the previous point here)
      ku <- uprev + dt * (c1 * duprev + dt * a21 * k1c)
      k2 <- f1(duprev, ku, p, t + c1 * dt)
      d2 <- k2 - cache$k2o
      if (four) {
        ku <- uprev + dt * (tab$c2 * duprev + dt * tab$a32 * k2)
        k3 <- f1(duprev, ku, p, t + tab$c2 * dt)
        d3 <- k3 - cache$k3o
        du <- duprev + dt * (tab$b1 * k1c + tab$bbar1 * k1c2 + tab$b2 * d2 + tab$b3 * d3)
        u <- uprev + dt * tab$bconst1 * duprev +
          dt * (tab$bconst2 * duprev2 + dt * (tab$bbar2 * d2 + tab$bbar3 * d3))
      } else {
        k3 <- NULL
        u <- uprev + tab$bconst1 * dt * duprev + dt * (tab$bconst2 * duprev2 + dt * tab$bbar2 * d2)
        du <- duprev + dt * (tab$b1 * k1c + tab$bbar1 * k1c2 + tab$b2 * d2)
      }
      fsal <- c(f1(du, u, p, t + dt), ctx$f2(du, u, p, t + dt))
      .count_nf2(ctx, 1L)
      list(u = c(du, u), err = NULL, k = list(cache$fsal, fsal),
           cache = list(fsal = fsal, k1c = fsal[i$iv], k1c2 = k1c, k2o = k2, k3o = k3),
           nf = if (four) 4L else 3L)
    })
}

#' @rdname rkn_methods
#' @export
IRKN3 <- function() .irkn_alg("IRKN3", 3L, FALSE)
#' @rdname rkn_methods
#' @export
IRKN4 <- function() .irkn_alg("IRKN4", 4L, TRUE)
