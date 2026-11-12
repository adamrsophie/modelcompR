# Symplectic integrators of OrdinaryDiffEqSymplecticRK (algorithms.jl,
# alg_utils.jl, symplectic_caches.jl, symplectic_perform_step.jl and the
# coefficients of symplectic_tableaus.jl in R/solvers-rkn-tableaus.R), for
# DynamicalODEProblem and SecondOrderODEProblem (R/scimlbase-problems.R).
#
# The state is c(v, u), velocities first as Julia's ArrayPartition((v, u)),
# with dv/dt = f1(v, u, p, t) and du/dt = f2(v, u, p, t). As the out of
# place caches of Julia, a step keeps the derivatives c(kdu, ku) it ends
# with (fsallast) for the next step and for the cubic Hermite dense output
# between fsalfirst and fsallast. Julia counts the evaluations of f1 in nf
# and those of f2 in nf2; so does sol$stats.

#' Symplectic integrators
#'
#' Fixed step symplectic methods of OrdinaryDiffEqSymplecticRK for
#' [DynamicalODEProblem()] and [SecondOrderODEProblem()] (pass `dt`). The
#' state is `c(v, u)` with `dv/dt = f1(v, u, p, t)` and `du/dt = f2(v, u, p,
#' t)`; separable Hamiltonian systems have `f1 = -dV/du` and `f2 = dT/dv`.
#'
#' * `SymplecticEuler()`: first order, semi implicit Euler.
#' * `VelocityVerlet()`: second order; needs `f2(v, u, p, t) = v` (a second
#'   order problem) and stops otherwise.
#' * `VerletLeapfrog()`: second order kick drift kick leapfrog, one `f1`
#'   evaluation per step; `LeapfrogDriftKickDrift()` is its drift kick
#'   drift form for `f1` depending on `v` (two evaluations). Both default to
#'   linear interpolation of the saved steps (`dense = FALSE`).
#' * `PseudoVerletLeapfrog()` and `McAte2()` (second order), `Ruth3()` and
#'   `McAte3()` (third), `CandyRoz4()`, `McAte4()`, `CalvoSanz4()` and
#'   `McAte42()` (fourth; McLachlan's method upstream marks as broken),
#'   `McAte5()` (fifth), `Yoshida6()` and `KahanLi6()` (sixth), `McAte8()`
#'   and `KahanLi8()` (eighth) and `SofSpa10()` (tenth order) are
#'   compositions of drifts and kicks with the coefficients of Julia.
#'
#' Dense output is cubic Hermite interpolation from the derivatives the
#' steps keep at both ends. `sol$stats$nf` counts evaluations of `f1` and
#' `sol$stats$nf2` those of `f2`, as Julia does. After an event that changes
#' the state, every method steps on from the derivatives of the new state
#' (upstream `VelocityVerlet()`, `VerletLeapfrog()` and the compositions
#' reuse those from before the event).
#'
#' @return An algorithm object.
#' @examples
#' prob <- SecondOrderODEProblem(function(v, u, p, t) -u, 0, 1, c(0, 1))
#' sol <- solve(prob, McAte5(), dt = 0.1)
#' sol$u[[length(sol$t)]] - c(-sin(1), cos(1))
#' solve(prob, VelocityVerlet(), dt = 0.1)$stats$nf
#' @name symplectic_methods
NULL

# The parts of a partitioned problem, or the error of Julia's solve for a
# partitioned method on another problem (OrdinaryDiffEqCore solve.jl).
.dyn_parts <- function(ctx, z) {
  if (is.null(ctx$nv) || is.null(ctx$f1) || is.null(ctx$f2))
    stop("You can not use a solver designed for partitioned ODE with this problem. ",
         "Please choose a solver suitable for your problem", call. = FALSE)
  iv <- seq_len(ctx$nv)
  list(iv = iv, iu = ctx$nv + seq_len(length(z) - ctx$nv), v = z[iv], u = z[-iv])
}

.dyn_idx <- function(ctx, z) {
  iv <- seq_len(ctx$nv)
  list(iv = iv, iu = ctx$nv + seq_len(length(z) - ctx$nv))
}

# stats.nf2: the f2 evaluations, kept in the counter environment of ctx
.count_nf2 <- function(ctx, n) {
  e <- ctx$counts
  e$nf2 <- (e$nf2 %||% 0L) + as.integer(n)
  invisible(NULL)
}

# After an event changed the state: f at the new state, counted once as
# Julia's reset_fsal! does. Upstream refreshes only fsalfirst, so the steps
# that read fsallast (VelocityVerlet, VerletLeapfrog, the compositions) start
# from the derivatives before the event; here every method restarts from the
# new state (and IRKN3/IRKN4 bootstrap again).
.part_on_event <- function(ctx, t, z, cache) list(cache = list(fsal = ctx$f(z, ctx$p, t)), nf = 1L)

.symp_alg <- function(name, order, step, init = .symp_init, linear = FALSE)
  .new_alg(name, order, adaptive = FALSE, step = step, init = init, interp = .hermite_interp,
           partitioned = TRUE, linear_interpolation = linear, on_event = .part_on_event)

# initialize! of the Hamilton caches: f1 and f2 at the start.
.symp_init <- function(ctx, t, z) {
  d <- .dyn_parts(ctx, z)
  .count_nf2(ctx, 1L)
  list(fsal = c(ctx$f1(d$v, d$u, ctx$p, t), ctx$f2(d$v, d$u, ctx$p, t)), nf0 = 1L)
}

# The end of a step: the new state, the Hermite data and the counts.
.symp_out <- function(ctx, du, u, kdu, ku, cache, nf, nf2) {
  .count_nf2(ctx, nf2)
  fsal <- c(kdu, ku)
  list(u = c(du, u), err = NULL, k = list(cache$fsal, fsal), cache = list(fsal = fsal), nf = nf)
}

#' @rdname symplectic_methods
#' @export
SymplecticEuler <- function() {
  .symp_alg("SymplecticEuler", 1L,
    init = function(ctx, t, z) {
      # fsalfirst = (f1, f2) at the start; Julia counts two f1 and one f2 call
      d <- .dyn_parts(ctx, z)
      .count_nf2(ctx, 1L)
      list(fsal = c(ctx$f1(d$v, d$u, ctx$p, t), ctx$f2(d$v, d$u, ctx$p, t)), nf0 = 2L)
    },
    step = function(ctx, t, z, dt, cache) {
      i <- .dyn_idx(ctx, z)
      duprev <- z[i$iv]; uprev <- z[i$iu]
      u <- uprev + dt * cache$fsal[i$iu]
      kdu <- ctx$f1(duprev, u, ctx$p, t)
      du <- duprev + dt * kdu
      ku <- ctx$f2(du, u, ctx$p, t)
      .symp_out(ctx, du, u, kdu, ku, cache, 1L, 1L)
    })
}

#' @rdname symplectic_methods
#' @export
VelocityVerlet <- function() {
  .symp_alg("VelocityVerlet", 2L,
    init = function(ctx, t, z) {
      d <- .dyn_parts(ctx, z)
      ku <- ctx$f2(d$v, d$u, ctx$p, t)
      # verify_f2: the method needs du/dt = v
      if (length(ku) != length(d$v) || !isTRUE(all(ku == d$v)))
        stop("Algorithm VelocityVerlet is not applicable if f2(p, q, t) != p", call. = FALSE)
      .count_nf2(ctx, 1L)
      list(fsal = c(ctx$f1(d$v, d$u, ctx$p, t), ku), nf0 = 1L)
    },
    step = function(ctx, t, z, dt, cache) {
      i <- .dyn_idx(ctx, z)
      duprev <- z[i$iv]; uprev <- z[i$iu]
      # x(t + dt) = x(t) + v(t) dt + a(t) dt^2 / 2
      a0 <- cache$fsal[i$iv]
      u <- uprev + dt * duprev + dt^2 * (0.5 * a0)
      kdu <- ctx$f1(duprev, u, ctx$p, t + dt)
      # v(t + dt) = v(t) + (a(t) + a(t + dt)) dt / 2
      du <- duprev + dt * (0.5 * a0 + 0.5 * kdu)
      .symp_out(ctx, du, u, kdu, du, cache, 1L, 0L)
    })
}

#' @rdname symplectic_methods
#' @export
VerletLeapfrog <- function() {
  .symp_alg("VerletLeapfrog", 2L, linear = TRUE,
    step = function(ctx, t, z, dt, cache) {
      i <- .dyn_idx(ctx, z)
      duprev <- z[i$iv]; uprev <- z[i$iu]
      # kick, drift, kick
      du <- duprev + dt * 0.5 * cache$fsal[i$iv]
      ku <- ctx$f2(du, uprev, ctx$p, t + 0.5 * dt)
      u <- uprev + dt * ku
      kdu <- ctx$f1(du, u, ctx$p, t + dt)
      du <- du + dt * 0.5 * kdu
      .symp_out(ctx, du, u, kdu, ku, cache, 1L, 1L)
    })
}

#' @rdname symplectic_methods
#' @export
LeapfrogDriftKickDrift <- function() {
  .symp_alg("LeapfrogDriftKickDrift", 2L, linear = TRUE,
    step = function(ctx, t, z, dt, cache) {
      i <- .dyn_idx(ctx, z)
      duprev <- z[i$iv]; uprev <- z[i$iu]
      f1 <- ctx$f1; f2 <- ctx$f2; p <- ctx$p
      # half drift and half kick, then the full kick and the second half drift
      ku <- f2(duprev, uprev, p, t)
      u <- uprev + dt * 0.5 * ku
      kdu <- f1(duprev, uprev, p, t)
      du <- duprev + dt * 0.5 * kdu
      kdu <- f1(du, u, p, t + 0.5 * dt)
      du <- duprev + dt * kdu
      ku <- f2(du, u, p, t + dt)
      u <- u + dt * 0.5 * ku
      .symp_out(ctx, du, u, kdu, ku, cache, 2L, 2L)
    })
}

# The generic drift kick step of the composition methods (perform_step! of
# SymplecticTableau): drifts with b, kicks with a; the clock advances with
# a, as upstream.
.symp_tableau_alg <- function(name, order) {
  tab <- .SYMP_TABLEAUS[[name]]
  a <- tab$a; b <- tab$b; N <- length(a)
  .symp_alg(name, order,
    step = function(ctx, t, z, dt, cache) {
      i <- .dyn_idx(ctx, z)
      f1 <- ctx$f1; f2 <- ctx$f2; p <- ctx$p
      duprev <- z[i$iv]; uprev <- z[i$iu]
      u <- uprev + dt * b[1L] * cache$fsal[i$iu]
      kdu <- f1(duprev, u, p, t)
      du <- duprev + dt * a[1L] * kdu
      tnew <- t + a[1L] * dt
      for (j in 2:N) {
        ku <- f2(du, u, p, tnew)
        u <- u + dt * b[j] * ku
        kdu <- f1(du, u, p, tnew)
        du <- du + dt * a[j] * kdu
        if (j < N) tnew <- tnew + a[j] * dt
      }
      # a last kick changed du: a fresh f1 for the next step
      if (a[N] != 0) kdu <- f1(du, u, p, tnew)
      ku <- f2(du, u, p, tnew)
      .symp_out(ctx, du, u, kdu, ku, cache, N + as.integer(a[N] != 0), N)
    })
}

#' @rdname symplectic_methods
#' @export
PseudoVerletLeapfrog <- function() .symp_tableau_alg("PseudoVerletLeapfrog", 2L)
#' @rdname symplectic_methods
#' @export
McAte2 <- function() .symp_tableau_alg("McAte2", 2L)
#' @rdname symplectic_methods
#' @export
Ruth3 <- function() .symp_tableau_alg("Ruth3", 3L)
#' @rdname symplectic_methods
#' @export
McAte3 <- function() .symp_tableau_alg("McAte3", 3L)
#' @rdname symplectic_methods
#' @export
CandyRoz4 <- function() .symp_tableau_alg("CandyRoz4", 4L)
#' @rdname symplectic_methods
#' @export
McAte4 <- function() .symp_tableau_alg("McAte4", 4L)
#' @rdname symplectic_methods
#' @export
CalvoSanz4 <- function() .symp_tableau_alg("CalvoSanz4", 4L)
#' @rdname symplectic_methods
#' @export
McAte42 <- function() .symp_tableau_alg("McAte42", 4L)
#' @rdname symplectic_methods
#' @export
McAte5 <- function() .symp_tableau_alg("McAte5", 5L)
#' @rdname symplectic_methods
#' @export
Yoshida6 <- function() .symp_tableau_alg("Yoshida6", 6L)
#' @rdname symplectic_methods
#' @export
KahanLi6 <- function() .symp_tableau_alg("KahanLi6", 6L)
#' @rdname symplectic_methods
#' @export
McAte8 <- function() .symp_tableau_alg("McAte8", 8L)
#' @rdname symplectic_methods
#' @export
KahanLi8 <- function() .symp_tableau_alg("KahanLi8", 8L)
#' @rdname symplectic_methods
#' @export
SofSpa10 <- function() .symp_tableau_alg("SofSpa10", 10L)
