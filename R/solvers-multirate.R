# Multirate methods for split problems du/dt = f1(u, p, t) + f2(u, p, t)
# with f1 the fast part and f2 the slow part, after OrdinaryDiffEqMultirate
# (algorithms.jl, alg_utils.jl, multirate_tableaus.jl,
# multirate_perform_step.jl: MREEF, MRAB, MRIGARKERK22a, MRIGARKERK22b,
# MRIGARKERK33a, MRIGARKERK45a, MIS) and OrdinaryDiffEqMultirateImplicit
# 1.0.1 (algorithms.jl, alg_utils.jl, multirate_implicit_tableaus.jl,
# multirate_implicit_caches.jl, multirate_implicit_perform_step.jl: MREIL,
# MRIGARKIRK21a, MRIGARKESDIRK34a, MRIGARKESDIRK46a). The steps follow the
# out of place perform_step! methods: those are @muladd, so their broadcast
# multiply adds and the scalar times t + c dt are fused (.fma), while the
# inner rate and substage helpers of the MRI-GARK methods are not.

#' Multirate methods
#'
#' Methods for a [SplitODEProblem()] `du/dt = f1(u, p, t) + f2(u, p, t)`
#' whose first part `f1` is fast (or stiff) and whose second part `f2` is
#' slow, after OrdinaryDiffEqMultirate and OrdinaryDiffEqMultirateImplicit.
#' The fast part is advanced with small substeps inside each step while the
#' slow part is evaluated only a few times per step. `sol$stats$nf` counts
#' the evaluations of `f1` (and those of `f2` inside the Newton iterations
#' of the implicit MRI-GARK methods), `sol$stats$nf2` the other evaluations
#' of `f2`.
#'
#' * `MREEF()`: multirate Euler with the slow rate frozen over each macro
#'   interval and `m` explicit Euler substeps of the fast part, raised to
#'   order `order` by Aitken-Neville (Richardson) extrapolation over the
#'   step counts `seq` (`"harmonic"`: 1, 2, 3, ...; `"romberg"`: 1, 2, 4,
#'   ...).
#' * `MRAB()`: multirate Adams-Bashforth: the slow rate is evaluated once per
#'   step and the `m` fast substeps use a `k` step Adams-Bashforth formula on
#'   the combined rate, started with lower orders inside each step. The error
#'   estimate compares orders `k` and `k - 1` (with `k = 1` there is none:
#'   the error estimate stays at 1, as in OrdinaryDiffEq).
#' * `MRIGARKERK22a()`, `MRIGARKERK22b()`, `MRIGARKERK33a()`,
#'   `MRIGARKERK45a()`: the explicit multirate infinitesimal GARK methods of
#'   Sandu (orders 2, 2, 3 and 4). Each stage integrates a modified fast
#'   problem with `m` explicit Runge-Kutta micro steps (midpoint, RK3 or RK4)
#'   and a slow forcing that is polynomial in time.
#' * `MIS()`: the multirate infinitesimal step method MIS2 of Wensch, Knoth
#'   and Galant (order 2), whose inner problems take `max(1, ceiling(m d_i))`
#'   explicit midpoint micro steps.
#' * `MREIL()`: as `MREEF()` with linearly implicit Euler substeps of the
#'   fast part, `(I - h J) delta = h (f1 + f2)` with `J` the Jacobian of `f1`
#'   (from `SplitFunction(f1, f2, jac = J1)`, else finite differences), for a
#'   stiff fast part. One Jacobian per step, one factorization per
#'   extrapolation column. A constant mass matrix enters as `M - h J`.
#' * `MRIGARKIRK21a()`, `MRIGARKESDIRK34a()`, `MRIGARKESDIRK46a()`: the solve
#'   decoupled implicit MRI-GARK methods of Sandu (orders 2, 3 and 4) for a
#'   stiff slow part: the fast stages are explicit as above, the slow
#'   corrections are diagonally implicit stages solved by Newton's method on
#'   `f2` (finite difference Jacobian of `f2`).
#'
#' @param m Number of fast substeps (or inner micro steps) per step or stage.
#' @param order Extrapolation order: the number of base solutions.
#' @param seq Step count sequence of the extrapolation, `"harmonic"` or
#'   `"romberg"`.
#' @param k Adams-Bashforth order, from 1 to 5.
#' @param autodiff As for [ImplicitEuler()]: how the Jacobian is formed when
#'   none is given.
#' @return An ODE algorithm.
#' @examples
#' # a fast linear relaxation coupled to a slow nonlinear oscillator
#' fast <- function(u, p, t) c(-20 * (u[1] - u[2]), 0)
#' slow <- function(u, p, t) c(0, -u[1] - 0.1 * u[2] * u[2] * u[2])
#' prob <- SplitODEProblem(fast, slow, c(1, 0.5), c(0, 1))
#' sol <- solve(prob, MRIGARKERK33a(m = 3), reltol = 1e-6, abstol = 1e-8)
#' sol$u[[length(sol$u)]]
#' solve(prob, MREEF(), dt = 1 / 16, adaptive = FALSE)$u[[17]]
#' @name multirate_methods
NULL

# ---------------------------------------------------------------------------
# Shared parts

.mr_check_split <- function(ctx, name) {
  if (is.null(ctx$f1) || is.null(ctx$f2))
    stop(name, "() integrates a SplitODEProblem (f1 fast, f2 slow)", call. = FALSE)
}

.mr_check_m <- function(m, name) {
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m != round(m) || m < 1)
    stop(name, ": `m` must be >= 1", call. = FALSE)
  as.integer(m)
}

# initialize!: fsalfirst = f1(uprev) + f2(uprev), counted in nf and nf2.
.mr_init_fn <- function(name) function(ctx, t, u) {
  .mr_check_split(ctx, name)
  list(fsal = ctx$f1(u, ctx$p, t) + ctx$f2(u, ctx$p, t), nf0 = 1L, nf20 = 1L)
}

# The methods are not FSAL: an event changes nothing they keep.
.mr_on_event <- function(ctx, t, u, cache) list(cache = cache, nf = 0L)

# _extrapolation_sequence: the step counts of the extrapolation columns.
.mr_ext_seq <- function(seq, order) {
  if (identical(seq, "harmonic")) return(seq_len(order))
  if (identical(seq, "romberg")) return(2L^(seq_len(order) - 1L))
  stop("unknown sequence `", seq, "`, choose \"harmonic\" or \"romberg\"", call. = FALSE)
}

.mr_check_seq <- function(seq) {
  seq <- as.character(seq)
  if (length(seq) != 1L || !seq %in% c("harmonic", "romberg"))
    stop("unknown sequence `", paste(seq, collapse = " "), "`, choose \"harmonic\" or \"romberg\"",
         call. = FALSE)
  seq
}

# ---------------------------------------------------------------------------
# MREEF

#' @rdname multirate_methods
#' @export
MREEF <- function(m = 4L, order = 4L, seq = "harmonic") {
  m <- as.integer(m)
  order <- as.integer(order)
  seq <- .mr_check_seq(seq)
  .new_alg("MREEF", order,
           step = function(ctx, t, u, dt, cache) .mr_mreef_step(ctx, t, u, dt, cache, m, order, seq),
           init = .mr_init_fn("MREEF"), on_event = .mr_on_event, m = m, extrapolation_order = order,
           seq = seq)
}

.mr_mreef_step <- function(ctx, t, u, dt, cache, m, order, seq) {
  f1 <- ctx$f1; f2 <- ctx$f2; p <- ctx$p
  ns <- .mr_ext_seq(seq, order)
  if (ctx$adaptive && order < 2L)
    stop("MREEF: an adaptive solve needs `order` >= 2 (the error estimate compares two columns)",
         call. = FALSE)
  nf <- 0L; nf2 <- 0L
  Tc <- vector("list", order)
  fsal0 <- NULL
  for (j in seq_len(order)) {
    nj <- ns[j]
    h_mac <- dt / nj
    h_fast <- h_mac / m
    uc <- u
    for (i_mac in seq_len(nj)) {
      t_mac <- .fma(i_mac - 1, h_mac, t)
      ks <- f2(uc, p, t_mac)
      nf2 <- nf2 + 1L
      for (i_fast in seq_len(m)) {
        t_fast <- .fma(i_fast - 1, h_fast, t_mac)
        kf <- f1(uc, p, t_fast)
        nf <- nf + 1L
        if (is.null(fsal0)) fsal0 <- ks + kf
        uc <- .fma(h_fast, kf, .fma(h_fast, ks, uc))
      }
    }
    Tc[[j]] <- uc
  }
  # Aitken-Neville extrapolation, rows in reverse order
  for (k in seq_len(order - 1L)) {
    for (j in order:(k + 1L)) {
      ratio <- ns[j] / ns[j - k]
      Tc[[j]] <- Tc[[j]] + (Tc[[j]] - Tc[[j - 1L]]) / (ratio - 1)
    }
  }
  unew <- Tc[[order]]
  EEst <- if (ctx$adaptive) .err_norm(Tc[[order]] - Tc[[order - 1L]], u, unew, ctx$abstol, ctx$reltol) else 0
  fsal1 <- f1(unew, p, t + dt) + f2(unew, p, t + dt)
  list(u = unew, EEst = EEst, k = list(fsal0, fsal1), cache = list(),
       nf = nf + 1L, nf2 = nf2 + 1L)
}

# ---------------------------------------------------------------------------
# MRAB

# MRABTableau: the Adams-Bashforth weights of orders 1 to 5, newest rate
# first.
.mr_ab_weights <- function() list(
  1,
  c(3 / 2, -1 / 2),
  c(23 / 12, -16 / 12, 5 / 12),
  c(55 / 24, -59 / 24, 37 / 24, -9 / 24),
  c(1901 / 720, -2774 / 720, 2616 / 720, -1274 / 720, 251 / 720)
)

#' @rdname multirate_methods
#' @export
MRAB <- function(k = 2L, m = 4L) {
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k != round(k) || k < 1 || k > 5)
    stop("MRAB: `k` must satisfy 1 <= k <= 5", call. = FALSE)
  k <- as.integer(k)
  m <- .mr_check_m(m, "MRAB")
  .new_alg("MRAB", k,
           step = function(ctx, t, u, dt, cache) .mr_mrab_step(ctx, t, u, dt, cache, k, m),
           init = .mr_init_fn("MRAB"), on_event = .mr_on_event, k = k, m = m)
}

.mr_mrab_step <- function(ctx, t, u, dt, cache, k, m) {
  f1 <- ctx$f1; f2 <- ctx$f2; p <- ctx$p
  beta <- .mr_ab_weights()
  h <- dt / m
  bk <- beta[[k]]
  uc <- u
  ks <- f2(uc, p, t)
  # the history restarts in every step; with m < k the oldest entries are
  # never filled (OrdinaryDiffEq leaves them undefined): they count as zero
  # in the error estimate
  hist <- rep(list(numeric(length(u))), k)
  n_hist <- 0L
  fsal0 <- NULL
  for (l in seq_len(m)) {
    t_fast <- .fma(l - 1, h, t)
    kf <- f1(uc, p, t_fast)
    Fl <- kf + ks
    if (l == 1L) fsal0 <- Fl
    n_hist <- min(n_hist + 1L, k)
    if (n_hist >= 2L) for (j in n_hist:2L) hist[[j]] <- hist[[j - 1L]]
    hist[[1L]] <- Fl
    kk <- min(l, k)
    bu <- beta[[kk]]
    du <- numeric(length(u))
    for (i in seq_len(kk)) du <- .fma(h * bu[i], hist[[i]], du)
    uc <- uc + du
  }
  EEst <- 1
  if (ctx$adaptive && k >= 2L) {
    blow <- beta[[k - 1L]]
    ut <- h * (bk[1L] - blow[1L]) * hist[[1L]]
    if (k >= 3L) for (i in 2:(k - 1L)) ut <- .fma(h * (bk[i] - blow[i]), hist[[i]], ut)
    ut <- .fma(h * bk[k], hist[[k]], ut)
    EEst <- .err_norm(ut, u, uc, ctx$abstol, ctx$reltol)
  }
  fsal1 <- f1(uc, p, t + dt) + f2(uc, p, t + dt)
  list(u = uc, EEst = EEst, k = list(fsal0, fsal1), cache = list(),
       nf = m + 1L, nf2 = 2L)
}

# ---------------------------------------------------------------------------
# MRI-GARK methods (Sandu 2019). Stage i integrates the fast problem
#   v' = dc_i dt f1(v) + dt sum_j omega_ij(tau) f2(z_j),  tau in [0, 1],
# omega_ij(tau) = W0[i, j] + W1[i, j] tau, with an explicit Runge-Kutta
# method of order q; a stage with dc_i = 0 is a pure slow quadrature (or,
# with gamma0[i] != 0, an implicit slow correction). Wemb0 and Wemb1 give
# the embedded last stage of the error estimate (NULL: the difference to
# the last stage value).

.mr_gark_tableau <- function(name) {
  z <- function(s) matrix(0, s, s)
  switch(name,
    MRIGARKERK22a = list(dc = c(1 / 2, 1 / 2), W0 = rbind(c(1 / 2, 0), c(-1 / 2, 1)), W1 = z(2L),
                         Wemb0 = c(1 / 2, 0), Wemb1 = NULL, g0 = c(0, 0), q = 2L, order = 2L),
    MRIGARKERK22b = list(dc = c(1, 0), W0 = rbind(c(1, 0), c(-1 / 2, 1 / 2)), W1 = z(2L),
                         Wemb0 = NULL, Wemb1 = NULL, g0 = c(0, 0), q = 2L, order = 2L),
    MRIGARKERK33a = list(dc = c(1 / 3, 1 / 3, 1 / 3),
                         W0 = rbind(c(1 / 3, 0, 0), c(-1 / 3, 2 / 3, 0), c(0, -2 / 3, 1)),
                         W1 = rbind(c(0, 0, 0), c(0, 0, 0), c(1 / 2, 0, -1 / 2)),
                         Wemb0 = c(1 / 12, -1 / 3, 7 / 12), Wemb1 = NULL, g0 = c(0, 0, 0), q = 3L,
                         order = 3L),
    MRIGARKERK45a = .mr_gark_erk45a(),
    MRIGARKIRK21a = list(dc = c(1, 0), W0 = rbind(c(1, 0), c(-1 / 2, 0)), W1 = z(2L),
                         Wemb0 = NULL, Wemb1 = NULL, g0 = c(0, 1 / 2), q = 2L, order = 2L),
    MRIGARKESDIRK34a = .mr_gark_esdirk34a(),
    MRIGARKESDIRK46a = .mr_gark_esdirk46a())
}

.mr_gark_erk45a <- function() {
  W0 <- matrix(0, 5L, 5L)
  W0[1, 1] <- 1 / 5
  W0[2, 1] <- -53 / 16; W0[2, 2] <- 281 / 80
  W0[3, 1] <- -36562993 / 71394880; W0[3, 2] <- 34903117 / 17848720
  W0[3, 3] <- -88770499 / 71394880
  W0[4, 1] <- -7631593 / 71394880; W0[4, 2] <- -166232021 / 35697440
  W0[4, 3] <- 6068517 / 1519040; W0[4, 4] <- 8644289 / 8924360
  W0[5, 1] <- 277061 / 303808; W0[5, 2] <- -209323 / 1139280
  W0[5, 3] <- -1360217 / 1139280; W0[5, 4] <- -148789 / 56964; W0[5, 5] <- 147889 / 45120
  W1 <- matrix(0, 5L, 5L)
  W1[2, 1] <- 503 / 80; W1[2, 2] <- -503 / 80
  W1[3, 1] <- -1365537 / 35697440; W1[3, 2] <- 4963773 / 7139488
  W1[3, 3] <- -1465833 / 2231090
  W1[4, 1] <- 66974357 / 35697440; W1[4, 2] <- 21445367 / 7139488
  W1[4, 3] <- -3; W1[4, 4] <- -8388609 / 4462180
  W1[5, 1] <- -18227 / 7520; W1[5, 2] <- 2; W1[5, 3] <- 1; W1[5, 4] <- 5
  W1[5, 5] <- -41933 / 7520
  list(dc = rep(1 / 5, 5L), W0 = W0, W1 = W1,
       Wemb0 = c(-88227 / 47470, 756870829 / 340217490, -713704111 / 1360869960,
                 -31967827 / 340217490, 129673 / 286680),
       Wemb1 = c(6213 / 1880, -6213 / 1880, 0, 0, 0), g0 = rep(0, 5L), q = 4L, order = 4L)
}

# The decimal coefficients of MRIGARKESDIRK34aTableau, as Julia parses them
.mr_gark_esdirk34a <- function() {
  b <- 0.435866521508459
  W0 <- matrix(0, 6L, 6L)
  W0[1, 1] <- 1 / 3
  W0[2, 1] <- -b
  W0[3, 1] <- -0.3045790611944505
  W0[3, 3] <- 0.6379123945277838
  W0[4, 1] <- 0.21169131056402665
  W0[4, 3] <- -0.6475578320724856
  W0[5, 1] <- 0.4454209388055495
  W0[5, 3] <- 0.8813784805616198
  W0[5, 5] <- -0.993466086033836
  W0[6, 1] <- -b
  list(dc = c(1 / 3, 0, 1 / 3, 0, 1 / 3, 0), W0 = W0, W1 = matrix(0, 6L, 6L), Wemb0 = NULL,
       Wemb1 = NULL, g0 = c(0, b, 0, b, 0, b), q = 3L, order = 3L)
}

.mr_gark_esdirk46a <- function() {
  W0 <- matrix(0, 11L, 11L)
  W0[1, 1] <- 1 / 5
  W0[2, 1] <- -1 / 4
  W0[3, 1] <- 1771023115159 / 1929363690800
  W0[3, 3] <- -1385150376999 / 1929363690800
  W0[4, 1] <- 914009 / 345800
  W0[4, 3] <- -1000459 / 345800
  W0[5, 1] <- 18386293581909 / 36657910125200
  W0[5, 3] <- 5506531089 / 80566835440
  W0[5, 5] <- -178423463189 / 482340922700
  W0[6, 1] <- 36036097 / 8299200
  W0[6, 3] <- 4621 / 118560
  W0[6, 5] <- -38434367 / 8299200
  W0[7, 1] <- -247809665162987 / 146631640500800
  W0[7, 3] <- 10604946373579 / 14663164050080
  W0[7, 5] <- 10838126175385 / 5865265620032
  W0[7, 7] <- -24966656214317 / 36657910125200
  W0[8, 1] <- 38519701 / 11618880
  W0[8, 3] <- 10517363 / 9682400
  W0[8, 5] <- -23284701 / 19364800
  W0[8, 7] <- -10018609 / 2904720
  W0[9, 1] <- -52907807977903 / 33838070884800
  W0[9, 3] <- 74846944529257 / 73315820250400
  W0[9, 5] <- 365022522318171 / 146631640500800
  W0[9, 7] <- -20513210406809 / 109973730375600
  W0[9, 9] <- -2918009798 / 1870301537
  W0[10, 1] <- 19 / 100
  W0[10, 3] <- -73 / 300
  W0[10, 5] <- 127 / 300
  W0[10, 7] <- 127 / 300
  W0[10, 9] <- -313 / 300
  W1 <- matrix(0, 11L, 11L)
  W1[3, 1] <- -1674554930619 / 964681845400
  W1[3, 3] <- 1674554930619 / 964681845400
  W1[4, 1] <- -1007739 / 172900
  W1[4, 3] <- 1007739 / 172900
  W1[5, 1] <- -8450070574289 / 18328955062600
  W1[5, 3] <- -39429409169 / 40283417720
  W1[5, 5] <- 173621393067 / 120585230675
  W1[6, 1] <- -122894383 / 16598400
  W1[6, 3] <- 14501 / 237120
  W1[6, 5] <- 121879313 / 16598400
  W1[7, 1] <- 32410002731287 / 15434909526400
  W1[7, 3] <- -46499276605921 / 29326328100160
  W1[7, 5] <- -34914135774643 / 11730531240064
  W1[7, 7] <- 45128506783177 / 18328955062600
  W1[8, 1] <- -128357303 / 23237760
  W1[8, 3] <- -35433927 / 19364800
  W1[8, 5] <- 71038479 / 38729600
  W1[8, 7] <- 8015933 / 1452360
  W1[9, 1] <- 136721604296777 / 67676141769600
  W1[9, 3] <- -349632444539303 / 146631640500800
  W1[9, 5] <- -1292744859249609 / 293263281001600
  W1[9, 7] <- 8356250416309 / 54986865187800
  W1[9, 9] <- 17282943803 / 3740603074
  W1[10, 1] <- 3 / 25
  W1[10, 3] <- -29 / 300
  W1[10, 5] <- 71 / 300
  W1[10, 7] <- 71 / 300
  W1[10, 9] <- -149 / 300
  Wemb0 <- numeric(11L)
  Wemb0[1] <- -1 / 4
  Wemb0[3] <- 5595 / 8804
  Wemb0[5] <- -2445 / 8804
  Wemb0[7] <- -4225 / 8804
  Wemb0[9] <- 2205 / 4402
  Wemb0[11] <- -567 / 4402
  list(dc = c(rep(c(1 / 5, 0), 5L), 0), W0 = W0, W1 = W1, Wemb0 = Wemb0, Wemb1 = NULL,
       g0 = c(rep(c(0, 1 / 4), 5L), 0), q = 4L, order = 4L)
}

# _mrigark_rate (not @muladd): the rate of the fast problem at tau.
.mr_gark_rate <- function(f1, w, p, t, dt, dci, cprev, tau, w0, w1, fS, upto) {
  t_phys <- t + (cprev + tau * dci) * dt
  g <- dt * dci * f1(w, p, t_phys)
  for (j in seq_len(upto)) {
    om <- if (is.null(w1)) w0[j] else w0[j] + w1[j] * tau
    if (om == 0) next
    g <- g + dt * om * fS[[j]]
  }
  g
}

# _mrigark_substage (not @muladd): one fast stage from vstart over tau in
# [0, 1] with m micro steps of the order q explicit method; returns the end
# value and the count of f1 evaluations.
.mr_gark_substage <- function(vstart, f1, p, t, dt, dci, cprev, m, q, w0, w1, fS, upto) {
  if (dci == 0) {
    v <- vstart
    for (j in seq_len(upto)) {
      om <- if (is.null(w1)) w0[j] else w0[j] + w1[j] / 2
      if (om == 0) next
      v <- v + dt * om * fS[[j]]
    }
    return(list(v = v, nf = 0L))
  }
  rate <- function(w, tau) .mr_gark_rate(f1, w, p, t, dt, dci, cprev, tau, w0, w1, fS, upto)
  h <- 1 / m
  v <- vstart
  for (ks in seq_len(m)) {
    tau0 <- (ks - 1) * h
    k1 <- rate(v, tau0)
    vmid <- v + (h / 2) * k1
    k2 <- rate(vmid, tau0 + h / 2)
    if (q == 2L) {
      v <- v + h * k2
    } else if (q == 3L) {
      vmid2 <- v - h * k1 + 2 * h * k2
      k3 <- rate(vmid2, tau0 + h)
      v <- v + (h / 6) * (k1 + 4 * k2 + k3)
    } else {
      vmid2 <- v + (h / 2) * k2
      k3 <- rate(vmid2, tau0 + h / 2)
      vmid3 <- v + h * k3
      k4 <- rate(vmid3, tau0 + h)
      v <- v + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    }
  }
  list(v = v, nf = m * q)
}

.mr_gark_alg <- function(name, m, implicit = FALSE) {
  m <- .mr_check_m(m, name)
  tab <- .mr_gark_tableau(name)
  if (!implicit) {
    return(.new_alg(name, tab$order,
                    step = function(ctx, t, u, dt, cache) .mr_gark_step(ctx, t, u, dt, cache, tab, m, name),
                    init = .mr_init_fn(name), on_event = .mr_on_event, m = m))
  }
  .new_alg(name, tab$order,
           step = function(ctx, t, u, dt, cache) .mr_gark_step(ctx, t, u, dt, cache, tab, m, name),
           init = .mr_gark_implicit_init(name), on_event = .mr_on_event, m = m, stiff = TRUE,
           qsteady_max = 6 / 5, on_reject = .nl_on_reject)
}

.mr_gark_implicit_init <- function(name) function(ctx, t, u) {
  .mr_check_split(ctx, name)
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  cache$fsal <- ctx$f1(u, ctx$p, t) + ctx$f2(u, ctx$p, t)
  cache$nf0 <- 1L
  cache$nf20 <- 1L
  cache
}

#' @rdname multirate_methods
#' @export
MRIGARKERK22a <- function(m) .mr_gark_alg("MRIGARKERK22a", m)

#' @rdname multirate_methods
#' @export
MRIGARKERK22b <- function(m) .mr_gark_alg("MRIGARKERK22b", m)

#' @rdname multirate_methods
#' @export
MRIGARKERK33a <- function(m) .mr_gark_alg("MRIGARKERK33a", m)

#' @rdname multirate_methods
#' @export
MRIGARKERK45a <- function(m) .mr_gark_alg("MRIGARKERK45a", m)

#' @rdname multirate_methods
#' @export
MRIGARKIRK21a <- function(m, autodiff = NULL) .with_autodiff(autodiff, .mr_gark_alg("MRIGARKIRK21a", m, TRUE))

#' @rdname multirate_methods
#' @export
MRIGARKESDIRK34a <- function(m, autodiff = NULL)
  .with_autodiff(autodiff, .mr_gark_alg("MRIGARKESDIRK34a", m, TRUE))

#' @rdname multirate_methods
#' @export
MRIGARKESDIRK46a <- function(m, autodiff = NULL)
  .with_autodiff(autodiff, .mr_gark_alg("MRIGARKESDIRK46a", m, TRUE))

# perform_step! of MRIGARKConstantCache and MRIGARKImplicitConstantCache.
# The implicit stages solve z = tmp + gamma dt f2(z) by the DIRK form of
# NLNewton on f2 (nlsolve_f of an MRI-GARK implicit method is f2), with the
# finite difference Jacobian of f2.
.mr_gark_step <- function(ctx, t, u, dt, cache, tab, m, name) {
  f1 <- ctx$f1; f2 <- ctx$f2; p <- ctx$p
  fsal0 <- f1(u, p, t) + f2(u, p, t)
  nf <- 1L; nf2 <- 1L
  cnt <- list(nf = 0L, njacs = 0L, nsolve = 0L, nw = 0L)
  s <- length(tab$dc)
  q <- tab$q
  z <- vector("list", s + 1L)
  fS <- vector("list", s)
  z[[1L]] <- u
  cprev <- 0
  implicit <- any(tab$g0 != 0)
  if (implicit) {
    nl <- cache$nl
    gamma <- tab$g0[tab$g0 != 0][1L]
    ctx2 <- ctx
    ctx2$f <- f2
    ctx2$jac <- NULL
    ctx2$tgrad <- NULL
    first <- TRUE
  }
  for (i in seq_len(s)) {
    fS[[i]] <- f2(z[[i]], p, .fma(cprev, dt, t))
    nf2 <- nf2 + 1L
    cnext <- cprev + tab$dc[i]
    if (tab$g0[i] == 0) {
      r <- .mr_gark_substage(z[[i]], f1, p, t, dt, tab$dc[i], cprev, m, q, tab$W0[i, ], tab$W1[i, ],
                             fS, i)
      z[[i + 1L]] <- r$v
      nf <- nf + r$nf
    } else {
      tmp <- z[[i]]
      for (j in seq_len(i)) {
        wb <- tab$W0[i, j] + tab$W1[i, j] / 2
        if (wb == 0) next
        tmp <- .fma(dt * wb, fS[[j]], tmp)
      }
      r <- .nlsolve(nl, ctx2, t, dt, dt * fS[[i]], tmp, gamma, 1, .fma(cnext, dt, t), u,
                    method = "dirk", firststage = first)
      first <- FALSE
      for (nm in names(cnt)) cnt[[nm]] <- cnt[[nm]] + as.integer(r[[nm]] %||% 0L)
      if (!r$ok)
        return(list(fail = TRUE, u = u, nf = nf + cnt$nf, nf2 = nf2, njacs = cnt$njacs,
                    nsolve = cnt$nsolve, nw = cnt$nw, cache = cache))
      z[[i + 1L]] <- .fma(gamma, r$z, tmp)
    }
    cprev <- cnext
  }
  unew <- z[[s + 1L]]
  EEst <- 0
  if (ctx$adaptive) {
    ut <- if (is.null(tab$Wemb0)) unew - z[[s]] else {
      r <- .mr_gark_substage(z[[s]], f1, p, t, dt, tab$dc[s], cprev - tab$dc[s], m, q, tab$Wemb0,
                             tab$Wemb1, fS, s)
      nf <- nf + r$nf
      unew - r$v
    }
    EEst <- .err_norm(ut, u, unew, ctx$abstol, ctx$reltol)
  }
  fsal1 <- f1(unew, p, t + dt) + f2(unew, p, t + dt)
  list(u = unew, EEst = EEst, k = list(fsal0, fsal1), cache = if (implicit) cache else list(),
       nf = nf + 1L + cnt$nf, nf2 = nf2 + 1L, njacs = cnt$njacs, nsolve = cnt$nsolve, nw = cnt$nw)
}

# ---------------------------------------------------------------------------
# MIS: the MIS2 tableau of Wensch, Knoth and Galant. d (row sums of beta),
# c = (I - alpha - gamma) \ d, ctilde = alpha c and bhat = b - btilde are
# the values MIS2Tableau(Float64) computes.

.mr_mis2_tableau <- function() {
  alpha <- matrix(0, 4L, 4L)
  alpha[3, 2] <- 0.53694656671
  alpha[4, 2] <- 0.480892968551; alpha[4, 3] <- 0.500561163566
  beta <- matrix(0, 4L, 4L)
  beta[2, 1] <- 0.126848494553
  beta[3, 1] <- -0.784838278826; beta[3, 2] <- 1.37442675268
  beta[4, 1] <- -0.0456727081749; beta[4, 2] <- -0.0087508227119; beta[4, 3] <- 0.524775788629
  gamma <- matrix(0, 4L, 4L)
  gamma[3, 2] <- 0.652465126004
  gamma[4, 2] <- -0.0732769849457; gamma[4, 3] <- 0.14490243042
  list(alpha = alpha, beta = beta, gamma = gamma,
       d = c(0, 0.126848494553, 0.589588473854, 0.4703522577422),
       c = c(0, 0.126848494553, 0.7404635564785064, 0.9999999999985448),
       ctilde = c(0, 0, 0.0681108636425655, 0.43164784851091725),
       bhat = c(-1.4031673973738952, 0.87839160874344, 0.524775788629, 0))
}

#' @rdname multirate_methods
#' @export
MIS <- function(m) {
  m <- .mr_check_m(m, "MIS")
  .new_alg("MIS", 2L, step = function(ctx, t, u, dt, cache) .mr_mis_step(ctx, t, u, dt, cache, m),
           init = .mr_init_fn("MIS"), on_event = .mr_on_event, m = m)
}

.mr_mis_step <- function(ctx, t, u, dt, cache, m) {
  f1 <- ctx$f1; f2 <- ctx$f2; p <- ctx$p
  tab <- .mr_mis2_tableau()
  fsal0 <- f1(u, p, t) + f2(u, p, t)
  nf <- 1L; nf2 <- 1L
  s <- length(tab$d)
  Y <- vector("list", s)
  fS <- vector("list", s)
  Y[[1L]] <- u
  fS[[1L]] <- f2(u, p, t)
  nf2 <- nf2 + 1L
  for (i in 2:s) {
    di <- tab$d[i]
    inv_di <- 1 / di
    v <- u
    for (j in seq_len(i - 1L)) if (tab$alpha[i, j] != 0) v <- .fma(tab$alpha[i, j], Y[[j]] - u, v)
    off <- numeric(length(u))
    for (j in seq_len(i - 1L)) {
      if (tab$gamma[i, j] != 0) off <- .fma(tab$gamma[i, j] * inv_di / dt, Y[[j]] - u, off)
      if (tab$beta[i, j] != 0) off <- .fma(tab$beta[i, j] * inv_di, fS[[j]], off)
    }
    M_in <- max(1, ceiling(m * di))
    h_in <- di * dt / M_in
    slope <- (tab$c[i] - tab$ctilde[i]) * inv_di
    t0 <- .fma(tab$ctilde[i], dt, t)
    for (ks in seq_len(M_in)) {
      tau <- (ks - 1) * h_in
      kf <- f1(v, p, .fma(slope, tau, t0))
      vmid <- .fma(h_in / 2, off + kf, v)
      kf <- f1(vmid, p, .fma(slope, tau + h_in / 2, t0))
      v <- .fma(h_in, off + kf, v)
    }
    nf <- nf + 2L * as.integer(M_in)
    Y[[i]] <- v
    fS[[i]] <- f2(v, p, .fma(tab$c[i], dt, t))
    nf2 <- nf2 + 1L
  }
  unew <- Y[[s]]
  EEst <- 0
  if (ctx$adaptive) {
    ut <- (dt * tab$bhat[1L]) * fS[[1L]]
    for (j in 2:s) if (tab$bhat[j] != 0) ut <- .fma(dt * tab$bhat[j], fS[[j]], ut)
    EEst <- .err_norm(ut, u, unew, ctx$abstol, ctx$reltol)
  }
  fsal1 <- f1(unew, p, t + dt) + f2(unew, p, t + dt)
  list(u = unew, EEst = EEst, k = list(fsal0, fsal1), cache = list(), nf = nf + 1L, nf2 = nf2 + 1L)
}

# ---------------------------------------------------------------------------
# MREIL: MREEF with linearly implicit Euler substeps of the fast part,
# (M - h J) delta = h (f2 + f1), solved as W delta = -(f2 + f1) with
# W = J - M / h (jacobian2W); J is the Jacobian of f1 at the start of the
# step (the out of place cache evaluates it every attempt).

#' @rdname multirate_methods
#' @export
MREIL <- function(m = 4L, order = 4L, seq = "harmonic", autodiff = NULL) {
  m <- .mr_check_m(m, "MREIL")
  if (!is.numeric(order) || length(order) != 1L || is.na(order) || order < 2)
    stop("MREIL: `order` must be >= 2", call. = FALSE)
  order <- as.integer(order)
  seq <- .mr_check_seq(seq)
  .with_autodiff(autodiff, .new_alg(
    "MREIL", order, step = function(ctx, t, u, dt, cache) .mr_mreil_step(ctx, t, u, dt, cache, m, order, seq),
    init = .mr_init_fn("MREIL"), on_event = .mr_on_event, stiff = TRUE, qsteady_max = 6 / 5, m = m,
    extrapolation_order = order, seq = seq, mass_matrix = TRUE))
}

.mr_mreil_step <- function(ctx, t, u, dt, cache, m, order, seq) {
  f1 <- ctx$f1; f2 <- ctx$f2; p <- ctx$p
  ns <- .mr_ext_seq(seq, order)
  ctx1 <- ctx
  ctx1$f <- f1
  ctx1$jac <- ctx$jac1
  jj <- .calc_J(ctx1, u, t)
  J <- jj$J
  nf <- jj$nf; nf2 <- 0L; nw <- 0L; nsolve <- 0L
  out <- function(...) list(nf = nf, nf2 = nf2, njacs = 1L, nw = nw, nsolve = nsolve, ...)
  Tc <- vector("list", order)
  fsal0 <- NULL
  for (j in seq_len(order)) {
    nj <- ns[j]
    h_mac <- dt / nj
    h_fast <- h_mac / m
    W <- if (is.null(ctx$mass)) J - diag(1 / h_fast, nrow(J)) else J - ctx$mass * (1 / h_fast)
    nw <- nw + 1L
    Wf <- .factorize(W)
    # a singular W fails the attempt (force_stepfail)
    if (is.null(Wf)) return(out(fail = TRUE, u = u, cache = list()))
    uc <- u
    for (i_mac in seq_len(nj)) {
      t_mac <- .fma(i_mac - 1, h_mac, t)
      ks <- f2(uc, p, t_mac)
      nf2 <- nf2 + 1L
      for (i_fast in seq_len(m)) {
        t_fast <- .fma(i_fast - 1, h_fast, t_mac)
        kf <- f1(uc, p, t_fast)
        nf <- nf + 1L
        if (is.null(fsal0)) fsal0 <- ks + kf
        dz <- .wsolve(Wf, -(ks + kf))
        nsolve <- nsolve + 1L
        uc <- uc + dz
      }
    }
    Tc[[j]] <- uc
  }
  # the ratio is a Rational in OrdinaryDiffEq: ratio - 1 rounds once
  for (k in seq_len(order - 1L)) {
    for (j in order:(k + 1L)) {
      Tc[[j]] <- Tc[[j]] + (Tc[[j]] - Tc[[j - 1L]]) / ((ns[j] - ns[j - k]) / ns[j - k])
    }
  }
  unew <- Tc[[order]]
  EEst <- if (ctx$adaptive) .err_norm(Tc[[order]] - Tc[[order - 1L]], u, unew, ctx$abstol, ctx$reltol) else 0
  fsal1 <- f1(unew, p, t + dt) + f2(unew, p, t + dt)
  nf <- nf + 1L
  nf2 <- nf2 + 1L
  out(u = unew, EEst = EEst, k = list(fsal0, fsal1), cache = list())
}
