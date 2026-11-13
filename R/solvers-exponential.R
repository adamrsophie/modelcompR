# Exponential Runge-Kutta, EPIRK and exponential Rosenbrock methods, after
# OrdinaryDiffEqExponentialRK.jl 2.4 (algorithms.jl, alg_utils.jl,
# exponential_rk_caches.jl, exponential_rk_perform_step.jl).
#
# The steps follow the out of place perform_step! (the ConstantCache
# methods). The linear part is the operator of a SplitODEProblem (whose
# second function is the nonlinear part) or, for any ODEProblem, the
# Jacobian at the start of the step. With krylov = FALSE the classical
# methods precompute exp and phi functions of dt times the operator once,
# from the operator as constructed and the dt passed to solve(), as
# expRK_operators does when the cache is built.

.exprk_is_split <- function(ctx) inherits(ctx$fobj, "SplitFunction")

.exprk_operator <- function(ctx, name) {
  L <- ctx$fobj$f1
  if (!inherits(L, "AbstractSciMLOperator"))
    stop(name, "() needs the linear part of a SplitODEProblem to be a SciML operator, ",
         "such as MatrixOperator(A)", call. = FALSE)
  L
}

# convert(Number or AbstractMatrix, f.f1.f): a number, a diagonal (class
# "eu_diagonal", whose phi functions are taken entrywise) or a matrix
.exprk_concrete <- function(L) {
  if (inherits(L, "ScalarOperator")) return(L$value)
  if (inherits(L, "DiagonalOperator")) return(structure(.as_num(L$diag), class = "eu_diagonal"))
  .sop_matrix(L)
}

# products of the precomputed operators (numbers or matrices)
.exprk_dense <- function(A) if (inherits(A, "eu_diagonal")) diag(unclass(A), length(A)) else A
.exprk_mv <- function(X, v) if (is.matrix(X)) .as_num(X %*% v) else X * v
.exprk_mm <- function(X, Y) if (is.matrix(X)) X %*% Y else X * Y

# expRK_operators(alg, dt, A)
.exprk_operators <- function(name, dt, A) {
  P <- function(c, k) {
    X <- if (inherits(A, "eu_diagonal")) structure(c * unclass(A), class = "eu_diagonal") else c * A
    lapply(.eu_phi(X, k), .exprk_dense)
  }
  Ad <- .exprk_dense(A)
  switch(name,
    LawsonEuler = list(if (is.matrix(Ad)) .eu_exp_base(dt * Ad) else exp(dt * Ad)),
    NorsettEuler = list(P(dt, 1L)[[2L]]),
    ETDRK2 = P(dt, 2L)[2:3],
    ETDRK3 = {
      Ph <- P(dt / 2, 1L)
      Pd <- P(dt, 3L)
      list(0.5 * Ph[[2L]], Pd[[2L]], 4 * Pd[[4L]] - 3 * Pd[[3L]] + Pd[[2L]],
           -8 * Pd[[4L]] + 4 * Pd[[3L]], 4 * Pd[[4L]] - Pd[[3L]])
    },
    ETDRK4 = {
      Pd <- P(dt, 3L)
      Ph <- P(dt / 2, 1L)
      list(0.5 * Ph[[2L]], .exprk_mm((dt / 4) * Ad, .exprk_mm(Ph[[2L]], Ph[[2L]])), Ph[[2L]],
           Pd[[2L]] - 3 * Pd[[3L]] + 4 * Pd[[4L]], 2 * Pd[[3L]] - 4 * Pd[[4L]], -Pd[[3L]] + 4 * Pd[[4L]])
    },
    HochOst4 = {
      Pd <- P(dt, 3L)
      Ph <- P(dt / 2, 3L)
      A21 <- 0.5 * Ph[[2L]]
      A52 <- 0.5 * Ph[[3L]] - Pd[[4L]] + 0.25 * Pd[[3L]] - 0.5 * Ph[[4L]]
      A54 <- 0.25 * Ph[[3L]] - A52
      list(A21, A21 - Ph[[3L]], Ph[[3L]], Pd[[2L]] - 2 * Pd[[3L]], Pd[[3L]],
           0.5 * Ph[[2L]] - 2 * A52 - A54, A52, A54,
           Pd[[2L]] - 3 * Pd[[3L]] + 4 * Pd[[4L]], -Pd[[3L]] + 4 * Pd[[4L]], 4 * Pd[[3L]] - 8 * Pd[[4L]])
    },
    Friedli = {
      Pd <- P(dt, 3L)
      Ph <- P(dt / 2, 2L)
      A21 <- 0.5 * Ph[[2L]]
      list(A21, A21 - 0.5 * Ph[[3L]], 0.5 * Ph[[3L]], Pd[[2L]] - 2 * Pd[[3L]],
           (-26 / 25) * Pd[[2L]] + (2 / 25) * Pd[[3L]], (26 / 25) * Pd[[2L]] + (48 / 25) * Pd[[3L]],
           Pd[[2L]] - 3 * Pd[[3L]] + 4 * Pd[[4L]], 4 * Pd[[3L]] - 8 * Pd[[4L]], -Pd[[3L]] + 4 * Pd[[4L]])
    },
    ETD2 = {
      Pd <- P(dt, 2L)
      list(Pd[[1L]], Pd[[2L]], Pd[[2L]] + Pd[[3L]], -Pd[[3L]])
    })
}

# The cached operators (expRK_operators) of the first step, and those for a
# step of dt: upstream cannot change dt (isdtchangeable is false) and steps
# past the final time; here a step shortened to reach a tstop gets operators
# for its own size.
.exprk_ops_cache <- function(name, cache, dt, L) {
  if (is.null(cache$ops)) {
    cache$ops <- .exprk_operators(name, dt, .exprk_concrete(L))
    cache$ops_dt <- dt
  }
  cache
}

.exprk_ops_for <- function(name, cache, dt, L)
  if (.eu_isapprox(dt, cache$ops_dt)) cache$ops else .exprk_operators(name, dt, .exprk_concrete(L))

# The step context of the classical methods: the linear operator A (the
# split operator, or the Jacobian at the step start) and the nonlinear part
# nl(u, t) (f2 of a split problem, otherwise f(u) - A u).
.exprk_setup <- function(ctx, name, t, u, f0) {
  if (.exprk_is_split(ctx)) {
    A <- .exprk_operator(ctx, name)
    f2 <- .as_rhs(ctx$fobj$f2)
    return(list(A = A, split = TRUE, njacs = 0L, nfj = 0L,
                nl = function(v, tt) f2(v, ctx$p, tt)))
  }
  # the out of place jacobian evaluates f at u again (n + 1 calls for
  # forward differences, as OrdinaryDiffEqDifferentiation counts them)
  J <- .calc_J(ctx, u, t)
  A <- J$J
  list(A = A, split = FALSE, njacs = 1L, nfj = J$nf,
       nl = function(v, tt) ctx$f(v, ctx$p, tt) - .as_num(A %*% v))
}

.exprk_alg <- function(name, order, core, krylov = FALSE, m = 30L, iop = 0L, cached = TRUE,
                       adaptive = FALSE, autodiff = NULL, ...) {
  krylov <- isTRUE(krylov)
  m <- as.integer(m)
  iop <- as.integer(iop)
  init <- function(ctx, t, u) {
    cache <- list(fsal = ctx$f(u, ctx$p, t))
    if (cached && !krylov) {
      if (!.exprk_is_split(ctx))
        stop(name, "(): caching can only be used with SplitFunction; use krylov = TRUE ",
             "or a SplitODEProblem", call. = FALSE)
    }
    cache
  }
  step <- function(ctx, t, u, dt, cache) {
    u <- .as_num(u)
    S <- .exprk_setup(ctx, name, t, u, cache$fsal)
    n <- length(u)
    if (cached && !krylov) {
      cache <- .exprk_ops_cache(name, cache, dt, S$A)
      S$ops <- .exprk_ops_for(name, cache, dt, S$A)
    }
    S$kw <- list(m = min(m, n), iop = iop)
    S$F1 <- cache$fsal
    r <- core(ctx, S, t, u, dt)
    fl <- ctx$f(r$u, ctx$p, t + dt)
    # evaluations of the nonlinear part count as nf2 for split problems
    nnl <- r$nnl %||% 0L
    res <- list(u = r$u, k = list(cache$fsal, fl), cache = list(fsal = fl, ops = cache$ops, ops_dt = cache$ops_dt),
                nf = 1L + S$nfj + (r$nf %||% 0L) + if (S$split) 0L else nnl,
                njacs = S$njacs)
    if (S$split) res$nf2 <- nnl
    res$err <- r$err
    res
  }
  alg <- .new_alg(name, order, adaptive = adaptive, step = step, init = init,
                  krylov = krylov, m = m, iop = iop, ...)
  .with_autodiff(autodiff, alg)
}

# phiv(t, A, b, k; m, iop) of the step context
.exprk_phiv <- function(S, t, b, k) .eu_phiv(t, S$A, b, k, m = S$kw$m, iop = S$kw$iop)
.exprk_arnoldi <- function(S, b) .eu_arnoldi(S$A, b, m = S$kw$m, iop = S$kw$iop)
.exprk_phivk <- function(t, Ks, k) .eu_phiv_ks(t, Ks, k)$w

#' Exponential Runge-Kutta methods
#'
#' Exponential integrators of OrdinaryDiffEqExponentialRK for semilinear
#' problems `u' = A u + f(u, t)`, given as `SplitODEProblem(L, f, u0,
#' tspan)` with a SciML operator `L` (such as `MatrixOperator(A)` or
#' `DiagonalOperator(d)`), or for any `ODEProblem`, whose Jacobian at the
#' start of each step is then the linear part (from the `jac_fn` of the
#' `ODEFunction` when given, otherwise by finite differences).
#'
#' `LawsonEuler()` is the first order Lawson Euler method, `NorsettEuler()`
#' (alias `ETD1()`) the first order exponential Euler method, `ETDRK2()`,
#' `ETDRK3()` and `ETDRK4()` the exponential time differencing Runge-Kutta
#' methods of orders 2 to 4, `HochOst4()` the fourth order method of
#' Hochbruck and Ostermann with stiff order 4, `Friedli()` Friedli's fourth
#' order method, and `ETD2()` the second order exponential time differencing
#' multistep method (split problems only). These are fixed step methods.
#' With `krylov = FALSE` (the default) the exponential and phi functions of
#' the operator are computed once (this needs a `SplitODEProblem`; a last
#' step shortened to reach the final time gets its own operators, where
#' upstream steps past it); with `krylov = TRUE` their actions are
#' approximated in Krylov subspaces of dimension `m`, with incomplete
#' orthogonalization of length `iop` (0 for full Arnoldi; Lanczos for
#' symmetric operators). Either way the steps use the operator as
#' constructed, as upstream does for out of place problems: an update
#' function of the operator only enters through the right-hand side. The
#' Krylov form of `HochOst4()` uses
#' the third stage in its fifth stage, as the precomputed form does
#' (upstream repeats the second there).
#'
#' @param krylov Whether to use Krylov approximations instead of precomputed
#'   operators.
#' @param m Krylov subspace dimension.
#' @param iop Incomplete orthogonalization length (0 for full Arnoldi).
#' @param autodiff An ADType for the Jacobian of an `ODEProblem` without
#'   `jac_fn`.
#' @return An algorithm object.
#' @examples
#' L <- MatrixOperator(matrix(c(-2, 1, 1, -2), 2))
#' prob <- SplitODEProblem(L, function(u, p, t) -u * u * u, c(1, 0.5), c(0, 1))
#' solve(prob, ETDRK4(), dt = 0.1)$u[[11]]
#' solve(prob, LawsonEuler(krylov = TRUE, m = 2), dt = 0.1)$u[[11]]
#' @name exponential_rk_methods
NULL

#' @rdname exponential_rk_methods
#' @export
LawsonEuler <- function(krylov = FALSE, m = 30, iop = 0, autodiff = NULL)
  .exprk_alg("LawsonEuler", 1, function(ctx, S, t, u, dt) {
    nl <- S$nl(u, t)
    v <- u + dt * nl
    list(u = if (krylov) .eu_expv(dt, S$A, v, m = S$kw$m, iop = S$kw$iop) else .exprk_mv(S$ops[[1L]], v),
         nnl = 1L)
  }, krylov, m, iop, autodiff = autodiff)

#' @rdname exponential_rk_methods
#' @export
NorsettEuler <- function(krylov = FALSE, m = 30, iop = 0, autodiff = NULL)
  .exprk_alg("NorsettEuler", 1, function(ctx, S, t, u, dt) {
    list(u = if (krylov) u + dt * .exprk_phiv(S, dt, S$F1, 1L)[, 2L]
             else u + dt * .exprk_mv(S$ops[[1L]], S$F1))
  }, krylov, m, iop, autodiff = autodiff)

#' @rdname exponential_rk_methods
#' @export
ETD1 <- function(krylov = FALSE, m = 30, iop = 0, autodiff = NULL) NorsettEuler(krylov, m, iop, autodiff)

#' @rdname exponential_rk_methods
#' @export
ETDRK2 <- function(krylov = FALSE, m = 30, iop = 0, autodiff = NULL)
  .exprk_alg("ETDRK2", 2, function(ctx, S, t, u, dt) {
    F1 <- S$F1
    if (krylov) {
      w1 <- .exprk_phiv(S, dt, F1, 2L)
      U2 <- u + dt * w1[, 2L]
      F2 <- S$nl(U2, t + dt) + .eu_mul(S$A, u)
      w2 <- .exprk_phiv(S, dt, F2, 2L)
      return(list(u = u + dt * (w1[, 2L] - w1[, 3L] + w2[, 3L]), nnl = 1L))
    }
    G1 <- S$nl(u, t)
    U2 <- u + dt * .exprk_mv(S$ops[[1L]], F1)
    G2 <- S$nl(U2, t + dt)
    list(u = U2 + dt * .exprk_mv(S$ops[[2L]], G2 - G1), nnl = 2L)
  }, krylov, m, iop, autodiff = autodiff)

#' @rdname exponential_rk_methods
#' @export
ETDRK3 <- function(krylov = FALSE, m = 30, iop = 0, autodiff = NULL)
  .exprk_alg("ETDRK3", 3, function(ctx, S, t, u, dt) {
    Au <- .eu_mul(S$A, u)
    F1 <- S$F1
    if (krylov) {
      Ks <- .exprk_arnoldi(S, F1)
      w1_half <- .exprk_phivk(dt / 2, Ks, 1L)
      w1 <- .exprk_phivk(dt, Ks, 3L)
      U2 <- u + dt / 2 * w1_half[, 2L]
      F2 <- S$nl(U2, t + dt / 2) + Au
      w2 <- .exprk_phiv(S, dt, F2, 3L)
      U3 <- u + dt * (2 * w2[, 2L] - w1[, 2L])
      F3 <- S$nl(U3, t + dt) + Au
      w3 <- .exprk_phiv(S, dt, F3, 3L)
      return(list(u = u + dt * (4 * w1[, 4L] - 3 * w1[, 3L] + w1[, 2L] - 8 * w2[, 4L] + 4 * w2[, 3L] +
                                  4 * w3[, 4L] - w3[, 3L]), nnl = 2L))
    }
    o <- S$ops
    U2 <- u + dt * .exprk_mv(o[[1L]], F1)
    F2 <- S$nl(U2, t + dt / 2) + Au
    U3 <- u + dt * .exprk_mv(o[[2L]], 2 * F2 - F1)
    F3 <- S$nl(U3, t + dt) + Au
    list(u = u + dt * (.exprk_mv(o[[3L]], F1) + .exprk_mv(o[[4L]], F2) + .exprk_mv(o[[5L]], F3)), nnl = 2L)
  }, krylov, m, iop, autodiff = autodiff)

#' @rdname exponential_rk_methods
#' @export
ETDRK4 <- function(krylov = FALSE, m = 30, iop = 0, autodiff = NULL)
  .exprk_alg("ETDRK4", 4, function(ctx, S, t, u, dt) {
    Au <- .eu_mul(S$A, u)
    F1 <- S$F1
    halfdt <- dt / 2
    if (krylov) {
      Ks <- .exprk_arnoldi(S, F1)
      w1_half <- .exprk_phivk(halfdt, Ks, 1L)
      w1 <- .exprk_phivk(dt, Ks, 3L)
      U2 <- u + halfdt * w1_half[, 2L]
      F2 <- S$nl(U2, t + halfdt) + Au
      Ks <- .exprk_arnoldi(S, F2)
      w2_half <- .exprk_phivk(halfdt, Ks, 1L)
      w2 <- .exprk_phivk(dt, Ks, 3L)
      U3 <- u + halfdt * w2_half[, 2L]
      F3 <- S$nl(U3, t + halfdt) + Au
      w3 <- .exprk_phiv(S, dt, F3, 3L)
      rtmp <- 2 * F3 - F1 - Au + .eu_mul(S$A, U2)
      wtmp <- .exprk_phiv(S, halfdt, rtmp, 1L)
      U4 <- U2 + halfdt * wtmp[, 2L]
      F4 <- S$nl(U4, t + dt) + Au
      w4 <- .exprk_phiv(S, dt, F4, 3L)
      return(list(u = u + dt * (w1[, 2L] - 3 * w1[, 3L] + 4 * w1[, 4L] + 2 * w2[, 3L] - 4 * w2[, 4L] +
                                  2 * w3[, 3L] - 4 * w3[, 4L] + 4 * w4[, 4L] - w4[, 3L]), nnl = 3L))
    }
    o <- S$ops
    U2 <- u + dt * .exprk_mv(o[[1L]], F1)
    F2 <- S$nl(U2, t + halfdt) + Au
    U3 <- u + dt * .exprk_mv(o[[1L]], F2)
    F3 <- S$nl(U3, t + halfdt) + Au
    U4 <- u + dt * (.exprk_mv(o[[2L]], F1) + .exprk_mv(o[[3L]], F3))
    F4 <- S$nl(U4, t + dt) + Au
    list(u = u + dt * (.exprk_mv(o[[4L]], F1) + .exprk_mv(o[[5L]], F2 + F3) + .exprk_mv(o[[6L]], F4)),
         nnl = 3L)
  }, krylov, m, iop, autodiff = autodiff)

#' @rdname exponential_rk_methods
#' @export
HochOst4 <- function(krylov = FALSE, m = 30, iop = 0, autodiff = NULL)
  .exprk_alg("HochOst4", 4, function(ctx, S, t, u, dt) {
    Au <- .eu_mul(S$A, u)
    F1 <- S$F1
    halfdt <- dt / 2
    if (krylov) {
      Ks <- .exprk_arnoldi(S, F1)
      w1_half <- .exprk_phivk(halfdt, Ks, 3L)
      w1 <- .exprk_phivk(dt, Ks, 3L)
      U2 <- u + halfdt * w1_half[, 2L]
      F2 <- S$nl(U2, t + halfdt) + Au
      Ks <- .exprk_arnoldi(S, F2)
      w2_half <- .exprk_phivk(halfdt, Ks, 3L)
      w2 <- .exprk_phivk(dt, Ks, 3L)
      U3 <- u + dt * (0.5 * w1_half[, 2L] - w1_half[, 3L] + w2_half[, 3L])
      F3 <- S$nl(U3, t + halfdt) + Au
      Ks <- .exprk_arnoldi(S, F3)
      w3_half <- .exprk_phivk(halfdt, Ks, 3L)
      w3 <- .exprk_phivk(dt, Ks, 3L)
      U4 <- u + dt * (w1[, 2L] - 2 * w1[, 3L] + w2[, 3L] + w3[, 3L])
      F4 <- S$nl(U4, t + dt) + Au
      Ks <- .exprk_arnoldi(S, F4)
      w4_half <- .exprk_phivk(halfdt, Ks, 3L)
      w4 <- .exprk_phivk(dt, Ks, 3L)
      # A51 F1 + A52 (F2 + F3) + A54 F4 of the cached operators; upstream
      # has w2[, 4] in place of w3[, 4] in the F3 group (an indexing slip)
      U5 <- u + dt * (0.5 * w1_half[, 2L] - 0.75 * w1_half[, 3L] + 0.5 * w1_half[, 4L] + w1[, 4L] -
                        0.25 * w1[, 3L] +
                        0.5 * w2_half[, 3L] - w2[, 4L] + 0.25 * w2[, 3L] - 0.5 * w2_half[, 4L] +
                        0.5 * w3_half[, 3L] - w3[, 4L] + 0.25 * w3[, 3L] - 0.5 * w3_half[, 4L] +
                        w4[, 4L] - 0.25 * w4[, 3L] - 0.25 * w4_half[, 3L] + 0.5 * w4_half[, 4L])
      F5 <- S$nl(U5, t + halfdt) + Au
      w5 <- .exprk_phiv(S, dt, F5, 3L)
      return(list(u = u + dt * (w1[, 2L] - 3 * w1[, 3L] + 4 * w1[, 4L] - w4[, 3L] + 4 * w4[, 4L] +
                                  4 * w5[, 3L] - 8 * w5[, 4L]), nnl = 4L))
    }
    o <- S$ops
    U2 <- u + dt * .exprk_mv(o[[1L]], F1)
    F2 <- S$nl(U2, t + halfdt) + Au
    U3 <- u + dt * (.exprk_mv(o[[2L]], F1) + .exprk_mv(o[[3L]], F2))
    F3 <- S$nl(U3, t + halfdt) + Au
    U4 <- u + dt * (.exprk_mv(o[[4L]], F1) + .exprk_mv(o[[5L]], F2 + F3))
    F4 <- S$nl(U4, t + dt) + Au
    U5 <- u + dt * (.exprk_mv(o[[6L]], F1) + .exprk_mv(o[[7L]], F2 + F3) + .exprk_mv(o[[8L]], F4))
    F5 <- S$nl(U5, t + halfdt) + Au
    list(u = u + dt * (.exprk_mv(o[[9L]], F1) + .exprk_mv(o[[10L]], F4) + .exprk_mv(o[[11L]], F5)),
         nnl = 4L)
  }, krylov, m, iop, autodiff = autodiff)

#' @rdname exponential_rk_methods
#' @export
Friedli <- function(krylov = FALSE, m = 30, iop = 0, autodiff = NULL)
  .exprk_alg("Friedli", 4, function(ctx, S, t, u, dt) {
    Au <- .eu_mul(S$A, u)
    F1 <- S$F1
    halfdt <- dt / 2
    if (krylov) {
      Ks <- .exprk_arnoldi(S, F1)
      w1_half <- .exprk_phivk(halfdt, Ks, 2L)
      w1 <- .exprk_phivk(dt, Ks, 3L)
      U2 <- u + halfdt * w1_half[, 2L]
      F2 <- S$nl(U2, t + halfdt) + Au
      Ks <- .exprk_arnoldi(S, F2)
      w2_half <- .exprk_phivk(halfdt, Ks, 2L)
      w2 <- .exprk_phivk(dt, Ks, 2L)
      U3 <- u + halfdt * (w1_half[, 2L] - w1_half[, 3L] + w2_half[, 3L])
      F3 <- S$nl(U3, t + halfdt) + Au
      w3 <- .exprk_phiv(S, dt, F3, 3L)
      U4 <- u + dt * (w1[, 2L] - 2 * w1[, 3L] - (26 / 25) * w2[, 2L] + (2 / 25) * w2[, 3L] +
                        (26 / 25) * w3[, 2L] + (48 / 25) * w3[, 3L])
      F4 <- S$nl(U4, t + dt) + Au
      w4 <- .exprk_phiv(S, dt, F4, 3L)
      return(list(u = u + dt * (w1[, 2L] - 3 * w1[, 3L] + 4 * w1[, 4L] + 4 * w3[, 3L] - 8 * w3[, 4L] -
                                  w4[, 3L] + 4 * w4[, 4L]), nnl = 3L))
    }
    o <- S$ops
    U2 <- u + dt * .exprk_mv(o[[1L]], F1)
    F2 <- S$nl(U2, t + halfdt) + Au
    U3 <- u + dt * (.exprk_mv(o[[2L]], F1) + .exprk_mv(o[[3L]], F2))
    F3 <- S$nl(U3, t + halfdt) + Au
    U4 <- u + dt * (.exprk_mv(o[[4L]], F1) + .exprk_mv(o[[5L]], F2) + .exprk_mv(o[[6L]], F3))
    F4 <- S$nl(U4, t + dt) + Au
    list(u = u + dt * (.exprk_mv(o[[7L]], F1) + .exprk_mv(o[[8L]], F3) + .exprk_mv(o[[9L]], F4)),
         nnl = 3L)
  }, krylov, m, iop, autodiff = autodiff)

#' @rdname exponential_rk_methods
#' @export
ETD2 <- function() {
  name <- "ETD2"
  init <- function(ctx, t, u) {
    if (!.exprk_is_split(ctx)) stop("ETD2() needs a SplitODEProblem", call. = FALSE)
    L <- .exprk_operator(ctx, name)
    lin <- .op_as_function(L)(u, ctx$p, t)
    nl <- .as_rhs(ctx$fobj$f2)(u, ctx$p, t)
    list(lin = lin, nl = nl, nlprev = 0 * nl, fsal = lin + nl, nf20 = 1L)
  }
  step <- function(ctx, t, u, dt, cache) {
    u <- .as_num(u)
    L <- .exprk_operator(ctx, name)
    cache <- .exprk_ops_cache(name, cache, dt, L)
    o <- .exprk_ops_for(name, cache, dt, L)
    k1 <- cache$lin + cache$nl
    unew <- if (ctx$iter == 1L) .exprk_mv(o[[1L]], u) + dt * .exprk_mv(o[[2L]], cache$nl)
            else .exprk_mv(o[[1L]], u) + dt * (.exprk_mv(o[[3L]], cache$nl) + .exprk_mv(o[[4L]], cache$nlprev))
    lin <- .op_as_function(L)(unew, ctx$p, t + dt)
    nl <- .as_rhs(ctx$fobj$f2)(unew, ctx$p, t + dt)
    fl <- lin + nl
    list(u = unew, k = list(k1, fl), nf = 1L, nf2 = 1L,
         cache = list(lin = lin, nl = nl, nlprev = cache$nl, fsal = fl, ops = cache$ops, ops_dt = cache$ops_dt))
  }
  .new_alg(name, 2, adaptive = FALSE, step = step, init = init, krylov = FALSE)
}

# ---------------------------------------------------------------------------
# EPIRK methods (Krylov with adaptive time stepping)

.epirk_alg <- function(name, order, core, adaptive_krylov = TRUE, m = 30L, iop = 0L, autodiff = NULL) {
  adaptive_krylov <- isTRUE(adaptive_krylov)
  m <- as.integer(m)
  iop <- as.integer(iop)
  step <- function(ctx, t, u, dt, cache) {
    u <- .as_num(u)
    f0 <- cache$fsal
    J <- .calc_J(ctx, u, t)
    Jm <- J$J
    n <- length(u)
    # phiv_timestep(ts, J, B; tol = reltol, iop, opnorm, adaptive, tau, m)
    pts <- function(ts, B) .eu_phiv_timestep(ts, Jm, B, tol = ctx$reltol, iop = iop,
                                             opnorm = .eu_opnorminf, adaptive = adaptive_krylov,
                                             tau = if (adaptive_krylov) 0 else dt, m = min(m, n))
    Jv <- function(v) .as_num(Jm %*% v)
    fe <- function(v, tt) ctx$f(v, ctx$p, tt)
    r <- core(t, u, dt, f0, pts, Jv, fe)
    fl <- ctx$f(r$u, ctx$p, t + dt)
    list(u = r$u, k = list(f0, fl), cache = list(fsal = fl), nf = 1L + J$nf + r$nf, njacs = 1L)
  }
  alg <- .new_alg(name, order, adaptive = FALSE, step = step, adaptive_krylov = adaptive_krylov,
                  m = m, iop = iop)
  .with_autodiff(autodiff, alg)
}

#' EPIRK and exponential Rosenbrock methods
#'
#' Exponential propagation iterative (EPIRK) methods of
#' OrdinaryDiffEqExponentialRK, which use the Jacobian at the start of each
#' step as the linear part (from the `jac_fn` of the `ODEFunction` when
#' given, otherwise by finite differences) and compute the phi function actions by
#' the adaptive Krylov time stepping of Niesen and Wright, with tolerance
#' `reltol`. `Exp4()` is the fourth order method of Hochbruck, Lubich and
#' Selhofer, `EPIRK4s3A()` and `EPIRK4s3B()` fourth order and `EPIRK5s3()`,
#' `EXPRB53s3()`, `EPIRK5P1()` and `EPIRK5P2()` fifth order EPIRK methods
#' (upstream marks `EPIRK5s3()` as broken: it does not reach fifth order).
#' These are fixed step methods. `Exprb32()` and `Exprb43()` are adaptive
#' exponential Rosenbrock methods of orders 3 and 4 with embedded error
#' estimates, using fixed Krylov subspaces of dimension `m`. For a split
#' problem they also linearize the whole right-hand side, with the remainder
#' `f(u) - J u` (upstream takes the second split function as the remainder,
#' which does not match that Jacobian).
#'
#' @param adaptive_krylov Whether to use the adaptive Krylov time stepping
#'   (otherwise one Krylov step of length `dt`).
#' @param m Krylov subspace dimension (the initial one when adaptive).
#' @param iop Incomplete orthogonalization length (0 for full Arnoldi).
#' @param autodiff An ADType for the Jacobian when the problem has no
#'   `jac_fn`.
#' @return An algorithm object.
#' @examples
#' f <- function(u, p, t) c(-2 * u[1] + u[2], u[1] - 2 * u[2] - u[2] * u[2] * u[2])
#' jac <- function(u, p, t) matrix(c(-2, 1, 1, -2 - 3 * u[2] * u[2]), 2)
#' prob <- ODEProblem(ODEFunction(f, jac_fn = jac), c(1, 0.5), c(0, 1))
#' solve(prob, EPIRK4s3A(), dt = 0.1)$u[[11]]
#' solve(prob, Exprb43())$retcode
#' @name epirk_methods
NULL

#' @rdname epirk_methods
#' @export
Exp4 <- function(adaptive_krylov = TRUE, m = 30, iop = 0, autodiff = NULL)
  .epirk_alg("Exp4", 4, function(t, u, dt, f0, pts, Jv, fe) {
    ts <- c(dt / 3, 2 * dt / 3, dt)
    K1 <- pts(ts, cbind(0 * f0, f0))
    for (i in 1:3) K1[, i] <- K1[, i] / ts[i]
    w4 <- .as_num(K1 %*% c(-7 / 300, 97 / 150, -37 / 300))
    u4 <- u + dt * w4
    d4 <- fe(u4, t + dt) - f0 - dt * Jv(w4)
    K2 <- pts(ts, cbind(0 * d4, d4))
    for (i in 1:3) K2[, i] <- K2[, i] / ts[i]
    w7 <- .as_num(K1 %*% c(59 / 300, -7 / 75, 269 / 300)) + .as_num(K2 %*% c(2 / 3, 2 / 3, 2 / 3))
    u7 <- u + dt * w7
    d7 <- fe(u7, t + dt) - f0 - dt * Jv(w7)
    k7 <- pts(ts[1L], cbind(0 * d7, d7))
    k7 <- k7 / ts[1L]
    list(u = u + dt * (K1[, 3L] + K2[, 1L] - 4 / 3 * K2[, 2L] + K2[, 3L] + k7 / 6), nf = 2L)
  }, adaptive_krylov, m, iop, autodiff)

#' @rdname epirk_methods
#' @export
EPIRK4s3A <- function(adaptive_krylov = TRUE, m = 30, iop = 0, autodiff = NULL)
  .epirk_alg("EPIRK4s3A", 4, function(t, u, dt, f0, pts, Jv, fe) {
    K <- pts(c(dt / 2, 2 * dt / 3), cbind(0 * f0, f0))
    U2 <- u + K[, 1L]
    U3 <- u + K[, 2L]
    R2 <- fe(U2, t + dt / 2) - f0 - Jv(K[, 1L])
    R3 <- fe(U3, t + 2 * dt / 3) - f0 - Jv(K[, 2L])
    B <- matrix(0, length(f0), 5L)
    B[, 2L] <- f0
    B[, 4L] <- (32 * R2 - 13.5 * R3) / (dt * dt)
    B[, 5L] <- (-144 * R2 + 81 * R3) / (dt * dt * dt)
    list(u = u + pts(dt, B), nf = 2L)
  }, adaptive_krylov, m, iop, autodiff)

#' @rdname epirk_methods
#' @export
EPIRK4s3B <- function(adaptive_krylov = TRUE, m = 30, iop = 0, autodiff = NULL)
  .epirk_alg("EPIRK4s3B", 4, function(t, u, dt, f0, pts, Jv, fe) {
    K <- pts(c(dt / 2, 3 * dt / 4), cbind(0 * f0, 0 * f0, f0))
    K[, 1L] <- K[, 1L] * (8 / (3 * dt))
    K[, 2L] <- K[, 2L] * (16 / (9 * dt))
    U2 <- u + K[, 1L]
    U3 <- u + K[, 2L]
    R2 <- fe(U2, t + dt / 2) - f0 - Jv(K[, 1L])
    R3 <- fe(U3, t + 3 * dt / 4) - f0 - Jv(K[, 2L])
    B <- matrix(0, length(f0), 5L)
    B[, 2L] <- f0
    B[, 4L] <- (54 * R2 - 16 * R3) / (dt * dt)
    B[, 5L] <- (-324 * R2 + 144 * R3) / (dt * dt * dt)
    list(u = u + pts(dt, B), nf = 2L)
  }, adaptive_krylov, m, iop, autodiff)

#' @rdname epirk_methods
#' @export
EPIRK5s3 <- function(adaptive_krylov = TRUE, m = 30, iop = 0, autodiff = NULL)
  .epirk_alg("EPIRK5s3", 5, function(t, u, dt, f0, pts, Jv, fe) {
    dt2 <- dt * dt
    dt3 <- dt * dt * dt
    B <- matrix(0, length(f0), 4L)
    B[, 3L] <- (55 / (8 * dt)) * f0
    B[, 4L] <- (-3025 / (192 * dt2)) * f0
    k <- pts(48 * dt / 55, B)
    U2 <- u + k
    R2 <- fe(U2, t + 48 * dt / 55) - f0 - Jv(k)
    B[, 2L] <- (53 / 5) * f0
    B[, 3L] <- (-648 / (5 * dt)) * f0
    B[, 4L] <- (2916 / (5 * dt2)) * f0 + (32065 / (1152 * dt2)) * R2
    k <- pts(4 * dt / 9, B)
    U3 <- u + k
    R3 <- fe(U3, t + 4 * dt / 9) - f0 - Jv(k)
    B <- matrix(0, length(f0), 5L)
    B[, 2L] <- f0
    B[, 4L] <- (-166375 / (61056 * dt2)) * R2 + (2187 / (106 * dt2)) * R3
    B[, 5L] <- (499125 / (27136 * dt3)) * R2 - (2187 / (106 * dt3)) * R3
    list(u = u + pts(dt, B), nf = 2L)
  }, adaptive_krylov, m, iop, autodiff)

#' @rdname epirk_methods
#' @export
EXPRB53s3 <- function(adaptive_krylov = TRUE, m = 30, iop = 0, autodiff = NULL)
  .epirk_alg("EXPRB53s3", 5, function(t, u, dt, f0, pts, Jv, fe) {
    dt2 <- dt * dt
    dt3 <- dt * dt * dt
    K <- pts(c(dt / 2, 9 * dt / 10), cbind(0 * f0, f0))
    U2 <- u + K[, 1L]
    R2 <- fe(U2, t + dt / 2) - f0 - Jv(K[, 1L])
    U3 <- u + K[, 2L]
    z <- 0 * R2
    K <- pts(c(dt / 2, 9 * dt / 10), cbind(z, z, z, R2))
    U3 <- U3 + ((216 / (25 * dt2)) * K[, 1L] + (8 / dt2) * K[, 2L])
    R3 <- fe(U3, t + 9 * dt / 10) - f0 - Jv(U3 - u)
    B <- matrix(0, length(f0), 5L)
    B[, 2L] <- f0
    B[, 4L] <- (18 / dt2) * R2 - (250 / (81 * dt2)) * R3
    B[, 5L] <- (-60 / dt3) * R2 + (500 / (27 * dt3)) * R3
    list(u = u + pts(dt, B), nf = 2L)
  }, adaptive_krylov, m, iop, autodiff)

#' @rdname epirk_methods
#' @export
EPIRK5P1 <- function(adaptive_krylov = TRUE, m = 30, iop = 0, autodiff = NULL)
  .epirk_alg("EPIRK5P1", 5, function(t, u, dt, f0, pts, Jv, fe) {
    g11 <- 0.3512959269505819 * dt
    g21 <- 0.8440547201165712 * dt
    g31 <- dt
    g32 <- 0.7111109536436687 * dt
    g33 <- 0.623781119533715 * dt
    a22 <- 2.3773915340418665
    b2 <- 0x1.ca2d3b3871873p+0
    b3 <- 0x1.2b7936858f2e0p+3 / (dt * dt)
    K1 <- pts(c(g11, g21, g31), cbind(0 * f0, f0))
    U1 <- u + K1[, 1L]
    R1 <- fe(U1, t + g11) - f0 - Jv(K1[, 1L])
    k2 <- pts(g32, cbind(0 * R1, R1))
    U2 <- u + K1[, 2L] + a22 * k2
    R2 <- fe(U2, t + g21) - f0 - Jv(U2 - u)
    B <- matrix(0, length(f0), 4L)
    B[, 4L] <- R2 - 2 * R1
    k3 <- pts(g33, B)
    list(u = u + K1[, 3L] + b2 * k2 + b3 * k3, nf = 2L)
  }, adaptive_krylov, m, iop, autodiff)

#' @rdname epirk_methods
#' @export
EPIRK5P2 <- function(adaptive_krylov = TRUE, m = 30, iop = 0, autodiff = NULL)
  .epirk_alg("EPIRK5P2", 5, function(t, u, dt, f0, pts, Jv, fe) {
    g11 <- 0x1.dd7c325a66039p-2 * dt
    g21 <- 0.8821791265336386 * dt
    g31 <- dt
    g32 <- 0.9207491648814004 * dt
    g33 <- 0.7979156183266451 * dt
    a22 <- 2.8062037328933127 / dt
    b2 <- 2.528063102562463 / dt
    b31 <- -0.12848678265700555
    b32 <- -0.16102803317280917 / dt
    b33 <- 5.267263316169096 / (dt * dt)
    K1 <- pts(c(g11, g21, g31), cbind(0 * f0, f0))
    U1 <- u + K1[, 1L]
    R1 <- fe(U1, t + g11) - f0 - Jv(K1[, 1L])
    k2 <- pts(g32, cbind(0 * R1, 0 * R1, R1))
    U2 <- u + K1[, 2L] + a22 * k2
    R2 <- fe(U2, t + g21) - f0 - Jv(U2 - u)
    dR <- R2 - 2 * R1
    k3 <- pts(g33, cbind(0 * dR, b31 * dR, b32 * dR, b33 * dR))
    list(u = u + K1[, 3L] + b2 * k2 + k3, nf = 2L)
  }, adaptive_krylov, m, iop, autodiff)

# ---------------------------------------------------------------------------
# Adaptive exponential Rosenbrock methods

.exprb_alg <- function(name, order, core, m = 30L, iop = 0L, autodiff = NULL) {
  m <- as.integer(m)
  iop <- as.integer(iop)
  step <- function(ctx, t, u, dt, cache) {
    u <- .as_num(u)
    # the linear part is the Jacobian of the whole right-hand side, also for
    # a split problem, so the remainder is f(u) - J u; upstream takes the
    # second split function as the remainder there, which is inconsistent
    # with that Jacobian
    ctx$fobj <- NULL
    S <- .exprk_setup(ctx, name, t, u, cache$fsal)
    S$kw <- list(m = min(m, length(u)), iop = iop)
    S$F1 <- cache$fsal
    r <- core(ctx, S, t, u, dt)
    fl <- ctx$f(r$u, ctx$p, t + dt)
    res <- list(u = r$u, k = list(cache$fsal, fl), cache = list(fsal = fl),
                nf = 1L + S$nfj + r$nnl, njacs = 1L)
    if (ctx$adaptive) res$err <- r$utilde
    res
  }
  alg <- .new_alg(name, order, adaptive = TRUE, step = step, m = m, iop = iop)
  .with_autodiff(autodiff, alg)
}

#' @rdname epirk_methods
#' @export
Exprb32 <- function(m = 30, iop = 0, autodiff = NULL)
  .exprb_alg("Exprb32", 3, function(ctx, S, t, u, dt) {
    F1 <- S$F1
    w1 <- .exprk_phiv(S, dt, F1, 3L)
    U2 <- u + dt * w1[, 2L]
    F2 <- S$nl(U2, t + dt) + .eu_mul(S$A, u)
    w2 <- .exprk_phiv(S, dt, F2, 3L)
    list(u = u + dt * (w1[, 2L] - 2 * w1[, 4L] + 2 * w2[, 4L]),
         utilde = 2 * dt * (-w1[, 4L] + w2[, 4L]), nnl = 1L)
  }, m, iop, autodiff)

#' @rdname epirk_methods
#' @export
Exprb43 <- function(m = 30, iop = 0, autodiff = NULL)
  .exprb_alg("Exprb43", 4, function(ctx, S, t, u, dt) {
    Au <- .eu_mul(S$A, u)
    F1 <- S$F1
    Ks <- .exprk_arnoldi(S, F1)
    w1_half <- .exprk_phivk(dt / 2, Ks, 1L)
    w1 <- .exprk_phivk(dt, Ks, 4L)
    U2 <- u + dt / 2 * w1_half[, 2L]
    F2 <- S$nl(U2, t + dt / 2) + Au
    w2 <- .exprk_phiv(S, dt, F2, 4L)
    U3 <- u + dt * w2[, 2L]
    F3 <- S$nl(U3, t + dt) + Au
    w3 <- .exprk_phiv(S, dt, F3, 4L)
    list(u = u + dt * (w1[, 2L] - 14 * w1[, 4L] + 36 * w1[, 5L] + 16 * w2[, 4L] - 48 * w2[, 5L] -
                         2 * w3[, 4L] + 12 * w3[, 5L]),
         utilde = dt * (36 * w1[, 5L] - 48 * w2[, 5L] + 12 * w3[, 5L]), nnl = 2L)
  }, m, iop, autodiff)
