# Lie group, Magnus and Crouch-Grossman methods for u' = A(u, t) u and the
# exact solution of linear problems, after OrdinaryDiffEqLinear.jl 2.3
# (algorithms.jl, alg_utils.jl, linear_perform_step.jl).
#
# The problem is ODEProblem(L, u0, tspan) for a SciML operator L (usually
# MatrixOperator(A, update_func)), or for CayleyEuler() also a
# SplitODEProblem whose first function is the operator. The steps follow the
# in-place perform_step! of OrdinaryDiffEqLinear (the only one most methods
# have): the operator is updated to A(u, t) before every use, and the step
# ends with the derivative at the new state for the Hermite interpolation.

# The right-hand side function of an operator: L(u, p, t) = A(u, p, t) u. A
# state of length n^2 for an n x n operator is a matrix state (CayleyEuler),
# passed to the coefficient update as a matrix and multiplied as a matrix.
.op_as_function <- function(L) {
  force(L)
  function(u, p, t) {
    n <- L$size[2L]
    if (!is.na(n) && length(u) != n && length(u) == n * n) {
      U <- matrix(u, n, n)
      return(.as_num(.sop_matrix(update_coefficients(L, U, p, t), n) %*% U))
    }
    op_apply(L, u, u, p, t)
  }
}

# A function or an operator as a right-hand side function
.as_rhs <- function(f) if (inherits(f, "AbstractSciMLOperator")) .op_as_function(f) else f

# The ODEFunction of an operator; the operator stays in $operator
.operator_odefunction <- function(L) {
  fn <- ODEFunction(.op_as_function(L))
  fn$operator <- L
  fn
}

.lin_operator <- function(ctx, name) {
  fo <- ctx$fobj
  L <- if (inherits(fo, "SplitFunction")) fo$f1 else fo$operator
  if (!inherits(L, "AbstractSciMLOperator"))
    stop(name, "() needs a problem defined by a SciML operator, such as ",
         "ODEProblem(MatrixOperator(A, update_func), u0, tspan)", call. = FALSE)
  L
}

# the matrix of the operator at (u, p, t) (update_coefficients! then convert)
.lin_A <- function(L, u, p, t) .sop_matrix(update_coefficients(L, u, p, t), length(u))

# exponential!(X, ExpMethodGeneric()) * v
.lin_expmul <- function(X, v) .as_num(.eu_exp_generic(X) %*% v)

.lie_alg <- function(name, order, core, krylov = FALSE, m = 30L, iop = 0L, adaptive = FALSE, ...) {
  krylov <- isTRUE(krylov)
  m <- as.integer(m)
  iop <- as.integer(iop)
  # expv(t, A, v; m = min(alg.m, size(L, 1)), alg.iop)
  kexpv <- function(t, A, v) .eu_expv(t, A, v, m = min(m, length(v)), iop = iop)
  step <- function(ctx, t, u, dt, cache) {
    L <- .lin_operator(ctx, name)
    r <- core(ctx, L, t, .as_num(u), dt, kexpv)
    fl <- ctx$f(r$u, ctx$p, t + dt)
    res <- list(u = r$u, k = list(cache$fsal, fl), cache = list(fsal = fl), nf = 1L)
    res$err <- r$err
    res
  }
  .new_alg(name, order, adaptive = adaptive, step = step, interp = .hermite_interp,
           krylov = krylov, m = m, iop = iop, ...)
}

#' Lie group and Magnus methods for linear problems
#'
#' Methods of OrdinaryDiffEqLinear for `u' = A(u, t) u`, where the problem
#' is `ODEProblem(L, u0, tspan)` with a SciML operator `L`, usually
#' `MatrixOperator(A, update_func)` whose `update_func(A, u, p, t)` returns
#' the matrix at `(u, p, t)`. All are fixed step methods (pass `dt`) except
#' `MagnusAdapt4()`, which is adaptive.
#'
#' `MagnusMidpoint()` and `MagnusLeapfrog()` are the second order Magnus
#' midpoint and leapfrog methods, `LieEuler()` the first order Lie Euler
#' method, `MagnusGauss4()`, `MagnusGL4()`, `MagnusGL6()`, `MagnusGL8()`,
#' `MagnusNC6()` and `MagnusNC8()` Magnus methods of orders 4 to 8 with
#' Gauss, Gauss-Legendre and Newton-Cotes quadratures, `RKMK2()` and
#' `RKMK4()` the Runge-Kutta-Munthe-Kaas methods of orders 2 and 4,
#' `LieRK4()` the fourth order Lie Runge-Kutta method, `CG2()`, `CG3()` and
#' `CG4a()` the Crouch-Grossman methods of orders 2, 3 and 4, and
#' `MagnusAdapt4()` the adaptive fourth order Magnus method of Li and Deng.
#' `CayleyEuler()` is the first order Cayley method for matrix states `U' =
#' A(U, t) U` (an `ODEProblem` of the operator, or a `SplitODEProblem` whose
#' first function is the operator), whose step is `U V'` multiplied by `V =
#' (I - dt A / 2)^-1 (I + dt A / 2)`; the state is the matrix in column
#' major order. `LinearExponential()` is the exact solution `exp(dt A) u` of
#' a linear problem with a constant operator (one step over the interval
#' when no `dt` is given).
#'
#' With `krylov = TRUE` the exponential acting on the state is approximated
#' in a Krylov subspace of dimension `m` (Arnoldi, or Lanczos for symmetric
#' matrices); `iop` sets the length of the incomplete orthogonalization.
#' `LinearExponential(krylov = "simple")` uses a fixed subspace and
#' `krylov = "adaptive"` the adaptive time stepping of Niesen and Wright
#' with tolerance `reltol`.
#'
#' @param krylov Whether (or for `LinearExponential()` how: `"off"`,
#'   `"simple"` or `"adaptive"`) to use Krylov approximations.
#' @param m Krylov subspace dimension.
#' @param iop Incomplete orthogonalization length (0 for full Arnoldi).
#' @return An algorithm object.
#' @examples
#' L <- MatrixOperator(diag(2), function(A, u, p, t) matrix(c(0, 1, -1, -0.1 * t), 2))
#' prob <- ODEProblem(L, c(1, 0), c(0, 1))
#' solve(prob, MagnusGL4(), dt = 0.1)$u[[11]]
#' solve(ODEProblem(MatrixOperator(matrix(c(-1, 1, 0, -2), 2)), c(1, 1), c(0, 1)),
#'       LinearExponential())$u[[2]]
#' @name linear_lie_methods
NULL

#' @rdname linear_lie_methods
#' @export
MagnusMidpoint <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("MagnusMidpoint", 2, function(ctx, L, t, u, dt, kexpv) {
    A <- .lin_A(L, u, ctx$p, t + dt / 2)
    list(u = if (krylov) kexpv(dt, A, u) else .lin_expmul(dt * A, u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
MagnusLeapfrog <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("MagnusLeapfrog", 2, function(ctx, L, t, u, dt, kexpv) {
    if (ctx$iter == 1L) {
      A <- .lin_A(L, u, ctx$p, t + dt / 2)
      list(u = if (krylov) kexpv(dt, A, u) else .lin_expmul(dt * A, u))
    } else {
      # u_{n+1} = exp(2 dt A(u_n, t_n)) u_{n-1}
      A <- .lin_A(L, u, ctx$p, t)
      u2 <- .as_num(ctx$uprev2)
      list(u = if (krylov) kexpv(2 * dt, A, u2) else .lin_expmul(2 * dt * A, u2))
    }
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
LieEuler <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("LieEuler", 1, function(ctx, L, t, u, dt, kexpv) {
    A <- .lin_A(L, u, ctx$p, t)
    list(u = if (krylov) kexpv(dt, A, u) else .lin_expmul(dt * A, u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
MagnusGauss4 <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("MagnusGauss4", 4, function(ctx, L, t, u, dt, kexpv) {
    A <- .lin_A(L, u, ctx$p, t + dt * (1 / 2 + sqrt(3) / 6))
    B <- .lin_A(L, u, ctx$p, t + dt * (1 / 2 - sqrt(3) / 6))
    C <- B %*% A - A %*% B
    list(u = if (krylov) kexpv(dt, (A + B) / 2 + (dt * sqrt(3)) * C / 12, u)
             else .lin_expmul((dt / 2) * (A + B) + ((dt * dt) * (sqrt(3) / 12)) * C, u))
  }, krylov, m, iop)

# the commutator X Y - Y X
.lin_comm <- function(X, Y) X %*% Y - Y %*% X

# the high order part of MagnusGL8 and MagnusNC8 from B0, ..., B3
.lin_magnus8 <- function(B0, B1, B2, B3, dt) {
  Q1 <- .lin_comm(-38 * B0 / 5 + 24 * B2, B3)
  Q2 <- .lin_comm(63 * B0 / 5 - 84 * B2, -5 * B1 / 28 + B3)
  Q3 <- .lin_comm(19 * B0 / 28 - 15 * B2 / 7, .lin_comm(B0, B2 + dt * (61 * Q1 / 588 - Q2 / 12)))
  Q4 <- .lin_comm(B3, 20 * Q1 / 7 + 10 * Q2)
  Q5 <- .lin_comm(-6025 * B0 / 4116 + 2875 * B2 / 343, .lin_comm(B2, Q1))
  Q6 <- .lin_comm(B3, 20 * (Q3 + Q4) / 7 + 820 * dt * Q5 / 189)
  Q7 <- (-1 / 42) * .lin_comm(B0, .lin_comm(B0, Q3 - Q4 / 3 + dt * Q5))
  O1 <- dt * B0
  O2 <- (dt * dt) * (Q1 + Q2)
  O3 <- (dt * dt * dt) * (Q3 + Q4) + (dt^4) * (Q5 + Q6) + (dt^5) * Q7
  O1 + O2 + O3
}

# the high order part of MagnusGL6 and MagnusNC6 from B0, B1, B2
.lin_magnus6 <- function(B0, B1, B2, dt) {
  O1 <- dt * B0
  O2 <- (dt * dt) * .lin_comm(B1, 3 * B0 / 2 - 6 * B2)
  O34 <- (dt * dt) * .lin_comm(B0, .lin_comm(B0, dt * B2 / 2 - O2 / 60)) + (3 * dt / 5) * .lin_comm(B1, O2)
  O1 + O2 + O34
}

#' @rdname linear_lie_methods
#' @export
MagnusNC6 <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("MagnusNC6", 6, function(ctx, L, t, u, dt, kexpv) {
    p <- ctx$p
    A0 <- .lin_A(L, u, p, t)
    A1 <- .lin_A(L, u, p, t + dt / 4)
    A2 <- .lin_A(L, u, p, t + dt / 2)
    A3 <- .lin_A(L, u, p, t + 3 * dt / 4)
    A4 <- .lin_A(L, u, p, t + dt)
    B0 <- (1 / 90) * (7 * (A0 + A4) + 32 * (A1 + A3) + 12 * (A2))
    B1 <- (1 / 90) * ((7 / 2) * (A4 - A0) + 8 * (A3 - A1))
    B2 <- (1 / 90) * ((7 / 4) * (A0 + A4) + 2 * (A1 + A3))
    O <- .lin_magnus6(B0, B1, B2, dt)
    list(u = if (krylov) kexpv(1, O, u) else .lin_expmul(O, u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
MagnusGL6 <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("MagnusGL6", 6, function(ctx, L, t, u, dt, kexpv) {
    p <- ctx$p
    A1 <- .lin_A(L, u, p, t - dt * (sqrt(3 / 20)) + dt / 2)
    A2 <- .lin_A(L, u, p, t + dt / 2)
    A3 <- .lin_A(L, u, p, t + dt * (sqrt(3 / 20)) + dt / 2)
    B0 <- (1 / 18) * (5 * (A1 + A3) + 8 * A2)
    B1 <- (sqrt(15) / 36) * (A3 - A1)
    B2 <- (1 / 24) * (A1 + A3)
    O <- .lin_magnus6(B0, B1, B2, dt)
    list(u = if (krylov) kexpv(1, O, u) else .lin_expmul(O, u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
MagnusGL8 <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("MagnusGL8", 8, function(ctx, L, t, u, dt, kexpv) {
    p <- ctx$p
    v1 <- (1 / 2) * sqrt((3 + 2 * sqrt(6 / 5)) / 7)
    v2 <- (1 / 2) * sqrt((3 - 2 * sqrt(6 / 5)) / 7)
    A1 <- .lin_A(L, u, p, t - dt * v1 + dt / 2)
    A4 <- .lin_A(L, u, p, t + dt * v1 + dt / 2)
    A2 <- .lin_A(L, u, p, t - dt * v2 + dt / 2)
    A3 <- .lin_A(L, u, p, t + dt * v2 + dt / 2)
    w1 <- (1 / 2) - (1 / 6) * sqrt(5 / 6)
    w2 <- (1 / 2) + (1 / 6) * sqrt(5 / 6)
    S1 <- A1 + A4
    S2 <- A2 + A3
    R1 <- A4 - A1
    R2 <- A3 - A2
    B0 <- (1 / 2) * (w1 * S1 + w2 * S2)
    B2 <- (1 / 2) * ((v1 * v1) * w1 * S1 + (v2 * v2) * w2 * S2)
    B1 <- (1 / 2) * (v1 * w1 * R1 + v2 * w2 * R2)
    B3 <- (1 / 2) * ((v1 * v1 * v1) * w1 * R1 + (v2 * v2 * v2) * w2 * R2)
    O <- .lin_magnus8(B0, B1, B2, B3, dt)
    list(u = if (krylov) kexpv(1, O, u) else .lin_expmul(O, u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
MagnusNC8 <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("MagnusNC8", 8, function(ctx, L, t, u, dt, kexpv) {
    p <- ctx$p
    A0 <- .lin_A(L, u, p, t)
    A1 <- .lin_A(L, u, p, t + dt / 6)
    A2 <- .lin_A(L, u, p, t + 2 * dt / 6)
    A3 <- .lin_A(L, u, p, t + 3 * dt / 6)
    A4 <- .lin_A(L, u, p, t + 4 * dt / 6)
    A5 <- .lin_A(L, u, p, t + 5 * dt / 6)
    A6 <- .lin_A(L, u, p, t + dt)
    S1 <- A0 + A6
    S2 <- A1 + A5
    S3 <- A2 + A4
    S4 <- A3
    R1 <- A6 - A0
    R2 <- A5 - A1
    R3 <- A4 - A2
    B0 <- (1 / 840) * (41 * S1 + 216 * S2 + 27 * S3 + 272 * S4)
    B2 <- (1 / 840) * (41 * S1 / 4 + 216 * S2 / 9 + 27 * S3 / 36)
    B1 <- (1 / 840) * (41 * R1 / 2 + 216 * R2 / 3 + 27 * R3 / 6)
    B3 <- (1 / 840) * (41 * R1 / 8 + 216 * R2 / 27 + 27 * R3 / 216)
    O <- .lin_magnus8(B0, B1, B2, B3, dt)
    list(u = if (krylov) kexpv(1, O, u) else .lin_expmul(O, u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
MagnusGL4 <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("MagnusGL4", 4, function(ctx, L, t, u, dt, kexpv) {
    A1 <- .lin_A(L, u, ctx$p, t + dt * (1 / 2 - sqrt(3) / 6))
    A2 <- .lin_A(L, u, ctx$p, t + dt * (1 / 2 + sqrt(3) / 6))
    O <- (dt / 2) * (A1 + A2) - (dt * dt) * (sqrt(3) / 12) * .lin_comm(A1, A2)
    list(u = if (krylov) kexpv(1, O, u) else .lin_expmul(O, u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
RKMK2 <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("RKMK2", 2, function(ctx, L, t, u, dt, kexpv) {
    k1 <- dt * .lin_A(L, u, ctx$p, t)
    k2 <- dt * .lin_A(L, .lin_expmul(k1, u), ctx$p, t)
    list(u = if (krylov) kexpv(1 / 2, k1 + k2, u) else .lin_expmul((1 / 2) * (k1 + k2), u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
RKMK4 <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("RKMK4", 4, function(ctx, L, t, u, dt, kexpv) {
    p <- ctx$p
    k1 <- dt * .lin_A(L, u, p, t)
    k2 <- dt * .lin_A(L, .lin_expmul(k1 / 2, u), p, t)
    k3 <- dt * .lin_A(L, .lin_expmul(k1 / 2 - .lin_comm(k1, k2) / 8, u), p, t)
    k4 <- dt * .lin_A(L, .lin_expmul(k3, u), p, t)
    O <- k1 + 2 * k2 + 2 * k3 + k4 - .lin_comm(k1, k4) / 2
    list(u = if (krylov) kexpv(1 / 6, O, u) else .lin_expmul((1 / 6) * O, u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
LieRK4 <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("LieRK4", 4, function(ctx, L, t, u, dt, kexpv) {
    p <- ctx$p
    k1 <- dt * .lin_A(L, u, p, t)
    Y2 <- .lin_expmul(k1 / 2, u)
    k2 <- dt * .lin_A(L, Y2, p, t)
    Y3 <- .lin_expmul(k2 / 2, u)
    k3 <- dt * .lin_A(L, Y3, p, t)
    Y4 <- .lin_expmul(k3 - k1 / 2, Y2)
    k4 <- dt * .lin_A(L, Y4, p, t)
    y12 <- .lin_expmul((3 * k1 + 2 * k2 + 2 * k3 - k4) / 12, u)
    O <- -k1 + 2 * k2 + 2 * k3 + 3 * k4
    list(u = if (krylov) kexpv(1 / 12, O, y12) else .lin_expmul((1 / 12) * O, y12))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
CG2 <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("CG2", 2, function(ctx, L, t, u, dt, kexpv) {
    k1 <- dt * .lin_A(L, u, ctx$p, t)
    k2 <- dt * .lin_A(L, .lin_expmul(k1, u), ctx$p, t)
    list(u = .as_num((.eu_exp_generic((1 / 2) * (k1)) %*% .eu_exp_generic((1 / 2) * k2)) %*% u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
CG3 <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("CG3", 3, function(ctx, L, t, u, dt, kexpv) {
    p <- ctx$p
    E <- .eu_exp_generic
    A <- .lin_A(L, u, p, t)
    v2 <- .lin_expmul((3 / 4) * dt * A, u)
    B <- .lin_A(L, v2, p, t + (3 * dt / 4))
    v3 <- .as_num(E((119 / 216) * dt * B) %*% (E((17 / 108) * dt * A) %*% u))
    C <- .lin_A(L, v3, p, t + (17 * dt / 24))
    list(u = .as_num(((E(dt * (24 / 17) * C) %*% E(dt * (-2 / 3) * B)) %*% E(dt * (13 / 51) * A)) %*% u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
CG4a <- function(krylov = FALSE, m = 30, iop = 0)
  .lie_alg("CG4a", 4, function(ctx, L, t, u, dt, kexpv) {
    p <- ctx$p
    E <- .eu_exp_generic
    A <- .lin_A(L, u, p, t)
    v2 <- .lin_expmul((0.8177227988124852) * dt * A, u)
    B <- .lin_A(L, v2, p, t + (0.8177227988124852 * dt))
    v3 <- .as_num(E((0.3199876375476427) * dt * B) %*% (E((0.0659864263556022) * dt * A) %*% u))
    C <- .lin_A(L, v3, p, t + (0.3859740639032449 * dt))
    v4 <- .as_num(E((0.9214417194464946) * dt * C) %*% (E((0.4997857776773573) * dt * B) %*%
                                                           (E((-1.0969984448371582) * dt * A) %*% u)))
    D <- .lin_A(L, v4, p, t + (0.3242290522866937 * dt))
    v5 <- .as_num((((E((0.3552358559023322) * dt * D) %*% E((0.2390958372307326) * dt * C)) %*%
                        E((1.3918565724203246) * dt * B)) %*% E((-1.1092979392113465) * dt * A)) %*% u)
    F <- .lin_A(L, v5, p, t + (0.8768903263420429 * dt))
    M <- (((E(dt * (0.3322195591068374) * F) %*% E(dt * (-0.1907142565505889) * D)) %*%
             E(dt * (0.739781398537078) * C)) %*% E(dt * (-0.018369853156402) * B)) %*%
      E(dt * (0.1370831520630755) * A)
    list(u = .as_num(M %*% u))
  }, krylov, m, iop)

#' @rdname linear_lie_methods
#' @export
MagnusAdapt4 <- function()
  .lie_alg("MagnusAdapt4", 4, function(ctx, L, t, u, dt, kexpv) {
    p <- ctx$p
    A0 <- .lin_A(L, u, p, t)
    k1 <- dt * A0
    Q1 <- k1
    y2 <- (1 / 2) * Q1
    A1 <- .lin_A(L, .lin_expmul(y2, u), p, t + dt / 2)
    k2 <- dt * A1
    Q2 <- k2 - k1
    y3 <- (1 / 2) * Q1 + (1 / 4) * Q2
    A2 <- .lin_A(L, .lin_expmul(y3, u), p, t + dt / 2)
    k3 <- dt * A2
    Q3 <- k3 - k2
    y4 <- Q1 + Q2
    A3 <- .lin_A(L, .lin_expmul(y4, u), p, t + dt)
    k4 <- dt * A3
    Q4 <- k4 - 2 * k2 + k1
    y5 <- (1 / 2) * Q1 + (1 / 4) * Q2 + (1 / 3) * Q3 - (1 / 24) * Q4 - (1 / 48) * .lin_comm(Q1, Q2)
    A4 <- .lin_A(L, .lin_expmul(y5, u), p, t + dt / 2)
    k5 <- dt * A4
    Q5 <- k5 - k2
    y6 <- Q1 + Q2 + (2 / 3) * Q3 + (1 / 6) * Q4 - (1 / 6) * .lin_comm(Q1, Q2)
    A5 <- .lin_A(L, .lin_expmul(y6, u), p, t + dt)
    k6 <- dt * A5
    Q6 <- k6 - 2 * k2 + k1
    X <- Q2 - Q3 + Q5 + (1 / 2) * Q6
    v4 <- Q1 + Q2 + (2 / 3) * Q5 + (1 / 6) * Q6 - (1 / 6) * (Q1 %*% X - X %*% Q1)
    unew <- .lin_expmul(v4, u)
    res <- list(u = unew)
    if (ctx$adaptive) res$err <- unew - .lin_expmul(y6, u)
    res
  }, adaptive = TRUE)

#' @rdname linear_lie_methods
#' @export
CayleyEuler <- function()
  .lie_alg("CayleyEuler", 2, function(ctx, L, t, u, dt, kexpv) {
    n <- L$size[1L]
    if (is.na(n) || length(u) != n * n)
      stop("CayleyEuler() needs a matrix state with as many rows and columns as the operator", call. = FALSE)
    U <- matrix(u, n, n)
    X <- .sop_matrix(update_coefficients(L, U, ctx$p, t), n) * dt
    I <- diag(n)
    # cay(A) = inv(I - A / 2) (I + A / 2); the step is V (U V') (in place)
    V <- solve(I - (1 / 2) * X) %*% (I + (1 / 2) * X)
    list(u = .as_num(V %*% (U %*% t(V))))
  })

#' @rdname linear_lie_methods
#' @export
LinearExponential <- function(krylov = c("off", "simple", "adaptive"), m = 10, iop = 0) {
  krylov <- match.arg(krylov)
  m <- as.integer(m)
  iop <- as.integer(iop)
  alg <- .lie_alg("LinearExponential", 1, function(ctx, L, t, u, dt, kexpv) {
    # A = convert(AbstractMatrix, f.f): the operator as left by the last
    # derivative evaluation, at the state and time of the step start
    A <- .lin_A(L, u, ctx$p, t)
    n <- length(u)
    list(u = switch(krylov,
      off = .lin_expmul(dt * A, u),
      simple = .eu_expv_ks(dt, .eu_arnoldi(A, u, m = min(m, n), iop = iop)),
      adaptive = .eu_expv_timestep(dt, A, u, adaptive = TRUE, m = min(m, n), iop = iop, tol = ctx$reltol)))
  })
  alg$krylov <- krylov
  alg$m <- m
  alg$iop <- iop
  alg$dt_required <- FALSE
  alg
}
