# Adjoints of optimization problems, after SciMLSensitivity.jl v7.119
# src/optimization_adjoint.jl (OptimizationAdjoint: implicit
# differentiation of the KKT conditions with active set detection) and
# src/concrete_solve.jl (UnconstrainedOptimizationAdjoint: the steady state
# adjoint of the stationarity condition grad f(u*, p) = 0).
#
# Julia takes the objective gradient, the Hessians and the parameter
# derivatives of the KKT residual from the OptimizationFunction and forward
# mode automatic differentiation. Here a problem built from a System gets
# them symbolically (the gradient of the cost and of each constraint, and
# their partial derivatives by the parameters); a problem built from R
# functions uses its `grad`, `hess`, `cons_j` and `cons_h` when given and
# central differences otherwise.

#' Adjoints of optimization problems
#'
#' `UnconstrainedOptimizationAdjoint()` differentiates the solution `u*` of
#' an unconstrained [OptimizationProblem()] with respect to the parameters
#' by the implicit function theorem on the optimality condition
#' `grad f(u*, p) = 0`: with the Hessian `H`, the adjoint `lambda = H^-1 g_u`
#' gives `dG/dp = g_p - lambda' d(grad f)/dp`, the [SteadyStateAdjoint()] of
#' the stationarity condition.
#'
#' `OptimizationAdjoint()` differentiates a constrained problem through its
#' KKT conditions. Equality constraints (`lcons == ucons`), inequality
#' constraints `lcons <= cons(u, p) <= ucons` within `active_tol` of a bound
#' and bounds `lb <= u <= ub` within `active_tol` are active; the
#' multipliers are the least squares solution of the stationarity
#' condition, and constraints whose multiplier has the wrong sign are
#' dropped from the active set. One solve of the symmetric KKT system
#' `[L_xx J'; J 0] lambda = [g_u; 0]` then gives `dG/dp = -lambda' dF/dp`
#' for the KKT residual `F = (grad_x L, active constraints)`.
#'
#' Both are chosen automatically by [solve_gradient()] for an
#' `OptimizationProblem` (the KKT adjoint when the problem has constraints or
#' bounds), and [adjoint_sensitivities()] applies them to an optimization
#' solution with `dgdu`. The derivatives come from the System the problem
#' was built from (symbolic), or from the `grad`, `hess`, `cons_j` and
#' `cons_h` of the `OptimizationFunction` and central differences. Julia
#' requires the stored gradient of the `OptimizationFunction`; here a
#' missing gradient is taken by central differences, the counterpart of its
#' automatic differentiation backend.
#'
#' @param chunk_size,autodiff,diff_type,objective_ad Options of the Julia
#'   signatures (the automatic differentiation choices).
#' @param autojacvec `FALSE` takes central differences also for problems
#'   built from a System.
#' @param linsolve,linsolve_kwargs Linear solver options (unused; the KKT
#'   system is solved directly).
#' @param active_tol Tolerance of the active set detection;
#'   `sqrt(.Machine$double.eps)` when `NULL`.
#' @return A sensitivity algorithm.
#' @examples
#' f <- function(x, p) (x[1] - p[1])^2 + p[2] * (x[2] - x[1]^2)^2
#' prob <- OptimizationProblem(f, c(0, 0), c(1.5, 3))
#' # u* = (p1, p1^2): the gradient of u1 + u2 is (1 + 2 p1, 0)
#' solve_gradient(function(sol) sum(sol$u), prob, BFGS(),
#'                sensealg = UnconstrainedOptimizationAdjoint())$p
#' @name optimization_adjoints
NULL

#' @rdname optimization_adjoints
#' @export
UnconstrainedOptimizationAdjoint <- function(chunk_size = 0, autodiff = TRUE, diff_type = "central",
                                             objective_ad = TRUE, autojacvec = NULL, linsolve = NULL,
                                             linsolve_kwargs = list()) {
  .sensealg("UnconstrainedOptimizationAdjoint", "AdjointSensitivityAlgorithm", autojacvec = autojacvec,
            linsolve = linsolve, linsolve_kwargs = if (length(linsolve_kwargs)) linsolve_kwargs,
            objective_ad = objective_ad)
}

#' @rdname optimization_adjoints
#' @export
OptimizationAdjoint <- function(chunk_size = 0, autodiff = NULL, diff_type = "forward", autojacvec = NULL,
                                linsolve = NULL, linsolve_kwargs = list(), active_tol = NULL) {
  .sensealg("OptimizationAdjoint", "AdjointSensitivityAlgorithm", autojacvec = autojacvec,
            linsolve = linsolve, linsolve_kwargs = if (length(linsolve_kwargs)) linsolve_kwargs,
            active_tol = active_tol, autodiff = autodiff)
}

# The derivatives of an optimization problem at (u, p): grad and Hessian of
# the objective, constraint values, Jacobian and Hessians, and the
# parameter derivatives gq = d(grad f)/dp, cq = dc/dp and cjq[[i]] =
# d(grad c_i)/dp.
.sens_opt_fns <- function(prob, sensealg) {
  of <- prob$f
  sys <- prob$sys
  m <- length(prob$lcons)
  pv <- .sens_pvec
  out <- list(m = m)
  out$cons <- if (m) function(u, p) as.numeric(of$cons(u, p))
  if (!is.null(sys) && inherits(sys, "System") && !identical(sensealg$autojacvec, FALSE)) {
    us <- unknowns(sys)
    ps <- .problem_params(sys)
    args <- .opt_args(sys)
    mat <- function(exprs) {
      if (!length(ps)) return(function(u, p) matrix(0, length(exprs), 0L))
      .build_closure(.partials(exprs, ps, us), args, out = "matrix")
    }
    gexpr <- .opt_gradient(sys, cost(sys))
    grad <- .build_closure(gexpr, args, out = "vector")
    out$grad <- function(u, p) as.numeric(grad(u, p))
    out$hess <- generate_cost_hessian(sys)
    out$gq <- mat(gexpr)
    if (m) {
      res <- canonical_constraints(sys)
      cj <- generate_constraint_jacobian(sys)
      out$cons_j <- function(u, p) matrix(cj(u, p), m, length(u))
      out$cons_h <- generate_constraint_hessian(sys)
      out$cq <- mat(.opt_full(sys, res))
      cjq <- lapply(res, function(r) mat(.opt_gradient(sys, r)))
      out$cjq <- function(u, p) lapply(cjq, function(f) f(u, p))
    }
    return(out)
  }
  out$grad <- if (is.function(of$grad)) function(u, p) as.numeric(of$grad(u, p)) else
    function(u, p) .sens_fd_grad(function(z) of$f(z, p), u)
  out$hess <- if (is.function(of$hess)) function(u, p) as.matrix(of$hess(u, p)) else
    function(u, p) .sens_fd_jac(function(z) out$grad(z, p), u)
  out$gq <- function(u, p) .sens_fd_jac(function(q) out$grad(u, pv(p, q)), as.numeric(p))
  if (m) {
    out$cons_j <- if (is.function(of$cons_j)) function(u, p) matrix(of$cons_j(u, p), m, length(u)) else
      function(u, p) .sens_fd_jac(function(z) out$cons(z, p), u)
    out$cons_h <- if (is.function(of$cons_h)) function(u, p) lapply(of$cons_h(u, p), as.matrix) else
      function(u, p) lapply(seq_len(m), function(i) .sens_fd_jac(function(z) out$cons_j(z, p)[i, ], u))
    out$cq <- function(u, p) .sens_fd_jac(function(q) out$cons(u, pv(p, q)), as.numeric(p))
    out$cjq <- function(u, p) lapply(seq_len(m), function(i)
      .sens_fd_jac(function(q) out$cons_j(u, pv(p, q))[i, ], as.numeric(p)))
  }
  out
}

# The least squares solution of A x = b by QR (LinearSolve's
# QRFactorization); coefficients of dependent columns are zero.
.sens_opt_lsq <- function(A, b) {
  x <- qr.coef(qr(A), b)
  x[is.na(x)] <- 0
  as.numeric(x)
}

# dG/dp of the solution x of an optimization problem for the cotangent du,
# by the stationarity condition (UnconstrainedOptimizationAdjoint) or the
# KKT system (OptimizationAdjointSensitivityFunction and
# OptimizationAdjointProblem).
.sens_opt_dp <- function(prob, x, p, du, sensealg) {
  fns <- .sens_opt_fns(prob, sensealg)
  n <- length(x)
  np <- length(p)
  if (!np)
    stop("Your model does not have parameters, and thus it is impossible to calculate the derivative of ",
         "the solution with respect to the parameters.", call. = FALSE)
  if (inherits(sensealg, "UnconstrainedOptimizationAdjoint"))
    return(.sens_ift_pullback(list(jac = fns$hess, pjac = fns$gq), x, p, du, numeric(np)))
  lcons <- prob$lcons
  ucons <- prob$ucons
  m <- fns$m
  atol <- sensealg$active_tol %||% sqrt(.Machine$double.eps)
  eq <- which(abs(lcons - ucons) <= 2 * atol)
  ineq <- which(abs(lcons - ucons) > 2 * atol)
  cval <- if (m) fns$cons(x, p) else numeric()
  alb <- ineq[which(abs(cval[ineq] - lcons[ineq]) <= atol)]
  aub <- ineq[which(abs(cval[ineq] - ucons[ineq]) <= atol)]
  alb_v <- if (!is.null(prob$lb)) which(abs(x - prob$lb) <= atol) else integer()
  aub_v <- if (!is.null(prob$ub)) which(abs(x - prob$ub) <= atol) else integer()
  gradf <- fns$grad(x, p)
  J <- if (m) fns$cons_j(x, p) else matrix(0, 0L, n)
  Jg <- J[eq, , drop = FALSE]
  neq <- length(eq)
  mtol <- sqrt(.Machine$double.eps)
  max_iters <- length(alb) + length(aub) + length(alb_v) + length(aub_v) + 1L
  negatives <- TRUE
  iter <- 0L
  while (negatives && iter < max_iters) {
    iter <- iter + 1L
    nact <- length(alb) + length(aub)
    nb <- length(alb_v) + length(aub_v)
    Jb <- matrix(0, nb, n)
    Jb[cbind(seq_along(alb_v), alb_v)] <- -1
    Jb[cbind(length(alb_v) + seq_along(aub_v), aub_v)] <- 1
    Jh <- rbind(-J[alb, , drop = FALSE], J[aub, , drop = FALSE], Jb)
    ntot <- nact + nb
    duals <- if (neq + ntot == 0L) numeric() else .sens_opt_lsq(t(rbind(Jg, Jh)), -gradf)
    ystar <- duals[seq_len(neq)]
    zI <- duals[neq + seq_len(nact)]
    zb <- duals[neq + nact + seq_len(nb)]
    negatives <- any(zI < -mtol) || any(zb < -mtol)
    if (negatives) {
      nlb <- length(alb)
      nlbv <- length(alb_v)
      keep_ub <- zI[nlb + seq_along(aub)] >= -mtol
      alb <- alb[zI[seq_len(nlb)] >= -mtol]
      aub <- aub[keep_ub]
      keep_ubv <- zb[nlbv + seq_along(aub_v)] >= -mtol
      alb_v <- alb_v[zb[seq_len(nlbv)] >= -mtol]
      aub_v <- aub_v[keep_ubv]
    }
  }
  mu <- numeric(m)
  mu[eq] <- ystar
  mu[alb] <- mu[alb] - zI[seq_along(alb)]
  mu[aub] <- mu[aub] + zI[length(alb) + seq_along(aub)]
  # the Lagrangian Hessian
  Lxx <- fns$hess(x, p)
  if (m) {
    Hc <- fns$cons_h(x, p)
    for (i in which(mu != 0)) Lxx <- Lxx + mu[i] * Hc[[i]]
  }
  N <- n + neq + ntot
  K <- matrix(0, N, N)
  K[seq_len(n), seq_len(n)] <- Lxx
  if (neq) {
    K[seq_len(n), n + seq_len(neq)] <- t(Jg)
    K[n + seq_len(neq), seq_len(n)] <- Jg
  }
  if (ntot) {
    K[seq_len(n), n + neq + seq_len(ntot)] <- t(Jh)
    K[n + neq + seq_len(ntot), seq_len(n)] <- Jh
  }
  lam <- solve(K, c(du, numeric(neq + ntot)))
  lx <- lam[seq_len(n)]
  ly <- lam[n + seq_len(neq)]
  lz <- lam[n + neq + seq_len(nact)]
  # dG/dp = -d/dq (lx' grad_x L + ly' g + lz' h_I); the bound rows do not
  # depend on the parameters
  phi <- as.numeric(crossprod(fns$gq(x, p), lx))
  if (m) {
    cjq <- fns$cjq(x, p)
    for (i in which(mu != 0)) phi <- phi + mu[i] * as.numeric(crossprod(cjq[[i]], lx))
    cq <- fns$cq(x, p)
    if (neq) phi <- phi + as.numeric(crossprod(cq[eq, , drop = FALSE], ly))
    if (nact) phi <- phi + as.numeric(crossprod(rbind(-cq[alb, , drop = FALSE], cq[aub, , drop = FALSE]), lz))
  }
  -phi
}

# adjoint_sensitivities() of an optimization solution.
.sens_opt_adjoint <- function(sol, sensealg, dgdu = NULL, dgdp = NULL, g = NULL) {
  env <- .sol_env(sol)
  prob <- env$prob
  if (is.null(dgdu) && is.null(g)) stop("dgdu must be specified for OptimizationAdjoint", call. = FALSE)
  x <- as.numeric(env$U[1L, ])
  p <- prob$p
  du <- if (!is.null(dgdu)) as.numeric(dgdu(x, p)) else .sens_fd_grad(function(z) g(z, p), x)
  dp <- .sens_opt_dp(prob, x, p, du, sensealg)
  if (!is.null(dgdp)) dp <- dp + as.numeric(dgdp(x, p)) else if (!is.null(g))
    dp <- dp + .sens_fd_grad(function(q) g(x, .sens_pvec(p, q)), as.numeric(p))
  list(du0 = stats::setNames(numeric(length(x)), names(prob$u0)), dp = stats::setNames(dp, names(prob$p)))
}

# solve_gradient() of an optimization problem: the optimum does not depend
# on the starting point.
.sens_opt_gradient <- function(loss, prob, alg, dots, sensealg) {
  if (!inherits(sensealg, c("OptimizationAdjoint", "UnconstrainedOptimizationAdjoint")))
    stop("optimization problems are differentiated with OptimizationAdjoint() or ",
         "UnconstrainedOptimizationAdjoint()", call. = FALSE)
  sol <- do.call(solve, c(list(prob, alg), dots))
  ct <- .sens_loss_cotangent(loss, sol)
  x <- as.numeric(.sol_env(sol)$U[1L, ])
  dp <- .sens_opt_dp(prob, x, prob$p, as.numeric(ct$dU[1L, ]), sensealg) + ct$dp
  list(u0 = stats::setNames(numeric(length(x)), names(prob$u0)), p = stats::setNames(dp, names(prob$p)),
       value = ct$value, sol = sol)
}
