# HybridExplicitImplicitRK and Tsit5DA of OrdinaryDiffEqRosenbrock 2.7.5
# (algorithms.jl; perform_step! of HybridExplicitImplicitConstantCache in
# rosenbrock_perform_step.jl; Tsit5DATableau in rosenbrock_tableaus.jl):
# explicit Runge-Kutta stages for the differential variables and linearly
# implicit stages for the algebraic variables of a semi explicit index 1
# DAE (diagonal mass matrix), so only the algebraic block of the Jacobian is
# factorized. Without a mass matrix the method is the explicit Runge-Kutta
# method of the tableau. The Jacobian of the algebraic block is evaluated by
# differentiation at each attempt, as upstream not counted in the
# statistics.

# Matrix times vector as OpenBLAS gemv: one fused multiply add per column.
.rosen_matvec <- function(A, x) {
  acc <- numeric(nrow(A))
  for (j in seq_len(ncol(A))) acc <- .fma(A[, j], x[j], acc)
  acc
}

# DI.jacobian(uf, autodiff, uprev): finite differences for AutoFiniteDiff
# (also with an analytic Jacobian, which upstream ignores here), else the
# analytic Jacobian or the differentiation of R/solvers-implicit.R.
.rosen_hybrid_jac <- function(ctx, u, t) {
  ad <- ctx$autodiff
  if (inherits(ad, "AutoSparse")) ad <- ad$dense_ad
  if (inherits(ad, "AutoFiniteDiff")) {
    fn <- function(x) ctx$f(x, ctx$p, t)
    return(if (ad$fdtype == "forward") .jac_fd(fn, u) else .jac_central(fn, u))
  }
  .calc_J(ctx, u, t)$J
}

.rosen_hybrid_step <- function(tab) {
  force(tab)
  A <- tab$A; C <- tab$C; b <- tab$b; bhat <- tab$bhat; cc <- tab$c; d <- tab$d; H <- tab$H
  s <- nrow(A)
  function(ctx, t, u, dt, cache) {
    f <- ctx$f; p <- ctx$p; M <- ctx$mass
    n <- length(u)
    nf <- 0L
    has_alg <- FALSE
    if (!is.null(M)) {
      dg <- diag(M)
      dv <- which(dg != 0)
      av <- which(dg == 0)
      has_alg <- length(av) > 0L
    }
    if (has_alg) {
      tg <- .rosen_tgrad(ctx, u, t)
      nf <- nf + tg$nf
      dT <- tg$dT
      J <- .rosen_hybrid_jac(ctx, u, t)
      g_z <- J[av, av, drop = FALSE]
      g_y <- J[av, dv, drop = FALSE]
      Wz <- .rosen_lu(-tab$gamma * g_z)
      if (is.null(Wz)) stop("the algebraic block of the Jacobian is singular at t = ", t, call. = FALSE)
    }
    du <- f(u, p, t)
    nf <- nf + 1L
    stage_k <- function(du, st, ks) {
      if (is.null(M)) return(dt * du)
      kn <- numeric(n)
      kn[dv] <- dt * du[dv]
      if (has_alg) {
        if (st == 1L) {
          rhs <- .fma(dt * d[1L], dT[av], du[av] + .rosen_matvec(g_y, C[1L, 1L] * kn[dv]))
        } else {
          gy <- C[st, st] * kn[dv]
          gz <- numeric(length(av))
          for (j in seq_len(st - 1L)) gy <- .fma(C[st, j], ks[[j]][dv], gy)
          for (j in seq_len(st - 1L)) gz <- .fma(C[st, j], ks[[j]][av], gz)
          rhs <- .fma(dt * d[st], dT[av], (du[av] + .rosen_matvec(g_y, gy)) + .rosen_matvec(g_z, gz))
        }
        kn[av] <- .rosen_lusolve(Wz, rhs)
      }
      kn
    }
    ks <- vector("list", s)
    ks[[1L]] <- stage_k(du, 1L, ks)
    for (st in 2:s) {
      us <- u
      for (j in seq_len(st - 1L)) us <- .fma(A[st, j], ks[[j]], us)
      du <- f(us, p, .fma(cc[st], dt, t))
      nf <- nf + 1L
      ks[[st]] <- stage_k(du, st, ks)
    }
    unew <- u
    for (i in seq_len(s)) unew <- .fma(b[i], ks[[i]], unew)
    EEst <- 0
    if (ctx$adaptive) {
      err <- numeric(n)
      for (i in seq_len(s)) err <- .fma(b[i] - bhat[i], ks[[i]], err)
      EEst <- .rosen_resid_norm(ctx, err, u, unew)
    }
    k <- NULL
    if (!isFALSE(ctx$calck)) {
      k <- lapply(seq_len(nrow(H)), function(j) {
        kj <- numeric(n)
        for (i in seq_len(s)) kj <- .fma(H[j, i], ks[[i]], kj)
        kj
      })
    }
    list(u = unew, EEst = EEst, k = k, cache = cache, nf = nf)
  }
}

# A tableau argument: a list with A, C, gamma, b, bhat, c, d and H, or a
# constructor function(T, T2) returning one.
.rosen_hybrid_tab <- function(tab) {
  if (is.function(tab)) tab <- tab("Float64", "Float64")
  need <- c("A", "C", "gamma", "b", "bhat", "c", "d", "H")
  if (!is.list(tab) || !all(need %in% names(tab)))
    stop("`tab` must be a tableau with entries ", paste(need, collapse = ", "), call. = FALSE)
  tab
}

#' Hybrid explicit and linearly implicit Runge-Kutta algorithms
#'
#' `HybridExplicitImplicitRK(tab, order)` advances the differential
#' variables of a semi explicit index 1 DAE in mass matrix form (a diagonal
#' mass matrix whose zero diagonal entries mark the algebraic variables)
#' with explicit Runge-Kutta stages and the algebraic variables with
#' Rosenbrock type linearly implicit stages, so only the algebraic block of
#' the Jacobian is factorized; without algebraic variables it is the
#' explicit method of the tableau. `Tsit5DA()` is the 12 stage order 5(4)
#' instance of Steinebach based on Tsit5, with a third order dense output.
#' As in OrdinaryDiffEq.jl the Jacobian of the algebraic block is evaluated
#' at every step (`max_jac_age` and `jac_reuse_gamma_tol` have no effect)
#' and is not counted in the statistics.
#'
#' @param tab The tableau: a list with entries `A`, `C`, `gamma`, `b`,
#'   `bhat`, `c`, `d` and `H` (the explicit coefficients, the linearly
#'   implicit Gamma coefficients, the weights of the solution and of the
#'   embedded method, the abscissae, the time derivative coefficients and
#'   the dense output matrix), or a function `(T, T2)` returning one.
#' @param order The order of the method, used by the step size controller.
#' @inheritParams rosenbrock_algorithms
#' @return An algorithm object.
#' @examples
#' f <- function(u, p, t) c(-u[1] + u[3], u[1] - 2 * u[2], u[3] - u[2] * u[2] - 0.5 * u[1])
#' M <- diag(c(1, 1, 0))
#' prob <- ODEProblem(ODEFunction(f, mass_matrix = M), c(1, 0.5, 0.75), c(0, 1))
#' solve(prob, Tsit5DA())(1)
#' @name hybrid_rosenbrock_algorithms
NULL

#' @rdname hybrid_rosenbrock_algorithms
#' @export
HybridExplicitImplicitRK <- function(tab, order, autodiff = NULL, concrete_jac = NULL, linsolve = NULL,
                                     step_limiter = NULL, stage_limiter = NULL, max_jac_age = 20L,
                                     jac_reuse_gamma_tol = 0.03) {
  tab <- .rosen_hybrid_tab(tab)
  .with_autodiff(autodiff, .new_alg("HybridExplicitImplicitRK", as.integer(order), step = .rosen_hybrid_step(tab),
                                    interp = .rosen_interp(nrow(tab$H)),
                                    init = function(ctx, t, u) list(nf0 = 0L), stiff = TRUE, mass_matrix = TRUE,
                                    qsteady_max = 6 / 5,
                                    step_limiter = .rosen_check_limiter(step_limiter, "step_limiter"),
                                    stage_limiter = .rosen_check_limiter(stage_limiter, "stage_limiter"),
                                    concrete_jac = concrete_jac, linsolve = linsolve,
                                    max_jac_age = as.integer(max_jac_age),
                                    jac_reuse_gamma_tol = .as_num(jac_reuse_gamma_tol), tab = tab))
}

#' @rdname hybrid_rosenbrock_algorithms
#' @export
Tsit5DA <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                    stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03) {
  alg <- HybridExplicitImplicitRK(.ROSEN_TABLEAUS$Tsit5DATableau, 5L, autodiff = autodiff,
                                  concrete_jac = concrete_jac, linsolve = linsolve, step_limiter = step_limiter,
                                  stage_limiter = stage_limiter, max_jac_age = max_jac_age,
                                  jac_reuse_gamma_tol = jac_reuse_gamma_tol)
  alg$name <- "Tsit5DA"
  class(alg) <- c("Tsit5DA", class(alg))
  alg
}
