# Rosenbrock methods with mass matrices, after OrdinaryDiffEq.jl
# (OrdinaryDiffEqRosenbrock 2.7.5, OrdinaryDiffEqDifferentiation 3.12.5):
#   lib/OrdinaryDiffEqRosenbrock/src/rosenbrock_perform_step.jl
#     perform_step! for Rosenbrock23ConstantCache, Rosenbrock32ConstantCache
#     and the unified Rodas loop (RosenbrockCombinedConstantCache)
#   lib/OrdinaryDiffEqRosenbrock/src/rosenbrock_interpolants.jl
#     Shampine's interpolant of Rosenbrock23/32, the H matrix dense output
#     and the Hermite fallback of the methods without an H matrix
#   lib/OrdinaryDiffEqRosenbrock/src/rosenbrock_caches.jl, alg_utils.jl,
#     algorithms.jl  tableau of each method, orders, W-method and FSAL traits
#   lib/OrdinaryDiffEqDifferentiation/src/derivative_utils.jl
#     _rosenbrock_jac_reuse_decision, calc_rosenbrock_differentiation (the
#     out of place path), calc_W (W = J - M / (dt gamma)), calc_tderivative
#
# The functions take u by value, as the out of place (ConstantCache) path of
# OrdinaryDiffEq, which is the reference for the steps, the counts and the
# rounding. Each attempt evaluates J afresh (counted twice: calc_W and the
# copy kept for reuse), except that the W-methods reuse J along accepted
# steps as CVODE does (max_jac_age, jac_reuse_gamma_tol). A singular W
# rejects the attempt (EEst = 2). The tableaus are in
# R/solvers-rosenbrock-tableaus.R (tools/rosenbrock-dump.jl).

# ---------------------------------------------------------------------------
# Jacobian, time derivative and W

# JacReuseState of a fresh cache (_make_jac_reuse_state)
.rosen_new_cache <- function() {
  cache <- new.env(parent = emptyenv())
  cache$last_dtgamma <- 0
  cache$pending_dtgamma <- 0
  cache$last_naccept <- 0L
  cache$last_step_iter <- 0L
  cache$last_u_length <- 0L
  cache$J <- NULL
  cache$dT <- NULL
  cache
}

# _rosenbrock_jac_reuse_decision: TRUE for a fresh Jacobian. The number of
# accepted steps is ctx$success_iter (the steps accepted since the cache was
# made, as the cache starts with last_naccept = 0).
.rosen_reuse_decision <- function(ctx, cache, spec, dtgamma, n) {
  if (!spec$wmethod) return(TRUE)
  if (ctx$iter <= 1L) return(TRUE)
  # a cache made again after set_u! (R restarts the caches) has no Jacobian
  if (is.null(cache$J)) return(TRUE)
  op <- ctx$fobj$operator
  if (inherits(op, "AbstractSciMLOperator") && islinear(op)) return(FALSE)
  if (!ctx$adaptive) return(TRUE)
  naccept <- ctx$success_iter
  if (!is.null(ctx$mass)) {
    if (naccept > cache$last_naccept) {
      cache$last_naccept <- naccept
      return(TRUE)
    }
    return(FALSE)
  }
  if (isTRUE(ctx$composite)) return(TRUE)
  if (naccept > cache$last_naccept) cache$last_dtgamma <- cache$pending_dtgamma
  if (cache$last_dtgamma == 0) return(TRUE)
  if (cache$last_step_iter != 0L && ctx$iter > cache$last_step_iter + 1L) return(TRUE)
  if (isTRUE(ctx$derivative_discontinuity)) return(TRUE)
  if (n != cache$last_u_length && cache$last_u_length != 0L) return(TRUE)
  if (ctx$EEst_prev > 1) return(TRUE)
  if (abs(dtgamma / cache$last_dtgamma - 1) > spec$jac_reuse_gamma_tol) return(TRUE)
  if (naccept - cache$last_naccept >= spec$max_jac_age) {
    cache$last_naccept <- naccept
    return(TRUE)
  }
  FALSE
}

# calc_tderivative: prob$f$tgrad, else finite differences in t (forward in
# the direction of diffdir, or central for AutoFiniteDiff("central")),
# counted as one evaluation. f0 is f(u, p, t), returned for the step to use
# as its first stage value.
.rosen_tgrad <- function(ctx, u, t) {
  if (!is.null(ctx$tgrad)) return(list(dT = .as_num(ctx$tgrad(u, ctx$p, t)), nf = 0L, f0 = NULL))
  f0 <- ctx$f(u, ctx$p, t)
  ad <- ctx$autodiff
  if (inherits(ad, "AutoSparse")) ad <- ad$dense_ad
  # AutoForwardDiff: ForwardDiff.derivative of f in t with a dual number,
  # or finite differences when f cannot run on duals
  if (inherits(ad, "AutoForwardDiff") && !inherits(ad, "AutoFiniteDiff")) {
    jr <- .ad_try_jacobian(function(tt) ctx$f(u, ctx$p, tt), t)
    if (!is.null(jr)) return(list(dT = as.numeric(jr$jacobian), nf = 1L, f0 = f0))
  }
  dT <- if (inherits(ad, "AutoFiniteDiff") && ad$fdtype == "central") {
    h <- .Machine$double.eps^(1 / 3) * max(1, abs(t))
    (ctx$f(u, ctx$p, t + h) - ctx$f(u, ctx$p, t - h)) / (2 * h)
  } else {
    h <- sqrt(.Machine$double.eps) * max(1, abs(t)) * .rosen_diffdir(ctx, u, t)
    (ctx$f(u, ctx$p, t + h) - f0) / h
  }
  list(dT = dT, nf = 1L, f0 = f0)
}

# diffdir (OrdinaryDiffEqCore misc_utils.jl): the side of t the difference
# uses, so that the stencil stays inside the time span near its end. Needs
# ctx$tspan; without it the difference is forward in t.
.rosen_diffdir <- function(ctx, u, t) {
  ts <- ctx$tspan
  if (is.null(ts)) return(1)
  difference <- max(abs(u)) * sqrt(.Machine$double.eps)
  if (ts[2L] > ts[1L]) {
    if (t > ts[2L] - difference) -1 else 1
  } else if (t < ts[2L] + difference) 1 else -1
}

# calc_W: W = J - M / (dt gamma), with M * inv(dtgamma) for a mass matrix.
.rosen_W <- function(ctx, J, dtgamma) {
  M <- ctx$mass
  if (is.null(M)) {
    W <- J
    diag(W) <- diag(W) - 1 / dtgamma
  } else W <- J - M * (1 / dtgamma)
  .rosen_lu(W)
}

# DiffEqBase.default_factorize, lu(W; check = false), and W \ b as Julia
# runs them through OpenBLAS: up to 16 equations OpenBLAS factorizes with its
# left looking dgetf2 (fused dot and gemv kernels) and solves with dgetrs
# (fused axpy), which .rosen_getf2 and .rosen_getrs reproduce operation by
# operation; larger systems use LAPACK through Matrix::lu (blocked OpenBLAS
# rounds differently). NULL when W is not finite or a pivot is exactly zero
# (issuccess_W).
.rosen_lu <- function(W) {
  W <- as.matrix(W)
  if (!all(is.finite(W))) return(NULL)
  n <- nrow(W)
  if (n > 16L) {
    F <- Matrix::lu(W, warnSing = FALSE)
    return(if (any(F@x[seq.int(1L, n * n, by = n + 1L)] == 0)) NULL else F)
  }
  F <- .rosen_getf2(W)
  if (any(diag(F$LU) == 0)) NULL else F
}

.rosen_lusolve <- function(F, b) {
  if (is.list(F)) .rosen_getrs(F, b) else .as_num(Matrix::solve(F, b))
}

# OpenBLAS lapack/getf2/getf2.c: column j gets the earlier row swaps, the
# unit lower solve by dot products (fused, accumulated from zero), the gemv
# update of the rows below (fused, column by column), then the pivot search,
# the row swap of columns 1 to j and the scaling by the reciprocal pivot.
.rosen_getf2 <- function(A) {
  n <- nrow(A)
  ipiv <- integer(n)
  for (j in seq_len(n)) {
    b <- A[, j]
    for (i in seq_len(j - 1L)) {
      jp <- ipiv[i]
      if (jp != i) b[c(i, jp)] <- b[c(jp, i)]
    }
    if (j > 2L) for (i in 2:(j - 1L)) {
      s <- 0
      for (k in seq_len(i - 1L)) s <- .fma(A[i, k], b[k], s)
      b[i] <- b[i] - s
    }
    rows <- j:n
    for (k in seq_len(j - 1L)) b[rows] <- .fma(A[rows, k], -b[k], b[rows])
    jp <- j - 1L + which.max(abs(b[rows]))
    ipiv[j] <- jp
    A[, j] <- b
    if (b[jp] != 0) {
      if (jp != j) A[c(j, jp), seq_len(j)] <- A[c(jp, j), seq_len(j)]
      if (j < n) A[(j + 1L):n, j] <- A[(j + 1L):n, j] * (1 / A[j, j])
    }
  }
  list(LU = A, ipiv = ipiv)
}

# OpenBLAS dgetrs for one right hand side: the row swaps, then trsv_NLU and
# trsv_NUN (division by the pivot, fused axpy updates).
.rosen_getrs <- function(F, b) {
  A <- F$LU
  n <- nrow(A)
  for (i in seq_len(n)) {
    jp <- F$ipiv[i]
    if (jp != i) b[c(i, jp)] <- b[c(jp, i)]
  }
  for (i in seq_len(n - 1L)) {
    r <- (i + 1L):n
    b[r] <- .fma(-b[i], A[r, i], b[r])
  }
  for (i in n:1) {
    b[i] <- b[i] / A[i, i]
    if (i > 1L) {
      r <- seq_len(i - 1L)
      b[r] <- .fma(-b[i], A[r, i], b[r])
    }
  }
  b
}

# calc_rosenbrock_differentiation (out of place): dT and the factorized W
# (NULL when singular) with the counts. A fresh Jacobian is evaluated by
# calc_W and again for the reuse cache, so njacs and its evaluations count
# twice; a reused one only builds W.
.rosen_differentiation <- function(ctx, t, u, dtgamma, cache, spec) {
  new_jac <- .rosen_reuse_decision(ctx, cache, spec, dtgamma, length(u))
  cache$last_step_iter <- ctx$iter
  out <- list(nf = 0L, njacs = 0L, nw = 1L, f0 = NULL, eigen_est = NULL)
  if (new_jac) {
    tg <- .rosen_tgrad(ctx, u, t)
    jj <- .calc_J(ctx, u, t)
    cache$J <- jj$J
    cache$dT <- tg$dT
    cache$pending_dtgamma <- dtgamma
    cache$last_u_length <- length(u)
    out$nf <- tg$nf + 2L * jj$nf
    out$njacs <- 2L
    out$f0 <- tg$f0
    if (isTRUE(ctx$composite)) out$eigen_est <- max(rowSums(abs(jj$J)))
  }
  out$W <- .rosen_W(ctx, cache$J, dtgamma)
  out$dT <- cache$dT
  out
}

# ---------------------------------------------------------------------------
# Error norms

# calculate_residuals (fused denominator) and ODE_DEFAULT_NORM
.rosen_resid_norm <- function(ctx, err, u0, u1) {
  .rosen_ode_norm(err / .fma(pmax(abs(u0), abs(u1)), ctx$reltol, ctx$abstol))
}

# DiffEqBase.ODE_DEFAULT_NORM for arrays: its @fastmath loop accumulates
# the squares with fused multiply adds, then sqrt(x / n). This is Julia's
# result bit for bit up to 8 components; LLVM vectorizes longer loops in an
# order not reproduced here, which changes the last bits.
.rosen_ode_norm <- function(x) {
  if (is.complex(x)) return(.rms_norm(x))
  s <- 0
  for (v in x) s <- .fma(v, v, s)
  sqrt(s / max(length(x), 1L))
}

# LinearAlgebra.norm (generic_norm2, vectors shorter than 32)
.rosen_norm2 <- function(x) {
  m <- max(abs(x))
  if (m == 0 || !is.finite(m)) return(m)
  if (is.finite(length(x) * m * m) && m * m != 0) return(sqrt(sum(x * x)))
  m * sqrt(sum((x / m)^2))
}

# get_differential_vars for the Hermite mask and the algebraic residual of
# Rosenbrock23/32: NULL (all differential), a logical vector, or FALSE
# (DifferentialVarsUndefined: no variable is treated as differential).
.rosen_dvars <- function(M) {
  if (is.null(M)) return(NULL)
  if (all(M != 0)) return(rep(TRUE, nrow(M)))
  if (all(M[row(M) != col(M)] == 0)) return(diag(M) != 0)
  FALSE
}

# ---------------------------------------------------------------------------
# Dense output

# The H matrix interpolant (interp_order 2, 3 or 4, fused as the @..
# broadcast under @muladd) and the Hermite fallback (interp_order -1) of
# the methods without an H matrix, linear in the algebraic variables.
.rosen_interp <- function(order) {
  force(order)
  function(theta, dt, y0, y1, k) {
    t1 <- 1 - theta
    if (order < 0L) {
      m <- attr(k, "dvars")
      if (is.null(m)) return(.hermite_interp_fused(theta, dt, y0, y1, k))
      lin <- .fma(theta, y1, t1 * y0)
      if (identical(m, FALSE)) return(lin)
      return(ifelse(m, .hermite_interp_fused(theta, dt, y0, y1, k), lin))
    }
    inner <- k[[order]]
    for (j in rev(seq_len(order - 1L))) inner <- .fma(theta, inner, k[[j]])
    .fma(theta, .fma(t1, inner, y1), t1 * y0)
  }
}

# Shampine's interpolant (MATLAB ode23s), rosenbrock2332pre0
.rosen_interp23 <- function(d) {
  force(d)
  function(theta, dt, y0, y1, k) {
    c1 <- theta * (1 - theta) / (1 - 2 * d)
    c2 <- theta * (theta - 2 * d) / (1 - 2 * d)
    y0 + dt * (c1 * k[[1L]] + c2 * k[[2L]])
  }
}

# ---------------------------------------------------------------------------
# Steps

# The unified Rodas loop (perform_step! of RosenbrockCombinedConstantCache).
.rosen_rodas_step <- function(spec) {
  force(spec)
  tab <- spec$tab
  A <- tab$A; C <- tab$C; cc <- tab$c; d <- tab$d; H <- tab$H; b <- tab$b; bt <- tab$btilde
  s <- nrow(A)
  nH <- nrow(H)
  # stages that repeat the abscissa and the row of the previous one reuse
  # its f value
  same <- vapply(seq_len(s), function(st) st > 2L && cc[st] == cc[st - 1L] &&
                   all(A[st, seq_len(st - 1L)] == A[st - 1L, seq_len(st - 1L)]), TRUE)
  function(ctx, t, u, dt, cache) {
    f <- ctx$f; p <- ctx$p; M <- ctx$mass
    n <- length(u)
    dtC <- C / dt
    dtd <- dt * d
    dtgamma <- dt * tab$gamma
    nf_reset <- .rosen_reset_fsal(ctx, t, u, cache)
    dw <- .rosen_differentiation(ctx, t, u, dtgamma, cache, spec)
    out <- list(nf = dw$nf + nf_reset, njacs = dw$njacs, nw = dw$nw, nsolve = 0L, eigen_est = dw$eigen_est,
                cache = cache)
    W <- dw$W
    if (is.null(W))
      return(c(out, list(u = u, EEst = 2, k = rep(list(numeric(n)), max(nH, 2L)))))
    dT <- dw$dT
    du <- dw$f0 %||% f(u, p, t)
    out$nf <- out$nf + 1L
    fsalfirst <- du
    ks <- vector("list", s)
    ks[[1L]] <- .rosen_lusolve(W, -.fma(dtd[1L], dT, du))
    for (st in 2:s) {
      ui <- u
      for (i in seq_len(st - 1L)) ui <- .fma(A[st, i], ks[[i]], ui)
      if (!same[st]) {
        du <- f(ui, p, .fma(cc[st], dt, t))
        out$nf <- out$nf + 1L
      }
      lt <- numeric(n)
      for (i in seq_len(st - 1L)) lt <- .fma(dtC[st, i], ks[[i]], lt)
      if (!is.null(M)) lt <- .mass_mul(M, lt)
      ks[[st]] <- .rosen_lusolve(W, -.fma(dtd[st], dT, du + lt))
      out$nsolve <- out$nsolve + 1L
    }
    unew <- u
    for (i in seq_len(s)) if (b[i] != 0) unew <- .fma(b[i], ks[[i]], unew)
    EEst <- 0
    if (ctx$adaptive && !is.null(bt)) {
      err <- numeric(n)
      for (i in seq_len(s)) if (bt[i] != 0) err <- .fma(bt[i], ks[[i]], err)
      EEst <- .rosen_resid_norm(ctx, err, u, unew)
    }
    k <- NULL
    if (!isFALSE(ctx$calck)) {
      if (nH > 0L) {
        k <- lapply(seq_len(nH), function(j) {
          kj <- numeric(n)
          for (i in seq_len(s)) kj <- .fma(H[j, i], ks[[i]], kj)
          kj
        })
        if (spec$residual_control && ctx$adaptive && isTRUE(EEst < 1)) {
          # Rodas5Pr: the residual of the interpolant at the midpoint
          k2 <- 0.5 * ((u + unew) + 0.5 * (k[[1L]] + 0.5 * (k[[2L]] + 0.5 * k[[3L]])))
          du1 <- ((0.25 * (k[[2L]] + k[[3L]]) - u) + unew) / dt
          dum <- f(k2, p, t + dt / 2)
          out$nf <- out$nf + 1L
          du2 <- if (is.null(M)) du1 - dum else .mass_mul(M, du1) - dum
          EEst <- max(.rosen_norm2(du2) / .rosen_norm2(.fma(ctx$reltol, k2, ctx$abstol)), EEst)
        }
      } else {
        # no H matrix: f at both ends for the Hermite interpolant
        k <- list(fsalfirst, f(unew, p, t + dt))
        out$nf <- out$nf + 1L
        attr(k, "dvars") <- .rosen_dvars(M)
      }
    }
    c(out, list(u = unew, EEst = EEst, k = k))
  }
}

# The algebraic residual added to EEst on mass matrix problems by
# Rosenbrock23/32: ifelse(differential_vars, 0, f(u)) / abstol.
.rosen_alg_residual <- function(ctx, fu) {
  dv <- .rosen_dvars(ctx$mass)
  .rosen_ode_norm(ifelse(dv, 0, fu) / ctx$abstol)
}

# perform_step! of Rosenbrock23ConstantCache (r32 = FALSE) and
# Rosenbrock32ConstantCache (r32 = TRUE). Upstream Rosenbrock32 keeps
# f(uprev + dt k2), the derivative at the second order point it needs for
# k3, as the FSAL derivative of the next step, although the step ends at the
# third order solution (rosenbrock_perform_step.jl, both caches). Here the
# next step evaluates f at the new u instead (deliberate difference).
.rosen_r2332_step <- function(spec, r32) {
  force(spec); force(r32)
  c32 <- spec$tab$c32; d <- spec$tab$d
  function(ctx, t, u, dt, cache) {
    f <- ctx$f; p <- ctx$p; M <- ctx$mass
    dtgamma <- dt * d
    neginvdtgamma <- -(1 / dtgamma)
    dto2 <- dt / 2
    dto6 <- dt / 6
    nf_reset <- .rosen_reset_fsal(ctx, t, u, cache)
    dw <- .rosen_differentiation(ctx, t, u, dtgamma, cache, spec)
    out <- list(nf = dw$nf + nf_reset, njacs = dw$njacs, nw = dw$nw, nsolve = 0L, eigen_est = dw$eigen_est,
                cache = cache)
    W <- dw$W
    if (is.null(W)) {
      # a fixed step is accepted all the same
      if (!ctx$adaptive) cache$reeval_fsal <- TRUE
      return(c(out, list(u = u, EEst = 2, k = list(numeric(length(u)), numeric(length(u))))))
    }
    dT <- dw$dT
    fsalfirst <- cache$fsal
    k1 <- if (r32) .rosen_lusolve(W, -(fsalfirst + dtgamma * dT)) / dtgamma else
      .rosen_lusolve(W, fsalfirst + dtgamma * dT) * neginvdtgamma
    f1 <- f(.fma(dto2, k1, u), p, t + dto2)
    k2 <- .fma(.rosen_lusolve(W, f1 - .mass_mul(M, k1)), neginvdtgamma, k1)
    out$nf <- out$nf + 1L
    out$nsolve <- 2L
    lin_tmp <- function(fsallast) {
      if (is.null(M)) return(.fma(dt, dT, .fma(-2, k1 - fsalfirst, .fma(-c32, k2 - f1, fsallast))))
      lt <- .mass_mul(M, .fma(2, k1, c32 * k2))
      .fma(dt, dT, .fma(2, fsalfirst, .fma(c32, f1, fsallast - lt)))
    }
    EEst <- 0
    fsallast <- NULL
    if (r32) {
      fsallast <- f(.fma(dt, k2, u), p, t + dt)
      k3 <- .rosen_lusolve(W, lin_tmp(fsallast)) * neginvdtgamma
      unew <- .fma(dto6, .fma(4, k2, k1 + k3), u)
      out$nf <- out$nf + 1L
      out$nsolve <- 3L
      if (ctx$adaptive) EEst <- .rosen_resid_norm(ctx, dto6 * ((k1 - 2 * k2) + k3), u, unew)
    } else {
      unew <- u + dt * k2
      if (ctx$adaptive) {
        fsallast <- f(unew, p, t + dt)
        k3 <- .rosen_lusolve(W, lin_tmp(fsallast)) * neginvdtgamma
        out$nf <- out$nf + 1L
        out$nsolve <- 3L
        EEst <- .rosen_resid_norm(ctx, .mass_mul(M, dto6 * ((k1 - 2 * k2) + k3)), u, unew)
      }
    }
    if (ctx$adaptive && !is.null(M)) EEst <- EEst + .rosen_alg_residual(ctx, fsallast)
    # update_fsal!: only_diagonal_mass_matrix methods evaluate f again
    # after every fixed step (reset_fsal!, at the start of the next step);
    # so does Rosenbrock32 after every accepted step (see above)
    if (!ctx$adaptive || (r32 && isTRUE(EEst <= 1))) cache$reeval_fsal <- TRUE else
      if (isTRUE(EEst <= 1)) cache$fsal <- fsallast
    c(out, list(u = unew, EEst = EEst, k = list(k1, k2)))
  }
}

# ---------------------------------------------------------------------------
# Caches and algorithm objects

# initialize!: Rosenbrock23/32 evaluate the FSAL derivative (counted); the
# methods that are FSAL in OrdinaryDiffEqCore's sense keep f(u) for the step
# limiter and the callbacks (reset_fsal!) but the unified loop evaluates
# nothing at initialization.
.rosen_init <- function(spec) {
  force(spec)
  function(ctx, t, u) {
    M <- ctx$mass
    if (spec$kind == "r2332" && !is.null(M) && any(M[row(M) != col(M)] != 0))
      stop(spec$name, " only works with diagonal mass matrices. Please choose a solver suitable for your ",
           "problem (e.g. Rodas5P)", call. = FALSE)
    if (identical(spec$name, "Rosenbrock32") && !is.null(M) && any(colSums(M != 0) == 0))
      stop("Rosenbrock32 is unstable on DAEs (singular mass matrix). Use Rosenbrock23 or Rodas5P instead.",
           call. = FALSE)
    cache <- .rosen_new_cache()
    cache$nf0 <- 0L
    if (spec$fsal) {
      # the step limiter of the driver asks for reset_fsal! at the next step
      cache$defer_fsal <- TRUE
      cache$reeval_fsal <- FALSE
      # the RodasTableau methods do not count it: it is the f0 of the
      # initial step estimate
      cache$fsal <- ctx$f(u, ctx$p, t)
      if (spec$kind == "r2332") cache$nf0 <- 1L
    }
    cache
  }
}

# reset_fsal! of update_fsal! (in loopheader! of OrdinaryDiffEqCore): after a
# step limiter, a callback or a fixed Rosenbrock23/32 step the FSAL methods
# evaluate f again at the start of the next step (so never after the last
# one). The RodasTableau methods evaluate f(uprev) in every step and only
# count this evaluation.
.rosen_reset_fsal <- function(ctx, t, u, cache) {
  if (!isTRUE(cache$reeval_fsal)) return(0L)
  cache$reeval_fsal <- FALSE
  cache$fsal <- ctx$f(u, ctx$p, t)
  1L
}

# A callback changed u: the cache and its Jacobian state carry over the
# event and the FSAL methods evaluate f again at the next step.
.rosen_on_event <- function(ctx, t, u, cache) {
  if (!is.null(cache$fsal)) cache$reeval_fsal <- TRUE
  list(cache = cache, nf = 0L)
}

.rosen_check_limiter <- function(x, arg) {
  if (!is.null(x) && !is.function(x)) stop("`", arg, "` must be a function or NULL", call. = FALSE)
  x
}

# A Rosenbrock algorithm object. kind is "rodas" (RodasTableau methods) or
# "r2332" (Rosenbrock23/32); tab_name names the tableau in .ROSEN_TABLEAUS.
.rosen_alg <- function(name, order, tab_name, kind = "rodas", wmethod = FALSE, fsal = TRUE, adaptive = TRUE,
                       autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                       stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03,
                       residual_control = FALSE, adaptive_order = NULL) {
  tab <- .ROSEN_TABLEAUS[[tab_name]]
  spec <- list(name = name, kind = kind, tab = tab, wmethod = wmethod, fsal = fsal,
               max_jac_age = as.integer(max_jac_age), jac_reuse_gamma_tol = .as_num(jac_reuse_gamma_tol),
               residual_control = residual_control)
  if (kind == "r2332") {
    step <- .rosen_r2332_step(spec, r32 = name == "Rosenbrock32")
    interp <- .rosen_interp23(tab$d)
  } else {
    step <- .rosen_rodas_step(spec)
    nH <- nrow(tab$H)
    interp <- .rosen_interp(if (nH == 0L) -1L else if (name %in% c("Rodas3P", "Rodas23W")) 2L else nH)
  }
  args <- list(name, order, adaptive = adaptive, step = step, interp = interp, init = .rosen_init(spec),
               stiff = TRUE, mass_matrix = TRUE, on_event = .rosen_on_event,
               step_limiter = .rosen_check_limiter(step_limiter, "step_limiter"),
               stage_limiter = .rosen_check_limiter(stage_limiter, "stage_limiter"),
               concrete_jac = concrete_jac, linsolve = linsolve, max_jac_age = spec$max_jac_age,
               jac_reuse_gamma_tol = spec$jac_reuse_gamma_tol, wmethod = wmethod)
  if (adaptive) args$qsteady_max <- 6 / 5
  if (!is.null(adaptive_order)) args$adaptive_order <- adaptive_order
  .with_autodiff(autodiff, do.call(.new_alg, args))
}

#' Rosenbrock algorithms
#'
#' Linearly implicit Rosenbrock methods for stiff problems and for
#' differential algebraic equations in mass matrix form (index 1, singular
#' mass matrices as produced by [mtkcompile()]), as in OrdinaryDiffEq.jl.
#' `Rosenbrock23()` is the order 2(3) W-method of MATLAB's ode23s and
#' `Rosenbrock32()` its order 3(2) counterpart, both with Shampine's
#' interpolant; `Rodas3()`, `Rodas3d()` and `ROS3P()` are third order,
#' `Rodas3P()` and the W-method `Rodas23W()` third order with a stiff aware
#' interpolant, `Rodas4()`, `Rodas42()`, `Rodas4P()`, `Rodas4P2()` and the
#' W-method `Rodas4PW()` fourth order, `Rodas5()`, `Rodas5P()`, `Rodas5Pe()`
#' (another embedded method) and `Rodas5Pr()` (with a residual test of the
#' interpolant) fifth order and `Rodas6P()` sixth order stiffly accurate
#' methods with their stiff aware dense output. Each step uses the Jacobian
#' `prob$f$jac` when present and differentiation otherwise, and
#' `prob$f$tgrad` for non autonomous problems when present. The W-methods
#' keep the Jacobian along accepted steps until it is `max_jac_age` steps
#' old or the step size changes by more than `jac_reuse_gamma_tol`.
#' After an accepted step `Rosenbrock32()` evaluates the derivative at the
#' new solution for the next step; OrdinaryDiffEq.jl reuses the derivative
#' at the second order point there, which costs one evaluation less per step.
#'
#' @param autodiff How Jacobians are computed when the problem has no
#'   symbolic Jacobian: `AutoForwardDiff()` (dual numbers, as ForwardDiff;
#'   central differences when the function cannot run on them),
#'   `AutoFiniteDiff(fdtype)` or `AutoSparse(ad)`; see [AutoForwardDiff()].
#'   The default uses forward differences.
#' @param concrete_jac,linsolve Accepted for compatibility with
#'   OrdinaryDiffEq.jl; the linear systems are always solved densely.
#' @param step_limiter A function `limiter(u, integrator, p, t)` returning
#'   the state after each accepted step, or `NULL`.
#' @param stage_limiter Accepted for compatibility: as in OrdinaryDiffEq.jl,
#'   the out of place steps do not apply a stage limiter.
#' @param max_jac_age Number of accepted steps after which a W-method
#'   evaluates the Jacobian again (1 for `Rosenbrock23()` and
#'   `Rosenbrock32()`, which then evaluate it at every step).
#' @param jac_reuse_gamma_tol Relative change of `dt * gamma` above which a
#'   W-method evaluates the Jacobian again.
#' @return An algorithm object.
#' @examples
#' vdp <- function(u, p, t) c(u[2], p * (1 - u[1]^2) * u[2] - u[1])
#' prob <- ODEProblem(vdp, c(2, 0), c(0, 1), p = 1000)
#' sol <- solve(prob, Rodas5P())
#' sol$retcode
#' solve(prob, Rosenbrock23())$stats$njacs
#' @name rosenbrock_algorithms
NULL

# The constructors pass their arguments (an environment) to .rosen_alg.
.rosen_from <- function(a, name, order, tab_name, ...)
  .rosen_alg(name, order, tab_name, autodiff = a$autodiff, concrete_jac = a$concrete_jac, linsolve = a$linsolve,
             step_limiter = a$step_limiter, stage_limiter = a$stage_limiter, max_jac_age = a$max_jac_age,
             jac_reuse_gamma_tol = a$jac_reuse_gamma_tol, ...)

#' @rdname rosenbrock_algorithms
#' @export
Rosenbrock23 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                         stage_limiter = NULL, max_jac_age = 1L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rosenbrock23", 2L, "Rosenbrock23Tableau", kind = "r2332", wmethod = TRUE, adaptive_order = 3L)

#' @rdname rosenbrock_algorithms
#' @export
Rosenbrock32 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                         stage_limiter = NULL, max_jac_age = 1L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rosenbrock32", 3L, "Rosenbrock32Tableau", kind = "r2332", wmethod = TRUE, adaptive_order = 2L)

#' @rdname rosenbrock_algorithms
#' @export
ROS3P <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                  stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS3P", 3L, "ROS3PRodasTableau", fsal = TRUE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas3 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                   stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas3", 3L, "Rodas3RodasTableau", fsal = TRUE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas3d <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                    stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas3d", 3L, "Rodas3dRodasTableau", fsal = TRUE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas23W <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                     stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas23W", 3L, "Rodas23WRodasTableau", wmethod = TRUE, fsal = FALSE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas3P <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                    stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas3P", 3L, "Rodas3PRodasTableau", fsal = FALSE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas4 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                   stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas4", 4L, "Rodas4Tableau", fsal = FALSE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas42 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                    stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas42", 4L, "Rodas42Tableau", fsal = FALSE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas4P <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                    stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas4P", 4L, "Rodas4PTableau", fsal = FALSE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas4P2 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                     stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas4P2", 4L, "Rodas4P2Tableau", fsal = FALSE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas4PW <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                     stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas4PW", 4L, "Rodas4PWTableau", wmethod = TRUE, fsal = FALSE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas5 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                   stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas5", 5L, "Rodas5Tableau", fsal = FALSE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas5P <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                    stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas5P", 5L, "Rodas5PTableau", fsal = FALSE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas5Pe <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                     stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas5Pe", 5L, "Rodas5PeTableau", fsal = FALSE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas5Pr <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                     stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas5Pr", 5L, "Rodas5PTableau", fsal = FALSE, residual_control = TRUE)

#' @rdname rosenbrock_algorithms
#' @export
Rodas6P <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, step_limiter = NULL,
                    stage_limiter = NULL, max_jac_age = 20L, jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Rodas6P", 6L, "Rodas6PTableau", fsal = FALSE)
