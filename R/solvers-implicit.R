# Implicit ODE methods and their shared modified Newton solver, after
# OrdinaryDiffEq.jl:
#   lib/OrdinaryDiffEqNonlinearSolve/src/nlsolve.jl   nlsolve!, initial_η
#   lib/OrdinaryDiffEqNonlinearSolve/src/newton.jl    compute_step!, _compute_rhs
#     (COEFFICIENT_MULTISTEP form and the DAEFunction form)
#   lib/OrdinaryDiffEqNonlinearSolve/src/type.jl      NLNewton defaults
#   lib/OrdinaryDiffEqDifferentiation/src/derivative_utils.jl  do_newJW,
#     update_W!, calc_W, calc_J, calc_tderivative
#   lib/OrdinaryDiffEqSDIRK/src/generic_imex_perform_step.jl  ImplicitEuler
#     (:ie_dd2 error estimate), Trapezoid (:trap_dd3)
#   lib/OrdinaryDiffEqBDF/src/bdf_perform_step.jl, bdf_utils.jl,
#     controllers.jl  FBDF (and DFBDF, R/solvers-dae.R), with the time
#     filters (_fbdf_time_filter, backdiff!, _fbdf_filter_weight,
#     _filt_bdf_coefficients!, bdf3stab_coeff, _fbdf_select_filter!)
#   lib/OrdinaryDiffEqBDF/src/stald.jl  stability limit detection (CVODE's
#     cvSLdet): StabilityLimitDetectionState, stald_collect_data!,
#     stald_check!, stald_reset!, _stald_detect
# (OrdinaryDiffEqBDF 2.4.14, OrdinaryDiffEq.jl commit 4460ff8)
#
# The Newton state lives in an environment, so it survives rejected steps:
# Jacobian J (and dF/d(du) for DAEs), the factorized W, the step factor it
# was built for (W_γdt), the time of the Jacobian (J_t), status, failure
# count and the convergence rate estimates.

# ---------------------------------------------------------------------------
# Derivatives

# Jacobian of f(u, p, t) with respect to u: the analytic one when the
# problem has it (prob$f$jac), else forward differences with the
# FiniteDiff.jl step sqrt(eps) * max(1, |u_j|).
.jac_fd <- function(fn, u, f0 = fn(u)) {
  n <- length(u)
  J <- matrix(0, length(f0), n)
  for (j in seq_len(n)) {
    h <- sqrt(.Machine$double.eps) * max(1, abs(u[j]))
    uh <- u
    uh[j] <- u[j] + h
    J[, j] <- (fn(uh) - f0) / h
  }
  J
}

.calc_J <- function(ctx, u, t, f0 = NULL) {
  if (!is.null(ctx$jac)) return(list(J = as.matrix(ctx$jac(u, ctx$p, t)), nf = 0L))
  fn <- function(x) ctx$f(x, ctx$p, t)
  ad <- ctx$autodiff
  if (inherits(ad, "AutoSparse")) ad <- ad$dense_ad
  # AutoForwardDiff without a symbolic Jacobian, and AutoFiniteDiff("central"),
  # use central differences; the default and "forward" forward differences
  if (inherits(ad, "AutoForwardDiff") || (inherits(ad, "AutoFiniteDiff") && ad$fdtype != "forward")) {
    J <- .jac_central(fn, u)
    return(list(J = J, nf = 2L * length(u)))
  }
  extra <- is.null(f0)
  f0 <- f0 %||% fn(u)
  list(J = .jac_fd(fn, u, f0), nf = length(u) + extra)
}

.jac_central <- function(fn, u) {
  n <- length(u)
  cols <- lapply(seq_len(n), function(i) {
    h <- .Machine$double.eps^(1 / 3) * max(1, abs(u[i]))
    up <- u; um <- u
    up[i] <- up[i] + h; um[i] <- um[i] - h
    (fn(up) - fn(um)) / (2 * h)
  })
  matrix(unlist(cols), ncol = n)
}

# The `autodiff` choice of an implicit method (an ADType, see ?AutoForwardDiff)
.with_autodiff <- function(autodiff, alg) {
  if (!is.null(autodiff)) {
    if (!inherits(autodiff, "ADType"))
      stop("`autodiff` must be AutoForwardDiff(), AutoFiniteDiff() or AutoSparse()", call. = FALSE)
    alg$autodiff <- autodiff
  }
  alg
}

# Time derivative of f, calc_tderivative: prob$f$tgrad or forward differences.
.calc_tgrad <- function(ctx, u, t, f0) {
  if (!is.null(ctx$tgrad)) return(list(dT = as.numeric(ctx$tgrad(u, ctx$p, t)), nf = 0L))
  h <- sqrt(.Machine$double.eps) * max(1, abs(t))
  list(dT = (ctx$f(u, ctx$p, t + h) - f0) / h, nf = 1L)
}

# Factorization of W; NULL when W is not finite or singular (issuccess_W
# false). The LU of Julia (default_factorize, R/solvers-rosenbrock.R): on
# the badly conditioned W of a stiff DAE a QR solve differs from it far
# beyond rounding.
.factorize <- function(W) .rosen_lu(W)

.wsolve <- function(Wf, b) .rosen_lusolve(Wf, b)

.mass_mul <- function(M, x) if (is.null(M)) x else as.numeric(M %*% x)

# ---------------------------------------------------------------------------
# NLNewton (type.jl): κ = 1/100, max_iter = 10, fast_convergence_cutoff =
# 1/5, new_W_dt_cutoff = 1/5, check_div = TRUE.

.new_nlsolver <- function() {
  nl <- new.env(parent = emptyenv())
  nl$kappa <- 1 / 100
  nl$maxiters <- 10L
  nl$new_W_dt_cutoff <- 1 / 5
  nl$J <- NULL
  nl$Jdu <- NULL
  nl$W <- NULL
  nl$W_gdt <- 1
  nl$J_t <- NA_real_
  nl$status <- "Convergence"
  nl$nfails <- 0L
  nl$eta_old <- 1
  nl$prev_theta <- 1
  nl$firstcall <- TRUE
  nl
}

# do_newJW (derivative_utils.jl); returns c(new_jac, new_W). `gdt` is
# nlsolver.gamma * dt, against which the step change is measured.
.do_newJW <- function(nl, ctx, t, gdt, cj = NULL) {
  if (ctx$iter <= 1L || isTRUE(nl$firstcall)) return(c(TRUE, TRUE))
  if (!isTRUE(ctx$adaptive)) return(c(TRUE, TRUE))
  errorfail <- ctx$EEst_prev > 1
  fresh <- isTRUE(nl$J_t == t) && !isTRUE(ctx$derivative_discontinuity)
  if (fresh) {
    jbad <- FALSE
    small <- TRUE
  } else {
    small <- if (!is.null(cj)) abs(cj / nl$W_gdt - 1) <= nl$new_W_dt_cutoff else
      abs((1 / gdt) / (1 / nl$W_gdt) - 1) <= nl$new_W_dt_cutoff
    jbad <- nl$status == "TryAgain" && small
  }
  # stages after the first of a multistage method (isfirststage false) keep W
  # after a rejected step
  wbad <- !small || (!isFALSE(nl$isfs) && errorfail) || nl$status == "Divergence"
  c(jbad, jbad || wbad)
}

# update_W!: Jacobian at (uprev, t) and W = J - M / gdt for ODEs; for DAEs
# J_u and J_du (calc_J_dae: an analytic Jacobian at (duprev, uprev), finite
# differences at (du_pred, uprev)) and W = J_u + cj * J_du. `gdt_test` is
# the step factor do_newJW compares with (gamma * dt of the nonlinear
# solver; IRKC's error estimate forms W with dt).
.update_W <- function(nl, ctx, t, uprev, gdt, dae = NULL, gdt_test = gdt) {
  nj <- .do_newJW(nl, ctx, t, gdt_test, cj = dae$cj)
  cnt <- list(nf = 0L, njacs = 0L, nw = 0L)
  # out of place, calc_W computes J afresh whenever it forms W
  if (nj[1L] || (nj[2L] && is.null(dae))) {
    if (is.null(dae)) {
      jj <- .calc_J(ctx, uprev, t)
      nl$J <- jj$J
    } else {
      jj <- .calc_J_dae(ctx, if (is.null(ctx$jac)) dae$du_pred else ctx$du, uprev, t)
      nl$J <- jj$J
      nl$Jdu <- jj$Jdu
    }
    cnt$nf <- jj$nf
    cnt$njacs <- 1L
    nl$J_t <- t
  }
  nl$new_W <- nj[2L]
  if (nj[2L]) {
    W <- if (is.null(dae)) {
      if (is.null(ctx$mass)) nl$J - diag(1 / gdt, nrow(nl$J)) else nl$J - ctx$mass / gdt
    } else nl$J + dae$cj * nl$Jdu
    nl$W <- .factorize(W)
    nl$W_gdt <- if (is.null(dae)) gdt else dae$cj
    cnt$nw <- 1L
  }
  if (isTRUE(ctx$composite) && nj[1L]) cnt$eigen_est <- max(rowSums(abs(nl$J)))
  cnt
}

# nlsolve! with compute_step! of the Newton method. Solves for z:
#   ODE (COEFFICIENT_MULTISTEP): tmp + f(z, p, tstep) - (alpha / (gamma dt)) M z = 0
#   ODE (DIRK): dt f(tmp + gamma z, p, tstep) - M z = 0
#   DAE: F((tmp + alpha z) / (gamma dt), u0 + z, p, tstep) = 0
# `z` is the initial guess, `uprev` the state at the start of the step (for
# the error weights). A retry with a fresh Jacobian
# (TryAgain) continues from the last iterate. `firststage` is
# markfirststage!: FALSE for the later stages of a multistage method, which
# keep W after a rejected step. `eta_first` is isadaptive(alg): FALSE skips
# the rate test on the first iteration (check_eta_convergence of an algorithm
# that is not adaptive, such as CNAB2 or a fixed step SDIRK method). Returns list(z, ok, nf, njacs, nsolve, nw,
# eigen_est).
.nlsolve <- function(nl, ctx, t, dt, z, tmp, gamma, alpha, tstep, uprev, u0 = NULL,
                     method = "multistep", firststage = TRUE, eta_first = TRUE) {
  nl$isfs <- firststage
  dae <- !is.null(u0)
  dirk <- identical(method, "dirk")
  cnt <- list(nf = 0L, njacs = 0L, nsolve = 0L, nw = 0L)
  add <- function(x) for (s in names(cnt)) cnt[[s]] <<- cnt[[s]] + as.integer(x[[s]] %||% 0L)
  invgdt <- 1 / (gamma * dt)
  eig <- NULL
  repeat {
    gdt <- if (dirk) gamma * dt else gamma * dt / alpha
    # do_newJW measures the step change by 1 / (gamma dt), without alpha
    w <- .update_W(nl, ctx, t, uprev, gdt,
                   dae = if (dae) list(cj = alpha * invgdt, du_pred = tmp * invgdt) else NULL,
                   gdt_test = gamma * dt)
    add(w)
    if (!is.null(w$eigen_est)) eig <- w$eigen_est
    nl$status <- "Convergence"
    eta <- if (isTRUE(nl$new_W)) max(nl$eta_old, .Machine$double.eps)^0.8 else nl$eta_old
    ndz <- 1
    converged <- FALSE
    for (it in seq_len(nl$maxiters)) {
      ndzprev <- ndz
      if (dae) {
        ustep <- u0 + z
        ztmp <- ctx$f((tmp + alpha * z) * invgdt, ustep, ctx$p, tstep)
      } else if (dirk) {
        ustep <- tmp + gamma * z
        ztmp <- (dt * ctx$f(ustep, ctx$p, tstep) - .mass_mul(ctx$mass, z)) * invgdt
      } else {
        ustep <- z
        ztmp <- tmp + ctx$f(z, ctx$p, tstep) - (alpha * invgdt) * .mass_mul(ctx$mass, z)
      }
      cnt$nf <- cnt$nf + 1L
      if (is.null(nl$W) || any(!is.finite(ztmp))) {
        ndz <- Inf
      } else {
        dz <- .wsolve(nl$W, ztmp)
        cnt$nsolve <- cnt$nsolve + 1L
        ndz <- .err_norm(dz, uprev, ustep, ctx$abstol, ctx$reltol)
        # has_special_newton_error (QNDF, NordsieckBDF): the increment norm
        # times the error constant
        if (!is.null(ctx$newton_scale)) ndz <- ndz * ctx$newton_scale
      }
      if (!is.finite(ndz)) {
        nl$status <- "Divergence"
        nl$nfails <- nl$nfails + 1L
        break
      }
      if (it > 1L) {
        theta <- nl$prev_theta <- max(0.3 * nl$prev_theta, ndz / ndzprev)
        if (abs(theta - 1) <= 100 * sqrt(.Machine$double.eps)) {
          if (ndz <= 1) {
            z <- z - dz
            converged <- TRUE
            break
          }
          nl$status <- "Divergence"
          nl$nfails <- nl$nfails + 1L
          break
        }
        if (theta > 2) {
          nl$status <- "Divergence"
          nl$nfails <- nl$nfails + 1L
          break
        }
      } else {
        if (!isTRUE(ctx$accept_prev)) nl$prev_theta <- 1
        theta <- nl$prev_theta
      }
      z <- z - dz
      eta <- theta / (1 - theta)
      if ((it == 1L && ndz < 1e-5) || ((it > 1L || eta_first) && eta >= 0 && eta * ndz < nl$kappa)) {
        converged <- TRUE
        break
      }
    }
    if (converged) {
      nl$status <- "Convergence"
      nl$nfails <- 0L
    } else if (nl$status != "Divergence") {
      nl$status <- "Divergence"
      nl$nfails <- nl$nfails + 1L
    }
    # a stale Jacobian gets one retry with a fresh one
    if (nl$status == "Divergence" && !isTRUE(nl$J_t == t)) {
      nl$status <- "TryAgain"
      nl$nfails <- nl$nfails + 1L
      next
    }
    break
  }
  nl$eta_old <- eta
  nl$firstcall <- FALSE
  c(list(z = z, ok = nl$status == "Convergence", eigen_est = eig), cnt)
}

# Reset the convergence rate after a rejected step (step_reject_controller!).
.nl_on_reject <- function(cache) if (is.environment(cache$nl)) cache$nl$prev_theta <- 1

# ---------------------------------------------------------------------------
# ImplicitEuler and Trapezoid

.implicit_init <- function(ctx, t, u) {
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  cache$fsal <- ctx$f(u, ctx$p, t)
  cache$uprev3 <- u
  cache$tprev2 <- NA_real_
  cache
}

# After an event: the Newton state and the error history are kept, the FSAL
# derivative is evaluated at the new state.
.implicit_on_event <- function(ctx, t, u, cache) {
  cache$fsal <- ctx$f(u, ctx$p, t)
  list(cache = cache, nf = 1L)
}

# Derivatives for the Hermite dense output: M^-1 f with a mass matrix, NA
# for algebraic variables (interpolated linearly, as _dae_hermite_mask).
.mm_deriv <- function(ctx, f) {
  M <- ctx$mass
  if (is.null(M)) return(f)
  d <- diag(M)
  if (all(M[row(M) != col(M)] == 0)) return(ifelse(d != 0, f / d, NA_real_))
  tryCatch(as.numeric(solve(M, f)), error = function(e) rep(NA_real_, length(f)))
}

.hermite_dae_interp <- function(theta, dt, y0, y1, k) {
  lin <- (1 - theta) * y0 + theta * y1
  out <- .hermite_interp(theta, dt, y0, y1, k)
  ifelse(is.na(out), lin, out)
}

# Norm of the algebraic residuals over abstol, added to EEst on singular
# mass matrices (:ie_dd2 DAE tail, Rosenbrock23).
.alg_penalty <- function(ctx, fu) {
  if (is.null(ctx$algebraic) || !any(ctx$algebraic)) return(0)
  x <- ifelse(ctx$algebraic, fu, 0) / ctx$abstol
  sqrt(mean(x^2))
}

.implicit_euler_step <- function(ctx, t, u, dt, cache) {
  fprev <- cache$fsal
  tmp <- .mass_mul(ctx$mass, u) / dt
  nl <- .nlsolve(cache$nl, ctx, t, dt, u, tmp, 1, 1, t + dt, u)
  out <- list(nf = nl$nf, njacs = nl$njacs, nsolve = nl$nsolve, nw = nl$nw, eigen_est = nl$eigen_est,
              cache = cache)
  if (!nl$ok) return(c(out, list(fail = TRUE, u = u)))
  unew <- nl$z
  EEst <- 1
  if (ctx$adaptive && ctx$success_iter > 0L) {
    tprev <- ctx$tprev
    dt1 <- dt * (t + dt - tprev)
    dt2 <- (t - tprev) * (t + dt - tprev)
    r <- 7 / 12 * dt^2
    est <- r * abs((unew - u) / dt1 - (u - ctx$uprev2) / dt2)
    EEst <- .err_norm(est, u, unew, ctx$abstol, ctx$reltol)
  }
  fnew <- ctx$f(unew, ctx$p, t + dt)
  if (ctx$adaptive) EEst <- EEst + .alg_penalty(ctx, fnew)
  if (isTRUE(EEst <= 1)) cache$fsal <- fnew
  out$nf <- out$nf + 1L
  c(out, list(u = unew, EEst = EEst, k = list(.mm_deriv(ctx, fprev), .mm_deriv(ctx, fnew))))
}

.trapezoid_step <- function(ctx, t, u, dt, cache) {
  gdt <- dt / 2
  fprev <- cache$fsal
  tmp <- .mass_mul(ctx$mass, u) * (1 / gdt) + fprev
  nl <- .nlsolve(cache$nl, ctx, t, dt, u, tmp, 1 / 2, 1, t + dt, u)
  out <- list(nf = nl$nf, njacs = nl$njacs, nsolve = nl$nsolve, nw = nl$nw, eigen_est = nl$eigen_est,
              cache = cache)
  if (!nl$ok) return(c(out, list(fail = TRUE, u = u)))
  unew <- nl$z
  EEst <- 1
  if (ctx$adaptive) {
    # the third divided difference needs three distinct past points since
    # the last (re)initialization
    if (ctx$iter > 2L && ctx$success_iter > 1L && is.finite(cache$tprev2) &&
        cache$tprev2 != ctx$tprev && ctx$tprev != t) {
      tprev <- ctx$tprev; uprev2 <- ctx$uprev2
      tprev2 <- cache$tprev2; uprev3 <- cache$uprev3
      dt1 <- dt * (t + dt - tprev)
      dt2 <- (t - tprev) * (t + dt - tprev)
      dt3 <- (t - tprev) * (t - tprev2)
      dt4 <- (tprev - tprev2) * (t - tprev2)
      dt5 <- t + dt - tprev2
      r <- 7 / 12 * (dt * dt * dt) / 2
      DD31 <- (unew - u) / dt1 - (u - uprev2) / dt2
      DD30 <- (u - uprev2) / dt3 - (uprev2 - uprev3) / dt4
      # as OrdinaryDiffEqSDIRK (:trap_dd3), the norm of the third divided
      # difference weighted componentwise
      EEst <- .err_norm(rep(r * .rms_norm((DD31 - DD30) / dt5), length(u)), u, unew, ctx$abstol, ctx$reltol)
      if (isTRUE(EEst <= 1)) {
        cache$uprev3 <- uprev2
        cache$tprev2 <- tprev
      }
    } else if (ctx$success_iter > 0L) {
      cache$uprev3 <- ctx$uprev2
      cache$tprev2 <- ctx$tprev
    }
  }
  fnew <- ctx$f(unew, ctx$p, t + dt)
  if (isTRUE(EEst <= 1)) cache$fsal <- fnew
  out$nf <- out$nf + 1L
  c(out, list(u = unew, EEst = EEst, k = list(.mm_deriv(ctx, fprev), .mm_deriv(ctx, fnew))))
}

#' Implicit ODE algorithms
#'
#' Implicit methods for stiff problems and differential algebraic
#' equations in mass matrix form, as in OrdinaryDiffEq.jl.
#' `ImplicitEuler()` is the backward Euler method with a divided difference
#' error estimate, `Trapezoid()` the trapezoidal rule, and `FBDF()` the
#' variable order (1 to 5) fixed leading coefficient backward
#' differentiation formula. Each step is solved by a modified Newton
#' iteration that reuses the Jacobian across steps until convergence
#' slows down; the Jacobian is `prob$f$jac` when the problem has one
#' (`ODEProblem(sys, ..., jac = TRUE)`) and finite differences otherwise.
#'
#' @param max_order Maximum order of `FBDF()`.
#' @param autodiff How Jacobians are computed when the problem has no
#'   symbolic Jacobian: `AutoForwardDiff()` (central differences, as R has no
#'   dual numbers), `AutoFiniteDiff(fdtype)` or `AutoSparse(ad)`; see
#'   [AutoForwardDiff()]. The default uses forward differences.
#' @return An algorithm object.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -1000 * (u - cos(t)), 0, c(0, 1))
#' sol <- solve(prob, FBDF(), reltol = 1e-6)
#' sol(1)
#' sol2 <- solve(prob, ImplicitEuler())
#' sol2$retcode
#' @name implicit_algorithms
NULL

#' @rdname implicit_algorithms
#' @export
ImplicitEuler <- function(autodiff = NULL) .with_autodiff(autodiff, {
  .new_alg("ImplicitEuler", 1L, step = .implicit_euler_step, interp = .hermite_dae_interp,
           init = .implicit_init, on_event = .implicit_on_event,
           stiff = TRUE, mass_matrix = TRUE, qsteady_max = 6 / 5, on_reject = .nl_on_reject)
})

#' @rdname implicit_algorithms
#' @export
Trapezoid <- function(autodiff = NULL) .with_autodiff(autodiff, {
  .new_alg("Trapezoid", 2L, step = .trapezoid_step, interp = .hermite_dae_interp,
           init = .implicit_init, on_event = .implicit_on_event,
           stiff = TRUE, mass_matrix = TRUE, qsteady_max = 6 / 5, on_reject = .nl_on_reject)
})

# ---------------------------------------------------------------------------
# FBDF: fixed leading coefficient BDF of variable order (bdf_perform_step.jl,
# FBDFConstantCache). The history ts / u_history holds past accepted points;
# the predictor and the corrector values come from the Lagrange polynomial
# through them. Orders 3 to 5 are guarded by stability limit detection
# (stald.jl); the optional time filters (time_filter = FALSE by default, as
# in Julia) add filtered candidates of neighboring orders.

# _make_bdf_coeffs_fbdf and _make_bdf_coeffs_dfbdf (bdf_utils.jl)
.BDF_COEFFS_FBDF <- rbind(
  c(1, -1, 0, 0, 0, 0),
  c(3 / 2, -2, 1 / 2, 0, 0, 0),
  c(11 / 6, -3, 3 / 2, -1 / 3, 0, 0),
  c(25 / 12, -4, 3, -4 / 3, 1 / 4, 0),
  c(137 / 60, -5, 5, -10 / 3, 5 / 4, -1 / 5))
.BDF_COEFFS_DFBDF <- rbind(
  c(1, -1, 0, 0, 0, 0),
  c(2 / 3, -4 / 3, 1 / 3, 0, 0, 0),
  c(6 / 11, -18 / 11, 9 / 11, -2 / 11, 0, 0),
  c(12 / 25, -48 / 25, 36 / 25, -16 / 25, 3 / 25, 0),
  c(60 / 137, -300 / 137, 300 / 137, -200 / 137, 75 / 137, -12 / 137))

# _eval_lagrange_oop: the interpolant through (thetas[j], values[[j]]) at xi
# (the basis as Julia forms it, L *= (xi - t_m) / (t_j - t_m))
.eval_lagrange <- function(xi, thetas, values, n = length(thetas)) {
  out <- 0
  for (j in seq_len(n)) {
    L <- 1
    for (m in seq_len(n)) if (m != j) L <- L * ((xi - thetas[m]) / (thetas[j] - thetas[m]))
    out <- out + L * values[[j]]
  }
  out
}

# calc_finite_difference_weights! (Fornberg's algorithm).
.fd_weights <- function(ts, t, order, size = 6L) {
  c <- matrix(0, size, size)
  c1 <- 1
  c4 <- ts[1L] - t
  c[1L, 1L] <- 1
  if (order >= 1L) for (i in 2:(order + 1L)) {
    c2 <- 1
    c5 <- c4
    c4 <- ts[i] - t
    for (j in 1:(i - 1L)) {
      c3 <- ts[i] - ts[j]
      c2 <- c2 * c3
      if (j == i - 1L) {
        for (k in i:2) c[i, k] <- c1 * ((k - 1) * c[i - 1L, k - 1L] - c5 * c[i - 1L, k]) / c2
        c[i, 1L] <- 0
      }
      for (k in i:2) c[j, k] <- (c4 * c[j, k] - (k - 1) * c[j, k - 1L]) / c3
      c[j, 1L] <- 0
    }
    c1 <- c2
  }
  c
}

# estimate_terk: h^(k-1) times the (k-1)th derivative through u and the
# first k-1 history values.
.estimate_terk <- function(cache, k, u, t, dt) {
  w <- .fd_weights(cache$ts_tmp, t + dt, k - 1L)
  terk <- w[1L, k] * u
  if (k >= 2L) for (i in 2:k) terk <- terk + w[i, k] * cache$u_history[[i - 1L]]
  terk * abs(dt^(k - 1L))
}

# ---------------------------------------------------------------------------
# Stability limit detection (stald.jl). BDF orders 3 to 5 are only alpha
# stable; when an eigenvalue of h J nears the stability boundary the scaled
# derivative norms of the last five steps follow a geometric pattern whose
# ratio rr (the characteristic root) exceeds rrcut, and the order is lowered.

.stald_new <- function(enabled = TRUE, rrcut = 0.98, vrrtol = 1e-4, vrrt2 = 5e-4, sqtol = 1e-3,
                       rrtol = 1e-2, tiny = 0x1.04bd984990e6fp-299) {
  st <- new.env(parent = emptyenv())
  st$enabled <- isTRUE(enabled)
  st$rrcut <- rrcut; st$vrrtol <- vrrtol; st$vrrt2 <- vrrt2; st$sqtol <- sqtol
  st$rrtol <- rrtol; st$tiny <- tiny
  .stald_reset(st)
  st
}

# stald_reset!
.stald_reset <- function(st) {
  if (is.null(st)) return(invisible())
  st$ssdat <- matrix(0, 5L, 3L)
  st$nscon <- 0L
  st$last_order <- 0L
  if (is.null(st$nor)) st$nor <- 0L
  invisible()
}

# stald_collect_data!: squared norms of the order q-2, q-1 and q estimates,
# newest in row 1.
.stald_collect <- function(st, order, terkm2, terkm1, terk) {
  if (!st$enabled || order < 3L) return(invisible())
  st$ssdat[2:5, ] <- st$ssdat[1:4, ]
  st$ssdat[1L, ] <- c(terkm2, terkm1, terk)^2
  invisible()
}

# stald_check!: TRUE when the order used for this step should be lowered.
.stald_check <- function(st, order) {
  if (!st$enabled) return(FALSE)
  if (order < 3L) {
    st$last_order <- order
    return(FALSE)
  }
  if (order != st$last_order) {
    st$nscon <- 0L
    st$last_order <- order
    return(FALSE)
  }
  st$nscon <- st$nscon + 1L
  if (st$nscon < order + 5L) return(FALSE)
  if (.stald_detect(st, order) > 3L) {
    st$nor <- st$nor + 1L
    st$nscon <- 0L
    return(TRUE)
  }
  FALSE
}

# _stald_detect (cvSLdet): negative when the data do not determine the
# characteristic root, 1 to 3 for a stable root found by the normal, quartic
# or Newton corrected quartic case, 4 to 6 for the same with rr > rrcut.
.stald_detect <- function(st, q) {
  ssdat <- st$ssdat
  tiny <- st$tiny
  smax <- ssmax <- rav <- vrat <- numeric(3L)
  qc <- matrix(0, 5L, 3L)
  for (k in 1:3) {
    smink <- min(ssdat[, k])
    smaxk <- max(0, ssdat[, k])
    if (isTRUE(smink < tiny * smaxk)) return(-1L)
    smax[k] <- smaxk
    ssmax[k] <- smaxk * smaxk
    rat <- ssdat[1:4, k] / ssdat[2:5, k]
    rav[k] <- 0.25 * sum(rat)
    vrat[k] <- abs(0.25 * sum(rat * rat) - rav[k] * rav[k])
    qc[5L, k] <- ssdat[1L, k] * ssdat[3L, k] - ssdat[2L, k] * ssdat[2L, k]
    qc[4L, k] <- ssdat[2L, k] * ssdat[3L, k] - ssdat[1L, k] * ssdat[4L, k]
    qc[3L, k] <- 0
    qc[2L, k] <- ssdat[2L, k] * ssdat[5L, k] - ssdat[3L, k] * ssdat[4L, k]
    qc[1L, k] <- ssdat[4L, k] * ssdat[4L, k] - ssdat[3L, k] * ssdat[5L, k]
  }
  quartic <- function(k, r) qc[5L, k] + r * (qc[4L, k] + r * r * (qc[2L, k] + r * qc[1L, k]))
  vmin <- min(vrat); vmax <- max(vrat)
  kflag <- 0L
  if (isTRUE(vmin < st$vrrtol * st$vrrtol)) {
    if (isTRUE(vmax > st$vrrt2 * st$vrrt2)) return(-2L)
    rr <- sum(rav) / 3
    if (isTRUE(max(abs(rav - rr)) > st$vrrt2)) return(-3L)
    kflag <- 1L
  } else {
    qco <- qc
    if (isTRUE(abs(qco[4L, 1L]) > tiny * smax[1L] * smax[1L])) {
      r1 <- qco[4L, 2L] / qco[4L, 1L]
      r2 <- qco[4L, 3L] / qco[4L, 1L]
      for (j in c(5L, 3L, 2L, 1L)) {
        qco[j, 2L] <- qco[j, 2L] - r1 * qco[j, 1L]
        qco[j, 3L] <- qco[j, 3L] - r2 * qco[j, 1L]
      }
    } else return(-4L)
    if (isTRUE(abs(qco[5L, 2L]) > tiny * smax[2L] * smax[2L])) {
      r1 <- qco[5L, 3L] / qco[5L, 2L]
      for (j in c(3L, 2L, 1L)) qco[j, 3L] <- qco[j, 3L] - r1 * qco[j, 2L]
    } else if (isTRUE(abs(qco[5L, 3L]) > tiny * smax[3L] * smax[3L])) return(-4L)
    if (isTRUE(abs(qco[4L, 3L]) < tiny)) return(-4L)
    rr <- -qco[5L, 3L] / qco[4L, 3L]
    if (isTRUE(rr < tiny) || isTRUE(rr > 100)) return(-5L)
    sqmax <- max(vapply(1:3, function(k) abs(quartic(k, rr)) / ssmax[k], 0))
    if (isTRUE(sqmax < st$sqtol)) {
      kflag <- 2L
    } else {
      # Newton corrections of rr
      for (it in 1:3) {
        rrc <- vapply(1:3, function(k) {
          qp <- qc[4L, k] + rr * rr * (3 * qc[2L, k] + rr * 4 * qc[1L, k])
          rr + if (isTRUE(abs(qp) > tiny * ssmax[k])) -quartic(k, rr) / qp else 0
        }, 0)
        sqmin <- Inf
        kmin <- 0L
        for (k in 1:3) if (isTRUE(rrc[k] > 0)) {
          saqk <- abs(quartic(k, rrc[k])) / ssmax[k]
          if (isTRUE(saqk < sqmin)) {
            sqmin <- saqk
            kmin <- k
          }
        }
        if (kmin == 0L) return(-6L)
        rr <- rrc[kmin]
        if (isTRUE(sqmin < st$sqtol)) {
          kflag <- 3L
          break
        }
      }
      if (kflag == 0L) return(-6L)
    }
  }
  sigsq <- numeric(3L)
  for (k in 1:3) {
    rsa <- ssdat[1L, k]
    rsb <- ssdat[2L, k] * rr
    rsc <- ssdat[3L, k] * rr * rr
    rsd <- ssdat[4L, k] * rr * rr * rr
    rd1a <- rsa - rsb; rd1b <- rsb - rsc; rd1c <- rsc - rsd
    rd2a <- rd1a - rd1b; rd2b <- rd1b - rd1c
    rd3a <- rd2a - rd2b
    if (isTRUE(abs(rd1b) < tiny * smax[k])) return(-7L)
    cest1 <- -rd3a / rd1b
    if (isTRUE(cest1 < tiny) || isTRUE(cest1 > 4)) return(-7L)
    sigsq[k] <- ssdat[3L, k] + (rd2b / cest1) / (rr * rr)
  }
  if (isTRUE(sigsq[2L] < tiny)) return(-8L)
  ratp <- sigsq[3L] / sigsq[2L]
  ratm <- sigsq[1L] / sigsq[2L]
  qfac1 <- 0.25 * (q * q - 1)
  qfac2 <- 2 / (q - 1)
  bb <- ratp * ratm - 1 - qfac1 * ratp
  tem <- 1 - qfac2 * bb
  if (isTRUE(abs(tem) < tiny)) return(-8L)
  if (isTRUE(abs(1 / tem - rr) > st$rrtol)) return(-9L)
  if (isTRUE(rr > st$rrcut)) return(kflag + 3L)
  kflag
}

# ---------------------------------------------------------------------------
# Time filters (bdf_utils.jl, after DeCaria et al., arXiv:1810.06670)

# backdiff!: c[q + 1, i] is the weight of y(ts[i]) in the q-th divided
# difference over the newest q + 1 of the ascending points ts[1:n].
.filt_backdiff <- function(ts, n) {
  D <- matrix(0, n, n)
  D[cbind(seq_len(n), n + 1L - seq_len(n))] <- 1
  cc <- matrix(0, n, n)
  cc[1L, ] <- D[1L, ]
  if (n >= 2L) for (q in 1:(n - 1L)) {
    for (j in seq_len(n - q)) D[j, ] <- (D[j, ] - D[j + 1L, ]) / (ts[n + 1L - j] - ts[n + 1L - j - q])
    cc[q + 1L, ] <- D[1L, ]
  }
  cc
}

# _filt_divided_diff: row `row` of the table applied to y (the newest value)
# and the history (newest first).
.filt_divided_diff <- function(cc, row, n, y, hist) {
  d <- cc[row, n] * y
  for (i in seq_len(n - 1L)) d <- d + cc[row, i] * hist[[n - i]]
  d
}

# _fbdf_filter_weight: cancels the fixed coefficient solve's error on a monic
# polynomial of degree k + 1.
.fbdf_filter_weight <- function(cc, ts, n, k, dt, coeffs) {
  defect <- 0
  for (i in seq_len(k)) {
    th <- -i
    prod <- 1
    for (j in seq_len(n - 1L)) prod <- prod * (th - (ts[j] - ts[n]) / dt)
    defect <- defect - coeffs[k, i + 1L] * (th^(k + 1L) - prod)
  }
  err <- defect * dt^(k + 1L) / coeffs[k, 1L]
  err / (1 + cc[n, n] * err)
}

# _filt_bdf_coefficients!: the order p variable step BDF coefficients.
.filt_bdf_coefficients <- function(cc, ts, n, p) {
  ab <- numeric(n)
  for (i in seq_len(n)) {
    pre <- 1
    for (j in seq_len(p)) {
      ab[i] <- ab[i] + pre * cc[j + 1L, i]
      pre <- pre * (ts[n] - ts[n - j])
    }
  }
  ab
}

.filt_step_ratio <- function(est, p) {
  if (isTRUE(est == 0)) return(Inf)
  if (!is.finite(est) || est <= 0) return(0)
  est^(-1 / (p + 1))
}

# _fbdf_time_filter: returns list(u, fsal, nf, EEst) and sets
# cache$filter_order (0 when no filter applies).
.fbdf_time_filter <- function(ctx, cache, t, dt, uprev, u_bdf, k, EEst) {
  cache$filter_order <- 0L
  tdt <- t + dt
  f <- function(y) ctx$f(y, ctx$p, tdt)
  hist <- cache$u_history
  maxo <- cache$max_order
  if (!(k <= 4L && cache$iters_from_event >= k)) return(list(u = u_bdf, fsal = f(u_bdf), nf = 1L, EEst = EEst))
  n <- k + 2L
  ts_asc <- c(cache$ts[(n - 1L):1], tdt)
  cc <- .filt_backdiff(ts_asc, n)
  eta <- .fbdf_filter_weight(cc, ts_asc, n, k, dt, cache$bdf_coeffs)
  est <- -eta * .filt_divided_diff(cc, n, n, u_bdf, hist)
  y_hi <- u_bdf + est
  if (!isTRUE(ctx$adaptive)) {
    cache$filter_order <- min(k + 1L, maxo)
    y <- if (k < maxo) y_hi else u_bdf
    return(list(u = y, fsal = f(y), nf = 1L, EEst = EEst))
  }
  nrm <- function(x, y) .err_norm(x, uprev, y, ctx$abstol, ctx$reltol)
  err <- nrm(est, u_bdf)
  f_hi <- f(y_hi)
  ab <- .filt_bdf_coefficients(cc, ts_asc, n, k + 1L)
  res <- ab[n] * y_hi
  for (i in seq_len(n - 1L)) res <- res + ab[i] * hist[[n - i]]
  err_hi <- nrm((res - f_hi) / ab[n], y_hi)
  y_lo <- u_bdf
  err_lo <- Inf
  if (k == 3L) {
    # BDF3-Stab (equation 3.16, mu = 9/125)
    est_lo <- (9 / 125) / cc[4L, n] * .filt_divided_diff(cc, 4L, n, u_bdf, hist)
    y_lo <- u_bdf + est_lo
    err_lo <- nrm(est_lo, y_lo)
  }
  # _fbdf_select_filter!
  order <- k
  estimate <- err
  ratio <- .filt_step_ratio(err, k)
  if (k == 3L && .filt_step_ratio(err_lo, 2L) > ratio) {
    order <- 2L
    estimate <- err_lo
    ratio <- .filt_step_ratio(err_lo, 2L)
  }
  if (k < maxo && cache$qwait == 0L && .filt_step_ratio(err_hi, k + 1L) > ratio) {
    order <- k + 1L
    estimate <- err_hi
  }
  cache$filter_order <- order
  if (order == k + 1L) return(list(u = y_hi, fsal = f_hi, nf = 1L, EEst = estimate))
  y <- if (order < k) y_lo else u_bdf
  list(u = y, fsal = f(y), nf = 2L, EEst = estimate)
}

.bdf_init <- function(ctx, t, u, max_order = 5L, dae = FALSE, opts = list()) {
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  cache$max_order <- max_order
  cache$dae <- dae
  cache$bdf_coeffs <- if (dae) .BDF_COEFFS_DFBDF else .BDF_COEFFS_FBDF
  cache$ts <- numeric(max_order + 2L)
  cache$ts_tmp <- numeric(max_order + 2L)
  cache$u_history <- rep(list(u * 0), max_order + 2L)
  cache$order <- 1L
  cache$prev_order <- 1L
  cache$nconsteps <- 0L
  cache$consfailcnt <- 0L
  cache$qwait <- 3L
  cache$terkm2 <- 1; cache$terkm1 <- 1; cache$terk <- 1; cache$terkp1 <- 1
  cache$iters_from_event <- 0L
  cache$filter_order <- 0L
  cache$time_filter <- isTRUE(opts$time_filter) && !dae
  if (cache$time_filter && !is.null(ctx$mass))
    stop("FBDF(time_filter = TRUE) requires the identity mass matrix; use time_filter = FALSE",
         call. = FALSE)
  # STALD runs for FBDF only (DFBDF has none in Julia)
  cache$stald <- if (dae) NULL else do.call(.stald_new, opts$stald %||% list())
  cache$fsal <- if (dae) NULL else ctx$f(u, ctx$p, t)
  # reinitFBDF! after a derivative discontinuity
  cache$ts[1L] <- t
  cache$u_history[[1L]] <- u
  cache
}

# reinitFBDF!: shift the history when the previous step was accepted.
.bdf_reinit <- function(cache, t, uprev) {
  ord <- cache$order
  if (cache$iters_from_event == 0L) {
    cache$ts[1L] <- t
    cache$u_history[[1L]] <- uprev
  } else if (cache$iters_from_event == 1L && t != cache$ts[1L]) {
    cache$ts[2L] <- cache$ts[1L]
    cache$ts[1L] <- t
    cache$u_history[[2L]] <- cache$u_history[[1L]]
    cache$u_history[[1L]] <- uprev
  } else if (cache$consfailcnt == 0L && t != cache$ts[1L]) {
    for (i in (ord + 2L):2L) {
      cache$ts[i] <- cache$ts[i - 1L]
      cache$u_history[[i]] <- cache$u_history[[i - 1L]]
    }
    cache$ts[1L] <- t
    cache$u_history[[1L]] <- uprev
  }
}

# The shared FBDF / DFBDF step: predictor, corrector stage solve, local
# error and the terk estimates for order selection.
.bdf_step <- function(ctx, t, u, dt, cache) {
  .bdf_reinit(cache, t, u)
  k <- cache$order
  ts <- cache$ts
  hist <- cache$u_history
  coeffs <- cache$bdf_coeffs
  tdt <- t + dt
  n_pred <- k + 1L
  pred_thetas <- (ts[seq_len(n_pred)] - t) / dt
  # without history the guess is integrator.u: the result of the rejected
  # attempt when the step is retried
  u0 <- if (cache$iters_from_event >= 1L) .eval_lagrange(1, pred_thetas, hist, n_pred) else
    if (isTRUE(cache$u_att_t == t)) cache$u_att else u
  u_corr <- list()
  if (cache$iters_from_event >= 1L && k >= 2L)
    for (i in 1:(k - 1L)) u_corr[[i]] <- .eval_lagrange(-i, pred_thetas, hist, n_pred)
  else if (k >= 2L) for (i in 1:(k - 1L)) u_corr[[i]] <- u * 0
  if (!cache$dae) {
    tmp <- -u * coeffs[k, 2L]
    if (k >= 2L) for (i in 1:(k - 1L)) tmp <- tmp - u_corr[[i]] * coeffs[k, i + 2L]
    tmp <- .mass_mul(ctx$mass, tmp) / dt
    gamma <- 1 / coeffs[k, 1L]
    nl <- .nlsolve(cache$nl, ctx, t, dt, u0, tmp, gamma, 1, tdt, u)
  } else {
    tmp <- u * coeffs[k, 2L]
    if (k >= 2L) for (i in 1:(k - 1L)) tmp <- tmp + u_corr[[i]] * coeffs[k, i + 2L]
    tmp <- tmp + u0
    gamma <- coeffs[k, 1L]
    nl <- .nlsolve(cache$nl, ctx, t, dt, u0 * 0, tmp, gamma, 1, tdt, u, u0 = u0)
  }
  out <- list(nf = nl$nf, njacs = nl$njacs, nsolve = nl$nsolve, nw = nl$nw, eigen_est = nl$eigen_est,
              cache = cache)
  if (!nl$ok) return(c(out, list(fail = TRUE, u = u)))
  z <- nl$z
  unew <- if (cache$dae) z + u0 else z
  r <- numeric(k)
  if (k >= 2L) for (j in 2:k) {
    r[j] <- 1 - j
    for (i in 2:(k + 1L)) r[j] <- r[j] * ((tdt - j * dt) - ts[i]) / (i * dt)
  }
  terkp1 <- if (cache$dae) z else unew - u0
  for (j in 1:(k + 1L)) terkp1 <- terkp1 * (j * dt / (tdt - ts[j]))
  lte <- -1 / (1 + k)
  if (k >= 2L) for (j in 2:k)
    lte <- lte - (if (cache$dae) coeffs[k, j] / coeffs[k, 1L] else coeffs[k, j]) * r[j]
  lte <- lte * terkp1
  EEst <- 0
  if (ctx$adaptive) {
    cache$ts_tmp[1L] <- tdt
    cache$ts_tmp[2:(k + 2L)] <- ts[1:(k + 1L)]
    nrm <- function(x) .err_norm(x, u, unew, ctx$abstol, ctx$reltol)
    EEst <- nrm(lte)
    cache$terk <- nrm(.estimate_terk(cache, k + 1L, unew, t, dt))
    if (k > 1L) cache$terkm1 <- nrm(.estimate_terk(cache, k, unew, t, dt))
    if (k > 2L) cache$terkm2 <- nrm(.estimate_terk(cache, k - 1L, unew, t, dt))
    cache$terkp1 <- if (cache$qwait == 0L && k < cache$max_order) nrm(terkp1) else 0
  }
  du <- NULL
  if (cache$dae) {
    du <- (z + tmp) / (gamma * dt)
  } else if (cache$time_filter) {
    fl <- .fbdf_time_filter(ctx, cache, t, dt, u, unew, k, EEst)
    unew <- fl$u
    EEst <- fl$EEst
    cache$fsal <- fl$fsal
    out$nf <- out$nf + fl$nf
    if (!ctx$adaptive) {
      # _fbdf_finish_fixed_step!
      cache$prev_order <- cache$order
      cache$order <- max(cache$order, cache$filter_order)
      cache$iters_from_event <- cache$iters_from_event + 1L
      cache$nconsteps <- cache$nconsteps + 1L
      cache$consfailcnt <- 0L
    }
  } else {
    cache$fsal <- ctx$f(unew, ctx$p, tdt)
    out$nf <- out$nf + 1L
  }
  # dense output: the Lagrange polynomial through u_new and the history
  cache$u_att <- unew
  cache$u_att_t <- t
  dense_order <- if (cache$time_filter) max(k, cache$filter_order) else k
  nk <- min(dense_order + 1L, cache$max_order + 1L)
  thetas <- c(1, (ts[seq_len(nk - 1L)] - t) / dt)
  vals <- c(list(unew), hist[seq_len(nk - 1L)])
  c(out, list(u = unew, EEst = EEst, du = du, k = list(thetas = thetas, values = vals)))
}

.bdf_interp <- function(theta, dt, y0, y1, k) .eval_lagrange(theta, k$thetas, k$values)

# choose_order! (controllers.jl, constant cache version)
.bdf_choose_order <- function(cache, t, dt, u, uprev, ctx) {
  k <- cache$order
  terkm2 <- cache$terkm2; terkm1 <- cache$terkm1; terk <- cache$terk; terkp1 <- cache$terkp1
  if (k < cache$max_order && cache$qwait == 0L &&
      ((k == 1L && terk > terkp1) ||
       (k == 2L && terkm1 > terk && terk > terkp1) ||
       (k > 2L && terkm2 > terkm1 && terkm1 > terk && terk > terkp1))) {
    k <- k + 1L
    terk <- terkp1
  } else {
    while (!(terkm2 > terkm1 && terkm1 > terk && terk > terkp1) && k > 2L) {
      terkp1 <- terk
      terk <- terkm1
      terkm1 <- terkm2
      w <- .fd_weights(cache$ts_tmp, t + dt, k - 2L)
      # fd_weights[k - 2, 1] as in choose_order! of controllers.jl
      tt <- w[k - 2L, 1L] * u
      if (k - 2L >= 2L) for (i in 2:(k - 2L)) tt <- tt + w[i, k - 2L] * cache$u_history[[i - 1L]]
      tt <- tt * abs(dt^(k - 2L))
      terkm2 <- .err_norm(tt, uprev, u, ctx$abstol, ctx$reltol)
      k <- k - 1L
    }
  }
  list(k = k, terk = terk)
}

# The driver calls these through the "BDF" controller type; `cache$step_info`
# holds what the Julia controller reads from the integrator.
.bdf_stepsize <- function(cache, EEst, co, first) {
  si <- cache$step_info
  cache$prev_order <- cache$order
  qmax <- if (first) co$qmax_first else co$qmax
  if (cache$filter_order > 0L) {
    # a filtered step: its candidate's order and error estimate
    .stald_reset(cache$stald)
    q <- co$gamma / .filt_step_ratio(EEst, cache$filter_order)
    return(min(max(q, 1 / qmax), 1 / co$qmin))
  }
  # a non finite error estimate only rejects the step (bdf_reject shrinks dt)
  if (!is.finite(EEst) || !all(is.finite(c(cache$terkm2, cache$terkm1, cache$terk, cache$terkp1))))
    return(1 / co$qmin)
  # STALD on the order of the step just taken, before the order selection
  step_order <- cache$prev_order
  stald_reduce <- FALSE
  if (!is.null(cache$stald) && step_order >= 3L) {
    .stald_collect(cache$stald, step_order, cache$terkm2, cache$terkm1, cache$terk)
    stald_reduce <- .stald_check(cache$stald, step_order)
  }
  ch <- .bdf_choose_order(cache, si$t, si$dt, si$u, si$uprev, si$ctx)
  k <- ch$k
  terk <- ch$terk
  if (k != cache$order) {
    cache$nconsteps <- 0L
    cache$order <- k
  }
  if (stald_reduce) {
    k <- min(k, step_order - 1L)
    cache$order <- k
    cache$nconsteps <- 0L
    terk <- cache$terkm1
  }
  if (terk == 0) return(1 / qmax)
  alpha0 <- cache$bdf_coeffs[k, 1L]
  (6 * terk / (alpha0 * (k + 1)))^(1 / (k + 1))
}

.bdf_accept <- function(cache, q, dt, co) {
  cache$consfailcnt <- 0L
  if (cache$filter_order > cache$order) {
    cache$order <- min(cache$filter_order, cache$max_order)
    cache$nconsteps <- 0L
  }
  if (q <= co$qsteady_max && q >= co$qsteady_min) q <- 1
  cache$nconsteps <- cache$nconsteps + 1L
  cache$iters_from_event <- cache$iters_from_event + 1L
  if (cache$order != cache$prev_order) cache$qwait <- cache$order + 2L
  else if (cache$qwait > 0L) cache$qwait <- cache$qwait - 1L
  dt / q
}

# bdf_step_reject_controller! with EEst1 = terkm1 (Inf with error order
# filter_order after a filtered step).
.bdf_reject <- function(cache, EEst, dt, co) {
  cache$nl$prev_theta <- 1
  k <- cache$order
  filtered <- cache$filter_order > 0L
  EEst1 <- if (filtered) Inf else cache$terkm1
  error_order <- if (filtered) cache$filter_order else k
  h <- abs(dt)
  cache$consfailcnt <- cache$consfailcnt + 1L
  cache$nconsteps <- 0L
  if (cache$consfailcnt > 1L) h <- h / 2
  z <- co$gamma * EEst^(1 / (error_order + 1))
  hk <- if (z <= 10) h / z else 0.1 * h
  hn <- hk
  kn <- k
  if (k > 1L) {
    zk1 <- 1.3 * EEst1^(1 / k)
    hk1 <- if (zk1 <= 10) h / zk1 else 0.1 * h
    if (cache$consfailcnt > 2L || abs(hk1) > abs(hk)) {
      hn <- min(abs(h), abs(hk1))
      kn <- k - 1L
    }
  }
  # restart from first order (clear the history) after repeated failures
  if (kn == 1L && cache$consfailcnt > 3L) .bdf_restart(cache)
  cache$order <- kn
  hn
}

# reinitFBDF! after derivative_discontinuity!(integrator, true)
.bdf_restart <- function(cache) {
  cache$order <- 1L
  cache$consfailcnt <- 0L
  cache$nconsteps <- 0L
  cache$qwait <- 3L
  cache$iters_from_event <- 0L
  .stald_reset(cache$stald)
  # the history is cleared and rebuilt from the current point
  cache$ts[] <- 0
  cache$u_history <- lapply(cache$u_history, function(h) h * 0)
}

# post_newton_controller!
.bdf_newton_fail <- function(cache, dt, co) {
  if (cache$order > 1L && cache$nl$nfails >= 3L) cache$order <- cache$order - 1L
  cache$consfailcnt <- cache$consfailcnt + 1L
  cache$nconsteps <- 0L
  dt / co$failfactor
}

# Wrap the step so the controller hooks see the step's t, dt, u and uprev.
.bdf_step_wrapped <- function(ctx, t, u, dt, cache) {
  res <- .bdf_step(ctx, t, u, dt, cache)
  cache$step_info <- list(t = t, dt = dt, u = res$u, uprev = u, ctx = ctx)
  res
}

# The FBDF and DFBDF constructors' controller defaults (algorithms.jl,
# default_controller of controllers.jl): qmax 10, qsteady_max 2, gamma 1.2.
.new_bdf_alg <- function(name, max_order, dae, qmax, qsteady_min, qsteady_max, opts = list(),
                         step_limiter = NULL) {
  max_order <- as.integer(max_order)
  if (length(max_order) != 1L || is.na(max_order) || max_order < 1L || max_order > 5L)
    stop("max_order must be between 1 and 5", call. = FALSE)
  .new_alg(name, 1L, step = .bdf_step_wrapped, interp = .bdf_interp,
           init = function(ctx, t, u) .bdf_init(ctx, t, u, max_order, dae, opts),
           on_event = if (!dae) function(ctx, t, u, cache) {
             # reinitFBDF! on the derivative discontinuity; the Newton state stays
             .bdf_restart(cache)
             .implicit_on_event(ctx, t, u, cache)
           },
           stiff = TRUE, mass_matrix = !dae, dae = dae, controller = "BDF",
           qmax = qmax, qsteady_min = qsteady_min, qsteady_max = qsteady_max, gamma = 12 / 10,
           bdf_stepsize = .bdf_stepsize, bdf_accept = .bdf_accept, bdf_reject = .bdf_reject,
           newton_fail = .bdf_newton_fail, max_order = max_order, step_limiter = step_limiter,
           time_filter = isTRUE(opts$time_filter), stald = !dae && !isFALSE(opts$stald$enabled))
}

#' @rdname implicit_algorithms
#' @param stald Stability limit detection for orders 3 to 5 (CVODE's
#'   algorithm): lowers the order when the step size nears the stability
#'   boundary of the BDF formula.
#' @param stald_rrcut,stald_vrrtol,stald_vrrt2,stald_sqtol,stald_rrtol,stald_tiny
#'   Tuning constants of the stability limit detection (CVODE's defaults).
#' @param qmax,qsteady_min,qsteady_max Largest step growth factor and the
#'   range of factors that keep the step size.
#' @param step_limiter A function `limiter(u, integ, p, t)` returning the
#'   state after each accepted step; prefer the `step_limiter` option of
#'   [solve()], which every method honors.
#' @param autodiff As for [ImplicitEuler()].
#' @param time_filter Apply the time filters of DeCaria et al. (adds
#'   filtered candidates of orders `k + 1` and, at order 3, BDF3-Stab;
#'   the step uses the candidate with the largest step ratio). Requires an
#'   identity mass matrix.
#' @export
FBDF <- function(max_order = 5L, stald = TRUE, stald_rrcut = 0.98, stald_vrrtol = 1e-4,
                 stald_vrrt2 = 5e-4, stald_sqtol = 1e-3, stald_rrtol = 1e-2, stald_tiny = 0x1.04bd984990e6fp-299,
                 qmax = 10, qsteady_min = 9 / 10, qsteady_max = 2, step_limiter = NULL,
                 time_filter = FALSE, autodiff = NULL) .with_autodiff(autodiff, {
  if (!is.null(step_limiter) && !is.function(step_limiter))
    stop("step_limiter must be a function(u, integ, p, t)", call. = FALSE)
  opts <- list(time_filter = isTRUE(time_filter),
               stald = list(enabled = isTRUE(stald), rrcut = stald_rrcut, vrrtol = stald_vrrtol,
                            vrrt2 = stald_vrrt2, sqtol = stald_sqtol, rrtol = stald_rrtol,
                            tiny = stald_tiny))
  .new_bdf_alg("FBDF", max_order, dae = FALSE, qmax = qmax, qsteady_min = qsteady_min,
               qsteady_max = qsteady_max, opts = opts, step_limiter = step_limiter)
})
