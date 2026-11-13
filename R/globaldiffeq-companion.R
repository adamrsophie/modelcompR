# Global error estimation and control with a defect driven companion
# equation, after GlobalDiffEq.jl v1.7.3: src/companion.jl (the
# GlobalErrorEquation and GlobalErrorMode enums, the companion right-hand
# sides, the InterpolatingMode and SimultaneousMode estimators and the
# endpoint error refinement loop) and src/estimation.jl
# (GlobalErrorEstimation, global_error_estimate(prob, alg)).
#
# After a dense forward solve with interpolant P(t), the global error
# eps(t) ~ y(t) - P(t) solves a companion ODE driven by the defect
# d(t) = f(P(t), p, t) - P'(t), with eps(t0) = 0. P' is the first derivative
# of the dense output (Julia's sol(t, Val{1})), ported below for the
# interpolants of the solvers.

# ---------------------------------------------------------------------------
# Enums. Each instance is a classed string; the enum type is the list of its
# instances (as ReturnCode), and `inherits(x, type)` is Julia's `x isa type`.

.gde_enum <- function(type, name, value) structure(name, value = value, class = c(type, "GlobalDiffEqEnum"))

#' @export
print.GlobalDiffEqEnum <- function(x, ...) {
  cat(unclass(x)[1L], "::", class(x)[1L], " = ", attr(x, "value"), "\n", sep = "")
  invisible(x)
}

#' Options of the global error estimators
#'
#' Enumerations that select how [GlobalErrorEstimation()] and
#' [GlobalAdjoint()] work, as in GlobalDiffEq.jl. Each enumeration is a list
#' of its instances (`GlobalErrorMode$SimultaneousMode`), and each instance
#' is also exported by name.
#'
#' `GlobalErrorEquation` selects the companion ODE of
#' [GlobalErrorEstimation()]:
#' * `DefectCorrection` integrates the nonlinear correction equation
#'   `eps' = f(P + eps, p, t) - P'(t)`, `eps(t0) = 0`, where `P` is the dense
#'   interpolant of the numerical solution. It needs no Jacobian, only extra
#'   evaluations of `f` (Zadunaisky 1976; Dormand, Duckers and Prince 1984).
#' * `ErrorTransport` integrates the linearized error transport (first
#'   variational) equation `eps' = J(P, p, t) eps + d(t)` with the defect
#'   `d(t) = f(P, p, t) - P'(t)` and `J = df/du` applied as a Jacobian
#'   vector product (Shampine 1986; Berzins 1988), the first order
#'   linearization of `DefectCorrection`.
#'
#' `GlobalErrorMode` selects how the companion ODE is integrated:
#' * `InterpolatingMode` solves the problem densely, saving every step, and
#'   integrates the companion against that interpolant in a second solve.
#'   It works with any dense solver but stores the whole solution.
#' * `SimultaneousMode` drives the forward integrator one accepted step at
#'   a time and advances the companion across each step with that step's
#'   dense output: no stored solution and no second solve. Its estimate
#'   closely agrees with `InterpolatingMode`.
#'
#' `GlobalAdjointScope` selects what [GlobalAdjoint()] estimates and
#' reports in `sol$global_error`:
#' * `EndpointError` (the default): the global error 2-norm at the final
#'   time only; `sol$global_error[n]` holds it and the earlier entries are
#'   zero. One adjoint solve per sample.
#' * `TrajectoryError`: the global error 2-norm at every saved time
#'   (`sol$global_error[1]` is 0). One adjoint solve per saved time per
#'   sample; the returned solution is dense and saved at every step.
#'
#' `GlobalAdjointControl` selects how [GlobalAdjoint()] meets `gtol`:
#' * `ToleranceRefinement` (the default) tightens the local `abstol` and
#'   `reltol` and solves again adaptively until the estimated endpoint
#'   error is at most `gtol`.
#' * `StepGridRefinement` tightens the tolerances in the same way to find a
#'   step grid that meets `gtol`, then returns the solution of a fixed step
#'   solve on exactly that grid (`adaptive = FALSE`), a reproducible fixed
#'   step schedule.
#'
#' @format Enumeration instances (classed character strings) and lists of
#'   them.
#' @references P. E. Zadunaisky, On the estimation of errors propagated in
#'   the numerical integration of ordinary differential equations,
#'   Numerische Mathematik 27, 1976. L. F. Shampine, Global error estimation
#'   with one-step methods, Computers and Mathematics with Applications 12A,
#'   1986.
#' @examples
#' GlobalErrorEquation$ErrorTransport
#' identical(GlobalAdjoint(Tsit5())$scope, EndpointError)
#' names(GlobalAdjointControl)
#' @name global_error_options
NULL

#' @rdname global_error_options
#' @export
DefectCorrection <- .gde_enum("GlobalErrorEquation", "DefectCorrection", 0L)
#' @rdname global_error_options
#' @export
ErrorTransport <- .gde_enum("GlobalErrorEquation", "ErrorTransport", 1L)
#' @rdname global_error_options
#' @export
GlobalErrorEquation <- list(DefectCorrection = DefectCorrection, ErrorTransport = ErrorTransport)

#' @rdname global_error_options
#' @export
InterpolatingMode <- .gde_enum("GlobalErrorMode", "InterpolatingMode", 0L)
#' @rdname global_error_options
#' @export
SimultaneousMode <- .gde_enum("GlobalErrorMode", "SimultaneousMode", 1L)
#' @rdname global_error_options
#' @export
GlobalErrorMode <- list(InterpolatingMode = InterpolatingMode, SimultaneousMode = SimultaneousMode)

#' @rdname global_error_options
#' @export
EndpointError <- .gde_enum("GlobalAdjointScope", "EndpointError", 0L)
#' @rdname global_error_options
#' @export
TrajectoryError <- .gde_enum("GlobalAdjointScope", "TrajectoryError", 1L)
#' @rdname global_error_options
#' @export
GlobalAdjointScope <- list(EndpointError = EndpointError, TrajectoryError = TrajectoryError)

#' @rdname global_error_options
#' @export
ToleranceRefinement <- .gde_enum("GlobalAdjointControl", "ToleranceRefinement", 0L)
#' @rdname global_error_options
#' @export
StepGridRefinement <- .gde_enum("GlobalAdjointControl", "StepGridRefinement", 1L)
#' @rdname global_error_options
#' @export
GlobalAdjointControl <- list(ToleranceRefinement = ToleranceRefinement, StepGridRefinement = StepGridRefinement)

# ---------------------------------------------------------------------------
# First derivatives of the dense output: the Val{1} _ode_interpolant methods
# of the ConstantCache (out of place) paths, by the interpolant a step
# carries. P' = d/dt = (1 / dt) d/dtheta; @evalpoly is a fused Horner.

.gde_horner <- function(x, coefs) {
  h <- coefs[length(coefs)]
  for (q in rev(seq_len(length(coefs) - 1L))) h <- .fma(x, h, coefs[q])
  h
}

# OrdinaryDiffEqCore hermite_interpolant, Val{1}
.gde_hermite_d <- function(theta, dt, y0, y1, k) {
  k1 <- k[[1L]]; k2 <- k[[2L]]
  k1 + theta * (-4 * dt * k1 - 2 * dt * k2 - 6 * y0 +
                  theta * (3 * dt * k1 + 3 * dt * k2 + 6 * y0 - 6 * y1) + 6 * y1) / dt
}

# .hermite_dae_interp: the Hermite derivative, the secant where the slopes
# are missing (algebraic components)
.gde_hermite_dae_d <- function(theta, dt, y0, y1, k) {
  out <- .gde_hermite_d(theta, dt, y0, y1, k)
  ifelse(is.na(out), (y1 - y0) / dt, out)
}

# .bdf_interp: the derivative of the Lagrange polynomial through the
# history values. OrdinaryDiffEqBDF resamples that polynomial at fixed
# Chebyshev nodes (_resample_at_chebyshev!, bdf_utils.jl) and differentiates
# the Lagrange basis on them (_fbdf_interp, Val{1}, bdf_interpolants.jl).
# It is the same polynomial, but the slope over a short early step is
# bound by cancellation, so Julia's order of evaluation is kept.
.gde_bdf_d <- function(theta, dt, y0, y1, k) {
  n <- length(k$thetas)
  nodes <- (1 + cospi((2 * seq_len(n) - 1) / (2 * n))) / 2
  kc <- lapply(nodes, function(xi) .eval_lagrange(xi, k$thetas, k$values, n))
  dbasis <- function(i) {
    result <- 0
    for (m in seq_len(n)) {
      if (m == i) next
      term <- 1 / (nodes[i] - nodes[m])
      for (l in seq_len(n)) if (l != i && l != m) term <- term * ((theta - nodes[l]) / (nodes[i] - nodes[l]))
      result <- result + term
    }
    result
  }
  invdt <- 1 / dt
  out <- (dbasis(1L) * invdt) * kc[[1L]]
  for (j in seq_len(n)[-1L]) out <- .fma(dbasis(j) * invdt, kc[[j]], out)
  out
}

# OrdinaryDiffEqTsit5 @tsit5pre1
.gde_tsit5_d <- function(theta, dt, y0, y1, k) {
  r <- .TSIT5$r
  b1 <- .gde_horner(theta, r[[1L]] * c(1, 2, 3, 4))
  b <- vapply(2:7, function(i) theta * .gde_horner(theta, r[[i]] * c(2, 3, 4)), 0)
  k[[1L]] * b1 + k[[2L]] * b[1L] + k[[3L]] * b[2L] + k[[4L]] * b[3L] + k[[5L]] * b[4L] +
    k[[6L]] * b[5L] + k[[7L]] * b[6L]
}

# OrdinaryDiffEqLowOrderRK @dp5pre1
.gde_dp5_d <- function(theta, dt, y0, y1, k) {
  b20 <- .gde_horner(theta, c(1, -2))
  b30 <- theta * .gde_horner(theta, c(2, -3))
  b40 <- theta * .gde_horner(theta, c(2, -6, 4))
  k[[1L]] + k[[2L]] * b20 + k[[3L]] * b30 + k[[4L]] * b40
}

# OrdinaryDiffEqVerner @vern6pre1 to @vern9pre1: stage 1 weighs
# evalpoly(theta, r1, 2 r2, ...), the others theta evalpoly(theta, 2 r2, ...).
.gde_verner_d <- function(theta, dt, y0, y1, k) {
  # the value interpolant computes the lazy extra stages
  if (!k$lazy$done) .verner_interp(theta, dt, y0, y1, k)
  R <- k$tab$R
  K <- k$lazy$K
  acc <- NULL
  for (j in which(rowSums(R != 0) > 0)) {
    m <- if (j == 1L) 1L else 2L
    q <- m:max(which(R[j, ] != 0))
    h <- .gde_horner(theta, q * R[j, q])
    term <- K[, j] * (if (m == 1L) h else theta * h)
    acc <- if (is.null(acc)) term else acc + term
  }
  acc
}

# OrdinaryDiffEqRosenbrock @rosenbrock2332pre1, for the interpolant
# .rosen_interp23(d)
.gde_rosenbrock23_d <- function(d) {
  force(d)
  function(theta, dt, y0, y1, k)
    ((1 - 2 * theta) / (1 - 2 * d)) * k[[1L]] + ((2 * theta - 2 * d) / (1 - 2 * d)) * k[[2L]]
}

# The derivative of .rosen_interp(order): the Rodas H matrix dense output of
# interp_order 2, 3 and 4, and the Hermite fallback (order -1), linear in the
# algebraic variables.
.gde_rosen_d <- function(order) {
  force(order)
  function(theta, dt, y0, y1, k) {
    if (order > 0L) return(.gde_rodas_d(order, theta, dt, y0, y1, k))
    m <- attr(k, "dvars")
    if (is.null(m)) return(.gde_hermite_d(theta, dt, y0, y1, k))
    lin <- (y1 - y0) / dt
    if (identical(m, FALSE)) return(lin)
    ifelse(m, .gde_hermite_d(theta, dt, y0, y1, k), lin)
  }
}

.gde_rodas_d <- function(order, theta, dt, y0, y1, k) {
  inner <- switch(order - 1L,
    -2 * k[[1L]] + 2 * k[[2L]] - 3 * k[[2L]] * theta,
    -2 * k[[1L]] + 2 * k[[2L]] + theta * (-3 * k[[2L]] + 3 * k[[3L]] - 4 * theta * k[[3L]]),
    -2 * k[[1L]] + 2 * k[[2L]] +
      theta * (-3 * k[[2L]] + 3 * k[[3L]] + theta * (-4 * k[[3L]] + 4 * k[[4L]] - 5 * theta * k[[4L]])))
  (k[[1L]] + theta * inner - y0 + y1) / dt
}

# The derivative of a step's interpolant, or NULL when there is none (the
# linear and constant fallbacks of Julia's _has_dense_derivative, and the
# interpolants not ported here).
.gde_dinterp_of <- function(interp) {
  if (!is.function(interp)) return(NULL)
  if (identical(interp, .tsit5_interp)) return(.gde_tsit5_d)
  if (identical(interp, .hermite_interp) || identical(interp, .hermite_interp_fused)) return(.gde_hermite_d)
  if (identical(interp, .dp5_interp)) return(.gde_dp5_d)
  if (identical(interp, .verner_interp)) return(.gde_verner_d)
  # the Rosenbrock interpolants are closures over their order or d
  if (identical(body(interp), body(.rosen_interp(2L)))) return(.gde_rosen_d(environment(interp)$order))
  if (identical(body(interp), body(.rosen_interp23(0)))) return(.gde_rosenbrock23_d(environment(interp)$d))
  if (identical(interp, .hermite_dae_interp)) return(.gde_hermite_dae_d)
  if (identical(interp, .bdf_interp)) return(.gde_bdf_d)
  NULL
}

.gde_alg_has_dense_derivative <- function(alg) {
  algs <- if (inherits(alg, "CompositeAlgorithm")) alg$algs else list(alg)
  all(vapply(algs, function(a) !is.null(.gde_dinterp_of(a$interp)), TRUE))
}

.gde_has_dense_derivative <- function(sol) {
  env <- .sol_env(sol)
  steps <- env$steps
  isTRUE(env$dense) && length(steps) == length(env$t) - 1L &&
    all(vapply(steps, function(s) !is.null(.gde_dinterp_of(s$interp)), TRUE))
}

# The dense output of a stored solution as OrdinaryDiffEqCore's
# ode_interpolation: `right` is continuity = :right (the step starting at a
# node), otherwise the step ending at it; theta and dt come from the saved
# times. Returns the value and the first derivative at t.
.gde_sol_dense <- function(sol, right) {
  env <- .sol_env(sol)
  ts <- env$t
  U <- unname(env$U)
  steps <- env$steps
  n <- length(ts)
  locate <- function(t) {
    if (right) {
      im <- max(1L, findInterval(t, ts))
      ip <- if (im < n) im + 1L else im
    } else {
      ip <- min(n, max(2L, findInterval(t, ts, left.open = TRUE) + 1L))
      im <- if (ip > 1L) ip - 1L else ip
    }
    c(im, ip)
  }
  eval_at <- function(t, deriv) {
    ix <- locate(t)
    im <- ix[1L]; ip <- ix[2L]
    if (im == ip) {
      if (!deriv) return(U[im, ])
      # at the end with continuity = :right Julia differentiates the last
      # step with dt = 0, which is NaN for the Hermite and Rodas forms; the
      # derivative is taken from the last step (its value at theta = 1)
      im <- n - 1L
    }
    dt <- ts[ip] - ts[im]
    theta <- if (dt == 0) 1 else (t - ts[im]) / dt
    st <- steps[[im]]
    fn <- if (deriv) .gde_dinterp_of(st$interp) else st$interp
    as.numeric(fn(theta, dt, U[im, ], U[ip, ], st$k))
  }
  list(value = function(t) eval_at(t, FALSE), deriv = function(t) eval_at(t, TRUE))
}

# The live dense output of an integrator inside its last step (Julia's
# integrator(t) and integrator(t, Val{1})).
.gde_integ_dense <- function(integ) {
  st <- integ$last_step
  dfn <- .gde_dinterp_of(st$interp)
  theta_of <- function(t) (t - st$t0) / st$dt
  list(value = function(t) as.numeric(st$interp(theta_of(t), st$dt, st$y0, st$y1, st$k)),
       deriv = function(t) as.numeric(dfn(theta_of(t), st$dt, st$y0, st$y1, st$k)))
}

# ---------------------------------------------------------------------------
# Validation

.gde_tol_ok <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0

.gde_validate_tolerances <- function(abstol, reltol, name) {
  if (!.gde_tol_ok(abstol)) stop(name, "_abstol must be a nonnegative finite real number", call. = FALSE)
  if (!.gde_tol_ok(reltol)) stop(name, "_reltol must be a nonnegative finite real number", call. = FALSE)
  if (abstol == 0 && reltol == 0) stop(name, "_abstol and ", name, "_reltol cannot both be zero", call. = FALSE)
  invisible()
}

.gde_positive_finite <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0

# _has_callbacks: an empty CallbackSet counts as none
.gde_has_callbacks <- function(cb) {
  !is.null(cb) && !(inherits(cb, "CallbackSet") && !length(cb$continuous) && !length(cb$discrete))
}

.gde_validate_state <- function(prob, name) {
  u0 <- prob$u0
  if (!is.double(u0) || !is.null(dim(u0))) stop(name, " requires a real floating-point vector state", call. = FALSE)
  if (!length(u0)) stop(name, " requires a nonempty state", call. = FALSE)
  if (!(prob$tspan[1L] < prob$tspan[2L])) stop(name, " requires a forward-time ODEProblem", call. = FALSE)
  if (!.is_identity_mass(prob$f$mass_matrix)) stop(name, " currently requires the standard mass matrix", call. = FALSE)
  invisible()
}

.gde_validate_problem <- function(prob, name) {
  .gde_validate_state(prob, name)
  if (.gde_has_callbacks(prob$kwargs$callback)) stop(name, " does not currently support callbacks", call. = FALSE)
  invisible()
}

.gde_validate_solution <- function(sol, name) {
  if (!successful_retcode(sol)) stop("the forward solve failed with retcode ", sol$retcode, call. = FALSE)
  if (length(sol$t) < 2L) stop("the forward solution must save at least its two endpoints", call. = FALSE)
  if (!.gde_has_dense_derivative(sol))
    stop(name, " requires a solver with dense (continuous) output that provides a first-derivative ",
         "interpolation; the wrapped solver produced none. Use a dense solver (e.g. Tsit5, Vern7) and do ",
         "not disable `dense`.", call. = FALSE)
  invisible()
}

.gde_companion_options <- function(gtol, maxiters, safety, companion_abstol, companion_reltol) {
  if (!(is.numeric(maxiters) && length(maxiters) == 1L && maxiters == round(maxiters) && maxiters > 0))
    stop("maxiters must be a positive integer", call. = FALSE)
  if (!(is.numeric(safety) && length(safety) == 1L && is.finite(safety) && safety > 0 && safety < 1))
    stop("safety must be between zero and one", call. = FALSE)
  if (!(is.null(gtol) || .gde_positive_finite(gtol))) stop("gtol must be a positive finite real number", call. = FALSE)
  .gde_validate_tolerances(companion_abstol, companion_reltol, "companion")
  list(maxiters = as.integer(maxiters), safety = safety, companion_abstol = companion_abstol,
       companion_reltol = companion_reltol)
}

# ---------------------------------------------------------------------------
# Companion right-hand sides: make_rhs(dense, prob) returns rhs(eps, t) for
# the dense output `dense` (value and derivative).

.gde_oop_rhs <- function(prob) {
  f <- prob$f$f
  function(u, p, t) as.numeric(f(u, p, t))
}

# J(u) v: the exact Jacobian of the problem (jac_fn, or symbolic for a
# problem from a System) for AutoForwardDiff, otherwise a finite difference
# along v (central, or forward for AutoFiniteDiff("forward")).
.gde_jvp <- function(prob, autodiff) {
  foop <- .gde_oop_rhs(prob)
  p <- prob$p
  ad <- if (inherits(autodiff, "AutoSparse")) autodiff$dense_ad else autodiff
  fd <- inherits(ad, "AutoFiniteDiff")
  if (!fd && (is.function(prob$f$jac) || (inherits(prob$sys, "System") && !is.null(prob$f$sys)))) {
    J <- .sens_fns(prob)$jac
    return(function(u, v, t) as.numeric(J(u, p, t) %*% v))
  }
  forward <- fd && identical(ad$fdtype, "forward")
  function(u, v, t) {
    nv <- sqrt(sum(v * v))
    if (nv == 0) return(numeric(length(u)))
    w <- v / nv
    nu <- max(1, sqrt(sum(u * u)))
    if (forward) {
      h <- sqrt(.Machine$double.eps) * nu
      return(nv * (foop(u + h * w, p, t) - foop(u, p, t)) / h)
    }
    h <- .Machine$double.eps^(1 / 3) * nu
    nv * (foop(u + h * w, p, t) - foop(u - h * w, p, t)) / (2 * h)
  }
}

.gde_rhs_builder <- function(equation, autodiff) {
  if (identical(equation, DefectCorrection)) {
    return(function(dense, prob) {
      foop <- .gde_oop_rhs(prob)
      p <- prob$p
      function(e, t) foop(dense$value(t) + e, p, t) - dense$deriv(t)
    })
  }
  function(dense, prob) {
    foop <- .gde_oop_rhs(prob)
    jvp <- .gde_jvp(prob, autodiff)
    p <- prob$p
    function(e, t) {
      u <- dense$value(t)
      jvp(u, e, t) + (foop(u, p, t) - dense$deriv(t))
    }
  }
}

# Solve eps' = rhs(eps, t), eps(t0) = e0 over [t0, t1] and return eps(t1).
.gde_advance <- function(rhs, e0, t0, t1, alg, abstol, reltol, tstops = NULL) {
  cprob <- ODEProblem(function(u, p, t) rhs(u, t), e0, c(t0, t1))
  csol <- solve(cprob, alg, abstol = abstol, reltol = reltol, tstops = tstops, dense = FALSE,
                save_everystep = FALSE, save_start = FALSE, save_end = TRUE)
  if (!successful_retcode(csol))
    stop("the companion error solve failed with retcode ", csol$retcode, call. = FALSE)
  u <- csol$u
  u[[length(u)]]
}

.GDE_DENSE_OPTS <- list(dense = TRUE, save_everystep = TRUE, save_start = TRUE, save_end = TRUE,
                        saveat = numeric(), save_idxs = NULL)

# InterpolatingMode: a dense forward solve, then the companion against its
# interpolant (right continuous at the step nodes, which the companion
# solver lands on). Returns the endpoint error 2-norm.
.gde_estimate_interpolating <- function(make_rhs, name, prob, inner_alg, companion_alg, opts, abstol, reltol,
                                        companion_abstol, companion_reltol) {
  if (.gde_has_callbacks(opts$callback)) stop(name, " does not currently support callbacks", call. = FALSE)
  .gde_validate_problem(prob, name)
  sopts <- utils::modifyList(opts, .GDE_DENSE_OPTS)
  sol <- do.call(solve, c(list(prob, inner_alg), list(abstol = abstol, reltol = reltol), sopts))
  .gde_validate_solution(sol, name)
  rhs <- make_rhs(.gde_sol_dense(sol, right = TRUE), sol$prob)
  e <- .gde_advance(rhs, numeric(length(prob$u0)), prob$tspan[1L], prob$tspan[2L], companion_alg,
                    companion_abstol, companion_reltol, tstops = sol$t)
  sqrt(sum(e * e))
}

# SimultaneousMode: step the forward integrator and advance the companion
# across each completed step with the step's live interpolant.
.gde_estimate_simultaneous <- function(make_rhs, name, prob, inner_alg, companion_alg, opts, abstol, reltol,
                                       companion_abstol, companion_reltol) {
  if (.gde_has_callbacks(opts$callback)) stop(name, " does not currently support callbacks", call. = FALSE)
  .gde_validate_problem(prob, name)
  iopts <- utils::modifyList(list(abstol = abstol, reltol = reltol, dense = TRUE, save_everystep = FALSE,
                                  save_start = TRUE, save_end = TRUE), opts)
  integ <- do.call(init, c(list(prob, inner_alg), iopts))
  I <- .integ_I(integ)
  if (!isTRUE(I$dense) || !.gde_alg_has_dense_derivative(I$alg))
    stop(name, " with SimultaneousMode requires a solver with dense (continuous) output that provides a ",
         "first-derivative interpolation; the wrapped solver produced none. Use a dense solver (e.g. Tsit5, ",
         "Vern7) and do not disable `dense`.", call. = FALSE)
  e <- numeric(length(prob$u0))
  # the driver ends on a call that takes no step (t unchanged)
  seen <- integ$t
  while (!I$done) {
    step_integrator(integ)
    tp <- integ$tprev
    tt <- integ$t
    if (!(tt > tp) || tt == seen) next
    seen <- tt
    e <- .gde_advance(make_rhs(.gde_integ_dense(integ), prob), e, tp, tt, companion_alg, companion_abstol,
                      companion_reltol)
  }
  if (!successful_retcode(integ$retcode))
    stop("the forward solve failed with retcode ", integ$retcode, call. = FALSE)
  sqrt(sum(e * e))
}

# The endpoint error refinement loop (_refine_to_gtol): estimate(la, lr)
# returns the estimate of a trial solve; finish(la, lr) the final solution.
.gde_refine <- function(estimate, finish, gtol, options, abstol, reltol, what = "global error estimate") {
  la <- abstol
  lr <- reltol
  last <- Inf
  for (it in seq_len(options$maxiters)) {
    last <- estimate(la, lr)
    if (!is.finite(last)) stop("the ", what, " is not finite", call. = FALSE)
    if (last <= gtol) return(finish(la, lr, last))
    scale <- min(0.5, options$safety * gtol / last)
    la <- la * scale
    lr <- lr * scale
    if (la == 0 && lr == 0) stop("local tolerances underflowed during global error control", call. = FALSE)
  }
  stop("failed to meet gtol=", .gde_jl_float(gtol), " after ", options$maxiters, " iterations; last estimated ",
       "global error was ", .gde_jl_float(last), call. = FALSE)
}

# Julia's string of a Float64: the shortest digits that read back, written
# in scientific form below 1e-4 and from 1e16 ("1.0e-9", "0.0026798526949479")
.gde_jl_float <- function(x) {
  if (!is.finite(x)) return(if (is.na(x)) "NaN" else if (x > 0) "Inf" else "-Inf")
  if (x == 0) return("0.0")
  for (dg in 1:17) {
    s <- formatC(abs(x), digits = dg - 1L, format = "e")
    if (as.numeric(s) == abs(x)) break
  }
  d <- sub("0+$", "", gsub(".", "", sub("e.*", "", s), fixed = TRUE))
  if (!nzchar(d)) d <- "0"
  e <- as.integer(sub(".*e", "", s))
  out <- if (e < -4L || e >= 16L) {
    paste0(substr(d, 1L, 1L), ".", if (nchar(d) > 1L) substring(d, 2L) else "0", "e", e)
  } else if (e < 0L) {
    paste0("0.", strrep("0", -e - 1L), d)
  } else {
    if (nchar(d) < e + 1L) d <- paste0(d, strrep("0", e + 1L - nchar(d)))
    frac <- substring(d, e + 2L)
    paste0(substr(d, 1L, e + 1L), ".", if (nzchar(frac)) frac else "0")
  }
  if (x < 0) paste0("-", out) else out
}

# ---------------------------------------------------------------------------
# GlobalErrorEstimation

#' Global error estimation with a companion equation
#'
#' `GlobalErrorEstimation(alg)` wraps an ODE algorithm with global error
#' estimation and control based on a defect driven companion equation, as
#' GlobalDiffEq.jl. After a dense solve with interpolant `P(t)`, the global
#' error `eps(t) ~ y(t) - P(t)` is estimated by integrating a companion ODE
#' driven by the defect `d(t) = f(P(t), p, t) - P'(t)` of the dense output,
#' from `eps(t0) = 0`.
#'
#' `equation` ([GlobalErrorEquation]) selects the companion ODE:
#' [DefectCorrection] (the default) or its linearization [ErrorTransport],
#' whose Jacobian vector product uses `autodiff`. `mode`
#' ([GlobalErrorMode]) selects how it is integrated: [InterpolatingMode]
#' (the default) or [SimultaneousMode].
#'
#' `global_error_estimate(prob, alg)` solves `prob` and returns the
#' estimated 2-norm of the global error at the final time (`abstol` and
#' `reltol` control the solve, `companion_abstol` and `companion_reltol` the
#' companion solve). `solve(prob, GlobalErrorEstimation(alg, gtol = ...))`
#' tightens the local tolerances (which default to `gtol`) until the
#' estimated endpoint error is at most `gtol` and returns the solution of
#' `alg` at those tolerances.
#'
#' The solver must have a dense output with a first derivative: Tsit5, DP5,
#' the Verner methods, Rosenbrock23, the Rodas methods, FBDF and the methods
#' with cubic Hermite dense output (ImplicitEuler and Trapezoid among them).
#' The low storage and SSP methods, whose dense output is linear in this
#' package (OrdinaryDiffEq has Hermite or quadratic dense output for them),
#' are rejected. The companion solver defaults to `alg`; a method that only
#' evaluates `f` at the end of each step (ImplicitEuler) samples the defect
#' where it vanishes, so give such a method an explicit `companion_alg`.
#' Forward time ODEs with the standard mass matrix, real vector states and
#' no callbacks are supported.
#'
#' @param alg The ODE algorithm.
#' @param equation [DefectCorrection] or [ErrorTransport].
#' @param mode [InterpolatingMode] or [SimultaneousMode].
#' @param companion_alg The algorithm of the companion solve.
#' @param autodiff How `ErrorTransport` applies the Jacobian:
#'   [AutoForwardDiff()] (the exact Jacobian of the problem when it has one,
#'   `jac_fn` or symbolic, else a central difference along the error) or
#'   [AutoFiniteDiff()].
#' @param gtol The requested absolute endpoint error, needed by `solve()`.
#' @param maxiters The maximum number of tolerance refinements.
#' @param safety The safety factor of each tolerance reduction.
#' @param companion_abstol,companion_reltol Tolerances of the companion
#'   solve.
#' @return An algorithm object.
#' @references P. E. Zadunaisky, Numerische Mathematik 27, 1976. J. R.
#'   Dormand, R. R. Duckers and P. J. Prince, Global error estimation with
#'   Runge-Kutta methods, IMA Journal of Numerical Analysis 4, 1984. J. Lang
#'   and J. Verwer, On global error estimation and control for initial
#'   value problems, SIAM Journal on Scientific Computing 29, 2007.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
#' alg <- GlobalErrorEstimation(Tsit5())
#' global_error_estimate(prob, alg, abstol = 1e-4, reltol = 1e-4)
#' sol <- solve(prob, Tsit5(), abstol = 1e-4, reltol = 1e-4)
#' abs(sol$u[[length(sol$t)]] - exp(-1))
#' sol <- solve(prob, GlobalErrorEstimation(Tsit5(), gtol = 1e-7))
#' abs(sol$u[[length(sol$t)]] - exp(-1))
#' @export
GlobalErrorEstimation <- function(alg, equation = DefectCorrection, mode = InterpolatingMode, companion_alg = alg,
                                  autodiff = AutoForwardDiff(), gtol = NULL, maxiters = 6, safety = 0.8,
                                  companion_abstol = if (is.null(gtol)) 1e-10 else min(gtol / 100, 1e-10),
                                  companion_reltol = if (is.null(gtol)) 1e-8 else min(gtol / 100, 1e-8)) {
  if (!inherits(alg, "ODEAlgorithm")) stop("GlobalErrorEstimation wraps an ODE algorithm", call. = FALSE)
  if (!(identical(equation, DefectCorrection) || identical(equation, ErrorTransport)))
    stop("`equation` must be DefectCorrection or ErrorTransport", call. = FALSE)
  if (!(identical(mode, InterpolatingMode) || identical(mode, SimultaneousMode)))
    stop("`mode` must be InterpolatingMode or SimultaneousMode", call. = FALSE)
  options <- .gde_companion_options(gtol, maxiters, safety, companion_abstol, companion_reltol)
  structure(list(name = "GlobalErrorEstimation", alg = alg, companion_alg = companion_alg, equation = equation,
                 mode = mode, autodiff = autodiff, gtol = gtol, options = options, order = alg$order,
                 adaptive = TRUE, solve_override = .gde_solve),
            class = c("GlobalErrorEstimation", "ODEAlgorithm"))
}

#' @rdname global_error_estimate
#' @export
global_error_estimate.ODEProblem <- function(x, alg, ..., abstol = 1e-6, reltol = 1e-3,
                                             companion_abstol = alg$options$companion_abstol,
                                             companion_reltol = alg$options$companion_reltol) {
  if (!inherits(alg, "GlobalErrorEstimation"))
    stop("global_error_estimate(prob, alg) needs a GlobalErrorEstimation algorithm", call. = FALSE)
  make_rhs <- .gde_rhs_builder(alg$equation, alg$autodiff)
  estimate <- if (identical(alg$mode, SimultaneousMode)) .gde_estimate_simultaneous else .gde_estimate_interpolating
  estimate(make_rhs, "GlobalErrorEstimation", x, alg$alg, alg$companion_alg, list(...), abstol, reltol,
           companion_abstol, companion_reltol)
}

# SciMLBase.__solve for GlobalErrorEstimation
.gde_solve <- function(prob, alg, opts) {
  if (is.null(alg$gtol)) stop("GlobalErrorEstimation requires a positive `gtol` constructor keyword", call. = FALSE)
  abstol <- opts$abstol %||% alg$gtol
  reltol <- opts$reltol %||% alg$gtol
  opts$abstol <- NULL
  opts$reltol <- NULL
  .gde_validate_tolerances(abstol, reltol, "local")
  if (.gde_has_callbacks(opts$callback))
    stop("GlobalErrorEstimation does not currently support callbacks", call. = FALSE)
  .gde_refine(function(la, lr) do.call(global_error_estimate, c(list(prob, alg), opts, list(abstol = la, reltol = lr))),
              function(la, lr, est) do.call(solve, c(list(prob, alg$alg), opts, list(abstol = la, reltol = lr))),
              alg$gtol, alg$options, abstol, reltol)
}
