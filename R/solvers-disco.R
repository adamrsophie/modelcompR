# Discontinuity detection of OrdinaryDiffEqCore (disco.jl, the
# zero_func_struct of solve.jl and handle_disco_accept! /
# handle_disco_reject! of integrators/controllers.jl), and the Predictor
# enumeration of OrdinaryDiffEqCore.jl.
#
# With `solve(..., discontinuity_detection = TRUE)` a rejected step looks for
# a sign change of the continuous callbacks (those with maybe_discontinuity)
# between the start of the step and the rejected end point. The root of the
# condition on the rejected step's interpolant, found by ModAB on theta in
# [0, 1], gives the next step size; the step after the accepted one is at
# most a quarter of the distance to the end of the rejected step.

#' Discontinuity detection
#'
#' `set_discontinuity(integrator)` returns the step size with which the
#' current step of `integrator` should be tried again so that it ends near
#' the first root of a continuous callback inside the step, or `-1` when no
#' callback changes sign in the step. Only callbacks with
#' `maybe_discontinuity = TRUE` (the default of [ContinuousCallback()]) are
#' considered. With the root at the fraction `theta` of the step, the new
#' step is `(0.5 + 1.2372 (theta - 0.5) - 1.7487 (theta - 0.5)^3) dt`, a
#' smooth version of `(0.5 + 0.4 sin(pi (theta - 0.5))) dt`.
#'
#' `solve(prob, alg, discontinuity_detection = TRUE)` applies it after every
#' rejected step of the integral, PI, PID and predictive controllers (as the
#' `discontinuity_detection` option of the controllers of OrdinaryDiffEq), so
#' that a step rejected because it jumped over a discontinuity of the right
#' hand side retries just before the discontinuity instead of shrinking by
#' the controller's factor.
#'
#' @param integrator An integrator from [init()]. Its current step is the
#'   last step taken (or attempted, inside the solver).
#' @return The new step size, or `-1`.
#' @examples
#' f <- function(u, p, t) c(if (u[2] < 0.5) u[1] else -u[1] * u[1], 1)
#' cb <- ContinuousCallback(function(u, t, integ) u[2] - 0.5, function(integ) NULL)
#' prob <- ODEProblem(f, c(1, 0), c(0, 1))
#' sol <- solve(prob, Tsit5(), callback = cb, discontinuity_detection = TRUE)
#' sol$stats$nreject
#' @export
set_discontinuity <- function(integrator) {
  I <- integrator$.I
  if (!is.environment(I)) stop("set_discontinuity() needs an integrator from init()", call. = FALSE)
  step <- I$disco_attempt %||% integrator$last_step
  if (is.null(step)) return(-1)
  fd <- .find_discontinuity(I$cbs, integrator, step)
  .disco_dt(fd$theta, step$dt)
}

# set_discontinuity: the retry step for the root at theta (1: none).
.disco_dt <- function(theta, dt) {
  if (theta < 1) {
    x <- theta - 0.5
    return((0.5 + 1.2372 * x - 1.7487 * (x * x * x)) * dt)
  }
  -1
}

# find_discontinuity: the smallest root theta in (0, 1) of the continuous
# callbacks on the step `step` (list(t0, dt, y0, y1, k, interp)), 1 if none.
# The zero function returns the end point values at theta = 0 and 1 and
# otherwise evaluates the condition on the interpolant. OrdinaryDiffEq
# brackets a later callback on [0, theta] with theta the root found so far
# without checking the sign at theta, so ModAB returns its left end 0 when
# the later crossing lies beyond theta; here such a callback is skipped (see
# PORTING.md).
.find_discontinuity <- function(cbs, integ, step) {
  interp <- .inplace_interp(step$interp)
  theta_bp <- 1
  found <- FALSE
  for (cb in cbs$continuous) {
    if (isFALSE(cb$maybe_discontinuity)) next
    lo <- .as_num(cb$condition(step$y0, step$t0, integ))
    hi <- .as_num(cb$condition(step$y1, step$t0 + step$dt, integ))
    for (j in seq_along(lo)) {
      if (!isTRUE(lo[j] * hi[j] < 0)) next
      g <- function(th) {
        if (th == 0) return(lo[j])
        if (th == 1) return(hi[j])
        uth <- interp(th, step$dt, step$y0, step$y1, step$k)
        .as_num(cb$condition(uth, step$t0 + th * step$dt, integ))[j]
      }
      if (theta_bp < 1 && !isTRUE(lo[j] * g(theta_bp) < 0)) next
      r <- .br_modab(g, 0, theta_bp, lo[j], g(theta_bp), .nl_default_tol(), 1000L)
      if (!is.nan(r$u) && r$u < theta_bp) {
        theta_bp <- r$u
        found <- TRUE
      }
    }
  }
  list(theta = theta_bp, found = found)
}

# The zero function interpolates in place (ode_interpolant!), where the
# multiply adds of the broadcast may be fused. The Hermite form is fused
# throughout; for the method specific interpolants LLVM contracts some of
# the muladd calls and not others, so those keep the out of place form.
.inplace_interp <- function(interp) {
  if (identical(interp, .hermite_interp)) return(.hermite_interp_fused)
  if (identical(interp, .lsrk_interp)) return(.lsrk_interp_fused)
  interp
}

#' Initial guess strategies of implicit Runge-Kutta stages
#'
#' `Predictor` lists the per stage initial guesses of the Newton iterations
#' of implicit Runge-Kutta methods (the `predictor` option of the SDIRK
#' methods): `Trivial` (a zero increment), `Linear` (linear extrapolation
#' with the derivative at the start of the step), `MaxOrder` (the full order
#' interpolant of the previous step), `VariableOrder` (the interpolant with
#' the order reduced as the stage extrapolates farther), `CutoffOrder` (full
#' order below a cutoff abscissa, else order one), `CopyPrev` (the previous
#' stage's derivative), `StageExtrap` (extrapolation of recent stage
#' derivatives) and `Tableau` (the guess derived from the tableau).
#'
#' @format A named list of strings.
#' @examples
#' Predictor$Linear
#' @export
Predictor <- list(Trivial = "Trivial", Linear = "Linear", MaxOrder = "MaxOrder",
                  VariableOrder = "VariableOrder", CutoffOrder = "CutoffOrder",
                  CopyPrev = "CopyPrev", StageExtrap = "StageExtrap", Tableau = "Tableau")
