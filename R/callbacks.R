# Numerical callbacks for the integrators, after DiffEqBase and
# DiffEqCallbacks. Symbolic events (SymbolicContinuousCallback and
# SymbolicDiscreteCallback) compile to these.
#
# An affect function receives the integrator, an environment with fields
# t, u, p, tprev, uprev, terminate, and changes them in place:
#   affect = function(integ) integ$u[1] <- -integ$u[1]

#' Numerical callbacks
#'
#' `ContinuousCallback(condition, affect)` fires when `condition(u, t, integ)`
#' crosses zero; the time is located by root finding on the solver's
#' interpolant. `affect` runs on upcrossings and `affect_neg` on
#' downcrossings (by default the same function). A vector valued condition
#' fires per component and the affect receives the component index as its
#' second argument. `DiscreteCallback(condition, affect)` fires after a step
#' when `condition(u, t, integ)` is `TRUE`. `PresetTimeCallback(times,
#' affect)` fires at given times, and `PeriodicCallback(affect, dt)` every
#' `dt`. `CallbackSet()` combines callbacks. Affects modify the integrator
#' `integ` in place: `integ$u`, `integ$p`, and `terminate(integ)` stops the
#' integration.
#'
#' @param condition A function `(u, t, integ)`.
#' @param affect,affect_neg Functions `(integ)` (or `(integ, idx)` for vector
#'   conditions); `affect_neg = NULL` ignores downcrossings.
#' @param initialize A function `(integ)` called at the start.
#' @param finalize A function `(integ)` called at the end.
#' @param rootfind How the event time is located: `"left"` (or `TRUE`, the
#'   default, as SciMLBase's LeftRootFind) applies the affect just before the
#'   crossing, `"right"` just after it, and `"none"` (or `FALSE`) at the end
#'   of the step.
#' @param save_positions Save the state before and after the affect.
#' @param dtrelax After a continuous event the next step starts from
#'   `dtrelax` times the step that contained the event (DiffEqBase's
#'   `apply_callback!`).
#' @param maybe_discontinuity Whether the event may introduce a
#'   discontinuity: only these callbacks take part in the discontinuity
#'   detection of `solve(..., discontinuity_detection = TRUE)`.
#' @param times Event times.
#' @param dt Period.
#' @param ... Callbacks.
#' @return A callback object.
#' @export
ContinuousCallback <- function(condition, affect, affect_neg = affect, initialize = NULL,
                               finalize = NULL, rootfind = "left", save_positions = c(TRUE, TRUE),
                               dtrelax = 1, maybe_discontinuity = TRUE) {
  structure(list(condition = condition, affect = affect, affect_neg = affect_neg,
                 initialize = initialize, finalize = finalize, rootfind = rootfind,
                 save_positions = save_positions, dtrelax = dtrelax,
                 maybe_discontinuity = maybe_discontinuity),
            class = c("ContinuousCallback", "DECallback"))
}

#' @rdname ContinuousCallback
#' @export
VectorContinuousCallback <- function(condition, affect, affect_neg = affect, ...) {
  cb <- ContinuousCallback(condition, affect, affect_neg, ...)
  cb$vector <- TRUE
  cb
}

#' @rdname ContinuousCallback
#' @export
DiscreteCallback <- function(condition, affect, initialize = NULL, finalize = NULL,
                             save_positions = c(TRUE, TRUE)) {
  structure(list(condition = condition, affect = affect, initialize = initialize,
                 finalize = finalize, save_positions = save_positions, tstops = numeric()),
            class = c("DiscreteCallback", "DECallback"))
}

#' @rdname ContinuousCallback
#' @param filter_tstops `PresetTimeCallback()`: with `TRUE` only the times
#'   inside the time span become stops; with `FALSE` every time is added
#'   with [add_tstop()], which is an error for a time before the start.
#' @export
PresetTimeCallback <- function(times, affect, initialize = NULL, filter_tstops = TRUE,
                               save_positions = c(TRUE, TRUE)) {
  if (!is.numeric(times)) stop("tstops must either be a number or a Vector", call. = FALSE)
  times <- sort(.as_num(times))
  # PresetTimeFunction of DiffEqCallbacks: the initialization adds the stops
  # and applies the affect when the initial time is one of the times
  init <- function(integ) {
    if (!is.null(initialize)) initialize(integ)
    if (!filter_tstops) for (tt in times) add_tstop(integ, tt)
    if (.near_any(times, integ$t)) affect(integ)
  }
  cb <- DiscreteCallback(function(u, t, integ) .near_any(times, t),
                         affect, initialize = init, save_positions = save_positions)
  cb$tstops <- if (filter_tstops) times else numeric()
  cb
}

#' @rdname ContinuousCallback
#' @param phase Time offset of the first event, `PeriodicCallback()`.
#' @param initial_affect,final_affect `PeriodicCallback()`: also apply the
#'   affect at the initial time, or when the integration reaches the final
#'   time (which is never a periodic stop).
#' @export
PeriodicCallback <- function(affect, dt, phase = 0, initial_affect = FALSE, final_affect = FALSE,
                             initialize = NULL, save_positions = c(TRUE, TRUE)) {
  if (phase < 0) stop("phase offset must be non-negative", call. = FALSE)
  # the stops t0 + phase + k dt strictly before the final time: the next stop
  # is only added while it lies before tstops_max (add_next_tstop!)
  times_for <- function(tsp) {
    dir <- if (tsp[2L] < tsp[1L]) -1 else 1
    first <- tsp[1L] + dir * phase
    if ((tsp[2L] - first) * dir < 0) return(numeric())
    tt <- seq(first, tsp[2L], by = dir * abs(dt))
    tt[tt != tsp[1L] & !vapply(tt, .near_any, TRUE, times = tsp[2L])]
  }
  st <- new.env(parent = emptyenv())
  # the condition of DiffEqCallbacks: a periodic time before the end, or the
  # end itself with final_affect
  cond <- function(u, t, integ) {
    fin <- .near_any(st$tf, t)
    (!fin && .near_any(st$times, t)) || (final_affect && fin)
  }
  init <- function(integ) {
    if (!is.null(initialize)) initialize(integ)
    st$times <- sort(times_for(integ$prob$tspan))
    st$tf <- integ$prob$tspan[2L]
    if (initial_affect) affect(integ)
  }
  cb <- DiscreteCallback(cond, affect, save_positions = save_positions, initialize = init)
  cb$tstops_fn <- times_for
  cb$periodic <- list(dt = dt, phase = phase)
  cb
}

# Whether t is within rounding of one of the sorted times.
.near_any <- function(times, t) {
  if (!length(times)) return(FALSE)
  i <- findInterval(t, times)
  tol <- 1e-12 * max(1, abs(t))
  (i >= 1L && abs(times[i] - t) <= tol) || (i < length(times) && abs(times[i + 1L] - t) <= tol)
}

#' @rdname ContinuousCallback
#' @export
CallbackSet <- function(...) {
  cbs <- list(...)
  out <- list(continuous = list(), discrete = list())
  for (cb in cbs) {
    if (is.null(cb)) next
    s <- .as_callback_set(cb)
    out$continuous <- c(out$continuous, s$continuous)
    out$discrete <- c(out$discrete, s$discrete)
  }
  structure(out, class = "CallbackSet")
}

.as_callback_set <- function(cb) {
  if (is.null(cb)) return(structure(list(continuous = list(), discrete = list()), class = "CallbackSet"))
  if (inherits(cb, "CallbackSet")) return(cb)
  if (inherits(cb, "ContinuousCallback")) return(structure(list(continuous = list(cb), discrete = list()), class = "CallbackSet"))
  if (inherits(cb, "DiscreteCallback")) return(structure(list(continuous = list(), discrete = list(cb)), class = "CallbackSet"))
  if (is.list(cb)) return(do.call(CallbackSet, cb))
  stop("not a callback: ", class(cb)[1L], call. = FALSE)
}

#' Stop an integration from a callback
#'
#' @param integ The integrator passed to an affect function.
#' @return Invisibly `NULL`.
#' @export
terminate <- function(integ) {
  integ$terminate <- TRUE
  invisible(NULL)
}

#' Reset the step size after an affect
#'
#' Calling `u_modified(integ)` tells the integrator that an affect changed
#' the state, so it restarts with a fresh step. Changes to `integ$u` made by
#' affects are detected automatically; this is for affects that change only
#' the parameters.
#'
#' @param integ The integrator.
#' @param value Logical.
#' @return Invisibly `NULL`.
#' @export
u_modified <- function(integ, value = TRUE) {
  integ$u_modified <- value
  integ$u_modified_flag <- isTRUE(value)
  invisible(NULL)
}

.init_callbacks <- function(cbs, integ) {
  for (cb in c(cbs$continuous, cbs$discrete)) {
    if (!is.null(cb$initialize)) {
      u_before <- integ$u; p_before <- integ$p
      cb$initialize(integ)
      if (!identical(u_before, integ$u) || !identical(p_before, integ$p)) integ$u_modified <- TRUE
    }
  }
  # continuous conditions start from their sign at t0
  integ$cb_last <- lapply(cbs$continuous, function(cb) cb$condition(integ$u, integ$t, integ))
}

.step_value <- function(step, t, dt = step$dt) {
  # dt: integrator.dt for the integrator interpolation; sol(t) passes the
  # difference of the saved times, as ode_interpolation of SciMLBase does
  theta <- if (dt == 0) 0 else (t - step$t0) / dt
  if (abs(theta) < 1e-15) return(step$y0)
  if (abs(theta - 1) < 1e-15) return(step$y1)
  step$interp(theta, dt, step$y0, step$y1, step$k)
}

# Find the first crossing of g on [ta, tb] with g(ta) = ga, g(tb) = gb of
# opposite signs. Returns the right end of the final bracket, so the
# condition has changed sign there, or the left end with `left = TRUE`
# (LeftRootFind), where it still has its old sign.
.find_root <- function(g, ta, tb, ga, gb, tol, left = FALSE) {
  # Illinois variant of regula falsi, with a bisection step whenever the
  # bracket shrinks too slowly, so it converges on steep functions too.
  side <- 0L
  for (i in 1:200) {
    w <- abs(tb - ta)
    if (w <= tol) return(if (left) ta else tb)
    tm <- (ta * gb - tb * ga) / (gb - ga)
    if (!is.finite(tm) || tm <= min(ta, tb) || tm >= max(ta, tb) || i %% 4L == 0L) tm <- (ta + tb) / 2
    gm <- g(tm)
    if (!is.finite(gm)) {
      tm <- (ta + tb) / 2
      gm <- g(tm)
    }
    if (gm == 0) return(tm)
    if (sign(gm) == sign(ga)) {
      ta <- tm; ga <- gm
      if (side == -1L) gb <- gb / 2
      side <- -1L
    } else {
      tb <- tm; gb <- gm
      if (side == 1L) ga <- ga / 2
      side <- 1L
    }
  }
  warning("event root finding did not converge; the event time may be inaccurate", call. = FALSE)
  tb
}

# SciMLBase's RootfindOpt: a Bool converts as `b ? LeftRootFind : NoRootFind`.
.rootfind_mode <- function(rf) {
  if (isTRUE(rf)) return("left")
  if (isFALSE(rf) || is.null(rf)) return("none")
  match.arg(as.character(rf), c("left", "right", "none"))
}

# Run the callbacks after an accepted step. Continuous events truncate the
# step at the first event time. Returns list(modified, save).
.apply_callbacks <- function(cbs, integ, step, ctx) {
  modified <- FALSE
  best <- NULL
  left <- NULL
  fired <- list()
  if (length(cbs$continuous)) {
    t0 <- step$t0; t1 <- integ$t
    best <- NULL
    tol <- 1e-12 * max(1, abs(t1))
    for (i in seq_along(cbs$continuous)) {
      cb <- cbs$continuous[[i]]
      gfun <- function(tt, j = NULL) {
        v <- cb$condition(.step_value(step, tt), tt, integ)
        if (is.null(j)) v else v[j]
      }
      g0 <- integ$cb_last[[i]]
      if (is.null(g0) || length(g0) == 0L) g0 <- cb$condition(step$y0, t0, integ)
      g1 <- cb$condition(integ$u, t1, integ)
      ncomp <- length(g1)
      for (j in seq_len(ncomp)) {
        ga0 <- g0[min(j, length(g0))]
        tb0 <- t0
        # after an event of this condition that left it at its event value,
        # start the search a little later (nudge_tprev of DiffEqBase
        # callbacks.jl: repeat_nudge = 1/100 of the step, abstol = 10 eps),
        # so a left root is not found again
        last <- integ$cb_last_event
        if (!is.null(last) && last$i == i && last$j == j && isTRUE(last$t == t0) &&
            abs(ga0 - last$err) <= (cb$abstol %||% (10 * .Machine$double.eps))) {
          tb0 <- t0 + (t1 - t0) * (cb$repeat_nudge %||% 0.01)
          ga0 <- gfun(tb0, j)
        }
        # sample inside the step so that pairs of roots are not missed
        tt <- seq(tb0, t1, length.out = 11L)
        gv <- c(ga0, vapply(tt[2:10], function(z) gfun(z, j), 0), g1[j])
        # a zero at the start is the previous event: use the sign just after it
        if (isTRUE(gv[1L] == 0)) gv[1L] <- gv[2L]
        for (s in 1:10) {
          ga <- gv[s]; gb <- gv[s + 1L]
          if (!is.finite(ga) || !is.finite(gb)) next
          up <- ga < 0 && gb >= 0
          down <- ga > 0 && gb <= 0
          if (!up && !down) next
          fn <- if (up) cb$affect else cb$affect_neg
          if (is.null(fn)) next
          rf <- .rootfind_mode(cb$rootfind)
          tr <- if (rf == "none") t1
                else if (gb != 0) .find_root(function(z) gfun(z, j), tt[s], tt[s + 1L], ga, gb, tol,
                                             left = rf == "left")
                else tt[s + 1L]
          if (is.null(best) || (tr - best$t) * sign(t1 - t0) < 0)
            best <- list(t = tr, cb = cb, fn = fn, idx = j, i = i, err = gfun(tr, j))
          break
        }
      }
    }
    if (!is.null(best)) {
      uev <- .step_value(step, best$t)
      integ$t <- best$t
      integ$u[] <- uev
      st <- integ$last_step
      st$tend <- best$t
      # the stages stay those of the full step (OrdinaryDiffEq recomputes them
      # for the shortened step), so its interpolant keeps the full dt
      st$truncated <- TRUE
      integ$last_step <- st
      left <- list(t = best$t, u = integ$u)
      u_before <- integ$u; p_before <- integ$p
      if (isTRUE(best$cb$vector)) best$fn(integ, best$idx) else best$fn(integ)
      fired <- c(fired, list(best$cb))
      integ$sol_events[[length(integ$sol_events) + 1L]] <- list(t = best$t, callback = best$i, index = best$idx)
      modified <- TRUE
    }
    integ$cb_last_event <- if (is.null(best)) NULL else list(i = best$i, j = best$idx, t = best$t, err = best$err)
    integ$cb_last <- lapply(cbs$continuous, function(cb) cb$condition(integ$u, integ$t, integ))
  }
  for (cb in cbs$discrete) {
    if (isTRUE(cb$condition(integ$u, integ$t, integ))) {
      if (is.null(left)) left <- list(t = integ$t, u = integ$u)
      integ$u_modified_flag <- NULL
      cb$affect(integ)
      fired <- c(fired, list(cb))
      # an affect that declares the state unchanged (u_modified(integ, FALSE),
      # derivative_discontinuity(integ, FALSE)) keeps the caches, as in SciMLBase
      if (!isFALSE(integ$u_modified_flag)) modified <- TRUE
    }
  }
  if (isTRUE(integ$u_modified)) modified <- TRUE
  if (modified && length(cbs$continuous))
    integ$cb_last <- lapply(cbs$continuous, function(cb) cb$condition(integ$u, integ$t, integ))
  # save_positions of the callbacks that fired; `force_save` callbacks (the
  # jump callbacks) save around the event also without save_everystep, as
  # DiffEqBase's apply_callback! does for every callback
  sp <- vapply(fired, function(cb) as.logical(cb$save_positions %||% c(TRUE, TRUE)), logical(2L))
  # a continuous event sets the proposed step to dtrelax times the step
  list(modified = modified, left = left, save_pre = any(sp[1L, ]), save_post = !length(fired) || any(sp[2L, ]),
       dtrelax = if (!is.null(best)) best$cb$dtrelax %||% 1,
       maybe_disco = !is.null(best) && !isFALSE(best$cb$maybe_discontinuity),
       force = any(vapply(fired, function(cb) isTRUE(cb$force_save), TRUE)))
}
