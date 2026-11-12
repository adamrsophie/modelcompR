# The callback library of DiffEqCallbacks.jl: saving.jl (SavedValues,
# SavingCallback, LinearizingSavingCallback with is_linear_enough! and
# linearize_period), independentlylinearizedutils.jl
# (IndependentlyLinearizedSolution), function_caller.jl
# (FunctionCallingCallback), domain.jl (PositiveDomain, GeneralDomain),
# manifold.jl (ManifoldProjection), autoabstol.jl (AutoAbstol),
# stepsizelimiters.jl (StepsizeLimiter), terminatesteadystate.jl
# (TerminateSteadyState), probints.jl (ProbIntsUncertainty,
# AdaptiveProbIntsUncertainty), iterative_and_periodic.jl
# (IterativeCallback), integrating.jl, integrating_sum.jl,
# integrating_GK_affect.jl and integrating_GK_sum.jl (the integrating
# callbacks with IntegrandValues and IntegrandValuesSum).
#
# All are discrete callbacks that act after every accepted step, through the
# integrator interface (R/integrator-interface.R): integ$interp(t) is the
# interpolant of the last step, from integ$tprev to integ$t.

.true_condition <- function(u, t, integ) TRUE

#' Callback library
#'
#' The callbacks of DiffEqCallbacks.jl, for [solve()] and [init()] of ODE,
#' DAE and SDE problems (pass them as `callback`, alone or in a
#' [CallbackSet()]).
#'
#' * `SavingCallback(save_func, saved_values)` calls `save_func(u, t, integ)`
#'   after every step (or at `saveat`) and stores the results in
#'   `saved_values`, a `SavedValues()` object with elements `t` and
#'   `saveval`. `FunctionCallingCallback(func)` calls `func(u, t, integ)` the
#'   same way without storing anything.
#' * `LinearizingSavingCallback(ils)` stores, in an
#'   `IndependentlyLinearizedSolution(prob)`, for each unknown separately the
#'   fewest time points such that linear interpolation is within `abstol` and
#'   `reltol` of the solver's interpolant. `ils(t)` interpolates, `ils$ts` and
#'   `ils$us` hold the points of each unknown.
#' * `PositiveDomain()` keeps the solution non negative: the next step is
#'   shortened until the extrapolated state is above `-abstol`, and negative
#'   values are set to zero. `GeneralDomain(g)` keeps `g(u, p, t)` (or `g(u,
#'   p)`) below `abstol` the same way and projects onto `g = 0`.
#' * `ManifoldProjection(g)` projects the state onto the manifold `g(u, p, t)
#'   = 0` (or `g(u, p)`) after every step, the closest point in the
#'   Euclidean norm (Newton on the Lagrange conditions).
#' * `AutoAbstol()` sets the absolute tolerance to the relative tolerance
#'   times the largest absolute value of each unknown so far.
#' * `StepsizeLimiter(dtFE)` limits each step to `safety_factor *
#'   dtFE(u, p, t)`.
#' * `TerminateSteadyState(abstol, reltol)` stops when every derivative is
#'   below `max(abstol, reltol |u|)` (after `min_t`), or when `test(integ,
#'   abstol, reltol, min_t)` is `TRUE`.
#' * `ProbIntsUncertainty(sigma, order)` adds Gaussian noise of standard
#'   deviation `sigma dt^order` after each step (the ProbInts method of
#'   Conrad et al.); `AdaptiveProbIntsUncertainty(order)` uses the error
#'   estimate of the step for `sigma`.
#' * `IterativeCallback(time_choice, affect)` applies `affect(integ)` at the
#'   times `time_choice(integ)` returns one after the other (`NULL` ends it).
#' * `IntegratingCallback(integrand_func, integrand_values)` integrates
#'   `integrand_func(u, t, integ)` over each step by Gauss-Legendre
#'   quadrature on the interpolant and stores one integral per step in
#'   `IntegrandValues()`; `IntegratingSumCallback()` accumulates the total
#'   in `IntegrandValuesSum()`. The `GK` versions use adaptive Gauss-Kronrod
#'   quadrature with `abstol` and `reltol`.
#'
#' @param save_func,func Functions `(u, t, integ)`.
#' @param saved_values A `SavedValues()` object.
#' @param saveat,funcat Times (or, as one number, an interval).
#' @param save_everystep,func_everystep Call after every step.
#' @param save_start,save_end,func_start Call at the start or the end.
#' @param tdir Direction of time.
#' @param ils An `IndependentlyLinearizedSolution()`.
#' @param prob An ODE problem.
#' @param num_derivatives Number of derivatives to store with each point (0 or
#'   1).
#' @param interpolate_mask Logical vector: the unknowns whose interpolant is
#'   refined (others are stored at the step ends).
#' @param abstol,reltol Tolerances (default: those of the integrator).
#' @param u Unused cache argument of the Julia signature.
#' @param save Save after the callback acts.
#' @param scalefactor Factor by which a step is shortened (default `1/2`).
#' @param g A residual function `(u, p, t)` or `(u, p)`.
#' @param autonomous Whether `g` takes no time argument (default: from the
#'   number of arguments).
#' @param manifold_jacobian,domain_jacobian Optional Jacobian of `g`.
#' @param nlsolve,... Unused options of the Julia signature.
#' @param init_curmax Initial largest absolute values for `AutoAbstol()`.
#' @param dtFE Function `(u, p, t)` giving the largest stable step.
#' @param safety_factor Fraction of `dtFE` allowed.
#' @param max_step Always propose the largest allowed step.
#' @param cached_dtcache Unused option of the Julia signature.
#' @param test Function `(integ, abstol, reltol, min_t)`.
#' @param min_t Earliest time to stop.
#' @param wrap_test Pass `test` as the condition `(u, t, integ)` itself when
#'   `FALSE`.
#' @param sigma,order Noise scale and order of the method.
#' @param time_choice Function `(integ)` returning the next time or `NULL`.
#' @param affect Function `(integ)`.
#' @param initial_affect Apply `affect` at the start.
#' @param integrand_func Function `(u, t, integ)`.
#' @param integrand_values An `IntegrandValues()` (or `IntegrandValuesSum()`)
#'   object.
#' @param integrand_prototype Unused (the integrand's shape is taken from its
#'   first value).
#' @return A callback, or a storage object.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, c(1, 2), c(0, 1))
#' sv <- SavedValues()
#' sol <- solve(prob, Tsit5(), callback = SavingCallback(function(u, t, integ) sum(u), sv,
#'                                                       saveat = 0.25))
#' unlist(sv$saveval)
#' iv <- IntegrandValuesSum()
#' sol <- solve(prob, Tsit5(), callback = IntegratingSumCallback(function(u, t, integ) u, iv))
#' iv$integrand # 1 - exp(-1) times c(1, 2)
#' @param save_positions Save the state before and after the affect.
#' @name callback_library
NULL

# ---------------------------------------------------------------------------
# Saving

#' @rdname callback_library
#' @export
SavedValues <- function() {
  e <- new.env(parent = emptyenv())
  e$t <- numeric()
  e$saveval <- list()
  class(e) <- "SavedValues"
  e
}

#' @export
print.SavedValues <- function(x, ...) {
  cat("SavedValues with", length(x$t), "saved values\n")
  invisible(x)
}

# the save or call times of a SavingCallback or FunctionCallingCallback
.cb_times <- function(at, tspan, start, end) {
  if (length(at) == 1L) {
    t0 <- tspan[1L] + if (!start) at else 0
    tf <- tspan[2L] - if (!end) at else 0
    return(seq(t0, tf, by = sign(tspan[2L] - tspan[1L]) * abs(at)))
  }
  sort(as.numeric(at), decreasing = tspan[2L] < tspan[1L])
}

#' @rdname callback_library
#' @export
SavingCallback <- function(save_func, saved_values, saveat = numeric(), save_everystep = !length(saveat),
                           save_start = save_everystep || !length(saveat) || length(saveat) == 1L,
                           save_end = save_everystep || !length(saveat) || length(saveat) == 1L, tdir = 1) {
  if (!inherits(saved_values, "SavedValues")) stop("`saved_values` must be SavedValues()", call. = FALSE)
  st <- new.env(parent = emptyenv())
  push <- function(t, val) {
    saved_values$t <- c(saved_values$t, t)
    saved_values$saveval[[length(saved_values$saveval) + 1L]] <- val
  }
  affect <- function(integ) {
    I <- .integ_I(integ)
    td <- I$tdir %||% tdir
    just_saved <- FALSE
    while (length(st$queue) && td * st$queue[1L] <= td * integ$t) {
      curt <- st$queue[1L]
      st$queue <- st$queue[-1L]
      if (curt != integ$t) {
        push(curt, save_func(integ$interp(curt), curt, integ))
      } else {
        just_saved <- TRUE
        push(integ$t, save_func(integ$u, integ$t, integ))
      }
    }
    if ((!just_saved && save_everystep) || (st$save_end && integ$t == integ$prob$tspan[2L]))
      push(integ$t, save_func(integ$u, integ$t, integ))
    derivative_discontinuity(integ, FALSE)
  }
  init <- function(integ) {
    tspan <- integ$prob$tspan
    st$save_end <- save_end
    st$queue <- if (length(saveat)) .cb_times(saveat, tspan, save_start, save_end) else numeric()
    if (length(saveat) == 1L && length(st$queue) && st$queue[length(st$queue)] == tspan[2L]) st$save_end <- FALSE
    saved_values$t <- numeric()
    saved_values$saveval <- list()
    if (save_start) {
      if (length(st$queue) && st$queue[1L] == integ$t) st$queue <- st$queue[-1L]
      push(integ$t, save_func(integ$u, integ$t, integ))
    }
  }
  DiscreteCallback(.true_condition, affect, initialize = init, save_positions = c(FALSE, FALSE))
}

#' @rdname callback_library
#' @export
FunctionCallingCallback <- function(func, funcat = numeric(), func_everystep = !length(funcat),
                                    func_start = TRUE, tdir = 1) {
  st <- new.env(parent = emptyenv())
  affect <- function(integ) {
    I <- .integ_I(integ)
    td <- I$tdir %||% tdir
    while (length(st$queue) && td * st$queue[1L] <= td * integ$t) {
      curt <- st$queue[1L]
      st$queue <- st$queue[-1L]
      if (curt != integ$t) func(integ$interp(curt), curt, integ) else func(integ$u, integ$t, integ)
    }
    if (func_everystep) func(integ$u, integ$t, integ)
    derivative_discontinuity(integ, FALSE)
  }
  init <- function(integ) {
    st$queue <- if (length(funcat)) .cb_times(funcat, integ$prob$tspan, func_start, TRUE) else numeric()
    if (func_start) {
      if (length(st$queue) && st$queue[1L] == integ$t) st$queue <- st$queue[-1L]
      func(integ$u, integ$t, integ)
    }
  }
  DiscreteCallback(.true_condition, affect, initialize = init, save_positions = c(FALSE, FALSE))
}

#' @rdname callback_library
#' @export
IndependentlyLinearizedSolution <- function(prob, num_derivatives = 0) {
  if (!num_derivatives %in% 0:1) stop("num_derivatives must be 0 or 1", call. = FALSE)
  e <- new.env(parent = emptyenv())
  n <- length(prob$u0)
  e$ts <- rep(list(numeric()), n)
  e$us <- rep(list(if (num_derivatives) matrix(numeric(), 0L, 2L) else numeric()), n)
  e$num_derivatives <- num_derivatives
  e$names <- names(prob$u0)
  f <- function(t, idxs = seq_len(n)) {
    vapply(idxs, function(i) {
      ts <- e$ts[[i]]
      us <- if (e$num_derivatives) e$us[[i]][, 1L] else e$us[[i]]
      if (!length(ts)) return(NA_real_)
      stats::approx(ts, us, xout = t, rule = 2, ties = "ordered")$y
    }, 0)
  }
  attr(f, "data") <- e
  class(f) <- c("IndependentlyLinearizedSolution", "function")
  f
}

#' @export
`$.IndependentlyLinearizedSolution` <- function(x, name) get(name, envir = attr(x, "data"))

#' @export
print.IndependentlyLinearizedSolution <- function(x, ...) {
  e <- attr(x, "data")
  cat("IndependentlyLinearizedSolution with", length(e$ts), "unknowns and",
      sum(lengths(e$ts)), "points\n")
  invisible(x)
}

# is_linear_enough!: the three quartiles of [t0, t1] against the line
.lin_check <- function(t0, t1, u0, u1, interp, abstol, reltol) {
  tq <- t0 + (t1 - t0) * (1:3) / 4
  yl <- sapply(tq, function(t) u0 + (t - t0) * (u1 - u0) / (t1 - t0))
  yi <- sapply(tq, interp)
  yl <- matrix(yl, length(u0)); yi <- matrix(yi, length(u0))
  err <- abs(yl - yi)
  tol <- abstol + reltol * pmax(abs(yl), abs(yi))
  ok <- rowSums(err <= tol) == 3L
  worst <- which(err == max(err), arr.ind = TRUE)[1L, 2L]
  list(ok = ok, tb = tq[worst])
}

#' @rdname callback_library
#' @export
LinearizingSavingCallback <- function(ils, interpolate_mask = NULL, abstol = NULL, reltol = NULL) {
  e <- attr(ils, "data")
  n <- length(e$ts)
  mask <- interpolate_mask %||% rep(TRUE, n)
  store <- function(integ, t, u, sel) {
    du <- if (e$num_derivatives) as.numeric(.integ_I(integ)$ctx$f(u, integ$p, t))
    for (i in which(sel)) {
      e$ts[[i]] <- c(e$ts[[i]], t)
      e$us[[i]] <- if (e$num_derivatives) rbind(e$us[[i]], c(u[i], du[i])) else c(e$us[[i]], u[i])
    }
  }
  linearize <- function(integ, t0, t1, u0, u1, sel, dtmin, at, rt) {
    if (t1 - t0 < dtmin) stop("Linearization failed, fell below linearization subdivision threshold", call. = FALSE)
    chk <- .lin_check(t0, t1, u0, u1, integ$interp, at, rt)
    nonlin <- !chk$ok & sel & mask
    if (any(nonlin)) {
      ub <- integ$interp(chk$tb)
      linearize(integ, t0, chk$tb, u0, ub, nonlin, dtmin, at, rt)
      linearize(integ, chk$tb, t1, ub, u1, sel, dtmin, at, rt)
    } else store(integ, t1, u1, sel)
  }
  affect <- function(integ) {
    t0 <- integ$tprev; t1 <- integ$t
    u0 <- integ$interp(t0); u1 <- integ$u
    derivative_discontinuity(integ, FALSE)
    if (!length(e$ts[[1L]])) store(integ, t0, u0, rep(TRUE, n))
    if (t1 == t0) return(invisible())
    I <- .integ_I(integ)
    linearize(integ, t0, t1, u0, u1, rep(TRUE, n), .Machine$double.eps * abs(t1 - t0) * 1000,
              abstol %||% I$abstol, reltol %||% I$reltol)
  }
  init <- function(integ) {
    e$ts <- rep(list(numeric()), n)
    e$us <- rep(list(if (e$num_derivatives) matrix(numeric(), 0L, 2L) else numeric()), n)
  }
  DiscreteCallback(.true_condition, affect, initialize = init, save_positions = c(FALSE, FALSE))
}

# ---------------------------------------------------------------------------
# Domains and manifolds

# the affect of the domain callbacks (affect! of domain.jl): shorten the
# proposed step until the extrapolated state is accepted
.domain_affect <- function(accepted, modify, abstol, scalefactor) function(integ) {
  I <- .integ_I(integ)
  if (!isTRUE(I$adaptive)) stop("domain callback can only be applied to adaptive algorithms", call. = FALSE)
  tol <- abstol %||% I$abstol
  sf <- scalefactor %||% (1 / 2)
  dt <- get_proposed_dt(integ)
  ft <- .integ_first_tstop(integ)
  if (!is.null(ft)) dt <- min(abs(dt), abs(ft - integ$t))
  td <- I$tdir
  t <- integ$t + td * dt
  modified <- FALSE
  while (dt > 0) {
    if (accepted(integ$interp(t), integ$p, t, tol)) break
    dtc <- dt
    dt <- dt * sf
    modified <- TRUE
    t <- integ$t + td * dt
    if (dtc == dt) {
      warning("Could not restrict values to domain. Iteration was canceled since proposed time step dt = ",
              dt, " could not be reduced.", call. = FALSE)
      break
    }
  }
  if (modified) set_proposed_dt(integ, 9 / 10 * dt)
  derivative_discontinuity(integ, modify(integ))
}

#' @rdname callback_library
#' @export
PositiveDomain <- function(u = NULL, save = TRUE, abstol = NULL, scalefactor = NULL) {
  accepted <- function(uu, p, t, tol) {
    if (length(tol) > 1L && length(tol) != length(uu))
      stop("numbers of states and tolerances do not match", call. = FALSE)
    all(uu > -tol)
  }
  modify <- function(integ) {
    neg <- integ$u < 0
    if (any(neg)) integ$u[neg] <- 0
    any(neg)
  }
  DiscreteCallback(.true_condition, .domain_affect(accepted, modify, abstol, scalefactor),
                   save_positions = c(FALSE, save))
}

.wrap_manifold <- function(g, autonomous) {
  if (is.null(g)) return(NULL)
  auto <- autonomous %||% (length(formals(g)) == 2L)
  if (auto) function(u, p, t) g(u, p) else g
}

#' @rdname callback_library
#' @export
GeneralDomain <- function(g, u = NULL, save = TRUE, abstol = NULL, scalefactor = NULL, autonomous = NULL,
                          domain_jacobian = NULL, ...) {
  gg <- .wrap_manifold(g, autonomous)
  accepted <- function(uu, p, t, tol) {
    r <- as.numeric(gg(uu, p, t))
    if (length(tol) > 1L && length(tol) != length(r))
      stop("numbers of residuals and tolerances do not match", call. = FALSE)
    all(r < tol)
  }
  domain <- DiscreteCallback(.true_condition, .domain_affect(accepted, function(integ) FALSE, abstol, scalefactor),
                             save_positions = c(FALSE, save))
  CallbackSet(ManifoldProjection(g, save = FALSE, autonomous = autonomous, manifold_jacobian = domain_jacobian),
              domain)
}

#' @rdname callback_library
#' @export
ManifoldProjection <- function(g, nlsolve = NULL, save = TRUE, autonomous = NULL, manifold_jacobian = NULL, ...) {
  gg <- .wrap_manifold(g, autonomous)
  jac <- .wrap_manifold(manifold_jacobian, autonomous)
  jacobian <- function(u, p, t) {
    if (!is.null(jac)) return(as.matrix(jac(u, p, t)))
    r0 <- as.numeric(gg(u, p, t))
    J <- matrix(0, length(r0), length(u))
    for (i in seq_along(u)) {
      h <- sqrt(.Machine$double.eps) * max(abs(u[i]), 1)
      up <- u
      up[i] <- up[i] + h
      J[, i] <- (as.numeric(gg(up, p, t)) - r0) / h
    }
    J
  }
  # the closest point of g = 0: z = (lambda, u) with u - u_hat + J0' lambda = 0
  # and g(u) = 0, J0 the Jacobian at u_hat (NonlinearSolveManifoldProjectionCache)
  affect <- function(integ) {
    uh <- as.numeric(integ$u); p <- integ$p; t <- integ$t
    J0 <- jacobian(uh, p, t)
    m <- nrow(J0); n <- length(uh)
    u <- uh; lam <- numeric(m)
    ok <- FALSE
    for (it in seq_len(50L)) {
      F <- c(u - uh + as.numeric(crossprod(J0, lam)), as.numeric(gg(u, p, t)))
      if (max(abs(F)) <= 10 * .Machine$double.eps * max(1, max(abs(u)))) {
        ok <- TRUE
        break
      }
      Ju <- jacobian(u, p, t)
      K <- rbind(cbind(diag(1, n), t(J0)), cbind(Ju, matrix(0, m, m)))
      d <- tryCatch(solve(K, -F), error = function(e) NULL)
      if (is.null(d)) break
      u <- u + d[seq_len(n)]
      lam <- lam + d[n + seq_len(m)]
      if (max(abs(d)) <= 1e-12 * max(1, max(abs(u)))) {
        ok <- TRUE
        break
      }
    }
    if (!ok) {
      terminate(integ)
      return(invisible())
    }
    integ$u[] <- u
  }
  DiscreteCallback(.true_condition, affect, save_positions = c(FALSE, save))
}

# ---------------------------------------------------------------------------
# Tolerances, step sizes, termination, uncertainty

#' @rdname callback_library
#' @export
AutoAbstol <- function(save = TRUE, init_curmax = 0) {
  st <- new.env(parent = emptyenv())
  st$curmax <- abs(init_curmax)
  affect <- function(integ) {
    I <- .integ_I(integ)
    st$curmax <- if (length(st$curmax) > 1L) pmax(st$curmax, abs(integ$u)) else max(st$curmax, abs(integ$u))
    set_abstol(integ, st$curmax * I$reltol)
    derivative_discontinuity(integ, FALSE)
  }
  init <- function(integ) {
    if (all(st$curmax == 0)) st$curmax <- .integ_I(integ)$abstol
  }
  DiscreteCallback(.true_condition, affect, initialize = init, save_positions = c(save, FALSE))
}

#' @rdname callback_library
#' @export
StepsizeLimiter <- function(dtFE, safety_factor = 9 / 10, max_step = FALSE, cached_dtcache = 0) {
  affect <- function(integ) {
    I <- .integ_I(integ)
    I$dtmax <- safety_factor * dtFE(integ$u, integ$p, integ$t)
    if (max_step) set_proposed_dt(integ, I$dtmax)
    derivative_discontinuity(integ, FALSE)
  }
  DiscreteCallback(.true_condition, affect, initialize = affect, save_positions = c(FALSE, FALSE))
}

# allDerivPass of terminatesteadystate.jl
.all_deriv_pass <- function(integ, abstol, reltol, min_t) {
  if (!is.null(min_t) && integ$t < min_t) return(FALSE)
  d <- if (inherits(integ$prob, "DiscreteProblem")) integ$u - integ$uprev else get_du(integ)
  all(abs(d) <= pmax(abstol, reltol * abs(integ$u)))
}

#' @rdname callback_library
#' @export
TerminateSteadyState <- function(abstol = 1e-8, reltol = 1e-6, test = .all_deriv_pass, min_t = NULL,
                                 wrap_test = TRUE) {
  condition <- if (wrap_test) function(u, t, integ) test(integ, abstol, reltol, min_t) else test
  DiscreteCallback(condition, function(integ) terminate(integ), save_positions = c(TRUE, FALSE))
}

#' @rdname callback_library
#' @export
ProbIntsUncertainty <- function(sigma, order, save = TRUE) {
  DiscreteCallback(.true_condition, function(integ) {
    integ$u[] <- integ$u + sigma * sqrt(integ$dt^(2 * order)) * stats::rnorm(length(integ$u))
  }, save_positions = c(save, FALSE))
}

#' @rdname callback_library
#' @export
AdaptiveProbIntsUncertainty <- function(order, save = TRUE) {
  DiscreteCallback(.true_condition, function(integ) {
    integ$u[] <- integ$u + integ$EEst * sqrt(integ$dt^(2 * order)) * stats::rnorm(length(integ$u))
  }, save_positions = c(save, FALSE))
}

#' @rdname callback_library
#' @export
IterativeCallback <- function(time_choice, affect, initial_affect = FALSE, save_positions = c(TRUE, TRUE)) {
  st <- new.env(parent = emptyenv())
  st$tnext <- Inf
  aff <- function(integ) {
    affect(integ)
    tnew <- time_choice(integ)
    st$tnext <- tnew %||% NA_real_
    if (!is.null(tnew)) add_tstop(integ, tnew)
  }
  init <- function(integ) {
    if (initial_affect) {
      st$tnext <- integ$t
      aff(integ)
    } else {
      tn <- time_choice(integ)
      st$tnext <- tn %||% NA_real_
      if (!is.null(tn)) add_tstop(integ, tn)
    }
  }
  DiscreteCallback(function(u, t, integ) isTRUE(abs(t - st$tnext) <= 1e-12 * max(1, abs(t))), aff,
                   initialize = init, save_positions = save_positions)
}

# ---------------------------------------------------------------------------
# Integrating callbacks

# Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch)
.gauss_legendre <- function(n) {
  if (n == 1L) return(list(x = 0, w = 2))
  b <- seq_len(n - 1L) / sqrt(4 * seq_len(n - 1L)^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1L), 2:n)] <- b
  J[cbind(2:n, seq_len(n - 1L))] <- b
  ev <- eigen(J, symmetric = TRUE)
  o <- order(ev$values)
  list(x = ev$values[o], w = 2 * ev$vectors[1L, o]^2)
}

# the integral of f(u(t), t, integ) over the last step
.step_integral <- function(integ, f, gk = FALSE, abstol = 1e-6, reltol = 1e-3) {
  a <- integ$tprev; b <- integ$t
  if (a == b) return(0 * f(integ$u, b, integ))
  g <- function(s) f(integ$interp(s), s, integ)
  if (gk) {
    v <- g(b)
    out <- vapply(seq_along(v), function(i) .sens_quadgk(function(s) g(s)[i], a, b, abstol, reltol), 0)
    return(if (length(v) == 1L) out[1L] else out)
  }
  alg <- .integ_I(integ)$alg
  n <- if (inherits(integ$prob, c("SDEProblem", "RODEProblem"))) 10L else max(1L, (as.integer(alg$order %||% 1) + 1L) %/% 2L)
  gl <- .gauss_legendre(min(n, 10L))
  acc <- 0
  for (k in seq_along(gl$x)) acc <- acc + gl$w[k] * g((b - a) / 2 * gl$x[k] + (a + b) / 2)
  acc * (b - a) / 2
}

#' @rdname callback_library
#' @export
IntegrandValues <- function() {
  e <- new.env(parent = emptyenv())
  e$ts <- numeric()
  e$integrand <- list()
  class(e) <- "IntegrandValues"
  e
}

#' @rdname callback_library
#' @export
IntegrandValuesSum <- function() {
  e <- new.env(parent = emptyenv())
  e$integrand <- 0
  class(e) <- "IntegrandValuesSum"
  e
}

#' @export
print.IntegrandValues <- function(x, ...) {
  cat("IntegrandValues with", length(x$ts), "steps\n")
  invisible(x)
}

#' @export
print.IntegrandValuesSum <- function(x, ...) {
  cat("IntegrandValuesSum:", format(x$integrand), "\n")
  invisible(x)
}

.integrating_cb <- function(integrand_func, integrand_values, sum, gk, abstol = 1e-6, reltol = 1e-3) {
  cls <- if (sum) "IntegrandValuesSum" else "IntegrandValues"
  if (!inherits(integrand_values, cls)) .stopf("`integrand_values` must be %s()", cls)
  affect <- function(integ) {
    v <- .step_integral(integ, integrand_func, gk, abstol, reltol)
    if (sum) integrand_values$integrand <- integrand_values$integrand + v else {
      integrand_values$ts <- c(integrand_values$ts, integ$t)
      integrand_values$integrand[[length(integrand_values$integrand) + 1L]] <- v
    }
    derivative_discontinuity(integ, FALSE)
  }
  init <- function(integ) {
    if (sum) integrand_values$integrand <- 0 else {
      integrand_values$ts <- numeric()
      integrand_values$integrand <- list()
    }
  }
  DiscreteCallback(.true_condition, affect, initialize = init, save_positions = c(FALSE, FALSE))
}

#' @rdname callback_library
#' @export
IntegratingCallback <- function(integrand_func, integrand_values, integrand_prototype = NULL)
  .integrating_cb(integrand_func, integrand_values, FALSE, FALSE)

#' @rdname callback_library
#' @export
IntegratingSumCallback <- function(integrand_func, integrand_values, integrand_prototype = NULL)
  .integrating_cb(integrand_func, integrand_values, TRUE, FALSE)

#' @rdname callback_library
#' @export
IntegratingGKCallback <- function(integrand_func, integrand_values, integrand_prototype = NULL,
                                  abstol = 1e-6, reltol = 1e-3)
  .integrating_cb(integrand_func, integrand_values, FALSE, TRUE, abstol, reltol)

#' @rdname callback_library
#' @export
IntegratingGKSumCallback <- function(integrand_func, integrand_values, integrand_prototype = NULL,
                                     abstol = 1e-6, reltol = 1e-3)
  .integrating_cb(integrand_func, integrand_values, TRUE, TRUE, abstol, reltol)
