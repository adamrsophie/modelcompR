# Forward sensitivity analysis, after SciMLSensitivity.jl v7.119
# src/forward_sensitivity.jl (ODEForwardSensitivityFunction,
# ODEForwardSensitivityProblem, extract_local_sensitivities, _extract) and
# src/callback_tracking.jl (the event corrections, here applied forward).
#
# The problem is the original ODE augmented with the sensitivities
# S = du/dtheta (n x m, stored column by column):
#   dS/dt = J S + f_p Q,     S(t0) = du0/dtheta,
# where Q = dp/dtheta. ForwardSensitivity takes theta = p as in Julia;
# ForwardDiffSensitivity (Julia's dual number method) takes the seed matrices
# du0 and dp, so theta can include the initial values. When the problem was
# built from a System whose initial values depend on the parameters through
# initialization, S(t0) includes du0/dp.
#
# Events: each affect is replaced by one that also maps (S, Q) through the
# linearization of the event (see .sens_event_jacobian), so sensitivities
# are correct across events that move with the parameters and that change
# states or parameters. Q is then part of the state (constant between
# events), since parameter changing events update it.

#' Forward sensitivity problems
#'
#' `ODEForwardSensitivityProblem(prob)` (or `(f, u0, tspan, p)`) builds an
#' ODE problem whose state is the original state followed by its
#' sensitivities, solved with [solve()] like any `ODEProblem`.
#' `extract_local_sensitivities(sol)` splits the solution into the values
#' and the sensitivities.
#'
#' With `sensealg = ForwardSensitivity()` the sensitivities are with respect
#' to the parameters `p`. `ForwardDiffSensitivity()` takes seed matrices:
#' `du0` (n x m) and `dp` (length(p) x m) give the derivatives of the
#' initial values and the parameters with respect to `m` directions
#' (`cbind(diag(n), 0)` and `cbind(0, diag(length(p)))` differentiate with
#' respect to both). For a problem built from a System, `du0` adds the
#' dependence of the initial values on the parameters through
#' initialization. Events (callbacks of the problem) are differentiated,
#' including the dependence of the event time on the parameters.
#'
#' Derivatives of the right-hand side come from the System (symbolic
#' Jacobians compiled to R) or, for R functions, from `jac_fn` or central
#' finite differences.
#'
#' @param f An `ODEProblem`, a complete `System`, an R function
#'   `f(u, p, t)` or an `ODEFunction`.
#' @param u0,tspan,p Initial values, time span and parameters (ignored when
#'   `f` is a problem; for a System `u0` is the operating point).
#' @param sensealg [ForwardSensitivity()] or [ForwardDiffSensitivity()].
#' @param du0,dp Seed matrices (see Details).
#' @param ... Further arguments for [ODEProblem()] when `f` is a System or
#'   function.
#' @return An `ODEProblem` with the augmented state.
#' @examples
#' prob <- ODEForwardSensitivityProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
#' sol <- solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-10)
#' s <- extract_local_sensitivities(sol, length(sol$t))
#' s$du[[1]] + exp(-2)   # d u(1) / dp = -exp(-2)
#' @export
ODEForwardSensitivityProblem <- function(f, u0 = NULL, tspan = NULL, p = NULL, sensealg = ForwardSensitivity(),
                                         du0 = NULL, dp = NULL, ...) {
  prob <- if (inherits(f, "ODEProblem")) f else if (inherits(f, "System"))
    ODEProblem(f, u0 %||% list(), tspan, ...) else ODEProblem(f, u0, tspan, p, ...)
  if (!.is_sensealg(sensealg, "ForwardSensitivityAlgorithm"))
    stop("sensealg must be ForwardSensitivity() or ForwardDiffSensitivity()", call. = FALSE)
  n <- length(prob$u0); np <- length(prob$p)
  if (is.null(du0) && is.null(dp) && !np)
    stop("You must have parameters to use parameter sensitivity calculations!", call. = FALSE)
  ini <- .sens_init_jacobian(prob)
  dp <- dp %||% diag(1, np)
  du0 <- du0 %||% matrix(0, n, ncol(dp))
  dp <- as.matrix(dp); du0 <- as.matrix(du0)
  if (ncol(du0) != ncol(dp)) stop("Same number of perturbations of initial conditions and parameters required",
                                  call. = FALSE)
  if (nrow(du0) != n) stop("Perturbations for all initial conditions required", call. = FALSE)
  if (nrow(dp) != np) stop("Perturbations for all parameters required", call. = FALSE)
  # theta moves the problem parameters, which move u0 and the effective
  # parameters through initialization
  S0 <- du0 + ini$du0 %*% dp
  Q0 <- ini$dp %*% dp
  .sens_forward_problem(prob, S0, Q0, sensealg, colnames(dp))
}

.sens_forward_problem <- function(prob, S0, Q0, sensealg, theta_names = NULL) {
  n <- length(prob$u0); np <- length(prob$p); m <- ncol(S0)
  fns <- .sens_fns(prob, sensealg)
  cbs <- .as_callback_set(prob$kwargs$callback)
  has_cb <- length(cbs$continuous) + length(cbs$discrete) > 0L
  if (has_cb && !is.null(fns$mass))
    stop("events in mass matrix problems are not supported by the sensitivity methods", call. = FALSE)
  if (has_cb && prob$tspan[2L] < prob$tspan[1L])
    stop("events with a decreasing time span are not supported by the sensitivity methods", call. = FALSE)
  nq <- if (has_cb) np * m else 0L
  iu <- seq_len(n); iS <- n + seq_len(n * m); iQ <- n + n * m + seq_len(nq)
  f <- fns$f
  rhs <- function(z, p, t) {
    u <- z[iu]
    S <- matrix(z[iS], n, m)
    Q <- if (has_cb) matrix(z[iQ], np, m) else Q0
    dS <- fns$jac(u, p, t) %*% S
    if (np) dS <- dS + fns$pjac(u, p, t) %*% Q
    c(f(u, p, t), as.numeric(dS), numeric(nq))
  }
  M <- fns$mass
  Maug <- if (is.null(M)) NULL else
    .sens_blockdiag(c(rep(list(M), m + 1L), if (nq) list(diag(1, nq))))
  un <- names(prob$u0) %||% paste0("u", iu)
  th <- theta_names %||% if (m == np && !is.null(names(prob$p))) names(prob$p) else paste0("theta", seq_len(m))
  z0 <- c(as.numeric(prob$u0), as.numeric(S0), if (nq) as.numeric(Q0))
  names(z0) <- c(un, paste0("d(", rep(un, m), ")/d(", rep(th, each = n), ")"),
                 if (nq) paste0("d(", rep(names(prob$p) %||% paste0("p", seq_len(np)), m), ")/d(",
                                rep(th, each = np), ")"))
  kw <- prob$kwargs %||% list()
  rec <- new.env(parent = emptyenv())
  rec$events <- list()
  if (has_cb) kw$callback <- .sens_forward_callbacks(cbs, prob, fns, n, np, m, rec)
  fn <- ODEFunction(rhs, mass_matrix = Maug)
  out <- .new_problem("ODEProblem", f = fn, u0 = z0, tspan = prob$tspan, p = prob$p, sys = prob$sys,
                      kwargs = kw)
  out$sensitivity <- list(n = n, m = m, np = np, sensealg = sensealg, theta = th, nq = nq,
                          original = prob, events = rec)
  out
}

.sens_blockdiag <- function(mats) {
  sz <- vapply(mats, nrow, 0L)
  out <- matrix(0, sum(sz), sum(sz))
  off <- 0L
  for (A in mats) {
    k <- nrow(A)
    out[off + seq_len(k), off + seq_len(k)] <- A
    off <- off + k
  }
  out
}

# Callbacks of the augmented problem: conditions see the original state,
# affects also map the sensitivities through the event linearization.
# Each event is recorded with the change of its time, dtau = tau_x' X-, and
# the derivatives f- and f+, for samples saved at the (moving) event time.
.sens_forward_callbacks <- function(cbs, prob, fns, n, np, m, rec) {
  iu <- seq_len(n)
  wrap_cond <- function(cond) {
    force(cond)
    function(u, t, integ) cond(u[iu], t, integ)
  }
  wrap_affect <- function(aff, cb, continuous) {
    if (is.null(aff)) return(NULL)
    force(aff)
    function(integ, ...) {
      z <- integ$u
      ev <- .sens_event_record(integ, aff, list(...), cb, continuous, n)
      fake <- .sens_fake_integ(ev$u, ev$p, ev$t, integ)
      do.call(aff, c(list(fake), ev$args))
      ev$u_right <- as.numeric(fake$u); ev$p_right <- fake$p
      D <- .sens_event_jacobian(ev, fns, prob)
      X <- rbind(matrix(z[n + seq_len(n * m)], n, m), matrix(z[n + n * m + seq_len(np * m)], np, m))
      tx <- attr(D, "tau_x")
      rec$events[[length(rec$events) + 1L]] <- list(
        t = ev$t, dtau = if (is.null(tx)) numeric(m) else as.numeric(crossprod(tx, X)),
        fm = fns$f(ev$u, ev$p, ev$t), fp = fns$f(ev$u_right, .sens_pvec(ev$p, ev$p_right), ev$t))
      X <- D %*% X
      integ$u <- c(ev$u_right, as.numeric(X[iu, , drop = FALSE]), as.numeric(X[n + seq_len(np), , drop = FALSE]))
      integ$p[] <- as.numeric(fake$p)
      if (isTRUE(fake$terminate)) integ$terminate <- TRUE
      invisible(NULL)
    }
  }
  conts <- lapply(cbs$continuous, function(cb) {
    cb2 <- cb
    cb2$condition <- wrap_cond(cb$condition)
    cb2$affect <- wrap_affect(cb$affect, cb, TRUE)
    cb2$affect_neg <- wrap_affect(cb$affect_neg, cb, TRUE)
    cb2
  })
  discs <- lapply(cbs$discrete, function(cb) {
    cb2 <- cb
    cb2$condition <- wrap_cond(cb$condition)
    cb2$affect <- wrap_affect(cb$affect, cb, FALSE)
    cb2
  })
  reset <- DiscreteCallback(function(u, t, integ) FALSE, function(integ) NULL,
                            initialize = function(integ) rec$events <- list())
  structure(list(continuous = conts, discrete = c(list(reset), discs)), class = "CallbackSet")
}

# ---------------------------------------------------------------------------
# Events. An event at time tau maps x- = (u-, p-) to x+ = R(x-, tau); for a
# continuous callback tau moves with the trajectory, c(u(tau), p, tau) = 0.
# With F = (f, 0), the linearization (the saltation matrix) is
#   D = R_x + b tau_x',  b = R_x F- + R_t - F+,
#   tau_x = -c_x / (c_t + c_u f-),
# forward X+ = D X-, adjoint z- = D' z+ (callback_tracking.jl, eq. 13 of
# arXiv:1905.10403, with the parameter changing and time dependent terms).
# Discrete and preset time callbacks do not move: tau_x = 0.

.sens_event_record <- function(integ, aff, args, cb, continuous, n) {
  list(t = integ$t, tprev = integ$tprev, u = as.numeric(integ$u)[seq_len(n)], p = integ$p,
       affect = aff, args = args, condition = cb$condition, continuous = continuous,
       idx = if (continuous && length(args)) args[[1L]] else 1L)
}

# An integrator stand-in for running affects and conditions on given values.
.sens_fake_integ <- function(u, p, t, integ = NULL) {
  e <- new.env(parent = emptyenv())
  e$u <- u; e$p <- p; e$t <- t
  e$tprev <- integ$tprev %||% t; e$uprev <- u
  e$prob <- integ$prob
  e$terminate <- FALSE; e$u_modified <- FALSE
  e$sol_events <- list()
  e
}

.sens_event_jacobian <- function(ev, fns, prob) {
  n <- length(ev$u); np <- length(ev$p)
  pn <- ev$p
  reset <- function(x, t) {
    q <- pn
    q[] <- x[n + seq_len(np)]
    e <- .sens_fake_integ(x[seq_len(n)], q, t)
    e$prob <- prob
    do.call(ev$affect, c(list(e), ev$args))
    c(as.numeric(e$u), as.numeric(e$p))
  }
  x <- c(ev$u, as.numeric(ev$p))
  Rx <- .sens_fd_jac(function(z) reset(z, ev$t), x)
  if (!isTRUE(ev$continuous)) return(Rx)
  Rt <- as.numeric(.sens_fd_jac(function(s) reset(x, s), ev$t))
  cfun <- function(xx, t) {
    q <- pn
    q[] <- xx[n + seq_len(np)]
    e <- .sens_fake_integ(xx[seq_len(n)], q, t)
    e$prob <- prob
    as.numeric(ev$condition(xx[seq_len(n)], t, e))[ev$idx]
  }
  cx <- .sens_fd_grad(function(z) cfun(z, ev$t), x)
  ct <- .sens_fd_grad(function(s) cfun(x, s), ev$t)
  fm <- fns$f(ev$u, ev$p, ev$t)
  pr <- pn
  pr[] <- as.numeric(ev$p_right)
  fp <- fns$f(as.numeric(ev$u_right), pr, ev$t)
  denom <- ct + sum(cx[seq_len(n)] * fm)
  if (!is.finite(denom) || abs(denom) < 1e-12)
    stop("the event at t = ", format(ev$t), " is grazing (the condition does not cross zero ",
         "transversally); its sensitivity is not defined", call. = FALSE)
  tau_x <- -cx / denom
  b <- Rx %*% c(fm, numeric(np)) + Rt - c(fp, numeric(np))
  D <- Rx + b %*% t(tau_x)
  attr(D, "tau_x") <- tau_x
  D
}

# ---------------------------------------------------------------------------
# Extraction

#' @rdname ODEForwardSensitivityProblem
#' @param sol A solution of an `ODEForwardSensitivityProblem`.
#' @param i Index of a saved time; `NULL` for all saved times.
#' @param asmatrix Return the sensitivities as one n x m matrix (per time)
#'   instead of a list with one vector (or time series) per direction.
#' @param t Times to interpolate at instead of `i`.
#' @export
extract_local_sensitivities <- function(sol, i = NULL, asmatrix = FALSE, t = NULL) {
  env <- .sol_env(sol)
  info <- env$prob$sensitivity
  if (is.null(info)) stop("the solution is not of an ODEForwardSensitivityProblem", call. = FALSE)
  n <- info$n; m <- info$m
  un <- colnames(env$U)[seq_len(n)]
  split <- function(z) {
    S <- matrix(z[n + seq_len(n * m)], n, m, dimnames = list(un, info$theta))
    u <- stats::setNames(as.numeric(z[seq_len(n)]), un)
    if (asmatrix) list(u = u, du = S) else
      list(u = u, du = stats::setNames(lapply(seq_len(m), function(j) stats::setNames(S[, j], un)), info$theta))
  }
  if (!is.null(i)) return(split(env$U[i, ]))
  if (!is.null(t)) {
    Z <- sol(t)
    if (!is.matrix(Z)) return(split(Z))
    return(lapply(seq_len(nrow(Z)), function(k) split(Z[k, ])))
  }
  U <- env$U
  u <- U[, seq_len(n), drop = FALSE]
  if (asmatrix) return(list(u = u, du = lapply(seq_len(nrow(U)), function(k) split(U[k, ])$du)))
  du <- lapply(seq_len(m), function(j) {
    D <- U[, n * j + seq_len(n), drop = FALSE]
    colnames(D) <- un
    D
  })
  list(u = u, du = stats::setNames(du, info$theta))
}
