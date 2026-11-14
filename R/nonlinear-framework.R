# The composable nonlinear solvers of NonlinearSolve.jl: descents
# (NonlinearSolveBase/src/descent: NewtonDescent, SteepestDescent,
# DampedNewtonDescent, Dogleg, GeodesicAcceleration, MoreTrustRegionDescent
# with MINPACK's lmpar), the trust region schemes
# (NonlinearSolveFirstOrder/src/trust_region.jl: GenericTrustRegionScheme
# with the Simple, NLsolve, NocedalWright, Hei, Yuan, More, Bastin and Fan
# radius updates), GeneralizedFirstOrderAlgorithm (solve.jl),
# PseudoTransient (pseudo_transient.jl, switched evolution relaxation),
# QuasiNewtonAlgorithm with the Broyden and Klement updates and
# LimitedMemoryBroyden (NonlinearSolveQuasiNewton), GeneralizedDFSane and
# DFSane (NonlinearSolveSpectralMethods), and the simple solvers of
# SimpleNonlinearSolve (SimpleBroyden, SimpleKlement, SimpleDFSane,
# SimpleGaussNewton, SimpleHalley, SimpleLimitedMemoryBroyden).

#' Composable nonlinear solvers
#'
#' The building blocks of NonlinearSolve's first order and quasi Newton
#' methods, and the methods built from them.
#'
#' Descents give the step `du` from the Jacobian `J` and residual `fu`:
#' `NewtonDescent()` solves `J du = -fu` (least squares when `J` is not
#' square), `SteepestDescent()` takes `-J' fu`, `DampedNewtonDescent(
#' initial_damping, damping_fn)` solves `(J'J + lambda D) du = -J' fu`,
#' `Dogleg()` the dogleg path inside a trust region,
#' `MoreTrustRegionDescent(scaling)` the exact trust region step of Moré
#' (MINPACK `lmpar`; `scaling = TrustRegionScaling$Jacobian` scales by the
#' column norms of `J`, `"Auto"` when they differ by more than 20), and
#' `GeodesicAcceleration(descent, finite_diff_step_geodesic, alpha)` adds the
#' second order geodesic correction to another descent.
#'
#' `GeneralizedFirstOrderAlgorithm(descent, linesearch, trustregion)` runs a
#' descent with an optional line search (see [line_searches]) or trust
#' region scheme (`GenericTrustRegionScheme(method)` with `method` one of
#' `RadiusUpdateSchemes`). `TrustRegionDogleg()` is the trust region with the
#' dogleg subproblem and `TrustRegionRobust()` with Moré's subproblem;
#' `PseudoTransient(alpha_initial)` is the damped Newton method with the
#' switched evolution relaxation of Kelley and Keyes.
#'
#' `QuasiNewtonAlgorithm(descent, linesearch, update_rule, initialization,
#' max_resets)` updates an approximate Jacobian (`"good_broyden"`,
#' `"bad_broyden"` or `"klement"`); `LimitedMemoryBroyden(threshold)` keeps
#' the last `threshold` Broyden updates of the inverse Jacobian.
#' `GeneralizedDFSane(linesearch, sigma_min, sigma_max)` and `DFSane()` are
#' the derivative free spectral methods of La Cruz, Martinez and Raydan.
#' `SimpleBroyden()`, `SimpleKlement()`, `SimpleDFSane()`,
#' `SimpleGaussNewton()`, `SimpleLimitedMemoryBroyden()` and `SimpleHalley()`
#' (Halley's third order method, with the second derivatives by finite
#' differences) are the methods of SimpleNonlinearSolve.
#'
#' All accept the `solve()` options `abstol`, `reltol`, `maxiters` and
#' `termination_condition` (see [termination_modes]).
#'
#' @param descent,linesearch,trustregion Components, see Details.
#' @param initial_damping,damping_fn,alpha_initial Damping settings.
#' @param scaling,min_damping_D Scaling of the Moré subproblem.
#' @param finite_diff_step_geodesic,alpha Geodesic acceleration settings.
#' @param method Radius update scheme.
#' @param max_trust_radius,initial_trust_radius,step_threshold,shrink_threshold,expand_threshold,shrink_factor,expand_factor,max_shrink_times
#'   Trust region settings (0 selects the scheme's default).
#' @param update_rule,initialization,max_resets,reset_tolerance,threshold
#'   Quasi Newton settings.
#' @param sigma_min,sigma_max,sigma_1,M,gamma,tau_min,tau_max,n_exp,max_inner_iterations,eta_strategy
#'   Spectral method settings.
#' @param name Name of the composed algorithm.
#' @param du,u,success,linsolve_success,extras Fields of a `DescentResult()`:
#'   the step, the new point, whether the descent and its linear solve
#'   succeeded, and further results.
#' @param autodiff,linsolve,concrete_jac,jvp_autodiff,vjp_autodiff,jacobian_reuse,mass_matrix
#'   Options of the Julia signatures; Jacobians come from the problem or
#'   from finite differences.
#' @param ... Options.
#' @return An algorithm object.
#' @examples
#' prob <- NonlinearProblem(function(u, p) c(u[1]^2 + u[2]^2 - 2, u[1] - u[2]), c(2, 0.5))
#' solve(prob, DFSane())$u
#' solve(prob, LimitedMemoryBroyden())$u
#' solve(prob, GeneralizedFirstOrderAlgorithm(descent = DampedNewtonDescent()))$u
#' solve(prob, PseudoTransient())$u
#' solve(prob, SimpleHalley())$u
#' @name composable_nonlinear
NULL

.desc <- function(name, ...) structure(list(name = name, ...), class = c(name, "AbstractDescentDirection"))

#' @rdname composable_nonlinear
#' @export
NewtonDescent <- function(linsolve = NULL) .desc("NewtonDescent")
#' @rdname composable_nonlinear
#' @export
SteepestDescent <- function(linsolve = NULL) .desc("SteepestDescent")
#' @rdname composable_nonlinear
#' @export
DampedNewtonDescent <- function(linsolve = NULL, initial_damping = 1, damping_fn = NULL)
  .desc("DampedNewtonDescent", initial_damping = initial_damping, damping_fn = damping_fn)
#' @rdname composable_nonlinear
#' @export
Dogleg <- function(linsolve = NULL) .desc("Dogleg")
#' @rdname composable_nonlinear
#' @export
TrustRegionScaling <- list(None = "None", Jacobian = "Jacobian", Auto = "Auto")
#' @rdname composable_nonlinear
#' @export
MoreTrustRegionDescent <- function(linsolve = NULL, scaling = TrustRegionScaling$None, min_damping_D = 1e-8)
  .desc("MoreTrustRegionDescent", scaling = scaling, min_damping_D = min_damping_D)
#' @rdname composable_nonlinear
#' @export
GeodesicAcceleration <- function(descent = NewtonDescent(), finite_diff_step_geodesic = 0.1, alpha = 0.75)
  .desc("GeodesicAcceleration", descent = descent, h = finite_diff_step_geodesic, alpha = alpha)
#' @rdname composable_nonlinear
#' @export
DescentResult <- function(du = NULL, u = NULL, success = TRUE, linsolve_success = TRUE, extras = list())
  structure(list(du = du, u = u, success = success, linsolve_success = linsolve_success, extras = extras),
            class = "DescentResult")

#' @rdname composable_nonlinear
#' @export
RadiusUpdateSchemes <- list(Simple = "Simple", NLsolve = "NLsolve", NocedalWright = "NocedalWright", Hei = "Hei",
                            Yuan = "Yuan", More = "More", Bastin = "Bastin", Fan = "Fan")
#' @rdname composable_nonlinear
#' @export
TrustRegionSubproblem <- list(More = "More", Dogleg = "Dogleg")

#' @rdname composable_nonlinear
#' @export
GenericTrustRegionScheme <- function(method = RadiusUpdateSchemes$Simple, max_trust_radius = 0, initial_trust_radius = 0,
                                     step_threshold = 0, shrink_threshold = 0, expand_threshold = 0,
                                     shrink_factor = 0, expand_factor = 0)
  structure(list(method = method, max_trust_radius = max_trust_radius, initial_trust_radius = initial_trust_radius,
                 step_threshold = step_threshold, shrink_threshold = shrink_threshold,
                 expand_threshold = expand_threshold, shrink_factor = shrink_factor, expand_factor = expand_factor),
            class = "GenericTrustRegionScheme")

#' @rdname composable_nonlinear
#' @export
GeneralizedFirstOrderAlgorithm <- function(descent = NewtonDescent(), linesearch = NULL, trustregion = NULL,
                                           max_shrink_times = 32L, autodiff = NULL, jvp_autodiff = NULL,
                                           vjp_autodiff = NULL, concrete_jac = NULL, jacobian_reuse = NULL,
                                           name = "GeneralizedFirstOrderAlgorithm") {
  a <- .nl_alg(name, descent = descent, linesearch = linesearch, trustregion = trustregion,
               max_shrink_times = as.integer(max_shrink_times))
  a$run <- function(F, J, x0, opts, lb, ub, ls) .nlf_first_order(a, F, J, x0, opts, lb, ub, ls)
  a
}

#' @rdname composable_nonlinear
#' @export
TrustRegionDogleg <- function(subproblem = TrustRegionSubproblem$Dogleg, ...) {
  # [TrustRegion()] with the dogleg subproblem (trust_region.jl)
  if (!(identical(subproblem, TrustRegionSubproblem$Dogleg) || inherits(subproblem, "Dogleg")))
    stop("`TrustRegionDogleg` fixes `subproblem = TrustRegionSubproblem.Dogleg`; use `TrustRegion(; subproblem)` ",
         "to select a different subproblem.", call. = FALSE)
  TrustRegion(subproblem = subproblem, ...)
}

#' @rdname composable_nonlinear
#' @param subproblem The trust region subproblem: a member of
#'   `TrustRegionSubproblem` or a descent.
#' @export
TrustRegionRobust <- function(subproblem = MoreTrustRegionDescent(), ...) {
  # [TrustRegion()] with More's subproblem; its robust sparse linear solver
  # (a column pivoted sparse QR of the augmented system) only changes sparse
  # Jacobians, and the dense path is MINPACK's lmpar either way
  if (!inherits(subproblem, "MoreTrustRegionDescent"))
    stop("`TrustRegionRobust` fixes `subproblem`'s `linsolve` to `RobustTrustRegionLinsolve()`; use ",
         "`TrustRegion(; subproblem)` to select a different subproblem solver.", call. = FALSE)
  TrustRegion(subproblem = subproblem, ...)
}

#' @rdname composable_nonlinear
#' @export
PseudoTransient <- function(alpha_initial = 1e-3, linesearch = NULL, mass_matrix = NULL, concrete_jac = NULL,
                            linsolve = NULL, autodiff = NULL, jvp_autodiff = NULL, vjp_autodiff = NULL,
                            jacobian_reuse = NULL)
{
  d <- DampedNewtonDescent(initial_damping = alpha_initial, damping_fn = "SER")
  d$mass_matrix <- mass_matrix
  GeneralizedFirstOrderAlgorithm(descent = d, linesearch = linesearch, name = "PseudoTransient")
}

# ---------------------------------------------------------------------------
# Descents

# Moré's lmpar: lambda >= 0 with |D du| <= Delta for (J'J + lambda D'D) du = -J'f
.nlf_lmpar <- function(Jx, fx, Delta, Dd) {
  dN <- .nl_lstsq(Jx, -fx)
  nN <- sqrt(sum((Dd * dN)^2))
  if (nN <= 1.1 * Delta) return(list(du = dN, lambda = 0, step_norm = nN))
  JtJ <- crossprod(Jx); Jtf <- as.numeric(crossprod(Jx, fx))
  solve_l <- function(l) as.numeric(solve(JtJ + l * diag(Dd^2, length(Dd)), -Jtf))
  gnorm <- sqrt(sum((Jtf / Dd)^2))
  hi <- gnorm / Delta
  lo <- 0
  l <- max(1e-3 * hi, sqrt(.Machine$double.eps))
  du <- dN * min(1, Delta / max(nN, .Machine$double.xmin))
  for (it in 1:10) {
    dl <- tryCatch(solve_l(l), error = function(e) NULL)
    if (is.null(dl) || !all(is.finite(dl))) { l <- 2 * l + 1e-8; next }
    du <- dl
    nd <- sqrt(sum((Dd * du)^2))
    phi <- nd - Delta
    if (abs(phi) <= 0.1 * Delta) break
    if (phi > 0) lo <- max(lo, l) else hi <- min(hi, l)
    # Newton step on 1/|D du| (Hebden): dphi/dl = -|q|^2 / |D du|
    q <- as.numeric(solve(JtJ + l * diag(Dd^2, length(Dd)), Dd^2 * du))
    dphi <- -sum(du * q) / nd
    ln <- l - (phi / dphi) * (nd / Delta)
    l <- if (is.finite(ln) && ln > lo && ln < hi) ln else max(0.001 * hi, sqrt(lo * hi))
  }
  list(du = du, lambda = l, step_norm = sqrt(sum((Dd * du)^2)))
}

.nlf_descent <- function(d, Jx, fx, x, state, Delta = NULL, F = NULL) {
  switch(class(d)[1L],
    NewtonDescent = list(du = .nl_lstsq(Jx, -fx)),
    SteepestDescent = list(du = -as.numeric(crossprod(Jx, fx))),
    DampedNewtonDescent = {
      lam <- state$lambda %||% d$initial_damping
      if (identical(d$damping_fn, "SER")) {
        # switched evolution relaxation (Kelley and Keyes): solve
        # (J + alpha^-1 M) du = -f (M = I without a mass matrix) with alpha^-1
        # scaled by the ratio of successive residual norms; the first ratio
        # is against the residual at init, so it is 0 / 0 at a root and the
        # step is not finite (Unstable), as upstream
        nf <- sqrt(sum(fx^2))
        inv <- state$inv_alpha %||% (1 / d$initial_damping)
        inv <- inv * nf / (state$res_norm %||% nf)
        state$res_norm <- nf
        state$inv_alpha <- inv
        A <- if (is.null(d$mass_matrix)) Jx + diag(inv, ncol(Jx)) else Jx + inv * as.matrix(d$mass_matrix)
        if (!all(is.finite(A))) return(list(du = rep(NaN, ncol(Jx))))
        return(list(du = .nl_lstsq(A, -fx)))
      }
      JtJ <- crossprod(Jx)
      list(du = as.numeric(solve(JtJ + lam * diag(pmax(diag(JtJ), 1e-8), ncol(Jx)), -crossprod(Jx, fx))))
    },
    Dogleg = list(du = .nl_dogleg(Jx, fx, as.numeric(crossprod(Jx, fx)), Delta %||% Inf)),
    MoreTrustRegionDescent = {
      cn <- sqrt(colSums(Jx^2))
      use <- switch(d$scaling, Jacobian = TRUE, Auto = {
        if (is.null(state$auto)) { nz <- cn[cn > 0]; state$auto <- length(nz) && max(nz) / min(nz) > 20 }
        state$auto
      }, FALSE)
      Dd <- if (use) { state$D <- pmax(state$D %||% 0, cn, sqrt(d$min_damping_D)); state$D } else rep(1, ncol(Jx))
      r <- .nlf_lmpar(Jx, fx, Delta %||% Inf, Dd)
      list(du = r$du, lambda = r$lambda, step_norm = r$step_norm)
    },
    GeodesicAcceleration = {
      v <- .nlf_descent(d$descent, Jx, fx, x, state, Delta, F)$du
      h <- d$h
      fh <- F(x + h * v)
      rv <- (2 / h) * ((fh - fx) / h - as.numeric(Jx %*% v))
      a <- .nl_lstsq(Jx, -rv)
      nv <- sqrt(sum(v^2)); na <- sqrt(sum(a^2))
      if (nv > 0 && 2 * na / nv > d$alpha) list(du = v) else list(du = v + a / 2)
    },
    stop("unknown descent ", class(d)[1L], call. = FALSE))
}

# ---------------------------------------------------------------------------
# Trust region schemes (GenericTrustRegionScheme)

.nlf_tr_init <- function(tr, x, fx, Jx) {
  m <- tr$method
  T <- function(v, dflt) if (is.null(v) || v == 0) dflt else v
  fn <- sqrt(sum(fx^2)); u0n <- sqrt(sum(x^2))
  mtr <- T(tr$max_trust_radius, if (m %in% c("Simple", "NocedalWright")) max(fn, max(x) - min(x)) else Inf)
  itr <- T(tr$initial_trust_radius, switch(m, NLsolve = if (u0n > 0) u0n else 1, More = if (u0n > 1e-4) u0n else 1,
                                           Hei = 1, Bastin = 1, Fan = fn^0.99 / 10, mtr / 11))
  p <- switch(m, NLsolve = c(1 / 2, 0, 0, 0), Hei = c(5, 1 / 10, 15 / 100, 15 / 100), Yuan = c(2, 1 / 6, 6, 0),
              Fan = c(1 / 10, 1 / 4, 12, 1e18), Bastin = c(5 / 2, 1 / 4, 0, 0), c(0, 0, 0, 0))
  if (m == "Yuan") itr <- p[1L] * sqrt(sum(crossprod(Jx, fx)^2))
  list(method = m, max = mtr, radius = itr, p = p,
       step = T(tr$step_threshold, switch(m, Hei = 0, Yuan = 1e-3, More = 1e-3, Bastin = 1 / 20, 1e-4)),
       shrink = T(tr$shrink_threshold, switch(m, Hei = 0, NLsolve = 1 / 20, Bastin = 1 / 20, 1 / 4)),
       expand = T(tr$expand_threshold, switch(m, NLsolve = 9 / 10, Hei = 0, Bastin = 9 / 10, 3 / 4)),
       shrink_factor = T(tr$shrink_factor, switch(m, NLsolve = 1 / 2, Hei = 0, Bastin = 1 / 20, 1 / 4)),
       expand_factor = T(tr$expand_factor, 2), shrink_counter = 0L)
}

.nlf_rfunc <- function(r, c2, M, g1, g2, b)
  if (r >= c2) (2 * (M - 1 - g2) * atan(r - c2) + (1 + g2)) / pi else (1 - g1 - b) * (exp(r - c2) + b / (1 - g1 - b))

# Accept or reject the step du and update the radius; returns the decision.
.nlf_tr_update <- function(st, F, J, x, fx, Jx, du, extras) {
  xn <- x + du
  fn <- F(xn)
  Jdu <- as.numeric(Jx %*% du)
  denom <- sum(du * crossprod(Jx, fx)) + sum(Jdu^2) / 2
  num <- (sum(fn^2) - sum(fx^2)) / 2
  finite <- is.finite(num) && all(is.finite(fn))
  rho <- if (denom < 0 && finite) num / denom else -1
  accepted <- rho > st$step
  nd <- sqrt(sum(du^2))
  m <- st$method
  if (m == "Simple") {
    if (rho < st$shrink) { st$radius <- st$radius * st$shrink_factor; st$shrink_counter <- st$shrink_counter + 1L }
    else { st$shrink_counter <- 0L; if (rho > st$expand && rho > st$step) st$radius <- st$expand_factor * st$radius }
  } else if (m == "NLsolve") {
    if (rho < st$shrink) { st$radius <- st$radius * st$shrink_factor; st$shrink_counter <- st$shrink_counter + 1L }
    else {
      st$shrink_counter <- 0L
      if (rho >= st$expand) st$radius <- st$expand_factor * nd
      else if (rho >= st$p[1L]) st$radius <- max(st$radius, st$expand_factor * nd)
    }
  } else if (m == "More") {
    sn <- extras$step_norm %||% nd; lam <- extras$lambda %||% sn
    if (rho < st$shrink) { st$radius <- st$shrink_factor * min(st$radius, 10 * sn); st$shrink_counter <- st$shrink_counter + 1L }
    else { st$shrink_counter <- 0L; if (rho >= st$expand || lam == 0) st$radius <- st$expand_factor * sn }
  } else if (m == "NocedalWright") {
    if (rho < st$shrink) { st$radius <- st$shrink_factor * nd; st$shrink_counter <- st$shrink_counter + 1L }
    else {
      st$shrink_counter <- 0L
      if (rho > st$expand && abs(nd - st$radius) < 1e-6 * st$radius) st$radius <- st$expand_factor * st$radius
    }
  } else if (m == "Hei") {
    trn <- .nlf_rfunc(rho, st$shrink, st$p[1L], st$p[3L], st$p[4L], st$p[2L]) * nd
    st$shrink_counter <- if (trn < st$radius) st$shrink_counter + 1L else 0L
    st$radius <- trn
  } else if (m == "Yuan") {
    if (rho < st$shrink) { st$p[1L] <- st$p[2L] * st$p[1L]; st$shrink_counter <- st$shrink_counter + 1L }
    else { if (rho >= st$expand && 2 * nd > st$radius) st$p[1L] <- st$p[3L] * st$p[1L]; st$shrink_counter <- 0L }
    xt <- if (finite) xn else x; ft <- if (finite) fn else fx
    st$radius <- st$p[1L] * sqrt(sum(crossprod(.nl_jacobian(F, J, xt, ft), ft)^2))
  } else if (m == "Fan") {
    if (rho < st$shrink) { st$p[1L] <- st$p[1L] * st$p[2L]; st$shrink_counter <- st$shrink_counter + 1L }
    else { st$shrink_counter <- 0L; if (rho > st$expand) st$p[1L] <- min(st$p[1L] * st$p[3L], st$p[4L]) }
    st$radius <- st$p[1L] * sqrt(sum((if (finite) fn else fx)^2))^0.99
  } else if (m == "Bastin") {
    if (rho > st$step) {
      Jn <- .nl_jacobian(F, J, xn, fn)
      d1 <- sum(crossprod(Jn, fn)^2); d2 <- sum(crossprod(Jn, Jn %*% du)^2)
      if (num / (d1 + d2 / 2) >= st$expand) st$radius <- st$p[1L] * nd
      st$shrink_counter <- 0L
    } else { st$radius <- st$radius * st$p[2L]; st$shrink_counter <- st$shrink_counter + 1L }
  }
  st$radius <- min(st$radius, st$max)
  list(accepted = accepted, x = xn, fx = fn, st = st)
}

# ---------------------------------------------------------------------------
# The first order loop

.nlf_term <- function(opts, fx, x, ls, simple = FALSE) {
  nrm <- if (ls) function(v) sqrt(sum(v^2)) else NULL
  # default_termination_mode: Val(:simple) for SimpleNonlinearSolve, Val(:regular) otherwise
  mode <- opts$termination_condition %||%
    (if (simple) AbsNormTerminationMode(nrm) else AbsNormSafeBestTerminationMode(internalnorm = nrm, max_stalled_steps = 32L))
  abstol <- opts$abstol %||% .nl_default_tol()
  reltol <- opts$reltol %||% .nl_default_tol()
  .term_init(mode, fx, x, abstol, reltol, ls)
}

.nlf_finish <- function(tc, x, fx, code, it) {
  r <- .term_result(tc, x, fx, code)
  .nl_result(r$u, r$fu, r$retcode, it)
}

.nlf_first_order <- function(alg, F, J, x0, opts, lb, ub, ls) {
  maxiters <- as.integer(opts$maxiters %||% 1000L)
  x <- .nl_project(as.numeric(x0), lb, ub)
  fx <- F(x)
  if (!all(is.finite(fx))) return(.nl_result(x, fx, "InitialFailure", 0L))
  tc <- .nlf_term(opts, fx, x, ls)
  # no check before the first step: solve! of NonlinearSolveBase steps first
  state <- new.env(parent = emptyenv())
  Jx <- .nl_jacobian(F, J, x, fx)
  st <- if (!is.null(alg$trustregion)) .nlf_tr_init(alg$trustregion, x, fx, Jx)
  for (it in seq_len(maxiters)) {
    if (it > 1L) Jx <- .nl_jacobian(F, J, x, fx)
    if (!all(is.finite(Jx))) return(.nlf_finish(tc, x, fx, "Unstable", it))
    d <- .nlf_descent(alg$descent, Jx, fx, x, state, if (!is.null(st)) st$radius, F)
    du <- d$du
    if (!all(is.finite(du))) return(.nlf_finish(tc, x, fx, "Unstable", it))
    xprev <- x
    if (!is.null(st)) {
      r <- .nlf_tr_update(st, F, J, x, fx, Jx, du, d)
      st <- r$st
      if (r$accepted) { x <- .nl_project(r$x, lb, ub); fx <- if (identical(x, r$x)) r$fx else F(x) }
      else {
        # a rejected step still reaches the termination check (solve.jl
        # step!), with the state unchanged; it overrides a shrink stop
        if (.term_check(tc, fx, x, xprev)) return(.nlf_finish(tc, x, fx, "Terminated", it))
        if (st$shrink_counter > alg$max_shrink_times) return(.nlf_finish(tc, x, fx, "ShrinkThresholdExceeded", it))
        next
      }
    } else {
      a <- 1
      if (!is.null(alg$linesearch)) {
        phi <- function(s) { v <- F(.nl_project(x + s * du, lb, ub)); if (all(is.finite(v))) sum(v^2) / 2 else Inf }
        g <- as.numeric(crossprod(Jx, fx))
        dphi0 <- sum(g * du)
        dphi <- function(s) {
          xs <- x + s * du; fs <- F(xs)
          sum(crossprod(.nl_jacobian(F, J, xs, fs), fs) * du)
        }
        lsol <- .ls_run(alg$linesearch, phi, dphi0, dphi, state, sqrt(sum(du^2)))
        a <- lsol$step_size
        if (identical(lsol$retcode, "Failure") && a == 0) return(.nlf_finish(tc, x, fx, "InternalLineSearchFailed", it))
      }
      x <- .nl_project(x + a * du, lb, ub)
      fx <- F(x)
      if (!all(is.finite(fx))) return(.nlf_finish(tc, xprev, F(xprev), "Unstable", it))
    }
    if (!is.null(opts$postcondition)) { x <- opts$postcondition(x, xprev); fx <- F(x) }
    if (.term_check(tc, fx, x, xprev)) return(.nlf_finish(tc, x, fx, "Terminated", it))
  }
  .nlf_finish(tc, x, fx, "MaxIters", maxiters)
}

# ---------------------------------------------------------------------------
# Quasi Newton methods

#' @rdname composable_nonlinear
#' @export
QuasiNewtonAlgorithm <- function(descent = NewtonDescent(), linesearch = NULL,
                                 update_rule = c("good_broyden", "bad_broyden", "klement"),
                                 initialization = c("identity", "true_jacobian"), max_resets = 100L,
                                 reset_tolerance = NULL, alpha = NULL, name = "QuasiNewtonAlgorithm") {
  a <- .nl_alg(name, descent = descent, linesearch = linesearch, update_rule = match.arg(update_rule),
               initialization = match.arg(initialization), max_resets = as.integer(max_resets),
               reset_tolerance = reset_tolerance, alpha = alpha)
  a$run <- function(F, J, x0, opts, lb, ub, ls) .nlf_quasi_newton(a, F, J, x0, opts, lb, ub, ls)
  a
}

#' @rdname composable_nonlinear
#' @export
LimitedMemoryBroyden <- function(max_resets = 3L, linesearch = NULL, threshold = 10L, reset_tolerance = NULL,
                                 alpha = NULL) {
  a <- .nl_alg("LimitedMemoryBroyden", linesearch = linesearch, threshold = as.integer(threshold),
               max_resets = as.integer(max_resets), reset_tolerance = reset_tolerance, alpha = alpha)
  # the quasi Newton kernel of R/nonlinear-solvers.R runs the default and
  # "BackTracking" line searches; other line search objects run here
  if (!is.null(linesearch) && !identical(linesearch, "BackTracking"))
    a$run <- function(F, J, x0, opts, lb, ub, ls) .nlf_lbroyden(a, F, J, x0, opts, lb, ub, ls)
  a
}

.nlf_alpha0 <- function(alpha, x, fx) alpha %||% max(2 * sqrt(sum(fx^2)) / max(sqrt(sum(x^2)), 1), .Machine$double.eps)

.nlf_quasi_newton <- function(alg, F, J, x0, opts, lb, ub, ls) {
  maxiters <- as.integer(opts$maxiters %||% 1000L)
  x <- as.numeric(x0); fx <- F(x)
  if (!all(is.finite(fx))) return(.nl_result(x, fx, "InitialFailure", 0L))
  tc <- .nlf_term(opts, fx, x, ls)
  # no check before the first step: solve! of NonlinearSolveBase steps first
  init_B <- function() {
    if (alg$initialization == "true_jacobian") return(.nl_jacobian(F, J, x, fx))
    diag(.nlf_alpha0(alg$alpha, x, fx), length(x))
  }
  B <- init_B(); resets <- 0L
  rtol <- alg$reset_tolerance %||% .Machine$double.eps^(3 / 4)
  state <- new.env(parent = emptyenv())
  for (it in seq_len(maxiters)) {
    du <- .nl_lstsq(B, -fx)
    a <- 1
    if (!is.null(alg$linesearch)) {
      phi <- function(s) { v <- F(x + s * du); if (all(is.finite(v))) sum(v^2) / 2 else Inf }
      lsol <- .ls_run(alg$linesearch, phi, sum(crossprod(B, fx) * du), NULL, state, sqrt(sum(du^2)))
      a <- lsol$step_size
    }
    xprev <- x; fprev <- fx
    s <- a * du
    x <- x + s; fx <- F(x)
    if (!all(is.finite(fx))) return(.nlf_finish(tc, xprev, fprev, "Unstable", it))
    if (.term_check(tc, fx, x, xprev)) return(.nlf_finish(tc, x, fx, "Terminated", it))
    y <- fx - fprev
    # NoChangeInStateReset / IllConditionedJacobianReset
    if (max(abs(s)) <= rtol || max(abs(y)) <= rtol || !all(is.finite(B))) {
      resets <- resets + 1L
      if (resets > alg$max_resets) return(.nlf_finish(tc, x, fx, "ConvergenceFailure", it))
      B <- init_B(); next
    }
    B <- switch(alg$update_rule,
      good_broyden = B + outer(y - as.numeric(B %*% s), s) / sum(s^2),
      bad_broyden = { Bi <- solve(B); Bi2 <- Bi + outer(s - as.numeric(Bi %*% y), y) / sum(y^2); solve(Bi2) },
      klement = {
        # Klement's diagonal-weighted update
        Js <- as.numeric(B %*% s)
        w <- colSums(B^2) * s^2
        B + outer(y - Js, s * colSums(B^2)) / max(sum(w), .Machine$double.eps)
      })
    if (kappa(B) > 1e12) { resets <- resets + 1L; if (resets > alg$max_resets) return(.nlf_finish(tc, x, fx, "ConvergenceFailure", it)); B <- init_B() }
  }
  .nlf_finish(tc, x, fx, "MaxIters", maxiters)
}

# Limited memory good Broyden on the inverse Jacobian: H = -alpha^-1 I + U V'.
.nlf_lbroyden <- function(alg, F, J, x0, opts, lb, ub, ls) {
  maxiters <- as.integer(opts$maxiters %||% 1000L)
  x <- as.numeric(x0); fx <- F(x)
  if (!all(is.finite(fx))) return(.nl_result(x, fx, "InitialFailure", 0L))
  tc <- .nlf_term(opts, fx, x, ls)
  # no check before the first step: solve! of NonlinearSolveBase steps first
  n <- length(x); mth <- alg$threshold
  a0 <- 1 / .nlf_alpha0(alg$alpha, x, fx)
  U <- matrix(0, n, 0); V <- matrix(0, n, 0)
  Hv <- function(v) a0 * v + if (ncol(U)) as.numeric(U %*% crossprod(V, v)) else 0
  Htv <- function(v) a0 * v + if (ncol(U)) as.numeric(V %*% crossprod(U, v)) else 0
  resets <- 0L
  rtol <- alg$reset_tolerance %||% .Machine$double.eps^(3 / 4)
  state <- new.env(parent = emptyenv())
  for (it in seq_len(maxiters)) {
    du <- -Hv(fx)
    a <- 1
    if (!is.null(alg$linesearch)) {
      phi <- function(s) { v <- F(x + s * du); if (all(is.finite(v))) sum(v^2) / 2 else Inf }
      a <- .ls_run(alg$linesearch, phi, -sum(fx^2), NULL, state, sqrt(sum(du^2)))$step_size
    }
    xprev <- x; fprev <- fx
    s <- a * du; x <- x + s; fx <- F(x)
    if (!all(is.finite(fx))) return(.nlf_finish(tc, xprev, fprev, "Unstable", it))
    if (.term_check(tc, fx, x, xprev)) return(.nlf_finish(tc, x, fx, "Terminated", it))
    y <- fx - fprev
    if (max(abs(s)) <= rtol || max(abs(y)) <= rtol) {
      resets <- resets + 1L
      if (resets > alg$max_resets) return(.nlf_finish(tc, x, fx, "ConvergenceFailure", it))
      U <- matrix(0, n, 0); V <- matrix(0, n, 0); next
    }
    Hy <- Hv(y)
    vt <- Htv(s)
    den <- sum(s * Hy)
    if (abs(den) < .Machine$double.eps) next
    if (ncol(U) >= mth) { U <- U[, -1L, drop = FALSE]; V <- V[, -1L, drop = FALSE] }
    U <- cbind(U, (s - Hy) / den); V <- cbind(V, vt)
  }
  .nlf_finish(tc, x, fx, "MaxIters", maxiters)
}

# ---------------------------------------------------------------------------
# Spectral methods

#' @rdname composable_nonlinear
#' @export
GeneralizedDFSane <- function(linesearch = RobustNonMonotoneLineSearch(), sigma_min = 1e-10, sigma_max = 1e10,
                              sigma_1 = NULL, name = "GeneralizedDFSane") {
  a <- .nl_alg(name, linesearch = linesearch, sigma_min = sigma_min, sigma_max = sigma_max, sigma_1 = sigma_1)
  a$run <- function(F, J, x0, opts, lb, ub, ls) .nlf_dfsane(a, F, x0, opts, ls)
  a
}

#' @rdname composable_nonlinear
#' @export
DFSane <- function(sigma_min = 1e-10, sigma_max = 1e10, sigma_1 = 1, M = 10L, gamma = 1e-4, tau_min = 0.1,
                   tau_max = 0.5, n_exp = 2L, max_inner_iterations = 100L,
                   eta_strategy = function(fn_1, n, x_n, f_n) fn_1 / n^2)
  GeneralizedDFSane(RobustNonMonotoneLineSearch(gamma = gamma, sigma_1 = sigma_1, M = M, tau_min = tau_min,
                                                tau_max = tau_max, n_exp = n_exp, maxiters = max_inner_iterations,
                                                eta_strategy = eta_strategy),
                    sigma_min, sigma_max, NULL, name = "DFSane")

.nlf_dfsane <- function(alg, F, x0, opts, ls) {
  maxiters <- as.integer(opts$maxiters %||% 1000L)
  x <- as.numeric(x0); fx <- F(x)
  if (!all(is.finite(fx))) return(.nl_result(x, fx, "InitialFailure", 0L))
  tc <- .nlf_term(opts, fx, x, ls)
  # no check before the first step: solve! of NonlinearSolveBase steps first
  lsr <- alg$linesearch
  nexp <- lsr$n_exp %||% 2L
  M <- lsr$M %||% 10L
  fnorm <- function(v) sqrt(sum(v^2))^nexp
  state <- new.env(parent = emptyenv())
  state$fn_1 <- fnorm(fx); state$history <- rep(fnorm(fx), M); state$nsteps <- 0L
  # safe_dot is BLAS ddot, accumulated with fma (R/nonlinear-wrappers.R)
  spectral <- function(s, y, fu) {
    sig <- .wr_dot(s, s) / .wr_dot(s, y)
    if (!is.finite(sig) || !(alg$sigma_min <= abs(sig) && abs(sig) <= alg$sigma_max))
      sig <- min(max(1 / sqrt(sum(fu^2)), 1), 1e5)
    sig
  }
  sigma <- if (is.null(alg$sigma_1)) spectral(x, fx, fx) else alg$sigma_1
  for (it in seq_len(maxiters)) {
    du <- -sigma * fx
    phi <- function(s) { v <- F(x + s * du); if (all(is.finite(v))) sum(v^2) / 2 else Inf }
    lsol <- .ls_run(lsr, phi, NULL, NULL, state, sqrt(sum(du^2)))
    if (!identical(lsol$retcode, "Success")) return(.nlf_finish(tc, x, fx, "InternalLineSearchFailed", it))
    xprev <- x; fprev <- fx
    # axpy!(alpha, du, u)
    x <- .wr_fma(rep_len(lsol$step_size, length(du)), du, x); fx <- F(x)
    if (!is.null(opts$postcondition)) { x <- opts$postcondition(x, xprev); fx <- F(x) }
    if (.term_check(tc, fx, x, xprev)) return(.nlf_finish(tc, x, fx, "Terminated", it))
    sigma <- spectral(x - xprev, fx - fprev, fx)
    # callback_into_cache!: the non monotone history
    state$history[(state$nsteps - 1L) %% M + 1L] <- fnorm(fx)
    state$nsteps <- state$nsteps + 1L
  }
  .nlf_finish(tc, x, fx, "MaxIters", maxiters)
}

# ---------------------------------------------------------------------------
# SimpleNonlinearSolve

#' @rdname composable_nonlinear
#' @export
SimpleBroyden <- function(linesearch = NULL, alpha = NULL) {
  a <- QuasiNewtonAlgorithm(linesearch = linesearch, update_rule = "good_broyden", alpha = alpha, name = "SimpleBroyden")
  a
}
#' @rdname composable_nonlinear
#' @export
SimpleKlement <- function() QuasiNewtonAlgorithm(update_rule = "klement", initialization = "true_jacobian",
                                                 name = "SimpleKlement")
#' @rdname composable_nonlinear
#' @export
SimpleDFSane <- function(...) {
  a <- DFSane(...)
  a$name <- "SimpleDFSane"
  a
}
#' @rdname composable_nonlinear
#' @export
SimpleGaussNewton <- function() {
  a <- GaussNewton()
  a$name <- "SimpleGaussNewton"
  a
}
#' @rdname composable_nonlinear
#' @export
SimpleLimitedMemoryBroyden <- function(threshold = 10L, linesearch = NULL, alpha = NULL) {
  a <- LimitedMemoryBroyden(threshold = threshold, linesearch = linesearch, alpha = alpha)
  a$name <- "SimpleLimitedMemoryBroyden"
  a
}

# Halley's method: u+ = u - (J - 0.5 H[a])^-1 f with a = J^-1 f and H[a]
# the directional second derivative of F along a.
#' @rdname composable_nonlinear
#' @export
SimpleHalley <- function(autodiff = NULL) {
  a <- .nl_alg("SimpleHalley")
  a$run <- function(F, J, x0, opts, lb, ub, ls) {
    maxiters <- as.integer(opts$maxiters %||% 1000L)
    x <- as.numeric(x0); fx <- F(x)
    if (!all(is.finite(fx))) return(.nl_result(x, fx, "InitialFailure", 0L))
    tc <- .nlf_term(opts, fx, x, ls, simple = TRUE)
    if (.term_check(tc, fx, x, x)) return(.nlf_finish(tc, x, fx, "Terminated", 0L))
    for (it in seq_len(maxiters)) {
      Jx <- .nl_jacobian(F, J, x, fx)
      aa <- .nl_lstsq(Jx, fx)
      h <- .Machine$double.eps^(1 / 4) * max(1, sqrt(sum(x^2))) / max(sqrt(sum(aa^2)), 0x1.56e1fc2f8f359p-997)
      Haa <- (F(x + h * aa) - 2 * fx + F(x - h * aa)) / h^2
      bb <- .nl_lstsq(Jx, Haa)
      cc <- aa * aa / (aa - bb / 2)
      cc[!is.finite(cc)] <- aa[!is.finite(cc)]
      xprev <- x
      x <- x - cc; fx <- F(x)
      if (!all(is.finite(fx))) return(.nlf_finish(tc, xprev, F(xprev), "Unstable", it))
      if (.term_check(tc, fx, x, xprev)) return(.nlf_finish(tc, x, fx, "Terminated", it))
    }
    .nlf_finish(tc, x, fx, "MaxIters", maxiters)
  }
  a
}
