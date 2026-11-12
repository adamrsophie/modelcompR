# Nonlinear solvers in pure R, after NonlinearSolve.jl: NewtonRaphson and
# TrustRegion (lib/NonlinearSolveFirstOrder raphson.jl and trust_region.jl:
# the dogleg subproblem with the radius update schemes Simple, NLsolve, Fan
# and Bastin of GenericTrustRegionScheme), LevenbergMarquardt, GaussNewton,
# the quasi Newton methods Broyden and Klement (lib/NonlinearSolveQuasiNewton
# broyden.jl, klement.jl, solve.jl, initialization.jl, reset_conditions.jl),
# the polyalgorithms of src/poly_algs.jl, lib/NonlinearSolveFirstOrder/src/
# poly_algs.jl and src/default.jl, the default termination of
# NonlinearSolveBase/src/termination_conditions.jl (AbsNormSafeBest), the
# bracketing methods of lib/BracketingNonlinearSolve (itp.jl, bisection.jl,
# brent.jl, falsi.jl), and SSRootfind / DynamicSS of SteadyStateDiffEq.jl.
#
# An algorithm is an S3 list of class c("<Name>", "NonlinearAlgorithm") (or
# "IntervalAlgorithm", "SteadyStateAlgorithm"). The numerical kernels take a
# residual F(x), an optional Jacobian J(x) and return
# list(u, resid, retcode, iters).

.nl_alg <- function(name, ..., cls = "NonlinearAlgorithm") {
  structure(list(name = name, ...), class = c(name, cls))
}

#' @export
print.NonlinearAlgorithm <- function(x, ...) {
  extra <- if (!is.null(x$radius_update_scheme) && x$radius_update_scheme != "Simple")
    paste0("radius_update_scheme = ", x$radius_update_scheme)
  else if (!is.null(x$linesearch)) paste0("linesearch = ", x$linesearch)
  cat(x$name, "(", extra, ")\n", sep = "")
  invisible(x)
}
#' @export
print.IntervalAlgorithm <- function(x, ...) print.NonlinearAlgorithm(x, ...)
#' @export
print.SteadyStateAlgorithm <- function(x, ...) print.NonlinearAlgorithm(x, ...)

#' Nonlinear solver algorithms
#'
#' Algorithms for [solve()] on nonlinear problems, named as in
#' NonlinearSolve.jl.
#'
#' * `NewtonRaphson()` is Newton's method; `linesearch = "BackTracking"`
#'   adds a backtracking line search on the residual norm.
#'   `SimpleNewtonRaphson()` is the same iteration without options.
#' * `TrustRegion()` is a dogleg trust region method. `radius_update_scheme`
#'   selects how the radius changes, with the constants of NonlinearSolve's
#'   `RadiusUpdateSchemes`: `"Simple"` (the default here), `"NLsolve"`,
#'   `"Fan"` and `"Bastin"`. `nlsolve_update_rule = TRUE` is
#'   `radius_update_scheme = "NLsolve"`. NonlinearSolve's own default solves
#'   the subproblem with More's method, which is not ported; the dogleg
#'   subproblem is its `TrustRegionDogleg()`.
#' * `LevenbergMarquardt()` is the damped least squares method, and
#'   `GaussNewton()` the undamped one; both also solve
#'   `NonlinearLeastSquaresProblem`s.
#' * `Broyden()` is Broyden's quasi Newton method: it starts from the
#'   Jacobian `alpha * I` (or the true Jacobian with `init_jacobian =
#'   "true_jacobian"`), updates the inverse with the good (or bad) Broyden
#'   rule and starts over when neither the state nor the residual changed for
#'   three steps. `Klement()` keeps a diagonal Jacobian with Klement's update
#'   and starts over when a diagonal entry vanishes. Neither needs Jacobian
#'   evaluations, and both can succeed where a singular Jacobian stops
#'   Newton's method.
#' * `FastShortcutNonlinearPolyalg()` tries `Broyden`, `Klement`, Newton,
#'   the trust region methods (Simple and Fan) and Levenberg Marquardt in
#'   turn, until one succeeds. With `u0_len <= 25` it starts at Newton, and
#'   with `must_use_jacobian = TRUE` (a problem with a Jacobian function) it
#'   leaves out the quasi Newton methods, as NonlinearSolve's default. It is
#'   the default for `NonlinearProblem`s, built with the problem size and
#'   whether the problem has a Jacobian.
#' * `RobustMultiNewton()` tries the trust region methods with the Simple,
#'   Bastin, NLsolve and Fan radius updates and Newton with and without line
#'   search.
#' * `FastShortcutNLLSPolyalg()` is the default for least squares problems
#'   (Gauss Newton, Levenberg Marquardt, trust region, Gauss Newton with line
#'   search).
#'
#' Solver options given to `solve()`: `abstol` (tolerance on the largest
#' absolute residual, default `eps^0.8`), `reltol`, `maxiters` (default
#' 1000).
#'
#' @param linesearch `NULL` or `"BackTracking"`.
#' @param nlsolve_update_rule Use the trust region radius update of NLsolve.
#' @param max_trust_radius,initial_trust_radius Trust region radii; chosen
#'   from the initial residual when `NULL`.
#' @param radius_update_scheme `"Simple"`, `"NLsolve"`, `"Fan"` or
#'   `"Bastin"`.
#' @param damping_initial Initial damping of Levenberg Marquardt.
#' @param max_resets Largest number of Jacobian resets of the quasi Newton
#'   methods before giving up.
#' @param reset_tolerance Tolerance of the reset test of `Broyden()`;
#'   `eps^(3/4)` when `NULL`.
#' @param init_jacobian `"identity"` (scaled by `alpha`), `"true_jacobian"`
#'   or, for `Klement()`, `"true_jacobian_diagonal"`.
#' @param alpha Scale of the initial identity Jacobian; `NULL` uses
#'   `2 |f(u0)| / max(|u0|, 1)`.
#' @param update_rule `"good_broyden"`, `"bad_broyden"` or `"diagonal"`.
#' @param u0_len Problem size; the polyalgorithm starts at Newton when it is
#'   at most 25.
#' @param must_use_jacobian Leave out the methods that do not use the
#'   Jacobian.
#' @return An algorithm object.
#' @examples
#' prob <- NonlinearProblem(function(u, p) c(u[1]^2 + u[2]^2 - p, u[1] - u[2]), c(1, 2), 2)
#' solve(prob, NewtonRaphson())$u
#' solve(prob, TrustRegion(radius_update_scheme = "Fan"))$u
#' solve(prob, LevenbergMarquardt())$u
#' solve(prob, Broyden())$u
#' solve(prob, Klement())$u
#' @name nonlinear_algorithms
NULL

#' @rdname nonlinear_algorithms
#' @export
NewtonRaphson <- function(linesearch = NULL) .nl_alg("NewtonRaphson", linesearch = linesearch)
#' @rdname nonlinear_algorithms
#' @export
SimpleNewtonRaphson <- function() .nl_alg("SimpleNewtonRaphson", linesearch = NULL)
#' @rdname nonlinear_algorithms
#' @export
TrustRegion <- function(nlsolve_update_rule = FALSE, max_trust_radius = NULL, initial_trust_radius = NULL,
                        radius_update_scheme = NULL) {
  scheme <- radius_update_scheme %||% if (isTRUE(nlsolve_update_rule)) "NLsolve" else "Simple"
  scheme <- match.arg(scheme, c("Simple", "NLsolve", "Fan", "Bastin"))
  .nl_alg("TrustRegion", radius_update_scheme = scheme, nlsolve_update_rule = scheme == "NLsolve",
          max_trust_radius = max_trust_radius, initial_trust_radius = initial_trust_radius)
}
#' @rdname nonlinear_algorithms
#' @export
SimpleTrustRegion <- function(nlsolve_update_rule = FALSE, max_trust_radius = NULL, initial_trust_radius = NULL) {
  alg <- TrustRegion(nlsolve_update_rule, max_trust_radius, initial_trust_radius)
  alg$name <- "SimpleTrustRegion"
  alg
}
#' @rdname nonlinear_algorithms
#' @export
LevenbergMarquardt <- function(damping_initial = 1e-3) .nl_alg("LevenbergMarquardt", damping_initial = damping_initial)
#' @rdname nonlinear_algorithms
#' @export
GaussNewton <- function(linesearch = NULL) .nl_alg("GaussNewton", linesearch = linesearch)
#' @rdname nonlinear_algorithms
#' @export
Broyden <- function(max_resets = 100L, linesearch = NULL, reset_tolerance = NULL, init_jacobian = "identity",
                    alpha = NULL, update_rule = "good_broyden") {
  init_jacobian <- match.arg(init_jacobian, c("identity", "true_jacobian"))
  update_rule <- match.arg(update_rule, c("good_broyden", "bad_broyden", "diagonal"))
  if (update_rule == "diagonal" && init_jacobian != "identity")
    stop("Unknown combination of `init_jacobian = \"", init_jacobian, "\"` and `update_rule = \"diagonal\"`. ",
         "Please choose a valid combination.", call. = FALSE)
  .nl_alg("Broyden", max_resets = as.integer(max_resets), linesearch = linesearch,
          reset_tolerance = reset_tolerance, init_jacobian = init_jacobian, alpha = alpha,
          update_rule = update_rule)
}
#' @rdname nonlinear_algorithms
#' @export
Klement <- function(max_resets = 100L, linesearch = NULL, alpha = NULL, init_jacobian = "identity") {
  init_jacobian <- match.arg(init_jacobian, c("identity", "true_jacobian", "true_jacobian_diagonal"))
  .nl_alg("Klement", max_resets = as.integer(max_resets), linesearch = linesearch, alpha = alpha,
          init_jacobian = init_jacobian)
}
#' @rdname nonlinear_algorithms
#' @export
FastShortcutNonlinearPolyalg <- function(u0_len = NULL, must_use_jacobian = FALSE) {
  if (isTRUE(must_use_jacobian)) {
    algs <- list(NewtonRaphson(), TrustRegion(), TrustRegion(radius_update_scheme = "Fan"), LevenbergMarquardt())
    start <- 1L
  } else {
    algs <- list(Broyden(), Klement(), NewtonRaphson(), TrustRegion(), TrustRegion(radius_update_scheme = "Fan"),
                 LevenbergMarquardt())
    start <- if (!is.null(u0_len) && u0_len <= 25) 3L else 1L
  }
  .nl_alg("FastShortcutNonlinearPolyalg", algs = algs, start_index = start)
}
#' @rdname nonlinear_algorithms
#' @export
RobustMultiNewton <- function() {
  .nl_alg("RobustMultiNewton",
          algs = list(TrustRegion(), TrustRegion(radius_update_scheme = "Bastin"), NewtonRaphson(),
                      NewtonRaphson(linesearch = "BackTracking"), TrustRegion(radius_update_scheme = "NLsolve"),
                      TrustRegion(radius_update_scheme = "Fan")),
          start_index = 1L)
}
#' @rdname nonlinear_algorithms
#' @export
FastShortcutNLLSPolyalg <- function() {
  .nl_alg("FastShortcutNLLSPolyalg",
          algs = list(GaussNewton(), LevenbergMarquardt(), TrustRegion(), GaussNewton(linesearch = "BackTracking"),
                      TrustRegion(radius_update_scheme = "Fan"), LevenbergMarquardt()),
          start_index = 1L)
}

# The algorithm solve() uses without one (src/default.jl): the polyalgorithm
# for the problem size, leaving out the quasi Newton methods when the problem
# has a Jacobian function.
.nl_default_alg <- function(u0, jac = NULL) {
  FastShortcutNonlinearPolyalg(u0_len = length(u0), must_use_jacobian = !is.null(jac))
}

# ---------------------------------------------------------------------------
# Numerical kernels

.nl_default_tol <- function() .Machine$double.eps^0.8

# Least squares solution of J dx = rhs: QR, and the minimum norm solution by
# SVD when J is rank deficient or wide.
.nl_lstsq <- function(J, rhs) {
  n <- ncol(J)
  if (!n) return(numeric())
  if (nrow(J) >= n) {
    q <- qr(J)
    if (q$rank == n) {
      dx <- tryCatch(qr.coef(q, rhs), error = function(e) NULL)
      if (!is.null(dx) && all(is.finite(dx))) return(dx)
    }
  }
  s <- svd(J)
  tol <- max(dim(J)) * .Machine$double.eps * max(s$d, 0)
  keep <- s$d > tol
  if (!any(keep)) return(numeric(n))
  as.numeric(s$v[, keep, drop = FALSE] %*% ((t(s$u[, keep, drop = FALSE]) %*% rhs) / s$d[keep]))
}

.nl_jacobian <- function(F, J, x, fx) {
  if (!is.null(J)) {
    Jx <- J(x)
    if (!is.matrix(Jx)) Jx <- matrix(Jx, nrow = length(fx))
    return(Jx)
  }
  .fd_jacobian(F, x, fx)
}

.nl_project <- function(x, lb, ub) {
  if (!is.null(lb)) x <- pmax(x, lb)
  if (!is.null(ub)) x <- pmin(x, ub)
  x
}

.nl_result <- function(x, fx, retcode, iters) list(u = x, resid = fx, retcode = retcode, iters = iters)

# Termination shared by the kernels: a root is found when the largest
# absolute residual is below abstol. For least squares problems a stationary
# point (small projected gradient) also ends the iteration successfully.
.nl_converged <- function(fx, abstol) all(is.finite(fx)) && max(abs(fx), 0) <= abstol

.nl_stationary <- function(g, fx, x, lb, ub, rtol = 1e-9) {
  # components pushing against an active bound do not count
  g[!.nl_free(x, g, lb, ub)] <- 0
  max(abs(g), 0) <= rtol * max(1, sqrt(sum(fx^2)))
}

# The free variables of a bound constrained step: a variable sitting on a
# bound that the descent direction -g pushes outward is held fixed, so the
# step for the others is computed with it at the bound.
.nl_free <- function(x, g, lb, ub) {
  free <- rep(TRUE, length(x))
  if (!is.null(lb)) free[x <= lb & g > 0] <- FALSE
  if (!is.null(ub)) free[x >= ub & g < 0] <- FALSE
  free
}

# Solve for a step in the free variables only.
.nl_free_step <- function(Jx, rhs, free) {
  if (all(free)) return(.nl_lstsq(Jx, rhs))
  dx <- numeric(ncol(Jx))
  if (any(free)) dx[free] <- .nl_lstsq(Jx[, free, drop = FALSE], rhs)
  dx
}

# At the floating point limit of the iterates: a root solve succeeded when
# the residual is small relative to where the iteration started; a least
# squares solve when the projected gradient vanishes. Otherwise a stall.
.nl_stall_code <- function(fx, f0norm, abstol, ls, F = NULL, J = NULL, x = NULL, lb = NULL, ub = NULL) {
  if (!ls) return(if (max(abs(fx)) <= max(abstol, 1e-9 * max(1, f0norm))) "Success" else "Stalled")
  if (.nl_converged(fx, abstol)) return("Success")
  Jx <- .nl_jacobian(F, J, x, fx)
  if (all(is.finite(Jx)) && .nl_stationary(as.numeric(crossprod(Jx, fx)), fx, x, lb, ub, rtol = 1e-6)) "Success"
  else "Stalled"
}

# Newton and Gauss Newton, with optional backtracking and a postcondition
# hook applied to every accepted iterate (see limited()).
.nl_newton <- function(F, J, x0, abstol, maxiters, lb = NULL, ub = NULL, linesearch = FALSE,
                       post = NULL, ls = FALSE) {
  x <- .nl_project(as.numeric(x0), lb, ub)
  fx <- F(x)
  if (!all(is.finite(fx))) return(.nl_result(x, fx, "InitialFailure", 0L))
  f0norm <- max(abs(fx), 0)
  best <- list(x = x, fx = fx, n = sqrt(sum(fx^2)))
  stalled <- 0L
  for (it in seq_len(maxiters)) {
    if (.nl_converged(fx, abstol)) return(.nl_result(x, fx, "Success", it - 1L))
    Jx <- .nl_jacobian(F, J, x, fx)
    if (!all(is.finite(Jx))) return(.nl_result(best$x, best$fx, "Unstable", it))
    g <- as.numeric(crossprod(Jx, fx))
    if (ls && .nl_stationary(g, fx, x, lb, ub)) return(.nl_result(x, fx, "Success", it - 1L))
    dx <- .nl_free_step(Jx, -fx, .nl_free(x, g, lb, ub))
    nf <- sqrt(sum(fx^2))
    lam <- 1
    repeat {
      xn <- .nl_project(x + lam * dx, lb, ub)
      if (!is.null(post)) xn <- post(xn, x)
      fn <- F(xn)
      ok <- all(is.finite(fn))
      if (ok && (!linesearch || sqrt(sum(fn^2)) <= (1 - 1e-4 * lam) * nf)) break
      lam <- lam / 2
      if (lam < 1e-10) break
    }
    if (!all(is.finite(fn))) return(.nl_result(best$x, best$fx, "Unstable", it))
    step <- max(abs(xn - x))
    x <- xn; fx <- fn
    nn <- sqrt(sum(fx^2))
    if (nn < best$n) {
      stalled <- if (nn > best$n * (1 - 1e-12)) stalled + 1L else 0L
      best <- list(x = x, fx = fx, n = nn)
    } else stalled <- stalled + 1L
    if (step <= 4 * .Machine$double.eps * max(1, max(abs(x)))) {
      if (.nl_converged(fx, abstol)) return(.nl_result(x, fx, "Success", it))
      return(.nl_result(best$x, best$fx, .nl_stall_code(best$fx, f0norm, abstol, ls, F, J, best$x, lb, ub), it))
    }
    if (stalled > 32L)
      return(.nl_result(best$x, best$fx, .nl_stall_code(best$fx, f0norm, abstol, ls, F, J, best$x, lb, ub), it))
  }
  if (.nl_converged(fx, abstol)) return(.nl_result(x, fx, "Success", maxiters))
  .nl_result(best$x, best$fx, "MaxIters", maxiters)
}

# The dogleg step for the model 0.5 |f + J d|^2 inside radius Delta.
.nl_dogleg <- function(Jx, fx, g, Delta) {
  dN <- .nl_lstsq(Jx, -fx)
  if (sqrt(sum(dN^2)) <= Delta) return(dN)
  lg <- sqrt(sum(g^2))
  if (lg == 0) return(dN * 0)
  Jg <- as.numeric(Jx %*% g)
  dc <- lg^3 / sum(Jg^2)
  if (!is.finite(dc) || dc >= Delta) return(-(Delta / lg) * g)
  dsd <- -(dc / lg) * g
  d <- dN - dsd
  a <- sum(d^2); b <- 2 * sum(dsd * d); cc <- sum(dsd^2) - Delta^2
  tau <- (-b + sqrt(max(b^2 - 4 * a * cc, 0))) / (2 * a)
  dsd + tau * d
}

# The trust region method with the dogleg subproblem and the radius updates
# of GenericTrustRegionScheme (trust_region.jl: max_trust_radius,
# initial_trust_radius, step_threshold, shrink_threshold, expand_threshold,
# shrink_factor, get_parameters and InternalAPI.solve!).
.nl_tr_constants <- function(scheme, x, fnorm) {
  simple <- scheme == "Simple"
  list(max_radius = if (simple) max(fnorm, max(x) - min(x)) else Inf,
       initial = switch(scheme, NLsolve = { nx <- sqrt(sum(x^2)); if (nx > 0) nx else 1 },
                        Bastin = 1, Fan = fnorm^0.99 / 10, NA_real_),
       step = if (scheme == "Bastin") 1 / 20 else 1e-4,
       shrink = if (scheme %in% c("NLsolve", "Bastin")) 1 / 20 else 1 / 4,
       expand = if (scheme %in% c("NLsolve", "Bastin")) 9 / 10 else 3 / 4,
       shrink_factor = switch(scheme, NLsolve = 1 / 2, Bastin = 1 / 20, 1 / 4),
       expand_factor = 2,
       p = switch(scheme, NLsolve = c(1 / 2, 0, 0, 0), Fan = c(1 / 10, 1 / 4, 12, 1e18),
                  Bastin = c(5 / 2, 1 / 4, 0, 0), c(0, 0, 0, 0)))
}

.nl_trust_region <- function(F, J, x0, abstol, maxiters, lb = NULL, ub = NULL, scheme = "Simple",
                             max_radius = NULL, initial_radius = NULL, post = NULL, ls = FALSE) {
  max_shrink <- 32L
  xo <- .nl_project(as.numeric(x0), lb, ub)
  fx <- F(xo)
  if (!all(is.finite(fx))) return(.nl_result(xo, fx, "InitialFailure", 0L))
  if (.nl_converged(fx, abstol)) return(.nl_result(xo, fx, "Success", 0L))
  f0norm <- max(abs(fx), 0)
  k <- .nl_tr_constants(scheme, xo, sqrt(sum(fx^2)))
  Dmax <- max_radius %||% k$max_radius
  Delta <- min(initial_radius %||% if (is.na(k$initial)) Dmax / 11 else k$initial, Dmax)
  p1 <- k$p[1L]
  Jx <- .nl_jacobian(F, J, xo, fx)
  fk <- sum(fx^2) / 2
  g <- as.numeric(crossprod(Jx, fx))
  shrink <- 0L
  for (it in seq_len(maxiters)) {
    if (ls && .nl_stationary(g, fx, xo, lb, ub)) return(.nl_result(xo, fx, "Success", it - 1L))
    free <- .nl_free(xo, g, lb, ub)
    delta <- numeric(length(xo))
    if (any(free)) delta[free] <- .nl_dogleg(Jx[, free, drop = FALSE], fx, g[free], Delta)
    x <- .nl_project(xo + delta, lb, ub)
    if (!is.null(post)) x <- post(x, xo)
    delta <- x - xo
    fnew <- F(x)
    finite <- all(is.finite(fnew))
    fk1 <- if (finite) sum(fnew^2) / 2 else Inf
    # the model predicts fk + pred; after projection on the bounds it may
    # predict an increase
    pred <- sum(delta * g) + sum((Jx %*% delta)^2) / 2
    num <- fk1 - fk
    r <- if (finite && pred < 0) num / pred else -1
    accepted <- r > k$step
    nd <- sqrt(sum(delta^2))
    Jnew <- NULL
    if (scheme == "Simple") {
      if (r < k$shrink) {
        Delta <- k$shrink_factor * Delta
        shrink <- shrink + 1L
      } else {
        shrink <- 0L
        if (r > k$expand && r > k$step) Delta <- k$expand_factor * Delta
      }
    } else if (scheme == "NLsolve") {
      if (r < k$shrink) {
        Delta <- k$shrink_factor * Delta
        shrink <- shrink + 1L
      } else {
        shrink <- 0L
        if (r >= k$expand) Delta <- k$expand_factor * nd
        else if (r >= p1) Delta <- max(Delta, k$expand_factor * nd)
      }
    } else if (scheme == "Fan") {
      if (r < k$shrink) {
        p1 <- p1 * k$p[2L]
        shrink <- shrink + 1L
      } else {
        shrink <- 0L
        if (r > k$expand) p1 <- min(p1 * k$p[3L], k$p[4L])
      }
      Delta <- p1 * sqrt(sum((if (finite) fnew else fx)^2))^0.99
    } else {
      # Bastin: the ratio is recomputed with the model at the trial point
      if (accepted) {
        Jnew <- .nl_jacobian(F, J, x, fnew)
        d1 <- sum(crossprod(Jnew, fnew)^2)
        d2 <- sum(crossprod(Jnew, Jnew %*% delta)^2)
        if (num / (d1 + d2 / 2) >= k$expand) Delta <- p1 * nd
        shrink <- 0L
      } else {
        Delta <- Delta * k$p[2L]
        shrink <- shrink + 1L
      }
    }
    Delta <- min(Delta, Dmax)
    if (accepted) {
      xo <- x; fx <- fnew
      if (.nl_converged(fx, abstol)) return(.nl_result(xo, fx, "Success", it))
      Jx <- Jnew %||% .nl_jacobian(F, J, xo, fx)
      fk <- fk1
      g <- as.numeric(crossprod(Jx, fx))
    }
    if (shrink > max_shrink) {
      code <- .nl_stall_code(fx, f0norm, abstol, ls, F, J, xo, lb, ub)
      if (code == "Stalled") code <- "ShrinkThresholdExceeded"
      return(.nl_result(xo, fx, code, it))
    }
  }
  .nl_result(xo, fx, if (.nl_converged(fx, abstol)) "Success" else "MaxIters", maxiters)
}

# ---------------------------------------------------------------------------
# Quasi Newton methods (NonlinearSolveQuasiNewton)

# initial_jacobian_scaling_alpha (NonlinearSolveBase/src/utils.jl)
.nl_qn_alpha <- function(alpha, x, fx) {
  if (!is.null(alpha)) return(as.numeric(alpha))
  fn <- sqrt(sum(fx^2))
  if (fn < 1e-5) return(1)
  2 * fn / max(sqrt(sum(x^2)), 1)
}

# The inverse used when a Jacobian is (re)built: the minimum norm solution
# for a singular one, as LinearSolve's rescue.
.nl_inverse <- function(A) {
  inv <- tryCatch(solve(A), error = function(e) NULL)
  if (!is.null(inv) && all(is.finite(inv))) return(inv)
  s <- svd(A)
  keep <- s$d > max(dim(A)) * .Machine$double.eps * max(s$d, 0)
  if (!any(keep)) return(matrix(0, ncol(A), nrow(A)))
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

# The default termination of NonlinearSolve (AbsNormSafeBestTerminationMode
# with the maximum norm, patience 100 steps, multiplier 3, min_max_factor
# 1.3 and 32 stalled steps); returns NULL to continue or a return code.
.nl_safe_best <- function(state, fx, x, xprev, abstol, ls = FALSE) {
  obj <- max(abs(fx), 0)
  if (!is.finite(obj)) return("Unstable")
  if (obj < state$best_obj) {
    state$best_obj <- obj
    state$best_x <- x
    state$best_f <- fx
  }
  if (obj <= abstol) return("Success")
  state$n <- state$n + 1L
  state$objs[(state$n - 1L) %% 100L + 1L] <- obj
  if (obj <= 3 * abstol && state$n > 100L) {
    o <- state$objs[seq_len(min(state$n, 100L))]
    if (min(o) < 1.3 * max(o)) return(if (ls) "StalledSuccess" else "Stalled")
  }
  state$steps[(state$n - 1L) %% 32L + 1L] <- sqrt(sum((x - xprev)^2))
  if (state$n > 32L && max(state$steps) <= abstol) return(if (ls) "StalledSuccess" else "Stalled")
  NULL
}

# Broyden and Klement: QuasiNewtonAlgorithm with NewtonDescent, no
# globalization unless a line search is given, and the reset rules
# NoChangeInStateReset (Broyden) and IllConditionedJacobianReset (Klement).
.nl_quasi_newton <- function(F, J, x0, abstol, maxiters, alg, lb = NULL, ub = NULL, post = NULL, ls = FALSE) {
  x <- .nl_project(as.numeric(x0), lb, ub)
  n <- length(x)
  fx <- F(x)
  if (!all(is.finite(fx))) return(.nl_result(x, fx, "InitialFailure", 0L))
  if (.nl_converged(fx, abstol)) return(.nl_result(x, fx, "Success", 0L))
  broyden <- inherits(alg, "Broyden")
  diag_rule <- broyden && alg$update_rule == "diagonal"
  # the Jacobian (Klement) or its inverse (Broyden): a vector when diagonal
  mf <- length(fx)
  if (mf != n && (!broyden || diag_rule || alg$init_jacobian == "true_jacobian_diagonal"))
    stop(alg$name, "() with a diagonal Jacobian needs as many residuals as unknowns", call. = FALSE)
  init <- function(x, fx) {
    if (alg$init_jacobian == "identity") {
      a <- .nl_qn_alpha(alg$alpha, x, fx)
      # the inverse of alpha I, n x length(fx) for a rectangular residual
      if (broyden && !diag_rule) return(diag(1 / a, n, mf))
      return(rep(if (broyden) 1 / a else a, n))
    }
    Jt <- .nl_jacobian(F, J, x, fx)
    if (alg$init_jacobian == "true_jacobian_diagonal") return(diag(Jt))
    if (broyden) .nl_inverse(Jt) else Jt
  }
  M <- init(x, fx)
  full <- is.matrix(M)
  tol <- alg$reset_tolerance %||% .Machine$double.eps^(3 / 4)
  cond_max <- 1 / sqrt(.Machine$double.eps)
  st <- new.env(parent = emptyenv())
  st$best_obj <- max(abs(fx)); st$best_x <- x; st$best_f <- fx
  st$n <- 0L; st$objs <- numeric(100L); st$steps <- numeric(32L)
  fprev_rule <- fx                 # the update rule's previous residual
  fprev_reset <- fx                # the reset rule's previous residual
  n_du <- 0L; n_dfu <- 0L
  nresets <- 0L; since_reset <- 0L
  force <- FALSE
  du <- numeric(n)
  linesearch <- identical(alg$linesearch, "BackTracking")
  # SafeBest termination returns the best iterate, except at a root
  finish <- function(code, it) {
    if (code == "Success") return(.nl_result(x, fx, code, it))
    .nl_result(st$best_x, st$best_f, code, it)
  }
  for (it in seq_len(maxiters)) {
    if (it > 1L) {
      reinit <- force
      if (!force) {
        if (broyden) {
          # NoChangeInStateReset with nsteps = 3, checking du and dfu
          if (all(abs(du) <= tol)) {
            n_du <- n_du + 1L
            if (n_du >= 3L) { n_du <- 0L; n_dfu <- 0L; reinit <- TRUE }
          } else { n_du <- 0L; n_dfu <- 0L }
          if (!reinit) {
            if (all(abs(fx - fprev_reset) <= tol)) {
              n_dfu <- n_dfu + 1L
              if (n_dfu >= 3L) { n_dfu <- 0L; n_du <- 0L; reinit <- TRUE }
            } else { n_dfu <- 0L; n_du <- 0L }
            fprev_reset <- fx
          }
        } else {
          reinit <- if (full) {
            k <- tryCatch(kappa(M, exact = TRUE), error = function(e) Inf)
            !is.finite(k) || k >= cond_max
          } else any(M == 0)
        }
      }
      force <- FALSE
      if (reinit) {
        nresets <- nresets + 1L
        if (nresets >= alg$max_resets) return(finish("ConvergenceFailure", it))
        M <- init(x, fx)
        since_reset <- 0L
      } else since_reset <- since_reset + 1L
    }
    # NewtonDescent: the step from the inverse, or a solve with the Jacobian
    du <- if (broyden) {
      if (full) -as.numeric(M %*% fx) else -M * fx
    } else if (full) {
      tryCatch(-as.numeric(solve(M, fx)), error = function(e) NULL)
    } else -fx / M
    if (is.null(du) || !all(is.finite(du))) {
      # the linear solve failed: retry with a fresh Jacobian
      if (since_reset == 0L) return(finish("InternalLinearSolveFailed", it))
      force <- TRUE
      du <- numeric(n)
      next
    }
    lam <- 1
    xn <- .nl_project(x + du, lb, ub)
    if (!is.null(post)) xn <- post(xn, x)
    fn <- F(xn)
    if (linesearch) {
      nf <- sqrt(sum(fx^2))
      while (!(all(is.finite(fn)) && sqrt(sum(fn^2)) <= (1 - 1e-4 * lam) * nf) && lam > 1e-10) {
        lam <- lam / 2
        xn <- .nl_project(x + lam * du, lb, ub)
        if (!is.null(post)) xn <- post(xn, x)
        fn <- F(xn)
      }
      if (lam <= 1e-10 && since_reset > 5L) {
        force <- TRUE
        next
      }
    }
    xprev <- x
    du <- xn - x
    x <- xn; fx <- fn
    code <- .nl_safe_best(st, fx, x, xprev, abstol, ls)
    if (!is.null(code)) return(finish(code, it))
    # the Jacobian update with the step actually taken
    dfu <- fx - fprev_rule
    if (broyden && full) {
      Jdfu <- as.numeric(M %*% dfu)
      if (alg$update_rule == "good_broyden") {
        rmul <- as.numeric(crossprod(M, du))
        denom <- sum(du * Jdfu)
      } else {
        rmul <- dfu
        denom <- sum(dfu^2)
      }
      M <- M + outer((du - Jdfu) / (if (denom == 0) 1e-5 else denom), rmul)
    } else if (broyden) {
      # diagonal Broyden on the inverse
      Jdfu <- M * dfu * du
      denom <- sum(Jdfu)
      M <- M + (du - Jdfu) * du * M / (if (denom == 0) 1e-5 else denom)
    } else if (full) {
      Jdu_cache <- as.numeric((t(M)^2) %*% (du^2))
      Jdu <- as.numeric(M %*% du)
      r <- (fx - fprev_rule - Jdu) / ifelse(Jdu_cache == 0, 1e-5, Jdu_cache)
      M <- M + (outer(r, du) * M) %*% M
    } else {
      Jdu <- M^2 * du^2
      M <- M + ((fx - fprev_rule - M * du) / ifelse(Jdu == 0, 1e-5, Jdu)) * du * M^2
    }
    fprev_rule <- fx
  }
  finish(if (.nl_converged(fx, abstol)) "Success" else "MaxIters", maxiters)
}

# Levenberg Marquardt with Nielsen's damping update and column scaling.
.nl_levenberg_marquardt <- function(F, J, x0, abstol, maxiters, lb = NULL, ub = NULL, lambda0 = 1e-3,
                                    post = NULL, ls = TRUE) {
  x <- .nl_project(as.numeric(x0), lb, ub)
  fx <- F(x)
  if (!all(is.finite(fx))) return(.nl_result(x, fx, "InitialFailure", 0L))
  f0norm <- max(abs(fx), 0)
  n <- length(x)
  Jx <- .nl_jacobian(F, J, x, fx)
  D <- pmax(sqrt(colSums(Jx^2)), 1e-8)
  lam <- lambda0
  nu <- 2
  small <- 0L
  for (it in seq_len(maxiters)) {
    if (.nl_converged(fx, abstol)) return(.nl_result(x, fx, "Success", it - 1L))
    g <- as.numeric(crossprod(Jx, fx))
    if (.nl_stationary(g, fx, x, lb, ub))
      return(.nl_result(x, fx, if (ls) "Success" else .nl_stall_code(fx, f0norm, abstol, FALSE), it - 1L))
    free <- .nl_free(x, g, lb, ub)
    dx <- .nl_free_step(rbind(Jx, diag(sqrt(lam) * D, n)), c(-fx, numeric(n)), free)
    xn <- .nl_project(x + dx, lb, ub)
    if (!is.null(post)) xn <- post(xn, x)
    dx <- xn - x
    fn <- F(xn)
    cost <- sum(fx^2); costn <- if (all(is.finite(fn))) sum(fn^2) else Inf
    pred <- cost - sum((fx + Jx %*% dx)^2)
    rho <- if (is.finite(costn) && pred > 0) (cost - costn) / pred else -1
    if (rho > 0) {
      rel <- (cost - costn) / max(cost, .Machine$double.xmin)
      step <- max(abs(dx))
      x <- xn; fx <- fn
      Jx <- .nl_jacobian(F, J, x, fx)
      D <- pmax(D, sqrt(colSums(Jx^2)))
      lam <- lam * max(1 / 3, 1 - (2 * rho - 1)^3)
      nu <- 2
      small <- if (rel < 1e-14 || step <= 4 * .Machine$double.eps * max(1, max(abs(x)))) small + 1L else 0L
      if (small >= 3L) return(.nl_result(x, fx, .nl_stall_code(fx, f0norm, abstol, ls, F, J, x, lb, ub), it))
    } else {
      lam <- lam * nu
      nu <- 2 * nu
      if (lam > 1e16) return(.nl_result(x, fx, .nl_stall_code(fx, f0norm, abstol, ls, F, J, x, lb, ub), it))
    }
  }
  .nl_result(x, fx, if (.nl_converged(fx, abstol)) "Success" else "MaxIters", maxiters)
}

# Dispatch on the algorithm object. `ls` marks a least squares problem.
.nl_run <- function(alg, F, J, x0, opts, lb = NULL, ub = NULL, ls = FALSE) {
  abstol <- opts$abstol %||% .nl_default_tol()
  maxiters <- as.integer(opts$maxiters %||% 1000L)
  post <- opts$postcondition
  if (!length(x0)) {
    fx <- F(numeric())
    return(.nl_result(numeric(), fx, if (ls || .nl_converged(fx, max(abstol, 1e-9))) "Success" else "Failure", 0L))
  }
  if (!is.null(alg$algs)) {
    best <- NULL
    algs <- alg$algs
    if (!is.null(alg$start_index)) algs <- algs[alg$start_index:length(algs)]
    for (a in algs) {
      res <- tryCatch(.nl_run(a, F, J, x0, opts, lb, ub, ls), error = function(e) NULL)
      if (is.null(res)) next
      res$alg <- a$name
      if (identical(res$retcode, "Success")) return(res)
      if (is.null(best) || (all(is.finite(res$resid)) && sum(res$resid^2) < sum(best$resid^2))) best <- res
    }
    if (is.null(best)) stop("every algorithm of ", alg$name, " failed", call. = FALSE)
    return(best)
  }
  switch(class(alg)[1L],
    NewtonRaphson = , SimpleNewtonRaphson = , GaussNewton =
      .nl_newton(F, J, x0, abstol, maxiters, lb, ub, linesearch = identical(alg$linesearch, "BackTracking"),
                 post = post, ls = ls),
    TrustRegion = , SimpleTrustRegion =
      .nl_trust_region(F, J, x0, abstol, maxiters, lb, ub,
                       scheme = alg$radius_update_scheme %||% (if (isTRUE(alg$nlsolve_update_rule)) "NLsolve" else "Simple"),
                       max_radius = alg$max_trust_radius, initial_radius = alg$initial_trust_radius,
                       post = post, ls = ls),
    Broyden = , Klement =
      .nl_quasi_newton(F, J, x0, abstol, maxiters, alg, lb, ub, post = post, ls = ls),
    LevenbergMarquardt =
      .nl_levenberg_marquardt(F, J, x0, abstol, maxiters, lb, ub, lambda0 = alg$damping_initial %||% 1e-3,
                              post = post, ls = ls),
    if (!is.null(alg$run)) alg$run(F, J, x0, opts, lb, ub, ls)
    else stop("not a nonlinear algorithm: ", class(alg)[1L], call. = FALSE))
}

.nl_problem_fns <- function(prob) {
  p <- prob$p
  f <- prob$f$f
  jac <- prob$f$jac
  list(F = function(x) as.numeric(f(x, p)), J = if (is.null(jac)) NULL else function(x) jac(x, p))
}

.nl_opts <- function(prob, dots) {
  opts <- utils::modifyList(prob$kwargs %||% list(), dots)
  # the corrector of limited() takes the current parameters as third argument
  pc <- opts$postcondition
  if (is.function(pc) && length(formals(pc)) >= 3L) {
    p <- prob$p
    opts$postcondition <- function(xnew, xold) pc(xnew, xold, p)
  }
  opts
}

#' Solve nonlinear problems
#'
#' Methods of [base::solve()] for the nonlinear problem types. The solution
#' has `sol$u` (the root, a named vector), `sol$resid` (the residual there),
#' `sol$retcode` (`"Success"`, `"MaxIters"`, `"Stalled"`,
#' `"InitialFailure"`, ...) and supports symbolic indexing: `sol[x]`,
#' `sol[x + y]`, `sol$ps[a]`.
#'
#' @param a A problem.
#' @param b An algorithm; see [nonlinear_algorithms], [interval_algorithms]
#'   and [steady_state_algorithms]. A default is chosen when missing.
#' @param ... Solver options: `abstol`, `reltol`, `maxiters`.
#' @return A solution object.
#' @examples
#' x <- Variable("x"); y <- Variable("y")
#' sys <- complete(System(list(0 ~ x^2 - y, 0 ~ x + y - 2), name = "ns"))
#' sol <- solve(NonlinearProblem(sys, list(x ~ 0.5, y ~ 0.5)), TrustRegion())
#' sol$retcode; sol$u; sol$resid
#' @name solve_nonlinear
#' @export
solve.NonlinearProblem <- function(a, b, ...) {
  alg <- if (missing(b) || is.null(b)) .nl_default_alg(a$u0, a$f$jac) else b
  if (!inherits(alg, "NonlinearAlgorithm")) stop("not a nonlinear algorithm: ", class(alg)[1L], call. = FALSE)
  opts <- .nl_opts(a, list(...))
  fns <- .nl_problem_fns(a)
  res <- .nl_run(alg, fns$F, fns$J, a$u0, opts, a$lb, a$ub, ls = FALSE)
  .nl_new_solution("NonlinearSolution", a, alg, res$u, res$retcode, resid = res$resid,
                   stats = list(nsteps = res$iters, alg = res$alg %||% alg$name))
}

#' @rdname solve_nonlinear
#' @export
solve.NonlinearLeastSquaresProblem <- function(a, b, ...) {
  alg <- if (missing(b) || is.null(b)) FastShortcutNLLSPolyalg() else b
  if (!inherits(alg, "NonlinearAlgorithm")) stop("not a nonlinear algorithm: ", class(alg)[1L], call. = FALSE)
  opts <- .nl_opts(a, list(...))
  fns <- .nl_problem_fns(a)
  res <- .nl_run(alg, fns$F, fns$J, a$u0, opts, a$lb, a$ub, ls = TRUE)
  .nl_new_solution("NonlinearSolution", a, alg, res$u, res$retcode, resid = res$resid,
                   stats = list(nsteps = res$iters, alg = res$alg %||% alg$name))
}

#' @rdname solve_nonlinear
#' @export
solve.SCCNonlinearProblem <- function(a, b, ...) {
  alg <- if (missing(b) || is.null(b)) NULL else b
  res <- .nl_scc_run(a, alg, .nl_opts(a, list(...)))
  .nl_new_solution("NonlinearSolution", a, alg %||% FastShortcutNonlinearPolyalg(), res$u, res$retcode,
                   resid = res$resid, stats = list(nsteps = res$iters, nblocks = length(a$var_sccs)))
}

# Solve the blocks of an SCC problem in order (SCCNonlinearProblem solve of
# SciMLBase/NonlinearSolve): linear blocks by a least squares solve, blocks
# with homotopy nodes by continuation, the others with `alg` or the default
# polyalgorithm for the block.
.nl_scc_run <- function(a, alg, opts) {
  fns <- .nl_problem_fns(a)
  p <- a$p
  u <- as.numeric(a$u0)
  abstol <- opts$abstol %||% .nl_default_tol()
  kinds <- a$scc_kinds %||% rep("nonlinear", length(a$var_sccs))
  retcode <- "Success"
  iters <- 0L
  for (k in seq_along(a$var_sccs)) {
    vi <- a$var_sccs[[k]]; ei <- a$eq_sccs[[k]]
    Fk <- function(z) {
      uu <- u
      uu[vi] <- z
      fns$F(uu)[ei]
    }
    Jk <- if (is.null(fns$J)) NULL else function(z) {
      uu <- u
      uu[vi] <- z
      fns$J(uu)[ei, vi, drop = FALSE]
    }
    res <- NULL
    # the shortcuts do not know bounds; bounded blocks use the nonlinear solver
    bounded <- any(is.finite(c(a$lb[vi], a$ub[vi])))
    if (!bounded && kinds[k] == "linear") {
      # one exact step for an affine block, and one refinement
      z <- u[vi]
      for (r in 1:2) {
        fz <- Fk(z)
        if (!all(is.finite(fz))) break
        z <- z + .nl_lstsq(.nl_jacobian(Fk, Jk, z, fz), -fz)
      }
      fz <- Fk(z)
      if (.nl_converged(fz, max(abstol, 1e-12 * max(1, abs(z)))))
        res <- .nl_result(z, fz, "Success", 1L)
    } else if (!bounded && kinds[k] == "homotopy" && !is.null(a$f$homotopy)) {
      hf <- a$f$homotopy
      at <- function(lam) list(
        F = function(z) {
          uu <- u
          uu[vi] <- z
          as.numeric(hf$f(uu, p, lam))[ei]
        },
        J = function(z) {
          uu <- u
          uu[vi] <- z
          hf$jac(uu, p, lam)[ei, vi, drop = FALSE]
        })
      sw <- if (inherits(alg, "HomotopySweep")) alg else HomotopySweep(inner = alg)
      res <- .nl_sweep(at, u[vi], sw, opts, 0, 1)
    }
    if (is.null(res))
      res <- .nl_run(if (is.null(alg) || inherits(alg, "HomotopySweep")) .nl_default_alg(vi, Jk) else alg,
                     Fk, Jk, u[vi], opts, a$lb[vi], a$ub[vi], ls = FALSE)
    u[vi] <- res$u
    iters <- iters + res$iters
    if (!identical(res$retcode, "Success")) {
      retcode <- res$retcode
      break
    }
  }
  list(u = u, resid = fns$F(u), retcode = retcode, iters = iters)
}

#' @rdname solve_nonlinear
#' @export
solve.LinearProblem <- function(a, b, ...) {
  opts <- utils::modifyList(a$kwargs %||% list(), list(...))
  if (!missing(b) && !is.null(b) || !is.null(opts$assumptions))
    return(.lin_solution(a, .lin_solve(a, if (missing(b)) NULL else b, opts)))
  A <- a$A; rhs <- a$b
  if (!is.matrix(A)) A <- matrix(A, nrow = length(rhs))
  u <- .nl_lstsq(A, rhs)
  resid <- as.numeric(A %*% u - rhs)
  ok <- all(is.finite(u)) && (nrow(A) > ncol(A) || max(abs(resid), 0) <= 1e-8 * max(1, max(abs(rhs), 0)))
  if (!is.null(a$sys)) a$u0 <- stats::setNames(u, .nl_names(unknowns(a$sys)))
  .nl_new_solution(c("LinearSolution", "NonlinearSolution"), a, .nl_alg("QRFactorization", cls = "LinearAlgorithm"),
                   u, if (ok) "Success" else "Failure", resid = resid)
}

# ---------------------------------------------------------------------------
# Steady state

#' Steady state algorithms
#'
#' `SSRootfind(alg)` finds the steady state of a [SteadyStateProblem()] by
#' solving `f(u, p, Inf) = 0` with a nonlinear algorithm. `DynamicSS(alg)`
#' integrates the differential equations with an ODE algorithm until every
#' derivative satisfies `|du| <= abstol + reltol * |u|`, or `tspan` is
#' reached.
#'
#' @param alg A nonlinear algorithm for `SSRootfind`, an ODE algorithm for
#'   `DynamicSS`; defaults are chosen when `NULL`.
#' @param abstol,reltol Tolerances on the derivatives for `DynamicSS`.
#' @param tspan The longest time to integrate.
#' @return An algorithm object.
#' @examples
#' prob <- SteadyStateProblem(function(u, p, t) p[1] - p[2] * u, u0 = 0, p = c(3, 2))
#' solve(prob, SSRootfind())$u
#' solve(prob, DynamicSS())$u
#' @name steady_state_algorithms
#' @export
SSRootfind <- function(alg = NULL) .nl_alg("SSRootfind", alg = alg, cls = "SteadyStateAlgorithm")
#' @rdname steady_state_algorithms
#' @export
DynamicSS <- function(alg = NULL, abstol = 1e-8, reltol = 1e-6, tspan = 1e6) {
  .nl_alg("DynamicSS", alg = alg, abstol = abstol, reltol = reltol, tspan = tspan, cls = "SteadyStateAlgorithm")
}

#' @rdname solve_nonlinear
#' @export
solve.SteadyStateProblem <- function(a, b, ...) {
  alg <- if (missing(b) || is.null(b)) SSRootfind() else b
  if (inherits(alg, "NonlinearAlgorithm")) alg <- SSRootfind(alg)
  if (inherits(alg, "ODEAlgorithm")) alg <- DynamicSS(alg)
  dots <- list(...)
  f <- a$f$f; p <- a$p
  u0 <- a$u0
  if (inherits(alg, "SSRootfind")) {
    F <- function(x) as.numeric(f(x, p, Inf))
    f0 <- F(as.numeric(u0))
    if (!all(is.finite(f0)))
      stop("the right-hand side is not finite at t = Inf; the system depends on time and has no ",
           "steady state to root find. Use DynamicSS().", call. = FALSE)
    jac <- a$f$jac
    J <- if (is.null(jac)) NULL else function(x) jac(x, p, Inf)
    res <- .nl_run(alg$alg %||% .nl_default_alg(u0, J), F, J, u0, .nl_opts(a, dots), ls = FALSE)
    return(.nl_new_solution(c("SteadyStateSolution", "NonlinearSolution"), a, alg, res$u, res$retcode,
                            resid = res$resid, stats = list(nsteps = res$iters), t = Inf))
  }
  if (!inherits(alg, "DynamicSS")) stop("not a steady state algorithm: ", class(alg)[1L], call. = FALSE)
  M <- a$f$mass_matrix
  drows <- if (is.matrix(M)) which(rowSums(abs(M)) != 0) else seq_along(u0)
  abstol <- alg$abstol; reltol <- alg$reltol
  cb <- DiscreteCallback(
    function(u, t, integ) {
      du <- f(u, integ$p, t)
      all(abs(du[drows]) <= abstol + reltol * abs(u[drows]))
    },
    function(integ) terminate(integ), save_positions = c(FALSE, FALSE))
  oprob <- a
  oprob$tspan <- c(a$tspan[1L], a$tspan[1L] + alg$tspan)
  class(oprob) <- c("ODEProblem", "SciMLProblem")
  if (!is.null(oprob$kwargs$callback)) cb <- CallbackSet(oprob$kwargs$callback, cb)
  oprob$kwargs$callback <- NULL
  args <- c(list(oprob, alg$alg %||% .default_ode_alg(oprob)),
            utils::modifyList(list(callback = cb, save_everystep = FALSE, dtmin = 1e-12, maxiters = 1e6), dots))
  osol <- do.call(solve, args)
  u <- utils::tail(osol$u, 1L)[[1L]]
  tend <- utils::tail(osol$t, 1L)
  du <- as.numeric(f(u, p, tend))
  retcode <- if (identical(osol$retcode, "Terminated")) "Success" else
    if (identical(osol$retcode, "Success")) "Failure" else osol$retcode
  .nl_new_solution(c("SteadyStateSolution", "NonlinearSolution"), a, alg, u, retcode, resid = du,
                   stats = osol$stats, extra = list(original = osol), t = tend)
}

# ---------------------------------------------------------------------------
# Bracketing methods

#' Bracketing algorithms for interval problems
#'
#' Root finders for [IntervalNonlinearProblem()], as in
#' BracketingNonlinearSolve.jl. `ITP()` (interpolate, truncate and project)
#' is the default; `Bisection()`, `Brent()` and `Falsi()` (regula falsi) are
#' the classical methods. The residual must have opposite signs at the ends
#' of the interval, otherwise the return code is `"InitialFailure"`.
#'
#' @param scaled_k1,k2,n0 Tuning parameters of ITP.
#' @param exact_left,exact_right Return the left or right end of the final
#'   bracket instead of its midpoint.
#' @return An algorithm object.
#' @examples
#' prob <- IntervalNonlinearProblem(function(u, p) u^3 - p, c(0, 3), p = 8)
#' solve(prob, ITP())$u
#' solve(prob, Bisection())$u
#' solve(prob, Falsi())$u
#' @name interval_algorithms
#' @export
ITP <- function(scaled_k1 = 0.2, k2 = 2, n0 = 10L) {
  if (scaled_k1 < 0) stop("Hyper-parameter k1 should not be negative", call. = FALSE)
  if (n0 < 0) stop("Hyper-parameter n0 should not be negative", call. = FALSE)
  if (!(k2 >= 1 && k2 <= 1.5 + sqrt(5) / 2))
    stop("Hyper-parameter k2 should be between 1 and 1 + phi where phi = 1.618... is the golden ratio",
         call. = FALSE)
  .nl_alg("ITP", scaled_k1 = scaled_k1, k2 = k2, n0 = as.integer(n0), cls = "IntervalAlgorithm")
}
#' @rdname interval_algorithms
#' @export
Bisection <- function(exact_left = FALSE, exact_right = FALSE) {
  .nl_alg("Bisection", exact_left = exact_left, exact_right = exact_right, cls = "IntervalAlgorithm")
}
#' @rdname interval_algorithms
#' @export
Brent <- function() .nl_alg("Brent", cls = "IntervalAlgorithm")
#' @rdname interval_algorithms
#' @export
Falsi <- function() .nl_alg("Falsi", cls = "IntervalAlgorithm")

.nl_bracket_result <- function(u, fu, left, right, retcode, iters = 0L) {
  list(u = u, resid = fu, left = left, right = right, retcode = retcode, iters = iters)
}

.nl_bisection <- function(f, left, right, fl, fr, abstol, maxiters, exact_left = FALSE, exact_right = FALSE) {
  it <- 0L
  while (it < maxiters) {
    it <- it + 1L
    mid <- (left + right) / 2
    if (mid == left || mid == right)
      return(.nl_bracket_result(left, fl, left, right, "FloatingPointLimit", it))
    fm <- f(mid)
    if (abs((right - left) / 2) < abstol) {
      u <- if (exact_left) left else if (exact_right) right else mid
      fu <- if (exact_left) fl else if (exact_right) fr else fm
      return(.nl_bracket_result(u, fu, left, right, "Success", it))
    }
    if (fm == 0) return(.nl_bracket_result(mid, fm, left, right, "Success", it))
    if (sign(fm) == sign(fl)) {
      left <- mid; fl <- fm
    } else {
      right <- mid; fr <- fm
    }
  }
  .nl_bracket_result(left, fl, left, right, "MaxIters", it)
}

.nl_itp <- function(f, left, right, fl, fr, abstol, maxiters, alg) {
  eps_ <- abstol
  k2 <- alg$k2
  span <- right - left
  k1 <- alg$scaled_k1 * span^(1 - k2)
  n0 <- alg$n0
  if (span / 2 > eps_ * .Machine$double.xmax) {
    eps_s <- span / 2 * 2^n0
  } else {
    n_h <- floor(log2(span / (2 * eps_)))
    eps_s <- eps_ * 2^n_h * 2^n0
  }
  for (i in seq_len(maxiters)) {
    span <- right - left
    mid <- (left + right) / 2
    r <- eps_s - span / 2
    x_f <- left + span * (fl / (fl - fr))
    delta <- max(k1 * span^k2, .Machine$double.eps * max(abs(x_f), .Machine$double.xmin))
    d <- mid - x_f
    sgn <- if (d >= 0) 1 else -1
    xt <- if (delta <= abs(d)) x_f + sgn * delta else mid
    xp <- if (abs(xt - mid) <= r) xt else mid - sgn * r
    if (span < 2 * eps_) return(.nl_bracket_result(xt, f(xt), left, right, "Success", i))
    yp <- f(xp)
    yps <- yp * sign(fr)
    if (yps > 0) {
      right <- xp; fr <- yp
    } else if (yps < 0) {
      left <- xp; fl <- yp
    } else return(.nl_bracket_result(xp, yp, left, right, "Success", i))
    eps_s <- eps_s / 2
    if ((left + right) / 2 == left || (left + right) / 2 == right)
      return(.nl_bracket_result(right, fr, left, right, "FloatingPointLimit", i))
  }
  .nl_bracket_result(left, fl, left, right, "MaxIters", maxiters)
}

.nl_falsi <- function(f, left, right, fl, fr, abstol, maxiters) {
  i <- 1L
  while (i <= maxiters) {
    if ((left + right) / 2 == left || (left + right) / 2 == right)
      return(.nl_bracket_result(left, fl, left, right, "FloatingPointLimit", i))
    mid <- (fr * left - fl * right) / (fr - fl)
    if (!is.finite(mid) || mid <= min(left, right) || mid >= max(left, right)) break
    fm <- f(mid)
    if (fm == 0) return(.nl_bracket_result(mid, fm, left, right, "Success", i))
    if (abs((right - left) / 2) < abstol) return(.nl_bracket_result(mid, fm, left, right, "Success", i))
    if (abs(fm) < abstol) {
      right <- mid; fr <- fm
      break
    }
    if (sign(fl) == sign(fm)) {
      fl <- fm; left <- mid
    } else {
      fr <- fm; right <- mid
    }
    i <- i + 1L
  }
  # fallback to bisection, as the Julia implementation
  if (fr == 0) return(.nl_bracket_result(right, fr, left, right, "Success", i))
  if (sign(fl) == sign(fr)) return(.nl_bracket_result(right, fr, left, right, "Success", i))
  res <- .nl_bisection(f, left, right, fl, fr, abstol, max(maxiters - i, 1L))
  res$iters <- res$iters + i
  res
}

# Brent's method (inverse quadratic interpolation, secant and bisection).
.nl_brent <- function(f, left, right, fl, fr, abstol, maxiters) {
  a <- left; b <- right; fa <- fl; fb <- fr
  if (abs(fa) < abs(fb)) {
    tmp <- a; a <- b; b <- tmp
    tmp <- fa; fa <- fb; fb <- tmp
  }
  cc <- a; fc <- fa
  d <- cc
  mflag <- TRUE
  eps <- .Machine$double.eps
  for (i in seq_len(maxiters)) {
    if (fb == 0) return(.nl_bracket_result(b, fb, min(a, b), max(a, b), "Success", i))
    if (abs(b - a) < 2 * abstol || abs(b - a) <= 2 * eps * abs(b))
      return(.nl_bracket_result(b, fb, min(a, b), max(a, b),
                                if (abs(b - a) < 2 * abstol) "Success" else "FloatingPointLimit", i))
    if (fa != fc && fb != fc) {
      s <- a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) +
        cc * fa * fb / ((fc - fa) * (fc - fb))
    } else {
      s <- b - fb * (b - a) / (fb - fa)
    }
    lo <- min((3 * a + b) / 4, b); hi <- max((3 * a + b) / 4, b)
    tol1 <- 2 * eps * abs(b) + abstol
    if (!(s > lo && s < hi) ||
        (mflag && abs(s - b) >= abs(b - cc) / 2) ||
        (!mflag && abs(s - b) >= abs(cc - d) / 2) ||
        (mflag && abs(b - cc) < tol1) ||
        (!mflag && abs(cc - d) < tol1)) {
      s <- (a + b) / 2
      mflag <- TRUE
    } else mflag <- FALSE
    fs <- f(s)
    d <- cc
    cc <- b; fc <- fb
    if (sign(fa) != sign(fs)) {
      b <- s; fb <- fs
    } else {
      a <- s; fa <- fs
    }
    if (abs(fa) < abs(fb)) {
      tmp <- a; a <- b; b <- tmp
      tmp <- fa; fa <- fb; fb <- tmp
    }
  }
  .nl_bracket_result(b, fb, min(a, b), max(a, b), "MaxIters", maxiters)
}

#' @rdname solve_nonlinear
#' @export
solve.IntervalNonlinearProblem <- function(a, b, ...) {
  alg <- if (missing(b) || is.null(b)) ITP() else b
  if (!inherits(alg, "IntervalAlgorithm")) stop("not a bracketing algorithm: ", class(alg)[1L], call. = FALSE)
  opts <- .nl_opts(a, list(...))
  abstol <- opts$abstol %||% .nl_default_tol()
  maxiters <- as.integer(opts$maxiters %||% 1000L)
  p <- a$p
  ff <- a$f$f
  f <- function(x) as.numeric(ff(x, p))
  left <- min(a$uspan); right <- max(a$uspan)
  fl <- f(left); fr <- f(right)
  res <- if (isTRUE(fl == 0)) .nl_bracket_result(left, fl, left, right, "ExactSolutionLeft")
  else if (isTRUE(fr == 0)) .nl_bracket_result(right, fr, left, right, "ExactSolutionRight")
  else if (!isTRUE(alg$no_bracket) && (!is.finite(fl) || !is.finite(fr) || sign(fl) == sign(fr))) {
    warning("The interval is not an enclosing interval, opposite signs at the boundaries are required.",
            call. = FALSE)
    .nl_bracket_result(left, fl, left, right, "InitialFailure")
  } else switch(class(alg)[1L],
    ITP = .nl_itp(f, left, right, fl, fr, abstol, maxiters, alg),
    Bisection = .nl_bisection(f, left, right, fl, fr, abstol, maxiters, alg$exact_left, alg$exact_right),
    Brent = .nl_brent(f, left, right, fl, fr, abstol, maxiters),
    Falsi = .nl_falsi(f, left, right, fl, fr, abstol, maxiters),
    if (!is.null(alg$run)) alg$run(f, left, right, fl, fr, abstol, maxiters)
    else stop("not a bracketing algorithm: ", class(alg)[1L], call. = FALSE))
  .nl_new_solution(c("IntervalNonlinearSolution", "NonlinearSolution"), a, alg, res$u, res$retcode,
                   resid = res$resid, stats = list(nsteps = res$iters),
                   extra = list(left = res$left, right = res$right))
}
