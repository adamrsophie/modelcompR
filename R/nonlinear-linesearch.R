# Termination modes of NonlinearSolveBase (termination_conditions.jl) and
# the line searches of LineSearch.jl (no_search.jl, backtracking.jl,
# armijo.jl, golden_section.jl, li_fukushima.jl, robust_non_monotone.jl,
# strong_wolfe.jl, more_thuente.jl, hager_zhang.jl, merit.jl).

# ---------------------------------------------------------------------------
# Termination modes

#' Termination modes of nonlinear solvers
#'
#' Values for the `termination_condition` option of `solve()` on nonlinear
#' problems, as in NonlinearSolve. With `f` the residual and `u` the iterate:
#' `AbsTerminationMode()` stops when every `|f_i| <= abstol`;
#' `RelTerminationMode()` when every `|f_i| <= reltol |u_i + f_i|`;
#' `NormTerminationMode(internalnorm)` when `norm(f) <= abstol` or `norm(f) <=
#' reltol norm(f + u)`; `AbsNormTerminationMode()` and
#' `RelNormTerminationMode()` the norm forms of the absolute and relative
#' tests. The safe variants (`AbsNormSafeTerminationMode()`,
#' `RelNormSafeTerminationMode()`) also stop with `"Unstable"` when the
#' objective grows past `protective_threshold` times its initial value, and
#' with `"Stalled"` when it stays within `min_max_factor` over
#' `patience_steps` steps near the tolerance (or the steps stay below the
#' tolerance for `max_stalled_steps` steps); the safe best variants
#' (`AbsNormSafeBestTerminationMode()`, `RelNormSafeBestTerminationMode()`)
#' return the best iterate seen. The default of the solvers is
#' `AbsNormSafeBestTerminationMode(max_stalled_steps = 32)` with the maximum
#' norm (the 2-norm for least squares problems).
#'
#' @param internalnorm Norm function; the maximum absolute value by default.
#' @param protective_threshold,patience_steps,patience_objective_multiplier,min_max_factor,max_stalled_steps
#'   Settings of the safe modes.
#' @return A termination mode.
#' @examples
#' prob <- NonlinearProblem(function(u, p) u^2 - 2, 1)
#' solve(prob, NewtonRaphson(), termination_condition = AbsTerminationMode())$u
#' @name termination_modes
NULL

.term_mode <- function(kind, internalnorm = NULL, safe = FALSE, best = FALSE, ...)
  structure(list(kind = kind, internalnorm = internalnorm %||% function(x) max(abs(x), 0), safe = safe,
                 best = best, ...), class = c(paste0(kind, "TerminationMode"), "NonlinearTerminationMode"))

.term_safe <- function(kind, best, internalnorm, protective_threshold, patience_steps,
                       patience_objective_multiplier, min_max_factor, max_stalled_steps)
  .term_mode(kind, internalnorm, safe = TRUE, best = best, protective_threshold = protective_threshold,
             patience_steps = as.integer(patience_steps), patience_objective_multiplier = patience_objective_multiplier,
             min_max_factor = min_max_factor, max_stalled_steps = max_stalled_steps)

#' @rdname termination_modes
#' @export
AbsTerminationMode <- function() .term_mode("Abs")
#' @rdname termination_modes
#' @export
RelTerminationMode <- function() .term_mode("Rel")
#' @rdname termination_modes
#' @export
NormTerminationMode <- function(internalnorm = NULL) .term_mode("Norm", internalnorm)
#' @rdname termination_modes
#' @export
AbsNormTerminationMode <- function(internalnorm = NULL) .term_mode("AbsNorm", internalnorm)
#' @rdname termination_modes
#' @export
RelNormTerminationMode <- function(internalnorm = NULL) .term_mode("RelNorm", internalnorm)
#' @rdname termination_modes
#' @export
AbsNormSafeTerminationMode <- function(internalnorm = NULL, protective_threshold = 1000, patience_steps = 100,
                                       patience_objective_multiplier = 3, min_max_factor = 1.3,
                                       max_stalled_steps = NULL)
  .term_safe("AbsNormSafe", FALSE, internalnorm, protective_threshold, patience_steps,
             patience_objective_multiplier, min_max_factor, max_stalled_steps)
#' @rdname termination_modes
#' @export
RelNormSafeTerminationMode <- function(internalnorm = NULL, protective_threshold = 1000, patience_steps = 100,
                                       patience_objective_multiplier = 3, min_max_factor = 1.3,
                                       max_stalled_steps = NULL)
  .term_safe("RelNormSafe", FALSE, internalnorm, protective_threshold, patience_steps,
             patience_objective_multiplier, min_max_factor, max_stalled_steps)
#' @rdname termination_modes
#' @export
AbsNormSafeBestTerminationMode <- function(internalnorm = NULL, protective_threshold = 1000, patience_steps = 100,
                                           patience_objective_multiplier = 3, min_max_factor = 1.3,
                                           max_stalled_steps = NULL)
  .term_safe("AbsNormSafeBest", TRUE, internalnorm, protective_threshold, patience_steps,
             patience_objective_multiplier, min_max_factor, max_stalled_steps)
#' @rdname termination_modes
#' @export
RelNormSafeBestTerminationMode <- function(internalnorm = NULL, protective_threshold = 1000, patience_steps = 100,
                                           patience_objective_multiplier = 3, min_max_factor = 1.3,
                                           max_stalled_steps = NULL)
  .term_safe("RelNormSafeBest", TRUE, internalnorm, protective_threshold, patience_steps,
             patience_objective_multiplier, min_max_factor, max_stalled_steps)

# The termination cache: init with the first residual and iterate.
.term_init <- function(mode, fu, u, abstol, reltol, leastsq = FALSE) {
  tc <- new.env(parent = emptyenv())
  tc$mode <- mode; tc$abstol <- abstol; tc$reltol <- reltol; tc$leastsq <- leastsq
  tc$retcode <- "Default"; tc$nsteps <- 0L; tc$u <- u; tc$best_fu <- fu
  nrm <- mode$internalnorm
  rel <- grepl("^RelNorm", mode$kind)
  if (isTRUE(mode$safe)) {
    tc$initial <- if (rel) nrm(fu) / (nrm(fu + u) + .Machine$double.eps) else nrm(fu)
    tc$best <- tc$initial
    tc$trace <- rep(NA_real_, mode$patience_steps)
    tc$step_trace <- if (!is.null(mode$max_stalled_steps)) rep(NA_real_, mode$max_stalled_steps)
    tc$u0_norm <- sqrt(sum(u^2))
  } else tc$best <- Inf
  tc
}

# TRUE when the iteration should stop; tc$retcode tells why.
.term_check <- function(tc, fu, u, uprev) {
  m <- tc$mode; nrm <- m$internalnorm
  abstol <- tc$abstol; reltol <- tc$reltol
  if (!isTRUE(m$safe)) {
    ok <- switch(m$kind,
      Abs = all(abs(fu) <= abstol),
      Rel = all(abs(fu) <= reltol * abs(u + fu)),
      Norm = { dn <- nrm(fu); dn <= abstol || dn <= reltol * nrm(fu + u) },
      AbsNorm = nrm(fu) <= abstol,
      RelNorm = nrm(fu) <= reltol * nrm(fu + u))
    if (isTRUE(ok)) tc$retcode <- "Success"
    return(isTRUE(ok))
  }
  rel <- grepl("^RelNorm", m$kind)
  obj <- if (rel) nrm(fu) / (nrm(fu + u) + .Machine$double.eps) else nrm(fu)
  crit <- if (rel) reltol else abstol
  if (!is.finite(obj)) { tc$retcode <- "Unstable"; return(TRUE) }
  if (!is.null(m$protective_threshold) && obj > tc$initial * m$protective_threshold * length(fu)) {
    tc$retcode <- "Unstable"; return(TRUE)
  }
  if (isTRUE(m$best) && obj < tc$best) { tc$best <- obj; tc$u <- u; tc$best_fu <- fu }
  if (obj <= crit) { tc$retcode <- "Success"; return(TRUE) }
  tc$nsteps <- tc$nsteps + 1L
  tc$trace[(tc$nsteps - 1L) %% length(tc$trace) + 1L] <- obj
  if (obj <= m$patience_objective_multiplier * crit && tc$nsteps > m$patience_steps) {
    r <- range(tc$trace, na.rm = TRUE)
    if (r[1L] < m$min_max_factor * r[2L]) {
      tc$retcode <- if (tc$leastsq) "StalledSuccess" else "Stalled"; return(TRUE)
    }
  }
  if (!is.null(tc$step_trace)) {
    tc$step_trace[(tc$nsteps - 1L) %% length(tc$step_trace) + 1L] <- sqrt(sum((u - uprev)^2))
    if (tc$nsteps > m$max_stalled_steps) {
      mx <- max(tc$step_trace, na.rm = TRUE)
      stalled <- if (rel) mx <= reltol * (mx + tc$u0_norm) else mx <= abstol
      if (stalled) { tc$retcode <- if (tc$leastsq) "StalledSuccess" else "Stalled"; return(TRUE) }
    }
  }
  tc$retcode <- "Failure"
  FALSE
}

# The final iterate and return code (termination_condition_result).
.term_result <- function(tc, u, fu, solver_retcode) {
  m <- tc$mode
  if (!isTRUE(m$safe)) return(list(u = u, fu = fu, retcode = if (solver_retcode == "Terminated") "Success" else solver_retcode))
  code <- if (solver_retcode == "Terminated") { if (tc$retcode != "Default") tc$retcode else "Success" }
          else if (solver_retcode == "Success") "Failure" else solver_retcode
  if (isTRUE(m$best) && !is.null(tc$u)) list(u = tc$u, fu = tc$best_fu, retcode = code)
  else list(u = u, fu = fu, retcode = code)
}

# ---------------------------------------------------------------------------
# Merits and line searches

#' Line searches
#'
#' The line searches of LineSearch.jl, for the `linesearch` argument of the
#' nonlinear solvers (and the gradient based optimizers). The merit is
#' `ResidualMerit()`, `0.5 |F(u + alpha du)|^2`, for nonlinear problems and
#' `ObjectiveMerit()`, the objective, for optimization.
#'
#' `NoLineSearch(alpha)` takes the step `alpha`. `BackTracking(c_1, rho_hi,
#' rho_lo, order, maxstep, initial_alpha, maxiters)` backtracks from
#' `initial_alpha` by quadratic (`order = 2`) or cubic interpolation until
#' the Armijo condition holds; `ArmijoLineSearch(c_1, contraction,
#' initial_alpha, maxiters)` contracts by a fixed factor, and
#' `ProjectedBackTracking()` does so on the projection of the trial point
#' onto the bounds. `GoldenSection(tol, maxiters)` minimizes the merit on
#' `[0, 1]` by golden section. `LiFukushimaLineSearch()` is the derivative
#' free search of Li and Fukushima and `RobustNonMonotoneLineSearch()` the
#' non monotone search of La Cruz, Martinez and Raydan used by `DFSane()`.
#' `StrongWolfeLineSearch(c1, c2, ...)`, `MoreThuenteLineSearch(ftol, gtol,
#' xtol, ...)` and `HagerZhangLineSearch(delta, sigma, ...)` find a step that
#' satisfies the strong Wolfe conditions (the latter two with the safeguarded
#' cubic steps of More and Thuente and the bracketing of Hager and Zhang).
#' `LineSearchesJL(method)` selects one of these by the name of its
#' LineSearches.jl method (`"BackTracking"`, `"StrongWolfe"`,
#' `"MoreThuente"`, `"HagerZhang"`, `"Static"`). `LineSearchSolution(step_size,
#' retcode)` is the result of a search; `get_trial()` and
#' `set_initial_step()` are its cache interface.
#'
#' @param alpha,initial_alpha,alpha_init,lambda_0 Initial steps.
#' @param c_1,c1,c2,ftol,gtol,xtol,delta,sigma,epsilon,theta,rho,beta,sigma_1,sigma_2,eta,gamma,tau_min,tau_max
#'   Line search constants.
#' @param rho_hi,rho_lo Bounds of the backtracking factor.
#' @param order Interpolation order of `BackTracking()`.
#' @param maxstep,alpha_max,alpha_min Step limits.
#' @param contraction Contraction factor.
#' @param tol Tolerance of the golden section.
#' @param maxiters,nan_maxiters,zoom_maxiters Iteration limits.
#' @param M,n_exp,eta_strategy Settings of the non monotone search.
#' @param method LineSearches.jl method name.
#' @param autodiff Ignored (derivatives are finite differences or exact).
#' @param step_size,retcode Fields of a solution.
#' @param cache A line search cache.
#' @return A line search algorithm or solution.
#' @examples
#' prob <- NonlinearProblem(function(u, p) c(u[1]^2 - 2, u[2] - u[1]), c(5, 5))
#' solve(prob, NewtonRaphson(linesearch = BackTracking()))$u
#' @name line_searches
NULL

.ls_alg <- function(name, ...) structure(list(name = name, ...), class = c(name, "AbstractLineSearchAlgorithm"))

#' @export
print.AbstractLineSearchAlgorithm <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

#' @rdname line_searches
#' @export
ResidualMerit <- function() structure(list(name = "ResidualMerit"), class = "ResidualMerit")
#' @rdname line_searches
#' @export
ObjectiveMerit <- function() structure(list(name = "ObjectiveMerit"), class = "ObjectiveMerit")

#' @rdname line_searches
#' @export
NoLineSearch <- function(alpha = 1) .ls_alg("NoLineSearch", alpha = alpha)
#' @rdname line_searches
#' @export
BackTracking <- function(autodiff = NULL, c_1 = 1e-4, rho_hi = 0.5, rho_lo = 0.1, order = 3L, maxstep = Inf,
                         initial_alpha = 1, maxiters = 1000L) {
  stopifnot(order %in% 2:3)
  .ls_alg("BackTracking", c_1 = c_1, rho_hi = rho_hi, rho_lo = rho_lo, order = as.integer(order), maxstep = maxstep,
          initial_alpha = initial_alpha, maxiters = as.integer(maxiters))
}
#' @rdname line_searches
#' @export
ArmijoLineSearch <- function(autodiff = NULL, c_1 = 1e-4, contraction = 0.5, initial_alpha = 1, maxiters = 40L)
  .ls_alg("ArmijoLineSearch", c_1 = c_1, contraction = contraction, initial_alpha = initial_alpha,
          maxiters = as.integer(maxiters))
#' @rdname line_searches
#' @export
ProjectedBackTracking <- function(c_1 = 1e-4, contraction = 0.5, initial_alpha = 1, maxiters = 40L)
  .ls_alg("ProjectedBackTracking", c_1 = c_1, contraction = contraction, initial_alpha = initial_alpha,
          maxiters = as.integer(maxiters))
#' @rdname line_searches
#' @export
GoldenSection <- function(tol = 1e-7, maxiters = 100L) .ls_alg("GoldenSection", tol = tol, maxiters = as.integer(maxiters))
#' @rdname line_searches
#' @export
LiFukushimaLineSearch <- function(lambda_0 = 1, beta = 0.5, sigma_1 = 0.001, sigma_2 = 0.001, eta = 0.1, rho = 0.9,
                                  nan_maxiters = 5L, maxiters = 100L)
  .ls_alg("LiFukushimaLineSearch", lambda_0 = lambda_0, beta = beta, sigma_1 = sigma_1, sigma_2 = sigma_2, eta = eta,
          rho = rho, nan_maxiters = nan_maxiters, maxiters = as.integer(maxiters))
#' @rdname line_searches
#' @export
RobustNonMonotoneLineSearch <- function(gamma = 1e-4, sigma_1 = 1, M = 10L, tau_min = 0.1, tau_max = 0.5, n_exp = 2L,
                                        maxiters = 100L, eta_strategy = function(fn_1, n, u, f) fn_1 / n^2)
  .ls_alg("RobustNonMonotoneLineSearch", gamma = gamma, sigma_1 = sigma_1, M = as.integer(M), tau_min = tau_min,
          tau_max = tau_max, n_exp = as.integer(n_exp), maxiters = as.integer(maxiters), eta_strategy = eta_strategy)
#' @rdname line_searches
#' @export
StrongWolfeLineSearch <- function(autodiff = NULL, c1 = 1e-4, c2 = 0.9, alpha_init = 1, alpha_max = 65536,
                                  maxiters = 10L, zoom_maxiters = 10L)
  .ls_alg("StrongWolfeLineSearch", c1 = c1, c2 = c2, alpha_init = alpha_init, alpha_max = alpha_max,
          maxiters = as.integer(maxiters), zoom_maxiters = as.integer(zoom_maxiters))
#' @rdname line_searches
#' @export
MoreThuenteLineSearch <- function(autodiff = NULL, ftol = 1e-4, gtol = 0.9, xtol = 1e-10, alpha_init = 1,
                                  alpha_min = 0, alpha_max = Inf, maxiters = 50L)
  .ls_alg("MoreThuenteLineSearch", ftol = ftol, gtol = gtol, xtol = xtol, alpha_init = alpha_init,
          alpha_min = alpha_min, alpha_max = alpha_max, maxiters = as.integer(maxiters))
#' @rdname line_searches
#' @export
HagerZhangLineSearch <- function(autodiff = NULL, delta = 0.1, sigma = 0.9, epsilon = 1e-6, theta = 0.5, rho = 5,
                                 maxiters = 50L, alpha_init = 1, alpha_max = Inf)
  .ls_alg("HagerZhangLineSearch", delta = delta, sigma = sigma, epsilon = epsilon, theta = theta, rho = rho,
          maxiters = as.integer(maxiters), alpha_init = alpha_init, alpha_max = alpha_max)
#' @rdname line_searches
#' @export
LineSearchesJL <- function(method = "BackTracking", initial_alpha = 1, autodiff = NULL) {
  m <- if (is.character(method)) method else method$name %||% "BackTracking"
  ls <- switch(m, BackTracking = BackTracking(initial_alpha = initial_alpha),
               StrongWolfe = StrongWolfeLineSearch(alpha_init = initial_alpha),
               MoreThuente = MoreThuenteLineSearch(alpha_init = initial_alpha),
               HagerZhang = HagerZhangLineSearch(alpha_init = initial_alpha),
               Static = NoLineSearch(initial_alpha),
               stop("unknown LineSearches.jl method ", m, call. = FALSE))
  ls$wrapper <- "LineSearchesJL"
  ls
}
#' @rdname line_searches
#' @export
LineSearchSolution <- function(step_size, retcode = "Success")
  structure(list(step_size = step_size, retcode = retcode), class = "LineSearchSolution")

#' @rdname line_searches
#' @export
get_trial <- function(cache) list(u = cache$u_trial, fu = cache$fu_trial)

#' @rdname line_searches
#' @export
set_initial_step <- function(cache, alpha) {
  if (!(is.finite(alpha) && alpha > 0)) stop("The initial step must be positive and finite.", call. = FALSE)
  cache$alpha <- alpha
  invisible(cache)
}

# One search: phi(a) is the merit at u + a du, dphi0 its slope at 0 (NULL
# when the search is derivative free), dphi(a) the slope at a; state holds
# what persists between iterations (non monotone history). Returns the step
# and a return code.
.ls_run <- function(ls, phi, dphi0, dphi = NULL, state = NULL, du_norm = 1, project = NULL) {
  name <- class(ls)[1L]
  phi0 <- phi(0)
  if (name == "NoLineSearch") return(LineSearchSolution(ls$alpha))
  if (name == "BackTracking") {
    a2 <- min(ls$initial_alpha, ls$maxstep / max(du_norm, .Machine$double.eps))
    f1 <- phi(a2); it <- 1L
    while (!is.finite(f1) && it <= 52L) { a2 <- a2 / 2; f1 <- phi(a2); it <- it + 1L }
    if (f1 <= phi0 + ls$c_1 * a2 * dphi0) return(LineSearchSolution(a2))
    at <- -(dphi0 * a2^2) / (2 * (f1 - phi0 - dphi0 * a2))
    a1 <- a2
    a2 <- max(min(at, a2 * ls$rho_hi), a2 * ls$rho_lo)
    f0 <- f1; f1 <- phi(a2)
    for (k in seq_len(max(0L, ls$maxiters - it))) {
      if (f1 <= phi0 + ls$c_1 * a2 * dphi0) return(LineSearchSolution(a2))
      at <- if (ls$order == 2L) -(dphi0 * a2^2) / (2 * (f1 - phi0 - dphi0 * a2)) else {
        d <- 1 / (a1^2 * a2^2 * (a2 - a1))
        aa <- (a1^2 * (f1 - phi0 - dphi0 * a2) - a2^2 * (f0 - phi0 - dphi0 * a1)) * d
        bb <- (-a1^3 * (f1 - phi0 - dphi0 * a2) + a2^3 * (f0 - phi0 - dphi0 * a1)) * d
        if (abs(aa) < .Machine$double.eps) dphi0 / (2 * bb) else {
          disc <- max(bb^2 - 3 * aa * dphi0, 0)
          (-bb + sqrt(disc)) / (3 * aa)
        }
      }
      a1 <- a2
      a2 <- max(min(at, a2 * ls$rho_hi), a2 * ls$rho_lo)
      if (!is.finite(a2)) a2 <- a1 * ls$rho_hi
      f0 <- f1; f1 <- phi(a2)
    }
    return(LineSearchSolution(a2, "Failure"))
  }
  if (name %in% c("ArmijoLineSearch", "ProjectedBackTracking")) {
    a <- state$alpha %||% ls$initial_alpha
    for (k in seq_len(ls$maxiters)) {
      slope <- a * dphi0
      v <- phi(a)
      if (is.finite(slope) && slope < 0 && is.finite(v) && v <= phi0 + ls$c_1 * slope) return(LineSearchSolution(a))
      a <- a * ls$contraction
    }
    return(LineSearchSolution(0, "Failure"))
  }
  if (name == "GoldenSection") {
    g <- (sqrt(5) - 1) / 2
    lo <- 0; hi <- 1
    x1 <- hi - g * (hi - lo); x2 <- lo + g * (hi - lo)
    f1 <- phi(x1); f2 <- phi(x2)
    for (k in seq_len(ls$maxiters)) {
      if (abs(hi - lo) < ls$tol) break
      if (f1 < f2) { hi <- x2; x2 <- x1; f2 <- f1; x1 <- hi - g * (hi - lo); f1 <- phi(x1) }
      else { lo <- x1; x1 <- x2; f1 <- f2; x2 <- lo + g * (hi - lo); f2 <- phi(x2) }
    }
    return(LineSearchSolution((lo + hi) / 2))
  }
  if (name == "LiFukushimaLineSearch") {
    nrm0 <- sqrt(2 * phi0)
    nrm <- function(a) sqrt(2 * phi(a))
    if (!is.finite(nrm0)) return(LineSearchSolution(1, "Failure"))
    if (nrm(1) <= ls$rho * nrm0 - ls$sigma_2 * du_norm^2) return(LineSearchSolution(1))
    l2 <- ls$lambda_0
    if (!is.finite(nrm(l2)) && !is.null(ls$nan_maxiters)) {
      ok <- FALSE
      for (k in seq_len(ls$nan_maxiters)) { l2 <- ls$beta * l2; if (is.finite(nrm(l2))) { ok <- TRUE; break } }
      if (!ok) return(LineSearchSolution(1, "Failure"))
    }
    for (k in seq_len(ls$maxiters)) {
      if (nrm(l2) <= (1 + ls$eta) * nrm0 - ls$sigma_1 * l2^2 * du_norm^2) return(LineSearchSolution(l2))
      l2 <- ls$beta * l2
    }
    return(LineSearchSolution(1, "Failure"))
  }
  if (name == "RobustNonMonotoneLineSearch") {
    fn <- function(a) (2 * phi(a))^(ls$n_exp / 2)
    f_old <- fn(0)
    ap <- am <- ls$sigma_1
    # nsteps is 0 on the first search, so eta is infinite and the first step is taken
    eta <- ls$eta_strategy(state$fn_1 %||% f_old, state$nsteps %||% 0L, NULL, f_old)
    fbar <- max(state$history %||% f_old)
    for (k in seq_len(ls$maxiters)) {
      f <- fn(ap)
      if (f <= fbar + eta - ls$gamma * ap * f_old) return(LineSearchSolution(ap))
      ap <- ap * min(max(ap * f_old / (f + (2 * ap - 1) * f_old), ls$tau_min), ls$tau_max)
      f <- fn(-am)
      if (f <= fbar + eta - ls$gamma * am * f_old) return(LineSearchSolution(-am))
      am <- am * min(max(am * f_old / (f + (2 * am - 1) * f_old), ls$tau_min), ls$tau_max)
    }
    return(LineSearchSolution(ls$sigma_1, "Failure"))
  }
  dphi <- dphi %||% function(a) { h <- 1e-7 * max(1, abs(a)); (phi(a + h) - phi(a - h)) / (2 * h) }
  if (name == "StrongWolfeLineSearch") return(.ls_strong_wolfe(ls, phi, phi0, dphi0, dphi))
  if (name == "MoreThuenteLineSearch") return(.ls_more_thuente(ls, phi, phi0, dphi0, dphi))
  if (name == "HagerZhangLineSearch") return(.ls_hager_zhang(ls, phi, phi0, dphi0, dphi))
  stop("unknown line search ", name, call. = FALSE)
}

# Strong Wolfe search with zoom (Nocedal and Wright, algorithms 3.5 and 3.6).
.ls_strong_wolfe <- function(ls, phi, phi0, dphi0, dphi) {
  zoom <- function(lo, hi, flo) {
    for (j in seq_len(ls$zoom_maxiters)) {
      a <- (lo + hi) / 2
      fa <- phi(a)
      if (fa > phi0 + ls$c1 * a * dphi0 || fa >= flo) hi <- a
      else {
        da <- dphi(a)
        if (abs(da) <= -ls$c2 * dphi0) return(a)
        if (da * (hi - lo) >= 0) hi <- lo
        lo <- a; flo <- fa
      }
    }
    lo
  }
  aprev <- 0; fprev <- phi0
  a <- ls$alpha_init
  for (i in seq_len(ls$maxiters)) {
    fa <- phi(a)
    if (!is.finite(fa)) { a <- (aprev + a) / 2; next }
    if (fa > phi0 + ls$c1 * a * dphi0 || (i > 1L && fa >= fprev)) return(LineSearchSolution(zoom(aprev, a, fprev)))
    da <- dphi(a)
    if (abs(da) <= -ls$c2 * dphi0) return(LineSearchSolution(a))
    if (da >= 0) return(LineSearchSolution(zoom(a, aprev, fa)))
    aprev <- a; fprev <- fa
    a <- min(2 * a, ls$alpha_max)
  }
  LineSearchSolution(a, "Failure")
}

# More and Thuente's search (dcsrch/dcstep of MINPACK-2).
.ls_more_thuente <- function(ls, phi, phi0, dphi0, dphi) {
  if (!(dphi0 < 0)) return(LineSearchSolution(0, "Failure"))
  ftol <- ls$ftol; gtol <- ls$gtol; xtol <- ls$xtol
  stpmin <- ls$alpha_min; stpmax <- ls$alpha_max
  stp <- ls$alpha_init
  ginit <- dphi0; gtest <- ftol * ginit
  brackt <- FALSE; stage <- 1L
  stx <- 0; fx <- phi0; gx <- ginit
  sty <- 0; fy <- phi0; gy <- ginit
  width <- stpmax - stpmin; width1 <- 2 * width
  stmin <- 0; stmax <- stp + 4 * stp
  cstep <- function(stx, fx, dx, sty, fy, dy, stp, fp, dp, brackt, stpmin, stpmax) {
    sgnd <- dp * sign(dx)
    if (fp > fx) {
      theta <- 3 * (fx - fp) / (stp - stx) + dx + dp
      s <- max(abs(c(theta, dx, dp)))
      gamma <- s * sqrt(max((theta / s)^2 - (dx / s) * (dp / s), 0)); if (stp < stx) gamma <- -gamma
      p <- (gamma - dx) + theta; q <- ((gamma - dx) + gamma) + dp
      stpc <- stx + (p / q) * (stp - stx)
      stpq <- stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2) * (stp - stx)
      stpf <- if (abs(stpc - stx) < abs(stpq - stx)) stpc else stpc + (stpq - stpc) / 2
      brackt <- TRUE
    } else if (sgnd < 0) {
      theta <- 3 * (fx - fp) / (stp - stx) + dx + dp
      s <- max(abs(c(theta, dx, dp)))
      gamma <- s * sqrt(max((theta / s)^2 - (dx / s) * (dp / s), 0)); if (stp > stx) gamma <- -gamma
      p <- (gamma - dp) + theta; q <- ((gamma - dp) + gamma) + dx
      stpc <- stp + (p / q) * (stx - stp)
      stpq <- stp + (dp / (dp - dx)) * (stx - stp)
      stpf <- if (abs(stpc - stp) > abs(stpq - stp)) stpc else stpq
      brackt <- TRUE
    } else if (abs(dp) < abs(dx)) {
      theta <- 3 * (fx - fp) / (stp - stx) + dx + dp
      s <- max(abs(c(theta, dx, dp)))
      gamma <- s * sqrt(max(0, (theta / s)^2 - (dx / s) * (dp / s))); if (stp > stx) gamma <- -gamma
      p <- (gamma - dp) + theta; q <- (gamma + (dx - dp)) + gamma
      r <- p / q
      stpc <- if (r < 0 && gamma != 0) stp + r * (stx - stp) else if (stp > stx) stpmax else stpmin
      stpq <- stp + (dp / (dp - dx)) * (stx - stp)
      if (brackt) {
        stpf <- if (abs(stpc - stp) < abs(stpq - stp)) stpc else stpq
        stpf <- if (stp > stx) min(stp + 0.66 * (sty - stp), stpf) else max(stp + 0.66 * (sty - stp), stpf)
      } else {
        stpf <- if (abs(stpc - stp) > abs(stpq - stp)) stpc else stpq
        stpf <- min(stpmax, max(stpmin, stpf))
      }
    } else {
      if (brackt) {
        theta <- 3 * (fp - fy) / (sty - stp) + dy + dp
        s <- max(abs(c(theta, dy, dp)))
        gamma <- s * sqrt(max((theta / s)^2 - (dy / s) * (dp / s), 0)); if (stp > sty) gamma <- -gamma
        p <- (gamma - dp) + theta; q <- ((gamma - dp) + gamma) + dy
        stpf <- stp + (p / q) * (sty - stp)
      } else stpf <- if (stp > stx) stpmax else stpmin
    }
    if (fp > fx) { sty <- stp; fy <- fp; dy <- dp }
    else {
      if (sgnd < 0) { sty <- stx; fy <- fx; dy <- dx }
      stx <- stp; fx <- fp; dx <- dp
    }
    list(stx = stx, fx = fx, dx = dx, sty = sty, fy = fy, dy = dy, stp = stpf, brackt = brackt)
  }
  for (it in seq_len(ls$maxiters)) {
    f <- phi(stp); g <- dphi(stp)
    ftest <- phi0 + stp * gtest
    if (stage == 1L && f <= ftest && g >= 0) stage <- 2L
    if (f <= ftest && abs(g) <= gtol * (-ginit)) return(LineSearchSolution(stp))
    if (brackt && (stp <= stmin || stp >= stmax)) return(LineSearchSolution(stp, "Failure"))
    if (brackt && stmax - stmin <= xtol * stmax) return(LineSearchSolution(stp))
    if (stp == stpmax && f <= ftest && g <= gtest) return(LineSearchSolution(stp))
    if (stp == stpmin && (f > ftest || g >= gtest)) return(LineSearchSolution(stp, "Failure"))
    if (stage == 1L && f <= fx && f > ftest) {
      fm <- f - stp * gtest; fxm <- fx - stx * gtest; fym <- fy - sty * gtest
      gm <- g - gtest; gxm <- gx - gtest; gym <- gy - gtest
      r <- cstep(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax)
      stx <- r$stx; sty <- r$sty; brackt <- r$brackt; stp <- r$stp
      fx <- r$fx + stx * gtest; fy <- r$fy + sty * gtest
      gx <- r$dx + gtest; gy <- r$dy + gtest
    } else {
      r <- cstep(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax)
      stx <- r$stx; fx <- r$fx; gx <- r$dx; sty <- r$sty; fy <- r$fy; gy <- r$dy; brackt <- r$brackt; stp <- r$stp
    }
    if (brackt) {
      if (abs(sty - stx) >= 0.66 * width1) stp <- stx + 0.5 * (sty - stx)
      width1 <- width; width <- abs(sty - stx)
      stmin <- min(stx, sty); stmax <- max(stx, sty)
    } else {
      stmin <- stp + 1.1 * (stp - stx); stmax <- stp + 4 * (stp - stx)
    }
    stp <- min(max(stp, stpmin), stpmax)
    if (brackt && (stp <= stmin || stp >= stmax || stmax - stmin <= xtol * stmax)) stp <- stx
  }
  LineSearchSolution(stp, "Failure")
}

# Hager and Zhang's search (CG_DESCENT): Wolfe or approximate Wolfe
# conditions with the secant2 bracketing steps.
.ls_hager_zhang <- function(ls, phi, phi0, dphi0, dphi) {
  d <- ls$delta; s <- ls$sigma; eps <- ls$epsilon * abs(phi0); th <- ls$theta
  wolfe <- function(a, fa, ga) (d * dphi0 >= (fa - phi0) / a && ga >= s * dphi0) ||
    ((2 * d - 1) * dphi0 >= ga && ga >= s * dphi0 && fa <= phi0 + eps)
  upd <- function(a, b, ga_b, c) {
    if (c <= a || c >= b) return(c(a, b))
    fc <- phi(c); gc <- dphi(c)
    if (gc >= 0) return(c(a, c))
    if (fc <= phi0 + eps) return(c(c, b))
    A <- a; B <- c
    for (k in 1:50) {
      dd <- (1 - th) * A + th * B
      fd <- phi(dd); gd <- dphi(dd)
      if (gd >= 0) return(c(A, dd))
      if (fd <= phi0 + eps) A <- dd else B <- dd
    }
    c(A, B)
  }
  secant <- function(a, b, ga, gb) (a * gb - b * ga) / (gb - ga)
  c0 <- ls$alpha_init
  # bracket
  a <- 0; ga <- dphi0; cc <- c0
  b <- NA
  for (k in seq_len(ls$maxiters)) {
    fc <- phi(cc); gc <- dphi(cc)
    if (!is.finite(fc)) { cc <- cc / 2; next }
    if (wolfe(cc, fc, gc)) return(LineSearchSolution(cc))
    if (gc >= 0) { b <- cc; break }
    if (fc > phi0 + eps) {
      ab <- upd(0, cc, gc, th * cc); a <- ab[1L]; b <- ab[2L]
      A <- 0; B <- cc
      for (j in 1:50) {
        dd <- (1 - th) * A + th * B
        fd <- phi(dd); gd <- dphi(dd)
        if (gd >= 0) { a <- A; b <- dd; break }
        if (fd <= phi0 + eps) A <- dd else B <- dd
      }
      break
    }
    a <- cc; cc <- min(ls$rho * cc, ls$alpha_max)
  }
  if (is.na(b)) return(LineSearchSolution(cc, "Failure"))
  for (k in seq_len(ls$maxiters)) {
    ga <- dphi(a); gb <- dphi(b)
    cs <- secant(a, b, ga, gb)
    ab <- upd(a, b, gb, cs)
    A <- ab[1L]; B <- ab[2L]
    if (is.finite(cs) && (cs == B || cs == A)) {
      gA <- dphi(A); gB <- dphi(B)
      c2 <- if (cs == B) secant(b, B, gb, gB) else secant(a, A, ga, gA)
      ab <- upd(A, B, gB, c2); A <- ab[1L]; B <- ab[2L]
    }
    for (x in c(A, B)) {
      fx <- phi(x)
      if (x > 0 && wolfe(x, fx, dphi(x))) return(LineSearchSolution(x))
    }
    if (B - A > 0.66 * (b - a)) {
      cc <- (A + B) / 2
      ab <- upd(A, B, dphi(B), cc); A <- ab[1L]; B <- ab[2L]
    }
    a <- A; b <- B
    if (b - a < 1e-14 * max(1, b)) break
  }
  LineSearchSolution(if (a > 0) a else b, "Failure")
}
