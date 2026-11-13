# Global error estimation, after GlobalDiffEq.jl v1.7.3: src/richardson.jl
# (GlobalRichardson), src/glee/tableaus.jl, algorithms.jl, solve.jl and
# perform_step.jl (the explicit general linear methods with global error
# estimation GLEE23, GLEE24, GLEE35 of Constantinescu 2016 and MM5GEE of
# Makazaga and Murua 2003; global_error_estimate) and Richardson.jl v1.4
# (extrapolate, the Neville tableau of Richardson extrapolation).
#
# A GLEE method integrates the partitioned state z = (y, eps): the solution
# y and an asymptotically correct estimate eps of its global error.
#   Y_i  = U1_i y + U2_i eps + dt sum_j A_ij f(Y_j)
#   y'   = y + dt sum_j B1_j f(Y_j),  eps' = eps + dt sum_j B2_j f(Y_j)
# As in Julia the problem is extended to the 2n state (f returns zero rates
# for eps), so the initial step, saving and dense output see the whole state,
# while the error estimate of the step controller is the increment of eps
# measured on y. The solution holds y; sol$global_error holds eps.
#
# GlobalErrorEstimation (the companion equation estimators) is in
# R/globaldiffeq-companion.R, GlobalAdjoint in R/globaldiffeq-adjoint.R.

# Tableaus in the y-eps form, the Float64 values of the exact rationals of
# GlobalDiffEq (tools/globaldiffeq-probe.jl prints them with repr).
.GLEE_TABS <- list(
  GLEE23 = list(
    order = 2L, stage1_fsal = TRUE, solution_stage = 0L,
    A = list(numeric(), c(1), c(0.25, 0.25)),
    U1 = c(1, 1, 1), U2 = c(0, 10, -1),
    B1 = c(0.08333333333333333, 0.08333333333333333, 0.8333333333333334),
    B2 = c(0.08333333333333333, 0.08333333333333333, -0.16666666666666666),
    c = c(0, 1, 0.5)),
  GLEE24 = list(
    order = 2L, stage1_fsal = FALSE, solution_stage = 0L,
    A = list(numeric(), c(0.75), c(0.25, 0.48333333333333334),
             c(-0.4772727272727273, 3.2954545454545454, -0x1.d1745d1745d17p+0)),
    U1 = c(1, 1, 1, 1), U2 = c(1, -0.29310344827586204, 1, 1),
    B1 = c(0x1.95e0595e0595ep-2, 0.7733333333333333, -0.33636363636363636, 0.16666666666666666),
    B2 = c(-0.12363636363636364, -0.7733333333333333, 0x1.304a7904a7905p+0, -0.2916666666666667),
    c = c(0, 0.75, 0.7333333333333333, 1)),
  GLEE35 = list(
    order = 3L, stage1_fsal = FALSE, solution_stage = 0L,
    A = list(numeric(), c(-0.08923467120428263),
             c(0.4943505136012235, -0.20930879618576095),
             c(0.2672542831101993, -0.5315983098317378, 1.0976653267020604),
             c(0.33695524969705265, -0.1092922590079333, -0.49456308711329705, 0.17355331181302344)),
    U1 = c(1, 1, 1, 1, 1),
    U2 = c(0.12420389705428309, -0.5227266959480462, 0x1.c1ea536248481p-4, 0.22674397647986055,
           0x1.d771358cb4293p-1),
    B1 = c(1.0800978502147018, -0.2696730454926488, 0.1515758333550661, 0.47080233376076, -0.4328029718378791),
    B2 = c(-1.1781678619971476, -0.2634261629453425, 0.4465694976326858, -0.06249847633365974,
           1.057523003643464),
    c = c(0, -0.08923467120428263, 0.2850417174154626, 0.8333212999805217, -0.09334678461115424)),
  MM5GEE = list(
    order = 5L, stage1_fsal = TRUE, solution_stage = 7L,
    A = list(numeric(), c(0.2), c(0.075, 0.225),
             c(0x1.f49f49f49f49fp-1, -3.7333333333333334, 3.5555555555555554),
             c(2.9525986892242035, -0x1.7310bd29520e4p+3, 9.822892851699436, -0.2908093278463649),
             c(2.8462752525252526, -10.757575757575758, 8.906422717743473, 0.2784090909090909,
               -0.2735313036020583),
             c(0x1.7555555555555p-4, 0, 0.44923629829290207, 0.6510416666666666, -0.322376179245283,
               0.13095238095238096),
             c(0.34865672575863554, -0.44326147133552585, 0.07403665054486801, 4.006711705879765,
               -5.6705632306057385, 4.724052718286655, -2.791759465478842),
             c(-1.0978863874778173, 3.11120755355011, -0x1.1d98b9d431204p+1, -10.433060281027627,
               15.144788740098694, -12.769884961024545, 7.535501566795239, 1.299437402369378),
             c(1.0811258559146544, -2.7774020191741435, 1.6732539147417829, 5.8141342443104955,
               -8.186456961426194, 6.903753869200552, -4.050168519594712, -0.07293520350314299,
               0.6146948195307068)),
    U1 = rep(1, 10),
    U2 = c(0, 0, 0, 0, 0, 0, 0, 0x1.f5f9f598c2fa8p-1, 1.0502232142857142, 0x1.f09c9f98a074ap-1),
    B1 = c(0x1.7555555555555p-4, 0, 0.44923629829290207, 0.6510416666666666, -0.322376179245283,
           0.13095238095238096, 0, 0, 0, 0),
    B2 = c(-0.019351583855225524, 0, -0.6188211010753616, -0.4476192487967594, 0.3724185995037892,
           -0x1.fce0c1bd42d9ap-3, -0.03927986906710311, 0.49173729184454695, 0.2965675298624644,
           0x1.b3dd0a84145d7p-3),
    c = c(0, 0.2, 0.3, 0.8, 0.8888888888888888, 1, 1, 0.24787363304981835, 0.5588803088803129,
          0.9999999999999988)))

# One step of GLEEConstantCache (the out of place perform_step!): stage sums
# accumulate term by term as Julia's unfused vector muladd; t + c dt is a
# fused scalar muladd.
.glee_step <- function(tab) {
  force(tab)
  function(ctx, t, u, dt, cache) {
    n <- length(u) %/% 2L
    iy <- seq_len(n)
    yprev <- u[iy]; eprev <- u[n + iy]
    f <- ctx$fobj$glee_inner
    p <- ctx$p
    s <- length(tab$c)
    ks <- vector("list", s)
    nf <- 0L
    if (tab$stage1_fsal) {
      ks[[1L]] <- cache$fsal[iy]
    } else {
      ks[[1L]] <- f(tab$U1[1L] * yprev + tab$U2[1L] * eprev, p, .fma(tab$c[1L], dt, t))
      nf <- nf + 1L
    }
    for (i in 2:s) {
      Y <- tab$U1[i] * yprev + tab$U2[i] * eprev
      a <- tab$A[[i]]
      for (j in seq_len(i - 1L)) if (a[j] != 0) Y <- Y + (dt * a[j]) * ks[[j]]
      ks[[i]] <- f(Y, p, .fma(tab$c[i], dt, t))
      nf <- nf + 1L
    }
    ynew <- yprev
    eloc <- numeric(n)
    for (j in seq_len(s)) {
      ynew <- ynew + (dt * tab$B1[j]) * ks[[j]]
      eloc <- eloc + (dt * tab$B2[j]) * ks[[j]]
    }
    fy <- if (tab$solution_stage == 0L) {
      nf <- nf + 1L
      f(ynew, p, t + dt)
    } else ks[[tab$solution_stage]]
    fsallast <- c(fy, eloc / dt)
    # calculate_residuals(eloc, yprev, ynew) on the y part
    EEst <- if (isTRUE(ctx$adaptive)) .err_norm(eloc, yprev, ynew, ctx$abstol, ctx$reltol)
    list(u = c(ynew, eprev + eloc), EEst = EEst, k = list(cache$fsal, fsallast),
         cache = list(fsal = fsallast), nf = nf)
  }
}

# __init of the GLEE methods: the problem is extended to the (y, eps) state
# and solved; then split into the solution and the global error estimate.
.glee_solve <- function(prob, alg, opts) {
  u0 <- prob$u0
  if (!is.numeric(u0) || !length(u0)) stop("GLEE methods require an array state", call. = FALSE)
  if (!.is_identity_mass(prob$f$mass_matrix)) stop("GLEE methods require the standard mass matrix", call. = FALSE)
  if (identical(prob$initialization_data$retcode, "InitialFailure")) return(.init_failure_solution(prob, alg))
  n <- length(u0)
  iy <- seq_len(n)
  p <- prob$p
  if (!is.null(prob$initialization_data %||% prob$f$initialization_data)) {
    # the initialization of the user's problem, then the extended problem
    iv <- get_initial_values(prob)
    if (isTRUE(iv$success)) {
      u0[] <- as.numeric(iv$u0)
      p <- iv$p
    }
  }
  f <- prob$f$f
  inner <- function(y, p, t) as.numeric(f(y, p, t))
  ext <- prob
  ext$f <- ODEFunction(function(z, p, t) c(inner(z[iy], p, t), numeric(n)))
  ext$f$glee_inner <- inner
  ext$initialization_data <- NULL
  un <- names(u0) %||% paste0("u", iy)
  ext$u0 <- stats::setNames(c(as.numeric(u0), numeric(n)), c(un, paste0("eps_", un)))
  ext$p <- p
  esol <- .integrate(ext, alg, opts)
  env <- .sol_env(esol)
  U <- env$U
  out <- .new_ode_solution(t = env$t, u = U[, iy, drop = FALSE], prob = prob, alg = alg,
                           retcode = env$retcode,
                           steps = if (length(env$steps)) lapply(env$steps, .glee_project_step, iy = iy),
                           stats = env$stats, p = env$p_final, events = env$events, dense = env$dense)
  oenv <- .sol_env(out)
  colnames(oenv$U) <- names(prob$u0)
  oenv$global_error <- lapply(seq_len(nrow(U)), function(i) unname(U[i, n + iy]))
  out
}

# The dense output of a step restricted to the solution partition (the
# interpolants are componentwise).
.glee_project_step <- function(step, iy) {
  step$y0 <- step$y0[iy]
  step$y1 <- step$y1[iy]
  step$k <- lapply(step$k, function(v) if (length(v) >= length(iy)) v[iy] else v)
  step
}

.glee_alg <- function(name) {
  tab <- .GLEE_TABS[[name]]
  .new_alg(name, tab$order, step = .glee_step(tab),
           init = function(ctx, t, u) list(fsal = ctx$f(u, ctx$p, t)),
           solve_override = .glee_solve, glee = TRUE)
}

#' Methods with global error estimation
#'
#' Explicit general linear methods that propagate the solution together
#' with an asymptotically correct estimate of its global (accumulated)
#' error, as GlobalDiffEq.jl. `GLEE23()` (3 stages, order 2), `GLEE24()`
#' (4 stages, order 2, both decoupling conditions; the recommended second
#' order method) and `GLEE35()` (5 stages, order 3) are the methods of
#' Constantinescu (2016), eqs. (4.6), (A.3) and (4.9) (PETSc's TSGLEE23,
#' TSGLEE24 and TSGLEE35); `MM5GEE()` is the scheme of Makazaga and Murua
#' (2003) on the Dormand-Prince 5(4) pair, which propagates an order 6
#' companion solution with three extra stages.
#'
#' Solving an `ODEProblem` with these methods returns an ordinary solution
#' whose `sol$global_error` holds the estimated global error at each saved
#' time (also [global_error_estimate()]). The per step increment of the
#' estimate is an asymptotically correct local error estimate and drives the
#' usual step size control, so the tolerances behave as for adaptive Runge
#' Kutta methods. Only explicit ODEs without a mass matrix are supported.
#'
#' @return An ODE algorithm object.
#' @references E. M. Constantinescu, Estimating global errors in time
#'   stepping, SIAM Journal on Numerical Analysis 54(6), 2016.
#'   J. Makazaga and A. Murua, New Runge-Kutta based schemes for ODEs with
#'   cheap global error estimation, BIT Numerical Mathematics 43, 2003.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
#' sol <- solve(prob, GLEE24(), reltol = 1e-6, abstol = 1e-8)
#' exp(-1) - sol$u[[length(sol$t)]]   # the true global error
#' global_error_estimate(sol, length(sol$t))
#' @name glee_algorithms
NULL

#' @rdname glee_algorithms
#' @export
GLEE23 <- function() .glee_alg("GLEE23")
#' @rdname glee_algorithms
#' @export
GLEE24 <- function() .glee_alg("GLEE24")
#' @rdname glee_algorithms
#' @export
GLEE35 <- function() .glee_alg("GLEE35")
#' @rdname glee_algorithms
#' @export
MM5GEE <- function() .glee_alg("MM5GEE")

#' Global error estimate of a solution
#'
#' `global_error_estimate(sol)` returns the estimated global errors of a
#' solution computed with a GLEE method ([GLEE23()], [GLEE24()],
#' [GLEE35()], [MM5GEE()]): a list matching `sol$u`, whose element `i` is
#' the estimated global error of `sol$u[[i]]` at `sol$t[i]`;
#' `global_error_estimate(sol, i)` returns its `i`th element. For a
#' [GlobalAdjoint()] solution it returns `sol$global_error`, the estimated
#' error 2-norms.
#'
#' `global_error_estimate(prob, alg)` with a [GlobalErrorEstimation()]
#' algorithm solves `prob` and returns the estimated 2-norm of its global
#' error at the final time.
#'
#' @param x A solution of a GLEE method, or an `ODEProblem`.
#' @param i Index of a saved time.
#' @param alg A [GlobalErrorEstimation()] algorithm.
#' @param ... Unused for solutions; options of the forward solve for
#'   problems.
#' @param abstol,reltol Tolerances of the forward solve.
#' @param companion_abstol,companion_reltol Tolerances of the companion
#'   solve (by default those of `alg`).
#' @return A list of vectors, one vector, or the estimated error norm.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
#' sol <- solve(prob, GLEE35(), reltol = 1e-6)
#' global_error_estimate(sol, 2)
#' global_error_estimate(prob, GlobalErrorEstimation(Tsit5()), abstol = 1e-5, reltol = 1e-5)
#' @export
global_error_estimate <- function(x, ...) UseMethod("global_error_estimate")

#' @rdname global_error_estimate
#' @export
global_error_estimate.default <- function(x, i = NULL, ...) {
  ge <- if (inherits(x, "SciMLSolution")) .sol_env(x)$global_error
  if (is.null(ge)) stop("global_error_estimate expects a solution produced by a GLEE method", call. = FALSE)
  if (is.null(i)) ge else ge[[i]]
}

# ---------------------------------------------------------------------------
# Richardson extrapolation

# Richardson.jl _extrapolate: the Neville tableau along h, h c, h c^2, ...
# (f values are numeric arrays; norm is the 2-norm of all entries).
.richardson_extrapolate <- function(f, h, contract = 0.125, x0 = 0, power = 1, atol = 0,
                                    rtol = if (atol > 0) 0 else sqrt(.Machine$double.eps),
                                    maxeval = Inf, breaktol = 2) {
  if (!(rtol >= 0 && atol >= 0)) stop("rtol and atol must be nonnegative", call. = FALSE)
  if (!(breaktol > 0)) stop("breaktol must be positive", call. = FALSE)
  if (!(abs(contract) > 0 && abs(contract) < 1)) stop("contract must be in (0,1)", call. = FALSE)
  nrm <- function(x) sqrt(sum(x * x))
  invcontract <- (1 / contract)^power
  neville <- list(f(x0 + h))
  f0 <- neville[[1L]]
  err <- Inf
  numeval <- 1
  while (numeval < maxeval) {
    numeval <- numeval + 1
    h <- h * contract
    neville[[length(neville) + 1L]] <- f(x0 + h)
    cc <- invcontract
    minerr <- Inf
    for (i in rev(seq_len(length(neville) - 1L))) {
      old <- neville[[i]]
      neville[[i]] <- neville[[i + 1L]] + (neville[[i + 1L]] - neville[[i]]) / (cc - 1)
      e <- nrm(neville[[i]] - old)
      minerr <- min(minerr, e)
      if (e < err) {
        f0 <- neville[[i]]
        err <- e
      }
      cc <- cc * invcontract
    }
    if (minerr > breaktol * err || !is.finite(minerr)) break
    if (err <= max(rtol * nrm(f0), atol)) break
  }
  list(value = f0, err = err)
}

#' Global Richardson extrapolation
#'
#' `GlobalRichardson(alg)` wraps an ODE (or DAE) algorithm with global
#' Richardson extrapolation, as GlobalDiffEq.jl: the problem is solved with
#' fixed steps `dt`, `dt / 2`, `dt / 4`, ... and the solutions at the times
#' `tstops` (default `seq(t0, tf, by = dt)`) are extrapolated to zero step
#' size with a Neville tableau until the estimated error is below
#' `max(reltol * norm, abstol)` (or stops decreasing). As in Julia the
#' returned solution is the last (finest) of these fixed step solves; the
#' extrapolated values at `tstops` and their error estimate are attached as
#' `sol$extrapolated` (one column per time) and `sol$extrapolation_error`.
#'
#' @param alg The algorithm of the fixed step solves.
#' @return An algorithm object; `solve()` needs `dt`.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
#' sol <- solve(prob, GlobalRichardson(Tsit5()), dt = 0.25, reltol = 1e-8)
#' length(sol$t)
#' sol$extrapolated[, 5] - exp(-1)
#' @export
GlobalRichardson <- function(alg) {
  if (!inherits(alg, "ODEAlgorithm")) stop("GlobalRichardson wraps an ODE algorithm", call. = FALSE)
  structure(list(name = "GlobalRichardson", alg = alg, order = alg$order, adaptive = FALSE,
                 dae = isTRUE(alg$dae), mass_matrix = isTRUE(alg$mass_matrix),
                 solve_override = .global_richardson_solve),
            class = c("GlobalRichardson", "ODEAlgorithm"))
}

#' @export
print.GlobalRichardson <- function(x, ...) {
  cat("GlobalRichardson(", x$alg$name, "())\n", sep = "")
  invisible(x)
}

.global_richardson_solve <- function(prob, alg, opts) {
  dt <- opts$dt
  if (is.null(dt)) stop("GlobalRichardson needs `dt`", call. = FALSE)
  other <- opts
  other$dt <- NULL
  t0 <- prob$tspan[1L]; tf <- prob$tspan[2L]
  # range(t0, stop = tf, step = dt), the Julia float range
  tstops <- as.numeric(opts$tstops %||% .sids_float_range(t0, dt, tf))
  last <- NULL
  inner <- alg$alg
  ex <- .richardson_extrapolate(function(h) {
    last <<- do.call(solve, c(list(prob, inner), other, list(dt = h, adaptive = FALSE)))
    t(.interp_states(.sol_env(last), tstops))
  }, dt, contract = 0.5, rtol = opts$reltol %||% 1e-3, atol = opts$abstol %||% 1e-6)
  env <- .sol_env(last)
  env$extrapolated <- ex$value
  env$extrapolation_error <- ex$err
  last
}
