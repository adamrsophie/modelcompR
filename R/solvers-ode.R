# Pure R ODE integrators, after OrdinaryDiffEq.jl (OrdinaryDiffEqCore
# solve.jl, integrators/integrator_utils.jl, integrators/controllers.jl,
# initdt.jl; OrdinaryDiffEqTsit5; OrdinaryDiffEqDefault/src/default_alg.jl).
#
# An algorithm is an S3 object of class c("<Name>", "ODEAlgorithm") with
#   order      the order used by the step size controller
#   adaptive   whether it estimates local errors
#   step(ctx, t, u, dt, cache)  one step from (t, u); returns
#              list(u = <new state>, err = <local error vector or NULL>,
#                   k = <interpolation data>, cache = <updated cache>,
#                   fail = <TRUE when an implicit solve failed>)
#              and optionally EEst (scalar error norm, replaces err),
#              eigen_est (stiffness estimate for automatic switching),
#              du (DAE problems) and the counts nf, njacs, nsolve, nw.
#              A cache that is an environment is kept after rejected steps
#              (implicit methods keep their Jacobians and failure counters).
#   interp(theta, dt, y0, y1, k)  dense output inside a step
#   init(ctx, t, u)  initial cache (for example the FSAL derivative)
# Optional fields set the controller ("PI", "I" or "BDF" with the hooks
# bdf_stepsize/bdf_accept/bdf_reject), its constants (beta1, beta2, qmin,
# qmax, gamma, qsteady_min, qsteady_max, failfactor), newton_fail and
# on_reject hooks, mass_matrix (supports singular mass matrices) and dae
# (solves DAEProblem). ctx carries f, p, mass, jac, tgrad, tolerances and
# the integrator state the implicit methods need (iter, success_iter,
# accept_prev, EEst_prev, tprev, uprev2, du).
# The driver .integrate() handles step size control, saving, tstops and
# callbacks, and returns an ODESolution.

# mass_matrix follows ... so that an option such as the Krylov dimension m
# of the exponential methods is not partially matched to it
.new_alg <- function(name, order, adaptive = TRUE, step, interp = NULL, init = NULL,
                     stiff = FALSE, ..., mass_matrix = FALSE) {
  structure(list(name = name, order = order, adaptive = adaptive, step = step,
                 interp = interp %||% .hermite_interp, init = init %||% .fsal_init,
                 stiff = stiff, mass_matrix = mass_matrix, ...),
            class = c(name, "ODEAlgorithm"))
}

#' @export
print.ODEAlgorithm <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

.fsal_init <- function(ctx, t, u) list(fsal = ctx$f(u, ctx$p, t))

# Cubic Hermite interpolation from the derivatives at both ends of a step;
# k = list(f0, f1).
.hermite_interp <- function(theta, dt, y0, y1, k) {
  f0 <- k[[1L]]; f1 <- k[[2L]]
  (1 - theta) * y0 + theta * y1 +
    theta * (theta - 1) * ((1 - 2 * theta) * (y1 - y0) + (theta - 1) * dt * f0 + theta * dt * f1)
}

# The in place form (hermite_interpolant!, a @.. broadcast under @muladd):
# every multiply add is fused element by element. OrdinaryDiffEq evaluates
# the root function of discontinuity detection in place.
.hermite_interp_fused <- function(theta, dt, y0, y1, k) {
  s <- (1 - 2 * theta) * (y1 - y0)
  s <- .fma((theta - 1) * dt, k[[1L]], s)
  s <- .fma(theta * dt, k[[2L]], s)
  .fma(theta * (theta - 1), s, .fma(theta, y1, (1 - theta) * y0))
}

# ---------------------------------------------------------------------------
# Explicit methods

#' ODE solver algorithms
#'
#' Algorithms for [solve()] on ODE problems, named as in OrdinaryDiffEq.jl.
#' `Tsit5()` is the Tsitouras 5(4) Runge-Kutta method with its free 4th
#' order interpolant, the default for nonstiff problems. `Euler()` and
#' `RK4()` are fixed step methods (pass `dt`). `BS3()` is the
#' Bogacki-Shampine 3(2) method, `DP5()` the Dormand-Prince 5(4) method
#' with its free dense output, and `Vern7()` and `Vern9()` are Verner's
#' efficient 7(6) and 9(8) methods with their high order interpolants
#' (whose extra stages are computed when the solution is first
#' interpolated inside a step).
#'
#' @return An algorithm object.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
#' sol <- solve(prob, Vern7(), reltol = 1e-10, abstol = 1e-12)
#' sol(0.5) - exp(-0.5)
#' solve(prob, RK4(), dt = 0.1)$retcode
#' @name ode_algorithms
NULL

.TSIT5 <- list(
  c = c(0.161, 0.327, 0.9, 0x1.f5c5e85edbf1dp-1, 1, 1),
  a = list(c(0.161),
           c(-0.008480655492356989, 0.335480655492357),
           c(2.8971530571054935, -6.359448489975075, 4.3622954328695815),
           c(5.325864828439257, -11.748883564062828, 7.4955393428898365, -0x1.7adc1b510a9a4p-4),
           c(5.86145544294642, -12.92096931784711, 8.159367898576159, -0.071584973281401, -0.028269050394068383),
           c(0x1.8b1a71ed4a6eep-4, 0.01, 0.4798896504144996, 1.379008574103742, -3.290069515436081, 2.324710524099774)),
  btilde = c(-0x1.d29e8505f7d7fp-10, -0.0008164344596567469, 0.007880878010261995,
             -0.1447110071732629, 0.5823571654525552, -0.45808210592918697, 0.015151515151515152),
  r = list(c(1.0, -2.763706197274826, 2.9132554618219126, -1.0530884977290216),
           c(0.13169999999999998, -0.2234, 0.1017),
           c(0x1.f713f27487049p+1, -5.941033872131505, 2.490627285651253),
           c(-12.411077166933676, 30.33818863028232, -16.548102889244902),
           c(37.50931341651104, -88.1789048947664, 47.37952196281928),
           c(-27.896526289197286, 65.09189467479366, -34.87065786149661),
           c(1.5, -4, 2.5)))

.tsit5_step <- function(ctx, t, u, dt, cache) {
  f <- ctx$f; p <- ctx$p
  A <- .TSIT5$a; cc <- .TSIT5$c
  k1 <- cache$fsal
  k2 <- f(u + dt * (A[[1]][1] * k1), p, t + cc[1] * dt)
  k3 <- f(u + dt * (A[[2]][1] * k1 + A[[2]][2] * k2), p, t + cc[2] * dt)
  k4 <- f(u + dt * (A[[3]][1] * k1 + A[[3]][2] * k2 + A[[3]][3] * k3), p, t + cc[3] * dt)
  k5 <- f(u + dt * (A[[4]][1] * k1 + A[[4]][2] * k2 + A[[4]][3] * k3 + A[[4]][4] * k4), p, t + cc[4] * dt)
  g6 <- u + dt * (A[[5]][1] * k1 + A[[5]][2] * k2 + A[[5]][3] * k3 + A[[5]][4] * k4 + A[[5]][5] * k5)
  k6 <- f(g6, p, t + dt)
  a7 <- A[[6]]
  unew <- u + dt * (a7[1] * k1 + a7[2] * k2 + a7[3] * k3 + a7[4] * k4 + a7[5] * k5 + a7[6] * k6)
  k7 <- f(unew, p, t + dt)
  bt <- .TSIT5$btilde
  err <- dt * (bt[1] * k1 + bt[2] * k2 + bt[3] * k3 + bt[4] * k4 + bt[5] * k5 + bt[6] * k6 + bt[7] * k7)
  # Hairer II, page 22, with the infinity norm (stiffness detection)
  eig <- if (isTRUE(ctx$composite)) .eigen_est(k7 - k6, unew - g6) else NULL
  list(u = unew, err = err, k = list(k1, k2, k3, k4, k5, k6, k7), cache = list(fsal = k7), nf = 6L,
       eigen_est = eig)
}

# norm((k7 - k6) ./ (g7 - g6), Inf). A non finite estimate counts as stiff,
# as in Julia, except that components with 0/0 (states whose stages
# coincide, such as constant states) are skipped: Julia's NaN there would
# flag every such problem as stiff and make the default method alternate.
.eigen_est <- function(df, du) {
  if (!length(df)) return(0)
  r <- abs(df / du)
  r[df == 0 & du == 0] <- 0
  max(r)
}

.tsit5_interp <- function(theta, dt, y0, y1, k) {
  r <- .TSIT5$r
  th2 <- theta * theta
  b1 <- theta * (r[[1]][1] + theta * (r[[1]][2] + theta * (r[[1]][3] + theta * r[[1]][4])))
  bs <- vapply(2:7, function(i) th2 * (r[[i]][1] + theta * (r[[i]][2] + theta * r[[i]][3])), 0)
  y0 + dt * (k[[1]] * b1 + k[[2]] * bs[1] + k[[3]] * bs[2] + k[[4]] * bs[3] + k[[5]] * bs[4] +
               k[[6]] * bs[5] + k[[7]] * bs[6])
}

#' @rdname ode_algorithms
#' @export
Tsit5 <- function() .new_alg("Tsit5", 5L, step = .tsit5_step, interp = .tsit5_interp,
                              stability_size = 3.5068)

#' @rdname ode_algorithms
#' @export
Euler <- function() {
  .new_alg("Euler", 1L, adaptive = FALSE,
    step = function(ctx, t, u, dt, cache) {
      k1 <- cache$fsal
      unew <- u + dt * k1
      k2 <- ctx$f(unew, ctx$p, t + dt)
      list(u = unew, err = NULL, k = list(k1, k2), cache = list(fsal = k2), nf = 1L)
    })
}

#' @rdname ode_algorithms
#' @export
RK4 <- function() {
  .new_alg("RK4", 4L, adaptive = FALSE,
    step = function(ctx, t, u, dt, cache) {
      f <- ctx$f; p <- ctx$p
      k1 <- cache$fsal
      k2 <- f(u + dt / 2 * k1, p, t + dt / 2)
      k3 <- f(u + dt / 2 * k2, p, t + dt / 2)
      k4 <- f(u + dt * k3, p, t + dt)
      unew <- u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      k5 <- f(unew, p, t + dt)
      list(u = unew, err = NULL, k = list(k1, k5), cache = list(fsal = k5), nf = 4L)
    })
}

# ---------------------------------------------------------------------------
# Step size control and initial step, after OrdinaryDiffEqCore
# (initdt.jl, integrators/controllers.jl, alg_utils.jl).

# ODE_DEFAULT_NORM: the root of the sum of squares over the length (mean()
# adds a correction pass that rounds differently).
.rms_norm <- function(x) sqrt(sum(x * x) / length(x))

.err_norm <- function(err, u0, u1, abstol, reltol) {
  if (!length(err)) return(0)
  .rms_norm(err / (abstol + pmax(abs(u0), abs(u1)) * reltol))
}

# _ode_initdt_oop of initdt.jl (Hairer and Wanner's estimate).
.initial_dt <- function(ctx, t0, u0, f0, tdir, order, abstol, reltol, dtmax, dtmin = 0) {
  if (!length(u0)) return(dtmax)
  smalldt <- max(dtmin, 1e-6)
  sc <- abstol + abs(u0) * reltol
  d0 <- .rms_norm(u0 / sc)
  d1 <- .rms_norm(f0 / sc)
  if (!is.finite(d1)) return(max(dtmin, 1e-12))
  dt0 <- if (d0 < 1e-5 || d1 < 1e-5) smalldt else (d0 / d1) / 100
  dt0 <- min(dt0, dtmax)
  u1 <- u0 + tdir * dt0 * f0
  f1 <- ctx$f(u1, ctx$p, t0 + tdir * dt0)
  if (identical(f0, f1)) return(max(dtmin, 100 * dt0))
  d2 <- .rms_norm((f1 - f0) / sc) / dt0
  if (!is.finite(d2)) d2 <- 1e10
  dt1 <- if (max(d1, d2) <= 1e-15) max(smalldt, dt0 * 1e-3) else 10^(-(2 + log10(max(d1, d2))) / order)
  max(dtmin, min(100 * dt0, dt1, dtmax))
}

# Controller settings of an algorithm: the *_default functions of
# OrdinaryDiffEqCore/src/alg_utils.jl with the per algorithm overrides the
# constructors store in the algorithm object. `type` is "PI" (PIController),
# "I" (IController), "Predictive" (PredictiveController) or "BDF"
# (BDFController, OrdinaryDiffEqBDF/src/controllers.jl).
.ctrl_opts <- function(alg, opts) {
  ord <- alg$order
  beta2 <- opts$beta2 %||% alg$beta2 %||% (2 / (5 * ord))
  list(type = alg$controller %||% "PI",
       pid_beta = opts$pid_beta %||% alg$pid_beta, accept_safety = alg$accept_safety %||% 0.81,
       pid_k = min(ord, alg$adaptive_order %||% ord) + 1,
       adaptive_order = alg$adaptive_order %||% (ord - 1),
       beta1 = opts$beta1 %||% alg$beta1 %||% (7 / (10 * ord)), beta2 = beta2,
       qmin = opts$qmin %||% alg$qmin %||% (1 / 5), qmax = opts$qmax %||% alg$qmax %||% 10,
       qmax_first = 10000, gamma = opts$gamma %||% alg$gamma %||% 0.9,
       qsteady_min = alg$qsteady_min %||% 1, qsteady_max = alg$qsteady_max %||% 1,
       failfactor = opts$failfactor %||% alg$failfactor %||% 2,
       # handle_disco_reject! / handle_disco_accept! (R/solvers-disco.R) of the
       # I, PI, PID and predictive controllers
       discontinuity_detection = isTRUE(opts$discontinuity_detection) &&
         (alg$controller %||% "PI") %in% c("I", "PI", "PID", "Predictive"))
}

# FastPower.fastpower, which the I and PI controllers of OrdinaryDiffEqCore
# use for EEst^expo: a Float32 approximation of x^y through a fast log2 (David
# Goldberg's rational formula) and Julia's Float32 exp2 kernel. Float32 values
# are emulated by rounding doubles with writeBin(size = 4); the products and
# sums of two Float32 values are exact in double precision, so every fused
# multiply add rounds once, as in Julia. Matching it bit for bit keeps the
# step sizes equal to Julia's.
.f32 <- function(x) readBin(writeBin(as.numeric(x), raw(), size = 4L), "double", n = length(x), size = 4L)

.f32_bits <- function(x) readBin(writeBin(as.numeric(x), raw(), size = 4L), "integer", n = length(x), size = 4L)

.f32_from_bits <- function(b) readBin(writeBin(as.integer(b), raw(), size = 4L), "double", n = length(b), size = 4L)

.fastlog2 <- function(x) {
  u <- .f32_bits(x)
  hi <- bitwAnd(u, 0x00400000L) == 0L
  e <- bitwShiftR(bitwAnd(u, 0x7F800000L), 23L) - ifelse(hi, 127L, 126L)
  sig <- .f32_from_bits(bitwOr(bitwAnd(u, 0x007FFFFFL), ifelse(hi, 0x3f800000L, 0x3f000000L)))
  quot <- .f32(.f32(sig * .f32(0.338953) + .f32(1.8596461)) / .f32(sig + .f32(0.523692)))
  .f32((sig - 1) * quot + e)
}

.exp2_f32 <- function(x) {
  cf <- .f32(c(1, 0.6931472, 0.2402265, 0.05550411, 0.009618025, 0.0013333423, 0.00015469732, 1.5316464e-5))
  N <- round(x)
  r <- x - N
  s <- rep(cf[8L], length(x))
  for (i in 7:1) s <- .f32(r * s + cf[i])
  out <- .f32(2^N * s)
  out[x >= 128] <- Inf
  out[x <= .f32(-126.00001)] <- 0
  out
}

.fastpower <- function(x, y) {
  out <- .exp2_f32(.f32(.f32(y) * .fastlog2(x)))
  out[x == 0] <- 0
  out[is.infinite(x) & is.infinite(y)] <- Inf
  out
}

.ctrl_state <- function() list(errold = 1e-4, q11 = 1, dtreject = NA_real_, err = c(1, 1, 1), dt_factor = 1,
                               dtacc = 1, erracc = 1, qold = 1, first = TRUE)

# PIDController (integrators/controllers.jl): the step factor from the last
# three error estimates, limited by 1 + atan(x - 1); a step is accepted when
# the factor is at least accept_safety.
.ctrl_pid_q <- function(co, cs, EEst) {
  EEst <- if (is.finite(EEst)) max(EEst, .Machine$double.eps) else Inf
  cs$err[1L] <- 1 / EEst
  b <- co$pid_beta; k <- co$pid_k
  f <- cs$err[1L]^(b[1L] / k) * cs$err[2L]^(b[2L] / k) * cs$err[3L]^(b[3L] / k)
  if (is.nan(f)) f <- 0
  cs$dt_factor <- 1 + atan(f - 1)
  list(q = 1 / cs$dt_factor, cs = cs)
}

# Whether the controller rejects the step.
.ctrl_rejects <- function(co, cs, EEst) {
  if (identical(co$type, "PID")) return(cs$dt_factor < co$accept_safety)
  EEst > 1
}

# stepsize_controller!: the factor q with dt_new = dt / q.
.ctrl_q <- function(co, cs, EEst, first, dt, alg, cache) {
  qmax <- if (first) co$qmax_first else co$qmax
  if (co$type == "BDF") return(list(q = alg$bdf_stepsize(cache, EEst, co, first), cs = cs))
  if (co$type == "PID") return(.ctrl_pid_q(co, cs, EEst))
  if (co$type == "Predictive") {
    # PredictiveController: the factor from the error of this step alone
    cs$first <- first
    if (EEst == 0) return(list(q = 1 / qmax, cs = cs))
    q <- if (is.finite(EEst)) .fastpower(EEst, 1 / (co$adaptive_order + 1)) / co$gamma else Inf
    q <- min(1 / co$qmin, max(1 / qmax, q))
    cs$qold <- q
    return(list(q = q, cs = cs))
  }
  if (!is.finite(EEst)) {
    cs$q11 <- Inf
    cs$dtreject <- dt * co$qmin
    return(list(q = 1 / co$qmin, cs = cs))
  }
  if (EEst == 0) return(list(q = 1 / qmax, cs = cs))
  if (co$type == "I") {
    # get_current_adaptive_order: alg_adaptive_order (alg_order - 1 unless
    # the method sets it), the current order of a variable order method
    ord <- if (is.environment(cache) && !is.null(cache$order)) cache$order else co$adaptive_order
    q <- min(1 / co$qmin, max(1 / qmax, .fastpower(EEst, 1 / (ord + 1)) / co$gamma))
    cs$dtreject <- dt / q
    return(list(q = q, cs = cs))
  }
  q11 <- .fastpower(EEst, co$beta1)
  q <- q11 / .fastpower(cs$errold, co$beta2)
  cs$q11 <- q11
  list(q = min(1 / co$qmin, max(1 / qmax, q / co$gamma)), cs = cs)
}

# step_accept_controller!: the next step size.
.ctrl_accept <- function(co, cs, q, EEst, dt, alg, cache) {
  if (co$type == "BDF") return(list(dt = alg$bdf_accept(cache, q, dt, co), cs = cs))
  if (co$type == "PID") {
    fac <- cs$dt_factor
    if (1 / fac >= co$qsteady_min && 1 / fac <= co$qsteady_max) fac <- 1
    cs$err[3L] <- cs$err[2L]; cs$err[2L] <- cs$err[1L]
    return(list(dt = dt * fac, cs = cs))
  }
  if (co$type == "Predictive") {
    # Gustafsson's predictive factor from the last accepted step
    qacc <- q
    if (!cs$first) {
      qgus <- (cs$dtacc / dt) * .fastpower((EEst^2) / cs$erracc, 1 / (co$adaptive_order + 1))
      qgus <- max(1 / co$qmax, min(1 / co$qmin, qgus / co$gamma))
      qacc <- max(q, qgus)
    }
    if (qacc >= co$qsteady_min && qacc <= co$qsteady_max) qacc <- 1
    cs$dtacc <- dt
    cs$erracc <- max(1e-2, EEst)
    return(list(dt = dt / qacc, cs = cs))
  }
  if (q >= co$qsteady_min && q <= co$qsteady_max) q <- 1
  if (co$type == "PI") cs$errold <- max(EEst, 1e-4)
  list(dt = dt / q, cs = cs)
}

# step_reject_controller!: the step size to retry with.
.ctrl_reject <- function(co, cs, dt, EEst, alg, cache) {
  if (co$type == "BDF") return(alg$bdf_reject(cache, EEst, dt, co))
  if (co$type == "PID") return(dt * cs$dt_factor)
  if (co$type == "I") return(cs$dtreject)
  if (co$type == "Predictive") return(if (cs$first) 0.1 * dt else dt / cs$qold)
  dt / min(1 / co$qmin, cs$q11 / co$gamma)
}

# ---------------------------------------------------------------------------
# The integration driver

# The default algorithm, OrdinaryDiffEqDefault/src/default_alg.jl: a
# composite of nonstiff and stiff methods with automatic switching
# (DefaultODEAlgorithm, R/solvers-default.R).
.default_ode_alg <- function(prob) DefaultODEAlgorithm()

.is_identity_mass <- function(M) {
  is.null(M) || !is.matrix(M) || (nrow(M) == ncol(M) && all(M == diag(nrow(M))))
}

#' Solve a problem
#'
#' `solve(prob, alg)` solves a problem built by [ODEProblem()] and the other
#' problem constructors. It is a method of the base R generic [base::solve()],
#' so `solve()` keeps working for matrices and for symcompR expressions.
#'
#' @param a A problem.
#' @param b An algorithm such as [Tsit5()]; [DefaultODEAlgorithm()] when
#'   missing.
#' @param ... Solver options: `saveat` (times to save at), `abstol`,
#'   `reltol`, `dt` (initial or fixed step), `dtmax`, `dtmin`, `maxiters`,
#'   `tstops`, `callback`, `save_everystep`, `save_start`, `save_end`,
#'   `dense`, `adaptive`, `step_limiter` (a function `limiter(u, integ, p,
#'   t)` returning the state to keep after each accepted step), and the
#'   controller settings `qmin`, `qmax`, `gamma`, `beta1`, `beta2`.
#' @return A solution object. Its `retcode` is one of "Success",
#'   "Terminated", "MaxIters", "DtLessThanMin", "Unstable",
#'   "InitialFailure" or "ConvergenceFailure"; `sol$stats` counts function
#'   evaluations (`nf`), Jacobians (`njacs`), linear solves (`nsolve`),
#'   factorizations (`nw`), accepted (`naccept`) and rejected (`nreject`)
#'   steps.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -p * u, 1, c(0, 1), p = 2)
#' sol <- solve(prob, Tsit5(), reltol = 1e-8)
#' sol(1) - exp(-2)
#' @export
solve.ODEProblem <- function(a, b, ...) {
  alg <- if (missing(b) || is.null(b)) .default_ode_alg(a) else b
  if (!inherits(alg, "ODEAlgorithm")) stop("not an ODE algorithm: ", class(alg)[1L], call. = FALSE)
  if (isTRUE(alg$dae)) stop(alg$name, "() solves DAEProblem objects; use a mass matrix method such as Rodas5P()",
                            call. = FALSE)
  if (!.is_identity_mass(a$f$mass_matrix) && !isTRUE(alg$mass_matrix))
    stop(alg$name, "() does not support mass matrices; use a mass matrix solver such as Rodas5P()",
         call. = FALSE)
  # a problem whose initialization failed is not integrated (R/initialization-problem.R)
  if (identical(a$initialization_data$retcode, "InitialFailure")) return(.init_failure_solution(a, alg))
  opts <- utils::modifyList(a$kwargs %||% list(), list(...))
  # solve(prob, alg, sensealg = ...) records events for the adjoint (R/sensitivity-interface.R)
  if (!is.null(opts$sensealg) && !is.null(.mtk$hooks$sensitivity_solve))
    return(.mtk$hooks$sensitivity_solve(a, alg, opts))
  .integrate(a, alg, opts)
}

.integrate <- function(prob, alg, opts) {
  I <- .integ_init(prob, alg, opts)
  while (!I$done) .integ_step(I)
  .integ_finalize(I)
}

# The integrator of .integrate in three parts, as SciMLBase's init, step!
# and the end of solve!: the state of the integration lives in the
# environment I, and the code runs in it (evalq), so its variables persist
# between steps. .integ_step makes attempts until one step is accepted or
# the integration ends (I$done).
.INTEG_INIT <- quote({
  t0 <- prob$tspan[1L]; tf <- prob$tspan[2L]
  tdir <- sign(tf - t0)
  if (tdir == 0) tdir <- 1
  abstol <- opts$abstol %||% 1e-6
  reltol <- opts$reltol %||% 1e-3
  maxiters <- opts$maxiters %||% 1e5
  dtmax <- opts$dtmax %||% abs(tf - t0)
  # dtmin defaults to 0 (OrdinaryDiffEqCore solve.jl); the step is bounded
  # below by the spacing of t, timedepentdtmin of DiffEqBase
  dtmin <- opts$dtmin %||% 0
  dtmin_at <- function(t) max(dtmin, .Machine$double.eps * abs(t), .Machine$double.xmin)
  is_dae <- inherits(prob, "DAEProblem")
  algs <- if (inherits(alg, "CompositeAlgorithm")) alg$algs else list(alg)
  adaptive <- (opts$adaptive %||% TRUE) && all(vapply(algs, function(a) isTRUE(a$adaptive), TRUE))
  saveat <- sort(unique(as.numeric(opts$saveat %||% numeric())), decreasing = tdir < 0)
  if (length(saveat) == 1L && !is.null(opts$saveat) && length(opts$saveat) == 1L) {
    # initialize_saveat: the range (t0 + s):s:tf of Julia, with the start
    # and the end saved as well
    ds <- tdir * abs(saveat)
    saveat <- c(t0, .jl_range(t0 + ds, ds, tf))
    if (saveat[length(saveat)] != tf) saveat <- c(saveat, tf)
  }
  save_everystep <- opts$save_everystep %||% !length(saveat)
  save_start <- opts$save_start %||% (!length(saveat) || any(abs(saveat - t0) < 1e-14))
  save_end <- opts$save_end %||% (!length(saveat) || any(abs(saveat - tf) < 1e-14))
  # default_linear_interpolation (VerletLeapfrog): no dense output by default
  dense <- opts$dense %||% (!length(saveat) && save_everystep && !isTRUE(alg$linear_interpolation))
  tstops <- as.numeric(c(unlist(opts$tstops), unlist(opts$sys_tstops)))
  tstops <- tstops[(tstops - t0) * tdir > 0 & (tstops - tf) * tdir <= 0]
  tstops <- sort(unique(c(tstops, tf)), decreasing = tdir < 0)

  integ <- new.env(parent = emptyenv())
  # the driver state, for the integrator interface (R/integrator-interface.R)
  integ$.I <- environment()
  integ$t <- t0
  integ$u <- as.numeric(prob$u0)
  names(integ$u) <- names(prob$u0)
  integ$p <- prob$p
  integ$tprev <- t0
  integ$uprev <- integ$u
  integ$terminate <- FALSE
  integ$u_modified <- FALSE
  integ$prob <- prob
  integ$sol_events <- list()
  if (is_dae) integ$du <- as.numeric(prob$du0)
  M <- prob$f$mass_matrix
  mm <- !is_dae && !.is_identity_mass(M)
  ctx <- list(f = prob$f$f, fobj = prob$f, f1 = prob$f$f1, f2 = prob$f$f2, jac1 = prob$f$jac1,
              p = integ$p, mass = if (mm) M else NULL, jac = prob$f$jac,
              tgrad = prob$f$tgrad, abstol = abstol, reltol = reltol, n = length(integ$u),
              adaptive = adaptive, composite = length(algs) > 1L, dae = is_dae,
              algebraic = if (mm) rowSums(abs(M)) == 0 else NULL,
              differential_vars = if (is_dae) prob$differential_vars else
                if (mm) colSums(abs(M)) != 0 else NULL,
              dtmin = dtmin, dtmax = dtmax, lazy_nf = new.env(parent = emptyenv()),
              autodiff = Reduce(function(a, b) a %||% b, lapply(algs, `[[`, "autodiff"), NULL),
              # the dimension of the position part of a DynamicalODEProblem and
              # the count of f2 calls (nf2) of the partitioned methods
              # (R/solvers-symplectic.R)
              nv = prob$nv, counts = new.env(parent = emptyenv()))
  ctx$lazy_nf$n <- 0L
  if (any(!is.finite(integ$u))) stop("initial values are not finite: ",
                                     paste(names(integ$u)[!is.finite(integ$u)], collapse = ", "), call. = FALSE)
  stats <- list(nf = 0L, nf2 = 0L, njacs = 0L, nsolve = 0L, nw = 0L, naccept = 0L, nreject = 0L)
  add_stats <- function(res) {
    for (s in c("nf", "nf2", "njacs", "nsolve", "nw")) stats[[s]] <<- stats[[s]] + as.integer(res[[s]] %||% 0L)
  }
  # the evaluations of a fresh cache (nf2 counts the second part of a split
  # function, as DEStats)
  add_init <- function(cache) {
    stats$nf <<- stats$nf + as.integer(cache$nf0 %||% 1L)
    stats$nf2 <<- stats$nf2 + as.integer(cache$nf20 %||% 0L)
  }
  finish <- function(retcode, ts, us, steps, alg_choice = NULL) {
    # a state resized during the integration (addat, deleteat) leaves NA in
    # the components that did not exist yet
    nmax <- max(c(length(prob$u0), lengths(us)))
    if (length(us) && any(lengths(us) != nmax))
      us <- lapply(us, function(v) c(v, rep(NA_real_, nmax - length(v))))
    U <- if (length(us)) do.call(rbind, us) else matrix(numeric(), 0L, length(prob$u0))
    if (ncol(U) == length(prob$u0)) colnames(U) <- names(prob$u0)
    # stages computed only for interpolation (Verner methods)
    stats$nf <- stats$nf + ctx$lazy_nf$n
    if (!is.null(ctx$counts$nf2)) stats$nf2 <- stats$nf2 + ctx$counts$nf2
    sol <- .new_ode_solution(t = ts, u = U, prob = prob, alg = alg, retcode = retcode,
                             steps = if (dense) steps else NULL, stats = stats, p = integ$p,
                             events = integ$sol_events, dense = dense)
    env <- .sol_env(sol)
    ctx$lazy_nf$sol_env <- env
    if (!is.null(opts$save_idxs)) env$save_idxs <- get_save_idxs_and_saved_subsystem(prob, opts$save_idxs)$save_idxs
    if (length(integ$param_p) > 1L)
      env$param_history <- list(t = integ$param_t, p = do.call(rbind, integ$param_p))
    if (!is.null(alg_choice)) env$alg_choice <- alg_choice
    if (is_dae) class(sol) <- c("DAESolution", class(sol))
    sol
  }

  # OverrideInit (and DefaultInit with an initialization problem): the
  # initialization problem gives u0 and p
  ialg <- .resolve_initializealg(opts$initializealg, prob)
  if (inherits(ialg, "OverrideInit")) {
    iv <- get_initial_values(prob)
    if (!isTRUE(iv$success)) {
      warning("the initialization problem was not solved", call. = FALSE)
      early <- finish("InitialFailure", t0, list(integ$u), list())
      done <- TRUE
    } else {
      integ$u[] <- as.numeric(iv$u0)
      integ$p <- iv$p
      ctx$p <- iv$p
    }
    ialg <- NoInit()
  }
  # consistent algebraic values (BrownFullBasicInit) for mass matrix and
  # fully implicit DAEs, R/solvers-dae.R
  if (!done && (mm || is_dae)) {
    ini <- .dae_initialize(ctx, t0, integ$u, integ$du, ialg, dtmax)
    add_stats(ini)
    if (!ini$ok) {
      warning("could not find consistent initial values for the algebraic variables", call. = FALSE)
      early <- finish("InitialFailure", t0, list(integ$u), list())
      done <- TRUE
    } else {
      integ$u[] <- ini$u
      if (is_dae) integ$du <- ini$du
    }
  }
  if (!done) {

  cbs <- .as_callback_set(opts$callback)
  for (cb in cbs$discrete) {
    cts <- if (!is.null(cb$tstops_fn)) cb$tstops_fn(prob$tspan) else cb$tstops
    tstops <- sort(unique(c(tstops, cts[(cts - t0) * tdir > 0 & (cts - tf) * tdir <= 0])),
                   decreasing = tdir < 0)
  }
  # parameter values in force from each time on, for observed expressions
  integ$param_t <- t0
  integ$param_p <- list(integ$p)
  ts <- numeric(); us <- list()
  steps <- list()
  push_save <- function(t, u) {
    ts[length(ts) + 1L] <<- t
    us[[length(us) + 1L]] <<- u
  }
  if (save_start) push_save(t0, integ$u)
  .init_callbacks(cbs, integ)
  if (integ$u_modified) {
    integ$u_modified <- FALSE
    if (save_start) {
      ts <- ts[-length(ts)]; us <- us[-length(us)]
      push_save(t0, integ$u)
    }
  }
  if (!identical(integ$p, integ$param_p[[1L]])) integ$param_p[[1L]] <- integ$p
  ctx$p <- integ$p
  ctx$du <- integ$du

  # composite algorithms choose the method before each step (CompositeAlgorithm)
  AS <- if (length(algs) > 1L) alg$choice_state() else NULL
  cur <- if (is.null(AS)) 1L else alg$choose(AS, list(init = TRUE, n = ctx$n, reltol = reltol, mass = M))
  caches <- vector("list", length(algs))
  cos <- lapply(algs, .ctrl_opts, opts = opts)
  css <- lapply(algs, function(a) .ctrl_state())
  caches[[cur]] <- algs[[cur]]$init(ctx, integ$t, integ$u)
  add_init(caches[[cur]])
  start_dt <- function(a, cache) {
    if (is_dae) return(1e-6 * abs(tf - t0))
    if (mm) return(max(1e-6, dtmin_at(integ$t)))
    f0 <- cache$fsal %||% ctx$f(integ$u, ctx$p, integ$t)
    .initial_dt(ctx, integ$t, integ$u, f0, tdir, a$order, abstol, reltol, dtmax, dtmin_at(integ$t))
  }
  dt <- opts$dt %||% NULL
  if (is.null(dt)) {
    if (!adaptive && isFALSE(algs[[cur]]$dt_required)) {
      # dt_required(alg) = false: the step reaches the next tstop
      dt <- abs(tf - t0)
    } else {
      if (!adaptive) stop(algs[[cur]]$name, "() is a fixed step method; give dt", call. = FALSE)
      dt <- start_dt(algs[[cur]], caches[[cur]])
      # auto_dt_reset! counts the two evaluations of the initial step estimate
      stats$nf <- stats$nf + 2L
    }
  }
  dt <- abs(dt)
  next_save <- 1L
  retcode <- "Success"
  iter <- 0L
  success_iter <- 0L
  accept_prev <- TRUE
  EEst_prev <- 1
  # discontinuity detection: is_disco_step, disco_checkpoint
  is_disco_step <- FALSE
  disco_checkpoint <- 0
  disco_attempt <- NULL
  eigen_est <- 1
  stop_idx <- 1L
  alg_choice <- if (is.null(AS)) NULL else integer()
  }
})

.integ_init <- function(prob, alg, opts) {
  I <- new.env(parent = environment(.integ_init))
  I$prob <- prob
  I$alg <- alg
  I$opts <- opts
  I$done <- FALSE
  I$disc_next <- FALSE
  I$early <- NULL
  # .INTEG_INIT is the initialization code above, run in the integrator
  eval(.INTEG_INIT, I)
  I$integ$.I <- I
  .integ_bindings(I)
  I
}

# .integ_step and .integ_finalize evaluate their bodies in the integrator
# state environment I (set up by .INTEG_INIT), so these names are bindings
# of I, not globals.
utils::globalVariables(c("AS", "M", "abstol", "adaptive", "add_stats", "alg", "alg_choice", "algs", "cbs",
                         "dense", "dtmax", "dtmin", "dtmin_at", "finish", "is_dae", "maxiters", "mm", "opts", "prob",
                         "push_save", "reltol", "disc_next", "add_init", "save_end", "save_everystep",
                         "save_start", "saveat", "steps", "t0", "tdir", "tf", "ts", "us"))

.integ_step <- function(I) {
  evalq(repeat {
    if (done) break
    if (!((integ$t - tf) * tdir < 0 && !integ$terminate)) {
      done <- TRUE
      break
    }
    if (isTRUE(integ$u_modified)) {
      # the state was changed between steps (set_u, set_t,
      # derivative_discontinuity): the caches and the controller restart
      ctx$p <- integ$p
      ctx$du <- integ$du
      caches[[cur]] <- algs[[cur]]$init(ctx, integ$t, integ$u)
      add_init(caches[[cur]])
      css[[cur]] <- .ctrl_state()
      success_iter <- 0L
      integ$u_modified <- FALSE
    }
    iter <- iter + 1L
    if (iter > maxiters) {
      retcode <- "MaxIters"
      warning("maximum number of iterations reached at t = ", integ$t, call. = FALSE)
      done <- TRUE
      break
    }
    if (!is.null(AS) && iter > 1L) {
      ch <- alg$choose(AS, list(init = FALSE, n = ctx$n, reltol = reltol, mass = M,
                                eigen_est = eigen_est, dt = dt))
      dt <- ch$dt %||% dt
      if (ch$current != cur) {
        cur <- ch$current
        caches[[cur]] <- algs[[cur]]$init(ctx, integ$t, integ$u)
        add_init(caches[[cur]])
      }
    }
    a <- algs[[cur]]; co <- cos[[cur]]
    if (length(integ$tstops_added)) {
      # stops added by callbacks (add_tstop)
      add <- integ$tstops_added
      integ$tstops_added <- NULL
      add <- add[(add - integ$t) * tdir > 0 & (add - tf) * tdir <= 0]
      if (length(add)) {
        tstops <- sort(unique(c(tstops[seq_along(tstops) >= stop_idx], add)), decreasing = tdir < 0)
        stop_idx <- 1L
      }
    }
    while (stop_idx <= length(tstops) && (tstops[stop_idx] - integ$t) * tdir <= 1e-14 * max(1, abs(integ$t)))
      stop_idx <- stop_idx + 1L
    next_stop <- tstops[min(stop_idx, length(tstops))]
    # a method that chooses before each step (StochasticCompositeAlgorithm)
    # sees the step size limited by dtmax and the next stop and may change it
    if (!is.null(a$pre_step))
      dt <- a$pre_step(caches[[cur]], iter, min(dt, dtmax, abs(next_stop - integ$t)), eigen_est, dt, integ)
    hit_stop <- FALSE
    dt <- min(dt, dtmax)
    dtcur <- dt
    # modify_dt_for_tstops!: a step that reaches the next stop to within 100
    # eps ends exactly at the stop, with the shorter of its size and the
    # distance to the stop
    dist <- abs(next_stop - integ$t)
    if (dtcur + 100 * .eps64(max(abs(integ$t), abs(next_stop))) >= dist) {
      dtcur <- min(dtcur, dist)
      hit_stop <- TRUE
    }
    ctx$p <- integ$p
    ctx$iter <- iter
    ctx$success_iter <- success_iter
    # integrator.accept_step: false until the first accepted step
    ctx$accept_prev <- accept_prev && stats$naccept > 0L
    ctx$EEst_prev <- EEst_prev
    ctx$tprev <- integ$tprev
    ctx$uprev2 <- integ$uprev
    ctx$du <- integ$du
    ctx$derivative_discontinuity <- isTRUE(disc_next)
    disc_next <- FALSE
    res <- a$step(ctx, integ$t, integ$u, tdir * dtcur, caches[[cur]])
    add_stats(res)
    if (!is.null(res$eigen_est)) eigen_est <- res$eigen_est
    if (is.environment(res$cache)) caches[[cur]] <- res$cache
    if (isTRUE(res$fail)) {
      # Newton failure (force_stepfail): post_newton_controller!
      if (!adaptive) {
        retcode <- "ConvergenceFailure"
        warning("the nonlinear solver failed at t = ", integ$t, call. = FALSE)
        done <- TRUE
        break
      }
      # res$nreject = 0L: not counted as a rejection (force_stepfail of the
      # Nordsieck methods)
      stats$nreject <- stats$nreject + (res$nreject %||% 1L)
      accept_prev <- FALSE
      # a method that finishes the failed step keeps its error estimate
      if (!is.null(res$EEst)) EEst_prev <- res$EEst
      dt <- if (!is.null(a$newton_fail)) a$newton_fail(caches[[cur]], dtcur, co) else dtcur / co$failfactor
      if (dt <= dtmin_at(integ$t)) {
        retcode <- "DtLessThanMin"
        warning("step size fell below dtmin at t = ", integ$t,
                "; the nonlinear solver does not converge", call. = FALSE)
        done <- TRUE
        break
      }
      next
    }
    EEst <- if (!adaptive) 0 else res$EEst %||% .err_norm(res$err, integ$u, res$u, abstol, reltol)
    if (is.nan(EEst)) EEst <- Inf
    if (adaptive) {
      qq <- .ctrl_q(co, css[[cur]], EEst, success_iter == 0L, dtcur, a, caches[[cur]])
      css[[cur]] <- qq$cs
    }
    EEst_prev <- EEst
    if (adaptive && .ctrl_rejects(co, css[[cur]], EEst)) {
      stats$nreject <- stats$nreject + 1L
      accept_prev <- FALSE
      ddt <- -1
      if (co$discontinuity_detection) {
        # handle_disco_reject!: retry just before the first callback root
        disco_attempt <- list(t0 = integ$t, dt = tdir * dtcur, y0 = integ$u, y1 = res$u, k = res$k,
                              interp = a$interp)
        fd <- .find_discontinuity(cbs, integ, disco_attempt)
        disco_attempt <- NULL
        is_disco_step <- fd$found
        if (fd$found) disco_checkpoint <- integ$t + tdir * dtcur
        ddt <- .disco_dt(fd$theta, dtcur)
      }
      dt <- if (ddt > 0) ddt else .ctrl_reject(co, css[[cur]], dtcur, EEst, a, caches[[cur]])
      if (!is.null(a$on_reject)) a$on_reject(caches[[cur]])
      if (dt <= dtmin_at(integ$t)) {
        retcode <- "DtLessThanMin"
        warning("step size fell below dtmin at t = ", integ$t, "; the problem may be stiff or singular",
                call. = FALSE)
        done <- TRUE
        break
      }
      next
    }
    if (any(!is.finite(res$u))) {
      retcode <- "Unstable"
      warning("the solution became non finite at t = ", integ$t, call. = FALSE)
      done <- TRUE
      break
    }
    stats$naccept <- stats$naccept + 1L
    success_iter <- success_iter + 1L
    after_reject <- !accept_prev
    # integrator.dt after the accepted step
    dt_int <- if (hit_stop && after_reject) max(dt, dtcur) else dtcur
    accept_prev <- TRUE
    if (!is.null(AS)) alg_choice <- c(alg_choice, cur)
    if (adaptive) {
      # a step shortened to hit a tstop proposes from the controller's dt
      # only after a rejected step: after an accepted step OrdinaryDiffEq
      # shortens dt twice (apply_step! and loopheader!), so the dt it
      # restores is the shortened one
      acc <- .ctrl_accept(co, css[[cur]], qq$q, EEst, dt_int,
                          a, caches[[cur]])
      css[[cur]] <- acc$cs
      dtnew <- acc$dt
      if (co$discontinuity_detection && is_disco_step) {
        # handle_disco_accept!: at most a quarter of the way to the end of
        # the rejected step
        is_disco_step <- FALSE
        dtnew <- min(tdir * (disco_checkpoint - (if (hit_stop) next_stop else integ$t + tdir * dtcur)) / 4, dtnew)
      }
      # dtnew_modification: the largest stable step of stabilized methods
      if (!is.null(a$dtnew_mod)) dtnew <- a$dtnew_mod(dtnew, eigen_est)
    } else dtnew <- dt
    tnew <- if (hit_stop) next_stop else integ$t + tdir * dtcur
    # calc_dt_propose!: the step is the interval it spans in floating point
    if (adaptive) dtnew <- abs((tnew + tdir * abs(dtnew)) - tnew)
    # apply_solve_step_limiter!: the limiter may change the accepted state
    # (and the parameters), so the FSAL derivative is evaluated again
    limiter <- opts$step_limiter %||% a$step_limiter
    if (!is.null(limiter)) {
      ul <- as.numeric(limiter(res$u, integ, integ$p, tnew))
      if (length(ul) != length(res$u)) stop("step_limiter must return the state vector", call. = FALSE)
      res$u <- ul
      if (!is_dae && !is.null(res$cache$fsal)) {
        # an environment cache is updated in place, a list cache is replaced
        res$cache$fsal <- ctx$f(ul, integ$p, tnew)
        stats$nf <- stats$nf + 1L
      }
    }
    # the interpolants take the step size of the step (integrator.dt), also
    # when the step ends at a stop
    step <- list(t0 = integ$t, dt = tdir * dtcur, y0 = integ$u, y1 = res$u, k = res$k, interp = a$interp)
    integ$tprev <- integ$t; integ$uprev <- integ$u
    integ$t <- tnew
    integ$u <- res$u
    names(integ$u) <- names(prob$u0)
    if (is_dae) integ$du <- res$du
    integ$last_step <- step
    if (!is.environment(res$cache)) caches[[cur]] <- res$cache
    # callbacks may move t back to an event time and change u
    ev <- .apply_callbacks(cbs, integ, step, ctx)
    step <- integ$last_step
    if (dense) steps[[length(steps) + 1L]] <- step
    # saving
    tend_step <- integ$t
    if (length(saveat)) {
      while (next_save <= length(saveat) && (saveat[next_save] - tend_step) * tdir <= 0) {
        ts_ <- saveat[next_save]
        if (!(abs(ts_ - t0) < 1e-14 && save_start) &&
            (ts_ - step$t0) * tdir >= -1e-14) {
          # savevalues!: a save time that is not the end of the step is
          # interpolated
          theta <- if (step$dt == 0) 1 else (ts_ - step$t0) / step$dt
          push_save(ts_, if (ts_ == tend_step || step$dt == 0) step$y1
                         else step$interp(theta, step$dt, step$y0, step$y1, step$k))
        }
        next_save <- next_save + 1L
      }
    }
    if (save_everystep && (integ$t - tf) * tdir < 0 && !isTRUE(ev$modified) &&
        !(length(ts) && ts[length(ts)] == integ$t)) {
      # after an event the points before and after the affect are saved below
      push_save(integ$t, integ$u)
    }
    if (!identical(integ$p, integ$param_p[[length(integ$param_p)]])) {
      integ$param_t <- c(integ$param_t, integ$t)
      integ$param_p[[length(integ$param_p) + 1L]] <- integ$p
    }
    if (isTRUE(ev$modified)) {
      ctx$p <- integ$p
      if (mm || is_dae) {
        # the affect may break the algebraic constraints: reinitialize
        ini <- .dae_initialize(ctx, integ$t, integ$u, integ$du, opts$initializealg, dtmax)
        add_stats(ini)
        if (!ini$ok) {
          retcode <- "InitialFailure"
          warning("could not find consistent algebraic values after an event at t = ", integ$t, call. = FALSE)
          done <- TRUE
          break
        }
        integ$u[] <- ini$u
        if (is_dae) integ$du <- ini$du
        integ$cb_last <- lapply(cbs$continuous, function(cb) cb$condition(integ$u, integ$t, integ))
      }
      if (isTRUE(ev$force) && !save_everystep && !is.null(ev$left)) {
        # forced saves around the event (DiffEqBase apply_callback!): a
        # saveat point at the event time already holds the state before it
        if (ev$save_pre && !(length(ts) && ts[length(ts)] == ev$left$t)) push_save(ev$left$t, ev$left$u)
        if (ev$save_post) push_save(integ$t, integ$u)
      } else {
        if (save_everystep && !is.null(ev$left)) push_save(ev$left$t, ev$left$u)
        if ((save_everystep && !isFALSE(ev$save_post)) ||
            (!save_everystep && length(saveat) && length(ts) && abs(ts[length(ts)] - integ$t) <= 1e-14 * max(1, abs(integ$t)))) {
          # a saved point at the event time takes the value after the affect
          if (!save_everystep) { ts <- ts[-length(ts)]; us <- us[-length(us)] }
          push_save(integ$t, integ$u)
        }
      }
      ctx$du <- integ$du
      if (!is.null(a$on_event)) {
        # methods with Newton state keep their caches across the event (as
        # OrdinaryDiffEq, which only evaluates the FSAL derivative again)
        # and see a derivative discontinuity on the next step
        evr <- a$on_event(ctx, integ$t, integ$u, caches[[cur]])
        caches[[cur]] <- evr$cache
        stats$nf <- stats$nf + as.integer(evr$nf)
        stats$nf2 <- stats$nf2 + as.integer(evr$nf2 %||% 0L)
        disc_next <- TRUE
      } else {
        caches[[cur]] <- a$init(ctx, integ$t, integ$u)
        add_init(caches[[cur]])
      }
      integ$u_modified <- FALSE
      # the step size carries over the event (as in OrdinaryDiffEq); the
      # caches, the controller history and multistep histories restart
      # a continuous event proposes dtrelax times the step (set_proposed_dt!
      # in DiffEqBase's apply_callback!)
      if (adaptive) dt <- if (!is.null(ev$dtrelax)) max(dtmin, ev$dtrelax * dt_int) else min(abs(dtnew), dtmax)
      if (co$discontinuity_detection && isTRUE(ev$maybe_disco)) {
        # reinit_controller! after an event that may be a discontinuity
        css[[cur]] <- utils::modifyList(css[[cur]], list(q11 = 1, errold = 1e-4, err = c(1, 1, 1), dt_factor = 1,
                                                         dtacc = 1, erracc = 1, qold = 1))
      }
      break
    }
    dt <- min(abs(dtnew), dtmax)
    break
  }, I)
  invisible(I)
}

.integ_finalize <- function(I) {
  if (!is.null(I$early)) return(I$early)
  evalq({
  if (save_end && (!length(ts) || abs(ts[length(ts)] - integ$t) > 1e-14 * max(1, abs(integ$t)) ||
                   integ$terminate)) {
    if (!length(ts) || ts[length(ts)] != integ$t) push_save(integ$t, integ$u)
  }
  if (integ$terminate && retcode == "Success") retcode <- "Terminated"
  # finalize affects run after the last save; parameter changes are recorded
  for (cb in c(cbs$continuous, cbs$discrete)) if (!is.null(cb$finalize)) cb$finalize(integ)
  if (!identical(integ$p, integ$param_p[[length(integ$param_p)]])) {
    integ$param_t <- c(integ$param_t, integ$t)
    integ$param_p[[length(integ$param_p) + 1L]] <- integ$p
  }
  finish(retcode, ts, us, steps, alg_choice)
  }, I)
}
