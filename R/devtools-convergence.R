# Convergence testing of DiffEqDevTools.jl (src/convergence.jl:
# ConvergenceSimulation, ConvergenceTrajectory, test_convergence,
# analyticless_test_convergence; src/test_solution.jl: TestSolution,
# appxtrue; the ConvergenceSimulation plot recipe of src/plotrecipes.jl),
# with the error measures of SciMLBase that they rely on
# (calculate_solution_errors! of src/solutions/ode_solutions.jl,
# dae_solutions.jl and rode_solutions.jl, calculate_ensemble_errors of
# src/ensemble/ensemble_solutions.jl).
#
# Errors are named lists. The Julia keys l\u221e, L\u221e, weak_l\u221e and
# weak_L\u221e are spelled l_inf, L_inf, weak_l_inf and weak_L_inf; the
# functions that take an error name accept either spelling (and linf).
# The order estimates (the Julia field \U0001d4aaest) are in the field Oest.

.dv_key <- function(k) {
  k <- as.character(k)
  switch(k,
    "l\u221E" = , linf = , l_inf = "l_inf",
    "L\u221E" = , Linf = , L_inf = "L_inf",
    "weak_l\u221E" = , weak_linf = , weak_l_inf = "weak_l_inf",
    "weak_L\u221E" = , weak_Linf = , weak_L_inf = "weak_L_inf",
    k)
}

# Julia's mean (sum / n, without R's refinement pass) and 2-norm.
.dv_mean <- function(x) sum(x) / length(x)
.dv_norm <- function(x) {
  if (length(x) == 1L) return(abs(x))
  m <- max(abs(x))
  if (is.na(m) || m == 0 || is.infinite(m)) return(m)
  sqrt(sum(x * x))
}
.dv_last <- function(x) x[[length(x)]]
.dv_rows <- function(u) {
  if (is.matrix(u)) return(u)
  do.call(rbind, lapply(u, function(v) as.numeric(unlist(v))))
}
.dv_as_rows_list <- function(M) lapply(seq_len(nrow(M)), function(i) M[i, ])
.dv_same_t <- function(a, b) length(a) == length(b) && all(a == b)

# collect(range(t0, stop = t1, length = 100)): the correctly rounded points
# for the usual endpoints
.dv_dense_times <- function(t0, t1, n = 100L) {
  i <- 0:(n - 1L)
  (t0 * (n - 1L - i) + t1 * i) / (n - 1L)
}

# Values at the times ts of a solution, a TestSolution or an interpolating
# function, one row per time.
.dv_at <- function(x, ts) {
  f <- if (inherits(x, "TestSolution")) x$interp else x
  v <- f(ts)
  if (is.list(v) && !is.data.frame(v)) return(.dv_rows(v))
  if (is.matrix(v)) return(v)
  matrix(as.numeric(v), nrow = length(ts))
}

.dv_ts_errors <- function(D) {
  list(l2 = sqrt(.dv_mean(as.vector(t(D))^2)), l_inf = max(abs(D)))
}
.dv_dense_errors <- function(D) {
  list(L2 = sqrt(.dv_mean(as.vector(t(D))^2)), L_inf = max(abs(D)))
}

# A copy of a solution with new analytic values and errors (build_solution).
.dv_with_errors <- function(sol, u_analytic, errors) {
  data <- as.list(.sol_env(sol), all.names = TRUE)
  data$obs_cache <- NULL
  data$u_analytic <- u_analytic
  data$errors <- errors
  .new_solution(setdiff(class(sol), c("SciMLSolution", "function")), data)
}

# The analytic solution f.analytic of a problem as a function of t (and of
# the Brownian value W for random problems), or NULL.
.dv_analytic <- function(prob) {
  an <- if (is.list(prob$f)) prob$f$analytic
  if (!is.function(an)) return(NULL)
  if (inherits(prob, "DAEProblem")) return(function(t, W = NULL) an(prob$du0, prob$u0, prob$p, t))
  if (inherits(prob, c("SDEProblem", "RODEProblem", "SDDEProblem")))
    return(function(t, W = NULL) an(prob$u0, prob$p, t, W))
  if (inherits(prob, "DDEProblem")) return(function(t, W = NULL) an(prob$u0, prob$h, prob$p, t))
  function(t, W = NULL) an(prob$u0, prob$p, t)
}

# The Brownian values of a random solution at the times ts: the saved path
# when it was saved at those times, from the noise of the problem when it is
# a fixed path (NoiseGrid, NoiseFunction), otherwise interpolated linearly in
# the saved path (where Julia samples a Brownian bridge). SciMLBase takes
# sol.W.u[i] for the ith saved time, which for a NoiseGrid finer than the
# steps is the grid value at another time; R uses W(t[i]).
.dv_W_rows <- function(env, ts, aligned = TRUE) {
  W <- env$W
  Wm <- if (!is.null(W)) as.matrix(W$W)
  if (aligned && !is.null(Wm) && nrow(Wm) == length(ts) && isTRUE(all(as.numeric(W$t) == ts)))
    return(.dv_as_rows_list(Wm))
  proc <- env$prob$noise_process %||% env$prob$noise
  if (inherits(proc, c("NoiseGrid", "NoiseFunction"))) {
    path <- .np_env(proc)$path
    return(lapply(ts, function(s) as.numeric(path(s, NULL, NULL))))
  }
  if (is.null(W)) stop("the solution has no saved noise path for its analytic solution", call. = FALSE)
  lapply(ts, function(s) {
    hit <- which(W$t == s)
    if (length(hit)) return(Wm[hit[1L], ])
    apply(Wm, 2L, function(col) .approx_hold(W$t, col, s))
  })
}

# calculate_solution_errors! of SciMLBase: the analytic values at the saved
# times, the final error and, with more than two saved points, the
# timeseries errors l2 and l_inf and (for a dense solution) the dense errors
# L2 and L_inf on 100 equally spaced times. For random problems the
# timeseries errors do not need more than two points and the dense errors
# do not need the timeseries errors.
.dv_solution_errors <- function(sol, timeseries_errors = TRUE, dense_errors = FALSE) {
  env <- .sol_env(sol)
  an <- .dv_analytic(env$prob)
  if (is.null(an)) return(sol)
  ts <- env$t
  U <- .dv_rows(sol$u)
  rode <- inherits(sol, "RODESolution")
  Ws <- if (rode) .dv_W_rows(env, ts) else vector("list", length(ts))
  ua <- lapply(seq_along(ts), function(i) as.numeric(an(ts[i], Ws[[i]])))
  A <- .dv_rows(ua)
  n <- nrow(U)
  errors <- list(final = abs(.dv_mean(abs(U[n, ] - A[n, ]))))
  dense <- function() {
    dts <- .dv_dense_times(ts[1L], ts[n])
    IU <- .dv_at(sol, dts)
    if (rode) {
      Wd <- .dv_W_rows(env, dts, aligned = FALSE)
      u1 <- U[1L, ]
      ia <- lapply(seq_along(dts), function(i) as.numeric(env$prob$f$analytic(u1, env$prob$p, dts[i], Wd[[i]])))
    } else {
      ia <- lapply(dts, function(s) as.numeric(an(s)))
    }
    .dv_dense_errors(IU - .dv_rows(ia))
  }
  if ((rode || n > 2L) && isTRUE(timeseries_errors)) {
    errors <- c(errors, .dv_ts_errors(U - A))
    if (!rode && isTRUE(env$dense) && isTRUE(dense_errors)) errors <- c(errors, dense())
  }
  if (rode && isTRUE(dense_errors)) errors <- c(errors, dense())
  .dv_with_errors(sol, ua, errors)
}

# Solve and attach the analytic errors (solve's calculate_error,
# timeseries_errors and dense_errors options).
.dv_solve <- function(prob, alg, opts) {
  ts_err <- opts$timeseries_errors %||% TRUE
  dense_err <- opts$dense_errors %||% FALSE
  calc <- opts$calculate_error %||% TRUE
  for (nm in c("timeseries_errors", "dense_errors", "calculate_error", "weak_timeseries_errors",
               "weak_dense_errors"))
    opts[[nm]] <- NULL
  sol <- if (is.null(alg)) do.call(solve, c(list(prob), opts)) else do.call(solve, c(list(prob, alg), opts))
  if (!isTRUE(calc) || !inherits(sol, "SciMLSolution")) return(sol)
  .dv_solution_errors(sol, ts_err, dense_err)
}

# ---------------------------------------------------------------------------
# Test solutions and approximate errors

#' Errors against a reference solution
#'
#' `TestSolution(t, u)` holds reference values `u` (a list of states, or a
#' matrix with one row per time) at the times `t`; `TestSolution(t, u,
#' interp)` adds an interpolating function `interp(ts)` (such as a
#' solution), and `TestSolution(sol)` wraps the interpolation of a solution
#' alone.
#'
#' `appxtrue(sol, ref)` computes the errors of the solution `sol` against
#' the reference `ref` (a `TestSolution` or a solution): the `final` error
#' (the mean absolute difference at the end), the timeseries errors `l2`
#' and `l_inf` at the times of `sol` and the dense errors `L2` and `L_inf`
#' at 100 equally spaced times. Without an interpolation of the reference
#' only the final error is computed, plus the timeseries errors when the
#' times agree. It returns a copy of `sol` whose `errors` and `u_analytic`
#' fields hold the errors and the reference values. For a reference
#' solution, `timeseries_errors` and `dense_errors` default to whether it
#' is dense.
#'
#' For an `EnsembleSolution` of SDE trajectories, `appxtrue(sim, setup)`
#' solves each trajectory again on its own Brownian path with the
#' algorithm `setup$alg` and the other options of `setup`, and returns the
#' ensemble errors of the trajectories against those solutions (with the
#' options `weak_timeseries_errors` and `weak_dense_errors` in `...`).
#'
#' @param t Times, or a solution for `TestSolution(sol)`.
#' @param u Reference values.
#' @param interp Interpolating function of the times.
#' @param sol A solution or an `EnsembleSolution`.
#' @param sol2 The reference: a `TestSolution`, a solution or, for an
#'   ensemble, a list of solver options with the algorithm `alg`.
#' @param timeseries_errors,dense_errors Which errors to compute against a
#'   reference solution.
#' @param ... Options of the ensemble errors.
#' @return `TestSolution()` returns a `TestSolution`; `appxtrue()` a solution
#'   with errors (an `EnsembleTestSolution` for an ensemble).
#' @examples
#' f <- function(u, p, t) -u
#' prob <- ODEProblem(f, 1, c(0, 1))
#' sol <- solve(prob, Euler(), dt = 1 / 8)
#' ref <- solve(prob, RK4(), dt = 1 / 64)
#' appxtrue(sol, ref)$errors
#' appxtrue(sol, TestSolution(ref))$errors$final
#' @export
TestSolution <- function(t, u = NULL, interp = NULL) {
  if (is.null(u) && is.null(interp) && is.function(t))
    return(structure(list(t = NULL, u = NULL, interp = t, dense = TRUE, retcode = "Success"),
                     class = "TestSolution"))
  if (is.matrix(u)) u <- .dv_as_rows_list(u)
  if (!is.list(u)) u <- as.list(u)
  structure(list(t = as.numeric(t), u = u, interp = interp, dense = !is.null(interp), retcode = "Success"),
            class = "TestSolution")
}

#' @export
print.TestSolution <- function(x, ...) {
  cat("TestSolution", if (!is.null(x$t)) paste0(" with ", length(x$t), " time points"),
      if (x$dense) " and an interpolation", "\n", sep = "")
  invisible(x)
}

#' @rdname TestSolution
#' @export
appxtrue <- function(sol, sol2, timeseries_errors = NULL, dense_errors = NULL, ...) {
  if (inherits(sol, "EnsembleSolution")) return(.dv_appxtrue_ensemble(sol, sol2, ...))
  ts <- sol$t
  U <- .dv_rows(sol$u)
  n <- nrow(U)
  if (inherits(sol2, "TestSolution")) {
    ref <- if (is.null(sol2$u) && !is.null(sol2$interp))
      TestSolution(ts, .dv_as_rows_list(.dv_at(sol2, ts)), sol2$interp) else sol2
    errors <- list(final = .dv_mean(abs(U[n, ] - as.numeric(.dv_last(ref$u)))))
    if (isTRUE(ref$dense)) {
      TA <- .dv_at(ref, ts)
      ua <- .dv_as_rows_list(TA)
      errors <- c(errors, .dv_ts_errors(U - TA))
      dts <- .dv_dense_times(ts[1L], ts[n])
      errors <- c(errors, .dv_dense_errors(.dv_at(sol, dts) - .dv_at(ref, dts)))
    } else {
      ua <- sol2$u
      if (.dv_same_t(ts, sol2$t)) errors <- c(errors, .dv_ts_errors(U - .dv_rows(sol2$u)))
    }
    return(.dv_with_errors(sol, ua, errors))
  }
  # solutions of random problems are not dense (the dense = false default of
  # StochasticDiffEq) though they interpolate
  dense2 <- !inherits(sol2, "RODESolution") && isTRUE(.sol_env(sol2)$dense)
  timeseries_errors <- timeseries_errors %||% dense2
  dense_errors <- dense_errors %||% dense2
  errors <- list(final = .dv_mean(abs(U[n, ] - as.numeric(.dv_last(sol2$u)))))
  if (dense2) {
    TA <- .dv_at(sol2, ts)
    ua <- .dv_as_rows_list(TA)
    errors <- c(errors, .dv_ts_errors(U - TA))
    if (isTRUE(dense_errors)) {
      dts <- .dv_dense_times(ts[1L], ts[n])
      errors <- c(errors, .dv_dense_errors(.dv_at(sol, dts) - .dv_at(sol2, dts)))
    }
  } else {
    ua <- sol2$u
    if (isTRUE(timeseries_errors) && .dv_same_t(ts, sol2$t))
      errors <- c(errors, .dv_ts_errors(U - .dv_rows(sol2$u)))
  }
  .dv_with_errors(sol, ua, errors)
}

# appxtrue(sim::EnsembleSolution, appx_setup): each trajectory against a
# solution on its own Brownian path. (The Julia method cannot run: it calls
# the removed constructor Vector{T}(n), and drops the parameters.)
.dv_appxtrue_ensemble <- function(sim, setup, ...) {
  setup <- as.list(setup)
  alg <- setup$alg
  setup$alg <- NULL
  new <- lapply(sim$u, function(s) {
    prob <- .sol_env(s)$prob
    nrp <- if (identical(prob$noise_type, "general")) matrix(0, length(prob$u0), prob$f$nbrownians)
    prob2 <- SDEProblem(prob$f, prob$f$g, prob$u0, prob$tspan, prob$p, noise = NoiseWrapper(s$W),
                        noise_rate_prototype = nrp)
    true_sol <- do.call(solve, c(list(prob2, alg), setup))
    appxtrue(s, true_sol)
  })
  .dv_ensemble_errors(new, elapsed = sim$elapsed %||% 0, converged = sim$converged %||% FALSE, ...)
}

# ---------------------------------------------------------------------------
# Ensemble errors (calculate_ensemble_errors)

.dv_vec_mean <- function(vs) Reduce(`+`, lapply(vs, as.numeric)) / length(vs)

.dv_ensemble_errors <- function(u, elapsed = 0, converged = FALSE, weak_timeseries_errors = FALSE,
                                weak_dense_errors = FALSE) {
  keys <- names(u[[1L]]$errors)
  errors <- stats::setNames(lapply(keys, function(k) vapply(u, function(s) s$errors[[k]], 0)), keys)
  error_means <- lapply(errors, .dv_mean)
  error_medians <- lapply(errors, stats::median)
  m_final <- .dv_vec_mean(lapply(u, function(s) .dv_last(s$u)))
  m_an <- .dv_vec_mean(lapply(u, function(s) .dv_last(s$u_analytic)))
  weak <- list(weak_final = .dv_norm(m_final - m_an))
  rms <- function(e) sqrt(sum(e * e) / length(e))
  if (isTRUE(weak_timeseries_errors)) {
    tsw <- lapply(seq_along(u[[1L]]$u), function(i)
      .dv_vec_mean(lapply(u, function(s) as.numeric(s$u[[i]]) - as.numeric(s$u_analytic[[i]]))))
    weak$weak_l2 <- rms(vapply(tsw, rms, 0))
    weak$weak_l_inf <- max(vapply(tsw, function(e) max(abs(e)), 0))
  }
  if (isTRUE(weak_dense_errors)) {
    t1 <- u[[1L]]$t
    dts <- .dv_dense_times(t1[1L], t1[length(t1)])
    an <- lapply(u, function(s) {
      env <- .sol_env(s)
      Wd <- .dv_W_rows(env, dts, aligned = FALSE)
      lapply(seq_along(dts), function(i) as.numeric(env$prob$f$analytic(env$prob$u0, env$prob$p, dts[i], Wd[[i]])))
    })
    ud <- lapply(u, function(s) .dv_at(s, dts))
    dw <- lapply(seq_along(dts), function(i)
      .dv_vec_mean(lapply(seq_along(u), function(j) ud[[j]][i, ] - an[[j]][[i]])))
    weak$weak_L2 <- rms(vapply(dw, rms, 0))
    weak$weak_L_inf <- max(vapply(dw, function(e) max(abs(e)), 0))
  }
  structure(list(u = u, errors = errors, weak_errors = weak, error_means = error_means,
                 error_medians = error_medians, elapsed = elapsed, converged = converged),
            class = c("EnsembleTestSolution", "EnsembleSolution"))
}

# ---------------------------------------------------------------------------
# Convergence simulations

#' Results of a convergence test
#'
#' `ConvergenceSimulation(solutions, convergence_axis)` collects the errors
#' of solutions computed at the points of a convergence axis (such as the
#' time steps): `errors[[k]]` holds the mean error `k` of each solution
#' (with `additional_errors`, a named list of error series, added), and
#' `Oest[[k]]` the estimated order of convergence, the absolute mean of
#' \eqn{\log_2(e_{i+1}/e_i)} over adjacent points. Give `expected_value`
#' when all error series come from `additional_errors`.
#'
#' `ConvergenceTrajectory(t, u, u_analytic, errors)` is the part of an SDE
#' trajectory a convergence test keeps (its final time, final state, final
#' analytic state and errors) when the solutions are not retained.
#'
#' `length(sim)` is the number of solutions and `sim[[i]]` the `i`th one;
#' `plot(sim)` draws the errors against the convergence axis on log scales.
#'
#' @param solutions List of solutions with errors.
#' @param convergence_axis Values of the convergence axis.
#' @param auxdata Extra data stored in the result.
#' @param additional_errors Named list of error series.
#' @param expected_value The value of the weak errors, which skips the errors
#'   of the solutions.
#' @param t,u,u_analytic Final time, state and analytic state.
#' @param errors Named list of errors.
#' @return A `ConvergenceSimulation` (a list with `solutions`, `errors`,
#'   `N`, `auxdata`, `Oest` and `convergence_axis`), or a
#'   `ConvergenceTrajectory`.
#' @examples
#' sols <- list(list(errors = list(final = 0.1)), list(errors = list(final = 0.025)))
#' sim <- ConvergenceSimulation(sols, c(0.1, 0.05))
#' sim$Oest$final
#' @export
ConvergenceSimulation <- function(solutions, convergence_axis, auxdata = NULL, additional_errors = NULL,
                                  expected_value = NULL) {
  errors <- list()
  if (is.null(expected_value)) {
    e1 <- solutions[[1L]]$errors
    if (!length(e1))
      stop("Errors dictionary is empty. No analytical solution set. If you used `test_convergence` ",
           "you may consider `analyticless_test_convergence` instead.", call. = FALSE)
    for (k in names(e1))
      errors[[k]] <- vapply(solutions, function(s) {
        e <- s$errors[[k]]
        if (is.null(e)) stop("a solution has no error ", k, call. = FALSE)
        .dv_mean(e)
      }, 0)
  }
  for (k in names(additional_errors)) errors[[.dv_key(k)]] <- as.numeric(unlist(additional_errors[[k]]))
  Oest <- lapply(errors, function(e) abs(.dv_mean(log2(e[-1L] / e[-length(e)]))))
  structure(list(solutions = solutions, errors = errors, N = length(solutions), auxdata = auxdata,
                 Oest = Oest, convergence_axis = convergence_axis),
            class = "ConvergenceSimulation")
}

#' @rdname ConvergenceSimulation
#' @export
ConvergenceTrajectory <- function(t, u, u_analytic, errors)
  structure(list(t = list(t), u = list(u), u_analytic = list(u_analytic), errors = as.list(errors)),
            class = "ConvergenceTrajectory")

#' @export
length.ConvergenceSimulation <- function(x) .subset2(x, "N")

#' @export
`[[.ConvergenceSimulation` <- function(x, i, ...) {
  if (is.numeric(i)) return(.subset2(x, "solutions")[[i]])
  .subset2(x, i)
}

#' @export
print.ConvergenceSimulation <- function(x, ...) {
  cat("ConvergenceSimulation with ", x$N, " solutions\n", sep = "")
  for (k in names(x$errors))
    cat("  ", k, ": order ", format(x$Oest[[k]], digits = 4), "; errors ",
        paste(format(x$errors[[k]], digits = 4), collapse = " "), "\n", sep = "")
  invisible(x)
}

#' @export
print.ConvergenceTrajectory <- function(x, ...) {
  cat("ConvergenceTrajectory at t = ", format(x$t[[1L]]), " with errors ",
      paste(names(x$errors), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @export
print.EnsembleTestSolution <- function(x, ...) {
  cat("EnsembleTestSolution with ", length(x$u), " trajectories\n", sep = "")
  for (k in names(x$error_means)) cat("  ", k, " mean ", format(x$error_means[[k]], digits = 4), "\n", sep = "")
  for (k in names(x$weak_errors)) cat("  ", k, " ", format(x$weak_errors[[k]], digits = 4), "\n", sep = "")
  invisible(x)
}

#' @rdname ConvergenceSimulation
#' @param x A `ConvergenceSimulation`.
#' @param xlab,ylab,log,type,lwd,pch Plot settings.
#' @param ... Passed to [graphics::matplot()].
#' @export
plot.ConvergenceSimulation <- function(x, xlab = "Convergence Axis", ylab = "Error", log = "xy", type = "b",
                                       lwd = 3, pch = seq_along(x$errors), ...) {
  vals <- do.call(cbind, x$errors)
  graphics::matplot(x$convergence_axis, vals, xlab = xlab, ylab = ylab, log = log, type = type, lwd = lwd,
                    pch = pch, lty = 1, col = seq_along(x$errors), ...)
  graphics::legend("topleft", legend = names(x$errors), col = seq_along(x$errors), lty = 1, pch = pch,
                   bty = "n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Convergence tests

#' Convergence tests
#'
#' `test_convergence(dts, prob, alg)` solves the problem with each time step
#' of `dts` (fixed steps: `adaptive = FALSE`, `save_everystep = TRUE` by
#' default) and returns the [ConvergenceSimulation()] of the errors
#' against the analytic solution of the problem (the `analytic` field of its
#' function, `analytic(u0, p, t)`, `analytic(du0, u0, p, t)` for DAEs,
#' `analytic(u0, h, p, t)` for DDEs and `analytic(u0, p, t, W)` for SDE,
#' RODE and SDDE problems). The options `timeseries_errors` (default `TRUE`)
#' and `dense_errors` (default `FALSE`) select the errors; other options go
#' to `solve()`.
#'
#' For SDE, RODE and SDDE problems (or an [EnsembleProblem()] of one),
#' `trajectories` paths are solved for each time step, with the ensemble
#' algorithm `ensemblealg` (the fourth argument), and the errors are the
#' means over the trajectories plus the weak errors `weak_final` (and
#' `weak_l2`, `weak_l_inf` with `weak_timeseries_errors = TRUE`, `weak_L2`,
#' `weak_L_inf` with `weak_dense_errors = TRUE`). With `expected_value`, the
#' weak error is that of the trajectory values (the outputs of the ensemble)
#' against it. Unless `retain_solutions = TRUE` each trajectory is reduced to
#' a [ConvergenceTrajectory()] and only the first one is kept.
#'
#' `test_convergence(probs, convergence_axis, alg)` solves each problem of
#' the list `probs` with the options in `...`, and
#' `test_convergence(setup, alg)` takes `probs` and `convergence_axis` from
#' a list `setup` with those fields.
#'
#' `analyticless_test_convergence(dts, prob, alg, setup)` compares an ODE
#' solution for each time step with a reference solution computed by
#' `solve(prob, setup$alg, ...)` with the other options of `setup`. For a
#' random problem, the fourth argument is the time step of the reference: each
#' of `trajectories` Brownian paths on that grid (a [NoiseGrid()], or with
#' `use_noise_grid = FALSE` the path of the reference solve through a
#' [NoiseWrapper()]) is solved with every step of `dts` and compared with the
#' solution on the reference step.
#'
#' @param dts Time steps, a list of problems, or a convergence setup.
#' @param prob The problem, the convergence axis for a list of problems, or
#'   the algorithm for a setup.
#' @param alg The algorithm.
#' @param ... Options, see Details, and solver options.
#' @param test_dt The reference setup (ODE) or reference time step (SDE).
#' @param dt,p The solver options `dt` and `p`, named in the signature so
#'   that they do not partially match `dts` and `prob`.
#' @return A [ConvergenceSimulation()].
#' @examples
#' f <- function(u, p, t) -u
#' fn <- ODEFunction(f, analytic = function(u0, p, t) u0 * exp(-t))
#' prob <- ODEProblem(fn, 1, c(0, 1))
#' sim <- test_convergence(1 / 2^(3:5), prob, Euler())
#' sim$Oest$final
#' sim2 <- analyticless_test_convergence(1 / 2^(3:5), ODEProblem(f, 1, c(0, 1)), Heun(),
#'   list(alg = RK4(), dt = 1 / 256))
#' sim2$Oest$final
#' @export
test_convergence <- function(dts, prob, alg, ..., dt = NULL, p = NULL) {
  # dt and p are solver options, named here so that they do not partially
  # match dts and prob
  opts <- .dv_dots(list(...), dt, p)
  if (inherits(dts, "ConvergenceSetup") || (is.list(dts) && !is.null(dts[["probs"]]) &&
                                             !is.null(dts[["convergence_axis"]])))
    return(.dv_test_convergence_probs(dts[["probs"]], dts[["convergence_axis"]], prob, opts))
  if (is.list(dts) && length(dts) && all(vapply(dts, inherits, TRUE, "SciMLProblem")))
    return(.dv_test_convergence_probs(dts, prob, alg, opts))
  if (inherits(prob, c("SDEProblem", "RODEProblem", "SDDEProblem", "EnsembleProblem")))
    return(.dv_test_convergence_random(dts, prob, alg, opts))
  .dv_test_convergence_ode(dts, prob, alg, opts)
}

# The options of a call plus the exactly named dt and p.
.dv_dots <- function(opts, dt = NULL, p = NULL) {
  if (!is.null(dt)) opts$dt <- dt
  if (!is.null(p)) opts$p <- p
  opts
}

# The settings of a function among its options (exact names, with their
# defaults) and the remaining options.
.dv_settings <- function(opts, defaults) {
  nm <- names(opts) %||% rep("", length(opts))
  list(s = .dv_merge(defaults, opts[nm %in% names(defaults)]), rest = opts[!nm %in% names(defaults)])
}

.dv_test_convergence_probs <- function(probs, convergence_axis, alg, opts)
  ConvergenceSimulation(lapply(probs, .dv_solve, alg = alg, opts = opts), convergence_axis)

.dv_test_convergence_ode <- function(dts, prob, alg, opts) {
  o <- .dv_merge(list(save_everystep = TRUE, adaptive = FALSE), opts)
  sols <- lapply(dts, function(dt) .dv_solve(prob, alg, .dv_merge(o, list(dt = dt))))
  ConvergenceSimulation(sols, dts, auxdata = list(dts = dts))
}

# _convergence_output_func: a trajectory reduced to its errors and endpoints
.dv_convergence_output_func <- function(sol, i) {
  ua <- sol$u_analytic
  ua <- if (is.null(ua)) .dv_last(sol$u) else .dv_last(ua)
  list(ConvergenceTrajectory(.dv_last(sol$t), .dv_last(sol$u), ua, sol$errors), FALSE)
}

.dv_test_convergence_random <- function(dts, prob, alg, opts) {
  # the ensemble algorithm is the fourth argument, as in Julia
  nm <- names(opts) %||% rep("", length(opts))
  pos <- which(!nzchar(nm))
  if (length(pos)) {
    ea <- opts[[pos[1L]]]
    opts <- opts[-pos[1L]]
    opts$ensemblealg <- ea
  }
  st <- .dv_settings(opts, list(ensemblealg = EnsembleThreads(), trajectories = NULL, save_start = TRUE,
                                save_everystep = TRUE, timeseries_errors = NULL, dense_errors = FALSE,
                                adaptive = FALSE, weak_timeseries_errors = FALSE, weak_dense_errors = FALSE,
                                expected_value = NULL, retain_solutions = FALSE))
  s <- st$s
  trajectories <- s$trajectories
  if (is.null(trajectories)) stop("test_convergence() of a random problem needs `trajectories`", call. = FALSE)
  save_everystep <- isTRUE(s$save_everystep)
  timeseries_errors <- s$timeseries_errors %||% save_everystep
  dense_errors <- s$dense_errors
  weak_timeseries_errors <- isTRUE(s$weak_timeseries_errors)
  weak_dense_errors <- isTRUE(s$weak_dense_errors)
  expected_value <- s$expected_value
  retain_solutions <- isTRUE(s$retain_solutions)
  N <- length(dts)
  is_ens <- inherits(prob, "EnsembleProblem")
  reduce <- !retain_solutions && !is_ens && is.null(expected_value) && !weak_timeseries_errors &&
    !weak_dense_errors
  ens <- if (is_ens) prob else if (reduce) EnsembleProblem(prob, output_func = .dv_convergence_output_func) else
    EnsembleProblem(prob)
  # the errors of a trajectory are computed by its solve, before output_func
  out <- ens$output_func
  ens$output_func <- function(sol, i) {
    if (inherits(sol, "SciMLSolution")) sol <- .dv_solution_errors(sol, timeseries_errors, dense_errors)
    out(sol, i)
  }
  sols <- vector("list", N)
  for (i in seq_len(N)) {
    sim <- do.call(solve, c(list(ens, alg, s$ensemblealg),
                            .dv_merge(st$rest, list(dt = dts[i], adaptive = s$adaptive, save_start = s$save_start,
                                                    save_everystep = save_everystep,
                                                    trajectories = as.integer(trajectories)))))
    message("dt: ", dts[i], " (", i, "/", N, ")")
    sols[[i]] <- if (is.null(expected_value)) {
      summ <- .dv_ensemble_errors(sim$u, sim$elapsed, sim$converged, weak_timeseries_errors, weak_dense_errors)
      if (reduce) summ$u <- summ$u[seq_len(min(1L, length(summ$u)))]
      summ
    } else sim
  }
  add <- list()
  if (is.null(expected_value)) {
    for (k in names(sols[[1L]]$weak_errors)) add[[k]] <- vapply(sols, function(s) s$weak_errors[[k]], 0)
  } else {
    flat <- function(v) if (inherits(v, "SciMLSolution")) unlist(v$u) else as.numeric(unlist(v))
    scalar <- length(expected_value) == 1L
    vals <- vapply(sols, function(s) {
      d <- lapply(s$u, function(v) flat(v) - expected_value)
      # norm(mean(sol.u .- ev)) for a number; norm(mean(sol .- ev)), the mean of
      # all entries, otherwise
      if (scalar) .dv_norm(.dv_vec_mean(d)) else abs(.dv_mean(unlist(d)))
    }, 0)
    add[[if (scalar) "weak_final" else "weak_l2"]] <- vals
  }
  ConvergenceSimulation(sols, dts, auxdata = list(dts = dts), additional_errors = add,
                        expected_value = expected_value)
}

#' @rdname test_convergence
#' @export
analyticless_test_convergence <- function(dts, prob, alg, test_dt, ..., dt = NULL, p = NULL) {
  opts <- .dv_dots(list(...), dt, p)
  if (inherits(prob, c("SDEProblem", "RODEProblem", "SDDEProblem")))
    return(.dv_analyticless_random(dts, prob, alg, test_dt, opts))
  .dv_analyticless_ode(dts, prob, alg, test_dt, opts)
}

.dv_analyticless_ode <- function(dts, prob, alg, appxsol_setup, opts) {
  setup <- as.list(appxsol_setup)
  ref_alg <- setup$alg
  setup$alg <- NULL
  true_sol <- do.call(solve, c(list(prob, ref_alg), setup))
  o <- .dv_merge(list(save_everystep = TRUE, adaptive = FALSE), opts)
  sols <- lapply(dts, function(dt) appxtrue(.dv_solve(prob, alg, .dv_merge(o, list(dt = dt))), true_sol))
  ConvergenceSimulation(sols, dts, auxdata = list(dts = dts))
}

.dv_analyticless_random <- function(dts, prob, alg, test_dt, opts) {
  s <- .dv_settings(opts, list(trajectories = 100L, save_everystep = TRUE, timeseries_errors = NULL,
                               adaptive = FALSE, weak_timeseries_errors = FALSE, weak_dense_errors = FALSE,
                               use_noise_grid = TRUE, verbose = TRUE))$s
  trajectories <- s$trajectories
  adaptive <- s$adaptive
  weak_timeseries_errors <- isTRUE(s$weak_timeseries_errors)
  weak_dense_errors <- isTRUE(s$weak_dense_errors)
  use_noise_grid <- isTRUE(s$use_noise_grid)
  verbose <- isTRUE(s$verbose)
  tmp <- matrix(list(), trajectories, length(dts))
  nrp <- prob$noise_rate_prototype
  m <- if (identical(prob$noise_type, "general")) prob$f$nbrownians else length(prob$u0)
  for (j in seq_len(trajectories)) {
    if (verbose) message("Monte Carlo iteration: ", j, "/", trajectories)
    if (use_noise_grid) {
      t <- seq(prob$tspan[1L], prob$tspan[2L], by = test_dt)
      path <- function() {
        incs <- lapply(seq_len(length(t) - 1L), function(i) sqrt(test_dt) * stats::rnorm(m))
        Reduce(`+`, incs, numeric(m), accumulate = TRUE)
      }
      W <- path()
      Z <- path()
      prob_j <- prob
      prob_j$noise_process <- NoiseGrid(t, W, Z)
      if (inherits(prob, "RODEProblem")) prob_j$noise <- prob_j$noise_process
      true_sol <- solve(prob_j, alg, adaptive = adaptive, dt = test_dt)
      for (i in seq_along(dts)) {
        sol <- solve(prob_j, alg, dt = dts[i], adaptive = adaptive)
        tmp[[j, i]] <- appxtrue(sol, true_sol)
      }
    } else {
      true_sol <- solve(prob, alg, adaptive = adaptive, dt = test_dt, save_noise = TRUE)
      for (i in seq_along(dts)) {
        prob_i <- prob
        prob_i$noise_process <- NoiseWrapper(true_sol$W)
        if (inherits(prob, "RODEProblem")) prob_i$noise <- prob_i$noise_process
        sol <- solve(prob_i, alg, dt = dts[i], adaptive = adaptive)
        tmp[[j, i]] <- appxtrue(sol, true_sol)
      }
    }
  }
  sols <- lapply(seq_along(dts), function(i)
    .dv_ensemble_errors(tmp[, i], 0, TRUE, weak_timeseries_errors, weak_dense_errors))
  add <- list()
  for (k in names(sols[[1L]]$weak_errors)) add[[k]] <- vapply(sols, function(s) s$weak_errors[[k]], 0)
  ConvergenceSimulation(sols, dts, auxdata = list(dts = dts), additional_errors = add)
}
