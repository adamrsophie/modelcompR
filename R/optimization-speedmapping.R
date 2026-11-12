# OptimizationSpeedMapping (lib/OptimizationSpeedMapping/src, v0.2.9) and the
# minimization path of SpeedMapping.jl (src/acx.jl and SpeedMapping_common.jl,
# v0.4): alternating cyclic extrapolations (ACX) of gradient descent steps,
# with orders (2, 3, 3), the initial learning rate found by a bracketing
# search, and box constraints. The dot products are those of BLAS (fused
# multiply adds, .ojl_dot), so the iterates agree with Julia to the last bit.

#' SpeedMapping
#'
#' `SpeedMappingOpt()` minimizes with the alternating cyclic extrapolation
#' (ACX) method of SpeedMapping.jl: cycles of two or three gradient descent
#' steps followed by a squared or cubic extrapolation, with an adaptive
#' learning rate and box constraints (`lb`, `ub`). It stops when the norm of
#' the gradient (of its components that do not push against a bound) is
#' below `1e-8`. `maxiters` limits the number of gradient evaluations and
#' `maxtime` the time in seconds; `abstol` and `reltol` are not used.
#'
#' The return code is `"Success"` when the gradient is small, `"MaxIters"`
#' when a limit was reached, `"MaxTime"` and otherwise `"Failure"`. As in
#' Julia the statistics only carry the time; `sol$original` has the counts
#' (`maps`, the gradient evaluations, `f_calls`, `iterations`), `status`,
#' `residual_norm` and `last_learning_rate`. The gradient of the problem is
#' used, finite differences when it has none.
#'
#' @return An optimization algorithm.
#' @examples
#' f <- function(u, p) (1 - u[1])^2 + 100 * (u[2] - u[1]^2)^2
#' g <- function(u, p) c(-2 * (1 - u[1]) - 400 * u[1] * (u[2] - u[1]^2), 200 * (u[2] - u[1]^2))
#' sol <- solve(OptimizationProblem(f, c(-1.2, 1), grad = g), SpeedMappingOpt())
#' sol$u
#' sol$original$maps
#' @export
SpeedMappingOpt <- function() {
  a <- .opt_alg("SpeedMappingOpt")
  a$solver <- .speedmapping_solve
  a$traits <- list(grad = TRUE, bounds = TRUE, callback = FALSE)
  a
}

.speedmapping_solve <- function(prob, alg, opts) {
  fs <- .opt_fns(prob)
  if (is.function(prob$f$cons))
    stop(IncompatibleOptimizerError(paste0("The algorithm SpeedMappingOpt does not support constraints. Either ",
                                           "remove the `cons` function passed to `OptimizationFunction` or use ",
                                           "a different algorithm.")))
  maxiters <- .opt_maxiters(opts$maxiters)
  if (!is.null(maxiters)) message("maxiters defines maximum gradient calls for SpeedMappingOpt()")
  if (!is.null(opts$abstol)) message("common abstol is currently not used by SpeedMappingOpt()")
  if (!is.null(opts$reltol)) message("common reltol is currently not used by SpeedMappingOpt()")
  x0 <- as.numeric(prob$u0)
  lower <- if (is.null(prob$lb)) NULL else as.numeric(prob$lb)
  upper <- if (is.null(prob$ub)) NULL else as.numeric(prob$ub)
  if (!is.null(lower) && !is.null(upper)) {
    bad <- which(!(lower <= upper))
    if (length(bad))
      stop("bound constraint #", bad[1L], " is infeasible because ", lower[bad[1L]], " > ", upper[bad[1L]], ".",
           call. = FALSE)
  }
  t0 <- proc.time()[["elapsed"]]
  res <- .speedmapping_acx(fs$f, fs$g, x0, lower, upper, maps_limit = maxiters %||% 1e9,
                           time_limit = opts$maxtime %||% Inf)
  retcode <- switch(res$status, first_order = "Success", max_time = "MaxTime", max_iter = , max_eval = "MaxIters",
                    "Failure")
  sol <- .opt_solution(prob, alg, res$minimizer, fs$f(res$minimizer), retcode, res)
  .opt_set_fields(sol, stats = OptimizationStats(time = proc.time()[["elapsed"]] - t0))
}

# isbad of SpeedMapping: a non finite sum
.speedmapping_isbad <- function(x) {
  s <- if (length(x) == 1L) x else sum(x)
  is.na(s) || is.infinite(s)
}

# do_after_map! of a gradient: the descent step x_in + alpha * r (alpha is
# negative), held inside the bounds, and the squared norm of the residual of
# the components that are not at a bound
.speedmapping_after_map <- function(r, alpha, x_in, lower, upper, abstol) {
  x_out <- x_in + alpha * r
  err <- r
  if (!is.null(lower)) {
    for (i in seq_along(x_out)) {
      l <- x_in[i] * 0 + lower[i] * 1
      if (x_out[i] <= l) {
        x_out[i] <- l
        r[i] <- (l - x_in[i]) / alpha
        if (x_in[i] - lower[i] < abstol) err[i] <- (lower[i] - x_in[i]) / alpha
      }
    }
  }
  if (!is.null(upper)) {
    for (i in seq_along(x_out)) {
      u <- x_in[i] * 0 + upper[i] * 1
      if (x_out[i] >= u) {
        x_out[i] <- u
        r[i] <- (u - x_in[i]) / alpha
        if (upper[i] - x_in[i] < abstol) err[i] <- (upper[i] - x_in[i]) / alpha
      }
    }
  }
  list(x = x_out, r = r, rr = .ojl_dot(err, err))
}

.speedmapping_acx <- function(f, G, x_in, lower, upper, orders = c(2L, 3L, 3L), abstol = 1e-8,
                              initial_learning_rate = 1, maps_limit = 1e9, iter_limit = 1e9, time_limit = Inf,
                              reltol_resid_grow = 1e5, abstol_obj_grow = sqrt(abstol), buffer = 0) {
  abstolsq <- abstol * abstol
  max_time <- if (time_limit < Inf) proc.time()[["elapsed"]] + time_limit else Inf
  maps_limit <- ceiling(min(maps_limit, .Machine$integer.max))
  x_now <- x_in
  if (!is.null(lower)) x_now <- pmax(lower, x_now)
  if (!is.null(upper)) x_now <- pmin(upper, x_now)
  has_cons <- !is.null(lower) || !is.null(upper)
  n <- length(x_now)
  x_next <- x_best <- r1 <- r2 <- r3 <- numeric(n)
  reltol_rr <- min(reltol_resid_grow^2, .Machine$double.xmax)
  maps <- f_calls <- iter <- i_ord <- 0L
  rr_now <- rr_best <- lpnormsq_best <- Inf
  alpha <- sigma <- 1
  am <- function(r, a, x) .speedmapping_after_map(r, a, x, lower, upper, abstol)
  result <- function(x, status) list(minimizer = x, residual_norm = sqrt(lpnormsq_best), maps = maps,
                                     f_calls = f_calls, iterations = iter, status = status, algo = "acx",
                                     last_learning_rate = -alpha)
  # the initial learning rate: the gradient norm must not grow too fast, the
  # extrapolation step must not be too large and the objective must decrease
  check_obj <- !is.null(f)
  if (check_obj) {
    obj_in <- obj_now <- f(x_now)
    f_calls <- f_calls + 1L
    max_new_obj_val <- obj_now + abstol_obj_grow
  }
  r1 <- G(x_now)
  maps <- maps + 1L
  if (has_cons) r3 <- r1
  log_a <- log(initial_learning_rate)
  low <- -Inf
  up <- Inf
  tries <- 0L
  while (tries < 100L && up - low > 0.001) {
    alpha <- -exp(log_a)
    if (has_cons) r1 <- r3
    m <- am(r1, alpha, x_now)
    x_next <- m$x
    r1 <- m$r
    rr_best <- m$rr
    if (.speedmapping_isbad(rr_best)) stop("Initial residual has NaN or Inf gradient", call. = FALSE)
    lpnormsq_best <- rr_best
    if (lpnormsq_best < abstolsq) {
      maps <- 1L
      f_calls <- 0L
      return(result(x_now, "first_order"))
    }
    if (check_obj) {
      obj_now <- f(x_next)
      f_calls <- f_calls + 1L
      if (!isTRUE(obj_now < max_new_obj_val)) {
        up <- log_a
        log_a <- if (low == -Inf) log_a - 5 else (log_a + low) * 0.5
      } else if (obj_now < obj_in + 0.8 * alpha * rr_best) {
        low <- log_a
        log_a <- if (up == Inf) log_a + 5 else (log_a + up) * 0.5
      } else {
        check_obj <- FALSE
      }
    }
    if (!check_obj) {
      m <- am(G(x_next), alpha, x_next)
      x_best <- m$x
      r2 <- m$r
      rr_now <- m$rr
      maps <- maps + 1L
      x_best <- r2 - r1
      if (.speedmapping_isbad(rr_now) || rr_now > rr_best * 4) {
        up <- log_a
        log_a <- if (low == -Inf) log_a - 5 else (log_a + low) * 0.5
      } else if (!isTRUE(.ojl_dot(x_best, r1) / .ojl_dot(x_best, x_best) < 100)) {
        low <- log_a
        log_a <- if (up == Inf) log_a + 5 else (log_a + up) * 0.5
      } else {
        break
      }
    }
    tries <- tries + 1L
  }
  if (rr_best < rr_now) {
    x_best <- x_now
  } else {
    rr_best <- rr_now
    x_best <- x_next
  }
  m <- am(r2, alpha, x_next)
  x_now <- m$x
  r2 <- m$r
  rr_now <- m$rr
  start_map <- 3L
  while (maps < maps_limit && (max_time == Inf || proc.time()[["elapsed"]] < max_time) && iter < iter_limit) {
    iter <- iter + 1L
    i_ord <- i_ord %% length(orders) + 1L
    p <- orders[i_ord]
    lpnormsq_best <- rr_best
    i <- start_map
    good_maps <- TRUE
    while (i <= p) {
      m <- am(G(x_now), alpha, x_now)
      if (i == 1L) r1 <- m$r else if (i == 2L) r2 <- m$r else r3 <- m$r
      x_next <- m$x
      rr_now <- m$rr
      maps <- maps + 1L
      if (!isTRUE(rr_now < rr_best * reltol_rr)) {
        good_maps <- FALSE
        break
      } else if (rr_now <= rr_best) {
        rr_best <- rr_now
        lpnormsq_best <- rr_best
        tmp <- x_best
        x_best <- x_now
        x_now <- tmp
        if (lpnormsq_best < abstolsq) break
      }
      tmp <- x_now
      x_now <- x_next
      x_next <- tmp
      i <- i + 1L
    }
    start_map <- 2L
    if (lpnormsq_best < abstolsq) break
    if (good_maps) {
      if (p == 3L) r3 <- r3 + (-2 * r2 + r1)
      r2 <- r2 - r1
      sigma <- -abs(if (p == 2L) .ojl_dot(r2, r1) / .ojl_dot(r2, r2) else .ojl_dot(r3, r2) / .ojl_dot(r3, r3))
    }
    if (.speedmapping_isbad(sigma) || !good_maps) {
      alpha <- alpha * (if (good_maps) 5 else 0.2)
      sigma <- 1
      m <- am(G(x_best), alpha, x_best)
      maps <- maps + 1L
      x_now <- m$x
      r1 <- m$r
      rr_best <- m$rr
    } else {
      # the squared or cubic extrapolation from x_now
      s2 <- sigma * sigma
      if (p == 2L) {
        x_now <- x_now + ((-2 * alpha * (1 + sigma)) * r1 + (alpha * (-1 + s2)) * r2)
      } else {
        x_now <- x_now + ((-3 * alpha * (1 + sigma)) * r1 + (-3 * alpha * (1 - s2)) * r2 +
                            (-alpha * (1 + s2 * sigma)) * r3)
      }
      if (!is.null(lower)) x_now <- pmax(x_now, (1 - buffer) * lower + buffer * x_best)
      if (!is.null(upper)) x_now <- pmin(x_now, (1 - buffer) * upper + buffer * x_best)
      alpha <- alpha * (if (abs(sigma) < 1) 2 / 3 else if (abs(sigma) < 2) 1.01 else 3 / 2)
      damp <- 1
      for (attempt in 1:100) {
        m <- am(G(x_now), alpha, x_now)
        x_next <- m$x
        r1 <- m$r
        rr_now <- m$rr
        maps <- maps + 1L
        if (isTRUE(rr_now < rr_best * reltol_rr)) break
        damp <- damp * 0.5
        x_now <- x_best + damp * (x_now - x_best)
      }
      rr_best <- rr_now
      tmp <- x_best
      x_best <- x_now
      x_now <- x_next
      x_next <- tmp
    }
  }
  status <- if (lpnormsq_best < abstolsq) "first_order" else if (iter_limit - iter <= 0) "max_iter" else
    if (maps_limit - maps <= 0) "max_eval" else if (proc.time()[["elapsed"]] >= max_time) "max_time" else "failure"
  result(x_best, status)
}
