# OptimizationPolyalgorithms (lib/OptimizationPolyalgorithms/src, v0.3.12) and
# OptimizationSophia (lib/OptimizationSophia/src, v1.4.2).

#' Polyalgorithm and Sophia
#'
#' `PolyOpt()` chains two optimizers. For a deterministic objective it runs
#' `Adam(0.01)` for `maxiters` iterations (300 by default) and continues
#' from its result with `Optim$BFGS(initial_stepnorm = 0.01)` (for
#' `maxiters` iterations as well, when given). An objective that gives
#' different values at the same point needs `maxiters` and is minimized with
#' `Adam(0.1)`. The solution and its statistics are those of the last
#' optimizer. As in Julia, problems with bounds are not accepted.
#'
#' `Sophia(eta, betas, epsilon, lambda, k, rho)` is the second order
#' optimizer of Liu et al. (2023): a moving average `m` of the gradients, a
#' moving average `h` of the diagonal of the Hessian estimated as
#' `v * (H v)` with a standard normal `v` (Hutchinson's estimator), weight
#' decay `x - eta * lambda * x` and the step `-eta * clamp(m / max(h,
#' epsilon), -rho, rho)`. `maxiters` defaults to 1000. As in Julia, where
#' the Hessian estimate is made when the index of the data batch `i`
#' satisfies `i %% k == 1`, it is made at every iteration for `k >= 2` (the
#' problem has a single batch) and never for `k = 1`. The product `H v`
#' uses the Hessian of the problem when it has one, else central finite
#' differences of the gradient; the random vectors come from [rnorm()]. The
#' objective of the solution is the one evaluated before the last update and
#' the return code is `"Success"`; the callback is called at every
#' iteration and stops the optimization when it returns `TRUE`.
#'
#' @param eta Learning rate.
#' @param betas Decay of the moving averages of the gradient and of the
#'   Hessian diagonal.
#' @param epsilon Lower bound of the Hessian diagonal in the step.
#' @param lambda Weight decay.
#' @param k Frequency of the Hessian estimates (see above).
#' @param rho Bound of the clipped step.
#' @return An optimization algorithm.
#' @examples
#' f <- function(u, p) (1 - u[1])^2 + 100 * (u[2] - u[1]^2)^2
#' g <- function(u, p) c(-2 * (1 - u[1]) - 400 * u[1] * (u[2] - u[1]^2), 200 * (u[2] - u[1]^2))
#' prob <- OptimizationProblem(f, c(-1.2, 1), grad = g)
#' solve(prob, PolyOpt())$u
#' set.seed(1)
#' solve(prob, Sophia(eta = 0.01), maxiters = 200)$objective
#' @name polyopt_sophia
#' @export
PolyOpt <- function() {
  a <- .opt_alg("PolyOpt")
  a$solver <- .polyopt_solve
  a$traits <- list(grad = TRUE)
  a
}

.polyopt_solve <- function(prob, alg, opts) {
  f <- function(x) as.numeric(prob$f$f(x, prob$p))[1L]
  deterministic <- isTRUE(f(prob$u0) == f(prob$u0))
  maxiters <- opts$maxiters
  # SciMLBase.allowsbounds is false for PolyOpt, so the bounded branch of its
  # __solve (BFGS alone) is never reached
  if (!is.null(prob$lb) || !is.null(prob$ub))
    stop(IncompatibleOptimizerError(paste0("The algorithm PolyOpt does not support box constraints. Either remove ",
                                           "the `lb` or `ub` bounds passed to `OptimizationProblem` or use a ",
                                           "different algorithm.")))
  if (!deterministic && is.null(maxiters))
    stop("Automatic optimizer determination requires deterministic loss functions (and no data) or maxiters ",
         "must be specified.", call. = FALSE)
  if (!deterministic) return(.opt_solve_core(prob, Adam(0.1), opts))
  res1 <- .opt_solve_core(prob, Adam(0.01), utils::modifyList(opts, list(maxiters = maxiters %||% 300L)))
  prob2 <- prob
  prob2$u0 <- res1$u
  .opt_solve_core(prob2, Optim$BFGS(initial_stepnorm = 0.01), opts)
}

#' @rdname polyopt_sophia
#' @export
Sophia <- function(eta = 1e-3, betas = c(0.9, 0.999), epsilon = 1e-8, lambda = 1e-1, k = 10L, rho = 0.04) {
  if (length(betas) != 2L) stop("`betas` must have two components", call. = FALSE)
  a <- .opt_alg("Sophia", eta = eta, betas = as.numeric(betas), epsilon = epsilon, lambda = lambda,
                k = as.integer(k), rho = rho)
  a$solver <- .sophia_solve
  a$traits <- list(grad = TRUE, hess = TRUE)
  a
}

.sophia_solve <- function(prob, alg, opts) {
  fs <- .opt_fns(prob)
  if (!is.null(prob$lb) || !is.null(prob$ub))
    stop(IncompatibleOptimizerError(paste0("The algorithm Sophia does not support box constraints. Either remove ",
                                           "the `lb` or `ub` bounds passed to `OptimizationProblem` or use a ",
                                           "different algorithm.")))
  if (is.function(prob$f$cons))
    stop(IncompatibleOptimizerError(paste0("The algorithm Sophia does not support constraints. Either remove ",
                                           "the `cons` function passed to `OptimizationFunction` or use a ",
                                           "different algorithm.")))
  maxiters <- .opt_maxiters(opts$maxiters %||% 1000L)
  n <- length(prob$u0)
  hv <- if (is.function(fs$fn$hess)) {
    function(x, v) as.numeric(matrix(as.numeric(fs$fn$hess(x, fs$p)), n, n) %*% v)
  } else {
    function(x, v) {
      h <- .Machine$double.eps^(1 / 3) * max(1, sqrt(sum(x * x))) / max(sqrt(sum(v * v)), .Machine$double.xmin)
      (fs$g(x + h * v) - fs$g(x - h * v)) / (2 * h)
    }
  }
  b1 <- alg$betas[1L]
  b2 <- alg$betas[2L]
  eta <- alg$eta
  theta <- as.numeric(prob$u0)
  m <- h <- numeric(n)
  x <- fs$f(theta)
  t0 <- proc.time()[["elapsed"]]
  iters <- 0L
  for (epoch in seq_len(maxiters)) {
    g <- fs$g(theta)
    x <- fs$f(theta)
    iters <- epoch
    m <- b1 * m + (1 - b1) * g
    # the index of the data batch is 1: an estimate at every iteration
    # unless k is 1
    if (1L %% alg$k == 1L) {
      v <- stats::rnorm(n)
      h <- b2 * h + (1 - b2) * (v * hv(theta, v))
    }
    theta <- theta - eta * alg$lambda * theta
    theta <- theta - eta * pmax(pmin(m / pmax(h, alg$epsilon), alg$rho), -alg$rho)
  }
  sol <- .opt_solution(prob, alg, theta, x, "Success", NULL)
  .opt_set_fields(sol, stats = OptimizationStats(iterations = iters, time = proc.time()[["elapsed"]] - t0))
}
