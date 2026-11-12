# Random ordinary differential equations, after SciMLBase (RODEProblem,
# RODEFunction, RODESolution) and StochasticDiffEqRODE (RandomEM, RandomHeun,
# RandomTamedEM, RandomTaylor15; perform_step/low_order.jl).

#' Random ordinary differential equations
#'
#' `RODEProblem(f, u0, tspan, p)` is `du/dt = f(u, p, t, W)` driven by a
#' noise process `W` (a [WienerProcess()] with the shape of
#' `rand_prototype`, or of `u0`, unless `noise` gives one; see
#' [noise_processes]). `RODEFunction(f)` wraps the right-hand side.
#' `solve(prob, alg, dt = )` integrates with fixed steps: `RandomEM()` is
#' the Euler method with the noise held at its value at the start of the
#' step, `RandomHeun()` the Heun method, `RandomTamedEM()` the tamed Euler
#' method `u + dt f / (1 + dt |f|)`, and `RandomTaylor15()` the derivative
#' free order 1.5 Taylor method of Asai, for a scalar [NoiseGrid()] that
#' resolves the path more finely than the steps. The solution is a
#' `RODESolution`; `sol$W` is the noise path.
#'
#' @param f Right-hand side `f(u, p, t, W)`, or a `RODEFunction`.
#' @param u0,tspan,p Initial state, time span and parameters.
#' @param rand_prototype Shape of the noise when `noise` is not given.
#' @param noise A noise process.
#' @param seed Seed for the noise.
#' @param ... Problem options passed to `solve()`.
#' @return A problem, a function, an algorithm or a solution.
#' @examples
#' f <- function(u, p, t, W) 1.01 * u + 0.87 * u * W
#' prob <- RODEProblem(f, 0.5, c(0, 1))
#' sol <- solve(prob, RandomEM(), dt = 0.01, seed = 1)
#' tail(sol$u, 1)
#' @name rode_problems
NULL

#' @rdname rode_problems
#' @export
RODEFunction <- function(f, ...) structure(list(f = f, ...), class = c("RODEFunction", "SciMLFunction"))

#' @rdname rode_problems
#' @export
RODEProblem <- function(f, u0, tspan, p = numeric(), rand_prototype = NULL, noise = NULL, seed = 0L, ...) {
  fn <- if (inherits(f, "RODEFunction")) f else RODEFunction(f)
  u0 <- as.numeric(unlist(u0))
  prob <- .new_problem("RODEProblem", f = fn, u0 = u0, tspan = .as_tspan(tspan), p = p, sys = NULL,
                       kwargs = list(...))
  prob$rand_prototype <- if (!is.null(rand_prototype)) as.numeric(rand_prototype) else 0 * u0
  prob$noise <- noise
  prob$seed <- seed
  prob
}

.rode_alg <- function(name, order, perform, ...)
  structure(list(name = name, order = order, perform = perform, ...), class = c(name, "RODEAlgorithm"))

#' @export
print.RODEAlgorithm <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

#' @rdname rode_problems
#' @export
RandomEM <- function() .rode_alg("RandomEM", 0.5, function(f, u, p, t, dt, W0, dW, ns) u + dt * f(u, p, t, W0))

#' @rdname rode_problems
#' @export
RandomTamedEM <- function() .rode_alg("RandomTamedEM", 0.5, function(f, u, p, t, dt, W0, dW, ns) {
  ft <- f(u, p, t, W0)
  u + dt * ft / (1 + dt * sqrt(sum(ft^2)))
})

#' @rdname rode_problems
#' @export
RandomHeun <- function() .rode_alg("RandomHeun", 0.5, function(f, u, p, t, dt, W0, dW, ns) {
  ft <- f(u, p, t, W0)
  tmp <- u + dt * ft
  u + (dt / 2) * (ft + f(tmp, p, t + dt, W0 + dW))
})

# path_integrals: the integrals of v = W - w0 and v^2 over the step from the
# saved path (trapezoid and Simpson type rules between saved points).
.rode_path_integrals <- function(ts, Ws, t, dt, w0, dW) {
  tend <- t + dt
  I1 <- 0; I2 <- 0; tprev <- t; vprev <- 0
  for (i in which(ts > t & ts <= tend)) {
    vi <- Ws[i] - w0
    I1 <- I1 + (ts[i] - tprev) * (vprev + vi) / 2
    I2 <- I2 + (ts[i] - tprev) * (vprev^2 + vprev * vi + vi^2) / 3
    tprev <- ts[i]; vprev <- vi
  }
  c(I1 + (tend - tprev) * (vprev + dW) / 2, I2 + (tend - tprev) * (vprev^2 + vprev * dW + dW^2) / 3)
}

#' @rdname rode_problems
#' @export
RandomTaylor15 <- function() .rode_alg("RandomTaylor15", 1.5, function(f, u, p, t, dt, W0, dW, ns) {
  adt <- abs(dt); h <- sqrt(adt)
  I <- .rode_path_integrals(ns$grid_t, ns$grid_W, t, dt, W0, dW)
  ft <- f(u, p, t, W0)
  ut <- u + dt * ft
  f0 <- f(ut, p, t + dt, W0)
  fp <- f(ut, p, t + dt, W0 + h)
  fm <- f(ut, p, t + dt, W0 - h)
  u + (dt / 2) * (ft + f0) + (I[1L] / (2 * h)) * (fp - fm) + (I[2L] / (2 * adt)) * (fp - 2 * f0 + fm)
}, needs_grid = TRUE)

#' @export
solve.RODEProblem <- function(a, b = RandomEM(), ...) {
  opts <- utils::modifyList(a$kwargs %||% list(), list(...))
  if (!inherits(b, "RODEAlgorithm")) stop("not a RODE algorithm: ", class(b)[1L], call. = FALSE)
  if (is.null(opts$dt)) stop("RODE solvers take fixed steps: pass dt", call. = FALSE)
  seed <- opts$seed %||% (if (!identical(a$seed, 0L) && !identical(a$seed, 0)) a$seed)
  opts$seed <- NULL
  noise <- a$noise %||% WienerProcess(a$tspan[1L], a$rand_prototype)
  if (isTRUE(b$needs_grid) && !(inherits(noise, "NoiseGrid") && ncol(.np_env(noise)$grid_W) == 1L))
    stop("RandomTaylor15 needs a scalar NoiseGrid", call. = FALSE)
  .with_seed(seed, {
    proc <- .np_copy(noise)
    pe <- .np_env(proc)
    if (isTRUE(pe$reset)) .np_reinit(pe, a$tspan[1L])
    ns <- .noise_new(a$tspan[1L], length(pe$W0), FALSE, proc = proc)
    ff <- a$f$f
    grid <- if (isTRUE(b$needs_grid)) list(grid_t = pe$grid_t, grid_W = as.numeric(pe$grid_W))
    runner <- .new_alg(b$name, b$order, adaptive = FALSE,
      step = function(ctx, t, u, dt, cache) {
        inc <- .noise_take(ns, t, dt, u, ctx$p)
        unew <- b$perform(ff, u, ctx$p, t, dt, inc$W0, inc$dW, grid)
        list(u = unew, err = NULL, k = NULL, cache = cache, nf = 1L)
      },
      interp = .sde_linear_interp, init = function(ctx, t, u) list())
    dprob <- a
    dprob$f <- list(f = function(u, p, t) 0 * u)
    opts$adaptive <- FALSE
    sol <- .integrate(dprob, runner, opts)
    env <- .sol_env(sol)
    .noise_settle(ns, env$t[length(env$t)])
    env$W <- list(t = ns$log_t, W = do.call(rbind, ns$log_W))
    env$alg <- b
    env$prob <- a
    .as_rode_solution(sol)
  })
}
