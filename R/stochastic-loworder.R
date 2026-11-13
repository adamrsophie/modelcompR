# Further low order SDE methods of StochasticDiffEqLowOrder (perform_step/
# low_order.jl, split.jl, predcorr.jl) and StochasticDiffEqMilstein
# (perform_step/explicit_3s_mil_methods.jl).

#' Further low order SDE methods
#'
#' * `SimplifiedEM()`: the simplified Euler-Maruyama method, weak order 1,
#'   with two point random variables in place of the Brownian increments.
#' * `SplitEM()`: Euler-Maruyama for a [SplitSDEProblem()], the drift being
#'   the sum of its two parts.
#' * `PCEuler(ggprime, theta, eta)`: the drift implicit predictor corrector
#'   Euler method of Kloeden and Platen; `ggprime(u, p, t)` is the Ito
#'   correction \eqn{\sum_k g_k' g_k} (for scalar noise \eqn{g g'}).
#' * `WangLi3SMil_A()` to `WangLi3SMil_F()`: Wang and Li's three stage
#'   Milstein methods (strong order 1) for diagonal and scalar noise. Each
#'   composes a drift step, a half step correction with the derivative free
#'   diffusion difference and a Milstein noise step, in six orders.
#'
#' @param ggprime The function `ggprime(u, p, t)` of `PCEuler()`.
#' @param theta,eta Implicitness of the drift and of the diffusion in
#'   `PCEuler()`.
#' @return An SDE algorithm.
#' @examples
#' prob <- SDEProblem(function(u, p, t) 1.01 * u, function(u, p, t) 0.87 * u, 0.5, c(0, 1))
#' sol <- solve(prob, WangLi3SMil_A(), dt = 0.01, seed = 1)
#' pc <- PCEuler(function(u, p, t) 0.87^2 * u)
#' sol2 <- solve(prob, pc, dt = 0.01, seed = 1)
#' @name sde_loworder_methods
NULL

#' @rdname sde_loworder_methods
#' @export
SimplifiedEM <- function() {
  .new_sde_alg("SimplifiedEM", 0.5, FALSE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    K <- u + dt * f(u, p, t)
    list(u = K + .gdw(g(u, p, t), calc_twopoint_random(.sqdt(dt), dW), type), nf = 1L)
  }, weak = TRUE)
}

#' @rdname sde_loworder_methods
#' @export
SplitEM <- function() {
  .new_sde_alg("SplitEM", 0.5, FALSE, function(f, g, u, p, t, dt, dW, dZ, type, delta)
    list(u = dt * f(u, p, t) + .gdw(g(u, p, t), dW, type) + u, nf = 1L))
}

#' @rdname sde_loworder_methods
#' @export
PCEuler <- function(ggprime, theta = 1 / 2, eta = 1 / 2) {
  if (missing(ggprime) || !is.function(ggprime))
    stop("PCEuler needs ggprime(u, p, t), the Ito correction sum_k g_k' g_k", call. = FALSE)
  .new_sde_alg("PCEuler", 0.5, FALSE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    gdW_n <- .gdw(g(u, p, t), dW, type)
    f_n <- f(u, p, t)
    ubar <- u + dt * f_n + gdW_n
    fbar_n <- f_n - eta * ggprime(u, p, t)
    tnp1 <- t + dt
    fbar_np1 <- f(ubar, p, tnp1) - eta * ggprime(ubar, p, tnp1)
    gdW_np1 <- .gdw(g(ubar, p, tnp1), dW, type)
    list(u = u + (theta * dt) * fbar_np1 + ((1 - theta) * dt) * fbar_n + eta * gdW_np1 + (1 - eta) * gdW_n,
         nf = 2L)
  }, theta = theta, eta = eta, ggprime = ggprime)
}

# One method of Wang and Li: D is the drift step, C the correction
# u - dt/2 (g(u + g sqdt) - g) / sqdt and N the noise step
# u + g dW + dW^2/2 (g(u + g sqdt) - g) / sqdt, all at time t.
.wangli_alg <- function(name, stages) {
  .new_sde_alg(name, 1, FALSE, function(f, g, u, p, t, dt, dW, dZ, type, delta) {
    sqdt <- .sqdt(dt)
    # upstream divides by sqdt in the first diffusion stage and by
    # abs(sqdt) in the second (explicit_3s_mil_methods.jl)
    div <- c(sqdt, abs(sqdt))
    ng <- 0L
    for (s in stages) {
      if (s == "D") {
        u <- u + dt * f(u, p, t)
        next
      }
      ng <- ng + 1L
      k <- as.numeric(g(u, p, t))
      k1 <- (as.numeric(g(u + k * sqdt, p, t)) - k) / div[ng]
      u <- if (s == "C") u - (dt / 2) * k1 else u + k * dW + 0.5 * dW^2 * k1
    }
    list(u = u, nf = 1L)
  }, noise = c("diagonal", "scalar"))
}

#' @rdname sde_loworder_methods
#' @export
WangLi3SMil_A <- function() .wangli_alg("WangLi3SMil_A", c("D", "C", "N"))
#' @rdname sde_loworder_methods
#' @export
WangLi3SMil_B <- function() .wangli_alg("WangLi3SMil_B", c("N", "D", "C"))
#' @rdname sde_loworder_methods
#' @export
WangLi3SMil_C <- function() .wangli_alg("WangLi3SMil_C", c("D", "N", "C"))
#' @rdname sde_loworder_methods
#' @export
WangLi3SMil_D <- function() .wangli_alg("WangLi3SMil_D", c("C", "D", "N"))
#' @rdname sde_loworder_methods
#' @export
WangLi3SMil_E <- function() .wangli_alg("WangLi3SMil_E", c("C", "N", "D"))
#' @rdname sde_loworder_methods
#' @export
WangLi3SMil_F <- function() .wangli_alg("WangLi3SMil_F", c("N", "C", "D"))
