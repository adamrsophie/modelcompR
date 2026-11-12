# The BAOAB Langevin integrator of StochasticDiffEqRODE (src/algorithms.jl,
# src/alg_utils.jl, src/caches/dynamical_caches.jl,
# src/perform_step/dynamical.jl), for a DynamicalSDEProblem with the
# velocities v first in the state c(v, u).

#' BAOAB Langevin dynamics integrator
#'
#' `BAOAB(gamma, scale_noise)` integrates the Langevin dynamics of a
#' [DynamicalSDEProblem()] with `dv = f1(v, u, p, t) dt - gamma v dt +
#' sigma dW` and `du = v dt` by the splitting kick (B), drift (A),
#' Ornstein Uhlenbeck (O), drift (A), kick (B) of Leimkuhler and Matthews
#' (2013). The O part damps the velocity by `c1 = exp(-gamma dt)` and adds
#' `c2 g dW` with `c2 = sqrt((1 - c1^2) / dt)` (or `1` with `scale_noise =
#' FALSE`), both fixed by the initial step size. The method needs `f2(v, u,
#' p, t) = v`, diagonal noise on the velocities, evaluated as `g(v, u, p,
#' t)` at the drifted position and half kicked velocity, and a fixed step
#' size `dt`. It has strong order 1.
#'
#' @param gamma The friction coefficient.
#' @param scale_noise Scale the noise by `sqrt((1 - exp(-2 gamma dt)) / dt)`.
#' @return An SDE algorithm.
#' @references Leimkuhler B., Matthews C. (2013). Robust and efficient
#'   configurational molecular sampling via Langevin dynamics. Journal of
#'   Chemical Physics 138, 174102.
#' @examples
#' prob <- DynamicalSDEProblem(function(v, u, p, t) -u, function(v, u, p, t) v,
#'                             function(v, u, p, t) 1, 0.5, 1, c(0, 1))
#' sol <- solve(prob, BAOAB(gamma = 0.5), dt = 0.1, seed = 1)
#' sol$u[[11]]
#' @export
BAOAB <- function(gamma = 1, scale_noise = TRUE) {
  alg <- .new_sde_alg("BAOAB", 1, FALSE, NULL, noise = "diagonal", gamma = gamma, scale_noise = scale_noise)
  alg$bind <- function(prob, opts) .baoab_bind(alg, prob, opts)
  alg
}

.baoab_bind <- function(alg, prob, opts) {
  if (!inherits(prob, "DynamicalSDEProblem"))
    stop("BAOAB() integrates a DynamicalSDEProblem", call. = FALSE)
  nv <- prob$nv
  f1 <- prob$f$f1; f2 <- prob$f$f2; gd <- prob$f$g_dyn
  iv <- seq_len(nv)
  a <- alg
  a$bind <- NULL
  a$noise_dim <- nv
  # initialize!: verify_f2 and the first kick derivative
  a$init <- function(ctx, t, u) {
    v <- u[iv]; q <- u[-iv]
    if (any(as.numeric(f2(v, q, ctx$p, t)) != v))
      stop("Algorithm BAOAB is not applicable if f2(p, q, t) != p", call. = FALSE)
    list(k = as.numeric(f1(v, q, ctx$p, t)), c = NULL)
  }
  a$step <- function(ctx, g, u, t, dt, dW, type, delta, cache, ...) {
    p <- ctx$p
    if (is.null(cache$c)) {
      # alg_cache: c1 and c2 from the step size at the start
      dt0 <- opts$dt %||% dt
      c1 <- exp(-alg$gamma * dt0)
      cache$c <- c(c1, if (alg$scale_noise) sqrt((1 - c1^2) / abs(dt0)) else 1)
    }
    half <- 1 / 2
    du1 <- u[iv]; u1 <- u[-iv]
    du2 <- du1 + half * dt * cache$k
    u2 <- u1 + half * dt * du2
    noise <- rep_len(as.numeric(gd(du2, u2, p, t + dt * half)), nv) * dW
    du3 <- cache$c[1L] * du2 + cache$c[2L] * noise
    unew <- u2 + half * dt * du3
    cache$k <- as.numeric(f1(du3, unew, p, t + dt))
    list(u = c(du3 + half * dt * cache$k, unew), err = NULL, k = NULL, cache = cache, nf = 2L)
  }
  a
}
