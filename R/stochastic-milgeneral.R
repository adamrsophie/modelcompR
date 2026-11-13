# RKMilGeneral of StochasticDiffEqMilstein (src/algorithms.jl, src/
# alg_utils.jl: _levyarea_terms and _z_prototype, src/caches/, src/
# perform_step/rkmilgeneral.jl: the out of place step with
# _compute_iterated_I). The second noise process carries the MronRoe
# Fourier coefficients of the Levy area of each step.
#
# Deliberate differences from upstream, each in its own helper:
# * upstream unpacks the increments dZ of the second Wiener process as
#   standard normal coefficients, but they have variance dt; .rkmg_scale
#   divides them by sqrt(dt) (for a step, a stored piece or a grid interval);
# * the grid strategy sums only the grid intervals inside [t, t + dt];
#   .rkmg_grid_ends adds the partial intervals at both ends, commutatively;
# * the out of place step for nondiagonal noise resets ggprime_norm for
#   every noise column, so the error estimate sees only the last one;
#   .rkmg_norm_acc sums the norms over the columns, as the in place step;
# * the error estimate for diagonal noise squares (du2 - L) / sqrt(dt), the
#   drift at the predictor minus the diffusion, where RKMil (and the
#   nondiagonal branch) square the Milstein derivative ggprime; that term
#   does not vanish with the step and stalls adaptive solves, so
#   .rkmg_diag_en uses ggprime as the RKMil step of StochasticDiffEqLowOrder
#   (src/perform_step/low_order.jl) does.

#' General Runge-Kutta Milstein method
#'
#' `RKMilGeneral()` is the Runge-Kutta Milstein method of strong order 1 for
#' any noise. For diagonal and scalar noise it uses `J = dW^2 / 2` (minus
#' `dt / 2` in the Ito interpretation); for nondiagonal noise it uses the
#' iterated integrals `J = dW dW' / 2 + dt A` with the Levy area `A` of
#' [MronRoe()], computed from the second noise process, whose increments over
#' a step carry the `coefficient_length(m, n, MronRoe())` Fourier
#' coefficients (scaled by `sqrt(dt)`). With `n` terms given by `p`, or by
#' `terms_needed(m, dt, c dt^(3/2), MronRoe(), MaxL2())`. After a rejected
#' step the iterated integrals are assembled from the stored noise pieces of
#' the step, and with a [NoiseGrid()] from the grid intervals inside the step
#' when there are at least two. The error estimate is that of the
#' Runge-Kutta Milstein method.
#'
#' @param interpretation `"Ito"` or `"Stratonovich"`.
#' @param ii_approx [IILevyArea()] or [IICommutative()], the evaluator of
#'   the iterated integrals when the step has no coefficients.
#' @param c The accuracy constant of the truncation, `eps = c dt^(3/2)`.
#' @param p The number of Fourier terms, `NULL` (or `TRUE`) to choose it.
#' @param dt Accepted for the signature; the truncation comes from the step
#'   size of the solve.
#' @return An SDE algorithm.
#' @examples
#' f <- function(u, p, t) -u
#' g <- function(u, p, t) matrix(c(0.2 * u[1], 0.1, 0.1, 0.2 * u[2]), 2, 2)
#' prob <- SDEProblem(f, g, c(1, 1), c(0, 1), noise_rate_prototype = matrix(0, 2, 2))
#' sol <- solve(prob, RKMilGeneral(p = 5), dt = 0.1, adaptive = FALSE, seed = 1)
#' sol$u[[11]]
#' @export
RKMilGeneral <- function(interpretation = c("Ito", "Stratonovich"), ii_approx = IILevyArea(), c = 1, p = NULL,
                         dt = NULL) {
  interpretation <- match.arg(interpretation)
  if (isTRUE(p)) p <- NULL
  alg <- .new_sde_alg("RKMilGeneral", 1, TRUE, perform = NULL, needs_dZ = TRUE, delta = 1,
                      interpretation = interpretation, ii_approx = ii_approx, c = as.integer(c), p = p)
  # _z_prototype: the length of the second process
  alg$z_length <- function(m, dt) coefficient_length(m, .rkmg_terms(alg, m, dt), MronRoe())
  alg$bind <- function(prob, opts) .rkmg_bind(alg, prob, opts)
  alg
}

# _levyarea_terms
.rkmg_terms <- function(alg, m, dt) {
  if (is.numeric(alg$p)) return(max(1L, as.integer(alg$p)))
  if (is.null(dt) || !is.finite(dt) || dt == 0) return(max(10L, m))
  h <- abs(dt)
  max(1L, terms_needed(m, h, alg$c * h^(1 + 1 / 2), MronRoe(), MaxL2()))
}

.rkmg_scale <- function(dZ, h) dZ / sqrt(abs(h))

.rkmg_norm_acc <- function(acc, x) acc + x

.rkmg_diag_en <- function(dW, ggprime, du2, L, sqdt) dW^3 * ggprime^2 / 6

# _unpack_dZ_to_coefficients after scaling: X (m x n), Y (n x m) and the
# tail ((m^2 + m) / 2), packed one after the other.
.rkmg_coeffs <- function(dZ, m, h) {
  z <- .rkmg_scale(as.numeric(dZ), h)
  ntail <- (m^2 + m) %/% 2
  n <- (length(z) - ntail) %/% (2L * m)
  if (n <= 0L) return(NULL)
  LevyAreaCoefficients(matrix(z[seq_len(m * n)], m, n), matrix(z[m * n + seq_len(n * m)], n, m),
                       z[2L * m * n + seq_len(ntail)], m, n)
}

# The Stratonovich iterated integrals of one piece (h, dW, dZ): with
# coefficients dW dW' / 2 + h A, the commutative part otherwise.
.rkmg_piece_J <- function(h, dW, dZ, m) {
  J <- 0.5 * outer(dW, dW)
  if (!is.null(dZ) && length(dZ) >= 2L * m && h > 0) {
    cf <- .rkmg_coeffs(dZ, m, h)
    if (!is.null(cf)) J <- J + h * levyarea(dW / sqrt(h), cf$n, MronRoe(), cf)
  }
  J
}

# Chen's relation over consecutive pieces: the pieces' own integrals plus
# the cross terms W_cumsum dW_i'.
.rkmg_chain <- function(pieces, m) {
  J <- matrix(0, m, m)
  Wc <- numeric(m)
  for (pc in pieces) {
    J <- J + .rkmg_piece_J(pc$h, pc$dW, pc$dZ, m) + outer(Wc, pc$dW)
    Wc <- Wc + pc$dW
  }
  J
}

# The partial grid intervals at both ends of the step, without Levy area.
.rkmg_grid_ends <- function(first, inner, last) c(first, inner, last)

# _compute_II_from_grid on the grid of a NoiseGrid
.rkmg_grid <- function(pe, t, dt, m) {
  tg <- pe$grid_t
  if (is.null(tg)) return(NULL)
  Wg <- pe$grid_W; Zg <- pe$grid_Z
  te <- t + dt
  tol <- 100 * .Machine$double.eps * max(1, abs(te))
  i0 <- which(tg >= t - tol)[1L]
  i1 <- utils::tail(which(tg <= te + tol), 1L)
  if (is.na(i0) || !length(i1) || i1 - i0 < 2L) return(NULL)
  inner <- lapply((i0 + 1L):i1, function(n)
    list(h = tg[n] - tg[n - 1L], dW = Wg[n, ] - Wg[n - 1L, ], dZ = if (!is.null(Zg)) Zg[n, ] - Zg[n - 1L, ]))
  W0 <- as.numeric(pe$path(t)); W1 <- as.numeric(pe$path(te))
  first <- if (tg[i0] > t + tol) list(list(h = tg[i0] - t, dW = Wg[i0, ] - W0, dZ = NULL))
  last <- if (tg[i1] < te - tol) list(list(h = te - tg[i1], dW = W1 - Wg[i1, ], dZ = NULL))
  .rkmg_chain(.rkmg_grid_ends(first, inner, last), m)
}

# _compute_iterated_I: the stored pieces of the step, the noise grid, the
# step's coefficients; NULL for diagonal noise or without coefficients.
.rkmg_J <- function(dt, dW, dZ, ns, t, m) {
  if (length(ns$pieces) >= 2L) return(.rkmg_chain(ns$pieces, m))
  if (identical(ns$kind, "path")) {
    J <- .rkmg_grid(ns$pe, t, dt, m)
    if (!is.null(J)) return(J)
  }
  if (is.null(dZ) || length(dZ) < 2L * m) return(NULL)
  cf <- .rkmg_coeffs(dZ, m, dt)
  if (is.null(cf)) return(NULL)
  0.5 * outer(dW, dW) + dt * levyarea(dW / sqrt(dt), cf$n, MronRoe(), cf)
}

.rkmg_bind <- function(alg, prob, opts) {
  type <- prob$noise_type %||% prob$f$noise_type %||% "diagonal"
  general <- type == "general"
  pe <- if (!is.null(prob$noise_process)) .np_env(prob$noise_process)
  a <- alg
  a$bind <- NULL
  # the coefficients are only read for nondiagonal noise, and a path
  # without a second process leaves the step to get_iterated_I
  a$needs_dZ <- general && (is.null(pe) || !identical(pe$kind, "path") || !is.null(pe$zpath))
  ito <- alg$interpretation == "Ito"
  a$init <- function(ctx, t, u) list(Jalg = NULL)
  a$step <- function(ctx, g, u, t, dt, dW, type, delta, cache, dZ = NULL, ns = NULL, ...) {
    p <- ctx$p
    f <- ctx$f
    m <- length(dW)
    if (is.null(cache$Jalg)) cache$Jalg <- get_Jalg(dW, dt, prob, alg)
    J <- if (general) .rkmg_J(dt, dW, dZ, ns, t, m)
    if (is.null(J)) J <- get_iterated_I(dt, dW, dZ, cache$Jalg, alg$p, alg$c, alg$order)
    if (ito) J <- if (general) J - diag(abs(dt) / 2, m) else J - abs(dt) / 2
    sqdt <- .sqdt(dt)
    du1 <- f(u, p, t)
    L <- g(u, p, t)
    K <- u + dt * du1
    gnorm <- 0
    if (!general) {
      utilde <- (if (ito) K else u) + L * sqdt
      ggprime <- (g(utilde, p, t) - L) / sqdt
      unew <- K + L * dW + ggprime * J
    } else {
      L <- as.matrix(L)
      mil <- numeric(length(u))
      for (i in seq_len(m)) {
        ggprime <- (as.matrix(g(K + sqdt * L[, i], p, t)) - L) / sqdt
        gnorm <- .rkmg_norm_acc(gnorm, .rms(ggprime))
        mil <- mil + as.numeric(ggprime %*% J[i, ])
      }
      unew <- K + as.numeric(L %*% dW) + mil
    }
    du2 <- f(K, p, t + dt)
    err <- if (!general) {
      delta * (dt * (du2 - du1) / 2) + .rkmg_diag_en(dW, ggprime, du2, L, sqdt)
    } else {
      rep_len(delta * .rms(dt * (du2 - du1) / 2) + .rms(dW)^3 * gnorm^2 / 6, length(u))
    }
    list(u = unew, err = err, k = NULL, cache = cache, nf = 2L + m)
  }
  a
}
