# Iterated stochastic integrals and Levy areas, after
# StochasticDiffEqLevyArea.jl (src/fourier.jl, milstein.jl, wiktorsson.jl,
# mronroe.jl: the four Fourier series algorithms and their error
# coefficients; error_norms.jl: MaxL2, FrobeniusL2 and the internal LpMax;
# coefficients.jl and generate.jl: LevyAreaCoefficients; alg_utils.jl:
# terms_needed, optimal_algorithm; iterated_integrals.jl; reconstruct.jl:
# reconstruct_path, iterated_integrals_subinterval) and the iterated integral
# evaluators of StochasticDiffEqCore (src/iterated_integrals.jl: JDiagonal,
# JCommute, IteratedIntegralAlgorithm_iip, get_iterated_I, get_Jalg;
# src/algorithms.jl: IteratedIntegralApprox, IICommutative, IILevyArea).
#
# Each algorithm draws its random numbers in the order of
# generate_coefficients, so levyarea() without coefficients is
# generate_coefficients() followed by the deterministic kernel. The `rng`
# argument is NULL (R's global generator) or a function of a count n that
# returns n standard normal numbers, which lets a caller feed the numbers of
# another generator.

#' Levy areas and iterated stochastic integrals
#'
#' The Levy area algorithms of StochasticDiffEqLevyArea: `Fourier()` (the
#' plain truncated Fourier series, error order 1/2 in the number of terms),
#' `Milstein()` (with Milstein's Gaussian tail approximation, order 1/2),
#' `Wiktorsson()` (Wiktorsson's tail sum correction, order 1) and `MronRoe()`
#' (the Mrongowius and Roessler correction, order 1). `MaxL2()` and
#' `FrobeniusL2()` are the error norms used to choose the number of terms.
#'
#' `levyarea(W, n, alg)` approximates the antisymmetric Levy area matrix of
#' the normalized Wiener increment `W` (the increment over the step divided by
#' the square root of the step) with `n` Fourier terms. The random numbers
#' come from `coeffs`, a [LevyAreaCoefficients()] object, when given (a
#' deterministic computation), and are otherwise drawn as
#' `generate_coefficients()` draws them.
#'
#' `iterated_integrals(W, h, eps)` simulates the matrix of iterated integrals
#' of the increment `W` over a step `h` with error at most `eps` (default
#' `h^(3/2)`) in the norm `error_norm`, with the algorithm of fewest random
#' numbers unless `alg` is given; `ito_correction = TRUE` gives Ito
#' integrals, `FALSE` Stratonovich ones. Passing coefficients as `eps` uses
#' them (with `MronRoe()` by default); `q_12` gives the square roots of the
#' covariance eigenvalues of a Q-Wiener process. A scalar `W` gives the exact
#' integral `W^2 / 2` (minus `h / 2` for Ito).
#'
#' `terms_needed(dim, stepsize, eps, alg, norm)` is the number of Fourier
#' terms for an error at most `eps`, `optimal_algorithm(dim, stepsize, eps,
#' norm)` the algorithm needing the fewest random numbers and
#' `coefficient_length(m, n, alg)` the number of random numbers of an
#' algorithm. `reconstruct_path(dW, h, coeffs, t_points)` evaluates the
#' Brownian path of the Fourier expansion at times inside the step and
#' `iterated_integrals_subinterval(dW, h, coeffs, t_start, t_end)` the
#' iterated integrals over a part of the step, by a Riemann sum over
#' `n_quadrature` points of the reconstructed path.
#'
#' `AbstractIteratedIntegralAlgorithm` and `AbstractErrorNorm` are the class
#' names of the algorithms and error norms.
#'
#' @param W A Wiener increment (normalized for `levyarea()`).
#' @param n Number of Fourier terms.
#' @param alg A Levy area algorithm.
#' @param coeffs A [LevyAreaCoefficients()] object.
#' @param rng `NULL` (R's generator) or a function of `n` returning `n`
#'   standard normal numbers.
#' @param h,stepsize Step size.
#' @param eps Error bound, or a [LevyAreaCoefficients()] object.
#' @param ito_correction Ito (`TRUE`) or Stratonovich (`FALSE`) integrals.
#' @param error_norm,norm An error norm.
#' @param q_12 Square roots of the covariance eigenvalues of a Q-Wiener
#'   process.
#' @param dim,m Dimension of the Wiener process.
#' @return An algorithm, an error norm, a matrix, a count or a list of
#'   vectors.
#' @examples
#' W <- c(0.3, -0.7, 1.1)
#' cf <- generate_coefficients(3, 4, MronRoe())
#' levyarea(W, 4, MronRoe(), cf)
#' iterated_integrals(W * sqrt(0.01), 0.01, cf)
#' terms_needed(3, 0.01, 0.01^1.5, Fourier(), MaxL2())
#' optimal_algorithm(3, 0.01)
#' @name levy_area
NULL

.ii_alg <- function(name, convorder, errc, norv)
  structure(list(name = name, convorder = convorder, errc = errc, norv = norv),
            class = c(name, "AbstractIteratedIntegralAlgorithm"))

#' @export
print.AbstractIteratedIntegralAlgorithm <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

#' @rdname levy_area
#' @export
AbstractIteratedIntegralAlgorithm <- "AbstractIteratedIntegralAlgorithm"

#' @rdname levy_area
#' @export
AbstractErrorNorm <- "AbstractErrorNorm"

#' @rdname levy_area
#' @export
Fourier <- function()
  .ii_alg("Fourier", 1 / 2, function(m, h) h * sqrt(3) / (sqrt(2) * pi), function(m, n) 2 * m * n)

#' @rdname levy_area
#' @export
Milstein <- function()
  .ii_alg("Milstein", 1 / 2, function(m, h) h / (sqrt(2) * pi), function(m, n) 2 * m * n + m)

#' @rdname levy_area
#' @export
Wiktorsson <- function()
  .ii_alg("Wiktorsson", 1, function(m, h) sqrt(5 * m) * h / (sqrt(12) * pi),
          function(m, n) 2 * m * n + (m^2 - m) %/% 2)

#' @rdname levy_area
#' @export
MronRoe <- function()
  .ii_alg("MronRoe", 1, function(m, h) sqrt(m) * h / (sqrt(12) * pi),
          function(m, n) 2 * m * n + (m^2 + m) %/% 2)

.ii_algs <- function() list(Fourier(), Milstein(), Wiktorsson(), MronRoe())

.err_norm_obj <- function(name, p, q = NULL)
  structure(list(name = name, p = p, q = q), class = c(name, "AbstractErrorNorm"))

#' @export
print.AbstractErrorNorm <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

#' @rdname levy_area
#' @export
MaxL2 <- function() .err_norm_obj("MaxL2", 2)

#' @rdname levy_area
#' @export
FrobeniusL2 <- function() .err_norm_obj("FrobeniusL2", 2, 2)

# LpMax{2}, the norm of the largest entry error (internal in Julia too)
.LpMax <- function() .err_norm_obj("LpMax", 2)

# errcoeff (error_norms.jl): the constant of the error bound for the norm.
.la_errcoeff <- function(m, h, alg, norm, q_12 = NULL) {
  if (norm$p != 2) stop("only the L2 error norms are implemented", call. = FALSE)
  base <- alg$errc(m, h)
  switch(norm$name,
    MaxL2 = if (is.null(q_12)) base else {
      prods <- unlist(lapply(seq_len(m), function(i) q_12[i] * q_12[seq_len(i - 1L)]))
      max(prods) * base
    },
    FrobeniusL2 = if (is.null(q_12)) (m^2 - m)^(1 / norm$q) * base else {
      qq <- norm$q
      (sum(q_12^qq)^2 - sum(q_12^(2 * qq)))^(1 / qq) * base
    },
    LpMax = 1 / 2^(1 / 2) * .la_errcoeff(m, h, alg, FrobeniusL2(), q_12),
    stop("unknown error norm", call. = FALSE))
}

#' @rdname levy_area
#' @export
terms_needed <- function(dim, stepsize, eps, alg, norm, q_12 = NULL) {
  if (!is.null(q_12) && length(q_12) != dim)
    stop("Length of q_12 must be equal to the dimension.", call. = FALSE)
  as.integer(ceiling((.la_errcoeff(dim, stepsize, alg, norm, q_12) / eps)^(1 / alg$convorder)))
}

.effective_cost <- function(dim, stepsize, eps, alg, norm, q_12 = NULL)
  alg$norv(dim, terms_needed(dim, stepsize, eps, alg, norm, q_12))

#' @rdname levy_area
#' @export
optimal_algorithm <- function(dim, stepsize, eps = stepsize^(3 / 2), norm = NULL, q_12 = NULL) {
  norm <- norm %||% if (is.null(q_12)) MaxL2() else FrobeniusL2()
  algs <- .ii_algs()
  algs[[which.min(vapply(algs, function(a) .effective_cost(dim, stepsize, eps, a, norm, q_12), 0))]]
}

#' @rdname levy_area
#' @export
coefficient_length <- function(m, n, alg) alg$norv(m, n)

#' @rdname levy_area
#' @param X,Y,tail The Fourier coefficient matrices and the tail random
#'   numbers.
#' @export
LevyAreaCoefficients <- function(X, Y, tail, m, n)
  structure(list(X = as.matrix(X), Y = as.matrix(Y), tail = .as_num(tail), m = as.integer(m), n = as.integer(n)),
            class = "LevyAreaCoefficients")

#' @export
print.LevyAreaCoefficients <- function(x, ...) {
  cat("LevyAreaCoefficients (m = ", x$m, ", n = ", x$n, ")\n", sep = "")
  invisible(x)
}

.la_randn <- function(rng, k) {
  if (is.null(rng)) return(stats::rnorm(k))
  z <- .as_num(rng(k))
  if (length(z) != k) stop("rng(n) must return n numbers", call. = FALSE)
  z
}

#' @rdname levy_area
#' @export
generate_coefficients <- function(m, n, alg, rng = NULL) {
  m <- as.integer(m); n <- as.integer(n)
  if (alg$name == "MronRoe") {
    X <- matrix(.la_randn(rng, m * n), m, n)
    Y <- matrix(.la_randn(rng, n * m), n, m)
    tail <- c(.la_randn(rng, m), .la_randn(rng, (m^2 - m) %/% 2))
  } else {
    X <- matrix(.la_randn(rng, n * m), n, m)
    Y <- matrix(.la_randn(rng, m * n), m, n)
    tail <- switch(alg$name, Fourier = numeric(), Milstein = .la_randn(rng, m),
                   Wiktorsson = .la_randn(rng, (m^2 - m) %/% 2), stop("unknown Levy area algorithm", call. = FALSE))
  }
  LevyAreaCoefficients(X, Y, tail, m, n)
}

# The antisymmetrization of the Fourier, Milstein and MronRoe kernels: the
# upper triangle (A[i, j] + extra - A[j, i]) / (2 pi), row by row.
.la_antisym <- function(A, m, extra = NULL) {
  idx <- 0L
  for (i in seq_len(m)) {
    A[i, i] <- 0
    j <- i + 1L
    while (j <= m) {
      idx <- idx + 1L
      A[i, j] <- if (is.null(extra)) (A[i, j] - A[j, i]) / (2 * pi) else (A[i, j] + extra[idx] - A[j, i]) / (2 * pi)
      A[j, i] <- -A[i, j]
      j <- j + 1L
    }
  }
  A
}

.la_kernel <- function(W, n, alg, cf) {
  m <- length(W)
  W <- .as_num(W)
  sq2 <- sqrt(2)
  if (alg$name == "MronRoe") {
    Y <- sweep(cf$Y - sq2 * matrix(W, n, m, byrow = TRUE), 1L, seq_len(n), "/")
    A <- cf$X %*% Y
    a <- sqrt(2 * trigamma(n + 1))
    A <- A + outer(a * W, cf$tail[seq_len(m)])
    return(.la_antisym(A, m, a * cf$tail[-seq_len(m)]))
  }
  Y <- sweep(cf$Y - sq2 * W, 2L, seq_len(n), "/")
  A <- Y %*% cf$X
  if (alg$name == "Fourier") return(.la_antisym(A, m))
  a <- sqrt(2 * trigamma(n + 1))
  if (alg$name == "Milstein") return(.la_antisym(A + outer(a * W, cf$tail[seq_len(m)]), m))
  # Wiktorsson
  G <- matrix(0, m, m)
  idx <- 0L
  for (j in seq_len(m)) {
    i <- j + 1L
    while (i <= m) {
      idx <- idx + 1L
      g <- a * cf$tail[idx]
      G[i, j] <- g; G[j, i] <- -g
      A[i, j] <- A[i, j] + g
      i <- i + 1L
    }
  }
  A <- A + outer((1 / (1 + sqrt(1 + sum(W * W)))) * .as_num(G %*% W), W)
  (1 / (2 * pi)) * (A - t(A))
}

#' @rdname levy_area
#' @export
levyarea <- function(W, n, alg, coeffs = NULL, rng = NULL) {
  n <- as.integer(n)
  cf <- coeffs %||% generate_coefficients(length(W), n, alg, rng)
  .la_kernel(W, n, alg, cf)
}

.ito_correction <- function(I, h = 1) {
  if (nrow(I) != ncol(I)) stop("Matrix is not square: dimensions are ", nrow(I), " x ", ncol(I), call. = FALSE)
  diag(I) <- diag(I) - h / 2
  I
}

#' @rdname levy_area
#' @export
iterated_integrals <- function(W, h, eps = NULL, ito_correction = TRUE, error_norm = NULL, alg = NULL,
                               rng = NULL, q_12 = NULL) {
  W <- .as_num(W)
  if (length(W) == 1L) {
    if (!ito_correction) return(0.5 * W^2)
    return(if (is.null(q_12)) 0.5 * W^2 - 0.5 * h else 0.5 * W^2 - 0.5 * h * q_12^2)
  }
  if (inherits(eps, "LevyAreaCoefficients")) {
    cf <- eps
    I <- levyarea(W / sqrt(h), cf$n, alg %||% MronRoe(), cf)
    if (ito_correction) I <- .ito_correction(I)
    return(0.5 * outer(W, W) + h * I)
  }
  m <- length(W)
  if (!is.null(q_12)) {
    if (is.null(eps)) stop("the Q-Wiener form needs eps", call. = FALSE)
    error_norm <- error_norm %||% FrobeniusL2()
    alg <- alg %||% optimal_algorithm(m, h, eps, error_norm, q_12 = q_12)
    n <- terms_needed(m, h, eps, alg, error_norm, q_12 = q_12)
    I <- levyarea(W / q_12 / sqrt(h), n, alg, rng = rng)
    if (ito_correction) I <- .ito_correction(I)
    return(0.5 * outer(W, W) + (matrix(h * q_12, m, m, byrow = TRUE) * I) * q_12)
  }
  eps <- eps %||% h^(3 / 2)
  error_norm <- error_norm %||% MaxL2()
  alg <- alg %||% optimal_algorithm(m, h, eps, error_norm)
  n <- terms_needed(m, h, eps, alg, error_norm)
  I <- levyarea(W / sqrt(h), n, alg, rng = rng)
  if (ito_correction) I <- .ito_correction(I)
  0.5 * outer(W, W) + h * I
}

#' @rdname levy_area
#' @param dW The Wiener increment over the step.
#' @param t_points Times in `[0, h]`.
#' @export
reconstruct_path <- function(dW, h, coeffs, t_points) {
  dW <- .as_num(dW)
  m <- length(dW)
  n <- coeffs$n
  X <- coeffs$X
  transposed <- nrow(X) == m && ncol(X) == n
  lapply(t_points, function(t) {
    Wt <- (t / h) * dW
    for (k in seq_len(n)) {
      cf <- sqrt(2 * h) / (k * pi)
      s <- sin(k * pi * t / h)
      a <- if (transposed) X[, k] else X[k, ]
      Wt <- Wt + a * cf * s
    }
    Wt
  })
}

#' @rdname levy_area
#' @param t_start,t_end The sub-interval of the step.
#' @param n_quadrature Number of Riemann sum points.
#' @export
iterated_integrals_subinterval <- function(dW, h, coeffs, t_start, t_end, n_quadrature = 64L,
                                           ito_correction = TRUE) {
  m <- length(dW)
  Wv <- reconstruct_path(dW, h, coeffs, seq(t_start, t_end, length.out = n_quadrature + 1L))
  J <- matrix(0, m, m)
  Wc <- numeric(m)
  for (s in 2:length(Wv)) {
    dWs <- Wv[[s]] - Wv[[s - 1L]]
    J <- J + outer(Wc, dWs)
    Wc <- Wc + dWs
  }
  if (ito_correction) diag(J) <- diag(J) - (t_end - t_start) / 2
  J
}

# ---------------------------------------------------------------------------
# Iterated integral evaluators of the SDE solvers (StochasticDiffEqCore)

#' Iterated integral evaluators of the SDE solvers
#'
#' Higher order SDE solvers need the Stratonovich iterated integrals `J` of
#' the Wiener increment. An algorithm chooses how to approximate them with an
#' `IteratedIntegralApprox`: `IICommutative()` assumes commutative noise, so
#' the integrals are `dW dW' / 2`, and `IILevyArea()` approximates the Levy
#' areas with the algorithm of [optimal_algorithm()]. `get_Jalg(dW, dt, prob,
#' alg)` turns that choice (the `ii_approx` of `alg`) and the noise of the
#' problem into an evaluator: `JDiagonal_oop()` (`dW^2 / 2` for diagonal or
#' scalar noise), `JCommute_oop()` (`dW dW' / 2`), a Levy area algorithm, or
#' their in place counterparts `JDiagonal_iip(dW)`, `JCommute_iip(dW)` and
#' `IteratedIntegralAlgorithm_iip(dW, levyalg)`, which keep the result in
#' their field `J`. `get_iterated_I(dt, dW, dZ, alg, p, c, gamma)` evaluates
#' it (`get_iterated_I!` in Julia, for the in place evaluators); with a Levy
#' area algorithm `p` terms are used, by default the number for the error
#' `c * dt^(gamma + 1/2)`, and the random numbers come from `coeffs` or
#' `rng` as in [levyarea()].
#'
#' @param dt The step size.
#' @param dW The Wiener increment.
#' @param dZ The second increment (unused by these evaluators).
#' @param alg An evaluator, or for `get_Jalg()` an SDE algorithm with an
#'   `ii_approx`.
#' @param p Number of Levy area terms, or `NULL`.
#' @param c,gamma Constant and exponent of the error target.
#' @param coeffs,rng Random numbers of a Levy area algorithm, see
#'   [levyarea()].
#' @param levyalg A Levy area algorithm.
#' @param prob An SDE problem.
#' @return An approximation choice, an evaluator, or the matrix (vector for
#'   diagonal noise) of iterated integrals.
#' @examples
#' dW <- c(0.1, -0.2)
#' get_iterated_I(0.01, dW, NULL, JCommute_oop())
#' J <- JDiagonal_iip(dW)
#' get_iterated_I(0.01, dW, NULL, J)
#' J$J
#' IILevyArea()
#' @name iterated_integral_evaluators
NULL

#' @rdname iterated_integral_evaluators
#' @export
IteratedIntegralApprox <- "IteratedIntegralApprox"

#' @rdname iterated_integral_evaluators
#' @export
IICommutative <- function() structure(list(name = "IICommutative"), class = c("IICommutative", "IteratedIntegralApprox"))

#' @rdname iterated_integral_evaluators
#' @export
IILevyArea <- function() structure(list(name = "IILevyArea"), class = c("IILevyArea", "IteratedIntegralApprox"))

#' @export
print.IteratedIntegralApprox <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

.jalg_iip <- function(name, J, cls, ...) {
  e <- new.env(parent = emptyenv())
  e$name <- name
  e$J <- J
  for (nm in names(list(...))) assign(nm, list(...)[[nm]], envir = e)
  class(e) <- c(name, cls)
  e
}

#' @export
print.AbstractJ <- function(x, ...) {
  cat(x$name, "\n", sep = "")
  invisible(x)
}

#' @rdname iterated_integral_evaluators
#' @export
JDiagonal_oop <- function() structure(list(name = "JDiagonal_oop"), class = c("JDiagonal_oop", "AbstractJDiagonal", "AbstractJ"))

#' @rdname iterated_integral_evaluators
#' @export
JDiagonal_iip <- function(dW) .jalg_iip("JDiagonal_iip", 0 * dW * dW, c("AbstractJDiagonal", "AbstractJ"))

#' @rdname iterated_integral_evaluators
#' @export
JCommute_oop <- function() structure(list(name = "JCommute_oop"), class = c("JCommute_oop", "AbstractJCommute", "AbstractJ"))

#' @rdname iterated_integral_evaluators
#' @export
JCommute_iip <- function(dW) .jalg_iip("JCommute_iip", 0 * outer(.as_num(dW), .as_num(dW)), c("AbstractJCommute", "AbstractJ"))

#' @rdname iterated_integral_evaluators
#' @export
IteratedIntegralAlgorithm_iip <- function(dW, levyalg)
  .jalg_iip("IteratedIntegralAlgorithm_iip", 0 * outer(.as_num(dW), .as_num(dW)),
            c("AbstractIteratedIntegralAlgorithm", "AbstractJ"), levyalg = levyalg)

.levy_J <- function(dt, dW, alg, p, c, gamma, coeffs, rng) {
  if (is.null(p)) p <- terms_needed(length(dW), dt, c * dt^(gamma + 1 / 2), alg, MaxL2())
  I <- levyarea(dW / sqrt(dt), p, alg, coeffs = coeffs, rng = rng)
  0.5 * outer(dW, dW) + dt * I
}

#' @rdname iterated_integral_evaluators
#' @export
get_iterated_I <- function(dt, dW, dZ, alg, p = NULL, c = 1, gamma = 1, coeffs = NULL, rng = NULL) {
  dW <- .as_num(dW)
  J <- if (inherits(alg, "AbstractJDiagonal")) 0.5 * dW * dW
  else if (inherits(alg, "AbstractJCommute")) 0.5 * outer(dW, dW)
  else if (inherits(alg, "IteratedIntegralAlgorithm_iip")) .levy_J(dt, dW, alg$levyalg, p, c, gamma, coeffs, rng)
  else if (inherits(alg, "AbstractIteratedIntegralAlgorithm")) .levy_J(dt, dW, alg, p, c, gamma, coeffs, rng)
  else stop("not an iterated integral evaluator", call. = FALSE)
  if (is.environment(alg)) alg$J <- J
  J
}

#' @rdname iterated_integral_evaluators
#' @export
get_Jalg <- function(dW, dt, prob, alg) {
  ii <- alg$ii_approx
  iip <- isTRUE(prob$iip)
  diag_like <- length(dW) == 1L || is_diagonal_noise(prob)
  if (inherits(ii, "IILevyArea")) {
    if (diag_like) return(if (iip) JDiagonal_iip(dW) else JDiagonal_oop())
    la <- optimal_algorithm(length(dW), dt)
    return(if (iip) IteratedIntegralAlgorithm_iip(dW, la) else la)
  }
  if (inherits(ii, "IICommutative")) {
    if (diag_like) return(if (iip) JDiagonal_iip(dW) else JDiagonal_oop())
    return(if (iip) JCommute_iip(dW) else JCommute_oop())
  }
  if (iip) IteratedIntegralAlgorithm_iip(dW, ii) else ii
}
