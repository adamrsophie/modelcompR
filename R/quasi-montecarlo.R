# Quasi-Monte Carlo point sets of QuasiMonteCarlo.jl (version 0.4.4): the
# deterministic low discrepancy samplers (Sobol, Halton, Faure, Kronecker,
# lattice rules, van der Corput and grids), the random samplers (plain Monte
# Carlo, Latin hypercube and randomized Halton), the randomizations (shift,
# Owen and Matousek scrambling and the digital shift) and the design matrices
# that GlobalSensitivity.jl draws from them. Points are the columns of a
# d x n matrix, as in Julia. QuasiMonteCarlo.jl is MIT licensed, Copyright (c)
# 2019 Chris Rackauckas; see inst/COPYRIGHTS.
#
# Random draws use R's generator (or an `rng` function of n returning n
# uniform numbers on [0, 1)); they cannot reproduce Julia's streams. The
# deterministic sequences match Julia bit for bit.

.qmc_unif <- function(rng, n) {
  if (is.null(rng)) return(stats::runif(n))
  u <- as.numeric(rng(n))
  if (length(u) != n) stop("rng(n) must return n numbers", call. = FALSE)
  u
}
# rand(rng, lo:hi, n)
.qmc_int <- function(rng, n, lo, hi) lo + pmin(floor(.qmc_unif(rng, n) * (hi - lo + 1)), hi - lo)
# shuffle(rng, x)
.qmc_shuffle <- function(rng, x) x[order(.qmc_unif(rng, length(x)))]

.qmc_check <- function(n, lb = NULL, ub = NULL) {
  if (!is.null(lb)) {
    if (length(lb) != length(ub)) stop("Dimensionality of lb and ub must match", call. = FALSE)
    if (any(lb > ub)) stop("Lower bound exceeds upper bound (lb > ub)", call. = FALSE)
  }
  if (!(n > 0)) stop("Number of samples must be greater than zero", call. = FALSE)
  invisible(NULL)
}

# log(b, n) of Julia: log(n) / log(b)
.qmc_log <- function(b, n) log(n) / log(b)

.qmc_logi <- function(b, n) {
  m <- round(.qmc_log(b, n))
  if (b^m != n) stop(n, " is not a power of ", b, call. = FALSE)
  m
}

.qmc_isprime <- function(n) {
  if (n < 2) return(FALSE)
  if (n < 4) return(TRUE)
  if (n %% 2 == 0) return(FALSE)
  k <- 3
  while (k * k <= n) {
    if (n %% k == 0) return(FALSE)
    k <- k + 2
  }
  TRUE
}
.qmc_nextprime <- function(n) {
  while (!.qmc_isprime(n)) n <- n + 1
  n
}
# nextprimes(start, k): the first k primes >= start
.qmc_nextprimes <- function(start, k) {
  out <- numeric(k)
  p <- start
  for (i in seq_len(k)) {
    p <- .qmc_nextprime(p)
    out[i] <- p
    p <- p + 1
  }
  out
}

.qmc_sampler <- function(name, ..., random = FALSE) {
  structure(list(...), class = c(name, if (random) "RandomSamplingAlgorithm" else "DeterministicSamplingAlgorithm",
                                 "SamplingAlgorithm"))
}
.qmc_rand_method <- function(name, ...) structure(list(...), class = c(name, "RandomizationMethod"))

#' Quasi-Monte Carlo samplers
#'
#' Ports of the point set generators of QuasiMonteCarlo.jl. Pass a sampler to
#' [qmc_sample()] or [generate_design_matrices()]. Deterministic samplers
#' (`SobolSample()`, `HaltonSample()`, `FaureSample()`, `KroneckerSample()`,
#' `GoldenSample()`, `LatticeRuleSample()`, `VanDerCorputSample()` and
#' `GridSample()`) carry a default randomization `R` (see
#' [randomization_methods]); random samplers (`RandomSample()`,
#' `LatinHypercubeSample()` and `RandomizedHaltonSample()`) draw a new point
#' set at every call.
#'
#' The deterministic sequences match Julia bit for bit. Random draws use R's
#' generator, so seeded results differ from Julia's.
#'
#' @param R Randomization method applied to the point set, `NoRand()` by
#'   default.
#' @param rng `NULL` for R's generator, or a function of `n` returning `n`
#'   uniform numbers on \[0, 1).
#' @param base Prime base of the van der Corput sequence.
#' @param generator Generating vector of the Kronecker sequence; `NULL` takes
#'   the powers of the harmonious number of the dimension (the generalized
#'   golden ratio).
#' @param ... Arguments of `KroneckerSample()`.
#' @return A sampler object.
#' @examples
#' qmc_sample(8, 2, SobolSample())
#' qmc_sample(9, c(0, -1), c(1, 1), HaltonSample())
#' @name qmc_samplers
NULL

#' @rdname qmc_samplers
#' @export
SobolSample <- function(R = NoRand()) .qmc_sampler("SobolSample", R = R)
#' @rdname qmc_samplers
#' @export
HaltonSample <- function(R = NoRand()) .qmc_sampler("HaltonSample", R = R)
#' @rdname qmc_samplers
#' @export
FaureSample <- function(R = NoRand()) .qmc_sampler("FaureSample", R = R)
#' @rdname qmc_samplers
#' @export
LatticeRuleSample <- function(R = NoRand()) .qmc_sampler("LatticeRuleSample", R = R)
#' @rdname qmc_samplers
#' @export
GridSample <- function(R = NoRand()) .qmc_sampler("GridSample", R = R)
#' @rdname qmc_samplers
#' @export
VanDerCorputSample <- function(base, R = NoRand()) .qmc_sampler("VanDerCorputSample", base = base, R = R)
#' @rdname qmc_samplers
#' @export
KroneckerSample <- function(generator = NULL, R = NoRand()) {
  .qmc_sampler("KroneckerSample", generator = if (!is.null(generator)) as.numeric(generator), R = R)
}
#' @rdname qmc_samplers
#' @export
GoldenSample <- function(...) KroneckerSample(...)
#' @rdname qmc_samplers
#' @export
RandomSample <- function(rng = NULL) .qmc_sampler("RandomSample", rng = rng, random = TRUE)
#' @rdname qmc_samplers
#' @export
LatinHypercubeSample <- function(rng = NULL) .qmc_sampler("LatinHypercubeSample", rng = rng, random = TRUE)
#' @rdname qmc_samplers
#' @export
RandomizedHaltonSample <- function(rng = NULL) .qmc_sampler("RandomizedHaltonSample", rng = rng, random = TRUE)

#' Randomization of quasi-Monte Carlo point sets
#'
#' `NoRand()` leaves the points as they are. `QMCShift()` adds one uniform
#' random shift modulo 1 to every point (QuasiMonteCarlo's `Shift`, renamed
#' because [Shift()] is the clock operator of ModelingToolkit).
#' `OwenScramble()`, `MatousekScramble()` and `DigitalShift()` scramble the
#' base `base` digits of point sets of `base^m` points: nested uniform
#' permutations, a random linear matrix scramble plus digital shift, and a
#' digital shift alone. Digits beyond the first `m` of the `pad` digits are
#' drawn at random.
#'
#' @param x A d x n matrix of points in the unit box.
#' @param R A randomization method.
#' @param base Base of the digital expansion.
#' @param pad Number of digits kept per coordinate.
#' @param rng `NULL` for R's generator, or a function of `n` returning `n`
#'   uniform numbers on \[0, 1).
#' @return `randomize()` returns the randomized matrix; the constructors return
#'   randomization methods.
#' @examples
#' x <- qmc_sample(16, 2, SobolSample())
#' set.seed(1)
#' randomize(x, OwenScramble(base = 2))
#' @name randomization_methods
NULL

#' @rdname randomization_methods
#' @export
NoRand <- function() .qmc_rand_method("NoRand")
#' @rdname randomization_methods
#' @export
QMCShift <- function(rng = NULL) .qmc_rand_method("QMCShift", rng = rng)
#' @rdname randomization_methods
#' @export
OwenScramble <- function(base, pad = 32L, rng = NULL)
  .qmc_rand_method(c("OwenScramble", "ScrambleMethod"), base = base, pad = as.integer(pad), rng = rng)
#' @rdname randomization_methods
#' @export
MatousekScramble <- function(base, pad = 32L, rng = NULL)
  .qmc_rand_method(c("MatousekScramble", "ScrambleMethod"), base = base, pad = as.integer(pad), rng = rng)
#' @rdname randomization_methods
#' @export
DigitalShift <- function(base, pad = 32L, rng = NULL)
  .qmc_rand_method(c("DigitalShift", "ScrambleMethod"), base = base, pad = as.integer(pad), rng = rng)

#' Draw quasi-Monte Carlo points
#'
#' `qmc_sample(n, d, sampler)` returns a d x n matrix of points in the unit
#' box; `qmc_sample(n, lb, ub, sampler)` maps them to the box with lower
#' bounds `lb` and upper bounds `ub`. This is `QuasiMonteCarlo.sample`, renamed
#' because [base::sample()] is not generic. A function of `k` returning `k`
#' random numbers (for example `function(k) rnorm(k)`) is also accepted as a
#' sampler, as the Distributions extension of QuasiMonteCarlo accepts a
#' distribution.
#'
#' @param n Number of points.
#' @param d Dimension, or the lower bounds `lb` when `ub` is given.
#' @param ub Upper bounds, or the sampler in the unit box form.
#' @param sampler A sampler from [qmc_samplers].
#' @return A matrix with one point per column.
#' @examples
#' qmc_sample(4, 1, VanDerCorputSample(2))
#' @export
qmc_sample <- function(n, d, ub, sampler) {
  if (missing(sampler)) return(.qmc_unit(n, d, ub))
  lb <- as.numeric(d); ub <- as.numeric(ub)
  .qmc_check(n, lb, ub)
  (ub - lb) * .qmc_unit(n, length(lb), sampler) + lb
}

.qmc_unit <- function(n, d, S) {
  if (is.function(S)) {
    .qmc_check(n)
    return(matrix(as.numeric(S(n * d)), d, n))
  }
  if (!inherits(S, "SamplingAlgorithm")) stop("`sampler` must be a sampling algorithm", call. = FALSE)
  switch(class(S)[1L],
    RandomSample = { .qmc_check(n); matrix(.qmc_unif(S$rng, d * n), d, n) },
    LatinHypercubeSample = {
      .qmc_check(n)
      seq <- (seq_len(n) - 0.5) / n
      out <- matrix(0, d, n)
      for (i in seq_len(d)) out[i, ] <- .qmc_shuffle(S$rng, seq)
      out
    },
    RandomizedHaltonSample = {
      .qmc_check(n)
      bases <- .qmc_nextprimes(1, d)
      out <- matrix(0, d, n)
      for (i in seq_len(d)) out[i, ] <- .qmc_randradinv(seq_len(n), S$rng, bases[i])
      out
    },
    SobolSample = {
      if (n < 0) stop("number of samples must be non-negative", call. = FALSE)
      seq <- matrix(0, d, n)
      if (n == 0) return(seq)
      s <- .sobol_new(d)
      for (k in seq_len(2^floor(log2(n + 1)) - 1)) .sobol_next(s)
      for (j in seq_len(n)) seq[, j] <- .sobol_next(s)
      randomize(seq, S$R)
    },
    HaltonSample = {
      bases <- .qmc_nextprimes(1, d)
      out <- matrix(0, d, n)
      for (i in seq_len(d)) out[i, ] <- .qmc_vdc(ceiling(.qmc_log(bases[i], n)), bases[i], n)
      randomize(out, S$R)
    },
    VanDerCorputSample = {
      if (d != 1) stop("Van der Corput sequence only supports 1D sampling", call. = FALSE)
      .qmc_check(n)
      randomize(matrix(.qmc_vdc(ceiling(.qmc_log(S$base, n)), S$base, n), 1L, n), S$R)
    },
    FaureSample = .qmc_faure(n, d, S),
    KroneckerSample = {
      g <- S$generator %||% .qmc_harmonious_powers(d)
      if (length(g) != d) stop("Dimensions of generator and sample must match.", call. = FALSE)
      randomize(outer(g, seq_len(n)) %% 1, S$R)
    },
    LatticeRuleSample = randomize(.qmc_lattice(n, d), S$R),
    GridSample = {
      out <- matrix(0, d, n)
      for (i in seq_len(d)) out[i, ] <- (.qmc_int(NULL, n, 0, n)) / n
      randomize(out, S$R)
    },
    stop("no sample method for ", class(S)[1L], call. = FALSE)
  )
}

# _vdc(lambda, n_digits, base; n): the n first points of the van der Corput
# sequence in base `base`, centered in their cells
.qmc_vdc <- function(n_digits, base, n) {
  inv_base <- 1 / base
  dgs <- rep(base - 1, n_digits)
  offset <- 1 / (2 * n)
  out <- numeric(n)
  for (k in seq_len(n)) {
    dgs <- .qmc_base_sum(dgs, base)
    out[k] <- .jl_horner(inv_base, dgs) * inv_base + offset
  }
  out
}

# _base_sum!: add one to the little endian digits, ignoring overflow
.qmc_base_sum <- function(dgs, base) {
  i <- 1L
  repeat {
    dgs[i] <- dgs[i] + 1
    if (dgs[i] != base) break
    dgs[i] <- 0
    if (i == length(dgs)) break
    i <- i + 1L
  }
  dgs
}

.qmc_faure <- function(n, d, S) {
  base <- .qmc_nextprime(d)
  n_digits <- ceiling(.qmc_log(base, n))
  # prevpow(base, n)
  pp <- base^floor(.qmc_log(base, n))
  if (pp > n) pp <- pp / base
  if (pp * base <= n) pp <- pp * base
  if (n %% pp != 0) {
    m <- n - n %% base^(n_digits - 1)
    stop("Invalid sample size: Faure sequences must be multiples of `base^power`. Try ", m, " or ",
         m + base^(n_digits - 1), " instead.", call. = FALSE)
  }
  inv_base <- 1 / base
  pascal <- matrix(0, n_digits, n_digits)
  pascal[upper.tri(pascal, diag = TRUE)] <- 1
  for (i in seq_len(n_digits)[-1L]) for (j in 2:i) pascal[j, i] <- (pascal[j - 1L, i - 1L] + pascal[j, i - 1L]) %% base
  faure <- matrix(0, d, n)
  half <- 1 / (2 * base^n_digits)
  for (dim in seq_len(d)) {
    counter <- rep(base - 1, n_digits)
    if (dim != 1L) {
      perm <- matrix(0, n_digits, n_digits)
      for (j in seq_len(n_digits)) for (i in seq_len(j)) {
        # powermod(power, j - i, base)
        pw <- 1
        for (k in seq_len(j - i)) pw <- (pw * (dim - 1)) %% base
        perm[i, j] <- pw * pascal[i, j]
      }
    }
    for (k in seq_len(n)) {
      counter <- .qmc_base_sum(counter, base)
      dgs <- if (dim != 1L) as.numeric(perm %*% counter) %% base else counter
      faure[dim, k] <- (.jl_horner(inv_base, dgs) + half) * inv_base
    }
  }
  randomize(faure, S$R)
}

# harmonious(d): the fixed point of y = (1 + y)^(1 / (d + 1)), inverted
.qmc_harmonious_powers <- function(d) {
  tol <- 2 * .Machine$double.eps
  y_old <- 1
  e <- 1 / (d + 1)
  y <- (1 + y_old)^e
  while (abs(y - y_old) > tol) {
    y_old <- y
    y <- (1 + y)^e
  }
  r <- 1 / y
  vapply(seq_len(d), function(i) .simple_powi(r, i), numeric(1))
}

.qmc_reversebits <- function(k) sum(((k %/% 2^(0:31)) %% 2) * 2^(31:0))

.qmc_lattice <- function(n, d) {
  if (!(d > 0)) stop("number of dimensions d must be larger than 0", call. = FALSE)
  if (d <= 250) { z <- .qmc_CKN_250_20[seq_len(d)]; nmax <- 2^20 } else {
    if (d > 3600) stop("number of dimensions s must be less than or equal to 3600", call. = FALSE)
    z <- .qmc_K_3600_32[seq_len(d)]; nmax <- 2^32
  }
  if (n > nmax) stop("lattice index out of bounds", call. = FALSE)
  out <- matrix(0, d, n)
  for (k in seq_len(n) - 1) {
    phi <- .qmc_reversebits(k) * 2^-32
    x <- phi * z
    out[, k + 1] <- x - floor(x)
  }
  out
}

.qmc_randradinv <- function(ind, rng, b) {
  b2r <- 1 / b
  ans <- numeric(length(ind))
  res <- ind
  while (1 - b2r < 1) {
    dig <- res %% b
    perm <- .qmc_shuffle(rng, seq_len(b)) - 1
    ans <- ans + perm[dig + 1] * b2r
    b2r <- b2r / b
    res <- (res - dig) / b
  }
  ans
}

#' @rdname randomization_methods
#' @export
randomize <- function(x, R) {
  if (!inherits(R, "RandomizationMethod")) stop("`R` must be a randomization method", call. = FALSE)
  if (inherits(R, "NoRand")) return(x)
  x <- as.matrix(x)
  if (inherits(R, "QMCShift")) {
    U <- .qmc_unif(R$rng, nrow(x))
    y <- x + U
    return(y - floor(y))
  }
  # ScrambleMethod: digits of the n x d transpose
  pts <- t(x)
  bits <- .qmc_unif2bits(pts, R$base, R$pad)
  rb <- if (inherits(R, "OwenScramble")) .qmc_owen_bits(bits, .qmc_which_perm(bits, R$base), R)
        else .qmc_scramble_bits(bits, R)
  t(.qmc_bits2unif(rb, R$base))
}

# unif2bits: the first pad digits in base b of every entry, as a
# pad x n x d array; divrem with an exact remainder (fmod)
.qmc_unif2bits <- function(x, b, pad = 32L) {
  dims <- dim(x) %||% length(x)
  y <- as.numeric(x)
  out <- array(0, c(pad, length(y)))
  ib <- 1 / b
  for (i in seq_len(pad)) {
    z <- .simple_powi(ib, i)
    q <- floor(y / z)
    r <- .fma(-q, z, y)
    lo <- r < 0
    if (any(lo)) { q[lo] <- q[lo] - 1; r[lo] <- .fma(-q[lo], z, y[lo]) }
    hi <- r >= z
    if (any(hi)) { q[hi] <- q[hi] + 1; r[hi] <- .fma(-q[hi], z, y[hi]) }
    out[i, ] <- q
    y <- r
  }
  array(out, c(pad, dims))
}

.qmc_bits2unif <- function(bits, b) {
  d <- dim(bits)
  m <- matrix(bits, d[1L])
  y <- numeric(ncol(m))
  for (j in rev(seq_len(d[1L]))) y <- (y + m[j, ]) / b
  array(y, d[-1L])
}

.qmc_which_perm <- function(bits, b) {
  n <- dim(bits)[2L]; d <- dim(bits)[3L]
  m <- .qmc_logi(b, n)
  idx <- array(0, c(m, n, d))
  for (k in seq_len(m)[-1L]) {
    v <- 0
    for (kk in (k - 1L):1L) v <- v * b + bits[kk, , , drop = FALSE]
    idx[k, , ] <- v
  }
  idx + 1
}

.qmc_pad_random <- function(rb, m, R) {
  pad <- dim(rb)[1L]
  if (pad > m) {
    k <- (pad - m) * prod(dim(rb)[-1L])
    rb[(m + 1L):pad, , ] <- .qmc_int(R$rng, k, 0, R$base - 1)
  }
  rb
}

.qmc_owen_bits <- function(bits, idx, R) {
  m <- dim(idx)[1L]; d <- dim(idx)[3L]; b <- R$base
  if (m < 1) stop("We need m >= 1", call. = FALSE)
  rb <- bits
  for (s in seq_len(d)) {
    perms <- matrix(0, m, b^(m - 1))
    for (k in seq_len(m)) perms[k, seq_len(b^(k - 1))] <- .qmc_int(R$rng, b^(k - 1), 0, b - 1)
    for (k in seq_len(m)) rb[k, , s] <- (bits[k, , s] + perms[k, idx[k, , s]]) %% b
  }
  .qmc_pad_random(rb, m, R)
}

.qmc_scramble_bits <- function(bits, R) {
  n <- dim(bits)[2L]; d <- dim(bits)[3L]; b <- R$base
  m <- .qmc_logi(b, n)
  if (m < 1) stop("We need m >= 1", call. = FALSE)
  rb <- bits
  for (s in seq_len(d)) {
    B <- matrix(bits[seq_len(m), , s], m)
    if (inherits(R, "MatousekScramble")) {
      M <- diag(.qmc_int(R$rng, m, 1, b - 1), m)
      C <- .qmc_int(R$rng, m, 0, b - 1)
      for (i in seq_len(m)[-1L]) for (j in seq_len(i - 1L)) M[i, j] <- .qmc_int(R$rng, 1, 0, b - 1)
      rb[seq_len(m), , s] <- (M %*% B + C) %% b
    } else {
      DS <- .qmc_int(R$rng, m, 0, b - 1)
      rb[seq_len(m), , s] <- (B + DS) %% b
    }
  }
  .qmc_pad_random(rb, m, R)
}

#' Design matrices for quasi-Monte Carlo error estimation
#'
#' `generate_design_matrices(n, lb, ub, sampler, num_mats)` returns a list of
#' `num_mats` matrices of `n` points in the box from `lb` to `ub`;
#' `generate_design_matrices(n, d, sampler, num_mats)` draws them in the unit
#' box. A deterministic sampler is randomized independently for every matrix
#' by its `R` field, or by `R` when given; with `NoRand()` the matrices are
#' consecutive blocks of dimensions of one `num_mats * d` dimensional
#' sequence, which are not independent (QuasiMonteCarlo warns). A random
#' sampler draws every matrix anew. `DesignMatrix()` is the iterator form:
#' `$next_matrix()` returns the next matrix.
#'
#' The randomized matrices scramble the same `n` base points, so each one is a
#' valid randomized QMC point set and the estimates from different matrices
#' are independent replicates, but the matrices are dependent column by
#' column (the i-th points of two matrices are strongly correlated). They
#' suit replicated error estimation, not designs that need independent pairs
#' of points such as the `A` and `B` matrices of Sobol indices; for those,
#' randomize one `2 d` dimensional point set and split its rows, as [gsa()]
#' does.
#'
#' @param n Number of points per matrix.
#' @param d Dimension, or the lower bounds when `ub` is given.
#' @param ub Upper bounds (bounds form only).
#' @param sampler A sampler from [qmc_samplers].
#' @param num_mats Number of matrices.
#' @param R Randomization method; the sampler's `R` by default.
#' @return A list of d x n matrices; `DesignMatrix()` returns an iterator.
#' @examples
#' set.seed(1)
#' AB <- generate_design_matrices(16, c(0, 0), c(1, 2), SobolSample(R = QMCShift()), 2)
#' lengths(AB)
#' @export
generate_design_matrices <- function(n, d, ub = NULL, sampler, num_mats = 2L, R = NULL) {
  if (!is.null(ub)) {
    if (!(n > 0)) stop("Number of samples must be greater than zero", call. = FALSE)
    lb <- as.numeric(d); ub <- as.numeric(ub)
    if (length(lb) != length(ub)) stop("Dimensionality of lb and ub must match", call. = FALSE)
    out <- generate_design_matrices(n, length(lb), NULL, sampler, num_mats, R)
    return(lapply(out, function(m) (ub - lb) * m + lb))
  }
  if (inherits(sampler, "RandomSamplingAlgorithm") || is.function(sampler))
    return(lapply(seq_len(num_mats), function(j) .qmc_unit(n, d, sampler)))
  R <- R %||% sampler$R
  if (inherits(R, "NoRand")) {
    out <- .qmc_unit(n, num_mats * d, sampler)
    warning("The `generate_design_matrices(n, d, sampler, R = NoRand(), num_mats)` method does not produce true ",
            "and independent QMC matrices. Prefer a randomization such as `R = QMCShift()` or ",
            "`R = MatousekScramble(base = 2)`.", call. = FALSE)
    return(lapply(seq_len(num_mats) - 1L, function(j) out[j * d + seq_len(d), , drop = FALSE]))
  }
  dm <- DesignMatrix(n, d, sampler, num_mats, R)
  lapply(seq_len(num_mats), function(j) dm$next_matrix())
}

#' @rdname generate_design_matrices
#' @export
DesignMatrix <- function(n, d, sampler, num_mats, R = NULL) {
  e <- new.env(parent = emptyenv())
  e$count <- num_mats
  if (inherits(sampler, "RandomSamplingAlgorithm") || is.function(sampler)) {
    e$next_matrix <- function() .qmc_unit(n, d, sampler)
  } else {
    R <- R %||% sampler$R
    s0 <- sampler; s0$R <- NoRand()
    X <- .qmc_unit(n, d, s0)
    if (inherits(R, "ScrambleMethod")) {
      bits <- .qmc_unif2bits(t(X), R$base, R$pad)
      idx <- if (inherits(R, "OwenScramble")) .qmc_which_perm(bits, R$base)
      e$next_matrix <- function() {
        rb <- if (!is.null(idx)) .qmc_owen_bits(bits, idx, R) else .qmc_scramble_bits(bits, R)
        t(.qmc_bits2unif(rb, R$base))
      }
    } else if (inherits(R, "QMCShift")) {
      # randomize!(X, Shift) shifts the same matrix again at every call
      e$next_matrix <- function() {
        X <<- randomize(X, R)
        X
      }
    } else {
      e$next_matrix <- function() X
    }
  }
  structure(e, class = "DesignMatrix")
}

#' @export
length.DesignMatrix <- function(x) x$count
