# Global sensitivity analysis of GlobalSensitivity.jl (version 2.12.9): the
# Sobol, Morris, regression, DGSM, eFAST, delta moment independent, EASI,
# fractional factorial, RBD-FAST, Shapley, regional (RSA) and mutual
# information methods behind gsa(). GlobalSensitivity.jl is MIT licensed,
# Copyright (c) 2020 Vaibhav Dixit; the kernel density estimate of
# KernelDensity.jl, the quadratic B-spline of Interpolations.jl, trapz of
# Trapz.jl and the histogram entropy of ComplexityMeasures.jl that the delta
# and mutual information methods use are ported here as well (MIT licensed;
# see inst/COPYRIGHTS).
#
# Parameter sets are the columns of d x n matrices, as in Julia. Random draws
# use R's generator (or an `rng` function of n returning n uniform numbers on
# [0, 1)) and so differ from Julia's streams; the analyses of given samples
# follow Julia's arithmetic.

.gsa_method <- function(name, ...) structure(list(...), class = c(name, "GSAMethod"))

#' Global sensitivity analysis methods
#'
#' Method objects for [gsa()], ports of the methods of GlobalSensitivity.jl.
#'
#' * `Sobol()`: variance based first, second and total order indices, with the
#'   Monte Carlo estimators of Jansen (1999, default), Homma and Saltelli
#'   (1996), Sobol (2007) and Janon (2014) for the total effects, and
#'   confidence intervals over `nboot` independent replicates.
#' * `Morris()`: elementary effects along one at a time trajectories on a grid
#'   of `p_steps` levels; the `num_trajectory` trajectories of largest spread
#'   among `total_num_trajectory` candidates are kept.
#' * `RegressionGSA()`: Pearson correlations, standardized regression and
#'   partial correlation coefficients, and their rank versions.
#' * `DGSM()`: derivative based global sensitivity measures from forward mode
#'   gradients (and Hessians when `crossed`): the means of the derivative,
#'   its absolute value and its square, and for uniform inputs the measures
#'   `tao` and `sigma` of Kucherenko and Iooss (2016), whose ratios to the
#'   output variance bound the total Sobol index (`NA` for other laws).
#' * `eFAST()`: the extended Fourier amplitude sensitivity test.
#' * `DeltaMoment()`: the moment independent delta measure of Borgonovo, with
#'   kernel density estimates and a bootstrap bias correction.
#' * `EASI()`: first order indices from spectral analysis of sorted samples.
#' * `FractionalFactorial()`: main effects of a resolution IV two level design.
#' * `RBDFAST()`: random balance designs Fourier amplitude sensitivity test,
#'   with the bias correction of Tissot and Prieur (2012).
#' * `Shapley()`: Shapley effects for dependent inputs, exact over all
#'   permutations (`n_perms = -1`) or over random permutations.
#' * `RSA()`: regional sensitivity analysis (Kolmogorov Smirnov distances
#'   between accepted and rejected samples), with dummy parameters.
#' * `MutualInformation()`: mutual information between each input and the
#'   output from histograms, against a permutation bound.
#'
#' @param order Orders of the Sobol indices: 0 (total), 1 (first), 2 (second).
#' @param nboot Number of bootstrap replicates.
#' @param conf_level Confidence level.
#' @param p_steps Number of grid levels per parameter (100 when missing).
#' @param relative_scale Scale the elementary effects by the output.
#' @param num_trajectory,total_num_trajectory Trajectories kept and drawn.
#' @param len_design_mat Points per trajectory.
#' @param rank Also compute the rank based coefficients.
#' @param crossed Also compute the crossed second derivative measures.
#' @param num_harmonics,max_harmonic Number of harmonics.
#' @param Ygrid_length Number of quadrature points of the output densities.
#' @param num_classes Number of classes of the inputs; `NULL` chooses it from
#'   the sample size.
#' @param dct_method Use the discrete cosine transform in EASI.
#' @param n_perms Number of random permutations, `-1` for all permutations.
#' @param n_var,n_outer,n_inner Sample sizes of the output variance, the
#'   outer loop and the inner conditional loop.
#' @param n_dummy_parameters Number of dummy parameters of RSA.
#' @param acceptance_threshold Threshold of acceptance, or a function of the
#'   outputs giving it.
#' @param n_bootstraps Number of permutations for the mutual information bound.
#' @return A method object.
#' @name gsa_methods
NULL

#' @rdname gsa_methods
#' @export
Sobol <- function(order = c(0L, 1L), nboot = 1L, conf_level = 0.95) {
  if (!(nboot >= 1)) stop("`nboot` must be positive", call. = FALSE)
  .gsa_method("Sobol", order = as.integer(order), nboot = as.integer(nboot), conf_level = conf_level)
}
#' @rdname gsa_methods
#' @export
Morris <- function(p_steps = integer(), relative_scale = FALSE, num_trajectory = 10L,
                   total_num_trajectory = 5L * num_trajectory, len_design_mat = 10L) {
  .gsa_method("Morris", p_steps = as.integer(p_steps), relative_scale = isTRUE(relative_scale),
              num_trajectory = as.integer(num_trajectory), total_num_trajectory = as.integer(total_num_trajectory),
              len_design_mat = as.integer(len_design_mat))
}
#' @rdname gsa_methods
#' @export
RegressionGSA <- function(rank = FALSE) .gsa_method("RegressionGSA", rank = isTRUE(rank))
#' @rdname gsa_methods
#' @export
DGSM <- function(crossed = FALSE) .gsa_method("DGSM", crossed = isTRUE(crossed))
#' @rdname gsa_methods
#' @export
eFAST <- function(num_harmonics = 4L) .gsa_method("eFAST", num_harmonics = as.integer(num_harmonics))
#' @rdname gsa_methods
#' @export
DeltaMoment <- function(nboot = 500L, conf_level = 0.95, Ygrid_length = 2048L, num_classes = NULL) {
  .gsa_method("DeltaMoment", nboot = as.integer(nboot), conf_level = conf_level,
              Ygrid_length = as.integer(Ygrid_length), num_classes = num_classes)
}
#' @rdname gsa_methods
#' @export
EASI <- function(max_harmonic = 4L, dct_method = FALSE)
  .gsa_method("EASI", max_harmonic = as.integer(max_harmonic), dct_method = isTRUE(dct_method))
#' @rdname gsa_methods
#' @export
FractionalFactorial <- function() .gsa_method("FractionalFactorial")
#' @rdname gsa_methods
#' @export
RBDFAST <- function(num_harmonics = 6L) .gsa_method("RBDFAST", num_harmonics = as.integer(num_harmonics))
#' @rdname gsa_methods
#' @export
Shapley <- function(n_perms = -1L, n_var, n_outer, n_inner = 3L) {
  .gsa_method("Shapley", n_perms = as.integer(n_perms), n_var = as.integer(n_var),
              n_outer = as.integer(n_outer), n_inner = as.integer(n_inner))
}
#' @rdname gsa_methods
#' @export
RSA <- function(n_dummy_parameters = 10L, acceptance_threshold = mean)
  .gsa_method("RSA", n_dummy_parameters = as.integer(n_dummy_parameters), acceptance_threshold = acceptance_threshold)
#' @rdname gsa_methods
#' @export
MutualInformation <- function(n_bootstraps = 1000L, conf_level = 0.95)
  .gsa_method("MutualInformation", n_bootstraps = as.integer(n_bootstraps), conf_level = conf_level)

#' Global sensitivity analysis
#'
#' `gsa(f, method, p_range, samples = n)` evaluates the model `f` on a design
#' over the parameter ranges and returns the sensitivity measures of
#' `method` (see [gsa_methods]). `p_range` is a list of `c(lower, upper)`
#' pairs, or of distributions for `eFAST()` and `DGSM()`, a [SklarDist()] for
#' `Shapley()`. `f` takes a parameter vector and returns a number, a vector
#' or an ODE solution; with `batch = TRUE` it takes a matrix with one
#' parameter set per column and returns one value (or one column) per set.
#'
#' The samples of the range forms are randomized quasi Monte Carlo (Owen
#' scrambled Sobol points), except for `EASI()`, whose bias correction
#' assumes independent uniform samples and which draws those. `Sobol()` evaluates `nboot` independent scrambled
#' designs of dimension `2 d`, whose first and last `d` rows are the `A` and
#' `B` matrices, so the estimators are applied to independent pairs of
#' points, repeated runs are independent replicates, and the intervals with
#' `nboot > 1` come from independent replicates. GlobalSensitivity.jl uses
#' the unscrambled sequence; pass matrices from [generate_design_matrices()]
#' with `NoRand()` to `gsa(f, Sobol(), A, B)` for its deterministic design.
#'
#' Other forms, as in GlobalSensitivity.jl:
#' * `gsa(f, Sobol(), A, B)` with design matrices `A` and `B` (one set per
#'   column), for example from [generate_design_matrices()];
#' * `gsa(X, Y, method)` for `RegressionGSA()`, `DeltaMoment()` and `EASI()`
#'   on samples `X` (d x n) and outputs `Y`;
#' * `gsa(f, FractionalFactorial(), num_params = k, p_range = levels)` and
#'   `gsa(f, RBDFAST(), num_params = k, samples = n)`.
#'
#' @param f The model, or the sample matrix `X`.
#' @param method A method from [gsa_methods], or the outputs `Y`.
#' @param p_range Parameter ranges or distributions, the design matrix `A`,
#'   or the method in the `gsa(X, Y, method)` form.
#' @param B Second design matrix of the Sobol method.
#' @param samples Number of samples.
#' @param batch Whether `f` evaluates a matrix of parameter sets.
#' @param rng `NULL` for R's generator, or a function of `n` returning `n`
#'   uniform numbers on \[0, 1).
#' @param Ei_estimator Estimator of the total effects of `Sobol()`:
#'   `"Jansen1999"`, `"Homma1996"`, `"Sobol2007"` or `"Janon2014"`.
#' @param num_params Number of parameters (fractional factorial and
#'   RBD-FAST).
#' @param ... Unused.
#' @return A result list whose fields follow the Julia result types, for
#'   example `S1`, `ST` and `S2` for Sobol, `means`, `means_star` and
#'   `variances` for Morris.
#' @examples
#' ishi <- function(x) sin(x[1]) + 7 * sin(x[2])^2 + 0.1 * x[3]^4 * sin(x[1])
#' rg <- rep(list(c(-pi, pi)), 3)
#' set.seed(1)
#' res <- gsa(ishi, Sobol(), rg, samples = 1000)
#' res$S1
#' set.seed(1)
#' gsa(ishi, eFAST(), rg, samples = 500)$S1
#' @export
gsa <- function(f, method, p_range = NULL, B = NULL, samples = NULL, batch = FALSE, rng = NULL,
                Ei_estimator = "Jansen1999", num_params = NULL, ...) {
  if (!is.function(f)) {
    X <- as.matrix(f); Y <- method; m <- p_range
    if (!inherits(m, "GSAMethod")) stop("`gsa(X, Y, method)` needs a method", call. = FALSE)
    return(switch(class(m)[1L],
      RegressionGSA = .gsa_regression(X, Y, m),
      DeltaMoment = .gsa_delta(X, as.numeric(Y), m, rng),
      EASI = .gsa_easi(X, Y, m),
      stop("the gsa(X, Y, method) form is not available for ", class(m)[1L], call. = FALSE)))
  }
  if (!inherits(method, "GSAMethod")) stop("`method` must be a GSA method", call. = FALSE)
  switch(class(method)[1L],
    Sobol = {
      if (is.matrix(p_range)) return(.gsa_sobol_AB(f, method, p_range, B, batch, Ei_estimator))
      .gsa_need(samples)
      r <- .gsa_bounds(p_range)
      # nboot independent scrambled Sobol designs of dimension 2 d, split
      # into A and B (see .gsa_rqmc)
      d <- length(r$lb)
      U <- lapply(seq_len(method$nboot), function(b) .gsa_rqmc(samples, c(r$lb, r$lb), c(r$ub, r$ub), rng))
      .gsa_sobol_AB(f, method, do.call(cbind, lapply(U, function(u) u[seq_len(d), , drop = FALSE])),
                    do.call(cbind, lapply(U, function(u) u[d + seq_len(d), , drop = FALSE])), batch, Ei_estimator)
    },
    Morris = .gsa_morris(f, method, p_range, batch, rng),
    RegressionGSA = {
      .gsa_need(samples)
      r <- .gsa_bounds(p_range)
      X <- .gsa_rqmc(samples, r$lb, r$ub, rng)
      ev <- .gsa_eval(f, X, batch)
      Y <- if (ev$multi) ev$Y else matrix(ev$Y, 1L)
      .gsa_regression(X, Y, method)
    },
    DGSM = { .gsa_need(samples); .gsa_dgsm(f, method, p_range, samples) },
    eFAST = { .gsa_need(samples); .gsa_efast(f, method, p_range, samples, batch, rng) },
    DeltaMoment = {
      .gsa_need(samples)
      r <- .gsa_bounds(p_range)
      X <- .gsa_rqmc(samples, r$lb, r$ub, rng)
      ev <- .gsa_eval(f, X, batch)
      if (ev$multi) stop("DeltaMoment sensitivity only supports scalar output functions", call. = FALSE)
      .gsa_delta(X, ev$Y, method, rng)
    },
    EASI = {
      .gsa_need(samples)
      r <- .gsa_bounds(p_range)
      # independent uniform samples: the S1_Corr correction assumes the
      # spectrum of the residual noise is flat, which low discrepancy
      # points break (it then overcorrects)
      X <- r$lb + (r$ub - r$lb) * matrix(.qmc_unif(rng, length(r$lb) * samples), length(r$lb))
      ev <- .gsa_eval(f, X, batch)
      .gsa_easi(X, ev$Y, method)
    },
    FractionalFactorial = .gsa_ff(f, num_params, p_range),
    RBDFAST = { .gsa_need(samples); .gsa_rbdfast(f, method, num_params, samples, rng, batch) },
    Shapley = .gsa_shapley(f, method, p_range, batch),
    RSA = { .gsa_need(samples); .gsa_rsa(f, method, p_range, samples, batch) },
    MutualInformation = { .gsa_need(samples); .gsa_mi(f, method, p_range, samples, batch, rng) },
    stop("unknown GSA method", call. = FALSE))
}

# A randomized quasi Monte Carlo sample: Owen scrambled Sobol points, each
# uniform on the box, so the Monte Carlo estimators keep their expectations
# and repeated runs give independent replicates. GlobalSensitivity.jl uses
# the unscrambled sequence here (deterministic, with the n points that
# follow its skip) and, for Sobol, the first and last d dimensions of a
# 2 d nboot dimensional sequence; the scrambled design is the change of the
# port (see PORTING.md). The deterministic design stays available through
# gsa(f, Sobol(), A, B) with generate_design_matrices().
.gsa_rqmc <- function(n, lb, ub, rng) {
  # Owen scrambling works on 2^m points; every scrambled point is uniform on
  # its own, so the first n of the next power of 2 keep that property
  N <- 2^ceiling(log2(max(n, 2)))
  U <- randomize(qmc_sample(N, length(lb), SobolSample()), OwenScramble(2L, rng = rng))
  lb + (ub - lb) * U[, seq_len(n), drop = FALSE]
}

.gsa_need <- function(samples) if (is.null(samples)) stop("`samples` is required", call. = FALSE)

.gsa_bounds <- function(p_range) {
  if (is.matrix(p_range) && ncol(p_range) == 2L) p_range <- lapply(seq_len(nrow(p_range)), function(i) p_range[i, ])
  list(lb = vapply(p_range, function(r) as.numeric(r[[1L]]), 0), ub = vapply(p_range, function(r) as.numeric(r[[2L]]), 0))
}

# A model output as a numeric vector, with the shape of a solution or array
.gsa_out <- function(y) {
  if (inherits(y, "SciMLSolution")) {
    m <- do.call(cbind, lapply(y$u, as.numeric))
    return(list(v = as.numeric(m), size = dim(m)))
  }
  if (is.array(y) && length(dim(y)) >= 2L) return(list(v = as.numeric(y), size = dim(y)))
  list(v = as.numeric(y), size = NULL)
}

# The outputs of f on the columns of X: a vector (one output) or a matrix
# with one column per parameter set; y_size is the shape of a solution
.gsa_eval <- function(f, X, batch) {
  if (batch) {
    Y <- f(X)
    return(list(Y = if (is.matrix(Y)) Y else as.numeric(Y), multi = is.matrix(Y), y_size = NULL))
  }
  ys <- lapply(seq_len(ncol(X)), function(j) .gsa_out(f(X[, j])))
  y_size <- ys[[1L]]$size
  multi <- !is.null(y_size) || length(ys[[1L]]$v) != 1L
  Y <- if (multi) do.call(cbind, lapply(ys, `[[`, "v")) else vapply(ys, `[[`, 0, "v")
  list(Y = Y, multi = multi, y_size = y_size)
}

# reshape the columns of x (outputs x parameters) to the output shape
.gsa_shape <- function(x, y_size) {
  if (is.null(y_size)) return(x)
  lapply(seq_len(ncol(x)), function(i) array(x[, i], y_size))
}

# std(x; mean) of Julia
.gsa_std_mean <- function(x, m) sqrt(sum((x - m)^2) / (length(x) - 1L))

# ---------------------------------------------------------------------------
# Sobol

.gsa_fuse <- function(A, B, second_order) {
  d <- nrow(A)
  Ab <- lapply(seq_len(d), function(i) { M <- A; M[i, ] <- B[i, ]; M })
  out <- cbind(A, B, do.call(cbind, Ab))
  if (second_order) {
    Ba <- lapply(seq_len(d), function(i) { M <- B; M[i, ] <- A[i, ]; M })
    out <- cbind(out, do.call(cbind, Ba))
  }
  out
}

.gsa_sobol_AB <- function(f, method, A, B, batch, Ei_estimator) {
  A <- as.matrix(A); B <- as.matrix(B)
  d <- nrow(A)
  nboot <- method$nboot
  n <- ncol(A) %/% nboot
  so <- 2L %in% method$order
  pts <- do.call(cbind, lapply(seq_len(nboot), function(i) {
    cols <- (n * (i - 1L) + 1L):(n * i)
    .gsa_fuse(A[, cols, drop = FALSE], B[, cols, drop = FALSE], so)
  }))
  ev <- .gsa_eval(f, pts, batch)
  .gsa_sobol_analysis(method, ev$Y, d, n, Ei_estimator, ev$y_size, ev$multi)
}

.gsa_sobol_analysis <- function(method, all_y, d, n, Ei_estimator, y_size, multi) {
  nboot <- method$nboot
  so <- 2L %in% method$order
  step <- if (so) 2L * d + 2L else d + 2L
  Ym <- if (multi) all_y else matrix(all_y, 1L)
  nout <- nrow(Ym)
  blk <- function(j) Ym[, (j * n + 1L):((j + 1L) * n), drop = FALSE]
  rs <- function(M) rowSums(M)
  est <- match.arg(Ei_estimator, c("Jansen1999", "Homma1996", "Sobol2007", "Janon2014"))
  Varys <- Vis <- Eis <- Vijs <- vector("list", nboot)
  b <- 0L
  for (i in seq.int(1L, step * nboot, by = step)) {
    b <- b + 1L
    AB <- Ym[, ((i - 1L) * n + 1L):((i + 1L) * n), drop = FALSE]
    Varys[[b]] <- apply(AB, 1L, stats::var)
    fA <- blk(i - 1L); fB <- blk(i)
    fAi <- lapply((i + 1L):(i + d), blk)
    if (so) fBi <- lapply((i + d + 1L):(i + 2L * d), blk)
    Vis[[b]] <- vapply(seq_len(d), function(k) rs(fB * (fAi[[k]] - fA)), numeric(nout)) / n
    if (so) {
      M <- array(0, c(d, d, nout))
      for (k in seq_len(d)) for (j in seq_len(d)) if (j > k) M[k, j, ] <- rs((fBi[[k]] * fAi[[j]]) - (fA * fB)) / n
      Vijs[[b]] <- M
    }
    Eis[[b]] <- vapply(seq_len(d), function(k) {
      switch(est,
        Homma1996 = rs((fA - rs(fA) / n)^2) / (n - 1) - rs(fA * fAi[[k]]) / n + (rs(fA) / n)^2,
        Sobol2007 = rs(fA * (fA - fAi[[k]])) / n,
        Jansen1999 = rs((fA - fAi[[k]])^2) / (2 * n),
        Janon2014 = {
          fk <- fAi[[k]]
          (rs(fA^2 + fk^2) / (2 * n) - (rs(fA + fk) / (2 * n))^2) *
            (1 - (1 / n * rs(fA * fk) - (1 / n * rs((fA + fk) / 2))^2) /
                 (1 / n * rs((fA^2 + fk^2) / 2) - (1 / n * rs((fA + fk) / 2))^2))
        })
    }, numeric(nout))
  }
  # nout x d matrices (vapply gives d columns; a single output gives a vector)
  as_mat <- function(v) matrix(v, nout, d)
  Sis <- lapply(seq_len(nboot), function(b) as_mat(Vis[[b]]) / Varys[[b]])
  Tis <- lapply(seq_len(nboot), function(b) as_mat(Eis[[b]]) / Varys[[b]])
  if (so) {
    Sijs <- lapply(seq_len(nboot), function(b) {
      Vi <- as_mat(Vis[[b]])
      out <- Vijs[[b]]
      for (l in seq_len(nout)) {
        M <- matrix(0, d, d)
        for (k in seq_len(d)) for (j in seq_len(d)) if (j > k) M[k, j] <- Vi[l, k] + Vi[l, j]
        out[, , l] <- (Vijs[[b]][, , l] - M) / Varys[[b]][l]
      }
      out
    })
  }
  S1_CI <- ST_CI <- S2_CI <- NULL
  if (nboot > 1L) {
    z <- stats::qnorm((1 + method$conf_level) / 2)
    agg <- function(L) {
      arr <- simplify2array(L)
      dm <- dim(L[[1L]])
      k <- length(dm) + 1L
      mn <- apply(arr, seq_len(k - 1L), mean)
      ci <- apply(arr, seq_len(k - 1L), function(x) z * .gsa_std_mean(x, mean(x)) / sqrt(length(x)))
      list(mean = array(mn, dm), ci = array(ci, dm))
    }
    a1 <- agg(Sis); aT <- agg(Tis)
    Si <- a1$mean; Ti <- aT$mean; S1_CI <- a1$ci; ST_CI <- aT$ci
    if (so) { a2 <- agg(Sijs); Sij <- a2$mean; S2_CI <- a2$ci }
  } else {
    Si <- Sis[[1L]]; Ti <- Tis[[1L]]
    if (so) Sij <- Sijs[[1L]]
  }
  simp <- function(x) if (!multi && !is.null(x)) {
    if (is.array(x) && length(dim(x)) == 3L) x[, , 1L] else as.numeric(x)
  } else x
  if (!is.null(y_size)) { Si <- .gsa_shape(Si, y_size); Ti <- .gsa_shape(Ti, y_size) }
  structure(list(S1 = if (is.null(y_size)) simp(Si) else Si, S1_Conf_Int = simp(S1_CI),
                 S2 = if (so) simp(Sij), S2_Conf_Int = simp(S2_CI),
                 ST = if (is.null(y_size)) simp(Ti) else Ti, ST_Conf_Int = simp(ST_CI)),
            class = "SobolResult")
}

# ---------------------------------------------------------------------------
# Morris

.gsa_morris_design <- function(lb, ub, p_steps, rng, len) {
  k <- length(lb)
  ps <- lapply(seq_len(k), function(i) .jl_linrange(lb[i], ub[i], p_steps[i]))
  idx <- vapply(p_steps, function(i) .qmc_int(rng, 1, 1, i), 0)
  M <- matrix(0, k, len)
  for (i in seq_len(len)) {
    j <- .qmc_int(rng, 1, 1, k)
    idx[j] <- idx[j] + (if (.qmc_unif(rng, 1) < 0.5) -1 else 1)
    if (idx[j] > p_steps[j]) idx[j] <- idx[j] - 2 else if (idx[j] < 1) idx[j] <- idx[j] + 2
    M[, i] <- vapply(seq_len(k), function(u) ps[[u]][idx[u]], 0)
  }
  M
}

.gsa_morris <- function(f, method, p_range, batch, rng) {
  r <- .gsa_bounds(p_range)
  k <- length(r$lb)
  p_steps <- method$p_steps
  if (length(p_steps) != k) p_steps <- c(p_steps, rep(100L, k - length(p_steps)))
  nt <- method$num_trajectory; len <- method$len_design_mat
  if (method$total_num_trajectory < nt)
    stop("total_num_trajectory should be greater than num_trajectory preferably atleast 3-4 times higher", call. = FALSE)
  mats <- lapply(seq_len(method$total_num_trajectory), function(i) .gsa_morris_design(r$lb, r$ub, p_steps, rng, len))
  spread <- vapply(mats, function(m) {
    s <- 0
    for (i in seq_len(ncol(m))[-1L]) s <- s + sqrt(sum((m[, i] - m[, i - 1L])^2))
    s
  }, 0)
  dm <- do.call(cbind, mats[order(spread, decreasing = TRUE)][seq_len(nt)])
  ev <- .gsa_eval(f, dm, batch)
  Y <- if (ev$multi) ev$Y else matrix(ev$Y, 1L)
  effects <- vector("list", k)
  for (i in seq_len(nt)) {
    y1 <- Y[, (i - 1L) * len + 1L]
    for (j in ((i - 1L) * len + 1L):(i * len - 1L)) {
      y2 <- y1
      del <- dm[, j + 1L] - dm[, j]
      ci <- which(abs(del) > 0)[1L]
      del <- sum(del)
      y1 <- Y[, j + 1L]
      e <- if (!method$relative_scale) (y1 - y2) / del else if (del > 0) (y1 - y2) / (y2 * del) else (y1 - y2) / (y1 * del)
      if (!is.na(ci)) effects[[ci]] <- c(effects[[ci]], list(e))
    }
  }
  nout <- nrow(Y)
  stat <- function(fun) vapply(effects, function(e) if (length(e)) fun(do.call(cbind, e)) else numeric(nout), numeric(nout))
  means <- stat(function(E) rowMeans(E))
  means_star <- stat(function(E) rowMeans(abs(E)))
  variances <- stat(function(E) apply(E, 1L, stats::var))
  mk <- function(x) {
    x <- matrix(x, nout, k)
    if (!is.null(ev$y_size)) lapply(seq_len(k), function(i) array(x[, i], ev$y_size)) else x
  }
  ee <- lapply(effects, function(e) if (ev$multi) e else unlist(e))
  structure(list(means = mk(means), means_star = mk(means_star), variances = mk(variances), elementary_effects = ee),
            class = "MorrisResult")
}

# ---------------------------------------------------------------------------
# Regression

# pinv of Julia: singular values below eps * min(m, n) * sigma_max are dropped
.gsa_pinv <- function(A) {
  s <- svd(A)
  tol <- .Machine$double.eps * min(dim(A)) * max(s$d)
  keep <- s$d > tol
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

.gsa_rowsd <- function(M) apply(M, 1L, stats::sd)

.gsa_reg_coefs <- function(X, Y) {
  beta <- qr.coef(qr(t(X)), t(Y))
  t(beta * .gsa_rowsd(X) / rep(.gsa_rowsd(Y), each = nrow(X)))
}
.gsa_reg_corr <- function(X, Y) t(stats::cov(t(X), t(Y)) / outer(.gsa_rowsd(X), .gsa_rowsd(Y)))
.gsa_reg_partial <- function(X, Y) {
  XY <- rbind(X, Y)
  sdv <- .gsa_rowsd(XY)
  prec <- .gsa_pinv(stats::cov(t(XY)) / outer(sdv, sdv))
  pcc <- -prec / sqrt(outer(diag(prec), diag(prec)))
  t(pcc[seq_len(nrow(X)), nrow(X) + seq_len(nrow(Y)), drop = FALSE])
}

.gsa_regression <- function(X, Y, method) {
  X <- as.matrix(X)
  Y <- if (is.matrix(Y)) Y else matrix(as.numeric(Y), 1L)
  res <- list(pearson = .gsa_reg_corr(X, Y), standard_regression = .gsa_reg_coefs(X, Y),
              partial_correlation = .gsa_reg_partial(X, Y),
              pearson_rank = NULL, standard_rank_regression = NULL, partial_rank_correlation = NULL)
  if (method$rank) {
    # ranks of the samples (GlobalSensitivity.jl uses sortperm, the inverse
    # permutation, here; see PORTING.md)
    rk <- function(M) t(apply(M, 1L, rank, ties.method = "first"))
    Xr <- rk(X); Yr <- rk(Y)
    if (nrow(Y) == 1L) Yr <- matrix(Yr, 1L)
    if (nrow(X) == 1L) Xr <- matrix(Xr, 1L)
    res$pearson_rank <- .gsa_reg_corr(Xr, Yr)
    res$standard_rank_regression <- .gsa_reg_coefs(Xr, Yr)
    res$partial_rank_correlation <- .gsa_reg_partial(Xr, Yr)
  }
  structure(res, class = "RegressionGSAResult")
}

# ---------------------------------------------------------------------------
# DGSM

.gsa_grad <- function(f, x) {
  r <- .ad_try_jacobian(f, x)
  if (!is.null(r)) return(as.numeric(r$jacobian))
  as.numeric(.fd_jacobian(function(z) f(z), x))
}
.gsa_hess <- function(f, x) {
  H <- tryCatch(.ad_hessian(f, x), error = function(e) NULL)
  if (!is.null(H) && !anyNA(H)) return(H)
  .fd_jacobian(function(z) .gsa_grad(f, z), x)
}

.gsa_dgsm <- function(f, method, distr, samples) {
  k <- length(distr)
  XX <- lapply(seq_len(samples), function(s) vapply(distr, function(D) .dist_r(D, 1), 0))
  dfdx <- do.call(rbind, lapply(XX, function(x) .gsa_grad(f, x)))
  dfdx <- matrix(dfdx, samples, k)
  a <- colMeans(dfdx); asq <- colMeans(dfdx^2); absa <- colMeans(abs(dfdx))
  # tau_i = E[(df/du_i)^2 (1 - 3 u_i + 3 u_i^2) / 6] and
  # sigma_i = E[(df/du_i)^2 u_i (1 - u_i) / 2] of Kucherenko and Iooss (2016,
  # eq. 32 and UB2) for an input uniform on [0, 1]; a uniform input on
  # [lo, hi] is rescaled to u = (x - lo) / (hi - lo), and the measures are NA
  # for other laws, where these bounds do not hold. GlobalSensitivity.jl
  # weights tau with 1 - 3 x + x^2 and uses the raw inputs (see PORTING.md).
  X <- matrix(unlist(XX), samples, k, byrow = TRUE)
  tao <- sigma <- rep(NA_real_, k)
  for (i in seq_len(k)) {
    D <- distr[[i]]
    if (!identical(D$family, "unif")) next
    lo <- .dist_q(D, 0); hi <- .dist_q(D, 1)
    u <- (X[, i] - lo) / (hi - lo)
    g2 <- (dfdx[, i] * (hi - lo))^2
    tao[i] <- mean(g2 * (1 - 3 * u + 3 * u^2) / 6)
    sigma[i] <- mean(g2 * u * (1 - u) / 2)
  }
  res <- list(a = a, absa = absa, asq = asq, sigma = sigma, tao = tao, crossed = NULL, abscrossed = NULL, crossedsq = NULL)
  if (method$crossed) {
    H <- lapply(XX, function(x) .gsa_hess(f, x))
    cr <- crsq <- abscr <- matrix(0, k, k)
    for (a1 in seq_len(k)) for (b1 in seq_len(k)) if (b1 > a1) {
      v <- vapply(H, function(h) h[b1, a1], 0)
      cr[b1, a1] <- cr[a1, b1] <- mean(v)
      crsq[b1, a1] <- crsq[a1, b1] <- mean(v^2)
      abscr[b1, a1] <- abscr[a1, b1] <- mean(abs(v))
    }
    res$crossed <- cr; res$abscrossed <- abscr; res$crossedsq <- crsq
  }
  structure(res, class = "DGSMResult")
}

# ---------------------------------------------------------------------------
# eFAST

.gsa_efast <- function(f, method, p_range, samples, batch, rng) {
  nh <- method$num_harmonics
  np <- length(p_range)
  omega1 <- (samples - 1) %/% (2 * nh)
  m <- omega1 %/% (2 * nh)
  isdist <- vapply(p_range, function(r) inherits(r, "Distribution"), NA)
  if (m >= np - 1) {
    omega <- c(omega1, floor(.jl_linrange(1, m, np - 1)))
  } else {
    if (m == 0) stop("eFAST needs more samples for this number of harmonics and parameters", call. = FALSE)
    omega <- c(omega1, (seq_len(np - 1) - 1) %% m + 1)
  }
  s <- 2 * (seq_len(samples) - 1) / samples
  ps <- matrix(0, np, samples * np)
  for (i in seq_len(np)) {
    om <- numeric(np)
    om[i] <- omega[1L]
    for (k in seq_len(i - 1L)) om[k] <- omega[k + 1L]
    for (k in seq_len(np)[-seq_len(i)]) om[k] <- omega[k]
    l <- ((i - 1L) * samples + 1L):(i * samples)
    phi <- 2 * .qmc_unif(rng, 1)
    for (j in seq_len(np)) {
      u <- 0.5 + (1 / pi) * asin(sinpi(om[j] * s + phi))
      if (!isdist[j] && p_range[[j]][[1L]] == p_range[[j]][[2L]]) ps[j, l] <- p_range[[j]][[1L]]
      else if (isdist[j]) ps[j, l] <- .dist_q(p_range[[j]], u)
      else ps[j, l] <- p_range[[j]][[1L]] + u * (p_range[[j]][[2L]] - p_range[[j]][[1L]])
    }
  }
  ev <- .gsa_eval(f, ps, batch)
  Y <- if (ev$multi) ev$Y else matrix(ev$Y, 1L)
  nout <- nrow(Y)
  # the positive frequencies below Nyquist; GlobalSensitivity.jl stops one
  # short for odd sample sizes, which drops a frequency of the spectrum and
  # puts the last harmonic out of bounds (see PORTING.md)
  top <- (samples - 1) %/% 2 + 1
  S1 <- ST <- matrix(0, nout, np)
  for (i in seq_len(np)) for (j in seq_len(nout)) {
    ft <- stats::fft(Y[j, ((i - 1L) * samples + 1L):(i * samples)])
    buf <- Mod(ft[2:top])^2
    z <- sum(buf)
    S1[j, i] <- sum(buf[seq_len(nh) * omega1]) / z
    ST[j, i] <- 1 - sum(buf[seq_len(omega1 %/% 2)]) / z
  }
  if (!is.null(ev$y_size)) {
    S1 <- do.call(cbind, lapply(seq_len(np), function(i) array(S1[, i], ev$y_size)))
    ST <- do.call(cbind, lapply(seq_len(np), function(i) array(ST[, i], ev$y_size)))
  }
  structure(list(S1 = S1, ST = ST), class = "eFASTResult")
}

# ---------------------------------------------------------------------------
# RBD-FAST

.gsa_rbdfast <- function(f, method, num_params, samples, rng, batch) {
  if (is.null(num_params)) stop("`num_params` is required", call. = FALSE)
  s0 <- .jl_linrange(-pi, pi, samples)
  s <- lapply(seq_len(num_params), function(i) s0[.qmc_shuffle(rng, seq_len(samples))])
  x <- vapply(s, function(si) 0.5 + asin(sin(si)) / pi, numeric(samples))
  x <- matrix(x, samples)
  Y <- if (batch) as.numeric(f(t(x))) else vapply(seq_len(samples), function(j) as.numeric(f(x[j, ])), 0)
  S <- .gsa_rbdfast_analysis(Y, s, method$num_harmonics)
  # bias correction of Tissot and Prieur (2012), lambda = 2 M / N: the
  # spectrum of the noise adds about lambda (1 - S) to each index
  # (GlobalSensitivity.jl returns S uncorrected; see PORTING.md)
  lam <- 2 * method$num_harmonics / samples
  S - lam / (1 - lam) * (1 - S)
}

.gsa_rbdfast_analysis <- function(Y, s, H) {
  n <- length(Y)
  vapply(s, function(si) {
    ft <- stats::fft(Y[order(si)])[seq_len(n %/% 2 + 1L)]
    a2 <- Mod(ft)^2
    L <- length(a2)
    V <- 2 * sum(a2[2:(L - 1L)]) + a2[L]
    Vi <- 2 * sum(a2[2:(1L + H)])
    Vi / V
  }, 0)
}

# ---------------------------------------------------------------------------
# EASI

.gsa_dct2 <- function(x) {
  # unnormalized DCT-II, sum_n x_n cos(pi k (2n + 1) / (2N)); the indices
  # are ratios of sums of squares over k >= 1, where FFTW's orthonormal
  # scaling is a common factor
  N <- length(x)
  V <- stats::fft(c(x, rev(x)))[seq_len(N)]
  Re(exp(-1i * pi * (seq_len(N) - 1) / (2 * N)) * V) / 2
}

.gsa_easi <- function(X, Y, method) {
  X <- as.matrix(X); Y <- as.numeric(Y)
  K <- nrow(X); n <- ncol(X); mh <- method$max_harmonic
  S <- vapply(seq_len(K), function(i) {
    perm <- order(X[i, ])
    if (method$dct_method) {
      d <- .gsa_dct2(Y[perm])
      sum(d[2:min(1L + mh, n)]^2) / sum(d[2:n]^2)
    } else {
      odd <- perm[seq(1L, n, by = 2L)]
      even <- if (n >= 2L) perm[rev(seq(2L, 2L * (n %/% 2L), by = 2L))] else integer()
      ft <- stats::fft(Y[c(odd, even)])[seq_len(n %/% 2 + 1L)]
      a2 <- Mod(ft)^2
      L <- length(a2)
      sum(a2[2:min(1L + mh, L - 1L)]) / sum(a2[2:(L - 1L)])
    }
  }, 0)
  lam <- (2 * mh) / n
  structure(list(S1 = S, S1_Corr = S - (lam / (1 - lam)) * (1 - S)), class = "EASIResult")
}

# ---------------------------------------------------------------------------
# Fractional factorial

.gsa_hadamard <- function(k) {
  if (k < 2 || bitwAnd(k, k - 1L) != 0) stop("the Hadamard matrix needs a power of 2 of at least 2", call. = FALSE)
  h <- matrix(1, k, k)
  h[2, 2] <- -1
  br <- 2L
  while (br < k) {
    cur <- h[seq_len(br), seq_len(br)]
    h[seq_len(br), br + seq_len(br)] <- cur
    h[br + seq_len(br), seq_len(br)] <- cur
    h[br + seq_len(br), br + seq_len(br)] <- -cur
    br <- 2L * br
  }
  h
}

.gsa_ff <- function(f, num_params, p_range) {
  if (is.null(num_params)) stop("`num_params` is required", call. = FALSE)
  k <- as.integer(round(2^ceiling(log2(num_params))))
  s <- .gsa_hadamard(k)
  design <- rbind(s, -s)
  smp <- design
  if (!is.null(p_range)) for (ci in seq_along(p_range)) {
    smp[design[, ci] == -1, ci] <- p_range[[ci]][[1L]]
    smp[design[, ci] == 1, ci] <- p_range[[ci]][[2L]]
  }
  y <- vapply(seq_len(nrow(smp)), function(i) as.numeric(f(smp[i, ])), 0)
  me <- as.numeric(crossprod(design, y)) / nrow(design)
  list(me, me^2)
}

# ---------------------------------------------------------------------------
# Delta moment independent measure

# KernelDensity.jl's default_bandwidth (Silverman's rule)
.kde_bandwidth <- function(x, alpha = 0.9) {
  n <- length(x)
  if (n <= 1) return(alpha)
  vw <- stats::sd(x)
  q <- stats::quantile(x, c(0.25, 0.75), names = FALSE)
  w <- min(vw, (q[2L] - q[1L]) / 1.34)
  if (w == 0) w <- if (vw == 0) 1 else vw
  alpha * w * n^(-0.2)
}

# solve tridiag(1, 6, 1) x = r: the Thomas algorithm, whose pivots reach the
# constant 3 - 2 sqrt(2) after a few rows, where the sweeps are recursive
# filters
.kde_tri6 <- function(r) {
  m <- length(r)
  if (m == 0L) return(numeric())
  L <- min(m, 40L)
  ci <- 3 - 2 * sqrt(2)
  cp <- numeric(L); dp <- numeric(m)
  cp[1L] <- 1 / 6; dp[1L] <- r[1L] / 6
  for (i in seq_len(L)[-1L]) {
    cp[i] <- 1 / (6 - cp[i - 1L])
    dp[i] <- (r[i] - dp[i - 1L]) * cp[i]
  }
  if (m > L) dp[(L + 1L):m] <- stats::filter(ci * r[(L + 1L):m], -ci, method = "recursive", init = dp[L])
  x <- numeric(m)
  x[m] <- dp[m]
  if (m - 1L > L) {
    idx <- (m - 1L):(L + 1L)
    x[idx] <- stats::filter(dp[idx], -ci, method = "recursive", init = x[m])
  }
  for (i in rev(seq_len(min(L, m - 1L)))) x[i] <- dp[i] - cp[i] * x[i + 1L]
  x
}

# pdf(kde(data), grid): the binned Gaussian kernel density estimate on 2048
# points, interpolated by a quadratic B-spline with Line(OnGrid()) boundary
# conditions (InterpKDE) and zero outside its grid
.kde_pdf <- function(data, grid, npoints = 2048L) {
  bw <- .kde_bandwidth(data)
  if (!(bw > 0)) stop("Bandwidth must be positive", call. = FALSE)
  lo <- min(data) - 4 * bw; hi <- max(data) + 4 * bw
  if (!(lo < hi)) stop("boundary (a,b) must have a < b", call. = FALSE)
  mid <- .jl_linrange(lo, hi, npoints)
  s <- attr(mid, "step")
  mid <- as.numeric(mid)
  K <- npoints
  dens <- numeric(K)
  ainc <- 1 / (s * s)
  w <- 1 / length(data)
  kk <- findInterval(data, mid, left.open = TRUE) + 1L
  for (i in seq_along(data)) {
    k <- kk[i]; j <- k - 1L
    if (j >= 1L && j <= K - 1L) {
      dens[j] <- dens[j] + (mid[k] - data[i]) * ainc * w
      dens[k] <- dens[k] + (data[i] - mid[j]) * ainc * w
    }
  }
  ft <- stats::fft(dens)
  jj <- seq_len(K) - 1
  t <- -2 * pi / (s * K) * pmin(jj, K - jj)
  dens <- pmax(0, Re(stats::fft(ft * exp(-bw^2 / 2 * t^2), inverse = TRUE)) / K)
  # coefficients c_0 .. c_{K+1}: c_1 = y_1, c_K = y_K, interior rows
  # c_{i-1} + 6 c_i + c_{i+1} = 8 y_i, ends extrapolated linearly
  cf <- numeric(K + 2L)
  cf[2L] <- dens[1L]; cf[K + 1L] <- dens[K]
  if (K > 2L) {
    r <- 8 * dens[2:(K - 1L)]
    r[1L] <- r[1L] - dens[1L]
    r[K - 2L] <- r[K - 2L] - dens[K]
    cf[3:K] <- .kde_tri6(r)
  }
  cf[1L] <- 2 * cf[2L] - cf[3L]
  cf[K + 2L] <- 2 * cf[K + 1L] - cf[K]
  u <- (grid - mid[1L]) / s + 1
  inb <- grid >= mid[1L] & grid <= mid[K]
  out <- numeric(length(grid))
  if (any(inb)) {
    ui <- u[inb]
    xm <- ifelse(ui < K + 0.5, floor(ui + 0.5), ceiling(ui + 0.5) - 1)
    dx <- ui - xm
    out[inb] <- (dx - 0.5)^2 / 2 * cf[xm] + (0.75 - dx^2) * cf[xm + 1] + (dx + 0.5)^2 / 2 * cf[xm + 2]
  }
  out
}

# trapz of Trapz.jl
.gsa_trapz <- function(x, y) {
  N <- length(x)
  if (N <= 1L) return(0)
  r <- (x[2L] - x[1L]) * y[1L]
  if (N > 2L) for (i in 2:(N - 1L)) r <- r + (x[i + 1L] - x[i - 1L]) * y[i]
  r <- r + (x[N] - x[N - 1L]) * y[N]
  r / 2
}

.gsa_calc_delta <- function(Xi, Y, Ygrid, cuts) {
  if (all(Y[1L] == Y)) return(0)
  n <- length(Y)
  fy <- .kde_pdf(Y, Ygrid)
  xr <- rank(Xi, ties.method = "min")
  seps <- numeric(length(cuts) - 1L)
  for (j in seq_len(length(cuts) - 1L)) {
    inc <- which((xr > cuts[j]) == (xr <= cuts[j + 1L]))
    if (!length(inc)) next
    fyc <- .kde_pdf(Y[inc], Ygrid)
    seps[j] <- length(inc) * .gsa_trapz(Ygrid, abs(fy - fyc))
  }
  sum(seps) / (2 * n)
}

.gsa_delta <- function(X, Y, method, rng) {
  X <- as.matrix(X)
  n <- ncol(X)
  M <- if (is.null(method$num_classes)) {
    ex <- 2 / (7 + tanh((1500 - n) / 500))
    round(min(ceiling(n^ex), 48))
  } else method$num_classes
  cuts <- as.numeric(.jl_linrange(0, n, M + 1))
  Ygrid <- as.numeric(.jl_linrange(min(Y), max(Y), method$Ygrid_length))
  nf <- nrow(X)
  deltas <- adj <- lo <- hi <- numeric(nf)
  z <- stats::qnorm(0.5 + method$conf_level / 2)
  for (fi in seq_len(nf)) {
    Xi <- X[fi, ]
    deltas[fi] <- .gsa_calc_delta(Xi, Y, Ygrid, cuts)
    r <- matrix(.qmc_int(rng, method$nboot * n, 1, n), method$nboot, n)
    bd <- vapply(seq_len(method$nboot), function(i) .gsa_calc_delta(Xi[r[i, ]], Y[r[i, ]], Ygrid, cuts), 0)
    adj[fi] <- 2 * deltas[fi] - mean(bd)
    band <- z * stats::sd(bd) / sqrt(method$nboot)
    lo[fi] <- adj[fi] - band; hi[fi] <- adj[fi] + band
  }
  structure(list(deltas = deltas, adjusted_deltas = adj, adjusted_deltas_low = lo, adjusted_deltas_hi = hi),
            class = "DeltaResult")
}

# ---------------------------------------------------------------------------
# Shapley effects

.gsa_permutations <- function(v) {
  if (length(v) <= 1L) return(list(v))
  out <- list()
  for (i in seq_along(v)) for (p in .gsa_permutations(v[-i])) out[[length(out) + 1L]] <- c(v[i], p)
  out
}

.gsa_cond_sampling <- function(D, n, idx, idx_c, x_cond) {
  if (inherits(D$C, "IndependentCopula")) return(.dist_rand_mat(.dist_subset(D, idx), n))
  z <- vapply(seq_along(idx_c), function(i) stats::qnorm(.dist_p(D$m[[idx_c[i]]], x_cond[i])), 0)
  S <- D$C$Sigma
  cdinv <- t(S[idx_c, idx, drop = FALSE]) %*% solve(S[idx_c, idx_c, drop = FALSE])
  mu <- as.numeric(cdinv %*% z)
  V <- S[idx, idx, drop = FALSE] - cdinv %*% S[idx_c, idx, drop = FALSE]
  k <- length(idx)
  smp <- if (k == 1L) matrix(stats::rnorm(n, mu, sqrt(V[1L, 1L])), 1L)
         else mu + t(chol((V + t(V)) / 2)) %*% matrix(stats::rnorm(k * n), k, n)
  for (i in seq_len(k)) smp[i, ] <- .dist_q(D$m[[idx[i]]], stats::pnorm(smp[i, ]))
  smp
}

.gsa_shapley <- function(f, method, D, batch) {
  if (!inherits(D, "SklarDist")) stop("Shapley needs a SklarDist input distribution", call. = FALSE)
  dim <- length(D)
  n_var <- method$n_var; n_outer <- method$n_outer; n_inner <- method$n_inner
  if (method$n_perms == -1L) {
    message("Since `n_perms` wasn't set the exact version of Shapley will be used")
    perms <- .gsa_permutations(seq_len(dim))
  } else {
    message("Since `n_perms` was set the random version of Shapley will be used")
    perms <- lapply(seq_len(method$n_perms), function(i) sample.int(dim))
  }
  m <- length(perms)
  A <- .dist_rand_mat(D, n_var)
  Bm <- matrix(0, dim, m * (dim - 1L) * n_outer * n_inner)
  for (ip in seq_len(m)) {
    perm <- perms[[ip]]
    srt <- order(perm)
    for (j in seq_len(dim - 1L)) {
      plus <- perm[seq_len(j)]
      minus <- perm[(j + 1L):dim]
      comp <- .dist_rand_mat(.dist_subset(D, minus), n_outer)
      for (l in seq_len(n_outer)) {
        cur <- comp[, l]
        xj <- .gsa_cond_sampling(D, n_inner, plus, minus, cur)
        xx <- rbind(xj, matrix(cur, length(cur), ncol(xj)))
        ind <- (ip - 1L) * (dim - 1L) * n_outer * n_inner + (j - 1L) * n_outer * n_inner + (l - 1L) * n_inner + 1L
        Bm[, ind:(ind + n_inner - 1L)] <- xx[srt, , drop = FALSE]
      }
    }
  }
  X <- cbind(A, Bm)
  ev <- .gsa_eval(f, X, batch)
  Ym <- if (ev$multi) ev$Y else matrix(ev$Y, 1L)
  nout <- nrow(Ym)
  Sh <- Sh2 <- matrix(0, nout, dim)
  Y <- Ym[, seq_len(n_var), drop = FALSE]
  y <- Ym[, -seq_len(n_var), drop = FALSE]
  VarY <- apply(Y, 1L, stats::var)
  pos <- 0L
  for (p in seq_len(m)) {
    perm <- perms[[p]]
    prevC <- numeric(nout)
    for (j in seq_len(dim)) {
      if (j == dim) {
        Chat <- VarY
        delta <- Chat - prevC
      } else {
        cVar <- vapply(seq_len(n_outer), function(l)
          apply(y[, pos + (l - 1L) * n_inner + seq_len(n_inner), drop = FALSE], 1L, stats::var), numeric(nout))
        cVar <- matrix(cVar, nout)
        pos <- pos + n_outer * n_inner
        Chat <- rowMeans(cVar)
        delta <- Chat - prevC
        Sh2[, perm[j]] <- Sh2[, perm[j]] + rowMeans((cVar - prevC)^2) - delta^2
      }
      Sh[, perm[j]] <- Sh[, perm[j]] + delta
      prevC <- Chat
    }
  }
  Sh <- Sh / m / VarY
  Sh2 <- Sh2 / m / VarY^2
  SE <- sqrt(Sh2 / n_outer)
  lo <- Sh - 1.96 * SE; hi <- Sh + 1.96 * SE
  fin <- function(x) if (!ev$multi) as.numeric(x) else if (!is.null(ev$y_size)) .gsa_shape(x, ev$y_size) else x
  structure(list(shapley_effects = fin(Sh), std_err = fin(SE), CI_lower = fin(lo), CI_upper = fin(hi)),
            class = "ShapleyResult")
}

# ---------------------------------------------------------------------------
# Regional sensitivity analysis

.gsa_rsa_ks <- function(Xi, flag) {
  o <- order(Xi)
  acc <- cumsum(flag[o]); rej <- cumsum(1 - flag[o])
  max(abs(acc / max(acc) - rej / max(rej)))
}

.gsa_rsa <- function(f, method, p_range, samples, batch) {
  r <- .gsa_bounds(p_range)
  d <- length(r$lb); nd <- method$n_dummy_parameters
  X <- qmc_sample(samples, c(r$lb, rep(0, nd)), c(r$ub, rep(1, nd)), LatinHypercubeSample())
  # the model sees its own parameters only (GlobalSensitivity.jl passes the
  # dummies too; see PORTING.md)
  ev <- .gsa_eval(f, X[seq_len(d), , drop = FALSE], batch)
  .gsa_rsa_analysis(X, as.numeric(ev$Y), method)
}

.gsa_rsa_analysis <- function(X, Y, method) {
  nd <- method$n_dummy_parameters
  K <- nrow(X)
  thr <- if (is.function(method$acceptance_threshold)) method$acceptance_threshold(Y) else method$acceptance_threshold
  flag <- as.numeric(Y > thr)
  sens <- vapply(seq_len(K), function(i) .gsa_rsa_ks(X[i, ], flag), 0)
  dummies <- sens[(K - nd + 1L):K]
  structure(list(S = sens[seq_len(K - nd)], Sd = c(mean(dummies), stats::sd(dummies))), class = "RSAResult")
}

# ---------------------------------------------------------------------------
# Mutual information

# Shannon entropy (base 2) of ComplexityMeasures' ValueBinning(nbins): nbins
# equal bins from the minimum to two floats past the maximum per dimension
.gsa_nextfloat2 <- function(x) .ojl_nextfloat(.ojl_nextfloat(x))
.gsa_entropy <- function(M, nbins) {
  M <- as.matrix(M)
  if (nrow(M) == 1L && ncol(M) > 1L) M <- t(M)
  lin <- rep(0, nrow(M)); mult <- 1
  ok <- rep(TRUE, nrow(M))
  for (j in seq_len(ncol(M))) {
    x <- M[, j]
    mini <- min(x)
    r <- .jl_linrange(mini, .gsa_nextfloat2(max(x)), nbins + 1L)
    b <- floor((x - r[1L]) / attr(r, "step")) + 1
    ok <- ok & b >= 1 & b <= nbins
    lin <- lin + (b - 1) * mult
    mult <- mult * nbins
  }
  cts <- as.numeric(table(lin[ok]))
  p <- cts / sum(cts)
  p <- p[p != 0]
  s <- 0
  for (v in p) s <- s + v * log2(v)
  -s
}

.gsa_mi_bounds <- function(Xi, Y, conf, nb, rng) {
  nbins <- round(sqrt(length(Y)))
  hX <- .gsa_entropy(Xi, nbins); hY <- .gsa_entropy(Y, nbins)
  Yp <- Y
  mi <- numeric(nb)
  for (i in seq_len(nb)) {
    Yp <- .qmc_shuffle(rng, Yp)
    mi[i] <- hX + hY - .gsa_entropy(cbind(Xi, Yp), nbins)
  }
  stats::quantile(mi, conf, names = FALSE)
}

.gsa_mi <- function(f, method, p_range, samples, batch, rng) {
  r <- .gsa_bounds(p_range)
  X <- .gsa_rqmc(samples, r$lb, r$ub, rng)
  .gsa_mi_analysis(X, as.numeric(.gsa_eval(f, X, batch)$Y), method, rng)
}

.gsa_mi_analysis <- function(X, Y, method, rng) {
  K <- nrow(X)
  if (method$n_bootstraps > ncol(X))
    stop("Number of bootstraps must be less than or equal to the number of samples", call. = FALSE)
  nbins <- round(sqrt(length(Y)))
  hY <- .gsa_entropy(Y, nbins)
  mi <- bounds <- numeric(K)
  for (i in seq_len(K)) {
    Xi <- X[i, ]
    mi[i] <- .gsa_entropy(Xi, nbins) + hY - .gsa_entropy(cbind(Xi, Y), nbins)
    bounds[i] <- .gsa_mi_bounds(Xi, Y, method$conf_level, method$n_bootstraps, rng)
  }
  structure(list(S = pmax(0, mi - bounds), mutual_information = mi, bounds = bounds), class = "MutualInformationResult")
}
