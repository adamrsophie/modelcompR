# Univariate distributions and the Sklar copula distributions of Copulas.jl
# (IndependentCopula and GaussianCopula) that GlobalSensitivity.jl samples
# from. A univariate distribution wraps the r, p and q functions of an R
# family; Copulas.jl is MIT licensed, Copyright (c) 2022 Oskar Laverny and
# contributors; see inst/COPYRIGHTS.

#' Distributions for global sensitivity analysis
#'
#' `univariate_distribution()` describes the law of one input, as the
#' distributions of Distributions.jl do in GlobalSensitivity.jl. Give a family
#' of the stats package by its suffix (`"unif"`, `"norm"`, `"gamma"`, ...) and
#' its parameters, or the random generator `r(n)`, the distribution function
#' `p(x)` and the quantile function `q(u)` directly. `SklarDist()` joins
#' marginal distributions with a copula: `IndependentCopula(d)` for
#' independent inputs or `GaussianCopula(Sigma)` with the correlation matrix
#' `Sigma`, as `SklarDist` of Copulas.jl.
#'
#' @param family Suffix of a stats distribution family, or `NULL`.
#' @param ... Parameters of the family, such as `min` and `max`.
#' @param r,p,q Random generator, distribution function and quantile function;
#'   they override the family.
#' @param d Dimension of the copula.
#' @param Sigma Correlation matrix of the Gaussian copula.
#' @param copula A copula.
#' @param marginals List of univariate distributions, one per input.
#' @return A distribution object.
#' @examples
#' x1 <- univariate_distribution("unif", min = -pi, max = pi)
#' x2 <- univariate_distribution("norm", mean = 0, sd = 2)
#' D <- SklarDist(GaussianCopula(matrix(c(1, 0.5, 0.5, 1), 2)), list(x1, x2))
#' set.seed(1)
#' dist_rand(D, 3)
#' @name gsa_distributions
NULL

#' @rdname gsa_distributions
#' @export
univariate_distribution <- function(family = NULL, ..., r = NULL, p = NULL, q = NULL) {
  args <- list(...)
  fam <- function(prefix) {
    if (is.null(family)) return(NULL)
    fn <- get(paste0(prefix, family), envir = asNamespace("stats"))
    function(x) do.call(fn, c(list(x), args))
  }
  structure(list(family = family, params = args, r = r %||% fam("r"), p = p %||% fam("p"), q = q %||% fam("q")),
            class = c("UnivariateDistribution", "Distribution"))
}

#' @rdname gsa_distributions
#' @export
IndependentCopula <- function(d) structure(list(d = as.integer(d)), class = c("IndependentCopula", "Copula"))

#' @rdname gsa_distributions
#' @export
GaussianCopula <- function(Sigma) {
  Sigma <- as.matrix(Sigma)
  if (nrow(Sigma) != ncol(Sigma)) stop("`Sigma` must be a square matrix", call. = FALSE)
  structure(list(d = nrow(Sigma), Sigma = Sigma), class = c("GaussianCopula", "Copula"))
}

#' @rdname gsa_distributions
#' @export
SklarDist <- function(copula, marginals) {
  if (!inherits(copula, "Copula")) stop("`copula` must be a copula", call. = FALSE)
  if (inherits(marginals, "UnivariateDistribution")) marginals <- list(marginals)
  if (length(marginals) != copula$d) stop("one marginal distribution per dimension of the copula is needed", call. = FALSE)
  structure(list(C = copula, m = marginals), class = c("SklarDist", "Distribution"))
}

#' @export
length.SklarDist <- function(x) length(x$m)

.dist_need <- function(D, what) {
  if (is.null(D[[what]])) stop("the distribution has no `", what, "` function", call. = FALSE)
  D[[what]]
}
.dist_r <- function(D, n) as.numeric(.dist_need(D, "r")(n))
.dist_p <- function(D, x) as.numeric(.dist_need(D, "p")(x))
.dist_q <- function(D, u) as.numeric(.dist_need(D, "q")(u))

#' @param D A distribution.
#' @param n Number of draws.
#' @return `dist_rand()` returns `n` draws of a univariate distribution, or a
#'   d x n matrix for a `SklarDist`.
#' @rdname gsa_distributions
#' @export
dist_rand <- function(D, n) {
  if (inherits(D, "UnivariateDistribution")) return(.dist_r(D, n))
  if (!inherits(D, "SklarDist")) stop("`D` must be a distribution", call. = FALSE)
  U <- .copula_rand(D$C, n)
  for (i in seq_along(D$m)) U[i, ] <- .dist_q(D$m[[i]], U[i, ])
  U
}

# d x n uniform draws of a copula
.copula_rand <- function(C, n) {
  if (inherits(C, "IndependentCopula")) return(matrix(stats::runif(C$d * n), C$d, n))
  L <- t(chol(C$Sigma))
  stats::pnorm(L %*% matrix(stats::rnorm(C$d * n), C$d, n))
}

# subsetdims: the law of the inputs idx, univariate for a single input
.dist_subset <- function(D, idx) {
  if (length(idx) == 1L) return(D$m[[idx]])
  C <- if (inherits(D$C, "IndependentCopula")) IndependentCopula(length(idx)) else GaussianCopula(D$C$Sigma[idx, idx, drop = FALSE])
  SklarDist(C, D$m[idx])
}

# rand of a distribution as a matrix with one row per dimension
.dist_rand_mat <- function(D, n) {
  x <- dist_rand(D, n)
  if (is.matrix(x)) x else matrix(x, 1L)
}
