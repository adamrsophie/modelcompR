# Covariance matrix adaptation evolution strategies: the
# OptimizationCMAEvolutionStrategy wrapper of Optimization.jl with the
# CMAEvolutionStrategy.jl (0.2) algorithm, and the OptimizationGCMAES
# wrapper with the GCMAES.jl (0.1) algorithm, ported to R.

# _check_opt_alg of OptimizationBase for the algorithm families of these
# files: bounds, constraints and callbacks against the traits.
.opt_check_traits <- function(a, alg) {
  t <- .opt_traits(alg)
  cl <- class(alg)[1L]
  has_bounds <- !is.null(a$lb) || !is.null(a$ub)
  if (!t$bounds && has_bounds)
    stop(IncompatibleOptimizerError(paste0("The algorithm ", cl, " does not support box constraints. Either remove ",
                                           "the `lb` or `ub` bounds passed to `OptimizationProblem` or use a ",
                                           "different algorithm.")))
  if (t$rbounds && is.null(a$lb))
    stop(IncompatibleOptimizerError(paste0("The algorithm ", cl, " requires box constraints. Either pass `lb` and ",
                                           "`ub` bounds to `OptimizationProblem` or use a different algorithm.")))
  if (!t$cons && is.function(a$f$cons))
    stop(IncompatibleOptimizerError(paste0("The algorithm ", cl, " does not support constraints. Either remove the ",
                                           "`cons` function passed to `OptimizationFunction` or use a different ",
                                           "algorithm.")))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# CMAEvolutionStrategy.jl

#' CMA-ES of CMAEvolutionStrategy.jl
#'
#' `CMAEvolutionStrategyOpt()` minimizes with the covariance matrix
#' adaptation evolution strategy of CMAEvolutionStrategy.jl, ported to R:
#' active covariance updates with negative recombination weights,
#' cumulative step size adaptation and, for box bounds, the
#' linear-quadratic transformation that keeps every evaluated point inside
#' `lb` and `ub`. No derivatives are used.
#'
#' The `solve` arguments map as in OptimizationCMAEvolutionStrategy:
#' `maxiters` to the maximal number of generations, `maxtime` (seconds),
#' `abstol` to `ftol` (the range of the last three generations' values,
#' default 1e-11) and `sigma0`, the initial step size (0.1). Further
#' arguments of `CMAEvolutionStrategy.minimize` are passed on: `popsize`
#' (default `4 + floor(3 log(n))`), `maxfevals`, `stagnation`, `ftarget`,
#' `xtol`, `seed` (the random numbers are R's) and `noise_handling` (`TRUE`
#' or a list of `r`, `epsilon`, `theta`, `c`, `alpha_sigma` and
#' `callback`, the uncertainty handling of Hansen et al. 2009). The return
#' code is `Success` for `ftarget`, `xtol` and `ftol`, `MaxIters` for
#' `maxiter` and `maxfevals`, `MaxTime` and `Stalled` for stagnation. The
#' solution's `original` holds the optimizer: `p` (with `lambda`, the
#' population size, and `mean`), `logger`, `stop` (with `reason`), `xbest`
#' and `fbest`.
#'
#' @return An algorithm object.
#' @examples
#' rosenbrock <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
#' prob <- OptimizationProblem(rosenbrock, c(0, 0), c(1, 100), lb = c(-1, -1), ub = c(0.8, 0.8))
#' set.seed(1)
#' sol <- solve(prob, CMAEvolutionStrategyOpt())
#' round(sol$u, 3); sol$retcode
#' @export
CMAEvolutionStrategyOpt <- function()
  structure(list(name = "CMAEvolutionStrategyOpt", traits = list(bounds = TRUE, callback = TRUE)),
            class = c("CMAEvolutionStrategyOpt", "OptimizationAlgorithm"))

.cmaes_box <- function(lb, ub) {
  if (is.null(lb) && is.null(ub)) return(NULL)
  n <- length(lb %||% ub)
  lb <- as.numeric(lb %||% rep(-Inf, n)); ub <- as.numeric(ub %||% rep(Inf, n))
  al <- ifelse(is.finite(lb), pmin((ub - lb) / 2, (1 + abs(lb)) / 20), 1)
  au <- ifelse(is.finite(ub), pmin((ub - lb) / 2, (1 + abs(ub)) / 20), 1)
  list(lb = lb, ub = ub, al = al, au = au)
}

.cmaes_transform <- function(b, x) {
  if (is.null(b)) return(x)
  for (i in seq_along(x)) {
    lb <- b$lb[i]; ub <- b$ub[i]; al <- b$al[i]; au <- b$au[i]; y <- x[i]
    if (y < lb - 2 * al - (ub - lb) / 2 || y > ub + 2 * au + (ub - lb) / 2) {
      r <- 2 * (ub - lb + al + au)
      s <- lb - 2 * al - (ub - lb) / 2
      y <- y - r * floor((y - s) / r)
    }
    if (y > ub + au) y <- y - 2 * (y - ub - au)
    if (y < lb - al) y <- y + 2 * (lb - al - y)
    x[i] <- if (y < lb + al) lb + (y - (lb - al))^2 / 4 / al else if (y < ub - au) y else
      if (y < ub + 3 * au) ub - (y - (ub + au))^2 / 4 / au else ub + au - (y - (ub + au))
  }
  x
}

.cmaes_backtransform <- function(b, y) {
  if (is.null(b)) return(y)
  for (i in seq_along(y)) {
    lb <- b$lb[i]; ub <- b$ub[i]; al <- b$al[i]; au <- b$au[i]; v <- y[i]
    if (!(lb <= v && v <= ub)) stop(v, " is not within [", lb, ", ", ub, "].", call. = FALSE)
    y[i] <- if (v < lb + al) (lb - al) + 2 * sqrt(al * (v - lb)) else if (v < ub - au) v else
      (ub + au) - 2 * sqrt(au * (ub - v))
  }
  y
}

.cmaes_meigen <- function(m, delta = 1e-15) {
  repeat {
    e <- eigen((m + t(m)) / 2, symmetric = TRUE)
    if (min(e$values) >= 1e-16) return(list(mat = m, vectors = e$vectors, sqrtvalues = sqrt(e$values)))
    m <- m + diag(delta, nrow(m))
    delta <- 10 * delta
  }
}

.cmaes_mueff <- function(w) sum(w)^2 / sum(w^2)

# The recombination weights (finalized) and the step size and covariance
# parameters of CMAEvolutionStrategy.Parameters
.cmaes_parameters <- function(n, lam) {
  w <- log((lam + 1) / 2) - log(seq_len(lam))
  mu <- sum(w > 0)
  w[seq_len(mu)] <- w[seq_len(mu)] / sum(w[seq_len(mu)])
  neg <- if (mu < lam) (mu + 1L):lam else integer()
  if (length(neg)) w[neg] <- w[neg] / -sum(w[neg])
  mueff <- .cmaes_mueff(w[seq_len(mu)])
  sg <- list(c = (mueff + 2) / (n + mueff + 3))
  sg$d <- 1 + max(0, sqrt((mueff - 1) / (n + 1)) - 1) + sg$c
  sg$e <- sqrt(sg$c * (2 - sg$c) * mueff)
  sg$chin <- sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n^2))
  cv <- list(c = (4 + mueff / n) / (n + 4 + 2 * mueff / n), c1 = 2 / ((n + 1.3)^2 + mueff))
  cv$cmu <- min(1 - cv$c1, 2 * (0.25 + mueff - 2 + 1 / mueff) / ((n + 2)^2 + mueff))
  cv$e <- sqrt(cv$c * (2 - cv$c) * mueff)
  if (length(neg)) {
    # finalize_negative_weights!
    w[neg] <- w[neg] * (1 + cv$c1 / cv$cmu)
    lim <- function(v) { sw <- -sum(w[neg]); if (sw > v) w[neg] <<- w[neg] * (v / sw) }
    lim((1 - cv$c1 - cv$cmu) / cv$cmu / n)
    lim(1 + 2 * .cmaes_mueff(w[neg]) / (mueff + 2))
  }
  list(w = w, mu = mu, mueff = mueff, sigma = sg, cov = cv)
}

# CMAEvolutionStrategy.minimize(f, x0, s0; kwargs...)
.cmaes_minimize <- function(f, x0, s0, lower = NULL, upper = NULL, popsize = NULL, maxiter = NULL, maxtime = NULL,
                            maxfevals = NULL, stagnation = NULL, ftarget = NULL, xtol = NULL, ftol = 1e-11,
                            noise_handling = NULL) {
  n <- length(x0)
  lam <- as.integer(popsize %||% (4L + as.integer(floor(3 * log(n)))))
  pr <- .cmaes_parameters(n, lam)
  w <- pr$w; mu <- pr$mu
  sg <- c(pr$sigma, list(sigma = s0, h = TRUE, g = 0L, p = numeric(n)))
  cv <- c(pr$cov, list(p = numeric(n), C = .cmaes_meigen(diag(exp(1e-4 / n * (0:(n - 1L))), n))))
  box <- .cmaes_box(lower, upper)
  mean <- .cmaes_backtransform(box, as.numeric(x0))
  if (isTRUE(noise_handling)) noise_handling <- list(alpha_sigma = 1 + 2 / (n + 10))
  if (is.list(noise_handling))
    noise_handling <- utils::modifyList(list(r = 0.3, epsilon = 0, theta = 0.5, c = 0.3, alpha_sigma = 1.1,
                                             sbar = 0, fevals = 0L, callback = function(s) s > 0), noise_handling)
  stagnation <- stagnation %||% floor(100 + 100 * n^1.5 / lam)
  log <- list(fbest = numeric(), xbest = list(), fmedian = numeric(), frange = numeric())
  it <- 0L
  reason <- "none"
  t0 <- proc.time()[["elapsed"]]
  evaluate <- function(y) vapply(seq_len(ncol(y)), function(i)
    as.numeric(f(.cmaes_transform(box, sg$sigma * y[, i] + mean)))[1L], 0)
  whiten <- function(x) as.numeric(cv$C$vectors %*% (crossprod(cv$C$vectors, x) / cv$C$sqrtvalues))
  repeat {
    y <- cv$C$vectors %*% (cv$C$sqrtvalues * matrix(stats::rnorm(n * lam), n, lam))
    fv <- evaluate(y)
    perm <- order(fv)
    if (!length(log$fbest) || fv[perm[1L]] < log$fbest[length(log$fbest)]) {
      log$fbest <- c(log$fbest, fv[perm[1L]])
      log$xbest <- c(log$xbest, list(mean + sg$sigma * y[, perm[1L]]))
    }
    log$fmedian <- c(log$fmedian, stats::median(fv))
    log$frange <- c(log$frange, max(fv) - min(fv))
    if (is.list(noise_handling)) {
      nh <- .cmaes_noise(noise_handling, lam, fv, y, evaluate)
      noise_handling <- nh$state
      if (nh$inflate) sg$sigma <- sg$sigma * noise_handling$alpha_sigma
      fv <- nh$fvals
      perm <- order(fv)
    }
    # update!(p, y, perm)
    sm <- as.numeric(y[, perm[seq_len(mu)], drop = FALSE] %*% w[seq_len(mu)])
    mean <- mean + sg$sigma * sm
    sg$p <- (1 - sg$c) * sg$p + sg$e * whiten(sm)
    np <- sqrt(sum(sg$p^2))
    sg$h <- np^2 / (1 - (1 - sg$c)^(2 * (sg$g + 1L))) / n - 1 < 1 + 4 / (n + 1)
    sg$sigma <- sg$sigma * exp(min(max(sg$c / sg$d * (np / sg$chin - 1), -1), 1))
    sg$g <- sg$g + 1L
    cv$p <- cv$p * (1 - cv$c)
    if (sg$h) cv$p <- cv$p + cv$e * sm
    wts <- cv$cmu * w
    M <- cv$C$mat * (1 - sg$h * cv$c1 - sum(wts)) + cv$c1 * tcrossprod(cv$p)
    for (i in seq_len(lam)) {
      wi <- wts[i]
      if (wi == 0) next
      v <- y[, perm[i]]
      if (wi < 0) wi <- wi * n / sum((crossprod(cv$C$vectors, v) / cv$C$sqrtvalues)^2)
      M <- M + wi * tcrossprod(v)
    }
    cv$C <- .cmaes_meigen(M)
    # terminate!
    it <- it + 1L
    nf <- length(log$frange)
    reason <- if (!is.null(maxiter) && it >= maxiter) "maxiter" else
      if (!is.null(maxtime) && proc.time()[["elapsed"]] - t0 > maxtime) "maxtime" else
      if (!is.null(maxfevals) && it * lam + (noise_handling$fevals %||% 0L) > maxfevals) "maxfevals" else
      if (!is.null(ftarget) && log$fbest[length(log$fbest)] < ftarget) "ftarget" else
      if (!is.null(xtol) && all(sg$sigma * cv$p < xtol)) "xtol" else
      if (!is.null(ftol) && nf > 10L && all(log$frange[(nf - 2L):nf] < ftol)) "ftol" else
      if (!is.null(stagnation) && it > n * (5 + 100 / lam) && length(log$fmedian) >= 2 * stagnation &&
          stats::median(log$fmedian[(nf - stagnation + 1L):nf]) >=
            stats::median(log$fmedian[(nf - 2L * stagnation + 1L):(nf - stagnation)])) "stagnation" else "none"
    if (reason != "none") break
  }
  xb <- .cmaes_transform(box, log$xbest[[length(log$xbest)]])
  list(p = list(n = n, lambda = lam, mean = mean, population_mean = .cmaes_transform(box, mean), sigma = sg,
                cov = cv, weights = w, mu = mu),
       logger = log, stop = list(it = it, reason = reason), xbest = xb, fbest = log$fbest[length(log$fbest)])
}

# the uncertainty handling of Hansen, Niederberger, Guzzella and Koumoutsakos
# (2009): reevaluate part of the population and inflate the step size when
# the rank changes exceed their limits
.cmaes_noise <- function(nh, lam, fv, y, evaluate) {
  rate <- lam * nh$r
  nre <- as.integer(floor(rate) + (stats::runif(1L) < rate - floor(rate)))
  nh$fevals <- nh$fevals + nre
  if (nre == 0L) return(list(state = nh, inflate = FALSE, fvals = fv))
  if (nh$epsilon != 0) stop("Not implemented.", call. = FALSE)
  newf <- evaluate(y[, seq_len(nre), drop = FALSE])
  joint <- c(fv, newf, fv[-seq_len(nre)])
  rank <- order(order(joint))
  dlim <- function(R) nh$theta * stats::median(abs(seq_len(2L * lam - 1L) - R))
  s <- mean(vapply(seq_len(nre), function(i) {
    d <- rank[i] - rank[lam + i]
    d <- d - sign(d)
    2 * abs(d) - dlim(rank[i] - (newf[i] < fv[i])) - dlim(rank[lam + i] - (fv[i] < newf[i]))
  }, 0))
  nh$sbar <- (1 - nh$c) * nh$sbar + nh$c * s
  fv[seq_len(nre)] <- 0.5 * fv[seq_len(nre)] + 0.5 * newf
  list(state = nh, inflate = isTRUE(nh$callback(nh$sbar)), fvals = fv)
}

.cmaes_retcode <- function(reason)
  switch(reason, ftarget = , xtol = , ftol = "Success", maxiter = , maxfevals = "MaxIters", maxtime = "MaxTime",
         stagnation = "Stalled", "Default")

.mtk$hooks$opt_solve_CMAEvolutionStrategyOpt <- function(a, alg, opts) {
  .opt_check_traits(a, alg)
  p <- a$p
  loss <- function(x) as.numeric(a$f$f(x, p))[1L]
  args <- list(lower = a$lb, upper = a$ub, maxiter = opts$maxiters, maxtime = opts$maxtime,
               popsize = opts$popsize, maxfevals = opts$maxfevals, stagnation = opts$stagnation,
               ftarget = opts$ftarget, xtol = opts$xtol, noise_handling = opts$noise_handling)
  args$ftol <- opts$abstol %||% opts$ftol %||% 1e-11
  args <- args[!vapply(args, is.null, TRUE)]
  r <- .with_seed(opts$seed, do.call(.cmaes_minimize, c(list(loss, as.numeric(a$u0), opts$sigma0 %||% 0.1), args)))
  .opt_solution(a, alg, r$xbest, r$fbest, .cmaes_retcode(r$stop$reason), r)
}

# ---------------------------------------------------------------------------
# GCMAES.jl

#' Gradient assisted CMA-ES of GCMAES.jl
#'
#' `GCMAESOpt()` minimizes with the CMA-ES of GCMAES.jl, ported to R. Each
#' generation samples the population from the covariance matrix, moves the
#' mean to the weighted average of the best half and, when the problem has
#' a gradient, also evaluates the mean itself (the line search of GCMAES
#' along the negative gradient, which in a single process tries only the
#' step length zero). Termination follows GCMAES: flat fitness inflates
#' the step size, and the run stops on the conditions MaxIter, TolHistFun,
#' EqualFunVals, TolX, TolUpSigma, Stagnation, ConditionCov or TolFun,
#' with return code `Success`, or after the evaluation budget with
#' `Failure`.
#'
#' The problem needs bounds (`lb` and `ub`), as in OptimizationGCMAES,
#' although GCMAES only keeps them and searches the whole space. The
#' `solve` arguments are `maxiters` (generations; without it the budget is
#' `1000 n^2` evaluations), `maxtime` (seconds) and `sigma0`, the initial
#' step size (0.2). Maximization (`sense = MaxSense`) is supported.
#'
#' @return An algorithm object.
#' @examples
#' f <- function(x, p) (p[1] - x[1])^2
#' prob <- OptimizationProblem(f, 0, 1, lb = -10, ub = 10)
#' set.seed(1)
#' round(solve(prob, GCMAESOpt())$u, 4)
#' @export
GCMAESOpt <- function()
  structure(list(name = "GCMAESOpt", traits = list(bounds = TRUE, rbounds = TRUE, callback = FALSE)),
            class = c("GCMAESOpt", "OptimizationAlgorithm"))

.gcmaes_ptp <- function(x) max(x) - min(x)

# The strategy parameters of GCMAES.CMAESOpt
.gcmaes_parameters <- function(N, lambda = 0L) {
  lam <- if (lambda == 0L) as.integer(round(4 + 3 * log(N))) else max(4L, as.integer(lambda))
  mu <- as.integer(ceiling(lam / 2))
  w <- log(mu + 1 / 2) - log(seq_len(mu))
  w <- w * (1 / sum(abs(w)))  # normalize!(w, 1) multiplies by the inverse norm
  mueff <- 1 / sum(w^2)
  cs <- (mueff + 2) / (N + mueff + 3)
  c1 <- 2 / ((N + 1.3)^2 + mueff)
  list(lam = lam, mu = mu, w = w, mueff = mueff, cc = 4 / (N + 4), cs = cs, c1 = c1,
       cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((N + 2)^2 + mueff)),
       ds = 1 + 2 * max(0, sqrt((mueff - 1) / (N + 1)) - 1) + cs, chin = sqrt(N) * (1 - 1 / (4 * N) + 1 / (21 * N^2)))
}

# GCMAES.minimize((f, g), x0, sigma0, lo, hi; maxiter, maxtime)
.gcmaes_minimize <- function(f, g, x0, sigma0, maxiter = 0L, maxtime = Inf, maxfevals = 0, lambda = 0L,
                             equal_best = Inf, lazydecomp = FALSE) {
  N <- length(x0)
  xbar <- as.numeric(x0)
  xmin <- xbar; fmin <- f(xbar); sigma <- sigma0
  pr <- .gcmaes_parameters(N, lambda)
  lam <- pr$lam; mu <- pr$mu; w <- pr$w; mueff <- pr$mueff
  cc <- pr$cc; cs <- pr$cs; c1 <- pr$c1; cmu <- pr$cmu; ds <- pr$ds; chin <- pr$chin
  pc <- numeric(N); ps <- numeric(N)
  D <- rep(sigma, N); B <- diag(N); BD <- B * rep(D, each = N); C <- diag(D^2, N)
  fmins <- numeric(); fmeds <- numeric(); feqls <- numeric()
  maxfevals <- if (maxfevals == 0) 1e3 * N^2 else maxfevals
  if (maxiter != 0) maxfevals <- maxiter * lam
  iter <- 0L; fcount <- 0; status <- 0L
  t0 <- proc.time()[["elapsed"]]; t_run <- 0
  histiter <- 10L + as.integer(round(30 * N / lam))
  while (fcount < maxfevals && t_run <= maxtime) {
    iter <- iter + 1L; fcount <- fcount + lam
    # update_candidates!
    arz <- matrix(stats::rnorm(N * lam), N, lam)
    ary <- BD %*% arz
    arx <- xbar + sigma * ary
    arf <- vapply(seq_len(lam), function(k) f(arx[, k]), 0)
    idx <- order(arf)
    arf <- arf[idx]
    if (arf[1L] < fmin) { xmin <- arx[, idx[1L]]; fmin <- arf[1L] }
    fmins <- c(fmins, arf[1L])
    fmeds <- c(fmeds, arf[ceiling(lam / 2)])
    feqls <- c(feqls, if (lam > 1L) as.numeric(arf[1L] == arf[ceiling(1.1 + lam / 4)]) else 0)
    # update_parameters!
    imu <- idx[seq_len(mu)]
    xold <- xbar
    xbar <- as.numeric(arx[, imu, drop = FALSE] %*% w)
    if (is.function(g)) {
      # the line search along the negative gradient: in one process the only
      # step length tried is zero, so the mean itself is evaluated
      d <- -g(xbar)
      if (sqrt(sum(d^2)) != 0) {
        fx <- f(xbar)
        if (fx < fmin) { xmin <- xbar; fmin <- fx }
      }
    }
    zbar <- as.numeric(arz[, imu, drop = FALSE] %*% w)
    ps <- (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * as.numeric(B %*% zbar)
    hsig <- sqrt(sum(ps^2)) / sqrt(1 - (1 - cs)^(2 * iter)) / chin < 1.4 + 2 / (N + 1)
    pc <- pc * (1 - cc) + (hsig * sqrt(cc * (2 - cc) * mueff) / sigma) * (xbar - xold)
    C <- C * (1 - c1 - cmu + (1 - hsig) * c1 * cc * (2 - cc)) + c1 * tcrossprod(pc)
    art <- ary[, imu, drop = FALSE]
    C <- C + cmu * ((art * rep(w, each = N)) %*% t(art))
    sigma <- sigma * exp((sqrt(sum(ps^2)) / chin - 1) * cs / ds)
    if (!lazydecomp || (iter %% (1 / (c1 + cmu) / N / 10)) < 1) {
      e <- eigen((C + t(C)) / 2, symmetric = TRUE)
      D <- sqrt(pmax(e$values, 0)); B <- e$vectors
      BD <- B * rep(D, each = N)
    }
    t_run <- proc.time()[["elapsed"]] - t0
    # terminate
    nf <- length(fmins)
    stagiter <- min(as.integer(round(0.2 * nf + 120 + 30 * N / lam)), 20000L)
    hist <- fmins[max(1L, nf - histiter):nf]
    if (arf[1L] == arf[ceiling(0.7 * lam)]) sigma <- sigma * exp(0.2 + cs / ds)
    stop_now <- nf >= 100 + 50 * (N + 3)^2 / sqrt(lam) ||
      (length(hist) > histiter && .gcmaes_ptp(hist) < 1e-13) ||
      (nf > N && mean(feqls[(nf - N):nf]) > 0.33) ||
      (sqrt(max(abs(diag(C)))) * sigma / sigma0 < 1e-11 && max(abs(pc)) * sigma / sigma0 < 1e-11) ||
      sigma / sigma0 > 1e20 * max(D) ||
      (nf > stagiter && stats::median(fmins[(nf - 20L):nf]) >= stats::median(fmins[(nf - stagiter):(nf - stagiter + 20L)]) &&
         stats::median(fmeds[(nf - 20L):nf]) >= stats::median(fmeds[(nf - stagiter):(nf - stagiter + 20L)])) ||
      max(D) > 1e7 * min(D) ||
      (length(hist) > histiter && .gcmaes_ptp(c(hist, arf)) < 1e-11) ||
      (nf > equal_best && .gcmaes_ptp(utils::tail(cummin(fmins), equal_best + 1)) < 1e-12)
    if (isTRUE(stop_now)) { status <- 1L; break }
  }
  list(xmin = xmin, fmin = fmin, status = status, iterations = iter, fevals = fcount, sigma = sigma, mean = xbar,
       C = C)
}

.mtk$hooks$opt_solve_GCMAESOpt <- function(a, alg, opts) {
  .opt_check_traits(a, alg)
  p <- a$p
  loss <- function(x) as.numeric(a$f$f(x, p))[1L]
  g <- if (is.function(a$f$grad)) function(x) as.numeric(a$f$grad(x, p)) else NULL
  r <- .gcmaes_minimize(loss, g, as.numeric(a$u0), opts$sigma0 %||% 0.2, maxiter = opts$maxiters %||% 0L,
                        maxtime = opts$maxtime %||% Inf)
  .opt_solution(a, alg, r$xmin, r$fmin, if (r$status == 1L) "Success" else "Failure", r)
}

# ---------------------------------------------------------------------------
# pycma (the Python package cma, version 4)

#' CMA-ES of the Python package cma
#'
#' `PyCMAOpt()` minimizes with the default CMA-ES of the Python package
#' `cma` (pycma), ported to R: active covariance matrix adaptation with
#' negative recombination weights, cumulative step size adaptation, lazy
#' eigendecomposition, the bound transformation `BoundTransform` (a
#' periodic linear-quadratic map into `lb` and `ub`, with the coordinate
#' standard deviations limited to a third of the bound range) and the
#' termination criteria of pycma. No derivatives are used. As in
#' OptimizationPyCMA, the initial step size is 1.
#'
#' The `solve` arguments map as in OptimizationPyCMA: `maxiters` to the
#' option `maxiter`, `maxtime` to `timeout` (seconds), `abstol` to
#' `tolfun` and `reltol` to `tolfunrel`. Further pycma options can be
#' given by name: `popsize`, `maxfevals`, `ftarget`, `tolx`,
#' `tolfunhist`, `tolstagnation`, `tolupsigma`, `tolconditioncov`,
#' `tolfacupx`, `tolflatfitness`, `tolxstagnation`, `maxstd`, `CMA_active`
#' and `seed` (the random numbers are R's). The return code is `Success`
#' for the criteria `ftarget`, `tolfun` and `tolx`, `MaxIters` for
#' `maxiter` and `maxfevals`, `MaxTime` for `timeout`, `Failure` for the
#' stagnation and numerical criteria and `Terminated` when the callback
#' halts. The solution's `original` holds `xbest`, `fbest`, `countiter`,
#' `countevals`, `stop` (the criteria met), `mean`, `sigma` and `C`.
#'
#' @return An algorithm object.
#' @examples
#' rosenbrock <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
#' prob <- OptimizationProblem(rosenbrock, c(0, 0), c(1, 100), lb = c(-1, -1), ub = c(0.8, 0.8))
#' sol <- solve(prob, PyCMAOpt(), seed = 42)
#' round(sol$u, 3); sol$retcode
#' @export
PyCMAOpt <- function()
  structure(list(name = "PyCMAOpt", traits = list(bounds = TRUE, callback = TRUE)),
            class = c("PyCMAOpt", "OptimizationAlgorithm"))

# Python's sum() of floats (CPython 3.12 and later): left to right with the
# compensated summation of Neumaier
.py_sum <- function(x) {
  s <- 0; c <- 0
  for (v in x) {
    t <- s + v
    c <- c + if (abs(s) >= abs(v)) (s - t) + v else (v - t) + s
    s <- t
  }
  if (c != 0 && is.finite(c)) s + c else s
}

# Mh.chiN: the expected norm of a standard normal vector, tabulated up to
# dimension 100 (the values of pycma, as exact hexadecimal constants)
.pycma_chin_table <- c(
  0x1.9884533d43652p-1, 0x1.40d931ff62706p+0, 0x1.9884533d43652p+0, 0x1.e145caff13a89p+0, 0x1.1058377e2cee1p+1,
  0x1.2ccb9edf6c496p+1, 0x1.46d04297691dap+1, 0x1.5eed8eaf53aafp+1, 0x1.75804c1ac1467p+1, 0x1.8acb40853e204p+1,
  0x1.9f00548f816acp+1, 0x1.b245fa2c2abd1p+1, 0x1.c4ba8acb18d19p+1, 0x1.d67679afd8f79p+1, 0x1.e78dd08bf357cp+1,
  0x1.f811393c6876fp+1, 0x1.04075e2881c87p+2, 0x1.0bc92668177f3p+2, 0x1.1353186720016p+2, 0x1.1aa9a88a51b0ep+2,
  0x1.21d0bb5f14376p+2, 0x1.28cbbdc46f603p+2, 0x1.2f9db81a76b3dp+2, 0x1.36495daa74703p+2, 0x1.3cd11926bea57p+2,
  0x1.433716e6e3f4ep+2, 0x1.497d4d5b7e978p+2, 0x1.4fa5841727cd1p+2, 0x1.55b159b4376ddp+2, 0x1.5ba248ced6efcp+2,
  0x1.6179ac3ed83cap+2, 0x1.6738c2b399d5bp+2, 0x1.6ce0b1cd424f1p+2, 0x1.727288c936a46p+2, 0x1.77ef42d37a9f0p+2,
  0x1.7d57c90b5de55p+2, 0x1.82acf4473c4bfp+2, 0x1.87ef8ea1040f5p+2, 0x1.8d2054d38a08bp+2, 0x1.923ff76f5bbedp+2,
  0x1.974f1beca7d48p+2, 0x1.9c4e5d9ef13d5p+2, 0x1.a13e4e8e86759p+2, 0x1.a61f783b1b8e1p+2, 0x1.aaf25c4a65dd9p+2,
  0x1.afb775252da2ap+2, 0x1.b46f3684f0a9cp+2, 0x1.b91a0df3e647bp+2, 0x1.bdb86340f5c86p+2, 0x1.c24a98e8fb13fp+2,
  0x1.c6d10c7687dbfp+2, 0x1.cb4c16d9290f3p+2, 0x1.cfbc0cb51712cp+2, 0x1.d4213eac1b144p+2, 0x1.d87bf9a05b218p+2,
  0x1.dccc86f1a9cd6p+2, 0x1.e1132cb5e3c52p+2, 0x1.e5502decd5fa8p+2, 0x1.e983cab0192b4p+2, 0x1.edae405f43990p+2,
  0x1.f1cfc9c8c7280p+2, 0x1.f5e89f4fc8fdbp+2, 0x1.f9f8f70f37891p+2, 0x1.fe0104fa5cbfdp+2, 0x1.01007d7d920cbp+3,
  0x1.02fc848723197p+3, 0x1.04f4aeaec38edp+3, 0x1.06e9122c17feap+3, 0x1.08d9c464f453ep+3, 0x1.0ac6d9f80949ep+3,
  0x1.0cb066c6e19b9p+3, 0x1.0e967dff3c9f0p+3, 0x1.10793223d2c1dp+3, 0x1.125895148f3dcp+3, 0x1.1434b816495ecp+3,
  0x1.160dabda06cc5p+3, 0x1.17e38083cf780p+3, 0x1.19b645b11b193p+3, 0x1.1b860a7edf762p+3, 0x1.1d52dd8f461cdp+3,
  0x1.1f1ccd0f0fa86p+3, 0x1.20e3e6baaa304p+3, 0x1.22a837e3000c7p+3, 0x1.2469cd7203adcp+3, 0x1.2628b3eefcf6fp+3,
  0x1.27e4f7829c1a8p+3, 0x1.299ea3fad5c68p+3, 0x1.2b55c4ce8c0eep+3, 0x1.2d0a652107504p+3, 0x1.2ebc8fc542038p+3,
  0x1.306c4f410a459p+3, 0x1.3219adcffba2dp+3, 0x1.33c4b56653872p+3, 0x1.356d6fb3a28b8p+3, 0x1.3713e6255caf1p+3,
  0x1.38b821e94a66ep+3, 0x1.3a5a2befdc4ffp+3, 0x1.3bfa0cee632d5p+3, 0x1.3d97cd612dc22p+3, 0x1.3f33758d8dfc1p+3)
.pycma_chin <- function(N)
  if (N <= 100) .pycma_chin_table[N] else sqrt(N) * (1 - 1 / (4 * N) + 1 / (32 * N^2))

# BoxConstraintsLinQuadTransformation with margin_width2
.pycma_box <- function(lb, ub, N) {
  if (is.null(lb) && is.null(ub)) return(NULL)
  lb <- rep_len(as.numeric(lb %||% -Inf), N); ub <- rep_len(as.numeric(ub %||% Inf), N)
  if (all(!is.finite(lb)) && all(!is.finite(ub))) return(NULL)
  if (any(lb >= ub)) stop("Lower bounds need to be smaller than upper bounds.", call. = FALSE)
  al <- ifelse(is.finite(lb), pmin((ub - lb) / 2, pmax(1, abs(lb)) / 20), 1)
  au <- ifelse(is.finite(ub), pmin((ub - lb) / 2, pmax(1, abs(ub)) / 20), 1)
  list(lb = lb, ub = ub, al = al, au = au)
}

.pycma_transform <- function(b, y) {
  if (is.null(b)) return(y)
  lb <- b$lb; ub <- b$ub; al <- b$al; au <- b$au
  i <- y < lb - 2 * al - (ub - lb) / 2 | y > ub + 2 * au + (ub - lb) / 2
  if (any(i)) {
    r <- 2 * (ub[i] - lb[i] + al[i] + au[i])
    s <- lb[i] - 2 * al[i] - (ub[i] - lb[i]) / 2
    y[i] <- y[i] - r * floor((y[i] - s) / r)
  }
  i <- y > ub + au; y[i] <- y[i] - 2 * (y[i] - ub[i] - au[i])
  i <- y < lb - al; y[i] <- y[i] + 2 * (lb[i] - al[i] - y[i])
  i <- y < lb + al; y[i] <- lb[i] + (y[i] - (lb[i] - al[i]))^2 / 4 / al[i]
  i <- y > ub - au; y[i] <- ub[i] - (y[i] - (ub[i] + au[i]))^2 / 4 / au[i]
  y
}

# The recombination weights and the strategy parameters of pycma
# (CMAParameters, RecombinationWeights, CMAAdaptSigmaCSA) and the default
# maxiter and tolstagnation
.pycma_parameters <- function(N, popsize = NULL, CMA_active = TRUE) {
  lam <- as.integer(floor(popsize %||% (4 + 3 * log(N))))
  lam_mirr <- if (lam < 6) as.integer(0.5 + 0.16 * min(lam, 2 * N + 2) + 0.29) else 0L
  # RecombinationWeights
  w <- log((lam + 1) / 2) - log(seq_len(lam))
  mu <- sum(w > 0)
  w <- w / .py_sum(w[seq_len(mu)])
  mueff <- 1 / .py_sum(w[seq_len(mu)]^2)
  neg <- if (mu < lam) (mu + 1L):lam else integer()
  sum_neg <- .py_sum(w[neg])
  if (length(neg) && sum_neg != 0) w[neg] <- w[neg] / -sum_neg
  # CMAParameters
  cc <- (4 + mueff / N) / (N + (4 + 2 * mueff / N))
  c1 <- min(1, lam / 6) * 2 / ((N + 1.3)^2 + mueff)
  cmu <- min(1 - c1, 2 * (0.25 + mueff + 1 / mueff - 2) / ((N + 2)^2 + 2 * mueff / 2))
  # the learning rates of GaussFullSampler.parameters, used for the update
  # and for the final negative weights
  cmu_sm <- min(1 - c1, 2 * (0.25 + mueff - 2 + 1 / mueff) / ((N + 2)^2 + 2 * mueff / 2))
  finalize <- function(c1, cmu) {
    if (!(length(neg) && w[lam] < 0)) return()
    if (!isTRUE(CMA_active)) { w[neg] <<- 0; return() }
    lim <- function(v) {
      s <- .py_sum(w[neg])
      if (s < -abs(v)) { fac <- abs(v / s); if (fac < 1) w[neg] <<- w[neg] * fac }
    }
    if (cmu > 0) {
      w[neg] <<- w[neg] * abs((1 + c1 / cmu) / .py_sum(w[neg]))
      lim((1 - c1 - cmu) / cmu / N)
    }
    lim(1 + 2 * (.py_sum(w[neg])^2 / .py_sum(w[neg]^2)) / (mueff + 2))
  }
  finalize(c1, cmu)
  finalize(c1, cmu_sm)
  # CMAAdaptSigmaCSA
  cs <- min(0.5, (mueff + 2) / (N + (mueff + 3)))
  die <- max(1, 3 * (1 - 0.5^(N / 9)))
  damps <- 0.5 + min(1, (lam_mirr / (0.159 * lam) - 1)^2) / 2 + 2 * max(0, die * ((mueff - 1) / (N + 1))^0.5 - 1) + cs
  chin <- .pycma_chin(N)
  # the defaults of maxiter and tolstagnation use the option value of
  # popsize, which is not rounded
  ps <- popsize %||% (4 + 3 * log(N))
  list(lam = lam, mu = mu, w = w, mueff = mueff, cc = cc, c1 = c1, cmu = cmu, cmu_sm = cmu_sm, cs = cs, damps = damps,
       chin = chin, lam_mirr = lam_mirr, maxiter = 100 + floor(150 * (N + 3)^2 / sqrt(ps)),
       tolstagnation = as.integer(100 + 100 * N^1.5 / ps))
}

# CMAEvolutionStrategy(x0, sigma0, opts).optimize(f) of pycma with the
# default options
.pycma_minimize <- function(f, x0, sigma0 = 1, lb = NULL, ub = NULL, popsize = NULL, maxiter = NULL,
                            maxfevals = Inf, ftarget = -Inf, tolfun = 1e-11, tolfunrel = 0, tolfunhist = 1e-12,
                            tolx = 1e-11, tolstagnation = NULL, tolupsigma = 1e20, tolconditioncov = 1e14,
                            tolfacupx = 1e3, tolflatfitness = 1, tolxstagnation = c(1e-9, 20, 0.1), timeout = Inf,
                            maxstd = NULL, CMA_active = TRUE) {
  N <- length(x0)
  pr <- .pycma_parameters(N, popsize, CMA_active)
  lam <- pr$lam; mu <- pr$mu; w <- pr$w; mueff <- pr$mueff; cc <- pr$cc; c1 <- pr$c1; cmu <- pr$cmu_sm
  cs <- pr$cs; damps <- pr$damps; chin <- pr$chin
  maxiter <- maxiter %||% pr$maxiter
  tolstagnation <- tolstagnation %||% pr$tolstagnation
  # GaussFullSampler
  C <- diag(exp(1e-4 / N * (0:(N - 1L))), N)
  B <- diag(N); D <- sqrt(diag(C))
  gap <- 1 / (c1 + pr$cmu + 0x1.82db34012b251p-77) / N / 10
  count_tell <- 0L; last_update <- 0L
  decompose <- function() {
    repeat {
      C <<- (C + t(C)) / 2
      e <- eigen(C, symmetric = TRUE)
      if (all(e$values > 0)) break
      C <<- C + diag(min(D)^2, N)
    }
    o <- order(e$values)
    D <<- sqrt(e$values[o]); B <<- e$vectors[, o, drop = FALSE]
  }
  box <- .pycma_box(lb, ub, N)
  mean <- .cmaes_backtransform(box, as.numeric(x0))
  sigma <- sigma0
  scaling <- 1
  if (is.null(maxstd) && !is.null(box)) maxstd <- (box$ub - box$lb) / 3
  maxstd <- rep_len(maxstd %||% Inf, N)
  countiter <- 0L; countevals <- 0
  into_limits <- function() {
    if (!any(is.finite(maxstd))) return()
    fix <- if (countiter <= 2L) exp(1e-4) else 1
    s <- sigma * scaling * sqrt(diag(C))
    for (i in which(s > maxstd * fix)) {
      if (length(scaling) == 1L) scaling <<- rep(scaling, N)
      if (s[i] != 0 && is.finite(scaling[i] * maxstd[i])) scaling[i] <<- scaling[i] * maxstd[i] / s[i]
    }
  }
  into_limits()
  ps <- numeric(N); pc <- numeric(N)
  best <- list(x = NULL, f = Inf, evals = 0)
  fit <- NULL; hist <- numeric(); histbest <- numeric(); histmedian <- numeric()
  median0 <- NULL; median_min <- Inf; flat <- 0L
  xs <- list(x = mean, count = 0L, count_x = 0L)
  t0 <- proc.time()[["elapsed"]]
  stopd <- character()
  repeat {
    if (countiter > 0L) {
      sx <- sigma * scaling * sqrt(diag(C))
      cur <- max(fit) - min(fit); his <- max(hist) - min(hist)
      stag <- FALSE
      if (tolstagnation) {
        l <- max(tolstagnation / 5 / 2, length(histbest) / 10)
        if (l <= countiter) {
          l <- as.integer(l)
          stag <- countiter > N * (5 + 100 / lam) && (countevals - best$evals) / lam > tolstagnation / 2 &&
            length(histbest) > 100 && 2 * l < length(histbest) &&
            stats::median(histmedian[seq_len(l)]) >= stats::median(histmedian[(l + 1L):(2L * l)]) &&
            stats::median(histbest[seq_len(l)]) >= stats::median(histbest[(l + 1L):(2L * l)])
        }
      }
      ia <- countiter %% N + 1L
      if (fit[1L] < fit[as.integer(0.75 * lam) + 1L]) flat <- 0L else flat <- flat + 1L
      conds <- c(ftarget = best$f <= ftarget, maxfevals = countevals - 1 >= maxfevals, maxiter = countiter >= maxiter,
                 tolfacupx = any(sx > sigma0 * tolfacupx),
                 tolx = all(sx < tolx) && all(sigma * scaling * pc < tolx),
                 tolfun = cur < tolfun && his < tolfun,
                 tolfunrel = cur < tolfunrel * (median0 - median_min),
                 tolfunhist = length(hist) > 9 && his < tolfunhist, tolstagnation = stag,
                 tolxstagnation = xs$count - xs$count_x > tolxstagnation[2L] + tolxstagnation[3L] * xs$count,
                 tolupsigma = tolupsigma > 0 && sigma / max(D) > sigma0 * tolupsigma,
                 timeout = proc.time()[["elapsed"]] - t0 > timeout,
                 noeffectcoord = any(mean == mean + 0.2 * sx),
                 noeffectaxis = all(mean == mean + 0.1 * sigma * D[ia] * scaling * B[, ia]),
                 tolconditioncov = tolconditioncov > 0 && D[N] > sqrt(tolconditioncov) * D[1L],
                 tolflatfitness = flat > tolflatfitness)
      stopd <- names(conds)[!is.na(conds) & conds]
      if (length(stopd)) break
    }
    # ask
    if (!(count_tell < last_update + gap || (gap == count_tell - last_update && count_tell == last_update))) {
      decompose(); last_update <- count_tell
    }
    ary <- scaling * (B %*% (D * matrix(stats::rnorm(N * lam), N, lam)))
    X <- mean + sigma * ary
    P <- matrix(apply(X, 2L, function(x) .pycma_transform(box, x)), N)
    fv <- vapply(seq_len(lam), function(k) as.numeric(f(P[, k]))[1L], 0)
    # tell
    countiter <- countiter + 1L; countevals <- countevals + lam
    if (anyNA(fv)) fv[is.na(fv)] <- stats::median(fv[!is.na(fv)])
    k <- which.min(fv)
    if (length(k) && fv[k] < best$f) best <- list(x = P[, k], f = fv[k], evals = countevals - lam + k)
    idx <- order(fv)
    fit <- fv[idx]
    hist <- c(fit[1L], hist)
    med <- stats::median(fit)
    if (countiter %% 5L == 0L) { histbest <- c(fit[1L], histbest); histmedian <- c(med, histmedian) }
    if (length(histbest) > 2e4) { histbest <- histbest[-length(histbest)]; histmedian <- histmedian[-length(histmedian)] }
    if (length(hist) > 10 + 30 * N / lam) hist <- hist[-length(hist)]
    if (is.null(median0)) median0 <- med
    median_min <- min(median_min, med)
    pop <- X[, idx, drop = FALSE]
    mold <- mean
    mean <- as.numeric(pop[, seq_len(mu), drop = FALSE] %*% w[seq_len(mu)])
    z <- as.numeric(B %*% (crossprod(B, (mean - mold) / scaling) / D)) * sqrt(mueff) / sigma
    ps <- (1 - cs) * ps + sqrt(cs * (2 - cs)) * z
    hsig <- sum(ps^2) / (1 - (1 - cs)^(2 * countiter)) / N - 1 < 1 + 4 / (N + 1)
    c1a <- c1 * (1 - (1 - hsig^2) * cc * (2 - cc))
    pc <- (1 - cc) * pc + hsig * (sqrt(cc * (2 - cc) * mueff) / sigma) * (mean - mold) / scaling
    if (c1a + cmu > 0) {
      V <- cbind(sqrt(c1 / (c1a + 0x1.82db34012b251p-77)) * pc, (pop - mold) / (sigma * scaling))
      wts <- c(c1a, cmu * w)
      C <- C * (1 - sum(wts))
      for (j in which(wts < 0)) wts[j] <- wts[j] * N / (sqrt(sum((crossprod(B, V[, j]) / D)^2)) + 1e-9)^2
      C <- C + (V * rep(wts, each = N)) %*% t(V)
      count_tell <- count_tell + 1L
    }
    sigma <- sigma * exp(min(max((sqrt(sum(ps^2)) / chin - 1) * cs / damps, -1), 1))
    into_limits()
    if (sigma > 1e9 * sigma0) {
      alpha <- sigma / sqrt(max(diag(C)))
      if (alpha > 1) {
        C <- C * alpha; D <- D * sqrt(alpha)
        sigma <- sigma / sqrt(alpha); tolupsigma <- tolupsigma / sqrt(alpha)
      }
    }
    # _StopTolXStagnation
    xs$count <- xs$count + 1L
    dl <- tolxstagnation[1L] * sqrt(max(1, (xs$count - xs$count_x) / (tolxstagnation[2L] + tolxstagnation[3L] * xs$count)))
    if (dl < 0 || sqrt(sum((xs$x - mean)^2)) > dl) { xs$x <- mean; xs$count_x <- xs$count }
  }
  list(xbest = best$x, fbest = best$f, countiter = countiter, countevals = countevals, stop = stopd, mean = mean,
       sigma = sigma, C = C, popsize = lam)
}

.pycma_retcode <- function(stopd) {
  if (any(c("ftarget", "tolfun", "tolx") %in% stopd)) "Success" else
    if (any(c("maxiter", "maxfevals") %in% stopd)) "MaxIters" else
    if ("timeout" %in% stopd) "MaxTime" else
    if (any(c("tolupsigma", "tolconditioncov", "noeffectcoord", "noeffectaxis", "tolxstagnation", "tolflatfitness",
              "tolfacupx", "tolstagnation") %in% stopd)) "Failure" else "Default"
}

.mtk$hooks$opt_solve_PyCMAOpt <- function(a, alg, opts) {
  .opt_check_traits(a, alg)
  p <- a$p
  loss <- function(x) as.numeric(a$f$f(x, p))[1L]
  keys <- c("popsize", "maxfevals", "ftarget", "tolx", "tolfunhist", "tolstagnation", "tolupsigma", "tolconditioncov",
            "tolfacupx", "tolflatfitness", "tolxstagnation", "maxstd", "CMA_active", "tolfun", "tolfunrel")
  args <- opts[intersect(names(opts), keys)]
  args$tolfun <- opts$abstol %||% args$tolfun
  args$tolfunrel <- opts$reltol %||% args$tolfunrel
  args$timeout <- opts$maxtime
  args$maxiter <- opts$maxiters
  args <- args[!vapply(args, is.null, TRUE)]
  r <- .with_seed(opts$seed, do.call(.pycma_minimize, c(list(loss, as.numeric(a$u0), 1, a$lb, a$ub), args)))
  .opt_solution(a, alg, r$xbest, r$fbest, .pycma_retcode(r$stop), r)
}
