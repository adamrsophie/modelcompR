# The derivative free and subgradient solvers of Manopt.jl (Nelder Mead,
# particle swarm, CMA-ES, convex bundle method), the OptimizationManopt
# optimizers and the `Manopt` module binding.

.mo_nelder_mead <- function(mp, population = NULL, stopping_criterion = NULL, alpha = 1, gamma = 2, rho = 0.5,
                            sigma = 0.5, retraction_method = "exp", ...) {
  M <- mp$M
  pts <- population %||% lapply(seq_len(.mf_dim(M) + 1L), function(i) .mf_rand(M))
  sc <- stopping_criterion %||% list(.mo_sc("StopAfterIteration", n = 2000L),
                                     .mo_sc("StopWhenPopulationConcentrated", tol_f = 1e-8, tol_p = 1e-8))
  st <- .mo_state(mp, pts[[1L]], sc)
  st$pts <- pts
  st$costs <- vapply(pts, mp$f, 0)
  st$p <- pts[[which.min(st$costs)]]
  .mo_run(st, function(k) {
    o <- order(st$costs)
    st$costs <- st$costs[o]
    st$pts <- st$pts[o]
    np <- length(st$pts)
    m <- .mf_mean(M, st$pts[-np])
    xi <- .mf_invretract(M, m, st$pts[[np]])
    xr <- .mf_retract(M, m, -alpha * xi, 1, retraction_method)
    cr <- mp$f(xr)
    go <- TRUE
    if (cr >= st$costs[1L] && cr <= st$costs[np - 1L]) {
      st$pts[[np]] <- xr; st$costs[np] <- cr; go <- FALSE
    }
    if (cr < st$costs[1L]) {
      xe <- .mf_retract(M, m, -gamma * alpha * xi, 1, retraction_method)
      ce <- mp$f(xe)
      st$pts[[np]] <- if (ce < cr) xe else xr
      st$costs[np] <- min(ce, cr)
      go <- FALSE
    }
    if (go && cr > st$costs[np - 1L]) {
      xc <- .mf_retract(M, m, (if (cr < st$costs[np]) -rho else rho) * xi, 1, retraction_method)
      cc <- mp$f(xc)
      if (cc < st$costs[np]) {
        st$pts[[np]] <- xc; st$costs[np] <- cc; go <- FALSE
      }
    }
    if (go) {
      for (i in 2:np) {
        st$pts[[i]] <- .mf_retract(M, st$pts[[1L]], .mf_invretract(M, st$pts[[1L]], st$pts[[i]]), sigma,
                                   retraction_method)
        st$costs[i] <- mp$f(st$pts[[i]])
      }
    }
    st$p <- st$pts[[which.min(st$costs)]]
  })
}

.mo_nm_concentrated <- function(M, st, c) {
  np <- length(st$pts)
  vf <- max(abs(st$costs[1L] - st$costs[-1L]))
  vp <- max(vapply(st$pts[-1L], function(q) .mf_distance(M, st$pts[[1L]], q), 0))
  vf < c$tol_f && vp < c$tol_p
}

.mo_particle_swarm <- function(mp, swarm, stopping_criterion = NULL, inertia = 0.65, social_weight = 1.4,
                               cognitive_weight = 1.4, retraction_method = "exp", vector_transport_method = "parallel",
                               ...) {
  M <- mp$M
  sc <- stopping_criterion %||% list(.mo_sc("StopAfterIteration", n = 500L),
                                     .mo_sc("StopWhenSwarmVelocityLess", tol = 1e-4))
  st <- .mo_state(mp, swarm[[1L]], sc)
  st$swarm <- swarm
  st$velocity <- lapply(swarm, function(y) .mf_rand_vector(M, y))
  st$best <- swarm
  st$p <- swarm[[which.min(vapply(swarm, mp$f, 0))]]
  .mo_run(st, function(k) {
    for (i in seq_along(st$swarm)) {
      x <- st$swarm[[i]]
      cog <- .mf_invretract(M, x, st$best[[i]])
      soc <- .mf_invretract(M, x, st$p)
      v <- inertia * st$velocity[[i]] + cognitive_weight * stats::runif(1L) * cog +
        social_weight * stats::runif(1L) * soc
      xn <- .mf_retract(M, x, v, 1, retraction_method)
      st$velocity[[i]] <- .mf_transport(M, x, v, xn, vector_transport_method)
      st$swarm[[i]] <- xn
      c <- mp$f(xn)
      if (c < mp$f(st$best[[i]])) {
        st$best[[i]] <- xn
        if (c < mp$f(st$p)) st$p <- xn
      }
    }
  })
}

# CMA-ES on manifolds (Colutto, Fruhauf, Fuchs and Scherzer 2010) as in
# Manopt.cma_es: the search distribution lives in the tangent space at the
# mean, expressed in an orthonormal basis.
.mo_cma_es <- function(mp, p, sigma = 1, lambda = NULL, tol_fun = 1e-12, tol_x = 1e-12, stopping_criterion = NULL,
                       retraction_method = "exp", vector_transport_method = "parallel", ...) {
  M <- mp$M
  n <- .mf_dim(M)
  lam <- as.integer(lambda %||% (4L + as.integer(floor(3 * log(n)))))
  wp <- log((lam + 1) / 2) - log(seq_len(lam))
  mu <- lam %/% 2L
  mueff <- sum(wp[seq_len(mu)])^2 / sum(wp[seq_len(mu)]^2)
  mueffm <- sum(wp[(mu + 1L):lam])^2 / sum(wp[(mu + 1L):lam]^2)
  a_mueffm <- 1 + 2 * mueffm / (mueff + 2)
  c1 <- 2 / ((n + 1.3)^2 + mueff)
  cmu <- min(1 - c1, 2 * (0.25 + mueff + 1 / mueff - 2) / ((n + 2)^2 + 2 * mueff / 2))
  a_mu <- 1 + c1 / cmu
  a_posdef <- (1 - c1 - cmu) / (n * cmu)
  wpos <- 1 / sum(wp[wp > 0])
  wneg <- -min(a_mu, a_mueffm, a_posdef) / sum(wp[wp < 0])
  w <- ifelse(wp >= 0, wp * wpos, wp * wneg)
  cs <- (mueff + 2) / (n + mueff + 5)
  ds <- 1 + 2 * max(0, sqrt((mueff - 1) / (n + 1)) - 1) + cs
  cc <- (4 + mueff / n) / (n + 4 + 2 * mueff / n)
  sc <- stopping_criterion %||% list(
    .mo_sc("StopAfterIteration", n = 50000L), .mo_sc("StopWhenCovarianceIllConditioned", tol = 1e14),
    .mo_sc("StopWhenBestCostInGenerationConstant", n = as.integer(10 + ceiling(30 * n / lam))),
    .mo_sc("StopWhenEvolutionStagnates", min_size = as.integer(120 + ceiling(30 * n / lam)), max_size = 20000L,
           fraction = 0.3),
    .mo_sc("StopWhenPopulationDiverges", tol = 1e4),
    .mo_sc("StopWhenPopulationCostConcentrated", tol = tol_fun, n = as.integer(10 + ceiling(30 * n / lam))),
    .mo_sc("StopWhenPopulationStronglyConcentrated", tol = tol_x))
  st <- .mo_state(mp, p, sc)
  st$p_m <- p
  st$p_obj <- Inf
  st$sigma <- sigma
  C <- diag(n)
  ev <- eigen(C, symmetric = TRUE)
  st$dev <- sqrt(ev$values)
  st$cond <- max(ev$values) / min(ev$values)
  st$ps <- numeric(n); st$pc <- numeric(n)
  emv <- sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n^2))
  last_var <- ev$values
  .mo_run(st, function(k) {
    D2 <- ev$values; B <- ev$vectors
    st$cond <- max(abs(D2)) / min(abs(D2))
    if (min(D2) <= 0) {
      warning("Covariance matrix has nonpositive eigenvalues; try reformulating the objective, modifying ",
              "stopping criteria or adjusting optimization parameters to avoid this.", call. = FALSE)
      D2[D2 <= 0] <- last_var[D2 <= 0]
    }
    st$dev <- sqrt(D2)
    cinv <- B %*% (t(B) / st$dev)
    ys <- lapply(seq_len(lam), function(i) as.numeric(B %*% (stats::rnorm(n) * st$dev)))
    pop <- lapply(ys, function(y) .mf_retract(M, st$p_m, .mf_get_vector(M, st$p_m, y), st$sigma, retraction_method))
    fit <- vapply(pop, mp$f, 0)
    st$best_fit <- min(fit); st$worst_fit <- max(fit); st$median_fit <- stats::median(fit)
    for (i in seq_along(fit)) if (fit[i] < st$p_obj) { st$p_obj <- fit[i]; st$p <- pop[[i]] }
    ys <- ys[order(fit)]
    yw <- Reduce(`+`, lapply(seq_len(mu), function(i) w[i] * ys[[i]]))
    new_m <- .mf_retract(M, st$p_m, .mf_get_vector(M, st$p_m, yw), st$sigma, retraction_method)
    cy <- as.numeric(cinv %*% yw)
    st$ps <- (1 - cs) * st$ps + sqrt(cs * (2 - cs) * mueff) * cy
    st$sigma <- st$sigma * exp(cs / ds * (.jl_norm(st$ps) / emv - 1))
    st$pc <- st$pc * (1 - cc)
    if (.jl_norm(st$ps) / sqrt(1 - (1 - cs)^(2 * k)) < (1.4 + 2 / (n + 1)) * emv) {
      st$pc <- st$pc + sqrt(cc * (2 - cc) * mueff) * yw
      dh <- 0
    } else {
      dh <- cc * (2 - cc)
    }
    C <<- C * (1 + c1 * dh - c1 - cmu * sum(w)) + c1 * tcrossprod(st$pc)
    for (i in seq_len(lam)) {
      wi <- w[i]
      if (wi < 0) wi <- wi * n / .jl_norm(as.numeric(cinv %*% ys[[i]]))^2
      C <<- C + (cmu * wi) * tcrossprod(ys[[i]])
    }
    last_var <<- D2
    e2 <- eigen((C + t(C)) / 2, symmetric = TRUE)
    V <- e2$vectors
    if (.mf_sphere(M)) {
      V <- vapply(seq_len(n), function(j) .mf_get_coordinates(M, new_m, .mf_transport(M, st$p_m,
        .mf_get_vector(M, st$p_m, V[, j]), new_m, vector_transport_method)), numeric(n))
      V <- matrix(V, n, n)
      st$ps <- .mf_get_coordinates(M, new_m, .mf_transport(M, st$p_m, .mf_get_vector(M, st$p_m, st$ps), new_m,
                                                          vector_transport_method))
      st$pc <- .mf_get_coordinates(M, new_m, .mf_transport(M, st$p_m, .mf_get_vector(M, st$p_m, st$pc), new_m,
                                                          vector_transport_method))
    }
    ev <<- list(values = e2$values, vectors = V)
    C <<- V %*% (e2$values * t(V))
    st$p_m <- new_m
  })
}

.mo_cma_stop <- function(c, e, st, k) {
  switch(c$type,
    StopWhenCovarianceIllConditioned = k > 0 && st$cond > c$tol,
    StopWhenBestCostInGenerationConstant = {
      if (k <= 0) { e$best <- Inf; e$since <- 0L; return(FALSE) }
      if (e$since >= c$n) return(TRUE)
      if (e$best != st$best_fit) { e$best <- st$best_fit; e$since <- 0L } else e$since <- e$since + 1L
      FALSE
    },
    StopWhenEvolutionStagnates = {
      if (k <= 0) { e$bh <- numeric(); e$mh <- numeric(); return(FALSE) }
      N <- length(e$bh)
      if (N >= c$min_size) {
        lo <- ceiling(N * c$fraction); hi <- N - lo + 1
        if (lo >= 1 && hi >= 1 && stats::median(e$bh[seq_len(lo)]) <= stats::median(e$bh[hi:N]) &&
            stats::median(e$mh[seq_len(lo)]) <= stats::median(e$mh[hi:N])) return(TRUE)
      }
      e$bh <- utils::tail(c(e$bh, st$best_fit), c$max_size)
      e$mh <- utils::tail(c(e$mh, st$median_fit), c$max_size)
      FALSE
    },
    StopWhenPopulationDiverges = {
      if (k <= 0) { e$last <- st$sigma * max(st$dev); return(FALSE) }
      st$sigma * max(st$dev) / e$last > c$tol
    },
    StopWhenPopulationCostConcentrated = {
      if (k <= 0) { e$h <- numeric(); return(FALSE) }
      e$h <- utils::tail(c(e$h, st$best_fit), c$n)
      length(e$h) == c$n && max(e$h) - min(e$h) < c$tol && st$worst_fit - st$best_fit < c$tol
    },
    StopWhenPopulationStronglyConcentrated = k > 0 && max(abs(st$dev)) < c$tol && st$sigma * max(abs(st$pc)) < c$tol,
    stop("unknown stopping criterion: ", c$type, call. = FALSE))
}

# min 0.5 l'Hl + e'l over the unit simplex, by a primal active set method
.mo_qp_simplex <- function(H, e) {
  d <- length(e)
  if (d == 1L) return(1)
  H <- (H + t(H)) / 2
  ridge <- 1e-13 * (1 + max(abs(diag(H))))
  lam <- numeric(d)
  S <- which.min(diag(H) / 2 + e)
  lam[S] <- 1
  for (it in seq_len(20L * d + 50L)) {
    s <- length(S)
    K <- rbind(cbind(H[S, S, drop = FALSE] + diag(ridge, s), 1), c(rep(1, s), 0))
    sol <- tryCatch(solve(K, c(-e[S], 1)), error = function(err) NULL)
    if (is.null(sol)) break
    z <- sol[seq_len(s)]
    if (all(z >= 0)) {
      lam[] <- 0
      lam[S] <- z
      mu <- as.numeric(H %*% lam) + e + sol[s + 1L]
      out <- setdiff(seq_len(d), S)
      if (!length(out)) break
      j <- out[which.min(mu[out])]
      if (mu[j] >= -1e-14 * (1 + max(abs(mu)))) break
      S <- c(S, j)
    } else {
      dir <- z - lam[S]
      neg <- dir < 0
      a <- min(lam[S][neg] / -dir[neg])
      lam[S] <- lam[S] + a * dir
      lam[abs(lam) < 1e-16] <- 0
      S <- S[lam[S] > 0]
      if (!length(S)) { S <- which.min(diag(H) / 2 + e); lam[] <- 0; lam[S] <- 1 }
    }
  }
  lam <- pmax(lam, 0)
  lam / sum(lam)
}

.mo_convex_bundle <- function(mp, p, stopping_criterion = NULL, atol_lambda = sqrt(.Machine$double.eps),
                              bundle_cap = 25L, contraction_factor = 0.975, m = 1e-3, retraction_method = "exp",
                              vector_transport_method = "parallel", domain = NULL, ...) {
  M <- mp$M
  dom <- domain %||% function(M, q) is.finite(mp$f(q))
  sc <- stopping_criterion %||% list(.mo_sc("StopWhenLagrangeMultiplierLess", tol = 1e-8),
                                     .mo_sc("StopAfterIteration", n = 5000L))
  st <- .mo_state(mp, p, sc)
  st$pls <- p
  st$X <- mp$grad(p)
  st$xi <- 0
  bundle <- list(list(q = p, X = st$X))
  lam <- 0; errs <- 0
  st$last_stepsize <- 0
  .mo_run(st, function(k) {
    Tg <- lapply(bundle, function(b) .mf_transport(M, b$q, b$X, st$pls, vector_transport_method))
    Hm <- outer(seq_along(Tg), seq_along(Tg), Vectorize(function(i, j) .mf_inner(M, st$pls, Tg[[i]], Tg[[j]])))
    lam <<- .mo_qp_simplex(Hm, errs)
    g <- Reduce(`+`, Map(`*`, lam, Tg))
    eps <- sum(lam * errs)
    ng <- .mf_norm(M, st$pls, g)
    st$xi <- -ng^2 - eps
    # DomainBackTracking
    t <- 1
    cand <- .mf_retract(M, st$pls, -t * g, 1, retraction_method)
    bad <- function(q, t) {
      d <- .mf_distance(M, st$pls, q)
      !dom(M, q) || (d < t * ng && !isTRUE(all.equal(d, t * ng, tolerance = sqrt(.Machine$double.eps))))
    }
    while (bad(cand, t) && t > 0x1.56e1fc2f8f359p-997) {
      t <- t * contraction_factor
      cand <- .mf_retract(M, st$pls, -t * g, 1, retraction_method)
    }
    st$last_stepsize <- t
    st$p <- cand
    fls <- mp$f(st$pls)
    if (mp$f(st$p) <= fls + t * m * st$xi) {
      st$pls <- st$p
      st$X <- mp$grad(st$p)
    } else {
      # NullStepBackTracking
      cand <- .mf_retract(M, st$pls, -t * g, 1, retraction_method)
      Xc <- mp$grad(cand)
      null_cond <- function(q, Xq, t) {
        v <- .mf_invretract(M, q, st$pls)
        .mf_inner(M, st$pls, .mf_transport(M, q, Xq, st$pls, vector_transport_method), t * g) >=
          -m * t * st$xi - (fls - mp$f(q) - .mf_inner(M, q, Xq, v))
      }
      while (null_cond(cand, Xc, t) && t > 0x1.56e1fc2f8f359p-997) {
        t <- t * contraction_factor
        cand <- .mf_retract(M, st$pls, -t * g, 1, retraction_method)
        Xc <- mp$grad(cand)
      }
      st$p <- cand
      st$X <- Xc
    }
    keep <- lam > atol_lambda
    bundle <<- bundle[keep]; lam <<- lam[keep]; errs <<- errs[keep]
    while (length(bundle) >= bundle_cap) {
      j <- which(vapply(bundle, function(b) !identical(b$q, st$pls), TRUE))[1L]
      if (is.na(j)) break
      bundle <<- bundle[-j]; lam <<- lam[-j]; errs <<- errs[-j]
    }
    bundle <<- c(bundle, list(list(q = st$p, X = st$X)))
    lam <<- c(lam, 0)
    fls <- mp$f(st$pls)
    errs <<- vapply(bundle, function(b) fls - mp$f(b$q) - .mf_inner(M, b$q, b$X, .mf_invretract(M, b$q, st$pls)), 0)
  })
  st$p <- st$pls
  st
}

# ---------------------------------------------------------------------------
# Optimization.jl interface

.mo_opt <- function(name, grad = FALSE, hess = FALSE) {
  structure(list(name = name, traits = list(grad = grad, hess = hess, bounds = FALSE, cons = FALSE,
                                            callback = TRUE)),
            class = c(name, "AbstractManoptOptimizer", "OptimizationAlgorithm"))
}

#' Manopt optimizers on manifolds
#'
#' The optimizers of OptimizationManopt: Riemannian optimization with the
#' Manopt.jl algorithms, ported to R, on the manifold given as the
#' `manifold` of the [OptimizationProblem()] ([Euclidean()] or
#' [Sphere()]). Gradient based methods turn the Euclidean gradient of the
#' problem (finite differences when it has none) into the Riemannian one by
#' projecting it onto the tangent space; the trust region and cubic
#' regularization methods use the Riemannian Hessian built the same way
#' from Hessian vector products.
#'
#' * `GradientDescentOptimizer()`: gradient descent with an Armijo line
#'   search (Manopt.gradient_descent).
#' * `ConjugateGradientDescentOptimizer()`: nonlinear conjugate gradients
#'   with the conjugate descent coefficient and restarts on non descent
#'   directions.
#' * `QuasiNewtonOptimizer()`: limited memory Riemannian BFGS with a Wolfe
#'   Powell line search.
#' * `NelderMeadOptimizer()`: the Riemannian Nelder Mead simplex method; as
#'   upstream, the simplex is random and the starting point is not used.
#' * `ParticleSwarmOptimizer()`: particle swarm optimization with
#'   `population_size` particles (100), one of them the starting point.
#' * `CMAESOptimizer()`: the Riemannian covariance matrix adaptation
#'   evolution strategy.
#' * `ConvexBundleOptimizer()`: the convex bundle method for (sub)gradient
#'   problems.
#' * `FrankWolfeOptimizer()`: the Frank Wolfe method with a gradient descent
#'   oracle.
#' * `TrustRegionsOptimizer()`: Riemannian trust regions with a truncated
#'   conjugate gradient subsolver.
#' * `AdaptiveRegularizationCubicOptimizer()`: adaptive regularization with
#'   cubics and a Lanczos subsolver.
#'
#' The common `solve` arguments map as in Optimization.jl: `maxiters` to
#' `StopAfterIteration`, `maxtime` (seconds) to `StopAfter`, `abstol` to
#' `StopWhenGradientNormLess` for gradient based methods and to
#' `StopWhenChangeLess` otherwise; a `stopping_criterion` from [Manopt] is
#' added to them. When `maxiters` or `maxtime` is given without `abstol` or
#' a `stopping_criterion`, a convergence criterion is added (gradient norm
#' below 1e-8 for gradient descent and conjugate gradients, 1e-6 for quasi
#' Newton and trust regions, 1e-9 for cubic regularization, gradient norm
#' or change below 1e-8 for Frank Wolfe). Without any of them the solver
#' runs with its Manopt defaults. Other `solve` arguments reach the solver:
#' `stepsize`, `retraction_method`, `vector_transport_method`,
#' `population_size`, `sigma` and so on. Constraints and bounds are not
#' supported.
#'
#' @return An algorithm object.
#' @examples
#' rosenbrock <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
#' prob <- OptimizationProblem(rosenbrock, c(0, 0), c(1, 100), manifold = Euclidean(2))
#' sol <- solve(prob, QuasiNewtonOptimizer())
#' round(sol$u, 4)
#' # the smallest eigenvector of a symmetric matrix on the sphere
#' A <- diag(c(3, 1, 2))
#' rq <- function(x, p) sum(x * (A %*% x))
#' sph <- OptimizationProblem(rq, c(1, 1, 1) / sqrt(3), manifold = Sphere(2))
#' round(abs(solve(sph, TrustRegionsOptimizer())$u), 6)
#' @name manopt_optimizers
#' @export
AbstractManoptOptimizer <- function() .mo_opt("AbstractManoptOptimizer")
#' @rdname manopt_optimizers
#' @export
GradientDescentOptimizer <- function() .mo_opt("GradientDescentOptimizer", grad = TRUE)
#' @rdname manopt_optimizers
#' @export
ConjugateGradientDescentOptimizer <- function() .mo_opt("ConjugateGradientDescentOptimizer", grad = TRUE)
#' @rdname manopt_optimizers
#' @export
QuasiNewtonOptimizer <- function() .mo_opt("QuasiNewtonOptimizer", grad = TRUE)
#' @rdname manopt_optimizers
#' @export
NelderMeadOptimizer <- function() .mo_opt("NelderMeadOptimizer")
#' @rdname manopt_optimizers
#' @export
ParticleSwarmOptimizer <- function() .mo_opt("ParticleSwarmOptimizer")
#' @rdname manopt_optimizers
#' @export
CMAESOptimizer <- function() .mo_opt("CMAESOptimizer")
#' @rdname manopt_optimizers
#' @export
ConvexBundleOptimizer <- function() .mo_opt("ConvexBundleOptimizer", grad = TRUE)
#' @rdname manopt_optimizers
#' @export
FrankWolfeOptimizer <- function() .mo_opt("FrankWolfeOptimizer", grad = TRUE)
#' @rdname manopt_optimizers
#' @export
TrustRegionsOptimizer <- function() .mo_opt("TrustRegionsOptimizer", grad = TRUE, hess = TRUE)
#' @rdname manopt_optimizers
#' @export
AdaptiveRegularizationCubicOptimizer <- function() .mo_opt("AdaptiveRegularizationCubicOptimizer", grad = TRUE,
                                                           hess = TRUE)

.mo_gradient_based <- c("GradientDescentOptimizer", "ConjugateGradientDescentOptimizer", "QuasiNewtonOptimizer",
                        "ConvexBundleOptimizer", "FrankWolfeOptimizer", "AdaptiveRegularizationCubicOptimizer",
                        "TrustRegionsOptimizer")

.mo_fallback_tol <- c(GradientDescentOptimizer = 1e-8, ConjugateGradientDescentOptimizer = 1e-8,
                      QuasiNewtonOptimizer = 1e-6, TrustRegionsOptimizer = 1e-6,
                      AdaptiveRegularizationCubicOptimizer = 1e-9, FrankWolfeOptimizer = 1e-8)

# the problem of a manifold optimization: cost, Riemannian gradient and
# Riemannian Hessian from the Euclidean ones
.mo_problem <- function(M, f, eg, ehv) {
  list(M = M, f = f,
       grad = function(x) .mf_project(M, x, eg(x)),
       hess = function(x, v) {
         Y <- .mf_project(M, x, ehv(x, v))
         if (!.mf_sphere(M)) return(Y)
         G <- eg(x)
         Y + (-.jl_dot(x, G - .mf_project(M, x, G))) * v
       })
}

.mo_solve <- function(a, alg, opts) {
  .opt_check_traits(a, alg)
  M <- opts$manifold
  if (is.null(M))
    stop("Manifold not specified in the problem, for example ",
         "`OptimizationProblem(f, x, p, manifold = Sphere(2))`.", call. = FALSE)
  if (!inherits(M, "AbstractManifold")) stop("`manifold` must be Euclidean(n) or Sphere(n)", call. = FALSE)
  x0 <- as.numeric(a$u0)
  if (length(x0) != .mf_rep(M))
    stop("the starting point has length ", length(x0), " but points of ", class(M)[1L], "(", M$n,
         ") have length ", .mf_rep(M), call. = FALSE)
  p <- a$p; fn <- a$f
  f <- function(x) as.numeric(fn$f(x, p))[1L]
  eg <- if (is.function(fn$grad)) function(x) as.numeric(fn$grad(x, p)) else function(x) .opt_fd_gradient(f, x)
  ehv <- if (is.function(fn$hess)) function(x, v) as.numeric(fn$hess(x, p) %*% v) else function(x, v) {
    nv <- sqrt(sum(v * v))
    if (nv == 0) return(numeric(length(x)))
    h <- .Machine$double.eps^(1 / 3) * max(1, sqrt(sum(x * x))) / nv
    (eg(x + h * v) - eg(x - h * v)) / (2 * h)
  }
  mp <- .mo_problem(M, f, eg, ehv)
  name <- alg$name
  crit <- list()
  if (!is.null(opts$maxiters)) crit <- c(crit, list(.mo_sc("StopAfterIteration", n = as.integer(opts$maxiters))))
  if (!is.null(opts$maxtime)) crit <- c(crit, list(.mo_sc("StopAfter", seconds = opts$maxtime)))
  if (!is.null(opts$stopping_criterion)) crit <- c(crit, .mo_flatten_sc(opts$stopping_criterion))
  if (!is.null(opts$abstol)) {
    crit <- c(crit, list(if (name %in% .mo_gradient_based) .mo_sc("StopWhenGradientNormLess", tol = opts$abstol) else
      .mo_sc("StopWhenChangeLess", tol = opts$abstol)))
  } else if (is.null(opts$stopping_criterion) && (!is.null(opts$maxiters) || !is.null(opts$maxtime)) &&
             name %in% names(.mo_fallback_tol)) {
    crit <- c(crit, list(.mo_sc("StopWhenGradientNormLess", tol = .mo_fallback_tol[[name]])))
    if (name == "FrankWolfeOptimizer") crit <- c(crit, list(.mo_sc("StopWhenChangeLess", tol = 1e-8)))
  }
  args <- opts[setdiff(names(opts), c("maxiters", "maxtime", "abstol", "reltol", "stopping_criterion", "manifold",
                                      "callback", "verbose", "population_size", "population"))]
  args$stopping_criterion <- if (length(crit)) crit
  st <- switch(name,
    GradientDescentOptimizer = do.call(.mo_gradient_descent, c(list(mp, x0), args)),
    ConjugateGradientDescentOptimizer = do.call(.mo_conjugate_gradient, c(list(mp, x0), args)),
    QuasiNewtonOptimizer = do.call(.mo_quasi_newton, c(list(mp, x0), args)),
    TrustRegionsOptimizer = do.call(.mo_trust_regions, c(list(mp, x0), args)),
    AdaptiveRegularizationCubicOptimizer = do.call(.mo_arc, c(list(mp, x0), args)),
    FrankWolfeOptimizer = do.call(.mo_frank_wolfe, c(list(mp, x0), args)),
    NelderMeadOptimizer = do.call(.mo_nelder_mead, c(list(mp, opts$population), args)),
    ParticleSwarmOptimizer = {
      np <- as.integer(opts$population_size %||% 100L)
      swarm <- c(list(x0), lapply(seq_len(np - 1L), function(i) .mf_rand(M)))
      do.call(.mo_particle_swarm, c(list(mp, swarm), args))
    },
    CMAESOptimizer = do.call(.mo_cma_es, c(list(mp, x0), args)),
    ConvexBundleOptimizer = do.call(.mo_convex_bundle, c(list(mp, x0), args)),
    stop("AbstractManoptOptimizer() is the abstract type of the Manopt optimizers; use one of them, ",
         "for example GradientDescentOptimizer()", call. = FALSE))
  .opt_solution(a, alg, st$p, f(st$p), .mo_retcode(st$stopper$active), st)
}

.mo_retcode <- function(act) {
  if (is.null(act)) return("Failure")
  if (isTRUE(.mo_sc_converges[act$type])) return("Success")
  switch(act$type, StopWhenChangeLess = "Success", StopAfterIteration = "MaxIters", StopAfter = "MaxTime",
         StopWhenStepsizeLess = "Stalled", "Failure")
}

.mtk$hooks$opt_solve_AbstractManoptOptimizer <- function(a, alg, opts) .mo_solve(a, alg, opts)

# ---------------------------------------------------------------------------
# The Manopt module

.mo_user_problem <- function(M, f, grad_f = NULL, Hess_f = NULL) {
  cost <- function(p) f(M, p)
  g <- if (is.function(grad_f)) function(p) grad_f(M, p) else function(p) .mf_project(M, p, .opt_fd_gradient(cost, p))
  h <- if (is.function(Hess_f)) function(p, X) Hess_f(M, p, X) else function(p, X) {
    nx <- .mf_norm(M, p, X)
    if (nx == 0) return(numeric(length(p)))
    c <- 2^-14 / nx
    q <- .mf_retract(M, p, X, c)
    (.mf_transport(M, q, g(q), p) - g(p)) / c
  }
  list(M = M, f = cost, grad = g, hess = h)
}

.mo_return <- function(st, return_state) {
  if (isTRUE(return_state)) {
    act <- st$stopper$active
    return(structure(list(p = st$p, iterations = st$iterations, stopped_by = act,
                          converged = !is.null(act) && isTRUE(.mo_sc_converges[act$type])),
                     class = "ManoptSolverState"))
  }
  st$p
}

#' @export
print.ManoptSolverState <- function(x, ...) {
  cat("Manopt solver state after ", x$iterations, " iterations\n", sep = "")
  if (!is.null(x$stopped_by)) cat("stopped by ", x$stopped_by$type, "\n", sep = "")
  invisible(x)
}

#' The Manopt module
#'
#' `Manopt` gathers the parts of Manopt.jl that OptimizationManopt
#' re-exports, ported to R: the solvers on a manifold `M` ([Euclidean()] or
#' [Sphere()]), called with the cost `f(M, p)` and, where needed, the
#' Riemannian gradient `grad_f(M, p)` and Hessian `Hess_f(M, p, X)`
#' (approximated by finite differences of the gradient when missing), the
#' stopping criteria and the step sizes. Reach them with `Manopt$name`.
#'
#' Solvers: `gradient_descent(M, f, grad_f, p, ...)`,
#' `conjugate_gradient_descent(M, f, grad_f, p, ...)`,
#' `quasi_Newton(M, f, grad_f, p, ...)`,
#' `trust_regions(M, f, grad_f, Hess_f, p, ...)`,
#' `adaptive_regularization_with_cubics(M, f, grad_f, Hess_f, p, ...)`,
#' `Frank_Wolfe_method(M, f, grad_f, p, ...)`,
#' `NelderMead(M, f, population, ...)`, `particle_swarm(M, f, swarm, ...)`,
#' `cma_es(M, f, p, ...)` and `convex_bundle_method(M, f, subgrad_f, p, ...)`.
#' Each returns the final point, or with `return_state = TRUE` a state
#' whose `p` is that point and `stopped_by` the criterion that stopped it;
#' `get_solver_result(state)` returns the point. Keyword arguments include
#' `stopping_criterion`, `stepsize`, `retraction_method` and
#' `vector_transport_method`.
#'
#' Stopping criteria: `StopAfterIteration(n)`, `StopAfter(seconds)`,
#' `StopWhenGradientNormLess(tol)`, `StopWhenChangeLess(tol)`,
#' `StopWhenStepsizeLess(tol)`, `StopWhenCostLess(tol)`,
#' `StopWhenPopulationConcentrated(tol_f, tol_p)`,
#' `StopWhenSwarmVelocityLess(tol)`, `StopWhenLagrangeMultiplierLess(tol)`
#' and `StopWhenAny(...)`, which stops when any of its criteria does.
#'
#' Step sizes: `ArmijoLinesearch(...)`, `WolfePowellLinesearch(...)`,
#' `DecreasingStepsize(...)` and `ConstantLength(...)`. Retractions and
#' transports: `ExponentialRetraction()`, `ProjectionRetraction()`,
#' `ParallelTransport()` and `ProjectionTransport()`.
#'
#' @format A list of functions.
#' @examples
#' M <- Sphere(2)
#' A <- diag(c(2, 1, 3))
#' f <- function(M, p) sum(p * (A %*% p))
#' grad_f <- function(M, p) { g <- 2 * as.numeric(A %*% p); g - sum(p * g) * p }
#' p <- Manopt$gradient_descent(M, f, grad_f, c(1, 1, 1) / sqrt(3),
#'                              stopping_criterion = Manopt$StopAfterIteration(50))
#' round(abs(p), 6)
#' @export
Manopt <- list(
  gradient_descent = function(M, f, grad_f, p, ..., return_state = FALSE)
    .mo_return(.mo_gradient_descent(.mo_user_problem(M, f, grad_f), p, ...), return_state),
  conjugate_gradient_descent = function(M, f, grad_f, p, ..., return_state = FALSE)
    .mo_return(.mo_conjugate_gradient(.mo_user_problem(M, f, grad_f), p, ...), return_state),
  quasi_Newton = function(M, f, grad_f, p, ..., return_state = FALSE)
    .mo_return(.mo_quasi_newton(.mo_user_problem(M, f, grad_f), p, ...), return_state),
  trust_regions = function(M, f, grad_f, Hess_f = NULL, p, ..., return_state = FALSE)
    .mo_return(.mo_trust_regions(.mo_user_problem(M, f, grad_f, Hess_f), p, ...), return_state),
  adaptive_regularization_with_cubics = function(M, f, grad_f, Hess_f = NULL, p, ..., return_state = FALSE)
    .mo_return(.mo_arc(.mo_user_problem(M, f, grad_f, Hess_f), p, ...), return_state),
  Frank_Wolfe_method = function(M, f, grad_f, p, ..., return_state = FALSE)
    .mo_return(.mo_frank_wolfe(.mo_user_problem(M, f, grad_f), p, ...), return_state),
  NelderMead = function(M, f, population = NULL, ..., return_state = FALSE)
    .mo_return(.mo_nelder_mead(.mo_user_problem(M, f), population, ...), return_state),
  particle_swarm = function(M, f, swarm = NULL, ..., return_state = FALSE)
    .mo_return(.mo_particle_swarm(.mo_user_problem(M, f), swarm %||% lapply(1:100, function(i) .mf_rand(M)), ...),
               return_state),
  cma_es = function(M, f, p = .mf_rand(M), ..., return_state = FALSE)
    .mo_return(.mo_cma_es(.mo_user_problem(M, f), p, ...), return_state),
  convex_bundle_method = function(M, f, subgrad_f, p, ..., return_state = FALSE)
    .mo_return(.mo_convex_bundle(.mo_user_problem(M, f, subgrad_f), p, ...), return_state),
  get_solver_result = function(state) if (inherits(state, "ManoptSolverState")) state$p else state,
  StopAfterIteration = function(n) .mo_sc("StopAfterIteration", n = as.integer(n)),
  StopAfter = function(seconds) .mo_sc("StopAfter", seconds = seconds),
  StopWhenGradientNormLess = function(tol) .mo_sc("StopWhenGradientNormLess", tol = tol),
  StopWhenChangeLess = function(tol) .mo_sc("StopWhenChangeLess", tol = tol),
  StopWhenStepsizeLess = function(tol) .mo_sc("StopWhenStepsizeLess", tol = tol),
  StopWhenCostLess = function(tol) .mo_sc("StopWhenCostLess", tol = tol),
  StopWhenPopulationConcentrated = function(tol_f = 1e-8, tol_p = 1e-8)
    .mo_sc("StopWhenPopulationConcentrated", tol_f = tol_f, tol_p = tol_p),
  StopWhenSwarmVelocityLess = function(tol) .mo_sc("StopWhenSwarmVelocityLess", tol = tol),
  StopWhenLagrangeMultiplierLess = function(tol = 1e-6) .mo_sc("StopWhenLagrangeMultiplierLess", tol = tol),
  StopWhenAny = function(...) .mo_sc("any", criteria = list(...)),
  ArmijoLinesearch = function(...) .manopt_armijo(...),
  WolfePowellLinesearch = function(...) .manopt_wolfe_powell(...),
  DecreasingStepsize = function(...) .manopt_decreasing(...),
  ConstantLength = function(length = NULL, type = c("relative", "absolute"))
    .mo_ss("Constant", length = length, kind = match.arg(type)),
  ExponentialRetraction = function() "exp",
  ProjectionRetraction = function() "project",
  ParallelTransport = function() "parallel",
  ProjectionTransport = function() "project"
)

.manopt_armijo <- function(contraction_factor = 0.95, sufficient_decrease = 0.1, initial_stepsize = 1,
                           retraction_method = NULL, stop_when_stepsize_less = 0, stop_when_stepsize_exceeds = NULL,
                           stop_increasing_at_step = 100L, stop_decreasing_at_step = 1000L) {
  .mo_ss("Armijo", contraction_factor = contraction_factor, sufficient_decrease = sufficient_decrease,
         initial_stepsize = initial_stepsize, retraction_method = retraction_method,
         stop_when_stepsize_less = stop_when_stepsize_less, stop_when_stepsize_exceeds = stop_when_stepsize_exceeds,
         stop_increasing_at_step = as.integer(stop_increasing_at_step),
         stop_decreasing_at_step = as.integer(stop_decreasing_at_step))
}

.manopt_wolfe_powell <- function(sufficient_decrease = 1e-4, sufficient_curvature = 0.999, initial_guess = 1,
                                 max_stepsize = NULL, retraction_method = NULL, vector_transport_method = NULL,
                                 stop_when_stepsize_less = 0, stop_increasing_at_step = 100L,
                                 stop_decreasing_at_step = 1000L) {
  .mo_ss("WolfePowell", sufficient_decrease = sufficient_decrease, sufficient_curvature = sufficient_curvature,
         initial_guess = initial_guess, max_stepsize = max_stepsize, retraction_method = retraction_method,
         vector_transport_method = vector_transport_method, stop_when_stepsize_less = stop_when_stepsize_less,
         stop_increasing_at_step = as.integer(stop_increasing_at_step),
         stop_decreasing_at_step = as.integer(stop_decreasing_at_step))
}

.manopt_decreasing <- function(length = 1, factor = 1, subtrahend = 0, exponent = 1, shift = 0,
                               type = c("relative", "absolute")) {
  .mo_ss("Decreasing", length = length, factor = factor, subtrahend = subtrahend, exponent = exponent, shift = shift,
         kind = match.arg(type))
}
