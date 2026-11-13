# Optimization on Riemannian manifolds: the OptimizationManopt wrapper of
# Optimization.jl and the Manopt.jl (0.6) solvers it calls, on the Euclidean
# and Sphere manifolds of Manifolds.jl.

# ---------------------------------------------------------------------------
# Floating point as in Julia: LinearAlgebra.dot of short vectors runs the BLAS
# kernel, a sequence of fused multiply adds; LinearAlgebra.norm of a short
# vector is the square root of a plain running sum of squares.

.jl_dot <- function(x, y) {
  n <- length(x)
  if (n > 16L) return(sum(x * y))
  p <- x * y
  sa <- 134217729 * x; ah <- sa - (sa - x); al <- x - ah
  sb <- 134217729 * y; bh <- sb - (sb - y); bl <- y - bh
  e <- ((ah * bh - p) + ah * bl + al * bh) + al * bl
  s <- 0
  for (i in seq_len(n)) {
    t <- p[i] + s
    bb <- t - p[i]
    s <- t + (((p[i] - (t - bb)) + (s - bb)) + e[i])
  }
  s
}

# eps(x) of Julia: the spacing of doubles at x
.jl_eps <- function(x) {
  x <- abs(x)
  if (x < 0x1.0000000000000p-1022) return(4.9406564584124654e-324)
  2^(floor(log2(x)) - 52)
}

.jl_norm <- function(x) {
  if (!length(x)) return(0)
  m <- max(abs(x))
  if (is.na(m) || m == 0 || is.infinite(m)) return(m)
  if (is.finite(length(x) * m * m) && m * m != 0) return(sqrt(sum(x * x)))
  m * sqrt(sum((x / m) * (x / m)))
}

# ---------------------------------------------------------------------------
# Manifolds

#' Manifolds for optimization
#'
#' `Euclidean(n)` is the space \eqn{R^n} and `Sphere(n)` the unit sphere
#' \eqn{S^n} embedded in \eqn{R^{n+1}}, the manifolds of Manifolds.jl that
#' the Manopt optimizers ([GradientDescentOptimizer()] and the others)
#' work on. Pass one as the `manifold` of an [OptimizationProblem()]; on the
#' sphere the starting point must have norm one. The sphere carries the
#' geometry of Manifolds.jl: the exponential map followed by a projection
#' back onto the sphere as its retraction, the parallel transport, the
#' logarithmic map and the projection onto tangent spaces, through which a
#' Euclidean gradient becomes a Riemannian one.
#'
#' `SymmetricPositiveDefinite(n)` is the manifold of symmetric positive
#' definite \eqn{n \times n} matrices with the affine invariant metric
#' \eqn{g_p(X, Y) = tr(p^{-1} X p^{-1} Y)}: its points are matrices (the
#' starting point and the solution of an [OptimizationProblem()]), its
#' retraction the exponential map \eqn{p^{1/2} Exp(p^{-1/2} X p^{-1/2})
#' p^{1/2}}, with the logarithmic map, the parallel transport, the
#' projection onto the symmetric matrices and the Riemannian gradient
#' \eqn{p sym(G) p} of a Euclidean gradient \eqn{G}. `Lorentz(n)` is
#' \eqn{R^n} with the Minkowski metric (timelike in the last coordinate), as
#' in Manifolds.jl: its exponential map, transport and retraction are those
#' of \eqn{R^n} and its Riemannian gradient flips the sign of the last
#' coordinate. Both are also the manifolds of the geodesic convexity
#' analysis, `analyze(ex, M)` ([analyze()]).
#'
#' @param n The dimension of the manifold.
#' @return A manifold object.
#' @examples
#' M <- Sphere(2)
#' M
#' f <- function(x, p) x[1] + 2 * x[2] + 3 * x[3]
#' prob <- OptimizationProblem(f, c(1, 0, 0), manifold = M)
#' sol <- solve(prob, GradientDescentOptimizer(), maxiters = 100)
#' round(sol$u, 4)
#' @name manifolds
#' @export
Euclidean <- function(n) structure(list(n = .mf_check_dim(n)), class = c("Euclidean", "AbstractManifold"))
#' @rdname manifolds
#' @export
Sphere <- function(n) structure(list(n = .mf_check_dim(n)), class = c("Sphere", "AbstractManifold"))

.mf_check_dim <- function(n) {
  n <- suppressWarnings(as.integer(n))
  if (length(n) != 1L || is.na(n) || n < 1L) stop("the dimension must be a positive integer", call. = FALSE)
  n
}

#' @export
print.AbstractManifold <- function(x, ...) {
  cat(class(x)[1L], "(", x$n, ")\n", sep = "")
  invisible(x)
}

.mf_sphere <- function(M) inherits(M, "Sphere")
.mf_spd <- function(M) inherits(M, "SymmetricPositiveDefinite")
.mf_lorentz <- function(M) inherits(M, "Lorentz")
# manifold_dimension
.mf_dim <- function(M) if (.mf_spd(M)) M$n * (M$n + 1L) %/% 2L else M$n
# the length of a point as a vector (a positive definite matrix by columns)
.mf_rep <- function(M) if (.mf_sphere(M)) M$n + 1L else if (.mf_spd(M)) M$n * M$n else M$n
# a point or tangent vector as the user sees it
.mf_shape <- function(M, x) if (.mf_spd(M)) matrix(x, M$n, M$n) else x
.mf_inj <- function(M) if (.mf_sphere(M)) pi else Inf
.mf_inner <- function(M, p, X, Y) {
  if (.mf_spd(M)) return(.spd_inner(M$n, p, X, Y))
  if (.mf_lorentz(M)) return(.mf_minkowski(X, Y))
  .jl_dot(X, Y)
}
.mf_norm <- function(M, p, X) {
  if (.mf_spd(M)) return(sqrt(.spd_inner(M$n, p, X, X)))
  if (.mf_lorentz(M)) return(.mf_sqrt(.mf_minkowski(X, X)))
  .jl_norm(X)
}
.mf_project <- function(M, p, X) {
  if (.mf_sphere(M)) return(X - .jl_dot(p, X) * p)
  if (.mf_spd(M)) return(as.numeric(.spd_sym(matrix(X, M$n))))
  X
}
.mf_project_point <- function(M, p) if (.mf_sphere(M)) p / .jl_norm(p) else p
.mf_cbrt_eps <- 6.055454452393343e-06

# the Minkowski metric of Manifolds.jl, timelike in the last coordinate
.mf_minkowski <- function(a, b) {
  n <- length(a)
  -a[n] * b[n] + .jl_dot(a[-n], b[-n])
}
# Julia's sqrt throws on a negative argument (a timelike vector has no norm)
.mf_sqrt <- function(x) {
  if (x < 0) stop("DomainError with ", .jl_show_float(x), ":\nsqrt was called with a negative real argument but ",
                  "will only return a complex result if called with a complex argument. Try sqrt(Complex(x)).",
                  call. = FALSE)
  sqrt(x)
}

# The affine invariant geometry of SymmetricPositiveDefinite(n) as Manifolds.jl
# computes it (HermitianPositiveDefiniteAffineInvariant.jl): a point and a
# tangent vector are n x n matrices, stored by columns.
.spd_sym <- function(X) (X + t(X)) / 2
# Hermitian(A): the matrix its upper triangle defines
.spd_herm <- function(A) {
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  A
}
.spd_eigen <- function(A) {
  e <- eigen(.spd_herm(A), symmetric = TRUE)
  o <- order(e$values)
  list(values = e$values[o], vectors = e$vectors[, o, drop = FALSE])
}
.spd_fun <- function(U, d) U %*% (d * t(U))
# hpd_sqrt_and_sqrt_inv
.spd_sqrt2 <- function(p) {
  e <- .spd_eigen(p)
  S <- pmax(e$values, .Machine$double.xmin)
  list(sqrt = .spd_herm(.spd_fun(e$vectors, sqrt(S))), inv = .spd_herm(.spd_fun(e$vectors, 1 / sqrt(S))))
}
.spd_inner <- function(n, p, X, Y) {
  p <- matrix(p, n); X <- matrix(X, n); Y <- matrix(Y, n)
  R <- chol(.spd_herm(p))
  sum(backsolve(R, forwardsolve(t(R), .spd_herm(X))) * t(backsolve(R, forwardsolve(t(R), t(.spd_herm(Y))))))
}
.spd_exp <- function(n, p, X) {
  r <- .spd_sqrt2(matrix(p, n))
  e <- .spd_eigen(r$inv %*% matrix(X, n) %*% r$inv)
  pU <- r$sqrt %*% e$vectors
  as.numeric(pU %*% (exp(e$values) * t(pU)))
}
.spd_log <- function(n, p, q) {
  r <- .spd_sqrt2(matrix(p, n))
  e <- .spd_eigen(r$inv %*% matrix(q, n) %*% r$inv)
  a <- max(abs(e$values))
  c <- if (a < .Machine$double.eps) max(.Machine$double.eps / 2 * a, .Machine$double.xmin) else .Machine$double.eps
  pU <- r$sqrt %*% e$vectors
  as.numeric(pU %*% (log(pmax(e$values, c)) * t(pU)))
}
.spd_distance <- function(n, p, q) {
  if (identical(as.numeric(p), as.numeric(q))) return(0)
  U <- chol(.spd_herm(matrix(q, n)))
  s <- .spd_eigen(t(forwardsolve(t(U), t(forwardsolve(t(U), matrix(p, n))))))$values
  sqrt(sum(abs(log(pmax(s, .Machine$double.xmin)))^2))
}
.spd_transport <- function(n, p, X, q) {
  if (.spd_distance(n, p, q) < 2 * .Machine$double.eps) return(X)
  r <- .spd_sqrt2(matrix(p, n))
  tv <- .spd_herm(r$inv %*% matrix(X, n) %*% r$inv)
  e2 <- .spd_eigen(r$inv %*% matrix(q, n) %*% r$inv)
  logty <- .spd_herm(.spd_fun(e2$vectors, log(pmax(e2$values, .Machine$double.xmin))))
  e3 <- .spd_eigen(logty)
  pU <- r$sqrt %*% .spd_fun(e3$vectors, exp(e3$values / 2))
  as.numeric(.spd_herm(pU %*% tv %*% t(pU)))
}
# the orthonormal basis of the tangent space at p (get_basis_orthonormal),
# and the basis that diagonalizes the curvature along the direction v with
# its sectional curvatures (get_basis_diagonalizing)
.spd_basis <- function(n, p, v = NULL) {
  r <- .spd_sqrt2(matrix(p, n))
  V <- diag(n); lam <- NULL
  if (!is.null(v)) {
    e <- .spd_eigen(r$inv %*% matrix(v, n) %*% r$inv)
    V <- e$vectors; lam <- e$values
  }
  B <- list(); k <- numeric()
  for (i in seq_len(n)) for (j in i:n) {
    E <- V[, i] %*% t(V[, j]) + V[, j] %*% t(V[, i])
    B[[length(B) + 1L]] <- as.numeric((if (i == j) 1 / 2 else 1 / sqrt(2)) * (r$sqrt %*% E %*% r$sqrt))
    k <- c(k, if (is.null(lam)) 0 else -1 / 4 * (lam[i] - lam[j])^2)
  }
  list(vectors = B, kappa = k)
}

.mf_cos_usinc_sq <- function(t2) {
  if (t2 < .mf_cbrt_eps) return(c(1 - t2 * (0.5 - t2 / 24), 1 - t2 * (1 / 6 - t2 / 120)))
  th <- sqrt(t2)
  c(cos(th), sin(th) / th)
}

.mf_exp <- function(M, p, X, t = 1) {
  if (.mf_spd(M)) return(.spd_exp(M$n, p, t * X))
  if (!.mf_sphere(M)) return(p + t * X)
  cs <- .mf_cos_usinc_sq((t * t) * .jl_dot(X, X))
  cs[1L] * p + (cs[2L] * t) * X
}

# retract_fused(M, p, X, t, method): "exp" is the default of each manifold
# (on the sphere the exponential map followed by a projection),
# "project" the projection retraction.
.mf_retract <- function(M, p, X, t = 1, method = "exp") {
  if (.mf_spd(M)) {
    if (!identical(method, "exp")) .mf_no_method(M, "retraction", method)
    return(.spd_exp(M$n, p, t * X))
  }
  if (!.mf_sphere(M)) return(p + t * X)
  q <- if (identical(method, "project")) p + t * X else .mf_exp(M, p, X, t)
  q / .jl_norm(q)
}

.mf_no_method <- function(M, what, method)
  stop("MethodError: no ", what, " method \"", method, "\" on ", class(M)[1L], "(", M$n, ")", call. = FALSE)

.mf_log <- function(M, p, q) {
  if (.mf_spd(M)) return(.spd_log(M$n, p, q))
  if (!.mf_sphere(M)) return(q - p)
  ct <- min(max(.jl_dot(p, q), -1), 1)
  X <- if (ct >= 0) q - p else q + p
  X <- X - .jl_dot(p, X) * p
  s2 <- .jl_dot(X, X)
  if (ct >= 0) {
    if (s2 < .mf_cbrt_eps) {
      X <- X * (1 + s2 * (1 / 6 + s2 * (3 / 40)))
    } else {
      st <- sqrt(s2)
      X <- X * (atan2(st, ct) / st)
    }
    return(.mf_project(M, p, X))
  }
  X <- X - .jl_dot(p, X) * p
  st <- .jl_norm(X)
  if (st == 0) {
    X <- numeric(length(p))
    if (isTRUE(all.equal(abs(p[1L]), 1, tolerance = sqrt(.Machine$double.eps)))) X[2L] <- 1 else X[1L] <- 1
    X <- X - .jl_dot(p, X) * p
    X <- X * (pi / .jl_norm(X))
  } else {
    X <- X * ((pi - atan2(st, -ct)) / st)
  }
  .mf_project(M, p, X)
}

.mf_invretract <- function(M, p, q, method = "log") {
  if (.mf_spd(M)) {
    if (!identical(method, "log")) .mf_no_method(M, "inverse retraction", method)
    return(.spd_log(M$n, p, q))
  }
  if (!.mf_sphere(M)) return(q - p)
  if (identical(method, "project")) return(q / .jl_dot(p, q) - p)
  .mf_log(M, p, q)
}

.mf_distance <- function(M, p, q) {
  if (.mf_spd(M)) return(.spd_distance(M$n, p, q))
  # the norm of the logarithmic map of the Minkowski metric
  if (.mf_lorentz(M)) return(.mf_sqrt(.mf_minkowski(q - p, q - p)))
  if (!.mf_sphere(M)) return(sqrt(sum((p - q) * (p - q))))
  ct <- .jl_dot(p, q)
  if (-1 < ct && ct < 1 - sqrt(.Machine$double.eps) / 8) return(acos(ct))
  2 * abs(atan2(.jl_norm(p - q), .jl_norm(p + q)))
}

# vector_transport_to(M, p, X, q, method): "parallel" (the default) or
# "project".
.mf_transport <- function(M, p, X, q, method = "parallel") {
  if (.mf_spd(M)) {
    if (identical(method, "project")) return(.mf_project(M, q, X))
    return(.spd_transport(M$n, p, X, q))
  }
  if (!.mf_sphere(M)) return(X)
  if (identical(method, "project")) return(X - .jl_dot(q, X) * q)
  m <- p + q
  X - m * (2 * .jl_dot(X, q) / .jl_dot(m, m))
}

.mf_rand <- function(M) {
  if (.mf_spd(M)) {
    # rand(M): the exponential at the identity of a Gaussian symmetric
    # tangent vector of deviation 1/sqrt(n)
    n <- M$n
    X <- .spd_sym(matrix(stats::rnorm(n * n), n)) / sqrt(n)
    return(.spd_exp(n, diag(n), X))
  }
  v <- stats::rnorm(.mf_rep(M))
  .mf_project_point(M, v)
}
.mf_rand_vector <- function(M, p, sigma = 1) {
  if (.mf_spd(M))
    return(.mf_get_vector(M, p, sigma * stats::rnorm(.mf_dim(M))))
  .mf_project(M, p, sigma * stats::rnorm(.mf_rep(M)))
}

# coordinates in the default orthonormal basis of the tangent space
.mf_get_vector <- function(M, p, c) {
  if (.mf_spd(M)) return(Reduce(`+`, Map(function(b, ck) ck * b, .spd_basis(M$n, p)$vectors, c)))
  if (!.mf_sphere(M)) return(c)
  n <- M$n
  ct <- abs(p[1L])
  lam <- if (ct == 0) 1 else p[1L] / ct
  pend <- p[2:(n + 1L)]
  pX <- .jl_dot(pend, c)
  c(-lam * pX, c - pend * (pX / (1 + ct)))
}
.mf_get_coordinates <- function(M, p, X) {
  if (.mf_spd(M)) return(vapply(.spd_basis(M$n, p)$vectors, function(b) .spd_inner(M$n, p, X, b), 0))
  if (!.mf_sphere(M)) return(X)
  n <- M$n
  ct <- abs(p[1L])
  lam <- if (ct == 0) 1 else p[1L] / ct
  X[2:(n + 1L)] - p[2:(n + 1L)] * (lam * X[1L] / (1 + ct))
}

# the Riemannian mean by geodesic interpolation
.mf_mean <- function(M, pts) {
  if (!.mf_sphere(M) && !.mf_spd(M)) return(Reduce(`+`, pts) / length(pts))
  q <- pts[[1L]]
  for (i in seq_along(pts)[-1L]) q <- .mf_retract(M, q, .mf_log(M, q, pts[[i]]), 1 / i)
  # points beyond the radius pi/2: gradient descent estimation from there
  if (.mf_spd(M) || !any(vapply(pts, function(x) .mf_distance(M, q, x) >= pi / 2, TRUE))) return(q)
  for (it in seq_len(100L)) {
    yold <- q
    v <- .mf_log(M, yold, pts[[1L]])
    for (j in seq_along(pts)[-1L]) v <- v + (1 / j) * (.mf_log(M, yold, pts[[j]]) - v)
    q <- .mf_retract(M, yold, v, 0.5)
    if (.jl_norm(q - yold) <= sqrt(.Machine$double.eps) * max(.jl_norm(q), .jl_norm(yold))) break
  }
  q
}

# riemannian_gradient and riemannian_Hessian of ManifoldDiff and Manifolds:
# the projection of the Euclidean gradient G, then the change of representer
# to the metric of M; the Hessian from G and the Euclidean Hessian H at X
.mf_rgrad <- function(M, p, G) {
  if (.mf_spd(M)) {
    P <- matrix(p, M$n)
    return(as.numeric(P %*% .spd_sym(matrix(G, M$n)) %*% P))
  }
  if (.mf_lorentz(M)) return(c(G[-length(G)], -G[length(G)]))
  .mf_project(M, p, G)
}
.mf_rhess <- function(M, p, G, H, X) {
  if (.mf_spd(M)) {
    n <- M$n
    P <- matrix(p, n)
    Y <- .spd_sym(matrix(X, n) %*% .spd_sym(matrix(G, n)) %*% P)
    return(as.numeric(Y + P %*% .spd_sym(matrix(H, n)) %*% P))
  }
  if (.mf_lorentz(M)) return(c(H[-length(H)], -H[length(H)]))
  Y <- .mf_project(M, p, H)
  if (!.mf_sphere(M)) return(Y)
  Y + (-.jl_dot(p, G - .mf_project(M, p, G))) * X
}

# ---------------------------------------------------------------------------
# Stopping criteria

.mo_sc <- function(type, ...) structure(list(type = type, ...), class = "ManoptStoppingCriterion")

#' @export
print.ManoptStoppingCriterion <- function(x, ...) {
  if (x$type == "any") {
    cat("StopWhenAny of\n")
    for (c in x$criteria) print(c)
  } else {
    a <- x[setdiff(names(x), "type")]
    cat(x$type, "(", paste(vapply(a, function(v) format(v), ""), collapse = ", "), ")\n", sep = "")
  }
  invisible(x)
}

.mo_flatten_sc <- function(sc) {
  if (is.null(sc)) return(list())
  if (inherits(sc, "ManoptStoppingCriterion") && sc$type != "any") return(list(sc))
  if (inherits(sc, "ManoptStoppingCriterion")) return(unlist(lapply(sc$criteria, .mo_flatten_sc), recursive = FALSE))
  unlist(lapply(sc, .mo_flatten_sc), recursive = FALSE)
}

# which criteria signal convergence (indicates_convergence of Manopt)
.mo_sc_converges <- c(StopWhenGradientNormLess = TRUE, StopWhenLagrangeMultiplierLess = TRUE,
                      StopWhenFirstOrderProgress = TRUE, StopWhenPopulationCostConcentrated = TRUE,
                      StopWhenPopulationStronglyConcentrated = TRUE, StopWhenEvolutionStagnates = TRUE,
                      StopWhenBestCostInGenerationConstant = TRUE)

# A stopper for one solve: stop(st, k) evaluates the criteria of `sc` on the
# solver state `st` (an environment) at iteration k; the criterion that
# fired is kept in the stopper's environment.
.mo_stopper <- function(sc, M) {
  crits <- .mo_flatten_sc(sc)
  envs <- lapply(crits, function(c) new.env(parent = emptyenv()))
  self <- new.env(parent = emptyenv())
  self$active <- NULL
  eval1 <- function(c, e, st, k) {
    switch(c$type,
      StopAfterIteration = k >= c$n,
      StopAfter = {
        if (k <= 0) e$start <- proc.time()[["elapsed"]]
        k > 0 && proc.time()[["elapsed"]] - e$start > c$seconds
      },
      StopWhenGradientNormLess = k > 0 && .mf_norm(M, st$p, st$X) < c$tol,
      StopWhenChangeLess = {
        hit <- FALSE
        if (!is.null(e$p_old)) {
          d <- .mf_distance(M, st$p, e$p_old)
          hit <- d < c$tol && k > 0
        }
        e$p_old <- st$p
        hit
      },
      StopWhenStepsizeLess = k > 0 && st$last_stepsize < c$tol,
      StopWhenCostLess = k > 0 && st$cost(st$p) < c$tol,
      StopWhenPopulationConcentrated = k > 0 && .mo_nm_concentrated(M, st, c),
      StopWhenSwarmVelocityLess = k > 0 && .jl_norm(vapply(st$velocity, function(v) .mf_norm(M, NULL, v), 0)) < c$tol,
      StopWhenLagrangeMultiplierLess = k > 0 && -st$xi < c$tol,
      StopWhenAllLanczosVectorsUsed = k > 0 && length(st$lanczos$V) >= c$n,
      StopWhenFirstOrderProgress = .mo_first_order_progress(st, c, k),
      .mo_cma_stop(c, e, st, k))
  }
  self$check <- function(st, k) {
    if (k <= 0) {
      hit <- vapply(seq_along(crits), function(i) isTRUE(eval1(crits[[i]], envs[[i]], st, k)), TRUE)
      if (any(hit)) { self$active <- crits[[which(hit)[1L]]]; return(TRUE) }
      return(FALSE)
    }
    for (i in seq_along(crits)) {
      if (isTRUE(eval1(crits[[i]], envs[[i]], st, k))) {
        self$active <- crits[[i]]
        return(TRUE)
      }
    }
    FALSE
  }
  self
}

# ---------------------------------------------------------------------------
# Step sizes

.mo_ss <- function(type, ...) structure(list(type = type, ...), class = "ManoptStepsize")

# A step size functor for one solve: step(st, k, eta, X) for the direction
# eta at the iterate st$p with gradient X.
.mo_stepsize <- function(spec, mp, retr, vt) {
  M <- mp$M
  self <- new.env(parent = emptyenv())
  self$spec <- spec
  self$failed <- FALSE
  if (spec$type == "Armijo") {
    self$last <- spec$initial_stepsize
    swse0 <- spec$stop_when_stepsize_exceeds %||% .mf_inj(M)
    rm <- spec$retraction_method %||% retr
    self$step <- function(st, k, eta, X) {
      gX <- st$X
      ms <- .mf_inj(M)
      guess <- if (is.finite(ms)) min(self$last, ms / .mf_norm(M, st$p, gX)) else self$last
      l <- .mf_norm(M, st$p, eta)
      self$last <- .mo_backtrack(mp, st$p, guess, spec$sufficient_decrease, spec$contraction_factor, eta, X, rm,
                                 spec$stop_when_stepsize_less / l, swse0 / l, spec$stop_increasing_at_step,
                                 spec$stop_decreasing_at_step)
      self$last
    }
  } else if (spec$type == "WolfePowell") {
    self$last <- 0
    rm <- spec$retraction_method %||% retr
    vm <- spec$vector_transport_method %||% vt
    self$step <- function(st, k, eta, X) {
      r <- .mo_wolfe_powell(mp, spec, st$p, eta, X, rm, vm)
      self$failed <- r$failed
      self$last <- r$step
      r$step
    }
  } else if (spec$type == "Decreasing") {
    self$step <- function(st, k, eta, X) {
      d <- (k + spec$shift)^spec$exponent
      ds <- (spec$length - k * spec$subtrahend) * (spec$factor^k) / (if (d == 0) 1 else d)
      if (spec$kind == "absolute") {
        ns <- .mf_norm(M, st$p, X %||% st$X)
        if (ns > .Machine$double.eps) ds <- ds / ns
      }
      self$last <- ds
      ds
    }
  } else if (spec$type == "Constant") {
    self$step <- function(st, k, eta, X) {
      s <- spec$length %||% (.mf_inj(M) / 2)
      if (spec$kind == "absolute") {
        ns <- .mf_norm(M, st$p, X %||% st$X)
        if (ns > .Machine$double.eps) s <- s / ns
      }
      self$last <- s
      s
    }
  } else {
    stop("unknown step size: ", spec$type, call. = FALSE)
  }
  self
}

# linesearch_backtrack! of Manopt
.mo_backtrack <- function(mp, p, s, decrease, contract, eta, X, retr, stop_less, stop_exceeds, inc_steps, dec_steps) {
  M <- mp$M
  lf0 <- mp$f(p)
  Dlf0 <- if (is.null(X)) 0 else .mf_inner(M, p, X, eta)
  if (s > stop_exceeds) s <- min(s, stop_exceeds)
  fq <- mp$f(.mf_retract(M, p, eta, s, retr))
  i <- 0L
  while (isTRUE(fq < lf0 + decrease * s * Dlf0)) {
    if (inc_steps == 0L) break
    i <- i + 1L
    s <- min(s / contract, stop_exceeds)
    fq <- mp$f(.mf_retract(M, p, eta, s, retr))
    if (i == inc_steps) break
    if (s >= stop_exceeds) break
  }
  i <- 0L
  while (isTRUE(fq > lf0 + decrease * s * Dlf0)) {
    if (i > 0L && i == dec_steps) break
    if (i > 0L && s < stop_less) break
    i <- i + 1L
    s <- contract * s
    fq <- mp$f(.mf_retract(M, p, eta, s, retr))
  }
  s
}

# WolfePowellLinesearchStepsize of Manopt
.mo_wolfe_powell <- function(mp, spec, p, eta, X, retr, vt) {
  M <- mp$M
  c1 <- spec$sufficient_decrease; c2 <- spec$sufficient_curvature
  l <- .mf_inner(M, p, X, eta)
  gn <- .mf_norm(M, p, eta)
  mx <- spec$max_stepsize %||% .mf_inj(M)
  max_inc <- if (is.finite(mx)) min(1e9, mx / gn) else 1e9
  step <- min(spec$initial_guess, max_inc)
  s_plus <- step; s_minus <- step
  failed <- FALSE
  f0 <- mp$f(p)
  cand <- .mf_retract(M, p, eta, step, retr)
  fNew <- mp$f(cand)
  slope <- function(q) .mf_inner(M, q, mp$grad(q), .mf_transport(M, p, eta, q, vt))
  if (isTRUE(fNew > f0 + c1 * step * l)) {
    i <- 0L
    while (isTRUE(fNew > f0 + c1 * step * l && s_minus > 1e-9)) {
      s_minus <- s_minus / 2
      step <- s_minus
      fNew <- mp$f(.mf_retract(M, p, eta, step, retr))
      i <- i + 1L
      if (i == spec$stop_decreasing_at_step) break
    }
    failed <- isTRUE(fNew > f0 + c1 * step * l)
    s_plus <- min(2 * s_minus, max_inc)
  } else if (isTRUE(slope(cand) < c2 * l)) {
    i <- 0L
    while (isTRUE(fNew <= f0 + c1 * step * l && s_plus < max_inc)) {
      s_minus <- s_plus
      s_plus <- min(2 * s_plus, max_inc)
      step <- s_plus
      fNew <- mp$f(.mf_retract(M, p, eta, step, retr))
      i <- i + 1L
      if (i == spec$stop_increasing_at_step) break
    }
    if (isTRUE(s_plus == max_inc && fNew <= f0 + c1 * s_plus * l)) s_minus <- s_plus
  }
  if (s_minus < s_plus) {
    cand <- .mf_retract(M, p, eta, s_minus, retr)
    i <- 0L
    while (isTRUE(slope(cand) < c2 * l)) {
      if (i > 0L && abs(s_plus - s_minus) <= spec$stop_when_stepsize_less) { failed <- TRUE; break }
      step <- (s_minus + s_plus) / 2
      if (step == s_minus || step == s_plus) { failed <- TRUE; break }
      i <- i + 1L
      fNew <- mp$f(.mf_retract(M, p, eta, step, retr))
      if (isTRUE(fNew <= f0 + c1 * step * l)) {
        s_minus <- step
        failed <- FALSE
      } else {
        s_plus <- step
      }
      cand <- .mf_retract(M, p, eta, s_minus, retr)
    }
  }
  list(step = s_minus, failed = failed)
}

# ---------------------------------------------------------------------------
# The solvers. `mp` is the problem: the manifold M, the cost f(p), the
# Riemannian gradient grad(p) and Hessian hess(p, X). Each returns the
# result point and the solver state; st$stopper$active is the criterion
# that stopped it.

.mo_state <- function(mp, p, sc) {
  st <- new.env(parent = emptyenv())
  # a positive definite matrix is held by columns
  st$p <- as.numeric(p)
  st$M <- mp$M
  st$cost <- mp$f
  st$stopper <- .mo_stopper(sc, mp$M)
  st$iterations <- 0L
  st
}

.mo_run <- function(st, step) {
  k <- 0L
  while (!st$stopper$check(st, k)) {
    k <- k + 1L
    step(k)
  }
  st$iterations <- k
  st
}

.mo_gradient_descent <- function(mp, p, stopping_criterion = NULL, stepsize = NULL, retraction_method = "exp", ...) {
  M <- mp$M
  sc <- stopping_criterion %||% list(.mo_sc("StopAfterIteration", n = 200L),
                                     .mo_sc("StopWhenGradientNormLess", tol = 1e-8))
  st <- .mo_state(mp, p, sc)
  ss <- .mo_stepsize(stepsize %||% .manopt_armijo(), mp, retraction_method, "parallel")
  st$X <- mp$grad(st$p)
  st$last_stepsize <- 0
  .mo_run(st, function(k) {
    st$X <- mp$grad(st$p)
    s <- ss$step(st, k, -st$X, st$X)
    st$last_stepsize <- s
    st$p <- .mf_retract(M, st$p, st$X, -s, retraction_method)
  })
}

.mo_conjugate_gradient <- function(mp, p, stopping_criterion = NULL, stepsize = NULL, retraction_method = "exp",
                                   vector_transport_method = "parallel", ...) {
  M <- mp$M
  sc <- stopping_criterion %||% list(.mo_sc("StopAfterIteration", n = 500L),
                                     .mo_sc("StopWhenGradientNormLess", tol = 1e-8))
  st <- .mo_state(mp, p, sc)
  ss <- .mo_stepsize(stepsize %||% .manopt_armijo(), mp, retraction_method, vector_transport_method)
  st$X <- mp$grad(st$p)
  st$delta <- -st$X
  store <- list(p = st$p, X = st$X)
  st$beta <- 0
  st$last_stepsize <- 0
  .mo_run(st, function(k) {
    p_old <- st$p
    s <- ss$step(st, k, st$delta, st$X)
    st$last_stepsize <- s
    st$p <- .mf_retract(M, p_old, st$delta, s, retraction_method)
    st$X <- mp$grad(st$p)
    # ConjugateDescentCoefficient
    st$beta <- .mf_inner(M, st$p, st$X, st$X) / .mf_inner(M, store$p, store$X, -st$delta)
    store <<- list(p = st$p, X = st$X)
    d <- .mf_transport(M, p_old, st$delta, st$p, vector_transport_method)
    d <- d * st$beta
    st$delta <- d - st$X
    # RestartOnNonDescent
    if (.mf_inner(M, st$p, st$X, st$delta) >= 0) {
      st$delta <- -st$X
      st$beta <- 0
    }
  })
}

# limited memory inverse BFGS: the two loop recursion of Manopt
.mo_lbfgs_direction <- function(M, mem, p, X) {
  r <- X
  m <- length(mem$S)
  if (m == 0L) return(-r)
  for (i in m:1L) {
    mem$xi[i] <- .mf_inner(M, p, mem$S[[i]], r) * mem$rho[i]
    r <- r - mem$xi[i] * mem$Y[[i]]
  }
  safe <- which(abs(mem$rho) > 0)
  if (!length(safe)) return(-r)
  ls <- max(safe)
  ny <- .mf_norm(M, p, mem$Y[[ls]])
  r <- r * (mem$initial_scale / (mem$rho[ls] * (ny * ny)))
  for (i in seq_along(mem$rho)) {
    if (abs(mem$rho[i]) > 0) {
      coeff <- mem$xi[i] - mem$rho[i] * .mf_inner(M, p, mem$Y[[i]], r)
      r <- r + coeff * mem$S[[i]]
    }
  }
  -r
}

.mo_lbfgs_update <- function(M, mem, p_old, p, sk, yk, vt) {
  cap <- length(mem$rho)
  if (cap == 0L) return(invisible(mem))
  m <- length(mem$S)
  start <- if (m == cap) 2L else 1L
  rho_of <- function(i) {
    v <- .mf_inner(M, p, mem$S[[i]], mem$Y[[i]])
    if (v == 0) 0 else 1 / v
  }
  if (m >= start) for (i in start:m) {
    mem$S[[i]] <- .mf_transport(M, p_old, mem$S[[i]], p, vt)
    mem$Y[[i]] <- .mf_transport(M, p_old, mem$Y[[i]], p, vt)
    mem$rho[i] <- rho_of(i)
  }
  if (m == cap) {
    mem$S <- c(mem$S[-1L], list(sk))
    mem$Y <- c(mem$Y[-1L], list(yk))
    mem$rho <- c(mem$rho[-1L], mem$rho[1L])
  } else {
    mem$S <- c(mem$S, list(sk))
    mem$Y <- c(mem$Y, list(yk))
  }
  mem$rho[length(mem$S)] <- rho_of(length(mem$S))
  invisible(mem)
}

.mo_quasi_newton <- function(mp, p, stopping_criterion = NULL, stepsize = NULL, retraction_method = "exp",
                             vector_transport_method = "parallel", memory_size = NULL, initial_scale = 1, ...) {
  M <- mp$M
  msize <- as.integer(memory_size %||% min(.mf_dim(M), 20L))
  sc <- stopping_criterion %||% list(.mo_sc("StopAfterIteration", n = max(1000L, msize)),
                                     .mo_sc("StopWhenGradientNormLess", tol = 1e-6))
  st <- .mo_state(mp, p, sc)
  ss <- .mo_stepsize(stepsize %||% .manopt_wolfe_powell(stop_when_stepsize_less = 1e-10), mp, retraction_method,
                     vector_transport_method)
  mem <- new.env(parent = emptyenv())
  reset <- function() {
    mem$S <- list(); mem$Y <- list()
    mem$rho <- numeric(msize); mem$xi <- numeric(msize)
    mem$initial_scale <- initial_scale
  }
  reset()
  st$X <- mp$grad(st$p)
  st$last_stepsize <- 0
  .mo_run(st, function(k) {
    eta <- .mo_lbfgs_direction(M, mem, st$p, st$X)
    if (.mf_inner(M, st$p, eta, st$X) >= 0) {
      reset()
      eta <- .mo_lbfgs_direction(M, mem, st$p, st$X)
    }
    alpha <- ss$step(st, k, eta, st$X)
    st$last_stepsize <- alpha
    p_old <- st$p
    st$p <- .mf_retract(M, p_old, eta, alpha, retraction_method)
    eta <- eta * alpha
    vanish <- alpha == 0 || .mf_norm(M, p_old, eta) == 0
    sk <- .mf_transport(M, p_old, eta, st$p, vector_transport_method)
    beta <- if (vanish) 1 else .mf_norm(M, p_old, eta) / .mf_norm(M, st$p, sk)
    X_old <- .mf_transport(M, p_old, st$X, st$p, vector_transport_method)
    st$X <- mp$grad(st$p)
    yk <- st$X / beta - X_old
    if (!vanish) .mo_lbfgs_update(M, mem, p_old, st$p, sk, yk, vector_transport_method)
  })
}

# the truncated conjugate gradient subsolver of the trust region method
.mo_tcg <- function(mp, p, X, radius, kappa = 0.1, theta = 1) {
  M <- mp$M
  Y <- numeric(length(p)); HY <- Y
  r <- X
  z <- r
  d <- -z
  dHd <- 0; YPd <- 0
  dPd <- .mf_inner(M, p, r, z)
  YPY <- 0
  model <- 0
  zr <- .mf_inner(M, p, z, r)
  r0 <- .mf_norm(M, p, r)
  r2 <- radius * radius
  model_seen <- Inf
  maxit <- .mf_dim(M)
  k <- 0L
  repeat {
    # StopAfterIteration | residual reduced | trust region exceeded |
    # negative curvature | model increased
    if (k >= maxit) break
    if (k > 0L && .mf_norm(M, p, r) <= r0 * min(kappa, r0^theta)) break
    if (YPY >= r2) break
    if (k > 0L && dHd <= 0) break
    if (k > 0L && model > model_seen) break
    model_seen <- model
    k <- k + 1L
    if (zr == 0) next
    Hd <- mp$hess(p, d)
    dHd <- .mf_inner(M, p, d, Hd)
    a <- zr / dHd
    YPY_new <- YPY + 2 * a * YPd + a^2 * dPd
    if (dHd <= 0 || YPY_new >= r2) {
      tau <- (-YPd + sqrt(YPd^2 + dPd * (r2 - YPY))) / dPd
      Y <- Y + tau * d
      HY <- HY + tau * Hd
      YPY <- r2
      next
    }
    YPY <- YPY_new
    nY <- Y + a * d
    nHY <- HY + a * Hd
    nmodel <- .mf_inner(M, p, nY, X) + 0.5 * .mf_inner(M, p, nY, nHY)
    if (nmodel >= model) {
      model <- nmodel
      next
    }
    Y <- nY; model <- nmodel; HY <- nHY
    r <- r + a * Hd
    z <- r
    zrn <- .mf_inner(M, p, z, r)
    b <- zrn / zr
    zr <- zrn
    d <- -z + b * d
    YPd <- b * (a * dPd + YPd)
    dPd <- zr + b^2 * dPd
  }
  list(Y = Y, HY = HY, exceeded = YPY >= r2)
}

.mo_trust_regions <- function(mp, p, stopping_criterion = NULL, retraction_method = "exp",
                              max_trust_region_radius = NULL, trust_region_radius = NULL, acceptance_rate = 0.1,
                              rho_regularization = 1e3, kappa = 0.1, theta = 1, reduction_threshold = 0.1,
                              reduction_factor = 0.25, augmentation_threshold = 0.75, augmentation_factor = 2, ...) {
  M <- mp$M
  dmax <- max_trust_region_radius %||% sqrt(.mf_dim(M))
  st <- .mo_state(mp, p, stopping_criterion %||% list(.mo_sc("StopAfterIteration", n = 1000L),
                                                     .mo_sc("StopWhenGradientNormLess", tol = 1e-6)))
  st$radius <- trust_region_radius %||% (dmax / 8)
  st$X <- mp$grad(st$p)
  st$last_stepsize <- 0
  .mo_run(st, function(k) {
    sub <- .mo_tcg(mp, st$p, mp$grad(st$p), st$radius, kappa, theta)
    Y <- sub$Y; HY <- sub$HY
    f <- mp$f(st$p)
    prop <- .mf_retract(M, st$p, Y, 1, retraction_method)
    reg <- max(1, abs(f)) * .Machine$double.eps * rho_regularization
    num <- f - mp$f(prop)
    den <- -.mf_inner(M, st$p, Y, st$X) - 0.5 * .mf_inner(M, st$p, Y, HY)
    num <- num + reg
    den <- den + reg
    rho <- num / den
    decreased <- den >= 0
    if (rho < reduction_threshold || !decreased || is.nan(rho)) {
      st$radius <- st$radius * reduction_factor
    } else if (rho > augmentation_threshold && sub$exceeded) {
      st$radius <- min(augmentation_factor * st$radius, dmax)
    }
    if (decreased && rho > acceptance_rate) {
      st$p <- prop
      st$X <- mp$grad(st$p)
    }
  })
}

# the Lanczos subsolver of the adaptive regularization with cubics
.mo_lanczos <- function(mp, p, X, sigma, ls, maxit, theta) {
  M <- mp$M
  nX <- .mf_norm(M, p, X)
  if (nX == 0) return(numeric(length(p)))
  Tm <- matrix(0, maxit, maxit)
  coef <- numeric(maxit)
  S <- numeric(length(p))
  Hp <- NULL; res <- NULL
  k <- 0L
  repeat {
    if (k >= maxit) break
    if (k > 0L) {
      # StopWhenFirstOrderProgress
      y <- coef[seq_len(k - 1L)]
      Ty <- if (k > 1L) Tm[seq_len(k), seq_len(k - 1L), drop = FALSE] %*% y else numeric(k)
      ny <- .jl_norm(y)
      mg <- .jl_norm(nX * c(1, numeric(k - 1L)) + as.numeric(Ty) + sigma * ny * c(y, 0))
      if (k > 1L && mg <= theta * ny^2) break
    }
    k <- k + 1L
    if (k == 1L) {
      v <- X / nX
      if (length(ls$V) == 0L) ls$V <- list(v) else ls$V[[1L]] <- v
      Hp <- mp$hess(p, v)
      a <- .mf_inner(M, p, v, Hp)
      Tm[1L, 1L] <- a
      res <- Hp - a * v
      coef[1L] <- (a - sqrt(a^2 + 4 * sigma * nX)) / (2 * sigma)
    } else {
      b <- .mf_norm(M, p, res)
      if (b > 1e-12) {
        v <- res / b
      } else {
        res <- .mf_rand_vector(M, p)
        for (j in seq_len(k - 1L)) res <- res - .mf_inner(M, p, ls$V[[j]], res) * ls$V[[j]]
        v <- res / .mf_norm(M, p, res)
      }
      if (length(ls$V) < k) ls$V[[k]] <- v else ls$V[[k]] <- v
      Hp <- mp$hess(p, v)
      res <- Hp - b * ls$V[[k - 1L]]
      a <- .mf_inner(M, p, res, v)
      res <- res - a * v
      Tm[k, k] <- a; Tm[k - 1L, k] <- b; Tm[k, k - 1L] <- b
      coef[seq_len(k)] <- .mo_min_cubic_newton(Tm[seq_len(k), seq_len(k), drop = FALSE], nX, sigma)
    }
    S <- Reduce(`+`, lapply(seq_len(k), function(j) ls$V[[j]] * coef[j]))
  }
  S
}

.mo_min_cubic_newton <- function(Tk, nX, sigma) {
  k <- nrow(Tk)
  g <- c(nX, numeric(k - 1L))
  lam <- max(svd(Tk, nu = 0, nv = 0)$d) + 2
  lmin <- min(eigen(Tk, symmetric = TRUE, only.values = TRUE)$values)
  lower <- max(0, -lmin)
  y <- numeric(k)
  for (j in seq_len(200L)) {
    Tl <- Tk + diag(lam, k)
    y <- -solve(Tl, g, tol = 0)
    yn <- sqrt(sum(y * y))
    phi <- 1 / yn - sigma / lam
    if (abs(phi) < 1e-16 * yn) break
    psi <- yn^2
    dy <- -solve(Tl, y, tol = 0)
    psi_p <- 2 * sum(y * dy)
    p0 <- 2 * sigma * psi^1.5
    p1 <- -2 * psi - lam * psi_p
    p2 <- psi_p
    disc <- sqrt(p1^2 - 4 * p2 * p0)
    dl <- max((-p1 + disc) / (2 * p2), (-p1 - disc) / (2 * p2)) - lam
    if (!is.finite(dl)) break
    if (lam + dl <= lower) dl <- -0.5 * (lam - lower)
    if (abs(dl) <= .jl_eps(lam)) break
    lam <- lam + dl
  }
  y
}

.mo_first_order_progress <- function(st, c, k) FALSE

.mo_arc <- function(mp, p, stopping_criterion = NULL, retraction_method = "exp", sigma = NULL, sigma_min = 1e-10,
                    eta1 = 0.1, eta2 = 0.9, gamma1 = 0.1, gamma2 = 2, theta = 0.5, rho_regularization = 1e3,
                    maxIterLanczos = NULL, ...) {
  M <- mp$M
  maxl <- as.integer(maxIterLanczos %||% min(300L, .mf_dim(M)))
  sc <- stopping_criterion %||% list(.mo_sc("StopAfterIteration", n = 40L),
                                     .mo_sc("StopWhenGradientNormLess", tol = 1e-9),
                                     .mo_sc("StopWhenAllLanczosVectorsUsed", n = maxl - 1L))
  st <- .mo_state(mp, p, sc)
  st$sigma <- sigma %||% (100 / sqrt(.mf_dim(M)))
  st$lanczos <- new.env(parent = emptyenv())
  st$lanczos$V <- list()
  st$X <- mp$grad(st$p)
  st$last_stepsize <- 0
  .mo_run(st, function(k) {
    s <- .mo_lanczos(mp, st$p, st$X, st$sigma, st$lanczos, maxl, theta)
    q <- .mf_retract(M, st$p, s, 1, retraction_method)
    cost <- mp$f(st$p)
    num <- cost - mp$f(q)
    H <- mp$hess(st$p, s)
    den <- -.mf_inner(M, st$p, s, st$X + 0.5 * H)
    reg <- rho_regularization * .Machine$double.eps * max(abs(cost), 1)
    rho <- (num + reg) / (den + reg)
    if (rho >= eta1) {
      st$p <- q
      st$X <- mp$grad(st$p)
    }
    if (rho >= eta2) st$sigma <- max(sigma_min, gamma1 * st$sigma) else if (rho < eta1) st$sigma <- gamma2 * st$sigma
  })
}

# the gradient of the Frank Wolfe oracle q -> <X, log_p q>: the adjoint
# differential of the logarithm in its argument, the adjoint Jacobi field of
# ManifoldDiff along the geodesic from q to p in the basis diagonalizing the
# curvature
.mo_adjoint_dlog <- function(M, p, q, X) {
  if (.mf_spd(M)) {
    n <- M$n
    d <- .spd_distance(n, q, p)
    if (.jl_norm(q - p) <= sqrt(.Machine$double.eps) * max(.jl_norm(q), .jl_norm(p))) return(X)
    v <- .spd_log(n, q, p)
    tX <- .spd_transport(n, p, X, q)
    B <- .spd_basis(n, q, v / d)
    Y <- numeric(n * n)
    for (k in seq_along(B$vectors)) {
      kap <- B$kappa[k]
      beta <- if (d == 0 || kap == 0) 1 else if (kap < 0) sqrt(-kap) * d / sinh(sqrt(-kap) * d)
              else sqrt(kap) * d / sin(sqrt(kap) * d)
      Y <- Y + beta * .spd_inner(n, q, B$vectors[[k]], tX) * B$vectors[[k]]
    }
    return(Y)
  }
  if (!.mf_sphere(M)) return(X)
  d <- .mf_distance(M, q, p)
  if (.jl_norm(q - p) <= sqrt(.Machine$double.eps) * max(.jl_norm(q), .jl_norm(p))) return(X)
  v <- .mf_log(M, q, p)
  x <- .mf_project_point(M, .mf_exp(M, q, v, 1))
  v <- v / d
  v <- v / .mf_norm(M, q, v)
  tX <- .mf_transport(M, x, X, q)
  along <- v * .jl_dot(v, tX)
  along + (d / sin(d)) * (tX - along)
}

.mo_frank_wolfe <- function(mp, p, stopping_criterion = NULL, stepsize = NULL, retraction_method = "exp",
                            sub_stopping_criterion = NULL, ...) {
  M <- mp$M
  sc <- stopping_criterion %||% list(.mo_sc("StopAfterIteration", n = 200L),
                                     .mo_sc("StopWhenGradientNormLess", tol = 1e-8),
                                     .mo_sc("StopWhenChangeLess", tol = 1e-8))
  st <- .mo_state(mp, p, sc)
  ss <- .mo_stepsize(stepsize %||% .manopt_decreasing(length = 2, shift = 2), mp, retraction_method, "parallel")
  st$X <- mp$grad(st$p)
  st$last_stepsize <- 0
  # the oracle: gradient descent on q -> <X, log_p q>, started from the
  # previous oracle solution with a step size that keeps its state
  sub_mp <- list(M = M, f = function(q) .mf_inner(M, st$p, st$X, .mf_log(M, st$p, q)),
                 grad = function(q) .mo_adjoint_dlog(M, st$p, q, st$X))
  sub_sc <- sub_stopping_criterion %||% list(.mo_sc("StopAfterIteration", n = 300L),
                                             .mo_sc("StopWhenStepsizeLess", tol = 1e-8))
  sub_ss <- .mo_stepsize(.manopt_armijo(), sub_mp, retraction_method, "parallel")
  q <- st$p
  .mo_run(st, function(k) {
    st$X <- mp$grad(st$p)
    sst <- .mo_state(sub_mp, q, sub_sc)
    sst$X <- sub_mp$grad(sst$p)
    sst$last_stepsize <- sub_ss$last
    .mo_run(sst, function(j) {
      sst$X <- sub_mp$grad(sst$p)
      s <- sub_ss$step(sst, j, -sst$X, sst$X)
      sst$last_stepsize <- s
      sst$p <- .mf_retract(M, sst$p, sst$X, -s, retraction_method)
    })
    q <<- sst$p
    s <- ss$step(st, k, NULL, st$X)
    st$last_stepsize <- s
    st$p <- .mf_retract(M, st$p, s * .mf_invretract(M, st$p, q), 1, retraction_method)
  })
}
