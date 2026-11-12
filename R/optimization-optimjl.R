# The optimizers of Optim.jl (v2.3.2) with the line searches of
# LineSearches.jl (v7.8.2) and the caching objectives of NLSolversBase.jl
# (v8.0.1), reached through OptimizationOptimJL (v0.4.23): first order
# methods (src/multivariate/solvers/first_order), Newton methods
# (second_order), the box constrained Fminbox (constrained/fminbox.jl), the
# interior point IPNewton (constrained/ipnewton), the zeroth order methods
# (zeroth_order) and SAMIN (constrained/samin.jl), with the driver loop of
# multivariate/optimize/optimize.jl.

# ---------------------------------------------------------------------------
# Floating point helpers

.ojl_ulp <- function(x) {
  a <- abs(x)
  if (!is.finite(a)) return(NaN)
  if (a < 0x1.0000000000000p-1022) return(4.9406564584124654e-324)
  e <- floor(log2(a))
  if (2^e > a) e <- e - 1 else if (2^(e + 1) <= a) e <- e + 1
  2^(e - 52)
}

.ojl_nextfloat <- function(x) {
  if (is.na(x) || x == Inf) return(x)
  if (x == -Inf) return(-.Machine$double.xmax)
  if (x == 0) return(4.9406564584124654e-324)
  if (x > 0) return(x + .ojl_ulp(x))
  a <- -x
  u <- .ojl_ulp(a)
  if (a >= 0x1.0000000000000p-1022 && a == 2^(log2(u) + 52)) u <- u / 2
  x + u
}

# Julia's `x != y` for arrays (true when any entry differs or is NaN)
.ojl_ne <- function(a, b) !isTRUE(all(a == b))

.ojl_nanmax <- function(a, b) if (is.nan(a)) b else if (is.nan(b)) a else max(a, b)

.ojl_clamp <- function(x, lo, hi) if (isTRUE(x > hi)) hi else if (isTRUE(x < lo)) lo else x

.ojl_mod1 <- function(i, m) ((i - 1L) %% m) + 1L

# muladd(a, b, c) of Julia, which is fused (one rounding) on hardware with a
# fused multiply add. The emulation of Boldo and Melquiond (2008): a*b is
# split exactly (Dekker), c + a*b is summed exactly (TwoSum), the small parts
# are added with rounding to odd, and the final sum rounds to nearest; this is
# the correctly rounded a*b + c, ties included
.ojl_fma <- function(a, b, c) {
  n <- max(length(a), length(b), length(c))
  a <- rep_len(a, n); b <- rep_len(b, n); c <- rep_len(c, n)
  p <- a * b
  sa <- 134217729 * a
  ah <- sa - (sa - a)
  al <- a - ah
  sb <- 134217729 * b
  bh <- sb - (sb - b)
  bl <- b - bh
  e <- al * bl - (((p - ah * bh) - al * bh) - ah * bl)
  th <- c + p
  bb <- th - c
  tl <- (c - (th - bb)) + (p - bb)
  v <- tl + e
  bb <- v - tl
  w <- (tl - (v - bb)) + (e - bb)
  k <- which(w != 0)
  if (length(k)) {
    vk <- v[k]
    av <- abs(vk)
    ex <- floor(log2(av))
    lo <- 2^ex > av
    ex[lo] <- ex[lo] - 1
    hi <- 2 * 2^ex <= av
    ex[hi] <- ex[hi] + 1
    u <- 2^(ex - 52)
    even <- (av / u) %% 2 == 0
    up <- sign(w[k]) == sign(vk)
    step <- ifelse(up, u, ifelse(av == 2^ex, u / 2, u))
    v[k] <- ifelse(even, vk + ifelse(up, 1, -1) * sign(vk) * step, vk)
  }
  out <- th + v
  bad <- !is.finite(out)
  out[bad] <- (p + c)[bad]
  out
}

# dot(x, y) and A * x of Julia's BLAS (OpenBLAS): for up to 16 entries the
# dot product is the chain of fused multiply adds x1 y1, fma(x2, y2, .), ...
# and every row of a matrix vector product is such a chain over the columns
.ojl_dot <- function(x, y) {
  n <- length(x)
  if (!n) return(0)
  if (n > 16L) return(sum(x * y))
  r <- x[1L] * y[1L]
  for (i in seq_len(n - 1L) + 1L) r <- .ojl_fma(x[i], y[i], r)
  r
}

# A \ b of Julia for a small dense square matrix: the triangular shortcuts of
# LinearAlgebra, else the LU factorization of OpenBLAS (the left looking getf2
# with fused dot products, reciprocal pivot scaling and the axpy based
# triangular solves of getrs), so the solution matches Julia bit for bit
.ojl_backslash <- function(A, b) {
  n <- nrow(A)
  low <- all(A[upper.tri(A)] == 0)
  upp <- all(A[lower.tri(A)] == 0)
  if (low && upp) return(b / diag(A))
  ipiv <- seq_len(n)
  if (!low && !upp) {
    for (j in seq_len(n)) {
      col <- A[, j]
      for (i in seq_len(j - 1L)) {
        jp <- ipiv[i]
        if (jp != i) {
          t <- col[i]; col[i] <- col[jp]; col[jp] <- t
        }
      }
      if (j > 2L) for (i in 2:(j - 1L)) col[i] <- col[i] - .ojl_dot(A[i, seq_len(i - 1L)], col[seq_len(i - 1L)])
      rows <- j:n
      for (k in seq_len(j - 1L)) col[rows] <- .ojl_fma(-col[k], A[rows, k], col[rows])
      jp <- j - 1L + which.max(abs(col[rows]))
      ipiv[j] <- jp
      A[, j] <- col
      piv <- A[jp, j]
      if (piv == 0) stop("SingularException(", j, ")", call. = FALSE)
      if (jp != j) {
        cs <- seq_len(j)
        t <- A[j, cs]; A[j, cs] <- A[jp, cs]; A[jp, cs] <- t
      }
      if (j < n) A[(j + 1L):n, j] <- A[(j + 1L):n, j] * (1 / piv)
    }
    for (i in seq_len(n)) {
      jp <- ipiv[i]
      if (jp != i) {
        t <- b[i]; b[i] <- b[jp]; b[jp] <- t
      }
    }
  }
  if (!upp) for (i in seq_len(n)) {
    if (low) b[i] <- b[i] / A[i, i]
    if (i < n) b[(i + 1L):n] <- .ojl_fma(-b[i], A[(i + 1L):n, i], b[(i + 1L):n])
  }
  if (!low) for (i in n:1) {
    b[i] <- b[i] / A[i, i]
    if (i > 1L) b[seq_len(i - 1L)] <- .ojl_fma(-b[i], A[seq_len(i - 1L), i], b[seq_len(i - 1L)])
  }
  b
}

.ojl_gemv <- function(A, x) {
  n <- ncol(A)
  if (!n) return(numeric(nrow(A)))
  r <- A[, 1L] * x[1L]
  for (j in seq_len(n - 1L) + 1L) r <- .ojl_fma(A[, j], x[j], r)
  r
}

.ojl_opt_fd_hessian <- function(g, x) {
  n <- length(x)
  H <- matrix(0, n, n)
  for (j in seq_len(n)) {
    h <- .Machine$double.eps^(1 / 3) * max(1, abs(x[j]))
    xp <- x; xm <- x
    xp[j] <- x[j] + h; xm[j] <- x[j] - h
    H[, j] <- (g(xp) - g(xm)) / (2 * h)
  }
  (H + t(H)) / 2
}

# ---------------------------------------------------------------------------
# NLSolversBase objectives: the value, gradient and Hessian of the last
# evaluation points are cached and counted

.ojl_objective <- function(f, g, h = NULL, n, hv = NULL) {
  d <- new.env(parent = emptyenv())
  d$F <- NaN; d$DF <- rep(NaN, n); d$H <- matrix(NaN, n, n); d$JVP <- NaN
  d$x_f <- d$x_df <- d$x_jvp <- d$v_jvp <- d$x_h <- rep(NaN, n)
  d$f_calls <- 0L; d$df_calls <- 0L; d$h_calls <- 0L
  vbb <- function(x) {
    d$F <- f(x); d$f_calls <- d$f_calls + 1L; d$x_f <- x
    d$F
  }
  gbb <- function(x) {
    d$DF <- g(x); d$df_calls <- d$df_calls + 1L; d$x_df <- x
    d$DF
  }
  vgbb <- function(x) {
    d$DF <- g(x); d$F <- f(x)
    d$f_calls <- d$f_calls + 1L; d$df_calls <- d$df_calls + 1L
    d$x_f <- x; d$x_df <- x
  }
  hbb <- function(x) {
    d$H <- h(x); d$h_calls <- d$h_calls + 1L; d$x_h <- x
    d$H
  }
  d$value <- function(x) {
    d$f_calls <- d$f_calls + 1L
    f(x)
  }
  d$value_bang <- function(x) {
    if (.ojl_ne(x, d$x_f)) vbb(x)
    d$F
  }
  d$gradient_bang <- function(x) {
    if (.ojl_ne(x, d$x_df)) gbb(x)
    d$DF
  }
  d$value_gradient <- function(x) {
    if (.ojl_ne(x, d$x_f)) {
      if (.ojl_ne(x, d$x_df)) vgbb(x) else vbb(x)
    } else if (.ojl_ne(x, d$x_df)) gbb(x)
    list(d$F, d$DF)
  }
  d$value_jvp <- function(x, v) {
    if (.ojl_ne(x, d$x_f)) {
      if (.ojl_ne(x, d$x_jvp) || .ojl_ne(v, d$v_jvp)) {
        d$value_gradient(x)
        d$JVP <- .ojl_dot(d$DF, v); d$x_jvp <- x; d$v_jvp <- v
      } else vbb(x)
    } else if (.ojl_ne(x, d$x_jvp) || .ojl_ne(v, d$v_jvp)) {
      d$gradient_bang(x)
      d$JVP <- .ojl_dot(d$DF, v); d$x_jvp <- x; d$v_jvp <- v
    }
    list(d$F, d$JVP)
  }
  d$hessian_bang <- function(x) {
    if (.ojl_ne(x, d$x_h)) hbb(x)
    d$H
  }
  d$value_gradient_hessian <- function(x) {
    if (.ojl_ne(x, d$x_f)) {
      if (.ojl_ne(x, d$x_h)) {
        vgbb(x); hbb(x)
      } else if (.ojl_ne(x, d$x_df)) vgbb(x) else vbb(x)
    } else if (.ojl_ne(x, d$x_df)) {
      gbb(x)
      if (.ojl_ne(x, d$x_h)) hbb(x)
    } else if (.ojl_ne(x, d$x_h)) hbb(x)
    list(d$F, d$DF, d$H)
  }
  d$gradient_hessian <- function(x) {
    if (.ojl_ne(x, d$x_df)) gbb(x)
    if (.ojl_ne(x, d$x_h)) hbb(x)
    list(d$DF, d$H)
  }
  # hvp! of a TwiceDifferentiable built with a Hessian vector product (the hv
  # of OptimizationOptimJL for KrylovTrustRegion), cached on x and v
  d$HVP <- rep(NaN, n); d$x_hvp <- d$v_hvp <- rep(NaN, n); d$hvp_calls <- 0L
  d$hvp <- function(x, v) {
    if (.ojl_ne(x, d$x_hvp) || .ojl_ne(v, d$v_hvp)) {
      d$HVP <- if (is.function(hv)) hv(x, v) else .ojl_gemv(h(x), v)
      d$hvp_calls <- d$hvp_calls + 1L
      d$x_hvp <- x; d$v_hvp <- v
    }
    d$HVP
  }
  d$obj <- d
  d
}

# ---------------------------------------------------------------------------
# Manifolds (Manifolds.jl): Flat (the default, where every operation is the
# identity), Sphere, Stiefel, PowerManifold and ProductManifold. retract
# maps a point back to the manifold and project the gradient onto the
# tangent space at a point

.ojl_manifold <- function(type, ...) structure(list(type = type, ...), class = c(type, "Manifold"))

.ojl_Flat <- function() .ojl_manifold("Flat")

.ojl_Sphere <- function() .ojl_manifold("Sphere")

.ojl_Stiefel <- function(retraction = "SVD") {
  retraction <- as.character(retraction)
  if (!(retraction %in% c("SVD", "CholQR")))
    stop("the retraction of Stiefel must be \"SVD\" or \"CholQR\"", call. = FALSE)
  .ojl_manifold("Stiefel", retraction = retraction)
}

.ojl_inner_manifold <- function(M, dims) {
  if (!inherits(M, "Manifold"))
    stop("the inner manifolds must be manifolds: Optim$Flat(), Optim$Sphere() or Optim$Stiefel()", call. = FALSE)
  M$dims <- as.integer(dims)
  M
}

# the product of prod(outer_dims) copies of a manifold, stored one after the
# other as blocks of prod(inner_dims) entries
.ojl_PowerManifold <- function(inner_manifold, inner_dims, outer_dims) {
  .ojl_manifold("PowerManifold", inner_manifold = .ojl_inner_manifold(inner_manifold, inner_dims),
                inner_dims = as.integer(inner_dims), outer_dims = as.integer(outer_dims))
}

# the product of two manifolds, x1 stored before x2
.ojl_ProductManifold <- function(m1, m2, dims1, dims2) {
  .ojl_manifold("ProductManifold", m1 = .ojl_inner_manifold(m1, dims1), m2 = .ojl_inner_manifold(m2, dims2),
                dims1 = as.integer(dims1), dims2 = as.integer(dims2))
}

# the index blocks of get_inner: one per inner manifold
.ojl_inner_blocks <- function(M, x) {
  if (M$type == "PowerManifold") {
    ni <- prod(M$inner_dims)
    no <- prod(M$outer_dims)
    if (length(x) < ni * no)
      stop("BoundsError: the point has ", length(x), " entries, the power manifold needs ", ni * no, call. = FALSE)
    lapply(seq_len(no), function(i) ((i - 1) * ni + 1):(i * ni))
  } else {
    n1 <- prod(M$dims1)
    n2 <- prod(M$dims2)
    if (length(x) != n1 + n2) stop("AssertionError: length(x) == N1 + N2", call. = FALSE)
    list(seq_len(n1), n1 + seq_len(n2))
  }
}

.ojl_inner_list <- function(M) {
  if (M$type == "PowerManifold") list(M$inner_manifold) else list(M$m1, M$m2)
}

# a manifold argument: NULL or Flat (Euclidean space), or a manifold
.ojl_as_manifold <- function(manifold) {
  if (is.null(manifold) || identical(manifold, "Flat")) return(NULL)
  if (!inherits(manifold, "Manifold"))
    stop("`manifold` must be NULL or a manifold: Optim$Flat(), Optim$Sphere(), Optim$Stiefel(), ",
         "Optim$PowerManifold() or Optim$ProductManifold()", call. = FALSE)
  if (manifold$type == "Flat") NULL else manifold
}

# norm(x) of LinearAlgebra (generic_norm2, the sum of squares in order)
.ojl_norm2 <- function(x) {
  if (!length(x)) return(0)
  maxabs <- max(abs(x))
  if (is.na(maxabs) || maxabs == 0 || is.infinite(maxabs)) return(maxabs)
  if (is.finite(length(x) * maxabs * maxabs) && maxabs * maxabs != 0) {
    s <- x[1L] * x[1L]
    for (v in x[-1L]) s <- s + v * v
    return(sqrt(s))
  }
  s <- 0
  for (v in x) {
    t <- v / maxabs
    s <- s + t * t
  }
  maxabs * sqrt(s)
}

.ojl_stiefel_dims <- function(M, x) M$dims %||% c(length(x), 1L)

.ojl_retract <- function(M, x) {
  if (is.null(M)) return(x)
  switch(M$type,
    Flat = x,
    Sphere = x / .ojl_norm2(x),
    Stiefel = {
      X <- matrix(x, .ojl_stiefel_dims(M, x)[1L])
      Y <- if (M$retraction == "SVD") {
        s <- svd(X)
        s$u %*% t(s$v)
      } else {
        U <- chol(crossprod(X))
        t(backsolve(U, t(X), transpose = TRUE))
      }
      as.numeric(Y)
    },
    PowerManifold = , ProductManifold = {
      ms <- .ojl_inner_list(M)
      bl <- .ojl_inner_blocks(M, x)
      for (i in seq_along(bl)) x[bl[[i]]] <- .ojl_retract(ms[[min(i, length(ms))]], x[bl[[i]]])
      x
    })
}

.ojl_project <- function(M, g, x) {
  if (is.null(M)) return(g)
  switch(M$type,
    Flat = g,
    Sphere = g - .ojl_dot(x, g) * x,
    Stiefel = {
      nr <- .ojl_stiefel_dims(M, x)[1L]
      X <- matrix(x, nr); G <- matrix(g, nr)
      XG <- crossprod(X, G)
      as.numeric(G - X %*% ((XG + t(XG)) / 2))
    },
    PowerManifold = , ProductManifold = {
      ms <- .ojl_inner_list(M)
      bl <- .ojl_inner_blocks(M, x)
      for (i in seq_along(bl)) g[bl[[i]]] <- .ojl_project(ms[[min(i, length(ms))]], g[bl[[i]]], x[bl[[i]]])
      g
    })
}

# ManifoldObjective: the objective seen by the line searches, evaluated at
# the retracted point with the projected gradient. As in Julia, the projection
# overwrites the cached gradient of the objective
.ojl_manifold_objective <- function(M, d) {
  if (is.null(M)) return(d)
  w <- new.env(parent = emptyenv())
  w$obj <- d$obj
  w$value_bang <- function(x) d$value_bang(.ojl_retract(M, x))
  w$value_jvp <- function(x, v) {
    xin <- .ojl_retract(M, x)
    r <- d$value_gradient(xin)
    g <- .ojl_project(M, r[[2L]], xin)
    d$DF <- g
    list(r[[1L]], .ojl_dot(g, v))
  }
  w
}

# Stiefel retractions need the shape of the point: the dimensions of u0 are
# set on every manifold of a method and its nested methods
.ojl_with_dims <- function(alg, dims) {
  if (!is.list(alg)) return(alg)
  if (inherits(alg$manifold, "Manifold")) alg$manifold$dims <- dims
  for (nm in c("nlprecon", "method")) if (is.list(alg[[nm]])) alg[[nm]] <- .ojl_with_dims(alg[[nm]], dims)
  alg
}

# ---------------------------------------------------------------------------
# LineSearches.jl

.ojl_ls_fail <- function(msg, alpha)
  stop(structure(class = c("ojl_linesearch_exception", "error", "condition"),
                 list(message = msg, call = NULL, alpha = alpha)))

.ojl_ls <- function(type, ...) structure(list(type = type, ...), class = c(type, "LineSearchesMethod"))

#' Line searches and initial step guesses of LineSearches.jl
#'
#' `HagerZhang()`, `Static()`, `BackTracking()`, `StrongWolfe()` (Nocedal
#' and Wright, Algorithms 3.5 and 3.6) and `MoreThuente()` (More and Thuente,
#' 1994) are the line searches of LineSearches.jl, and `InitialStatic()`, `InitialPrevious()`,
#' `InitialHagerZhang()`, `InitialQuadratic()` and `InitialConstantChange()`
#' its initial step length guesses. They are reached as
#' `Optim$LineSearches$HagerZhang()` and so on (Julia's
#' `Optim.LineSearches.HagerZhang()`), and passed as the `linesearch` and
#' `alphaguess` arguments of the Optim methods.
#'
#' @param delta,sigma The Wolfe conditions of the Hager Zhang search.
#' @param alphamax The largest step.
#' @param rho The expansion factor of the bracketing phase.
#' @param epsilon The approximate Wolfe tolerance.
#' @param gamma The interval shrinking factor.
#' @param linesearchmax The largest number of iterations.
#' @param psi3 The factor that shrinks a step with a non finite value.
#' @param check_flatness Stop when the function is flat on the bracket.
#' @param c_1 The sufficient decrease constant.
#' @param c_2 The curvature constant of the strong Wolfe conditions.
#' @param f_tol,gtol The sufficient decrease and curvature constants of the
#'   More Thuente search.
#' @param x_tol The relative width of the interval of uncertainty at which
#'   the More Thuente search stops.
#' @param maxfev The largest number of evaluations of the More Thuente search.
#' @param rho_hi,rho_lo The safeguards of the backtracking step.
#' @param iterations The largest number of backtracking steps.
#' @param order The order of the interpolation (2 or 3).
#' @param maxstep The largest step norm.
#' @param alpha The step length.
#' @param scaled Scale the step by the norm of the direction.
#' @param alphamin The smallest step.
#' @param psi0,psi1,psi2 The parameters of the Hager Zhang initial guess.
#' @param alpha0 The first initial step.
#' @param quadstep Use the quadratic step of the Hager Zhang guess.
#' @param alpha_min,alpha_max,rho_0 The bounds and the largest decrease of
#'   the quadratic and constant change guesses.
#' @param snap2one Step lengths in this closed interval become 1.
#' @return A line search or initial guess object.
#' @examples
#' ls <- Optim$LineSearches$HagerZhang()
#' ag <- Optim$LineSearches$InitialStatic(alpha = 0.5)
#' prob <- OptimizationProblem(function(x, p) sum((x - 1)^2), c(0, 0),
#'                             grad = function(x, p) 2 * (x - 1))
#' solve(prob, GradientDescent(alphaguess = ag, linesearch = ls))$u
#' @name optim_linesearches
NULL

.ojl_HagerZhang <- function(delta = 0.1, sigma = 0.9, alphamax = Inf, rho = 5, epsilon = 1e-6, gamma = 2 / 3,
                            linesearchmax = 50L, psi3 = 0.1, check_flatness = FALSE)
  .ojl_ls("HagerZhang", delta = delta, sigma = sigma, alphamax = alphamax, rho = rho, epsilon = epsilon,
          gamma = gamma, linesearchmax = as.integer(linesearchmax), psi3 = psi3, check_flatness = check_flatness)

.ojl_Static <- function() .ojl_ls("Static")

.ojl_BackTracking <- function(c_1 = 1e-4, rho_hi = 0.5, rho_lo = 0.1, iterations = 1000L, order = 3L,
                              maxstep = Inf)
  .ojl_ls("BackTracking", c_1 = c_1, rho_hi = rho_hi, rho_lo = rho_lo, iterations = as.integer(iterations),
          order = as.integer(order), maxstep = maxstep)

.ojl_StrongWolfe <- function(c_1 = 1e-4, c_2 = 0.9, rho = 2) {
  if (!isTRUE(0 < c_1 && c_1 < c_2 && c_2 < 1))
    stop("The Wolfe constants must satisfy 0 < c_1 < c_2 < 1. Got c_1 = ", c_1, " and c_2 = ", c_2, call. = FALSE)
  if (!isTRUE(rho > 1))
    stop("The bracket growth factor must be larger than one, otherwise bracketing never terminates. Got \u03C1 = ",
         rho, call. = FALSE)
  .ojl_ls("StrongWolfe", c_1 = c_1, c_2 = c_2, rho = rho)
}

.ojl_MoreThuente <- function(f_tol = 1e-4, gtol = 0.9, x_tol = 1e-8, alphamin = 1e-16, alphamax = 65536,
                             maxfev = 100L)
  .ojl_ls("MoreThuente", f_tol = f_tol, gtol = gtol, x_tol = x_tol, alphamin = alphamin, alphamax = alphamax,
          maxfev = as.integer(maxfev))

.ojl_InitialStatic <- function(alpha = 1, scaled = FALSE) .ojl_ls("InitialStatic", alpha = alpha, scaled = scaled)

.ojl_InitialPrevious <- function(alpha = 1, alphamin = 0, alphamax = Inf)
  .ojl_ls("InitialPrevious", alpha = alpha, alphamin = alphamin, alphamax = alphamax)

.ojl_InitialHagerZhang <- function(psi0 = 0.01, psi1 = 0.2, psi2 = 2, psi3 = 0.1, alphamax = Inf, alpha0 = 1,
                                   quadstep = TRUE)
  .ojl_ls("InitialHagerZhang", psi0 = psi0, psi1 = psi1, psi2 = psi2, psi3 = psi3, alphamax = alphamax,
          alpha0 = alpha0, quadstep = quadstep)

.ojl_InitialQuadratic <- function(alpha_min = 1e-12, alpha_max = 1, alpha0 = 1, rho_0 = 0.25, snap2one = c(0.75, Inf))
  .ojl_ls("InitialQuadratic", amin = alpha_min, amax = alpha_max, a0 = alpha0, rho = rho_0, snap2one = snap2one)

.ojl_InitialConstantChange <- function(alpha_min = 1e-12, alpha_max = 1, alpha0 = 1, rho_0 = 0.25,
                                       snap2one = c(0.75, Inf)) {
  r <- .ojl_ls("InitialConstantChange", amin = alpha_min, amax = alpha_max, a0 = alpha0, rho = rho_0,
               snap2one = snap2one)
  r$prev <- new.env(parent = emptyenv())
  r$prev$dphi <- NaN
  r
}

.ojl_alphaguess <- function(a) if (is.numeric(a)) .ojl_InitialStatic(alpha = a) else a

# HagerZhang (hagerzhang.jl): returns c(step, value)
.ojl_hz <- function(ls, phidphi, c, phi_0, dphi_0) {
  delta <- ls$delta; sigma <- ls$sigma; rho <- ls$rho; epsilon <- ls$epsilon; gamma <- ls$gamma
  lsmax <- ls$linesearchmax; psi3 <- ls$psi3; alphamax <- ls$alphamax
  if (!(is.finite(phi_0) && is.finite(dphi_0))) .ojl_ls_fail("Value and slope at step length = 0 must be finite.", 0)
  if (dphi_0 >= 0) .ojl_ls_fail("Search direction is not a direction of descent.", 0)
  itfmax <- 52L
  alphas <- 0; values <- phi_0; slopes <- dphi_0
  phi_lim <- phi_0 + epsilon * abs(phi_0)
  if (c < 0) .ojl_ls_fail("Initial step length must be non-negative.", 0)
  if (c <= .Machine$double.eps) .ojl_ls_fail("Initial step length must exceed machine epsilon.", c)
  if (!is.finite(c)) .ojl_ls_fail("Initial step length must be finite.", 0)
  if (c > alphamax) .ojl_ls_fail("Initial step length exceeds alphamax.", 0)
  wolfe <- function(c, phi_c, dphi_c) {
    w1 <- isTRUE(delta * dphi_0 >= (phi_c - phi_0) / c) && isTRUE(dphi_c >= sigma * dphi_0)
    w2 <- isTRUE((2 * delta - 1) * dphi_0 >= dphi_c) && isTRUE(dphi_c >= sigma * dphi_0) && isTRUE(phi_c <= phi_lim)
    w1 || w2
  }
  push <- function(a, r) {
    alphas <<- c(alphas, a); values <<- c(values, r[1L]); slopes <<- c(slopes, r[2L])
    length(alphas)
  }
  bisect <- function(ia, ib) {
    a <- alphas[ia]; b <- alphas[ib]
    while (b - a > .ojl_ulp(b)) {
      d <- (a + b) / 2
      r <- phidphi(d)
      id <- push(d, r)
      if (wolfe(d, r[1L], r[2L])) return(c(ia, id, 1L))
      if (isTRUE(r[2L] >= 0)) return(c(ia, id, 0L))
      if (isTRUE(r[1L] <= phi_lim) && is.finite(r[2L])) {
        a <- d; ia <- id
      } else {
        b <- d; ib <- id
      }
    }
    c(ia, ib, 0L)
  }
  update <- function(ia, ib, ic) {
    a <- alphas[ia]; b <- alphas[ib]; c <- alphas[ic]
    if (!isTRUE(a < c && c < b)) return(c(ia, ib, 0L))
    if (isTRUE(slopes[ic] >= 0)) return(c(ia, ic, 0L))
    if (isTRUE(values[ic] <= phi_lim) && is.finite(slopes[ic])) return(c(ic, ib, 0L))
    bisect(ia, ic)
  }
  secant <- function(a, b, da, db) (a * db - b * da) / (db - da)
  secant2 <- function(ia, ib) {
    a <- alphas[ia]; b <- alphas[ib]; da <- slopes[ia]; db <- slopes[ib]
    if (!isTRUE(da < 0))
      stop("Search direction is not a direction of descent; this error may indicate that user-provided ",
           "derivatives are inaccurate.", call. = FALSE)
    c <- secant(a, b, da, db)
    if (!isTRUE(a < c && c < b)) c <- (a + b) / 2
    r <- phidphi(c)
    ic <- push(c, r)
    if (wolfe(c, r[1L], r[2L])) return(c(1L, ic, ic))
    u <- update(ia, ib, ic)
    if (u[3L] == 1L) return(c(1L, u[2L], u[2L]))
    iA <- u[1L]; iB <- u[2L]
    a <- alphas[iA]; b <- alphas[iB]
    if (iB == ic) c <- secant(alphas[ib], alphas[iB], slopes[ib], slopes[iB])
    else if (iA == ic) c <- secant(alphas[ia], alphas[iA], slopes[ia], slopes[iA])
    if ((iA == ic || iB == ic) && isTRUE(a < c && c < b)) {
      r <- phidphi(c)
      ic <- push(c, r)
      if (wolfe(c, r[1L], r[2L])) return(c(1L, ic, ic))
      u <- update(iA, iB, ic)
      if (u[3L] == 1L) return(c(1L, u[2L], u[2L]))
      iA <- u[1L]; iB <- u[2L]
    }
    c(0L, iA, iB)
  }
  r <- phidphi(c)
  iterfinite <- 1L
  while (!(is.finite(r[1L]) && is.finite(r[2L])) && iterfinite < itfmax) {
    iterfinite <- iterfinite + 1L
    c <- c * psi3
    r <- phidphi(c)
  }
  if (!(is.finite(r[1L]) && is.finite(r[2L]))) .ojl_ls_fail("Failed to achieve finite new evaluation point.", c)
  phi_c <- r[1L]; dphi_c <- r[2L]
  push(c, r)
  if (wolfe(c, phi_c, dphi_c)) return(c(c, phi_c))
  isbracketed <- FALSE
  ia <- 1L; ib <- 2L
  iter <- 1L
  while (!isbracketed && iter < lsmax) {
    if (dphi_c >= 0) {
      ib <- length(alphas)
      for (i in (ib - 1L):1L) if (values[i] <= phi_lim) {
        ia <- i
        break
      }
      isbracketed <- TRUE
    } else if (values[length(values)] > phi_lim) {
      ib <- length(alphas); ia <- 1L
      b <- bisect(ia, ib)
      ia <- b[1L]; ib <- b[2L]
      if (b[3L] == 1L) return(c(alphas[ib], values[ib]))
      isbracketed <- TRUE
    } else {
      cold <- c; phi_cold <- phi_c
      if (.ojl_nextfloat(cold) >= alphamax) {
        if (alphamax < ls$alphamax)
          .ojl_ls_fail("Failed to bracket: pinned at non-finite feasibility boundary.", cold)
        return(c(cold, phi_cold))
      }
      c <- c * rho
      if (c > alphamax) c <- alphamax
      r <- phidphi(c)
      iterfinite <- 1L
      while (!(is.finite(r[1L]) && is.finite(r[2L])) && c > .ojl_nextfloat(cold) && iterfinite < itfmax) {
        alphamax <- c
        iterfinite <- iterfinite + 1L
        c <- (cold + c) / 2
        r <- phidphi(c)
      }
      if (!(is.finite(r[1L]) && is.finite(r[2L]))) {
        if (alphamax < ls$alphamax)
          .ojl_ls_fail("Failed to expand bracket to finite values after alphamax was shrunk.", cold)
        return(c(cold, phi_cold))
      }
      phi_c <- r[1L]; dphi_c <- r[2L]
      push(c, r)
      if (wolfe(c, phi_c, dphi_c)) return(c(c, phi_c))
    }
    iter <- iter + 1L
  }
  while (iter < lsmax) {
    a <- alphas[ia]; b <- alphas[ib]
    if (b - a <= .ojl_ulp(b)) return(c(a, values[ia]))
    s2 <- secant2(ia, ib)
    if (s2[1L] == 1L) return(c(alphas[s2[2L]], values[s2[2L]]))
    iA <- s2[2L]; iB <- s2[3L]
    A <- alphas[iA]; B <- alphas[iB]
    if (B - A < gamma * (b - a)) {
      if (ls$check_flatness && .ojl_nextfloat(values[ia]) >= values[ib] &&
          .ojl_nextfloat(values[iA]) >= values[iB]) return(c(A, values[iA]))
      ia <- iA; ib <- iB
    } else {
      c <- (A + B) / 2
      r <- phidphi(c)
      ic <- push(c, r)
      if (wolfe(c, r[1L], r[2L])) return(c(c, r[1L]))
      u <- update(iA, iB, ic)
      ia <- u[1L]; ib <- u[2L]
      if (u[3L] == 1L) return(c(alphas[ib], values[ib]))
    }
    iter <- iter + 1L
  }
  best <- alphas[which.min(values)]
  .ojl_ls_fail(sprintf("Linesearch failed to converge, reached maximum iterations %d.", lsmax), best)
}

# Static (static.jl)
.ojl_static <- function(phi, alpha) {
  if (!(alpha > 0)) stop("the step of the Static line search must be positive", call. = FALSE)
  pa <- phi(alpha)
  it <- 0L
  while (!is.finite(pa) && it < 52L) {
    it <- it + 1L
    alpha <- alpha / 2
    pa <- phi(alpha)
  }
  if (!is.finite(pa)) .ojl_ls_fail("Static: failed to achieve finite new evaluation point.", alpha)
  c(alpha, pa)
}

# BackTracking (backtracking.jl)
.ojl_backtracking <- function(ls, phi, s, alpha_0, phi_0, dphi_0) {
  c_1 <- ls$c_1; rho_hi <- ls$rho_hi; rho_lo <- ls$rho_lo; order <- ls$order; iterations <- ls$iterations
  alpha_0 <- min(alpha_0, min(Inf, ls$maxstep / max(abs(s))))
  if (!(is.finite(phi_0) && is.finite(dphi_0)))
    .ojl_ls_fail("Value and slope at step length = 0 must be finite.", 0)
  if (dphi_0 > 0) .ojl_ls_fail("Search direction is not a direction of descent.", 0)
  iteration <- 0L
  phix_0 <- phi_0
  alpha_1 <- alpha_2 <- alpha_0
  phix_1 <- phi(alpha_2)
  iterfinite <- 0L
  while (!is.finite(phix_1) && iterfinite < 52L) {
    iterfinite <- iterfinite + 1L
    alpha_2 <- alpha_2 / 2
    phix_1 <- phi(alpha_2)
  }
  if (!is.finite(phix_1)) .ojl_ls_fail("Backtracking: failed to achieve finite new evaluation point.", alpha_2)
  while (phix_1 > phi_0 + c_1 * alpha_2 * dphi_0) {
    iteration <- iteration + 1L
    if (iteration > iterations)
      .ojl_ls_fail(sprintf("Linesearch failed to converge, reached maximum iterations %d.", iterations), alpha_2)
    if (order == 2L || iteration == 1L) {
      alpha_tmp <- -(dphi_0 * (alpha_2 * alpha_2)) / (2 * (phix_1 - phi_0 - dphi_0 * alpha_2))
    } else {
      a1s <- alpha_1 * alpha_1; a2s <- alpha_2 * alpha_2
      div <- 1 / (a1s * a2s * (alpha_2 - alpha_1))
      a <- (a1s * (phix_1 - phi_0 - dphi_0 * alpha_2) - a2s * (phix_0 - phi_0 - dphi_0 * alpha_1)) * div
      b <- (-(a1s * alpha_1) * (phix_1 - phi_0 - dphi_0 * alpha_2) +
              (a2s * alpha_2) * (phix_0 - phi_0 - dphi_0 * alpha_1)) * div
      if (a == 0) {
        alpha_tmp <- -dphi_0 / (2 * b)
      } else {
        dd <- max(b * b - 3 * a * dphi_0, 0)
        alpha_tmp <- -dphi_0 / (b + sqrt(dd))
      }
    }
    alpha_1 <- alpha_2
    alpha_tmp <- .ojl_nanmin(alpha_tmp, alpha_2 * rho_hi)
    alpha_2 <- .ojl_nanmax(alpha_tmp, alpha_2 * rho_lo)
    phix_new <- phi(alpha_2)
    while (!is.finite(phix_new) && iteration < iterations) {
      iteration <- iteration + 1L
      alpha_2 <- alpha_2 / 2
      phix_new <- phi(alpha_2)
    }
    if (!is.finite(phix_new)) .ojl_ls_fail("Backtracking: failed to achieve finite new evaluation point.", alpha_2)
    phix_0 <- phix_1
    phix_1 <- phix_new
  }
  c(alpha_2, phix_1)
}

.ojl_nanmin <- function(a, b) if (is.nan(a)) b else if (is.nan(b)) a else min(a, b)

# cubic_d1_d2 (LineSearches.jl): the terms d1 and d2 of the Hermite cubic
# through two points (Nocedal and Wright, eq. 3.59)
.ojl_cubic_d1_d2 <- function(a1, p1, dp1, a2, p2, dp2) {
  d1 <- 3 * (p1 - p2) / (a2 - a1) + dp1 + dp2
  s <- max(abs(d1), abs(dp1), abs(dp2))
  t <- d1 / s
  v <- t * t - (dp1 / s) * (dp2 / s)
  c(d1, s * sqrt(if (is.nan(v)) 0 else max(0, v)))
}

# StrongWolfe (strongwolfe.jl): Algorithms 3.5 and 3.6 of Nocedal and Wright
.ojl_strongwolfe <- function(ls, phidphi, alpha0, phi_0, dphi_0) {
  c_1 <- ls$c_1; c_2 <- ls$c_2
  if (!(is.finite(phi_0) && is.finite(dphi_0)))
    .ojl_ls_fail("Value and slope at step length = 0 must be finite.", 0)
  if (dphi_0 > 0) .ojl_ls_fail("Search direction is not a direction of descent.", 0)
  a_prev <- 0
  a_i <- alpha0
  a_max <- 65536
  p_prev <- phi_0
  dp_prev <- dphi_0
  i <- 1L
  while (isTRUE(a_i < a_max)) {
    v <- phidphi(a_i)
    if (isTRUE(v[1L] > phi_0 + c_1 * a_i * dphi_0) || (isTRUE(v[1L] >= p_prev) && i > 1L))
      return(.ojl_sw_zoom(a_prev, p_prev, dp_prev, a_i, v[1L], v[2L], dphi_0, phi_0, phidphi, c_1, c_2))
    if (isTRUE(abs(v[2L]) <= -c_2 * dphi_0)) return(c(a_i, v[1L]))
    if (isTRUE(v[2L] >= 0))
      return(.ojl_sw_zoom(a_i, v[1L], v[2L], a_prev, p_prev, dp_prev, dphi_0, phi_0, phidphi, c_1, c_2))
    a_prev <- a_i
    p_prev <- v[1L]
    dp_prev <- v[2L]
    a_i <- a_i * ls$rho
    i <- i + 1L
  }
  .ojl_ls_fail("StrongWolfe: bracketing reached a_max=65536.0 without satisfying Wolfe conditions.", a_max)
}

.ojl_sw_zoom <- function(a_lo, p_lo, dp_lo, a_hi, p_hi, dp_hi, dphi_0, phi_0, phidphi, c_1, c_2) {
  a_j <- NaN
  for (iteration in 1:10) {
    a_j <- if (isTRUE(a_lo < a_hi)) .ojl_sw_interpolate(a_lo, a_hi, p_lo, p_hi, dp_lo, dp_hi)
           else .ojl_sw_interpolate(a_hi, a_lo, p_hi, p_lo, dp_hi, dp_lo)
    v <- phidphi(a_j)
    if (isTRUE(v[1L] > phi_0 + c_1 * a_j * dphi_0) || isTRUE(v[1L] > p_lo)) {
      a_hi <- a_j; p_hi <- v[1L]; dp_hi <- v[2L]
    } else {
      if (isTRUE(abs(v[2L]) <= -c_2 * dphi_0)) return(c(a_j, v[1L]))
      if (isTRUE(v[2L] * (a_hi - a_lo) >= 0)) {
        a_hi <- a_lo; p_hi <- p_lo; dp_hi <- dp_lo
      }
      a_lo <- a_j; p_lo <- v[1L]; dp_lo <- v[2L]
    }
  }
  .ojl_ls_fail("StrongWolfe: zoom reached maximum iterations (10) without satisfying Wolfe conditions.", a_j)
}

# interpolate of strongwolfe.jl: the minimizer of the cubic, or the midpoint
.ojl_sw_interpolate <- function(a1, a2, p1, p2, dp1, dp2) {
  d <- .ojl_cubic_d1_d2(a1, p1, dp1, a2, p2, dp2)
  a_j <- a2 - (a2 - a1) * ((dp2 + d[2L] - d[1L]) / (dp2 - dp1 + 2 * d[2L]))
  if (isTRUE(min(a1, a2) < a_j) && isTRUE(a_j < max(a1, a2))) a_j else (a1 + a2) / 2
}

# MoreThuente (morethuente.jl): the line search of More and Thuente (1994),
# from MINPACK's cvsrch and cstep
.ojl_morethuente <- function(ls, phidphi, alpha, phi_0, dphi_0) {
  f_tol <- ls$f_tol; gtol <- ls$gtol; x_tol <- ls$x_tol; alphamin <- ls$alphamin; maxfev <- ls$maxfev
  alphamax <- ls$alphamax
  iterfinitemax <- -log2(.Machine$double.eps)
  fail <- .ojl_ls_fail
  if (!(is.finite(phi_0) && is.finite(dphi_0))) fail("Value and slope at step length = 0 must be finite.", 0)
  if (isTRUE(alpha <= 0) || f_tol < 0 || gtol < 0 || x_tol < 0 || alphamin < 0 || alphamax < alphamin || maxfev <= 0)
    fail("Invalid parameters to MoreThuente.", 0)
  if (dphi_0 >= 0) fail("Search direction is not a direction of descent.", 0)
  if (!is.finite(alpha)) fail("Initial step length must be finite.", alpha)
  bracketed <- FALSE
  stage1 <- TRUE
  nfev <- 0L
  dgtest <- f_tol * dphi_0
  width <- alphamax - alphamin
  width1 <- 2 * width
  stx <- 0; fx <- phi_0; dgx <- dphi_0
  sty <- 0; fy <- phi_0; dgy <- dphi_0
  stmin <- stx
  stmax <- alpha + 4 * (alpha - stx)
  alpha <- min(max(alpha, alphamin), alphamax)
  v <- phidphi(alpha)
  f <- v[1L]; dg <- v[2L]
  nfev <- nfev + 1L
  iterfinite <- 0L
  while ((!is.finite(f) || !is.finite(dg)) && iterfinite < iterfinitemax) {
    iterfinite <- iterfinite + 1L
    alpha <- alpha / 2
    v <- phidphi(alpha)
    f <- v[1L]; dg <- v[2L]
    nfev <- nfev + 1L
  }
  if (!(is.finite(f) && is.finite(dg))) fail("Failed to achieve finite new evaluation point.", alpha)
  alpha_checked <- alpha
  firstpass <- TRUE
  repeat {
    if (bracketed) {
      stmin <- min(stx, sty)
      stmax <- max(stx, sty)
    } else {
      stmin <- stx
      stmax <- alpha + 4 * (alpha - stx)
    }
    stmin <- max(alphamin, stmin)
    stmax <- min(alphamax, stmax)
    alpha <- min(max(alpha, alphamin), alphamax)
    if ((bracketed && (isTRUE(alpha <= stmin) || isTRUE(alpha >= stmax))) || nfev >= maxfev - 1L ||
        (bracketed && isTRUE(stmax - stmin <= x_tol * stmax))) alpha <- stx
    if (!(firstpass && isTRUE(alpha == alpha_checked))) {
      v <- phidphi(alpha)
      f <- v[1L]; dg <- v[2L]
      nfev <- nfev + 1L
    }
    firstpass <- FALSE
    iterfinite <- 0L
    while (!(is.finite(f) && is.finite(dg))) {
      if (nfev >= maxfev)
        fail(sprintf(paste0("MoreThuente: reached maximum function evaluations (%d) while recovering from ",
                            "non-finite evaluation."), maxfev), alpha)
      if (iterfinite >= iterfinitemax || isTRUE(alpha == stx))
        fail("MoreThuente: failed to recover finite evaluation by bisecting toward stx.", alpha)
      alphamax <- min(alphamax, alpha)
      iterfinite <- iterfinite + 1L
      alpha <- (alpha + stx) / 2
      v <- phidphi(alpha)
      f <- v[1L]; dg <- v[2L]
      nfev <- nfev + 1L
    }
    ftest1 <- phi_0 + alpha * dgtest
    if (f <= ftest1 && abs(dg) <= -gtol * dphi_0) return(c(alpha, f))
    if (bracketed && isTRUE(stmax - stmin <= x_tol * stmax))
      fail(sprintf("MoreThuente: bracket width fell below x_tol=%s without satisfying Wolfe conditions.",
                   format(x_tol)), alpha)
    if (nfev >= maxfev)
      fail(sprintf("MoreThuente: reached maximum function evaluations (%d) without satisfying Wolfe conditions.",
                   maxfev), alpha)
    if (alpha == alphamin && (f > ftest1 || dg >= dgtest))
      fail("MoreThuente: step pinned at alphamin without satisfying Wolfe conditions.", alpha)
    if (alpha == alphamax && f <= ftest1 && dg <= dgtest) return(c(alpha, f))
    if (bracketed && (alpha <= stmin || alpha >= stmax))
      fail("MoreThuente: rounding errors prevent further progress; tolerances may be too small.", alpha)
    if (stage1 && f <= ftest1 && dg >= min(f_tol, gtol) * dphi_0) stage1 <- FALSE
    if (stage1 && f <= fx && f > ftest1) {
      # the modified function psi(alpha) = phi(alpha) - phi_0 - alpha * dgtest
      cs <- .ojl_cstep(stx, fx - stx * dgtest, dgx - dgtest, sty, fy - sty * dgtest, dgy - dgtest,
                       alpha, f - alpha * dgtest, dg - dgtest, bracketed, stmin, stmax)
      stx <- cs[[1L]]; sty <- cs[[4L]]; alpha <- cs[[7L]]; bracketed <- cs[[10L]]
      fx <- cs[[2L]] + stx * dgtest
      fy <- cs[[5L]] + sty * dgtest
      dgx <- cs[[3L]] + dgtest
      dgy <- cs[[6L]] + dgtest
    } else {
      cs <- .ojl_cstep(stx, fx, dgx, sty, fy, dgy, alpha, f, dg, bracketed, stmin, stmax)
      stx <- cs[[1L]]; fx <- cs[[2L]]; dgx <- cs[[3L]]
      sty <- cs[[4L]]; fy <- cs[[5L]]; dgy <- cs[[6L]]
      alpha <- cs[[7L]]; bracketed <- cs[[10L]]
    }
    if (bracketed) {
      if (isTRUE(abs(sty - stx) >= (2 / 3) * width1)) alpha <- stx + (sty - stx) / 2
      width1 <- width
      width <- abs(sty - stx)
    }
  }
}

# cstep: a safeguarded step of the More and Thuente search and the update of
# the interval of uncertainty
.ojl_cstep <- function(stx, fx, dgx, sty, fy, dgy, alpha, f, dg, bracketed, alphamin, alphamax) {
  if ((bracketed && (isTRUE(alpha <= min(stx, sty)) || isTRUE(alpha >= max(stx, sty)))) ||
      isTRUE(dgx * (alpha - stx) >= 0) || alphamax < alphamin)
    .ojl_ls_fail(paste0("MoreThuente: cstep called with alpha = ", format(alpha), " outside the bracket [",
                        format(min(stx, sty)), ", ", format(max(stx, sty)), "]."), alpha)
  sgnd <- dg * (dgx / abs(dgx))
  if (f > fx) {
    bound <- TRUE
    tg <- .ojl_cubic_d1_d2(stx, fx, dgx, alpha, f, dg)
    theta <- tg[1L]; gamma <- tg[2L]
    if (alpha < stx) gamma <- -gamma
    p <- gamma - dgx + theta
    q <- gamma - dgx + gamma + dg
    r <- p / q
    alphac <- stx + r * (alpha - stx)
    alphaq <- stx + (dgx / ((fx - f) / (alpha - stx) + dgx)) / 2 * (alpha - stx)
    alphaf <- if (isTRUE(abs(alphac - stx) < abs(alphaq - stx))) alphac else (alphac + alphaq) / 2
    bracketed <- TRUE
  } else if (isTRUE(sgnd < 0)) {
    bound <- FALSE
    tg <- .ojl_cubic_d1_d2(stx, fx, dgx, alpha, f, dg)
    theta <- tg[1L]; gamma <- tg[2L]
    if (alpha > stx) gamma <- -gamma
    p <- gamma - dg + theta
    q <- gamma - dg + gamma + dgx
    r <- p / q
    alphac <- alpha + r * (stx - alpha)
    alphaq <- alpha + (dg / (dg - dgx)) * (stx - alpha)
    alphaf <- if (isTRUE(abs(alphac - alpha) > abs(alphaq - alpha))) alphac else alphaq
    bracketed <- TRUE
  } else if (isTRUE(abs(dg) < abs(dgx))) {
    bound <- TRUE
    tg <- .ojl_cubic_d1_d2(stx, fx, dgx, alpha, f, dg)
    theta <- tg[1L]; gamma <- tg[2L]
    if (alpha > stx) gamma <- -gamma
    p <- gamma - dg + theta
    q <- gamma + dgx - dg + gamma
    r <- p / q
    alphac <- if (isTRUE(r < 0) && isTRUE(gamma != 0)) alpha + r * (stx - alpha)
              else if (alpha > stx) alphamax else alphamin
    alphaq <- alpha + (dg / (dg - dgx)) * (stx - alpha)
    closer <- isTRUE(abs(alpha - alphac) < abs(alpha - alphaq))
    farther <- isTRUE(abs(alpha - alphac) > abs(alpha - alphaq))
    alphaf <- if (if (bracketed) closer else farther) alphac else alphaq
  } else {
    bound <- FALSE
    if (bracketed) {
      tg <- .ojl_cubic_d1_d2(sty, fy, dgy, alpha, f, dg)
      theta <- tg[1L]; gamma <- tg[2L]
      if (alpha > sty) gamma <- -gamma
      p <- gamma - dg + theta
      q <- gamma - dg + gamma + dgy
      r <- p / q
      alphaf <- alpha + r * (sty - alpha)
    } else {
      alphaf <- if (alpha > stx) alphamax else alphamin
    }
  }
  if (f > fx) {
    sty <- alpha; fy <- f; dgy <- dg
  } else {
    if (isTRUE(sgnd < 0)) {
      sty <- stx; fy <- fx; dgy <- dgx
    }
    stx <- alpha; fx <- f; dgx <- dg
  }
  alpha <- max(alphamin, min(alphamax, alphaf))
  if (bracketed && bound) {
    alpha <- if (sty > stx) min(stx + (2 / 3) * (sty - stx), alpha) else max(stx + (2 / 3) * (sty - stx), alpha)
  }
  list(stx, fx, dgx, sty, fy, dgy, alpha, f, dg, bracketed)
}

# the initial step guesses (initialguess.jl) set st$alpha
.ojl_alphaguess_apply <- function(ig, ls, st, phi_0, dphi_0, d) {
  switch(ig$type,
    InitialStatic = {
      ns <- sqrt(sum(st$s * st$s))
      st$alpha <- if (isTRUE(ig$scaled) && ns > 0) min(ig$alpha, ns) / ns else ig$alpha
    },
    InitialPrevious = {
      if (is.nan(st$alpha)) st$alpha <- ig$alpha
      st$alpha <- min(ig$alphamax, st$alpha)
      st$alpha <- max(ig$alphamin, st$alpha)
    },
    InitialQuadratic = {
      if (!is.finite(st$f_x_previous) || dphi_0 == 0) {
        ag <- ig$a0
      } else {
        ag <- 2 * (phi_0 - st$f_x_previous) / dphi_0
        ag <- .ojl_nanmax(.ojl_nanmax(ig$amin, st$alpha * ig$rho), ag)
        ag <- .ojl_nanmin(ig$amax, ag)
        if (isTRUE(ig$snap2one[1L] <= ag && ag <= ig$snap2one[2L])) ag <- 1
      }
      st$alpha <- ag
    },
    InitialConstantChange = {
      prev <- ig$prev$dphi
      if (!is.finite(prev) || !is.finite(st$alpha) || dphi_0 == 0) {
        ag <- ig$a0
      } else {
        ag <- st$alpha * prev / dphi_0
        ag <- .ojl_nanmax(.ojl_nanmax(ig$amin, st$alpha * ig$rho), ag)
        ag <- .ojl_nanmin(ig$amax, ag)
        if (isTRUE(ig$snap2one[1L] <= ag && ag <= ig$snap2one[2L])) ag <- 1
      }
      ig$prev$dphi <- dphi_0
      st$alpha <- ag
    },
    InitialHagerZhang = {
      if (is.nan(st$f_x_previous) && is.nan(ig$alpha0)) {
        st$alpha <- .ojl_hzI0(st$x, st$g_x, phi_0, ig$alphamax, ig$psi0)
      } else {
        st$alpha <- .ojl_hzI12(st$alpha, d, st, phi_0, dphi_0, ig$psi1, ig$psi2, ig$psi3, ig$alphamax,
                               ig$quadstep)
      }
    },
    stop("unknown initial step guess ", ig$type, call. = FALSE))
  invisible(NULL)
}

.ojl_hzI0 <- function(x, gr, f_x, alphamax, psi0) {
  alpha <- 1
  gr_max <- max(abs(gr))
  if (gr_max != 0) {
    x_max <- max(abs(x))
    if (x_max != 0) {
      alpha <- psi0 * x_max / gr_max
    } else if (f_x != 0) {
      ng <- sqrt(sum(gr * gr))
      alpha <- psi0 * abs(f_x) / (ng * ng)
    }
  }
  min(alpha, alphamax)
}

.ojl_hzI12 <- function(alpha, d, st, phi_0, dphi_0, psi1, psi2, psi3, alphamax, quadstep) {
  x <- st$x; s <- st$s
  phi <- function(a) {
    xn <- .ojl_fma(a, s, x)
    st$x_ls <- xn
    d$value_bang(xn)
  }
  get_finite <- function(a, ph) {
    it <- 1L
    while (!is.finite(ph)) {
      if (it >= 52L) return(c(0, phi_0))
      a <- psi3 * a
      ph <- phi(a)
      it <- it + 1L
    }
    c(a, ph)
  }
  alphatest <- psi1 * alpha
  alphatest <- min(alphatest, alphamax)
  if (alphatest == 0) return(alphatest)
  phitest <- phi(alphatest)
  r <- get_finite(alphatest, phitest)
  alphatest <- r[1L]; phitest <- r[2L]
  if (alphatest == 0) return(alphatest)
  quadstep_success <- FALSE
  if (quadstep) {
    a <- ((phitest - phi_0) / alphatest - dphi_0) / alphatest
    if (is.finite(a) && a > 0 && phitest <= phi_0) {
      alphatest2 <- -dphi_0 / 2 / a
      alphatest2 <- min(alphatest2, alphamax)
      phitest2 <- phi(alphatest2)
      if (is.finite(phitest2)) {
        quadstep_success <- TRUE
        alphatest <- alphatest2
        phitest <- phitest2
      }
    }
  }
  if ((!quadstep || !quadstep_success) && phitest <= phi_0) {
    alphatest <- psi2 * alpha
    alphatest <- min(alphatest, alphamax)
    phitest <- phi(alphatest)
    r <- get_finite(alphatest, phitest)
    alphatest <- r[1L]
  }
  alphatest
}

# one call of the line search of a method from st$x along st$s; sets st$x_ls
.ojl_linesearch_call <- function(ls, d, st, phi_0, dphi_0) {
  x <- st$x; s <- st$s
  # make_phidphi: the value and the slope along s
  phidphi <- function(a) {
    xn <- .ojl_fma(a, s, x)
    st$x_ls <- xn
    v <- d$value_jvp(xn, s)
    c(v[[1L]], v[[2L]])
  }
  r <- switch(ls$type,
    HagerZhang = .ojl_hz(ls, phidphi, st$alpha, phi_0, dphi_0),
    StrongWolfe = .ojl_strongwolfe(ls, phidphi, st$alpha, phi_0, dphi_0),
    MoreThuente = .ojl_morethuente(ls, phidphi, st$alpha, phi_0, dphi_0),
    Static = .ojl_static(function(a) {
      xn <- .ojl_fma(a, s, x)
      st$x_ls <- xn
      d$value_bang(xn)
    }, st$alpha),
    BackTracking = .ojl_backtracking(ls, function(a) {
      xn <- .ojl_fma(a, s, x)
      st$x_ls <- xn
      d$value_bang(xn)
    }, s, st$alpha, phi_0, dphi_0),
    stop("unknown line search ", ls$type, call. = FALSE))
  st$x_ls <- .ojl_fma(r[1L], s, x)
  r
}

# perform_linesearch! (utilities/perform_linesearch.jl)
.ojl_perform_linesearch <- function(st, m, d) {
  # the methods on a manifold search along ManifoldObjective(manifold, d)
  d <- .ojl_manifold_objective(m$manifold, d)
  dphi_0 <- .ojl_dot(st$g_x, st$s)
  if (dphi_0 >= 0 && .ojl_reset_direction(st, m)) dphi_0 <- .ojl_dot(st$g_x, st$s)
  phi_0 <- st$f_x
  attempt <- function(dphi_0) {
    .ojl_alphaguess_apply(m$alphaguess, m$linesearch, st, phi_0, dphi_0, d)
    tryCatch({
      r <- .ojl_linesearch_call(m$linesearch, d, st, phi_0, dphi_0)
      st$alpha <- r[1L]
      if (!(st$alpha > 0)) .ojl_ls_fail("Line search returned non-positive alpha.", st$alpha)
      TRUE
    }, ojl_linesearch_exception = function(e) {
      st$alpha <- e$alpha
      FALSE
    })
  }
  ok <- attempt(dphi_0)
  if (!ok && .ojl_reset_direction(st, m)) ok <- attempt(.ojl_dot(st$g_x, st$s))
  st$f_x_previous <- phi_0
  st$x_previous <- st$x
  if (!is.null(st$g_x_previous)) st$g_x_previous <- st$g_x
  ok
}

# reset_search_direction!
.ojl_reset_direction <- function(st, m) {
  switch(m$name,
    BFGS = {
      st$invH <- .ojl_bfgs_invH0(m, st$x, st$g_x)
      st$s <- -st$g_x
      TRUE
    },
    LBFGS = {
      st$pseudo_iteration <- 1L
      st$s <- -st$g_x
      TRUE
    },
    ConjugateGradient = {
      st$s <- -st$pg
      TRUE
    },
    FALSE)
}

# ---------------------------------------------------------------------------
# Preconditioners (multivariate/precon.jl)

.ojl_precondprep <- function(m, x) {
  if (is.null(m$P) || !is.function(m$precondprep)) return(invisible(NULL))
  r <- m$precondprep(m$P, x)
  if (!inherits(m$P, "InverseDiagonal") && !is.null(r)) m$P <- r
  invisible(NULL)
}

.ojl_precondition_apply <- function(P, g) {
  if (is.null(P)) g
  else if (inherits(P, "InverseDiagonal")) P$diag * g
  else if (is.matrix(P)) as.numeric(solve(P, g))
  else g / P
}

.ojl_inverse_precondition <- function(P, s) {
  if (is.null(P)) .ojl_dot(s, s)
  else if (inherits(P, "InverseDiagonal")) .ojl_dot(s, s / P$diag)
  else if (is.matrix(P)) sum(s * as.numeric(P %*% s))
  else sum(s * P * s)
}

# ---------------------------------------------------------------------------
# Options and the optimization loop (types.jl, optimize.jl)

.ojl_default_options <- list(
  x_abstol = 0, x_reltol = 0, f_abstol = 0, f_reltol = 0, g_abstol = 1e-8,
  outer_x_abstol = 0, outer_x_reltol = 0, outer_f_abstol = 0, outer_f_reltol = 0, outer_g_abstol = 1e-8,
  f_calls_limit = 0L, g_calls_limit = 0L, h_calls_limit = 0L, allow_f_increases = TRUE,
  allow_outer_f_increases = TRUE, successive_f_tol = 1L, iterations = 1000L, outer_iterations = 1000L,
  show_trace = FALSE, show_warnings = TRUE, callback = NULL, time_limit = NaN)

# OptimizationOptimJL builds Optim.Options from the solver arguments alone, so
# the per method default_options of Optim (successive_f_tol = 2 for IPNewton)
# do not apply on this path
.ojl_options <- function(..., method = NULL) {
  o <- .ojl_default_options
  args <- list(...)
  for (nm in names(args)) if (!is.null(args[[nm]])) o[[nm]] <- args[[nm]]
  o$iterations <- as.integer(o$iterations); o$outer_iterations <- as.integer(o$outer_iterations)
  o
}

.ojl_g_residual <- function(st) if (is.null(st$g_x)) NaN else max(abs(st$g_x))

.ojl_assess <- function(x, x_previous, f_x, f_x_previous, g_x, x_abstol, x_reltol, f_abstol, f_reltol, g_abstol) {
  dx <- max(abs(x - x_previous))
  df <- abs(f_x - f_x_previous)
  c(x = isTRUE(dx <= x_abstol) || isTRUE(dx <= x_reltol * max(abs(x))),
    f = isTRUE(df <= f_abstol) || isTRUE(df <= f_reltol * abs(f_x)),
    g = isTRUE(max(abs(g_x)) <= g_abstol),
    inc = isTRUE(f_x > f_x_previous))
}

.ojl_assess_state <- function(st, d, o, m) {
  if (m$name == "NewtonTrustRegion" && !(st$rho > st$eta)) return(c(x = FALSE, f = FALSE, g = FALSE, inc = FALSE))
  if (!is.null(m$assess)) return(m$assess(st, d, o, m))
  # assess_convergence(::ZerothOrderState): SimulatedAnnealing and
  # ParticleSwarm run until a limit
  if (isTRUE(m$zeroth)) return(c(x = FALSE, f = FALSE, g = FALSE, inc = FALSE))
  .ojl_assess(st$x, st$x_previous, st$f_x, st$f_x_previous, st$g_x, o$x_abstol, o$x_reltol, o$f_abstol,
              o$f_reltol, o$g_abstol)
}

.ojl_converged <- function(r) {
  sb <- r$stopped_by
  flags <- sb$x_converged || sb$f_converged || sb$g_converged
  x_fin <- is.finite(r$x_abschange) || is.nan(r$x_relchange)
  # Optim tests the Bool stopped_by.iterations here (> 0), not the count
  f_fin <- if (isTRUE(sb$iterations)) is.finite(r$f_abschange) || is.nan(r$f_relchange) else TRUE
  flags && x_fin && f_fin && is.finite(r$g_residual)
}

.ojl_retcode <- function(r) {
  if (isTRUE(r$samin_converge > 0L)) return("Success")
  if (.ojl_converged(r)) "Success" else if (isTRUE(r$stopped_by$iterations)) "MaxIters" else "Failure"
}

.ojl_changes <- function(x, x_previous, f_x, f_x_previous) {
  dx <- if (length(x)) max(abs(x - x_previous)) else 0
  list(x_abschange = dx, x_relchange = dx / max(abs(x)), f_abschange = abs(f_x - f_x_previous),
       f_relchange = abs(f_x - f_x_previous) / abs(f_x))
}

# optimize(d, initial_x, method, options, state)
.ojl_optimize <- function(d, x0, m, o, st = NULL) {
  if (is.null(st)) st <- .ojl_init(m, o, d, x0)
  t0 <- proc.time()[["elapsed"]]
  stopped_by_time_limit <- FALSE
  f_limit <- g_limit <- h_limit <- FALSE
  conv <- c(x = FALSE, f = FALSE, g = FALSE, inc = FALSE)
  counter_f_tol <- 0L
  small_tr <- FALSE
  zeroth <- isTRUE(m$zeroth)
  if (zeroth) {
    conv[["g"]] <- FALSE; stopped <- FALSE
  } else {
    conv[["g"]] <- isTRUE(.ojl_g_residual(st) <= o$g_abstol)
    stopped <- !is.finite(st$f_x) || any(!is.finite(st$g_x))
  }
  if (!is.null(m$initial_convergence)) {
    ic <- m$initial_convergence(st, o)
    conv[["g"]] <- ic[[1L]]; stopped <- ic[[2L]]
  }
  converged <- conv[["g"]] || stopped
  iteration <- 0L
  stopped_by_callback <- is.function(o$callback) && isTRUE(o$callback(st))
  stopped <- stopped || stopped_by_callback
  ls_success <- TRUE
  while (!converged && !stopped && iteration < o$iterations) {
    iteration <- iteration + 1L
    ls_success <- !.ojl_update(d, st, m)
    if (!ls_success) break
    .ojl_fgh(d, st, m)
    if (!.ojl_accept(d, st, m, o)) {
      if (isTRUE(o$show_warnings)) warning("Terminated early: trial iterate had non-finite values.", call. = FALSE)
      ls_success <- FALSE
      break
    }
    conv <- .ojl_assess_state(st, d, o, m)
    counter_f_tol <- if (conv[["f"]]) counter_f_tol + 1L else 0L
    converged <- conv[["x"]] || conv[["g"]] || (counter_f_tol > o$successive_f_tol)
    if (is.function(o$callback)) stopped_by_callback <- isTRUE(o$callback(st))
    stopped_by_time_limit <- isTRUE(proc.time()[["elapsed"]] - t0 > o$time_limit)
    f_limit <- o$f_calls_limit > 0 && d$obj$f_calls >= o$f_calls_limit
    g_limit <- o$g_calls_limit > 0 && d$obj$df_calls >= o$g_calls_limit
    h_limit <- o$h_calls_limit > 0 && d$obj$h_calls >= o$h_calls_limit
    if ((conv[["inc"]] && !o$allow_f_increases) || stopped_by_callback || stopped_by_time_limit || f_limit ||
        g_limit || h_limit) stopped <- TRUE
    if (m$name == "NewtonTrustRegion" && st$delta <= m$delta_min) {
      small_tr <- TRUE
      stopped <- TRUE
    }
    if (!is.null(st$g_x) && !all(is.finite(st$g_x))) {
      if (isTRUE(o$show_warnings)) warning("Terminated early due to NaN in gradient.", call. = FALSE)
      break
    }
    if (!is.null(st$H_x) && !all(is.finite(st$H_x))) {
      if (isTRUE(o$show_warnings)) warning("Terminated early due to NaN in Hessian.", call. = FALSE)
      break
    }
  }
  if (!is.null(m$after_while)) m$after_while(d, st, m, o)
  f_incr_pick <- conv[["inc"]] && !o$allow_f_increases
  best_x <- if (f_incr_pick && !isTRUE(m$keep_x)) st$x_previous else st$x
  best_f <- if (f_incr_pick && !isTRUE(m$keep_x)) st$f_x_previous else st$f_x
  ch <- if (zeroth) list(x_abschange = NaN, x_relchange = NaN, f_abschange = NaN, f_relchange = NaN)
        else .ojl_changes(st$x, st$x_previous, st$f_x, st$f_x_previous)
  gres <- if (!is.null(m$g_residual)) m$g_residual(st) else if (zeroth) NaN else .ojl_g_residual(st)
  c(list(minimizer = best_x, minimum = best_f, iterations = iteration, g_residual = gres,
         f_calls = d$obj$f_calls, g_calls = d$obj$df_calls, h_calls = d$obj$h_calls, state = st,
         stopped_by = list(x_converged = conv[["x"]], f_converged = conv[["f"]], g_converged = conv[["g"]],
                           f_limit_reached = f_limit, g_limit_reached = g_limit, h_limit_reached = h_limit,
                           time_limit = stopped_by_time_limit, callback = stopped_by_callback,
                           f_increased = f_incr_pick, ls_failed = !ls_success,
                           iterations = iteration == o$iterations, small_trustregion_radius = small_tr)),
    ch)
}

# ---------------------------------------------------------------------------
# Method dispatch: initial_state, update_state!, update_fgh!, accept_step!

.ojl_state <- function(...) list2env(list(...), parent = emptyenv())

.ojl_first_order_state <- function(d, x0, ..., M = NULL) {
  r <- d$value_gradient(x0)
  na <- rep(NaN, length(x0))
  st <- .ojl_state(x = x0, g_x = .ojl_project(M, r[[2L]], x0), f_x = r[[1L]], x_previous = na, f_x_previous = NaN,
                   s = na,
                   x_candidate = na, g_candidate = na, f_candidate = NaN, x_ls = na, alpha = 1)
  extra <- list(...)
  for (nm in names(extra)) st[[nm]] <- extra[[nm]]
  st
}

.ojl_bfgs_invH0 <- function(m, x, g) {
  n <- length(x)
  if (is.function(m$initial_invH)) return(m$initial_invH(x))
  if (is.null(m$initial_stepnorm)) return(diag(1, n))
  diag(m$initial_stepnorm * (1 / max(abs(g))), n)
}

.ojl_init <- function(m, o, d, x0) {
  x0 <- as.numeric(x0)
  n <- length(x0)
  na <- rep(NaN, n)
  M <- m$manifold
  if (m$name %in% c("GradientDescent", "AcceleratedGradientDescent", "MomentumGradientDescent", "BFGS", "LBFGS",
                    "ConjugateGradient")) x0 <- .ojl_retract(M, x0)
  switch(m$name,
    GradientDescent = , AcceleratedGradientDescent = , MomentumGradientDescent = {
      st <- .ojl_first_order_state(d, x0, M = M)
      if (m$name == "AcceleratedGradientDescent") {
        st$iteration <- 0L; st$y <- x0; st$y_previous <- na; st$y_candidate <- na
      }
      if (m$name == "MomentumGradientDescent") {
        st$x_previous <- x0; st$x_momentum <- na
      }
      st
    },
    BFGS = {
      st <- .ojl_first_order_state(d, x0, g_x_previous = na, dx = na, dg = na, u = na, M = M)
      st$invH <- .ojl_bfgs_invH0(m, x0, st$g_x)
      st
    },
    LBFGS = .ojl_first_order_state(d, x0, g_x_previous = na, rho = rep(NaN, m$m),
                                   dx_history = rep(list(na), m$m), dg_history = rep(list(na), m$m), dx = na,
                                   dg = na, u = na, twoloop_q = na, twoloop_alpha = rep(NaN, m$m),
                                   pseudo_iteration = 0L, M = M),
    ConjugateGradient = {
      st <- .ojl_first_order_state(d, x0, g_x_previous = na, y = 0 * x0, py = 0 * x0, beta = 0, M = M)
      .ojl_precondprep(m, x0)
      st$pg <- .ojl_precondition_apply(m$P, st$g_x)
      if (!is.null(m$P)) st$pg <- .ojl_project(M, st$pg, x0)
      st$s <- -st$pg
      st
    },
    Newton = {
      r <- d$value_gradient_hessian(x0)
      .ojl_state(x = x0, g_x = r[[2L]], H_x = r[[3L]], f_x = r[[1L]], x_previous = na, f_x_previous = NaN,
                 s = na, x_candidate = na, g_candidate = na, H_candidate = matrix(NaN, n, n), f_candidate = NaN,
                 x_ls = na, alpha = 1)
    },
    NewtonTrustRegion = {
      r <- d$value_gradient_hessian(x0)
      .ojl_state(x = x0, g_x = r[[2L]], H_x = r[[3L]], f_x = r[[1L]], x_previous = na, f_x_previous = NaN,
                 s = na, x_cache = na, g_cache = na, hard_case = FALSE, reached_subproblem_solution = TRUE,
                 interior = TRUE, delta = m$initial_delta, lambda = NaN, eta = m$eta, rho = 0)
    },
    NGMRES = , OACCEL = .ojl_ngmres_init(m, o, d, x0),
    if (is.function(m$init)) m$init(m, o, d, x0) else stop("unknown Optim method ", m$name, call. = FALSE))
}

.ojl_update <- function(d, st, m) {
  M <- m$manifold
  switch(m$name,
    GradientDescent = {
      .ojl_precondprep(m, st$x)
      st$s <- -.ojl_precondition_apply(m$P, st$g_x)
      if (!is.null(m$P)) st$s <- .ojl_project(M, st$s, st$x)
      ok <- .ojl_perform_linesearch(st, m, d)
      st$x_candidate <- .ojl_retract(M, st$x_ls)
      !ok
    },
    MomentumGradientDescent = {
      st$s <- -st$g_x
      st$x_momentum <- st$x - st$x_previous
      ok <- .ojl_perform_linesearch(st, m, d)
      st$x_candidate <- .ojl_retract(M, st$x + st$alpha * st$s + m$mu * st$x_momentum)
      !ok
    },
    AcceleratedGradientDescent = {
      st$s <- -st$g_x
      ok <- .ojl_perform_linesearch(st, m, d)
      st$y_candidate <- .ojl_retract(M, st$x_ls)
      k <- st$iteration + 1L
      scaling <- (k - 1) / (k + 2)
      st$x_candidate <- .ojl_retract(M, st$y_candidate + scaling * (st$y_candidate - st$y))
      !ok
    },
    BFGS = {
      st$s <- .ojl_project(M, -.ojl_gemv(st$invH, st$g_x), st$x)
      ok <- .ojl_perform_linesearch(st, m, d)
      st$x_candidate <- .ojl_retract(M, st$x_ls)
      st$dx <- st$x_candidate - st$x
      !ok
    },
    LBFGS = {
      st$pseudo_iteration <- st$pseudo_iteration + 1L
      .ojl_precondprep(m, st$x)
      .ojl_twoloop(st, m)
      st$s <- .ojl_project(M, st$s, st$x)
      ok <- .ojl_perform_linesearch(st, m, d)
      st$x_candidate <- .ojl_retract(M, st$x_ls)
      st$dx <- st$x_candidate - st$x
      !ok
    },
    ConjugateGradient = {
      ok <- .ojl_perform_linesearch(st, m, d)
      st$x_candidate <- .ojl_retract(M, st$x_ls)
      !ok
    },
    Newton = {
      L <- .ojl_positive_cholesky(st$H_x)
      st$s <- -backsolve(t(L), forwardsolve(L, st$g_x))
      ok <- .ojl_perform_linesearch(st, m, d)
      st$x_candidate <- st$x_ls
      !ok
    },
    NewtonTrustRegion = .ojl_ntr_update(d, st, m),
    NGMRES = , OACCEL = .ojl_ngmres_update(d, st, m),
    m$update(d, st, m))
}

.ojl_fgh <- function(d, st, m) {
  switch(m$name,
    GradientDescent = , MomentumGradientDescent = , AcceleratedGradientDescent = , BFGS = , LBFGS = ,
    ConjugateGradient = {
      r <- d$value_gradient(st$x_candidate)
      st$g_candidate <- .ojl_project(m$manifold, r[[2L]], st$x_candidate); st$f_candidate <- r[[1L]]
    },
    Newton = {
      r <- d$value_gradient_hessian(st$x_candidate)
      st$f_candidate <- r[[1L]]; st$g_candidate <- r[[2L]]; st$H_candidate <- r[[3L]]
    },
    NewtonTrustRegion = NULL,
    NGMRES = , OACCEL = .ojl_ngmres_fgh(d, st, m),
    if (is.function(m$fgh)) m$fgh(d, st, m))
  invisible(NULL)
}

.ojl_candidate_finite <- function(st)
  is.finite(st$f_candidate) && all(is.finite(st$g_candidate)) && all(is.finite(st$x_candidate))

.ojl_commit <- function(st) {
  st$x <- st$x_candidate; st$g_x <- st$g_candidate; st$f_x <- st$f_candidate
}

.ojl_accept <- function(d, st, m, o) {
  switch(m$name,
    GradientDescent = , MomentumGradientDescent = {
      if (!.ojl_candidate_finite(st)) return(FALSE)
      .ojl_commit(st)
      TRUE
    },
    AcceleratedGradientDescent = {
      if (!.ojl_candidate_finite(st) || !all(is.finite(st$y_candidate))) return(FALSE)
      st$y_previous <- st$y; st$y <- st$y_candidate
      .ojl_commit(st)
      st$iteration <- st$iteration + 1L
      TRUE
    },
    BFGS = {
      if (!.ojl_candidate_finite(st)) return(FALSE)
      .ojl_commit(st)
      .ojl_bfgs_update(st)
      TRUE
    },
    LBFGS = {
      if (!.ojl_candidate_finite(st)) {
        st$pseudo_iteration <- st$pseudo_iteration - 1L
        return(FALSE)
      }
      .ojl_commit(st)
      .ojl_lbfgs_store(st, m)
      TRUE
    },
    ConjugateGradient = {
      if (!.ojl_candidate_finite(st)) return(FALSE)
      .ojl_commit(st)
      .ojl_precondprep(m, st$x)
      dPd <- .ojl_inverse_precondition(m$P, st$s)
      etak <- m$eta * .ojl_dot(st$s, st$g_x_previous) / dPd
      st$y <- st$g_x - st$g_x_previous
      ydots <- .ojl_dot(st$y, st$s)
      pg_old <- st$pg
      st$pg <- .ojl_precondition_apply(m$P, st$g_x)
      st$py <- st$pg - pg_old
      betak <- (.ojl_dot(st$y, st$pg) - .ojl_dot(st$y, st$py) * .ojl_dot(st$g_x, st$s) / ydots * m$theta) / ydots
      beta <- .ojl_nanmax(betak, etak)
      beta <- .ojl_clamp(beta, -m$betamax, m$betamax)
      st$beta <- beta
      st$s <- .ojl_project(m$manifold, beta * st$s - st$pg, st$x)
      TRUE
    },
    Newton = {
      if (!.ojl_candidate_finite(st) || !all(is.finite(st$H_candidate))) return(FALSE)
      .ojl_commit(st)
      st$H_x <- st$H_candidate
      TRUE
    },
    if (is.function(m$accept)) m$accept(d, st, m, o) else TRUE)
}

.ojl_bfgs_update <- function(st) {
  st$dg <- st$g_x - st$g_x_previous
  dx <- st$dx
  dx_dg <- .ojl_dot(dx, st$dg)
  if (isTRUE(dx_dg > 0)) {
    u <- .ojl_gemv(st$invH, st$dg)
    st$u <- u
    c1 <- (dx_dg + .ojl_dot(st$dg, u)) / (dx_dg * dx_dg)
    c2 <- 1 / dx_dg
    n <- length(dx)
    inner <- .ojl_fma(rep(c2 * u, each = n), rep(-dx, n), as.numeric(st$invH))
    mid <- .ojl_fma(rep(-u, n), rep(c2 * dx, each = n), inner)
    st$invH <- matrix(.ojl_fma(rep(dx, n), rep(c1 * dx, each = n), mid), n, n)
  }
  invisible(NULL)
}

.ojl_lbfgs_store <- function(st, m) {
  st$dg <- st$g_x - st$g_x_previous
  rho_it <- 1 / .ojl_dot(st$dx, st$dg)
  if (is.infinite(rho_it)) {
    st$pseudo_iteration <- 0L
    return(invisible(NULL))
  }
  idx <- .ojl_mod1(st$pseudo_iteration, m$m)
  st$dx_history[[idx]] <- st$dx
  st$dg_history[[idx]] <- st$dg
  st$rho[idx] <- rho_it
  invisible(NULL)
}

# twoloop! of l_bfgs.jl
.ojl_twoloop <- function(st, m) {
  M <- m$m
  pit <- st$pseudo_iteration
  lower <- pit - M; upper <- pit - 1L
  q <- st$g_x
  alpha <- st$twoloop_alpha
  for (index in seq.int(upper, lower, by = -1L)) {
    if (index < 1L) next
    i <- .ojl_mod1(index, M)
    alpha[i] <- st$rho[i] * .ojl_dot(st$dx_history[[i]], q)
    q <- q - alpha[i] * st$dg_history[[i]]
  }
  if (isTRUE(m$scaleinvH0) && pit > 1L) {
    i <- .ojl_mod1(upper, M)
    dgi <- st$dg_history[[i]]
    scaling <- .ojl_dot(st$dx_history[[i]], dgi) / sum(dgi * dgi)
    s <- scaling * q
  } else {
    s <- .ojl_precondition_apply(m$P, q)
  }
  for (index in seq.int(lower, upper)) {
    if (index < 1L) next
    i <- .ojl_mod1(index, M)
    beta <- st$rho[i] * .ojl_dot(st$dg_history[[i]], s)
    s <- s + st$dx_history[[i]] * (alpha[i] - beta)
  }
  st$twoloop_alpha <- alpha
  st$twoloop_q <- q
  st$s <- -s
  invisible(NULL)
}

# reset! of the methods that Fminbox restarts
.ojl_reset <- function(m, st, d, x) {
  if (!(m$name %in% c("GradientDescent", "ConjugateGradient", "BFGS", "LBFGS")) && !is.function(m$reset))
    stop("Fminbox cannot restart ", m$name, "(): Optim defines no reset! for it.", call. = FALSE)
  if (is.function(m$reset)) return(m$reset(m, st, d, x))
  st$x <- .ojl_retract(m$manifold, x)
  r <- d$value_gradient(st$x)
  st$g_x <- .ojl_project(m$manifold, r[[2L]], st$x); st$f_x <- r[[1L]]
  na <- rep(NaN, length(x))
  st$x_previous <- na; st$f_x_previous <- NaN
  switch(m$name,
    GradientDescent = st$s <- na,
    ConjugateGradient = {
      st$g_x_previous <- na
      .ojl_precondprep(m, st$x)
      st$pg <- .ojl_precondition_apply(m$P, st$g_x)
      if (!is.null(m$P)) st$pg <- .ojl_project(m$manifold, st$pg, st$x)
      st$s <- -st$pg
    },
    BFGS = {
      st$g_x_previous <- na
      st$invH <- .ojl_bfgs_invH0(m, x, r[[2L]])
    },
    LBFGS = {
      st$g_x_previous <- na
      st$pseudo_iteration <- 0L
    })
  invisible(NULL)
}

# cholesky!(Positive, H) of PositiveFactorizations.jl: the lower factor of a
# positive definite modification of H (the unblocked algorithm, which is
# the blocked one for up to 32 unknowns)
.ojl_positive_cholesky <- function(A) {
  K <- nrow(A)
  if (!K) return(A)
  tol <- 10 * K * .Machine$double.eps * max(abs(A))
  B <- A
  for (j in seq_len(K)) {
    Bjj <- B[j, j]
    if (abs(Bjj) > tol) {
      dj <- sign(Bjj)
      s <- sqrt(abs(Bjj))
      B[j, j] <- s
      if (j < K) {
        idx <- (j + 1L):K
        B[idx, j] <- B[idx, j] * (dj / s)
        cc <- B[idx, j]
        for (k in idx) {
          ii <- k:K
          B[ii, k] <- B[ii, k] + cc[ii - j] * (-dj * cc[k - j])
        }
      }
    } else {
      B[j, j] <- 1
      if (j < K) B[(j + 1L):K, j] <- 0
    }
  }
  B[upper.tri(B)] <- 0
  B
}

# ---------------------------------------------------------------------------
# NewtonTrustRegion (second_order/newton_trust_region.jl)

.ojl_tr_subproblem <- function(gr, H, delta, max_iters = 5L) {
  n <- length(gr)
  eps <- .Machine$double.eps
  delta_sq <- delta * delta
  out <- function(m, interior, lambda, hard_case, reached, s)
    list(m = m, interior = interior, lambda = lambda, hard_case = hard_case, reached = reached, s = s)
  if (any(!is.finite(H))) return(out(Inf, FALSE, 0, FALSE, FALSE, rep(0, n)))
  e <- eigen(H, symmetric = TRUE)
  vals <- rev(e$values)
  vecs <- e$vectors[, n:1L, drop = FALSE]
  min_ev <- vals[1L]; max_ev <- vals[n]
  H_scale <- max(abs(min_ev), abs(max_ev))
  gr_norm <- sqrt(sum(gr * gr))
  if (H_scale == 0) {
    if (gr_norm == 0) return(out(0, TRUE, 0, FALSE, TRUE, rep(0, n)))
    return(out(-delta * gr_norm, FALSE, gr_norm / delta, FALSE, TRUE, (-delta / gr_norm) * gr))
  }
  scale_tol <- sqrt(eps) * H_scale
  gr_tol <- sqrt(eps) * gr_norm
  qg <- as.numeric(crossprod(vecs, gr))
  calc_p <- function(lambda, min_i) {
    p <- rep(0, n)
    if (min_i <= n) for (i in min_i:n) p <- p + (-qg[i] / (vals[i] + lambda)) * vecs[, i]
    p
  }
  s <- rep(0, n)
  hard_case <- FALSE
  reached <- TRUE
  lambda <- 0
  pd <- min_ev > scale_tol
  if (pd) s <- calc_p(0, 1L)
  if (pd && sum(s * s) <= delta_sq) {
    interior <- TRUE
  } else {
    interior <- FALSE
    qg_tol <- max(gr_tol, scale_tol * delta)
    hc <- FALSE; min_i <- 1L
    if (vals[1L] < 0) {
      hc <- TRUE; li <- 1L
      repeat {
        if (li > n) break
        if (abs(vals[1L] - vals[li]) > scale_tol) break
        if (abs(qg[li]) > qg_tol) {
          hc <- FALSE
          break
        }
        li <- li + 1L
      }
      min_i <- li
    }
    lambda_lb <- max(0, .ojl_nextfloat(-min_ev))
    lambda_ub <- gr_norm / delta - min_ev
    lambda <- lambda_lb
    lambda_tol <- max((eps^(1 / 3))^2 * (lambda_ub - lambda_lb), 4 * eps * lambda_ub)
    if (hc) {
      s <- calc_p(lambda, min_i)
      p2 <- sum(s * s)
      if (!(p2 > delta_sq)) {
        hard_case <- TRUE
        reached <- TRUE
        tau <- sqrt(delta_sq - p2)
        s <- s + tau * vecs[, 1L]
      }
    }
    norm2_lb <- 0
    first_nz <- 1L
    boundary <- FALSE
    for (i in seq_len(n)) {
      dd <- vals[i] + lambda_lb
      if (dd <= scale_tol) {
        first_nz <- i + 1L
        if (abs(qg[i]) > gr_tol) {
          boundary <- TRUE
          break
        }
      } else {
        norm2_lb <- norm2_lb + (qg[i] / dd)^2
      }
    }
    interior_at_lb <- !boundary && norm2_lb <= delta_sq
    if (!hard_case && interior_at_lb) {
      s <- calc_p(lambda_lb, first_nz)
      interior <- TRUE
      reached <- TRUE
      lambda <- lambda_lb
    } else if (!hard_case) {
      lambda <- .ojl_initial_safeguards(H, gr, delta, lambda)
      reached <- FALSE
      for (iter in seq_len(max_iters)) {
        lambda_previous <- lambda
        Hr <- H
        diag(Hr) <- diag(H) + lambda
        R <- tryCatch(chol(Hr), error = function(e) NULL)
        if (is.null(R)) {
          lambda <- max(2 * lambda, sqrt(eps) * H_scale)
          next
        }
        s <- -backsolve(R, forwardsolve(t(R), gr))
        q_l <- forwardsolve(t(R), s)
        norm2_s <- sum(s * s)
        lambda_update <- norm2_s * (sqrt(norm2_s) - delta) / (delta * sum(q_l * q_l))
        lambda <- lambda + lambda_update
        if (lambda < lambda_lb) lambda <- (lambda_previous + lambda_lb) / 2
        else if (lambda > lambda_ub) lambda <- (lambda_previous + lambda_ub) / 2
        if (abs(lambda - lambda_previous) < lambda_tol) {
          reached <- abs(sqrt(norm2_s) - delta) <= eps^(1 / 3) * delta
          break
        }
      }
    }
  }
  m <- sum(gr * s) + sum(colSums(s * H) * s) / 2
  out(m, interior, lambda, hard_case, reached, s)
}

.ojl_initial_safeguards <- function(H, gr, delta, lambda) {
  lS <- -min(diag(H))
  gn <- sqrt(sum(gr * gr))
  Hn <- max(colSums(abs(H)))
  lL <- max(0, lS, gn / delta - Hn)
  lU <- gn / delta + Hn
  lambda <- .ojl_clamp(lambda, lL, lU)
  if (lambda <= lS) lambda <- max(1 / 1000 * lU, sqrt(lL * lU))
  lambda
}

.ojl_ntr_update <- function(d, st, m) {
  r <- .ojl_tr_subproblem(st$g_x, st$H_x, st$delta)
  st$s <- r$s
  st$interior <- r$interior; st$lambda <- r$lambda; st$hard_case <- r$hard_case
  st$reached_subproblem_solution <- r$reached
  mm <- r$m
  st$x_cache <- st$x
  f_cache <- st$f_x
  st$x <- st$x + st$s
  if (m$use_fg) {
    st$g_cache <- st$g_x
    v <- d$value_gradient(st$x)
    st$g_x <- v[[2L]]; st$f_x <- v[[1L]]
  } else {
    st$f_x <- d$value_bang(st$x)
  }
  f_x_diff <- f_cache - st$f_x
  st$rho <- if (abs(mm) <= .Machine$double.eps) 1 else if (mm > 0) -1 else f_x_diff / (-mm)
  accept <- isTRUE(st$rho > st$eta)
  if (!accept) {
    s_norm <- sqrt(sum(st$s * st$s))
    st$delta <- (if (is.finite(s_norm)) min(st$delta, s_norm) else st$delta) / 4
  } else if (st$rho < m$rho_lower) {
    st$delta <- st$delta / 4
  } else if (st$rho > m$rho_upper && !st$interior) {
    st$delta <- min(2 * st$delta, m$delta_hat)
  }
  if (accept) {
    if (m$use_fg) {
      st$H_x <- d$hessian_bang(st$x)
    } else {
      v <- d$gradient_hessian(st$x)
      st$g_x <- v[[1L]]; st$H_x <- v[[2L]]
    }
    st$x_previous <- st$x_cache
    st$f_x_previous <- f_cache
  } else {
    st$x <- st$x_cache
    st$f_x <- f_cache
    if (m$use_fg) st$g_x <- st$g_cache
  }
  FALSE
}

# ---------------------------------------------------------------------------
# NGMRES and OACCEL (first_order/ngmres.jl): the state shares x, g_x and
# x_previous with the state of the nonlinear preconditioner

.ojl_ngmres_init <- function(m, o, d, x0) {
  inner <- .ojl_init(m$nlprecon, m$nlpreconopts, d, x0)
  n <- length(x0); w <- m$wmax
  st <- new.env(parent = emptyenv())
  for (nm in c("x", "g_x", "x_previous")) local({
    field <- nm
    makeActiveBinding(field, function(v) if (missing(v)) inner[[field]] else assign(field, v, envir = inner), st)
  })
  st$inner <- inner
  st$f_x <- inner$f_x
  st$x_previous_0 <- rep(NaN, n)
  st$f_x_previous <- NaN; st$f_x_previous_0 <- NaN; st$f_xP <- NaN; st$grnorm_xP <- NaN
  st$s <- rep(NaN, n)
  st$X <- matrix(NaN, n, w); st$R <- matrix(NaN, n, w); st$Q <- matrix(NaN, w, w)
  st$xi <- if (m$name == "OACCEL") matrix(NaN, w, 2L) else rep(NaN, w)
  st$curw <- 1L
  st$A <- matrix(NaN, w, w); st$b <- rep(NaN, w)
  st$k <- 0L; st$restart <- FALSE; st$g_abstol <- o$g_abstol
  st$subspacealpha <- rep(NaN, w)
  st$x_ls <- rep(NaN, n); st$alpha <- 1
  st$X[, 1L] <- inner$x
  st$R[, 1L] <- inner$g_x
  .ojl_ngmres_Q(st, 1L, 1L, m)
  st
}

.ojl_ngmres_Q <- function(st, i, j, m) {
  if (m$name == "NGMRES") {
    st$Q[j, i] <- .ojl_dot(st$R[, j], st$R[, i])
    if (i != j) st$Q[i, j] <- st$Q[j, i]
  } else {
    st$Q[i, j] <- .ojl_dot(st$X[, i], st$R[, j])
    if (i != j) st$Q[j, i] <- .ojl_dot(st$X[, j], st$R[, i])
  }
}

.ojl_refresh_precon <- function(d, st, p) {
  inner <- st$inner
  r <- d$value_gradient(inner$x)
  inner$f_x <- r[[1L]]
  inner$g_x <- .ojl_project(p$manifold, r[[2L]], inner$x)
}

.ojl_ngmres_post_accelerate <- function(d, st, p) {
  inner <- st$inner
  if (p$name == "LBFGS") {
    inner$pseudo_iteration <- inner$pseudo_iteration + 1L
    .ojl_refresh_precon(d, st, p)
    .ojl_lbfgs_store(inner, p)
  } else if (p$name == "BFGS") {
    .ojl_refresh_precon(d, st, p)
    .ojl_bfgs_update(inner)
  } else {
    .ojl_refresh_precon(d, st, p)
  }
}

.ojl_ngmres_update <- function(d, st, m) {
  st$x_previous_0 <- st$x
  st$f_x_previous_0 <- st$f_x
  st$k <- st$k + 1L
  curw <- st$curw
  .ojl_optimize(d, st$x, m$nlprecon, m$nlpreconopts, st$inner)
  st$g_x <- .ojl_project(m$manifold, st$g_x, st$x)
  st$f_xP <- st$inner$f_x
  st$f_x <- st$inner$f_x
  if (!all(is.finite(st$x)) || !all(is.finite(st$g_x)) || !is.finite(st$f_xP)) {
    warning("Non-finite values attained from preconditioner ", m$nlprecon$name, ".", call. = FALSE)
    return(TRUE)
  }
  st$grnorm_xP <- max(abs(st$g_x))
  if (st$grnorm_xP <= st$g_abstol) return(FALSE)
  .ojl_refresh_precon(d, st, m$nlprecon)
  x <- st$x; g <- st$g_x
  oaccel <- m$name == "OACCEL"
  eta <- if (oaccel) .ojl_dot(x, g) else .ojl_dot(g, g)
  for (i in seq_len(curw)) {
    if (oaccel) {
      st$xi[i, 1L] <- .ojl_dot(st$X[, i], g)
      st$xi[i, 2L] <- .ojl_dot(x, st$R[, i])
      st$b[i] <- eta - st$xi[i, 1L]
    } else {
      st$xi[i] <- .ojl_dot(g, st$R[, i])
      st$b[i] <- eta - st$xi[i]
    }
  }
  for (i in seq_len(curw)) for (j in seq_len(curw)) {
    st$A[i, j] <- if (oaccel) st$Q[i, j] - st$xi[i, 1L] - st$xi[j, 2L] + eta
                  else st$Q[i, j] - st$xi[i] - st$xi[j] + eta
  }
  Av <- st$A[seq_len(curw), seq_len(curw), drop = FALSE]
  bv <- st$b[seq_len(curw)]
  dl <- m$epsilon0 * max(max(diag(Av)), m$epsilon0)
  alpha <- tryCatch(.ojl_backslash(Av + dl * diag(1, curw), bv), error = function(e) {
    warning("Calculating alpha failed in ", m$name, ".", call. = FALSE)
    rep(NaN, curw)
  })
  st$subspacealpha[seq_len(curw)] <- alpha
  if (any(is.nan(alpha))) {
    warning("Calculated alpha is NaN in ", m$name, ". Restarting ...", call. = FALSE)
    st$s <- rep(0, length(x))
    st$restart <- TRUE
  } else {
    acc <- st$X[, 1L] * alpha[1L]
    if (curw > 1L) for (k in 2:curw) acc <- acc + st$X[, k] * alpha[k]
    xA <- (1 - sum(alpha)) * x + acc
    st$s <- xA - x
  }
  sTg <- .ojl_dot(st$s, st$g_x)
  if (sTg >= 0 || !is.finite(sTg)) {
    st$restart <- TRUE
    ok <- TRUE
    st$alpha <- 0
  } else {
    st$restart <- FALSE
    st$f_x_previous <- st$inner$f_x_previous
    if (m$alphaguess$type == "InitialConstantChange" && identical(m$nlprecon$alphaguess$type, "InitialConstantChange"))
      m$alphaguess$prev$dphi <- m$nlprecon$alphaguess$prev$dphi
    ok <- .ojl_perform_linesearch(st, m, d)
    if (ok) {
      st$x <- .ojl_retract(m$manifold, st$x_ls)
      st$inner$f_x_previous <- st$f_x_previous
      if (m$alphaguess$type == "InitialConstantChange" &&
          identical(m$nlprecon$alphaguess$type, "InitialConstantChange"))
        m$nlprecon$alphaguess$prev$dphi <- m$alphaguess$prev$dphi
      .ojl_ngmres_post_accelerate(d, st, m$nlprecon)
    } else {
      st$restart <- TRUE
      st$alpha <- 0
      ok <- TRUE
    }
  }
  st$x_previous <- st$x_previous_0
  st$f_x_previous <- st$f_x_previous_0
  !ok
}

.ojl_ngmres_fgh <- function(d, st, m) {
  r <- d$value_gradient(st$x)
  st$g_x <- .ojl_project(m$manifold, r[[2L]], st$x)
  st$f_x <- r[[1L]]
  if (st$restart) {
    st$k <- 0L
    st$curw <- 1L
  } else {
    st$curw <- min(st$curw + 1L, m$wmax)
  }
  j <- (st$k %% m$wmax) + 1L
  st$X[, j] <- st$x
  st$R[, j] <- r[[2L]]
  for (i in seq_len(st$curw)) .ojl_ngmres_Q(st, i, j, m)
}

# ---------------------------------------------------------------------------
# Fminbox (constrained/fminbox.jl)

.ojl_barrier_value <- function(x, l, u) {
  v <- 0
  for (i in seq_along(x)) {
    dxl <- x[i] - l[i]; dxu <- u[i] - x[i]
    t <- if (dxl <= 0 || dxu <= 0) Inf else {
      (if (is.finite(dxl)) -log(dxl) else 0) + (if (is.finite(dxu)) -log(dxu) else 0)
    }
    v <- if (i == 1L) t else v + t
  }
  v
}

.ojl_barrier_gradient <- function(x, l, u, transform = identity) {
  g <- numeric(length(x))
  for (i in seq_along(x)) {
    gi <- 0
    if (is.finite(l[i])) gi <- gi + transform(-1 / (x[i] - l[i]))
    if (is.finite(u[i])) gi <- gi + transform(1 / (u[i] - x[i]))
    g[i] <- gi
  }
  g
}

.ojl_barrier_wrapper <- function(obj, mu, l, u) {
  b <- new.env(parent = emptyenv())
  b$obj <- obj; b$mu <- mu
  in_box <- function(x) all(l <= x & x <= u)
  b$value_gradient <- function(x) {
    DFb <- .ojl_barrier_gradient(x, l, u)
    Fb <- .ojl_barrier_value(x, l, u)
    if (in_box(x)) {
      r <- obj$value_gradient(x)
      list(.ojl_fma(mu, Fb, r[[1L]]), .ojl_fma(mu, DFb, r[[2L]]))
    } else list(mu * Fb, mu * DFb)
  }
  b$value_bang <- function(x) {
    Fb <- .ojl_barrier_value(x, l, u)
    if (in_box(x)) .ojl_fma(mu, Fb, obj$value_bang(x)) else mu * Fb
  }
  b$value <- function(x) {
    Fb <- .ojl_barrier_value(x, l, u)
    if (in_box(x)) .ojl_fma(mu, Fb, obj$value(x)) else mu * Fb
  }
  b$gradient_bang <- function(x) {
    DFb <- .ojl_barrier_gradient(x, l, u)
    if (in_box(x)) .ojl_fma(mu, DFb, obj$gradient_bang(x)) else mu * DFb
  }
  b$value_jvp <- function(x, v) {
    JVPb <- sum(.ojl_barrier_gradient(x, l, u) * v)
    Fb <- .ojl_barrier_value(x, l, u)
    if (in_box(x)) {
      r <- obj$value_jvp(x, v)
      list(.ojl_fma(mu, Fb, r[[1L]]), .ojl_fma(mu, JVPb, r[[2L]]))
    } else list(mu * Fb, mu * JVPb)
  }
  b
}

.ojl_barrier_method <- function(m, P, precondprep) {
  if (m$name %in% c("ConjugateGradient", "LBFGS", "GradientDescent")) {
    mm <- m
    mm$P <- P
    mm$precondprep <- precondprep
    if (m$name == "LBFGS") {
      mm$m <- 10L
      mm$scaleinvH0 <- FALSE
    }
    if (m$name == "ConjugateGradient") {
      mm$theta <- 1
      mm$betamax <- Inf
    }
    return(mm)
  }
  if (m$name %in% c("NelderMead", "SimulatedAnnealing", "ParticleSwarm", "BFGS", "NGMRES", "OACCEL")) return(m)
  stop("You need to specify a valid inner optimizer for Fminbox, ", m$name,
       "() is not supported. Please consult the documentation.", call. = FALSE)
}

.ojl_fminbox <- function(df, l, u, x0, Fm, o) {
  t0 <- proc.time()[["elapsed"]]
  x <- as.numeric(x0)
  n <- length(x)
  P <- new.env(parent = emptyenv())
  P$diag <- x
  class(P) <- "InverseDiagonal"
  moved <- integer()
  for (i in seq_len(n)) {
    if (x[i] == l[i]) {
      x[i] <- 99 / 100 * l[i] + 1 / 100 * u[i]
      moved <- c(moved, i)
    } else if (x[i] == u[i]) {
      x[i] <- 1 / 100 * l[i] + 99 / 100 * u[i]
      moved <- c(moved, i)
    } else if (x[i] < l[i] || x[i] > u[i]) {
      stop(sprintf("Initial x[%d]=%s is outside of [%s, %s]", i, format(x[i]), format(l[i]), format(u[i])),
           call. = FALSE)
    }
  }
  if (length(moved))
    warning("Initial position cannot be on the boundary of the box. Moving elements to the interior.\n",
            "Element indices affected: [", paste(moved, collapse = ", "), "]", call. = FALSE)
  dfbox <- .ojl_barrier_wrapper(df, 0, l, u)
  prep <- function(P, x) {
    P$diag <- 1 / (dfbox$mu * (1 / (x - l)^2 + 1 / (u - x)^2) + 1)
    P
  }
  inner <- .ojl_method_env(.ojl_barrier_method(Fm$method, P, prep))
  st <- .ojl_init(inner, o, df, x)
  f_x <- st$f_x
  g_x <- if (!is.null(st$g_x)) st$g_x else df$gradient_bang(x)
  mu <- if (is.nan(Fm$mu0)) {
    gb <- sum(.ojl_barrier_gradient(x, l, u, abs))
    if (gb > 0) Fm$mufactor * sum(abs(g_x)) / gb else 0
  } else Fm$mu0
  if (!is.finite(mu))
    stop("Could not compute finite initial mu for Fminbox, got: ", Fm$mu0, ". Please provide a finite mu0 value.",
         call. = FALSE)
  f_calls <- 0L; g_calls <- 0L; h_calls <- NULL
  iteration <- 0L
  last <- NULL
  x_previous <- x
  f_x_previous <- f_x
  run_inner <- function() {
    dfbox <<- .ojl_barrier_wrapper(df, mu, l, u)
    .ojl_reset(inner, st, dfbox, x)
    res <- .ojl_optimize(dfbox, x, inner, o, st)
    f_calls <<- f_calls + df$f_calls
    g_calls <<- g_calls + df$df_calls
    if (is.null(h_calls)) h_calls <<- df$h_calls
    df$f_calls <- 0L; df$df_calls <- 0L; df$h_calls <- 0L
    res
  }
  converged <- FALSE
  repeat {
    iteration <- iteration + 1L
    if (iteration > 1L) mu <- mu * Fm$mufactor
    x_previous <- x
    f_x_previous <- f_x
    last <- run_inner()
    x <- last$minimizer
    r <- df$value_gradient(x)
    f_x <- r[[1L]]; g_x <- r[[2L]]
    g <- x - pmin(pmax(x - g_x, l), u)
    conv <- .ojl_assess(x, x_previous, f_x, f_x_previous, g, o$outer_x_abstol, o$outer_x_reltol,
                        o$outer_f_abstol, o$outer_f_reltol, o$outer_g_abstol)
    converged <- conv[["x"]] || conv[["f"]] || conv[["g"]] || last$stopped_by$callback
    stopped <- isTRUE(proc.time()[["elapsed"]] - t0 > o$time_limit)
    if (conv[["inc"]] && !o$allow_outer_f_increases) {
      warning("f(x) increased: stopping optimization", call. = FALSE)
      break
    }
    if (converged || stopped || iteration >= o$outer_iterations) break
  }
  ch <- .ojl_changes(x, x_previous, f_x, f_x_previous)
  c(list(minimizer = x, minimum = f_x, iterations = iteration, g_residual = max(abs(g)), f_calls = f_calls,
         g_calls = g_calls, h_calls = h_calls %||% 0L, state = st,
         stopped_by = list(x_converged = conv[["x"]], f_converged = conv[["f"]], g_converged = conv[["g"]],
                           callback = last$stopped_by$callback, time_limit = stopped,
                           f_increased = conv[["inc"]] && !o$allow_f_increases, ls_failed = FALSE,
                           iterations = last$stopped_by$iterations)),
    ch)
}

# ---------------------------------------------------------------------------
# Method objects and the OptimizationOptimJL solver

.ojl_alg <- function(name, ..., grad = TRUE, hess = FALSE) {
  a <- .opt_alg(name, ...)
  class(a) <- c(name, "OptimAlgorithm", "OptimizationAlgorithm")
  a$solver <- .ojl_solve
  a$traits <- list(grad = grad, hess = hess, bounds = name != "SimulatedAnnealing",
                   rbounds = name %in% c("Fminbox", "SAMIN"), cons = name == "IPNewton",
                   consjac = name == "IPNewton", conshess = name == "IPNewton")
  a
}

# a mutable copy of a method (its preconditioner and its nested methods)
.ojl_method_env <- function(alg) {
  if (is.environment(alg)) return(alg)
  m <- list2env(unclass(alg), parent = emptyenv())
  if (!is.null(m$nlprecon)) m$nlprecon <- .ojl_method_env(m$nlprecon)
  if (!is.null(m$method)) m$method <- unclass(m$method)
  m
}

#' Optim.jl optimizers
#'
#' The optimizers of Optim.jl, through OptimizationOptimJL: first order
#' methods (`GradientDescent()`, `ConjugateGradient()`,
#' `AcceleratedGradientDescent()`, `MomentumGradientDescent()`, `NGMRES()`,
#' `OACCEL()`, and Optim's `BFGS` and `LBFGS` as `Optim$BFGS()` and
#' `Optim$LBFGS()`), second order methods (`Newton()`,
#' `NewtonTrustRegion()`), the box constrained `Fminbox()` (which wraps a
#' first order method in a logarithmic barrier) and the interior point
#' `IPNewton()` for bounds and nonlinear constraints, and the derivative free
#' `ParticleSwarm()`, `SimulatedAnnealing()`, `SAMIN()` and Optim's Nelder
#' Mead simplex `Optim$NelderMead()`. Optim's own `Optim$Adam()` and
#' `Optim$AdaMax()` (which OptimizationOptimJL does not export, `Adam` being
#' the rule of Optimisers) and the Hessian vector product trust region
#' `Optim$KrylovTrustRegion()` (it uses the `hv` of the function, `(x, v, p)`,
#' when there is one, else the Hessian times the vector) are solved too.
#' The first order methods keep their iterates on a manifold:
#' `Optim$Flat()` (the default), `Optim$Sphere()`, `Optim$Stiefel()`,
#' `Optim$PowerManifold()` or `Optim$ProductManifold()`.
#' `Optim` is the module: a list of these constructors and of
#' `Optim$LineSearches`.
#'
#' `solve(prob, alg, maxiters =, abstol =, reltol =, maxtime =)` maps
#' `maxiters` to Optim's `iterations` (`outer_iterations` for `Fminbox`,
#' with `local_maxiters` for its inner iterations), `abstol` to `f_abstol`,
#' `reltol` to `f_reltol` and `maxtime` to `time_limit`; Optim's other
#' options (`g_abstol`, `x_abstol`, `x_reltol`, `allow_f_increases`,
#' `successive_f_tol`, `outer_g_abstol`, ...) pass through. A problem with
#' bounds and a method that does not handle them is wrapped in `Fminbox()`
#' (bounds are given to `ParticleSwarm()`). Without a gradient function a
#' central difference gradient is used, and without a Hessian a difference
#' of the gradient. The return code is `"Success"` when a convergence test
#' holds, `"MaxIters"` at the iteration limit and `"Failure"` otherwise.
#'
#' @param alphaguess The initial step guess of the line search (see
#'   [optim_linesearches]) or a number, a static first step.
#' @param linesearch The line search (for `IPNewton()`, a function
#'   `(state, constraints, ...)`; the default is its backtracking search on
#'   the barrier merit function).
#' @param P A preconditioner: `NULL`, a vector (a diagonal matrix) or a
#'   matrix; the search direction solves `P s = g`.
#' @param precondprep A function `(P, x)` returning the preconditioner at
#'   `x`.
#' @param manifold `NULL` or `Optim$Flat()` (Euclidean space), or the
#'   manifold the iterates are kept on: `Optim$Sphere()` (unit norm),
#'   `Optim$Stiefel()` (orthonormal columns, of the shape of `u0`),
#'   `Optim$PowerManifold()` or `Optim$ProductManifold()`.
#' @param eta,theta,betamax The parameters of the Hager Zhang conjugate
#'   gradient update; for `NewtonTrustRegion()`, `eta` is the threshold of
#'   the ratio of actual and predicted reductions for accepting a step.
#' @param mu The momentum of `MomentumGradientDescent()`.
#' @param initial_delta,delta_hat,delta_min The initial, largest and
#'   smallest trust region radii.
#' @param rho_lower,rho_upper The ratios below which the trust region
#'   shrinks and above which it grows.
#' @param use_fg Evaluate the gradient with the value after each step.
#' @param nlprecon,nlpreconopts The nonlinear preconditioner of `NGMRES()`
#'   and `OACCEL()` and its options (a list of Optim options).
#' @param epsilon0 The regularization of the subspace system.
#' @param wmax The window size.
#' @param method The inner method of `Fminbox()`.
#' @param mu0,mufactor The initial barrier multiplier of `Fminbox()` (`NaN`
#'   computes one) or of `IPNewton()` (`"auto"` computes one, or a number)
#'   and the decrease factor of `Fminbox()`.
#' @param lower,upper The bounds of the particles.
#' @param n_particles The number of particles (0 chooses the default).
#' @param neighbor A function of the current point returning the proposed
#'   point (the default adds standard normal noise).
#' @param temperature A function of the iteration giving the temperature
#'   (the default is `1 / log(t)`).
#' @param keep_best Return the best point visited.
#' @param nt,ns,t0,rt,r_expand,bounds_ratio,neps,coverage_ok,verbosity The
#'   parameters of SAMIN: the temperature is reduced every `nt * ns * n`
#'   evaluations and the search bounds adjusted every `ns * n`; the initial
#'   temperature, its reduction factor, the factor that raises the
#'   temperature while the parameter space is not covered, the acceptance
#'   ratio above which the search bounds grow (below one minus it they
#'   shrink), the number of previous values compared in the convergence
#'   test, whether the parameter space counts as covered from the start and
#'   the printing level (0 to 3).
#' @param show_linesearch Print the line search of `IPNewton()`.
#' @section The Optim module:
#' Constructors reached through `Optim$` only (their names are taken by the
#' `stats::optim` wrappers or by the rules of Optimisers):
#' \describe{
#'   \item{`Optim$BFGS(alphaguess, linesearch, initial_invH, initial_stepnorm, manifold)`}{
#'     `initial_invH` is a function of `x` giving the initial inverse
#'     Hessian and `initial_stepnorm` the norm of the first step.}
#'   \item{`Optim$LBFGS(m = 10, alphaguess, linesearch, P, precondprep, manifold, scaleinvH0)`}{
#'     `m` is the memory; `scaleinvH0` scales the initial inverse Hessian.}
#'   \item{`Optim$NelderMead(parameters, initial_simplex)`}{`parameters` is a
#'     list with `alpha`, `beta`, `gamma` and `delta` (`NULL` for the
#'     adaptive ones) and `initial_simplex` the initial simplex (`NULL` for
#'     the affine simplex).}
#'   \item{`Optim$Adam(alpha = 1e-4, beta_mean = 0.9, beta_var = 0.999, epsilon = 1e-8, manifold)`,
#'     `Optim$AdaMax(alpha = 0.002, beta_mean, beta_var, epsilon, manifold)`}{
#'     `alpha` is the step or a function of the iteration giving it,
#'     `beta_mean` and `beta_var` the decay rates of the moment estimates and
#'     `epsilon` keeps the denominators positive.}
#'   \item{`Optim$KrylovTrustRegion(initial_radius = 1, max_radius = 100, eta = 0.1, rho_lower = 0.25, rho_upper = 0.75, cg_tol = 0.01)`}{
#'     the initial and largest trust region radii, the acceptance threshold,
#'     the ratios that shrink and grow the region and the relative residual
#'     at which the conjugate gradient iteration stops.}
#'   \item{`Optim$Flat()`, `Optim$Sphere()`, `Optim$Stiefel(retraction = "SVD")`}{
#'     the manifolds; the retraction of Stiefel is `"SVD"` (the polar factor)
#'     or `"CholQR"` (the Cholesky orthogonalization).}
#'   \item{`Optim$PowerManifold(inner_manifold, inner_dims, outer_dims)`}{
#'     `prod(outer_dims)` copies of `inner_manifold`, each a block of
#'     `prod(inner_dims)` entries (of shape `inner_dims`) of the point.}
#'   \item{`Optim$ProductManifold(m1, m2, dims1, dims2)`}{the product of two
#'     manifolds: the first `prod(dims1)` entries lie on `m1`, the next
#'     `prod(dims2)` on `m2`.}
#' }
#' @return An optimization algorithm.
#' @examples
#' rosen <- function(x, p) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
#' rosg <- function(x, p) c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2), 200 * (x[2] - x[1]^2))
#' prob <- OptimizationProblem(rosen, c(-1.2, 1), grad = rosg)
#' solve(prob, ConjugateGradient())$u
#' solve(prob, Optim$BFGS())$u
#' solve(prob, Optim$KrylovTrustRegion(), maxiters = 40)$u
#' solve(prob, Optim$LBFGS(manifold = Optim$Sphere()))$u
#' bprob <- OptimizationProblem(rosen, c(0.5, 0.5), grad = rosg, lb = c(0, 0), ub = c(0.8, 0.8))
#' solve(bprob, Fminbox(), maxiters = 5)$u
#' @name optim_methods
NULL

#' @rdname optim_methods
#' @export
GradientDescent <- function(alphaguess = NULL, linesearch = NULL, P = NULL, precondprep = NULL, manifold = NULL) {
  .ojl_alg("GradientDescent", alphaguess = .ojl_alphaguess(alphaguess %||% .ojl_InitialPrevious()),
           linesearch = linesearch %||% .ojl_HagerZhang(), P = P, precondprep = precondprep,
           manifold = .ojl_as_manifold(manifold))
}

#' @rdname optim_methods
#' @export
ConjugateGradient <- function(alphaguess = NULL, linesearch = NULL, eta = 0.4, theta = 1, betamax = Inf, P = NULL,
                              precondprep = NULL, manifold = NULL) {
  .ojl_alg("ConjugateGradient", alphaguess = .ojl_alphaguess(alphaguess %||% .ojl_InitialHagerZhang()),
           linesearch = linesearch %||% .ojl_HagerZhang(), eta = eta, theta = theta, betamax = betamax, P = P,
           precondprep = precondprep, manifold = .ojl_as_manifold(manifold))
}

#' @rdname optim_methods
#' @export
AcceleratedGradientDescent <- function(alphaguess = NULL, linesearch = NULL, manifold = NULL) {
  .ojl_alg("AcceleratedGradientDescent", alphaguess = .ojl_alphaguess(alphaguess %||% .ojl_InitialPrevious()),
           linesearch = linesearch %||% .ojl_HagerZhang(), manifold = .ojl_as_manifold(manifold))
}

#' @rdname optim_methods
#' @export
MomentumGradientDescent <- function(mu = 0.01, alphaguess = NULL, linesearch = NULL, manifold = NULL) {
  .ojl_alg("MomentumGradientDescent", mu = mu, alphaguess = .ojl_alphaguess(alphaguess %||% .ojl_InitialPrevious()),
           linesearch = linesearch %||% .ojl_HagerZhang(), manifold = .ojl_as_manifold(manifold))
}

.ojl_BFGS <- function(alphaguess = NULL, linesearch = NULL, initial_invH = NULL, initial_stepnorm = NULL,
                      manifold = NULL) {
  .ojl_alg("BFGS", alphaguess = .ojl_alphaguess(alphaguess %||% .ojl_InitialStatic()),
           linesearch = linesearch %||% .ojl_HagerZhang(), initial_invH = initial_invH,
           initial_stepnorm = initial_stepnorm, manifold = .ojl_as_manifold(manifold))
}

.ojl_LBFGS <- function(m = 10L, alphaguess = NULL, linesearch = NULL, P = NULL, precondprep = NULL,
                       manifold = NULL, scaleinvH0 = is.null(P)) {
  .ojl_alg("LBFGS", m = as.integer(m), alphaguess = .ojl_alphaguess(alphaguess %||% .ojl_InitialStatic()),
           linesearch = linesearch %||% .ojl_HagerZhang(), P = P, precondprep = precondprep,
           scaleinvH0 = scaleinvH0, manifold = .ojl_as_manifold(manifold))
}

#' @rdname optim_methods
#' @export
Newton <- function(alphaguess = NULL, linesearch = NULL)
  .ojl_alg("Newton", alphaguess = .ojl_alphaguess(alphaguess %||% .ojl_InitialStatic()),
           linesearch = linesearch %||% .ojl_HagerZhang(), hess = TRUE)

#' @rdname optim_methods
#' @export
NewtonTrustRegion <- function(initial_delta = 1, delta_hat = 100, delta_min = 0, eta = 0.1, rho_lower = 0.25,
                              rho_upper = 0.75, use_fg = TRUE) {
  if (!(delta_hat > 0)) stop("maximum trust region radius must be positive", call. = FALSE)
  if (!(0 < initial_delta && initial_delta < delta_hat))
    stop("initial trust region radius must be positive and below the maximum trust region radius (", delta_hat,
         ")", call. = FALSE)
  if (!(delta_min >= 0)) stop("smallest allowable trust region radius must be non-negative", call. = FALSE)
  if (!(eta >= 0)) stop("the threshold eta must be non-negative", call. = FALSE)
  if (!(rho_lower > eta))
    stop("maximum threshold of actual and predicted reduction for shrinking the trust region must be greater ",
         "than the minimum threshold for accepting a step (", eta, ")", call. = FALSE)
  if (!(rho_upper > rho_lower))
    stop("minimum threshold of actual and predicted reduction for growing the trust region must be greater ",
         "than the minimum threshold for shrinking it (", rho_lower, ")", call. = FALSE)
  .ojl_alg("NewtonTrustRegion", initial_delta = initial_delta, delta_hat = delta_hat, delta_min = delta_min,
           eta = eta, rho_lower = rho_lower, rho_upper = rho_upper, use_fg = use_fg, hess = TRUE)
}

.ojl_ngmres_alg <- function(name, alphaguess, linesearch, nlprecon, nlpreconopts, epsilon0, wmax, manifold) {
  manifold <- .ojl_as_manifold(manifold)
  nlprecon <- nlprecon %||% GradientDescent(alphaguess = .ojl_InitialStatic(alpha = 1e-4, scaled = TRUE),
                                            linesearch = .ojl_Static(), manifold = manifold)
  if (!identical(manifold, nlprecon$manifold)) stop("AssertionError: manifold == nlprecon.manifold", call. = FALSE)
  if (!(nlprecon$name %in% c("GradientDescent", "LBFGS")))
    warning("Use caution. N-GMRES/O-ACCEL has only been tested with Gradient Descent and L-BFGS preconditioning.",
            call. = FALSE)
  nlpreconopts <- do.call(.ojl_options, utils::modifyList(list(iterations = 1L, allow_f_increases = TRUE),
                                                          as.list(nlpreconopts)))
  .ojl_alg(name, alphaguess = .ojl_alphaguess(alphaguess %||% .ojl_InitialStatic()),
           linesearch = linesearch %||% .ojl_HagerZhang(), nlprecon = nlprecon, nlpreconopts = nlpreconopts,
           epsilon0 = epsilon0, wmax = as.integer(wmax), manifold = manifold)
}

#' @rdname optim_methods
#' @export
NGMRES <- function(alphaguess = NULL, linesearch = NULL, nlprecon = NULL, nlpreconopts = list(), epsilon0 = 1e-12,
                   wmax = 10L, manifold = NULL)
  .ojl_ngmres_alg("NGMRES", alphaguess, linesearch, nlprecon, nlpreconopts, epsilon0, wmax, manifold)

#' @rdname optim_methods
#' @export
OACCEL <- function(alphaguess = NULL, linesearch = NULL, nlprecon = NULL, nlpreconopts = list(), epsilon0 = 1e-12,
                   wmax = 10L, manifold = NULL)
  .ojl_ngmres_alg("OACCEL", alphaguess, linesearch, nlprecon, nlpreconopts, epsilon0, wmax, manifold)

#' @rdname optim_methods
#' @export
Fminbox <- function(method = NULL, mu0 = NaN, mufactor = 0.001) {
  method <- method %||% .ojl_LBFGS()
  if (!inherits(method, "OptimAlgorithm"))
    stop("the inner method of Fminbox must be an Optim method", call. = FALSE)
  if (method$name %in% c("Newton", "NewtonTrustRegion"))
    stop("Newton is not supported as the Fminbox optimizer.", call. = FALSE)
  .ojl_alg("Fminbox", method = method, mu0 = mu0, mufactor = mufactor)
}

# ---------------------------------------------------------------------------
# Optim's Adam and AdaMax (first_order/adam.jl, adamax.jl), which
# OptimizationOptimJL does not export (Adam is the rule of Optimisers) but
# solves as Optim.Adam() and Optim.AdaMax(). alpha is a step or a function of
# the iteration (a schedule)

.ojl_adam_family <- function(name, alpha, beta_mean, beta_var, epsilon, manifold) {
  a <- .ojl_alg(name, alpha = alpha, beta_mean = beta_mean, beta_var = beta_var, epsilon = epsilon,
                manifold = .ojl_as_manifold(manifold))
  a$init <- function(m, o, d, x0) {
    x0 <- .ojl_retract(m$manifold, as.numeric(x0))
    st <- .ojl_first_order_state(d, x0, M = m$manifold)
    st$m <- st$g_x
    st$u <- rep(0, length(x0))
    st$alpha <- if (is.function(m$alpha)) m$alpha(1L) else m$alpha
    st$iter <- 0L
    st
  }
  a$update <- function(d, st, m) {
    st$iter <- st$iter + 1L
    if (is.function(m$alpha)) st$alpha <- m$alpha(st$iter)
    b1 <- m$beta_mean; b2 <- m$beta_var
    g <- st$g_x
    st$m <- b1 * st$m + (1 - b1) * g
    st$x_previous <- st$x
    st$f_x_previous <- st$f_x
    if (m$name == "Adam") {
      st$u <- b2 * st$u + (1 - b2) * (g * g)
      at <- st$alpha * sqrt(1 - .simple_powi(b2, st$iter)) / (1 - .simple_powi(b1, st$iter))
      st$x <- st$x - at * st$m / (sqrt(st$u) + m$epsilon)
    } else {
      st$u <- pmax(m$epsilon, pmax(b2 * st$u, abs(g)))
      st$x <- st$x - (st$alpha / (1 - .simple_powi(b1, st$iter))) * st$m / st$u
    }
    FALSE
  }
  # update_fgh! of the first order methods
  a$fgh <- function(d, st, m) {
    r <- d$value_gradient(st$x)
    st$g_x <- .ojl_project(m$manifold, r[[2L]], st$x)
    st$f_x <- r[[1L]]
  }
  a$reset <- function(m, st, d, x) {
    st$x <- .ojl_retract(m$manifold, x)
    r <- d$value_gradient(st$x)
    st$g_x <- .ojl_project(m$manifold, r[[2L]], st$x)
    st$f_x <- r[[1L]]
    na <- rep(NaN, length(x))
    st$x_previous <- na; st$f_x_previous <- NaN; st$s <- na
    st$m <- st$g_x
    st$u <- rep(0, length(x))
  }
  a
}

.ojl_Adam <- function(alpha = 1e-4, beta_mean = 0.9, beta_var = 0.999, epsilon = 1e-8, manifold = NULL)
  .ojl_adam_family("Adam", alpha, beta_mean, beta_var, epsilon, manifold)

.ojl_AdaMax <- function(alpha = 0.002, beta_mean = 0.9, beta_var = 0.999, epsilon = sqrt(.Machine$double.eps),
                        manifold = NULL)
  .ojl_adam_family("AdaMax", alpha, beta_mean, beta_var, epsilon, manifold)

# ---------------------------------------------------------------------------
# KrylovTrustRegion (second_order/krylov_trust_region.jl): a trust region
# method whose subproblem is solved by the Steihaug conjugate gradient
# iteration on Hessian vector products

.ojl_KrylovTrustRegion <- function(initial_radius = 1, max_radius = 100, eta = 0.1, rho_lower = 0.25,
                                   rho_upper = 0.75, cg_tol = 0.01) {
  a <- .ojl_alg("KrylovTrustRegion", initial_radius = initial_radius, max_radius = max_radius, eta = eta,
                rho_lower = rho_lower, rho_upper = rho_upper, cg_tol = cg_tol, hess = TRUE)
  a$init <- .ojl_ktr_init
  a$update <- .ojl_ktr_update
  a$assess <- .ojl_ktr_assess
  a
}

.ojl_ktr_init <- function(m, o, d, x0) {
  chk <- function(ok, what) if (!isTRUE(ok)) stop("AssertionError: ", what, call. = FALSE)
  chk(m$max_radius > 0, "method.max_radius > 0")
  chk(0 < m$initial_radius && m$initial_radius < m$max_radius, "0 < method.initial_radius < method.max_radius")
  chk(0 <= m$eta && m$eta < m$rho_lower, "0 <= method.eta < method.rho_lower")
  chk(m$rho_lower < m$rho_upper, "method.rho_lower < method.rho_upper")
  chk(m$rho_lower >= 0, "method.rho_lower >= 0")
  x0 <- as.numeric(x0)
  na <- rep(NaN, length(x0))
  r <- d$value_gradient(x0)
  .ojl_state(x = x0, f_x = r[[1L]], g_x = r[[2L]], x_previous = x0, f_x_previous = NaN, x_cache = na, s = na,
             interior = TRUE, accept_step = TRUE, radius = m$initial_radius, m_diff = 0, f_diff = 0, rho = 0,
             r = na, d = na, cg_iters = 0L)
}

# cg_steihaug!: the step in st$s and the change of the quadratic model
.ojl_cg_steihaug <- function(obj, st, m) {
  x <- st$x
  g <- st$g_x
  z <- rep(0, length(x))
  r <- g
  dd <- -r
  rho0 <- 0x1.249ad2594c37dp+332
  st$cg_iters <- 0L
  for (i in seq_along(x)) {
    st$cg_iters <- st$cg_iters + 1L
    Hd <- obj$hvp(x, dd)
    dHd <- .ojl_dot(dd, Hd)
    if (isTRUE(-1e-15 < dHd && dHd < 1e-15)) break
    alpha <- .ojl_dot(r, r) / dHd
    if (isTRUE(dHd < 0) || isTRUE(.ojl_norm2(z + alpha * dd) >= st$radius)) {
      a_ <- .ojl_dot(dd, dd)
      b_ <- 2 * .ojl_dot(z, dd)
      c_ <- .ojl_dot(z, z) - st$radius * st$radius
      tau <- (-b_ + sqrt(b_ * b_ - 4 * a_ * c_)) / (2 * a_)
      z <- z + tau * dd
      break
    }
    z <- z + alpha * dd
    rho_prev <- .ojl_dot(r, r)
    if (i == 1L) rho0 <- rho_prev
    r <- r + alpha * Hd
    rho_next <- .ojl_dot(r, r)
    dd <- -r + (rho_next / rho_prev) * dd
    if (isTRUE(rho_next / rho0 < m$cg_tol * m$cg_tol)) break
  }
  st$s <- z; st$r <- r; st$d <- dd
  Hd <- obj$hvp(x, z)
  .ojl_dot(g, z) + 0.5 * .ojl_dot(z, Hd)
}

.ojl_ktr_update <- function(d, st, m) {
  st$m_diff <- .ojl_cg_steihaug(d, st, m)
  if (!isTRUE(st$m_diff <= 0)) stop("AssertionError: state.m_diff <= 0", call. = FALSE)
  st$x_cache <- st$x + st$s
  f_x_cache <- d$value_bang(st$x_cache)
  st$f_diff <- f_x_cache - st$f_x
  st$rho <- st$f_diff / st$m_diff
  st$interior <- isTRUE(.ojl_norm2(st$s) < 0.9 * st$radius)
  if (isTRUE(st$rho < m$rho_lower)) {
    st$radius <- st$radius * 0.25
  } else if (isTRUE(st$rho > m$rho_upper) && !st$interior) {
    st$radius <- min(2 * st$radius, m$max_radius)
  }
  st$accept_step <- isTRUE(st$rho > m$eta)
  if (st$accept_step) {
    st$x_previous <- st$x
    st$f_x_previous <- st$f_x
    st$x <- st$x_cache
    st$f_x <- f_x_cache
    st$g_x <- d$gradient_bang(st$x)
  }
  FALSE
}

.ojl_ktr_assess <- function(st, d, o, m) {
  if (!st$accept_step) return(c(x = isTRUE(st$radius < o$x_abstol), f = FALSE, g = FALSE, inc = FALSE))
  ftol <- max(o$f_reltol * (abs(st$f_x) + o$f_reltol), .ojl_ulp(abs(st$f_x) + abs(st$f_x_previous)))
  c(x = isTRUE(max(abs(st$s)) < o$x_abstol), f = isTRUE(abs(st$f_diff) < ftol),
    g = isTRUE(max(abs(st$g_x)) < o$g_abstol), inc = FALSE)
}

#' @export
print.OptimAlgorithm <- function(x, ...) {
  if (x$name == "Fminbox") cat("Fminbox(", x$method$name, "())\n", sep = "") else cat(x$name, "()\n", sep = "")
  invisible(x)
}

# the options of a solve call (__map_optimizer_args of OptimizationOptimJL)
.ojl_solver_options <- function(alg, opts) {
  pass <- opts[intersect(names(opts), names(.ojl_default_options))]
  pass$callback <- NULL
  maxiters <- .opt_maxiters(opts$maxiters)
  local_maxiters <- .opt_maxiters(opts$local_maxiters)
  if (!is.null(opts$verbose)) pass$show_trace <- isTRUE(opts$verbose)
  if (!is.null(opts$abstol)) pass$f_abstol <- opts$abstol
  if (!is.null(maxiters)) {
    if (alg$name == "Fminbox") pass$outer_iterations <- maxiters else pass$iterations <- maxiters
  }
  if (!is.null(local_maxiters) && alg$name == "Fminbox") pass$iterations <- local_maxiters
  if (!is.null(opts$maxtime)) pass$time_limit <- opts$maxtime
  if (!is.null(opts$reltol)) pass$f_reltol <- opts$reltol
  do.call(.ojl_options, c(pass, list(method = if (alg$name == "Fminbox") alg$method else alg)))
}

.ojl_solve <- function(prob, alg, opts) {
  fs <- .opt_fns(prob)
  n <- length(prob$u0)
  bounded <- !is.null(prob$lb) || !is.null(prob$ub)
  lb <- if (is.null(prob$lb)) rep(-Inf, n) else as.numeric(prob$lb)
  ub <- if (is.null(prob$ub)) rep(Inf, n) else as.numeric(prob$ub)
  if (is.function(prob$f$cons) && alg$name != "IPNewton")
    stop(alg$name, "() does not handle constraints; use IPNewton()", call. = FALSE)
  if (bounded && alg$name == "SimulatedAnnealing")
    stop(IncompatibleOptimizerError("SimulatedAnnealing() does not handle bounds"))
  if (bounded && !(alg$name %in% c("Fminbox", "SAMIN", "IPNewton"))) {
    alg <- if (alg$name == "ParticleSwarm") ParticleSwarm(lower = lb, upper = ub, n_particles = alg$n_particles)
           else Fminbox(alg)
  }
  if (alg$name %in% c("Fminbox", "SAMIN") && !bounded)
    stop(alg$name, "() requires bounds: give `lb` and `ub` to the problem", call. = FALSE)
  o <- .ojl_solver_options(alg, opts)
  alg <- .ojl_with_dims(alg, dim(prob$u0) %||% c(n, 1L))
  hfun <- if (is.function(fs$fn$hess)) function(x) {
    H <- fs$fn$hess(x, fs$p)
    matrix(as.numeric(H), n, n)
  } else function(x) .ojl_opt_fd_hessian(fs$g, x)
  # the Hessian vector product of KrylovTrustRegion: the hv of the function
  # (x, v, p), else the Hessian times v
  hvfun <- if (is.function(fs$fn$hv)) {
    sg <- if (identical(unclass(prob$sense %||% MinSense), "MaxSense")) -1 else 1
    function(x, v) sg * as.numeric(fs$fn$hv(x, v, fs$p))
  }
  d <- .ojl_objective(fs$f, fs$g, hfun, n, hv = hvfun)
  t0 <- proc.time()[["elapsed"]]
  x0 <- as.numeric(prob$u0)
  res <- if (alg$name == "Fminbox") {
    .ojl_fminbox(d, lb, ub, x0, alg, o)
  } else if (alg$name == "SAMIN") {
    .ojl_samin(d, lb, ub, x0, .ojl_method_env(alg), o)
  } else if (alg$name == "IPNewton") {
    .ojl_ipnewton_solve(prob, fs, d, x0, lb, ub, .ojl_method_env(alg), o)
  } else {
    .ojl_optimize(d, x0, .ojl_method_env(alg), o)
  }
  sol <- .opt_solution(prob, alg, res$minimizer, res$minimum, .ojl_retcode(res), res)
  .opt_set_fields(sol, stats = OptimizationStats(iterations = res$iterations,
                                                time = proc.time()[["elapsed"]] - t0, fevals = res$f_calls,
                                                gevals = res$g_calls, hevals = res$h_calls))
}

# ---------------------------------------------------------------------------
# Zeroth order methods: SimulatedAnnealing, ParticleSwarm and Nelder Mead
# (zeroth_order/), whose iterations end with value!(d, state.x)

.ojl_zeroth_fgh <- function(d, st, m) st$f_x <- d$value_bang(st$x)

#' @rdname optim_methods
#' @export
SimulatedAnnealing <- function(neighbor = NULL, temperature = NULL, keep_best = TRUE) {
  a <- .ojl_alg("SimulatedAnnealing", neighbor = neighbor %||% function(x) x + stats::rnorm(length(x)),
                temperature = temperature %||% function(t) 1 / log(t), keep_best = keep_best, grad = FALSE)
  a$zeroth <- TRUE
  a$keep_x <- TRUE
  a$init <- function(m, o, d, x0) {
    f_x <- d$value_bang(x0)
    .ojl_state(x = x0, f_x = f_x, iteration = 1L, x_current = x0, f_x_current = f_x,
               x_proposal = rep(NaN, length(x0)), f_proposal = NaN)
  }
  a$update <- function(d, st, m) {
    t <- m$temperature(st$iteration)
    st$x_proposal <- as.numeric(m$neighbor(st$x_current))
    st$f_proposal <- d$value_bang(st$x_proposal)
    if (st$f_proposal <= st$f_x_current) {
      st$x_current <- st$x_proposal
      st$f_x_current <- st$f_proposal
      if (st$f_proposal < st$f_x) {
        st$x <- st$x_proposal
        st$f_x <- st$f_proposal
      }
    } else {
      p <- exp(-(st$f_proposal - st$f_x_current) / t)
      if (stats::runif(1L) <= p) {
        st$x_current <- st$x_proposal
        st$f_x_current <- st$f_proposal
      }
    }
    st$iteration <- st$iteration + 1L
    FALSE
  }
  a$fgh <- .ojl_zeroth_fgh
  a
}

#' @rdname optim_methods
#' @export
ParticleSwarm <- function(lower = NULL, upper = NULL, n_particles = 0L) {
  a <- .ojl_alg("ParticleSwarm", lower = lower, upper = upper, n_particles = as.integer(n_particles), grad = FALSE)
  a$zeroth <- TRUE
  a$init <- .ojl_ps_init
  a$update <- .ojl_ps_update
  a$fgh <- .ojl_zeroth_fgh
  a
}

.ojl_ps_init <- function(m, o, d, x0) {
  n <- length(x0)
  limit <- FALSE
  lower <- if (length(m$lower)) {
    limit <- TRUE
    as.numeric(m$lower)
  } else rep(-Inf, n)
  upper <- if (length(m$upper)) {
    limit <- TRUE
    as.numeric(m$upper)
  } else rep(Inf, n)
  if (length(lower) != n) stop("limits must be of same length as x_initial.", call. = FALSE)
  if (!all(upper > lower)) stop("upper must be greater than lower", call. = FALSE)
  np <- if (m$n_particles > 0L) {
    if (m$n_particles < 3L) {
      warning("Number of particles is set to 3 (minimum required)", call. = FALSE)
      3L
    } else m$n_particles
  } else max(3L, n)
  X <- V <- matrix(0, n, np)
  dx <- numeric(n)
  score <- numeric(np)
  f_x <- d$value_bang(x0)
  score[1L] <- f_x
  for (i in seq_len(np)) for (j in seq_len(n)) {
    if (limit) {
      ww <- upper[j] - lower[j]
      X[j, i] <- lower[j] + ww * stats::runif(1L)
      V[j, i] <- ww * (stats::runif(1L) * 2 - 1) / 10
    } else {
      # Optim reads x0[i] here (the particle index); the unknown j is meant
      if (i == 1L) dx[j] <- if (abs(x0[j]) > 0) abs(x0[j]) else 1
      X[j, i] <- x0[j] + dx[j] * stats::runif(1L)
      V[j, i] <- abs(X[j, i]) * (stats::runif(1L) * 2 - 1)
    }
  }
  X[, 1L] <- x0
  X_best <- X
  if (np > 1L) for (i in 2:np) score[i] <- d$value(X[, i])
  .ojl_state(x = x0, f_x = f_x, iteration = 0L, lower = lower, upper = upper, c1 = 2, c2 = 2, w = 1,
             limit = limit, n_particles = np, X = X, V = V, X_best = X_best, score = score,
             best_score = numeric(np), x_learn = x0, current_state = 0L, iterations = o$iterations)
}

.ojl_ps_mu <- function(f) {
  mu1 <- if (0 <= f && f <= 4 / 10) 0 else if (4 / 10 < f && f <= 6 / 10) 5 * f - 2 else
    if (6 / 10 < f && f <= 7 / 10) 1 else if (7 / 10 < f && f <= 8 / 10) -10 * f + 8 else 0
  mu2 <- if (0 <= f && f <= 2 / 10) 0 else if (2 / 10 < f && f <= 3 / 10) 10 * f - 2 else
    if (3 / 10 < f && f <= 4 / 10) 1 else if (4 / 10 < f && f <= 6 / 10) -5 * f + 3 else 0
  mu3 <- if (0 <= f && f <= 1 / 10) 1 else if (1 / 10 < f && f <= 3 / 10) -5 * f + 3 / 2 else 0
  mu4 <- if (0 <= f && f <= 7 / 10) 0 else if (7 / 10 < f && f <= 9 / 10) 5 * f - 7 / 2 else 1
  c(mu1, mu2, mu3, mu4)
}

.ojl_ps_next_state <- function(prev, mu) {
  if (prev == 0L) return(which.max(mu))
  if (prev == 1L) return(if (mu[1L] > 0) 1L else if (mu[2L] > 0) 2L else if (mu[4L] > 0) 4L else 3L)
  if (prev == 2L) return(if (mu[2L] > 0) 2L else if (mu[3L] > 0) 3L else if (mu[1L] > 0) 1L else 4L)
  if (prev == 3L) return(if (mu[3L] > 0) 3L else if (mu[4L] > 0) 4L else if (mu[2L] > 0) 2L else 1L)
  if (mu[4L] > 0) 4L else if (mu[1L] > 0) 1L else if (mu[2L] > 0) 2L else 3L
}

.ojl_ps_update <- function(d, st, m) {
  n <- length(st$x)
  np <- st$n_particles
  if (st$iteration == 0L) st$best_score <- st$score
  # housekeeping!
  Fb <- st$f_x
  for (i in seq_len(np)) {
    if (st$score[i] <= st$best_score[i]) {
      st$best_score[i] <- st$score[i]
      st$X_best[, i] <- st$X[, i]
      if (st$score[i] <= Fb) {
        st$x <- st$X[, i]
        Fb <- st$score[i]
      }
    }
  }
  st$f_x <- Fb
  i_worst <- which.max(st$score)
  st$x_learn <- st$x
  ri <- sample.int(n, 1L)
  stats::rnorm(1L)
  sigma_learn <- 1 - (1 - 0.1) * st$iteration / st$iterations
  r3 <- stats::rnorm(1L) * sigma_learn
  if (st$limit) {
    st$x_learn[ri] <- st$x_learn[ri] + (st$upper[ri] - st$lower[ri]) / 3 * r3
    st$x_learn[ri] <- min(max(st$x_learn[ri], st$lower[ri]), st$upper[ri])
  } else {
    st$x_learn[ri] <- st$x_learn[ri] + st$x_learn[ri] * r3
  }
  score_learn <- d$value_bang(st$x_learn)
  if (score_learn < st$f_x) {
    st$x <- st$x_learn
    st$X_best[, i_worst] <- st$x_learn
    st$X[, i_worst] <- st$x_learn
    st$f_x <- score_learn
    st$score[i_worst] <- score_learn
    st$best_score[i_worst] <- score_learn
  }
  # get_swarm_state
  X <- st$X
  i_best <- which.min(st$score)
  XtX <- crossprod(X)
  trc <- sum(diag(XtX))
  dd <- colSums(XtX)
  dd <- sqrt(pmax(np * diag(XtX) + trc - 2 * dd, 0))
  f <- (dd[i_best] - min(dd)) / max(max(dd) - min(dd), sqrt(.Machine$double.eps))
  cs <- .ojl_ps_next_state(st$current_state, .ojl_ps_mu(f))
  st$current_state <- cs
  # update_swarm_params!
  dc1 <- 5 / 100 + stats::runif(1L) / 20
  dc2 <- 5 / 100 + stats::runif(1L) / 20
  c1 <- st$c1; c2 <- st$c2
  if (cs == 1L) {
    c1 <- c1 + dc1; c2 <- c2 - dc2
  } else if (cs == 2L) {
    c1 <- c1 + dc1 / 2; c2 <- c2 - dc2 / 2
  } else if (cs == 3L) {
    c1 <- c1 + dc1 / 2; c2 <- c2 + dc2 / 2
  } else if (cs == 4L) {
    c1 <- c1 - dc1; c2 <- c2 - dc2
  }
  c1 <- min(max(c1, 3 / 2), 5 / 2)
  # Optim sets c2 to 5/2 when it falls below 3/2; the lower bound 3/2 is meant
  c2 <- min(max(c2, 3 / 2), 5 / 2)
  if (c1 + c2 > 4) {
    ct <- c1 + c2
    c1 <- c1 / ct * 4; c2 <- c2 / ct * 4
  }
  st$w <- 1 / (1 + 3 / 2 * exp(-26 / 10 * f))
  st$c1 <- c1; st$c2 <- c2
  # update_swarm!
  for (i in seq_len(np)) for (j in seq_len(n)) {
    r1 <- stats::runif(1L); r2 <- stats::runif(1L)
    vx <- st$X_best[j, i] - st$X[j, i]
    vg <- st$x[j] - st$X[j, i]
    st$V[j, i] <- st$V[j, i] * st$w + c1 * r1 * vx + c2 * r2 * vg
    st$X[j, i] <- st$X[j, i] + st$V[j, i]
  }
  if (st$limit) st$X <- pmin(pmax(st$X, st$lower), st$upper)
  for (i in seq_len(np)) st$score[i] <- d$value(st$X[, i])
  st$iteration <- st$iteration + 1L
  FALSE
}

.ojl_nm_objective <- function(y, m, n) {
  mu <- sum(y) / length(y)
  sqrt(sum((y - mu)^2) / (length(y) - 1L) * (m / n))
}

.ojl_nm_centroid <- function(simplex, h) {
  n <- length(simplex[[1L]])
  cc <- numeric(n)
  for (i in seq_along(simplex)) if (i != h) cc <- cc + simplex[[i]]
  cc * (1 / n)
}

.ojl_NelderMead <- function(parameters = NULL, initial_simplex = NULL) {
  a <- .ojl_alg("NelderMead", parameters = parameters, initial_simplex = initial_simplex, grad = FALSE)
  a$zeroth <- TRUE
  a$keep_x <- TRUE
  a$init <- function(m, o, d, x0) {
    n <- length(x0); mm <- n + 1L
    simplex <- .ojl_nm_simplex(m, x0)
    f_x <- d$value_bang(simplex[[1L]])
    fs <- numeric(mm)
    fs[1L] <- f_x
    for (i in 2:mm) fs[i] <- d$value(simplex[[i]])
    io <- order(fs)
    pr <- m$parameters
    par <- if (is.null(pr)) c(1, 1 + 2 / n, 0.75 - 1 / (2 * n), 1 - 1 / n)
           else c(pr$alpha %||% 1, pr$beta %||% 2, pr$gamma %||% 0.5, pr$delta %||% 0.5)
    .ojl_state(x = x0, f_x = f_x, m = mm, simplex = simplex, x_centroid = .ojl_nm_centroid(simplex, io[mm]),
               f_simplex = fs, nm_x = .ojl_nm_objective(fs, n, mm), f_lowest = fs[io[1L]], i_order = io,
               a = par[1L], b = par[2L], g = par[3L], dl = par[4L], step_type = "initial")
  }
  a$initial_convergence <- function(st, o) {
    nmo <- .ojl_nm_objective(st$f_simplex, st$m, length(st$x))
    list(isTRUE(nmo <= o$g_abstol), !is.finite(st$f_x) || !is.finite(nmo))
  }
  a$assess <- function(st, d, o, m) c(x = FALSE, f = FALSE, g = isTRUE(st$nm_x <= o$g_abstol), inc = FALSE)
  a$g_residual <- function(st) st$nm_x
  a$update <- .ojl_nm_update
  a$fgh <- .ojl_zeroth_fgh
  a$reset <- function(m, st, d, x) {
    st$simplex <- .ojl_nm_simplex(m, x)
    st$f_x <- d$value_bang(st$simplex[[1L]])
    st$f_simplex[1L] <- st$f_x
    for (i in 2:length(st$simplex)) st$f_simplex[i] <- d$value(st$simplex[[i]])
    st$i_order <- order(st$f_simplex)
  }
  a$after_while <- function(d, st, m, o) {
    st$i_order <- order(st$f_simplex)
    xc <- .ojl_nm_centroid(st$simplex, st$i_order[st$m])
    fc <- d$value(xc)
    i <- which.min(st$f_simplex)
    x_min <- st$simplex[[i]]; f_min <- st$f_simplex[i]
    if (isTRUE(fc < f_min)) {
      x_min <- xc; f_min <- fc
    }
    st$x <- x_min
    st$f_x <- f_min
  }
  a
}

.ojl_nm_simplex <- function(m, x0) {
  if (is.function(m$initial_simplex)) return(lapply(m$initial_simplex(x0), as.numeric))
  a <- 0.025; b <- 0.5
  if (is.list(m$initial_simplex)) {
    a <- m$initial_simplex$a %||% a; b <- m$initial_simplex$b %||% b
  }
  s <- rep(list(x0), length(x0) + 1L)
  for (j in seq_along(x0)) s[[j + 1L]][j] <- (1 + b) * s[[j + 1L]][j] + a
  s
}

.ojl_nm_update <- function(d, st, m) {
  shrink <- FALSE
  n <- length(st$x); mm <- st$m
  io <- st$i_order
  st$x_centroid <- .ojl_nm_centroid(st$simplex, io[mm])
  x_lowest <- st$simplex[[io[1L]]]
  x_highest <- st$simplex[[io[mm]]]
  st$f_lowest <- st$f_simplex[io[1L]]
  f_second <- st$f_simplex[io[n]]
  f_highest <- st$f_simplex[io[mm]]
  xc <- st$x_centroid
  x_reflect <- xc + st$a * (xc - x_highest)
  f_reflect <- d$value(x_reflect)
  if (f_reflect < st$f_lowest) {
    x_cache <- xc + st$b * (x_reflect - xc)
    f_expand <- d$value(x_cache)
    if (f_expand < f_reflect) {
      st$simplex[[io[mm]]] <- x_cache
      st$f_simplex[io[mm]] <- f_expand
      st$step_type <- "expansion"
    } else {
      st$simplex[[io[mm]]] <- x_reflect
      st$f_simplex[io[mm]] <- f_reflect
      st$step_type <- "reflection"
    }
    st$i_order <- c(io[mm], io[-mm])
  } else if (f_reflect < f_second) {
    st$simplex[[io[mm]]] <- x_reflect
    st$f_simplex[io[mm]] <- f_reflect
    st$step_type <- "reflection"
    st$i_order <- order(st$f_simplex)
  } else {
    if (f_reflect < f_highest) {
      x_cache <- xc + st$g * (x_reflect - xc)
      f_oc <- d$value(x_cache)
      if (f_oc < f_reflect) {
        st$simplex[[io[mm]]] <- x_cache
        st$f_simplex[io[mm]] <- f_oc
        st$step_type <- "outside contraction"
        st$i_order <- order(st$f_simplex)
      } else shrink <- TRUE
    } else {
      x_cache <- xc - st$g * (x_reflect - xc)
      f_ic <- d$value(x_cache)
      if (f_ic < f_highest) {
        st$simplex[[io[mm]]] <- x_cache
        st$f_simplex[io[mm]] <- f_ic
        st$step_type <- "inside contraction"
        st$i_order <- order(st$f_simplex)
      } else shrink <- TRUE
    }
  }
  if (shrink) {
    io <- st$i_order
    for (i in 2:mm) {
      ord <- io[i]
      st$simplex[[ord]] <- x_lowest + st$dl * (st$simplex[[ord]] - x_lowest)
      st$f_simplex[ord] <- d$value(st$simplex[[ord]])
    }
    st$step_type <- "shrink"
    st$i_order <- order(st$f_simplex)
  }
  st$nm_x <- .ojl_nm_objective(st$f_simplex, n, mm)
  FALSE
}

#' @rdname optim_methods
#' @export
SAMIN <- function(nt = 5L, ns = 5L, t0 = 2, rt = 0.9, r_expand = 10, bounds_ratio = 0.6, neps = 5L,
                  coverage_ok = FALSE, verbosity = 0L)
  .ojl_alg("SAMIN", nt = as.integer(nt), ns = as.integer(ns), t0 = t0, rt = rt, r_expand = r_expand,
           bounds_ratio = bounds_ratio, neps = as.integer(neps), coverage_ok = coverage_ok,
           verbosity = as.integer(verbosity), grad = FALSE)

# optimize(obj_fn, lb, ub, x0, ::SAMIN, options) of constrained/samin.jl.
# Optim falls off the end of the function (it returns nothing) when the
# annealing converges or stops on wide bounds; here the best point is
# returned with the convergence flags.
.ojl_samin <- function(d, lb, ub, x0, m, o) {
  t0 <- proc.time()[["elapsed"]]
  x <- as.numeric(x0)
  n <- length(x)
  f_x <- d$value_bang(x)
  x_tol <- o$x_abstol; f_tol <- o$f_abstol
  t <- m$t0
  converge <- 0L
  x_converged <- FALSE; f_converged <- FALSE
  x_abs <- Inf; f_abs <- Inf
  fstar <- rep(Inf, m$neps)
  f_opt <- f_x
  x_opt <- x
  bounds <- ub - lb
  coverage_ok <- m$coverage_ok
  for (i in seq_len(n)) if (x[i] > ub[i] || x[i] < lb[i])
    stop("samin: initial parameter ", i, " out of bounds", call. = FALSE)
  result <- function(stopped) {
    list(minimizer = x_opt, minimum = f_opt, iterations = d$f_calls, g_residual = NaN, f_calls = d$f_calls,
         g_calls = d$df_calls, h_calls = d$h_calls, x_abschange = x_abs, x_relchange = NaN, f_abschange = f_abs,
         f_relchange = NaN, samin_converge = converge,
         stopped_by = list(x_converged = x_converged, f_converged = f_converged, g_converged = FALSE,
                           iterations = stopped && d$f_calls >= o$iterations, callback = FALSE,
                           f_increased = FALSE, ls_failed = FALSE))
  }
  while (converge == 0L) {
    nup <- 0L; nrej <- 0L; ndown <- 0L; lnobds <- 0L
    for (mi in seq_len(m$nt)) {
      nacp <- numeric(n)
      for (j in seq_len(m$ns)) {
        for (h in seq_len(n)) {
          if (lb[h] != ub[h]) {
            xp <- x
            xp[h] <- xp[h] + (2 * stats::runif(1L) - 1) * bounds[h]
            if (xp[h] < lb[h] || xp[h] > ub[h]) {
              xp[h] <- lb[h] + (ub[h] - lb[h]) * stats::runif(1L)
              lnobds <- lnobds + 1L
            }
            f_prop <- d$value_bang(xp)
            if (f_prop <= f_x) {
              x <- xp; f_x <- f_prop
              nacp[h] <- nacp[h] + 1; nup <- nup + 1L
              if (f_prop < f_opt) {
                x_opt <- xp; f_opt <- f_prop
              }
            } else {
              p <- exp(-(f_prop - f_x) / t)
              if (stats::runif(1L) < p) {
                x <- xp; f_x <- f_prop
                nacp[h] <- nacp[h] + 1; ndown <- ndown + 1L
              } else nrej <- nrej + 1L
            }
          }
          if (d$f_calls >= o$iterations || isTRUE(proc.time()[["elapsed"]] - t0 > o$time_limit)) {
            x_converged <- x_abs < 0
            f_converged <- f_abs <= 0
            return(result(TRUE))
          }
        }
      }
      test <- 0L
      for (i in seq_len(n)) {
        if (lb[i] != ub[i]) {
          ratio <- nacp[i] / m$ns
          if (ratio > m$bounds_ratio) bounds[i] <- bounds[i] * (1 + 2 * (ratio - 0.6) / 0.4)
          if (ratio < 1 - m$bounds_ratio) bounds[i] <- bounds[i] / (1 + 2 * ((0.4 - ratio) / 0.4))
          if (bounds[i] > ub[i] - lb[i]) {
            bounds[i] <- ub[i] - lb[i]
            test <- test + 1L
          }
        } else test <- test + 1L
      }
      if (!coverage_ok) coverage_ok <- test == n
    }
    if (coverage_ok) {
      fstar[1L] <- f_x
      f_abs <- abs(f_opt - f_x)
      if (all(abs(f_opt - fstar) < f_tol)) {
        f_converged <- TRUE
        if (any(bounds > x_tol)) {
          x_converged <- FALSE
          converge <- 0L
          break
        } else {
          converge <- 1L
          x_converged <- TRUE
          x_abs <- max(bounds)
        }
      } else f_converged <- FALSE
      if (converge == 1L && lnobds > 0L) converge <- 2L
      t <- t * m$rt
      fstar <- c(f_x, fstar)[seq_len(m$neps)]
    } else {
      t <- t * m$r_expand
      if (m$neps >= 2L) for (i in m$neps:2) fstar[i] <- fstar[i - 1L]
    }
    f_x <- f_opt
    x <- x_opt
  }
  result(FALSE)
}


# ---------------------------------------------------------------------------
# IPNewton (constrained/ipnewton): the interior point primal dual Newton
# method of Optim, with the constraint bounds of NLSolversBase

# parse_constraints: l <= v <= u as equalities and one sided inequalities
# sigma (v - b) >= 0
.ojl_parse_constraints <- function(l, u) {
  eq <- integer(); val <- numeric(); ineq <- integer(); sg <- numeric(); b <- numeric()
  if (length(l) != length(u))
    stop("Size of lower bounds (", length(l), ") and number of upper bounds (", length(u), ") must be equal.",
         call. = FALSE)
  for (i in seq_along(l)) {
    li <- l[i]; ui <- u[i]
    if (li == ui) {
      eq <- c(eq, i); val <- c(val, ui)
    } else if (li < ui) {
      if (is.finite(li)) {
        ineq <- c(ineq, i); sg <- c(sg, 1); b <- c(b, li)
      }
      if (is.finite(ui)) {
        ineq <- c(ineq, i); sg <- c(sg, -1); b <- c(b, ui)
      }
    } else stop("Lower bound (", li, ") must not be greater than upper bound (", ui, ").", call. = FALSE)
  }
  list(eq = eq, val = val, ineq = ineq, sigma = sg, b = b)
}

.ojl_symmetrize <- function(l, u) {
  if (!length(l) && length(u)) l <- rep(-Inf, length(u))
  if (length(l) && !length(u)) u <- rep(Inf, length(l))
  list(as.numeric(l), as.numeric(u))
}

# ConstraintBounds(lx, ux, lc, uc)
.ojl_constraint_bounds <- function(lx, ux, lc, uc) {
  x <- do.call(.ojl_parse_constraints, .ojl_symmetrize(lx, ux))
  cs <- .ojl_symmetrize(lc, uc)
  cc <- do.call(.ojl_parse_constraints, cs)
  list(nc = length(cs[[1L]]), eqx = x$eq, valx = x$val, ineqx = x$ineq, sx = x$sigma, bx = x$b,
       eqc = cc$eq, valc = cc$val, ineqc = cc$ineq, sc = cc$sigma, bc = cc$b)
}

# the pivoted cholesky(Positive, A, Val{true}) of PositiveFactorizations.jl:
# the lower factor L and the permutation p with A[p, p] = L L' after the
# negative and tiny pivots are replaced
.ojl_positive_cholesky_pivoted <- function(A) {
  KA <- nrow(A)
  if (!KA) return(list(L = A, piv = integer()))
  tol <- 10 * KA * .Machine$double.eps * max(abs(A))
  d <- numeric(KA)
  piv <- seq_len(KA)
  Ad <- diag(A)
  j0 <- 1L
  while (j0 <= KA) {
    jend <- min(KA, j0 + 31L)
    for (j in j0:jend) {
      Amax <- 0
      jmax <- j - 1L
      for (jj in j:KA) {
        tmp <- abs(Ad[jj])
        if (isTRUE(tmp > Amax)) {
          Amax <- tmp
          jmax <- jj
        }
      }
      if (jmax > j) {
        A <- .ojl_pivot(A, j, jmax)
        Ad[c(j, jmax)] <- Ad[c(jmax, j)]
        piv[c(j, jmax)] <- piv[c(jmax, j)]
      }
      Ajj <- A[j, j]
      ks <- if (j > j0) j0:(j - 1L) else integer()
      for (k in ks) {
        tmp <- A[j, k]
        Ajj <- Ajj - d[k] * tmp * tmp
      }
      if (isTRUE(abs(Ajj) > tol)) {
        d[j] <- sign(Ajj)
        s <- sqrt(abs(Ajj))
        A[j, j] <- s
        f <- d[j] / s
      } else {
        d[j] <- 0
        A[j, j] <- 1
        f <- 0
      }
      if (j < KA) {
        ii <- (j + 1L):KA
        for (k in ks) A[ii, j] <- A[ii, j] - (d[k] * A[j, k]) * A[ii, k]
        tmp <- A[ii, j] * f
        A[ii, j] <- tmp
        Ad[ii] <- Ad[ii] - d[j] * tmp * tmp
      }
    }
    if (jend < KA) {
      # the block update B22 -= B21 diag(d) B21' (a BLAS syrk in Julia)
      ii <- (jend + 1L):KA
      C <- A[ii, j0:jend, drop = FALSE]
      U <- C %*% (d[j0:jend] * t(C))
      B22 <- A[ii, ii, drop = FALSE]
      low <- lower.tri(B22, diag = TRUE)
      B22[low] <- B22[low] - U[low]
      A[ii, ii] <- B22
    }
    j0 <- jend + 1L
  }
  A[upper.tri(A)] <- 0
  list(L = A, piv = piv)
}

# pivot! of PositiveFactorizations.jl: swap rows and columns i and j of the
# lower triangle
.ojl_pivot <- function(A, i, j) {
  i0 <- min(i, j); j <- max(i, j); i <- i0
  K <- nrow(A)
  for (k in seq_len(i - 1L)) A[c(i, j), k] <- A[c(j, i), k]
  t1 <- A[i, i]; A[i, i] <- A[j, j]; A[j, j] <- t1
  if (j - i > 1L) for (k in (i + 1L):(j - 1L)) {
    t1 <- A[k, i]; A[k, i] <- A[j, k]; A[j, k] <- t1
  }
  if (j < K) for (k in (j + 1L):K) A[k, c(i, j)] <- A[k, c(j, i)]
  A
}

# F \ b for a factor list(L, piv) (piv NULL: no pivoting)
.ojl_chol_solve <- function(Fc, b) {
  vec <- is.null(dim(b))
  B <- if (vec) matrix(b, ncol = 1L) else b
  if (!nrow(Fc$L) || !ncol(B)) return(if (vec) as.numeric(B) else B)
  P <- Fc$piv
  if (length(P)) B <- B[P, , drop = FALSE]
  X <- backsolve(Fc$L, forwardsolve(Fc$L, B), upper.tri = FALSE, transpose = TRUE)
  if (length(P)) X[P, ] <- X
  if (vec) as.numeric(X) else X
}

# Matrix(F) of the pivoted factor
.ojl_chol_matrix <- function(Fc) {
  M <- tcrossprod(Fc$L)
  if (!length(Fc$piv)) return(M)
  ip <- order(Fc$piv)
  M[ip, ip, drop = FALSE]
}

.ojl_ip_jx <- function(idx, n) {
  M <- matrix(0, length(idx), n)
  if (length(idx)) M[cbind(seq_along(idx), idx)] <- 1
  M
}

.ojl_ip_JE <- function(B, st)
  rbind(.ojl_ip_jx(B$eqx, length(st$x)), st$constr_J[B$eqc, , drop = FALSE])

.ojl_ip_JI <- function(B, st)
  rbind(.ojl_ip_jx(B$ineqx, length(st$x)), st$constr_J[B$ineqc, , drop = FALSE])

.ojl_ip_userlambda <- function(lcI, B) {
  lam <- numeric(B$nc)
  for (i in seq_along(B$ineqc)) lam[B$ineqc[i]] <- lam[B$ineqc[i]] - lcI[i] * B$sc[i]
  lam
}

.ojl_ip_bv <- function(v) {
  if (!length(v)) return(0)
  lv <- rep(-Inf, length(v))
  pos <- !is.na(v) & v > 0
  lv[pos] <- log(v[pos])
  -sum(lv)
}

.ojl_ip_ev <- function(B, x, cc, bs) {
  (sum(bs$lx * (bs$slack_x - B$sx * (x[B$ineqx] - B$bx))) +
     sum(bs$lc * (bs$slack_c - B$sc * (cc[B$ineqc] - B$bc)))) +
    sum(bs$lxE * (B$valx - x[B$eqx])) + sum(bs$lcE * (B$valc - cc[B$eqc]))
}

# lagrangian_fg!: the value of the Lagrangian with the barrier and the
# gradients in x and in the slack and multiplier variables
.ojl_ip_lagrangian_fg <- function(B, x, f_x, g_x, cc, J, bs, mu) {
  ev <- .ojl_ip_ev(B, x, cc, bs)
  L <- f_x + mu * (.ojl_ip_bv(bs$slack_x) + .ojl_ip_bv(bs$slack_c)) + ev
  gx <- g_x
  bg <- lapply(bs, function(v) numeric(length(v)))
  bg$slack_x <- bg$slack_x - mu / bs$slack_x
  bg$slack_c <- bg$slack_c - mu / bs$slack_c
  for (i in seq_along(B$ineqx)) {
    li <- bs$lx[i]
    bg$slack_x[i] <- bg$slack_x[i] + li
    j <- B$ineqx[i]
    gx[j] <- gx[j] - li * B$sx[i]
  }
  bg$slack_c <- bg$slack_c + bs$lc
  if (length(B$ineqc)) gx <- gx - as.numeric(crossprod(J[B$ineqc, , drop = FALSE], bs$lc * B$sc))
  gx[B$eqx] <- gx[B$eqx] - bs$lxE
  if (length(B$eqc)) gx <- gx - as.numeric(crossprod(J[B$eqc, , drop = FALSE], bs$lcE))
  bg$lx <- bg$lx + (bs$slack_x - B$sx * (x[B$ineqx] - B$bx))
  bg$lc <- bg$lc + (bs$slack_c - B$sc * (cc[B$ineqc] - B$bc))
  bg$lxE <- bg$lxE + (B$valx - x[B$eqx])
  bg$lcE <- bg$lcE + (B$valc - cc[B$eqc])
  list(L = L, ev = ev, g = gx, bgrad = bg)
}

.ojl_ip_slopealpha <- function(sx, gx, bstep, bgrad) {
  r <- sum(sx * gx)
  for (nm in c("slack_x", "slack_c", "lx", "lc", "lxE", "lcE")) r <- r + sum(bstep[[nm]] * bgrad[[nm]])
  r
}

# update_fgh!: the objective, the Lagrangian, the modified gradient and the
# factored Hessian of the Lagrangian at state.x
.ojl_ip_fgh <- function(d, cons, B, st) {
  r <- d$value_gradient_hessian(st$x)
  st$f_x <- r[[1L]]; st$g_x <- r[[2L]]; st$H_x <- r[[3L]]
  lf <- .ojl_ip_lagrangian_fg(B, st$x, st$f_x, st$g_x, st$constr_c, st$constr_J, st$bstate, st$mu)
  st$L_x <- lf$L; st$ev <- lf$ev; st$g_L_x <- lf$g; st$bgrad <- lf$bgrad
  bs <- st$bstate; bg <- st$bgrad
  gt <- st$g_L_x
  if (length(B$ineqc)) {
    Hssc <- bs$lc / bs$slack_c
    gt <- gt + as.numeric(crossprod(st$constr_J[B$ineqc, , drop = FALSE], B$sc * (bs$lc - Hssc * bg$lc)))
  }
  for (i in seq_along(B$ineqx)) {
    j <- B$ineqx[i]
    gt[j] <- gt[j] + B$sx[i] * (bs$lx[i] - bg$lx[i] * bs$lx[i] / bs$slack_x[i])
  }
  st$gtilde <- gt
  lam <- .ojl_ip_userlambda(bs$lc, B)
  lam[B$eqc] <- -bs$lcE
  Hxx <- cons$h(st$H_x, st$x, lam)
  JIc <- st$constr_J[B$ineqc, , drop = FALSE]
  Hxx <- Hxx + crossprod(JIc * (bs$lc / bs$slack_c), JIc)
  for (i in seq_along(B$ineqx)) {
    j <- B$ineqx[i]
    Hxx[j, j] <- Hxx[j, j] + bs$lx[i] / bs$slack_x[i]
  }
  st$H_L_x <- Hxx
  st$Htilde <- .ojl_positive_cholesky_pivoted(Hxx)
  invisible(NULL)
}

# initialize_mu_lambda!
.ojl_ip_init_mu_lambda <- function(B, st, HcI, mu0, beta = 1 / 100) {
  bs <- st$bstate
  if (B$nc == 0L && !length(unique(c(B$ineqx, B$eqx)))) {
    st$mu <- 0
    st$bstate <- lapply(bs, function(v) numeric(length(v)))
    return(invisible(NULL))
  }
  gf <- st$g_x
  st$g_L_x <- gf
  JE <- .ojl_ip_JE(B, st)
  Cc <- list(L = .ojl_positive_cholesky(tcrossprod(JE)), piv = NULL)
  perp <- function(v) v - as.numeric(crossprod(JE, .ojl_chol_solve(Cc, as.numeric(JE %*% v))))
  Pg <- perp(gf)
  JI <- .ojl_ip_JI(B, st)
  s <- c(bs$slack_x, bs$slack_c)
  if (is.null(mu0) || identical(mu0, "auto")) {
    Dperp <- sum(Pg * Pg)
    PDg <- perp(as.numeric(crossprod(JI, c(B$sx, B$sc) / s)))
    DI <- sum(PDg * PDg)
    y <- as.numeric(JI %*% Pg) / s
    kperp <- sum(as.numeric(st$H_x %*% Pg) * Pg)
    kI <- sum(as.numeric(HcI %*% Pg) * Pg) + sum(y * y)
    mu <- beta * (if (isTRUE(kperp == 0)) sqrt(Dperp / DI) else min(sqrt(Dperp / DI), abs(kperp / kI)))
    if (!is.finite(mu)) {
      PDgt <- perp(as.numeric(crossprod(JI, 1 / s)))
      mu <- beta * sqrt(Dperp / sum(PDgt * PDgt))
    }
    if (!is.finite(mu) || mu == 0) mu <- 1
  } else mu <- as.numeric(mu0)
  st$mu <- mu
  bs$lx <- mu / bs$slack_x
  bs$lc <- mu / bs$slack_c
  dbI <- gf - as.numeric(crossprod(JI, c(bs$lx, bs$lc)))
  lE <- .ojl_chol_solve(Cc, as.numeric(JE %*% dbI)) +
    (c(B$valx, B$valc) - c(st$x[B$eqx], st$constr_c[B$eqc])) / mu
  bs$lxE <- lE[seq_along(B$eqx)]
  bs$lcE <- lE[length(B$eqx) + seq_along(B$eqc)]
  st$bstate <- bs
  invisible(NULL)
}

.ojl_ip_maxstep <- function(amax, x, s) {
  for (i in seq_along(s)) if (isTRUE(s[i] < 0)) amax <- min(amax, -x[i] / s[i])
  amax
}

.ojl_ip_small <- function(ref, step)
  isTRUE(all(step == 0 | abs(step) < vapply(ref, .ojl_ulp, numeric(1))))

.ojl_ip_solve_slack <- function(bstep, s, B, bs, bg, J, mu) {
  for (i in seq_along(B$ineqx)) {
    bstep$slack_x[i] <- -bg$lx[i] + B$sx[i] * s[B$ineqx[i]]
    bstep$lx[i] <- -(-mu / bs$slack_x[i] + bs$lx[i]) - bs$lx[i] * bstep$slack_x[i] / bs$slack_x[i]
  }
  SJ <- B$sc * as.numeric(J[B$ineqc, , drop = FALSE] %*% s)
  for (i in seq_along(bstep$lc)) {
    bstep$slack_c[i] <- -bg$lc[i] + SJ[i]
    bstep$lc[i] <- -(-mu / bs$slack_c[i] + bs$lc[i]) - bs$lc[i] * bstep$slack_c[i] / bs$slack_c[i]
  }
  bstep
}

# solve_step!: the Newton step of the primal dual system; NULL when the
# step is not a descent step (the step is then zero)
.ojl_ip_solve_step <- function(st, B) {
  mu <- st$mu
  bs <- st$bstate; bg <- st$bgrad
  J <- st$constr_J; Ht <- st$Htilde
  JE <- .ojl_ip_JE(B, st)
  gE <- c(bg$lxE, bg$lcE)
  MF <- .ojl_positive_cholesky_pivoted(JE %*% .ojl_chol_solve(Ht, t(JE)))
  dlE0 <- .ojl_chol_solve(MF, gE + as.numeric(JE %*% .ojl_chol_solve(Ht, st$gtilde)))
  JEt_dlE0 <- as.numeric(crossprod(JE, dlE0))
  dx0 <- .ojl_chol_solve(Ht, JEt_dlE0 - st$gtilde)
  Hpstepx <- as.numeric(.ojl_chol_matrix(Ht) %*% dx0) - JEt_dlE0
  HsteplE <- -as.numeric(JE %*% dx0)
  nrm <- function(v) sqrt(sum(v * v))
  if (isTRUE(nrm(gE) + nrm(st$gtilde) < max(nrm(HsteplE + gE), nrm(Hpstepx + st$gtilde)))) {
    st$s <- numeric(length(st$x))
    st$bstep <- lapply(st$bstep, function(v) numeric(length(v)))
    return(NULL)
  }
  bstep <- .ojl_ip_solve_slack(st$bstep, dx0, B, bs, bg, J, 0)
  as1 <- .ojl_ip_maxstep(.ojl_ip_maxstep(1, bs$slack_x, bstep$slack_x), bs$slack_c, bstep$slack_c)
  al <- .ojl_ip_maxstep(.ojl_ip_maxstep(1, bs$lx, bstep$lx), bs$lc, bstep$lc)
  mm <- max(1, length(bs$slack_x) + length(bs$slack_c))
  muaff <- (sum((bs$slack_x + as1 * bstep$slack_x) * (bs$lx + al * bstep$lx)) +
              sum((bs$slack_c + as1 * bstep$slack_c) * (bs$lc + al * bstep$lc))) / mm
  mumean <- (sum(bs$slack_x * bs$lx) + sum(bs$slack_c * bs$lc)) / mm
  r <- muaff / mumean
  st$munext <- .ojl_nanmax(r * r * r * mumean, mumean / 10)
  musinv <- mu * c(B$sx / bs$slack_x, B$sc / bs$slack_c)
  gtmu <- st$gtilde - as.numeric(crossprod(.ojl_ip_JI(B, st), musinv))
  dlE <- .ojl_chol_solve(MF, gE + as.numeric(JE %*% .ojl_chol_solve(Ht, gtmu)))
  JEt_dlE <- as.numeric(crossprod(JE, dlE))
  dx <- .ojl_chol_solve(Ht, JEt_dlE - gtmu)
  st$s <- dx
  bstep$lxE <- dlE[seq_along(B$eqx)]
  bstep$lcE <- dlE[length(B$eqx) + seq_along(B$eqc)]
  bstep <- .ojl_ip_solve_slack(bstep, dx, B, bs, bg, J, mu)
  st$bstep <- bstep
  Hstepx <- as.numeric(st$H_L_x %*% dx) - JEt_dlE
  HsteplE <- -as.numeric(JE %*% dx)
  c(st$L_x, .ojl_ip_slopealpha(st$s, st$g_L_x, bstep, bg), sum(dx * Hstepx) + sum(dlE * HsteplE))
}

# backtrack_constrained_grad: a backtracking line search on the Lagrangian
# with a condition on its slope
.ojl_backtrack_constrained_grad <- function(phi, alpha, alphamax, qp, c1 = 0.9, c2 = 0.9, rho = 0.5,
                                            alphaminfrac = sqrt(.Machine$double.eps), show_linesearch = FALSE) {
  alpha <- min(alpha, 0.999 * alphamax)
  amin <- alphaminfrac * alpha
  L0 <- qp[1L]; L1 <- qp[2L]; L2 <- qp[3L]
  if (show_linesearch) message("L0 = ", L0, ", L1 = ", L1, ", L2 = ", L2)
  while (isTRUE(alpha >= amin) && alpha > 0) {
    vs <- phi(alpha)
    val <- vs[1L]; sl <- vs[2L]
    dval <- L1 * alpha
    dslope <- L2 * alpha
    if (show_linesearch)
      message("alpha = ", alpha, ", value: (", L0, ", ", val, ", ", L0 + dval, "), slope: (", L1, ", ", sl, ", ",
              L1 + dslope, ")")
    if (is.finite(val) && isTRUE(val - (L0 + dval) <= c1 * abs(val - L0)) &&
        (isTRUE(sl < c2 * abs(L1)) || isTRUE(sl - (L1 + dslope) <= c2 * abs(sl - L1)))) return(alpha)
    alpha <- alpha * rho
  }
  phi(0)
  0
}

# the Lagrangian and its slope along the step (lagrangian_lineslope!),
# which also sets the objective and the Lagrangian of the state
.ojl_ip_lineslope <- function(a, d, cons, B, st) {
  st$x_ls <- st$x + a * st$s
  bls <- Map(function(b, s) b + a * s, st$bstate, st$bstep)
  cls <- cons$c(st$x_ls)
  Jls <- cons$J(st$x_ls)
  r <- d$value_gradient(st$x_ls)
  lf <- .ojl_ip_lagrangian_fg(B, st$x_ls, r[[1L]], r[[2L]], cls, Jls, bls, st$mu)
  st$f_x <- r[[1L]]; st$L_x <- lf$L; st$ev <- lf$ev
  c(lf$L, .ojl_ip_slopealpha(st$s, lf$g, st$bstep, lf$bgrad))
}

# update_state! of IPNewton
.ojl_ip_update <- function(d, cons, B, st, m) {
  st$f_x_previous <- st$f_x
  st$L_x_previous <- st$L_x
  qp <- .ojl_ip_solve_step(st, B)
  if (is.null(qp)) return(invisible(NULL))
  bs <- st$bstate; bst <- st$bstep
  if (.ojl_ip_small(st$x, st$s) && .ojl_ip_small(bs$slack_x, bst$slack_x) &&
      .ojl_ip_small(bs$slack_c, bst$slack_c) && .ojl_ip_small(bs$lx, bst$lx) && .ojl_ip_small(bs$lc, bst$lc))
    return(invisible(NULL))
  amax <- Inf
  for (nm in c("slack_x", "slack_c", "lx", "lc")) amax <- .ojl_ip_maxstep(amax, bs[[nm]], bst[[nm]])
  phi <- function(a) .ojl_ip_lineslope(a, d, cons, B, st)
  st$alpha <- m$linesearch(phi, 1, amax, qp, show_linesearch = isTRUE(m$show_linesearch))
  st$x_previous <- st$x
  st$x <- st$x + st$alpha * st$s
  bs <- Map(function(b, s) b + st$alpha * s, bs, bst)
  p <- st$mu / bs$slack_x
  bs$lx <- pmax(pmin(bs$lx, 1e10 * p), 1e-10 * p)
  p <- st$mu / bs$slack_c
  bs$lc <- pmax(pmin(bs$lc, 1e10 * p), 1e-10 * p)
  st$bstate <- bs
  st$mu <- st$munext
  st$constr_c <- cons$c(st$x)
  st$constr_J <- cons$J(st$x)
  st$ev <- .ojl_ip_ev(B, st$x, st$constr_c, bs)
  invisible(NULL)
}

.ojl_ip_state <- function(d, cons, B, x0, m) {
  n <- length(x0)
  cc <- cons$c(x0)
  if (!isTRUE(all(sign(x0[B$ineqx] - B$bx) == B$sx)) || !isTRUE(all(sign(cc[B$ineqc] - B$bc) == B$sc)))
    warning("Initial guess is not an interior point", call. = FALSE)
  r <- d$value_gradient_hessian(x0)
  J <- cons$J(x0)
  na <- function(k) rep(NaN, k)
  bs <- list(slack_x = abs(B$sx * (x0[B$ineqx] - B$bx)), slack_c = abs(B$sc * (cc[B$ineqc] - B$bc)),
             lx = na(length(B$ineqx)), lc = na(length(B$ineqc)), lxE = na(length(B$eqx)), lcE = na(length(B$eqc)))
  st <- .ojl_state(x = x0, g_x = r[[2L]], H_x = r[[3L]], f_x = r[[1L]], x_previous = na(n), f_x_previous = NaN,
                   s = na(n), mu = 1, munext = 1, L_x = NaN, g_L_x = na(n), H_L_x = matrix(NaN, n, n),
                   L_x_previous = NaN, bstate = bs, bgrad = bs, bstep = bs, constr_c = cc, constr_J = J, ev = NaN,
                   x_ls = na(n), alpha = 1, gtilde = na(n), Htilde = NULL)
  HcI <- cons$h(matrix(0, n, n), x0, .ojl_ip_userlambda(1 / bs$slack_c, B))
  .ojl_ip_init_mu_lambda(B, st, HcI, m$mu0)
  .ojl_ip_fgh(d, cons, B, st)
  st
}

# optimize(d, constraints, x0, ::IPNewton, options) of interior.jl
.ojl_ipnewton_optimize <- function(d, cons, B, x0, m, o) {
  st <- .ojl_ip_state(d, cons, B, x0, m)
  t0 <- proc.time()[["elapsed"]]
  stopped_by_time_limit <- f_limit <- g_limit <- h_limit <- FALSE
  conv <- c(x = FALSE, f = FALSE, g = FALSE, inc = FALSE)
  counter_f_tol <- 0L
  stopped <- !is.finite(st$f_x) || any(!is.finite(st$g_x))
  ninf <- function(v) if (length(v)) max(abs(v)) else 0
  conv[["g"]] <- isTRUE(ninf(st$g_x) + sum(vapply(st$bgrad, ninf, 0)) < o$g_abstol)
  converged <- conv[["g"]]
  iteration <- 0L
  while (!converged && !stopped && iteration < o$iterations) {
    iteration <- iteration + 1L
    .ojl_ip_update(d, cons, B, st, m)
    .ojl_ip_fgh(d, cons, B, st)
    conv <- c(x = isTRUE(max(abs(st$x - st$x_previous)) <= o$x_abstol),
              f = isTRUE(abs(st$L_x - st$L_x_previous) <= o$f_reltol * abs(st$L_x)),
              g = isTRUE(ninf(c(st$g_x, unlist(st$bgrad, use.names = FALSE))) <= o$g_abstol),
              inc = isTRUE(st$L_x > st$L_x_previous))
    counter_f_tol <- if (conv[["f"]]) counter_f_tol + 1L else 0L
    converged <- conv[["x"]] || conv[["g"]] || counter_f_tol > o$successive_f_tol
    stopped_by_time_limit <- isTRUE(proc.time()[["elapsed"]] - t0 > o$time_limit)
    f_limit <- o$f_calls_limit > 0 && d$f_calls >= o$f_calls_limit
    g_limit <- o$g_calls_limit > 0 && d$df_calls >= o$g_calls_limit
    h_limit <- o$h_calls_limit > 0 && d$h_calls >= o$h_calls_limit
    if ((conv[["inc"]] && !o$allow_f_increases) || stopped_by_time_limit || f_limit || g_limit || h_limit)
      stopped <- TRUE
  }
  pick <- conv[["inc"]] && !o$allow_f_increases
  c(list(minimizer = if (pick) st$x_previous else st$x, minimum = if (pick) st$f_x_previous else st$f_x,
         iterations = iteration, g_residual = ninf(st$g_x), f_calls = d$f_calls, g_calls = d$df_calls,
         h_calls = d$h_calls, state = st,
         stopped_by = list(x_converged = conv[["x"]], f_converged = conv[["f"]], g_converged = conv[["g"]],
                           f_increased = conv[["inc"]], iterations = iteration == o$iterations, ls_failed = FALSE)),
    .ojl_changes(st$x, st$x_previous, st$f_x, st$f_x_previous))
}

#' @rdname optim_methods
#' @export
IPNewton <- function(linesearch = NULL, mu0 = "auto", show_linesearch = FALSE) {
  if (!(identical(mu0, "auto") || is.numeric(mu0))) stop("`mu0` must be \"auto\" or a number", call. = FALSE)
  .ojl_alg("IPNewton", linesearch = linesearch %||% .ojl_backtrack_constrained_grad, mu0 = mu0,
           show_linesearch = show_linesearch, hess = TRUE)
}

# the constraints of OptimizationOptimJL's IPNewton path
# (TwiceDifferentiableConstraints of the problem's cons, cons_j and cons_h,
# with central differences when cons_j or cons_h is missing)
.ojl_ipnewton_solve <- function(prob, fs, d, x0, lb, ub, m, o) {
  n <- length(x0)
  fn <- fs$fn; p <- fs$p
  if (is.function(fn$cons)) {
    cfun <- function(x) as.numeric(fn$cons(x, p))
    mc <- length(prob$lcons %||% cfun(x0))
    Jfun <- if (is.function(fn$cons_j)) function(x) matrix(as.numeric(fn$cons_j(x, p)), mc, n)
            else function(x) .fd_jacobian(cfun, x)
    Hs <- if (is.function(fn$cons_h)) function(x) fn$cons_h(x, p)
          else function(x) lapply(seq_len(mc), function(i) .ojl_opt_fd_hessian(function(z) Jfun(z)[i, ], x))
    hfun <- function(H, x, lam) {
      hs <- Hs(x)
      for (i in seq_along(lam)) H <- H + lam[i] * matrix(as.numeric(hs[[i]]), n, n)
      H
    }
    lc <- prob$lcons %||% rep(0, mc)
    uc <- prob$ucons %||% rep(0, mc)
  } else {
    cfun <- function(x) numeric()
    Jfun <- function(x) matrix(0, 0L, n)
    hfun <- function(H, x, lam) H
    lc <- uc <- numeric()
  }
  B <- .ojl_constraint_bounds(lb, ub, lc, uc)
  .ojl_ipnewton_optimize(d, list(c = cfun, J = Jfun, h = hfun), B, x0, m, o)
}

# ---------------------------------------------------------------------------
# The module

#' @rdname optim_methods
#' @export
Optim <- list(
  GradientDescent = GradientDescent, ConjugateGradient = ConjugateGradient,
  AcceleratedGradientDescent = AcceleratedGradientDescent, MomentumGradientDescent = MomentumGradientDescent,
  BFGS = .ojl_BFGS, LBFGS = .ojl_LBFGS, NGMRES = NGMRES, OACCEL = OACCEL, Newton = Newton,
  NewtonTrustRegion = NewtonTrustRegion, IPNewton = IPNewton, Fminbox = Fminbox, ParticleSwarm = ParticleSwarm,
  SimulatedAnnealing = SimulatedAnnealing, SAMIN = SAMIN, NelderMead = .ojl_NelderMead, Adam = .ojl_Adam,
  AdaMax = .ojl_AdaMax, KrylovTrustRegion = .ojl_KrylovTrustRegion, Flat = .ojl_Flat, Sphere = .ojl_Sphere,
  Stiefel = .ojl_Stiefel, PowerManifold = .ojl_PowerManifold, ProductManifold = .ojl_ProductManifold,
  LineSearches = list(HagerZhang = .ojl_HagerZhang, Static = .ojl_Static, BackTracking = .ojl_BackTracking,
                      StrongWolfe = .ojl_StrongWolfe, MoreThuente = .ojl_MoreThuente,
                      InitialStatic = .ojl_InitialStatic, InitialPrevious = .ojl_InitialPrevious,
                      InitialHagerZhang = .ojl_InitialHagerZhang, InitialQuadratic = .ojl_InitialQuadratic,
                      InitialConstantChange = .ojl_InitialConstantChange))
