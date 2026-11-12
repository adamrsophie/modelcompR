# The solvers of SimpleOptimization (lib/SimpleOptimization/src/SimpleOptimization.jl,
# v2.0.2): SimpleLBFGS with the static Strong Wolfe line search of LineSearch.jl
# (src/strong_wolfe.jl, v0.1.20); SimpleBFGS and SimpleNewton, which run the
# SimpleBroyden and SimpleNewtonRaphson iterations of SimpleNonlinearSolve
# (src/broyden.jl and src/raphson.jl, v2) on the gradient; SimpleGradientDescent
# and SimpleSOAP. Julia differentiates the objective with ForwardDiff; here the
# gradient (and for SimpleNewton the Hessian) of the problem is used when it is
# given, else finite differences. The fused multiply adds of the BLAS kernels
# and of Julia's integer power are repeated with .ojl_fma, so that with the same
# derivatives the iterates agree with Julia to the last bit.

.simple_alg <- function(name, ..., bounds = FALSE) {
  a <- .opt_alg(name, ...)
  class(a) <- c(name, "SimpleOptimizationAlgorithm", "OptimizationAlgorithm")
  a$solver <- .simple_solve
  a$traits <- list(bounds = bounds)
  a
}

#' SimpleOptimization solvers
#'
#' Small, low overhead solvers of SimpleOptimization for
#' `solve(prob, SimpleLBFGS())`. Julia differentiates the objective with
#' ForwardDiff; these use the gradient (and Hessian) of the problem when it
#' has them, else finite differences.
#'
#' * `SimpleLBFGS(threshold, linesearch)`: limited memory BFGS with
#'   `threshold` stored pairs and a Strong Wolfe line search
#'   ([StrongWolfeLineSearch()], 50 bracketing and 50 zoom iterations by
#'   default). The only one of these that accepts box constraints (`lb`,
#'   `ub`), by projection. It stops with `"Success"` when the largest
#'   component of the projected gradient is at most
#'   `max(abstol, reltol * initial value)` (`abstol` at least
#'   `sqrt(.Machine$double.eps)`), `"MaxIters"` after `maxiters` iterations
#'   (default 100) and `"Failure"` on a non finite point or a failed line
#'   search.
#' * `SimpleBFGS()`: the Broyden iteration (SimpleBroyden without line
#'   search) on the gradient, which converges to a stationary point;
#'   `abstol` (default `3e-13`) bounds the largest gradient component and
#'   `maxiters` defaults to 100.
#' * `SimpleNewton()`: the Newton Raphson iteration on the gradient, with the
#'   Hessian of the problem (finite differences of the gradient when it has
#'   none). Same tolerances as `SimpleBFGS()`.
#' * `SimpleGradientDescent(eta)`: `u - eta * gradient` until the Euclidean
#'   norm of the gradient is below `abstol` (default `1e-8`); `maxiters`
#'   defaults to 1000.
#' * `SimpleSOAP(...)`: SOAP (Shampoo with Adam in the preconditioner's
#'   eigenbasis). For a matrix `u0` it runs AdamW in the eigenbasis of the
#'   Shampoo preconditioners, refreshed every `freq` iterations by a power
#'   iteration and a QR factorization; for a vector it is AdamW. It always
#'   runs `maxiters` iterations (default 1000) and returns `"MaxIters"`.
#'
#' The statistics follow Julia: `SimpleBFGS()` and `SimpleNewton()` report
#' `maxiters` iterations and no evaluations, `SimpleLBFGS()` the iterations
#' only.
#'
#' @param threshold Number of stored correction pairs.
#' @param linesearch A [StrongWolfeLineSearch()].
#' @param eta Learning rate.
#' @param beta Decay of the first and second moments.
#' @param shampoo_beta Decay of the preconditioners; a negative value uses
#'   `beta[2]`.
#' @param epsilon Constant against division by zero.
#' @param freq Iterations between updates of the eigenbases.
#' @param max_dim Dimensions larger than this keep the identity basis.
#' @param weight_decay Decoupled weight decay, applied as `eta * weight_decay`.
#' @return An optimization algorithm.
#' @examples
#' f <- function(u, p) (1 - u[1])^2 + 100 * (u[2] - u[1]^2)^2
#' g <- function(u, p) c(-2 * (1 - u[1]) - 400 * u[1] * (u[2] - u[1]^2), 200 * (u[2] - u[1]^2))
#' prob <- OptimizationProblem(f, c(-1.2, 1), grad = g)
#' solve(prob, SimpleLBFGS())$u
#' solve(prob, SimpleGradientDescent(eta = 1e-3), maxiters = 200)$objective
#' @name simple_optimization
#' @export
SimpleLBFGS <- function(threshold = 10L, linesearch = StrongWolfeLineSearch(maxiters = 50L, zoom_maxiters = 50L)) {
  threshold <- as.integer(threshold)
  if (length(threshold) != 1L || is.na(threshold) || threshold < 1L)
    stop("`threshold` must be a positive integer", call. = FALSE)
  .simple_alg("SimpleLBFGS", threshold = threshold, linesearch = linesearch, bounds = TRUE)
}
#' @rdname simple_optimization
#' @export
SimpleBFGS <- function() .simple_alg("SimpleBFGS")
#' @rdname simple_optimization
#' @export
SimpleGradientDescent <- function(eta = 0.01) .simple_alg("SimpleGradientDescent", eta = eta)
#' @rdname simple_optimization
#' @export
SimpleNewton <- function() .simple_alg("SimpleNewton")
#' @rdname simple_optimization
#' @export
SimpleSOAP <- function(eta = 3e-3, beta = c(0.95, 0.95), shampoo_beta = -1, epsilon = 1e-8, freq = 10L,
                       max_dim = 10000L, weight_decay = 0.01) {
  if (length(beta) != 2L) stop("`beta` must have two components", call. = FALSE)
  .simple_alg("SimpleSOAP", eta = eta, beta = as.numeric(beta), shampoo_beta = shampoo_beta, epsilon = epsilon,
              freq = as.integer(freq), max_dim = as.integer(max_dim), weight_decay = weight_decay)
}

.simple_solve <- function(prob, alg, opts) {
  fs <- .opt_fns(prob)
  if (is.function(prob$f$cons))
    stop(IncompatibleOptimizerError(paste0("The algorithm ", alg$name, " does not support constraints. Either ",
                                           "remove the `cons` function passed to `OptimizationFunction` or use ",
                                           "a different algorithm.")))
  bounded <- !is.null(prob$lb) || !is.null(prob$ub)
  if (bounded && alg$name != "SimpleLBFGS")
    stop(IncompatibleOptimizerError(paste0("The algorithm ", alg$name, " does not support box constraints. ",
                                           "Either remove the `lb` or `ub` bounds passed to ",
                                           "`OptimizationProblem` or use a different algorithm.")))
  maxiters <- .opt_maxiters(opts$maxiters)
  x0 <- prob$u0
  n <- length(x0)
  t0 <- proc.time()[["elapsed"]]
  G <- function(x) fs$g(x)
  r <- switch(alg$name,
    SimpleLBFGS = {
      ls <- alg$linesearch
      if (!inherits(ls, "StrongWolfeLineSearch"))
        stop("SimpleLBFGS requires a StrongWolfeLineSearch", call. = FALSE)
      lb <- ub <- NULL
      if (bounded) {
        lb <- if (is.null(prob$lb)) rep(-Inf, n) else as.numeric(prob$lb)
        ub <- if (is.null(prob$ub)) rep(Inf, n) else as.numeric(prob$ub)
      }
      .simple_lbfgs(fs$f, G, as.numeric(x0), lb, ub, maxiters %||% 100L, opts$abstol, opts$reltol,
                    ls, alg$threshold)
    },
    SimpleBFGS = , SimpleNewton = {
      tol <- opts$abstol %||% 3e-13
      nl <- if (alg$name == "SimpleBFGS") {
        .simple_broyden(G, as.numeric(x0), maxiters %||% 100L, tol)
      } else {
        H <- if (is.function(fs$fn$hess)) function(x) matrix(as.numeric(fs$fn$hess(x, fs$p)), n, n)
             else function(x) .ojl_opt_fd_hessian(fs$g, x)
        .simple_newton_raphson(G, H, as.numeric(x0), maxiters %||% 100L, tol)
      }
      list(u = nl$u, objective = fs$f(nl$u), retcode = nl$retcode, iterations = maxiters %||% 100L,
           fevals = 0L, original = nl)
    },
    SimpleGradientDescent = .simple_gradient_descent(fs$f, G, as.numeric(x0), alg$eta, maxiters %||% 1000L,
                                                     opts$abstol %||% 1e-8),
    SimpleSOAP = {
      m <- maxiters %||% 1000L
      if (is.matrix(x0)) {
        d <- dim(x0)
        th <- .simple_soap_matrix(function(x) matrix(fs$g(x), d[1L], d[2L]), matrix(as.numeric(x0), d[1L], d[2L]),
                                  alg, m)
        ge <- m + 1L
      } else {
        th <- .simple_soap_vector(G, as.numeric(x0), alg, m)
        ge <- m
      }
      list(u = th, objective = fs$f(th), retcode = "MaxIters", iterations = m, fevals = ge)
    })
  sol <- .opt_solution(prob, alg, r$u, r$objective, r$retcode, r$original)
  .opt_set_fields(sol, stats = OptimizationStats(iterations = r$iterations, time = proc.time()[["elapsed"]] - t0,
                                                fevals = r$fevals %||% 0L, gevals = 0L, hevals = 0L))
}

# ---------------------------------------------------------------------------
# Arithmetic of Julia

# sum of a short vector, accumulated left to right as Julia's sum does for
# fewer than 16 elements
.simple_sum <- function(v) Reduce(`+`, v, 0)

# LinearAlgebra.generic_norm2, the norm of vectors shorter than 32
.simple_norm <- function(v) {
  if (!length(v)) return(0)
  m <- max(abs(v))
  if (is.na(m) || m == 0 || is.infinite(m)) return(m)
  if (is.finite(length(v) * m * m) && m * m != 0) return(sqrt(.simple_sum(v * v)))
  m * sqrt(.simple_sum((v / m) * (v / m)))
}

# NonlinearSolveBase.L2_NORM: a simd sum of squares, a chain of fused
# multiply adds for short vectors
.simple_l2 <- function(v) sqrt(.ojl_dot(v, v))

# transpose(A) * x of OpenBLAS gemv for small matrices: a fused dot product
# per column
.simple_gemv_t <- function(A, x) vapply(seq_len(ncol(A)), function(j) .ojl_dot(A[, j], x), numeric(1))

# x^n of Julia for a Float64 and an integer (Base.Math.pow_body): the
# compensated power by squaring, reliable for |n| below 2^20; beyond 3 * 2^13
# Julia uses its exp and log based power, approximated by R's
.simple_two_mul <- function(a, b) {
  p <- a * b
  c(p, .ojl_fma(a, b, -p))
}
.simple_powi <- function(x, n) {
  if (n == 0) return(1)
  if (n < -4096 || n > 24576) return(x^n)
  if (n == 3) return(x * x * x)
  y <- 1
  xnlo <- -0
  ynlo <- 0
  if (n < 0) {
    rx <- 1 / x
    if (n == -2) return(rx * rx)
    if (is.finite(x)) xnlo <- -.ojl_fma(x, rx, -1) * rx
    x <- rx
    n <- -n
  }
  while (n > 1) {
    if (n %% 2 == 1) {
      err <- .ojl_fma(y, xnlo, x * ynlo)
      ym <- .simple_two_mul(x, y)
      y <- ym[1L]
      ynlo <- ym[2L] + err
    }
    err <- x * 2 * xnlo
    xm <- .simple_two_mul(x, x)
    x <- xm[1L]
    xnlo <- xm[2L] + err
    n <- n %/% 2
  }
  err <- .ojl_fma(y, xnlo, x * ynlo)
  if (is.finite(x) && is.finite(err) && err != 0) .ojl_fma(x, y, err) else x * y
}

# ---------------------------------------------------------------------------
# SimpleLBFGS and the static Strong Wolfe line search of LineSearch.jl

.simple_clamp <- function(x, lo, hi) if (isTRUE(x > hi)) hi else if (isTRUE(x < lo)) lo else x

# _sw_interpolate: the minimizer of the cubic through both ends of the bracket
.simple_sw_interpolate <- function(alo, ahi, plo, phi, dlo, dhi) {
  d1 <- dlo + dhi - 3 * (plo - phi) / (alo - ahi)
  desc <- d1 * d1 - dlo * dhi
  d2 <- sqrt(max(desc, 0))
  cand <- ahi - (ahi - alo) * ((dhi + d2 - d1) / (dhi - dlo + 2 * d2))
  if (!is.finite(cand)) cand <- (alo + ahi) / 2
  if (isTRUE(desc < 0)) (alo + ahi) / 2 else cand
}

# _sw_zoom (Nocedal and Wright, algorithm 3.6)
.simple_sw_zoom <- function(ev, alo, ahi, p0, d0, plo, dlo, phi, dhi, c1, c2, maxiters) {
  for (k in seq_len(maxiters)) {
    aj <- .simple_sw_interpolate(alo, ahi, plo, phi, dlo, dhi)
    br <- 0.01 * abs(ahi - alo)
    aj <- .simple_clamp(aj, min(alo, ahi) + br, max(alo, ahi) - br)
    e <- ev(aj)
    if (isTRUE(e[1L] > p0 + c1 * aj * d0) || isTRUE(e[1L] >= plo)) {
      ahi <- aj
      phi <- e[1L]
      dhi <- e[2L]
    } else if (isTRUE(abs(e[2L]) <= -c2 * d0)) {
      return(list(alpha = aj, ok = TRUE))
    } else {
      if (isTRUE(e[2L] * (ahi - alo) >= 0)) {
        ahi <- alo
        phi <- plo
        dhi <- dlo
      }
      alo <- aj
      plo <- e[1L]
      dlo <- e[2L]
    }
  }
  list(alpha = alo, ok = FALSE)
}

# _sw_search (Nocedal and Wright, algorithm 3.5)
.simple_sw_search <- function(ev, p0, d0, c1, c2, ainit, amax, maxiters, zmax) {
  if (d0 >= 0) return(list(alpha = 0, ok = FALSE))
  aprev <- 0
  ai <- ainit
  pprev <- p0
  dprev <- d0
  for (i in seq_len(maxiters)) {
    e <- ev(ai)
    if (isTRUE(e[1L] > p0 + c1 * ai * d0) || (isTRUE(e[1L] >= pprev) && i > 1L))
      return(.simple_sw_zoom(ev, aprev, ai, p0, d0, pprev, dprev, e[1L], e[2L], c1, c2, zmax))
    if (isTRUE(abs(e[2L]) <= -c2 * d0)) return(list(alpha = ai, ok = TRUE))
    if (isTRUE(e[2L] >= 0))
      return(.simple_sw_zoom(ev, ai, aprev, p0, d0, e[1L], e[2L], pprev, dprev, c1, c2, zmax))
    # blocked at the cap while still descending: accepted on Armijo alone
    if (ai >= amax) return(list(alpha = ai, ok = TRUE))
    aprev <- ai
    pprev <- e[1L]
    dprev <- e[2L]
    ai <- min(ai * 2, amax)
  }
  list(alpha = ai, ok = FALSE)
}

# _objective_line_search: the step length along d, NULL when it fails
.simple_sw_step <- function(f, G, x, d, lb, ub, ls) {
  amax <- ls$alpha_max
  if (!is.null(lb)) {
    for (i in seq_along(x)) {
      if (d[i] > 0) amax <- min(amax, (ub[i] - x[i]) / d[i])
      else if (d[i] < 0) amax <- min(amax, (lb[i] - x[i]) / d[i])
    }
    amax <- max(amax, 0)
  }
  if (isTRUE(amax <= 0)) return(NULL)
  amaxv <- max(0, min(ls$alpha_max, amax))
  ainit <- min(ls$alpha_init, amaxv)
  if (!isTRUE(amaxv > 0)) return(NULL)
  ev <- function(a) {
    u <- x + a * d
    c(f(u), .ojl_dot(G(u), d))
  }
  e0 <- ev(0)
  if (!(is.finite(e0[1L]) && is.finite(e0[2L]) && e0[2L] < 0)) return(NULL)
  r <- .simple_sw_search(ev, e0[1L], e0[2L], ls$c1, ls$c2, ainit, amaxv, ls$maxiters, ls$zoom_maxiters)
  if (r$ok) r$alpha else NULL
}

# _lbfgs_direction: the two loop recursion over the stored pairs
.simple_lbfgs_direction <- function(g, S, Y, pit, M) {
  al <- numeric(M)
  q <- g
  lower <- pit - M
  upper <- pit - 1L
  for (index in upper:lower) {
    if (index < 1L) next
    j <- (index - 1L) %% M + 1L
    rho <- 1 / .ojl_dot(Y[[j]], S[[j]])
    al[j] <- rho * .ojl_dot(S[[j]], q)
    q <- q - al[j] * Y[[j]]
  }
  r <- q
  if (pit > 1L) {
    j <- (upper - 1L) %% M + 1L
    r <- (.ojl_dot(S[[j]], Y[[j]]) / .simple_sum(Y[[j]] * Y[[j]])) * q
  }
  for (index in lower:upper) {
    if (index < 1L) next
    j <- (index - 1L) %% M + 1L
    rho <- 1 / .ojl_dot(Y[[j]], S[[j]])
    r <- r + S[[j]] * (al[j] - rho * .ojl_dot(Y[[j]], r))
  }
  -r
}

.simple_lbfgs <- function(f, G, x0, lb, ub, maxiters, abstol, reltol, ls, M) {
  bounded <- !is.null(lb)
  proj <- function(x) if (bounded) pmin(pmax(x, lb), ub) else x
  pgrad <- function(x, g) if (bounded) x - proj(x - g) else g
  feas <- function(x, d) {
    if (bounded) d[which((x <= lb & d < 0) | (x >= ub & d > 0))] <- 0
    d
  }
  eps <- .Machine$double.eps
  x <- proj(x0)
  fx <- f(x)
  g <- G(x)
  dtol <- sqrt(eps)
  abs_tol <- if (is.null(abstol)) dtol else max(abstol, dtol)
  rel_tol <- if (is.null(reltol)) 0 else max(reltol, 0)
  tol <- max(abs_tol, rel_tol * max(abs(pgrad(x, g))))
  S <- Y <- rep(list(numeric(length(x))), M)
  pit <- 0L
  retcode <- "MaxIters"
  iters <- 0L
  done <- function() {
    if (!is.finite(fx) || !all(is.finite(g))) return("Failure")
    if (isTRUE(max(abs(pgrad(x, g))) <= tol)) return("Success")
    NULL
  }
  for (i in seq_len(maxiters)) {
    rc <- done()
    if (!is.null(rc)) {
      retcode <- rc
      break
    }
    iters <- i
    pit <- pit + 1L
    d <- feas(x, .simple_lbfgs_direction(g, S, Y, pit, M))
    if (isTRUE(.ojl_dot(d, g) >= 0)) {
      pit <- 1L
      d <- feas(x, -g)
    }
    if (isTRUE(max(abs(d)) <= eps)) {
      retcode <- "Failure"
      break
    }
    a <- .simple_sw_step(f, G, x, d, lb, ub, ls)
    if (is.null(a)) {
      retcode <- "Failure"
      break
    }
    xn <- proj(x + a * d)
    fxn <- f(xn)
    gn <- G(xn)
    s <- xn - x
    y <- gn - g
    curv <- .ojl_dot(s, y)
    if (!is.finite(curv) || curv <= dtol * max(.simple_norm(s) * .simple_norm(y), eps)) {
      pit <- 0L
    } else {
      j <- (pit - 1L) %% M + 1L
      S[[j]] <- s
      Y[[j]] <- y
    }
    x <- xn
    fx <- fxn
    g <- gn
  }
  if (retcode == "MaxIters") retcode <- done() %||% "MaxIters"
  list(u = x, objective = fx, retcode = retcode, iterations = iters, fevals = 0L)
}

# ---------------------------------------------------------------------------
# SimpleBroyden and SimpleNewtonRaphson of SimpleNonlinearSolve on the
# gradient, with their default termination (the largest component of the
# residual at most abstol)

.simple_broyden <- function(G, x, maxiters, abstol) {
  fx <- G(x)
  if (all(fx == 0)) return(list(u = x, resid = fx, retcode = "Success"))
  fnorm <- .simple_l2(fx)
  a0 <- if (fnorm >= 1e-5) max(.simple_l2(x), 1) / (2 * fnorm) else 1
  Ji <- diag(a0, length(x))
  xo <- x
  fprev <- fx
  for (k in seq_len(maxiters)) {
    dx <- -.ojl_gemv(Ji, fprev)
    x <- xo + dx
    fx <- G(x)
    if (isTRUE(max(abs(fx)) <= abstol)) return(list(u = x, resid = fx, retcode = "Success"))
    Jdf <- .ojl_gemv(Ji, fx - fprev)
    d <- .ojl_dot(dx, Jdf)
    Ji <- Ji + outer((dx - Jdf) / d, .simple_gemv_t(Ji, dx))
    xo <- x
    fprev <- fx
  }
  list(u = x, resid = fx, retcode = "MaxIters")
}

# the residual tested is the one before the step, as in SimpleNewtonRaphson
.simple_newton_raphson <- function(G, H, x, maxiters, abstol) {
  fx <- G(x)
  if (all(fx == 0)) return(list(u = x, resid = fx, retcode = "Success"))
  J <- H(x)
  for (k in seq_len(maxiters)) {
    x <- x - .ojl_backslash(J, fx)
    if (isTRUE(max(abs(fx)) <= abstol)) return(list(u = x, resid = fx, retcode = "Success"))
    fx <- G(x)
    J <- H(x)
  }
  list(u = x, resid = fx, retcode = "MaxIters")
}

# ---------------------------------------------------------------------------
# SimpleGradientDescent and SimpleSOAP

.simple_gradient_descent <- function(f, G, theta, eta, maxiters, abstol) {
  g <- G(theta)
  retcode <- "MaxIters"
  iters <- maxiters
  for (i in seq_len(maxiters)) {
    theta <- theta - eta * g
    g <- G(theta)
    if (isTRUE(.simple_norm(g) < abstol)) {
      retcode <- "Success"
      iters <- i
      break
    }
  }
  list(u = theta, objective = f(theta), retcode = retcode, iterations = iters, fevals = iters + 1L)
}

# AdamW, the vector path of SimpleSOAP
.simple_soap_vector <- function(G, theta, alg, maxiters) {
  b1 <- alg$beta[1L]
  b2 <- alg$beta[2L]
  eta <- alg$eta
  wd <- alg$weight_decay
  ea <- eas <- numeric(length(theta))
  for (t in seq_len(maxiters)) {
    g <- G(theta)
    ea <- b1 * ea + (1 - b1) * g
    eas <- b2 * eas + (1 - b2) * (g * g)
    s <- eta * sqrt(1 - .simple_powi(b2, t)) / (1 - .simple_powi(b1, t))
    theta <- theta - s * ea / (sqrt(eas) + alg$epsilon) - eta * wd * theta
  }
  theta
}

# the matrix path of SimpleSOAP: Adam in the eigenbases QL, QR of the
# left and right Shampoo preconditioners L = E[G G'] and R = E[G' G]
.simple_soap_matrix <- function(G, theta, alg, maxiters) {
  m <- nrow(theta)
  n <- ncol(theta)
  b1 <- alg$beta[1L]
  b2 <- alg$beta[2L]
  sb <- if (alg$shampoo_beta >= 0) alg$shampoo_beta else b2
  eta <- alg$eta
  wd <- alg$weight_decay
  ea <- eas <- matrix(0, m, n)
  uL <- m <= alg$max_dim
  uR <- n <= alg$max_dim
  L <- matrix(0, m, m)
  R <- matrix(0, n, n)
  QL <- diag(m)
  QR <- diag(n)
  fwd <- function(X) {
    if (uL) X <- crossprod(QL, X)
    if (uR) X <- X %*% QR
    X
  }
  bwd <- function(X) {
    if (uL) X <- QL %*% X
    if (uR) X <- tcrossprod(X, QR)
    X
  }
  accum <- function(g) {
    if (uL) L <<- (1 - sb) * tcrossprod(g) + sb * L
    if (uR) R <<- (1 - sb) * crossprod(g) + sb * R
  }
  # eigenvectors by decreasing eigenvalue
  eigh <- function(P) eigen((P + t(P)) / 2 + 0x1.4484bfeebc2a0p-100 * diag(nrow(P)), symmetric = TRUE)$vectors
  # one power iteration and a QR factorization; the second moment follows
  # the new order of the eigenvectors
  pqr <- function(Q, GG, row) {
    perm <- order(diag(crossprod(Q, GG %*% Q)), decreasing = TRUE)
    Qn <- qr.Q(qr(GG %*% Q[, perm, drop = FALSE], tol = 0))
    eas <<- if (row) eas[perm, , drop = FALSE] else eas[, perm, drop = FALSE]
    Qn
  }
  accum(G(theta))
  if (uL) QL <- eigh(L)
  if (uR) QR <- eigh(R)
  for (t in seq_len(maxiters)) {
    g <- G(theta)
    buf <- fwd(g)
    ea <- b1 * ea + (1 - b1) * buf
    eas <- b2 * eas + (1 - b2) * (buf * buf)
    buf <- bwd(ea / (sqrt(eas) + alg$epsilon))
    s <- eta * sqrt(1 - .simple_powi(b2, t)) / (1 - .simple_powi(b1, t))
    theta <- theta - s * buf - eta * wd * theta
    ea_orig <- bwd(ea)
    accum(g)
    if (t %% alg$freq == 0L) {
      if (uL) QL <- pqr(QL, L, TRUE)
      if (uR) QR <- pqr(QR, R, FALSE)
    }
    ea <- fwd(ea_orig)
  }
  theta
}
