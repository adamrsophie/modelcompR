# Runge-Kutta tableau tools of DiffEqDevTools.jl (src/tableau_info.jl:
# stability_region, imaginary_stability_interval, residual_order_condition,
# check_tableau; src/ode_tableaus.jl: deduce_Butcher_tableau; the tableau
# plot recipe of src/plotrecipes.jl), with the rooted tree machinery of
# RootedTrees.jl they use (RootedTreeIterator, SubtreeIterator, symmetry,
# density, derivative_weight, elementary_weight, residual_order_condition of
# src/RootedTrees.jl and src/time_integration_methods.jl) and the scalar
# SimpleTrustRegion of SimpleNonlinearSolve.jl (src/trust_region.jl) with
# AutoFiniteDiff derivatives (FiniteDiff forward differences), which the
# stability interval searches run.
#
# A tableau is any list with fields A, c, alpha (the weights of the method),
# alphaEEst (those of the embedded method), order and adaptiveorder, such as
# ExplicitRKTableau() (R/solvers-erk-tableau.R).

.dt_is_tableau <- function(x) is.list(x) && !is.null(x$A) && !is.null(x$alpha) && !inherits(x, "ODEAlgorithm")

# ---------------------------------------------------------------------------
# Rooted trees as level sequences (integer vectors)

# RootedTreeIterator(order): the canonical level sequences of all rooted
# trees with `order` vertices, in the order of the iterator.
.rt_trees <- function(order) {
  ls <- seq_len(order)
  out <- list(ls)
  if (order <= 0L) return(out)
  repeat {
    two <- ls[1L] + 1L
    p <- 1L
    for (i in seq_along(ls)[-1L]) if (ls[i] > two) p <- i
    if (p == 1L) break
    level_q <- ls[p] - 1L
    q <- 1L
    for (i in seq_len(p)) if (ls[i] == level_q) q <- i
    for (i in p:length(ls)) ls[i] <- ls[i - (p - q)]
    out[[length(out) + 1L]] <- ls
  }
  out
}

# SubtreeIterator: the child subtrees of a tree, as level sequences.
.rt_subtrees <- function(ls) {
  out <- list()
  n <- length(ls)
  i <- 2L
  while (i <= n) {
    j <- i
    while (j < n && ls[j + 1L] > ls[i]) j <- j + 1L
    out[[length(out) + 1L]] <- ls[i:j]
    i <- j + 1L
  }
  out
}

.rt_equal <- function(a, b) length(a) == length(b) && all(a - a[1L] == b - b[1L])

# symmetry (sigma) of a canonical tree: identical subtrees are adjacent.
.rt_symmetry <- function(ls) {
  if (length(ls) <= 2L) return(1)
  sub <- .rt_subtrees(ls)
  result <- 1
  same <- 1
  prev <- sub[[1L]]
  for (s in sub[-1L]) {
    if (.rt_equal(s, prev)) {
      same <- same + 1
    } else {
      result <- result * factorial(same) * .rt_symmetry(prev)^same
      same <- 1
    }
    prev <- s
  }
  result * factorial(same) * .rt_symmetry(prev)^same
}

# density (gamma): the order times the densities of the subtrees.
.rt_density <- function(ls) {
  if (!length(ls)) return(1)
  Reduce(`*`, lapply(.rt_subtrees(ls), .rt_density), length(ls))
}

.rt_derivative_weight <- function(ls, A, c) {
  result <- 0 * c + 1
  for (s in .rt_subtrees(ls)) result <- result * as.numeric(A %*% .rt_derivative_weight(s, A, c))
  result
}

# (Phi(t) - 1 / gamma(t)) / sigma(t)
.rt_residual <- function(ls, A, b, c)
  (sum(b * .rt_derivative_weight(ls, A, c)) - 1 / .rt_density(ls)) * (1 / .rt_symmetry(ls))

# ---------------------------------------------------------------------------
# Order conditions

#' Order conditions of Runge-Kutta tableaus
#'
#' `residual_order_condition(tab, order)` evaluates the residuals
#' \eqn{(\Phi(t) - 1/\gamma(t)) / \sigma(t)} of the order conditions of a
#' Runge-Kutta tableau for every rooted tree `t` with `order` vertices (in
#' the order of RootedTrees.jl's tree iterator), with the elementary weight
#' \eqn{\Phi}, the density \eqn{\gamma} and the symmetry \eqn{\sigma} of the
#' tree. With a `reducer` (a function of two arguments such as `` `+` ``) the
#' residuals are first transformed by `mapper` and then combined from left to
#' right. `embedded = TRUE` uses the weights of the embedded method.
#'
#' `check_tableau(tab)` checks the conditions of every order up to the
#' declared order of the method (and of the embedded method, when it has
#' one): the sum of the absolute residuals of each order must be below `tol`.
#' It returns `TRUE` or signals an error.
#'
#' @param tab A Runge-Kutta tableau, such as [ExplicitRKTableau()].
#' @param order Number of vertices of the rooted trees.
#' @param reducer `NULL` for the vector of residuals, or a function combining
#'   two values.
#' @param mapper Function applied to each residual before reducing.
#' @param embedded Use the embedded method.
#' @param tol Tolerance on the sum of absolute residuals.
#' @return A numeric vector or value; `check_tableau()` returns `TRUE`.
#' @examples
#' tab <- constructTsit5ExplicitRK()
#' residual_order_condition(tab, 3)
#' residual_order_condition(tab, 5, `+`, abs)
#' # the order 5 conditions of Tsit5 hold only approximately, as in Julia
#' tryCatch(check_tableau(tab), error = conditionMessage)
#' @export
residual_order_condition <- function(tab, order, reducer = NULL, mapper = function(x) x^2, embedded = FALSE) {
  A <- as.matrix(tab$A); cc <- as.numeric(tab$c)
  b <- as.numeric(if (embedded) tab$alphaEEst else tab$alpha)
  res <- vapply(.rt_trees(as.integer(order)), .rt_residual, 0, A = A, b = b, c = cc)
  if (is.null(reducer)) return(res)
  Reduce(match.fun(reducer), lapply(res, match.fun(mapper)))
}

#' @rdname residual_order_condition
#' @export
check_tableau <- function(tab, tol = 10 * .Machine$double.eps) {
  ok <- all(vapply(seq_len(tab$order), function(i) residual_order_condition(tab, i, `+`, abs) < tol, TRUE))
  if (!ok) stop("Tableau's order is not correct.", call. = FALSE)
  if (!is.null(tab$adaptiveorder) && tab$adaptiveorder != 0) {
    # every order up to the embedded order (Julia iterates over the integer
    # adaptiveorder itself, which checks that order only)
    ok <- all(vapply(seq_len(tab$adaptiveorder), function(i)
      residual_order_condition(tab, i, `+`, abs, embedded = TRUE) < tol, TRUE))
    if (!ok) stop("Tableau's embedded order is not correct.", call. = FALSE)
  }
  TRUE
}

# ---------------------------------------------------------------------------
# Stability

# (I - z A) \ e: forward substitution for lower triangular A, as Julia's `\`
# takes the triangular path (BLAS trsv, column oriented); LU otherwise.
.dt_stage_solve <- function(M) {
  s <- nrow(M)
  e <- rep(1, s)
  # a non finite z makes z * 0 = NaN above the diagonal, so Julia takes the
  # LU path, whose LAPACK call rejects the matrix
  if (!all(is.finite(M))) stop("matrix contains Infs or NaNs", call. = FALSE)
  if (all(M[upper.tri(M)] == 0)) {
    x <- if (is.complex(M)) complex(real = e) else e
    for (j in seq_len(s)) {
      x[j] <- x[j] / M[j, j]
      if (j < s) x[(j + 1L):s] <- x[(j + 1L):s] - x[j] * M[(j + 1L):s, j]
    }
    return(x)
  }
  solve(M, e)
}

.dt_tableau_r <- function(z, tab, embedded = FALSE) {
  A <- as.matrix(tab$A)
  b <- as.numeric(if (embedded) tab$alphaEEst else tab$alpha)
  stages <- .dt_stage_solve(diag(nrow(A)) - z * A)
  1 + z * sum(b * stages)
}

# One step of size 1 of u' = z u from u = 1. A complex z runs the real
# system of the real and imaginary parts, whose arithmetic is that of the
# complex product.
.dt_alg_r <- function(z, alg) {
  if (is.complex(z) && Im(z) != 0) {
    a <- Re(z); b <- Im(z)
    prob <- ODEProblem(function(u, p, t) c(a * u[1L] - b * u[2L], a * u[2L] + b * u[1L]), c(1, 0), c(0, 1))
    u <- solve(prob, alg, dt = 1, adaptive = FALSE)$u
    u <- u[[length(u)]]
    return(complex(real = u[1L], imaginary = u[2L]))
  }
  z <- Re(z)
  prob <- ODEProblem(function(u, p, t) z * u, 1, c(0, 1))
  u <- solve(prob, alg, dt = 1, adaptive = FALSE)$u
  u[[length(u)]][[1L]]
}

# SimpleTrustRegion of SimpleNonlinearSolve for a scalar problem f(x) = 0,
# with AutoFiniteDiff (forward difference) derivatives and the default
# absolute termination |f(x)| <= eps^(4/5).
.dt_scalar_trust_region <- function(f, x, abstol = NULL, maxiters = 1000L) {
  tol <- abstol %||% .Machine$double.eps^0.8
  h0 <- sqrt(.Machine$double.eps)
  deriv <- function(x) {
    h <- max(h0 * abs(x), h0)
    (f(x + h) - f(x)) / h
  }
  dogleg <- function(J, fx, g, D) {
    dN <- -(fx / J)
    if (abs(dN) <= D) return(dN)
    dsd <- -g
    nsd <- abs(dsd)
    if (nsd >= D) return(dsd * (D / nsd))
    dNsd <- dN - dsd
    a <- dNsd * dNsd
    bb <- dsd * dNsd
    fact <- bb^2 - a * (dsd * dsd - D^2)
    tau <- (-bb + sqrt(fact)) / a
    dsd + tau * dNsd
  }
  fx <- f(x)
  xo <- x
  J <- deriv(x)
  Dmax <- max(abs(fx), 0)
  D <- Dmax / 11
  fk <- abs(fx)^2 / 2
  H <- J * J
  g <- J * fx
  if (abs(fx) <= tol) return(x)
  shrink <- 0L
  for (it in seq_len(maxiters)) {
    d <- dogleg(J, fx, g, D)
    x <- xo + d
    fx <- f(x)
    fk1 <- abs(fx)^2 / 2
    r <- (fk1 - fk) / (d * g + (d * (H * d)) / 2)
    if (r >= 0.25) {
      shrink <- 0L
    } else {
      D <- 0.25 * D
      shrink <- shrink + 1L
      if (shrink > 32L) return(x)
    }
    if (r >= 1e-4) {
      if (abs(fx) <= tol) return(x)
      xo <- x
      J <- deriv(x)
      fx <- f(x)
      if (r > 0.75) D <- min(2 * D, Dmax)
      fk <- fk1
      H <- J * J
      g <- J * fx
    }
  }
  x
}

#' Stability regions of Runge-Kutta methods
#'
#' `stability_region(z, tab)` evaluates the stability function
#' \eqn{r(z) = 1 + z b^T (I - zA)^{-1} e} of a Runge-Kutta tableau at a real
#' or complex `z` (`embedded = TRUE` for the embedded method); the method is
#' stable where \eqn{|r(z)| < 1}. `stability_region(z, alg)` computes it for
#' an ODE algorithm by taking one step of size 1 of \eqn{u' = zu} from
#' \eqn{u = 1}.
#'
#' `stability_region(tab)` and `stability_region(alg)` (one argument) return
#' the length of the stability region on the negative real axis: the root of
#' \eqn{|r(x)| = 1} found from `initial_guess` by a trust region method with
#' finite difference derivatives (SimpleNonlinearSolve's
#' `SimpleTrustRegion(autodiff = AutoFiniteDiff())`).
#' `imaginary_stability_interval()` returns the length of the stability
#' interval on the imaginary axis in the same way, from `initial_guess` (by
#' default the number of stages minus one for a tableau, 20 for an
#' algorithm).
#'
#' @param z A real or complex number, or (alone) a tableau or algorithm.
#' @param tab A Runge-Kutta tableau ([ExplicitRKTableau()]) or an ODE
#'   algorithm.
#' @param embedded Use the embedded method of the tableau.
#' @param initial_guess Start of the root search.
#' @param abstol,maxiters Options of the root search.
#' @param ... Unused.
#' @return A number (complex for complex `z`).
#' @examples
#' tab <- constructTsit5ExplicitRK()
#' stability_region(-1, tab)
#' stability_region(1i, RK4())
#' stability_region(tab)
#' imaginary_stability_interval(RK4())
#' @export
stability_region <- function(z, tab, ..., embedded = FALSE, initial_guess = -3, abstol = NULL, maxiters = 1000L) {
  if (missing(tab)) {
    tab <- z
    f <- function(x) abs(stability_region(x, tab)) - 1
    return(.dt_scalar_trust_region(f, initial_guess, abstol, maxiters))
  }
  if (.dt_is_tableau(tab)) return(.dt_tableau_r(z, tab, embedded))
  if (inherits(tab, "ODEAlgorithm")) return(.dt_alg_r(z, tab))
  stop("stability_region() needs a Runge-Kutta tableau or an ODE algorithm", call. = FALSE)
}

#' @rdname stability_region
#' @export
imaginary_stability_interval <- function(tab, initial_guess = NULL, abstol = NULL, maxiters = 1000L) {
  initial_guess <- initial_guess %||% if (.dt_is_tableau(tab)) nrow(as.matrix(tab$A)) - 1 else 20
  f <- function(x) abs(stability_region(complex(real = 0, imaginary = x), tab)) - 1
  .dt_scalar_trust_region(f, initial_guess, abstol, maxiters)
}

# ---------------------------------------------------------------------------
# Deducing a tableau from an explicit method

#' Deduce the Butcher tableau of an explicit Runge-Kutta method
#'
#' `deduce_Butcher_tableau(erk)` recovers the coefficients `A`, `b` and `c`
#' of an explicit Runge-Kutta method by taking one step of size 1 with it on
#' special problems: the number of stages is the number of function
#' evaluations of the step (including the evaluation at the start and, for
#' methods with the first same as last property, the one at the end), the
#' nodes are the times of the evaluations, the weights are the step of a
#' problem whose `i`th evaluation returns the `i`th unit vector, and row `i`
#' of `A` is the state seen by the `i`th evaluation of that problem.
#' Computations are in double precision.
#'
#' @param erk An explicit ODE algorithm.
#' @return A list with the matrix `A` and the vectors `b` and `c`.
#' @examples
#' deduce_Butcher_tableau(RK4())
#' @export
deduce_Butcher_tableau <- function(erk) {
  one_step <- function(f, u0) {
    u <- solve(ODEProblem(f, u0, c(0, 1)), erk, dt = 1, adaptive = FALSE)$u
    u[[length(u)]]
  }
  n <- 0L
  one_step(function(u, p, t) {
    n <<- n + 1L
    0 * u
  }, 0)
  s <- n
  cc <- numeric(s)
  n <- 0L
  one_step(function(u, p, t) {
    n <<- n + 1L
    if (n <= s) cc[n] <<- t
    0 * u
  }, 0)
  n <- 0L
  b <- one_step(function(u, p, t) {
    n <<- n + 1L
    as.numeric(seq_along(u) == n)
  }, numeric(s))
  A <- matrix(0, s, s)
  n <- 0L
  one_step(function(u, p, t) {
    n <<- n + 1L
    A[n, ] <<- u
    as.numeric(seq_along(u) == n)
  }, numeric(s))
  list(A = A, b = as.numeric(b), c = cc)
}

# ---------------------------------------------------------------------------
# Plot of the stability region

#' Plot the stability region of a tableau
#'
#' Draws the region \eqn{|r(z)| < 1} of the stability function of a
#' Runge-Kutta tableau in the complex plane with base graphics (the
#' DiffEqDevTools plot recipe). With `order_star = TRUE` it draws the order
#' star \eqn{|r(z) / e^z| < 1} instead.
#'
#' @param x A Runge-Kutta tableau.
#' @param dx,dy Grid spacing.
#' @param xlim,ylim Plot ranges of the real and imaginary parts.
#' @param order_star Plot the order star.
#' @param embedded Use the embedded method.
#' @param col Fill colors of the outside and inside of the region.
#' @param ... Passed to [graphics::image()].
#' @return Invisibly `x`.
#' @examples
#' plot(constructTsit5ExplicitRK(), dx = 0.1, dy = 0.1)
#' @export
plot.Tableau <- function(x, dx = 1 / 100, dy = 1 / 100, xlim = c(-6, 1), ylim = c(-5, 5), order_star = FALSE,
                         embedded = FALSE, col = c("white", "gray60"), ...) {
  xs <- seq(xlim[1L], xlim[2L], by = dx)
  ys <- seq(ylim[1L], ylim[2L], by = dy)
  inside <- outer(xs, ys, Vectorize(function(u, v) {
    z <- complex(real = u, imaginary = v)
    r <- .dt_tableau_r(z, x, embedded)
    abs(if (order_star) r / exp(z) else r) < 1
  }))
  graphics::image(xs, ys, inside + 0, col = col, asp = 1, xlab = "Re(z)", ylab = "Im(z)", ...)
  graphics::contour(xs, ys, inside + 0, levels = 0.5, drawlabels = FALSE, add = TRUE)
  invisible(x)
}
