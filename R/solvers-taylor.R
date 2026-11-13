# Taylor series methods of OrdinaryDiffEqTaylorSeries (algorithms.jl,
# alg_utils.jl, TaylorSeries_caches.jl, TaylorSeries_perform_step.jl,
# interpolants.jl) with the Taylor mode arithmetic of TaylorDiff
# (primitive.jl, utils.jl).
#
# Julia builds the Taylor jet of ExplicitTaylor and
# ExplicitTaylorAdaptiveOrder once, symbolically (Symbolics and TaylorDiff,
# then build_function with common subexpression elimination). Here the jet
# is computed at every step by evaluating the user's f on truncated Taylor
# polynomials (class mcr_taylor below) with TaylorDiff's recurrences; the
# coefficients agree with Julia's up to the rounding of the different
# operation order. ExplicitTaylor2 evaluates f on first order Taylor numbers
# in Julia too, and matches it bit for bit.

# ---------------------------------------------------------------------------
# Taylor mode arithmetic. A traced value is a list (class mcr_taylor) whose
# elements are numeric vectors: a number (length 1) or the coefficients
# (value, then the normalized derivatives) of a truncated Taylor polynomial
# (length K + 1, K >= 1). Indexing, c(), arithmetic, math and summary
# methods return mcr_taylor values, so a right-hand side written for
# numeric vectors (u[1], c(...), sum(u)) computes the Taylor expansion.
# Comparisons use the values, as TaylorDiff.

.taylor_vec <- function(xs, nms = NULL) structure(xs, names = nms, class = "mcr_taylor")

# Elements of a value as a list of numeric vectors.
.taylor_elems <- function(z) {
  if (is.null(z)) return(list())
  if (inherits(z, "mcr_taylor")) return(unclass(z))
  if (is.numeric(z) || is.logical(z)) return(as.list(as.numeric(z)))
  if (is.list(z)) return(do.call(c, c(list(list()), lapply(unname(z), .taylor_elems))))
  stop("cannot evaluate a value of class ", class(z)[1L], " in Taylor mode", call. = FALSE)
}

#' @export
`[.mcr_taylor` <- function(x, i, ...) {
  if (missing(i)) return(x)
  out <- unclass(x)[i]
  if (any(vapply(out, is.null, TRUE))) stop("index out of bounds in a Taylor vector", call. = FALSE)
  structure(out, class = "mcr_taylor")
}

#' @export
`[[.mcr_taylor` <- function(x, i, ...) .taylor_vec(list(unclass(x)[[i]]), if (is.character(i)) i else NULL)

#' @export
`$.mcr_taylor` <- function(x, name) x[[name]]

#' @export
`[<-.mcr_taylor` <- function(x, i, value) {
  y <- unclass(x)
  y[i] <- .taylor_elems(value)
  .taylor_vec(y, names(y))
}

#' @export
`[[<-.mcr_taylor` <- function(x, i, value) {
  y <- unclass(x)
  v <- .taylor_elems(value)
  if (length(v) != 1L) stop("replacement has length ", length(v), call. = FALSE)
  y[[i]] <- v[[1L]]
  .taylor_vec(y, names(y))
}

#' @export
c.mcr_taylor <- function(...) {
  parts <- list(...)
  out <- list()
  nms <- character()
  argn <- names(parts) %||% rep("", length(parts))
  for (k in seq_along(parts)) {
    el <- .taylor_elems(parts[[k]])
    pn <- names(parts[[k]]) %||% rep("", length(el))
    if (length(pn) != length(el)) pn <- rep("", length(el))
    if (nzchar(argn[k])) pn <- if (length(el) == 1L) argn[k] else paste0(argn[k], seq_along(el))
    out <- c(out, el)
    nms <- c(nms, pn)
  }
  .taylor_vec(out, if (any(nzchar(nms))) nms else NULL)
}

#' @export
as.list.mcr_taylor <- function(x, ...) {
  out <- lapply(seq_along(x), function(i) .taylor_vec(list(unclass(x)[[i]])))
  names(out) <- names(x)
  out
}

#' @export
Ops.mcr_taylor <- function(e1, e2) {
  a <- .taylor_elems(e1)
  if (missing(e2)) {
    if (.Generic == "-") return(.taylor_vec(lapply(a, function(x) -x), names(e1)))
    if (.Generic == "+") return(e1)
    stop("unary ", .Generic, " is not defined in Taylor mode", call. = FALSE)
  }
  b <- .taylor_elems(e2)
  na <- length(a); nb <- length(b)
  if (!na || !nb) return(.taylor_vec(list()))
  n <- max(na, nb)
  cmp <- .Generic %in% c("==", "!=", "<", ">", "<=", ">=")
  if (!cmp && !.Generic %in% c("+", "-", "*", "/", "^"))
    stop("operator ", .Generic, " is not defined in Taylor mode", call. = FALSE)
  res <- lapply(seq_len(n), function(i) {
    x <- a[[(i - 1L) %% na + 1L]]; y <- b[[(i - 1L) %% nb + 1L]]
    if (cmp) return(get(.Generic, envir = baseenv())(x[1L], y[1L]))
    .taylor_op(.Generic, x, y)
  })
  if (cmp) return(unlist(res))
  nms <- if (na >= nb && inherits(e1, "mcr_taylor")) names(e1) else names(e2)
  .taylor_vec(res, if (length(nms) == n) nms else NULL)
}

#' @export
chooseOpsMethod.mcr_taylor <- function(x, y, mx, my, cl, reverse) TRUE

#' @export
Math.mcr_taylor <- function(x, ...) {
  .taylor_vec(lapply(unclass(x), function(z) .taylor_math(.Generic, z)), names(x))
}

#' @export
Summary.mcr_taylor <- function(..., na.rm = FALSE) {
  el <- .taylor_elems(list(...))
  r <- switch(.Generic,
    sum = if (length(el)) Reduce(function(x, y) .taylor_op("+", x, y), el) else 0,
    prod = if (length(el)) Reduce(function(x, y) .taylor_op("*", x, y), el) else 1,
    stop(.Generic, "() is not defined in Taylor mode", call. = FALSE))
  .taylor_vec(list(r))
}

# A number as a Taylor polynomial of order K (promote).
.taylor_const <- function(x, K) c(x, numeric(K))

# The binary operations of TaylorDiff on two elements (numbers have length
# one). Each loop accumulates in the order of the unrolled @immutable code.
.taylor_op <- function(op, x, y) {
  nx <- length(x); ny <- length(y)
  if (nx == 1L && ny == 1L) return(get(op, envir = baseenv())(x, y))
  switch(op,
    "+" = if (nx == 1L) c(x + y[1L], y[-1L]) else if (ny == 1L) c(x[1L] + y, x[-1L]) else x + y,
    "-" = if (nx == 1L) c(x - y[1L], -y[-1L]) else if (ny == 1L) c(x[1L] - y, x[-1L]) else x - y,
    "*" = if (nx == 1L || ny == 1L) x * y else .taylor_mul(x, y),
    "/" = if (ny == 1L) x / y else .taylor_div(if (nx == 1L) .taylor_const(x, ny - 1L) else x, y),
    "^" = if (ny == 1L) .taylor_pow(x, y)
          else if (nx == 1L) .taylor_math("exp", y * log(x))
          else .taylor_math("exp", .taylor_mul(y, .taylor_math("log", x))))
}

.taylor_mul <- function(a, b) {
  P <- length(a) - 1L
  v <- numeric(P + 1L)
  for (i in 0:P) {
    s <- a[1L] * b[i + 1L]
    if (i > 0L) for (j in 1:i) s <- s + a[j + 1L] * b[i - j + 1L]
    v[i + 1L] <- s
  }
  v
}

.taylor_div <- function(a, b) {
  P <- length(a) - 1L
  v <- numeric(P + 1L)
  v[1L] <- a[1L] / b[1L]
  for (i in seq_len(P)) {
    s <- a[i + 1L]
    for (j in 0:(i - 1L)) s <- s - b[i - j + 1L] * v[j + 1L]
    v[i + 1L] <- s / b[1L]
  }
  v
}

# t^n: literal_pow for the integer powers 0, 1, 2, 3, -1, -2 (Julia's
# x^2 in source code), else the power recurrence.
.taylor_pow <- function(f, n) {
  if (n == 2) return(.taylor_mul(f, f))
  if (n == 3) return(.taylor_mul(.taylor_mul(f, f), f))
  if (n == 1) return(f)
  if (n == 0) return(.taylor_const(1, length(f) - 1L))
  if (n == -1) return(.taylor_div(.taylor_const(1, length(f) - 1L), f))
  if (n == -2) {
    i <- .taylor_div(.taylor_const(1, length(f) - 1L), f)
    return(.taylor_mul(i, i))
  }
  P <- length(f) - 1L
  v <- numeric(P + 1L)
  v[1L] <- f[1L]^n
  for (i in seq_len(P)) {
    s <- 0
    for (j in 0:(i - 1L)) {
      term <- ((n * (i - j) - j) * v[j + 1L]) * f[i - j + 1L]
      s <- if (j == 0L) term else s + term
    }
    v[i + 1L] <- s / (i * f[1L])
  }
  v
}

# exp, sin and cos (hand written recurrences of TaylorDiff), sqrt (t^0.5),
# log (first order rule, then raiseinv), abs (first order rule, sign).
.taylor_math <- function(fn, f) {
  if (length(f) == 1L) return(get(fn, envir = baseenv())(f))
  P <- length(f) - 1L
  switch(fn,
    exp = {
      v <- numeric(P + 1L)
      v[1L] <- exp(f[1L])
      for (i in seq_len(P)) {
        s <- 0
        for (j in 0:(i - 1L)) {
          term <- ((i - j) * v[j + 1L]) * f[i - j + 1L]
          s <- if (j == 0L) term else s + term
        }
        v[i + 1L] <- s / i
      }
      v
    },
    sin = , cos = {
      s <- numeric(P + 1L); co <- numeric(P + 1L)
      s[1L] <- sin(f[1L]); co[1L] <- cos(f[1L])
      for (i in seq_len(P)) {
        for (j in 0:(i - 1L)) {
          ts <- ((i - j) * co[j + 1L]) * f[i - j + 1L]
          tc <- ((i - j) * s[j + 1L]) * f[i - j + 1L]
          if (j == 0L) { s[i + 1L] <- ts; co[i + 1L] <- -tc } else {
            s[i + 1L] <- s[i + 1L] + ts; co[i + 1L] <- co[i + 1L] - tc
          }
        }
        s[i + 1L] <- s[i + 1L] / i; co[i + 1L] <- co[i + 1L] / i
      }
      if (fn == "sin") s else co
    },
    sqrt = .taylor_pow(f, 0.5),
    log = {
      if (P == 1L) return(c(log(f[1L]), 0 + f[2L] / f[1L]))
      # raiseinv(log(f0), z, t) = integrate(differentiate(t) / z, f0), z the
      # polynomial truncated to order P - 1
      d <- .taylor_div(f[-1L] * seq_len(P), f[seq_len(P)])
      c(log(f[1L]), d / seq_len(P))
    },
    abs = if (P == 1L) c(abs(f[1L]), 0 + sign(f[1L]) * f[2L]) else
      stop("abs() of order above one is not supported in Taylor mode", call. = FALSE),
    stop(fn, "() is not supported in Taylor mode; use polynomial, rational, exp, log, ",
         "sqrt, sin or cos terms", call. = FALSE))
}

# Coefficient `i` (0 for the value) of each component of f's result.
.taylor_coef <- function(fu, i, n) {
  el <- .taylor_elems(fu)
  if (length(el) != n) stop("f returned ", length(el), " components for a state of length ", n, call. = FALSE)
  vapply(el, function(x) if (length(x) == 1L) (if (i == 0L) x else 0) else x[i + 1L], 0)
}

# The Taylor jet of order P at (t, u) as an n x (P + 1) coefficient matrix
# (build_jet): u0 + f0 t, then each order from coefficient index - 1 of f on
# the polynomial so far, divided by index.
.taylor_jet <- function(ctx, t, u, P) {
  f <- ctx$f; p <- ctx$p; n <- length(u)
  C <- cbind(as.numeric(u), as.numeric(f(u, p, t)))
  if (P >= 2L) for (index in 2:P) {
    tT <- .taylor_vec(list(c(t, 1, numeric(index - 2L))))
    uT <- .taylor_vec(lapply(seq_len(n), function(i) C[i, ]), names(u))
    C <- cbind(C, .taylor_coef(f(uT, p, tT), index - 1L, n) / index)
  }
  C
}

# evaluate_polynomial: Horner's rule (Qin Jiushao), not fused.
.taylor_horner <- function(C, z) {
  ex <- C[, ncol(C)]
  for (i in rev(seq_len(ncol(C) - 1L))) ex <- C[, i] + z * ex
  ex
}

# ---------------------------------------------------------------------------

#' Explicit Taylor series methods
#'
#' `ExplicitTaylor2()` is the fixed step second order Taylor method: the
#' second derivative comes from evaluating `f` on first order Taylor numbers.
#' `ExplicitTaylor(order)` is the adaptive Taylor method of the given order
#' with the error estimate of its last term and the PI controller of that
#' order. `ExplicitTaylorAdaptiveOrder(min_order, max_order)` computes the
#' jet one order above its current order every step and moves to the order,
#' among the current one and its neighbours, that minimizes the work
#' (order squared over the proposed step); it starts at `max_order - 1`.
#'
#' The right-hand side `f(u, p, t)` is evaluated in Taylor mode arithmetic
#' (TaylorDiff's recurrences): it may use `+`, `-`, `*`, `/`, `^`, `exp`,
#' `log`, `sqrt`, `sin`, `cos`, `sum`, indexing and `c()`, and must return its
#' result built from those (as `c(...)`). Julia builds the jet once,
#' symbolically, so the higher order coefficients agree with it up to
#' rounding; `ExplicitTaylor2()` matches it exactly. The dense output is the
#' Taylor polynomial of each step. Function evaluations are counted as in
#' Julia: 3 per step for `ExplicitTaylor2()`, `order + 1` for
#' `ExplicitTaylor()`, the current order plus one for the adaptive order
#' method.
#'
#' @param order Order of the Taylor expansion (a positive integer).
#' @param min_order,max_order Bounds of the order of
#'   `ExplicitTaylorAdaptiveOrder()`; `max_order` must exceed `min_order`.
#' @param stage_limiter Stored only (not called, as in Julia).
#' @param step_limiter Optional `function(u, integrator, p, t)` applied after
#'   every accepted step.
#' @param thread Ignored.
#' @return An ODE algorithm.
#' @examples
#' prob <- ODEProblem(function(u, p, t) c(u[2], -u[1]), c(0, 1), c(0, 1))
#' solve(prob, ExplicitTaylor2(), dt = 0.1)$u[[11]]
#' solve(prob, ExplicitTaylor(order = 6), reltol = 1e-8, abstol = 1e-10)$u[[2]]
#' sol <- solve(prob, ExplicitTaylorAdaptiveOrder(), reltol = 1e-8, abstol = 1e-10)
#' sol$u[[length(sol$u)]]
#' @name taylor_methods
NULL

#' @rdname taylor_methods
#' @export
ExplicitTaylor2 <- function(stage_limiter = NULL, step_limiter = NULL, thread = NULL) {
  .new_alg("ExplicitTaylor2", 2L, adaptive = FALSE, step = .taylor2_step, interp = .taylor2_interp,
           init = function(ctx, t, u) list(nf0 = 0L), on_event = .taylor2_on_event,
           step_limiter = step_limiter, stage_limiter = stage_limiter)
}

# ExplicitTaylor2 keeps the FSAL default of OrdinaryDiffEq, so an event
# evaluates the (unused) FSAL derivative once more.
.taylor2_on_event <- function(ctx, t, u, cache) list(cache = cache, nf = 1L)

# perform_step! of ExplicitTaylor2ConstantCache:
# u = uprev + dt k1 + dt^2 / 2 k2 under @muladd (fused, element by element).
.taylor2_step <- function(ctx, t, u, dt, cache) {
  f <- ctx$f; p <- ctx$p
  k1 <- as.numeric(f(u, p, t))
  u1 <- .taylor_vec(lapply(seq_along(u), function(i) c(u[i], k1[i])), names(u))
  k2 <- .taylor_coef(f(u1, p, .taylor_vec(list(c(t, 1)))), 1L, length(u))
  unew <- .fma(dt * dt / 2, k2, .fma(dt, k1, u))
  list(u = unew, k = list(k1, k2), cache = cache, nf = 3L)
}

.taylor2_interp <- function(theta, dt, y0, y1, k) {
  h <- theta * dt
  .fma(k[[2L]] / 2, h * h, .fma(k[[1L]], h, y0))
}

# _ode_interpolant of the ExplicitTaylor caches: Horner's rule over the
# saved coefficients, y0 + h (k1 + h (k2 + ... + h kP)), not fused.
.taylor_interp <- function(theta, dt, y0, y1, k) {
  h <- theta * dt
  P <- length(k)
  res <- k[[P]]
  if (P > 1L) for (i in (P - 1L):1L) res <- k[[i]] + h * res
  y0 + h * res
}

#' @rdname taylor_methods
#' @export
ExplicitTaylor <- function(order = 1L, stage_limiter = NULL, step_limiter = NULL, thread = NULL) {
  order <- as.integer(order)
  if (length(order) != 1L || is.na(order) || order < 1L) stop("order must be a positive integer", call. = FALSE)
  .new_alg("ExplicitTaylor", order, step = function(ctx, t, u, dt, cache) .taylor_step(ctx, t, u, dt, cache, order),
           interp = .taylor_interp, init = function(ctx, t, u) list(nf0 = 0L), on_event = .taylor_on_event,
           step_limiter = step_limiter, stage_limiter = stage_limiter)
}

# perform_step! of ExplicitTaylorConstantCache{P}.
.taylor_step <- function(ctx, t, u, dt, cache, P) {
  C <- .taylor_jet(ctx, t, u, P)
  unew <- .taylor_horner(C, dt)
  EEst <- if (isTRUE(ctx$adaptive)) .firk_err(C[, P + 1L] * dt^(P + 1L), u, unew, ctx$abstol, ctx$reltol) else 0
  list(u = unew, EEst = EEst, k = lapply(seq_len(P), function(i) C[, i + 1L]), cache = cache, nf = P + 1L)
}

#' @rdname taylor_methods
#' @export
ExplicitTaylorAdaptiveOrder <- function(min_order = 1L, max_order = 10L, stage_limiter = NULL,
                                        step_limiter = NULL, thread = NULL) {
  min_order <- as.integer(min_order); max_order <- as.integer(max_order)
  if (min_order < 1L) stop("min_order must be a positive integer", call. = FALSE)
  if (max_order <= min_order)
    stop("ExplicitTaylorAdaptiveOrder requires max_order > min_order, got min_order = ", min_order,
         ", max_order = ", max_order, call. = FALSE)
  # alg_order is min_order (the PI controller's beta1 and beta2); the initial
  # step uses the current order, max_order - 1
  .new_alg("ExplicitTaylorAdaptiveOrder", max_order - 1L,
           step = function(ctx, t, u, dt, cache) .taylor_ao_step(ctx, t, u, dt, cache, min_order, max_order),
           interp = .taylor_interp, init = .taylor_ao_init, on_event = .taylor_on_event,
           beta1 = 7 / (10 * min_order), beta2 = 2 / (5 * min_order),
           step_limiter = step_limiter, stage_limiter = stage_limiter)
}

# The caches are not rebuilt after a callback changes u (and these methods
# are not FSAL, so nothing is evaluated again).
.taylor_on_event <- function(ctx, t, u, cache) list(cache = cache, nf = 0L)

# current_order = max_order - 1 when the cache is built; it is kept when a
# callback changes u (the cache is not rebuilt).
.taylor_ao_init <- function(ctx, t, u) {
  cache <- new.env(parent = emptyenv())
  cache$current_order <- NA_integer_
  cache$nf0 <- 0L
  cache
}

# perform_step! of ExplicitTaylorAdaptiveOrderConstantCache: the jet one
# order above the current one, then for each candidate order i the error
# estimate of its last term and a pseudo step of the PI controller (whose
# state is restored), keeping the order of least work i^2 / dtpropose.
.taylor_ao_step <- function(ctx, t, u, dt, cache, min_order, max_order) {
  if (is.na(cache$current_order)) cache$current_order <- max_order - 1L
  cur <- cache$current_order
  C <- .taylor_jet(ctx, t, u, cur + 1L)
  unew <- .taylor_horner(C, dt)
  k <- lapply(seq_len(max_order), function(i) if (i <= cur + 1L) C[, i + 1L] else numeric(length(u)))
  EEst <- 0
  if (isTRUE(ctx$adaptive)) {
    min_work <- Inf
    EEst <- ctx$EEst_prev
    co <- ctx$ctrl$co; cs <- ctx$ctrl$cs
    for (i in max(min_order, cur - 1L):min(max_order - 1L, cur + 1L)) {
      A <- i * i
      EEi <- .firk_err(C[, i + 1L] * dt^i, u, unew, ctx$abstol, ctx$reltol)
      qq <- .ctrl_q(co, cs, EEi, ctx$success_iter == 0L, abs(dt), NULL, NULL)
      dtpropose <- .ctrl_accept(co, qq$cs, qq$q, EEi, dt, NULL, NULL)$dt
      work <- A / dtpropose
      if (work < min_work) {
        cache$current_order <- i
        min_work <- work
        EEst <- EEi
      }
    }
  }
  list(u = unew, EEst = EEst, k = k, cache = cache, nf = cur + 1L)
}
