# Forward mode automatic differentiation with dual numbers, after
# ForwardDiff.jl 1.x (src/dual.jl, src/partials.jl, src/prelude.jl) with the
# derivative rules of DiffRules.jl (src/rules.jl), and the Hessian that
# DifferentiationInterface 0.7 computes for SecondOrder(AutoForwardDiff(),
# AutoForwardDiff()) (src/second_order/hessian.jl and hvp.jl: a pushforward
# of the ForwardDiff gradient for every batch of basis seeds), which is how
# OptimizationBase differentiates an OptimizationFunction with
# AutoForwardDiff().
#
# A dual number (class "fdual") holds a vector of values `v` and a list `p`
# of partials, one vector per perturbation, each as long as `v`; `tag` is
# the nesting level (ForwardDiff's tag order). The values and partials of a
# nested dual are themselves duals of a lower tag, so a Hessian is a dual
# of duals, as in ForwardDiff. Every operation applies ForwardDiff's rule
# with the same arithmetic in the same order, so the derivatives agree with
# Julia bit for bit when the R function computes its values as the Julia
# function does.

.ad_tag <- function(x) if (inherits(x, "fdual")) .subset2(x, "tag") else 0L
.ad_v <- function(x) .subset2(x, "v")
.ad_p <- function(x) .subset2(x, "p")

# recycle a value (numeric or dual) to length n
.ad_recycle <- function(a, n) {
  la <- length(a)
  if (la == n) a else a[rep_len(seq_len(la), n)]
}

.ad_new <- function(v, p, tag) {
  n <- length(v)
  structure(list(v = v, p = lapply(p, .ad_recycle, n), tag = tag), class = "fdual")
}

# the innermost (plain numeric) values
.ad_primal <- function(x) {
  while (inherits(x, "fdual")) x <- .ad_v(x)
  x
}

# elementwise iszero of a value or a partial (a dual is zero when its value
# and all its partials are zero)
.ad_iszero <- function(x) {
  if (!inherits(x, "fdual")) return(x == 0)
  z <- .ad_iszero(.ad_v(x))
  for (q in .ad_p(x)) z <- z & .ad_iszero(q)
  z
}

# elementwise iszero of all partials of x at tag `tag` (ForwardDiff's
# iszero(partials(x)) and isconstant)
.ad_const_at <- function(x, tag, n) {
  if (.ad_tag(x) != tag) return(rep(TRUE, n))
  z <- rep(TRUE, n)
  for (q in .ad_p(x)) z <- z & rep_len(.ad_iszero(q), n)
  z
}

# x at the elements where `keep` is FALSE, zero where it is TRUE
.ad_zero_where <- function(x, keep) {
  if (!any(keep)) return(x)
  .ad_assign(x, which(keep), 0)
}

# x as a dual of tag `tag` with K partials (a constant gets zero partials)
.ad_promote <- function(x, tag, K) {
  if (.ad_tag(x) == tag) return(x)
  z <- rep(0, length(x))
  .ad_new(x, rep(list(z), K), tag)
}

# x[i] <- value for numeric or dual x and value of any tags
.ad_assign <- function(x, i, value) {
  tag <- max(.ad_tag(x), .ad_tag(value))
  if (tag == 0L) {
    x[i] <- value
    return(x)
  }
  K <- length(.ad_p(if (.ad_tag(x) == tag) x else value))
  x <- .ad_promote(x, tag, K)
  value <- .ad_promote(value, tag, K)
  p <- Map(function(a, b) .ad_assign(a, i, b), .ad_p(x), .ad_p(value))
  .ad_new(.ad_assign(.ad_v(x), i, .ad_v(value)), p, tag)
}

# c() of numeric and dual pieces of any tags
.ad_c <- function(args) {
  args <- args[vapply(args, length, 0L) > 0L]
  tags <- vapply(args, .ad_tag, 0L)
  if (!length(args)) return(numeric())
  tag <- max(tags)
  if (tag == 0L) return(do.call(c, args))
  K <- length(.ad_p(args[[which(tags == tag)[1L]]]))
  args <- lapply(args, .ad_promote, tag, K)
  v <- .ad_c(lapply(args, .ad_v))
  p <- lapply(seq_len(K), function(k) .ad_c(lapply(args, function(a) .ad_p(a)[[k]])))
  .ad_new(v, p, tag)
}

#' @export
length.fdual <- function(x) length(.ad_v(x))

#' @export
`[.fdual` <- function(x, i) {
  if (missing(i)) return(x)
  .ad_new(.ad_v(x)[i], lapply(.ad_p(x), function(q) q[i]), .ad_tag(x))
}

#' @export
`[[.fdual` <- function(x, i) {
  if (length(i) != 1L) stop("subscript out of bounds", call. = FALSE)
  x[i]
}

#' @export
`[<-.fdual` <- function(x, i, value) .ad_assign(x, if (missing(i)) seq_len(length(x)) else i, value)

#' @export
`[[<-.fdual` <- function(x, i, value) .ad_assign(x, i, value)

#' @export
c.fdual <- function(...) .ad_c(list(...))

#' @export
rep.fdual <- function(x, ...) x[rep(seq_len(length(x)), ...)]

#' @export
rev.fdual <- function(x) x[rev(seq_len(length(x)))]

#' @export
as.double.fdual <- function(x, ...) {
  stop("cannot convert a dual number of forward mode differentiation to a double", call. = FALSE)
}

#' @export
is.na.fdual <- function(x) is.na(.ad_primal(x))

#' @export
mean.fdual <- function(x, ...) sum(x) / length(x)

#' @export
print.fdual <- function(x, ...) {
  cat("Dual (tag ", .ad_tag(x), ", ", length(.ad_p(x)), " partials)\n", sep = "")
  print(.ad_primal(x))
  invisible(x)
}

# x^k for the literal powers 0 to 3 (Base.literal_pow of ForwardDiff and
# the x*x, x*x*x of Julia's Float64)
.ad_litpow <- function(v, k) {
  if (inherits(v, "fdual")) return(v^k)
  switch(k + 1L, rep(1, length(v)), v, v * v, v * v * v)
}

.ad_inv <- function(x) if (inherits(x, "fdual")) .ad_unary("inv", x) else 1 / x
.ad_abs2 <- function(x) if (inherits(x, "fdual")) .ad_unary("abs2", x) else x * x
.ad_sech <- function(x) if (inherits(x, "fdual")) .ad_unary("sech", x) else 1 / cosh(x)

#' @export
Ops.fdual <- function(e1, e2) {
  g <- .Generic
  if (missing(e2)) {
    if (g == "-") return(.ad_new(-.ad_v(e1), lapply(.ad_p(e1), function(q) -q), .ad_tag(e1)))
    if (g == "+") return(e1)
    stop("operation ", g, " is not defined for dual numbers", call. = FALSE)
  }
  if (g %in% c("==", "!=", "<", ">", "<=", ">="))
    return(get(g, envir = baseenv())(.ad_primal(e1), .ad_primal(e2)))
  tag <- max(.ad_tag(e1), .ad_tag(e2))
  xd <- .ad_tag(e1) == tag
  yd <- .ad_tag(e2) == tag
  vx <- if (xd) .ad_v(e1) else e1
  vy <- if (yd) .ad_v(e2) else e2
  px <- if (xd) .ad_p(e1)
  py <- if (yd) .ad_p(e2)
  switch(g,
    "+" = {
      if (xd && yd) .ad_new(vx + vy, Map(`+`, px, py), tag)
      else if (xd) .ad_new(vx + vy, px, tag)
      else .ad_new(vx + vy, py, tag)
    },
    "-" = {
      if (xd && yd) .ad_new(vx - vy, Map(`-`, px, py), tag)
      else if (xd) .ad_new(vx - vy, px, tag)
      else .ad_new(vx - vy, lapply(py, function(q) -q), tag)
    },
    "*" = {
      if (xd && yd) .ad_new(vx * vy, Map(function(a, b) a * vy + b * vx, px, py), tag)
      else if (xd) .ad_new(vx * vy, lapply(px, function(a) a * vy), tag)
      else .ad_new(vx * vy, lapply(py, function(b) b * vx), tag)
    },
    "/" = {
      if (xd && yd) {
        val <- vx / vy
        .ad_new(val, Map(function(a, b) (a - val * b) / vy, px, py), tag)
      } else if (xd) {
        .ad_new(vx / vy, lapply(px, function(a) a / vy), tag)
      } else {
        divv <- vx / vy
        s <- -(divv / vy)
        .ad_new(divv, lapply(py, function(b) s * b), tag)
      }
    },
    "^" = .ad_pow(vx, vy, px, py, xd, yd, tag),
    stop("operation ", g, " is not defined for dual numbers", call. = FALSE))
}

.ad_pow <- function(vx, vy, px, py, xd, yd, tag) {
  if (xd && !yd) {
    if (!inherits(vy, "fdual") && length(vy) == 1L && vy %in% 0:3) {
      # Base.literal_pow
      k <- as.integer(vy)
      if (k == 0L) return(.ad_new(.ad_litpow(vx, 0L), lapply(px, function(a) 0 * a), tag))
      expv <- .ad_litpow(vx, k)
      deriv <- k * .ad_litpow(vx, k - 1L)
      return(.ad_new(expv, lapply(px, function(a) deriv * a), tag))
    }
    expv <- vx^vy
    n <- length(expv)
    pw <- vx^(vy - 1)
    zero <- rep_len(.ad_primal(vy) == 0, n) & rep_len(!inherits(vy, "fdual"), n)
    zp <- rep(TRUE, n)
    for (a in px) zp <- zp & rep_len(.ad_iszero(a), n)
    p <- lapply(px, function(a) .ad_zero_where(.ad_recycle(a * vy * pw, n), zero | zp))
    return(.ad_new(expv, p, tag))
  }
  if (!xd && yd) {
    expv <- vx^vy
    n <- length(expv)
    deriv <- .ad_zero_where(.ad_recycle(expv * log(vx), n),
                            rep_len(.ad_primal(vx) == 0, n) & rep_len(.ad_primal(vy) > 0, n))
    return(.ad_new(expv, lapply(py, function(b) deriv * b), tag))
  }
  expv <- vx^vy
  n <- length(expv)
  powval <- vy * vx^(vy - 1)
  ycon <- rep(TRUE, n)
  for (b in py) ycon <- ycon & rep_len(.ad_iszero(b), n)
  logval <- .ad_recycle(expv * log(vx), n)
  zx <- rep_len(.ad_primal(vx) == 0, n) & rep_len(.ad_primal(vy) > 0, n)
  logval <- .ad_zero_where(logval, zx & !ycon)
  if (any(ycon)) logval <- .ad_assign(logval, which(ycon), 1)
  .ad_new(expv, Map(function(a, b) a * powval + b * logval, px, py), tag)
}

# the unary rules: value and derivative from the value v of the dual (a
# dual of a lower tag inside a nested dual)
.ad_unary <- function(g, x, ...) {
  v <- .ad_v(x)
  tag <- .ad_tag(x)
  logtwo <- 0.6931471805599453
  logten <- 2.302585092994046
  r <- switch(g,
    abs = list(abs(v), ifelse(.ad_primal(v) < 0 | (1 / .ad_primal(v)) < 0, -1, 1)),
    abs2 = list(v * v, v + v),
    inv = {
      val <- .ad_inv(v)
      list(val, -.ad_abs2(val))
    },
    sqrt = {
      val <- sqrt(v)
      list(val, .ad_inv(2 * val))
    },
    exp = {
      val <- exp(v)
      list(val, val)
    },
    expm1 = list(expm1(v), exp(v)),
    log = {
      base <- list(...)$base
      if (is.null(base)) list(log(v), .ad_inv(v)) else list(log(v) / log(base), .ad_inv(v) / log(base))
    },
    log2 = list(log2(v), .ad_inv(logtwo * v)),
    log10 = list(log10(v), .ad_inv(logten * v)),
    log1p = list(log1p(v), .ad_inv(v + 1)),
    sin = list(sin(v), cos(v)),
    cos = list(cos(v), -sin(v)),
    tan = {
      val <- tan(v)
      list(val, 1 + .ad_litpow(val, 2L))
    },
    asin = list(asin(v), .ad_inv(sqrt((1 - v) * (1 + v)))),
    acos = list(acos(v), -.ad_inv(sqrt((1 - v) * (1 + v)))),
    atan = list(atan(v), .ad_inv(1 + .ad_litpow(v, 2L))),
    sinh = list(sinh(v), cosh(v)),
    cosh = list(cosh(v), sinh(v)),
    tanh = list(tanh(v), .ad_litpow(.ad_sech(v), 2L)),
    sech = {
      val <- .ad_sech(v)
      list(val, -tanh(v) * val)
    },
    asinh = list(asinh(v), .ad_inv(sqrt(.ad_litpow(v, 2L) + 1))),
    acosh = list(acosh(v), .ad_inv(sqrt(v - 1) * sqrt(v + 1))),
    atanh = list(atanh(v), .ad_inv((1 - v) * (1 + v))),
    sinpi = list(sinpi(v), pi * cospi(v)),
    cospi = list(cospi(v), -(pi * sinpi(v))),
    gamma = {
      val <- gamma(v)
      list(val, digamma(v) * val)
    },
    lgamma = list(lgamma(v), digamma(v)),
    digamma = list(digamma(v), trigamma(v)),
    trigamma = list(trigamma(v), .ad_psigamma(v, 2L)),
    stop("the function ", g, " is not defined for dual numbers", call. = FALSE))
  deriv <- r[[2L]]
  .ad_new(r[[1L]], lapply(.ad_p(x), function(a) deriv * a), tag)
}

.ad_psigamma <- function(x, k) {
  if (!inherits(x, "fdual")) return(psigamma(x, k))
  deriv <- .ad_psigamma(.ad_v(x), k + 1L)
  .ad_new(.ad_psigamma(.ad_v(x), k), lapply(.ad_p(x), function(a) deriv * a), .ad_tag(x))
}

#' @export
Math.fdual <- function(x, ...) {
  g <- .Generic
  if (g %in% c("sign", "floor", "ceiling", "round", "signif", "trunc"))
    return(get(g, envir = baseenv())(.ad_primal(x), ...))
  if (g == "cumsum" || g == "cumprod") {
    n <- length(x)
    if (n == 0L) return(x)
    out <- x[1L]
    acc <- x[1L]
    for (i in seq_len(n)[-1L]) {
      acc <- if (g == "cumsum") acc + x[i] else acc * x[i]
      out <- .ad_c(list(out, acc))
    }
    return(out)
  }
  .ad_unary(g, x, ...)
}

#' @export
Summary.fdual <- function(..., na.rm = FALSE) {
  g <- .Generic
  x <- .ad_c(list(...))
  n <- length(x)
  switch(g,
    sum = {
      if (n == 0L) return(0)
      # Julia's sum over a vector of duals: a sequential fold
      acc <- x[1L]
      for (i in seq_len(n)[-1L]) acc <- acc + x[i]
      acc
    },
    prod = {
      if (n == 0L) return(1)
      acc <- x[1L]
      for (i in seq_len(n)[-1L]) acc <- acc * x[i]
      acc
    },
    max = x[which.max(.ad_primal(x))],
    min = x[which.min(.ad_primal(x))],
    range = .ad_c(list(x[which.min(.ad_primal(x))], x[which.max(.ad_primal(x))])),
    stop("the function ", g, " is not defined for dual numbers", call. = FALSE))
}

# ForwardDiff's pickchunksize with the default threshold 12
.ad_chunksize <- function(n, threshold = 12L) {
  if (n <= threshold) return(as.integer(n))
  nchunks <- ceiling(n / threshold)
  as.integer(ceiling(n / nchunks))
}

.ad_unit <- function(n, k) {
  e <- numeric(n)
  e[k] <- 1
  e
}

# ForwardDiff.jacobian (and gradient) of f at x in chunks: the value and the
# matrix of derivatives, one row per output and one column per input
.ad_jacobian <- function(f, x, chunk = .ad_chunksize(length(x))) {
  n <- length(x)
  J <- NULL
  val <- NULL
  for (start in seq.int(1L, n, by = max(chunk, 1L))) {
    idx <- start:min(n, start + chunk - 1L)
    xd <- .ad_new(x, lapply(idx, function(k) .ad_unit(n, k)), 1L)
    y <- f(xd)
    if (inherits(y, "fdual")) {
      if (.ad_tag(y) != 1L) stop("unexpected dual number tag", call. = FALSE)
      yv <- .ad_v(y)
      m <- length(yv)
      cols <- vapply(.ad_p(y), function(q) as.numeric(.ad_recycle(q, m)), numeric(m))
    } else if (is.numeric(y)) {
      yv <- as.numeric(y)
      m <- length(yv)
      cols <- matrix(0, m, length(idx))
    } else {
      stop("the function does not return numbers for dual inputs", call. = FALSE)
    }
    if (is.null(J)) {
      J <- matrix(0, m, n)
      val <- as.numeric(yv)
    }
    J[, idx] <- cols
  }
  if (is.null(J)) J <- matrix(0, length(f(x)), 0L)
  list(value = val, jacobian = J)
}

.ad_gradient <- function(f, x) as.numeric(.ad_jacobian(f, x)$jacobian)

# The Hessian of DifferentiationInterface for SecondOrder(AutoForwardDiff(),
# AutoForwardDiff()): for each batch of basis seeds (batch size
# pickchunksize(n), the last batch padded with seeds from the start), the
# pushforward of the ForwardDiff gradient, column j = d(grad f)/dx_j. The
# outer perturbation has tag 1 and the gradient's tag 2.
.ad_hessian <- function(f, x) {
  n <- length(x)
  B <- .ad_chunksize(n)
  A <- ceiling(n / B)
  H <- matrix(0, n, n)
  zeros <- rep(list(numeric(n)), B)
  for (a in seq_len(A)) {
    cols <- 1L + ((a - 1L) * B + seq_len(B) - 1L) %% n
    xo <- .ad_new(x, lapply(cols, function(j) .ad_unit(n, j)), 1L)
    for (start in seq.int(1L, n, by = B)) {
      idx <- start:min(n, start + B - 1L)
      xi <- .ad_new(xo, lapply(idx, function(k) .ad_new(.ad_unit(n, k), zeros, 1L)), 2L)
      y <- f(xi)
      if (!inherits(y, "fdual") && !is.numeric(y))
        stop("the function does not return numbers for dual inputs", call. = FALSE)
      if (!inherits(y, "fdual") || .ad_tag(y) != 2L) next
      for (k in seq_along(idx)) {
        gk <- .ad_p(y)[[k]]
        if (.ad_tag(gk) != 1L) next
        for (b in seq_len(B)) H[idx[k], cols[b]] <- .ad_primal(.ad_p(gk)[[b]])
      }
    }
  }
  H
}

# Whether f runs on dual numbers at x (forward mode needs R code that only
# uses the operations above; a function that writes duals into a numeric
# vector, for example, does not).
.ad_works <- function(f, x) {
  ok <- tryCatch({
    xd <- .ad_new(as.numeric(x), list(.ad_unit(length(x), 1L)), 1L)
    y <- suppressWarnings(f(xd))
    inherits(y, "fdual") || (is.numeric(y) && !is.list(y))
  }, error = function(e) FALSE)
  isTRUE(ok)
}

# The derivatives of OptimizationBase.instantiate_function for
# AutoForwardDiff() from the objective f(x) and the constraints cons(x)
# (vector of length m): gradient, constraint Jacobian, Hessian and the
# Hessian of the Lagrangian sigma f + dot(lambda, cons) (lag_h; with
# `weighted_zero`, the matrix form that sums lambda_i times the constraint
# Hessians when sigma is 0). Hessians come back with their lower triangle
# mirrored, the triangle the solvers read.
.ad_opt_derivatives <- function(f, cons, m) {
  lower <- function(H) {
    H[upper.tri(H)] <- t(H)[upper.tri(H)]
    H
  }
  lagrangian <- function(sigma, lambda) function(z) {
    s <- sigma * f(z)
    if (!m) return(s)
    cz <- cons(z)
    acc <- 0
    for (i in seq_len(m)) acc <- acc + lambda[i] * cz[i]
    s + acc
  }
  list(
    grad = function(x) .ad_gradient(f, x),
    jac = function(x) .ad_jacobian(cons, x)$jacobian,
    hess = function(x) lower(.ad_hessian(f, x)),
    cons_h = function(x) lapply(seq_len(m), function(i) lower(.ad_hessian(function(z) cons(z)[i], x))),
    lag_h = function(x, sigma, lambda, weighted_zero = FALSE) {
      if (weighted_zero && sigma == 0) {
        H <- matrix(0, length(x), length(x))
        for (i in seq_len(m)) if (lambda[i] != 0)
          H <- H + lambda[i] * .ad_hessian(function(z) cons(z)[i], x)
        return(lower(H))
      }
      lower(.ad_hessian(lagrangian(sigma, lambda), x))
    })
}

# ForwardDiff.jacobian of f at x with the number of chunks (the f calls
# OrdinaryDiffEq counts), or NULL when f cannot run on dual numbers
.ad_try_jacobian <- function(f, x) {
  x <- as.numeric(x)
  r <- tryCatch(suppressWarnings(.ad_jacobian(function(z) {
    y <- f(z)
    if (!inherits(y, "fdual") && !(is.numeric(y) && !is.object(y))) stop("not a dual result")
    y
  }, x)), error = function(e) NULL)
  if (is.null(r) || !is.matrix(r$jacobian) || ncol(r$jacobian) != length(x) || anyNA(r$jacobian)) return(NULL)
  r$chunks <- as.integer(ceiling(length(x) / max(1L, .ad_chunksize(length(x)))))
  r
}
