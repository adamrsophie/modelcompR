# More bracketing solvers, after BracketingNonlinearSolve.jl: alefeld.jl
# (Alefeld, Potra and Shi's algorithm 4.2 with common.jl's bracket,
# newton_quadratic and ipzero), modAB.jl (bisection switching to the
# Anderson-Bjork method), muller.jl and ridder.jl.

#' More bracketing algorithms
#'
#' Algorithms for [solve()] on an [IntervalNonlinearProblem()], beside
#' [ITP()], [Bisection()], [Brent()] and [Falsi()]: `Alefeld()` is the
#' enclosing method of Alefeld, Potra and Shi with inverse cubic
#' interpolation; `ModAB()` starts with bisection and switches to the
#' Anderson-Bjork regula falsi once the function looks linear on the
#' bracket; `Ridder()` is Ridder's exponential fit; `Muller(middle)` is
#' Muller's quadratic interpolation from the two ends and a middle point
#' (the midpoint by default), which needs no sign change.
#'
#' @param middle The third starting point of `Muller()`.
#' @return An algorithm object.
#' @examples
#' prob <- IntervalNonlinearProblem(function(u, p) u^3 - 2, c(0, 2))
#' solve(prob, Alefeld())$u
#' solve(prob, ModAB())$u
#' solve(prob, Ridder())$u
#' solve(prob, Muller())$u
#' @name interval_algorithms_more
NULL

.br_alg <- function(name, run, ...) {
  a <- .nl_alg(name, cls = "IntervalAlgorithm", ...)
  a$run <- run
  a
}

.br_bracket <- function(f, a, b, c) {
  fc <- f(c)
  if (fc == 0) return(c(a, b, c))
  if (f(a) * fc < 0) return(c(a, c, b))
  if (f(b) * fc < 0) return(c(c, b, a))
  c(a, b, c)
}

.br_newton_quadratic <- function(f, a, b, d, k) {
  A <- ((f(d) - f(b)) / (d - b) - (f(b) - f(a)) / (b - a)) / (d - a)
  B <- (f(b) - f(a)) / (b - a)
  if (A == 0) return(a - (1 / B) * f(a))
  r <- if (A * f(a) > 0) a else b
  for (i in seq_len(k)) r <- r - (f(a) + B * (r - a) + A * (r - a) * (r - b)) / (B + A * (2 * r - a - b))
  r
}

.br_ipzero <- function(f, a, b, c, d) {
  Q11 <- (c - d) * f(c) / (f(d) - f(c))
  Q21 <- (b - c) * f(b) / (f(c) - f(b))
  Q31 <- (a - b) * f(a) / (f(b) - f(a))
  D21 <- (b - c) * f(c) / (f(c) - f(b))
  D31 <- (a - b) * f(b) / (f(b) - f(a))
  Q22 <- (D21 - Q11) * f(b) / (f(d) - f(b))
  Q32 <- (D31 - Q21) * f(a) / (f(c) - f(a))
  D32 <- (D31 - Q21) * f(c) / (f(c) - f(a))
  Q33 <- (D32 - Q22) * f(a) / (f(d) - f(a))
  a + Q31 + Q32 + Q33
}

.br_res <- function(u, fu, l, r, code, it) list(u = u, resid = fu, left = l, right = r, retcode = code, iters = it)

.br_alefeld <- function(f, a, b, fa, fb, abstol, maxiters) {
  c <- a - (b - a) / (fb - fa) * fa
  fc <- f(c)
  if (a == c || b == c) return(.br_res(c, fc, a, b, "FloatingPointLimit", 1L))
  if (fc == 0) return(.br_res(c, fc, c, c, "Success", 1L))
  abd <- .br_bracket(f, a, b, c); a <- abd[1L]; b <- abd[2L]; d <- abd[3L]
  if (abs((b - a) / 2) < abstol) return(.br_res(c, fc, a, b, "Success", 1L))
  e <- d
  distinct <- function(v) length(unique(v)) == length(v)
  for (i in 2:max(2L, maxiters)) {
    if (i == 2L) c <- .br_newton_quadratic(f, a, b, d, 2L)
    else {
      if (!distinct(c(f(a), f(b), f(d), f(e)))) c <- .br_newton_quadratic(f, a, b, d, 2L)
      else {
        c <- .br_ipzero(f, a, b, d, e)
        if ((c - a) * (c - b) >= 0) c <- .br_newton_quadratic(f, a, b, d, 2L)
      }
    }
    ebar <- d; fc <- f(c)
    if (a == c || b == c) return(.br_res(c, fc, a, b, "FloatingPointLimit", i))
    if (fc == 0) return(.br_res(c, fc, c, c, "Success", i))
    z <- .br_bracket(f, a, b, c); ab <- z[1L]; bb <- z[2L]; db <- z[3L]
    if (abs((bb - ab) / 2) < abstol) return(.br_res(c, fc, ab, bb, "Success", i))
    if (!distinct(c(f(ab), f(bb), f(db), f(ebar)))) c <- .br_newton_quadratic(f, ab, bb, db, 3L)
    else {
      c <- .br_ipzero(f, ab, bb, db, ebar)
      if ((c - ab) * (c - bb) >= 0) c <- .br_newton_quadratic(f, ab, bb, db, 3L)
    }
    fc <- f(c)
    if (ab == c || bb == c) return(.br_res(c, fc, ab, bb, "FloatingPointLimit", i))
    if (fc == 0) return(.br_res(c, fc, c, c, "Success", i))
    z <- .br_bracket(f, ab, bb, c); ab <- z[1L]; bb <- z[2L]; db <- z[3L]
    if (abs((bb - ab) / 2) < abstol) return(.br_res(c, fc, ab, bb, "Success", i))
    u <- if (abs(f(ab)) < abs(f(bb))) ab else bb
    c <- u - 2 * (bb - ab) / (f(bb) - f(ab)) * f(u)
    if (abs(c - u) > 0.5 * (bb - ab)) c <- 0.5 * (ab + bb)
    fc <- f(c)
    if (ab == c || bb == c) return(.br_res(c, fc, ab, bb, "FloatingPointLimit", i))
    if (fc == 0) return(.br_res(c, fc, c, c, "Success", i))
    z <- .br_bracket(f, ab, bb, c); ab <- z[1L]; bb <- z[2L]; d <- z[3L]
    if (abs((bb - ab) / 2) < abstol) return(.br_res(c, fc, ab, bb, "Success", i))
    if (bb - ab < 0.5 * (b - a)) {
      a <- ab; b <- bb; e <- db
    } else {
      e <- d
      c <- 0.5 * (ab + bb); fc <- f(c)
      if (ab == c || bb == c) return(.br_res(c, fc, ab, bb, "FloatingPointLimit", i))
      if (fc == 0) return(.br_res(c, fc, c, c, "Success", i))
      z <- .br_bracket(f, ab, bb, c); a <- z[1L]; b <- z[2L]; d <- z[3L]
    }
    if (abs((b - a) / 2) < abstol) return(.br_res(c, fc, a, b, "Success", i))
  }
  if (b == c) b <- d else if (a == c) a <- d
  .br_res(c, f(c), a, b, "MaxIters", maxiters)
}

#' @rdname interval_algorithms_more
#' @export
Alefeld <- function() .br_alg("Alefeld", .br_alefeld)

.br_modab <- function(f, x1, x2, y1, y2, abstol, maxiters) {
  same <- function(x, y) (x < 0 && y < 0) || (x > 0 && y > 0)
  secant <- function(x1, y1, x2, y2) {
    a <- abs(y1); b <- abs(y2); den <- a + b
    if (is.infinite(den)) {
      if (is.infinite(a) || is.infinite(b)) return(x1 / 2 + x2 / 2)
      a <- a / 2; b <- b / 2; den <- a + b
    }
    min(max((b / den) * x1 + (a / den) * x2, x1), x2)
  }
  abf <- function(y3, y) { m <- 1 - y3 / y; if (m > 0) m else 1 / 2 }
  bisecting <- TRUE; side <- 0L
  threshold <- x2 - x1; C <- 2
  f1 <- y1; f2 <- y2; yMin <- 0
  i <- 1L
  while (i < maxiters) {
    if (bisecting) {
      x3 <- x1 / 2 + x2 / 2; y3 <- f(x3)
      if (is.finite(f2 - f1)) {
        ym <- (f1 + f2) / 2
        r <- 1 - abs(ym / (f2 - f1)); k <- r * r
        if (abs(ym - y3) < k * abs(y3) + k * abs(ym)) {
          bisecting <- FALSE
          threshold <- C * (x2 - x1)
          y1 <- f1; y2 <- f2
        }
      }
    } else {
      x3 <- secant(x1, y1, x2, y2)
      y3 <- if (x3 == x1) f1 else if (x3 == x2) f2 else f(x3)
      threshold <- threshold / 2
      yMin <- min(abs(f1), abs(f2))
    }
    if (y3 == 0) return(.br_res(x3, y3, x3, x3, "Success", i))
    if (is.nan(y3)) return(.br_res(x3, y3, x1, x2, "Failure", i))
    if ((x2 - x1) < 2 * abstol) return(.br_res(x3, y3, x1, x2, "Success", i))
    if (bisecting) {
      if (same(f1, y3)) { x1 <- x3; f1 <- y3 } else { x2 <- x3; f2 <- y3 }
    } else {
      if (same(f1, y3)) {
        if (side == 1L) y2 <- y2 * abf(y3, y1)
        x1 <- x3; y1 <- y3; f1 <- y3; side <- 1L
      } else {
        if (side == -1L) y1 <- y1 * abf(y3, y2)
        x2 <- x3; y2 <- y3; f2 <- y3; side <- -1L
      }
      if (x2 - x1 > threshold && abs(y3) > yMin / 2) { bisecting <- TRUE; side <- 0L }
    }
    nx <- x1 + max(abs(x1) * .Machine$double.eps, .Machine$double.xmin)
    if (nx >= x2) return(.br_res(x2, f(x2), x1, x2, "FloatingPointLimit", i))
    i <- i + 1L
  }
  .br_res(x1, f1, x1, x2, "MaxIters", maxiters)
}

#' @rdname interval_algorithms_more
#' @export
ModAB <- function() .br_alg("ModAB", .br_modab)

.br_ridder <- function(f, left, right, fl, fr, abstol, maxiters) {
  for (i in seq_len(maxiters)) {
    mid <- (left + right) / 2
    if (mid == left || mid == right) return(.br_res(left, fl, left, right, "FloatingPointLimit", i))
    fm <- f(mid)
    if (fm == 0) return(.br_res(mid, fm, mid, mid, "Success", i))
    s <- sqrt(fm^2 - fl * fr)
    if (s == 0) return(.br_res(left, fl, left, right, "Failure", i))
    x <- mid + (mid - left) * sign(fl - fm) * fm / s
    fx <- f(x)
    if (fx == 0) return(.br_res(x, fx, x, x, "Success", i))
    if (abs((right - left) / 2) < abstol) return(.br_res(mid, fm, left, right, "Success", i))
    if (sign(fx) != sign(fm)) { left <- mid; fl <- fm; right <- x; fr <- fx }
    else if (sign(fx) != sign(fl)) { right <- x; fr <- fx }
    else { left <- x; fl <- fx }
    if (left > right) { tmp <- left; left <- right; right <- tmp; tf <- fl; fl <- fr; fr <- tf }
  }
  .br_res(left, fl, left, right, "MaxIters", maxiters)
}

#' @rdname interval_algorithms_more
#' @export
Ridder <- function() .br_alg("Ridder", .br_ridder)

#' @rdname interval_algorithms_more
#' @export
Muller <- function(middle = NULL) {
  run <- function(f, a, b, fa, fb, abstol, maxiters) {
    x2 <- a; x0 <- b
    x1 <- middle %||% ((x2 + x0) / 2)
    f2 <- fa; f1 <- f(x1); f0 <- fb
    xn <- x2; fn <- f2
    for (i in seq_len(maxiters)) {
      q <- (x0 - x1) / (x1 - x2)
      A <- q * f0 - q * (1 + q) * f1 + q^2 * f2
      B <- (2 * q + 1) * f0 - (1 + q)^2 * f1 + q^2 * f2
      C <- (1 + q) * f0
      disc <- sqrt(as.complex(B^2 - 4 * A * C))
      dp <- B + disc; dm <- B - disc
      step <- if (Mod(dp) >= Mod(dm)) (x0 - x1) * 2 * C / dp else (x0 - x1) * 2 * C / dm
      xn <- Re(x0 - step)
      fn <- f(xn)
      if (abstol >= abs(fn)) return(.br_res(xn, fn, xn, xn, "Success", i))
      x2 <- x1; x1 <- x0; x0 <- xn
      f2 <- f1; f1 <- f0; f0 <- fn
    }
    .br_res(xn, fn, xn, xn, "MaxIters", maxiters)
  }
  a <- .br_alg("Muller", run, middle = middle)
  a$no_bracket <- TRUE
  a
}
