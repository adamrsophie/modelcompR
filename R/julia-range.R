# Julia's floating point ranges start:step:stop (the (:) method for IEEE
# floats and the TwicePrecision StepRangeLen of base/twiceprecision.jl), as
# OrdinaryDiffEqCore builds the save times of a numeric saveat. Inputs with
# exact rational forms give the multiples of the rational step in twice
# precision, so 0.37:0.37:3 holds 1.11 rather than 3 * 0.37.

# rat: a continued fraction approximation n / d of x, with terms up to the
# largest integer of single precision
.jlr_rat <- function(x) {
  y <- x
  a <- 1; d <- 1; b <- 0; c <- 0
  m <- 16777216
  while (abs(y) <= m) {
    f <- trunc(y)
    y <- y - f
    a0 <- a; b0 <- b
    a <- f * a + c; c <- a0
    b <- f * b + d; d <- b0
    if (max(abs(a), abs(b)) > m) return(c(c, d))
    if (a / b == x) break
    y <- 1 / y
  }
  c(a, b)
}

.jlr_gcd <- function(a, b) {
  a <- abs(a); b <- abs(b)
  while (b != 0) {
    r <- a %% b
    a <- b; b <- r
  }
  a
}

# div of Julia: integer division rounding toward zero
.jlr_div <- function(a, b) sign(a) * sign(b) * (abs(a) %/% abs(b))

.jlr_between <- function(a, x, b) (a <= x && x <= b) || (b <= x && x <= a)

# truncbits: x with the nb lowest bits of its representation cleared
.jlr_truncbits <- function(x, nb) {
  if (nb <= 0) return(x)
  r <- writeBin(x, raw(), size = 8L, endian = "little")
  bits <- rawToBits(r)
  bits[seq_len(nb)] <- as.raw(0L)
  readBin(packBits(bits, "raw"), "double", size = 8L, endian = "little")
}

.jlr_canon <- function(big, little) {
  h <- big + little
  c(h, (big - h) + little)
}

# mul12: the product and its rounding error
.jlr_mul12 <- function(x, y) {
  p <- x * y
  if (p == 0 || !is.finite(p)) return(c(p, p))
  tx <- x * 134217729; xh <- tx - (tx - x); xl <- x - xh
  ty <- y * 134217729; yh <- ty - (ty - y); yl <- y - yh
  .jlr_canon(p, ((xh * yh - p) + xh * yl + xl * yh) + xl * yl)
}

# TwicePrecision{Float64}((n, d)): the integer n split into halves, divided
# by d in twice precision
.jlr_tp_ratio <- function(n, d) {
  hi <- .jlr_truncbits(n, 27L)
  x <- .jlr_canon(hi, n - hi)
  h <- x[1L] / d
  if (h == 0 || !is.finite(h)) return(c(h, h))
  u <- .jlr_mul12(h, d)
  .jlr_canon(h, (((x[1L] - u[1L]) - u[2L]) + x[2L]) / d)
}

.jlr_top_set_bit <- function(x) if (x <= 0) 0 else floor(log2(x)) + 1

.jlr_elements <- function(ref, step, len, offset) {
  vapply(seq_len(len), function(i) {
    u <- i - offset
    sh <- u * step[1L]; sl <- u * step[2L]
    xy <- if (abs(sh) > abs(ref[1L])) .jlr_canon(sh, ref[1L]) else .jlr_canon(ref[1L], sh)
    xy[1L] + (xy[2L] + (sl + ref[2L]))
  }, 0)
}

.jl_range <- function(start, step, stop) {
  if (step == 0) stop("range step cannot be zero", call. = FALSE)
  sr <- .jlr_rat(step)
  if (sr[2L] != 0 && sr[1L] / sr[2L] == step) {
    ar <- .jlr_rat(start)
    zr <- .jlr_rat(stop)
    if (ar[2L] != 0 && zr[2L] != 0 && ar[1L] / ar[2L] == start && zr[1L] / zr[2L] == stop) {
      den <- ar[2L] * .jlr_div(sr[2L], .jlr_gcd(ar[2L], sr[2L]))
      m <- 2^53
      if (den != 0 && abs(start * den) <= m && abs(step * den) <= m &&
          den %% ar[2L] == 0 && den %% sr[2L] == 0) {
        start_n <- ar[1L] * .jlr_div(den, ar[2L])
        step_n <- sr[1L] * .jlr_div(den, sr[2L])
        len <- max(0, .jlr_div(den * zr[1L] - zr[2L] * start_n + step_n * zr[2L], step_n * zr[2L]))
        if (.jlr_between(start, start + (len - 1) * step, stop + step / 2) &&
            !.jlr_between(start, start + len * step, stop)) {
          # floatrange: the reference is the element of smallest magnitude
          if (len < 2) {
            ref <- .jlr_tp_ratio(start_n, den)
            return(.jlr_elements(ref, c(0, 0), len, 1))
          }
          imin <- min(max(round(-start_n / step_n + 1), 1), len)
          ref <- .jlr_tp_ratio(start_n + (imin - 1) * step_n, den)
          nb <- min(27, if (len < 2) 0 else .jlr_top_set_bit(max(imin - 1, len - imin) - 1) + 1)
          st <- .jlr_tp_ratio(step_n, den)
          hi <- .jlr_truncbits(st[1L], nb)
          return(.jlr_elements(ref, c(hi, (st[1L] - hi) + st[2L]), len, imin))
        }
      }
    }
  }
  # start and step taken literally
  lf <- (stop - start) / step
  len <- if (lf < 0) 0 else if (lf == 0) 1 else {
    n <- round(lf) + 1
    s2 <- start + (n - 1) * step
    n - ((start < stop && stop < s2) + (start > stop && stop > s2))
  }
  .jlr_elements(c(start, 0), c(step, 0), len, 1)
}

# range(start, stop; length = len) of Julia for Float64 endpoints
# (range_start_stop_length and _linspace of base/twiceprecision.jl): exact
# rational endpoints give the ratios in twice precision, other endpoints a
# twice precision reference at the element of smallest magnitude. Returns
# the elements, with the step(r) of the range as attribute "step".
.jl_linrange <- function(start, stop, len) {
  out <- function(ref, step, offset) {
    x <- .jlr_elements(ref, step, len, offset)
    attr(x, "step") <- step[1L] + step[2L]
    x
  }
  if (len < 2) {
    if (len == 1 && start != stop) stop("range: endpoints differ", call. = FALSE)
    return(out(c(start, 0), c(start, -stop), 1))
  }
  if (start == stop) return(out(c(start, 0), c(0, 0), 1))
  nb_len <- function(imin) min(27, .jlr_top_set_bit(max(imin - 1, len - imin) - 1) + 1)
  ar <- .jlr_rat(start); zr <- .jlr_rat(stop)
  if (ar[2L] != 0 && zr[2L] != 0) {
    den <- ar[2L] * .jlr_div(zr[2L], .jlr_gcd(ar[2L], zr[2L]))
    m <- 2^53
    if (den != 0 && abs(den * start) <= m && abs(den * stop) <= m) {
      sn <- round(den * start); en <- round(den * stop)
      if (sn / den == start && en / den == stop) {
        tmin <- -sn / (en - sn)
        imin <- min(max(round(tmin * (len - 1) + 1), 1), len)
        ref_num <- (len - imin) * sn + (imin - 1) * en
        ref_den <- (len - 1) * den
        ref <- .jlr_tp_ratio(ref_num, ref_den)
        st <- .jlr_tp_ratio(en - sn, ref_den)
        hi <- .jlr_truncbits(st[1L], nb_len(imin))
        return(out(ref, c(hi, (st[1L] - hi) + st[2L]), imin))
      }
    }
  }
  if (!is.finite(start) || !is.finite(stop)) stop("start and stop must be finite", call. = FALSE)
  delta <- stop - start; dfac <- 1
  if (!is.finite(delta)) { delta <- stop / len - start / len; dfac <- len }
  tmin <- -(start / delta) / dfac
  lenn1 <- len - 1
  imin <- round(tmin * lenn1 + 1)
  if (1 < imin && imin < len) {
    t <- (imin - 1) / lenn1
    ref <- (1 - t) * start + t * stop
    step <- if (imin - 1 < len - imin) (ref - start) / (imin - 1) else (stop - ref) / (len - imin)
  } else if (imin <= 1) {
    imin <- 1; ref <- start; step <- (delta / lenn1) * dfac
  } else {
    imin <- len; ref <- stop; step <- (delta / lenn1) * dfac
  }
  mx <- 1.7976931348623155e308
  k <- max(imin - 1, len - imin)
  step_pre <- min(max(step, max(-(mx + ref) / k, (-mx + ref) / k)), min((mx - ref) / k, (mx + ref) / k))
  step_hi <- .jlr_truncbits(step_pre, nb_len(imin))
  x1 <- .jlr_add12((1 - imin) * step_hi, ref)
  x2 <- .jlr_add12((len - imin) * step_hi, ref)
  a <- (start - x1[1L]) - x1[2L]
  b <- (stop - x2[1L]) - x2[2L]
  step_lo <- (b - a) / (len - 1)
  ref_lo <- a - (1 - imin) * step_lo
  out(c(ref, ref_lo), c(step_hi, step_lo), imin)
}

.jlr_add12 <- function(x, y) if (abs(y) > abs(x)) .jlr_canon(y, x) else .jlr_canon(x, y)
