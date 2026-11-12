# Elementary functions as Julia's Base computes them, for the places where a
# solver of OrdinaryDiffEq calls them itself: the C library of the platform
# rounds differently in the last bit at some arguments, and a step size
# controller carries such a bit into every later step.

# evalpoly of Julia (the @horner macro): Horner's rule with fused multiply adds
.jl_horner <- function(x, cf) {
  ex <- cf[length(cf)]
  for (i in rev(seq_len(length(cf) - 1L))) ex <- .fma(x, ex, cf[i])
  ex
}

# atan(::Float64) of Base/special/trig.jl (openlibm's s_atan.c)
.jl_atan <- function(x) {
  vapply(x, function(x) {
    absx <- abs(x)
    if (is.nan(x)) return(x)
    if (absx >= 0x1p66) return(sign(x) * 0x1.921fb54442d18p+0)
    pq <- function(x) {
      x2 <- x * x
      x4 <- x2 * x2
      p <- x2 * .jl_horner(x4, c(0x1.555555555550dp-2, 0x1.24924920083ffp-3, 0x1.745cdc54c206ep-4,
                                 0x1.10d66a0d03d51p-4, 0x1.97b4b24760debp-5, 0x1.0ad3ae322da11p-6))
      q <- x4 * .jl_horner(x4, c(-0x1.999999998ebc4p-3, -0x1.c71c6fe231671p-4, -0x1.3b0f2af749a6dp-4,
                                 -0x1.dde2d52defd9ap-5, -0x1.2b4442c6a6c2fp-5))
      p + q
    }
    if (absx < 0.4375) {
      if (absx < 0x1p-27) return(x)
      return(x - x * pq(x))
    }
    xsign <- sign(x)
    if (absx < 1.1875) {
      if (absx < 0.6875) {
        hi <- 0x1.dac670561bb4fp-2
        lo <- 0x1.a2b7f222f65e2p-56
        x <- (2 * absx - 1) / (2 + absx)
      } else {
        hi <- 0x1.921fb54442d18p-1
        lo <- 0x1.1a62633145c07p-55
        x <- (absx - 1) / (absx + 1)
      }
    } else if (absx < 2.4375) {
      hi <- 0x1.f730bd281f69bp-1
      lo <- 0x1.007887af0cbbdp-56
      x <- (absx - 1.5) / (1 + 1.5 * absx)
    } else {
      hi <- 0x1.921fb54442d18p+0
      lo <- 0x1.1a62633145c07p-54
      x <- -1 / absx
    }
    xsign * (hi - ((x * pq(x) - lo) - x))
  }, 0)
}
