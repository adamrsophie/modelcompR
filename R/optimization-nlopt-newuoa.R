# Powell's NEWUOA derivative-free trust region method as modified for NLopt
# 2.11.0 by S. G. Johnson (src/algs/newuoa/newuoa.c), ported to R: LN_NEWUOA
# (unconstrained) and LN_NEWUOA_BOUND (the bound constrained variant, whose
# trust region and Lagrange function subproblems are solved by a nested
# NLopt LD_MMA optimizer, as in the C code).
#
# The NLopt library that NLopt.jl loads (NLopt_jll) rounds every product
# and sum separately (no fused multiply add), so plain R arithmetic in the
# order of the C expressions reproduces it. Every sum is accumulated in the
# loop order of the C code.
#
# Copyright (c) 2004 M. J. D. Powell (mjdp@cam.ac.uk); C translation and
# NLopt modifications Copyright (c) 2008 Massachusetts Institute of
# Technology, released under the MIT license (see
# R/optimization-nlopt-util.R).

.nnu_min2 <- function(a, b) if (a <= b) a else b
.nnu_max2 <- function(a, b) if (a >= b) a else b
.nnu_failed <- function(ret) ret %in% c("FAILURE", "INVALID_ARGS", "OUT_OF_MEMORY", "ROUNDOFF_LIMITED", "FORCED_STOP")

# cos(angle) and sin(angle), which the C code takes in pairs: clang merges
# each pair into one __sincos_stret of the macOS arm64 libm. For |x| <= pi/4
# its sine rounds the product x * (x^2 p(x^2)) before adding x, where sin()
# fuses that multiply add, so the two differ in the last bit for some x;
# beyond pi/4 (after argument reduction) they agree. The cosine there is
# 1 + x^2 q(x^2), which R's cos() (unlike that libm) misses in the last bit
# for some tiny x. The polynomials p and q and their Horner schemes of fused
# multiply adds are those of the libm.
.nnu_sin <- function(x) {
  if (!(abs(x) <= pi / 4) || x == 0) return(sin(x))
  z <- x * x
  p <- .fma(z, 0x1.5d8fd1fd19ccdp-33, -0x1.ae5e5a9291f5dp-26)
  p <- .fma(z, p, 0x1.71de3567d48a1p-19)
  p <- .fma(z, p, -0x1.a01a019bfdf03p-13)
  p <- .fma(z, p, 0x1.111111110f7d0p-7)
  p <- .fma(z, p, -0x1.5555555555548p-3)
  x * (z * p) + x
}
.nnu_cos <- function(x) {
  if (!(abs(x) <= pi / 4)) return(cos(x))
  z <- x * x
  q <- .fma(z, -0x1.8fa49a0861a9bp-37, 0x1.1ee9d7b4e3f05p-29)
  q <- .fma(z, q, -0x1.27e4f7eac4bc6p-22)
  q <- .fma(z, q, 0x1.a01a019c844f5p-16)
  q <- .fma(z, q, -0x1.6c16c16c14f91p-10)
  q <- .fma(z, q, 0x1.555555555554bp-5)
  q <- .fma(z, q, -0.5)
  .fma(z, q, 1)
}

# the product of the model Hessian (implicit part from xpt and pq, explicit
# part from the packed lower triangle hq) with the vector v
.nnu_hessmul <- function(n, npt, xpt, pq, hq, v) {
  hd <- numeric(n)
  for (k in seq_len(npt)) {
    temp <- 0
    for (j in seq_len(n)) temp <- temp + xpt[k, j] * v[j]
    temp <- temp * pq[k]
    hd <- hd + temp * xpt[k, ]
  }
  ih <- 0L
  for (j in seq_len(n)) for (i in seq_len(j)) {
    ih <- ih + 1L
    if (i < j) hd[j] <- hd[j] + hq[ih] * v[i]
    hd[i] <- hd[i] + hq[ih] * v[j]
  }
  hd
}

# quad_model: the quadratic model at xopt + x (value and gradient)
.nnu_quad_model <- function(n, npt, xpt, pq, hq, gq, xopt, x) {
  y <- xopt + x
  hd <- numeric(n)
  for (k in seq_len(npt)) {
    temp <- 0
    for (j in seq_len(n)) temp <- temp + xpt[k, j] * y[j]
    temp <- temp * pq[k]
    hd <- hd + temp * xpt[k, ]
  }
  k <- 0L
  for (j in seq_len(n)) {
    for (i in seq_len(j - 1L)) {
      k <- k + 1L
      hd[j] <- hd[j] + hq[k] * y[i]
      hd[i] <- hd[i] + hq[k] * y[j]
    }
    k <- k + 1L
    hd[j] <- hd[j] + hq[k] * y[j]
  }
  val <- 0
  for (i in seq_len(n)) val <- val + (gq[i] + 0.5 * hd[i]) * y[i]
  list(val = val, grad = gq + hd)
}

# rho_constraint: |x|^2 <= rho^2
.nnu_rho_constraint <- function(x, rho) {
  val <- -rho * rho
  for (i in seq_along(x)) val <- val + x[i] * x[i]
  val
}

# a nested LD_MMA run of nlopt_optimize from x0 with the trust region
# constraint |x| <= delta, the bounds [lo, hi], the absolute x tolerance
# xtol and at most 1000 evaluations
.nnu_mma <- function(fg, n, delta, lo, hi, xtol, x0) {
  o <- .nlopt_create("LD_MMA", n)
  o$f <- function(x) fg(x)$val
  o$fg <- fg
  o$fc <- list(list(m = 1L, f = function(x) .nnu_rho_constraint(x, delta),
                    fg = function(x) list(val = .nnu_rho_constraint(x, delta), grad = matrix(2 * x, 1L)), tol = 0))
  o$lb <- lo; o$ub <- hi
  o$xtol_abs <- xtol
  o$maxeval <- 1000L
  .nlopt_optimize(o, x0)
}

# the box [lo, hi] of steps from xbase + xopt that stay within the trust
# region and the bounds, containing 0
.nnu_step_box <- function(n, delta, xbase, xopt, lb, ub) {
  lo <- numeric(n); hi <- numeric(n)
  for (j in seq_len(n)) {
    hi[j] <- delta
    lo[j] <- -hi[j]
    if (lo[j] < lb[j] - xbase[j] - xopt[j]) lo[j] <- lb[j] - xbase[j] - xopt[j]
    if (hi[j] > ub[j] - xbase[j] - xopt[j]) hi[j] <- ub[j] - xbase[j] - xopt[j]
    if (lo[j] > 0) lo[j] <- 0
    if (hi[j] < 0) hi[j] <- 0
  }
  list(lo = lo, hi = hi)
}

# trsapp_: an approximate minimizer of the quadratic model within the trust
# region (truncated conjugate gradients, or LD_MMA with bounds); returns
# list(ret, step, crvmin)
.nnu_trsapp <- function(n, npt, xopt, xpt, gq, hq, pq, delta, xbase, lb, ub, bounded) {
  if (bounded) {
    qm <- function(x) .nnu_quad_model(n, npt, xpt, pq, hq, gq, xopt, x)
    box <- .nnu_step_box(n, delta, xbase, xopt, lb, ub)
    r <- .nnu_mma(qm, n, delta, box$lo, box$hi, rep(1e-7 * delta, n), numeric(n))
    step <- r$x
    if (.nnu_rho_constraint(step, delta) > -1e-6 * delta * delta) {
      crv <- 0
    } else {
      g <- qm(step - xopt)$grad
      crv <- 0; gg <- 0
      for (j in seq_len(n)) {
        crv <- crv + step[j] * (g[j] - gq[j])
        gg <- gg + step[j] * step[j]
      }
      crv <- if (gg <= 1e-16 * crv) 1e16 else crv / gg
    }
    return(list(ret = r$ret, step = step, crvmin = crv))
  }
  twopi <- atan(1) * 8
  delsq <- delta * delta
  iterc <- 0L; itermax <- n; itersw <- itermax
  step <- numeric(n); hs <- numeric(n); g <- numeric(n); hd <- numeric(n)
  d <- xopt
  dd <- 0; gg <- 0; sg <- 0; shs <- 0; qred <- 0; ggbeg <- 0; tempa <- 0; tempb <- 0; bstep <- 0
  ds <- 0; ss <- 0; crvmin <- 0
  lab <- "L170"
  repeat {
    if (lab == "L170") {
      hd <- .nnu_hessmul(n, npt, xpt, pq, hq, d)
      lab <- if (iterc == 0L) "L20" else if (iterc <= itersw) "L50" else "L120"
    } else if (lab == "L20") {
      qred <- 0; dd <- 0
      step[] <- 0; hs[] <- 0
      g <- gq + hd
      d <- -g
      for (i in seq_len(n)) dd <- dd + d[i] * d[i]
      crvmin <- 0
      if (dd == 0) break
      ds <- 0; ss <- 0; gg <- dd; ggbeg <- gg
      lab <- "L40"
    } else if (lab == "L40") {
      iterc <- iterc + 1L
      temp <- delsq - ss
      bstep <- temp / (ds + sqrt(ds * ds + dd * temp))
      lab <- "L170"
    } else if (lab == "L50") {
      dhd <- 0
      for (j in seq_len(n)) dhd <- dhd + d[j] * hd[j]
      alpha <- bstep
      if (dhd > 0) {
        temp <- dhd / dd
        if (iterc == 1L) crvmin <- temp
        crvmin <- .nnu_min2(crvmin, temp)
        alpha <- .nnu_min2(alpha, gg / dhd)
      }
      qadd <- alpha * (gg - 0.5 * alpha * dhd)
      qred <- qred + qadd
      ggsav <- gg; gg <- 0
      for (i in seq_len(n)) {
        step[i] <- step[i] + alpha * d[i]
        hs[i] <- hs[i] + alpha * hd[i]
        t <- g[i] + hs[i]
        gg <- gg + t * t
      }
      if (alpha < bstep) {
        if (qadd <= qred * 0.01) break
        if (gg <= ggbeg * 1e-4) break
        if (iterc == itermax) break
        temp <- gg / ggsav
        dd <- 0; ds <- 0; ss <- 0
        for (i in seq_len(n)) {
          d[i] <- temp * d[i] - g[i] - hs[i]
          dd <- dd + d[i] * d[i]
          ds <- ds + d[i] * step[i]
          ss <- ss + step[i] * step[i]
        }
        if (ds <= 0) break
        if (ss < delsq) { lab <- "L40"; next }
      }
      crvmin <- 0
      itersw <- iterc
      lab <- "L90"
    } else if (lab == "L90") {
      if (gg <= ggbeg * 1e-4) break
      sg <- 0; shs <- 0
      for (i in seq_len(n)) {
        sg <- sg + step[i] * g[i]
        shs <- shs + step[i] * hs[i]
      }
      sgk <- sg + shs
      angtest <- sgk / sqrt(gg * delsq)
      if (angtest <= -0.99) break
      iterc <- iterc + 1L
      temp <- sqrt(delsq * gg - sgk * sgk)
      tempa <- delsq / temp
      tempb <- sgk / temp
      d <- tempa * (g + hs) - tempb * step
      lab <- "L170"
    } else {
      # L120
      dg <- 0; dhd <- 0; dhs <- 0
      for (i in seq_len(n)) {
        dg <- dg + d[i] * g[i]
        dhd <- dhd + hd[i] * d[i]
        dhs <- dhs + hd[i] * step[i]
      }
      cf <- 0.5 * (shs - dhd)
      qbeg <- sg + cf
      qsav <- qbeg; qmin <- qbeg
      isave <- 0L; iu <- 49L
      temp <- twopi / (iu + 1)
      qnew <- 0
      for (i in seq_len(iu)) {
        angle <- i * temp
        cth <- .nnu_cos(angle); sth <- .nnu_sin(angle)
        qnew <- (sg + cf * cth) * cth + (dg + dhs * cth) * sth
        if (qnew < qmin) {
          qmin <- qnew; isave <- i; tempa <- qsav
        } else if (i == isave + 1L) {
          tempb <- qnew
        }
        qsav <- qnew
      }
      if (isave == 0L) tempa <- qnew
      if (isave == iu) tempb <- qbeg
      angle <- 0
      if (tempa != tempb) {
        tempa <- tempa - qmin
        tempb <- tempb - qmin
        angle <- 0.5 * (tempa - tempb) / (tempa + tempb)
      }
      angle <- temp * (isave + angle)
      cth <- .nnu_cos(angle); sth <- .nnu_sin(angle)
      reduc <- qbeg - (sg + cf * cth) * cth - (dg + dhs * cth) * sth
      gg <- 0
      for (i in seq_len(n)) {
        step[i] <- cth * step[i] + sth * d[i]
        hs[i] <- cth * hs[i] + sth * hd[i]
        t <- g[i] + hs[i]
        gg <- gg + t * t
      }
      qred <- qred + reduc
      ratio <- reduc / qred
      if (iterc < itermax && ratio > 0.01) { lab <- "L90"; next }
      break
    }
  }
  list(ret = "SUCCESS", step = step, crvmin = crvmin)
}

# biglag_: a step d of length delta from xopt that makes the KNEW-th
# Lagrange function large (LD_MMA with bounds); returns list(ret, d, alpha)
.nnu_biglag <- function(n, npt, xopt, xpt, bmat, zmat, idz, knew, delta, xbase, lb, ub, bounded) {
  twopi <- atan(1) * 8
  delsq <- delta * delta
  nptm <- npt - n - 1L
  iterc <- 0L
  hcol <- numeric(npt)
  for (j in seq_len(nptm)) {
    temp <- zmat[knew, j]
    if (j < idz) temp <- -temp
    for (k in seq_len(npt)) {
      hcol[k] <- hcol[k] + temp * zmat[k, j]
      if (.nlopt_isinf(hcol[k])) return(list(ret = "ROUNDOFF_LIMITED"))
    }
  }
  alpha <- hcol[knew]
  dd <- 0
  d <- xpt[knew, ] - xopt
  gc <- bmat[knew, ]
  gd <- numeric(n)
  for (i in seq_len(n)) dd <- dd + d[i] * d[i]
  for (k in seq_len(npt)) {
    temp <- 0; sum <- 0
    for (j in seq_len(n)) {
      temp <- temp + xpt[k, j] * xopt[j]
      sum <- sum + xpt[k, j] * d[j]
    }
    temp <- hcol[k] * temp
    sum <- hcol[k] * sum
    gc <- gc + temp * xpt[k, ]
    gd <- gd + sum * xpt[k, ]
  }
  gg <- 0; sp <- 0; dhd <- 0
  for (i in seq_len(n)) {
    gg <- gg + gc[i] * gc[i]
    sp <- sp + d[i] * gc[i]
    dhd <- dhd + d[i] * gd[i]
  }
  scale <- delta / sqrt(dd)
  if (sp * dhd < 0) scale <- -scale
  temp <- 0
  if (sp * sp > dd * 0.99 * gg) temp <- 1
  tau <- scale * (abs(sp) + 0.5 * scale * abs(dhd))
  if (gg * delsq < tau * 0.01 * tau) temp <- 1
  d <- scale * d
  gd <- scale * gd
  s <- gc + temp * gd
  if (bounded) {
    brow <- bmat[knew, ]
    flip <- FALSE
    lag <- function(dx) {
      y <- xopt + dx
      val <- 0
      for (j in seq_len(n)) val <- val + brow[j] * y[j]
      grad <- brow
      for (i in seq_len(npt)) {
        dot <- 0
        for (j in seq_len(n)) dot <- dot + xpt[i, j] * y[j]
        val <- val + 0.5 * hcol[i] * (dot * dot)
        dot <- dot * hcol[i]
        grad <- grad + dot * xpt[i, ]
      }
      if (flip) list(val = -val, grad = -grad) else list(val = val, grad = grad)
    }
    d <- d * 0.99999
    box <- .nnu_step_box(n, delta, xbase, xopt, lb, ub)
    for (j in seq_len(n)) {
      if (d[j] < box$lo[j]) d[j] <- box$lo[j] else if (d[j] > box$hi[j]) d[j] <- box$hi[j]
    }
    flip <- lag(d)$val > 0
    r <- .nnu_mma(lag, n, delta, box$lo, box$hi, rep(1e-5 * delta, n), d)
    return(list(ret = r$ret, d = r$x, alpha = alpha))
  }
  tempa <- 0; tempb <- 0
  repeat {
    # L80
    iterc <- iterc + 1L
    dd <- 0; sp <- 0; ss <- 0
    for (i in seq_len(n)) {
      dd <- dd + d[i] * d[i]
      sp <- sp + d[i] * s[i]
      ss <- ss + s[i] * s[i]
    }
    temp <- dd * ss - sp * sp
    if (temp <= dd * 1e-8 * ss) break
    denom <- sqrt(temp)
    s <- (dd * s - sp * d) / denom
    w <- numeric(n)
    for (k in seq_len(npt)) {
      sum <- 0
      for (j in seq_len(n)) sum <- sum + xpt[k, j] * s[j]
      sum <- hcol[k] * sum
      w <- w + sum * xpt[k, ]
    }
    cf1 <- 0; cf2 <- 0; cf3 <- 0; cf4 <- 0; cf5 <- 0
    for (i in seq_len(n)) {
      cf1 <- cf1 + s[i] * w[i]
      cf2 <- cf2 + d[i] * gc[i]
      cf3 <- cf3 + s[i] * gc[i]
      cf4 <- cf4 + d[i] * gd[i]
      cf5 <- cf5 + s[i] * gd[i]
    }
    cf1 <- 0.5 * cf1
    cf4 <- 0.5 * cf4 - cf1
    taubeg <- cf1 + cf2 + cf4
    taumax <- taubeg; tauold <- taubeg
    isave <- 0L; iu <- 49L
    temp <- twopi / (iu + 1)
    for (i in seq_len(iu)) {
      angle <- i * temp
      cth <- .nnu_cos(angle); sth <- .nnu_sin(angle)
      tau <- cf1 + (cf2 + cf4 * cth) * cth + (cf3 + cf5 * cth) * sth
      if (abs(tau) > abs(taumax)) {
        taumax <- tau; isave <- i; tempa <- tauold
      } else if (i == isave + 1L) {
        tempb <- tau
      }
      tauold <- tau
    }
    if (isave == 0L) tempa <- tau
    if (isave == iu) tempb <- taubeg
    step <- 0
    if (tempa != tempb) {
      tempa <- tempa - taumax
      tempb <- tempb - taumax
      step <- 0.5 * (tempa - tempb) / (tempa + tempb)
    }
    angle <- temp * (isave + step)
    cth <- .nnu_cos(angle); sth <- .nnu_sin(angle)
    tau <- cf1 + (cf2 + cf4 * cth) * cth + (cf3 + cf5 * cth) * sth
    d <- cth * d + sth * s
    gd <- cth * gd + sth * w
    s <- gc + gd
    if (abs(tau) <= abs(taubeg) * 1.1) break
    if (iterc >= n) break
  }
  list(ret = "SUCCESS", d = d, alpha = alpha)
}

# bigden_: a revision of the step d that makes the denominator of the
# updating formula large; returns list(ret, d, w, vlag, beta)
.nnu_bigden <- function(n, npt, xopt, xpt, bmat, zmat, idz, kopt, knew, d, xbase, lb, ub, bounded) {
  ndim <- npt + n
  twopi <- atan(1) * 8
  nptm <- npt - n - 1L
  roundoff <- list(ret = "ROUNDOFF_LIMITED")
  w <- numeric(ndim)
  for (j in seq_len(nptm)) {
    temp <- zmat[knew, j]
    if (j < idz) temp <- -temp
    w[n + seq_len(npt)] <- w[n + seq_len(npt)] + temp * zmat[, j]
  }
  alpha <- w[n + knew]
  dd <- 0; ds <- 0; ss <- 0; xoptsq <- 0
  s <- numeric(n)
  for (i in seq_len(n)) {
    dd <- dd + d[i] * d[i]
    s[i] <- xpt[knew, i] - xopt[i]
    ds <- ds + d[i] * s[i]
    ss <- ss + s[i] * s[i]
    xoptsq <- xoptsq + xopt[i] * xopt[i]
  }
  if (ds * ds > dd * 0.99 * ss) {
    ksav <- knew
    dtest <- ds * ds / ss
    for (k in seq_len(npt)) if (k != kopt) {
      dstemp <- 0; sstemp <- 0
      for (i in seq_len(n)) {
        diff <- xpt[k, i] - xopt[i]
        dstemp <- dstemp + d[i] * diff
        sstemp <- sstemp + diff * diff
      }
      if (sstemp == 0) return(roundoff)
      if (dstemp * dstemp / sstemp < dtest) {
        ksav <- k
        dtest <- dstemp * dstemp / sstemp
        ds <- dstemp
        ss <- sstemp
      }
    }
    s <- xpt[ksav, ] - xopt
  }
  ssden <- dd * ss - ds * ds
  iterc <- 0L
  densav <- 0
  den <- numeric(9); denex <- numeric(9); par <- numeric(9)
  wvec <- matrix(0, ndim, 5L); prod <- matrix(0, ndim, 5L)
  vlag <- numeric(ndim)
  tempa <- 0; tempb <- 0; tempc <- 0
  rotate <- function(par) {
    for (j in c(4L, 6L, 8L)) {
      par[j] <- par[2] * par[j - 2L] - par[3] * par[j - 1L]
      par[j + 1L] <- par[2] * par[j - 1L] + par[3] * par[j - 2L]
    }
    par
  }
  repeat {
    # L70
    iterc <- iterc + 1L
    if (ssden < 0) return(roundoff)
    temp <- 1 / sqrt(ssden)
    xoptd <- 0; xopts <- 0
    for (i in seq_len(n)) {
      s[i] <- temp * (dd * s[i] - ds * d[i])
      xoptd <- xoptd + xopt[i] * d[i]
      if (.nlopt_isinf(s[i])) return(roundoff)
      xopts <- xopts + xopt[i] * s[i]
    }
    tempa <- 0.5 * xoptd * xoptd
    tempb <- 0.5 * xopts * xopts
    den[1] <- dd * (xoptsq + 0.5 * dd) + tempa + tempb
    den[2] <- 2 * xoptd * dd
    den[3] <- 2 * xopts * dd
    den[4] <- tempa - tempb
    den[5] <- xoptd * xopts
    den[6:9] <- 0
    for (k in seq_len(npt)) {
      tempa <- 0; tempb <- 0; tempc <- 0
      for (i in seq_len(n)) {
        tempa <- tempa + xpt[k, i] * d[i]
        tempb <- tempb + xpt[k, i] * s[i]
        tempc <- tempc + xpt[k, i] * xopt[i]
      }
      wvec[k, 1] <- 0.25 * (tempa * tempa + tempb * tempb)
      wvec[k, 2] <- tempa * tempc
      wvec[k, 3] <- tempb * tempc
      wvec[k, 4] <- 0.25 * (tempa * tempa - tempb * tempb)
      wvec[k, 5] <- 0.5 * tempa * tempb
    }
    for (i in seq_len(n)) {
      ip <- i + npt
      wvec[ip, ] <- c(0, d[i], s[i], 0, 0)
    }
    for (jc in 1:5) {
      nw <- if (jc == 2L || jc == 3L) ndim else npt
      prod[seq_len(npt), jc] <- 0
      for (j in seq_len(nptm)) {
        sum <- 0
        for (k in seq_len(npt)) sum <- sum + zmat[k, j] * wvec[k, jc]
        if (j < idz) sum <- -sum
        prod[seq_len(npt), jc] <- prod[seq_len(npt), jc] + sum * zmat[, j]
      }
      if (nw == ndim) for (k in seq_len(npt)) {
        sum <- 0
        for (j in seq_len(n)) sum <- sum + bmat[k, j] * wvec[npt + j, jc]
        prod[k, jc] <- prod[k, jc] + sum
      }
      for (j in seq_len(n)) {
        sum <- 0
        for (i in seq_len(nw)) sum <- sum + bmat[i, j] * wvec[i, jc]
        prod[npt + j, jc] <- sum
      }
    }
    for (k in seq_len(ndim)) {
      sum <- 0
      for (i in 1:5) {
        par[i] <- 0.5 * prod[k, i] * wvec[k, i]
        sum <- sum + par[i]
      }
      den[1] <- den[1] - par[1] - sum
      tempa <- prod[k, 1] * wvec[k, 2] + prod[k, 2] * wvec[k, 1]
      tempb <- prod[k, 2] * wvec[k, 4] + prod[k, 4] * wvec[k, 2]
      tempc <- prod[k, 3] * wvec[k, 5] + prod[k, 5] * wvec[k, 3]
      den[2] <- den[2] - tempa - 0.5 * (tempb + tempc)
      den[6] <- den[6] - 0.5 * (tempb - tempc)
      tempa <- prod[k, 1] * wvec[k, 3] + prod[k, 3] * wvec[k, 1]
      tempb <- prod[k, 2] * wvec[k, 5] + prod[k, 5] * wvec[k, 2]
      tempc <- prod[k, 3] * wvec[k, 4] + prod[k, 4] * wvec[k, 3]
      den[3] <- den[3] - tempa - 0.5 * (tempb - tempc)
      den[7] <- den[7] - 0.5 * (tempb + tempc)
      tempa <- prod[k, 1] * wvec[k, 4] + prod[k, 4] * wvec[k, 1]
      den[4] <- den[4] - tempa - par[2] + par[3]
      tempa <- prod[k, 1] * wvec[k, 5] + prod[k, 5] * wvec[k, 1]
      tempb <- prod[k, 2] * wvec[k, 3] + prod[k, 3] * wvec[k, 2]
      den[5] <- den[5] - tempa - 0.5 * tempb
      den[8] <- den[8] - par[4] + par[5]
      tempa <- prod[k, 4] * wvec[k, 5] + prod[k, 5] * wvec[k, 4]
      den[9] <- den[9] - 0.5 * tempa
    }
    sum <- 0
    for (i in 1:5) {
      par[i] <- 0.5 * (prod[knew, i] * prod[knew, i])
      sum <- sum + par[i]
    }
    denex[1] <- alpha * den[1] + par[1] + sum
    tempa <- 2 * prod[knew, 1] * prod[knew, 2]
    tempb <- prod[knew, 2] * prod[knew, 4]
    tempc <- prod[knew, 3] * prod[knew, 5]
    denex[2] <- alpha * den[2] + tempa + tempb + tempc
    denex[6] <- alpha * den[6] + tempb - tempc
    tempa <- 2 * prod[knew, 1] * prod[knew, 3]
    tempb <- prod[knew, 2] * prod[knew, 5]
    tempc <- prod[knew, 3] * prod[knew, 4]
    denex[3] <- alpha * den[3] + tempa + tempb - tempc
    denex[7] <- alpha * den[7] + tempb + tempc
    tempa <- 2 * prod[knew, 1] * prod[knew, 4]
    denex[4] <- alpha * den[4] + tempa + par[2] - par[3]
    tempa <- 2 * prod[knew, 1] * prod[knew, 5]
    denex[5] <- alpha * den[5] + tempa + prod[knew, 2] * prod[knew, 3]
    denex[8] <- alpha * den[8] + par[4] - par[5]
    denex[9] <- alpha * den[9] + prod[knew, 4] * prod[knew, 5]
    sum <- denex[1] + denex[2] + denex[4] + denex[6] + denex[8]
    denold <- sum; denmax <- sum
    isave <- 0L; iu <- 49L
    temp <- twopi / (iu + 1)
    par[1] <- 1
    for (i in seq_len(iu)) {
      angle <- i * temp
      par[2] <- .nnu_cos(angle); par[3] <- .nnu_sin(angle)
      par <- rotate(par)
      sumold <- sum
      sum <- 0
      for (j in 1:9) sum <- sum + denex[j] * par[j]
      if (abs(sum) > abs(denmax)) {
        denmax <- sum; isave <- i; tempa <- sumold
      } else if (i == isave + 1L) {
        tempb <- sum
      }
    }
    if (isave == 0L) tempa <- sum
    if (isave == iu) tempb <- denold
    step <- 0
    if (tempa != tempb) {
      tempa <- tempa - denmax
      tempb <- tempb - denmax
      step <- 0.5 * (tempa - tempb) / (tempa + tempb)
    }
    angle <- temp * (isave + step)
    par[2] <- .nnu_cos(angle); par[3] <- .nnu_sin(angle)
    par <- rotate(par)
    beta <- 0; denmax <- 0
    for (j in 1:9) {
      beta <- beta + den[j] * par[j]
      denmax <- denmax + denex[j] * par[j]
    }
    for (k in seq_len(ndim)) {
      vlag[k] <- 0
      for (j in 1:5) vlag[k] <- vlag[k] + prod[k, j] * par[j]
    }
    tau <- vlag[knew]
    dd <- 0; tempa <- 0; tempb <- 0
    for (i in seq_len(n)) {
      d[i] <- par[2] * d[i] + par[3] * s[i]
      w[i] <- xopt[i] + d[i]
      dd <- dd + d[i] * d[i]
      tempa <- tempa + d[i] * w[i]
      tempb <- tempb + w[i] * w[i]
    }
    if (iterc >= n) break
    if (iterc > 1L) densav <- .nnu_max2(densav, denold)
    if (abs(denmax) <= abs(densav) * 1.1) break
    densav <- denmax
    for (i in seq_len(n)) {
      temp <- tempa * xopt[i] + tempb * d[i] - vlag[npt + i]
      s[i] <- tau * bmat[knew, i] + alpha * temp
    }
    for (k in seq_len(npt)) {
      sum <- 0
      for (j in seq_len(n)) sum <- sum + xpt[k, j] * w[j]
      if (.nlopt_isinf(tau * w[n + k]) || .nlopt_isinf(alpha * vlag[k])) return(roundoff)
      temp <- (tau * w[n + k] - alpha * vlag[k]) * sum
      s <- s + temp * xpt[k, ]
    }
    ss <- 0; ds <- 0
    for (i in seq_len(n)) {
      ss <- ss + s[i] * s[i]
      ds <- ds + d[i] * s[i]
    }
    ssden <- dd * ss - ds * ds
    if (!(ssden >= dd * 1e-8 * ss)) break
  }
  # L340
  if (bounded) for (k in seq_len(n)) {
    if (d[k] > ub[k] - xbase[k] - xopt[k]) d[k] <- ub[k] - xbase[k] - xopt[k]
    else if (d[k] < lb[k] - xbase[k] - xopt[k]) d[k] <- lb[k] - xbase[k] - xopt[k]
  }
  for (k in seq_len(ndim)) {
    w[k] <- 0
    for (j in 1:5) w[k] <- w[k] + wvec[k, j] * par[j]
  }
  vlag[kopt] <- vlag[kopt] + 1
  list(ret = "SUCCESS", d = d, w = w, vlag = vlag, beta = beta)
}

# update_: the updating of bmat and zmat when xpt[knew, ] is replaced;
# returns list(bmat, zmat, idz, vlag, w)
.nnu_update <- function(n, npt, bmat, zmat, idz, vlag, beta, knew) {
  ndim <- npt + n
  nptm <- npt - n - 1L
  jl <- 1L
  for (j in seq_len(nptm)[-1L]) {
    if (j == idz) {
      jl <- idz
    } else if (zmat[knew, j] != 0) {
      temp <- sqrt(zmat[knew, jl] * zmat[knew, jl] + zmat[knew, j] * zmat[knew, j])
      tempa <- zmat[knew, jl] / temp
      tempb <- zmat[knew, j] / temp
      t <- tempa * zmat[, jl] + tempb * zmat[, j]
      zmat[, j] <- tempa * zmat[, j] - tempb * zmat[, jl]
      zmat[, jl] <- t
      zmat[knew, j] <- 0
    }
  }
  tempa <- zmat[knew, 1]
  if (idz >= 2L) tempa <- -tempa
  if (jl > 1L) tempb <- zmat[knew, jl]
  w <- numeric(ndim)
  w[seq_len(npt)] <- tempa * zmat[, 1]
  if (jl > 1L) w[seq_len(npt)] <- w[seq_len(npt)] + tempb * zmat[, jl]
  alpha <- w[knew]
  tau <- vlag[knew]
  tausq <- tau * tau
  denom <- alpha * beta + tausq
  vlag[knew] <- vlag[knew] - 1
  iflag <- 0L
  ip <- seq_len(npt)
  if (jl == 1L) {
    temp <- sqrt(abs(denom))
    tempb <- tempa / temp
    tempa <- tau / temp
    zmat[, 1] <- tempa * zmat[, 1] - tempb * vlag[ip]
    if (idz == 1L && temp < 0) idz <- 2L
    if (idz >= 2L && temp >= 0) iflag <- 1L
  } else {
    ja <- 1L
    if (beta >= 0) ja <- jl
    jb <- jl + 1L - ja
    temp <- zmat[knew, jb] / denom
    tempa <- temp * beta
    tempb <- temp * tau
    temp <- zmat[knew, ja]
    scala <- 1 / sqrt(abs(beta) * temp * temp + tausq)
    scalb <- scala * sqrt(abs(denom))
    zmat[, ja] <- scala * (tau * zmat[, ja] - temp * vlag[ip])
    zmat[, jb] <- scalb * (zmat[, jb] - tempa * w[ip] - tempb * vlag[ip])
    if (denom <= 0) {
      if (beta < 0) idz <- idz + 1L
      if (beta >= 0) iflag <- 1L
    }
  }
  if (iflag == 1L) {
    idz <- idz - 1L
    temp <- zmat[, 1]
    zmat[, 1] <- zmat[, idz]
    zmat[, idz] <- temp
  }
  for (j in seq_len(n)) {
    jp <- npt + j
    w[jp] <- bmat[knew, j]
    tempa <- (alpha * vlag[jp] - tau * w[jp]) / denom
    tempb <- (-beta * w[jp] - tau * vlag[jp]) / denom
    for (i in seq_len(jp)) {
      bmat[i, j] <- bmat[i, j] + tempa * vlag[i] + tempb * w[i]
      if (i > npt) bmat[jp, i - npt] <- bmat[i, j]
    }
  }
  list(bmat = bmat, zmat = zmat, idz = idz, vlag = vlag, w = w)
}

# newuoa: minimize fn from x with npt interpolation conditions and initial
# trust region radius rhobeg; lb and ub are NULL for LN_NEWUOA. Returns
# list(ret, x, minf).
.nlopt_newuoa <- function(n, npt, x, lb, ub, rhobeg, stop, fn) {
  bounded <- !is.null(lb) && !is.null(ub)
  np <- n + 1L
  nptm <- npt - np
  if (n < 2L) {
    stop$msg <- sprintf("dimension %d must be >= 2", n)
    return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  }
  if (npt < n + 2L || npt > ((n + 2L) * np) %/% 2L) {
    stop$msg <- sprintf("invalid # of interpolation conditions %d", npt)
    return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  }
  ndim <- npt + n
  rhoend <- stop$xtol_rel * rhobeg
  if (!is.null(stop$xtol_abs)) {
    xa <- rep_len(stop$xtol_abs, n)
    for (j in seq_len(n)) if (rhoend < xa[j]) rhoend <- xa[j]
  }
  half <- 0.5; tenth <- 0.1
  nh <- (n * np) %/% 2L
  xbase <- x
  xpt <- matrix(0, npt, n); bmat <- matrix(0, ndim, n); zmat <- matrix(0, npt, nptm)
  hq <- numeric(nh); pq <- numeric(npt); gq <- numeric(n); fval <- numeric(npt)
  xopt <- numeric(n); xnew <- numeric(n); d <- numeric(n); vlag <- numeric(ndim)
  w <- numeric(2L * npt + n)
  rhosq <- rhobeg * rhobeg
  recip <- 1 / rhosq
  reciq <- sqrt(half) / rhosq
  nf <- 0L; nfm <- 0L; nfmm <- 0L
  rc <- "SUCCESS"
  f <- Inf; fopt <- Inf; fbeg <- 0; fsave <- 0; kopt <- 1L; knew <- 0L; ksave <- 0L
  ipt <- 0L; jpt <- 0L; xipt <- 0; xjpt <- 0
  rho <- 0; delta <- 0; idz <- 1L; diffa <- 0; diffb <- 0; diffc <- 0; itest <- 0L; xoptsq <- 0
  nfsav <- 0L; dsq <- 0; dnorm <- 0; crvmin <- 0; ratio <- 0; alpha <- 0; beta <- 0; dstep <- 0
  vquad <- 0; diff <- 0
  clamp <- function(x) {
    if (bounded) for (j in seq_len(n)) {
      if (x[j] < lb[j]) x[j] <- lb[j] else if (x[j] > ub[j]) x[j] <- ub[j]
    }
    x
  }
  lab <- "L50"
  repeat {
    switch(lab,
      L50 = {
        nfm <- nf
        nfmm <- nf - n
        nf <- nf + 1L
        if (nfm <= 2L * n) {
          if (nfm >= 1L && nfm <= n) xpt[nf, nfm] <- rhobeg
          else if (nfm > n) xpt[nf, nfmm] <- -rhobeg
        } else {
          itemp <- (nfmm - 1L) %/% n
          jpt <- nfm - itemp * n - n
          ipt <- jpt + itemp
          if (ipt > n) {
            itemp <- jpt
            jpt <- ipt - n
            ipt <- itemp
          }
          xipt <- rhobeg
          if (fval[ipt + np] < fval[ipt + 1L]) xipt <- -xipt
          xjpt <- rhobeg
          if (fval[jpt + np] < fval[jpt + 1L]) xjpt <- -xjpt
          xpt[nf, ipt] <- xipt
          xpt[nf, jpt] <- xjpt
        }
        x <- clamp(xpt[nf, ] + xbase)
        lab <- "L310"
      },
      L70 = {
        fval[nf] <- f
        if (nf == 1L) {
          fbeg <- f; fopt <- f; kopt <- 1L
        } else if (f < fopt) {
          fopt <- f; kopt <- nf
        }
        if (nfm <= 2L * n) {
          if (nfm >= 1L && nfm <= n) {
            gq[nfm] <- (f - fbeg) / rhobeg
            if (npt < nf + n) {
              bmat[1, nfm] <- -1 / rhobeg
              bmat[nf, nfm] <- 1 / rhobeg
              bmat[npt + nfm, nfm] <- -half * rhosq
            }
          } else if (nfm > n) {
            bmat[nf - n, nfmm] <- half / rhobeg
            bmat[nf, nfmm] <- -half / rhobeg
            zmat[1, nfmm] <- -reciq - reciq
            zmat[nf - n, nfmm] <- reciq
            zmat[nf, nfmm] <- reciq
            ih <- (nfmm * (nfmm + 1L)) %/% 2L
            temp <- (fbeg - f) / rhobeg
            hq[ih] <- (gq[nfmm] - temp) / rhobeg
            gq[nfmm] <- half * (gq[nfmm] + temp)
          }
        } else {
          ih <- (ipt * (ipt - 1L)) %/% 2L + jpt
          if (xipt < 0) ipt <- ipt + n
          if (xjpt < 0) jpt <- jpt + n
          zmat[1, nfmm] <- recip
          zmat[nf, nfmm] <- recip
          zmat[ipt + 1L, nfmm] <- -recip
          zmat[jpt + 1L, nfmm] <- -recip
          hq[ih] <- (fbeg - fval[ipt + 1L] - fval[jpt + 1L] + f) / (xipt * xjpt)
        }
        if (nf < npt) {
          lab <- "L50"
        } else {
          rho <- rhobeg; delta <- rho; idz <- 1L
          diffa <- 0; diffb <- 0; itest <- 0L; xoptsq <- 0
          for (i in seq_len(n)) {
            xopt[i] <- xpt[kopt, i]
            xoptsq <- xoptsq + xopt[i] * xopt[i]
          }
          lab <- "L90"
        }
      },
      L90 = {
        nfsav <- nf
        lab <- "L100"
      },
      L100 = {
        knew <- 0L
        tr <- .nnu_trsapp(n, npt, xopt, xpt, gq, hq, pq, delta, xbase, lb, ub, bounded)
        if (.nnu_failed(tr$ret)) {
          rc <- tr$ret; lab <- "L530"
        } else {
          d <- tr$step; crvmin <- tr$crvmin
          dsq <- 0
          for (i in seq_len(n)) dsq <- dsq + d[i] * d[i]
          dnorm <- .nnu_min2(delta, sqrt(dsq))
          lab <- "L120"
          if (dnorm < half * rho) {
            knew <- -1L
            delta <- tenth * delta
            ratio <- -1
            if (delta <= rho * 1.5) delta <- rho
            lab <- if (nf <= nfsav + 2L) "L460" else
              if (crvmin * 0.125 * rho * rho <= .nnu_max2(.nnu_max2(diffa, diffb), diffc)) "L460" else "L490"
          }
        }
      },
      L120 = {
        if (dsq <= xoptsq * 0.001) {
          tempq <- xoptsq * 0.25
          for (k in seq_len(npt)) {
            sum <- 0
            for (i in seq_len(n)) sum <- sum + xpt[k, i] * xopt[i]
            temp <- pq[k] * sum
            sum <- sum - half * xoptsq
            w[npt + k] <- sum
            for (i in seq_len(n)) {
              gq[i] <- gq[i] + temp * xpt[k, i]
              xpt[k, i] <- xpt[k, i] - half * xopt[i]
              vlag[i] <- bmat[k, i]
              w[i] <- sum * xpt[k, i] + tempq * xopt[i]
              ip <- npt + i
              for (j in seq_len(i)) bmat[ip, j] <- bmat[ip, j] + vlag[i] * w[j] + w[i] * vlag[j]
            }
          }
          for (k in seq_len(nptm)) {
            sumz <- 0
            for (i in seq_len(npt)) {
              sumz <- sumz + zmat[i, k]
              w[i] <- w[npt + i] * zmat[i, k]
            }
            for (j in seq_len(n)) {
              sum <- tempq * sumz * xopt[j]
              for (i in seq_len(npt)) sum <- sum + w[i] * xpt[i, j]
              vlag[j] <- sum
              if (k < idz) sum <- -sum
              bmat[seq_len(npt), j] <- bmat[seq_len(npt), j] + sum * zmat[, k]
            }
            for (i in seq_len(n)) {
              ip <- i + npt
              temp <- vlag[i]
              if (k < idz) temp <- -temp
              for (j in seq_len(i)) bmat[ip, j] <- bmat[ip, j] + temp * vlag[j]
            }
          }
          ih <- 0L
          for (j in seq_len(n)) {
            w[j] <- 0
            for (k in seq_len(npt)) {
              w[j] <- w[j] + pq[k] * xpt[k, j]
              xpt[k, j] <- xpt[k, j] - half * xopt[j]
            }
            for (i in seq_len(j)) {
              ih <- ih + 1L
              if (i < j) gq[j] <- gq[j] + hq[ih] * xopt[i]
              gq[i] <- gq[i] + hq[ih] * xopt[j]
              hq[ih] <- hq[ih] + w[i] * xopt[j] + xopt[i] * w[j]
              bmat[npt + i, j] <- bmat[npt + j, i]
            }
          }
          xbase <- xbase + xopt
          xopt[] <- 0
          xoptsq <- 0
        }
        lab <- "L290"
        if (knew > 0L) {
          bl <- .nnu_biglag(n, npt, xopt, xpt, bmat, zmat, idz, knew, dstep, xbase, lb, ub, bounded)
          if (.nnu_failed(bl$ret)) {
            rc <- bl$ret; lab <- "L530"; next
          }
          d <- bl$d; alpha <- bl$alpha
        }
        for (k in seq_len(npt)) {
          suma <- 0; sumb <- 0; sum <- 0
          for (j in seq_len(n)) {
            suma <- suma + xpt[k, j] * d[j]
            sumb <- sumb + xpt[k, j] * xopt[j]
            sum <- sum + bmat[k, j] * d[j]
          }
          w[k] <- suma * (half * suma + sumb)
          vlag[k] <- sum
        }
        beta <- 0
        for (k in seq_len(nptm)) {
          sum <- 0
          for (i in seq_len(npt)) sum <- sum + zmat[i, k] * w[i]
          if (k < idz) {
            beta <- beta + sum * sum
            sum <- -sum
          } else {
            beta <- beta - sum * sum
          }
          vlag[seq_len(npt)] <- vlag[seq_len(npt)] + sum * zmat[, k]
        }
        bsum <- 0; dx <- 0
        for (j in seq_len(n)) {
          sum <- 0
          for (i in seq_len(npt)) sum <- sum + w[i] * bmat[i, j]
          bsum <- bsum + sum * d[j]
          jp <- npt + j
          for (k in seq_len(n)) sum <- sum + bmat[jp, k] * d[k]
          vlag[jp] <- sum
          bsum <- bsum + sum * d[j]
          dx <- dx + d[j] * xopt[j]
        }
        beta <- dx * dx + dsq * (xoptsq + dx + dx + half * dsq) + beta - bsum
        vlag[kopt] <- vlag[kopt] + 1
        if (knew > 0L) {
          d1 <- vlag[knew]
          if (d1 == 0) {
            rc <- "ROUNDOFF_LIMITED"; lab <- "L530"; next
          }
          temp <- 1 + alpha * beta / (d1 * d1)
          if (abs(temp) <= 0.8) {
            bd <- .nnu_bigden(n, npt, xopt, xpt, bmat, zmat, idz, kopt, knew, d, xbase, lb, ub, bounded)
            if (.nnu_failed(bd$ret)) {
              rc <- bd$ret; lab <- "L530"; next
            }
            d <- bd$d; w[seq_len(ndim)] <- bd$w; vlag <- bd$vlag; beta <- bd$beta
          }
        }
      },
      L290 = {
        xnew <- xopt + d
        x <- clamp(xbase + xnew)
        nf <- nf + 1L
        lab <- "L310"
      },
      L310 = {
        if (.nlopt_stop_forced(stop)) {
          rc <- "FORCED_STOP"
        } else if (stop$nevals > 0L) {
          if (.nlopt_stop_evals(stop)) rc <- "MAXEVAL_REACHED"
          else if (.nlopt_stop_time(stop)) rc <- "MAXTIME_REACHED"
        }
        if (rc != "SUCCESS") { lab <- "L530"; next }
        stop$nevals <- stop$nevals + 1L
        f <- fn(x)
        if (f < stop$minf_max) {
          rc <- "STOPVAL_REACHED"; lab <- "L530"; next
        }
        if (nf <= npt) { lab <- "L70"; next }
        if (knew == -1L) { lab <- "L530"; next }
        vquad <- 0
        ih <- 0L
        for (j in seq_len(n)) {
          vquad <- vquad + d[j] * gq[j]
          for (i in seq_len(j)) {
            ih <- ih + 1L
            temp <- d[i] * xnew[j] + d[j] * xopt[i]
            if (i == j) temp <- half * temp
            vquad <- vquad + temp * hq[ih]
          }
        }
        for (k in seq_len(npt)) vquad <- vquad + pq[k] * w[k]
        diff <- f - fopt - vquad
        diffc <- diffb
        diffb <- diffa
        diffa <- abs(diff)
        if (dnorm > rho) nfsav <- nf
        fsave <- fopt
        if (f < fopt) {
          fopt <- f
          xoptsq <- 0
          for (i in seq_len(n)) {
            xopt[i] <- xnew[i]
            xoptsq <- xoptsq + xopt[i] * xopt[i]
          }
          if (.nlopt_stop_ftol(stop, fopt, fsave)) {
            rc <- "FTOL_REACHED"; lab <- "L530"; next
          }
        }
        ksave <- knew
        if (knew > 0L) { lab <- "L410"; next }
        if (vquad >= 0) { lab <- "L530"; next }
        ratio <- (f - fsave) / vquad
        if (ratio <= tenth) {
          delta <- half * dnorm
        } else if (ratio <= 0.7) {
          delta <- .nnu_max2(half * delta, dnorm)
        } else {
          delta <- .nnu_max2(half * delta, dnorm + dnorm)
        }
        if (delta <= rho * 1.5) delta <- rho
        t <- .nnu_max2(tenth * delta, rho)
        rhosq <- t * t
        ktemp <- 0L
        detrat <- 0
        if (f >= fsave) {
          ktemp <- kopt
          detrat <- 1
        }
        for (k in seq_len(npt)) {
          hdiag <- 0
          for (j in seq_len(nptm)) {
            temp <- if (j < idz) -1 else 1
            hdiag <- hdiag + temp * (zmat[k, j] * zmat[k, j])
          }
          temp <- abs(beta * hdiag + vlag[k] * vlag[k])
          distsq <- 0
          for (j in seq_len(n)) {
            t <- xpt[k, j] - xopt[j]
            distsq <- distsq + t * t
          }
          if (distsq > rhosq) {
            t <- distsq / rhosq
            temp <- temp * (t * (t * t))
          }
          if (temp > detrat && k != ktemp) {
            detrat <- temp
            knew <- k
          }
        }
        lab <- if (knew == 0L) "L460" else "L410"
      },
      L410 = {
        up <- .nnu_update(n, npt, bmat, zmat, idz, vlag, beta, knew)
        bmat <- up$bmat; zmat <- up$zmat; idz <- up$idz; vlag <- up$vlag; w[seq_len(ndim)] <- up$w
        fval[knew] <- f
        ih <- 0L
        for (i in seq_len(n)) {
          temp <- pq[knew] * xpt[knew, i]
          for (j in seq_len(i)) {
            ih <- ih + 1L
            hq[ih] <- hq[ih] + temp * xpt[knew, j]
          }
        }
        pq[knew] <- 0
        for (j in seq_len(nptm)) {
          temp <- diff * zmat[knew, j]
          if (j < idz) temp <- -temp
          pq <- pq + temp * zmat[, j]
        }
        gqsq <- 0
        for (i in seq_len(n)) {
          gq[i] <- gq[i] + diff * bmat[knew, i]
          gqsq <- gqsq + gq[i] * gq[i]
          xpt[knew, i] <- xnew[i]
        }
        if (ksave == 0L && delta == rho) {
          if (abs(ratio) > 0.01) {
            itest <- 0L
          } else {
            vlag[seq_len(npt)] <- fval - fval[kopt]
            gisq <- 0
            for (i in seq_len(n)) {
              sum <- 0
              for (k in seq_len(npt)) sum <- sum + bmat[k, i] * vlag[k]
              gisq <- gisq + sum * sum
              w[i] <- sum
            }
            itest <- itest + 1L
            if (gqsq < gisq * 100) itest <- 0L
            if (itest >= 3L) {
              gq <- w[seq_len(n)]
              hq[] <- 0
              for (j in seq_len(nptm)) {
                w[j] <- 0
                for (k in seq_len(npt)) w[j] <- w[j] + vlag[k] * zmat[k, j]
                if (j < idz) w[j] <- -w[j]
              }
              for (k in seq_len(npt)) {
                pq[k] <- 0
                for (j in seq_len(nptm)) pq[k] <- pq[k] + zmat[k, j] * w[j]
              }
              itest <- 0L
            }
          }
        }
        if (f < fsave) kopt <- knew
        if (f <= fsave + tenth * vquad || ksave > 0L) {
          lab <- "L100"
        } else {
          knew <- 0L
          lab <- "L460"
        }
      },
      L460 = {
        distsq <- delta * 4 * delta
        for (k in seq_len(npt)) {
          sum <- 0
          for (j in seq_len(n)) {
            t <- xpt[k, j] - xopt[j]
            sum <- sum + t * t
          }
          if (sum > distsq) {
            knew <- k
            distsq <- sum
          }
        }
        if (knew > 0L) {
          dstep <- .nnu_max2(.nnu_min2(tenth * sqrt(distsq), half * delta), rho)
          dsq <- dstep * dstep
          lab <- "L120"
        } else if (ratio > 0 || .nnu_max2(delta, dnorm) > rho) {
          lab <- "L100"
        } else {
          lab <- "L490"
        }
      },
      L490 = {
        if (rho > rhoend) {
          delta <- half * rho
          ratio <- rho / rhoend
          if (ratio <= 16) rho <- rhoend
          else if (ratio <= 250) rho <- sqrt(ratio) * rhoend
          else rho <- tenth * rho
          delta <- .nnu_max2(delta, rho)
          lab <- "L90"
        } else if (knew == -1L) {
          lab <- "L290"
        } else {
          rc <- "XTOL_REACHED"
          lab <- "L530"
        }
      },
      L530 = {
        if (fopt <= f) {
          x <- xbase + xopt
          f <- fopt
        }
        return(list(ret = rc, x = x, minf = f))
      })
  }
}

# the LN_NEWUOA and LN_NEWUOA_BOUND cases of nlopt_optimize_: 2n + 1
# interpolation conditions and the smallest initial step as rhobeg
.nlopt_newuoa_minimize <- function(o, x, stop, bounded) {
  n <- length(x)
  step <- min(abs(.nlopt_dx(o, x)))
  .nlopt_newuoa(n, 2L * n + 1L, x, if (bounded) o$lb, if (bounded) o$ub, step, stop, o$f)
}
