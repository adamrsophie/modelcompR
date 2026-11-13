# SLSQP (sequential least squares programming) of Dieter Kraft as in NLopt
# 2.11.0 (src/algs/slsqp/slsqp.c, the f2c translation of ACM TOMS Algorithm
# 733 cleaned up by S. G. Johnson), ported to R: LD_SLSQP.
#
# The C code passes pointers into one workspace and its subroutines (nnls,
# ldp, lsi, hfti, lsei, lsq, ldl, slsqpb) write through them, often with
# aliasing. The port keeps every array of the C code in one numeric vector
# W (the integer index arrays too, as exact doubles) and passes 1-based
# offsets into it: a C pointer p maps to the offset of p[0], so p[i] (after
# the f2c `--p`) is W[p + i - 1] and a column-major a[i + j * lda] is
# W[p + (i - 1) + (j - 1) * lda]. Each routine returns the updated W.
#
# The NLopt library that NLopt.jl loads (NLopt_jll) rounds every product
# and sum separately (no fused multiply add), so plain R arithmetic in the
# order of the C expressions reproduces it; the BLAS-like loops (ddot,
# dnrm2) accumulate in their loop order.
#
# The Fortran code of SLSQP is Copyright (c) 1988 Dieter Kraft and was
# released by the ACM under the BSD license; the C translation is Copyright
# (c) 2009 Massachusetts Institute of Technology, released under the MIT
# license (see R/optimization-nlopt-util.R).

.nsq_min2 <- function(a, b) if (a <= b) a else b
.nsq_max2 <- function(a, b) if (a >= b) a else b
.nsq_up <- function(a, b) if (a <= b) a:b else integer(0)
.nsq_down <- function(a, b) if (a >= b) a:b else integer(0)

# the BLAS-like helpers of slsqp.c
.nsq_dcopy <- function(W, n, px, incx, py, incy) {
  if (n <= 0) return(W)
  k <- 0:(n - 1)
  W[py + k * incy] <- if (incx == 0) W[px] else W[px + k * incx]
  W
}
.nsq_daxpy <- function(W, n, da, px, incx, py, incy) {
  if (n <= 0 || da == 0) return(W)
  k <- 0:(n - 1)
  iy <- py + k * incy
  W[iy] <- W[iy] + da * W[px + k * incx]
  W
}
.nsq_ddot <- function(W, n, px, incx, py, incy) {
  if (n <= 0) return(0)
  k <- 0:(n - 1)
  p <- W[px + k * incx] * W[py + k * incy]
  s <- 0
  for (v in p) s <- s + v
  s
}
.nsq_dnrm2 <- function(W, n, px, incx) {
  if (n <= 0) return(0)
  v <- W[px + (0:(n - 1)) * incx]
  xmax <- 0
  for (a in abs(v)) if (xmax < a) xmax <- a
  if (xmax == 0) return(0)
  scale <- 1 / xmax
  s <- 0
  for (a in scale * v) s <- s + a * a
  xmax * sqrt(s)
}
.nsq_dsrot <- function(W, n, px, incx, py, incy, c, s) {
  k <- 0:(n - 1)
  ix <- px + k * incx; iy <- py + k * incy
  x <- W[ix]; y <- W[iy]
  W[ix] <- c * x + s * y
  W[iy] <- c * y - s * x
  W
}
.nsq_dsrotg <- function(da, db) {
  absa <- abs(da); absb <- abs(db)
  if (absa > absb) { roe <- da; scale <- absa } else { roe <- db; scale <- absb }
  if (scale != 0) {
    iscale <- 1 / scale
    tmpa <- da * iscale; tmpb <- db * iscale
    r <- (if (roe < 0) -scale else scale) * sqrt(tmpa * tmpa + tmpb * tmpb)
    c <- da / r; s <- db / r
    da <- r
    db <- if (c != 0 && abs(c) <= s) 1 / c else s
  } else {
    c <- 1; s <- 0; da <- 0; db <- 0
  }
  list(da = da, db = db, c = c, s = s)
}
.nsq_dscal <- function(W, n, da, px, incx) {
  if (n <= 0) return(W)
  ix <- px + (0:(n - 1)) * incx
  W[ix] <- W[ix] * da
  W
}

# h12: construction (mode 1) and application (mode 2) of a Householder
# transformation; u (stride iue) holds the pivot vector, W[pup] is up
.nsq_h12 <- function(W, mode, lpivot, l1, m, pu, iue, pup, pc, ice, icv, ncv) {
  if (0 >= lpivot || lpivot >= l1 || l1 > m) return(W)
  U <- function(j) pu + (j - 1) * iue
  cl <- abs(W[U(lpivot)])
  if (mode != 2) {
    for (j in .nsq_up(l1, m)) {
      sm <- abs(W[U(j)])
      cl <- .nsq_max2(sm, cl)
    }
    if (cl <= 0) return(W)
    clinv <- 1 / cl
    d <- W[U(lpivot)] * clinv
    sm <- d * d
    for (j in .nsq_up(l1, m)) {
      d <- W[U(j)] * clinv
      sm <- sm + d * d
    }
    cl <- cl * sqrt(sm)
    if (W[U(lpivot)] > 0) cl <- -cl
    W[pup] <- W[U(lpivot)] - cl
    W[U(lpivot)] <- cl
  } else if (cl <= 0) {
    return(W)
  }
  if (ncv <= 0) return(W)
  b <- W[pup] * W[U(lpivot)]
  if (b >= 0) return(W)
  b <- 1 / b
  i2 <- 1 - icv + ice * (lpivot - 1)
  incr <- ice * (l1 - lpivot)
  for (j in .nsq_up(1, ncv)) {
    i2 <- i2 + icv
    i3 <- i2 + incr
    i4 <- i3
    sm <- W[pc + i2 - 1] * W[pup]
    for (i in .nsq_up(l1, m)) {
      sm <- sm + W[pc + i3 - 1] * W[U(i)]
      i3 <- i3 + ice
    }
    if (sm == 0) next
    sm <- sm * b
    W[pc + i2 - 1] <- W[pc + i2 - 1] + sm * W[pup]
    for (i in .nsq_up(l1, m)) {
      W[pc + i4 - 1] <- W[pc + i4 - 1] + sm * W[U(i)]
      i4 <- i4 + ice
    }
  }
  W
}

# nnls: nonnegative least squares of Lawson and Hanson; pup is a scratch
# slot for the Householder up
.nsq_nnls <- function(W, pa, mda, m, n, pb, px, pw, pz, pind, pup) {
  A <- function(i, j) pa + (i - 1) + (j - 1) * mda
  B <- function(i) pb + i - 1
  X <- function(i) px + i - 1
  Wd <- function(i) pw + i - 1
  Z <- function(i) pz + i - 1
  IND <- function(i) pind + i - 1
  factor <- 0.01
  mode <- 2; rnorm <- 0
  if (m <= 0 || n <= 0) return(list(W = W, rnorm = rnorm, mode = mode))
  mode <- 1
  iter <- 0; itmax <- n * 3
  for (i in .nsq_up(1, n)) W[IND(i)] <- i
  iz1 <- 1; iz2 <- n; nsetp <- 0; npp1 <- 1
  W[X(1)] <- 0
  W <- .nsq_dcopy(W, n, X(1), 0, X(1), 1)
  jj <- 0; i <- 0; j <- 0; iz <- 0; izmax <- 0
  lab <- "L110"
  repeat switch(lab,
    L110 = {
      if (iz1 > iz2 || nsetp >= m) { lab <- "L280"; next }
      for (iz in .nsq_up(iz1, iz2)) {
        j <- W[IND(iz)]
        W[Wd(j)] <- .nsq_ddot(W, m - nsetp, A(npp1, j), 1, B(npp1), 1)
      }
      lab <- "L130"
    },
    L130 = {
      wmax <- 0
      for (iz in .nsq_up(iz1, iz2)) {
        j <- W[IND(iz)]
        if (W[Wd(j)] <= wmax) next
        wmax <- W[Wd(j)]
        izmax <- iz
      }
      if (wmax <= 0) { lab <- "L280"; next }
      iz <- izmax
      j <- W[IND(iz)]
      asave <- W[A(npp1, j)]
      W <- .nsq_h12(W, 1, npp1, npp1 + 1, m, A(1, j), 1, pup, Z(1), 1, 1, 0)
      unorm <- .nsq_dnrm2(W, nsetp, A(1, j), 1)
      t <- factor * abs(W[A(npp1, j)])
      d <- unorm + t
      if (d - unorm > 0) {
        W <- .nsq_dcopy(W, m, B(1), 1, Z(1), 1)
        W <- .nsq_h12(W, 2, npp1, npp1 + 1, m, A(1, j), 1, pup, Z(1), 1, 1, 1)
        if (W[Z(npp1)] / W[A(npp1, j)] > 0) { lab <- "L160"; next }
      }
      W[A(npp1, j)] <- asave
      W[Wd(j)] <- 0
      lab <- "L130"
    },
    L160 = {
      W <- .nsq_dcopy(W, m, Z(1), 1, B(1), 1)
      W[IND(iz)] <- W[IND(iz1)]
      W[IND(iz1)] <- j
      iz1 <- iz1 + 1
      nsetp <- npp1
      npp1 <- npp1 + 1
      for (jz in .nsq_up(iz1, iz2)) {
        jj <- W[IND(jz)]
        W <- .nsq_h12(W, 2, nsetp, npp1, m, A(1, j), 1, pup, A(1, jj), 1, mda, 1)
      }
      k <- .nsq_min2(npp1, mda)
      W[Wd(j)] <- 0
      W <- .nsq_dcopy(W, m - nsetp, Wd(j), 0, A(k, j), 1)
      lab <- "L180"
    },
    L180 = {
      for (ip in .nsq_down(nsetp, 1)) {
        if (ip != nsetp) W <- .nsq_daxpy(W, ip, -W[Z(ip + 1)], A(1, jj), 1, Z(1), 1)
        jj <- W[IND(ip)]
        W[Z(ip)] <- W[Z(ip)] / W[A(ip, jj)]
      }
      iter <- iter + 1
      if (iter > itmax) { mode <- 3; lab <- "L280"; next }
      alpha <- 1
      jj <- 0
      for (ip in .nsq_up(1, nsetp)) {
        if (W[Z(ip)] > 0) next
        l <- W[IND(ip)]
        t <- -W[X(l)] / (W[Z(ip)] - W[X(l)])
        if (alpha < t) next
        alpha <- t
        jj <- ip
      }
      for (ip in .nsq_up(1, nsetp)) {
        l <- W[IND(ip)]
        W[X(l)] <- (1 - alpha) * W[X(l)] + alpha * W[Z(ip)]
      }
      if (jj == 0) { lab <- "L110"; next }
      i <- W[IND(jj)]
      lab <- "L250"
    },
    L250 = {
      W[X(i)] <- 0
      jj <- jj + 1
      for (j in .nsq_up(jj, nsetp)) {
        ii <- W[IND(j)]
        W[IND(j - 1)] <- ii
        rot <- .nsq_dsrotg(W[A(j - 1, ii)], W[A(j, ii)])
        W[A(j - 1, ii)] <- rot$da; W[A(j, ii)] <- rot$db
        t <- W[A(j - 1, ii)]
        W <- .nsq_dsrot(W, n, A(j - 1, 1), mda, A(j, 1), mda, rot$c, rot$s)
        W[A(j - 1, ii)] <- t
        W[A(j, ii)] <- 0
        W <- .nsq_dsrot(W, 1, B(j - 1), 1, B(j), 1, rot$c, rot$s)
      }
      npp1 <- nsetp
      nsetp <- nsetp - 1
      iz1 <- iz1 - 1
      W[IND(iz1)] <- i
      if (nsetp <= 0) { mode <- 3; lab <- "L280"; next }
      again <- FALSE
      for (jj in .nsq_up(1, nsetp)) {
        i <- W[IND(jj)]
        if (W[X(i)] <= 0) { again <- TRUE; break }
      }
      if (again) { lab <- "L250"; next }
      W <- .nsq_dcopy(W, m, B(1), 1, Z(1), 1)
      lab <- "L180"
    },
    L280 = {
      k <- .nsq_min2(npp1, m)
      rnorm <- .nsq_dnrm2(W, m - nsetp, B(k), 1)
      if (npp1 > m) {
        W[Wd(1)] <- 0
        W <- .nsq_dcopy(W, n, Wd(1), 0, Wd(1), 1)
      }
      break
    })
  list(W = W, rnorm = rnorm, mode = mode)
}

# ldp: least distance programming
.nsq_ldp <- function(W, pg, mg, m, n, ph, px, xnorm, pw, pind, pup) {
  G <- function(i, j) pg + (i - 1) + (j - 1) * mg
  H <- function(i) ph + i - 1
  X <- function(i) px + i - 1
  Wd <- function(i) pw + i - 1
  out <- function(mode) list(W = W, xnorm = xnorm, mode = mode)
  if (n <= 0) return(out(2))
  W[X(1)] <- 0
  W <- .nsq_dcopy(W, n, X(1), 0, X(1), 1)
  xnorm <- 0
  if (m == 0) return(out(1))
  iw <- 0
  for (j in .nsq_up(1, m)) {
    for (i in .nsq_up(1, n)) {
      iw <- iw + 1
      W[Wd(iw)] <- W[G(j, i)]
    }
    iw <- iw + 1
    W[Wd(iw)] <- W[H(j)]
  }
  if_ <- iw + 1
  for (i in .nsq_up(1, n)) {
    iw <- iw + 1
    W[Wd(iw)] <- 0
  }
  W[Wd(iw + 1)] <- 1
  n1 <- n + 1
  iz <- iw + 2
  iy <- iz + n1
  iwdual <- iy + m
  r <- .nsq_nnls(W, Wd(1), n1, n1, m, Wd(if_), Wd(iy), Wd(iwdual), Wd(iz), pind, pup)
  W <- r$W
  if (r$mode != 1) return(out(r$mode))
  if (r$rnorm <= 0) return(out(4))
  fac <- 1 - .nsq_ddot(W, m, H(1), 1, Wd(iy), 1)
  d <- 1 + fac
  if (d - 1 <= 0) return(out(4))
  fac <- 1 / fac
  for (j in .nsq_up(1, n)) W[X(j)] <- fac * .nsq_ddot(W, m, G(1, j), 1, Wd(iy), 1)
  xnorm <- .nsq_dnrm2(W, n, X(1), 1)
  W[Wd(1)] <- 0
  W <- .nsq_dcopy(W, m, Wd(1), 0, Wd(1), 1)
  W <- .nsq_daxpy(W, m, fac, Wd(iy), 1, Wd(1), 1)
  out(1)
}

# lsi: least squares with linear inequality constraints; pt is a scratch
# slot for the Householder up
.nsq_lsi <- function(W, pe, pf, pg, ph, le, me, lg, mg, n, px, xnorm, pw, pjw, pt, pup) {
  epmach <- 2.22e-16
  E <- function(i, j) pe + (i - 1) + (j - 1) * le
  Fv <- function(i) pf + i - 1
  G <- function(i, j) pg + (i - 1) + (j - 1) * lg
  H <- function(i) ph + i - 1
  X <- function(i) px + i - 1
  for (i in .nsq_up(1, n)) {
    j <- .nsq_min2(i + 1, n)
    W <- .nsq_h12(W, 1, i, i + 1, me, E(1, i), 1, pt, E(1, j), 1, le, n - i)
    W <- .nsq_h12(W, 2, i, i + 1, me, E(1, i), 1, pt, Fv(1), 1, 1, 1)
  }
  for (i in .nsq_up(1, mg)) {
    for (j in .nsq_up(1, n)) {
      if (abs(W[E(j, j)]) < epmach) return(list(W = W, xnorm = xnorm, mode = 5))
      W[G(i, j)] <- (W[G(i, j)] - .nsq_ddot(W, j - 1, G(i, 1), lg, E(1, j), 1)) / W[E(j, j)]
    }
    W[H(i)] <- W[H(i)] - .nsq_ddot(W, n, G(i, 1), lg, Fv(1), 1)
  }
  r <- .nsq_ldp(W, pg, lg, mg, n, ph, px, xnorm, pw, pjw, pup)
  W <- r$W; xnorm <- r$xnorm
  if (r$mode != 1) return(list(W = W, xnorm = xnorm, mode = r$mode))
  W <- .nsq_daxpy(W, n, 1, Fv(1), 1, X(1), 1)
  for (i in .nsq_down(n, 1)) {
    j <- .nsq_min2(i + 1, n)
    W[X(i)] <- (W[X(i)] - .nsq_ddot(W, n - i, E(i, j), le, X(j), 1)) / W[E(i, i)]
  }
  j <- .nsq_min2(n + 1, me)
  t <- .nsq_dnrm2(W, me - n, Fv(j), 1)
  xnorm <- sqrt(xnorm * xnorm + t * t)
  list(W = W, xnorm = xnorm, mode = 1)
}

# hfti: rank deficient least squares by Householder triangularization with
# column pivoting (nb right hand sides); returns rnorm[1] and the rank
.nsq_hfti <- function(W, pa, mda, m, n, pb, mdb, nb, tau, rnorm1, ph, pg, pip) {
  factor <- 0.001
  A <- function(i, j) pa + (i - 1) + (j - 1) * mda
  B <- function(i, j) pb + (i - 1) + (j - 1) * mdb
  H <- function(i) ph + i - 1
  G <- function(i) pg + i - 1
  IP <- function(i) pip + i - 1
  k <- 0
  ldiag <- .nsq_min2(m, n)
  if (ldiag <= 0) return(list(W = W, rnorm1 = rnorm1, krank = 0))
  hmax <- 0
  for (j in .nsq_up(1, ldiag)) {
    recompute <- TRUE
    if (j != 1) {
      lmax <- j
      for (l in .nsq_up(j, n)) {
        d <- W[A(j - 1, l)]
        W[H(l)] <- W[H(l)] - d * d
        if (W[H(l)] > W[H(lmax)]) lmax <- l
      }
      d <- hmax + factor * W[H(lmax)]
      if (d - hmax > 0) recompute <- FALSE
    }
    if (recompute) {
      lmax <- j
      for (l in .nsq_up(j, n)) {
        W[H(l)] <- 0
        for (i in .nsq_up(j, m)) {
          d <- W[A(i, l)]
          W[H(l)] <- W[H(l)] + d * d
        }
        if (W[H(l)] > W[H(lmax)]) lmax <- l
      }
      hmax <- W[H(lmax)]
    }
    W[IP(j)] <- lmax
    if (lmax != j) {
      for (i in .nsq_up(1, m)) {
        tmp <- W[A(i, j)]
        W[A(i, j)] <- W[A(i, lmax)]
        W[A(i, lmax)] <- tmp
      }
      W[H(lmax)] <- W[H(j)]
    }
    i <- .nsq_min2(j + 1, n)
    W <- .nsq_h12(W, 1, j, j + 1, m, A(1, j), 1, H(j), A(1, i), 1, mda, n - j)
    W <- .nsq_h12(W, 2, j, j + 1, m, A(1, j), 1, H(j), B(1, 1), 1, mdb, nb)
  }
  k <- ldiag
  for (j in .nsq_up(1, ldiag)) if (abs(W[A(j, j)]) <= tau) { k <- j - 1; break }
  kp1 <- k + 1
  for (jb in .nsq_up(1, nb)) {
    v <- .nsq_dnrm2(W, m - k, B(kp1, jb), 1)
    if (jb == 1) rnorm1 <- v
  }
  if (k <= 0) {
    for (jb in .nsq_up(1, nb)) for (i in .nsq_up(1, n)) W[B(i, jb)] <- 0
    return(list(W = W, rnorm1 = rnorm1, krank = k))
  }
  if (k != n) {
    for (i in .nsq_down(k, 1)) W <- .nsq_h12(W, 1, i, kp1, n, A(i, 1), mda, G(i), A(1, 1), mda, 1, i - 1)
  }
  for (jb in .nsq_up(1, nb)) {
    for (i in .nsq_down(k, 1)) {
      j <- .nsq_min2(i + 1, n)
      W[B(i, jb)] <- (W[B(i, jb)] - .nsq_ddot(W, k - i, A(i, j), mda, B(j, jb), 1)) / W[A(i, i)]
    }
    if (k != n) {
      for (j in .nsq_up(kp1, n)) W[B(j, jb)] <- 0
      for (i in .nsq_up(1, k)) W <- .nsq_h12(W, 2, i, kp1, n, A(i, 1), mda, G(i), B(1, jb), 1, mdb, 1)
    }
    for (j in .nsq_down(ldiag, 1)) {
      l <- W[IP(j)]
      if (l == j) next
      tmp <- W[B(l, jb)]
      W[B(l, jb)] <- W[B(j, jb)]
      W[B(j, jb)] <- tmp
    }
  }
  list(W = W, rnorm1 = rnorm1, krank = k)
}

# lsei: least squares with linear equality and inequality constraints
.nsq_lsei <- function(W, pc, pd, pe, pf, pg, ph, lc, mc, le, me, lg, mg, n, px, xnrm, pw, pjw, pt, pup) {
  epmach <- 2.22e-16
  C <- function(i, j) pc + (i - 1) + (j - 1) * lc
  D <- function(i) pd + i - 1
  E <- function(i, j) pe + (i - 1) + (j - 1) * le
  Fv <- function(i) pf + i - 1
  G <- function(i, j) pg + (i - 1) + (j - 1) * lg
  H <- function(i) ph + i - 1
  X <- function(i) px + i - 1
  Wd <- function(i) pw + i - 1
  out <- function(mode) list(W = W, xnrm = xnrm, mode = mode)
  if (mc > n) return(out(2))
  l <- n - mc
  mc1 <- mc + 1
  iw <- (l + 1) * (mg + 2) + 2 * mg + mc
  ie <- iw + mc + 1
  if_ <- ie + me * l
  ig <- if_ + me
  for (i in .nsq_up(1, mc)) {
    j <- .nsq_min2(i + 1, lc)
    W <- .nsq_h12(W, 1, i, i + 1, n, C(i, 1), lc, Wd(iw + i), C(j, 1), lc, 1, mc - i)
    W <- .nsq_h12(W, 2, i, i + 1, n, C(i, 1), lc, Wd(iw + i), E(1, 1), le, 1, me)
    W <- .nsq_h12(W, 2, i, i + 1, n, C(i, 1), lc, Wd(iw + i), G(1, 1), lg, 1, mg)
  }
  for (i in .nsq_up(1, mc)) {
    if (abs(W[C(i, i)]) < epmach) return(out(6))
    W[X(i)] <- (W[D(i)] - .nsq_ddot(W, i - 1, C(i, 1), lc, X(1), 1)) / W[C(i, i)]
  }
  mode <- 1
  W[Wd(mc1)] <- 0
  W <- .nsq_dcopy(W, mg, Wd(mc1), 0, Wd(mc1), 1)
  if (mc != n) {
    for (i in .nsq_up(1, me)) W[Wd(if_ - 1 + i)] <- W[Fv(i)] - .nsq_ddot(W, mc, E(i, 1), le, X(1), 1)
    for (i in .nsq_up(1, me)) W <- .nsq_dcopy(W, l, E(i, mc1), le, Wd(ie - 1 + i), me)
    for (i in .nsq_up(1, mg)) W <- .nsq_dcopy(W, l, G(i, mc1), lg, Wd(ig - 1 + i), mg)
    if (mg > 0) {
      for (i in .nsq_up(1, mg)) W[H(i)] <- W[H(i)] - .nsq_ddot(W, mc, G(i, 1), lg, X(1), 1)
      r <- .nsq_lsi(W, Wd(ie), Wd(if_), Wd(ig), H(1), me, me, mg, mg, l, X(mc1), xnrm, Wd(mc1), pjw, pt, pup)
      W <- r$W; xnrm <- r$xnorm; mode <- r$mode
      if (mc == 0) return(out(mode))
      t <- .nsq_dnrm2(W, mc, X(1), 1)
      xnrm <- sqrt(xnrm * xnrm + t * t)
      if (mode != 1) return(out(mode))
    } else {
      k <- .nsq_max2(le, n)
      t <- sqrt(epmach)
      r <- .nsq_hfti(W, Wd(ie), me, me, l, Wd(if_), k, 1, t, xnrm, Wd(1), Wd(l + 1), pjw)
      W <- r$W; xnrm <- r$rnorm1
      W <- .nsq_dcopy(W, l, Wd(if_), 1, X(mc1), 1)
      if (r$krank != l) return(out(7))
      mode <- 1
    }
  }
  for (i in .nsq_up(1, me)) W[Fv(i)] <- .nsq_ddot(W, n, E(i, 1), le, X(1), 1) - W[Fv(i)]
  for (i in .nsq_up(1, mc)) {
    W[D(i)] <- .nsq_ddot(W, me, E(1, i), 1, Fv(1), 1) - .nsq_ddot(W, mg, G(1, i), 1, Wd(mc1), 1)
  }
  for (i in .nsq_down(mc, 1)) W <- .nsq_h12(W, 2, i, i + 1, n, C(i, 1), lc, Wd(iw + i), X(1), 1, 1, 1)
  for (i in .nsq_down(mc, 1)) {
    j <- .nsq_min2(i + 1, lc)
    W[Wd(i)] <- (W[D(i)] - .nsq_ddot(W, mc - i, C(j, i), 1, Wd(j), 1)) / W[C(i, i)]
  }
  out(mode)
}

# lsq: the quadratic programming subproblem of SLSQP as a least squares
# problem (l the packed LDL' factor, xl and xu the bounds of the step)
.nsq_lsq <- function(W, m, meq, n, nl, la, pl, pg, pa, pb, pxl, pxu, px, py, pw, pjw, pt, pup) {
  A <- function(i, j) pa + (i - 1) + (j - 1) * la
  Wd <- function(i) pw + i - 1
  Lv <- function(i) pl + i - 1
  n1 <- n + 1
  mineq <- m - meq
  m1 <- mineq + n + n
  n2 <- (n1 * n) %/% 2 + 1
  n2 <- if (n2 == nl) 0 else 1
  n3 <- n - n2
  i2 <- 1; i3 <- 1; i4 <- 1; ie <- 1
  if_ <- n * n + 1
  for (i in .nsq_up(1, n3)) {
    i1 <- n1 - i
    diag <- sqrt(W[Lv(i2)])
    W[Wd(i3)] <- 0
    W <- .nsq_dcopy(W, i1, Wd(i3), 0, Wd(i3), 1)
    W <- .nsq_dcopy(W, i1 - n2, Lv(i2), 1, Wd(i3), n)
    W <- .nsq_dscal(W, i1 - n2, diag, Wd(i3), n)
    W[Wd(i3)] <- diag
    W[Wd(if_ - 1 + i)] <- (W[pg + i - 1] - .nsq_ddot(W, i - 1, Wd(i4), 1, Wd(if_), 1)) / diag
    i2 <- i2 + i1 - n2
    i3 <- i3 + n1
    i4 <- i4 + n
  }
  if (n2 == 1) {
    W[Wd(i3)] <- W[Lv(nl)]
    W[Wd(i4)] <- 0
    W <- .nsq_dcopy(W, n3, Wd(i4), 0, Wd(i4), 1)
    W[Wd(if_ - 1 + n)] <- 0
  }
  W <- .nsq_dscal(W, n, -1, Wd(if_), 1)
  ic <- if_ + n
  id <- ic + meq * n
  if (meq > 0) {
    for (i in .nsq_up(1, meq)) W <- .nsq_dcopy(W, n, A(i, 1), la, Wd(ic - 1 + i), meq)
    W <- .nsq_dcopy(W, meq, pb, 1, Wd(id), 1)
    W <- .nsq_dscal(W, meq, -1, Wd(id), 1)
  }
  ig <- id + meq
  if (mineq > 0) for (i in .nsq_up(1, mineq)) W <- .nsq_dcopy(W, n, A(meq + i, 1), la, Wd(ig - 1 + i), m1)
  ip <- ig + mineq
  for (i in .nsq_up(1, n)) {
    W[Wd(ip - 1 + i)] <- 0
    W <- .nsq_dcopy(W, n, Wd(ip - 1 + i), 0, Wd(ip - 1 + i), m1)
  }
  for (i in .nsq_up(1, n)) if (!.nlopt_isinf(W[pxl + i - 1])) W[Wd((ip - (m1 + 1)) + i * (m1 + 1))] <- 1
  im <- ip + n
  for (i in .nsq_up(1, n)) {
    W[Wd(im - 1 + i)] <- 0
    W <- .nsq_dcopy(W, n, Wd(im - 1 + i), 0, Wd(im - 1 + i), m1)
  }
  for (i in .nsq_up(1, n)) if (!.nlopt_isinf(W[pxu + i - 1])) W[Wd((im - (m1 + 1)) + i * (m1 + 1))] <- -1
  ih <- ig + m1 * n
  if (mineq > 0) {
    W <- .nsq_dcopy(W, mineq, pb + meq, 1, Wd(ih), 1)
    W <- .nsq_dscal(W, mineq, -1, Wd(ih), 1)
  }
  il <- ih + mineq
  iu <- il + n
  for (i in .nsq_up(1, n)) {
    xl <- W[pxl + i - 1]; xu <- W[pxu + i - 1]
    W[Wd(il - 1 + i)] <- if (.nlopt_isinf(xl)) 0 else xl
    W[Wd(iu - 1 + i)] <- if (.nlopt_isinf(xu)) 0 else -xu
  }
  iw <- iu + n
  r <- .nsq_lsei(W, Wd(ic), Wd(id), Wd(ie), Wd(if_), Wd(ig), Wd(ih), .nsq_max2(1, meq), meq, n, n,
                 m1, m1, n, px, 0, Wd(iw), pjw, pt, pup)
  W <- r$W
  if (r$mode == 1) {
    W <- .nsq_dcopy(W, m, Wd(iw), 1, py, 1)
    W <- .nsq_dcopy(W, n3, Wd(iw + m), 1, py + m, 1)
    W <- .nsq_dcopy(W, n3, Wd(iw + m + n), 1, py + m + n3, 1)
    for (i in .nsq_up(1, n)) {
      xi <- W[px + i - 1]
      if (xi < W[pxl + i - 1]) W[px + i - 1] <- W[pxl + i - 1]
      else if (xi > W[pxu + i - 1]) W[px + i - 1] <- W[pxu + i - 1]
    }
  }
  list(W = W, mode = r$mode)
}

# ldl: the rank one update a + sigma z z' of the packed LDL' factor a
.nsq_ldl <- function(W, n, pa, pz, sigma, pw) {
  epmach <- 2.22e-16
  Av <- function(i) pa + i - 1
  Z <- function(i) pz + i - 1
  Wd <- function(i) pw + i - 1
  if (sigma == 0) return(W)
  ij <- 1
  t <- 1 / sigma
  if (!(sigma > 0)) {
    for (i in .nsq_up(1, n)) W[Wd(i)] <- W[Z(i)]
    for (i in .nsq_up(1, n)) {
      v <- W[Wd(i)]
      t <- t + v * v / W[Av(ij)]
      for (j in .nsq_up(i + 1, n)) {
        ij <- ij + 1
        W[Wd(j)] <- W[Wd(j)] - v * W[Av(ij)]
      }
      ij <- ij + 1
    }
    if (t >= 0) t <- epmach / sigma
    for (i in .nsq_up(1, n)) {
      j <- n + 1 - i
      ij <- ij - i
      u <- W[Wd(j)]
      W[Wd(j)] <- t
      t <- t - u * u / W[Av(ij)]
    }
  }
  for (i in .nsq_up(1, n)) {
    v <- W[Z(i)]
    delta <- v / W[Av(ij)]
    tp <- if (sigma < 0) W[Wd(i)] else t + delta * v
    alpha <- tp / t
    W[Av(ij)] <- alpha * W[Av(ij)]
    if (i == n) break
    beta <- delta / tp
    if (alpha > 4) {
      gamma <- t / tp
      for (j in .nsq_up(i + 1, n)) {
        ij <- ij + 1
        u <- W[Av(ij)]
        W[Av(ij)] <- gamma * u + beta * W[Z(j)]
        W[Z(j)] <- W[Z(j)] - v * u
      }
    } else {
      for (j in .nsq_up(i + 1, n)) {
        ij <- ij + 1
        W[Z(j)] <- W[Z(j)] - v * W[Av(ij)]
        W[Av(ij)] <- W[Av(ij)] + beta * W[Z(j)]
      }
    }
    ij <- ij + 1
    t <- tp
  }
  W
}

# slsqpb: one reverse communication step of SLSQP. `st` carries the saved
# state of the C code (slsqpb_state) plus iter and acc; p holds the offsets
# of the arrays in W.
.nsq_slsqpb <- function(W, st, m, meq, la, n, p, f, mode) {
  for (nm in names(st)) assign(nm, st[[nm]])
  X <- function(i) p$x + i - 1
  XL <- function(i) p$xl + i - 1
  XU <- function(i) p$xu + i - 1
  Cv <- function(i) p$c + i - 1
  G <- function(i) p$g + i - 1
  A <- function(i, j) p$a + (i - 1) + (j - 1) * la
  R <- function(i) p$r + i - 1
  L <- function(i) p$l + i - 1
  X0 <- function(i) p$x0 + i - 1
  MU <- function(i) p$mu + i - 1
  S <- function(i) p$s + i - 1
  U <- function(i) p$u + i - 1
  V <- function(i) p$v + i - 1
  ccon <- function(j) if (j <= meq) W[Cv(j)] else 0
  lsq <- function(nn) .nsq_lsq(W, m, meq, nn, n3, la, L(1), G(1), A(1, 1), Cv(1), U(1), V(1), S(1), R(1),
                               p$w, p$iw, p$t, p$up)
  lab <- if (mode == -1) "L260" else if (mode == 0) "L100" else "L220"
  repeat switch(lab,
    L100 = {
      itermx <- iter
      iexact <- if (acc >= 0) 0 else 1
      acc <- abs(acc)
      tol <- 10 * acc
      iter <- 0
      ireset <- 0
      n1 <- n + 1
      n2 <- (n1 * n) %/% 2
      n3 <- n2 + 1
      W[S(1)] <- 0
      W[MU(1)] <- 0
      W <- .nsq_dcopy(W, n, S(1), 0, S(1), 1)
      W <- .nsq_dcopy(W, m, MU(1), 0, MU(1), 1)
      lab <- "L110"
    },
    L110 = {
      ireset <- ireset + 1
      if (ireset > 5) { lab <- "L255"; next }
      W[L(1)] <- 0
      W <- .nsq_dcopy(W, n2, L(1), 0, L(1), 1)
      j <- 1
      for (i in .nsq_up(1, n)) {
        W[L(j)] <- 1
        j <- j + n1 - i
      }
      lab <- "L130"
    },
    L130 = {
      iter <- iter + 1
      mode <- 9
      if (iter > itermx && itermx > 0) { lab <- "L330"; next }
      W <- .nsq_dcopy(W, n, XL(1), 1, U(1), 1)
      W <- .nsq_dcopy(W, n, XU(1), 1, V(1), 1)
      W <- .nsq_daxpy(W, n, -1, X(1), 1, U(1), 1)
      W <- .nsq_daxpy(W, n, -1, X(1), 1, V(1), 1)
      h4 <- 1
      r <- lsq(n)
      W <- r$W; mode <- r$mode
      if (mode == 6 && n == meq) mode <- 4
      if (mode == 4) {
        for (j in .nsq_up(1, m)) W[A(j, n1)] <- if (j <= meq) -W[Cv(j)] else .nsq_max2(-W[Cv(j)], 0)
        W[S(1)] <- 0
        W <- .nsq_dcopy(W, n, S(1), 0, S(1), 1)
        h3 <- 0
        W[G(n1)] <- 0
        W[L(n3)] <- 100
        W[S(n1)] <- 1
        W[U(n1)] <- 0
        W[V(n1)] <- 1
        incons <- 0
        lab <- "L150"
      } else if (mode != 1) {
        lab <- "L330"
      } else {
        lab <- "L160"
      }
    },
    L150 = {
      r <- lsq(n1)
      W <- r$W; mode <- r$mode
      h4 <- 1 - W[S(n1)]
      if (mode == 4) {
        W[L(n3)] <- 10 * W[L(n3)]
        incons <- incons + 1
        lab <- if (incons > 5) "L330" else "L150"
      } else if (mode != 1) {
        lab <- "L330"
      } else {
        lab <- "L160"
      }
    },
    L160 = {
      for (i in .nsq_up(1, n)) W[V(i)] <- W[G(i)] - .nsq_ddot(W, m, A(1, i), 1, R(1), 1)
      f0 <- f
      W <- .nsq_dcopy(W, n, X(1), 1, X0(1), 1)
      gs <- .nsq_ddot(W, n, G(1), 1, S(1), 1)
      h1 <- abs(gs)
      h2 <- 0
      for (j in .nsq_up(1, m)) {
        h3 <- ccon(j)
        h2 <- h2 + .nsq_max2(-W[Cv(j)], h3)
        h3 <- abs(W[R(j)])
        W[MU(j)] <- .nsq_max2(h3, (W[MU(j)] + h3) / 2)
        h1 <- h1 + h3 * abs(W[Cv(j)])
      }
      mode <- 0
      if (h1 < acc && h2 < acc) { lab <- "L330"; next }
      h1 <- 0
      for (j in .nsq_up(1, m)) {
        h3 <- ccon(j)
        h1 <- h1 + W[MU(j)] * .nsq_max2(-W[Cv(j)], h3)
      }
      t0 <- f + h1
      h3 <- gs - h1 * h4
      mode <- 8
      if (h3 >= 0) { lab <- "L110"; next }
      line <- 0
      alpha <- 1
      lab <- if (iexact == 1) "L210" else "L190"
    },
    L190 = {
      line <- line + 1
      h3 <- alpha * h3
      W <- .nsq_dscal(W, n, alpha, S(1), 1)
      W <- .nsq_dcopy(W, n, X0(1), 1, X(1), 1)
      W <- .nsq_daxpy(W, n, 1, S(1), 1, X(1), 1)
      for (i in .nsq_up(1, n)) {
        if (W[X(i)] < W[XL(i)]) W[X(i)] <- W[XL(i)]
        else if (W[X(i)] > W[XU(i)]) W[X(i)] <- W[XU(i)]
      }
      mode <- if (line == 1) -2 else 1
      lab <- "L330"
    },
    L200 = {
      if (is.finite(h1)) {
        if (h1 <= h3 / 10 || line > 10) { lab <- "L240"; next }
        alpha <- .nsq_max2(h3 / (2 * (h3 - h1)), 0.1)
      } else {
        alpha <- .nsq_max2(alpha * 0.5, 0.1)
      }
      lab <- "L190"
    },
    L210 = {
      # the exact line search is compiled out of NLopt
      mode <- 9
      lab <- "L330"
    },
    L220 = {
      t <- f
      for (j in .nsq_up(1, m)) {
        h1 <- ccon(j)
        t <- t + W[MU(j)] * .nsq_max2(-W[Cv(j)], h1)
      }
      h1 <- t - t0
      lab <- if (iexact == 0) "L200" else "L210"
    },
    L240 = {
      h3 <- 0
      for (j in .nsq_up(1, m)) {
        h1 <- ccon(j)
        h3 <- h3 + .nsq_max2(-W[Cv(j)], h1)
      }
      mode <- if ((abs(f - f0) < acc || .nsq_dnrm2(W, n, S(1), 1) < acc) && h3 < acc) 0 else -1
      lab <- "L330"
    },
    L255 = {
      mode <- if ((abs(f - f0) < tol || .nsq_dnrm2(W, n, S(1), 1) < tol) && h3 < tol) 0 else 8
      lab <- "L330"
    },
    L260 = {
      for (i in .nsq_up(1, n)) W[U(i)] <- W[G(i)] - .nsq_ddot(W, m, A(1, i), 1, R(1), 1) - W[V(i)]
      k <- 0
      for (i in .nsq_up(1, n)) {
        h1 <- 0
        k <- k + 1
        for (j in .nsq_up(i + 1, n)) {
          k <- k + 1
          h1 <- h1 + W[L(k)] * W[S(j)]
        }
        W[V(i)] <- W[S(i)] + h1
      }
      k <- 1
      for (i in .nsq_up(1, n)) {
        W[V(i)] <- W[L(k)] * W[V(i)]
        k <- k + n1 - i
      }
      for (i in .nsq_down(n, 1)) {
        h1 <- 0
        k <- i
        for (j in .nsq_up(1, i - 1)) {
          h1 <- h1 + W[L(k)] * W[V(j)]
          k <- k + n - j
        }
        W[V(i)] <- W[V(i)] + h1
      }
      h1 <- .nsq_ddot(W, n, S(1), 1, U(1), 1)
      h2 <- .nsq_ddot(W, n, S(1), 1, V(1), 1)
      h3 <- h2 * 0.2
      if (h1 < h3) {
        h4 <- (h2 - h3) / (h2 - h1)
        h1 <- h3
        W <- .nsq_dscal(W, n, h4, U(1), 1)
        W <- .nsq_daxpy(W, n, 1 - h4, V(1), 1, U(1), 1)
      }
      W <- .nsq_ldl(W, n, L(1), U(1), 1 / h1, V(1))
      W <- .nsq_ldl(W, n, L(1), V(1), -1 / h2, U(1))
      lab <- "L130"
    },
    L330 = break)
  for (nm in names(st)) st[[nm]] <- get(nm)
  list(W = W, st = st, mode = mode)
}

# nlopt_slsqp: the driver of slsqp.c
.nlopt_slsqp <- function(o, x, stop) {
  n <- length(x)
  fc <- o$fc; h <- o$h
  mtot <- .nlopt_count_constraints(fc)
  ptot <- .nlopt_count_constraints(h)
  if (ptot > n) {
    stop$msg <- "slsqp: more equality constraints than variables"
    return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  }
  mpi <- mtot + ptot; meq <- ptot; mpi1 <- if (mpi > 0) mpi else 1
  n1 <- n + 1; mineq <- mpi - meq + n1 + n1
  len_w <- (3 * n1 + mpi) * (n1 + 1) + (n1 - meq + 1) * (mineq + 2) + 2 * mineq +
    (n1 + mineq) * (n1 - meq) + 2 * meq + n1 + ((n + 1) * n) %/% 2 + 2 * mpi + 3 * n + 3 * n1 + 1
  len_jw <- mineq
  # the layout of W: xcur, lb, ub, c, grad, cgrad, w, jw, two scratch slots
  p <- list(x = 1)
  p$xl <- p$x + n; p$xu <- p$xl + n; p$c <- p$xu + n; p$g <- p$c + mpi1
  p$a <- p$g + n1; p$wk <- p$a + mpi1 * n1; p$iw <- p$wk + len_w
  p$t <- p$iw + len_jw; p$up <- p$t + 1
  W <- numeric(p$up)
  W[p$x + 0:(n - 1)] <- x
  W[p$xl + 0:(n - 1)] <- o$lb
  W[p$xu + 0:(n - 1)] <- o$ub
  # the workspace offsets of slsqp()
  la <- mpi1
  im <- 1; il <- im + la; ix <- il + (n1 * n) %/% 2 + 1; ir <- ix + n
  is <- ir + n + n + la; iu <- is + n1; iv <- iu + n1; iw <- iv + n1
  wk <- function(i) p$wk + i - 1
  p$r <- wk(ir); p$l <- wk(il); p$x0 <- wk(ix); p$mu <- wk(im); p$s <- wk(is)
  p$u <- wk(iu); p$v <- wk(iv); p$w <- wk(iw)
  st <- list(t = 0, f0 = 0, h1 = 0, h2 = 0, h3 = 0, h4 = 0, n1 = 0, n2 = 0, n3 = 0, t0 = 0, gs = 0,
             tol = 0, line = 0, alpha = 0, iexact = 0, incons = 0, ireset = 0, itermx = 0,
             iter = 0, acc = 0)
  mode <- 0; prev_mode <- 0
  ret <- "SUCCESS"
  feasible <- FALSE; feasible_cur <- FALSE
  infeasibility <- Inf; infeasibility_cur <- Inf
  want_grad <- TRUE
  xcur <- x; xprev <- x; xbest <- x
  fprev <- Inf; fcur <- Inf; minf <- Inf
  ix_x <- p$x + 0:(n - 1)
  first <- TRUE
  repeat {
    do_eval <- FALSE
    if (first) {
      first <- FALSE
      want_grad <- TRUE; do_eval <- TRUE
    } else {
      r <- .nsq_slsqpb(W, st, mpi, meq, la, n, p, fcur, mode)
      W <- r$W; st <- r$st; mode <- r$mode
      xcur <- W[ix_x]
      if (mode == -1) {
        if (!(prev_mode == -2 && !want_grad)) { want_grad <- TRUE; do_eval <- TRUE }
      } else if (mode == -2) {
        want_grad <- TRUE; do_eval <- TRUE
      } else if (mode == 1) {
        do_eval <- TRUE
      } else if (mode == 0) {
        break
      } else if (mode == 8) {
        ret <- "ROUNDOFF_LIMITED"
        if (feasible_cur) {
          save <- c(stop$ftol_rel, stop$xtol_rel, stop$ftol_abs)
          stop$ftol_rel <- stop$ftol_rel * 10
          stop$ftol_abs <- stop$ftol_abs * 10
          stop$xtol_rel <- stop$xtol_rel * 10
          if (.nlopt_stop_ftol(stop, fcur, st$f0)) ret <- "FTOL_REACHED"
          else if (.nlopt_stop_x(stop, xcur, W[p$x0 + 0:(n - 1)])) ret <- "XTOL_REACHED"
          stop$ftol_rel <- save[1]; stop$xtol_rel <- save[2]; stop$ftol_abs <- save[3]
        }
        break
      } else if (mode %in% c(5, 6, 7)) {
        ret <- "ROUNDOFF_LIMITED"
        break
      } else if (mode %in% c(3, 4, 9)) {
        stop$msg <- "bug: more than iter SQP iterations"
        ret <- "FAILURE"
        break
      } else {
        stop$msg <- "bug: workspace is too small"
        ret <- "INVALID_ARGS"
        break
      }
    }
    if (do_eval) {
      newcgrad <- want_grad
      feasible_cur <- TRUE; infeasibility_cur <- 0
      if (want_grad) {
        fr <- o$fg(xcur)
        fcur <- fr$val
        W[p$g + 0:(n - 1)] <- fr$grad
      } else {
        fcur <- o$f(xcur)
      }
      stop$nevals <- stop$nevals + 1L
      if (.nlopt_stop_forced(stop)) { fcur <- Inf; ret <- "FORCED_STOP"; break }
      if (is.finite(fcur)) {
        want_grad <- FALSE
        ii <- 0
        forced <- FALSE
        for (pass in 1:2) {
          cons <- if (pass == 1) h else fc
          for (cc in cons) {
            if (newcgrad) {
              cr <- cc$fg(xcur)
              val <- as.numeric(cr$val)
              grad <- matrix(cr$grad, cc$m, n)
            } else {
              val <- as.numeric(cc$f(xcur))
            }
            if (.nlopt_stop_forced(stop)) { forced <- TRUE; break }
            tol <- rep_len(cc$tol, cc$m)
            for (k in seq_len(cc$m)) {
              cv <- val[k]
              if (pass == 1) {
                infeasibility_cur <- .nsq_max2(infeasibility_cur, abs(cv))
                feasible_cur <- feasible_cur && abs(cv) <= tol[k]
                if (newcgrad) W[p$a + ii + (0:(n - 1)) * mpi1] <- grad[k, ]
                W[p$c + ii] <- cv
              } else {
                infeasibility_cur <- .nsq_max2(infeasibility_cur, cv)
                feasible_cur <- feasible_cur && cv <= tol[k]
                if (newcgrad) W[p$a + ii + (0:(n - 1)) * mpi1] <- -grad[k, ]
                W[p$c + ii] <- -cv
              }
              ii <- ii + 1
            }
          }
          if (forced) break
        }
        if (forced) { ret <- "FORCED_STOP"; break }
      }
    }
    prev_mode <- mode
    if (is.finite(fcur) && ((fcur < minf && (feasible_cur || !feasible)) ||
                            (!feasible && infeasibility_cur < infeasibility))) {
      minf <- fcur
      feasible <- feasible_cur
      infeasibility <- infeasibility_cur
      xbest <- xcur
    }
    if (mode == -1) {
      if (!.nlopt_isinf(fprev) && feasible_cur) {
        if (.nlopt_stop_ftol(stop, fcur, fprev)) ret <- "FTOL_REACHED"
        else if (.nlopt_stop_x(stop, xcur, xprev)) ret <- "XTOL_REACHED"
      }
      fprev <- fcur
      xprev <- xcur
    }
    if (.nlopt_stop_evals(stop)) ret <- "MAXEVAL_REACHED"
    else if (.nlopt_stop_time(stop)) ret <- "MAXTIME_REACHED"
    else if (feasible && minf < stop$minf_max) ret <- "STOPVAL_REACHED"
    if (ret != "SUCCESS") break
  }
  if (.nlopt_isinf(minf)) {
    if (.nlopt_isinf(fcur)) { minf <- fprev; xbest <- xprev } else { minf <- fcur; xbest <- xcur }
  }
  list(ret = ret, x = xbest, minf = minf)
}
