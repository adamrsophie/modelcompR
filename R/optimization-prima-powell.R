# Kernels shared by the UOBYQA, NEWUOA, BOBYQA and LINCOA ports of PRIMA.
#
# No Fortran source of PRIMA is at hand. These routines are reconstructed
# from the compiled library that PRIMA.jl calls (libprimaf of PRIMA_jll,
# PRIMA 0.7.1, built by gfortran with fused multiply add contraction) and
# checked routine by routine against that library, calling its module
# procedures directly on random inputs, until the results agree in every
# bit. gfortran contracts a*b + c into fma(a, b, c); `.pw_fma()` (the exact
# emulation `.fma()` of R/solvers-abm.R) reproduces each contraction in the order of
# the compiled code. Where the compiled code runs a two lane vector loop
# and a scalar remainder that contract differently, `.pw_isvec()` tells the
# lanes apart. Indices are 1 based as in Fortran.

.pw_twopi <- 6.283185307179586

# z rounded to odd when its exact value lies on the side sign(ze) of z: an
# inexact z with an even last bit moves one unit in the last place toward
# the exact value
.pw_round_odd <- function(z, ze) {
  k <- which(ze != 0)
  if (length(k)) {
    zk <- z[k]
    az <- abs(zk)
    ex <- floor(log2(az))
    ex <- ex - (2^ex > az)
    ex <- ex + (2^(ex + 1) <= az)
    ulp <- 2^(pmax(ex, -1022) - 52)
    even <- (az / ulp) %% 2 == 0
    if (any(even)) {
      dir <- sign(ze[k]) * sign(zk)
      step <- ifelse(dir < 0 & az == 2^ex & ex > -1022, ulp / 2, ulp)
      z[k[even]] <- (sign(zk) * (az + dir * step))[even]
    }
  }
  z
}

# fma(a, b, c) correctly rounded, the operation that gfortran emits.
# `.fma()` (R/solvers-abm.R) rounds twice when the low parts of the exact
# sum fall on a rounding tie (as for squares of 1 + 2^-52); this version
# (Boldo and Melquiond, Emulation of FMA and correctly rounded sums: proved
# algorithms using rounding to odd, IEEE Trans. Comput. 57, 2008) sums the
# low parts rounded to odd, which leaves one rounding. A factor above 2^995,
# whose splitting would overflow, trades 2^60 with the other factor, which
# leaves the product unchanged. A product above 2^1000, which may overflow
# although the result does not, is formed scaled by 2^-64 (with c). A
# product below 2^-960 (or one that underflows to zero), whose low parts
# would underflow, is formed scaled by 2^120 (with c), and the result is
# scaled back; a subnormal result is first rounded to odd, so that the
# scaling rounds only once. Where the
# splitting still overflows (a non finite operand or product) the plain
# a * b + c gives the IEEE result, except that an infinite c with finite
# factors is returned as it is.
.pw_fma <- function(a, b, c) {
  a0 <- a
  b0 <- b
  uh <- a * b
  big <- which(abs(a) > 0x1p995 | abs(b) > 0x1p995)
  small <- which(abs(uh) < 0x1p-960 & abs(c) < 0x1p900 & (uh != 0 | (a != 0 & b != 0)))
  large <- which(!(abs(uh) <= 0x1p1000))
  back <- NULL
  down <- NULL
  if (length(big) || length(small) || length(large)) {
    m <- max(length(a), length(b), length(c))
    a <- rep_len(a, m)
    b <- rep_len(b, m)
    c <- rep_len(c, m)
    uh <- rep_len(uh, m)
    big <- which(abs(a) > 0x1p995 | abs(b) > 0x1p995)
    up <- abs(a[big]) > abs(b[big])
    up[is.na(up)] <- TRUE
    a[big] <- a[big] * ifelse(up, 0x1p-60, 0x1p60)
    b[big] <- b[big] * ifelse(up, 0x1p60, 0x1p-60)
    back <- which(abs(uh) < 0x1p-960 & abs(c) < 0x1p900 & (uh != 0 | (a != 0 & b != 0)))
    lo <- abs(a[back]) < abs(b[back])
    a[back] <- a[back] * ifelse(lo, 0x1p120, 1)
    b[back] <- b[back] * ifelse(lo, 1, 0x1p120)
    c[back] <- c[back] * 0x1p120
    uh <- a * b
    down <- which(!(abs(uh) <= 0x1p1000) & is.finite(a) & is.finite(b) & is.finite(c))
    a[down] <- a[down] * 0x1p-64
    cd <- c[down]
    c[down] <- cd * 0x1p-64
    # a tiny c that the scaling flushes to zero still breaks a rounding tie
    lost <- c[down] == 0 & cd != 0
    c[down[lost]] <- sign(cd[lost]) * 2^-1074
    uh <- a * b
  }
  ta <- a * 134217729
  ah <- ta - (ta - a)
  al <- a - ah
  tb <- b * 134217729
  bh <- tb - (tb - b)
  bl <- b - bh
  ul <- ((ah * bh - uh) + ah * bl + al * bh) + al * bl
  th <- c + ul
  bv <- th - c
  tl <- (c - (th - bv)) + (ul - bv)
  vh <- uh + th
  bv <- vh - uh
  vl <- (uh - (vh - bv)) + (th - bv)
  z <- vl + tl
  bv <- z - vl
  ze <- (vl - (z - bv)) + (tl - bv)
  zn <- z
  z <- .pw_round_odd(z, ze)
  r <- vh + z
  if (anyNA(r)) {
    bad <- is.na(r)
    r[bad] <- (a * b + c)[bad]
    # an infinite c with finite factors is the result, even when the product
    # (or a traded factor) overflows
    inf <- bad & is.infinite(c) & is.finite(a0) & is.finite(b0)
    r[inf] <- rep_len(c, length(r))[inf]
  }
  if (length(back)) {
    # a subnormal result: scaling back rounds once more, so r (rounded to
    # nearest) gives way to r rounded to odd, which rounds correctly to
    # 51 bits or fewer; on the 52 bit grid of [2^-1023, 2^-1022) an odd r
    # is a midpoint, decided by the side of the exact sum
    sb <- back[is.finite(r[back]) & abs(r[back]) < 0x1p-902]
    if (length(sb)) {
      bv <- r[sb] - vh[sb]
      e <- (vh[sb] - (r[sb] - bv)) + (z[sb] - bv)
      side <- ifelse(e != 0, sign(e), sign(ze[sb] - (z[sb] - zn[sb])))
      top <- abs(r[sb]) >= 0x1p-903
      ulp <- 0x1p-955
      odd <- (abs(r[sb]) / ulp) %% 2 == 1
      r[sb] <- ifelse(top, ifelse(odd & side != 0, r[sb] + side * ulp, r[sb]), .pw_round_odd(r[sb], side))
    }
    r[back] <- r[back] * 0x1p-120
  }
  if (length(down)) r[down] <- r[down] * 0x1p64
  # a zero result with a zero product takes the sign of a * b + c, and with
  # an underflowing product (and a zero c) the sign of the product
  zr <- which(r == 0)
  if (length(zr)) {
    p <- rep_len(a0 * b0, length(r))
    zr <- zr[p[zr] == 0]
    tiny <- rep_len(a0 != 0 & b0 != 0, length(r))[zr]
    r[zr] <- ifelse(tiny, p[zr], p[zr] + rep_len(c, length(r))[zr])
  }
  r
}

.pw_fmaxnm <- function(a, b) if (is.nan(a)) b else if (is.nan(b)) a else max(a, b)
.pw_fminnm <- function(a, b) if (is.nan(a)) b else if (is.nan(b)) a else min(a, b)

# elementwise fmaxnm with a scalar b that is not NaN
.pw_vfmaxnm <- function(a, b) {
  r <- pmax(a, b)
  r[is.nan(a)] <- b
  r
}

# elements 1..m of a loop that gfortran runs in two lane vector form when
# m - 1 > t (all but the last element of an odd m)
.pw_isvec <- function(m, t) {
  idx <- seq_len(m)
  (m - 1 > t) & (idx <= 2 * (m %/% 2))
}

# first index of the largest (smallest) element that is not NaN, 1 when all
# elements are NaN (maxloc and minloc of gfortran)
.pw_maxloc <- function(v) {
  k <- which.max(v)
  if (length(k)) k else 1L
}
.pw_minloc <- function(v) {
  k <- which.min(v)
  if (length(k)) k else 1L
}

# linalg: inner products and matrix products as fused chains
.pw_inprod <- function(x, y) {
  z <- 0
  for (i in seq_along(x)) z <- .pw_fma(x[i], y[i], z)
  z
}
.pw_matprod12 <- function(x, Y) {
  z <- numeric(ncol(Y))
  for (i in seq_along(x)) z <- .pw_fma(x[i], Y[i, ], z)
  z
}
.pw_matprod21 <- function(X, y) {
  z <- numeric(nrow(X))
  for (j in seq_along(y)) z <- .pw_fma(X[, j], y[j], z)
  z
}
.pw_matprod22 <- function(X, Y) {
  nr <- nrow(X)
  nc <- ncol(Y)
  Z <- numeric(nr * nc)
  for (i in seq_len(ncol(X))) Z <- .pw_fma(rep(X[, i], times = nc), rep(Y[i, ], each = nr), Z)
  matrix(Z, nr, nc)
}
.pw_sumabs <- function(x) {
  z <- 0
  for (v in x) z <- z + abs(v)
  z
}
# norm(x) of linalg (2 norm without scaling)
.pw_norm <- function(x) {
  if (!length(x)) return(0)
  if (!all(is.finite(x))) return(.pw_sumabs(x))
  if (!any(abs(x) > 0)) return(0)
  sqrt(.pw_inprod(x, x))
}
# projection of x onto the span of v (project1 of linalg)
.pw_project1 <- function(x, v) {
  if (isTRUE(all(x == 0)) || isTRUE(all(v == 0))) return(numeric(length(x)))
  u <- v / sqrt(.pw_inprod(v, v))
  sc <- max(abs(x))
  (.pw_inprod(x / sc, u) * sc) * u
}
.pw_linspace <- function(a, b, n) {
  h <- (b - a) / (n - 1)
  r <- .pw_fma(seq_len(n) - 1, h, a)
  r[n] <- b
  r
}
# plane rotation G with G %*% x = c(r, 0) (planerot of linalg)
.pw_planerot <- function(x) {
  if (is.nan(x[1]) || is.nan(x[2])) {
    cs <- 1; sn <- 0
  } else if (is.infinite(x[1]) && is.infinite(x[2])) {
    cs <- sign(x[1]) / sqrt(2); sn <- sign(x[2]) / sqrt(2)
  } else if (abs(x[1]) <= 0 && abs(x[2]) <= 0) {
    cs <- 1; sn <- 0
  } else if (abs(x[2]) <= .pr_eps * abs(x[1])) {
    cs <- if (x[1] >= 0) 1 else -1; sn <- 0
  } else if (abs(x[1]) <= .pr_eps * abs(x[2])) {
    cs <- 0; sn <- if (x[2] >= 0) 1 else -1
  } else if (all(abs(x) > sqrt(.Machine$double.xmin)) && all(abs(x) < sqrt(.pr_realmax / 2.1))) {
    r <- sqrt(.pw_fma(x[2], x[2], x[1] * x[1]))
    cs <- x[1] / r; sn <- x[2] / r
  } else if (abs(x[1]) > abs(x[2])) {
    t <- x[2] / x[1]
    u <- max(1, abs(t), sqrt(.pw_fma(t, t, 1))) * (if (x[1] >= 0) 1 else -1)
    cs <- 1 / u; sn <- t / u
  } else {
    t <- x[1] / x[2]
    u <- max(1, abs(t), sqrt(.pw_fma(t, t, 1))) * (if (x[2] >= 0) 1 else -1)
    cs <- t / u; sn <- 1 / u
  }
  matrix(c(cs, -sn, sn, cs), 2L, 2L)
}

# angle in [0, 2 pi) minimizing (maximizing the modulus of) fun(angle, args)
# on a grid of gs points, refined by a parabola through the neighbors
# (circle_min and circle_maxabs of univariate); fun takes a vector of angles
.pw_circle_pick <- function(fun, args, gs, maxabs) {
  ang <- .pw_linspace(0, .pw_twopi, gs + 1L)
  fval <- fun(ang[seq_len(gs)], args)
  if (all(is.nan(fval))) return(0)
  i <- if (maxabs) which.max(abs(fval)) else which.min(fval)
  fm <- fval[i]
  fprev <- fval[(i - 2L) %% gs + 1L]
  fnext <- fval[i %% gs + 1L]
  if (isTRUE(abs(fprev - fnext) > 0)) {
    a <- fprev - fm
    b <- fnext - fm
    step <- ((a - b) * 0.5) / (a + b)
    if (is.finite(step) && abs(step) > 0) return(((i - 1) + step) * (.pw_twopi / gs))
  }
  ang[i]
}

# powalg: the H matrix [BMAT, ZMAT, IDZ] of the interpolation system
.pw_omega_col <- function(idz, zmat, k) {
  zk <- zmat[k, ]
  if (idz > 1L) zk[seq_len(idz - 1L)] <- -zk[seq_len(idz - 1L)]
  .pw_matprod21(zmat, zk)
}
.pw_omega_mul <- function(idz, zmat, x) {
  xz <- .pw_matprod12(x, zmat)
  if (idz > 1L) xz[seq_len(idz - 1L)] <- -xz[seq_len(idz - 1L)]
  .pw_matprod21(zmat, xz)
}
# Hessian of the model (HQ plus sum of PQ(k) XPT(:, k) XPT(:, k)') times x
.pw_hess_mul <- function(x, xpt, pq, hq = NULL) {
  y <- .pw_matprod21(xpt, .pw_matprod12(x, xpt) * pq)
  if (!is.null(hq)) for (j in seq_along(x)) y <- .pw_fma(hq[, j], x[j], y)
  y
}
# Q(xopt + d) - Q(xopt) of the quadratic model [gq, pq, hq] (quadinc)
.pw_quadinc <- function(d, xpt, gq, pq, hq) {
  w <- .pw_matprod12(d, xpt)
  .pw_fma(0.5, .pw_inprod(w, pq * w), .pw_inprod(d, .pw_fma(0.5, .pw_matprod21(hq, d), gq)))
}
.pw_wcheck <- function(kref, d, xpt) {
  w <- .pw_matprod12(d, xpt)
  w * .pw_fma(0.5, w, .pw_matprod12(xpt[, kref], xpt))
}
.pw_calvlag <- function(kref, bmat, d, xpt, zmat, idz) {
  npt <- ncol(xpt)
  wcheck <- .pw_wcheck(kref, d, xpt)
  v1 <- .pw_omega_mul(idz, zmat, wcheck) + .pw_matprod12(d, bmat[, seq_len(npt), drop = FALSE])
  v1[kref] <- v1[kref] + 1
  c(v1, .pw_matprod21(bmat, c(wcheck, d)))
}
.pw_calbeta <- function(kref, bmat, d, xpt, zmat, idz) {
  npt <- ncol(xpt)
  xref <- xpt[, kref]
  wcheck <- .pw_wcheck(kref, d, xpt)
  vlag1 <- .pw_omega_mul(idz, zmat, wcheck) + .pw_matprod12(d, bmat[, seq_len(npt), drop = FALSE])
  vlag2 <- .pw_matprod21(bmat, c(wcheck, d))
  dxref <- .pw_inprod(d, xref)
  dsq <- .pw_inprod(d, d)
  xrefsq <- .pw_inprod(xref, xref)
  .pw_fma(dxref, dxref, dsq * .pw_fma(dsq, 0.5, xrefsq + dxref + dxref)) - .pw_inprod(vlag2, d) - .pw_inprod(vlag1, wcheck)
}
.pw_calden <- function(kref, bmat, d, xpt, zmat, idz) {
  npt <- ncol(xpt)
  vlag <- .pw_calvlag(kref, bmat, d, xpt, zmat, idz)
  beta <- .pw_calbeta(kref, bmat, d, xpt, zmat, idz)
  m <- ncol(zmat)
  hdiag <- vapply(seq_len(npt), function(k) {
    zp <- zmat[k, seq.int(idz, length.out = m - idz + 1L)]
    zn <- zmat[k, seq_len(idz - 1L)]
    .pw_inprod(zp, zp) - .pw_inprod(zn, zn)
  }, 0)
  .pw_fma(vlag[seq_len(npt)], vlag[seq_len(npt)], hdiag * beta)
}

# update of [BMAT, ZMAT, IDZ] when XPT(:, knew) moves to XPT(:, kref) + d
# (updateh of powalg); returns list(idz, bmat, zmat, info)
.pw_updateh <- function(knew, kref, d, xpt, idz, bmat, zmat) {
  n <- nrow(xpt)
  npt <- ncol(xpt)
  hcol <- .pw_omega_col(idz, zmat, knew)
  alpha <- hcol[knew]
  vlag <- .pw_calvlag(kref, bmat, d, xpt, zmat, idz)
  beta <- .pw_calbeta(kref, bmat, d, xpt, zmat, idz)
  tau <- vlag[knew]
  vlag[knew] <- tau - 1
  denom <- .pw_fma(alpha, beta, tau * tau)
  if (!(is.finite(sum(abs(hcol)) + sum(abs(vlag)) + abs(beta)) && isTRUE(abs(denom) > 0)))
    return(list(idz = idz, bmat = bmat, zmat = zmat, info = 7L))
  w <- c(hcol, bmat[, knew])
  jj <- seq_len(n)
  wn <- w[npt + jj]
  vn <- vlag[npt + jj]
  v1 <- ifelse(.pw_isvec(n, 5), .pw_fma(-tau, wn, alpha * vn), .pw_fma(alpha, vn, -(tau * wn))) / denom
  v2 <- .pw_fma(-beta, wn, -(tau * vn)) / denom
  bmat <- (bmat + outer(v1, vlag)) + outer(v2, w)
  B <- bmat[, npt + jj, drop = FALSE]
  B[upper.tri(B)] <- t(B)[upper.tri(B)]
  bmat[, npt + jj] <- B
  ztest <- 1e-20 * max(abs(zmat))
  jl <- 1L
  for (j in seq_len(npt - n - 1L)[-1L]) {
    if (j == idz) {
      jl <- idz
      next
    }
    if (isTRUE(abs(zmat[knew, j]) > ztest)) {
      grot <- .pw_planerot(c(zmat[knew, jl], zmat[knew, j]))
      zmat[, c(jl, j)] <- .pw_matprod22(zmat[, c(jl, j), drop = FALSE], t(grot))
    }
    zmat[knew, j] <- 0
  }
  if (jl == 1L) {
    temp <- sqrt(abs(denom))
    tempb <- zmat[knew, 1L] / temp
    tempa <- tau / temp
    zmat[, 1L] <- ifelse(.pw_isvec(npt, 3), .pw_fma(-tempb, vlag[seq_len(npt)], zmat[, 1L] * tempa),
                         .pw_fma(tempa, zmat[, 1L], -(tempb * vlag[seq_len(npt)])))
  } else {
    ja <- if (beta >= 0) jl else 1L
    jb <- jl + 1L - ja
    temp <- zmat[knew, jb] / denom
    tempa <- temp * beta
    tempb <- temp * tau
    temp <- zmat[knew, ja]
    scala <- 1 / sqrt(.pw_fma(abs(beta), temp * temp, tau * tau))
    scalb <- scala * sqrt(abs(denom))
    vk <- vlag[seq_len(npt)]
    zja <- ifelse(.pw_isvec(npt, 5), .pw_fma(-temp, vk, zmat[, ja] * tau), .pw_fma(tau, zmat[, ja], -(temp * vk))) * scala
    zjb <- .pw_fma(-tempb, vk, .pw_fma(-tempa, w[seq_len(npt)], zmat[, jb])) * scalb
    zmat[, ja] <- zja
    zmat[, jb] <- zjb
  }
  nzc <- npt - n - 1L
  if (denom < 0) {
    dec <- FALSE
    if (idz == 1L) {
      idz <- 2L
    } else if (nzc > idz) {
      if (beta < 0) idz <- idz + 1L else dec <- TRUE
    } else if (nzc < idz) {
      dec <- beta >= 0
    } else {
      dec <- TRUE
    }
    if (dec) {
      idz <- idz - 1L
      if (idz > 1L) zmat[, c(1L, idz)] <- zmat[, c(idz, 1L)]
    }
  }
  list(idz = idz, bmat = bmat, zmat = zmat, info = 0L)
}

# shift of XBASE to XBASE + XOPT with the matching revision of BMAT and HQ
# (shiftbase_lfqint); returns list(xbase, xpt, bmat, hq)
.pw_shiftbase <- function(kopt, xbase, xpt, zmat, bmat, pq, hq, idz) {
  n <- nrow(xpt)
  npt <- ncol(xpt)
  xopt <- xpt[, kopt]
  xoptsq <- .pw_inprod(xopt, xopt)
  yh <- xpt - xopt * 0.5
  sxpt <- .pw_matprod12(xopt, xpt) - 0.5 * xoptsq
  qxoptq <- 0.25 * xoptsq
  ymat <- .pw_fma(rep(sxpt, each = n), yh, qxoptq * xopt)
  ymat <- matrix(ymat, n, npt)
  jj <- seq_len(n)
  kk <- seq_len(npt)
  M1 <- .pw_matprod22(bmat[, kk, drop = FALSE], t(ymat))
  bmat[, npt + jj] <- bmat[, npt + jj] + (M1 + t(M1))
  yz <- .pw_matprod22(ymat, zmat)
  yzc <- yz
  if (idz > 1L) yzc[, seq_len(idz - 1L)] <- -yzc[, seq_len(idz - 1L)]
  bmat[, npt + jj] <- bmat[, npt + jj] + .pw_matprod22(yz, t(yzc))
  bmat[, kk] <- bmat[, kk] + .pw_matprod22(yzc, t(zmat))
  spq <- 0
  for (pqk in pq) spq <- spq + pqk
  v <- .pw_fma(-xopt, spq * 0.5, .pw_matprod21(xpt, pq))
  M4 <- outer(v, xopt)
  hq <- hq + (M4 + t(M4))
  xpt <- xpt - xopt
  xpt[, kopt] <- 0
  list(xbase = xbase + xopt, xpt = xpt, bmat = bmat, hq = hq)
}

# squared distances of the columns of xpt from x, each an fma chain
.pw_distsq <- function(xpt, x) {
  z <- numeric(ncol(xpt))
  for (i in seq_len(nrow(xpt))) {
    dv <- xpt[i, ] - x[i]
    z <- .pw_fma(dv, dv, z)
  }
  z
}
# A + alpha x x' with the lower triangle computed and mirrored (r1 of linalg)
.pw_r1sym <- function(A, alpha, x) {
  n <- length(x)
  B <- matrix(.pw_fma(rep(alpha * x, times = n), rep(x, each = n), A), n, n)
  B[upper.tri(B)] <- t(B)[upper.tri(B)]
  B
}

# evaluate of PRIMA: f at x, moderated to [-REALMAX, FUNCMAX] (NaN and NA
# become FUNCMAX, which is 2^100 in libprima); x containing NaN is not
# evaluated
.pw_funcmax <- 2^100
.pw_evaluatef <- function(fun, x) {
  if (any(is.nan(x))) return(sum(x))
  f <- fun(pmin(pmax(x, -.pr_realmax), .pr_realmax))
  if (is.na(f)) .pw_funcmax else min(max(f, -.pr_realmax), .pw_funcmax)
}
.pw_checkexit <- function(maxfun, nf, f, ftarget, x) {
  info <- 0L
  if (any(is.nan(x) | is.infinite(x))) info <- -1L
  if (is.nan(f) || f == Inf) info <- -2L
  if (isTRUE(f <= ftarget)) info <- 1L
  if (nf >= maxfun) info <- 3L
  info
}
# the messages of PRIMA (the fmsg, rhomsg and retmsg routines of its
# message module) for iprint = 1 (exit), 2 (each new rho) and 3 (each
# evaluation); reals are written in the Fortran format ES23.15E3
.pw_ffmt <- function(v) {
  s <- sprintf("%.15e", v)
  e <- suppressWarnings(as.integer(sub("^.*e", "", s)))
  s <- ifelse(is.finite(v), paste0(sub("e.*$", "", s), "E", ifelse(e < 0L, "-", "+"), sprintf("%03d", abs(e))),
              ifelse(is.nan(v), "NaN", ifelse(v > 0, "Infinity", "-Infinity")))
  sprintf("%23s", s)
}
.pw_xlines <- function(x) {
  f <- .pw_ffmt(x)
  if (length(f) <= 2L) return(paste0("The corresponding X is: ", paste(f, collapse = "  ")))
  c("The corresponding X is:", vapply(split(f, (seq_along(f) - 1L) %/% 4L), paste, "", collapse = "  "))
}
.pw_say <- function(lines) cat(paste0(lines, "\n"), sep = "")
.pw_fmsg <- function(iprint, state, nf, delta, f, x) {
  if (iprint >= 3L) {
    .pw_say(c("", paste0(state, " step with radius = ", .pw_ffmt(delta)),
              paste0("Function number ", nf, "   F = ", .pw_ffmt(f)), .pw_xlines(x)))
  }
}
.pw_rhomsg <- function(iprint, nf, delta, f, rho, x) {
  if (iprint >= 2L) {
    .pw_say(c(if (iprint >= 3L) "", "", paste0("New RHO = ", .pw_ffmt(rho), "   Delta = ", .pw_ffmt(delta)),
              paste0("Number of function values = ", nf, "   Least value of F = ", .pw_ffmt(f)), .pw_xlines(x)))
  }
}
.pw_retmsg <- function(solver, iprint, info, nf, f, x) {
  if (iprint >= 1L) {
    why <- switch(as.character(info),
                  "0" = "the trust region radius reaches its lower bound.",
                  "1" = "the target function value is achieved.",
                  "2" = "a trust region step has failed to reduce the quadratic model.",
                  "3" = "the maximal number of function evaluations has been reached.",
                  "20" = "the maximal number of trust region iterations has been reached.",
                  "-1" = "NaN or Inf occurs in x.",
                  "-2" = "the objective or constraint functions return NaN or +Inf.",
                  "-3" = "NaN or Inf occurs in the models.",
                  "7" = "rounding errors are becoming damaging.",
                  "UNKNOWN EXIT FLAG")
    .pw_say(c(if (iprint >= 2L) "", "", paste0("Return from ", solver, " because ", why),
              paste0("Number of function values = ", nf, "   Least value of F = ", .pw_ffmt(f)), .pw_xlines(x)))
  }
}

# moderatex: NaN becomes 0, the rest is clipped to [-REALMAX, REALMAX]
.pw_moderatex <- function(x) {
  x[is.nan(x)] <- 0
  pmin(pmax(x, -.pr_realmax), .pr_realmax)
}

# preproc of PRIMA for UOBYQA, NEWUOA, BOBYQA and LINCOA: revises invalid
# options (the parameters eta1, eta2, gamma1 and gamma2 keep their defaults,
# which PRIMA.jl cannot change)
.pw_preproc <- function(solver, n, maxfun, ftarget, rhobeg, rhoend, npt = NULL) {
  min_maxfun <- if (solver == "UOBYQA") ((n + 1L) * (n + 2L)) %/% 2L + 1L else n + 3L
  warn <- function(...) warning(solver, ": Invalid ", ..., ".", call. = FALSE)
  if (maxfun < min_maxfun) {
    maxfun <- min_maxfun
    warn("MAXFUN; it should be at least ", if (solver == "UOBYQA") "(N+1)(N+2)/2 + 1" else "N + 3",
         "; it is set to ", maxfun)
  }
  if (is.nan(ftarget)) {
    ftarget <- -.pr_realmax
    warn("FTARGET; it should be a real number; it is set to ", .pw_ffmt(ftarget))
  }
  if (!is.null(npt) && (npt < n + 2L || npt > min(maxfun - 1L, ((n + 2L) * (n + 1L)) %/% 2L))) {
    npt <- min(maxfun - 1L, 2L * n + 1L)
    warn("NPT; it should be an integer in the interval [N+2, (N+1)(N+2)/2] and less than MAXFUN; it is set to ",
         npt)
  }
  if (isTRUE(abs(rhobeg - rhoend) < 1e2 * .pr_eps * max(abs(rhobeg), 1))) rhoend <- rhobeg
  if (!isTRUE(rhobeg > 0) || !is.finite(rhobeg)) {
    rhobeg <- if (is.finite(rhoend) && rhoend > 0) max(10 * rhoend, 1) else 1
    warn("RHOBEG; it should be a positive number; it is set to ", .pw_ffmt(rhobeg))
  }
  if (!isTRUE(rhoend > 0) || rhobeg < rhoend || !is.finite(rhoend)) {
    rhoend <- max(.pr_eps, min(0.1 * rhobeg, 1e-6))
    warn("RHOEND; it should be a positive number and RHOEND <= RHOBEG; it is set to ", .pw_ffmt(rhoend))
  }
  list(maxfun = as.integer(maxfun), ftarget = ftarget, rhobeg = rhobeg, rhoend = rhoend,
       npt = if (!is.null(npt)) as.integer(npt))
}
