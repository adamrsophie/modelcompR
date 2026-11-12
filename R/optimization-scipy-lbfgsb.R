# L-BFGS-B 3.0 as translated to C in SciPy 1.18.1 (scipy/optimize/src/lbfgsb.c)
# with the driver _minimize_lbfgsb of _lbfgsb_py.py. The reverse communication
# of setulb is unrolled into a direct loop with the same evaluation sequence.
# Variable indices held in index arrays are 1-based; positions such as head,
# itail, nbreak and nleft keep the 0-based values of the C code.

.lb_rng <- function(a, b) if (a <= b) a:b else integer()

# solve with the upper triangle of U (dtrtrs): transpose = TRUE solves U'x = b
.lb_trsolve <- function(U, b, transpose = FALSE) {
  if (any(diag(U) == 0)) return(NULL)
  backsolve(U, b, transpose = transpose)
}

.lb_bmv <- function(sy, wt, col, v) {
  p <- numeric(2L * col)
  if (col == 0L) return(list(p = p, info = 0L))
  p[col + 1L] <- v[col + 1L]
  for (i in .lb_rng(1L, col - 1L)) {
    i2 <- col + i
    ssum <- 0.0
    for (k in 0:(i - 1L)) ssum <- ssum + sy[i + 1L, k + 1L] * v[k + 1L] / sy[k + 1L, k + 1L]
    p[i2 + 1L] <- v[i2 + 1L] + ssum
  }
  U <- wt[seq_len(col), seq_len(col), drop = FALSE]
  q <- .lb_trsolve(U, p[col + seq_len(col)], transpose = TRUE)
  if (is.null(q)) return(list(p = p, info = 1L))
  p[col + seq_len(col)] <- q
  for (i in seq_len(col)) p[i] <- v[i] / sqrt(sy[i, i])
  q <- .lb_trsolve(U, p[col + seq_len(col)])
  if (is.null(q)) return(list(p = p, info = 1L))
  p[col + seq_len(col)] <- q
  for (i in seq_len(col)) p[i] <- -p[i] / sqrt(sy[i, i])
  for (i in seq_len(col)) {
    ssum <- 0.0
    for (k in .lb_rng(i + 1L, col)) ssum <- ssum + sy[k, i] * p[col + k] / sy[i, i]
    p[i] <- p[i] + ssum
  }
  list(p = p, info = 0L)
}

.lb_hpsolb <- function(n, t, iorder, iheap) {
  if (iheap == 0L) {
    for (k in .lb_rng(2L, n + 1L)) {
      ddum <- t[k]
      indxin <- iorder[k]
      i <- k
      while (i > 1L) {
        j <- i %/% 2L
        if (ddum < t[j]) {
          t[i] <- t[j]
          iorder[i] <- iorder[j]
          i <- j
        } else {
          break
        }
      }
      t[i] <- ddum
      iorder[i] <- indxin
    }
  }
  if (n > 0L) {
    i <- 1L
    out <- t[1L]
    indxout <- iorder[1L]
    ddum <- t[n + 1L]
    indxin <- iorder[n + 1L]
    repeat {
      j <- i + i
      if (j > n) break
      if (t[j + 1L] < t[j]) j <- j + 1L
      if (t[j] < ddum) {
        t[i] <- t[j]
        iorder[i] <- iorder[j]
        i <- j
      } else {
        break
      }
    }
    t[i] <- ddum
    iorder[i] <- indxin
    t[n + 1L] <- out
    iorder[n + 1L] <- indxout
  }
  list(t = t, iorder = iorder)
}

.sp_lbfgsb <- function(fun, x0, jac = NULL, bounds = NULL, maxcor = 10L, ftol = 2.2204460492503131e-09, gtol = 1e-5,
                       eps = 1e-8, maxfun = 15000L, maxiter = 15000L, maxls = 20L, finite_diff_rel_step = NULL, ...) {
  m <- as.integer(maxcor)
  pgtol <- gtol
  epsmach <- 2.220446049250313e-016
  factr <- ftol / .Machine$double.eps
  x0 <- as.numeric(x0)
  n <- length(x0)
  lbv <- rep(-Inf, n)
  ubv <- rep(Inf, n)
  if (!is.null(bounds)) {
    lbv <- rep_len(as.numeric(bounds$lb), n)
    ubv <- rep_len(as.numeric(bounds$ub), n)
    if (any(lbv > ubv)) stop("LBFGSB - one of the lower bounds is greater than an upper bound.", call. = FALSE)
    x0 <- pmin(pmax(x0, lbv), ubv)
  }
  if (!(maxls > 0)) stop("maxls must be positive.", call. = FALSE)
  sf <- .sp_prepare_scalar_function(fun, x0, jac, lb = lbv, ub = ubv, epsilon = eps,
                                    finite_diff_rel_step = finite_diff_rel_step)
  nbd <- integer(n)
  l <- numeric(n)
  u <- numeric(n)
  for (i in seq_len(n)) {
    L <- is.finite(lbv[i]); U <- is.finite(ubv[i])
    if (L) l[i] <- lbv[i]
    if (U) u[i] <- ubv[i]
    nbd[i] <- if (!L && !U) 0L else if (L && !U) 1L else if (L && U) 2L else 3L
  }
  if (m <= 0L) stop("maxcor must be positive", call. = FALSE)
  if (factr < 0) stop("factr must be nonnegative", call. = FALSE)

  # workspace
  ws <- matrix(0, n, m); wy <- matrix(0, n, m)
  sy <- matrix(0, m, m); ss <- matrix(0, m, m); wt <- matrix(0, m, m)
  wn <- matrix(0, 2L * m, 2L * m); wn1 <- matrix(0, 2L * m, 2L * m)
  z <- numeric(n); r <- numeric(n); d <- numeric(n); t <- numeric(n)
  wa <- numeric(8L * m)
  index <- integer(n); iwhere <- integer(n); indx2 <- integer(n)

  x <- x0
  col <- 0L; head <- 0L; theta <- 1.0; iupdat <- 0L; updatd <- FALSE; iback <- 0L; itail <- 0L
  iword <- 0L; nact <- 0L; ileave <- 0L; nenter <- 0L; fold <- 0.0; dnorm <- 0.0; gd <- 0.0; sbgnrm <- 0.0
  stp <- 0.0; gdold <- 0.0; dtd <- 0.0; iter <- 0L; nfgv <- 0L; nseg <- 0L; nintol <- 0L; nskip <- 0L
  nfree <- n; ifun <- 0L; tol <- factr * epsmach; info <- 0L; stpmx <- 0.0

  # active: project x0 and initialize iwhere
  prjctd <- FALSE; cnstnd <- FALSE; boxed <- TRUE
  for (i in seq_len(n)) {
    if (nbd[i] > 0L) {
      if (nbd[i] <= 2L && x[i] <= l[i]) {
        if (x[i] < l[i]) { prjctd <- TRUE; x[i] <- l[i] }
      } else if (nbd[i] >= 2L && x[i] >= u[i]) {
        if (x[i] > u[i]) { prjctd <- TRUE; x[i] <- u[i] }
      }
    }
  }
  for (i in seq_len(n)) {
    if (nbd[i] != 2L) boxed <- FALSE
    if (nbd[i] == 0L) {
      iwhere[i] <- -1L
    } else {
      cnstnd <- TRUE
      iwhere[i] <- if (nbd[i] == 2L && u[i] - l[i] <= 0.0) 3L else 0L
    }
  }

  projgr <- function() {
    s <- 0.0
    for (i in seq_len(n)) {
      gi <- g[i]
      if (is.nan(gi)) return(gi)
      if (nbd[i] != 0L) {
        if (gi < 0.0) {
          if (nbd[i] >= 2L) gi <- max(x[i] - u[i], gi)
        } else {
          if (nbd[i] <= 2L) gi <- min(x[i] - l[i], gi)
        }
      }
      s <- max(s, abs(gi))
    }
    s
  }

  reset_memory <- function() {
    info <<- 0L; col <<- 0L; head <<- 0L; theta <<- 1.0; iupdat <<- 0L; updatd <<- FALSE
  }

  cauchy <- function() {
    f1 <- 0.0; tl <- 0.0; tu <- 0.0
    if (sbgnrm <= 0.0) {
      z <<- x
      return(invisible())
    }
    bnded <- TRUE
    nfree_c <- n
    nbreak <- -1L
    ibkmin <- 0L
    bkmin <- 0.0
    col2 <- 2L * col
    P <- seq_len(col2)
    C <- 2L * m + P
    wa[P] <<- 0.0
    for (i in seq_len(n)) {
      neggi <- -g[i]
      if (iwhere[i] != 3L && iwhere[i] != -1L) {
        if (nbd[i] <= 2L) tl <- x[i] - l[i]
        if (nbd[i] >= 2L) tu <- u[i] - x[i]
        xlower <- nbd[i] <= 2L && tl <= 0.0
        xupper <- nbd[i] >= 2L && tu <= 0.0
        iwhere[i] <<- 0L
        if (xlower) {
          if (neggi <= 0.0) iwhere[i] <<- 1L
        } else if (xupper) {
          if (neggi >= 0.0) iwhere[i] <<- 2L
        } else {
          if (abs(neggi) <= 0.0) iwhere[i] <<- -3L
        }
      }
      pointr <- head
      if (iwhere[i] != 0L && iwhere[i] != -1L) {
        d[i] <<- 0.0
      } else {
        d[i] <<- neggi
        f1 <- f1 - neggi * neggi
        for (j in seq_len(col)) {
          wa[j] <<- wa[j] + wy[i, pointr + 1L] * neggi
          wa[col + j] <<- wa[col + j] + ws[i, pointr + 1L] * neggi
          pointr <- (pointr + 1L) %% m
        }
        if (nbd[i] <= 2L && nbd[i] != 0L && neggi < 0.0) {
          nbreak <- nbreak + 1L
          indx2[nbreak + 1L] <<- i
          t[nbreak + 1L] <<- tl / (-neggi)
          if (nbreak == 0L || t[nbreak + 1L] < bkmin) {
            bkmin <- t[nbreak + 1L]
            ibkmin <- nbreak
          }
        } else if (nbd[i] >= 2L && neggi > 0.0) {
          nbreak <- nbreak + 1L
          indx2[nbreak + 1L] <<- i
          t[nbreak + 1L] <<- tu / neggi
          if (nbreak == 0L || t[nbreak + 1L] < bkmin) {
            bkmin <- t[nbreak + 1L]
            ibkmin <- nbreak
          }
        } else {
          nfree_c <- nfree_c - 1L
          indx2[nfree_c + 1L] <<- i
          if (abs(neggi) > 0.0) bnded <- FALSE
        }
      }
    }
    if (theta != 1.0 && col > 0L) wa[col + seq_len(col)] <<- theta * wa[col + seq_len(col)]
    z <<- x
    if (nbreak == -1L && nfree_c == n) return(invisible())
    wa[C] <<- 0.0
    f2 <- -theta * f1
    f2_org <- f2
    if (col > 0L) {
      bm <- .lb_bmv(sy, wt, col, wa[P])
      wa[6L * m + P] <<- bm$p
      if (bm$info != 0L) {
        info <<- bm$info
        return(invisible())
      }
      f2 <- f2 - sum(bm$p * wa[P])
    }
    dtm <- -f1 / f2
    tsum <- 0.0
    nseg <<- 1L
    nleft <- nbreak
    it <- 0L
    tj <- 0.0
    while (nbreak > -1L) {
      tj0 <- tj
      if (it == 0L) {
        tj <- bkmin
        ibp <- indx2[ibkmin + 1L]
      } else {
        if (it == 1L && ibkmin != nbreak) {
          t[ibkmin + 1L] <<- t[nbreak + 1L]
          indx2[ibkmin + 1L] <<- indx2[nbreak + 1L]
        }
        hp <- .lb_hpsolb(nleft, t, indx2, it - 1L)
        t <<- hp$t
        indx2 <<- hp$iorder
        tj <- t[nleft + 1L]
        ibp <- indx2[nleft + 1L]
      }
      dt <- tj - tj0
      if (dtm < dt) {
        if (dtm <= 0.0) dtm <- 0.0
        tsum <- tsum + dtm
        z <<- z + tsum * d
        if (col > 0L) wa[C] <<- wa[C] + dtm * wa[P]
        return(invisible())
      }
      tsum <- tsum + dt
      nleft <- nleft - 1L
      it <- it + 1L
      dibp <- d[ibp]
      d[ibp] <<- 0.0
      if (dibp > 0.0) {
        zibp <- u[ibp] - x[ibp]
        z[ibp] <<- u[ibp]
        iwhere[ibp] <<- 2L
      } else {
        zibp <- l[ibp] - x[ibp]
        z[ibp] <<- l[ibp]
        iwhere[ibp] <<- 1L
      }
      # the C code compares the 0-based nbreak with n here, so the test never holds
      if (nleft == -1L && nbreak == n) {
        dtm <- dt
        if (col > 0L) wa[C] <<- wa[C] + dtm * wa[P]
        return(invisible())
      }
      nseg <<- nseg + 1L
      dibp2 <- dibp * dibp
      f1 <- f1 + dt * f2 + dibp2 - theta * dibp * zibp
      f2 <- f2 - theta * dibp2
      if (col > 0L) {
        wa[C] <<- wa[C] + dt * wa[P]
        pointr <- head
        wbp <- numeric(col2)
        for (j in seq_len(col)) {
          wbp[j] <- wy[ibp, pointr + 1L]
          wbp[col + j] <- theta * ws[ibp, pointr + 1L]
          pointr <- (pointr + 1L) %% m
        }
        wa[4L * m + P] <<- wbp
        bm <- .lb_bmv(sy, wt, col, wbp)
        wa[6L * m + P] <<- bm$p
        if (bm$info != 0L) {
          info <<- bm$info
          return(invisible())
        }
        v <- bm$p
        wmc <- sum(wa[C] * v)
        wmp <- sum(wa[P] * v)
        wmw <- sum(wbp * v)
        wa[P] <<- wa[P] + (-dibp) * wbp
        f1 <- f1 + dibp * wmc
        f2 <- f2 + dibp * 2.0 * wmp - dibp2 * wmw
      }
      f2 <- max(epsmach * f2_org, f2)
      if (nleft >= 0L) {
        dtm <- -f1 / f2
      } else if (bnded) {
        dtm <- 0.0
        break
      } else {
        dtm <- -f1 / f2
        break
      }
    }
    if (dtm <= 0.0) dtm <- 0.0
    tsum <- tsum + dtm
    z <<- z + tsum * d
    if (col > 0L) wa[C] <<- wa[C] + dtm * wa[P]
    invisible()
  }

  freev <- function() {
    nenter <<- 0L
    ileave <<- n
    if (iter > 0L && cnstnd) {
      for (i in seq_len(nfree)) {
        k <- index[i]
        if (iwhere[k] > 0L) {
          ileave <<- ileave - 1L
          indx2[ileave + 1L] <<- k
        }
      }
      for (i in .lb_rng(nfree + 1L, n)) {
        k <- index[i]
        if (iwhere[k] <= 0L) {
          indx2[nenter + 1L] <<- k
          nenter <<- nenter + 1L
        }
      }
    }
    wrk <- ileave < n || nenter > 0L || updatd
    nfree <<- 0L
    iact <- n
    for (i in seq_len(n)) {
      if (iwhere[i] <= 0L) {
        index[nfree + 1L] <<- i
        nfree <<- nfree + 1L
      } else {
        iact <- iact - 1L
        index[iact + 1L] <<- i
      }
    }
    wrk
  }

  formk <- function() {
    m2 <- 2L * m
    if (updatd) {
      if (iupdat > m) {
        for (J in .lb_rng(1L, m - 1L)) {
          wn1[J:(m - 1L), J] <<- wn1[(J + 1L):m, J + 1L]
          JS <- m + J
          wn1[JS:(m2 - 1L), JS] <<- wn1[(JS + 1L):m2, JS + 1L]
          wn1[(m + 1L):(m2 - 1L), J] <<- wn1[(m + 2L):m2, J + 1L]
        }
      }
      fr <- index[seq_len(nfree)]
      ac <- index[.lb_rng(nfree + 1L, n)]
      iy <- col - 1L
      is <- m + col - 1L
      ipntr <- head + col - 1L
      if (ipntr >= m) ipntr <- ipntr - m
      jpntr <- head
      for (jy in 0:(col - 1L)) {
        js <- m + jy
        temp1 <- 0.0; temp2 <- 0.0; temp3 <- 0.0
        for (k1 in fr) temp1 <- temp1 + wy[k1, ipntr + 1L] * wy[k1, jpntr + 1L]
        for (k1 in ac) {
          temp2 <- temp2 + ws[k1, ipntr + 1L] * ws[k1, jpntr + 1L]
          temp3 <- temp3 + ws[k1, ipntr + 1L] * wy[k1, jpntr + 1L]
        }
        wn1[iy + 1L, jy + 1L] <<- temp1
        wn1[is + 1L, js + 1L] <<- temp2
        wn1[is + 1L, jy + 1L] <<- temp3
        jpntr <- (jpntr + 1L) %% m
      }
      jy <- col - 1L
      jpntr <- head + col - 1L
      if (jpntr >= m) jpntr <- jpntr - m
      ipntr <- head
      for (i in 0:(col - 1L)) {
        is <- m + i
        temp3 <- 0.0
        for (k1 in fr) temp3 <- temp3 + ws[k1, ipntr + 1L] * wy[k1, jpntr + 1L]
        ipntr <- (ipntr + 1L) %% m
        wn1[is + 1L, jy + 1L] <<- temp3
      }
      upcl <- col - 1L
    } else {
      upcl <- col
    }
    ent <- indx2[seq_len(nenter)]
    lv <- indx2[.lb_rng(ileave + 1L, n)]
    ipntr <- head
    for (iy in .lb_rng(0L, upcl - 1L)) {
      is <- m + iy
      jpntr <- head
      for (jy in 0:iy) {
        js <- m + jy
        temp1 <- 0.0; temp2 <- 0.0; temp3 <- 0.0; temp4 <- 0.0
        for (k1 in ent) {
          temp1 <- temp1 + wy[k1, ipntr + 1L] * wy[k1, jpntr + 1L]
          temp2 <- temp2 + ws[k1, ipntr + 1L] * ws[k1, jpntr + 1L]
        }
        for (k1 in lv) {
          temp3 <- temp3 + wy[k1, ipntr + 1L] * wy[k1, jpntr + 1L]
          temp4 <- temp4 + ws[k1, ipntr + 1L] * ws[k1, jpntr + 1L]
        }
        wn1[iy + 1L, jy + 1L] <<- wn1[iy + 1L, jy + 1L] + temp1 - temp3
        wn1[is + 1L, js + 1L] <<- wn1[is + 1L, js + 1L] - temp2 + temp4
        jpntr <- (jpntr + 1L) %% m
      }
      ipntr <- (ipntr + 1L) %% m
    }
    ipntr <- head
    for (is in .lb_rng(m, m + upcl - 1L)) {
      jpntr <- head
      for (jy in .lb_rng(0L, upcl - 1L)) {
        temp1 <- 0.0; temp3 <- 0.0
        for (k1 in ent) temp1 <- temp1 + ws[k1, ipntr + 1L] * wy[k1, jpntr + 1L]
        for (k1 in lv) temp3 <- temp3 + ws[k1, ipntr + 1L] * wy[k1, jpntr + 1L]
        if (is <= jy + m) {
          wn1[is + 1L, jy + 1L] <<- wn1[is + 1L, jy + 1L] + temp1 - temp3
        } else {
          wn1[is + 1L, jy + 1L] <<- wn1[is + 1L, jy + 1L] - temp1 + temp3
        }
        jpntr <- (jpntr + 1L) %% m
      }
      ipntr <- (ipntr + 1L) %% m
    }
    for (iy in 0:(col - 1L)) {
      is <- col + iy
      is1 <- m + iy
      for (jy in 0:iy) {
        js <- col + jy
        js1 <- m + jy
        wn[jy + 1L, iy + 1L] <<- wn1[iy + 1L, jy + 1L] / theta
        wn[js + 1L, is + 1L] <<- wn1[is1 + 1L, js1 + 1L] * theta
      }
      for (jy in .lb_rng(0L, iy - 1L)) wn[jy + 1L, is + 1L] <<- -wn1[is1 + 1L, jy + 1L]
      for (jy in iy:(col - 1L)) wn[jy + 1L, is + 1L] <<- wn1[is1 + 1L, jy + 1L]
      wn[iy + 1L, iy + 1L] <<- wn[iy + 1L, iy + 1L] + sy[iy + 1L, iy + 1L]
    }
    A <- seq_len(col)
    B <- col + seq_len(col)
    ch <- .sp_potrf(wn[A, A, drop = FALSE])
    if (ch$info != 0L) {
      info <<- -1L
      return(invisible())
    }
    wn[A, A] <<- wn[A, A] * lower.tri(wn[A, A]) + ch$U
    wn[A, B] <<- backsolve(ch$U, wn[A, B, drop = FALSE], transpose = TRUE)
    for (is in B) for (js in is:(2L * col)) wn[is, js] <<- wn[is, js] + sum(wn[A, is] * wn[A, js])
    ch <- .sp_potrf(wn[B, B, drop = FALSE])
    if (ch$info != 0L) {
      info <<- -2L
      return(invisible())
    }
    wn[B, B] <<- wn[B, B] * lower.tri(wn[B, B]) + ch$U
    invisible()
  }

  cmprlb <- function() {
    if (!cnstnd && col > 0L) {
      r <<- -g
      return(invisible())
    }
    fr <- index[seq_len(nfree)]
    r[seq_len(nfree)] <<- -theta * (z[fr] - x[fr]) - g[fr]
    bm <- .lb_bmv(sy, wt, col, wa[2L * m + seq_len(2L * col)])
    wa[seq_len(2L * col)] <<- bm$p
    if (bm$info != 0L) {
      info <<- -8L
      return(invisible())
    }
    pointr <- head
    for (j in seq_len(col)) {
      a1 <- wa[j]
      a2 <- theta * wa[col + j]
      r[seq_len(nfree)] <<- r[seq_len(nfree)] + wy[fr, pointr + 1L] * a1 + ws[fr, pointr + 1L] * a2
      pointr <- (pointr + 1L) %% m
    }
    invisible()
  }

  subsm <- function() {
    nsub <- nfree
    if (nsub <= 0L) return(invisible())
    ind <- index[seq_len(nsub)]
    dd <- r[seq_len(nsub)]
    col2 <- 2L * col
    wv <- numeric(col2)
    pointr <- head
    for (i in seq_len(col)) {
      temp1 <- 0.0; temp2 <- 0.0
      for (j in seq_len(nsub)) {
        temp1 <- temp1 + wy[ind[j], pointr + 1L] * dd[j]
        temp2 <- temp2 + ws[ind[j], pointr + 1L] * dd[j]
      }
      wv[i] <- temp1
      wv[col + i] <- theta * temp2
      pointr <- (pointr + 1L) %% m
    }
    U <- wn[seq_len(col2), seq_len(col2), drop = FALSE]
    q <- .lb_trsolve(U, wv, transpose = TRUE)
    if (is.null(q)) {
      info <<- 1L
      return(invisible())
    }
    wv <- q
    wv[seq_len(col)] <- -wv[seq_len(col)]
    q <- .lb_trsolve(U, wv)
    if (is.null(q)) {
      info <<- 1L
      return(invisible())
    }
    wv <- q
    pointr <- head
    for (jy in seq_len(col)) {
      js <- col + jy
      dd <- dd + (wy[ind, pointr + 1L] * wv[jy] / theta) + (ws[ind, pointr + 1L] * wv[js])
      pointr <- (pointr + 1L) %% m
    }
    dd <- (1.0 / theta) * dd
    iword <<- 0L
    xp <- z
    for (i in seq_len(nsub)) {
      k <- ind[i]
      dk <- dd[i]
      xk <- z[k]
      if (nbd[k] != 0L) {
        if (nbd[k] == 1L) {
          z[k] <<- max(l[k], xk + dk)
          if (z[k] == l[k]) iword <<- 1L
        } else if (nbd[k] == 2L) {
          xk <- max(l[k], xk + dk)
          z[k] <<- min(u[k], xk)
          if (z[k] == l[k] || z[k] == u[k]) iword <<- 1L
        } else if (nbd[k] == 3L) {
          z[k] <<- min(u[k], xk + dk)
          if (z[k] == u[k]) iword <<- 1L
        }
      } else {
        z[k] <<- xk + dk
      }
    }
    r[seq_len(nsub)] <<- dd
    if (iword == 0L) return(invisible())
    dd_p <- 0.0
    for (i in seq_len(n)) dd_p <- dd_p + (z[i] - x[i]) * g[i]
    if (dd_p > 0.0) z <<- xp else return(invisible())
    alpha <- 1.0
    temp1 <- 1.0
    ibd <- 1L
    for (i in seq_len(nsub)) {
      k <- ind[i]
      dk <- dd[i]
      if (nbd[k] != 0L) {
        if (dk < 0.0 && nbd[k] <= 2L) {
          temp2 <- l[k] - z[k]
          if (temp2 >= 0.0) temp1 <- 0.0 else if (dk * alpha < temp2) temp1 <- temp2 / dk
        } else if (dk > 0.0 && nbd[k] >= 2L) {
          temp2 <- u[k] - z[k]
          if (temp2 <= 0.0) temp1 <- 0.0 else if (dk * alpha > temp2) temp1 <- temp2 / dk
        }
        if (temp1 < alpha) {
          alpha <- temp1
          ibd <- i
        }
      }
    }
    if (alpha < 1.0) {
      dk <- dd[ibd]
      k <- ind[ibd]
      if (dk > 0.0) {
        z[k] <<- u[k]
        dd[ibd] <- 0.0
      } else if (dk < 0.0) {
        z[k] <<- l[k]
        dd[ibd] <- 0.0
      }
    }
    z[ind] <<- z[ind] + alpha * dd
    r[seq_len(nsub)] <<- dd
    invisible()
  }

  matupd <- function(rr, dr) {
    if (iupdat <= m) {
      col <<- iupdat
      itail <<- (head + iupdat - 1L) %% m
    } else {
      itail <<- (itail + 1L) %% m
      head <<- (head + 1L) %% m
    }
    ws[, itail + 1L] <<- d
    wy[, itail + 1L] <<- r
    theta <<- rr / dr
    if (iupdat > m) {
      for (j in .lb_rng(1L, col - 1L)) {
        ss[1:j, j] <<- ss[2:(j + 1L), j + 1L]
        sy[j:(col - 1L), j] <<- sy[(j + 1L):col, j + 1L]
      }
    }
    pointr <- head
    for (j in .lb_rng(1L, col - 1L)) {
      sy[col, j] <<- sum(d * wy[, pointr + 1L])
      ss[j, col] <<- sum(ws[, pointr + 1L] * d)
      pointr <- (pointr + 1L) %% m
    }
    ss[col, col] <<- if (stp == 1.0) dtd else stp * stp * dtd
    sy[col, col] <<- dr
  }

  formt <- function() {
    W <- matrix(0, col, col)
    for (j in seq_len(col)) W[1L, j] <- theta * ss[1L, j]
    for (i in .lb_rng(2L, col)) {
      for (j in i:col) {
        ddum <- 0.0
        for (k in seq_len(i - 1L)) ddum <- ddum + ((sy[i, k] * sy[j, k]) / sy[k, k])
        W[i, j] <- ddum + theta * ss[i, j]
      }
    }
    ch <- .sp_potrf(W)
    if (ch$info != 0L) {
      info <<- -3L
      return(invisible())
    }
    wt[seq_len(col), seq_len(col)] <<- ch$U
    invisible()
  }

  fg <- sf$fun_and_grad(x)
  f <- fg$f
  g <- as.numeric(fg$g)
  n_iterations <- 0L
  task <- "FG_START"
  nfgv <- 1L
  sbgnrm <- projgr()
  if (sbgnrm <= pgtol) task <- "CONVERGENCE"
  while (task == "FG_START") {
    iword <- -1L
    if (cnstnd || col == 0L) {
      cauchy()
      if (info != 0L) {
        reset_memory()
        next
      }
      nintol <- nintol + nseg
      wrk <- freev()
      nact <- n - nfree
    } else {
      z <- x
      wrk <- updatd
      nseg <- 0L
    }
    if (nfree != 0L && col != 0L) {
      if (wrk) formk()
      if (info != 0L) {
        reset_memory()
        next
      }
      cmprlb()
      if (info == 0L) subsm()
      if (info != 0L) {
        reset_memory()
        next
      }
    }
    d <- z - x
    # lnsrlb
    dnorm <- .sp_norm2(d)
    dtd <- dnorm^2
    stpmx <- 1.0e+10
    if (cnstnd) {
      if (iter == 0L) {
        stpmx <- 1.0
      } else {
        for (i in seq_len(n)) {
          a1 <- d[i]
          if (nbd[i] != 0L) {
            if (a1 < 0.0 && nbd[i] <= 2L) {
              a2 <- l[i] - x[i]
              if (a2 >= 0.0) stpmx <- 0.0 else if (a1 * stpmx < a2) stpmx <- a2 / a1
            } else if (a1 > 0.0 && nbd[i] >= 2L) {
              a2 <- u[i] - x[i]
              if (a2 <= 0.0) stpmx <- 0.0 else if (a1 * stpmx > a2) stpmx <- a2 / a1
            }
          }
        }
      }
    }
    stp <- if (iter == 0L && !boxed) min(1.0 / dnorm, stpmx) else 1.0
    t <- x
    r <- g
    fold <- f
    ifun <- 0L
    iback <- 0L
    ls_iterate <- .sp_dcsrch_state(1e-3, 0.9, 0.1, 0.0, stpmx)
    ls_task <- "START"
    repeat {
      gd <- sum(g * d)
      if (ifun == 0L) {
        gdold <- gd
        if (gd >= 0.0) {
          info <- -4L
          break
        }
      }
      ls <- ls_iterate(stp, f, gd, ls_task)
      stp <- ls$stp
      ls_task <- ls$task
      if (startsWith(ls_task, "CONV") || startsWith(ls_task, "WARN")) break
      ifun <- ifun + 1L
      nfgv <- nfgv + 1L
      iback <- ifun - 1L
      if (stp == 1.0) {
        x <- z
      } else {
        x <- stp * d + t
        x <- ifelse(nbd == 1L | nbd == 2L, pmax(x, l), x)
        x <- ifelse(nbd == 2L | nbd == 3L, pmin(x, u), x)
      }
      if (iback >= maxls) break
      fg <- sf$fun_and_grad(x)
      f <- fg$f
      g <- as.numeric(fg$g)
    }
    if (info != 0L || iback >= maxls) {
      x <- t
      g <- r
      f <- fold
      if (col == 0L) {
        if (info == 0L) {
          info <- -9L
          nfgv <- nfgv - 1L
          ifun <- ifun - 1L
          iback <- iback - 1L
        }
        task <- "ABNORMAL"
        iter <- iter + 1L
        break
      } else {
        if (info == 0L) nfgv <- nfgv - 1L
        reset_memory()
        next
      }
    }
    iter <- iter + 1L
    sbgnrm <- projgr()
    n_iterations <- n_iterations + 1L
    if (n_iterations >= maxiter) {
      task <- "STOP"
      break
    } else if (sf$nfev() > maxfun) {
      task <- "STOP"
      break
    }
    if (sbgnrm <= pgtol) {
      task <- "CONVERGENCE"
      break
    }
    ddum <- max(abs(fold), abs(f), 1.0)
    if ((fold - f) <= tol * ddum) {
      task <- "CONVERGENCE"
      if (iback >= 10L) info <- -5L
      break
    }
    r <- g - r
    rr <- .sp_norm2(r)^2
    if (stp == 1.0) {
      dr <- gd - gdold
      ddum <- -gdold
    } else {
      dr <- (gd - gdold) * stp
      d <- stp * d
      ddum <- -gdold * stp
    }
    if (dr <= epsmach * ddum) {
      nskip <- nskip + 1L
      updatd <- FALSE
      next
    }
    updatd <- TRUE
    iupdat <- iupdat + 1L
    matupd(rr, dr)
    formt()
    if (info != 0L) {
      reset_memory()
      next
    }
  }
  warnflag <- if (task == "CONVERGENCE") 0L else if (sf$nfev() > maxfun || n_iterations >= maxiter) 1L else 2L
  msg <- switch(as.character(warnflag), "0" = "CONVERGENCE", "1" = "STOP: TOTAL NO. OF ITERATIONS OR FUNCTION EVALUATIONS REACHED LIMIT",
                "ABNORMAL: the line search failed")
  list(fun = f, jac = g, nfev = sf$nfev(), njev = sf$ngev, nit = n_iterations, status = warnflag, message = msg,
       x = x, success = warnflag == 0L)
}
