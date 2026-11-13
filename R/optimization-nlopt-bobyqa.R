# Powell's BOBYQA bound constrained derivative-free trust region method as
# translated to C for NLopt 2.11.0 by S. G. Johnson
# (src/algs/bobyqa/bobyqa.c), ported to R: LN_BOBYQA.
#
# The NLopt library that NLopt.jl loads (NLopt_jll) rounds every product
# and sum separately (no fused multiply add), so plain R arithmetic in the
# order of the C expressions reproduces it. Every sum is accumulated in the
# loop order of the C code.
#
# Copyright (c) 2009 M. J. D. Powell (mjdp@cam.ac.uk); C translation and
# NLopt modifications Copyright (c) 2009 Massachusetts Institute of
# Technology, released under the MIT license (see
# R/optimization-nlopt-util.R).

.nbq_min2 <- function(a, b) if (a <= b) a else b
.nbq_max2 <- function(a, b) if (a >= b) a else b

# update_: the updating of bmat and zmat when xpt[knew, ] is replaced;
# returns list(bmat, zmat, vlag)
.nbq_update <- function(n, npt, bmat, zmat, vlag, beta, denom, knew) {
  nptm <- npt - n - 1L
  ztest <- 0
  for (k in seq_len(npt)) for (j in seq_len(nptm)) ztest <- .nbq_max2(ztest, abs(zmat[k, j]))
  ztest <- ztest * 1e-20
  for (j in seq_len(nptm)[-1L]) {
    if (abs(zmat[knew, j]) > ztest) {
      temp <- sqrt(zmat[knew, 1] * zmat[knew, 1] + zmat[knew, j] * zmat[knew, j])
      tempa <- zmat[knew, 1] / temp
      tempb <- zmat[knew, j] / temp
      t <- tempa * zmat[, 1] + tempb * zmat[, j]
      zmat[, j] <- tempa * zmat[, j] - tempb * zmat[, 1]
      zmat[, 1] <- t
    }
    zmat[knew, j] <- 0
  }
  ip <- seq_len(npt)
  w <- numeric(npt + n)
  w[ip] <- zmat[knew, 1] * zmat[, 1]
  alpha <- w[knew]
  tau <- vlag[knew]
  vlag[knew] <- vlag[knew] - 1
  temp <- sqrt(denom)
  tempb <- zmat[knew, 1] / temp
  tempa <- tau / temp
  zmat[, 1] <- tempa * zmat[, 1] - tempb * vlag[ip]
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
  list(bmat = bmat, zmat = zmat, vlag = vlag)
}

# the evaluation point of xbase + v within [xl, xu], exactly on a bound
# when v is on the shifted bound sl or su
.nbq_point <- function(xl, xu, xbase, v, sl, su) {
  x <- numeric(length(v))
  for (i in seq_along(v)) {
    x[i] <- .nbq_min2(.nbq_max2(xl[i], xbase[i] + v[i]), xu[i])
    if (v[i] == sl[i]) x[i] <- xl[i]
    if (v[i] == su[i]) x[i] <- xu[i]
  }
  x
}

# rescue_: a new set of interpolation points when the denominators of the
# updating formula are too small; s is the state list of bobyqb_ (xbase,
# xpt, fval, xopt, gopt, hq, pq, bmat, zmat, sl, su, kopt, vlag), returned
# updated with the result code ret
.nbq_rescue <- function(n, npt, xl, xu, stop, calfun, s, delta) {
  ndim <- npt + n
  xbase <- s$xbase; xpt <- s$xpt; fval <- s$fval; xopt <- s$xopt; gopt <- s$gopt; hq <- s$hq
  pq <- s$pq; bmat <- s$bmat; zmat <- s$zmat; sl <- s$sl; su <- s$su; kopt <- s$kopt; vlag <- s$vlag
  out <- function(ret) list(ret = ret, xbase = xbase, xpt = xpt, fval = fval, xopt = xopt, gopt = gopt,
                            hq = hq, pq = pq, bmat = bmat, zmat = zmat, sl = sl, su = su, kopt = kopt, vlag = vlag)
  np <- n + 1L
  sfrac <- 0.5 / np
  nptm <- npt - np
  w <- numeric(ndim + npt)
  ptsaux <- matrix(0, 2L, n)
  ptsid <- numeric(npt)
  xp <- 0; xq <- 0; ihp <- 0L; beta <- 0; denom <- 0
  sumpq <- 0; winc <- 0
  for (k in seq_len(npt)) {
    distsq <- 0
    for (j in seq_len(n)) {
      xpt[k, j] <- xpt[k, j] - xopt[j]
      distsq <- distsq + xpt[k, j] * xpt[k, j]
    }
    sumpq <- sumpq + pq[k]
    w[ndim + k] <- distsq
    winc <- .nbq_max2(winc, distsq)
    zmat[k, ] <- 0
  }
  ih <- 0L
  for (j in seq_len(n)) {
    w[j] <- 0.5 * sumpq * xopt[j]
    for (k in seq_len(npt)) w[j] <- w[j] + pq[k] * xpt[k, j]
    for (i in seq_len(j)) {
      ih <- ih + 1L
      hq[ih] <- hq[ih] + w[i] * xopt[j] + w[j] * xopt[i]
    }
  }
  for (j in seq_len(n)) {
    xbase[j] <- xbase[j] + xopt[j]
    sl[j] <- sl[j] - xopt[j]
    su[j] <- su[j] - xopt[j]
    xopt[j] <- 0
    ptsaux[1, j] <- .nbq_min2(delta, su[j])
    ptsaux[2, j] <- .nbq_max2(-delta, sl[j])
    if (ptsaux[1, j] + ptsaux[2, j] < 0) {
      temp <- ptsaux[1, j]
      ptsaux[1, j] <- ptsaux[2, j]
      ptsaux[2, j] <- temp
    }
    if (abs(ptsaux[2, j]) < 0.5 * abs(ptsaux[1, j])) ptsaux[2, j] <- 0.5 * ptsaux[1, j]
    bmat[, j] <- 0
  }
  fbase <- fval[kopt]
  ptsid[1] <- sfrac
  for (j in seq_len(n)) {
    jp <- j + 1L
    jpn <- jp + n
    ptsid[jp] <- j + sfrac
    if (jpn <= npt) {
      ptsid[jpn] <- j / np + sfrac
      temp <- 1 / (ptsaux[1, j] - ptsaux[2, j])
      bmat[jp, j] <- -temp + 1 / ptsaux[1, j]
      bmat[jpn, j] <- temp + 1 / ptsaux[2, j]
      bmat[1, j] <- -bmat[jp, j] - bmat[jpn, j]
      zmat[1, j] <- sqrt(2) / abs(ptsaux[1, j] * ptsaux[2, j])
      zmat[jp, j] <- zmat[1, j] * ptsaux[2, j] * temp
      zmat[jpn, j] <- -zmat[1, j] * ptsaux[1, j] * temp
    } else {
      bmat[1, j] <- -1 / ptsaux[1, j]
      bmat[jp, j] <- 1 / ptsaux[1, j]
      bmat[j + npt, j] <- -0.5 * (ptsaux[1, j] * ptsaux[1, j])
    }
  }
  if (npt >= n + np) for (k in seq.int(2L * np, npt)) {
    iw <- as.integer(((k - np) - 0.5) / n)
    ip <- k - np - iw * n
    iq <- ip + iw
    if (iq > n) iq <- iq - n
    ptsid[k] <- ip + iq / np + sfrac
    temp <- 1 / (ptsaux[1, ip] * ptsaux[1, iq])
    zmat[1, k - np] <- temp
    zmat[ip + 1L, k - np] <- -temp
    zmat[iq + 1L, k - np] <- -temp
    zmat[k, k - np] <- temp
  }
  nrem <- npt
  kold <- 1L
  knew <- kopt
  lab <- "L80"
  repeat {
    if (lab == "L80") {
      temp <- bmat[kold, ]; bmat[kold, ] <- bmat[knew, ]; bmat[knew, ] <- temp
      temp <- zmat[kold, ]; zmat[kold, ] <- zmat[knew, ]; zmat[knew, ] <- temp
      ptsid[kold] <- ptsid[knew]
      ptsid[knew] <- 0
      w[ndim + knew] <- 0
      nrem <- nrem - 1L
      if (knew != kopt) {
        temp <- vlag[kold]; vlag[kold] <- vlag[knew]; vlag[knew] <- temp
        up <- .nbq_update(n, npt, bmat, zmat, vlag, beta, denom, knew)
        bmat <- up$bmat; zmat <- up$zmat; vlag <- up$vlag
        if (nrem == 0L) return(out("SUCCESS"))
        w[ndim + seq_len(npt)] <- abs(w[ndim + seq_len(npt)])
      }
      lab <- "L120"
    } else if (lab == "L120") {
      dsqmin <- 0
      for (k in seq_len(npt)) if (w[ndim + k] > 0) {
        if (dsqmin == 0 || w[ndim + k] < dsqmin) {
          knew <- k
          dsqmin <- w[ndim + k]
        }
      }
      if (dsqmin == 0) break
      w[npt + seq_len(n)] <- xpt[knew, ]
      for (k in seq_len(npt)) {
        sum <- 0
        if (k == kopt) {
        } else if (ptsid[k] == 0) {
          for (j in seq_len(n)) sum <- sum + w[npt + j] * xpt[k, j]
        } else {
          ip <- as.integer(ptsid[k])
          if (ip > 0L) sum <- w[npt + ip] * ptsaux[1, ip]
          iq <- as.integer(np * ptsid[k] - ip * np)
          if (iq > 0L) {
            iw <- if (ip == 0L) 2L else 1L
            sum <- sum + w[npt + iq] * ptsaux[iw, iq]
          }
        }
        w[k] <- 0.5 * sum * sum
      }
      for (k in seq_len(npt)) {
        sum <- 0
        for (j in seq_len(n)) sum <- sum + bmat[k, j] * w[npt + j]
        vlag[k] <- sum
      }
      beta <- 0
      for (j in seq_len(nptm)) {
        sum <- 0
        for (k in seq_len(npt)) sum <- sum + zmat[k, j] * w[k]
        beta <- beta - sum * sum
        vlag[seq_len(npt)] <- vlag[seq_len(npt)] + sum * zmat[, j]
      }
      bsum <- 0; distsq <- 0
      for (j in seq_len(n)) {
        sum <- 0
        for (k in seq_len(npt)) sum <- sum + bmat[k, j] * w[k]
        jp <- j + npt
        bsum <- bsum + sum * w[jp]
        for (ip in seq.int(npt + 1L, ndim)) sum <- sum + bmat[ip, j] * w[ip]
        bsum <- bsum + sum * w[jp]
        vlag[jp] <- sum
        distsq <- distsq + xpt[knew, j] * xpt[knew, j]
      }
      beta <- 0.5 * distsq * distsq + beta - bsum
      vlag[kopt] <- vlag[kopt] + 1
      denom <- 0; vlmxsq <- 0
      for (k in seq_len(npt)) {
        if (ptsid[k] != 0) {
          hdiag <- 0
          for (j in seq_len(nptm)) hdiag <- hdiag + zmat[k, j] * zmat[k, j]
          den <- beta * hdiag + vlag[k] * vlag[k]
          if (den > denom) {
            kold <- k
            denom <- den
          }
        }
        vlmxsq <- .nbq_max2(vlmxsq, vlag[k] * vlag[k])
      }
      if (denom <= vlmxsq * 0.01) {
        w[ndim + knew] <- -w[ndim + knew] - winc
        lab <- "L120"
      } else {
        lab <- "L80"
      }
    }
  }
  # L260
  for (kpt in seq_len(npt)) {
    if (ptsid[kpt] == 0) next
    if (.nlopt_stop_forced(stop)) return(out("FORCED_STOP"))
    else if (.nlopt_stop_evals(stop)) return(out("MAXEVAL_REACHED"))
    else if (.nlopt_stop_time(stop)) return(out("MAXTIME_REACHED"))
    ih <- 0L
    for (j in seq_len(n)) {
      w[j] <- xpt[kpt, j]
      xpt[kpt, j] <- 0
      temp <- pq[kpt] * w[j]
      for (i in seq_len(j)) {
        ih <- ih + 1L
        hq[ih] <- hq[ih] + temp * w[i]
      }
    }
    pq[kpt] <- 0
    ip <- as.integer(ptsid[kpt])
    iq <- as.integer(np * ptsid[kpt] - ip * np)
    if (ip > 0L) {
      xp <- ptsaux[1, ip]
      xpt[kpt, ip] <- xp
    }
    if (iq > 0L) {
      xq <- ptsaux[1, iq]
      if (ip == 0L) xq <- ptsaux[2, iq]
      xpt[kpt, iq] <- xq
    }
    vquad <- fbase
    if (ip > 0L) {
      ihp <- (ip + ip * ip) %/% 2L
      vquad <- vquad + xp * (gopt[ip] + 0.5 * xp * hq[ihp])
    }
    if (iq > 0L) {
      ihq <- (iq + iq * iq) %/% 2L
      vquad <- vquad + xq * (gopt[iq] + 0.5 * xq * hq[ihq])
      if (ip > 0L) {
        iw <- .nbq_max2(ihp, ihq) - abs(ip - iq)
        vquad <- vquad + xp * xq * hq[iw]
      }
    }
    for (k in seq_len(npt)) {
      temp <- 0
      if (ip > 0L) temp <- temp + xp * xpt[k, ip]
      if (iq > 0L) temp <- temp + xq * xpt[k, iq]
      vquad <- vquad + 0.5 * pq[k] * temp * temp
    }
    w[seq_len(n)] <- .nbq_point(xl, xu, xbase, xpt[kpt, ], sl, su)
    stop$nevals <- stop$nevals + 1L
    f <- calfun(w[seq_len(n)])
    fval[kpt] <- f
    if (f < fval[kopt]) kopt <- kpt
    if (.nlopt_stop_forced(stop)) return(out("FORCED_STOP"))
    else if (f < stop$minf_max) return(out("STOPVAL_REACHED"))
    else if (.nlopt_stop_evals(stop)) return(out("MAXEVAL_REACHED"))
    else if (.nlopt_stop_time(stop)) return(out("MAXTIME_REACHED"))
    diff <- f - vquad
    gopt <- gopt + diff * bmat[kpt, ]
    for (k in seq_len(npt)) {
      sum <- 0
      for (j in seq_len(nptm)) sum <- sum + zmat[k, j] * zmat[kpt, j]
      temp <- diff * sum
      if (ptsid[k] == 0) {
        pq[k] <- pq[k] + temp
      } else {
        ip <- as.integer(ptsid[k])
        iq <- as.integer(np * ptsid[k] - ip * np)
        ihq <- (iq * iq + iq) %/% 2L
        if (ip == 0L) {
          hq[ihq] <- hq[ihq] + temp * (ptsaux[2, iq] * ptsaux[2, iq])
        } else {
          ihp <- (ip * ip + ip) %/% 2L
          hq[ihp] <- hq[ihp] + temp * (ptsaux[1, ip] * ptsaux[1, ip])
          if (iq > 0L) {
            hq[ihq] <- hq[ihq] + temp * (ptsaux[1, iq] * ptsaux[1, iq])
            iw <- .nbq_max2(ihp, ihq) - abs(iq - ip)
            hq[iw] <- hq[iw] + temp * ptsaux[1, ip] * ptsaux[1, iq]
          }
        }
      }
    }
    ptsid[kpt] <- 0
  }
  out("SUCCESS")
}

# altmov_: an alternative to the trust region step that makes the KNEW-th
# Lagrange function large; returns list(xnew, xalt, alpha, cauchy)
.nbq_altmov <- function(n, npt, xpt, xopt, bmat, zmat, sl, su, kopt, knew, adelt) {
  const <- 1 + sqrt(2)
  hcol <- numeric(npt)
  for (j in seq_len(npt - n - 1L)) {
    temp <- zmat[knew, j]
    hcol <- hcol + temp * zmat[, j]
  }
  alpha <- hcol[knew]
  ha <- 0.5 * alpha
  glag <- bmat[knew, ]
  for (k in seq_len(npt)) {
    temp <- 0
    for (j in seq_len(n)) temp <- temp + xpt[k, j] * xopt[j]
    temp <- hcol[k] * temp
    glag <- glag + temp * xpt[k, ]
  }
  presav <- 0; ksav <- 0L; stpsav <- 0; ibdsav <- 0L; step <- 0
  for (k in seq_len(npt)) {
    if (k == kopt) next
    dderiv <- 0; distsq <- 0
    for (i in seq_len(n)) {
      temp <- xpt[k, i] - xopt[i]
      dderiv <- dderiv + glag[i] * temp
      distsq <- distsq + temp * temp
    }
    subd <- adelt / sqrt(distsq)
    slbd <- -subd
    ilbd <- 0L; iubd <- 0L
    sumin <- .nbq_min2(1, subd)
    for (i in seq_len(n)) {
      temp <- xpt[k, i] - xopt[i]
      if (temp > 0) {
        if (slbd * temp < sl[i] - xopt[i]) {
          slbd <- (sl[i] - xopt[i]) / temp
          ilbd <- -i
        }
        if (subd * temp > su[i] - xopt[i]) {
          subd <- .nbq_max2(sumin, (su[i] - xopt[i]) / temp)
          iubd <- i
        }
      } else if (temp < 0) {
        if (slbd * temp > su[i] - xopt[i]) {
          slbd <- (su[i] - xopt[i]) / temp
          ilbd <- i
        }
        if (subd * temp < sl[i] - xopt[i]) {
          subd <- .nbq_max2(sumin, (sl[i] - xopt[i]) / temp)
          iubd <- -i
        }
      }
    }
    if (k == knew) {
      diff <- dderiv - 1
      step <- slbd
      vlag <- slbd * (dderiv - slbd * diff)
      isbd <- ilbd
      temp <- subd * (dderiv - subd * diff)
      if (abs(temp) > abs(vlag)) {
        step <- subd; vlag <- temp; isbd <- iubd
      }
      tempd <- 0.5 * dderiv
      tempa <- tempd - diff * slbd
      tempb <- tempd - diff * subd
      if (tempa * tempb < 0) {
        temp <- tempd * tempd / diff
        if (abs(temp) > abs(vlag)) {
          step <- tempd / diff; vlag <- temp; isbd <- 0L
        }
      }
    } else {
      step <- slbd
      vlag <- slbd * (1 - slbd)
      isbd <- ilbd
      temp <- subd * (1 - subd)
      if (abs(temp) > abs(vlag)) {
        step <- subd; vlag <- temp; isbd <- iubd
      }
      if (subd > 0.5 && abs(vlag) < 0.25) {
        step <- 0.5; vlag <- 0.25; isbd <- 0L
      }
      vlag <- vlag * dderiv
    }
    temp <- step * (1 - step) * distsq
    predsq <- vlag * vlag * (vlag * vlag + ha * temp * temp)
    if (predsq > presav) {
      presav <- predsq; ksav <- k; stpsav <- step; ibdsav <- isbd
    }
  }
  xnew <- numeric(n)
  for (i in seq_len(n)) {
    temp <- xopt[i] + stpsav * (xpt[ksav, i] - xopt[i])
    xnew[i] <- .nbq_max2(sl[i], .nbq_min2(su[i], temp))
  }
  if (ibdsav < 0L) xnew[-ibdsav] <- sl[-ibdsav]
  if (ibdsav > 0L) xnew[ibdsav] <- su[ibdsav]
  bigstp <- adelt + adelt
  iflag <- 0L
  xalt <- numeric(n); w <- numeric(2L * n); cauchy <- 0; csave <- 0
  repeat {
    # L100
    wfixsq <- 0; ggfree <- 0
    for (i in seq_len(n)) {
      w[i] <- 0
      tempa <- .nbq_min2(xopt[i] - sl[i], glag[i])
      tempb <- .nbq_max2(xopt[i] - su[i], glag[i])
      if (tempa > 0 || tempb < 0) {
        w[i] <- bigstp
        ggfree <- ggfree + glag[i] * glag[i]
      }
    }
    if (ggfree == 0) {
      cauchy <- 0
      break
    }
    repeat {
      # L120
      temp <- adelt * adelt - wfixsq
      if (!(temp > 0)) break
      wsqsav <- wfixsq
      step <- sqrt(temp / ggfree)
      ggfree <- 0
      for (i in seq_len(n)) if (w[i] == bigstp) {
        temp <- xopt[i] - step * glag[i]
        if (temp <= sl[i]) {
          w[i] <- sl[i] - xopt[i]
          wfixsq <- wfixsq + w[i] * w[i]
        } else if (temp >= su[i]) {
          w[i] <- su[i] - xopt[i]
          wfixsq <- wfixsq + w[i] * w[i]
        } else {
          ggfree <- ggfree + glag[i] * glag[i]
        }
      }
      if (!(wfixsq > wsqsav && ggfree > 0)) break
    }
    gw <- 0
    for (i in seq_len(n)) {
      if (w[i] == bigstp) {
        w[i] <- -step * glag[i]
        xalt[i] <- .nbq_max2(sl[i], .nbq_min2(su[i], xopt[i] + w[i]))
      } else if (w[i] == 0) {
        xalt[i] <- xopt[i]
      } else if (glag[i] > 0) {
        xalt[i] <- sl[i]
      } else {
        xalt[i] <- su[i]
      }
      gw <- gw + glag[i] * w[i]
    }
    curv <- 0
    for (k in seq_len(npt)) {
      temp <- 0
      for (j in seq_len(n)) temp <- temp + xpt[k, j] * w[j]
      curv <- curv + hcol[k] * temp * temp
    }
    if (iflag == 1L) curv <- -curv
    if (curv > -gw && curv < -const * gw) {
      scale <- -gw / curv
      for (i in seq_len(n)) {
        temp <- xopt[i] + scale * w[i]
        xalt[i] <- .nbq_max2(sl[i], .nbq_min2(su[i], temp))
      }
      t <- 0.5 * gw * scale
      cauchy <- t * t
    } else {
      t <- gw + 0.5 * curv
      cauchy <- t * t
    }
    if (iflag == 0L) {
      glag <- -glag
      w[n + seq_len(n)] <- xalt
      csave <- cauchy
      iflag <- 1L
      next
    }
    if (csave > cauchy) {
      xalt <- w[n + seq_len(n)]
      cauchy <- csave
    }
    break
  }
  list(xnew = xnew, xalt = xalt, alpha = alpha, cauchy = cauchy)
}

# trsbox_: an approximate minimizer of the quadratic model within the trust
# region and the bounds; returns list(xnew, d, dsq, crvmin, gnew)
.nbq_trsbox <- function(n, npt, xpt, xopt, gopt, hq, pq, sl, su, delta) {
  xbdi <- numeric(n); d <- numeric(n); gnew <- gopt; s <- numeric(n); hs <- numeric(n); hred <- numeric(n)
  iterc <- 0L; nact <- 0L
  iact <- 0L; xsav <- 0; angbd <- 0; dredg <- 0; sredg <- 0; ggsav <- 0; dredsq <- 0; gredsq <- 0
  itcsav <- 0L; rdnext <- 0; rdprev <- 0; itermax <- 0L; beta <- 0; stepsq <- 0; angt <- 0
  for (i in seq_len(n)) {
    if (xopt[i] <= sl[i]) {
      if (gopt[i] >= 0) xbdi[i] <- -1
    } else if (xopt[i] >= su[i]) {
      if (gopt[i] <= 0) xbdi[i] <- 1
    }
    if (xbdi[i] != 0) nact <- nact + 1L
  }
  delsq <- delta * delta
  qred <- 0
  crvmin <- -1
  lab <- "L20"
  repeat {
    if (lab == "L20") {
      beta <- 0
      lab <- "L30"
    } else if (lab == "L30") {
      stepsq <- 0
      for (i in seq_len(n)) {
        if (xbdi[i] != 0) s[i] <- 0
        else if (beta == 0) s[i] <- -gnew[i]
        else s[i] <- beta * s[i] - gnew[i]
        stepsq <- stepsq + s[i] * s[i]
      }
      if (stepsq == 0) break
      if (beta == 0) {
        gredsq <- stepsq
        itermax <- iterc + n - nact
      }
      lab <- if (gredsq * delsq <= qred * 1e-4 * qred) "L190" else "L210"
      if (lab == "L190") break
    } else if (lab == "L50") {
      resid <- delsq; ds <- 0; shs <- 0
      for (i in seq_len(n)) if (xbdi[i] == 0) {
        resid <- resid - d[i] * d[i]
        ds <- ds + s[i] * d[i]
        shs <- shs + s[i] * hs[i]
      }
      if (resid <= 0) { lab <- "L90"; next }
      temp <- sqrt(stepsq * resid + ds * ds)
      blen <- if (ds < 0) (temp - ds) / stepsq else resid / (temp + ds)
      stplen <- blen
      if (shs > 0) stplen <- .nbq_min2(blen, gredsq / shs)
      iact <- 0L
      for (i in seq_len(n)) if (s[i] != 0) {
        xsum <- xopt[i] + d[i]
        temp <- if (s[i] > 0) (su[i] - xsum) / s[i] else (sl[i] - xsum) / s[i]
        if (temp < stplen) {
          stplen <- temp
          iact <- i
        }
      }
      sdec <- 0
      if (stplen > 0) {
        iterc <- iterc + 1L
        temp <- shs / stepsq
        if (iact == 0L && temp > 0) {
          crvmin <- .nbq_min2(crvmin, temp)
          if (crvmin == -1) crvmin <- temp
        }
        ggsav <- gredsq
        gredsq <- 0
        for (i in seq_len(n)) {
          gnew[i] <- gnew[i] + stplen * hs[i]
          if (xbdi[i] == 0) gredsq <- gredsq + gnew[i] * gnew[i]
          d[i] <- d[i] + stplen * s[i]
        }
        sdec <- .nbq_max2(stplen * (ggsav - 0.5 * stplen * shs), 0)
        qred <- qred + sdec
      }
      if (iact > 0L) {
        nact <- nact + 1L
        xbdi[iact] <- 1
        if (s[iact] < 0) xbdi[iact] <- -1
        delsq <- delsq - d[iact] * d[iact]
        lab <- if (delsq <= 0) "L90" else "L20"
        next
      }
      if (stplen < blen) {
        if (iterc == itermax) break
        if (sdec <= qred * 0.01) break
        beta <- gredsq / ggsav
        lab <- "L30"
        next
      }
      lab <- "L90"
    } else if (lab == "L90") {
      crvmin <- 0
      lab <- "L100"
    } else if (lab == "L100") {
      if (nact >= n - 1L) break
      dredsq <- 0; dredg <- 0; gredsq <- 0
      for (i in seq_len(n)) {
        if (xbdi[i] == 0) {
          dredsq <- dredsq + d[i] * d[i]
          dredg <- dredg + d[i] * gnew[i]
          gredsq <- gredsq + gnew[i] * gnew[i]
          s[i] <- d[i]
        } else {
          s[i] <- 0
        }
      }
      itcsav <- iterc
      lab <- "L210"
    } else if (lab == "L120") {
      iterc <- iterc + 1L
      temp <- gredsq * dredsq - dredg * dredg
      if (temp <= qred * 1e-4 * qred) break
      temp <- sqrt(temp)
      for (i in seq_len(n)) s[i] <- if (xbdi[i] == 0) (dredg * d[i] - dredsq * gnew[i]) / temp else 0
      sredg <- -temp
      angbd <- 1
      iact <- 0L
      lab <- "L210"
      for (i in seq_len(n)) if (xbdi[i] == 0) {
        tempa <- xopt[i] + d[i] - sl[i]
        tempb <- su[i] - xopt[i] - d[i]
        if (tempa <= 0) {
          nact <- nact + 1L
          xbdi[i] <- -1
          lab <- "L100"
          break
        } else if (tempb <= 0) {
          nact <- nact + 1L
          xbdi[i] <- 1
          lab <- "L100"
          break
        }
        ssq <- d[i] * d[i] + s[i] * s[i]
        t <- xopt[i] - sl[i]
        temp <- ssq - t * t
        if (temp > 0) {
          temp <- sqrt(temp) - s[i]
          if (angbd * temp > tempa) {
            angbd <- tempa / temp
            iact <- i
            xsav <- -1
          }
        }
        t <- su[i] - xopt[i]
        temp <- ssq - t * t
        if (temp > 0) {
          temp <- sqrt(temp) + s[i]
          if (angbd * temp > tempb) {
            angbd <- tempb / temp
            iact <- i
            xsav <- 1
          }
        }
      }
    } else if (lab == "L150") {
      shs <- 0; dhs <- 0; dhd <- 0
      for (i in seq_len(n)) if (xbdi[i] == 0) {
        shs <- shs + s[i] * hs[i]
        dhs <- dhs + d[i] * hs[i]
        dhd <- dhd + d[i] * hred[i]
      }
      redmax <- 0; isav <- 0L; redsav <- 0
      iu <- as.integer(angbd * 17 + 3.1)
      for (i in seq_len(iu)) {
        angt <- angbd * i / iu
        sth <- (angt + angt) / (1 + angt * angt)
        temp <- shs + angt * (angt * dhd - dhs - dhs)
        rednew <- sth * (angt * dredg - sredg - 0.5 * sth * temp)
        if (rednew > redmax) {
          redmax <- rednew
          isav <- i
          rdprev <- redsav
        } else if (i == isav + 1L) {
          rdnext <- rednew
        }
        redsav <- rednew
      }
      if (isav == 0L) break
      if (isav < iu) {
        temp <- (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext)
        angt <- angbd * (isav + 0.5 * temp) / iu
      }
      cth <- (1 - angt * angt) / (1 + angt * angt)
      sth <- (angt + angt) / (1 + angt * angt)
      temp <- shs + angt * (angt * dhd - dhs - dhs)
      sdec <- sth * (angt * dredg - sredg - 0.5 * sth * temp)
      if (sdec <= 0) break
      dredg <- 0; gredsq <- 0
      for (i in seq_len(n)) {
        gnew[i] <- gnew[i] + (cth - 1) * hred[i] + sth * hs[i]
        if (xbdi[i] == 0) {
          d[i] <- cth * d[i] + sth * s[i]
          dredg <- dredg + d[i] * gnew[i]
          gredsq <- gredsq + gnew[i] * gnew[i]
        }
        hred[i] <- cth * hred[i] + sth * hs[i]
      }
      qred <- qred + sdec
      if (iact > 0L && isav == iu) {
        nact <- nact + 1L
        xbdi[iact] <- xsav
        lab <- "L100"
        next
      }
      if (sdec > qred * 0.01) { lab <- "L120"; next }
      break
    } else {
      # L210: hs = H s
      ih <- 0L
      for (j in seq_len(n)) {
        hs[j] <- 0
        for (i in seq_len(j)) {
          ih <- ih + 1L
          if (i < j) hs[j] <- hs[j] + hq[ih] * s[i]
          hs[i] <- hs[i] + hq[ih] * s[j]
        }
      }
      for (k in seq_len(npt)) if (pq[k] != 0) {
        temp <- 0
        for (j in seq_len(n)) temp <- temp + xpt[k, j] * s[j]
        temp <- temp * pq[k]
        hs <- hs + temp * xpt[k, ]
      }
      if (crvmin != 0) {
        lab <- "L50"
      } else if (iterc > itcsav) {
        lab <- "L150"
      } else {
        hred <- hs
        lab <- "L120"
      }
    }
  }
  # L190
  dsq <- 0
  xnew <- numeric(n)
  for (i in seq_len(n)) {
    xnew[i] <- .nbq_max2(.nbq_min2(xopt[i] + d[i], su[i]), sl[i])
    if (xbdi[i] == -1) xnew[i] <- sl[i]
    if (xbdi[i] == 1) xnew[i] <- su[i]
    d[i] <- xnew[i] - xopt[i]
    dsq <- dsq + d[i] * d[i]
  }
  list(xnew = xnew, d = d, dsq = dsq, crvmin = crvmin, gnew = gnew)
}

# prelim_: the initial interpolation points and quadratic model; returns
# the state list of bobyqb_ with ret, x and kopt
.nbq_prelim <- function(n, npt, x, xl, xu, rhobeg, stop, calfun, sl, su) {
  ndim <- npt + n
  np <- n + 1L
  rhosq <- rhobeg * rhobeg
  recip <- 1 / rhosq
  xbase <- x
  xpt <- matrix(0, npt, n); bmat <- matrix(0, ndim, n); zmat <- matrix(0, npt, npt - np)
  hq <- numeric((n * np) %/% 2L); pq <- numeric(npt); fval <- numeric(npt); gopt <- numeric(n)
  kopt <- 1L
  fbeg <- 0; stepa <- 0; stepb <- 0; ipt <- 0L; jpt <- 0L
  out <- function(ret) list(ret = ret, x = x, xbase = xbase, xpt = xpt, fval = fval, gopt = gopt, hq = hq,
                            pq = pq, bmat = bmat, zmat = zmat, kopt = kopt)
  nf <- 0L
  repeat {
    nfm <- nf
    nfx <- nf - n
    nf <- nf + 1L
    if (nfm <= 2L * n) {
      if (nfm >= 1L && nfm <= n) {
        stepa <- rhobeg
        if (su[nfm] == 0) stepa <- -stepa
        xpt[nf, nfm] <- stepa
      } else if (nfm > n) {
        stepa <- xpt[nf - n, nfx]
        stepb <- -rhobeg
        if (sl[nfx] == 0) stepb <- .nbq_min2(2 * rhobeg, su[nfx])
        if (su[nfx] == 0) stepb <- .nbq_max2(-2 * rhobeg, sl[nfx])
        xpt[nf, nfx] <- stepb
      }
    } else {
      itemp <- (nfm - np) %/% n
      jpt <- nfm - itemp * n - n
      ipt <- jpt + itemp
      if (ipt > n) {
        itemp <- jpt
        jpt <- ipt - n
        ipt <- itemp
      }
      xpt[nf, ipt] <- xpt[ipt + 1L, ipt]
      xpt[nf, jpt] <- xpt[jpt + 1L, jpt]
    }
    x <- .nbq_point(xl, xu, xbase, xpt[nf, ], sl, su)
    stop$nevals <- stop$nevals + 1L
    f <- calfun(x)
    fval[nf] <- f
    if (nf == 1L) {
      fbeg <- f
      kopt <- 1L
    } else if (f < fval[kopt]) {
      kopt <- nf
    }
    if (nf <= 2L * n + 1L) {
      if (nf >= 2L && nf <= n + 1L) {
        gopt[nfm] <- (f - fbeg) / stepa
        if (npt < nf + n) {
          bmat[1, nfm] <- -1 / stepa
          bmat[nf, nfm] <- 1 / stepa
          bmat[npt + nfm, nfm] <- -0.5 * rhosq
        }
      } else if (nf >= n + 2L) {
        ih <- (nfx * (nfx + 1L)) %/% 2L
        temp <- (f - fbeg) / stepb
        diff <- stepb - stepa
        hq[ih] <- 2 * (temp - gopt[nfx]) / diff
        gopt[nfx] <- (gopt[nfx] * stepb - temp * stepa) / diff
        if (stepa * stepb < 0 && f < fval[nf - n]) {
          fval[nf] <- fval[nf - n]
          fval[nf - n] <- f
          if (kopt == nf) kopt <- nf - n
          xpt[nf - n, nfx] <- stepb
          xpt[nf, nfx] <- stepa
        }
        bmat[1, nfx] <- -(stepa + stepb) / (stepa * stepb)
        bmat[nf, nfx] <- -0.5 / xpt[nf - n, nfx]
        bmat[nf - n, nfx] <- -bmat[1, nfx] - bmat[nf, nfx]
        zmat[1, nfx] <- sqrt(2) / (stepa * stepb)
        zmat[nf, nfx] <- sqrt(0.5) / rhosq
        zmat[nf - n, nfx] <- -zmat[1, nfx] - zmat[nf, nfx]
      }
    } else {
      ih <- (ipt * (ipt - 1L)) %/% 2L + jpt
      zmat[1, nfx] <- recip
      zmat[nf, nfx] <- recip
      zmat[ipt + 1L, nfx] <- -recip
      zmat[jpt + 1L, nfx] <- -recip
      temp <- xpt[nf, ipt] * xpt[nf, jpt]
      hq[ih] <- (fbeg - fval[ipt + 1L] - fval[jpt + 1L] + f) / temp
    }
    if (.nlopt_stop_forced(stop)) return(out("FORCED_STOP"))
    else if (f < stop$minf_max) return(out("STOPVAL_REACHED"))
    else if (.nlopt_stop_evals(stop)) return(out("MAXEVAL_REACHED"))
    else if (.nlopt_stop_time(stop)) return(out("MAXTIME_REACHED"))
    if (nf >= npt) return(out("SUCCESS"))
  }
}

# bobyqb_: the main iteration of BOBYQA on the scaled variables; returns
# list(ret, x, minf)
.nbq_bobyqb <- function(n, npt, x, xl, xu, rhobeg, rhoend, stop, calfun, sl, su) {
  ndim <- npt + n
  np <- n + 1L
  nptm <- npt - np
  nh <- (n * np) %/% 2L
  pr <- .nbq_prelim(n, npt, x, xl, xu, rhobeg, stop, calfun, sl, su)
  x <- pr$x; xbase <- pr$xbase; xpt <- pr$xpt; fval <- pr$fval; gopt <- pr$gopt; hq <- pr$hq; pq <- pr$pq
  bmat <- pr$bmat; zmat <- pr$zmat; kopt <- pr$kopt
  vlag <- numeric(ndim); d <- numeric(n); xnew <- numeric(n); xalt <- numeric(n)
  gnew <- numeric(n); wsum <- numeric(npt)
  xopt <- xpt[kopt, ]
  xoptsq <- 0
  for (i in seq_len(n)) xoptsq <- xoptsq + xopt[i] * xopt[i]
  fsave <- fval[1]
  rc <- "SUCCESS"
  knew <- 0L; kbase <- 1L; rho <- rhobeg; delta <- rho; nresc <- 0L; ntrits <- 0L
  diffa <- 0; diffb <- 0; diffc <- 0; itest <- 0L; nfsav <- 0L
  alpha <- 0; adelt <- 0; ratio <- 0; cauchy <- 0; dsq <- 0; crvmin <- 0; dnorm <- 0
  beta <- 0; denom <- 0; distsq <- 0; f <- 0; fopt <- 0; vquad <- 0; diff <- 0
  lab <- "L20"
  if (pr$ret != "SUCCESS") {
    rc <- pr$ret
    lab <- "L720"
  } else {
    nresc <- stop$nevals
    nfsav <- stop$nevals
  }
  state <- function() list(xbase = xbase, xpt = xpt, fval = fval, xopt = xopt, gopt = gopt, hq = hq, pq = pq,
                           bmat = bmat, zmat = zmat, sl = sl, su = su, kopt = kopt, vlag = vlag)
  repeat {
    switch(lab,
      L20 = {
        if (kopt != kbase) {
          ih <- 0L
          for (j in seq_len(n)) for (i in seq_len(j)) {
            ih <- ih + 1L
            if (i < j) gopt[j] <- gopt[j] + hq[ih] * xopt[i]
            gopt[i] <- gopt[i] + hq[ih] * xopt[j]
          }
          if (stop$nevals > npt) for (k in seq_len(npt)) {
            temp <- 0
            for (j in seq_len(n)) temp <- temp + xpt[k, j] * xopt[j]
            temp <- pq[k] * temp
            gopt <- gopt + temp * xpt[k, ]
          }
        }
        lab <- "L60"
      },
      L60 = {
        tb <- .nbq_trsbox(n, npt, xpt, xopt, gopt, hq, pq, sl, su, delta)
        xnew <- tb$xnew; d <- tb$d; dsq <- tb$dsq; crvmin <- tb$crvmin; gnew <- tb$gnew
        dnorm <- .nbq_min2(delta, sqrt(dsq))
        if (dnorm < 0.5 * rho) {
          ntrits <- -1L
          t <- 10 * rho
          distsq <- t * t
          lab <- "L680"
          if (stop$nevals <= nfsav + 2L) {
            lab <- "L650"
          } else {
            errbig <- .nbq_max2(.nbq_max2(diffa, diffb), diffc)
            frhosq <- rho * 0.125 * rho
            if (crvmin > 0 && errbig > frhosq * crvmin) {
              lab <- "L650"
            } else {
              bdtol <- errbig / rho
              for (j in seq_len(n)) {
                bdtest <- bdtol
                if (xnew[j] == sl[j]) bdtest <- gnew[j]
                if (xnew[j] == su[j]) bdtest <- -gnew[j]
                if (bdtest < bdtol) {
                  curv <- hq[(j + j * j) %/% 2L]
                  for (k in seq_len(npt)) curv <- curv + pq[k] * (xpt[k, j] * xpt[k, j])
                  bdtest <- bdtest + 0.5 * curv * rho
                  if (bdtest < bdtol) {
                    lab <- "L650"
                    break
                  }
                }
              }
            }
          }
        } else {
          ntrits <- ntrits + 1L
          lab <- "L90"
        }
      },
      L90 = {
        if (dsq <= xoptsq * 0.001) {
          fracsq <- xoptsq * 0.25
          sumpq <- 0
          wa <- numeric(n)
          for (k in seq_len(npt)) {
            sumpq <- sumpq + pq[k]
            sum <- -0.5 * xoptsq
            for (i in seq_len(n)) sum <- sum + xpt[k, i] * xopt[i]
            wsum[k] <- sum
            temp <- fracsq - 0.5 * sum
            for (i in seq_len(n)) {
              wa[i] <- bmat[k, i]
              vlag[i] <- sum * xpt[k, i] + temp * xopt[i]
              ip <- npt + i
              for (j in seq_len(i)) bmat[ip, j] <- bmat[ip, j] + wa[i] * vlag[j] + vlag[i] * wa[j]
            }
          }
          for (jj in seq_len(nptm)) {
            sumz <- 0; sumw <- 0
            for (k in seq_len(npt)) {
              sumz <- sumz + zmat[k, jj]
              vlag[k] <- wsum[k] * zmat[k, jj]
              sumw <- sumw + vlag[k]
            }
            for (j in seq_len(n)) {
              sum <- (fracsq * sumz - 0.5 * sumw) * xopt[j]
              for (k in seq_len(npt)) sum <- sum + vlag[k] * xpt[k, j]
              wa[j] <- sum
              bmat[seq_len(npt), j] <- bmat[seq_len(npt), j] + sum * zmat[, jj]
            }
            for (i in seq_len(n)) {
              ip <- i + npt
              temp <- wa[i]
              for (j in seq_len(i)) bmat[ip, j] <- bmat[ip, j] + temp * wa[j]
            }
          }
          ih <- 0L
          for (j in seq_len(n)) {
            wa[j] <- -0.5 * sumpq * xopt[j]
            for (k in seq_len(npt)) {
              wa[j] <- wa[j] + pq[k] * xpt[k, j]
              xpt[k, j] <- xpt[k, j] - xopt[j]
            }
            for (i in seq_len(j)) {
              ih <- ih + 1L
              hq[ih] <- hq[ih] + wa[i] * xopt[j] + xopt[i] * wa[j]
              bmat[npt + i, j] <- bmat[npt + j, i]
            }
          }
          xbase <- xbase + xopt
          xnew <- xnew - xopt
          sl <- sl - xopt
          su <- su - xopt
          xopt[] <- 0
          xoptsq <- 0
        }
        lab <- if (ntrits == 0L) "L210" else "L230"
      },
      L190 = {
        nfsav <- stop$nevals
        kbase <- kopt
        rs <- .nbq_rescue(n, npt, xl, xu, stop, calfun, state(), delta)
        xbase <- rs$xbase; xpt <- rs$xpt; fval <- rs$fval; xopt <- rs$xopt; gopt <- rs$gopt; hq <- rs$hq
        pq <- rs$pq; bmat <- rs$bmat; zmat <- rs$zmat; sl <- rs$sl; su <- rs$su; kopt <- rs$kopt; vlag <- rs$vlag
        xoptsq <- 0
        if (kopt != kbase) for (i in seq_len(n)) {
          xopt[i] <- xpt[kopt, i]
          xoptsq <- xoptsq + xopt[i] * xopt[i]
        }
        if (rs$ret != "SUCCESS") {
          rc <- rs$ret
          lab <- "L720"
        } else {
          nresc <- stop$nevals
          if (nfsav < stop$nevals) {
            nfsav <- stop$nevals
            lab <- "L20"
          } else {
            lab <- if (ntrits > 0L) "L60" else "L210"
          }
        }
      },
      L210 = {
        am <- .nbq_altmov(n, npt, xpt, xopt, bmat, zmat, sl, su, kopt, knew, adelt)
        xnew <- am$xnew; xalt <- am$xalt; alpha <- am$alpha; cauchy <- am$cauchy
        d <- xnew - xopt
        lab <- "L230"
      },
      L230 = {
        w1 <- numeric(npt)
        for (k in seq_len(npt)) {
          suma <- 0; sumb <- 0; sum <- 0
          for (j in seq_len(n)) {
            suma <- suma + xpt[k, j] * d[j]
            sumb <- sumb + xpt[k, j] * xopt[j]
            sum <- sum + bmat[k, j] * d[j]
          }
          w1[k] <- suma * (0.5 * suma + sumb)
          vlag[k] <- sum
          wsum[k] <- suma
        }
        beta <- 0
        for (jj in seq_len(nptm)) {
          sum <- 0
          for (k in seq_len(npt)) sum <- sum + zmat[k, jj] * w1[k]
          beta <- beta - sum * sum
          vlag[seq_len(npt)] <- vlag[seq_len(npt)] + sum * zmat[, jj]
        }
        dsq <- 0; bsum <- 0; dx <- 0
        for (j in seq_len(n)) {
          dsq <- dsq + d[j] * d[j]
          sum <- 0
          for (k in seq_len(npt)) sum <- sum + w1[k] * bmat[k, j]
          bsum <- bsum + sum * d[j]
          jp <- npt + j
          for (i in seq_len(n)) sum <- sum + bmat[jp, i] * d[i]
          vlag[jp] <- sum
          bsum <- bsum + sum * d[j]
          dx <- dx + d[j] * xopt[j]
        }
        beta <- dx * dx + dsq * (xoptsq + dx + dx + 0.5 * dsq) + beta - bsum
        vlag[kopt] <- vlag[kopt] + 1
        lab <- "L360"
        if (ntrits == 0L) {
          denom <- vlag[knew] * vlag[knew] + alpha * beta
          if (denom < cauchy && cauchy > 0) {
            xnew <- xalt
            d <- xnew - xopt
            cauchy <- 0
            lab <- "L230"
          } else if (denom <= 0.5 * (vlag[knew] * vlag[knew])) {
            if (stop$nevals > nresc) {
              lab <- "L190"
            } else {
              rc <- "ROUNDOFF_LIMITED"
              lab <- "L720"
            }
          }
        } else {
          delsq <- delta * delta
          scaden <- 0; biglsq <- 0
          knew <- 0L
          for (k in seq_len(npt)) {
            if (k == kopt) next
            hdiag <- 0
            for (jj in seq_len(nptm)) hdiag <- hdiag + zmat[k, jj] * zmat[k, jj]
            den <- beta * hdiag + vlag[k] * vlag[k]
            distsq <- 0
            for (j in seq_len(n)) {
              t <- xpt[k, j] - xopt[j]
              distsq <- distsq + t * t
            }
            t <- distsq / delsq
            temp <- .nbq_max2(1, t * t)
            if (temp * den > scaden) {
              scaden <- temp * den
              knew <- k
              denom <- den
            }
            biglsq <- .nbq_max2(biglsq, temp * (vlag[k] * vlag[k]))
          }
          if (scaden <= 0.5 * biglsq) {
            if (stop$nevals > nresc) {
              lab <- "L190"
            } else {
              rc <- "ROUNDOFF_LIMITED"
              lab <- "L720"
            }
          }
        }
      },
      L360 = {
        x <- .nbq_point(xl, xu, xbase, xnew, sl, su)
        if (.nlopt_stop_forced(stop)) rc <- "FORCED_STOP"
        else if (.nlopt_stop_evals(stop)) rc <- "MAXEVAL_REACHED"
        else if (.nlopt_stop_time(stop)) rc <- "MAXTIME_REACHED"
        if (rc != "SUCCESS") { lab <- "L720"; next }
        stop$nevals <- stop$nevals + 1L
        f <- calfun(x)
        if (ntrits == -1L) {
          fsave <- f
          rc <- "XTOL_REACHED"
          if (fsave < fval[kopt]) return(list(ret = rc, x = x, minf = f))
          lab <- "L720"
          next
        }
        if (f < stop$minf_max) return(list(ret = "STOPVAL_REACHED", x = x, minf = f))
        fopt <- fval[kopt]
        vquad <- 0
        ih <- 0L
        for (j in seq_len(n)) {
          vquad <- vquad + d[j] * gopt[j]
          for (i in seq_len(j)) {
            ih <- ih + 1L
            temp <- d[i] * d[j]
            if (i == j) temp <- 0.5 * temp
            vquad <- vquad + hq[ih] * temp
          }
        }
        for (k in seq_len(npt)) vquad <- vquad + 0.5 * pq[k] * (wsum[k] * wsum[k])
        diff <- f - fopt - vquad
        diffc <- diffb
        diffb <- diffa
        diffa <- abs(diff)
        if (dnorm > rho) nfsav <- stop$nevals
        if (ntrits > 0L) {
          if (vquad >= 0) {
            rc <- "ROUNDOFF_LIMITED"; lab <- "L720"; next
          }
          ratio <- (f - fopt) / vquad
          if (ratio <= 0.1) {
            delta <- .nbq_min2(0.5 * delta, dnorm)
          } else if (ratio <= 0.7) {
            delta <- .nbq_max2(0.5 * delta, dnorm)
          } else {
            delta <- .nbq_max2(0.5 * delta, dnorm + dnorm)
          }
          if (delta <= rho * 1.5) delta <- rho
          if (f < fopt) {
            ksav <- knew
            densav <- denom
            delsq <- delta * delta
            scaden <- 0; biglsq <- 0
            knew <- 0L
            for (k in seq_len(npt)) {
              hdiag <- 0
              for (jj in seq_len(nptm)) hdiag <- hdiag + zmat[k, jj] * zmat[k, jj]
              den <- beta * hdiag + vlag[k] * vlag[k]
              distsq <- 0
              for (j in seq_len(n)) {
                t <- xpt[k, j] - xnew[j]
                distsq <- distsq + t * t
              }
              t <- distsq / delsq
              temp <- .nbq_max2(1, t * t)
              if (temp * den > scaden) {
                scaden <- temp * den
                knew <- k
                denom <- den
              }
              biglsq <- .nbq_max2(biglsq, temp * (vlag[k] * vlag[k]))
            }
            if (scaden <= 0.5 * biglsq) {
              knew <- ksav
              denom <- densav
            }
          }
        }
        up <- .nbq_update(n, npt, bmat, zmat, vlag, beta, denom, knew)
        bmat <- up$bmat; zmat <- up$zmat; vlag <- up$vlag
        ih <- 0L
        pqold <- pq[knew]
        pq[knew] <- 0
        for (i in seq_len(n)) {
          temp <- pqold * xpt[knew, i]
          for (j in seq_len(i)) {
            ih <- ih + 1L
            hq[ih] <- hq[ih] + temp * xpt[knew, j]
          }
        }
        for (jj in seq_len(nptm)) {
          temp <- diff * zmat[knew, jj]
          pq <- pq + temp * zmat[, jj]
        }
        fval[knew] <- f
        xpt[knew, ] <- xnew
        wv <- bmat[knew, ]
        for (k in seq_len(npt)) {
          suma <- 0
          for (jj in seq_len(nptm)) suma <- suma + zmat[knew, jj] * zmat[k, jj]
          if (.nlopt_isinf(suma)) {
            rc <- "ROUNDOFF_LIMITED"
            break
          }
          sumb <- 0
          for (j in seq_len(n)) sumb <- sumb + xpt[k, j] * xopt[j]
          temp <- suma * sumb
          wv <- wv + temp * xpt[k, ]
        }
        if (rc != "SUCCESS") { lab <- "L720"; next }
        gopt <- gopt + diff * wv
        if (f < fopt) {
          kopt <- knew
          xoptsq <- 0
          ih <- 0L
          for (j in seq_len(n)) {
            xopt[j] <- xnew[j]
            xoptsq <- xoptsq + xopt[j] * xopt[j]
            for (i in seq_len(j)) {
              ih <- ih + 1L
              if (i < j) gopt[j] <- gopt[j] + hq[ih] * d[i]
              gopt[i] <- gopt[i] + hq[ih] * d[j]
            }
          }
          for (k in seq_len(npt)) {
            temp <- 0
            for (j in seq_len(n)) temp <- temp + xpt[k, j] * d[j]
            temp <- pq[k] * temp
            gopt <- gopt + temp * xpt[k, ]
          }
          if (.nlopt_stop_ftol(stop, f, fopt)) {
            rc <- "FTOL_REACHED"; lab <- "L720"; next
          }
        }
        if (ntrits > 0L) {
          vlag[seq_len(npt)] <- fval - fval[kopt]
          w1 <- numeric(npt)
          for (j in seq_len(nptm)) {
            sum <- 0
            for (k in seq_len(npt)) sum <- sum + zmat[k, j] * vlag[k]
            w1 <- w1 + sum * zmat[, j]
          }
          w2 <- numeric(npt)
          for (k in seq_len(npt)) {
            sum <- 0
            for (j in seq_len(n)) sum <- sum + xpt[k, j] * xopt[j]
            w2[k] <- w1[k]
            w1[k] <- sum * w1[k]
          }
          gqsq <- 0; gisq <- 0
          for (i in seq_len(n)) {
            sum <- 0
            for (k in seq_len(npt)) sum <- sum + bmat[k, i] * vlag[k] + xpt[k, i] * w1[k]
            if (xopt[i] == sl[i]) {
              t <- .nbq_min2(0, gopt[i]); gqsq <- gqsq + t * t
              t <- .nbq_min2(0, sum); gisq <- gisq + t * t
            } else if (xopt[i] == su[i]) {
              t <- .nbq_max2(0, gopt[i]); gqsq <- gqsq + t * t
              t <- .nbq_max2(0, sum); gisq <- gisq + t * t
            } else {
              gqsq <- gqsq + gopt[i] * gopt[i]
              gisq <- gisq + sum * sum
            }
            vlag[npt + i] <- sum
          }
          itest <- itest + 1L
          if (gqsq < 10 * gisq) itest <- 0L
          if (itest >= 3L) {
            gopt <- vlag[npt + seq_len(n)]
            pq <- w2
            hq[] <- 0
            itest <- 0L
          }
        }
        if (ntrits == 0L || f <= fopt + 0.1 * vquad) {
          lab <- "L60"
        } else {
          t2 <- 2 * delta; t10 <- 10 * rho
          distsq <- .nbq_max2(t2 * t2, t10 * t10)
          lab <- "L650"
        }
      },
      L650 = {
        knew <- 0L
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
          dist <- sqrt(distsq)
          if (ntrits == -1L) {
            delta <- .nbq_min2(0.1 * delta, 0.5 * dist)
            if (delta <= rho * 1.5) delta <- rho
          }
          ntrits <- 0L
          adelt <- .nbq_max2(.nbq_min2(0.1 * dist, delta), rho)
          dsq <- adelt * adelt
          lab <- "L90"
        } else if (ntrits == -1L) {
          lab <- "L680"
        } else if (ratio > 0 || .nbq_max2(delta, dnorm) > rho) {
          lab <- "L60"
        } else {
          lab <- "L680"
        }
      },
      L680 = {
        if (rho > rhoend) {
          delta <- 0.5 * rho
          ratio <- rho / rhoend
          if (ratio <= 16) rho <- rhoend
          else if (ratio <= 250) rho <- sqrt(ratio) * rhoend
          else rho <- 0.1 * rho
          delta <- .nbq_max2(delta, rho)
          ntrits <- 0L
          nfsav <- stop$nevals
          lab <- "L60"
        } else {
          lab <- if (ntrits == -1L) "L360" else "L720"
        }
      },
      L720 = {
        x <- .nbq_point(xl, xu, xbase, xopt, sl, su)
        return(list(ret = rc, x = x, minf = fval[kopt]))
      })
  }
}

# bobyqa: LN_BOBYQA from x with initial steps dx (scaled so that the first
# is rhobeg) and 2n + 1 interpolation conditions; returns list(ret, x, minf)
.nlopt_bobyqa <- function(f, lb, ub, x, dx, stop) {
  n <- length(x)
  npt <- 2L * n + 1L
  s <- .nlopt_compute_rescaling(n, dx)
  bad <- function(msg) {
    stop$msg <- msg
    list(ret = "INVALID_ARGS", x = x * s, minf = Inf)
  }
  for (j in seq_len(n)) if (s[j] == 0 || !is.finite(s[j]))
    return(bad(sprintf("invalid scaling %g of dimension %d: possible over/underflow?", s[j], j - 1L)))
  x <- x / s
  bnd <- .nlopt_reorder_bounds(lb / s, ub / s)
  xl <- bnd$lb; xu <- bnd$ub
  rhobeg <- abs(dx[1] / s[1])
  rhoend <- stop$xtol_rel * rhobeg
  if (!is.null(stop$xtol_abs)) {
    xa <- rep_len(stop$xtol_abs, n)
    for (j in seq_len(n)) if (rhoend < xa[j] / abs(s[j])) rhoend <- xa[j] / abs(s[j])
  }
  np <- n + 1L
  if (npt < n + 2L || npt > ((n + 2L) * np) %/% 2L) return(bad("invalid number of sampling points"))
  sl <- numeric(n); su <- numeric(n)
  for (j in seq_len(n)) {
    temp <- xu[j] - xl[j]
    if (temp < rhobeg + rhobeg)
      return(bad(sprintf("insufficient space between the bounds: %g - %g < %g", xu[j], xl[j], rhobeg + rhobeg)))
    sl[j] <- xl[j] - x[j]
    su[j] <- xu[j] - x[j]
    if (sl[j] >= -rhobeg) {
      if (sl[j] >= 0) {
        x[j] <- xl[j]; sl[j] <- 0; su[j] <- temp
      } else {
        x[j] <- xl[j] + rhobeg
        sl[j] <- -rhobeg
        su[j] <- .nbq_max2(xu[j] - x[j], rhobeg)
      }
    } else if (su[j] <= rhobeg) {
      if (su[j] <= 0) {
        x[j] <- xu[j]; sl[j] <- -temp; su[j] <- 0
      } else {
        x[j] <- xu[j] - rhobeg
        sl[j] <- .nbq_min2(xl[j] - x[j], -rhobeg)
        su[j] <- rhobeg
      }
    }
  }
  r <- .nbq_bobyqb(n, npt, x, xl, xu, rhobeg, rhoend, stop, function(xx) f(xx * s), sl, su)
  r$x <- r$x * s
  r
}
