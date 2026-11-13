# LINCOA of PRIMA (libprima 0.7.1, fortran/lincoa), the solver that
# OptimizationPRIMA.jl calls through PRIMA.jl: Powell's method for
# minimization without derivatives subject to bounds and linear
# constraints. Translated from the Fortran reference with the arithmetic of
# the libprima build, as R/optimization-prima-bobyqa.R; the common routines
# of powalg are those of R/optimization-prima-powell.R and the filter of the
# selectx module is that of libprima (selectx.f90).

.lc_tinycv <- 2^-200
.lc_funcmax <- 2^100

# maximum of linalg: NaN when an element is NaN
.lc_maximum <- function(v) {
  s <- sum(abs(v))
  if (is.nan(s)) s else max(v)
}

# isbetter, savefilt and selectx of the selectx module of libprima
.lc_isbetter <- function(f1, c1, f2, c2, ctol) {
  ib <- (is.nan(f1) || is.nan(c1)) && !(is.nan(f2) || is.nan(c2))
  ib <- ib || isTRUE(f1 < f2 && c1 <= c2)
  ib <- ib || isTRUE(f1 <= f2 && c1 < c2)
  cref <- 10 * max(.pr_eps, min(ctol, 1e-2 * .lc_funcmax))
  ib || isTRUE(f1 < .pr_realmax && c1 <= ctol && (c2 > max(ctol, cref) || is.nan(c2)))
}
.lc_savefilt <- function(filt, cstrv, ctol, cweight, f, x) {
  nf <- filt$nfilt
  idx <- seq_len(nf)
  for (i in idx) if (.lc_isbetter(filt$ffilt[i], filt$cfilt[i], f, cstrv, ctol)) return(filt)
  keep <- !vapply(idx, function(i) .lc_isbetter(f, cstrv, filt$ffilt[i], filt$cfilt[i], ctol), TRUE)
  maxfilt <- length(filt$ffilt)
  if (sum(keep) == maxfilt) {
    cfs <- .bq_pmax(filt$cfilt - ctol, 0)
    phi <- if (cweight <= 0) filt$ffilt else if (cweight == Inf) cfs else
      .bq_pmax(filt$ffilt, -.pr_realmax) + cweight * cfs
    phimax <- .uo_maxval(phi)
    cref <- .uo_maxval(cfs[which(phi >= phimax)])
    fref <- .uo_maxval(filt$ffilt[which(cfs >= cref)])
    cand <- which(filt$ffilt >= fref)
    kworst <- if (length(cand)) cand[.pw_maxloc(filt$cfilt[cand])] else 1L
    keep[kworst] <- FALSE
  }
  kept <- which(keep)
  nk <- length(kept)
  if (nk) {
    filt$xfilt[, seq_len(nk)] <- filt$xfilt[, kept, drop = FALSE]
    filt$ffilt[seq_len(nk)] <- filt$ffilt[kept]
    filt$cfilt[seq_len(nk)] <- filt$cfilt[kept]
  }
  filt$xfilt[, nk + 1L] <- x
  filt$ffilt[nk + 1L] <- f
  filt$cfilt[nk + 1L] <- cstrv
  filt$nfilt <- nk + 1L
  filt
}
.lc_selectx <- function(fh, ch, cweight, ctol) {
  anyw <- function(v) any(v, na.rm = TRUE)
  fm <- .lc_funcmax
  if (anyw(fh < fm & ch < fm)) {
    fref <- fm; cref <- fm
  } else if (anyw(fh < .pr_realmax & ch < fm)) {
    fref <- .pr_realmax; cref <- fm
  } else if (anyw(fh < fm & ch < .pr_realmax)) {
    fref <- fm; cref <- .pr_realmax
  } else {
    fref <- .pr_realmax; cref <- .pr_realmax
  }
  if (!anyw(fh < fref & ch < cref)) return(length(fh))
  chs <- .bq_pmax(ch - ctol, 0)
  cmin <- .uo_minval(chs[which(fh < fref)])
  cref <- max(.pr_eps, 2 * cmin)
  phi <- if (cweight <= 0) fh else if (cweight == Inf) chs else .bq_pmax(fh, -.pr_realmax) + cweight * chs
  phimin <- .uo_minval(phi[which(fh < fref & chs <= cref)])
  cref <- .uo_minval(chs[which(fh < fref & phi <= phimin)])
  fref <- .uo_minval(fh[which(chs <= cref)])
  cand <- which(fh <= fref)
  cand[.pw_minloc(ch[cand])]
}

# hypotenuse of linalg
.lc_hypot <- function(x1, x2) {
  if (!is.finite(x1)) return(abs(x1))
  if (!is.finite(x2)) return(abs(x2))
  y <- sort(abs(c(x1, x2)))
  r <- if (y[1L] > sqrt(.Machine$double.xmin) && y[2L] < sqrt(.pr_realmax / 2.1)) sqrt(.pw_inprod(y, y))
    else if (y[2L] > 0) { t <- y[1L] / y[2L]; y[2L] * sqrt(.pw_fma(t, t, 1)) }
    else 0
  min(y[1L] + y[2L], max(y[2L], r))
}

# solve of linalg for a triangular A (forward substitution when A is lower
# triangular, back substitution otherwise)
.lc_solve <- function(A, b) {
  n <- nrow(A)
  x <- numeric(n)
  if (!n) return(x)
  if (all(A[upper.tri(A)] == 0)) {
    for (i in seq_len(n)) x[i] <- (b[i] - .pw_inprod(A[i, seq_len(i - 1L)], x[seq_len(i - 1L)])) / A[i, i]
  } else {
    for (i in n:1) {
      j <- seq_len(n)[seq_len(n) > i]
      x[i] <- (b[i] - .pw_inprod(A[i, j], x[j])) / A[i, i]
    }
  }
  x
}

# lsqr_Rfull of linalg: the least squares solution of A x = b with A = Q R
.lc_lsqr <- function(b, Q, R) {
  n <- ncol(R)
  x <- .pw_matprod12(b, Q)
  for (i in rev(seq_len(n))) {
    for (j in seq_len(n)[seq_len(n) > i]) x[i] <- .pw_fma(-R[i, j], x[j], x[i])
    x[i] <- x[i] / R[i, i]
  }
  x
}

# qradd_Rfull of powalg: adds the column cvec to the QR factorization whose
# first nact columns are active; returns list(Q, R, nact)
.lc_qradd <- function(cvec, Q, R, nact) {
  m <- nrow(Q)
  cq <- .pw_matprod12(cvec, Q)
  for (k in rev(seq_len(m - 1L))) {
    if (k < nact + 1L) break
    if (abs(cq[k + 1L]) > 0) {
      G <- .pw_planerot(cq[c(k, k + 1L)])
      Q[, c(k, k + 1L)] <- .pw_matprod22(Q[, c(k, k + 1L), drop = FALSE], t(G))
      cq[k] <- sqrt(.pw_fma(cq[k], cq[k], cq[k + 1L] * cq[k + 1L]))
    }
  }
  if (nact) R[seq_len(nact), nact + 1L] <- .pw_matprod12(cvec, Q[, seq_len(nact), drop = FALSE])
  if (cq[nact + 1L] < 0) Q[, nact + 1L] <- -Q[, nact + 1L]
  R[nact + 1L, nact + 1L] <- abs(cq[nact + 1L])
  list(Q = Q, R = R, nact = nact + 1L)
}

# qrexc_Rfull of powalg: moves active column i to position nact; R holds the
# first nact columns of rfac
.lc_qrexc <- function(Q, R, i, nact) {
  if (i <= 0L || i >= nact) return(list(Q = Q, R = R))
  for (k in i:(nact - 1L)) {
    G <- .pw_planerot(R[c(k + 1L, k), k + 1L])
    hypt <- .lc_hypot(R[k + 1L, k + 1L], R[k, k + 1L])
    Q[, c(k, k + 1L)] <- .pw_matprod22(Q[, c(k + 1L, k), drop = FALSE], t(G))
    R[c(k, k + 1L), k:nact] <- .pw_matprod22(G, R[c(k + 1L, k), k:nact, drop = FALSE])
    R[seq_len(k + 1L), c(k, k + 1L)] <- R[seq_len(k + 1L), c(k + 1L, k)]
    R[c(k, k + 1L), k] <- c(hypt, 0)
  }
  list(Q = Q, R = R)
}

# addact and delact of getact; st holds iact, nact, qfac, rfac, resact,
# resnew and vlam
.lc_addact <- function(l, cvec, st) {
  qr <- .lc_qradd(cvec, st$qfac, st$rfac, st$nact)
  st$qfac <- qr$Q
  st$rfac <- qr$R
  st$nact <- qr$nact
  st$iact[st$nact] <- l
  st$resact[st$nact] <- st$resnew[l]
  st$resnew[l] <- 0
  st$vlam[st$nact] <- 0
  st
}
.lc_delact <- function(icon, st) {
  nact <- st$nact
  qr <- .lc_qrexc(st$qfac, st$rfac, icon, nact)
  st$qfac <- qr$Q
  st$rfac <- qr$R
  rot <- c(seq_len(nact)[seq_len(nact) > icon], icon)
  st$iact[icon:nact] <- st$iact[rot]
  st$resact[icon:nact] <- st$resact[rot]
  st$resnew[st$iact[nact]] <- .pw_fmaxnm(st$resact[nact], .lc_tinycv)
  st$vlam[icon:nact] <- st$vlam[rot]
  st$nact <- nact - 1L
  st
}

# the projection of v onto the complement of the active gradients
.lc_project <- function(Q, nact, v) {
  n <- nrow(Q)
  Q2 <- Q[, seq_len(n)[seq_len(n) > nact], drop = FALSE]
  .pw_matprod21(Q2, .pw_matprod12(v, Q2))
}

# getact: the active set and the steepest descent direction psd projected
# onto it
.lc_getact <- function(amat, delta, g, st) {
  m <- ncol(amat)
  n <- length(g)
  if (m <= 0L) {
    st$nact <- 0L
    st$qfac <- diag(n)
    st$psd <- -g
    return(st)
  }
  gg <- .pw_inprod(g, g)
  tdel <- 0.2 * delta
  if (st$nact == 0L) st$qfac <- diag(n)
  st$vlam <- numeric(n)
  for (icon in rev(seq_len(st$nact))) if (st$resact[icon] > tdel) st <- .lc_delact(icon, st)
  while (st$nact > 0L) {
    a <- seq_len(st$nact)
    st$vlam[a] <- .lc_lsqr(g, st$qfac[, a, drop = FALSE], st$rfac[a, a, drop = FALSE])
    if (!any(st$vlam[a] >= 0, na.rm = TRUE)) break
    st <- .lc_delact(max(which(st$vlam[a] >= 0)), st)
  }
  psdsav <- numeric(n)
  ddsav <- 2 * gg
  psd <- numeric(n)
  for (iter in seq_len(max(2L * (m + n), m))) {
    if (st$nact >= n) {
      psd <- numeric(n)
      break
    }
    psd <- -.lc_project(st$qfac, st$nact, g)
    dd <- .pw_inprod(psd, psd)
    dnorm <- sqrt(dd)
    if (is.nan(dnorm) || dnorm <= .pr_eps) break
    if (dd >= ddsav) {
      psd <- numeric(n)
      break
    }
    if (isTRUE(.pw_inprod(psd, g) > 0) || !is.finite(sum(abs(psd)))) {
      psd <- psdsav
      break
    }
    psdsav <- psd
    ddsav <- dd
    apsd <- .pw_matprod12(psd, amat)
    mask <- st$resnew > 0 & st$resnew <= tdel & apsd > (dnorm / delta) * st$resnew
    mask[is.na(mask)] <- FALSE
    if (any(mask)) {
      l <- which(mask)[.pw_maxloc(apsd[mask])]
      violmx <- apsd[l]
    } else {
      l <- 0L
      violmx <- -.pr_realmax
    }
    ainf <- if (st$nact) .uo_maxval(abs(apsd[st$iact[seq_len(st$nact)]])) else 0
    if (!any(mask) || !isTRUE(violmx > .pw_fmaxnm(.pr_eps * dnorm, 10 * ainf))) break
    st <- .lc_addact(l, amat[, l], st)
    while (isTRUE(violmx > 0) && st$nact > 0L) {
      nact <- st$nact
      a <- seq_len(nact)
      v <- numeric(nact)
      v[nact] <- 1 / st$rfac[nact, nact]
      vmu <- .lc_solve(st$rfac[a, a, drop = FALSE], v)
      frac <- rep(.pr_realmax, nact)
      w <- which(vmu < 0 & st$vlam[a] < 0)
      frac[w] <- st$vlam[w] / vmu[w]
      vmult <- .uo_minval(c(violmx, frac))
      icon <- max(0L, which(frac <= vmult))
      violmx <- .pw_fmaxnm(violmx - vmult, 0)
      st$vlam[a] <- .pw_fma(-vmult, vmu, st$vlam[a])
      if (icon > 0L && icon <= nact) st$vlam[icon] <- 0
      for (ic in rev(a)) if (isTRUE(st$vlam[ic] >= 0)) st <- .lc_delact(ic, st)
    }
    if (st$nact == 0L) break
  }
  if (st$nact == 0L) {
    st$qfac <- diag(n)
    psd <- -g
  }
  st$psd <- psd
  st
}

# trstep: the trust region step subject to the linear constraints by
# truncated conjugate gradients on the active set; returns list(s, ngetact,
# iact, nact, qfac, rfac)
.lc_trstep <- function(amat, delta, gopt_in, hq_in, pq_in, rescon, xpt, iact, nact, qfac, rfac) {
  m <- ncol(amat)
  n <- length(gopt_in)
  ctest <- 0.01
  gmax <- max(abs(gopt_in))
  if (isTRUE(gmax > 999999995904)) {
    modscal <- .pw_fmaxnm(2 * .Machine$double.xmin, 1 / gmax)
    gopt <- gopt_in * modscal
    pq <- pq_in * modscal
    hq <- hq_in * modscal
  } else {
    gopt <- gopt_in
    pq <- pq_in
    hq <- hq_in
  }
  out <- function(s, ngetact) list(s = s, ngetact = ngetact, iact = st$iact, nact = st$nact, qfac = st$qfac,
                                   rfac = st$rfac)
  st <- list(iact = iact, nact = nact, qfac = qfac, rfac = rfac)
  if (!is.finite(sum(abs(gopt)))) return(out(numeric(n), 0L))
  resnew <- rescon
  nn <- which(rescon >= 0)
  resnew[nn] <- .bq_pmax(.lc_tinycv, rescon[nn])
  resnew[which(rescon >= delta)] <- -1
  a <- seq_len(nact)
  resnew[iact[a]] <- 0
  st$resact <- numeric(n)
  st$resact[a] <- rescon[iact[a]]
  st$resnew <- resnew
  st$vlam <- numeric(n)
  g <- gopt
  delsq <- delta * delta
  s <- numeric(n)
  ss <- 0
  reduct <- 0
  ngetact <- 0L
  newact <- TRUE
  itercg <- -1L
  gamma <- 0
  d <- numeric(n)
  dproj <- numeric(n)
  sqrtd3 <- function(ds, dd, resid) .uo_maxval(c(sqrt(.pw_fma(ds, ds, dd * resid)), abs(ds), sqrt(dd * resid)))
  for (iter in seq_len(min(10000L, 10L * (m + n)))) {
    if (newact) {
      ngetact <- ngetact + 1L
      st <- .lc_getact(amat, delta, g, st)
      psd <- st$psd
      dd <- .pw_inprod(psd, psd)
      if (is.nan(dd) || dd <= .pr_eps * delsq) break
      psd <- ((0.2 * delta) / sqrt(dd)) * psd
      gamma <- 0
      a <- seq_len(st$nact)
      if (any(st$resact[a] > 1e-4 * delta)) {
        dproj <- .pw_matprod21(st$qfac[, a, drop = FALSE], .lc_solve(t(st$rfac[a, a, drop = FALSE]), st$resact[a]))
        sp <- s + psd
        ds <- .pw_inprod(dproj, sp)
        dd <- .pw_inprod(dproj, dproj)
        resid <- delsq - .pw_inprod(sp, sp)
        if (resid > 0 && dd > .pr_eps * delsq && !is.nan(ds)) {
          sqrtd <- sqrtd3(ds, dd, resid)
          gamma <- if (ds <= 0) (sqrtd - ds) / dd else resid / (sqrtd + ds)
          if (!is.finite(gamma) || gamma < 0) gamma <- 0
          ad <- rep(-1, m)
          pos <- which(st$resnew > 0)
          ad[pos] <- .pw_matprod12(dproj, amat[, pos, drop = FALSE])
          frac <- rep(1, m)
          pa <- which(ad > 0)
          frac[pa] <- (st$resnew[pa] - .pw_matprod12(psd, amat[, pa, drop = FALSE])) / ad[pa]
          gamma <- .uo_minval(c(gamma, 1, frac))
        }
      }
      if (gamma > 0) {
        d <- .pw_fma(gamma, dproj, psd)
        itercg <- -1L
      } else {
        d <- psd
        itercg <- 0L
      }
    }
    itercg <- itercg + 1L
    resid <- delsq - ss
    dg <- .pw_inprod(d, g)
    ds <- .pw_inprod(d, s)
    dd <- .pw_inprod(d, d)
    if (resid <= 0 || is.nan(ds) || dd <= .pr_eps * delsq) break
    sqrtd <- sqrtd3(ds, dd, resid)
    alpha <- if (ds <= 0) (sqrtd - ds) / dd else resid / (sqrtd + ds)
    if (!is.finite(alpha) || alpha <= 0) break
    if (is.nan(alpha * dg) || -alpha * dg <= ctest * reduct) break
    hd <- .pw_hess_mul(d, xpt, pq, hq)
    dhd <- .pw_inprod(d, hd)
    alpht <- alpha
    if (isTRUE(.pw_fma(alpha, dhd, dg) > 0)) alpha <- -dg / dhd
    alphm <- alpha
    ad <- rep(-1, m)
    pos <- which(st$resnew > 0)
    ad[pos] <- .pw_matprod12(d, amat[, pos, drop = FALSE])
    frac <- rep(alpha, m)
    pa <- which(ad > 0)
    frac[pa] <- st$resnew[pa] / ad[pa]
    frac[is.nan(frac)] <- alpha
    jsav <- 0L
    if (any(frac < alpha)) {
      jsav <- .pw_minloc(frac)
      alpha <- frac[jsav]
    }
    alpha <- if (itercg == 0L) 1 else if (itercg == 1L && gamma <= 0) .pw_fmaxnm(alpha, 1) else .pw_fmaxnm(alpha, 0)
    alpha <- .pw_fminnm(alpha, alphm)
    sold <- s
    s <- .pw_fma(alpha, d, s)
    ss <- .pw_inprod(s, s)
    if (!is.finite(ss)) {
      s <- sold
      break
    }
    g <- .pw_fma(alpha, hd, g)
    if (!is.finite(sum(abs(g)))) break
    restmp <- .pw_fma(-alpha, ad, st$resnew)
    pos <- which(st$resnew > 0)
    st$resnew[pos] <- .bq_pmax(.lc_tinycv, restmp[pos])
    if (itercg == 0L) {
      a <- seq_len(st$nact)
      st$resact[a] <- .pw_fma(-alpha, gamma, 1) * st$resact[a]
    }
    reduct <- .pw_fma(-alpha, .pw_fma(0.5 * alpha, dhd, dg), reduct)
    if (is.nan(reduct) || reduct <= 0) {
      s <- sold
      break
    }
    if (isTRUE(alpha >= alpht) || isTRUE(-alphm * .pw_fma(0.5 * alphm, dhd, dg) <= ctest * reduct)) break
    newact <- jsav > 0L
    if (newact) next
    if (itercg >= n - st$nact) break
    pg <- if (st$nact <= 0L) g else .lc_project(st$qfac, st$nact, g)
    beta <- if (itercg == 0L) 0 else .pw_inprod(pg, hd) / dhd
    d <- .pw_fma(beta, d, -pg)
  }
  out(s, ngetact)
}

# trrad of LINCOA
.lc_trrad <- function(delta_in, dnorm, eta1, eta2, gamma1, gamma2, ratio) {
  if (isTRUE(ratio <= eta1)) gamma1 * dnorm
  else if (isTRUE(ratio <= eta2)) .pw_fmaxnm(gamma1 * delta_in, dnorm)
  else .pw_fmaxnm(gamma1 * delta_in, gamma2 * dnorm)
}

# setdrop_tr of LINCOA
.lc_setdrop_tr <- function(idz, kopt, ximproved, bmat, d, delta, rho, xpt, zmat) {
  distsq <- .pw_distsq(xpt, if (ximproved) xpt[, kopt] + d else xpt[, kopt])
  r <- .pw_fmaxnm(0.1 * delta, rho)
  w <- .pw_vfmaxnm(distsq / (r * r), 1)
  score <- ((w * w) * w) * abs(.pw_calden(kopt, bmat, d, xpt, zmat, idz))
  if (!ximproved) score[kopt] <- -1
  if (any(score > 0, na.rm = TRUE)) return(.pw_maxloc(score))
  if (ximproved) return(.pw_maxloc(distsq))
  0L
}

# geostep of LINCOA; returns list(s, feasible)
.lc_geostep <- function(iact, idz, knew, kopt, nact, amat, bmat, delbar, qfac, rescon, xpt, zmat) {
  n <- nrow(xpt)
  xopt <- xpt[, kopt]
  pqlag <- .pw_omega_col(idz, zmat, knew)
  glag <- bmat[, knew] + .pw_hess_mul(xopt, xpt, pqlag)
  dderiv <- .pw_matprod12(glag, xpt) - .pw_inprod(glag, xopt)
  distsq <- .pw_distsq(xpt, xopt)
  distsq[kopt] <- 1
  stplen <- -delbar / sqrt(distsq)
  vlagabs <- abs((stplen * (1 - stplen)) * dderiv)
  dk <- dderiv[knew]
  if (isTRUE(dk * (dk - 1) < 0)) stplen[knew] <- -stplen[knew]
  sk <- stplen[knew]
  vlagabs[knew] <- .pw_fma(sk * sk, abs(dk - 1), abs(sk * dk))
  vlagabs[kopt] <- -1
  k <- knew
  if (any(vlagabs > vlagabs[knew], na.rm = TRUE)) k <- .pw_maxloc(vlagabs)
  s <- stplen[k] * (xpt[, k] - xopt)
  denabs <- abs(.pw_calden(kopt, bmat, s, xpt, zmat, idz)[knew])
  curv <- function(v) .pw_inprod(v, .pw_hess_mul(v, xpt, pqlag))
  gnorm <- .pw_norm(glag)
  if (is.finite(gnorm) && gnorm > .pr_eps) {
    gstp <- (delbar / gnorm) * glag
    if (isTRUE(curv(gstp) < 0)) gstp <- -gstp
    dg <- abs(.pw_calden(kopt, bmat, gstp, xpt, zmat, idz)[knew])
    if (isTRUE(dg > denabs) || is.nan(denabs)) {
      denabs <- dg
      s <- gstp
    }
  }
  rstat <- rep(1L, length(rescon))
  rstat[abs(rescon) >= delbar] <- -1L
  rstat[iact[seq_len(nact)]] <- 0L
  viol <- function(v, cols) .lc_maximum(c(0, .pw_matprod12(v, amat[, cols, drop = FALSE]) - rescon[cols]))
  cstrv <- viol(s, which(rstat >= 0L))
  feasible <- isTRUE(cstrv <= 0)
  pglag <- .lc_project(qfac, nact, glag)
  gnorm <- .pw_norm(pglag)
  if (nact > 0L && is.finite(gnorm) && gnorm > .pr_eps) {
    pgstp <- (delbar / gnorm) * pglag
    if (isTRUE(curv(pgstp) < 0)) pgstp <- -pgstp
    cstrv <- viol(pgstp, which(rstat == 1L))
    cvtol <- .pw_fmaxnm(.pr_eps * .pw_norm(pgstp),
                        10 * .uo_maxval(abs(.pw_matprod12(pgstp, amat[, iact[seq_len(nact)], drop = FALSE]))))
    take <- FALSE
    if (isTRUE(cstrv <= cvtol)) take <- isTRUE(abs(.pw_calden(kopt, bmat, pgstp, xpt, zmat, idz)[knew]) > 0.1 * denabs)
    if (take || is.nan(denabs)) {
      s <- pgstp
      feasible <- isTRUE(cstrv <= cvtol)
    }
  }
  if (is.nan(sum(abs(s)))) {
    s <- xpt[, knew] - xopt
    scaling <- delbar / .pw_norm(s)
    s <- .pw_fmaxnm(0.6 * scaling, .pw_fminnm(0.5, scaling)) * s
    cstrv <- viol(s, which(rstat >= 0L))
    feasible <- isTRUE(cstrv <= 0)
  }
  list(s = s, feasible = feasible)
}

# updateres: the residuals of the constraints at the new XOPT, negative for
# the constraints that are far from it
.lc_updateres <- function(ximproved, amat, b, delta, dnorm, xopt, rescon) {
  if (!ximproved) return(rescon)
  mask <- abs(rescon) < dnorm + delta
  mask[is.na(mask)] <- FALSE
  out <- .bq_pmin(-abs(rescon) + dnorm, -delta)
  w <- which(mask)
  out[w] <- .bq_pmax(b[w] - .pw_matprod12(xopt, amat[, w, drop = FALSE]), 0)
  far <- which(out >= delta)
  out[far] <- -out[far]
  out
}

# inith of LINCOA: the initial [BMAT, ZMAT] with IDZ = 1
.lc_inith <- function(ij, xpt) {
  n <- nrow(xpt)
  npt <- ncol(xpt)
  rhobeg <- max(abs(xpt[, 2L]))
  rhosq <- rhobeg * rhobeg
  recip <- 1 / rhobeg
  reciq <- 0.5 / rhobeg
  bmat <- matrix(0, n, npt + n)
  if (npt <= 2L * n + 1L) {
    k1 <- seq_len(npt - n - 1L)
    bmat[cbind(k1, k1 + 1L)] <- reciq
    bmat[cbind(k1, k1 + n + 1L)] <- -reciq
    k2 <- seq_len(n)[seq_len(n) >= npt - n]
    bmat[k2, 1L] <- -recip
    bmat[cbind(k2, k2 + 1L)] <- recip
    bmat[cbind(k2, k2 + npt)] <- -(0.5 * rhosq)
  } else {
    k1 <- seq_len(n)
    bmat[cbind(k1, k1 + 1L)] <- reciq
    bmat[cbind(k1, k1 + n + 1L)] <- -reciq
  }
  recip <- 1 / rhosq
  reciq <- sqrt(0.5) / rhosq
  zmat <- matrix(0, npt, npt - n - 1L)
  if (npt <= 2L * n + 1L) {
    k1 <- seq_len(npt - n - 1L)
    zmat[1L, ] <- -reciq - reciq
    zmat[cbind(k1 + 1L, k1)] <- reciq
    zmat[cbind(k1 + n + 1L, k1)] <- reciq
  } else {
    k1 <- seq_len(n)
    zmat[1L, k1] <- -reciq - reciq
    zmat[cbind(k1 + 1L, k1)] <- reciq
    zmat[cbind(k1 + n + 1L, k1)] <- reciq
    k3 <- seq_len(npt - 2L * n - 1L)
    zmat[1L, n + k3] <- recip
    zmat[cbind(2L * n + 1L + k3, n + k3)] <- recip
    for (k in k3) zmat[ij[, k] + 1L, k + n] <- -recip
  }
  list(bmat = bmat, zmat = zmat)
}

# the constraint values [xl - x, x - xu, -(Aeq x - beq), Aeq x - beq,
# Aineq x - bineq] of lincob at x and their violation
.lc_constr <- function(x, lc) {
  leq <- .pw_matprod21(lc$Aeq, x) - lc$beq
  constr <- c(lc$xl[lc$ixl] - x[lc$ixl], x[lc$ixu] - lc$xu[lc$ixu], -leq, leq,
              .pw_matprod21(lc$Aineq, x) - lc$bineq)
  list(constr = constr, cstrv = .lc_maximum(c(0, constr)))
}

# get_lincon: the constraints as amat' x <= bvec with normalized gradients,
# the right hand sides relaxed so that x0 is feasible
.lc_get_lincon <- function(lc, rhoend, x0) {
  n <- length(x0)
  ixl <- which(lc$xl > -.pr_realmax)
  ixu <- which(lc$xu < .pr_realmax)
  rnorm <- function(A) sqrt(vapply(seq_len(nrow(A)), function(i) .pw_inprod(A[i, ], A[i, ]), 0))
  Aeq_norm <- rnorm(lc$Aeq)
  Aineq_norm <- rnorm(lc$Aineq)
  ieq <- which(Aeq_norm > 0)
  iineq <- which(Aineq_norm > 0)
  if (length(ieq) < nrow(lc$Aeq) || length(iineq) < nrow(lc$Aineq))
    warning("LINCOA: Some linear constraints have zero gradients; they are ignored", call. = FALSE)
  I <- diag(n)
  Ae <- t(lc$Aeq[ieq, , drop = FALSE])
  amat <- cbind(-I[, ixl, drop = FALSE], I[, ixu, drop = FALSE], -Ae, Ae, t(lc$Aineq[iineq, , drop = FALSE]))
  bvec <- c(-lc$xl[ixl], lc$xu[ixu], -lc$beq[ieq], lc$beq[ieq], lc$bineq[iineq])
  Aeqx0 <- .pw_matprod21(lc$Aeq, x0)
  Aineqx0 <- .pw_matprod21(lc$Aineq, x0)
  bvec <- .bq_pmax(bvec, c(-x0[ixl], x0[ixu], -Aeqx0[ieq], Aeqx0[ieq], Aineqx0[iineq]))
  Anorm <- c(Aeq_norm[ieq], Aeq_norm[ieq], Aineq_norm[iineq])
  lin <- length(ixl) + length(ixu) + seq_along(Anorm)
  amat[, lin] <- amat[, lin, drop = FALSE] / rep(Anorm, each = n)
  bvec[lin] <- bvec[lin] / Anorm
  smallx <- 1e-6 * rhoend
  if (any(x0 + smallx < lc$xl) || any(x0 - smallx > lc$xu) || any(abs(Aeqx0 - lc$beq) > smallx * Aeq_norm) ||
      any(Aineqx0 - lc$bineq > smallx * Aineq_norm))
    warning("LINCOA: The starting point is infeasible. LINCOA modified the right-hand sides of the constraints ",
            "to make it feasible", call. = FALSE)
  list(amat = amat, bvec = bvec, ixl = ixl, ixu = ixu)
}

# lincob: the main loop of LINCOA (with initxf)
.lc_lincob <- function(fun, iprint, maxfilt, maxfun, npt, lc, amat, bvec, ctol, cweight, ftarget, rhobeg, rhoend,
                       x0) {
  eta1 <- 0.1
  eta2 <- 0.7
  gamma1 <- 0.5
  gamma2 <- 2
  n <- length(x0)
  m <- ncol(amat)
  filt <- list(nfilt = 0L, xfilt = matrix(0, n, maxfilt), ffilt = numeric(maxfilt), cfilt = numeric(maxfilt))
  # initxf
  xbase <- x0
  xpt <- matrix(0, n, npt)
  for (k in seq_len(n)) xpt[k, k + 1L] <- rhobeg
  for (k in seq_len(min(n, npt - n - 1L))) xpt[k, k + n + 1L] <- -rhobeg
  ij <- .nw_setij(n, npt)
  if (npt > 2L * n + 1L)
    for (k in (2L * n + 2L):npt) xpt[, k] <- xpt[, ij[1L, k - 2L * n - 1L] + 1L] + xpt[, ij[2L, k - 2L * n - 1L] + 1L]
  b <- bvec - .pw_matprod12(xbase, amat)
  cval <- vapply(seq_len(npt), function(k) {
    v <- .lc_maximum(c(0, .pw_matprod12(xpt[, k], amat) - b))
    if (is.nan(v)) .pr_realmax else v
  }, 0)
  feasible <- cval <= 0
  evaluated <- logical(npt)
  fval <- rep(.pr_realmax, npt)
  info <- 0L
  for (k in seq_len(npt)) {
    x <- xbase + xpt[, k]
    f <- .pw_evaluatef(fun, x)
    cc <- .lc_constr(x, lc)
    .pw_fmsg(iprint, "Initialization", k, rhobeg, f, x, cc$cstrv, cc$constr)
    evaluated[k] <- TRUE
    cval[k] <- cc$cstrv
    fval[k] <- f
    info <- .pr_checkbreak_con(maxfun, k, cc$cstrv, ctol, f, ftarget, x)
    if (info != 0L) break
  }
  nf <- sum(evaluated)
  feasible[.bq_minloc_mask(cval, evaluated)] <- TRUE
  kopt <- .bq_minloc_mask(fval, evaluated & feasible)
  for (k in which(evaluated)) filt <- .lc_savefilt(filt, cval[k], ctol, cweight, fval[k], xbase + xpt[, k])
  finish <- function(info) {
    nfl <- seq_len(filt$nfilt)
    kf <- .lc_selectx(filt$ffilt[nfl], filt$cfilt[nfl], cweight, ctol)
    x <- filt$xfilt[, kf]
    f <- filt$ffilt[kf]
    cc <- .lc_constr(x, lc)
    .pw_retmsg("LINCOA", iprint, info, nf, f, x, cc$cstrv, cc$constr)
    list(x = x, f = f, cstrv = cc$cstrv, nf = nf, info = info)
  }
  if (info != 0L) return(finish(info))
  h <- .lc_inith(ij, xpt)
  idz <- 1L
  bmat <- h$bmat
  zmat <- h$zmat
  hq <- matrix(0, n, n)
  pq <- .pw_omega_mul(idz, zmat, fval)
  gopt <- .pw_matprod21(bmat[, seq_len(npt), drop = FALSE], fval) + .pw_hess_mul(xpt[, kopt], xpt, pq)
  pqalt <- pq
  galt <- gopt
  rescon <- .bq_pmax(b - .pw_matprod12(xpt[, kopt], amat), 0)
  rescon[which(rescon >= rhobeg)] <- -rescon[which(rescon >= rhobeg)]
  rho <- rhobeg
  delta <- rho
  ratio <- -1
  dnorm_rec <- rep(.pr_realmax, 4L)
  shortd <- FALSE
  trfail <- FALSE
  qalt_better <- logical(3L)
  knew_tr <- 0L
  qfac <- diag(n)
  rfac <- matrix(0, n, n)
  nact <- 0L
  iact <- seq_len(m)
  gamma3 <- .pw_fmaxnm(1, .pw_fminnm(0.75 * gamma2, 1.5))
  maxtr <- if (maxfun > .Machine$integer.max %/% 2L) maxfun else 2L * maxfun
  info <- 20L
  d <- numeric(n)
  cc <- .lc_constr(xbase + xpt[, kopt], lc)
  # evaluates f at XBASE + XOPT + d and saves it in the filter
  evalx <- function(state, rad) {
    x <- xbase + (xpt[, kopt] + d)
    f <<- .pw_evaluatef(fun, x)
    nf <<- nf + 1L
    cc <<- .lc_constr(x, lc)
    .pw_fmsg(iprint, state, nf, rad, f, x, cc$cstrv, cc$constr)
    filt <<- .lc_savefilt(filt, cc$cstrv, ctol, cweight, f, x)
    .pr_checkbreak_con(maxfun, nf, cc$cstrv, ctol, f, ftarget, x)
  }
  # replaces XPT(:, knew) by XOPT + d and updates the models
  swap_point <- function(knew, ximproved, moderr) {
    xdrop <- xpt[, knew]
    xosav <- xpt[, kopt]
    uh <- .pw_updateh(knew, kopt, d, xpt, idz, bmat, zmat)
    idz <<- uh$idz
    bmat <<- uh$bmat
    zmat <<- uh$zmat
    xpt[, knew] <<- xosav + d
    fval[knew] <<- f
    if (ximproved) kopt <<- knew
    q <- .nw_updateq(idz, knew, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
    fo <- fval - fval[kopt]
    pqalt <<- .pw_omega_mul(idz, zmat, fo)
    galt <<- .pw_matprod21(bmat[, seq_len(npt), drop = FALSE], fo) + .pw_hess_mul(xpt[, kopt], xpt, pqalt)
    if (all(qalt_better)) {
      q <- list(gopt = galt, hq = matrix(0, n, n), pq = pqalt)
      qalt_better <<- logical(3L)
    }
    gopt <<- q$gopt
    hq <<- q$hq
    pq <<- q$pq
    rescon <<- .lc_updateres(ximproved, amat, b, delta, .pw_norm(d), xpt[, kopt], rescon)
  }
  f <- fval[kopt]
  for (tr in seq_len(maxtr)) {
    ts <- .lc_trstep(amat, delta, gopt, hq, pq, rescon, xpt, iact, nact, qfac, rfac)
    d <- ts$s
    iact <- ts$iact
    nact <- ts$nact
    qfac <- ts$qfac
    rfac <- ts$rfac
    dnorm <- .pw_fminnm(delta, .pw_norm(d))
    shortd <- isTRUE((dnorm < 0.5 * delta && ts$ngetact < 2L) || dnorm < 0.1999 * delta)
    dnorm_rec <- c(dnorm_rec[-1L], dnorm)
    if (delta > rho || !shortd) dnorm_rec[] <- .pr_realmax
    qred <- -.pw_quadinc(d, xpt, gopt, pq, hq)
    trfail <- !isTRUE(qred > 9.99999974737875e-06 * (rho * rho))
    if (shortd || trfail) {
      delta <- 0.5 * delta
      if (delta <= gamma3 * rho) delta <- rho
    } else {
      sub <- evalx("Trust region", delta)
      if (sub != 0L) {
        info <- sub
        break
      }
      fopt <- fval[kopt]
      moderr <- (f - fopt) + qred
      moderr_alt <- (f - fopt) - .pw_quadinc(d, xpt, galt, pqalt, NULL)
      qalt_better <- c(qalt_better[-1L], isTRUE(abs(moderr_alt) < 0.1 * abs(moderr)))
      ratio <- .pr_redrat(fopt - f, qred, eta1)
      delta <- .lc_trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
      if (delta <= gamma3 * rho) delta <- rho
      ximproved <- isTRUE(f < fopt)
      knew_tr <- .lc_setdrop_tr(idz, kopt, ximproved, bmat, d, delta, rho, xpt, zmat)
      if (knew_tr > 0L) swap_point(knew_tr, ximproved, moderr)
    }
    accurate_mod <- all(dnorm_rec <= 0.5 * rho) || all(dnorm_rec[3:4] <= 0.20000000298023224 * rho)
    distsq <- .pw_distsq(xpt, xpt[, kopt])
    close_itpset <- all(distsq <= 4 * (delta * delta))
    adequate_geo <- isTRUE((shortd && accurate_mod) || close_itpset)
    small_trrad <- isTRUE(.pw_fmaxnm(delta, dnorm) <= rho)
    improve_geo <- (shortd || trfail || isTRUE(ratio <= eta1) || knew_tr == 0L) && !adequate_geo
    reduce_rho <- (shortd || trfail || isTRUE(ratio <= 0) || knew_tr == 0L) && adequate_geo && small_trrad
    if (improve_geo) {
      knew_geo <- .pw_maxloc(distsq)
      delbar <- .pw_fmaxnm(0.1 * delta, rho)
      gs <- .lc_geostep(iact, idz, knew_geo, kopt, nact, amat, bmat, delbar, qfac, rescon, xpt, zmat)
      d <- gs$s
      sub <- evalx("Geometry", delbar)
      if (sub != 0L) {
        info <- sub
        break
      }
      fopt <- fval[kopt]
      moderr <- (f - fopt) - .pw_quadinc(d, xpt, gopt, pq, hq)
      moderr_alt <- (f - fopt) - .pw_quadinc(d, xpt, galt, pqalt, NULL)
      qalt_better <- c(qalt_better[-1L], isTRUE(abs(moderr_alt) < 0.1 * abs(moderr)))
      swap_point(knew_geo, isTRUE(f < fopt) && gs$feasible, moderr)
    }
    if (reduce_rho) {
      if (rho <= rhoend) {
        info <- 0L
        break
      }
      rhonew <- .pr_redrho(rho, rhoend)
      delta <- .pw_fmaxnm(0.5 * rho, rhonew)
      rho <- rhonew
      .pw_rhomsg(iprint, nf, delta, fval[kopt], rho, xbase + xpt[, kopt], cc$cstrv, cc$constr)
      dnorm_rec[] <- .pr_realmax
    }
    xopt <- xpt[, kopt]
    if (isTRUE(.pw_inprod(xopt, xopt) >= 1000 * (delta * delta))) {
      b <- b - .pw_matprod12(xopt, amat)
      sb <- .pw_shiftbase(kopt, xbase, xpt, zmat, bmat, pq, hq, idz)
      xbase <- sb$xbase
      xpt <- sb$xpt
      bmat <- sb$bmat
      hq <- sb$hq
      fo <- fval - fval[kopt]
      pqalt <- .pw_omega_mul(idz, zmat, fo)
      galt <- .pw_matprod21(bmat[, seq_len(npt), drop = FALSE], fo) + .pw_hess_mul(xpt[, kopt], xpt, pqalt)
    }
  }
  if (info == 0L && shortd && nf < maxfun) evalx("Trust region", rho)
  finish(info)
}

# lincoa: the constraints Aineq x <= bineq, Aeq x = beq and xl <= x <= xu as
# libprima prepares them, then lincob
.lc_lincoa <- function(fun, x0, Aineq = NULL, bineq = NULL, Aeq = NULL, beq = NULL, xl = NULL, xu = NULL,
                       rhobeg = 1, rhoend = 1e-6 * rhobeg, ftarget = -Inf, maxfun = 500L * length(x0),
                       npt = 2L * length(x0) + 1L, iprint = 0L) {
  n <- length(x0)
  x0 <- .pw_moderatex(as.numeric(x0))
  mat <- function(A) if (is.null(A)) matrix(0, 0L, n) else matrix(as.numeric(A), ncol = n)
  lc <- list(Aineq = mat(Aineq), bineq = as.numeric(bineq), Aeq = mat(Aeq), beq = as.numeric(beq),
             xl = if (is.null(xl)) rep(-.pr_realmax, n) else as.numeric(xl),
             xu = if (is.null(xu)) rep(.pr_realmax, n) else as.numeric(xu))
  op <- .pw_preproc("LINCOA", n, maxfun, ftarget, rhobeg, rhoend, npt)
  gl <- .lc_get_lincon(lc, op$rhoend, x0)
  lc$ixl <- gl$ixl
  lc$ixu <- gl$ixu
  .lc_lincob(fun, as.integer(iprint), min(op$maxfun, 2000L), op$maxfun, op$npt, lc, gl$amat, gl$bvec, sqrt(.pr_eps), 1e8,
             op$ftarget, op$rhobeg, op$rhoend, x0)
}
