# TNC of SciPy 1.18.1 (BSD 3-clause license): the driver _minimize_tnc of
# scipy/optimize/_tnc.py and scipy/optimize/tnc/tnc.c, the C version by
# Jean-Sebastien Roy (MIT license) of the TNBC Fortran code of Stephen G.
# Nash. SciPy compiles tnc.c with clang, which fuses a * b + c within one
# expression into a multiply add (-ffp-contract=on); .fma() stands at each
# such site.

.sptn_eps <- 2.220446049250313e-16

# RCSTRINGS of _tnc.py, for the return codes -1 to 7
.sptn_messages <- c("Infeasible (lower bound > upper bound)", "Local minimum reached (|pg| ~= 0)",
                    "Converged (|f_n-f_(n-1)| ~= 0)", "Converged (|x_n-x_(n-1)| ~= 0)",
                    "Max. number of function evaluations reached", "Linear search failed",
                    "All lower bounds are equal to the upper bounds", "Unable to progress",
                    "User requested end of minimization")

# tnc_rc_string of tnc.c, for the return codes -3 to 7
.sptn_c_messages <- c("Memory allocation failed", "Invalid parameters (n<0)", "Infeasible (low bound > up bound)",
                      "Local minima reach (|pg| ~= 0)", "Converged (|f_n-f_(n-1)| ~= 0)",
                      "Converged (|x_n-x_(n-1)| ~= 0)", "Maximum number of function evaluations reached",
                      "Linear search failed", "All lower bounds are equal to the upper bounds", "Unable to progress",
                      "User requested end of minimization")

# ddot1 (dtemp += dy[i] * dx[i], fused)
.sptn_dot <- function(x, y) {
  s <- 0
  for (i in seq_along(x)) s <- .fma(y[i], x[i], s)
  s
}

# ddot1 as the SciPy wheel compiles it inside msolve: clang vectorizes the
# loop there, so blocks of four products are formed first and added in order
# without fusion, and only the last n %% 4 terms are fused (n < 4: all fused).
# The ddot1 calls inlined into tnc_minimize stay scalar and fused.
.sptn_dot4 <- function(x, y) {
  s <- 0
  n <- length(x)
  m <- if (n >= 4L) n - n %% 4L else 0L
  for (i in seq_len(m)) s <- s + x[i] * y[i]
  for (i in seq_len(n - m) + m) s <- .fma(y[i], x[i], s)
  s
}

# dnrm21: the scaled sum of squares
.sptn_nrm2 <- function(x) {
  ssq <- 1
  scale <- 0
  for (xi in x) {
    if (!isTRUE(xi == 0)) {
      a <- abs(xi)
      if (isTRUE(scale < a)) {
        r <- scale / a
        ssq <- .fma(ssq * r, r, 1)
        scale <- a
      } else {
        r <- a / scale
        ssq <- .fma(r, r, ssq)
      }
    }
  }
  scale * sqrt(ssq)
}

.sptn_coerce <- function(x, low, up) {
  i <- which(x < low)
  x[i] <- low[i]
  i <- which(x > up)
  x[i] <- up[i]
  x
}

.sptn_unscale <- function(x, xscale, xoffset) .fma(x, xscale, xoffset)

.sptn_scale <- function(x, xscale, xoffset) {
  i <- which(xscale > 0)
  x[i] <- (x[i] - xoffset[i]) / xscale[i]
  x
}

.sptn_project <- function(x, pivot) {
  x[pivot != 0L] <- 0
  x
}

# setConstraints; the test at the upper bound is the one of tnc.c, which
# marks x active only when it lies beyond up by the tolerance
.sptn_set_constraints <- function(x, xscale, xoffset, low, up) {
  n <- length(x)
  pivot <- integer(n)
  for (i in seq_len(n)) {
    if (xscale[i] == 0) {
      pivot[i] <- 2L
    } else if (!isTRUE(low[i] == -Inf) &&
               isTRUE(.fma(x[i], xscale[i], xoffset[i]) - low[i] <= .sptn_eps * 10 * (abs(low[i]) + 1))) {
      pivot[i] <- -1L
    } else if (!isTRUE(up[i] == Inf) &&
               isTRUE(.fma(x[i], xscale[i], xoffset[i]) - up[i] >= .sptn_eps * 10 * (abs(up[i]) + 1))) {
      pivot[i] <- 1L
    }
  }
  pivot
}

.sptn_step_max <- function(step, x, dir, pivot, low, up, xscale, xoffset) {
  for (i in seq_along(x)) {
    if (pivot[i] == 0L && !isTRUE(dir[i] == 0)) {
      if (isTRUE(dir[i] < 0)) {
        t <- (low[i] - xoffset[i]) / xscale[i] - x[i]
        if (isTRUE(t > step * dir[i])) step <- t / dir[i]
      } else {
        t <- (up[i] - xoffset[i]) / xscale[i] - x[i]
        if (isTRUE(t < step * dir[i])) step <- t / dir[i]
      }
    }
  }
  step
}

.sptn_add_constraint <- function(x, p, pivot, low, up, xscale, xoffset) {
  newcon <- FALSE
  for (i in seq_along(x)) {
    if (pivot[i] == 0L && !isTRUE(p[i] == 0)) {
      if (isTRUE(p[i] < 0) && !isTRUE(low[i] == -Inf)) {
        tol <- .sptn_eps * 10 * (abs(low[i]) + 1)
        if (isTRUE(.fma(x[i], xscale[i], xoffset[i]) - low[i] <= tol)) {
          pivot[i] <- -1L
          x[i] <- (low[i] - xoffset[i]) / xscale[i]
          newcon <- TRUE
        }
      } else if (!isTRUE(up[i] == Inf)) {
        tol <- .sptn_eps * 10 * (abs(up[i]) + 1)
        if (isTRUE(up[i] - .fma(x[i], xscale[i], xoffset[i]) <= tol)) {
          pivot[i] <- 1L
          x[i] <- (up[i] - xoffset[i]) / xscale[i]
          newcon <- TRUE
        }
      }
    }
  }
  list(x = x, pivot = pivot, newcon = newcon)
}

.sptn_remove_constraint <- function(gtpnew, gnorm, pgtolfs, f, flast, g, pivot) {
  if (isTRUE((flast - f) <= (gtpnew * -0.5)) && isTRUE(gnorm > pgtolfs)) return(list(pivot = pivot, removed = FALSE))
  imax <- 0L
  cmax <- 0
  for (i in seq_along(g)) {
    if (pivot[i] == 2L) next
    t <- -pivot[i] * g[i]
    if (isTRUE(t < cmax)) {
      cmax <- t
      imax <- i
    }
  }
  if (imax > 0L) pivot[imax] <- 0L
  list(pivot = pivot, removed = imax > 0L)
}

.sptn_initial_step <- function(fnew, fmin, gtp, smax) {
  d <- abs(fnew - fmin)
  alpha <- 1
  if (isTRUE(d * 2 <= -gtp) && d >= .sptn_eps) alpha <- (d * -2) / gtp
  if (isTRUE(alpha >= smax)) alpha <- smax
  alpha
}

# self scaled BFGS
.sptn_ssbfgs <- function(gamma, sj, hjv, hjyj, yjsj, yjhyj, vsj, vhyj) {
  if (isTRUE(yjsj == 0)) {
    delta <- 0
    beta <- 0
  } else {
    delta <- (gamma * yjhyj / yjsj + 1) * vsj / yjsj - gamma * vhyj / yjsj
    beta <- -gamma * vsj / yjsj
  }
  .fma(beta, hjyj, .fma(gamma, hjv, delta * sj))
}

# msolve: the two step self scaled BFGS preconditioner applied to g
.sptn_msolve <- function(g, P) {
  if (P$upd1) return(g / P$diagb)
  gsk <- .sptn_dot4(g, P$sk)
  rdiagb <- 1 / P$diagb
  hg <- g * rdiagb
  hyk <- P$yk * rdiagb
  if (P$lreset) {
    ykhyk <- .sptn_dot4(P$yk, hyk)
    ghyk <- .sptn_dot4(g, hyk)
    return(.sptn_ssbfgs(1, P$sk, hg, hyk, P$yksk, ykhyk, gsk, ghyk))
  }
  hyr <- P$yr * rdiagb
  gsr <- .sptn_dot4(g, P$sr)
  ghyr <- .sptn_dot4(g, hyr)
  yrhyr <- .sptn_dot4(P$yr, hyr)
  hg <- .sptn_ssbfgs(1, P$sr, hg, hyr, P$yrsr, yrhyr, gsr, ghyr)
  yksr <- .sptn_dot4(P$yk, P$sr)
  ykhyr <- .sptn_dot4(P$yk, hyr)
  hyk <- .sptn_ssbfgs(1, P$sr, hyk, hyr, P$yrsr, yrhyr, yksr, ykhyr)
  ykhyk <- .sptn_dot4(hyk, P$yk)
  ghyk <- .sptn_dot4(hyk, g)
  .sptn_ssbfgs(1, P$sk, hg, hyk, P$yksk, ykhyk, gsk, ghyk)
}

.sptn_init_preconditioner <- function(P) {
  if (P$upd1) return(P$diagb)
  td <- P$diagb
  yksk <- P$yksk
  if (P$lreset) {
    bsk <- td * P$sk
    sds <- .sptn_dot(P$sk, bsk)
    if (yksk == 0) yksk <- 1
    if (sds == 0) sds <- 1
    return(td - td * td * P$sk * P$sk / sds + P$yk * P$yk / yksk)
  }
  yrsr <- P$yrsr
  bsk <- td * P$sr
  sds <- .sptn_dot(P$sr, bsk)
  srds <- .sptn_dot(P$sk, bsk)
  yrsk <- .sptn_dot(P$yr, P$sk)
  if (yrsr == 0) yrsr <- 1
  if (sds == 0) sds <- 1
  bsk <- .fma(td, P$sk, -(bsk * srds / sds)) + P$yr * yrsk / yrsr
  emat <- td - td * td * P$sr * P$sr / sds + P$yr * P$yr / yrsr
  sds <- .sptn_dot(P$sk, bsk)
  if (yksk == 0) yksk <- 1
  if (sds == 0) sds <- 1
  emat - (bsk * bsk / sds + P$yk * P$yk / yksk)
}

.sptn_diagonal_scaling <- function(e, v, gv, r) {
  vr <- 1 / .sptn_dot(v, r)
  vgv <- 1 / .sptn_dot(v, gv)
  e <- e + .fma(-r * r, vr, gv * gv * vgv)
  e[which(e <= 1e-6)] <- 1
  e
}

# the Hessian times v by a forward difference of the gradient
.sptn_hessian_times_vector <- function(v, x, g, Q) {
  delta <- Q$accuracy * (Q$xnorm + 1)
  xv <- .sptn_coerce(.sptn_unscale(.fma(delta, v, x), Q$xscale, Q$xoffset), Q$low, Q$up)
  gv <- Q$fg(xv)$g * (Q$xscale * Q$fscale)
  gv <- (gv - g) * (1 / delta)
  gv[Q$xscale == 0] <- 0
  gv
}

# tnc_direction: preconditioned conjugate gradients on the Newton equations;
# returns the direction and the new diagonal preconditioner
.sptn_direction <- function(x, g, maxCGit, maxnfeval, P, Q, S, pivot, gnorm) {
  free <- pivot == 0L
  if (maxCGit == 0L) {
    zsol <- -g
    zsol[!free] <- 0
    return(list(zsol = zsol, diagb = P$diagb))
  }
  rhsnrm <- gnorm
  tol <- 1e-12
  qold <- 0
  rzold <- 0
  n <- length(x)
  emat <- .sptn_init_preconditioner(P)
  r <- -g
  v <- numeric(n)
  zsol <- numeric(n)
  for (k in seq_len(maxCGit) - 1L) {
    r[!free] <- 0
    zk <- .sptn_msolve(r, P)
    zk[!free] <- 0
    rz <- .sptn_dot(r, zk)
    if (isTRUE(rz / rhsnrm < tol) || S$nfeval >= maxnfeval - 1L) {
      if (k == 0L) {
        zsol <- -g
        zsol[!free] <- 0
      }
      break
    }
    beta <- if (k == 0L) 0 else rz / rzold
    v <- .fma(beta, v, zk)
    v[!free] <- 0
    gv <- .sptn_hessian_times_vector(v, x, g, Q)
    S$nfeval <- S$nfeval + 1L
    gv[!free] <- 0
    vgv <- .sptn_dot(v, gv)
    if (isTRUE(vgv / rhsnrm < tol)) {
      if (k == 0L) {
        zsol <- -.sptn_msolve(g, P)
        zsol[!free] <- 0
      }
      break
    }
    emat <- .sptn_diagonal_scaling(emat, v, gv, r)
    alpha <- rz / vgv
    zsol <- .fma(alpha, v, zsol)
    r <- .fma(-alpha, gv, r)
    gtp <- .sptn_dot(zsol, g)
    pr <- .sptn_dot(r, zsol)
    qnew <- (gtp + pr) * 0.5
    qtest <- (k + 1) * (1 - qold / qnew)
    if (isTRUE(qtest <= 0.5)) break
    if (isTRUE(gtp > 0)) {
      zsol <- .fma(-alpha, v, zsol)
      break
    }
    qold <- qnew
    rzold <- rz
  }
  list(zsol = zsol, diagb = emat)
}

# getptcInit and getptcIter: the safeguarded cubic interpolation of the line
# search; L holds the variables that tnc.c passes by pointer, L$xmin being
# the step alpha. Return codes: 0 a point found, 1 evaluate f at xmin + u,
# 2 bad input, 3 no lower point
.sptn_getptc_init <- function(L, tnytol, eta, rmu, xbnd) {
  if (!isTRUE(L$u > 0) || xbnd <= tnytol || isTRUE(L$gu > 0)) return(2L)
  if (xbnd < L$abstol) L$abstol <- xbnd
  L$tol <- L$abstol
  L$a <- 0
  L$xw <- 0
  L$xmin <- 0
  L$oldf <- L$fu
  L$fmin <- L$fu
  L$fw <- L$fu
  L$gw <- L$gu
  L$gmin <- L$gu
  L$step <- L$u
  L$factor <- 5
  L$braktd <- FALSE
  L$scxbnd <- xbnd
  L$b <- .fma(L$reltol, abs(L$scxbnd), L$scxbnd) + L$abstol
  L$e <- L$b + L$b
  L$b1 <- L$b
  L$gtest1 <- -rmu * L$gu
  L$gtest2 <- -eta * L$gu
  .sptn_getptc_bound(L, xbnd)
}

# the end of both routines: the step capped by the scaled bound and kept
# at least tol away from the origin
.sptn_getptc_bound <- function(L, xbnd) {
  if (isTRUE(L$step >= L$scxbnd)) {
    L$step <- L$scxbnd
    L$scxbnd <- L$scxbnd - .fma(L$reltol, abs(xbnd), L$abstol) / (1 + L$reltol)
  }
  L$u <- L$step
  if (isTRUE(abs(L$step) < L$tol) && isTRUE(L$step < 0)) L$u <- -L$tol
  if (isTRUE(abs(L$step) < L$tol) && isTRUE(L$step >= 0)) L$u <- L$tol
  1L
}

.sptn_getptc_iter <- function(L, big, rtsmll, tnytol, fpresn, xbnd) {
  moved <- FALSE
  if (isTRUE(L$fu <= L$fmin)) {
    chordu <- .fma(-(L$xmin + L$u), L$gtest1, L$oldf)
    if (isTRUE(L$fu > chordu)) {
      # no sufficient decrease: the upper bound moves to this point
      chordm <- .fma(-L$xmin, L$gtest1, L$oldf)
      L$gu <- -L$gmin
      denom <- chordm - L$fmin
      if (isTRUE(abs(denom) < 1e-15)) {
        denom <- 1e-15
        if (isTRUE(chordm - L$fmin < 0)) denom <- -denom
      }
      if (L$xmin != 0) L$gu <- L$gmin * (chordu - L$fu) / denom
      L$fu <- .fma(0.5 * L$u, L$gmin + L$gu, L$fmin)
      if (isTRUE(L$fu < L$fmin)) L$fu <- L$fmin
    } else {
      L$fw <- L$fmin
      L$fmin <- L$fu
      L$gw <- L$gmin
      L$gmin <- L$gu
      L$xmin <- L$xmin + L$u
      L$a <- L$a - L$u
      L$b <- L$b - L$u
      L$xw <- -L$u
      L$scxbnd <- L$scxbnd - L$u
      if (isTRUE(L$gu <= 0)) {
        L$a <- 0
      } else {
        L$b <- 0
        L$braktd <- TRUE
      }
      L$tol <- .fma(abs(L$xmin), L$reltol, L$abstol)
      moved <- TRUE
    }
  }
  if (!moved) {
    if (isTRUE(L$u < 0)) {
      L$a <- L$u
    } else {
      L$b <- L$u
      L$braktd <- TRUE
    }
    L$xw <- L$u
    L$fw <- L$fu
    L$gw <- L$gu
  }
  twotol <- L$tol + L$tol
  xmidpt <- 0.5 * (L$a + L$b)
  convrg <- isTRUE(abs(xmidpt) <= .fma(-0.5, L$b - L$a, twotol)) ||
    (isTRUE(abs(L$gmin) <= L$gtest2) && isTRUE(L$fmin < L$oldf) &&
       (isTRUE(abs(L$xmin - xbnd) > L$tol) || !L$braktd))
  if (convrg) {
    if (L$xmin != 0) return(0L)
    # no reduction: check the change of f against tol and shrink tol
    if (isTRUE(abs(L$oldf - L$fw) <= fpresn)) return(3L)
    L$tol <- 0.1 * L$tol
    if (L$tol < tnytol) return(3L)
    L$reltol <- 0.1 * L$reltol
    L$abstol <- 0.1 * L$abstol
    twotol <- 0.1 * twotol
  }
  r <- 0
  q <- 0
  s <- 0
  if (isTRUE(abs(L$e) > L$tol)) {
    # the cubic through xmin and xw
    r <- 3 * (L$fmin - L$fw) / L$xw + L$gmin + L$gw
    absr <- abs(r)
    q <- absr
    found <- FALSE
    if (!isTRUE(L$gw == 0) && !isTRUE(L$gmin == 0)) {
      abgw <- abs(L$gw)
      abgmin <- abs(L$gmin)
      s <- sqrt(abgmin) * sqrt(abgw)
      if (isTRUE(L$gw / abgw * L$gmin > 0)) {
        if (isTRUE(r >= s) || isTRUE(r <= -s)) {
          q <- sqrt(abs(r + s)) * sqrt(abs(r - s))
        } else {
          r <- 0
          q <- 0
          found <- TRUE
        }
      } else {
        sumsq <- 1
        p <- 0
        if (isTRUE(absr >= s)) {
          if (absr > rtsmll) p <- absr * rtsmll
          if (isTRUE(s >= p)) {
            value <- s / absr
            sumsq <- .fma(value, value, 1)
          }
          scale <- absr
        } else {
          if (s > rtsmll) p <- s * rtsmll
          if (isTRUE(absr >= p)) {
            value <- absr / s
            sumsq <- .fma(value, value, 1)
          }
          scale <- s
        }
        sumsq <- sqrt(sumsq)
        q <- big
        if (isTRUE(scale < big / sumsq)) q <- scale * sumsq
      }
    }
    if (!found) {
      if (isTRUE(L$xw < 0)) q <- -q
      s <- L$xw * (L$gmin - r - q)
      q <- L$gw - L$gmin + q + q
      if (isTRUE(q > 0)) s <- -s
      if (isTRUE(q <= 0)) q <- -q
      r <- L$e
      if (!isTRUE(L$b1 == L$step) || L$braktd) L$e <- L$step
    }
  }
  # an artificial bound on the step
  a1 <- L$a
  L$b1 <- L$b
  L$step <- xmidpt
  if (!L$braktd || (L$a == 0 && isTRUE(L$xw < 0)) || (L$b == 0 && isTRUE(L$xw > 0))) {
    if (L$braktd) {
      d1 <- L$xw
      d2 <- L$a
      if (L$a == 0) d2 <- L$b
      L$u <- -d1 / d2
      L$step <- 5 * d2 * (0.1 + 1 / L$u) / 11
      if (isTRUE(L$u < 1)) L$step <- 0.5 * d2 * sqrt(L$u)
    } else {
      L$step <- -L$factor * L$xw
      if (isTRUE(L$step > L$scxbnd)) L$step <- L$scxbnd
      if (!isTRUE(L$step == L$scxbnd)) L$factor <- 5 * L$factor
    }
    if (isTRUE(L$step <= 0)) a1 <- L$step
    if (isTRUE(L$step > 0)) L$b1 <- L$step
  }
  # reject the interpolated step outside the interval or larger than half the
  # step of the last but one iteration
  if (isTRUE(abs(s) <= abs(0.5 * q * r)) || isTRUE(s <= q * a1) || isTRUE(s >= q * L$b1)) {
    L$e <- L$b - L$a
  } else {
    L$step <- s / q
    if (isTRUE(L$step - L$a < twotol) || isTRUE(L$b - L$step < twotol)) {
      L$step <- if (isTRUE(xmidpt <= 0)) -L$tol else L$tol
    }
  }
  .sptn_getptc_bound(L, xbnd)
}

# linearSearch (Gill and Murray): returns the code (0 a lower point, 1 too
# many evaluations, 2 failure), the step, and x, f and the gradient moved to
# the new point on success
.sptn_line_search <- function(Q, pivot, eta, ftol, xbnd, p, x, f, alpha, gfull, maxnfeval, S) {
  eps <- .sptn_eps
  maxlsit <- 64L
  gu <- .sptn_dot(gfull * (Q$xscale * Q$fscale), p)
  xnorm <- .sptn_nrm2(.sptn_project(x, pivot))
  rteps <- sqrt(eps)
  pe <- .sptn_nrm2(p) + eps
  L <- new.env(parent = emptyenv())
  L$reltol <- rteps * (xnorm + 1) / pe
  L$abstol <- -eps * (1 + abs(f)) / (gu - eps)
  tnytol <- eps * (xnorm + 1) / pe
  rtsmll <- eps
  big <- 1 / (eps * eps)
  itcnt <- 0L
  L$u <- alpha
  L$fu <- f
  L$fmin <- f
  L$gu <- gu
  L$xmin <- alpha
  itest <- .sptn_getptc_init(L, tnytol, eta, 1e-4, xbnd)
  newgfull <- gfull
  while (itest == 1L) {
    itcnt <- itcnt + 1L
    if (itcnt > maxlsit || S$nfeval >= maxnfeval) break
    ualpha <- L$xmin + L$u
    temp <- .sptn_coerce(.sptn_unscale(.fma(ualpha, p, x), Q$xscale, Q$xoffset), Q$low, Q$up)
    r <- Q$fg(temp)
    S$nfeval <- S$nfeval + 1L
    L$fu <- r$f * Q$fscale
    L$gu <- .sptn_dot(r$g * (Q$xscale * Q$fscale), p)
    itest <- .sptn_getptc_iter(L, big, rtsmll, tnytol, ftol, xbnd)
    if (isTRUE(L$xmin == ualpha)) newgfull <- r$g
  }
  if (itest == 0L) return(list(rc = 0L, alpha = L$xmin, x = .fma(L$xmin, p, x), f = L$fmin, gfull = newgfull))
  rc <- if (itcnt > maxlsit || itest != 1L) 2L else 1L
  list(rc = rc, alpha = L$xmin, x = x, f = f, gfull = gfull)
}

.sptn_print_iter <- function(f, g, niter, nfeval, pivot) {
  gtg <- 0
  for (i in which(pivot == 0L)) gtg <- .fma(g[i], g[i], gtg)
  message(sprintf(" %4d %4d %22.15E  %15.8E", niter, nfeval, f, gtg))
}

# tnc() and tnc_minimize() of tnc.c: fg(x) returns list(f, g); low and up
# are the bounds, scale and offset numeric(0) for the defaults
.sptn_tnc <- function(fg, x, low, up, scale, offset, messages, maxCGit, maxnfeval, eta, stepmx, accuracy, fmin,
                      ftol, xtol, pgtol, rescale, callback) {
  eps <- .sptn_eps
  n <- length(x)
  S <- new.env(parent = emptyenv())
  S$nfeval <- 0L
  niter <- 0L
  f <- Inf
  gfull <- numeric(n)
  done <- function(rc) {
    if (bitwAnd(messages, 8L) != 0L) message("tnc: ", .sptn_c_messages[rc + 4L])
    list(rc = rc, nfeval = S$nfeval, niter = niter, x = x, f = f, g = gfull)
  }
  if (n == 0L) return(done(5L))
  if (any(low > up, na.rm = TRUE)) return(done(-1L))
  x <- .sptn_coerce(x, low, up)
  if (maxnfeval < 1L) return(done(3L))
  r <- fg(x)
  S$nfeval <- 1L
  f <- r$f
  gfull <- r$g
  fixed <- (low == up) %in% TRUE
  if (length(scale)) fixed <- fixed | (scale == 0) %in% TRUE
  if (sum(fixed) == n) return(done(5L))
  xscale <- numeric(n)
  xoffset <- numeric(n)
  for (i in seq_len(n)) {
    if (length(scale)) {
      xscale[i] <- abs(scale[i])
      xoffset[i] <- x[i]
      if (xscale[i] == 0) {
        low[i] <- x[i]
        up[i] <- x[i]
      }
    } else if (!isTRUE(low[i] == -Inf) && !isTRUE(up[i] == Inf)) {
      xscale[i] <- up[i] - low[i]
      xoffset[i] <- (up[i] + low[i]) * 0.5
    } else {
      xscale[i] <- 1 + abs(x[i])
      xoffset[i] <- x[i]
    }
    if (length(offset)) xoffset[i] <- offset[i]
  }
  rteps <- sqrt(eps)
  if (stepmx < rteps * 10) stepmx <- 10
  if (eta < 0 || eta >= 1) eta <- 0.25
  if (rescale < 0) rescale <- 1.3
  if (maxCGit < 0L) maxCGit <- min(max(n %/% 2L, 1L), 50L)
  if (maxCGit > n) maxCGit <- n
  if (accuracy <= eps) accuracy <- rteps
  if (ftol < 0) ftol <- accuracy
  if (pgtol < 0) pgtol <- 1e-2 * sqrt(accuracy)
  if (xtol < 0) xtol <- rteps

  # tnc_minimize
  fscale <- 1
  difnew <- 0
  epsred <- 0.05
  icycle <- n - 1L
  newcon <- TRUE
  alpha <- 0
  P <- list(upd1 = TRUE, lreset = FALSE, yrsr = 0, yksk = 0, diagb = rep(1, n), sk = numeric(n), yk = numeric(n),
            sr = numeric(n), yr = numeric(n))
  Q <- list(fg = fg, xscale = xscale, xoffset = xoffset, low = low, up = up, accuracy = accuracy)
  x <- .sptn_scale(x, xscale, xoffset)
  f <- f * fscale
  pivot <- .sptn_set_constraints(x, xscale, xoffset, low, up)
  g <- gfull * (xscale * fscale)
  # a bound whose multiplier has the wrong sign is released
  pivot[which(-pivot * g < 0)] <- 0L
  g <- .sptn_project(g, pivot)
  gnorm <- .sptn_nrm2(g)
  flast_con <- f
  flast_reset <- f
  if (bitwAnd(messages, 1L) != 0L) {
    message("  NIT   NF   F                       GTG")
    .sptn_print_iter(f / fscale, gfull, niter, S$nfeval, pivot)
  }
  repeat {
    if (isTRUE(.sptn_nrm2(g) <= pgtol * fscale)) {
      g <- .sptn_project(gfull, pivot)
      if (bitwAnd(messages, 2L) != 0L) message(sprintf("tnc: |pg| = %g -> local minimum", .sptn_nrm2(g) / fscale))
      rc <- 0L
      break
    }
    if (S$nfeval >= maxnfeval) {
      rc <- 3L
      break
    }
    newscale <- .sptn_nrm2(g)
    if (isTRUE(newscale > eps) && isTRUE(abs(log10(newscale)) > rescale)) {
      newscale <- 1 / newscale
      f <- f * newscale
      fscale <- fscale * newscale
      gnorm <- gnorm * newscale
      flast_con <- flast_con * newscale
      flast_reset <- flast_reset * newscale
      difnew <- difnew * newscale
      g <- g * newscale
      P$diagb <- rep(1, n)
      P$upd1 <- TRUE
      icycle <- n - 1L
      newcon <- TRUE
      if (bitwAnd(messages, 2L) != 0L) message(sprintf("tnc: fscale = %g", fscale))
    }
    xnorm <- .sptn_nrm2(.sptn_project(x, pivot))
    oldnfeval <- S$nfeval
    Q$fscale <- fscale
    Q$xnorm <- xnorm
    d <- .sptn_direction(x, g, maxCGit, maxnfeval, P, Q, S, pivot, gnorm)
    pk <- d$zsol
    P$diagb <- d$diagb
    if (!newcon) {
      if (!P$lreset) {
        # the accumulated step and gradient change
        P$sr <- P$sr + P$sk
        P$yr <- P$yr + P$yk
        icycle <- icycle + 1L
      } else {
        P$sr <- P$sk
        P$yr <- P$yk
        flast_reset <- f
        icycle <- 1L
      }
    }
    oldg <- g
    oldf <- f
    oldgtp <- .sptn_dot(pk, g)
    ustpmax <- stepmx / (.sptn_nrm2(pk) + eps)
    spe <- .sptn_step_max(ustpmax, x, pk, pivot, low, up, xscale, xoffset)
    if (isTRUE(spe > 0)) {
      alpha <- .sptn_initial_step(f, fmin / fscale, oldgtp, spe)
      ls <- .sptn_line_search(Q, pivot, eta, ftol, spe, pk, x, f, alpha, gfull, maxnfeval, S)
      alpha <- ls$alpha
      x <- ls$x
      f <- ls$f
      gfull <- ls$gfull
      if (ls$rc == 2L) {
        rc <- 4L
        break
      }
      # the step reached the maximum unconstrained step: increase it
      if (isTRUE(alpha >= 0.9 * ustpmax)) {
        stepmx <- stepmx * 1e2
        if (bitwAnd(messages, 2L) != 0L) message(sprintf("tnc: stepmx = %g", stepmx))
      }
      # the step reached the maximum constrained step: a new constraint
      if (isTRUE(alpha - spe >= -eps * 10)) {
        newcon <- TRUE
      } else {
        if (ls$rc != 0L) {
          rc <- if (ls$rc == 1L) 3L else 4L
          break
        }
        newcon <- FALSE
      }
    } else {
      newcon <- TRUE
    }
    if (newcon) {
      ac <- .sptn_add_constraint(x, pk, pivot, low, up, xscale, xoffset)
      x <- ac$x
      pivot <- ac$pivot
      if (!ac$newcon && S$nfeval == oldnfeval) {
        rc <- 6L
        break
      }
      flast_con <- f
    }
    niter <- niter + 1L
    if (is.function(callback)) callback(.sptn_unscale(x, xscale, xoffset))
    difold <- difnew
    difnew <- oldf - f
    # the first iteration of a cycle sets the reduction factor of the reset test
    if (icycle == 1L) {
      if (isTRUE(difnew > difold * 2)) epsred <- epsred + epsred
      if (isTRUE(difnew < difold * 0.5)) epsred <- epsred * 0.5
    }
    g <- gfull * (xscale * fscale)
    gnorm <- .sptn_nrm2(.sptn_project(g, pivot))
    rmc <- .sptn_remove_constraint(oldgtp, gnorm, pgtol * fscale, f, flast_con, g, pivot)
    pivot <- rmc$pivot
    if (rmc$removed) {
      gnorm <- .sptn_nrm2(.sptn_project(g, pivot))
      flast_con <- f
    }
    if (!rmc$removed && !newcon) {
      if (isTRUE(abs(difnew) <= ftol * fscale)) {
        if (bitwAnd(messages, 2L) != 0L) message(sprintf("tnc: |fn-fn-1] = %g -> convergence", abs(difnew) / fscale))
        rc <- 1L
        break
      }
      if (isTRUE(alpha * .sptn_nrm2(pk) <= xtol)) {
        if (bitwAnd(messages, 2L) != 0L) message(sprintf("tnc: |xn-xn-1] = %g -> convergence", alpha * .sptn_nrm2(pk)))
        rc <- 2L
        break
      }
    }
    g <- .sptn_project(g, pivot)
    if (bitwAnd(messages, 1L) != 0L) .sptn_print_iter(f / fscale, gfull, niter, S$nfeval, pivot)
    if (!newcon) {
      P$yk <- g - oldg
      P$sk <- alpha * pk
      P$yksk <- .sptn_dot(P$yk, P$sk)
      if (icycle == n - 1L || isTRUE(difnew < epsred * (flast_reset - f))) {
        P$lreset <- TRUE
      } else {
        P$yrsr <- .sptn_dot(P$yr, P$sr)
        P$lreset <- isTRUE(P$yrsr <= 0)
      }
      P$upd1 <- FALSE
    }
  }
  if (bitwAnd(messages, 1L) != 0L) .sptn_print_iter(f / fscale, gfull, niter, S$nfeval, pivot)
  x <- .sptn_coerce(.sptn_unscale(x, xscale, xoffset), low, up)
  f <- f / fscale
  done(rc)
}

# _minimize_tnc
.sptn_minimize <- function(fun, x0, jac, bounds, eps, scale, offset, mesg_num, maxCGit, eta, stepmx, accuracy,
                           minfev, ftol, xtol, gtol, rescale, disp, callback, finite_diff_rel_step, maxfun) {
  n <- length(x0)
  low <- rep(-Inf, n)
  up <- rep(Inf, n)
  if (!is.null(bounds)) {
    low <- rep_len(as.numeric(bounds$lb), n)
    up <- rep_len(as.numeric(bounds$ub), n)
  }
  messages <- if (!is.null(mesg_num)) {
    m <- c(0L, 1L, 2L, 4L, 8L, 15L)[as.integer(mesg_num) + 1L]
    if (length(m) != 1L || is.na(m)) 15L else m
  } else if (isTRUE(disp)) 15L else 0L
  scale <- as.numeric(scale %||% numeric())
  offset <- as.numeric(offset %||% numeric())
  if ((length(scale) && length(scale) != n) || (length(offset) && length(offset) != n))
    stop("tnc: vector sizes must be equal", call. = FALSE)
  sf <- .sp_prepare_scalar_function(fun, x0, jac = jac, lb = low, ub = up, epsilon = eps,
                                    finite_diff_rel_step = finite_diff_rel_step)
  if (is.null(maxfun)) maxfun <- max(100L, 10L * n)
  res <- .sptn_tnc(sf$fun_and_grad, x0, low, up, scale, offset, messages, as.integer(maxCGit), as.integer(maxfun),
                   eta, stepmx, accuracy, minfev, ftol, xtol, gtol, rescale, callback)
  # x, f and g of tnc.c may be slightly out of sync through the scaling, so
  # SciPy evaluates them again
  fg <- sf$fun_and_grad(res$x)
  rc <- res$rc
  list(x = res$x, fun = fg$f, jac = fg$g, nfev = sf$nfev(), nit = res$niter, status = rc,
       message = .sptn_messages[rc + 2L], success = rc > -1L && rc < 3L)
}

# SciPy's minimize(method = "TNC"): the variables fixed by equal bounds are
# always removed (_minimize, gh-14565), then _minimize_tnc runs on the rest
.sp_tnc <- function(fun, x0, jac = NULL, bounds = NULL, eps = 1e-8, scale = NULL, offset = NULL, mesg_num = NULL,
                    maxCGit = -1L, eta = -1, stepmx = 0, accuracy = 0, minfev = 0, ftol = -1, xtol = -1, gtol = -1,
                    rescale = -1, disp = FALSE, callback = NULL, finite_diff_rel_step = NULL, maxfun = NULL,
                    workers = NULL, ...) {
  extra <- names(list(...))
  if (length(extra)) warning("Unknown solver options: ", paste(extra, collapse = ", "), call. = FALSE)
  x0 <- as.numeric(x0)
  run <- function(fun, x0, jac, bounds, callback) {
    .sptn_minimize(fun, x0, jac, bounds, eps, scale, offset, mesg_num, maxCGit, eta, stepmx, accuracy, minfev, ftol,
                   xtol, gtol, rescale, disp, callback, finite_diff_rel_step, maxfun)
  }
  if (!is.null(bounds)) {
    lb <- rep_len(as.numeric(bounds$lb), length(x0))
    ub <- rep_len(as.numeric(bounds$ub), length(x0))
    i_fixed <- lb == ub
    if (any(i_fixed)) {
      x_fixed <- lb[i_fixed]
      full <- function(x) {
        out <- numeric(length(i_fixed))
        out[i_fixed] <- x_fixed
        out[!i_fixed] <- x
        out
      }
      fun0 <- fun
      jac0 <- jac
      cb0 <- callback
      res <- run(function(x) fun0(full(x)), x0[!i_fixed],
                 if (is.function(jac0)) function(x) as.numeric(jac0(full(x)))[!i_fixed] else jac0,
                 list(lb = lb[!i_fixed], ub = ub[!i_fixed]), if (is.function(cb0)) function(x) cb0(full(x)))
      res$x <- full(res$x)
      g <- rep(NaN, length(i_fixed))
      g[!i_fixed] <- res$jac
      res$jac <- g
      return(res)
    }
  }
  run(fun, x0, jac, bounds, callback)
}
