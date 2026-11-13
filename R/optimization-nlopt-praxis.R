# Brent's principal axis method of NLopt 2.11.0 (LN_PRAXIS,
# src/algs/praxis/praxis.c, the f2c translation of Brent's PRAXIS modified
# by S. G. Johnson for the NLopt stopping criteria), ported to R, with the
# bound handling of the LN_PRAXIS case of nlopt_optimize_ (optimize.c): the
# bounded coordinates are mapped to an unbounded space by tanh (two finite
# bounds) or a square (one finite bound).
#
# Copyright (c) 2007-2014 Massachusetts Institute of Technology, released
# under the MIT license (see R/optimization-nlopt-util.R). PRAXIS is by
# Richard P. Brent, from "Algorithms for Minimization without Derivatives"
# (1973), as distributed with NLopt.

# x_bound_inv and x_bound of optimize.c
.nlopt_x_bound_inv <- function(x, lb, ub) {
  for (i in seq_along(x)) {
    if (!.nlopt_isinf(lb[i]) && !.nlopt_isinf(ub[i])) {
      mid <- (lb[i] + ub[i]) * 0.5
      width <- (ub[i] - lb[i]) * 0.5
      v <- (x[i] - mid) / width
      if (v <= -1.0) v <- -1.0 + 1e-12
      if (v >= 1.0) v <- 1.0 - 1e-12
      x[i] <- atanh(v)
    } else if (!.nlopt_isinf(lb[i])) {
      x[i] <- if (x[i] > lb[i]) sqrt(x[i] - lb[i]) else 0.0
    } else if (!.nlopt_isinf(ub[i])) {
      x[i] <- if (x[i] < ub[i]) sqrt(ub[i] - x[i]) else 0.0
    }
  }
  x
}
.nlopt_x_bound <- function(x, lb, ub) {
  for (i in seq_along(x)) {
    if (!.nlopt_isinf(lb[i]) && !.nlopt_isinf(ub[i])) {
      mid <- (lb[i] + ub[i]) * 0.5
      width <- (ub[i] - lb[i]) * 0.5
      x[i] <- mid + tanh(x[i]) * width
    } else if (!.nlopt_isinf(lb[i])) {
      x[i] <- lb[i] + x[i] * x[i]
    } else if (!.nlopt_isinf(ub[i])) {
      x[i] <- ub[i] - x[i] * x[i]
    }
  }
  x
}

# minfit_: the singular value decomposition of ab (n by n) by Golub and
# Reinsch; returns list(ab, q, e) with `e` the work vector (which praxis
# shares with its random step vector z)
.nlopt_praxis_minfit <- function(n, machep, tol, ab, q, e) {
  if (n == 1L) {
    q[1L] <- ab[1L, 1L]
    ab[1L, 1L] <- 1
    return(list(ab = ab, q = q, e = e))
  }
  eps <- machep
  g <- 0; x <- 0; f <- 0; h <- 0
  for (i in seq_len(n)) {
    e[i] <- g
    s <- 0
    l <- i + 1L
    for (j in i:n) s <- s + ab[j, i] * ab[j, i]
    g <- 0
    if (!(s < tol)) {
      f <- ab[i, i]
      g <- sqrt(s)
      if (f >= 0) g <- -g
      h <- f * g - s
      ab[i, i] <- f - g
      if (l <= n) {
        for (j in l:n) {
          f <- 0
          for (k in i:n) f <- f + ab[k, i] * ab[k, j]
          f <- f / h
          for (k in i:n) ab[k, j] <- ab[k, j] + f * ab[k, i]
        }
      }
    }
    q[i] <- g
    s <- 0
    if (i != n) for (j in l:n) s <- s + ab[i, j] * ab[i, j]
    g <- 0
    if (!(s < tol)) {
      if (i != n) f <- ab[i, i + 1L]
      g <- sqrt(s)
      if (f >= 0) g <- -g
      h <- f * g - s
      if (i != n) {
        ab[i, i + 1L] <- f - g
        for (j in l:n) e[j] <- ab[i, j] / h
        for (j in l:n) {
          s <- 0
          for (k in l:n) s <- s + ab[j, k] * ab[i, k]
          for (k in l:n) ab[j, k] <- ab[j, k] + s * e[k]
        }
      }
    }
    y <- abs(q[i]) + abs(e[i])
    if (y > x) x <- y
  }
  # accumulation of the right hand transformations
  ab[n, n] <- 1
  g <- e[n]
  l <- n
  for (ii in 2L:n) {
    i <- n - ii + 1L
    if (g != 0) {
      h <- ab[i, i + 1L] * g
      for (j in l:n) ab[j, i] <- ab[i, j] / h
      for (j in l:n) {
        s <- 0
        for (k in l:n) s <- s + ab[i, k] * ab[k, j]
        for (k in l:n) ab[k, j] <- ab[k, j] + s * ab[k, i]
      }
    }
    for (j in l:n) {
      ab[i, j] <- 0
      ab[j, i] <- 0
    }
    ab[i, i] <- 1
    g <- e[i]
    l <- i
  }
  # diagonalization of the bidiagonal form
  eps <- eps * x
  hyp <- function(a, b) { d1 <- a / b; abs(b) * sqrt(d1 * d1 + 1) }
  for (kk in seq_len(n)) {
    k <- n - kk + 1L
    kt <- 0L
    repeat {
      kt <- kt + 1L
      if (kt > 30L) e[k] <- 0
      converged <- FALSE
      for (ll2 in seq_len(k)) {
        l <- k - ll2 + 1L
        if (abs(e[l]) <= eps) { converged <- TRUE; break }
        if (l == 1L) next
        if (abs(q[l - 1L]) <= eps) break
      }
      if (!converged) {
        # cancellation of e(l) if l > 1
        cc <- 0; s <- 1
        for (i in l:k) {
          f <- s * e[i]
          e[i] <- cc * e[i]
          if (abs(f) <= eps) break
          g <- q[i]
          h <- if (abs(f) < abs(g)) hyp(f, g) else if (f != 0) hyp(g, f) else 0
          q[i] <- h
          if (h == 0) { g <- 1; h <- 1 }
          cc <- g / h
          s <- -f / h
        }
      }
      z <- q[k]
      if (l == k) break
      # shift from the bottom 2 by 2 minor
      x <- q[l]; y <- q[k - 1L]; g <- e[k - 1L]; h <- e[k]
      f <- ((y - z) * (y + z) + (g - h) * (g + h)) / (h * 2 * y)
      g <- sqrt(f * f + 1)
      temp <- f - g
      if (f >= 0) temp <- f + g
      f <- ((x - z) * (x + z) + h * (y / temp - h)) / x
      # next QR transformation
      cc <- 1; s <- 1
      if (l + 1L <= k) {
        for (i in (l + 1L):k) {
          g <- e[i]
          y <- q[i]
          h <- s * g
          g <- g * cc
          z <- if (abs(f) < abs(h)) hyp(f, h) else if (f != 0) hyp(h, f) else 0
          e[i - 1L] <- z
          if (z == 0) { f <- 1; z <- 1 }
          cc <- f / z
          s <- h / z
          f <- x * cc + g * s
          g <- -x * s + g * cc
          h <- y * s
          y <- y * cc
          for (j in seq_len(n)) {
            x <- ab[j, i - 1L]
            z <- ab[j, i]
            ab[j, i - 1L] <- x * cc + z * s
            ab[j, i] <- -x * s + z * cc
          }
          z <- if (abs(f) < abs(h)) hyp(f, h) else if (f != 0) hyp(h, f) else 0
          q[i - 1L] <- z
          if (z == 0) { f <- 1; z <- 1 }
          cc <- f / z
          s <- h / z
          f <- cc * g + s * y
          x <- -s * g + cc * y
        }
      }
      e[l] <- 0
      e[k] <- f
      q[k] <- x
    }
    # convergence: q(k) is made non negative
    if (z < 0) {
      q[k] <- -z
      ab[, k] <- -ab[, k]
    }
  }
  list(ab = ab, q = q, e = e)
}

# praxis_: minimizes f (of the unbounded variables) from x with the maximal
# step h0; returns list(ret, x, minf)
.nlopt_praxis <- function(t0, machep, h0, x, f, stop) {
  n <- length(x)
  small <- machep * machep
  vsmall <- small * small
  large <- 1 / small
  vlarge <- 1 / vsmall
  m2 <- sqrt(machep)
  m4 <- sqrt(m2)
  ldfac <- 0.01
  ktm <- 1L
  illc <- FALSE
  kt <- 0L
  nl <- 0L
  nf <- 1L
  fx <- f(x)
  fbest <- fx; prev_fbest <- fx
  xbest <- x; prev_xbest <- x
  stop$nevals <- stop$nevals + 1L
  qf1 <- fx
  if (t0 > 0) {
    t_old <- small + t0
  } else {
    t_old <- 0
    if (!is.null(stop$xtol_abs)) for (i in seq_len(n)) if (stop$xtol_abs[i] > t_old) t_old <- stop$xtol_abs[i]
    t_old <- t_old + small
  }
  t2_old <- t_old
  dmin <- small
  h <- h0
  if (h < t_old * 100) h <- t_old * 100
  ldt <- h
  V <- diag(1, n)
  d <- numeric(n)
  qd0 <- 0; qd1 <- 0; qa <- 0; qb <- 0; qc <- 0
  q0 <- x; q1 <- x
  y <- numeric(n)
  z <- numeric(n) # also minfit's work vector e
  fret <- "SUCCESS"

  flin <- function(j, l) {
    if (j != 0L) {
      t <- x + l * V[, j]
    } else {
      qa <<- l * (l - qd1) / (qd0 * (qd0 + qd1))
      qb <<- (l + qd0) * (qd1 - l) / (qd0 * qd1)
      qc <<- l * (l + qd0) / (qd1 * (qd0 + qd1))
      t <- qa * q0 + qb * x + qc * q1
    }
    nf <<- nf + 1L
    v <- f(t)
    stop$nevals <- stop$nevals + 1L
    if (v < fbest) { fbest <<- v; xbest <<- t }
    if (.nlopt_stop_forced(stop)) fret <<- "FORCED_STOP"
    else if (.nlopt_stop_evals(stop)) fret <<- "MAXEVAL_REACHED"
    else if (.nlopt_stop_time(stop)) fret <<- "MAXTIME_REACHED"
    else if (v <= stop$minf_max) fret <<- "STOPVAL_REACHED"
    v
  }

  # min_: returns list(ret, d2, x1); updates x and fx
  minj <- function(j, nits, d2, x1, f1, fk) {
    fret <<- "SUCCESS"
    out <- function() list(ret = fret, d2 = d2, x1 = x1)
    sf1 <- f1; sx1 <- x1
    k <- 0L
    xm <- 0
    fm <- fx; f0 <- fx
    dz <- d2 < machep
    s <- 0
    for (i in seq_len(n)) s <- s + x[i] * x[i]
    s <- sqrt(s)
    temp <- d2
    if (dz) temp <- dmin
    t2 <- m4 * sqrt(abs(fx) / temp + s * ldt) + m2 * ldt
    s <- m4 * s + t_old
    if (dz && t2 > s) t2 <- s
    t2 <- if (t2 > small) t2 else small
    hh <- h * 0.01
    t2 <- if (t2 < hh) t2 else hh
    if (fk && !(f1 > fm)) { xm <- x1; fm <- f1 }
    if (!(fk && abs(x1) >= t2)) {
      temp <- 1
      if (x1 < 0) temp <- -1
      x1 <- temp * t2
      f1 <- flin(j, x1)
      if (fret != "SUCCESS") return(out())
    }
    if (!(f1 > fm)) { xm <- x1; fm <- f1 }
    at4 <- TRUE
    repeat {
      if (at4) {
        if (dz) {
          x2 <- -x1
          if (f0 >= f1) x2 <- x1 * 2
          f2 <- flin(j, x2)
          if (fret != "SUCCESS") return(out())
          if (!(f2 > fm)) { xm <- x2; fm <- f2 }
          d2 <- (x2 * (f1 - f0) - x1 * (f2 - f0)) / (x1 * x2 * (x1 - x2))
        }
        d1 <- (f1 - f0) / x1 - x1 * d2
        dz <- TRUE
        if (d2 > small) {
          x2 <- d1 * -0.5 / d2
        } else {
          x2 <- h
          if (d1 >= 0) x2 <- -x2
        }
        if (!(abs(x2) <= h)) x2 <- if (x2 <= 0) -h else h
      }
      f2 <- flin(j, x2)
      if (fret != "SUCCESS") return(out())
      if (k >= nits || f2 <= f0) break
      k <- k + 1L
      if (f0 < f1 && x1 * x2 > 0) { at4 <- TRUE; next }
      x2 <- x2 * 0.5
      at4 <- FALSE
    }
    nl <<- nl + 1L
    if (f2 <= fm) fm <- f2 else x2 <- xm
    if (abs(x2 * (x2 - x1)) <= small) {
      if (k > 0L) d2 <- 0
    } else {
      d2 <- (x2 * (f1 - f0) - x1 * (fm - f0)) / (x1 * x2 * (x1 - x2))
    }
    if (d2 <= small) d2 <- small
    x1 <- x2
    fx <<- fm
    if (!(sf1 >= fx)) { fx <<- sf1; x1 <- sx1 }
    if (j != 0L) x <<- x + x1 * V[, j]
    list(ret = "SUCCESS", d2 = d2, x1 = x1)
  }

  quad <- function() {
    s <- fx
    fx <<- qf1
    qf1 <<- s
    qd1 <<- 0
    for (i in seq_len(n)) {
      s <- x[i]
      l <- q1[i]
      x[i] <<- l
      q1[i] <<- s
      qd1 <<- qd1 + (s - l) * (s - l)
    }
    qd1 <<- sqrt(qd1)
    l <- qd1
    s <- 0
    if (qd0 <= 0 || qd1 <= 0 || nl < n * 3L * n) {
      fx <<- qf1
      qa <<- 0; qb <<- 0; qc <<- 1
    } else {
      # the return code of this line search is ignored, as in praxis.c
      r <- minj(0L, 2L, s, l, qf1, TRUE)
      l <- r$x1
      qa <<- l * (l - qd1) / (qd0 * (qd0 + qd1))
      qb <<- (l + qd0) * (qd1 - l) / (qd0 * qd1)
      qc <<- l * (l + qd0) / (qd1 * (qd0 + qd1))
    }
    qd0 <<- qd1
    for (i in seq_len(n)) {
      s <- q0[i]
      q0[i] <<- x[i]
      x[i] <<- qa * s + qb * x[i] + qc * q1[i]
    }
  }

  ret <- "SUCCESS"
  repeat { # the main loop (L40)
    sf <- d[1L]
    d[1L] <- 0
    r <- minj(1L, 2L, d[1L], 0, fx, FALSE)
    d[1L] <- r$d2; s <- r$x1
    if (r$ret != "SUCCESS") { ret <- r$ret; break }
    if (!(s > 0)) V[, 1L] <- -V[, 1L]
    if (!(sf > d[1L] * 0.9 && sf * 0.9 < d[1L])) if (n >= 2L) d[2L:n] <- 0
    done <- FALSE
    for (k in seq_len(n)[-1L]) { # the inner loop
      y <- x
      sf <- fx
      if (kt > 0L) illc <- TRUE
      repeat { # L80
        kl <- k
        df <- 0
        if (illc) {
          # a random step (to avoid resolution valleys)
          for (i in seq_len(n)) {
            s <- (ldt * 0.1 + t2_old * 10^kt) * .nlopt_urand(-0.5, 0.5)
            z[i] <- s
            x <- x + s * V[, i]
          }
          fx <- f(x)
          nf <- nf + 1L
          stop$nevals <- stop$nevals + 1L
        }
        for (k2 in k:n) {
          sl <- fx
          r <- minj(k2, 2L, d[k2], 0, fx, FALSE)
          d[k2] <- r$d2; s <- r$x1
          if (r$ret != "SUCCESS") { ret <- r$ret; done <- TRUE; break }
          if (illc) {
            d1 <- s + z[k2]
            s <- d[k2] * (d1 * d1)
          } else s <- sl - fx
          if (!(df > s)) { df <- s; kl <- k2 }
        }
        if (done) break
        if (illc || df >= abs(machep * 100 * fx)) break
        illc <- TRUE
      }
      if (done) break
      # minimize along the conjugate directions v(*,1), ..., v(*,k-1)
      for (k2 in seq_len(k - 1L)) {
        r <- minj(k2, 2L, d[k2], 0, fx, FALSE)
        d[k2] <- r$d2
        if (r$ret != "SUCCESS") { ret <- r$ret; done <- TRUE; break }
      }
      if (done) break
      f1 <- fx
      fx <- sf
      lds <- 0
      for (i in seq_len(n)) {
        sl <- x[i]
        x[i] <- y[i]
        sl <- sl - y[i]
        y[i] <- sl
        lds <- lds + sl * sl
      }
      lds <- sqrt(lds)
      if (!(lds <= small)) {
        # discard direction v(*,kl)
        klmk <- kl - k
        if (klmk >= 1L) {
          for (ii in seq_len(klmk)) {
            i <- kl - ii
            V[, i + 1L] <- V[, i]
            d[i + 1L] <- d[i]
          }
        }
        d[k] <- 0
        V[, k] <- y / lds
        r <- minj(k, 4L, d[k], lds, f1, TRUE)
        d[k] <- r$d2; lds <- r$x1
        if (r$ret != "SUCCESS") { ret <- r$ret; break }
        if (!(lds > 0)) {
          lds <- -lds
          V[, k] <- -V[, k]
        }
      }
      ldt <- ldfac * ldt
      if (ldt < lds) ldt <- lds
      t2_old <- 0
      for (i in seq_len(n)) t2_old <- t2_old + x[i] * x[i]
      t2_old <- m2 * sqrt(t2_old) + t_old
      # the length of the step taken since starting the inner loop
      if (ldt > t2_old * 0.5 && !.nlopt_stop_f(stop, fbest, prev_fbest) &&
          !.nlopt_stop_x(stop, xbest, prev_xbest)) kt <- -1L
      kt <- kt + 1L
      if (kt > ktm) {
        if (.nlopt_stop_f(stop, fbest, prev_fbest)) ret <- "FTOL_REACHED"
        else if (.nlopt_stop_x(stop, xbest, prev_xbest)) ret <- "XTOL_REACHED"
        done <- TRUE
        break
      }
      prev_fbest <- fbest
      prev_xbest <- xbest
    }
    if (done) break
    # quadratic extrapolation in case we are in a curved valley
    quad()
    dn <- 0
    for (i in seq_len(n)) {
      d[i] <- 1 / sqrt(d[i])
      if (dn < d[i]) dn <- d[i]
    }
    for (j in seq_len(n)) {
      s <- d[j] / dn
      V[, j] <- s * V[, j]
    }
    # (scbd = 1: no scaling of the axes) transpose V for minfit
    V <- t(V)
    mf <- .nlopt_praxis_minfit(n, machep, vsmall, V, d, z)
    V <- mf$ab; d <- mf$q; z <- mf$e
    for (i in seq_len(n)) {
      dni <- dn * d[i]
      d[i] <- if (dni > large) vsmall else if (dni < small) vlarge else 1 / (dni * dni)
    }
    # sort the eigenvalues and eigenvectors (descending)
    if (n > 1L) {
      for (i in seq_len(n - 1L)) {
        k <- i
        s <- d[i]
        for (j in (i + 1L):n) if (d[j] > s) { k <- j; s <- d[j] }
        if (k > i) {
          d[k] <- d[i]
          d[i] <- s
          tmp <- V[, i]; V[, i] <- V[, k]; V[, k] <- tmp
        }
      }
    }
    dmin <- d[n]
    if (dmin < small) dmin <- small
    illc <- m2 * d[1L] > dmin
  }
  list(ret = ret, x = xbest, minf = fbest)
}

# the LN_PRAXIS case of nlopt_optimize_
.nlopt_praxis_dispatch <- function(o, x, stop) {
  dx <- .nlopt_dx(o, x)
  step <- Inf
  for (i in seq_along(dx)) if (step > abs(dx[i])) step <- abs(dx[i])
  lb <- o$lb; ub <- o$ub
  xu <- .nlopt_x_bound_inv(x, lb, ub)
  fb <- function(xx) o$f(.nlopt_x_bound(xx, lb, ub))
  res <- .nlopt_praxis(.nlopt_get_param(o, "t0_tol", 0), .Machine$double.eps, step, xu, fb, stop)
  res$x <- .nlopt_x_bound(res$x, lb, ub)
  res
}
