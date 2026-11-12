# The nonlinear program solver of the native dynamic optimization backend: a
# primal dual interior point method for
#
#   min f(z)  subject to  cE(z) = 0,  cI(z) <= 0,  lo <= z <= hi,
#
# after Waechter and Biegler, "On the implementation of an interior-point
# filter line-search algorithm for large-scale nonlinear programming" (2006),
# the algorithm of Ipopt, which the JuMP, InfiniteOpt, CasADi and Pyomo
# backends of ModelingToolkit (ext/MTKInfiniteOptExt.jl and the others) call.
# Inequalities get slacks s >= 0 and bounds a logarithmic barrier with
# primal dual bound multipliers (equations 13 to 16 of the paper). The Newton
# step solves the KKT system with the inertia free regularization of Chiang
# and Zavala (a curvature test along the step), steps stay inside the bounds
# by the fraction to the boundary rule, and the filter line search with
# second order corrections of section 2 of the paper accepts them; the
# barrier parameter follows the monotone update of Ipopt and the problem is
# scaled by gradients at the start as Ipopt's default. The restoration phase
# of Ipopt is replaced by a minimum norm step toward feasibility.
#
# The KKT matrix of a collocation problem is banded when unknowns and rows
# are ordered in time: `keys` gives the time position of every unknown and
# row (Inf for global ones such as tuned parameters), the banded part is
# factorized by Gaussian elimination with partial pivoting inside the band,
# and the global unknowns enter through a Schur complement.

# Solve K x = B for a matrix K whose nonzeros lie within `bw` of the
# diagonal; partial pivoting widens the upper band to 2 bw.
.nlp_band_solve <- function(K, B, bw) {
  n <- nrow(K)
  B <- as.matrix(B)
  if (!n) return(B)
  bu <- 2L * bw
  for (j in seq_len(n)) {
    rl <- min(n, j + bw)
    cu <- min(n, j + bu)
    if (rl > j) {
      cand <- j:rl
      p <- cand[which.max(abs(K[cand, j]))]
      if (p != j) {
        cols <- j:cu
        tmp <- K[j, cols]; K[j, cols] <- K[p, cols]; K[p, cols] <- tmp
        tmp <- B[j, ]; B[j, ] <- B[p, ]; B[p, ] <- tmp
      }
    }
    if (K[j, j] == 0) K[j, j] <- .Machine$double.eps * max(1, abs(K[j, j:cu]))
    if (rl > j) {
      r <- (j + 1L):rl
      l <- K[r, j] / K[j, j]
      nzl <- l != 0
      if (any(nzl)) {
        r <- r[nzl]; l <- l[nzl]
        if (cu > j) {
          cols <- (j + 1L):cu
          K[r, cols] <- K[r, cols, drop = FALSE] - outer(l, K[j, cols])
        }
        B[r, ] <- B[r, , drop = FALSE] - outer(l, B[j, ])
      }
    }
  }
  for (i in n:1) {
    cu <- min(n, i + bu)
    if (cu > i) {
      cols <- (i + 1L):cu
      B[i, ] <- B[i, ] - colSums(K[i, cols] * B[cols, , drop = FALSE])
    }
    B[i, ] <- B[i, ] / K[i, i]
  }
  B
}

# Solve the KKT system given as triplets (i, j, v) of a square matrix of
# dimension `dim` (duplicates are summed) and right-hand side r. `key` is
# the time position of each index; indices with key Inf form the border.
.nlp_kkt_solve <- function(ti, tj, tv, dim, r, key) {
  ord <- order(key, seq_len(dim))
  pos <- integer(dim)
  pos[ord] <- seq_len(dim)
  border <- is.infinite(key)
  nb <- sum(border)
  nc <- dim - nb
  pi <- pos[ti]; pj <- pos[tj]
  core <- pi <= nc & pj <= nc
  bw <- if (any(core)) max(abs(pi[core] - pj[core]), 1L) else 1L
  K <- matrix(0, dim, dim)
  lin <- (pj - 1) * dim + pi
  s <- rowsum(tv, lin)
  K[as.numeric(rownames(s))] <- s[, 1L]
  rr <- numeric(dim)
  rr[pos] <- r
  if (!nb) {
    x <- .nlp_band_solve(K, rr, bw)[, 1L]
  } else {
    ci <- seq_len(nc); bi <- nc + seq_len(nb)
    Y <- .nlp_band_solve(K[ci, ci, drop = FALSE], cbind(rr[ci], K[ci, bi, drop = FALSE]), bw)
    x1 <- Y[, 1L]
    Yb <- Y[, -1L, drop = FALSE]
    S <- K[bi, bi, drop = FALSE] - K[bi, ci, drop = FALSE] %*% Yb
    x2 <- tryCatch(solve(S, rr[bi] - K[bi, ci, drop = FALSE] %*% x1), error = function(e) NULL)
    if (is.null(x2)) x2 <- qr.solve(S + diag(1e-12 * max(1, abs(S)), nb), rr[bi] - K[bi, ci, drop = FALSE] %*% x1)
    x <- c(x1 - as.numeric(Yb %*% x2), as.numeric(x2))
  }
  x[pos]
}

# sum(v[idx == k]) for k in 1:n.
.nlp_accum <- function(idx, v, n) {
  out <- numeric(n)
  if (length(idx)) {
    a <- rowsum(v, idx)
    out[as.integer(rownames(a))] <- a[, 1L]
  }
  out
}

# Start strictly inside the bounds (Ipopt's bound_push and bound_frac).
.nlp_push_inside <- function(z, lo, hi, k1 = 1e-2, k2 = 1e-2) {
  pl <- ifelse(is.finite(lo), pmin(k1 * pmax(1, abs(lo)), ifelse(is.finite(hi), k2 * (hi - lo), Inf)), 0)
  pu <- ifelse(is.finite(hi), pmin(k1 * pmax(1, abs(hi)), ifelse(is.finite(lo), k2 * (hi - lo), Inf)), 0)
  z <- ifelse(is.finite(lo), pmax(z, lo + pl), z)
  ifelse(is.finite(hi), pmin(z, hi - pu), z)
}

# The largest step in (0, 1] keeping x + a dx at least (1 - tau) of the way
# from its bound (the fraction to the boundary rule).
.nlp_ftb <- function(gap, dgap, tau) {
  neg <- dgap < 0 & is.finite(gap)
  if (!any(neg)) return(1)
  min(1, min(-tau * gap[neg] / dgap[neg]))
}

.nlp_interior_point <- function(evalf, z0, lo, hi, var_key, rowE_key, rowI_key, tol = 1e-8, maxiter = 300L,
                                verbose = FALSE) {
  nz <- length(z0)
  # a fixed unknown keeps a tiny interval so the barrier is defined
  fixed <- is.finite(lo) & is.finite(hi) & hi - lo < 1e-12 * (1 + abs(lo))
  lo[fixed] <- lo[fixed] - 1e-9 * (1 + abs(lo[fixed]))
  hi[fixed] <- hi[fixed] + 1e-9 * (1 + abs(hi[fixed]))
  z <- .nlp_push_inside(z0, lo, hi)
  # gradient based scaling of the objective and of each constraint so that
  # no gradient entry exceeds 100 at the start (nlp_scaling_method =
  # gradient-based)
  e1 <- evalf(z, 1L)
  nE <- length(e1$cE); nI <- length(e1$cI)
  rowmax <- function(J, nr) {
    out <- numeric(nr)
    if (length(J$i)) {
      a <- tapply(abs(J$v), J$i, max)
      out[as.integer(names(a))] <- a
    }
    out
  }
  gmax <- max(abs(e1$gf), 0)
  sf <- if (gmax > 100) 100 / gmax else 1
  sE <- pmin(1, 100 / pmax(rowmax(e1$JE, nE), 1e-300))
  sI <- pmin(1, 100 / pmax(rowmax(e1$JI, nI), 1e-300))
  sev <- function(z, level, lamE = NULL, lamI = NULL) {
    e <- evalf(z, level, if (!is.null(lamE)) lamE * sE, if (!is.null(lamI)) lamI * sI, sf)
    e$f <- sf * e$f; e$cE <- sE * e$cE; e$cI <- sI * e$cI
    if (level >= 1L) {
      e$gf <- sf * e$gf
      e$JE$v <- sE[e$JE$i] * e$JE$v
      e$JI$v <- sI[e$JI$i] * e$JI$v
    }
    e
  }
  ev <- sev(z, 0L)
  m <- nE + nI
  ny <- nz + nI
  ylo <- c(lo, rep(0, nI)); yhi <- c(hi, rep(Inf, nI))
  hl <- is.finite(ylo); hu <- is.finite(yhi)
  s <- .nlp_push_inside(pmax(-ev$cI, 0), rep(0, nI), rep(Inf, nI))
  y <- c(z, s)
  mu <- 0.1
  lam <- numeric(m)
  zl <- ifelse(hl, 1, 0); zu <- ifelse(hu, 1, 0)
  dw_last <- 0
  dim <- ny + m
  key <- c(var_key, rowI_key, rowE_key, rowI_key)
  status <- "MaxIters"
  barrier <- function(yy) {
    a <- if (any(hl)) sum(log(yy[hl] - ylo[hl])) else 0
    b <- if (any(hu)) sum(log(yhi[hu] - yy[hu])) else 0
    -mu * (a + b)
  }
  cons <- function(e, yy) c(e$cE, e$cI + yy[nz + seq_len(nI)])
  # the filter of (constraint violation, barrier objective) pairs
  theta0 <- sum(abs(cons(ev, y)))
  th_max <- 1e4 * max(1, theta0); th_min <- 1e-4 * max(1, theta0)
  filt <- matrix(numeric(), 0L, 2L)
  gth <- 1e-5; gph <- 1e-8; sth <- 1.1; sph <- 2.3; eta <- 1e-4
  acc <- 0L
  fails <- 0L
  # the watchdog: after ten shortened steps in a row the full step is taken
  # and kept when one of the next three iterates is acceptable to the
  # reference point (watchdog_shortened_iter_trigger = 10,
  # watchdog_trial_iter_max = 3); otherwise the reference point is restored
  shortened <- 0L
  wd <- NULL
  it <- 0L
  for (it in seq_len(maxiter)) {
    z <- y[seq_len(nz)]
    ev <- sev(z, 2L, lam[seq_len(nE)], lam[nE + seq_len(nI)])
    c0 <- cons(ev, y)
    # gradient of the Lagrangian with respect to y
    gl <- c(ev$gf, numeric(nI))
    JT <- c(.nlp_accum(ev$JE$j, lam[ev$JE$i] * ev$JE$v, nz) + .nlp_accum(ev$JI$j, lam[nE + ev$JI$i] * ev$JI$v, nz),
            lam[nE + seq_len(nI)])
    gapl <- ifelse(hl, y - ylo, Inf); gapu <- ifelse(hu, yhi - y, Inf)
    dual <- gl + JT - zl + zu
    sd <- max(100, (sum(abs(lam)) + sum(abs(zl)) + sum(abs(zu))) / max(1, m + 2 * ny)) / 100
    sc <- max(100, (sum(abs(zl)) + sum(abs(zu))) / max(1, 2 * ny)) / 100
    compl <- max(0, abs((gapl * zl)[hl]), abs((gapu * zu)[hu]))
    err0 <- max(max(abs(dual), 0) / sd, max(abs(c0), 0), compl / sc)
    if (is.finite(err0) && err0 <= tol) {
      status <- "Success"
      break
    }
    # an acceptable point fifteen times in a row also ends the iteration
    acc <- if (is.finite(err0) && err0 <= 1e-6) acc + 1L else 0L
    if (acc >= 15L) {
      status <- "Success"
      break
    }
    # barrier parameter update (monotone Fiacco McCormick); a new barrier
    # problem starts with an empty filter
    repeat {
      complm <- max(0, abs((gapl * zl)[hl] - mu), abs((gapu * zu)[hu] - mu))
      errm <- max(max(abs(dual), 0) / sd, max(abs(c0), 0), complm / sc)
      if (errm > 10 * mu || mu <= tol / 10) break
      mu <- max(tol / 10, min(0.2 * mu, mu^1.5))
      filt <- matrix(numeric(), 0L, 2L)
    }
    # the Newton system: [W + Sigma, J'; J, -dc] [dy; dlam] = -[grad phi + J' lam; c]
    sig <- ifelse(hl, zl / gapl, 0) + ifelse(hu, zu / gapu, 0)
    gphi <- gl - ifelse(hl, mu / gapl, 0) + ifelse(hu, mu / gapu, 0)
    rhs <- -c(gphi + JT, c0)
    Hi <- ev$H$i; Hj <- ev$H$j; Hv <- ev$H$v
    ji <- c(ev$JE$i, nE + ev$JI$i, nE + seq_len(nI))
    jj <- c(ev$JE$j, ev$JI$j, nz + seq_len(nI))
    jv <- c(ev$JE$v, ev$JI$v, rep(1, nI))
    dc <- 1e-9
    dw <- 0
    for (attempt in 1:25) {
      ti <- c(Hi, seq_len(ny), ny + ji, jj)
      tj <- c(Hj, seq_len(ny), jj, ny + ji)
      tv <- c(Hv, sig + dw, jv, jv)
      if (m) {
        ti <- c(ti, ny + seq_len(m)); tj <- c(tj, ny + seq_len(m)); tv <- c(tv, rep(-dc, m))
      }
      d <- tryCatch(.nlp_kkt_solve(ti, tj, tv, dim, rhs, key), error = function(e) NULL)
      ok <- !is.null(d) && all(is.finite(d))
      if (ok) {
        dy <- d[seq_len(ny)]
        hdy <- .nlp_accum(Hi, Hv * dy[Hj], ny)
        curv <- sum(dy * hdy) + sum((sig + dw) * dy^2)
        if (curv >= 1e-10 * sum(dy^2)) break
      }
      dw <- if (dw == 0) (if (dw_last == 0) 1e-4 else max(1e-20, dw_last / 3)) else 8 * dw
      if (dw > 1e40) break
    }
    if (dw > 0) dw_last <- dw
    if (is.null(d) || !all(is.finite(d))) {
      status <- "Failure"
      break
    }
    solve_kkt <- function(r) .nlp_kkt_solve(ti, tj, tv, dim, r, key)
    tau <- max(0.99, 1 - mu)
    step_of <- function(d) {
      dy <- d[seq_len(ny)]
      list(dy = dy, dlam = d[ny + seq_len(m)],
           dzl = ifelse(hl, mu / gapl - zl - zl / gapl * dy, 0),
           dzu = ifelse(hu, mu / gapu - zu + zu / gapu * dy, 0),
           amax = min(.nlp_ftb(gapl, dy, tau), .nlp_ftb(gapu, -dy, tau)))
    }
    st <- step_of(d)
    th0 <- sum(abs(c0))
    phi0 <- ev$f + barrier(y)
    D <- sum(gphi * st$dy)
    amin <- 0.05 * if (D < 0) {
      if (th0 <= th_min) min(gth, gph * th0 / (-D), th0^sth / (-D)^sph) else min(gth, gph * th0 / (-D))
    } else gth
    # the filter acceptance test: "f" for an Armijo step on the barrier
    # objective, "h" for a step that is accepted by the filter
    judge <- function(tht, pht, a) {
      if (!is.finite(tht) || !is.finite(pht) || tht > th_max) return(NA_character_)
      if (nrow(filt) && any(tht >= filt[, 1L] & pht >= filt[, 2L])) return(NA_character_)
      if (D < 0 && a * (-D)^sph > th0^sth && th0 <= th_min)
        return(if (pht <= phi0 + eta * a * D) "f" else NA_character_)
      if (tht <= (1 - gth) * th0 || pht <= phi0 - gph * th0) "h" else NA_character_
    }
    trial <- function(st, a) {
      yt <- y + a * st$dy
      et <- tryCatch(sev(yt[seq_len(nz)], 0L), error = function(e) NULL)
      if (is.null(et) || !is.finite(et$f)) return(list(y = yt, th = Inf, ph = Inf, c = NULL))
      ct <- cons(et, yt)
      list(y = yt, th = sum(abs(ct)), ph = et$f + barrier(yt), c = ct)
    }
    if (!is.null(wd)) {
      # judge the full step against the reference point of the watchdog
      tr <- trial(st, st$amax)
      ok <- is.finite(tr$th) && tr$th <= th_max &&
        !(nrow(wd$filt) && any(tr$th >= wd$filt[, 1L] & tr$ph >= wd$filt[, 2L])) &&
        (tr$th <= (1 - gth) * wd$th0 || tr$ph <= wd$phi0 - gph * wd$th0)
      wd$n <- wd$n + 1L
      if (!ok && wd$n >= 3L) {
        y <- wd$y; lam <- wd$lam; zl <- wd$zl; zu <- wd$zu; filt <- wd$filt
        wd <- NULL
        shortened <- 0L
        next
      }
      if (ok) {
        filt <- rbind(wd$filt, c((1 - gth) * wd$th0, wd$phi0 - gph * wd$th0))
        wd <- NULL
        shortened <- 0L
      }
      a <- st$amax
      az <- min(.nlp_ftb(zl, st$dzl, tau), .nlp_ftb(zu, st$dzu, tau))
      y <- tr$y
      lam <- lam + a * st$dlam
      zl <- zl + az * st$dzl; zu <- zu + az * st$dzu
      gapl <- ifelse(hl, y - ylo, Inf); gapu <- ifelse(hu, yhi - y, Inf)
      zl <- ifelse(hl, pmin(pmax(zl, mu / (1e10 * gapl)), 1e10 * mu / gapl), 0)
      zu <- ifelse(hu, pmin(pmax(zu, mu / (1e10 * gapu)), 1e10 * mu / gapu), 0)
      if (verbose) message(sprintf("iter %3d  f %.8g  inf_pr %.2e  inf_du %.2e  mu %.1e  watchdog step",
                                   it, ev$f / sf, max(abs(c0), 0), max(abs(dual), 0), mu))
      next
    }
    if (shortened >= 10L && st$amax > 0) {
      wd <- list(y = y, lam = lam, zl = zl, zu = zu, filt = filt, th0 = th0, phi0 = phi0, n = 0L)
      shortened <- 0L
      tr <- trial(st, st$amax)
      if (is.finite(tr$th)) {
        a <- st$amax
        az <- min(.nlp_ftb(zl, st$dzl, tau), .nlp_ftb(zu, st$dzu, tau))
        y <- tr$y
        lam <- lam + a * st$dlam
        zl <- zl + az * st$dzl; zu <- zu + az * st$dzu
        gapl <- ifelse(hl, y - ylo, Inf); gapu <- ifelse(hu, yhi - y, Inf)
        zl <- ifelse(hl, pmin(pmax(zl, mu / (1e10 * gapl)), 1e10 * mu / gapl), 0)
        zu <- ifelse(hu, pmin(pmax(zu, mu / (1e10 * gapu)), 1e10 * mu / gapu), 0)
        if (verbose) message(sprintf("iter %3d  f %.8g  inf_pr %.2e  inf_du %.2e  mu %.1e  watchdog start",
                                     it, ev$f / sf, max(abs(c0), 0), max(abs(dual), 0), mu))
        next
      }
      wd <- NULL
    }
    a <- st$amax
    kind <- NA_character_
    first <- TRUE
    repeat {
      tr <- trial(st, a)
      kind <- judge(tr$th, tr$ph, a)
      if (!is.na(kind)) break
      if (first && tr$th >= th0 && !is.null(tr$c)) {
        # second order corrections against the Maratos effect
        csoc <- a * c0 + tr$c
        th_old <- th0
        for (p in 1:4) {
          ds <- tryCatch(solve_kkt(-c(gphi + JT, csoc)), error = function(e) NULL)
          if (is.null(ds) || !all(is.finite(ds))) break
          ss <- step_of(ds)
          ts <- trial(ss, ss$amax)
          kind <- judge(ts$th, ts$ph, a)
          if (!is.na(kind)) {
            st <- ss; tr <- ts; a <- ss$amax
            break
          }
          if (is.null(ts$c) || ts$th > 0.99 * th_old) break
          th_old <- ts$th
          csoc <- ss$amax * csoc + ts$c
        }
        if (!is.na(kind)) break
      }
      first <- FALSE
      a <- a / 2
      if (a < amin) break
    }
    if (is.na(kind)) {
      # the line search failed: a step toward feasibility (in place of the
      # restoration phase of Ipopt), after adding the iterate to the filter
      filt <- rbind(filt, c((1 - gth) * th0, phi0 - gph * th0))
      fails <- fails + 1L
      if (fails >= 10L) {
        status <- "Stalled"
        break
      }
      dfz <- tryCatch(solve_kkt(c(numeric(ny), -c0)), error = function(e) NULL)
      if (is.null(dfz) || !all(is.finite(dfz))) {
        status <- "Stalled"
        break
      }
      st <- step_of(dfz)
      a <- st$amax
      repeat {
        tr <- trial(st, a)
        if (tr$th <= (1 - 1e-4 * a) * th0 || a < 1e-10) break
        a <- a / 2
      }
      if (a < 1e-10) {
        status <- "Stalled"
        break
      }
      dw_last <- max(1e-4, 10 * dw_last)
    } else {
      fails <- 0L
      if (kind == "h") filt <- rbind(filt, c((1 - gth) * th0, phi0 - gph * th0))
    }
    shortened <- if (a < st$amax) shortened + 1L else 0L
    if (verbose) message(sprintf("iter %3d  f %.8g  inf_pr %.2e  inf_du %.2e  mu %.1e  dw %.1e  alpha %.2e %s",
                                 it, ev$f / sf, max(abs(c0), 0), max(abs(dual), 0), mu, dw, a,
                                 if (is.na(kind)) "r" else kind))
    az <- min(.nlp_ftb(zl, st$dzl, tau), .nlp_ftb(zu, st$dzu, tau))
    y <- tr$y
    lam <- lam + a * st$dlam
    zl <- zl + az * st$dzl
    zu <- zu + az * st$dzu
    gapl <- ifelse(hl, y - ylo, Inf); gapu <- ifelse(hu, yhi - y, Inf)
    # keep the bound multipliers near the central path (kappa_Sigma = 1e10)
    zl <- ifelse(hl, pmin(pmax(zl, mu / (1e10 * gapl)), 1e10 * mu / gapl), 0)
    zu <- ifelse(hu, pmin(pmax(zu, mu / (1e10 * gapu)), 1e10 * mu / gapu), 0)
  }
  z <- y[seq_len(nz)]
  ev <- evalf(z, 0L)
  viol <- max(0, abs(ev$cE), ev$cI)
  if (status != "Success" && viol > sqrt(tol)) status <- "Infeasible"
  list(z = z, f = ev$f, lamE = lam[seq_len(nE)] * sE / sf, lamI = lam[nE + seq_len(nI)] * sI / sf,
       violation = viol, iterations = it, status = status)
}
