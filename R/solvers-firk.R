# Fully implicit Runge-Kutta methods of OrdinaryDiffEqFIRK (algorithms.jl,
# alg_utils.jl, controllers.jl, firk_caches.jl, firk_perform_step.jl,
# firk_interpolants.jl, firk_tableaus.jl): RadauIIA3, RadauIIA5, RadauIIA9,
# the adaptive order AdaptiveRadau and GaussLegendre, on the out of place
# ConstantCache paths.
#
# The Radau methods solve the stage system in the eigenbasis of the inverse
# Runge-Kutta matrix (Hairer and Wanner's transformation): w = TI z
# decouples it into one real block J - gamma / dt M and complex blocks
# J - (alpha + i beta) / dt M, each factorized once per step (a fresh
# Jacobian every step, as Julia's out of place code). The simplified Newton
# iteration stops on the rate estimate eta * ||dw|| < kappa; the error
# estimate is Hairer's (smoothed through the real block, with a second
# evaluation after a first rejected or discontinuous step). AdaptiveRadau
# moves between 3, 5, 7, ... stages by the average Newton iteration count.
# GaussLegendre solves the full s n by s n stage system with the PI
# controller and an embedded error estimate of order s - 1.
#
# The stage arithmetic follows Julia's expression shapes under @muladd: the
# broadcast sums of products are fused multiply adds (.fma), the vector
# sums outside broadcasts are not. Linear systems are solved with LAPACK's
# LU (solve()), the factorization Julia's lu uses; tol = 0 drops R's
# reciprocal condition check, so only an exactly singular matrix fails as
# in lu(A) with check = true.

# calculate_residuals with internalnorm (DiffEqBase): the scalar residual
# err / (abstol + max(|u0|, |u1|) * reltol) is a fused multiply add under
# @muladd; the norm is the root mean square.
.firk_err <- function(err, u0, u1, abstol, reltol)
  .rms_norm(err / .fma(pmax(abs(u0), abs(u1)), reltol, abstol))

# A broadcast sum of products A[i, 1] * v1 + A[i, 2] * v2 + ... under
# @muladd: the first product, then one fused multiply add per term.
.firk_chain <- function(A, i, v) {
  acc <- A[i, 1L] * v[[1L]]
  for (j in seq_along(v)[-1L]) acc <- .fma(A[i, j], v[[j]], acc)
  acc
}

# ---------------------------------------------------------------------------
# Tableaus

# The square matrix nm11, nm12, ... of a dumped tableau; `extra` holds the
# entries the hand written methods leave out (T32 = 1, T33 = 0, ...).
.firk_mat <- function(d, nm, s, extra = list()) {
  v <- c(d, extra)
  keys <- outer(seq_len(s), seq_len(s), function(i, j) paste0(nm, i, j))
  matrix(vapply(keys, function(k) v[[k]], 0), s, s)
}

# RadauIIA3Tableau, RadauIIA5Tableau and RadauIIA9Tableau (firk_tableaus.jl)
# as matrices, with c ending in 1.
.firk_iia_tab <- function(name) {
  d <- .FIRK_TAB[[name]]
  switch(name,
    RadauIIA3 = list(T = .firk_mat(d, "T", 2L), TI = .firk_mat(d, "TI", 2L), c = c(d$c1, d$c2),
                     alpha = d$alpha, beta = d$beta, e = c(d$e1, d$e2)),
    RadauIIA5 = list(T = .firk_mat(d, "T", 3L, list(T32 = 1, T33 = 0)), TI = .firk_mat(d, "TI", 3L),
                     c = c(d$c1, d$c2, 1), gamma = d$gamma, alpha = d$alpha, beta = d$beta,
                     e = c(d$e1, d$e2, d$e3)),
    RadauIIA9 = list(T = .firk_mat(d, "T", 5L, list(T52 = 1, T53 = 0, T54 = 1, T55 = 0)),
                     TI = .firk_mat(d, "TI", 5L), c = c(d$c1, d$c2, d$c3, d$c4, 1), gamma = d$gamma,
                     alpha = c(d$alpha1, d$alpha2), beta = c(d$beta1, d$beta2),
                     e = c(d$e1, d$e2, d$e3, d$e4, d$e5)))
}

# Nodes and weights of the s point Gauss-Legendre rule on [-1, 1]
# (gausslegendre of FastGaussQuadrature): Golub and Welsch's eigenvalues,
# polished by Newton's method on the Legendre polynomial; the weights are
# 2 / ((1 - x^2) P_s'(x)^2).
.firk_legendre <- function(s, x) {
  # P_s and P_s' at x by the three term recurrence
  p0 <- rep(1, length(x)); p1 <- x
  if (s == 0L) return(list(p = p0, dp = 0 * x))
  for (k in seq_len(s - 1L)) {
    p2 <- ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
    p0 <- p1; p1 <- p2
  }
  list(p = p1, p_prev = p0, dp = s * (x * p1 - p0) / (x * x - 1))
}

.firk_gauss_rule <- function(s) {
  if (s == 1L) return(list(x = 0, w = 2))
  k <- seq_len(s - 1L)
  Jm <- matrix(0, s, s)
  Jm[cbind(k, k + 1L)] <- Jm[cbind(k + 1L, k)] <- k / sqrt(4 * k * k - 1)
  x <- sort(eigen(Jm, symmetric = TRUE, only.values = TRUE)$values)
  for (it in 1:3) {
    lg <- .firk_legendre(s, x)
    x <- x - lg$p / lg$dp
  }
  lg <- .firk_legendre(s, x)
  list(x = x, w = 2 / ((1 - x * x) * lg$dp * lg$dp))
}

# Left Gauss-Radau nodes on [-1, 1] (gaussradau): -1 and the zeros of the
# Jacobi polynomial P^(0, 1)_(s - 1), polished by Newton's method on
# P_s + P_(s - 1).
.firk_radau_nodes <- function(s) {
  if (s == 1L) return(-1)
  m <- s - 1L
  k <- seq_len(m) - 1L
  a <- 1 / ((2 * k + 1) * (2 * k + 3))
  Jm <- diag(a, m)
  if (m > 1L) {
    kk <- seq_len(m - 1L)
    b <- sqrt(kk * (kk + 1)) / (2 * kk + 1)
    Jm[cbind(kk, kk + 1L)] <- Jm[cbind(kk + 1L, kk)] <- b
  }
  x <- sort(eigen(Jm, symmetric = TRUE, only.values = TRUE)$values)
  for (it in 1:3) {
    lg <- .firk_legendre(s, x)
    dprev <- if (s > 1L) .firk_legendre(s - 1L, x)$dp else 0
    x <- x - (lg$p + lg$p_prev) / (lg$dp + dprev)
  }
  c(-1, x)
}

# generateRadauTableau: the Radau IIA collocation matrix from the nodes,
# its eigen decomposition (eigenvalues sorted by real then imaginary part,
# as LinearAlgebra.eigen) and the error weights. Julia caches the 3, 5 and 7
# stage tableaus at load time; those come from the dump in .FIRK_TAB and
# this builds the others (their coefficients agree with Julia's to rounding).
.firk_radau_tab <- function(s) {
  key <- paste0("Radau", s)
  if (!is.null(.FIRK_TAB[[key]])) return(.FIRK_TAB[[key]])
  cc <- rev((1 - .firk_radau_nodes(s)) / 2)
  cp <- outer(cc, seq_len(s) - 1L, `^`)
  cq <- outer(seq_len(s), seq_len(s), function(i, j) cp[cbind(i, j)] * cc[i] / j)
  a <- t(solve(t(cp), t(cq)))
  ev <- eigen(a)
  ord <- order(Re(ev$values), Im(ev$values))
  val <- 1 / ev$values[ord]
  vec <- ev$vectors[, ord, drop = FALSE]
  odd <- seq(1L, s - 1L, by = 2L)
  Tm <- matrix(0, s, s)
  for (i in seq(2L, s - 1L, by = 2L)) {
    v <- vec[, i] / vec[s, i]
    Tm[, i] <- Re(v)
    Tm[, i + 1L] <- Im(v)
  }
  Tm[, 1L] <- Re(vec[, s])
  A <- t(cp) / seq_len(s)
  list(T = Tm, TI = solve(Tm), c = cc, gamma = Re(val[s]), alpha = Re(val[odd]), beta = Im(val[odd]),
       e = solve(A, c(-(s^2), -0.5, numeric(s - 2L))))
}

# generateGaussLegendreTableau: nodes, weights, the collocation matrix by a
# 2s point quadrature of the Lagrange basis on [0, c_i], and the embedded
# weights of order s - 1 (e = b - bhat). The 2 to 5 stage tableaus come from
# the dump of Julia's load time cache.
.firk_gauss_tab <- function(s) {
  key <- paste0("Gauss", s)
  if (!is.null(.FIRK_TAB[[key]])) return(.FIRK_TAB[[key]])
  gr <- .firk_gauss_rule(s)
  cc <- (gr$x + 1) / 2
  b <- gr$w / 2
  q <- .firk_gauss_rule(2L * s)
  A <- matrix(0, s, s)
  for (i in seq_len(s)) for (j in seq_len(s)) {
    nodes <- (q$x + 1) / 2 * cc[i]
    wts <- q$w / 2 * cc[i]
    L <- rep(1, length(nodes))
    for (k in seq_len(s)) if (k != j) L <- L * ((nodes - cc[k]) / (cc[j] - cc[k]))
    A[i, j] <- sum(wts * L)
  }
  e <- if (s > 1L) {
    V <- outer(seq_len(s), cc, function(i, cj) cj^(i - 1L))
    b - solve(V, c(1 / seq_len(s - 1L), 0))
  } else b
  list(A = A, b = b, c = cc, e = e)
}

# ---------------------------------------------------------------------------
# Shared pieces

# The step matrices -scale M + J (UniformScaling mass: the diagonal), real
# or complex.
.firk_shift <- function(J, M, scale) {
  if (is.null(M)) {
    A <- if (is.complex(scale)) matrix(complex(real = J, imaginary = 0), nrow(J)) else J
    diag(A) <- diag(A) - scale
    A
  } else J - scale * M
}

# Collocation interpolation through (0, 0) and (c_i, z_i): the Lagrange
# basis of the fixed stage methods (l_i in node order, products left to
# right) summed under @muladd, and the loop of AdaptiveRadau.
.firk_radau_interp <- function(theta, dt, y0, y1, k) {
  cc <- k$c; z <- k$z; s <- length(cc)
  if (isTRUE(k$adaptive)) {
    for (i in seq_len(s)) {
      tmp <- z[[i]]
      for (j in seq_len(s)) if (j != i) tmp <- tmp * ((theta - cc[j]) / (cc[i] - cc[j]))
      y0 <- y0 + tmp * theta / cc[i]
    }
    return(y0)
  }
  out <- y0
  for (i in seq_len(s)) {
    num <- theta - 0; den <- cc[i] - 0
    for (j in seq_len(s)) if (j != i) {
      num <- num * (theta - cc[j])
      den <- den * (cc[i] - cc[j])
    }
    out <- .fma(z[[i]], num / den, out)
  }
  out
}

# After an event OrdinaryDiffEq evaluates the FSAL derivative again; the
# Newton state stays and the next step sees a derivative discontinuity.
.firk_on_event <- function(ctx, t, u, cache) {
  cache$fsal <- as.numeric(ctx$f(u, ctx$p, t))
  cache$ulast <- NULL
  list(cache = cache, nf = 1L)
}

# The step factor of the predictive controller for the FIRK methods
# (fac_default_gamma false, isfirk): gamma scaled down by the Newton
# iterations of the step.
.firk_pred_fac <- function(maxiters) function(cache, gamma)
  min(gamma, (1 + 2 * maxiters) * gamma / (cache$iter + 2 * maxiters))

# ---------------------------------------------------------------------------
# The Radau steps (RadauIIA3, RadauIIA5, RadauIIA9 and AdaptiveRadau
# ConstantCache perform_step!)

.firk_radau_init <- function(ctx, t, u, max_stages, adaptive_tabs = NULL, min_stages = NULL) {
  cache <- new.env(parent = emptyenv())
  cache$fsal <- as.numeric(ctx$f(u, ctx$p, t))
  cache$etaold <- 1
  cache$iter <- 10000L
  cache$dtprev <- 1
  cache$status <- "Convergence"
  cache$kz <- rep(list(numeric(length(u))), max_stages)
  if (!is.null(adaptive_tabs)) {
    cache$tabs <- adaptive_tabs
    cache$num_stages <- min_stages
    cache$index <- 1L
    cache$step <- 1L
    cache$hist_iter <- 0
  }
  cache$nf0 <- 1L
  cache
}

# The extrapolated starting values from the collocation polynomial of the
# last accepted step, in each method's own expression shapes.
.firk_extrapolate <- function(kind, zs, cc, ratio) {
  s <- length(cc)
  if (kind == "IIA3") {
    c1 <- cc[1L]; c1m1 <- c1 - 1
    c2p <- ratio; c1p <- c1 * c2p
    k3 <- (zs[[1L]] - zs[[2L]]) / c1m1
    k4 <- k3 - (zs[[1L]] / c1)
    return(list(c1p * .fma(c1p - c1m1, k4, k3), c2p * .fma(c2p - c1m1, k4, k3)))
  }
  if (kind == "IIA5") {
    c1 <- cc[1L]; c2 <- cc[2L]
    c1m1 <- c1 - 1; c2m1 <- c2 - 1; c1mc2 <- c1 - c2
    c3p <- ratio; cp <- c(c1 * c3p, c2 * c3p, c3p)
    k1 <- (zs[[2L]] - zs[[3L]]) / c2m1
    tmp <- (zs[[1L]] - zs[[2L]]) / c1mc2
    k2 <- (tmp - k1) / c1m1
    k3 <- k2 - (tmp - zs[[1L]] / c1) / c2
    return(lapply(cp, function(x) x * .fma(x - c2m1, .fma(x - c1m1, k3, k2), k1)))
  }
  if (kind == "IIA9") {
    c1 <- cc[1L]; c2 <- cc[2L]; c3 <- cc[3L]; c4 <- cc[4L]
    c1m1 <- c1 - 1; c2m1 <- c2 - 1; c3m1 <- c3 - 1; c4m1 <- c4 - 1
    c1mc2 <- c1 - c2; c1mc3 <- c1 - c3; c1mc4 <- c1 - c4
    c2mc3 <- c2 - c3; c2mc4 <- c2 - c4; c3mc4 <- c3 - c4
    c5p <- ratio; cp <- c(c1 * c5p, c2 * c5p, c3 * c5p, c4 * c5p, c5p)
    z1 <- zs[[1L]]; z2 <- zs[[2L]]; z3 <- zs[[3L]]; z4 <- zs[[4L]]; z5 <- zs[[5L]]
    k3 <- (z4 - z5) / c4m1
    tmp1 <- (z3 - z4) / c3mc4
    k4 <- (tmp1 - k3) / c3m1
    tmp2 <- (z2 - z3) / c2mc3
    tmp3 <- (tmp2 - tmp1) / c2mc4
    k5 <- (tmp3 - k4) / c2m1
    tmp4 <- (z1 - z2) / c1mc2
    tmp5 <- (tmp4 - tmp2) / c1mc3
    tmp6 <- (tmp5 - tmp3) / c1mc4
    k6 <- (tmp6 - k5) / c1m1
    tmp7 <- z1 / c1
    tmp8 <- (tmp4 - tmp7) / c2
    tmp9 <- (tmp5 - tmp8) / c3
    tmp10 <- (tmp6 - tmp9) / c4
    k7 <- k6 - tmp10
    return(lapply(cp, function(x)
      x * .fma(x - c4m1, .fma(x - c3m1, .fma(x - c2m1, .fma(x - c1m1, k7, k6), k5), k4), k3)))
  }
  # AdaptiveRadau: the divided differences of the stage values
  D <- matrix(list(), s, s)
  D[[1L, 1L]] <- zs[[1L]] / cc[1L]
  for (j in seq_len(s)[-1L]) D[[1L, j]] <- (zs[[j - 1L]] - zs[[j]]) / (cc[j - 1L] - cc[j])
  for (i in seq_len(s)[-1L]) {
    D[[i, i]] <- (D[[i - 1L, i]] - D[[i - 1L, i - 1L]]) / cc[i]
    for (j in seq_len(s)[seq_len(s) > i]) D[[i, j]] <- (D[[i - 1L, j - 1L]] - D[[i - 1L, j]]) / (cc[j - i] - cc[j])
  }
  cp <- c(cc[-s] * ratio, ratio)
  lapply(seq_len(s), function(i) {
    zi <- D[[s, s]]
    for (j in rev(seq_len(s - 1L))) zi <- .fma(zi, cp[i] - cc[s - j] + 1, D[[j, s]])
    zi * cp[i]
  })
}

# z = T w: fused sums, with the last row of each method written out
# (T32 = 1, T33 = 0 of RadauIIA5; T52 = T54 = 1 of RadauIIA9; the unit
# entries of AdaptiveRadau's last row added without fusing).
.firk_t_times <- function(kind, Tm, w) {
  s <- length(w)
  z <- lapply(seq_len(s - 1L), function(i) .firk_chain(Tm, i, w))
  z[[s]] <- switch(kind,
    IIA3 = .firk_chain(Tm, s, w),
    IIA5 = .fma(Tm[3L, 1L], w[[1L]], w[[2L]]),
    IIA9 = .fma(Tm[5L, 1L], w[[1L]], w[[2L]] + w[[4L]]),
    {
      acc <- Tm[s, 1L] * w[[1L]]
      for (i in seq(2L, s - 1L, by = 2L)) acc <- acc + w[[i]]
      acc
    })
  z
}

.firk_radau_step <- function(ctx, t, u, dt, cache, alg) {
  f <- ctx$f; p <- ctx$p; M <- ctx$mass
  kind <- alg$firk_kind
  adaptive_order <- kind == "adaptive"
  tab <- if (adaptive_order) cache$tabs[[cache$index]] else alg$tab
  cc <- tab$c; s <- length(cc); n <- length(u)
  maxiters <- alg$maxiters
  nf <- 0L; nsolve <- 0L
  # the transformed tolerances of the Newton norm and the error estimate
  # (radau5.f)
  ex <- switch(kind, IIA3 = 3 / 4, IIA5 = 2 / 3, IIA9 = 3 / 5, (s + 1) / (s * 2))
  rtol <- ctx$reltol^ex / 10
  atol <- rtol * (ctx$abstol / ctx$reltol)
  kappa <- if (kind == "IIA3") 1 / 100 else if (!is.null(alg$kappa)) alg$kappa else if (adaptive_order) {
    # newton_tolerance
    r <- min(rtol)
    max(10 * .Machine$double.eps / r, min(3 / 100, r^((s - 1) / (s + 1))))
  } else 1 / 100
  jj <- .calc_J(ctx, u, t)
  nf <- nf + jj$nf
  J <- jj$J
  # the real block (odd s) and the complex pairs (stages 2k, 2k + 1; stages
  # 1 and 2 of the two stage RadauIIA3)
  real_block <- !is.null(tab$gamma)
  pairs <- if (real_block) lapply(seq_along(tab$alpha), function(k) c(2L * k, 2L * k + 1L)) else list(c(1L, 2L))
  gdt <- if (real_block) tab$gamma / dt
  adt <- tab$alpha / dt; bdt <- tab$beta / dt
  LU1 <- if (real_block) .firk_shift(J, M, gdt)
  LU2 <- lapply(seq_along(pairs), function(k) .firk_shift(J, M, complex(real = adt[k], imaginary = bdt[k])))
  first <- ctx$iter == 1L || isTRUE(ctx$derivative_discontinuity) || alg$extrapolant == "constant"
  if (first) {
    cache$dtprev <- 1
    z <- rep(list(numeric(n)), s)
    w <- z
    cache$kz[seq_len(s)] <- z
  } else {
    z <- .firk_extrapolate(kind, cache$kz[seq_len(s)], cc, dt / cache$dtprev)
    w <- lapply(seq_len(s), function(i) .firk_chain(tab$TI, i, z))
  }
  # integrator.u: the state of the last attempt (a rejected candidate)
  ucur <- if (isTRUE(ctx$accept_prev)) u else cache$ulast %||% u
  ndw <- 1
  eta <- max(cache$etaold, .Machine$double.eps)^0.8
  fail <- TRUE
  iter <- 0L
  ff <- rep(list(numeric(n)), s)
  while (iter < maxiters) {
    iter <- iter + 1L
    ff <- lapply(seq_len(s), function(i) as.numeric(f(u + z[[i]], p, .fma(cc[i], dt, t))))
    nf <- nf + s
    fw <- lapply(seq_len(s), function(i) .firk_chain(tab$TI, i, ff))
    Mw <- if (is.null(M)) w else lapply(w, function(x) as.numeric(M %*% x))
    dw <- vector("list", s)
    if (real_block) {
      dw[[1L]] <- as.numeric(solve(LU1, .fma(-gdt, Mw[[1L]], fw[[1L]]), tol = 0))
      nsolve <- nsolve + 1L
    }
    for (k in seq_along(pairs)) {
      a <- pairs[[k]][1L]; b <- pairs[[k]][2L]
      rhs_a <- .fma(bdt[k], Mw[[b]], .fma(-adt[k], Mw[[a]], fw[[a]]))
      rhs_b <- .fma(-adt[k], Mw[[b]], .fma(-bdt[k], Mw[[a]], fw[[b]]))
      x <- solve(LU2[[k]], complex(real = rhs_a, imaginary = rhs_b), tol = 0)
      dw[[a]] <- Re(x); dw[[b]] <- Im(x)
      nsolve <- nsolve + 1L
    }
    ndwprev <- ndw
    ndw <- 0
    for (i in seq_len(s)) ndw <- ndw + .firk_err(dw[[i]], u, ucur, atol, rtol)
    if (iter > 1L) {
      theta <- ndw / ndwprev
      diverge <- isTRUE(theta > 1)
      if (diverge) cache$status <- "Divergence"
      slow <- isTRUE(ndw * theta^(maxiters - iter) > kappa * (1 - theta))
      if (slow) cache$status <- "VerySlowConvergence"
      if (diverge || slow) break
      if (!(adaptive_order && ndw == 0)) eta <- theta / (1 - theta)
    }
    w <- lapply(seq_len(s), function(i) w[[i]] - dw[[i]])
    z <- .firk_t_times(kind, tab$T, w)
    if (isTRUE(eta * ndw < kappa) && (iter > 1L || ndw == 0 || ctx$success_iter != 0L)) {
      cache$status <- if (eta < alg$fast_convergence_cutoff) "FastConvergence" else "Convergence"
      fail <- FALSE
      break
    }
  }
  out <- list(nf = nf, njacs = 1L, nw = 1L, nsolve = nsolve, cache = cache)
  # the out of place RadauIIA3 of Julia has no convergence check
  # (firk_iia3_check = FALSE reproduces it); the in place one fails the step
  if (fail && (kind != "IIA3" || isTRUE(alg$firk_iia3_check)))
    return(c(out, list(fail = TRUE, nreject = 0L, u = u)))
  cache$etaold <- eta
  cache$iter <- iter
  unew <- u + z[[s]]
  EEst <- 1
  if (ctx$adaptive) {
    if (kind == "IIA3") {
      utilde <- dt * .fma(tab$e[2L], ff[[2L]], tab$e[1L] * ff[[1L]])
      EEst <- .firk_err(utilde, u, unew, atol, rtol)
    } else {
      edt <- tab$e / dt
      tmp <- edt[1L] * z[[1L]]
      for (i in seq_len(s)[-1L]) tmp <- .fma(edt[i], z[[i]], tmp)
      if (!is.null(M)) tmp <- as.numeric(M %*% tmp)
      utilde <- cache$fsal + tmp
      if (alg$smooth_est) {
        utilde <- as.numeric(solve(LU1, utilde, tol = 0))
        nsolve <- nsolve + 1L
      }
      EEst <- .firk_err(utilde, u, unew, atol, rtol)
      if ((!isTRUE(EEst < 1) && ctx$iter == 1L) || isTRUE(ctx$derivative_discontinuity)) {
        f0 <- as.numeric(f(u + utilde, p, t))
        nf <- nf + 1L
        utilde <- f0 + tmp
        if (alg$smooth_est) {
          utilde <- as.numeric(solve(LU1, utilde, tol = 0))
          nsolve <- nsolve + 1L
        }
        EEst <- .firk_err(utilde, u, unew, atol, rtol)
      }
    }
  }
  accepted <- !ctx$adaptive || isTRUE(EEst <= 1)
  if (accepted) {
    cache$dtprev <- dt
    if (alg$extrapolant != "constant") cache$kz[seq_len(s)] <- z
  }
  fnew <- as.numeric(f(unew, p, t + dt))
  # the out of place RadauIIA3 does not count this evaluation
  if (kind != "IIA3") nf <- nf + 1L
  k <- list(cache$fsal, fnew, z = z, c = cc, adaptive = adaptive_order)
  if (accepted) cache$fsal <- fnew
  cache$ulast <- unew
  out$nf <- nf; out$nsolve <- nsolve
  c(out, list(u = unew, EEst = EEst, k = k))
}

# AdaptiveRadau's step_accept_controller! and step_reject_controller!
# bookkeeping: the average Newton iteration count moves the method up or
# down by two stages after more than ten steps at the current count.
.firk_ar_update <- function(cache, min_stages, max_stages, accept) {
  step <- cache$step; ns <- cache$num_stages; it <- cache$iter
  cache$step <- step + 1L
  hist <- cache$hist_iter * 0.8 + it * 0.2
  cache$hist_iter <- hist
  if (step > 10L) {
    slow <- hist > 8 || cache$status %in% c("VerySlowConvergence", "Divergence")
    if (accept && hist < 2.6 && ns < max_stages) {
      cache$num_stages <- ns + 2L
      cache$index <- cache$index + 1L
      cache$step <- 1L
      cache$hist_iter <- it
    } else if (slow && ns > min_stages) {
      cache$num_stages <- ns - 2L
      cache$index <- cache$index - 1L
      cache$step <- 1L
      cache$hist_iter <- it
    }
  }
  invisible(cache)
}

# ---------------------------------------------------------------------------
# GaussLegendre (GaussLegendreConstantCache perform_step! and
# _gausslegendre_substep_constant)

.firk_gauss_init <- function(ctx, t, u, s) {
  cache <- new.env(parent = emptyenv())
  cache$fsal <- as.numeric(ctx$f(u, ctx$p, t))
  cache$etaold <- 1
  cache$iter <- 10000L
  cache$dtprev <- 1
  cache$status <- "Convergence"
  cache$kz <- rep(list(numeric(length(u))), s)
  cache$nf0 <- 1L
  cache
}

.firk_gauss_step <- function(ctx, t, u, dt, cache, alg) {
  f <- ctx$f; p <- ctx$p
  tab <- alg$tab
  A <- tab$A; b <- tab$b; cc <- tab$c; e <- tab$e
  s <- length(cc); n <- length(u)
  maxiters <- alg$maxiters
  kappa <- alg$kappa %||% (1 / 100)
  rtol <- ctx$reltol^((s + 1) / (s * 2)) / 10
  atol <- rtol * (ctx$abstol / ctx$reltol)
  jj <- .calc_J(ctx, u, t)
  nf <- jj$nf
  if (ctx$iter == 1L || isTRUE(ctx$derivative_discontinuity) || alg$extrapolant == "constant") {
    cache$dtprev <- 1
    z <- rep(list(numeric(n)), s)
    cache$kz <- z
  } else {
    # the collocation polynomial through (0, 0) and (c_i, z_i)
    zs <- cache$kz
    ratio <- dt / cache$dtprev
    cp <- cc * ratio
    D <- matrix(list(), s, s)
    D[[1L, 1L]] <- zs[[1L]] / cc[1L]
    for (j in seq_len(s)[-1L]) D[[1L, j]] <- (zs[[j - 1L]] - zs[[j]]) / (cc[j - 1L] - cc[j])
    for (i in seq_len(s)[-1L]) {
      D[[i, i]] <- (D[[i - 1L, i]] - D[[i - 1L, i - 1L]]) / cc[i]
      for (j in seq_len(s)[seq_len(s) > i]) D[[i, j]] <- (D[[i - 1L, j - 1L]] - D[[i - 1L, j]]) / (cc[j - i] - cc[j])
    }
    z <- lapply(seq_len(s), function(i) {
      zi <- D[[s, s]]
      for (j in rev(seq_len(s - 1L))) zi <- .fma(zi, cp[i] - cc[j], D[[j, j]])
      zi * cp[i]
    })
  }
  # W = I - dt (A kron J), one LU per step
  W <- diag(s * n) - dt * kronecker(A, jj$J)
  ndw <- 1
  eta <- max(cache$etaold, .Machine$double.eps)^0.8
  ok <- FALSE
  iter <- 0L
  nsolve <- 0L
  while (iter < maxiters) {
    iter <- iter + 1L
    ff <- lapply(seq_len(s), function(i) as.numeric(f(u + z[[i]], p, .fma(cc[i], dt, t))))
    nf <- nf + s
    r <- unlist(lapply(seq_len(s), function(i) .fma(-dt, .firk_chain(A, i, ff), z[[i]])))
    dwf <- as.numeric(solve(W, r, tol = 0))
    nsolve <- nsolve + 1L
    dw <- lapply(seq_len(s), function(i) dwf[(i - 1L) * n + seq_len(n)])
    ndwprev <- ndw
    ndw <- .firk_err(dw[[1L]], u, u, atol, rtol)
    for (i in seq_len(s)[-1L]) ndw <- ndw + .firk_err(dw[[i]], u, u, atol, rtol)
    if (iter > 1L) {
      theta <- ndw / ndwprev
      diverge <- isTRUE(theta > 1)
      if (diverge) cache$status <- "Divergence"
      slow <- isTRUE(ndw * theta^(maxiters - iter) > kappa * (1 - theta))
      if (slow) cache$status <- "VerySlowConvergence"
      if (diverge || slow) break
      eta <- theta / (1 - theta)
    }
    z <- lapply(seq_len(s), function(i) z[[i]] - dw[[i]])
    if (isTRUE(eta * ndw < kappa) && (iter > 1L || ndw == 0 || ctx$success_iter != 0L)) {
      cache$status <- if (eta < alg$fast_convergence_cutoff) "FastConvergence" else "Convergence"
      ok <- TRUE
      break
    }
  }
  cache$etaold <- eta
  cache$iter <- iter
  out <- list(njacs = 1L, nw = 1L, cache = cache)
  if (!ok) return(c(out, list(nf = nf, nsolve = nsolve, fail = TRUE, nreject = 0L, u = u)))
  ff <- lapply(seq_len(s), function(i) as.numeric(f(u + z[[i]], p, .fma(cc[i], dt, t))))
  nf <- nf + s
  unew <- u
  for (i in seq_len(s)) unew <- .fma(dt * b[i], ff[[i]], unew)
  EEst <- 1
  if (ctx$adaptive) {
    utilde <- (dt * e[1L]) * ff[[1L]]
    for (i in seq_len(s)[-1L]) utilde <- .fma(dt * e[i], ff[[i]], utilde)
    EEst <- .firk_err(utilde, u, unew, atol, rtol)
  }
  accepted <- !ctx$adaptive || isTRUE(EEst <= 1)
  if (accepted) {
    cache$dtprev <- dt
    if (alg$extrapolant != "constant") cache$kz <- z
  }
  fnew <- as.numeric(f(unew, p, t + dt))
  nf <- nf + 1L
  k <- list(cache$fsal, fnew)
  if (accepted) cache$fsal <- fnew
  c(out, list(nf = nf, nsolve = nsolve, u = unew, EEst = EEst, k = k))
}

# ---------------------------------------------------------------------------
# Constructors

.firk_check_extrapolant <- function(extrapolant) {
  extrapolant <- as.character(extrapolant)
  if (!extrapolant %in% c("dense", "constant"))
    stop("`extrapolant` must be \"dense\" or \"constant\"", call. = FALSE)
  extrapolant
}

.firk_radau_alg <- function(name, kind, order, adaptive_order, autodiff, extrapolant, fast_convergence_cutoff,
                            new_W_gamma_dt_cutoff, kappa, maxiters, smooth_est, step_limiter) {
  tab <- .firk_iia_tab(name)
  alg <- NULL
  alg <- .with_autodiff(autodiff, .new_alg(
    name, order, step = function(ctx, t, u, dt, cache) .firk_radau_step(ctx, t, u, dt, cache, alg),
    interp = .firk_radau_interp,
    init = function(ctx, t, u) .firk_radau_init(ctx, t, u, length(tab$c)),
    stiff = TRUE, controller = "Predictive", adaptive_order = adaptive_order, qmax = 8,
    qsteady_max = 6 / 5, pred_fac = .firk_pred_fac(maxiters), on_event = .firk_on_event,
    firk_kind = kind, tab = tab, extrapolant = .firk_check_extrapolant(extrapolant),
    fast_convergence_cutoff = fast_convergence_cutoff, new_W_gamma_dt_cutoff = new_W_gamma_dt_cutoff,
    kappa = kappa, maxiters = as.integer(maxiters), smooth_est = isTRUE(smooth_est), firk_iia3_check = TRUE,
    step_limiter = step_limiter, mass_matrix = TRUE))
  alg
}

#' Radau IIA methods
#'
#' Fully implicit Radau IIA collocation methods of orders 3, 5 and 9 for
#' stiff problems, after OrdinaryDiffEqFIRK. The stage system is solved in
#' the eigenbasis of the inverse Runge-Kutta matrix: one real and one or two
#' complex linear systems per Newton iteration, factorized once per step with
#' a fresh Jacobian. The simplified Newton iteration starts from the
#' collocation polynomial of the last step and stops on its convergence rate
#' estimate; the error estimate is Hairer's, smoothed through the real
#' block (`smooth_est`). The step size follows the predictive (Gustafsson)
#' controller with the factor reduced by the Newton iteration count.
#' `RadauIIA5()` is the method of Hairer and Wanner's RADAU5.
#'
#' Two corrections to Julia: the out of place `RadauIIA3` of Julia accepts a
#' step whose Newton iteration diverged; here a failed iteration fails the
#' step and retries with half the step, as Julia's in place method does. The
#' dense output always uses the collocation polynomial of its own step (with
#' `extrapolant = "constant"` Julia interpolates with zero stage values).
#'
#' @param autodiff Differentiation choice for the Jacobian when the problem
#'   has none (see [AutoForwardDiff()]).
#' @param extrapolant `"dense"` (the default: the previous collocation
#'   polynomial gives the starting values) or `"constant"` (zero starting
#'   values).
#' @param fast_convergence_cutoff Rate below which the Newton iteration counts
#'   as fast.
#' @param new_W_gamma_dt_cutoff Stored only; the out of place methods
#'   factorize every step.
#' @param kappa Tolerance of the Newton stopping test (default 1/100;
#'   `RadauIIA3()` always uses 1/100, as Julia).
#' @param maxiters Maximum number of Newton iterations.
#' @param smooth_est Whether the error estimate is smoothed through the real
#'   block.
#' @param step_limiter Optional `function(u, integrator, p, t)` applied after
#'   every accepted step.
#' @param concrete_jac,linsolve Accepted for compatibility; dense LU
#'   factorizations are used.
#' @return An ODE algorithm.
#' @examples
#' vdp <- function(u, p, t) c(u[2], p[1] * ((1 - u[1] * u[1]) * u[2]) - u[1])
#' prob <- ODEProblem(vdp, c(2, 0), c(0, 1), 10)
#' sol <- solve(prob, RadauIIA5(), reltol = 1e-6, abstol = 1e-8)
#' sol$u[[length(sol$u)]]
#' @name radau_methods
NULL

#' @rdname radau_methods
#' @export
RadauIIA3 <- function(autodiff = NULL, extrapolant = "dense", fast_convergence_cutoff = 1 / 5,
                      new_W_gamma_dt_cutoff = 1 / 5, kappa = NULL, maxiters = 10L, step_limiter = NULL,
                      concrete_jac = NULL, linsolve = NULL) {
  .firk_radau_alg("RadauIIA3", "IIA3", 3L, 1L, autodiff, extrapolant, fast_convergence_cutoff,
                  new_W_gamma_dt_cutoff, kappa, maxiters, FALSE, step_limiter)
}

#' @rdname radau_methods
#' @export
RadauIIA5 <- function(autodiff = NULL, extrapolant = "dense", fast_convergence_cutoff = 1 / 5,
                      new_W_gamma_dt_cutoff = 1 / 5, kappa = NULL, maxiters = 10L, smooth_est = TRUE,
                      step_limiter = NULL, concrete_jac = NULL, linsolve = NULL) {
  .firk_radau_alg("RadauIIA5", "IIA5", 5L, 3L, autodiff, extrapolant, fast_convergence_cutoff,
                  new_W_gamma_dt_cutoff, kappa, maxiters, smooth_est, step_limiter)
}

#' @rdname radau_methods
#' @export
RadauIIA9 <- function(autodiff = NULL, extrapolant = "dense", fast_convergence_cutoff = 1 / 5,
                      new_W_gamma_dt_cutoff = 1 / 5, kappa = NULL, maxiters = 10L, smooth_est = TRUE,
                      step_limiter = NULL, concrete_jac = NULL, linsolve = NULL) {
  .firk_radau_alg("RadauIIA9", "IIA9", 9L, 5L, autodiff, extrapolant, fast_convergence_cutoff,
                  new_W_gamma_dt_cutoff, kappa, maxiters, smooth_est, step_limiter)
}

#' Adaptive order Radau IIA method
#'
#' The Radau IIA collocation methods with 3, 5, 7, ... stages (orders 5, 9,
#' 13, ...) under one controller, after OrdinaryDiffEqFIRK's
#' `AdaptiveRadau`: the method starts at the stage count of `min_order` and
#' moves up two stages when the Newton iteration has converged in fewer than
#' 2.6 iterations on average over more than ten steps, down when it needed
#' more than 8. The tableaus of 3, 5 and 7 stages are Julia's; larger ones
#' are computed (Gauss-Radau nodes and the eigen decomposition of the
#' collocation matrix) and agree with Julia's to rounding. The Newton
#' tolerance follows the transformed relative tolerance unless `kappa` is
#' given. The dense output uses the stage count of each step (Julia
#' interpolates with the stage count in force at the end).
#'
#' A correction to Julia: with a mass matrix, Julia's out of place step
#' computes `mass_matrix * w` with `w` the vector of stage vectors, which
#' mixes the stages (and fails unless the state has as many components as
#' the method has stages); here the mass matrix multiplies each stage.
#'
#' @inheritParams radau_methods
#' @param threading Accepted for compatibility; the stage systems are solved
#'   serially.
#' @param min_order,max_order Order range (at least 5); the stage counts are
#'   `(order - 1) %/% 4 * 2 + 1`.
#' @return An ODE algorithm.
#' @examples
#' vdp <- function(u, p, t) c(u[2], p[1] * ((1 - u[1] * u[1]) * u[2]) - u[1])
#' prob <- ODEProblem(vdp, c(2, 0), c(0, 1), 10)
#' sol <- solve(prob, AdaptiveRadau(), reltol = 1e-8, abstol = 1e-10)
#' sol$u[[length(sol$u)]]
#' @export
AdaptiveRadau <- function(autodiff = NULL, min_order = 5L, max_order = 13L, extrapolant = "dense",
                          fast_convergence_cutoff = 1 / 5, new_W_gamma_dt_cutoff = 1 / 5, kappa = NULL,
                          maxiters = 10L, smooth_est = TRUE, step_limiter = NULL, concrete_jac = NULL,
                          linsolve = NULL, threading = FALSE) {
  max_stages <- (as.integer(max_order) - 1L) %/% 4L * 2L + 1L
  min_stages <- (as.integer(min_order) - 1L) %/% 4L * 2L + 1L
  if (min_order < 5) stop("min_order choice ", min_order, " below 5 is not compatible with the algorithm",
                          call. = FALSE)
  if (max_stages < min_stages) stop("max_order ", max_order, " is below min_order ", min_order, call. = FALSE)
  tabs <- lapply(seq(min_stages, max_stages, by = 2L), .firk_radau_tab)
  alg <- NULL
  alg <- .with_autodiff(autodiff, .new_alg(
    "AdaptiveRadau", 2L * min_stages - 1L,
    step = function(ctx, t, u, dt, cache) .firk_radau_step(ctx, t, u, dt, cache, alg),
    interp = .firk_radau_interp,
    init = function(ctx, t, u) .firk_radau_init(ctx, t, u, max_stages, tabs, min_stages),
    stiff = TRUE, controller = "Predictive", adaptive_order = min_stages, qmax = 8, qsteady_max = 6 / 5,
    pred_fac = .firk_pred_fac(maxiters), cur_adaptive_order = function(cache) cache$num_stages,
    on_accept = function(cache) .firk_ar_update(cache, min_stages, max_stages, TRUE),
    on_reject = function(cache) .firk_ar_update(cache, min_stages, max_stages, FALSE),
    on_event = .firk_on_event, firk_kind = "adaptive", extrapolant = .firk_check_extrapolant(extrapolant),
    fast_convergence_cutoff = fast_convergence_cutoff, new_W_gamma_dt_cutoff = new_W_gamma_dt_cutoff,
    kappa = kappa, maxiters = as.integer(maxiters), smooth_est = isTRUE(smooth_est),
    min_order = as.integer(min_order), max_order = as.integer(max_order), step_limiter = step_limiter,
    mass_matrix = TRUE))
  alg
}

#' Gauss-Legendre collocation method
#'
#' The `num_stages` stage Gauss-Legendre implicit Runge-Kutta method of
#' order `2 num_stages` (symplectic and A-stable), after OrdinaryDiffEqFIRK's
#' `GaussLegendre`: the full stage system is solved by a simplified Newton
#' iteration with one LU factorization per step, the error estimate comes
#' from the embedded weights of order `num_stages - 1` (so adaptive stepping
#' needs at least two stages) and the step size from the PI controller. The
#' tableaus of 2 to 5 stages are Julia's; others are computed and agree with
#' Julia's to rounding. The dense output is cubic Hermite.
#'
#' @inheritParams radau_methods
#' @param num_stages Number of stages (a positive integer).
#' @param precs Accepted for compatibility (no preconditioner is used).
#' @param controller Accepted for compatibility (the PI controller is used,
#'   as in Julia).
#' @return An ODE algorithm.
#' @examples
#' osc <- function(u, p, t) c(u[2], -u[1])
#' prob <- ODEProblem(osc, c(1, 0), c(0, 1))
#' sol <- solve(prob, GaussLegendre(num_stages = 3), reltol = 1e-8, abstol = 1e-10)
#' sol$u[[length(sol$u)]]
#' @export
GaussLegendre <- function(num_stages = 2L, autodiff = NULL, extrapolant = "dense", fast_convergence_cutoff = 1 / 5,
                          new_W_gamma_dt_cutoff = 1 / 5, controller = "PI", kappa = NULL, maxiters = 10L,
                          step_limiter = NULL, concrete_jac = NULL, linsolve = NULL, precs = NULL) {
  s <- as.integer(num_stages)
  if (length(s) != 1L || is.na(s) || s < 1L) stop("`num_stages` must be a positive integer", call. = FALSE)
  tab <- .firk_gauss_tab(s)
  alg <- NULL
  alg <- .with_autodiff(autodiff, .new_alg(
    "GaussLegendre", 2L * s, adaptive = s >= 2L,
    step = function(ctx, t, u, dt, cache) .firk_gauss_step(ctx, t, u, dt, cache, alg),
    init = function(ctx, t, u) .firk_gauss_init(ctx, t, u, s),
    stiff = TRUE, adaptive_order = s - 1L, qmax = 8, qsteady_max = 6 / 5, on_event = .firk_on_event,
    tab = tab, num_stages = s, extrapolant = .firk_check_extrapolant(extrapolant),
    fast_convergence_cutoff = fast_convergence_cutoff, new_W_gamma_dt_cutoff = new_W_gamma_dt_cutoff,
    firk_controller = controller, kappa = kappa, maxiters = as.integer(maxiters), step_limiter = step_limiter))
  alg
}

# ---------------------------------------------------------------------------
# Tableaus of Julia (OrdinaryDiffEqFIRK firk_tableaus.jl): the hand written
# RadauIIA3, RadauIIA5 and RadauIIA9 tableaus, and the 3, 5 and 7 stage
# Radau and 2 to 5 stage Gauss-Legendre tableaus Julia computes and caches at
# load time, dumped with repr().
.FIRK_TAB <- list(
  RadauIIA3 = list(
    T11 = 0.10540925533894596,
    T12 = -0.29814239699997197,
    T21 = 0.9486832980505138,
    T22 = 0.0,
    TI11 = 0.0,
    TI12 = 1.0540925533894598,
    TI21 = -3.3541019662496843,
    TI22 = 0.3726779962499649,
    c1 = 0.3333333333333333,
    c2 = 1.0,
    alpha = 2.0,
    beta = -1.4142135623730951,
    e1 = 0.25,
    e2 = -0.25
  ),
  RadauIIA5 = list(
    T11 = 0x1.75b01984e65f7p-4,
    T12 = -0.1412552950209542,
    T13 = -0.030029194105147424,
    T21 = 0.241717932707107,
    T22 = 0x1.a20e91e20b551p-3,
    T23 = 0.3829421127572619,
    T31 = 0.966048182615093,
    TI11 = 4.325579890063155,
    TI12 = 0.33919925181580984,
    TI13 = 0.5417705399358749,
    TI21 = -4.178718591551905,
    TI22 = -0.32768282076106237,
    TI23 = 0x1.e810014c55221p-2,
    TI31 = -0.5028726349457868,
    TI32 = 2.571926949855605,
    TI33 = -0.5960392048282249,
    c1 = 0.15505102572168222,
    c2 = 0.6449489742783178,
    gamma = 0x1.d1a48d83e731ep+1,
    alpha = 2.6810828736277523,
    beta = 3.0504301992474105,
    e1 = -10.048809399827414,
    e2 = 1.382142733160748,
    e3 = -0.3333333333333333
  ),
  RadauIIA9 = list(
    T11 = -0x1.9a2d1fb890948p-7,
    T12 = -0.010242047817908826,
    T13 = 0.04767387729029572,
    T14 = -0x1.7820c3d05d519p-7,
    T15 = -0.01401985889287541,
    T21 = -0.0014916701518953824,
    T22 = 0.05017286451737106,
    T23 = -0x1.826214e5b71e4p-4,
    T24 = -0.007668830749180163,
    T25 = 0.024708578426518527,
    T31 = -0.07298187638808715,
    T32 = -0.23053953404341795,
    T33 = 0.1027030453801259,
    T34 = 0.01939846399882895,
    T35 = 0.08180035370375117,
    T41 = -0.3800914400035681,
    T42 = 0.37789390224886127,
    T43 = 0.46674413033249434,
    T44 = 0x1.a1623391d6c89p-2,
    T45 = 0.19968242788680252,
    T51 = -0x1.d80da0c71ef28p-1,
    TI11 = -30.0415677215444,
    TI12 = -0x1.bbaef6abc6593p+3,
    TI13 = -3.4800027747951856,
    TI14 = 1.0320087978252634,
    TI15 = -0.8043030450739899,
    TI21 = 5.344186437834912,
    TI22 = 4.593615567759161,
    TI23 = -3.0363603234594243,
    TI24 = 1.050660190231459,
    TI25 = -0.27277861186429625,
    TI31 = 3.748059807439805,
    TI32 = -3.984965736343885,
    TI33 = -1.0444156416080188,
    TI34 = 1.1840985681379486,
    TI35 = -0x1.ccb73e72d9b82p-2,
    TI41 = -33.041880213519,
    TI42 = -17.376953479063566,
    TI43 = -0.17212906325400557,
    TI44 = -0.09916977798254265,
    TI45 = 0.5312281158383066,
    TI51 = -8.611443979875292,
    TI52 = 9.699991409528808,
    TI53 = 1.9147286396968743,
    TI54 = 2.41869200608494,
    TI55 = -1.0474634879353375,
    c1 = 0.05710419611451768,
    c2 = 0.2768430136381238,
    c3 = 0.5835904323689168,
    c4 = 0.8602401356562195,
    gamma = 6.2867047517292765,
    alpha1 = 3.655694325463572,
    beta1 = 6.543736899360077,
    alpha2 = 5.70095329867179,
    beta2 = 3.2102656003085497,
    e1 = -27.78093394406464,
    e2 = 3.641478498049213,
    e3 = -1.2525477211691187,
    e4 = 0.5920031671845428,
    e5 = -0.2
  ),
  Radau3 = list(
    T = matrix(c(0.09123239487089282, 0.2417179327071075, 0x1.ee9dde0d62559p-1, -0.1412552950209543, 0.2041293522938006, 1.0, -0.030029194105147327, 0x1.8821fa2a36defp-2, -0.0), 3, 3),
    TI = matrix(c(4.325579890063157, -4.178718591551905, -0.5028726349457857, 0.3391992518158088, -0.32768282076106137, 2.571926949855605, 0.5417705399358754, 0.47662355450054994, -0.5960392048282265), 3, 3),
    c = c(0.15505102572168217, 0.6449489742783178, 1.0),
    gamma = 3.6378342527444967,
    alpha = c(2.6810828736277528),
    beta = c(3.0504301992474083),
    e = c(-10.048809399827416, 0x1.61d41b2d5457fp+0, -0.333333333333333)
  ),
  Radau5 = list(
    T = matrix(c(0.01251758622050091, 0.001491670151895388, 0.07298187638808792, 0.3800914400035683, 0.9219789736812104, -0.010242047817909406, 0.050172864517372094, -0.2305395340434175, 0x1.82f69e7e1110fp-2, 1.0, 0.04767387729029541, -0.09433181918161004, 0x1.a4abf2d1b73afp-4, 0.46674413033249373, -0.0, -0.011478515255229471, -0.007668830749180366, 0.019398463998828844, 0.4076011712801977, 1.0, -0x1.cb6719b83f683p-7, 0x1.94d34a11716ebp-6, 0.08180035370375069, 0.1996824278868023, -0.0), 5, 5),
    TI = matrix(c(30.041567721544386, 5.344186437834917, 0x1.dfc06c7c3404ap+1, -33.04188021351899, -8.611443979875666, 0x1.bbaef6abc65acp+3, 4.59361556775922, -0x1.fe135b74a3815p+1, -17.376953479063662, 9.699991409528666, 3.4800027747952047, -3.036360323459413, -1.044415641608052, -0.17212906325403302, 0x1.ea2ba7f836764p+0, -1.0320087978252603, 1.050660190231443, 1.1840985681379648, -0.09916977798252875, 2.4186920060849233, 0.8043030450739889, -0.2727786118642877, -0.4499177701567845, 0.5312281158382987, -1.0474634879353304), 5, 5),
    c = c(0.05710419611451767, 0.2768430136381238, 0.5835904323689168, 0.8602401356562195, 1.0),
    gamma = 6.28670475172926,
    alpha = c(0x1.d3edcaaa04fa7p+1, 5.700953298671795),
    beta = c(6.543736899360087, 3.2102656003085626),
    e = c(-27.780933944064632, 3.641478498049206, -1.2525477211691116, 0.592003167184537, -0x1.9999999999944p-3)
  ),
  Radau7 = list(
    T = matrix(c(-0.0021537546273123417, 0.001600025077882338, -0.004059107301950501, -0.01575048807937698, -0x1.cec1a9f471647p-4, -0.4583810431839396, -0.8813915783538118, 0.02156755135132616, -0.038131648134425954, 0.05739650893942264, -0.03821469359705627, -0.24917421246511323, 0.5315846490836865, 1.0, 0.008783567925124847, -0.021525560593974338, 0.05885052920838343, -0x1.536dd24148bcep-3, 0.2735633057987673, 0.48632283661737585, 0.0, -0.004055161452330414, 0.008415568276561277, -0.00856043106161273, -0.037371242302372834, 0.005366761379201571, 0.5265742264585006, 1.0, 0.004427232753273333, -0.0040319495702327785, -0.006923212665015953, 0.008239007298518909, 0.19321111610127575, 0.27553439498960397, 0.0, -0.001238646187954856, -6.666635339249368e-5, -0.0023521809829451716, 0.00311507115235087, 0.10171773248173102, 0.5217519452747852, 0x1.fffffffffffffp-1, -0.0027606174805457897, 0x1.a186cc6485d10p-9, 0.0004169077725283774, 0.02511660491344946, 0x1.854ded613ad12p-4, 0.12807194463554966, -0.0), 7, 7),
    TI = matrix(c(-258.13192631980695, -3.0073901694642187, -8.441963188321433, 74.67833223504702, 58.35652885182899, -299.18624802810734, -93.07650289727287, -0x1.7a25c4464b852p+7, -0x1.6081f9942df31p+3, -0.6505252740447515, 87.40858897986983, -10.068773957859182, -243.04074536860531, 23.881631056347, -49.08731481790715, 1.4877994561383148, 6.940670730369991, 4.024158737361161, -30.36638884256207, -48.77710407800511, 39.27888073080811, -4.110647469658472, 2.1303881595573904, -3.2050475255994755, -3.714806315157185, -1.020020865182681, -2.0386719057385965, 14.388915685487172, -4.053447889316159, -0x1.d0ee9f4ccd13ep+0, 1.0712809435471853, -3.4300939859826443, -0.11241750037813132, 1.6735602398603917, -3.510438399399342, 3.112755366608249, 1.1343255878949288, -0.3548507491218917, 2.6966048097658537, 1.890640831000311, -1.0873740320566434, 4.863284885565753, -1.6467749135588885, -0.41469904594322143, 0x1.78c5c74f5c2b7p-4, -0.9386927436077704, -0x1.f17bee32b51b5p-1, 0x1.cdcad9a77975fp-1, -2.246482729591072), 7, 7),
    c = c(0x1.e05201f8c2060p-6, 0.1480785996684843, 0.3369846902811543, 0.5586715187715501, 0.7692338620300545, 0.9269456713197411, 1.0),
    gamma = 8.93683278840461,
    alpha = c(4.378693561507815, 7.141055219187439, 8.511834825102486),
    beta = c(0x1.456e209c07d99p+3, 6.6230459226385, 3.2810136243247583),
    e = c(-54.37443689412865, 7.0000240042593145, -2.355661091987765, 1.1322890661063771, -0.646891326767591, 0.3875333853755371, -0.14285714285721904)
  ),
  Gauss2 = list(
    A = matrix(c(0.25, 0.5386751345948129, -0.03867513459481287, 0.25000000000000006), 2, 2),
    b = c(0.5, 0.5),
    c = c(0.21132486540518708, 0.7886751345948129),
    e = c(-0.8660254037844386, 0.8660254037844386)
  ),
  Gauss3 = list(
    A = matrix(c(0.1388888888888889, 0.3002631949808647, 0.26798833376246955, -0.03597666752493891, 0.22222222222222235, 0x1.ebf38310dda6ep-2, 0.009789444015308329, -0.02248541720308683, 0.1388888888888889), 3, 3),
    b = c(0.2777777777777778, 0.4444444444444444, 0.2777777777777778),
    c = c(0.1127016653792583, 0.5, 0.8872983346207417),
    e = c(1.1111111111111107, -2.2222222222222214, 1.1111111111111112)
  ),
  Gauss4 = list(
    A = matrix(c(0.0869637112843634, 0.18811811749986793, 0.1671919219741885, 0.17748257225452227, -0.026604180084998784, 0.16303628871563652, 0.3539530060337441, 0.3134451147418685, 0.01262746268940472, -0x1.c8cafe42f2e72p-6, 0.16303628871563633, 0.3526767575162717, -0.0035551496857956825, 0.006735500594538166, -0.014190694931141137, 0.0869637112843633), 4, 4),
    b = c(0.17392742256872692, 0.3260725774312731, 0.3260725774312731, 0.17392742256872692),
    c = c(0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262),
    e = c(-1.855135017009738, 4.698862351888767, -4.698862351888765, 1.8551350170097367)
  ),
  Gauss5 = list(
    A = matrix(c(0.05923172126404729, 0.12815100567004528, 0x1.d20715c09d603p-4, 0.12123243692686413, 0x1.deb8aa62b9f92p-4, -0.019570364359076053, 0.1196571676248416, 0.2600046516806416, 0x1.d4fbe22aa380fp-3, 0.2449081289104954, 0.011254400818642962, -0.024592114619642193, 0.1422222222222222, 0.30903655906408667, 0.2731900436258016, -0.005593793660812189, 0.01031828067068336, -0.020690316430958304, 0x1.ea1da25ae4155p-4, 0.25888469960875904, 0.0015881129678660005, -0.0027689943987696036, 0.004687154523869945, -0.009687563141950718, 0.05923172126404734), 5, 5),
    b = c(0.11846344252809454, 0.23931433524968324, 0.28444444444444444, 0.23931433524968324, 0.11846344252809454),
    c = c(0.04691007703066802, 0.23076534494715845, 0.5, 0.7692346550528415, 0x1.e7fb670299d34p-1),
    e = c(3.667944222886206, -0x1.4c6a0a00befbdp+3, 13.439999999999944, -10.38794422288621, 0x1.d57f32406c9b3p+1)
  )
)
