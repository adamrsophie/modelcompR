# The GLTR trust region subproblem solver of trlib (Felix Lenders, MIT
# license), as SciPy 1.18.1 vendors it in scipy/optimize/_trlib and drives it
# from _trlib.pyx for method = "trust-krylov". The reverse communication state
# machine of trlib_krylov.c, the tridiagonal subproblem of trlib_tri_factor.c,
# trlib_leftmost.c, trlib_eigen_inverse.c, trlib_quadratic_zero.c and the
# LAPACK routines dpttrf, dpttrs, dptrfs and dlagtm are ported. SciPy fixes
# equality = 0, ctl_invariant = 0, convexify = 1, earlyterm = 1, refine = 1 and
# itmax_lanczos = 100, so only one irreducible block ever occurs; the branches
# that need other settings are not reachable from SciPy and are not ported.

.sptl_eps <- 2.2204460492503131e-16
.sptl_eps_pow_4 <- 0x1.2611186bae670p-21
.sptl_eps_pow_5 <- 1.4901161193847656e-08
.sptl_eps_pow_75 <- 1.8189894035458565e-12

.sptl_nrm2 <- function(x) sqrt(sum(x * x))

# comparisons follow C, where any comparison with NaN is false and fmin and
# fmax ignore a NaN argument
.sptl_sign1 <- function(a) if (isTRUE(a < 0) || (isTRUE(a == 0) && 1 / a < 0)) -1 else 1

.sptl_fmax <- function(a, b) if (is.na(a)) b else if (is.na(b)) a else max(a, b)

.sptl_fmin <- function(a, b) if (is.na(a)) b else if (is.na(b)) a else min(a, b)

# LAPACK dpttrf: L D L' factorization of a symmetric positive definite
# tridiagonal matrix; returns the partial factors on failure as LAPACK does
.sptl_pttrf <- function(d, e) {
  n <- length(d)
  for (i in seq_len(n - 1L)) {
    if (isTRUE(d[i] <= 0)) return(list(d = d, e = e, info = i))
    ei <- e[i]
    e[i] <- ei / d[i]
    d[i + 1L] <- d[i + 1L] - e[i] * ei
  }
  if (isTRUE(d[n] <= 0)) return(list(d = d, e = e, info = n))
  list(d = d, e = e, info = 0L)
}

# LAPACK dpttrs (dptts2) with one right hand side
.sptl_pttrs <- function(d, e, b) {
  n <- length(d)
  if (n == 1L) return(b * (1 / d[1L]))
  for (i in 2:n) b[i] <- b[i] - b[i - 1L] * e[i - 1L]
  b[n] <- b[n] / d[n]
  for (i in (n - 1L):1L) b[i] <- b[i] / d[i] - b[i + 1L] * e[i]
  b
}

# LAPACK dptrfs: the iterative refinement of x (the error bounds it also
# computes are discarded by trlib)
.sptl_ptrfs <- function(d, e, df, ef, b, x) {
  n <- length(d)
  eps <- 0.5 * .sptl_eps
  safe1 <- 4 * .Machine$double.xmin
  safe2 <- safe1 / eps
  count <- 1L
  lstres <- 3
  repeat {
    if (n == 1L) {
      dx <- d * x
      r <- b - dx
      w <- abs(b) + abs(dx)
    } else {
      cx <- c(0, e * x[-n])
      dx <- d * x
      ex <- c(e * x[-1L], 0)
      r <- b - cx - dx - ex
      w <- abs(b) + abs(cx) + abs(dx) + abs(ex)
    }
    s <- 0
    for (i in seq_len(n))
      s <- .sptl_fmax(s, if (isTRUE(w[i] > safe2)) abs(r[i]) / w[i] else (abs(r[i]) + safe1) / (w[i] + safe1))
    if (!(isTRUE(s > eps) && isTRUE(2 * s <= lstres) && count <= 5L)) break
    x <- x + .sptl_pttrs(df, ef, r)
    lstres <- s
    count <- count + 1L
  }
  x
}

# LAPACK dlagtm with alpha = beta = 1: T x + b for symmetric tridiagonal T
.sptl_lagtm <- function(d, e, x, b) {
  n <- length(d)
  if (n == 1L) return(b + d * x)
  b + c(0, e * x[-n]) + d * x + c(e * x[-1L], 0)
}

# trlib_quadratic_zero: the zeros of t^2 + c_lin t + c_abs with one Newton
# polish each
.sptl_quadratic_zero <- function(c_abs, c_lin, tol) {
  lin_sq <- c_lin * c_lin
  t1 <- 0
  t2 <- 0
  if (isTRUE(abs(c_abs) > tol * lin_sq)) {
    q <- lin_sq - 4 * c_abs
    if (isTRUE(abs(q) <= (.sptl_eps * c_lin) * (.sptl_eps * c_lin))) {
      t1 <- -0.5 * c_lin
      t2 <- t1
    } else if (isTRUE(q < 0)) {
      return(c(0, 0))
    } else {
      t1 <- -0.5 * (c_lin + .sptl_sign1(c_lin) * sqrt(q))
      t2 <- c_abs / t1
      if (isTRUE(t2 < t1)) {
        q <- t2
        t2 <- t1
        t1 <- q
      }
    }
  } else if (isTRUE(c_lin < 0)) {
    t2 <- -c_lin
  } else {
    t1 <- -c_lin
  }
  q <- (t1 + c_lin) * t1 + c_abs
  dq <- 2 * t1 + c_lin
  if (!isTRUE(dq == 0)) t1 <- t1 - q / dq
  q <- (t2 + c_lin) * t2 + c_abs
  dq <- 2 * t2 + c_lin
  if (!isTRUE(dq == 0)) t2 <- t2 - q / dq
  c(t1, t2)
}

# trlib_leftmost_irreducible, cold start (the only start SciPy reaches): the
# smallest eigenvalue of an irreducible tridiagonal matrix by a safeguarded
# Newton iteration on the last Parlett Reid pivot
.sptl_leftmost <- function(diag, offdiag, itmax, tol_abs) {
  n <- length(diag)
  if (n == 1L) return(diag[1L])
  oabs0 <- abs(offdiag[1L])
  oabs1 <- abs(offdiag[n - 1L])
  low <- .sptl_fmin(diag[1L] - oabs0, diag[n] - oabs1)
  # trlib writes diag[n-1] - oabs1 for the upper bound; it is still an upper
  # bound of the leftmost eigenvalue, kept as upstream
  up <- .sptl_fmax(diag[1L] + oabs0, diag[n] - oabs1)
  for (ii in seq_len(n - 2L) + 1L) {
    oabs1 <- abs(offdiag[ii])
    low <- .sptl_fmin(low, diag[ii] - oabs0 - oabs1)
    up <- .sptl_fmax(up, diag[ii] + oabs0 + oabs1)
    oabs0 <- oabs1
  }
  leftmost <- low + 0.1 * (up - low)
  itmax <- itmax * n
  iter_pr <- 0L
  repeat {
    iter_pr <- iter_pr + 1L
    if (iter_pr > itmax) return(leftmost)
    n_neg_piv <- 0L
    restart <- FALSE
    for (jj in seq_len(n)) {
      if (jj == 1L) {
        prlp <- diag[1L] - leftmost
        dprlp <- -1
        ddprlp <- 0
      } else {
        o <- offdiag[jj - 1L]
        obyprlp <- o / prlp
        ddprlp <- obyprlp * obyprlp * (ddprlp - 2 * dprlp * dprlp / prlp)
        dprlp <- -1 + o * o * dprlp / (prlp * prlp)
        prlp <- diag[jj] - leftmost - o * obyprlp
      }
      if (isTRUE(prlp == 0)) {
        if (n_neg_piv == 0L && jj == n) return(leftmost)
        up <- leftmost
        leftmost <- 0.5 * (low + up)
        restart <- TRUE
        break
      } else if (isTRUE(prlp < 0)) {
        n_neg_piv <- n_neg_piv + 1L
        up <- leftmost
        if (n_neg_piv > 1L) {
          leftmost <- 0.5 * (low + up)
          restart <- TRUE
          break
        }
      }
    }
    if (restart) next
    if (n_neg_piv == 0L) low <- leftmost else up <- leftmost
    tol <- tol_abs * .sptl_fmax(1, .sptl_fmax(abs(low), abs(up)))
    if (isTRUE(up - low <= tol) || isTRUE(abs(prlp) <= tol_abs)) return(leftmost)
    dleftmost <- -prlp / dprlp
    attempt <- leftmost + dleftmost
    if (isTRUE(low <= attempt) && isTRUE(attempt <= up)) {
      if (isTRUE(abs(dleftmost) <= tol)) return(leftmost)
      leftmost <- attempt
    } else {
      leftmost <- 0.5 * (low + up)
    }
  }
}

# trlib_eigen_inverse: inverse iteration for the eigenvector of the leftmost
# eigenvalue. trlib draws its start vectors from rand() seeded by the clock;
# the port uses fixed start vectors, which changes at most the sign of the
# eigenvector and not the trust region step built from it
.sptl_eigen_inverse <- function(S, diag, offdiag, lam_init, itmax = 10L, tol_abs = .sptl_eps_pow_5) {
  n <- length(diag)
  pert <- 0
  minuslam <- -lam_init
  repeat {
    f <- .sptl_pttrf(diag + minuslam, offdiag)
    if (f$info == 0L || pert > 1 / .sptl_eps) break
    pert <- if (pert == 0) .sptl_eps_pow_4 * .sptl_fmax(1, -lam_init) else 10 * pert
    minuslam <- pert - lam_init
    if (pert > 1 / .sptl_eps) break
  }
  S$delta_fac[seq_len(n)] <- f$d
  S$gamma_fac[seq_len(n - 1L)] <- f$e
  start <- function(jj) {
    eig <- 1 + ((seq_len(n) * (2L * jj - 1L)) %% 7L) / 7
    eig * (1 / .sptl_nrm2(eig))
  }
  iterate <- function(eig) {
    res <- Inf
    for (k in seq_len(itmax)) {
      eig <- .sptl_pttrs(f$d, f$e, eig)
      invnorm <- 1 / .sptl_nrm2(eig)
      eig <- eig * invnorm
      res <- abs(invnorm - pert)
      if (isTRUE(res <= tol_abs)) return(list(eig = eig, res = res, conv = TRUE))
    }
    list(eig = eig, res = res, conv = FALSE)
  }
  residuals <- numeric(5L)
  for (jj in 1:5) {
    r <- iterate(start(jj))
    if (r$conv) return(r$eig)
    residuals[jj] <- r$res
  }
  iterate(start(which.min(residuals)))$eig
}

# trlib_tri_factor_regularize_posdef
.sptl_regularize_posdef <- function(diag, offdiag, tol_away, security_step) {
  n <- length(diag)
  regdiag <- numeric(n)
  if (isTRUE(diag[1L] <= tol_away)) regdiag[1L] <- security_step * tol_away
  diag_fac <- diag[1L] + regdiag[1L]
  for (pivot in seq_len(n - 1L)) {
    o2 <- offdiag[pivot] * offdiag[pivot]
    if (isTRUE(diag[pivot + 1L] - o2 / diag_fac <= tol_away * diag_fac))
      regdiag[pivot + 1L] <- security_step * abs(o2 / diag_fac - diag[pivot + 1L])
    diag_fac <- diag[pivot + 1L] + regdiag[pivot + 1L] - o2 / diag_fac
  }
  regdiag
}

# trlib_tri_factor_min for one irreducible block: the Newton iteration on
# 1/||h(lam)|| - 1/radius for the tridiagonal subproblem of dimension n0,
# with the hard case handled by the leftmost eigenvector
.sptl_tri_factor_min <- function(S, n0, radius, itmax, warm_fac0) {
  tol_rel <- .sptl_eps
  tol_newton_tiny <- 1e-11
  idx <- seq_len(n0)
  idm <- seq_len(n0 - 1L)
  diag <- S$delta[idx]
  offdiag <- S$gamma[idm]
  neglin <- S$neglin[idx]
  S$iter_tri <- 0L
  S$h[idx] <- 0
  ret <- 0L
  norm_sol0 <- 0
  info_fac <- 0L
  sol0 <- S$h0[idx]
  setfac <- function(f) {
    S$delta_fac0[idx] <- f$d
    S$gamma_fac0[idm] <- f$e
  }
  if (S$warm_lam0) {
    if (!warm_fac0) {
      f <- .sptl_pttrf(diag + S$lam0, offdiag)
      setfac(f)
      if (f$info != 0L) S$warm_lam0 <- FALSE else warm_fac0 <- TRUE
    }
    if (warm_fac0) {
      sol0 <- .sptl_pttrs(S$delta_fac0[idx], S$gamma_fac0[idm], neglin)
      norm_sol0 <- .sptl_nrm2(sol0)
      S$warm_lam0 <- isTRUE(norm_sol0 >= radius)
    }
  }
  if (!S$warm_lam0) {
    S$lam0 <- 0
    sol0 <- neglin
    f <- .sptl_pttrf(diag, offdiag)
    setfac(f)
    info_fac <- f$info
    if (info_fac == 0L) {
      sol0 <- .sptl_pttrs(f$d, f$e, neglin)
      norm_sol0 <- .sptl_nrm2(sol0)
    }
    if (info_fac != 0L) {
      S$leftmost <- .sptl_leftmost(diag, offdiag, 1000L, .sptl_eps_pow_75)
      lm <- S$leftmost
      pert_low <- -.sptl_eps_pow_5 * abs(lm)
      pert_up <- 1 / .sptl_eps
      lam_pert <- 0
      for (jj in 1:50) {
        S$lam0 <- -lm + lam_pert
        f <- .sptl_pttrf(diag + S$lam0, offdiag)
        setfac(f)
        info_fac <- f$info
        if (info_fac == 0L) {
          pert_up <- lam_pert
          sol0 <- .sptl_pttrs(f$d, f$e, neglin)
          norm_sol0 <- .sptl_nrm2(sol0)
          if (isTRUE(norm_sol0 >= radius)) break
          lam_pert <- 0.5 * (pert_low + pert_up)
        } else {
          pert_low <- lam_pert
          if (isTRUE(pert_up == 1 / .sptl_eps)) {
            lam_pert <- if (isTRUE(lam_pert == 0)) (1 + abs(lm)) * .sptl_eps_pow_75 else 2 * lam_pert
          } else {
            lam_pert <- 0.5 * (pert_low + pert_up)
          }
        }
      }
      if (info_fac != 0L) {
        S$lam0 <- -lm + pert_up
        f <- .sptl_pttrf(diag + S$lam0, offdiag)
        setfac(f)
        if (f$info != 0L) return(-2L)
        # trlib solves with the previous contents of sol0 here, as upstream
        sol0 <- .sptl_pttrs(f$d, f$e, sol0)
        norm_sol0 <- .sptl_nrm2(sol0)
      }
    }
  }
  if (isTRUE(norm_sol0 >= radius)) {
    repeat {
      dfac <- S$delta_fac0[idx]
      efac <- S$gamma_fac0[idm]
      w <- sol0 * (1 / norm_sol0)
      inv_d <- w[1L] * w[1L] / dfac[1L]
      for (jj in seq_len(n0 - 1L) + 1L) {
        w[jj] <- w[jj] - efac[jj - 1L] * w[jj - 1L]
        inv_d <- inv_d + w[jj] * w[jj] / dfac[jj]
      }
      dlam <- (norm_sol0 - radius) / (radius * inv_d)
      S$iter_tri <- S$iter_tri + 1L
      if (isTRUE(S$lam0 + dlam == S$lam0) || isTRUE(abs(dlam) <= tol_newton_tiny * .sptl_fmax(1, abs(S$lam0))) ||
          S$iter_tri > itmax) {
        ret <- if (S$iter_tri > itmax) -1L else 3L
        break
      }
      S$lam0 <- S$lam0 + dlam
      diag_lam <- diag + S$lam0
      f <- .sptl_pttrf(diag_lam, offdiag)
      setfac(f)
      if (f$info != 0L) {
        S$h0[idx] <- sol0
        return(-2L)
      }
      sol0 <- .sptl_ptrfs(diag_lam, offdiag, f$d, f$e, neglin, .sptl_pttrs(f$d, f$e, neglin))
      norm_sol0 <- .sptl_nrm2(sol0)
      if (isTRUE(norm_sol0 - radius <= tol_rel * radius)) {
        ret <- 0L
        break
      }
    }
  }
  S$warm_lam0 <- TRUE
  if (isTRUE(abs(radius - norm_sol0) >= .sptl_eps_pow_5 * radius)) {
    if (isTRUE(S$lam0 == 0)) {
      ret <- 1L
    } else {
      lm <- S$leftmost
      eig <- .sptl_eigen_inverse(S, diag, offdiag, lm)
      dot <- sum(sol0 * eig)
      tt <- .sptl_quadratic_zero(norm_sol0 * norm_sol0 - radius * radius, 2 * dot, .sptl_eps_pow_75)
      dot2 <- sum(eig * neglin)
      lin <- lm * dot - dot2
      t <- if (isTRUE(0.5 * lm * tt[1L] * tt[1L] + lin * tt[1L] <= 0.5 * lm * tt[2L] * tt[2L] + lin * tt[2L]))
        tt[1L] else tt[2L]
      sol0 <- sol0 + t * eig
      ret <- 4L
    }
  }
  S$h0[idx] <- sol0
  S$lam <- S$lam0
  S$warm_lam <- TRUE
  S$h[idx] <- sol0
  w <- .sptl_lagtm(diag, offdiag, sol0, -2 * neglin)
  S$obj <- 0.5 * sum(sol0 * w)
  ret
}

# one transition of trlib_krylov_min_internal. trlib compares the pointer
# action with TRLIB_CLA_TRIVIAL at the end of its loop, so every call makes
# exactly one transition; the outer trlib_krylov_min loops over them.
# Returns codes: 0 conv bound, 1 conv interior, 2 approx hard, 4 hard init
# lam, 5 fail hard, 6 unbounded below, 7 unlikely convergence, 10 continue,
# -1 itmax, -5 fail numeric, -7 fail ttr, -8 pc indefinite, -9 unexpected
# interior. Actions: 0 trivial, 1 init, 2 retransform, 3 update stationary
# point, 4 update gradient, 5 update direction, 8 objective value. Types:
# 1 conjugate gradient, 2 Lanczos.
.sptl_internal <- function(S, io, init, radius, g_dot_g, v_dot_g, p_dot_Hp) {
  zero <- 2e-16
  obj_lo <- -1e20
  ret <- 10L
  io$iter <- S$ii
  if (init %in% c(1L, 2L, 4L)) S$status <- init
  ii <- S$ii
  i <- ii + 1L
  set <- function(ityp = io$ityp, action, status = S$status, flt1 = io$flt1, flt2 = io$flt2, flt3 = io$flt3) {
    io$ityp <- ityp
    io$action <- action
    S$status <- status
    io$flt1 <- flt1
    io$flt2 <- flt2
    io$flt3 <- flt3
  }
  cg_tri <- function() {
    warm_fac0 <- FALSE
    if (S$warm_lam0 && ii >= 1L) {
      warm_fac0 <- !isTRUE(S$delta_fac0[i - 1L] == 0)
      if (warm_fac0) {
        S$gamma_fac0[i - 1L] <- S$gamma[i - 1L] / S$delta_fac0[i - 1L]
        S$delta_fac0[i] <- S$delta[i] + S$lam0 - S$gamma[i - 1L] * S$gamma[i - 1L] / S$delta_fac0[i - 1L]
        warm_fac0 <- isTRUE(S$delta_fac0[i] > 0)
      }
    }
    S$irblk1 <- ii + 1L
    S$exit_tri <- .sptl_tri_factor_min(S, ii + 1L, radius, 100L + 3L * ii, warm_fac0)
  }
  early <- function() {
    ii > 10L && isTRUE(S$convhist[i - 10L] > 0.1 * S$convhist[i])
  }
  doit <- function() {
    k <- (ii - 10L):(ii - 1L) + 1L
    !any(S$convhist[k + 1L] > S$convhist[k], na.rm = TRUE)
  }
  st <- S$status
  if (st == 1L) {
    S$ii <- 0L
    io$iter <- 0L
    S$pos_def <- TRUE
    S$interior <- TRUE
    S$warm_lam0 <- FALSE
    S$warm_lam <- FALSE
    S$lanczos_switch <- -1L
    S$exit_tri <- 0L
    S$iter_tri <- 0L
    S$gamma[seq_len(S$itmax)] <- 0
    set(1L, 1L, 7L)
  } else if (st == 7L) {
    if (isTRUE(v_dot_g <= 0) && isTRUE(g_dot_g > 0)) {
      set(1L, 0L)
      ret <- -8L
    } else {
      eta_i <- 0
      if (S$tol_rel_i >= 0) eta_i <- S$tol_rel_i
      if (abs(S$tol_rel_i + 1) < 1e-8) eta_i <- .sptl_fmin(0.5, sqrt(sqrt(v_dot_g)))
      if (abs(S$tol_rel_i + 2) < 1e-8) eta_i <- .sptl_fmin(0.5, sqrt(v_dot_g))
      S$stop_i <- .sptl_fmax(0, eta_i * sqrt(v_dot_g))
      S$stop_i <- S$stop_i * S$stop_i
      eta_b <- 0
      if (S$tol_rel_b >= 0) eta_b <- S$tol_rel_b
      if (abs(S$tol_rel_b + 1) < 1e-8 || abs(S$tol_rel_b + 3) < 1e-8) eta_b <- .sptl_fmin(0.5, sqrt(sqrt(v_dot_g)))
      if (abs(S$tol_rel_b + 2) < 1e-8 || abs(S$tol_rel_b + 4) < 1e-8) eta_b <- .sptl_fmin(0.5, sqrt(v_dot_g))
      if (S$tol_rel_b < -2.5) eta_b <- max(1e-6, eta_b)
      S$stop_b <- .sptl_fmax(0, eta_b * sqrt(v_dot_g))
      S$v_g <- v_dot_g
      S$p_Hp <- p_dot_Hp
      S$neglin[1L] <- -sqrt(v_dot_g)
      S$cgl <- 1
      S$sigma <- 1
      S$leftmost <- 0
      S$lam <- 0
      S$obj <- 0
      S$s_Mp <- 0
      S$s_Ms <- 0
      S$p_Mp <- S$v_g
      S$raymax <- S$p_Hp / S$p_Mp
      S$raymin <- S$raymax
      S$delta[1L] <- 0
      set(1L, 0L, 8L)
    }
  } else if (st == 8L) {
    if (isTRUE(abs(S$p_Hp) <= zero)) {
      set(action = 0L, status = 12L)
    } else {
      S$alpha[i] <- S$v_g / S$p_Hp
      S$delta[i] <- S$delta[i] + S$p_Hp / S$v_g
      S$pos_def <- S$pos_def && isTRUE(S$alpha[i] > 0)
      if (ii > 0L) S$sigma <- -.sptl_sign1(S$alpha[i - 1L]) * S$sigma
      S$cgl <- S$sigma / sqrt(S$v_g)
      if (ii > 0L) {
        S$s_Mp <- S$beta[i - 1L] * (S$s_Mp + S$alpha[i - 1L] * S$p_Mp)
        S$p_Mp <- S$v_g + S$beta[i - 1L] * S$beta[i - 1L] * S$p_Mp
      }
      sp_Msp <- S$s_Ms + S$alpha[i] * (2 * S$s_Mp + S$alpha[i] * S$p_Mp)
      S$raymax <- .sptl_fmax(S$raymax, S$p_Hp / S$p_Mp)
      S$raymin <- .sptl_fmin(S$raymin, S$p_Hp / S$p_Mp)
      S$interior <- S$interior && S$pos_def && isTRUE(sp_Msp < radius * radius)
      if (S$interior) {
        S$s_Ms <- sp_Msp
        S$obj <- S$obj - 0.5 * S$alpha[i] * S$alpha[i] * S$p_Hp
        set(1L, 3L, 9L, flt1 = S$alpha[i])
      } else {
        cg_tri()
        if (S$exit_tri < 0L) {
          set(1L, 2L)
          ret <- -7L
        } else if (S$exit_tri == 1L && S$pos_def) {
          set(1L, 2L)
          ret <- -9L
        } else {
          set(1L, 0L, 9L)
        }
      }
    }
  } else if (st == 9L) {
    set(1L, 4L, 10L, flt1 = S$alpha[i], flt2 = S$cgl)
  } else if (st == 10L) {
    bad <- if (S$interior) 0L else 2L
    if (is.nan(v_dot_g)) {
      set(1L, bad)
      ret <- -5L
    } else if (isTRUE(g_dot_g > 0) && isTRUE(v_dot_g <= 0)) {
      set(1L, bad)
      ret <- -8L
    } else if (isTRUE(g_dot_g <= zero)) {
      set(1L, bad)
      ret <- 5L
      S$gamma[i] <- 0
    } else {
      S$beta[i] <- v_dot_g / S$v_g
      S$delta[i + 1L] <- (v_dot_g * S$p_Hp) / (S$v_g * S$v_g)
      S$gamma[i] <- abs(sqrt(v_dot_g) * S$p_Hp / (sqrt(S$v_g) * S$v_g))
      S$v_g <- v_dot_g
      S$convhist[i] <- if (S$interior) sqrt(S$v_g) / S$stop_i else S$gamma[i] * abs(S$h[i]) / S$stop_b
      if (S$interior && isTRUE(S$v_g <= S$stop_i)) {
        set(1L, 0L)
        ret <- 1L
      } else if (!S$interior && isTRUE(S$gamma[i] * abs(S$h[i]) <= S$stop_b)) {
        set(1L, 2L)
        ret <- 0L
      } else if (!S$interior && early()) {
        if (doit()) {
          set(1L, 2L)
          ret <- 7L
        } else {
          set(1L, 5L, 11L, flt1 = -1, flt2 = S$beta[i])
        }
      } else if (is.nan(S$obj) || isTRUE(S$obj < obj_lo)) {
        set(1L, bad)
        ret <- 6L
      } else {
        set(1L, 5L, 11L, flt1 = -1, flt2 = S$beta[i])
      }
    }
  } else if (st == 11L) {
    S$p_Hp <- p_dot_Hp
    S$ii <- ii + 1L
    if (S$ii >= S$itmax) {
      set(1L, if (S$interior) 0L else 2L)
      ret <- -1L
    } else {
      set(1L, 0L, 8L)
    }
  } else if (st == 2L) {
    S$irblk1 <- ii + 1L
    S$exit_tri <- .sptl_tri_factor_min(S, ii + 1L, radius, 100L + 3L * ii, FALSE)
    set(1L, 2L)
    # earlyterm: a hotstart only resolves the tridiagonal subproblem
    ret <- if (S$exit_tri < 0L) -7L else if (S$exit_tri == 4L) 4L else 7L
  } else if (st == 4L) {
    k <- seq_len(S$irblk1)
    S$regdiag[k] <- .sptl_regularize_posdef(S$delta[k], S$gamma[seq_len(S$irblk1 - 1L)], 1e-12, 10)
    S$delta[k] <- S$delta[k] + S$regdiag[k]
    S$h[seq_len(ii + 1L)] <- 0
    S$warm_lam0 <- FALSE
    S$warm_lam <- FALSE
    S$irblk1 <- ii + 1L
    S$exit_tri <- .sptl_tri_factor_min(S, ii + 1L, radius, 100L + 3L * ii, FALSE)
    # trlib removes the shift over the new block length
    k <- seq_len(S$irblk1)
    S$delta[k] <- S$delta[k] - S$regdiag[k]
    set(action = 2L)
    ret <- S$exit_tri
  } else if (st == 12L) {
    S$lanczos_switch <- ii
    if (ii > 0L) S$sigma <- -.sptl_sign1(S$alpha[i - 1L]) * S$sigma
    set(2L, 5L, 13L, flt1 = S$sigma / sqrt(S$v_g), flt2 = 0)
  } else if (st == 13L) {
    if (isTRUE(abs(p_dot_Hp) <= zero)) {
      set(1L, 2L)
      ret <- 5L
    } else {
      S$delta[i] <- p_dot_Hp
      S$raymax <- .sptl_fmax(S$raymax, p_dot_Hp)
      S$raymin <- .sptl_fmin(S$raymin, p_dot_Hp)
      cg_tri()
      if (S$exit_tri < 0L) {
        set(2L, 2L)
        ret <- -7L
      } else if (S$exit_tri == 4L) {
        set(2L, 2L)
        ret <- 4L
      } else {
        # gamma[ii - 1] at ii = 0 reads the word before gamma in trlib's work
        # array, the last entry of delta_fac
        gm1 <- if (ii >= 1L) S$gamma[i - 1L] else S$delta_fac[S$itmax + 1L]
        f1 <- -S$delta[i] / gm1
        f2 <- if (ii < 2L) -gm1 else -gm1 / S$gamma[i - 2L]
        f3 <- 1
        if (ii == S$lanczos_switch) {
          f1 <- -S$delta[i] / sqrt(S$v_g) * S$sigma
          f2 <- -gm1 * S$cgl
          f3 <- gm1 / sqrt(S$v_g)
          S$cgl <- 1
        }
        if (ii == 0L) {
          f1 <- -S$delta[i] / gm1
          f2 <- 0
          f3 <- 1
        }
        set(2L, 4L, 14L, flt1 = f1, flt2 = f2, flt3 = f3)
      }
    }
  } else if (st == 14L) {
    bad <- if (S$interior) 0L else 2L
    if (is.nan(v_dot_g)) {
      set(2L, bad)
      ret <- -5L
    } else if (isTRUE(v_dot_g <= 0) && isTRUE(g_dot_g > 0)) {
      set(2L, bad)
      ret <- -8L
    } else {
      if (isTRUE(g_dot_g <= zero)) ret <- 2L
      S$gamma[i] <- sqrt(v_dot_g)
      S$v_g <- v_dot_g * S$h[i] * S$h[i]
      set(2L, 0L, 16L)
    }
  } else if (st == 16L) {
    S$convhist[i] <- if (S$interior) sqrt(S$v_g) / sqrt(S$stop_i) else S$v_g / sqrt(S$stop_b)
    if (S$exit_tri != 1L && isTRUE(S$v_g <= sqrt(S$stop_b))) {
      set(2L, 2L)
      ret <- S$exit_tri
    } else if (S$exit_tri == 1L && isTRUE(S$v_g <= sqrt(S$stop_i))) {
      set(2L, 2L)
      ret <- S$exit_tri
    } else if (!S$interior && early()) {
      if (doit()) {
        set(2L, 2L)
        ret <- 7L
      } else {
        set(2L, 0L, 17L)
      }
    } else if (is.nan(S$obj) || isTRUE(S$obj < obj_lo)) {
      set(2L, 2L)
      ret <- 6L
    } else {
      set(2L, 0L, 17L)
    }
  } else if (st == 17L) {
    S$ii <- ii + 1L
    if (S$ii >= S$itmax || S$ii - S$lanczos_switch >= 100L) {
      set(2L, 2L)
      ret <- -1L
    } else {
      io$iter <- S$ii
      set(2L, 5L, 13L, flt1 = 1 / S$gamma[S$ii], flt2 = 0)
    }
  } else {
    io$action <- 0L
  }
  ret
}

# trlib_krylov_min: the outer loop over transitions with the deferred return
# codes and the convexification check
.sptl_krylov_min <- function(S, io, init, radius, g_dot_g, v_dot_g, p_dot_Hp) {
  ret <- -1000L
  io$iter <- S$ii
  if (init == 1L || init == 2L) S$outerstatus <- 0L
  if (S$outerstatus < 100L || S$outerstatus == 300L) {
    repeat {
      ret <- .sptl_internal(S, io, init, radius, g_dot_g, v_dot_g, p_dot_Hp)
      if (init > 0L || ret < 10L || io$action != 0L) break
    }
  }
  if (ret >= 0L || ret == -1000L) {
    if (S$outerstatus < 100L && ret < 10L && io$action != 0L) {
      S$outerstatus <- 100L + ret
      return(10L)
    }
    if (S$outerstatus >= 100L && S$outerstatus < 200L) {
      ret <- S$outerstatus - 100L
      S$outerstatus <- 0L
      io$action <- 0L
    }
    if (ret < 10L && S$outerstatus < 100L &&
        isTRUE(S$lam > 1e-2 * .sptl_fmax(1, S$raymax)) && isTRUE(S$raymin < 0) &&
        isTRUE(abs(S$raymin) < 1e-8 * S$raymax)) {
      S$outerstatus <- 200L + ret
      io$action <- 8L
      return(10L)
    }
    if (S$outerstatus > 200L && S$outerstatus < 300L) {
      if (isTRUE(abs(S$obj - g_dot_g) > .sptl_fmax(1e-6, 1e-1 * abs(g_dot_g))) || isTRUE(g_dot_g > 0)) {
        ret <- .sptl_internal(S, io, 4L, radius, g_dot_g, v_dot_g, p_dot_Hp)
        S$outerstatus <- 300L
        return(ret)
      }
      ret <- S$outerstatus - 200L
      S$outerstatus <- 0L
      return(ret)
    }
    if (S$outerstatus == 300L && ret < 10L) {
      S$outerstatus <- 0L
      return(ret)
    }
  }
  ret
}

.sptl_new <- function(n, trlib) {
  itmax <- as.integer(min(1e9 / n, 2 * n))
  S <- new.env(parent = emptyenv())
  S$itmax <- itmax
  S$tol_rel_i <- trlib$tol_rel_i
  S$tol_rel_b <- trlib$tol_rel_b
  len <- itmax + 2L
  for (nm in c("alpha", "beta", "neglin", "h0", "h", "delta", "delta_fac0", "delta_fac", "gamma", "gamma_fac0",
               "gamma_fac", "convhist", "regdiag"))
    assign(nm, numeric(len), envir = S)
  S$status <- 0L
  S$ii <- 0L
  S$outerstatus <- 0L
  S$irblk1 <- 0L
  S$warm_lam0 <- FALSE
  S$warm_lam <- FALSE
  S$lam0 <- 0
  S$lam <- 0
  S$leftmost <- 0
  S$raymax <- 0
  S$raymin <- 0
  S$obj <- 0
  S$init <- 1L
  S$Q <- matrix(0, itmax + 1L, n)
  S$s <- numeric(n)
  S$g <- numeric(n)
  S$v <- numeric(n)
  S$gm <- numeric(n)
  S$p <- numeric(n)
  S$Hp <- numeric(n)
  S
}

# TRLIBQuadraticSubproblem.solve of _trlib.pyx: answers the reverse
# communication requests of trlib with vectors of the subproblem m. The trlib
# state lives in m, so a rejected step hotstarts from the previous solve.
.sp_trlib_solve <- function(m, trust_radius, trlib) {
  jac <- m$jac()
  S <- m$trlib_state
  if (is.null(S)) {
    S <- .sptl_new(length(jac), trlib)
    m$trlib_state <- S
  }
  io <- new.env(parent = emptyenv())
  io$action <- 0L
  io$iter <- 0L
  io$ityp <- 0L
  io$flt1 <- 0
  io$flt2 <- 0
  io$flt3 <- 0
  g_dot_g <- 0
  v_dot_g <- 0
  p_dot_Hp <- 0
  init <- S$init
  repeat {
    ret <- .sptl_krylov_min(S, io, init, trust_radius, g_dot_g, v_dot_g, p_dot_Hp)
    init <- 0L
    a <- io$action
    it <- io$iter
    if (a == 1L) {
      S$s[] <- 0
      S$gm[] <- 0
      S$g <- jac
      S$v <- S$g
      g_dot_g <- sum(S$g * S$g)
      v_dot_g <- sum(S$v * S$g)
      S$p <- -S$v
      S$Hp <- m$hessp(S$p)
      p_dot_Hp <- sum(S$p * S$Hp)
      S$Q[1L, ] <- S$v / sqrt(v_dot_g)
    } else if (a == 2L) {
      s <- numeric(length(jac))
      for (k in seq_len(it + 1L)) s <- s + S$h[k] * S$Q[k, ]
      S$s <- s
    } else if (a == 3L) {
      if (io$ityp == 1L) S$s <- S$s + io$flt1 * S$p
    } else if (a == 4L) {
      if (io$ityp == 1L) {
        S$Q[it + 1L, ] <- io$flt2 * S$v
        S$gm <- S$g
        S$g <- S$g + io$flt1 * S$Hp
      }
      if (io$ityp == 2L) {
        S$s <- S$Hp + io$flt1 * S$g + io$flt2 * S$gm
        S$gm <- io$flt3 * S$g
        S$g <- S$s
      }
      S$v <- S$g
      g_dot_g <- sum(S$g * S$g)
      v_dot_g <- sum(S$v * S$g)
    } else if (a == 5L) {
      S$p <- io$flt1 * S$v + io$flt2 * S$p
      S$Hp <- m$hessp(S$p)
      p_dot_Hp <- sum(S$p * S$Hp)
      if (io$ityp == 2L) S$Q[it + 1L, ] <- S$p
    } else if (a == 8L) {
      g_dot_g <- 0.5 * sum(S$s * m$hessp(S$s)) + sum(S$s * jac)
    }
    if (ret < 10L) break
    S$init <- 2L
  }
  list(p = S$s, hits = isTRUE(S$lam > 0))
}
