# The Krylov solvers that the SciPy ports use (scipy/sparse/linalg/_isolve of
# SciPy 1.18.1): lsmr (lsmr.py) for tr_solver = "lsmr" of least_squares, and
# lgmres with the flexible GMRES Arnoldi process of _gcrotmk.py for the
# KrylovJacobian of root(method = "krylov") (_nonlin.py). Linear operators are
# passed as functions computing A v (and A' v for lsmr).

.spk_norm <- function(x) sqrt(sum(x * x))

# _sym_ortho of lsqr.py: a stable Givens rotation, c(c, s, r)
.spk_sym_ortho <- function(a, b) {
  if (b == 0) return(c(sign(a), 0, abs(a)))
  if (a == 0) return(c(0, sign(b), abs(b)))
  if (abs(b) > abs(a)) {
    tau <- a / b
    s <- sign(b) / sqrt(1 + tau * tau)
    cc <- s * tau
    r <- b / s
  } else {
    tau <- b / a
    cc <- sign(a) / sqrt(1 + tau * tau)
    s <- cc * tau
    r <- a / cc
  }
  c(cc, s, r)
}

# lsmr(A, b, damp, atol, btol, conlim, maxiter) with A of shape m x n given by
# mv(v) = A v and rmv(u) = A' u; returns list(x, istop, itn, normr, normar,
# normA, condA, normx)
.spk_lsmr <- function(mv, rmv, m, n, b, damp = 0, atol = 1e-6, btol = 1e-6, conlim = 1e8, maxiter = NULL,
                      show = FALSE) {
  if (is.null(maxiter)) maxiter <- min(m, n)
  u <- b
  normb <- .spk_norm(b)
  x <- numeric(n)
  beta <- normb
  if (beta > 0) {
    u <- (1 / beta) * u
    v <- rmv(u)
    alpha <- .spk_norm(v)
  } else {
    v <- numeric(n)
    alpha <- 0
  }
  if (alpha > 0) v <- (1 / alpha) * v
  itn <- 0L
  zetabar <- alpha * beta
  alphabar <- alpha
  rho <- 1; rhobar <- 1; cbar <- 1; sbar <- 0
  h <- v
  hbar <- numeric(n)
  betadd <- beta; betad <- 0; rhodold <- 1; tautildeold <- 0; thetatilde <- 0; zeta <- 0; d <- 0
  normA2 <- alpha * alpha
  maxrbar <- 0; minrbar <- 0x1.249ad2594c37dp+332
  normA <- sqrt(normA2)
  condA <- 1; normx <- 0
  istop <- 0L
  ctol <- if (conlim > 0) 1 / conlim else 0
  normr <- beta
  normar <- alpha * beta
  out <- function() list(x = x, istop = istop, itn = itn, normr = normr, normar = normar, normA = normA,
                         condA = condA, normx = normx)
  if (normar == 0) return(out())
  if (normb == 0) {
    x[] <- 0
    return(out())
  }
  while (itn < maxiter) {
    itn <- itn + 1L
    u <- u * (-alpha) + mv(v)
    beta <- .spk_norm(u)
    if (beta > 0) {
      u <- u * (1 / beta)
      v <- v * (-beta) + rmv(u)
      alpha <- .spk_norm(v)
      if (alpha > 0) v <- v * (1 / alpha)
    }
    so <- .spk_sym_ortho(alphabar, damp)
    chat <- so[1L]; shat <- so[2L]; alphahat <- so[3L]
    rhoold <- rho
    so <- .spk_sym_ortho(alphahat, beta)
    cc <- so[1L]; s <- so[2L]; rho <- so[3L]
    thetanew <- s * alpha
    alphabar <- cc * alpha
    rhobarold <- rhobar
    zetaold <- zeta
    thetabar <- sbar * rho
    rhotemp <- cbar * rho
    so <- .spk_sym_ortho(cbar * rho, thetanew)
    cbar <- so[1L]; sbar <- so[2L]; rhobar <- so[3L]
    zeta <- cbar * zetabar
    zetabar <- -sbar * zetabar
    hbar <- hbar * (-(thetabar * rho / (rhoold * rhobarold))) + h
    x <- x + (zeta / (rho * rhobar)) * hbar
    h <- h * (-(thetanew / rho)) + v
    betaacute <- chat * betadd
    betacheck <- -shat * betadd
    betahat <- cc * betaacute
    betadd <- -s * betaacute
    thetatildeold <- thetatilde
    so <- .spk_sym_ortho(rhodold, thetabar)
    ctildeold <- so[1L]; stildeold <- so[2L]; rhotildeold <- so[3L]
    thetatilde <- stildeold * rhobar
    rhodold <- ctildeold * rhobar
    betad <- -stildeold * betad + ctildeold * betahat
    tautildeold <- (zetaold - thetatildeold * tautildeold) / rhotildeold
    taud <- (zeta - thetatilde * tautildeold) / rhodold
    d <- d + betacheck * betacheck
    normr <- sqrt(d + (betad - taud)^2 + betadd * betadd)
    normA2 <- normA2 + beta * beta
    normA <- sqrt(normA2)
    normA2 <- normA2 + alpha * alpha
    maxrbar <- max(maxrbar, rhobarold)
    if (itn > 1L) minrbar <- min(minrbar, rhobarold)
    condA <- max(maxrbar, rhotemp) / min(minrbar, rhotemp)
    normar <- abs(zetabar)
    normx <- .spk_norm(x)
    test1 <- normr / normb
    test2 <- if ((normA * normr) != 0) normar / (normA * normr) else Inf
    test3 <- 1 / condA
    t1 <- test1 / (1 + normA * normx / normb)
    rtol <- btol + atol * normA * normx / normb
    if (itn >= maxiter) istop <- 7L
    if (1 + test3 <= 1) istop <- 6L
    if (1 + test2 <= 1) istop <- 5L
    if (1 + t1 <= 1) istop <- 4L
    if (test3 <= ctol) istop <- 3L
    if (test2 <= atol) istop <- 2L
    if (test1 <= rtol) istop <- 1L
    if (istop > 0L) break
  }
  out()
}

# dlartg of LAPACK 3.10 (the convention of the BLAS SciPy links): c(c, s, r)
.spk_lartg <- function(f, g) {
  if (g == 0) return(c(1, 0, f))
  if (f == 0) return(c(0, sign(g), abs(g)))
  d <- sqrt(f * f + g * g)
  r <- if (f < 0) -d else d
  c(abs(f) / d, g / r, r)
}

# _fgmres of _gcrotmk.py (no C/U columns, no right preconditioner): the QR
# factorization of the Hessenberg matrix is updated by qr_insert, that is a
# Givens rotation of the last two rows per new column
.spk_fgmres <- function(matvec, v0, m, atol, lpsolve = identity, outer_v = list(), prepend_outer_v = FALSE) {
  vs <- list(v0)
  zs <- list()
  nout <- length(outer_v)
  m <- m + nout
  Q <- matrix(1, 1, 1)
  R <- matrix(0, 1, 0)
  breakdown <- FALSE
  res <- NaN
  for (j in seq_len(m) - 1L) {
    w <- NULL
    if (prepend_outer_v && j < nout) {
      z <- outer_v[[j + 1L]]$z
      w <- outer_v[[j + 1L]]$w
    } else if (prepend_outer_v && j == nout) {
      z <- v0
    } else if (!prepend_outer_v && j >= m - nout) {
      z <- outer_v[[j - (m - nout) + 1L]]$z
      w <- outer_v[[j - (m - nout) + 1L]]$w
    } else {
      z <- vs[[length(vs)]]
    }
    if (is.null(w)) w <- lpsolve(matvec(z))
    w_norm <- .spk_norm(w)
    hcur <- numeric(j + 2L)
    for (i in seq_along(vs)) {
      alpha <- sum(vs[[i]] * w)
      hcur[i] <- alpha
      w <- w - alpha * vs[[i]]
    }
    hcur[j + 2L] <- .spk_norm(w)
    alpha <- 1 / hcur[j + 2L]
    if (is.finite(alpha)) w <- alpha * w
    if (!(hcur[j + 2L] > .Machine$double.eps * w_norm)) breakdown <- TRUE
    vs[[length(vs) + 1L]] <- w
    zs[[length(zs) + 1L]] <- z
    # qr_insert(Q2, R2, hcur, j, which = "col")
    Q2 <- matrix(0, j + 2L, j + 2L)
    Q2[seq_len(j + 1L), seq_len(j + 1L)] <- Q
    Q2[j + 2L, j + 2L] <- 1
    R2 <- matrix(0, j + 2L, j + 1L)
    R2[seq_len(j + 1L), seq_len(j)] <- R
    col <- as.numeric(crossprod(Q2, hcur))
    g <- .spk_lartg(col[j + 1L], col[j + 2L])
    col[j + 1L] <- g[3L]
    col[j + 2L] <- 0
    R2[, j + 1L] <- col
    qa <- Q2[, j + 1L]; qb <- Q2[, j + 2L]
    Q2[, j + 1L] <- g[1L] * qa + g[2L] * qb
    Q2[, j + 2L] <- g[1L] * qb - g[2L] * qa
    Q <- Q2
    R <- R2
    res <- abs(Q[1L, j + 2L])
    if (res < atol || breakdown) break
  }
  if (!is.finite(R[j + 1L, j + 1L])) .sp_stop("sp_linalg", "LinAlgError")
  k <- j + 1L
  # lstsq (gelsd): the minimum norm solution, also for a singular R
  y <- .sp_lstsq(R[seq_len(k), seq_len(k), drop = FALSE], Q[1L, seq_len(k)])
  list(Q = Q, R = R, vs = vs, zs = zs, y = y, res = res)
}

# lgmres(A, b, x0 = 0, rtol, atol, maxiter, inner_m, outer_k, outer_v,
# store_outer_Av, prepend_outer_v) without preconditioner; the augmentation
# vectors live in the environment `store` (store$outer_v) and are kept
# between calls, as the list SciPy mutates
.spk_lgmres <- function(matvec, b, rtol = 1e-5, atol = 0, maxiter = 1000L, inner_m = 30L, outer_k = 3L,
                        store = NULL, store_outer_Av = TRUE, prepend_outer_v = FALSE) {
  if (!all(is.finite(b))) stop("RHS must contain only finite numbers", call. = FALSE)
  if (is.null(store)) {
    store <- new.env(parent = emptyenv())
    store$outer_v <- list()
  }
  x <- numeric(length(b))
  b_norm <- .spk_norm(b)
  atol <- max(atol, rtol * b_norm)
  if (b_norm == 0) return(list(x = b, info = 0L))
  ptol_max_factor <- 1
  for (k_outer in seq_len(maxiter)) {
    r_outer <- matvec(x) - b
    r_norm <- .spk_norm(r_outer)
    if (r_norm <= max(atol, rtol * b_norm)) return(list(x = x, info = 0L))
    v0 <- -r_outer
    inner_res_0 <- .spk_norm(v0)
    if (inner_res_0 == 0)
      stop("Preconditioner returned a zero vector; |v| ~ ", signif(r_norm, 1), ", |M v| = 0", call. = FALSE)
    v0 <- (1 / inner_res_0) * v0
    ptol <- min(ptol_max_factor, max(atol, rtol * b_norm) / r_norm)
    fg <- tryCatch({
      r <- .spk_fgmres(matvec, v0, inner_m, atol = ptol, outer_v = store$outer_v, prepend_outer_v = prepend_outer_v)
      r$y <- r$y * inner_res_0
      if (!all(is.finite(r$y))) .sp_stop("sp_linalg", "LinAlgError")
      r
    }, sp_linalg = function(e) NULL)
    if (is.null(fg)) return(list(x = x, info = k_outer))
    if (fg$res > ptol) ptol_max_factor <- min(1, 1.5 * ptol_max_factor)
    else ptol_max_factor <- max(1e-16, 0.25 * ptol_max_factor)
    y <- fg$y
    dx <- fg$zs[[1L]] * y[1L]
    for (i in seq_along(y)[-1L]) dx <- dx + y[i] * fg$zs[[i]]
    nx <- .spk_norm(dx)
    if (nx > 0) {
      if (store_outer_Av) {
        q <- as.numeric(fg$Q %*% (fg$R %*% y))
        ax <- fg$vs[[1L]] * q[1L]
        for (i in seq_along(q)[-1L]) ax <- ax + q[i] * fg$vs[[i]]
        store$outer_v[[length(store$outer_v) + 1L]] <- list(z = dx / nx, w = ax / nx)
      } else {
        store$outer_v[[length(store$outer_v) + 1L]] <- list(z = dx / nx, w = NULL)
      }
    }
    while (length(store$outer_v) > outer_k) store$outer_v[[1L]] <- NULL
    x <- x + dx
  }
  list(x = x, info = as.integer(maxiter))
}

# KrylovJacobian of _nonlin.py with the lgmres inner solver: the Jacobian
# vector products are forward differences of the function
.spk_krylov_jacobian <- function(opts) {
  known <- c("rdiff", "method", "inner_maxiter", "inner_M", "outer_k")
  inner <- setdiff(names(opts), known)
  for (key in inner) {
    if (!startsWith(key, "inner_")) stop("Unknown parameter ", key, call. = FALSE)
    if (!(substring(key, 7L) %in% c("rtol", "atol", "maxiter", "inner_m", "outer_k", "outer_v", "store_outer_Av",
                                    "prepend_outer_v", "callback", "x0", "M"))) {
      warning("Option '", key, "' is invalid for the inner method: lgmres. It will be ignored.", call. = FALSE)
      opts[[key]] <- NULL
    }
  }
  method <- opts$method %||% "lgmres"
  if (!identical(method, "lgmres"))
    stop("KrylovJacobian(method = '", method, "') is not available; the port uses lgmres", call. = FALSE)
  if (!is.null(opts$inner_M)) stop("KrylovJacobian: a preconditioner inner_M is not supported", call. = FALSE)
  kw <- list(maxiter = 1L, outer_k = opts$outer_k %||% 10L, prepend_outer_v = TRUE, store_outer_Av = FALSE,
             atol = 0, inner_m = 30L)
  for (key in intersect(inner, names(opts))) kw[[substring(key, 7L)]] <- opts[[key]]
  J <- new.env(parent = emptyenv())
  J$store <- new.env(parent = emptyenv())
  J$store$outer_v <- list()
  J$rdiff <- opts$rdiff
  diff_step <- function() {
    J$omega <- J$rdiff * max(1, max(abs(J$x0))) / max(1, max(abs(J$f0)))
  }
  J$matvec <- function(v) {
    nv <- .spk_norm(v)
    if (nv == 0) return(0 * v)
    sc <- J$omega / nv
    r <- (J$func(J$x0 + sc * v) - J$f0) / sc
    if (!all(is.finite(r)) && all(is.finite(v))) stop("Function returned non-finite results", call. = FALSE)
    r
  }
  J$setup <- function(x, F, func) {
    J$func <- func
    J$x0 <- x
    J$f0 <- F
    if (is.null(J$rdiff)) J$rdiff <- .Machine$double.eps^(1 / 2)
    diff_step()
  }
  J$solve <- function(f, tol = 0) {
    rtol <- if (!is.null(kw$rtol)) kw$rtol else tol
    .spk_lgmres(J$matvec, f, rtol = rtol, atol = kw$atol, maxiter = kw$maxiter, inner_m = kw$inner_m,
                outer_k = kw$outer_k, store = J$store, store_outer_Av = kw$store_outer_Av,
                prepend_outer_v = kw$prepend_outer_v)$x
  }
  J$update <- function(x, f) {
    J$x0 <- x
    J$f0 <- f
    diff_step()
  }
  J
}

# numpy.roots: the eigenvalues of the companion matrix
.spk_roots <- function(p) {
  nz <- which(p != 0)
  if (!length(nz)) return(complex())
  trailing <- length(p) - nz[length(nz)]
  p <- p[nz[1L]:nz[length(nz)]]
  N <- length(p)
  r <- complex()
  if (N > 1L) {
    A <- matrix(0, N - 1L, N - 1L)
    if (N > 2L) A[cbind(2:(N - 1L), 1:(N - 2L))] <- 1
    A[1L, ] <- -p[-1L] / p[1L]
    r <- as.complex(eigen(A, symmetric = FALSE, only.values = TRUE)$values)
  }
  c(r, complex(trailing))
}

# solve_trust_region_2d of _lsq/common.py
.spk_tr_2d <- function(B, g, Delta) {
  R <- tryCatch(chol(B), error = function(e) NULL)
  if (!is.null(R)) {
    p <- -backsolve(R, forwardsolve(t(R), g))
    if (sum(p * p) <= Delta^2) return(list(p = p, newton = TRUE))
  }
  a <- B[1L, 1L] * Delta^2
  b <- B[1L, 2L] * Delta^2
  cc <- B[2L, 2L] * Delta^2
  d <- g[1L] * Delta
  f <- g[2L] * Delta
  t <- .spk_roots(c(-b + d, 2 * (a - cc + f), 6 * b, 2 * (-a + cc + f), -b - d))
  t <- Re(t[Im(t) == 0])
  P <- Delta * rbind(2 * t / (1 + t^2), (1 - t^2) / (1 + t^2))
  value <- 0.5 * colSums(P * (B %*% P)) + as.numeric(crossprod(g, P))
  list(p = P[, which.min(value)], newton = FALSE)
}
