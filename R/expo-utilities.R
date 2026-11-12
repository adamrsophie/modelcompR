# Matrix exponentials, phi functions and Krylov subspace methods, after
# ExponentialUtilities.jl 1.37 (exp_generic.jl, exp_baseexp.jl, phi.jl,
# phi_almohy.jl, arnoldi.jl, krylov_phiv.jl, krylov_phiv_adaptive.jl),
# PureGebal.jl (balancing) and the dense exp! of Julia's LinearAlgebra.
#
# These are the numerical kernels of the exponential integrators
# (R/solvers-linear-lie.R, R/solvers-exponential.R, R/solvers-rkip.R). A
# linear operator `A` is a matrix or a SciML operator (R/sciml-operators.R).

# A %*% x for a matrix or an operator
.eu_mul <- function(A, x) {
  if (is.matrix(A)) return(as.numeric(A %*% x))
  as.numeric(A$apply(as.numeric(x)))
}

# isapprox(x, y) with the default relative tolerance sqrt(eps)
.eu_isapprox <- function(x, y)
  x == y || (is.finite(x) && is.finite(y) && abs(x - y) <= sqrt(.Machine$double.eps) * max(abs(x), abs(y)))

.eu_matrix <- function(A) if (is.matrix(A)) A else .sop_matrix(A)

# LinearAlgebra.ishermitian: exact symmetry
.eu_ishermitian <- function(A) {
  if (!is.matrix(A)) {
    if (inherits(A, "DiagonalOperator")) return(TRUE)
    if (is.null(A$as_matrix)) return(FALSE)
    A <- A$as_matrix()
  }
  nrow(A) == ncol(A) && all(A == t(A))
}

.eu_opnorm1 <- function(A) max(colSums(abs(A)))
.eu_opnorminf <- function(A) max(rowSums(abs(A)))
.eu_norm2 <- function(x) sqrt(sum(x * x))

# intlog2(x) of exp_generic.jl: the bit length of ceil(x) - 1
.eu_intlog2 <- function(x) {
  v <- ceiling(x) - 1
  s <- 0
  while (v >= 1) {
    v <- floor(v / 2)
    s <- s + 1
  }
  s
}

# right division A / B (rdiv! with the LU factorization of B)
.eu_rdiv <- function(A, B) t(solve(t(B), t(A)))

# ---------------------------------------------------------------------------
# ExpMethodGeneric: scaling and squaring with the [13/13] Pade approximant in
# Horner form, p(x) / p(-x)

.EU_GENERIC_C <- c(2.8101705462199623e-15, 0x1.1cc1e2df80824p-42, 1.48377004840414e-11,
                   0x1.5aad6134c4c94p-31, 2.0431513566525008e-8, 5.175983436853002e-7,
                   1.0351966873706003e-5, 0.00016304347826086958, 0.0019927536231884057,
                   0.018333333333333333, 0.12, 0.5)

.eu_exp_generic <- function(x) {
  x <- as.matrix(x)
  nx <- .eu_opnorm1(x)
  if (is.nan(nx) || nx > 4611686018427387904) return(x * sum(x * exp(nx)))
  s <- if (nx > 1) .eu_intlog2(nx) else 0
  padep <- function(z) {
    y <- z * 1.5440497506703088e-17
    for (cf in .EU_GENERIC_C) {
      diag(y) <- diag(y) + cf
      y <- y %*% z
    }
    diag(y) <- diag(y) + 1
    y
  }
  if (s > 0) x <- x * (1 / 2^s)
  q <- padep(-x)
  y <- .eu_rdiv(padep(x), q)
  for (i in seq_len(s)) y <- y %*% y
  y
}

# ---------------------------------------------------------------------------
# Balancing (xGEBAL with job 'B': permutation, then scaling by powers of 2)

.eu_gebal <- function(A) {
  n <- nrow(A)
  s <- rep(1, n)
  if (n == 0L) return(list(A = A, scale = s, ilo = 1L, ihi = 0L))
  if (any(!is.finite(A))) stop("matrix contains Infs or NaNs", call. = FALSE)
  # isolate eigenvalues in the trailing rows
  ihi <- n + 1L
  while (ihi > 1L) {
    ihi <- ihi - 1L
    exchange <- FALSE
    source <- 0L
    for (j in ihi:1L) {
      exchange <- TRUE
      source <- j
      for (i in seq_len(ihi)) if (i != j && A[j, i] != 0) {
        exchange <- FALSE
        break
      }
      if (exchange) break
    }
    if (!exchange) break
    s[ihi] <- source
    if (source != ihi) {
      r <- seq_len(ihi)
      tmp <- A[r, source]; A[r, source] <- A[r, ihi]; A[r, ihi] <- tmp
      tmp <- A[source, ]; A[source, ] <- A[ihi, ]; A[ihi, ] <- tmp
    }
  }
  ilo <- 1L
  if (ihi > 1L) {
    # isolate eigenvalues in the leading columns
    ilo <- 0L
    while (ilo < n) {
      ilo <- ilo + 1L
      exchange <- FALSE
      source <- 0L
      if (ilo <= ihi) for (j in ilo:ihi) {
        exchange <- TRUE
        source <- j
        for (i in ilo:ihi) if (i != j && A[i, j] != 0) {
          exchange <- FALSE
          break
        }
        if (exchange) break
      }
      if (!exchange) break
      s[ilo] <- source
      if (source != ilo) {
        r <- seq_len(ihi)
        tmp <- A[r, source]; A[r, source] <- A[r, ilo]; A[r, ilo] <- tmp
        cc <- ilo:n
        tmp <- A[source, cc]; A[source, cc] <- A[ilo, cc]; A[ilo, cc] <- tmp
      }
    }
  }
  # scale the remaining block
  sfmin2 <- 2 * (0x1.0000000000000p-1022 / 2.220446049250313e-16)
  sfmax2 <- 1 / sfmin2
  repair <- function(nrm, a, cm) {
    if (!is.finite(nrm) || (nrm == 0 && cm != 0)) {
      e <- floor(log2(cm))
      nrm <- sqrt(sum((a / 2^e)^2)) * 2^e
    }
    nrm
  }
  converged <- FALSE
  while (!converged && ilo <= ihi) {
    converged <- TRUE
    blk <- ilo:ihi
    for (i in blk) {
      col <- A[blk, i]
      row <- A[i, blk]
      ccmax <- max(abs(col))
      rcmax <- max(abs(row))
      cnorm <- repair(sqrt(sum(col * col)), col, ccmax)
      rnorm <- repair(sqrt(sum(row * row)), row, rcmax)
      cmax <- sqrt(max(col * col))
      rmax <- sqrt(max(row * row))
      if (!is.finite(cmax) || (cmax == 0 && ccmax != 0)) cmax <- ccmax
      if (!is.finite(rmax) || (rmax == 0 && rcmax != 0)) rmax <- rcmax
      if (cnorm == 0 || rnorm == 0) next
      g <- rnorm / 2
      osum <- cnorm + rnorm
      f <- 1
      while (cnorm < rnorm / 2) {
        if (cnorm >= g || max(f, cnorm, cmax) >= sfmax2 || min(rnorm, g, rmax) <= sfmin2) break
        f <- f * 2
        cnorm <- cnorm * 2
        cmax <- cmax * 2
        rnorm <- rnorm / 2
        g <- g / 2
        rmax <- rmax / 2
      }
      g <- cnorm / 2
      while (rnorm <= cnorm / 2) {
        if (g < rnorm || max(rnorm, rmax) >= sfmax2 || min(f, cnorm, g, cmax) <= sfmin2) break
        f <- f / 2
        cnorm <- cnorm / 2
        g <- g / 2
        cmax <- cmax / 2
        rnorm <- rnorm * 2
        rmax <- rmax * 2
      }
      if (cnorm + rnorm >= 0.95 * osum) next
      converged <- FALSE
      s[i] <- s[i] * f
      A[i, ilo:n] <- A[i, ilo:n] * (1 / f)
      A[seq_len(ihi), i] <- A[seq_len(ihi), i] * f
    }
  }
  list(A = A, scale = s, ilo = ilo, ihi = ihi)
}

.eu_rcswap <- function(X, i, j) {
  if (i == j) return(X)
  tmp <- X[i, ]; X[i, ] <- X[j, ]; X[j, ] <- tmp
  tmp <- X[, i]; X[, i] <- X[, j]; X[, j] <- tmp
  X
}

.eu_gebak <- function(X, bal) {
  n <- nrow(X)
  s <- bal$scale
  ilo <- bal$ilo
  ihi <- bal$ihi
  if (ilo <= ihi) for (j in ilo:ihi) {
    X[j, ] <- X[j, ] * s[j]
    X[, j] <- X[, j] / s[j]
  }
  if (ilo > 1L) for (j in (ilo - 1L):1L) X <- .eu_rcswap(X, j, as.integer(s[j]))
  if (ihi < n) for (j in (ihi + 1L):n) X <- .eu_rcswap(X, j, as.integer(s[j]))
  X
}

# ---------------------------------------------------------------------------
# Pade approximants of Higham (2005)

.EU_PADE <- list(
  c3 = c(120, 60, 12, 1),
  c5 = c(30240, 15120, 3360, 420, 30, 1),
  c7 = c(17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1),
  c9 = c(17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1),
  c13 = c(64764752532480000, 32382376266240000, 7771770303897600, 1187353796428800, 129060195264000,
          10559470521600, 670442572800, 33522128640, 1323241920, 40840800, 960960, 16380, 182, 1))

# U and V by powers of A^2 (_pade_evaluate! of exp_baseexp.jl), then (V - U) \ (V + U)
.eu_pade_loop <- function(A, C) {
  n <- nrow(A)
  A2 <- A %*% A
  P <- diag(n)
  U <- C[2L] * P
  V <- C[1L] * P
  for (k in seq_len(length(C) %/% 2L - 1L)) {
    P <- P %*% A2
    U <- U + C[2L * k + 2L] * P
    V <- V + C[2L * k + 1L] * P
  }
  U <- A %*% U
  solve(V - U, V + U)
}

.eu_pade_degree <- function(nA) {
  if (nA > 0.95) .EU_PADE$c9 else if (nA > 0.25) .EU_PADE$c7 else if (nA > 0.015) .EU_PADE$c5 else .EU_PADE$c3
}

# ExpMethodHigham2005Base: the reduced matrices of the Krylov methods, the
# augmented matrices of the phi functions
.eu_exp_hb <- function(A) {
  A <- as.matrix(A)
  if (!length(A)) return(A)
  bal <- .eu_gebal(A)
  A <- bal$A
  nA <- .eu_opnorm1(A)
  if (nA <= 2.1) {
    X <- .eu_pade_loop(A, .eu_pade_degree(nA))
  } else {
    s <- log2(nA / 5.4)
    si <- 0
    if (s > 0) {
      si <- ceiling(s)
      A <- A / 2^si
    }
    X <- .eu_pade_loop(A, .EU_PADE$c13)
    for (i in seq_len(si)) X <- X %*% X
  }
  .eu_gebak(X, bal)
}

# exp(A) of Julia's LinearAlgebra (exp! for dense matrices)
.eu_exp_base <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  if (!n) return(A)
  if (all(A[row(A) != col(A)] == 0)) return(diag(exp(diag(A)), n))
  if (all(A == t(A))) {
    e <- eigen(A, symmetric = TRUE)
    Uv <- e$vectors %*% diag(exp(e$values / 2), n)
    return(tcrossprod(Uv))
  }
  bal <- .eu_gebal(A)
  A <- bal$A
  nA <- .eu_opnorm1(A)
  if (nA <= 2.1) {
    X <- .eu_pade_loop(A, .eu_pade_degree(nA))
  } else {
    s <- log2(nA / 5.4)
    si <- 0
    if (s > 0) {
      si <- ceiling(s)
      A <- A / 2^si
    }
    CC <- .EU_PADE$c13
    I <- diag(n)
    A2 <- A %*% A
    A4 <- A2 %*% A2
    A6 <- A2 %*% A4
    U <- A %*% (A6 %*% (CC[14] * A6 + CC[12] * A4 + CC[10] * A2) + (CC[8] * A6 + CC[6] * A4 + CC[4] * A2 + CC[2] * I))
    V <- A6 %*% (CC[13] * A6 + CC[11] * A4 + CC[9] * A2) + (CC[7] * A6 + CC[5] * A4 + CC[3] * A2 + CC[1] * I)
    X <- solve(V - U, V + U)
    for (i in seq_len(si)) X <- X %*% X
  }
  .eu_gebak(X, bal)
}

# ---------------------------------------------------------------------------
# phi functions

# phi(z, k) for a number: the first row of the exponential of the
# (k + 1) x (k + 1) matrix with z in the corner and ones on the superdiagonal
.eu_phi_scalar <- function(z, k) {
  M <- matrix(0, k + 1L, k + 1L)
  M[1L, 1L] <- z
  for (i in seq_len(k)) M[i, i + 1L] <- 1
  .eu_exp_hb(M)[1L, ]
}

# phiv_dense(A, v, k): [phi_0(A) v, ..., phi_k(A) v] through the exponential
# of the (m + k) x (m + k) augmented matrix
.eu_phiv_dense <- function(A, v, k) {
  m <- length(v)
  C <- matrix(0, m + k, m + k)
  C[seq_len(m), seq_len(m)] <- A
  C[seq_len(m), m + 1L] <- v
  if (k > 1L) for (i in (m + 1L):(m + k - 1L)) C[i, i + 1L] <- 1
  P <- .eu_exp_hb(C)
  w <- matrix(0, m, k + 1L)
  w[, 1L] <- P[seq_len(m), seq_len(m), drop = FALSE] %*% v
  for (i in seq_len(k)) w[, i + 1L] <- P[seq_len(m), m + i]
  w
}

# Matrix phi functions phi_0(A), ..., phi_p(A): the scaling and recovering
# algorithm of Al-Mohy and Liu (phi_almohy.jl, _phi_almohy!)
.EU_PHI_THETA <- matrix(c(
  1.999463452408407e-5, 0x1.3bac6a3574299p-15, 7.366006416045163e-5, 0.00014973317297025854, 0.0003152443333771182, 0.0006855209983764435, 0.0015357294906993542,
  0.0038062018282832713, 0.006090206286125726, 0x1.4368e64f45ad8p-7, 0.016211831146383013, 0x1.ba1c1baea9a28p-6, 0.0454757968381803, 0.07749273259855331,
  0.03971636005661334, 0.05806968886880692, 0.08534220076759817, 0.12612151517169362, 0.18736524835661617, 0.27955116495245524, 0x1.ac43caa219960p-2,
  0.15442675548312682, 0x1.b3c69d3eb7bc3p-3, 0.29371996708854947, 0.40617647304246707, 0.5623843002320996, 0.7787883505754265, 1.0464245027100287,
  0.37980898016147974, 0.5014624976587007, 0.6621033824456818, 0.8739858777744354, 1.110828442901451, 1.3375805277521873, 1.5770955744229305,
  0.726177195703321, 0.9281910159646274, 1.1591052927815408, 1.4012982671152012, 0x1.a833aeef17afep+0, 1.924026602416338, 2.20029216725632,
  1.1898666361923196, 1.4469917284172122, 1.7187282676344413, 2.0024009817133357, 2.295729523138923, 2.5968019781751672, 2.9040336511997444,
  0x1.c2b573a2560d9p+0, 2.060907194742016, 2.371480030315257, 2.690078917436836, 3.014877598333601, 3.3443898815994957, 3.677415494646134,
  2.425818958547233, 2.7623053687512256, 3.105174730511504, 3.4527057639100476, 0x1.e6d9daf98e155p+1, 4.156538184296475, 4.510892686686127,
  3.173113456793749, 3.5393251225873685, 3.9086764034913664, 4.279911937822335, 4.652058062541427, 5.024366624424222, 5.396268333440519,
  3.991025201329815, 4.3813599805888455, 4.772204432206606, 5.162703837030561, 5.552220615863406, 5.940286915538899, 6.326566605567779,
  4.869485489784578, 5.279199870248916, 5.687344501175957, 6.09339280009268, 6.496977936247796, 6.8978550856394785, 7.295871912058635),
  nrow = 12L, byrow = TRUE)

.eu_phi_theta <- function(m, p) .EU_PHI_THETA[m, min(p, 7L)]

.eu_factf <- function(k) {
  r <- 1
  for (i in seq_len(max(k - 1L, 0L)) + 1L) r <- r * i
  r
}

.eu_phi_log_be_coeff <- function(m, p) {
  r <- 0
  for (j in seq_len(m)) r <- r - log(m + p + j)
  for (j in seq_len(m + p + 1L)) r <- r - log(m + j)
  r
}

.eu_phi_ell <- function(A, normT, m, p, phat) {
  if (normT == 0) return(0)
  t0 <- if (normT > 1) log2(normT) else 0
  scalefac <- 2^t0
  normTs <- normT / scalefac
  delta <- (p - 1) * (p - phat) / p + 1
  K <- 2L * m + p + 1L
  cc <- exp((.eu_phi_log_be_coeff(m, p) - delta * log(normTs)) / K)
  B <- cc * abs(A) / scalefac
  v <- rep(1, nrow(A))
  for (i in seq_len(K)) v <- as.numeric(crossprod(B, v))
  t <- log2(2 * max(v) / 2.220446049250313e-16) / (K - delta) + t0
  if (!is.finite(t)) return(0)
  max(ceiling(t), 0)
}

.eu_phi_select <- function(A, p) {
  phat <- if (.eu_phi_theta(12L, p) < 1) 0L else p
  r_max <- floor((1 + sqrt(1 + 4 * (2 * 12 + phat + 1))) / 2)
  eta <- numeric(r_max)
  cur <- A %*% A
  eta[1L] <- .eu_opnorm1(cur)^(1 / 2)
  for (j in seq_len(r_max - 1L) + 1L) {
    cur <- cur %*% A
    eta[j] <- .eu_opnorm1(cur)^(1 / (j + 1))
  }
  alpha <- pmax(eta[-r_max], eta[-1L])
  normT <- .eu_opnorm1(A)
  imax <- 7L
  tvals <- vapply(0:imax, function(i) {
    m_i <- (i + 3L)^2 %/% 8L
    .eu_phi_ell(A, normT, m_i, p, if (.eu_phi_theta(m_i, p) < 1) 0L else p)
  }, 0)
  C <- matrix(0, imax + 1L, r_max - 1L)
  for (i in 0:imax) {
    m_i <- (i + 3L)^2 %/% 8L
    th <- .eu_phi_theta(m_i, p)
    phat_i <- if (th < 1) 0L else p
    t <- tvals[i + 1L]
    for (r in seq_len(r_max - 1L) + 1L) {
      if (2 * m_i + phat_i + 1 >= r * (r - 1)) {
        a <- alpha[r - 1L]
        s0 <- if (a > 0 && is.finite(a)) max(ceiling(log2(a / th)), t) else t
        C[i + 1L, r - 1L] <- i + p + s0 * (p + 1)
      }
    }
  }
  minval <- min(C[C > 0])
  # the last minimum in column major order
  hit <- which(C == minval)
  i_star <- (hit[length(hit)] - 1L) %% (imax + 1L)
  m <- (i_star + 3L)^2 %/% 8L
  s <- round((minval - i_star - p) / (p + 1))
  tau <- floor(sqrt(2 * m))
  if (tau - 1 + 2 * (m %/% tau) - 2 * (m %% tau == 0) != i_star) tau <- ceiling(sqrt(2 * m))
  list(m = m, s = s, tau = tau)
}

.eu_phi_pade_coef <- function(m, p) {
  n1 <- 1
  for (i in (m + 1L):(2L * m + 1L)) n1 <- n1 * i
  d1 <- n1
  Nc <- Dc <- numeric(m + 1L)
  for (k in seq_len(p)) {
    if (k == p) {
      Am <- matrix(0, m + 1L, m + 1L)
      cp <- 1
      Am[1L, 1L] <- n1
      for (r in seq_len(m) + 1L) {
        cp <- cp / (k + (r - 1))
        Am[r, 1L] <- n1 * cp
      }
      for (cc in seq_len(m) + 1L) for (r in seq_len(m) + 1L) {
        I <- r - 1
        J <- cc - 1
        Am[r, cc] <- (-(m + 1 - J) * (k + 1 + I - J)) / ((2 * m + k + 1 - J) * J)
      }
      for (r in seq_len(m + 1L)) {
        s <- 0
        pr <- 1
        for (cc in seq_len(r)) {
          pr <- pr * Am[r, cc]
          s <- s + pr
        }
        Nc[r] <- s
      }
      Dc[1L] <- d1
      cpd <- 1
      for (r in seq_len(m) + 1L) {
        i <- r - 1
        cpd <- cpd * (-(m + 1 - i) / (i * (2 * m + k + 1 - i)))
        Dc[r] <- d1 * cpd
      }
    }
    n1 <- n1 * (2 * m + k + 1) / (k + 1)
    d1 <- d1 * (2 * m + k + 1)
  }
  list(N = Nc, D = Dc)
}

.eu_paterson_stockmeyer <- function(As, Nc, Dc, m, tau) {
  n <- nrow(As)
  Apow <- vector("list", tau + 1L)
  Apow[[1L]] <- diag(n)
  Apow[[2L]] <- As
  for (i in seq_len(tau - 1L) + 2L) Apow[[i]] <- Apow[[i - 1L]] %*% As
  Atau <- Apow[[tau + 1L]]
  addblock <- function(M, coef, i) {
    start <- i * tau + 1L
    stop <- min((i + 1L) * tau, m + 1L)
    for (l in seq_len(stop - start + 1L)) M <- M + coef[start + l - 1L] * Apow[[l]]
    M
  }
  nu <- m %/% tau
  if (m %% tau == 0L) {
    N <- addblock(Nc[m + 1L] * Atau, Nc, nu - 1L)
    D <- addblock(Dc[m + 1L] * Atau, Dc, nu - 1L)
    inext <- nu - 2L
  } else {
    N <- addblock(matrix(0, n, n), Nc, nu)
    D <- addblock(matrix(0, n, n), Dc, nu)
    inext <- nu - 1L
  }
  if (inext >= 0L) for (i in inext:0L) {
    N <- addblock(N %*% Atau, Nc, i)
    D <- addblock(D %*% Atau, Dc, i)
  }
  list(N = N, D = D)
}

.eu_phi_matrix <- function(A, p) {
  A <- as.matrix(A)
  n <- nrow(A)
  sel <- .eu_phi_select(A, p)
  m <- sel$m
  s <- sel$s
  As <- A / 2^s
  cf <- .eu_phi_pade_coef(m, p)
  ND <- .eu_paterson_stockmeyer(As, cf$N, cf$D, m, sel$tau)
  out <- vector("list", p + 1L)
  X <- tryCatch(solve(ND$D, ND$N), error = function(e) matrix(NaN, n, n))
  if (!all(is.finite(X))) return(rep(list(matrix(NaN, n, n)), p + 1L))
  out[[p + 1L]] <- X
  for (k in p:1L) {
    out[[k]] <- As %*% out[[k + 1L]]
    diag(out[[k]]) <- diag(out[[k]]) + 1 / .eu_factf(k - 1L)
  }
  for (it in seq_len(s)) {
    for (j in p:1L) {
      tmp <- out[[1L]] %*% out[[j + 1L]]
      for (k in seq_len(j)) tmp <- tmp + (1 / .eu_factf(j - k)) * out[[k + 1L]]
      out[[j + 1L]] <- tmp / 2^j
    }
    out[[1L]] <- out[[1L]] %*% out[[1L]]
  }
  out
}

# phi(A, k) for a matrix, a diagonal operator or a number: list of phi_0..phi_k
.eu_phi <- function(A, k) {
  if (is.numeric(A) && !is.matrix(A) && length(A) == 1L) return(as.list(.eu_phi_scalar(A, k)))
  if (inherits(A, "DiagonalOperator") || inherits(A, "eu_diagonal")) {
    d <- if (inherits(A, "eu_diagonal")) unclass(A) else A$diag
    vals <- vapply(d, .eu_phi_scalar, numeric(k + 1L), k = k)
    vals <- matrix(vals, nrow = k + 1L)
    return(lapply(seq_len(k + 1L), function(j) diag(vals[j, ], length(d))))
  }
  .eu_phi_matrix(.eu_matrix(A), k)
}

# ---------------------------------------------------------------------------
# Krylov subspaces (arnoldi.jl)

# KrylovSubspace: m, beta, V (n x (m + 1)), H ((m + 1) x m), wasbreakdown
.eu_arnoldi <- function(A, b, m, tol = 1e-7, ishermitian = .eu_ishermitian(A), iop = 0L) {
  n <- length(b)
  V <- matrix(0, n, m + 1L)
  H <- matrix(0, m + 1L, m)
  beta <- .eu_norm2(b)
  Ks <- list(m = m, beta = beta, V = V, H = H, wasbreakdown = FALSE)
  if (beta == 0) return(Ks)
  V[, 1L] <- b * (1 / beta)
  if (ishermitian) {
    # lanczos!
    for (j in seq_len(m)) {
      x <- V[, j]
      y <- .eu_mul(A, x)
      al <- sum(x * y)
      H[j, j] <- al
      y <- y - al * x
      if (j > 1L) y <- y - H[j, j - 1L] * V[, j - 1L]
      bt <- .eu_norm2(y)
      H[j + 1L, j] <- bt
      V[, j + 1L] <- y / bt
      if (tol > bt) {
        Ks$m <- j
        Ks$wasbreakdown <- TRUE
        break
      }
    }
    for (j in seq_len(m - 1L)) H[j, j + 1L] <- H[j + 1L, j]
  } else {
    if (iop == 0L) iop <- m
    for (j in seq_len(m)) {
      y <- .eu_mul(A, V[, j])
      for (i in max(1L, j - iop + 1L):j) {
        al <- sum(V[, i] * y)
        H[i, j] <- al
        y <- y - al * V[, i]
      }
      bt <- .eu_norm2(y)
      H[j + 1L, j] <- bt
      V[, j + 1L] <- y / bt
      if (bt < tol) {
        Ks$m <- j
        Ks$wasbreakdown <- TRUE
        break
      }
    }
  }
  Ks$V <- V
  Ks$H <- H
  Ks
}

.eu_getV <- function(Ks) Ks$V[, seq_len(Ks$m + 1L), drop = FALSE]
.eu_getH <- function(Ks) Ks$H[seq_len(Ks$m + 1L), seq_len(Ks$m), drop = FALSE]

# expv!(w, t, Ks)
.eu_expv_ks <- function(t, Ks) {
  m <- Ks$m
  n <- nrow(Ks$V)
  if (Ks$beta == 0) return(numeric(n))
  Hm <- Ks$H[seq_len(m), seq_len(m), drop = FALSE]
  Vm <- Ks$V[, seq_len(m), drop = FALSE]
  if (all(Hm == t(Hm))) {
    Tm <- matrix(0, m, m)
    diag(Tm) <- diag(Hm)
    if (m > 1L) for (j in seq_len(m - 1L)) Tm[j + 1L, j] <- Tm[j, j + 1L] <- Hm[j + 1L, j]
    e <- eigen(Tm, symmetric = TRUE)
    ord <- rev(seq_len(m))
    vec <- e$vectors[, ord, drop = FALSE]
    expHe <- vec %*% (exp(t * e$values[ord]) * vec[1L, ])
    return(Ks$beta * as.numeric(Vm %*% expHe))
  }
  E <- .eu_exp_hb(t * Hm)
  Ks$beta * as.numeric(Vm %*% E[, 1L])
}

# _phiv!(w, t, Ks, k): list(w = n x (k + 1), err)
.eu_phiv_ks <- function(t, Ks, k, correct = FALSE) {
  m <- Ks$m
  n <- nrow(Ks$V)
  if (Ks$beta == 0) return(list(w = matrix(0, n, k + 1L), err = 0))
  H <- .eu_getH(Ks)
  Hm <- t * Ks$H[seq_len(m), seq_len(m), drop = FALSE]
  e <- numeric(m)
  e[1L] <- 1
  C2 <- .eu_phiv_dense(Hm, e, k)
  w <- Ks$beta * (Ks$V[, seq_len(m), drop = FALSE] %*% C2)
  hend <- H[m + 1L, m]
  if (correct) {
    betah <- Ks$beta * hend * t
    vlast <- Ks$V[, m + 1L]
    for (i in seq_len(k)) w[, i] <- w[, i] + (betah * C2[m, i + 1L]) * vlast
  }
  list(w = w, err = abs(Ks$beta * hend * t * C2[m, k + 1L]))
}

# expv(t, A, b; m, iop): the Krylov approximation of exp(t A) b
.eu_expv <- function(t, A, b, m = min(30L, length(b)), iop = 0L, tol = 1e-7,
                     ishermitian = .eu_ishermitian(A)) {
  .eu_expv_ks(t, .eu_arnoldi(A, b, m = m, tol = tol, ishermitian = ishermitian, iop = iop))
}

# phiv(t, A, b, k; m, iop): [phi_0(tA) b, ..., phi_k(tA) b]
.eu_phiv <- function(t, A, b, k, m = min(30L, length(b)), iop = 0L, tol = 1e-7,
                     ishermitian = .eu_ishermitian(A), correct = FALSE) {
  .eu_phiv_ks(t, .eu_arnoldi(A, b, m = m, tol = tol, ishermitian = ishermitian, iop = iop), k,
              correct = correct)$w
}

# ---------------------------------------------------------------------------
# Adaptive time stepping of the Krylov phi functions (krylov_phiv_adaptive.jl,
# Niesen and Wright, Algorithm 919)

.eu_phiv_timestep_flops <- function(m, tau, n, p, nnz, iop, Hnorm, maxtau) {
  flops_W <- 2 * (p - 1) * (nnz + n)
  flops_u <- (2 * p + 1) * n
  if (iop == 0) iop <- m
  flops_matvec <- 2 * m * nnz
  flops_vecvec <- 0
  for (i in seq_len(m)) flops_vecvec <- flops_vecvec + 3 * min(i, iop)
  MH <- 44 / 3 + 2 * ceiling(max(0, log2(Hnorm / 5.37)))
  flops_phiv <- round(MH * (m + p)^3)
  (flops_W + flops_u + flops_matvec + flops_vecvec + flops_phiv) * ceiling(maxtau / tau)
}

.eu_phiv_timestep_adapt <- function(m, tau, epsilon, m_old, tau_old, epsilon_old, q, kappa,
                                    gamma, omega, maxtau, n, p, nnz, iop, Hnorm) {
  if (tau_old > tau) q <- log(tau / tau_old) / log(epsilon / epsilon_old) - 1
  tau_new <- tau * (gamma / omega)^(1 / (q + 1))
  tau_new <- min(max(tau_new, tau / 5), 2 * tau, maxtau)
  if (m_old < m) kappa <- (epsilon / epsilon_old)^(1 / (m_old - m))
  m_new <- m + ceiling(log(omega / gamma) / log(kappa))
  m_new <- min(max(m_new, (3 * m) %/% 4, 1), ceiling(4 * m / 3))
  cost_tau <- .eu_phiv_timestep_flops(m, tau_new, n, p, nnz, iop, Hnorm, maxtau)
  cost_m <- .eu_phiv_timestep_flops(m_new, tau, n, p, nnz, iop, Hnorm, maxtau)
  if (cost_tau < cost_m) m_new <- m else tau_new <- tau
  list(m = as.integer(m_new), tau = tau_new, q = q, kappa = kappa)
}

# phiv_timestep(ts, A, B): u(t) = phi_0(t A) b_0 + t phi_1(t A) b_1 + ... +
# t^p phi_p(t A) b_p at the times ts (a column per time; a vector when ts
# is a single time)
.eu_phiv_timestep <- function(ts, A, B, tau = 0, m = min(10L, NROW(B)), tol = 1e-7, opnorm = NULL,
                              iop = 0L, correct = FALSE, adaptive = FALSE, delta = 1.2,
                              ishermitian = .eu_ishermitian(A), gamma = 0.8, nnz = 0L) {
  B <- as.matrix(B)
  single <- length(ts) == 1L
  arnoldi_scale <- is.null(opnorm)
  abstol <- NULL
  m <- as.integer(m)
  iop <- as.integer(iop)
  if (!arnoldi_scale) {
    opn <- if (is.numeric(opnorm)) opnorm else opnorm(A)
    abstol <- tol * opn
    if (tau == 0) {
      b0norm <- max(abs(B[, 1L]))
      tau <- 10 / opn * (abstol * ((m + 1) / exp(1))^(m + 1) * sqrt(2 * pi * (m + 1)) /
                           (4 * opn * b0norm))^(1 / m)
    }
  }
  n <- nrow(B)
  ts <- sort(ts)
  tend <- ts[length(ts)]
  seed_arnoldi_tau <- arnoldi_scale && tau == 0
  if (seed_arnoldi_tau) tau <- tend
  p <- ncol(B) - 1L
  U <- matrix(0, n, length(ts))
  u <- B[, 1L]
  coeffs <- rep(1, max(p, 1L))
  if (adaptive) {
    if (ishermitian) iop <- 2L
    if (nnz == 0L) nnz <- sum(.eu_matrix(A) != 0)
  }
  W <- matrix(0, n, p + 1L)
  t <- 0
  snapshot <- 1L
  while (t < tend) {
    if (t + tau > tend) tau <- tend - t
    W[, 1L] <- u
    for (l in seq_len(max(p - 1L, 0L))) coeffs[l + 1L] <- coeffs[l] * t / l
    for (j in seq_len(p)) {
      W[, j + 1L] <- .eu_mul(A, W[, j])
      for (l in 0:(p - j)) W[, j + 1L] <- W[, j + 1L] + coeffs[l + 1L] * B[, j + l + 1L]
    }
    Ks <- .eu_arnoldi(A, W[, p + 1L], m = m, tol = tol, iop = iop)
    if (is.null(abstol)) {
      opn <- .eu_opnorm1(.eu_getH(Ks))
      abstol <- tol * opn
      if (seed_arnoldi_tau) {
        b0norm <- max(abs(B[, 1L]))
        tau <- min(tend - t, gamma * 10 / opn * (abstol * ((m + 1) / exp(1))^(m + 1) * sqrt(2 * pi * (m + 1)) /
                                                   (4 * opn * b0norm))^(1 / m))
      }
    }
    breakdown <- Ks$wasbreakdown || Ks$beta == 0
    if (breakdown) tau <- tend - t
    ph <- .eu_phiv_ks(tau, Ks, p + 1L, correct)
    P <- ph$w
    epsilon <- ph$err
    if (adaptive && !breakdown) {
      omega <- (tend / tau) * (epsilon / abstol)
      epsilon_old <- epsilon
      m_old <- m
      tau_old <- tau
      q <- m / 4
      kappa <- 2
      maxtau <- tend - t
      while (omega > delta) {
        ad <- .eu_phiv_timestep_adapt(m, tau, epsilon, m_old, tau_old, epsilon_old, q, kappa, gamma, omega,
                                      maxtau, n, p, nnz, iop, .eu_opnorm1(.eu_getH(Ks)))
        q <- ad$q
        kappa <- ad$kappa
        m_old <- m
        m <- ad$m
        tau_old <- tau
        tau <- ad$tau
        Ks <- .eu_arnoldi(A, W[, p + 1L], m = m, tol = tol, iop = iop)
        ph <- .eu_phiv_ks(tau, Ks, p + 1L, correct)
        P <- ph$w
        epsilon_old <- epsilon
        epsilon <- ph$err
        omega <- (tend / tau) * (epsilon / abstol)
      }
    }
    u <- tau^p * P[, p + 1L]
    for (l in seq_len(max(p - 1L, 0L))) coeffs[l + 1L] <- coeffs[l] * tau / l
    for (j in seq_len(p) - 1L) u <- u + coeffs[j + 1L] * W[, j + 1L]
    while (snapshot <= length(ts) && t + tau >= ts[snapshot]) {
      tau_s <- ts[snapshot] - t
      Ps <- .eu_phiv_ks(tau_s, Ks, p + 1L, correct)$w
      us <- tau_s^p * Ps[, p + 1L]
      for (l in seq_len(max(p - 1L, 0L))) coeffs[l + 1L] <- coeffs[l] * tau_s / l
      for (j in seq_len(p) - 1L) us <- us + coeffs[j + 1L] * W[, j + 1L]
      U[, snapshot] <- us
      snapshot <- snapshot + 1L
    }
    t <- t + tau
  }
  if (single) U[, 1L] else U
}

# expv_timestep(ts, A, b) = phiv_timestep(ts, A, b) with p = 0
.eu_expv_timestep <- function(ts, A, b, ...) .eu_phiv_timestep(ts, A, matrix(b, ncol = 1L), ...)
