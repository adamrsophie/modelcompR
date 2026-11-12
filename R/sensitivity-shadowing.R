# Shadowing sensitivities of long time averages of chaotic systems, after
# SciMLSensitivity.jl v7.119 src/lss.jl (ForwardLSSProblem,
# AdjointLSSProblem, LSSSchur, wB!, wE!, B!, E!, SchurLU, b!,
# wBcorrect!, h!, shadow_forward, shadow_adjoint), src/nilss.jl
# (NILSSProblem, forward_sense, perp!, renormalize!, compute_Cinv!,
# compute_d!, compute_B!, compute_a!, compute_v!, compute_xi,
# shadow_forward) and src/nilsas.jl (NILSASProblem, terminate_conditions,
# adjoint_sense, renormalize!, nilsas_min, shadow_adjoint).
#
# Differences from Julia: the LSS Schur complement B W B' + E W E' is block
# tridiagonal and is solved by block elimination instead of a dense LU; the
# Jacobians at the midpoints use the midpoint times; the NILSAS adjoint
# reads the forward state from the dense forward solution (Julia integrates
# it backward, BacksolveAdjoint, with checkpoints at the saved points).

.lss_check_g <- function(sensealg, g) {
  if (inherits(sensealg$LSSregularizer, "TimeDilation") && is.null(g))
    stop("Time dilation needs explicit knowledge of g. Either pass `g` as an argument to `ForwardLSS(g = g)` ",
         "or `AdjointLSS(g = g)` or use ForwardLSS/AdjointLSS with windowing.", call. = FALSE)
}

# Gradients of the objective at a saved point: list(u, p), p NULL when the
# objective's parameter gradient is not used.
.lss_dg <- function(g, dgdu, dgdp, np) {
  function(u, p, t, i) {
    if (is.null(dgdu)) {
      return(list(u = .sens_fd_grad(function(z) g(z, p, t), u),
                  p = if (np) .sens_fd_grad(function(q) g(u, .sens_pvec(p, q), t), as.numeric(p)) else numeric()))
    }
    list(u = as.numeric(dgdu(u, p, t, i)), p = if (!is.null(dgdp)) as.numeric(dgdp(u, p, t, i)))
  }
}

.lss_costs <- function(sol, t, dgdu_discrete, dgdp_discrete, dgdu_continuous, dgdp_continuous) {
  if (!is.null(t)) {
    if (length(t) != length(sol$t) || any(abs(t - sol$t) > 1e-12 * max(1, abs(t))))
      stop("the discrete cost times must be the saved times of the solution", call. = FALSE)
    return(list(dgdu = dgdu_discrete, dgdp = dgdp_discrete))
  }
  list(dgdu = if (!is.null(dgdu_continuous)) function(u, p, tt, i) dgdu_continuous(u, p, tt),
       dgdp = if (!is.null(dgdp_continuous)) function(u, p, tt, i) dgdp_continuous(u, p, tt))
}

#' Least squares shadowing problems
#'
#' `ForwardLSSProblem(sol, ForwardLSS(g = g))` and `AdjointLSSProblem(sol,
#' AdjointLSS(g = g))` set up least squares shadowing on the saved points of
#' a solution `sol` on the attractor; `shadow_forward()` and
#' `shadow_adjoint()` return the derivative of the long time average of
#' `g(u, p, t)` with respect to the parameters. `NILSSProblem(prob,
#' NILSS(nseg, nstep, g = g))` and `NILSASProblem(sol, NILSAS(nseg, nstep,
#' M, g = g), alg)` are the non intrusive forward and adjoint versions; they
#' integrate tangent (or adjoint) solutions over `nseg` segments with
#' `nstep` saved points each, so the cost scales with the number of unstable
#' directions instead of the number of states.
#'
#' The objective's gradients are taken by finite differences of `g` unless
#' given: `dgdu_continuous(u, p, t)` and `dgdp_continuous(u, p, t)`, or
#' with `t = sol$t` the discrete forms `dgdu_discrete(u, p, t, i)` and
#' `dgdp_discrete(u, p, t, i)`.
#'
#' @param sol A solution on the attractor (`saveat` gives the time grid of
#'   LSS); for `NILSASProblem` it must be dense.
#' @param prob An `ODEProblem` whose initial value is on the attractor.
#' @param sensealg The shadowing algorithm, see [shadowing_algorithms].
#' @param alg ODE algorithm for the segment solves of NILSS and NILSAS.
#' @param t,dgdu_discrete,dgdp_discrete,dgdu_continuous,dgdp_continuous
#'   Gradients of the objective (see Details).
#' @param g The objective; by default `sensealg$g`.
#' @param ... Options for the segment solves ([solve()]).
#' @return A problem object; `shadow_forward()` and `shadow_adjoint()`
#'   return the gradient with respect to the parameters.
#' @examples
#' lorenz <- function(u, p, t) c(10 * (u[2] - u[1]), u[1] * (p[1] - u[3]) - u[2],
#'                               u[1] * u[2] - 8 / 3 * u[3])
#' prob <- ODEProblem(lorenz, c(-2.26, -2.47, 18.94), c(0, 5), 28)
#' sol <- solve(prob, Tsit5(), saveat = 0.01, reltol = 1e-8, abstol = 1e-8)
#' g <- function(u, p, t) u[3]
#' shadow_forward(ForwardLSSProblem(sol, ForwardLSS(g = g)))
#' @export
ForwardLSSProblem <- function(sol, sensealg = ForwardLSS(), t = NULL, dgdu_discrete = NULL, dgdp_discrete = NULL,
                              dgdu_continuous = NULL, dgdp_continuous = NULL, g = sensealg$g) {
  .lss_problem(sol, sensealg, "ForwardLSSProblem", t, dgdu_discrete, dgdp_discrete, dgdu_continuous,
               dgdp_continuous, g)
}

#' @rdname ForwardLSSProblem
#' @export
AdjointLSSProblem <- function(sol, sensealg = AdjointLSS(), t = NULL, dgdu_discrete = NULL, dgdp_discrete = NULL,
                              dgdu_continuous = NULL, dgdp_continuous = NULL, g = sensealg$g) {
  .lss_problem(sol, sensealg, "AdjointLSSProblem", t, dgdu_discrete, dgdp_discrete, dgdu_continuous,
               dgdp_continuous, g)
}

.lss_problem <- function(sol, sensealg, cls, t, dgdu_discrete, dgdp_discrete, dgdu_continuous, dgdp_continuous, g) {
  env <- .sol_env(sol)
  prob <- env$prob
  .lss_check_g(sensealg, g)
  if (!length(prob$p)) stop("You must have parameters to use parameter sensitivity calculations!", call. = FALSE)
  fns <- .sens_fns(prob, sensealg)
  if (!is.null(fns$mass)) stop("shadowing methods are only compatible with ODEs, not with DAEs", call. = FALSE)
  costs <- .lss_costs(sol, t, dgdu_discrete, dgdp_discrete, dgdu_continuous, dgdp_continuous)
  if (is.null(g) && is.null(costs$dgdu)) stop("give the objective `g` or its gradient", call. = FALSE)
  ts <- env$t
  U <- env$U
  n <- ncol(U); Nt <- length(ts); Ndt <- Nt - 1L
  p <- env$p_final %||% prob$p
  dt <- diff(ts)
  umid <- (U[-1L, , drop = FALSE] + U[-Nt, , drop = FALSE]) / 2
  dudt <- (U[-1L, , drop = FALSE] - U[-Nt, , drop = FALSE]) / dt
  tmid <- (ts[-1L] + ts[-Nt]) / 2
  Dt <- ts[Nt] - ts[1L]
  # wB!: inverse weights of v, wE!: of the time dilation
  wB <- matrix(1, Nt, n)
  wB[1L, ] <- 1 / dt[1L]
  wB[Nt, ] <- 1 / dt[Ndt]
  if (Nt > 2L) wB[2:(Nt - 1L), ] <- 1 / (dt[-1L] + dt[-Ndt])
  wB <- wB * 2 * Dt
  reg <- sensealg$LSSregularizer
  td <- inherits(reg, "TimeDilation")
  wE <- if (td) Dt / (reg$alpha^2 * dt) else NULL
  # B!: row block i holds B1_i at column block i and B0_i at block i + 1
  Js <- lapply(seq_len(Ndt), function(i) fns$jac(umid[i, ], p, tmid[i]))
  I <- diag(1, n)
  B1 <- lapply(seq_len(Ndt), function(i) -I / dt[i] - Js[[i]] / 2)
  B0 <- lapply(seq_len(Ndt), function(i) I / dt[i] - Js[[i]] / 2)
  # the Schur complement B W B' + E W E', block tridiagonal
  Sd <- lapply(seq_len(Ndt), function(i) {
    S <- B1[[i]] %*% (wB[i, ] * t(B1[[i]])) + B0[[i]] %*% (wB[i + 1L, ] * t(B0[[i]]))
    if (td) S <- S + wE[i] * tcrossprod(dudt[i, ])
    S
  })
  So <- if (Ndt > 1L) lapply(seq_len(Ndt - 1L), function(i) B0[[i]] %*% (wB[i + 1L, ] * t(B1[[i + 1L]]))) else list()
  out <- list(sensealg = sensealg, sol = sol, prob = prob, fns = fns, ts = ts, U = U, p = p, n = n, Nt = Nt,
              Ndt = Ndt, dt = dt, umid = umid, dudt = dudt, tmid = tmid, Dt = Dt, wB = wB, wE = wE,
              B1 = B1, B0 = B0, schur = .lss_block_factor(Sd, So), g = g,
              dg = .lss_dg(g, costs$dgdu, costs$dgdp, length(p)), discrete = !is.null(t))
  structure(out, class = c(cls, "LSSProblem"))
}

# Block LDL' style elimination of a symmetric block tridiagonal matrix with
# diagonal blocks Sd and upper blocks So (Sd[[i]], So[[i]] = S[i, i + 1]).
.lss_block_factor <- function(Sd, So) {
  N <- length(Sd)
  Dinv <- vector("list", N)
  Dinv[[1L]] <- solve(Sd[[1L]])
  for (i in seq_len(N - 1L) + 1L) {
    L <- t(So[[i - 1L]]) %*% Dinv[[i - 1L]]
    Dinv[[i]] <- solve(Sd[[i]] - L %*% So[[i - 1L]])
  }
  list(Dinv = Dinv, So = So, N = N)
}

# Solve S X = R for R with n N rows.
.lss_block_solve <- function(F, R) {
  R <- as.matrix(R)
  N <- F$N; n <- nrow(F$Dinv[[1L]])
  blk <- function(i) (i - 1L) * n + seq_len(n)
  Y <- R
  for (i in seq_len(N - 1L) + 1L)
    Y[blk(i), ] <- R[blk(i), , drop = FALSE] - t(F$So[[i - 1L]]) %*% (F$Dinv[[i - 1L]] %*% Y[blk(i - 1L), , drop = FALSE])
  X <- Y
  X[blk(N), ] <- F$Dinv[[N]] %*% Y[blk(N), , drop = FALSE]
  for (i in rev(seq_len(N - 1L)))
    X[blk(i), ] <- F$Dinv[[i]] %*% (Y[blk(i), , drop = FALSE] - F$So[[i]] %*% X[blk(i + 1L), , drop = FALSE])
  X
}

# B' w (one block per saved point) and B v.
.lss_Bt <- function(P, w) {
  n <- P$n
  W <- matrix(w, n)
  out <- matrix(0, n, P$Nt)
  for (i in seq_len(P$Ndt)) {
    out[, i] <- out[, i] + crossprod(P$B1[[i]], W[, i])
    out[, i + 1L] <- out[, i + 1L] + crossprod(P$B0[[i]], W[, i])
  }
  out
}

.lss_B <- function(P, v) {
  V <- matrix(v, P$n)
  vapply(seq_len(P$Ndt), function(i) as.numeric(P$B1[[i]] %*% V[, i] + P$B0[[i]] %*% V[, i + 1L]),
         numeric(P$n))
}

.lss_window <- function(P) {
  w <- 1 - cos((P$ts - P$ts[1L]) * 2 * pi / P$Dt)
  if (inherits(P$sensealg$LSSregularizer, "Cos2Windowing")) w <- w^2
  w / sum(w)
}

#' @rdname ForwardLSSProblem
#' @param lssprob A problem from `ForwardLSSProblem()`, `NILSSProblem()`,
#'   `AdjointLSSProblem()` or `NILSASProblem()`.
#' @export
shadow_forward <- function(lssprob, alg = NULL, ...) {
  if (inherits(lssprob, "NILSSProblem")) return(.nilss_forward(lssprob, alg, list(...)))
  if (!inherits(lssprob, "ForwardLSSProblem")) stop("shadow_forward() needs a ForwardLSSProblem or NILSSProblem",
                                                    call. = FALSE)
  P <- lssprob
  n <- P$n; np <- length(P$p)
  # b!: the parameter Jacobians at the midpoints
  b <- do.call(rbind, lapply(seq_len(P$Ndt), function(i) P$fns$pjac(P$umid[i, ], P$p, P$tmid[i])))
  W <- .lss_block_solve(P$schur, b)
  dgs <- lapply(seq_len(P$Nt), function(j) P$dg(P$U[j, ], P$p, P$ts[j], j))
  res <- numeric(np)
  reg <- P$sensealg$LSSregularizer
  if (inherits(reg, "TimeDilation")) {
    n0 <- which(P$ts >= P$ts[1L] + reg$t0skip)[1L]
    n1 <- which(P$ts >= P$ts[P$Nt] - reg$t1skip)[1L]
    gmid <- vapply(n0:(n1 - 1L), function(j) P$g(P$umid[j, ], P$p, P$tmid[j]), 0)
    for (i in seq_len(np)) {
      v <- t(P$wB) * .lss_Bt(P, W[, i])
      eta <- P$wE * vapply(seq_len(P$Ndt), function(j) sum(P$dudt[j, ] * W[(j - 1L) * n + seq_len(n), i]), 0)
      acc <- 0
      for (j in n0:n1) acc <- acc + sum(dgs[[j]]$u * v[, j]) + if (length(dgs[[j]]$p)) dgs[[j]]$p[i] else 0
      r <- acc / (n1 - n0 + 1)
      er <- eta[n0:(n1 - 1L)]
      r <- r - sum(er * gmid) / (n1 - n0) + sum(er) * sum(gmid) / (n1 - n0)^2
      res[i] <- r
    }
  } else {
    win <- .lss_window(P)
    for (i in seq_len(np)) {
      v <- t(P$wB) * .lss_Bt(P, W[, i])
      res[i] <- sum(vapply(seq_len(P$Nt), function(j)
        (sum(dgs[[j]]$u * v[, j]) + if (length(dgs[[j]]$p)) dgs[[j]]$p[i] else 0) * win[j], 0))
    }
  }
  stats::setNames(res, names(P$p))
}

#' @rdname ForwardLSSProblem
#' @export
shadow_adjoint <- function(lssprob, alg = NULL, ...) {
  if (inherits(lssprob, "NILSASProblem")) return(.nilsas_adjoint(lssprob, alg %||% lssprob$alg, list(...)))
  if (!inherits(lssprob, "AdjointLSSProblem")) stop("shadow_adjoint() needs an AdjointLSSProblem or NILSASProblem",
                                                    call. = FALSE)
  P <- lssprob
  reg <- P$sensealg$LSSregularizer
  if (!inherits(reg, "TimeDilation")) stop("AdjointLSS supports the TimeDilation regularizer", call. = FALSE)
  n <- P$n; np <- length(P$p)
  dgs <- lapply(seq_len(P$Nt), function(j) P$dg(P$U[j, ], P$p, P$ts[j], j))
  # wBcorrect!, h!
  wb <- t(P$wB * do.call(rbind, lapply(dgs, `[[`, "u")) / P$Nt)
  h <- vapply(seq_len(P$Ndt), function(j) P$g(P$umid[j, ], P$p, P$tmid[j]), 0)
  h <- P$wE * (-(h - mean(h)) / P$Ndt)
  bvec <- as.numeric(.lss_B(P, as.numeric(wb)) + t(P$dudt * h))
  wa <- matrix(.lss_block_solve(P$schur, bvec), n)
  n0 <- which(P$ts >= P$ts[1L] + reg$t0skip)[1L]
  n1 <- which(P$ts >= P$ts[P$Nt] - reg$t1skip)[1L]
  res <- numeric(np)
  if (length(dgs[[1L]]$p)) {
    gp <- vapply(n0:(n1 - 1L), function(j) {
      d <- P$dg(P$umid[j, ], P$p, P$tmid[j], j)$p
      if (length(d)) d else numeric(np)
    }, numeric(np))
    res <- rowMeans(matrix(gp, np))
  }
  for (j in n0:(n1 - 1L)) res <- res + as.numeric(crossprod(P$fns$pjac(P$umid[j, ], P$p, P$tmid[j]), wa[, j]))
  stats::setNames(res, names(P$p))
}

# ---------------------------------------------------------------------------
# NILSS

#' @rdname ForwardLSSProblem
#' @export
NILSSProblem <- function(prob, sensealg, t = NULL, dgdu_discrete = NULL, dgdp_discrete = NULL,
                         dgdu_continuous = NULL, dgdp_continuous = NULL, g = sensealg$g) {
  if (is.null(g)) stop("To use NILSS, g must be passed as an argument to `NILSS(g = g)`.", call. = FALSE)
  n <- length(prob$u0); np <- length(prob$p)
  if (!np) stop("You must have parameters to use parameter sensitivity calculations!", call. = FALSE)
  nus <- sensealg$nus %||% (n - 1L)
  if (nus >= n) stop("`nus` must be smaller than `numindvar`.", call. = FALSE)
  T_seg <- (prob$tspan[2L] - prob$tspan[1L]) / sensealg$nseg
  dtsave <- T_seg / (sensealg$nstep - 1L)
  if (!is.null(t)) {
    jevery <- round(stats::median(diff(t)) / dtsave)
    if (jevery < 1L || abs(jevery * dtsave - stats::median(diff(t))) > 1e-8 * dtsave)
      stop("the discrete cost times must be multiples of the NILSS save step", call. = FALSE)
  }
  fns <- .sens_fns(prob, sensealg)
  if (!is.null(fns$mass)) stop("shadowing methods are only compatible with ODEs, not with DAEs", call. = FALSE)
  w0 <- .with_seed(sensealg$seed, matrix(stats::runif(n * nus), n, nus))
  w0 <- sweep(w0, 2L, sqrt(colSums(w0^2)), "/")
  structure(list(sensealg = sensealg, prob = prob, fns = fns, nus = nus, T_seg = T_seg, dtsave = dtsave,
                 w0 = w0, g = g, t = t,
                 dgu = if (!is.null(t)) dgdu_discrete else if (!is.null(dgdu_continuous))
                   function(u, p, tt, i) dgdu_continuous(u, p, tt),
                 dgp = if (is.null(t) && !is.null(dgdp_continuous))
                   function(u, p, tt, i) dgdp_continuous(u, p, tt)),
            class = "NILSSProblem")
}

.nilss_forward <- function(P, alg, dots) {
  sa <- P$sensealg
  prob <- P$prob; fns <- P$fns
  n <- length(prob$u0); np <- length(prob$p); nus <- P$nus
  nseg <- sa$nseg; nstep <- sa$nstep
  p <- prob$p
  alg <- alg %||% Tsit5()
  # the tangent problem: y, nus homogeneous solutions, np inhomogeneous ones
  m <- nus + np
  rhs <- function(z, pp, tt) {
    y <- z[seq_len(n)]
    S <- matrix(z[n + seq_len(n * m)], n, m)
    dS <- fns$jac(y, pp, tt) %*% S
    dS[, nus + seq_len(np)] <- dS[, nus + seq_len(np)] + fns$pjac(y, pp, tt)
    c(fns$f(y, pp, tt), as.numeric(dS))
  }
  y0 <- as.numeric(prob$u0)
  W0 <- P$w0
  V0 <- matrix(0, n, np)
  ys <- vector("list", nseg); ws <- ys; vs <- ys
  R <- array(0, c(np, nseg - 1L, nus, nus)); bb <- array(0, c(np, (nseg - 1L) * nus))
  dudt <- ys; gsave <- matrix(0, nstep, nseg); dgdu <- ys
  wperp <- ys; vperp <- ys
  perp <- function(X, d) X - d %*% (crossprod(d, X) / sum(d * d))
  kcost <- 1L
  for (iseg in seq_len(nseg)) {
    t1 <- prob$tspan[1L] + (iseg - 1L) * P$T_seg
    t2 <- t1 + P$T_seg
    sp <- ODEProblem(rhs, c(y0, as.numeric(W0), as.numeric(V0)), c(t1, t2), p)
    sv <- do.call(solve, c(list(sp, alg), dots, list(saveat = seq(t1, t2, length.out = nstep))))
    Z <- as.matrix(sv)
    if (nrow(Z) != nstep) stop("a NILSS segment solve failed (retcode ", sv$retcode, ")", call. = FALSE)
    ys[[iseg]] <- Z[, seq_len(n), drop = FALSE]
    ws[[iseg]] <- lapply(seq_len(nstep), function(j) matrix(Z[j, n + seq_len(n * nus)], n, nus))
    vs[[iseg]] <- lapply(seq_len(nstep), function(j) matrix(Z[j, n + n * nus + seq_len(n * np)], n, np))
    dudt[[iseg]] <- t(vapply(seq_len(nstep), function(j) fns$f(ys[[iseg]][j, ], p, sv$t[j]), numeric(n)))
    gsave[, iseg] <- vapply(seq_len(nstep), function(j) P$g(ys[[iseg]][j, ], p, sv$t[j]), 0)
    dgdu[[iseg]] <- t(vapply(seq_len(nstep), function(j) {
      if (!is.null(P$t)) {
        jev <- round(stats::median(diff(P$t)) / P$dtsave)
        if ((j - 1L) %% jev != 0L) return(numeric(n))
        i <- kcost + (j - 1L) %/% jev
        return(as.numeric(P$dgu(ys[[iseg]][j, ], p, sv$t[j], i)))
      }
      if (!is.null(P$dgu)) as.numeric(P$dgu(ys[[iseg]][j, ], p, sv$t[j], j)) else
        .sens_fd_grad(function(z) P$g(z, p, sv$t[j]), ys[[iseg]][j, ])
    }, numeric(n)))
    if (!is.null(P$t)) kcost <- kcost + (nstep - 1L) %/% round(stats::median(diff(P$t)) / P$dtsave)
    wperp[[iseg]] <- lapply(seq_len(nstep), function(j) perp(ws[[iseg]][[j]], dudt[[iseg]][j, ]))
    vperp[[iseg]] <- lapply(seq_len(nstep), function(j) perp(vs[[iseg]][[j]], dudt[[iseg]][j, ]))
    if (iseg < nseg) {
      qrw <- qr(wperp[[iseg]][[nstep]])
      Q <- qr.Q(qrw)[, seq_len(nus), drop = FALSE]
      Rm <- qr.R(qrw)[seq_len(nus), seq_len(nus), drop = FALSE]
      # qr() pivots columns; undo it so that W_end = Q R
      Rm[, qrw$pivot] <- Rm
      vend <- vperp[[iseg]][[nstep]]
      for (i in seq_len(np)) {
        bi <- as.numeric(crossprod(Q, vend[, i]))
        bb[i, (iseg - 1L) * nus + seq_len(nus)] <- bi
        R[i, iseg, , ] <- Rm
        V0[, i] <- vend[, i] - Q %*% bi
      }
      W0 <- Q
      y0 <- ys[[iseg]][nstep, ]
    }
  }
  weight <- rep(1, nstep)
  weight[c(1L, nstep)] <- 0.5
  gavg <- sum(weight * gsave) / ((nstep - 1L) * nseg)
  # the time average of the direct parameter derivative of g (Julia's NILSS
  # leaves it out)
  res <- numeric(np)
  if (is.null(P$dgu) || !is.null(P$dgp)) {
    for (iseg in seq_len(nseg)) for (j in seq_len(nstep)) {
      tj <- prob$tspan[1L] + (iseg - 1L) * P$T_seg + (j - 1L) * P$dtsave
      y <- ys[[iseg]][j, ]
      gp <- if (!is.null(P$dgp)) as.numeric(P$dgp(y, p, tj, j)) else
        .sens_fd_grad(function(q) P$g(y, .sens_pvec(p, q), tj), as.numeric(p))
      res <- res + weight[j] * gp / ((nstep - 1L) * nseg)
    }
  }
  for (i in seq_len(np)) {
    Cinv <- matrix(0, nseg * nus, nseg * nus)
    d <- numeric(nseg * nus)
    for (iseg in seq_len(nseg)) {
      C <- Reduce(`+`, Map(function(w, wt) wt * crossprod(w), wperp[[iseg]], weight))
      blk <- (iseg - 1L) * nus + seq_len(nus)
      Cinv[blk, blk] <- solve(C)
      d[blk] <- Reduce(`+`, Map(function(w, v, wt) wt * as.numeric(crossprod(w, v[, i])),
                                wperp[[iseg]], vperp[[iseg]], weight))
    }
    B <- matrix(0, (nseg - 1L) * nus, nseg * nus)
    for (iseg in seq_len(nseg - 1L)) {
      blk <- (iseg - 1L) * nus + seq_len(nus)
      B[blk, blk] <- -R[i, iseg, , ]
      B[cbind(blk, iseg * nus + seq_len(nus))] <- 1
    }
    a <- if (nseg > 1L) {
      lbd <- solve(-B %*% Cinv %*% t(B), B %*% Cinv %*% d + bb[i, ])
      -Cinv %*% (t(B) %*% lbd + d)
    } else -Cinv %*% d
    for (iseg in seq_len(nseg)) {
      blk <- (iseg - 1L) * nus + seq_len(nus)
      v <- lapply(seq_len(nstep), function(j) vs[[iseg]][[j]][, i] + ws[[iseg]][[j]] %*% a[blk])
      vterm <- sum(vapply(seq_len(nstep), function(j) weight[j] * sum(v[[j]] * dgdu[[iseg]][j, ]), 0))
      de <- dudt[[iseg]][nstep, ]
      xi_end <- sum(v[[nstep]] * de) / sum(de * de)
      res[i] <- res[i] + vterm / ((nstep - 1L) * nseg) +
        xi_end * (gavg - gsave[nstep, iseg]) / (P$dtsave * (nstep - 1L) * nseg)
    }
  }
  stats::setNames(res, names(prob$p))
}

# ---------------------------------------------------------------------------
# NILSAS

#' @rdname ForwardLSSProblem
#' @export
NILSASProblem <- function(sol, sensealg, alg = NULL, t = NULL, dgdu_discrete = NULL, dgdp_discrete = NULL,
                          dgdu_continuous = NULL, dgdp_continuous = NULL, g = sensealg$g, ...) {
  if (is.null(g)) stop("To use NILSAS, g must be passed as an argument to `NILSAS(g = g)`.", call. = FALSE)
  if (is.null(sensealg$M))
    stop("Please provide an `M` with `M >= nus + 1`, where nus is the number of unstable covariant ",
         "Lyapunov vectors.", call. = FALSE)
  if (!is.null(t)) stop("`NILSAS` can currently not be used with discrete cost functions.", call. = FALSE)
  env <- .sol_env(sol)
  prob <- env$prob
  n <- length(prob$u0); np <- length(prob$p); M <- sensealg$M
  if (!np) stop("You must have parameters to use parameter sensitivity calculations!", call. = FALSE)
  fns <- .sens_fns(prob, sensealg)
  if (!is.null(fns$mass)) stop("shadowing methods are only compatible with ODEs, not with DAEs", call. = FALSE)
  tr <- .sens_trajectory(sol)
  tf <- prob$tspan[2L]
  # terminate_conditions: M - 1 random directions orthogonal to f, and f
  y <- env$U[nrow(env$U), ]
  fu <- fns$f(y, prob$p, tf)
  fu <- fu / sqrt(sum(fu^2))
  W <- if (M > 1L) {
    Wr <- .with_seed(sensealg$seed, matrix(stats::runif(n * (M - 1L)), n, M - 1L))
    Wr <- Wr - fu %*% crossprod(fu, Wr)
    cbind(qr.Q(qr(Wr))[, seq_len(M - 1L), drop = FALSE], fu)
  } else matrix(fu, n, 1L)
  gu <- dgdu_continuous %||% function(u, p, tt) .sens_fd_grad(function(z) g(z, p, tt), u)
  gp <- dgdp_continuous %||% if (is.null(dgdu_continuous)) function(u, p, tt)
    .sens_fd_grad(function(q) g(u, .sens_pvec(p, q), tt), as.numeric(p))
  structure(list(sensealg = sensealg, sol = sol, prob = prob, fns = fns, tr = tr, W0 = W, alg = alg,
                 T_seg = (tf - prob$tspan[1L]) / sensealg$nseg, gu = gu, gp = gp, dots = list(...)),
            class = "NILSASProblem")
}

.nilsas_adjoint <- function(P, alg, dots) {
  sa <- P$sensealg
  prob <- P$prob; fns <- P$fns; tr <- P$tr
  n <- length(prob$u0); np <- length(prob$p); M <- sa$M
  nseg <- sa$nseg; nstep <- sa$nstep
  p <- prob$p
  alg <- alg %||% Tsit5()
  dots <- c(P$dots, dots)
  # state: lambda_v, lambda_w (n x M), their parameter gradients, and the
  # quadratures C, dwv, dwf, dvf, dJs
  il <- seq_len(n * (M + 1L))
  ig <- max(il) + seq_len(np * (M + 1L))
  iC <- max(ig) + seq_len(M * M)
  iwv <- max(iC) + seq_len(M); iwf <- max(iwv) + seq_len(M); ivf <- max(iwf) + 1L
  iJ <- ivf + seq_len(np)
  # BacksolveAdjoint carries the state backward from each segment end
  backsolve <- inherits(sa$adjoint_sensealg, "BacksolveAdjoint")
  iy <- if (backsolve) max(iJ) + seq_len(n) else integer()
  rhs <- function(z, pp, tt) {
    y <- if (backsolve) z[iy] else tr$at(tt)
    L <- matrix(z[il], n, M + 1L)
    dL <- -crossprod(fns$jac(y, pp, tt), L)
    dG <- -crossprod(fns$pjac(y, pp, tt), L)
    dL[, 1L] <- dL[, 1L] - as.numeric(P$gu(y, pp, tt))
    fy <- fns$f(y, pp, tt)
    Lw <- L[, -1L, drop = FALSE]
    dC <- -crossprod(Lw)
    dwv <- -as.numeric(crossprod(Lw, L[, 1L]))
    dwf <- -as.numeric(crossprod(Lw, fy))
    dvf <- -sum(L[, 1L] * fy)
    dJ <- if (!is.null(P$gp)) -as.numeric(P$gp(y, pp, tt)) else numeric(np)
    c(as.numeric(dL), as.numeric(dG), as.numeric(dC), dwv, dwf, dvf, dJ, if (backsolve) fy)
  }
  nz <- max(iJ) + length(iy)
  z0 <- numeric(nz)
  z0[il] <- c(numeric(n), as.numeric(P$W0))
  K <- nseg
  Rs <- array(0, c(M, M, K)); bs <- matrix(0, M, K)
  Cs <- array(0, c(M, M, K)); dwv <- matrix(0, M, K); dwf <- matrix(0, M, K); dvf <- numeric(K)
  dvfs <- matrix(0, np, K); dwfs <- matrix(0, np * M, K); dJs <- matrix(0, np, K)
  tf <- prob$tspan[2L]
  for (iseg in nseg:1L) {
    t1 <- tf - (nseg - iseg + 1L) * P$T_seg
    t2 <- tf - (nseg - iseg) * P$T_seg
    ckpts <- seq(t1, t2, length.out = nstep)
    # BacksolveAdjoint resets the state at the checkpoints (the saved points
    # of the segment): the backward flow of a dissipative system is unstable
    reset <- if (backsolve) {
      z0[iy] <- tr$at(t2)
      DiscreteCallback(function(u, tt, integ) .near_any(ckpts, tt),
                       function(integ) integ$u[iy] <- tr$at(integ$t), save_positions = c(FALSE, FALSE))
    }
    ap <- ODEProblem(rhs, z0, c(t2, t1), p)
    sv <- do.call(solve, c(list(ap, alg), dots, list(save_everystep = FALSE, save_start = FALSE,
                                                     tstops = ckpts, callback = reset)))
    x <- sv$u[[length(sv$t)]]
    L <- matrix(x[il], n, M + 1L)
    vstar <- L[, 1L]
    qrw <- qr(L[, -1L, drop = FALSE])
    Q <- qr.Q(qrw)[, seq_len(M), drop = FALSE]
    Rm <- qr.R(qrw)[seq_len(M), seq_len(M), drop = FALSE]
    Rm[, qrw$pivot] <- Rm
    b <- as.numeric(crossprod(Q, vstar))
    Rs[, , iseg] <- Rm; bs[, iseg] <- b
    G <- matrix(x[ig], np, M + 1L)
    dvfs[, iseg] <- G[, 1L]; dwfs[, iseg] <- as.numeric(G[, -1L])
    Cs[, , iseg] <- matrix(x[iC], M, M); dwv[, iseg] <- x[iwv]; dwf[, iseg] <- x[iwf]; dvf[iseg] <- x[ivf]
    dJs[, iseg] <- x[iJ]
    z0 <- numeric(nz)
    z0[il] <- c(vstar - as.numeric(Q %*% b), as.numeric(Q))
  }
  a <- .nilsas_min(Cs, Rs, bs, dwv, dwf, dvf)
  res <- rowSums(dvfs) + rowSums(dJs)
  for (i in seq_len(K)) res <- res + as.numeric(matrix(dwfs[, i], np, M) %*% a[, i])
  stats::setNames(res / (nseg * P$T_seg), names(prob$p))
}

# nilsas_min: the minimization of NILSAS as one linear solve.
.nilsas_min <- function(C, R, b, dwv, dwf, dvf) {
  M <- dim(C)[1L]; K <- dim(C)[3L]
  Cinv <- matrix(0, M * K, M * K)
  for (i in seq_len(K)) {
    blk <- (i - 1L) * M + seq_len(M)
    Cinv[blk, blk] <- solve(C[, , i])
  }
  B <- matrix(0, M * K - M + 1L, M * K)
  for (i in seq_len(K)) {
    if (i < K) {
      blk <- (i - 1L) * M + seq_len(M)
      B[blk, i * M + seq_len(M)] <- -R[, , i + 1L]
      B[cbind(blk, blk)] <- 1
    }
    B[nrow(B), (i - 1L) * M + seq_len(M)] <- dwf[, i]
  }
  d <- as.numeric(dwv)
  bv <- c(as.numeric(b)[-seq_len(M)], -sum(dvf))
  lam <- solve(-B %*% Cinv %*% t(B), B %*% Cinv %*% d + bv)
  matrix(-Cinv %*% (t(B) %*% lam + d), M, K)
}
