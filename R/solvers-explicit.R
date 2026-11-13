# Explicit Runge-Kutta methods, after OrdinaryDiffEq.jl:
#   lib/OrdinaryDiffEqLowOrderRK/src/low_order_rk_perform_step.jl
#     perform_step! for BS3ConstantCache and DP5ConstantCache
#   lib/OrdinaryDiffEqLowOrderRK/src/interpolants.jl  DP5 dense output (dp5pre0)
#   lib/OrdinaryDiffEqLowOrderRK/src/alg_utils.jl     DP5 controller gains
#   lib/OrdinaryDiffEqVerner/src/verner_rk_perform_step.jl  Vern7, Vern9
#   lib/OrdinaryDiffEqVerner/src/verner_addsteps.jl, interpolants.jl
#     the lazy extra stages and the high order interpolants
#   lib/OrdinaryDiffEqVerner/src/alg_utils.jl  stability sizes
# Tableau constants are in R/solvers-tableaus.R, extracted from the Julia
# sources by a script.

# ---------------------------------------------------------------------------
# BS3: Bogacki-Shampine 3(2), FSAL, cubic Hermite dense output.

.bs3_step <- function(ctx, t, u, dt, cache) {
  T <- .BS3_TAB
  f <- ctx$f; p <- ctx$p
  k1 <- cache$fsal
  k2 <- f(u + dt * T[["a21"]] * k1, p, t + T[["c1"]] * dt)
  k3 <- f(u + dt * T[["a32"]] * k2, p, t + T[["c2"]] * dt)
  unew <- u + dt * (T[["a41"]] * k1 + T[["a42"]] * k2 + T[["a43"]] * k3)
  k4 <- f(unew, p, t + dt)
  err <- dt * (T[["btilde1"]] * k1 + T[["btilde2"]] * k2 + T[["btilde3"]] * k3 + T[["btilde4"]] * k4)
  list(u = unew, err = err, k = list(k1, k4), cache = list(fsal = k4), nf = 3L)
}

# ---------------------------------------------------------------------------
# DP5: Dormand-Prince 5(4), FSAL, with its free 4th order dense output.

.dp5_step <- function(ctx, t, u, dt, cache) {
  T <- .DP5_TAB
  D <- .DP5_DENSE
  f <- ctx$f; p <- ctx$p
  k1 <- cache$fsal
  k2 <- f(u + dt * T[["a21"]] * k1, p, t + T[["c1"]] * dt)
  k3 <- f(u + dt * (T[["a31"]] * k1 + T[["a32"]] * k2), p, t + T[["c2"]] * dt)
  k4 <- f(u + dt * (T[["a41"]] * k1 + T[["a42"]] * k2 + T[["a43"]] * k3), p, t + T[["c3"]] * dt)
  k5 <- f(u + dt * (T[["a51"]] * k1 + T[["a52"]] * k2 + T[["a53"]] * k3 + T[["a54"]] * k4), p,
          t + T[["c4"]] * dt)
  g6 <- u + dt * (T[["a61"]] * k1 + T[["a62"]] * k2 + T[["a63"]] * k3 + T[["a64"]] * k4 + T[["a65"]] * k5)
  k6 <- f(g6, p, t + dt)
  update <- T[["a71"]] * k1 + T[["a73"]] * k3 + T[["a74"]] * k4 + T[["a75"]] * k5 + T[["a76"]] * k6
  unew <- u + dt * update
  k7 <- f(unew, p, t + dt)
  err <- dt * (T[["btilde1"]] * k1 + T[["btilde3"]] * k3 + T[["btilde4"]] * k4 + T[["btilde5"]] * k5 +
                 T[["btilde6"]] * k6 + T[["btilde7"]] * k7)
  bspl <- k1 - update
  kd <- list(update, bspl, update - k7 - bspl,
             D[["d1"]] * k1 + D[["d3"]] * k3 + D[["d4"]] * k4 + D[["d5"]] * k5 + D[["d6"]] * k6 + D[["d7"]] * k7)
  list(u = unew, err = err, k = kd, cache = list(fsal = k7), nf = 6L,
       eigen_est = if (isTRUE(ctx$composite)) .eigen_est(k7 - k6, unew - g6) else NULL)
}

.dp5_interp <- function(theta, dt, y0, y1, k) {
  b20 <- theta * (1 - theta)
  y0 + dt * (k[[1L]] * theta + k[[2L]] * b20 + k[[3L]] * theta * b20 + k[[4L]] * b20^2)
}

# ---------------------------------------------------------------------------
# Verner methods. The Julia field names encode the stage numbers: aIJ with
# three digits is row II, column J (Vern7 main tableau), with four digits
# row II, column JJ; bI, btildeI weight stage I; rIIP is the coefficient of
# theta^P in the interpolation weight of stage II. The abscissae follow the
# perform_step! code: Vern7 stage i uses c_i, Vern9 stage i uses c_(i-1),
# and the last two main stages use t + dt.

.verner_method <- function(main, extra, interp, nmain, c_main) {
  all <- c(main, extra)
  nstage <- nmain + max(as.integer(sub("^c", "", grep("^c[0-9]+$", names(extra), value = TRUE))) - nmain)
  A <- matrix(0, nstage, nstage)
  for (nm in grep("^a[0-9]+$", names(all), value = TRUE)) {
    dg <- sub("^a", "", nm)
    i <- as.integer(substr(dg, 1L, 2L))
    j <- as.integer(substr(dg, 3L, nchar(dg)))
    A[i, j] <- all[[nm]]
  }
  cc <- numeric(nstage)
  cc[2:nmain] <- c_main
  for (i in (nmain + 1L):nstage) cc[i] <- extra[[paste0("c", i)]]
  vec <- function(prefix) {
    v <- numeric(nstage)
    for (nm in grep(paste0("^", prefix, "[0-9]+$"), names(main), value = TRUE))
      v[as.integer(sub(prefix, "", nm))] <- main[[nm]]
    v
  }
  R <- matrix(0, nstage, 9L)
  for (nm in grep("^r[0-9]+$", names(interp), value = TRUE)) {
    R[as.integer(substr(nm, 2L, 3L)), as.integer(substr(nm, 4L, 4L))] <- interp[[nm]]
  }
  list(A = A, c = cc, b = vec("b"), btilde = vec("btilde"), R = R, nmain = nmain, nstage = nstage)
}

# Built on first use: R/solvers-tableaus.R is collated after this file.
.verner_tabs <- new.env(parent = emptyenv())
.verner_tab <- function(name) {
  tab <- .verner_tabs[[name]]
  if (is.null(tab)) {
    tab <- switch(name,
      Vern6 = .verner_method(.vern6_main(), .VERN6_EXTRA, .VERN6_INTERP, 9L,
                             c(.VERN6_TAB[paste0("c", 1:6)], 1, 1)),
      Vern7 = .verner_method(.VERN7_TAB, .VERN7_EXTRA, .VERN7_INTERP, 10L,
                             c(.VERN7_TAB[paste0("c", 2:8)], 1, 1)),
      Vern8 = .verner_method(.VERN8_TAB, .VERN8_EXTRA, .VERN8_INTERP, 13L,
                             c(.VERN8_TAB[paste0("c", 2:11)], 1, 1)),
      Vern9 = .verner_method(.VERN9_TAB, .VERN9_EXTRA, .VERN9_INTERP, 16L,
                             c(.VERN9_TAB[paste0("c", 1:13)], 1, 1)))
    tab$name <- name
    assign(name, tab, envir = .verner_tabs)
  }
  tab
}

# Vern6 names its main coefficients aIJ with one digit each and computes
# the solution from row 9 (a9J), whose stage is the FSAL derivative at the
# new state (perform_step! for Vern6ConstantCache): rewritten to the aIIJ
# names of Vern7 with row 9 as the weights bJ.
.vern6_main <- function() {
  v <- .VERN6_TAB
  nm <- names(v)
  a2 <- grepl("^a[0-9]{2}$", nm)
  nm[a2] <- sub("^a([0-9])([0-9])$", "a0\\1\\2", nm[a2])
  v <- stats::setNames(v, nm)
  w <- v[grepl("^a09[0-9]$", names(v))]
  c(v, stats::setNames(w, sub("^a09", "b", names(w))))
}

# The sum of w[j] * K[, j] over the nonzero weights, term by term in stage
# order. OrdinaryDiffEq writes these sums out (a1 * k1 + a3 * k3 + ...) and
# its out of place muladd on vectors does not fuse, whereas R's BLAS (%*%)
# fuses multiply and add; the plain sum gives Julia's bits.
.lincomb <- function(K, w, idx = which(w != 0)) {
  if (!length(idx)) return(numeric(nrow(K)))
  acc <- K[, idx[1L]] * w[idx[1L]]
  for (j in idx[-1L]) acc <- acc + K[, j] * w[j]
  acc
}

# The input of stage i: uprev + (dt * a21) * k1 for the second stage, else
# uprev + dt * (sum of a_ij k_j), as the Verner perform_step! methods.
.rk_stage_input <- function(tab, K, i, u, dt) {
  nz <- tab$rows[[i]]
  if (i == 2L && identical(nz, 1L)) return(u + (dt * tab$A[2L, 1L]) * K[, 1L])
  u + dt * .lincomb(K, tab$A[i, ], nz)
}

# Stages from..to of an explicit tableau, given the earlier columns of K;
# the stage times t + c dt are muladd (fused) on scalars.
.rk_stages <- function(tab, K, from, to, f, p, t, u, dt) {
  for (i in from:to) K[, i] <- f(.rk_stage_input(tab, K, i, u, dt), p, .fma(tab$c[i], dt, t))
  K
}

.verner_step <- function(tab, fsal = FALSE) {
  force(tab)
  tab$rows <- lapply(seq_len(tab$nstage), function(i) which(tab$A[i, seq_len(i - 1L)] != 0))
  ib <- which(tab$b[seq_len(tab$nmain)] != 0)
  ibt <- which(tab$btilde[seq_len(tab$nmain)] != 0)
  function(ctx, t, u, dt, cache) {
    f <- ctx$f; p <- ctx$p
    nm <- tab$nmain
    K <- matrix(0, length(u), tab$nstage)
    # Vern6 is FSAL: its first stage is the last stage of the previous step
    K[, 1L] <- if (fsal && !is.null(cache$fsal)) cache$fsal else f(u, p, t)
    K <- .rk_stages(tab, K, 2L, nm, f, p, t, u, dt)
    unew <- u + dt * .lincomb(K, tab$b, ib)
    err <- dt * .lincomb(K, tab$btilde, ibt)
    eig <- NULL
    if (isTRUE(ctx$composite))
      eig <- .eigen_est(K[, nm] - K[, nm - 1L], .rk_stage_input(tab, K, nm, u, dt) - .rk_stage_input(tab, K, nm - 1L, u, dt))
    # the extra interpolation stages are computed on first use (lazy = TRUE)
    lazy <- new.env(parent = emptyenv())
    lazy$K <- K
    lazy$done <- FALSE
    lazy$args <- list(f = f, p = p, t = t, u = u, dt = dt, counter = ctx$lazy_nf)
    list(u = unew, err = err, k = list(tab = tab, lazy = lazy),
         cache = if (fsal) list(fsal = K[, nm]) else list(), nf = if (fsal) nm - 1L else nm, eigen_est = eig)
  }
}

# The dense output: the weight of stage 1 is theta times a polynomial with
# the coefficients of the powers 1 to p, the others theta^2 times one from
# power 2 (@evalpoly: Horner with fused muladd), then
# y0 + dt * (sum of k_j b_j(theta)) term by term.
.verner_interp <- function(theta, dt, y0, y1, k) {
  tab <- k$tab
  lz <- k$lazy
  if (!lz$done) {
    a <- lz$args
    lz$K <- .rk_stages(tab, lz$K, tab$nmain + 1L, tab$nstage, a$f, a$p, a$t, a$u, a$dt)
    lz$done <- TRUE
    .count_lazy_nf(a$counter, tab$nstage - tab$nmain)
  }
  th2 <- theta * theta
  acc <- NULL
  for (j in which(rowSums(tab$R != 0) > 0)) {
    m <- if (j == 1L) 1L else 2L
    pmax <- max(which(tab$R[j, ] != 0))
    h <- tab$R[j, pmax]
    if (pmax > m) for (q in (pmax - 1L):m) h <- .fma(theta, h, tab$R[j, q])
    term <- lz$K[, j] * ((if (m == 1L) theta else th2) * h)
    acc <- if (is.null(acc)) term else acc + term
  }
  y0 + dt * acc
}

# Vern7, Vern8 and Vern9 evaluate nothing in initialize! (no FSAL).
.verner_init <- function(ctx, t, u) list(nf0 = 0L)

# The extra stages count as function evaluations: during the solve they are
# added to stats$nf when the solution is built, afterwards to the
# solution's stats (counter set up by .integrate). Julia's _ode_addsteps!
# calls f without counting.
.count_lazy_nf <- function(counter, n) {
  if (is.null(counter)) return(invisible())
  env <- counter$sol_env
  if (is.null(env)) counter$n <- counter$n + as.integer(n)
  else env$stats$nf <- env$stats$nf + as.integer(n)
  invisible()
}

#' @rdname ode_algorithms
#' @export
BS3 <- function() {
  .new_alg("BS3", 3L, step = .bs3_step, stability_size = NULL)
}

#' @rdname ode_algorithms
#' @export
DP5 <- function() {
  .new_alg("DP5", 5L, step = .dp5_step, interp = .dp5_interp, beta2 = 4 / 100,
           beta1 = 1 / 5 - 3 * (4 / 100) / 4, stability_size = 3.3066)
}

#' @rdname ode_algorithms
#' @export
Vern7 <- function() {
  .new_alg("Vern7", 7L, step = .verner_step(.verner_tab("Vern7")), interp = .verner_interp,
           init = .verner_init, stability_size = 4.64)
}

#' @rdname ode_algorithms
#' @export
Vern6 <- function() {
  .new_alg("Vern6", 6L, step = .verner_step(.verner_tab("Vern6"), fsal = TRUE), interp = .verner_interp,
           init = function(ctx, t, u) list(fsal = ctx$f(u, ctx$p, t)), stability_size = 4.8553)
}

#' @rdname ode_algorithms
#' @export
Vern8 <- function() {
  .new_alg("Vern8", 8L, step = .verner_step(.verner_tab("Vern8")), interp = .verner_interp,
           init = .verner_init, stability_size = 5.8641)
}

#' @rdname ode_algorithms
#' @export
Vern9 <- function() {
  .new_alg("Vern9", 9L, step = .verner_step(.verner_tab("Vern9")), interp = .verner_interp,
           init = .verner_init, stability_size = 4.4762)
}
