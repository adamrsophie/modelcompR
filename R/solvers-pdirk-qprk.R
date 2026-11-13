# KuttaPRK2p5 of OrdinaryDiffEqPRK, QPRK98 of OrdinaryDiffEqQPRK and
# PDIRK44 of OrdinaryDiffEqPDIRK (algorithms.jl, alg_utils.jl, *_caches.jl,
# *_perform_step.jl, qprk_tableaus.jl): the out of place ConstantCache
# paths. Stage sums follow Julia's @muladd expression shapes: a vector
# muladd is not fused, a scalar one (the stage times t + c * dt) is.

#' Parallel explicit Runge-Kutta methods
#'
#' `KuttaPRK2p5()` is the fixed step fifth order method of Jackson and
#' Norsett with five stages arranged for two processors (the last two stages
#' are independent). `QPRK98()` is the adaptive 9(8) pair of Kovalnogov,
#' Fedorov, Karpukhina, Simos and Tsitouras, designed for quadruple precision
#' (here in double precision), with the PI controller of order 9 and FSAL.
#'
#' The stages are evaluated one after the other; `threading` and `thread`
#' are accepted for compatibility and do not change the results. As in
#' OrdinaryDiffEq's out of place method, `KuttaPRK2p5()` counts no function
#' evaluations in `sol$stats$nf`.
#'
#' @param threading Ignored (Julia's threading option of the parallel stages).
#' @param stage_limiter Stored only: as in Julia, the stages of `QPRK98()`
#'   do not call it.
#' @param step_limiter Optional `function(u, integrator, p, t)` applied after
#'   every accepted step.
#' @param thread Ignored.
#' @return An ODE algorithm.
#' @examples
#' prob <- ODEProblem(function(u, p, t) c(u[2], -u[1]), c(0, 1), c(0, 1))
#' solve(prob, KuttaPRK2p5(), dt = 0.1)$u[[11]]
#' solve(prob, QPRK98(), reltol = 1e-10, abstol = 1e-10)$u[[2]]
#' @name parallel_rk_methods
NULL

#' @rdname parallel_rk_methods
#' @export
KuttaPRK2p5 <- function(threading = TRUE) {
  .new_alg("KuttaPRK2p5", 5L, adaptive = FALSE, step = .prk_step,
           init = function(ctx, t, u) list(fsal = ctx$f(u, ctx$p, t), nf0 = 0L),
           # a callback that changes u evaluates the FSAL derivative again
           # (counted, unlike the evaluations of the steps)
           on_event = function(ctx, t, u, cache) list(cache = list(fsal = ctx$f(u, ctx$p, t), nf0 = 0L), nf = 1L))
}

# perform_step! of KuttaPRK2p5ConstantCache: no evaluation is counted.
.prk_step <- function(ctx, t, u, dt, cache) {
  f <- ctx$f; p <- ctx$p
  k1 <- f(u, p, t)
  k2 <- f(u + (dt * (1 / 3)) * k1, p, .fma(1 / 3, dt, t))
  k3 <- f(u + dt * ((4 / 25) * k1 + (6 / 25) * k2), p, .fma(2 / 5, dt, t))
  k4 <- f(u + dt * ((1 / 4) * k1 + (-3) * k2 + (15 / 4) * k3), p, .fma(1, dt, t))
  k5 <- f(u + dt * ((6 / 81) * k1 + (90 / 81) * k2 + (-50 / 81) * k3 + (8 / 81) * k4), p,
          .fma(2 / 3, dt, t))
  k6 <- f(u + dt * ((6 / 75) * k1 + (36 / 75) * k2 + (10 / 75) * k3 + (8 / 75) * k4), p,
          .fma(4 / 5, dt, t))
  unew <- u + dt * ((23 / 192) * k1 + (125 / 192) * k3 + (-81 / 192) * k5 + (125 / 192) * k6)
  knew <- f(unew, p, t + dt)
  list(u = unew, k = list(cache$fsal, knew), cache = list(fsal = knew, nf0 = 0L), nf = 0L)
}

#' @rdname parallel_rk_methods
#' @export
QPRK98 <- function(stage_limiter = NULL, step_limiter = NULL, thread = NULL) {
  .new_alg("QPRK98", 9L, step = .qprk_step, step_limiter = step_limiter,
           stage_limiter = stage_limiter)
}

# The stages of QPRK98: the indices of the earlier stages each one uses, in
# Julia's summation order.
.QPRK98_ROWS <- list(NULL, 1L, 1:2, c(1L, 3L), c(1L, 3L, 4L), c(1L, 4L, 5L), c(1L, 4L, 5L, 6L),
                     c(1L, 6L, 7L), c(1L, 6L, 7L, 8L), c(1L, 6:9), c(1L, 6:10), c(1L, 6:11),
                     c(1L, 6:12), c(1L, 6:13), c(1L, 6:14), c(1L, 6:14))

.qprk_coef <- function(i, j) {
  tb <- .QPRK98_TAB
  tb[[if (i < 10L) paste0("b", i, j) else paste0("b", i, "_", j)]]
}

# perform_step! of QPRK98ConstantCache.
.qprk_step <- function(ctx, t, u, dt, cache) {
  f <- ctx$f; p <- ctx$p; tb <- .QPRK98_TAB
  K <- matrix(0, length(u), 16L)
  K[, 1L] <- cache$fsal
  K[, 2L] <- f(u + (tb$b21 * K[, 1L]) * dt, p, .fma(tb$d2, dt, t))
  for (i in 3:16) {
    idx <- .QPRK98_ROWS[[i]]
    w <- numeric(16L)
    w[idx] <- vapply(idx, function(j) .qprk_coef(i, j), 0)
    ti <- if (i <= 14L) .fma(tb[[paste0("d", i)]], dt, t) else t + dt
    K[, i] <- f(u + dt * .lincomb(K, w, idx), p, ti)
  }
  wi <- c(1L, 8:16)
  w <- numeric(16L); w[wi] <- vapply(wi, function(j) tb[[paste0("w", j)]], 0)
  unew <- u + dt * .lincomb(K, w, wi)
  ei <- c(1L, 8:15)
  e <- numeric(16L); e[ei] <- vapply(ei, function(j) tb[[paste0("eps", j)]], 0)
  utilde <- dt * .lincomb(K, e, ei)
  EEst <- .firk_err(utilde, u, unew, ctx$abstol, ctx$reltol)
  fsal <- f(unew, p, t + dt)
  list(u = unew, EEst = EEst, k = list(K[, 1L], fsal), cache = list(fsal = fsal), nf = 16L)
}

#' Parallel diagonally implicit Runge-Kutta method
#'
#' `PDIRK44()` is the fixed step fourth order parallel DIRK method of
#' OrdinaryDiffEqPDIRK: two rounds of two independent implicit stages
#' (diagonal coefficients 1/2 and 2/3), each solved by the simplified Newton
#' method with a fresh Jacobian and factorization at every stage. With
#' `threading = TRUE` (the default) each of the two stage columns keeps its
#' own nonlinear solver state, as Julia's two `NLSolver`s; with `FALSE` one
#' solver serves all four stages. The stages are evaluated in order either
#' way.
#'
#' @param autodiff As for [ImplicitEuler()]: how the Jacobian is formed
#'   when the problem gives none.
#' @param threading Whether the two stage columns have separate nonlinear
#'   solvers (Julia's threaded branch).
#' @return An ODE algorithm.
#' @examples
#' prob <- ODEProblem(function(u, p, t) c(-50 * (u[1] - cos(t)), -u[2]), c(0, 1), c(0, 1))
#' solve(prob, PDIRK44(), dt = 0.05)$u[[21]]
#' @export
PDIRK44 <- function(autodiff = NULL, threading = TRUE) .with_autodiff(autodiff, {
  .new_alg("PDIRK44", 4L, adaptive = FALSE, step = .pdirk_step, interp = .newmark_interp,
           init = function(ctx, t, u) .pdirk_init(ctx, t, u, threading),
           on_event = function(ctx, t, u, cache) list(cache = cache, nf = 0L), stiff = TRUE)
})

.pdirk_init <- function(ctx, t, u, threading) {
  cache <- new.env(parent = emptyenv())
  cache$nl <- if (isTRUE(threading)) list(.new_nlsolver(), .new_nlsolver()) else list(.new_nlsolver())
  cache$nf0 <- 0L
  cache
}

# perform_step! of PDIRK44ConstantCache. Each stage solves
# z = dt f(tmp + gamma z, t + c dt) (DIRK form, zero initial guess).
.pdirk_step <- function(ctx, t, u, dt, cache) {
  gs <- c(1 / 2, 2 / 3); cs <- c(1 / 2, 2 / 3, 1 / 2, 1 / 3)
  a1 <- c(-5 / 2, -5 / 3); a2 <- c(5 / 2, 4 / 3)
  cnt <- list(nf = 0L, njacs = 0L, nsolve = 0L, nw = 0L)
  # isadaptive(PDIRK44) is false: nlsolve! does not trust the convergence
  # rate on the first Newton iteration
  solve1 <- function(nl, tmp, g, c, fs) {
    r <- .nlsolve(nl, ctx, t, dt, 0 * u, tmp, g, 1, .fma(c, dt, t), u, method = "dirk", firststage = fs, eta_first = FALSE)
    for (s in names(cnt)) cnt[[s]] <<- cnt[[s]] + as.integer(r[[s]] %||% 0L)
    r
  }
  failed <- function() c(list(u = u, fail = TRUE, nreject = 0L, cache = cache), cnt)
  if (length(cache$nl) == 2L) {
    # threaded: stage column i keeps nlsolver[i]; the second round is not
    # marked as a first stage
    r1 <- lapply(1:2, function(i) solve1(cache$nl[[i]], u, gs[i], cs[i], TRUE))
    if (!r1[[1L]]$ok || !r1[[2L]]$ok) return(failed())
    k1 <- lapply(r1, `[[`, "z")
    r2 <- lapply(1:2, function(i) solve1(cache$nl[[i]], u + a1[i] * k1[[1L]] + a2[i] * k1[[2L]], gs[i], cs[2L + i], FALSE))
    if (!r2[[1L]]$ok || !r2[[2L]]$ok) return(failed())
    k2 <- lapply(r2, `[[`, "z")
  } else {
    nl <- cache$nl[[1L]]
    k1 <- list(); k2 <- list()
    for (i in 1:2) {
      r <- solve1(nl, u, gs[i], cs[i], TRUE)
      if (!r$ok) return(failed())
      k1[[i]] <- r$z
    }
    for (i in 1:2) {
      r <- solve1(nl, u + a1[i] * k1[[1L]] + a2[i] * k1[[2L]], gs[i], cs[2L + i], TRUE)
      if (!r$ok) return(failed())
      k2[[i]] <- r$z
    }
  }
  unew <- u + (-1) * k1[[1L]] + (-1) * k2[[1L]] + (3 / 2) * k1[[2L]] + (3 / 2) * k2[[2L]]
  c(list(u = unew, k = list(f = ctx$f, p = ctx$p, t = t), cache = cache), cnt)
}

# QPRK98Tableau(Float64, Float64), printed by Julia with repr (the BigInt
# rationals of qprk_tableaus.jl rounded to Float64).
.QPRK98_TAB <- list(
  d2 = 0.016304347826086956,
  d3 = 0.03865889883769841,
  d4 = 0.057988348256547605,
  d5 = 0.9880697460997246,
  d6 = 0.08767153884170796,
  d7 = 0.2086403684291034,
  d8 = 0.24692658939234283,
  d9 = 0.1851949420442571,
  d10 = 0.23412654896370216,
  d11 = 0.6160965794768611,
  d12 = 0.8330182682406116,
  d13 = 0.9995,
  d14 = 0x1.ffec202405bebp-1,
  eps1 = 0x1.b29bfc2d63797p-14,
  eps8 = -0.031585807132642474,
  eps9 = -0.006741120724439175,
  eps10 = 0.037789967808180114,
  eps11 = 0.0012404893989179156,
  eps12 = -0.002985857841439465,
  eps13 = 0x1.b4a87a783fe0dp+6,
  eps14 = -358.28478752540417,
  eps15 = 249.12243655003593,
  w1 = 0.03875543106489954,
  w8 = 5.863446352480979,
  w9 = 1.4351778532846822,
  w10 = -6.869189080199178,
  w11 = 0.28048709689195883,
  w12 = 0.14896062923887027,
  w13 = 1097.4353761326668,
  w14 = -3588.2929966912766,
  w15 = 2490.498318858044,
  w16 = 0.4616634178037687,
  b21 = 0.016304347826086956,
  b31 = -0.007172755248832501,
  b32 = 0.04583165408653091,
  b41 = 0x1.db0a6169f324bp-7,
  b43 = 0.04349126119241071,
  b51 = 0x1.ed8213d52c345p+6,
  b53 = -392.4206182220937,
  b54 = 270.031659187425,
  b61 = 0.021427839322007335,
  b64 = 0.0662417785653732,
  b65 = 1.9209543274403406e-6,
  b71 = 0.15197913719183012,
  b74 = -0.5606462658457096,
  b75 = 0.00017331706423201162,
  b76 = 0.6171341800187509,
  b81 = 0.027436287710260313,
  b86 = 0.12654637717262396,
  b87 = 0.09294392450945856,
  b91 = 0.027489874209694416,
  b96 = 0.12611939341068437,
  b97 = 0x1.6d6c10cf521b5p-5,
  b98 = -0.013021519362486828,
  b10_1 = 0.02745986119504744,
  b10_6 = 0.1263572516595024,
  b10_7 = 0.08601756786180088,
  b10_8 = -0.009070455381866497,
  b10_9 = 0.0033623236292179364,
  b11_1 = 2.873495572211196,
  b11_6 = -15.861151603498543,
  b11_7 = 2.406702715992152,
  b11_8 = 162.5864745953057,
  b11_9 = 74.07106779134126,
  b11_10 = -225.4604924918749,
  b12_1 = -0x1.f17851133bf8fp+5,
  b12_6 = 343.8261357134057,
  b12_7 = -58.181454546277664,
  b12_8 = -3261.213133835117,
  b12_9 = -0x1.807e2b726fc0fp+10,
  b12_10 = 4574.0602299311595,
  b12_11 = 0x1.3f89bf19abd70p+1,
  b13_1 = 760.8157537655514,
  b13_6 = -4201.7573404025625,
  b13_7 = 716.6942980313362,
  b13_8 = 39888.22388083239,
  b13_9 = 18799.064323744275,
  b13_10 = -55936.71343730077,
  b13_11 = -0x1.a894dbf709a0bp+4,
  b13_12 = 1.2083639172512697,
  b14_1 = 764.9891298171596,
  b14_6 = -4224.816316208864,
  b14_7 = 720.6210614061454,
  b14_8 = 40107.36167774771,
  b14_9 = 18902.29593416604,
  b14_10 = -56243.98051773594,
  b14_11 = -26.68576221023212,
  b14_12 = 0x1.36fe44ef8f17ep+0,
  b14_13 = -0.00017595335744664912,
  b15_1 = 766.8104383772252,
  b15_6 = -4234.879525577625,
  b15_7 = 722.3347429905238,
  b15_8 = 40202.99629756439,
  b15_9 = 18947.34750830305,
  b15_10 = -56378.07587008917,
  b15_11 = -26.75097212246604,
  b15_12 = 1.2176325822640046,
  b15_13 = -3.3955779181668334e-5,
  b15_14 = -0.00021807241143683643,
  b16_1 = 709.1541214095432,
  b16_6 = -3916.568062852313,
  b16_7 = 667.8932620361838,
  b16_8 = 37184.914892338886,
  b16_9 = 17524.222234648252,
  b16_10 = -52145.022470848264,
  b16_11 = -0x1.8bfb2f6b51a2ap+4,
  b16_12 = 1.1528369625257169,
  b16_13 = 0.004142192536637572,
  b16_14 = -0.002131367987008277
)
