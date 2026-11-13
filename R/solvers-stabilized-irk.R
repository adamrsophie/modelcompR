# IRKC of OrdinaryDiffEqStabilizedIRK (algorithms.jl, alg_utils.jl,
# irkc_caches.jl, irkc_perform_step.jl, irkc_utils.jl): the implicit
# Runge-Kutta-Chebyshev method for split problems du/dt = f1(u) + f2(u),
# with the stiff part f1 solved by Newton's method at every stage (DIRK
# form, Jacobian of f1) and f2 treated explicitly through the shifted
# Chebyshev recursion.

#' Implicit Runge-Kutta-Chebyshev method
#'
#' `IRKC()` integrates a [SplitODEProblem()] `du/dt = f1(u, p, t) + f2(u, p,
#' t)`: each of its Chebyshev stages solves the stiff part `f1` implicitly
#' (modified Newton, with the Jacobian of `f1`: `SplitFunction(f1, f2, jac
#' = J1)`, else finite differences) and takes the mildly stiff part `f2`
#' explicitly. Second order, with an error estimate of first order and the
#' integral controller (`gamma = 0.8`).
#'
#' As in OrdinaryDiffEq, the cache fixes the number of stages at 50 (its
#' smallest and largest degree are both 50), so every step takes 50 stages
#' whatever the spectral radius of `f2`; the radius, by the power iteration
#' on `f2` every step or from `eigen_est`, only enters the stage count
#' formula. `sol$stats$nf2` counts the evaluations of `f2`.
#'
#' @param eigen_est An upper bound on the spectral radius of the Jacobian
#'   of `f2`: a number or a function `function(u, p, t)`. By default the
#'   nonlinear power iteration on `f2` estimates it every step.
#' @param autodiff As for [ImplicitEuler()]: how the Jacobian of `f1` is
#'   formed when none is given.
#' @return An ODE algorithm.
#' @examples
#' # a stiff linear decay (implicit) with a mild explicit forcing
#' prob <- SplitODEProblem(function(u, p, t) -50 * u, function(u, p, t) -u + cos(t),
#'                         c(1, 2), c(0, 1))
#' sol <- solve(prob, IRKC(), reltol = 1e-4, abstol = 1e-6)
#' sol$u[[length(sol$u)]]
#' @export
IRKC <- function(eigen_est = NULL, autodiff = NULL) .with_autodiff(autodiff, {
  if (!is.null(eigen_est) && !is.numeric(eigen_est) && !is.function(eigen_est))
    stop("eigen_est must be a number or a function(u, p, t)", call. = FALSE)
  .new_alg("IRKC", 2L, step = function(ctx, t, u, dt, cache) .irkc_step(ctx, t, u, dt, cache, eigen_est),
           init = .irkc_init, on_event = .irkc_on_event,
           stiff = TRUE, controller = "I", adaptive_order = 1L, gamma = 8 / 10,
           eigen_est = eigen_est, on_reject = .nl_on_reject)
})

# initialize!: du1 = f1(u), du2 = f2(u), fsalfirst = du1 + du2.
.irkc_init <- function(ctx, t, u) {
  if (is.null(ctx$f1) || is.null(ctx$f2))
    stop("IRKC() integrates a SplitODEProblem (f1 implicit, f2 explicit)", call. = FALSE)
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  cache$minm <- 50L
  cache$st <- new.env(parent = emptyenv())
  cache$st$zprev <- u
  cache$st$eig <- 0
  .irkc_derivs(ctx, cache, t, u)
  cache$nf0 <- 1L
  cache$nf20 <- 1L
  cache
}

.irkc_derivs <- function(ctx, cache, t, u) {
  cache$du1 <- ctx$f1(u, ctx$p, t)
  cache$du2 <- ctx$f2(u, ctx$p, t)
  cache$du_u <- .as_num(u)
  cache$du_stale <- FALSE
  cache$fsal <- cache$du1 + cache$du2
}

# After an event the Newton state stays; the derivatives are evaluated at
# the new state when the next step starts.
.irkc_on_event <- function(ctx, t, u, cache) {
  cache$du_stale <- TRUE
  list(cache = cache, nf = 0L)
}

.irkc_step <- function(ctx, t, u, dt, cache, eigen_est = NULL) {
  f1 <- ctx$f1; f2 <- ctx$f2; p <- ctx$p
  nl <- cache$nl
  nf <- 0L; nf2 <- 0L
  cnt <- list(nf = 0L, njacs = 0L, nsolve = 0L, nw = 0L)
  add <- function(x) for (s in names(cnt)) cnt[[s]] <<- cnt[[s]] + as.integer(x[[s]] %||% 0L)
  # du1 and du2 belong to the start of the step: after an event or a step
  # limiter they are evaluated again (OrdinaryDiffEq keeps the values of the
  # last attempt, see PORTING.md)
  if (isTRUE(cache$du_stale) || !identical(cache$du_u, .as_num(u))) {
    .irkc_derivs(ctx, cache, t, u)
    nf <- nf + 1L
    nf2 <- nf2 + 1L
  }
  du1 <- cache$du1; du2 <- cache$du2
  # maxeig! on f2 every step, or the user's bound
  st <- cache$st
  e <- eigen_est
  if (!is.null(e)) {
    st$eig <- .as_num(if (is.function(e)) e(u, p, t) else e)
  } else {
    st$disc <- isTRUE(ctx$derivative_discontinuity)
    ctx2 <- ctx
    ctx2$f <- f2
    nf2 <- nf2 + .stab_maxeig(ctx2, st, t, u, du2, TRUE, 50L)
  }
  maxm <- 50L
  mdeg <- 1 + floor(sqrt(1.54 * abs(dt) * st$eig + 1))
  mdeg <- if (is.finite(mdeg)) min(maxm, max(cache$minm, mdeg)) else maxm
  w0 <- 1 + 2 / (13 * (mdeg * mdeg))
  temp1 <- w0 * w0 - 1
  temp2 <- sqrt(temp1)
  theta <- mdeg * log(w0 + temp2)
  w1 <- (sinh(theta) * temp1) / (cosh(theta) * mdeg * temp2 - w0 * sinh(theta))
  Bjm2 <- 1 / (4 * (w0 * w0))
  Bjm1 <- 1 / w0
  # the implicit part: Newton on f1 with its Jacobian
  ctx1 <- ctx
  ctx1$f <- f1
  ctx1$jac <- ctx$jac1
  ctx1$tgrad <- NULL
  # stage 1
  f1jm2 <- du1
  gprev2 <- u
  mus <- w1 * Bjm1
  mus1 <- mus
  tmp <- u + dt * mus1 * du2
  r <- .nlsolve(nl, ctx1, t, dt, dt * du1, tmp, mus1, 1, t + mus * dt, u, method = "dirk")
  add(r)
  gprev <- tmp + mus1 * r$z
  ok <- r$ok
  Cjm2 <- 0; Cjm1 <- mus
  Tjm1 <- w0; Tjm2 <- 1
  Tjm1p <- 1; Tjm2p <- 0
  Tjm1pp <- 0; Tjm2pp <- 0
  unew <- gprev
  for (j in seq_len(mdeg - 1L) + 1L) {
    Tj <- 2 * w0 * Tjm1 - Tjm2
    Tjp <- 2 * w0 * Tjm1p + 2 * Tjm1 - Tjm2p
    Tjpp <- 2 * w0 * Tjm1pp + 4 * Tjm1p - Tjm2pp
    Bj <- Tjpp / (Tjp * Tjp)
    mu <- (2 * w0 * Bj) / Bjm1
    nu <- -Bj / Bjm2
    mus <- (mu * w1) / w0
    nus <- -(1 - Tjm1 * Bjm1) * mus
    Cj <- mu * Cjm1 + nu * Cjm2 + mus + nus
    f1jm1 <- f1(gprev, p, t + Cjm1 * dt)
    f2jm1 <- f2(gprev, p, t + Cjm1 * dt)
    nf <- nf + 1L
    nf2 <- nf2 + 1L
    tmp <- (1 - mu - nu) * u + mu * gprev + nu * gprev2 + dt * mus * f2jm1 + dt * nus * du2 +
      (nus - (1 - mu - nu) * mus1) * dt * du1 - nu * mus1 * dt * f1jm2
    # the Newton failures of the stages are ignored; the last one decides
    # whether the step fails (force_stepfail)
    r <- .nlsolve(nl, ctx1, t, dt, dt * f1jm1, tmp, mus1, 1, t + Cj * dt, u, method = "dirk",
                  firststage = FALSE)
    add(r)
    ok <- r$ok
    unew <- tmp + mus1 * r$z
    if (j < mdeg) {
      f1jm2 <- f1jm1
      gprev2 <- gprev
      gprev <- unew
      Cjm2 <- Cjm1; Cjm1 <- Cj
      Bjm2 <- Bjm1; Bjm1 <- Bj
      Tjm2 <- Tjm1; Tjm1 <- Tj
      Tjm2p <- Tjm1p; Tjm1p <- Tjp
      Tjm2pp <- Tjm1pp; Tjm1pp <- Tjpp
    }
  }
  du1n <- f1(unew, p, t + dt)
  du2n <- f2(unew, p, t + dt)
  nf <- nf + 1L
  nf2 <- nf2 + 1L
  EEst <- 1
  if (ctx$adaptive) {
    # W = J - I / dt, with do_newJW measuring against mus1 * dt
    w <- .update_W(nl, ctx1, t, u, dt, gdt_test = mus1 * dt)
    add(w)
    est <- dt * (0.5 * (du2n - du2) + (0.5 - mus1) * (du1n - du1))
    EEst <- if (is.null(nl$W)) Inf else .err_norm(.wsolve(nl$W, est), u, unew, ctx$abstol, ctx$reltol)
  }
  fnew <- du1n + du2n
  out <- list(nf = cnt$nf + nf, nf2 = nf2, njacs = cnt$njacs, nsolve = cnt$nsolve, nw = cnt$nw,
              cache = cache)
  if (!ok) return(c(out, list(fail = TRUE, u = u, EEst = EEst)))
  if (isTRUE(EEst <= 1)) {
    cache$du1 <- du1n
    cache$du2 <- du2n
    cache$du_u <- .as_num(unew)
    cache$fsal <- fnew
  }
  c(out, list(u = unew, EEst = EEst, k = list(du1 + du2, fnew)))
}
