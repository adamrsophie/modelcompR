# IMEX multistep methods of OrdinaryDiffEqIMEXMultistep (algorithms.jl,
# alg_utils.jl, imex_multistep_caches.jl, imex_multistep_perform_step.jl):
# CNAB2 and CNLF2 for split problems du/dt = f1(u, p, t) + f2(u, p, t), with
# f1 treated implicitly (Crank-Nicolson, NLNewton in the DIRK form on f1)
# and f2 explicitly. Both are fixed step two step methods; the out of place
# perform_step! methods are not @muladd, so no multiply add is fused.

#' IMEX multistep methods
#'
#' Fixed step second order implicit explicit two step methods for a
#' [SplitODEProblem()] `du/dt = f1(u, p, t) + f2(u, p, t)`: the first part
#' `f1` is treated implicitly by the Crank-Nicolson rule (Newton's method
#' with the Jacobian of `f1`: `SplitFunction(f1, f2, jac = J1)`, else finite
#' differences), the second part `f2` explicitly.
#'
#' * `CNAB2()`: Crank-Nicolson Adams-Bashforth (He and Li 2010).
#' * `CNLF2()`: Crank-Nicolson leapfrog (Jiang, Kubacki, Layton, Moraiti and
#'   Tran 2015).
#'
#' The first step, and the first step after an event that changes the
#' state, use the one step startup formula. Give `dt`; the methods are not
#' adaptive. As in OrdinaryDiffEq, `f2` enters through the difference of
#' the derivative `f1 + f2` carried from the previous step and `f1`.
#'
#' @param autodiff As for [ImplicitEuler()]: how the Jacobian of `f1` is
#'   formed when none is given.
#' @param extrapolant The initial guess of the Newton iteration; only
#'   `"linear"` (the default of OrdinaryDiffEq, which starts from
#'   `dt * f1(u)`) is used by these methods.
#' @return An ODE algorithm.
#' @examples
#' prob <- SplitODEProblem(function(u, p, t) -9 * u, function(u, p, t) -u, 1, c(0, 1))
#' sol <- solve(prob, CNAB2(), dt = 0.05)
#' sol$u[[length(sol$u)]]
#' solve(prob, CNLF2(), dt = 0.05)$u[[21]]
#' @name imex_multistep_methods
NULL

#' @rdname imex_multistep_methods
#' @export
CNAB2 <- function(autodiff = NULL, extrapolant = "linear") .imex_ms_alg("CNAB2", autodiff, extrapolant)

#' @rdname imex_multistep_methods
#' @export
CNLF2 <- function(autodiff = NULL, extrapolant = "linear") .imex_ms_alg("CNLF2", autodiff, extrapolant)

.imex_ms_alg <- function(name, autodiff, extrapolant) {
  extrapolant <- as.character(extrapolant)
  if (length(extrapolant) != 1L) stop("`extrapolant` must be a single string", call. = FALSE)
  .with_autodiff(autodiff, .new_alg(name, 2L, adaptive = FALSE,
                                    step = function(ctx, t, u, dt, cache) .imex_ms_step(ctx, t, u, dt, cache, name),
                                    init = .imex_ms_init_fn(name), on_event = .imex_ms_on_event, stiff = TRUE,
                                    extrapolant = extrapolant))
}

# initialize!: the pre start derivative f1(u) + f2(u).
.imex_ms_init_fn <- function(name) function(ctx, t, u) {
  if (is.null(ctx$f1) || is.null(ctx$f2))
    stop(name, "() integrates a SplitODEProblem (f1 implicit, f2 explicit)", call. = FALSE)
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  cache$fsal <- ctx$f1(u, ctx$p, t) + ctx$f2(u, ctx$p, t)
  cache$k2 <- NULL
  cache$uprev2 <- u
  cache$nf0 <- 1L
  cache$nf20 <- 1L
  cache
}

# The methods are FSAL: after an event the derivative is evaluated again
# (reset_fsal!, the whole split function counted once in nf), and the next
# step is a startup step (derivative_discontinuity).
.imex_ms_on_event <- function(ctx, t, u, cache) {
  cache$fsal <- ctx$f(u, ctx$p, t)
  list(cache = cache, nf = 1L)
}

.imex_ms_step <- function(ctx, t, u, dt, cache, name) {
  f1 <- ctx$f1; f2 <- ctx$f2; p <- ctx$p
  ctx1 <- ctx
  ctx1$f <- f1
  ctx1$jac <- ctx$jac1
  ctx1$tgrad <- NULL
  startup <- ctx$iter == 1L || isTRUE(ctx$derivative_discontinuity)
  du1 <- f1(u, p, t)
  if (identical(name, "CNAB2")) {
    k1 <- cache$fsal - du1
    tmp <- if (startup) u + dt * k1 else u + dt * (1.5 * k1 - 0.5 * cache$k2)
    gamma <- 1 / 2
    z0 <- dt * du1
    tmp <- tmp + gamma * z0
  } else {
    tmp <- if (startup) u + dt * (cache$fsal - du1) else cache$uprev2 + 2 * dt * (cache$fsal - du1)
    gamma <- 1
    if (!startup) tmp <- tmp + dt * cache$k2
    z0 <- dt * du1
  }
  r <- .nlsolve(cache$nl, ctx1, t, dt, z0, tmp, gamma, 1, t + dt, u, method = "dirk", eta_first = FALSE)
  out <- list(nf = 1L + r$nf, njacs = r$njacs, nsolve = r$nsolve, nw = r$nw, cache = cache)
  if (!r$ok) return(c(out, list(fail = TRUE, u = u)))
  unew <- tmp + gamma * r$z
  if (identical(name, "CNAB2")) {
    cache$k2 <- k1
  } else {
    cache$uprev2 <- u
    cache$k2 <- du1
  }
  fsal0 <- cache$fsal
  cache$fsal <- f1(unew, p, t + dt) + f2(unew, p, t + dt)
  out$nf <- out$nf + 1L
  c(out, list(nf2 = 1L, u = unew, k = list(fsal0, cache$fsal)))
}
