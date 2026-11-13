# Implicit integrating factor methods of StochasticDiffEqIIF (src/
# algorithms.jl, src/alg_utils.jl, src/caches/iif_caches.jl,
# src/perform_step/iif.jl) for a SplitSDEProblem whose f1(u, p, t) gives
# the matrix A of the linear part A u, with the nonlinear part f2 treated
# implicitly by SimpleTrustRegion (NLSOLVEJL_SETUP of StochasticDiffEqCore
# src/misc_utils.jl). IIF1Mil adds the Milstein correction of the in place
# step, which upstream calls with the old f1(t, u, du) signature; the out
# of place step stops with "Milstein correction does not work.".

#' Implicit integrating factor methods for semilinear SDEs
#'
#' For `du = (A(u, p, t) u + f2(u, p, t)) dt + g(u, p, t) dW` given as a
#' [SplitSDEProblem()] whose `f1(u, p, t)` returns the matrix `A` (a number
#' for scalar problems), `IIF1M()` solves `u = exp(A dt) (uprev + g dW) + dt
#' f2(u, p, t + dt)` and `IIF2M()` solves `u = exp(A dt) (uprev + dt / 2
#' f2(uprev, p, t) + g dW) + dt / 2 f2(u, p, t + dt)`, both of strong order
#' 1/2. `IIF1Mil()` adds the Milstein correction `sum_j (M g(uprev +
#' sqrt(dt) M g_j) - M g(uprev)) / sqrt(dt) I_j` with `M = exp(A dt)` and
#' the iterated integrals `I = (dW dW' - dt I) / 2` of commutative noise, for
#' strong order 1. `A` is evaluated at the start of the step. The implicit
#' equation is solved by the trust region method of [SimpleTrustRegion()]
#' from the current state; `nlsolve` is a [NLSOLVEJL_SETUP()] object, kept
#' for the signature. These methods take a fixed step size `dt`.
#'
#' `IIFNLSolveFunc(f)` wraps the residual function `f(u)` of such an implicit
#' equation, as the `Val{:init}` call of `NLSOLVEJL_SETUP()` does upstream,
#' and `NLSOLVEJL_SETUP()` solves `f(u) = 0` from `u0` when called as
#' `setup(wrapped, u0)`.
#'
#' @param nlsolve The nonlinear solver setup.
#' @param autodiff The differentiation choice of the inner solve; R uses
#'   finite differences.
#' @param f A residual function of `u`.
#' @return An SDE algorithm, an `NLSOLVEJL_SETUP` object (a function of a
#'   wrapped residual and a start value) or an `IIFNLSolveFunc` object.
#' @examples
#' prob <- SplitSDEProblem(function(u, p, t) -1, function(u, p, t) 0.5 * u,
#'                         function(u, p, t) 0.2 * u, 1, c(0, 1))
#' sol <- solve(prob, IIF1M(), dt = 0.1, seed = 1)
#' sol$u[[11]]
#' setup <- NLSOLVEJL_SETUP()
#' setup(IIFNLSolveFunc(function(u) u^2 - 2), 1)
#' @name iif_methods
NULL

#' @rdname iif_methods
#' @export
NLSOLVEJL_SETUP <- function(autodiff = AutoForwardDiff()) {
  structure(function(f, u0) {
    if (inherits(f, "IIFNLSolveFunc")) f <- f$f
    .iif_root(f, u0)
  }, autodiff = autodiff, class = c("NLSOLVEJL_SETUP", "function"))
}

#' @rdname iif_methods
#' @export
IIFNLSolveFunc <- function(f) structure(list(f = f), class = "IIFNLSolveFunc")

# SimpleTrustRegion from u0, run until the residual is at rounding level.
.iif_root <- function(F, u0) {
  u0 <- .as_num(u0)
  tol <- 8 * .Machine$double.eps * max(1, abs(u0))
  .nl_trust_region(function(x) .as_num(F(x)), NULL, u0, tol, 100L)$u
}

#' @rdname iif_methods
#' @export
IIF1M <- function(nlsolve = NLSOLVEJL_SETUP()) .iif_alg("IIF1M", 1 / 2, nlsolve)

#' @rdname iif_methods
#' @export
IIF2M <- function(nlsolve = NLSOLVEJL_SETUP()) .iif_alg("IIF2M", 1 / 2, nlsolve)

#' @rdname iif_methods
#' @export
IIF1Mil <- function(nlsolve = NLSOLVEJL_SETUP()) .iif_alg("IIF1Mil", 1, nlsolve)

.iif_alg <- function(name, order, nlsolve) {
  alg <- .new_sde_alg(name, order, FALSE, NULL, nlsolve = nlsolve)
  alg$bind <- function(prob, opts) .iif_bind(alg, prob)
  alg
}

.iif_bind <- function(alg, prob) {
  f1 <- prob$f$f1; f2 <- prob$f$f2
  if (is.null(f1) || is.null(f2))
    .stopf("%s() needs a SplitSDEProblem whose f1(u, p, t) gives the matrix of the linear part", alg$name)
  name <- alg$name
  a <- alg
  a$bind <- NULL
  a$init <- function(ctx, t, u) list()
  a$step <- function(ctx, g, u, t, dt, dW, type, delta, cache, ...) {
    p <- ctx$p
    A <- f1(u, p, t)
    scalar <- length(A) == 1L && length(u) == 1L
    M <- if (scalar) exp(.as_num(A) * dt) else .expm(as.matrix(A) * dt)
    Mv <- function(v) if (scalar) M * v else .as_num(M %*% v)
    G <- g(u, p, t)
    gdw <- .gdw(G, dW, type)
    tmp <- switch(name,
      IIF1M = Mv(u + gdw),
      IIF2M = Mv(u + (dt / 2) * f2(u, p, t) + gdw),
      IIF1Mil = Mv(u) + Mv(gdw) + .iif_milstein(G, g, u, p, t, dt, dW, type, Mv))
    c <- if (name == "IIF2M") dt / 2 else dt
    unew <- .iif_root(function(x) x - tmp - c * f2(x, p, t + dt), u)
    list(u = unew, err = NULL, k = NULL, cache = cache, nf = 2L)
  }
  a
}

# The Milstein correction of the in place IIF1Mil step: Dg_j = (M g(uprev +
# sqrt(dt) (M G)_j) - M G) / sqrt(dt) applied to the column j of I.
.iif_milstein <- function(G, g, u, p, t, dt, dW, type, Mv) {
  n <- length(u)
  m <- length(dW)
  Gm <- switch(type, diagonal = diag(.as_num(G), n), scalar = matrix(.as_num(G), n, 1L), as.matrix(G))
  MG <- apply(Gm, 2L, Mv)
  if (!is.matrix(MG)) MG <- matrix(MG, n)
  I <- 0.5 * outer(dW, dW)
  diag(I) <- diag(I) - 0.5 * dt
  sqdt <- .sqdt(dt)
  mil <- numeric(n)
  for (j in seq_len(m)) {
    Kj <- u + sqdt * MG[, j]
    Gj <- g(Kj, p, t)
    Gjm <- switch(type, diagonal = diag(.as_num(Gj), n), scalar = matrix(.as_num(Gj), n, 1L), as.matrix(Gj))
    MGj <- apply(Gjm, 2L, Mv)
    if (!is.matrix(MGj)) MGj <- matrix(MGj, n)
    mil <- mil + .as_num(((MGj - MG) / sqdt) %*% I[, j])
  }
  mil
}
