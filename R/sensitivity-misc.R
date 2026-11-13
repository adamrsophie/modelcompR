# Smaller names of SciMLSensitivity.jl v7.119: SensitivityFunction (the
# abstract type of the sensitivity right-hand sides), ODEForwardSensitivityFunction
# (src/forward_sensitivity.jl), AdjointSensitivityIntegrand
# (src/quadrature_adjoint.jl), StochasticTransformedFunction
# (src/sde_tools.jl), SundialsAdjoint, SensitivityAlg and
# supports_functor_params (src/sensitivity_algorithms.jl) and
# ReactantDualTag (src/adjoint_common.jl).

#' Abstract sensitivity function type
#'
#' `SensitivityFunction` is the class name shared by the right-hand sides of
#' the adjoint equations, SciMLSensitivity's abstract type: the right-hand
#' side of an [ODEAdjointProblem()] has the class
#' `ODEInterpolatingAdjointSensitivityFunction`,
#' `ODEBacksolveSensitivityFunction`,
#' `ODEQuadratureAdjointSensitivityFunction` or
#' `ODEGaussAdjointSensitivityFunction` (by method), and that of a
#' [DAEAdjointProblem()] `DAEAdjointSensitivityFunction`, each also a
#' `SensitivityFunction`. As in Julia, [ODEForwardSensitivityFunction()] is
#' not one.
#'
#' @format A character string.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
#' sol <- solve(prob, Tsit5())
#' ap <- ODEAdjointProblem(sol, t = 1, dgdu_discrete = function(u, p, t, i) 1)
#' inherits(ap$f$f, SensitivityFunction)
#' @export
SensitivityFunction <- "SensitivityFunction"

# The class of the adjoint right-hand side of each method.
.sens_function_class <- function(method) {
  switch(method, BacksolveAdjoint = "ODEBacksolveSensitivityFunction",
         QuadratureAdjoint = "ODEQuadratureAdjointSensitivityFunction",
         GaussAdjoint = , GaussKronrodAdjoint = "ODEGaussAdjointSensitivityFunction",
         "ODEInterpolatingAdjointSensitivityFunction")
}

#' Forward sensitivity right-hand side
#'
#' `ODEForwardSensitivityFunction(prob, sensealg, Q)` is the right-hand side
#' `F(z, p, t)` of the forward sensitivity equations of an `ODEProblem`, the
#' function [ODEForwardSensitivityProblem()] integrates: for the state
#' `z = (u, S)` with the sensitivities `S` (`n` x `m`, stored column by
#' column), `F` returns `(f(u, p, t), J S + f_p Q)`, where `Q` (`length(p)`
#' x `m`) holds the derivatives of the parameters with respect to the `m`
#' directions (the identity for the parameter sensitivities). With `tracked =
#' TRUE`, `Q` is part of the state (after `S`) and has zero rate, so events
#' can change it.
#'
#' Derivatives come from the System (symbolic Jacobians) or, for R
#' functions, from `jac_fn` or central differences, as for
#' [ODEForwardSensitivityProblem()].
#'
#' @param prob An `ODEProblem`.
#' @param sensealg [ForwardSensitivity()] or [ForwardDiffSensitivity()].
#' @param Q The parameter seed matrix; the identity when `NULL`.
#' @param tracked Carry `Q` in the state.
#' @return A function `(z, p, t)` of class `ODEForwardSensitivityFunction`
#'   with attributes `numindvar` and `numparams`.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
#' F <- ODEForwardSensitivityFunction(prob)
#' F(c(1, 0), 2, 0)   # (-p u, -p S - u)
#' @export
ODEForwardSensitivityFunction <- function(prob, sensealg = ForwardSensitivity(), Q = NULL, tracked = FALSE) {
  n <- length(prob$u0); np <- length(prob$p)
  Q0 <- as.matrix(Q %||% diag(1, np))
  if (nrow(Q0) != np) stop("Q needs one row per parameter", call. = FALSE)
  m <- ncol(Q0)
  fns <- .sens_fns(prob, sensealg)
  nq <- if (tracked) np * m else 0L
  iu <- seq_len(n); iS <- n + seq_len(n * m); iQ <- n + n * m + seq_len(nq)
  f <- fns$f
  rhs <- function(z, p, t) {
    u <- z[iu]
    S <- matrix(z[iS], n, m)
    Qt <- if (tracked) matrix(z[iQ], np, m) else Q0
    dS <- fns$jac(u, p, t) %*% S
    if (np) dS <- dS + fns$pjac(u, p, t) %*% Qt
    c(f(u, p, t), as.numeric(dS), numeric(nq))
  }
  structure(rhs, class = c("ODEForwardSensitivityFunction", "function"),
            numindvar = n, numparams = np)
}

#' @export
print.ODEForwardSensitivityFunction <- function(x, ...) {
  cat("ODEForwardSensitivityFunction(numindvar = ", attr(x, "numindvar"), ", numparams = ",
      attr(x, "numparams"), ")\n", sep = "")
  invisible(x)
}

# The parameter gradient integrand lambda' f_p + g_p of the quadrature and
# Gauss adjoints at the adjoint state l (forward state from tr$at).
.sens_mu_integrand <- function(fns, tr, gp_c = NULL) {
  function(l, p, s) {
    y <- tr$at(s)
    v <- fns$vjp(l, y, p, s)$p
    if (!is.null(gp_c)) v <- v + as.numeric(gp_c(y, p, s))
    v
  }
}

#' Integrand of the quadrature adjoint
#'
#' `AdjointSensitivityIntegrand(sol, adj_sol, sensealg, dgdp)` is the
#' integrand of the parameter gradient of [QuadratureAdjoint()]: a function
#' of `t` returning `lambda(t)' f_p(u(t), p, t) + dgdp(u(t), p, t)`, with
#' `u` from the dense forward solution `sol` and `lambda` the adjoint state
#' interpolated from the backward solution `adj_sol` (the first `length(u0)`
#' components; for a [DAEAdjointProblem()] solution the multiplier block, and
#' the integrand is then `dgdp - lambda' F_p` of the residual `F`). Its
#' integral from `t0` to `T` is the gradient `dp` of
#' [adjoint_sensitivities()].
#'
#' @param sol The forward solution (dense).
#' @param adj_sol The solution of the backward adjoint problem.
#' @param sensealg An adjoint method, [QuadratureAdjoint()] by default.
#' @param dgdp `NULL` or a function `(u, p, t)` returning the gradient of
#'   the continuous cost with respect to the parameters.
#' @return A function of `t` of class `AdjointSensitivityIntegrand`.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
#' sol <- solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-12)
#' ap <- ODEAdjointProblem(sol, QuadratureAdjoint(), Tsit5(), t = 1,
#'                         dgdu_discrete = function(u, p, t, i) 1)
#' adj <- solve(ap, Tsit5(), reltol = 1e-10, abstol = 1e-12)
#' integrand <- AdjointSensitivityIntegrand(sol, adj)
#' integrand(0.5)   # lambda(0.5) * (-u(0.5))
#' @export
AdjointSensitivityIntegrand <- function(sol, adj_sol, sensealg = QuadratureAdjoint(), dgdp = NULL) {
  prob <- .sol_env(sol)$prob
  dae <- inherits(prob, "DAEProblem")
  fns <- if (dae) .sens_dae_fns(prob, sol) else .sens_fns(prob, sensealg)
  n <- length(prob$u0)
  tr <- if (dae) fns$tr else .sens_trajectory(sol)
  p <- prob$p
  aenv <- .sol_env(adj_sol)
  # the adjoint state of a DAEProblem is (w, lambda)
  il <- if (dae) n + seq_len(n) else seq_len(n)
  vjp_p <- .sens_mu_integrand(fns, tr)
  integrand <- function(t) {
    v <- vjp_p(.interp_states(aenv, t)[1L, il], p, t)
    # for a residual F(du, u, p, t) the integrand is g_p - lambda' F_p
    if (dae) v <- -v
    if (!is.null(dgdp)) v <- v + as.numeric(dgdp(tr$at(t), p, t))
    v
  }
  structure(integrand, class = c("AdjointSensitivityIntegrand", "function"))
}

#' Drift of the Stratonovich form of an Ito SDE
#'
#' `StochasticTransformedFunction(sol, f, g, corfunc_analytical)` wraps the
#' drift `f(u, p, t)` of an Ito SDE solved in `sol` for the stochastic
#' adjoints: it returns `f(u, p, t) - c(u, p, t)` with the correction `c =
#' J' vec(g(u, p, t))`, `J` the Jacobian of `vec(g)` with respect to `u` (as
#' SciMLSensitivity, which computes it by reverse mode differentiation; here
#' central differences), or `c = corfunc_analytical(u, p, t)` when given.
#'
#' @param sol The forward SDE solution (gives the problem).
#' @param f,g Drift and diffusion functions `(u, p, t)`; those of the
#'   problem when `NULL`.
#' @param corfunc_analytical `NULL` or a function `(u, p, t)` returning the
#'   correction.
#' @return A function `(u, p, t)` of class `StochasticTransformedFunction`.
#' @examples
#' prob <- SDEProblem(function(u, p, t) p[1] * u, function(u, p, t) p[2] * u, 1, c(0, 1), c(1.01, 0.87))
#' sol <- solve(prob, EM(), dt = 0.1, seed = 1)
#' tf <- StochasticTransformedFunction(sol)
#' tf(2, c(1.01, 0.87), 0)   # 1.01 * 2 - 0.87^2 * 2
#' @export
StochasticTransformedFunction <- function(sol, f = NULL, g = NULL, corfunc_analytical = NULL) {
  prob <- .sol_env(sol)$prob
  f <- f %||% prob$f$f
  g <- g %||% prob$f$g
  if (inherits(f, "SciMLFunction")) f <- f$f
  if (inherits(g, "SDEFunction")) g <- g$g
  out <- function(u, p, t) {
    ducor <- if (!is.null(corfunc_analytical)) as.numeric(corfunc_analytical(u, p, t)) else {
      gv <- as.numeric(g(u, p, t))
      as.numeric(crossprod(.sens_fd_jac(function(x) as.numeric(g(x, p, t)), u), gv))
    }
    as.numeric(f(u, p, t)) - ducor
  }
  structure(out, class = c("StochasticTransformedFunction", "function"))
}

#' SUNDIALS adjoint
#'
#' `SundialsAdjoint()` names the adjoint of SciMLSensitivity that runs the
#' CVODES (and IDAS) adjoint: the forward solve is checkpointed, the adjoint
#' state and the parameter gradient are integrated backward together, with
#' the forward state interpolated between checkpoints. In R it is that same
#' method, [InterpolatingAdjoint()] with `checkpointing = TRUE`, solved with
#' the backward algorithm given to [adjoint_sensitivities()] (the forward
#' algorithm by default). The SUNDIALS options are kept for the Julia
#' signature.
#'
#' @param chunk_size,autodiff,diff_type Options of the Julia signature.
#' @param autojacvec As for [InterpolatingAdjoint()].
#' @param steps Forward steps between checkpoints of CVODES (unused: the
#'   checkpoints are the saved times).
#' @param interp `"hermite"` or `"polynomial"`, the CVODES interpolation
#'   (checked, unused).
#' @param quad_error_control Whether the gradient quadrature is error
#'   controlled (it always is here).
#' @return A sensitivity algorithm.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -p[1] * u, 1, c(0, 1), 2)
#' sol <- solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-12)
#' adjoint_sensitivities(sol, Tsit5(), t = 1, dgdu_discrete = function(u, p, t, i) 1,
#'                       sensealg = SundialsAdjoint(), abstol = 1e-10, reltol = 1e-10)$dp
#' @export
SundialsAdjoint <- function(chunk_size = 0, autodiff = TRUE, diff_type = "central", autojacvec = NULL,
                            steps = 150, interp = "hermite", quad_error_control = TRUE) {
  interp <- as.character(interp)
  if (!interp %in% c("hermite", "polynomial"))
    stop("SundialsAdjoint `interp` must be `:hermite` or `:polynomial`, got `", interp, "`.", call. = FALSE)
  a <- InterpolatingAdjoint(autojacvec = autojacvec, checkpointing = TRUE)
  a$name <- "SundialsAdjoint"
  a$steps <- steps
  a$interp <- interp
  a$quad_error_control <- quad_error_control
  class(a) <- c("SundialsAdjoint", class(a))
  a
}

#' Deprecated sensitivity algorithm selector
#'
#' `SensitivityAlg()` is the selector of old versions of SciMLSensitivity;
#' as in Julia it only reports (as an error message, without stopping) that
#' the choice mechanism was replaced by the sensitivity algorithms
#' ([sensitivity_algorithms]) and returns `NULL`.
#'
#' @param ... Ignored.
#' @return `NULL`, invisibly.
#' @examples
#' SensitivityAlg()
#' @export
SensitivityAlg <- function(...) {
  message("Error: The SensitivityAlg choice mechanism was completely overhauled. Please consult the local ",
          "sensitivity documentation for more information")
  invisible(NULL)
}

#' Sensitivity algorithms that accept structured parameters
#'
#' `supports_functor_params(sensealg)` is `TRUE` for the sensitivity
#' algorithms that differentiate structured (non vector) parameters
#' directly, [GaussAdjoint()] and [GaussKronrodAdjoint()], as in
#' SciMLSensitivity.
#'
#' @param sensealg A sensitivity algorithm.
#' @return `TRUE` or `FALSE`.
#' @examples
#' supports_functor_params(GaussAdjoint())
#' supports_functor_params(QuadratureAdjoint())
#' @export
supports_functor_params <- function(sensealg) {
  if (!.is_sensealg(sensealg)) stop("`sensealg` must be a sensitivity algorithm", call. = FALSE)
  inherits(sensealg, c("GaussAdjoint", "GaussKronrodAdjoint"))
}

#' Reactant dual tag
#'
#' `ReactantDualTag()` is the tag SciMLSensitivity's Reactant extension uses
#' to tell dual number compilation from floating point compilation. R has no
#' Reactant backend; the object only exists for the Julia name.
#'
#' @return An object of class `ReactantDualTag`.
#' @examples
#' ReactantDualTag()
#' @export
ReactantDualTag <- function() structure(list(), class = "ReactantDualTag")

#' @export
print.ReactantDualTag <- function(x, ...) {
  cat("ReactantDualTag()\n")
  invisible(x)
}
