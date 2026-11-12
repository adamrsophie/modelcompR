# Sensitivity algorithms (sensealg objects), after SciMLSensitivity.jl v7.119
# src/sensitivity_algorithms.jl (ForwardSensitivity, ForwardDiffSensitivity,
# BacksolveAdjoint, InterpolatingAdjoint, QuadratureAdjoint, GaussAdjoint,
# GaussKronrodAdjoint, SteadyStateAdjoint, ForwardLSS, AdjointLSS,
# TimeDilation, CosWindowing, Cos2Windowing, NILSS, NILSAS,
# ForwardDiffOverAdjoint).
#
# Julia chooses how vector Jacobian products are computed with an automatic
# differentiation backend (`autojacvec = ZygoteVJP()`, `EnzymeVJP()`, ...).
# R has no automatic differentiation: `autojacvec = TRUE` (or NULL, the
# default) uses the symbolic derivatives of the System the problem was built
# from (Jacobians and the gradient of lambda' f, compiled to R code), and
# `autojacvec = FALSE`, or a problem built from a plain R function, uses
# central finite differences (or the user Jacobian `jac_fn` when given).

.sensealg <- function(name, kind, ...) {
  structure(list(name = name, ...), class = c(name, kind, "SensitivityAlgorithm"))
}

#' @export
print.SensitivityAlgorithm <- function(x, ...) {
  flds <- x[setdiff(names(x), c("name", "g", "adjalg"))]
  flds <- flds[!vapply(flds, is.null, TRUE)]
  args <- vapply(names(flds), function(n) paste0(n, " = ", format(flds[[n]])[1L]), "")
  cat(x$name, "(", paste(args, collapse = ", "), ")\n", sep = "")
  invisible(x)
}

#' Sensitivity algorithms
#'
#' Choices of `sensealg` for [solve_gradient()], [adjoint_sensitivities()]
#' and [ODEForwardSensitivityProblem()], named as in SciMLSensitivity.jl.
#'
#' Forward methods integrate the sensitivity equations
#' `dS/dt = J S + f_p` together with the state. `ForwardSensitivity()`
#' differentiates with respect to the parameters; `ForwardDiffSensitivity()`
#' (Julia's dual number method) differentiates with respect to the initial
#' values and the parameters and is the forward method that supports
#' events.
#'
#' Adjoint methods solve the adjoint equation backward in time.
#' `InterpolatingAdjoint()` integrates the adjoint state and the parameter
#' gradient together, reading the forward state from the dense forward
#' solution, or with `checkpointing = TRUE` from short re-solves between the
#' saved points. `BacksolveAdjoint()` integrates the state backward along
#' with the adjoint (with `checkpointing = TRUE` the state is reset to the
#' saved forward values). `QuadratureAdjoint()` solves the adjoint state
#' alone and integrates the parameter gradient afterwards by adaptive
#' Gauss-Kronrod quadrature; `GaussAdjoint()` and `GaussKronrodAdjoint()`
#' integrate it step by step with a Gauss-Legendre or Gauss-Kronrod rule on
#' each step of the adjoint solve. `SteadyStateAdjoint()` differentiates
#' steady states and roots of nonlinear problems by the implicit function
#' theorem.
#'
#' `ForwardLSS()`, `AdjointLSS()`, `NILSS()` and `NILSAS()` are the
#' shadowing methods for long time averages of chaotic systems (see
#' [ForwardLSSProblem()]). `ForwardDiffOverAdjoint()` names the second
#' order method of [second_order_sensitivities()].
#'
#' @param autojacvec `NULL` or `TRUE` use the symbolic derivatives of the
#'   System the problem was built from; `FALSE` uses finite differences.
#'   Problems built from R functions always use finite differences (or
#'   their `jac_fn`). The automatic differentiation backends of Julia
#'   (`ZygoteVJP()`, `EnzymeVJP()`, `ReverseDiffVJP()`, ...) have no R
#'   counterpart.
#' @param checkpointing For `InterpolatingAdjoint()`, `GaussAdjoint()` and
#'   `GaussKronrodAdjoint()`: re-solve the forward problem between saved
#'   points instead of reading a dense forward solution. For
#'   `BacksolveAdjoint()`: reset the backward state to the saved forward
#'   values.
#' @param abstol,reltol Tolerances of the quadrature of
#'   `QuadratureAdjoint()`.
#' @param convert_tspan Unused; kept for the Julia signature.
#' @return A sensitivity algorithm object.
#' @examples
#' InterpolatingAdjoint(checkpointing = TRUE)
#' GaussAdjoint()
#' ForwardDiffSensitivity()
#' @name sensitivity_algorithms
NULL

#' @rdname sensitivity_algorithms
#' @export
ForwardSensitivity <- function(autojacvec = NULL) {
  .sensealg("ForwardSensitivity", "ForwardSensitivityAlgorithm", autojacvec = autojacvec)
}

#' @rdname sensitivity_algorithms
#' @export
ForwardDiffSensitivity <- function(autojacvec = NULL, convert_tspan = NULL) {
  .sensealg("ForwardDiffSensitivity", "ForwardSensitivityAlgorithm", autojacvec = autojacvec,
            convert_tspan = convert_tspan)
}

#' @rdname sensitivity_algorithms
#' @export
InterpolatingAdjoint <- function(autojacvec = NULL, checkpointing = FALSE) {
  .sensealg("InterpolatingAdjoint", "AdjointSensitivityAlgorithm", autojacvec = autojacvec,
            checkpointing = isTRUE(checkpointing))
}

#' @rdname sensitivity_algorithms
#' @export
BacksolveAdjoint <- function(autojacvec = NULL, checkpointing = TRUE) {
  .sensealg("BacksolveAdjoint", "AdjointSensitivityAlgorithm", autojacvec = autojacvec,
            checkpointing = isTRUE(checkpointing))
}

#' @rdname sensitivity_algorithms
#' @export
QuadratureAdjoint <- function(autojacvec = NULL, abstol = 1e-6, reltol = 1e-3) {
  .sensealg("QuadratureAdjoint", "AdjointSensitivityAlgorithm", autojacvec = autojacvec,
            abstol = abstol, reltol = reltol, checkpointing = FALSE)
}

#' @rdname sensitivity_algorithms
#' @export
GaussAdjoint <- function(autojacvec = NULL, checkpointing = FALSE) {
  .sensealg("GaussAdjoint", "AdjointSensitivityAlgorithm", autojacvec = autojacvec,
            checkpointing = isTRUE(checkpointing))
}

#' @rdname sensitivity_algorithms
#' @export
GaussKronrodAdjoint <- function(autojacvec = NULL, checkpointing = FALSE) {
  .sensealg("GaussKronrodAdjoint", "AdjointSensitivityAlgorithm", autojacvec = autojacvec,
            checkpointing = isTRUE(checkpointing))
}

#' @rdname sensitivity_algorithms
#' @export
SteadyStateAdjoint <- function(autojacvec = NULL) {
  .sensealg("SteadyStateAdjoint", "AdjointSensitivityAlgorithm", autojacvec = autojacvec)
}

#' @rdname sensitivity_algorithms
#' @param adjalg The first order adjoint method used inside.
#' @export
ForwardDiffOverAdjoint <- function(adjalg = InterpolatingAdjoint()) {
  .sensealg("ForwardDiffOverAdjoint", "SecondOrderSensitivityAlgorithm", adjalg = adjalg)
}

#' Shadowing sensitivity algorithms for chaotic systems
#'
#' Least squares shadowing computes derivatives of long time averages
#' `<g> = (1/T) int g(u, p, t) dt` of ergodic chaotic systems, whose
#' ordinary sensitivities grow without bound. `ForwardLSS()` and
#' `AdjointLSS()` solve the discretized least squares problem on the saved
#' points of a solution; `NILSS()` (forward) and `NILSAS()` (adjoint)
#' restrict the computation to the unstable subspace, integrating tangent or
#' adjoint solutions over `nseg` segments of `nstep` saved points each.
#' `TimeDilation()`, `CosWindowing()` and `Cos2Windowing()` choose how
#' LSS treats the time average.
#'
#' @param LSSregularizer `TimeDilation(alpha, t0skip, t1skip)` (the
#'   default, weight `alpha` of the time dilation, times skipped at both
#'   ends), `CosWindowing()` or `Cos2Windowing()`.
#' @param g The instantaneous objective `g(u, p, t)`.
#' @param nseg,nstep Number of segments and of saved points per segment.
#' @param nus Dimension of the unstable subspace for `NILSS()` (default
#'   `length(u0) - 1`).
#' @param M Number of homogeneous adjoint solutions for `NILSAS()` (at least
#'   the number of positive Lyapunov exponents plus one).
#' @param seed Seed of the random initial tangent or adjoint directions.
#' @param adjoint_sensealg For `NILSAS()`: the continuous adjoint method of
#'   the segment solves. [BacksolveAdjoint()] integrates the state backward
#'   with the adjoint from the state saved at each segment end (Julia's
#'   default); the other methods (the default `NULL` is
#'   [InterpolatingAdjoint()]) read the state from the dense forward solution.
#' @param autojacvec As in [sensitivity_algorithms].
#' @param alpha,t0skip,t1skip Time dilation weight and skipped times.
#' @return An algorithm object.
#' @examples
#' ForwardLSS(g = function(u, p, t) u[3])
#' NILSS(10, 101, g = function(u, p, t) u[3])
#' @name shadowing_algorithms
NULL

#' @rdname shadowing_algorithms
#' @export
ForwardLSS <- function(LSSregularizer = TimeDilation(10), g = NULL, autojacvec = NULL) {
  .sensealg("ForwardLSS", "ShadowingSensitivityAlgorithm", LSSregularizer = LSSregularizer, g = g,
            autojacvec = autojacvec)
}

#' @rdname shadowing_algorithms
#' @export
AdjointLSS <- function(LSSregularizer = TimeDilation(10), g = NULL, autojacvec = NULL) {
  .sensealg("AdjointLSS", "ShadowingSensitivityAlgorithm", LSSregularizer = LSSregularizer, g = g,
            autojacvec = autojacvec)
}

#' @rdname shadowing_algorithms
#' @export
TimeDilation <- function(alpha, t0skip = 0, t1skip = 0) {
  structure(list(alpha = alpha, t0skip = t0skip, t1skip = t1skip),
            class = c("TimeDilation", "LSSRegularizer"))
}

#' @rdname shadowing_algorithms
#' @export
CosWindowing <- function() structure(list(power = 1), class = c("CosWindowing", "LSSRegularizer"))

#' @rdname shadowing_algorithms
#' @export
Cos2Windowing <- function() structure(list(power = 2), class = c("Cos2Windowing", "LSSRegularizer"))

#' @export
format.LSSRegularizer <- function(x, ...) {
  if (inherits(x, "TimeDilation")) sprintf("TimeDilation(%g, %g, %g)", x$alpha, x$t0skip, x$t1skip)
  else paste0(class(x)[1L], "()")
}

#' @rdname shadowing_algorithms
#' @export
NILSS <- function(nseg, nstep, nus = NULL, g = NULL, seed = NULL, autojacvec = NULL) {
  .sensealg("NILSS", "ShadowingSensitivityAlgorithm", nseg = as.integer(nseg), nstep = as.integer(nstep),
            nus = nus, g = g, seed = seed, autojacvec = autojacvec)
}

#' @rdname shadowing_algorithms
#' @export
NILSAS <- function(nseg, nstep, M = NULL, g = NULL, seed = NULL, autojacvec = NULL,
                   adjoint_sensealg = NULL) {
  if (nstep <= 1) stop("At least the start and the end point of each segment must be stored. ",
                       "Please use `nstep >= 2`.", call. = FALSE)
  if (!is.null(adjoint_sensealg) && !.is_sensealg(adjoint_sensealg, "AdjointSensitivityAlgorithm"))
    stop("`adjoint_sensealg` must be a continuous adjoint method such as BacksolveAdjoint()", call. = FALSE)
  .sensealg("NILSAS", "ShadowingSensitivityAlgorithm", nseg = as.integer(nseg), nstep = as.integer(nstep),
            M = M, g = g, seed = seed, autojacvec = autojacvec, adjoint_sensealg = adjoint_sensealg)
}

.is_sensealg <- function(x, kind = "SensitivityAlgorithm") inherits(x, kind)
