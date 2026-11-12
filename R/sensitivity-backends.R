# Automatic differentiation choices of SciMLSensitivity.jl
# (src/sensitivity_algorithms.jl: ZygoteVJP, EnzymeVJP, ReverseDiffVJP,
# TrackerVJP, MooncakeVJP, ReactantVJP; the discretize then optimize
# adjoints ZygoteAdjoint, EnzymeAdjoint, ReverseDiffAdjoint, TrackerAdjoint,
# MooncakeAdjoint).

#' Automatic differentiation backends of the sensitivity methods
#'
#' Julia computes the vector Jacobian products of the adjoint methods with
#' an automatic differentiation backend chosen by `autojacvec`. R has no
#' automatic differentiation, so the backend objects `ZygoteVJP()`,
#' `EnzymeVJP()`, `ReverseDiffVJP(compile)`, `TrackerVJP()`, `MooncakeVJP()`
#' and `ReactantVJP()` all select the exact vector Jacobian products: the
#' symbolic derivatives of the System the problem was built from, or the
#' problem's `jac_fn` and central finite differences for problems built from
#' R functions (as `autojacvec = TRUE`). `ReactantVJPConfig()` holds the
#' options of the Reactant backend (unused).
#'
#' The discretize then optimize adjoints `ZygoteAdjoint()`,
#' `EnzymeAdjoint()`, `ReverseDiffAdjoint()`, `TrackerAdjoint()` and
#' `MooncakeAdjoint()` differentiate through the solver in Julia; in R they
#' are the forward sensitivities of the discretized solve,
#' [ForwardDiffSensitivity()], which give the derivatives of the computed
#' solution.
#'
#' @param compile,mode,allow_nothing Options of the Julia signatures.
#' @param ... Options of the Julia signatures.
#' @return A backend or a sensitivity algorithm.
#' @examples
#' InterpolatingAdjoint(autojacvec = ZygoteVJP())
#' ZygoteAdjoint()
#' @name sensitivity_backends
NULL

.vjp_choice <- function(name, ...) structure(list(name = name, ...), class = c(name, "VJPChoice"))

#' @export
print.VJPChoice <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

#' @rdname sensitivity_backends
#' @export
ZygoteVJP <- function(allow_nothing = FALSE) .vjp_choice("ZygoteVJP")
#' @rdname sensitivity_backends
#' @export
EnzymeVJP <- function(...) .vjp_choice("EnzymeVJP")
#' @rdname sensitivity_backends
#' @export
ReverseDiffVJP <- function(compile = FALSE) .vjp_choice("ReverseDiffVJP", compile = compile)
#' @rdname sensitivity_backends
#' @export
TrackerVJP <- function(allow_nothing = FALSE) .vjp_choice("TrackerVJP")
#' @rdname sensitivity_backends
#' @export
MooncakeVJP <- function(...) .vjp_choice("MooncakeVJP")
#' @rdname sensitivity_backends
#' @export
ReactantVJP <- function(...) .vjp_choice("ReactantVJP")
#' @rdname sensitivity_backends
#' @export
ReactantVJPConfig <- function(...) structure(list(...), class = "ReactantVJPConfig")

.dto_adjoint <- function(name) {
  a <- ForwardDiffSensitivity()
  a$name <- name
  class(a) <- c(name, class(a))
  a
}

#' @rdname sensitivity_backends
#' @export
ZygoteAdjoint <- function() .dto_adjoint("ZygoteAdjoint")
#' @rdname sensitivity_backends
#' @export
EnzymeAdjoint <- function(mode = NULL) .dto_adjoint("EnzymeAdjoint")
#' @rdname sensitivity_backends
#' @export
ReverseDiffAdjoint <- function(compile = FALSE) .dto_adjoint("ReverseDiffAdjoint")
#' @rdname sensitivity_backends
#' @export
TrackerAdjoint <- function() .dto_adjoint("TrackerAdjoint")
#' @rdname sensitivity_backends
#' @export
MooncakeAdjoint <- function() .dto_adjoint("MooncakeAdjoint")
