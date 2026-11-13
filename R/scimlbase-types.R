# Small SciMLBase.jl and OrdinaryDiffEq.jl objects: the ReturnCode enum
# (src/retcodes.jl), NullParameters, the ODESolution constructor, the
# automatic differentiation choices of ADTypes.jl reexported by
# OrdinaryDiffEq (AutoForwardDiff, AutoFiniteDiff, AutoSparse) and the
# threading options of the explicit Runge-Kutta methods (Sequential,
# BaseThreads, PolyesterThreads).

#' SciMLBase return codes and small types
#'
#' `ReturnCode` lists the return codes of SciMLBase (`ReturnCode$Success`,
#' `ReturnCode$MaxIters`, ...); solutions carry them as strings in
#' `sol$retcode`, and [successful_retcode()] tells which ones mean success.
#' `NullParameters()` is the empty parameter object of problems without
#' parameters. `ODESolution(t, u, prob, alg)` builds a solution from saved
#' times and states (a matrix with one row per time, or a list of states),
#' interpolated linearly.
#'
#' The Jacobians of the implicit solvers come from the problem's symbolic
#' Jacobian (from a `System`, or `jac_fn`) or from differentiation.
#' `AutoForwardDiff()` uses the exact Jacobian when the problem has one,
#' else dual numbers as ForwardDiff does (central differences when the
#' function cannot run on them), `AutoFiniteDiff(fdtype)` finite
#' differences of the given type,
#' and `AutoSparse(ad)` the method of `ad` (Jacobians are dense in R).
#' `Sequential()`, `BaseThreads()` and `PolyesterThreads()` are the threading
#' options of the Runge-Kutta stages; R evaluates the stages in turn.
#'
#' @param t,u Saved times and states.
#' @param prob The problem the solution belongs to.
#' @param alg The algorithm (any object with a `name`).
#' @param retcode The return code.
#' @param fdtype `"forward"`, `"central"` or `"complex"`.
#' @param chunksize,tag Options of the Julia signature.
#' @param ad An automatic differentiation choice.
#' @param sparsity_detector,coloring_algorithm Options of the Julia
#'   signature.
#' @return An object of the named kind.
#' @examples
#' ReturnCode$Success
#' successful_retcode(ReturnCode$Terminated)
#' sol <- ODESolution(c(0, 1), rbind(0, 1), ODEProblem(function(u, p, t) 1, 0, c(0, 1)))
#' sol(0.5)
#' @name scimlbase_types
NULL

#' @rdname scimlbase_types
#' @export
ReturnCode <- local({
  codes <- c("Default", "Success", "Terminated", "DtNaN", "MaxIters", "MaxNumSub", "DtLessThanMin", "Unstable",
             "InitialFailure", "ConvergenceFailure", "Failure", "ExactSolutionLeft", "ExactSolutionRight",
             "FloatingPointLimit", "Infeasible", "MaxTime", "InternalLineSearchFailed",
             "ShrinkThresholdExceeded", "Stalled", "StalledSuccess", "InternalLinearSolveFailed",
             "APosterioriSafetyFailure")
  as.list(stats::setNames(codes, codes))
})

#' @rdname scimlbase_types
#' @export
NullParameters <- function() structure(numeric(), class = c("NullParameters", "numeric"))

#' @rdname scimlbase_types
#' @export
ODESolution <- function(t, u, prob, alg = list(name = "ODESolution"), retcode = "Success") {
  U <- if (is.list(u)) do.call(rbind, lapply(u, as.numeric)) else as.matrix(u)
  t <- as.numeric(t)
  if (nrow(U) != length(t)) stop("give one state per time", call. = FALSE)
  lin <- function(theta, dt, y0, y1, k) (1 - theta) * y0 + theta * y1
  steps <- lapply(seq_len(max(0L, length(t) - 1L)), function(i)
    list(t0 = t[i], dt = t[i + 1L] - t[i], tend = t[i + 1L], y0 = U[i, ], y1 = U[i + 1L, ], k = NULL, interp = lin))
  colnames(U) <- names(prob$u0)
  .new_ode_solution(t = t, u = U, prob = prob, alg = alg, retcode = retcode, steps = steps, p = prob$p, dense = TRUE)
}

.adtype <- function(name, ...) structure(list(name = name, ...), class = c(name, "ADType"))

#' @export
print.ADType <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

#' @rdname scimlbase_types
#' @export
AutoForwardDiff <- function(chunksize = NULL, tag = NULL) .adtype("AutoForwardDiff")

#' @rdname scimlbase_types
#' @export
AutoFiniteDiff <- function(fdtype = c("forward", "central", "complex")) .adtype("AutoFiniteDiff", fdtype = match.arg(fdtype))

#' @rdname scimlbase_types
#' @export
AutoSparse <- function(ad, sparsity_detector = NULL, coloring_algorithm = NULL) .adtype("AutoSparse", dense_ad = ad)

#' @rdname scimlbase_types
#' @export
Sequential <- function() structure(list(name = "Sequential"), class = "ThreadingOption")

#' @rdname scimlbase_types
#' @export
BaseThreads <- function() structure(list(name = "BaseThreads"), class = "ThreadingOption")

#' @rdname scimlbase_types
#' @export
PolyesterThreads <- function() structure(list(name = "PolyesterThreads"), class = "ThreadingOption")
