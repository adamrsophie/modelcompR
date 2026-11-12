# The automatic differentiation choices of ADTypes.jl (src/dense.jl,
# src/legacy.jl) reexported by NonlinearSolve, and NonlinearSolveBase's
# pickchunksize (ForwardDiff.pickchunksize, src/prelude.jl).
#
# R has no automatic differentiation. Every backend that differentiates
# exactly (forward, reverse or symbolic mode) is a subclass of
# AutoForwardDiff here, so the solvers use the problem's Jacobian when it has
# one and central differences otherwise, as for AutoForwardDiff(). The
# finite difference backends of other packages are subclasses of
# AutoFiniteDiff with central differences.

#' Automatic differentiation choices
#'
#' The backends of ADTypes.jl, for the `autodiff` arguments of the solvers.
#' R has no automatic differentiation: every backend that differentiates
#' exactly (`AutoZygote()`, `AutoEnzyme()`, `AutoReverseDiff()`,
#' `AutoSymbolics()`, ...) uses the problem's Jacobian when it has one and
#' central differences otherwise, as [AutoForwardDiff()] does, and
#' `AutoFiniteDifferences()` uses central differences (the bounded solvers
#' use their bound aware differences only for `AutoFiniteDiff()` itself).
#' `NoAutoDiff()` selects no differentiation: the bounded solvers
#' ([TrustRegionReflective()] and the others) signal
#' `NoAutoDiffSelectedError()` when the problem has no Jacobian; the other
#' solvers fall back to central differences.
#'
#' The deprecated constructors of ADTypes return their replacements:
#' `AutoSparseForwardDiff()` is `AutoSparse(AutoForwardDiff())` and so on,
#' `AutoModelingToolkit()` is `AutoSymbolics()` (sparse when either flag is
#' set) and `AutoTapir()` is kept for `AutoMooncake()`.
#'
#' `pickchunksize(x)` is the chunk size ForwardDiff would use for an input of
#' length `x` (a single whole number) or of the length of the vector `x`:
#' the length itself up to 12, otherwise the length divided evenly into
#' chunks of at most 12.
#'
#' @param ruleconfig,mode,function_annotation,fdm,descriptor,config,chunksize,tag,compile,safe_mode,order,simd,jet
#'   Options of the Julia constructors, stored in the object.
#' @param obj_sparse,cons_sparse Sparsity flags of `AutoModelingToolkit()`.
#' @param fdtype Difference type of `AutoSparseFiniteDiff()`.
#' @param msg The message of the condition.
#' @param x An input length, or a vector.
#' @param threshold The largest chunk size.
#' @return An ADType object; `NoAutoDiffSelectedError()` a condition;
#'   `pickchunksize()` an integer.
#' @examples
#' f <- function(u, p) c(u[1] * u[1] - 2, u[2] - 1)
#' prob <- NonlinearProblem(f, u0 = c(1, 0), lb = c(0, 0), ub = c(2, 2))
#' solve(prob, TrustRegionReflective(autodiff = AutoZygote()))$u
#' AutoSparseForwardDiff()
#' pickchunksize(30)
#' tryCatch(stop(NoAutoDiffSelectedError()), NoAutoDiffSelectedError = conditionMessage)
#' @name adtypes
NULL

.ad_exact <- function(name, mode, ...) {
  a <- .adtype(name, mode = mode, ...)
  class(a) <- c(name, "AutoForwardDiff", "ADType")
  a
}

#' @rdname adtypes
#' @export
AutoChainRules <- function(ruleconfig = NULL) .ad_exact("AutoChainRules", "ForwardOrReverseMode", ruleconfig = ruleconfig)

#' @rdname adtypes
#' @export
AutoDiffractor <- function() .ad_exact("AutoDiffractor", "ForwardOrReverseMode")

#' @rdname adtypes
#' @export
AutoEnzyme <- function(mode = NULL, function_annotation = NULL)
  .ad_exact("AutoEnzyme", "ForwardOrReverseMode", enzyme_mode = mode, function_annotation = function_annotation)

#' @rdname adtypes
#' @export
AutoReactant <- function(mode = NULL) {
  if (is.null(mode)) mode <- AutoEnzyme()
  if (!inherits(mode, "AutoEnzyme")) stop("`mode` must be AutoEnzyme() or NULL", call. = FALSE)
  .ad_exact("AutoReactant", "ForwardOrReverseMode", enzyme = mode)
}

#' @rdname adtypes
#' @export
AutoFastDifferentiation <- function() .ad_exact("AutoFastDifferentiation", "SymbolicMode")

#' @rdname adtypes
#' @export
AutoFiniteDifferences <- function(fdm = NULL) {
  a <- .adtype("AutoFiniteDifferences", mode = "ForwardMode", fdtype = "central", fdm = fdm)
  class(a) <- c("AutoFiniteDifferences", "AutoFiniteDiff", "ADType")
  a
}

#' @rdname adtypes
#' @export
AutoGTPSA <- function(descriptor = NULL) .ad_exact("AutoGTPSA", "ForwardMode", descriptor = descriptor)

#' @rdname adtypes
#' @export
AutoHyperHessians <- function(chunksize = NULL, simd = FALSE, jet = FALSE) {
  if (!is.null(chunksize) && !(chunksize > 0)) stop("chunksize must be positive, got ", chunksize, call. = FALSE)
  .ad_exact("AutoHyperHessians", "ForwardMode", chunksize = chunksize, simd = simd, jet = jet)
}

#' @rdname adtypes
#' @export
AutoMooncake <- function(config = NULL) .ad_exact("AutoMooncake", "ReverseMode", config = config)

#' @rdname adtypes
#' @export
AutoMooncakeForward <- function(config = NULL) .ad_exact("AutoMooncakeForward", "ForwardMode", config = config)

#' @rdname adtypes
#' @export
AutoPolyesterForwardDiff <- function(chunksize = NULL, tag = NULL)
  .ad_exact("AutoPolyesterForwardDiff", "ForwardMode", chunksize = chunksize, tag = tag)

#' @rdname adtypes
#' @export
AutoReverseDiff <- function(compile = FALSE) .ad_exact("AutoReverseDiff", "ReverseMode", compile = isTRUE(compile))

#' @rdname adtypes
#' @export
AutoSymbolics <- function() .ad_exact("AutoSymbolics", "SymbolicMode")

#' @rdname adtypes
#' @export
AutoTapir <- function(safe_mode = TRUE) .ad_exact("AutoTapir", "ReverseMode", safe_mode = isTRUE(safe_mode))

#' @rdname adtypes
#' @export
AutoTaylorDiff <- function(order = 1L) .ad_exact("AutoTaylorDiff", "ForwardMode", order = as.integer(order))

#' @rdname adtypes
#' @export
AutoTracker <- function() .ad_exact("AutoTracker", "ReverseMode")

#' @rdname adtypes
#' @export
AutoZygote <- function() .ad_exact("AutoZygote", "ReverseMode")

#' @rdname adtypes
#' @export
AutoModelingToolkit <- function(obj_sparse = FALSE, cons_sparse = FALSE)
  if (isTRUE(obj_sparse) || isTRUE(cons_sparse)) AutoSparse(AutoSymbolics()) else AutoSymbolics()

#' @rdname adtypes
#' @export
AutoSparseFastDifferentiation <- function() AutoSparse(AutoFastDifferentiation())

#' @rdname adtypes
#' @export
AutoSparseFiniteDiff <- function(fdtype = c("forward", "central", "complex")) AutoSparse(AutoFiniteDiff(match.arg(fdtype)))

#' @rdname adtypes
#' @export
AutoSparseForwardDiff <- function(chunksize = NULL, tag = NULL) AutoSparse(AutoForwardDiff(chunksize, tag))

#' @rdname adtypes
#' @export
AutoSparsePolyesterForwardDiff <- function(chunksize = NULL, tag = NULL)
  AutoSparse(AutoPolyesterForwardDiff(chunksize, tag))

#' @rdname adtypes
#' @export
AutoSparseReverseDiff <- function(compile = FALSE) AutoSparse(AutoReverseDiff(compile))

#' @rdname adtypes
#' @export
AutoSparseZygote <- function() AutoSparse(AutoZygote())

#' @rdname adtypes
#' @export
NoAutoDiff <- function() .adtype("NoAutoDiff")

#' @rdname adtypes
#' @export
NoAutoDiffSelectedError <- function(msg = "Automatic differentiation can not be used with NoAutoDiff()")
  structure(class = c("NoAutoDiffSelectedError", "error", "condition"), list(message = msg, call = NULL))

#' @rdname adtypes
#' @export
pickchunksize <- function(x, threshold = 12L) {
  n <- if (is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)) as.integer(x) else length(x)
  if (n <= threshold) return(as.integer(n))
  nchunks <- ceiling(n / threshold)
  as.integer(ceiling(n / nchunks))
}
