# The OptimizationBase.jl interface reexported by Optimization.jl
# (src/OptimizationBase.jl, cache.jl, solve.jl and SciMLBase's optimizer
# traits): the objective sense, the default callback, optimizer traits,
# OptimizationStats, the cache of init/solve!, verbosity and the error
# types.

#' Optimization interface
#'
#' `OptimizationProblem(..., sense = MaxSense)` maximizes instead of
#' minimizing (`MinSense`, the default; `ObjSense` lists both); the solution's
#' `objective` is the value of the original function. `solve(prob, alg,
#' callback = cb)` calls `cb(state, loss)` after every evaluation of the
#' objective, and stops when it returns `TRUE` (`DEFAULT_CALLBACK` never
#' stops; `DEFAULT_DATA` is the data iterator of Julia's data driven
#' optimizers). Solutions carry `sol$stats`, an `OptimizationStats` with the
#' numbers of iterations and evaluations and the time.
#'
#' `init(prob, alg, ...)` returns an `OptimizationCache`, which
#' `solve_integrator()` (Julia's `solve!`) solves; `reinit(cache, u0, p)`
#' changes its initial point or parameters.
#'
#' The traits `requiresgradient`, `requireshessian`, `requiresconsjac`,
#' `requiresconshess`, `requiresbounds`, `requiresconstraints`,
#' `allowsbounds`, `allowsconstraints`, `allowscallback` describe an
#' optimizer, and `lag_hess_structure` is the sparsity pattern of the
#' Hessian of the Lagrangian of an `OptimizationFunction` (dense here).
#' `NoAD()` builds no derivatives (the user's or none). `OptimizationVerbosity()`
#' holds verbosity settings. A solver refuses problems it cannot handle with
#' an `IncompatibleOptimizerError` condition; `OptimizerMissingError` is
#' signaled when no optimizer is given and no default applies.
#'
#' @param alg An optimization algorithm.
#' @param prob An `OptimizationProblem`.
#' @param iterations,time,fevals,gevals,hevals Fields of `OptimizationStats()`.
#' @param f An `OptimizationFunction`.
#' @param msg A message.
#' @param state,loss Arguments of a callback.
#' @param ... Solver options.
#' @param level Verbosity: `"none"`, `"warn"`, `"info"` or `"all"`.
#' @return A trait value, a cache, a solution, a condition or a settings
#'   object.
#' @examples
#' prob <- OptimizationProblem(function(x, p) -(x[1] - 1)^2, 0, sense = MaxSense)
#' sol <- solve(prob, NelderMead())
#' sol$u; sol$objective
#' requiresgradient(LBFGSB())
#' cache <- init(prob, NelderMead())
#' solve_integrator(cache)$u
#' @name optimization_interface
NULL

#' @rdname optimization_interface
#' @export
MinSense <- structure("MinSense", class = "ObjSense")

#' @rdname optimization_interface
#' @export
MaxSense <- structure("MaxSense", class = "ObjSense")

#' @rdname optimization_interface
#' @export
ObjSense <- list(MinSense = MinSense, MaxSense = MaxSense)

#' @export
print.ObjSense <- function(x, ...) {
  cat(unclass(x), "\n")
  invisible(x)
}

#' @rdname optimization_interface
#' @export
DEFAULT_CALLBACK <- function(state, loss) FALSE

#' @rdname optimization_interface
#' @export
DEFAULT_DATA <- structure(list(), class = "OptimizationDefaultData")

#' @rdname optimization_interface
#' @export
NoAD <- function() .adtype("NoAD")

#' @rdname optimization_interface
#' @export
OptimizationVerbosity <- function(level = c("warn", "none", "info", "all"))
  structure(list(level = match.arg(level)), class = "OptimizationVerbosity")

#' @rdname optimization_interface
#' @export
OptimizationStats <- function(iterations = 0L, time = 0, fevals = 0L, gevals = 0L, hevals = 0L)
  structure(list(iterations = as.integer(iterations), time = time, fevals = as.integer(fevals),
                 gevals = as.integer(gevals), hevals = as.integer(hevals)), class = "OptimizationStats")

#' @export
print.OptimizationStats <- function(x, ...) {
  cat(sprintf("OptimizationStats: %d iterations, %d f, %d gradient, %d Hessian evaluations, %.3g s\n",
              x$iterations, x$fevals, x$gevals, x$hevals, x$time))
  invisible(x)
}

#' @rdname optimization_interface
#' @export
IncompatibleOptimizerError <- function(msg)
  structure(class = c("IncompatibleOptimizerError", "error", "condition"), list(message = msg, call = NULL))

#' @rdname optimization_interface
#' @export
OptimizerMissingError <- function(msg = "no optimization algorithm was given")
  structure(class = c("OptimizerMissingError", "error", "condition"), list(message = msg, call = NULL))

# traits of the R optimizers: which derivatives they use and what they accept
.opt_traits <- function(alg) {
  name <- alg$name %||% class(alg)[1L]
  t <- list(grad = FALSE, hess = FALSE, consjac = FALSE, conshess = FALSE, rbounds = FALSE, rcons = FALSE,
            bounds = FALSE, cons = FALSE, callback = TRUE)
  if (is.list(alg$traits)) return(utils::modifyList(t, alg$traits))
  switch(name,
    LBFGSB = , LBFGS = { t$grad <- TRUE; t$bounds <- TRUE },
    Nlminb = { t$grad <- TRUE; t$bounds <- TRUE },
    AugLag = { t$grad <- TRUE; t$bounds <- TRUE; t$cons <- TRUE; t$consjac <- TRUE },
    BranchAndBound = { t$bounds <- TRUE; t$cons <- TRUE; t$rbounds <- TRUE },
    NSGA2 = { t$bounds <- TRUE; t$rbounds <- TRUE; t$cons <- TRUE },
    WeightedSumScan = { t$grad <- TRUE; t$bounds <- TRUE; t$cons <- TRUE })
  t
}

#' @rdname optimization_interface
#' @export
requiresgradient <- function(alg) .opt_traits(alg)$grad
#' @rdname optimization_interface
#' @export
requireshessian <- function(alg) .opt_traits(alg)$hess
#' @rdname optimization_interface
#' @export
requiresconsjac <- function(alg) .opt_traits(alg)$consjac
#' @rdname optimization_interface
#' @export
requiresconshess <- function(alg) .opt_traits(alg)$conshess
#' @rdname optimization_interface
#' @export
requiresbounds <- function(alg) .opt_traits(alg)$rbounds
#' @rdname optimization_interface
#' @export
requiresconstraints <- function(alg) .opt_traits(alg)$rcons
#' @rdname optimization_interface
#' @export
allowsbounds <- function(alg) .opt_traits(alg)$bounds
#' @rdname optimization_interface
#' @export
allowsconstraints <- function(alg) .opt_traits(alg)$cons
#' @rdname optimization_interface
#' @export
allowscallback <- function(alg) .opt_traits(alg)$callback

#' @rdname optimization_interface
#' @export
lag_hess_structure <- function(f) {
  n <- f$n %||% length(f$u0 %||% numeric())
  if (!is.null(f$lag_hess_prototype)) return(f$lag_hess_prototype != 0)
  if (!n) return(NULL)
  matrix(TRUE, n, n)
}

#' @export
init.OptimizationProblem <- function(prob, alg = NULL, ...)
  structure(list2env(list(prob = prob, alg = alg, opts = list(...), sol = NULL)), class = "OptimizationCache")

#' @rdname optimization_interface
#' @export
OptimizationCache <- function(prob, alg = NULL, ...) init.OptimizationProblem(prob, alg, ...)

#' @export
print.OptimizationCache <- function(x, ...) {
  cat("OptimizationCache (", x$alg$name %||% "default", ")\n", sep = "")
  invisible(x)
}
