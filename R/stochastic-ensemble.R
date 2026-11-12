# Ensemble simulations, after SciMLBase's EnsembleProblem, EnsembleSolution
# and EnsembleSummary (src/ensemble/ensemble_problems.jl,
# ensemble_solutions.jl and ensemble_analysis.jl), used by the SDE and jump
# tutorials of ModelingToolkit. Trajectories run serially (EnsembleSerial).

#' Ensemble problems
#'
#' `EnsembleProblem(prob)` repeats a problem (usually stochastic) many times.
#' `prob_func(prob, i, repeat)` returns the problem of trajectory `i` (by
#' default `prob` itself), `output_func(sol, i)` returns `list(value,
#' rerun)` with the value kept for trajectory `i` and whether to run it again
#' (by default the solution, no rerun), and `reduction(u, batch, I)` combines
#' the values of a batch, returning `list(u, done)`. `solve(ensemble, alg,
#' trajectories = n)` returns an `EnsembleSolution`, a list of trajectory
#' values with fields `u` and `elapsed`. A `seed` makes the whole ensemble
#' reproducible; trajectories draw from one random stream in turn.
#'
#' @param prob A problem.
#' @param prob_func,output_func,reduction Functions as described above.
#' @param u_init Initial value of the reduction.
#' @param safetycopy Unused; kept for compatibility.
#' @return An `EnsembleProblem`.
#' @examples
#' prob <- SDEProblem(function(u, p, t) 1.5 * u, function(u, p, t) 0.3 * u, 1, c(0, 1))
#' ens <- EnsembleProblem(prob)
#' sim <- solve(ens, EM(), trajectories = 20, dt = 0.05, seed = 1)
#' summ <- EnsembleSummary(sim)
#' tail(summ$u, 1)
#' @export
EnsembleProblem <- function(prob, prob_func = function(prob, i, repeat_) prob,
                            output_func = function(sol, i) list(sol, FALSE),
                            reduction = function(u, data, I) list(c(u, data), FALSE),
                            u_init = NULL, safetycopy = TRUE) {
  structure(list(prob = prob, prob_func = prob_func, output_func = output_func,
                 reduction = reduction, u_init = u_init),
            class = c("EnsembleProblem", "SciMLProblem"))
}

#' @rdname EnsembleProblem
#' @export
EnsembleSerial <- function() structure(list(name = "EnsembleSerial"), class = "EnsembleAlgorithm")

#' @rdname EnsembleProblem
#' @param cores Number of worker processes (default `getOption("mc.cores", 2)`).
#' @export
EnsembleThreads <- function(cores = NULL)
  structure(list(name = "EnsembleThreads", cores = cores), class = "EnsembleAlgorithm")

#' @rdname EnsembleProblem
#' @param cl A cluster of the parallel package; by default one is created
#'   for the solve.
#' @export
EnsembleDistributed <- function(cl = NULL, cores = NULL)
  structure(list(name = "EnsembleDistributed", cl = cl, cores = cores), class = "EnsembleAlgorithm")

#' @rdname EnsembleProblem
#' @export
EnsembleSplitThreads <- function(cl = NULL, cores = NULL)
  structure(list(name = "EnsembleSplitThreads", cl = cl, cores = cores), class = "EnsembleAlgorithm")

#' @rdname EnsembleProblem
#' @param backend,cpu_offload Options of the Julia signature; R has no GPU
#'   arrays, so the trajectories run on CPU workers as with `EnsembleThreads()`.
#' @export
EnsembleGPUKernel <- function(backend = NULL, cpu_offload = 0, cores = NULL)
  structure(list(name = "EnsembleGPUKernel", cores = cores), class = "EnsembleAlgorithm")

#' @rdname EnsembleProblem
#' @export
EnsembleGPUArray <- function(backend = NULL, cpu_offload = 0, cores = NULL)
  structure(list(name = "EnsembleGPUArray", cores = cores), class = "EnsembleAlgorithm")

#' @rdname EnsembleProblem
#' @param u,elapsed,converged Values of an ensemble solution.
#' @export
EnsembleSolution <- function(u, elapsed = 0, converged = TRUE)
  structure(list(u = u, elapsed = elapsed, converged = converged), class = "EnsembleSolution")

# How a batch of trajectories is run: in turn (EnsembleSerial, one random
# stream), or by workers (forked processes for the thread and GPU kernel
# variants, a socket cluster for the distributed ones), each trajectory with
# its own L'Ecuyer-CMRG stream, so the results do not depend on the number
# of workers.
.ensemble_mapper <- function(alg) {
  name <- alg$name %||% "EnsembleSerial"
  if (name == "EnsembleSerial") return(list(map = function(I, FUN) lapply(I, FUN), close = function() NULL))
  cores <- alg$cores %||% getOption("mc.cores", 2L)
  streams <- function(I) {
    old <- RNGkind()[1L]
    on.exit(RNGkind(old))
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) stats::runif(1L)
    RNGkind("L'Ecuyer-CMRG")
    s <- .Random.seed
    lapply(I, function(i) s <<- parallel::nextRNGStream(s))
  }
  with_stream <- function(FUN) function(args) {
    RNGkind("L'Ecuyer-CMRG")
    assign(".Random.seed", args$seed, envir = globalenv())
    FUN(args$i)
  }
  if (name %in% c("EnsembleThreads", "EnsembleGPUKernel", "EnsembleGPUArray")) {
    return(list(map = function(I, FUN) {
      args <- Map(function(i, sd) list(i = i, seed = sd), I, streams(I))
      cores_ok <- if (.Platform$OS.type == "windows") 1L else cores
      parallel::mclapply(args, with_stream(FUN), mc.cores = cores_ok, mc.set.seed = FALSE)
    }, close = function() NULL))
  }
  cl <- alg$cl
  own <- is.null(cl)
  if (own) {
    cl <- parallel::makePSOCKcluster(cores)
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(modelcompR)))
  }
  list(map = function(I, FUN) {
    args <- Map(function(i, sd) list(i = i, seed = sd), I, streams(I))
    parallel::parLapply(cl, args, with_stream(FUN))
  }, close = function() if (own) parallel::stopCluster(cl))
}

#' @rdname EnsembleProblem
#' @param a An `EnsembleProblem`.
#' @param b The solver algorithm for each trajectory.
#' @param ensemblealg How trajectories run: `EnsembleSerial()` (the default,
#'   one random stream used in turn), `EnsembleThreads()` (forked workers;
#'   serial on Windows), `EnsembleDistributed()` and
#'   `EnsembleSplitThreads()` (a socket cluster of the parallel package),
#'   or `EnsembleGPUKernel()` and `EnsembleGPUArray()` (CPU workers). The
#'   parallel ones give each trajectory its own L'Ecuyer-CMRG stream.
#' @param trajectories Number of trajectories.
#' @param batch_size Trajectories per reduction batch.
#' @param seed Seed for the random number generator.
#' @param ... Solver options for each trajectory.
#' @export
solve.EnsembleProblem <- function(a, b = NULL, ensemblealg = EnsembleSerial(), ...,
                                  trajectories, batch_size = trajectories, seed = NULL) {
  if (inherits(b, "EnsembleAlgorithm")) {
    ensemblealg <- b
    b <- NULL
  }
  if (missing(trajectories)) stop("give the number of trajectories", call. = FALSE)
  start <- proc.time()[["elapsed"]]
  dots <- list(...)
  one <- function(i) {
    rep_ <- 1L
    repeat {
      prob_i <- a$prob_func(a$prob, i, rep_)
      sol <- if (is.null(b)) do.call(solve, c(list(prob_i), dots)) else do.call(solve, c(list(prob_i, b), dots))
      out <- a$output_func(sol, i)
      rerun <- is.list(out) && length(out) == 2L && is.logical(out[[2L]]) && length(out[[2L]]) == 1L
      val <- if (rerun) out[[1L]] else out
      if (!(rerun && isTRUE(out[[2L]]))) break
      rep_ <- rep_ + 1L
    }
    val
  }
  run <- function() {
    mapper <- .ensemble_mapper(ensemblealg)
    on.exit(mapper$close())
    u <- a$u_init
    done <- 0L
    while (done < trajectories) {
      I <- seq.int(done + 1L, min(trajectories, done + batch_size))
      batch <- mapper$map(I, one)
      red <- a$reduction(u, batch, I)
      u <- red[[1L]]
      done <- max(I)
      if (isTRUE(red[[2L]])) break
    }
    u
  }
  u <- .with_seed(seed, run())
  structure(list(u = u, elapsed = proc.time()[["elapsed"]] - start, converged = TRUE),
            class = "EnsembleSolution")
}

#' @export
`$.EnsembleSolution` <- function(x, name) .subset2(x, name)

#' @export
`[[.EnsembleSolution` <- function(x, i) .subset2(x, "u")[[i]]

#' @export
length.EnsembleSolution <- function(x) length(.subset2(x, "u"))

#' @export
print.EnsembleSolution <- function(x, ...) {
  cat("EnsembleSolution with ", length(x), " trajectories\n", sep = "")
  invisible(x)
}

#' Summary statistics of an ensemble
#'
#' `EnsembleSummary(sim, ts)` computes, at the times `ts` (by default the
#' saved times of the first trajectory), the componentwise mean `u`,
#' variance `v` and the `quantiles` `qlow` and `qhigh` over the trajectories
#' of an `EnsembleSolution` whose values are solutions. Each is a matrix with
#' one row per time. Solutions are interpolated at `ts` (piecewise constant
#' for jump processes, where the value at a jump time is the one after it).
#'
#' @param sim An `EnsembleSolution`.
#' @param ts Times.
#' @param quantiles Lower and upper quantile levels.
#' @param idxs Optional variables or expressions to summarize instead of the
#'   states.
#' @return An `EnsembleSummary` with fields `t`, `u`, `v`, `qlow`, `qhigh`
#'   and `num_monte`.
#' @examples
#' prob <- SDEProblem(function(u, p, t) u, function(u, p, t) 0.2 * u, 1, c(0, 1))
#' sim <- solve(EnsembleProblem(prob), EM(), trajectories = 10, dt = 0.1, seed = 3)
#' EnsembleSummary(sim, ts = c(0, 0.5, 1))$u
#' @export
EnsembleSummary <- function(sim, ts = NULL, quantiles = c(0.05, 0.95), idxs = NULL) {
  sols <- .subset2(sim, "u")
  if (!length(sols) || !inherits(sols[[1L]], "SciMLSolution"))
    stop("EnsembleSummary needs an ensemble of solutions", call. = FALSE)
  if (is.null(ts)) ts <- sols[[1L]]$t
  vals <- lapply(sols, function(s) {
    out <- if (is.null(idxs)) s(ts) else s(ts, idxs = idxs)
    if (!is.matrix(out)) out <- matrix(out, nrow = length(ts))
    out
  })
  arr <- array(unlist(vals), dim = c(length(ts), ncol(vals[[1L]]), length(vals)),
               dimnames = list(NULL, colnames(vals[[1L]]), NULL))
  stat <- function(fn) apply(arr, c(1L, 2L), fn)
  u <- stat(mean)
  structure(list(t = ts, u = u, v = stat(stats::var),
                 qlow = stat(function(z) stats::quantile(z, quantiles[1L], names = FALSE)),
                 qhigh = stat(function(z) stats::quantile(z, quantiles[2L], names = FALSE)),
                 num_monte = length(sols)),
            class = "EnsembleSummary")
}

#' @export
print.EnsembleSummary <- function(x, ...) {
  cat("EnsembleSummary of ", x$num_monte, " trajectories at ", length(x$t), " times\n", sep = "")
  invisible(x)
}

#' @export
plot.EnsembleSummary <- function(x, ...) {
  graphics::matplot(x$t, x$u, type = "l", lty = 1, xlab = "t", ylab = "", ...)
  graphics::matlines(x$t, x$qlow, lty = 2)
  graphics::matlines(x$t, x$qhigh, lty = 2)
  invisible(x)
}
