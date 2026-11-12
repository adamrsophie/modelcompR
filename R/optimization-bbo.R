# Multiobjective global optimization: the BBO_borg_moea wrapper of
# OptimizationBBO with the Borg MOEA of BlackBoxOptim.jl (borg_moea.jl and
# the parts of the BlackBoxOptim framework it runs on: the epsilon box
# dominance fitness, the epsilon box archive, the population with one
# transient member, the tournament selector, the seven recombination
# operators, the mutation clocks, the random bound embedding and the run
# controller), ported to R.

#' Borg MOEA through BlackBoxOptim
#'
#' `BBO_borg_moea()` solves a multiobjective problem (an
#' [OptimizationProblem()] built on a [MultiObjectiveOptimizationFunction()])
#' with the Borg MOEA of BlackBoxOptim.jl (Hadka and Reed 2013), ported to
#' R, as the BBO_borg_moea wrapper of OptimizationBBO does. Borg keeps an
#' archive of the points that are not dominated in epsilon boxes (boxes of
#' side `epsilon` in the space of the objectives) and a population. Each
#' step picks one of seven recombination operators (differential evolution,
#' simulated binary, simplex, two parent centric and two unimodal normal
#' distribution crossovers) with probabilities that follow how many archive
#' points each operator made, takes the parents by tournament from the
#' population and one parent from the archive, and puts the child into the
#' population in place of a point it dominates. When the population no
#' longer fits the archive, or the archive stops making progress, the
#' population is restarted from the archive and mutated archive points.
#'
#' The problem needs finite bounds (`lb` and `ub`); `u0` is put into the
#' initial population. The `solve` arguments map as in OptimizationBBO:
#' `maxiters` to `MaxSteps` (default 10000; the run stops after one step
#' more), `maxtime` to `MaxTime` (seconds) and `abstol` to
#' `MinDeltaFitnessTolerance` (not used by Borg); `reltol` is not used. A
#' `callback(state, objective)` is called after every step with
#' `state$iter`, `state$u` and `state$objective` (the best point so far,
#' the one with the least sum of the objectives on the archive) and
#' `state$original` (the state of the run); `TRUE` stops the run. Other
#' arguments are BlackBoxOptim parameters: `FitnessScheme` (a
#' [ParetoFitnessScheme()]; by default one with as many objectives as the
#' function returns), `PopulationSize` (50), `MaxFuncEvals`,
#' `MaxStepsWithoutProgress` (10000 calls to the archive without epsilon
#' progress), `MaxNumStepsWithoutFuncEvals` (100), `epsilon` (0.1, the side
#' of the epsilon boxes, one value or one per objective), `tau` (0.02, the
#' tournament size as a fraction of the population), `gamma` (4, the
#' population to archive ratio), `gamma_delta` (0.25), `theta` (0.9),
#' `zeta` (1), `RestartCheckPeriod` (1000), `OperatorsUpdatePeriod` (100),
#' `MaxStepsWithoutEpsProgress` (100), `MaxArchiveSize` (100000), the
#' operator parameters `DEX_f` (0.6), `DEX_cr` (0.7), `SBX_p` (0.2),
#' `SBX_eta` (3), `SPX_epsilon` (1.1), `PCX_zeta` (0.1), `PCX_eta` (0.5),
#' `UNDX_zeta` (0.5), `UNDX_eta` (0.1) and `PM_eta` (50), and `Population`
#' (a matrix with one initial point per column; the last column is the
#' transient member). Greek parameter names are spelled out.
#'
#' The solution's `u` and `objective` are the best point and its
#' objectives, as in OptimizationBBO; `original` holds the stop reason and
#' the whole archive (`frontier$X`, one point per row, and `frontier$F`,
#' their objectives). The return code is `MaxIters` when the run reaches
#' `MaxSteps`, `MaxTime` at `MaxTime` and `Default` otherwise (with a
#' warning that names the stop reason). The run uses R's random numbers:
#' call [set.seed()] for a repeatable run.
#'
#' @return An algorithm object.
#' @examples
#' sphere_rastrigin <- function(x, p) c(sum(x^2), sum(x^2 - 10 * cos(2 * pi * x) + 10))
#' prob <- OptimizationProblem(MultiObjectiveOptimizationFunction(sphere_rastrigin), c(0.25, 0.25),
#'                             lb = c(0, 0), ub = c(2, 2))
#' set.seed(1)
#' sol <- solve(prob, BBO_borg_moea(), maxiters = 300)
#' sol$u; sol$objective
#' @export
BBO_borg_moea <- function()
  structure(list(name = "BBO_borg_moea", method = "borg_moea",
                 traits = list(bounds = TRUE, rbounds = TRUE, callback = TRUE)),
            class = c("BBO_borg_moea", "OptimizationAlgorithm"))

#' Pareto dominance fitness scheme
#'
#' `ParetoFitnessScheme()` is the BlackBoxOptim fitness scheme for `N`
#' objectives compared by Pareto dominance, the `FitnessScheme` argument of
#' `solve` with [BBO_borg_moea()]. Borg turns it into epsilon box
#' dominance.
#'
#' @param N The number of objectives.
#' @param is_minimizing Whether the objectives are minimized (`TRUE`) or
#'   maximized.
#' @param aggregator The function that maps the objectives of a point to
#'   one number; the best point of a run is the one with the least (or the
#'   greatest, when maximizing) aggregate on the archive.
#' @return A fitness scheme object.
#' @examples
#' ParetoFitnessScheme(2)
#' @export
ParetoFitnessScheme <- function(N, is_minimizing = TRUE, aggregator = sum)
  structure(list(N = as.integer(N), is_minimizing = isTRUE(is_minimizing), aggregator = aggregator),
            class = c("ParetoFitnessScheme", "TupleFitnessScheme", "FitnessScheme"))

#' @export
print.ParetoFitnessScheme <- function(x, ...) {
  cat("ParetoFitnessScheme{", x$N, ", Float64, ", if (x$is_minimizing) "true" else "false", "}\n", sep = "")
  invisible(x)
}

.mtk$hooks$opt_solve_BBO_borg_moea <- function(a, alg, opts) {
  .opt_check_traits(a, alg)
  maxiters <- opts$maxiters
  if (!is.null(maxiters)) {
    if (maxiters <= 0) stop("The number of maxiters has to be a non-negative and non-zero number.", call. = FALSE)
    maxiters <- round(maxiters)
  }
  maxtime <- opts$maxtime
  if (!is.null(maxtime) && maxtime <= 0)
    stop("The maximum time has to be a non-negative and non-zero number.", call. = FALSE)
  if (!is.null(opts$reltol)) message("common reltol is currently not used by BBO_borg_moea(:borg_moea)")
  p <- a$p
  cb <- opts$callback
  # the bbsetup parameters of OptimizationBBO
  params <- opts[setdiff(names(opts), c("maxiters", "maxtime", "abstol", "reltol", "callback", "verbose",
                                        "progress"))]
  params$Method <- "borg_moea"
  params$SearchRange <- list(lower = as.numeric(a$lb), upper = as.numeric(a$ub))
  if (!is.null(maxiters)) params$MaxSteps <- maxiters
  if (!is.null(maxtime)) params$MaxTime <- maxtime
  if (!is.null(opts$abstol)) params$MinDeltaFitnessTolerance <- opts$abstol
  callback <- if (is.function(cb)) function(b) {
    best <- b$archive$best
    state <- list(iter = b$num_steps, u = best$x, p = p, objective = best$fit$orig, original = b)
    halt <- cb(state, best$fit$orig)
    if (!is.logical(halt) || length(halt) != 1L || is.na(halt))
      stop("The callback should return a boolean `halt` for whether to stop the optimization process.",
           call. = FALSE)
    if (halt) b$shutdown <- TRUE
  }
  res <- .bbo_optimize(function(x) as.numeric(a$f$f(x, p)), params, as.numeric(a$u0), callback)
  sol <- .opt_solution(a, alg, res$best_candidate, res$best_fitness, .bbo_retcode(res$stop_reason), res)
  .opt_set_fields(sol, stats = OptimizationStats(iterations = res$iterations, time = res$elapsed_time,
                                                 fevals = res$f_calls))
}

# deduce_retcode of OptimizationBase on a BlackBoxOptim stop reason
.bbo_retcode <- function(reason) {
  map <- c("Delta fitness .* below tolerance .*" = "Success", "Fitness .* within tolerance .* of optimum" = "Success",
           "Terminated" = "Terminated", "MaxIters|MAXITERS_EXCEED|Max number of steps .* reached" = "MaxIters",
           "MaxTime|TIME_LIMIT" = "MaxTime", "Max time" = "MaxTime", "Optimization completed" = "Success",
           "Convergence achieved" = "Success")
  for (pat in names(map)) if (grepl(pat, reason)) return(map[[pat]])
  warning("Unrecognized stop reason: ", reason, ". Defaulting to ReturnCode.Default.", call. = FALSE)
  "Default"
}

.bbo_jlnum <- function(x) sub("^\\[(.*)\\]$", "\\1", .qd_jlvec(x))

# ---------------------------------------------------------------------------
# Random numbers. Every draw goes through these functions, which the tests
# can feed with the numbers BlackBoxOptim drew (tools/bbo-probe.jl).

.bbo_rng <- new.env(parent = emptyenv())

.bbo_draw <- function(n, gen) {
  q <- .bbo_rng$queue
  if (is.null(q)) return(gen(n))
  .bbo_rng$queue <- q[-seq_len(n)]
  q[seq_len(n)]
}
.bbo_runif <- function(n = 1L) .bbo_draw(n, stats::runif)
.bbo_rnorm <- function(n = 1L) .bbo_draw(n, stats::rnorm)
# rand(1:n)
.bbo_randint <- function(n) .bbo_draw(1L, function(k) sample.int(n, 1L))
# rand(Bool)
.bbo_coin <- function() as.logical(.bbo_draw(1L, function(k) stats::runif(1L) < 0.5))
# the index drawn by the weights w (Categorical, and StatsBase sample)
.bbo_wsample <- function(w) {
  t <- .bbo_runif() * sum(w)
  i <- 1L
  cw <- w[1L]
  while (cw < t && i < length(w)) {
    i <- i + 1L
    cw <- cw + w[i]
  }
  i
}

# ---------------------------------------------------------------------------
# Epsilon box dominance (ntuple_fitness.jl)

.bbo_imax <- 0x1.0000000000000p+63  # typemax(Int)
.bbo_maxuu <- 9223372036854774784  # prevfloat(Float64(typemax(Int)))
.bbo_tol <- 10 * .Machine$double.eps

.bbo_scheme <- function(fs, eps) {
  if (!inherits(fs, "TupleFitnessScheme"))
    stop("BorgMOEA can only solve problems with `TupleFitnessScheme`", call. = FALSE)
  N <- fs$N
  eps <- as.numeric(eps)
  if (length(eps) == 1L) {
    if (!(eps > 0)) stop("\u03F5 must be positive", call. = FALSE)
    eps <- rep(eps, N)
  } else {
    if (length(eps) != N)
      stop("The length of \u03F5 vector (", length(eps), ") does not match the specified fitness dimensions (",
           N, ")", call. = FALSE)
    if (!all(eps > 0)) stop("\u03F5 must be positive", call. = FALSE)
  }
  agg <- fs$aggregator
  # Julia's sum adds left to right
  if (identical(agg, sum)) agg <- function(f) Reduce(`+`, f)
  list(N = N, eps = eps, is_minimizing = fs$is_minimizing, agg = agg)
}

# IndexedTupleFitness: the objectives u, their aggregate, the epsilon box
# index and the distance of u from the corner of its box
.bbo_indexed <- function(u, sch) {
  n <- length(u)
  ix <- numeric(n)
  rem <- numeric(n)
  for (i in seq_len(n)) {
    if (is.na(u[i])) {
      ix[i] <- if (sch$is_minimizing) .bbo_imax else -.bbo_imax
      next
    }
    uu <- u[i] / sch$eps[i]
    if (uu > .bbo_maxuu) {
      warning("Clamping the epsilon-box index (value too large). Probably you need to increase the \u03F5 value ",
              "of your fitness scheme.", call. = FALSE)
      ix[i] <- .bbo_imax
    } else if (uu < -.bbo_maxuu) {
      warning("Clamping the epsilon-box index (value too small). Probably you need to increase the \u03F5 value ",
              "of your fitness scheme.", call. = FALSE)
      ix[i] <- -.bbo_imax
    } else if (sch$is_minimizing) {
      ix[i] <- floor(min(uu + .bbo_tol, .bbo_maxuu))
      rem[i] <- uu - ix[i]
    } else {
      ix[i] <- ceiling(uu)
      rem[i] <- ix[i] - uu
    }
  }
  s <- 0
  for (r in rem) s <- s + r^2
  list(orig = u, agg = sch$agg(u), index = ix, dist = sqrt(s))
}

.bbo_nafitness <- function(N) list(orig = rep(NaN, N), agg = NaN, index = rep(.bbo_imax, N), dist = 0)

# hat_compare_eps_box of the fitness with index ui and distance ud against
# each row of the index matrix VI (distances vd): -1 when it dominates the
# row, 1 when the row dominates it, 0 otherwise
.bbo_hat_rows <- function(ui, ud, VI, vd, minimizing) {
  gt <- lt <- logical(nrow(VI))
  for (j in seq_along(ui)) {
    gt <- gt | ui[j] > VI[, j]
    lt <- lt | ui[j] < VI[, j]
  }
  comp <- gt - lt
  if (!minimizing) comp <- -comp
  eq <- which(!gt & !lt)
  if (length(eq)) {
    d <- ud - vd[eq]
    comp[eq] <- (d > .bbo_tol) - (d < -.bbo_tol)
  }
  comp
}

.bbo_hat <- function(u, v, minimizing)
  .bbo_hat_rows(u$index, u$dist, matrix(v$index, 1L), v$dist, minimizing)

# ---------------------------------------------------------------------------
# The epsilon box archive (archives/epsbox_archive.jl); the frontier is
# held by rows

.bbo_archive <- function(sch, D, max_size) {
  a <- new.env(parent = emptyenv())
  a$sch <- sch
  a$num_candidates <- 0L
  a$last_progress <- 0L
  a$last_restart <- 0L
  a$n_restarts <- 0L
  a$n_oversize_inserts <- 0L
  a$max_size <- max_size
  N <- sch$N
  a$FO <- matrix(0, 0L, N)
  a$FI <- matrix(0, 0L, N)
  a$FD <- numeric(0)
  a$FA <- numeric(0)
  a$X <- matrix(0, 0L, D)
  a$tag <- integer(0)
  a$nfe <- integer(0)
  a$nrs <- integer(0)
  a$best <- NULL
  a
}

.bbo_front_elem <- function(a, i)
  list(fit = list(orig = a$FO[i, ], agg = a$FA[i], index = a$FI[i, ], dist = a$FD[i]), x = a$X[i, ],
       tag = a$tag[i], nfe = a$nfe[i], nrs = a$nrs[i])

.bbo_noprogress <- function(a, since_restart = FALSE)
  a$num_candidates - (if (since_restart) max(a$last_progress, a$last_restart) else a$last_progress)

.bbo_add_candidate <- function(a, fit, x, tag = 0L, nfe = -1L) {
  a$num_candidates <- a$num_candidates + 1L
  if (nfe == -1L) nfe <- a$num_candidates
  minim <- a$sch$is_minimizing
  k <- length(a$FD)
  h <- if (k) .bbo_hat_rows(fit$index, fit$dist, a$FI, a$FD, minim) else numeric(0)
  if (any(h == 1)) {
    # an element of the frontier dominates the candidate
    if (k <= a$max_size) a$n_oversize_inserts <- 0L
    return(invisible(a))
  }
  same <- if (k) which(rowSums(a$FI == matrix(fit$index, k, length(fit$index), byrow = TRUE)) ==
                         length(fit$index)) else integer(0)
  elem <- list(fit = fit, x = x, tag = tag, nfe = nfe, nrs = a$n_restarts)
  if (length(same)) {
    # the box is on the frontier: no epsilon progress
    i <- same[1L]
    hat <- h[i]
    if (hat < 0) {
      a$FO[i, ] <- fit$orig; a$FI[i, ] <- fit$index; a$FD[i] <- fit$dist; a$FA[i] <- fit$agg
      a$X[i, ] <- x; a$tag[i] <- tag; a$nfe[i] <- nfe; a$nrs[i] <- a$n_restarts
    } else {
      elem <- .bbo_front_elem(a, i)
    }
  } else {
    a$last_progress <- a$num_candidates
    hat <- -1
    keep <- h != -1
    a$FO <- rbind(a$FO[keep, , drop = FALSE], fit$orig)
    a$FI <- rbind(a$FI[keep, , drop = FALSE], fit$index)
    a$FD <- c(a$FD[keep], fit$dist)
    a$FA <- c(a$FA[keep], fit$agg)
    a$X <- rbind(a$X[keep, , drop = FALSE], x)
    a$tag <- c(a$tag[keep], tag)
    a$nfe <- c(a$nfe[keep], nfe)
    a$nrs <- c(a$nrs[keep], a$n_restarts)
    if (length(a$FD) > a$max_size) a$n_oversize_inserts <- a$n_oversize_inserts + 1L
  }
  if (is.null(a$best)) {
    a$best <- elem
  } else if (hat < 0) {
    d <- a$best$fit$agg - elem$fit$agg
    if (identical(a$best$fit$index, elem$fit$index) || (d > 0 && minim) || (d < 0 && !minim)) a$best <- elem
  }
  if (length(a$FD) <= a$max_size) a$n_oversize_inserts <- 0L
  invisible(a)
}

# the tag (operator) counts of the frontier, each discounted by theta for
# every restart it survived
.bbo_tagcounts <- function(a, theta, nops) {
  if (!(theta > 0 && theta <= 1)) stop("\u03B8 (", theta, ") should be in (0.0, 1.0] range", call. = FALSE)
  cnt <- numeric(nops)
  for (i in seq_along(a$tag)) {
    t <- a$tag[i]
    if (t > 0) cnt[t] <- cnt[t] + theta^(a$n_restarts - a$nrs[i])
  }
  cnt
}

# ---------------------------------------------------------------------------
# The population: one column per individual, the last ntransient columns
# transient, and the fitness by rows

.bbo_population <- function(X, N, ntransient) {
  pop <- new.env(parent = emptyenv())
  pop$X <- X
  n <- ncol(X)
  pop$ntransient <- ntransient
  pop$popsize <- n - ntransient
  pop$FO <- matrix(NaN, n, N)
  pop$FI <- matrix(.bbo_imax, n, N)
  pop$FD <- numeric(n)
  pop$FA <- rep(NaN, n)
  pop
}

.bbo_pop_fit <- function(pop, i)
  list(orig = pop$FO[i, ], agg = pop$FA[i], index = pop$FI[i, ], dist = pop$FD[i])

.bbo_pop_set <- function(pop, i, x, fit) {
  if (!is.null(x)) pop$X[, i] <- x
  pop$FO[i, ] <- fit$orig; pop$FI[i, ] <- fit$index; pop$FD[i] <- fit$dist; pop$FA[i] <- fit$agg
}

.bbo_lhs <- function(lo, hi, n) {
  D <- length(lo)
  R <- matrix(0, D, n)
  for (i in seq_len(D)) {
    step <- (hi[i] - lo[i]) / n
    cube <- lo[i] + step * (seq_len(n) - 1 + .bbo_runif(n))
    R[i, ] <- cube[sample.int(n)]
  }
  R
}

# ---------------------------------------------------------------------------
# Genetic operators (genetic_operators/): each crossover maps the
# population matrix and the parent columns to its children

.bbo_de <- function(cr, f) function(X, pix) {
  D <- nrow(X)
  m <- .bbo_randint(D)
  u <- numeric(D)
  u[-m] <- .bbo_runif(D - 1L)
  cross <- seq_len(D) == m | u <= cr
  t <- X[, pix[3L]]
  t[cross] <- X[cross, pix[3L]] + f * (X[cross, pix[1L]] - X[cross, pix[2L]])
  list(t)
}

.bbo_sbx <- function(p, eta) {
  if (!(p >= 0 && p <= 1)) stop("p must be within [0,1]", call. = FALSE)
  if (!(eta > 0)) stop("\u03B7 must be positive", call. = FALSE)
  eta_exp <- 1 / (eta + 1)
  function(X, pix) {
    a <- X[, pix[1L]]
    b <- X[, pix[2L]]
    t1 <- a
    t2 <- b
    for (i in seq_along(a)) {
      if (.bbo_runif() <= p) {
        u <- .bbo_runif()
        beta <- (if (u < 0.5) 2 * u else 1 / (2 * (1 - u)))^eta_exp
        avg <- 0.5 * (a[i] + b[i])
        delta <- 0.5 * (a[i] - b[i])
        if (.bbo_coin()) delta <- -delta
        t1[i] <- avg + beta * delta
        t2[i] <- avg - beta * delta
      }
    }
    list(t1, t2)
  }
}

.bbo_spx <- function(np, eps) function(X, pix) {
  offset <- numeric(nrow(X))
  psum <- X[, pix[1L]]
  for (i in seq_len(np - 1L)) {
    s <- .bbo_runif()^(1 / i)
    nxt <- X[, pix[i + 1L]]
    offset <- s * (eps * (X[, pix[i]] - nxt) + offset)
    psum <- psum + nxt
  }
  s <- (1 - eps) / np
  list(s * psum + eps * X[, pix[np]] + offset)
}

.bbo_center <- function(P) Reduce(`+`, lapply(seq_len(ncol(P)), function(j) P[, j])) / ncol(P)

.bbo_pcx <- function(np, zeta, eta) function(X, pix) {
  P <- X[, pix, drop = FALSE]
  center <- .bbo_center(P)
  p1 <- P[, 1L] - center
  others <- P[, -1L, drop = FALSE] - center
  sq <- Reduce(`+`, p1^2)
  if (isTRUE(sq > 0)) {
    dots <- colSums(p1 * others) / sq
    others <- others - matrix(dots, nrow(others), ncol(others), byrow = TRUE) * p1
  }
  sd <- mean(sqrt(rowSums(others^2)))
  if (isTRUE(sd > 1e-8)) {
    s <- svd(others)
    sv <- s$d * sd * zeta / sqrt(Reduce(`+`, s$d^2))
    t <- drop(s$u %*% (sv * .bbo_rnorm(length(sv))))
  } else {
    t <- numeric(nrow(X))
  }
  list(t + (P[, 1L] + .bbo_rnorm() * eta * p1))
}

.bbo_undx <- function(np, zeta, eta) function(X, pix) {
  P <- X[, pix, drop = FALSE]
  center <- .bbo_center(P)
  P <- P - center
  sd <- mean(sqrt(rowSums(P^2)))
  s <- svd(P)
  sv <- s$d / sqrt(Reduce(`+`, s$d^2))
  for (i in seq_along(sv)) sv[i] <- sd * .bbo_rnorm() * (if (isTRUE(abs(sv[i]) > 1e-8)) sv[i] * zeta else eta)
  list(drop(s$u %*% sv) + center)
}

# the mutation clock over a polynomial (pm = TRUE) or a uniform mutation
.bbo_clock <- function(lo, hi, pm, rate, eta = 50) {
  mc <- new.env(parent = emptyenv())
  mc$nextvar <- 1 + ceiling(-log(.bbo_runif()) / rate)
  mc$apply <- function(v) {
    n <- length(v)
    i <- mc$nextvar
    while (i <= n) {
      if (pm) {
        u <- 2 * .bbo_runif()
        v[i] <- if (u <= 1) v[i] + (u^(1 / (1 + eta)) - 1) * (v[i] - lo[i]) else
          v[i] + (1 - (2 - u)^(1 / (1 + eta))) * (hi[i] - v[i])
      } else {
        v[i] <- lo[i] + .bbo_runif() * (hi[i] - lo[i])
      }
      i <- i + ceiling(-log(.bbo_runif()) / rate)
    }
    mc$nextvar <- i - n
    v
  }
  mc
}

# RandomBound: a coordinate out of bounds is drawn between the bound and
# the reference point
.bbo_bound <- function(t, ref, lo, hi) {
  for (i in seq_along(t)) {
    if (is.na(t[i])) next
    if (t[i] < lo[i]) t[i] <- lo[i] + .bbo_runif() * (ref[i] - lo[i]) else
      if (t[i] > hi[i]) t[i] <- hi[i] + .bbo_runif() * (ref[i] - hi[i])
  }
  t
}

# ---------------------------------------------------------------------------
# Borg MOEA (borg_moea.jl)

.bbo_borg <- function(fn, sch, X, prm) {
  b <- new.env(parent = emptyenv())
  b$fn <- fn
  b$sch <- sch
  b$num_evals <- 0L
  D <- nrow(X)
  b$lo <- prm$lo
  b$hi <- prm$hi
  b$archive <- .bbo_archive(sch, D, prm$MaxArchiveSize)
  b$pop <- .bbo_population(X, sch$N, 1L)
  b$order <- integer(0)
  b$n_restarts <- 0L
  b$n_steps <- 0L
  b$last_restart_check <- 0L
  b$last_restart <- 0L
  b$last_wupdate <- 0L
  b$tau <- prm$tau
  b$gamma <- prm$gamma
  b$gamma_delta <- prm$gamma_delta
  b$min_popsize <- prm$PopulationSize
  b$theta <- prm$theta
  b$zeta <- prm$zeta
  b$wperiod <- prm$OperatorsUpdatePeriod
  b$rperiod <- prm$RestartCheckPeriod
  b$max_noprogress <- prm$MaxStepsWithoutEpsProgress
  b$ops <- list(
    list(np = 3L, f = .bbo_de(prm$DEX_cr, prm$DEX_f)),
    list(np = 2L, f = .bbo_sbx(prm$SBX_p, prm$SBX_eta)),
    list(np = 3L, f = .bbo_spx(3L, prm$SPX_epsilon)),
    list(np = 2L, f = .bbo_pcx(2L, prm$PCX_zeta, prm$PCX_eta)),
    list(np = 3L, f = .bbo_pcx(3L, prm$PCX_zeta, prm$PCX_eta)),
    list(np = 2L, f = .bbo_undx(2L, prm$UNDX_zeta, prm$UNDX_eta)),
    list(np = 3L, f = .bbo_undx(3L, prm$UNDX_zeta, prm$UNDX_eta)))
  b$weights <- rep(1 / 7, 7L)
  b$mutators <- list(.bbo_clock(b$lo, b$hi, TRUE, 1 / D, prm$PM_eta), .bbo_clock(b$lo, b$hi, FALSE, 1 / D))
  b$select_size <- ceiling(b$tau * b$pop$popsize)
  b
}

# fitness(params, ::ProblemEvaluator): evaluate and add to the archive
.bbo_eval <- function(b, x, tag) {
  fit <- .bbo_indexed(b$fn(x), b$sch)
  b$num_evals <- b$num_evals + 1L
  .bbo_add_candidate(b$archive, fit, x, tag, b$num_evals)
  fit
}

.bbo_step <- function(b) {
  pop <- b$pop
  a <- b$archive
  if (b$n_steps == 0L) {
    # the fitness of every member of the population
    for (i in seq_len(pop$popsize))
      if (anyNA(pop$FO[i, ])) .bbo_pop_set(pop, i, NULL, .bbo_eval(b, pop$X[, i], 0L))
  }
  b$n_steps <- b$n_steps + 1L
  if (b$n_steps >= b$last_restart_check + b$rperiod) {
    b$last_restart_check <- b$n_steps
    na <- length(a$FD)
    if ((na > 0 && abs(pop$popsize - b$gamma * na) >= b$gamma_delta * na) ||
        .bbo_noprogress(a, TRUE) >= b$max_noprogress) .bbo_restart(b)
  }
  if (b$n_steps >= b$last_wupdate + b$wperiod) {
    cnt <- .bbo_tagcounts(a, b$theta, length(b$ops))
    b$weights <- (cnt + b$zeta) / (sum(cnt) + length(b$ops) * b$zeta)
    b$last_wupdate <- b$n_steps
  }
  k <- .bbo_wsample(b$weights)
  .bbo_recombine(b, k)
}

.bbo_select <- function(b, n_tours) {
  pop <- b$pop
  size <- b$select_size
  cand <- sample.int(pop$popsize, min(pop$popsize, size * n_tours))
  minim <- b$sch$is_minimizing
  vapply(seq_len(n_tours), function(i) {
    tc <- cand[(i - 1L) * size + seq_len(size)]
    w <- tc[1L]
    for (c in tc[-1L])
      if (.bbo_hat(.bbo_pop_fit(pop, c), .bbo_pop_fit(pop, w), minim) == -1) w <- c
    w
  }, 0L)
}

.bbo_recombine <- function(b, k) {
  pop <- b$pop
  a <- b$archive
  op <- b$ops[[k]]
  na <- length(a$FD)
  pix <- .bbo_select(b, if (na == 0L) op$np else op$np - 1L)
  if (na > 0L) {
    # one parent from the archive, in the transient member (BlackBoxOptim's
    # transient_range starts one column early, at the last persistent member)
    ix <- pop$popsize + 1L
    e <- .bbo_front_elem(a, .bbo_randint(na))
    .bbo_pop_set(pop, ix, e$x, e$fit)
    pix <- c(pix, ix)
  }
  children <- op$f(pop$X, pix)
  ref <- pop$X[, pix[1L]]
  children <- lapply(children, .bbo_bound, ref = ref, lo = b$lo, hi = b$hi)
  fits <- lapply(children, function(x) .bbo_eval(b, x, k))
  for (j in seq_along(children)) .bbo_process(b, children[[j]], fits[[j]])
}

# process_candidate!: the child replaces the first member it dominates, or
# a random checked member when it is dominated by none of them
.bbo_process <- function(b, x, fit) {
  pop <- b$pop
  ps <- pop$popsize
  if (length(b$order) != ps) b$order <- sample.int(ps)
  lo <- min(ps, 2L * b$select_size)
  n_checks <- lo - 1L + .bbo_randint(ps - lo + 1L)
  ord <- b$order
  i <- seq_len(n_checks)
  js <- i + floor(.bbo_runif(n_checks) * (ps - i + 1L))
  swap <- function(o, m) {
    for (k in seq_len(m)) {
      j <- js[k]
      t <- o[j]; o[j] <- o[k]; o[k] <- t
    }
    o
  }
  o <- swap(ord, n_checks)
  checked <- o[i]
  comps <- .bbo_hat_rows(fit$index, fit$dist, pop$FI[checked, , drop = FALSE], pop$FD[checked],
                         b$sch$is_minimizing)
  first <- which(comps != 0)[1L]
  if (is.na(first)) {
    b$order <- o
    .bbo_pop_set(pop, o[.bbo_randint(n_checks)], x, fit)
  } else {
    b$order <- swap(ord, first)
    if (comps[first] < 0) .bbo_pop_set(pop, checked[first], x, fit)
  }
}

.bbo_restart <- function(b) {
  a <- b$archive
  pop <- b$pop
  a$n_restarts <- a$n_restarts + 1L
  a$last_restart <- a$num_candidates
  na <- length(a$FD)
  ps <- max(b$min_popsize, ceiling(b$gamma * na))
  # resize!: keep the first members
  keep <- seq_len(min(ps, pop$popsize))
  n <- ps + pop$ntransient
  X <- matrix(NaN, nrow(pop$X), n)
  X[, keep] <- pop$X[, keep]
  pop$X <- X
  grow <- function(M) rbind(M[keep, , drop = FALSE], matrix(NaN, n - length(keep), ncol(M)))
  pop$FO <- grow(pop$FO)
  pop$FI <- grow(pop$FI)
  pop$FD <- c(pop$FD[keep], numeric(n - length(keep)))
  pop$FA <- c(pop$FA[keep], rep(NaN, n - length(keep)))
  pop$popsize <- ps
  nd <- min(na, ps)
  dix <- sort(sample.int(na, nd))
  for (j in seq_len(nd)) {
    e <- .bbo_front_elem(a, dix[j])
    .bbo_pop_set(pop, j, e$x, e$fit)
  }
  # the rest of the population: mutated delegates
  for (ix in nd + seq_len(ps - nd)) {
    x <- pop$X[, .bbo_randint(nd)]
    x <- b$mutators[[.bbo_wsample(c(0.75, 0.25))]]$apply(x)
    x <- .bbo_bound(x, pop$X[, .bbo_randint(nd)], b$lo, b$hi)
    .bbo_pop_set(pop, ix, x, .bbo_eval(b, x, 0L))
  }
  b$select_size <- min(ps, max(3, ceiling(b$tau * ps)))
  b$last_restart <- b$n_steps
  b$n_restarts <- b$n_restarts + 1L
}

# ---------------------------------------------------------------------------
# bbsetup, bboptimize and the OptRunController

.bbo_defaults <- list(
  MaxTime = 0, MaxFuncEvals = 0L, MaxSteps = 10000L, MaxStepsWithoutProgress = 10000L,
  MinDeltaFitnessTolerance = 0x1.dee7a4ad4b81fp-167, FitnessTolerance = 1e-8, MaxNumStepsWithoutFuncEvals = 100L,
  PopulationSize = 50L,
  # Borg and its archive
  MaxArchiveSize = 100000L, epsilon = 0.1, tau = 0.02, gamma = 4, gamma_delta = 0.25, theta = 0.9, zeta = 1,
  RestartCheckPeriod = 1000L, OperatorsUpdatePeriod = 100L, MaxStepsWithoutEpsProgress = 100L,
  # the operators
  DEX_f = 0.6, DEX_cr = 0.7, SBX_p = 0.2, SBX_eta = 3, SPX_epsilon = 1.1, PCX_zeta = 0.1, PCX_eta = 0.5,
  UNDX_zeta = 0.5, UNDX_eta = 0.1, PM_eta = 50)

.bbo_optimize <- function(fn, params, x0, callback = NULL) {
  explicit <- names(params)
  prm <- utils::modifyList(.bbo_defaults, params)
  lo <- prm$SearchRange$lower
  hi <- prm$SearchRange$upper
  for (i in seq_along(lo)) if (lo[i] > hi[i])
    stop("at index ", i, " of search space, the low bound (", .bbo_jlnum(lo[i]), ") exceeds the high bound (",
         .bbo_jlnum(hi[i]), ")", call. = FALSE)
  D <- length(lo)
  # setup_problem: one evaluation at a random point checks the objectives
  ind <- .bbo_runif(D) * (hi - lo) + lo
  res <- fn(ind)
  fs <- prm$FitnessScheme %||% ParetoFitnessScheme(length(res))
  if (inherits(fs, "TupleFitnessScheme") && (!is.numeric(res) || length(res) != fs$N))
    stop("The supplied fitness function does NOT return the expected fitness type NTuple{", fs$N,
         ", Float64}when called with a potential solution (when called with ", .qd_jlvec(ind), " it returned ",
         paste(format(res), collapse = ", "), " so we cannot optimize it!", call. = FALSE)
  # check_valid!
  if (!is.numeric(prm$MaxTime) || prm$MaxTime < 0) stop("MaxTime parameter must be a non-negative number",
                                                        call. = FALSE)
  if (prm$MaxTime > 0) {
    if (!"MaxFuncEvals" %in% explicit) prm$MaxFuncEvals <- 0L
    if (!"MaxSteps" %in% explicit) prm$MaxSteps <- 0L
  }
  if (!is.numeric(prm$MaxFuncEvals) || prm$MaxFuncEvals < 0 || prm$MaxFuncEvals != round(prm$MaxFuncEvals))
    stop("MaxFuncEvals parameter MUST be a non-negative integer", call. = FALSE)
  if (prm$MaxFuncEvals > 0 && !"MaxSteps" %in% explicit) prm$MaxSteps <- 0L
  if (!is.numeric(prm$MaxSteps) || prm$MaxSteps < 0)
    stop("The number of iterations (MaxSteps) MUST be a non-negative number", call. = FALSE)
  if (prm$PopulationSize < 2) stop("The population size MUST be at least 2", call. = FALSE)
  if (!identical(prm$Method, "borg_moea"))
    stop("The method specified, ", prm$Method, ", is NOT among the valid methods", call. = FALSE)
  # borg_moea(problem, options)
  sch <- .bbo_scheme(fs, prm$epsilon)
  X <- prm[["Population"]] %||% .bbo_lhs(lo, hi, prm$PopulationSize + 1L)
  prm$lo <- lo
  prm$hi <- hi
  b <- .bbo_borg(fn, sch, X, prm)
  if (length(x0)) {
    if (length(x0) != D || any(x0 < lo | x0 > hi))
      stop("Starting point ", .qd_jlvec(x0), " is not in the search space/range.", call. = FALSE)
    i <- .bbo_randint(b$pop$popsize)
    b$pop$X[, i] <- x0
    .bbo_pop_set(b$pop, i, NULL, .bbo_nafitness(sch$N))
  }
  # run!(::OptRunController)
  b$shutdown <- FALSE
  b$num_steps <- 0L
  start <- proc.time()[["elapsed"]]
  last_fe <- 0L
  nswf <- 0L
  reason <- ""
  while (!nzchar(reason)) {
    .bbo_step(b)
    b$num_steps <- b$num_steps + 1L
    if (b$num_evals == last_fe) nswf <- nswf + 1L else nswf <- 0L
    last_fe <- b$num_evals
    if (!is.null(callback)) callback(b)
    if (prm$MaxTime > 0) elapsed <- proc.time()[["elapsed"]] - start
    a <- b$archive
    reason <- if (b$shutdown) "Run explicitly stopped via shutdown method" else
      if (prm$MaxTime > 0 && elapsed > prm$MaxTime) paste0("Max time (", .bbo_jlnum(prm$MaxTime), " s) reached") else
      if (prm$MaxFuncEvals > 0 && b$num_evals > prm$MaxFuncEvals)
        paste0("Max number of function evaluations (", prm$MaxFuncEvals, ") reached") else
      if (prm$MaxNumStepsWithoutFuncEvals > 0 && nswf > prm$MaxNumStepsWithoutFuncEvals)
        paste0("Too many steps (", nswf, ") without any function evaluations (probably search has converged)") else
      if (prm$MaxSteps > 0 && b$num_steps > prm$MaxSteps)
        paste0("Max number of steps (", prm$MaxSteps, ") reached") else
      if (prm$MaxStepsWithoutProgress > 0 && .bbo_noprogress(a) > prm$MaxStepsWithoutProgress)
        paste0("No epsilon-progress for more than ", prm$MaxStepsWithoutProgress, " iterations") else
      if (a$n_oversize_inserts >= 10)
        paste0("Pareto frontier size (", length(a$FD), ") exceeded maximum (", a$max_size, ")") else ""
  }
  elapsed <- proc.time()[["elapsed"]] - start
  a <- b$archive
  frontier <- list(X = a$X, F = a$FO)
  list(method = "borg_moea", stop_reason = reason, iterations = b$num_steps, elapsed_time = elapsed,
       f_calls = b$num_evals, best_candidate = a$best$x, best_fitness = a$best$fit$orig, frontier = frontier,
       parameters = params, state = b)
}
