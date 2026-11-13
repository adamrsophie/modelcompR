# OptimizationNLopt: the NLopt algorithm identifiers (Algorithm), the
# optimizer object (Opt) and the NLopt module of algorithm constructors, after
# lib/OptimizationNLopt/src/OptimizationNLopt.jl and NLopt.jl, and the driver
# of the NLopt C library (src/api/optimize.c of NLopt 2.11.0: the dimension
# elimination for fixed variables, the memoization of the best point and the
# dispatch to the algorithms). The algorithms are native R ports of the NLopt
# C sources, one file per family (R/optimization-nlopt-*.R).
#
# NLopt is Copyright (c) 2007-2014 Massachusetts Institute of Technology,
# released under the MIT license (see R/optimization-nlopt-util.R).

# the algorithm enumeration of NLopt (nlopt_algorithm, values 0 to 43)
.nlopt_algorithms <- c(
  "GN_DIRECT", "GN_DIRECT_L", "GN_DIRECT_L_RAND", "GN_DIRECT_NOSCAL", "GN_DIRECT_L_NOSCAL",
  "GN_DIRECT_L_RAND_NOSCAL", "GN_ORIG_DIRECT", "GN_ORIG_DIRECT_L", "GD_STOGO", "GD_STOGO_RAND",
  "LD_LBFGS_NOCEDAL", "LD_LBFGS", "LN_PRAXIS", "LD_VAR1", "LD_VAR2", "LD_TNEWTON", "LD_TNEWTON_RESTART",
  "LD_TNEWTON_PRECOND", "LD_TNEWTON_PRECOND_RESTART", "GN_CRS2_LM", "GN_MLSL", "GD_MLSL", "GN_MLSL_LDS",
  "GD_MLSL_LDS", "LD_MMA", "LN_COBYLA", "LN_NEWUOA", "LN_NEWUOA_BOUND", "LN_NELDERMEAD", "LN_SBPLX",
  "LN_AUGLAG", "LD_AUGLAG", "LN_AUGLAG_EQ", "LD_AUGLAG_EQ", "LN_BOBYQA", "GN_ISRES", "AUGLAG", "AUGLAG_EQ",
  "G_MLSL", "G_MLSL_LDS", "LD_SLSQP", "LD_CCSAQ", "GN_ESCH", "GN_AGS")

# the identifiers whose NLopt algorithm is ported (the others raise on solve)
.nlopt_ported <- c("LN_NELDERMEAD", "LN_SBPLX", "LD_LBFGS_NOCEDAL", "GN_CRS2_LM", "GN_ESCH", "GN_ISRES",
                   "LD_MMA", "LD_CCSAQ", "AUGLAG", "AUGLAG_EQ", "LN_AUGLAG", "LN_AUGLAG_EQ", "LD_AUGLAG",
                   "LD_AUGLAG_EQ", "G_MLSL", "G_MLSL_LDS", "GN_MLSL", "GD_MLSL", "GN_MLSL_LDS", "GD_MLSL_LDS",
                   "GN_DIRECT", "GN_DIRECT_L", "GN_DIRECT_L_RAND", "GN_DIRECT_NOSCAL", "GN_DIRECT_L_NOSCAL",
                   "GN_DIRECT_L_RAND_NOSCAL", "GN_ORIG_DIRECT", "GN_ORIG_DIRECT_L", "LN_PRAXIS", "GN_AGS", "GD_STOGO", "GD_STOGO_RAND",
                   "LD_LBFGS", "LD_VAR1", "LD_VAR2", "LD_TNEWTON", "LD_TNEWTON_RESTART", "LD_TNEWTON_PRECOND",
                   "LD_TNEWTON_PRECOND_RESTART", "LN_COBYLA", "LN_NEWUOA", "LN_NEWUOA_BOUND", "LN_BOBYQA",
                   "LD_SLSQP")

#' NLopt algorithms and optimizers
#'
#' The interface of OptimizationNLopt.jl: `Algorithm(name)` is an NLopt
#' algorithm identifier such as `"LD_LBFGS"` (also given by its number in
#' the NLopt enumeration, from 0), `NLopt` holds a constructor for each
#' identifier (`NLopt$LN_COBYLA()`) and NLopt's `srand(seed)` and
#' `srand_time()`, and `Opt(algorithm, n)` is an optimizer for `n` unknowns
#' whose settings (`xtol_rel`, `xtol_abs`, `ftol_rel`, `ftol_abs`,
#' `maxeval`, `maxtime`, `stopval`, `initial_step`, `local_optimizer`,
#' `lower_bounds`, `upper_bounds`, `population`, `vector_storage` and the
#' algorithm parameters `params`) are fields set with `$<-`. An algorithm
#' or an `Opt` is passed to [solve()].
#'
#' The algorithms are native R ports of the C code of NLopt 2.11.0, so they
#' take the steps of NLopt: `LN_NELDERMEAD` (the Nelder Mead simplex with
#' bounds), `LN_SBPLX` (Rowan's subplex), `GN_CRS2_LM` (controlled random
#' search), `GN_ESCH` (an evolutionary strategy), `GN_ISRES` (stochastic
#' ranking evolution), `LD_MMA` and `LD_CCSAQ` (conservative convex
#' separable approximations), the augmented Lagrangian family (`AUGLAG`,
#' `AUGLAG_EQ`, `LN_AUGLAG`, `LN_AUGLAG_EQ`, `LD_AUGLAG`, `LD_AUGLAG_EQ`) and
#' the multi level single linkage family (`G_MLSL`, `G_MLSL_LDS`, `GN_MLSL`,
#' `GD_MLSL`, `GN_MLSL_LDS`, `GD_MLSL_LDS`), the DIRECT family (`GN_DIRECT`,
#' `GN_DIRECT_L`, `GN_DIRECT_L_RAND` and their `_NOSCAL` variants, and
#' Gablonsky's `GN_ORIG_DIRECT` and `GN_ORIG_DIRECT_L`), `LN_PRAXIS` (Brent's
#' principal axis method), `GD_STOGO` and `GD_STOGO_RAND` (branch and bound
#' with trust region BFGS searches), `GN_AGS` (the index method on a Peano
#' evolvent), Luksan's limited memory and
#' truncated Newton methods (`LD_LBFGS`, `LD_VAR1`, `LD_VAR2`, `LD_TNEWTON`,
#' `LD_TNEWTON_RESTART`, `LD_TNEWTON_PRECOND`, `LD_TNEWTON_PRECOND_RESTART`),
#' Powell's `LN_COBYLA`, `LN_NEWUOA`, `LN_NEWUOA_BOUND` and `LN_BOBYQA`, and
#' Kraft's `LD_SLSQP`. `LD_LBFGS_NOCEDAL` has no implementation in NLopt and
#' returns `INVALID_ARGS` (a `Failure`). The
#' stochastic algorithms draw from NLopt's Mersenne Twister, seeded by
#' `NLopt$srand(seed)` (from the clock otherwise).
#'
#' The constraints `lcons <= cons(u) <= ucons` become NLopt's equality
#' constraints `cons(u) - lcons = 0` (when `lcons == ucons`) and inequality
#' constraints `cons(u) - ucons <= 0` and `lcons - cons(u) <= 0`, with the
#' tolerance `cons_tol` (default `1e-6`). `maxiters` sets `maxeval`,
#' `maxtime` sets `maxtime`, `abstol` sets `ftol_abs` and `reltol` sets
#' `ftol_rel`; `local_method`, `local_maxiters`, `local_maxtime` and
#' `local_options` configure the local optimizer of the AUGLAG and MLSL
#' algorithms, and any other option of [solve()] becomes an algorithm
#' parameter (NLopt's `nlopt_set_param`). The `stopval` of an `Opt` stops
#' the search at the first point at least that good, and `maxtime` after
#' that many seconds (return codes `STOPVAL_REACHED`, a `Success`, and
#' `MAXTIME_REACHED`, `MaxTime`).
#'
#' @param name An algorithm name (with or without the `NLOPT_` prefix) or
#'   its number.
#' @param algorithm An `Algorithm` or its name.
#' @param n The number of unknowns.
#' @return `Algorithm()` an algorithm object, `Opt()` an optimizer object
#'   (an environment).
#' @examples
#' prob <- OptimizationProblem(function(u, p) (u[1] - 1) * (u[1] - 1) + 2 * (u[2] + 1) * (u[2] + 1),
#'                             c(0, 0), numeric())
#' solve(prob, NLopt$LN_NELDERMEAD(), maxiters = 400)$u
#' opt <- Opt("LN_SBPLX", 2)
#' opt$xtol_rel <- 1e-8
#' solve(prob, opt)$u
#' @name nlopt_optimizers
NULL

#' @rdname nlopt_optimizers
#' @export
Algorithm <- function(name) {
  if (is.numeric(name)) name <- .nlopt_algorithms[as.integer(name) + 1L]
  name <- sub("^NLOPT_", "", toupper(as.character(name)))
  if (length(name) != 1L || is.na(name) || !(name %in% .nlopt_algorithms))
    stop("unknown NLopt algorithm ", name, call. = FALSE)
  structure(list(name = name, value = match(name, .nlopt_algorithms) - 1L,
                 traits = list(callback = TRUE, bounds = TRUE,
                               cons = grepl("AUGLAG|CCSA|MMA|COBYLA|ISRES|AGS|ORIG_DIRECT|SLSQP", name)),
                 solve = .nlopt_solve),
            class = c("NLoptAlgorithm", "OptimizationAlgorithm"))
}

#' @export
print.NLoptAlgorithm <- function(x, ...) {
  cat(x$name, "::Algorithm = ", x$value, "\n", sep = "")
  invisible(x)
}

#' @rdname nlopt_optimizers
#' @export
NLopt <- c(stats::setNames(lapply(.nlopt_algorithms, function(nm) {
  force(nm)
  function() Algorithm(nm)
}), .nlopt_algorithms), list(srand = function(seed) invisible(.nlopt_srand(seed)),
                             srand_time = function() invisible(.nlopt_srand_time())))

#' @rdname nlopt_optimizers
#' @export
Opt <- function(algorithm, n) {
  alg <- if (inherits(algorithm, "NLoptAlgorithm")) algorithm else Algorithm(algorithm)
  e <- new.env(parent = emptyenv())
  e$algorithm <- alg
  e$n <- as.integer(n)
  e$lower_bounds <- rep(-Inf, n)
  e$upper_bounds <- rep(Inf, n)
  e$stopval <- -Inf
  e$ftol_rel <- 0
  e$ftol_abs <- 0
  e$xtol_rel <- 0
  e$xtol_abs <- rep(0, n)
  e$x_weights <- NULL
  e$maxeval <- 0L
  e$maxtime <- 0
  e$initial_step <- NULL
  e$local_optimizer <- NULL
  e$population <- 0L
  e$vector_storage <- 0L
  e$params <- list()
  e$numevals <- 0L
  e$traits <- alg$traits
  e$solve <- .nlopt_solve
  class(e) <- c("NLoptOpt", "OptimizationAlgorithm")
  e
}

#' @export
print.NLoptOpt <- function(x, ...) {
  cat("Opt(", x$algorithm$name, ", ", x$n, ")\n", sep = "")
  invisible(x)
}

# nlopt_set_lower_bounds / nlopt_set_upper_bounds: a bound pair closer than
# a subnormal number collapses to a fixed variable
.nlopt_set_bounds <- function(e, lb = NULL, ub = NULL) {
  if (!is.null(lb)) {
    lb <- rep_len(as.numeric(lb), e$n)
    tiny <- lb < e$upper_bounds & vapply(e$upper_bounds - lb, .nlopt_istiny, TRUE)
    lb[tiny] <- e$upper_bounds[tiny]
    e$lower_bounds <- lb
  }
  if (!is.null(ub)) {
    ub <- rep_len(as.numeric(ub), e$n)
    tiny <- e$lower_bounds < ub & vapply(ub - e$lower_bounds, .nlopt_istiny, TRUE)
    ub[tiny] <- e$lower_bounds[tiny]
    e$upper_bounds <- ub
  }
  invisible(e)
}

# the plain settings of an Opt (the optimizer description of the driver)
.nlopt_settings <- function(e) {
  lo <- e$local_optimizer
  list(algorithm = e$algorithm$name, n = e$n, lb = e$lower_bounds, ub = e$upper_bounds, stopval = e$stopval,
       ftol_rel = e$ftol_rel, ftol_abs = e$ftol_abs, xtol_rel = e$xtol_rel, xtol_abs = e$xtol_abs,
       x_weights = e$x_weights, maxeval = as.integer(e$maxeval), maxtime = e$maxtime, dx = e$initial_step,
       local_opt = if (!is.null(lo)) .nlopt_settings(lo), params = e$params %||% list(),
       population = as.integer(e$population %||% 0L), vector_storage = as.integer(e$vector_storage %||% 0L))
}

.nlopt_get_param <- function(o, name, default) {
  v <- o$params[[name]]
  if (is.null(v)) default else as.numeric(v)
}

.nlopt_solve <- function(a, alg, opts) {
  is_opt <- inherits(alg, "NLoptOpt")
  name <- if (is_opt) alg$algorithm$name else alg$name
  fn <- a$f
  p <- a$p
  x0 <- as.numeric(a$u0)
  n <- length(x0)
  if (substr(name, 2L, 2L) != "N" && !is.function(fn$grad))
    stop(IncompatibleOptimizerError(paste0(
      "The NLopt algorithm ", name, " requires gradients, but no gradient function is available. Please use ",
      "`OptimizationFunction` with an automatic differentiation backend, or provide gradients manually via ",
      "the `grad` kwarg.")))
  if (is_opt && alg$n != n)
    stop("Passed NLopt.Opt optimization dimension does not match OptimizationProblem dimension.", call. = FALSE)
  if (!(name %in% .nlopt_ported))
    stop("the NLopt algorithm ", name, " is not yet ported to native R in modelcompR", call. = FALSE)
  e <- if (is_opt) alg else Opt(alg, n)
  # the objective, already in the minimizing sense (solve() negated it for
  # MaxSense); NLopt.jl's max_objective! and nlopt_optimize flip stopval
  f <- function(x) as.numeric(fn$f(x, p))[1L]
  fg <- if (is.function(fn$grad)) function(x) list(val = f(x), grad = as.numeric(fn$grad(x, p)))
  maximize <- identical(unclass(a$sense %||% MinSense), "MaxSense")
  # constraints: an inequality group and an equality group, as
  # OptimizationNLopt registers them
  m <- length(a$lcons %||% numeric())
  if (m && !alg$traits$cons)
    stop("The algorithm ", name, " does not support constraints", call. = FALSE)
  fc <- list(); h <- list()
  if (m && is.function(fn$cons)) {
    lc <- rep_len(as.numeric(a$lcons), m); uc <- rep_len(as.numeric(a$ucons %||% a$lcons), m)
    cv <- function(x) as.numeric(fn$cons(x, p))
    cj <- function(x) {
      J <- if (is.function(fn$cons_j)) fn$cons_j(x, p) else .fd_jacobian(cv, x)
      matrix(as.numeric(J), m, n)
    }
    ie <- which(lc == uc); iu <- which(lc != uc & is.finite(uc)); il <- which(lc != uc & is.finite(lc))
    tol <- opts$cons_tol %||% 1e-6
    mk <- function(vfun, gfun, k) list(m = k, f = vfun, fg = function(x) list(val = vfun(x), grad = gfun(x)),
                                       tol = rep(tol, k))
    if (length(iu) + length(il)) {
      fc <- list(mk(function(x) { v <- cv(x); c(v[iu] - uc[iu], lc[il] - v[il]) },
                    function(x) { J <- cj(x); rbind(J[iu, , drop = FALSE], -J[il, , drop = FALSE]) },
                    length(iu) + length(il)))
    }
    if (length(ie)) {
      h <- list(mk(function(x) cv(x)[ie] - lc[ie], function(x) cj(x)[ie, , drop = FALSE], length(ie)))
    }
  }
  # __map_optimizer_args!
  if (name %in% c("LN_AUGLAG", "LD_AUGLAG", "AUGLAG") && is.null(opts$local_method))
    stop("NLopt.", name, " requires a local optimization method. Please specify a local_method, e.g., ",
         "solve(prob, NLopt.", name, "(); local_method = NLopt.LN_NELDERMEAD())", call. = FALSE)
  if (!is.null(opts$local_method)) {
    lm <- opts$local_method
    if (inherits(lm, "NLoptOpt")) {
      if (lm$n != n)
        stop("Passed local NLopt.Opt optimization dimension does not match OptimizationProblem dimension.",
             call. = FALSE)
      local <- lm
    } else local <- Opt(lm, n)
    # as OptimizationNLopt, local_options are set on the main optimizer
    for (nm in names(opts$local_options)) e$params[[nm]] <- opts$local_options[[nm]]
    if (!is.null(opts$local_maxiters)) local$maxeval <- as.integer(opts$local_maxiters)
    if (!is.null(opts$local_maxtime)) local$maxtime <- opts$local_maxtime
    e$local_optimizer <- local
  }
  known <- c("maxiters", "maxtime", "abstol", "reltol", "local_method", "local_maxiters", "local_maxtime",
             "local_options", "cons_tol", "callback", "progress", "verbose")
  for (nm in setdiff(names(opts), known)) e$params[[nm]] <- opts[[nm]]
  .nlopt_set_bounds(e, lb = a$lb, ub = a$ub)
  if (!is.null(opts$maxiters)) e$maxeval <- as.integer(opts$maxiters)
  if (!is.null(opts$maxtime)) e$maxtime <- as.numeric(opts$maxtime)
  if (!is.null(opts$abstol)) e$ftol_abs <- opts$abstol
  if (!is.null(opts$reltol)) e$ftol_rel <- opts$reltol
  o <- .nlopt_settings(e)
  if (maximize) o$stopval <- if (is.infinite(o$stopval) && o$stopval < 0) -Inf else -o$stopval
  o$f <- f; o$fg <- fg; o$fc <- fc; o$h <- h
  t0 <- .nlopt_seconds()
  res <- .nlopt_optimize(o, x0)
  e$numevals <- res$numevals
  retcode <- switch(res$ret, SUCCESS = , STOPVAL_REACHED = , FTOL_REACHED = , XTOL_REACHED = "Success",
                    MAXEVAL_REACHED = "MaxIters", MAXTIME_REACHED = "MaxTime", "Failure")
  sol <- .opt_solution(a, alg, res$x, res$minf, retcode, list(ret = res$ret, numevals = res$numevals))
  .opt_set_fields(sol, stats = OptimizationStats(time = .nlopt_seconds() - t0))
}

# ---------------------------------------------------------------------------
# The driver (optimize.c)

# nlopt_set_default_initial_step
.nlopt_default_initial_step <- function(lb, ub, x) {
  vapply(seq_along(x), function(i) {
    step <- Inf
    if (!.nlopt_isinf(ub[i]) && !.nlopt_isinf(lb[i]) && (ub[i] - lb[i]) * 0.25 < step && ub[i] > lb[i])
      step <- (ub[i] - lb[i]) * 0.25
    if (!.nlopt_isinf(ub[i]) && ub[i] - x[i] < step && ub[i] > x[i]) step <- (ub[i] - x[i]) * 0.75
    if (!.nlopt_isinf(lb[i]) && x[i] - lb[i] < step && x[i] > lb[i]) step <- (x[i] - lb[i]) * 0.75
    if (.nlopt_isinf(step)) {
      if (!.nlopt_isinf(ub[i]) && abs(ub[i] - x[i]) < abs(step)) step <- (ub[i] - x[i]) * 1.1
      if (!.nlopt_isinf(lb[i]) && abs(x[i] - lb[i]) < abs(step)) step <- (x[i] - lb[i]) * 1.1
    }
    if (.nlopt_isinf(step) || .nlopt_istiny(step)) step <- x[i]
    if (.nlopt_isinf(step) || step == 0) step <- 1
    step
  }, 0)
}

# the initial step of an optimizer at x: its own, else the default one
.nlopt_dx <- function(o, x) if (!is.null(o$dx)) rep_len(as.numeric(o$dx), length(x)) else
  .nlopt_default_initial_step(o$lb, o$ub, x)

.nlopt_elimdim_algs <- c("GN_DIRECT", "GN_DIRECT_L", "GN_DIRECT_L_RAND", "GN_DIRECT_NOSCAL", "GN_DIRECT_L_NOSCAL",
                         "GN_DIRECT_L_RAND_NOSCAL", "GN_ORIG_DIRECT", "GN_ORIG_DIRECT_L", "GN_CRS2_LM", "LN_PRAXIS",
                         "LN_COBYLA", "LN_NEWUOA", "LN_NEWUOA_BOUND", "LN_BOBYQA", "LN_NELDERMEAD", "LN_SBPLX",
                         "GN_ISRES", "GN_ESCH", "GN_AGS", "GD_STOGO", "GD_STOGO_RAND")
.nlopt_memoize_algs <- c("LN_COBYLA", "LD_TNEWTON", "LD_TNEWTON_RESTART", "LD_TNEWTON_PRECOND",
                         "LD_TNEWTON_PRECOND_RESTART")

# nlopt_optimize for a minimization: `o` holds the settings (see
# .nlopt_settings) and the objective `f(x)` (with `fg(x)` giving
# list(val, grad)), the inequality constraints `fc` and the equality
# constraints `h` (lists of list(m, f, fg, tol)). Returns list(ret, x,
# minf, numevals).
.nlopt_optimize <- function(o, x) {
  x <- as.numeric(x)
  o$force_stop <- o$force_stop %||% local({ e <- new.env(parent = emptyenv()); e$value <- FALSE; e })
  memo <- NULL
  if (!length(o$fc) && !length(o$h) && o$algorithm %in% .nlopt_memoize_algs) {
    memo <- new.env(parent = emptyenv())
    memo$minf <- .Machine$double.xmax; memo$bestx <- x
    keep <- function(xx, val) {
      if (all(xx >= o$lb) && all(xx <= o$ub) && val < memo$minf) { memo$minf <- val; memo$bestx <- xx }
      val
    }
    f0 <- o$f; fg0 <- o$fg
    o$f <- function(xx) keep(xx, f0(xx))
    if (is.function(fg0)) o$fg <- function(xx) { r <- fg0(xx); keep(xx, r$val); r }
  }
  fixed <- o$lb == o$ub
  if (any(fixed) && o$algorithm %in% .nlopt_elimdim_algs) {
    # elimdim_create: the free variables only
    full <- function(x0) { xx <- o$lb; xx[!fixed] <- x0; xx }
    oe <- o
    oe$n <- sum(!fixed)
    oe$lb <- o$lb[!fixed]; oe$ub <- o$ub[!fixed]
    if (!is.null(o$xtol_abs)) oe$xtol_abs <- rep_len(o$xtol_abs, o$n)[!fixed]
    if (!is.null(o$dx)) oe$dx <- rep_len(o$dx, o$n)[!fixed]
    oe$f <- function(x0) o$f(full(x0))
    if (is.function(o$fg)) oe$fg <- function(x0) { r <- o$fg(full(x0)); r$grad <- r$grad[!fixed]; r }
    elim <- function(c) {
      c$f <- local({ g <- c$f; function(x0) g(full(x0)) })
      c$fg <- local({ g <- c$fg; function(x0) { r <- g(full(x0)); r$grad <- r$grad[, !fixed, drop = FALSE]; r } })
      c
    }
    oe$fc <- lapply(o$fc, elim); oe$h <- lapply(o$h, elim)
    res <- .nlopt_optimize_(oe, x[!fixed])
    res$x <- full(res$x)
  } else {
    res <- .nlopt_optimize_(o, x)
  }
  if (!is.null(memo)) { res$x <- memo$bestx; res$minf <- memo$minf }
  res
}

.nlopt_finite_domain <- function(lb, ub) !any(vapply(ub - lb, .nlopt_isinf, TRUE))

# nlopt_optimize_: the checks and the dispatch to the algorithms
.nlopt_optimize_ <- function(o, x) {
  n <- length(x)
  out <- function(ret, xx = x, minf = Inf, nevals = 0L) list(ret = ret, x = xx, minf = minf, numevals = nevals)
  if (n == 0L) return(out("SUCCESS", x, o$f(x)))
  .nlopt_srand_time_default()
  if (any(o$lb > o$ub | x < o$lb | x > o$ub)) return(out("INVALID_ARGS"))
  stop <- .nlopt_stopping(n, minf_max = o$stopval, ftol_rel = o$ftol_rel, ftol_abs = o$ftol_abs,
                          xtol_rel = o$xtol_rel, xtol_abs = o$xtol_abs, x_weights = o$x_weights,
                          maxeval = o$maxeval, maxtime = o$maxtime, force_stop = o$force_stop)
  alg <- o$algorithm
  pop <- as.integer(o$population %||% 0L) # POP(0)
  global_ok <- .nlopt_finite_domain(o$lb, o$ub)
  bad_domain <- list(ret = "INVALID_ARGS", x = x, minf = Inf)
  res <- switch(alg,
    LN_NELDERMEAD = .nlopt_nldrmd_minimize(o$f, o$lb, o$ub, x, .nlopt_dx(o, x), stop),
    LN_SBPLX = .nlopt_sbplx_minimize(o$f, o$lb, o$ub, x, .nlopt_dx(o, x), stop),
    AUGLAG = , AUGLAG_EQ = , LN_AUGLAG = , LN_AUGLAG_EQ = , LD_AUGLAG = , LD_AUGLAG_EQ =
      .nlopt_auglag_dispatch(o, x, stop),
    G_MLSL = , G_MLSL_LDS = , GN_MLSL = , GD_MLSL = , GN_MLSL_LDS = , GD_MLSL_LDS =
      .nlopt_mlsl_dispatch(o, x, stop, pop),
    LD_MMA = .nlopt_ccsa_minimize(o, x, stop, quadratic = FALSE),
    LD_CCSAQ = .nlopt_ccsa_minimize(o, x, stop, quadratic = TRUE),
    GN_CRS2_LM = if (!global_ok) bad_domain else .nlopt_crs_minimize(o$f, o$lb, o$ub, x, stop, pop),
    GN_ESCH = if (!global_ok) bad_domain else
      .nlopt_esch_minimize(o$f, o$lb, o$ub, x, stop, pop, as.integer(floor(pop * 1.5))),
    GN_ISRES = if (!global_ok) bad_domain else .nlopt_isres_minimize(o$f, o$fc, o$h, o$lb, o$ub, x, stop, pop),
    GN_DIRECT = , GN_DIRECT_L = , GN_DIRECT_L_RAND = , GN_DIRECT_NOSCAL = , GN_DIRECT_L_NOSCAL = ,
    GN_DIRECT_L_RAND_NOSCAL = .nlopt_cdirect_dispatch(o, x, stop),
    GN_ORIG_DIRECT = , GN_ORIG_DIRECT_L = .nlopt_orig_direct_dispatch(o, x, stop),
    GN_AGS = .nlopt_ags_dispatch(o, x, stop),
    GD_STOGO = , GD_STOGO_RAND = .nlopt_stogo_dispatch(o, x, stop),
    LN_PRAXIS = .nlopt_praxis_dispatch(o, x, stop),
    LD_LBFGS = .nlopt_luksan_plis(o$fg, o$lb, o$ub, x, stop, o$vector_storage, .nlopt_get_param(o, "tolg", 0)),
    LD_VAR1 = , LD_VAR2 = .nlopt_luksan_plip(o$fg, o$lb, o$ub, x, stop, o$vector_storage,
                                            if (alg == "LD_VAR1") 1L else 2L, .nlopt_get_param(o, "tolg", 0)),
    LD_TNEWTON = , LD_TNEWTON_RESTART = , LD_TNEWTON_PRECOND = , LD_TNEWTON_PRECOND_RESTART = {
      k <- match(alg, c("LD_TNEWTON", "LD_TNEWTON_RESTART", "LD_TNEWTON_PRECOND", "LD_TNEWTON_PRECOND_RESTART")) - 1L
      .nlopt_luksan_pnet(o$fg, o$lb, o$ub, x, stop, o$vector_storage, 1L + k %% 2L, 1L + k %/% 2L,
                         .nlopt_get_param(o, "tolg", 0))
    },
    LN_COBYLA = .nlopt_cobyla_minimize(o, x, stop, .nlopt_dx(o, x)),
    LN_NEWUOA = .nlopt_newuoa_minimize(o, x, stop, FALSE),
    LN_NEWUOA_BOUND = .nlopt_newuoa_minimize(o, x, stop, TRUE),
    LN_BOBYQA = .nlopt_bobyqa(o$f, o$lb, o$ub, x, .nlopt_dx(o, x), stop),
    LD_SLSQP = .nlopt_slsqp(o, x, stop),
    list(ret = "INVALID_ARGS", x = x, minf = Inf))
  out(res$ret, res$x, res$minf, stop$nevals)
}
