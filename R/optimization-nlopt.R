# OptimizationNLopt: the NLopt algorithm identifiers (Algorithm), the
# optimizer object (Opt) and the NLopt module of algorithm constructors, after
# lib/OptimizationNLopt/src/OptimizationNLopt.jl and NLopt.jl. The NLopt C
# library is not available here: the algorithms that have a native kernel in
# modelcompR run that kernel (see .nlopt_kernels), the others raise an error.

# the algorithm enumeration of NLopt (nlopt_algorithm, values 0 to 43)
.nlopt_algorithms <- c(
  "GN_DIRECT", "GN_DIRECT_L", "GN_DIRECT_L_RAND", "GN_DIRECT_NOSCAL", "GN_DIRECT_L_NOSCAL",
  "GN_DIRECT_L_RAND_NOSCAL", "GN_ORIG_DIRECT", "GN_ORIG_DIRECT_L", "GD_STOGO", "GD_STOGO_RAND",
  "LD_LBFGS_NOCEDAL", "LD_LBFGS", "LN_PRAXIS", "LD_VAR1", "LD_VAR2", "LD_TNEWTON", "LD_TNEWTON_RESTART",
  "LD_TNEWTON_PRECOND", "LD_TNEWTON_PRECOND_RESTART", "GN_CRS2_LM", "GN_MLSL", "GD_MLSL", "GN_MLSL_LDS",
  "GD_MLSL_LDS", "LD_MMA", "LN_COBYLA", "LN_NEWUOA", "LN_NEWUOA_BOUND", "LN_NELDERMEAD", "LN_SBPLX",
  "LN_AUGLAG", "LD_AUGLAG", "LN_AUGLAG_EQ", "LD_AUGLAG_EQ", "LN_BOBYQA", "GN_ISRES", "AUGLAG", "AUGLAG_EQ",
  "G_MLSL", "G_MLSL_LDS", "LD_SLSQP", "LD_CCSAQ", "GN_ESCH", "GN_AGS")

# algorithms with a native kernel: the kernel standing in for NLopt's code
.nlopt_kernels <- c(LN_NELDERMEAD = "nelder-mead", LN_COBYLA = "cobyla", LD_LBFGS = "l-bfgs-b",
                    LD_LBFGS_NOCEDAL = "l-bfgs-b", LD_SLSQP = "slsqp", LN_BOBYQA = "bobyqa", LN_NEWUOA = "newuoa",
                    LN_NEWUOA_BOUND = "bobyqa", GN_DIRECT = "direct", GN_DIRECT_L = "direct_l",
                    GN_ORIG_DIRECT = "direct", GN_ORIG_DIRECT_L = "direct_l")

#' NLopt algorithms and optimizers
#'
#' The interface of OptimizationNLopt.jl: `Algorithm(name)` is an NLopt
#' algorithm identifier such as `"LD_LBFGS"` (also given by its number in
#' the NLopt enumeration, from 0), `NLopt` holds a constructor for each
#' identifier (`NLopt$LN_COBYLA()`), and `Opt(algorithm, n)` is an optimizer
#' for `n` unknowns whose settings (`xtol_rel`, `xtol_abs`, `ftol_rel`,
#' `ftol_abs`, `maxeval`, `maxtime`, `stopval`, `initial_step`,
#' `local_optimizer`, `lower_bounds`, `upper_bounds`) are fields set with
#' `$<-`. An algorithm or an `Opt` is passed to [solve()].
#'
#' The NLopt library is not available to modelcompR, so each algorithm runs a
#' native kernel of the package that implements the same method:
#' `LN_NELDERMEAD` the Nelder Mead simplex, `LN_COBYLA` COBYLA (PRIMA),
#' `LD_LBFGS` and `LD_LBFGS_NOCEDAL` limited memory BFGS (L-BFGS-B). Results
#' agree with NLopt at the solution, not step by step. The other identifiers
#' raise an error on [solve()] until their kernels are ported.
#'
#' The constraints `lcons <= cons(u) <= ucons` become NLopt's equality
#' constraints `cons(u) - lcons = 0` (when `lcons == ucons`) and inequality
#' constraints `cons(u) - ucons <= 0` and `lcons - cons(u) <= 0`. `maxiters`
#' sets `maxeval`, `abstol` sets `ftol_abs` and `reltol` sets `ftol_rel`.
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
#' opt <- Opt("LN_COBYLA", 2)
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
NLopt <- stats::setNames(lapply(.nlopt_algorithms, function(nm) {
  force(nm)
  function() Algorithm(nm)
}), .nlopt_algorithms)

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
  e$maxeval <- 0L
  e$maxtime <- 0
  e$initial_step <- NULL
  e$local_optimizer <- NULL
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

.nlopt_solve <- function(a, alg, opts) {
  is_opt <- inherits(alg, "NLoptOpt")
  name <- if (is_opt) alg$algorithm$name else alg$name
  fn <- a$f
  p <- a$p
  x0 <- as.numeric(a$u0)
  n <- length(x0)
  if (is_opt && alg$n != n)
    stop("Passed NLopt.Opt optimization dimension does not match OptimizationProblem dimension.", call. = FALSE)
  if (substr(name, 2L, 2L) != "N" && !is.function(fn$grad))
    stop(IncompatibleOptimizerError(paste0(
      "The NLopt algorithm ", name, " requires gradients, but no gradient function is available. Please use ",
      "`OptimizationFunction` with an automatic differentiation backend, or provide gradients manually via ",
      "the `grad` kwarg.")))
  if (name %in% c("LN_AUGLAG", "LD_AUGLAG", "AUGLAG") && is.null(opts$local_method))
    stop("NLopt.", name, " requires a local optimization method. Please specify a local_method, e.g., ",
         "solve(prob, NLopt.", name, "(); local_method = NLopt.LN_NELDERMEAD())", call. = FALSE)
  kernel <- unname(.nlopt_kernels[name])
  # only the kernels with a branch below run; the others are listed for when
  # their ports land
  if (is.na(kernel) || !(kernel %in% c("nelder-mead", "cobyla", "l-bfgs-b")))
    stop("the NLopt algorithm ", name, " is not yet ported to native R in modelcompR", call. = FALSE)
  # settings of the Opt, then the common arguments of solve()
  st <- if (is_opt) as.list.environment(alg) else list(stopval = -Inf, ftol_rel = 0, ftol_abs = 0, xtol_rel = 0,
                                           xtol_abs = rep(0, n), maxeval = 0L, maxtime = 0)
  if (!is.null(opts$maxiters)) st$maxeval <- as.integer(opts$maxiters)
  if (!is.null(opts$abstol)) st$ftol_abs <- opts$abstol
  if (!is.null(opts$reltol)) st$ftol_rel <- opts$reltol
  for (nm in intersect(names(opts), c("xtol_rel", "xtol_abs", "ftol_rel", "ftol_abs", "stopval", "initial_step")))
    st[[nm]] <- opts[[nm]]
  lb <- if (!is.null(a$lb)) as.numeric(a$lb) else if (is_opt) alg$lower_bounds else rep(-Inf, n)
  ub <- if (!is.null(a$ub)) as.numeric(a$ub) else if (is_opt) alg$upper_bounds else rep(Inf, n)
  f <- function(x) as.numeric(fn$f(x, p))[1L]
  grad <- if (is.function(fn$grad)) function(x) as.numeric(fn$grad(x, p))
  m <- length(a$lcons %||% numeric())
  maxeval <- if (st$maxeval > 0) st$maxeval else NULL
  if (m && !grepl("AUGLAG|CCSA|MMA|COBYLA|ISRES|AGS|ORIG_DIRECT|SLSQP", name))
    stop("The algorithm ", name, " does not support constraints", call. = FALSE)
  res <- switch(kernel,
    "nelder-mead" = {
      # NLopt stops when any set tolerance holds, the SciPy simplex when both
      # of its tolerances hold: an unset one is made infinite, and with none
      # set the SciPy defaults apply
      xt <- max(st$xtol_abs, st$xtol_rel)
      ft <- max(st$ftol_abs, st$ftol_rel)
      o <- list(maxfev = maxeval %||% (200L * n), maxiter = .Machine$integer.max,
                xatol = if (xt > 0) xt else if (ft > 0) Inf else 1e-4,
                fatol = if (ft > 0) ft else if (xt > 0) Inf else 1e-4)
      bnd <- if (any(is.finite(c(lb, ub)))) list(lb = lb, ub = ub)
      r <- do.call(.sp_neldermead, c(list(f, x0, bounds = bnd), o))
      list(x = r$x, f = r$fun, ret = if (identical(r$status, 1L) || identical(r$status, 2L)) "MAXEVAL_REACHED" else
        "XTOL_REACHED", nevals = r$nfev)
    },
    "l-bfgs-b" = {
      r <- .sp_lbfgsb(f, x0, jac = grad, bounds = list(lb = lb, ub = ub),
                      ftol = if (st$ftol_rel > 0) st$ftol_rel else 2.220446049250313e-09,
                      gtol = if (st$ftol_rel > 0 || st$ftol_abs > 0) 0 else 1e-5,
                      maxfun = maxeval %||% 15000L, maxiter = 15000L)
      list(x = r$x, f = r$fun, ret = if (isTRUE(r$success)) "SUCCESS" else if (identical(r$status, 1L))
        "MAXEVAL_REACHED" else "FAILURE", nevals = r$nfev)
    },
    "cobyla" = {
      lc <- rep_len(as.numeric(a$lcons %||% numeric()), m)
      uc <- rep_len(as.numeric(a$ucons %||% numeric()), m)
      cons <- if (m) function(x) as.numeric(fn$cons(x, p))
      ie <- which(lc == uc); iu <- which(lc != uc & is.finite(uc)); il <- which(lc != uc & is.finite(lc))
      calcfc <- function(x) {
        cv <- if (m) cons(x) else numeric()
        list(f = f(x), c = c(cv[ie] - lc[ie], lc[ie] - cv[ie], cv[iu] - uc[iu], lc[il] - cv[il]))
      }
      step <- st$initial_step
      rhobeg <- if (!is.null(step)) max(step) else {
        w <- ub - lb
        s <- ifelse(is.finite(w), 0.25 * w, ifelse(x0 != 0, 0.25 * abs(x0), 1))
        max(min(s), 1e-4)
      }
      rhoend <- max(st$xtol_rel * max(abs(x0), 1), max(st$xtol_abs), 1e-8 * rhobeg)
      r0 <- calcfc(x0)
      r <- .pr_cobyla(calcfc, length(r0$c), x0, xl = lb, xu = ub, f0 = r0$f, nlconstr0 = r0$c, rhobeg = rhobeg,
                      rhoend = min(rhoend, rhobeg), maxfun = maxeval %||% (1000L * n), fortran = TRUE)
      list(x = r$x, f = r$f, ret = switch(as.character(r$info), "0" = "XTOL_REACHED", "1" = "STOPVAL_REACHED",
                                          "3" = "MAXEVAL_REACHED", "FAILURE"), nevals = r$nf)
    })
  if (is_opt) alg$numevals <- res$nevals
  retcode <- switch(res$ret, SUCCESS = , STOPVAL_REACHED = , FTOL_REACHED = , XTOL_REACHED = "Success",
                    MAXEVAL_REACHED = "MaxIters", MAXTIME_REACHED = "MaxTime", FORCED_STOP = "Terminated", "Failure")
  .opt_solution(a, alg, res$x, res$f, retcode, list(ret = res$ret, numevals = res$nevals))
}
