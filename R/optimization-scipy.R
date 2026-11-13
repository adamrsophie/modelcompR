# OptimizationSciPy (lib/OptimizationSciPy/src/OptimizationSciPy.jl): the
# SciPy optimizers driven through solve(). The algorithms themselves are
# native R ports of SciPy 1.18.1 (scipy/optimize), see the
# optimization-scipy-*.R files; this file holds the algorithm constructors,
# their traits and the translation of an OptimizationProblem into the
# arguments of the SciPy functions, as the Julia wrapper does.

.scipy_minimize_methods <- c("Nelder-Mead", "Powell", "CG", "BFGS", "Newton-CG", "L-BFGS-B", "TNC", "COBYLA",
                             "COBYQA", "SLSQP", "trust-constr", "dogleg", "trust-ncg", "trust-krylov",
                             "trust-exact")

.scipy_check_method <- function(method, valid) {
  if (!(method %in% valid))
    stop("Invalid method: ", method, ". Valid methods are: ", paste(valid, collapse = ", "), call. = FALSE)
  method
}

.scipy_alg <- function(name, method = NULL, ..., traits = list(), vector_objective = FALSE) {
  t <- utils::modifyList(list(callback = TRUE), traits)
  structure(list(name = name, method = method, ..., traits = t, vector_objective = vector_objective,
                 solve = .scipy_solve), class = c(name, "ScipyOptimizer", "OptimizationAlgorithm"))
}

#' @export
print.ScipyOptimizer <- function(x, ...) {
  cat(x$name, "(", if (!is.null(x$method)) paste0("\"", x$method, "\""), ")\n", sep = "")
  invisible(x)
}

#' SciPy optimizers
#'
#' Native R ports of the optimizers of SciPy 1.18.1 (`scipy.optimize`), with
#' the interface of OptimizationSciPy.jl: pass them to [solve()] with an
#' [OptimizationProblem()].
#'
#' * `ScipyMinimize(method)` runs `scipy.optimize.minimize` with one of the
#'   methods `"Nelder-Mead"`, `"Powell"`, `"CG"`, `"BFGS"` (the default),
#'   `"Newton-CG"`, `"L-BFGS-B"`, `"TNC"`, `"COBYLA"`, `"COBYQA"`, `"SLSQP"`,
#'   `"trust-constr"`, `"dogleg"`, `"trust-ncg"`, `"trust-krylov"` and
#'   `"trust-exact"`; `ScipyNelderMead()` ... `ScipyTrustExact()` are
#'   shorthands. Bounds are used by Nelder-Mead, Powell, L-BFGS-B, TNC,
#'   COBYLA, COBYQA, SLSQP and trust-constr, constraints by SLSQP,
#'   trust-constr, COBYLA and COBYQA. Without a gradient function, the
#'   gradient based methods use the forward differences of SciPy.
#' * `ScipyMinimizeScalar(method)` minimizes a function of one variable by
#'   `"brent"` (the default; `ScipyBrent()`), `"bounded"` (`ScipyBounded()`,
#'   on the interval of the bounds) or `"golden"` (`ScipyGolden()`).
#' * `ScipyLeastSquares(method, loss)` minimizes the sum of squares of the
#'   vector returned by the objective, by `"trf"` (`ScipyLeastSquaresTRF()`),
#'   `"dogbox"` (`ScipyLeastSquaresDogbox()`) or `"lm"`
#'   (`ScipyLeastSquaresLM()`, MINPACK's Levenberg Marquardt), with a robust
#'   `loss` (`"linear"`, `"soft_l1"`, `"huber"`, `"cauchy"`, `"arctan"`).
#'   The solution's objective is the cost `0.5 * sum(rho(r^2))`.
#' * `ScipyRootScalar(method)` finds a zero of a function of one variable:
#'   the bracketing methods `"brentq"`, `"brenth"`, `"bisect"`, `"ridder"`
#'   and `"toms748"` use the bounds as the bracket, `"newton"`, `"secant"`
#'   and `"halley"` start from `u0`.
#' * `ScipyRoot(method)` finds a zero of the vector function returned by the
#'   objective: `"hybr"` (MINPACK's Powell hybrid method), `"lm"`,
#'   `"broyden1"`, `"broyden2"`, `"anderson"`, `"linearmixing"`,
#'   `"diagbroyden"`, `"excitingmixing"`, `"krylov"` and `"df-sane"`.
#' * `ScipyLinprog(method)` solves the linear program whose cost vector is
#'   the value of the objective, with the options `A_ub`, `b_ub`, `A_eq`,
#'   `b_eq` given to `solve()` and the bounds of the problem.
#'   `ScipyMilp()` solves a mixed integer linear program with the options
#'   `integrality`, `A`, `lb_con` and `ub_con`.
#' * `ScipyDifferentialEvolution()`, `ScipyBasinhopping()`,
#'   `ScipyDualAnnealing()`, `ScipyShgo()`, `ScipyDirect()` and
#'   `ScipyBrute()` are the global optimizers of SciPy; all but
#'   basinhopping need finite bounds. The randomized ones take the option
#'   `seed` (an integer for [set.seed()]; the caller's random state is
#'   restored).
#'
#' SciPy solves linear and mixed integer programs with the HiGHS library;
#' here a native port of the HiGHS LP presolve and postsolve and of the
#' HiGHS simplex solver (dual simplex method with scaling, cost
#' perturbation, dual steepest edge pricing and the bound flipping ratio
#' test, primal simplex cleanup, the HFactor LU factorization with
#' Forrest-Tomlin updates) solves them, with
#' `presolve = FALSE` passed on as in SciPy; `"highs-ipm"` runs IPX with
#' crossover (`run_crossover = "off"` skips it) and the simplex clean up of
#' HiGHS, and `ScipyMilp()`
#' branches and bounds depth first over its LPs without the cuts and
#' heuristics of the HiGHS MIP solver, so node counts differ. `ScipyCOBYQA()` solves its
#' interpolation systems with R's LAPACK (SciPy uses Accelerate or
#' OpenBLAS), so its points agree to rounding and a long run can take a
#' slightly different path to the same solution; it evaluates the nonlinear
#' constraints at the point of the user also when variables are fixed by the
#' bounds or scaled (`scale = TRUE`), where SciPy 1.18.1 passes them the
#' reduced point.
#' `ScipyTrustConstr()` factorizes its augmented systems with a dense LU
#' (SciPy uses SuperLU), so a last bit difference can change its path while
#' the solution agrees. `ScipyShgo()` samples with the simplicial method only
#' (the other methods triangulate with Qhull).
#' `ScipyDirect()` drops the rectangles off the convex hull, as SciPy's
#' DIRECT code does only with its log switched on.
#'
#' Common options of `solve()`: `maxiters`, `abstol` and `reltol`, mapped to
#' the SciPy options as OptimizationSciPy.jl does; every other option is
#' passed to the SciPy function (for example `xatol`, `maxfev`, `adaptive`
#' for Nelder-Mead, `popsize` for differential evolution). The solution's
#' `original` is the SciPy result (`x`, `fun`, `nit`, `nfev`, `status`,
#' `message`, ...).
#'
#' @param method The SciPy method.
#' @param loss The robust loss of least squares.
#' @return An algorithm object.
#' @examples
#' rosen <- function(u, p) (p[1] - u[1])^2 + p[2] * (u[2] - u[1] * u[1])^2
#' prob <- OptimizationProblem(rosen, c(-1, 2), c(1, 100))
#' solve(prob, ScipyNelderMead())$u
#' solve(prob, ScipyBFGS())$u
#' @name scipy_optimizers
#' @export
ScipyMinimize <- function(method = "BFGS") {
  .scipy_check_method(method, .scipy_minimize_methods)
  gradient_free <- method %in% c("Nelder-Mead", "Powell", "COBYLA", "COBYQA")
  .scipy_alg("ScipyMinimize", method,
             traits = list(grad = !gradient_free,
                           hess = method %in% c("Newton-CG", "dogleg", "trust-ncg", "trust-exact", "trust-krylov"),
                           bounds = TRUE,
                           cons = method %in% c("SLSQP", "trust-constr", "COBYLA", "COBYQA"),
                           consjac = method %in% c("SLSQP", "trust-constr")))
}
#' @rdname scipy_optimizers
#' @export
ScipyNelderMead <- function() ScipyMinimize("Nelder-Mead")
#' @rdname scipy_optimizers
#' @export
ScipyPowell <- function() ScipyMinimize("Powell")
#' @rdname scipy_optimizers
#' @export
ScipyCG <- function() ScipyMinimize("CG")
#' @rdname scipy_optimizers
#' @export
ScipyBFGS <- function() ScipyMinimize("BFGS")
#' @rdname scipy_optimizers
#' @export
ScipyNewtonCG <- function() ScipyMinimize("Newton-CG")
#' @rdname scipy_optimizers
#' @export
ScipyLBFGSB <- function() ScipyMinimize("L-BFGS-B")
#' @rdname scipy_optimizers
#' @export
ScipyTNC <- function() ScipyMinimize("TNC")
#' @rdname scipy_optimizers
#' @export
ScipyCOBYLA <- function() ScipyMinimize("COBYLA")
#' @rdname scipy_optimizers
#' @export
ScipyCOBYQA <- function() ScipyMinimize("COBYQA")
#' @rdname scipy_optimizers
#' @export
ScipySLSQP <- function() ScipyMinimize("SLSQP")
#' @rdname scipy_optimizers
#' @export
ScipyTrustConstr <- function() ScipyMinimize("trust-constr")
#' @rdname scipy_optimizers
#' @export
ScipyDogleg <- function() ScipyMinimize("dogleg")
#' @rdname scipy_optimizers
#' @export
ScipyTrustNCG <- function() ScipyMinimize("trust-ncg")
#' @rdname scipy_optimizers
#' @export
ScipyTrustKrylov <- function() ScipyMinimize("trust-krylov")
#' @rdname scipy_optimizers
#' @export
ScipyTrustExact <- function() ScipyMinimize("trust-exact")

#' @rdname scipy_optimizers
#' @export
ScipyMinimizeScalar <- function(method = "brent") {
  .scipy_check_method(method, c("brent", "bounded", "golden"))
  .scipy_alg("ScipyMinimizeScalar", method, traits = list(bounds = method == "bounded"))
}
#' @rdname scipy_optimizers
#' @export
ScipyBrent <- function() ScipyMinimizeScalar("brent")
#' @rdname scipy_optimizers
#' @export
ScipyBounded <- function() ScipyMinimizeScalar("bounded")
#' @rdname scipy_optimizers
#' @export
ScipyGolden <- function() ScipyMinimizeScalar("golden")

#' @rdname scipy_optimizers
#' @export
ScipyLeastSquares <- function(method = "trf", loss = "linear") {
  .scipy_check_method(method, c("trf", "dogbox", "lm"))
  if (!(loss %in% c("linear", "soft_l1", "huber", "cauchy", "arctan")))
    stop("Invalid loss: ", loss, ". Valid loss functions are: linear, soft_l1, huber, cauchy, arctan", call. = FALSE)
  .scipy_alg("ScipyLeastSquares", method, loss = loss, traits = list(bounds = method %in% c("trf", "dogbox")),
             vector_objective = TRUE)
}
#' @rdname scipy_optimizers
#' @export
ScipyLeastSquaresTRF <- function() ScipyLeastSquares(method = "trf")
#' @rdname scipy_optimizers
#' @export
ScipyLeastSquaresDogbox <- function() ScipyLeastSquares(method = "dogbox")
#' @rdname scipy_optimizers
#' @export
ScipyLeastSquaresLM <- function() ScipyLeastSquares(method = "lm")

#' @rdname scipy_optimizers
#' @export
ScipyRootScalar <- function(method = "brentq") {
  .scipy_check_method(method, c("brentq", "brenth", "bisect", "ridder", "newton", "secant", "halley", "toms748"))
  .scipy_alg("ScipyRootScalar", method, traits = list(hess = method == "halley",
                                                      bounds = method %in% c("brentq", "brenth", "bisect", "ridder")))
}

#' @rdname scipy_optimizers
#' @export
ScipyRoot <- function(method = "hybr") {
  .scipy_check_method(method, c("hybr", "lm", "broyden1", "broyden2", "anderson", "linearmixing", "diagbroyden",
                                "excitingmixing", "krylov", "df-sane"))
  .scipy_alg("ScipyRoot", method, vector_objective = TRUE)
}

#' @rdname scipy_optimizers
#' @export
ScipyLinprog <- function(method = "highs") {
  .scipy_check_method(method, c("highs", "highs-ds", "highs-ipm", "interior-point", "revised simplex", "simplex"))
  .scipy_alg("ScipyLinprog", method, traits = list(bounds = TRUE, cons = TRUE))
}
#' @rdname scipy_optimizers
#' @export
ScipyMilp <- function() .scipy_alg("ScipyMilp", traits = list(bounds = TRUE, cons = TRUE))
#' @rdname scipy_optimizers
#' @export
ScipyDifferentialEvolution <- function() .scipy_alg("ScipyDifferentialEvolution", traits = list(bounds = TRUE))
#' @rdname scipy_optimizers
#' @export
ScipyBasinhopping <- function() .scipy_alg("ScipyBasinhopping", traits = list(bounds = TRUE))
#' @rdname scipy_optimizers
#' @export
ScipyDualAnnealing <- function() .scipy_alg("ScipyDualAnnealing", traits = list(bounds = TRUE))
#' @rdname scipy_optimizers
#' @export
ScipyShgo <- function() .scipy_alg("ScipyShgo", traits = list(bounds = TRUE, cons = TRUE))
#' @rdname scipy_optimizers
#' @export
ScipyDirect <- function() .scipy_alg("ScipyDirect", traits = list(bounds = TRUE))
#' @rdname scipy_optimizers
#' @export
ScipyBrute <- function() .scipy_alg("ScipyBrute", traits = list(bounds = TRUE))

# ---------------------------------------------------------------------------
# solve(): the __init and __solve methods of OptimizationSciPy.jl

# options consumed by the wrapper; the rest goes to the SciPy function
.scipy_exclude <- c("maxiters", "maxtime", "abstol", "reltol", "callback", "progress", "cons_tol", "jac_sparsity",
                    "keep_feasible", "hess_update", "verbose", "seed")

.scipy_kwargs <- function(opts) opts[setdiff(names(opts), .scipy_exclude)]

.scipy_retcode <- function(status, success) {
  if (isTRUE(success)) return("Success")
  switch(as.character(status %||% -1L), "0" = "Success", "1" = "MaxIters", "2" = "Infeasible", "3" = "Unstable",
         "4" = "Terminated", "9" = "MaxIters", "Failure")
}

# run `expr` with set.seed(seed), restoring the caller's random state
.scipy_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  expr
}

.scipy_solve <- function(a, alg, opts) {
  fn <- a$f
  p <- a$p
  x0 <- as.numeric(a$u0)
  n <- length(x0)
  lb <- a$lb; ub <- a$ub
  cls <- class(alg)[1L]
  # ScipyShgo too: upstream leaves its missing bounds to scipy.optimize.shgo,
  # which raises its own error
  if (cls %in% c("ScipyDifferentialEvolution", "ScipyDirect", "ScipyDualAnnealing", "ScipyBrute", "ScipyShgo") &&
      (is.null(lb) || is.null(ub)))
    stop(IncompatibleOptimizerError(paste0(cls, " requires bounds")))
  if (cls %in% c("ScipyMinimizeScalar", "ScipyRootScalar") && n != 1L)
    stop(cls, " requires exactly 1 variable, got ", n, ". Use ",
         if (cls == "ScipyMinimizeScalar") "ScipyMinimize" else "ScipyRoot", " for multivariate problems.",
         call. = FALSE)
  if (identical(alg$method, "bounded") && (is.null(lb) || is.null(ub)))
    stop("ScipyMinimizeScalar with method='bounded' requires bounds", call. = FALSE)
  if (cls == "ScipyRootScalar" && alg$method %in% c("brentq", "brenth", "bisect", "ridder") &&
      (is.null(lb) || is.null(ub)))
    stop("ScipyRootScalar with method='", alg$method, "' requires bracket (bounds)", call. = FALSE)
  if (!is.null(lb) && !is.null(ub) && any(lb > ub))
    stop("Lower bounds must be less than or equal to upper bounds", call. = FALSE)
  f <- function(x) fn$f(x, p)
  grad <- if (is.function(fn$grad)) function(x) as.numeric(fn$grad(x, p))
  hess <- if (is.function(fn$hess)) function(x) as.matrix(fn$hess(x, p))
  m <- length(a$lcons %||% numeric())
  cons <- if (is.function(fn$cons) && m) function(x) as.numeric(fn$cons(x, p))
  cons_j <- if (is.function(fn$cons_j) && m) function(x) {
    J <- fn$cons_j(x, p)
    if (!is.matrix(J)) J <- matrix(J, ncol = length(x))
    J
  }
  kw <- .scipy_kwargs(opts)
  runner <- c(ScipyMinimize = ".scipy_run_minimize", ScipyMinimizeScalar = ".scipy_run_scalar",
              ScipyLeastSquares = ".scipy_run_lsq", ScipyRootScalar = ".scipy_run_root_scalar",
              ScipyRoot = ".scipy_run_root", ScipyLinprog = ".scipy_run_linprog", ScipyMilp = ".scipy_run_milp",
              ScipyDifferentialEvolution = ".scipy_run_de", ScipyBasinhopping = ".scipy_run_basinhopping",
              ScipyDualAnnealing = ".scipy_run_dual_annealing", ScipyShgo = ".scipy_run_shgo",
              ScipyDirect = ".scipy_run_direct", ScipyBrute = ".scipy_run_brute")[[cls]]
  .scipy_need(runner, cls)
  res <- .scipy_with_seed(opts$seed, switch(cls,
    ScipyMinimize = .scipy_run_minimize(alg, f, grad, hess, cons, cons_j, a, x0, lb, ub, opts, kw),
    ScipyMinimizeScalar = .scipy_run_scalar(alg, f, lb, ub, kw),
    ScipyLeastSquares = .scipy_run_lsq(alg, f, x0, lb, ub, opts, kw),
    ScipyRootScalar = .scipy_run_root_scalar(alg, f, grad, hess, x0, lb, ub, kw),
    ScipyRoot = .scipy_run_root(alg, f, grad, x0, kw),
    ScipyLinprog = .scipy_run_linprog(alg, f, x0, lb, ub, opts, kw),
    ScipyMilp = .scipy_run_milp(alg, f, x0, lb, ub, a, kw),
    ScipyDifferentialEvolution = .scipy_run_de(f, lb, ub, cons, a, opts, kw),
    ScipyBasinhopping = .scipy_run_basinhopping(f, x0, opts, kw),
    ScipyDualAnnealing = .scipy_run_dual_annealing(f, lb, ub, opts, kw),
    ScipyShgo = .scipy_run_shgo(f, lb, ub, cons, a, opts, kw),
    ScipyDirect = .scipy_run_direct(f, lb, ub, opts, kw),
    ScipyBrute = .scipy_run_brute(f, lb, ub, kw)))
  retcode <- res$retcode %||% .scipy_retcode(res$status, res$success)
  .opt_solution(a, alg, res$x, res$objective %||% res$fun, retcode, res)
}

# ScipyMinimize: the option mapping of the Julia wrapper, then minimize()
.scipy_run_minimize <- function(alg, f, grad, hess, cons, cons_j, a, x0, lb, ub, opts, kw) {
  method <- alg$method
  options <- list()
  if (method == "trust-constr") {
    options <- list(initial_tr_radius = 1, verbose = 0, finite_diff_rel_step = 1e-8, gtol = 1e-10, maxiter = 50000L)
  } else if (method %in% c("dogleg", "trust-ncg", "trust-krylov", "trust-exact")) {
    options <- list(gtol = 1e-10, maxiter = 50000L)
  }
  if (!is.null(opts$maxiters)) options$maxiter <- as.integer(opts$maxiters)
  if (!is.null(opts$abstol)) {
    if (method %in% c("Nelder-Mead", "Powell")) options$xatol <- opts$abstol
    else if (method %in% c("L-BFGS-B", "TNC", "SLSQP", "trust-constr")) options$ftol <- opts$abstol
    else if (method == "COBYQA") options$feasibility_tol <- opts$abstol
  }
  if (!is.null(opts$reltol) && method %in% c("CG", "BFGS", "Newton-CG", "L-BFGS-B", "TNC", "SLSQP", "trust-constr"))
    options$gtol <- opts$reltol
  options <- utils::modifyList(options, kw)
  # the gradient is used by the gradient based methods only
  jac <- if (alg$traits$grad) grad
  # the Hessian: without one, the trust region methods use the BFGS update
  # (OptimizationSciPy.jl asks for it through options["hess"], which SciPy
  # rejects; see the porting notes)
  hs <- NULL
  if (alg$traits$hess || method == "trust-constr") {
    hs <- hess
    if (is.null(hs)) {
      if (method %in% c("trust-constr", "dogleg", "trust-ncg", "trust-krylov", "trust-exact")) {
        hs <- .sp_hessian_update("BFGS")
      } else {
        stop("Method ", method, " requires Hessian but none was provided", call. = FALSE)
      }
    }
    if (method == "trust-constr" && !is.null(opts$hess_update)) hs <- .sp_hessian_update(opts$hess_update)
  }
  bounds <- if (!is.null(lb) && !is.null(ub) &&
                method %in% c("L-BFGS-B", "TNC", "SLSQP", "trust-constr", "COBYLA", "COBYQA")) list(lb = lb, ub = ub)
  constraints <- list()
  if (alg$traits$cons) {
    if (!is.null(cons)) {
      cj <- if (alg$traits$consjac && !is.null(cons_j)) cons_j else "2-point"
      constraints <- list(.sp_nonlinear_constraint(cons, a$lcons, a$ucons, jac = cj,
                                                   keep_feasible = isTRUE(opts$keep_feasible),
                                                   finite_diff_rel_step = opts$cons_tol %||% 1e-8))
    }
  } else if (!is.null(cons)) {
    stop("Method ", method, " does not support constraints. Use SLSQP, trust-constr, COBYLA, or COBYQA instead.",
         call. = FALSE)
  }
  .sp_minimize(f, x0, method = method, jac = jac, hess = hs, bounds = bounds, constraints = constraints,
               options = options)
}

# ScipyMinimizeScalar: minimize_scalar() without a bracket, as the wrapper
.scipy_run_scalar <- function(alg, f, lb, ub, kw) {
  fs <- function(x) as.numeric(f(x))
  bounds <- if (alg$method == "bounded") c(lb[1L], ub[1L])
  # the solver keywords are minimize_scalar() arguments: bracket, tol and options
  bad <- setdiff(names(kw), c("bracket", "tol", "options"))
  if (length(bad)) stop("minimize_scalar() got an unexpected keyword argument '", bad[1L], "'", call. = FALSE)
  .sp_minimize_scalar(fs, bracket = kw$bracket, bounds = bounds, method = alg$method, tol = kw$tol,
                      options = kw$options %||% list())
}
