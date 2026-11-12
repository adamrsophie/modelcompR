# Optimizers for OptimizationProblem, named as in Optimization.jl and
# OptimizationOptimJL (LBFGS, BFGS, NelderMead; Optimization.LBFGS is the
# L-BFGS-B method) and Optimization.AugLag. The unconstrained and box
# constrained methods wrap stats::optim and stats::nlminb; general
# constraints are handled by an augmented Lagrangian method written here
# (the Powell, Hestenes and Rockafellar form for inequalities, as
# Optimization.jl's AugLag), with L-BFGS-B as the inner solver.
#
# Multiobjective problems (MultiObjectiveOptimizationFunction) are solved by
# the NSGA-II method that Optimization.jl reaches through
# OptimizationEvolutionary and OptimizationMetaheuristics
# (lib/OptimizationEvolutionary/src/OptimizationEvolutionary.jl: the
# solution's u is the Pareto set) and by a weighted sum scan; problems with
# integer and binary unknowns by branch and bound over the continuous
# relaxation, the method of the Juniper solver of OptimizationMOI.

.opt_alg <- function(name, ...) structure(list(name = name, ...), class = c(name, "OptimizationAlgorithm"))

#' @export
print.OptimizationAlgorithm <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

#' Optimization algorithms
#'
#' Algorithms for [solve()] on an [OptimizationProblem()].
#'
#' * `LBFGSB()` is the limited memory BFGS method with box constraints
#'   ([stats::optim()] method `"L-BFGS-B"`); `LBFGS()` is another name for
#'   it. It is the default for problems without constraints.
#' * `BFGS()` and `NelderMead()` are the `stats::optim()` methods of those
#'   names; they do not support bounds.
#' * `Nlminb()` is the PORT quasi Newton method of [stats::nlminb()], with
#'   box constraints; it uses the symbolic Hessian when the problem has one.
#' * `AugLag()` is an augmented Lagrangian method for equality and
#'   inequality constraints, with bounds handled by the inner L-BFGS-B
#'   solver. It is the default for problems with constraints.
#'
#' Gradients are the symbolic ones of the problem when present, and finite
#' differences otherwise.
#'
#' Solver options given to `solve()`: `maxiters`, `abstol`, `reltol`.
#'
#' @param m Number of corrections kept by L-BFGS-B.
#' @param rho0 Initial penalty parameter.
#' @param gamma Factor by which the penalty grows when the constraint
#'   violation does not decrease by `tau`.
#' @param tau Required decrease of the constraint violation.
#' @param ctol Tolerance on the constraint violation.
#' @param max_outer Maximum number of multiplier updates.
#' @return An algorithm object.
#' @examples
#' rosenbrock <- function(u, p) (p[1] - u[1])^2 + p[2] * (u[2] - u[1]^2)^2
#' prob <- OptimizationProblem(rosenbrock, c(-1, 2), c(1, 100))
#' solve(prob, LBFGSB())$u
#' solve(prob, BFGS())$u
#' solve(prob, NelderMead())$u
#' solve(prob, Nlminb())$u
#' @name optimization_algorithms
#' @export
LBFGSB <- function(m = 10L) .opt_alg("LBFGSB", m = as.integer(m))
#' @rdname optimization_algorithms
#' @export
LBFGS <- function(m = 10L) .opt_alg("LBFGS", m = as.integer(m))
#' @rdname optimization_algorithms
#' @export
BFGS <- function() .opt_alg("BFGS")
#' @rdname optimization_algorithms
#' @export
NelderMead <- function() .opt_alg("NelderMead")
#' @rdname optimization_algorithms
#' @export
Nlminb <- function() .opt_alg("Nlminb")
#' @rdname optimization_algorithms
#' @export
AugLag <- function(rho0 = 10, gamma = 10, tau = 0.25, ctol = 1e-8, max_outer = 50L, m = 10L) {
  .opt_alg("AugLag", rho0 = rho0, gamma = gamma, tau = tau, ctol = ctol, max_outer = as.integer(max_outer),
           m = as.integer(m))
}

.opt_fd_gradient <- function(f, x, fx = f(x)) {
  g <- numeric(length(x))
  for (j in seq_along(x)) {
    h <- .Machine$double.eps^(1 / 3) * max(1, abs(x[j]))
    xp <- x; xm <- x
    xp[j] <- x[j] + h; xm[j] <- x[j] - h
    g[j] <- (f(xp) - f(xm)) / (2 * h)
  }
  g
}

# L-BFGS-B through stats::optim with tight tolerances.
.opt_lbfgsb <- function(f, g, x0, lb, ub, maxiters, m = 10L, factr = 10, pgtol = 0) {
  n <- length(x0)
  lo <- lb %||% rep(-Inf, n); hi <- ub %||% rep(Inf, n)
  x0 <- pmin(pmax(x0, lo), hi)
  big <- 1e100
  safe <- function(x) {
    v <- f(x)
    if (is.finite(v)) v else big
  }
  # a trial point where the gradient overflows is rejected by the line
  # search through the huge value of `safe`; its gradient only has to be finite
  if (!is.null(g)) {
    g0 <- g
    g <- function(x) {
      gx <- g0(x)
      gx[!is.finite(gx)] <- 0
      gx
    }
  }
  res <- stats::optim(x0, safe, g, method = "L-BFGS-B", lower = lo, upper = hi,
                      control = list(maxit = maxiters, factr = factr, pgtol = pgtol, lmm = m))
  if (res$convergence %in% c(51L, 52L)) {
    # the line search gave up at the tight tolerance; confirm convergence
    # from that point at the default tolerance of optim()
    res2 <- stats::optim(res$par, safe, g, method = "L-BFGS-B", lower = lo, upper = hi,
                         control = list(maxit = maxiters, lmm = m))
    if (res2$convergence == 0L && res2$value <= res$value) {
      res2$counts <- res$counts + res2$counts
      res <- res2
    }
  }
  res
}

.opt_retcode <- function(code) {
  if (code == 0L) "Success" else if (code == 1L) "MaxIters" else "Failure"
}

# Augmented Lagrangian for min f(x) subject to ceq(x) = 0, cin(x) <= 0 and
# bounds. `cons(x)` returns the constraint function g and (lcons, ucons)
# give its bounds; rows with lcons == ucons are equalities.
.opt_auglag <- function(f, g, cons, cons_j, lcons, ucons, x0, lb, ub, alg, maxiters) {
  eq <- which(lcons == ucons)
  up <- which(lcons != ucons & is.finite(ucons))
  lo <- which(lcons != ucons & is.finite(lcons))
  # c(x): equalities first, then inequalities written as c <= 0
  cval <- function(x) {
    gx <- cons(x)
    list(e = gx[eq] - ucons[eq], i = c(gx[up] - ucons[up], lcons[lo] - gx[lo]))
  }
  cjac <- function(x) {
    J <- if (is.null(cons_j)) .fd_jacobian(cons, x) else cons_j(x)
    if (!is.matrix(J)) J <- matrix(J, ncol = length(x))
    list(e = J[eq, , drop = FALSE], i = rbind(J[up, , drop = FALSE], -J[lo, , drop = FALSE]))
  }
  x <- as.numeric(x0)
  lam <- numeric(length(eq)); mu <- numeric(length(up) + length(lo))
  rho <- alg$rho0
  viol <- function(cv) max(0, abs(cv$e), cv$i)
  cv <- cval(x)
  vprev <- viol(cv)
  retcode <- "MaxIters"
  inner_evals <- 0L
  for (outer in seq_len(alg$max_outer)) {
    L <- function(z) {
      cz <- cval(z)
      f(z) + sum(lam * cz$e) + rho / 2 * sum(cz$e^2) + sum(pmax(0, mu + rho * cz$i)^2 - mu^2) / (2 * rho)
    }
    gL <- function(z) {
      cz <- cval(z); Jz <- cjac(z)
      gf <- if (is.null(g)) .opt_fd_gradient(f, z) else g(z)
      gf + as.numeric(crossprod(Jz$e, lam + rho * cz$e)) + as.numeric(crossprod(Jz$i, pmax(0, mu + rho * cz$i)))
    }
    res <- .opt_lbfgsb(L, gL, x, lb, ub, maxiters, alg$m)
    inner_evals <- inner_evals + sum(res$counts, na.rm = TRUE)
    step <- max(abs(res$par - x), 0)
    x <- res$par
    cv <- cval(x)
    v <- viol(cv)
    lam <- lam + rho * cv$e
    mu <- pmax(0, mu + rho * cv$i)
    if (v <= alg$ctol && (step <= 1e-9 * max(1, abs(x)) || outer > 1L && abs(v - vprev) <= alg$ctol)) {
      # stationarity of the Lagrangian with the updated multipliers
      Jz <- cjac(x)
      gf <- if (is.null(g)) .opt_fd_gradient(f, x) else g(x)
      kkt <- gf + as.numeric(crossprod(Jz$e, lam)) + as.numeric(crossprod(Jz$i, mu))
      if (!is.null(lb)) kkt[x <= lb & kkt > 0] <- 0
      if (!is.null(ub)) kkt[x >= ub & kkt < 0] <- 0
      # complementary slackness: a multiplier is positive only on an active
      # inequality
      slack <- max(abs(mu * cv$i), 0) <= 1e-6 * max(1, abs(mu))
      if (slack && max(abs(kkt), 0) <= 1e-5 * max(1, abs(gf))) {
        retcode <- "Success"
        break
      }
    }
    if (v > alg$tau * vprev) rho <- min(rho * alg$gamma, 1e12)
    vprev <- v
  }
  if (retcode != "Success" && viol(cv) > 1e-6) retcode <- "Infeasible"
  list(u = x, objective = f(x), retcode = retcode, lambda = lam, mu = mu, rho = rho, viol = viol(cv),
       outer = outer, evals = inner_evals)
}

#' Solve optimization problems
#'
#' `solve(prob, alg)` minimizes an [OptimizationProblem()]. The solution
#' has `sol$u` (the minimizer, a named vector), `sol$objective` (the cost
#' there), `sol$retcode` and `sol$original` (the result of the underlying
#' optimizer), and supports symbolic indexing: `sol[x]`, `sol$ps[a]`.
#'
#' Problems with integer or binary unknowns are solved by
#' [BranchAndBound()] (the default for them), and multiobjective problems
#' (built with `multiobjective = TRUE`) by [NSGA2()] or [WeightedSumScan()]
#' (the default); their solution holds the Pareto set in `sol$u` and the
#' Pareto front in `sol$objective`, see [multiobjective_algorithms].
#'
#' @param a An `OptimizationProblem`.
#' @param b An algorithm, see [optimization_algorithms]; by default
#'   `AugLag()` for problems with constraints and `LBFGSB()` otherwise,
#'   `BranchAndBound()` for problems with integer unknowns and
#'   `WeightedSumScan()` for multiobjective problems.
#' @param ... Solver options: `maxiters`, `abstol`, `reltol`.
#' @return An `OptimizationSolution`.
#' @examples
#' x <- Variable("x"); y <- Variable("y")
#' sys <- complete(System(list(), costs = list((1 - x)^2 + 100 * (y - x^2)^2), name = "rosenbrock",
#'                        constraints = list(Inequality(x^2 + y^2, 1))))
#' sol <- solve(OptimizationProblem(sys, list(x ~ 0, y ~ 0)), AugLag())
#' sol$u; sol$objective; sol$retcode
#' @export
solve.OptimizationProblem <- function(a, b, ...) {
  opts <- utils::modifyList(a$kwargs %||% list(), list(...))
  b <- if (missing(b)) NULL else b
  if (!is.null(b) && !inherits(b, "OptimizationAlgorithm"))
    stop(IncompatibleOptimizerError(paste0("not an optimization algorithm: ", class(b)[1L])))
  if (inherits(a$f, "MultiObjectiveOptimizationFunction")) return(.opt_solve_core(a, b, opts))
  # the objective sense, the callback and the statistics (OptimizationBase's
  # OptimizationCache and __solve): wrap the functions of the problem
  start <- proc.time()[["elapsed"]]
  cnt <- new.env(parent = emptyenv())
  cnt$f <- 0L; cnt$g <- 0L; cnt$h <- 0L; cnt$best <- NULL
  sgn <- if (identical(unclass(a$sense %||% MinSense), "MaxSense")) -1 else 1
  cb <- opts$callback
  opts$callback <- NULL
  fn <- a$f
  w <- a
  w$f$f <- function(x, p) {
    v <- sgn * as.numeric(fn$f(x, p))
    cnt$f <- cnt$f + 1L
    if (is.null(cnt$best) || v < cnt$best$v) cnt$best <- list(x = x, v = v)
    if (is.function(cb) && isTRUE(cb(list(u = x, objective = sgn * v, iter = cnt$f, p = p), sgn * v)))
      stop(structure(class = c("opt_callback_stop", "condition"), list(message = "stopped", call = NULL)))
    v
  }
  if (is.function(fn$grad)) w$f$grad <- function(x, p) {
    cnt$g <- cnt$g + 1L
    sgn * fn$grad(x, p)
  }
  if (is.function(fn$hess)) w$f$hess <- function(x, p) {
    cnt$h <- cnt$h + 1L
    sgn * fn$hess(x, p)
  }
  sol <- tryCatch(.opt_solve_core(w, b, opts), opt_callback_stop = function(e) NULL)
  if (is.null(sol)) {
    # the callback stopped the optimization: the best point so far
    alg <- b %||% (if (is.function(fn$cons)) AugLag() else LBFGSB())
    sol <- .opt_solution(w, alg, cnt$best$x, cnt$best$v, "Terminated", NULL)
  }
  obj <- sgn * (sol$objective %||% 0)
  stats <- OptimizationStats(iterations = cnt$f, time = proc.time()[["elapsed"]] - start, fevals = cnt$f,
                             gevals = cnt$g, hevals = cnt$h)
  .opt_set_fields(sol, objective = obj, stats = stats, prob = a)
}

.opt_set_fields <- function(sol, ...) {
  env <- .sol_env(sol)
  vals <- list(...)
  for (nm in names(vals)) env[[nm]] <- vals[[nm]]
  sol
}

# solve!(cache) of OptimizationBase
.opt_solve_cache <- function(cache) {
  if (is.null(cache$sol)) cache$sol <- do.call(solve, c(list(cache$prob, cache$alg), cache$opts))
  cache$sol
}

.opt_solve_core <- function(a, b, opts) {
  if (inherits(a$f, "MultiObjectiveOptimizationFunction")) return(.opt_solve_multiobjective(a, b, opts))
  if (any(a$int %||% FALSE) || inherits(b, "BranchAndBound")) {
    if (!is.null(b) && !inherits(b, "BranchAndBound"))
      stop("the problem has integer or binary unknowns, which ", b$name, "() does not handle; use ",
           "BranchAndBound(relaxation = ", b$name, "())", call. = FALSE)
    return(.opt_branch_and_bound(a, b %||% BranchAndBound(), opts))
  }
  p <- a$p
  fn <- a$f
  f <- function(x) as.numeric(fn$f(x, p))
  g <- if (is.function(fn$grad)) function(x) as.numeric(fn$grad(x, p)) else NULL
  has_cons <- is.function(fn$cons)
  alg <- b %||% (if (has_cons) AugLag() else LBFGSB())
  if (!inherits(alg, "OptimizationAlgorithm")) stop("not an optimization algorithm: ", class(alg)[1L], call. = FALSE)
  if (has_cons && !inherits(alg, "AugLag"))
    stop(alg$name, "() does not handle constraints; use AugLag(), or constraints_to_penalties()", call. = FALSE)
  lb <- a$lb; ub <- a$ub
  bounded <- !is.null(lb) && any(is.finite(lb)) || !is.null(ub) && any(is.finite(ub))
  if (bounded && inherits(alg, c("BFGS", "NelderMead")))
    stop(alg$name, "() does not handle bounds; use LBFGSB(), Nlminb() or AugLag()", call. = FALSE)
  maxiters <- as.integer(opts$maxiters %||% 1000L)
  x0 <- as.numeric(a$u0)
  if (!length(x0)) stop("the problem has no optimization variables", call. = FALSE)
  if (inherits(alg, "AugLag")) {
    if (has_cons) {
      cons <- function(x) as.numeric(fn$cons(x, p))
      cj <- if (is.function(fn$cons_j)) function(x) fn$cons_j(x, p) else NULL
      m <- length(cons(x0))
      res <- .opt_auglag(f, g, cons, cj, a$lcons %||% rep(0, m), a$ucons %||% rep(0, m), x0, lb, ub, alg, maxiters)
    } else {
      r <- .opt_lbfgsb(f, g, x0, lb, ub, maxiters, alg$m)
      res <- list(u = r$par, objective = r$value, retcode = .opt_retcode(r$convergence))
    }
    return(.opt_solution(a, alg, res$u, res$objective, res$retcode, res))
  }
  if (inherits(alg, c("LBFGSB", "LBFGS"))) {
    factr <- if (is.null(opts$reltol)) 10 else opts$reltol / .Machine$double.eps
    r <- .opt_lbfgsb(f, g, x0, lb, ub, maxiters, alg$m, factr = factr, pgtol = opts$abstol %||% 0)
    return(.opt_solution(a, alg, r$par, r$value, .opt_retcode(r$convergence), r))
  }
  if (inherits(alg, "Nlminb")) {
    h <- if (is.function(fn$hess)) function(x) fn$hess(x, p) else NULL
    n <- length(x0)
    r <- stats::nlminb(x0, f, gradient = g, hessian = h, lower = lb %||% rep(-Inf, n), upper = ub %||% rep(Inf, n),
                       control = list(iter.max = maxiters, eval.max = 2L * maxiters,
                                      rel.tol = opts$reltol %||% 1e-10, abs.tol = opts$abstol %||% 0))
    rc <- if (r$convergence == 0L) "Success" else if (r$iterations >= maxiters) "MaxIters" else "Failure"
    return(.opt_solution(a, alg, r$par, r$objective, rc, r))
  }
  method <- if (inherits(alg, "BFGS")) "BFGS" else if (inherits(alg, "NelderMead")) "Nelder-Mead" else
    stop("not an optimization algorithm: ", class(alg)[1L], call. = FALSE)
  if (method == "Nelder-Mead" && length(x0) == 1L) {
    # stats::optim() warns that Nelder-Mead is unreliable in one dimension
    r <- stats::optim(x0, f, g, method = "BFGS", control = list(maxit = maxiters, reltol = opts$reltol %||% 1e-14))
  } else {
    ctrl <- list(maxit = if (method == "Nelder-Mead") max(maxiters, 5000L) else maxiters,
                 reltol = opts$reltol %||% 1e-14)
    if (!is.null(opts$abstol)) ctrl$abstol <- opts$abstol
    r <- stats::optim(x0, f, if (method == "BFGS") g, method = method, control = ctrl)
  }
  .opt_solution(a, alg, r$par, r$value, .opt_retcode(r$convergence), r)
}

.opt_solution <- function(prob, alg, u, objective, retcode, original) {
  .nl_new_solution(c("OptimizationSolution", "NonlinearSolution"), prob, alg, u, retcode,
                   extra = list(objective = objective, original = original))
}

# ---------------------------------------------------------------------------
# Integer and binary unknowns

#' Branch and bound for integer unknowns
#'
#' `BranchAndBound()` solves an [OptimizationProblem()] whose unknowns
#' include integers (`Variable("n", integer = TRUE)`) or binaries
#' (`binary = TRUE`, bounds 0 and 1) by branch and bound over the
#' continuous relaxation, the method of the Juniper solver that
#' Optimization.jl reaches through OptimizationMOI. Each node solves the
#' relaxation with the bounds of its branch; a node whose relaxed
#' objective is no better than the best integer point found is pruned, and
#' otherwise the most fractional integer unknown is branched on. It is the
#' default for problems with integer unknowns.
#'
#' The relaxations are solved by local methods, so the result is the global
#' optimum when the relaxed problem is convex (for example a convex cost
#' with linear constraints); otherwise it is a good integer point but not
#' guaranteed optimal.
#'
#' @param relaxation The algorithm for the relaxations; `AugLag()` for
#'   problems with constraints and `LBFGSB()` otherwise when `NULL`.
#' @param max_nodes Largest number of relaxations to solve.
#' @param int_tol Distance to the nearest integer below which a value
#'   counts as integral.
#' @param ctol Largest constraint violation of an accepted point.
#' @return An algorithm object.
#' @examples
#' x <- Variable("x", integer = TRUE, bounds = c(-10, 10)); y <- Variable("y")
#' sys <- complete(System(list(), costs = list((x - 2.6)^2 + (y - 1.3)^2), name = "intq"))
#' sol <- solve(OptimizationProblem(sys, list(x ~ 0, y ~ 0)), BranchAndBound())
#' sol$u
#' @export
BranchAndBound <- function(relaxation = NULL, max_nodes = 1000L, int_tol = 1e-6, ctol = 1e-6) {
  .opt_alg("BranchAndBound", relaxation = relaxation, max_nodes = as.integer(max_nodes), int_tol = int_tol,
           ctol = ctol)
}

.opt_violation <- function(fn, x, p, lcons, ucons) {
  if (!is.function(fn$cons)) return(0)
  cx <- as.numeric(fn$cons(x, p))
  max(0, lcons - cx, cx - ucons)
}

.opt_branch_and_bound <- function(a, alg, opts) {
  int <- a$int %||% rep(FALSE, length(a$u0))
  n <- length(a$u0)
  fn <- a$f
  p <- a$p
  lb0 <- a$lb %||% rep(-Inf, n); ub0 <- a$ub %||% rep(Inf, n)
  lb0[int] <- ceiling(lb0[int] - alg$int_tol); ub0[int] <- floor(ub0[int] + alg$int_tol)
  m <- if (is.function(fn$cons)) length(fn$cons(as.numeric(a$u0), p)) else 0L
  lcons <- a$lcons %||% rep(0, m); ucons <- a$ucons %||% rep(0, m)
  relax <- alg$relaxation
  f <- function(x) as.numeric(fn$f(x, p))
  # a node is a box; nodes are explored depth first
  stack <- list(list(lb = lb0, ub = ub0, x = pmin(pmax(as.numeric(a$u0), lb0), ub0)))
  best_x <- NULL; best_f <- Inf
  nodes <- 0L
  unsure <- FALSE
  while (length(stack) && nodes < alg$max_nodes) {
    node <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL
    if (any(node$lb > node$ub)) next
    nodes <- nodes + 1L
    sub <- a
    sub$int <- rep(FALSE, n)
    sub$lb <- node$lb; sub$ub <- node$ub
    sub$u0 <- stats::setNames(node$x, names(a$u0))
    sol <- tryCatch(do.call(solve, c(list(sub, relax), opts)), error = function(e) NULL)
    if (is.null(sol)) {
      unsure <- TRUE
      next
    }
    x <- as.numeric(sol$u)
    fx <- f(x)
    if (!is.finite(fx) || .opt_violation(fn, x, p, lcons, ucons) > alg$ctol) {
      # an infeasible relaxation closes the node; an unfinished one leaves
      # the result uncertain
      if (!identical(sol$retcode, "Infeasible")) unsure <- TRUE
      next
    }
    solved <- identical(sol$retcode, "Success")
    if (!solved) unsure <- TRUE
    # ponytail: the bound is the local optimum of the relaxation, exact only
    # for convex relaxations; a global NLP solver would make it certified
    if (solved && is.finite(best_f) && fx >= best_f - 1e-9 * (1 + abs(best_f))) next
    frac <- abs(x - round(x))
    frac[!int] <- 0
    if (all(frac <= alg$int_tol)) {
      xr <- x
      xr[int] <- round(x[int])
      fr <- f(xr)
      if (is.finite(fr) && fr < best_f && .opt_violation(fn, xr, p, lcons, ucons) <= alg$ctol) {
        best_x <- xr; best_f <- fr
      }
      next
    }
    j <- which.max(frac)
    down <- node; down$ub[j] <- floor(x[j]); down$x <- x; down$x[j] <- floor(x[j])
    up <- node; up$lb[j] <- ceiling(x[j]); up$x <- x; up$x[j] <- ceiling(x[j])
    # the branch toward the nearer integer is explored first
    stack <- c(stack, if (x[j] - floor(x[j]) < 0.5) list(up, down) else list(down, up))
  }
  if (is.null(best_x)) {
    retcode <- if (length(stack)) "MaxIters" else "Infeasible"
    return(.opt_solution(a, alg, as.numeric(a$u0), f(as.numeric(a$u0)), retcode, list(nodes = nodes)))
  }
  retcode <- if (length(stack)) "MaxIters" else if (unsure) "ConvergenceFailure" else "Success"
  .opt_solution(a, alg, best_x, best_f, retcode, list(nodes = nodes, open_nodes = length(stack)))
}

# ---------------------------------------------------------------------------
# Multiobjective problems

#' Multiobjective optimization algorithms
#'
#' Algorithms for an [OptimizationProblem()] built with `multiobjective =
#' TRUE` (or from a [MultiObjectiveOptimizationFunction()]), whose objective
#' is the vector of the costs of the system. Both return an approximation of
#' the Pareto front: the solution's `sol$u` is the Pareto set, a matrix with
#' one row per point and one column per unknown, `sol$objective` the matrix
#' of the cost vectors at those points, and `sol[x]` the values of `x` along
#' the front.
#'
#' * `NSGA2()` is the elitist nondominated sorting genetic algorithm of Deb
#'   and coauthors, as `Evolutionary.NSGA2` and `Metaheuristics.NSGA2` used
#'   through Optimization.jl: simulated binary crossover, polynomial
#'   mutation, selection by nondomination rank and crowding distance, and
#'   constraints handled by Deb's feasibility rules. It needs finite bounds on
#'   every unknown and finds nonconvex fronts.
#' * `WeightedSumScan()` minimizes `sum(w * (f - ideal) / (nadir - ideal))`
#'   for weights `w` on a regular grid of the simplex with a gradient based
#'   method, starting each solve from the previous point. The ideal and
#'   nadir points come from minimizing each cost alone. A weighted sum
#'   reaches only the convex parts of a front. It is the default.
#'
#' @param population Population size.
#' @param generations Number of generations.
#' @param crossover_rate Probability of crossover of a pair of parents.
#' @param mutation_rate Probability of mutating each unknown; `1 / n` when
#'   `NULL`.
#' @param eta_c,eta_m Distribution indices of the crossover and mutation.
#' @param seed Seed of the random numbers, or `NULL` to use the current
#'   stream. With a seed the caller's random number state is restored.
#' @param n Number of weights along each edge of the simplex.
#' @param inner The algorithm for each scalarized problem; `AugLag()` with
#'   constraints and `LBFGSB()` otherwise when `NULL`.
#' @return An algorithm object.
#' @examples
#' x <- Variable("x", bounds = c(-1, 3))
#' sys <- complete(System(list(), costs = list(x^2, (x - 2)^2), name = "mo"))
#' prob <- OptimizationProblem(sys, list(x ~ 1), multiobjective = TRUE)
#' sol <- solve(prob, WeightedSumScan(n = 5))
#' sol$u
#' sol$objective
#' @name multiobjective_algorithms
#' @export
NSGA2 <- function(population = 100L, generations = 200L, crossover_rate = 0.9, mutation_rate = NULL,
                  eta_c = 20, eta_m = 20, seed = NULL) {
  .opt_alg("NSGA2", population = 2L * ((as.integer(population) + 1L) %/% 2L),
           generations = as.integer(generations), crossover_rate = crossover_rate, mutation_rate = mutation_rate,
           eta_c = eta_c, eta_m = eta_m, seed = seed)
}
#' @rdname multiobjective_algorithms
#' @export
WeightedSumScan <- function(n = 21L, inner = NULL) .opt_alg("WeightedSumScan", points = as.integer(n), inner = inner)

.opt_solve_multiobjective <- function(a, alg, opts) {
  alg <- alg %||% WeightedSumScan()
  if (!inherits(alg, c("NSGA2", "WeightedSumScan")))
    stop("this is a multiobjective problem: its objective is a vector. Use NSGA2() or WeightedSumScan(), or ",
         "build the problem with `weights` to minimize a weighted sum of the costs.", call. = FALSE)
  if (any(a$int %||% FALSE)) stop("integer unknowns are not supported in multiobjective problems", call. = FALSE)
  res <- if (inherits(alg, "NSGA2")) .opt_nsga2(a, alg, opts) else .opt_weighted_scan(a, alg, opts)
  X <- res$X
  colnames(X) <- names(a$u0)
  Fm <- res$F
  colnames(Fm) <- paste0("cost", seq_len(ncol(Fm)))
  .new_solution(c("MultiObjectiveOptimizationSolution", "OptimizationSolution", "NonlinearSolution"),
                list(t = seq_len(nrow(X)), U = X, prob = a, alg = alg, retcode = res$retcode, objective = Fm,
                     stats = res$stats %||% list(), p_final = a$p, sys = a$sys, original = res))
}

#' @export
`$.MultiObjectiveOptimizationSolution` <- function(x, name) {
  if (name == "u") return(.sol_env(x)$U)
  `$.SciMLSolution`(x, name)
}

#' @export
print.MultiObjectiveOptimizationSolution <- function(x, ...) {
  env <- .sol_env(x)
  cat("MultiObjectiveOptimizationSolution with retcode ", env$retcode, " (", env$alg$name, ")\n", sep = "")
  cat("Pareto front of ", nrow(env$U), " points; the first:\n", sep = "")
  print(utils::head(cbind(env$U, env$objective)))
  invisible(x)
}

# Pareto dominance under Deb's constraint rules: d[i, j] is TRUE when point
# i dominates point j.
.opt_dominance <- function(Fm, cv) {
  np <- nrow(Fm)
  le <- matrix(TRUE, np, np); lt <- matrix(FALSE, np, np)
  for (k in seq_len(ncol(Fm))) {
    le <- le & outer(Fm[, k], Fm[, k], "<=")
    lt <- lt | outer(Fm[, k], Fm[, k], "<")
  }
  d <- le & lt
  feas <- cv <= 0
  both <- outer(feas, feas, "&")
  d[!both] <- FALSE
  d[outer(feas, !feas, "&")] <- TRUE
  inf <- outer(!feas, !feas, "&")
  d[inf] <- outer(cv, cv, "<")[inf]
  d
}

# Fast nondominated sorting: the rank (front number) of each point.
.opt_ranks <- function(d) {
  np <- nrow(d)
  rank <- integer(np)
  count <- colSums(d)
  r <- 0L
  left <- rep(TRUE, np)
  while (any(left)) {
    r <- r + 1L
    front <- which(left & count == 0L)
    if (!length(front)) front <- which(left)
    rank[front] <- r
    left[front] <- FALSE
    count <- count - colSums(d[front, , drop = FALSE])
  }
  rank
}

.opt_crowding <- function(Fm) {
  np <- nrow(Fm)
  cd <- numeric(np)
  if (np <= 2L) return(rep(Inf, np))
  for (k in seq_len(ncol(Fm))) {
    o <- order(Fm[, k])
    rng <- Fm[o[np], k] - Fm[o[1L], k]
    cd[o[c(1L, np)]] <- Inf
    if (rng > 0) cd[o[2:(np - 1L)]] <- cd[o[2:(np - 1L)]] + (Fm[o[3:np], k] - Fm[o[1:(np - 2L)], k]) / rng
  }
  cd
}

.opt_mo_fns <- function(a) {
  fn <- a$f
  p <- a$p
  m <- if (is.function(fn$cons)) length(fn$cons(as.numeric(a$u0), p)) else 0L
  lcons <- a$lcons %||% rep(0, m); ucons <- a$ucons %||% rep(0, m)
  list(F = function(x) as.numeric(fn$f(x, p)),
       cv = function(x) .opt_violation(fn, x, p, lcons, ucons),
       J = if (is.function(fn$jac)) function(x) fn$jac(x, p) else NULL)
}

.opt_nsga2 <- function(a, alg, opts) {
  n <- length(a$u0)
  lo <- a$lb; hi <- a$ub
  if (is.null(lo) || is.null(hi) || !all(is.finite(c(lo, hi))))
    stop("NSGA2() needs finite bounds on every unknown; give them with `bounds` metadata or `lb` and `ub`",
         call. = FALSE)
  if (!is.null(alg$seed)) {
    old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv())
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()),
            add = TRUE)
    set.seed(alg$seed)
  }
  fns <- .opt_mo_fns(a)
  N <- alg$population
  pm <- alg$mutation_rate %||% (1 / n)
  span <- hi - lo
  evalp <- function(X) {
    Fm <- do.call(rbind, lapply(seq_len(nrow(X)), function(i) fns$F(X[i, ])))
    Fm[!is.finite(Fm)] <- .Machine$double.xmax
    list(F = Fm, cv = apply(X, 1L, fns$cv))
  }
  X <- matrix(stats::runif(N * n), N, n) * rep(span, each = N) + rep(lo, each = N)
  X[1L, ] <- pmin(pmax(as.numeric(a$u0), lo), hi)
  ev <- evalp(X)
  Fm <- ev$F; cv <- ev$cv
  sbx <- function(x1, x2) {
    c1 <- x1; c2 <- x2
    for (j in seq_len(n)) {
      if (stats::runif(1L) > 0.5 || abs(x1[j] - x2[j]) < 1e-14) next
      y1 <- min(x1[j], x2[j]); y2 <- max(x1[j], x2[j])
      r <- stats::runif(1L)
      beta <- 1 + 2 * (y1 - lo[j]) / (y2 - y1)
      alpha <- 2 - beta^-(alg$eta_c + 1)
      bq <- if (r <= 1 / alpha) (r * alpha)^(1 / (alg$eta_c + 1)) else (1 / (2 - r * alpha))^(1 / (alg$eta_c + 1))
      a1 <- 0.5 * ((y1 + y2) - bq * (y2 - y1))
      beta <- 1 + 2 * (hi[j] - y2) / (y2 - y1)
      alpha <- 2 - beta^-(alg$eta_c + 1)
      bq <- if (r <= 1 / alpha) (r * alpha)^(1 / (alg$eta_c + 1)) else (1 / (2 - r * alpha))^(1 / (alg$eta_c + 1))
      a2 <- 0.5 * ((y1 + y2) + bq * (y2 - y1))
      a1 <- min(max(a1, lo[j]), hi[j]); a2 <- min(max(a2, lo[j]), hi[j])
      if (stats::runif(1L) <= 0.5) { c1[j] <- a2; c2[j] <- a1 } else { c1[j] <- a1; c2[j] <- a2 }
    }
    list(c1, c2)
  }
  mutate <- function(x) {
    for (j in seq_len(n)) {
      if (stats::runif(1L) > pm || span[j] <= 0) next
      d1 <- (x[j] - lo[j]) / span[j]; d2 <- (hi[j] - x[j]) / span[j]
      r <- stats::runif(1L)
      e <- 1 / (alg$eta_m + 1)
      dq <- if (r < 0.5) (2 * r + (1 - 2 * r) * (1 - d1)^(alg$eta_m + 1))^e - 1
            else 1 - (2 * (1 - r) + 2 * (r - 0.5) * (1 - d2)^(alg$eta_m + 1))^e
      x[j] <- min(max(x[j] + dq * span[j], lo[j]), hi[j])
    }
    x
  }
  for (g in seq_len(alg$generations)) {
    rank <- .opt_ranks(.opt_dominance(Fm, cv))
    crowd <- numeric(N)
    for (r in unique(rank)) crowd[rank == r] <- .opt_crowding(Fm[rank == r, , drop = FALSE])
    better <- function(i, j) rank[i] < rank[j] || (rank[i] == rank[j] && crowd[i] > crowd[j])
    pick <- function() {
      ij <- sample.int(N, 2L)
      if (better(ij[1L], ij[2L])) ij[1L] else ij[2L]
    }
    Q <- matrix(0, N, n)
    for (k in seq(1L, N, by = 2L)) {
      x1 <- X[pick(), ]; x2 <- X[pick(), ]
      ch <- if (stats::runif(1L) <= alg$crossover_rate) sbx(x1, x2) else list(x1, x2)
      Q[k, ] <- mutate(ch[[1L]])
      Q[k + 1L, ] <- mutate(ch[[2L]])
    }
    eq <- evalp(Q)
    RX <- rbind(X, Q); RF <- rbind(Fm, eq$F); Rcv <- c(cv, eq$cv)
    rr <- .opt_ranks(.opt_dominance(RF, Rcv))
    keep <- integer()
    for (r in sort(unique(rr))) {
      idx <- which(rr == r)
      if (length(keep) + length(idx) <= N) {
        keep <- c(keep, idx)
      } else {
        cdr <- .opt_crowding(RF[idx, , drop = FALSE])
        keep <- c(keep, idx[order(-cdr)][seq_len(N - length(keep))])
      }
      if (length(keep) >= N) break
    }
    X <- RX[keep, , drop = FALSE]; Fm <- RF[keep, , drop = FALSE]; cv <- Rcv[keep]
  }
  rank <- .opt_ranks(.opt_dominance(Fm, cv))
  sel <- which(rank == 1L & cv <= 0)
  retcode <- if (length(sel)) "Success" else "Infeasible"
  if (!length(sel)) sel <- which(rank == 1L)
  sel <- sel[!duplicated(round(Fm[sel, , drop = FALSE], 12))]
  sel <- sel[order(Fm[sel, 1L])]
  list(X = X[sel, , drop = FALSE], F = Fm[sel, , drop = FALSE], retcode = retcode,
       stats = list(generations = alg$generations, fevals = N * (alg$generations + 1L)))
}

.opt_weighted_scan <- function(a, alg, opts) {
  fns <- .opt_mo_fns(a)
  fn <- a$f
  x0 <- as.numeric(a$u0)
  k <- length(fns$F(x0))
  inner <- alg$inner
  # one scalarized problem: minimize sum(w * (f - shift) / scale)
  scalar <- function(w, shift, scale, start) {
    obj <- function(u, p) sum(w * (fns$F(u) - shift) / scale)
    grad <- if (!is.null(fns$J)) function(u, p) as.numeric(crossprod(fns$J(u), w / scale)) else NULL
    sp <- OptimizationProblem(obj, start, a$p, grad = grad, cons = fn$cons, cons_j = fn$cons_j,
                              lcons = a$lcons, ucons = a$ucons)
    sp$lb <- a$lb; sp$ub <- a$ub
    sol <- tryCatch(do.call(solve, c(list(sp, inner), opts)), error = function(e) NULL)
    if (is.null(sol) || !identical(sol$retcode, "Success")) return(NULL)
    u <- as.numeric(sol$u)
    list(u = u, f = fns$F(u), ok = fns$cv(u) <= 1e-6)
  }
  # the ideal point and an estimate of the nadir point from the anchors
  anchors <- lapply(seq_len(k), function(i) scalar(as.numeric(seq_len(k) == i), 0, 1, x0))
  ok <- !vapply(anchors, is.null, TRUE)
  if (!any(ok)) stop("WeightedSumScan(): every scalarized problem failed", call. = FALSE)
  AF <- do.call(rbind, lapply(anchors[ok], `[[`, "f"))
  ideal <- apply(AF, 2L, min); nadir <- apply(AF, 2L, max)
  scale <- nadir - ideal
  scale[!(scale > 0)] <- 1
  # weights on the simplex lattice with n points along each edge
  steps <- max(alg$points, 2L) - 1L
  grid <- as.matrix(expand.grid(rep(list(0:steps), k - 1L)))
  grid <- grid[rowSums(grid) <= steps, , drop = FALSE]
  W <- cbind(grid, steps - rowSums(grid)) / steps
  pts <- list(); start <- x0
  for (r in seq_len(nrow(W))) {
    s <- scalar(W[r, ], ideal, scale, start)
    if (is.null(s) || !s$ok) next
    pts[[length(pts) + 1L]] <- s
    start <- s$u
  }
  pts <- c(pts, Filter(function(s) !is.null(s) && s$ok, anchors))
  if (!length(pts)) stop("WeightedSumScan(): no feasible point was found", call. = FALSE)
  X <- do.call(rbind, lapply(pts, `[[`, "u")); Fm <- do.call(rbind, lapply(pts, `[[`, "f"))
  keep <- which(colSums(.opt_dominance(Fm, rep(0, nrow(Fm)))) == 0L)
  keep <- keep[!duplicated(round(Fm[keep, , drop = FALSE], 8))]
  keep <- keep[order(Fm[keep, 1L])]
  list(X = X[keep, , drop = FALSE], F = Fm[keep, , drop = FALSE], retcode = "Success",
       stats = list(nsolves = nrow(W) + k))
}
