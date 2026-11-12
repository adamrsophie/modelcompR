# OptimizationIpopt and OptimizationMadNLP: IpoptOptimizer() and
# MadNLPOptimizer() run the primal dual interior point method of
# R/dynamic-opt-nlp.R (Waechter and Biegler's filter line search algorithm,
# the algorithm of Ipopt and of MadNLP) on a general OptimizationProblem.
# Derivatives the problem does not supply are taken by finite differences.

.ipm_alg <- function(name, fields) {
  structure(c(list(name = name), fields, list(traits = list(callback = TRUE, bounds = TRUE, cons = TRUE),
                                              solve = .ipm_solve)),
            class = c(name, "OptimizationAlgorithm"))
}

#' Interior point optimizers
#'
#' `IpoptOptimizer()` (OptimizationIpopt.jl) and `MadNLPOptimizer()`
#' (OptimizationMadNLP.jl) solve
#' \deqn{\min f(u) \quad \textrm{subject to} \quad lcons \le cons(u) \le ucons,
#' \; lb \le u \le ub}{min f(u) subject to lcons <= cons(u) <= ucons, lb <= u <= ub}
#' with the primal dual interior point method and filter line search of
#' Waechter and Biegler (2006), the algorithm of both solvers: a logarithmic
#' barrier for the bounds and for slacks of the inequalities, the monotone
#' barrier update, gradient based scaling of the problem, inertia free
#' regularization of the Newton system and second order corrections. The
#' restoration phase of the solvers is replaced by a minimum norm step toward
#' feasibility, so iterates and counts differ from Ipopt and MadNLP while
#' the solutions agree.
#'
#' Gradients, Hessians and constraint Jacobians and Hessians come from the
#' problem when it has them and from finite differences otherwise.
#' `reltol` (Ipopt) or `abstol` (MadNLP) sets the convergence tolerance
#' (default `1e-8`) and `maxiters` the iteration limit (default 3000); for
#' Ipopt, `tol` and `max_iter` in `additional_options` take precedence, as in
#' OptimizationIpopt.jl. The other fields of the Julia structs (linear
#' solvers, printing, scaling, barrier strategy, acceptable point
#' tolerances) are accepted and stored but do not change the native method,
#' which uses Ipopt's defaults for them.
#'
#' @param ... Fields of the Julia structs `IpoptOptimizer` (for example
#'   `acceptable_tol`, `mu_init`, `linear_solver`, `hessian_approximation`,
#'   `additional_options`) and `MadNLPOptimizer` (for example
#'   `acceptable_tol`, `linear_solver`, `barrier`, `additional_options`),
#'   with the Julia defaults.
#' @return An algorithm object for [solve()].
#' @examples
#' prob <- OptimizationProblem(function(u, p) (u[1] - 2) * (u[1] - 2) + (u[2] - 1) * (u[2] - 1),
#'                             c(0, 0), numeric(), cons = function(u, p) u[1] * u[1] + u[2] * u[2],
#'                             lcons = -Inf, ucons = 1)
#' solve(prob, IpoptOptimizer())$u
#' @name interior_point_optimizers
NULL

.ipopt_fields <- list(
  acceptable_tol = 1e-6, acceptable_iter = 15L, dual_inf_tol = 1, constr_viol_tol = 1e-4, compl_inf_tol = 1e-4,
  print_timing_statistics = "no", print_info_string = "no", linear_solver = "mumps", linear_system_scaling = "none",
  hsllib = "", pardisolib = "", linear_scaling_on_demand = "yes", nlp_scaling_method = "gradient-based",
  nlp_scaling_max_gradient = 100, honor_original_bounds = "no", check_derivatives_for_naninf = "no",
  mu_strategy = "monotone", mu_oracle = "quality-function", mu_init = 0.1,
  adaptive_mu_globalization = "obj-constr-filter", warm_start_init_point = "no", hessian_approximation = "exact",
  limited_memory_max_history = 6L, limited_memory_update_type = "bfgs", accept_every_trial_step = "no",
  line_search_method = "filter", expect_infeasible_problem = "no", additional_options = list())

.madnlp_fields <- list(
  rethrow_error = TRUE, disable_garbage_collector = FALSE, blas_num_threads = 1L, output_file = "",
  file_print_level = "INFO", acceptable_tol = 1e-6, acceptable_iter = 15L, jacobian_constant = FALSE,
  hessian_constant = FALSE, hessian_approximation = "ExactHessian", nlp_scaling = TRUE,
  nlp_scaling_max_gradient = 100, linear_solver = NULL, kkt_system = NULL, barrier = NULL,
  quasi_newton_options = NULL, additional_options = list())

.ipm_fields <- function(cls, defaults, args) {
  bad <- setdiff(names(args), names(defaults))
  if (length(bad) || (length(args) && any(!nzchar(names(args) %||% ""))))
    stop(cls, "() has no field ", paste(bad, collapse = ", "), call. = FALSE)
  utils::modifyList(defaults, args, keep.null = TRUE)
}

#' @rdname interior_point_optimizers
#' @export
IpoptOptimizer <- function(...) .ipm_alg("IpoptOptimizer", .ipm_fields("IpoptOptimizer", .ipopt_fields, list(...)))

#' @rdname interior_point_optimizers
#' @export
MadNLPOptimizer <- function(...) .ipm_alg("MadNLPOptimizer", .ipm_fields("MadNLPOptimizer", .madnlp_fields, list(...)))

# central difference Jacobian of a vector function
.ipm_fd_jac <- function(g, x, gx = g(x)) {
  vapply(seq_along(x), function(k) {
    h <- .Machine$double.eps^(1 / 3) * max(1, abs(x[k]))
    (g(replace(x, k, x[k] + h)) - g(replace(x, k, x[k] - h))) / (2 * h)
  }, numeric(length(gx)))
}

.ipm_triplets <- function(M) {
  M <- as.matrix(M)
  list(i = as.integer(row(M)), j = as.integer(col(M)), v = as.numeric(M))
}

.ipm_solve <- function(a, alg, opts) {
  fn <- a$f
  p <- a$p
  x0 <- as.numeric(a$u0)
  n <- length(x0)
  ipopt <- inherits(alg, "IpoptOptimizer")
  ao <- alg$additional_options %||% list()
  tol <- if (ipopt) ao$tol %||% opts$reltol %||% 1e-8 else opts$abstol %||% 1e-8
  maxiter <- as.integer(if (ipopt) ao$max_iter %||% opts$maxiters %||% 3000L else opts$maxiters %||% 3000L)
  f <- function(x) as.numeric(fn$f(x, p))[1L]
  grad <- if (is.function(fn$grad)) function(x) as.numeric(fn$grad(x, p)) else
    function(x) as.numeric(.ipm_fd_jac(function(z) f(z), x, f(x)))
  hess <- if (is.function(fn$hess)) function(x) as.matrix(fn$hess(x, p)) else function(x) {
    H <- .ipm_fd_jac(grad, x)
    (H + t(H)) / 2
  }
  m <- length(a$lcons %||% numeric())
  cons <- if (m) function(x) as.numeric(fn$cons(x, p)) else function(x) numeric()
  cjac <- if (m && is.function(fn$cons_j)) function(x) matrix(fn$cons_j(x, p), m, n) else
    function(x) matrix(.ipm_fd_jac(cons, x), m, n)
  chess <- if (m && is.function(fn$cons_h)) function(x) lapply(fn$cons_h(x, p), as.matrix) else function(x) {
    J <- function(z) as.numeric(cjac(z))
    D <- .ipm_fd_jac(J, x)
    lapply(seq_len(m), function(i) {
      Hi <- D[i + m * (seq_len(n) - 1L), , drop = FALSE]
      (Hi + t(Hi)) / 2
    })
  }
  lc <- rep_len(as.numeric(a$lcons %||% numeric()), m)
  uc <- rep_len(as.numeric(a$ucons %||% numeric()), m)
  # equalities cE = c - lcons; inequalities c - ucons <= 0 and lcons - c <= 0
  ie <- which(lc == uc)
  iu <- which(lc != uc & is.finite(uc))
  il <- which(lc != uc & is.finite(lc))
  sgnI <- c(rep(1, length(iu)), rep(-1, length(il)))
  rowsI <- c(iu, il)
  evalf <- function(z, level, lamE = NULL, lamI = NULL, obj = 1) {
    cv <- cons(z)
    out <- list(f = f(z), cE = cv[ie] - lc[ie], cI = sgnI * (cv[rowsI] - c(uc[iu], lc[il])))
    if (level >= 1L) {
      J <- cjac(z)
      out$gf <- grad(z)
      out$JE <- .ipm_triplets(J[ie, , drop = FALSE])
      out$JI <- .ipm_triplets(sgnI * J[rowsI, , drop = FALSE])
    }
    if (level >= 2L) {
      H <- obj * hess(z)
      if (m) {
        w <- numeric(m)
        w[ie] <- w[ie] + (lamE %||% numeric(length(ie)))
        if (length(rowsI)) {
          lI <- lamI %||% numeric(length(rowsI))
          for (k in seq_along(rowsI)) w[rowsI[k]] <- w[rowsI[k]] + sgnI[k] * lI[k]
        }
        if (any(w != 0)) {
          Hc <- chess(z)
          for (i in which(w != 0)) H <- H + w[i] * Hc[[i]]
        }
      }
      out$H <- .ipm_triplets(H)
    }
    out
  }
  lo <- if (is.null(a$lb)) rep(-Inf, n) else rep_len(as.numeric(a$lb), n)
  hi <- if (is.null(a$ub)) rep(Inf, n) else rep_len(as.numeric(a$ub), n)
  # one dense block: every unknown and row carries the same key
  r <- .nlp_interior_point(evalf, x0, lo, hi, rep(0, n), rep(0, length(ie)), rep(0, length(rowsI)), tol = tol,
                           maxiter = maxiter, verbose = FALSE)
  retcode <- switch(r$status, Success = "Success", MaxIters = "MaxIters", Infeasible = , Stalled = "Infeasible",
                    "Failure")
  lambda <- numeric(m)
  lambda[ie] <- r$lamE
  if (length(rowsI)) for (k in seq_along(rowsI)) lambda[rowsI[k]] <- lambda[rowsI[k]] + sgnI[k] * r$lamI[k]
  original <- list(status = r$status, iterations = r$iterations, nit = r$iterations, lambda = lambda,
                   constraint_violation = r$violation)
  .opt_solution(a, alg, r$z, r$f, retcode, original)
}
