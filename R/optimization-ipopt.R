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
#' Waechter and Biegler (2006), ported from Ipopt 3.14: a logarithmic
#' barrier for the bounds and for slacks of the inequalities, the monotone
#' or adaptive barrier update, gradient based scaling of the problem,
#' Ipopt's inertia correction of the Newton system, the backtracking filter
#' line search with second order corrections, the watchdog and the soft
#' restoration phase, and Ipopt's restoration phase, which minimizes the
#' 1-norm of the infeasibility when the line search fails.
#'
#' Gradients, Hessians and constraint Jacobians and Hessians come from the
#' problem when it has them and from finite differences otherwise.
#' `reltol` (Ipopt) or `abstol` (MadNLP) sets the convergence tolerance
#' (default `1e-8`) and `maxiters` the iteration limit (default 3000); for
#' Ipopt, `tol` and `max_iter` in `additional_options` take precedence, as in
#' OptimizationIpopt.jl. The fields of `IpoptOptimizer` and Ipopt options in
#' `additional_options` act as in Ipopt (tolerances, scaling,
#' `mu_strategy`, `mu_oracle`, `hessian_approximation`,
#' `warm_start_init_point` and the others); a linear solver or scaling
#' library that Ipopt_jll lacks (the HSL solvers, MC19) ends the solve
#' before the first iteration, as Ipopt's Invalid_Option does. Options that
#' only change printing or the linear algebra library are accepted without
#' effect on the iterates. The penalty line searches of Ipopt
#' (`line_search_method = "penalty"` or `"cg-penalty"`) are not ported: they
#' run the filter line search. `MadNLPOptimizer()` runs the same Ipopt
#' algorithm, not MadNLP's own implementation.
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

# Ipopt's options that only change printing, timing or the linear algebra
# library: accepted with no effect on the iterates.
.ipopt_quiet_options <- c(
  "print_level", "print_user_options", "print_options_documentation", "print_options_mode", "output_file",
  "file_print_level", "file_append", "print_frequency_iter", "print_frequency_time", "print_timing_statistics",
  "print_info_string", "print_advanced_options", "print_ls_trace", "sb", "timing_statistics", "hsllib",
  "pardisolib", "linear_scaling_on_demand", "linear_system_scaling", "mumps_pivtol", "mumps_pivtolmax",
  "mumps_mem_percent", "mumps_permuting_scaling", "mumps_pivot_order", "mumps_scaling", "mumps_dep_tol",
  "max_refinement_steps", "min_refinement_steps", "residual_ratio_max", "residual_ratio_singular",
  "residual_improvement_factor", "check_derivatives_for_naninf", "nlp_lower_bound_inf",
  "nlp_upper_bound_inf", "max_cpu_time", "derivative_test", "skip_finalize_solution_call", "inf_pr_output",
  "mu_oracle", "adaptive_mu_globalization", "limited_memory_max_history", "limited_memory_update_type",
  "mu_strategy", "hessian_approximation", "line_search_method", "warm_start_init_point", "linear_solver",
  "fixed_variable_treatment", "alpha_for_y", "recalc_y")

# Linear solvers of an Ipopt_jll build without HSL, Pardiso or WSMP: the
# others make Ipopt stop with Invalid_Option.
.ipopt_linear_solvers <- c("mumps", "spral")

# The options of the native method from the fields of IpoptOptimizer or
# MadNLPOptimizer and their additional_options; NULL for an option that
# Ipopt rejects (an invalid value or a library that is not available).
.ipm_options <- function(alg, ipopt) {
  ao <- alg$additional_options %||% list()
  get <- function(k) if (!is.null(ao[[k]])) ao[[k]] else alg[[k]]
  yes <- function(v) identical(v, "yes") || isTRUE(v)
  defaults <- .nlp_ip_defaults()
  o <- list(variant = if (ipopt) "ipopt" else "madnlp")
  for (k in c("acceptable_tol", "acceptable_iter", "dual_inf_tol", "constr_viol_tol", "compl_inf_tol",
              "mu_init", "nlp_scaling_max_gradient"))
    if (!is.null(get(k))) o[[k]] <- get(k)
  if (!ipopt) {
    o$nlp_scaling_method <- if (isTRUE(alg$nlp_scaling)) "gradient-based" else "none"
    for (k in intersect(names(ao), names(defaults))) o[[k]] <- ao[[k]]
    return(o)
  }
  sm <- get("nlp_scaling_method")
  # equilibration-based scaling needs MC19 of HSL
  if (!sm %in% c("none", "user-scaling", "gradient-based")) return(NULL)
  o$nlp_scaling_method <- sm
  o$honor_original_bounds <- yes(get("honor_original_bounds"))
  o$accept_every_trial_step <- yes(get("accept_every_trial_step"))
  o$expect_infeasible_problem <- yes(get("expect_infeasible_problem"))
  o$warm_start_init_point <- yes(get("warm_start_init_point"))
  if (!get("linear_solver") %in% .ipopt_linear_solvers) return(NULL)
  if (!get("linear_system_scaling") %in% c("none", "slack-based")) return(NULL)
  if (!get("mu_strategy") %in% c("monotone", "adaptive")) return(NULL)
  if (!get("hessian_approximation") %in% c("exact", "limited-memory")) return(NULL)
  if (!get("line_search_method") %in% c("filter", "penalty", "cg-penalty")) return(NULL)
  if (!get("mu_oracle") %in% c("quality-function", "probing", "loqo")) return(NULL)
  if (!get("adaptive_mu_globalization") %in% c("obj-constr-filter", "kkt-error", "never-monotone-mode"))
    return(NULL)
  o$mu_strategy <- get("mu_strategy")
  o$mu_oracle <- get("mu_oracle")
  o$adaptive_mu_globalization <- get("adaptive_mu_globalization")
  o$hessian_approximation <- get("hessian_approximation")
  o$limited_memory_max_history <- as.integer(get("limited_memory_max_history"))
  o$limited_memory_update_type <- get("limited_memory_update_type")
  o$obj_scaling_factor <- as.numeric(ao$obj_scaling_factor %||% 1)
  # the restoration phase starts from zero constraint multipliers unless
  # constr_mult_init_max is set (RestoIterateInitializer)
  if (!is.null(ao$constr_mult_init_max)) o$resto_constr_mult_init_max <- ao$constr_mult_init_max
  for (k in names(ao)) {
    if (k %in% c("tol", "max_iter", "obj_scaling_factor") || k %in% names(alg)) next
    if (k %in% names(defaults)) {
      v <- ao[[k]]
      o[[k]] <- if (is.character(v) && v %in% c("yes", "no")) v == "yes" else v
    } else if (!k %in% .ipopt_quiet_options) {
      return(NULL)
    }
  }
  o
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
  io <- .ipm_options(alg, ipopt)
  if (is.null(io)) {
    # Ipopt stops before the first iteration (Invalid_Option)
    original <- list(status = "InvalidOption", iterations = 0L, nit = 0L,
                     lambda = rep(0, length(a$lcons %||% numeric())), constraint_violation = NA_real_)
    return(.opt_solution(a, alg, x0, as.numeric(fn$f(x0, p))[1L], "Failure", original))
  }
  if (!is.null(opts$maxtime)) io$max_wall_time <- as.numeric(opts$maxtime)
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
                           maxiter = maxiter, verbose = isTRUE(opts$verbose), opts = io,
                           rowI_bound = c(uc[iu], lc[il]))
  # the ReturnCode of OptimizationIpopt.jl for Ipopt's status, and of
  # OptimizationMadNLP.jl (map_madnlp_status) for MadNLP's
  retcode <- if (ipopt) {
    switch(r$status, Success = , Acceptable = , FeasiblePointFound = "Success",
           Infeasible = , TinyStep = , Diverging = "Infeasible", MaxIters = "MaxIters",
           MaxTime = "MaxTime", "Failure")
  } else {
    switch(r$status, Success = , Acceptable = "Success",
           Infeasible = , TinyStep = , Diverging = , RestorationFailed = , TooFewDof = "Infeasible",
           MaxIters = "MaxIters", MaxTime = "MaxTime", "Failure")
  }
  lambda <- numeric(m)
  lambda[ie] <- r$lamE
  if (length(rowsI)) for (k in seq_along(rowsI)) lambda[rowsI[k]] <- lambda[rowsI[k]] + sgnI[k] * r$lamI[k]
  original <- list(status = r$status, iterations = r$iterations, nit = r$iterations, lambda = lambda,
                   constraint_violation = r$violation, restorations = r$restorations)
  .opt_solution(a, alg, r$z, r$f, retcode, original)
}
