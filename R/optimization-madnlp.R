# MadNLP's interior point method (MadNLP.jl 0.11.0, src/IPM and
# src/Callbacks) for MadNLPOptimizer (OptimizationMadNLP.jl): the
# reformulation of the problem by SparseCallback (gradient based scaling,
# slacks for the inequalities, fixed variables made parameters, relaxed
# bounds), the regular iterations of `regular!` with MadNLP's filter line
# search, the soft restoration `restore!` and the robust restoration
# `robust!`, the barrier updates MonotoneUpdate, QualityFunctionUpdate and
# LOQOUpdate, the inertia correctors InertiaBased, InertiaFree and
# InertiaIgnore, Richardson's iterative refinement and the compact limited
# memory BFGS approximation of the Hessian. The KKT system is the
# augmented system of SparseKKTSystem, factorized densely; its inertia
# comes from the Bunch-Kaufman factorization (.nlp_inertia), as MUMPS
# reports it.

# MadNLPOptions (src/IPM/options.jl) without the primary options.
.madnlp_ipm_defaults <- function() {
  list(rethrow_error = TRUE, disable_garbage_collector = FALSE, blas_num_threads = 1L,
       iterator = "RichardsonIterator", intermediate_callback = NULL, output_file = "", print_level = "INFO",
       file_print_level = "INFO", acceptable_tol = 1e-6, acceptable_iter = 15L, diverging_iterates_tol = 1e20,
       max_iter = 3000L, max_wall_time = 1e6, s_max = 100, kappa_d = 1e-5, fixed_variable_treatment = NULL,
       equality_treatment = NULL, bound_relax_factor = 1e-8, jacobian_constant = FALSE,
       hessian_constant = FALSE, hessian_approximation = "ExactHessian", quasi_newton_options = NULL,
       inertia_correction_method = "InertiaAuto", inertia_free_tol = 0, default_primal_regularization = 0,
       default_dual_regularization = 0, dual_initialized = FALSE, dual_initialization_method = NULL,
       constr_mult_init_max = 1e3, bound_push = 1e-2, bound_fac = 1e-2, nlp_scaling = TRUE,
       nlp_scaling_max_gradient = 100, min_hessian_perturbation = 1e-20, first_hessian_perturbation = 1e-4,
       max_hessian_perturbation = 1e20, perturb_inc_fact_first = 1e2, perturb_inc_fact = 8,
       perturb_dec_fact = 1 / 3, jacobian_regularization_exponent = 1 / 4,
       jacobian_regularization_value = 1e-8, soft_resto_pderror_reduction_factor = 0.9999,
       required_infeasibility_reduction = 0.9, obj_max_inc = 5, max_soc = 4L, alpha_min_frac = 0.05,
       s_theta = 1.1, s_phi = 2.3, eta_phi = 1e-4, kappa_soc = 0.99, gamma_theta = 1e-5, gamma_phi = 1e-5,
       delta = 1, kappa_sigma = 1e10, barrier_tol_factor = 10, rho = 1000, barrier = NULL, tau_min = 0.99,
       kkt_options = list())
}

# The options of the linear solvers and of the iterative refinement
# (MumpsOptions, LapackOptions, UmfpackOptions, LDLFactorizationsOptions,
# RichardsonOptions); only the Richardson options change the iterates.
.madnlp_solver_options <- c(
  "mumps_mem_percent", "mumps_permuting_scaling", "mumps_pivot_order", "mumps_pivtol", "mumps_pivtolmax",
  "mumps_scaling", "lapack_algorithm", "umfpack_pivtol", "umfpack_pivtolmax", "umfpack_sym_pivtol",
  "umfpack_block_size", "umfpack_strategy", "ldl_tol", "cholmod_algorithm", "richardson_max_iter",
  "richardson_tol", "richardson_acceptable_tol")

# The barrier updates of src/IPM/types.jl: a name ("MonotoneUpdate",
# "QualityFunctionUpdate", "LOQOUpdate") or a list with `type` and fields.
# NULL is MadNLP's default, MonotoneUpdate(tol, barrier_tol_factor).
.madnlp_barrier <- function(b, tol, barrier_tol_factor) {
  defaults <- list(
    MonotoneUpdate = list(mu_init = 1e-1, mu_min = 1e-11, mu_superlinear_decrease_power = 1.5,
                          mu_linear_decrease_factor = 0.2),
    QualityFunctionUpdate = list(mu_init = 1e-1, mu_min = 1e-11, mu_max = 1e5, sigma_min = 1e-6, sigma_max = 1e2,
                                 sigma_tol = 1e-2, gamma = 1, max_gs_iter = 8L,
                                 mu_superlinear_decrease_power = 1.5, mu_linear_decrease_factor = 0.2,
                                 free_mode = TRUE, globalization = TRUE, n_update = 0L),
    LOQOUpdate = list(mu_init = 1e-1, mu_min = 1e-11, mu_max = 1e5, gamma = 0.1, r = 0.95,
                      mu_superlinear_decrease_power = 1.5, mu_linear_decrease_factor = 0.2, free_mode = TRUE,
                      globalization = TRUE))
  if (is.null(b)) {
    out <- defaults$MonotoneUpdate
    out$mu_min <- min(1e-4, tol) / (barrier_tol_factor + 1)
    out$type <- "MonotoneUpdate"
    return(out)
  }
  if (is.character(b)) b <- list(type = b)
  type <- sub("^MadNLP\\.", "", b$type %||% "")
  if (!type %in% names(defaults))
    stop("unknown MadNLP barrier update ", b$type, "; use MonotoneUpdate, QualityFunctionUpdate or LOQOUpdate",
         call. = FALSE)
  out <- utils::modifyList(defaults[[type]], b[setdiff(names(b), "type")])
  out$type <- type
  out
}

# MadNLP's Status values and their ReturnCode (map_madnlp_status of
# OptimizationMadNLP.jl).
.madnlp_retcode <- function(status) {
  switch(status, SOLVE_SUCCEEDED = , SOLVED_TO_ACCEPTABLE_LEVEL = , USER_REQUESTED_STOP = "Success",
         INFEASIBLE_PROBLEM_DETECTED = , SEARCH_DIRECTION_BECOMES_TOO_SMALL = , DIVERGING_ITERATES = ,
         RESTORATION_FAILED = , NOT_ENOUGH_DEGREES_OF_FREEDOM = "Infeasible",
         MAXIMUM_ITERATIONS_EXCEEDED = "MaxIters", MAXIMUM_WALLTIME_EXCEEDED = "MaxTime", "Failure")
}

# MadNLPOptimizer on an OptimizationProblem: the options of
# __map_optimizer_args (tol from abstol, max_iter from maxiters, max_wall_time
# from maxtime, the fields and additional_options) and MadNLP's solve.
.madnlp_solve <- function(a, alg, opts) {
  fn <- a$f
  p <- a$p
  x0 <- as.numeric(a$u0)
  n <- length(x0)
  ao <- alg$additional_options %||% list()
  o <- .madnlp_ipm_defaults()
  ls_opt <- list()
  ignored <- character()
  for (k in names(ao)) {
    if (k %in% c("tol", "linear_solver", "callback", "kkt_system")) {
      o[[k]] <- ao[[k]]
    } else if (k %in% names(o)) {
      o[k] <- list(ao[[k]])
    } else if (k %in% .madnlp_solver_options) {
      ls_opt[[k]] <- ao[[k]]
    } else {
      ignored <- c(ignored, k)
    }
  }
  if (length(ignored))
    warning("The following options are ignored: ", paste(ignored, collapse = ", "), call. = FALSE)
  if (!is.null(alg$barrier)) o$barrier <- alg$barrier
  if (!is.null(alg$quasi_newton_options)) o$quasi_newton_options <- alg$quasi_newton_options
  if (!is.null(alg$linear_solver)) o$linear_solver <- alg$linear_solver
  if (!is.null(alg$kkt_system)) o$kkt_system <- alg$kkt_system
  for (k in c("rethrow_error", "disable_garbage_collector", "blas_num_threads", "output_file",
              "file_print_level", "acceptable_tol", "acceptable_iter", "jacobian_constant", "hessian_constant",
              "hessian_approximation", "nlp_scaling", "nlp_scaling_max_gradient"))
    o[[k]] <- alg[[k]]
  o$tol <- opts$abstol %||% 1e-8
  o$max_iter <- as.integer(opts$maxiters %||% 3000L)
  o$max_wall_time <- as.numeric(opts$maxtime %||% 1e6)
  o$ls_opt <- ls_opt
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
  # the Hessian of the Lagrangian obj_weight f + y' c (hess_coord! of
  # OptimizationNLPModels)
  lag_hess <- function(x, w, y) {
    H <- w * hess(x)
    if (m && any(y != 0)) {
      Hc <- chess(x)
      for (i in seq_len(m)) H <- H + y[i] * Hc[[i]]
    }
    H
  }
  pb <- list(n = n, m = m, x0 = x0, f = f, grad = grad, cons = cons, jac = cjac, lag_hess = lag_hess,
             lvar = if (is.null(a$lb)) rep(-Inf, n) else rep_len(as.numeric(a$lb), n),
             uvar = if (is.null(a$ub)) rep(Inf, n) else rep_len(as.numeric(a$ub), n),
             lcon = rep_len(as.numeric(a$lcons %||% numeric()), m),
             ucon = rep_len(as.numeric(a$ucons %||% numeric()), m))
  r <- .madnlp_ipm(pb, o, verbose = isTRUE(opts$verbose))
  original <- list(status = r$status, iterations = r$iterations, nit = r$iterations, lambda = r$lambda,
                   counters = list(k = r$iterations), multipliers = r$lambda,
                   multipliers_L = r$zl, multipliers_U = r$zu, constraint_violation = r$violation)
  .opt_solution(a, alg, r$x, r$objective, .madnlp_retcode(r$status), original)
}

# MadNLPSolver and solve! for the problem `pb` (dense callbacks of an
# AbstractNLPModel: f, grad, cons, jac, lag_hess, bounds) with the options
# `o`. Returns the status, the solution and the multipliers in the
# unscaled problem (MadNLPExecutionStats).
.madnlp_ipm <- function(pb, o, verbose = FALSE) {
  start_time <- proc.time()[["elapsed"]]
  eps <- .Machine$double.eps
  tol <- as.numeric(o$tol)
  # the presets of MadNLPOptions{T}(nlp) for a sparse model
  kkt_system <- sub("^MadNLP\\.", "", o$kkt_system %||% "SparseKKTSystem")
  linear_solver <- sub("^MadNLP\\.", "", o$linear_solver %||% "MumpsSolver")
  known_kkt <- c("SparseKKTSystem", "SparseUnreducedKKTSystem", "ScaledSparseKKTSystem",
                 "SparseCondensedKKTSystem", "DenseKKTSystem", "DenseCondensedKKTSystem")
  if (!kkt_system %in% known_kkt) stop("unknown MadNLP KKT system ", kkt_system, call. = FALSE)
  known_ls <- c("MumpsSolver", "LapackCPUSolver", "UmfpackSolver", "LDLSolver")
  if (!linear_solver %in% known_ls)
    stop("the MadNLP linear solver ", linear_solver, " is not available; use one of ",
         paste(known_ls, collapse = ", "), call. = FALSE)
  dense_kkt <- kkt_system %in% c("DenseKKTSystem", "DenseCondensedKKTSystem")
  if (dense_kkt && linear_solver != "LapackCPUSolver")
    stop("[options] Sparse Linear solver is not supported in dense mode.\n",
         "Please use a dense linear solver or change `kkt_system` ", call. = FALSE)
  condensed <- kkt_system == "SparseCondensedKKTSystem"
  hess_approx <- sub("^MadNLP\\.", "", o$hessian_approximation %||% "ExactHessian")
  if (!hess_approx %in% c("ExactHessian", "CompactLBFGS", "BFGS", "DampedBFGS"))
    stop("unknown MadNLP Hessian approximation ", hess_approx, call. = FALSE)
  if (hess_approx %in% c("BFGS", "DampedBFGS") && !dense_kkt)
    stop("[options] BFGS and DampedBFGS quasi-Newton approximations\n",
         "require a dense KKT system (DenseKKTSystem or DenseCondensedKKTSystem).", call. = FALSE)
  if (hess_approx != "ExactHessian" && hess_approx != "CompactLBFGS")
    stop("MadNLP's ", hess_approx, " approximation is not ported; use ExactHessian or CompactLBFGS",
         call. = FALSE)
  fixed_treatment <- sub("^MadNLP\\.", "", o$fixed_variable_treatment %||%
                           if (condensed) "RelaxBound" else "MakeParameter")
  eq_treatment <- sub("^MadNLP\\.", "", o$equality_treatment %||%
                        if (condensed) "RelaxEquality" else "EnforceEquality")
  dual_init <- sub("^MadNLP\\.", "", o$dual_initialization_method %||%
                     if (condensed) "DualInitializeSetZero" else "DualInitializeLeastSquares")
  barrier <- .madnlp_barrier(o$barrier, tol, o$barrier_tol_factor)
  adaptive <- barrier$type != "MonotoneUpdate"
  ic_method <- sub("^MadNLP\\.", "", o$inertia_correction_method %||% "InertiaAuto")
  if (ic_method == "InertiaAuto") ic_method <- if (linear_solver == "UmfpackSolver") "InertiaFree" else "InertiaBased"
  rich_max <- as.integer(o$ls_opt$richardson_max_iter %||% 10L)
  rich_tol <- o$ls_opt$richardson_tol %||% tol^(5 / 4)
  rich_acc <- o$ls_opt$richardson_acceptable_tol %||% tol^(5 / 8)
  s_max <- o$s_max

  # SparseCallback: fixed variables, equalities and inequalities
  nvar <- pb$n
  m <- pb$m
  lvar0 <- pb$lvar; uvar0 <- pb$uvar
  isfixed <- lvar0 == uvar0
  make_param <- fixed_treatment == "MakeParameter" && any(isfixed)
  free <- if (make_param) which(lvar0 < uvar0) else seq_len(nvar)
  n <- length(free)
  x_full <- pb$x0
  if (make_param) x_full[isfixed] <- lvar0[isfixed]
  lcon <- pb$lcon; ucon <- pb$ucon
  if (m > 0 && eq_treatment == "EnforceEquality") {
    ind_ineq <- which(lcon != ucon)
  } else {
    ind_ineq <- seq_len(m)
  }
  ns <- length(ind_ineq)
  N <- n + ns
  lvar <- lvar0[free]; uvar <- uvar0[free]
  xl <- c(lvar, lcon[ind_ineq]); xu <- c(uvar, ucon[ind_ineq])
  ind_llb <- which(lvar != -Inf & uvar == Inf)
  ind_uub <- which(lvar == -Inf & uvar != Inf)
  ind_lb <- which(xl != -Inf)
  ind_ub <- which(xu != Inf)
  nlb <- length(ind_lb); nub <- length(ind_ub)
  n_eq <- if (eq_treatment == "EnforceEquality") sum(lcon == ucon) else 0L
  iP <- seq_len(N); iY <- N + seq_len(m); iL <- N + m + seq_len(nlb); iU <- N + m + nlb + seq_len(nub)
  nd <- N + m + nlb + nub
  slk <- n + seq_len(ns)

  status <- "INITIAL"
  out_x <- function(x) {
    xf <- x_full
    xf[free] <- x[seq_len(n)]
    xf
  }
  if (n == 0L) {
    return(list(status = status, iterations = 0L, x = x_full, objective = pb$f(x_full), lambda = numeric(m),
                zl = numeric(nvar), zu = numeric(nvar), violation = NA_real_))
  }

  obj_scale <- 1
  con_scale <- rep(1, m)
  cnt <- list(k = 0L, l = 0L, t = 0L, ir = 0L, acceptable_cnt = 0L, unsuccessful_iterate = 0L,
              restoration_fail_count = 0L, obj = 0L, grad = 0L, con = 0L, jac = 0L, hess = 0L)
  invalid <- function(what) stop(structure(class = c("madnlp_invalid", "error", "condition"),
                                           list(message = what, call = NULL, what = what)))
  # the callbacks of the scaled problem (src/IPM/callbacks.jl)
  rhs <- numeric(m)
  eval_f <- function(x) {
    v <- pb$f(out_x(x)) * obj_scale
    cnt$obj <<- cnt$obj + 1L
    if (cnt$obj == 1L && !is.finite(v)) invalid("INVALID_NUMBER_OBJECTIVE")
    v
  }
  eval_grad <- function(x) {
    g <- pb$grad(out_x(x))[free] * obj_scale
    cnt$grad <<- cnt$grad + 1L
    if (cnt$grad == 1L && !all(is.finite(g))) invalid("INVALID_NUMBER_GRADIENT")
    c(g, numeric(ns))
  }
  eval_cons <- function(x) {
    cv <- pb$cons(out_x(x)) * con_scale
    cv[ind_ineq] <- cv[ind_ineq] - x[slk]
    cv <- cv - rhs
    cnt$con <<- cnt$con + 1L
    if (cnt$con == 1L && !all(is.finite(cv))) invalid("INVALID_NUMBER_CONSTRAINTS")
    cv
  }
  eval_jac <- function(x) {
    J <- matrix(0, m, N)
    if (m) {
      J[, seq_len(n)] <- pb$jac(out_x(x))[, free, drop = FALSE] * con_scale
      if (ns) J[cbind(ind_ineq, slk)] <- -1
    }
    cnt$jac <<- cnt$jac + 1L
    if (cnt$jac == 1L && !all(is.finite(J))) invalid("INVALID_NUMBER_JACOBIAN")
    J
  }
  eval_hess <- function(x, y, w = 1) {
    H <- matrix(0, N, N)
    H[seq_len(n), seq_len(n)] <- pb$lag_hess(out_x(x), w * obj_scale, y * con_scale)[free, free, drop = FALSE]
    cnt$hess <<- cnt$hess + 1L
    if (cnt$hess == 1L && !all(is.finite(H))) invalid("INVALID_NUMBER_HESSIAN_LAGRANGIAN")
    H
  }

  # the compact limited memory BFGS approximation (src/quasi_newton.jl)
  qn <- if (hess_approx == "CompactLBFGS") .madnlp_lbfgs_new(n, o$quasi_newton_options) else NULL

  # solver state
  x <- numeric(N); zl <- numeric(N); zu <- numeric(N)
  y <- numeric(m); fgrad <- numeric(N); cval <- numeric(m); jacl <- numeric(N)
  obj_val <- 0
  J <- matrix(0, m, N); H <- matrix(0, N, N)
  reg <- rep(1, N); pr_diag <- rep(1, N); du_diag <- numeric(m)
  l_diag <- rep(1, nlb); u_diag <- rep(1, nub); l_lower <- numeric(nlb); u_lower <- numeric(nub)
  d <- numeric(nd); pv <- numeric(nd)
  x_trial <- numeric(N); c_trial <- numeric(m); obj_val_trial <- 0
  inf_pr <- 0; inf_du <- 0; inf_compl <- 0; inf_compl_mu <- 0
  theta_min <- 0; theta_max <- 0; mu <- 0; tau <- 0
  alpha <- 0; alpha_z <- 0; ftype <- " "
  del_w <- 0; del_w_last <- 0; del_c <- 0
  filter <- matrix(numeric(), 0L, 2L)
  RR <- NULL
  fac <- NULL

  x_lr <- function(v = x) v[ind_lb]
  x_ur <- function(v = x) v[ind_ub]
  reset_filter <- function(tm = theta_max) filter <<- matrix(c(tm, -Inf), 1L, 2L)
  filter_ok <- function(flt, th, ph) {
    if (is.na(th) || is.infinite(th) || is.na(ph) || is.infinite(ph)) return(FALSE)
    all(th <= flt[, 1L] | ph <= flt[, 2L])
  }
  augment <- function(flt, th, ph) rbind(flt, c((1 - o$gamma_theta) * th, ph - o$gamma_theta * th))

  # kernels (src/IPM/kernels.jl)
  varphi_of <- function(obj, xv, muv) {
    v <- obj
    dl <- xv[ind_lb] - xl[ind_lb]
    du <- xu[ind_ub] - xv[ind_ub]
    for (t in dl) v <- v + if (t < 0) Inf else -muv * log(t)
    for (t in du) v <- v + if (t < 0) Inf else -muv * log(t)
    v
  }
  varphi_d_of <- function(fv, xv, dx, muv) sum((fv - muv / (xv - xl) + muv / (xu - xv)) * dx)
  alpha_max_of <- function(xv, dx, tv) {
    a <- 1
    for (i in seq_len(N)) {
      if (dx[i] < 0) a <- min(a, (-xv[i] + xl[i]) * tv / dx[i])
      if (dx[i] > 0) a <- min(a, (-xv[i] + xu[i]) * tv / dx[i])
    }
    a
  }
  alpha_z_of <- function(zlr, zur, dzl, dzu, tv) {
    a <- 1
    for (i in seq_along(zlr)) if (dzl[i] < 0) a <- min(a, -zlr[i] * tv / dzl[i])
    for (i in seq_along(zur)) if (dzu[i] < 0) a <- min(a, -zur[i] * tv / dzu[i])
    a
  }
  get_sd <- function() max(s_max, (sum(abs(y)) + sum(abs(zl[ind_lb])) + sum(abs(zu[ind_ub]))) /
                             max(1, m + nlb + nub)) / s_max
  get_sc <- function() max(s_max, (sum(abs(zl[ind_lb])) + sum(abs(zu[ind_ub]))) / max(1, nlb + nub)) / s_max
  get_inf_du <- function(fv, sd) max(0, abs(fv - zl + zu + jacl)) / sd
  compl_of <- function(muv, sc) {
    v <- 0
    if (nlb) v <- max(v, abs((x[ind_lb] - xl[ind_lb]) * zl[ind_lb] - muv))
    if (nub) v <- max(v, abs((xu[ind_ub] - x[ind_ub]) * zu[ind_ub] - muv))
    v / sc
  }
  get_mu_next <- function(muv, mu_min, lin, sup) {
    a <- min(99 * mu_min / tol, 0.01)
    max(mu_min, a * tol, min(lin * muv, muv^sup))
  }
  get_tau <- function(muv) max(o$tau_min, 1 - muv)
  alpha_min_of <- function(theta, varphi_d) {
    if (varphi_d < 0) {
      if (theta <= theta_min) {
        o$alpha_min_frac * min(o$gamma_theta, o$gamma_phi * theta / (-varphi_d),
                               o$delta * theta^o$s_theta / (-varphi_d)^o$s_phi)
      } else {
        o$alpha_min_frac * min(o$gamma_theta, -o$gamma_phi * theta / varphi_d)
      }
    } else {
      o$alpha_min_frac * o$gamma_theta
    }
  }
  is_switching <- function(varphi_d, a, theta) varphi_d < 0 && a * (-varphi_d)^o$s_phi > o$delta * theta^o$s_theta
  is_armijo <- function(vt, v, a, varphi_d) vt <= v + o$eta_phi * a * varphi_d
  is_progress <- function(tt, th, vt, v) {
    (m > 0 && tt <= (1 - o$gamma_theta) * th + 10 * eps * abs(th)) ||
      vt <= v - o$gamma_phi * th + 10 * eps * abs(v)
  }
  rapid_increase <- function(v, vt) vt >= v && log10(vt - v) > o$obj_max_inc + max(1, log10(abs(v)))
  get_ftype <- function(flt, th, tt, v, vt, sw, arm) {
    if (!filter_ok(flt, tt, vt)) return(" ")
    if (rapid_increase(v, vt)) return(" ")
    if (th <= theta_min && sw) {
      if (arm) return("f")
    } else if (is_progress(tt, th, vt, v)) {
      return("h")
    }
    " "
  }
  reset_bound_dual <- function(z, gap, muv) {
    k <- o$kappa_sigma
    pmax(pmin(z, (k * muv) / gap), (muv / k) / gap)
  }
  adjust_boundary <- function(muv) {
    c1 <- eps * muv
    c2 <- eps^(3 / 4)
    if (nlb) {
      xr <- x[ind_lb]; lr <- xl[ind_lb]
      xl[ind_lb] <<- ifelse(xr - lr < c1, lr - c2 * pmax(1, abs(xr)), lr)
    }
    if (nub) {
      xr <- x[ind_ub]; ur <- xu[ind_ub]
      xu[ind_ub] <<- ifelse(ur - xr < c1, ur + c2 * pmax(1, abs(xr)), ur)
    }
  }
  avg_compl <- function() {
    if (nlb + nub == 0) return(0)
    (sum(x[ind_lb] * zl[ind_lb]) - sum(xl[ind_lb] * zl[ind_lb]) +
       sum(xu[ind_ub] * zu[ind_ub]) - sum(x[ind_ub] * zu[ind_ub])) / (nlb + nub)
  }
  min_compl <- function() {
    v <- Inf
    if (nlb) v <- min(v, (x[ind_lb] - xl[ind_lb]) * zl[ind_lb])
    if (nub) v <- min(v, (xu[ind_ub] - x[ind_ub]) * zu[ind_ub])
    v
  }

  # the augmented KKT system of SparseKKTSystem
  set_aug_diagonal <- function() {
    reg <<- rep(o$default_primal_regularization, N)
    du_diag <<- rep(-o$default_dual_regularization, m)
    l_diag <<- xl[ind_lb] - x[ind_lb]
    u_diag <<- x[ind_ub] - xu[ind_ub]
    l_lower <<- zl[ind_lb]
    u_lower <<- zu[ind_ub]
    pr_diag <<- reg
    pr_diag[ind_lb] <<- pr_diag[ind_lb] - l_lower / l_diag
    pr_diag[ind_ub] <<- pr_diag[ind_ub] - u_lower / u_diag
  }
  initialize_kkt <- function() {
    reg <<- rep(1, N); pr_diag <<- rep(1, N); du_diag <<- numeric(m)
    H <<- matrix(0, N, N)
    l_lower <<- numeric(nlb); u_lower <<- numeric(nub); l_diag <<- rep(1, nlb); u_diag <<- rep(1, nub)
    if (!is.null(qn)) qn$Bk <<- 0
  }
  # the Hessian block of the factorized matrix (the diagonal xi I of the
  # compact LBFGS) and of the KKT system (xi I - U U' + V V')
  hess_C <- function() {
    if (is.null(qn)) return(H)
    diag(c(rep(qn$Bk, n), numeric(ns)), N)
  }
  hess_matrix <- function() {
    if (is.null(qn)) return(H)
    B <- hess_C()
    k <- ncol(qn$U)
    if (k) {
      ii <- seq_len(n)
      B[ii, ii] <- B[ii, ii] - tcrossprod(qn$U) + tcrossprod(qn$V)
    }
    B
  }
  # C^-1 b for the factorized matrix; a singular matrix leaves b (MUMPS)
  csolve <- function(b) {
    if (is.null(fac)) factorize()
    if (fac$sing) b else tryCatch(.nlp_kkt_apply(fac$lu, b), error = function(e) b)
  }
  # factorize_wrapper!: the matrix, its inertia (as MUMPS reports it) and
  # the factors
  factorize <- function() {
    dim <- N + m
    K <- matrix(0, dim, dim)
    K[iP, iP] <- hess_C()
    diag(K)[iP] <- diag(K)[iP] + pr_diag
    if (m) {
      K[iY, iP] <- J
      K[iP, iY] <- t(J)
      diag(K)[iY] <- du_diag
    }
    tr <- .ipm_triplets(K)
    inr <- .nlp_inertia(tr$i, tr$j, tr$v, dim)
    sing <- inr$zero > 0
    lu <- if (sing) NULL else tryCatch(.nlp_kkt_factor(tr$i, tr$j, tr$v, dim, rep(0, dim)),
                                       error = function(e) NULL)
    sing <- is.null(lu)
    neg <- inr$neg
    fac <<- list(pos = dim - sing - neg, zero = as.integer(sing), neg = neg, sing = sing, lu = lu)
  }
  # solve_kkt! of an AbstractReducedKKTSystem on an unreduced vector
  solve_kkt <- function(w) {
    if (nlb) w[ind_lb] <- w[ind_lb] - w[iL] / l_diag
    if (nub) w[ind_ub] <- w[ind_ub] - w[iU] / u_diag
    pd <- c(iP, iY)
    w[pd] <- csolve(w[pd])
    k <- if (is.null(qn)) 0L else ncol(qn$U)
    if (k) {
      # the Sherman-Morrison-Woodbury correction of the low rank part
      E <- matrix(0, N + m, 2L * k)
      E[seq_len(n), ] <- cbind(qn$U, qn$V)
      Hh <- apply(E, 2L, csolve)
      Tk <- diag(c(rep(-1, k), rep(1, k)), 2L * k) + crossprod(E, Hh)
      xr <- solve(Tk, as.numeric(crossprod(E, w[pd])))
      w[pd] <- w[pd] - as.numeric(Hh %*% xr)
    }
    if (nlb) w[iL] <- (-w[iL] + l_lower * w[ind_lb]) / l_diag
    if (nub) w[iU] <- (w[iU] - u_lower * w[ind_ub]) / u_diag
    w
  }
  # mul! of the unreduced KKT system
  kkt_mul <- function(v) {
    vp <- v[iP]; vy <- v[iY]
    wp <- as.numeric(hess_matrix() %*% vp) + reg * vp
    if (m) wp <- wp + as.numeric(crossprod(J, vy))
    wy <- if (m) as.numeric(J %*% vp) + du_diag * vy else numeric()
    if (nlb) wp[ind_lb] <- wp[ind_lb] - v[iL]
    if (nub) wp[ind_ub] <- wp[ind_ub] + v[iU]
    c(wp, wy, vp[ind_lb] * l_lower - v[iL] * l_diag, vp[ind_ub] * u_lower + v[iU] * u_diag)
  }
  # solve_refine! of RichardsonIterator: the solution and whether the
  # residual ratio is acceptable
  solve_refine <- function(b) {
    nb <- max(abs(b), 0)
    xs <- numeric(nd)
    ratio <- 0
    if (nb != 0) {
      w <- b
      cnt$ir <<- 0L
      repeat {
        w <- solve_kkt(w)
        xs <- xs + w
        w <- b - kkt_mul(xs)
        nw <- max(abs(w), 0); nx <- max(abs(xs), 0)
        ratio <- nw / (min(nx, 1e6 * nb) + nb)
        cnt$ir <<- cnt$ir + 1L
        if (cnt$ir >= rich_max || ratio < rich_tol) break
      }
    }
    list(x = xs, ok = isTRUE(ratio < rich_acc))
  }
  set_aug_rhs <- function(cv, muv) {
    w <- numeric(nd)
    w[iP] <- -fgrad + zl - zu - jacl
    w[iY] <- -cv
    if (nlb) w[iL] <- (xl[ind_lb] - x[ind_lb]) * zl[ind_lb] + muv
    if (nub) w[iU] <- (xu[ind_ub] - x[ind_ub]) * zu[ind_ub] - muv
    w
  }
  dual_inf_perturbation <- function(w, muv) {
    w[ind_llb] <- w[ind_llb] - muv * o$kappa_d
    w[ind_uub] <- w[ind_uub] + muv * o$kappa_d
    w
  }
  regularize <- function(dw, dc) {
    reg <<- reg + dw
    pr_diag <<- pr_diag + dw
    du_diag <<- du_diag - dc
  }
  inertia_correct <- function() fac$zero == 0 && fac$pos == N
  regularize_dual <- function() condensed || fac$zero != 0
  next_del_w <- function(n_trial) {
    if (n_trial == 0L) {
      del_w <<- if (del_w_last == 0) o$first_hessian_perturbation else
        max(o$min_hessian_perturbation, o$perturb_dec_fact * del_w_last)
    } else {
      del_w <<- del_w * (if (del_w_last == 0) o$perturb_inc_fact_first else o$perturb_inc_fact)
      if (del_w > o$max_hessian_perturbation) {
        cnt$k <<- cnt$k + 1L
        return(FALSE)
      }
    }
    TRUE
  }
  jac_reg <- function() o$jacobian_regularization_value * mu^o$jacobian_regularization_exponent
  # inertia_correction! (InertiaBased, InertiaIgnore, InertiaFree): the
  # step in d from the right-hand side pv
  inertia_correction <- function() {
    n_trial <- 0L
    dwp <- 0; dcp <- 0
    del_w <<- 0; del_c <<- 0
    if (ic_method == "InertiaFree") {
      g <- fgrad - mu / (x - xl) + mu / (xu - x) + jacl
      p0 <- numeric(nd)
      p0[iY] <- -cval
      attempt <- function() {
        factorize()
        r0 <- solve_refine(p0)
        if (!r0$ok) return(FALSE)
        r1 <- solve_refine(pv)
        d <<- r1$x
        if (!r1$ok) return(FALSE)
        nvec <- r0$x[iP]
        tt <- d[iP] - nvec
        wx <- as.numeric(hess_C() %*% tt) + tt * pr_diag
        sum(wx * tt) + max(sum(wx * nvec) - sum(g * nvec), 0) - o$inertia_free_tol * sum(tt * tt) >= 0
      }
      ok <- attempt()
      while (!ok) {
        if (!next_del_w(n_trial)) return(FALSE)
        del_c <<- jac_reg()
        regularize(del_w - dwp, del_c - dcp)
        dwp <- del_w; dcp <- del_c
        ok <- attempt()
        n_trial <- n_trial + 1L
      }
    } else {
      attempt <- function() {
        factorize()
        if (ic_method == "InertiaBased" && !inertia_correct()) return(FALSE)
        r <- solve_refine(pv)
        d <<- r$x
        r$ok
      }
      ok <- attempt()
      while (!ok) {
        if (!next_del_w(n_trial)) return(FALSE)
        del_c <<- if (ic_method == "InertiaIgnore" || regularize_dual()) jac_reg() else 0
        regularize(del_w - dwp, del_c - dcp)
        dwp <- del_w; dcp <- del_c
        ok <- attempt()
        n_trial <- n_trial + 1L
      }
    }
    if (del_w != 0) del_w_last <<- del_w
    TRUE
  }

  print_iter <- function(resto = FALSE) {
    if (!verbose) return(invisible())
    if (cnt$k %% 10L == 0L)
      cat("iter    objective    inf_pr   inf_du inf_compl lg(mu) lg(rg) alpha_pr ir ls\n")
    if (resto) {
      ip <- RR$inf_pr_R; id <- RR$inf_du_R; ic <- RR$inf_compl_R; lm <- log10(RR$mu_R)
    } else {
      ip <- inf_pr; id <- inf_du; ic <- inf_compl; lm <- log10(mu)
    }
    cat(sprintf("%4d%s% 10.7e %6.2e %6.2e %7.2e %5.1f  %s  %6.2e %2d %2d%s\n", cnt$k, if (resto) "r" else " ",
                obj_val / obj_scale, ip, id, ic, lm,
                if (del_w == 0) "   - " else sprintf("%5.1f", log10(del_w)), alpha, cnt$ir, cnt$l, ftype))
  }

  # update_barrier! (src/IPM/barrier.jl)
  update_monotone <- function(sc) {
    inf_compl_mu <<- compl_of(mu, sc)
    while (mu > max(barrier$mu_min, tol / 10) &&
           max(inf_pr, inf_du, inf_compl_mu) <= o$barrier_tol_factor * mu) {
      mu <<- get_mu_next(mu, barrier$mu_min, barrier$mu_linear_decrease_factor,
                         barrier$mu_superlinear_decrease_power)
      tau <<- get_tau(mu)
      inf_compl_mu <<- compl_of(mu, sc)
      reset_filter()
    }
  }
  quality <- function(sigma, aff, cen, res_primal, res_dual) {
    dd <- aff + sigma * cen
    d <<- dd
    a_pr <- alpha_max_of(x, dd[iP], tau)
    a_du <- alpha_z_of(zl[ind_lb], zu[ind_ub], dd[iL], dd[iU], tau)
    cl <- sum(((x[ind_lb] + a_pr * dd[ind_lb] - xl[ind_lb]) * (zl[ind_lb] + a_du * dd[iL]))^2)
    cu <- sum(((xu[ind_ub] - x[ind_ub] - a_pr * dd[ind_ub]) * (zu[ind_ub] + a_du * dd[iU]))^2)
    ipr <- if (m > 0) (1 - a_pr)^2 * res_primal^2 / m else 0
    idu <- (1 - a_du)^2 * res_dual^2 / N
    idu + ipr + (cl + cu) / (nlb + nub)
  }
  golden <- function(s_lb, s_ub, aff, cen, rp, rd) {
    gfac <- (3 - sqrt(5)) / 2
    s1 <- s_lb; s2 <- s_ub
    q <- function(s) quality(s, aff, cen, rp, rd)
    p1 <- q(s1); p2 <- q(s2)
    s1_in <- s1; s2_in <- s2; p1_in <- p1; p2_in <- p2
    sm1 <- s_lb + gfac * (s_ub - s_lb)
    sm2 <- s_lb + (1 - gfac) * (s_ub - s_lb)
    pm1 <- q(sm1); pm2 <- q(sm2)
    for (i in seq_len(barrier$max_gs_iter)) {
      if (pm1 > pm2) {
        s1 <- sm1; p1 <- pm1
        sm1 <- sm2
        sm2 <- s1 + (1 - gfac) * (s2 - s1)
        pm1 <- pm2
        pm2 <- q(sm2)
      } else {
        s2 <- sm2; p2 <- pm2
        sm2 <- sm1
        sm1 <- s1 + gfac * (s2 - s1)
        pm2 <- pm1
        pm1 <- q(sm1)
      }
      if (s2 - s1 < barrier$sigma_tol * s2) break
    }
    if (pm1 < pm2) { sigma <- sm1; phi <- pm1 } else { sigma <- sm2; phi <- pm2 }
    if (s2 == s2_in && p2_in < phi) {
      sigma <- s2_in
    } else if (s1 == s1_in && p1_in < phi) {
      sigma <- s1_in
    }
    sigma
  }
  adaptive_mu <- function() {
    if (nlb + nub == 0) return(barrier$mu_min)
    if (barrier$type == "LOQOUpdate") {
      mua <- avg_compl()
      xi <- min_compl() / mua
      sigma <- barrier$gamma * min((1 - barrier$r) * ((1 - xi) / xi), 2)^3
      return(min(max(sigma * mua, barrier$mu_min), barrier$mu_max))
    }
    pv <<- set_aug_rhs(cval, 0)
    res_primal <- sqrt(sum(pv[iY]^2))
    res_dual <- sqrt(sum(pv[iP]^2))
    aff <- solve_kkt(pv)
    mua <- avg_compl()
    pc <- numeric(nd)
    if (nlb) pc[iL] <- mua
    if (nub) pc[iU] <- -mua
    pc <- dual_inf_perturbation(pc, mua)
    pv <<- pc
    cen <- solve_kkt(pc)
    phi1 <- quality(1, aff, cen, res_primal, res_dual)
    s1m <- 1 - 1e-4
    phi1m <- quality(s1m, aff, cen, res_primal, res_dual)
    if (phi1m > phi1) {
      smin <- 1
      smax <- min(barrier$sigma_max, barrier$mu_max / mua)
    } else {
      smin <- max(barrier$sigma_min, barrier$mu_min / mua)
      smax <- min(max(smin, s1m), barrier$mu_max / mua)
    }
    sopt <- golden(smin, smax, aff, cen, res_primal, res_dual)
    barrier$n_update <<- barrier$n_update + 1L
    min(max(sopt * mua, barrier$mu_min), barrier$mu_max)
  }
  update_barrier <- function(sc) {
    if (!adaptive) return(update_monotone(sc))
    old_mu <- mu
    progress <- TRUE
    if (isTRUE(barrier$globalization)) {
      theta <- sum(abs(cval))
      varphi <- varphi_of(obj_val, x, mu)
      kkt_error <- max(inf_pr, inf_du, inf_compl)
      dl <- 1e-5 * min(1, kkt_error)
      progress <- filter_ok(filter, theta + dl, varphi + dl)
    }
    if (!barrier$free_mode) {
      if (progress) barrier$free_mode <<- TRUE else update_monotone(sc)
    } else if (!progress) {
      barrier$free_mode <<- FALSE
      mu <<- min(max(0.8 * avg_compl(), barrier$mu_min), barrier$mu_max)
    }
    if (barrier$free_mode) mu <<- adaptive_mu()
    if (mu != old_mu) {
      tau <<- get_tau(mu)
      reset_filter()
    }
  }

  # eval_lag_hess_wrapper!: the exact Hessian or the update of the
  # quasi-Newton approximation
  eval_lag_hess <- function(resto = FALSE) {
    if (is.null(qn)) {
      H <<- eval_hess(x, y, if (resto) 0 else 1)
      return(invisible())
    }
    if (cnt$grad >= 2L) {
      sk <- x[seq_len(n)] - qn$last_x
      yk <- fgrad[seq_len(n)] - qn$last_g
      if (m > 0) {
        jacl <<- as.numeric(crossprod(J, y))
        yk <- yk + jacl[seq_len(n)]
        Jold <- pb$jac(out_x(c(qn$last_x, numeric(ns))))[, free, drop = FALSE] * con_scale
        yk <- yk - as.numeric(crossprod(Jold, y))
      }
      qn <<- .madnlp_lbfgs_update(qn, sk, yk)
    } else {
      qn <<- .madnlp_lbfgs_init(qn, fgrad[seq_len(n)], obj_val)
    }
    qn$last_x <<- x[seq_len(n)]
    qn$last_g <<- fgrad[seq_len(n)]
  }

  # initialize! (src/IPM/solver.jl and the callback's initialize!)
  initialize <- function() {
    lv <- lvar; uv <- uvar
    relax <- function(lo, hi) {
      if (o$bound_relax_factor > 0) {
        lo <- lo - pmax(1, abs(lo)) * o$bound_relax_factor
        hi <- hi + pmax(1, abs(hi)) * o$bound_relax_factor
      }
      list(lo = lo, hi = hi)
    }
    push <- function(v, lo, hi) {
      bp <- o$bound_push; bf <- o$bound_fac
      for (i in seq_along(v)) {
        if (lo[i] != -Inf && hi[i] != Inf) {
          v[i] <- min(hi[i] - min(bp * max(1, abs(hi[i])), bf * (hi[i] - lo[i])),
                      max(lo[i] + min(bp * max(1, abs(lo[i])), bf * (hi[i] - lo[i])), v[i]))
        } else if (lo[i] != -Inf) {
          v[i] <- max(lo[i] + bp * max(1, abs(lo[i])), v[i])
        } else if (hi[i] != Inf) {
          v[i] <- min(hi[i] - bp * max(1, abs(hi[i])), v[i])
        }
      }
      v
    }
    rb <- relax(lv, uv)
    xv <- push(x_full[free], rb$lo, rb$hi)
    x[seq_len(n)] <<- xv
    cv <- pb$cons(out_x(c(xv, numeric(ns))))
    rhs <<- if (eq_treatment == "EnforceEquality") ifelse(lcon == ucon, lcon, 0) else numeric(m)
    rs <- relax(lcon[ind_ineq], ucon[ind_ineq])
    sv <- push(cv[ind_ineq], rs$lo, rs$hi)
    x[slk] <<- sv
    xl <<- c(rb$lo, rs$lo); xu <<- c(rb$hi, rs$hi)
    jacl <<- numeric(N)
    zl[ind_lb] <<- 1; zu[ind_ub] <<- 1
    if (isTRUE(o$nlp_scaling)) {
      Jx <- if (m) pb$jac(out_x(x))[, free, drop = FALSE] else matrix(0, 0, n)
      cs <- rep(1, m)
      if (m) cs <- pmax(1, apply(abs(Jx), 1L, max))
      con_scale <<- pmin(1, o$nlp_scaling_max_gradient / cs)
      g0 <- pb$grad(out_x(x))[free]
      obj_scale <<- min(1, o$nlp_scaling_max_gradient / max(abs(g0), 0))
      y <<- y / con_scale
      rhs <<- rhs * con_scale
      x[slk] <<- x[slk] * con_scale[ind_ineq]
      xl[slk] <<- xl[slk] * con_scale[ind_ineq]
      xu[slk] <<- xu[slk] * con_scale[ind_ineq]
    }
    initialize_kkt()
    J <<- eval_jac(x)
    fgrad <<- eval_grad(x)
    if (!isTRUE(o$dual_initialized)) {
      if (dual_init == "DualInitializeLeastSquares") {
        pv <<- numeric(nd)
        pv[iP] <<- -fgrad + zl - zu
        factorize()
        r <- solve_refine(pv)
        d <<- r$x
        y <<- if (!r$ok || max(abs(d[iY]), 0) > o$constr_mult_init_max) numeric(m) else d[iY]
      } else {
        y <<- numeric(m)
      }
    }
    obj_val <<- eval_f(x)
    cval <<- eval_cons(x)
    eval_lag_hess()
    theta <- sum(abs(cval))
    theta_max <<- 1e4 * max(1, theta)
    theta_min <<- 1e-4 * max(1, theta)
    mu <<- barrier$mu_init
    tau <<- max(o$tau_min, 1 - barrier$mu_init)
    reset_filter()
    "REGULAR"
  }

  # second_order_correction
  soc <- function(alpha_max, theta, varphi, theta_trial, varphi_d, sw) {
    wy <- c_trial + alpha_max * cval
    theta_soc_old <- theta_trial
    for (p in seq_len(o$max_soc)) {
      pv <<- dual_inf_perturbation(set_aug_rhs(wy, mu), mu)
      r <- solve_refine(pv)
      wx <- r$x[iP]
      alpha_soc <- alpha_max_of(x, wx, tau)
      x_trial <<- x + alpha_soc * wx
      c_trial <<- eval_cons(x_trial)
      obj_val_trial <<- eval_f(x_trial)
      theta_soc <- sum(abs(c_trial))
      varphi_soc <- varphi_of(obj_val_trial, x_trial, mu)
      if (!filter_ok(filter, theta_soc, varphi_soc)) break
      if (theta <= theta_min && sw) {
        if (is_armijo(varphi_soc, varphi, alpha, varphi_d)) {
          ftype <<- "F"
          alpha <<- alpha_soc
          return(TRUE)
        }
      } else if (is_progress(theta_soc, theta, varphi_soc, varphi)) {
        ftype <<- "H"
        alpha <<- alpha_soc
        return(TRUE)
      }
      if (theta_soc > o$kappa_soc * theta_soc_old) break
      theta_soc_old <- theta_soc
    }
    FALSE
  }
  reset_to_regular <- function() {
    y <<- numeric(m)
    zl[ind_lb] <<- 1; zu[ind_ub] <<- 1
    reset_filter()
  }
  # filter_line_search!
  filter_line_search <- function() {
    theta <- sum(abs(cval))
    varphi <- varphi_of(obj_val, x, mu)
    dx <- d[iP]
    varphi_d <- varphi_d_of(fgrad, x, dx, mu)
    alpha_max <- alpha_max_of(x, dx, tau)
    alpha_z <<- alpha_z_of(zl[ind_lb], zu[ind_ub], d[iL], d[iU], tau)
    alpha_min <- alpha_min_of(theta, varphi_d)
    cnt$l <<- 1L
    alpha <<- alpha_max
    varphi_trial <- 0; theta_trial <- 0
    small <- max(abs(dx) / (1 + abs(x)), 0) < 10 * eps
    sw <- is_switching(varphi_d, alpha, theta)
    arm <- FALSE
    unsuccessful <- FALSE
    repeat {
      x_trial <<- x
      x_trial <<- x + alpha * dx
      obj_val_trial <<- eval_f(x_trial)
      c_trial <<- eval_cons(x_trial)
      theta_trial <- sum(abs(c_trial))
      varphi_trial <- varphi_of(obj_val_trial, x_trial, mu)
      arm <- is_armijo(varphi_trial, varphi, alpha, varphi_d)
      if (small) break
      ftype <<- get_ftype(filter, theta, theta_trial, varphi, varphi_trial, sw, arm)
      if (ftype %in% c("f", "h")) break
      if (cnt$l == 1L && theta_trial >= theta &&
          soc(alpha_max, theta, varphi, theta_trial, varphi_d, sw)) break
      unsuccessful <- TRUE
      alpha <<- alpha / 2
      cnt$l <<- cnt$l + 1L
      if (alpha < alpha_min) {
        cnt$k <<- cnt$k + 1L
        return("RESTORE")
      } else if (alpha * sqrt(sum(dx^2)) < eps * 10) {
        cnt$restoration_fail_count <<- cnt$restoration_fail_count + 1L
        if (cnt$restoration_fail_count >= 4L)
          return(if (cnt$acceptable_cnt > 0) "SOLVED_TO_ACCEPTABLE_LEVEL" else "SEARCH_DIRECTION_BECOMES_TOO_SMALL")
        reset_to_regular()
        cnt$k <<- cnt$k + 1L
        return("REGULAR")
      }
    }
    if (unsuccessful) {
      cnt$unsuccessful_iterate <<- cnt$unsuccessful_iterate + 1L
      if (cnt$unsuccessful_iterate >= 4L) {
        if (theta_max / 10 > theta_trial) {
          theta_max <<- theta_max / 10
          reset_filter()
        }
        cnt$unsuccessful_iterate <<- 0L
      }
    } else {
      cnt$unsuccessful_iterate <<- 0L
    }
    if (!sw || !arm) filter <<- augment(filter, theta_trial, varphi_trial)
    "LINESEARCH_SUCCEEDED"
  }

  wall <- function() proc.time()[["elapsed"]] - start_time >= o$max_wall_time
  regular <- function() {
    repeat {
      if (cnt$k != 0L && !isTRUE(o$jacobian_constant)) J <<- eval_jac(x)
      jacl <<- if (m) as.numeric(crossprod(J, y)) else numeric(N)
      sd <- get_sd(); sc <- get_sc()
      inf_pr <<- max(abs(cval), 0)
      inf_du <<- get_inf_du(fgrad, sd)
      inf_compl <<- compl_of(0, sc)
      print_iter()
      inf_total <- max(inf_pr, inf_du, inf_compl)
      if (inf_total <= tol) return("SOLVE_SUCCEEDED")
      if (inf_total <= o$acceptable_tol) {
        if (cnt$acceptable_cnt < o$acceptable_iter) cnt$acceptable_cnt <<- cnt$acceptable_cnt + 1L else
          return("SOLVED_TO_ACCEPTABLE_LEVEL")
      } else {
        cnt$acceptable_cnt <<- 0L
      }
      if (inf_total >= o$diverging_iterates_tol) return("DIVERGING_ITERATES")
      if (cnt$k >= o$max_iter) return("MAXIMUM_ITERATIONS_EXCEEDED")
      if (wall()) return("MAXIMUM_WALLTIME_EXCEEDED")
      if (cnt$k != 0L && !isTRUE(o$hessian_constant)) eval_lag_hess()
      update_barrier(sc)
      set_aug_diagonal()
      pv <<- dual_inf_perturbation(set_aug_rhs(cval, mu), mu)
      if (!inertia_correction()) return("ROBUST")
      st <- filter_line_search()
      if (st != "LINESEARCH_SUCCEEDED") return(st)
      x <<- x_trial
      cval <<- c_trial
      obj_val <<- obj_val_trial
      adjust_boundary(mu)
      y <<- y + alpha * d[iY]
      zl[ind_lb] <<- zl[ind_lb] + alpha_z * d[iL]
      zu[ind_ub] <<- zu[ind_ub] + alpha_z * d[iU]
      zl[ind_lb] <<- reset_bound_dual(zl[ind_lb], x[ind_lb] - xl[ind_lb], mu)
      zu[ind_ub] <<- reset_bound_dual(zu[ind_ub], xu[ind_ub] - x[ind_ub], mu)
      fgrad <<- eval_grad(x)
      cnt$k <<- cnt$k + 1L
    }
  }
  get_F <- function(muv) {
    F1 <- sum(abs(cval))
    F2 <- sum(abs(fgrad - zl + zu + jacl))
    F3 <- 0
    for (i in ind_lb) F3 <- F3 + if (x[i] >= xl[i] && zl[i] >= 0) abs((x[i] - xl[i]) * zl[i] - muv) else Inf
    F4 <- 0
    for (i in ind_ub) F4 <- F4 + if (xu[i] >= x[i] && zu[i] >= 0) abs((xu[i] - x[i]) * zu[i] - muv) else Inf
    F1 + F2 + F3 + F4
  }
  restore <- function() {
    del_w <<- 0
    x_bak <- x; y_bak <- y; c_bak <- cval
    Fv <- get_F(mu)
    alpha_z <<- 0
    ftype <<- "R"
    repeat {
      a_max <- alpha_max_of(x, d[iP], tau)
      alpha <<- min(a_max, alpha_z_of(zl[ind_lb], zu[ind_ub], d[iL], d[iU], tau))
      x <<- x + alpha * d[iP]
      y <<- y + alpha * d[iY]
      zl[ind_lb] <<- zl[ind_lb] + alpha * d[iL]
      zu[ind_ub] <<- zu[ind_ub] + alpha * d[iU]
      cval <<- eval_cons(x)
      fgrad <<- eval_grad(x)
      obj_val <<- eval_f(x)
      if (!isTRUE(o$jacobian_constant)) J <<- eval_jac(x)
      jacl <<- if (m) as.numeric(crossprod(J, y)) else numeric(N)
      F_trial <- get_F(mu)
      if (F_trial > o$soft_resto_pderror_reduction_factor * Fv) {
        x <<- x_bak
        y <<- y_bak
        cval <<- c_bak
        return("ROBUST")
      }
      adjust_boundary(mu)
      Fv <- F_trial
      theta <- sum(abs(cval))
      varphi <- varphi_of(obj_val, x, mu)
      cnt$k <<- cnt$k + 1L
      if (filter_ok(filter, theta, varphi)) return("REGULAR") else cnt$t <<- cnt$t + 1L
      if (cnt$k >= o$max_iter) return("MAXIMUM_ITERATIONS_EXCEEDED")
      if (wall()) return("MAXIMUM_WALLTIME_EXCEEDED")
      sd <- get_sd(); sc <- get_sc()
      inf_pr <<- max(abs(cval), 0)
      inf_du <<- get_inf_du(fgrad, sd)
      inf_compl <<- compl_of(0, sc)
      inf_compl_mu <<- compl_of(mu, sc)
      print_iter()
      if (!isTRUE(o$hessian_constant)) eval_lag_hess()
      set_aug_diagonal()
      pv <<- dual_inf_perturbation(set_aug_rhs(cval, mu), mu)
      factorize()
      d <<- solve_refine(pv)$x
      ftype <<- "f"
    }
  }

  # the robust restoration phase (RobustRestorer)
  rho <- o$rho
  init_rr <- function() {
    if (is.null(RR)) RR <<- list()
    RR$x_ref <<- x
    RR$theta_ref <<- sum(abs(cval))
    RR$D_R <<- pmin(1, 1 / abs(RR$x_ref))
    RR$mu_R <<- max(mu, max(abs(cval), 0))
    RR$tau_R <<- max(o$tau_min, 1 - RR$mu_R)
    RR$zeta <<- sqrt(RR$mu_R)
    mR <- RR$mu_R
    a <- (mR - rho * cval) / (2 * rho)
    RR$nn <<- a + sqrt(a^2 + mR * cval / (2 * rho))
    RR$pp <<- cval + RR$nn
    RR$zp <<- RR$mu_R / RR$pp
    RR$zn <<- RR$mu_R / RR$nn
    RR$obj_val_R <<- obj_R(RR$pp, RR$nn, x)
    RR$f_R <<- numeric(N)
    RR$filter <<- matrix(c(theta_max, -Inf), 1L, 2L)
    y <<- numeric(m)
    zl[ind_lb] <<- pmin(rho, zl[ind_lb])
    zu[ind_ub] <<- pmin(rho, zu[ind_ub])
    cnt$t <<- 0L
    del_w <<- 0
  }
  obj_R <- function(pp, nn, xv) rho * sum(pp + nn) + sum(RR$zeta / 2 * RR$D_R^2 * (xv - RR$x_ref)^2)
  compl_R <- function(muv, sc) {
    v <- compl_of(muv, 1)
    if (m) v <- max(v, abs(RR$pp * RR$zp - muv), abs(RR$nn * RR$zn - muv))
    v / sc
  }
  varphi_R_of <- function(objR, xv, pp, nn) {
    v <- varphi_of(objR, xv, RR$mu_R)
    for (t in pp) v <- v - if (t < 0) Inf else RR$mu_R * log(t)
    for (t in nn) v <- v - if (t < 0) Inf else RR$mu_R * log(t)
    v
  }
  update_monotone_RR <- function(sc) {
    icm <- compl_R(RR$mu_R, sc)
    while (RR$mu_R >= barrier$mu_min && max(RR$inf_pr_R, RR$inf_du_R, icm) <= o$barrier_tol_factor * RR$mu_R) {
      RR$mu_R <<- get_mu_next(RR$mu_R, barrier$mu_min, barrier$mu_linear_decrease_factor,
                              barrier$mu_superlinear_decrease_power)
      icm <- compl_R(RR$mu_R, sc)
      RR$tau_R <<- max(o$tau_min, 1 - RR$mu_R)
      RR$zeta <<- sqrt(RR$mu_R)
      RR$filter <<- matrix(c(theta_max, -Inf), 1L, 2L)
    }
  }
  line_search_RR <- function() {
    theta_R <- sum(abs(cval - RR$pp + RR$nn))
    varphi_R <- varphi_R_of(RR$obj_val_R, x, RR$pp, RR$nn)
    dx <- d[iP]
    mR <- RR$mu_R
    varphi_d_R <- sum((RR$f_R - mR / (x - xl) + mR / (xu - x)) * dx) + sum((rho - mR / RR$pp) * RR$dpp) +
      sum((rho - mR / RR$nn) * RR$dnn)
    a_max <- alpha_max_of(x, dx, RR$tau_R)
    for (i in seq_len(m)) {
      if (RR$dpp[i] < 0) a_max <- min(a_max, -RR$pp[i] * RR$tau_R / RR$dpp[i])
      if (RR$dnn[i] < 0) a_max <- min(a_max, -RR$nn[i] * RR$tau_R / RR$dnn[i])
    }
    az <- alpha_z_of(zl[ind_lb], zu[ind_ub], d[iL], d[iU], RR$tau_R)
    for (i in seq_len(m)) {
      if (RR$dzp[i] < 0) az <- min(az, -RR$zp[i] * RR$tau_R / RR$dzp[i])
      if (RR$dzn[i] < 0) az <- min(az, -RR$zn[i] * RR$tau_R / RR$dzn[i])
    }
    alpha_z <<- az
    alpha_min <- alpha_min_of(theta_R, varphi_d_R)
    alpha <<- a_max
    cnt$l <<- 1L
    theta_R_trial <- 0; varphi_R_trial <- 0
    small <- max(abs(dx) / (1 + abs(x)), 0) < 10 * eps
    sw <- is_switching(varphi_d_R, alpha, theta_R)
    arm <- FALSE
    repeat {
      x_trial <<- x + alpha * dx
      RR$pp_trial <<- RR$pp + alpha * RR$dpp
      RR$nn_trial <<- RR$nn + alpha * RR$dnn
      RR$obj_val_R_trial <<- obj_R(RR$pp_trial, RR$nn_trial, x_trial)
      c_trial <<- eval_cons(x_trial)
      theta_R_trial <- sum(abs(c_trial - RR$pp_trial + RR$nn_trial))
      varphi_R_trial <- varphi_R_of(RR$obj_val_R_trial, x_trial, RR$pp_trial, RR$nn_trial)
      arm <- is_armijo(varphi_R_trial, varphi_R, alpha, varphi_d_R)
      if (small) break
      ftype <<- get_ftype(RR$filter, theta_R, theta_R_trial, varphi_R, varphi_R_trial, sw, arm)
      if (ftype %in% c("f", "h")) break
      alpha <<- alpha / 2
      cnt$l <<- cnt$l + 1L
      if (alpha < alpha_min) {
        cnt$restoration_fail_count <<- cnt$restoration_fail_count + 1L
        if (cnt$restoration_fail_count >= 4L) return("RESTORATION_FAILED")
        reset_to_regular()
        cnt$k <<- cnt$k + 1L
        cnt$t <<- cnt$t + 1L
        return("REGULAR")
      } else if (alpha < eps * 10) {
        return(if (cnt$acceptable_cnt > 0) "SOLVED_TO_ACCEPTABLE_LEVEL" else "SEARCH_DIRECTION_BECOMES_TOO_SMALL")
      }
    }
    if (!sw || !arm) RR$filter <<- augment(RR$filter, theta_R_trial, varphi_R_trial)
    "LINESEARCH_SUCCEEDED"
  }
  robust <- function() {
    init_rr()
    repeat {
      if (!isTRUE(o$jacobian_constant)) J <<- eval_jac(x)
      jacl <<- if (m) as.numeric(crossprod(J, y)) else numeric(N)
      sd <- get_sd(); sc <- get_sc()
      inf_pr <<- max(abs(cval), 0)
      inf_du <<- get_inf_du(fgrad, sd)
      inf_compl <<- compl_of(0, sc)
      RR$inf_pr_R <<- max(abs(cval - RR$pp + RR$nn), 0)
      RR$inf_du_R <<- max(abs(RR$f_R - zl + zu + jacl), abs(rho - y - RR$zp), abs(rho + y - RR$zn), 0) / sd
      RR$inf_compl_R <<- compl_R(0, sc)
      print_iter(resto = TRUE)
      if (max(RR$inf_pr_R, RR$inf_du_R, RR$inf_compl_R) <= tol) return("INFEASIBLE_PROBLEM_DETECTED")
      if (cnt$k >= o$max_iter) return("MAXIMUM_ITERATIONS_EXCEEDED")
      if (wall()) return("MAXIMUM_WALLTIME_EXCEEDED")
      update_monotone_RR(sc)
      if (!isTRUE(o$hessian_constant)) eval_lag_hess(resto = TRUE)
      # set_aug_RR! and set_aug_rhs_RR!
      reg <<- o$default_primal_regularization + RR$zeta * RR$D_R^2
      du_diag <<- -o$default_dual_regularization - RR$pp / RR$zp - RR$nn / RR$zn
      l_lower <<- zl[ind_lb]; u_lower <<- zu[ind_ub]
      l_diag <<- xl[ind_lb] - x[ind_lb]; u_diag <<- x[ind_ub] - xu[ind_ub]
      pr_diag <<- reg
      pr_diag[ind_lb] <<- pr_diag[ind_lb] - l_lower / l_diag
      pr_diag[ind_ub] <<- pr_diag[ind_ub] - u_lower / u_diag
      mR <- RR$mu_R
      w <- numeric(nd)
      w[iP] <- -RR$f_R + zl - zu - jacl
      w[iY] <- -cval + RR$pp - RR$nn + (mR - (rho - y) * RR$pp) / RR$zp - (mR - (rho + y) * RR$nn) / RR$zn
      if (nlb) w[iL] <- (xl[ind_lb] - x[ind_lb]) * zl[ind_lb] + mR
      if (nub) w[iU] <- (xu[ind_ub] - x[ind_ub]) * zu[ind_ub] - mR
      pv <<- w
      if (!inertia_correction()) return("RESTORATION_FAILED")
      dl <- d[iY]
      RR$dzp <<- rho - y - dl - RR$zp
      RR$dzn <<- rho + y + dl - RR$zn
      RR$dpp <<- -RR$pp + mR / RR$zp - (RR$pp / RR$zp) * RR$dzp
      RR$dnn <<- -RR$nn + mR / RR$zn - (RR$nn / RR$zn) * RR$dzn
      st <- line_search_RR()
      if (st != "LINESEARCH_SUCCEEDED") return(st)
      x <<- x_trial
      cval <<- c_trial
      RR$pp <<- RR$pp_trial
      RR$nn <<- RR$nn_trial
      RR$obj_val_R <<- RR$obj_val_R_trial
      RR$f_R <<- RR$zeta * RR$D_R^2 * (x - RR$x_ref)
      y <<- y + alpha * d[iY]
      RR$zp <<- RR$zp + alpha_z * RR$dzp
      RR$zn <<- RR$zn + alpha_z * RR$dzn
      zl[ind_lb] <<- zl[ind_lb] + alpha_z * d[iL]
      zu[ind_ub] <<- zu[ind_ub] + alpha_z * d[iU]
      zl[ind_lb] <<- reset_bound_dual(zl[ind_lb], x[ind_lb] - xl[ind_lb], RR$mu_R)
      zu[ind_ub] <<- reset_bound_dual(zu[ind_ub], xu[ind_ub] - x[ind_ub], RR$mu_R)
      RR$zp <<- reset_bound_dual(RR$zp, RR$pp, RR$mu_R)
      RR$zn <<- reset_bound_dual(RR$zn, RR$nn, RR$mu_R)
      adjust_boundary(mu)
      obj_val <<- eval_f(x)
      fgrad <<- eval_grad(x)
      theta <- sum(abs(cval))
      varphi <- varphi_of(obj_val, x, mu)
      if (filter_ok(filter, theta, varphi) && theta <= o$required_infeasibility_reduction * RR$theta_ref) {
        pv <<- numeric(nd)
        pv[iP] <<- -fgrad + zl - zu
        initialize_kkt()
        factorize()
        d <<- solve_refine(pv)$x
        y <<- if (max(abs(d[iY]), 0) > o$constr_mult_init_max) numeric(m) else d[iY]
        cnt$k <<- cnt$k + 1L
        cnt$t <<- cnt$t + 1L
        return("REGULAR")
      }
      if (cnt$k >= o$max_iter) return("MAXIMUM_ITERATIONS_EXCEEDED")
      if (wall()) return("MAXIMUM_WALLTIME_EXCEEDED")
      cnt$k <<- cnt$k + 1L
      cnt$t <<- cnt$t + 1L
    }
  }

  # solve!
  status <- tryCatch({
    if (nvar < n_eq) {
      "NOT_ENOUGH_DEGREES_OF_FREEDOM"
    } else {
      st <- initialize()
      while (st %in% c("REGULAR", "RESTORE", "ROBUST")) {
        if (st == "REGULAR") st <- regular()
        if (st == "RESTORE") st <- restore()
        if (st == "ROBUST") st <- robust()
      }
      st
    }
  }, madnlp_invalid = function(e) e$what)
  # MadNLPExecutionStats: the unscaled solution
  xs <- out_x(x)
  zl_full <- numeric(nvar); zu_full <- numeric(nvar)
  zl_full[free] <- zl[seq_len(n)] / obj_scale
  zu_full[free] <- zu[seq_len(n)] / obj_scale
  cv <- if (m) pb$cons(xs) else numeric()
  viol <- if (m) max(pmax(pb$lcon - cv, cv - pb$ucon, 0), 0) else 0
  list(status = status, iterations = cnt$k, x = xs, objective = obj_val / obj_scale,
       lambda = y * con_scale / obj_scale, zl = zl_full, zu = zu_full, violation = viol)
}

# CompactLBFGS (src/quasi_newton.jl): the compact representation
# B = sigma I - U U' + V V' of the limited memory BFGS approximation, with
# the options of QuasiNewtonOptions.
.madnlp_lbfgs_new <- function(n, opts = NULL) {
  qo <- utils::modifyList(list(init_strategy = "SCALAR1", max_history = 6L, init_value = 1, sigma_min = 1e-8,
                               sigma_max = 1e8), as.list(opts %||% list()))
  list(init_strategy = sub("^MadNLP\\.", "", qo$init_strategy), max_mem = as.integer(qo$max_history),
       init_value = qo$init_value, sigma_min = qo$sigma_min, sigma_max = qo$sigma_max, n = n,
       current_mem = 0L, skipped_iter = 0L, Sk = matrix(0, n, 0L), Yk = matrix(0, n, 0L),
       Lk = matrix(0, 0L, 0L), Dk = numeric(), U = matrix(0, n, 0L), V = matrix(0, n, 0L),
       max_mem_reached = FALSE, Bk = 0, last_x = numeric(n), last_g = numeric(n))
}

.madnlp_lbfgs_reset <- function(qn) {
  qn$current_mem <- 0L
  qn$skipped_iter <- 0L
  qn$Dk <- numeric()
  qn$Sk <- matrix(0, qn$n, 0L); qn$Yk <- matrix(0, qn$n, 0L)
  qn$Lk <- matrix(0, 0L, 0L)
  qn$U <- matrix(0, qn$n, 0L); qn$V <- matrix(0, qn$n, 0L)
  qn$max_mem_reached <- FALSE
  qn
}

.madnlp_lbfgs_init <- function(qn, g0, f0) {
  ng <- sum(g0 * g0)
  rho0 <- if (ng < sqrt(.Machine$double.eps)) 1 else if (f0 == 0) 1 / ng else abs(f0) / ng
  qn$Bk <- 2 * rho0 * qn$init_value
  qn
}

.madnlp_lbfgs_update <- function(qn, sk, yk) {
  eps <- .Machine$double.eps
  ns <- sqrt(sum(sk * sk)); ny <- sqrt(sum(yk * yk))
  if (ns < 100 * eps || ny < 100 * eps || sum(sk * yk) < sqrt(eps) * ns * ny) {
    qn$skipped_iter <- qn$skipped_iter + 1L
    if (qn$skipped_iter >= 2L) qn <- .madnlp_lbfgs_reset(qn)
    return(qn)
  }
  n <- qn$n
  if (qn$current_mem < qn$max_mem) {
    qn$current_mem <- qn$current_mem + 1L
    qn$Sk <- cbind(qn$Sk, sk); qn$Yk <- cbind(qn$Yk, yk)
    k <- qn$current_mem
    qn$Lk <- matrix(0, k, k)
  } else {
    k <- qn$current_mem
    if (k > 1L) {
      qn$Sk[, seq_len(k - 1L)] <- qn$Sk[, 2:k]
      qn$Yk[, seq_len(k - 1L)] <- qn$Yk[, 2:k]
    }
    qn$Sk[, k] <- sk; qn$Yk[, k] <- yk
  }
  dimnames(qn$Sk) <- NULL; dimnames(qn$Yk) <- NULL
  if (qn$max_mem_reached) {
    for (i in seq_len(k - 1L)) {
      qn$Dk[i] <- qn$Dk[i + 1L]
      for (j in seq_len(i - 1L)) qn$Lk[i, j] <- qn$Lk[i + 1L, j + 1L]
    }
    qn$Lk[k, ] <- as.numeric(crossprod(qn$Yk, qn$Sk[, k]))
    qn$Dk[k] <- qn$Lk[k, k]
    qn$Lk[k, k] <- 0
  } else {
    L <- crossprod(qn$Sk, qn$Yk)
    qn$Dk <- c(qn$Dk, L[k, k])
    L[upper.tri(L, diag = TRUE)] <- 0
    qn$Lk <- L
    if (qn$current_mem == qn$max_mem) qn$max_mem_reached <- TRUE
  }
  sty <- sum(sk * yk); sts <- sum(sk * sk); yty <- sum(yk * yk)
  sigma <- switch(qn$init_strategy, SCALAR1 = sty / sts, SCALAR2 = yty / sty,
                  SCALAR3 = ((sty / sts) + (yty / sty)) / 2, SCALAR4 = sqrt((sty / sts) * (yty / sty)),
                  stop("unknown BFGS init_strategy ", qn$init_strategy, call. = FALSE))
  sigma <- min(max(sigma, qn$sigma_min), qn$sigma_max)
  qn$Bk <- sigma
  delta <- 1 / sqrt(qn$Dk)
  DkLk <- delta * t(qn$Lk)
  Mk <- sigma * crossprod(qn$Sk) + crossprod(DkLk)
  R <- chol(Mk)
  qn$V <- sweep(qn$Yk, 2L, delta, `*`)
  U <- sigma * qn$Sk + qn$V %*% DkLk
  qn$U <- t(forwardsolve(t(R), t(U)))
  qn
}
