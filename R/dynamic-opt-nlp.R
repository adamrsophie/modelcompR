# The nonlinear program solver of the native dynamic optimization backend: a
# primal dual interior point method for
#
#   min f(z)  subject to  cE(z) = 0,  cI(z) <= 0,  lo <= z <= hi,
#
# after Waechter and Biegler, "On the implementation of an interior-point
# filter line-search algorithm for large-scale nonlinear programming" (2006),
# the algorithm of Ipopt, which the JuMP, InfiniteOpt, CasADi and Pyomo
# backends of ModelingToolkit (ext/MTKInfiniteOptExt.jl and the others) call.
# Inequalities get slacks s >= 0 and bounds a logarithmic barrier with
# primal dual bound multipliers (equations 13 to 16 of the paper). The Newton
# step solves the KKT system with Ipopt's inertia correction (the inertia
# from LAPACK's Bunch-Kaufman factorization; the curvature test of Chiang
# and Zavala above .nlp_inertia_max_dim) and iterative refinement, steps
# stay inside the bounds by the fraction to the boundary rule, and the
# backtracking filter line search of Ipopt (second order corrections,
# watchdog, soft restoration, filter reset, tiny steps) accepts them. The
# barrier parameter follows Ipopt's monotone or adaptive update (quality
# function, probing or LOQO oracle), the Hessian is exact or Ipopt's limited
# memory BFGS or SR1, and the problem is scaled by gradients at the start.
# When the line search fails, Ipopt's restoration phase (RestoMinC_1Nrm)
# minimizes the 1-norm of the infeasibility with a proximity term as an
# interior point problem of its own, with its own filter line search.
#
# The KKT matrix of a collocation problem is banded when unknowns and rows
# are ordered in time: `keys` gives the time position of every unknown and
# row (Inf for global ones such as tuned parameters), the banded part is
# factorized by Gaussian elimination with partial pivoting inside the band,
# and the global unknowns enter through a Schur complement.

# Solve K x = B for a matrix K whose nonzeros lie within `bw` of the
# diagonal; partial pivoting widens the upper band to 2 bw.
.nlp_band_solve <- function(K, B, bw) {
  n <- nrow(K)
  B <- as.matrix(B)
  if (!n) return(B)
  bu <- 2L * bw
  for (j in seq_len(n)) {
    rl <- min(n, j + bw)
    cu <- min(n, j + bu)
    if (rl > j) {
      cand <- j:rl
      p <- cand[which.max(abs(K[cand, j]))]
      if (p != j) {
        cols <- j:cu
        tmp <- K[j, cols]; K[j, cols] <- K[p, cols]; K[p, cols] <- tmp
        tmp <- B[j, ]; B[j, ] <- B[p, ]; B[p, ] <- tmp
      }
    }
    # a zero pivot means a singular matrix; the caller regularizes it
    if (K[j, j] == 0) stop("singular KKT matrix", call. = FALSE)
    if (rl > j) {
      r <- (j + 1L):rl
      l <- K[r, j] / K[j, j]
      nzl <- l != 0
      if (any(nzl)) {
        r <- r[nzl]; l <- l[nzl]
        if (cu > j) {
          cols <- (j + 1L):cu
          K[r, cols] <- K[r, cols, drop = FALSE] - outer(l, K[j, cols])
        }
        B[r, ] <- B[r, , drop = FALSE] - outer(l, B[j, ])
      }
    }
  }
  for (i in n:1) {
    cu <- min(n, i + bu)
    if (cu > i) {
      cols <- (i + 1L):cu
      B[i, ] <- B[i, ] - colSums(K[i, cols] * B[cols, , drop = FALSE])
    }
    B[i, ] <- B[i, ] / K[i, i]
  }
  B
}

# Solve the KKT system given as triplets (i, j, v) of a square matrix of
# dimension `dim` (duplicates are summed) and right-hand side r. `key` is
# the time position of each index; indices with key Inf form the border.
.nlp_kkt_solve <- function(ti, tj, tv, dim, r, key) {
  ord <- order(key, seq_len(dim))
  pos <- integer(dim)
  pos[ord] <- seq_len(dim)
  border <- is.infinite(key)
  nb <- sum(border)
  nc <- dim - nb
  pi <- pos[ti]; pj <- pos[tj]
  core <- pi <= nc & pj <= nc
  bw <- if (any(core)) max(abs(pi[core] - pj[core]), 1L) else 1L
  K <- matrix(0, dim, dim)
  lin <- (pj - 1) * dim + pi
  s <- rowsum(tv, lin)
  K[as.numeric(rownames(s))] <- s[, 1L]
  rr <- numeric(dim)
  rr[pos] <- r
  if (!nb) {
    x <- .nlp_band_solve(K, rr, bw)[, 1L]
  } else {
    ci <- seq_len(nc); bi <- nc + seq_len(nb)
    Y <- .nlp_band_solve(K[ci, ci, drop = FALSE], cbind(rr[ci], K[ci, bi, drop = FALSE]), bw)
    x1 <- Y[, 1L]
    Yb <- Y[, -1L, drop = FALSE]
    S <- K[bi, bi, drop = FALSE] - K[bi, ci, drop = FALSE] %*% Yb
    x2 <- tryCatch(solve(S, rr[bi] - K[bi, ci, drop = FALSE] %*% x1), error = function(e) NULL)
    if (is.null(x2)) stop("singular KKT matrix", call. = FALSE)
    x <- c(x1 - as.numeric(Yb %*% x2), as.numeric(x2))
  }
  x[pos]
}

# The inertia of the symmetric matrix given as triplets (duplicates are
# summed) from the Bunch-Kaufman factorization of LAPACK: the counts of
# positive, negative and zero pivots, 2 by 2 blocks counted by their
# eigenvalue signs. Only an exactly zero pivot is zero, as MUMPS reports a
# singular matrix to Ipopt (no null pivot detection in the factorization).
.nlp_inertia_max_dim <- 2000L
.nlp_inertia <- function(ti, tj, tv, dim) {
  K <- matrix(0, dim, dim)
  s <- rowsum(tv, (tj - 1) * dim + ti)
  K[as.numeric(rownames(s))] <- s[, 1L]
  K <- (K + t(K)) / 2
  tol <- 0
  # a singular matrix gives zero pivots, which LAPACK reports by a warning
  bk <- suppressWarnings(Matrix::BunchKaufman(methods::new("dsyMatrix", Dim = c(dim, dim), x = as.numeric(K),
                                                           uplo = "U")))
  X <- matrix(bk@x, dim)
  piv <- bk@perm
  pos <- 0L; neg <- 0L; zero <- 0L
  k <- dim
  while (k >= 1L) {
    if (piv[k] > 0L || k == 1L) {
      dk <- X[k, k]
      if (abs(dk) <= tol) zero <- zero + 1L else if (dk > 0) pos <- pos + 1L else neg <- neg + 1L
      k <- k - 1L
    } else {
      a <- X[k - 1L, k - 1L]; b <- X[k - 1L, k]; cc <- X[k, k]
      dt <- a * cc - b * b
      if (abs(dt) <= tol * tol) {
        zero <- zero + 1L
        if (a + cc > 0) pos <- pos + 1L else neg <- neg + 1L
      } else if (dt < 0) {
        pos <- pos + 1L; neg <- neg + 1L
      } else if (a + cc > 0) {
        pos <- pos + 2L
      } else {
        neg <- neg + 2L
      }
      k <- k - 2L
    }
  }
  list(pos = pos, neg = neg, zero = zero)
}

# sum(v[idx == k]) for k in 1:n.
.nlp_accum <- function(idx, v, n) {
  out <- numeric(n)
  if (length(idx)) {
    a <- rowsum(v, idx)
    out[as.integer(rownames(a))] <- a[, 1L]
  }
  out
}

# Start strictly inside the bounds (Ipopt's bound_push and bound_frac).
# Start strictly inside the bounds (Ipopt's bound_push and bound_frac). The
# push is relative to the size of `ref_lo` and `ref_hi`, by default the
# bounds themselves.
.nlp_push_inside <- function(z, lo, hi, k1 = 1e-2, k2 = 1e-2, ref_lo = lo, ref_hi = hi) {
  pl <- ifelse(is.finite(lo), pmin(k1 * pmax(1, abs(ref_lo)), ifelse(is.finite(hi), k2 * (hi - lo), Inf)), 0)
  pu <- ifelse(is.finite(hi), pmin(k1 * pmax(1, abs(ref_hi)), ifelse(is.finite(lo), k2 * (hi - lo), Inf)), 0)
  z <- ifelse(is.finite(lo), pmax(z, lo + pl), z)
  ifelse(is.finite(hi), pmin(z, hi - pu), z)
}

# The largest step in (0, 1] keeping x + a dx at least (1 - tau) of the way
# from its bound (the fraction to the boundary rule).
.nlp_ftb <- function(gap, dgap, tau) {
  neg <- dgap < 0 & is.finite(gap)
  if (!any(neg)) return(1)
  min(1, min(-tau * gap[neg] / dgap[neg]))
}

# Ipopt's options of the algorithm with their default values (the names and
# defaults of Ipopt 3.14's options; `variant` selects the restoration phase
# of Ipopt or of MadNLP).
.nlp_ip_defaults <- function() {
  list(acceptable_tol = 1e-6, acceptable_iter = 15L, dual_inf_tol = 1, constr_viol_tol = 1e-4,
       compl_inf_tol = 1e-4, acceptable_dual_inf_tol = 1e10, acceptable_constr_viol_tol = 1e-2,
       acceptable_compl_inf_tol = 1e-2, acceptable_obj_change_tol = 1e20, diverging_iterates_tol = 1e20,
       mu_init = 0.1, mu_target = 0, mu_linear_decrease_factor = 0.2, mu_superlinear_decrease_power = 1.5,
       mu_allow_fast_monotone_decrease = TRUE, barrier_tol_factor = 10, tau_min = 0.99, bound_push = 1e-2,
       bound_frac = 1e-2, slack_bound_push = 1e-2, slack_bound_frac = 1e-2, bound_relax_factor = 1e-8,
       honor_original_bounds = FALSE, bound_mult_init_val = 1, constr_mult_init_max = 1e3,
       nlp_scaling_method = "gradient-based", nlp_scaling_max_gradient = 100, nlp_scaling_min_value = 1e-8,
       obj_scaling_factor = 1, kappa_d = 1e-5, kappa_sigma = 1e10, s_max = 100, theta_max_fact = 1e4,
       theta_min_fact = 1e-4, eta_phi = 1e-8, delta = 1, s_phi = 2.3, s_theta = 1.1, gamma_phi = 1e-8,
       gamma_theta = 1e-5, alpha_min_frac = 0.05, max_soc = 4L, kappa_soc = 0.99, obj_max_inc = 5,
       max_filter_resets = 5L, filter_reset_trigger = 5L, alpha_red_factor = 0.5,
       watchdog_shortened_iter_trigger = 10L, watchdog_trial_iter_max = 3L, accept_every_trial_step = FALSE,
       accept_after_max_steps = -1L, tiny_step_tol = 10 * .Machine$double.eps, tiny_step_y_tol = 1e-2,
       soft_resto_pderror_reduction_factor = 0.9999, max_soft_resto_iters = 10L,
       expect_infeasible_problem = FALSE, expect_infeasible_problem_ctol = 1e-3,
       expect_infeasible_problem_ytol = 1e8, start_with_resto = FALSE, required_infeasibility_reduction = 0.9,
       max_resto_iter = 3000000L, resto_penalty_parameter = 1000, resto_proximity_weight = 1,
       bound_mult_reset_threshold = 1000, constr_mult_reset_threshold = 0,
       resto_failure_feasibility_threshold = NULL, first_hessian_perturbation = 1e-4,
       min_hessian_perturbation = 1e-20, max_hessian_perturbation = 1e20, perturb_inc_fact_first = 100,
       perturb_inc_fact = 8, perturb_dec_fact = 1 / 3, jacobian_regularization_value = 1e-8,
       jacobian_regularization_exponent = 0.25, perturb_always_cd = FALSE, hessian_approximation = "exact",
       limited_memory_max_history = 6L, limited_memory_update_type = "bfgs", max_wall_time = 1e20,
       mu_strategy = "monotone", mu_oracle = "quality-function", fixed_mu_oracle = "average_compl",
       adaptive_mu_globalization = "obj-constr-filter", mu_max_fact = 1e3, mu_max = NA_real_,
       mu_min = NA_real_, adaptive_mu_safeguard_factor = 0, adaptive_mu_kkterror_red_iters = 4L,
       adaptive_mu_kkterror_red_fact = 0.9999, filter_margin_fact = 1e-5, filter_max_margin = 1,
       adaptive_mu_restore_previous_iterate = FALSE, adaptive_mu_monotone_init_factor = 0.8,
       sigma_max = 100, sigma_min = 1e-6, quality_function_norm_type = "2-norm-squared",
       quality_function_centrality = "none", quality_function_balancing_term = "none",
       quality_function_max_section_steps = 8L, quality_function_section_sigma_tol = 1e-2,
       quality_function_section_qf_tol = 0, variant = "ipopt")
}

# The primal dual interior point method with a filter line search of Ipopt
# (Waechter and Biegler 2006) on
#
#   min f(z)  subject to  cE(z) = 0,  cI(z) <= 0,  lo <= z <= hi.
#
# `evalf(z, level, lamE, lamI, obj)` returns f, cE and cI, with level >= 1
# also the gradient gf and the Jacobians JE and JI as triplets and with
# level 2 the Hessian H of obj f + lamE' cE + lamI' cI as triplets. The keys
# give the time position of every unknown and row for the banded KKT solve.
# `rowI_bound` holds the original bound of every inequality row (the size of
# the slack's bound for the bound push and relaxation). With `resto` the
# call is the restoration phase of an outer call: it skips the scaling, the
# bound relaxation and the push and starts from the iterate in `resto`.
#
# The status is one of Success, Acceptable, FeasiblePointFound, Infeasible
# (local infeasibility), TinyStep, Diverging, MaxIters, MaxTime,
# RestorationFailed, StepFailed, InvalidNumber and TooFewDof (Ipopt's
# SolverReturn values).
.nlp_interior_point <- function(evalf, z0, lo, hi, var_key, rowE_key, rowI_key, tol = 1e-8, maxiter = 3000L,
                                verbose = FALSE, opts = list(), rowI_bound = NULL, resto = NULL) {
  o <- utils::modifyList(.nlp_ip_defaults(), opts, keep.null = TRUE)
  o$tol <- tol
  o$max_iter <- as.integer(maxiter)
  start_time <- proc.time()[["elapsed"]]
  inner <- !is.null(resto)
  madnlp <- identical(o$variant, "madnlp")
  nz <- length(z0)
  lo <- rep_len(as.numeric(lo), nz); hi <- rep_len(as.numeric(hi), nz)
  eps <- .Machine$double.eps
  cmp_le <- function(lhs, rhs, bas) lhs - rhs <= 10 * eps * abs(bas)
  if (!inner) {
    # fixed variables become parameters (fixed_variable_treatment =
    # make_parameter); the other bounds are relaxed by bound_relax_factor
    fixed <- is.finite(lo) & is.finite(hi) & lo == hi
    z0[fixed] <- lo[fixed]
    lo_orig <- lo; hi_orig <- hi
    brf <- o$bound_relax_factor
    rl <- is.finite(lo) & !fixed; ru <- is.finite(hi) & !fixed
    lo[rl] <- lo[rl] - pmin(brf * pmax(1, abs(lo[rl])), o$constr_viol_tol)
    hi[ru] <- hi[ru] + pmin(brf * pmax(1, abs(hi[ru])), o$constr_viol_tol)
    # gradient based scaling at the starting point (nlp_scaling_method)
    e1 <- evalf(z0, 1L)
    nE <- length(e1$cE); nI <- length(e1$cI)
    if (!all(is.finite(c(e1$f, e1$cE, e1$cI, e1$gf, e1$JE$v, e1$JI$v))))
      return(list(z = z0, f = e1$f, lamE = numeric(nE), lamI = numeric(nI), violation = Inf, iterations = 0L,
                  status = "InvalidNumber"))
    rowmax <- function(J, nr) {
      out <- numeric(nr)
      keep <- !fixed[J$j]
      if (any(keep)) {
        a <- tapply(abs(J$v[keep]), J$i[keep], max)
        out[as.integer(names(a))] <- a
      }
      out
    }
    sf <- o$obj_scaling_factor; sE <- rep(1, nE); sI <- rep(1, nI)
    if (identical(o$nlp_scaling_method, "gradient-based")) {
      mg <- o$nlp_scaling_max_gradient; mv <- o$nlp_scaling_min_value
      gmax <- max(abs(e1$gf[!fixed]), 0)
      if (gmax > mg) sf <- sf * max(mv, mg / gmax)
      rE <- rowmax(e1$JE, nE); rI <- rowmax(e1$JI, nI)
      sE <- ifelse(rE > mg, pmax(mv, mg / rE), 1)
      sI <- ifelse(rI > mg, pmax(mv, mg / rI), 1)
    }
    bI <- if (is.null(rowI_bound)) numeric(nI) else rep_len(as.numeric(rowI_bound), nI)
    slack_lo <- -pmin(brf * pmax(1, abs(sI * bI)), o$constr_viol_tol)
  } else {
    fixed <- resto$fixed
    lo_orig <- lo; hi_orig <- hi
    sf <- 1
    e1 <- evalf(z0, 0L)
    nE <- length(e1$cE); nI <- length(e1$cI)
    sE <- rep(1, nE); sI <- rep(1, nI)
    slack_lo <- resto$slack_lo
    bI <- numeric(nI)
  }
  m <- nE + nI
  ny <- nz + nI
  nfree <- sum(!fixed)
  if (!inner && nE > nfree)
    return(list(z = z0, f = e1$f, lamE = numeric(nE), lamI = numeric(nI), violation = Inf, iterations = 0L,
                status = "TooFewDof"))
  sev <- function(z, level, lamE = NULL, lamI = NULL, obj = 1) {
    e <- evalf(z, level, if (!is.null(lamE)) lamE * sE, if (!is.null(lamI)) lamI * sI, sf * obj)
    e$f <- sf * e$f; e$cE <- sE * e$cE; e$cI <- sI * e$cI
    if (level >= 1L) {
      e$gf <- sf * e$gf
      e$gf[fixed] <- 0
      e$JE$v <- sE[e$JE$i] * e$JE$v
      e$JI$v <- sI[e$JI$i] * e$JI$v
      e$JE$v[fixed[e$JE$j]] <- 0
      e$JI$v[fixed[e$JI$j]] <- 0
    }
    if (level >= 2L) e$H$v[fixed[e$H$i] | fixed[e$H$j]] <- 0
    e
  }
  ylo <- c(lo, slack_lo); yhi <- c(hi, rep(Inf, nI))
  free <- c(!fixed, rep(TRUE, nI))
  hl <- is.finite(ylo) & free; hu <- is.finite(yhi) & free
  lonly <- hl & !hu; uonly <- hu & !hl
  nbnd <- sum(hl) + sum(hu)
  key <- c(var_key, rowI_key, rowE_key, rowI_key)
  dim <- ny + m
  ix <- seq_len(nz); isl <- nz + seq_len(nI)
  cons_of <- function(e, yy) c(e$cE, e$cI + yy[isl])
  gaps <- function(yy) list(l = ifelse(hl, yy - ylo, Inf), u = ifelse(hu, yhi - yy, Inf))
  # the barrier objective with Ipopt's damping of one sided bounds (kappa_d)
  barrier_of <- function(f, yy, mu) {
    g <- gaps(yy)
    if (any(hl[] & g$l <= 0) || any(hu & g$u <= 0)) return(Inf)
    f - mu * (sum(log(g$l[hl])) + sum(log(g$u[hu]))) +
      o$kappa_d * mu * (sum(g$l[lonly]) + sum(g$u[uonly]))
  }
  jac_t <- function(e, lam) {
    c(.nlp_accum(e$JE$j, lam[e$JE$i] * e$JE$v, nz) + .nlp_accum(e$JI$j, lam[nE + e$JI$i] * e$JI$v, nz),
      lam[nE + seq_len(nI)])
  }
  ji_of <- function(e) list(i = c(e$JE$i, nE + e$JI$i, nE + seq_len(nI)), j = c(e$JE$j, e$JI$j, isl),
                            v = c(e$JE$v, e$JI$v, rep(1, nI)))
  # the KKT matrix [W + Sigma + dw, J'; J, -dc] as triplets; a fixed unknown
  # keeps the unit row and column of a parameter
  kkt_triplets <- function(H, sig, dw, J, dc) {
    d <- ifelse(free, sig + dw, 1)
    ti <- c(H$i, seq_len(ny), ny + J$i, J$j)
    tj <- c(H$j, seq_len(ny), J$j, ny + J$i)
    tv <- c(H$v, d, J$v, J$v)
    if (m) {
      ti <- c(ti, ny + seq_len(m)); tj <- c(tj, ny + seq_len(m)); tv <- c(tv, rep(-dc, m))
    }
    list(i = ti, j = tj, v = tv)
  }
  # iterative refinement as in Ipopt's PDFullSpaceSolver: at least one and at
  # most ten corrections, until the residual ratio is below 1e-10
  ksolve <- function(K, r) {
    x <- .nlp_kkt_solve(K$i, K$j, K$v, dim, r, key)
    nr <- max(abs(r), 0)
    for (k in seq_len(10L)) {
      res <- r - .nlp_accum(K$i, K$v * x[K$j], dim)
      ratio <- max(abs(res), 0) / (min(max(abs(x), 0), 1e6 * nr) + nr)
      if (k > 1L && ratio <= 1e-10) break
      if (!is.finite(ratio)) break
      x <- x + .nlp_kkt_solve(K$i, K$j, K$v, dim, res, key)
    }
    x
  }
  # least squares multipliers of the constraints (constr_mult_init_max)
  ls_mult <- function(e, yy, zl, zu, cap) {
    if (!m || cap <= 0) return(numeric(m))
    J <- ji_of(e)
    gl <- c(e$gf, numeric(nI))
    K <- kkt_triplets(list(i = integer(), j = integer(), v = numeric()), rep(1, ny), 0, J, 0)
    r <- c(ifelse(free, -(gl - zl + zu), 0), numeric(m))
    d <- tryCatch(ksolve(K, r), error = function(err) NULL)
    if (is.null(d) || !all(is.finite(d))) return(numeric(m))
    lam <- d[ny + seq_len(m)]
    if (max(abs(lam)) > cap) numeric(m) else lam
  }
  # starting point
  if (!inner) {
    # the warm start initializer starts from multipliers given by the
    # caller, zero here (as Ipopt.jl passes them), pushed to 1e-3; Ipopt
    # 3.14 through Ipopt.jl starts the primal variables at the point of the
    # default initializer (its iteration 0 equals the cold start's)
    warm <- isTRUE(o$warm_start_init_point)
    z <- z0
    z[!fixed] <- .nlp_push_inside(z0, lo, hi, o$bound_push, o$bound_frac)[!fixed]
    ev <- sev(z, 1L)
    s <- .nlp_push_inside(-ev$cI, slack_lo, rep(Inf, nI), o$slack_bound_push, o$slack_bound_frac,
                          ref_lo = pmax(abs(sI * bI), abs(slack_lo)))
    y <- c(z, s)
    zinit <- if (warm) 1e-3 else o$bound_mult_init_val
    zl <- ifelse(hl, zinit, 0); zu <- ifelse(hu, zinit, 0)
    lam <- if (madnlp || warm) numeric(m) else ls_mult(ev, y, zl, zu, o$constr_mult_init_max)
    mu <- o$mu_init
  } else {
    y <- resto$y0
    zl <- resto$zl0; zu <- resto$zu0
    mu <- resto$mu0
    ev <- sev(y[ix], 1L)
    lam <- if (madnlp) numeric(m) else ls_mult(ev, y, zl, zu, o$constr_mult_init_max)
  }
  tau <- max(o$tau_min, 1 - mu)
  # the restoration objective depends on the barrier parameter
  set_resto_mu <- function(mu) if (inner) resto$env$mu <- mu
  set_resto_mu(mu)
  # line search state
  filt <- matrix(numeric(), 0L, 2L)
  theta_max <- -1; theta_min <- -1
  ref <- NULL
  last_rej_filter <- FALSE; n_filter_rej <- 0L; n_filter_resets <- 0L
  filter_reset <- function() {
    filt <<- matrix(numeric(), 0L, 2L)
    last_rej_filter <<- FALSE
    n_filter_rej <<- 0L
  }
  add_filter <- function(theta, phi) {
    th <- (1 - o$gamma_theta) * theta; ph <- phi - o$gamma_phi * theta
    if (nrow(filt)) filt <<- filt[!(filt[, 1L] >= th & filt[, 2L] >= ph), , drop = FALSE]
    filt <<- rbind(filt, c(th, ph))
  }
  filter_ok <- function(theta, phi) !nrow(filt) || all(theta <= filt[, 1L] | phi <= filt[, 2L])
  is_ftype <- function(a) {
    gbd <- ref$gbd
    if (ref$theta == 0 && gbd > 0 && gbd < 100 * eps) gbd <- -eps
    gbd < 0 && a * (-gbd)^o$s_phi > o$delta * ref$theta^o$s_theta
  }
  armijo <- function(barr, a) cmp_le(barr - ref$barr, o$eta_phi * a * ref$gbd, ref$barr)
  ok_to_iterate <- function(barr, theta, from_resto = FALSE) {
    if (!from_resto && barr > ref$barr) {
      basval <- if (abs(ref$barr) > 10) log10(abs(ref$barr)) else 1
      if (log10(barr - ref$barr) > o$obj_max_inc + basval) return(FALSE)
    }
    cmp_le(theta, (1 - o$gamma_theta) * ref$theta, ref$theta) ||
      cmp_le(barr - ref$barr, -o$gamma_phi * ref$theta, ref$barr)
  }
  # CheckAcceptabilityOfTrialPoint of the filter acceptor
  check_trial <- function(tr, a) {
    if (theta_max < 0) theta_max <<- o$theta_max_fact * max(1, ref$theta)
    if (theta_min < 0) theta_min <<- o$theta_min_fact * max(1, ref$theta)
    if (!is.finite(tr$theta) || !is.finite(tr$barr)) return(FALSE)
    if (theta_max > 0 && tr$theta > theta_max) return(FALSE)
    acc <- if (a > 0 && is_ftype(a) && ref$theta <= theta_min) armijo(tr$barr, a) else
      ok_to_iterate(tr$barr, tr$theta)
    if (!acc) {
      last_rej_filter <<- FALSE
      return(FALSE)
    }
    if (!filter_ok(tr$theta, tr$barr)) {
      last_rej_filter <<- TRUE
      return(FALSE)
    }
    if (o$max_filter_resets > 0 && n_filter_resets < o$max_filter_resets) {
      if (last_rej_filter) {
        n_filter_rej <<- n_filter_rej + 1L
        if (n_filter_rej >= o$filter_reset_trigger) {
          filter_reset()
          n_filter_resets <<- n_filter_resets + 1L
        }
      } else {
        n_filter_rej <<- 0L
      }
    }
    last_rej_filter <<- FALSE
    TRUE
  }
  trial_at <- function(yt, level = 0L) {
    et <- tryCatch(sev(yt[ix], level), error = function(err) NULL)
    if (is.null(et) || !is.finite(et$f) || !all(is.finite(c(et$cE, et$cI))))
      return(list(y = yt, theta = Inf, barr = Inf, e = NULL, c = NULL))
    ct <- cons_of(et, yt)
    list(y = yt, theta = sum(abs(ct)), barr = barrier_of(et$f, yt, mu), e = et, c = ct)
  }
  # quantities of the optimality error at an iterate
  measures <- function(e, yy, lam, zl, zu, mu_c = 0) {
    g <- gaps(yy)
    dual <- c(e$gf, numeric(nI)) + jac_t(e, lam) - zl + zu
    dual[!free] <- 0
    cc <- cons_of(e, yy)
    cl <- (g$l * zl)[hl] - mu_c; cu <- (g$u * zu)[hu] - mu_c
    sd <- max(o$s_max, (sum(abs(lam)) + sum(abs(zl)) + sum(abs(zu))) / max(1, m + nbnd)) / o$s_max
    sc <- max(o$s_max, (sum(abs(zl)) + sum(abs(zu))) / max(1, nbnd)) / o$s_max
    list(dual = dual, c = cc, theta = sum(abs(cc)), inf_pr = max(abs(cc), 0), inf_du = max(abs(dual), 0),
         compl = max(abs(cl), abs(cu), 0), compl1 = sum(abs(cl)) + sum(abs(cu)), sd = sd, sc = sc)
  }
  pd_error <- function(e, yy, lam, zl, zu) {
    q <- measures(e, yy, lam, zl, zu, mu)
    sum(abs(q$dual)) / max(1, nfree + nI) + (if (m) sum(abs(q$c)) / m else 0) +
      (if (nbnd) q$compl1 / nbnd else 0)
  }
  unscaled <- function(e, q, yy) {
    du <- c(abs(q$dual[ix]) / sf, abs(q$dual[isl]) * sI / sf)
    # Ipopt's unscaled constraint violation measures the inequalities
    # against their relaxed bounds and leaves out the variable bounds
    viol <- max(abs(e$cE / sE), pmax((e$cI + slack_lo) / sI, 0), 0, na.rm = TRUE)
    list(dual = max(du, 0), viol = viol, compl = q$compl / sf)
  }
  step_of <- function(d, yy, zl, zu) {
    g <- gaps(yy)
    dy <- d[seq_len(ny)]
    list(dy = dy, dlam = d[ny + seq_len(m)],
         dzl = ifelse(hl, mu / g$l - zl - zl / g$l * dy, 0),
         dzu = ifelse(hu, mu / g$u - zu + zu / g$u * dy, 0))
  }
  primal_ftb <- function(yy, dy) {
    g <- gaps(yy)
    min(.nlp_ftb(g$l, dy, tau), .nlp_ftb(g$u, -dy, tau))
  }
  dual_ftb <- function(st, zl, zu) min(.nlp_ftb(zl, st$dzl, tau), .nlp_ftb(zu, st$dzu, tau))
  # keep the bound multipliers near the central path (kappa_sigma)
  clip_z <- function(yy, zl, zu) {
    g <- gaps(yy)
    ks <- o$kappa_sigma
    list(l = ifelse(hl, pmin(pmax(zl, mu / (ks * g$l)), ks * mu / g$l), 0),
         u = ifelse(hu, pmin(pmax(zu, mu / (ks * g$u)), ks * mu / g$u), 0))
  }
  # counters and states of Ipopt's algorithm objects
  status <- "MaxIters"
  it <- if (inner) resto$it0 else 0L
  acc_count <- 0L
  last_f <- NA_real_
  # the state of Ipopt's PDPerturbationHandler: degeneracy of the Hessian
  # and of the Jacobian ("undetermined", "degenerate", "not"), the test in
  # progress and the current and last perturbations
  pt <- new.env(parent = emptyenv())
  pt$always_cd <- isTRUE(o$perturb_always_cd)
  pt$hess <- "undetermined"
  pt$jac <- if (pt$always_cd) "not" else "undetermined"
  pt$degen_iters <- 0L
  pt$test <- "none"
  pt$dx <- 0; pt$dc <- 0; pt$dx_last <- 0; pt$dc_last <- 0
  tiny_flag <- FALSE; tiny_last <- FALSE
  in_wd <- FALSE; wd <- NULL; wd_shortened <- 0L; wd_trial <- 0L
  count_shortened <- 0L
  soft <- FALSE; soft_count <- 0L
  last_mu <- mu
  expect_inf <- isTRUE(o$expect_infeasible_problem)
  start_resto <- isTRUE(o$start_with_resto)
  acceptable_point <- NULL
  first_iter <- TRUE
  mu_initialized <- FALSE
  n_resto <- 0L
  # the first restoration iteration shows the line search information of
  # the failed iteration of the original problem (RestoMinC_1Nrm)
  info <- if (inner && !is.null(resto$info0)) resto$info0 else
    list(alpha_pr = 0, alpha_du = 0, char = " ", ls = 0L, dw = 0, dnorm = 0)
  skip_out <- FALSE
  out_line <- function(q, uq, e) {
    if (skip_out) {
      skip_out <<- FALSE
      return(invisible())
    }
    if (!verbose) return(invisible())
    # a restoration iteration shows the objective and the constraint
    # violation of the original problem (RestoIterationOutput)
    op <- if (inner) resto$orig_print(y) else list(f = e$f / sf, viol = uq$viol)
    message(sprintf("%4d%s %14.7e %7.2e %7.2e %5.1f %7.2e %5s %7.2e %7.2e%s %2d", it, if (inner) "r" else " ",
                    op$f, op$viol, q$inf_du, log10(mu), info$dnorm,
                    if (info$dw > 0) sprintf("%5.1f", log10(info$dw)) else "-", info$alpha_du, info$alpha_pr,
                    info$char, info$ls))
  }
  # the limited memory quasi Newton approximation of the Hessian of the
  # Lagrangian in the unknowns `lm_idx` (LimMemQuasiNewtonUpdater with the
  # scalar1 initialization), as a dense matrix
  lm <- identical(o$hessian_approximation, "limited-memory")
  lm_idx <- if (inner) resto$lm_idx else which(!fixed)
  lmS <- list(); lmY <- list(); lm_sigma <- 1; lm_last <- NULL; lm_skipped <- 0L
  lm_update <- function(e, yy, lam) {
    gnow <- e$gf + jac_t(e, lam)[ix]
    xnow <- yy[lm_idx]
    if (is.null(lm_last) || lm_skipped >= 2L) {
      lmS <<- list(); lmY <<- list(); lm_sigma <<- 1; lm_skipped <<- 0L
    } else {
      s <- xnow - lm_last$x
      glast <- lm_last$gf + jac_t(lm_last, lam)[ix]
      yv <- (gnow - glast)[lm_idx]
      sty <- sum(s * yv)
      skip <- max(abs(s), 0) < 100 * eps
      if (!skip && identical(o$limited_memory_update_type, "bfgs"))
        skip <- sty <= sqrt(eps) * sqrt(sum(s^2)) * sqrt(sum(yv^2))
      if (skip) {
        lm_skipped <<- lm_skipped + 1L
      } else {
        lm_skipped <<- 0L
        lmS <<- c(lmS, list(s)); lmY <<- c(lmY, list(yv))
        if (length(lmS) > o$limited_memory_max_history) {
          lmS <<- lmS[-1L]; lmY <<- lmY[-1L]
        }
        lm_sigma <<- if (identical(o$limited_memory_update_type, "sr1") && length(lmS) == 1L) 1 else
          max(1e-8, min(1e8, max(if (identical(o$limited_memory_update_type, "sr1")) 1e-8 else -Inf, abs(sty)) /
                                sum(s^2)))
      }
    }
    lm_last <<- list(x = xnow, gf = e$gf, JE = e$JE, JI = e$JI)
    k <- length(lm_idx)
    B <- diag(lm_sigma, k)
    if (length(lmS) && o$limited_memory_max_history > 0) {
      S <- do.call(cbind, lmS); Y <- do.call(cbind, lmY)
      SY <- crossprod(S, Y)
      Dd <- diag(diag(SY), ncol(S))
      L <- SY; L[upper.tri(L, diag = TRUE)] <- 0
      if (identical(o$limited_memory_update_type, "sr1")) {
        V <- Y - lm_sigma * S
        Z <- Dd + L + t(L) - lm_sigma * crossprod(S)
        Zi <- tryCatch(solve(Z), error = function(err) NULL)
        if (!is.null(Zi)) B <- B + V %*% Zi %*% t(V)
      } else {
        W <- cbind(lm_sigma * S, Y)
        M <- rbind(cbind(lm_sigma * crossprod(S), L), cbind(t(L), -Dd))
        Mi <- tryCatch(solve(M), error = function(err) NULL)
        if (!is.null(Mi)) B <- B - W %*% Mi %*% t(W)
      }
    }
    list(i = rep(lm_idx, times = k), j = rep(lm_idx, each = k), v = as.numeric(B))
  }
  # the right-hand side of the reduced KKT system for barrier parameter mm
  rhs_at <- function(mm) {
    g <- gaps(y)
    gphi <- c(ev$gf, numeric(nI)) - ifelse(hl, mm / g$l, 0) + ifelse(hu, mm / g$u, 0) +
      o$kappa_d * mm * (as.numeric(lonly) - as.numeric(uonly))
    gphi[!free] <- 0
    -c(ifelse(free, gphi + jac_t(ev, lam), 0), q$c)
  }
  # solve with the regularization (delta_w, delta_c) of Ipopt's perturbation
  # handler, raised until the matrix has the inertia (n, m, 0); above
  # .nlp_inertia_max_dim the curvature test of Chiang and Zavala stands in
  # PDPerturbationHandler of Ipopt 3.14: delta_cd, the deltas for a wrong
  # inertia, the degeneracy tests (finalize_test) and the perturbations for
  # a new system, a singular one and one with the wrong inertia
  pt_cd <- function() o$jacobian_regularization_value * mu^o$jacobian_regularization_exponent
  pt_wrong_deltas <- function() {
    pt$dx <- if (pt$dx == 0) {
      if (pt$dx_last == 0) o$first_hessian_perturbation else
        max(o$min_hessian_perturbation, pt$dx_last * o$perturb_dec_fact)
    } else if (pt$dx_last == 0 || 1e5 * pt$dx_last < pt$dx) {
      o$perturb_inc_fact_first * pt$dx
    } else {
      o$perturb_inc_fact * pt$dx
    }
    if (pt$dx > o$max_hessian_perturbation) {
      pt$dx_last <- 0
      return(FALSE)
    }
    TRUE
  }
  pt_finalize <- function() {
    switch(pt$test,
      none = NULL,
      c0x0 = {
        if (pt$hess == "undetermined") pt$hess <- "not"
        if (pt$jac == "undetermined") pt$jac <- "not"
      },
      c1x0 = {
        if (pt$hess == "undetermined") pt$hess <- "not"
        if (pt$jac == "undetermined") {
          pt$degen_iters <- pt$degen_iters + 1L
          if (pt$degen_iters >= 3L) pt$jac <- "degenerate"
        }
      },
      c0x1 = {
        if (pt$jac == "undetermined") pt$jac <- "not"
        if (pt$hess == "undetermined") {
          pt$degen_iters <- pt$degen_iters + 1L
          if (pt$degen_iters >= 3L) pt$hess <- "degenerate"
        }
      },
      c1x1 = {
        pt$degen_iters <- pt$degen_iters + 1L
        if (pt$degen_iters >= 3L) {
          pt$hess <- "degenerate"; pt$jac <- "degenerate"
        }
      })
    invisible()
  }
  pt_new_system <- function() {
    pt_finalize()
    if (pt$dx > 0) pt$dx_last <- pt$dx
    if (pt$dc > 0) pt$dc_last <- pt$dc
    pt$test <- if (pt$hess == "undetermined" || pt$jac == "undetermined") {
      if (pt$always_cd) "c1x0" else "c0x0"
    } else "none"
    pt$dc <- if (pt$jac == "degenerate" || pt$always_cd) pt_cd() else 0
    pt$dx <- 0
    if (pt$hess == "degenerate" && !pt_wrong_deltas()) return(FALSE)
    TRUE
  }
  pt_singular <- function() {
    if (pt$hess == "undetermined" || pt$jac == "undetermined") {
      switch(pt$test,
        c0x0 = {
          if (pt$jac == "undetermined") {
            pt$dc <- pt_cd(); pt$test <- "c1x0"
          } else {
            if (!pt_wrong_deltas()) return(FALSE)
            pt$test <- "c0x1"
          }
        },
        c1x0 = {
          if (!pt$always_cd) {
            pt$dc <- 0
            if (!pt_wrong_deltas()) return(FALSE)
            pt$test <- "c0x1"
          } else {
            if (!pt_wrong_deltas()) return(FALSE)
            pt$test <- "c1x1"
          }
        },
        c0x1 = {
          pt$dc <- pt_cd()
          if (!pt_wrong_deltas()) return(FALSE)
          pt$test <- "c1x1"
        },
        c1x1 = {
          if (!pt_wrong_deltas()) return(FALSE)
        })
    } else if (pt$dc > 0) {
      if (!pt_wrong_deltas()) return(FALSE)
    } else {
      pt$dc <- pt_cd()
    }
    TRUE
  }
  pt_wrong_inertia <- function() {
    pt_finalize()
    ok <- pt_wrong_deltas()
    if (!ok && pt$dc == 0) {
      pt$dc <- pt_cd(); pt$dx <- 0; pt$test <- "none"
      if (pt$hess == "degenerate") pt$hess <- "undetermined"
      ok <- pt_wrong_deltas()
    }
    ok
  }
  # PDFullSpaceSolver::SolveOnce: factorize with the current perturbations
  # until the inertia is right (m negative eigenvalues); a singular matrix
  # or too few negative eigenvalues ask for the perturbation for
  # singularity, too many for the one for a wrong inertia
  reg_solve <- function(rhs) {
    g <- gaps(y)
    sig <- ifelse(hl, zl / g$l, 0) + ifelse(hu, zu / g$u, 0)
    J <- ji_of(ev)
    d <- NULL
    K <- NULL
    ok <- pt_new_system()
    while (ok) {
      K <- kkt_triplets(ev$H, sig, pt$dx, J, pt$dc)
      d <- tryCatch(ksolve(K, rhs), error = function(err) NULL)
      if (!is.null(d) && !all(is.finite(d))) d <- NULL
      singular <- is.null(d)
      wrong <- FALSE
      if (!singular) {
        if (dim <= .nlp_inertia_max_dim) {
          # the inertia of the matrix, as MUMPS reports it to Ipopt
          ine <- .nlp_inertia(K$i, K$j, K$v, dim)
          if (ine$zero > 0 || ine$neg < m) singular <- TRUE else if (ine$neg > m) wrong <- TRUE
        } else {
          # above the dense limit, the curvature test of the inertia free
          # regularization stands in for the inertia
          dy <- d[seq_len(ny)]
          hdy <- .nlp_accum(ev$H$i, ev$H$v * dy[ev$H$j], ny)
          curv <- sum(dy * hdy) + sum(ifelse(free, sig + pt$dx, 0) * dy^2)
          wrong <- curv < 1e-10 * sum(dy^2)
        }
      }
      if (!singular && !wrong) break
      d <- NULL
      ok <- if (singular && m > 0) pt_singular() else pt_wrong_inertia()
    }
    list(d = d, K = K, dw = pt$dx)
  }
  # AdaptiveMuUpdate: a free mode that takes the barrier parameter from an
  # oracle (quality function, Mehrotra's probing or LOQO's rule) and a fixed
  # mode with the monotone decrease, switched on the progress of the iterates
  adaptive <- identical(o$mu_strategy, "adaptive")
  free_mode <- TRUE
  ad_filter <- matrix(numeric(), 0L, 2L)
  ad_refs <- numeric()
  ad_accepted <- NULL
  ad_kkt <- NULL
  ad_mu_max <- o$mu_max
  ad_mu_min_default <- is.na(o$mu_min)
  ad_mu_min <- if (ad_mu_min_default) (if (inner) 1e-9 else 1e-11) else o$mu_min
  ad_init_inf <- c(-1, -1)
  if (adaptive && !inner) mu <- 1
  compl_vec <- function(yy, zl, zu) {
    g <- gaps(yy)
    c((g$l * zl)[hl], (g$u * zu)[hu])
  }
  ad_norm <- function(v, n) {
    n <- max(n, 1)
    switch(o$quality_function_norm_type, "1-norm" = sum(abs(v)) / n, "2-norm-squared" = sum(v^2) / n,
           "max-norm" = max(abs(v), 0), "2-norm" = sqrt(sum(v^2)) / sqrt(n))
  }
  ad_extra <- function(du, pr, cm, xi) {
    cen <- switch(o$quality_function_centrality, none = 0, log = -cm * log(xi), reciprocal = cm / xi,
                  "cubed-reciprocal" = cm / xi^3)
    bal <- if (identical(o$quality_function_balancing_term, "cubic")) max(0, max(du, pr) - cm)^3 else 0
    cen + bal
  }
  ad_kkt_error <- function() {
    cv <- compl_vec(y, zl, zu)
    du <- ad_norm(q$dual[free], sum(free)); pr <- ad_norm(q$c, m); cm <- ad_norm(cv, nbnd)
    xi <- if (length(cv)) min(cv) / mean(cv) else 0
    du + pr + cm + ad_extra(du, pr, cm, xi)
  }
  ad_safeguard <- function() {
    if (o$adaptive_mu_safeguard_factor == 0) return(0)
    du <- sum(abs(q$dual[free])) / sum(free)
    pr <- sum(abs(q$c)) / max(m, 1)
    if (ad_init_inf[1L] < 0) ad_init_inf[1L] <<- max(1, du)
    if (ad_init_inf[2L] < 0) ad_init_inf[2L] <<- max(1, pr)
    v <- o$adaptive_mu_safeguard_factor * max(du / ad_init_inf[1L], pr / ad_init_inf[2L])
    if (identical(o$adaptive_mu_globalization, "kkt-error") && length(ad_refs)) v <- min(v, min(ad_refs))
    v
  }
  ad_progress <- function(nlp_error) {
    switch(o$adaptive_mu_globalization,
      "kkt-error" = length(ad_refs) < o$adaptive_mu_kkterror_red_iters ||
        any(ad_kkt_error() <= o$adaptive_mu_kkterror_red_fact * ad_refs),
      "obj-constr-filter" = {
        mg <- o$filter_margin_fact * min(o$filter_max_margin, nlp_error)
        all(ev$f + mg <= ad_filter[, 1L] | q$theta + mg <= ad_filter[, 2L])
      },
      TRUE)
  }
  ad_remember <- function() {
    if (identical(o$adaptive_mu_globalization, "kkt-error")) {
      if (length(ad_refs) >= o$adaptive_mu_kkterror_red_iters) ad_refs <<- ad_refs[-1L]
      ad_refs <<- c(ad_refs, ad_kkt_error())
    } else if (identical(o$adaptive_mu_globalization, "obj-constr-filter")) {
      v <- c(ev$f, q$theta)
      keep <- !(v[1L] <= ad_filter[, 1L] & v[2L] <= ad_filter[, 2L])
      ad_filter <<- rbind(ad_filter[keep, , drop = FALSE], v)
    }
    if (isTRUE(o$adaptive_mu_restore_previous_iterate)) ad_accepted <<- list(y = y, lam = lam, zl = zl, zu = zu)
  }
  # the barrier parameter of an oracle; NULL when the linear system fails
  ad_oracle <- function(mu_lo, mu_hi) {
    avg <- mean(compl_vec(y, zl, zu))
    if (identical(o$mu_oracle, "loqo")) {
      cv <- compl_vec(y, zl, zu)
      xi <- min(cv) / avg
      sigma <- 0.1 * min(0.05 * (1 - xi) / xi, 2)^3
      return(max(min(mu_hi, sigma * avg), mu_lo))
    }
    rs <- reg_solve(rhs_at(0))
    if (is.null(rs$d)) return(NULL)
    ad_kkt <<- rs
    g <- gaps(y)
    dir_of <- function(d, mm) {
      dy <- d[seq_len(ny)]
      list(dy = dy, dzl = ifelse(hl, mm / g$l - zl - zl / g$l * dy, 0),
           dzu = ifelse(hu, mm / g$u - zu + zu / g$u * dy, 0))
    }
    a0 <- dir_of(rs$d, 0)
    if (identical(o$mu_oracle, "probing")) {
      ap <- min(.nlp_ftb(g$l, a0$dy, 1), .nlp_ftb(g$u, -a0$dy, 1))
      ad <- min(.nlp_ftb(zl, a0$dzl, 1), .nlp_ftb(zu, a0$dzu, 1))
      mu_aff <- mean(c(((g$l + ap * a0$dy) * (zl + ad * a0$dzl))[hl],
                       ((g$u - ap * a0$dy) * (zu + ad * a0$dzu))[hu]))
      sigma <- min((mu_aff / avg)^3, o$sigma_max)
      return(max(min(sigma * avg, mu_hi), mu_lo))
    }
    d1 <- tryCatch(ksolve(rs$K, rhs_at(avg)), error = function(err) NULL)
    if (is.null(d1) || !all(is.finite(d1))) return(NULL)
    a1 <- dir_of(d1, avg)
    # the centering direction is the change of the step per unit of mu
    cen <- list(dy = (a1$dy - a0$dy) / avg, dzl = (a1$dzl - a0$dzl) / avg, dzu = (a1$dzu - a0$dzu) / avg)
    dnorm2 <- sum(q$dual[free]^2); pnorm2 <- sum(q$c^2)
    nd <- sum(free)
    qf <- function(sigma) {
      s_mu <- sigma * avg
      dy <- a0$dy + s_mu * cen$dy; dzl <- a0$dzl + s_mu * cen$dzl; dzu <- a0$dzu + s_mu * cen$dzu
      ap <- min(.nlp_ftb(g$l, dy, tau), .nlp_ftb(g$u, -dy, tau))
      ad <- min(.nlp_ftb(zl, dzl, tau), .nlp_ftb(zu, dzu, tau))
      cv <- c(((g$l + ap * dy) * (zl + ad * dzl))[hl], ((g$u - ap * dy) * (zu + ad * dzu))[hu])
      switch(o$quality_function_norm_type,
        "1-norm" = {
          du <- (1 - ad) * sum(abs(q$dual[free])) / nd
          pr <- (1 - ap) * sum(abs(q$c)) / max(m, 1); cm <- sum(abs(cv)) / nbnd
        },
        "2-norm-squared" = {
          du <- (1 - ad)^2 * dnorm2 / nd; pr <- (1 - ap)^2 * pnorm2 / max(m, 1); cm <- sum(cv^2) / nbnd
        },
        "max-norm" = {
          du <- (1 - ad) * max(abs(q$dual[free]), 0); pr <- (1 - ap) * max(abs(q$c), 0); cm <- max(abs(cv), 0)
        },
        "2-norm" = {
          du <- (1 - ad) * sqrt(dnorm2) / sqrt(nd); pr <- (1 - ap) * sqrt(pnorm2) / sqrt(max(m, 1))
          cm <- sqrt(sum(cv^2)) / sqrt(nbnd)
        })
      du + pr + cm + ad_extra(du, pr, cm, min(cv) / mean(cv))
    }
    golden <- function(up_in, q_up, lo_in, q_lo) {
      up <- log(up_in); lo <- log(lo_in)
      gf <- (3 - sqrt(5)) / 2
      m1 <- lo + gf * (up - lo); m2 <- lo + (1 - gf) * (up - lo)
      q1 <- qf(exp(m1)); q2 <- qf(exp(m2))
      rel <- function() 1 - min(q_lo, q_up, q1, q2) / max(q_lo, q_up, q1, q2)
      ns <- 0L
      while ((up - lo) >= o$quality_function_section_sigma_tol * up && rel() >= o$quality_function_section_qf_tol &&
             ns < o$quality_function_max_section_steps) {
        ns <- ns + 1L
        if (q1 > q2) {
          lo <- m1; q_lo <- q1; m1 <- m2; q1 <- q2
          m2 <- lo + (1 - gf) * (up - lo); q2 <- qf(exp(m2))
        } else {
          up <- m2; q_up <- q2; m2 <- m1; q2 <- q1
          m1 <- lo + gf * (up - lo); q1 <- qf(exp(m1))
        }
      }
      if ((up - lo) >= o$quality_function_section_sigma_tol * up && rel() < o$quality_function_section_qf_tol) {
        qmin <- min(q_lo, q_up, q1, q2)
        s <- if (qmin == q_lo) lo else if (qmin == q1) m1 else if (qmin == q2) m2 else up
      } else {
        if (q1 < q2) {
          s <- m1; qq <- q1
        } else {
          s <- m2; qq <- q2
        }
        if (up == log(up_in)) {
          qt <- if (q_up < 0) qf(exp(up)) else q_up
          if (qt < qq) s <- up
        } else if (lo == log(lo_in)) {
          qt <- if (q_lo < 0) qf(exp(lo)) else q_lo
          if (qt < qq) s <- lo
        }
      }
      exp(s)
    }
    q_1 <- qf(1)
    s1m <- 1 - max(1e-4, o$quality_function_section_sigma_tol)
    q_1m <- qf(s1m)
    if (q_1m > q_1) {
      up <- min(o$sigma_max, mu_hi / avg)
      sigma <- if (1 >= up) up else golden(up, -100, 1, q_1)
    } else {
      lo <- max(o$sigma_min, mu_lo / avg)
      up <- min(max(lo, s1m), mu_hi / avg)
      sigma <- if (lo >= up) lo else golden(up, q_1m, lo, -100)
    }
    sigma * avg
  }
  # one call of UpdateBarrierParameter; NA when a tiny step leaves the
  # barrier parameter unchanged (Ipopt's TINY_STEP_DETECTED)
  adaptive_update <- function(nlp_error) {
    if (ad_mu_min_default) ad_mu_min <<- min(ad_mu_min, 0.5 * min(o$tol, o$compl_inf_tol * sf))
    if (is.na(ad_mu_max)) ad_mu_max <<- o$mu_max_fact * (if (nbnd) mean(compl_vec(y, zl, zu)) else 0)
    if (!nbnd) {
      mu <<- ad_mu_min; tau <<- o$tau_min
      set_resto_mu(mu)
      return(TRUE)
    }
    tflag <- tiny_flag
    tiny_flag <<- FALSE
    reset_ls <- function() {
      soft <<- FALSE; soft_count <<- 0L
      filter_reset()
    }
    if (!free_mode) {
      if (ad_progress(nlp_error) && !tflag) {
        free_mode <<- TRUE
        ad_remember()
      } else {
        qm <- measures(ev, y, lam, zl, zu, mu)
        if (max(qm$inf_du / qm$sd, qm$inf_pr, qm$compl / qm$sc) <= o$barrier_tol_factor * mu || tflag) {
          new_mu <- min(o$mu_linear_decrease_factor * mu, mu^o$mu_superlinear_decrease_power)
          new_mu <- max(new_mu, min(o$compl_inf_tol * sf, o$tol) / (o$barrier_tol_factor + 1))
          if (tflag && new_mu == mu) return(NA)
          tau <<- max(o$tau_min, 1 - mu)
          mu <<- new_mu
          reset_ls()
        }
      }
    } else {
      prog <- ad_progress(nlp_error)
      if (!identical(o$adaptive_mu_globalization, "never-monotone-mode") && tflag) prog <- FALSE
      if (prog) {
        ad_remember()
      } else {
        free_mode <<- FALSE
        if (!is.null(ad_accepted)) {
          y <<- ad_accepted$y; lam <<- ad_accepted$lam; zl <<- ad_accepted$zl; zu <<- ad_accepted$zu
          ev <<- sev(y[ix], if (lm) 1L else 2L, lam[seq_len(nE)], lam[nE + seq_len(nI)])
          if (lm) ev$H <<- lm_update(ev, y, lam)
          q <<- measures(ev, y, lam, zl, zu)
        }
        new_mu <- o$adaptive_mu_monotone_init_factor * mean(compl_vec(y, zl, zu))
        new_mu <- max(new_mu, ad_safeguard())
        new_mu <- min(max(min(new_mu, 1e19), ad_mu_min), ad_mu_max)
        if (tflag && new_mu == mu) return(NA)
        tau <<- max(o$tau_min, 1 - new_mu)
        mu <<- new_mu
        reset_ls()
      }
    }
    if (free_mode) {
      tau <<- max(o$tau_min, 1 - nlp_error)
      new_mu <- ad_oracle(max(ad_mu_min, o$mu_target), ad_mu_max)
      if (!is.null(new_mu)) {
        new_mu <- max(new_mu, ad_mu_min)
        sg <- ad_safeguard()
        if (new_mu < sg) new_mu <- sg
        mu <<- min(new_mu, ad_mu_max)
      }
      reset_ls()
    }
    set_resto_mu(mu)
    TRUE
  }
  restore_acceptable <- function() {
    if (is.null(acceptable_point)) return(FALSE)
    y <<- acceptable_point$y; lam <<- acceptable_point$lam
    zl <<- acceptable_point$zl; zu <<- acceptable_point$zu
    TRUE
  }
  repeat {
    set_resto_mu(mu)
    ev <- tryCatch(sev(y[ix], if (lm) 1L else 2L, lam[seq_len(nE)], lam[nE + seq_len(nI)]),
                   error = function(err) NULL)
    if (is.null(ev) || !all(is.finite(c(ev$f, ev$gf, ev$H$v, ev$JE$v, ev$JI$v)))) {
      status <- "InvalidNumber"
      break
    }
    if (lm) ev$H <- lm_update(ev, y, lam)
    q <- measures(ev, y, lam, zl, zu)
    uq <- unscaled(ev, q, y)
    # convergence (OptimalityErrorConvergenceCheck and, in the
    # restoration phase, RestoConvergenceCheck); Ipopt shows an iterate
    # after its convergence check, so a restoration phase that ends in
    # local infeasibility (an exception) does not show its last iterate
    if (!inner) out_line(q, uq, ev)
    if (inner) {
      if (it >= o$max_iter) { out_line(q, uq, ev); status <- "MaxIters"; break }
      rs <- resto$check(y, first_iter)
      if (identical(rs, "converged")) { out_line(q, uq, ev); status <- "Success"; break }
    }
    overall <- max(q$inf_du / q$sd, q$inf_pr, q$compl / q$sc)
    conv <- "continue"
    if (is.finite(overall) && overall <= o$tol && uq$dual <= o$dual_inf_tol && uq$viol <= o$constr_viol_tol &&
        uq$compl <= o$compl_inf_tol) {
      conv <- "converged"
    } else if (o$acceptable_iter > 0 && is.finite(overall) && overall <= o$acceptable_tol &&
               uq$dual <= o$acceptable_dual_inf_tol && uq$viol <= o$acceptable_constr_viol_tol &&
               uq$compl <= o$acceptable_compl_inf_tol &&
               (is.na(last_f) || abs(ev$f - last_f) / max(1, abs(ev$f)) <= o$acceptable_obj_change_tol)) {
      acc_count <- acc_count + 1L
      if (acc_count >= o$acceptable_iter) conv <- "acceptable"
    } else {
      acc_count <- 0L
    }
    last_f <- ev$f
    if (inner && conv != "continue") {
      orig_inf <- resto$orig_inf(y)
      if (orig_inf <= 1e2 * o$tol && o$tol > 1e-1 * resto$orig_tol) {
        o$tol <- 1e-2 * o$tol
        conv <- "continue"
      } else {
        status <- if (orig_inf <= 1e2 * o$tol) "FeasiblePointFound" else "Infeasible"
        break
      }
    }
    if (inner) out_line(q, uq, ev)
    if (conv == "converged") { status <- "Success"; break }
    if (conv == "acceptable") { status <- "Acceptable"; break }
    if (it >= o$max_iter) { status <- "MaxIters"; break }
    if (max(abs(y[ix][!fixed]), 0) > o$diverging_iterates_tol) { status <- "Diverging"; break }
    if (proc.time()[["elapsed"]] - start_time >= o$max_wall_time) { status <- "MaxTime"; break }
    # barrier parameter (MonotoneMuUpdate); no update in the first
    # iteration of a restoration phase
    if (adaptive) {
      if (is.na(adaptive_update(overall))) {
        status <- "TinyStep"
        break
      }
    } else if (!(inner && first_iter)) {
      berr <- function() {
        qm <- measures(ev, y, lam, zl, zu, mu)
        max(qm$inf_du / qm$sd, qm$inf_pr, qm$compl / qm$sc)
      }
      sub_err <- berr()
      done <- FALSE
      while ((sub_err <= o$barrier_tol_factor * mu || tiny_flag) && !done) {
        new_mu <- min(o$mu_linear_decrease_factor * mu, mu^o$mu_superlinear_decrease_power)
        new_mu <- max(new_mu, o$mu_target, min(o$tol, o$compl_inf_tol * sf) / (o$barrier_tol_factor + 1))
        changed <- new_mu != mu
        if (!changed && tiny_flag) {
          status <- "TinyStep"
          done <- NA
          break
        }
        mu <- new_mu
        tau <- max(o$tau_min, 1 - mu)
        set_resto_mu(mu)
        if (mu_initialized && !o$mu_allow_fast_monotone_decrease) {
          done <- TRUE
        } else if (!changed) {
          done <- TRUE
        } else {
          sub_err <- berr()
          done <- sub_err > o$barrier_tol_factor * mu
        }
        if (done && changed) {
          soft <- FALSE; soft_count <- 0L
          filter_reset()
        }
        tiny_flag <- FALSE
      }
      if (is.na(done)) break
      mu_initialized <- TRUE
    }
    first_iter <- FALSE
    # the search direction with Ipopt's inertia free regularization; a free
    # mode update of the adaptive barrier parameter has already regularized
    # the matrix of this iteration
    gl <- c(ev$gf, numeric(nI))
    JT <- jac_t(ev, lam)
    g <- gaps(y)
    gphi <- gl - ifelse(hl, mu / g$l, 0) + ifelse(hu, mu / g$u, 0) +
      o$kappa_d * mu * (as.numeric(lonly) - as.numeric(uonly))
    gphi[!free] <- 0
    c0 <- q$c
    rhs <- rhs_at(mu)
    fallback <- FALSE
    if (!is.null(ad_kkt)) {
      rs <- ad_kkt
      ad_kkt <- NULL
      rs$d <- tryCatch(ksolve(rs$K, rhs), error = function(err) NULL)
      if (!is.null(rs$d) && !all(is.finite(rs$d))) rs$d <- NULL
    } else {
      rs <- reg_solve(rhs)
    }
    d <- rs$d; dw <- rs$dw; K <- rs$K
    info$dw <- dw
    if (is.null(d)) fallback <- TRUE
    if (!fallback) st <- step_of(d, y, zl, zu)
    # FindAcceptableTrialPoint of the backtracking line search
    if (mu != last_mu) {
      in_wd <- FALSE; wd <- NULL; wd_shortened <- 0L
      last_mu <- mu
    }
    phi0 <- barrier_of(ev$f, y, mu)
    th0 <- q$theta
    is_acceptable_now <- conv == "continue" && acc_count > 0
    if (is_acceptable_now) acceptable_point <- list(y = y, lam = lam, zl = zl, zu = zu)
    goto_resto <- fallback || start_resto ||
      (expect_inf && m && max(abs(lam)) > o$expect_infeasible_problem_ytol)
    start_resto <- FALSE
    ref <- if (in_wd) wd$ref else list(theta = th0, barr = phi0, gbd = if (fallback) 0 else sum(gphi * st$dy))
    tiny <- !goto_resto && !fallback &&
      max(abs(st$dy[free]) / (1 + abs(y[free])), 0) <= o$tiny_step_tol && th0 <= 1e-4
    if (in_wd && (goto_resto || tiny)) {
      y <- wd$y; lam <- wd$lam; zl <- wd$zl; zu <- wd$zu; st <- wd$st
      ref <- wd$ref
      in_wd <- FALSE; wd <- NULL; wd_shortened <- 0L
      goto_resto <- FALSE; tiny <- FALSE
      ev <- sev(y[ix], 1L); g <- gaps(y)
      q <- measures(ev, y, lam, zl, zu); th0 <- q$theta; c0 <- q$c; phi0 <- ref$barr
    }
    if (o$watchdog_shortened_iter_trigger > 0 && !in_wd && !goto_resto && !tiny && !soft && !expect_inf &&
        wd_shortened >= o$watchdog_shortened_iter_trigger) {
      in_wd <- TRUE
      wd <- list(y = y, lam = lam, zl = zl, zu = zu, st = st, ref = ref, alpha_test = primal_ftb(y, st$dy))
      wd_trial <- 0L
    }
    accept <- FALSE
    n_steps <- 0L
    alpha <- 0
    tr <- NULL
    dual_done <- FALSE
    try_soft <- function() {
      if (o$soft_resto_pderror_reduction_factor == 0) return(list(ok = FALSE))
      a <- min(primal_ftb(y, st$dy), dual_ftb(st, zl, zu))
      tt <- trial_at(y + a * st$dy, 1L)
      if (is.null(tt$e)) return(list(ok = FALSE))
      lt <- lam + a * st$dlam; zlt <- zl + a * st$dzl; zut <- zu + a * st$dzu
      res <- list(ok = TRUE, orig = FALSE, a = a, tr = tt, lam = lt, zl = zlt, zu = zut)
      if (check_trial(tt, 0)) {
        res$orig <- TRUE
        return(res)
      }
      if (pd_error(tt$e, tt$y, lt, zlt, zut) <= o$soft_resto_pderror_reduction_factor *
          pd_error(ev, y, lam, zl, zu)) return(res)
      list(ok = FALSE)
    }
    if (tiny) {
      alpha <- primal_ftb(y, st$dy)
      tr <- trial_at(y + alpha * st$dy)
      if (tiny_last) tiny_flag <- TRUE
      tiny_last <- max(abs(st$dlam), 0) < o$tiny_step_y_tol
      accept <- TRUE
      info$char <- if (tiny_last) "T" else "t"
    } else {
      tiny_last <- FALSE
    }
    if (!goto_resto && !tiny) {
      if (soft) {
        soft_count <- soft_count + 1L
        if (soft_count <= o$max_soft_resto_iters) {
          sr <- try_soft()
          if (sr$ok) {
            accept <- TRUE; dual_done <- TRUE
            alpha <- sr$a; tr <- sr$tr; lam <- sr$lam; zl <- sr$zl; zu <- sr$zu
            info$alpha_du <- sr$a
            info$char <- "s"
            if (sr$orig) {
              soft <- FALSE; soft_count <- 0L
              info$char <- "S"
            }
          }
        }
      } else {
        skip_first <- FALSE
        repeat {
          # DoBacktrackingLineSearch
          a_max <- primal_ftb(y, st$dy)
          a_min <- if (in_wd) a_max else {
            gbd <- sum(gphi * st$dy)
            am <- o$gamma_theta
            if (gbd < 0) {
              am <- min(o$gamma_theta, o$gamma_phi * th0 / (-gbd))
              if (theta_min >= 0 && th0 <= theta_min) am <- min(am, o$delta * th0^o$s_theta / (-gbd)^o$s_phi)
            }
            o$alpha_min_frac * am
          }
          alpha <- a_max
          a_test <- if (in_wd) wd$alpha_test else alpha
          if (skip_first) alpha <- alpha * o$alpha_red_factor
          eval_err <- FALSE
          soc_taken <- FALSE
          while (alpha > a_min || n_steps == 0L) {
            tr <- trial_at(y + alpha * st$dy)
            eval_err <- is.null(tr$e)
            a_test <- alpha
            if (!eval_err && (isTRUE(o$accept_every_trial_step) ||
                              (o$accept_after_max_steps != -1 && n_steps >= o$accept_after_max_steps))) {
              filter_reset()
              accept <- TRUE
            } else {
              accept <- !eval_err && check_trial(tr, a_test)
            }
            if (accept || in_wd) break
            if (expect_inf && count_shortened >= 5L) break
            if (!eval_err && alpha == a_max && th0 <= tr$theta && o$max_soc > 0) {
              # second order corrections
              csoc <- c0
              a_soc <- alpha
              th_trial <- tr$theta
              trs <- tr
              k <- 0L
              th_old <- 0
              while (k < o$max_soc && !accept && (k == 0L || th_trial <= o$kappa_soc * th_old)) {
                th_old <- th_trial
                csoc <- a_soc * csoc + trs$c
                ds <- tryCatch(ksolve(K, -c(ifelse(free, gphi + JT, 0), csoc)), error = function(err) NULL)
                if (is.null(ds) || !all(is.finite(ds))) break
                ss <- step_of(ds, y, zl, zu)
                a_soc <- primal_ftb(y, ss$dy)
                trs <- trial_at(y + a_soc * ss$dy)
                if (is.null(trs$e)) break
                accept <- check_trial(trs, a_test)
                if (accept) {
                  alpha <- a_soc; st <- ss; tr <- trs
                } else {
                  k <- k + 1L
                  th_trial <- trs$theta
                }
              }
              if (accept) {
                soc_taken <- TRUE
                break
              }
            }
            alpha <- alpha * o$alpha_red_factor
            n_steps <- n_steps + 1L
          }
          if (accept) {
            info$char <- if (!is_ftype(a_test) || !armijo(tr$barr, a_test)) {
              add_filter(ref$theta, ref$barr)
              "h"
            } else "f"
            if (soc_taken) info$char <- toupper(info$char)
          } else if (in_wd) {
            info$char <- "w"
          }
          if (in_wd) {
            if (accept) {
              in_wd <- FALSE
              wd <- NULL
              break
            }
            wd_trial <- wd_trial + 1L
            if (eval_err || wd_trial > o$watchdog_trial_iter_max) {
              # StopWatchDog: back to the reference point and its direction
              y <- wd$y; lam <- wd$lam; zl <- wd$zl; zu <- wd$zu; st <- wd$st
              ref <- wd$ref
              in_wd <- FALSE; wd <- NULL; wd_shortened <- 0L
              ev <- sev(y[ix], 1L); g <- gaps(y)
              q <- measures(ev, y, lam, zl, zu); th0 <- q$theta; c0 <- q$c; phi0 <- ref$barr
              JT <- jac_t(ev, lam)
              gl <- c(ev$gf, numeric(nI))
              gphi <- gl - ifelse(hl, mu / g$l, 0) + ifelse(hu, mu / g$u, 0) +
                o$kappa_d * mu * (as.numeric(lonly) - as.numeric(uonly))
              gphi[!free] <- 0
              skip_first <- TRUE
              next
            }
            accept <- TRUE
            break
          }
          break
        }
      }
    }
    info$ls <- n_steps + 1L
    if (!accept) {
      if (!soft && !goto_resto && !expect_inf) {
        add_filter(ref$theta, ref$barr)
        sr <- try_soft()
        if (sr$ok) {
          accept <- TRUE; dual_done <- TRUE
          alpha <- sr$a; tr <- sr$tr; lam <- sr$lam; zl <- sr$zl; zu <- sr$zu
          info$alpha_du <- sr$a
          if (sr$orig) info$char <- "S" else {
            soft <- TRUE
            info$char <- "s"
          }
        }
      }
      if (!accept) {
        if (!soft) add_filter(ref$theta, ref$barr)
        if (is_acceptable_now) {
          status <- "Acceptable"
          break
        }
        if (th0 <= 1e-2 * o$tol && uq$viol <= 1e-1 * o$constr_viol_tol) {
          status <- if (restore_acceptable()) "Acceptable" else "StepFailed"
          break
        }
        info$char <- "R"
        # the information of the failed line search: its step after the last
        # reduction and no dual step
        info$alpha_pr <- if (exists("alpha", inherits = FALSE) && !fallback) alpha else 0
        info$alpha_du <- 0
        if (inner) {
          # the restoration phase of the restoration phase resets p and n
          tr <- list(y = resto$reset_pn(y, mu))
          lam_keep <- lam
          accept <- TRUE
          dual_done <- TRUE
        } else {
          n_resto <- n_resto + 1L
          rr <- .nlp_restoration(environment())
          if (!identical(rr$status, "Success")) {
            if (rr$status %in% c("Infeasible", "MaxIters", "MaxTime")) {
              status <- rr$status
              if (!is.null(rr$y)) y <- rr$y
            } else if (isTRUE(rr$direct)) {
              status <- "RestorationFailed"
              if (!is.null(rr$y)) y <- rr$y
            } else {
              status <- if (restore_acceptable()) "Acceptable" else "RestorationFailed"
            }
            if (!is.null(rr$it)) it <- rr$it
            break
          }
          tr <- list(y = rr$y)
          zl <- rr$zl; zu <- rr$zu; lam <- rr$lam
          it <- rr$it - 1L
          accept <- TRUE
          dual_done <- TRUE
          # the iterate the restoration phase returns was already shown
          # (Set_info_skip_output)
          skip_out <- TRUE
        }
        count_shortened <- 0L
        expect_inf <- FALSE
        soft <- FALSE; soft_count <- 0L
        wd_shortened <- 0L
      }
    } else if (!dual_done) {
      ad <- dual_ftb(st, zl, zu)
      zl <- zl + ad * st$dzl; zu <- zu + ad * st$dzu
      lam <- lam + alpha * st$dlam
      info$alpha_du <- ad
      if (n_steps == 0L) {
        count_shortened <- 0L
        wd_shortened <- 0L
      } else {
        count_shortened <- count_shortened + 1L
        wd_shortened <- wd_shortened + 1L
      }
    }
    info$alpha_pr <- alpha
    info$dnorm <- if (exists("st", inherits = FALSE) && !fallback) max(abs(st$dy), 0) else 0
    y <- tr$y
    zc <- clip_z(y, zl, zu)
    zl <- zc$l; zu <- zc$u
    if (expect_inf && !is.null(tr$theta) && tr$theta <= o$expect_infeasible_problem_ctol) expect_inf <- FALSE
    it <- it + 1L
  }
  z <- y[ix]
  # the objective reported is that of the final iterate, before the
  # projection onto the original bounds (as Ipopt computes obj_value)
  f_final <- if (inner) NA_real_ else evalf(z, 0L)$f
  if (!inner && isTRUE(o$honor_original_bounds)) z <- pmin(pmax(z, lo_orig), hi_orig)
  if (inner) {
    return(list(y = y, lam = lam, zl = zl, zu = zu, mu = mu, it = it, status = status))
  }
  ev <- evalf(z, 0L)
  viol <- max(0, abs(ev$cE), ev$cI)
  list(z = z, f = f_final, lamE = lam[seq_len(nE)] * sE / sf, lamI = lam[nE + seq_len(nI)] * sI / sf,
       zl = zl[ix] / sf, zu = zu[ix] / sf, violation = viol, iterations = it, status = status,
       restorations = n_resto)
}

# Ipopt's restoration phase (MinC_1NrmRestorationPhase): minimize
#
#   rho sum(p + n) + eta(mu) / 2 ||D_R (z - z_R)||^2
#   subject to  cE(z) - pE + nE = 0,  cI(z) + s - pI + nI = 0,  p, n >= 0
#
# with the interior point method itself, eta(mu) = sqrt(mu) and
# D_R = min(1, 1 / |z_R|), until the original filter accepts the iterate and
# its infeasibility fell by required_infeasibility_reduction. `env` is the
# environment of the outer call at the failed line search.
.nlp_restoration <- function(env) {
  o <- env$o
  nz <- env$nz; nE <- env$nE; nI <- env$nI; m <- env$m; nfix <- env$fixed
  y <- env$y; mu <- env$mu
  ix <- seq_len(nz); isl <- nz + seq_len(nI)
  rho <- o$resto_penalty_parameter
  zr <- y[ix]
  DR2 <- pmin(1, 1 / abs(zr))^2
  DR2[nfix] <- 0
  ev <- env$sev(zr, 0L)
  cvec <- c(ev$cE, ev$cI + y[isl])
  mu0 <- max(mu, abs(cvec), 0)
  quad_n <- function(cv, mu) {
    a <- mu / (2 * rho) - cv / 2
    b <- cv * mu / (2 * rho)
    sqrt(a * a + b) + a
  }
  nn <- quad_n(cvec, mu0)
  pp <- cvec + nn
  rv <- new.env(parent = emptyenv())
  rv$mu <- mu0
  pE <- nz + seq_len(nE); nEi <- nz + nE + seq_len(nE)
  pI <- nz + 2L * nE + seq_len(nI); nIi <- nz + 2L * nE + nI + seq_len(nI)
  nzR <- nz + 2L * m
  evalfR <- function(zR, level, lamE = NULL, lamI = NULL, obj = 1) {
    z <- zR[ix]
    e <- env$sev(z, level, lamE, lamI, obj = 0)
    eta <- o$resto_proximity_weight * sqrt(rv$mu)
    dz <- z - zr
    out <- list(f = rho * sum(zR[-ix]) + eta / 2 * sum(DR2 * dz^2),
                cE = e$cE - zR[pE] + zR[nEi], cI = e$cI - zR[pI] + zR[nIi])
    if (level >= 1L) {
      out$gf <- c(eta * DR2 * dz, rep(rho, 2L * m))
      out$JE <- list(i = c(e$JE$i, seq_len(nE), seq_len(nE)), j = c(e$JE$j, pE, nEi),
                     v = c(e$JE$v, rep(-1, nE), rep(1, nE)))
      out$JI <- list(i = c(e$JI$i, seq_len(nI), seq_len(nI)), j = c(e$JI$j, pI, nIi),
                     v = c(e$JI$v, rep(-1, nI), rep(1, nI)))
    }
    if (level >= 2L) {
      out$H <- list(i = c(e$H$i, ix), j = c(e$H$j, ix), v = c(e$H$v, obj * eta * DR2))
    }
    out
  }
  ylo <- env$ylo; yhi <- env$yhi
  loR <- c(ylo[ix], rep(0, 2L * m)); hiR <- c(yhi[ix], rep(Inf, 2L * m))
  hl <- env$hl; hu <- env$hu
  # the layout of the restoration unknowns is (z, pE, nE, pI, nI)
  iE <- seq_len(nE); iI <- nE + seq_len(nI)
  pn <- function(p, n) c(p[iE], n[iE], p[iI], n[iI])
  zl0 <- c(pmin(rho, env$zl[ix]), pn(mu0 / pp, mu0 / nn), pmin(rho, env$zl[isl]))
  zu0 <- c(pmin(rho, env$zu[ix]), rep(0, 2L * m), env$zu[isl])
  y0 <- c(zr, pn(pp, nn), y[isl])
  # the original problem's view of a restoration iterate
  ref <- env$ref
  orig_tol <- o$tol
  orig_inf_start <- max(abs(cvec), 0)
  kappa <- o$required_infeasibility_reduction
  square <- nE == env$nfree && nI == 0L
  if (square) kappa <- 0
  if (!square && env$expect_inf && env$n_resto == 0L && orig_inf_start > 1e-3) kappa <- 1e-3
  orig_at <- function(yR) {
    z <- yR[ix]; s <- yR[nzR + seq_len(nI)]
    e <- tryCatch(env$sev(z, 0L), error = function(err) NULL)
    if (is.null(e)) return(list(theta = Inf, inf = Inf, barr = Inf))
    cc <- c(e$cE, e$cI + s)
    yy <- c(z, s)
    list(theta = sum(abs(cc)), inf = max(abs(cc), 0), barr = env$barrier_of(e$f, yy, mu))
  }
  orig_print <- function(yR) {
    e <- env$sev(yR[ix], 0L)
    list(f = e$f / env$sf,
         viol = max(abs(e$cE / env$sE), pmax((e$cI + env$slack_lo) / env$sI, 0), 0, na.rm = TRUE))
  }
  successive <- 0L
  check <- function(yR, first) {
    successive <<- successive + 1L
    if (successive > o$max_resto_iter) return("continue")
    if (first) return("continue")
    oa <- orig_at(yR)
    inf_max <- if (kappa == 0) 0 else max(kappa * orig_inf_start, min(orig_tol, o$constr_viol_tol))
    if (square && oa$inf <= min(orig_tol, o$constr_viol_tol)) return("converged")
    if (oa$inf > inf_max) return("continue")
    if (!env$filter_ok(oa$theta, oa$barr)) return("continue")
    if (!env$ok_to_iterate(oa$barr, oa$theta, TRUE)) return("continue")
    "converged"
  }
  reset_pn <- function(yR, muR) {
    z <- yR[ix]; s <- yR[nzR + seq_len(nI)]
    e <- env$sev(z, 0L)
    cv <- c(e$cE, e$cI + s)
    n2 <- quad_n(cv, muR)
    yR[nz + seq_len(2L * m)] <- pn(cv + n2, n2)
    yR
  }
  ro <- o
  ro$expect_infeasible_problem <- FALSE
  ro$start_with_resto <- FALSE
  ro$constr_mult_init_max <- o$resto_constr_mult_init_max %||% 0
  keyE <- env$key[env$ny + seq_len(nE)]; keyI <- env$key[nz + seq_len(nI)]
  r <- .nlp_interior_point(evalfR, y0[seq_len(nzR)], loR, hiR,
                           c(env$key[ix], keyE, keyE, keyI, keyI), keyE, keyI, tol = o$tol,
                           maxiter = o$max_iter, verbose = env$verbose, opts = ro,
                           resto = list(y0 = y0, zl0 = zl0, zu0 = zu0, mu0 = mu0, it0 = env$it + 1L,
                                        fixed = c(nfix, rep(FALSE, 2L * m)), slack_lo = ylo[isl], env = rv,
                                        check = check, orig_inf = function(yR) orig_at(yR)$inf,
                                        orig_tol = orig_tol, reset_pn = reset_pn, lm_idx = which(!nfix),
                                        info0 = utils::modifyList(env$info, list(dnorm = 0)),
                                        orig_print = orig_print))
  yR <- r$y
  ynew <- c(yR[ix], yR[nzR + seq_len(nI)])
  if (!identical(r$status, "Success")) {
    st <- switch(r$status, Infeasible = "Infeasible", MaxIters = "MaxIters", MaxTime = "MaxTime",
                 "RestorationFailed")
    # a restoration phase that stops at a tiny step or an acceptable point
    # declares local infeasibility unless the point is nearly feasible
    direct <- FALSE
    if (r$status %in% c("TinyStep", "Acceptable")) {
      thr <- o$resto_failure_feasibility_threshold %||% (1e2 * orig_tol)
      st <- if (orig_at(yR)$inf <= thr) "RestorationFailed" else "Infeasible"
      direct <- TRUE
    }
    return(list(status = st, y = ynew, it = r$it, direct = direct))
  }
  # bound multipliers from a Newton step for complementarity over the
  # whole restoration phase, reset to 1 when too large
  g0 <- env$gaps(y); g1 <- env$gaps(ynew)
  dzl <- ifelse(hl, ((g0$l - g1$l) * env$zl + mu) / g0$l - env$zl, 0)
  dzu <- ifelse(hu, ((g0$u - g1$u) * env$zu + mu) / g0$u - env$zu, 0)
  tau <- env$tau
  ad <- min(.nlp_ftb(env$zl, dzl, tau), .nlp_ftb(env$zu, dzu, tau))
  zl <- env$zl + ad * dzl; zu <- env$zu + ad * dzu
  if (max(abs(zl), abs(zu), 0) > o$bound_mult_reset_threshold) {
    zl <- ifelse(hl, 1, 0); zu <- ifelse(hu, 1, 0)
  }
  lam <- numeric(m)
  if (o$constr_mult_reset_threshold > 0) {
    e1 <- env$sev(ynew[ix], 1L)
    lam <- env$ls_mult(e1, ynew, zl, zu, o$constr_mult_reset_threshold)
  }
  list(status = "Success", y = ynew, zl = zl, zu = zu, lam = lam, it = r$it)
}
