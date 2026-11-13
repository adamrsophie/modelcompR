# scipy.optimize.least_squares (scipy/optimize/_lsq): the trust region
# reflective method "trf" and the dogleg method in a rectangular trust region
# "dogbox" with the exact (SVD and lstsq) subproblem solvers, the robust loss
# functions, and the option mapping of OptimizationSciPy.jl. Method "lm"
# runs MINPACK (R/optimization-scipy-minpack.R).

.sp_nrm <- function(x) sqrt(sum(x * x))

.sp_lsq_loss <- function(loss, f_scale) {
  if (loss == "linear") return(NULL)
  rhofun <- switch(loss,
    huber = function(z) {
      mask <- z <= 1
      list(ifelse(mask, z, 2 * z^0.5 - 1), ifelse(mask, 1, z^-0.5), ifelse(mask, 0, -0.5 * z^-1.5))
    },
    soft_l1 = function(z) {
      t <- 1 + z
      list(2 * (t^0.5 - 1), t^-0.5, -0.5 * t^-1.5)
    },
    cauchy = function(z) {
      t <- 1 + z
      list(log1p(z), 1 / t, -1 / t^2)
    },
    arctan = function(z) {
      t <- 1 + z^2
      list(atan(z), 1 / t, -2 * z / t^2)
    })
  function(f, cost_only = FALSE) {
    rho <- rhofun((f / f_scale)^2)
    if (cost_only) return(0.5 * f_scale^2 * sum(rho[[1L]]))
    rho[[1L]] <- rho[[1L]] * f_scale^2
    rho[[3L]] <- rho[[3L]] / f_scale^2
    rho
  }
}

# scale_for_robust_loss_function: J and f of the equivalent linear problem
.sp_lsq_robust <- function(J, f, loss_function) {
  rho <- loss_function(f)
  js <- rho[[2L]] + 2 * rho[[3L]] * f^2
  js[js < .spmp_epsmch] <- .spmp_epsmch
  js <- js^0.5
  list(J = J * js, f = f * (rho[[2L]] / js), cost = 0.5 * sum(rho[[1L]]))
}

# J is a matrix or a function computing J v (a LinearOperator)
.sp_lsq_mv <- function(J, s) if (is.function(J)) J(s) else as.numeric(J %*% s)

.sp_lsq_eval_quad <- function(J, g, s, diag = NULL) {
  Js <- .sp_lsq_mv(J, s)
  q <- sum(Js * Js)
  if (!is.null(diag)) q <- q + sum(s * diag * s)
  0.5 * q + sum(s * g)
}

.sp_lsq_quad_1d <- function(J, g, s, diag = NULL, s0 = NULL) {
  v <- .sp_lsq_mv(J, s)
  a <- sum(v * v)
  if (!is.null(diag)) a <- a + sum(s * diag * s)
  a <- a * 0.5
  b <- sum(g * s)
  if (is.null(s0)) return(list(a = a, b = b))
  u <- .sp_lsq_mv(J, s0)
  b <- b + sum(u * v)
  c <- 0.5 * sum(u * u) + sum(g * s0)
  if (!is.null(diag)) {
    b <- b + sum(s0 * diag * s)
    c <- c + 0.5 * sum(s0 * diag * s0)
  }
  list(a = a, b = b, c = c)
}

.sp_lsq_min_quad_1d <- function(a, b, lb, ub, c = 0) {
  t <- c(lb, ub)
  if (a != 0) {
    extremum <- -0.5 * b / a
    if (lb < extremum && extremum < ub) t <- c(t, extremum)
  }
  y <- t * (a * t + b) + c
  i <- which.min(y)
  list(t = t[i], y = y[i])
}

.sp_lsq_step_to_bound <- function(x, s, lb, ub) {
  nz <- s != 0
  steps <- rep(Inf, length(x))
  steps[nz] <- pmax((lb - x)[nz] / s[nz], (ub - x)[nz] / s[nz])
  mn <- min(steps)
  list(step = mn, hits = as.integer(steps == mn) * as.integer(sign(s)))
}

.sp_lsq_intersect_tr <- function(x, s, Delta) {
  a <- sum(s * s)
  if (a == 0) stop("`s` is zero.", call. = FALSE)
  b <- sum(x * s)
  c <- sum(x * x) - Delta^2
  if (c > 0) stop("`x` is not within the trust region.", call. = FALSE)
  d <- sqrt(b * b - a * c)
  q <- -(b + if (b < 0 || identical(1 / b, -Inf)) -d else d)
  t1 <- q / a
  t2 <- c / q
  if (t1 < t2) c(t1, t2) else c(t2, t1)
}

.sp_lsq_active <- function(x, lb, ub, rtol = 1e-10) {
  active <- integer(length(x))
  if (rtol == 0) {
    active[x <= lb] <- -1L
    active[x >= ub] <- 1L
    return(active)
  }
  lower_dist <- x - lb
  upper_dist <- ub - x
  lower_threshold <- rtol * pmax(1, abs(lb))
  upper_threshold <- rtol * pmax(1, abs(ub))
  active[is.finite(lb) & lower_dist <= pmin(upper_dist, lower_threshold)] <- -1L
  active[is.finite(ub) & upper_dist <= pmin(lower_dist, upper_threshold)] <- 1L
  active
}

# numpy.nextafter for normal numbers and zero
.sp_nextafter <- function(a, b) {
  vapply(seq_along(a), function(i) {
    x <- a[i]; y <- b[i]
    if (x == y) return(y)
    if (x == 0) return(if (y > 0) 4.9406564584124654e-324 else -4.9406564584124654e-324)
    ax <- abs(x)
    e <- floor(log2(ax))
    if (2^e > ax) e <- e - 1
    if (2^(e + 1) <= ax) e <- e + 1
    sp <- 2^(max(e, -1022) - 52)
    away <- (y > x) == (x > 0)
    if (!away && ax == 2^e && e > -1022) sp <- sp / 2
    sign(x) * (if (away) ax + sp else ax - sp)
  }, 0)
}

.sp_lsq_strictly_feasible <- function(x, lb, ub, rstep = 1e-10) {
  active <- .sp_lsq_active(x, lb, ub, rstep)
  lo <- active == -1L
  up <- active == 1L
  if (rstep == 0) {
    x[lo] <- .sp_nextafter(lb[lo], ub[lo])
    x[up] <- .sp_nextafter(ub[up], lb[up])
  } else {
    x[lo] <- lb[lo] + rstep * pmax(1, abs(lb[lo]))
    x[up] <- ub[up] - rstep * pmax(1, abs(ub[up]))
  }
  tight <- x < lb | x > ub
  x[tight] <- 0.5 * (lb[tight] + ub[tight])
  x
}

.sp_lsq_update_tr <- function(Delta, actual_reduction, predicted_reduction, step_norm, bound_hit) {
  ratio <- if (predicted_reduction > 0) actual_reduction / predicted_reduction
           else if (predicted_reduction == 0 && actual_reduction == 0) 1 else 0
  if (ratio < 0.25) Delta <- 0.25 * step_norm
  else if (ratio > 0.75 && bound_hit) Delta <- Delta * 2
  list(Delta = Delta, ratio = ratio)
}

.sp_lsq_check_termination <- function(dF, F, dx_norm, x_norm, ratio, ftol, xtol) {
  fs <- dF < ftol * F && ratio > 0.25
  xs <- dx_norm < xtol * (xtol + x_norm)
  if (fs && xs) 4L else if (fs) 2L else if (xs) 3L else NULL
}

# solve_lsq_trust_region: the Levenberg Marquardt parameter by Newton
# iterations on the secular equation, from the SVD of the Jacobian
.sp_lsq_trust_region <- function(n, m, uf, s, V, Delta, initial_alpha = NULL, rtol = 0.01, max_iter = 10L) {
  suf <- s * uf
  phid <- function(alpha) {
    denom <- s^2 + alpha
    p_norm <- .sp_nrm(suf / denom)
    c(p_norm - Delta, -sum(suf^2 / denom^3) / p_norm)
  }
  full_rank <- m >= n && s[length(s)] > .spmp_epsmch * m * s[1L]
  if (full_rank) {
    p <- -as.numeric(V %*% (uf / s))
    if (.sp_nrm(p) <= Delta) return(list(p = p, alpha = 0, it = 0L))
  }
  alpha_upper <- .sp_nrm(suf) / Delta
  alpha_lower <- if (full_rank) { pd <- phid(0); -pd[1L] / pd[2L] } else 0
  alpha <- if (is.null(initial_alpha) || (!full_rank && initial_alpha == 0))
    max(0.001 * alpha_upper, (alpha_lower * alpha_upper)^0.5) else initial_alpha
  for (it in seq_len(max_iter)) {
    if (alpha < alpha_lower || alpha > alpha_upper) alpha <- max(0.001 * alpha_upper, (alpha_lower * alpha_upper)^0.5)
    pd <- phid(alpha)
    if (pd[1L] < 0) alpha_upper <- alpha
    ratio <- pd[1L] / pd[2L]
    alpha_lower <- max(alpha_lower, alpha - ratio)
    alpha <- alpha - (pd[1L] + Delta) * ratio / Delta
    if (abs(pd[1L]) < rtol * Delta) break
  }
  p <- -as.numeric(V %*% (suf / (s^2 + alpha)))
  list(p = p * (Delta / .sp_nrm(p)), alpha = alpha, it = it)
}

# select_step of trf: the trust region step, its reflection at the bound and
# the scaled gradient step, whichever reduces the quadratic model most
.sp_lsq_select_step <- function(x, J_h, diag_h, g_h, p, p_h, d, Delta, lb, ub, theta) {
  if (all(x + p >= lb & x + p <= ub))
    return(list(step = p, step_h = p_h, pred = -.sp_lsq_eval_quad(J_h, g_h, p_h, diag_h)))
  sb <- .sp_lsq_step_to_bound(x, p, lb, ub)
  p_stride <- sb$step
  r_h <- p_h
  r_h[sb$hits != 0L] <- -r_h[sb$hits != 0L]
  r <- d * r_h
  p <- p * p_stride
  p_h <- p_h * p_stride
  to_tr <- .sp_lsq_intersect_tr(p_h, r_h, Delta)[2L]
  to_bound <- .sp_lsq_step_to_bound(x + p, r, lb, ub)$step
  r_stride <- min(to_bound, to_tr)
  if (r_stride > 0) {
    r_stride_l <- (1 - theta) * p_stride / r_stride
    r_stride_u <- if (r_stride == to_bound) theta * to_bound else to_tr
  } else {
    r_stride_l <- 0
    r_stride_u <- -1
  }
  if (r_stride_l <= r_stride_u) {
    q <- .sp_lsq_quad_1d(J_h, g_h, r_h, s0 = p_h, diag = diag_h)
    mq <- .sp_lsq_min_quad_1d(q$a, q$b, r_stride_l, r_stride_u, c = q$c)
    r_h <- r_h * mq$t + p_h
    r <- r_h * d
    r_value <- mq$y
  } else {
    r_value <- Inf
  }
  p <- p * theta
  p_h <- p_h * theta
  p_value <- .sp_lsq_eval_quad(J_h, g_h, p_h, diag_h)
  ag_h <- -g_h
  ag <- d * ag_h
  to_tr <- Delta / .sp_nrm(ag_h)
  to_bound <- .sp_lsq_step_to_bound(x, ag, lb, ub)$step
  ag_stride <- if (to_bound < to_tr) theta * to_bound else to_tr
  q <- .sp_lsq_quad_1d(J_h, g_h, ag_h, diag = diag_h)
  mq <- .sp_lsq_min_quad_1d(q$a, q$b, 0, ag_stride)
  ag_h <- ag_h * mq$t
  ag <- ag * mq$t
  if (p_value < r_value && p_value < mq$y) return(list(step = p, step_h = p_h, pred = -p_value))
  if (r_value < p_value && r_value < mq$y) return(list(step = r, step_h = r_h, pred = -r_value))
  list(step = ag, step_h = ag_h, pred = -mq$y)
}

.sp_lsq_scale <- function(J, x_scale, scale_inv_old = NULL) {
  if (!identical(x_scale, "jac")) return(list(scale = x_scale, inv = 1 / x_scale))
  inv <- colSums(J^2)^0.5
  if (is.null(scale_inv_old)) inv[inv == 0] <- 1 else inv <- pmax(inv, scale_inv_old)
  list(scale = 1 / inv, inv = inv)
}

# trf: trust region reflective (Branch, Coleman and Li 1999) with the
# Coleman and Li scaling when there are bounds
.sp_lsq_trf <- function(fun, jac, x0, f0, J0, lb, ub, ftol, xtol, gtol, max_nfev, x_scale, loss_function,
                        tr_solver = "exact", tr_options = list()) {
  bounded <- !(all(lb == -Inf) && all(ub == Inf))
  x <- x0; f <- f0; f_true <- f
  nfev <- 1L; njev <- 1L
  J <- J0
  m <- nrow(J); n <- ncol(J)
  if (!is.null(loss_function)) {
    rb <- .sp_lsq_robust(J, f, loss_function); J <- rb$J; f <- rb$f; cost <- rb$cost
  } else {
    cost <- 0.5 * sum(f * f)
  }
  g <- as.numeric(crossprod(J, f))
  sc <- .sp_lsq_scale(J, x_scale)
  cl <- function(x, g) {
    v <- rep(1, n); dv <- numeric(n)
    mk <- g < 0 & is.finite(ub); v[mk] <- ub[mk] - x[mk]; dv[mk] <- -1
    mk <- g > 0 & is.finite(lb); v[mk] <- x[mk] - lb[mk]; dv[mk] <- 1
    list(v = v, dv = dv)
  }
  if (bounded) {
    cv <- cl(x, g)
    v <- cv$v
    v[cv$dv != 0] <- v[cv$dv != 0] * sc$inv[cv$dv != 0]
    Delta <- .sp_nrm(x0 * sc$inv / v^0.5)
  } else {
    Delta <- .sp_nrm(x0 * sc$inv)
  }
  if (Delta == 0) Delta <- 1
  lsmr <- tr_solver == "lsmr"
  if (lsmr) {
    # the options of the trf driver; the rest goes to lsmr
    reg_term <- 0
    damp <- 0
    if (!bounded) {
      damp <- tr_options$damp %||% 0
      tr_options$damp <- NULL
    }
    regularize <- tr_options$regularize %||% TRUE
    tr_options$regularize <- NULL
  }
  if (is.null(max_nfev)) max_nfev <- n * 100L
  alpha <- 0
  status <- NULL
  repeat {
    if (bounded) {
      cv <- cl(x, g)
      v <- cv$v; dv <- cv$dv
      g_norm <- max(abs(g * v))
    } else {
      g_norm <- max(abs(g))
    }
    if (g_norm < gtol) status <- 1L
    if (!is.null(status) || nfev == max_nfev) break
    if (bounded) {
      v[dv != 0] <- v[dv != 0] * sc$inv[dv != 0]
      d <- v^0.5 * sc$scale
      diag_h <- g * dv * sc$scale
    } else {
      d <- sc$scale
      diag_h <- NULL
    }
    g_h <- d * g
    if (!lsmr) {
      J_h <- J * rep(d, each = m)
      if (bounded) {
        sv <- svd(rbind(J_h, diag(diag_h^0.5, n, n)))
        uf <- as.numeric(crossprod(sv$u, c(f, numeric(n))))
      } else {
        sv <- svd(J_h)
        uf <- as.numeric(crossprod(sv$u, f))
      }
    } else {
      # the Gauss Newton step by lsmr and the two dimensional subspace
      # spanned by it and the gradient (right_multiplied_operator)
      J_h <- local({
        Jc <- J; dc <- d
        function(x) as.numeric(Jc %*% (x * dc))
      })
      J_hT <- local({
        Jc <- J; dc <- d
        function(y) dc * as.numeric(crossprod(Jc, y))
      })
      if (regularize) {
        q <- .sp_lsq_quad_1d(J_h, g_h, -g_h, diag = diag_h)
        ag_value <- .sp_lsq_min_quad_1d(q$a, q$b, 0, Delta / .sp_nrm(g_h))$y
        reg_term <- -ag_value / Delta^2
      }
      if (bounded) {
        dg <- (diag_h + reg_term)^0.5
        gn_h <- do.call(.spk_lsmr, c(list(function(x) c(J_h(x), dg * x),
                                          function(y) J_hT(y[seq_len(m)]) + dg * y[m + seq_len(n)], m + n, n,
                                          c(f, numeric(n))), tr_options))$x
      } else {
        gn_h <- do.call(.spk_lsmr, c(list(J_h, J_hT, m, n, f, damp = (damp^2 + reg_term)^0.5), tr_options))$x
      }
      S <- qr.Q(qr(cbind(g_h, gn_h), tol = 0))
      JS <- J %*% (S * d)
      B_S <- crossprod(JS)
      if (bounded) B_S <- B_S + crossprod(S * diag_h, S)
      g_S <- as.numeric(crossprod(S, g_h))
    }
    if (bounded) theta <- max(0.995, 1 - g_norm)
    actual_reduction <- -1
    while (actual_reduction <= 0 && nfev < max_nfev) {
      if (lsmr) {
        p_h <- as.numeric(S %*% .spk_tr_2d(B_S, g_S, Delta)$p)
      } else {
        tr <- .sp_lsq_trust_region(n, m, uf, sv$d, sv$v, Delta, initial_alpha = alpha)
        alpha <- tr$alpha
        p_h <- tr$p
      }
      if (bounded) {
        ss <- .sp_lsq_select_step(x, J_h, diag_h, g_h, d * p_h, p_h, d, Delta, lb, ub, theta)
        step <- ss$step; step_h <- ss$step_h; predicted_reduction <- ss$pred
        x_new <- .sp_lsq_strictly_feasible(x + step, lb, ub, rstep = 0)
      } else {
        step_h <- p_h
        predicted_reduction <- -.sp_lsq_eval_quad(J_h, g_h, step_h)
        step <- d * step_h
        x_new <- x + step
      }
      f_new <- fun(x_new)
      nfev <- nfev + 1L
      step_h_norm <- .sp_nrm(step_h)
      if (!all(is.finite(f_new))) {
        Delta <- 0.25 * step_h_norm
        next
      }
      cost_new <- if (!is.null(loss_function)) loss_function(f_new, TRUE) else 0.5 * sum(f_new * f_new)
      actual_reduction <- cost - cost_new
      ur <- .sp_lsq_update_tr(Delta, actual_reduction, predicted_reduction, step_h_norm, step_h_norm > 0.95 * Delta)
      status <- .sp_lsq_check_termination(actual_reduction, cost, .sp_nrm(step), .sp_nrm(x), ur$ratio, ftol, xtol)
      if (!is.null(status)) break
      alpha <- alpha * Delta / ur$Delta
      Delta <- ur$Delta
    }
    if (actual_reduction > 0) {
      x <- x_new; f <- f_new; f_true <- f; cost <- cost_new
      J <- jac(x)
      njev <- njev + 1L
      if (!is.null(loss_function)) {
        rb <- .sp_lsq_robust(J, f, loss_function); J <- rb$J; f <- rb$f
      }
      g <- as.numeric(crossprod(J, f))
      if (identical(x_scale, "jac")) sc <- .sp_lsq_scale(J, x_scale, sc$inv)
    }
  }
  list(x = x, cost = cost, fun = f_true, jac = J, grad = g, optimality = g_norm,
       active_mask = if (bounded) .sp_lsq_active(x, lb, ub, rtol = xtol) else integer(n), nfev = nfev, njev = njev,
       status = status %||% 0L)
}

# numpy.linalg.lstsq(A, b, rcond = -1): the minimum norm solution from the SVD
.sp_lstsq <- function(A, b) {
  sv <- svd(A)
  keep <- sv$d > .spmp_epsmch * max(sv$d)
  as.numeric(sv$v[, keep, drop = FALSE] %*% (crossprod(sv$u[, keep, drop = FALSE], b) / sv$d[keep]))
}

.sp_lsq_dogleg_step <- function(x, newton_step, g, a, b, tr_bounds, lb, ub) {
  lbc <- lb - x; ubc <- ub - x
  lbt <- pmax(lbc, -tr_bounds); ubt <- pmin(ubc, tr_bounds)
  orig_l <- lbt == lbc; orig_u <- ubt == ubc
  tr_l <- lbt == -tr_bounds; tr_u <- ubt == tr_bounds
  hits <- integer(length(x))
  if (all(newton_step >= lbt & newton_step <= ubt)) return(list(step = newton_step, hits = hits, tr_hit = FALSE))
  to_bounds <- .sp_lsq_step_to_bound(numeric(length(x)), -g, lbt, ubt)$step
  cauchy <- -.sp_lsq_min_quad_1d(a, b, 0, to_bounds)$t * g
  step_diff <- newton_step - cauchy
  sb <- .sp_lsq_step_to_bound(cauchy, step_diff, lbt, ubt)
  hits[sb$hits < 0L & orig_l] <- -1L
  hits[sb$hits > 0L & orig_u] <- 1L
  list(step = cauchy + sb$step * step_diff, hits = hits,
       tr_hit = any((sb$hits < 0L & tr_l) | (sb$hits > 0L & tr_u)))
}

# dogbox: Powell's dogleg in a rectangular trust region (Voglis and Lagaris)
.sp_lsq_dogbox <- function(fun, jac, x0, f0, J0, lb, ub, ftol, xtol, gtol, max_nfev, x_scale, loss_function,
                           tr_solver = "exact", tr_options = list()) {
  f <- f0; f_true <- f
  nfev <- 1L; njev <- 1L
  J <- J0
  n <- ncol(J)
  if (!is.null(loss_function)) {
    rb <- .sp_lsq_robust(J, f, loss_function); J <- rb$J; f <- rb$f; cost <- rb$cost
  } else {
    cost <- 0.5 * sum(f * f)
  }
  g <- as.numeric(crossprod(J, f))
  sc <- .sp_lsq_scale(J, x_scale)
  Delta <- max(abs(x0 * sc$inv))
  if (Delta == 0) Delta <- 1
  on_bound <- integer(n)
  on_bound[x0 == lb] <- -1L
  on_bound[x0 == ub] <- 1L
  x <- x0
  if (is.null(max_nfev)) max_nfev <- n * 100L
  status <- NULL
  repeat {
    active_set <- on_bound * g < 0
    free_set <- !active_set
    g_free <- g[free_set]
    g_full <- g
    g[active_set] <- 0
    g_norm <- max(abs(g))
    if (g_norm < gtol) status <- 1L
    if (!is.null(status) || nfev == max_nfev) break
    x_free <- x[free_set]
    lb_free <- lb[free_set]; ub_free <- ub[free_set]
    scale_free <- sc$scale[free_set]
    if (tr_solver == "lsmr") {
      # lsmr in scaled variables on the free columns (lsmr_operator)
      Jop <- local({
        Jc <- J
        function(x) as.numeric(Jc %*% x)
      })
      scl <- sc$scale
      act <- active_set
      lop <- function(x) as.numeric(J %*% (x * scl))
      lrop <- function(y) {
        r <- scl * as.numeric(crossprod(J, y))
        r[act] <- 0
        r
      }
      newton_step <- -do.call(.spk_lsmr, c(list(lop, lrop, nrow(J), n, f), tr_options))$x[free_set]
      newton_step <- newton_step * scale_free
      q <- .sp_lsq_quad_1d(Jop, g, -g)
    } else {
      J_free <- J[, free_set, drop = FALSE]
      newton_step <- .sp_lstsq(J_free, -f)
      q <- .sp_lsq_quad_1d(J_free, g_free, -g_free)
    }
    actual_reduction <- -1
    while (actual_reduction <= 0 && nfev < max_nfev) {
      ds <- .sp_lsq_dogleg_step(x_free, newton_step, g_free, q$a, q$b, Delta * scale_free, lb_free, ub_free)
      step <- numeric(n)
      step[free_set] <- ds$step
      predicted_reduction <- if (tr_solver == "lsmr") -.sp_lsq_eval_quad(Jop, g, step) else
        -.sp_lsq_eval_quad(J_free, g_free, ds$step)
      x_new <- pmin(pmax(x + step, lb), ub)
      f_new <- fun(x_new)
      nfev <- nfev + 1L
      step_h_norm <- max(abs(step * sc$inv))
      if (!all(is.finite(f_new))) {
        Delta <- 0.25 * step_h_norm
        next
      }
      cost_new <- if (!is.null(loss_function)) loss_function(f_new, TRUE) else 0.5 * sum(f_new * f_new)
      actual_reduction <- cost - cost_new
      ur <- .sp_lsq_update_tr(Delta, actual_reduction, predicted_reduction, step_h_norm, ds$tr_hit)
      Delta <- ur$Delta
      status <- .sp_lsq_check_termination(actual_reduction, cost, .sp_nrm(step), .sp_nrm(x), ur$ratio, ftol, xtol)
      if (!is.null(status)) break
    }
    if (actual_reduction > 0) {
      on_bound[free_set] <- ds$hits
      x <- x_new
      x[on_bound == -1L] <- lb[on_bound == -1L]
      x[on_bound == 1L] <- ub[on_bound == 1L]
      f <- f_new; f_true <- f; cost <- cost_new
      J <- jac(x)
      njev <- njev + 1L
      if (!is.null(loss_function)) {
        rb <- .sp_lsq_robust(J, f, loss_function); J <- rb$J; f <- rb$f
      }
      g <- as.numeric(crossprod(J, f))
      if (identical(x_scale, "jac")) sc <- .sp_lsq_scale(J, x_scale, sc$inv)
    }
  }
  list(x = x, cost = cost, fun = f_true, jac = J, grad = g_full, optimality = g_norm, active_mask = on_bound,
       nfev = nfev, njev = njev, status = status %||% 0L)
}

# ScipyLeastSquares: the option mapping of the Julia wrapper (bounds for
# trf and dogbox when both are given, jac "2-point", max_nfev = maxiters,
# ftol = abstol, gtol = reltol, then the solver keywords), then
# scipy.optimize.least_squares
.scipy_run_lsq <- function(alg, f, x0, lb, ub, opts, kw) {
  method <- alg$method
  loss <- alg$loss %||% "linear"
  n <- length(x0)
  ftol <- kw$ftol %||% opts$abstol %||% 1e-8
  xtol <- kw$xtol %||% 1e-8
  gtol <- kw$gtol %||% opts$reltol %||% 1e-8
  max_nfev <- kw$max_nfev %||% opts$maxiters
  if (!is.null(max_nfev)) max_nfev <- as.integer(max_nfev)
  if (!is.null(max_nfev) && max_nfev <= 0L) stop("`max_nfev` must be None or positive integer.", call. = FALSE)
  if (method == "lm" && loss != "linear") stop("method='lm' supports only 'linear' loss function.", call. = FALSE)
  use_bounds <- !is.null(lb) && !is.null(ub) && method %in% c("trf", "dogbox")
  lo <- if (use_bounds) rep_len(as.numeric(lb), n) else rep(-Inf, n)
  hi <- if (use_bounds) rep_len(as.numeric(ub), n) else rep(Inf, n)
  if (any(lo >= hi)) stop("Each lower bound must be strictly less than each upper bound.", call. = FALSE)
  if (!all(x0 >= lo & x0 <= hi)) stop("Initial guess is outside of provided bounds", call. = FALSE)
  x_scale <- kw$x_scale %||% (if (method == "lm") "jac" else 1)
  if (!identical(x_scale, "jac")) {
    if (!is.numeric(x_scale) || !all(is.finite(x_scale) & x_scale > 0))
      stop("`x_scale` must be 'jac' or array_like with positive numbers.", call. = FALSE)
    x_scale <- rep_len(as.numeric(x_scale), n)
  }
  if (method == "lm" && (ftol < .spmp_epsmch || xtol < .spmp_epsmch || gtol < .spmp_epsmch))
    stop("All tolerances must be higher than machine epsilon (2.22e-16) for method 'lm'.", call. = FALSE)
  if (ftol < .spmp_epsmch && xtol < .spmp_epsmch && gtol < .spmp_epsmch)
    stop("At least one of the tolerances must be higher than machine epsilon (2.22e-16).", call. = FALSE)
  if (method == "trf") x0 <- .sp_lsq_strictly_feasible(x0, lo, hi)
  tr_solver <- kw$tr_solver %||% "exact"
  if (!(tr_solver %in% c("exact", "lsmr"))) stop("`tr_solver` must be None, 'exact' or 'lsmr'.", call. = FALSE)
  tr_options <- kw$tr_options %||% list()
  if (method == "dogbox" && tr_solver == "lsmr" && "regularize" %in% names(tr_options)) {
    warning("The keyword 'regularize' in `tr_options` is not relevant for 'dogbox' method.", call. = FALSE)
    tr_options$regularize <- NULL
  }
  bad <- setdiff(names(tr_options), c("damp", "atol", "btol", "conlim", "maxiter", "show", "regularize"))
  if (tr_solver == "lsmr" && length(bad))
    stop("lsmr() got an unexpected keyword argument '", bad[1L], "'", call. = FALSE)
  vf <- .sp_vector_function(f, kw$diff_step, lo, hi)
  f0 <- vf$fun(x0)
  J0 <- vf$jac(x0)
  if (!all(is.finite(f0))) stop("Residuals are not finite in the initial point.", call. = FALSE)
  m <- length(f0)
  if (method == "lm") {
    if (m < n)
      stop("Method 'lm' doesn't work when the number of residuals is less than the number of variables.",
           call. = FALSE)
    r <- .spmp_lmdif(vf$fun, x0, jac = vf$jac, ftol = ftol, xtol = xtol, gtol = gtol,
                   maxfev = max_nfev %||% (100L * n), diag = if (identical(x_scale, "jac")) NULL else 1 / x_scale,
                   factor = 100)
    status <- c("0" = -1L, "1" = 2L, "2" = 3L, "3" = 4L, "4" = 1L, "5" = 0L)[[as.character(r$info)]]
    J <- vf$jac(r$x)
    g <- as.numeric(crossprod(J, r$fvec))
    res <- list(x = r$x, cost = 0.5 * .spmp_ssum(r$fvec * r$fvec), fun = r$fvec, jac = J, grad = g,
                optimality = max(abs(g)), active_mask = integer(n), nfev = r$nfev, njev = NULL, status = status)
  } else {
    lf <- .sp_lsq_loss(loss, kw$f_scale %||% 1)
    solver <- if (method == "trf") .sp_lsq_trf else .sp_lsq_dogbox
    res <- solver(vf$fun, vf$jac, x0, f0, J0, lo, hi, ftol, xtol, gtol, max_nfev, x_scale, lf, tr_solver = tr_solver,
                  tr_options = tr_options)
  }
  res$message <- .sp_lsq_messages[[as.character(res$status)]]
  res$success <- res$status > 0L
  res$objective <- res$cost
  res
}
