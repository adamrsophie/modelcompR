# minimize(method = "COBYQA") of SciPy 1.18.1: the COBYQA solver of Ragonneau
# and Zhang that SciPy vendors (scipy/_external/cobyqa: main.py,
# framework.py, models.py, problem.py, settings.py, subsolvers/geometry.py,
# subsolvers/optim.py and utils/math.py), the wrapper _minimize_cobyqa
# (scipy/optimize/_cobyqa_py.py) and the bounded variable least squares of
# lsq_linear(method = "bvls") (scipy/optimize/_lsq/bvls.py) that COBYQA uses
# for its Lagrange multipliers.
#
# Differences from SciPy:
# - the symmetric eigendecomposition of the interpolation system, the
#   pivoted QR factorizations, the SVD of the least squares solves and the
#   matrix products run in R's LAPACK and BLAS, while the arm64 SciPy wheel
#   uses Accelerate. The iterates agree with SciPy to rounding, not bit for
#   bit, so a last bit difference can change the path of a long run while
#   the solution agrees.
# - two code defects of the vendored COBYQA are fixed (see PORTING.md):
#   Problem.violation evaluates the nonlinear constraints at the reduced
#   point (fixed variables removed, scaled when scale = TRUE) instead of the
#   point of the user, and Models.update_interpolation stops updating the
#   constraint models once one model reports an ill conditioned system
#   ("or" short circuits). The internal argument .upstream_bugs = TRUE
#   restores both defects, for the comparison with SciPy.
# - with disp = TRUE, the objective is printed under the name "fun" (SciPy
#   prints the __name__ of the Python function) and the constraint values
#   are not printed (SciPy prints them only for functions with a __name__).

.spcq_eps <- .Machine$double.eps
.spcq_tiny <- .Machine$double.xmin
.spcq_barrier <- 2^100

.spcq_messages <- c("0" = "The lower bound for the trust-region radius has been reached",
                    "1" = "The target objective function value has been reached",
                    "2" = "All variables are fixed by the bound constraints",
                    "3" = "The callback requested to stop the optimization procedure",
                    "4" = "The feasibility problem received has been solved successfully",
                    "5" = "The maximum number of function evaluations has been exceeded",
                    "6" = "The maximum number of iterations has been exceeded",
                    "-1" = "The bound constraints are infeasible",
                    "-2" = "A linear algebra error occurred")

# the exceptions of COBYQA (MaxEvalError, TargetSuccess, CallbackSuccess,
# FeasibleSuccess and numpy's LinAlgError) as R conditions
.spcq_raise <- function(cls) stop(structure(class = c(cls, "error", "condition"), list(message = cls, call = NULL)))

# the Python built-ins max and min (a NaN argument after the first is skipped)
.spcq_max <- function(...) Reduce(function(a, b) if (isTRUE(b > a)) b else a, list(...))
.spcq_min <- function(...) Reduce(function(a, b) if (isTRUE(b < a)) b else a, list(...))

.spcq_dot <- function(a, b) sum(a * b)
.spcq_norm <- function(a) sqrt(sum(a * a))
# A x and A' v as sequential sums (numpy's products are these sequential
# sums up to three terms)
.spcq_mv <- function(A, x) rowSums(A * rep(x, each = nrow(A)))
.spcq_tmv <- function(A, v) colSums(A * v)
.spcq_clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

# numpy.argmax and argmin: the first NaN, else the first extremum
.spcq_argmax <- function(x) if (anyNA(x)) which(is.na(x))[1L] else which.max(x)
.spcq_argmin <- function(x) if (anyNA(x)) which(is.na(x))[1L] else which.min(x)

# the ties of COBYQA's _argmin
.spcq_argmin_all <- function(x) which(x <= min(x))

# numpy.linspace(start, stop, num)
.spcq_linspace <- function(start, stop, num) {
  div <- num - 1
  delta <- stop - start
  y <- seq.int(0, num - 1)
  step <- delta / div
  y <- if (step == 0) y / div * delta else y * step
  y <- y + start
  y[num] <- stop
  y
}

# numpy.linalg.lstsq(A, b, rcond): the minimum norm solution of the SVD
# (gelsd), singular values at or below rcond times the largest dropped;
# rcond = NULL is eps times max(m, n)
.spcq_lstsq <- function(A, b, rcond = NULL) {
  if (ncol(A) == 0L) return(numeric())
  if (is.null(rcond)) rcond <- .spcq_eps * max(dim(A))
  if (rcond < 0) rcond <- .spcq_eps
  sv <- svd(A)
  keep <- sv$d > rcond * max(sv$d)
  as.numeric(sv$v[, keep, drop = FALSE] %*% (crossprod(sv$u[, keep, drop = FALSE], b) / sv$d[keep]))
}

# utils.get_arrays_tol
.spcq_arrays_tol <- function(...) {
  arrays <- list(...)
  size <- max(vapply(arrays, length, 0L))
  weight <- max(vapply(arrays, function(a) max(c(1, abs(a[is.finite(a)]))), 0))
  10 * .spcq_eps * max(size, 1) * weight
}

# np.eye(n)[mask, ]
.spcq_eye_rows <- function(n, mask) diag(n)[mask, , drop = FALSE]

# ---------------------------------------------------------------------------
# lsq_linear(A, b, bounds = (lb, ub), method = "bvls") of SciPy: the solution
# only

.spcq_kkt <- function(g, on_bound) {
  g_kkt <- g * on_bound
  free <- on_bound == 0
  g_kkt[free] <- abs(g[free])
  max(g_kkt)
}

.spcq_lsq_linear_bvls <- function(A, b, lb, ub, tol = 1e-10) {
  n <- ncol(A)
  lb <- rep_len(lb, n)
  ub <- rep_len(ub, n)
  if (any(lb >= ub)) stop("Each lower bound must be strictly less than each upper bound.", call. = FALSE)
  x_lsq <- .spcq_lstsq(A, b, rcond = -1)
  if (all(x_lsq >= lb & x_lsq <= ub)) return(x_lsq)
  # bvls(A, b, x_lsq, lb, ub, tol, max_iter = None, verbose = 0)
  x <- x_lsq
  on_bound <- numeric(n)
  mask <- x <= lb
  x[mask] <- lb[mask]
  on_bound[mask] <- -1
  mask <- x >= ub
  x[mask] <- ub[mask]
  on_bound[mask] <- 1
  free_mask <- on_bound == 0
  active <- !free_mask
  free_set <- which(free_mask)
  r <- .spcq_mv(A, x) - b
  cost <- 0.5 * .spcq_dot(r, r)
  g <- .spcq_tmv(A, r)
  iteration <- 0L
  while (length(free_set) > 0L) {
    iteration <- iteration + 1L
    A_free <- A[, free_set, drop = FALSE]
    b_free <- b - .spcq_mv(A, x * active)
    z <- .spcq_lstsq(A_free, b_free)
    lbv <- z < lb[free_set]
    ubv <- z > ub[free_set]
    v <- lbv | ubv
    if (any(lbv)) {
      ind <- free_set[lbv]
      x[ind] <- lb[ind]
      active[ind] <- TRUE
      on_bound[ind] <- -1
    }
    if (any(ubv)) {
      ind <- free_set[ubv]
      x[ind] <- ub[ind]
      active[ind] <- TRUE
      on_bound[ind] <- 1
    }
    ind <- free_set[!v]
    x[ind] <- z[!v]
    r <- .spcq_mv(A, x) - b
    cost_new <- 0.5 * .spcq_dot(r, r)
    cost <- cost_new
    g <- .spcq_tmv(A, r)
    if (any(v)) free_set <- free_set[!v] else break
  }
  max_iter <- n + iteration
  status <- NULL
  optimality <- .spcq_kkt(g, on_bound)
  it <- iteration
  while (it < max_iter) {
    if (optimality < tol) status <- 1L
    if (!is.null(status)) break
    move_to_free <- .spcq_argmax(g * on_bound)
    on_bound[move_to_free] <- 0
    repeat {
      free_mask <- on_bound == 0
      active <- !free_mask
      free_set <- which(free_mask)
      x_free <- x[free_set]
      lb_free <- lb[free_set]
      ub_free <- ub[free_set]
      A_free <- A[, free_set, drop = FALSE]
      b_free <- b - .spcq_mv(A, x * active)
      z <- .spcq_lstsq(A_free, b_free)
      lbv <- which(z < lb_free)
      ubv <- which(z > ub_free)
      v <- c(lbv, ubv)
      if (length(v)) {
        alphas <- c(lb_free[lbv] - x_free[lbv], ub_free[ubv] - x_free[ubv]) / (z[v] - x_free[v])
        i <- .spcq_argmin(alphas)
        i_free <- v[i]
        alpha <- alphas[i]
        x_free <- x_free * (1 - alpha)
        x_free <- x_free + alpha * z
        x[free_set] <- x_free
        on_bound[free_set[i_free]] <- if (i <= length(lbv)) -1 else 1
      } else {
        x[free_set] <- z
        break
      }
    }
    r <- .spcq_mv(A, x) - b
    cost_new <- 0.5 * .spcq_dot(r, r)
    cost_change <- cost - cost_new
    if (isTRUE(cost_change < tol * cost)) status <- 2L
    cost <- cost_new
    g <- .spcq_tmv(A, r)
    optimality <- .spcq_kkt(g, on_bound)
    it <- it + 1L
  }
  x
}

# ---------------------------------------------------------------------------
# The problem (problem.py)

# LinearConstraints: the inequalities a_ub x <= b_ub and the equalities
# a_eq x = b_eq of a list of list(A, lb, ub), with the constraints kept for
# the violation
.spcq_linear <- function(cons, n) {
  a_ub <- matrix(0, 0L, n); b_ub <- numeric()
  a_eq <- matrix(0, 0L, n); b_eq <- numeric()
  for (con in cons) {
    is_eq <- abs(con$ub - con$lb) <= .spcq_arrays_tol(con$lb, con$ub)
    is_eq <- is_eq %in% TRUE
    if (any(is_eq)) {
      a_eq <- rbind(a_eq, con$A[is_eq, , drop = FALSE])
      b_eq <- c(b_eq, 0.5 * (con$lb[is_eq] + con$ub[is_eq]))
    }
    if (!all(is_eq)) {
      a_ub <- rbind(a_ub, con$A[!is_eq, , drop = FALSE], -con$A[!is_eq, , drop = FALSE])
      b_ub <- c(b_ub, con$ub[!is_eq], -con$lb[!is_eq])
    }
  }
  a_ub[is.na(a_ub)] <- 0
  a_eq[is.na(a_eq)] <- 0
  undef_ub <- is.na(b_ub) | is.infinite(b_ub)
  undef_eq <- is.na(b_eq)
  list(a_ub = a_ub[!undef_ub, , drop = FALSE], b_ub = b_ub[!undef_ub],
       a_eq = a_eq[!undef_eq, , drop = FALSE], b_eq = b_eq[!undef_eq],
       pcs = Filter(function(con) length(con$A) > 0L, cons))
}

# PreparedConstraint.violation: the excess over the bounds of each component
.spcq_excess <- function(ev, lb, ub) pmax(lb - ev, 0) + pmax(ev - ub, 0)

# BoundConstraints
.spcq_bounds <- function(xl, xu) {
  xl[is.na(xl)] <- -Inf
  xu[is.na(xu)] <- Inf
  list(xl = xl, xu = xu, is_feasible = all(xl <= xu) && all(xl < Inf) && all(xu > -Inf),
       m = sum(xl > -Inf) + sum(xu < Inf))
}

.spcq_project <- function(bd, x) if (bd$is_feasible) .spcq_clip(x, bd$xl, bd$xu) else x

# the value of a nonlinear constraint through the VectorFunction of its
# PreparedConstraint: memoized at the last point only
.spcq_nl_fun <- function(pc, x) {
  if (is.null(pc$x) || length(x) != length(pc$x) || !isTRUE(all(x == pc$x))) {
    pc$x <- x
    pc$f <- as.numeric(pc$fun(x))
    pc$nfev <- pc$nfev + 1L
  }
  pc$f
}

# NonlinearConstraints.__call__: the inequalities c_ub(x) <= 0 and the
# equalities c_eq(x) = 0
.spcq_nl_call <- function(pb, x) {
  if (length(pb$nl_cons) == 0L) {
    pb$m_nl_ub <- 0L; pb$m_nl_eq <- 0L
    return(list(cub = numeric(), ceq = numeric()))
  }
  if (length(pb$pcs) == 0L) {
    pb$pcs <- lapply(pb$nl_cons, function(con) {
      pc <- new.env(parent = emptyenv())
      pc$fun <- con$fun
      pc$nfev <- 0L
      f <- .spcq_nl_fun(pc, x)
      m <- length(f)
      pc$lb <- rep_len(as.numeric(con$lb), m)
      pc$ub <- rep_len(as.numeric(con$ub), m)
      is_eq <- (abs(pc$ub - pc$lb) <= .spcq_arrays_tol(pc$lb, pc$ub)) %in% TRUE
      pc$map_eq <- which(is_eq)
      pc$map_ub <- which(!is_eq)
      pc
    })
  }
  c_ub <- numeric(); c_eq <- numeric()
  for (pc in pb$pcs) {
    val <- .spcq_nl_fun(pc, x)
    ub_idx <- pc$map_ub
    if (length(ub_idx)) {
      ub_val <- val[ub_idx]
      xl <- pc$lb[ub_idx]
      xu <- pc$ub[ub_idx]
      fl <- xl > -Inf
      fu <- xu < Inf
      c_ub <- c(c_ub, xl[fl] - ub_val[fl], ub_val[fu] - xu[fu])
    }
    eq_idx <- pc$map_eq
    if (length(eq_idx)) c_eq <- c(c_eq, val[eq_idx] - 0.5 * (pc$ub[eq_idx] + pc$lb[eq_idx]))
  }
  pb$m_nl_ub <- length(c_ub)
  pb$m_nl_eq <- length(c_eq)
  list(cub = c_ub, ceq = c_eq)
}

.spcq_build_x <- function(pb, x) {
  x_full <- numeric(pb$n_orig)
  x_full[pb$fixed_idx] <- pb$fixed_val
  x_full[!pb$fixed_idx] <- x * pb$scaling_factor + pb$scaling_shift
  .spcq_project(pb$orig_bounds, x_full)
}

# Problem.violation (the values cub_val and ceq_val are not used: the
# constraints are evaluated again, through their memoization)
.spcq_violation <- function(pb, x) {
  v <- numeric()
  if (!pb$bounds$is_feasible) v <- c(v, .spcq_excess(x, pb$bounds$xl, pb$bounds$xu))
  for (pc in pb$linear$pcs) v <- c(v, .spcq_excess(.spcq_mv(pc$A, x), pc$lb, pc$ub))
  if (length(pb$pcs)) {
    xn <- if (pb$upstream_bugs) x else .spcq_build_x(pb, x)
    for (pc in pb$pcs) v <- c(v, .spcq_excess(.spcq_nl_fun(pc, xn), pc$lb, pc$ub))
  }
  v
}

.spcq_maxcv <- function(pb, x) {
  v <- .spcq_violation(pb, x)
  if (any(v != 0 | is.na(v))) max(v) else 0
}

.spcq_type <- function(pb) {
  if (is.null(pb$m_nl_ub)) return("nonlinearly constrained")
  if (pb$m_nl_ub > 0L || pb$m_nl_eq > 0L) "nonlinearly constrained"
  else if (nrow(pb$linear$a_ub) > 0L || nrow(pb$linear$a_eq) > 0L) "linearly constrained"
  else if (pb$bounds$m > 0L) "bound-constrained" else "unconstrained"
}

# numpy's printing of a float array under PRINT_OPTIONS of COBYQA
.spcq_format_float <- function(x) {
  if (is.na(x)) return(" nan")
  if (is.infinite(x)) return(if (x > 0) " inf" else "-inf")
  s <- sprintf("%.3e", x)
  if (x >= 0 && !grepl("^-", s)) s <- paste0(" ", s)
  s
}

.spcq_format_array <- function(x) {
  s <- vapply(x, .spcq_format_float, "")
  if (length(s) > 6L) s <- c(s[1:2], "...", s[length(s) - 1:0])
  paste0("[", paste(s, collapse = " "), "]")
}

# ObjectiveFunction
.spcq_obj <- function(pb, x) {
  if (is.null(pb$fun)) return(0)
  f <- as.numeric(pb$fun(x))
  if (length(f) != 1L) stop("The objective function must return a scalar.", call. = FALSE)
  pb$n_eval <- pb$n_eval + 1L
  if (pb$verbose) cat("fun(", .spcq_format_array(x), ") = ", .spsq_repr(f), "\n", sep = "")
  f
}

# Problem.__call__: the values at x (reduced), the filter and the callback
.spcq_eval <- function(pb, x, penalty = 0) {
  x_full <- .spcq_build_x(pb, x)
  fun_val <- .spcq_obj(pb, x_full)
  cv <- .spcq_nl_call(pb, x_full)
  cub_val <- cv$cub
  ceq_val <- cv$ceq
  maxcv_val <- .spcq_maxcv(pb, x)
  ff <- pb$fun_filter
  mf <- pb$maxcv_filter
  if (is.nan(fun_val) && is.nan(maxcv_val)) {
    include <- length(ff) == 0L
  } else if (is.nan(fun_val)) {
    include <- all((is.nan(ff) & maxcv_val < mf) | is.nan(mf))
  } else if (is.nan(maxcv_val)) {
    include <- all((is.nan(mf) & fun_val < ff) | is.nan(ff))
  } else {
    include <- all(fun_val < ff | maxcv_val < mf)
  }
  if (isTRUE(include)) {
    pb$fun_filter <- c(ff, fun_val)
    pb$maxcv_filter <- c(mf, maxcv_val)
    pb$x_filter[[length(pb$x_filter) + 1L]] <- x
    nf <- length(pb$fun_filter)
    if (nf >= 2L) {
      for (k in (nf - 1L):1L) {
        fk <- pb$fun_filter[k]; mk <- pb$maxcv_filter[k]
        remove <- if (is.nan(fun_val)) is.nan(fk) else if (is.nan(maxcv_val)) is.nan(mk) else
          is.nan(fk) || is.nan(mk) || (fun_val <= fk && maxcv_val <= mk)
        if (remove) {
          pb$fun_filter <- pb$fun_filter[-k]
          pb$maxcv_filter <- pb$maxcv_filter[-k]
          pb$x_filter <- pb$x_filter[-k]
        }
      }
    }
  }
  if (is.function(pb$callback)) {
    be <- .spcq_best_eval(pb, penalty)
    xb <- .spcq_build_x(pb, be$x)
    cb <- pb$callback
    stop_now <- tryCatch({
      if (identical(names(formals(cb)), "intermediate_result")) cb(intermediate_result = list(x = xb, fun = be$fun))
      else cb(xb)
      FALSE
    }, StopIteration = function(e) TRUE)
    if (stop_now) .spcq_raise("spcq_callback")
  }
  if (is.nan(fun_val)) fun_val <- .spcq_barrier
  cub_val[is.nan(cub_val)] <- .spcq_barrier
  ceq_val[is.nan(ceq_val)] <- .spcq_barrier
  fun_val <- max(min(fun_val, .spcq_barrier), -.spcq_barrier)
  cub_val <- pmax(pmin(cub_val, .spcq_barrier), -.spcq_barrier)
  ceq_val <- pmax(pmin(ceq_val, .spcq_barrier), -.spcq_barrier)
  list(fun = fun_val, cub = cub_val, ceq = ceq_val)
}

# Problem.best_eval
.spcq_best_eval <- function(pb, penalty) {
  if (length(pb$fun_filter) == 0L) .spcq_eval(pb, pb$x0)
  ff <- pb$fun_filter
  mf <- pb$maxcv_filter
  last <- function(m) max(which(m))
  finite <- is.finite(mf)
  if (any(finite)) {
    feas <- mf <= pb$feasibility_tol
    feas <- feas %in% TRUE
    if (any(feas) && !all(is.nan(ff[feas]))) {
      idx <- feas & (ff <= min(ff[feas], na.rm = TRUE)) %in% TRUE
      if (sum(idx) > 1L) idx <- idx & (mf <= min(mf[idx])) %in% TRUE
      i <- last(idx)
    } else if (any(feas)) {
      i <- last(feas)
    } else {
      merit <- rep(NaN, length(ff))
      merit[finite] <- ff[finite] + penalty * mf[finite]
      if (all(is.nan(merit))) {
        i <- last((mf <= min(mf, na.rm = TRUE)) %in% TRUE)
      } else {
        idx <- (merit <= min(merit, na.rm = TRUE)) %in% TRUE
        if (sum(idx) > 1L) idx <- idx & (mf <= min(mf[idx])) %in% TRUE
        if (sum(idx) > 1L) idx <- idx & (ff <= min(ff[idx])) %in% TRUE
        i <- last(idx)
      }
    }
  } else if (!all(is.nan(ff))) {
    i <- last((ff <= min(ff, na.rm = TRUE)) %in% TRUE)
  } else {
    i <- length(ff)
  }
  list(x = .spcq_project(pb$bounds, pb$x_filter[[i]]), fun = ff[i], maxcv = mf[i])
}

# Problem.__init__ with the bounds, the linear constraints and the nonlinear
# constraints of the user
.spcq_problem <- function(fun, x0, xl, xu, lin_cons, nl_cons, callback, feasibility_tol, scale, verbose,
                          upstream_bugs) {
  pb <- new.env(parent = emptyenv())
  pb$fun <- fun
  pb$n_eval <- 0L
  pb$verbose <- verbose
  pb$callback <- callback
  pb$nl_cons <- nl_cons
  pb$pcs <- list()
  pb$m_nl_ub <- NULL
  pb$m_nl_eq <- NULL
  pb$upstream_bugs <- upstream_bugs
  n0 <- length(x0)
  orig <- .spcq_bounds(xl, xu)
  linear <- .spcq_linear(lin_cons, n0)
  tol <- .spcq_arrays_tol(orig$xl, orig$xu)
  fixed <- (orig$xl <= orig$xu & abs(orig$xl - orig$xu) < tol) %in% TRUE
  fv <- 0.5 * (orig$xl[fixed] + orig$xu[fixed])
  pb$fixed_idx <- fixed
  pb$fixed_val <- .spcq_clip(fv, orig$xl[fixed], orig$xu[fixed])
  pb$n_orig <- n0
  pb$orig_bounds <- orig
  # Bounds() of SciPy rejects the empty bounds when every variable is fixed
  if (all(fixed)) stop("`lb` and `ub` must be non-empty.", call. = FALSE)
  pb$bounds <- .spcq_bounds(orig$xl[!fixed], orig$xu[!fixed])
  pb$x0 <- .spcq_project(pb$bounds, x0[!fixed])
  n <- length(pb$x0)
  b_eq <- linear$b_eq - .spcq_mv(linear$a_eq[, fixed, drop = FALSE], pb$fixed_val)
  pb$linear <- .spcq_linear(list(
    list(A = linear$a_ub[, !fixed, drop = FALSE], lb = rep(-Inf, length(linear$b_ub)),
         ub = linear$b_ub - .spcq_mv(linear$a_ub[, fixed, drop = FALSE], pb$fixed_val)),
    list(A = linear$a_eq[, !fixed, drop = FALSE], lb = b_eq, ub = b_eq)), n)
  scale <- scale && pb$bounds$is_feasible && all(is.finite(pb$bounds$xl)) && all(is.finite(pb$bounds$xu))
  if (scale) {
    pb$scaling_factor <- 0.5 * (pb$bounds$xu - pb$bounds$xl)
    pb$scaling_shift <- 0.5 * (pb$bounds$xu + pb$bounds$xl)
    pb$bounds <- .spcq_bounds(rep(-1, n), rep(1, n))
    lin <- pb$linear
    b_eq <- lin$b_eq - .spcq_mv(lin$a_eq, pb$scaling_shift)
    sf <- pb$scaling_factor
    pb$linear <- .spcq_linear(list(
      list(A = sweep(lin$a_ub, 2L, sf, "*"), lb = rep(-Inf, length(lin$b_ub)),
           ub = lin$b_ub - .spcq_mv(lin$a_ub, pb$scaling_shift)),
      list(A = sweep(lin$a_eq, 2L, sf, "*"), lb = b_eq, ub = b_eq)), n)
    pb$x0 <- (pb$x0 - pb$scaling_shift) / pb$scaling_factor
  } else {
    pb$scaling_factor <- rep(1, n)
    pb$scaling_shift <- numeric(n)
  }
  pb$feasibility_tol <- feasibility_tol
  pb$fun_filter <- numeric()
  pb$maxcv_filter <- numeric()
  pb$x_filter <- list()
  pb
}

# ---------------------------------------------------------------------------
# The interpolation and the quadratic models (models.py)

.spcq_interpolation <- function(pb, opt) {
  it <- new.env(parent = emptyenv())
  xl <- pb$bounds$xl
  xu <- pb$bounds$xu
  max_radius <- 0.5 * min(xu - xl)
  if (opt$radius_init > max_radius) {
    opt$radius_init <- max_radius
    opt$radius_final <- min(opt$radius_final, max_radius)
  }
  rb <- opt$radius_init
  n <- pb$n
  npt <- opt$npt
  xb <- pb$x0
  very_close_xl <- xb <= xl + 0.5 * rb
  xb[very_close_xl] <- xl[very_close_xl]
  close_xl <- (xl + 0.5 * rb < xb) & (xb <= xl + rb)
  xb[close_xl] <- pmin(xl[close_xl] + rb, xu[close_xl])
  very_close_xu <- xb >= xu - 0.5 * rb
  xb[very_close_xu] <- xu[very_close_xu]
  close_xu <- (xb < xu - 0.5 * rb) & (xu - rb <= xb)
  xb[close_xu] <- pmax(xu[close_xu] - rb, xl[close_xu])
  xpt <- matrix(0, n, npt)
  if (npt >= 2L) for (k in 1:(npt - 1L)) {
    if (k <= n) {
      xpt[k, k + 1L] <- if (very_close_xu[k]) -rb else rb
    } else if (k <= 2L * n) {
      i <- k - n
      xpt[i, k + 1L] <- if (very_close_xl[i]) 2 * rb else if (very_close_xu[i]) -2 * rb else -rb
    } else {
      spread <- (k - n - 1L) %/% n
      k1 <- k - (1L + spread) * n - 1L
      k2 <- (k1 + spread) %% n
      xpt[k1 + 1L, k + 1L] <- xpt[k1 + 1L, k1 + 2L]
      xpt[k2 + 1L, k + 1L] <- xpt[k2 + 1L, k2 + 2L]
    }
  }
  it$x_base <- xb
  it$xpt <- xpt
  it$cache <- NULL
  it
}

.spcq_point <- function(it, k) it$x_base + it$xpt[, k]

# build_system: the scaled interpolation matrix and its eigendecomposition
# (eigh, ascending), cached for the current points
.spcq_build_system <- function(it) {
  cc <- it$cache
  if (!is.null(cc) && identical(cc$xpt, it$xpt)) return(cc)
  xpt <- it$xpt
  n <- nrow(xpt); npt <- ncol(xpt)
  scale <- max(c(.spcq_eps, sqrt(colSums(xpt * xpt))))
  xs <- xpt / scale
  a <- matrix(0, npt + n + 1L, npt + n + 1L)
  a[1:npt, 1:npt] <- 0.5 * crossprod(xs)^2
  a[1:npt, npt + 1L] <- 1
  a[npt + 1L, 1:npt] <- 1
  if (n > 0L) {
    a[1:npt, npt + 1L + seq_len(n)] <- t(xs)
    a[npt + 1L + seq_len(n), 1:npt] <- xs
  }
  rs <- c(rep(1 / scale^2, npt), scale^2, rep(scale, n))
  e <- eigen(a, symmetric = TRUE)
  ord <- rev(seq_along(e$values))
  cc <- list(xpt = xpt, a = a, rs = rs, values = e$values[ord], vectors = e$vectors[, ord, drop = FALSE])
  it$cache <- cc
  cc
}

# Quadratic.solve_systems for the columns of rhs
.spcq_solve_systems <- function(it, rhs) {
  cc <- .spcq_build_system(it)
  rhs_scaled <- rhs * cc$rs
  if (!(all(is.finite(cc$a)) && all(is.finite(rhs_scaled)))) .spcq_raise("spcq_linalg")
  large <- abs(cc$values) > .spcq_eps
  V <- cc$vectors[, large, drop = FALSE]
  inv <- 1 / cc$values[large]
  sol <- V %*% (crossprod(V, rhs_scaled) * inv)
  list(x = sol * cc$rs, ill = !all(large))
}

.spcq_get_model <- function(it, values) {
  n <- nrow(it$xpt); npt <- ncol(it$xpt)
  s <- .spcq_solve_systems(it, matrix(c(values, numeric(n + 1L)), ncol = 1L))
  x <- s$x[, 1L]
  list(const = x[npt + 1L], grad = x[npt + 1L + seq_len(n)], i_hess = x[1:npt], ill = s$ill)
}

.spcq_quad_new <- function(it, values) {
  m <- .spcq_get_model(it, values)
  n <- nrow(it$xpt)
  list(const = m$const, grad = m$grad, i_hess = m$i_hess, e_hess = matrix(0, n, n))
}

.spcq_quad_hess_prod <- function(q, v, it) {
  .spcq_mv(q$e_hess, v) + .spcq_mv(it$xpt, q$i_hess * .spcq_tmv(it$xpt, v))
}

.spcq_quad_value <- function(q, x, it) {
  xd <- x - it$x_base
  q$const + .spcq_dot(q$grad, xd) +
    0.5 * (.spcq_dot(q$i_hess, .spcq_tmv(it$xpt, xd)^2) + .spcq_dot(.spcq_tmv(q$e_hess, xd), xd))
}

.spcq_quad_grad <- function(q, x, it) q$grad + .spcq_quad_hess_prod(q, x - it$x_base, it)

.spcq_quad_curv <- function(q, v, it) {
  .spcq_dot(.spcq_tmv(q$e_hess, v), v) + .spcq_dot(q$i_hess, .spcq_tmv(it$xpt, v)^2)
}

# Quadratic.update: the model after the point k_new moved (the new point
# already stored in it), with the change of the values
.spcq_quad_update <- function(q, it, k_new, dir_old, values_diff) {
  q$e_hess <- q$e_hess + q$i_hess[k_new] * outer(dir_old, dir_old)
  q$i_hess[k_new] <- 0
  m <- .spcq_get_model(it, values_diff)
  q$const <- q$const + m$const
  q$grad <- q$grad + m$grad
  q$i_hess <- q$i_hess + m$i_hess
  list(q = q, ill = m$ill)
}

.spcq_quad_shift <- function(q, it, new_x_base) {
  const <- .spcq_quad_value(q, new_x_base, it)
  grad <- .spcq_quad_grad(q, new_x_base, it)
  shift <- new_x_base - it$x_base
  upd <- outer(shift, .spcq_mv(it$xpt - 0.5 * shift, q$i_hess))
  q$const <- const
  q$grad <- grad
  q$e_hess <- q$e_hess + upd + t(upd)
  q
}

# Models: the interpolation points, the values and the models of the
# objective and the constraints
.spcq_models <- function(pb, opt, penalty) {
  md <- new.env(parent = emptyenv())
  md$it <- .spcq_interpolation(pb, opt)
  npt <- opt$npt
  r0 <- .spcq_eval(pb, .spcq_point(md$it, 1L), penalty)
  md$fun_val <- rep(NaN, npt)
  md$cub_val <- matrix(NaN, npt, length(r0$cub))
  md$ceq_val <- matrix(NaN, npt, length(r0$ceq))
  for (k in seq_len(npt)) {
    if (k - 1L >= opt$maxfev) .spcq_raise("spcq_maxeval")
    r <- if (k == 1L) r0 else .spcq_eval(pb, .spcq_point(md$it, k), penalty)
    md$fun_val[k] <- r$fun
    md$cub_val[k, ] <- r$cub
    md$ceq_val[k, ] <- r$ceq
    if (is.null(pb$fun) && .spcq_maxcv(pb, .spcq_point(md$it, k)) <= opt$feasibility_tol)
      .spcq_raise("spcq_feasible")
    if (md$fun_val[k] <= opt$target && .spcq_maxcv(pb, .spcq_point(md$it, k)) <= opt$feasibility_tol)
      .spcq_raise("spcq_target")
  }
  .spcq_reset_models(md)
  md
}

.spcq_reset_models <- function(md) {
  md$fun <- .spcq_quad_new(md$it, md$fun_val)
  md$cub <- lapply(seq_len(ncol(md$cub_val)), function(i) .spcq_quad_new(md$it, md$cub_val[, i]))
  md$ceq <- lapply(seq_len(ncol(md$ceq_val)), function(i) .spcq_quad_new(md$it, md$ceq_val[, i]))
}

.spcq_cub <- function(md, x, mask = NULL) {
  qs <- if (is.null(mask)) md$cub else md$cub[mask]
  vapply(qs, function(q) .spcq_quad_value(q, x, md$it), 0)
}

.spcq_ceq <- function(md, x) vapply(md$ceq, function(q) .spcq_quad_value(q, x, md$it), 0)

.spcq_grads <- function(qs, x, it, n) {
  if (length(qs) == 0L) return(matrix(0, 0L, n))
  do.call(rbind, lapply(qs, function(q) .spcq_quad_grad(q, x, it)))
}

.spcq_cub_grad <- function(md, x, mask = NULL) {
  qs <- if (is.null(mask)) md$cub else md$cub[mask]
  .spcq_grads(qs, x, md$it, nrow(md$it$xpt))
}

.spcq_ceq_grad <- function(md, x) .spcq_grads(md$ceq, x, md$it, nrow(md$it$xpt))

.spcq_update_interpolation <- function(md, k_new, x_new, fun_val, cub_val, ceq_val, upstream_bugs) {
  it <- md$it
  npt <- ncol(it$xpt)
  fun_diff <- numeric(npt)
  cub_diff <- matrix(0, npt, ncol(md$cub_val))
  ceq_diff <- matrix(0, npt, ncol(md$ceq_val))
  fun_diff[k_new] <- fun_val - .spcq_quad_value(md$fun, x_new, it)
  cub_diff[k_new, ] <- cub_val - .spcq_cub(md, x_new)
  ceq_diff[k_new, ] <- ceq_val - .spcq_ceq(md, x_new)
  md$fun_val[k_new] <- fun_val
  md$cub_val[k_new, ] <- cub_val
  md$ceq_val[k_new, ] <- ceq_val
  dir_old <- it$xpt[, k_new]
  it$xpt[, k_new] <- x_new - it$x_base
  u <- .spcq_quad_update(md$fun, it, k_new, dir_old, fun_diff)
  md$fun <- u$q
  ill <- u$ill
  # SciPy writes "ill_conditioned or model.update(...)", which leaves the
  # remaining constraint models stale once a system is ill conditioned
  for (i in seq_along(md$cub)) {
    if (upstream_bugs && ill) next
    u <- .spcq_quad_update(md$cub[[i]], it, k_new, dir_old, cub_diff[, i])
    md$cub[[i]] <- u$q
    ill <- ill || u$ill
  }
  for (i in seq_along(md$ceq)) {
    if (upstream_bugs && ill) next
    u <- .spcq_quad_update(md$ceq[[i]], it, k_new, dir_old, ceq_diff[, i])
    md$ceq[[i]] <- u$q
    ill <- ill || u$ill
  }
  ill
}

# Models.determinants: the values of sigma of the new point (for k_new, or
# for every point)
.spcq_determinants <- function(md, x_new, k_new = NULL) {
  it <- md$it
  n <- nrow(it$xpt); npt <- ncol(it$xpt)
  shift <- x_new - it$x_base
  new_col <- c(0.5 * .spcq_tmv(it$xpt, shift)^2, 1, shift)
  inv <- .spcq_solve_systems(it, matrix(new_col, ncol = 1L))$x[, 1L]
  beta <- 0.5 * .spcq_dot(shift, shift)^2 - .spcq_dot(new_col, inv)
  if (is.null(k_new)) {
    coord <- diag(1, npt + n + 1L, npt)
    alpha <- diag(.spcq_solve_systems(it, coord)$x)
    tau <- inv[1:npt]
  } else {
    coord <- matrix(0, npt + n + 1L, 1L)
    coord[k_new, 1L] <- 1
    alpha <- .spcq_solve_systems(it, coord)$x[k_new, 1L]
    tau <- inv[k_new]
  }
  alpha * beta + tau^2
}

.spcq_shift_x_base <- function(md, new_x_base) {
  it <- md$it
  md$fun <- .spcq_quad_shift(md$fun, it, new_x_base)
  md$cub <- lapply(md$cub, .spcq_quad_shift, it = it, new_x_base = new_x_base)
  md$ceq <- lapply(md$ceq, .spcq_quad_shift, it = it, new_x_base = new_x_base)
  shift <- new_x_base - it$x_base
  it$x_base <- it$x_base + shift
  it$xpt <- it$xpt - shift
}

# ---------------------------------------------------------------------------
# The subproblem solvers (subsolvers/optim.py and subsolvers/geometry.py)

# _alpha_tr: the step length to the trust region boundary, NULL for the
# ZeroDivisionError of SciPy
.spcq_alpha_tr <- function(step, sd, delta) {
  step_sd <- .spcq_dot(step, sd)
  sd_sq <- .spcq_dot(sd, sd)
  dist_tr_sq <- delta^2 - .spcq_dot(step, step)
  temp <- sqrt(.spcq_max(step_sd^2 + sd_sq * dist_tr_sq, 0))
  if (step_sd <= 0 && sd_sq > .spcq_tiny * abs(temp - step_sd)) {
    .spcq_max((temp - step_sd) / sd_sq, 0)
  } else if (abs(temp + step_sd) > .spcq_tiny * dist_tr_sq) {
    .spcq_max(dist_tr_sq / (temp + step_sd), 0)
  } else {
    NULL
  }
}

# QR factorization with column pivoting of M' (scipy.linalg.qr, mode full):
# the number of active constraints and the orthogonal factor
.spcq_qr_active <- function(M, n, factor_n) {
  if (nrow(M) == 0L) return(list(n_act = 0L, q = diag(n)))
  f <- qr(t(M), LAPACK = TRUE)
  q <- qr.Q(f, complete = TRUE)
  r <- qr.R(f)
  mn <- min(dim(t(M)))
  rr <- r[seq_len(mn), seq_len(mn), drop = FALSE]
  d <- abs(diag(rr))
  n_act <- sum(d >= 10 * .spcq_eps * factor_n * sqrt(colSums(rr * rr)))
  list(n_act = n_act, q = q)
}

.spcq_qr_tangential <- function(aub, aeq, free_xl, free_xu, free_ub) {
  n <- length(free_xl)
  M <- rbind(aeq, aub[!free_ub, , drop = FALSE], -.spcq_eye_rows(n, !free_xl), .spcq_eye_rows(n, !free_xu))
  .spcq_qr_active(M, n, n)
}

.spcq_qr_normal <- function(aub, free_xl, free_xu, free_slack, free_ub) {
  m <- nrow(aub); n <- ncol(aub)
  M <- rbind(cbind(aub[!free_ub, , drop = FALSE], -.spcq_eye_rows(m, !free_ub)),
             cbind(matrix(0, m - sum(free_slack), n), -.spcq_eye_rows(m, !free_slack)),
             cbind(-.spcq_eye_rows(n, !free_xl), matrix(0, n - sum(free_xl), m)),
             cbind(.spcq_eye_rows(n, !free_xu), matrix(0, n - sum(free_xu), m)))
  .spcq_qr_active(M, n + m, n + m)
}

# the projection q[, n_act:] q[, n_act:]' v
.spcq_proj <- function(q, n_act, v) {
  if (n_act >= ncol(q)) return(numeric(nrow(q)))
  Z <- q[, (n_act + 1L):ncol(q), drop = FALSE]
  .spcq_mv(Z, .spcq_tmv(Z, v))
}

# tangential_byrd_omojokun: the truncated conjugate gradient step of the
# bound constrained trust region subproblem
.spcq_tangential <- function(grad, hess_prod, xl, xu, delta, improve_tcg = TRUE) {
  xl <- pmin(xl, 0); xu <- pmax(xu, 0)
  n <- length(grad)
  grad_orig <- grad
  free_bd <- ((xl < 0) | (grad < 0)) & ((xu > 0) | (grad > 0))
  step <- numeric(n)
  sd <- numeric(n)
  sd[free_bd] <- -grad[free_bd]
  k <- 0L
  reduct <- 0
  boundary_reached <- FALSE
  while (k < sum(free_bd)) {
    grad_sd <- .spcq_dot(grad, sd)
    if (grad_sd >= -10 * .spcq_eps * n * .spcq_max(1, .spcq_norm(grad))) break
    alpha_tr <- .spcq_alpha_tr(step, sd, delta)
    if (is.null(alpha_tr)) break
    if (-alpha_tr * grad_sd <= 1e-08 * reduct) break
    hess_sd <- hess_prod(sd)
    curv_sd <- .spcq_dot(sd, hess_sd)
    alpha_quad <- if (curv_sd > .spcq_tiny * abs(grad_sd)) .spcq_max(-grad_sd / curv_sd, 0) else Inf
    alpha <- .spcq_min(alpha_tr, alpha_quad)
    if (-alpha * (grad_sd + 0.5 * alpha * curv_sd) <= 1e-08 * reduct) break
    i_xl <- (xl > -Inf) & (sd < -.spcq_tiny * abs(xl - step))
    i_xu <- (xu < Inf) & (sd > .spcq_tiny * abs(xu - step))
    all_xl <- rep(Inf, n); all_xu <- rep(Inf, n)
    all_xl[i_xl] <- pmax((xl[i_xl] - step[i_xl]) / sd[i_xl], 0)
    all_xu[i_xu] <- pmax((xu[i_xu] - step[i_xu]) / sd[i_xu], 0)
    alpha_xl <- min(all_xl); alpha_xu <- min(all_xu)
    alpha_bd <- .spcq_min(alpha_xl, alpha_xu)
    alpha <- .spcq_min(alpha, alpha_bd)
    if (alpha > 0) {
      step[free_bd] <- .spcq_clip(step[free_bd] + alpha * sd[free_bd], xl[free_bd], xu[free_bd])
      grad <- grad + alpha * hess_sd
      reduct <- reduct - alpha * (grad_sd + 0.5 * alpha * curv_sd)
    }
    if (alpha < .spcq_min(alpha_tr, alpha_bd)) {
      beta <- .spcq_dot(grad[free_bd], hess_sd[free_bd]) / curv_sd
      sd[free_bd] <- beta * sd[free_bd] - grad[free_bd]
      sd[!free_bd] <- 0
      k <- k + 1L
    } else if (alpha < alpha_tr) {
      if (alpha_xl <= alpha) {
        i_new <- .spcq_argmin(all_xl)
        step[i_new] <- xl[i_new]
      } else {
        i_new <- .spcq_argmin(all_xu)
        step[i_new] <- xu[i_new]
      }
      free_bd[i_new] <- FALSE
      sd[free_bd] <- -grad[free_bd]
      sd[!free_bd] <- 0
      k <- 0L
    } else {
      if (alpha_xl <= alpha) {
        i_new <- .spcq_argmin_all(all_xl)
        step[i_new] <- xl[i_new]
        free_bd[i_new] <- FALSE
      }
      if (alpha_xu <= alpha) {
        i_new <- .spcq_argmin_all(all_xu)
        step[i_new] <- xu[i_new]
        free_bd[i_new] <- FALSE
      }
      boundary_reached <- TRUE
      break
    }
  }
  if (improve_tcg && boundary_reached) {
    step_base <- step
    step_comparator <- .spcq_dot(grad_orig, step_base) + 0.5 * .spcq_dot(step_base, hess_prod(step_base))
    while (sum(free_bd) > 0L) {
      f <- free_bd
      step_sq <- .spcq_dot(step[f], step[f])
      grad_sq <- .spcq_dot(grad[f], grad[f])
      grad_step <- .spcq_dot(grad[f], step[f])
      grad_sd <- -sqrt(.spcq_max(step_sq * grad_sq - grad_step^2, 0))
      sd[f] <- grad_step * step[f] - step_sq * grad[f]
      sd[!f] <- 0
      if (grad_sd >= -1e-08 * reduct || any(grad_sd >= -.spcq_tiny * abs(sd[f]))) break
      sd[f] <- sd[f] / -grad_sd
      temp_xl <- numeric(n); temp_xu <- numeric(n)
      temp_xl[f] <- step[f]^2 + sd[f]^2 - xl[f]^2
      temp_xu[f] <- step[f]^2 + sd[f]^2 - xu[f]^2
      p <- temp_xl > 0
      temp_xl[p] <- sqrt(temp_xl[p]) - sd[p]
      p <- temp_xu > 0
      temp_xu[p] <- sqrt(temp_xu[p]) + sd[p]
      dist_xl <- pmax(step - xl, 0)
      dist_xu <- pmax(xu - step, 0)
      i_xl <- temp_xl > .spcq_tiny * dist_xl
      i_xu <- temp_xu > .spcq_tiny * dist_xu
      all_t_xl <- rep(1, n); all_t_xu <- rep(1, n)
      all_t_xl[i_xl] <- pmin(all_t_xl[i_xl], dist_xl[i_xl] / temp_xl[i_xl])
      all_t_xu[i_xu] <- pmin(all_t_xu[i_xu], dist_xu[i_xu] / temp_xu[i_xu])
      t_xl <- min(all_t_xl); t_xu <- min(all_t_xu)
      t_bd <- .spcq_min(t_xl, t_xu)
      hess_step <- hess_prod(step)
      hess_sd <- hess_prod(sd)
      curv_step <- .spcq_dot(step, hess_step)
      curv_sd <- .spcq_dot(sd, hess_sd)
      curv_step_sd <- .spcq_dot(step, hess_sd)
      n_samples <- as.integer((20 - 3) * t_bd + 3)
      t_s <- .spcq_linspace(t_bd / n_samples, t_bd, n_samples)
      sin_v <- 2 * t_s / (1 + t_s^2)
      all_reduct <- sin_v * (grad_step * t_s - grad_sd - t_s * curv_step +
                               sin_v * (t_s * curv_step_sd - 0.5 * (curv_sd - curv_step)))
      if (all(all_reduct <= 0, na.rm = FALSE) %in% TRUE) break
      i_max <- .spcq_argmax(all_reduct)
      cos_v <- (1 - t_s[i_max]^2) / (1 + t_s[i_max]^2)
      step[f] <- cos_v * step[f] + sin_v[i_max] * sd[f]
      grad <- grad + (cos_v - 1) * hess_step + sin_v[i_max] * hess_sd
      reduct <- reduct + all_reduct[i_max]
      if (t_bd < 1 && i_max == n_samples) {
        if (t_xl <= t_bd) {
          i_new <- .spcq_argmin_all(all_t_xl)
          step[i_new] <- xl[i_new]
          free_bd[i_new] <- FALSE
        }
        if (t_xu <= t_bd) {
          i_new <- .spcq_argmin_all(all_t_xu)
          step[i_new] <- xu[i_new]
          free_bd[i_new] <- FALSE
        }
      } else {
        break
      }
    }
    if (.spcq_dot(grad_orig, step) + 0.5 * .spcq_dot(step, hess_prod(step)) > step_comparator) step <- step_base
  }
  step
}

# constrained_tangential_byrd_omojokun: the same with the linearized
# inequality constraints aub s <= bub and the equalities aeq s = 0
.spcq_constrained_tangential <- function(grad, hess_prod, xl, xu, aub, bub, aeq, delta, improve_tcg = TRUE) {
  xl <- pmin(xl, 0); xu <- pmax(xu, 0)
  bub <- pmax(bub, 0)
  n <- length(grad)
  grad_orig <- grad
  free_xl <- (xl < 0) | (grad < 0)
  free_xu <- (xu > 0) | (grad > 0)
  free_ub <- (bub > 0) | (.spcq_mv(aub, grad) > 0)
  qa <- .spcq_qr_tangential(aub, aeq, free_xl, free_xu, free_ub)
  n_act <- qa$n_act; q <- qa$q
  step <- numeric(n)
  sd <- -.spcq_proj(q, n_act, grad)
  resid <- bub
  k <- 0L
  reduct <- 0
  boundary_reached <- FALSE
  while (k < n - n_act) {
    grad_sd <- .spcq_dot(grad, sd)
    if (grad_sd >= -10 * .spcq_eps * n * .spcq_max(1, .spcq_norm(grad))) break
    alpha_tr <- .spcq_alpha_tr(step, sd, delta)
    if (is.null(alpha_tr)) break
    if (-alpha_tr * grad_sd <= 1e-08 * reduct) break
    hess_sd <- hess_prod(sd)
    curv_sd <- .spcq_dot(sd, hess_sd)
    alpha_quad <- if (curv_sd > .spcq_tiny * abs(grad_sd)) .spcq_max(-grad_sd / curv_sd, 0) else Inf
    alpha <- .spcq_min(alpha_tr, alpha_quad)
    if (-alpha * (grad_sd + 0.5 * alpha * curv_sd) <= 1e-08 * reduct) break
    i_xl <- free_xl & (xl > -Inf) & (sd < -.spcq_tiny * abs(xl - step))
    i_xu <- free_xu & (xu < Inf) & (sd > .spcq_tiny * abs(xu - step))
    all_xl <- rep(Inf, n); all_xu <- rep(Inf, n)
    all_xl[i_xl] <- pmax((xl[i_xl] - step[i_xl]) / sd[i_xl], 0)
    all_xu[i_xu] <- pmax((xu[i_xu] - step[i_xu]) / sd[i_xu], 0)
    alpha_xl <- min(all_xl); alpha_xu <- min(all_xu)
    alpha_bd <- .spcq_min(alpha_xl, alpha_xu)
    aub_sd <- .spcq_mv(aub, sd)
    i_ub <- free_ub & (aub_sd > .spcq_tiny * abs(resid))
    all_ub <- rep(Inf, length(bub))
    all_ub[i_ub] <- resid[i_ub] / aub_sd[i_ub]
    alpha_ub <- min(c(Inf, all_ub))
    alpha <- .spcq_min(alpha, alpha_bd, alpha_ub)
    if (alpha > 0) {
      step <- .spcq_clip(step + alpha * sd, xl, xu)
      grad <- grad + alpha * hess_sd
      resid <- pmax(0, resid - alpha * aub_sd)
      reduct <- reduct - alpha * (grad_sd + 0.5 * alpha * curv_sd)
    }
    if (alpha < .spcq_min(alpha_tr, alpha_bd, alpha_ub)) {
      grad_proj <- .spcq_proj(q, n_act, grad)
      beta <- .spcq_dot(grad_proj, hess_sd) / curv_sd
      sd <- beta * sd - grad_proj
      k <- k + 1L
    } else if (alpha < alpha_tr) {
      if (alpha_xl <= alpha) {
        i_new <- .spcq_argmin(all_xl)
        step[i_new] <- xl[i_new]
        free_xl[i_new] <- FALSE
      } else if (alpha_xu <= alpha) {
        i_new <- .spcq_argmin(all_xu)
        step[i_new] <- xu[i_new]
        free_xu[i_new] <- FALSE
      } else {
        i_new <- .spcq_argmin(all_ub)
        free_ub[i_new] <- FALSE
      }
      qa <- .spcq_qr_tangential(aub, aeq, free_xl, free_xu, free_ub)
      n_act <- qa$n_act; q <- qa$q
      sd <- -.spcq_proj(q, n_act, grad)
      k <- 0L
    } else {
      if (alpha_xl <= alpha) {
        i_new <- .spcq_argmin_all(all_xl)
        step[i_new] <- xl[i_new]
        free_xl[i_new] <- FALSE
      }
      if (alpha_xu <= alpha) {
        i_new <- .spcq_argmin_all(all_xu)
        step[i_new] <- xu[i_new]
        free_xu[i_new] <- FALSE
      }
      if (alpha_ub <= alpha) {
        i_new <- .spcq_argmin_all(all_ub)
        free_ub[i_new] <- FALSE
      }
      qa <- .spcq_qr_tangential(aub, aeq, free_xl, free_xu, free_ub)
      n_act <- qa$n_act; q <- qa$q
      boundary_reached <- TRUE
      break
    }
  }
  if (improve_tcg && boundary_reached && n_act < n) {
    step_base <- step
    while (n_act < n) {
      step_proj <- .spcq_proj(q, n_act, step)
      grad_proj <- .spcq_proj(q, n_act, grad)
      step_sq <- .spcq_dot(step_proj, step_proj)
      grad_sq <- .spcq_dot(grad_proj, grad_proj)
      grad_step <- .spcq_dot(grad_proj, step_proj)
      grad_sd <- -sqrt(.spcq_max(step_sq * grad_sq - grad_step^2, 0))
      sd <- .spcq_proj(q, n_act, grad_step * step - step_sq * grad)
      if (grad_sd >= -1e-08 * reduct || any(grad_sd >= -.spcq_tiny * abs(sd))) break
      sd <- sd / -grad_sd
      temp_xl <- numeric(n); temp_xu <- numeric(n)
      dist_xl <- pmax(step - xl, 0)
      dist_xu <- pmax(xu - step, 0)
      temp_xl[free_xl] <- sd[free_xl]^2 - dist_xl[free_xl] * (dist_xl[free_xl] - 2 * step_proj[free_xl])
      temp_xu[free_xu] <- sd[free_xu]^2 - dist_xu[free_xu] * (dist_xu[free_xu] + 2 * step_proj[free_xu])
      p <- temp_xl > 0
      temp_xl[p] <- sqrt(temp_xl[p]) - sd[p]
      p <- temp_xu > 0
      temp_xu[p] <- sqrt(temp_xu[p]) + sd[p]
      i_xl <- temp_xl > .spcq_tiny * dist_xl
      i_xu <- temp_xu > .spcq_tiny * dist_xu
      all_t_xl <- rep(1, n); all_t_xu <- rep(1, n)
      all_t_xl[i_xl] <- pmin(all_t_xl[i_xl], dist_xl[i_xl] / temp_xl[i_xl])
      all_t_xu[i_xu] <- pmin(all_t_xu[i_xu], dist_xu[i_xu] / temp_xu[i_xu])
      t_xl <- min(all_t_xl); t_xu <- min(all_t_xu)
      t_bd <- .spcq_min(t_xl, t_xu)
      temp_ub <- numeric(length(resid))
      aub_step <- .spcq_mv(aub, step_proj)
      aub_sd <- .spcq_mv(aub, sd)
      temp_ub[free_ub] <- aub_sd[free_ub]^2 - resid[free_ub] * (resid[free_ub] + 2 * aub_step[free_ub])
      p <- temp_ub > 0
      temp_ub[p] <- sqrt(temp_ub[p]) + aub_sd[p]
      i_ub <- temp_ub > .spcq_tiny * resid
      all_t_ub <- rep(1, length(resid))
      all_t_ub[i_ub] <- pmin(all_t_ub[i_ub], resid[i_ub] / temp_ub[i_ub])
      t_ub <- min(c(1, all_t_ub))
      t_min <- .spcq_min(t_bd, t_ub)
      hess_step <- hess_prod(step_proj)
      hess_sd <- hess_prod(sd)
      curv_step <- .spcq_dot(step_proj, hess_step)
      curv_sd <- .spcq_dot(sd, hess_sd)
      curv_step_sd <- .spcq_dot(step_proj, hess_sd)
      n_samples <- as.integer((20 - 3) * t_min + 3)
      t_s <- .spcq_linspace(t_min / n_samples, t_min, n_samples)
      sin_v <- 2 * t_s / (1 + t_s^2)
      all_reduct <- sin_v * (grad_step * t_s - grad_sd -
                               sin_v * (0.5 * t_s^2 * curv_step - 2 * t_s * curv_step_sd + 0.5 * curv_sd))
      if (all(all_reduct <= 0) %in% TRUE) break
      i_max <- .spcq_argmax(all_reduct)
      cos_v <- (1 - t_s[i_max]^2) / (1 + t_s[i_max]^2)
      step <- .spcq_clip(step + (cos_v - 1) * step_proj + sin_v[i_max] * sd, xl, xu)
      grad <- grad + (cos_v - 1) * hess_step + sin_v[i_max] * hess_sd
      resid <- pmax(0, resid - (cos_v - 1) * aub_step - sin_v[i_max] * aub_sd)
      reduct <- reduct + all_reduct[i_max]
      if (t_min < 1 && i_max == n_samples) {
        if (t_xl <= t_min) {
          i_new <- .spcq_argmin_all(all_t_xl)
          step[i_new] <- xl[i_new]
          free_xl[i_new] <- FALSE
        }
        if (t_xu <= t_min) {
          i_new <- .spcq_argmin_all(all_t_xu)
          step[i_new] <- xu[i_new]
          # SciPy frees free_xl here (sic)
          free_xl[i_new] <- FALSE
        }
        if (t_ub <= t_min) {
          i_new <- .spcq_argmin_all(all_t_ub)
          free_ub[i_new] <- FALSE
        }
        qa <- .spcq_qr_tangential(aub, aeq, free_xl, free_xu, free_ub)
        n_act <- qa$n_act; q <- qa$q
      } else {
        break
      }
    }
    if (.spcq_dot(grad_orig, step) + 0.5 * .spcq_dot(step, hess_prod(step)) >
        .spcq_dot(grad_orig, step_base) + 0.5 * .spcq_dot(step_base, hess_prod(step_base))) step <- step_base
  }
  step
}

# normal_byrd_omojokun: the step that reduces the violation of the
# linearized constraints within the trust region
.spcq_normal <- function(aub, bub, aeq, beq, xl, xu, delta, improve_tcg = TRUE) {
  xl <- pmin(xl, 0); xu <- pmax(xu, 0)
  m_ub <- nrow(aub); n <- ncol(aub)
  ix <- seq_len(n); is <- n + seq_len(m_ub)
  grad <- c(.spcq_tmv(aeq, -beq), pmax(0, -bub))
  free_xl <- (xl < 0) | (grad[ix] < 0)
  free_xu <- (xu > 0) | (grad[ix] > 0)
  free_slack <- bub < 0
  free_ub <- (bub > 0) | (.spcq_mv(aub, grad[ix]) - grad[is] > 0)
  qa <- .spcq_qr_normal(aub, free_xl, free_xu, free_slack, free_ub)
  n_act <- qa$n_act; q <- qa$q
  delta_slack <- sqrt(.spcq_dot(beq, beq) + .spcq_dot(grad[is], grad[is]))
  step <- numeric(n)
  sd <- -.spcq_proj(q, n_act, grad)
  resid <- bub + grad[is]
  k <- 0L
  reduct <- 0
  boundary_reached <- FALSE
  while (k < n + m_ub - n_act) {
    grad_sd <- .spcq_dot(grad, sd)
    if (grad_sd >= -10 * .spcq_eps * n * .spcq_max(1, .spcq_norm(grad))) break
    alpha_tr <- .spcq_alpha_tr(step, sd[ix], delta)
    if (is.null(alpha_tr)) alpha_tr <- Inf
    a2 <- .spcq_alpha_tr(grad[is], sd[is], delta_slack)
    if (!is.null(a2)) alpha_tr <- .spcq_min(alpha_tr, a2)
    if (-alpha_tr * grad_sd <= 1e-08 * reduct) break
    hess_sd <- c(.spcq_tmv(aeq, .spcq_mv(aeq, sd[ix])), sd[is])
    curv_sd <- .spcq_dot(sd, hess_sd)
    alpha_quad <- if (curv_sd > .spcq_tiny * abs(grad_sd)) .spcq_max(-grad_sd / curv_sd, 0) else Inf
    alpha <- .spcq_min(alpha_tr, alpha_quad)
    if (-alpha * (grad_sd + 0.5 * alpha * curv_sd) <= 1e-08 * reduct) break
    sdx <- sd[ix]; sds <- sd[is]; gs <- grad[is]
    i_xl <- free_xl & (xl > -Inf) & (sdx < -.spcq_tiny * abs(xl - step))
    i_xu <- free_xu & (xu < Inf) & (sdx > .spcq_tiny * abs(xu - step))
    i_slack <- free_slack & (sds < -.spcq_tiny * abs(gs))
    all_xl <- rep(Inf, n); all_xu <- rep(Inf, n); all_slack <- rep(Inf, m_ub)
    all_xl[i_xl] <- pmax((xl[i_xl] - step[i_xl]) / sdx[i_xl], 0)
    all_xu[i_xu] <- pmax((xu[i_xu] - step[i_xu]) / sdx[i_xu], 0)
    all_slack[i_slack] <- pmax(-gs[i_slack] / sds[i_slack], 0)
    alpha_xl <- min(all_xl); alpha_xu <- min(all_xu)
    alpha_slack <- min(c(Inf, all_slack))
    alpha_bd <- .spcq_min(alpha_xl, alpha_xu, alpha_slack)
    aub_sd <- .spcq_mv(aub, sdx) - sds
    i_ub <- free_ub & (aub_sd > .spcq_tiny * abs(resid))
    all_ub <- rep(Inf, m_ub)
    all_ub[i_ub] <- resid[i_ub] / aub_sd[i_ub]
    alpha_ub <- min(c(Inf, all_ub))
    alpha <- .spcq_min(alpha, alpha_bd, alpha_ub)
    if (alpha > 0) {
      step <- .spcq_clip(step + alpha * sdx, xl, xu)
      grad <- grad + alpha * hess_sd
      resid <- pmax(0, resid - alpha * aub_sd)
      reduct <- reduct - alpha * (grad_sd + 0.5 * alpha * curv_sd)
    }
    if (alpha < .spcq_min(alpha_tr, alpha_bd, alpha_ub)) {
      grad_proj <- .spcq_proj(q, n_act, grad)
      beta <- .spcq_dot(grad_proj, hess_sd) / curv_sd
      sd <- beta * sd - grad_proj
      k <- k + 1L
    } else if (alpha < alpha_tr) {
      if (alpha_xl <= alpha) {
        i_new <- .spcq_argmin(all_xl)
        step[i_new] <- xl[i_new]
        free_xl[i_new] <- FALSE
      } else if (alpha_xu <= alpha) {
        i_new <- .spcq_argmin(all_xu)
        step[i_new] <- xu[i_new]
        free_xu[i_new] <- FALSE
      } else if (alpha_slack <= alpha) {
        i_new <- .spcq_argmin(all_slack)
        free_slack[i_new] <- FALSE
      } else {
        i_new <- .spcq_argmin(all_ub)
        free_ub[i_new] <- FALSE
      }
      qa <- .spcq_qr_normal(aub, free_xl, free_xu, free_slack, free_ub)
      n_act <- qa$n_act; q <- qa$q
      sd <- -.spcq_proj(q, n_act, grad)
      k <- 0L
    } else {
      if (alpha_xl <= alpha) {
        i_new <- .spcq_argmin_all(all_xl)
        step[i_new] <- xl[i_new]
        free_xl[i_new] <- FALSE
      }
      if (alpha_xu <= alpha) {
        i_new <- .spcq_argmin_all(all_xu)
        step[i_new] <- xu[i_new]
        free_xu[i_new] <- FALSE
      }
      boundary_reached <- TRUE
      break
    }
  }
  resid_fun <- function(s) list(ub = pmax(.spcq_mv(aub, s) - bub, 0), eq = .spcq_mv(aeq, s) - beq)
  if (improve_tcg && boundary_reached) {
    step_base <- step
    free_bd <- free_xl & free_xu
    rr <- resid_fun(step)
    grad <- .spcq_tmv(aub, rr$ub) + .spcq_tmv(aeq, rr$eq)
    sd <- numeric(n)
    while (sum(free_bd) > 0L) {
      f <- free_bd
      step_sq <- .spcq_dot(step[f], step[f])
      grad_sq <- .spcq_dot(grad[f], grad[f])
      grad_step <- .spcq_dot(grad[f], step[f])
      grad_sd <- -sqrt(.spcq_max(step_sq * grad_sq - grad_step^2, 0))
      sd[f] <- grad_step * step[f] - step_sq * grad[f]
      sd[!f] <- 0
      if (grad_sd >= -1e-08 * reduct || any(grad_sd >= -.spcq_tiny * abs(sd[f]))) break
      sd[f] <- sd[f] / -grad_sd
      temp_xl <- numeric(n); temp_xu <- numeric(n)
      temp_xl[f] <- step[f]^2 + sd[f]^2 - xl[f]^2
      temp_xu[f] <- step[f]^2 + sd[f]^2 - xu[f]^2
      p <- temp_xl > 0
      temp_xl[p] <- sqrt(temp_xl[p]) - sd[p]
      p <- temp_xu > 0
      temp_xu[p] <- sqrt(temp_xu[p]) + sd[p]
      dist_xl <- pmax(step - xl, 0)
      dist_xu <- pmax(xu - step, 0)
      i_xl <- temp_xl > .spcq_tiny * dist_xl
      i_xu <- temp_xu > .spcq_tiny * dist_xu
      all_t_xl <- rep(1, n); all_t_xu <- rep(1, n)
      all_t_xl[i_xl] <- pmin(all_t_xl[i_xl], dist_xl[i_xl] / temp_xl[i_xl])
      all_t_xu[i_xu] <- pmin(all_t_xu[i_xu], dist_xu[i_xu] / temp_xu[i_xu])
      t_xl <- min(all_t_xl); t_xu <- min(all_t_xu)
      t_bd <- .spcq_min(t_xl, t_xu)
      n_samples <- as.integer((20 - 3) * t_bd + 3)
      t_s <- .spcq_linspace(t_bd / n_samples, t_bd, n_samples)
      rr <- resid_fun(step)
      step_proj <- step
      step_proj[!f] <- 0
      base <- .spcq_dot(rr$ub, rr$ub) + .spcq_dot(rr$eq, rr$eq)
      all_reduct <- vapply(seq_len(n_samples), function(i) {
        sin_value <- 2 * t_s[i] / (1 + t_s[i]^2)
        step_alt <- .spcq_clip(step + sin_value * (sd - t_s[i] * step_proj), xl, xu)
        ra <- resid_fun(step_alt)
        0.5 * (.spcq_dot(rr$ub, rr$ub) + .spcq_dot(rr$eq, rr$eq) - .spcq_dot(ra$ub, ra$ub) - .spcq_dot(ra$eq, ra$eq))
      }, 0)
      if (all(all_reduct <= 0) %in% TRUE) break
      i_max <- .spcq_argmax(all_reduct)
      cos_v <- (1 - t_s[i_max]^2) / (1 + t_s[i_max]^2)
      sin_value <- 2 * t_s[i_max] / (1 + t_s[i_max]^2)
      step[f] <- cos_v * step[f] + sin_value * sd[f]
      rr <- resid_fun(step)
      grad <- .spcq_tmv(aub, rr$ub) + .spcq_tmv(aeq, rr$eq)
      reduct <- reduct + all_reduct[i_max]
      if (t_bd < 1 && i_max == n_samples) {
        if (t_xl <= t_bd) {
          i_new <- .spcq_argmin_all(all_t_xl)
          step[i_new] <- xl[i_new]
          free_bd[i_new] <- FALSE
        }
        if (t_xu <= t_bd) {
          i_new <- .spcq_argmin_all(all_t_xu)
          step[i_new] <- xu[i_new]
          free_bd[i_new] <- FALSE
        }
      } else {
        break
      }
    }
    rs <- resid_fun(step)
    rb <- resid_fun(step_base)
    if (.spcq_dot(rs$ub, rs$ub) + .spcq_dot(rs$eq, rs$eq) > .spcq_dot(rb$ub, rb$ub) + .spcq_dot(rb$eq, rb$eq))
      step <- step_base
  }
  step
}

# cauchy_geometry and spider_geometry: steps that make |q| large, for the
# Lagrange polynomial q of the point to replace
.spcq_cauchy_geom1 <- function(const, grad, curv, xl, xu, delta) {
  fixed_xl <- (xl < 0) & (grad > 0)
  fixed_xu <- (xu > 0) & (grad < 0)
  cs <- numeric(length(grad))
  cs[fixed_xl] <- xl[fixed_xl]
  cs[fixed_xu] <- xu[fixed_xu]
  if (.spcq_norm(cs) > delta) {
    working <- fixed_xl | fixed_xu
    repeat {
      g_norm <- .spcq_norm(grad[working])
      delta_reduced <- sqrt(delta^2 - .spcq_dot(cs[!working], cs[!working]))
      if (g_norm > .spcq_tiny * abs(delta_reduced)) mu <- .spcq_max(delta_reduced / g_norm, 0) else break
      cs[working] <- mu * grad[working]
      fixed_xl <- working & (cs < xl)
      fixed_xu <- working & (cs > xu)
      if (!any(fixed_xl) && !any(fixed_xu)) break
      cs[fixed_xl] <- xl[fixed_xl]
      cs[fixed_xu] <- xu[fixed_xu]
      working <- working & !(fixed_xl | fixed_xu)
    }
  }
  grad_step <- .spcq_dot(grad, cs)
  if (grad_step >= 0) {
    s_norm <- .spcq_norm(cs)
    alpha_tr <- if (s_norm > .spcq_tiny * delta) .spcq_max(delta / s_norm, 0) else 0
    curv_step <- curv(cs)
    alpha_quad <- if (curv_step < -.spcq_tiny * grad_step) .spcq_max(-grad_step / curv_step, 0) else Inf
    i_xl <- (xl > -Inf) & (cs < .spcq_tiny * xl)
    i_xu <- (xu < Inf) & (cs > .spcq_tiny * xu)
    alpha_xl <- min(c(Inf, xl[i_xl] / cs[i_xl]))
    alpha_xu <- min(c(Inf, xu[i_xu] / cs[i_xu]))
    alpha_bd <- .spcq_min(alpha_xl, alpha_xu)
    alpha <- .spcq_min(alpha_tr, alpha_quad, alpha_bd)
    list(step = .spcq_clip(alpha * cs, xl, xu), q = const + alpha * grad_step + 0.5 * alpha^2 * curv_step)
  } else {
    list(step = numeric(length(grad)), q = const)
  }
}

.spcq_cauchy_geometry <- function(const, grad, curv, xl, xu, delta) {
  xl <- pmin(xl, 0); xu <- pmax(xu, 0)
  r1 <- .spcq_cauchy_geom1(const, grad, curv, xl, xu, delta)
  r2 <- .spcq_cauchy_geom1(-const, -grad, function(x) -curv(x), xl, xu, delta)
  if (abs(r1$q) >= abs(r2$q)) r1$step else r2$step
}

.spcq_spider_geometry <- function(const, grad, curv, xpt, xl, xu, delta) {
  xl <- pmin(xl, 0); xu <- pmax(xu, 0)
  step <- numeric(length(grad))
  q_val <- const
  s_norm <- sqrt(colSums(xpt * xpt))
  XT <- t(xpt)
  XL <- matrix(xl, nrow(XT), ncol(XT), byrow = TRUE)
  XU <- matrix(xu, nrow(XT), ncol(XT), byrow = TRUE)
  # as in SciPy, each of these is one number for all the points (the
  # maximum over every point and coordinate)
  glob <- function(mask, B, init) max(c(init, B[mask] / XT[mask]))
  a_xl_pos <- glob((XL > -Inf) & (XT > -.spcq_tiny * XL), XL, -Inf)
  a_xl_neg <- glob((XL > -Inf) & (XT < .spcq_tiny * XL), XL, Inf)
  a_xu_neg <- glob((XU < Inf) & (XT < -.spcq_tiny * XU), XU, -Inf)
  a_xu_pos <- glob((XU < Inf) & (XT > .spcq_tiny * XU), XU, Inf)
  for (k in seq_len(ncol(xpt))) {
    if (s_norm[k] > .spcq_tiny * delta) alpha_tr <- .spcq_max(delta / s_norm[k], 0) else next
    alpha_bd_pos <- .spcq_max(.spcq_min(a_xu_pos, a_xl_neg), 0)
    alpha_bd_neg <- .spcq_min(.spcq_max(a_xl_pos, a_xu_neg), 0)
    d <- xpt[, k]
    gs <- .spcq_dot(grad, d)
    cs <- curv(d)
    aqp <- if ((gs >= 0 && cs < -.spcq_tiny * gs) || (gs <= 0 && cs > -.spcq_tiny * gs)) .spcq_max(-gs / cs, 0) else Inf
    aqn <- if ((gs >= 0 && cs > .spcq_tiny * gs) || (gs <= 0 && cs < .spcq_tiny * gs)) .spcq_min(-gs / cs, 0) else -Inf
    alpha_pos <- .spcq_min(alpha_tr, alpha_bd_pos)
    alpha_neg <- .spcq_max(-alpha_tr, alpha_bd_neg)
    q_pos <- const + alpha_pos * gs + 0.5 * alpha_pos^2 * cs
    q_neg <- const + alpha_neg * gs + 0.5 * alpha_neg^2 * cs
    if (aqp < alpha_pos) {
      qq <- const + aqp * gs + 0.5 * aqp^2 * cs
      if (abs(qq) > abs(q_pos)) {
        alpha_pos <- aqp
        q_pos <- qq
      }
    }
    if (aqn > alpha_neg) {
      qq <- const + aqn * gs + 0.5 * aqn^2 * cs
      if (abs(qq) > abs(q_neg)) {
        alpha_neg <- aqn
        q_neg <- qq
      }
    }
    if (abs(q_pos) >= abs(q_neg) && abs(q_pos) > abs(q_val)) {
      step <- .spcq_clip(alpha_pos * d, xl, xu)
      q_val <- q_pos
    } else if (abs(q_neg) > abs(q_pos) && abs(q_neg) > abs(q_val)) {
      step <- .spcq_clip(alpha_neg * d, xl, xu)
      q_val <- q_neg
    }
  }
  step
}

# ---------------------------------------------------------------------------
# The trust region framework (framework.py)

.spcq_x_best <- function(fw) .spcq_point(fw$md$it, fw$best)

.spcq_set_radius <- function(fw, radius) {
  fw$radius <- radius
  if (radius <= fw$cst$decrease_radius_threshold * fw$resolution) fw$radius <- fw$resolution
}

.spcq_lag_hess_prod <- function(fw, v) {
  md <- fw$md
  out <- .spcq_quad_hess_prod(md$fun, v, md$it)
  for (i in seq_along(md$cub)) out <- out + fw$lm_nl_ub[i] * .spcq_quad_hess_prod(md$cub[[i]], v, md$it)
  for (i in seq_along(md$ceq)) out <- out + fw$lm_nl_eq[i] * .spcq_quad_hess_prod(md$ceq[[i]], v, md$it)
  out
}

.spcq_fun_grad <- function(fw, x) .spcq_quad_grad(fw$md$fun, x, fw$md$it)

.spcq_sqp_fun <- function(fw, step) {
  .spcq_dot(step, .spcq_fun_grad(fw, .spcq_x_best(fw)) + 0.5 * .spcq_lag_hess_prod(fw, step))
}

.spcq_merit <- function(fw, x, fun_val = NULL, cub_val = NULL, ceq_val = NULL) {
  if (is.null(fun_val) || is.null(cub_val) || is.null(ceq_val)) fun_val <- .spcq_eval(fw$pb, x, fw$penalty)$fun
  m_val <- fun_val
  if (fw$penalty > 0) {
    c_val <- .spcq_violation(fw$pb, x)
    if (any(c_val != 0 | is.na(c_val))) m_val <- m_val + fw$penalty * .spcq_norm(c_val)
  }
  m_val
}

.spcq_linearizations <- function(fw, x) {
  lin <- fw$pb$linear
  md <- fw$md
  list(aub = rbind(lin$a_ub, .spcq_cub_grad(md, x)), bub = c(lin$b_ub - .spcq_mv(lin$a_ub, x), -.spcq_cub(md, x)),
       aeq = rbind(lin$a_eq, .spcq_ceq_grad(md, x)), beq = c(lin$b_eq - .spcq_mv(lin$a_eq, x), -.spcq_ceq(md, x)))
}

.spcq_trust_region_step <- function(fw) {
  xb <- .spcq_x_best(fw)
  L <- .spcq_linearizations(fw, xb)
  xl <- fw$pb$bounds$xl - xb
  xu <- fw$pb$bounds$xu - xb
  improve <- fw$cst$improve_tcg
  normal <- .spcq_normal(L$aub, L$bub, L$aeq, L$beq, xl, xu, fw$cst$byrd_omojokun_factor * fw$radius, improve)
  radius <- sqrt(fw$radius^2 - .spcq_dot(normal, normal))
  xl <- xl - normal
  xu <- xu - normal
  bub <- pmax(L$bub - .spcq_mv(L$aub, normal), 0)
  g_best <- .spcq_fun_grad(fw, xb) + .spcq_lag_hess_prod(fw, normal)
  hp <- function(v) .spcq_lag_hess_prod(fw, v)
  tangential <- if (.spcq_type(fw$pb) %in% c("unconstrained", "bound-constrained")) {
    .spcq_tangential(g_best, hp, xl, xu, radius, improve)
  } else {
    .spcq_constrained_tangential(g_best, hp, xl, xu, L$aub, bub, L$aeq, radius, improve)
  }
  list(normal = normal, tangential = tangential)
}

.spcq_geometry_step <- function(fw, k_new) {
  md <- fw$md
  it <- md$it
  npt <- ncol(it$xpt)
  xb <- .spcq_x_best(fw)
  coord <- numeric(npt)
  coord[k_new] <- 1
  lag <- .spcq_quad_new(it, coord)
  g_lag <- .spcq_quad_grad(lag, xb, it)
  curv <- function(v) .spcq_quad_curv(lag, v, it)
  xl <- fw$pb$bounds$xl - xb
  xu <- fw$pb$bounds$xu - xb
  step <- .spcq_cauchy_geometry(0, g_lag, curv, xl, xu, fw$radius)
  sigma <- .spcq_determinants(md, xb + step, k_new)
  xpt <- it$xpt - it$xpt[, fw$best]
  xpt[, c(1L, fw$best)] <- xpt[, c(fw$best, 1L)]
  step_alt <- .spcq_spider_geometry(0, g_lag, curv, xpt[, -1L, drop = FALSE], xl, xu, fw$radius)
  sigma_alt <- .spcq_determinants(md, xb + step_alt, k_new)
  if (abs(sigma_alt) > abs(sigma)) {
    step <- step_alt
    sigma <- sigma_alt
  }
  if (.spcq_type(fw$pb) %in% c("linearly constrained", "nonlinearly constrained")) {
    L <- .spcq_linearizations(fw, xb)
    tol_bd <- .spcq_arrays_tol(xl, xu)
    tol_ub <- .spcq_arrays_tol(L$bub)
    free_xl <- xl <= -tol_bd
    free_xu <- xu >= tol_bd
    free_ub <- L$bub >= tol_ub
    qa <- .spcq_qr_tangential(L$aub, L$aeq, free_xl, free_xu, free_ub)
    g_proj <- .spcq_proj(qa$q, qa$n_act, g_lag)
    norm_g <- .spcq_norm(g_proj)
    if (qa$n_act > 0L && qa$n_act < fw$pb$n && norm_g > .spcq_tiny * fw$radius) {
      step_alt <- fw$radius / norm_g * g_proj
      if (curv(step_alt) < 0) step_alt <- -step_alt
      cbd <- c(xl - step_alt, step_alt - xu)
      cub <- .spcq_mv(L$aub, step_alt) - L$bub
      ceq <- .spcq_mv(L$aeq, step_alt) - L$beq
      maxcv_val <- max(max(c(0, cbd)), max(c(0, cub)), max(c(0, abs(ceq))))
      tol <- max(c(0, abs(step_alt[!free_xl])))
      tol <- max(c(tol, abs(step_alt[!free_xu])))
      tol <- max(c(tol, abs(.spcq_mv(L$aub[!free_ub, , drop = FALSE], step_alt))))
      tol <- .spcq_min(10 * tol, 0.01 * .spcq_norm(step_alt))
      if (maxcv_val <= tol) {
        sigma_alt <- .spcq_determinants(md, xb + step_alt, k_new)
        if (abs(sigma_alt) >= 0.1 * abs(sigma)) step <- .spcq_clip(step_alt, xl, xu)
      }
    }
  }
  step
}

.spcq_soc_step <- function(fw, step) {
  xb <- .spcq_x_best(fw)
  L <- .spcq_linearizations(fw, xb)
  .spcq_normal(L$aub, L$bub, L$aeq, L$beq, fw$pb$bounds$xl - xb, fw$pb$bounds$xu - xb, .spcq_norm(step),
               fw$cst$improve_tcg)
}

.spcq_best_vals <- function(fw) {
  md <- fw$md
  list(fun = md$fun_val[fw$best], cub = md$cub_val[fw$best, ], ceq = md$ceq_val[fw$best, ])
}

.spcq_reduction_ratio <- function(fw, step, fun_val, cub_val, ceq_val) {
  xb <- .spcq_x_best(fw)
  b <- .spcq_best_vals(fw)
  md <- fw$md
  merit_old <- .spcq_merit(fw, xb, b$fun, b$cub, b$ceq)
  merit_new <- .spcq_merit(fw, xb + step, fun_val, cub_val, ceq_val)
  merit_model_old <- .spcq_merit(fw, xb, 0, .spcq_cub(md, xb), .spcq_ceq(md, xb))
  sqp_cub <- .spcq_cub(md, xb) + .spcq_mv(.spcq_cub_grad(md, xb), step)
  sqp_ceq <- .spcq_ceq(md, xb) + .spcq_mv(.spcq_ceq_grad(md, xb), step)
  merit_model_new <- .spcq_merit(fw, xb + step, .spcq_sqp_fun(fw, step), sqp_cub, sqp_ceq)
  if (isTRUE(abs(merit_model_old - merit_model_new) > .spcq_tiny * abs(merit_old - merit_new))) {
    (merit_old - merit_new) / abs(merit_model_old - merit_model_new)
  } else {
    -1
  }
}

.spcq_increase_penalty <- function(fw, step) {
  L <- .spcq_linearizations(fw, .spcq_x_best(fw))
  viol_diff <- .spcq_max(.spcq_norm(c(pmax(0, -L$bub), L$beq)) -
                           .spcq_norm(c(pmax(0, .spcq_mv(L$aub, step) - L$bub), .spcq_mv(L$aeq, step) - L$beq)), 0)
  sqp_val <- .spcq_sqp_fun(fw, step)
  threshold <- .spcq_norm(c(fw$lm_lin_ub, fw$lm_lin_eq, fw$lm_nl_ub, fw$lm_nl_eq))
  if (abs(viol_diff) > .spcq_tiny * abs(sqp_val)) threshold <- .spcq_max(threshold, sqp_val / viol_diff)
  best_save <- fw$best
  if (fw$penalty <= fw$cst$penalty_increase_threshold * threshold) {
    fw$penalty <- .spcq_max(fw$cst$penalty_increase_factor * threshold, 1)
    .spcq_set_best_index(fw)
  }
  best_save == fw$best
}

.spcq_set_best_index <- function(fw) {
  md <- fw$md
  pb <- fw$pb
  best <- fw$best
  xb <- .spcq_x_best(fw)
  m_best <- .spcq_merit(fw, xb, md$fun_val[best], md$cub_val[best, ], md$ceq_val[best, ])
  r_best <- .spcq_maxcv(pb, xb)
  tol <- 10 * .spcq_eps * max(nrow(md$it$xpt), ncol(md$it$xpt)) * .spcq_max(abs(m_best), 1)
  for (k in seq_len(ncol(md$it$xpt))) {
    if (k != fw$best) {
      x_val <- .spcq_point(md$it, k)
      m_val <- .spcq_merit(fw, x_val, md$fun_val[k], md$cub_val[k, ], md$ceq_val[k, ])
      r_val <- .spcq_maxcv(pb, x_val)
      if (isTRUE(m_val < m_best) || (isTRUE(m_val < m_best + tol) && isTRUE(r_val < r_best))) {
        best <- k
        m_best <- m_val
        r_best <- r_val
      }
    }
  }
  fw$best <- best
}

.spcq_index_to_remove <- function(fw, x_new = NULL) {
  it <- fw$md$it
  dd <- it$xpt - it$xpt[, fw$best]
  dist_sq <- colSums(dd^2)
  if (is.null(x_new)) {
    sigma <- 1
    weights <- dist_sq
  } else {
    sigma <- .spcq_determinants(fw$md, x_new)
    weights <- pmax(1, dist_sq / .spcq_max(fw$cst$low_radius_factor * fw$radius, fw$resolution)^2)^3
    weights[fw$best] <- -1
  }
  k_max <- .spcq_argmax(weights * abs(sigma))
  list(k = k_max, dist = sqrt(dist_sq[k_max]))
}

.spcq_update_radius <- function(fw, step, ratio) {
  s_norm <- .spcq_norm(step)
  cst <- fw$cst
  if (isTRUE(ratio <= cst$low_ratio)) {
    .spcq_set_radius(fw, fw$radius * cst$decrease_radius_factor)
  } else if (isTRUE(ratio <= cst$high_ratio)) {
    .spcq_set_radius(fw, .spcq_max(cst$decrease_radius_factor * fw$radius, s_norm))
  } else {
    .spcq_set_radius(fw, .spcq_min(cst$increase_radius_factor * fw$radius,
                                   .spcq_max(cst$decrease_radius_factor * fw$radius,
                                             cst$increase_radius_threshold * s_norm)))
  }
}

.spcq_enhance_resolution <- function(fw, opt) {
  cst <- fw$cst
  if (cst$large_resolution_threshold * opt$radius_final < fw$resolution) {
    fw$resolution <- fw$resolution * cst$decrease_resolution_factor
  } else if (cst$moderate_resolution_threshold * opt$radius_final < fw$resolution) {
    fw$resolution <- sqrt(fw$resolution * opt$radius_final)
  } else {
    fw$resolution <- opt$radius_final
  }
  fw$radius <- .spcq_max(cst$decrease_radius_factor * fw$radius, fw$resolution)
}

.spcq_set_multipliers <- function(fw, x) {
  pb <- fw$pb
  md <- fw$md
  lin <- pb$linear
  n <- pb$n
  incl_lin_ub <- .spcq_mv(lin$a_ub, x) >= lin$b_ub
  incl_nl_ub <- md$cub_val[fw$best, ] >= 0
  incl_xl <- pb$bounds$xl >= x
  incl_xu <- pb$bounds$xu <= x
  m_lin_ub <- sum(incl_lin_ub); m_nl_ub <- sum(incl_nl_ub)
  m_xl <- sum(incl_xl); m_xu <- sum(incl_xu)
  m_lin_eq <- nrow(lin$a_eq); m_nl_eq <- ncol(md$ceq_val)
  if (m_lin_ub + m_nl_ub + m_lin_eq + m_nl_eq > 0) {
    c_jac <- rbind(-.spcq_eye_rows(n, incl_xl), .spcq_eye_rows(n, incl_xu), lin$a_ub[incl_lin_ub, , drop = FALSE],
                   .spcq_cub_grad(md, x, incl_nl_ub), lin$a_eq, .spcq_ceq_grad(md, x))
    g_best <- .spcq_fun_grad(fw, x)
    xl_lm <- rep(-Inf, nrow(c_jac))
    xl_lm[seq_len(m_xl + m_xu + m_lin_ub + m_nl_ub)] <- 0
    res <- .spcq_lsq_linear_bvls(t(c_jac), -g_best, xl_lm, Inf)
    o <- m_xl + m_xu
    fw$lm_lin_ub[incl_lin_ub] <- res[o + seq_len(m_lin_ub)]
    fw$lm_lin_ub[!incl_lin_ub] <- 0
    o <- o + m_lin_ub
    fw$lm_nl_ub[incl_nl_ub] <- res[o + seq_len(m_nl_ub)]
    fw$lm_nl_ub[!incl_nl_ub] <- 0
    o <- o + m_nl_ub
    fw$lm_lin_eq[] <- res[o + seq_len(m_lin_eq)]
    o <- o + m_lin_eq
    fw$lm_nl_eq[] <- res[o + seq_len(m_nl_eq)]
  }
}

.spcq_low_penalty <- function(fw) {
  md <- fw$md
  lin <- fw$pb$linear
  P <- t(md$it$x_base + md$it$xpt)
  r_ub <- cbind(sweep(P %*% t(lin$a_ub), 2L, lin$b_ub), md$cub_val)
  r_eq <- sweep(P %*% t(lin$a_eq), 2L, lin$b_eq)
  r_val <- cbind(r_ub, r_eq, -r_eq, md$ceq_val, -md$ceq_val)
  nanext <- function(v, f) if (all(is.na(v))) NaN else f(v, na.rm = TRUE)
  c_min <- apply(r_val, 2L, nanext, f = min)
  c_max <- apply(r_val, 2L, nanext, f = max)
  if (ncol(r_val) == 0L) c_min <- c_max <- numeric()
  idx <- (c_min < fw$cst$threshold_ratio_constraints * c_max) %in% TRUE
  if (any(idx)) {
    f_min <- nanext(md$fun_val, min)
    f_max <- nanext(md$fun_val, max)
    c_min_neg <- pmin(0, c_min[idx])
    c_diff <- min(c_max[idx] - c_min_neg)
    if (c_diff > .spcq_tiny * (f_max - f_min)) (f_max - f_min) / c_diff else Inf
  } else {
    0
  }
}

.spcq_trust_region <- function(pb, opt, cst) {
  fw <- new.env(parent = emptyenv())
  fw$pb <- pb
  fw$cst <- cst
  fw$penalty <- 0
  fw$md <- .spcq_models(pb, opt, 0)
  fw$best <- 1L
  .spcq_set_best_index(fw)
  fw$lm_lin_ub <- numeric(nrow(pb$linear$a_ub))
  fw$lm_lin_eq <- numeric(nrow(pb$linear$a_eq))
  fw$lm_nl_ub <- numeric(ncol(fw$md$cub_val))
  fw$lm_nl_eq <- numeric(ncol(fw$md$ceq_val))
  .spcq_set_multipliers(fw, .spcq_x_best(fw))
  fw$resolution <- opt$radius_init
  fw$radius <- fw$resolution
  fw
}

# ---------------------------------------------------------------------------
# minimize (main.py) and the SciPy wrapper

.spcq_print_step <- function(message, pb, x, fun_val, r_val, n_eval, n_iter) {
  cat("\n", message, ".\n", sep = "")
  cat("Number of function evaluations: ", n_eval, ".\n", sep = "")
  cat("Number of iterations: ", n_iter, ".\n", sep = "")
  if (!is.null(pb$fun)) cat("Least value of fun: ", .spsq_repr(fun_val), ".\n", sep = "")
  cat("Maximum constraint violation: ", .spsq_repr(r_val), ".\n", sep = "")
  cat("Corresponding point: ", .spcq_format_array(x), ".\n", sep = "")
}

.spcq_result <- function(pb, penalty, success, status, n_iter, opt) {
  be <- .spcq_best_eval(pb, penalty)
  success <- success && is.finite(be$fun) && is.finite(be$maxcv)
  if (!(status %in% c(1L, 4L))) success <- success && be$maxcv <= opt$feasibility_tol
  res <- list(message = .spcq_messages[[as.character(status)]], success = success, status = status,
              x = .spcq_build_x(pb, be$x), fun = be$fun, maxcv = be$maxcv, nfev = pb$n_eval, nit = n_iter)
  if (opt$disp) .spcq_print_step(res$message, pb, res$x, res$fun, res$maxcv, res$nfev, res$nit)
  res
}

# the constants of COBYQA (_set_default_constants with the defaults)
.spcq_constants <- list(decrease_radius_factor = 0.5, increase_radius_factor = sqrt(2), increase_radius_threshold = 2,
                        decrease_radius_threshold = 1.4, decrease_resolution_factor = 0.1,
                        large_resolution_threshold = 250, moderate_resolution_threshold = 16, low_ratio = 0.1,
                        high_ratio = 0.7, very_low_ratio = 0.01, penalty_increase_threshold = 1.5,
                        penalty_increase_factor = 2, short_step_threshold = 0.5, low_radius_factor = 0.1,
                        byrd_omojokun_factor = 0.8, threshold_ratio_constraints = 2, large_shift_factor = 10,
                        large_gradient_factor = 10, resolution_factor = 2, improve_tcg = TRUE)

# _minimize_cobyqa of SciPy: the options of minimize(method = "COBYQA");
# bounds is list(lb, ub) or NULL, constraints a list of
# .sp_nonlinear_constraint() and list(type = "linear", A, lb, ub)
.sp_cobyqa <- function(fun, x0, callback = NULL, disp = FALSE, maxfev = NULL, maxiter = NULL, f_target = -Inf,
                       feasibility_tol = 1e-8, initial_tr_radius = 1.0, final_tr_radius = 1e-6, scale = FALSE,
                       bounds = NULL, constraints = list(), .upstream_bugs = FALSE, ...) {
  extra <- list(...)
  if (length(extra)) warning("Unknown solver options: ", paste(names(extra), collapse = ", "), call. = FALSE)
  x0 <- as.numeric(x0)
  n0 <- length(x0)
  opt <- new.env(parent = emptyenv())
  opt$disp <- isTRUE(as.logical(disp))
  opt$maxfev <- if (is.null(maxfev)) 500L * n0 else as.integer(maxfev)
  opt$maxiter <- if (is.null(maxiter)) 1000L * n0 else as.integer(maxiter)
  opt$target <- as.numeric(f_target)
  opt$feasibility_tol <- as.numeric(feasibility_tol)
  opt$radius_init <- as.numeric(initial_tr_radius)
  opt$radius_final <- as.numeric(final_tr_radius)
  opt$scale <- isTRUE(as.logical(scale))
  xl <- if (is.null(bounds)) rep(-Inf, n0) else rep_len(as.numeric(bounds$lb), n0)
  xu <- if (is.null(bounds)) rep(Inf, n0) else rep_len(as.numeric(bounds$ub), n0)
  lin_cons <- list()
  nl_cons <- list()
  for (con in constraints) {
    if (identical(con$type, "linear")) {
      A <- con$A
      if (!is.matrix(A)) A <- matrix(as.numeric(A), nrow = 1L)
      m <- nrow(A)
      lin_cons[[length(lin_cons) + 1L]] <- list(A = A, lb = rep_len(as.numeric(con$lb), m),
                                                ub = rep_len(as.numeric(con$ub), m))
    } else {
      len <- max(length(con$lb), length(con$ub))
      nl_cons[[length(nl_cons) + 1L]] <- list(fun = con$fun, lb = rep_len(as.numeric(con$lb), len),
                                              ub = rep_len(as.numeric(con$ub), len))
    }
  }
  pb <- .spcq_problem(fun, x0, xl, xu, lin_cons, nl_cons, callback, opt$feasibility_tol, opt$scale, opt$disp,
                      isTRUE(.upstream_bugs))
  # _set_default_options
  if (opt$radius_init <= 0) stop("The initial trust-region radius must be positive.", call. = FALSE)
  if (opt$radius_final < 0) stop("The final trust-region radius must be nonnegative.", call. = FALSE)
  if (opt$radius_init < opt$radius_final)
    stop("The initial trust-region radius must be greater than or equal to the final trust-region radius.",
         call. = FALSE)
  pb$n <- length(pb$x0)
  opt$npt <- 2L * pb$n + 1L
  if (opt$maxfev <= 0) stop("The maximum number of function evaluations must be positive.", call. = FALSE)
  if (opt$maxiter <= 0) stop("The maximum number of iterations must be positive.", call. = FALSE)
  cst <- .spcq_constants
  if (!pb$bounds$is_feasible) return(.spcq_result(pb, 0, FALSE, -1L, 0L, opt))
  if (pb$n == 0L) return(.spcq_result(pb, 0, TRUE, 2L, 0L, opt))
  if (opt$disp) {
    cat("Starting the optimization procedure.\n")
    cat("Initial trust-region radius: ", .spsq_repr(opt$radius_init), ".\n", sep = "")
    cat("Final trust-region radius: ", .spsq_repr(opt$radius_final), ".\n", sep = "")
    cat("Maximum number of function evaluations: ", opt$maxfev, ".\n", sep = "")
    cat("Maximum number of iterations: ", opt$maxiter, ".\n", sep = "")
    cat("\n")
  }
  init_status <- function(status, success) function(e) list(status = status, success = success)
  fw <- tryCatch(.spcq_trust_region(pb, opt, cst),
                 spcq_target = init_status(1L, TRUE), spcq_callback = init_status(3L, TRUE),
                 spcq_feasible = init_status(4L, TRUE), spcq_maxeval = init_status(6L, FALSE),
                 spcq_linalg = init_status(-2L, FALSE))
  if (!is.environment(fw)) return(.spcq_result(pb, 0, fw$success, fw$status, 0L, opt))
  # _eval: a new point, NULL with the exit when it ends the run
  st <- new.env(parent = emptyenv())
  st$success <- FALSE
  eval_step <- function(step) {
    tryCatch({
      if (pb$n_eval >= opt$maxfev) .spcq_raise("spcq_maxeval")
      x_eval <- .spcq_x_best(fw) + step
      r <- .spcq_eval(pb, x_eval, fw$penalty)
      r_val <- .spcq_maxcv(pb, x_eval)
      if (r$fun <= opt$target && r_val <= opt$feasibility_tol) .spcq_raise("spcq_target")
      if (is.null(pb$fun) && r_val <= opt$feasibility_tol) .spcq_raise("spcq_feasible")
      r
    }, spcq_target = function(e) { st$status <- 1L; st$success <- TRUE; NULL },
       spcq_feasible = function(e) { st$status <- 4L; st$success <- TRUE; NULL },
       spcq_callback = function(e) { st$status <- 3L; st$success <- TRUE; NULL },
       spcq_maxeval = function(e) { st$status <- 5L; NULL })
  }
  linalg <- function(expr) tryCatch(expr, spcq_linalg = function(e) { st$status <- -2L; NULL })
  n_iter <- 0L
  k_new <- NULL
  n_short <- 0L
  n_very_short <- 0L
  n_alt <- 0L
  repeat {
    if (n_iter >= opt$maxiter) {
      st$status <- 6L
      break
    }
    n_iter <- n_iter + 1L
    if (.spcq_norm(.spcq_x_best(fw) - fw$md$it$x_base) >= cst$large_shift_factor * fw$radius)
      .spcq_shift_x_base(fw$md, .spcq_x_best(fw))
    radius_save <- fw$radius
    ts <- .spcq_trust_region_step(fw)
    step <- ts$normal + ts$tangential
    s_norm <- .spcq_norm(step)
    if (s_norm <= cst$short_step_threshold * fw$resolution) {
      .spcq_set_radius(fw, fw$radius * cst$decrease_resolution_factor)
      if (radius_save > fw$resolution) {
        n_short <- 0L
        n_very_short <- 0L
      } else {
        n_short <- n_short + 1L
        n_very_short <- n_very_short + 1L
        if (s_norm > 0.1 * fw$resolution) n_very_short <- 0L
      }
      enhance <- n_short >= 5L || n_very_short >= 3L
      if (enhance) {
        n_short <- 0L
        n_very_short <- 0L
        improve_geometry <- FALSE
      } else {
        rm <- linalg(.spcq_index_to_remove(fw))
        if (is.null(rm)) break
        k_new <- rm$k
        improve_geometry <- rm$dist > .spcq_max(fw$radius, cst$resolution_factor * fw$resolution)
      }
    } else {
      same_best <- .spcq_increase_penalty(fw, step)
      if (same_best) {
        r <- eval_step(step)
        if (is.null(r)) break
        b <- .spcq_best_vals(fw)
        xb <- .spcq_x_best(fw)
        merit_old <- .spcq_merit(fw, xb, b$fun, b$cub, b$ceq)
        merit_new <- .spcq_merit(fw, xb + step, r$fun, r$cub, r$ceq)
        if (.spcq_type(pb) == "nonlinearly constrained" && isTRUE(merit_new > merit_old) &&
            .spcq_norm(ts$normal) > cst$byrd_omojokun_factor^2 * fw$radius) {
          soc <- .spcq_soc_step(fw, step)
          if (.spcq_norm(soc) > 0) {
            step <- step + soc
            r <- eval_step(step)
            if (is.null(r)) break
          }
        }
        ratio <- .spcq_reduction_ratio(fw, step, r$fun, r$cub, r$ceq)
        rm <- linalg(.spcq_index_to_remove(fw, .spcq_x_best(fw) + step))
        if (is.null(rm)) break
        k_new <- rm$k
        ill <- linalg(.spcq_update_interpolation(fw$md, k_new, .spcq_x_best(fw) + step, r$fun, r$cub, r$ceq,
                                                 pb$upstream_bugs))
        if (is.null(ill)) break
        .spcq_set_best_index(fw)
        .spcq_update_radius(fw, step, ratio)
        if (fw$radius <= fw$resolution) {
          if (isTRUE(ratio >= cst$very_low_ratio)) {
            n_alt <- 0L
          } else {
            n_alt <- n_alt + 1L
            grad <- .spcq_fun_grad(fw, .spcq_x_best(fw))
            grad_alt <- linalg(.spcq_quad_grad(.spcq_quad_new(fw$md$it, fw$md$fun_val), .spcq_x_best(fw), fw$md$it))
            if (is.null(grad_alt)) break
            if (.spcq_norm(grad) < cst$large_gradient_factor * .spcq_norm(grad_alt)) n_alt <- 0L
            if (n_alt >= 3L) {
              if (is.null(linalg({ .spcq_reset_models(fw$md); TRUE }))) break
              n_alt <- 0L
            }
          }
        }
        .spcq_set_multipliers(fw, .spcq_x_best(fw) + step)
        rm <- linalg(.spcq_index_to_remove(fw))
        if (is.null(rm)) break
        k_new <- rm$k
        improve_geometry <- ill || (isTRUE(ratio <= cst$low_ratio) &&
                                      rm$dist > .spcq_max(fw$radius, cst$resolution_factor * fw$resolution))
        enhance <- radius_save <= fw$resolution && isTRUE(ratio <= cst$low_ratio) && !improve_geometry
      } else {
        enhance <- FALSE
        improve_geometry <- FALSE
      }
    }
    if (enhance) {
      if (fw$resolution <= opt$radius_final) {
        st$success <- TRUE
        st$status <- 0L
        break
      }
      .spcq_enhance_resolution(fw, opt)
      fw$penalty <- .spcq_min(fw$penalty, .spcq_low_penalty(fw))
      .spcq_set_best_index(fw)
      if (opt$disp) {
        xb <- .spcq_x_best(fw)
        .spcq_print_step(paste0("New trust-region radius: ", .spsq_repr(fw$resolution)), pb, .spcq_build_x(pb, xb),
                         .spcq_best_vals(fw)$fun, .spcq_maxcv(pb, xb), pb$n_eval, n_iter)
        cat("\n")
      }
    }
    if (improve_geometry) {
      step <- linalg(.spcq_geometry_step(fw, k_new))
      if (is.null(step)) break
      r <- eval_step(step)
      if (is.null(r)) break
      if (is.null(linalg(.spcq_update_interpolation(fw$md, k_new, .spcq_x_best(fw) + step, r$fun, r$cub, r$ceq,
                                                    pb$upstream_bugs)))) break
      .spcq_set_best_index(fw)
    }
  }
  .spcq_result(pb, fw$penalty, st$success, st$status, n_iter, opt)
}
