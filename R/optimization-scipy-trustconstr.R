# trust-constr, the trust region method for constrained minimization of
# SciPy 1.18.1 (scipy/optimize/_trustregion_constr: minimize_trustregion_constr.py,
# equality_constrained_sqp.py, tr_interior_point.py, qp_subproblem.py,
# projections.py, canonical_constraint.py and report.py), with the pieces of
# scipy/optimize it builds on: VectorFunction, LinearVectorFunction and
# IdentityVectorFunction (_differentiable_functions.py), PreparedConstraint and
# strict_bounds (_constraints.py). The objective is the ScalarFunction of
# optimization-scipy-common.R.
#
# Differences from SciPy:
# - the matrices are dense. SciPy stores the Jacobian of the bounds (and of
#   every constraint when sparse_jacobian = TRUE) as a sparse matrix and
#   factorizes the augmented system [I A'; A 0] with SuperLU. Here the same
#   system is solved by a dense LU factorization (LAPACK dgesv), with the same
#   iterative refinement. "NormalEquation" needs scikit-sparse in SciPy; without
#   it SciPy warns and uses the augmented system, and so does this port.
# - the QR factorization with column pivoting, the SVD and the triangular
#   solves run in R's LAPACK, and the products are sequential sums, while the
#   arm64 SciPy wheel uses Accelerate (the BFGS and SR1 updates model its
#   dsyr and dsymv, see .sp_hessian_update). The iterates agree with SciPy to
#   rounding, not bit for bit, and where a last bit difference flips a branch
#   (a step that lands on a bound, a zero quasi Newton difference) the
#   iteration counts differ while the solution agrees.
# - warnings are signaled each time (Python shows a repeated warning once).

.sptc_fd <- c("2-point", "3-point", "cs")
.sptc_fd_names <- "('2-point', '3-point', 'cs')"

.sptc_messages <- c("0" = "The maximum number of function evaluations is exceeded.",
                    "1" = "`gtol` termination condition is satisfied.",
                    "2" = "`xtol` termination condition is satisfied.",
                    "3" = "`callback` raised `StopIteration`.",
                    "4" = "Constraint violation exceeds 'gtol'")

.sptc_dot <- function(a, b) sum(a * b)

.sptc_norm <- function(a) sqrt(sum(a * a))

# A x and A' v by rows and columns
.sptc_mv <- function(A, x) {
  if (nrow(A) == 0L) return(numeric())
  vapply(seq_len(nrow(A)), function(i) sum(A[i, ] * x), 0)
}

.sptc_tmv <- function(A, v) {
  if (ncol(A) == 0L) return(numeric())
  vapply(seq_len(ncol(A)), function(j) sum(A[, j] * v), 0)
}

# the Python built-ins max(a, b) and min(a, b), NaN included
.sptc_max2 <- function(a, b) if (isTRUE(b > a)) b else a

.sptc_min2 <- function(a, b) if (isTRUE(b < a)) b else a

.sptc_true <- function(v) any(v %in% TRUE)

# H.dot(p) for a Hessian: a matrix, a quasi Newton update, an operator (a
# function of p) or NULL for the zero Hessian of a linear function
.sptc_hmv <- function(H, p) {
  if (is.null(H)) return(numeric(length(p)))
  if (is.function(H)) return(H(p))
  if (inherits(H, "SpHessianUpdate")) return(as.numeric(H$dot(p)))
  .sptc_mv(H, p)
}

# ---------------------------------------------------------------------------
# Constraint functions

# VectorFunction: memoized constraint values, Jacobian and Hessian of the
# Lagrangian part, with SciPy's evaluation counts
.sptc_vector_function <- function(fun, x0, jac, hess, rel_step = NULL, lb = -Inf, ub = Inf,
                                  sparse_jacobian = NULL) {
  fdj <- is.character(jac)
  if (!is.function(jac) && !(fdj && jac %in% .sptc_fd))
    stop("`jac` must be either callable or one of ", .sptc_fd_names, ".", call. = FALSE)
  quasi <- inherits(hess, "SpHessianUpdate")
  if (!(is.function(hess) || quasi || (is.character(hess) && hess %in% .sptc_fd)))
    stop("`hess` must be either callable,HessianUpdateStrategy or one of ", .sptc_fd_names, ".", call. = FALSE)
  fdh <- is.character(hess)
  if (fdj && fdh)
    stop("Whenever the Jacobian is estimated via finite-differences, we require the Hessian to be estimated ",
         "using one of the quasi-Newton strategies.", call. = FALSE)
  vf <- new.env(parent = emptyenv())
  vf$x <- as.numeric(x0)
  vf$n <- length(vf$x)
  vf$nfev0 <- 0L; vf$jnfev <- 0L; vf$njev0 <- 0L; vf$nhev0 <- 0L; vf$hnjev <- 0L
  vf$f_updated <- FALSE; vf$J_updated <- FALSE; vf$H_updated <- FALSE
  fw <- function(x) as.numeric(fun(x))
  as_2d <- function(J) if (is.matrix(J)) J else matrix(as.numeric(J), nrow = 1L)
  fd_jac <- function(x, f0) {
    r <- .sp_approx_derivative(fw, x, method = jac, rel_step = rel_step, f0 = f0, lb = lb, ub = ub)
    list(J = matrix(r$J, nrow = length(f0)), nfev = r$nfev)
  }
  update_fun <- function() {
    if (!vf$f_updated) {
      vf$f <- fw(vf$x)
      vf$nfev0 <- vf$nfev0 + 1L
      vf$f_updated <- TRUE
    }
  }
  update_jac <- function() {
    if (!vf$J_updated) {
      if (fdj) {
        update_fun()
        r <- fd_jac(vf$x, vf$f)
        vf$jnfev <- vf$jnfev + r$nfev
        vf$J <- r$J
      } else {
        vf$njev0 <- vf$njev0 + 1L
        vf$J <- as_2d(jac(vf$x))
      }
      vf$J_updated <- TRUE
    }
  }
  # the finite difference Hessian of v'f at x: the linear operator of
  # approx_derivative(jac_dot_v, x, as_linear_operator = TRUE); every product
  # evaluates the Jacobian (once, twice for "3-point")
  fd_hess <- function(x, v, J0) {
    f0 <- .sptc_tmv(J0, v)
    h <- if (is.null(rel_step)) .sp_eps_for_method(hess) else rel_step
    jdv <- function(z) {
      vf$hnjev <- vf$hnjev + 1L
      J <- jac(z)
      if (!is.matrix(J)) J <- matrix(J, nrow = 1L)
      if (is.complex(J)) colSums(J * v) else .sptc_tmv(J, v)
    }
    function(p) {
      if (all(p == 0)) return(numeric(length(f0)))
      pn <- .sptc_norm(p)
      if (hess == "2-point") {
        dx <- h / pn
        (jdv(x + dx * p) - f0) / dx
      } else if (hess == "3-point") {
        dx <- 2 * h / pn
        f1 <- jdv(x - (dx / 2) * p)
        f2 <- jdv(x + (dx / 2) * p)
        (f2 - f1) / dx
      } else {
        dx <- h / pn
        Im(jdv(x + dx * p * 1i)) / dx
      }
    }
  }
  update_hess <- function() {
    if (!vf$H_updated) {
      if (fdh) {
        update_jac()
        vf$H <- fd_hess(vf$x, vf$v, vf$J)
      } else if (quasi) {
        update_jac()
        if (!is.null(vf$x_prev) && !is.null(vf$J_prev))
          hess$update(vf$x - vf$x_prev, .sptc_tmv(vf$J, vf$v) - .sptc_tmv(vf$J_prev, vf$v))
      } else {
        H <- hess(vf$x, vf$v)
        vf$H <- if (is.function(H)) H else as.matrix(H)
        vf$nhev0 <- vf$nhev0 + 1L
      }
      vf$H_updated <- TRUE
    }
  }
  update_x <- function(x) {
    if (!.sp_same_x(x, vf$x)) {
      if (quasi) {
        update_jac()
        vf$x_prev <- vf$x
        vf$J_prev <- vf$J
      }
      vf$x <- as.numeric(x)
      vf$f_updated <- FALSE; vf$J_updated <- FALSE; vf$H_updated <- FALSE
      if (quasi) update_hess()
    }
  }
  update_fun()
  vf$v <- numeric(length(vf$f))
  vf$m <- length(vf$f)
  if (fdj) {
    r <- fd_jac(vf$x, vf$f)
    vf$J <- r$J
    vf$nfev0 <- vf$nfev0 + r$nfev
  } else {
    vf$J <- as_2d(jac(vf$x))
    vf$njev0 <- vf$njev0 + 1L
  }
  vf$J_updated <- TRUE
  vf$sparse <- isTRUE(sparse_jacobian)
  if (fdh) {
    vf$H <- fd_hess(vf$x, vf$v, vf$J)
  } else if (quasi) {
    hess$initialize(vf$n, "hess")
    vf$H <- hess
    vf$x_prev <- NULL; vf$J_prev <- NULL
  } else {
    H <- hess(vf$x, vf$v)
    vf$H <- if (is.function(H)) H else as.matrix(H)
    vf$nhev0 <- vf$nhev0 + 1L
  }
  vf$H_updated <- TRUE
  vf$fun <- function(x) {
    update_x(x)
    update_fun()
    vf$f
  }
  vf$jac <- function(x) {
    update_x(x)
    update_jac()
    vf$J
  }
  vf$hess <- function(x, v) {
    if (!.sp_same_x(v, vf$v)) {
      vf$v <- v
      vf$H_updated <- FALSE
    }
    update_x(x)
    update_hess()
    vf$H
  }
  vf$counts <- function() c(nfev = vf$nfev0 + vf$jnfev, njev = vf$njev0 + vf$hnjev, nhev = vf$nhev0)
  vf
}

# LinearVectorFunction (A x) and IdentityVectorFunction (x, A = I)
.sptc_linear_function <- function(A, x0, sparse, identity = FALSE) {
  lf <- new.env(parent = emptyenv())
  lf$J <- A
  lf$m <- nrow(A); lf$n <- ncol(A)
  lf$x <- as.numeric(x0)
  value <- function(x) if (identity) x else .sptc_mv(lf$J, x)
  lf$f <- value(lf$x)
  lf$f_updated <- TRUE
  lf$v <- numeric(lf$m)
  lf$sparse <- sparse
  update_x <- function(x) {
    if (!.sp_same_x(x, lf$x)) {
      lf$x <- as.numeric(x)
      lf$f_updated <- FALSE
    }
  }
  lf$fun <- function(x) {
    update_x(x)
    if (!lf$f_updated) {
      lf$f <- value(lf$x)
      lf$f_updated <- TRUE
    }
    lf$f
  }
  lf$jac <- function(x) {
    update_x(x)
    lf$J
  }
  lf$hess <- function(x, v) {
    update_x(x)
    lf$v <- v
    NULL
  }
  lf$counts <- NULL
  lf
}

# PreparedConstraint of a NonlinearConstraint, a LinearConstraint or the bounds
.sptc_prepare <- function(con, x0, sparse_jacobian = NULL, lb = -Inf, ub = Inf, bounds = FALSE) {
  n <- length(x0)
  if (bounds) {
    sparse <- is.null(sparse_jacobian) || isTRUE(sparse_jacobian)
    fn <- .sptc_linear_function(diag(1, n), x0, sparse, identity = TRUE)
  } else {
    type <- tolower(con$type %||% "nonlinear")
    if (type == "nonlinear") {
      fn <- .sptc_vector_function(con$fun, x0, con$jac %||% "2-point", con$hess %||% .sp_hessian_update("BFGS"),
                                  con$finite_diff_rel_step, lb, ub, sparse_jacobian)
    } else if (type == "linear") {
      A <- con$A
      if (!is.matrix(A)) A <- matrix(as.numeric(A), nrow = 1L)
      fn <- .sptc_linear_function(A, x0, isTRUE(sparse_jacobian))
    } else {
      stop("`constraint` of an unknown type is passed.", call. = FALSE)
    }
  }
  m <- fn$m
  widen <- function(v) {
    if (!(length(v) %in% c(1L, m))) stop("operands could not be broadcast together with remapped shapes", call. = FALSE)
    rep_len(v, m)
  }
  clb <- widen(as.numeric(con$lb %||% -Inf))
  cub <- widen(as.numeric(con$ub %||% Inf))
  kf <- widen(as.logical(con$keep_feasible %||% FALSE))
  mask <- kf & clb != cub
  if (.sptc_true(fn$f[mask] < clb[mask]) || .sptc_true(fn$f[mask] > cub[mask]))
    stop("`x0` is infeasible with respect to some inequality constraint with `keep_feasible` set to True.",
         call. = FALSE)
  list(fun = fn, lb = clb, ub = cub, keep_feasible = kf)
}

# ---------------------------------------------------------------------------
# Canonical constraints: c_eq(x) = 0 and c_ineq(x) <= 0

.sptc_canonical_empty <- function(n) {
  list(n_eq = 0L, n_ineq = 0L,
       fun = function(x) list(eq = numeric(), ineq = numeric()),
       jac = function(x) list(eq = matrix(0, 0L, n), ineq = matrix(0, 0L, n)),
       hess = function(x, v_eq, v_ineq) NULL,
       keep_feasible = logical())
}

# CanonicalConstraint.from_PreparedConstraint
.sptc_canonical <- function(pc) {
  lb <- pc$lb; ub <- pc$ub; cf <- pc$fun; kf <- pc$keep_feasible
  n <- cf$n
  empty_J <- matrix(0, 0L, n)
  rows <- function(J, i) J[i, , drop = FALSE]
  if (all(lb == -Inf) && all(ub == Inf)) return(.sptc_canonical_empty(n))
  if (all(lb == ub)) {
    return(list(n_eq = length(lb), n_ineq = 0L,
                fun = function(x) list(eq = cf$fun(x) - lb, ineq = numeric()),
                jac = function(x) list(eq = cf$jac(x), ineq = empty_J),
                hess = function(x, v_eq, v_ineq) cf$hess(x, v_eq),
                keep_feasible = logical()))
  }
  if (all(lb == -Inf)) {
    fin <- ub < Inf
    if (all(fin)) {
      return(list(n_eq = 0L, n_ineq = length(ub),
                  fun = function(x) list(eq = numeric(), ineq = cf$fun(x) - ub),
                  jac = function(x) list(eq = empty_J, ineq = cf$jac(x)),
                  hess = function(x, v_eq, v_ineq) cf$hess(x, v_ineq),
                  keep_feasible = kf))
    }
    fin <- which(fin)
    ubf <- ub[fin]
    return(list(n_eq = 0L, n_ineq = length(fin),
                fun = function(x) list(eq = numeric(), ineq = cf$fun(x)[fin] - ubf),
                jac = function(x) list(eq = empty_J, ineq = rows(cf$jac(x), fin)),
                hess = function(x, v_eq, v_ineq) {
                  v <- numeric(cf$m)
                  v[fin] <- v_ineq
                  cf$hess(x, v)
                },
                keep_feasible = kf[fin]))
  }
  if (all(ub == Inf)) {
    fin <- lb > -Inf
    if (all(fin)) {
      return(list(n_eq = 0L, n_ineq = length(lb),
                  fun = function(x) list(eq = numeric(), ineq = lb - cf$fun(x)),
                  jac = function(x) list(eq = empty_J, ineq = -cf$jac(x)),
                  hess = function(x, v_eq, v_ineq) cf$hess(x, -v_ineq),
                  keep_feasible = kf))
    }
    fin <- which(fin)
    lbf <- lb[fin]
    return(list(n_eq = 0L, n_ineq = length(fin),
                fun = function(x) list(eq = numeric(), ineq = lbf - cf$fun(x)[fin]),
                jac = function(x) list(eq = empty_J, ineq = -rows(cf$jac(x), fin)),
                hess = function(x, v_eq, v_ineq) {
                  v <- numeric(cf$m)
                  v[fin] <- -v_ineq
                  cf$hess(x, v)
                },
                keep_feasible = kf[fin]))
  }
  # _interval_to_canonical
  lb_inf <- lb == -Inf
  ub_inf <- ub == Inf
  equal <- lb == ub
  less <- which(lb_inf & !ub_inf)
  greater <- which(ub_inf & !lb_inf)
  interval <- which(!equal & !lb_inf & !ub_inf)
  equal <- which(equal)
  n_less <- length(less); n_greater <- length(greater); n_interval <- length(interval)
  list(n_eq = length(equal), n_ineq = n_less + n_greater + 2L * n_interval,
       fun = function(x) {
         f <- cf$fun(x)
         list(eq = f[equal] - lb[equal],
              ineq = c(f[less] - ub[less], lb[greater] - f[greater], f[interval] - ub[interval],
                       lb[interval] - f[interval]))
       },
       jac = function(x) {
         J <- cf$jac(x)
         il <- rows(J, interval)
         list(eq = rows(J, equal), ineq = rbind(rows(J, less), -rows(J, greater), il, -il))
       },
       hess = function(x, v_eq, v_ineq) {
         v_l <- v_ineq[seq_len(n_less)]
         v_g <- v_ineq[n_less + seq_len(n_greater)]
         v_il <- v_ineq[n_less + n_greater + seq_len(n_interval)]
         v_ig <- v_ineq[n_less + n_greater + n_interval + seq_len(n_interval)]
         v <- numeric(length(lb))
         v[equal] <- v_eq
         v[less] <- v_l
         v[greater] <- -v_g
         v[interval] <- v_il - v_ig
         cf$hess(x, v)
       },
       keep_feasible = c(kf[less], kf[greater], kf[interval], kf[interval]))
}

# CanonicalConstraint.concatenate
.sptc_concatenate <- function(cl, n) {
  part <- function(lst, key) lapply(lst, `[[`, key)
  stack <- function(ms) if (length(ms)) do.call(rbind, ms) else matrix(0, 0L, n)
  list(n_eq = sum(vapply(cl, `[[`, 0L, "n_eq")), n_ineq = sum(vapply(cl, `[[`, 0L, "n_ineq")),
       fun = function(x) {
         r <- lapply(cl, function(cc) cc$fun(x))
         list(eq = as.numeric(unlist(part(r, "eq"))), ineq = as.numeric(unlist(part(r, "ineq"))))
       },
       jac = function(x) {
         r <- lapply(cl, function(cc) cc$jac(x))
         list(eq = stack(part(r, "eq")), ineq = stack(part(r, "ineq")))
       },
       hess = function(x, v_eq, v_ineq) {
         hess_all <- list()
         ie <- 0L; ii <- 0L
         for (cc in cl) {
           hess_all[length(hess_all) + 1L] <- list(cc$hess(x, v_eq[ie + seq_len(cc$n_eq)],
                                                           v_ineq[ii + seq_len(cc$n_ineq)]))
           ie <- ie + cc$n_eq
           ii <- ii + cc$n_ineq
         }
         function(p) {
           r <- numeric(length(p))
           for (h in hess_all) r <- r + .sptc_hmv(h, p)
           r
         }
       },
       keep_feasible = as.logical(unlist(part(cl, "keep_feasible"))))
}

# initial_constraints_as_canonical
.sptc_initial_canonical <- function(n, prepared) {
  c_eq <- list(); c_ineq <- list(); J_eq <- list(); J_ineq <- list()
  add <- function(lst, v) c(lst, list(v))
  for (pc in prepared) {
    f <- pc$fun$f; J <- pc$fun$J; lb <- pc$lb; ub <- pc$ub
    if (all(lb == ub)) {
      c_eq <- add(c_eq, f - lb); J_eq <- add(J_eq, J)
    } else if (all(lb == -Inf)) {
      fin <- ub < Inf
      c_ineq <- add(c_ineq, f[fin] - ub[fin]); J_ineq <- add(J_ineq, J[fin, , drop = FALSE])
    } else if (all(ub == Inf)) {
      fin <- lb > -Inf
      c_ineq <- add(c_ineq, lb[fin] - f[fin]); J_ineq <- add(J_ineq, -J[fin, , drop = FALSE])
    } else {
      lb_inf <- lb == -Inf
      ub_inf <- ub == Inf
      equal <- lb == ub
      less <- lb_inf & !ub_inf
      greater <- ub_inf & !lb_inf
      interval <- !equal & !lb_inf & !ub_inf
      c_eq <- add(c_eq, f[equal] - lb[equal])
      c_ineq <- c(c_ineq, list(f[less] - ub[less], lb[greater] - f[greater], f[interval] - ub[interval],
                               lb[interval] - f[interval]))
      J_eq <- add(J_eq, J[equal, , drop = FALSE])
      J_ineq <- c(J_ineq, list(J[less, , drop = FALSE], -J[greater, , drop = FALSE], J[interval, , drop = FALSE],
                               -J[interval, , drop = FALSE]))
    }
  }
  stack <- function(ms) if (length(ms)) do.call(rbind, ms) else matrix(0, 0L, n)
  list(c_eq = as.numeric(unlist(c_eq)), c_ineq = as.numeric(unlist(c_ineq)), J_eq = stack(J_eq),
       J_ineq = stack(J_ineq))
}

# ---------------------------------------------------------------------------
# Projections onto the null space and the row space of A (projections.py):
# Z, LS and Y as functions of a vector

.sptc_orthogonality <- function(A, g) {
  norm_g <- .sptc_norm(g)
  norm_A <- sqrt(sum(A * A))
  if (norm_g == 0 || norm_A == 0) return(0)
  .sptc_norm(.sptc_mv(A, g)) / (norm_A * norm_g)
}

.sptc_proj_augmented <- function(A, m, n, orth_tol, max_refin, tol) {
  # K = I without constraints: SuperLU solves it exactly
  if (m == 0L) return(list(Z = function(x) x, LS = function(x) numeric(), Y = function(x) numeric(n)))
  K <- rbind(cbind(diag(1, n), t(A)), cbind(A, matrix(0, m, m)))
  solve_K <- function(v) as.numeric(solve(K, v, tol = 0))
  if (inherits(tryCatch(solve_K(numeric(n + m)), error = function(e) e), "error")) {
    warning("Singular Jacobian matrix. Using dense SVD decomposition to perform the factorizations.", call. = FALSE)
    return(.sptc_proj_svd(A, m, n, orth_tol, max_refin, tol))
  }
  idx_x <- seq_len(n)
  list(Z = function(x) {
         v <- c(x, numeric(m))
         lu_sol <- solve_K(v)
         z <- lu_sol[idx_x]
         k <- 0L
         while (isTRUE(.sptc_orthogonality(A, z) > orth_tol)) {
           if (k >= max_refin) break
           lu_sol <- lu_sol + solve_K(v - .sptc_mv(K, lu_sol))
           z <- lu_sol[idx_x]
           k <- k + 1L
         }
         z
       },
       LS = function(x) solve_K(c(x, numeric(m)))[n + seq_len(m)],
       Y = function(x) solve_K(c(numeric(n), x))[idx_x])
}

.sptc_proj_qr <- function(A, m, n, orth_tol, max_refin, tol) {
  q <- qr(t(A), LAPACK = TRUE)
  Q <- qr.Q(q)
  R <- qr.R(q)
  P <- q$pivot
  if (max(abs(R[nrow(R), ])) < tol) {
    warning("Singular Jacobian matrix. Using SVD decomposition to perform the factorizations.", call. = FALSE)
    return(.sptc_proj_svd(A, m, n, orth_tol, max_refin, tol))
  }
  if (nrow(R) != ncol(R))
    stop("The 'expected square matrix' error can occur if there are more equality constraints than independent ",
         "variables. Consider how your constraints are set up, or use factorization_method='SVDFactorization'.",
         call. = FALSE)
  lsq <- function(x) {
    z <- numeric(m)
    z[P] <- backsolve(R, .sptc_tmv(Q, x))
    z
  }
  list(Z = function(x) {
         z <- x - .sptc_tmv(A, lsq(x))
         k <- 0L
         while (isTRUE(.sptc_orthogonality(A, z) > orth_tol)) {
           if (k >= max_refin) break
           z <- z - .sptc_tmv(A, lsq(z))
           k <- k + 1L
         }
         z
       },
       LS = lsq,
       Y = function(x) .sptc_mv(Q, backsolve(R, x[P], transpose = TRUE)))
}

.sptc_proj_svd <- function(A, m, n, orth_tol, max_refin, tol) {
  sv <- svd(A)
  keep <- sv$d > tol
  U <- sv$u[, keep, drop = FALSE]
  Vt <- t(sv$v[, keep, drop = FALSE])
  s <- sv$d[keep]
  lsq <- function(x) .sptc_mv(U, 1 / s * .sptc_mv(Vt, x))
  list(Z = function(x) {
         z <- x - .sptc_tmv(A, lsq(x))
         k <- 0L
         while (isTRUE(.sptc_orthogonality(A, z) > orth_tol)) {
           if (k >= max_refin) break
           z <- z - .sptc_tmv(A, lsq(z))
           k <- k + 1L
         }
         z
       },
       LS = lsq,
       Y = function(x) .sptc_tmv(Vt, 1 / s * .sptc_tmv(U, x)))
}

.sptc_projections <- function(A, sparse, method = NULL, orth_tol = 1e-12, max_refin = 3L, tol = 1e-15) {
  m <- nrow(A); n <- ncol(A)
  if (m * n == 0) sparse <- TRUE
  if (sparse) {
    if (is.null(method)) method <- "AugmentedSystem"
    if (!(method %in% c("NormalEquation", "AugmentedSystem"))) stop("Method not allowed for sparse array.", call. = FALSE)
    if (method == "NormalEquation") {
      warning("Only accepts 'NormalEquation' option when scikit-sparse is available. Using 'AugmentedSystem' ",
              "option instead.", call. = FALSE)
      method <- "AugmentedSystem"
    }
  } else {
    if (is.null(method)) method <- "QRFactorization"
    if (!(method %in% c("QRFactorization", "SVDFactorization"))) stop("Method not allowed for dense array.", call. = FALSE)
  }
  switch(method,
         AugmentedSystem = .sptc_proj_augmented(A, m, n, orth_tol, max_refin, tol),
         QRFactorization = .sptc_proj_qr(A, m, n, orth_tol, max_refin, tol),
         SVDFactorization = .sptc_proj_svd(A, m, n, orth_tol, max_refin, tol))
}

# ---------------------------------------------------------------------------
# Quadratic subproblems (qp_subproblem.py)

.sptc_copysign <- function(a, b) if (b < 0 || (b == 0 && 1 / b < 0)) -abs(a) else abs(a)

.sptc_sphere_intersections <- function(z, d, trust_radius, entire_line = FALSE) {
  if (.sptc_norm(d) == 0) return(list(ta = 0, tb = 0, intersect = FALSE))
  if (is.infinite(trust_radius)) {
    if (entire_line) return(list(ta = -Inf, tb = Inf, intersect = TRUE))
    return(list(ta = 0, tb = 1, intersect = TRUE))
  }
  a <- .sptc_dot(d, d)
  b <- 2 * .sptc_dot(z, d)
  cc <- .sptc_dot(z, z) - trust_radius^2
  discriminant <- b * b - 4 * a * cc
  if (isTRUE(discriminant < 0)) return(list(ta = 0, tb = 0, intersect = FALSE))
  aux <- b + .sptc_copysign(sqrt(discriminant), b)
  ta <- -aux / (2 * a)
  tb <- -2 * cc / aux
  if (isTRUE(tb < ta)) {
    tmp <- ta; ta <- tb; tb <- tmp
  }
  if (entire_line) return(list(ta = ta, tb = tb, intersect = TRUE))
  if (isTRUE(tb < 0) || isTRUE(ta > 1)) return(list(ta = 0, tb = 0, intersect = FALSE))
  list(ta = .sptc_max2(0, ta), tb = .sptc_min2(1, tb), intersect = TRUE)
}

.sptc_box_intersections <- function(z, d, lb, ub, entire_line = FALSE) {
  if (.sptc_norm(d) == 0) return(list(ta = 0, tb = 0, intersect = FALSE))
  zero_d <- d %in% 0
  if (.sptc_true(z[zero_d] < lb[zero_d]) || .sptc_true(z[zero_d] > ub[zero_d]))
    return(list(ta = 0, tb = 0, intersect = FALSE))
  keep <- !zero_d
  z <- z[keep]; d <- d[keep]; lb <- lb[keep]; ub <- ub[keep]
  t_lb <- (lb - z) / d
  t_ub <- (ub - z) / d
  ta <- Reduce(.sptc_max2, pmin(t_lb, t_ub))
  tb <- Reduce(.sptc_min2, pmax(t_lb, t_ub))
  intersect <- isTRUE(ta <= tb)
  if (!entire_line) {
    if (isTRUE(tb < 0) || isTRUE(ta > 1)) return(list(ta = 0, tb = 0, intersect = FALSE))
    ta <- .sptc_max2(0, ta)
    tb <- .sptc_min2(1, tb)
  }
  list(ta = ta, tb = tb, intersect = intersect)
}

.sptc_box_sphere <- function(z, d, lb, ub, trust_radius, entire_line = FALSE) {
  b <- .sptc_box_intersections(z, d, lb, ub, entire_line)
  s <- .sptc_sphere_intersections(z, d, trust_radius, entire_line)
  ta <- max(b$ta, s$ta)
  tb <- min(b$tb, s$tb)
  list(ta = ta, tb = tb, intersect = b$intersect && s$intersect && isTRUE(ta <= tb))
}

.sptc_inside <- function(x, lb, ub) isTRUE(all(lb <= x)) && isTRUE(all(x <= ub))

.sptc_reinforce <- function(x, lb, ub) pmin(pmax(x, lb), ub)

.sptc_modified_dogleg <- function(A, Y, b, trust_radius, lb, ub) {
  newton_point <- -Y(b)
  if (.sptc_inside(newton_point, lb, ub) && isTRUE(.sptc_norm(newton_point) <= trust_radius)) return(newton_point)
  g <- .sptc_tmv(A, b)
  A_g <- .sptc_mv(A, g)
  cauchy_point <- -.sptc_dot(g, g) / .sptc_dot(A_g, A_g) * g
  origin_point <- numeric(length(cauchy_point))
  z <- cauchy_point
  p <- newton_point - cauchy_point
  r <- .sptc_box_sphere(z, p, lb, ub, trust_radius)
  if (r$intersect) {
    x1 <- z + r$tb * p
  } else {
    z <- origin_point
    p <- cauchy_point
    r <- .sptc_box_sphere(z, p, lb, ub, trust_radius)
    x1 <- z + r$tb * p
  }
  z <- origin_point
  p <- newton_point
  r <- .sptc_box_sphere(z, p, lb, ub, trust_radius)
  x2 <- z + r$tb * p
  if (isTRUE(.sptc_norm(.sptc_mv(A, x1) + b) < .sptc_norm(.sptc_mv(A, x2) + b))) x1 else x2
}

.sptc_projected_cg <- function(H, cvec, Z, Y, b, trust_radius = Inf, lb = NULL, ub = NULL) {
  close_to_zero <- 0x1.ef2d0f5da7dd9p-84
  n <- length(cvec); m <- length(b)
  x <- Y(-b)
  r <- Z(H(x) + cvec)
  g <- Z(r)
  p <- -g
  H_p <- H(p)
  rt_g <- .sptc_norm(g)^2
  tr_distance <- trust_radius - .sptc_norm(x)
  if (isTRUE(tr_distance < 0)) stop("Trust region problem does not have a solution.", call. = FALSE)
  if (isTRUE(tr_distance < close_to_zero))
    return(list(x = x, info = list(niter = 0L, stop_cond = 2L, hits_boundary = TRUE)))
  tol <- .sptc_max2(.sptc_min2(0.01 * sqrt(rt_g), 0.1 * rt_g), close_to_zero)
  if (is.null(lb)) lb <- rep(-Inf, n)
  if (is.null(ub)) ub <- rep(Inf, n)
  max_iter <- n - m
  max_infeasible_iter <- n - m
  hits_boundary <- FALSE
  stop_cond <- 1L
  counter <- 0L
  last_feasible_x <- numeric(length(x))
  k <- 0L
  for (i in seq_len(max(0L, max_iter))) {
    if (isTRUE(rt_g < tol)) {
      stop_cond <- 4L
      break
    }
    k <- k + 1L
    pt_H_p <- .sptc_dot(H_p, p)
    if (isTRUE(pt_H_p <= 0)) {
      if (is.infinite(trust_radius)) stop("Negative curvature not allowed for unrestricted problems.", call. = FALSE)
      s <- .sptc_box_sphere(x, p, lb, ub, trust_radius, entire_line = TRUE)
      if (s$intersect) x <- x + s$tb * p
      x <- .sptc_reinforce(x, lb, ub)
      stop_cond <- 3L
      hits_boundary <- TRUE
      break
    }
    alpha <- rt_g / pt_H_p
    x_next <- x + alpha * p
    if (isTRUE(.sptc_norm(x_next) >= trust_radius)) {
      s <- .sptc_box_sphere(x, alpha * p, lb, ub, trust_radius)
      if (s$intersect) x <- x + s$tb * alpha * p
      x <- .sptc_reinforce(x, lb, ub)
      stop_cond <- 2L
      hits_boundary <- TRUE
      break
    }
    counter <- if (.sptc_inside(x_next, lb, ub)) 0L else counter + 1L
    if (counter > 0L) {
      s <- .sptc_box_sphere(x, alpha * p, lb, ub, trust_radius)
      if (s$intersect) {
        last_feasible_x <- .sptc_reinforce(x + s$tb * alpha * p, lb, ub)
        counter <- 0L
      }
    }
    if (counter > max_infeasible_iter) break
    r_next <- r + alpha * H_p
    g_next <- Z(r_next)
    rt_g_next <- .sptc_norm(g_next)^2
    beta <- rt_g_next / rt_g
    p <- -g_next + beta * p
    x <- x_next
    g <- g_next
    r <- g_next
    rt_g <- .sptc_norm(g)^2
    H_p <- H(p)
  }
  if (!.sptc_inside(x, lb, ub)) {
    x <- last_feasible_x
    hits_boundary <- TRUE
  }
  list(x = x, info = list(niter = k, stop_cond = stop_cond, hits_boundary = hits_boundary))
}

# ---------------------------------------------------------------------------
# Byrd-Omojokun SQP for equality constraints (equality_constrained_sqp.py).
# fun_and_constr(x) returns list(f, b, x): the barrier subproblem sets the
# slacks of keep_feasible constraints in place, so the point comes back too.

.sptc_eqsqp <- function(fun_and_constr, grad_and_jac, lagr_hess, x0, fun0, grad0, constr0, jac0, stop_criteria,
                        state, initial_penalty, initial_trust_radius, factorization_method, sparse,
                        trust_lb = NULL, trust_ub = NULL, scaling = NULL) {
  n <- length(x0)
  if (is.null(trust_lb)) trust_lb <- rep(-Inf, n)
  if (is.null(trust_ub)) trust_ub <- rep(Inf, n)
  if (is.null(scaling)) scaling <- function(x) function(d) d
  x <- x0
  trust_radius <- initial_trust_radius
  penalty <- initial_penalty
  f <- fun0; cvec <- grad0; b <- constr0; A <- jac0
  S <- scaling(x)
  pr <- .sptc_projections(A, sparse, factorization_method)
  v <- -pr$LS(cvec)
  H <- lagr_hess(x, v)
  optimality <- max(abs(cvec + .sptc_tmv(A, v)))
  constr_violation <- if (length(b)) max(abs(b)) else 0
  cg_info <- list(niter = 0L, stop_cond = 0L, hits_boundary = FALSE)
  last_iteration_failed <- FALSE
  while (!stop_criteria(state, x, last_iteration_failed, optimality, constr_violation, trust_radius, penalty,
                        cg_info)) {
    dn <- .sptc_modified_dogleg(A, pr$Y, b, 0.8 * trust_radius, 0.5 * trust_lb, 0.5 * trust_ub)
    c_t <- H(dn) + cvec
    b_t <- numeric(length(b))
    trust_radius_t <- sqrt(trust_radius^2 - .sptc_norm(dn)^2)
    cg <- .sptc_projected_cg(H, c_t, pr$Z, pr$Y, b_t, trust_radius_t, trust_lb - dn, trust_ub - dn)
    dt <- cg$x
    cg_info <- cg$info
    d <- dn + dt
    quadratic_model <- 1 / 2 * .sptc_dot(H(d), d) + .sptc_dot(cvec, d)
    linearized_constr <- .sptc_mv(A, d) + b
    vpred <- .sptc_max2(1e-16, .sptc_norm(b) - .sptc_norm(linearized_constr))
    previous_penalty <- penalty
    if (isTRUE(quadratic_model > 0)) {
      new_penalty <- quadratic_model / ((1 - 0.3) * vpred)
      penalty <- .sptc_max2(penalty, new_penalty)
    }
    predicted_reduction <- -quadratic_model + penalty * vpred
    merit_function <- f + penalty * .sptc_norm(b)
    nxt <- fun_and_constr(x + S(d))
    x_next <- nxt$x; f_next <- nxt$f; b_next <- nxt$b
    merit_function_next <- f_next + penalty * .sptc_norm(b_next)
    actual_reduction <- merit_function - merit_function_next
    reduction_ratio <- actual_reduction / predicted_reduction
    if (isTRUE(reduction_ratio < 1e-8) && isTRUE(.sptc_norm(dn) <= 0.1 * .sptc_norm(dt))) {
      y <- -pr$Y(b_next)
      bi <- .sptc_box_intersections(d, y, trust_lb, trust_ub)
      soc <- fun_and_constr(x + S(d + bi$tb * y))
      merit_function_soc <- soc$f + penalty * .sptc_norm(soc$b)
      actual_reduction_soc <- merit_function - merit_function_soc
      reduction_ratio_soc <- actual_reduction_soc / predicted_reduction
      if (bi$intersect && isTRUE(reduction_ratio_soc >= 1e-8)) {
        x_next <- soc$x; f_next <- soc$f; b_next <- soc$b
        reduction_ratio <- reduction_ratio_soc
      }
    }
    if (isTRUE(reduction_ratio >= 0.9)) {
      trust_radius <- .sptc_max2(7.0 * .sptc_norm(d), trust_radius)
    } else if (isTRUE(reduction_ratio >= 0.3)) {
      trust_radius <- .sptc_max2(2.0 * .sptc_norm(d), trust_radius)
    } else if (isTRUE(reduction_ratio < 1e-8)) {
      trust_reduction <- (1 - 1e-8) / (1 - reduction_ratio)
      new_trust_radius <- trust_reduction * .sptc_norm(d)
      if (isTRUE(new_trust_radius >= 0.5 * trust_radius)) {
        trust_radius <- trust_radius * 0.5
      } else if (isTRUE(new_trust_radius >= 0.1 * trust_radius)) {
        trust_radius <- new_trust_radius
      } else {
        trust_radius <- trust_radius * 0.1
      }
    }
    if (isTRUE(reduction_ratio >= 1e-8)) {
      x <- x_next
      f <- f_next; b <- b_next
      gj <- grad_and_jac(x)
      cvec <- gj$g; A <- gj$J
      S <- scaling(x)
      pr <- .sptc_projections(A, sparse, factorization_method)
      v <- -pr$LS(cvec)
      H <- lagr_hess(x, v)
      last_iteration_failed <- FALSE
      optimality <- max(abs(cvec + .sptc_tmv(A, v)))
      constr_violation <- if (length(b)) max(abs(b)) else 0
    } else {
      penalty <- previous_penalty
      last_iteration_failed <- TRUE
    }
  }
  x
}

# ---------------------------------------------------------------------------
# Barrier method for inequality constraints (tr_interior_point.py)

.sptc_tr_interior_point <- function(fun, grad, lagr_hess, n_vars, n_ineq, n_eq, constr, jac, x0, fun0, grad0,
                                    constr_ineq0, jac_ineq0, constr_eq0, jac_eq0, stop_criteria,
                                    enforce_feasibility, xtol, state, initial_barrier_parameter, initial_tolerance,
                                    initial_penalty, initial_trust_radius, factorization_method, fd_lb, fd_ub,
                                    sparse) {
  sb <- new.env(parent = emptyenv())
  sb$mu <- initial_barrier_parameter
  sb$tol <- initial_tolerance
  sb$terminate <- FALSE
  enforce <- if (length(enforce_feasibility)) enforce_feasibility else rep(FALSE, n_ineq)
  iv <- seq_len(n_vars)
  is <- n_vars + seq_len(n_ineq)
  compute_function <- function(f, c_ineq, s) {
    s[enforce] <- -c_ineq[enforce]
    log_s <- rep(-Inf, length(s))
    pos <- which(s > 0)
    log_s[pos] <- log(s[pos])
    list(value = f - sb$mu * sum(log_s), s = s)
  }
  compute_constr <- function(c_ineq, c_eq, s) c(c_eq, c_ineq + s)
  compute_gradient <- function(g) c(g, -sb$mu * rep(1, n_ineq))
  compute_jacobian <- function(J_eq, J_ineq, s) {
    if (n_ineq == 0L) return(J_eq)
    rbind(cbind(J_eq, matrix(0, n_eq, n_ineq)), cbind(J_ineq, diag(s, nrow = n_ineq)))
  }
  function_and_constraints <- function(z) {
    x <- z[iv]
    if (.sptc_true(x < fd_lb | x > fd_ub)) {
      f <- Inf
      c_eq <- numeric(n_eq)
      c_ineq <- numeric(n_ineq)
    } else {
      cc <- constr(x)
      c_eq <- cc$eq
      c_ineq <- cc$ineq
      f <- if (.sptc_true(c_ineq[enforce] >= 0)) Inf else fun(x)
    }
    cf <- compute_function(f, c_ineq, z[is])
    z[is] <- cf$s
    list(f = cf$value, b = compute_constr(c_ineq, c_eq, cf$s), x = z)
  }
  gradient_and_jacobian <- function(z) {
    g <- grad(z[iv])
    J <- jac(z[iv])
    list(g = compute_gradient(g), J = compute_jacobian(J$eq, J$ineq, z[is]))
  }
  lagrangian_hessian <- function(z, v) {
    Hx <- lagr_hess(z[iv], v[seq_len(n_eq)], v[n_eq + seq_len(n_ineq)])
    if (n_ineq > 0L) {
      vi <- v[length(v) - n_ineq + seq_len(n_ineq)]
      S_Hs_S <- ifelse(vi %in% NaN | !(vi > 0), sb$mu, vi * z[is])
    }
    function(vec) if (n_ineq > 0L) c(Hx(vec[iv]), S_Hs_S * vec[is]) else Hx(vec[iv])
  }
  scaling <- function(z) {
    diag_elements <- c(rep(1, n_vars), z[is])
    function(vec) diag_elements * vec
  }
  sub_stop <- function(state, z, last_iteration_failed, optimality, constr_violation, trust_radius, penalty,
                       cg_info) {
    if (stop_criteria(state, z[iv], last_iteration_failed, trust_radius, penalty, cg_info, sb$mu, sb$tol)) {
      sb$terminate <- TRUE
      return(TRUE)
    }
    (isTRUE(optimality < sb$tol) && isTRUE(constr_violation < sb$tol)) || isTRUE(trust_radius < xtol)
  }
  trust_radius <- initial_trust_radius
  s0 <- pmax(-1.5 * constr_ineq0, rep(1, n_ineq))
  cf <- compute_function(fun0, constr_ineq0, s0)
  s0 <- cf$s
  fun0_sub <- cf$value
  constr0_sub <- compute_constr(constr_ineq0, constr_eq0, s0)
  grad0_sub <- compute_gradient(grad0)
  jac0_sub <- compute_jacobian(jac_eq0, jac_ineq0, s0)
  z <- c(x0, s0)
  trust_lb <- c(rep(-Inf, n_vars), rep(-0.995, n_ineq))
  trust_ub <- rep(Inf, n_vars + n_ineq)
  repeat {
    z <- .sptc_eqsqp(function_and_constraints, gradient_and_jacobian, lagrangian_hessian, z, fun0_sub, grad0_sub,
                     constr0_sub, jac0_sub, sub_stop, state, initial_penalty, trust_radius, factorization_method,
                     sparse, trust_lb, trust_ub, scaling)
    if (sb$terminate) break
    trust_radius <- .sptc_max2(initial_trust_radius, 5 * state$tr_radius)
    sb$mu <- sb$mu * 0.2
    sb$tol <- sb$tol * 0.2
    fc <- function_and_constraints(z)
    z <- fc$x
    fun0_sub <- fc$f
    constr0_sub <- fc$b
    gj <- gradient_and_jacobian(z)
    grad0_sub <- gj$g
    jac0_sub <- gj$J
  }
  z[iv]
}

# ---------------------------------------------------------------------------
# Reports (report.py) and the driver (minimize_trustregion_constr.py)

.sptc_center <- function(s, w) {
  pad <- w - nchar(s)
  if (pad <= 0) return(s)
  left <- pad %/% 2
  paste0(strrep(" ", left), s, strrep(" ", pad - left))
}

# Python's format(x, '.<digits>e'), with '+' for the signed column
.sptc_fmt_e <- function(x, digits, plus = FALSE) {
  s <- if (is.nan(x)) "nan" else if (is.infinite(x)) (if (x > 0) "inf" else "-inf") else
    sprintf(paste0("%.", digits, "e"), x)
  if (plus && !startsWith(s, "-")) s <- paste0("+", s)
  s
}

.sptc_report_columns <- list(
  basic = list(names = c("niter", "f evals", "CG iter", "obj func", "tr radius", "opt", "c viol"),
               widths = c(7L, 7L, 7L, 13L, 10L, 10L, 10L)),
  sqp = list(names = c("niter", "f evals", "CG iter", "obj func", "tr radius", "opt", "c viol", "penalty", "CG stop"),
             widths = c(7L, 7L, 7L, 13L, 10L, 10L, 10L, 10L, 7L)),
  ip = list(names = c("niter", "f evals", "CG iter", "obj func", "tr radius", "opt", "c viol", "penalty",
                      "barrier param", "CG stop"),
            widths = c(7L, 7L, 7L, 13L, 10L, 10L, 10L, 10L, 13L, 7L)))

.sptc_print_header <- function(kind) {
  cols <- .sptc_report_columns[[kind]]
  cat("|", paste(mapply(.sptc_center, cols$names, cols$widths), collapse = "|"), "|\n", sep = "")
  cat("|", paste(vapply(cols$widths, function(w) strrep("-", w), ""), collapse = "|"), "|\n", sep = "")
}

.sptc_print_iteration <- function(kind, state) {
  cells <- c(as.character(state$nit), as.character(state$nfev), as.character(state$cg_niter),
             .sptc_fmt_e(state$fun, 4L, plus = TRUE), .sptc_fmt_e(state$tr_radius, 2L),
             .sptc_fmt_e(state$optimality, 2L), .sptc_fmt_e(state$constr_violation, 2L))
  if (kind != "basic") cells <- c(cells, .sptc_fmt_e(state$constr_penalty, 2L))
  if (kind == "ip") cells <- c(cells, .sptc_fmt_e(state$barrier_parameter, 2L))
  if (kind != "basic") cells <- c(cells, as.character(state$cg_stop_cond))
  cat("|", paste(mapply(.sptc_center, cells, .sptc_report_columns[[kind]]$widths), collapse = "|"), "|\n", sep = "")
}

# Python's format(x, '4.2'): the 'g' style with at least one digit after the point
.sptc_fmt_time <- function(x) {
  s <- sprintf("%.2g", x)
  if (!grepl("[.e]", s)) s <- paste0(s, ".0")
  formatC(s, width = 4L)
}

# minimize(method = "trust-constr"); jac is a function or a finite difference
# scheme ("2-point" when NULL), hess a function or a quasi Newton update
# (BFGS when NULL), bounds list(lb, ub, keep_feasible) and constraints a list
# of .sp_nonlinear_constraint() and list(type = "linear", A, lb, ub) entries.
.sp_trust_constr <- function(fun, x0, xtol = 1e-8, gtol = 1e-8, barrier_tol = 1e-8, sparse_jacobian = NULL,
                             callback = NULL, maxiter = 1000L, verbose = 0L, finite_diff_rel_step = NULL,
                             initial_constr_penalty = 1.0, initial_tr_radius = 1.0, initial_barrier_parameter = 0.1,
                             initial_barrier_tolerance = 0.1, factorization_method = NULL, disp = FALSE,
                             workers = NULL, jac = NULL, hess = NULL, hessp = NULL, bounds = NULL,
                             constraints = list(), ...) {
  extra <- list(...)
  if (length(extra))
    stop("_minimize_trustregion_constr() got an unexpected keyword argument '", names(extra)[1L], "'", call. = FALSE)
  x0 <- as.numeric(x0)
  n_vars <- length(x0)
  if (is.null(jac)) jac <- "2-point"
  if (is.null(hess)) {
    # HessianLinearOperator: hess(x) is the operator p -> hessp(x, p)
    hess <- if (is.function(hessp)) function(x) function(p) hessp(x, p) else .sp_hessian_update("BFGS")
  }
  if (isTRUE(disp) && verbose == 0) verbose <- 1L
  fd_lb <- -Inf; fd_ub <- Inf
  if (!is.null(bounds)) {
    lb <- rep_len(as.numeric(bounds$lb), n_vars)
    ub <- rep_len(as.numeric(bounds$ub), n_vars)
    fl <- is.finite(lb); fu <- is.finite(ub)
    lb[fl] <- .sp_nextafter(lb[fl], rep(-Inf, sum(fl)))
    ub[fu] <- .sp_nextafter(ub[fu], rep(Inf, sum(fu)))
    kf <- rep_len(as.logical(bounds$keep_feasible %||% FALSE), n_vars)
    bounds <- list(lb = lb, ub = ub, keep_feasible = kf)
    fd_lb <- ifelse(kf, lb, -Inf)
    fd_ub <- ifelse(kf, ub, Inf)
  }
  objective <- .sp_scalar_function(fun, x0, jac, hess, finite_diff_rel_step, fd_lb, fd_ub)
  prepared <- lapply(constraints, function(con) .sptc_prepare(con, x0, sparse_jacobian, fd_lb, fd_ub))
  if (length(prepared)) sparse_jacobian <- any(vapply(prepared, function(pc) pc$fun$sparse, TRUE))
  if (!is.null(bounds)) {
    if (is.null(sparse_jacobian)) sparse_jacobian <- TRUE
    prepared[[length(prepared) + 1L]] <- .sptc_prepare(bounds, x0, sparse_jacobian, bounds = TRUE)
  }
  sparse <- isTRUE(sparse_jacobian)
  ic <- .sptc_initial_canonical(n_vars, prepared)
  canonical_all <- lapply(prepared, .sptc_canonical)
  canonical <- if (length(canonical_all) == 0L) .sptc_canonical_empty(n_vars) else
    if (length(canonical_all) == 1L) canonical_all[[1L]] else .sptc_concatenate(canonical_all, n_vars)
  lagrangian_hess <- function(x, v_eq, v_ineq = numeric()) {
    H_objective <- objective$hess(x)
    H_constraints <- canonical$hess(x, v_eq, v_ineq)
    function(p) .sptc_hmv(H_objective, p) + .sptc_hmv(H_constraints, p)
  }
  method <- if (canonical$n_ineq == 0L) "equality_constrained_sqp" else "tr_interior_point"
  counts <- function(i) vapply(prepared, function(pc) if (is.null(pc$fun$counts)) 0L else pc$fun$counts()[[i]], 0L)
  state <- new.env(parent = emptyenv())
  state$nit <- 0L; state$nfev <- 0L; state$njev <- 0L; state$nhev <- 0L
  state$cg_niter <- 0L; state$cg_stop_cond <- 0L
  state$fun <- objective$f; state$grad <- objective$g; state$lagrangian_grad <- objective$g
  state$constr <- lapply(prepared, function(pc) pc$fun$f)
  state$jac <- lapply(prepared, function(pc) pc$fun$J)
  state$constr_nfev <- rep(0L, length(prepared)); state$constr_njev <- rep(0L, length(prepared))
  state$constr_nhev <- rep(0L, length(prepared))
  state$v <- lapply(prepared, function(pc) pc$fun$v)
  state$method <- method
  start_time <- proc.time()[["elapsed"]]
  update_state <- function(state, x, last_iteration_failed, tr_radius, constr_penalty, cg_info) {
    state$nit <- state$nit + 1L
    state$nfev <- objective$nfev()
    state$njev <- objective$ngev
    state$nhev <- objective$nhev
    state$constr_nfev <- counts("nfev")
    state$constr_njev <- counts("njev")
    state$constr_nhev <- counts("nhev")
    if (!last_iteration_failed) {
      state$x <- x
      state$fun <- objective$f
      state$grad <- objective$g
      state$v <- lapply(prepared, function(pc) pc$fun$v)
      state$constr <- lapply(prepared, function(pc) pc$fun$f)
      state$jac <- lapply(prepared, function(pc) pc$fun$J)
      lg <- state$grad
      for (pc in prepared) lg <- lg + .sptc_tmv(pc$fun$J, pc$fun$v)
      state$lagrangian_grad <- lg
      state$optimality <- max(abs(lg))
      cv <- 0
      for (i in seq_along(prepared)) cv <- max(cv, max(prepared[[i]]$lb - state$constr[[i]]),
                                               max(state$constr[[i]] - prepared[[i]]$ub))
      state$constr_violation <- cv
    }
    state$execution_time <- proc.time()[["elapsed"]] - start_time
    state$tr_radius <- tr_radius
    state$constr_penalty <- constr_penalty
    state$cg_niter <- state$cg_niter + cg_info$niter
    state$cg_stop_cond <- cg_info$stop_cond
  }
  run_callback <- function(state) {
    if (!is.function(callback)) return(FALSE)
    snap <- as.list(state)
    out <- tryCatch(if (identical(names(formals(callback)), "intermediate_result"))
                      callback(intermediate_result = snap) else callback(state$x, snap),
                    StopIteration = function(e) TRUE)
    isTRUE(as.logical(out)[1L])
  }
  finish_check <- function(state, xtol_met) {
    state$status <- NULL
    state$niter <- state$nit
    if (run_callback(state)) {
      state$status <- 3L
      return(TRUE)
    }
    if (isTRUE(state$optimality < gtol) && isTRUE(state$constr_violation < gtol)) state$status <- 1L
    else if (xtol_met) state$status <- 2L
    else if (state$nit >= maxiter) state$status <- 0L
    !is.null(state$status)
  }
  if (method == "equality_constrained_sqp") {
    stop_criteria <- function(state, x, last_iteration_failed, optimality, constr_violation, tr_radius,
                              constr_penalty, cg_info) {
      update_state(state, x, last_iteration_failed, tr_radius, constr_penalty, cg_info)
      if (verbose == 2) .sptc_print_iteration("basic", state)
      else if (verbose > 2) .sptc_print_iteration("sqp", state)
      finish_check(state, isTRUE(state$tr_radius < xtol))
    }
  } else {
    stop_criteria <- function(state, x, last_iteration_failed, tr_radius, constr_penalty, cg_info,
                              barrier_parameter, barrier_tolerance) {
      update_state(state, x, last_iteration_failed, tr_radius, constr_penalty, cg_info)
      state$barrier_parameter <- barrier_parameter
      state$barrier_tolerance <- barrier_tolerance
      if (verbose == 2) .sptc_print_iteration("basic", state)
      else if (verbose > 2) .sptc_print_iteration("ip", state)
      finish_check(state, isTRUE(state$tr_radius < xtol) && isTRUE(state$barrier_parameter < barrier_tol))
    }
  }
  kind <- if (verbose == 2) "basic" else if (method == "equality_constrained_sqp") "sqp" else "ip"
  if (verbose >= 2) .sptc_print_header(kind)
  if (method == "equality_constrained_sqp") {
    .sptc_eqsqp(function(x) list(f = objective$fun(x), b = canonical$fun(x)$eq, x = x),
                function(x) list(g = objective$grad(x), J = canonical$jac(x)$eq),
                function(x, v) lagrangian_hess(x, v), x0, objective$f, objective$g, ic$c_eq, ic$J_eq, stop_criteria,
                state, initial_constr_penalty, initial_tr_radius, factorization_method, sparse)
  } else {
    .sptc_tr_interior_point(objective$fun, objective$grad, lagrangian_hess, n_vars, canonical$n_ineq, canonical$n_eq,
                            canonical$fun, canonical$jac, x0, objective$f, objective$g, ic$c_ineq, ic$J_ineq, ic$c_eq,
                            ic$J_eq, stop_criteria, canonical$keep_feasible, xtol, state, initial_barrier_parameter,
                            initial_barrier_tolerance, initial_constr_penalty, initial_tr_radius,
                            factorization_method, fd_lb, fd_ub, sparse)
  }
  if (state$status %in% c(1L, 2L) && isTRUE(state$constr_violation > gtol)) state$status <- 4L
  state$success <- state$status %in% c(1L, 2L)
  state$message <- .sptc_messages[[as.character(state$status)]]
  state$niter <- state$nit
  if (verbose >= 2) cat("\n")
  if (verbose >= 1) {
    cat(state$message, "\n", sep = "")
    cat("Number of iterations: ", state$nit, ", function evaluations: ", state$nfev, ", CG iterations: ",
        state$cg_niter, ", optimality: ", .sptc_fmt_e(state$optimality, 2L), ", constraint violation: ",
        .sptc_fmt_e(state$constr_violation, 2L), ", execution time: ", .sptc_fmt_time(state$execution_time),
        " s.\n", sep = "")
  }
  res <- as.list(state)
  first <- c("x", "fun", "success", "status", "message", "nit", "nfev", "njev", "nhev")
  res[c(first, setdiff(names(res), first))]
}
