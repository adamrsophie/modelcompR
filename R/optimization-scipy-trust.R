# The trust region methods of scipy.optimize.minimize (SciPy 1.18.1,
# _trustregion.py, _trustregion_dogleg.py, _trustregion_ncg.py,
# _trustregion_exact.py and _trustregion_krylov.py with the GLTR solver of
# the trlib library, scipy/optimize/_trlib).

.sp_norm2 <- function(x) sqrt(sum(x * x))

# Cholesky factorization A = U'U column by column (LAPACK dpotf2, upper):
# list(U, info) with info the order of the first leading minor that is not
# positive definite (0 on success).
.sp_potrf <- function(A) {
  n <- nrow(A)
  U <- matrix(0, n, n)
  for (j in seq_len(n)) {
    ajj <- A[j, j] - sum(U[seq_len(j - 1L), j]^2)
    if (ajj <= 0 || is.nan(ajj)) {
      U[j, j] <- ajj
      return(list(U = U, info = j))
    }
    ajj <- sqrt(ajj)
    U[j, j] <- ajj
    if (j < n) {
      rest <- (j + 1L):n
      v <- A[j, rest] - as.numeric(crossprod(U[seq_len(j - 1L), j], U[seq_len(j - 1L), rest, drop = FALSE]))
      U[j, rest] <- v * (1 / ajj)
    }
  }
  list(U = U, info = 0L)
}

.sp_boundaries_intersections <- function(z, d, trust_radius) {
  a <- sum(d * d)
  b <- 2 * sum(z * d)
  c <- sum(z * z) - trust_radius^2
  sqrt_discriminant <- sqrt(b * b - 4 * a * c)
  aux <- b + abs(sqrt_discriminant) * (if (b < 0 || (b == 0 && 1 / b < 0)) -1 else 1)
  ta <- -aux / (2 * a)
  tb <- -2 * c / aux
  sort(c(ta, tb))
}

.sp_estimate_smallest_singular_value <- function(U) {
  n <- nrow(U)
  p <- numeric(n)
  w <- numeric(n)
  UT <- t(U)
  for (k in seq_len(n)) {
    wp <- (1 - p[k]) / UT[k, k]
    wm <- (-1 - p[k]) / UT[k, k]
    rest <- if (k < n) (k + 1L):n else integer()
    pp <- p[rest] + UT[rest, k] * wp
    pm <- p[rest] + UT[rest, k] * wm
    if (abs(wp) + sum(abs(pp)) >= abs(wm) + sum(abs(pm))) {
      w[k] <- wp
      p[rest] <- pp
    } else {
      w[k] <- wm
      p[rest] <- pm
    }
  }
  v <- backsolve(U, w)
  v_norm <- .sp_norm2(v)
  list(s_min = .sp_norm2(w) / v_norm, z_min = v / v_norm)
}

.sp_gershgorin_bounds <- function(H) {
  d <- diag(H)
  rs <- rowSums(abs(H))
  c(min(d + abs(d) - rs), max(d - abs(d) + rs))
}

.sp_singular_leading_submatrix <- function(A, U, k) {
  delta <- sum(U[seq_len(k - 1L), k]^2) - A[k, k]
  n <- nrow(A)
  v <- numeric(n)
  v[k] <- 1
  if (k != 1L) v[seq_len(k - 1L)] <- backsolve(U[seq_len(k - 1L), seq_len(k - 1L), drop = FALSE], -U[seq_len(k - 1L), k])
  list(delta = delta, v = v)
}

# A quadratic subproblem at x (BaseQuadraticSubproblem and its children):
# values are evaluated lazily and cached in the object.
.sp_subproblem <- function(type, x, fun, jac, hess, hessp, maxiter = NULL, trlib = NULL) {
  m <- new.env(parent = emptyenv())
  m$x <- x
  m$fun <- function() {
    if (is.null(m$f_)) m$f_ <- fun(x)
    m$f_
  }
  m$jac <- function() {
    if (is.null(m$g_)) m$g_ <- jac(x)
    m$g_
  }
  m$hess <- function() {
    if (is.null(m$h_)) m$h_ <- hess(x)
    m$h_
  }
  m$hessp <- function(p) if (!is.null(hessp)) hessp(x, p) else .sp_hdot(m$hess(), p)
  m$jac_mag <- function() {
    if (is.null(m$gm_)) m$gm_ <- .sp_norm2(m$jac())
    m$gm_
  }
  m$model <- function(p) m$fun() + sum(m$jac() * p) + 0.5 * sum(p * m$hessp(p))
  if (type == "dogleg") {
    m$solve <- function(trust_radius) {
      if (is.null(m$newton_)) {
        g <- m$jac()
        R <- tryCatch(chol(m$hess()), error = function(e) .sp_stop("sp_linalg", conditionMessage(e)))
        m$newton_ <- -backsolve(R, forwardsolve(t(R), g))
      }
      p_best <- m$newton_
      if (.sp_norm2(p_best) < trust_radius) return(list(p = p_best, hits = FALSE))
      if (is.null(m$cauchy_)) {
        g <- m$jac()
        Bg <- m$hessp(g)
        m$cauchy_ <- -(sum(g * g) / sum(g * Bg)) * g
      }
      p_u <- m$cauchy_
      p_u_norm <- .sp_norm2(p_u)
      if (p_u_norm >= trust_radius) return(list(p = p_u * (trust_radius / p_u_norm), hits = TRUE))
      tb <- .sp_boundaries_intersections(p_u, p_best - p_u, trust_radius)[2L]
      list(p = p_u + tb * (p_best - p_u), hits = TRUE)
    }
  } else if (type == "ncg") {
    m$solve <- function(trust_radius) {
      p_origin <- numeric(length(m$jac()))
      tolerance <- min(0.5, sqrt(m$jac_mag())) * m$jac_mag()
      if (m$jac_mag() < tolerance) return(list(p = p_origin, hits = FALSE))
      z <- p_origin
      r <- m$jac()
      d <- -r
      repeat {
        Bd <- m$hessp(d)
        dBd <- sum(d * Bd)
        if (dBd <= 0) {
          t <- .sp_boundaries_intersections(z, d, trust_radius)
          pa <- z + t[1L] * d
          pb <- z + t[2L] * d
          return(list(p = if (m$model(pa) < m$model(pb)) pa else pb, hits = TRUE))
        }
        r_squared <- sum(r * r)
        alpha <- r_squared / dBd
        z_next <- z + alpha * d
        if (.sp_norm2(z_next) >= trust_radius) {
          tb <- .sp_boundaries_intersections(z, d, trust_radius)[2L]
          return(list(p = z + tb * d, hits = TRUE))
        }
        r_next <- r + alpha * Bd
        r_next_squared <- sum(r_next * r_next)
        if (sqrt(r_next_squared) < tolerance) return(list(p = z_next, hits = FALSE))
        beta_next <- r_next_squared / r_squared
        d <- -r_next + beta_next * d
        z <- z_next
        r <- r_next
      }
    }
  } else if (type == "exact") {
    m$maxiter <- maxiter %||% 25
    H0 <- m$hess()
    m$dimension <- nrow(H0)
    gb <- .sp_gershgorin_bounds(H0)
    m$hess_inf <- max(rowSums(abs(H0)))
    m$hess_fro <- sqrt(sum(H0 * H0))
    m$close_to_zero <- m$dimension * .sp_eps * m$hess_inf
    m$previous_tr_radius <- -1
    m$lambda_lb <- NULL
    m$solve <- function(tr_radius) {
      H0 <- m$hess()
      n <- m$dimension
      lambda_ub <- max(0, m$jac_mag() / tr_radius + min(-gb[1L], m$hess_fro, m$hess_inf))
      lambda_lb <- max(0, -min(diag(H0)), m$jac_mag() / tr_radius - min(gb[2L], m$hess_fro, m$hess_inf))
      if (tr_radius < m$previous_tr_radius) lambda_lb <- max(m$lambda_lb, lambda_lb)
      lambda_current <- if (lambda_lb == 0) 0 else
        max(sqrt(lambda_lb * lambda_ub), lambda_lb + 0.01 * (lambda_ub - lambda_lb))
      hits <- TRUE
      already <- FALSE
      niter <- 0
      p <- NULL
      while (niter < m$maxiter) {
        if (already) {
          already <- FALSE
        } else {
          H <- H0 + lambda_current * diag(1, n)
          ch <- .sp_potrf(H)
          U <- ch$U; info <- ch$info
        }
        niter <- niter + 1
        if (info == 0 && m$jac_mag() > m$close_to_zero) {
          p <- -backsolve(U, forwardsolve(t(U), m$jac()))
          p_norm <- .sp_norm2(p)
          if (p_norm <= tr_radius && lambda_current == 0) {
            hits <- FALSE
            break
          }
          w <- forwardsolve(t(U), p)
          w_norm <- .sp_norm2(w)
          delta_lambda <- (p_norm / w_norm)^2 * (p_norm - tr_radius) / tr_radius
          lambda_new <- lambda_current + delta_lambda
          if (p_norm < tr_radius) {
            ss <- .sp_estimate_smallest_singular_value(U)
            t <- .sp_boundaries_intersections(p, ss$z_min, tr_radius)
            step_len <- t[which.min(abs(t))]
            quadratic_term <- sum(p * as.numeric(H %*% p))
            relative_error <- (step_len^2 * ss$s_min^2) / (quadratic_term + lambda_current * tr_radius^2)
            if (relative_error <= 0.2) {
              p <- p + step_len * ss$z_min
              break
            }
            lambda_ub <- lambda_current
            lambda_lb <- max(lambda_lb, lambda_current - ss$s_min^2)
            H <- H0 + lambda_new * diag(1, n)
            ch <- .sp_potrf(H)
            if (ch$info == 0) {
              lambda_current <- lambda_new
              already <- TRUE
              U <- ch$U
              info <- 0L
            } else {
              lambda_lb <- max(lambda_lb, lambda_new)
              lambda_current <- max(sqrt(abs(lambda_lb * lambda_ub)), lambda_lb + 0.01 * (lambda_ub - lambda_lb))
            }
          } else {
            relative_error <- abs(p_norm - tr_radius) / tr_radius
            if (relative_error <= 0.1) break
            lambda_lb <- lambda_current
            lambda_current <- lambda_new
          }
        } else if (info == 0 && m$jac_mag() <= m$close_to_zero) {
          if (lambda_current == 0) {
            p <- numeric(n)
            hits <- FALSE
            break
          }
          ss <- .sp_estimate_smallest_singular_value(U)
          step_len <- tr_radius
          p <- step_len * ss$z_min
          if (step_len^2 * ss$s_min^2 <= 0.2 * lambda_current * tr_radius^2) break
          lambda_ub <- lambda_current
          lambda_lb <- max(lambda_lb, lambda_current - ss$s_min^2)
          lambda_current <- max(sqrt(lambda_lb * lambda_ub), lambda_lb + 0.01 * (lambda_ub - lambda_lb))
        } else {
          sl <- .sp_singular_leading_submatrix(H, U, info)
          v_norm <- .sp_norm2(sl$v)
          lambda_lb <- max(lambda_lb, lambda_current + sl$delta / v_norm^2)
          lambda_current <- max(sqrt(abs(lambda_lb * lambda_ub)), lambda_lb + 0.01 * (lambda_ub - lambda_lb))
        }
      }
      m$lambda_lb <- lambda_lb
      m$previous_tr_radius <- tr_radius
      list(p = p, hits = hits)
    }
  } else if (type == "krylov") {
    # TRLIBQuadraticSubproblem sizes its work arrays from the gradient, so
    # every subproblem evaluates it on construction
    m$jac()
    m$solve <- function(trust_radius) .sp_trlib_solve(m, trust_radius, trlib)
  }
  m
}

.sp_trust_region <- function(fun, x0, jac = NULL, hess = NULL, hessp = NULL, subproblem = "ncg",
                             initial_trust_radius = 1.0, max_trust_radius = 1000.0, eta = 0.15, gtol = 1e-4,
                             maxiter = NULL, inexact = TRUE, subproblem_maxiter = NULL, ...) {
  if (is.null(jac)) stop("Jacobian is currently required for trust-region methods", call. = FALSE)
  if (subproblem %in% c("dogleg", "exact") && !is.function(hess))
    stop("Method ", if (subproblem == "dogleg") "dogleg" else "trust-exact", " needs the Hessian matrix; build the ",
         "problem with hess = TRUE or pass a Hessian function", call. = FALSE)
  if (is.null(hess) && is.null(hessp))
    stop("Either the Hessian or the Hessian-vector product is currently required for trust-region methods",
         call. = FALSE)
  if (!(0 <= eta && eta < 0.25)) stop("invalid acceptance stringency", call. = FALSE)
  if (max_trust_radius <= 0) stop("the max trust radius must be positive", call. = FALSE)
  if (initial_trust_radius <= 0) stop("the initial trust radius must be positive", call. = FALSE)
  if (initial_trust_radius >= max_trust_radius)
    stop("the initial trust radius must be less than the max trust radius", call. = FALSE)
  x0 <- as.numeric(x0)
  sf <- .sp_prepare_scalar_function(fun, x0, jac = jac, hess = hess)
  fn <- sf$fun
  jc <- sf$grad
  nh <- new.env(parent = emptyenv())
  nh$n <- 0L
  if (is.function(hess)) {
    hs <- sf$hess
    hp <- NULL
  } else {
    hs <- NULL
    hp <- function(x, p) {
      nh$n <- nh$n + 1L
      .sp_hdot(sf$hess(x), p)
    }
  }
  if (is.null(maxiter)) maxiter <- length(x0) * 200
  trlib <- if (subproblem == "krylov") (if (inexact) list(tol_rel_i = -2.0, tol_rel_b = -3.0) else
    list(tol_rel_i = 1e-8, tol_rel_b = 1e-6))
  mk <- function(x) .sp_subproblem(subproblem, x, fn, jc, hs, hp, subproblem_maxiter, trlib)
  warnflag <- 0L
  trust_radius <- initial_trust_radius
  x <- x0
  m <- mk(x)
  k <- 0L
  # comparisons with NaN are false, as in Python
  while (isTRUE(m$jac_mag() >= gtol)) {
    st <- tryCatch(m$solve(trust_radius), sp_linalg = function(e) NULL)
    if (is.null(st)) {
      warnflag <- 3L
      break
    }
    p <- st$p
    predicted_value <- m$model(p)
    x_proposed <- x + p
    m_proposed <- mk(x_proposed)
    actual_reduction <- m$fun() - m_proposed$fun()
    predicted_reduction <- m$fun() - predicted_value
    if (isTRUE(predicted_reduction <= 0)) {
      warnflag <- 2L
      break
    }
    rho <- actual_reduction / predicted_reduction
    if (isTRUE(rho < 0.25)) {
      trust_radius <- trust_radius * 0.25
    } else if (isTRUE(rho > 0.75) && st$hits) {
      trust_radius <- min(2 * trust_radius, max_trust_radius)
    }
    if (isTRUE(rho > eta)) {
      x <- x_proposed
      m <- m_proposed
    }
    k <- k + 1L
    if (isTRUE(m$jac_mag() < gtol)) {
      warnflag <- 0L
      break
    }
    if (k >= maxiter) {
      warnflag <- 1L
      break
    }
  }
  msgs <- c(.sp_status_message$success, .sp_status_message$maxiter,
            "A bad approximation caused failure to predict improvement.",
            "A linalg error occurred, such as a non-psd Hessian.")
  list(x = x, success = warnflag == 0L, status = warnflag, fun = m$fun(), jac = m$jac(), nfev = sf$nfev(),
       njev = sf$ngev, nhev = sf$nhev + nh$n, nit = k, message = msgs[warnflag + 1L])
}
