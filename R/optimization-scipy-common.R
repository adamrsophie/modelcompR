# Shared machinery of the SciPy ports (scipy/optimize of SciPy 1.18.1):
# finite differences (_numdiff.py approx_derivative), the memoizing
# ScalarFunction (_differentiable_functions.py), the quasi Newton Hessian
# updates (_hessian_update_strategy.py) and the Wolfe line searches
# (_linesearch.py, _dcsrch.py), and the minimize() dispatcher (_minimize.py).

.sp_eps <- .Machine$double.eps

.sp_eps_for_method <- function(method) {
  if (method %in% c("2-point", "cs")) sqrt(.sp_eps) else .sp_eps^(1 / 3)
}

.sp_stop <- function(cls, msg = cls, data = NULL)
  stop(structure(class = c(cls, "error", "condition"), list(message = msg, call = NULL, data = data)))

# _adjust_scheme_to_bounds
.sp_adjust_scheme <- function(x0, h, num_steps, scheme, lb, ub) {
  if (scheme == "1-sided") {
    one <- rep(TRUE, length(h))
  } else {
    h <- abs(h)
    one <- rep(FALSE, length(h))
  }
  if (all(lb == -Inf & ub == Inf)) return(list(h = h, one = one))
  h_total <- h * num_steps
  ha <- h
  lower <- x0 - lb
  upper <- ub - x0
  if (scheme == "1-sided") {
    x <- x0 + h_total
    violated <- x < lb | x > ub
    fitting <- abs(h_total) <= pmax(lower, upper)
    ha[violated & fitting] <- -ha[violated & fitting]
    fw <- upper >= lower & !fitting
    ha[fw] <- upper[fw] / num_steps
    bw <- upper < lower & !fitting
    ha[bw] <- -lower[bw] / num_steps
  } else {
    central <- lower >= h_total & upper >= h_total
    fw <- upper >= lower & !central
    ha[fw] <- pmin(h[fw], 0.5 * upper[fw] / num_steps)
    one[fw] <- TRUE
    bw <- upper < lower & !central
    ha[bw] <- -pmin(h[bw], 0.5 * lower[bw] / num_steps)
    one[bw] <- TRUE
    min_dist <- pmin(upper, lower) / num_steps
    adj <- !central & abs(ha) <= min_dist
    ha[adj] <- min_dist[adj]
    one[adj] <- FALSE
  }
  list(h = ha, one = one)
}

# approx_derivative for dense Jacobians: returns list(J, nfev); J is a
# vector for a scalar function and an m x n matrix otherwise.
.sp_approx_derivative <- function(fun, x0, method = "3-point", rel_step = NULL, abs_step = NULL, f0 = NULL,
                                  lb = -Inf, ub = Inf) {
  n <- length(x0)
  lb <- rep_len(lb, n); ub <- rep_len(ub, n)
  nfev <- 0L
  if (is.null(f0)) {
    f0 <- fun(x0)
    nfev <- 1L
  }
  if (any(x0 < lb | x0 > ub)) stop("`x0` violates bound constraints.", call. = FALSE)
  sgn <- ifelse(x0 >= 0, 1, -1)
  rstep <- .sp_eps_for_method(method)
  if (is.null(abs_step)) {
    default <- rstep * sgn * pmax(1, abs(x0))
    if (is.null(rel_step)) {
      h <- default
    } else {
      h <- rel_step * sgn * abs(x0)
      dx <- (x0 + h) - x0
      h <- ifelse(dx == 0, default, h)
    }
  } else {
    h <- rep_len(abs_step, n)
    dx <- (x0 + h) - x0
    h <- ifelse(dx == 0, rstep * sgn * pmax(1, abs(x0)), h)
  }
  if (method == "2-point") {
    adj <- .sp_adjust_scheme(x0, h, 1, "1-sided", lb, ub)
  } else if (method == "3-point") {
    adj <- .sp_adjust_scheme(x0, h, 1, "2-sided", lb, ub)
  } else {
    stop("finite difference method '", method, "' is not available", call. = FALSE)
  }
  h <- adj$h
  m <- length(f0)
  JT <- matrix(0, n, m)
  for (i in seq_len(n)) {
    if (method == "2-point") {
      x1 <- x0
      x1[i] <- x0[i] + h[i]
      f1 <- fun(x1)
      JT[i, ] <- (f1 - f0) / ((x0[i] + h[i]) - x0[i])
      nfev <- nfev + 1L
    } else {
      x1 <- x0; x2 <- x0
      if (adj$one[i]) {
        x1[i] <- x0[i] + h[i]
        x2[i] <- x0[i] + 2 * h[i]
        f1 <- fun(x1); f2 <- fun(x2)
        JT[i, ] <- (-3.0 * f0 + 4 * f1 - f2) / (x2[i] - x0[i])
      } else {
        x1[i] <- x0[i] - h[i]
        x2[i] <- x0[i] + h[i]
        f1 <- fun(x1); f2 <- fun(x2)
        JT[i, ] <- (f2 - f1) / (x2[i] - x1[i])
      }
      nfev <- nfev + 2L
    }
  }
  J <- if (m == 1L) as.numeric(JT) else t(JT)
  list(J = J, nfev = nfev)
}

.sp_same_x <- function(x, y) length(x) == length(y) && isTRUE(all(x == y))

# ScalarFunction: memoized objective, gradient and Hessian with SciPy's
# evaluation counts. `grad` is a function or "2-point"/"3-point"; `hess` a
# function, a Hessian update strategy (.sp_hessian_update) or NULL.
.sp_scalar_function <- function(fun, x0, grad, hess = NULL, rel_step = NULL, lb = -Inf, ub = Inf,
                                epsilon = NULL) {
  sf <- new.env(parent = emptyenv())
  sf$x <- as.numeric(x0)
  sf$n <- length(x0)
  sf$nfev0 <- 0L; sf$gnfev <- 0L; sf$ngev <- 0L; sf$nhev <- 0L
  sf$f_updated <- FALSE; sf$g_updated <- FALSE; sf$H_updated <- FALSE
  sf$lowest_x <- NULL; sf$lowest_f <- Inf
  fd <- is.character(grad)
  quasi <- inherits(hess, "SpHessianUpdate")
  wrapped_fun <- function(x) {
    fx <- fun(x)
    if (length(fx) != 1L) stop("The user-provided objective function must return a scalar value.", call. = FALSE)
    as.numeric(fx)
  }
  update_fun <- function() {
    if (!sf$f_updated) {
      fx <- wrapped_fun(sf$x)
      sf$nfev0 <- sf$nfev0 + 1L
      if (isTRUE(fx < sf$lowest_f)) {
        sf$lowest_x <- sf$x
        sf$lowest_f <- fx
      }
      sf$f <- fx
      sf$f_updated <- TRUE
    }
  }
  update_grad <- function() {
    if (!sf$g_updated) {
      if (fd) {
        update_fun()
        r <- .sp_approx_derivative(wrapped_fun, sf$x, method = grad, rel_step = if (is.null(epsilon)) rel_step,
                                   abs_step = epsilon, f0 = sf$f, lb = lb, ub = ub)
        sf$gnfev <- sf$gnfev + r$nfev
        sf$g <- r$J
      } else {
        sf$g <- as.numeric(grad(sf$x))
      }
      sf$ngev <- sf$ngev + 1L
      sf$g_updated <- TRUE
    }
  }
  update_hess <- function() {
    if (!sf$H_updated) {
      if (quasi) {
        update_grad()
        hess$update(sf$x - sf$x_prev, sf$g - sf$g_prev)
      } else if (is.function(hess)) {
        sf$nhev <- sf$nhev + 1L
        sf$H <- as.matrix(hess(sf$x))
      }
      sf$H_updated <- TRUE
    }
  }
  update_x <- function(x) {
    if (quasi) {
      update_grad()
      sf$x_prev <- sf$x
      sf$g_prev <- sf$g
      sf$x <- as.numeric(x)
      sf$f_updated <- FALSE; sf$g_updated <- FALSE; sf$H_updated <- FALSE
      update_hess()
    } else {
      sf$x <- as.numeric(x)
      sf$f_updated <- FALSE; sf$g_updated <- FALSE; sf$H_updated <- FALSE
    }
  }
  update_fun()
  update_grad()
  if (quasi) {
    hess$initialize(sf$n, "hess")
    sf$H <- hess
    sf$H_updated <- TRUE
    sf$x_prev <- NULL; sf$g_prev <- NULL
  } else if (is.function(hess)) {
    sf$H <- as.matrix(hess(sf$x))
    sf$nhev <- 1L
    sf$H_updated <- TRUE
  }
  sf$fun <- function(x) {
    if (!.sp_same_x(x, sf$x)) update_x(x)
    update_fun()
    sf$f
  }
  sf$grad <- function(x) {
    if (!.sp_same_x(x, sf$x)) update_x(x)
    update_grad()
    sf$g
  }
  sf$hess <- function(x) {
    if (!.sp_same_x(x, sf$x)) update_x(x)
    update_hess()
    sf$H
  }
  sf$fun_and_grad <- function(x) {
    if (!.sp_same_x(x, sf$x)) update_x(x)
    update_fun()
    update_grad()
    list(f = sf$f, g = sf$g)
  }
  sf$nfev <- function() sf$nfev0 + sf$gnfev
  sf
}

# _prepare_scalar_function of _optimize.py
.sp_prepare_scalar_function <- function(fun, x0, jac = NULL, lb = -Inf, ub = Inf, epsilon = NULL,
                                        finite_diff_rel_step = NULL, hess = NULL) {
  if (is.function(jac)) {
    grad <- jac
  } else if (is.character(jac)) {
    epsilon <- NULL
    grad <- jac
  } else {
    grad <- "2-point"
  }
  .sp_scalar_function(fun, x0, grad, hess, finite_diff_rel_step, lb, ub, epsilon = epsilon)
}

# ---------------------------------------------------------------------------
# Quasi Newton Hessian updates (BFGS and SR1 of _hessian_update_strategy.py)

.sp_hessian_update <- function(kind = c("BFGS", "SR1"), exception_strategy = "skip_update", min_curvature = NULL,
                               init_scale = "auto", min_denominator = 1e-8) {
  kind <- match.arg(kind)
  if (kind == "BFGS") {
    if (is.null(min_curvature)) min_curvature <- if (exception_strategy == "skip_update") 1e-8 else 0.2
  }
  st <- new.env(parent = emptyenv())
  st$first_iteration <- NULL
  st$initialize <- function(n, approx_type) {
    st$first_iteration <- TRUE
    st$n <- n
    st$approx_type <- approx_type
    st$B <- diag(1, n)
    st$H <- diag(1, n)
  }
  auto_scale <- function(delta_x, delta_grad) {
    s_norm2 <- sum(delta_x * delta_x)
    y_norm2 <- sum(delta_grad * delta_grad)
    ys <- abs(sum(delta_grad * delta_x))
    if (ys == 0 || y_norm2 == 0 || s_norm2 == 0) return(1)
    if (st$approx_type == "hess") y_norm2 / ys else ys / y_norm2
  }
  # BFGS: _update_inverse_hessian and _update_hessian
  bfgs_inv <- function(ys, Hy, yHy, s) {
    st$H <- st$H + (outer(s, s) * ((ys + yHy) / ys^2)) - (outer(Hy, s) + outer(s, Hy)) / ys
  }
  bfgs_hess <- function(ys, Bs, sBs, y) {
    st$B <- st$B + outer(y, y) / ys - outer(Bs, Bs) / sBs
  }
  update_bfgs <- function(delta_x, delta_grad) {
    if (st$approx_type == "hess") { w <- delta_x; z <- delta_grad } else { w <- delta_grad; z <- delta_x }
    wz <- sum(w * z)
    Mw <- as.numeric(st$dot(w))
    wMw <- sum(Mw * w)
    if (wMw <= 0) {
      scale <- auto_scale(delta_x, delta_grad)
      if (st$approx_type == "hess") st$B <- scale * diag(1, st$n) else st$H <- scale * diag(1, st$n)
      Mw <- as.numeric(st$dot(w))
      wMw <- sum(Mw * w)
    }
    if (wz <= min_curvature * wMw) {
      if (exception_strategy == "skip_update") return(invisible())
      update_factor <- (1 - min_curvature) / (1 - wz / wMw)
      z <- update_factor * z + (1 - update_factor) * Mw
      wz <- sum(w * z)
    }
    if (st$approx_type == "hess") bfgs_hess(wz, Mw, wMw, z) else bfgs_inv(wz, Mw, wMw, z)
  }
  update_sr1 <- function(delta_x, delta_grad) {
    if (st$approx_type == "hess") { w <- delta_x; z <- delta_grad } else { w <- delta_grad; z <- delta_x }
    Mw <- as.numeric(st$dot(w))
    z_minus_Mw <- z - Mw
    denominator <- sum(w * z_minus_Mw)
    if (abs(denominator) <= min_denominator * sqrt(sum(w * w)) * sqrt(sum(z_minus_Mw * z_minus_Mw))) return(invisible())
    if (st$approx_type == "hess") st$B <- st$B + outer(z_minus_Mw, z_minus_Mw) / denominator
    else st$H <- st$H + outer(z_minus_Mw, z_minus_Mw) / denominator
  }
  st$update <- function(delta_x, delta_grad) {
    if (all(delta_x == 0)) return(invisible())
    if (all(delta_grad == 0)) return(invisible())
    if (isTRUE(st$first_iteration)) {
      scale <- if (identical(init_scale, "auto")) auto_scale(delta_x, delta_grad) else init_scale
      if (is.matrix(scale)) {
        if (st$approx_type == "hess") st$B <- scale else st$H <- solve(scale)
      } else {
        if (st$approx_type == "hess") st$B <- st$B * scale else st$H <- st$H * scale
      }
      st$first_iteration <- FALSE
    }
    if (kind == "BFGS") update_bfgs(delta_x, delta_grad) else update_sr1(delta_x, delta_grad)
  }
  st$dot <- function(p) if (st$approx_type == "hess") st$B %*% p else st$H %*% p
  st$get_matrix <- function() if (st$approx_type == "hess") st$B else st$H
  st$kind <- kind
  class(st) <- c("SpHessianUpdate", "environment")
  st
}

# ---------------------------------------------------------------------------
# Line searches

# dcstep of _dcsrch.py (MINPACK-2)
.sp_dcstep <- function(stx, fx, dx, sty, fy, dy, stp, fp, dp, brackt, stpmin, stpmax) {
  sgnd <- sign(dp) * sign(dx)
  if (fp > fx) {
    theta <- 3.0 * (fx - fp) / (stp - stx) + dx + dp
    s <- max(abs(theta), abs(dx), abs(dp))
    gamma <- s * sqrt((theta / s)^2 - (dx / s) * (dp / s))
    if (stp < stx) gamma <- -gamma
    p <- (gamma - dx) + theta
    q <- ((gamma - dx) + gamma) + dp
    r <- p / q
    stpc <- stx + r * (stp - stx)
    stpq <- stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.0) * (stp - stx)
    stpf <- if (abs(stpc - stx) <= abs(stpq - stx)) stpc else stpc + (stpq - stpc) / 2.0
    brackt <- TRUE
  } else if (sgnd < 0.0) {
    theta <- 3 * (fx - fp) / (stp - stx) + dx + dp
    s <- max(abs(theta), abs(dx), abs(dp))
    gamma <- s * sqrt((theta / s)^2 - (dx / s) * (dp / s))
    if (stp > stx) gamma <- -gamma
    p <- (gamma - dp) + theta
    q <- ((gamma - dp) + gamma) + dx
    r <- p / q
    stpc <- stp + r * (stx - stp)
    stpq <- stp + (dp / (dp - dx)) * (stx - stp)
    stpf <- if (abs(stpc - stp) > abs(stpq - stp)) stpc else stpq
    brackt <- TRUE
  } else if (abs(dp) < abs(dx)) {
    theta <- 3 * (fx - fp) / (stp - stx) + dx + dp
    s <- max(abs(theta), abs(dx), abs(dp))
    gamma <- s * sqrt(max(0, (theta / s)^2 - (dx / s) * (dp / s)))
    if (stp > stx) gamma <- -gamma
    p <- (gamma - dp) + theta
    q <- (gamma + (dx - dp)) + gamma
    r <- p / q
    stpc <- if (r < 0 && gamma != 0) stp + r * (stx - stp) else if (stp > stx) stpmax else stpmin
    stpq <- stp + (dp / (dp - dx)) * (stx - stp)
    if (brackt) {
      stpf <- if (abs(stpc - stp) < abs(stpq - stp)) stpc else stpq
      stpf <- if (stp > stx) min(stp + 0.66 * (sty - stp), stpf) else max(stp + 0.66 * (sty - stp), stpf)
    } else {
      stpf <- if (abs(stpc - stp) > abs(stpq - stp)) stpc else stpq
      stpf <- min(max(stpf, stpmin), stpmax)
    }
  } else {
    if (brackt) {
      theta <- 3.0 * (fp - fy) / (sty - stp) + dy + dp
      s <- max(abs(theta), abs(dy), abs(dp))
      gamma <- s * sqrt((theta / s)^2 - (dy / s) * (dp / s))
      if (stp > sty) gamma <- -gamma
      p <- (gamma - dp) + theta
      q <- ((gamma - dp) + gamma) + dy
      r <- p / q
      stpf <- stp + r * (sty - stp)
    } else if (stp > stx) {
      stpf <- stpmax
    } else {
      stpf <- stpmin
    }
  }
  if (fp > fx) {
    sty <- stp; fy <- fp; dy <- dp
  } else {
    if (sgnd < 0) {
      sty <- stx; fy <- fx; dy <- dx
    }
    stx <- stp; fx <- fp; dx <- dp
  }
  list(stx = stx, fx = fx, dx = dx, sty = sty, fy = fy, dy = dy, stp = stpf, brackt = brackt)
}

# DCSRCH of _dcsrch.py: the More and Thuente line search. Returns
# list(stp, phi1, phi0, task); stp is NULL on failure.
# the reverse communication iteration of MINPACK-2 dcsrch: returns a function
# iterate(stp, f, g, task) that keeps the search state between calls
.sp_dcsrch_state <- function(ftol, gtol, xtol, stpmin, stpmax) {
  s <- new.env(parent = emptyenv())
  iterate <- function(stp, f, g, task) {
    p5 <- 0.5; p66 <- 0.66; xtrapl <- 1.1; xtrapu <- 4.0
    if (startsWith(task, "START")) {
      if (stp < stpmin) task <- "ERROR: STP .LT. STPMIN"
      if (stp > stpmax) task <- "ERROR: STP .GT. STPMAX"
      if (g >= 0) task <- "ERROR: INITIAL G .GE. ZERO"
      if (ftol < 0) task <- "ERROR: FTOL .LT. ZERO"
      if (gtol < 0) task <- "ERROR: GTOL .LT. ZERO"
      if (xtol < 0) task <- "ERROR: XTOL .LT. ZERO"
      if (stpmin < 0) task <- "ERROR: STPMIN .LT. ZERO"
      if (stpmax < stpmin) task <- "ERROR: STPMAX .LT. STPMIN"
      if (startsWith(task, "ERROR")) return(list(stp = stp, f = f, g = g, task = task))
      s$brackt <- FALSE
      s$stage <- 1L
      s$finit <- f; s$ginit <- g
      s$gtest <- ftol * s$ginit
      s$width <- stpmax - stpmin
      s$width1 <- s$width / p5
      s$stx <- 0.0; s$fx <- s$finit; s$gx <- s$ginit
      s$sty <- 0.0; s$fy <- s$finit; s$gy <- s$ginit
      s$stmin <- 0
      s$stmax <- stp + xtrapu * stp
      return(list(stp = stp, f = f, g = g, task = "FG"))
    }
    ftest <- s$finit + stp * s$gtest
    if (s$stage == 1L && f <= ftest && g >= 0) s$stage <- 2L
    if (s$brackt && (stp <= s$stmin || stp >= s$stmax)) task <- "WARNING: ROUNDING ERRORS PREVENT PROGRESS"
    if (s$brackt && s$stmax - s$stmin <= xtol * s$stmax) task <- "WARNING: XTOL TEST SATISFIED"
    if (stp == stpmax && f <= ftest && g <= s$gtest) task <- "WARNING: STP = STPMAX"
    if (stp == stpmin && (f > ftest || g >= s$gtest)) task <- "WARNING: STP = STPMIN"
    if (f <= ftest && abs(g) <= gtol * -s$ginit) task <- "CONVERGENCE"
    if (startsWith(task, "WARN") || startsWith(task, "CONV")) return(list(stp = stp, f = f, g = g, task = task))
    if (s$stage == 1L && f <= s$fx && f > ftest) {
      fm <- f - stp * s$gtest
      fxm <- s$fx - s$stx * s$gtest
      fym <- s$fy - s$sty * s$gtest
      gm <- g - s$gtest
      gxm <- s$gx - s$gtest
      gym <- s$gy - s$gtest
      r <- .sp_dcstep(s$stx, fxm, gxm, s$sty, fym, gym, stp, fm, gm, s$brackt, s$stmin, s$stmax)
      s$stx <- r$stx; s$sty <- r$sty; stp <- r$stp; s$brackt <- r$brackt
      s$fx <- r$fx + s$stx * s$gtest
      s$fy <- r$fy + s$sty * s$gtest
      s$gx <- r$dx + s$gtest
      s$gy <- r$dy + s$gtest
    } else {
      r <- .sp_dcstep(s$stx, s$fx, s$gx, s$sty, s$fy, s$gy, stp, f, g, s$brackt, s$stmin, s$stmax)
      s$stx <- r$stx; s$fx <- r$fx; s$gx <- r$dx
      s$sty <- r$sty; s$fy <- r$fy; s$gy <- r$dy
      stp <- r$stp; s$brackt <- r$brackt
    }
    if (s$brackt) {
      if (abs(s$sty - s$stx) >= p66 * s$width1) stp <- s$stx + p5 * (s$sty - s$stx)
      s$width1 <- s$width
      s$width <- abs(s$sty - s$stx)
    }
    if (s$brackt) {
      s$stmin <- min(s$stx, s$sty)
      s$stmax <- max(s$stx, s$sty)
    } else {
      s$stmin <- stp + xtrapl * (stp - s$stx)
      s$stmax <- stp + xtrapu * (stp - s$stx)
    }
    stp <- min(max(stp, stpmin), stpmax)
    if (s$brackt && (stp <= s$stmin || stp >= s$stmax) || (s$brackt && s$stmax - s$stmin <= xtol * s$stmax))
      stp <- s$stx
    list(stp = stp, f = f, g = g, task = "FG")
  }
  iterate
}

.sp_dcsrch <- function(phi, derphi, ftol, gtol, xtol, stpmin, stpmax, alpha1, phi0, derphi0, maxiter = 100L) {
  iterate <- .sp_dcsrch_state(ftol, gtol, xtol, stpmin, stpmax)
  phi1 <- phi0
  derphi1 <- derphi0
  task <- "START"
  stp <- alpha1
  done <- FALSE
  for (i in seq_len(maxiter)) {
    r <- iterate(alpha1, phi1, derphi1, task)
    stp <- r$stp; phi1 <- r$f; derphi1 <- r$g; task <- r$task
    if (!is.finite(stp)) {
      task <- "WARN"
      stp <- NULL
      done <- TRUE
      break
    }
    if (startsWith(task, "FG")) {
      alpha1 <- stp
      phi1 <- phi(stp)
      derphi1 <- derphi(stp)
    } else {
      done <- TRUE
      break
    }
  }
  if (!done) {
    stp <- NULL
    task <- "WARNING: dcsrch did not converge within max iterations"
  }
  if (startsWith(task, "ERROR") || startsWith(task, "WARN")) stp <- NULL
  list(stp = stp, phi1 = phi1, phi0 = phi0, task = task)
}

# line_search_wolfe1
.sp_line_search_wolfe1 <- function(f, fprime, xk, pk, gfk, old_fval, old_old_fval, c1 = 1e-4, c2 = 0.9,
                                   amax = 50, amin = 1e-8, xtol = 1e-14) {
  st <- new.env(parent = emptyenv())
  st$gval <- gfk; st$fc <- 0L; st$gc <- 0L
  phi <- function(s) {
    st$fc <- st$fc + 1L
    f(xk + s * pk)
  }
  derphi <- function(s) {
    st$gval <- fprime(xk + s * pk)
    st$gc <- st$gc + 1L
    sum(st$gval * pk)
  }
  derphi0 <- sum(gfk * pk)
  phi0 <- old_fval
  if (!is.null(old_old_fval) && derphi0 != 0) {
    alpha1 <- min(1.0, 1.01 * 2 * (phi0 - old_old_fval) / derphi0)
    if (alpha1 < 0) alpha1 <- 1.0
  } else {
    alpha1 <- 1.0
  }
  r <- .sp_dcsrch(phi, derphi, c1, c2, xtol, amin, amax, alpha1, phi0, derphi0, 100L)
  list(alpha = r$stp, fc = st$fc, gc = st$gc, fval = r$phi1, old_fval = r$phi0, gval = st$gval)
}

.sp_cubicmin <- function(a, fa, fpa, b, fb, c, fc) {
  C <- fpa
  db <- b - a
  dc <- c - a
  denom <- (db * dc)^2 * (db - dc)
  v1 <- fb - fa - C * db
  v2 <- fc - fa - C * dc
  A <- dc^2 * v1 + (-db^2) * v2
  B <- (-dc^3) * v1 + db^3 * v2
  if (!is.finite(denom) || denom == 0) return(NULL)
  A <- A / denom
  B <- B / denom
  radical <- B * B - 3 * A * C
  if (!is.finite(radical) || radical < 0 || A == 0) return(NULL)
  xmin <- a + (-B + sqrt(radical)) / (3 * A)
  if (!is.finite(xmin)) return(NULL)
  xmin
}

.sp_quadmin <- function(a, fa, fpa, b, fb) {
  D <- fa
  C <- fpa
  db <- b - a * 1.0
  if (db == 0) return(NULL)
  B <- (fb - D - C * db) / (db * db)
  if (!is.finite(B) || B == 0) return(NULL)
  xmin <- a - C / (2.0 * B)
  if (!is.finite(xmin)) return(NULL)
  xmin
}

.sp_zoom <- function(a_lo, a_hi, phi_lo, phi_hi, derphi_lo, phi, derphi, phi0, derphi0, c1, c2, extra_condition) {
  maxiter <- 10L
  i <- 0L
  delta1 <- 0.2; delta2 <- 0.1
  phi_rec <- phi0
  a_rec <- 0
  repeat {
    dalpha <- a_hi - a_lo
    if (dalpha < 0) { a <- a_hi; b <- a_lo } else { a <- a_lo; b <- a_hi }
    a_j <- NULL
    if (i > 0L) {
      cchk <- delta1 * dalpha
      a_j <- .sp_cubicmin(a_lo, phi_lo, derphi_lo, a_hi, phi_hi, a_rec, phi_rec)
    }
    if (i == 0L || is.null(a_j) || a_j > b - cchk || a_j < a + cchk) {
      qchk <- delta2 * dalpha
      a_j <- .sp_quadmin(a_lo, phi_lo, derphi_lo, a_hi, phi_hi)
      if (is.null(a_j) || a_j > b - qchk || a_j < a + qchk) a_j <- a_lo + 0.5 * dalpha
    }
    phi_aj <- phi(a_j)
    if (phi_aj > phi0 + c1 * a_j * derphi0 || phi_aj >= phi_lo) {
      phi_rec <- phi_hi
      a_rec <- a_hi
      a_hi <- a_j
      phi_hi <- phi_aj
    } else {
      derphi_aj <- derphi(a_j)
      if (abs(derphi_aj) <= -c2 * derphi0 && extra_condition(a_j, phi_aj))
        return(list(a = a_j, phi = phi_aj, derphi = derphi_aj))
      if (derphi_aj * (a_hi - a_lo) >= 0) {
        phi_rec <- phi_hi
        a_rec <- a_hi
        a_hi <- a_lo
        phi_hi <- phi_lo
      } else {
        phi_rec <- phi_lo
        a_rec <- a_lo
      }
      a_lo <- a_j
      phi_lo <- phi_aj
      derphi_lo <- derphi_aj
    }
    i <- i + 1L
    if (i > maxiter) return(list(a = NULL, phi = NULL, derphi = NULL))
  }
}

# line_search_wolfe2
.sp_line_search_wolfe2 <- function(f, fprime, xk, pk, gfk, old_fval, old_old_fval, c1 = 1e-4, c2 = 0.9,
                                   amax = NULL, extra_condition = NULL, maxiter = 10L) {
  st <- new.env(parent = emptyenv())
  st$fc <- 0L; st$gc <- 0L; st$gval <- NULL; st$gval_alpha <- NULL
  phi <- function(alpha) {
    st$fc <- st$fc + 1L
    f(xk + alpha * pk)
  }
  derphi <- function(alpha) {
    st$gc <- st$gc + 1L
    st$gval <- fprime(xk + alpha * pk)
    st$gval_alpha <- alpha
    sum(st$gval * pk)
  }
  derphi0 <- sum(gfk * pk)
  extra2 <- if (!is.null(extra_condition)) function(alpha, ph) {
    if (!identical(st$gval_alpha, alpha)) derphi(alpha)
    extra_condition(alpha, xk + alpha * pk, ph, st$gval)
  } else function(alpha, ph) TRUE
  phi0 <- old_fval
  old_phi0 <- old_old_fval
  alpha0 <- 0
  alpha1 <- if (!is.null(old_phi0) && derphi0 != 0) min(1.0, 1.01 * 2 * (phi0 - old_phi0) / derphi0) else 1.0
  if (alpha1 < 0) alpha1 <- 1.0
  if (!is.null(amax)) alpha1 <- min(alpha1, amax)
  phi_a1 <- phi(alpha1)
  phi_a0 <- phi0
  derphi_a0 <- derphi0
  res <- NULL
  for (i in seq_len(maxiter) - 1L) {
    if (alpha1 == 0 || (!is.null(amax) && alpha0 > amax)) {
      res <- list(a = NULL, phi = phi0, derphi = NULL)
      phi0 <- old_phi0
      break
    }
    if (phi_a1 > phi0 + c1 * alpha1 * derphi0 || (phi_a1 >= phi_a0 && i > 0L)) {
      res <- .sp_zoom(alpha0, alpha1, phi_a0, phi_a1, derphi_a0, phi, derphi, phi0, derphi0, c1, c2, extra2)
      break
    }
    derphi_a1 <- derphi(alpha1)
    if (abs(derphi_a1) <= -c2 * derphi0 && extra2(alpha1, phi_a1)) {
      res <- list(a = alpha1, phi = phi_a1, derphi = derphi_a1)
      break
    }
    if (derphi_a1 >= 0) {
      res <- .sp_zoom(alpha1, alpha0, phi_a1, phi_a0, derphi_a1, phi, derphi, phi0, derphi0, c1, c2, extra2)
      break
    }
    alpha2 <- 2 * alpha1
    if (!is.null(amax)) alpha2 <- min(alpha2, amax)
    alpha0 <- alpha1
    alpha1 <- alpha2
    phi_a0 <- phi_a1
    phi_a1 <- phi(alpha1)
    derphi_a0 <- derphi_a1
  }
  if (is.null(res)) res <- list(a = alpha1, phi = phi_a1, derphi = NULL)
  list(alpha = res$a, fc = st$fc, gc = st$gc, fval = res$phi, old_fval = phi0,
       gval = if (is.null(res$derphi)) NULL else st$gval)
}

# _line_search_wolfe12: wolfe1 then wolfe2; NULL alpha signals failure
.sp_line_search_wolfe12 <- function(f, fprime, xk, pk, gfk, old_fval, old_old_fval, c1 = 1e-4, c2 = 0.9,
                                    amax = 0x1.249ad2594c37dp+332, amin = 0x1.bff2ee48e0530p-333, extra_condition = NULL) {
  ret <- .sp_line_search_wolfe1(f, fprime, xk, pk, gfk, old_fval, old_old_fval, c1 = c1, c2 = c2, amax = amax,
                                amin = amin)
  if (!is.null(ret$alpha) && !is.null(extra_condition)) {
    xp1 <- xk + ret$alpha * pk
    if (!extra_condition(ret$alpha, xp1, ret$fval, ret$gval)) ret$alpha <- NULL
  }
  if (is.null(ret$alpha)) {
    ret <- .sp_line_search_wolfe2(f, fprime, xk, pk, gfk, old_fval, old_old_fval, c1 = c1, c2 = c2, amax = amax,
                                  extra_condition = extra_condition)
  }
  ret
}

# ---------------------------------------------------------------------------
# minimize()

.sp_status_message <- list(success = "Optimization terminated successfully.",
                           maxfev = "Maximum number of function evaluations has been exceeded.",
                           maxiter = "Maximum number of iterations has been exceeded.",
                           pr_loss = "Desired error not necessarily achieved due to precision loss.",
                           nan = "NaN result encountered.",
                           out_of_bounds = "The result is outside of the provided bounds.")

.sp_minimize <- function(fun, x0, method = "BFGS", jac = NULL, hess = NULL, bounds = NULL, constraints = list(),
                         options = list()) {
  meth <- tolower(method)
  fs <- function(x) as.numeric(fun(x))
  if (!is.null(bounds)) {
    lb <- rep_len(as.numeric(bounds$lb), length(x0)); ub <- rep_len(as.numeric(bounds$ub), length(x0))
    if (any(lb > ub)) stop("An upper bound is less than the corresponding lower bound.", call. = FALSE)
    bounds <- list(lb = lb, ub = ub)
    if (meth %in% c("tnc", "slsqp", "l-bfgs-b") && all(lb == ub)) {
      x <- lb
      return(list(x = x, fun = fs(x), nit = 0L, nfev = 1L, njev = 0L, status = 0L, success = TRUE,
                  message = "All independent variables were fixed by bounds."))
    }
  }
  args <- c(list(fs, x0), options)
  kernel <- c("nelder-mead" = ".sp_neldermead", powell = ".sp_powell", cg = ".sp_cg", bfgs = ".sp_bfgs",
              "newton-cg" = ".sp_newtoncg", "l-bfgs-b" = ".sp_lbfgsb", tnc = ".sp_tnc", cobyla = ".sp_cobyla",
              cobyqa = ".sp_cobyqa", slsqp = ".sp_slsqp", "trust-constr" = ".sp_trust_constr",
              "trust-krylov" = ".sp_trlib_solve")[meth]
  if (!is.na(kernel)) .scipy_need(kernel, paste0("minimize(method = '", method, "')"))
  res <- switch(meth,
    "nelder-mead" = do.call(.sp_neldermead, c(args, list(bounds = bounds))),
    "powell" = do.call(.sp_powell, c(args, list(bounds = bounds))),
    "cg" = do.call(.sp_cg, c(args, list(jac = jac))),
    "bfgs" = do.call(.sp_bfgs, c(args, list(jac = jac))),
    "newton-cg" = do.call(.sp_newtoncg, c(args, list(jac = jac, hess = hess))),
    "l-bfgs-b" = do.call(.sp_lbfgsb, c(args, list(jac = jac, bounds = bounds))),
    "tnc" = do.call(.sp_tnc, c(args, list(jac = jac, bounds = bounds))),
    "cobyla" = do.call(.sp_cobyla, c(args, list(constraints = constraints, bounds = bounds))),
    "cobyqa" = do.call(.sp_cobyqa, c(args, list(bounds = bounds, constraints = constraints))),
    "slsqp" = do.call(.sp_slsqp, c(args, list(jac = jac, bounds = bounds, constraints = constraints))),
    "trust-constr" = do.call(.sp_trust_constr, c(args, list(jac = jac, hess = hess, bounds = bounds,
                                                           constraints = constraints))),
    "dogleg" = do.call(.sp_trust_region, c(args, list(jac = jac, hess = hess, subproblem = "dogleg"))),
    "trust-ncg" = do.call(.sp_trust_region, c(args, list(jac = jac, hess = hess, subproblem = "ncg"))),
    "trust-krylov" = do.call(.sp_trust_region, c(args, list(jac = jac, hess = hess, subproblem = "krylov"))),
    "trust-exact" = do.call(.sp_trust_region, c(args, list(jac = jac, hess = hess, subproblem = "exact"))),
    stop("Unknown solver ", method, call. = FALSE))
  res
}

# a nonlinear constraint lb <= fun(x) <= ub (NonlinearConstraint)
.sp_nonlinear_constraint <- function(fun, lb, ub, jac = "2-point", keep_feasible = FALSE,
                                     finite_diff_rel_step = NULL) {
  list(type = "nonlinear", fun = fun, lb = as.numeric(lb), ub = as.numeric(ub), jac = jac,
       keep_feasible = keep_feasible, finite_diff_rel_step = finite_diff_rel_step)
}

# product of a Hessian (matrix or quasi Newton update) with a vector
.sp_hdot <- function(H, p) if (inherits(H, "SpHessianUpdate")) as.numeric(H$dot(p)) else as.numeric(H %*% p)

# the SciPy algorithms whose native port is still missing raise this error
.scipy_need <- function(fname, label) {
  if (!exists(fname, mode = "function", envir = asNamespace("modelcompR")))
    stop(label, " is not yet ported to native R in modelcompR", call. = FALSE)
  invisible(TRUE)
}
