# The unconstrained methods of scipy.optimize.minimize and the scalar
# minimizers, ported from scipy/optimize/_optimize.py (SciPy 1.18.1):
# Nelder-Mead, Powell, CG, BFGS, Newton-CG, Brent, golden section,
# bounded (Brent's fminbound), bracket() and brute().

# function wrapper that raises the "maxfun" condition after maxfun calls
.sp_maxfun_wrapper <- function(func, maxfun) {
  cnt <- new.env(parent = emptyenv())
  cnt$n <- 0L
  list(count = cnt, f = function(x) {
    if (cnt$n >= maxfun) .sp_stop("sp_maxfun", "Too many function calls")
    cnt$n <- cnt$n + 1L
    as.numeric(func(x))
  })
}

.sp_clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

.sp_neldermead <- function(func, x0, maxiter = NULL, maxfev = NULL, initial_simplex = NULL, xatol = 1e-4,
                           fatol = 1e-4, adaptive = FALSE, bounds = NULL, ...) {
  maxfun <- maxfev
  x0 <- as.numeric(x0)
  if (adaptive) {
    dim <- length(x0)
    rho <- 1; chi <- 1 + 2 / dim; psi <- 0.75 - 1 / (2 * dim); sigma <- 1 - 1 / dim
  } else {
    rho <- 1; chi <- 2; psi <- 0.5; sigma <- 0.5
  }
  nonzdelt <- 0.05
  zdelt <- 0.00025
  if (!is.null(bounds)) {
    lo <- bounds$lb; hi <- bounds$ub
    if (any(lo > hi)) stop("Nelder Mead - one of the lower bounds is greater than an upper bound.", call. = FALSE)
    x0 <- .sp_clip(x0, lo, hi)
  }
  if (is.null(initial_simplex)) {
    N <- length(x0)
    sim <- matrix(x0, N + 1L, N, byrow = TRUE)
    for (k in seq_len(N)) {
      y <- x0
      y[k] <- if (y[k] != 0) (1 + nonzdelt) * y[k] else zdelt
      sim[k + 1L, ] <- y
    }
  } else {
    sim <- as.matrix(initial_simplex)
    N <- ncol(sim)
  }
  if (is.null(maxiter) && is.null(maxfun)) {
    maxiter <- N * 200; maxfun <- N * 200
  } else if (is.null(maxiter)) {
    maxiter <- if (maxfun == Inf) N * 200 else Inf
  } else if (is.null(maxfun)) {
    maxfun <- if (maxiter == Inf) N * 200 else Inf
  }
  if (!is.null(bounds)) {
    up <- matrix(hi, N + 1L, N, byrow = TRUE)
    lw <- matrix(lo, N + 1L, N, byrow = TRUE)
    sim <- ifelse(sim > up, 2 * up - sim, sim)
    sim <- pmin(pmax(sim, lw), up)
  }
  w <- .sp_maxfun_wrapper(func, maxfun)
  fn <- w$f
  fcalls <- w$count
  fsim <- rep(Inf, N + 1L)
  tryCatch(for (k in seq_len(N + 1L)) fsim[k] <- fn(sim[k, ]), sp_maxfun = function(e) NULL)
  ord <- order(fsim, method = "radix")
  sim <- sim[ord, , drop = FALSE]
  fsim <- fsim[ord]
  iterations <- 1
  cl <- function(x) if (is.null(bounds)) x else .sp_clip(x, lo, hi)
  while (fcalls$n < maxfun && iterations < maxiter) {
    stop_now <- tryCatch({
      if (max(abs(sweep(sim[-1L, , drop = FALSE], 2L, sim[1L, ]))) <= xatol &&
          max(abs(fsim[1L] - fsim[-1L])) <= fatol) {
        TRUE
      } else {
        xbar <- colSums(sim[-(N + 1L), , drop = FALSE]) / N
        xr <- cl((1 + rho) * xbar - rho * sim[N + 1L, ])
        fxr <- fn(xr)
        doshrink <- FALSE
        if (fxr < fsim[1L]) {
          xe <- cl((1 + rho * chi) * xbar - rho * chi * sim[N + 1L, ])
          fxe <- fn(xe)
          if (fxe < fxr) {
            sim[N + 1L, ] <- xe; fsim[N + 1L] <- fxe
          } else {
            sim[N + 1L, ] <- xr; fsim[N + 1L] <- fxr
          }
        } else {
          if (fxr < fsim[N]) {
            sim[N + 1L, ] <- xr; fsim[N + 1L] <- fxr
          } else {
            if (fxr < fsim[N + 1L]) {
              xc <- cl((1 + psi * rho) * xbar - psi * rho * sim[N + 1L, ])
              fxc <- fn(xc)
              if (fxc <= fxr) {
                sim[N + 1L, ] <- xc; fsim[N + 1L] <- fxc
              } else {
                doshrink <- TRUE
              }
            } else {
              xcc <- cl((1 - psi) * xbar + psi * sim[N + 1L, ])
              fxcc <- fn(xcc)
              if (fxcc < fsim[N + 1L]) {
                sim[N + 1L, ] <- xcc; fsim[N + 1L] <- fxcc
              } else {
                doshrink <- TRUE
              }
            }
            if (doshrink) {
              for (j in 2:(N + 1L)) {
                sim[j, ] <- cl(sim[1L, ] + sigma * (sim[j, ] - sim[1L, ]))
                fsim[j] <- fn(sim[j, ])
              }
            }
          }
        }
        iterations <- iterations + 1
        FALSE
      }
    }, sp_maxfun = function(e) FALSE)
    if (isTRUE(stop_now)) break
    ord <- order(fsim, method = "radix")
    sim <- sim[ord, , drop = FALSE]
    fsim <- fsim[ord]
  }
  x <- sim[1L, ]
  fval <- min(fsim)
  if (fcalls$n >= maxfun) {
    warnflag <- 1L; msg <- .sp_status_message$maxfev
  } else if (iterations >= maxiter) {
    warnflag <- 2L; msg <- .sp_status_message$maxiter
  } else {
    warnflag <- 0L; msg <- .sp_status_message$success
  }
  list(x = x, fun = fval, nit = as.integer(iterations), nfev = fcalls$n, status = warnflag, success = warnflag == 0L,
       message = msg, final_simplex = list(sim, fsim))
}

# _line_for_search
.sp_line_for_search <- function(x0, alpha, lower_bound, upper_bound) {
  nz <- which(alpha != 0)
  lower_bound <- lower_bound[nz]; upper_bound <- upper_bound[nz]
  x0 <- x0[nz]; alpha <- alpha[nz]
  low <- (lower_bound - x0) / alpha
  high <- (upper_bound - x0) / alpha
  pos <- alpha > 0
  lmin <- max(ifelse(pos, low, 0) + ifelse(pos, 0, high))
  lmax <- min(ifelse(pos, high, 0) + ifelse(pos, 0, low))
  if (lmax >= lmin) c(lmin, lmax) else c(0, 0)
}

.sp_linesearch_powell <- function(func, p, xi, tol = 1e-3, lower_bound = NULL, upper_bound = NULL, fval = NULL) {
  myfunc <- function(alpha) func(p + alpha * xi)
  if (!any(xi != 0)) return(list(fval = if (!is.null(fval)) fval else func(p), x = p, xi = xi))
  if (is.null(lower_bound) && is.null(upper_bound)) {
    res <- .sp_recover_from_bracket_error(.sp_scalar_brent, myfunc, NULL, xtol = tol)
    xi <- res$x * xi
    return(list(fval = res$fun, x = p + xi, xi = xi))
  }
  bound <- .sp_line_for_search(p, xi, lower_bound, upper_bound)
  if (bound[1L] == -Inf && bound[2L] == Inf) {
    return(.sp_linesearch_powell(func, p, xi, fval = fval, tol = tol))
  } else if (bound[1L] != -Inf && bound[2L] != Inf) {
    res <- .sp_scalar_bounded(myfunc, bound, xatol = tol / 100)
    xi <- res$x * xi
    return(list(fval = res$fun, x = p + xi, xi = xi))
  }
  bound <- atan(bound)
  res <- .sp_scalar_bounded(function(x) myfunc(tan(x)), bound, xatol = tol / 100)
  xi <- tan(res$x) * xi
  list(fval = res$fun, x = p + xi, xi = xi)
}

.sp_powell <- function(func, x0, bounds = NULL, xtol = 1e-4, ftol = 1e-4, maxiter = NULL, maxfev = NULL,
                       direc = NULL, ...) {
  maxfun <- maxfev
  x <- as.numeric(x0)
  N <- length(x)
  if (is.null(maxiter) && is.null(maxfun)) {
    maxiter <- N * 1000; maxfun <- N * 1000
  } else if (is.null(maxiter)) {
    maxiter <- if (maxfun == Inf) N * 1000 else Inf
  } else if (is.null(maxfun)) {
    maxfun <- if (maxiter == Inf) N * 1000 else Inf
  }
  w <- .sp_maxfun_wrapper(func, maxfun)
  fn <- w$f
  fcalls <- w$count
  direc <- if (is.null(direc)) diag(1, N) else as.matrix(direc)
  lower_bound <- bounds$lb; upper_bound <- bounds$ub
  fval <- fn(x)
  x1 <- x
  iter <- 0
  repeat {
    brk <- tryCatch({
      fx <- fval
      bigind <- 1L
      delta <- 0.0
      for (i in seq_len(N)) {
        fx2 <- fval
        r <- .sp_linesearch_powell(fn, x, direc[i, ], tol = xtol * 100, lower_bound = lower_bound,
                                   upper_bound = upper_bound, fval = fval)
        fval <- r$fval; x <- r$x
        if ((fx2 - fval) > delta) {
          delta <- fx2 - fval
          bigind <- i
        }
      }
      iter <- iter + 1
      bnd <- ftol * (abs(fx) + abs(fval)) + 1e-20
      if (2.0 * (fx - fval) <= bnd || fcalls$n >= maxfun || iter >= maxiter || (is.nan(fx) && is.nan(fval))) {
        TRUE
      } else {
        direc1 <- x - x1
        x1 <- x
        lmax <- if (is.null(lower_bound) && is.null(upper_bound)) 1 else
          .sp_line_for_search(x, direc1, lower_bound, upper_bound)[2L]
        x2 <- x + min(lmax, 1) * direc1
        fx2 <- fn(x2)
        if (fx > fx2) {
          t <- 2.0 * (fx + fx2 - 2.0 * fval)
          temp <- (fx - fval - delta)
          t <- t * (temp * temp)
          temp <- fx - fx2
          t <- t - delta * temp * temp
          if (t < 0.0) {
            r <- .sp_linesearch_powell(fn, x, direc1, tol = xtol * 100, lower_bound = lower_bound,
                                       upper_bound = upper_bound, fval = fval)
            fval <- r$fval; x <- r$x
            if (any(r$xi != 0)) {
              direc[bigind, ] <- direc[N, ]
              direc[N, ] <- r$xi
            }
          }
        }
        FALSE
      }
    }, sp_maxfun = function(e) TRUE)
    if (brk) break
  }
  warnflag <- 0L
  msg <- .sp_status_message$success
  if (!is.null(bounds) && (any(lower_bound > x) || any(x > upper_bound))) {
    warnflag <- 4L; msg <- .sp_status_message$out_of_bounds
  } else if (fcalls$n >= maxfun) {
    warnflag <- 1L; msg <- .sp_status_message$maxfev
  } else if (iter >= maxiter) {
    warnflag <- 2L; msg <- .sp_status_message$maxiter
  } else if (is.nan(fval) || any(is.nan(x))) {
    warnflag <- 3L; msg <- .sp_status_message$nan
  }
  list(x = x, fun = fval, direc = direc, nit = as.integer(iter), nfev = fcalls$n, status = warnflag,
       success = warnflag == 0L, message = msg)
}

.sp_vecnorm <- function(x, ord = 2) {
  if (ord == Inf) max(abs(x)) else if (ord == -Inf) min(abs(x)) else sum(abs(x)^ord)^(1.0 / ord)
}

.sp_bfgs <- function(fun, x0, jac = NULL, gtol = 1e-5, norm = Inf, eps = sqrt(.Machine$double.eps), maxiter = NULL,
                     finite_diff_rel_step = NULL, xrtol = 0, c1 = 1e-4, c2 = 0.9, hess_inv0 = NULL, ...) {
  x0 <- as.numeric(x0)
  if (is.null(maxiter)) maxiter <- length(x0) * 200
  sf <- .sp_prepare_scalar_function(fun, x0, jac, epsilon = eps, finite_diff_rel_step = finite_diff_rel_step)
  f <- sf$fun
  myfprime <- sf$grad
  old_fval <- f(x0)
  gfk <- myfprime(x0)
  k <- 0L
  N <- length(x0)
  I <- diag(1, N)
  Hk <- if (is.null(hess_inv0)) I else hess_inv0
  old_old_fval <- old_fval + sqrt(sum(gfk * gfk)) / 2
  xk <- x0
  warnflag <- 0L
  gnorm <- .sp_vecnorm(gfk, norm)
  while (gnorm > gtol && k < maxiter) {
    pk <- -as.numeric(Hk %*% gfk)
    ls <- .sp_line_search_wolfe12(f, myfprime, xk, pk, gfk, old_fval, old_old_fval, amin = 0x1.bff2ee48e0530p-333, amax = 0x1.249ad2594c37dp+332,
                                  c1 = c1, c2 = c2)
    if (is.null(ls$alpha)) {
      warnflag <- 2L
      break
    }
    alpha_k <- ls$alpha
    old_fval <- ls$fval; old_old_fval <- ls$old_fval
    gfkp1 <- ls$gval
    sk <- alpha_k * pk
    xkp1 <- xk + sk
    xk <- xkp1
    if (is.null(gfkp1)) gfkp1 <- myfprime(xkp1)
    yk <- gfkp1 - gfk
    gfk <- gfkp1
    k <- k + 1L
    gnorm <- .sp_vecnorm(gfk, norm)
    if (gnorm <= gtol) break
    if (alpha_k * .sp_vecnorm(pk) <= xrtol * (xrtol + .sp_vecnorm(xk))) break
    if (!is.finite(old_fval)) {
      warnflag <- 2L
      break
    }
    rhok_inv <- sum(yk * sk)
    rhok <- if (rhok_inv == 0) 1000.0 else 1 / rhok_inv
    A1 <- I - outer(sk, yk) * rhok
    A2 <- I - outer(yk, sk) * rhok
    Hk <- A1 %*% (Hk %*% A2) + (rhok * outer(sk, sk))
  }
  fval <- old_fval
  if (warnflag == 2L) {
    msg <- .sp_status_message$pr_loss
  } else if (k >= maxiter) {
    warnflag <- 1L; msg <- .sp_status_message$maxiter
  } else if (is.nan(gnorm) || is.nan(fval) || any(is.nan(xk))) {
    warnflag <- 3L; msg <- .sp_status_message$nan
  } else {
    msg <- .sp_status_message$success
  }
  list(x = xk, fun = fval, jac = gfk, hess_inv = Hk, nfev = sf$nfev(), njev = sf$ngev, status = warnflag,
       success = warnflag == 0L, message = msg, nit = k)
}

.sp_cg <- function(fun, x0, jac = NULL, gtol = 1e-5, norm = Inf, eps = sqrt(.Machine$double.eps), maxiter = NULL,
                   finite_diff_rel_step = NULL, c1 = 1e-4, c2 = 0.4, ...) {
  x0 <- as.numeric(x0)
  if (is.null(maxiter)) maxiter <- length(x0) * 200
  sf <- .sp_prepare_scalar_function(fun, x0, jac, epsilon = eps, finite_diff_rel_step = finite_diff_rel_step)
  f <- sf$fun
  myfprime <- sf$grad
  old_fval <- f(x0)
  gfk <- myfprime(x0)
  k <- 0L
  xk <- x0
  old_old_fval <- old_fval + sqrt(sum(gfk * gfk)) / 2
  warnflag <- 0L
  pk <- -gfk
  gnorm <- .sp_vecnorm(gfk, norm)
  sigma_3 <- 0.01
  while (gnorm > gtol && k < maxiter) {
    deltak <- sum(gfk * gfk)
    cached <- new.env(parent = emptyenv())
    cached$step <- NULL
    prp_step <- function(alpha, gfkp1 = NULL) {
      xkp1 <- xk + alpha * pk
      if (is.null(gfkp1)) gfkp1 <- myfprime(xkp1)
      yk <- gfkp1 - gfk
      beta_k <- max(0, sum(yk * gfkp1) / deltak)
      pkp1 <- -gfkp1 + beta_k * pk
      list(alpha = alpha, xk = xkp1, pk = pkp1, gfk = gfkp1, gnorm = .sp_vecnorm(gfkp1, norm))
    }
    descent_condition <- function(alpha, xkp1, fp1, gfkp1) {
      cached$step <- prp_step(alpha, gfkp1)
      st <- cached$step
      if (st$gnorm <= gtol) return(TRUE)
      sum(st$pk * st$gfk) <= -sigma_3 * sum(st$gfk * st$gfk)
    }
    ls <- .sp_line_search_wolfe12(f, myfprime, xk, pk, gfk, old_fval, old_old_fval, c1 = c1, c2 = c2,
                                  amin = 0x1.bff2ee48e0530p-333, amax = 0x1.249ad2594c37dp+332, extra_condition = descent_condition)
    if (is.null(ls$alpha)) {
      warnflag <- 2L
      break
    }
    old_fval <- ls$fval; old_old_fval <- ls$old_fval
    st <- if (!is.null(cached$step) && ls$alpha == cached$step$alpha) cached$step else prp_step(ls$alpha, ls$gval)
    xk <- st$xk; pk <- st$pk; gfk <- st$gfk; gnorm <- st$gnorm
    k <- k + 1L
  }
  fval <- old_fval
  if (warnflag == 2L) {
    msg <- .sp_status_message$pr_loss
  } else if (k >= maxiter) {
    warnflag <- 1L; msg <- .sp_status_message$maxiter
  } else if (is.nan(gnorm) || is.nan(fval) || any(is.nan(xk))) {
    warnflag <- 3L; msg <- .sp_status_message$nan
  } else {
    msg <- .sp_status_message$success
  }
  list(x = xk, fun = fval, jac = gfk, nfev = sf$nfev(), njev = sf$ngev, status = warnflag, success = warnflag == 0L,
       message = msg, nit = k)
}

.sp_newtoncg <- function(fun, x0, jac = NULL, hess = NULL, hessp = NULL, xtol = 1e-5,
                         eps = sqrt(.Machine$double.eps), maxiter = NULL, c1 = 1e-4, c2 = 0.9, ...) {
  if (is.null(jac)) stop("Jacobian is required for Newton-CG method", call. = FALSE)
  x0 <- as.numeric(x0)
  sf <- .sp_prepare_scalar_function(fun, x0, jac, epsilon = eps, hess = hess)
  f <- sf$fun
  fprime <- sf$grad
  sf$hess(x0)
  fhess <- hess
  fhess_p <- hessp
  hcalls <- 0L
  if (is.null(maxiter)) maxiter <- length(x0) * 200
  cg_maxiter <- 20 * length(x0)
  xtol <- length(x0) * xtol
  update_l1norm <- .Machine$double.xmax
  xk <- x0
  k <- 0L
  gfk <- NULL
  old_fval <- f(x0)
  old_old_fval <- NULL
  terminate <- function(warnflag, msg)
    list(x = xk, fun = old_fval, jac = gfk, nfev = sf$nfev(), njev = sf$ngev, nhev = hcalls, status = warnflag,
         success = warnflag == 0L, message = msg, nit = k)
  while (update_l1norm > xtol) {
    if (k >= maxiter) return(terminate(1L, paste("Warning:", .sp_status_message$maxiter)))
    b <- -fprime(xk)
    maggrad <- sum(abs(b))
    eta <- min(0.5, sqrt(maggrad))
    termcond <- eta * maggrad
    xsupi <- numeric(length(x0))
    ri <- -b
    psupi <- -ri
    i <- 0L
    dri0 <- sum(ri * ri)
    if (!is.null(fhess)) {
      A <- sf$hess(xk)
      hcalls <- hcalls + 1L
    }
    converged <- FALSE
    for (k2 in seq_len(cg_maxiter)) {
      if (sum(abs(ri)) <= termcond) {
        converged <- TRUE
        break
      }
      if (is.null(fhess)) {
        if (is.null(fhess_p)) {
          stop("Newton-CG needs a Hessian or a Hessian vector product", call. = FALSE)
        } else {
          Ap <- fhess_p(xk, psupi)
          hcalls <- hcalls + 1L
        }
      } else {
        Ap <- .sp_hdot(A, psupi)
      }
      curv <- sum(psupi * Ap)
      if (0 <= curv && curv <= 3 * .sp_eps) {
        converged <- TRUE
        break
      } else if (curv < 0) {
        if (i > 0L) {
          converged <- TRUE
          break
        } else {
          xsupi <- dri0 / (-curv) * b
          converged <- TRUE
          break
        }
      }
      alphai <- dri0 / curv
      xsupi <- xsupi + alphai * psupi
      ri <- ri + alphai * Ap
      dri1 <- sum(ri * ri)
      betai <- dri1 / dri0
      psupi <- -ri + betai * psupi
      i <- i + 1L
      dri0 <- dri1
    }
    if (!converged) return(terminate(3L, "Warning: CG iterations didn't converge. The Hessian is not positive definite."))
    pk <- xsupi
    gfk <- -b
    ls <- .sp_line_search_wolfe12(f, fprime, xk, pk, gfk, old_fval, old_old_fval, c1 = c1, c2 = c2)
    if (is.null(ls$alpha)) return(terminate(2L, paste("Warning:", .sp_status_message$pr_loss)))
    old_fval <- ls$fval; old_old_fval <- ls$old_fval
    update <- ls$alpha * pk
    xk <- xk + update
    k <- k + 1L
    update_l1norm <- sum(abs(update))
  }
  if (is.nan(old_fval) || is.nan(update_l1norm)) return(terminate(3L, .sp_status_message$nan))
  terminate(0L, .sp_status_message$success)
}

# ---------------------------------------------------------------------------
# Scalar minimization

.sp_bracket <- function(func, xa = 0.0, xb = 1.0, grow_limit = 110.0, maxiter = 1000L) {
  gold <- 1.618034
  verysmall <- 1e-21
  fa <- func(xa)
  fb <- func(xb)
  if (fa < fb) {
    tmp <- xa; xa <- xb; xb <- tmp
    tmp <- fa; fa <- fb; fb <- tmp
  }
  xc <- xb + gold * (xb - xa)
  fc <- func(xc)
  funcalls <- 3L
  iter <- 0L
  while (fc < fb) {
    tmp1 <- (xb - xa) * (fb - fc)
    tmp2 <- (xb - xc) * (fb - fa)
    val <- tmp2 - tmp1
    denom <- if (abs(val) < verysmall) 2.0 * verysmall else 2.0 * val
    w <- xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom
    wlim <- xb + grow_limit * (xc - xb)
    if (iter > maxiter)
      stop("No valid bracket was found before the iteration limit was reached. Consider trying different ",
           "initial points or increasing `maxiter`.", call. = FALSE)
    iter <- iter + 1L
    done <- FALSE
    if ((w - xc) * (xb - w) > 0.0) {
      fw <- func(w)
      funcalls <- funcalls + 1L
      if (fw < fc) {
        xa <- xb; xb <- w; fa <- fb; fb <- fw
        done <- TRUE
      } else if (fw > fb) {
        xc <- w; fc <- fw
        done <- TRUE
      } else {
        w <- xc + gold * (xc - xb)
        fw <- func(w)
        funcalls <- funcalls + 1L
      }
    } else if ((w - wlim) * (wlim - xc) >= 0.0) {
      w <- wlim
      fw <- func(w)
      funcalls <- funcalls + 1L
    } else if ((w - wlim) * (xc - w) > 0.0) {
      fw <- func(w)
      funcalls <- funcalls + 1L
      if (fw < fc) {
        xb <- xc; xc <- w
        w <- xc + gold * (xc - xb)
        fb <- fc; fc <- fw
        fw <- func(w)
        funcalls <- funcalls + 1L
      }
    } else {
      w <- xc + gold * (xc - xb)
      fw <- func(w)
      funcalls <- funcalls + 1L
    }
    if (done) break
    xa <- xb; xb <- xc; xc <- w
    fa <- fb; fb <- fc; fc <- fw
  }
  cond1 <- (fb < fc && fb <= fa) || (fb < fa && fb <= fc)
  cond2 <- (xa < xb && xb < xc) || (xc < xb && xb < xa)
  cond3 <- is.finite(xa) && is.finite(xb) && is.finite(xc)
  out <- list(xa = xa, xb = xb, xc = xc, fa = fa, fb = fb, fc = fc, funcalls = funcalls)
  if (!(cond1 && cond2 && cond3))
    .sp_stop("sp_bracket_error", paste("The algorithm terminated without finding a valid bracket.",
                                       "Consider trying different initial points."), out)
  out
}

.sp_get_bracket <- function(func, brack) {
  if (is.null(brack)) return(.sp_bracket(func))
  if (length(brack) == 2L) return(.sp_bracket(func, xa = brack[1L], xb = brack[2L]))
  if (length(brack) == 3L) {
    xa <- brack[1L]; xb <- brack[2L]; xc <- brack[3L]
    if (xa > xc) { tmp <- xa; xa <- xc; xc <- tmp }
    if (!(xa < xb && xb < xc))
      stop("Bracketing values (xa, xb, xc) do not fulfill this requirement: (xa < xb) and (xb < xc)", call. = FALSE)
    fa <- func(xa); fb <- func(xb); fc <- func(xc)
    if (!(fb < fa && fb < fc))
      stop("Bracketing values (xa, xb, xc) do not fulfill this requirement: (f(xb) < f(xa)) and (f(xb) < f(xc))",
           call. = FALSE)
    return(list(xa = xa, xb = xb, xc = xc, fa = fa, fb = fb, fc = fc, funcalls = 3L))
  }
  stop("Bracketing interval must be length 2 or 3 sequence.", call. = FALSE)
}

.sp_recover_from_bracket_error <- function(solver, fun, bracket, ...) {
  tryCatch(solver(fun, bracket, ...), sp_bracket_error = function(e) {
    d <- e$data
    xs <- c(d$xa, d$xb, d$xc); fs <- c(d$fa, d$fb, d$fc)
    if (anyNA(c(xs, fs))) {
      x <- NaN; fv <- NaN
    } else {
      i <- which.min(fs)
      x <- xs[i]; fv <- fs[i]
    }
    list(fun = fv, nfev = d$funcalls, x = x, nit = 0L, success = FALSE, message = conditionMessage(e))
  })
}

# Brent.optimize of _optimize.py
.sp_scalar_brent <- function(func, brack = NULL, xtol = 1.48e-8, maxiter = 500L, ...) {
  tol <- xtol
  if (tol < 0) stop("tolerance should be >= 0, got ", tol, call. = FALSE)
  br <- .sp_get_bracket(func, brack)
  xa <- br$xa; xb <- br$xb; xc <- br$xc; funcalls <- br$funcalls
  mintol <- 1.0e-11
  cg <- 0.3819660
  x <- w <- v <- xb
  fw <- fv <- fx <- br$fb
  if (xa < xc) { a <- xa; b <- xc } else { a <- xc; b <- xa }
  deltax <- 0.0
  rat <- 0
  iter <- 0L
  while (iter < maxiter) {
    tol1 <- tol * abs(x) + mintol
    tol2 <- 2.0 * tol1
    xmid <- 0.5 * (a + b)
    if (abs(x - xmid) < (tol2 - 0.5 * (b - a))) break
    if (abs(deltax) <= tol1) {
      deltax <- if (x >= xmid) a - x else b - x
      rat <- cg * deltax
    } else {
      tmp1 <- (x - w) * (fx - fv)
      tmp2 <- (x - v) * (fx - fw)
      p <- (x - v) * tmp2 - (x - w) * tmp1
      tmp2 <- 2.0 * (tmp2 - tmp1)
      if (tmp2 > 0.0) p <- -p
      tmp2 <- abs(tmp2)
      dx_temp <- deltax
      deltax <- rat
      if (p > tmp2 * (a - x) && p < tmp2 * (b - x) && abs(p) < abs(0.5 * tmp2 * dx_temp)) {
        rat <- p * 1.0 / tmp2
        u <- x + rat
        if ((u - a) < tol2 || (b - u) < tol2) rat <- if (xmid - x >= 0) tol1 else -tol1
      } else {
        deltax <- if (x >= xmid) a - x else b - x
        rat <- cg * deltax
      }
    }
    u <- if (abs(rat) < tol1) (if (rat >= 0) x + tol1 else x - tol1) else x + rat
    fu <- func(u)
    funcalls <- funcalls + 1L
    if (fu > fx) {
      if (u < x) a <- u else b <- u
      if (fu <= fw || w == x) {
        v <- w; w <- u; fv <- fw; fw <- fu
      } else if (fu <= fv || v == x || v == w) {
        v <- u; fv <- fu
      }
    } else {
      if (u >= x) a <- x else b <- x
      v <- w; w <- x; x <- u
      fv <- fw; fw <- fx; fx <- fu
    }
    iter <- iter + 1L
  }
  success <- iter < maxiter && !(is.nan(x) || is.nan(fx))
  msg <- if (success) "Optimization terminated successfully" else if (iter >= maxiter)
    "Maximum number of iterations exceeded" else .sp_status_message$nan
  list(fun = fx, x = x, nit = iter, nfev = funcalls, success = success, message = msg)
}

.sp_scalar_golden <- function(func, brack = NULL, xtol = sqrt(.Machine$double.eps), maxiter = 5000L, ...) {
  tol <- xtol
  br <- .sp_get_bracket(func, brack)
  xa <- br$xa; xb <- br$xb; xc <- br$xc; funcalls <- br$funcalls
  gR <- 0.61803399
  gC <- 1.0 - gR
  x3 <- xc
  x0 <- xa
  if (abs(xc - xb) > abs(xb - xa)) {
    x1 <- xb
    x2 <- xb + gC * (xc - xb)
  } else {
    x2 <- xb
    x1 <- xb - gC * (xb - xa)
  }
  f1 <- func(x1)
  f2 <- func(x2)
  funcalls <- funcalls + 2L
  nit <- 0L
  for (i in seq_len(maxiter)) {
    if (abs(x3 - x0) <= tol * (abs(x1) + abs(x2))) break
    if (f2 < f1) {
      x0 <- x1; x1 <- x2
      x2 <- gR * x1 + gC * x3
      f1 <- f2
      f2 <- func(x2)
    } else {
      x3 <- x2; x2 <- x1
      x1 <- gR * x2 + gC * x0
      f2 <- f1
      f1 <- func(x1)
    }
    funcalls <- funcalls + 1L
    nit <- nit + 1L
  }
  if (f1 < f2) { xmin <- x1; fval <- f1 } else { xmin <- x2; fval <- f2 }
  success <- nit < maxiter && !(is.nan(fval) || is.nan(xmin))
  msg <- if (success) "Optimization terminated successfully" else if (nit >= maxiter)
    "Maximum number of iterations exceeded" else .sp_status_message$nan
  list(fun = fval, nfev = funcalls, x = xmin, nit = nit, success = success, message = msg)
}

# _minimize_scalar_bounded (fminbound)
.sp_scalar_bounded <- function(func, bounds, xatol = 1e-5, maxiter = 500L, ...) {
  maxfun <- maxiter
  x1 <- bounds[1L]; x2 <- bounds[2L]
  if (!(is.finite(x1) && is.finite(x2))) stop("Optimization bounds must be finite scalars.", call. = FALSE)
  if (x1 > x2) stop("The lower bound exceeds the upper bound.", call. = FALSE)
  flag <- 0L
  sqrt_eps <- sqrt(2.2e-16)
  golden_mean <- 0.5 * (3.0 - sqrt(5.0))
  a <- x1; b <- x2
  fulc <- a + golden_mean * (b - a)
  nfc <- xf <- fulc
  rat <- e <- 0.0
  x <- xf
  fx <- func(x)
  num <- 1L
  fu <- Inf
  ffulc <- fnfc <- fx
  xm <- 0.5 * (a + b)
  tol1 <- sqrt_eps * abs(xf) + xatol / 3.0
  tol2 <- 2.0 * tol1
  while (abs(xf - xm) > (tol2 - 0.5 * (b - a))) {
    golden <- TRUE
    if (abs(e) > tol1) {
      golden <- FALSE
      r <- (xf - nfc) * (fx - ffulc)
      q <- (xf - fulc) * (fx - fnfc)
      p <- (xf - fulc) * q - (xf - nfc) * r
      q <- 2.0 * (q - r)
      if (q > 0.0) p <- -p
      q <- abs(q)
      r <- e
      e <- rat
      if (abs(p) < abs(0.5 * q * r) && p > q * (a - xf) && p < q * (b - xf)) {
        rat <- (p + 0.0) / q
        x <- xf + rat
        if ((x - a) < tol2 || (b - x) < tol2) {
          si <- sign(xm - xf) + ((xm - xf) == 0)
          rat <- tol1 * si
        }
      } else {
        golden <- TRUE
      }
    }
    if (golden) {
      e <- if (xf >= xm) a - xf else b - xf
      rat <- golden_mean * e
    }
    si <- sign(rat) + (rat == 0)
    x <- xf + si * max(abs(rat), tol1)
    fu <- func(x)
    num <- num + 1L
    if (fu <= fx) {
      if (x >= xf) a <- xf else b <- xf
      fulc <- nfc; ffulc <- fnfc
      nfc <- xf; fnfc <- fx
      xf <- x; fx <- fu
    } else {
      if (x < xf) a <- x else b <- x
      if (fu <= fnfc || nfc == xf) {
        fulc <- nfc; ffulc <- fnfc
        nfc <- x; fnfc <- fu
      } else if (fu <= ffulc || fulc == xf || fulc == nfc) {
        fulc <- x; ffulc <- fu
      }
    }
    xm <- 0.5 * (a + b)
    tol1 <- sqrt_eps * abs(xf) + xatol / 3.0
    tol2 <- 2.0 * tol1
    if (num >= maxfun) {
      flag <- 1L
      break
    }
  }
  if (is.nan(xf) || is.nan(fx) || is.nan(fu)) flag <- 2L
  msg <- c("Solution found.", "Maximum number of function calls reached.", .sp_status_message$nan)[flag + 1L]
  list(fun = fx, status = flag, success = flag == 0L, message = msg, x = xf, nfev = num, nit = num)
}

.sp_minimize_scalar <- function(fun, bracket = NULL, bounds = NULL, method = NULL, tol = NULL, options = list()) {
  meth <- if (is.null(method)) (if (is.null(bounds)) "brent" else "bounded") else tolower(method)
  if (!is.null(bounds) && meth %in% c("brent", "golden"))
    stop("Use of `bounds` is incompatible with 'method=", method, "'.", call. = FALSE)
  if (!is.null(tol)) {
    if (meth == "bounded" && is.null(options$xatol)) options$xatol <- tol
    else if (is.null(options$xtol)) options$xtol <- tol
  }
  res <- switch(meth,
    brent = do.call(.sp_recover_from_bracket_error, c(list(.sp_scalar_brent, fun, bracket), options)),
    bounded = {
      if (is.null(bounds)) stop("The `bounds` parameter is mandatory for method `bounded`.", call. = FALSE)
      do.call(.sp_scalar_bounded, c(list(fun, bounds), options))
    },
    golden = do.call(.sp_recover_from_bracket_error, c(list(.sp_scalar_golden, fun, bracket), options)),
    stop("Unknown solver ", method, call. = FALSE))
  res$status <- res$status %||% if (isTRUE(res$success)) 0L else 1L
  res
}

# brute() with finish = fmin (Nelder-Mead with xatol = fatol = 1e-4)
.sp_brute <- function(func, lb, ub, Ns = 20L, finish = TRUE, ...) {
  N <- length(lb)
  if (N > 40L) stop("Brute Force not possible with more than 40 variables.", call. = FALSE)
  axes <- lapply(seq_len(N), function(k) {
    step <- (ub[k] - lb[k]) / (Ns - 1)
    (0:(Ns - 1)) * step + lb[k]
  })
  # grid points in C order (the last coordinate varies fastest)
  G <- as.matrix(expand.grid(rev(axes)))[, N:1, drop = FALSE]
  Jout <- vapply(seq_len(nrow(G)), function(i) as.numeric(func(G[i, ])), 0)
  i <- which.min(Jout)
  xmin <- as.numeric(G[i, ])
  Jmin <- Jout[i]
  success <- TRUE
  nfev <- length(Jout)
  if (isTRUE(finish)) {
    res <- .sp_neldermead(func, xmin, xatol = 1e-4, fatol = 1e-4)
    xmin <- res$x
    Jmin <- res$fun
    success <- res$status == 0L
    nfev <- nfev + res$nfev
  }
  list(x = xmin, fun = Jmin, grid = G, Jout = Jout, nfev = nfev, success = success, retcode = "Success")
}

.scipy_run_brute <- function(f, lb, ub, kw) {
  kw$finish <- kw$finish %||% TRUE
  do.call(.sp_brute, c(list(function(x) as.numeric(f(x)), lb, ub), kw))
}
