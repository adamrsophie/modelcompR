# Scalar root finding of scipy.optimize (SciPy 1.18.1): root_scalar
# (_root_scalar.py), newton, secant, halley and toms748 (_zeros_py.py) and the
# C bracketing solvers bisect, brentq, brenth and ridder (Zeros/*.c).

.sp_xtol <- 2e-12
.sp_rtol <- 4 * .Machine$double.eps
.sp_signbit <- function(x) x < 0 || (x == 0 && 1 / x < 0)

# result tuple of the C solvers: list(root, funcalls, iterations, flag)
.sp_zeros_c <- function(kind, f, xa, xb, xtol, rtol, iter) {
  calls <- 0L
  fv <- function(x) {
    calls <<- calls + 1L
    v <- as.numeric(f(x))
    if (is.nan(v)) {
      e <- structure(class = c("sp_nan_root", "error", "condition"),
                     list(message = sprintf("The function value at x=%s is NaN; solver cannot continue.", format(x)),
                          call = NULL, x = x, calls = calls))
      stop(e)
    }
    v
  }
  out <- function(x, flag, its) list(root = x, function_calls = calls, iterations = its, flag = flag)
  signerr <- function() stop("f(a) and f(b) must have different signs", call. = FALSE)
  its <- 0L
  if (kind == "bisect") {
    fa <- fv(xa); fb <- fv(xb)
    if (fa == 0) return(out(xa, 0L, its))
    if (fb == 0) return(out(xb, 0L, its))
    if (.sp_signbit(fa) == .sp_signbit(fb)) signerr()
    dm <- xb - xa
    for (i in seq_len(iter)) {
      its <- its + 1L
      dm <- dm * 0.5
      xm <- xa + dm
      fm <- fv(xm)
      if (.sp_signbit(fm) == .sp_signbit(fa)) xa <- xm
      if (fm == 0 || abs(dm) < xtol + rtol * abs(xm)) return(out(xm, 0L, its))
    }
    return(out(xa, -2L, its))
  }
  if (kind == "ridder") {
    tol <- xtol + rtol * min(abs(xa), abs(xb))
    fa <- fv(xa); fb <- fv(xb)
    if (fa == 0) return(out(xa, 0L, its))
    if (fb == 0) return(out(xb, 0L, its))
    if (.sp_signbit(fa) == .sp_signbit(fb)) signerr()
    xn <- 0.0
    for (i in seq_len(iter)) {
      its <- its + 1L
      dm <- 0.5 * (xb - xa)
      xm <- xa + dm
      fm <- fv(xm)
      if (fm == 0.0) return(out(xm, 0L, its))
      ratio <- fm / fa
      dn <- dm * ratio / sqrt(ratio * ratio - fb / fa)
      xn <- xm + (if (dn > 0) 1 else -1) * min(abs(dn), abs(dm) - 0.5 * tol)
      fn <- fv(xn)
      if (.sp_signbit(fn) != .sp_signbit(fm)) {
        xa <- xn; fa <- fn; xb <- xm; fb <- fm
      } else if (.sp_signbit(fn) != .sp_signbit(fa)) {
        xb <- xn; fb <- fn
      } else {
        xa <- xn; fa <- fn
      }
      # the C code uses rtol * xn here, without an absolute value
      tol <- xtol + rtol * xn
      if (fn == 0.0 || abs(xb - xa) < tol) return(out(xn, 0L, its))
    }
    return(out(xn, -2L, its))
  }
  # brentq and brenth
  xpre <- xa; xcur <- xb
  xblk <- 0; fblk <- 0; spre <- 0; scur <- 0
  fpre <- fv(xpre); fcur <- fv(xcur)
  if (fpre == 0) return(out(xpre, 0L, its))
  if (fcur == 0) return(out(xcur, 0L, its))
  if (.sp_signbit(fpre) == .sp_signbit(fcur)) signerr()
  for (i in seq_len(iter)) {
    its <- its + 1L
    if (fpre != 0 && fcur != 0 && .sp_signbit(fpre) != .sp_signbit(fcur)) {
      xblk <- xpre
      fblk <- fpre
      spre <- scur <- xcur - xpre
    }
    if (abs(fblk) < abs(fcur)) {
      xpre <- xcur; xcur <- xblk; xblk <- xpre
      fpre <- fcur; fcur <- fblk; fblk <- fpre
    }
    delta <- (xtol + rtol * abs(xcur)) / 2
    sbis <- (xblk - xcur) / 2
    if (fcur == 0 || abs(sbis) < delta) return(out(xcur, 0L, its))
    if (abs(spre) > delta && abs(fcur) < abs(fpre)) {
      if (xpre == xblk) {
        stry <- -fcur * (xcur - xpre) / (fcur - fpre)
      } else {
        dpre <- (fpre - fcur) / (xpre - xcur)
        dblk <- (fblk - fcur) / (xblk - xcur)
        stry <- if (kind == "brentq") -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre)) else
          -fcur * (fblk - fpre) / (fblk * dpre - fpre * dblk)
      }
      if (2 * abs(stry) < min(abs(spre), 3 * abs(sbis) - delta)) {
        spre <- scur
        scur <- stry
      } else {
        spre <- sbis
        scur <- sbis
      }
    } else {
      spre <- sbis
      scur <- sbis
    }
    xpre <- xcur; fpre <- fcur
    if (abs(scur) > delta) xcur <- xcur + scur else xcur <- xcur + (if (sbis > 0) delta else -delta)
    fcur <- fv(xcur)
  }
  out(xcur, -2L, its)
}

.sp_isclose <- function(a, b, rtol, atol) abs(a - b) <= atol + rtol * abs(b)

# newton(), secant (fprime NULL) and halley (fprime2 given), full_output = TRUE, disp = FALSE
.sp_newton <- function(func, x0, fprime = NULL, tol = 1.48e-8, maxiter = 50L, fprime2 = NULL, x1 = NULL, rtol = 0.0) {
  if (tol <= 0) stop(sprintf("tol too small (%g <= 0)", tol), call. = FALSE)
  if (maxiter < 1) stop("maxiter must be greater than 0", call. = FALSE)
  p0 <- as.numeric(x0) * 1.0
  funcalls <- 0L
  res <- function(x, its, flag, method) list(root = x, function_calls = funcalls, iterations = its, flag = flag,
                                             method = method)
  if (!is.null(fprime)) {
    method <- "newton"
    for (itr in 0:(maxiter - 1L)) {
      fval <- func(p0)
      funcalls <- funcalls + 1L
      if (fval == 0) return(res(p0, itr, 0L, method))
      fder <- fprime(p0)
      funcalls <- funcalls + 1L
      if (fder == 0) return(res(p0, itr + 1L, -2L, method))
      newton_step <- fval / fder
      if (!is.null(fprime2)) {
        fder2 <- fprime2(p0)
        funcalls <- funcalls + 1L
        method <- "halley"
        adj <- newton_step * fder2 / fder / 2
        if (abs(adj) < 1) newton_step <- newton_step / (1.0 - adj)
      }
      p <- p0 - newton_step
      if (.sp_isclose(p, p0, rtol, tol)) return(res(p, itr + 1L, 0L, method))
      p0 <- p
    }
  } else {
    method <- "secant"
    if (!is.null(x1)) {
      if (x1 == x0) stop("x1 and x0 must be different", call. = FALSE)
      p1 <- x1
    } else {
      eps <- 1e-4
      p1 <- x0 * (1 + eps)
      p1 <- p1 + (if (p1 >= 0) eps else -eps)
    }
    q0 <- func(p0); q1 <- func(p1)
    funcalls <- funcalls + 2L
    if (abs(q1) < abs(q0)) {
      tmp <- p0; p0 <- p1; p1 <- tmp
      tmp <- q0; q0 <- q1; q1 <- tmp
    }
    for (itr in 0:(maxiter - 1L)) {
      if (q1 == q0) {
        p <- (p1 + p0) / 2.0
        return(res(p, itr + 1L, -2L, method))
      }
      p <- if (abs(q1) > abs(q0)) (-q0 / q1 * p1 + p0) / (1 - q0 / q1) else (-q1 / q0 * p0 + p1) / (1 - q1 / q0)
      if (.sp_isclose(p, p1, rtol, tol)) return(res(p, itr + 1L, 0L, method))
      p0 <- p1; q0 <- q1
      p1 <- p
      q1 <- func(p1)
      funcalls <- funcalls + 1L
    }
  }
  res(p, itr + 1L, -2L, method)
}

# --------------------------------------------------------------------------
# toms748 (Alefeld, Potra and Shi, TOMS 748)

.sp_notclose <- function(fs, rtol, atol) {
  if (any(vapply(fs, is.null, TRUE))) return(FALSE)
  fs <- unlist(fs)
  if (!all(fs != 0) || !all(is.finite(fs))) return(FALSE)
  for (i in seq_len(length(fs) - 1L)) if (any(.sp_isclose(fs[i], fs[(i + 1L):length(fs)], rtol, atol))) return(FALSE)
  TRUE
}

.sp_secant2 <- function(xv, fv) {
  if (fv[1L] == fv[2L]) return(NaN)
  if (abs(fv[2L]) > abs(fv[1L])) (-fv[1L] / fv[2L] * xv[2L] + xv[1L]) / (1 - fv[1L] / fv[2L]) else
    (-fv[2L] / fv[1L] * xv[1L] + xv[2L]) / (1 - fv[2L] / fv[1L])
}

.sp_divided_differences <- function(xvals, fvals, forward = TRUE) {
  dd <- fvals
  row <- fvals
  m <- length(xvals)
  dd[1L] <- if (forward) fvals[1L] else fvals[m]
  for (i in .lb_rng(2L, m)) {
    k <- length(row)
    denom <- xvals[i:(i + k - 2L)] - xvals[1:(k - 1L)]
    row <- diff(row) / denom
    dd[i] <- if (forward) row[1L] else row[length(row)]
  }
  dd
}

.sp_interpolated_poly <- function(xvals, fvals, x) {
  N <- length(xvals)
  Q <- matrix(0, N, N); D <- matrix(0, N, N)
  Q[, 1L] <- fvals; D[, 1L] <- fvals
  for (k in .lb_rng(1L, N - 1L)) {
    alpha <- D[(k + 1L):N, k] - Q[k:(N - 1L), k]
    diffik <- xvals[1:(N - k)] - xvals[(k + 1L):N]
    Q[(k + 1L):N, k + 1L] <- (xvals[(k + 1L):N] - x) / diffik * alpha
    D[(k + 1L):N, k + 1L] <- (xvals[1:(N - k)] - x) / diffik * alpha
  }
  sum(Q[N, 2:N]) + Q[N, 1L]
}

.sp_newton_quadratic <- function(ab, fab, d, fd, k) {
  a <- ab[1L]; b <- ab[2L]; fa <- fab[1L]
  dd <- .sp_divided_differences(c(a, b, d), c(fa, fab[2L], fd))
  B <- dd[2L]; A <- dd[3L]
  P <- function(x) (A * (x - b) + B) * (x - a) + fa
  if (A == 0) return(a - fa / B)
  r <- if (sign(A) * sign(fa) > 0) a else b
  for (i in seq_len(k)) {
    r1 <- r - P(r) / (B + A * (2 * r - a - b))
    if (!(ab[1L] < r1 && r1 < ab[2L])) {
      if (ab[1L] < r && r < ab[2L]) return(r)
      r <- sum(ab) / 2.0
      break
    }
    r <- r1
  }
  r
}

.sp_toms748 <- function(f, a, b, k = 1L, xtol = .sp_xtol, rtol = .sp_rtol, maxiter = 100L) {
  if (xtol <= 0) stop(sprintf("xtol too small (%g <= 0)", xtol), call. = FALSE)
  if (rtol < .sp_rtol / 4) stop(sprintf("rtol too small (%g < %g)", rtol, .sp_rtol / 4), call. = FALSE)
  if (maxiter < 1) stop("maxiter must be greater than 0", call. = FALSE)
  if (!is.finite(a)) stop("a is not finite ", a, call. = FALSE)
  if (!is.finite(b)) stop("b is not finite ", b, call. = FALSE)
  if (a >= b) stop(sprintf("a and b are not an interval [%s, %s]", format(a), format(b)), call. = FALSE)
  if (!(k >= 1)) stop(sprintf("k too small (%s < 1)", format(k)), call. = FALSE)
  kk <- min(max(k, 1L), 100L)
  calls <- 0L
  its <- 0L
  callf <- function(x) {
    fx <- as.numeric(f(x))
    calls <<- calls + 1L
    if (!is.finite(fx)) stop(sprintf("Invalid function value: f(%f) -> %s ", x, format(fx)), call. = FALSE)
    fx
  }
  result <- function(x, flag = 0L) list(root = x, function_calls = calls, iterations = its, flag = flag)
  eps <- .Machine$double.eps
  ab <- c(a, b)
  fa <- callf(a)
  if (fa == 0) return(result(a))
  fb <- callf(b)
  if (fb == 0) return(result(b))
  if (sign(fb) * sign(fa) > 0)
    stop(sprintf("f(a) and f(b) must have different signs, but f(%e)=%e, f(%e)=%e ", a, fa, b, fb), call. = FALSE)
  fab <- c(fa, fb)
  upd <- function(c, fc) {
    idx <- if (sign(fab[1L]) * sign(fc) > 0) 1L else 2L
    out <- list(x = ab[idx], f = fab[idx])
    fab[idx] <<- fc
    ab[idx] <<- c
    out
  }
  c <- .sp_secant2(ab, fab)
  if (!isTRUE(ab[1L] < c && c < ab[2L])) c <- sum(ab) / 2.0
  fc <- callf(c)
  if (fc == 0) return(result(c))
  r <- upd(c, fc)
  d <- r$x; fd <- r$f
  e <- NULL; fe <- NULL
  its <- its + 1L
  repeat {
    its <- its + 1L
    ab_width <- ab[2L] - ab[1L]
    c <- NULL
    for (nsteps in 2:(kk + 1L)) {
      if (.sp_notclose(list(fab[1L], fab[2L], fd, fe), 0, 32 * eps)) {
        c0 <- .sp_interpolated_poly(c(fab[1L], fab[2L], fd, fe), c(ab[1L], ab[2L], d, e), 0)
        if (ab[1L] < c0 && c0 < ab[2L]) c <- c0
      }
      if (is.null(c)) c <- .sp_newton_quadratic(ab, fab, d, fd, nsteps)
      fc <- callf(c)
      if (fc == 0) return(result(c))
      e <- d; fe <- fd
      r <- upd(c, fc)
      d <- r$x; fd <- r$f
    }
    uix <- if (abs(fab[1L]) < abs(fab[2L])) 1L else 2L
    u <- ab[uix]; fu <- fab[uix]
    A <- .sp_divided_differences(ab, fab, forward = uix == 1L)[2L]
    c <- u - 2 * fu / A
    if (abs(c - u) > 0.5 * (ab[2L] - ab[1L])) {
      c <- sum(ab) / 2.0
    } else if (.sp_isclose(c, u, eps, 0)) {
      frs <- floor(log2(abs(fab))) + 1
      if (frs[uix] < frs[3L - uix] - 50) {
        c <- (31 * ab[uix] + ab[3L - uix]) / 32
      } else {
        mm <- if (uix == 1L) 1 else -1
        c <- u + (mm * abs(c) * rtol + mm * xtol)
      }
      if (!(ab[1L] < c && c < ab[2L])) c <- sum(ab) / 2.0
    }
    fc <- callf(c)
    if (fc == 0) return(result(c))
    e <- d; fe <- fd
    r <- upd(c, fc)
    d <- r$x; fd <- r$f
    if (ab[2L] - ab[1L] > 0.5 * ab_width) {
      e <- d; fe <- fd
      z <- sum(ab) / 2.0
      fz <- callf(z)
      if (fz == 0) return(result(z))
      r <- upd(z, fz)
      d <- r$x; fd <- r$f
    }
    if (.sp_isclose(ab[1L], ab[2L], rtol, xtol)) return(result(sum(ab) / 2.0))
    if (its >= maxiter) return(result(sum(ab) / 2.0, -2L))
  }
}

# root_scalar(); the result has root, iterations, function_calls, converged and flag
.sp_root_scalar <- function(f, method = NULL, bracket = NULL, fprime = NULL, fprime2 = NULL, x0 = NULL, x1 = NULL,
                            xtol = NULL, rtol = NULL, maxiter = NULL, options = list()) {
  kw <- list(xtol = xtol, rtol = rtol, maxiter = maxiter)
  kw <- c(kw[!vapply(kw, is.null, TRUE)], options)
  if (is.null(method)) {
    method <- if (!is.null(bracket)) "brentq" else if (!is.null(x0)) {
      if (!is.null(fprime)) (if (!is.null(fprime2)) "halley" else "newton") else if (!is.null(x1)) "secant" else
        "newton"
    }
  }
  if (is.null(method)) stop("Unable to select a solver as neither bracket nor starting point provided.", call. = FALSE)
  meth <- tolower(method)
  flag_name <- function(flag) c("0" = "converged", "-1" = "sign error", "-2" = "convergence error")[[as.character(flag)]]
  if (meth %in% c("bisect", "ridder", "brentq", "brenth", "toms748")) {
    if (is.null(bracket)) stop("Bracket needed for ", method, call. = FALSE)
    a <- bracket[1L]; b <- bracket[2L]
    sol <- if (meth == "toms748") {
      do.call(.sp_toms748, c(list(f, a, b), kw))
    } else {
      xt <- kw$xtol %||% .sp_xtol; rt <- kw$rtol %||% .sp_rtol; mi <- as.integer(kw$maxiter %||% 100L)
      if (xt <= 0) stop(sprintf("xtol too small (%g <= 0)", xt), call. = FALSE)
      if (rt < .sp_rtol) stop(sprintf("rtol too small (%g < %g)", rt, .sp_rtol), call. = FALSE)
      tryCatch(.sp_zeros_c(meth, f, a, b, xt, rt, mi), sp_nan_root = function(e)
        list(root = e$x, iterations = NaN, function_calls = e$calls, flag = conditionMessage(e)))
    }
  } else if (meth %in% c("secant", "newton", "halley")) {
    if (is.null(x0)) stop("x0 must not be None for ", method, call. = FALSE)
    if (!is.null(kw$xtol)) { kw$tol <- kw$xtol; kw$xtol <- NULL }
    if (meth == "secant") {
      sol <- do.call(.sp_newton, c(list(f, x0, fprime = NULL, fprime2 = NULL, x1 = x1), kw))
    } else if (meth == "newton") {
      fp <- fprime %||% function(x) .sp_approx_derivative(function(z) f(z[1L]), x, "2-point")$J[1L]
      sol <- do.call(.sp_newton, c(list(f, x0, fprime = fp, fprime2 = NULL), kw))
    } else {
      if (is.null(fprime)) stop("fprime must be specified for ", method, call. = FALSE)
      if (is.null(fprime2)) stop("fprime2 must be specified for ", method, call. = FALSE)
      sol <- do.call(.sp_newton, c(list(f, x0, fprime = fprime, fprime2 = fprime2), kw))
    }
  } else {
    stop("Unknown solver ", method, call. = FALSE)
  }
  sol$converged <- identical(sol$flag, 0L)
  if (is.numeric(sol$flag)) sol$flag <- flag_name(sol$flag)
  sol$method <- meth
  sol
}

# ScipyRootScalar: a bracket from the bounds for the bracketing methods
# (toms748 included, see the porting notes), else x0 = u0
.scipy_run_root_scalar <- function(alg, f, grad, hess, x0, lb, ub, kw) {
  fs <- function(x) as.numeric(f(x))
  method <- alg$method
  args <- list(fs, method = method)
  if (method %in% c("brentq", "brenth", "bisect", "ridder", "toms748")) {
    if (is.null(lb) || is.null(ub)) stop("Method ", method, " requires bracket (bounds)", call. = FALSE)
    args$bracket <- c(lb[1L], ub[1L])
  } else {
    args$x0 <- x0[1L]
  }
  if (method == "newton" && !is.null(grad)) args$fprime <- function(x) grad(x)[1L]
  if (method == "halley") {
    if (is.null(grad) || is.null(hess)) stop("Method 'halley' requires both gradient and Hessian", call. = FALSE)
    args$fprime <- function(x) grad(x)[1L]
    args$fprime2 <- function(x) hess(x)[1L, 1L]
  }
  args <- utils::modifyList(args, kw)
  sol <- do.call(.sp_root_scalar, args)
  root <- sol$root
  list(x = root, fun = if (is.null(root) || is.nan(root)) NaN else abs(fs(root)), success = sol$converged,
       retcode = if (isTRUE(sol$converged)) "Success" else "Failure", nit = sol$iterations,
       nfev = sol$function_calls, root = root, converged = sol$converged, flag = sol$flag)
}
