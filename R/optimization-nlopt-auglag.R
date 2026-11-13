# The augmented Lagrangian method of NLopt 2.11.0 (AUGLAG, AUGLAG_EQ,
# LN_AUGLAG, LN_AUGLAG_EQ, LD_AUGLAG, LD_AUGLAG_EQ; src/algs/auglag/auglag.c,
# after Birgin and Martinez), ported to R. Each subproblem is solved by the
# local optimizer through the R port of nlopt_optimize.
#
# Copyright (c) 2007-2014 Massachusetts Institute of Technology, released
# under the MIT license (see R/optimization-nlopt-util.R).

# the AUGLAG cases of nlopt_optimize_
.nlopt_auglag_dispatch <- function(o, x, stop) {
  alg <- o$algorithm
  sub <- o$local_opt
  if (alg %in% c("AUGLAG", "AUGLAG_EQ") && is.null(sub)) {
    stop$msg <- "local optimizer must be specified for AUGLAG"
    return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  }
  if (is.null(sub)) { # default
    sub <- .nlopt_create(if (alg %in% c("LN_AUGLAG", "LN_AUGLAG_EQ")) "LN_COBYLA" else "LD_MMA", length(x))
    sub$ftol_rel <- o$ftol_rel; sub$ftol_abs <- o$ftol_abs
    sub$xtol_rel <- o$xtol_rel; sub$xtol_abs <- o$xtol_abs
    sub$maxeval <- -1L # nlopt_local_search_maxeval
  }
  if (!is.null(o$dx)) sub$dx <- o$dx
  sub$force_stop <- o$force_stop
  sub_has_fc <- alg %in% c("AUGLAG_EQ", "LN_AUGLAG_EQ", "LD_AUGLAG_EQ")
  # nlopt_add_inequality_constraint refuses a subsidiary algorithm without
  # inequality constraints (inequality_ok)
  if (sub_has_fc && length(o$fc) && !(sub$algorithm %in% c("LD_MMA", "LD_CCSAQ", "LN_COBYLA", "GN_ORIG_DIRECT",
      "GN_ORIG_DIRECT_L", "GN_ISRES", "AUGLAG", "AUGLAG_EQ", "LN_AUGLAG", "LN_AUGLAG_EQ", "LD_AUGLAG",
      "LD_AUGLAG_EQ", "LD_SLSQP", "GN_AGS")))
    return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  .nlopt_auglag_minimize(o, x, stop, sub, sub_has_fc)
}

.nlopt_auglag_minimize <- function(o, x, stop, sub, sub_has_fc) {
  n <- length(x)
  f <- o$f
  tau <- 0.5; gam <- 10
  lam_min <- -1e20; lam_max <- 1e20; mu_max <- 1e20
  fc <- if (sub_has_fc) list() else o$fc
  h <- o$h
  lambda <- numeric(.nlopt_count_constraints(h))
  mu <- numeric(.nlopt_count_constraints(fc))
  rho <- 1
  forced <- function() .nlopt_stop_forced(stop)
  # the augmented Lagrangian objective
  aug <- function(xx, grad) {
    r <- if (grad) o$fg(xx) else list(val = f(xx))
    L <- r$val; g <- if (grad) as.numeric(r$grad)
    stop$nevals <- stop$nevals + 1L
    out <- function() if (grad) list(val = L, grad = g) else L
    if (forced()) return(out())
    ii <- 0L
    for (cc in h) {
      cr <- if (grad) cc$fg(xx) else list(val = cc$f(xx))
      if (forced()) return(out())
      for (k in seq_len(cc$m)) {
        ii <- ii + 1L
        hk <- cr$val[k] + lambda[ii] / rho
        L <- L + 0.5 * rho * hk * hk
        if (grad) g <- g + (rho * hk) * cr$grad[k, ]
      }
    }
    ii <- 0L
    for (cc in fc) {
      cr <- if (grad) cc$fg(xx) else list(val = cc$f(xx))
      if (forced()) return(out())
      for (k in seq_len(cc$m)) {
        ii <- ii + 1L
        fk <- cr$val[k] + mu[ii] / rho
        if (fk > 0) {
          L <- L + 0.5 * rho * fk * fk
          if (grad) g <- g + (rho * fk) * cr$grad[k, ]
        }
      }
    }
    out()
  }
  sub$f <- function(xx) aug(xx, FALSE)
  sub$fg <- function(xx) aug(xx, TRUE)
  sub$lb <- o$lb; sub$ub <- o$ub
  sub$stopval <- if (!length(fc) && !length(h)) stop$minf_max else -Inf
  sub$fc <- if (sub_has_fc) o$fc else list()
  sub$h <- list()
  xcur <- x
  minf <- Inf; minf_penalty <- Inf; minf_feasible <- FALSE
  ICM <- Inf
  done <- function(ret) list(ret = ret, x = x, minf = minf)
  # the constraint values at xx: list(penalty, feasible, con2, h values, fc values)
  measure <- function(xx) {
    penalty <- 0; feasible <- TRUE; con2 <- 0; hv <- list(); fv <- list()
    for (i in seq_along(h)) {
      v <- h[[i]]$f(xx)
      if (forced()) return(NULL)
      for (k in seq_len(h[[i]]$m)) {
        penalty <- penalty + abs(v[k])
        feasible <- feasible && abs(v[k]) <= h[[i]]$tol[k]
        con2 <- con2 + v[k] * v[k]
      }
      hv[[i]] <- v
    }
    for (i in seq_along(fc)) {
      v <- fc[[i]]$f(xx)
      if (forced()) return(NULL)
      for (k in seq_len(fc[[i]]$m)) {
        penalty <- penalty + (if (v[k] > 0) v[k] else 0)
        feasible <- feasible && v[k] <= fc[[i]]$tol[k]
        if (v[k] > 0) con2 <- con2 + v[k] * v[k]
      }
      fv[[i]] <- v
    }
    list(penalty = penalty, feasible = feasible, con2 = con2, hv = hv, fv = fv)
  }
  if (length(h) || length(fc)) {
    stop$nevals <- stop$nevals + 1L
    fcur <- f(xcur)
    if (forced()) return(done("FORCED_STOP"))
    ms <- measure(xcur)
    if (is.null(ms)) return(done("FORCED_STOP"))
    minf <- fcur
    minf_penalty <- ms$penalty
    minf_feasible <- ms$feasible
    rho <- if (ms$con2 > 0) max(1e-6, min(10, 2 * abs(minf) / ms$con2)) else 10
  }
  repeat {
    prev_ICM <- ICM
    rs <- .nlopt_optimize_limited(sub, xcur, stop$maxeval - stop$nevals, stop$maxtime - (.nlopt_seconds() - stop$start))
    ret <- rs$ret
    if (ret %in% c("FAILURE", "INVALID_ARGS", "OUT_OF_MEMORY", "ROUNDOFF_LIMITED", "FORCED_STOP")) break
    xcur <- rs$x
    stop$nevals <- stop$nevals + 1L
    fcur <- f(xcur)
    if (forced()) return(done("FORCED_STOP"))
    ms <- measure(xcur)
    if (is.null(ms)) return(done("FORCED_STOP"))
    ICM <- 0
    ii <- 0L
    for (i in seq_along(h)) for (k in seq_len(h[[i]]$m)) {
      ii <- ii + 1L
      hi <- ms$hv[[i]][k]
      newlam <- lambda[ii] + rho * hi
      ICM <- max(ICM, abs(hi))
      lambda[ii] <- min(max(lam_min, newlam), lam_max)
    }
    ii <- 0L
    for (i in seq_along(fc)) for (k in seq_len(fc[[i]]$m)) {
      ii <- ii + 1L
      fci <- ms$fv[[i]][k]
      newmu <- mu[ii] + rho * fci
      ICM <- max(ICM, abs(max(fci, -mu[ii] / rho)))
      mu[ii] <- min(max(0.0, newmu), mu_max)
    }
    if (ICM > tau * prev_ICM) rho <- rho * gam
    feasible <- ms$feasible; penalty <- ms$penalty
    if ((feasible && (!minf_feasible || penalty < minf_penalty || fcur < minf)) ||
        (!minf_feasible && penalty < minf_penalty)) {
      ret <- "SUCCESS"
      if (feasible) {
        if (fcur < stop$minf_max) ret <- "STOPVAL_REACHED"
        else if (.nlopt_stop_ftol(stop, fcur, minf)) ret <- "FTOL_REACHED"
        else if (.nlopt_stop_x(stop, xcur, x)) ret <- "XTOL_REACHED"
      }
      minf <- fcur
      minf_penalty <- penalty
      minf_feasible <- feasible
      x <- xcur
      if (ret != "SUCCESS") break
    }
    if (forced()) { ret <- "FORCED_STOP"; break }
    if (.nlopt_stop_evals(stop)) { ret <- "MAXEVAL_REACHED"; break }
    if (.nlopt_stop_time(stop)) { ret <- "MAXTIME_REACHED"; break }
    if (ICM == 0) { ret <- "FTOL_REACHED"; break }
  }
  done(ret)
}
