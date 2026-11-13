# The conservative convex separable approximation algorithms of NLopt
# 2.11.0, ported to R: Svanberg's method of moving asymptotes (LD_MMA,
# src/algs/mma/mma.c) and the CCSA variant with quadratic penalties
# (LD_CCSAQ, ccsa_quadratic.c, without preconditioners, which NLopt.jl does
# not set). The dual problem of each inner iteration is solved by a nested
# NLopt optimizer (LD_MMA by default) over the constraint multipliers.
#
# Copyright (c) 2007-2014 Massachusetts Institute of Technology, released
# under the MIT license (see R/optimization-nlopt-util.R).

# nlopt_optimize_limited: nlopt_optimize with maxeval and maxtime replaced
# by the given limits when the optimizer has none or they are more stringent
.nlopt_optimize_limited <- function(o, x, maxeval, maxtime) {
  if (o$maxeval <= 0 || (maxeval > 0 && maxeval < o$maxeval)) o$maxeval <- as.integer(maxeval)
  if (o$maxtime <= 0 || (maxtime > 0 && maxtime < o$maxtime)) o$maxtime <- maxtime
  .nlopt_optimize(o, x)
}

# a bare optimizer description (nlopt_create)
.nlopt_create <- function(algorithm, n) {
  list(algorithm = algorithm, n = n, lb = rep(-Inf, n), ub = rep(Inf, n), stopval = -Inf, ftol_rel = 0,
       ftol_abs = 0, xtol_rel = 0, xtol_abs = rep(0, n), x_weights = NULL, maxeval = 0L, maxtime = 0,
       dx = NULL, local_opt = NULL, params = list(), population = 0L, vector_storage = 0L,
       fc = list(), h = list())
}

# the dual optimizer of the LD_MMA and LD_CCSAQ cases of nlopt_optimize_
.nlopt_ccsa_dual_opt <- function(o, m) {
  lo <- o$local_opt
  LO <- function(field, def) if (!is.null(lo)) lo[[field]] else def
  da <- o$params[["dual_algorithm"]]
  alg <- if (is.null(da)) LO("algorithm", "LD_MMA") else if (is.character(da)) da else .nlopt_algorithms[as.integer(da) + 1L]
  d <- .nlopt_create(alg, m)
  d$ftol_rel <- .nlopt_get_param(o, "dual_ftol_rel", LO("ftol_rel", 1e-14))
  d$ftol_abs <- .nlopt_get_param(o, "dual_ftol_abs", LO("ftol_abs", 0))
  d$xtol_rel <- .nlopt_get_param(o, "dual_xtol_rel", 0)
  d$xtol_abs <- rep(.nlopt_get_param(o, "dual_xtol_abs", 0), m)
  d$maxeval <- as.integer(.nlopt_get_param(o, "dual_maxeval", LO("maxeval", 100000L)))
  d
}

# mma_minimize (quadratic = FALSE) and ccsa_quadratic_minimize (TRUE).
# `fg(x)` gives list(val, grad); the constraints are a list of list(m, f,
# fg, tol) with fg giving list(val, grad = m by n matrix).
.nlopt_ccsa_minimize <- function(o, x, stop, quadratic) {
  n <- length(x); lb <- o$lb; ub <- o$ub
  fc <- o$fc
  inner_maxeval <- as.integer(.nlopt_get_param(o, "inner_maxeval", 0))
  rho_init <- .nlopt_get_param(o, "rho_init", 1.0)
  inner_gradients <- as.integer(.nlopt_get_param(o, "inner_gradients", 1))
  always_improve <- as.integer(.nlopt_get_param(o, "always_improve", 1))
  sigma_min <- .nlopt_get_param(o, "sigma_min", 0.0)
  bad <- function() list(ret = "INVALID_ARGS", x = x, minf = Inf)
  if (!(rho_init > 0) && !.nlopt_isinf(rho_init)) return(bad())
  if (!(inner_gradients %in% 0:1) || !(always_improve %in% 0:1) || sigma_min < 0) return(bad())
  m <- .nlopt_count_constraints(fc)
  dual <- .nlopt_ccsa_dual_opt(o, m)
  nan <- function(v) is.nan(v) | is.na(v)
  if (quadratic) nan <- function(v) rep(FALSE, length(v))
  tols <- unlist(lapply(fc, function(cc) cc$tol))
  # all constraints at x: list(val, grad) with gradients when wanted
  cons_at <- function(xx, grad) {
    val <- numeric(m); J <- matrix(0, m, n); i <- 0L
    for (cc in fc) {
      idx <- i + seq_len(cc$m)
      if (grad) { r <- cc$fg(xx); val[idx] <- r$val; J[idx, ] <- r$grad } else val[idx] <- cc$f(xx)
      i <- i + cc$m
      if (.nlopt_stop_forced(stop)) return(NULL)
    }
    list(val = val, grad = J)
  }
  sigma <- vapply(seq_len(n), function(j) {
    s <- if (!is.null(o$dx) && o$dx[j] > 0) o$dx[j] else if (.nlopt_isinf(ub[j]) || .nlopt_isinf(lb[j])) 1.0
      else 0.5 * (ub[j] - lb[j])
    max(s, sigma_min)
  }, 0)
  rho <- rho_init
  rhoc <- rep(rho_init, m)
  y <- numeric(m)
  dual_ub <- rep(Inf, m)
  # the dual data (dd)
  dd <- new.env(parent = emptyenv())
  # dual_func: minus the dual function at y, setting xcur, gval, wval, gcval
  dual_func <- function(yy) {
    dd$count <- dd$count + 1L
    val <- dd$gval <- dd$fval
    dd$wval <- 0
    gcval <- if (quadratic) fcval else ifelse(nan(fcval), 0, fcval)
    for (i in seq_len(m)) val <- val + yy[i] * gcval[i]
    xcur <- numeric(n)
    for (j in seq_len(n)) {
      if (sigma[j] == 0) { xcur[j] <- x[j]; next }
      if (quadratic) {
        u <- dd$rho; v <- dfdx[j]
        for (i in seq_len(m)) { u <- u + rhoc[i] * yy[i]; v <- v + dfcdx[i, j] * yy[i] }
        sigma2 <- sigma[j] * sigma[j]
        dx <- -sigma2 * v / u
        if (abs(dx) > sigma[j]) dx <- if (dx < 0 || (dx == 0 && 1 / dx < 0)) -sigma[j] else sigma[j]
        xcur[j] <- x[j] + dx
        if (xcur[j] > ub[j]) xcur[j] <- ub[j] else if (xcur[j] < lb[j]) xcur[j] <- lb[j]
        dx <- xcur[j] - x[j]
        dx2 <- dx * dx
        val <- val + v * dx + 0.5 * u * dx2 / sigma2
        dx2sig <- 0.5 * dx2 / sigma2
        dd$gval <- dd$gval + dfdx[j] * dx + dd$rho * dx2sig
        dd$wval <- dd$wval + dx2sig
        for (i in seq_len(m)) gcval[i] <- gcval[i] + dfcdx[i, j] * dx + rhoc[i] * dx2sig
      } else {
        u <- dfdx[j]
        v <- abs(dfdx[j]) * sigma[j] + 0.5 * dd$rho
        for (i in seq_len(m)) if (!nan(fcval[i])) {
          u <- u + dfcdx[i, j] * yy[i]
          v <- v + (abs(dfcdx[i, j]) * sigma[j] + 0.5 * rhoc[i]) * yy[i]
        }
        sigma2 <- sigma[j] * sigma[j]
        u <- u * sigma2
        q <- u / (v * sigma[j])
        dx <- (u / v) / (-1 - sqrt(abs(1 - q * q)))
        xcur[j] <- x[j] + dx
        if (xcur[j] > ub[j]) xcur[j] <- ub[j] else if (xcur[j] < lb[j]) xcur[j] <- lb[j]
        if (xcur[j] > x[j] + 0.9 * sigma[j]) xcur[j] <- x[j] + 0.9 * sigma[j]
        else if (xcur[j] < x[j] - 0.9 * sigma[j]) xcur[j] <- x[j] - 0.9 * sigma[j]
        dx <- xcur[j] - x[j]
        dx2 <- dx * dx
        denominv <- 1.0 / (sigma2 - dx2)
        val <- val + (u * dx + v * dx2) * denominv
        cc <- sigma2 * dx
        dd$gval <- dd$gval + (dfdx[j] * cc + (abs(dfdx[j]) * sigma[j] + 0.5 * dd$rho) * dx2) * denominv
        dd$wval <- dd$wval + 0.5 * dx2 * denominv
        for (i in seq_len(m)) if (!nan(fcval[i]))
          gcval[i] <- gcval[i] + (dfcdx[i, j] * cc + (abs(dfcdx[i, j]) * sigma[j] + 0.5 * rhoc[i]) * dx2) * denominv
      }
    }
    dd$xcur <- xcur
    dd$gcval <- gcval
    list(val = -val, grad = -gcval)
  }
  dual$f <- function(yy) dual_func(yy)$val
  dual$fg <- dual_func
  dual$lb <- rep(0, m)
  done <- function(ret) list(ret = ret, x = x, minf = minf)
  r0 <- o$fg(x)
  dd$fval <- fcur <- minf <- r0$val
  dfdx <- as.numeric(r0$grad)
  stop$nevals <- stop$nevals + 1L
  xcur <- x
  if (.nlopt_stop_forced(stop)) return(done("FORCED_STOP"))
  cs <- cons_at(x, TRUE)
  if (is.null(cs)) return(done("FORCED_STOP"))
  fcval <- cs$val; dfcdx <- cs$grad
  feasible <- TRUE; infeasibility <- 0
  for (i in seq_len(m)) {
    feasible <- feasible && (fcval[i] <= 0 || nan(fcval[i]))
    if (!nan(fcval[i]) && fcval[i] > infeasibility) infeasibility <- fcval[i]
  }
  if (!feasible) dual_ub <- rep(1e40, m)
  dual$ub <- dual_ub
  k <- 0L
  xprev <- xprevprev <- numeric(n)
  rhomin <- 1e-5
  repeat { # outer iterations
    inner_nevals <- 0L
    fprev <- fcur
    if (.nlopt_stop_forced(stop)) return(done("FORCED_STOP"))
    if (.nlopt_stop_evals(stop)) return(done("MAXEVAL_REACHED"))
    if (.nlopt_stop_time(stop)) return(done("MAXTIME_REACHED"))
    if (feasible && minf < stop$minf_max) return(done("STOPVAL_REACHED"))
    k <- k + 1L
    if (k > 1L) xprevprev <- xprev
    xprev <- xcur
    repeat { # inner iterations
      dd$rho <- rho; dd$count <- 0L
      rd <- .nlopt_optimize_limited(dual, y, 0, stop$maxtime - (.nlopt_seconds() - stop$start))
      if (rd$ret %in% c("FAILURE", "INVALID_ARGS", "OUT_OF_MEMORY", "ROUNDOFF_LIMITED", "FORCED_STOP",
                        "MAXTIME_REACHED")) return(done(rd$ret))
      y <- rd$x
      dual_func(y)
      xcur <- dd$xcur
      if (inner_gradients == 1L) { r <- o$fg(xcur); fcur <- r$val; dfdx_cur <- as.numeric(r$grad) }
      else fcur <- o$f(xcur)
      stop$nevals <- stop$nevals + 1L
      inner_nevals <- inner_nevals + 1L
      if (.nlopt_stop_forced(stop)) return(done("FORCED_STOP"))
      feasible_cur <- TRUE; infeasibility_cur <- 0
      new_infeasible_constraint <- FALSE
      inner_done <- dd$gval >= fcur
      cs <- cons_at(xcur, inner_gradients == 1L)
      if (is.null(cs)) return(done("FORCED_STOP"))
      fcval_cur <- cs$val
      if (inner_gradients == 1L) dfcdx_cur <- cs$grad
      for (i in seq_len(m)) if (!nan(fcval_cur[i])) {
        feasible_cur <- feasible_cur && (fcval_cur[i] <= tols[i])
        if (!nan(fcval[i])) inner_done <- inner_done && (dd$gcval[i] >= fcval_cur[i])
        else if (fcval_cur[i] > 0) new_infeasible_constraint <- TRUE
        if (fcval_cur[i] > infeasibility_cur) infeasibility_cur <- fcval_cur[i]
      }
      inner_done <- inner_done || (inner_maxeval > 0L && inner_nevals == inner_maxeval)
      accept <- if (always_improve == 1L) (fcur < minf && (inner_done || feasible_cur || !feasible)) ||
        (!feasible && infeasibility_cur < infeasibility) else inner_done
      if (accept) {
        if (inner_gradients == 0L) {
          r <- o$fg(xcur); fcur <- r$val; dfdx_cur <- as.numeric(r$grad)
          if (.nlopt_stop_forced(stop)) return(done("FORCED_STOP"))
          feasible_cur <- TRUE; infeasibility_cur <- 0
          new_infeasible_constraint <- FALSE
          if (!quadratic) inner_done <- dd$gval >= fcur
          cs <- cons_at(xcur, TRUE)
          if (is.null(cs)) return(done("FORCED_STOP"))
          fcval_cur <- cs$val; dfcdx_cur <- cs$grad
          for (i in seq_len(m)) if (!nan(fcval_cur[i])) {
            feasible_cur <- feasible_cur && (fcval_cur[i] <= tols[i])
            if (nan(fcval[i]) && fcval_cur[i] > 0) new_infeasible_constraint <- TRUE
            if (fcval_cur[i] > infeasibility_cur) infeasibility_cur <- fcval_cur[i]
          }
        }
        dd$fval <- minf <- fcur
        infeasibility <- infeasibility_cur
        fcval <- fcval_cur
        x <- xcur
        dfdx <- dfdx_cur
        if (m) dfcdx <- dfcdx_cur
        if (infeasibility_cur == 0) {
          if (!feasible) { dual_ub <- rep(Inf, m); dual$ub <- dual_ub }
          feasible <- TRUE
        } else if (new_infeasible_constraint) feasible <- FALSE
      }
      if (.nlopt_stop_forced(stop)) return(done("FORCED_STOP"))
      if (.nlopt_stop_evals(stop)) return(done("MAXEVAL_REACHED"))
      if (.nlopt_stop_time(stop)) return(done("MAXTIME_REACHED"))
      if (feasible && minf < stop$minf_max) return(done("STOPVAL_REACHED"))
      if (inner_done) break
      if (fcur > dd$gval) rho <- min(10 * rho, 1.1 * (rho + (fcur - dd$gval) / dd$wval))
      for (i in seq_len(m))
        if (!nan(fcval_cur[i]) && fcval_cur[i] > dd$gcval[i])
          rhoc[i] <- min(10 * rhoc[i], 1.1 * (rhoc[i] + (fcval_cur[i] - dd$gcval[i]) / dd$wval))
    }
    ret <- "SUCCESS"
    if (.nlopt_stop_ftol(stop, fcur, fprev)) ret <- "FTOL_REACHED"
    if (.nlopt_stop_x(stop, xcur, xprev)) ret <- "XTOL_REACHED"
    if (ret != "SUCCESS") return(done(ret))
    rho <- max(0.1 * rho, rhomin)
    rhoc <- pmax(0.1 * rhoc, rhomin)
    if (k > 1L) {
      for (j in seq_len(n)) {
        dx2 <- (xcur[j] - xprev[j]) * (xprev[j] - xprevprev[j])
        gam <- if (dx2 < 0) 0.7 else if (dx2 > 0) 1.2 else 1
        sigma[j] <- sigma[j] * gam
        if (!.nlopt_isinf(ub[j]) && !.nlopt_isinf(lb[j])) {
          sigma[j] <- min(sigma[j], 10 * (ub[j] - lb[j]))
          sigma[j] <- max(sigma[j], (if (quadratic) 1e-8 else 0.01) * (ub[j] - lb[j]))
        }
        sigma[j] <- max(sigma[j], sigma_min)
      }
    }
  }
}
