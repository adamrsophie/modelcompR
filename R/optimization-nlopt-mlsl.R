# Multi level single linkage of NLopt 2.11.0 (G_MLSL, G_MLSL_LDS, GN_MLSL,
# GD_MLSL, GN_MLSL_LDS, GD_MLSL_LDS; src/algs/mlsl/mlsl.c, after Rinnooy Kan
# and Timmer and Kucherenko and Sytsko), ported to R. The sample points
# come from NLopt's Mersenne Twister or, for the LDS variants, from the
# Sobol sequence of sobolseq.c (.sobol_new, the same Joe and Kuo
# construction).
#
# Copyright (c) 2007-2014 Massachusetts Institute of Technology, released
# under the MIT license (see R/optimization-nlopt-util.R).
#
# The red black trees of mlsl.c compare by f only and insert an equal key
# to the left, so among equal values the later point comes first: the
# trees are kept as vectors ordered by (f, minus the insertion number).

# the MLSL cases of nlopt_optimize_
.nlopt_mlsl_dispatch <- function(o, x, stop, pop) {
  alg <- o$algorithm
  n <- length(x)
  if (!.nlopt_finite_domain(o$lb, o$ub)) return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  local <- o$local_opt
  if (is.null(local) && alg %in% c("G_MLSL", "G_MLSL_LDS")) return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  if (is.null(local)) {
    local <- .nlopt_create(if (alg %in% c("GN_MLSL", "GN_MLSL_LDS")) "LN_COBYLA" else "LD_MMA", n)
    local$ftol_rel <- o$ftol_rel; local$ftol_abs <- o$ftol_abs
    local$xtol_rel <- o$xtol_rel; local$xtol_abs <- o$xtol_abs
    local$maxeval <- -1L
  }
  if (!is.null(o$dx)) local$dx <- o$dx
  xa <- stop$xtol_abs
  if (local$ftol_rel <= 0 && local$ftol_abs <= 0 && local$xtol_rel <= 0 && !(length(xa) && all(xa > 0))) {
    # not sensible to run MLSL without some tolerance for the local search
    local$ftol_rel <- 1e-15
    local$xtol_rel <- 1e-7
  }
  local$force_stop <- o$force_stop
  .nlopt_mlsl_minimize(o, x, stop, local, pop, alg %in% c("GN_MLSL_LDS", "GD_MLSL_LDS", "G_MLSL_LDS"))
}

.nlopt_mlsl_minimize <- function(o, x, stop, local, Nsamples, lds) {
  n <- length(x); lb <- o$lb; ub <- o$ub; f <- o$f
  N <- if (Nsamples == 0L) 4L else Nsamples
  if (N < 1L) return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  K2PI <- 0x1.921fb54442d18p+2
  s <- if (lds) .sobol_new(n)
  # fcount: the local searches count against the MLSL evaluations
  local$f <- function(xx) { stop$nevals <- stop$nevals + 1L; f(xx) }
  if (is.function(o$fg)) local$fg <- function(xx) { stop$nevals <- stop$nevals + 1L; o$fg(xx) }
  local$lb <- lb; local$ub <- ub
  local$stopval <- stop$minf_max
  local$fc <- list(); local$h <- list()
  gamma <- 0.3
  z <- n %/% 2 # integer division, as in gam()
  gam <- sqrt((K2PI * z)^(1.0 / n) * z) * exp(-0.5)
  R_prefactor <- sqrt(2. / K2PI) * (gam * 2.)^(1.0 / n)
  for (i in seq_len(n)) R_prefactor <- R_prefactor * (ub[i] - lb[i])^(1.0 / n)
  dlm <- 1.0; dbound <- 1e-6
  # the points: x (rows), f, minimized, closest distances, insertion number
  px <- matrix(0, 0, n); pf <- numeric(); pmin_ <- logical(); cpd <- numeric(); cld <- numeric(); pid <- integer()
  lmx <- matrix(0, 0, n); lmf <- numeric(); lmid <- integer()
  ins <- 0L
  d2 <- function(a, b) { d <- 0; for (i in seq_len(n)) { dx <- a[i] - b[i]; d <- d + dx * dx }; d }
  pt_order <- function() order(pf, -pid)
  if (lds) { k <- 1; while (k * 2 < 10 * n + N) k <- k * 2; for (q in seq_len(k)) .sobol_next(s) }
  add_pt <- function(xx, val) {
    ins <<- ins + 1L
    px <<- rbind(px, xx); pf <<- c(pf, val); pmin_ <<- c(pmin_, FALSE)
    cpd <<- c(cpd, Inf); cld <<- c(cld, Inf); pid <<- c(pid, ins)
    length(pf)
  }
  minf <- Inf
  get_minf <- function() {
    if (length(pf)) { b <- pt_order()[1L]; minf <<- pf[b]; x <<- px[b, ] }
    if (length(lmf)) {
      b <- order(lmf, -lmid)[1L]
      if (lmf[b] < minf) { minf <<- lmf[b]; x <<- lmx[b, ] }
    }
  }
  status <- function(val) {
    if (.nlopt_stop_forced(stop)) "FORCED_STOP" else if (.nlopt_stop_evals(stop)) "MAXEVAL_REACHED"
    else if (.nlopt_stop_time(stop)) "MAXTIME_REACHED" else if (val < stop$minf_max) "STOPVAL_REACHED" else "SUCCESS"
  }
  v <- f(x)
  stop$nevals <- stop$nevals + 1L
  add_pt(x, v)
  ret <- status(v)
  while (ret == "SUCCESS") {
    get_minf()
    i <- 0L
    while (i < N && ret == "SUCCESS") { # sampling
      i <- i + 1L
      xx <- if (lds) lb + (ub - lb) * .sobol_next(s) else vapply(seq_len(n), function(j) .nlopt_urand(lb[j], ub[j]), 0)
      val <- f(xx)
      stop$nevals <- stop$nevals + 1L
      k <- add_pt(xx, val)
      ret <- status(val)
      if (ret == "SUCCESS") {
        # find_closest_pt, find_closest_lm, pts_update_newpt
        for (q in which(pf < val)) cpd[k] <- min(cpd[k], d2(xx, px[q, ]))
        for (q in which(lmf < val)) cld[k] <- min(cld[k], d2(xx, lmx[q, ]))
        for (q in which(pf > val & !pmin_)) { dd <- d2(xx, px[q, ]); if (dd < cpd[q]) cpd[q] <- dd }
      }
    }
    Np <- length(pf)
    R <- R_prefactor * (log(Np) / Np)^(1.0 / n)
    ord <- pt_order()
    cnt <- as.integer(ceiling(gamma * Np) + 0.5)
    for (q in ord) {
      if (cnt <= 0L || ret != "SUCCESS") break
      cnt <- cnt - 1L
      xq <- px[q, ]
      potential <- !pmin_[q] && cpd[q] > R * R && cld[q] > (dlm * R) * (dlm * R) &&
        !any((xq - lb <= dbound * R | ub - xq <= dbound * R) & ub - lb > dbound * R)
      if (!potential) next
      t <- .nlopt_seconds()
      if (.nlopt_stop_forced(stop)) { ret <- "FORCED_STOP"; break }
      if (.nlopt_stop_evals(stop)) { ret <- "MAXEVAL_REACHED"; break }
      if (stop$maxtime > 0 && t - stop$start >= stop$maxtime) { ret <- "MAXTIME_REACHED"; break }
      lr <- .nlopt_optimize_limited(local, xq, stop$maxeval - stop$nevals, stop$maxtime - (t - stop$start))
      pmin_[q] <- TRUE
      if (lr$ret %in% c("FAILURE", "INVALID_ARGS", "OUT_OF_MEMORY", "ROUNDOFF_LIMITED", "FORCED_STOP"))
        return(list(ret = lr$ret, x = x, minf = minf))
      ins <- ins + 1L
      lmx <- rbind(lmx, lr$x); lmf <- c(lmf, lr$minf); lmid <- c(lmid, ins)
      if (.nlopt_stop_forced(stop)) ret <- "FORCED_STOP"
      else if (lr$minf < stop$minf_max) ret <- "STOPVAL_REACHED"
      else if (.nlopt_stop_evals(stop)) ret <- "MAXEVAL_REACHED"
      else if (.nlopt_stop_time(stop)) ret <- "MAXTIME_REACHED"
      else for (r in which(pf > lr$minf & !pmin_)) { dd <- d2(lr$x, px[r, ]); if (dd < cld[r]) cld[r] <- dd }
    }
  }
  get_minf()
  list(ret = ret, x = x, minf = minf)
}
