# The stochastic global algorithms of NLopt 2.11.0, ported to R: controlled
# random search with local mutation (GN_CRS2_LM, src/algs/crs/crs.c), the
# evolutionary strategy of Carlos Henrique da Silva Santos (GN_ESCH,
# esch/esch.c) and the improved stochastic ranking evolution strategy
# (GN_ISRES, isres/isres.c). Each draws from the Mersenne Twister of NLopt,
# so a seed set by NLopt$srand() reproduces the run of NLopt.
#
# Copyright (c) 2007-2014 Massachusetts Institute of Technology and
# Copyright (c) 2008-2013 Carlos Henrique da Silva Santos, released under
# the MIT license (see R/optimization-nlopt-util.R).

# ---------------------------------------------------------------------------
# CRS2 with local mutation (Kaelo and Ali 2006). The red black tree of crs.c
# sorts the population by f with ties broken by the row; order() is stable,
# so its first and last entries are the tree minimum and maximum.

.nlopt_crs_minimize <- function(f, lb, ub, x, stop, population) {
  n <- length(x)
  N <- if (population == 0L) 10L * (n + 1L) else as.integer(population)
  if (N < n + 1L) {
    stop$msg <- sprintf("population %d should be >= dimension + 1 = %d", N, n + 1L)
    return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  }
  ps <- matrix(0, N, n)
  fv <- rep(Inf, N)
  filled <- 0L
  best_row <- function() { o <- order(fv[seq_len(filled)]); o[1L] }
  worst_row <- function() { o <- order(fv[seq_len(filled)]); o[filled] }
  # random_trial: x = 2 G - x_n over the best point and n other random
  # points (Vitter's method A)
  random_trial <- function(b) {
    xt <- ps[b, ]
    i0 <- b - 1L
    jn <- .nlopt_iurand(n)
    add <- function(i) {
      xi <- ps[i + 1L, ]
      if (jn == 0L) xt <<- xt - xi * (0.5 * n) else xt <<- xt + xi
      jn <<- jn - 1L
    }
    Nleft <- N - 1L; nleft <- n
    Nfree <- Nleft - nleft
    i <- 0L; i <- i + (i == i0)
    while (nleft > 1L) {
      q <- Nfree / Nleft
      v <- .nlopt_urand(0, 1)
      while (q > v) {
        i <- i + 1L; i <- i + (i == i0)
        Nfree <- Nfree - 1L; Nleft <- Nleft - 1L
        q <- (q * Nfree) / Nleft
      }
      add(i)
      i <- i + 1L; i <- i + (i == i0)
      Nleft <- Nleft - 1L; nleft <- nleft - 1L
    }
    i <- i + .nlopt_iurand(Nleft); i <- i + (i == i0)
    add(i)
    xt <- xt * (2 / n)
    ifelse(xt > ub, ub, ifelse(xt < lb, lb, xt))
  }
  trial <- function() {
    b <- best_row(); w <- worst_row()
    mutation <- 1L
    p <- random_trial(b)
    repeat {
      pf <- f(p)
      stop$nevals <- stop$nevals + 1L
      if (.nlopt_stop_forced(stop)) return("FORCED_STOP")
      if (pf < fv[w]) break
      if (.nlopt_stop_evals(stop)) return("MAXEVAL_REACHED")
      if (.nlopt_stop_time(stop)) return("MAXTIME_REACHED")
      if (mutation > 0L) {
        for (i in seq_len(n)) {
          wt <- .nlopt_urand(0, 1)
          p[i] <- ps[b, i] * (1 + wt) - wt * p[i]
          if (p[i] > ub[i]) p[i] <- ub[i] else if (p[i] < lb[i]) p[i] <- lb[i]
        }
        mutation <- mutation - 1L
      } else {
        p <- random_trial(b)
        mutation <- 1L
      }
    }
    ps[w, ] <<- p
    fv[w] <<- pf
    "SUCCESS"
  }
  # crs_init
  init <- function() {
    for (k in seq_len(N)) {
      if (k == 1L) ps[1L, ] <<- x
      else for (j in seq_len(n)) ps[k, j] <<- .nlopt_urand(lb[j], ub[j])
      fv[k] <<- f(ps[k, ])
      stop$nevals <- stop$nevals + 1L
      filled <<- k
      if (fv[k] < stop$minf_max) return("STOPVAL_REACHED")
      if (.nlopt_stop_evals(stop)) return("MAXEVAL_REACHED")
      if (.nlopt_stop_time(stop)) return("MAXTIME_REACHED")
    }
    "SUCCESS"
  }
  ret <- init()
  b <- best_row()
  minf <- fv[b]; x <- ps[b, ]
  while (ret == "SUCCESS") {
    ret <- trial()
    if (ret == "SUCCESS") {
      b <- best_row()
      if (fv[b] < minf) {
        if (fv[b] < stop$minf_max) ret <- "STOPVAL_REACHED"
        else if (.nlopt_stop_f(stop, fv[b], minf)) ret <- "FTOL_REACHED"
        else if (.nlopt_stop_x(stop, ps[b, ], x)) ret <- "XTOL_REACHED"
        minf <- fv[b]; x <- ps[b, ]
      }
      if (ret != "SUCCESS") {
        if (.nlopt_stop_evals(stop)) ret <- "MAXEVAL_REACHED"
        else if (.nlopt_stop_time(stop)) ret <- "MAXTIME_REACHED"
      }
    }
  }
  list(ret = ret, x = x, minf = minf)
}

# ---------------------------------------------------------------------------
# ESCH: np parents, no offspring, one point crossover and Cauchy mutation

.nlopt_esch_randcauchy <- function(lo, hi) {
  # randcauchy with mi = 0, t = 1, band = 10
  repeat {
    u <- .nlopt_urand(0, 1)
    cm <- tan((u - 0.5) * 0x1.921fb54442d18p+1)
    if (!(cm < -5 || cm > 5)) break
  }
  cm <- if (cm < 0) -cm else cm + 5
  lo + (hi - lo) * (cm / 10)
}

.nlopt_esch_minimize <- function(f, lb, ub, x, stop, np, no) {
  n <- length(x)
  if (np == 0L) np <- 40L
  if (no == 0L) no <- 60L
  minf <- Inf
  cauchy_row <- function() vapply(seq_len(n), function(j) .nlopt_esch_randcauchy(lb[j], ub[j]), 0)
  par <- vector("list", np); off <- vector("list", no)
  for (id in seq_len(np)) par[[id]] <- cauchy_row()
  par[[1L]] <- x
  for (id in seq_len(no)) off[[id]] <- cauchy_row()
  pfit <- numeric(np); ofit <- numeric(no)
  eval_one <- function(xi) {
    fi <- f(xi)
    stop$nevals <- stop$nevals + 1L
    if (minf > fi) { minf <<- fi; x <<- xi }
    ret <- if (.nlopt_stop_forced(stop)) "FORCED_STOP"
      else if (minf < stop$minf_max) "STOPVAL_REACHED"
      else if (.nlopt_stop_evals(stop)) "MAXEVAL_REACHED"
      else if (.nlopt_stop_time(stop)) "MAXTIME_REACHED"
      else "SUCCESS"
    list(f = fi, ret = ret)
  }
  for (id in seq_len(np)) {
    r <- eval_one(par[[id]]); pfit[id] <- r$f
    if (r$ret != "SUCCESS") return(list(ret = r$ret, x = x, minf = minf))
  }
  totalmutation <- max((no * n) %/% 10L, 1L)
  repeat {
    for (id in seq_len(no)) { # crossover
      p1 <- .nlopt_iurand(np) + 1L
      p2 <- .nlopt_iurand(np) + 1L
      cp <- .nlopt_iurand(n)
      off[[id]] <- c(par[[p1]][seq_len(cp)], par[[p2]][seq_len(n - cp) + cp])
    }
    for (k in seq_len(totalmutation)) { # mutation
      io <- .nlopt_iurand(no) + 1L
      ip <- .nlopt_iurand(n) + 1L
      off[[io]][ip] <- .nlopt_esch_randcauchy(lb[ip], ub[ip])
    }
    for (id in seq_len(no)) {
      r <- eval_one(off[[id]]); ofit[id] <- r$f
      if (r$ret != "SUCCESS") return(list(ret = r$ret, x = x, minf = minf))
    }
    # selection: the best np of parents and offspring become the parents
    all <- c(par, off); fit <- c(pfit, ofit)
    o <- .nlopt_qsort_order(fit)
    all <- all[o]; fit <- fit[o]
    par <- all[seq_len(np)]; pfit <- fit[seq_len(np)]
    off <- all[np + seq_len(no)]; ofit <- fit[np + seq_len(no)]
  }
}

# ---------------------------------------------------------------------------
# ISRES (Runarsson and Yao 2005): stochastic ranking of objective and
# squared constraint violation, then differential variation of the best
# survivors and log normal mutation of the rest

.nlopt_isres_minimize <- function(f, fc, h, lb, ub, x, stop, population) {
  ALPHA <- 0.2; GAMMA <- 0.85; PHI <- 1.0; PF <- 0.45; SURVIVOR <- 1.0 / 7.0
  n <- length(x)
  minf <- Inf
  minf_penalty <- Inf; minf_gpenalty <- Inf
  if (population == 0L) population <- 20L * (n + 1L)
  if (population < 1L) {
    stop$msg <- sprintf("population %d is too small", population)
    return(list(ret = "INVALID_ARGS", x = x, minf = minf))
  }
  survivors <- as.integer(ceiling(population * SURVIVOR))
  taup <- PHI / sqrt(2 * n)
  tau <- PHI / sqrt(2 * sqrt(n))
  if (any(vapply(c(lb, ub), .nlopt_isinf, TRUE))) {
    stop$msg <- "isres requires a finite search region"
    return(list(ret = "INVALID_ARGS", x = x, minf = minf))
  }
  sigmax <- (ub - lb) / sqrt(n)
  sigmas <- matrix(0, population, n); xs <- matrix(0, population, n)
  for (k in seq_len(population)) for (j in seq_len(n)) {
    sigmas[k, j] <- (ub[j] - lb[j]) / sqrt(n)
    xs[k, j] <- .nlopt_urand(lb[j], ub[j])
  }
  xs[1L, ] <- x
  fval <- numeric(population); penalty <- numeric(population)
  done <- function(ret) list(ret = ret, x = x, minf = minf)
  repeat { # one generation
    all_feasible <- TRUE
    for (k in seq_len(population)) {
      feasible <- TRUE
      stop$nevals <- stop$nevals + 1L
      xk <- xs[k, ]
      fval[k] <- f(xk)
      if (.nlopt_stop_forced(stop)) return(done("FORCED_STOP"))
      penalty[k] <- 0
      for (cc in fc) {
        g <- cc$f(xk)
        if (.nlopt_stop_forced(stop)) return(done("FORCED_STOP"))
        for (r in seq_len(cc$m)) {
          gval <- g[r]
          if (gval > cc$tol[r]) feasible <- FALSE
          if (gval < 0) gval <- 0
          penalty[k] <- penalty[k] + gval * gval
        }
      }
      gpenalty <- penalty[k]
      for (cc in h) {
        hv <- cc$f(xk)
        if (.nlopt_stop_forced(stop)) return(done("FORCED_STOP"))
        for (r in seq_len(cc$m)) {
          hval <- hv[r]
          if (abs(hval) > cc$tol[r]) feasible <- FALSE
          penalty[k] <- penalty[k] + hval * hval
        }
      }
      if (penalty[k] > 0) all_feasible <- FALSE
      ret <- "SUCCESS"
      pk <- if (feasible) 0 else penalty[k]
      if ((penalty[k] <= minf_penalty || feasible) && (fval[k] <= minf || minf_gpenalty > 0) &&
          (pk != minf_penalty || fval[k] != minf)) {
        if (fval[k] < stop$minf_max && feasible) ret <- "STOPVAL_REACHED"
        else if (!.nlopt_isinf(minf)) {
          if (.nlopt_stop_f(stop, fval[k], minf) && .nlopt_stop_f(stop, pk, minf_penalty)) ret <- "FTOL_REACHED"
          else if (.nlopt_stop_x(stop, xk, x)) ret <- "XTOL_REACHED"
        }
        x <- xk
        minf <- fval[k]
        minf_penalty <- pk
        minf_gpenalty <- if (feasible) 0 else gpenalty
        if (ret != "SUCCESS") return(done(ret))
      }
      if (.nlopt_stop_forced(stop)) return(done("FORCED_STOP"))
      if (.nlopt_stop_evals(stop)) return(done("MAXEVAL_REACHED"))
      if (.nlopt_stop_time(stop)) return(done("MAXTIME_REACHED"))
    }
    # selection: rank the population (irank holds 1 based rows)
    irank <- seq_len(population)
    if (all_feasible) irank <- .nlopt_qsort_order(fval)
    else {
      for (i in seq_len(population)) {
        swapped <- FALSE
        for (j in seq_len(population - 1L)) {
          u <- .nlopt_urand(0, 1)
          a <- irank[j]; b <- irank[j + 1L]
          if (u < PF || (penalty[a] == 0 && penalty[b] == 0)) {
            if (fval[a] > fval[b]) { irank[j] <- b; irank[j + 1L] <- a; swapped <- TRUE }
          } else if (penalty[a] > penalty[b]) {
            irank[j] <- b; irank[j + 1L] <- a; swapped <- TRUE
          }
        }
        if (!swapped) break
      }
    }
    # standard mutation of the rest
    if (survivors < population) for (k in (survivors + 1L):population) {
      taup_rand <- taup * .nlopt_nrand(0, 1)
      rk <- irank[k]
      ri <- irank[(k - 1L) %% survivors + 1L]
      for (j in seq_len(n)) {
        sigmas[rk, j] <- sigmas[ri, j] * exp(taup_rand + tau * .nlopt_nrand(0, 1))
        if (sigmas[rk, j] > sigmax[j]) sigmas[rk, j] <- sigmax[j]
        repeat {
          xs[rk, j] <- xs[ri, j] + sigmas[rk, j] * .nlopt_nrand(0, 1)
          if (!(xs[rk, j] < lb[j] || xs[rk, j] > ub[j])) break
        }
        sigmas[rk, j] <- sigmas[ri, j] + ALPHA * (sigmas[rk, j] - sigmas[ri, j])
      }
    }
    # differential variation of the survivors; as in isres.c the base point
    # is row 1 and the difference uses row k + 1 of the unranked population
    x0 <- xs[1L, ]
    for (k in seq_len(survivors)) {
      taup_rand <- taup * .nlopt_nrand(0, 1)
      rk <- irank[k]
      for (j in seq_len(n)) {
        xi <- xs[rk, j]
        if (k < survivors) xs[rk, j] <- xs[rk, j] + GAMMA * (x0[j] - xs[k + 1L, j])
        if (k == survivors || xs[rk, j] < lb[j] || xs[rk, j] > ub[j]) {
          sigi <- sigmas[rk, j]
          sigmas[rk, j] <- sigmas[rk, j] * exp(taup_rand + tau * .nlopt_nrand(0, 1))
          if (sigmas[rk, j] > sigmax[j]) sigmas[rk, j] <- sigmax[j]
          repeat {
            xs[rk, j] <- xi + sigmas[rk, j] * .nlopt_nrand(0, 1)
            if (!(xs[rk, j] < lb[j] || xs[rk, j] > ub[j])) break
          }
          sigmas[rk, j] <- sigi + ALPHA * (sigmas[rk, j] - sigi)
        }
      }
    }
  }
}
