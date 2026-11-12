# The stochastic global optimizers of scipy.optimize (SciPy 1.18.1):
# differential_evolution (_differentialevolution.py), basinhopping
# (_basinhopping.py) and dual_annealing (_dual_annealing.py). Random numbers
# come from R's generator, so runs are reproducible with set.seed() but do
# not repeat the NumPy streams.

.sp_fmod <- function(a, b) a - b * trunc(a / b)

# --------------------------------------------------------------------------
# differential_evolution

.sp_de_strategies <- c(best1bin = "best1", randtobest1bin = "randtobest1", currenttobest1bin = "currenttobest1",
                       best2bin = "best2", rand2bin = "rand2", rand1bin = "rand1", best1exp = "best1",
                       rand1exp = "rand1", randtobest1exp = "randtobest1", currenttobest1exp = "currenttobest1",
                       best2exp = "best2", rand2exp = "rand2")

.sp_differential_evolution <- function(func, lb, ub, strategy = "best1bin", maxiter = 1000L, popsize = 15L,
                                       tol = 0.01, mutation = c(0.5, 1), recombination = 0.7, maxfun = Inf,
                                       callback = NULL, disp = FALSE, polish = TRUE, init = "latinhypercube",
                                       atol = 0, updating = "immediate", workers = 1L, x0 = NULL,
                                       integrality = NULL, vectorized = FALSE, ...) {
  if (!strategy %in% names(.sp_de_strategies)) stop("Please select a valid mutation strategy", call. = FALSE)
  if (any(!is.finite(mutation)) || any(mutation >= 2) || any(mutation < 0))
    stop("The mutation constant must be a float in U[0, 2), or specified as a tuple(min, max) where min < max ",
         "and min, max are in U[0, 2).", call. = FALSE)
  if (!all(is.finite(c(lb, ub))))
    stop("bounds should be a sequence containing finite real valued (min, max) pairs for each value in x",
         call. = FALSE)
  if (is.null(maxiter)) maxiter <- 1000L
  if (is.null(maxfun)) maxfun <- Inf
  dither <- if (length(mutation) > 1L) sort(mutation[1:2])
  scale <- mutation[1L]
  N <- length(lb)
  arg1 <- 0.5 * (lb + ub)
  arg2 <- abs(lb - ub)
  rarg2 <- 1 / arg2
  rarg2[!is.finite(rarg2)] <- 0
  integ <- if (!is.null(integrality) && any(integrality)) rep_len(as.logical(integrality), N) else rep(FALSE, N)
  if (any(integ)) {
    l2 <- ceiling(lb); u2 <- floor(ub)
    if (!all(l2[integ] <= u2[integ]))
      stop("One of the integrality constraints does not have any possible integer values between the lower/upper ",
           "bounds.", call. = FALSE)
    lb[integ] <- l2[integ] - 0.5 + 2^-52 * abs(l2[integ] - 0.5)
    ub[integ] <- u2[integ] + 0.5 - 2^-52 * abs(u2[integ] + 0.5)
    arg1 <- 0.5 * (lb + ub); arg2 <- abs(lb - ub); rarg2 <- 1 / arg2; rarg2[!is.finite(rarg2)] <- 0
  }
  scale_p <- function(trial) {
    s <- arg1 + (trial - 0.5) * arg2
    if (any(integ)) s[integ] <- round(s[integ])
    s
  }
  unscale_p <- function(par) (par - arg1) * rarg2 + 0.5
  S <- max(5L, popsize * max(1L, N - sum(lb == ub)))
  if (is.character(init) && init == "latinhypercube") {
    segsize <- 1.0 / S
    samples <- segsize * matrix(stats::runif(S * N), S, N, byrow = TRUE) + (0:(S - 1L)) / S
    pop <- matrix(0, S, N)
    for (j in seq_len(N)) pop[, j] <- samples[sample.int(S), j]
  } else if (is.character(init) && init == "random") {
    pop <- matrix(stats::runif(S * N), S, N, byrow = TRUE)
  } else if (is.matrix(init)) {
    if (nrow(init) < 5L || ncol(init) != N)
      stop("The population supplied needs to have shape (S, len(x)), where S > 4.", call. = FALSE)
    pop <- pmin(pmax(t(apply(init, 1L, unscale_p)), 0), 1)
    if (N == 1L) pop <- matrix(pop, ncol = 1L)
    S <- nrow(pop)
  } else {
    stop("The population initialization method must be one of 'latinhypercube' or 'random', or an array of shape ",
         "(S, N) where N is the number of parameters and S>5", call. = FALSE)
  }
  if (!is.null(x0)) {
    xs <- unscale_p(as.numeric(x0))
    if (any(xs > 1 | xs < 0)) stop("Some entries in x0 lay outside the specified bounds", call. = FALSE)
    pop[1L, ] <- xs
  }
  energies <- rep(Inf, S)
  nfev <- 0
  fn <- function(x) {
    nfev <<- nfev + 1
    as.numeric(func(x))
  }
  promote <- function() {
    l <- which.min(energies)
    if (length(l) && l != 1L) {
      energies[c(1L, l)] <<- energies[c(l, 1L)]
      pop[c(1L, l), ] <<- pop[c(l, 1L), , drop = FALSE]
    }
  }
  select <- function(candidate, k) {
    idx <- sample.int(S)[seq_len(k + 1L)]
    idx[idx != candidate][seq_len(k)]
  }
  kind <- .sp_de_strategies[[strategy]]
  binomial <- endsWith(strategy, "bin")
  mutate <- function(candidate) {
    fill_point <- sample.int(N, 1L)
    r <- select(candidate, 5L)
    trial <- pop[candidate, ]
    bprime <- switch(kind,
      best1 = pop[1L, ] + scale * (pop[r[1L], ] - pop[r[2L], ]),
      rand1 = pop[r[1L], ] + scale * (pop[r[2L], ] - pop[r[3L], ]),
      randtobest1 = {
        b <- pop[r[1L], ]
        b <- b + scale * (pop[1L, ] - b)
        b + scale * (pop[r[2L], ] - pop[r[3L], ])
      },
      currenttobest1 = pop[candidate, ] + scale * (pop[1L, ] - pop[candidate, ] + pop[r[1L], ] - pop[r[2L], ]),
      best2 = pop[1L, ] + scale * (pop[r[1L], ] + pop[r[2L], ] - pop[r[3L], ] - pop[r[4L], ]),
      rand2 = pop[r[1L], ] + scale * (pop[r[2L], ] + pop[r[3L], ] - pop[r[4L], ] - pop[r[5L], ]))
    crossovers <- stats::runif(N) < recombination
    if (binomial) {
      crossovers[fill_point] <- TRUE
      trial <- ifelse(crossovers, bprime, trial)
    } else {
      i <- 0L
      crossovers[1L] <- TRUE
      while (i < N && crossovers[i + 1L]) {
        trial[fill_point] <- bprime[fill_point]
        fill_point <- fill_point %% N + 1L
        i <- i + 1L
      }
    }
    oob <- trial > 1 | trial < 0
    if (any(oob)) trial[oob] <- stats::runif(sum(oob))
    trial
  }
  converged <- function() {
    if (any(is.infinite(energies))) return(FALSE)
    sd0 <- sqrt(mean((energies - mean(energies))^2))
    sd0 <= atol + tol * abs(mean(energies))
  }
  # the initial population
  Sv <- min(S, maxfun - nfev)
  for (i in seq_len(Sv)) energies[i] <- fn(scale_p(pop[i, ]))
  promote()
  warning_flag <- FALSE
  message <- .sp_status_message$success
  nit <- 0L
  stopped <- FALSE
  for (nit in seq_len(maxiter)) {
    if (!is.null(dither)) scale <- stats::runif(1L, dither[1L], dither[2L])
    if (updating == "immediate") {
      for (candidate in seq_len(S)) {
        if (nfev > maxfun) {
          stopped <- TRUE
          break
        }
        trial <- mutate(candidate)
        energy <- fn(scale_p(trial))
        if (energy <= energies[candidate]) {
          pop[candidate, ] <- trial
          energies[candidate] <- energy
          if (energy <= energies[1L]) promote()
        }
      }
    } else {
      if (nfev >= maxfun) {
        stopped <- TRUE
      } else {
        trials <- t(vapply(seq_len(S), mutate, numeric(N)))
        if (N == 1L) trials <- matrix(trials, ncol = 1L)
        k <- min(S, maxfun - nfev)
        te <- rep(Inf, S)
        for (i in seq_len(k)) te[i] <- fn(scale_p(trials[i, ]))
        loc <- te <= energies
        pop[loc, ] <- trials[loc, , drop = FALSE]
        energies[loc] <- te[loc]
        promote()
      }
    }
    if (stopped) {
      warning_flag <- TRUE
      message <- if (nfev > maxfun) .sp_status_message$maxfev else
        "Maximum number of function evaluations has been reached."
      break
    }
    if (is.function(callback)) {
      if (isTRUE(callback(list(x = scale_p(pop[1L, ]), fun = energies[1L], nit = nit)))) {
        warning_flag <- TRUE
        message <- "callback function requested stop early"
      }
    }
    if (warning_flag || converged()) break
    if (nit == maxiter) {
      message <- .sp_status_message$maxiter
      warning_flag <- TRUE
    }
  }
  x <- scale_p(pop[1L, ])
  fun <- energies[1L]
  jac <- NULL
  if (!isFALSE(polish) && !all(integ)) {
    plb <- lb; pub <- ub
    if (any(integ)) { plb[integ] <- x[integ]; pub[integ] <- x[integ] }
    res <- .sp_minimize(fn, x, "L-BFGS-B", bounds = list(lb = plb, ub = pub))
    if (res$fun < fun && isTRUE(res$success) && all(res$x <= pub) && all(plb <= res$x)) {
      fun <- res$fun
      x <- res$x
      jac <- res$jac
    }
  }
  list(x = x, fun = fun, nfev = as.integer(nfev), nit = as.integer(nit), message = message,
       success = !warning_flag, jac = jac, population = t(apply(pop, 1L, scale_p)), population_energies = energies)
}

# --------------------------------------------------------------------------
# basinhopping

.sp_basinhopping <- function(func, x0, niter = 100L, T = 1.0, stepsize = 0.5, minimizer_kwargs = NULL,
                             take_step = NULL, accept_test = NULL, callback = NULL, interval = 50L, disp = FALSE,
                             niter_success = NULL, target_accept_rate = 0.5, stepwise_factor = 0.9, ...) {
  if (target_accept_rate <= 0 || target_accept_rate >= 1)
    stop("target_accept_rate has to be in range (0, 1)", call. = FALSE)
  if (stepwise_factor <= 0 || stepwise_factor >= 1) stop("stepwise_factor has to be in range (0, 1)", call. = FALSE)
  mk <- minimizer_kwargs %||% list()
  minimizer <- function(x) {
    method <- mk$method %||% (if (length(mk$constraints)) "SLSQP" else if (!is.null(mk$bounds)) "L-BFGS-B" else "BFGS")
    .sp_minimize(func, x, method, jac = mk$jac, hess = mk$hess, bounds = mk$bounds,
                 constraints = mk$constraints %||% list(), options = mk$options %||% list())
  }
  st <- new.env(parent = emptyenv())
  st$stepsize <- stepsize; st$nstep <- 0L; st$naccept <- 0L
  step <- function(x) {
    st$nstep <- st$nstep + 1L
    if (st$nstep %% interval == 0L) {
      accept_rate <- st$naccept / st$nstep
      st$stepsize <- if (accept_rate > target_accept_rate) st$stepsize / stepwise_factor else
        st$stepsize * stepwise_factor
    }
    if (is.function(take_step)) take_step(x) else x + stats::runif(length(x), -st$stepsize, st$stepsize)
  }
  beta <- if (T != 0) 1.0 / T else Inf
  metropolis <- function(new, old) {
    prod <- -(new$fun - old$fun) * beta
    w <- exp(min(0, prod))
    rand <- stats::runif(1L)
    isTRUE(w >= rand) && (isTRUE(new$success) || !isTRUE(old$success))
  }
  if (is.null(niter_success)) niter_success <- niter + 2L
  minres <- minimizer(as.numeric(x0))
  failures <- if (isTRUE(minres$success)) 0L else 1L
  nfev <- minres$nfev %||% 0L; njev <- minres$njev %||% 0L; nhev <- minres$nhev %||% 0L
  x <- minres$x; energy <- minres$fun; incumbent <- minres; lowest <- minres
  if (is.function(callback)) callback(lowest$x, lowest$fun, TRUE)
  count <- 0L
  message <- "requested number of basinhopping iterations completed successfully"
  i <- 0L
  for (i in seq_len(niter)) {
    res <- minimizer(step(x))
    if (!isTRUE(res$success)) failures <- failures + 1L
    nfev <- nfev + (res$nfev %||% 0L); njev <- njev + (res$njev %||% 0L); nhev <- nhev + (res$nhev %||% 0L)
    accept <- TRUE
    if (is.function(accept_test)) {
      tr <- accept_test(f_new = res$fun, x_new = res$x, f_old = energy, x_old = x)
      if (identical(tr, "force accept")) accept <- TRUE else if (!isTRUE(tr)) accept <- FALSE
      if (!identical(tr, "force accept") && !metropolis(res, incumbent)) accept <- FALSE
    } else if (!metropolis(res, incumbent)) {
      accept <- FALSE
    }
    if (accept) st$naccept <- st$naccept + 1L
    new_global_min <- FALSE
    if (accept) {
      energy <- res$fun
      x <- res$x
      incumbent <- res
      if (isTRUE(res$success) && (res$fun < lowest$fun || !isTRUE(lowest$success))) {
        lowest <- res
        new_global_min <- TRUE
      }
    }
    if (is.function(callback) && isTRUE(callback(res$x, res$fun, accept))) {
      message <- "callback function requested stop early byreturning True"
      break
    }
    count <- count + 1L
    if (new_global_min) {
      count <- 0L
    } else if (count > niter_success) {
      message <- "success condition satisfied"
      break
    }
  }
  list(x = lowest$x, fun = lowest$fun, message = message, nit = as.integer(i), success = isTRUE(lowest$success),
       nfev = nfev, njev = njev, nhev = nhev, minimization_failures = failures, lowest_optimization_result = lowest)
}

# --------------------------------------------------------------------------
# dual_annealing

.sp_dual_annealing <- function(func, lb, ub, maxiter = 1000L, minimizer_kwargs = NULL, initial_temp = 5230,
                               restart_temp_ratio = 2e-5, visit = 2.62, accept = -5.0, maxfun = 1e7,
                               no_local_search = FALSE, callback = NULL, x0 = NULL, ...) {
  lower <- as.numeric(lb); upper <- as.numeric(ub)
  if (!is.null(x0) && length(x0) != length(lower)) stop("Bounds size does not match x0", call. = FALSE)
  if (restart_temp_ratio <= 0 || restart_temp_ratio >= 1)
    stop("Restart temperature ratio has to be in range (0, 1)", call. = FALSE)
  if (any(!is.finite(c(lower, upper)))) stop("Some bounds values are inf values or nan values", call. = FALSE)
  if (!all(lower < upper)) stop("Bounds are not consistent min < max", call. = FALSE)
  n <- length(lower)
  fw <- new.env(parent = emptyenv())
  fw$nfev <- 0; fw$ngev <- 0; fw$nhev <- 0
  fun <- function(x) {
    fw$nfev <- fw$nfev + 1
    as.numeric(func(x))
  }
  # LocalSearchWrapper
  mk <- minimizer_kwargs %||% list()
  if (!length(mk)) {
    mk <- list(method = "L-BFGS-B", options = list(maxiter = min(max(n * 6L, 100L), 1000L)),
               bounds = list(lb = lower, ub = upper))
  }
  local_search <- function(x, e) {
    mres <- .sp_minimize(fun, x, mk$method %||% "BFGS", jac = mk$jac, hess = mk$hess, bounds = mk$bounds,
                         options = mk$options %||% list())
    if (!is.null(mres$njev)) fw$ngev <- fw$ngev + mres$njev
    if (!is.null(mres$nhev)) fw$nhev <- fw$nhev + mres$nhev
    ok <- all(is.finite(mres$x)) && is.finite(mres$fun) && all(mres$x >= lower) && all(mres$x <= upper)
    if (ok && mres$fun < e) list(e = mres$fun, x = mres$x) else list(e = e, x = x)
  }
  # VisitingDistribution
  bound_range <- upper - lower
  factor2 <- exp((4.0 - visit) * log(visit - 1.0))
  factor3 <- exp((2.0 - visit) * log(2.0) / (visit - 1.0))
  factor4_p <- sqrt(pi) * factor2 / (factor3 * (3.0 - visit))
  factor5 <- 1.0 / (visit - 1.0) - 0.5
  d1 <- 2.0 - factor5
  factor6 <- pi * (1.0 - factor5) / sin(pi * (1.0 - factor5)) / exp(lgamma(d1))
  visit_fn <- function(temperature, dim) {
    z <- matrix(stats::rnorm(2L * dim), dim, 2L, byrow = TRUE)
    xx <- z[, 1L]; yy <- z[, 2L]
    factor1 <- exp(log(temperature) / (visit - 1.0))
    factor4 <- factor4_p * factor1
    xx <- xx * exp(-(visit - 1.0) * log(factor6 / factor4) / (3.0 - visit))
    den <- exp((visit - 1.0) * log(abs(yy)) / (3.0 - visit))
    xx / den
  }
  tail_limit <- 1e8
  min_visit <- 1e-10
  visiting <- function(x, step, temperature) {
    if (step < n) {
      visits <- visit_fn(temperature, n)
      us <- stats::runif(2L)
      visits[visits > tail_limit] <- tail_limit * us[1L]
      visits[visits < -tail_limit] <- -tail_limit * us[2L]
      xv <- visits + x
      a <- xv - lower
      b <- .sp_fmod(a, bound_range) + bound_range
      xv <- .sp_fmod(b, bound_range) + lower
      small <- abs(xv - lower) < min_visit
      xv[small] <- xv[small] + 1e-10
    } else {
      xv <- x
      v <- visit_fn(temperature, 1L)[1L]
      if (v > tail_limit) v <- tail_limit * stats::runif(1L) else if (v < -tail_limit) v <- -tail_limit * stats::runif(1L)
      k <- step - n + 1L
      xv[k] <- v + x[k]
      a <- xv[k] - lower[k]
      b <- .sp_fmod(a, bound_range[k]) + bound_range[k]
      xv[k] <- .sp_fmod(b, bound_range[k]) + lower[k]
      if (abs(xv[k] - lower[k]) < min_visit) xv[k] <- xv[k] + min_visit
    }
    xv
  }
  # EnergyState
  es <- new.env(parent = emptyenv())
  es$ebest <- NULL; es$xbest <- NULL
  reset <- function(x0 = NULL) {
    es$current_location <- if (is.null(x0)) stats::runif(n, lower, upper) else as.numeric(x0)
    reinit <- 0L
    repeat {
      es$current_energy <- fun(es$current_location)
      if (is.finite(es$current_energy)) break
      if (reinit >= 1000L)
        stop("Stopping algorithm because function create NaN or (+/-) infinity values even with trying new ",
             "random parameters", call. = FALSE)
      es$current_location <- stats::runif(n, lower, upper)
      reinit <- reinit + 1L
    }
    if (is.null(es$ebest) && is.null(es$xbest)) {
      es$ebest <- es$current_energy
      es$xbest <- es$current_location
    }
  }
  update_best <- function(e, x, context) {
    es$ebest <- e
    es$xbest <- x
    if (is.function(callback) && isTRUE(callback(x, e, context)))
      return("Callback function requested to stop early by returning True")
    NULL
  }
  update_current <- function(e, x) {
    es$current_energy <- e
    es$current_location <- x
  }
  reset(x0)
  # StrategyChain
  sc <- new.env(parent = emptyenv())
  sc$emin <- es$current_energy; sc$xmin <- es$current_location
  sc$not_improved_idx <- 0L; sc$not_improved_max_idx <- 1000L; sc$temperature_step <- 0
  K <- 100 * n
  accept_reject <- function(j, e, xv) {
    r <- stats::runif(1L)
    pqv_temp <- 1.0 - ((1.0 - accept) * (e - es$current_energy) / sc$temperature_step)
    pqv <- if (pqv_temp <= 0) 0 else exp(log(pqv_temp) / (1.0 - accept))
    if (r <= pqv) {
      update_current(e, xv)
      sc$xmin <- es$current_location
    }
    if (sc$not_improved_idx >= sc$not_improved_max_idx) {
      if (j == 0L || es$current_energy < sc$emin) {
        sc$emin <- es$current_energy
        sc$xmin <- es$current_location
      }
    }
  }
  run_chain <- function(step, temperature) {
    sc$temperature_step <- temperature / (step + 1)
    sc$not_improved_idx <- sc$not_improved_idx + 1L
    for (j in 0:(2L * n - 1L)) {
      if (j == 0L) sc$improved <- step == 0L
      xv <- visiting(es$current_location, j, temperature)
      e <- fun(xv)
      if (e < es$current_energy) {
        update_current(e, xv)
        if (e < es$ebest) {
          val <- update_best(e, xv, 0L)
          if (!is.null(val)) return(val)
          sc$improved <- TRUE
          sc$not_improved_idx <- 0L
        }
      } else {
        accept_reject(j, e, xv)
      }
      if (fw$nfev >= maxfun) return("Maximum number of function call reached during annealing")
    }
    NULL
  }
  run_local <- function() {
    if (isTRUE(sc$improved)) {
      ls <- local_search(es$xbest, es$ebest)
      if (ls$e < es$ebest) {
        sc$not_improved_idx <- 0L
        val <- update_best(ls$e, ls$x, 1L)
        if (!is.null(val)) return(val)
        update_current(ls$e, ls$x)
      }
      if (fw$nfev >= maxfun) return("Maximum number of function call reached during local search")
    }
    do_ls <- FALSE
    if (K < 90 * n) {
      pls <- exp(K * (es$ebest - es$current_energy) / sc$temperature_step)
      if (pls >= stats::runif(1L)) do_ls <- TRUE
    }
    if (sc$not_improved_idx >= sc$not_improved_max_idx) do_ls <- TRUE
    if (do_ls) {
      ls <- local_search(sc$xmin, sc$emin)
      sc$xmin <- ls$x
      sc$emin <- ls$e
      sc$not_improved_idx <- 0L
      sc$not_improved_max_idx <- n
      if (ls$e < es$ebest) {
        val <- update_best(sc$emin, sc$xmin, 2L)
        if (!is.null(val)) return(val)
        update_current(ls$e, ls$x)
      }
      if (fw$nfev >= maxfun) return("Maximum number of function call reached during dual annealing")
    }
    NULL
  }
  temperature_restart <- initial_temp * restart_temp_ratio
  need_to_stop <- FALSE
  iteration <- 0L
  message <- character()
  success <- TRUE
  t1 <- exp((visit - 1) * log(2.0)) - 1.0
  while (!need_to_stop) {
    for (i in 0:(maxiter - 1L)) {
      s <- i + 2.0
      t2 <- exp((visit - 1) * log(s)) - 1.0
      temperature <- initial_temp * t1 / t2
      if (iteration >= maxiter) {
        message <- c(message, "Maximum number of iteration reached")
        need_to_stop <- TRUE
        break
      }
      if (temperature < temperature_restart) {
        reset()
        break
      }
      val <- run_chain(i, temperature)
      if (!is.null(val)) {
        message <- c(message, val)
        need_to_stop <- TRUE
        success <- FALSE
        break
      }
      if (!no_local_search) {
        val <- run_local()
        if (!is.null(val)) {
          message <- c(message, val)
          need_to_stop <- TRUE
          success <- FALSE
          break
        }
      }
      iteration <- iteration + 1L
    }
  }
  list(x = es$xbest, fun = es$ebest, nit = iteration, nfev = as.integer(fw$nfev), njev = as.integer(fw$ngev),
       nhev = as.integer(fw$nhev), message = message, success = success, status = 0L)
}

# --------------------------------------------------------------------------
# the OptimizationSciPy wrappers

# NumPy generator arguments (rng, seed) become a set.seed() call
.scipy_rng_kw <- function(kw) {
  seed <- kw$rng %||% kw$seed
  kw$rng <- NULL
  kw$seed <- NULL
  list(kw = kw, seed = if (is.numeric(seed)) seed)
}

.scipy_run_de <- function(f, lb, ub, cons, a, opts, kw) {
  r <- .scipy_rng_kw(kw)
  args <- utils::modifyList(list(maxiter = as.integer(opts$maxiters %||% 1000L), popsize = 15L, atol = 0.0,
                                 tol = 0.01, mutation = c(0.5, 1), recombination = 0.7, polish = TRUE,
                                 init = "latinhypercube", updating = "immediate", workers = 1L), r$kw)
  .scipy_with_seed(r$seed, do.call(.sp_differential_evolution, c(list(function(x) as.numeric(f(x)), lb, ub), args)))
}

.scipy_run_basinhopping <- function(f, x0, opts, kw) {
  r <- .scipy_rng_kw(kw)
  args <- utils::modifyList(list(niter = as.integer(opts$maxiters %||% 100L), T = 1.0, stepsize = 0.5,
                                 interval = 50L), r$kw)
  res <- .scipy_with_seed(r$seed, do.call(.sp_basinhopping, c(list(function(x) as.numeric(f(x)), x0), args)))
  low <- res$lowest_optimization_result
  res$success <- isTRUE(low$success)
  res$nit <- low$nit
  res
}

.scipy_run_dual_annealing <- function(f, lb, ub, opts, kw) {
  r <- .scipy_rng_kw(kw)
  args <- utils::modifyList(list(maxiter = as.integer(opts$maxiters %||% 1000L), initial_temp = 5230.0,
                                 restart_temp_ratio = 2.0e-5, visit = 2.62, accept = -5.0, maxfun = 1.0e7,
                                 no_local_search = FALSE), r$kw)
  .scipy_with_seed(r$seed, do.call(.sp_dual_annealing, c(list(function(x) as.numeric(f(x)), lb, ub), args)))
}
