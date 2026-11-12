# The single objective and combinatorial algorithms of Metaheuristics.jl
# (3.6), ported to R on the framework of optimization-metaheuristics-base.R. Each
# algorithm registers its initialize!, update_state!, final_stage! and (when
# it is not the template) stop_criteria! in .mh_impls.

# The implementation of an algorithm: initialize, update_state, final_stage
# and stop_criteria (the template by default). This file is collated before
# the other files of the family, so the registry lives here.
.mh_impls <- new.env(parent = emptyenv())

.mh_register <- function(name, initialize, update, final = NULL, stop_criteria = NULL)
  assign(name, list(initialize = initialize, update = update, final = final, stop_criteria = stop_criteria),
         envir = .mh_impls)

.mh_set_budget <- function(op, calls, iterations_from_calls) {
  if (op$f_calls_limit == 0) op$f_calls_limit <- calls
  if (op$iterations == 0L) op$iterations <- as.integer(iterations_from_calls(op$f_calls_limit))
}

# ---------------------------------------------------------------------------
# ECA: evolutionary centers algorithm

.mh_fitness_to_mass <- function(f) {
  m <- min(f); M <- max(f)
  if (m == 0 && M == 0) return(rep(1, length(f)))
  if (m < 0) f <- 2 * abs(m) + f
  2 * max(f) - f
}

.mh_get_mass <- function(U) {
  if (is.null(U[[1L]]$sv)) return(.mh_fitness_to_mass(.mh_fvals(U)))
  fit <- .mh_fvals(U)
  if (is.matrix(fit)) fit <- rowSums(fit)
  M <- 2 * max(abs(fit))
  .mh_fitness_to_mass(fit + max(M, 100) * .mh_svs(U))
}

.mh_center <- function(U, mass) {
  cc <- mass[1L] * U[[1L]]$x
  for (k in seq_along(U)[-1L]) cc <- cc + mass[k] * U[[k]]$x
  cc / .mh_seqsum(mass)
}

.mh_getU_ids <- function(K, I, i, N) if (i <= N - K) I[i:(K + i)] else I[(i:(K + i)) %% N + 1L]

.mh_adapt_crossover <- function(p_cr, M) {
  pb <- p_cr[which.min(M)]
  for (i in seq_along(p_cr)) {
    if (M[i] > 0.3) {
      pn <- abs(pb + 0.3 * stats::rnorm(1L))
      if (pn > 1) pn <- 1
      if (pn < 0) pn <- pb
      p_cr[i] <- pn
    }
  }
  p_cr
}

.mh_eca_solution <- function(st, par, op, pr, I, i) {
  p <- st$f_calls / op$f_calls_limit
  U <- st$population[.mh_getU_ids(par$K, I, i, par$N)]
  cc <- .mh_center(U, .mh_get_mass(U))
  uw <- .mh_argworst(U); ub <- .mh_argbest(U)
  eta <- par$eta_max * stats::runif(1L)
  xi <- st$population[[i]]$x
  y <- if (p < par$p_exploit) xi + eta * (cc - U[[uw]]$x) else xi + eta * (st$best_sol$x - cc)
  xb <- U[[ub]]$x
  M <- numeric(length(y))
  for (j in seq_along(y)) if (stats::runif(1L) < par$p_cr[j]) { y[j] <- xb[j]; M[j] <- 1 }
  list(y = .mh_evo_boundary_repairer(y, cc, pr$lb, pr$ub), M = M)
}

.mh_eca_env_selection <- function(pop, N) {
  new <- pop[(N + 1L):length(pop)]
  P <- pop[seq_len(N)]
  surv <- seq_len(N)
  for (i in seq_along(new)) {
    if (!.mh_is_better(new[[i]], P[[i]])) next
    wi <- .mh_argworst(P)
    surv[wi] <- N + i
    P[[wi]] <- new[[i]]
  }
  pop[sort(surv)]
}

.mh_register("ECA",
  initialize = function(st, par, pr, info, op) {
    D <- pr$D
    if (par$N <= par$K) par$N <- as.integer(min(par$K * D, 1000))
    .mh_set_budget(op, 10000 * D, function(fc) fc %/% par$N + 1)
    # Metaheuristics.jl assigns N_init to a local variable here, so that a
    # resized population shrinks to 2K at once; the field is set instead
    par$N_init <- par$N
    par$p_cr <- if (par$adaptive) stats::runif(D) else par$p_bin * rep(1, D)
    .mh_gen_initial_state(pr, par, info, op, st)
  },
  update = function(st, par, pr, info, op) {
    I <- sample.int(par$N)
    D <- pr$D
    Mcr_fail <- numeric(D)
    if (op$parallel_evaluation) Xn <- matrix(0, par$N, D)
    for (i in seq_len(par$N)) {
      r <- .mh_eca_solution(st, par, op, pr, I, i)
      if (op$parallel_evaluation) { Xn[i, ] <- r$y; next }
      sol <- .mh_create_solution(r$y, pr)
      if (.mh_is_better(sol, st$population[[i]])) {
        wi <- .mh_argworst(st$population)
        st$population[[wi]] <- sol
        if (.mh_is_better(sol, st$best_sol)) st$best_sol <- sol
      } else if (par$adaptive) {
        Mcr_fail <- Mcr_fail + r$M
      }
    }
    if (op$parallel_evaluation) {
      st$population <- .mh_eca_env_selection(c(st$population, .mh_create_solutions(Xn, pr)), par$N)
      st$best_sol <- .mh_get_best(st$population)
    }
    if (par$adaptive) par$p_cr <- .mh_adapt_crossover(par$p_cr, Mcr_fail / par$N)
    if (par$resize_population) {
      p <- st$f_calls / op$f_calls_limit
      K <- par$K
      par$N <- as.integer(max(2 * K + round((1 - p) * (par$N_init - 2 * K)), 2 * K))
      if (length(st$population) != par$N) st$population <- st$population[.mh_sortperm(st$population)[seq_len(par$N)]]
    }
  })

# ---------------------------------------------------------------------------
# DE, epsilonDE and SHADE: differential evolution

.mh_de_reproduction <- function(st, de, pr) {
  pop <- st$population
  N <- de$N
  xb <- st$best_sol$x
  X <- matrix(0, N, pr$D)
  for (i in seq_len(N)) {
    # Metaheuristics.jl passes the index 1 as the best solution, so that the
    # best strategies use an arbitrary member; the best one is used here
    v <- .mh_de_mutation(pop, de$F, de$strategy, 0L)
    v <- .mh_de_crossover(v, pop[[i]]$x, de$CR)
    X[i, ] <- .mh_evo_boundary_repairer(v, xb, pr$lb, pr$ub)
  }
  X
}

# environmental_selection! of the DE family: each trial against its target,
# the survivors kept in index order (old ones first)
.mh_de_selection <- function(pop, N, better) {
  surv <- vapply(seq_len(N), function(i) if (better(pop[[N + i]], pop[[i]])) N + i else i, 0L)
  pop[sort(surv)]
}

.mh_de_initialize <- function(st, par, pr, info, op) {
  D <- pr$D
  if (par$N <= 5L) par$N <- as.integer(10 * D)
  if (par$CR < 0 || par$CR > 1) par$CR <- 0.5
  .mh_set_budget(op, 10000 * D, function(fc) fc %/% par$N + 1)
  .mh_gen_initial_state(pr, par, info, op, st)
}

.mh_register("DE",
  initialize = .mh_de_initialize,
  update = function(st, par, pr, info, op) {
    if (par$F_min < par$F_max) par$F <- par$F_min + (par$F_max - par$F_min) * stats::runif(1L)
    if (par$CR_min < par$CR_max) par$CR <- par$CR_min + (par$CR_max - par$CR_min) * stats::runif(1L)
    new <- .mh_create_solutions(.mh_de_reproduction(st, par, pr), pr, op$h_tol)
    st$population <- .mh_de_selection(c(st$population, new), par$N, .mh_is_better)
    st$best_sol <- .mh_get_best(st$population)
  })

.mh_epsilon_better <- function(a, b, eps) {
  pa <- .mh_sv(a); pb <- .mh_sv(b)
  if ((pa <= eps && pb <= eps) || pa == pb) return(isTRUE(a$f < b$f))
  pa < pb
}

.mh_register("epsilonDE",
  initialize = function(st, par, pr, info, op) {
    if (!is.environment(par$de)) par$de <- list2env(par$de, parent = emptyenv())
    st <- .mh_de_initialize(st, par$de, pr, info, op)
    elite <- .mh_sortperm(st$population)
    theta <- round(0.2 * length(elite))
    cand <- elite[seq_len(theta)]
    s <- st$population[[cand[sample.int(length(cand), 1L)]]]
    par$eps0 <- .mh_sv(s)
    par$Tc <- round(0.2 * op$iterations)
    par$N <- par$de$N
    st
  },
  update = function(st, par, pr, info, op) {
    t <- st$iteration
    par$eps <- if (t < par$Tc) par$eps0 * (1 - t / par$Tc)^par$cp else 0
    new <- .mh_create_solutions(.mh_de_reproduction(st, par$de, pr), pr, op$h_tol)
    eps <- par$eps
    st$population <- .mh_de_selection(c(st$population, new), par$N, function(a, b) .mh_epsilon_better(a, b, eps))
    st$best_sol <- .mh_get_best(st$population)
  })

.mh_randc <- function(n) stats::rnorm(n) / stats::rnorm(n)

.mh_register("SHADE",
  initialize = function(st, par, pr, info, op) {
    .mh_set_budget(op, 10000 * pr$D, function(fc) fc %/% par$N + 1)
    par$MCR <- rep(0.5, par$N)
    par$MF <- rep(0.5, par$N)
    .mh_gen_initial_state(pr, par, info, op, st)
  },
  update = function(st, par, pr, info, op) {
    N <- par$N; D <- pr$D
    pop <- st$population[.mh_sortperm(st$population)]
    # the p best member: a random one among the best max(2, ceil(0.2 N U));
    # Metaheuristics.jl takes the member at that rank instead
    pb <- vapply(seq_len(N), function(i) sample.int(max(2L, ceiling(stats::runif(1L) * 0.2 * N)), 1L), 0L)
    Xpb <- .mh_positions(pop[pb])
    Xr1 <- .mh_positions(pop[sample.int(N, N, replace = TRUE)])
    # r2 from the population and the archive (Metaheuristics.jl draws it
    # from the population only)
    P <- c(pop, par$archive)
    Xr2 <- .mh_positions(P[sample.int(length(P), N, replace = TRUE)])
    CR <- pmin(pmax(stats::rnorm(N) * sqrt(0.1) + par$MCR[sample.int(N, N, replace = TRUE)], 0), 1)
    F <- .mh_randc(N) * sqrt(0.1) + par$MF[sample.int(N, N, replace = TRUE)]
    for (k in 1:100) {
      m <- F <= 0
      if (!any(m)) break
      F[m] <- pmin(1, .mh_randc(sum(m)) * sqrt(0.1) + par$MF[sample.int(N, sum(m), replace = TRUE)])
    }
    F <- pmin(pmax(F, 0), 1)
    par$CR <- CR
    par$F <- F
    mask <- matrix(stats::runif(N * D), N, D) < CR
    X <- .mh_positions(pop)
    Fm <- matrix(F, N, D)
    X[mask] <- X[mask] + Fm[mask] * (Xpb[mask] - X[mask] + Xr1[mask] - Xr2[mask])
    for (i in seq_len(N)) X[i, ] <- .mh_reset_to_violated_bounds(X[i, ], pr$lb, pr$ub)
    new <- .mh_create_solutions(X, pr, op$h_tol)
    surv <- integer(N)
    delta <- numeric(0)
    for (i in seq_len(N)) {
      h <- new[[i]]
      if (.mh_is_better(h, pop[[i]])) {
        surv[i] <- N + i
        vp <- .mh_sv(pop[[i]]); vc <- .mh_sv(h)
        delta <- c(delta, if (vp != vc) vp - vc else pop[[i]]$f - h$f)
      } else {
        surv[i] <- i
      }
    }
    if (length(delta)) {
      win <- surv[surv > N] - N
      valid <- is.finite(delta) & delta > 0 & F[win] > 0
      if (any(valid)) {
        w <- delta[valid] / sum(delta[valid])
        SCR <- CR[win[valid]]; SF <- F[win[valid]]
        k <- par$k
        par$MCR[k] <- sum(w * SCR)
        par$MF[k] <- sum(w * SF^2) / sum(w * SF)
        par$k <- k %% length(par$MCR) + 1L
      }
      # the archive keeps the replaced parents and at most N of them
      # (Metaheuristics.jl stores the trials and deletes N of them)
      arch <- c(par$archive, pop[win])
      if (length(arch) > N) arch <- arch[-sort(sample.int(length(arch))[seq_len(length(arch) - N)])]
      par$archive <- arch
    }
    st$population <- c(pop, new)[sort(surv)]
    st$best_sol <- .mh_get_best(st$population)
  })

# ---------------------------------------------------------------------------
# PSO: particle swarm optimization

.mh_register("PSO",
  initialize = function(st, par, pr, info, op) {
    D <- pr$D
    if (par$N == 0L) par$N <- as.integer(10 * D)
    .mh_set_budget(op, 10000 * D, function(fc) fc %/% par$N + 1)
    st <- .mh_gen_initial_state(pr, par, info, op, st)
    par$v <- matrix(0, par$N, D)
    par$flock <- st$population
    st
  },
  update = function(st, par, pr, info, op) {
    xg <- st$best_sol$x
    Xn <- matrix(0, par$N, pr$D)
    for (i in seq_len(par$N)) {
      x <- par$flock[[i]]$x
      r1 <- par$C1 * stats::runif(1L)
      r2 <- par$C2 * stats::runif(1L)
      v <- par$omega * par$v[i, ] + r1 * (st$population[[i]]$x - x) + r2 * (xg - x)
      par$v[i, ] <- v
      Xn[i, ] <- .mh_reset_to_violated_bounds(x + v, pr$lb, pr$ub)
    }
    new <- .mh_create_solutions(Xn, pr, op$h_tol)
    for (i in seq_along(new)) {
      sol <- new[[i]]
      if (.mh_is_better(sol, st$population[[i]])) {
        st$population[[i]] <- sol
        if (.mh_is_better(sol, st$best_sol)) st$best_sol <- sol
      }
      par$flock[[i]] <- sol
    }
  })

# ---------------------------------------------------------------------------
# ABC: artificial bee colony. A bee is a solution with its fitness `fit`
# and trial counter `t`.

.mh_abc_fit <- function(fx) if (fx >= 0) 1 / (1 + fx) else 1 + abs(fx)

.mh_bee <- function(sol, t = 0L) { sol$fit <- .mh_abc_fit(sol$f); sol$t <- t; sol }

.mh_abc_getK <- function(i, N) {
  j <- sample.int(N, 1L)
  while (j == i) j <- sample.int(N, 1L)
  j
}

.mh_abc_update_bee <- function(st, i, k, pr) {
  bee <- st$population[[i]]
  phi <- -1 + 2 * stats::runif(1L)
  xn <- bee$x + phi * (bee$x - st$population[[k]]$x)
  xn <- .mh_replace_with_random_in_bounds(xn, pr$lb, pr$ub)
  new <- .mh_create_solution(xn, pr)
  st$population[[i]] <- if (.mh_is_better(new, bee)) .mh_bee(new) else { bee$t <- bee$t + 1L; bee }
}

.mh_abc_best_bee <- function(bees) {
  best <- bees[[1L]]
  for (b in bees) if (isTRUE(b$f < best$f)) best <- b
  best
}

.mh_register("ABC",
  initialize = function(st, par, pr, info, op) {
    if (op$parallel_evaluation)
      stop("ABC is not supporting parallel evaluation. Put `options.parallel_evaluation=false`", call. = FALSE)
    D <- pr$D
    if (op$f_calls_limit == 0) {
      op$f_calls_limit <- 10000 * D
      op$iterations <- as.integer(par$N + (10000 * D) %/% par$N)
    }
    s0 <- .mh_gen_initial_state(pr, par, info, op, st)
    bees <- lapply(s0$population, .mh_bee)
    .mh_state(.mh_abc_best_bee(bees), bees, f_calls = length(bees))
  },
  update = function(st, par, pr, info, op) {
    N <- length(st$population)
    for (i in sample.int(N)[seq_len(par$Ne)]) .mh_abc_update_bee(st, i, .mh_abc_getK(i, N), pr)
    sum_f <- sum(vapply(st$population, function(b) b$fit, 0))
    for (k in seq_len(par$No)) {
      r <- stats::runif(1L)
      rs <- 0
      j <- N
      for (i in seq_len(N)) {
        rs <- rs + st$population[[i]]$fit / sum_f
        if (r <= rs) { j <- i; break }
      }
      .mh_abc_update_bee(st, j, .mh_abc_getK(j, N), pr)
    }
    cand <- .mh_abc_best_bee(st$population)
    best <- if (isTRUE(cand$f < st$best_sol$f)) cand else st$best_sol
    scouts <- which(vapply(st$population, function(b) b$t >= par$limit, TRUE))
    for (i in scouts) {
      b <- st$population[[i]]
      s <- .mh_create_solution(pr$lb + (pr$ub - pr$lb) * stats::runif(pr$D), pr)
      # the scout keeps its old fitness, as in Metaheuristics.jl
      s$fit <- b$fit; s$t <- 0L
      st$population[[i]] <- s
    }
    st$f_calls <- st$f_calls + par$Ne + par$No + length(scouts)
    st$best_sol <- best
    .mh_template_stop(st, par, pr, info, op)
  })

# ---------------------------------------------------------------------------
# SA: simulated annealing

.mh_sa_newsol <- function(y, mu) ((1 + mu)^abs(y) - 1) / mu * sign(y)

.mh_register("SA",
  initialize = function(st, par, pr, info, op) {
    if (!length(par$x_initial)) par$x_initial <- pr$lb + (pr$ub - pr$lb) * stats::runif(pr$D)
    .mh_set_budget(op, 10000 * pr$D, function(fc) fc %/% par$N)
    x <- par$x_initial
    best <- .mh_create_solution(x, pr)
    if (!is.null(best$sv))
      stop("SA has no method for constrained solutions (xfgh_solution); Metaheuristics.jl fails the same way",
           call. = FALSE)
    st <- .mh_state(best, list(best))
    par$x <- x
    par$fx <- best$f
    st$f_calls <- 1L
    st
  },
  update = function(st, par, pr, info, op) {
    T <- st$f_calls / op$f_calls_limit
    mu <- 10^(100 * T)
    l <- pr$lb; u <- pr$ub
    for (i in seq_len(par$N)) {
      dx <- .mh_sa_newsol(2 * stats::runif(length(par$x)) - 1, mu) * (u - l)
      x1 <- .mh_reset_to_violated_bounds(par$x + dx, l, u)
      fx1 <- as.numeric(.mh_evaluate(x1, pr))
      st$f_calls <- st$f_calls + 1L
      df <- fx1 - par$fx
      if (isTRUE(df < 0) || isTRUE(stats::runif(1L) < exp(-T * df / abs(par$fx) / par$tol_fun))) {
        par$x <- x1
        par$fx <- fx1
      }
      if (isTRUE(fx1 < st$best_sol$f)) st$best_sol <- list(x = x1, f = fx1)
      st$population <- list(st$best_sol)
      .mh_template_stop(st, par, pr, info, op)
      if (st$stop) break
    }
  },
  final = function(st, par, pr, info, op) st$population <- list(st$best_sol))

# ---------------------------------------------------------------------------
# WOA: whale optimization algorithm

.mh_register("WOA",
  initialize = function(st, par, pr, info, op) {
    if (op$iterations == 0L) op$iterations <- 500L
    if (op$f_calls_limit == 0) op$f_calls_limit <- as.numeric(op$iterations) * par$N
    .mh_gen_initial_state(pr, par, info, op, st)
  },
  update = function(st, par, pr, info, op) {
    lead <- st$best_sol$x
    Max_iter <- op$iterations
    N <- par$N; D <- pr$D
    t <- st$iteration
    a <- 2 - t * (2 / Max_iter)
    a2 <- -1 + t * (-1 / Max_iter)
    Xn <- matrix(0, N, D)
    for (i in seq_len(N)) {
      r1 <- stats::runif(1L); r2 <- stats::runif(1L)
      A <- 2 * a * r1 - a
      C <- 2 * r2
      b <- 1
      l <- (a2 - 1) * stats::runif(1L) + 1
      p <- stats::runif(1L)
      # the whale moves from its own position (Metaheuristics.jl starts
      # from the zero vector, which drops the position from the update)
      x <- st$population[[i]]$x
      for (j in seq_len(D)) {
        if (p < 0.5) {
          if (abs(A) >= 1) {
            xr <- st$population[[floor(N * stats::runif(1L) + 1)]]$x
            x[j] <- xr[j] - A * abs(C * xr[j] - x[j])
          } else {
            x[j] <- lead[j] - A * abs(C * lead[j] - x[j])
          }
        } else {
          x[j] <- abs(lead[j] - x[j]) * exp(b * l) * cos(l * 2 * pi) + lead[j]
        }
      }
      Xn[i, ] <- .mh_reset_to_violated_bounds(x, pr$lb, pr$ub)
    }
    st$population <- .mh_create_solutions(Xn, pr)
    sol <- .mh_get_best(st$population)
    if (.mh_is_better(sol, st$best_sol)) st$best_sol <- sol
  })

# ---------------------------------------------------------------------------
# CGSA: chaotic gravitational search algorithm

# chaos(index, curr_iter, max_iter, Value): the chaotic map value at
# curr_iter (the sequence is computed up to that point only)
.mh_chaos <- function(index, curr_iter, max_iter, Value) {
  x <- numeric(curr_iter + 1L)
  x[1L] <- if (index == 10) 0.6 else 0.7
  G <- 0
  for (i in seq_len(curr_iter)) {
    xi <- x[i]
    switch(as.character(index),
      "1" = { x[i + 1L] <- cos(i * acos(xi)); G <- ((xi + 1) * Value) / 2 },
      "2" = { x[i + 1L] <- (xi + 0.2 - (0.5 / (2 * pi)) * sin(2 * pi * xi)) %% 1; G <- xi * Value },
      "3" = { x[i + 1L] <- if (xi == 0) 0 else (1 / xi) %% 1; G <- xi * Value },
      "4" = { x[i + 1L] <- sin((0.7 * pi) / xi); G <- ((xi + 1) * Value) / 2 },
      "5" = { x[i + 1L] <- 4 * xi * (1 - xi); G <- xi * Value },
      "6" = {
        P <- 0.4
        nx <- 0
        if (xi >= 0 && xi < P) nx <- xi / P
        if (xi >= P && xi < 0.5) nx <- (xi - P) / (0.5 - P)
        if (xi >= 0.5 && xi < 1 - P) nx <- (1 - P - xi) / (0.5 - P)
        if (xi >= 1 - P && xi < 1) nx <- (1 - xi) / P
        x[i + 1L] <- nx
        G <- xi * Value
      },
      "7" = { x[i + 1L] <- sin(pi * xi); G <- xi * Value },
      "8" = {
        x[i + 1L] <- 1.07 * (7.86 * xi - 23.31 * (xi * xi) + 28.75 * (xi * xi * xi) - 13.302875 * .mh_jl_pow4(xi))
        G <- xi * Value
      },
      "9" = { x[i + 1L] <- 2.3 * (xi * xi) * sin(pi * xi); G <- xi * Value },
      "10" = {
        nx <- 0
        if (xi < 0.7) nx <- xi / 0.7
        if (xi >= 0.7) nx <- (10 / 3) * (1 - xi)
        x[i + 1L] <- nx
        G <- xi * Value
      },
      G <- 0)
  }
  if (curr_iter > max_iter) stop("chaos: iteration beyond max_iter", call. = FALSE)
  G
}

# x^4 as Julia computes it (pow_body, a compensated squaring), with the
# exact products of Dekker's algorithm
.mh_two_mul <- function(a, b) {
  p <- a * b
  ca <- 134217729 * a; ah <- ca - (ca - a); al <- a - ah
  cb <- 134217729 * b; bh <- cb - (cb - b); bl <- b - bh
  c(p, ((ah * bh - p) + ah * bl + al * bh) + al * bl)
}

.mh_jl_pow4 <- function(x) {
  s <- .mh_two_mul(x, x)
  err <- s[1L] * 2 * s[2L]
  q <- .mh_two_mul(s[1L], s[1L])
  q[1L] + (q[2L] + err)
}

.mh_mass_calculation <- function(fit) {
  Fmax <- max(fit); Fmin <- min(fit)
  M <- if (Fmax == Fmin) rep(1, length(fit)) else (fit - Fmax) / (Fmin - Fmax)
  M / sum(M)
}

.mh_gconstant <- function(iteration, max_it) 100 * exp(-20 * iteration / max_it)

.mh_gfield <- function(M, X, G, Rnorm, Rpower, ElitistCheck, iteration, max_it) {
  N <- nrow(X); D <- ncol(X)
  kbest <- if (ElitistCheck == 1) round(N * (2 + (1 - iteration / max_it) * (100 - 2)) / 100) else N
  ds <- order(M, decreasing = TRUE)
  E <- matrix(0, N, D)
  for (i in seq_len(N)) {
    for (ii in seq_len(kbest)) {
      j <- ds[ii]
      if (j == i) next
      dv <- X[i, ] - X[j, ]
      R <- if (Rnorm == 2) sqrt(sum(dv * dv)) else sum(abs(dv)^Rnorm)^(1 / Rnorm)
      E[i, ] <- E[i, ] + stats::runif(D) * M[j] * ((X[j, ] - X[i, ]) / (R^Rpower + .Machine$double.eps))
    }
  }
  E * G
}

.mh_register("CGSA",
  initialize = function(st, par, pr, info, op) {
    if (op$iterations == 0L) op$iterations <- 500L
    st <- .mh_gen_initial_state(pr, par, info, op, st)
    if (op$f_calls_limit == 0) op$f_calls_limit <- as.numeric(op$iterations) * par$N
    if (!length(par$V)) par$V <- matrix(0, par$N, pr$D)
    if (!length(par$X)) par$X <- .mh_positions(st$population)
    if (!length(par$fitness)) par$fitness <- .mh_fvals(st$population)
    st
  },
  update = function(st, par, pr, info, op) {
    it <- st$iteration
    max_it <- op$iterations
    chValue <- par$wMax - it * ((par$wMax - par$wMin) / max_it)
    M <- .mh_mass_calculation(par$fitness)
    G <- .mh_gconstant(it, max_it)
    if (1 <= par$chaosIndex && par$chaosIndex <= 10) G <- G + .mh_chaos(par$chaosIndex, it, max_it, chValue)
    a <- .mh_gfield(M, par$X, G, par$Rnorm, par$Rpower, par$ElitistCheck, it, max_it)
    V <- matrix(stats::runif(length(par$V)), nrow(par$V)) * par$V + a
    X <- par$X + V
    for (i in seq_len(nrow(X))) X[i, ] <- .mh_reset_to_violated_bounds(X[i, ], pr$lb, pr$ub)
    st$population <- .mh_create_solutions(X, pr)
    par$X <- X
    par$fitness <- .mh_fvals(st$population)
    par$V <- V
    cb <- .mh_get_best(st$population)
    if (.mh_is_better(cb, st$best_sol)) st$best_sol <- cb
    .mh_template_stop(st, par, pr, info, op)
  })

# ---------------------------------------------------------------------------
# CSO: competitive swarm optimizer

.mh_register("CSO",
  initialize = function(st, par, pr, info, op) {
    D <- pr$D
    if (par$N == 0L) par$N <- as.integer(min(max(10 * D, 10), 500))
    .mh_set_budget(op, 10000 * D, function(fc) fc %/% (par$N %/% 2L) + 1)
    par$V_loser <- matrix(0, par$N %/% 2L, D)
    .mh_gen_initial_state(pr, par, info, op, st)
  },
  update = function(st, par, pr, info, op) {
    N <- par$N; D <- pr$D; h <- N %/% 2L
    pop <- st$population
    rk <- sample.int(N)
    loser <- rk[seq_len(h)]
    winner <- rk[(length(rk) - h + 1L):length(rk)]
    m <- vapply(seq_len(h), function(i) .mh_is_better(pop[[loser[i]]], pop[[winner[i]]]), TRUE)
    tmp <- loser[m]; loser[m] <- winner[m]; winner[m] <- tmp
    Xl <- .mh_positions(pop[loser]); Xw <- .mh_positions(pop[winner])
    R1 <- matrix(stats::runif(h * D), h, D); R2 <- matrix(stats::runif(h * D), h, D); R3 <- matrix(stats::runif(h * D), h, D)
    xm <- colMeans(.mh_positions(pop))
    par$V_loser <- R1 * par$V_loser + R2 * (Xw - Xl) + par$phi * R3 * (matrix(xm, h, D, byrow = TRUE) - Xl)
    Xl <- Xl + par$V_loser
    for (i in seq_len(h)) Xl[i, ] <- .mh_reset_to_violated_bounds(Xl[i, ], pr$lb, pr$ub)
    pop[loser] <- .mh_create_solutions(Xl, pr, op$h_tol)
    st$population <- pop
    st$best_sol <- .mh_get_best(pop)
  })

# ---------------------------------------------------------------------------
# GA and BRKGA: the genetic algorithm with operator objects

.mh_op <- function(cl, ...) structure(list(...), class = c(cl, "MetaheuristicsOperator"))

.mh_rank_solutions <- function(pop) {
  if (.mh_is_mo(pop[[1L]])) return(.mh_non_dominated_sort(pop))
  order(.mh_sortperm(pop))
}

.mh_selection <- function(pop, op) {
  n <- length(pop)
  switch(class(op)[1L],
    TournamentSelection = {
      if (op$K < 1L) stop("TournamentSelection requires K >= 1.", call. = FALSE)
      fitness <- .mh_rank_solutions(pop)
      vapply(seq_len(op$N), function(j) {
        best <- sample.int(n, 1L)
        for (k in seq_len(op$K - 1L)) {
          i <- sample.int(n, 1L)
          if (fitness[best] > fitness[i]) best <- i
        }
        best
      }, 0L)
    },
    RouletteWheelSelection = {
      fitness <- .mh_fvals(pop)
      fmin <- min(fitness)
      shift <- if (fmin > 0) .Machine$double.eps else fmin - 1
      p <- cumsum(1 / (fitness - shift))
      p <- p / max(p)
      if (all(is.finite(p)) && !is.unsorted(p)) {
        vapply(seq_len(op$N), function(j) findInterval(stats::runif(1L), p) + 1L, 0L)
      } else {
        vapply(seq_len(op$N), function(j) which(p > stats::runif(1L))[1L], 0L)
      }
    },
    BiasedSelection = {
      ne <- op$num_elites; no <- op$num_offsprings
      parent <- integer(2L * no)
      parent[seq(1L, 2L * no, 2L)] <- sample.int(ne, no, replace = TRUE)
      parent[seq(2L, 2L * no, 2L)] <- ne + sample.int(n - ne, no, replace = TRUE)
      parent
    },
    stop("unknown selection ", class(op)[1L], call. = FALSE))
}

.mh_crossover <- function(parents, op, pr) {
  Q <- .mh_positions(parents)
  switch(class(op)[1L],
    SBX = {
      lb <- op$lb %mh_or% pr$lb; ub <- op$ub %mh_or% pr$ub
      n <- nrow(Q)
      if (n >= 2L) for (i in seq(1L, n - 1L, 2L)) {
        cc <- .mh_sbx(Q[i, ], Q[i + 1L, ], lb, ub, op$eta, op$p)
        Q[i, ] <- cc[[1L]]; Q[i + 1L, ] <- cc[[2L]]
      }
      Q
    },
    UniformCrossover = {
      n <- nrow(Q) %/% 2L
      O <- Q[seq_len(2L * n), , drop = FALSE]
      R <- matrix(stats::runif(n * ncol(O)), n, ncol(O))
      sw <- R <= op$p
      A <- O[seq_len(n), , drop = FALSE]; B <- O[n + seq_len(n), , drop = FALSE]
      O[seq_len(n), ] <- ifelse(sw, B, A)
      O[n + seq_len(n), ] <- ifelse(sw, A, B)
      O
    },
    BinomialCrossover = {
      n <- nrow(Q) %/% 2L
      R <- matrix(stats::runif(n * ncol(Q)), n, ncol(Q))
      first <- R <= op$p
      A <- Q[2L * seq_len(n) - 1L, , drop = FALSE]; B <- Q[2L * seq_len(n), , drop = FALSE]
      if (op$n_offsprings == 1L) ifelse(first, A, B) else rbind(ifelse(first, A, B), ifelse(first, B, A))
    },
    OrderCrossover = {
      N <- nrow(Q); D <- ncol(Q)
      s <- sample.int(D, N, replace = TRUE)
      if (N >= 2L) for (i in seq(1L, N - 1L, 2L)) {
        PA <- Q[i, ]; PB <- Q[i + 1L, ]
        fill <- function(P1, P2, k) {
          tail <- setdiff(P2, P1[seq_len(k)])
          if (length(tail) != D - k)
            stop("OrderCrossover requires each row to be a permutation of the same values.", call. = FALSE)
          c(P1[seq_len(k)], tail)
        }
        Q[i, ] <- fill(PA, PB, s[i])
        Q[i + 1L, ] <- fill(PB, PA, s[i])
      }
      Q
    },
    stop("unknown crossover ", class(op)[1L], call. = FALSE))
}

.mh_mutation <- function(Q, op, pr) {
  switch(class(op)[1L],
    PolynomialMutation = {
      lb <- op$lb %mh_or% pr$lb; ub <- op$ub %mh_or% pr$ub
      p <- op$p %mh_or% (1 / length(lb))
      if (!is.finite(p)) p <- 1e-2
      for (i in seq_len(nrow(Q))) Q[i, ] <- .mh_polynomial_mutation(Q[i, ], lb, ub, op$eta, p)
      Q
    },
    BitFlipMutation = {
      R <- matrix(stats::runif(length(Q)), nrow(Q))
      Q[] <- as.numeric(xor(Q != 0, R < op$p))
      Q
    },
    InsertRandomMutation = rbind(Q, matrix(stats::runif(op$num_mutants * ncol(Q)), op$num_mutants, ncol(Q))),
    SlightMutation = {
      N <- nrow(Q); D <- ncol(Q)
      k <- sample.int(D, N, replace = TRUE); s <- sample.int(D, N, replace = TRUE)
      for (i in seq_len(N)) {
        ki <- k[i]; si <- s[i]
        if (ki == si) next
        t <- Q[i, ki]
        if (si < ki) {
          Q[i, (si + 1L):ki] <- Q[i, si:(ki - 1L)]
          Q[i, si] <- t
        } else {
          if (ki <= si - 2L) Q[i, ki:(si - 2L)] <- Q[i, (ki + 1L):(si - 1L)]
          Q[i, si - 1L] <- t
        }
      }
      Q
    },
    stop("unknown mutation ", class(op)[1L], call. = FALSE))
}

.mh_environmental_selection <- function(pop, off, op) {
  switch(class(op)[1L],
    GenerationalReplacement = {
      if (length(pop) != length(off)) stop("population and offsprings must have the same length.", call. = FALSE)
      off
    },
    ElitistReplacement = {
      N <- length(pop)
      all <- c(pop, off)
      all[.mh_sortperm(all)][seq_len(N)]
    },
    RankAndCrowding = {
      # Metaheuristics.jl has no nondominated sorting of single objective
      # solutions either (a MethodError)
      if (!.mh_is_mo(pop[[1L]]))
        stop("RankAndCrowding needs a multiobjective population.", call. = FALSE)
      .mh_truncate_population(c(pop, off), length(pop))
    },
    stop("unknown environmental selection ", class(op)[1L], call. = FALSE))
}

.mh_register("GA",
  initialize = function(st, par, pr, info, op) {
    N <- par$initializer$N
    if (N <= 0L) stop("the initializer needs N > 0", call. = FALSE)
    if (op$iterations == 0L) op$iterations <- 500L
    if (op$f_calls_limit == 0) op$f_calls_limit <- 2 * op$iterations * N
    par$N <- N
    init <- par$initializer
    if (!is.null(st$best_sol) || inherits(init, "RandomInBounds")) return(.mh_gen_initial_state(pr, par, info, op, st))
    X <- if (inherits(init, "RandomPermutation")) {
      t(vapply(seq_len(N), function(i) as.numeric(sample.int(pr$D)), numeric(pr$D)))
    } else {
      matrix(as.numeric(stats::runif(N * pr$D) < 0.5), N, pr$D)
    }
    pop <- .mh_create_solutions(matrix(X, N, pr$D), pr, op$h_tol)
    .mh_state(.mh_get_best(pop), pop, f_calls = N)
  },
  update = function(st, par, pr, info, op) {
    pop <- st$population
    mask <- .mh_selection(pop, par$selection)
    Q <- .mh_crossover(pop[mask], par$crossover, pr)
    # the mutation's value is used, so that InsertRandomMutation adds its
    # mutants (Metaheuristics.jl drops them)
    Q <- .mh_mutation(Q, par$mutation, pr)
    off <- .mh_create_solutions(Q, pr)
    st$population <- .mh_environmental_selection(pop, off, par$environmental_selection)
    # the best of the new population (Metaheuristics.jl takes it before the
    # replacement, one generation late)
    cb <- .mh_get_best(st$population)
    if (.mh_is_better(cb, st$best_sol)) st$best_sol <- cb
  },
  stop_criteria = function(st, par, pr, info, op) invisible(NULL))

# ---------------------------------------------------------------------------
# MCCGA: machine coded compact genetic algorithm

# bits: the 32 bits of the single precision value, most significant first
.mh_mccga_bits <- function(x) {
  raw <- writeBin(as.numeric(x), raw(), size = 4L, endian = "big")
  b <- matrix(as.integer(rawToBits(raw)), 8L)
  as.integer(b[8:1, ])
}

.mh_mccga_floats <- function(bits) {
  n <- length(bits) %/% 32L
  m <- matrix(as.integer(bits)[seq_len(32L * n)], 32L)
  readBin(packBits(as.raw(m[32:1, , drop = FALSE]), "raw"), "double", n = n, size = 4L, endian = "little")
}

.mh_mccga_initialprobs <- function(lower, upper, maxsamples = 10000L) {
  p <- length(lower)
  vals <- matrix(stats::runif(maxsamples * p), maxsamples, p, byrow = TRUE)
  vals <- sweep(sweep(vals, 2L, upper - lower, "*"), 2L, lower, "+")
  raw <- writeBin(as.numeric(t(vals)), raw(), size = 4L, endian = "big")
  b <- matrix(as.integer(rawToBits(raw)), 8L)[8:1, , drop = FALSE]
  colMeans(matrix(as.integer(b), ncol = 32L * p, byrow = TRUE))
}

.mh_mccga_exploratory <- function(f, base, step) {
  x <- base
  fx <- f(x)
  for (i in seq_along(x)) {
    xl <- x; xl[i] <- xl[i] - step
    fl <- f(xl)
    if (fl < fx) { x <- xl; fx <- fl; next }
    xr <- x; xr[i] <- xr[i] + step
    fr <- f(xr)
    if (fr < fx) { x <- xr; fx <- fr }
  }
  list(x = x, fx = fx)
}

# hj: the Hooke and Jeeves pattern search of MCCGA
.mh_mccga_hj <- function(f, par, maxiter = 1000L, startstep = 5, endstep = 1e-4) {
  step <- startstep
  base <- par
  fbase <- f(base)
  iter <- 0L
  while (iter < maxiter) {
    e <- .mh_mccga_exploratory(f, base, step)
    xnew <- e$x; fnew <- e$fx
    if (fnew < fbase) {
      repeat {
        e2 <- .mh_mccga_exploratory(f, xnew + (xnew - base), step)
        if (e2$fx < fnew) {
          base <- xnew; fbase <- fnew; xnew <- e2$x; fnew <- e2$fx
        } else {
          base <- xnew; fbase <- fnew
          break
        }
      }
    } else {
      step <- step / 2
    }
    if (step < endstep) break
    iter <- iter + 1L
  }
  list(par = base, iter = iter, step = step)
}

.mh_register("MCCGA",
  initialize = function(st, par, pr, info, op) {
    if (op$iterations == 0L) op$iterations <- 100000000L
    if (op$f_calls_limit == 0) op$f_calls_limit <- as.numeric(op$iterations) * par$N + 1
    par$probvector <- .mh_mccga_initialprobs(pr$lb, pr$ub, par$maxsamples)
    .mh_gen_initial_state(pr, par, info, op, st)
  },
  update = function(st, par, pr, info, op) {
    pv <- par$probvector
    ch1 <- as.integer(stats::runif(length(pv)) < pv)
    ch2 <- as.integer(stats::runif(length(pv)) < pv)
    s1 <- .mh_create_solution(.mh_mccga_floats(ch1), pr)
    s2 <- .mh_create_solution(.mh_mccga_floats(ch2), pr)
    win <- s1; w <- ch1; l <- ch2
    if (.mh_is_better(s2, s1)) { win <- s2; w <- ch2; l <- ch1 }
    st$population[[1L + st$iteration %% par$N]] <- win
    if (.mh_is_better(win, st$best_sol)) st$best_sol <- win
    d <- w != l
    pv[d] <- pv[d] + ifelse(w[d] == 1L, par$mutation, -par$mutation)
    pv[d] <- pmin(pmax(pv[d], 0), 1)
    par$probvector <- pv
  },
  final = function(st, par, pr, info, op) {
    if (!par$use_local_search || !is.null(st$best_sol$sv)) return(invisible(NULL))
    if (st$stop && st$termination_status_code %in% c("EVALUATIONS_LIMIT", "TIME_LIMIT", "ACCURACY_LIMIT")) return(invisible(NULL))
    cost <- function(x) as.numeric(.mh_evaluate(x, pr))
    x0 <- .mh_mccga_floats(as.integer(stats::runif(length(par$probvector)) < par$probvector))
    bp <- .mh_mccga_hj(cost, x0)$par
    sol <- .mh_child(bp, cost(bp))
    st$f_calls <- pr$f_calls
    if (.mh_is_better(sol, st$best_sol)) st$best_sol <- sol
  },
  stop_criteria = function(st, par, pr, info, op) {
    if (st$stop) return(invisible(NULL))
    m <- par$mutation
    st$stop <- all(par$probvector <= m | par$probvector >= 1 - m)
    if (st$stop) st$termination_status_code <- "OTHER_LIMIT"
  })

# ---------------------------------------------------------------------------
# RDEx: the reconstructed differential evolution of Metaheuristics.jl
# (indices are kept 0 based as in the Julia module)

.mh_weighted_sample <- function(w) {
  total <- sum(w)
  if (total == 0 || !length(w)) return(sample.int(length(w), 1L) - 1L)
  r <- stats::runif(1L) * total
  cs <- 0
  for (i in seq_along(w)) {
    cs <- cs + w[i]
    if (r <= cs) return(i - 1L)
  }
  length(w) - 1L
}

.mh_rdex_weighted_mean <- function(values, weights, n) {
  sw <- .mh_seqsum(weights[seq_len(n)])
  if (sw == 0) return(1)
  w <- weights[seq_len(n)] / sw
  ssq <- .mh_seqsum(w * values[seq_len(n)] * values[seq_len(n)])
  sv <- .mh_seqsum(w * values[seq_len(n)])
  if (abs(sv) > 1e-8) ssq / sv else 1
}

.mh_rdex_trial <- function(xc, b1, b2, b3, F, Cr, lb, ub) {
  D <- length(xc)
  trial <- numeric(D)
  wc <- sample.int(D, 1L) - 1L
  pert <- stats::runif(1L) < 0.4
  acr <- 0
  for (j in seq_len(D)) {
    if (stats::runif(1L) < Cr || wc == j - 1L) {
      trial[j] <- xc[j] + F * (b1[j] - xc[j]) + F * (b2[j] - b3[j])
      if (trial[j] < lb[j] || trial[j] > ub[j]) trial[j] <- stats::runif(1L) * (ub[j] - lb[j]) + lb[j]
      acr <- acr + 1
    } else {
      trial[j] <- if (pert) stats::rcauchy(1L, xc[j], 0.1) else xc[j]
    }
  }
  list(trial = trial, cr = acr / D)
}

.mh_register("RDEx",
  initialize = function(st, par, pr, info, op) {
    D <- pr$D
    if (par$N <= 4L) par$N <- as.integer(max(4L, par$population_size * D))
    n0 <- par$N
    par$n_inds_front <- n0
    par$n_inds_front_max <- n0
    if (op$f_calls_limit == 0) op$f_calls_limit <- max(10000 * D, par$N * op$iterations)
    if (op$iterations == 0L) op$iterations <- as.integer(min(max(4 * op$f_calls_limit, 500), .Machine$integer.max))
    par$memory_size <- 5L
    par$memory_iter <- 0L
    par$success_filled <- 0L
    par$success_rate <- 0.5
    par$eb_hybrid_rate <- 0.7
    ps <- 2L * n0
    par$memory_cr <- rep(1, par$memory_size)
    par$memory_f <- rep(1, par$memory_size)
    par$temp_success_cr <- numeric(ps); par$temp_success_f <- numeric(ps)
    par$fit_delta <- numeric(ps); par$weights <- numeric(ps)
    par$fit_temp <- numeric(n0); par$eb_hybrid_flag <- numeric(n0)
    st <- .mh_gen_initial_state(pr, par, info, op, st)
    st$population <- st$population[.mh_sortperm(st$population)]
    st$best_sol <- st$population[[1L]]
    st$iteration <- 0L
    st
  },
  update = function(st, par, pr, info, op) {
    lb <- pr$lb; ub <- pr$ub
    st$population <- st$population[.mh_sortperm(st$population)]
    st$best_sol <- st$population[[1L]]
    X <- .mh_positions(st$population)
    fits <- .mh_fvals(st$population)
    nf <- length(st$population)
    par$n_inds_front <- nf
    if (length(par$fit_temp) < nf) length(par$fit_temp) <- nf
    if (length(par$eb_hybrid_flag) < nf) length(par$eb_hybrid_flag) <- nf
    ft2 <- exp(-(0:(nf - 1L)) / nf * 3)
    psize <- max(2L, as.integer(round(nf * 0.7 * exp(-par$success_rate * 7))))
    ftp <- 3 * (nf - 0:(nf - 1L))
    psize2 <- as.integer(round(nf * 0.17 * (1 - 0.5 * st$f_calls / op$f_calls_limit)))
    if (psize2 <= 1L) psize2 <- 2L
    idx <- order(fits) - 1L
    mean_f <- 0.4 + tanh(par$success_rate * 5) * 0.25
    par$success_filled <- 0L
    max_success <- length(par$temp_success_cr)
    prog <- st$f_calls / op$f_calls_limit
    for (it in seq_len(nf)) {
      if (st$f_calls >= op$f_calls_limit) break
      chosen <- sample.int(nf, 1L) - 1L
      mi <- sample.int(par$memory_size, 1L) - 1L
      mi2 <- sample.int(par$memory_size + 1L, 1L) - 1L
      reb <- stats::runif(1L)
      if (prog < 0.7) reb <- 2
      xc <- X[chosen + 1L, ]
      fc <- fits[chosen + 1L]
      if (reb * (1 - prog) < par$eb_hybrid_rate) {
        par$eb_hybrid_flag[it] <- 1
        repeat { prand <- idx[.mh_weighted_sample(ftp[seq_len(psize2)]) + 1L]; if (prand != chosen) break }
        repeat { r1 <- idx[.mh_weighted_sample(ftp) + 1L]; if (r1 != prand) break }
        repeat { r2 <- idx[.mh_weighted_sample(ftp) + 1L]; if (r2 != prand && r2 != r1) break }
        cand <- list(list(X[prand + 1L, ], fits[prand + 1L]), list(X[r1 + 1L, ], fits[r1 + 1L]), list(X[r2 + 1L, ], fits[r2 + 1L]))
        f1 <- cand[[1L]][[2L]]; f3 <- cand[[2L]][[2L]]; f4 <- cand[[3L]][[2L]]
        ord <- if (f1 <= f3 && f1 <= f4) { if (f3 <= f4) c(1L, 2L, 3L) else c(1L, 3L, 2L) } else
          if (f3 <= f1 && f3 <= f4) { if (f1 <= f4) c(2L, 1L, 3L) else c(2L, 3L, 1L) } else
          { if (f1 <= f3) c(3L, 1L, 2L) else c(3L, 2L, 1L) }
        loc <- if (mi2 < par$memory_size) par$memory_f[mi2 + 1L] else 0.9
        F <- stats::rcauchy(1L, loc, 0.1)
        while (F < 0) F <- stats::rcauchy(1L, loc, 0.1)
        if (F > 1) F <- 1
        if (prog < 0.6 && F > 0.7) F <- 0.7
        Cr <- if (mi2 < par$memory_size) {
          if (par$memory_cr[mi2 + 1L] < 0) 0 else stats::rnorm(1L, par$memory_cr[mi2 + 1L], 0.1)
        } else stats::rnorm(1L, 0.9, 0.1)
        Cr <- min(max(Cr, 0), 1)
        if (prog < 0.25) Cr <- max(Cr, 0.7)
        if (prog < 0.5) Cr <- max(Cr, 0.6)
        tr <- .mh_rdex_trial(xc, cand[[ord[1L]]][[1L]], cand[[ord[2L]]][[1L]], cand[[ord[3L]]][[1L]], F, Cr, lb, ub)
      } else {
        par$eb_hybrid_flag[it] <- 0
        repeat { prand <- idx[sample.int(psize, 1L)]; if (prand != chosen) break }
        repeat { r1 <- idx[.mh_weighted_sample(ft2) + 1L]; if (r1 != prand) break }
        repeat { r2 <- idx[sample.int(nf, 1L)]; if (r2 != prand && r2 != r1) break }
        F <- stats::rnorm(1L, mean_f, 0.02)
        while (F < 0 || F > 1) F <- stats::rnorm(1L, mean_f, 0.02)
        Cr <- min(max(stats::rnorm(1L, par$memory_cr[mi + 1L], 0.05), 0), 1)
        tr <- .mh_rdex_trial(xc, X[prand + 1L, ], X[r1 + 1L, ], X[r2 + 1L, ], F, Cr, lb, ub)
      }
      ts <- .mh_create_solution(tr$trial, pr)
      par$fit_temp[it] <- ts$f
      if (ts$f <= fc) {
        st$population[[chosen + 1L]] <- ts
        if (par$success_filled < max_success) {
          k <- par$success_filled + 1L
          par$temp_success_cr[k] <- tr$cr
          par$temp_success_f[k] <- F
          par$fit_delta[k] <- abs(fc - ts$f)
        }
        par$success_filled <- par$success_filled + 1L
      }
    }
    # update_eb_hybrid_param!
    eb <- par$eb_hybrid_flag[seq_len(nf)] == 1
    imp <- par$fit_temp[seq_len(nf)] <= fits
    gain <- fits - par$fit_temp[seq_len(nf)]
    s_eb <- sum(gain[eb & imp]); s_or <- sum(gain[!eb & imp])
    par$eb_hybrid_rate <- if (s_eb != 0 && s_or != 0) min(max(s_eb / (s_eb + s_or), 0), 1) else 0.7
    par$success_rate <- par$success_filled / nf
    nn <- as.integer(round((4 - par$n_inds_front_max) / op$f_calls_limit * st$f_calls + par$n_inds_front_max))
    nn <- min(max(nn, 4L), par$n_inds_front_max)
    st$population <- st$population[.mh_sortperm(st$population)]
    if (length(st$population) > nn) st$population <- st$population[seq_len(nn)]
    par$n_inds_front <- length(st$population)
    if (par$success_filled != 0L) {
      n <- min(par$success_filled, max_success)
      k <- par$memory_iter + 1L
      par$memory_cr[k] <- 0.5 * (.mh_rdex_weighted_mean(par$temp_success_cr, par$fit_delta, n) + par$memory_cr[k])
      par$memory_f[k] <- .mh_rdex_weighted_mean(par$temp_success_f, par$fit_delta, n)
      par$memory_iter <- (par$memory_iter + 1L) %% par$memory_size
    }
    st$best_sol <- st$population[[1L]]
    st$success_rate <- par$success_rate
  },
  final = function(st, par, pr, info, op) {
    if (length(st$population)) {
      st$population <- st$population[.mh_sortperm(st$population)]
      st$best_sol <- st$population[[1L]]
    }
  })

# ---------------------------------------------------------------------------
# GRASP, VND and VNS: the combinatorial local search methods

# the i-th neighbor of x (NULL past the last one), or a random neighbor
# when i is NULL
.mh_neighbor <- function(x, nb, i) {
  if (is.function(nb)) return(nb(x, i))
  if (inherits(nb, "TwoOptNeighborhood")) {
    k <- nb$k - 1L
    n <- length(x)
    if (is.null(i)) {
      i <- sample.int(n - k, 1L)
    } else if (!(1L <= i && i <= n - k && 1L <= k && k <= n)) {
      return(NULL)
    }
    x[i:(i + k)] <- rev(x[i:(i + k)])
    return(x)
  }
  if (is.function(nb$structure)) return(nb$structure(x, i))
  stop("define the neighborhood structure (a function(x, i))", call. = FALSE)
}

.mh_local_search_nb <- function(x, nb, ls, pr) {
  first <- inherits(ls, "FirstImproveSearch")
  best <- .mh_create_solution(x, pr)
  init <- best
  i <- 1L
  repeat {
    xn <- .mh_neighbor(x, nb, i)
    if (is.null(xn)) break
    sol <- .mh_create_solution(xn, pr)
    if (first) {
      if (.mh_is_better(sol, init)) return(sol)
    } else if (.mh_is_better(sol, best)) {
      best <- sol
    }
    i <- i + 1L
  }
  best
}

# local_search(x, ls, problem) of GRASP: a function local search, or a
# best or first improvement search over the neighborhood of `ls`
.mh_local_search <- function(x, ls, pr) {
  if (is.function(ls)) return(ls(x))
  nb <- ls$neighborhood
  if (is.null(nb)) {
    # Metaheuristics.jl fails here on box constrained spaces (its test names
    # an undefined type); the unsupported space branch is taken instead
    warning("Definition of local_search for BoxConstrainedSpace is required.", call. = FALSE)
    return(NULL)
  }
  .mh_local_search_nb(x, nb, ls, pr)
}

.mh_construct <- function(cons) {
  if (is.function(cons)) return(cons())
  if (inherits(cons, "GreedyRandomizedConstructor")) {
    cand <- cons$candidates
    S <- cand[0]
    while (length(cand)) {
      cost <- cons$compute_cost(cand, cons, cons$instance)
      cmin <- min(cost); cmax <- max(cost)
      RCL <- which(cost <= cmin + cons$alpha * (cmax - cmin))
      if (!length(RCL)) stop("RCL is empty. Try increasing alpha or check your `compute_cost` method.", call. = FALSE)
      s <- RCL[sample.int(length(RCL), 1L)]
      S <- c(S, cand[s])
      cand <- cand[-s]
    }
    return(S)
  }
  stop("Define your own randomized greedy constructor", call. = FALSE)
}

.mh_vns_initialize <- function(st, par, pr, info, op) {
  x0 <- par$initial %mh_or% .mh_rand_x(pr)
  op$f_calls_limit <- Inf
  if (op$iterations <= 0L) op$iterations <- 500L
  sol <- .mh_create_solution(as.numeric(x0), pr)
  .mh_state(sol, list(sol))
}

.mh_register("GRASP",
  initialize = .mh_vns_initialize,
  update = function(st, par, pr, info, op) {
    x <- .mh_construct(par$constructor)
    if (is.null(x) || !length(x)) {
      st$stop <- TRUE
      return(invisible(NULL))
    }
    xi <- .mh_local_search(x, par$local_search, pr)
    sol <- if (is.list(xi) && !is.null(xi$f)) xi else if (is.numeric(xi)) .mh_create_solution(xi, pr) else NULL
    if (is.null(sol)) return(invisible(NULL))
    if (.mh_is_better(sol, st$best_sol)) st$best_sol <- sol
  })

.mh_vnd_update <- function(st, par, pr, info, op) {
  improvement <- FALSE
  l <- 1L
  while (l <= length(par$neighborhood)) {
    sol <- .mh_local_search_nb(st$best_sol$x, par$neighborhood[[l]], par$local_search, pr)
    if (.mh_is_better(sol, st$best_sol)) {
      st$best_sol <- sol
      l <- 1L
      improvement <- TRUE
    } else {
      l <- l + 1L
    }
  }
  st$stop <- !improvement
}

.mh_register("VND", initialize = .mh_vns_initialize, update = .mh_vnd_update)

.mh_register("VNS",
  initialize = .mh_vns_initialize,
  update = function(st, par, pr, info, op) {
    sol <- st$population[[1L]]
    k <- 1L
    while (k <= length(par$neighborhood_shaking)) {
      xp <- .mh_neighbor(sol$x, par$neighborhood_shaking[[k]], NULL)
      vnd <- VND(initial = xp, neighborhood = par$neighborhood_local, local_search = par$local_search)
      res <- .mh_optimize(pr$f, pr$lb, pr$ub, vnd, pr$mo)
      new <- res$best_sol
      sequential <- inherits(par$neighborhood_change, "SequentialChange")
      if (sequential) {
        if (.mh_is_better(new, sol)) { k <- 1L; sol <- new } else k <- k + 1L
      } else {
        k <- k + 1L
        if (.mh_is_better(new, sol)) sol <- new
      }
      if (.mh_is_better(new, st$best_sol)) st$best_sol <- new
      pr$f_calls <- pr$f_calls + res$f_calls
    }
    st$population <- list(sol)
  })
