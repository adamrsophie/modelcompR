# The multiobjective algorithms of Metaheuristics.jl (3.6), ported to R on
# the framework of optimization-metaheuristics-base.R: NSGA2 (the base of CCMO
# and Metaheuristics$NSGA2), NSGA3, SPEA2, SMS_EMOA, MOEAD_DE and CCMO.

.mh_nsga_initialize <- function(st, par, pr, info, op) {
  if (par$p_m < 0) par$p_m <- 1 / pr$D
  if (op$iterations == 0L) op$iterations <- 500L
  if (op$f_calls_limit == 0) op$f_calls_limit <- as.numeric(op$iterations) * par$N + 1
  .mh_gen_initial_state(pr, par, info, op, st)
}

# tournament_selection by rank, then crowding
.mh_tournament <- function(P, a = -1L) {
  n <- length(P)
  if (a < 0L) a <- sample.int(n, 1L)
  b <- sample.int(n, 1L)
  while (a == b) b <- sample.int(n, 1L)
  pa <- P[[a]]; pb <- P[[b]]
  if (pa$rank < pb$rank || (pa$rank == pb$rank && pa$crowding > pb$crowding)) pa else pb
}

.mh_nsga_reproduction <- function(st, par, pr) {
  N <- par$N
  Q <- matrix(0, 2L * N, pr$D)
  for (i in seq_len(N)) {
    pa <- .mh_tournament(st$population)
    pb <- .mh_tournament(st$population)
    cc <- .mh_ga_reproduction(pa$x, pb$x, pr$lb, pr$ub, par$eta_cr, par$eta_m, par$p_cr, par$p_m)
    Q[2L * i - 1L, ] <- cc[[1L]]
    Q[2L * i, ] <- cc[[2L]]
  }
  Q
}

.mh_register("NSGA2",
  initialize = .mh_nsga_initialize,
  update = function(st, par, pr, info, op) {
    Q <- .mh_nsga_reproduction(st, par, pr)
    st$population <- .mh_truncate_population(c(st$population, .mh_create_solutions(Q, pr)), par$N)
  })

# ---------------------------------------------------------------------------
# SPEA2

# pairwise_distances: the Euclidean distances of the objective vectors as
# Distances.jl computes them (from the Gram matrix), the diagonal set to Inf
.mh_pairwise_distances <- function(pop) {
  Fm <- .mh_fvals(pop)
  if (!is.matrix(Fm)) Fm <- matrix(Fm, ncol = 1L)
  n <- nrow(Fm)
  G <- matrix(0, n, n)
  sa2 <- numeric(n)
  for (m in seq_len(ncol(Fm))) {
    G <- G + outer(Fm[, m], Fm[, m])
    sa2 <- sa2 + Fm[, m] * Fm[, m]
  }
  D <- sqrt(pmax(outer(sa2, sa2, "+") - 2 * G, 0))
  diag(D) <- Inf
  D
}

.mh_spea2_fitness <- function(pop, Dist = .mh_pairwise_distances(pop)) {
  N <- length(pop)
  dom <- .mh_dominance(.mh_fvals(pop), .mh_svs(pop))
  S <- rowSums(dom)
  R <- vapply(seq_len(N), function(i) sum(S[dom[, i]]), 0)
  srt <- t(apply(Dist, 1L, sort))
  if (N == 1L) srt <- matrix(srt, 1L)
  R + 1 / (srt[, floor(sqrt(N))] + 2)
}

.mh_spea2_truncation <- function(K, distance) {
  del <- rep(FALSE, nrow(distance))
  while (sum(del) < K) {
    remain <- which(!del)
    sub <- distance[remain, remain, drop = FALSE]
    nn <- which.min(sub)
    del[remain[(nn - 1L) %% length(remain) + 1L]] <- TRUE
  }
  del
}

.mh_spea2_selection <- function(pop, N) {
  Dist <- .mh_pairwise_distances(pop)
  fitness <- .mh_spea2_fitness(pop, Dist)
  nxt <- fitness < 1
  K <- sum(nxt)
  if (K < N) {
    nxt[order(fitness)[seq_len(N)]] <- TRUE
  } else if (K > N) {
    del <- .mh_spea2_truncation(K - N, Dist[nxt, nxt, drop = FALSE])
    nxt[which(nxt)[del]] <- FALSE
  }
  list(pop = pop[nxt], fitness = fitness[nxt], keep = nxt)
}

.mh_register("SPEA2",
  initialize = .mh_nsga_initialize,
  update = function(st, par, pr, info, op) {
    if (!length(par$fitness)) par$fitness <- .mh_spea2_fitness(st$population)
    Q <- .mh_nsga_reproduction(st, par, pr)
    r <- .mh_spea2_selection(c(st$population, .mh_create_solutions(Q, pr)), par$N)
    st$population <- r$pop
    par$fitness <- r$fitness
  })

# ---------------------------------------------------------------------------
# NSGA3

.mh_nsga3_normalize <- function(pop) {
  F <- t(.mh_fvals(pop))
  M <- nrow(F); n <- ncol(F)
  rk <- vapply(pop, function(s) s$rank, 0L)
  first <- rk == min(rk)
  ideal <- apply(F, 1L, min)
  worst <- apply(F, 1L, max)
  ext <- vapply(seq_len(M), function(i) {
    w <- (seq_len(M) == i) + 1e-6
    which.min(apply((F - ideal) / w, 2L, max))
  }, 0L)
  S <- t(F[, ext, drop = FALSE] - ideal)
  nadir <- .mh_nsga3_intercepts_nadir(S, ideal, worst)
  if (is.null(nadir)) {
    nadir <- apply(F[, first, drop = FALSE], 1L, max)
    bad <- !(nadir - ideal > 1e-6 * (worst - ideal))
    nadir[bad] <- worst[bad]
  }
  b <- nadir - ideal
  b[b < .Machine$double.eps] <- .Machine$double.eps
  (F - ideal) / b
}

.mh_nsga3_intercepts_nadir <- function(S, ideal, worst) {
  M <- length(ideal)
  if (det(S) == 0) return(NULL)
  hp <- tryCatch(solve(S, rep(1, M)), error = function(e) NULL)
  if (is.null(hp) || !all(is.finite(hp))) return(NULL)
  if (any(abs(drop(S %*% hp) - 1) > 1e-5)) return(NULL)
  ic <- 1 / hp
  if (!all(is.finite(ic) & ic > 1e-6 * (worst - ideal))) return(NULL)
  ideal + ic
}

# associate!: the closest reference line of each solution (the last one on
# ties), the squared distance to it and the niche counts of the first l - 1
# solutions; the sums run over the objectives in order, as in Julia
.mh_nsga3_associate <- function(F, refs, l) {
  M <- nrow(F); n <- ncol(F); H <- length(refs)
  W <- matrix(unlist(refs), ncol = H)
  ww <- vapply(refs, function(w) .mh_seqsum(w * w), 0)
  nich <- integer(n); dist <- rep(Inf, n); freq <- integer(H)
  for (i in seq_len(n)) {
    ws <- numeric(H)
    for (m in seq_len(M)) ws <- ws + W[m, ] * F[m, i]
    cfac <- ws / ww
    d <- numeric(H)
    for (m in seq_len(M)) {
      r <- F[m, i] - cfac * W[m, ]
      d <- d + r * r
    }
    jb <- max(which(d == min(d)))
    dist[i] <- d[jb]
    nich[i] <- jb
    if (i < l) freq[jb] <- freq[jb] + 1L
  }
  list(nich = nich, dist = dist, freq = freq)
}

.mh_nsga3_niching <- function(pop, refs, N, l) {
  n <- length(pop)
  if (n == N) return(pop)
  F <- .mh_nsga3_normalize(pop)
  H <- length(refs)
  a <- .mh_nsga3_associate(F, refs, l)
  nich <- a$nich; dist <- a$dist; freq <- a$freq
  members <- lapply(seq_len(H), function(j) { m <- which(nich == j); m[m >= l] })
  K <- N - l + 1L
  selected <- integer(K)
  for (t in seq_len(K)) {
    avail <- which(lengths(members) > 0L)
    mini <- min(freq[avail])
    cand <- avail[freq[avail] == mini]
    jh <- cand[sample.int(length(cand), 1L)]
    mem <- members[[jh]]
    p <- if (mini == 0L) which.min(dist[mem]) else sample.int(length(mem), 1L)
    selected[t] <- mem[p]
    members[[jh]] <- mem[-p]
    freq[jh] <- freq[jh] + 1L
  }
  c(pop[seq_len(l - 1L)], pop[selected])
}

.mh_nsga3_truncate <- function(pop, refs, N) {
  pop <- .mh_fast_nds(pop)
  if (length(pop) <= N) return(pop)
  k <- pop[[N]]$rank
  nk <- N
  while (nk < length(pop) && pop[[nk + 1L]]$rank == k) nk <- nk + 1L
  pop <- pop[seq_len(nk)]
  l <- N
  while (l > 1L && pop[[l - 1L]]$rank == k) l <- l - 1L
  .mh_nsga3_niching(pop, refs, N, l)
}

.mh_register("NSGA3",
  initialize = function(st, par, pr, info, op) {
    if (op$iterations == 0L) op$iterations <- 500L
    if (op$f_calls_limit == 0) op$f_calls_limit <- as.numeric(op$iterations) * par$N + 1
    st <- .mh_gen_initial_state(pr, par, info, op, st)
    if (!length(par$reference_points)) par$reference_points <- .mh_gen_weights(length(st$population[[1L]]$f), par$partitions)
    st
  },
  update = function(st, par, pr, info, op) {
    N <- par$N; D <- pr$D
    p_m <- if (par$p_m < 0) 1 / D else par$p_m
    pop <- st$population
    n <- length(pop)
    I <- sample.int(n)
    Q <- matrix(0, N, D)
    for (i in seq_len(ceiling(N / 2))) {
      pa <- pop[[I[(2L * i - 2L) %% n + 1L]]]
      pb <- pop[[I[(2L * i - 1L) %% n + 1L]]]
      cc <- .mh_ga_reproduction(pa$x, pb$x, pr$lb, pr$ub, par$eta_cr, par$eta_m, par$p_cr, p_m)
      Q[2L * i - 1L, ] <- cc[[1L]]
      if (2L * i <= N) Q[2L * i, ] <- cc[[2L]]
    }
    st$population <- .mh_nsga3_truncate(c(pop, .mh_create_solutions(Q, pr)), par$reference_points, N)
  })

# ---------------------------------------------------------------------------
# SMS_EMOA

.mh_sms_update_population <- function(pop, off, n_samples) {
  pop <- .mh_fast_nds(c(pop, list(off)))
  n <- length(pop)
  r <- pop[[n]]$rank
  i <- n
  while (i > 1L && pop[[i - 1L]]$rank == r) i <- i - 1L
  last <- i:n
  for (k in seq_len(n)) pop[[k]]$crowding <- Inf
  dS <- .mh_compute_contribution(pop[last], n_samples)
  for (k in seq_along(last)) pop[[last[k]]]$crowding <- dS[k]
  cmin <- Inf; worst <- last[length(last)]; nties <- 0L
  for (k in last) {
    cc <- pop[[k]]$crowding
    if (cc < cmin) {
      cmin <- cc; worst <- k; nties <- 1L
    } else if (cc == cmin) {
      nties <- nties + 1L
      if (sample.int(nties, 1L) == 1L) worst <- k
    }
  }
  pop[-worst]
}

.mh_register("SMS_EMOA",
  initialize = function(st, par, pr, info, op) {
    if (op$iterations == 0L) op$iterations <- 500L
    if (op$f_calls_limit == 0) op$f_calls_limit <- as.numeric(op$iterations) * par$N + 1
    .mh_gen_initial_state(pr, par, info, op, st)
  },
  update = function(st, par, pr, info, op) {
    N <- par$N
    p_m <- if (par$p_m < 0) 1 / pr$D else par$p_m
    I <- sample.int(N); J <- sample.int(N)
    if (op$parallel_evaluation) Q <- matrix(0, N, pr$D)
    for (i in seq_len(N)) {
      pa <- .mh_tournament(st$population, I[i])
      pb <- .mh_tournament(st$population, J[i])
      cc <- .mh_sbx(pa$x, pb$x, pr$lb, pr$ub, par$eta_cr, par$p_cr)[[2L]]
      cc <- .mh_reset_to_violated_bounds(.mh_polynomial_mutation(cc, pr$lb, pr$ub, par$eta_m, p_m), pr$lb, pr$ub)
      if (op$parallel_evaluation) { Q[i, ] <- cc; next }
      st$population <- .mh_sms_update_population(st$population, .mh_create_solution(cc, pr), par$n_samples)
    }
    if (op$parallel_evaluation) {
      pop <- st$population
      for (s in .mh_create_solutions(Q, pr)) pop <- .mh_sms_update_population(pop, s, par$n_samples)
      st$population <- pop
    }
  })

# ---------------------------------------------------------------------------
# MOEAD_DE

.mh_moead_neighbors <- function(lambda, T) {
  N <- length(lambda)
  L <- do.call(rbind, lambda)
  lapply(seq_len(N), function(i) {
    d <- sqrt(colSums((t(L) - L[i, ])^2))
    d[i] <- -1
    order(d)[2:(T + 1L)]
  })
}

.mh_moead_g <- function(fx, lam, z) max(lam * abs(fx - z))

.mh_moead_better <- function(g1, g2, cv1, cv2, tau, s1, s2) {
  if (cv1 == 0 && cv2 == 0) return(g1 <= g2)
  gte <- function(gx, V) if (V < tau) gx + s1 * V^2 else gx + s1 * tau^2 + s2 * (V - tau)
  gte(g1, cv1) <= gte(g2, cv2)
}

.mh_moead_reproduction <- function(i, P, pop, par, pr) {
  pick <- function() P[sample.int(length(P), 1L)]
  r2 <- pick(); while (r2 == i) r2 <- pick()
  r3 <- pick(); while (r3 == i || r3 == r2) r3 <- pick()
  a <- pop[[i]]$x; b <- pop[[r2]]$x; cc <- pop[[r3]]$x
  p_m <- if (par$p_m < 0) 1 / pr$D else par$p_m
  e <- 1 / (par$eta + 1)
  v <- numeric(length(a))
  for (j in seq_along(a)) {
    v[j] <- if (stats::runif(1L) < par$CR) a[j] + par$F * (b[j] - cc[j]) else a[j]
    if (stats::runif(1L) < p_m) {
      r <- stats::runif(1L)
      sk <- if (r < 0.5) (2 * r)^e - 1 else 1 - (2 - 2 * r)^e
      v[j] <- v[j] + sk * (pr$ub[j] - pr$lb[j])
    }
  }
  v
}

.mh_register("MOEAD_DE",
  initialize = function(st, par, pr, info, op) {
    if (op$iterations == 0L) op$iterations <- 500L
    if (op$f_calls_limit == 0) op$f_calls_limit <- as.numeric(op$iterations) * par$N + 1
    st <- .mh_gen_initial_state(pr, par, info, op, st)
    n <- length(st$population)
    if (n != length(par$lambda))
      stop("MOEAD_DE: the population size (", n, ") must equal the number of weight vectors (", length(par$lambda), ").",
           call. = FALSE)
    par$nobjectives <- length(st$population[[1L]]$f)
    for (s in st$population) par$z <- pmin(par$z, s$f)
    st
  },
  update = function(st, par, pr, info, op) {
    N <- par$N
    pop <- st$population
    Fpop <- do.call(cbind, lapply(pop, function(s) s$f))
    CV <- .mh_svs(pop)
    t0 <- .mh_now()
    for (i in seq_len(N)) {
      P <- if (stats::runif(1L) < par$delta) par$B[[i]] else seq_len(N)
      v <- .mh_replace_with_random_in_bounds(.mh_moead_reproduction(i, P, pop, par, pr), pr$lb, pr$ub)
      h <- .mh_create_solution(v, pr)
      par$z <- pmin(par$z, h$f)
      cvh <- .mh_sv(h)
      tau <- min(CV[P]) + 0.3 * (max(CV[P]) - min(CV[P]))
      par$tau <- tau
      c <- 0L
      P <- P[sample.int(length(P))]
      while (c < par$n_r && length(P)) {
        j <- P[length(P)]
        P <- P[-length(P)]
        g1 <- .mh_moead_g(h$f, par$lambda[[j]], par$z)
        g2 <- .mh_moead_g(Fpop[, j], par$lambda[[j]], par$z)
        if (.mh_moead_better(g1, g2, cvh, CV[j], tau, par$s1, par$s2)) {
          pop[[j]] <- h
          Fpop[, j] <- h$f
          CV[j] <- cvh
          c <- c + 1L
        }
      }
      st$f_calls <- pr$f_calls
      if (.mh_call_limit_stop(st, info, op)) { st$stop <- TRUE; break }
      if (is.finite(op$time_limit) && st$overall_time + (.mh_now() - t0) >= op$time_limit) {
        st$termination_status_code <- "TIME_LIMIT"
        st$stop <- TRUE
        break
      }
    }
    st$population <- pop
  })

# ---------------------------------------------------------------------------
# CCMO: coevolutionary constrained multiobjective optimization, with the
# NSGA2 operators and SPEA2 environmental selection

.mh_ccmo_selection <- function(pop, N, constraints = TRUE) {
  if (constraints) return(.mh_spea2_selection(pop, N))
  # the helper population ignores the constraints
  r <- .mh_spea2_selection(lapply(pop, function(s) { s$sv <- 0; s }), N)
  r$pop <- pop[r$keep]
  r
}

.mh_register("CCMO",
  initialize = function(st, par, pr, info, op) {
    if (!is.environment(par$base)) par$base <- list2env(par$base, parent = emptyenv())
    base <- par$base
    if (op$iterations == 0L) op$iterations <- 500L
    if (op$f_calls_limit == 0) op$f_calls_limit <- op$iterations * base$N + 1
    st <- .mh_gen_initial_state(pr, base, info, op, st)
    helper <- .mh_generate_population(base$N, pr, op$h_tol)
    r <- .mh_ccmo_selection(st$population, base$N)
    st$population <- r$pop
    par$fitness <- r$fitness
    r <- .mh_ccmo_selection(helper, base$N, FALSE)
    par$phelper <- r$pop
    par$fitness_helper <- r$fitness
    st
  },
  update = function(st, par, pr, info, op) {
    base <- par$base
    if (base$p_m < 0) base$p_m <- 1 / pr$D
    sel <- .mh_op("TournamentSelection", K = 2L, N = base$N %/% 2L)
    tour <- function(fitness, n) vapply(seq_len(sel$N), function(j) {
      best <- sample.int(n, 1L)
      i <- sample.int(n, 1L)
      if (fitness[best] > fitness[i]) best <- i
      best
    }, 0L)
    cx <- .mh_op("SBX", eta = base$eta_cr, p = base$p_cr)
    mu <- .mh_op("PolynomialMutation", eta = base$eta_m, p = base$p_m)
    pa <- tour(par$fitness, length(st$population))
    Qa <- .mh_mutation(.mh_crossover(st$population[pa], cx, pr), mu, pr)
    pb <- tour(par$fitness_helper, length(par$phelper))
    # the helper's offspring come from the helper population (Metaheuristics.jl
    # indexes the main population with the helper's selection)
    Qb <- .mh_mutation(.mh_crossover(par$phelper[pb], cx, pr), mu, pr)
    off <- c(.mh_create_solutions(Qa, pr), .mh_create_solutions(Qb, pr))
    r <- .mh_ccmo_selection(c(st$population, off), base$N)
    st$population <- r$pop
    par$fitness <- r$fitness
    r <- .mh_ccmo_selection(c(par$phelper, off), base$N, FALSE)
    par$phelper <- r$pop
    par$fitness_helper <- r$fitness
  })
