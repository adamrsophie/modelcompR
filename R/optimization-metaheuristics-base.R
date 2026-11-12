# The OptimizationMetaheuristics wrapper of Optimization.jl with the
# Metaheuristics.jl (3.6) framework, ported to R: solutions, problems,
# options, termination criteria, the optimize loop and the shared operators.
# The algorithms are in optimization-metaheuristics-algorithms.R (single
# objective and combinatorial) and optimization-metaheuristics-mo.R
# (multiobjective).
#
# Julia's mutable structs (the state, the parameters and the options of a
# run) are environments here; a solution is a list with the fields `x` and
# `f` (`xf_solution`), plus `g`, `h`, `sv` (the sum of violations) and
# `feasible` (`xfgh_solution`), plus `rank` and `crowding` when `f` is the
# vector of the objectives (`xFgh_solution`).

`%mh_or%` <- function(a, b) if (is.null(a)) b else a

.mh_now <- function() proc.time()[["elapsed"]]

# ---------------------------------------------------------------------------
# Solutions

# violationsSum: the mean positive part of g plus the mean of |h| over the
# equality constraints violated by more than eps.
.mh_violations_sum <- function(g, h, eps = 0) {
  if (!length(g) && !length(h)) return(0)
  if (!length(g)) g <- 0
  if (!length(h)) h <- 0
  ok <- !is.na(h) & (h == 0 | (is.finite(h) & abs(h) <= eps))
  .mh_seqsum(pmax(0, g)) / length(g) + .mh_seqsum(abs(h[!ok])) / length(h)
}

# a sum in double precision from left to right, as Julia's loops (R's sum
# accumulates in extended precision)
.mh_seqsum <- function(v) {
  s <- 0
  for (x in v) s <- s + x
  s
}

.mh_xFgh <- function(x, f, g, h, eps = 0) {
  sv <- .mh_violations_sum(g, h, eps)
  list(x = x, f = f, g = g, h = h, rank = 0L, crowding = 0, sv = sv, feasible = sv == 0)
}

# create_child: the solution of the value `r` of the objective at `x`, a
# number, or a list (f, g, h) of the objective and the inequality and
# equality constraints.
.mh_child <- function(x, r, eps = 0, mo = FALSE) {
  if (is.list(r)) {
    f <- as.numeric(r[[1L]])
    g <- as.numeric(r[[2L]] %mh_or% numeric(0))
    h <- as.numeric(r[[3L]] %mh_or% numeric(0))
    if (mo || length(f) > 1L) return(.mh_xFgh(x, f, g, h, eps))
    if (!length(g) && !length(h)) return(list(x = x, f = f))
    sv <- .mh_violations_sum(g, h, eps)
    return(list(x = x, f = f, g = g, h = h, sv = sv, feasible = sv == 0))
  }
  r <- as.numeric(r)
  if (mo || length(r) > 1L) return(.mh_xFgh(x, r, numeric(0), numeric(0), eps))
  list(x = x, f = r)
}

.mh_is_mo <- function(s) !is.null(s$rank)
.mh_sv <- function(s) s$sv %mh_or% 0
.mh_feasible <- function(s) if (is.null(s$feasible)) TRUE else s$feasible

# compare of two objective vectors: 0 equal, 1 a dominates b, 2 b dominates
# a, 3 incomparable
.mh_compare_vec <- function(a, b) {
  k <- length(a)
  ne <- which(is.na(a == b) | a != b)
  if (!length(ne)) return(0L)
  i <- ne[1L]
  rest <- seq_len(k - i) + i
  if (isTRUE(a[i] < b[i])) {
    if (any(b[rest] < a[rest], na.rm = TRUE)) 3L else 1L
  } else {
    if (any(a[rest] < b[rest], na.rm = TRUE)) 3L else 2L
  }
}

.mh_compare <- function(a, b) {
  sa <- .mh_sv(a); sb <- .mh_sv(b)
  if (isTRUE(sa < sb)) return(1L)
  if (isTRUE(sb < sa)) return(2L)
  .mh_compare_vec(a$f, b$f)
}

# is_better: the feasibility rules (violations first, then the objective,
# or Pareto dominance for several objectives)
.mh_is_better <- function(a, b) {
  sa <- a$sv
  if (is.null(sa)) return(isTRUE(a$f < b$f))
  sb <- b$sv %mh_or% 0
  if (isTRUE(sa < sb)) return(TRUE)
  if (isTRUE(sb < sa)) return(FALSE)
  if (is.null(a$rank)) isTRUE(a$f < b$f) else .mh_compare_vec(a$f, b$f) == 1L
}

.mh_get_best <- function(pop) {
  if (!length(pop)) return(NULL)
  best <- pop[[1L]]
  for (i in seq_along(pop)[-1L]) if (.mh_is_better(pop[[i]], best)) best <- pop[[i]]
  best
}

.mh_argbest <- function(pop) {
  if (!length(pop)) return(-1L)
  ib <- 1L
  for (i in seq_along(pop)[-1L]) if (.mh_is_better(pop[[i]], pop[[ib]])) ib <- i
  ib
}

.mh_argworst <- function(pop) {
  if (!length(pop)) return(-1L)
  iw <- 1L
  for (i in seq_along(pop)[-1L]) if (.mh_is_better(pop[[iw]], pop[[i]])) iw <- i
  iw
}

# sortperm(population, lt = is_better), a stable merge sort as Julia's
.mh_sortperm <- function(pop, lt = .mh_is_better) {
  idx <- seq_along(pop)
  msort <- function(v) {
    n <- length(v)
    if (n <= 1L) return(v)
    m <- n %/% 2L
    a <- msort(v[seq_len(m)]); b <- msort(v[(m + 1L):n])
    out <- integer(n); i <- 1L; j <- 1L; k <- 1L
    while (i <= length(a) && j <= length(b)) {
      if (lt(pop[[b[j]]], pop[[a[i]]])) { out[k] <- b[j]; j <- j + 1L } else { out[k] <- a[i]; i <- i + 1L }
      k <- k + 1L
    }
    while (i <= length(a)) { out[k] <- a[i]; i <- i + 1L; k <- k + 1L }
    while (j <= length(b)) { out[k] <- b[j]; j <- j + 1L; k <- k + 1L }
    out
  }
  msort(idx)
}

.mh_positions <- function(pop) {
  if (!length(pop)) return(matrix(0, 0, 0))
  do.call(rbind, lapply(pop, function(s) s$x))
}

.mh_fvals <- function(pop) {
  if (!length(pop)) return(numeric(0))
  if (.mh_is_mo(pop[[1L]])) return(do.call(rbind, lapply(pop, function(s) s$f)))
  vapply(pop, function(s) s$f, 0)
}

.mh_svs <- function(pop) vapply(pop, .mh_sv, 0)

# ---------------------------------------------------------------------------
# Problems (a box constrained search space) and evaluation

.mh_problem <- function(f, lb, ub, mo = FALSE) {
  pr <- new.env(parent = emptyenv())
  pr$f <- f
  pr$lb <- as.numeric(lb)
  pr$ub <- as.numeric(ub)
  pr$D <- length(lb)
  pr$Delta <- pr$ub - pr$lb
  pr$f_calls <- 0L
  pr$mo <- mo
  pr
}

.mh_create_solution <- function(x, pr, eps = 0) {
  pr$f_calls <- pr$f_calls + 1L
  .mh_child(x, pr$f(x), eps, pr$mo)
}

.mh_create_solutions <- function(X, pr, eps = 0)
  lapply(seq_len(nrow(X)), function(i) .mh_create_solution(X[i, ], pr, eps))

.mh_evaluate <- function(x, pr) {
  pr$f_calls <- pr$f_calls + 1L
  pr$f(x)
}

.mh_rand_x <- function(pr) pr$lb + pr$Delta * stats::runif(pr$D)

.mh_generate_population <- function(N, pr, eps = 0)
  lapply(seq_len(N), function(i) .mh_create_solution(.mh_rand_x(pr), pr, eps))

# ---------------------------------------------------------------------------
# Repair of points outside the bounds

.mh_reset_to_violated_bounds <- function(x, lb, ub) {
  lo <- which(lb > x); x[lo] <- lb[lo]
  hi <- which(x > ub); x[hi] <- ub[hi]
  x
}

.mh_replace_with_random_in_bounds <- function(x, lb, ub) {
  for (i in seq_along(x)) if (!isTRUE(lb[i] <= x[i] && x[i] <= ub[i])) x[i] <- lb[i] + stats::runif(1L) * (ub[i] - lb[i])
  x
}

.mh_evo_boundary_repairer <- function(x, xbest, lb, ub) {
  for (i in seq_along(x)) {
    if (isTRUE(lb[i] > x[i])) {
      a <- stats::runif(1L)
      x[i] <- a * lb[i] + (1 - a) * xbest[i]
    } else if (isTRUE(x[i] > ub[i])) {
      b <- stats::runif(1L)
      x[i] <- b * ub[i] + (1 - b) * xbest[i]
    }
  }
  x
}

.mh_reflect_far <- function(xi, l, u) {
  rho <- u - l
  if (is.nan(xi) || rho == 0) return(l)
  if (is.infinite(xi)) return(if (xi < l) l else u)
  t <- (xi - l) %% (2 * rho)
  min(max(if (t <= rho) l + t else l + (2 * rho - t), l), u)
}

.mh_reflected_back_to_bound <- function(x, lb, ub) {
  for (i in seq_along(x)) {
    l <- lb[i]; u <- ub[i]; xi <- x[i]
    if (isTRUE(l <= xi && xi <= u)) next
    y <- if (isTRUE(xi < l)) 2 * l - xi else 2 * u - xi
    x[i] <- if (isTRUE(l <= y && y <= u)) y else .mh_reflect_far(xi, l, u)
  }
  x
}

.mh_wrap_to_bounds <- function(x, lb, ub) {
  for (i in seq_along(x)) {
    l <- lb[i]; u <- ub[i]
    if (!isTRUE(l <= x[i] && x[i] <= u)) {
      rho <- u - l
      # Julia's % is the truncated remainder
      rem <- function(a, b) a - b * trunc(a / b)
      x[i] <- if (x[i] < l) u - rem(l - x[i], rho) else l + rem(x[i] - u, rho)
    }
  }
  x
}

# ---------------------------------------------------------------------------
# Variation operators

.mh_sbx_bq <- function(beta, eta, R) {
  alpha <- 2 - beta^(-eta - 1)
  s <- 1 / (eta + 1)
  if (isTRUE(R <= 1 / alpha)) (R * alpha)^s else (1 / (2 - R * alpha))^s
}

# SBX_crossover! with the random numbers given: u (crossover of each
# coordinate), R (spread) and swap (exchange of the children)
.mh_sbx_core <- function(c1, c2, lb, ub, eta, p, u, R, swap) {
  eps <- .Machine$double.eps
  for (j in seq_along(lb)) {
    a <- c1[j]; b <- c2[j]
    if (isTRUE(u[j] > p) || isTRUE(abs(b - a) <= eps)) next
    y1 <- min(a, b); y2 <- max(a, b)
    Dl <- max(eps, y2 - y1)
    bq <- .mh_sbx_bq(1 + (2 * (y1 - lb[j]) / Dl), eta, R[j])
    o1 <- 0.5 * (y1 + y2 - bq * Dl)
    bq <- .mh_sbx_bq(1 + (2 * (ub[j] - y2) / Dl), eta, R[j])
    o2 <- 0.5 * (y1 + y2 + bq * Dl)
    if (swap[j]) { t <- o1; o1 <- o2; o2 <- t }
    c1[j] <- o1; c2[j] <- o2
  }
  list(.mh_reset_to_violated_bounds(c1, lb, ub), .mh_reset_to_violated_bounds(c2, lb, ub))
}

.mh_sbx <- function(c1, c2, lb, ub, eta = 15, p = 0.9) {
  D <- length(lb)
  u <- stats::runif(D); R <- stats::runif(D); swap <- stats::runif(D) < 0.5
  .mh_sbx_core(c1, c2, lb, ub, eta, p, u, R, swap)
}

# polynomial_mutation! with the random numbers given: u (which coordinates
# mutate) and R (one number per mutated coordinate)
.mh_polynomial_core <- function(v, lb, ub, eta, prob, u, R) {
  s <- eta + 1
  mp <- 1 / (eta + 1)
  k <- 0L
  for (j in seq_along(v)) {
    if (!isTRUE(u[j] < prob)) next
    k <- k + 1L
    x <- v[j]
    Dl <- ub[j] - lb[j]
    d1 <- (x - lb[j]) / Dl
    d2 <- (ub[j] - x) / Dl
    r <- R[k]
    dq <- if (isTRUE(r < 0.5)) (2 * r + (1 - 2 * r) * (1 - d1)^s)^mp - 1 else 1 - (2 * (1 - r) + 2 * (r - 0.5) * (1 - d2)^s)^mp
    v[j] <- x + dq * Dl
  }
  v
}

.mh_polynomial_mutation <- function(v, lb, ub, eta = 20, prob = 1 / length(v)) {
  u <- stats::runif(length(v))
  R <- stats::runif(sum(u < prob))
  .mh_polynomial_core(v, lb, ub, eta, prob, u, R)
}

# GA_reproduction of the NSGA family: SBX, polynomial mutation and repair
.mh_ga_reproduction <- function(pa, pb, lb, ub, eta_cr = 20, eta_m = 15, p_cr = 0.9, p_m = 0.1) {
  cc <- .mh_sbx(pa, pb, lb, ub, eta_cr, p_cr)
  c1 <- .mh_reset_to_violated_bounds(.mh_polynomial_mutation(cc[[1L]], lb, ub, eta_m, p_m), lb, ub)
  c2 <- .mh_reset_to_violated_bounds(.mh_polynomial_mutation(cc[[2L]], lb, ub, eta_m, p_m), lb, ub)
  list(c1, c2)
}

# DE_mutation! (base vector, differences and strategy)
.mh_rand_distinct <- function(N, excluded) {
  if (length(unique(excluded)) >= N) return(sample.int(N, 1L))
  r <- sample.int(N, 1L)
  while (r %in% excluded) r <- sample.int(N, 1L)
  r
}

.mh_de_mutation <- function(pop, F = 1, strategy = "rand1", best_ind = 0L) {
  N <- length(pop)
  i <- sample.int(N, 1L)
  r1 <- .mh_rand_distinct(N, i)
  r2 <- .mh_rand_distinct(N, c(i, r1))
  x <- pop[[i]]$x; a <- pop[[r1]]$x; b <- pop[[r2]]$x
  if (strategy == "rand1") return(x + F * (a - b))
  r3 <- .mh_rand_distinct(N, c(i, r1, r2))
  cc <- pop[[r3]]$x
  if (strategy == "rand2") {
    r4 <- .mh_rand_distinct(N, c(i, r1, r2, r3))
    r5 <- .mh_rand_distinct(N, c(i, r1, r2, r3, r4))
    return(pop[[r5]]$x + F * (a - b + cc - pop[[r4]]$x))
  }
  if (best_ind == 0L) best_ind <- .mh_argbest(pop)
  xb <- pop[[best_ind]]$x
  switch(strategy,
         best1 = xb + F * (b - cc),
         randToBest1 = x + F * (xb - x + a - b),
         best2 = {
           r4 <- .mh_rand_distinct(N, c(i, r1, r2, r3, best_ind))
           xb + F * (a - b + cc - pop[[r4]]$x)
         },
         stop("Unknown strategy ", strategy, call. = FALSE))
}

# DE_crossover!: binomial crossover of the mutant u with the target x
.mh_de_crossover <- function(u, x, CR) {
  D <- length(x)
  jr <- sample.int(D, 1L)
  for (j in seq_len(D)) {
    if (stats::runif(1L) < CR || j == jr) next
    u[j] <- x[j]
  }
  u
}

# ---------------------------------------------------------------------------
# Nondominated sorting, crowding distance and the hypervolume contributions

# d[i, j] is TRUE when solution i dominates solution j (compare == 1):
# smaller violation first, then Pareto dominance of the objectives
.mh_dominance <- function(Fm, sv) {
  np <- nrow(Fm)
  le <- matrix(TRUE, np, np); lt <- matrix(FALSE, np, np)
  for (k in seq_len(ncol(Fm))) {
    le <- le & outer(Fm[, k], Fm[, k], "<=")
    lt <- lt | outer(Fm[, k], Fm[, k], "<")
  }
  d <- le & lt
  d[is.na(d)] <- FALSE
  less <- outer(sv, sv, "<")
  same <- outer(sv, sv, "==")
  less | (same & d)
}

.mh_non_dominated_sort <- function(pop) {
  if (!length(pop)) return(integer(0))
  .opt_ranks(.mh_dominance(.mh_fvals(pop), .mh_svs(pop)))
}

# fast_non_dominated_sort!: the ranks stored and the population ordered by
# rank (a stable sort)
.mh_fast_nds <- function(pop) {
  r <- .mh_non_dominated_sort(pop)
  for (i in seq_along(pop)) pop[[i]]$rank <- r[i]
  pop[order(r)]
}

# The crowding distance of a front. Metaheuristics.jl accumulates the
# distances of more than two objectives by position while it re-sorts the
# front for each objective, which mixes up the solutions; here every
# solution gets its own sum. The front is returned in the order that
# compute_crowding_distance leaves it (stably sorted by each objective in
# turn).
.mh_crowding <- function(pop) {
  n <- length(pop)
  Fm <- .mh_fvals(pop)
  M <- ncol(Fm)
  o1 <- order(Fm[, 1L])
  cd <- numeric(n)
  if (M == 2L) {
    Fs <- Fm[o1, , drop = FALSE]
    if (Fs[1L, 1L] == Fs[n, 1L]) return(list(pop = pop[o1], cd = cd))
    cd[c(1L, n)] <- Inf
    w1 <- Fs[n, 1L] - Fs[1L, 1L]
    w2 <- Fs[1L, 2L] - Fs[n, 2L]
    i <- 2:(n - 1L)
    cd[i] <- (Fs[i + 1L, 1L] - Fs[i - 1L, 1L]) / w1 + (Fs[i - 1L, 2L] - Fs[i + 1L, 2L]) / w2
    return(list(pop = pop[o1], cd = cd))
  }
  o <- seq_len(n)
  for (j in seq_len(M)) {
    o <- o[order(Fm[o, j])]
    cd[o[c(1L, n)]] <- Inf
    if (Fm[o[1L], j] != Fm[o[n], j]) {
      i <- 2:(n - 1L)
      cd[o[i]] <- cd[o[i]] + (Fm[o[i + 1L], j] - Fm[o[i - 1L], j]) / (Fm[o[n], j] - Fm[o[1L], j])
    }
  }
  list(pop = pop[o], cd = cd[o])
}

.mh_update_crowding <- function(pop) {
  if (length(pop) < 3L) return(pop)
  r <- .mh_crowding(pop)
  pop <- r$pop
  for (i in seq_along(pop)) pop[[i]]$crowding <- r$cd[i]
  pop
}

.mh_front_end <- function(pop, start) {
  r <- pop[[start]]$rank
  e <- start
  while (e < length(pop) && pop[[e + 1L]]$rank == r) e <- e + 1L
  e
}

# truncate_population!: the best N by rank, the last front cut by crowding
.mh_truncate_population <- function(pop, N) {
  pop <- .mh_fast_nds(pop)
  L <- length(pop)
  if (L == 0L) return(pop)
  ind <- 0L
  indnext <- .mh_front_end(pop, 1L)
  while (indnext <= N && indnext < L) {
    ind <- indnext
    indnext <- .mh_front_end(pop, ind + 1L)
  }
  if (indnext <= N) ind <- indnext
  if (ind + 1L <= indnext) {
    seg <- .mh_update_crowding(pop[(ind + 1L):indnext])
    cr <- vapply(seg, function(s) s$crowding, 0)
    pop[(ind + 1L):indnext] <- seg[order(-cr)]
  }
  pop[seq_len(min(N, length(pop)))]
}

.mh_get_non_dominated_perm <- function(pop) {
  n <- length(pop)
  if (n == 0L) return(integer(0))
  ids <- 1L
  for (i in seq_len(n)[-1L]) {
    j <- 1L
    dominated <- FALSE
    while (j <= length(ids)) {
      rel <- .mh_compare(pop[[i]], pop[[ids[j]]])
      if (rel == 2L) { dominated <- TRUE; break }
      if (rel == 1L) { ids <- ids[-j]; next }
      j <- j + 1L
    }
    if (!dominated) ids <- c(ids, i)
  }
  ids
}

# hypervolume contributions of the points (columns of P)
.mh_hv_contrib_2d <- function(P) {
  N <- ncol(P)
  dS <- rep(Inf, N)
  if (N < 3L) return(dS)
  p <- order(P[1L, ], P[2L, ])
  for (t in 2:(N - 1L)) {
    i <- p[t]
    dS[i] <- (P[1L, p[t + 1L]] - P[1L, i]) * (P[2L, p[t - 1L]] - P[2L, i])
  }
  dS
}

.mh_hv3d <- function(T3, ref) {
  n <- ncol(T3)
  if (n == 0L) return(0)
  r1 <- ref[1L]; r2 <- ref[2L]; r3 <- ref[3L]
  p <- order(T3[3L, ])
  xs <- numeric(n); ys <- numeric(n)
  k <- 0L; area <- 0; vol <- 0
  zprev <- T3[3L, p[1L]]
  for (t in seq_len(n)) {
    i <- p[t]
    x <- T3[1L, i]; y <- T3[2L, i]; z <- T3[3L, i]
    vol <- vol + area * (z - zprev)
    zprev <- z
    j <- if (k == 0L) 1L else sum(xs[seq_len(k)] < x) + 1L
    if (j > 1L && ys[j - 1L] <= y) next
    if (j <= k && xs[j] == x && ys[j] <= y) next
    e <- j
    while (e <= k && ys[e] >= y) e <- e + 1L
    xe <- if (e <= k) xs[e] else r1
    hgt <- if (j > 1L) ys[j - 1L] else r2
    area <- area + ((if (j <= k) xs[j] else r1) - x) * (hgt - y)
    if (j <= e - 1L) for (m in j:(e - 1L)) area <- area + ((if (m + 1L < e) xs[m + 1L] else xe) - xs[m]) * (ys[m] - y)
    nrem <- e - j
    if (nrem == 0L) {
      if (k >= j) for (m in k:j) { xs[m + 1L] <- xs[m]; ys[m + 1L] <- ys[m] }
    } else if (nrem > 1L) {
      if (e <= k) for (m in e:k) { xs[m - nrem + 1L] <- xs[m]; ys[m - nrem + 1L] <- ys[m] }
    }
    k <- k + 1L - nrem
    xs[j] <- x; ys[j] <- y
  }
  vol + area * (r3 - zprev)
}

.mh_hv_contrib_3d <- function(P, ref) {
  N <- ncol(P)
  dS <- numeric(N)
  for (i in seq_len(N)) {
    a <- P[1L, i]; b <- P[2L, i]; cc <- P[3L, i]
    q <- seq_len(N)[-i]
    Tm <- rbind(pmax(P[1L, q], a), pmax(P[2L, q], b), pmax(P[3L, q], cc))
    keep <- Tm[1L, ] < ref[1L] & Tm[2L, ] < ref[2L] & Tm[3L, ] < ref[3L]
    box <- (ref[1L] - a) * (ref[2L] - b) * (ref[3L] - cc)
    dS[i] <- max(0, box - .mh_hv3d(Tm[, keep, drop = FALSE], ref))
  }
  dS
}

.mh_hv_contrib_mc <- function(P, ref, n_samples) {
  M <- nrow(P); N <- ncol(P)
  lo <- apply(P, 1L, min)
  w <- ref - lo
  S <- matrix(stats::runif(n_samples * M), n_samples, M, byrow = TRUE)
  S <- sweep(sweep(S, 2L, w, "*"), 2L, lo, "+")
  dom <- matrix(TRUE, n_samples, N)
  for (m in seq_len(M)) dom <- dom & outer(S[, m], P[m, ], ">=")
  cnt <- rowSums(dom)
  owner <- max.col(dom, ties.method = "first")
  hits <- tabulate(owner[cnt == 1L], N)
  hits * (prod(w) / n_samples)
}

.mh_compute_contribution <- function(pop, n_samples = 10000L) {
  if (!length(pop)) return(numeric(0))
  P <- t(.mh_fvals(pop))
  M <- nrow(P); N <- ncol(P)
  if (M == 2L) return(.mh_hv_contrib_2d(P))
  if (N == 1L) return(Inf)
  ref <- vapply(seq_len(M), function(m) {
    lo <- min(P[m, ]); hi <- max(P[m, ])
    hi + if (hi > lo) 0.1 * (hi - lo) else 1
  }, 0)
  if (M == 3L) return(.mh_hv_contrib_3d(P, ref))
  .mh_hv_contrib_mc(P, ref, n_samples)
}

# gen_weights / gen_ref_dirs: the points of the simplex with H partitions
.mh_gen_weights <- function(nobj, H) {
  w <- list()
  a <- numeric(nobj)
  rec <- function(a, i, d) {
    for (k in 0:d) {
      if (i < nobj) {
        a[i] <- k
        rec(a, i + 1L, d - k)
      } else {
        a[i] <- d
        w[[length(w) + 1L]] <<- a / H
        break
      }
    }
  }
  rec(a, 1L, H)
  w
}

# Inverted generational distance of `front` to `reference`, as
# PerformanceIndicators.igd (rows are points)
.mh_generational_distance <- function(front, true_front, p = 1, inverted = FALSE, plus = FALSE) {
  if (inverted) { tmp <- true_front; true_front <- front; front <- tmp }
  s <- 0
  for (i in seq_len(nrow(front))) {
    x <- front[i, ]
    dif <- if (plus && inverted) pmax(sweep(true_front, 2L, x, "-"), 0) else if (plus) pmax(-sweep(true_front, 2L, x, "-"), 0) else
      sweep(true_front, 2L, x, "-")
    s <- s + min(sqrt(rowSums(dif * dif)))^p
  }
  s^(1 / p) / nrow(front)
}

# ---------------------------------------------------------------------------
# Options, information and termination

.mh_option_fields <- c("x_tol", "f_tol", "g_tol", "h_tol", "f_tol_rel", "f_calls_limit", "time_limit", "iterations",
                       "store_convergence", "debug", "seed", "parallel_evaluation", "verbose", "rng")

.mh_Options <- function(x_tol = 1e-8, f_tol = 1e-12, g_tol = 0, h_tol = 0, f_tol_rel = .Machine$double.eps, f_tol_abs = 0,
                        f_calls_limit = 0, g_calls_limit = 0, h_calls_limit = 0, iterations = 0L, time_limit = Inf,
                        store_convergence = FALSE, show_results = TRUE, debug = FALSE, search_type = "minimize",
                        parallel_evaluation = FALSE, verbose = FALSE, seed = NULL)
  structure(list(x_tol = as.numeric(x_tol), f_tol = as.numeric(f_tol), g_tol = as.numeric(g_tol),
                 h_tol = as.numeric(h_tol), f_tol_rel = as.numeric(f_tol_rel), f_calls_limit = as.numeric(f_calls_limit),
                 time_limit = as.numeric(time_limit), iterations = as.integer(iterations),
                 store_convergence = isTRUE(store_convergence), debug = isTRUE(debug), seed = seed,
                 parallel_evaluation = isTRUE(parallel_evaluation), verbose = isTRUE(verbose)),
            class = "MetaheuristicsOptions")

.mh_Information <- function(f_optimum = NaN, x_optimum = numeric(0))
  structure(list(f_optimum = as.numeric(f_optimum), x_optimum = as.numeric(x_optimum)), class = "MetaheuristicsInformation")

.mh_criterion <- function(cl, ...) structure(list(...), class = c(cl, "MetaheuristicsTermination"))

.mh_CheckConvergence <- function(f_tol_rel = .Machine$double.eps, f_tol_abs = 0, x_tol = 1e-8, criteria = NULL)
  .mh_criterion("CheckConvergence", criteria = criteria %mh_or% list(
    .mh_criterion("AbsoluteFunctionConvergence", ftol = f_tol_abs),
    .mh_criterion("RelativeFunctionConvergence", ftol = f_tol_rel),
    .mh_criterion("SmallStandardDeviation", ftol = 1e-6),
    .mh_criterion("RelativeParameterConvergence", xtol = x_tol)))

.mh_RobustConvergence <- function(ftol = 1e-3, period = 10L)
  .mh_criterion("RobustConvergence", ftol = ftol, period = as.integer(period))

.mh_Termination <- function(checkall = list(), checkany = list())
  structure(list(checkall = checkall, checkany = checkany), class = "MetaheuristicsTerminationSet")

# stop_check of a convergence criterion on the feasible solutions of a
# population
.mh_criterion_check <- function(pop, cr, mem) {
  switch(class(cr)[1L],
    CheckConvergence = {
      for (c1 in cr$criteria) if (!.mh_criterion_check(pop, c1, mem)) return(FALSE)
      TRUE
    },
    RelativeFunctionConvergence = {
      fx <- as.numeric(.mh_fvals(pop))
      m <- min(fx); M <- max(fx)
      isTRUE(abs(M - m) / max(abs(M), .Machine$double.eps) <= cr$ftol)
    },
    SmallStandardDeviation = isTRUE(stats::sd(as.numeric(.mh_fvals(pop))) <= cr$ftol),
    AbsoluteFunctionConvergence = {
      fx <- as.numeric(.mh_fvals(pop))
      isTRUE(abs(max(fx) - min(fx)) <= cr$ftol)
    },
    RelativeParameterConvergence = {
      fx <- as.numeric(.mh_fvals(pop))
      n <- length(pop)
      xhi <- pop[[(which.max(fx) - 1L) %% n + 1L]]$x
      xlo <- pop[[(which.min(fx) - 1L) %% n + 1L]]$x
      isTRUE(max(abs(xhi - xlo)) / max(max(pmax(abs(xlo), abs(xhi))), .Machine$double.eps) <= cr$xtol)
    },
    RobustConvergence = {
      if (is.null(mem$last)) mem$last <- pop
      mem$iter <- (mem$iter %mh_or% 0L) + 1L
      if (mem$iter %% cr$period != 0L) return(FALSE)
      stop <- .mh_robust_check(.mh_fvals(pop), .mh_fvals(mem$last), cr$ftol)
      mem$last <- pop
      stop
    },
    stop("unknown termination criterion ", class(cr)[1L], call. = FALSE))
}

.mh_robust_check <- function(fc, fl, tol) {
  fc <- as.matrix(fc); fl <- as.matrix(fl)
  fmin <- apply(fc, 2L, min); fmax <- apply(fc, 2L, max)
  dl <- fmax - fmin
  dl[dl == 0] <- 1
  if (!isTRUE(max(abs(apply(fl, 2L, min) / dl - fmin / dl)) <= tol)) return(FALSE)
  if (!isTRUE(max(abs(apply(fl, 2L, max) / dl - fmax / dl)) <= tol)) return(FALSE)
  fc <- sweep(sweep(fc, 2L, fmin, "-"), 2L, dl, "/")
  fl <- sweep(sweep(fl, 2L, fmin, "-"), 2L, dl, "/")
  isTRUE(.mh_generational_distance(fl, fc, inverted = TRUE) <= tol)
}

.mh_convergence_stop <- function(st, cr, mem, report = TRUE) {
  pop <- st$population
  if (!length(pop)) return(FALSE)
  feas <- vapply(pop, .mh_feasible, TRUE)
  if (sum(feas) / length(pop) < 0.3) return(FALSE)
  stop <- .mh_criterion_check(pop[feas], cr, mem)
  if (report && stop) st$termination_status_code <- class(cr)[1L]
  stop
}

.mh_term_stop <- function(st, term) {
  ca <- term$checkall
  if (length(ca)) {
    for (k in seq_along(ca)) if (!.mh_convergence_stop(st, ca[[k]], term$mem_all[[k]])) return(FALSE)
    return(TRUE)
  }
  for (k in seq_along(term$checkany)) if (.mh_convergence_stop(st, term$checkany[[k]], term$mem_any[[k]])) return(TRUE)
  FALSE
}

.mh_call_limit_stop <- function(st, info, op) {
  cond <- st$f_calls >= op$f_calls_limit
  if (cond) st$termination_status_code <- "EVALUATIONS_LIMIT"
  cond
}

.mh_iteration_stop <- function(st, info, op) {
  cond <- st$iteration >= op$iterations
  if (cond) st$termination_status_code <- "ITERATION_LIMIT"
  cond
}

.mh_time_stop <- function(st, info, op) {
  cond <- st$overall_time >= op$time_limit
  if (cond) st$termination_status_code <- "TIME_LIMIT"
  cond
}

.mh_accuracy_stop <- function(st, info, op) {
  b <- st$best_sol
  if (is.nan(info$f_optimum) || .mh_is_mo(b)) return(FALSE)
  if (!is.null(b$sv) && .mh_violations_sum(b$g, b$h, op$h_tol) != 0) return(FALSE)
  cond <- isTRUE(abs(b$f - info$f_optimum) < op$f_tol)
  if (cond) st$termination_status_code <- "ACCURACY_LIMIT"
  cond
}

.mh_default_stop_check <- function(st, info, op)
  .mh_call_limit_stop(st, info, op) || .mh_iteration_stop(st, info, op) || .mh_time_stop(st, info, op) ||
    .mh_accuracy_stop(st, info, op) ||
    .mh_convergence_stop(st, .mh_CheckConvergence(f_tol_abs = op$f_tol, f_tol_rel = op$f_tol_rel, x_tol = op$x_tol),
                         new.env())

# stop_criteria! of the template: the budget, iteration, time and accuracy
.mh_template_stop <- function(st, par, pr, info, op) {
  st$stop <- st$stop || .mh_call_limit_stop(st, info, op) || .mh_iteration_stop(st, info, op) ||
    .mh_time_stop(st, info, op) || .mh_accuracy_stop(st, info, op)
  invisible(NULL)
}

.mh_status_messages <- c(
  ITERATION_LIMIT = "Maximum number of iterations exceeded.",
  TIME_LIMIT = "Maximum time exceeded.",
  EVALUATIONS_LIMIT = "Maximum objective function calls exceeded.",
  ACCURACY_LIMIT = "The desired accuracy was obtained.",
  OBJECTIVE_VARIANCE_LIMIT = "Small variance of the objective function.",
  OBJECTIVE_DIFFERENCE_LIMIT = "Small difference of objective function values.",
  OTHER_LIMIT = "Other stopping criteria.",
  UNKNOWN_STOP_REASON = "Unknown stop reason.")

.mh_termination_message <- function(code)
  if (code %in% names(.mh_status_messages)) .mh_status_messages[[code]] else "Due to Convergence Termination criterion."

# ---------------------------------------------------------------------------
# State and the optimize loop

.mh_state <- function(best_sol, population, f_calls = 0L, iteration = 1L) {
  st <- new.env(parent = emptyenv())
  st$best_sol <- best_sol
  st$population <- population
  st$f_calls <- as.integer(f_calls)
  st$iteration <- as.integer(iteration)
  st$success_rate <- 0
  st$convergence <- list()
  st$start_time <- 0
  st$final_time <- 0
  st$overall_time <- 0
  st$stop <- FALSE
  st$termination_status_code <- "UNKNOWN_STOP_REASON"
  st
}

# gen_initial_state: a random population, or the population of the state
# completed with random solutions
.mh_gen_initial_state <- function(pr, par, info, op, st) {
  if (is.null(st$best_sol)) {
    pop <- .mh_generate_population(par$N, pr, op$h_tol)
    return(.mh_state(.mh_get_best(pop), pop, f_calls = length(pop), iteration = 1L))
  }
  if (pr$D != length(st$best_sol$x)) stop("Invalid population (dimension does not match with bounds)", call. = FALSE)
  pop <- st$population
  if (par$N < length(pop)) {
    par$N <- length(pop)
  } else if (par$N > length(pop)) {
    pop <- c(pop, .mh_generate_population(par$N - length(pop), pr, op$h_tol))
  }
  best <- .mh_get_best(pop)
  if (.mh_is_better(st$best_sol, best)) best <- st$best_sol
  .mh_state(best, pop, f_calls = length(pop))
}

.mh_snapshot <- function(st) {
  out <- mget(ls(st, all.names = TRUE), envir = st)
  out$convergence <- NULL
  structure(out, class = "MetaheuristicsState")
}

# optimize(f, bounds, method) of Metaheuristics.jl
.mh_optimize <- function(f, lb, ub, alg, mo = FALSE) {
  impl <- .mh_impls[[alg$impl %mh_or% alg$name]]
  op <- list2env(unclass(alg$options), parent = emptyenv())
  par <- list2env(alg$parameters, parent = emptyenv())
  info <- alg$information
  pr <- .mh_problem(f, lb, ub, mo)
  run <- function() {
    st <- if (is.null(alg$initial_state)) .mh_state(NULL, list()) else
      .mh_state(alg$initial_state$best_sol, alg$initial_state$population, iteration = 0L)
    t0 <- .mh_now()
    st <- impl$initialize(st, par, pr, info, op)
    st$start_time <- t0
    st$overall_time <- st$overall_time + .mh_now() - t0
    term <- alg$termination %mh_or% .mh_Termination()
    if (!length(term$checkall) && !length(term$checkany)) {
      term$checkany <- list(if (length(st$population) && .mh_is_mo(st$population[[1L]])) .mh_RobustConvergence(ftol = op$f_tol) else
        .mh_CheckConvergence(f_tol_abs = op$f_tol, f_tol_rel = op$f_tol_rel, x_tol = op$x_tol))
    }
    term$mem_all <- lapply(term$checkall, function(i) new.env())
    term$mem_any <- lapply(term$checkany, function(i) new.env())
    st$iteration <- 1L
    conv <- list()
    if (op$store_convergence) conv <- list(.mh_snapshot(st))
    st$stop <- st$stop || .mh_default_stop_check(st, info, op) || .mh_term_stop(st, term)
    stop_fn <- impl$stop_criteria %mh_or% .mh_template_stop
    while (!st$stop) {
      st$iteration <- st$iteration + 1L
      t <- .mh_now()
      impl$update(st, par, pr, info, op)
      st$overall_time <- st$overall_time + .mh_now() - t
      st$f_calls <- pr$f_calls
      if (op$store_convergence) conv[[length(conv) + 1L]] <- .mh_snapshot(st)
      st$stop <- st$stop || .mh_default_stop_check(st, info, op) || .mh_term_stop(st, term)
      if (!st$stop) stop_fn(st, par, pr, info, op)
    }
    t <- .mh_now()
    if (!is.null(impl$final)) impl$final(st, par, pr, info, op)
    st$overall_time <- st$overall_time + .mh_now() - t
    st$final_time <- .mh_now()
    out <- .mh_snapshot(st)
    out$convergence <- conv
    out$parameters <- as.list(par)
    out$options <- as.list(op)
    out
  }
  .with_seed(op$seed, run())
}

# ---------------------------------------------------------------------------
# Algorithm objects

.mh_algorithm <- function(name, parameters, information = NULL, options = NULL, termination = NULL,
                          initial_state = NULL, cons = FALSE, impl = NULL) {
  if (!is.null(termination) && !inherits(termination, "MetaheuristicsTerminationSet"))
    termination <- .mh_Termination(checkany = list(termination))
  structure(list(name = name, impl = impl, parameters = parameters, information = information %mh_or% .mh_Information(),
                 options = options %mh_or% .mh_Options(), termination = termination, initial_state = initial_state,
                 traits = list(bounds = TRUE, rbounds = TRUE, cons = cons, callback = FALSE)),
            class = c(name, "MetaheuristicsAlgorithm", "OptimizationAlgorithm"))
}

#' @export
print.MetaheuristicsAlgorithm <- function(x, ...) {
  pr <- x$parameters[!vapply(x$parameters, function(v) is.list(v) || is.function(v) || length(v) != 1L, TRUE)]
  cat(x$name, "(", paste0(names(pr), "=", vapply(pr, format, ""), collapse = ", "), ")\n", sep = "")
  invisible(x)
}

#' @export
print.MetaheuristicsState <- function(x, ...) {
  cat("Metaheuristics state after ", x$iteration, " iterations and ", x$f_calls, " evaluations\n", sep = "")
  if (.mh_is_mo(x$best_sol)) {
    nd <- .mh_get_non_dominated_perm(x$population)
    cat("Nondominated solutions: ", length(nd), " of ", length(x$population), "\n", sep = "")
  } else {
    cat("Minimum: ", format(x$best_sol$f, digits = 10), "\nMinimizer: ", paste(format(x$best_sol$x, digits = 10), collapse = " "),
        "\n", sep = "")
  }
  cat("Stop reason: ", .mh_termination_message(x$termination_status_code), "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# The OptimizationMetaheuristics wrapper

.mh_info_fields <- c("f_optimum", "x_optimum")

# Map the solve arguments: Options fields are set, Information fields and
# unknown keywords are errors, and maxiters, maxtime and abstol become
# iterations, time_limit and f_tol.
.mh_map_args <- function(alg, opts) {
  o <- alg$options
  jl <- paste0("Metaheuristics.Algorithm{", alg$jl_name %mh_or% alg$name, "}")
  for (nm in setdiff(names(opts), c("maxiters", "maxtime", "abstol", "reltol", "use_initial", "progress"))) {
    if (nm %in% .mh_info_fields)
      stop("Set ", nm, " by directly passing it to Information Structure which is passed to ", jl,
           " algorithms when calling solve().", call. = FALSE)
    if (nm %in% .mh_option_fields) {
      o[nm] <- list(opts[[nm]])
    } else {
      stop(nm, " keyword is not a valid option for Metaheuristics.AbstractAlgorithm algorithm.", call. = FALSE)
    }
  }
  if (!is.null(opts$maxiters)) o$iterations <- as.integer(opts$maxiters)
  if (!is.null(opts$maxtime)) o$time_limit <- as.numeric(opts$maxtime)
  if (!is.null(opts$abstol)) o$f_tol <- as.numeric(opts$abstol)
  alg$options <- o
  alg
}

# constraint_violations: g from the finite one sided bounds, h from the
# equalities
.mh_constraint_violations <- function(a) {
  fn <- a$f
  p <- a$p
  if (!is.function(fn$cons)) return(NULL)
  lcons <- a$lcons; ucons <- a$ucons
  m <- length(lcons %mh_or% ucons %mh_or% fn$cons(as.numeric(a$u0), p))
  lcons <- lcons %mh_or% rep(0, m); ucons <- ucons %mh_or% rep(0, m)
  eq <- which(lcons == ucons)
  up <- which(lcons != ucons & is.finite(ucons))
  lo <- which(lcons != ucons & is.finite(lcons))
  function(x) {
    res <- as.numeric(fn$cons(x, p))
    list(g = c(res[up] - ucons[up], lcons[lo] - res[lo]), h = res[eq] - ucons[eq])
  }
}

# initial_population!: run a copy of the algorithm for two iterations to
# learn the population size, then start from u0 and random points
.mh_initial_population <- function(alg, a, loss, mo) {
  ai <- alg
  ai$options$iterations <- 2L
  r <- .mh_optimize(loss, a$lb, a$ub, ai, mo)
  N <- r$parameters$N
  u0 <- as.numeric(a$u0)
  X <- lapply(seq_len(N - 1L), function(i) a$lb + stats::runif(length(u0)) * (a$ub - a$lb))
  X[[N]] <- u0
  pop <- lapply(X, function(x) .mh_child(x, loss(x), 0, mo))
  alg$parameters$N <- N
  alg$initial_state <- list(best_sol = .mh_get_best(pop), population = pop)
  alg
}

.mtk$hooks$opt_solve_MetaheuristicsAlgorithm <- function(a, alg, opts) {
  .opt_check_traits(a, alg)
  mo <- inherits(a$f, "MultiObjectiveOptimizationFunction")
  p <- a$p
  viol <- .mh_constraint_violations(a)
  loss <- if (mo) {
    function(x) {
      r <- a$f$f(x, p)
      if (is.list(r)) { F <- r[[1L]]; g <- r[[2L]]; h <- r[[3L]] } else { F <- r; g <- numeric(0); h <- numeric(0) }
      if (is.null(viol)) return(list(as.numeric(F), as.numeric(g), as.numeric(h)))
      v <- viol(x)
      list(as.numeric(F), c(as.numeric(g), v$g), c(as.numeric(h), v$h))
    }
  } else {
    function(x) {
      fx <- as.numeric(a$f$f(x, p))[1L]
      if (is.null(viol)) return(fx)
      v <- viol(x)
      list(fx, v$g, v$h)
    }
  }
  alg <- .mh_map_args(alg, opts)
  t0 <- .mh_now()
  if (isTRUE(opts$use_initial)) alg <- .mh_initial_population(alg, a, loss, mo)
  res <- .mh_optimize(loss, a$lb, a$ub, alg, mo)
  stats <- OptimizationStats(time = .mh_now() - t0)
  if (mo) {
    nd <- res$population[.mh_get_non_dominated_perm(res$population)]
    res$X <- .mh_positions(nd)
    res$F <- .mh_fvals(nd)
    res$retcode <- "Default"
    res$stats <- stats
    return(res)
  }
  sol <- .opt_solution(a, alg, res$best_sol$x, res$best_sol$f, "Default", res)
  .opt_set_fields(sol, stats = stats)
}
