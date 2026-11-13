# Jump aggregators, after JumpProcesses.jl: aggregators/aggregators.jl
# (the algorithm types, needs_depgraph, needs_vartojumps_map,
# supports_variablerates, select_aggregator), aggregators/ssajump.jl (the
# shared routines: fill_rates_and_sum!, calculate_jump_rate,
# update_dependent_rates!, update_state!, nomorejumps!, the rounding error
# bounds of the running rate sums, linear_search, rejectrx,
# no_jump_before_end!, execute_tied_jumps!), aggregators/direct.jl,
# sortingdirect.jl, nrm.jl, directcr.jl, rssa.jl, rssacr.jl, coevolve.jl,
# bracketing.jl (BracketData), prioritytable.jl (PriorityTable),
# massaction_rates.jl (evalrxrate, executerx!, var_to_jumps_map,
# jump_to_vars_map, make_dependency_graph) and jump_dependencies.jl
# (normalize_map).
#
# Jumps are numbered as in JumpProcesses: mass action jumps first, then
# constant rate jumps, then (for Coevolve) bounded variable rate jumps. An
# aggregator is an environment holding its state; .agg_initialize,
# .agg_execute and .agg_generate are initialize!, execute_jumps! and
# generate_jumps!. The integrator `integ` is an environment with u, p and t;
# with jumps coupled to an ODE or SDE the state also holds the integrated
# intensities of the variable rate jumps after the n unknowns.

#' Jump aggregators and the SSA stepper
#'
#' Aggregators choose the next jump of the constant rate and mass action jumps
#' of a [JumpProblem()], as in JumpProcesses.jl; all simulate the exact law
#' of the jump process. `Direct()` is Gillespie's direct method and
#' `SortingDirect()` the direct method with the jumps searched in order of
#' how often they fire. `NRM()` is Gibson and Bruck's next reaction method
#' and `DirectCR()` the composition rejection direct method. `RSSA()` and
#' `RSSACR()` are the rejection SSA methods, which bracket each population in
#' an interval (see [BracketData()]) and assume that the rates are monotone in
#' the populations. `Coevolve()` is a next reaction method that also handles
#' variable rate jumps with rate bounds (`urate`, `rateinterval` and
#' optionally `lrate` of [VariableRateJump()]) by thinning.
#' `FRM()` is Gillespie's first reaction method, `RDirect()` the direct method
#' with rejection sampling of the jump, and `CCNRM()` the next reaction method
#' with a constant complexity priority queue. `DirectFW()` and `FRMFW()` are
#' Julia's function wrapper variants of `Direct()` and `FRM()`; in R they are
#' the same algorithms.
#' `NullAggregator()`, the default of `JumpProblem()`, selects one of them
#' from the number and kind of jumps as JumpProcesses' `select_aggregator`.
#'
#' `SortingDirect`, `NRM`, `CCNRM`, `RDirect`, `DirectCR` and `Coevolve` need a
#' dependency graph
#' (for each jump, the jumps whose rates change when it fires) and `RSSA` and
#' `RSSACR` the maps from unknowns to the jumps whose rates depend on them
#' and from jumps to the unknowns they change. For systems they are built as
#' ModelingToolkit builds them ([asgraph()], [variable_dependencies()],
#' [eqeq_dependencies()]); problems without a system pass them as
#' `dependencies = list(dep_graph = , vartojumps_map = , jumptovars_map = )`.
#'
#' `SSAStepper()` steps from jump to jump in pure jump problems and is their
#' default solver. `reset_aggregated_jumps(integ)` recomputes the rates and
#' next jump times inside a callback that changed the state or parameters
#' (affects of the system's own events call it). `add_tstop(integ, t)` adds
#' a time the solver must stop at; a time behind the current time is an
#' error, as in OrdinaryDiffEq.
#'
#' @param integ The integrator passed to a callback.
#' @param t A time.
#' @return An aggregator or solver object; `reset_aggregated_jumps()` and
#'   `add_tstop()` return `NULL` invisibly.
#' @examples
#' t <- t_nounits
#' variables("X(t)"); parameters("lam = 10", "mu = 1")
#' birth <- SymbolicMassActionJump(lam, list(), list(X ~ 1))
#' death <- SymbolicMassActionJump(mu, list(X ~ 1), list(X ~ -1))
#' js <- complete(JumpSystem(list(birth, death), t, list(X), list(lam, mu), name = "bd"))
#' jprob <- JumpProblem(js, list(X ~ 0), c(0, 5), aggregator = NRM())
#' sol <- solve(jprob, SSAStepper(), seed = 1)
#' tail(sol[X], 1)
#' @name jump_aggregators
NULL

.new_aggregator <- function(name) structure(list(name = name), class = c(name, "JumpAggregator"))

#' @rdname jump_aggregators
#' @export
Direct <- function() .new_aggregator("Direct")
#' @rdname jump_aggregators
#' @export
SortingDirect <- function() .new_aggregator("SortingDirect")
#' @rdname jump_aggregators
#' @export
NRM <- function() .new_aggregator("NRM")
#' @rdname jump_aggregators
#' @export
DirectCR <- function() .new_aggregator("DirectCR")
#' @rdname jump_aggregators
#' @export
RSSA <- function() .new_aggregator("RSSA")
#' @rdname jump_aggregators
#' @export
RSSACR <- function() .new_aggregator("RSSACR")
#' @rdname jump_aggregators
#' @export
Coevolve <- function() .new_aggregator("Coevolve")
#' @rdname jump_aggregators
#' @export
FRM <- function() .new_aggregator("FRM")
#' @rdname jump_aggregators
#' @export
FRMFW <- function() .new_aggregator("FRMFW")
#' @rdname jump_aggregators
#' @export
DirectFW <- function() .new_aggregator("DirectFW")
#' @rdname jump_aggregators
#' @export
RDirect <- function() .new_aggregator("RDirect")
#' @rdname jump_aggregators
#' @export
CCNRM <- function() .new_aggregator("CCNRM")
#' @rdname jump_aggregators
#' @export
NullAggregator <- function() .new_aggregator("NullAggregator")
#' @rdname jump_aggregators
#' @export
SSAStepper <- function() structure(list(name = "SSAStepper"), class = "SSAStepper")

#' @export
print.JumpAggregator <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

.needs_depgraph <- function(agg)
  agg$name %in% c("DirectCR", "SortingDirect", "NRM", "CCNRM", "RDirect", "Coevolve")
.needs_vartojumps_map <- function(agg) agg$name %in% c("RSSA", "RSSACR")
.supports_variablerates <- function(agg) agg$name == "Coevolve"

# select_aggregator with USE_DIRECT_THRESHOLD = 20, USE_RSSA_THRESHOLD = 100
# and USE_SORTINGDIRECT_THRESHOLD = 200 (JumpProcesses.jl).
.select_aggregator <- function(nmaj, ncrj, nvrj, nbvrj, deps) {
  if (nvrj + nbvrj > 0) return(if (nbvrj > 0) Coevolve() else Direct())
  nj <- nmaj + ncrj
  if (nj < 20) return(Direct())
  have_maps <- !is.null(deps$vartojumps_map) && !is.null(deps$jumptovars_map)
  if (ncrj == 0 || have_maps) return(if (nj < 100) RSSA() else RSSACR())
  if (!is.null(deps$dep_graph)) return(if (nj < 200) SortingDirect() else DirectCR())
  Direct()
}

#' Population brackets for the rejection SSA methods
#'
#' `BracketData(fluctrate, threshold, delta_u)` sets the intervals in which
#' [RSSA()] and [RSSACR()] bracket each population `u`: `[(1 - fluctrate) u,
#' (1 + fluctrate) u]` when `u >= threshold`, and `[max(u - delta_u, 0), u +
#' delta_u]` below it. Each argument is a number or one value per unknown.
#' Pass it to [JumpProblem()] as `bracket_data`.
#'
#' @param fluctrate Relative half width for large populations.
#' @param threshold Population below which the absolute width is used.
#' @param delta_u Absolute half width for small populations.
#' @return A `BracketData` object.
#' @examples
#' BracketData(0.1, 25, 4)
#' @export
BracketData <- function(fluctrate = 0.1, threshold = 25, delta_u = 4)
  structure(list(fluctrate = fluctrate, threshold = threshold, delta_u = delta_u), class = "BracketData")

.bd_get <- function(x, i) if (length(x) > 1L) x[[i]] else x

# get_spec_brackets for real valued populations
.spec_brackets <- function(bd, i, u) {
  if (u == 0) return(c(0, 0))
  if (u < .bd_get(bd$threshold, i)) {
    du <- .bd_get(bd$delta_u, i)
    return(c(max(0, u - du), u + du))
  }
  d <- .bd_get(bd$fluctrate, i)
  c((1 - d) * u, (1 + d) * u)
}

# ---------------------------------------------------------------------------
# Dependency maps (jump_dependencies.jl normalize_map, massaction_rates.jl)

.normalize_map <- function(map, addself) {
  lapply(seq_along(map), function(i) {
    e <- as.integer(unlist(map[[i]]))
    if (addself) e <- c(e, i)
    sort(unique(e))
  })
}

.check_map <- function(map, name, nentries, nindices) {
  if (is.null(map)) return(invisible())
  if (length(map) != nentries)
    .stopf("`%s` has %d entries, but needs %d", name, length(map), nentries)
  for (i in seq_along(map)) {
    bad <- setdiff(unlist(map[[i]]), seq_len(nindices))
    if (length(bad)) .stopf("Entry %d of `%s` lists %s, which is not a valid index", i, name, format(bad[1L]))
  }
  invisible()
}

# The maps an aggregator uses, in normal form; generated from mass action
# jumps when every jump is one (normalize_dependencies).
.agg_dependencies <- function(agg, deps, js) {
  nj <- js$nmaj + js$ncrj + js$nbvrj
  if (nj == 0L) return(list())
  can_generate <- js$nmaj > 0L && nj == js$nmaj
  need <- function(name) {
    m <- deps[[name]]
    if (!is.null(m)) return(m)
    if (!can_generate) {
      if (name == "dep_graph")
        .stopf("To use ConstantRateJumps with the %s algorithm a dependency graph must be supplied.", agg$name)
      .stopf("To use the %s algorithm a map from %s to dependent %s must be supplied.", agg$name,
             if (name == "vartojumps_map") "variables" else "jumps",
             if (name == "vartojumps_map") "jumps" else "variables")
    }
    .maj_generated_map(name, js)
  }
  out <- list()
  if (.needs_depgraph(agg)) {
    out$dep_graph <- .normalize_map(need("dep_graph"), TRUE)
    .check_map(out$dep_graph, "dep_graph", nj, nj)
  }
  if (.needs_vartojumps_map(agg)) {
    out$vartojumps_map <- .normalize_map(need("vartojumps_map"), FALSE)
    out$jumptovars_map <- .normalize_map(need("jumptovars_map"), FALSE)
    .check_map(out$vartojumps_map, "vartojumps_map", js$n, nj)
    .check_map(out$jumptovars_map, "jumptovars_map", nj, js$n)
  }
  out
}

# var_to_jumps_map, jump_to_vars_map and make_dependency_graph
.maj_generated_map <- function(name, js) {
  v2j <- rep(list(integer()), js$n)
  for (k in seq_len(js$nmaj)) for (s in js$maj_rs[[k]][, 1L]) v2j[[s]] <- c(v2j[[s]], k)
  if (name == "vartojumps_map") return(v2j)
  j2v <- lapply(seq_len(js$nmaj), function(k) js$maj_ns[[k]][, 1L])
  if (name == "jumptovars_map") return(j2v)
  lapply(seq_len(js$nmaj), function(k) unlist(lapply(j2v[[k]], function(s) v2j[[s]])))
}

# ---------------------------------------------------------------------------
# Mass action rates (massaction_rates.jl), the real valued evalrxrate: the
# rate is zero when the smallest factor of a reactant is not positive.

.evalrxrate <- function(u, rs, rate) {
  val <- 1
  for (r in seq_len(nrow(rs))) {
    pop <- u[rs[r, 1L]]
    val <- val * pop
    s <- rs[r, 2L]
    while (s > 1) {
      pop <- pop - 1
      val <- val * pop
      s <- s - 1
    }
    if (pop <= 0) return(0)
  }
  val * rate
}

# ---------------------------------------------------------------------------
# Rounding error bounds of running rate sums (ssajump.jl)

.update_error_bound <- function(before, after, nops) nops * .Machine$double.eps * (abs(before) + abs(after))
.sum_error_bound <- function(s, n) n * .Machine$double.eps * abs(s)
.sum_rate_inaccurate <- function(s, err, n)
  err > max(sqrt(.Machine$double.eps), 4 * n * .Machine$double.eps) * s

# index j with sum(v[1:(j-1)]) <= r < sum(v[1:j]); 0 when r is not below the sum
.linear_search <- function(v, r) {
  j <- findInterval(r, cumsum(v)) + 1L
  if (j > length(v)) 0L else j
}

# ---------------------------------------------------------------------------
# Aggregator state

# A new aggregator for the numeric jump set `js` (see .jmp_numeric).
.agg_new <- function(agg, js, deps, save_positions, bracket_data = NULL) {
  a <- new.env(parent = emptyenv())
  a$type <- switch(agg$name, DirectFW = "Direct", FRMFW = "FRM", agg$name)
  a$js <- js
  a$nmaj <- js$nmaj
  a$ncrj <- js$ncrj
  a$nbvrj <- if (agg$name == "Coevolve") js$nbvrj else 0L
  a$J <- a$nmaj + a$ncrj + a$nbvrj
  a$deps <- deps
  a$save_positions <- save_positions
  a$next_jump <- 0L
  a$prev_jump <- 0L
  a$next_jump_time <- Inf
  a$end_time <- Inf
  a$sum_rate <- 0
  a$sum_rate_err <- 0
  a$bd <- bracket_data %||% BracketData()
  a$jump_search_order <- seq_len(a$J)
  a$jump_search_idx <- 1L
  a$minexponent <- -40L  # exponent(MINJUMPRATE = 2^exponent(1e-12))
  a
}

.agg_u <- function(a, integ) if (length(integ$u) == a$js$n) integ$u else integ$u[seq_len(a$js$n)]

# calculate_jump_rate
.agg_rate <- function(a, k, u, p, t) {
  if (k <= a$nmaj) return(.evalrxrate(u, a$js$maj_rs[[k]], a$maj_rates[k]))
  a$js$crj_rate[[k - a$nmaj]](u, p, t)
}

.agg_all_rates <- function(a, u, p, t) vapply(seq_len(a$J), function(k) .agg_rate(a, k, u, p, t), 0)

# fill_scaled_rates!
.agg_fill_scaled <- function(a, u, p, t) {
  a$maj_rates <- if (a$nmaj) as.numeric(a$js$maj_ratefun(u, p, t)) else numeric()
}

# update_state!: execute jump next_jump
.agg_update_state <- function(a, integ) {
  k <- a$next_jump
  if (k <= a$nmaj) {
    ns <- a$js$maj_ns[[k]]
    integ$u[ns[, 1L]] <- integ$u[ns[, 1L]] + ns[, 2L]
  } else if (k <= a$nmaj + a$ncrj) {
    a$js$crj_affect[[k - a$nmaj]](integ)
  } else {
    a$js$bvrj[[k - a$nmaj - a$ncrj]]$affect(integ)
  }
  a$prev_jump <- k
  invisible()
}

.agg_nomorejumps <- function(a, s) {
  if (s <= 0) {
    a$next_jump <- 0L
    a$next_jump_time <- Inf
    return(TRUE)
  }
  FALSE
}

# update_dependent_rates! of ssajump.jl (SortingDirect)
.agg_update_dependent <- function(a, deps, u, p, t) {
  s <- a$sum_rate
  for (rx in deps) {
    s <- s - a$cur_rates[rx]
    a$cur_rates[rx] <- .agg_rate(a, rx, u, p, t)
    s <- s + a$cur_rates[rx]
  }
  a$sum_rate_err <- a$sum_rate_err + .update_error_bound(a$sum_rate, s, 2 * length(deps))
  a$sum_rate <- s
}

.agg_resum <- function(a, rates) {
  a$sum_rate <- sum(rates)
  a$sum_rate_err <- .sum_error_bound(a$sum_rate, length(rates))
}

# initialize!
.agg_initialize <- function(a, integ) {
  if (!is.null(a$spatial)) return(.sp_initialize(a, integ))
  u <- .agg_u(a, integ); p <- integ$p; t <- integ$t
  a$end_time <- integ$prob$tspan[2L]
  .agg_fill_scaled(a, u, p, t)
  switch(a$type,
    Direct = , FRM = NULL,
    SortingDirect = {
      a$cur_rates <- .agg_all_rates(a, u, p, t)
      a$sum_rate <- sum(a$cur_rates)
      a$sum_rate_err <- .sum_error_bound(a$sum_rate, a$J)
    },
    RDirect = {
      a$cur_rates <- .agg_all_rates(a, u, p, t)
      a$sum_rate <- sum(a$cur_rates)
      a$sum_rate_err <- .sum_error_bound(a$sum_rate, a$J)
      a$max_rate <- if (a$J) max(a$cur_rates) else 0
      a$counter <- 0L
    },
    NRM = , CCNRM = {
      a$cur_rates <- numeric(a$J)
      a$pq <- numeric(a$J)
      for (i in seq_len(a$J)) {
        a$cur_rates[i] <- .agg_rate(a, i, u, p, t)
        a$pq[i] <- t + stats::rexp(1L) / a$cur_rates[i]
      }
    },
    DirectCR = {
      a$cur_rates <- .agg_all_rates(a, u, p, t)
      a$sum_rate <- sum(a$cur_rates)
      a$rt <- .pt_new(a$minexponent, a$cur_rates)
      a$sum_rate_err <- .sum_error_bound(a$sum_rate, .pt_numterms(a$rt))
    },
    RSSA = {
      .agg_set_bracketing(a, u, p, t)
      a$sum_rate_err <- .sum_error_bound(a$sum_rate, a$J)
    },
    RSSACR = {
      .agg_set_bracketing(a, u, p, t)
      a$rt <- .pt_new(a$minexponent, a$cur_rate_high)
      a$sum_rate_err <- .sum_error_bound(a$sum_rate, .pt_numterms(a$rt))
    },
    Coevolve = .coev_fill_rates_and_get_times(a, u, p, t),
    .stopf("unknown aggregator %s", a$type))
  .agg_generate(a, integ)
}

# execute_jumps!
.agg_execute <- function(a, integ) {
  if (!is.null(a$spatial)) return(.sp_execute(a, integ))
  .agg_update_state(a, integ)
  u <- .agg_u(a, integ); p <- integ$p; t <- integ$t
  k <- a$next_jump
  switch(a$type,
    Direct = , FRM = NULL,
    RDirect = {
      deps <- a$deps$dep_graph[[k]]
      s <- a$sum_rate
      increased <- FALSE
      for (rx in deps) {
        new <- .agg_rate(a, rx, u, p, t)
        s <- s + new - a$cur_rates[rx]
        if (new > a$max_rate) {
          a$max_rate <- new
          increased <- TRUE
        }
        a$cur_rates[rx] <- new
      }
      # counter_threshold = number of jumps (rdirect.jl)
      if (!increased && a$counter > a$J) a$max_rate <- max(a$cur_rates)
      a$sum_rate_err <- a$sum_rate_err + .update_error_bound(a$sum_rate, s, 2 * length(deps))
      a$sum_rate <- s
    },
    SortingDirect = {
      jsidx <- a$jump_search_idx
      if (jsidx != 1L) a$jump_search_order[c(jsidx - 1L, jsidx)] <- a$jump_search_order[c(jsidx, jsidx - 1L)]
      .agg_update_dependent(a, a$deps$dep_graph[[k]], u, p, t)
    },
    NRM = , CCNRM = {
      for (rx in a$deps$dep_graph[[k]]) {
        old <- a$cur_rates[rx]
        new <- a$cur_rates[rx] <- .agg_rate(a, rx, u, p, t)
        a$pq[rx] <- if (rx != k && old > 0) {
          if (new > 0) t + old / new * (a$pq[rx] - t) else Inf
        } else {
          if (new > 0) t + stats::rexp(1L) / new else Inf
        }
      }
    },
    DirectCR = {
      deps <- a$deps$dep_graph[[k]]
      for (rx in deps) {
        old <- a$cur_rates[rx]
        a$cur_rates[rx] <- .agg_rate(a, rx, u, p, t)
        .pt_update(a$rt, rx, old, a$cur_rates[rx])
      }
      s <- a$rt$gsum
      a$sum_rate_err <- a$sum_rate_err + .update_error_bound(a$sum_rate, s, 4 * length(deps))
      a$sum_rate <- s
    },
    RSSA = , RSSACR = .rssa_update_rates(a, a$deps$jumptovars_map[[k]], u, p, t),
    Coevolve = {
      for (i in a$deps$dep_graph[[k]]) {
        nt <- .coev_next_time(a, u, p, t, i)
        a$pq[i] <- nt[1L]
        a$cur_rates[i] <- nt[2L]
      }
    })
  invisible()
}

# generate_jumps!
.agg_generate <- function(a, integ) {
  if (!is.null(a$spatial)) return(.sp_generate(a, integ))
  u <- .agg_u(a, integ); p <- integ$p; t <- integ$t
  switch(a$type,
    Direct = {
      cum <- cumsum(.agg_all_rates(a, u, p, t))
      a$sum_rate <- if (a$J) cum[a$J] else 0
      a$next_jump_time <- t + stats::rexp(1L) / a$sum_rate
      a$next_jump <- min(findInterval(stats::runif(1L) * a$sum_rate, cum, left.open = TRUE) + 1L, a$J)
    },
    FRM = {
      # first reaction: a putative time for every jump, the earliest fires
      if (a$J) {
        dt <- stats::rexp(a$J) / .agg_all_rates(a, u, p, t)
        k <- which.min(dt)
        a$next_jump <- k
        a$next_jump_time <- t + dt[k]
      }
    },
    RDirect = {
      if (.sum_rate_inaccurate(a$sum_rate, a$sum_rate_err, a$J)) {
        .agg_resum(a, a$cur_rates)
        a$max_rate <- max(a$cur_rates)
      }
      if (.agg_nomorejumps(a, a$sum_rate)) return(invisible())
      counter <- 0L
      rx <- sample.int(a$J, 1L)
      while (a$cur_rates[rx] < stats::runif(1L) * a$max_rate) {
        rx <- sample.int(a$J, 1L)
        counter <- counter + 1L
      }
      a$counter <- counter
      a$next_jump <- rx
      a$next_jump_time <- t + stats::rexp(1L) / a$sum_rate
    },
    SortingDirect = {
      if (.sum_rate_inaccurate(a$sum_rate, a$sum_rate_err, a$J)) .agg_resum(a, a$cur_rates)
      if (.agg_nomorejumps(a, a$sum_rate)) return(invisible())
      a$next_jump_time <- t + stats::rexp(1L) / a$sum_rate
      if (a$next_jump_time < a$end_time) {
        jso <- a$jump_search_order
        idx <- .linear_search(a$cur_rates[jso], a$sum_rate * stats::runif(1L))
        while (idx == 0L) {
          # the search finds no jump only if the sum exceeds the sum of the rates
          .agg_resum(a, a$cur_rates)
          idx <- .linear_search(a$cur_rates[jso], a$sum_rate * stats::runif(1L))
        }
        a$jump_search_idx <- idx
        a$next_jump <- jso[idx]
      }
    },
    NRM = , CCNRM = , Coevolve = {
      # ponytail: CCNRM's binned priority time table is a which.min scan here;
      # same law as Julia, O(J) per jump instead of O(1).
      if (a$J) {
        k <- which.min(a$pq)
        a$next_jump <- k
        a$next_jump_time <- a$pq[k]
      }
    },
    DirectCR = {
      if (.sum_rate_inaccurate(a$sum_rate, a$sum_rate_err, .pt_numterms(a$rt))) .pt_resum_agg(a, a$cur_rates)
      if (.agg_nomorejumps(a, a$sum_rate)) return(invisible())
      a$next_jump_time <- t + stats::rexp(1L) / a$sum_rate
      if (a$next_jump_time < a$end_time) {
        k <- .pt_sample(a$rt, a$cur_rates)
        if (k == 0L) {
          .pt_resum_agg(a, a$cur_rates)
          while (k == 0L) k <- .pt_sample(a$rt, a$cur_rates)
        }
        a$next_jump <- k
      }
    },
    RSSA = , RSSACR = .rssa_generate(a, u, p, t))
  invisible()
}

# register_next_jump_time!: a stop at the next jump time
.agg_register <- function(a, integ) {
  if (isTRUE(integ$is_ssa)) integ$tstop <- a$next_jump_time
  else if (a$next_jump_time < a$end_time) add_tstop(integ, a$next_jump_time)
  invisible()
}

# The affect of the aggregator's callback for ODE and SDE integrators
# (execute_tied_jumps!): the jump, then the next jumps tied with it in time.
# ponytail: tied jumps are not saved in between, ties have probability zero
# except at Coevolve rate interval ends.
.agg_execute_tied <- function(a, integ) {
  repeat {
    .agg_execute(a, integ)
    if (isTRUE(integ$terminate)) return(invisible())
    .agg_generate(a, integ)
    if (!(a$next_jump_time == integ$t && integ$t < a$end_time)) break
    if (a$type == "Coevolve" && !.coev_accept_tied(a, integ)) break
  }
  .agg_register(a, integ)
  invisible()
}

#' @rdname jump_aggregators
#' @export
reset_aggregated_jumps <- function(integ) {
  a <- integ$jump_agg
  if (!is.null(a)) {
    .agg_initialize(a, integ)
    .agg_register(a, integ)
  }
  invisible(NULL)
}

#' @rdname jump_aggregators
#' @export
add_tstop <- function(integ, t) {
  t <- as.numeric(t)
  # add_tstop! of OrdinaryDiffEqCore refuses a stop behind the current time
  tsp <- integ$prob$tspan
  tdir <- if (length(tsp) == 2L && tsp[2L] < tsp[1L]) -1 else 1
  if (length(integ$t) == 1L && any(tdir * (t - integ$t) < 0))
    stop("Tried to add a tstop that is behind the current time. This is strictly forbidden", call. = FALSE)
  integ$tstops_added <- c(integ$tstops_added, t)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# RSSA and RSSACR (rssa.jl, rssacr.jl, bracketing.jl)

.agg_jump_brackets <- function(a, k, u, p, t) {
  if (k <= a$nmaj) {
    rs <- a$js$maj_rs[[k]]
    return(c(.evalrxrate(a$ulow, rs, a$maj_rates[k]), .evalrxrate(a$uhigh, rs, a$maj_rates[k])))
  }
  i <- k - a$nmaj
  b <- a$js$crj_bounds[[i]]
  if (!is.null(b)) return(as.numeric(b(a$ulow, a$uhigh, u, p, t))[1:2])
  r <- a$js$crj_rate[[i]]
  rl <- r(a$ulow, p, t); rh <- r(a$uhigh, p, t)
  if (rl <= rh) c(rl, rh) else c(rh, rl)
}

# set_bracketing! and validate_brackets
.agg_set_bracketing <- function(a, u, p, t) {
  n <- length(u)
  a$ulow <- numeric(n); a$uhigh <- numeric(n)
  for (i in seq_len(n)) {
    b <- .spec_brackets(a$bd, i, u[i])
    a$ulow[i] <- b[1L]; a$uhigh[i] <- b[2L]
  }
  a$cur_rate_low <- numeric(a$J); a$cur_rate_high <- numeric(a$J)
  for (k in seq_len(a$J)) {
    b <- .agg_jump_brackets(a, k, u, p, t)
    a$cur_rate_low[k] <- b[1L]; a$cur_rate_high[k] <- b[2L]
  }
  a$sum_rate <- sum(a$cur_rate_high)
  for (k in seq_len(a$J)) {
    lr <- a$cur_rate_low[k]; ur <- a$cur_rate_high[k]
    if (!(0 <= lr && lr <= ur && ur < Inf))
      .stopf("Invalid rate bounds for jump %d: expected 0 <= lrate <= urate < Inf, got (%g, %g).", k, lr, ur)
    r <- .agg_rate(a, k, u, p, t)
    if (!(lr <= r && r <= ur))
      .stopf("Rate bounds for jump %d do not bracket the rate at the initial state: %g is not in [%g, %g].",
             k, r, lr, ur)
  }
}

.rssa_update_rates <- function(a, changed, u, p, t) {
  s <- a$sum_rate
  err <- a$sum_rate_err
  for (uidx in changed) {
    uval <- u[uidx]
    if (uval == 0 || uval < a$ulow[uidx] || uval > a$uhigh[uidx]) {
      b <- .spec_brackets(a$bd, uidx, uval)
      a$ulow[uidx] <- b[1L]; a$uhigh[uidx] <- b[2L]
      deps <- a$deps$vartojumps_map[[uidx]]
      before <- s
      for (j in deps) {
        old <- a$cur_rate_high[j]
        if (a$type == "RSSA") s <- s - old
        jb <- .agg_jump_brackets(a, j, u, p, t)
        a$cur_rate_low[j] <- jb[1L]; a$cur_rate_high[j] <- jb[2L]
        if (a$type == "RSSA") s <- s + jb[2L] else .pt_update(a$rt, j, old, jb[2L])
      }
      if (a$type == "RSSA") {
        err <- err + .update_error_bound(before, s, 2 * length(deps))
      } else {
        after <- a$rt$gsum
        err <- err + .update_error_bound(s, after, 4 * length(deps))
        s <- after
      }
    }
  }
  a$sum_rate_err <- err
  a$sum_rate <- s
}

# rejectrx: TRUE when candidate jump j is rejected
.rssa_reject <- function(a, j, u, p, t) {
  r2 <- stats::runif(1L) * a$cur_rate_high[j]
  lo <- a$cur_rate_low[j]
  if (lo > 0 && r2 <= lo) return(FALSE)
  rate <- .agg_rate(a, j, u, p, t)
  !(rate > 0 && r2 <= rate)
}

.rssa_sample <- function(a) {
  if (a$type == "RSSA") return(.linear_search(a$cur_rate_high, stats::runif(1L) * a$sum_rate))
  .pt_sample(a$rt, a$cur_rate_high)
}

# resum_and_search! / resum_and_sample!
.rssa_resample <- function(a) {
  if (a$type == "RSSA") .agg_resum(a, a$cur_rate_high) else .pt_resum_agg(a, a$cur_rate_high)
  j <- 0L
  while (j == 0L) j <- .rssa_sample(a)
  j
}

.rssa_generate <- function(a, u, p, t) {
  nterms <- if (a$type == "RSSA") a$J else .pt_numterms(a$rt)
  if (.sum_rate_inaccurate(a$sum_rate, a$sum_rate_err, nterms)) {
    if (a$type == "RSSA") .agg_resum(a, a$cur_rate_high) else .pt_resum_agg(a, a$cur_rate_high)
  }
  if (.agg_nomorejumps(a, a$sum_rate)) return(invisible())
  rerl <- 0
  j <- .rssa_sample(a)
  if (j == 0L) j <- .rssa_resample(a)
  rerl <- rerl + stats::rexp(1L)
  nrejected <- 0L
  while (.rssa_reject(a, j, u, p, t)) {
    # after a long run of rejections, check if a jump can still come before the end time
    nrejected <- nrejected + 1L
    if (nrejected == max(16L, a$J)) {
      nrejected <- 0L
      if (t + rerl / a$sum_rate >= a$end_time ||
          all(.agg_all_rates(a, u, p, t) <= 0)) {
        a$next_jump_time <- Inf
        return(invisible())
      }
    }
    j <- .rssa_sample(a)
    if (j == 0L) j <- .rssa_resample(a)
    rerl <- rerl + stats::rexp(1L)
  }
  a$next_jump <- j
  a$next_jump_time <- t + rerl / a$sum_rate
  invisible()
}

# ---------------------------------------------------------------------------
# PriorityTable (prioritytable.jl): group 1 holds the zero priorities, group
# 2 those below minpriority, and group g > 2 the priorities in
# [maxpriority_g / 2, maxpriority_g), with powers of two bounds.

# exponent(x) of Julia, floor(log2(x)) for positive finite x
.exponent <- function(x) {
  e <- floor(log2(x))
  if (2^e > x) e <- e - 1 else if (2^(e + 1) <= x) e <- e + 1
  as.integer(e)
}

.priortogid <- function(priority, mingid) {
  if (priority <= .Machine$double.eps) return(1L)
  gid <- .exponent(priority) + 1L
  if (gid <= mingid) return(2L)
  gid - mingid + 2L
}

.pt_new <- function(minexponent, priorities) {
  pt <- new.env(parent = emptyenv())
  minp <- 2^minexponent
  pt$mingid <- minexponent
  numgroups <- .priortogid(2 * minp, minexponent) - 1L
  pt$gmax <- c(0, minp * 2^(seq_len(numgroups - 1L) - 1L))
  pt$maxpriority <- pt$gmax[numgroups]
  pt$gnum <- integer(numgroups)
  pt$gpids <- rep(list(integer()), numgroups)
  pt$gsums <- numeric(numgroups)
  pt$gsum <- 0
  pt$pid_gid <- integer(length(priorities))
  pt$pid_idx <- integer(length(priorities))
  for (pid in seq_along(priorities)) .pt_insert(pt, pid, priorities[pid])
  pt
}

.pt_numterms <- function(pt) length(pt$pid_gid) + length(pt$gmax)

.pt_pad <- function(pt, priority) {
  while (priority >= pt$maxpriority) {
    pt$maxpriority <- 2 * pt$maxpriority
    pt$gmax <- c(pt$gmax, pt$maxpriority)
    pt$gnum <- c(pt$gnum, 0L)
    pt$gpids[[length(pt$gpids) + 1L]] <- integer()
    pt$gsums <- c(pt$gsums, 0)
  }
}

.pt_group_insert <- function(pt, gid, pid) {
  n <- pt$gnum[gid] + 1L
  pt$gnum[gid] <- n
  pt$gpids[[gid]][n] <- pid
  pt$pid_gid[pid] <- gid
  pt$pid_idx[pid] <- n
}

.pt_insert <- function(pt, pid, priority) {
  gid <- .priortogid(priority, pt$mingid)
  if (priority >= pt$maxpriority) .pt_pad(pt, priority)
  pt$gsum <- pt$gsum + priority
  pt$gsums[gid] <- pt$gsums[gid] + priority
  .pt_group_insert(pt, gid, pid)
}

.pt_update <- function(pt, pid, old, new) {
  oldgid <- .priortogid(old, pt$mingid)
  newgid <- .priortogid(new, pt$mingid)
  if (new >= pt$maxpriority) .pt_pad(pt, new)
  pt$gsum <- pt$gsum + (new - old)
  if (oldgid == newgid) {
    pt$gsums[newgid] <- pt$gsums[newgid] + (new - old)
    return(invisible())
  }
  # remove pid from its old group by swapping the last id into its place
  idx <- pt$pid_idx[pid]
  n <- pt$gnum[oldgid]
  moved <- pt$gpids[[oldgid]][n]
  pt$gpids[[oldgid]][idx] <- moved
  pt$gnum[oldgid] <- n - 1L
  pt$pid_gid[moved] <- oldgid
  pt$pid_idx[moved] <- idx
  .pt_group_insert(pt, newgid, pid)
  pt$gsums[oldgid] <- if (pt$gnum[oldgid] == 0L) 0 else pt$gsums[oldgid] - old
  pt$gsums[newgid] <- pt$gsums[newgid] + new
  invisible()
}

.pt_resum <- function(pt, priorities) {
  for (g in seq_along(pt$gsums))
    pt$gsums[g] <- sum(priorities[pt$gpids[[g]][seq_len(pt$gnum[g])]])
  pt$gsum <- sum(pt$gsums)
}

.pt_resum_agg <- function(a, priorities) {
  .pt_resum(a$rt, priorities)
  a$sum_rate <- a$rt$gsum
  a$sum_rate_err <- .sum_error_bound(a$sum_rate, .pt_numterms(a$rt))
}

# sample(pt, priorities): a group by its sum (largest first), then an id
# within it by rejection, or by weight in the group of zero priorities
.pt_sample <- function(pt, priorities) {
  if (pt$gsum <= 0) return(0L)
  r <- stats::runif(1L) * pt$gsum
  ng <- length(pt$gsums)
  hit <- which(cumsum(pt$gsums[ng:1L]) >= r)
  if (!length(hit)) return(0L)
  gid <- ng - hit[1L] + 1L
  pids <- pt$gpids[[gid]]
  n <- pt$gnum[gid]
  if (pt$gmax[gid] == 0) {
    pr <- priorities[pids[seq_len(n)]]
    total <- sum(pr)
    if (total == 0) return(0L)
    j <- .linear_search(pr, stats::runif(1L) * total)
    if (j == 0L) j <- max(which(pr > 0))
    return(pids[j])
  }
  repeat {
    r <- stats::runif(1L) * n
    i <- trunc(r)
    pid <- pids[i + 1L]
    if ((r - i) * pt$gmax[gid] < priorities[pid]) return(pid)
  }
}

# ---------------------------------------------------------------------------
# Coevolve (coevolve.jl). Jump k is a mass action jump (k <= nmaj), a
# constant rate jump (uidx = k - nmaj <= ncrj, always accepted) or a bounded
# variable rate jump (lidx = uidx - ncrj > 0).

.coev_urate <- function(a, k, u, p, t) {
  uidx <- k - a$nmaj
  if (uidx <= 0L) return(.evalrxrate(u, a$js$maj_rs[[k]], a$maj_rates[k]))
  if (uidx <= a$ncrj) return(a$js$crj_rate[[uidx]](u, p, t))
  a$js$bvrj[[uidx - a$ncrj]]$urate(u, p, t)
}

.coev_next_candidate_time <- function(a, u, p, t, s, lidx) {
  if (lidx <= 0L) return(t + s)
  j <- a$js$bvrj[[lidx]]
  ri <- j$rateinterval(u, p, t)
  if (s > ri) {
    # a rejection when the candidate lies past the interval
    a$cur_lrates[lidx] <- Inf
    return(t + ri)
  }
  t <- t + s
  a$cur_lrates[lidx] <- if (t < a$end_time) {
    if (is.null(j$lrate)) 0 else j$lrate(u, p, t)
  } else Inf
  t
}

.coev_next_time <- function(a, u, p, t, i) {
  lidx <- i - a$nmaj - a$ncrj
  urate <- .coev_urate(a, i, u, p, t)
  if (i != a$next_jump && lidx <= 0L) {
    last <- a$cur_rates[i]
    if (last > 0) {
      s <- if (urate == 0) Inf else last / urate * (a$pq[i] - t)
      return(c(.coev_next_candidate_time(a, u, p, t, s, lidx), urate))
    }
  }
  s <- if (urate == 0) Inf else stats::rexp(1L) / urate
  c(.coev_next_candidate_time(a, u, p, t, s, lidx), urate)
}

.coev_fill_rates_and_get_times <- function(a, u, p, t) {
  a$cur_rates <- numeric(a$J)
  a$cur_lrates <- numeric(a$nbvrj)
  a$pq <- rep(Inf, a$J)
  for (i in seq_len(a$J)) {
    nt <- .coev_next_time(a, u, p, t, i)
    a$pq[i] <- nt[1L]
    a$cur_rates[i] <- nt[2L]
  }
}

# accept_next_jump!: thinning of a bounded variable rate jump; a rejected
# candidate gets a new time and the next jump is chosen again.
.coev_accept <- function(a, integ) {
  k <- a$next_jump
  lidx <- k - a$nmaj - a$ncrj
  if (lidx <= 0L) return(TRUE)
  u <- .agg_u(a, integ); p <- integ$p; t <- integ$t
  j <- a$js$bvrj[[lidx]]
  urate <- a$cur_rates[k]
  lrate <- a$cur_lrates[lidx]
  s <- -1
  if (lrate == Inf) {
    urate <- j$urate(u, p, t)
    s <- if (urate == 0) Inf else stats::rexp(1L) / urate
  } else if (lrate < urate) {
    v <- stats::runif(1L) * urate
    if (v > lrate) {
      rate <- j$rate(u, p, t)
      if (v > rate) {
        urate <- j$urate(u, p, t)
        s <- if (urate == 0) Inf else stats::rexp(1L) / urate
      }
    }
  }
  if (s < 0) return(TRUE)
  a$pq[k] <- .coev_next_candidate_time(a, u, p, t, s, lidx)
  a$cur_rates[k] <- urate
  a$prev_jump <- k
  .agg_generate(a, integ)
  .agg_register(a, integ)
  integ$u_modified <- FALSE
  FALSE
}

# candidates tied with the current time, after a jump
.coev_accept_tied <- function(a, integ) {
  um <- integ$u_modified
  tied <- FALSE
  while (a$next_jump_time == integ$t && integ$t < a$end_time) {
    if (.coev_accept(a, integ)) {
      tied <- TRUE
      break
    }
  }
  integ$u_modified <- um
  tied
}
