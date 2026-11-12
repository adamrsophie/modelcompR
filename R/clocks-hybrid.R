# Sampled data (hybrid) systems, after split_system of
# StateSelection.jl/lib/ModelingToolkitTearing/src/clock_inference/
# clock_inference.jl, the hybrid branch of mtkcompile! in
# ModelingToolkit.jl/src/systems/systemstructure.jl and
# docs/src/tutorials/SampledData.md. Clock inference itself is in
# structural-clocks.R.
#
# ModelingToolkit 11 rejects hybrid systems unless an external discrete
# compilation pass is supplied. modelcompR supplies one, as ModelingToolkit 9
# did with generate_discrete_affect: the continuous partition is compiled as
# usual, with every discrete variable (and its past values) as a parameter,
# and each discrete partition becomes a callback that, at every tick of its
# clock, reads the sampled continuous values, computes the new discrete values
# and shifts the history. This contract is an R extension: the update runs at
# every tick from the initial time (when the clock ticks there) and never at
# the final time; implicit partitions are solved with Newton's method; a
# value given for a discrete variable x is its value before the first tick,
# which is also x(k - 1) at that tick unless (k - 1)(x) is given.

.has_clock_ops <- function(e) hasshift(e) || hassample(e) || hashold(e) || .has_op(e, "SampleTime")

# Wraps the compile function of mtkcompile(): systems without clock operators
# compile as before, purely discrete ones through the shift branch, hybrid
# ones through .compile_hybrid().
.compile_clocked <- function(compile_fn) {
  force(compile_fn)
  function(sys) {
    eqs <- Filter(function(e) inherits(e, "Equation"), .f(sys, "eqs"))
    machines <- length(.f(sys, "metadata")$mtk_statemachines) > 0L
    if (is.null(.f(sys, "iv")) || ((!length(eqs) || !any(vapply(eqs, .has_clock_ops, TRUE))) && !machines))
      return(compile_fn(sys))
    hybrid <- machines || any(vapply(eqs, function(e) hassample(e) || hashold(e) || is_diff_equation(e), TRUE))
    if (!hybrid) return(compile_fn(.update_system(sys, eqs = .substitute_sample_time(eqs, NULL),
                                                  is_discrete = TRUE)))
    .compile_hybrid(sys, eqs, compile_fn)
  }
}

# SampleTime() is the sample interval of the clock of its equation
# (substitute_sample_time); the clock of a purely discrete system is the one
# of its shifts, and integer steps have interval 1.
.substitute_sample_time <- function(eqs, clock_id) {
  nodes <- .find_ops(.eq_sides(eqs), "SampleTime")
  if (!length(nodes)) return(eqs)
  if (is.null(clock_id)) {
    ids <- unique(vapply(.shift_nodes(.eq_sides(eqs)), function(n) .shift_parts(n)$clock, 0L))
    clock_id <- if (length(ids[ids > 0L]) == 1L) ids[ids > 0L] else 0L
  }
  to <- lapply(nodes, function(n) {
    cid <- as.integer(as.numeric(.args(n)[[1L]]))
    clk <- .clock_of(if (cid > 0L) cid else clock_id)
    if (is.null(clk) || identical(clk$type, "integer")) return(.S1())
    if (is.null(clk$dt)) stop("SampleTime() needs a clock with a known sample interval", call. = FALSE)
    as_sym(clk$dt)
  })
  rule <- .rule(nodes, to)
  lapply(eqs, function(e) Equation(.subst(e$lhs, rule), .subst(e$rhs, rule)))
}

# ---------------------------------------------------------------------------
# Partitioning (split_system) and the discrete compile pass

.strip_sample <- function(e) {
  nodes <- .find_ops(list(e$lhs, e$rhs), "Sample")
  if (!length(nodes)) return(e)
  rule <- .rule(nodes, lapply(nodes, function(n) .args(n)[[1L]]))
  Equation(.subst(e$lhs, rule), .subst(e$rhs, rule))
}

.strip_hold <- function(e) {
  nodes <- .find_ops(list(e$lhs, e$rhs), "Hold")
  if (!length(nodes)) return(e)
  rule <- .rule(nodes, lapply(nodes, function(n) .args(n)[[1L]]))
  Equation(.subst(e$lhs, rule), .subst(e$rhs, rule))
}

.compile_hybrid <- function(sys, eqs, compile_fn) {
  iv <- .f(sys, "iv")
  us <- .f(sys, "unknowns")
  # state machines (structural-statemachines.R) take part in clock inference
  # through the equations of their states, their transition conditions and
  # an equation tying all their variables to one clock
  machines <- .f(sys, "metadata")$mtk_statemachines %||% list()
  sm_eqs <- list(); sm_of <- integer(); sm_def <- logical()
  for (i in seq_along(machines)) {
    es <- unlist(machines[[i]]$state_eqs, recursive = FALSE)
    defined <- .unique_syms(Filter(.is_declared, lapply(es, `[[`, "lhs")))
    new <- c(es, lapply(machines[[i]]$transitions, function(tr) Equation(.S0(), tr$cond)),
             if (length(defined) > 1L) list(Equation(.S0(), Reduce(`+`, defined))))
    sm_eqs <- c(sm_eqs, new)
    sm_of <- c(sm_of, rep(i, length(new)))
    sm_def <- c(sm_def, rep(TRUE, length(es)), rep(FALSE, length(new) - length(es)))
  }
  inf <- .ci_infer_clocks(c(eqs, sm_eqs), .keys(us), .f(sys, "initialization_eqs"))
  sm_dom <- inf$eq[length(eqs) + seq_along(sm_eqs)]
  inf$eq <- inf$eq[seq_along(eqs)]
  vdom <- vapply(us, inf$var, 0L)
  disc_vars <- us[vdom > 0L]
  parts <- list()
  hist_all <- list()
  for (cid in sort(unique(inf$eq[inf$eq > 0L]))) {
    clk <- .clock_of(cid)
    if (identical(clk$type, "integer"))
      stop("a hybrid system cannot use ShiftIndex(t); give its discrete partitions a Clock()", call. = FALSE)
    if (identical(clk$type, "periodic") && is.null(clk$dt))
      stop("the clock of a discrete partition needs a sample interval", call. = FALSE)
    peqs <- lapply(.substitute_sample_time(eqs[inf$eq == cid], cid), .strip_sample)
    pvars <- us[vdom == cid]
    low <- .lower_shifts(peqs, pvars, iv)
    parts[[length(parts) + 1L]] <- c(list(clock = clk, vars = pvars, hist = low$hist),
                                     .explicit_partition(low$eqs, pvars))
    for (h in low$hist) hist_all <- c(hist_all, h$vars)
  }
  if (length(machines)) {
    ids <- vapply(seq_along(machines), function(i) {
      d <- unique(sm_dom[sm_of == i & sm_def])
      if (length(d) != 1L || d == 0L)
        stop("the states of state machine ", machines[[i]]$name, " must be on one discrete clock", call. = FALSE)
      d
    }, 0L)
    smp <- .sm_partitions(machines, ids, iv)
    parts <- c(parts, smp$parts)
    hist_all <- c(hist_all, smp$hist)
    disc_vars <- .unique_syms(c(disc_vars, smp$vars))
  }
  if (length(.find_ops(.eq_sides(eqs[inf$eq == 0L]), "Sample")))
    stop("Sample() in a continuous equation", call. = FALSE)
  ceqs <- lapply(eqs[inf$eq == 0L], .strip_hold)
  # discrete variables and their past values are parameters of the
  # continuous system, with values before the first tick
  binds <- .f(sys, "bindings"); ics <- .f(sys, "initial_conditions")
  for (v in disc_vars) {
    b <- smap_get(binds, v)
    if (!is.null(b)) {
      binds <- smap_del(binds, v)
      if (!smap_has(ics, v)) ics <- smap_set(ics, v, b)
    }
    if (!smap_has(ics, v)) ics <- smap_set(ics, v, 0)
  }
  for (h in hist_all) if (!smap_has(ics, h)) ics <- smap_set(ics, h, .meta_get(h, "unshifted"))
  csys <- .update_system(sys, eqs = ceqs, unknowns = us[vdom == 0L],
                         ps = .unique_syms(c(.f(sys, "ps"), disc_vars, hist_all)),
                         bindings = binds, initial_conditions = ics)
  newsys <- compile_fn(csys)
  md <- .f(newsys, "metadata") %||% list()
  md$mtk_clock_partitions <- parts
  .update_system(newsys, metadata = md)
}

# Explicit assignments of a discrete partition: an equation defines its
# left-hand side variable, or its right-hand side when that is a lone
# variable (connections such as y ~ r); the rest are residuals solved by
# Newton's method for the variables left undetermined.
.explicit_partition <- function(eqs, pvars) {
  vk <- .keys(pvars)
  classify <- function(allow) {
    assigned <- character(); asg <- list(); res <- list()
    for (i in seq_along(eqs)) {
      e <- eqs[[i]]
      lk <- if (.is_declared(e$lhs)) .key(e$lhs) else ""
      rk <- if (.is_declared(e$rhs)) .key(e$rhs) else ""
      if (allow[i] && nzchar(lk) && lk %in% vk && !(lk %in% assigned) && !(lk %in% .keys(.collect_vars(e$rhs)))) {
        asg[[length(asg) + 1L]] <- e
        assigned <- c(assigned, lk)
      } else if (allow[i] && nzchar(rk) && rk %in% vk && !(rk %in% assigned) &&
                 !(rk %in% .keys(.collect_vars(e$lhs)))) {
        asg[[length(asg) + 1L]] <- Equation(e$rhs, e$lhs)
        assigned <- c(assigned, rk)
      } else {
        res[[length(res) + 1L]] <- e$rhs - e$lhs
      }
    }
    list(asg = asg, res = res, assigned = assigned)
  }
  cl <- classify(rep(TRUE, length(eqs)))
  if (inherits(tryCatch(.toposort_observed(cl$asg), error = function(e) e), "error"))
    cl <- classify(rep(FALSE, length(eqs)))
  z <- Filter(function(v) !(.key(v) %in% cl$assigned), pvars)
  if (length(cl$res) != length(z))
    stop("a discrete partition has ", length(cl$res), " implicit equations for ", length(z),
         " undetermined variables", call. = FALSE)
  list(assigns = .toposort_observed(cl$asg), residuals = cl$res, newton = z)
}

# Callbacks of the discrete partitions, built when a problem is made.
.clock_partition_callbacks <- function(parts, ctx) {
  lapply(parts, function(pt) {
    if (identical(pt$type, "statemachine")) return(.sm_callback(pt, ctx))
    hist_new <- list(); hist_old <- list()
    for (h in pt$hist) {
      n <- length(h$vars)
      hist_old <- c(hist_old, h$vars)
      hist_new <- c(hist_new, list(h$var), if (n > 1L) h$vars[-n])
    }
    for (v in c(pt$vars, hist_old))
      if (!(.key(v) %in% ctx$pk)) stop("discrete variable ", format(v), " is not a parameter of the system",
                                     call. = FALSE)
    args <- list(list(name = "u", vars = ctx$us), list(name = "p", vars = ctx$ps),
                 list(name = "z", vars = pt$newton), list(name = "t", vars = ctx$iv))
    obs <- c(pt$assigns, ctx$obs)
    val_fn <- .build_closure(c(pt$vars, hist_new), args, observed = obs, out = "vector")
    res_fn <- if (length(pt$residuals)) .build_closure(pt$residuals, args, observed = obs, out = "vector")
    z_idx <- match(.keys(pt$newton), ctx$pk)
    idx <- match(.keys(c(pt$vars, hist_old)), ctx$pk)
    aff <- function(integ) {
      u <- integ$u; p <- integ$p; t <- integ$t
      z <- numeric()
      if (!is.null(res_fn)) {
        sol <- .newton_solve(function(zz) res_fn(u, p, zz, t), p[z_idx])
        if (!sol$converged) stop("the discrete partition could not be solved at t = ", t, call. = FALSE)
        z <- sol$x
      }
      integ$p[idx] <- val_fn(u, p, z, t)
      invisible(NULL)
    }
    clk <- pt$clock
    if (identical(clk$type, "solverstep")) {
      DiscreteCallback(function(u, t, integ) TRUE, aff, initialize = aff)
    } else {
      .periodic_event_callback(aff, clk$dt, clk$phase)
    }
  })
}

.mtk$hooks$compile_clocked <- function(compile_fn) .compile_clocked(compile_fn)
