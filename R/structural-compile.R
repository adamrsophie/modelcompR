# The structural simplification pipeline of mtkcompile, after
# ModelingToolkit.jl v11.47 (MIT): src/systems/systems.jl (__mtkcompile,
# map_variables_to_equations), src/systems/systemstructure.jl (mtkcompile!,
# _mtkcompile!, _mtkcompile_worker!, inputs_to_parameters!, validate_io!),
# src/systems/alias_elimination.jl (alias_elimination!, union_with_sign!,
# pick_alias_target, find_perfect_aliases!, eliminate_zero_variables!),
# src/structural_transformation/{pantelides,symbolics_tearing,utils}.jl
# (pantelides_reassemble, dae_index_lowering, tearing, dummy_derivative,
# sorted_incidence_matrix, but_ordered_incidence, variable masks),
# src/systems/substitute_component.jl and isolate_subsystem from
# src/systems/analysis_points.jl.
#
# The graph algorithms are in structural-graph.R, the StateSelection passes in
# structural-stateselection.R and the tearing state in structural-tearing.R.
# Clock inference is in structural-clocks.R and clocks-hybrid.R, state
# machines in structural-statemachines.R, IfLifting in
# structural-iflifting.R, inline solves of torn blocks in
# structural-inline.R, array arithmetic in structural-arrays.R and the stage
# system of solver_nlprob.jl in structural-nlstep.R.

# ---------------------------------------------------------------------------
# Inputs

# inputs_to_parameters!: inputs leave the unknowns and become parameters.
# As in .base_compile(), their defaults are kept.
.st_inputs_to_parameters <- function(ts, inputs, outputs) {
  sys <- ts$sys
  if (length(inputs)) {
    ik <- .keys(inputs)
    .ss_rm_eqs_vars(ts, integer(), which(.keys(ts$fullvars) %in% ik))
    sys <- .update_system(sys, unknowns = Filter(function(u) !(.key(u) %in% ik), .f(sys, "unknowns")),
                          ps = .unique_syms(c(.f(sys, "ps"), inputs)))
  }
  ts$sys <- .update_system(sys, inputs = inputs, outputs = outputs)
  invisible(ts)
}

# validate_io!
.st_validate_io <- function(ts, what, vars) {
  miss <- vars[!(.keys(vars) %in% .keys(ts$fullvars))]
  if (length(miss))
    .stopf("The %s %s of system %s do not occur in its equations.", what,
           paste(vapply(miss, .fmt, ""), collapse = ", "), nameof(ts$sys))
  invisible()
}

# ---------------------------------------------------------------------------
# Substitution with incidence update

# __substitute_and_update_incidence!: substitution never introduces
# variables, so the new incidence of an equation is the subset of the old one
# still present in the simplified expression.
.st_substitute_and_update <- function(ts, eq_idxs, rule) {
  s <- ts$structure
  for (e in unique(eq_idxs)) {
    ts$eqs[[e]] <- .subst(ts$eqs[[e]], rule)
    ts$original_eqs[[e]] <- .subst(ts$original_eqs[[e]], rule)
    present <- .keys(.ts_search_vars(ts$eqs[[e]]$rhs - ts$eqs[[e]]$lhs))
    nb <- s$graph$fadjlist[[e]]
    row <- nb[.keys(ts$fullvars[nb]) %in% present]
    .bg_set_neighbors(s$graph, e, row)
    if (!is.null(s$solvable_graph))
      .bg_set_neighbors(s$solvable_graph, e, intersect(row, s$solvable_graph$fadjlist[[e]]))
  }
  invisible(ts)
}

# ---------------------------------------------------------------------------
# Perfect aliases

# pick_alias_target: an irreducible member if there is one, else the member
# with the highest state priority; ties with a positive priority go to the
# variable in the most equations, then to the canonical order.
.st_pick_alias_target <- function(ts, group) {
  s <- ts$structure
  irr <- Filter(function(v) .ts_is_irreducible(ts, ts$fullvars[[v]]), group)
  if (length(irr)) return(irr[1L])
  sp <- s$state_priorities[group]
  cand <- group[sp == max(sp)]
  if (length(cand) > 1L && max(sp) > 0L) {
    deg <- lengths(s$graph$badjlist[cand])
    cand <- cand[deg == max(deg)]
    cand <- cand[order(s$canonical_ranks[cand])]
  }
  cand[1L]
}

# eliminate_perfect_aliases! (find_perfect_aliases!): equations v1 ~ v2 and
# v1 ~ -v2 between undifferentiated variables define alias groups. One
# variable per group survives; the others, and their derivatives, are
# substituted and become observed. A group with a sign contradiction
# (x ~ y and x ~ -y) is zero. Returns the index maps and the alias of each
# removed variable.
.st_eliminate_perfect_aliases <- function(ts) {
  s <- ts$structure; g <- s$graph
  nv <- length(ts$fullvars)
  # union find with the sign of each variable relative to its root
  parent <- integer(nv); parity <- integer(nv); members <- vector("list", nv)
  union_with_sign <- function(v1, v2, edge_sign) {
    for (v in c(v1, v2)) if (parent[v] == 0L) {
      parent[v] <<- v; parity[v] <<- 1L; members[[v]] <<- v
    }
    r1 <- parent[v1]; s1 <- parity[v1]; r2 <- parent[v2]; s2 <- parity[v2]
    if (r1 == r2) {
      if (s1 != edge_sign * s2) parity[members[[r1]]] <<- 0L
      return(invisible())
    }
    if (length(members[[r1]]) < length(members[[r2]])) {
      tmp <- r1; r1 <- r2; r2 <- tmp
      tmp <- s1; s1 <- s2; s2 <- tmp
    }
    r2_to_r1 <- s1 * edge_sign * s2
    mem <- members[[r2]]
    parent[mem] <<- r1
    parity[mem] <<- parity[mem] * r2_to_r1
    if (r2_to_r1 == 0L) parity[members[[r1]]] <<- 0L
    members[[r1]] <<- c(members[[r1]], mem)
    members[r2] <<- list(NULL)
    invisible()
  }
  cands <- list()
  for (ieq in seq_len(.bg_nsrcs(g))) {
    nb <- g$fadjlist[[ieq]]
    if (length(nb) != 2L || any(!is.na(s$diff_to_var[nb]))) next
    eq <- ts$eqs[[ieq]]
    ph <- list(Symbol(".sv1"), Symbol(".sv2"))
    low <- .subst(eq$rhs - eq$lhs, .rule(ts$fullvars[nb], ph))
    a1 <- diff(low, ph[[1L]]); a2 <- diff(low, ph[[2L]])
    if (!.is_num(a1) || !.is_num(a2)) next
    if (!.is_zero(.subst(low, .rule(ph, list(.S0(), .S0()))))) next
    a1 <- as.numeric(a1); a2 <- as.numeric(a2)
    edge_sign <- if (a1 + a2 == 0) 1L else if (a1 - a2 == 0) -1L else next
    cands[[length(cands) + 1L]] <- c(ieq, nb)
    union_with_sign(nb[1L], nb[2L], edge_sign)
  }
  eqs_rm <- integer(); vars_rm <- integer(); eqs_sub <- integer()
  sub_from <- list(); sub_to <- list()
  alias_target <- integer(); alias_sign <- numeric()
  if (length(cands)) {
    is_irr <- function(v) .ts_is_irreducible(ts, ts$fullvars[[v]])
    signed <- function(sg, x) if (sg == 1L) x else -x
    add_sub <- function(v, to) {
      sub_from[[length(sub_from) + 1L]] <<- ts$fullvars[[v]]
      sub_to[[length(sub_to) + 1L]] <<- to
    }
    move_edges <- function(v, target) {
      for (e in g$badjlist[[v]]) {
        eqs_sub <<- c(eqs_sub, e)
        .bg_rem_edge(g, e, v)
        .bg_add_edge(g, e, target)
        sg <- s$solvable_graph
        if (!is.null(sg) && .bg_has_edge(sg, e, v)) {
          .bg_rem_edge(sg, e, v)
          .bg_add_edge(sg, e, target)
        }
      }
    }
    clear_var <- function(v) {
      eqs_sub <<- c(eqs_sub, g$badjlist[[v]])
      .bg_set_dneighbors(g, v, integer())
      if (!is.null(s$solvable_graph)) .bg_set_dneighbors(s$solvable_graph, v, integer())
    }
    group_target <- integer(nv)
    irrs_by_root <- vector("list", nv)
    for (root in which(lengths(members) > 0L)) {
      group <- members[[root]]
      if (parity[root] == 0L) {
        # contradiction: every variable is zero, only irreducibles survive
        irrs_by_root[[root]] <- Filter(is_irr, group)
        for (v in group) {
          if (is_irr(v)) {
            ts$always_present[v] <- TRUE
            next
          }
          ts$additional_observed[[length(ts$additional_observed) + 1L]] <- Equation(ts$fullvars[[v]], .S0())
          dv <- v
          while (!is.na(dv)) {
            vars_rm <- c(vars_rm, dv)
            add_sub(dv, .S0())
            clear_var(dv)
            dv <- s$var_to_diff[dv]
          }
        }
        next
      }
      target <- .st_pick_alias_target(ts, group)
      group_target[root] <- target
      for (v in group) {
        if (is_irr(v) || v == target) {
          ts$always_present[v] <- TRUE
          next
        }
        sgn <- parity[v] * parity[target]
        rhs <- signed(sgn, ts$fullvars[[target]])
        vars_rm <- c(vars_rm, v)
        add_sub(v, rhs)
        ts$additional_observed[[length(ts$additional_observed) + 1L]] <- Equation(ts$fullvars[[v]], rhs)
        alias_target[v] <- target; alias_sign[v] <- sgn
        move_edges(v, target)
        # the derivatives of an alias are aliases of the target's derivatives
        dv <- s$var_to_diff[v]
        prev <- target
        dtarget <- s$var_to_diff[target]
        while (!is.na(dv)) {
          if (is.na(dtarget)) dtarget <- .ts_var_derivative(ts, prev)
          vars_rm <- c(vars_rm, dv)
          add_sub(dv, signed(sgn, ts$fullvars[[dtarget]]))
          alias_target[dv] <- dtarget; alias_sign[dv] <- sgn
          move_edges(dv, dtarget)
          dv <- s$var_to_diff[dv]
          prev <- dtarget
          dtarget <- s$var_to_diff[dtarget]
        }
      }
    }
    for (cd in cands) {
      ieq <- cd[1L]; v1 <- cd[2L]; v2 <- cd[3L]
      if (parity[v1] == 0L) {
        irrs <- irrs_by_root[[parent[v1]]]
        if (!length(irrs)) eqs_rm <- c(eqs_rm, ieq)
        else {
          # pin one irreducible of the contradictory group to zero
          v_pin <- irrs[length(irrs)]
          irrs_by_root[[parent[v1]]] <- irrs[-length(irrs)]
          .bg_set_neighbors(g, ieq, v_pin)
          if (!is.null(s$solvable_graph)) .bg_set_neighbors(s$solvable_graph, ieq, v_pin)
          ts$eqs[[ieq]] <- Equation(.S0(), ts$fullvars[[v_pin]])
          ts$original_eqs[[ieq]] <- Equation(ts$fullvars[[v_pin]], .S0())
        }
      } else {
        target <- group_target[parent[v1]]
        c1 <- if (is_irr(v1)) v1 else target
        c2 <- if (is_irr(v2)) v2 else target
        if (c1 == c2) eqs_rm <- c(eqs_rm, ieq)
      }
    }
    .st_substitute_and_update(ts, eqs_sub, .rule(sub_from, sub_to))
    # alias equations that became copies of one another
    seen <- character()
    for (cd in cands) {
      ieq <- cd[1L]
      if (ieq %in% eqs_rm) next
      nb <- g$fadjlist[[ieq]]
      if (length(nb) != 2L) next
      pair <- paste(nb, collapse = "-")
      if (pair %in% seen) eqs_rm <- c(eqs_rm, ieq) else seen <- c(seen, pair)
    }
  }
  n <- length(ts$fullvars)
  target <- integer(n); sgn <- rep(1, n)
  hit <- which(!is.na(alias_target[seq_len(n)]) & alias_target[seq_len(n)] > 0L)
  target[hit] <- alias_target[hit]; sgn[hit] <- alias_sign[hit]
  o2n <- .ss_rm_eqs_vars(ts, eqs_rm, vars_rm)
  list(eq = o2n$eq, var = o2n$var, alias_target = target, alias_sign = sgn)
}

# ---------------------------------------------------------------------------
# Alias elimination on the integer linear subsystem

# alias_elimination!: Bareiss on the equations that are integer linear
# combinations of variables. Redundant equations (rows that reduce to 0 ~ 0)
# are removed and the others are replaced by their reduced rows.
.st_alias_elimination <- function(ts) {
  mm <- .ts_linear_subsys_adjmat(ts)
  r <- .structural_singularity_removal(ts, mm)
  mm <- r$mm
  # Julia reports these through its verbosity system; they are kept for
  # inspection and show up in the consistency check when they matter
  ts$underconstrained <- ts$fullvars[r$underconstrained]
  fk <- .keys(ts$fullvars)
  eqs_rm <- integer()
  for (i in seq_along(mm$nzrows)) {
    eq <- mm$nzrows[i]; rc <- mm$row_cols[[i]]; rv <- mm$row_vals[[i]]
    if (!length(rc)) {
      eqs_rm <- c(eqs_rm, eq)
      next
    }
    rhs <- .ts_linear_combination(ts, rv, rc)
    ts$eqs[[eq]] <- Equation(.S0(), rhs)
    oeq <- ts$original_eqs[[eq]]
    col <- match(match(.key(oeq$lhs), fk), rc)
    ts$original_eqs[[eq]] <- if (!is.na(col) && rv[col] == -1) Equation(oeq$lhs, rhs + oeq$lhs) else ts$eqs[[eq]]
  }
  o2n <- .ss_rm_eqs_vars(ts, eqs_rm, integer())
  .ss_get_new_mm(mm, o2n$eq, o2n$var)
}

# ---------------------------------------------------------------------------
# Zero variables

# __eliminate_zero_variables!: a row of the linear subsystem with a single
# variable makes that variable, and all its derivatives, zero.
.st_find_zero_variables <- function(ts, mm) {
  s <- ts$structure; g <- s$graph
  zero_vars <- integer(); zero_eqs <- integer(); eqs_sub <- integer()
  zero_var_eq <- integer()
  queue <- which(lengths(mm$row_cols) == 1L)
  process_neighbors <- function(zvar) {
    nb <- g$badjlist[[zvar]]
    eqs_sub <<- union(eqs_sub, nb)
    rows <- match(nb, mm$nzrows)
    queue <<- c(queue, rows[!is.na(rows)])
  }
  while (length(queue)) {
    row_i <- queue[1L]; queue <- queue[-1L]
    eqidx <- mm$nzrows[row_i]
    if (eqidx %in% zero_eqs) next
    keep <- !(mm$row_cols[[row_i]] %in% zero_vars)
    mm$row_cols[[row_i]] <- mm$row_cols[[row_i]][keep]
    mm$row_vals[[row_i]] <- mm$row_vals[[row_i]][keep]
    if (sum(keep) != 1L) next
    zvar <- mm$row_cols[[row_i]][1L]
    zero_vars <- c(zero_vars, zvar); zero_eqs <- c(zero_eqs, eqidx)
    zero_var_eq[zvar] <- eqidx
    process_neighbors(zvar)
    dz <- s$var_to_diff[zvar]
    while (!is.na(dz)) {
      zero_vars <- c(zero_vars, dz)
      process_neighbors(dz)
      dz <- s$var_to_diff[dz]
    }
  }
  zero_vars <- unique(zero_vars)
  # a zero variable whose chain starts at an irreducible variable is kept
  protected <- integer()
  if (length(ts$irreducible_keys)) for (v in zero_vars) {
    root <- v
    while (!is.na(lower <- s$diff_to_var[root])) root <- lower
    if (.ts_is_irreducible(ts, ts$fullvars[[root]])) protected <- c(protected, v)
  }
  protected_eqs <- zero_var_eq[protected]
  protected_eqs <- protected_eqs[!is.na(protected_eqs)]
  # Remove the `var ~ 0` rows, as the Julia docstring says. (Julia's
  # __filter_mm_eqs! is called with the negated predicate and keeps only those
  # rows; the matrix is a cache, so only later reuse of rows differs.)
  drop <- mm$nzrows %in% setdiff(zero_eqs, protected_eqs)
  mm$nzrows <- mm$nzrows[!drop]; mm$row_cols <- mm$row_cols[!drop]; mm$row_vals <- mm$row_vals[!drop]
  list(mm = mm, zero_vars = zero_vars, zero_eqs = zero_eqs, eqs_sub = eqs_sub,
       protected_vars = protected, protected_eqs = protected_eqs)
}

# eliminate_zero_variables!: substitute the zero variables everywhere. A
# variable whose derivative is zero is integrated analytically: x ~ x#0 for
# D(x) ~ 0, x ~ xˍt#0 * t + x#0 for D(D(x)) ~ 0, with new parameters for the
# values at t = 0 that initialization solves for.
# Departure from Julia: each level is the integral of the one above,
# sum(p_i * t^(m - i) / (m - i)!). Julia's `rhs * t + p` is the same up to
# second order but for D(D(D(x))) ~ 0 gives x = a t^2 + b t + c next to
# D(x) = a t + b, which contradict each other.
# The new parameters are plain symbols; Julia declares them as functions of
# time only for the benefit of callbacks (see the FIXME in
# __create_t0_parameter_for).
.st_eliminate_zero_variables <- function(ts, mm) {
  z <- .st_find_zero_variables(ts, mm)
  mm <- z$mm
  if (!length(z$zero_vars)) return(list(mm = mm, modified = FALSE))
  unprotected <- setdiff(z$zero_vars, z$protected_vars)
  eqs_sub <- setdiff(z$eqs_sub, z$protected_eqs)
  if (!length(unprotected) && !length(eqs_sub)) return(list(mm = mm, modified = FALSE))
  s <- ts$structure; g <- s$graph
  iv <- ts$iv
  sys <- ts$sys
  sub_from <- list(); sub_to <- list()
  vars_rm <- integer(); mm_eqs_rm <- integer()
  new_ps <- list(); integ <- .f(sys, "analytically_integrated") %||% list()
  guesses <- .f(sys, "guesses"); binds <- .f(sys, "bindings")
  observe <- function(lhs, rhs) ts$additional_observed[[length(ts$additional_observed) + 1L]] <- Equation(lhs, rhs)
  analytic <- function(dvar, rhs) {
    ts$ad_from[[length(ts$ad_from) + 1L]] <- dvar
    ts$ad_to[[length(ts$ad_to) + 1L]] <- rhs
  }
  for (v in z$zero_vars) {
    sym <- ts$fullvars[[v]]
    ttsym <- default_toterm(sym)
    sub_from <- c(sub_from, list(sym, ttsym)); sub_to <- c(sub_to, list(.S0(), .S0()))
    if (v %in% z$protected_vars) next
    observe(ttsym, .S0())
    lower <- s$diff_to_var[v]
    if (is.na(lower)) next
    if (lower %in% z$zero_vars) {
      analytic(.st_D(default_toterm(ts$fullvars[[lower]]), iv), .S0())
      next
    }
    rhs <- .S0()
    consts <- list()
    while (!is.na(lower)) {
      vars_rm <- c(vars_rm, lower)
      sym <- ts$fullvars[[lower]]
      ttsym <- default_toterm(sym)
      analytic(.st_D(ttsym, iv), rhs)
      ts$pdm_from <- c(ts$pdm_from, list(.st_D(sym, iv), .st_D(ttsym, iv)))
      ts$pdm_to <- c(ts$pdm_to, list(rhs, rhs))
      t0 <- .make_var(paste0(.name(ttsym), "#0"), "parameter", NULL,
                      list(description = paste0("value of ", format(ttsym), " at ", format(iv), " = 0")))
      integ <- smap_set(integ, sym, t0)
      consts[[length(consts) + 1L]] <- t0
      m <- length(consts)
      rhs <- Reduce(`+`, lapply(seq_len(m), function(i) consts[[i]] * iv^(m - i) / factorial(m - i)))
      new_ps[[length(new_ps) + 1L]] <- t0
      observe(ttsym, rhs)
      guesses <- smap_set(guesses %||% symmap(), t0, sym)
      binds <- smap_set(binds %||% symmap(), t0, NA)
      mm_eqs_rm <- union(mm_eqs_rm, g$badjlist[[lower]])
      lower <- s$diff_to_var[lower]
    }
  }
  drop <- mm$nzrows %in% mm_eqs_rm
  mm$nzrows <- mm$nzrows[!drop]; mm$row_cols <- mm$row_cols[!drop]; mm$row_vals <- mm$row_vals[!drop]
  vars_rm <- c(vars_rm, unprotected)
  eqs_rm <- setdiff(z$zero_eqs, z$protected_eqs)
  .st_substitute_and_update(ts, eqs_sub, .rule(sub_from, sub_to))
  # a higher derivative may have dropped out with its zero cofactor
  for (v in seq_len(.bg_ndsts(g))) if (!.ss_is_present(s, v)) vars_rm <- c(vars_rm, v)
  ts$sys <- .update_system(sys, analytically_integrated = integ,
                           ps = .unique_syms(c(.f(sys, "ps"), new_ps)), bindings = binds, guesses = guesses)
  o2n <- .ss_rm_eqs_vars(ts, eqs_rm, vars_rm)
  mm <- .ss_get_new_mm(mm, o2n$eq, o2n$var)
  # substituted equations may have become integer linear combinations
  for (e in eqs_sub) {
    e2 <- o2n$eq[e]
    if (e2 == 0L || e2 %in% mm$nzrows) next
    r <- .ts_find_eq_solvables(ts, e2)
    if (!r$homogeneous) next
    k <- length(mm$nzrows) + 1L
    mm$nzrows[k] <- e2; mm$row_cols[[k]] <- r$cols; mm$row_vals[[k]] <- r$coeffs
  }
  list(mm = mm, modified = TRUE)
}

# AffectSystem (ModelingToolkitBase/src/systems/callbacks.jl) adds the
# parameters of the analytically integrated variables to the discrete
# parameters of every affect: an affect that changes x then solves the
# observed equation x ~ x#0 for x#0. Departure from Julia: only the constant
# of each variable the affect assigns is added (x#0 for x, xˍt#0 for xˍt),
# so that after D(D(x)) ~ 0 an event on x moves x#0 and keeps the slope
# xˍt#0; with all constants the affect would be underdetermined.
.st_affects_update_integrated <- function(sys) {
  integ <- unname(unclass(.f(sys, "analytically_integrated") %||% list()))
  if (!length(integ)) return(sys)
  ikeys <- vapply(integ, function(el) .key(default_toterm(el$var)), "")
  add <- function(a) {
    if (!inherits(a, "SymbolicAffect")) return(a)
    lk <- vapply(a$affect, function(e) .key(default_toterm(e$lhs)), "")
    hit <- lapply(integ[ikeys %in% lk], `[[`, "value")
    if (length(hit)) a$discrete_parameters <- .unique_syms(c(a$discrete_parameters, hit))
    a
  }
  fix <- function(cb) {
    for (fld in c("affect", "affect_neg", "initialize", "finalize")) if (!is.null(cb[[fld]])) cb[[fld]] <- add(cb[[fld]])
    cb
  }
  .update_system(sys, continuous_events = lapply(.f(sys, "continuous_events"), fix),
                 discrete_events = lapply(.f(sys, "discrete_events"), fix))
}

# eliminate_zero_variables_fixpoint!
.st_eliminate_zero_variables_fixpoint <- function(ts, mm, maxiters = 4L) {
  for (i in seq_len(maxiters)) {
    r <- .st_eliminate_zero_variables(ts, mm)
    mm <- r$mm
    if (!r$modified) break
  }
  mm
}

# ---------------------------------------------------------------------------
# The pipeline

# The integer Jacobian of equations with respect to variables, or NULL when
# an entry is not a small integer (the `jac` closure of dummy_derivative).
.st_integer_jacobian <- function(ts) {
  function(eqs, vars) {
    J <- matrix(0, length(eqs), length(vars))
    for (i in seq_along(eqs)) {
      eq <- ts$eqs[[eqs[i]]]
      nb <- ts$structure$graph$fadjlist[[eqs[i]]]
      cols <- which(vars %in% nb)
      if (!length(cols)) next
      ph <- lapply(nb, function(j) Symbol(paste0(".sv", j)))
      low <- .subst(eq$rhs - eq$lhs, .rule(ts$fullvars[nb], ph))
      for (j in cols) {
        a <- diff(low, ph[[match(vars[j], nb)]])
        if (!.is_num(a)) return(NULL)
        av <- as.numeric(a)
        if (av != round(av) || abs(av) > 127) return(NULL)
        J[i, j] <- av
      }
    }
    J
  }
}

# dummy_derivative on a tearing state: index reduction, state selection,
# tearing and reassembly.
# The solvable graph is computed here, from the equations as they stand after
# the alias passes, instead of being kept in step through those passes as in
# Julia; the passes above only need the incidence graph.
.st_dummy_derivative <- function(ts, fully_determined = TRUE, simplify = FALSE, tearing_alg = "carpanzano",
                                 maxiters = 8000L) {
  if (is.null(ts$structure$solvable_graph)) .ts_find_solvables(ts)
  m <- .pantelides(ts, maxiters = maxiters)
  tr <- .dummy_derivative_graph(ts, m, .st_integer_jacobian(ts), tearing_alg)
  .ts_reassemble(ts, tr, fully_determined = fully_determined, simplify = simplify)
}

# tearing on a tearing state: no index reduction; every differentiated
# variable is a state.
.st_tearing <- function(ts, fully_determined = TRUE, simplify = FALSE, tearing_alg = "carpanzano") {
  if (is.null(ts$structure$solvable_graph)) .ts_find_solvables(ts)
  tr <- .dummy_derivative_tearing(ts$structure, integer(), tearing_alg)
  .ts_reassemble(ts, tr, fully_determined = fully_determined, simplify = simplify)
}

.st_options <- function(ts, allow_symbolic = FALSE, allow_parameter = TRUE, conservative = FALSE) {
  ts$opts <- list(allow_symbolic = isTRUE(allow_symbolic), allow_parameter = isTRUE(allow_parameter),
                  conservative = isTRUE(conservative))
  invisible(ts)
}

#' Options of the structural pipeline of mtkcompile
#'
#' [mtkcompile()] passes its `...` to the structural simplification pipeline
#' (alias elimination, index reduction with the Pantelides algorithm, dummy
#' derivative state selection and tearing), which accepts these options.
#'
#' * `fully_determined`: `TRUE` requires as many equations as unknowns;
#'   `FALSE` skips the consistency check and only tears, keeping surplus
#'   variables as unknowns; `NULL` checks and decides.
#' * `check_consistency`: check that the system is balanced and not
#'   structurally singular (default `TRUE`).
#' * `dummy_derivative`: select states with the dummy derivative method
#'   (default `TRUE`). With `FALSE` the index reduced system keeps every
#'   differentiated variable as a state and drops the original constraints.
#' * `eliminate_mm_zeros`: eliminate variables that are identically zero and
#'   integrate their antiderivatives analytically, so that `D(x) ~ 0` gives
#'   the observed `x ~ x#0` with a parameter `x#0` solved at initialization.
#'   The default is `FALSE` (ModelingToolkit uses `TRUE`), which keeps such
#'   variables as states. With `TRUE`, an event that changes `x` updates
#'   `x#0`, as in ModelingToolkit.
#' * `tearing_alg`: `"carpanzano"` (default) or `"modia"`.
#' * `allow_symbolic`, `allow_parameter`, `conservative`: which coefficients
#'   an equation may be divided by when it is solved for a variable: ones
#'   containing unknowns (default `FALSE`) and ones containing parameters
#'   (default `TRUE`); with `conservative = TRUE` (default `FALSE`) a
#'   numeric coefficient must be 1 or -1, while symbolic coefficients still
#'   follow the first two options. A parameter whose default value is zero
#'   is never divided by.
#' * `maxiters`: iteration limit of the Pantelides algorithm.
#' * `inline_linear_sccs`: solve each block of algebraic equations that is
#'   linear in its variables inside the generated right-hand side (default
#'   `FALSE`), symbolically for blocks of at most
#'   `analytical_linear_scc_limit` variables (default 2) and numerically at
#'   run time for larger ones, so that the block's variables become observed.
#' * `inline_nonlinear_sccs`: an extension of modelcompR (default `FALSE`):
#'   solve each torn nonlinear algebraic block by a nested Newton iteration in
#'   the generated right-hand side, which turns such a DAE into an ODE.
#'
#' State selection honors the `state_priority` and `irreducible` metadata of
#' variables and the `state_priorities` and `irreducibles` arguments of
#' [System()]. Unbalanced systems raise conditions of class
#' `ExtraEquationsSystemException` or `ExtraVariablesSystemException`, and
#' structurally singular ones `InvalidSystemException`.
#'
#' @name mtkcompile_options
#' @seealso [tearing()], [TearingState()]
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)", "y(t)", "lam(t)"); parameters("g")
#' pend <- System(list(D(D(x)) ~ lam * x, D(D(y)) ~ lam * y - g, 0 ~ x^2 + y^2 - 1),
#'   t, name = "pendulum")
#' unknowns(mtkcompile(pend))
#' unknowns(mtkcompile(pend, dummy_derivative = FALSE))
NULL

# _mtkcompile! and _mtkcompile_worker!, registered as the hook
# `structural_compile` of mtkcompile(). `sys` arrives flattened with its
# connections expanded.
# Departure from Julia: `eliminate_mm_zeros` defaults to FALSE. Julia's
# default turns a state with D(x) ~ 0 into the observed x ~ x#0 with a new
# parameter x#0 that initialization solves for and that event affects update
# (.st_affects_update_integrated); here such states are kept as states.
# Pass eliminate_mm_zeros = TRUE to mtkcompile() for the Julia behavior.
.structural_compile <- function(sys, inputs = list(), outputs = list(), disturbance_inputs = list(),
                                fully_determined = TRUE, simplify = FALSE, check_consistency = TRUE,
                                dummy_derivative = TRUE, eliminate_mm_zeros = FALSE,
                                tearing_alg = c("carpanzano", "modia"), allow_symbolic = FALSE,
                                allow_parameter = TRUE, conservative = FALSE, maxiters = 8000L,
                                inline_linear_sccs = FALSE, analytical_linear_scc_limit = 2L,
                                inline_nonlinear_sccs = FALSE, ...) {
  tearing_alg <- match.arg(tearing_alg)
  ts <- TearingState(sys)
  .st_options(ts, allow_symbolic, allow_parameter, conservative)
  ts$opts$inline_linear_sccs <- isTRUE(inline_linear_sccs)
  ts$opts$analytical_linear_scc_limit <- as.integer(analytical_linear_scc_limit)
  ts$opts$inline_nonlinear_sccs <- isTRUE(inline_nonlinear_sccs)
  if (isTRUE(fully_determined) || identical(fully_determined, FALSE)) {
    check_consistency <- check_consistency && fully_determined
  } else {
    # NULL: find out whether the system is fully determined
    check_consistency <- TRUE
    fully_determined <- NULL
  }
  orig_inputs <- Filter(isinput, ts$fullvars)
  .st_validate_io(ts, "inputs", inputs)
  .st_validate_io(ts, "outputs", outputs)
  .st_validate_io(ts, "disturbance inputs", disturbance_inputs)
  .st_inputs_to_parameters(ts, .unique_syms(c(inputs, disturbance_inputs)), outputs)
  .st_eliminate_perfect_aliases(ts)
  .trivial_tearing(ts)
  mm <- .st_alias_elimination(ts)
  al <- .st_eliminate_perfect_aliases(ts)
  mm <- .ss_get_new_mm(mm, al$eq, al$var, al$alias_target, al$alias_sign)
  # after the second alias pass, so that a zero that is an alias eliminates
  # the right variable
  if (eliminate_mm_zeros) mm <- .st_eliminate_zero_variables_fixpoint(ts, mm)
  ts$mm <- mm
  if (check_consistency)
    fully_determined <- .check_consistency(ts, orig_inputs, nothrow = is.null(fully_determined))
  fully_determined <- isTRUE(fully_determined)

  sys <- if (fully_determined && dummy_derivative) {
    .st_dummy_derivative(ts, TRUE, simplify, tearing_alg, maxiters)
  } else if (fully_determined) {
    m <- .pantelides(ts, finalize = FALSE, maxiters = maxiters)
    ts2 <- TearingState(pantelides_reassemble(ts, m))
    ts2$opts <- ts$opts
    # keep what the passes above already eliminated (Julia starts from a
    # fresh state here and loses these observed equations)
    ts2$additional_observed <- ts$additional_observed
    ts2$pdm_from <- ts$pdm_from; ts2$pdm_to <- ts$pdm_to
    ts2$ad_from <- ts$ad_from; ts2$ad_to <- ts$ad_to
    ts2$mm <- .st_alias_elimination(ts2)
    .st_dummy_derivative(ts2, TRUE, simplify, tearing_alg, maxiters)
  } else {
    .st_tearing(ts, FALSE, simplify, tearing_alg)
  }
  sys <- .update_system(sys, observed = .toposort_observed(.f(sys, "observed")))
  .st_affects_update_integrated(sys)
}

# ---------------------------------------------------------------------------
# Exported structural transformations

#' Structural transformations of a system
#'
#' The passes that [mtkcompile()] combines, as the functions of the same
#' names in ModelingToolkit. Most users should call `mtkcompile()`.
#'
#' * `alias_elimination(sys)` removes the redundant equations among those
#'   that are integer linear combinations of variables and replaces the
#'   others by their reduced forms.
#' * `dae_index_lowering(sys)` runs the Pantelides algorithm and returns the
#'   system extended by the differentiated equations, which has index one.
#' * `dummy_derivative(sys)` reduces the index and selects states with the
#'   dummy derivative method, then tears and reassembles the system.
#' * `tearing(sys)` solves the equations that are linear in a variable for
#'   that variable, leaving the remaining (torn) equations as residuals. It
#'   does not reduce the index.
#' * `tearing_substitution(sys)` substitutes the observed equations into the
#'   equations of a simplified system.
#' * `pantelides_reassemble(state, var_eq_matching)` builds the system from a
#'   state after the Pantelides algorithm.
#' * `find_solvables(state)` (`find_solvables!` in Julia) fills
#'   `state$structure$solvable_graph` with the pairs of equations and
#'   variables such that the equation is linear in the variable with a
#'   coefficient that can be divided by.
#'
#' @param sys A `System`.
#' @param state A [TearingState()]; it is modified.
#' @param fully_determined Whether the system has as many equations as
#'   unknowns.
#' @param simplify Simplify the solved expressions.
#' @param tearing_alg `"carpanzano"` (the default) or `"modia"`.
#' @param allow_symbolic Allow solving for a variable whose coefficient
#'   contains other unknowns.
#' @param allow_parameter Allow solving for a variable whose coefficient
#'   contains parameters.
#' @param conservative Only divide by the numeric coefficients 1 and -1;
#'   symbolic coefficients follow `allow_symbolic` and `allow_parameter`.
#' @param var_eq_matching The matching returned by the Pantelides pass.
#' @param ... Options passed on: `maxiters` for the Pantelides algorithm.
#' @return A `System`; `find_solvables()` returns the state invisibly.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)", "y(t)", "z(t)"); parameters("k")
#' sys <- System(list(D(x) ~ -k * y, 0 ~ x - y + z, 0 ~ z - sin(x)), t, name = "sys")
#' torn <- tearing(sys)
#' equations(torn); observed(torn)
#' equations(tearing_substitution(torn))
#' ts <- TearingState(sys)
#' find_solvables(ts)
#' ts$structure$solvable_graph
#' # an index 2 system
#' sys2 <- System(list(D(x) ~ y, 0 ~ x - sin(t)), t, name = "idx")
#' equations(dae_index_lowering(sys2))
#' unknowns(dummy_derivative(sys2))
#' @export
tearing <- function(sys, state = TearingState(sys), fully_determined = TRUE, simplify = FALSE,
                    tearing_alg = c("carpanzano", "modia"), allow_symbolic = FALSE,
                    allow_parameter = TRUE, conservative = FALSE, ...) {
  .st_options(state, allow_symbolic, allow_parameter, conservative)
  .st_tearing(state, fully_determined, simplify, match.arg(tearing_alg))
}

#' @rdname tearing
#' @export
dummy_derivative <- function(sys, state = TearingState(sys), fully_determined = TRUE, simplify = FALSE,
                             tearing_alg = c("carpanzano", "modia"), allow_symbolic = FALSE,
                             allow_parameter = TRUE, conservative = FALSE, ...) {
  .st_options(state, allow_symbolic, allow_parameter, conservative)
  .st_dummy_derivative(state, fully_determined, simplify, match.arg(tearing_alg),
                       maxiters = list(...)$maxiters %||% 8000L)
}

#' @rdname tearing
#' @export
dae_index_lowering <- function(sys, ...) {
  state <- TearingState(sys)
  m <- .pantelides(state, finalize = FALSE, maxiters = list(...)$maxiters %||% 8000L)
  pantelides_reassemble(state, m)
}

#' @rdname tearing
#' @export
pantelides_reassemble <- function(state, var_eq_matching) {
  keep <- sort(unique(var_eq_matching[var_eq_matching > 0L]))
  .update_system(state$sys, eqs = state$original_eqs[keep],
                 unknowns = Filter(function(v) !.is_derivative(v), state$fullvars))
}

#' @rdname tearing
#' @export
alias_elimination <- function(sys) {
  state <- TearingState(sys)
  .st_alias_elimination(state)
  .update_system(state$sys, eqs = state$original_eqs)
}

#' @rdname tearing
#' @export
find_solvables <- function(state, allow_symbolic = FALSE, allow_parameter = TRUE, conservative = FALSE) {
  .st_options(state, allow_symbolic, allow_parameter, conservative)
  .ts_find_solvables(state)
  invisible(state)
}

#' @rdname tearing
#' @export
tearing_substitution <- function(sys) {
  .update_system(sys, eqs = full_equations(sys), schedule = NULL)
}

.bg_incidence_matrix <- function(g) {
  M <- matrix(0L, .bg_nsrcs(g), .bg_ndsts(g))
  for (e in seq_along(g$fadjlist)) M[e, g$fadjlist[[e]]] <- 1L
  M
}

#' Incidence matrices and variable masks
#'
#' `sorted_incidence_matrix(sys)` is the incidence matrix (equations by
#' unknowns, 1 where the unknown occurs in the equation) of a system
#' returned by [mtkcompile()], with the columns ordered by the strongly
#' connected components of its schedule. `but_ordered_incidence(ts)`
#' computes the block upper triangular ordering of a [TearingState()] and
#' returns the ordered incidence matrix and the block boundaries.
#' `highest_order_variable_mask(ts)` and `lowest_order_variable_mask(ts)`
#' return predicates on variable indices that select the variables without a
#' derivative, or that are not a derivative, in the state.
#'
#' @param sys A system returned by `mtkcompile()`.
#' @param ts A `TearingState`.
#' @param varmask A predicate on variable indices selecting the variables
#'   that take part in the ordering.
#' @return A 0/1 matrix; `but_ordered_incidence()` returns a list with the
#'   matrix and the block boundaries; the masks are functions.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)", "y(t)"); parameters("k")
#' eqs <- list(D(x) ~ -k * y, D(y) ~ x)
#' sorted_incidence_matrix(mtkcompile(System(eqs, t, name = "sys")))
#' ts <- TearingState(System(eqs, t, name = "sys"))
#' highest_order_variable_mask(ts)(seq_along(ts$fullvars))
#' but_ordered_incidence(ts)
#' @export
sorted_incidence_matrix <- function(sys) {
  ts <- if (inherits(sys, "System")) .f(sys, "tearing_state") else NULL
  sched <- if (inherits(sys, "System")) .f(sys, "schedule") else NULL
  if (is.null(ts) || is.null(sched) || !iscomplete(sys))
    stop("A simplified `System` is required. Call `mtkcompile` on the system first.", call. = FALSE)
  imat <- .bg_incidence_matrix(ts$structure$graph)
  q <- unlist(sched$var_sccs)
  imat[, c(q, setdiff(seq_len(ncol(imat)), q)), drop = FALSE]
}

#' @rdname sorted_incidence_matrix
#' @export
highest_order_variable_mask <- function(ts) {
  function(v) is.na(ts$structure$var_to_diff[v])
}

#' @rdname sorted_incidence_matrix
#' @export
lowest_order_variable_mask <- function(ts) {
  function(v) is.na(ts$structure$diff_to_var[v])
}

#' @rdname sorted_incidence_matrix
#' @export
but_ordered_incidence <- function(ts, varmask = highest_order_variable_mask(ts)) {
  g <- ts$structure$graph
  ok <- varmask(seq_len(.bg_ndsts(g)))
  m <- .maximal_matching(g, dstfilter = ok)
  vordering <- integer(); bb <- 1L; l <- 1L
  for (cmp in .find_var_sccs(g, m)) {
    vs <- cmp[ok[cmp]]
    if (!length(vs)) next
    vordering <- c(vordering, vs)
    l <- l + length(vs)
    bb <- c(bb, l)
  }
  vordering <- rev(vordering)
  rows <- m[vordering]
  list(.bg_incidence_matrix(g)[rows[rows > 0L], vordering, drop = FALSE], bb)
}

#' Which equation determines each variable
#'
#' For a system returned by [mtkcompile()], `map_variables_to_equations(sys)`
#' gives, for every unknown and observed variable, the equation that is used
#' to solve for it: the differential equation of a state, the matched
#' algebraic equation of an algebraic unknown, or the observed equation.
#'
#' @param sys A system returned by `mtkcompile()`.
#' @param rename_dummy_derivatives Name dummy derivative variables such as
#'   \ifelse{latex}{\verb{y_t(t)}}{\verb{yˍt(t)}} by their derivative form `Derivative(y(t), t)`.
#' @return A list of equations named by the printed form of the variables;
#'   the attribute `"vars"` holds the variables.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)", "y(t)")
#' sys <- mtkcompile(System(list(D(x) ~ -x, y ~ x + 1), t, name = "sys"))
#' map_variables_to_equations(sys)
#' @export
map_variables_to_equations <- function(sys, rename_dummy_derivatives = TRUE) {
  ts <- if (inherits(sys, "System")) .f(sys, "tearing_state") else NULL
  if (is.null(ts))
    stop("`map_variables_to_equations` requires a simplified system. Call `mtkcompile` on the ",
         "system before calling this function.", call. = FALSE)
  back <- list()
  sub <- .f(sys, "schedule")$dummy_sub
  if (rename_dummy_derivatives) for (k in names(sub)) {
    from <- .rule_lhs(sub, k)
    if (identical(.key(default_toterm(from)), .key(sub[[k]]))) back[[.key(sub[[k]])]] <- from
  }
  rename <- function(v) back[[.key(v)]] %||% v
  out <- list(); vars <- list()
  put <- function(v, eq) {
    v <- rename(v)
    out[[format(v)]] <<- eq
    vars[[format(v)]] <<- v
  }
  eqs <- equations(sys)
  for (eq in eqs) if (.is_derivative(eq$lhs)) put(.args(eq$lhs)[[1L]], eq)
  s <- ts$structure; g <- s$graph
  nd <- .bg_ndsts(g)
  algvars <- .ss_isalgvar(s, seq_len(nd))
  algeqs <- vapply(g$fadjlist, function(nb) !any(.ss_isdervar(s, nb)), TRUE)
  m <- .maximal_matching(g, algeqs, algvars)
  for (i in seq_len(nd)) if (m[i] > 0L) put(ts$fullvars[[i]], eqs[[m[i]]])
  for (eq in observed(sys)) put(eq$lhs, eq)
  structure(out, vars = vars)
}

# ---------------------------------------------------------------------------
# substitute_component (src/systems/substitute_component.jl)

# validate_replacement_rule: `rhs` must offer every unknown, parameter and
# subsystem of `lhs`, with the same connector and causality metadata.
.st_validate_replacement <- function(lhs, rhs, namespace = character()) {
  if (iscomplete(lhs)) stop("LHS of replacement rule cannot be completed.", call. = FALSE)
  if (iscomplete(rhs)) stop("RHS of replacement rule cannot be completed.", call. = FALSE)
  rhs_h <- .st_hierarchy(nameof(rhs))
  if (length(rhs_h) != 1L) stop("RHS of replacement rule must not be namespaced.", call. = FALSE)
  lhs_h <- .st_hierarchy(nameof(lhs))
  if (rhs_h != lhs_h[length(lhs_h)]) stop("LHS and RHS must have the same name.", call. = FALSE)
  if (!identical(.keys(independent_variables(lhs)), .keys(independent_variables(rhs))))
    stop("LHS and RHS of replacement rule must have the same independent variable.", call. = FALSE)
  where <- paste(c(namespace, nameof(lhs)), collapse = NAMESPACE_SEPARATOR)
  rhs_u <- .f(rhs, "unknowns")
  for (u in .f(lhs, "unknowns")) {
    i <- match(.key(u), .keys(rhs_u))
    if (is.na(i)) {
      if (!length(namespace)) .stopf("RHS of replacement rule does not contain unknown %s.", .fmt(u))
      .stopf("Subsystem %s of RHS does not contain unknown %s.", where, .fmt(u))
    }
    ru <- rhs_u[[i]]
    name <- paste(c(namespace, nameof(lhs), .name(u)), collapse = NAMESPACE_SEPARATOR)
    lc <- getconnect(u) %||% "Equality"
    if (!identical(lc, getconnect(ru) %||% "Equality"))
      .stopf("Variable %s should have connection metadata %s,", name, lc)
    if (isinput(u) != isinput(ru))
      .stopf("Variable %s has differing causality. Marked as `input = %s` in LHS and `input = %s` in RHS.",
             name, isinput(u), isinput(ru))
    if (isoutput(u) != isoutput(ru))
      .stopf("Variable %s has differing causality. Marked as `output = %s` in LHS and `output = %s` in RHS.",
             name, isoutput(u), isoutput(ru))
  }
  for (p in .f(lhs, "ps")) {
    if (.key(p) %in% .keys(.f(rhs, "ps"))) next
    if (!length(namespace)) .stopf("RHS of replacement rule does not contain parameter %s", .fmt(p))
    .stopf("Subsystem %s of RHS does not contain parameter %s.", where, .fmt(p))
  }
  rhs_s <- .f(rhs, "systems")
  for (s in .f(lhs, "systems")) {
    i <- match(nameof(s), vapply(rhs_s, nameof, ""))
    if (is.na(i)) .stopf("%s of replacement rule does not contain subsystem %s.", where, nameof(s))
    rs <- rhs_s[[i]]
    ct <- .f(s, "connector_type")
    if (!is.null(ct)) {
      name <- paste(c(namespace, nameof(lhs), nameof(s)), collapse = NAMESPACE_SEPARATOR)
      rct <- .f(rs, "connector_type")
      if (is.null(rct)) .stopf("Subsystem %s of RHS is not a connector.", name)
      if (!identical(ct, rct))
        .stopf("Subsystem %s of RHS has connection type %s but LHS has %s.", name, format(rct), format(ct))
    }
    .st_validate_replacement(s, rs, c(namespace, nameof(rhs)))
  }
  invisible(TRUE)
}

.st_hierarchy <- function(name) strsplit(name, paste0("[.", NAMESPACE_SEPARATOR, "]"))[[1L]]

# recursive_getproperty: the object at `hierarchy` below `sys`, namespaced
# from the first level on.
.st_recursive_getproperty <- function(sys, hierarchy) {
  cur <- toggle_namespacing(sys, FALSE)
  for (name in hierarchy) cur <- .getproperty(cur, name)
  cur
}

# recreate_connections: connections hold copies of the connectors they join;
# look them up again in the modified hierarchy.
.st_recreate_connections <- function(sys) {
  fresh <- function(arg) {
    if (is.null(arg)) return(NULL)
    name <- if (inherits(arg, "System")) nameof(arg) else .name(arg)
    .st_recursive_getproperty(sys, .st_hierarchy(name))
  }
  eqs <- lapply(.f(sys, "eqs"), function(eq) {
    if (inherits(eq, "AnalysisPoint")) {
      eq$input <- fresh(eq$input)
      eq$outputs <- lapply(eq$outputs, fresh)
    } else if (inherits(eq, "Connection")) {
      eq$systems <- lapply(eq$systems, fresh)
    }
    eq
  })
  .update_system(sys, eqs = eqs, systems = lapply(.f(sys, "systems"), .st_recreate_connections))
}

# Replace the subsystem of `sys` at `path` (names below the root).
.st_replace_subsystem <- function(sys, path, new) {
  if (!length(path)) return(new)
  subs <- .f(sys, "systems")
  i <- match(path[1L], vapply(subs, nameof, ""))
  if (is.na(i)) .stopf("System %s does not have a subsystem named %s.", nameof(sys), path[1L])
  subs[[i]] <- .st_replace_subsystem(subs[[i]], path[-1L], new)
  .update_system(sys, systems = subs)
}

#' Replace a component of a hierarchical model
#'
#' `substitute_component(sys, old, new)` returns `sys` with the subsystem
#' `old` (obtained as `sys$name`, or `sys$outer$name` for a nested one)
#' replaced by `new`, as `substitute_component(sys, old => new)` in
#' ModelingToolkit. `new` must have the same name as the component it
#' replaces and offer all of its unknowns, parameters and subsystems, with
#' the same connector types and input or output markings, so that every
#' connection to the component stays valid.
#'
#' @param sys A hierarchical `System` that is not complete.
#' @param old The subsystem to replace, namespaced by `sys`; or a list
#'   `list(old, new)`.
#' @param new The replacement system.
#' @return The modified `System`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); k <- Parameter("k")
#' decay <- System(list(D(x) ~ -k * x), t, name = "part")
#' growth <- System(list(D(x) ~ k * x), t, name = "part")
#' model <- System(list(), t, systems = list(decay), name = "model")
#' equations(substitute_component(model, model$part, growth))
#' @export
substitute_component <- function(sys, old, new = NULL) {
  if (!inherits(old, "System")) {
    new <- old[[2L]]
    old <- old[[1L]]
  }
  if (iscomplete(sys)) stop("Cannot replace subsystems of completed systems", call. = FALSE)
  .st_validate_replacement(old, new)
  hierarchy <- .st_hierarchy(nameof(old))
  if (nameof(sys) != hierarchy[1L])
    .stopf("The component %s is not namespaced by the system %s; obtain it as `sys$name`.",
           nameof(old), nameof(sys))
  .st_recreate_connections(.st_replace_subsystem(sys, hierarchy[-1L], new))
}

# ---------------------------------------------------------------------------
# isolate_subsystem (src/systems/analysis_points.jl)

#' Isolate the part of a model between analysis points
#'
#' `isolate_subsystem(sys, input_aps, output_aps)` returns the part of the
#' unsimplified hierarchical model `sys` bounded by the input and output
#' analysis points: the components reachable from the inside faces of those
#' points without crossing them, with the connections among them. All other
#' components and all equations involving them are removed. The analysis
#' points are given as `sys$name` or by name, singly or in a list.
#'
#' @param sys A hierarchical `System` that has not been completed or
#'   compiled.
#' @param input_aps,output_aps Analysis points bounding the subsystem.
#' @return A list with the isolated system `sys`, the variables `input_vars`
#'   at the inside face of each input point and `output_vars` at the inside
#'   face of each output point.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' block <- function(name, eq) {
#'   u <- Variable("u", t, input = TRUE); y <- Variable("y", t, output = TRUE)
#'   System(list(eq(u, y)), t, name = name)
#' }
#' P <- block("P", function(u, y) D(y) ~ u - y)
#' C <- block("C", function(u, y) y ~ -u)
#' cl <- System(list(connect(C$y, "u", P$u), connect(P$y, "y", C$u)), t,
#'   systems = list(P, C), name = "cl")
#' iso <- isolate_subsystem(cl, "u", "y")
#' vapply(get_systems(iso$sys), nameof, "")
#' iso$input_vars; iso$output_vars
#' @export
isolate_subsystem <- function(sys, input_aps, output_aps) {
  canon <- function(aps) {
    if (inherits(aps, "AnalysisPoint") || is.character(aps)) aps <- as.list(aps)
    if (inherits(aps, "AnalysisPoint")) aps <- list(aps)
    vapply(aps, function(a) {
      if (is.character(a)) return(paste(c(nameof(sys), .st_hierarchy(a)), collapse = NAMESPACE_SEPARATOR))
      nm <- .st_hierarchy(a$name)
      if (nm[1L] != nameof(sys)) nm <- c(nameof(sys), nm)
      paste(nm, collapse = NAMESPACE_SEPARATOR)
    }, "")
  }
  if (inherits(input_aps, "AnalysisPoint")) input_aps <- list(input_aps)
  if (inherits(output_aps, "AnalysisPoint")) output_aps <- list(output_aps)
  in_names <- canon(input_aps); out_names <- canon(output_aps)
  boundary <- union(in_names, out_names)
  pkey <- function(p) paste(p, collapse = NAMESPACE_SEPARATOR)

  # 1. index every subsystem of the hierarchy by its path of names
  paths <- list()
  collect <- function(cur, path) for (s in .f(cur, "systems")) {
    p <- c(path, nameof(s))
    paths[[length(paths) + 1L]] <<- p
    collect(s, p)
  }
  collect(sys, character())
  keys <- vapply(paths, pkey, "")
  n <- length(paths)
  adj <- matrix(FALSE, n, n)
  idx_of <- function(p) if (is.null(p)) NA_integer_ else match(pkey(p), keys)
  full_ap_name <- function(parent, local) pkey(c(nameof(sys), parent, local))
  # the component a connector or causal variable belongs to
  conn_path <- function(conn, parent) {
    segs <- .st_hierarchy(if (inherits(conn, "System")) nameof(conn) else .name(conn))
    if (length(segs) <= 1L) NULL else c(parent, segs[-length(segs)])
  }
  add_edge <- function(p1, p2) {
    i <- idx_of(p1); j <- idx_of(p2)
    if (!is.na(i) && !is.na(j) && i != j) adj[i, j] <<- adj[j, i] <<- TRUE
  }

  # 2. connections join components; boundary analysis points are cut
  seeds <- integer(); input_vars <- list(); output_vars <- list(); cuts <- list()
  build <- function(cur, parent) {
    for (eq in .f(cur, "eqs")) {
      if (inherits(eq, "AnalysisPoint")) {
        fname <- full_ap_name(parent, eq$name)
        is_boundary <- fname %in% boundary
        outs <- eq$outputs %||% list()
        if (is_boundary && fname %in% in_names) for (cn in outs) {
          input_vars[[length(input_vars) + 1L]] <<- ap_var(cn)
          seeds <<- c(seeds, idx_of(conn_path(cn, parent)))
        }
        if (is_boundary && fname %in% out_names && !is.null(eq$input)) {
          output_vars[[length(output_vars) + 1L]] <<- ap_var(eq$input)
          seeds <<- c(seeds, idx_of(conn_path(eq$input, parent)))
        }
        in_path <- if (is.null(eq$input)) NULL else conn_path(eq$input, parent)
        for (cn in outs) {
          out_path <- conn_path(cn, parent)
          if (is.null(out_path) || is.null(in_path)) next
          add_edge(in_path, out_path)
          i <- idx_of(in_path); j <- idx_of(out_path)
          if (is_boundary && !is.na(i) && !is.na(j) && i != j) cuts[[length(cuts) + 1L]] <<- c(i, j)
        }
      } else if (inherits(eq, "Connection")) {
        cps <- Filter(Negate(is.null), lapply(eq$systems, conn_path, parent = parent))
        for (i in seq_along(cps)) for (j in seq_along(cps)) if (i < j) add_edge(cps[[i]], cps[[j]])
      }
    }
    for (s in .f(cur, "systems")) build(s, c(parent, nameof(s)))
  }
  build(sys, character())
  for (ct in cuts) adj[ct[1L], ct[2L]] <- adj[ct[2L], ct[1L]] <- FALSE

  # 3. everything reachable from the inside faces is inside
  inside <- logical(n)
  queue <- unique(seeds[!is.na(seeds)])
  inside[queue] <- TRUE
  while (length(queue)) {
    v <- queue[1L]; queue <- queue[-1L]
    new <- which(adj[v, ] & !inside)
    inside[new] <- TRUE
    queue <- c(queue, new)
  }
  prefixes <- unique(unlist(lapply(paths[inside], function(p)
    vapply(seq_along(p), function(k) pkey(p[seq_len(k)]), ""))))
  has_inside <- function(p) pkey(p) %in% prefixes

  # 4. rebuild the hierarchy: inside components are kept whole, containers
  # keep only their inside subsystems and the connections among them
  rebuild <- function(cur, parent) {
    i <- idx_of(if (length(parent)) parent else NULL)
    if (!is.na(i) && inside[i]) return(cur)
    subs <- Filter(function(s) has_inside(c(parent, nameof(s))), .f(cur, "systems"))
    subs <- lapply(subs, function(s) rebuild(s, c(parent, nameof(s))))
    eqs <- Filter(function(eq) {
      if (inherits(eq, "AnalysisPoint")) {
        if (full_ap_name(parent, eq$name) %in% boundary) return(FALSE)
        conns <- c(if (!is.null(eq$input)) list(eq$input), eq$outputs)
      } else if (inherits(eq, "Connection")) {
        conns <- eq$systems
      } else return(FALSE)
      ps <- Filter(Negate(is.null), lapply(conns, conn_path, parent = parent))
      all(vapply(ps, has_inside, TRUE))
    }, .f(cur, "eqs"))
    System(eqs, .f(cur, "iv"), list(), list(), name = nameof(cur), systems = subs)
  }
  list(sys = rebuild(sys, character()), input_vars = input_vars, output_vars = output_vars)
}

# Register the pipeline with mtkcompile() (see tools/ARCHITECTURE.md).
.mtk$hooks$structural_compile <- .structural_compile
