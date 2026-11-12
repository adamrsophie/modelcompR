# Ported from StateSelection.jl (https://github.com/JuliaComputing/StateSelection.jl),
# commit 29b3c90. Portions are licensed under the GNU AGPL, version 3; see
# inst/COPYRIGHTS and inst/STATESELECTION-AGPL-LICENSE.
#
# The tearing state of a system and its reassembly into a simplified system,
# after StateSelection.jl/lib/ModelingToolkitTearing/src: tearingstate.jl
# (TearingState, canonical_sort_key, build_canonical_ranks,
# build_state_priorities, sort_fullvars, build_var_to_diff),
# stateselection_interface.jl (var_derivative!, eq_derivative!,
# linear_subsys_adjmat!, find_eq_solvables!, _check_allow_symbolic_parameter),
# trivial_tearing_interface.jl (possibly_explicit_equations) and reassemble.jl
# (substitute_derivatives_algevars!, generate_derivative_variables!,
# generate_system_equations!, blt_reorder_generated_equations!,
# reorder_vars!, update_simplified_system!, DefaultReassembleAlgorithm).
#
# A TearingState is an environment with
#   sys            the flattened system (metadata; the equations live in eqs)
#   iv             the independent variable
#   fullvars       list of variables: unknowns and their derivative nodes
#   structure      see structural-stateselection.R
#   eqs            equations in the form 0 ~ rhs, in graph order
#   original_eqs   the same equations as written by the user
#   additional_observed  observed equations found by the alias passes
#   always_present logical, variables that count as used without incidence
#   mm             integer linear subsystem matrix, or NULL
#   pdm_from, pdm_to     derivative rules D(v) => expr (param_derivative_map)
#   opts           solvability options: allow_symbolic, allow_parameter,
#                  conservative
#   ad_from, ad_to       analytical_derivatives of eliminated variables
# Every equation is stored as 0 ~ rhs - lhs (Julia keeps differential
# equations as written); the form the user wrote stays in original_eqs. A
# symbol declared as a parameter that is missing from the parameter list of
# the system is treated as a parameter, where Julia raises an error.
# Array equations arrive scalarized (structural-arrays.R); clock inference
# runs before the tearing state is built (structural-clocks.R), and inline
# solves of blocks are in structural-inline.R. The equation sort key
# (get_equation_sort_key!) orders the terms of a sum with a port of
# SymbolicUtils' degree ordering that is exact for polynomial terms and
# approximate for other expressions.

# D(e) with nested derivatives merged.
.st_D <- function(e, iv) {
  if (.is_derivative(e)) {
    di <- .deriv_info(e)
    return(.make_derivative(di$var, di$iv, di$order + 1L))
  }
  .make_derivative(e, iv, 1L)
}

# Variable leaves and derivative nodes of an expression. Operator
# applications such as Initial(x) and Pre(x) are atomic and skipped, and a
# derivative node hides its inner variable (OperatorIsAtomic in Julia).
.ts_search_vars <- function(e) {
  out <- list()
  seen <- new.env(parent = emptyenv())
  walk <- function(z) {
    k <- .key(z)
    if (!is.null(seen[[k]])) return(invisible())
    seen[[k]] <- TRUE
    ty <- .type(z)
    if (ty == "Derivative") {
      out[[length(out) + 1L]] <<- z
      return(invisible())
    }
    if (ty %in% c("Symbol", "Function")) {
      kind <- .var_kind(z)
      if (identical(kind, "operator")) return(invisible())
      if (!is.null(kind)) {
        out[[length(out) + 1L]] <<- z
        return(invisible())
      }
    }
    for (a in .args(z)) walk(a)
  }
  walk(e)
  out
}

.ts_has_subexpr <- function(e, key) {
  if (identical(.key(e), key)) return(TRUE)
  for (a in .args(e)) if (.ts_has_subexpr(a, key)) return(TRUE)
  FALSE
}

# canonical_sort_key: (operator signature, name, array index).
.ts_canonical_order <- function(vars) {
  if (!length(vars)) return(integer())
  base <- lapply(vars, function(v) if (.is_derivative(v)) .args(v)[[1L]] else v)
  opsig <- vapply(vars, function(v) as.integer(.is_derivative(v)), 0L)
  nm <- vapply(base, function(b) .meta_get(b, "array_parent") %||% .name(b) %||% format(b), "")
  idx <- t(vapply(base, function(b) {
    ai <- as.integer(.meta_get(b, "array_index") %||% 0L)
    c(ai, 0L, 0L)[1:2]
  }, integer(2)))
  # byte order of the names, as Julia compares Symbols
  order(opsig, xtfrm(factor(nm, levels = sort(unique(nm), method = "radix"))), idx[, 1L], idx[, 2L])
}

# get_degrees and monomial_lt of SymbolicUtils (src/ordering.jl, MIT): the
# degree lexicographic order of terms used by `<ₑ`.
.ts_degrees <- function(e) {
  none <- list(names = character(), degs = numeric())
  if (.is_num(e)) return(none)
  ty <- .type(e)
  if (ty == "Derivative") return(list(names = "Differential", degs = 1))
  if (ty %in% c("Symbol", "Function")) return(list(names = .name(e), degs = 1))
  a <- .args(e)
  if (ty == "Pow") {
    if (!.is_num(a[[2L]])) return(list(names = "^", degs = 1))
    d <- .ts_degrees(a[[1L]])
    d$degs <- d$degs * as.numeric(a[[2L]])
    return(d)
  }
  if (ty == "Add") {
    best <- none; bsum <- 0
    for (z in a) {
      d <- .ts_degrees(z)
      ds <- sum(d$degs)
      if (!length(best$names) || ds > bsum || (ds == bsum && .ts_lexlt(d, best))) {
        best <- d; bsum <- ds
      }
    }
    return(best)
  }
  if (ty == "Mul") {
    ds <- lapply(a, .ts_degrees)
    names <- unlist(lapply(ds, `[[`, "names")); degs <- unlist(lapply(ds, `[[`, "degs"))
    if (!length(names)) return(none)
    # a quotient of two symbolic expressions is a single opaque term in Julia
    if (any(degs < 0) && any(degs > 0)) return(list(names = "zzz__ZZZ_unknown", degs = 1))
    o <- order(xtfrm(factor(names, levels = sort(unique(names), method = "radix"))), degs)
    return(list(names = names[o], degs = degs[o]))
  }
  list(names = "zzz__ZZZ_unknown", degs = 1)
}

.ts_lexlt <- function(d1, d2) {
  for (i in seq_len(min(length(d1$names), length(d2$names)))) {
    if (d1$names[i] == d2$names[i]) {
      a <- d1$degs[i]; b <- d2$degs[i]
      if (a != b) return(abs(a) > abs(b) || (abs(a) == abs(b) && a > b))
    } else {
      return(identical(sort(c(d1$names[i], d2$names[i]), method = "radix")[1L], d1$names[i]))
    }
  }
  FALSE
}

.ts_monomial_lt <- function(d1, d2) {
  s1 <- sum(d1$degs); s2 <- sum(d2$degs)
  if (s1 != s2) return(abs(s1) < abs(s2) || (abs(s1) == abs(s2) && s1 > s2))
  .ts_lexlt(d1, d2)
}

# Terms in the order of `<ₑ`; ties keep the canonical order of their keys.
.ts_sort_terms <- function(terms) {
  n <- length(terms)
  if (n <= 1L) return(seq_len(n))
  degs <- lapply(terms, .ts_degrees)
  idx <- order(.keys(terms), method = "radix")
  out <- integer()
  for (i in idx) {
    pos <- length(out)
    while (pos > 0L && .ts_monomial_lt(degs[[i]], degs[[out[pos]]])) pos <- pos - 1L
    out <- append(out, i, pos)
  }
  out
}

# get_expression_sort_key!: rows of (canonical rank, coefficient, exponent)
# for the variables of an expression. `ranks` maps keys of fullvars to
# canonical ranks.
.ts_expr_sort_key <- function(e, ranks, cache) {
  k <- .key(e)
  hit <- cache[[k]]
  if (!is.null(hit)) return(hit)
  empty <- matrix(numeric(), 0L, 3L)
  r <- ranks[[k]]
  ty <- .type(e)
  a <- .args(e)
  bind <- function(parts) {
    out <- empty
    for (p in parts) {
      out <- rbind(out, p())
      if (nrow(out) > 100L) break
    }
    out
  }
  out <- if (.is_num(e)) matrix(c(0, as.numeric(e), 1), 1L)
  else if (!is.null(r)) matrix(c(r, 1, 1), 1L)
  else if (ty == "Add") {
    a <- Filter(function(z) !.is_num(z), a)
    coef <- vapply(a, function(z) {
      if (.type(z) == "Mul" && .is_num(.args(z)[[1L]])) as.numeric(.args(z)[[1L]]) else 1
    }, 0)
    terms <- lapply(a, function(z) {
      if (.type(z) == "Mul" && .is_num(.args(z)[[1L]])) Reduce(`*`, .args(z)[-1L]) else z
    })
    bind(lapply(.ts_sort_terms(terms), function(i) function() {
      kk <- .ts_expr_sort_key(terms[[i]], ranks, cache)
      kk[, 2L] <- kk[, 2L] * coef[i]
      kk
    }))
  } else if (ty == "Mul") {
    num <- vapply(a, .is_num, TRUE)
    cf <- prod(vapply(a[num], as.numeric, 0))
    fac <- a[!num]
    expo <- vapply(fac, function(z) if (.type(z) == "Pow" && .is_num(.args(z)[[2L]])) as.numeric(.args(z)[[2L]]) else 1, 0)
    base <- lapply(fac, function(z) if (.type(z) == "Pow" && .is_num(.args(z)[[2L]])) .args(z)[[1L]] else z)
    bind(lapply(.ts_sort_terms(base), function(i) function() {
      kk <- .ts_expr_sort_key(base[[i]], ranks, cache)
      cbind(kk[, 1L], abs(kk[, 2L])^expo[i] * cf, kk[, 3L] + expo[i], deparse.level = 0)
    }))
  } else if (ty == "Pow" && .is_num(a[[2L]])) {
    kk <- .ts_expr_sort_key(a[[1L]], ranks, cache)
    v <- as.numeric(a[[2L]])
    cbind(kk[, 1L], abs(kk[, 2L])^v, kk[, 3L] + v, deparse.level = 0)
  } else if (ty == "Symbol" || !length(a)) empty
  else bind(lapply(a, function(z) function() .ts_expr_sort_key(z, ranks, cache)))
  dim(out) <- c(length(out) %/% 3L, 3L)
  assign(k, out, envir = cache)
  out
}

# sortperm of the equation sort keys: lexicographic, a prefix sorts first.
.ts_equation_order <- function(exprs, ranks) {
  if (length(exprs) <= 1L) return(seq_along(exprs))
  cache <- new.env(parent = emptyenv())
  keys <- lapply(exprs, function(e) as.vector(t(.ts_expr_sort_key(e, ranks, cache))))
  w <- max(lengths(keys), 1L)
  cols <- lapply(seq_len(w), function(j) vapply(keys, function(k) if (j <= length(k)) k[j] else -Inf, 0))
  do.call(order, cols)
}

.ts_is_irreducible <- function(ts, v) {
  if (.is_derivative(v)) return(FALSE)
  .key(v) %in% ts$irreducible_keys || isirreducible(v)
}

#' The structural state of a system
#'
#' `TearingState(sys)` builds the structure that the simplification passes
#' of [mtkcompile()] work on, as the type of the same name in
#' ModelingToolkit: the list `fullvars` of unknowns and their derivatives,
#' the bipartite incidence graph between equations and `fullvars`
#' (`structure$graph`), and the derivative relations between variables
#' (`structure$var_to_diff`, `NA` for none). It is a mutable environment.
#'
#' @param sys A `System` with an independent variable or without one.
#' @param sort_eqs Order the equations by their canonical sort key, so that
#'   the result does not depend on the order in which they were written.
#' @param ... Unused.
#' @return An object of class `TearingState`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)", "y(t)"); parameters("k")
#' sys <- System(list(D(x) ~ -k * y, 0 ~ x - y), t, name = "sys")
#' ts <- TearingState(sys)
#' ts$fullvars
#' ts$structure$graph
#' ts$structure$var_to_diff
#' @export
TearingState <- function(sys, sort_eqs = TRUE, ...) {
  if (!inherits(sys, "System")) stop("TearingState() needs a System", call. = FALSE)
  sys <- flatten(sys)
  iv <- .f(sys, "iv")
  eqs0 <- lapply(Filter(function(e) inherits(e, "Equation") && is_sym(e$lhs) && is_sym(e$rhs),
                        .f(sys, "eqs")), expand_derivatives)
  dv_keys <- .keys(.f(sys, "unknowns"))
  ps <- .f(sys, "ps")
  ps_keys <- .keys(ps)
  obs_keys <- .keys(lapply(.f(sys, "observed"), `[[`, "lhs"))
  ivk <- if (is.null(iv)) "" else .key(iv)

  # derivatives of time dependent parameters are zero between callbacks
  pdm_from <- list(); pdm_to <- list()
  fullvars <- list()
  var2idx <- new.env(parent = emptyenv())
  addvar <- function(v) {
    k <- .key(v)
    if (is.null(var2idx[[k]])) {
      fullvars[[length(fullvars) + 1L]] <<- v
      var2idx[[k]] <- length(fullvars)
    }
    invisible()
  }
  eqs <- vector("list", length(eqs0))
  incidence <- vector("list", length(eqs0))
  for (i in seq_along(eqs0)) {
    eq <- eqs0[[i]]
    term <- if (.is_zero(eq$lhs)) eq$rhs else eq$rhs - eq$lhs
    inc <- list()
    zero_from <- list()
    for (v in .ts_search_vars(term)) {
      k <- .key(v)
      if (k == ivk || k %in% ps_keys || k %in% obs_keys) next
      if (.is_derivative(v)) {
        base <- .args(v)[[1L]]
        bk <- .key(base)
        if (bk %in% ps_keys || identical(.var_kind(base), "parameter") || identical(.var_kind(base), "discrete")) {
          zero_from[[length(zero_from) + 1L]] <- v
          next
        }
        if (!(bk %in% dv_keys))
          .stopf("Variable %s occurs differentiated as %s but is not an unknown of the system.",
                 .fmt(base), .fmt(v))
        if (.deriv_info(v)$n_specs > 1L || !identical(.key(.deriv_info(v)$iv), ivk))
          .stopf("%s is not a derivative with respect to the independent variable.", .fmt(v))
      } else if (!(k %in% dv_keys)) {
        if (!is.null(iv) && .type(v) == "Function" && .is_delay_or_evalat(v, iv)) next
        if (!identical(.var_kind(v), "unknown")) next
        .stopf(paste0("`%s` occurs in the system but is neither an unknown nor a parameter of it.\n",
                      "Declare it in the system, in whichever list it belongs to:\n",
                      "    System(eqs, t, list(unknowns...), list(parameters..., %s))"), .fmt(v), .fmt(v))
      }
      addvar(v)
      inc[[length(inc) + 1L]] <- v
    }
    if (length(zero_from)) {
      zeros <- rep(list(.S0()), length(zero_from))
      term <- .subst(term, .rule(zero_from, zeros))
      eqs0[[i]] <- .subst(eq, .rule(zero_from, zeros))
      pdm_from <- c(pdm_from, zero_from); pdm_to <- c(pdm_to, zeros)
    }
    eqs[[i]] <- Equation(.S0(), term)
    incidence[[i]] <- inc
  }
  # intermediate derivative orders: D(D(x)) brings D(x) and x
  is_der <- logical()
  i <- 1L
  while (i <= length(fullvars)) {
    v <- fullvars[[i]]
    is_der[i] <- .is_derivative(v)
    if (is_der[i]) {
      di <- .deriv_info(v)
      addvar(di$var)
      for (j in seq_len(di$order - 1L)) addvar(.make_derivative(di$var, di$iv, j))
    }
    i <- i + 1L
  }
  # sort_fullvars: derivative variables first, then the variables they
  # differentiate, then the rest, so the mass matrix is as diagonal as possible
  if (!is.null(iv) && length(fullvars)) {
    lower <- function(v) {
      di <- .deriv_info(v)
      .make_derivative(di$var, di$iv, di$order - 1L)
    }
    ders <- fullvars[is_der]
    ders <- ders[.ts_canonical_order(ders)]
    sorted <- .unique_syms(c(ders, lapply(ders, lower), fullvars[.ts_canonical_order(fullvars)]))
    fullvars <- sorted
  } else ders <- list()
  fk <- .keys(fullvars)
  nvars <- length(fullvars)
  cranks <- order(.ts_canonical_order(fullvars)) * 100L
  if (sort_eqs && length(eqs) > 1L) {
    ranks <- list2env(stats::setNames(as.list(cranks), fk), parent = emptyenv())
    # a differential equation is keyed by its right-hand side as written
    exprs <- lapply(seq_along(eqs), function(i) {
      if (any(vapply(incidence[[i]], .is_derivative, TRUE))) eqs0[[i]]$rhs else eqs[[i]]$rhs
    })
    o <- .ts_equation_order(exprs, ranks)
    eqs <- eqs[o]; eqs0 <- eqs0[o]; incidence <- incidence[o]
  }
  graph <- BipartiteGraph(fadjlist = lapply(incidence, function(inc) match(.keys(inc), fk)), ndsts = nvars)
  # state priorities: the highest priority along each derivative chain
  spm <- .f(sys, "state_priorities")
  sp <- vapply(fullvars, function(v) {
    p <- if (.is_derivative(v)) NULL else smap_get(spm, v)
    as.integer(round(as.numeric(p %||% 0)))
  }, 0L)
  s <- .ss_structure(graph, nvars, sp, cranks)
  for (d in seq_along(ders)) {
    lk <- .key(lower(fullvars[[d]]))
    .dg_set(s, "var", match(lk, fk), d)
  }
  for (i in which(is.na(s$diff_to_var))) {
    p <- s$state_priorities[i]; var <- i
    repeat {
      p <- max(p, s$state_priorities[var])
      s$state_priorities[var] <- p
      var <- s$var_to_diff[var]
      if (is.na(var)) break
    }
  }
  # unknowns that occur in no equation lose their bindings and initial values
  unused <- setdiff(dv_keys, .keys(lapply(fullvars, function(v) if (.is_derivative(v)) .args(v)[[1L]] else v)))
  binds <- .f(sys, "bindings"); ics <- .f(sys, "initial_conditions")
  for (k in unused) {
    binds[[k]] <- NULL
    ics[[k]] <- NULL
  }
  ts <- new.env(parent = emptyenv())
  ts$sys <- .update_system(sys, eqs = eqs, bindings = binds, initial_conditions = ics)
  ts$iv <- iv
  ts$fullvars <- fullvars
  ts$structure <- s
  ts$eqs <- eqs
  ts$original_eqs <- eqs0
  ts$additional_observed <- list()
  ts$always_present <- logical(nvars)
  ts$mm <- NULL
  ts$pdm_from <- pdm_from; ts$pdm_to <- pdm_to
  ts$ad_from <- list(); ts$ad_to <- list()
  ts$irreducible_keys <- .keys(.f(sys, "irreducibles"))
  ts$maybe_zero_keys <- .ts_discover_maybe_zeros(sys)
  class(ts) <- "TearingState"
  ts
}

#' @export
print.TearingState <- function(x, ...) {
  g <- x$structure$graph
  cat("TearingState of ", nameof(x$sys), " with ", .bg_nsrcs(g), " equations and ", .bg_ndsts(g),
      " variables\n", sep = "")
  for (i in seq_along(x$fullvars)) {
    d <- x$structure$var_to_diff[i]
    cat(sprintf("  %3d  %s%s\n", i, format(x$fullvars[[i]]),
                if (is.na(d)) "" else paste0("   (derivative: ", d, ")")))
  }
  invisible(x)
}

# discover_maybe_zeros: parameters and unknowns whose default value is zero
# must not be divided by when solving equations.
.ts_discover_maybe_zeros <- function(sys) {
  mz <- .keys(.f(sys, "maybe_zeros"))
  defs <- smap_merge(.f(sys, "bindings"), .f(sys, "initial_conditions"))
  defs <- defs[!vapply(unclass(defs), function(el) is.null(el$var), TRUE)]
  if (!length(defs)) return(mz)
  class(defs) <- "mtk_symmap"
  vals <- tryCatch(.evaluate_varmap(defs), error = function(e) list())
  zero <- names(vals)[vapply(vals, function(z) is.numeric(z) && length(z) == 1L && !is.na(z) && z == 0, TRUE)]
  union(mz, zero)
}

# query_maybe_zeros: whether an expression can vanish because one of its
# factors is marked as possibly zero.
.ts_maybe_zero <- function(e, mz) {
  if (!length(mz)) return(FALSE)
  if (.key(e) %in% mz) return(TRUE)
  if (.is_num(e)) return(.is_zero(e))
  ty <- .type(e)
  a <- .args(e)
  if (ty == "Mul") return(any(vapply(a, .ts_maybe_zero, TRUE, mz = mz)))
  if (ty == "Add") return(all(vapply(a, .ts_maybe_zero, TRUE, mz = mz)))
  if (ty == "Pow") return(.ts_maybe_zero(a[[1L]], mz))
  if (ty == "Function" && .name(e) %in% c("sin", "Max", "Min"))
    return(any(vapply(a, .ts_maybe_zero, TRUE, mz = mz)))
  FALSE
}

# ---------------------------------------------------------------------------
# Linear expansion and solvability

# expr = a * var + b, holding every derivative node independent.
.st_linear_expansion <- function(expr, var) {
  ders <- Filter(.is_derivative, .ts_search_vars(expr))
  hold <- .unique_syms(c(list(var), ders))
  ph <- lapply(seq_along(hold), function(i) Symbol(paste0(".sv", i)))
  low <- .subst(expr, .rule(hold, ph))
  back <- .rule(ph, hold)
  p <- ph[[1L]]
  a <- diff(low, p)
  list(a = .subst(a, back), b = .subst(.subst(low, .rule(list(p), list(.S0()))), back),
       islinear = !.has_symbol(a, p))
}

# Whether symbol `s` occurs in expression `e`, by structure rather than by
# the spelling of symcompR keys.
.has_symbol <- function(e, s) .key(s) %in% .keys(free_symbols(e))

# _check_allow_symbolic_parameter on a coefficient in the lowered form, where
# every variable of the system is a placeholder symbol named ".sv<index>".
.ts_allow_division <- function(ts, a, allow_symbolic, allow_parameter) {
  if (allow_symbolic) return(TRUE)
  if (!allow_parameter) return(FALSE)
  if (.ts_maybe_zero(a, ts$maybe_zero_keys)) return(FALSE)
  # the coefficient must not contain any variable placeholder
  !any(startsWith(vapply(free_symbols(a), function(z) .name(z) %||% "", ""), ".sv"))
}

# find_eq_solvables!: for each variable of equation `ieq`, whether the
# equation is linear in it with a coefficient that may be divided by. Edges
# to variables with a zero coefficient are removed from the graph. Returns
# whether the equation is an integer linear combination of its variables
# (`all_int_vars` with a zero remainder, `homogeneous`), with the columns and
# coefficients of that combination.
.ts_find_eq_solvables <- function(ts, ieq, may_be_zero = TRUE) {
  s <- ts$structure; g <- s$graph; sg <- s$solvable_graph
  allow_symbolic <- isTRUE(ts$opts$allow_symbolic)
  allow_parameter <- !identical(ts$opts$allow_parameter, FALSE)
  conservative <- isTRUE(ts$opts$conservative)
  eq <- ts$eqs[[ieq]]
  term <- eq$rhs - eq$lhs
  nb <- g$fadjlist[[ieq]]
  ph <- lapply(nb, function(j) Symbol(paste0(".sv", j)))
  low <- if (length(nb)) .subst(term, .rule(ts$fullvars[nb], ph)) else term
  all_int <- TRUE
  cols <- integer(); coeffs <- numeric(); to_rm <- integer()
  for (i in seq_along(nb)) {
    j <- nb[i]
    if (.ts_is_irreducible(ts, ts$fullvars[[j]])) {
      all_int <- FALSE
      next
    }
    a <- diff(low, ph[[i]])
    if (.has_symbol(a, ph[[i]])) {
      all_int <- FALSE
      next
    }
    if (!.is_num(a)) {
      all_int <- FALSE
      if (.ts_allow_division(ts, a, allow_symbolic, allow_parameter) && !is.null(sg)) .bg_add_edge(sg, ieq, j)
      next
    }
    av <- as.numeric(a)
    if (av != round(av) || abs(av) > 127) {
      all_int <- FALSE
      if (conservative) next
    } else if (conservative && abs(av) > 1) {
      # the row is incomplete without this coefficient, so it cannot join the
      # integer linear subsystem (Julia leaves the row inconsistent here)
      all_int <- FALSE
      next
    } else if (av != 0 || !may_be_zero) {
      cols <- c(cols, j); coeffs <- c(coeffs, av)
    }
    if (av != 0) {
      if (!is.null(sg)) .bg_add_edge(sg, ieq, j)
    } else if (may_be_zero) to_rm <- c(to_rm, j)
  }
  for (j in to_rm) .bg_rem_edge(g, ieq, j)
  homogeneous <- all_int &&
    .is_zero(if (length(nb)) .subst(low, .rule(ph, rep(list(.S0()), length(ph)))) else low)
  list(all_int_vars = all_int, homogeneous = homogeneous, cols = cols, coeffs = coeffs, to_rm = to_rm)
}

# find_solvables!: populate the solvable graph.
.ts_find_solvables <- function(ts) {
  s <- ts$structure
  s$solvable_graph <- BipartiteGraph(.bg_nsrcs(s$graph), .bg_ndsts(s$graph))
  for (ieq in seq_len(.bg_nsrcs(s$graph))) .ts_find_eq_solvables(ts, ieq)
  invisible(ts)
}

# linear_subsys_adjmat!: the equations that are integer linear combinations
# of variables, as a sparse matrix of coefficients.
.ts_linear_subsys_adjmat <- function(ts) {
  g <- ts$structure$graph
  nz <- integer(); cols <- list(); vals <- list()
  for (i in seq_len(.bg_nsrcs(g))) {
    r <- .ts_find_eq_solvables(ts, i)
    if (!r$homogeneous) next
    k <- length(nz) + 1L
    nz[k] <- i; cols[[k]] <- r$cols; vals[[k]] <- r$coeffs
  }
  list(nparentrows = .bg_nsrcs(g), ncols = .bg_ndsts(g), nzrows = nz, row_cols = cols, row_vals = vals)
}

# possibly_explicit_equations: (equation, variable) pairs where the user
# wrote `variable ~ expression` with the variable absent from the expression.
.ts_possibly_explicit_equations <- function(ts) {
  fk <- .keys(ts$fullvars)
  out <- list()
  for (i in seq_along(ts$original_eqs)) {
    eq <- ts$original_eqs[[i]]
    j <- match(.key(eq$lhs), fk)
    if (is.na(j) || .ts_is_irreducible(ts, eq$lhs)) next
    if (.is_zero(ts$eqs[[i]]$lhs) && .is_zero(ts$eqs[[i]]$rhs)) next
    if (.ts_has_subexpr(eq$rhs, .key(eq$lhs))) next
    out[[length(out) + 1L]] <- c(i, j)
  }
  out
}

# ---------------------------------------------------------------------------
# Differentiation of variables and equations (Pantelides)

# sum(coeffs * fullvars[vars])
.ts_linear_combination <- function(ts, coeffs, vars) {
  if (!length(vars)) return(.S0())
  Reduce(`+`, Map(function(cf, v) cf * ts$fullvars[[v]], coeffs, vars))
}

# var_derivative!: add the derivative of variable v to the state.
.ts_var_derivative <- function(ts, v) {
  s <- ts$structure
  vd <- .bg_add_dst(s$graph)
  if (!is.null(s$solvable_graph)) .bg_add_dst(s$solvable_graph)
  s$var_to_diff[vd] <- NA_integer_
  s$diff_to_var[vd] <- NA_integer_
  .dg_set(s, "var", v, vd)
  ts$fullvars[[vd]] <- .st_D(ts$fullvars[[v]], ts$iv)
  s$state_priorities[vd] <- s$state_priorities[v]
  s$canonical_ranks[vd] <- s$canonical_ranks[v] + 1L
  ts$always_present[vd] <- ts$always_present[v]
  if (!is.null(ts$mm)) ts$mm$ncols <- ts$mm$ncols + 1L
  vd
}

# eq_derivative!: add the time derivative of equation ieq to the state.
.ts_eq_derivative <- function(ts, ieq) {
  s <- ts$structure; g <- s$graph
  ed <- .bg_add_src(g)
  if (!is.null(s$solvable_graph)) .bg_add_src(s$solvable_graph)
  s$eq_to_diff[ed] <- NA_integer_
  s$diff_to_eq[ed] <- NA_integer_
  .dg_set(s, "eq", ieq, ed)
  mm <- ts$mm
  idx <- if (is.null(mm)) NA_integer_ else match(ieq, mm$nzrows)
  if (!is.null(mm)) ts$mm$nparentrows <- mm$nparentrows + 1L
  if (!is.na(idx)) {
    # the derivative of an integer linear row is the same row in the derivatives
    rcol <- s$var_to_diff[mm$row_cols[[idx]]]
    rval <- mm$row_vals[[idx]]
    o <- order(rcol)
    rcol <- rcol[o]; rval <- rval[o]
    k <- length(mm$nzrows) + 1L
    ts$mm$nzrows[k] <- ed
    ts$mm$row_cols[[k]] <- rcol
    ts$mm$row_vals[[k]] <- rval
    ts$eqs[[ed]] <- Equation(.S0(), .ts_linear_combination(ts, rval, rcol))
    ts$original_eqs[[ed]] <- ts$eqs[[ed]]
    for (v in rcol) {
      .bg_add_edge(g, ed, v)
      if (!is.null(s$solvable_graph)) .bg_add_edge(s$solvable_graph, ed, v)
    }
    return(ed)
  }
  eq <- ts$eqs[[ieq]]
  d <- diff(eq$rhs - eq$lhs, ts$iv)
  if (length(ts$pdm_from)) d <- .subst(d, .rule(ts$pdm_from, ts$pdm_to))
  # derivatives of parameters and of eliminated observed variables that have
  # no rule are zero: parameters stay constant between callbacks
  stray <- Filter(function(z) .is_derivative(z) && !(.key(z) %in% .keys(ts$fullvars)) &&
                    isparameter(.args(z)[[1L]]), .ts_search_vars(d))
  if (length(stray)) d <- .subst(d, .rule(stray, rep(list(.S0()), length(stray))))
  ts$eqs[[ed]] <- Equation(.S0(), d)
  ts$original_eqs[[ed]] <- ts$eqs[[ed]]
  for (v in g$fadjlist[[ieq]]) {
    .bg_add_edge(g, ed, v)
    .bg_add_edge(g, ed, s$var_to_diff[v])
  }
  if (!is.null(s$solvable_graph)) {
    r <- .ts_find_eq_solvables(ts, ed)
    if (!is.null(mm) && r$homogeneous) {
      k <- length(ts$mm$nzrows) + 1L
      ts$mm$nzrows[k] <- ed
      ts$mm$row_cols[[k]] <- r$cols
      ts$mm$row_vals[[k]] <- r$coeffs
    }
  } else {
    # keep the incidence exact: drop variables whose coefficient vanished
    present <- .keys(.ts_search_vars(d))
    nb <- g$fadjlist[[ed]]
    .bg_set_neighbors(g, ed, nb[.keys(ts$fullvars[nb]) %in% present])
  }
  ed
}

# ---------------------------------------------------------------------------
# Reassembly

.match_set <- function(m, i, v) {
  if (v > 0L) m[m == v] <- .UNASSIGNED
  m[i] <- as.integer(v)
  m
}

.ts_var_order <- function(s, dv) {
  order <- 0L
  while (!is.na(lower <- s$diff_to_var[dv])) {
    order <- order + 1L
    dv <- lower
  }
  list(order = order, lv = dv)
}

# Fixed point substitution of derivative rules D(x) => expr (total_sub, an
# ExpandDerivativeDict): a higher derivative without a rule of its own is
# obtained by differentiating the rule of the first derivative.
.ts_fixpoint_sub <- function(e, sub, iv, maxiter = 100L) {
  if (!length(ls(sub))) return(e)
  for (it in seq_len(maxiter)) {
    from <- list(); to <- list()
    for (d in Filter(.is_derivative, .ts_search_vars(e))) {
      hit <- sub[[.key(d)]]
      if (is.null(hit)) {
        di <- .deriv_info(d)
        if (di$order > 1L) {
          first <- sub[[.key(.make_derivative(di$var, di$iv, 1L))]]
          if (!is.null(first)) {
            hit <- expand_derivatives(.make_derivative(first, di$iv, di$order - 1L))
            assign(.key(d), hit, envir = sub)
          }
        }
      }
      if (!is.null(hit)) {
        from[[length(from) + 1L]] <- d
        to[[length(to) + 1L]] <- hit
      }
    }
    if (!length(from)) break
    e <- .subst(e, .rule(from, to))
  }
  e
}

# find_duplicate_dd: an existing equation 0 ~ D(x) - v, so that v can serve
# as the dummy derivative variable.
.ts_find_duplicate_dd <- function(ts, dv, mm) {
  if (is.null(mm)) return(NULL)
  s <- ts$structure
  for (eq in s$solvable_graph$badjlist[[dv]]) {
    mi <- match(eq, mm$nzrows)
    if (is.na(mi)) next
    cs <- mm$row_cols[[mi]]; vs <- mm$row_vals[[mi]]
    if (length(cs) == 2L && abs(vs[1L]) == 1 && vs[1L] == -vs[2L] && dv %in% cs) {
      v_t <- setdiff(cs, dv)
      if (is.na(s$diff_to_var[v_t])) return(c(eq, v_t))
    }
  }
  NULL
}

# contract_variables: the graph on the kept equations and variables, where a
# solved variable is replaced by the variables it is computed from.
.contract_variables <- function(g, m, var_rename, eq_rename, n_eq, n_var) {
  deps <- new.env(parent = emptyenv())
  var_deps <- function(v) {
    hit <- deps[[as.character(v)]]
    if (!is.null(hit)) return(hit)
    seen <- logical(length(var_rename)); stack <- v
    while (length(stack)) {
      x <- stack[length(stack)]; stack <- stack[-length(stack)]
      if (seen[x]) next
      seen[x] <- TRUE
      e <- m[x]
      if (e > 0L) stack <- c(stack, setdiff(g$fadjlist[[e]], x))
    }
    out <- var_rename[which(seen)]
    out <- out[out != 0L]
    assign(as.character(v), out, envir = deps)
    out
  }
  newg <- BipartiteGraph(n_eq, n_var)
  for (e in seq_along(g$fadjlist)) {
    ne <- eq_rename[e]
    if (ne == 0L) next
    row <- integer()
    for (v in g$fadjlist[[e]]) row <- c(row, if (var_rename[v] != 0L) var_rename[v] else var_deps(v))
    .bg_set_neighbors(newg, ne, row)
  }
  newg
}

# DefaultReassembleAlgorithm: turn the torn structure into the simplified
# system: differential equations D(x) ~ f, algebraic residuals 0 ~ g and
# observed equations for the solved variables, in block lower triangular
# order.
.ts_reassemble <- function(ts, tr, fully_determined = TRUE, simplify = FALSE, ...) {
  s <- ts$structure; g <- s$graph; sg <- s$solvable_graph
  iv <- ts$iv
  m <- tr$var_eq_matching; fullm <- tr$full_var_eq_matching; var_sccs <- tr$var_sccs
  nd0 <- .bg_ndsts(g)
  m <- m[seq_len(nd0)]; fullm <- fullm[seq_len(nd0)]
  m[is.na(m)] <- .UNASSIGNED; fullm[is.na(fullm)] <- .UNASSIGNED

  # get_extra_eqs_vars
  extra_eqs <- integer(); extra_vars <- integer()
  if (!fully_determined) {
    extra_vars <- which(fullm <= 0L & m == .UNASSIGNED)
    extra_eqs <- setdiff(seq_len(.bg_nsrcs(g)), fullm[fullm > 0L])
  }
  extra_unknowns <- ts$fullvars[extra_vars]
  neweqs <- ts$eqs
  dummy_from <- list(); dummy_to <- list()

  if (!is.null(iv)) {
    # substitute_derivatives_algevars!: derivatives of variables that are not
    # selected states become algebraic dummy derivative variables
    for (var in seq_along(ts$fullvars)) {
      dv <- s$var_to_diff[var]
      if (is.na(dv) || m[var] == .SELECTED_STATE) next
      dd <- ts$fullvars[[dv]]
      v_t <- default_toterm(dd)
      nbrs <- g$badjlist[[dv]]
      if (length(nbrs)) {
        dummy_from[[length(dummy_from) + 1L]] <- dd
        dummy_to[[length(dummy_to) + 1L]] <- v_t
        rule <- .rule(list(dd), list(v_t))
        for (eq in nbrs) neweqs[[eq]] <- .subst(neweqs[[eq]], rule)
      }
      ts$fullvars[[dv]] <- v_t
      dx <- dv; x_t <- v_t
      while (!is.na(ddx <- s$var_to_diff[dx])) {
        dx_t <- .st_D(x_t, iv)
        rule <- .rule(list(ts$fullvars[[ddx]]), list(dx_t))
        for (eq in g$badjlist[[ddx]]) neweqs[[eq]] <- .subst(neweqs[[eq]], rule)
        ts$fullvars[[ddx]] <- dx_t
        dx <- ddx; x_t <- dx_t
      }
      .dg_cut(s, "var", dv)
    }

    # generate_derivative_variables!: a state derivative D(x) that no
    # equation is solved for gets a variable x_t and the equation
    # 0 ~ D(x) - x_t, which lowers the system to first order
    scc_i <- integer(nd0); scc_j <- integer(nd0)
    for (i in seq_along(var_sccs)) {
      scc_i[var_sccs[[i]]] <- i
      scc_j[var_sccs[[i]]] <- seq_along(var_sccs[[i]])
    }
    pairs <- list()
    for (v in seq_len(length(s$var_to_diff))) {
      dv <- s$var_to_diff[v]
      if (is.na(dv) || m[dv] > 0L) next
      dd <- .ts_find_duplicate_dd(ts, dv, ts$mm)
      if (is.null(dd)) {
        dx <- ts$fullvars[[dv]]
        vo <- .ts_var_order(s, dv)
        x_t <- .diff2term(ts$fullvars[[vo$lv]], iv, vo$order)
        v_t <- .bg_add_dst(g)
        .bg_add_dst(sg)
        ts$fullvars[[v_t]] <- x_t
        s$state_priorities[v_t] <- s$state_priorities[dv]
        s$canonical_ranks[v_t] <- s$canonical_ranks[dv]
        ts$always_present[v_t] <- FALSE
        s$var_to_diff[v_t] <- NA_integer_
        s$diff_to_var[v_t] <- NA_integer_
        dummy_eq <- .bg_add_src(g)
        .bg_add_src(sg)
        s$eq_to_diff[dummy_eq] <- NA_integer_
        s$diff_to_eq[dummy_eq] <- NA_integer_
        neweqs[[dummy_eq]] <- Equation(.S0(), dx - x_t)
        .bg_add_edge(g, dummy_eq, dv)
        .bg_add_edge(g, dummy_eq, v_t)
        .bg_add_edge(sg, dummy_eq, dv)
        for (e in g$badjlist[[dv]]) .bg_add_edge(g, e, v_t)
        m[v_t] <- .UNASSIGNED; fullm[v_t] <- .UNASSIGNED
        dd <- c(dummy_eq, v_t)
      }
      dummy_eq <- dd[1L]; v_t <- dd[2L]
      .dg_set(s, "var", v_t, s$var_to_diff[dv])
      old <- fullm[dv]
      fullm <- .match_set(fullm, dv, dummy_eq)
      m <- .match_set(m, dv, dummy_eq)
      fullm <- .match_set(fullm, v_t, old)
      pairs[[length(pairs) + 1L]] <- c(v_t, dv)
    }
    if (length(pairs)) {
      chain_height <- function(dv) {
        h <- 0L; v <- dv
        repeat {
          if (v > length(s$var_to_diff)) break
          v <- s$var_to_diff[v]
          if (is.na(v)) break
          h <- h + 1L
        }
        h
      }
      ins_pos <- integer(); ins_var <- integer()
      remove <- vector("list", length(var_sccs))
      for (pr in pairs) {
        v_t <- pr[1L]; dv <- pr[2L]
        i <- scc_i[dv]
        var_sccs[[i]][scc_j[dv]] <- v_t
        if (v_t <= nd0) remove[[scc_i[v_t]]] <- c(remove[[scc_i[v_t]]], scc_j[v_t])
        ddv <- s$var_to_diff[dv]
        ins_pos <- c(ins_pos, if (!is.na(ddv) && ddv <= nd0) min(i, scc_i[ddv]) else i)
        ins_var <- c(ins_var, dv)
      }
      for (i in seq_along(var_sccs)) if (length(remove[[i]])) var_sccs[[i]] <- var_sccs[[i]][-remove[[i]]]
      o <- order(ins_pos, -vapply(ins_var, chain_height, 0L))
      ins_pos <- ins_pos[o]; ins_var <- ins_var[o]
      new_sccs <- list()
      for (i in seq_along(var_sccs)) {
        for (dv in ins_var[ins_pos == i]) new_sccs[[length(new_sccs) + 1L]] <- dv
        new_sccs[[length(new_sccs) + 1L]] <- var_sccs[[i]]
      }
      var_sccs <- Filter(length, new_sccs)
    }
  }

  # generate_system_equations!
  total_sub <- new.env(parent = emptyenv())
  out_eqs <- list(); eq_ordering <- integer(); var_ordering <- integer(); eq_scc <- integer()
  solved_eqs <- list(); solved_vars <- integer()
  NO_SCC <- .Machine$integer.max
  inv <- .inv_matching(m, .bg_nsrcs(g))
  done <- logical(.bg_nsrcs(g))
  codegen <- function(ieq, v, scc_idx) {
    done[ieq] <<- TRUE
    eq <- neweqs[[ieq]]
    issolvable <- v > 0L && .bg_has_edge(sg, ieq, v)
    if (issolvable && !is.na(s$diff_to_var[v])) {
      if (is.null(iv)) stop("Differential found in a non-differential system: ", .fmt(eq), call. = FALSE)
      le <- .st_linear_expansion(eq$rhs - eq$lhs, ts$fullvars[[v]])
      dx <- .st_D(ts$fullvars[[.ts_var_order(s, v)$lv]], iv)
      rhs <- .ts_fixpoint_sub(-le$b / le$a, total_sub, iv)
      for (e in g$badjlist[[v]]) {
        if (e == ieq) next
        for (w in g$fadjlist[[ieq]]) .bg_add_edge(g, e, w)
        .bg_rem_edge(g, e, v)
      }
      assign(.key(dx), rhs, envir = total_sub)
      k <- length(out_eqs) + 1L
      out_eqs[[k]] <<- Equation(dx, rhs)
      eq_ordering[k] <<- ieq; var_ordering[k] <<- s$diff_to_var[v]; eq_scc[k] <<- scc_idx
    } else if (issolvable) {
      var <- ts$fullvars[[v]]
      le <- .st_linear_expansion(eq$lhs - eq$rhs, var)
      if (.is_zero(le$a)) {
        warning("Tearing: solving ", .fmt(eq), " for ", .fmt(var), " is singular!", call. = FALSE)
        return(invisible())
      }
      rhs <- -le$b / le$a
      if (simplify) rhs <- symcompR::simplify(rhs)
      solved_eqs[[length(solved_eqs) + 1L]] <<- Equation(var, .ts_fixpoint_sub(rhs, total_sub, iv))
      solved_vars <<- c(solved_vars, v)
    } else {
      k <- length(out_eqs) + 1L
      out_eqs[[k]] <<- Equation(.S0(), .ts_fixpoint_sub(eq$rhs - eq$lhs, total_sub, iv))
      eq_ordering[k] <<- ieq; var_ordering[k] <<- 0L; eq_scc[k] <<- scc_idx
    }
    invisible()
  }
  for (eq in extra_eqs) if (inv[eq] > 0L) codegen(eq, inv[eq], NO_SCC)
  for (i in seq_along(var_sccs)) {
    # get_sorted_scc: solved equations in dependency order, then the residuals
    scc <- var_sccs[[i]]
    scc_eqs <- fullm[scc]; scc_eqs <- scc_eqs[scc_eqs > 0L]
    solved <- m[scc]; solved <- solved[solved > 0L]
    if (length(solved) > 1L) {
      succ <- function(k) {
        w <- m[g$fadjlist[[solved[k]]]]
        w <- match(w[w > 0L & w != solved[k]], solved)
        w[!is.na(w)]
      }
      solved <- solved[rev(.toposort_dfs(length(solved), succ))]
    }
    escc <- c(solved, setdiff(scc_eqs, solved))
    vscc <- inv[escc]; vscc <- vscc[vscc > 0L]
    vscc <- c(vscc, setdiff(scc, vscc))
    var_sccs[[i]] <- vscc
    if (length(escc) != length(vscc)) {
      escc <- setdiff(escc, extra_eqs)
      if (!length(escc) || !length(setdiff(vscc, extra_vars))) next
    }
    # inline linear (and nonlinear) solves of the block, structural-inline.R
    inl <- .ts_inline_scc(ts, escc, vscc, neweqs, m, fullm, inv, total_sub)
    if (!is.null(inl)) {
      fk <- .keys(ts$fullvars)
      new_inc <- unique(unlist(lapply(inl$sol, function(z) match(.keys(.ts_search_vars(z)), fk))))
      new_inc <- new_inc[!is.na(new_inc)]
      for (k in seq_along(inl$vars)) {
        ieq <- inl$eqs[k]; v <- inl$vars[k]; rhs <- inl$sol[[k]]
        done[ieq] <- TRUE
        if (!is.na(s$diff_to_var[v])) {
          dx <- .st_D(ts$fullvars[[.ts_var_order(s, v)$lv]], iv)
          for (e in g$badjlist[[v]]) {
            if (e == ieq) next
            for (w in new_inc) .bg_add_edge(g, e, w)
            .bg_rem_edge(g, e, v)
          }
          assign(.key(dx), rhs, envir = total_sub)
          kk <- length(out_eqs) + 1L
          out_eqs[[kk]] <- Equation(dx, rhs)
          eq_ordering[kk] <- ieq; var_ordering[kk] <- s$diff_to_var[v]; eq_scc[kk] <- i
        } else {
          solved_eqs[[length(solved_eqs) + 1L]] <- Equation(ts$fullvars[[v]], rhs)
          solved_vars <- c(solved_vars, v)
          # v is now computed by ieq from the inputs of the solve, so the
          # contraction of the graph below reaches them through v
          m[v] <- ieq
          .bg_set_neighbors(g, ieq, unique(c(v, new_inc)))
          .bg_set_neighbors(sg, ieq, unique(c(v, new_inc)))
        }
      }
      for (ieq in setdiff(escc, inl$eqs)) codegen(ieq, inv[ieq], i)
      next
    }
    for (ieq in escc) codegen(ieq, inv[ieq], i)
  }
  for (eq in extra_eqs) if (inv[eq] <= 0L) codegen(eq, .UNASSIGNED, NO_SCC)
  missed <- which(!done & lengths(g$fadjlist) > 0L)
  if (length(missed))
    stop("internal error: equation(s) not assigned to any block during reassembly:\n",
         .ss_print_array(neweqs[missed]), call. = FALSE)
  # an equation without unknowns determines nothing: it is dropped when it
  # holds identically, is an error when it cannot hold, and is reported when
  # it only constrains parameters
  for (ieq in which(!done)) {
    res <- neweqs[[ieq]]$rhs - neweqs[[ieq]]$lhs
    if (.is_zero(res)) next
    if (.is_num(res))
      .ss_error("InvalidSystemException", paste0(
        "InvalidSystemException: the equation ", .fmt(neweqs[[ieq]]),
        " contains no unknowns and cannot hold."))
    .warnf("The equation %s contains no unknowns of the system; it is not part of the simplified system.",
           .fmt(neweqs[[ieq]]))
  }

  diff_vars <- var_ordering[var_ordering != 0L]
  if (anyDuplicated(diff_vars)) stop("Tearing internal error: lowering DAE into semi-implicit ODE failed!", call. = FALSE)
  ispresent <- function(i) length(g$badjlist[[i]]) > 0L ||
    (!is.na(s$var_to_diff[i]) && length(g$badjlist[[s$var_to_diff[i]]]) > 0L)
  is_alg_unknown <- function(j) !(j %in% diff_vars || j %in% solved_vars || j %in% extra_vars) &&
    is.na(s$diff_to_var[j]) && ispresent(j)

  # blt_reorder_generated_equations!
  nd <- .bg_ndsts(g)
  scc_pos <- integer(nd)
  for (i in seq_along(var_sccs)) scc_pos[var_sccs[[i]]] <- i
  used <- integer(length(var_sccs))
  for (k in seq_along(out_eqs)) {
    if (var_ordering[k] != 0L || eq_scc[k] == NO_SCC) next
    i <- eq_scc[k]
    pool <- Filter(is_alg_unknown, var_sccs[[i]])
    used[i] <- used[i] + 1L
    if (used[i] > length(pool))
      stop("BLT reordering: a block has more algebraic equations than unknowns: ", .fmt(out_eqs[[k]]), call. = FALSE)
    var_ordering[k] <- pool[used[i]]
  }
  key <- ifelse(eq_scc == NO_SCC, NO_SCC, scc_pos[pmax(var_ordering, 1L)])
  perm <- order(key)
  out_eqs <- out_eqs[perm]; eq_ordering <- eq_ordering[perm]; var_ordering <- var_ordering[perm]
  var_ordering <- var_ordering[var_ordering != 0L]
  var_ordering <- c(var_ordering, setdiff(seq_len(nd), c(var_ordering, solved_vars)))

  # reorder_vars!
  eqsperm <- integer(.bg_nsrcs(g)); eqsperm[eq_ordering] <- seq_along(eq_ordering)
  varsperm <- integer(nd); varsperm[var_ordering] <- seq_along(var_ordering)
  new_s <- .ss_structure(
    .contract_variables(g, m, varsperm, eqsperm, length(eq_ordering), length(var_ordering)),
    length(var_ordering), s$state_priorities[var_ordering], s$canonical_ranks[var_ordering])
  new_s$solvable_graph <- .contract_variables(sg, m, varsperm, eqsperm, length(eq_ordering), length(var_ordering))
  for (v in seq_along(s$var_to_diff)) {
    d <- s$var_to_diff[v]
    if (varsperm[v] == 0L || is.na(d)) next
    .dg_set(new_s, "var", varsperm[v], if (varsperm[d] > 0L) varsperm[d] else NA_integer_)
  }
  for (e in seq_along(s$eq_to_diff)) {
    d <- s$eq_to_diff[e]
    if (eqsperm[e] == 0L || is.na(d)) next
    .dg_set(new_s, "eq", eqsperm[e], if (eqsperm[d] > 0L) eqsperm[d] else NA_integer_)
  }
  var_sccs <- Filter(length, lapply(var_sccs, function(scc) {
    x <- varsperm[scc]
    x[x != 0L]
  }))
  ts$structure <- new_s
  ts$fullvars <- ts$fullvars[var_ordering]
  ts$always_present <- ts$always_present[var_ordering]
  ts$eqs <- out_eqs
  ts$original_eqs <- out_eqs

  # update_simplified_system!
  g2 <- new_s$graph
  present2 <- function(i) length(g2$badjlist[[i]]) > 0L ||
    (!is.na(new_s$var_to_diff[i]) && length(g2$badjlist[[new_s$var_to_diff[i]]]) > 0L)
  for (eq in out_eqs) if (.is_derivative(eq$lhs)) {
    dummy_from[[length(dummy_from) + 1L]] <- eq$lhs
    dummy_to[[length(dummy_to) + 1L]] <- eq$rhs
  }
  obs_rule <- .rule(dummy_from, dummy_to)
  sys <- ts$sys
  obs <- c(.subst(.f(sys, "observed"), obs_rule), solved_eqs, .subst(ts$additional_observed, obs_rule))
  unknown_idxs <- Filter(function(i) is.na(new_s$diff_to_var[i]) && present2(i), seq_along(ts$fullvars))
  unknowns <- .unique_syms(c(ts$fullvars[unknown_idxs], extra_unknowns))
  var_sccs <- Filter(length, lapply(var_sccs, function(scc) scc[scc %in% unknown_idxs]))
  sys <- .update_system(sys, eqs = out_eqs, unknowns = unknowns, observed = obs,
                        schedule = list(var_sccs = var_sccs,
                                        dummy_sub = .rule(c(dummy_from, ts$ad_from), c(dummy_to, ts$ad_to))),
                        isscheduled = TRUE)
  ts$sys <- sys
  .update_system(sys, tearing_state = ts)
}
