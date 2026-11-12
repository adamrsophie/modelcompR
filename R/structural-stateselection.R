# Ported from StateSelection.jl (https://github.com/JuliaComputing/StateSelection.jl),
# commit 29b3c90. Portions are licensed under the GNU AGPL, version 3; see
# inst/COPYRIGHTS and inst/STATESELECTION-AGPL-LICENSE.
#
# Structural algorithms on the incidence graph of a system, after
# StateSelection.jl/src: graph/diff.jl (DiffGraph), interface.jl
# (SystemStructure, rm_eqs_vars!, isdervar and friends), utils.jl
# (check_consistency, singular_check, error_reporting, get_new_mm),
# math/bareiss.jl and math/sparsematrixclil.jl (fraction free Bareiss),
# singularity_removal.jl (aag_bareiss!, structural_singularity_removal!),
# pantelides.jl (computed_highest_diff_variables, pantelides!),
# partial_state_selection.jl (dummy_derivative_graph!, DummyDerivativeTearing),
# tearing.jl (trivial_tearing!), carpanzano_tearing.jl and modia_tearing.jl.
#
# Data layout. A structure is an environment with
#   graph, solvable_graph   BipartiteGraph (equations to variables)
#   var_to_diff, diff_to_var  integer vectors, NA for "no derivative"
#   eq_to_diff, diff_to_eq    the same for equations
#   state_priorities, canonical_ranks  integer vectors over variables
# The linear subsystem matrix (SparseMatrixCLIL) is a list with nzrows
# (equation index of each row), row_cols (sorted variable indices) and
# row_vals (integer coefficients stored as doubles).

.ss_structure <- function(graph, nvars, state_priorities = integer(nvars),
                          canonical_ranks = seq_len(nvars) * 100L) {
  s <- new.env(parent = emptyenv())
  s$graph <- graph
  s$solvable_graph <- NULL
  s$var_to_diff <- rep(NA_integer_, nvars)
  s$diff_to_var <- rep(NA_integer_, nvars)
  ne <- .bg_nsrcs(graph)
  s$eq_to_diff <- rep(NA_integer_, ne)
  s$diff_to_eq <- rep(NA_integer_, ne)
  s$state_priorities <- as.integer(state_priorities)
  s$canonical_ranks <- as.integer(canonical_ranks)
  s
}

# DiffGraph setindex!: var_to_diff[v] = d keeps the inverse in step and, as
# in Julia, leaves the forward edge of a previous preimage of d in place.
.dg_set <- function(s, kind, v, d) {
  fwd <- paste0(kind, "_to_diff"); bwd <- paste0("diff_to_", kind)
  old <- s[[fwd]][v]
  if (!is.na(old)) s[[bwd]][old] <- NA_integer_
  if (!is.na(d)) s[[bwd]][d] <- as.integer(v)
  s[[fwd]][v] <- as.integer(d)
  invisible(s)
}

# invview(var_to_diff)[d] = nothing
.dg_cut <- function(s, kind, d) {
  fwd <- paste0(kind, "_to_diff"); bwd <- paste0("diff_to_", kind)
  old <- s[[bwd]][d]
  if (!is.na(old)) s[[fwd]][old] <- NA_integer_
  s[[bwd]][d] <- NA_integer_
  invisible(s)
}

.ss_isdervar <- function(s, i) !is.na(s$diff_to_var[i])
.ss_isalgvar <- function(s, i) is.na(s$var_to_diff[i]) & is.na(s$diff_to_var[i])
.ss_isdiffvar <- function(s, i) !is.na(s$var_to_diff[i]) & is.na(s$diff_to_var[i])

# is_present: the variable or one of its derivatives occurs in an equation.
.ss_is_present <- function(s, v) {
  repeat {
    if (length(s$graph$badjlist[[v]])) return(TRUE)
    v <- s$var_to_diff[v]
    if (is.na(v)) return(FALSE)
  }
}

.ss_error <- function(class, msg) {
  stop(structure(class = c(class, "error", "condition"), list(message = msg, call = NULL)))
}

# ---------------------------------------------------------------------------
# rm_eqs_vars! and get_new_mm

# Remove equations and variables from the structure and the tearing state.
# Returns the old to new index maps (0 for removed entries). The linear
# subsystem matrix is not touched: callers apply .ss_get_new_mm().
.ss_rm_eqs_vars <- function(ts, eqs_rm, vars_rm) {
  s <- ts$structure; g <- s$graph
  ne <- .bg_nsrcs(g); nv <- .bg_ndsts(g)
  keep_e <- setdiff(seq_len(ne), eqs_rm)
  keep_v <- setdiff(seq_len(nv), vars_rm)
  o2n_eq <- integer(ne); o2n_eq[keep_e] <- seq_along(keep_e)
  o2n_var <- integer(nv); o2n_var[keep_v] <- seq_along(keep_v)
  remap <- function(h) {
    if (is.null(h)) return(NULL)
    BipartiteGraph(fadjlist = lapply(keep_e, function(e) {
      x <- o2n_var[h$fadjlist[[e]]]
      x[x > 0L]
    }), ndsts = length(keep_v))
  }
  s$graph <- remap(g)
  s$solvable_graph <- remap(s$solvable_graph)
  e2d <- s$eq_to_diff[keep_e]
  e2d[!is.na(e2d)] <- o2n_eq[e2d[!is.na(e2d)]]
  if (any(e2d == 0L, na.rm = TRUE)) stop("internal error: removed a differentiated equation", call. = FALSE)
  v2d <- s$var_to_diff[keep_v]
  v2d[!is.na(v2d)] <- o2n_var[v2d[!is.na(v2d)]]
  v2d[!is.na(v2d) & v2d == 0L] <- NA_integer_
  inverse <- function(x) {
    out <- rep(NA_integer_, length(x))
    hit <- which(!is.na(x))
    out[x[hit]] <- hit
    out
  }
  s$eq_to_diff <- e2d; s$diff_to_eq <- inverse(e2d)
  s$var_to_diff <- v2d; s$diff_to_var <- inverse(v2d)
  s$state_priorities <- s$state_priorities[keep_v]
  s$canonical_ranks <- s$canonical_ranks[keep_v]
  ts$fullvars <- ts$fullvars[keep_v]
  ts$always_present <- ts$always_present[keep_v]
  ts$eqs <- ts$eqs[keep_e]
  ts$original_eqs <- ts$original_eqs[keep_e]
  list(eq = o2n_eq, var = o2n_var)
}

# The coefficient matrix after removing equations and variables. A removed
# variable with an alias (target index and sign, 0 for none) contributes its
# coefficient to the target; a row with a removed variable that has no alias
# is no longer an integer linear equation and is dropped.
# Departure from get_new_mm: Julia carries the coefficient to the alias
# target without the sign of the alias (v ~ -target); the sign is applied
# here so that the row stays equal to the substituted equation.
.ss_get_new_mm <- function(mm, o2n_eq, o2n_var, alias_target = integer(length(o2n_var)),
                           alias_sign = rep(1, length(o2n_var))) {
  nz <- integer(); cols <- list(); vals <- list()
  for (i in seq_along(mm$nzrows)) {
    eq <- mm$nzrows[i]
    if (o2n_eq[eq] == 0L) next
    rc <- mm$row_cols[[i]]; rv <- mm$row_vals[[i]]
    nc <- o2n_var[rc]
    gone <- nc == 0L
    if (any(gone)) {
      tg <- alias_target[rc[gone]]
      if (any(tg == 0L)) next
      nc[gone] <- o2n_var[tg]
      rv[gone] <- rv[gone] * alias_sign[rc[gone]]
      if (any(nc == 0L)) stop("internal error: a variable is aliased to a removed variable", call. = FALSE)
    }
    if (length(nc)) {
      u <- sort.int(unique(nc))
      rv <- vapply(u, function(cc) sum(rv[nc == cc]), 0)
      nc <- u[rv != 0]; rv <- rv[rv != 0]
    }
    k <- length(nz) + 1L
    nz[k] <- o2n_eq[eq]; cols[[k]] <- as.integer(nc); vals[[k]] <- rv
  }
  list(nparentrows = sum(o2n_eq > 0L), ncols = sum(o2n_var > 0L), nzrows = nz,
       row_cols = cols, row_vals = vals)
}

# ---------------------------------------------------------------------------
# Fraction free Bareiss elimination

.BAREISS_MAX <- 2^53

.bareiss_guard <- function(x) {
  if (any(abs(x) >= .BAREISS_MAX))
    stop("integer overflow in the Bareiss factorization of the linear subsystem", call. = FALSE)
  x
}

# find_first_linear_variable: the pivot among rows k..n whose column passes
# `mask` (NULL for any column), preferring lower priority, then sparser rows.
.bareiss_find_pivot <- function(M, k, n, mask, prio) {
  ci <- 0L; cv <- 0L; cval <- 0; cnnz <- 0L
  for (i in k:n) {
    vs <- M$row_cols[[i]]
    nnz <- length(vs)
    if (cv != 0L && is.null(prio) && nnz >= cnnz) next
    for (j in seq_along(vs)) {
      v <- vs[j]
      if (!is.null(mask) && !mask[v]) next
      if (cv == 0L || (is.null(prio) && nnz < cnnz) ||
          (!is.null(prio) && (prio[v] < prio[cv] || (prio[v] == prio[cv] && nnz < cnnz)))) {
        ci <- i; cv <- v; cval <- M$row_vals[[i]][j]; cnnz <- nnz
      }
    }
  }
  if (cv == 0L) NULL else list(i = ci, v = cv, val = cval)
}

# bareiss_update_virtual_colswap_mtk!: eliminate column `vpivot` from the
# rows below k. Rows without an entry in that column are only rescaled by
# pivot / last_pivot, which is skipped when the two agree up to sign.
.bareiss_update <- function(M, k, n, vpivot, pivot, last_pivot) {
  if (k >= n) return(M)
  pivot_equal <- abs(pivot) == abs(last_pivot)
  kvars <- M$row_cols[[k]]; kcoeffs <- M$row_vals[[k]]
  for (ei in (k + 1L):n) {
    ivars <- M$row_cols[[ei]]; icoeffs <- M$row_vals[[ei]]
    pos <- match(vpivot, ivars)
    if (is.na(pos)) {
      if (pivot_equal) next
      coeff <- 0
    } else coeff <- icoeffs[pos]
    allv <- sort.int(union(ivars, kvars))
    allv <- allv[allv != vpivot]
    a <- icoeffs[match(allv, ivars)]; a[is.na(a)] <- 0
    b <- kcoeffs[match(allv, kvars)]; b[is.na(b)] <- 0
    num <- .bareiss_guard(.bareiss_guard(pivot * a) - .bareiss_guard(coeff * b))
    ci <- num / last_pivot
    if (any(ci != round(ci))) stop("internal error: inexact division in the Bareiss factorization", call. = FALSE)
    keep <- ci != 0
    M$row_cols[[ei]] <- as.integer(allv[keep])
    M$row_vals[[ei]] <- ci[keep]
  }
  M
}

# do_bareiss!: tiered pivoting. Pivots are taken first among the linear
# variables (rank1), then among the highest derivatives (rank2), then
# anywhere (rank3).
.do_bareiss <- function(M, lin, high, prio) {
  n <- length(M$nzrows)
  prev <- 1; rank1 <- NA_integer_; rank2 <- NA_integer_; rank3 <- n
  pivots <- integer()
  for (k in seq_len(n)) {
    r <- NULL
    if (is.na(rank1)) {
      r <- .bareiss_find_pivot(M, k, n, lin, prio)
      if (is.null(r)) rank1 <- k - 1L
    }
    if (is.null(r) && is.na(rank2)) {
      r <- .bareiss_find_pivot(M, k, n, high, prio)
      if (is.null(r)) rank2 <- k - 1L
    }
    if (is.null(r)) r <- .bareiss_find_pivot(M, k, n, NULL, prio)
    if (is.null(r)) {
      rank3 <- k - 1L
      break
    }
    pivots <- c(pivots, r$v)
    if (r$i != k) {
      sw <- c(k, r$i); ws <- c(r$i, k)
      M$nzrows[sw] <- M$nzrows[ws]
      M$row_cols[sw] <- M$row_cols[ws]
      M$row_vals[sw] <- M$row_vals[ws]
    }
    M <- .bareiss_update(M, k, n, r$v, r$val, prev)
    prev <- r$val
  }
  if (is.na(rank2)) rank2 <- rank3
  if (is.na(rank1)) rank1 <- rank2
  list(M = M, rank1 = rank1, rank2 = rank2, rank3 = rank3, pivots = pivots)
}

# Dense Bareiss with full pivoting (bareiss! with find_pivot_any) as used by
# bareiss.nullspace: returns the rank and the column order, whose first
# `rank` entries are the pivot columns.
.bareiss_rank <- function(A) {
  n <- nrow(A); mc <- ncol(A)
  col_order <- seq_len(mc); prev <- 1; rank <- 0L
  for (k in seq_len(min(n, mc))) {
    sub <- A[k:n, k:mc, drop = FALSE]
    p <- which(sub != 0)[1L]
    if (is.na(p)) break
    pi <- (p - 1L) %% nrow(sub) + k
    pj <- (p - 1L) %/% nrow(sub) + k
    pivot <- A[pi, pj]
    if (pj != k) {
      A[, c(k, pj)] <- A[, c(pj, k)]
      col_order[c(k, pj)] <- col_order[c(pj, k)]
    }
    if (pi != k) A[c(k, pi), ] <- A[c(pi, k), ]
    if (k < n) {
      if (k < mc) {
        r <- (k + 1L):n; cc <- (k + 1L):mc
        num <- .bareiss_guard(.bareiss_guard(A[r, cc, drop = FALSE] * pivot) -
                                .bareiss_guard(outer(A[r, k], A[k, cc])))
        A[r, cc] <- num / prev
      }
      A[(k + 1L):n, k] <- 0
    }
    prev <- pivot; rank <- k
  }
  list(rank = rank, col_order = col_order)
}

# Rows ordered by sparsity, then columns, then coefficients, then equation.
.mm_row_order <- function(mm, rows) {
  if (length(rows) <= 1L) return(rows)
  len <- lengths(mm$row_cols[rows])
  w <- max(len)
  if (w == 0L) return(rows[order(mm$nzrows[rows])])
  pad <- function(x) vapply(x, function(z) c(z, rep(-Inf, w - length(z))), numeric(w))
  cm <- matrix(pad(mm$row_cols[rows]), nrow = w)
  vm <- matrix(pad(mm$row_vals[rows]), nrow = w)
  keys <- c(list(len), lapply(seq_len(w), function(i) cm[i, ]), lapply(seq_len(w), function(i) vm[i, ]),
            list(mm$nzrows[rows]))
  rows[do.call(order, keys)]
}

# aag_bareiss!: Bareiss on each connected block of the linear subsystem.
.aag_bareiss <- function(s, mm) {
  g <- s$graph
  nv <- length(s$var_to_diff)
  is_alg <- is.na(s$var_to_diff) & is.na(s$diff_to_var)
  lin <- is_alg
  high <- .highest_diff_variables(s)
  lineq <- logical(.bg_nsrcs(g)); lineq[mm$nzrows] <- TRUE
  for (i in seq_len(.bg_nsrcs(g))) {
    nb <- g$fadjlist[[i]]
    if (lineq[i] && all(is_alg[nb])) next
    lin[nb] <- FALSE
  }
  solvable_variables <- which(lin)
  big <- max(s$canonical_ranks, 0L) + 1
  prio <- s$state_priorities * big + s$canonical_ranks

  # rows sharing a column belong to the same block
  nrows <- length(mm$nzrows)
  parent <- seq_len(nrows)
  find <- function(x) {
    while (parent[x] != x) x <- parent[x]
    x
  }
  col_row <- integer(nv)
  for (i in seq_len(nrows)) for (cc in mm$row_cols[[i]]) {
    j <- col_row[cc]
    if (j == 0L) col_row[cc] <- i
    else {
      ra <- find(i); rb <- find(j)
      if (ra != rb) parent[max(ra, rb)] <- min(ra, rb)
    }
  }
  roots <- vapply(seq_len(nrows), find, 1L)
  out <- list(nparentrows = mm$nparentrows, ncols = mm$ncols, nzrows = integer(),
              row_cols = list(), row_vals = list())
  t1 <- integer(); t2 <- integer(); t3 <- integer()
  for (r in unique(roots)) {
    rows <- .mm_row_order(mm, which(roots == r))
    sub <- list(nzrows = mm$nzrows[rows], row_cols = mm$row_cols[rows], row_vals = mm$row_vals[rows])
    b <- .do_bareiss(sub, lin, high, prio)
    out$nzrows <- c(out$nzrows, b$M$nzrows)
    out$row_cols <- c(out$row_cols, b$M$row_cols)
    out$row_vals <- c(out$row_vals, b$M$row_vals)
    t1 <- c(t1, b$pivots[seq_len(b$rank1)])
    t2 <- c(t2, b$pivots[b$rank1 + seq_len(b$rank2 - b$rank1)])
    t3 <- c(t3, b$pivots[b$rank2 + seq_len(b$rank3 - b$rank2)])
  }
  list(mm = out, solvable_variables = solvable_variables, rank1 = length(t1),
       rank2 = length(t1) + length(t2), rank3 = length(t1) + length(t2) + length(t3),
       pivots = c(t1, t2, t3))
}

# structural_singularity_removal!: replace the integer linear equations by
# their Bareiss reduced rows. Variables that occur only in such equations
# and are not pivots are underconstrained; they are returned for reporting.
.structural_singularity_removal <- function(ts, mm) {
  s <- ts$structure
  if (!length(mm$nzrows)) return(list(mm = mm, underconstrained = integer()))
  b <- .aag_bareiss(s, mm)
  mm <- b$mm
  for (i in seq_along(mm$nzrows)) {
    .bg_set_neighbors(s$graph, mm$nzrows[i], mm$row_cols[[i]])
    if (!is.null(s$solvable_graph)) .bg_set_neighbors(s$solvable_graph, mm$nzrows[i], mm$row_cols[[i]])
  }
  list(mm = mm, underconstrained = setdiff(b$solvable_variables, b$pivots[seq_len(b$rank1)]))
}

# ---------------------------------------------------------------------------
# Consistency check

.ss_is_unused_var <- function(ts, v) !ts$always_present[v] && !length(ts$structure$graph$badjlist[[v]])

.ss_print_array <- function(xs) paste0(" ", vapply(xs, .fmt, ""), collapse = "\n")

# singular_check: variables that stay unmatched when every derivative
# relation counts as one more equation.
.singular_check <- function(ts) {
  s <- ts$structure; g <- s$graph
  edges <- which(!is.na(s$var_to_diff))
  ext <- BipartiteGraph(fadjlist = c(g$fadjlist, lapply(edges, function(v) c(v, s$var_to_diff[v]))),
                        ndsts = .bg_ndsts(g))
  m <- .maximal_matching(ext)
  bad <- Filter(function(v) m[v] == .UNASSIGNED && !.ss_is_unused_var(ts, v), seq_len(.bg_ndsts(g)))
  ts$fullvars[bad]
}

# check_consistency: error for unbalanced or structurally singular systems
# (or return FALSE when `nothrow`).
.check_consistency <- function(ts, orig_inputs = list(), nothrow = FALSE) {
  s <- ts$structure; g <- s$graph
  neqs <- sum(lengths(g$fadjlist) > 0L)
  highest <- which(.highest_diff_variables(s))
  n_highest <- sum(!vapply(highest, function(v) .ss_is_unused_var(ts, v), TRUE))
  balanced <- n_highest == neqs
  if (neqs > 0L && !balanced) {
    if (nothrow) return(FALSE)
    m <- .maximal_matching(g, dstfilter = is.na(s$var_to_diff))
    head <- sprintf("The system is unbalanced. There are %d highest order derivative variables and %d equations.\n",
                    n_highest, neqs)
    if (n_highest < neqs) {
      bad <- which(.inv_matching(m, .bg_nsrcs(g)) == 0L)
      .ss_error("ExtraEquationsSystemException", paste0(
        "ExtraEquationsSystemException: ", head,
        "More equations than variables, here are the potential extra equation(s):\n",
        .ss_print_array(ts$eqs[bad]),
        "\nNote that the process of determining extra equations is a best-effort heuristic. ",
        "The true extra equations are dependent on the model and may not be in this list."))
    }
    bad <- ts$fullvars[which(m[seq_len(.bg_ndsts(g))] == .UNASSIGNED)]
    unset <- bad[.keys(bad) %in% .keys(orig_inputs)]
    extra <- if (length(unset)) paste0("In particular, the unset input(s) are:\n", .ss_print_array(unset),
                                       "\nThe rest of potentially unset variable(s) are:\n") else ""
    .ss_error("ExtraVariablesSystemException", paste0(
      "ExtraVariablesSystemException: ", head,
      "More variables than equations, here are the potential extra variable(s):\n", extra,
      .ss_print_array(bad),
      "\nNote that the process of determining extra variables is a best-effort heuristic. ",
      "The true extra variables are dependent on the model and may not be in this list."))
  }
  bad <- .singular_check(ts)
  if (length(bad) || !balanced) {
    if (nothrow) return(FALSE)
    .ss_error("InvalidSystemException", paste0(
      "InvalidSystemException: The system is structurally singular! ",
      "Here are the problematic variables: \n", .ss_print_array(bad)))
  }
  TRUE
}

# ---------------------------------------------------------------------------
# Pantelides index reduction

# computed_highest_diff_variables: the highest differentiated variable of
# each chain that occurs in the system. `varfilter` is a logical vector.
.highest_diff_variables <- function(s, varfilter = NULL) {
  g <- s$graph; v2d <- s$var_to_diff; d2v <- s$diff_to_var
  n <- length(v2d)
  wl <- logical(n)
  for (var in seq_len(n)) {
    if (!is.null(varfilter) && !varfilter[var]) next
    if (is.na(v2d[var]) && !wl[var]) {
      while (!length(g$badjlist[[var]])) {
        lower <- d2v[var]
        if (is.na(lower)) break
        var <- lower
      }
      wl[var] <- TRUE
    }
  }
  for (var in seq_len(n)) {
    if (!wl[var]) next
    hi <- var
    while (!is.na(hi <- v2d[hi])) {
      if (wl[hi]) {
        wl[var] <- FALSE
        break
      }
    }
  }
  wl
}

# pantelides!: differentiate equations until every equation can be matched
# to a highest order derivative. Returns the variable to equation matching.
.pantelides <- function(ts, finalize = TRUE, maxiters = 8000L, eqfilter = NULL, varfilter = NULL) {
  s <- ts$structure
  neqs <- .bg_nsrcs(s$graph)
  nvars <- length(s$var_to_diff)
  st <- new.env(parent = emptyenv())
  st$m <- integer(nvars)
  neqs0 <- neqs
  eq_ok <- function(e) is.null(eqfilter) || eqfilter[e]
  nnonempty <- sum(vapply(seq_len(neqs0), function(e)
    length(s$graph$fadjlist[[e]]) > 0L && is.na(s$eq_to_diff[e]) && eq_ok(e), TRUE))
  wl <- .highest_diff_variables(s, varfilter)
  if (nnonempty > sum(wl))
    .ss_error("InvalidSystemException", "InvalidSystemException: System is structurally singular")
  for (k in seq_len(neqs0)) {
    eq <- k
    if (!eq_ok(k) || !is.na(s$eq_to_diff[eq]) || !length(s$graph$fadjlist[[eq]])) next
    pathfound <- FALSE
    for (iter in seq_len(maxiters)) {
      st$dcolor <- logical(nvars)
      st$scolor <- logical(neqs)
      pathfound <- .augment(st, s$graph$fadjlist, eq, wl)
      if (pathfound) break
      vcol <- which(st$dcolor); ecol <- which(st$scolor)
      for (var in vcol) {
        if (is.na(s$var_to_diff[var])) {
          nvars <- nvars + 1L
          .ts_var_derivative(ts, var)
          st$m <- c(st$m, .UNASSIGNED)
          wl <- c(wl, FALSE)
        }
        wl[var] <- FALSE
        wl[s$var_to_diff[var]] <- TRUE
      }
      for (e in ecol) {
        neqs <- neqs + 1L
        .ts_eq_derivative(ts, e)
      }
      for (var in vcol) st$m[s$var_to_diff[var]] <- s$eq_to_diff[st$m[var]]
      eq <- s$eq_to_diff[eq]
    }
    if (!pathfound)
      stop("maxiters = ", maxiters, " reached in the Pantelides algorithm; the index of the system is ",
           "too high or the system is singular", call. = FALSE)
  }
  if (finalize) st$m[which(!wl)] <- .UNASSIGNED
  st$m
}

# ---------------------------------------------------------------------------
# Tearing

# update_full_var_eq_matching!: copy the torn matching of one component into
# the full matching and give the unsolved variables the remaining equations.
.update_full_matching <- function(g, fullm, m, vars, eqs, vf) {
  for (var in vars) {
    if (!vf[var]) next
    fullm[var] <- m[var]
    if (m[var] > 0L) eqs <- setdiff(eqs, m[var])
  }
  for (var in vars) {
    if (!length(eqs)) break
    if (!vf[var] || fullm[var] > 0L) next
    hit <- eqs[vapply(eqs, function(e) .bg_has_edge(g, e, var), TRUE)]
    eq <- if (length(hit)) hit[1L] else eqs[1L]
    fullm[var] <- eq
    eqs <- setdiff(eqs, eq)
  }
  fullm
}

.lex_less <- function(a, b) {
  d <- which(a != b)
  length(d) && a[d[1L]] < b[d[1L]]
}

# carpanzano_tear_scc!: algorithm A1 of Carpanzano (2000) on one strongly
# connected component. Variables left unassigned are the tear variables.
.carpanzano_tear_scc <- function(s, m, av, ae) {
  g <- s$graph; sg <- s$solvable_graph
  tearkey <- function(v) c(-s$state_priorities[v], s$canonical_ranks[v])
  av <- av[vapply(av, function(v) any(sg$badjlist[[v]] %in% ae), TRUE)]
  while (length(av)) {
    # an equation incident on, and solvable for, exactly one active variable.
    # (Julia first looks for one whose variable is a kept derivative, but its
    # find_single_solvable_eq! ignores that condition, so one scan suffices.)
    hit <- 0L
    for (ieq in ae) {
      nb <- intersect(g$fadjlist[[ieq]], av)
      if (length(nb) != 1L || !(nb %in% sg$fadjlist[[ieq]])) next
      m[nb] <- ieq
      hit <- ieq
      av <- setdiff(av, nb)
      break
    }
    if (hit) {
      ae <- setdiff(ae, hit)
      next
    }
    # Heuristic 1: among the equations with the fewest active variables, a
    # variable that is present but not solvable becomes a tear variable.
    cand <- 0L
    if (length(ae)) {
      cnt <- vapply(ae, function(ieq) sum(g$fadjlist[[ieq]] %in% av), 0L)
      for (ieq in ae[cnt == min(cnt)]) {
        ns <- intersect(setdiff(g$fadjlist[[ieq]], sg$fadjlist[[ieq]]), av)
        for (ivar in ns) if (cand == 0L || .lex_less(tearkey(ivar), tearkey(cand))) cand <- ivar
      }
    }
    if (cand != 0L) {
      av <- setdiff(av, cand)
      next
    }
    # Heuristic 2: the variable in the most equations, solvable in the fewest.
    best <- 0L; bcnt <- 0L; bsol <- 0L; bkey <- NULL
    for (ivar in av) {
      cnt <- sum(g$badjlist[[ivar]] %in% ae)
      sol <- sum(sg$badjlist[[ivar]] %in% ae)
      key <- tearkey(ivar)
      if (best == 0L || cnt > bcnt || (cnt == bcnt && (sol < bsol || (sol == bsol && .lex_less(key, bkey))))) {
        best <- ivar; bcnt <- cnt; bsol <- sol; bkey <- key
      }
    }
    av <- setdiff(av, best)
  }
  m
}

# CarpanzanoTearing: a maximal matching gives the strongly connected
# components; each is torn separately.
# ponytail: the exact Bareiss matching of pure integer linear components
# (exact_scc_matching!) is not ported. In the pinned StateSelection commit
# its `linear_rewrite` result has no consumer, so the graph and the symbolic
# equations would diverge. Port it together with a consumer if the structural
# heuristics ever pick a singular assignment for such a block.
.carpanzano_tearing <- function(s, varfilter = NULL, eqfilter = NULL) {
  g <- s$graph
  vf <- varfilter %||% rep(TRUE, .bg_ndsts(g))
  m <- .maximal_matching(g, eqfilter, vf)
  fullm <- m
  var_sccs <- .find_var_sccs(g, m)
  for (vars in var_sccs) {
    av <- integer(); ae <- integer()
    for (var in vars) {
      if (vf[var]) {
        av <- c(av, var)
        if (m[var] > 0L) ae <- c(ae, m[var])
      }
      m[var] <- .UNASSIGNED
    }
    m <- .carpanzano_tear_scc(s, m, av, ae)
    fullm <- .update_full_matching(g, fullm, m, vars, ae, vf)
  }
  list(var_eq_matching = m, full_var_eq_matching = fullm, var_sccs = var_sccs)
}

# ModiaTearing (tear_graph_modia, tearEquations!): assign each equation of a
# component to a solvable variable unless that closes an algebraic loop.
.modia_tearing <- function(s, isder = NULL, varfilter = NULL, eqfilter = NULL) {
  g <- s$graph; sg <- s$solvable_graph
  nd <- .bg_ndsts(g)
  vf <- varfilter %||% rep(TRUE, nd)
  m <- .maximal_matching(g, eqfilter, vf)
  fullm <- m
  var_sccs <- .find_var_sccs(g, m)
  matched <- unlist(lapply(var_sccs, function(vs) m[vs[vf[vs]]]))
  overdetermined <- length(setdiff(seq_len(.bg_nsrcs(g)), matched[matched > 0L])) > 0L
  ict <- integer(nd)     # the matching of the cycle tracker
  # adding the edges u -> vj for the other variables u of `eq` closes a
  # cycle when vj already reaches one of them
  # ponytail: a full reachability search per attempt instead of Julia's
  # incremental cycle tracker; fine for the component sizes R can handle.
  try_assign <- function(vj, eq) {
    others <- setdiff(g$fadjlist[[eq]], vj)
    inv <- .inv_matching(ict, .bg_nsrcs(g))
    seen <- logical(nd); stack <- vj
    while (length(stack)) {
      x <- stack[length(stack)]; stack <- stack[-length(stack)]
      if (seen[x]) next
      seen[x] <- TRUE
      w <- inv[g$badjlist[[x]]]
      w <- w[w > 0L]
      stack <- c(stack, w[!seen[w]])
    }
    if (any(seen[others])) return(FALSE)
    ict[vj] <<- eq
    TRUE
  }
  try_vars <- function(vs, active, eq, cond = NULL) {
    for (vj in vs) {
      if (!active[vj] || ict[vj] != .UNASSIGNED || (!is.null(cond) && !cond(vj))) next
      if (try_assign(vj, eq)) return(TRUE)
    }
    FALSE
  }
  for (vars in var_sccs) {
    ieqs <- integer(); active <- logical(nd)
    for (var in vars) {
      if (vf[var]) {
        active[var] <- TRUE
        if (m[var] > 0L) ieqs <- c(ieqs, m[var])
      }
      m[var] <- .UNASSIGNED
    }
    for (single in c(TRUE, FALSE)) for (eq in ieqs) {
      vs <- sg$fadjlist[[eq]]
      if (xor(length(vs) == 1L, single)) next
      if (length(vs) > 1L) vs <- vs[order(s$state_priorities[vs])]
      if (!is.null(isder)) {
        has_der <- FALSE
        cond <- function(v) {
          r <- isder(v)
          if (r) has_der <<- TRUE
          r
        }
        try_vars(vs, active, eq, cond)
        if (has_der) next
      }
      try_vars(vs, active, eq)
    }
    for (var in which(active)) m[var] <- ict[var]
    fullm <- .update_full_matching(g, fullm, m, vars, ieqs, vf)
    if (!overdetermined) ict[vars] <- .UNASSIGNED
  }
  list(var_eq_matching = m, full_var_eq_matching = fullm, var_sccs = var_sccs)
}

# DummyDerivativeTearing: tear with the derivative aware filters and mark the
# selected states. Variables whose derivative stays a real derivative are
# states and cannot be eliminated.
.dummy_derivative_tearing <- function(s, dummy_derivatives = integer(), alg = "carpanzano") {
  nd <- .bg_ndsts(s$graph)
  is_dd <- logical(nd); is_dd[dummy_derivatives] <- TRUE
  is_some_diff <- function(v) !is_dd[v] && .ss_is_present(s, v)
  can_eliminate <- vapply(seq_len(nd), function(v) {
    dv <- s$var_to_diff[v]
    is.na(dv) || !is_some_diff(dv)
  }, TRUE)
  isder <- function(v) !is.na(s$diff_to_var[v]) && is_some_diff(v)
  res <- if (identical(alg, "modia")) .modia_tearing(s, isder, can_eliminate)
         else .carpanzano_tearing(s, can_eliminate)
  for (v in seq_len(nd)) {
    if (!.ss_is_present(s, v)) next
    dv <- s$var_to_diff[v]
    if (is.na(dv) || !is_some_diff(dv)) next
    res$var_eq_matching[v] <- .SELECTED_STATE
  }
  res$can_eliminate <- can_eliminate
  res
}

# dummy_derivative_graph!: choose, level by level, which of the highest
# derivatives are dummy derivatives (algebraic variables), preferring the
# variables with the lowest state priority, then tear.
# `jac(eqs, vars)` returns the integer Jacobian of the equations with respect
# to the variables, or NULL when it is not an integer matrix.
.dummy_derivative_graph <- function(ts, var_eq_matching, jac = NULL, tearing_alg = "carpanzano") {
  s <- ts$structure; g <- s$graph
  d2e <- s$diff_to_eq; d2v <- s$diff_to_var; v2d <- s$var_to_diff
  sp <- s$state_priorities
  extended_sp <- function(var) {
    while (!is.na(v2d[var])) var <- v2d[var]
    minp <- 0; maxp <- 0
    repeat {
      p <- sp[var]
      maxp <- max(maxp, p); minp <- min(minp, p)
      var <- d2v[var]
      if (is.na(var)) break
    }
    if (minp < 0) minp else maxp
  }
  neqs <- .bg_nsrcs(g); nvars <- .bg_ndsts(g)
  dummy <- integer()
  dd_sccs <- list(); dd_prio <- list()
  for (vars0 in .find_var_sccs(g, var_eq_matching)) {
    eqs <- integer(); vars <- integer()
    for (var in vars0) {
      eq <- var_eq_matching[var]
      if (eq <= 0L) next
      if (!is.na(d2e[eq])) eqs <- c(eqs, eq)
      if (!is.na(v2d[var])) stop("internal error: invalid component in dummy derivative selection", call. = FALSE)
      if (!is.na(d2v[var]) && .ss_is_present(s, var)) vars <- c(vars, var)
    }
    if (!length(eqs)) next
    J <- if (is.null(jac)) NULL else jac(eqs, vars)
    isfirst <- TRUE
    next_eq <- integer(); next_var <- integer()
    repeat {
      nrows <- length(eqs)
      if (!nrows) break
      if (isfirst) {
        spv <- vapply(vars, extended_sp, 0)
        perm <- order(spv, s$canonical_ranks[vars])
        vars <- vars[perm]
        if (!is.null(J)) J <- J[, perm, drop = FALSE]
        dd_sccs[[length(dd_sccs) + 1L]] <- vars
        dd_prio[[length(dd_prio) + 1L]] <- spv[perm]
      }
      if (!is.null(J)) {
        if (!isfirst) J <- J[next_eq, next_var, drop = FALSE]
        b <- .bareiss_rank(J)
        rank <- b$rank
        dummy <- c(dummy, vars[b$col_order[seq_len(rank)]])
      } else {
        st <- new.env(parent = emptyenv())
        st$m <- integer(max(nvars, neqs))
        in_eqs <- logical(neqs); in_eqs[eqs] <- TRUE
        rank <- 0L
        for (var in vars) {
          st$dcolor <- logical(neqs)
          if (!.augment(st, g$badjlist, var, in_eqs)) next
          dummy <- c(dummy, var)
          rank <- rank + 1L
          if (rank == nrows) break
        }
      }
      if (rank != nrows) warning("The DAE system is singular!", call. = FALSE)
      # descend one differentiation level; only equations and variables that
      # are still differentiated take part in the next round
      next_eq <- integer(); next_var <- integer()
      new_eqs <- integer(); new_vars <- integer()
      for (i in seq_along(eqs)) {
        ieq <- d2e[eqs[i]]
        if (is.na(ieq) || is.na(d2e[ieq])) next
        next_eq <- c(next_eq, i); new_eqs <- c(new_eqs, ieq)
      }
      for (i in seq_along(vars)) {
        ivar <- d2v[vars[i]]
        if (is.na(ivar) || is.na(d2v[ivar])) next
        next_var <- c(next_var, i); new_vars <- c(new_vars, ivar)
      }
      eqs <- new_eqs; vars <- new_vars
      isfirst <- FALSE
    }
  }
  n_diff_eqs <- sum(!is.na(d2e))
  if (n_diff_eqs != length(dummy))
    warning("The number of dummy derivatives (", length(dummy),
            ") does not match the number of differentiated equations (", n_diff_eqs, ").", call. = FALSE)
  res <- .dummy_derivative_tearing(s, dummy, tearing_alg)
  res$dummy_derivatives <- dummy
  res$ddsummary <- list(var_sccs = dd_sccs, state_priority = dd_prio)
  res
}

# ---------------------------------------------------------------------------
# Trivial tearing

# trivial_tearing!: equations already written as `v ~ f(...)` whose left-hand
# side occurs nowhere else become observed equations before any other pass.
.trivial_tearing <- function(ts) {
  s <- ts$structure; g <- s$graph
  cand <- .ts_possibly_explicit_equations(ts)
  trivial <- integer(); matched <- integer()
  repeat {
    added <- FALSE
    for (cd in cand) {
      i <- cd[1L]; vari <- cd[2L]
      if (i %in% trivial) next
      if (s$state_priorities[vari] > 0L) next
      if (vari %in% matched) next
      if (!is.na(s$var_to_diff[vari]) || !is.na(s$diff_to_var[vari])) next
      eqidxs <- setdiff(g$badjlist[[vari]], trivial)
      if (length(eqidxs) != 1L || eqidxs != i) next
      valid <- TRUE
      for (v in g$fadjlist[[i]]) {
        if (v == vari || v %in% matched) next
        if (sum(!(g$badjlist[[v]] %in% trivial)) <= 1L || !is.na(s$diff_to_var[v])) {
          valid <- FALSE
          break
        }
      }
      if (!valid) next
      added <- TRUE
      trivial <- c(trivial, i); matched <- c(matched, vari)
    }
    if (!added) break
  }
  # the torn order is a topological order of the new observed equations
  ts$additional_observed <- c(ts$additional_observed, ts$original_eqs[trivial])
  .ss_rm_eqs_vars(ts, trivial, matched)
  invisible(ts)
}
