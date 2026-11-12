# Ported from StateSelection.jl (https://github.com/JuliaComputing/StateSelection.jl),
# commit 29b3c90. Portions are licensed under the GNU AGPL, version 3; see
# inst/COPYRIGHTS and inst/STATESELECTION-AGPL-LICENSE.
#
# Inline solves of torn algebraic blocks, after
# StateSelection.jl/lib/ModelingToolkitTearing/src/reassemble.jl
# (get_linear_scc_linsol, __reduce_linear_system!, safe_ldiv, numeric_ldiv!,
# the inline branch of generate_system_equations! and the
# `inline_linear_sccs` and `analytical_linear_scc_limit` options of
# DefaultReassembleAlgorithm).
#
# With `inline_linear_sccs = TRUE` a strongly connected block of algebraic
# equations that is linear in its variables is solved inside the generated
# right-hand side instead of being left to the solver as residuals: blocks of
# at most `analytical_linear_scc_limit` variables (and every block of one)
# are solved symbolically by Gaussian elimination, larger ones by a numeric
# solve at run time (`inline_linsolve`, the counterpart of safe_ldiv). The
# rows that tearing already solves for a variable are eliminated first, so
# only the torn part of the block enters the linear system.
#
# Extension (not in ModelingToolkit): with `inline_nonlinear_sccs = TRUE` a
# torn block that is not linear becomes a nested Newton solve: its tearing
# variables are observed through a generated function `inline_nlsolve_<n>`
# that solves the block residuals for them, warm started from its last
# solution, so the system compiles to an ODE. The block must have a locally
# unique solution with a nonsingular Jacobian; blocks that contain
# derivatives stay residual equations.

# A one entry memo: the N calls for the components of one block share a
# single solve of A x = b. Arguments: component j, size n, A by columns, b.
.inline_linsolve_memo <- new.env(parent = emptyenv())

.inline_linsolve_fn <- function(j, n, ...) {
  v <- c(...)
  if (!identical(v, .inline_linsolve_memo$args)) {
    n <- as.integer(n)
    A <- matrix(v[seq_len(n * n)], n, n)
    b <- v[n * n + seq_len(n)]
    x <- tryCatch(solve(A, b), error = function(e) rep(NaN, n))
    assign("args", v, envir = .inline_linsolve_memo)
    assign("x", x, envir = .inline_linsolve_memo)
  }
  .inline_linsolve_memo$x[as.integer(j)]
}

.INLINE_LINSOLVE <- register_symbolic("inline_linsolve", .inline_linsolve_fn)

# Symbolic Gaussian elimination with pivoting on entries that are not zero
# (sym_lu followed by a solve). NULL when no pivot is found.
.ts_sym_linsolve <- function(A, b) {
  n <- length(b)
  for (k in seq_len(n)) {
    cand <- which(!vapply(k:n, function(r) .is_zero(A[[r, k]]), TRUE)) + k - 1L
    if (!length(cand)) return(NULL)
    num <- cand[vapply(cand, function(r) .is_num(A[[r, k]]), TRUE)]
    p <- if (length(num)) num[1L] else cand[1L]
    if (p != k) {
      tmp <- A[k, ]; A[k, ] <- A[p, ]; A[p, ] <- tmp
      tmp <- b[[k]]; b[[k]] <- b[[p]]; b[[p]] <- tmp
    }
    for (r in setdiff(k:n, k)) {
      if (.is_zero(A[[r, k]])) next
      f <- A[[r, k]] / A[[k, k]]
      for (c in k:n) A[[r, c]] <- A[[r, c]] - f * A[[k, c]]
      b[[r]] <- b[[r]] - f * b[[k]]
    }
  }
  x <- vector("list", n)
  for (k in rev(seq_len(n))) {
    acc <- b[[k]]
    for (c in setdiff(seq_len(n), seq_len(k))) acc <- acc - A[[k, c]] * x[[c]]
    x[[k]] <- acc / A[[k, k]]
  }
  x
}

# _check_allow_symbolic_parameter on a matrix: every entry may be divided by.
.ts_allow_matrix <- function(ts, A) {
  if (isTRUE(ts$opts$allow_symbolic)) return(TRUE)
  fk <- .keys(ts$fullvars)
  for (a in A) {
    if (.is_num(a)) next
    if (identical(ts$opts$allow_parameter, FALSE)) return(FALSE)
    if (.ts_maybe_zero(a, ts$maybe_zero_keys)) return(FALSE)
    if (any(.keys(.ts_search_vars(a)) %in% fk)) return(FALSE)
  }
  TRUE
}

# __reduce_linear_system!: rows solved by tearing are eliminated, their
# variables written in terms of the others. Returns the reduced A and b and
# the masks of the rows and columns kept.
.ts_reduce_linear_system <- function(A, b, eqm, varm, solved_var) {
  N <- length(b)
  alias <- list(); const <- list()
  for (i in seq_len(N)) {
    if (!eqm[i] || is.na(solved_var[i])) next
    v <- solved_var[i]
    eqm[i] <- FALSE; varm[v] <- FALSE
    cv <- A[[i, v]]
    coeffs <- lapply(seq_len(N), function(j) if (j == v) .S0() else -A[[i, j]] / cv)
    alias[[as.character(v)]] <- coeffs
    const[[as.character(v)]] <- b[[i]] / cv
  }
  # an eliminated variable written in terms of other eliminated ones is
  # expanded in dependency order (tearing leaves no cycles)
  for (it in seq_along(alias)) for (v in names(alias)) {
    co <- alias[[v]]; cst <- const[[v]]
    for (w in names(alias)) {
      wi <- as.integer(w)
      if (.is_zero(co[[wi]])) next
      cw <- co[[wi]]; co[[wi]] <- .S0()
      for (j in seq_len(N)) co[[j]] <- co[[j]] + cw * alias[[w]][[j]]
      cst <- cst + cw * const[[w]]
    }
    alias[[v]] <- co; const[[v]] <- cst
  }
  keep_r <- which(eqm); keep_c <- which(varm)
  Ar <- matrix(list(.S0()), length(keep_r), length(keep_c))
  br <- vector("list", length(keep_r))
  for (ii in seq_along(keep_r)) {
    i <- keep_r[ii]
    bi <- b[[i]]
    row <- A[i, ]
    for (v in names(alias)) {
      vi <- as.integer(v)
      if (.is_zero(row[[vi]])) next
      bi <- bi - row[[vi]] * const[[v]]
      for (j in keep_c) row[[j]] <- row[[j]] + row[[vi]] * alias[[v]][[j]]
    }
    br[[ii]] <- bi
    for (jj in seq_along(keep_c)) Ar[[ii, jj]] <- expand(row[[keep_c[jj]]])
  }
  list(A = Ar, b = br, eqm = eqm, varm = varm)
}

# get_linear_scc_linsol, and the nonlinear extension. Returns NULL (the block
# is generated as usual) or list(eqs, vars, sol): the solution sol[[k]] of
# variable vars[k], attributed to equation eqs[k].
.ts_inline_scc <- function(ts, escc, vscc, neweqs, m, fullm, inv, total_sub) {
  opts <- ts$opts
  lin <- isTRUE(opts$inline_linear_sccs); nl <- isTRUE(opts$inline_nonlinear_sccs)
  iv <- ts$iv
  if (is.null(iv) || !(lin || nl) || length(escc) != length(vscc) || !length(escc)) return(NULL)
  s <- ts$structure; g <- s$graph
  # a fully torn block needs no solve
  if (all(m[vscc] > 0L & !.ss_isdervar(s, vscc))) return(NULL)
  N <- length(escc)
  eqm <- rep(TRUE, N); varm <- rep(TRUE, N)
  for (k in seq_len(N)) if (.ts_is_irreducible(ts, ts$fullvars[[vscc[k]]])) {
    varm[k] <- FALSE
    j <- match(fullm[vscc[k]], escc)
    if (!is.na(j)) eqm[j] <- FALSE
  }
  vars <- ts$fullvars[vscc]
  resid <- lapply(escc, function(e) .ts_fixpoint_sub(neweqs[[e]]$rhs - neweqs[[e]]$lhs, total_sub, iv))
  out <- if (lin) .ts_inline_linear(ts, escc, vscc, vars, resid, eqm, varm, inv) else NULL
  if (is.null(out) && nl) out <- .ts_inline_nonlinear(ts, escc, vscc, vars, resid, eqm, varm, inv)
  out
}

.ts_inline_linear <- function(ts, escc, vscc, vars, b, eqm, varm, inv) {
  g <- ts$structure$graph
  N <- length(escc)
  A <- matrix(list(.S0()), N, N)
  bk <- .keys(vars[varm])
  for (vi in which(varm)) for (ei in which(eqm)) {
    if (!.bg_has_edge(g, escc[ei], vscc[vi])) next
    le <- .st_linear_expansion(b[[ei]], vars[[vi]])
    # a coefficient containing a variable of the block is not linear
    if (!le$islinear || any(.keys(.ts_search_vars(le$a)) %in% bk)) return(NULL)
    if (!.is_zero(le$a)) A[[ei, vi]] <- le$a
    b[[ei]] <- le$b
  }
  # the remainder must be free of the block's variables
  if (any(vapply(b[eqm], function(z) any(.keys(.ts_search_vars(z)) %in% bk), TRUE))) return(NULL)
  b <- lapply(b, function(z) -z)
  solved_var <- match(inv[escc], vscc)
  red <- .ts_reduce_linear_system(A, b, eqm, varm, solved_var)
  Nr <- length(red$b)
  if (!Nr) return(NULL)
  limit <- ts$opts$analytical_linear_scc_limit %||% 2L
  sol <- NULL
  if (Nr == 1L || (Nr <= limit && .ts_allow_matrix(ts, red$A))) sol <- .ts_sym_linsolve(red$A, red$b)
  if (is.null(sol)) {
    args <- c(unclass(red$A), red$b)
    sol <- lapply(seq_len(Nr), function(j) do.call(.INLINE_LINSOLVE, c(list(j, Nr), args)))
  }
  list(eqs = escc[red$eqm], vars = vscc[red$varm], sol = sol)
}

.mtk$inline_nl_counter <- 0L

.ts_inline_nonlinear <- function(ts, escc, vscc, vars, resid, eqm, varm, inv) {
  s <- ts$structure
  if (any(!eqm) || any(.ss_isdervar(s, vscc))) return(NULL)
  zi <- which(!(vscc %in% inv[escc][inv[escc] > 0L]))
  ri <- which(inv[escc] <= 0L)
  if (!length(zi) || length(zi) != length(ri)) return(NULL)
  Z <- vars[zi]
  # the equations tearing solves become local assignments of the residual
  solved <- list()
  for (k in which(inv[escc] > 0L)) {
    v <- ts$fullvars[[inv[escc[k]]]]
    le <- .st_linear_expansion(-resid[[k]], v)
    solved[[length(solved) + 1L]] <- Equation(v, -le$b / le$a)
  }
  R <- resid[ri]
  exprs <- c(R, lapply(solved, `[[`, "rhs"))
  if (any(vapply(exprs, .has_derivative, TRUE))) return(NULL)
  local <- c(.keys(Z), .keys(lapply(solved, `[[`, "lhs")))
  inputs <- Filter(function(v) !(.key(v) %in% local), .collect_vars(exprs))
  ivk <- .key(ts$iv)
  if (!(ivk %in% .keys(inputs)) && ivk %in% .keys(free_symbols(Reduce(`+`, exprs)))) inputs <- c(inputs, list(ts$iv))
  args <- list(list(name = "z", vars = Z), list(name = "q", vars = inputs))
  F <- .build_closure(R, args, observed = solved, out = "vector")
  gs <- smap_merge(.f(ts$sys, "initial_conditions"), .f(ts$sys, "guesses"))
  z0 <- vapply(Z, function(z) {
    gv <- smap_get(gs, z)
    if (is.numeric(gv) && length(gv) == 1L && is.finite(gv)) gv else 0
  }, 0)
  # ponytail: the warm start belongs to the compiled system and is shared by
  # every problem built from it; a block with several roots can then follow
  # the branch of an earlier solve. Per integrator state needs the integrator
  # passed to generated functions.
  st <- new.env(parent = emptyenv())
  st$z <- z0
  st$q <- NULL; st$valid <- FALSE
  fn <- function(j, ...) {
    q <- as.numeric(c(...))
    if (!st$valid || !identical(q, st$q)) {
      sol <- .newton_solve(function(z) F(z, q), st$z, tol = 1e-13)
      if (!sol$converged || any(!is.finite(sol$x))) {
        sol <- .newton_solve(function(z) F(z, q), z0, tol = 1e-13)
        if (!sol$converged) stop("the inline nonlinear solve of ", paste(vapply(Z, format, ""), collapse = ", "),
                                 " did not converge", call. = FALSE)
      }
      st$q <- q; st$z <- sol$x; st$valid <- TRUE
    }
    st$z[as.integer(j)]
  }
  .mtk$inline_nl_counter <- .mtk$inline_nl_counter + 1L
  head <- register_symbolic(paste0("inline_nlsolve_", .mtk$inline_nl_counter), fn)
  sol <- lapply(seq_along(Z), function(j) do.call(head, c(list(j), inputs)))
  list(eqs = escc[ri], vars = vscc[zi], sol = sol)
}
