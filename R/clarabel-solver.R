# The Clarabel interior point solver for conic programs
#
#   minimize   x'Px/2 + q'x   subject to   Ax + s = b,  s in K
#
# with K a product of zero, nonnegative, second order, exponential, power,
# generalized power and positive semidefinite (triangle) cones: the homogeneous
# embedding, Ruiz equilibration, presolve of infinite bounds, the direct LDL
# KKT solver with static and dynamic regularization and iterative refinement,
# and the termination rules of Clarabel.jl.
#
# Ported from Clarabel.jl 0.11.1 (https://github.com/oxfordcontrol/Clarabel.jl),
# Copyright 2022 University of Oxford Control Group, Apache License 2.0
# (inst/CLARABEL-LICENSE). Floating point operations follow the order of the
# Julia code; BLAS dot products are emulated by .jl_dot.

.cl_rng <- function(a, b) if (b >= a) a:b else integer(0)
.cl_ssum <- function(v) {
  s <- 0
  for (x in v) s <- s + x
  s
}
.cl_norm2 <- function(v) {
  if (!length(v)) return(0)
  m <- max(abs(v))
  if (is.na(m) || m == 0 || is.infinite(m)) return(m)
  if (length(v) >= 32L) return(sqrt(sum(v * v)))
  s <- 0
  for (x in v) s <- s + x * x
  sqrt(s)
}
.cl_norminf <- function(v) if (length(v)) max(abs(v)) else 0
.cl_mean <- function(v) .cl_ssum(v) / length(v)
.cl_clip <- function(s, lo, hi) ifelse(s < lo, lo, ifelse(s > hi, hi, s))
.cl_logsafe <- function(v) if (v < 0) -.Machine$double.xmax else log(v)
.cl_floatmax <- .Machine$double.xmax
.cl_eps <- .Machine$double.eps
.cl_triangular_number <- function(k) (k * (k + 1)) %/% 2

# The scaled two norm of x .* v (norm_scaled) and the shifted norm of x + a dx
.cl_norm_scaled_vec <- function(xi) {
  scale <- 0; sumsq <- 1
  for (v in xi) {
    if (v != 0) {
      absxi <- abs(v)
      if (scale < absxi) {
        r <- scale / absxi
        scale <- absxi
        sumsq <- 1 + sumsq * r * r
      } else {
        r <- absxi / scale
        sumsq <- sumsq + r * r
      }
    }
  }
  scale * sqrt(sumsq)
}
.cl_norm_scaled <- function(x, v) .cl_norm_scaled_vec(x * v)
.cl_norm_shifted <- function(x, dx, a) .cl_norm_scaled_vec(x + a * dx)
.cl_norm_inf_scaled <- function(m, v) if (length(v)) max(0, abs(m * v)) else 0
.cl_dot_shifted <- function(z, s, dz, ds, a) {
  out <- 0
  zi <- z + a * dz; si <- s + a * ds
  for (i in seq_along(z)) out <- out + zi[i] * si[i]
  out
}

# ---------------------------------------------------------------------------
# Sparse matrices (1 based CSC lists made by .csc)

# sparse(I, J, V, m, n) of Julia: entries sorted by column then row, repeated
# entries summed, explicit zeros kept
.cl_sparse <- function(I, J, V, m, n) {
  if (!length(I)) return(.csc(m, n, rep(1L, n + 1L), integer(), numeric()))
  o <- order(J, I)
  I <- I[o]; J <- J[o]; V <- V[o]
  key <- (J - 1) * m + I
  first <- c(TRUE, key[-1L] != key[-length(key)])
  grp <- cumsum(first)
  vals <- numeric(sum(first))
  for (k in seq_along(V)) vals[grp[k]] <- vals[grp[k]] + V[k]
  rows <- I[first]; cols <- J[first]
  colptr <- c(1L, 1L + cumsum(tabulate(cols, n)))
  .csc(m, n, as.integer(colptr), as.integer(rows), vals)
}
.cl_nnz <- function(M) length(M$nzval)
.cl_csc_cols <- function(M) rep(seq_len(M$n), diff(M$colptr))

# rows of a CSC matrix kept by a logical mask (A[keep, :])
.cl_csc_keep_rows <- function(M, keep) {
  newrow <- cumsum(keep)
  cols <- .cl_csc_cols(M)
  sel <- keep[M$rowval]
  .csc(sum(keep), M$n, as.integer(c(1L, 1L + cumsum(tabulate(cols[sel], M$n)))), as.integer(newrow[M$rowval[sel]]),
       M$nzval[sel])
}

# y = alpha * Symmetric(M) x + y, as the sparse symmetric _mul! of SparseArrays
.cl_symv <- function(M, x, y, alpha = 1) {
  rv <- M$rowval; nzv <- M$nzval; cp <- M$colptr
  for (col in seq_len(M$n)) {
    axj <- x[col] * alpha
    sumcol <- 0
    for (j in .cl_rng(cp[col], cp[col + 1L] - 1L)) {
      row <- rv[j]
      if (row > col) next
      a <- nzv[j]
      if (row == col) {
        sumcol <- sumcol + a * x[row]
      } else {
        y[row] <- y[row] + a * axj
        sumcol <- sumcol + a * x[row]
      }
    }
    y[col] <- y[col] + alpha * sumcol
  }
  y
}
# -A'z (mul!(r, A', z, -1, 0))
.cl_neg_Atz <- function(A, z) {
  out <- numeric(A$n)
  for (col in seq_len(A$n)) {
    tmp <- 0
    for (j in .cl_rng(A$colptr[col], A$colptr[col + 1L] - 1L)) tmp <- tmp + A$nzval[j] * z[A$rowval[j]]
    out[col] <- 0 + tmp * -1
  }
  out
}
# y + A x (mul!(y, A, x, 1, 1))
.cl_Ax_add <- function(A, x, y) {
  for (col in seq_len(A$n)) {
    axj <- x[col] * 1
    for (j in .cl_rng(A$colptr[col], A$colptr[col + 1L] - 1L)) {
      r <- A$rowval[j]
      y[r] <- y[r] + A$nzval[j] * axj
    }
  }
  y
}
.cl_quad_form <- function(x, M, y) {
  if (M$n == 0) return(0)
  out <- 0
  for (j in seq_len(M$n)) {
    tmp1 <- 0; tmp2 <- 0
    for (p in .cl_rng(M$colptr[j], M$colptr[j + 1L] - 1L)) {
      i <- M$rowval[p]
      if (i < j) {
        tmp1 <- tmp1 + M$nzval[p] * x[i]
        tmp2 <- tmp2 + M$nzval[p] * y[i]
      } else if (i == j) {
        out <- out + M$nzval[p] * x[i] * y[i]
      }
    }
    out <- out + tmp1 * y[j] + tmp2 * x[j]
  }
  out
}

# ---------------------------------------------------------------------------
# Settings

.cl_default_settings <- function() list(
  max_iter = 200L, time_limit = Inf, verbose = TRUE, max_step_fraction = 0.99,
  tol_gap_abs = 1e-8, tol_gap_rel = 1e-8, tol_feas = 1e-8, tol_infeas_abs = 1e-8, tol_infeas_rel = 1e-8,
  tol_ktratio = 1e-6, reduced_tol_gap_abs = 5e-5, reduced_tol_gap_rel = 5e-5, reduced_tol_feas = 1e-4,
  reduced_tol_infeas_abs = 5e-12, reduced_tol_infeas_rel = 5e-5, reduced_tol_ktratio = 1e-4,
  equilibrate_enable = TRUE, equilibrate_max_iter = 10L, equilibrate_min_scaling = 1e-4,
  equilibrate_max_scaling = 1e+4, linesearch_backtrack_step = 0.8, min_switch_step_length = 1e-1,
  min_terminate_step_length = 1e-4, max_threads = 0L, direct_kkt_solver = TRUE, direct_solve_method = "auto",
  static_regularization_enable = TRUE, static_regularization_constant = 1e-8,
  static_regularization_proportional = .Machine$double.eps^2, dynamic_regularization_enable = TRUE,
  dynamic_regularization_eps = 1e-13, dynamic_regularization_delta = 2e-7, iterative_refinement_enable = TRUE,
  iterative_refinement_reltol = 1e-13, iterative_refinement_abstol = 1e-12, iterative_refinement_max_iter = 10L,
  iterative_refinement_stop_ratio = 5, presolve_enable = TRUE, input_sparse_dropzeros = FALSE,
  chordal_decomposition_enable = TRUE, chordal_decomposition_merge_method = "clique_graph",
  chordal_decomposition_compact = TRUE, chordal_decomposition_complete_dual = TRUE)

.CL_STATUS <- c(UNSOLVED = "unsolved", SOLVED = "solved", PRIMAL_INFEASIBLE = "primal infeasible",
                DUAL_INFEASIBLE = "dual infeasible", ALMOST_SOLVED = "solved (reduced accuracy)",
                ALMOST_PRIMAL_INFEASIBLE = "primal infeasible (reduced accuracy)",
                ALMOST_DUAL_INFEASIBLE = "dual infeasible (reduced accuracy)", MAX_ITERATIONS = "iteration limit",
                MAX_TIME = "time limit", NUMERICAL_ERROR = "numerical error",
                INSUFFICIENT_PROGRESS = "insufficient progress")
.cl_status_is_infeasible <- function(st)
  st %in% c("PRIMAL_INFEASIBLE", "DUAL_INFEASIBLE", "ALMOST_PRIMAL_INFEASIBLE", "ALMOST_DUAL_INFEASIBLE")
.cl_status_is_errored <- function(st) st %in% c("NUMERICAL_ERROR", "INSUFFICIENT_PROGRESS")

# ---------------------------------------------------------------------------
# Cone specifications (SupportedCone): list(type, dim) and friends

.cl_spec <- function(type, dim = NULL, alpha = NULL, dim2 = NULL)
  list(type = type, dim = dim, alpha = alpha, dim2 = dim2)
.cl_nvars <- function(spec) switch(spec$type,
  PSDTriangleConeT = .cl_triangular_number(spec$dim), ExponentialConeT = 3L, PowerConeT = 3L,
  GenPowerConeT = length(spec$alpha) + spec$dim2, spec$dim)

# Consecutive nonnegative cones (and one dimensional second order and PSD
# cones) merge into one nonnegative cone; empty cones are dropped.
.cl_cones_collapsed <- function(cones) {
  out <- list()
  merges <- function(c) c$type == "NonnegativeConeT" || (c$type %in% c("SecondOrderConeT", "PSDTriangleConeT") &&
                                                           c$dim == 1)
  i <- 1L
  while (i <= length(cones)) {
    cone <- cones[[i]]; i <- i + 1L
    if (.cl_nvars(cone) == 0) next
    if (merges(cone)) {
      total <- cone$dim
      while (i <= length(cones)) {
        nx <- cones[[i]]
        if (.cl_nvars(nx) == 0) {
        } else if (nx$type == "NonnegativeConeT") {
          total <- total + nx$dim
        } else if (merges(nx)) {
          total <- total + 1
        } else break
        i <- i + 1L
      }
      out[[length(out) + 1L]] <- .cl_spec("NonnegativeConeT", total)
    } else out[[length(out) + 1L]] <- cone
  }
  out
}

# ---------------------------------------------------------------------------
# Problem data: presolve, equilibration

.cl_problem_data <- function(P, q, A, b, cones, settings) {
  cones <- .cl_cones_collapsed(cones)
  # P is upper triangular here: the MOI wrapper passes an empty or triu P
  infbound <- 1e20
  presolver <- NULL
  if (settings$presolve_enable) {
    keep <- rep(TRUE, length(b))
    lim <- infbound * (1 - 10 * .cl_eps)
    idx <- 1L
    for (cone in cones) {
      nv <- .cl_nvars(cone)
      if (cone$type == "NonnegativeConeT") {
        for (k in .cl_rng(idx, idx + nv - 1L)) if (b[k] > lim) keep[k] <- FALSE
      }
      idx <- idx + nv
    }
    if (!all(keep)) {
      presolver <- list(keep = keep, infbound = infbound, mfull = length(b))
      A <- .cl_csc_keep_rows(A, keep)
      b <- b[keep]
      newcones <- list()
      idx <- 1L
      for (cone in cones) {
        nv <- .cl_nvars(cone)
        if (cone$type == "NonnegativeConeT") {
          nkeep <- sum(keep[.cl_rng(idx, idx + nv - 1L)])
          if (nkeep > 0) newcones[[length(newcones) + 1L]] <- .cl_spec("NonnegativeConeT", nkeep)
        } else newcones[[length(newcones) + 1L]] <- cone
        idx <- idx + nv
      }
      cones <- newcones
    }
  }
  # a PSD cone with a sparse aggregate pattern splits into its cliques
  chordal <- .chd_try_info(A, b, cones, settings)
  if (!is.null(chordal)) {
    r <- .chd_augment(chordal, P, q, A, b, settings)
    chordal <- r$info; P <- r$P; q <- r$q; A <- r$A; b <- r$b; cones <- r$cones
  }
  b <- pmin(b, infbound)
  n <- A$n; m <- A$m
  data <- new.env(parent = emptyenv())
  data$P <- P; data$q <- q; data$A <- A; data$b <- b; data$cones <- cones
  data$n <- n; data$m <- m
  data$d <- rep(1, n); data$dinv <- rep(1, n); data$e <- rep(1, m); data$einv <- rep(1, m); data$c <- 1
  data$normq <- .cl_norminf(q); data$normb <- .cl_norminf(b)
  data$presolver <- presolver
  data$chordal_info <- chordal
  data
}

.cl_col_norms_sym <- function(norms, M) {
  for (i in seq_len(M$n)) for (j in .cl_rng(M$colptr[i], M$colptr[i + 1L] - 1L)) {
    tmp <- abs(M$nzval[j]); r <- M$rowval[j]
    norms[i] <- if (norms[i] > tmp) norms[i] else tmp
    norms[r] <- if (norms[r] > tmp) norms[r] else tmp
  }
  norms
}
.cl_col_norms <- function(norms, M) {
  for (i in seq_along(norms)) for (j in .cl_rng(M$colptr[i], M$colptr[i + 1L] - 1L)) {
    tmp <- abs(M$nzval[j])
    norms[i] <- if (norms[i] > tmp) norms[i] else tmp
  }
  norms
}
.cl_row_norms <- function(norms, M) {
  for (i in seq_along(M$nzval)) {
    idx <- M$rowval[i]; tmp <- abs(M$nzval[i])
    norms[idx] <- if (norms[idx] > tmp) norms[idx] else tmp
  }
  norms
}
.cl_lrscale <- function(L, M, R) {
  M$nzval <- M$nzval * (L[M$rowval] * R[.cl_csc_cols(M)])
  M
}

.cl_equilibrate <- function(data, cones, settings) {
  if (!settings$equilibrate_enable) return(invisible())
  d <- data$d; e <- data$e; cc <- data$c
  P <- data$P; A <- data$A; q <- data$q; b <- data$b
  smin <- settings$equilibrate_min_scaling; smax <- settings$equilibrate_max_scaling
  for (it in seq_len(settings$equilibrate_max_iter)) {
    dwork <- .cl_col_norms_sym(numeric(data$n), P)
    dwork <- .cl_col_norms(dwork, A)
    ework <- .cl_row_norms(numeric(data$m), A)
    dwork[dwork == 0] <- 1
    ework[ework == 0] <- 1
    dwork <- 1 / sqrt(dwork)
    ework <- 1 / sqrt(ework)
    dwork <- .cl_clip(dwork, smin / d, smax / d)
    ework <- .cl_clip(ework, smin / e, smax / e)
    P <- .cl_lrscale(dwork, P, dwork)
    A <- .cl_lrscale(ework, A, dwork)
    q <- q * dwork
    b <- b * ework
    d <- d * dwork
    e <- e * ework
    pn <- .cl_col_norms(numeric(data$n), P)
    mean_col_norm_P <- .cl_mean(pn)
    inf_norm_q <- .cl_norminf(q)
    if (mean_col_norm_P != 0 && inf_norm_q != 0) {
      scale_cost <- max(inf_norm_q, mean_col_norm_P)
      ctmp <- 1 / scale_cost
      ctmp <- .cl_clip(ctmp, smin / cc, smax / cc)
      P$nzval <- P$nzval * ctmp
      q <- q * ctmp
      cc <- cc * ctmp
    }
  }
  # rectify_equilibration: cones other than zero and nonnegative cones need a
  # uniform scaling of their rows
  delta <- rep(1, data$m)
  changed <- FALSE
  for (k in seq_along(cones$cones)) {
    K <- cones$cones[[k]]
    if (!(K$type %in% c("Zero", "Nonnegative"))) {
      rng <- cones$rng_cones[[k]]
      tmp <- .cl_mean(e[rng])
      delta[rng] <- tmp / e[rng]
      changed <- TRUE
    }
  }
  if (changed) {
    A$nzval <- A$nzval * delta[A$rowval]
    b <- b * delta
    e <- e * delta
  }
  data$P <- P; data$A <- A; data$q <- q; data$b <- b
  data$d <- d; data$e <- e; data$c <- cc
  data$dinv <- 1 / d; data$einv <- 1 / e
  invisible()
}

# ---------------------------------------------------------------------------
# The composite cone

.cl_composite <- function(specs) {
  cones <- lapply(specs, .clk_make)
  numel <- vapply(cones, function(K) K$numel, 0)
  ends <- cumsum(numel)
  rng_cones <- lapply(seq_along(cones), function(i) .cl_rng(ends[i] - numel[i] + 1L, ends[i]))
  blen <- vapply(cones, function(K) if (.clk_Hs_is_diagonal(K)) K$numel else .cl_triangular_number(K$numel), 0)
  bends <- cumsum(blen)
  rng_blocks <- lapply(seq_along(cones), function(i) .cl_rng(bends[i] - blen[i] + 1L, bends[i]))
  list(cones = cones, numel = sum(numel), degree = sum(vapply(cones, function(K) K$degree, 0)),
       rng_cones = rng_cones, rng_blocks = rng_blocks,
       is_symmetric = all(vapply(cones, .clk_is_symmetric, TRUE)))
}

.cl_margins <- function(cones, z, pd) {
  a <- .cl_floatmax; bsum <- 0
  for (k in seq_along(cones$cones)) {
    mb <- .clk_margins(cones$cones[[k]], z[cones$rng_cones[[k]]], pd)
    a <- min(a, mb[1L]); bsum <- bsum + mb[2L]
  }
  c(a, bsum)
}
.cl_scaled_unit_shift <- function(cones, z, a, pd) {
  for (k in seq_along(cones$cones)) {
    r <- cones$rng_cones[[k]]
    z[r] <- .clk_scaled_unit_shift(cones$cones[[k]], z[r], a, pd)
  }
  z
}
.cl_get_Hs <- function(cones) {
  out <- numeric(if (length(cones$rng_blocks)) max(0, unlist(cones$rng_blocks)) else 0)
  for (k in seq_along(cones$cones)) out[cones$rng_blocks[[k]]] <- .clk_get_Hs(cones$cones[[k]])
  out
}
.cl_mul_Hs <- function(cones, x) {
  y <- numeric(length(x))
  for (k in seq_along(cones$cones)) {
    r <- cones$rng_cones[[k]]
    y[r] <- .clk_mul_Hs(cones$cones[[k]], x[r])
  }
  y
}
.cl_affine_ds <- function(cones, s) {
  ds <- numeric(length(s))
  for (k in seq_along(cones$cones)) {
    r <- cones$rng_cones[[k]]
    ds[r] <- .clk_affine_ds(cones$cones[[k]], s[r])
  }
  ds
}
.cl_combined_ds_shift <- function(cones, step_z, step_s, sigmamu) {
  shift <- numeric(length(step_z))
  for (k in seq_along(cones$cones)) {
    r <- cones$rng_cones[[k]]
    shift[r] <- .clk_combined_ds_shift(cones$cones[[k]], step_z[r], step_s[r], sigmamu)
  }
  shift
}
.cl_ds_from_dz_offset <- function(cones, ds, z) {
  out <- numeric(length(ds))
  for (k in seq_along(cones$cones)) {
    r <- cones$rng_cones[[k]]
    out[r] <- .clk_ds_from_dz_offset(cones$cones[[k]], ds[r], z[r])
  }
  out
}
.cl_step_length <- function(cones, dz, ds, z, s, settings, amax) {
  a <- amax
  inner <- function(a, symcond) {
    for (k in seq_along(cones$cones)) {
      K <- cones$cones[[k]]
      if (.clk_is_symmetric(K) != symcond) next
      r <- cones$rng_cones[[k]]
      ab <- .clk_step_length(K, dz[r], ds[r], z[r], s[r], settings, a)
      a <- min(a, ab[1L], ab[2L])
    }
    a
  }
  a <- inner(a, TRUE)
  if (!cones$is_symmetric) a <- min(a, 1 - sqrt(.cl_eps))
  inner(a, FALSE)
}
.cl_compute_barrier <- function(cones, z, s, dz, ds, a) {
  barrier <- 0
  for (k in seq_along(cones$cones)) {
    r <- cones$rng_cones[[k]]
    barrier <- barrier + .clk_compute_barrier(cones$cones[[k]], z[r], s[r], dz[r], ds[r], a)
  }
  barrier
}
.cl_update_scaling <- function(cones, s, z, mu, strategy) {
  for (k in seq_along(cones$cones)) {
    r <- cones$rng_cones[[k]]
    if (!.clk_update_scaling(cones$cones[[k]], s[r], z[r], mu, strategy)) return(FALSE)
  }
  TRUE
}

# ---------------------------------------------------------------------------
# The KKT matrix (upper triangle) and the direct LDL solver

.cl_kkt_assemble <- function(P, A, cones) {
  m <- A$m; n <- P$n
  smaps <- list()
  for (K in cones$cones) if (.clk_sparse_expandable(K)) {
    smaps[[length(smaps) + 1L]] <- if (K$type == "SecondOrder")
      list(type = "SOC", u = integer(K$numel), v = integer(K$numel), D = integer(2L), pdim = 2L, Dsigns = c(-1, 1))
    else list(type = "GenPow", p = integer(K$numel), q = integer(K$dim1), r = integer(K$dim2), D = integer(3L),
              pdim = 3L, Dsigns = c(-1, -1, 1))
  }
  pd <- sum(vapply(smaps, function(s) s$pdim, 0L))
  nnz_vec <- sum(vapply(smaps, function(s) if (s$type == "SOC") 2L * length(s$u) else
    length(s$p) + length(s$q) + length(s$r), 0L))
  N <- m + n + pd
  Pcols <- P$colptr
  has_diag <- vapply(seq_len(n), function(i) Pcols[i] != Pcols[i + 1L] && P$rowval[Pcols[i + 1L] - 1L] == i, TRUE)
  Hsmap <- integer(if (length(cones$rng_blocks)) max(0, unlist(cones$rng_blocks)) else 0)
  nnzKKT <- .cl_nnz(P) + n - sum(has_diag) + .cl_nnz(A) + length(Hsmap) + nnz_vec + pd
  colptr <- integer(N + 1L); rowval <- integer(nnzKKT); nzval <- numeric(nnzKKT)
  # column counts
  colptr[seq_len(n)] <- colptr[seq_len(n)] + diff(Pcols)
  colptr[seq_len(n)] <- colptr[seq_len(n)] + !has_diag
  for (r in A$rowval) colptr[r + n] <- colptr[r + n] + 1L
  pcol <- m + n + 1L
  si <- 0L
  for (k in seq_along(cones$cones)) {
    K <- cones$cones[[k]]
    row <- cones$rng_cones[[k]][1L] + n
    bd <- K$numel
    cols <- .cl_rng(row, row + bd - 1L)
    colptr[cols] <- colptr[cols] + if (.clk_Hs_is_diagonal(K)) 1L else seq_len(bd)
    if (.clk_sparse_expandable(K)) {
      si <- si + 1L
      sm <- smaps[[si]]
      if (sm$type == "SOC") {
        colptr[pcol] <- colptr[pcol] + bd
        colptr[pcol + 1L] <- colptr[pcol + 1L] + bd
      } else {
        colptr[pcol] <- colptr[pcol] + K$dim1
        colptr[pcol + 1L] <- colptr[pcol + 1L] + K$dim2
        colptr[pcol + 2L] <- colptr[pcol + 2L] + bd
      }
      dc <- pcol + seq_len(sm$pdim) - 1L
      colptr[dc] <- colptr[dc] + 1L
      pcol <- pcol + sm$pdim
    }
  }
  cur <- 1L
  for (i in seq_len(N + 1L)) {
    cnt <- colptr[i]; colptr[i] <- cur; cur <- cur + cnt
  }
  # fill; colptr[col] is the next free slot of column col
  put <- function(col, row, val) {
    dest <- colptr[col]
    rowval[dest] <<- row; nzval[dest] <<- val
    colptr[col] <<- dest + 1L
    dest
  }
  Pmap <- integer(.cl_nnz(P))
  for (i in seq_len(n)) for (j in .cl_rng(Pcols[i], Pcols[i + 1L] - 1L)) Pmap[j] <- put(i, P$rowval[j], P$nzval[j])
  for (i in seq_len(n)) if (!has_diag[i]) put(i, i, 0)
  Amap <- integer(.cl_nnz(A))
  for (i in seq_len(A$n)) for (j in .cl_rng(A$colptr[i], A$colptr[i + 1L] - 1L))
    Amap[j] <- put(A$rowval[j] + n, i, A$nzval[j])
  pcol <- m + n + 1L
  si <- 0L
  for (k in seq_along(cones$cones)) {
    K <- cones$cones[[k]]
    row <- cones$rng_cones[[k]][1L] + n
    bd <- K$numel
    blk <- cones$rng_blocks[[k]]
    if (.clk_Hs_is_diagonal(K)) {
      for (i in seq_len(bd)) Hsmap[blk[i]] <- put(row + i - 1L, row + i - 1L, 0)
    } else {
      kidx <- 1L
      for (col in .cl_rng(row, row + bd - 1L)) for (r in row:col) {
        Hsmap[blk[kidx]] <- put(col, r, 0)
        kidx <- kidx + 1L
      }
    }
    if (.clk_sparse_expandable(K)) {
      si <- si + 1L
      sm <- smaps[[si]]
      colvec <- function(len, initrow, col) vapply(seq_len(len), function(i) put(col, initrow + i - 1L, 0), 0L)
      if (sm$type == "SOC") {
        sm$v <- colvec(bd, row, pcol)
        sm$u <- colvec(bd, row, pcol + 1L)
      } else {
        sm$q <- colvec(K$dim1, row, pcol)
        sm$r <- colvec(K$dim2, row + K$dim1, pcol + 1L)
        sm$p <- colvec(bd, row, pcol + 2L)
      }
      sm$D <- vapply(seq_len(sm$pdim), function(i) put(pcol + i - 1L, pcol + i - 1L, 0), 0L)
      smaps[[si]] <- sm
      pcol <- pcol + sm$pdim
    }
  }
  for (i in rev(seq_len(N))) colptr[i + 1L] <- colptr[i]
  colptr[1L] <- 1L
  KKT <- .csc(N, N, colptr, rowval, nzval)
  list(KKT = KKT, map = list(P = Pmap, A = Amap, Hsblocks = Hsmap, sparse_maps = smaps,
                             diag_full = colptr[-1L] - 1L, diagP = colptr[2:(n + 1L)] - 1L), p = pd)
}

.cl_kkt_new <- function(data, cones, settings) {
  asm <- .cl_kkt_assemble(data$P, data$A, cones)
  m <- data$m; n <- data$n; p <- asm$p
  Dsigns <- c(rep(1, n), rep(-1, m), unlist(lapply(asm$map$sparse_maps, function(s) s$Dsigns)))
  ks <- new.env(parent = emptyenv())
  ks$m <- m; ks$n <- n; ks$p <- p
  ks$KKT <- asm$KKT; ks$map <- asm$map; ks$Dsigns <- Dsigns
  ks$settings <- settings
  ks$x <- numeric(n + m + p); ks$b <- numeric(n + m + p)
  ks$ldl <- .qdldl(asm$KKT, amd_dense_scale = 1.5, logical = TRUE, Dsigns = Dsigns,
                   regularize_eps = settings$dynamic_regularization_eps,
                   regularize_delta = settings$dynamic_regularization_delta)
  kkt <- new.env(parent = emptyenv())
  kkt$solver <- ks
  kkt$x1 <- numeric(n); kkt$z1 <- numeric(m); kkt$x2 <- numeric(n); kkt$z2 <- numeric(m)
  kkt
}

.cl_ks_update_values <- function(ks, idx, values) {
  ks$KKT$nzval[idx] <- values
  .qdldl_update_values(ks$ldl, idx, values)
}
.cl_ks_scale_values <- function(ks, idx, scale) {
  ks$KKT$nzval[idx] <- ks$KKT$nzval[idx] * scale
  w <- ks$ldl$AtoPAPt[idx]
  ks$ldl$triuA$nzval[w] <- ks$ldl$triuA$nzval[w] * scale
}

.cl_kktsolver_update <- function(ks, cones) {
  map <- ks$map
  Hs <- .cl_get_Hs(cones)
  .cl_ks_update_values(ks, map$Hsblocks, Hs * -1)
  si <- 0L
  for (K in cones$cones) if (.clk_sparse_expandable(K)) {
    si <- si + 1L
    sm <- map$sparse_maps[[si]]
    if (sm$type == "SOC") {
      eta2 <- K$eta^2
      .cl_ks_update_values(ks, sm$u, K$sp_u)
      .cl_ks_update_values(ks, sm$v, K$sp_v)
      .cl_ks_scale_values(ks, sm$u, -eta2)
      .cl_ks_scale_values(ks, sm$v, -eta2)
      .cl_ks_update_values(ks, sm$D, c(-eta2, eta2))
    } else {
      sqrtmu <- sqrt(K$mu)
      .cl_ks_update_values(ks, sm$q, K$q)
      .cl_ks_update_values(ks, sm$r, K$r)
      .cl_ks_update_values(ks, sm$p, K$p)
      .cl_ks_scale_values(ks, sm$q, -sqrtmu)
      .cl_ks_scale_values(ks, sm$r, -sqrtmu)
      .cl_ks_scale_values(ks, sm$p, -sqrtmu)
      .cl_ks_update_values(ks, sm$D, c(-1, -1, 1))
    }
  }
  # static regularization, factorization, then the unregularized diagonal back
  st <- ks$settings
  if (st$static_regularization_enable) {
    diag_kkt <- ks$KKT$nzval[map$diag_full]
    eps <- st$static_regularization_constant + st$static_regularization_proportional * .cl_norminf(diag_kkt)
    shifted <- ifelse(ks$Dsigns == 1, diag_kkt + eps, diag_kkt - eps)
    .cl_ks_update_values(ks, map$diag_full, shifted)
    ks$diagonal_regularizer <- eps
  }
  .qdldl_refactor(ks$ldl)
  ok <- all(is.finite(ks$ldl$Dinv))
  if (st$static_regularization_enable) ks$KKT$nzval[map$diag_full] <- diag_kkt
  ok
}

.cl_kktsolver_setrhs <- function(ks, rhsx, rhsz) ks$b <- c(rhsx, rhsz, numeric(ks$p))

# Solve with iterative refinement; the solution is in ks$x
.cl_kktsolver_solve <- function(ks) {
  b <- ks$b
  x <- .qdldl_solve(ks$ldl, b)
  st <- ks$settings
  ok <- TRUE
  if (st$iterative_refinement_enable) {
    refine_err <- function(xi) .cl_symv(ks$KKT, xi, b, -1)
    normb <- .cl_norminf(b)
    e <- refine_err(x)
    norme <- .cl_norminf(e)
    if (!is.finite(norme)) ok <- FALSE else {
      for (i in seq_len(st$iterative_refinement_max_iter)) {
        if (norme <= st$iterative_refinement_abstol + st$iterative_refinement_reltol * normb) break
        lastnorme <- norme
        dx <- .qdldl_solve(ks$ldl, e)
        dx <- dx + x
        e <- refine_err(dx)
        norme <- .cl_norminf(e)
        if (!is.finite(norme)) {
          ok <- FALSE
          break
        }
        ratio <- lastnorme / norme
        if (ratio < st$iterative_refinement_stop_ratio) {
          if (ratio > 1) x <- dx
          break
        }
        x <- dx
      }
    }
  } else ok <- all(is.finite(x))
  ks$x <- x
  ok
}
.cl_lhs_x <- function(ks) ks$x[seq_len(ks$n)]
.cl_lhs_z <- function(ks) ks$x[ks$n + seq_len(ks$m)]

.cl_kkt_update <- function(kkt, data, cones) {
  ks <- kkt$solver
  if (!.cl_kktsolver_update(ks, cones)) return(FALSE)
  .cl_kktsolver_setrhs(ks, -data$q, data$b)
  ok <- .cl_kktsolver_solve(ks)
  if (ok) {
    kkt$x2 <- .cl_lhs_x(ks); kkt$z2 <- .cl_lhs_z(ks)
  }
  ok
}

.cl_kkt_solve_initial_point <- function(kkt, vars, data) {
  ks <- kkt$solver
  if (.cl_nnz(data$P) == 0) {
    .cl_kktsolver_setrhs(ks, numeric(data$n), data$b)
    ok <- .cl_kktsolver_solve(ks)
    if (ok) {
      vars$x <- .cl_lhs_x(ks); vars$s <- .cl_lhs_z(ks)
    }
    vars$s <- -vars$s
    if (!ok) return(list(ok = ok, vars = vars))
    .cl_kktsolver_setrhs(ks, -data$q, numeric(data$m))
    ok <- .cl_kktsolver_solve(ks)
    if (ok) vars$z <- .cl_lhs_z(ks)
  } else {
    .cl_kktsolver_setrhs(ks, -data$q, data$b)
    ok <- .cl_kktsolver_solve(ks)
    if (ok) {
      vars$x <- .cl_lhs_x(ks); vars$z <- .cl_lhs_z(ks)
    }
    vars$s <- -vars$z
  }
  list(ok = ok, vars = vars)
}

.cl_kkt_solve <- function(kkt, lhs, rhs, data, vars, cones, steptype) {
  ks <- kkt$solver
  workx <- rhs$x
  dsconst <- if (steptype == "affine") vars$s else .cl_ds_from_dz_offset(cones, rhs$s, vars$z)
  workz <- dsconst - rhs$z
  .cl_kktsolver_setrhs(ks, workx, workz)
  if (!.cl_kktsolver_solve(ks)) return(list(ok = FALSE, lhs = lhs))
  x1 <- .cl_lhs_x(ks); z1 <- .cl_lhs_z(ks)
  x2 <- kkt$x2; z2 <- kkt$z2
  xi <- vars$x / vars$tau
  P <- data$P
  tau_num <- rhs$tau - rhs$kappa / vars$tau + .jl_dot(data$q, x1) + .jl_dot(data$b, z1) +
    2 * .cl_quad_form(xi, P, x1)
  xi <- xi - x2
  tau_den <- vars$kappa / vars$tau - .jl_dot(data$q, x2) - .jl_dot(data$b, z2)
  tau_den <- tau_den + (.cl_quad_form(xi, P, xi) - .cl_quad_form(x2, P, x2))
  lhs$tau <- tau_num / tau_den
  lhs$x <- x1 + lhs$tau * x2
  lhs$z <- z1 + lhs$tau * z2
  lhs$s <- -(.cl_mul_Hs(cones, lhs$z) + dsconst)
  lhs$kappa <- -(rhs$kappa + vars$kappa * lhs$tau) / vars$tau
  list(ok = TRUE, lhs = lhs)
}

# ---------------------------------------------------------------------------
# Variables, residuals and information

.cl_vars <- function(n, m) list(x = numeric(n), s = numeric(m), z = numeric(m), tau = 1, kappa = 1)

.cl_residuals_update <- function(vars, data) {
  qx <- .jl_dot(data$q, vars$x)
  bz <- .jl_dot(data$b, vars$z)
  sz <- .jl_dot(vars$s, vars$z)
  Px <- .cl_symv(data$P, vars$x, numeric(data$n), 1)
  xPx <- .jl_dot(vars$x, Px)
  rx_inf <- .cl_neg_Atz(data$A, vars$z)
  rz_inf <- .cl_Ax_add(data$A, vars$x, vars$s)
  list(rx = rx_inf - Px - data$q * vars$tau, rz = rz_inf - data$b * vars$tau,
       rtau = qx + bz + vars$kappa + xPx / vars$tau, rx_inf = rx_inf, rz_inf = rz_inf,
       dot_qx = qx, dot_bz = bz, dot_sz = sz, dot_xPx = xPx, Px = Px)
}

.cl_info_update <- function(info, data, vars, res) {
  tauinv <- 1 / vars$tau
  normb <- data$normb; normq <- data$normq
  d <- data$d; dinv <- data$dinv; e <- data$e; einv <- data$einv
  cinv <- 1 / data$c
  xPx2 <- res$dot_xPx * tauinv * tauinv / 2
  info$cost_primal <- (+res$dot_qx * tauinv + xPx2) * cinv
  info$cost_dual <- (-res$dot_bz * tauinv - xPx2) * cinv
  normx <- .cl_norm_scaled(d, vars$x)
  normz <- .cl_norm_scaled(e, vars$z) * cinv
  norms <- .cl_norm_scaled(einv, vars$s)
  info$res_primal_inf <- (.cl_norm_scaled(dinv, res$rx_inf) * cinv) / max(1, normz)
  info$res_dual_inf <- max(.cl_norm_scaled(dinv, res$Px) / max(1, normx),
                           .cl_norm_scaled(einv, res$rz_inf) / max(1, normx + norms))
  normx <- normx * tauinv; normz <- normz * tauinv; norms <- norms * tauinv
  info$res_primal <- .cl_norm_scaled(einv, res$rz) * tauinv / max(1, normb + normx + norms)
  info$res_dual <- .cl_norm_scaled(dinv, res$rx) * tauinv * cinv / max(1, normq + normx + normz)
  info$gap_abs <- abs(info$cost_primal - info$cost_dual)
  info$gap_rel <- info$gap_abs / max(1, min(abs(info$cost_primal), abs(info$cost_dual)))
  info$ktratio <- vars$kappa * tauinv
}

.cl_check_convergence <- function(info, res, st, reduced) {
  pre <- if (reduced) "reduced_" else ""
  g <- function(nm) st[[paste0(pre, nm)]]
  solved <- if (reduced) "ALMOST_SOLVED" else "SOLVED"
  pinf <- if (reduced) "ALMOST_PRIMAL_INFEASIBLE" else "PRIMAL_INFEASIBLE"
  dinf <- if (reduced) "ALMOST_DUAL_INFEASIBLE" else "DUAL_INFEASIBLE"
  if (info$ktratio <= 1 && (info$gap_abs < g("tol_gap_abs") || info$gap_rel < g("tol_gap_rel")) &&
      info$res_primal < g("tol_feas") && info$res_dual < g("tol_feas")) {
    info$status <- solved
  } else if (info$ktratio > 1000 / g("tol_ktratio")) {
    if (res$dot_bz < -g("tol_infeas_abs") && info$res_primal_inf < -g("tol_infeas_rel") * res$dot_bz) {
      info$status <- pinf
    } else if (res$dot_qx < -g("tol_infeas_abs") && info$res_dual_inf < -g("tol_infeas_rel") * res$dot_qx) {
      info$status <- dinf
    }
  }
}

.cl_info_check_termination <- function(info, res, st, iter) {
  info$status <- "UNSOLVED"
  .cl_check_convergence(info, res, st, FALSE)
  if (info$status == "UNSOLVED" && iter > 1 &&
      (info$res_dual > info$prev_res_dual || info$res_primal > info$prev_res_primal)) {
    if (info$ktratio < 100 * .cl_eps && (info$prev_gap_abs < st$tol_gap_abs || info$prev_gap_rel < st$tol_gap_rel))
      info$status <- "INSUFFICIENT_PROGRESS"
    if (info$ktratio < 1) {
      if ((info$res_dual > 100 * st$tol_feas && info$res_dual > 100 * info$prev_res_dual) ||
          (info$res_primal > 100 * st$tol_feas && info$res_primal > 100 * info$prev_res_primal))
        info$status <- "INSUFFICIENT_PROGRESS"
    }
  }
  if (info$status == "UNSOLVED") {
    if (st$max_iter == info$iterations) info$status <- "MAX_ITERATIONS"
    else if (info$solve_time > st$time_limit) info$status <- "MAX_TIME"
  }
  info$status != "UNSOLVED"
}

.cl_info_new <- function() {
  info <- new.env(parent = emptyenv())
  for (nm in c("mu", "sigma", "step_length", "iterations", "cost_primal", "cost_dual", "res_primal", "res_dual",
               "res_primal_inf", "res_dual_inf", "gap_abs", "gap_rel", "ktratio", "prev_cost_primal"))
    assign(nm, 0, envir = info)
  # DefaultInfo: the five fields after prev_cost_primal, and solve_time, start at floatmax
  for (nm in c("prev_cost_dual", "prev_res_primal", "prev_res_dual", "prev_gap_abs", "prev_gap_rel", "solve_time"))
    assign(nm, .cl_floatmax, envir = info)
  info$status <- "UNSOLVED"
  info
}
.cl_info_save_prev <- function(info) {
  info$prev_cost_primal <- info$cost_primal; info$prev_cost_dual <- info$cost_dual
  info$prev_res_primal <- info$res_primal; info$prev_res_dual <- info$res_dual
  info$prev_gap_abs <- info$gap_abs; info$prev_gap_rel <- info$gap_rel
}
.cl_info_reset_to_prev <- function(info) {
  info$cost_primal <- info$prev_cost_primal; info$cost_dual <- info$prev_cost_dual
  info$res_primal <- info$prev_res_primal; info$res_dual <- info$prev_res_dual
  info$gap_abs <- info$prev_gap_abs; info$gap_rel <- info$prev_gap_rel
}

# ---------------------------------------------------------------------------
# Solver

# Solver(P, q, A, b, cones, settings): P and A are CSC lists, cones a list of
# cone specifications. Returns the solver environment.
.cl_setup <- function(P, q, A, b, cones, settings) {
  n <- length(q); m <- length(b)
  pv <- sum(vapply(cones, .cl_nvars, 0))
  if (m != A$m) stop("A and b incompatible dimensions.", call. = FALSE)
  if (pv != m) stop("Constraint dimensions inconsistent with size of cones.", call. = FALSE)
  if (n != A$n) stop("A and q incompatible dimensions.", call. = FALSE)
  if (n != P$m || P$m != P$n) stop("P and q incompatible dimensions.", call. = FALSE)
  S <- new.env(parent = emptyenv())
  S$settings <- settings
  S$info <- .cl_info_new()
  S$solution <- list(x = numeric(A$n), z = numeric(A$m), s = numeric(A$m), status = "UNSOLVED", obj_val = NaN,
                     obj_val_dual = NaN, iterations = 0L, r_prim = NaN, r_dual = NaN)
  S$data <- .cl_problem_data(P, q, A, b, cones, settings)
  S$cones <- .cl_composite(S$data$cones)
  if (S$data$m != S$cones$numel) stop("DimensionMismatch", call. = FALSE)
  S$variables <- .cl_vars(S$data$n, S$data$m)
  .cl_equilibrate(S$data, S$cones, settings)
  S$kkt <- .cl_kkt_new(S$data, S$cones, settings)
  S$step_rhs <- .cl_vars(S$data$n, S$data$m)
  S$step_lhs <- .cl_vars(S$data$n, S$data$m)
  S$prev_vars <- .cl_vars(S$data$n, S$data$m)
  S
}

.cl_default_start <- function(S) {
  cones <- S$cones
  if (cones$is_symmetric) {
    for (K in cones$cones) .clk_set_identity_scaling(K)
    .cl_kkt_update(S$kkt, S$data, cones)
    r <- .cl_kkt_solve_initial_point(S$kkt, S$variables, S$data)
    v <- r$vars
    shift <- function(z, pd) {
      mm <- .cl_margins(cones, z, pd)
      min_margin <- mm[1L]; pos_margin <- mm[2L]
      target <- max(1, 0.1 * pos_margin / cones$degree)
      if (min_margin <= 0) {
        z <- .cl_scaled_unit_shift(cones, z, -min_margin, pd)
        .cl_scaled_unit_shift(cones, z, target, pd)
      } else if (min_margin < target) {
        .cl_scaled_unit_shift(cones, z, target - min_margin, pd)
      } else .cl_scaled_unit_shift(cones, z, 0, pd)
    }
    v$s <- shift(v$s, "primal")
    v$z <- shift(v$z, "dual")
    v$tau <- 1; v$kappa <- 1
    S$variables <- v
  } else {
    v <- S$variables
    for (k in seq_along(cones$cones)) {
      r <- cones$rng_cones[[k]]
      zs <- .clk_unit_initialization(cones$cones[[k]])
      v$z[r] <- zs$z; v$s[r] <- zs$s
    }
    v$x[] <- 0; v$tau <- 1; v$kappa <- 1
    S$variables <- v
  }
}

.cl_calc_step_length <- function(S, steptype) {
  v <- S$variables; st <- S$step_lhs
  atau <- if (st$tau < 0) -v$tau / st$tau else .cl_floatmax
  akap <- if (st$kappa < 0) -v$kappa / st$kappa else .cl_floatmax
  a <- min(atau, akap, 1)
  a <- .cl_step_length(S$cones, st$z, st$s, v$z, v$s, S$settings, a)
  if (steptype == "combined") a <- a * S$settings$max_step_fraction
  a
}

.cl_barrier <- function(S, a) {
  v <- S$variables; st <- S$step_lhs
  central <- S$cones$degree + 1
  cur_tau <- v$tau + a * st$tau
  cur_kappa <- v$kappa + a * st$kappa
  sz <- .cl_dot_shifted(v$z, v$s, st$z, st$s, a)
  mu <- (sz + cur_tau * cur_kappa) / central
  barrier <- central * .cl_logsafe(mu) - .cl_logsafe(cur_tau) - .cl_logsafe(cur_kappa)
  barrier + .cl_compute_barrier(S$cones, v$z, v$s, st$z, st$s, a)
}

.cl_get_step_length <- function(S, steptype, scaling) {
  a <- .cl_calc_step_length(S, steptype)
  if (!S$cones$is_symmetric && steptype == "combined" && scaling == "Dual") {
    step <- S$settings$linesearch_backtrack_step
    for (j in 1:50) {
      if (.cl_barrier(S, a) < 1) return(a)
      a <- step * a
    }
  }
  a
}

.cl_solve <- function(S) {
  iter <- 0L; sigma <- 1; a <- 0; mu <- .cl_floatmax
  info <- S$info
  info$status <- "UNSOLVED"; info$iterations <- 0L; info$solve_time <- 0
  t0 <- proc.time()[["elapsed"]]
  .cl_default_start(S)
  cones <- S$cones; data <- S$data; st <- S$settings
  scaling <- if (all(vapply(cones$cones, .clk_allows_primal_dual_scaling, TRUE))) "PrimalDual" else "Dual"
  repeat {
    res <- .cl_residuals_update(S$variables, data)
    S$residuals <- res
    mu <- (res$dot_sz + S$variables$tau * S$variables$kappa) / (cones$degree + 1)
    info$mu <- mu; info$step_length <- a; info$sigma <- sigma; info$iterations <- iter
    .cl_info_update(info, data, S$variables, res)
    info$solve_time <- proc.time()[["elapsed"]] - t0
    if (.cl_info_check_termination(info, res, st, iter)) {
      if (info$status != "INSUFFICIENT_PROGRESS") break
      .cl_info_reset_to_prev(info)
      S$variables <- S$prev_vars
      if (!cones$is_symmetric && scaling == "PrimalDual") {
        info$status <- "UNSOLVED"
        scaling <- "Dual"
        next
      }
      break
    }
    if (!.cl_update_scaling(cones, S$variables$s, S$variables$z, mu, scaling)) {
      info$status <- "NUMERICAL_ERROR"
      break
    }
    iter <- iter + 1L
    ok <- .cl_kkt_update(S$kkt, data, cones)
    v <- S$variables
    rhs <- list(x = res$rx, z = res$rz, s = .cl_affine_ds(cones, v$s), tau = res$rtau, kappa = v$tau * v$kappa)
    S$step_rhs <- rhs
    if (ok) {
      r <- .cl_kkt_solve(S$kkt, S$step_lhs, rhs, data, v, cones, "affine")
      ok <- r$ok
      S$step_lhs <- r$lhs
    }
    if (ok) {
      a <- .cl_get_step_length(S, "affine", scaling)
      sigma <- (1 - a)^3
      mstep <- if (iter > 1L) 1 else a
      lhs <- S$step_lhs
      sigmamu <- sigma * mu
      rhs$x <- (1 - sigma) * res$rx
      rhs$tau <- (1 - sigma) * res$rtau
      rhs$kappa <- -sigmamu + mstep * lhs$tau * lhs$kappa + v$tau * v$kappa
      stz <- if (mstep != 1) lhs$z * mstep else lhs$z
      rhs$s <- rhs$s + .cl_combined_ds_shift(cones, stz, lhs$s, sigmamu)
      rhs$z <- (1 - sigma) * res$rz
      S$step_rhs <- rhs
      r <- .cl_kkt_solve(S$kkt, S$step_lhs, rhs, data, v, cones, "combined")
      ok <- r$ok
      S$step_lhs <- r$lhs
    }
    if (!ok) {
      if (!cones$is_symmetric && scaling == "PrimalDual") {
        scaling <- "Dual"
        a <- 0
        next
      }
      info$status <- "NUMERICAL_ERROR"
      a <- 0
      break
    }
    a <- .cl_get_step_length(S, "combined", scaling)
    if (!cones$is_symmetric && scaling == "PrimalDual" && a < st$min_switch_step_length) {
      scaling <- "Dual"
      a <- 0
      next
    } else if (a <= max(0, st$min_terminate_step_length)) {
      info$status <- "INSUFFICIENT_PROGRESS"
      a <- 0
      break
    }
    .cl_info_save_prev(info)
    S$prev_vars <- S$variables
    lhs <- S$step_lhs
    v <- S$variables
    v$x <- v$x + a * lhs$x; v$s <- v$s + a * lhs$s; v$z <- v$z + a * lhs$z
    v$tau <- v$tau + a * lhs$tau; v$kappa <- v$kappa + a * lhs$kappa
    S$variables <- v
  }
  if (a == 0) {
    info$mu <- mu; info$step_length <- a; info$sigma <- sigma; info$iterations <- iter
  }
  if (.cl_status_is_errored(info$status) || info$status %in% c("MAX_ITERATIONS", "MAX_TIME"))
    .cl_check_convergence(info, S$residuals, st, TRUE)
  .cl_solution_post_process(S)
  info$solve_time <- proc.time()[["elapsed"]] - t0
  S$solution$solve_time <- info$solve_time
  S$solution
}

.cl_solution_post_process <- function(S) {
  info <- S$info; data <- S$data; v <- S$variables
  sol <- S$solution
  sol$status <- info$status
  infeas <- .cl_status_is_infeasible(info$status)
  if (infeas) {
    sol$obj_val <- NaN; sol$obj_val_dual <- NaN
  } else {
    sol$obj_val <- info$cost_primal; sol$obj_val_dual <- info$cost_dual
  }
  sol$iterations <- info$iterations
  sol$r_prim <- info$res_primal; sol$r_dual <- info$res_dual
  scaleinv <- if (infeas) 1 / v$kappa else 1 / v$tau
  cinv <- 1 / data$c
  v$x <- v$x * (data$d * scaleinv)
  v$z <- v$z * (data$e * (scaleinv * cinv))
  v$s <- v$s * (data$einv * scaleinv)
  v$tau <- v$tau * scaleinv; v$kappa <- v$kappa * scaleinv
  S$variables <- v
  # the chordal decomposition unwinds before the presolve
  if (!is.null(data$chordal_info)) v <- .chd_reverse(data$chordal_info, v$x, v$s, v$z, data$cones, S$settings)
  if (!is.null(data$presolver)) {
    keep <- data$presolver$keep
    sol$x <- v$x
    sol$s <- rep(data$presolver$infbound, length(keep)); sol$z <- numeric(length(keep))
    sol$s[keep] <- v$s; sol$z[keep] <- v$z
  } else {
    sol$x <- v$x; sol$z <- v$z; sol$s <- v$s
  }
  S$solution <- sol
}
