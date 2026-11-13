# Chordal decomposition of the positive semidefinite cones of Clarabel: the
# aggregate sparsity pattern of each PSD cone, its chordal extension by a
# logical LDL factorization, the supernodal elimination tree, the merge of
# cliques (clique graph, parent child or none), the standard or compact
# augmentation of the problem data, and the reverse map of the solution
# with the optional positive semidefinite completion of the dual.
#
# Ported from Clarabel.jl 0.11.1 (https://github.com/oxfordcontrol/Clarabel.jl),
# src/chordal, Copyright 2022 University of Oxford Control Group, Apache
# License 2.0 (inst/CLARABEL-LICENSE).

.CHD_NO_PARENT <- .Machine$integer.max
.CHD_INACTIVE <- .Machine$integer.max - 1L

# OrderedSet operations on integer vectors
.chd_push <- function(s, v) if (v %in% s) s else c(s, v)
.chd_union <- function(a, b) c(a, b[!(b %in% a)])

# Base.searchsortedfirst and searchsortedlast on v[lo:hi], also on data
# that is not sorted (the result then follows the same bisection)
.chd_ssf <- function(v, x, lo = 1L, hi = length(v)) {
  lo <- lo - 1L; hi <- hi + 1L
  while (lo < hi - 1L) {
    m <- lo + (hi - lo) %/% 2L
    if (v[m] < x) lo <- m else hi <- m
  }
  hi
}
.chd_ssl <- function(v, x, lo = 1L, hi = length(v)) {
  lo <- lo - 1L; hi <- hi + 1L
  while (lo < hi - 1L) {
    m <- lo + (hi - lo) %/% 2L
    if (x < v[m]) hi <- m else lo <- m
  }
  lo
}
# Base.insorted: true when the bisection of searchsorted meets x
.chd_insorted <- function(x, v) {
  lo <- 0L; hi <- length(v) + 1L
  while (lo < hi - 1L) {
    m <- lo + (hi - lo) %/% 2L
    if (v[m] < x) lo <- m else if (x < v[m]) hi <- m else return(TRUE)
  }
  FALSE
}

.chd_trinum <- function(k) (k * (k + 1L)) %/% 2L
.chd_coord_to_upper <- function(i, j) if (i <= j) .chd_trinum(j - 1L) + i else .chd_trinum(i - 1L) + j
.chd_upper_to_coord <- function(k) {
  # isqrt(8k) exactly
  r <- as.integer(floor(sqrt(8 * k)))
  while (r * r > 8 * k) r <- r - 1L
  while ((r + 1) * (r + 1) <= 8 * k) r <- r + 1L
  col <- (r + 1L) %/% 2L
  c(k - .chd_trinum(col - 1L), col)
}

.chd_rng_cones <- function(cones) {
  out <- vector("list", length(cones))
  start <- 1L
  for (k in seq_along(cones)) {
    nv <- .cl_nvars(cones[[k]])
    out[[k]] <- .cl_rng(start, start + nv - 1L)
    start <- start + nv
  }
  out
}

# ---------------------------------------------------------------------------
# Chordal info

.chd_try_info <- function(A, b, cones, settings) {
  if (!isTRUE(settings$chordal_decomposition_enable)) return(NULL)
  if (!any(vapply(cones, function(c) c$type == "PSDTriangleConeT" && c$dim > 2, TRUE))) return(NULL)
  info <- list(init_dims = c(A$n, A$m), init_cones = list(), spatterns = list(), H = NULL, cone_maps = NULL)
  rng <- .chd_rng_cones(cones)
  # rows touched by A or with a nonzero b
  active <- logical(length(b))
  active[A$rowval] <- TRUE
  active <- active | (b != 0)
  for (k in seq_along(cones)) {
    cone <- cones[[k]]
    if (cone$type != "PSDTriangleConeT") next
    sp <- .chd_analyse_psd(active[rng[[k]]], cone$dim, k, settings$chordal_decomposition_merge_method)
    if (!is.null(sp)) info$spatterns[[length(info$spatterns) + 1L]] <- sp
  }
  if (!length(info$spatterns)) return(NULL)
  info$init_cones <- cones
  info
}

.chd_analyse_psd <- function(mask, conedim, coneidx, merge_method) {
  for (i in seq_len(conedim)) mask[.chd_trinum(i)] <- TRUE
  if (all(mask)) return(NULL)
  g <- .chd_find_graph(mask)
  sp <- .chd_sparsity_pattern(g$L, g$ordering, coneidx, merge_method)
  if (sp$sntree$n_cliques == 1L) return(NULL)
  sp
}

.chd_find_graph <- function(mask) {
  n <- .chd_upper_to_coord(length(mask))[2L]
  idx <- which(mask)
  rc <- vapply(idx, .chd_upper_to_coord, integer(2))
  pattern <- .cl_sparse(rc[1L, ], rc[2L, ], rep(1, length(idx)), n, n)
  ws <- .qdldl(pattern, logical = TRUE)
  L <- list(n = n, colptr = ws$Lp, rowval = ws$Li)
  # connect_graph!: an empty column j takes the entry (j + 1, j)
  for (j in seq_len(n - 1L)) {
    if (L$colptr[j + 1L] > L$colptr[j]) next
    pos <- L$colptr[j]
    L$rowval <- append(L$rowval, j + 1L, after = pos - 1L)
    L$colptr[(j + 1L):(n + 1L)] <- L$colptr[(j + 1L):(n + 1L)] + 1L
  }
  list(L = L, ordering = ws$perm)
}

.chd_sparsity_pattern <- function(L, ordering, orig_index, merge_method) {
  t <- .chd_supernode_tree(L)
  if (t$n_cliques > 1L) {
    switch(merge_method,
      none = NULL,
      parent_child = .chd_merge_parent_child(t),
      clique_graph = .chd_merge_clique_graph(t),
      stop("Unknown merge strategy: ", merge_method, call. = FALSE))
  }
  ordering <- .chd_reorder_snode_consecutively(t, ordering)
  t$nblk <- vapply(seq_len(t$n_cliques), function(i) {
    cq <- t$snode_post[i]
    length(t$separators[[cq]]) + length(t$snode[[cq]])
  }, 0L)
  list(sntree = t, ordering = ordering, orig_index = orig_index)
}

# ---------------------------------------------------------------------------
# Supernode tree

.chd_children_from_parent <- function(parent) {
  ch <- rep(list(integer()), length(parent))
  for (i in seq_along(parent)) if (parent[i] != .CHD_NO_PARENT) ch[[parent[i]]] <- .chd_push(ch[[parent[i]]], i)
  ch
}

# post_order!: returns the new post order and the children with their sets sorted
.chd_post_order <- function(parent, children, nc) {
  np <- length(parent)
  rank <- rep(nc + 1L, np)
  root <- which(parent == .CHD_NO_PARENT)[1L]
  stack <- root
  i <- nc
  while (length(stack)) {
    v <- stack[length(stack)]; stack <- stack[-length(stack)]
    rank[v] <- i
    i <- i - 1L
    children[[v]] <- sort(children[[v]])
    stack <- c(stack, children[[v]])
  }
  post <- order(rank)
  list(post = post[seq_len(min(nc, np))], children = children)
}

.chd_supernode_tree <- function(L) {
  n <- L$n
  parent <- vapply(seq_len(n), function(v) if (v == n) .CHD_NO_PARENT else L$rowval[L$colptr[v]], 0L)
  po <- .chd_post_order(parent, .chd_children_from_parent(parent), n)
  post <- po$post
  degree <- integer(n)
  for (v in seq_len(n - 1L)) degree[v] <- L$colptr[v + 1L] - L$colptr[v]
  ps <- .chd_pothen_sun(parent, post, degree)
  snode <- rep(list(integer()), n)
  for (i in seq_len(n)) {
    f <- ps$snode_index[i]
    if (f < 0) snode[[i]] <- .chd_push(snode[[i]], i) else snode[[f]] <- .chd_push(snode[[f]], i)
  }
  snode <- Filter(length, snode)
  t <- new.env(parent = emptyenv())
  t$snode <- snode
  t$snode_parent <- ps$snode_parent
  sp <- .chd_post_order(t$snode_parent, .chd_children_from_parent(t$snode_parent), length(t$snode_parent))
  t$snode_post <- sp$post
  t$snode_children <- sp$children
  t$post <- post
  t$separators <- lapply(snode, function(sn) {
    vrep <- min(sn)
    adj <- L$rowval[.cl_rng(L$colptr[vrep], L$colptr[vrep + 1L] - 1L)]
    sep <- integer()
    for (nb in adj) if (!(nb %in% sn)) sep <- .chd_push(sep, nb)
    sep
  })
  t$nblk <- NULL
  t$n_cliques <- length(snode)
  t
}

.chd_pothen_sun <- function(parent, post, degree) {
  n <- length(parent)
  snode_index <- rep(-1L, n)
  snode_parent <- rep(.CHD_NO_PARENT, n)
  children <- rep(list(integer()), n)
  root_index <- which(parent == .CHD_NO_PARENT)[1L]
  for (v in post) {
    pv <- parent[v]
    if (pv == .CHD_NO_PARENT) children[[root_index]] <- .chd_push(children[[root_index]], v)
    else children[[pv]] <- .chd_push(children[[pv]], v)
    if (pv != .CHD_NO_PARENT) {
      if (degree[v] - 1L == degree[pv] && snode_index[pv] == -1L) {
        if (snode_index[v] < 0) {
          snode_index[pv] <- v
          snode_index[v] <- snode_index[v] - 1L
        } else {
          snode_index[pv] <- snode_index[v]
          snode_index[snode_index[v]] <- snode_index[snode_index[v]] - 1L
        }
      } else {
        if (snode_index[v] < 0) snode_parent[v] <- v
        else snode_parent[snode_index[v]] <- snode_index[v]
      }
    }
    k <- if (snode_index[v] < 0) v else snode_index[v]
    for (w in children[[v]]) {
      l <- if (snode_index[w] < 0) w else snode_index[w]
      if (l != k) snode_parent[l] <- k
    }
  }
  repr <- which(snode_index < 0)
  rp <- snode_parent[repr]
  out <- vapply(rp, function(r) { i <- match(r, repr); if (is.na(i)) .CHD_NO_PARENT else i }, 0L)
  list(snode_parent = out, snode_index = snode_index)
}

.chd_reorder_snode_consecutively <- function(t, ordering) {
  p <- integer(length(t$post))
  k <- 1L
  for (i in t$snode_post) {
    sn <- t$snode[[i]]
    m <- length(sn)
    p[.cl_rng(k, k + m - 1L)] <- sort(sn)
    t$snode[[i]] <- .cl_rng(k, k + m - 1L)
    k <- k + m
  }
  pinv <- integer(length(p)); pinv[p] <- seq_along(p)
  for (s in seq_along(t$separators)) t$separators[[s]] <- pinv[t$separators[[s]]]
  # invpermute!(ordering, pinv)
  out <- integer(length(ordering))
  out[pinv] <- ordering
  out
}

# ---------------------------------------------------------------------------
# Merge strategies

.chd_merge_parent_child <- function(t, t_fill = 8L, t_size = 8L) {
  ci <- length(t$snode) - 1L
  stop_ <- FALSE
  while (!stop_) {
    cq <- t$snode_post[ci]
    cand <- c(t$snode_parent[cq], cq)
    ok <- FALSE
    if (!stop_) {
      par <- cand[1L]; ch <- cand[2L]
      dps <- length(t$snode[[par]]); dpsep <- length(t$separators[[par]])
      dcs <- length(t$snode[[ch]]); dcsep <- length(t$separators[[ch]])
      fill <- ((dps + dpsep) - dcsep) * ((dcs + dcsep) - dcsep)
      ok <- fill <= t_fill || max(dcs, dps) <= t_size
    }
    if (ok) {
      pc <- if (cand[2L] %in% t$snode_children[[cand[1L]]]) cand else rev(cand)
      p <- pc[1L]; ch <- pc[2L]
      t$snode[[p]] <- .chd_union(t$snode[[p]], t$snode[[ch]])
      t$snode[[ch]] <- integer()
      t$separators[[ch]] <- integer()
      for (g in t$snode_children[[ch]]) t$snode_parent[g] <- p
      t$snode_parent[ch] <- .CHD_INACTIVE
      t$snode_children[[p]] <- .chd_union(setdiff(t$snode_children[[p]], ch), t$snode_children[[ch]])
      t$snode_children[[ch]] <- integer()
      t$n_cliques <- t$n_cliques - 1L
    }
    if (ci == 1L) stop_ <- TRUE else ci <- ci - 1L
    if (t$n_cliques == 1L) break
  }
  po <- .chd_post_order(t$snode_parent, t$snode_children, t$n_cliques)
  t$snode_post <- po$post; t$snode_children <- po$children
  invisible()
}

# The reduced clique graph as a dense weight matrix W and a mask S of the
# stored entries of its strictly lower triangular sparse matrix
.chd_merge_clique_graph <- function(t) {
  nsn <- length(t$snode)
  for (i in seq_len(nsn)) t$snode[[i]] <- .chd_union(t$snode[[i]], t$separators[[i]])
  t$snode_parent <- rep(.CHD_INACTIVE, length(t$snode_parent))
  t$snode_children <- rep(list(integer()), length(t$snode_children))
  # compute_reduced_clique_graph! sorts the separators by length, longest first
  t$separators <- t$separators[order(-vapply(t$separators, length, 0L))]
  rows <- integer(); cols <- integer()
  for (sep in t$separators) {
    ci <- which(vapply(t$snode, function(x) all(sep %in% x), TRUE))
    comp <- .chd_components(ci, sep, t$snode)
    nc <- length(ci)
    for (i in seq_len(nc)) for (j in .cl_rng(i + 1L, nc)) {
      a <- ci[i]; b <- ci[j]
      if (comp[match(a, ci)] != comp[match(b, ci)]) {
        rows <- c(rows, max(a, b)); cols <- c(cols, min(a, b))
      }
    }
  }
  N <- t$n_cliques
  W <- matrix(0L, N, N); S <- matrix(FALSE, N, N)
  for (k in seq_along(rows)) {
    w <- .chd_edge_metric(t$snode[[rows[k]]], t$snode[[cols[k]]])
    W[rows[k], cols[k]] <- W[rows[k], cols[k]] + w
    S[rows[k], cols[k]] <- TRUE
  }
  adj <- rep(list(integer()), N)
  for (cc in seq_len(N)) for (r in which(S[, cc])) {
    adj[[r]] <- .chd_push(adj[[r]], cc); adj[[cc]] <- .chd_push(adj[[cc]], r)
  }
  alive <- rep(TRUE, N)
  stop_ <- FALSE
  while (!stop_) {
    # stored entries in column major order, the nzval of the sparse matrix
    nz <- which(S)
    if (!length(nz)) stop("DomainError: Sparse matrix A doesn't contain any entries", call. = FALSE)
    vals <- W[nz]
    perm <- order(-vals)
    cand <- NULL
    for (k in perm) {
      e <- c((nz[k] - 1L) %% N + 1L, (nz[k] - 1L) %/% N + 1L)
      if (.chd_permissible(e, adj, t$snode)) { cand <- e; break }
    }
    if (is.null(cand)) stop("MethodError: no permissible merge candidate", call. = FALSE)
    c1 <- cand[1L]; c2 <- cand[2L]
    if (!(W[c1, c2] >= 0)) { stop_ <- TRUE; next }
    t$snode[[c1]] <- .chd_union(t$snode[[c1]], t$snode[[c2]])
    t$snode[[c2]] <- integer()
    t$n_cliques <- t$n_cliques - 1L
    # update_strategy!
    nbrs <- adj[[c1]]
    newn <- setdiff(setdiff(adj[[c2]], nbrs), c1)
    set_edge <- function(r, cc, v) {
      if (S[r, cc]) W[r, cc] <<- v else if (v != 0) { S[r, cc] <<- TRUE; W[r, cc] <<- v }
    }
    for (ni in nbrs) if (ni != c2)
      set_edge(max(c1, ni), min(c1, ni), .chd_edge_metric(t$snode[[c1]], t$snode[[ni]]))
    for (ni in newn) set_edge(max(c1, ni), min(c1, ni), .chd_edge_metric(t$snode[[c1]], t$snode[[ni]]))
    W[.cl_rng(c2 + 1L, N), c2] <- 0L
    W[c2, seq_len(c2)] <- 0L
    S[S & W == 0L] <- FALSE
    adj[[c1]] <- .chd_union(adj[[c1]], newn)
    for (ni in newn) adj[[ni]] <- .chd_push(adj[[ni]], c1)
    adj[[c2]] <- integer(); alive[c2] <- FALSE
    for (i in seq_len(N)) adj[[i]] <- setdiff(adj[[i]], c2)
    if (t$n_cliques == 1L) break
  }
  # post_process_merge!
  t$snode_post <- which(vapply(t$snode, length, 0L) > 0L)
  t$snode_parent <- rep(.CHD_INACTIVE, length(t$snode))
  if (t$n_cliques > 1L) .chd_clique_tree(t, W, S)
  t$snode <- lapply(t$snode, sort)
  t$separators <- lapply(t$separators, sort)
  invisible()
}

.chd_edge_metric <- function(a, b) {
  n1 <- length(a); n2 <- length(b)
  nm <- n1 + n2 - sum(a %in% b)
  as.integer(n1^3 + n2^3 - nm^3)
}

.chd_permissible <- function(e, adj, snode) {
  common <- intersect(adj[[e[1L]]], adj[[e[2L]]])
  for (nb in common) {
    i1 <- intersect(snode[[e[1L]]], snode[[nb]]); i2 <- intersect(snode[[e[2L]]], snode[[nb]])
    if (length(i1) != length(i2) || !all(i1 %in% i2)) return(FALSE)
  }
  TRUE
}

# components of the separator graph: cliques joined when their intersection
# is not exactly the separator
.chd_components <- function(ci, sep, snode) {
  n <- length(ci)
  H <- rep(list(integer()), n)
  for (i in seq_len(n)) for (j in .cl_rng(i + 1L, n)) {
    s1 <- snode[[ci[i]]]; s2 <- snode[[ci[j]]]
    inter <- s1[s1 %in% s2]
    equal <- length(inter) == length(sep) && all(inter %in% sep)
    if (!equal) { H[[i]] <- c(H[[i]], j); H[[j]] <- c(H[[j]], i) }
  }
  comp <- integer(n); k <- 0L
  for (v in seq_len(n)) if (!comp[v]) {
    k <- k + 1L
    stack <- v
    while (length(stack)) {
      u <- stack[length(stack)]; stack <- stack[-length(stack)]
      if (comp[u]) next
      comp[u] <- k
      stack <- c(stack, H[[u]][!comp[H[[u]]]])
    }
  }
  comp
}

.chd_clique_tree <- function(t, W, S) {
  N <- nrow(W)
  # clique_intersections!
  for (k in which(S)) {
    r <- (k - 1L) %% N + 1L; cc <- (k - 1L) %/% N + 1L
    W[r, cc] <- sum(t$snode[[r]] %in% t$snode[[cc]])
  }
  # kruskal!: the heaviest edges of a spanning forest become -1
  nz <- which(S)
  perm <- order(-W[nz])
  set <- seq_len(N)
  find <- function(x) { while (set[x] != x) x <- set[x]; x }
  found <- 0L
  for (k in perm) {
    r <- (nz[k] - 1L) %% N + 1L; cc <- (nz[k] - 1L) %/% N + 1L
    a <- find(r); b <- find(cc)
    if (a != b) {
      set[b] <- a
      W[r, cc] <- -1L
      found <- found + 1L
      if (found >= t$n_cliques - 1L) break
    }
  }
  # determine_parent_cliques!
  v <- t$post[length(t$post)]
  root <- 0L
  for (k in seq_along(t$snode)) if (v %in% t$snode[[k]]) {
    t$snode_parent[k] <- .CHD_NO_PARENT
    root <- k
    break
  }
  stack <- root
  while (length(stack)) {
    cq <- stack[length(stack)]; stack <- stack[-length(stack)]
    nbrs <- c(if (cq > 1L) which(W[cq, seq_len(cq - 1L)] != 0L), if (cq < N) which(S[, cq]))
    for (nb in nbrs) {
      if (W[max(cq, nb), min(cq, nb)] == -1L && t$snode_parent[cq] != nb) {
        t$snode_parent[nb] <- cq
        t$snode_children[[cq]] <- .chd_push(t$snode_children[[cq]], nb)
        stack <- c(stack, nb)
      }
    }
  }
  po <- .chd_post_order(t$snode_parent, t$snode_children, t$n_cliques)
  t$snode_post <- po$post; t$snode_children <- po$children
  t$separators <- rep(list(integer()), length(t$separators))
  for (j in seq_len(t$n_cliques - 1L)) {
    cq <- t$snode_post[j]
    p <- t$snode_parent[cq]
    t$separators[[cq]] <- intersect(t$snode[[cq]], t$snode[[p]])
    t$snode[[cq]] <- t$snode[[cq]][!(t$snode[[cq]] %in% t$separators[[cq]])]
  }
  invisible()
}

# ---------------------------------------------------------------------------
# Augmentation

.chd_get_clique <- function(t, i) { cq <- t$snode_post[i]; .chd_union(t$snode[[cq]], t$separators[[cq]]) }

.chd_decomposed_dim_and_overlaps <- function(info) {
  cols <- 0L; ov <- 0L; it <- 1L
  for (k in seq_along(info$init_cones)) {
    sp <- if (it <= length(info$spatterns)) info$spatterns[[it]] else NULL
    if (!is.null(sp) && sp$orig_index == k) {
      t <- sp$sntree
      for (i in seq_len(t$n_cliques)) {
        cols <- cols + .chd_trinum(t$nblk[i])
        ov <- ov + .chd_trinum(length(t$separators[[t$snode_post[i]]]))
      }
      it <- it + 1L
    } else cols <- cols + .cl_nvars(info$init_cones[[k]])
  }
  c(cols, ov)
}

.chd_blockdiag_P <- function(P, nadd) .csc(P$m + nadd, P$n + nadd, c(P$colptr, rep(P$colptr[P$n + 1L], nadd)),
                                          P$rowval, P$nzval)

.chd_augment <- function(info, P, q, A, b, settings) {
  if (isTRUE(settings$chordal_decomposition_compact)) .chd_augment_compact(info, P, q, A, b)
  else .chd_augment_standard(info, P, q, A, b)
}

.chd_augment_standard <- function(info, P, q, A, b) {
  m <- info$init_dims[2L]
  HI <- integer()
  cones <- list(.cl_spec("ZeroConeT", m))
  row <- 1L; it <- 1L
  for (k in seq_along(info$init_cones)) {
    cone <- info$init_cones[[k]]
    sp <- if (it <= length(info$spatterns)) info$spatterns[[it]] else NULL
    if (!is.null(sp) && sp$orig_index == k) {
      t <- sp$sntree
      for (i in seq_len(t$n_cliques)) {
        v <- sort(sp$ordering[.chd_get_clique(t, i)])
        for (j in seq_along(v)) for (ii in seq_len(j)) HI <- c(HI, row + .chd_coord_to_upper(v[ii], v[j]) - 1L)
        cones[[length(cones) + 1L]] <- .cl_spec("PSDTriangleConeT", t$nblk[i])
      }
      it <- it + 1L
    } else {
      HI <- c(HI, row + seq_len(.cl_nvars(cone)) - 1L)
      cones[[length(cones) + 1L]] <- cone
    }
    row <- row + .cl_nvars(cone)
  }
  lenH <- length(HI)
  H <- .cl_sparse(HI, seq_len(lenH), rep(1, lenH), max(HI), lenH)
  info$H <- H
  # [A H; 0 -I]
  I <- c(A$rowval, H$rowval, A$m + seq_len(lenH))
  J <- c(.cl_csc_cols(A), A$n + .cl_csc_cols(H), A$n + seq_len(lenH))
  V <- c(A$nzval, H$nzval, rep(-1, lenH))
  Anew <- .cl_sparse(I, J, V, A$m + lenH, A$n + lenH)
  list(info = info, P = .chd_blockdiag_P(P, lenH), q = c(q, numeric(lenH)), A = Anew, b = c(b, numeric(lenH)),
       cones = cones)
}

.chd_rows_subset <- function(rows, lo, hi) {
  if (!length(rows)) return(NULL)
  s <- .chd_ssf(rows, lo)
  if (s == length(rows) + 1L) return(NULL)
  if (rows[s] > hi || s == 0L) return(NULL)
  c(s, .chd_ssl(rows, hi))
}
.chd_rows_mat <- function(A, col, lo, hi) {
  a <- A$colptr[col]; z <- A$colptr[col + 1L] - 1L
  if (z < a) return(NULL)
  se <- .chd_rows_subset(A$rowval[a:z], lo, hi)
  if (is.null(se)) NULL else c(a + se[1L] - 1L, a + se[2L] - 1L)
}

.chd_row_index <- function(k, rowval, rlo, rc) {
  if (is.null(rc) || rc[2L] < rc[1L]) return(0L)
  ks <- rlo + k - 1L
  u <- min(rc[2L], rc[1L] + ks - 1L)
  r <- .chd_ssf(rowval, ks, rc[1L], u)
  if (r > u || rowval[r] != ks) 0L else r
}

.chd_block_indices <- function(snode, sep, nv) {
  bi <- list()
  for (j in sep) for (i in sep) if (i <= j) bi[[length(bi) + 1L]] <- c(i, j, 1L)
  for (j in snode) for (i in snode) if (i <= j) bi[[length(bi) + 1L]] <- c(i, j, 0L)
  for (i in snode) for (j in sep) bi[[length(bi) + 1L]] <- c(min(i, j), max(i, j), 0L)
  bi[order(vapply(bi, function(x) x[2L] * nv + x[1L], 0))]
}

.chd_augment_compact <- function(info, P, q, A, b) {
  cones0 <- info$init_cones
  do <- .chd_decomposed_dim_and_overlaps(info)
  Am <- do[1L]; nov <- do[2L]; An <- A$n + nov
  nnzA <- length(A$nzval)
  tot <- nnzA + 2L * nov
  AI <- integer(tot)
  AJ <- integer(tot); AJ[seq_len(nnzA)] <- .cl_csc_cols(A)
  AV <- rep(1, tot); AV[seq_len(nnzA)] <- A$nzval
  sv <- A$n + 1L
  for (i in seq.int(nnzA + 1L, length.out = nov)) {
    pos <- nnzA + 2L * (i - nnzA) - 1L
    AJ[pos] <- sv; AJ[pos + 1L] <- sv; AV[pos + 1L] <- -1
    sv <- sv + 1L
  }
  bnz <- which(b != 0); bval <- b[bnz]
  bI <- integer(length(bnz))
  rng <- .chd_rng_cones(cones0)
  cones <- list(); maps <- list()
  row_ptr <- 1L; ov_ptr <- nnzA + 1L; it <- 1L
  for (k in seq_along(cones0)) {
    cone <- cones0[[k]]; rr <- rng[[k]]
    rlo <- rr[1L]; rhi <- rr[length(rr)]
    sp <- if (it <= length(info$spatterns)) info$spatterns[[it]] else NULL
    if (!is.null(sp) && sp$orig_index == k) {
      t <- sp$sntree; ord <- sp$ordering
      # clique_rows_map
      crows <- list(); rs <- row_ptr
      for (i in rev(seq_len(t$n_cliques))) {
        nr <- .chd_trinum(t$nblk[i])
        crows[[as.character(t$snode_post[i])]] <- rs
        rs <- rs + nr
      }
      for (i in rev(seq_len(t$n_cliques))) {
        cq <- t$snode_post[i]
        sep <- sort(ord[t$separators[[cq]]])
        sn <- sort(ord[t$snode[[cq]]])
        bi <- .chd_block_indices(sn, sep, length(ord))
        if (i == t$n_cliques) { prow <- 0L; pcl <- integer() }
        else {
          pidx <- t$snode_parent[cq]
          prow <- crows[[as.character(pidx)]]
          pcl <- sort(ord[.chd_union(t$snode[[pidx]], t$separators[[pidx]])])
        }
        for (col in seq_len(A$n)) {
          rc <- .chd_rows_mat(A, col, rlo, rhi)
          rb <- if (col == 1L) .chd_rows_subset(bnz, rlo, rhi) else NULL
          counter <- 0L
          for (x in bi) {
            nr <- row_ptr + counter
            if (x[3L] == 1L) {
              if (col == 1L) {
                AI[ov_ptr] <- nr
                AI[ov_ptr + 1L] <- prow + .chd_coord_to_upper(.chd_ssf(pcl, x[1L]), .chd_ssf(pcl, x[2L])) - 1L
                ov_ptr <- ov_ptr + 2L
              }
            } else {
              kk <- .chd_coord_to_upper(x[1L], x[2L])
              r0 <- .chd_row_index(kk, A$rowval, rlo, rc)
              if (r0 != 0L) AI[r0] <- nr
              if (col == 1L) {
                r0 <- .chd_row_index(kk, bnz, rlo, rb)
                if (r0 != 0L) bI[r0] <- nr
              }
            }
            counter <- counter + 1L
          }
        }
        cones[[length(cones) + 1L]] <- .cl_spec("PSDTriangleConeT", t$nblk[i])
        maps[[length(maps) + 1L]] <- list(orig_index = sp$orig_index, tree = it, clique = i)
        row_ptr <- row_ptr + .chd_trinum(t$nblk[i])
      }
      it <- it + 1L
    } else {
      off <- row_ptr - rlo
      rb <- .chd_rows_subset(bnz, rlo, rhi)
      if (!is.null(rb)) for (j in .cl_rng(rb[1L], rb[2L])) bI[j] <- bnz[j] + off
      for (col in seq_len(A$n)) {
        rc <- .chd_rows_mat(A, col, rlo, rhi)
        if (!is.null(rc)) for (j in .cl_rng(rc[1L], rc[2L])) AI[j] <- A$rowval[j] + off
      }
      cones[[length(cones) + 1L]] <- cone
      oi <- if (!length(maps)) 1L else maps[[length(maps)]]$orig_index + 1L
      maps[[length(maps) + 1L]] <- list(orig_index = oi, tree = NULL, clique = NULL)
      row_ptr <- row_ptr + .cl_nvars(cone)
    }
  }
  info$cone_maps <- maps
  Anew <- .cl_sparse(AI, AJ, AV, Am, An)
  bnew <- numeric(Am)
  for (j in seq_along(bI)) bnew[bI[j]] <- bnew[bI[j]] + bval[j]
  nadd <- An - A$n
  list(info = info, P = .chd_blockdiag_P(P, nadd), q = c(q, numeric(nadd)), A = Anew, b = bnew, cones = cones)
}

# ---------------------------------------------------------------------------
# Reverse map of the solution

.chd_reverse <- function(info, x, s, z, cones, settings) {
  n <- info$init_dims[1L]; m <- info$init_dims[2L]
  nx <- x[seq_len(n)]
  ns <- numeric(m); nz <- numeric(m)
  if (isTRUE(settings$chordal_decomposition_compact)) {
    rng <- .chd_rng_cones(info$init_cones)
    row_ptr <- 1L
    for (k in seq_along(cones)) {
      mp <- info$cone_maps[[k]]
      rr <- rng[[mp$orig_index]]
      if (is.null(mp$tree)) {
        src <- row_ptr + seq_len(.cl_nvars(cones[[k]])) - 1L
        ns[rr] <- s[src]; nz[rr] <- z[src]
        row_ptr <- row_ptr + .cl_nvars(cones[[k]])
      } else {
        sp <- info$spatterns[[mp$tree]]
        cb <- sort(sp$ordering[.chd_get_clique(sp$sntree, mp$clique)])
        counter <- 0L
        for (j in cb) for (i in cb) if (i <= j) {
          off <- .chd_coord_to_upper(i, j) - 1L
          ns[rr[1L] + off] <- ns[rr[1L] + off] + s[row_ptr + counter]
          nz[rr[1L] + off] <- z[row_ptr + counter]
          counter <- counter + 1L
        }
        row_ptr <- row_ptr + .chd_trinum(length(cb))
      }
    }
  } else {
    H <- info$H
    sd <- s[-seq_len(m)]; zd <- z[-seq_len(m)]
    for (col in seq_len(H$n)) for (k in .cl_rng(H$colptr[col], H$colptr[col + 1L] - 1L)) {
      r <- H$rowval[k]
      ns[r] <- ns[r] + H$nzval[k] * sd[col]
      nz[r] <- nz[r] + H$nzval[k] * zd[col]
    }
    cnt <- tabulate(H$rowval, H$m)
    over <- which(cnt > 1)
    nz[over] <- nz[over] / cnt[over]
  }
  if (isTRUE(settings$chordal_decomposition_complete_dual)) {
    rng <- .chd_rng_cones(info$init_cones)
    for (sp in info$spatterns) {
      rr <- rng[[sp$orig_index]]
      nz[rr] <- .chd_complete(nz[rr], sp)
    }
  }
  list(x = nx, s = ns, z = nz)
}

.chd_complete <- function(z, sp) {
  n <- length(sp$ordering)
  Z <- .clk_svec_to_mat(z, n)
  t <- sp$sntree
  p <- sp$ordering
  ip <- integer(n); ip[p] <- seq_len(n)
  W <- Z[p, p, drop = FALSE]
  for (j in rev(seq_len(t$n_cliques - 1L))) {
    cq <- t$snode_post[j]
    nu <- t$snode[[cq]]; al <- t$separators[[cq]]
    i <- nu[1L]
    eta <- .cl_rng(i + 1L, n)
    eta <- eta[!vapply(eta, function(x) .chd_insorted(x, al) || .chd_insorted(x, nu), TRUE)]
    Waa <- W[al, al, drop = FALSE]
    Wan <- W[al, nu, drop = FALSE]
    # Julia's backslash fails only on an exactly singular block, then pinv
    Y <- tryCatch(solve(Waa, Wan, tol = 0), error = function(e) {
      sv <- svd(Waa)
      tol <- min(dim(Waa)) * .Machine$double.eps * max(sv$d, 0)
      dinv <- ifelse(sv$d > tol, 1 / sv$d, 0)
      (sv$v %*% (dinv * t(sv$u))) %*% Wan
    })
    W[eta, nu] <- W[eta, al, drop = FALSE] %*% Y
    W[nu, eta] <- t(W[eta, nu, drop = FALSE])
  }
  .clk_mat_to_svec(W[ip, ip, drop = FALSE])
}
