# Sparsity detection and coloring, as reexported by NonlinearSolve from
# ADTypes.jl (jacobian_sparsity, hessian_sparsity, column_coloring,
# row_coloring, symmetric_coloring, src/sparse.jl) with the methods of
# SparseConnectivityTracer.jl (the global gradient and Hessian tracers of
# TracerSparsityDetector: src/operators.jl, src/overloads/gradient_tracer.jl,
# src/overloads/hessian_tracer.jl) and SparseMatrixColorings.jl
# (GreedyColoringAlgorithm with the natural order: partial distance-2
# coloring and star coloring, src/coloring.jl and src/graph.jl).

#' Sparsity patterns and colorings
#'
#' `jacobian_sparsity(f, x, sd)` is the pattern of the Jacobian of `f` at
#' inputs like `x`, and `hessian_sparsity(f, x, sd)` the pattern of the
#' Hessian of the scalar function `f`, as logical sparse matrices. With the
#' default `sd = NULL` the pattern is found by tracing `f`, as
#' SparseConnectivityTracer's `TracerSparsityDetector()`: `f` is called with
#' a vector whose entries record which inputs they depend on, through the
#' arithmetic, math and summary functions of R. The pattern is global (it
#' does not depend on the values in `x`), and code whose control flow depends
#' on the values cannot be traced. `sd = "none"` gives a dense pattern
#' (`NoSparsityDetector()`), and a matrix gives that known pattern
#' (`KnownJacobianSparsityDetector()`, `KnownHessianSparsityDetector()`).
#'
#' `column_coloring(M)` partitions the columns of `M` into structurally
#' orthogonal groups (no two columns of a group have a nonzero in the same
#' row), `row_coloring(M)` the rows, and `symmetric_coloring(M)` the columns
#' of a symmetric matrix by a star coloring, so that a Jacobian or Hessian
#' can be computed from one directional derivative per color. They are the
#' greedy algorithms of SparseMatrixColorings with the natural vertex order
#' (`GreedyColoringAlgorithm()`); `ca = "none"` gives every column its own
#' color (`NoColoringAlgorithm()`).
#'
#' @param f A function of a numeric vector.
#' @param x The input vector.
#' @param sd The sparsity detector: `NULL`, `"none"` or a known pattern.
#' @param M A matrix (dense or from the Matrix package); its nonzero entries
#'   are the pattern.
#' @param ca The coloring algorithm: `NULL` (greedy) or `"none"`.
#' @return A logical sparse matrix (`jacobian_sparsity()`,
#'   `hessian_sparsity()`) or an integer vector of colors.
#' @examples
#' f <- function(x) c(x[1] * x[2], sin(x[3]), x[1] + x[4] * x[4])
#' S <- jacobian_sparsity(f, c(1, 2, 3, 4))
#' S
#' column_coloring(S)
#' row_coloring(S)
#' g <- function(x) x[1] * x[2] + exp(x[3]) + x[4] / x[1]
#' H <- hessian_sparsity(g, c(1, 2, 3, 4))
#' symmetric_coloring(H)
#' @name sparsity_coloring
NULL

# ---------------------------------------------------------------------------
# The tracer: a vector whose elements carry a gradient pattern (input
# indices) and a Hessian pattern (index pairs "i j").

.trel <- function(g = integer(), h = character()) structure(list(g = g, h = h), class = "nl_trel")
.trc <- function(els) structure(els, class = "nl_tracer")

.trc_elems <- function(z) {
  if (is.null(z)) return(list())
  if (inherits(z, "nl_tracer")) return(unclass(z))
  if (inherits(z, "nl_trel")) return(list(z))
  if (is.numeric(z) || is.logical(z)) return(rep(list(.trel()), length(z)))
  if (is.list(z)) return(do.call(c, c(list(list()), lapply(unname(z), .trc_elems))))
  stop("cannot trace a value of class ", class(z)[1L], call. = FALSE)
}

.trc_prod <- function(a, b) if (length(a) && length(b)) as.vector(outer(a, b, paste)) else character()

# 1 to 1 operator with the derivative flags of operators.jl
.trc_unary <- function(t, d1, d2) {
  g <- if (d1) t$g else integer()
  h <- if (d1 && d2) union(t$h, .trc_prod(t$g, t$g))
       else if (d1) t$h
       else if (d2) .trc_prod(t$g, t$g)
       else character()
  .trel(g, h)
}

# 2 to 1 operator: d1x, d2x, d1y, d2y, cross are TRUE when nonzero
.trc_binary <- function(x, y, d1x, d2x, d1y, d2y, cross) {
  g <- sort(union(if (d1x) x$g, if (d1y) y$g))
  h <- c(if (d1x) x$h, if (d1y) y$h, if (d2x) .trc_prod(x$g, x$g), if (d2y) .trc_prod(y$g, y$g),
         if (cross) c(.trc_prod(x$g, y$g), .trc_prod(y$g, x$g)))
  .trel(as.integer(g), unique(h))
}

# operator classes: s (f' and f'' nonzero), f (only f'), z (none)
.trc_math_class <- function(name) {
  if (name %in% c("abs")) return("f")
  if (name %in% c("sign", "floor", "ceiling", "trunc", "round", "signif")) return("z")
  "s"
}

# binary flags c(d1x, d2x, d1y, d2y, cross)
.trc_ops_flags <- function(op) switch(op,
  "^" = c(TRUE, TRUE, TRUE, TRUE, TRUE),
  "/" = c(TRUE, FALSE, TRUE, TRUE, TRUE),
  "*" = c(TRUE, FALSE, TRUE, FALSE, TRUE),
  "+" = , "-" = , "%%" = c(TRUE, FALSE, TRUE, FALSE, FALSE),
  c(FALSE, FALSE, FALSE, FALSE, FALSE))

#' @export
`[.nl_tracer` <- function(x, i, ...) {
  if (missing(i)) return(x)
  out <- unclass(x)[i]
  if (any(vapply(out, is.null, TRUE))) stop("index out of bounds in a traced vector", call. = FALSE)
  .trc(out)
}

#' @export
`[[.nl_tracer` <- function(x, i, ...) .trc(list(unclass(x)[[i]]))

#' @export
`[<-.nl_tracer` <- function(x, i, value) {
  y <- unclass(x)
  y[i] <- .trc_elems(value)
  .trc(y)
}

#' @export
`[[<-.nl_tracer` <- function(x, i, value) {
  y <- unclass(x)
  v <- .trc_elems(value)
  if (length(v) != 1L) stop("replacement has length ", length(v), call. = FALSE)
  y[[i]] <- v[[1L]]
  .trc(y)
}

#' @export
c.nl_tracer <- function(...) .trc(.trc_elems(list(...)))

#' @export
as.list.nl_tracer <- function(x, ...) lapply(unclass(x), function(e) .trc(list(e)))

#' @export
rev.nl_tracer <- function(x) .trc(rev(unclass(x)))

#' @export
Ops.nl_tracer <- function(e1, e2) {
  a <- .trc_elems(e1)
  if (missing(e2)) {
    if (.Generic %in% c("+", "-")) return(.trc(lapply(a, .trc_unary, TRUE, FALSE)))
    return(.trc(lapply(a, .trc_unary, FALSE, FALSE)))
  }
  b <- .trc_elems(e2)
  na <- length(a); nb <- length(b)
  if (!na || !nb) return(.trc(list()))
  fl <- .trc_ops_flags(.Generic)
  .trc(lapply(seq_len(max(na, nb)), function(i)
    .trc_binary(a[[(i - 1L) %% na + 1L]], b[[(i - 1L) %% nb + 1L]], fl[1], fl[2], fl[3], fl[4], fl[5])))
}

#' @export
chooseOpsMethod.nl_tracer <- function(x, y, mx, my, cl, reverse) TRUE

#' @export
Math.nl_tracer <- function(x, ...) {
  a <- unclass(x)
  switch(.Generic,
    cumsum = , cummax = , cummin = {
      acc <- .trel()
      return(.trc(lapply(a, function(e) acc <<- .trc_binary(acc, e, TRUE, FALSE, TRUE, FALSE, FALSE))))
    },
    cumprod = {
      acc <- NULL
      return(.trc(lapply(a, function(e)
        acc <<- if (is.null(acc)) e else .trc_binary(acc, e, TRUE, FALSE, TRUE, FALSE, TRUE))))
    })
  cl <- .trc_math_class(.Generic)
  .trc(lapply(a, .trc_unary, cl != "z", cl == "s"))
}

#' @export
Summary.nl_tracer <- function(..., na.rm = FALSE) {
  a <- .trc_elems(list(...))
  fold <- function(fl) {
    if (!length(a)) return(.trel())
    Reduce(function(x, y) .trc_binary(x, y, fl[1], fl[2], fl[3], fl[4], fl[5]), a)
  }
  switch(.Generic,
    sum = , max = , min = .trc(list(fold(c(TRUE, FALSE, TRUE, FALSE, FALSE)))),
    prod = .trc(list(fold(c(TRUE, FALSE, TRUE, FALSE, TRUE)))),
    range = { r <- fold(c(TRUE, FALSE, TRUE, FALSE, FALSE)); .trc(list(r, r)) },
    .trc(list(.trel())))
}

.trc_input <- function(x) .trc(lapply(seq_along(x), function(i) .trel(i)))

.trc_pattern <- function(rows, cols, m, n) {
  Matrix::sparseMatrix(i = as.integer(rows), j = as.integer(cols), x = rep(TRUE, length(rows)), dims = c(m, n))
}

.trc_check_known <- function(sd, m, n, what) {
  if (!all(dim(sd) == c(m, n)))
    stop(sprintf("%s size (%d, %d) of the known pattern doesn't match expected size (%d, %d).",
                 what, nrow(sd), ncol(sd), m, n), call. = FALSE)
  sd
}

#' @rdname sparsity_coloring
#' @export
jacobian_sparsity <- function(f, x, sd = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  if (identical(sd, "none")) return(.trc_pattern(rep(seq_len(length(f(x))), n), rep(seq_len(n), each = length(f(x))),
                                                 length(f(x)), n))
  if (!is.null(sd)) return(.trc_check_known(sd, length(f(x)), n, "Jacobian"))
  y <- .trc_elems(f(.trc_input(x)))
  rows <- rep(seq_along(y), vapply(y, function(e) length(e$g), 0L))
  cols <- unlist(lapply(y, `[[`, "g"))
  .trc_pattern(rows, cols, length(y), n)
}

#' @rdname sparsity_coloring
#' @export
hessian_sparsity <- function(f, x, sd = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  if (identical(sd, "none")) return(.trc_pattern(rep(seq_len(n), n), rep(seq_len(n), each = n), n, n))
  if (!is.null(sd)) return(.trc_check_known(sd, n, n, "Hessian"))
  y <- .trc_elems(f(.trc_input(x)))
  if (length(y) != 1L) stop("hessian_sparsity needs a scalar function", call. = FALSE)
  keys <- y[[1L]]$h
  if (!length(keys)) return(.trc_pattern(integer(), integer(), n, n))
  ij <- matrix(as.integer(unlist(strsplit(keys, " "))), ncol = 2L, byrow = TRUE)
  .trc_pattern(c(ij[, 1L], ij[, 2L]), c(ij[, 2L], ij[, 1L]), n, n)
}

# ---------------------------------------------------------------------------
# Colorings

# column -> sorted rows of the nonzeros (a CSC pattern)
.col_pattern <- function(M) {
  P <- if (inherits(M, "Matrix")) as.matrix(M != 0) else (as.matrix(M) != 0)
  P[is.na(P)] <- TRUE
  list(cols = lapply(seq_len(ncol(P)), function(j) which(P[, j])),
       rows = lapply(seq_len(nrow(P)), function(i) which(P[i, ])), m = nrow(P), n = ncol(P))
}

.partial_distance2 <- function(own, other, nv) {
  color <- integer(nv)
  forbidden <- integer(nv)
  for (v in seq_len(nv)) {
    for (w in own[[v]]) for (x in other[[w]]) if (color[x] != 0L) forbidden[color[x]] <- v
    for (i in seq_len(nv)) if (forbidden[i] != v) { color[v] <- i; break }
  }
  color
}

#' @rdname sparsity_coloring
#' @export
column_coloring <- function(M, ca = NULL) {
  if (identical(ca, "none")) return(seq_len(ncol(M)))
  p <- .col_pattern(M)
  .partial_distance2(p$cols, p$rows, p$n)
}

#' @rdname sparsity_coloring
#' @export
row_coloring <- function(M, ca = NULL) {
  if (identical(ca, "none")) return(seq_len(nrow(M)))
  p <- .col_pattern(M)
  .partial_distance2(p$rows, p$cols, p$m)
}

#' @rdname sparsity_coloring
#' @export
symmetric_coloring <- function(M, ca = NULL) {
  if (nrow(M) != ncol(M)) stop("symmetric_coloring needs a square matrix", call. = FALSE)
  if (identical(ca, "none")) return(seq_len(nrow(M)))
  p <- .col_pattern(M)
  nv <- p$n
  # CSC positions and edge indices (build_edge_to_index)
  colptr <- c(1L, 1L + cumsum(lengths(p$cols)))
  rowval <- unlist(p$cols)
  edge <- integer(length(rowval))
  offsets <- integer(nv)
  counter <- 0L
  for (j in seq_len(nv)) for (k in seq.int(colptr[j], length.out = length(p$cols[[j]]))) {
    i <- rowval[k]
    if (i > j) {
      counter <- counter + 1L
      edge[k] <- counter
      k2 <- colptr[i] + offsets[i]
      edge[k2] <- counter
      offsets[i] <- offsets[i] + 1L
    }
  }
  nbr <- function(v) {
    r <- seq.int(colptr[v], length.out = length(p$cols[[v]]))
    list(w = rowval[r], e = edge[r])
  }
  color <- integer(nv); forbidden <- integer(nv); treated <- integer(nv)
  fn_p <- integer(nv); fn_q <- integer(nv); fn_e <- integer(nv)
  star <- integer(length(rowval)); hub <- integer()
  treat <- function(v, w) {
    for (x in nbr(w)$w) {
      if (w == x || color[x] == 0L) next
      forbidden[color[x]] <<- v
    }
    treated[w] <<- v
  }
  for (v in seq_len(nv)) {
    nv_ <- nbr(v)
    for (k in seq_along(nv_$w)) {
      w <- nv_$w[k]; ivw <- nv_$e[k]
      if (v == w) next
      cw <- color[w]
      if (cw == 0L) next
      forbidden[cw] <- v
      pp <- fn_p[cw]; qq <- fn_q[cw]
      if (pp == v) {
        if (treated[qq] != v) treat(v, qq)
        treat(v, w)
      } else {
        fn_p[cw] <- v; fn_q[cw] <- w; fn_e[cw] <- ivw
        nw <- nbr(w)
        for (kk in seq_along(nw$w)) {
          x <- nw$w[kk]; iwx <- nw$e[kk]
          if (w == x) next
          cx <- color[x]
          if (x == v || cx == 0L) next
          if (star[iwx] > 0L && x == hub[star[iwx]]) forbidden[cx] <- v
        }
      }
    }
    for (i in seq_len(nv)) if (forbidden[i] != v) { color[v] <- i; break }
    # _update_stars!
    for (k in seq_along(nv_$w)) {
      w <- nv_$w[k]; ivw <- nv_$e[k]
      if (v == w) next
      cw <- color[w]
      if (cw == 0L) next
      x_exists <- FALSE
      nw <- nbr(w)
      for (kk in seq_along(nw$w)) {
        x <- nw$w[kk]; iwx <- nw$e[kk]
        if (w == x) next
        if (x != v && color[x] == color[v]) {
          swx <- star[iwx]
          hub[swx] <- w
          star[ivw] <- swx
          x_exists <- TRUE
          break
        }
      }
      if (!x_exists) {
        pp <- fn_p[cw]; qq <- fn_q[cw]; ipq <- fn_e[cw]
        if (pp == v && qq != w) {
          svq <- star[ipq]
          hub[svq] <- v
          star[ivw] <- svq
        } else {
          hub <- c(hub, -max(v, w))
          star[ivw] <- length(hub)
        }
      }
    }
  }
  color
}
