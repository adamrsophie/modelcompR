# Deterministic global optimization: the OptimizationQuadDIRECT wrapper of
# Optimization.jl with the QuadDIRECT.jl algorithm (the box tree, the
# minimum edge sweeps and the quasi Newton steps on regression models),
# ported to R.

#' QuadDIRECT global optimization
#'
#' `QuadDirect()` minimizes with QuadDIRECT.jl, ported to R: a deterministic
#' global method that splits the search box into a tree of boxes, each
#' holding one function value. Every sweep splits, along each coordinate,
#' the boxes on the lower edge of the width versus value plane (the
#' minimum edge of DIRECT, with values extrapolated by one dimensional
#' quadratic fits); then, in the boxes that look convex, a quadratic model
#' of the function is fitted by regression on the points of the tree and a
#' quasi Newton step with a bounded least squares solve and a backtracking
#' line search is tried. No derivatives are used.
#'
#' The problem needs bounds (`lb` and `ub`; infinite bounds are allowed)
#' and `solve` needs `splits`, a list with one increasing triple per
#' unknown: the initial values along that coordinate, inside the bounds.
#' The first evaluation point is made of the middle values. The `solve`
#' arguments map as in OptimizationQuadDIRECT: `maxiters` to the maximal
#' number of evaluations (default 2500; the count may run a little over),
#' `abstol` to `atol` (0) and `reltol` to `rtol` (0.001), the changes of
#' the minimum below which a sweep counts as converged; the run stops after
#' more such sweeps in a row than there are unknowns. `maxtime` is not
#' used. The return code is `Default` and the solution's `original` is the
#' root of the box tree (an environment with the fields `splitdim`,
#' `xvalues`, `fvalues` and `children`).
#'
#' @return An algorithm object.
#' @examples
#' rosenbrock <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
#' prob <- OptimizationProblem(rosenbrock, c(0, 0), c(1, 100), lb = c(-1, -1), ub = c(0.8, 0.8))
#' sol <- solve(prob, QuadDirect(), splits = list(c(-0.5, 0, 0.5), c(-0.5, 0, 0.5)),
#'              maxiters = 200)
#' round(sol$u, 3); sol$objective
#' @export
QuadDirect <- function()
  structure(list(name = "QuadDirect", traits = list(bounds = TRUE, rbounds = TRUE, callback = FALSE)),
            class = c("QuadDirect", "OptimizationAlgorithm"))

.mtk$hooks$opt_solve_QuadDirect <- function(a, alg, opts) {
  .opt_check_traits(a, alg)
  if (is.null(opts$splits))
    stop("You must provide the initial locations at which to evaluate the function in `splits` (a list of ",
         "3-vectors with values in strictly increasing order and within the specified bounds).", call. = FALSE)
  maxiters <- opts$maxiters
  if (!is.null(maxiters)) {
    if (maxiters <= 0) stop("The number of maxiters has to be a non-negative and non-zero number.", call. = FALSE)
    maxiters <- round(maxiters)
  }
  if (!is.null(opts$maxtime)) warning("common maxtime is currently not used by QuadDirect()", call. = FALSE)
  p <- a$p
  loss <- function(x) as.numeric(a$f$f(x, p))[1L]
  r <- .qd_analyze(loss, lapply(opts$splits, as.numeric), as.numeric(a$lb), as.numeric(a$ub),
                   rtol = opts$reltol %||% 1e-3, atol = opts$abstol %||% 0, maxevals = maxiters %||% 2500)
  box <- .qd_minimum(r$root)
  .opt_solution(a, alg, .qd_position(box, r$x0), .qd_value(box), "Default", r$root)
}

# ---------------------------------------------------------------------------
# Floating point as in Julia and its OpenBLAS on this platform: a fused
# multiply add, the positive definite Cholesky factor of
# PositiveFactorizations (its rank one updates are fused), the LAPACK solve
# with that factor (OpenBLAS trsm: blocks of 8 rows, reciprocal diagonals,
# fused updates inside a block and fused dot products between blocks) and
# the product L L'.

.qd_eps <- function(x) {
  x <- abs(x)
  if (is.na(x) || is.infinite(x)) return(NaN)
  if (x < 0x1.0000000000000p-1022) return(4.9406564584124654e-324)
  e <- floor(log2(x))
  if (2^e > x) e <- e - 1 else if (2^(e + 1) <= x) e <- e + 1
  2^(e - 52)
}

.qd_fma <- function(a, b, c) {
  p <- a * b
  sa <- 134217729 * a; ah <- sa - (sa - a); al <- a - ah
  sb <- 134217729 * b; bh <- sb - (sb - b); bl <- b - bh
  e <- ((ah * bh - p) + ah * bl + al * bh) + al * bl
  t <- p + c
  bb <- t - p
  r <- t + (((p - (t - bb)) + (c - bb)) + e)
  bad <- !is.finite(r)
  r[bad] <- t[bad]
  r
}

.qd_poschol <- function(B) {
  K <- nrow(B)
  tol <- 10 * K * 2.220446049250313e-16 * max(abs(B))
  for (j in seq_len(K)) {
    Bjj <- B[j, j]
    below <- if (j < K) (j + 1L):K else integer(0)
    if (isTRUE(abs(Bjj) > tol)) {
      d <- sign(Bjj)
      s <- sqrt(abs(Bjj))
      B[j, j] <- s
      B[below, j] <- B[below, j] * (d / s)
      for (jj in below) {
        ii <- jj:K
        B[ii, jj] <- .qd_fma(B[ii, j], -d * B[jj, j], B[ii, jj])
      }
    } else {
      B[j, j] <- 1
      B[below, j] <- 0
    }
  }
  B[upper.tri(B)] <- 0
  B
}

.qd_trsm_blocks <- function(m) {
  bl <- list()
  s <- 1L
  for (k in seq_len(m %/% 8L)) {
    bl[[length(bl) + 1L]] <- s:(s + 7L)
    s <- s + 8L
  }
  for (i in c(4L, 2L, 1L)) if (bitwAnd(m, i)) {
    bl[[length(bl) + 1L]] <- s:(s + i - 1L)
    s <- s + i
  }
  bl
}

# the solution of (L L') z = g
.qd_cholsolve <- function(L, g) {
  n <- length(g)
  bl <- .qd_trsm_blocks(n)
  y <- g
  done <- 0L
  for (b in bl) {
    if (done > 0L) for (r in b) y[r] <- y[r] - .jl_dot(L[r, seq_len(done)], y[seq_len(done)])
    for (j in b) {
      y[j] <- y[j] * (1 / L[j, j])
      ii <- b[b > j]
      if (length(ii)) y[ii] <- .qd_fma(-y[j], L[ii, j], y[ii])
    }
    done <- b[length(b)]
  }
  for (b in rev(bl)) {
    hi <- b[length(b)]
    if (hi < n) {
      k <- (hi + 1L):n
      for (r in b) y[r] <- y[r] - .jl_dot(L[k, r], y[k])
    }
    for (j in rev(b)) {
      y[j] <- y[j] * (1 / L[j, j])
      ii <- b[b < j]
      if (length(ii)) y[ii] <- .qd_fma(-y[j], L[j, ii], y[ii])
    }
  }
  y
}

.qd_llt <- function(L) {
  n <- nrow(L)
  M <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(i)) M[i, j] <- M[j, i] <- .jl_dot(L[i, seq_len(j)], L[j, seq_len(j)])
  M
}

# x' B x of LinearAlgebra.dot(x, B, x)
.qd_quadform <- function(x, B) .jl_dot(vapply(seq_along(x), function(j) .jl_dot(B[, j], x), 0), x)

.qd_clamp <- function(x, lo, hi) ifelse(x > hi, hi, ifelse(x < lo, lo, x))

# ---------------------------------------------------------------------------
# The box tree (Box of QuadDIRECT): environments, the root being its own
# parent

.qd_box <- function(parent = NULL, cindex = 0L, n = parent$n) {
  b <- new.env(parent = emptyenv())
  if (is.null(parent)) {
    b$parent <- b
    b$parent_cindex <- 0L
  } else {
    b$parent <- parent
    b$parent_cindex <- cindex
    parent$children[[cindex]] <- b
  }
  b$n <- n
  b$splitdim <- 0L
  b$qnconverged <- FALSE
  b$minmax <- c(Inf, -Inf)
  b$xvalues <- c(Inf, 0, -Inf)
  b$fvalues <- c(0, 0, 0)
  b$children <- list(b, b, b)
  b
}

# only the root has parent_cindex 0 (it is its own parent)
.qd_isroot <- function(b) b$parent_cindex == 0L

.qd_value <- function(b) {
  if (.qd_isroot(b)) stop("root box does not have a unique value", call. = FALSE)
  b$parent$fvalues[b$parent_cindex]
}

.qd_get_root <- function(b) {
  while (!.qd_isroot(b)) b <- b$parent
  b
}

# Base.isless on doubles: NaN last, -0.0 before 0.0
.qd_isless <- function(a, b) {
  if (is.nan(a)) return(FALSE)
  if (is.nan(b)) return(TRUE)
  a < b || (a == 0 && b == 0 && 1 / a < 0 && 1 / b > 0)
}

.qd_preorder <- function(root) {
  out <- vector("list", 64L)
  k <- 0L
  stack <- list(root)
  while (length(stack)) {
    b <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL
    k <- k + 1L
    if (k > length(out)) length(out) <- 2L * length(out)
    out[[k]] <- b
    if (b$splitdim != 0L) stack <- c(stack, list(b$children[[3L]], b$children[[2L]], b$children[[1L]]))
  }
  out[seq_len(k)]
}

# leaves(box): empty for a leaf, as in QuadDIRECT
.qd_leaves <- function(b) {
  if (b$splitdim == 0L) return(list())
  Filter(function(x) x$splitdim == 0L, .qd_preorder(b))
}

# minimum(root): the first box in depth first order of smallest value
.qd_minimum <- function(root) {
  nodes <- .qd_preorder(root)
  best <- root
  bv <- Inf
  for (b in nodes[-1L]) {
    v <- .qd_value(b)
    if (.qd_isless(v, bv)) {
      best <- b
      bv <- v
    }
  }
  best
}

.qd_in <- function(b, items) {
  for (x in items) if (identical(x, b)) return(TRUE)
  FALSE
}

.qd_set_add <- function(items, b) if (.qd_in(b, items)) items else c(items, list(b))

.qd_fpws <- function(b, sd) {
  while (!.qd_isroot(b)) {
    p <- b$parent
    if (p$splitdim == sd) return(b)
    b <- p
  }
  b
}

.qd_bb1 <- function(b) {
  if (.qd_isroot(b)) stop("cannot compute bounds on root Box", call. = FALSE)
  p <- b$parent
  xv <- p$xvalues
  switch(b$parent_cindex,
         c(p$minmax[1L], (xv[1L] + xv[2L]) / 2),
         c((xv[1L] + xv[2L]) / 2, (xv[2L] + xv[3L]) / 2),
         c((xv[2L] + xv[3L]) / 2, p$minmax[2L]))
}

.qd_bb <- function(b, lo, hi) if (.qd_isroot(b)) c(lo, hi) else .qd_bb1(b)

# boxbounds!(bb, flag, box): bb is an n x 2 matrix; rows not on the path
# keep their values
.qd_bb_fill <- function(bb, b) {
  n <- b$n
  flag <- logical(n)
  if (b$splitdim == 0L) {
    i <- b$parent$splitdim
    bb[i, ] <- .qd_bb1(b)
  } else {
    i <- b$splitdim
    bb[i, ] <- b$minmax
  }
  flag[i] <- TRUE
  nfilled <- 1L
  while (!.qd_isroot(b) && nfilled < n) {
    i <- b$parent$splitdim
    if (!flag[i]) {
      bb[i, ] <- .qd_bb1(b)
      flag[i] <- TRUE
      nfilled <- nfilled + 1L
    }
    b <- b$parent
  }
  bb
}

.qd_position <- function(b, x0) {
  x <- x0
  n <- length(x)
  flag <- logical(n)
  nf <- 0L
  while (!.qd_isroot(b) && nf < n) {
    i <- b$parent$splitdim
    if (!flag[i]) {
      x[i] <- b$parent$xvalues[b$parent_cindex]
      flag[i] <- TRUE
      nf <- nf + 1L
    }
    b <- b$parent
  }
  x
}

.qd_count_splits <- function(b) {
  ns <- integer(b$n)
  if (b$splitdim == 0L) {
    b <- b$parent
    if (b$splitdim == 0L) return(ns)
  }
  ns[b$splitdim] <- ns[b$splitdim] + 1L
  while (!.qd_isroot(b)) {
    b <- b$parent
    ns[b$splitdim] <- ns[b$splitdim] + 1L
  }
  ns
}

.qd_width <- function(b, i, x0, lower, upper) {
  p <- .qd_fpws(b, i)
  bb <- .qd_bb(p, lower[i], upper[i])
  x <- if (.qd_isroot(p)) x0[i] else p$parent$xvalues[p$parent_cindex]
  max(x - bb[1L], bb[2L] - x)
}

.qd_epswidth <- function(bb) {
  w1 <- if (is.finite(bb[1L])) .qd_eps(bb[1L]) else 0
  w2 <- if (is.finite(bb[2L])) .qd_eps(bb[2L]) else 0
  10 * min(w1, w2)
}

.qd_boxscale <- function(b, splits) {
  bx <- function(s) if (s[1L] == s[2L]) s[3L] - s[2L] else if (s[2L] == s[3L]) s[2L] - s[1L] else min(s[2L] - s[1L], s[3L] - s[2L])
  vapply(seq_len(b$n), function(i) {
    p <- .qd_fpws(b, i)
    if (.qd_isroot(p)) bx(splits[[i]]) else bx(p$parent$xvalues)
  }, 0)
}

.qd_add_children <- function(parent, sd, xv, fv, u, v) {
  if (parent$splitdim != 0L) stop("cannot add children to non-leaf node", call. = FALSE)
  if (!isTRUE(xv[1L] < xv[2L] && xv[2L] < xv[3L]))
    stop("xvalues must be monotonic, got ", .qd_jlvec(xv), call. = FALSE)
  parent$splitdim <- sd
  p <- .qd_fpws(parent, sd)
  parent$minmax <- if (.qd_isroot(p)) c(u, v) else .qd_bb1(p)
  mm <- parent$minmax
  minsep <- 2.220446049250313e-16 * (xv[3L] - xv[1L])
  if (!isTRUE(all(mm[1L] <= xv & xv <= mm[2L])) || !isTRUE(xv[2L] - xv[1L] > minsep) ||
      !isTRUE(xv[3L] - xv[2L] > minsep))
    stop("AssertionError: split values outside the box or too close", call. = FALSE)
  parent$xvalues <- xv
  parent$fvalues <- fv
  for (i in 1:3) .qd_box(parent, i)
  invisible(parent)
}

.qd_find_smallest_child_leaf <- function(b) {
  vmin <- .qd_value(b)
  bmin <- b
  for (leaf in .qd_leaves(b)) {
    v <- .qd_value(leaf)
    if (isTRUE(v <= vmin)) {
      vmin <- v
      bmin <- leaf
    }
  }
  bmin
}

.qd_unique_smallest_leaves <- function(boxes) {
  out <- list()
  for (b in boxes) out <- .qd_set_add(out, .qd_find_smallest_child_leaf(b))
  out
}

.qd_find_leaf_at <- function(root, x) {
  while (root$splitdim != 0L) {
    i <- root$splitdim
    found <- FALSE
    for (b in root$children) {
      bb <- .qd_bb1(b)
      if (isTRUE(bb[1L] <= x[i] && x[i] <= bb[2L])) {
        root <- b
        found <- TRUE
        break
      }
    }
    if (!found) stop(x[i], " not within ", .qd_jlvec(root$minmax), call. = FALSE)
  }
  root
}

.qd_within <- function(x, bb, dir) {
  if (!isTRUE(bb[1L] <= x && x <= bb[2L])) return(FALSE)
  if (x == bb[1L] && dir < 0) return(FALSE)
  if (x == bb[2L] && dir > 0) return(FALSE)
  TRUE
}

.qd_find_leaf_at_edge <- function(root, x, sd, dir) {
  if (root$splitdim == 0L) return(list(root, FALSE))
  while (root$splitdim != 0L) {
    i <- root$splitdim
    found <- FALSE
    for (b in root$children) if (.qd_within(x[i], .qd_bb1(b), dir)) {
      root <- b
      found <- TRUE
      break
    }
    if (!found) return(list(root, FALSE))
  }
  list(root, TRUE)
}

# ---------------------------------------------------------------------------
# One dimensional quadratic fits

.qd_lagrangecoefs <- function(xv, fv) {
  xm <- xv[1L]; x0 <- xv[2L]; xp <- xv[3L]
  if (!isTRUE(xp > x0 && x0 > xm && is.finite(xm) && is.finite(xp)))
    stop("AssertionError: xp > x0 && x0 > xm && isfinite(xm) && isfinite(xp)", call. = FALSE)
  c(fv[1L] / ((xm - x0) * (xm - xp)), fv[2L] / ((x0 - xm) * (x0 - xp)), fv[3L] / ((xp - xm) * (xp - x0)))
}

# qfit: c(xvert, fvert, qcoef)
.qd_qfit <- function(xv, fv) {
  cf <- .qd_lagrangecoefs(xv, fv)
  xm <- xv[1L]; x0 <- xv[2L]; xp <- xv[3L]
  qcoef <- cf[1L] + cf[2L] + cf[3L]
  if (isTRUE(fv[1L] == fv[2L] && fv[2L] == fv[3L])) return(c(x0, fv[2L], 0))
  xvert <- (cf[1L] * (x0 + xp) + cf[2L] * (xm + xp) + cf[3L] * (xm + x0)) / (2 * qcoef)
  qv <- cf[1L] * (xvert - x0) * (xvert - xp) + cf[2L] * (xvert - xm) * (xvert - xp) + cf[3L] * (xvert - xm) * (xvert - x0)
  c(xvert, qv, qcoef)
}

.qd_qdelta <- function(b) {
  xv <- b$parent$xvalues
  fv <- b$parent$fvalues
  xm <- xv[1L]; x0 <- xv[2L]; xp <- xv[3L]
  fbox <- fv[b$parent_cindex]
  cf <- .qd_lagrangecoefs(xv, fv)
  cm <- cf[1L]; c0 <- cf[2L]; cp <- cf[3L]
  qvalue <- function(x) cm * (x - x0) * (x - xp) + c0 * (x - xm) * (x - xp) + cp * (x - xm) * (x - x0)
  bb <- .qd_bb1(b)
  qcoef <- cm + c0 + cp
  xvert <- (cm * (x0 + xp) + c0 * (xm + xp) + cp * (xm + x0)) / (2 * qcoef)
  if (isTRUE(qcoef > 0 && bb[1L] <= xvert && xvert <= bb[2L])) return(qvalue(xvert) - fbox)
  if (is.infinite(bb[1L]) || is.infinite(bb[2L])) {
    lcoef <- -cm * (x0 + xp) - c0 * (xm + xp) - cp * (xm + x0)
    ccoef <- cm * x0 * xp + c0 * xm * xp + cp * xm * x0
    if (isTRUE(qcoef == 0)) return(min(lcoef * bb[1L] + ccoef, lcoef * bb[2L] + ccoef) - fbox)
    qinf <- function(x) if (is.infinite(x)) abs(x) * sign(qcoef) else qvalue(x)
    return(min(qinf(bb[1L]), qinf(bb[2L])) - fbox)
  }
  min(qvalue(bb[1L]), qvalue(bb[2L])) - fbox
}

.qd_qdelta_sd <- function(b, sd) {
  p <- .qd_fpws(b, sd)
  if (p$parent$splitdim == sd) .qd_qdelta(p) else 0
}

.qd_is_diag_convex <- function(b) {
  for (i in seq_len(b$n)) {
    p <- .qd_fpws(b, i)
    if (.qd_isroot(p)) return(FALSE)
    if (!isTRUE(.qd_qfit(p$parent$xvalues, p$parent$fvalues)[3L] > 0)) return(FALSE)
  }
  TRUE
}

# ---------------------------------------------------------------------------
# Small utilities

.qd_lohi <- function(x, y) if (x <= y) c(x, y) else c(y, x)

.qd_lohi3 <- function(x, y, z) {
  if (z <= x) c(z, x, y) else if (z <= y) c(x, z, y) else c(x, y, z)
}

# order_pairs of three (x, f) pairs: a list of the sorted x and f
.qd_order_pairs3 <- function(p1, p2, p3) {
  op <- function(a, b) if (a[1L] <= b[1L]) list(a, b) else list(b, a)
  r <- op(p1, p2); p1 <- r[[1L]]; p2 <- r[[2L]]
  r <- op(p2, p3); p2 <- r[[1L]]; p3 <- r[[2L]]
  r <- op(p1, p2); p1 <- r[[1L]]; p2 <- r[[2L]]
  list(x = c(p1[1L], p2[1L], p3[1L]), f = c(p1[2L], p2[2L], p3[2L]))
}

# biggest_interval of QuadDIRECT, which returns the smallest of the three
.qd_biggest_interval <- function(a, b, c, d) {
  ab <- b - a; bc <- c - b; cd <- d - c
  if (ab <= bc && ab <= cd) c(a, b) else if (bc <= ab && bc <= cd) c(b, c) else c(c, d)
}

.qd_ensure_distinct <- function(x, x1, x2, bb, minfrac = 0.1) {
  lh <- .qd_lohi(x1, x2)
  x1 <- lh[1L]; x2 <- lh[2L]
  dxmin <- minfrac * min(x2 - x1, if (bb[2L] == x2) Inf else bb[2L] - x2, if (bb[1L] == x1) Inf else x1 - bb[1L])
  if (!isTRUE(dxmin > 0)) stop("AssertionError: dxmin > 0", call. = FALSE)
  if (abs(x - x1) < dxmin) {
    s <- if (x == x1) 1 else sign(x - x1)
    x <- max(bb[1L], x1 + dxmin * s)
  } else if (abs(x - x2) < dxmin) {
    s <- if (x == x2) -1 else sign(x - x2)
    x <- min(bb[2L], x2 + dxmin * s)
  }
  x
}

.qd_pick3 <- function(a, b, bb) {
  lh <- .qd_lohi(a, b)
  a <- lh[1L]; b <- lh[2L]
  iv <- .qd_biggest_interval(bb[1L], a, b, bb[2L])
  if (is.infinite(iv[1L])) return(c(a - 2 * (b - a), a, b))
  if (is.infinite(iv[2L])) return(c(a, b, b + 2 * (b - a)))
  .qd_lohi3(a, b, (iv[1L] + iv[2L]) / 2)
}

.qd_issame <- function(x1, x2, scale) isTRUE(all(abs(x1 - x2) < 1.4901161193847656e-08 * scale))

.qd_jlvec <- function(x) {
  num <- function(v) {
    if (is.nan(v)) return("NaN")
    if (is.infinite(v)) return(if (v > 0) "Inf" else "-Inf")
    for (d in 1:17) {
      s <- trimws(formatC(v, digits = d, format = "g"))
      if (as.numeric(s) == v) break
    }
    if (!grepl("[.e]", s)) s <- paste0(s, ".0")
    s
  }
  paste0("[", paste(vapply(x, num, ""), collapse = ", "), "]")
}

# ---------------------------------------------------------------------------
# Minimum edge lists (MELink): a linked list headed by a dummy item; the
# last item points to itself

.qd_mel <- function(l, w, f, nxt = NULL) {
  m <- new.env(parent = emptyenv())
  m$l <- l
  m$w <- w
  m$f <- f
  m$nxt <- if (is.null(nxt)) m else nxt
  m
}

.qd_mel_dropnext <- function(prev, nxt) {
  if (!identical(nxt, nxt$nxt)) {
    nxt <- nxt$nxt
    prev$nxt <- nxt
  } else {
    prev$nxt <- prev
    nxt <- prev
  }
  nxt
}

.qd_mel_trim <- function(mel, w, f) {
  prev <- mel
  nxt <- mel$nxt
  while (!identical(prev, nxt) && isTRUE(w > nxt$w)) {
    if (isTRUE(f <= nxt$f)) {
      nxt <- .qd_mel_dropnext(prev, nxt)
    } else {
      prev <- nxt
      nxt <- nxt$nxt
    }
  }
  if (isTRUE(w == nxt$w && f < nxt$f)) nxt <- .qd_mel_dropnext(prev, nxt)
  list(prev, nxt)
}

.qd_mel_insert <- function(mel, w, l, f) {
  r <- .qd_mel_trim(mel, w, f)
  prev <- r[[1L]]
  nxt <- r[[2L]]
  if (identical(prev, nxt)) {
    prev$nxt <- .qd_mel(l, w, f)
  } else if (isTRUE(f < nxt$f)) {
    prev$nxt <- .qd_mel(l, w, f, nxt)
  }
  invisible(mel)
}

.qd_mel_pop <- function(mel) {
  item <- mel$nxt
  mel$nxt <- if (identical(item$nxt, item)) mel else item$nxt
  item
}

.qd_mel_length <- function(mel) {
  k <- 0L
  s <- mel
  while (!identical(s, s$nxt)) {
    k <- k + 1L
    s <- s$nxt
  }
  k
}

# ---------------------------------------------------------------------------
# Splitting

# split!: evaluates the three points (reusing fcur at xcur), adds the
# children and returns the child of smallest value and xtmp at its point
.qd_split <- function(st, b, xtmp, sd, xsplit, lower, upper, xcur, fcur) {
  fsplit <- c(Inf, Inf, Inf)
  fmin <- Inf
  idxmin <- 0L
  for (l in 1:3) {
    if (!is.finite(xsplit[l])) stop("AssertionError: isfinite(xsplit[l])", call. = FALSE)
    if (xsplit[l] == xcur && !is.nan(fcur)) {
      ftmp <- fcur
    } else {
      xtmp[sd] <- xsplit[l]
      ftmp <- st$f(xtmp)
    }
    if (isTRUE(ftmp < fmin)) {
      fmin <- ftmp
      idxmin <- l
    }
    fsplit[l] <- ftmp
  }
  if (idxmin == 0L) stop("function was not finite at any evaluation point ", .qd_jlvec(xsplit), call. = FALSE)
  if (idxmin == 1L && isTRUE(fsplit[1L] == fsplit[2L])) idxmin <- 2L
  .qd_add_children(b, sd, xsplit, fsplit, lower, upper)
  xtmp[sd] <- xsplit[idxmin]
  list(box = b$children[[idxmin]], xtmp = xtmp)
}

.qd_trimschedule <- function(mes, b, x0, lower, upper) {
  if (b$splitdim == 0L) return(invisible())
  for (child in b$children) {
    fval <- child$parent$fvalues[child$parent_cindex]
    for (i in seq_len(b$n)) .qd_mel_trim(mes[[i]], .qd_width(child, i, x0, lower, upper), fval)
  }
  invisible()
}

.qd_trisect <- function(st, b, xtmp, sd, bb, dx, xcur, fcur) {
  l <- bb[1L]; r <- bb[2L]
  if (is.infinite(l)) l <- r - 6 * dx else if (is.infinite(r)) r <- l + 6 * dx
  a <- 5 * l / 6 + r / 6
  m <- (l + r) / 2
  cc <- l / 6 + 5 * r / 6
  if (xcur < (a + m) / 2) a <- xcur else if (xcur < (m + cc) / 2) m <- xcur else cc <- xcur
  if (a < m && m < cc) .qd_split(st, b, xtmp, sd, c(a, m, cc), bb[1L], bb[2L], xcur, fcur)
  invisible()
}

.qd_autosplit <- function(st, b, mes, x0, xtmp, sd, splits, lower, upper, minwidth, visited) {
  if (.qd_in(b, visited$items)) stop("already visited box", call. = FALSE)
  if (b$qnconverged) return(list(b, .qd_value(b)))
  if (b$splitdim != 0L) {
    b <- .qd_find_smallest_child_leaf(b)
    xtmp <- .qd_position(b, x0)
  }
  n <- b$n
  if (sd == 0L) {
    ns <- .qd_count_splits(b)
    sd <- which.min(ns)
    if (ns[sd] > 0L) {
      bbs <- .qd_bb_fill(cbind(lower, upper), b)
      if (all(is.finite(bbs[sd, ])))
        for (i in seq_len(n)) if (!all(is.finite(bbs[i, ]))) {
          sd <- i
          break
        }
    }
  }
  lwr <- lower[sd]; upr <- upper[sd]
  p <- .qd_fpws(b, sd)
  if (.qd_isroot(p) && p$splitdim != sd) {
    .qd_split(st, b, xtmp, sd, splits[[sd]], lwr, upr, x0[sd], b$parent$fvalues[b$parent_cindex])
    .qd_trimschedule(mes, b, x0, lower, upper)
    return(list(b, min(b$fvalues)))
  }
  bb <- .qd_bb1(p)
  if (!isTRUE(bb[2L] - bb[1L] > max(minwidth[sd], .qd_epswidth(bb)))) return(list(b, .qd_value(b)))
  xp <- p$parent$xvalues
  fp <- p$parent$fvalues
  dx <- max(xp[2L] - xp[1L], xp[3L] - xp[2L])
  xcur <- xp[p$parent_cindex]
  if (xcur != xtmp[sd]) stop("AssertionError: xcur == xtmp[splitdim]", call. = FALSE)
  fcur <- b$parent$fvalues[b$parent_cindex]
  qcoef <- .qd_qfit(xp, fp)[3L]
  if (isTRUE(qcoef > 0)) {
    xnew <- if (is.infinite(bb[1L])) xcur - 8 * (bb[2L] - xcur) else if (is.infinite(bb[2L])) xcur + 8 * (xcur - bb[1L])
            else if (bb[2L] - xcur > xcur - bb[1L]) (xcur + 2 * bb[2L]) / 3 else (2 * bb[1L] + xcur) / 3
    if (!is.finite(xnew)) stop("AssertionError: isfinite(xnew)", call. = FALSE)
    xtmp[sd] <- xnew
    fnew <- st$f(xtmp)
    xvert <- (xcur + xnew) / 2 - (fnew - fcur) / (2 * qcoef * (xnew - xcur))
    if (isTRUE(bb[1L] <= xvert && xvert <= bb[2L])) {
      xvert <- .qd_ensure_distinct(xvert, xcur, xnew, bb)
      xtmp[sd] <- xvert
      fvert <- st$f(xtmp)
      o <- .qd_order_pairs3(c(xcur, fcur), c(xnew, fnew), c(xvert, fvert))
      .qd_add_children(b, sd, o$x, o$f, lwr, upr)
      .qd_trimschedule(mes, b, x0, lower, upper)
      return(list(b, min(b$fvalues)))
    }
    if (isTRUE(xvert < bb[1L])) {
      xtmp[sd] <- bb[1L]
      dir <- -1
    } else {
      xtmp[sd] <- bb[2L]
      dir <- 1
    }
    r <- .qd_find_leaf_at_edge(.qd_get_root(b), xtmp, sd, dir)
    nbr <- r[[1L]]
    success <- r[[2L]]
    xmid <- (xcur + xnew) / 2
    xtmp[sd] <- xmid
    fmid <- st$f(xtmp)
    o <- .qd_order_pairs3(c(xcur, fcur), c(xnew, fnew), c(xmid, fmid))
    .qd_add_children(b, sd, o$x, o$f, lwr, upr)
    .qd_trimschedule(mes, b, x0, lower, upper)
    minbox <- min(b$fvalues)
    if (!success || .qd_in(nbr, visited$items) || !is.finite(.qd_value(nbr)) || length(visited$items) > n)
      return(list(b, minbox))
    xn <- .qd_position(nbr, x0)
    visited$items <- .qd_set_add(visited$items, b)
    r <- .qd_autosplit(st, nbr, mes, x0, xn, 0L, splits, lower, upper, minwidth, visited)
    if (isTRUE(r[[2L]] < minbox)) return(r)
    return(list(b, minbox))
  }
  .qd_trisect(st, b, xtmp, sd, bb, dx, xcur, fcur)
  .qd_trimschedule(mes, b, x0, lower, upper)
  list(b, if (b$splitdim == 0L) .qd_value(b) else min(b$fvalues))
}

.qd_minimum_edges <- function(root, x0, lower, upper, minwidth, extrapolate = TRUE) {
  n <- root$n
  mes <- lapply(seq_len(n), function(i) .qd_mel(root, -Inf, -Inf))
  nthresh <- ceiling((((n + 1) * (n + 2)) %/% 2) / n)
  for (b in .qd_leaves(root)) {
    if (b$qnconverged) next
    fval <- .qd_value(b)
    if (!is.finite(fval)) next
    ns <- .qd_count_splits(b)
    nmin <- min(ns)
    for (i in seq_len(n)) {
      if (nmin < nthresh && ns[i] != nmin) next
      bb <- .qd_bb(.qd_fpws(b, i), lower[i], upper[i])
      if (isTRUE(bb[2L] - bb[1L] < minwidth[i])) next
      .qd_mel_insert(mes[[i]], .qd_width(b, i, x0, lower, upper), b,
                     if (extrapolate) fval + .qd_qdelta_sd(b, i) else fval)
    }
  }
  mes
}

# sweep!: the boxes split, as a list (QuadDIRECT's Set)
.qd_sweep <- function(st, root, x0, splits, lower, upper, extrapolate, fvalue, minwidth) {
  mes <- .qd_minimum_edges(root, x0, lower, upper, minwidth, extrapolate)
  visited <- new.env(parent = emptyenv())
  splitboxes <- list()
  for (i in order(vapply(mes, .qd_mel_length, 0L))) {
    me <- mes[[i]]
    while (!identical(me$nxt, me)) {
      b <- .qd_mel_pop(me)$l
      if (b$qnconverged) next
      bb <- .qd_bb(.qd_fpws(b, i), lower[i], upper[i])
      if (isTRUE(bb[2L] - bb[1L] < minwidth[i])) next
      xtmp <- .qd_position(b, x0)
      ns <- .qd_count_splits(b)
      if (ns[i] > 0L && any(ns == 0L)) next
      visited$items <- list()
      r <- .qd_autosplit(st, b, mes, x0, xtmp, i, splits, lower, upper, minwidth, visited)
      splitboxes <- .qd_set_add(splitboxes, r[[1L]])
      if (isTRUE(r[[2L]] <= fvalue)) return(.qd_unique_smallest_leaves(splitboxes))
    }
  }
  .qd_unique_smallest_leaves(splitboxes)
}

# ---------------------------------------------------------------------------
# Quadratic models by incremental Gaussian elimination (QmIGE)

.qd_qmige <- function(n) {
  m <- ((n + 1L) * (n + 2L)) %/% 2L - 1L
  Q <- new.env(parent = emptyenv())
  Q$n <- n
  Q$coefs <- matrix(0, m, m)
  Q$rhs <- numeric(m)
  Q$dimpiv <- integer(n)
  Q$rowzero <- rep(TRUE, m)
  Q$rowbox <- vector("list", m)
  Q$ndims <- 0L
  Q$nzrows <- 0L
  Q$rowtmp <- numeric(m)
  Q
}

.qd_setrow <- function(Q, dx, ndims, sd) {
  rt <- numeric(length(Q$rowtmp))
  k <- 0L
  have <- FALSE
  for (j in seq_len(ndims)) {
    d <- Q$dimpiv[j]
    have <- have || sd == d
    dxd <- dx[d]
    if (isTRUE(dxd == 0)) {
      k <- k + j + 1L
      next
    }
    k <- k + 1L
    rt[k] <- dxd
    for (i in seq_len(j - 1L)) {
      k <- k + 1L
      rt[k] <- dxd * dx[Q$dimpiv[i]]
    }
    k <- k + 1L
    rt[k] <- (dxd * dxd) / 2
  }
  if (!have) {
    ndims <- ndims + 1L
    Q$dimpiv[ndims] <- sd
    dxd <- dx[sd]
    k <- k + 1L
    rt[k] <- dxd
    for (i in seq_len(ndims - 1L)) {
      k <- k + 1L
      rt[k] <- dxd * dx[Q$dimpiv[i]]
    }
    k <- k + 1L
    rt[k] <- (dxd * dxd) / 2
  }
  while (ndims > 0L && isTRUE(dx[Q$dimpiv[ndims]] == 0)) ndims <- ndims - 1L
  Q$rowtmp <- rt
  ndims
}

# incremental_gaussian_elimination! on Q$rowtmp over the rows rng
.qd_ige <- function(Q, rng, newrhs = NULL, meta = NULL, canswap = TRUE) {
  rtol <- 1000 * 2.220446049250313e-16
  nr <- Q$rowtmp
  co <- Q$coefs
  lo <- rng[1L]
  i <- rng[length(rng)]
  while (i >= lo) {
    if (isTRUE(nr[i] == 0)) {
      i <- i - 1L
      next
    }
    if (Q$rowzero[i]) break
    if (canswap && isTRUE(abs(nr[i]) > abs(co[i, i]))) {
      js <- seq_len(i)
      tmp <- nr[js]
      nr[js] <- co[i, js]
      co[i, js] <- tmp
      if (!is.null(newrhs)) {
        tmp <- Q$rhs[i]
        Q$rhs[i] <- newrhs
        newrhs <- tmp
      }
      if (!is.null(meta)) {
        tmp <- Q$rowbox[[i]]
        Q$rowbox[i] <- list(meta)
        meta <- tmp
      }
    }
    cc <- nr[i] / co[i, i]
    if (i > 1L) {
      js <- seq_len(i - 1L)
      a <- nr[js]
      sub <- cc * co[i, js]
      same <- (a == sub) | (abs(a - sub) < rtol * (abs(a) + abs(sub)))
      same[is.na(same)] <- FALSE
      nr[js] <- ifelse(same, 0, a - sub)
    }
    nr[i] <- 0
    if (!is.null(newrhs)) newrhs <- newrhs - cc * Q$rhs[i]
    i <- i - 1L
  }
  if (i < lo) i <- 0L
  Q$rowtmp <- nr
  Q$coefs <- co
  list(i = i, newrhs = newrhs, meta = meta)
}

.qd_insert <- function(Q, dx, df, sd, b, canswap = TRUE) {
  ndims_old <- Q$ndims
  ndims <- .qd_setrow(Q, dx, ndims_old, sd)
  if (ndims > ndims_old) {
    i <- ndims + (ndims * (ndims + 1L)) %/% 2L
    Q$coefs[i, seq_len(i)] <- Q$rowtmp[seq_len(i)]
    Q$rowzero[i] <- FALSE
    Q$rowbox[[i]] <- b
    Q$rhs[i] <- df
    Q$ndims <- ndims
    Q$nzrows <- Q$nzrows + 1L
    return(i)
  }
  if (ndims == 0L) return(0L)
  bs <- ndims + ((ndims - 1L) * ndims) %/% 2L
  r <- .qd_ige(Q, bs:(bs + ndims), df, b, canswap)
  fr <- r$i
  if (fr > 0L) {
    Q$rowzero[fr] <- FALSE
    Q$coefs[fr, seq_len(fr)] <- Q$rowtmp[seq_len(fr)]
    Q$rhs[fr] <- r$newrhs
    Q$rowbox[fr] <- list(r$meta)
    Q$nzrows <- Q$nzrows + 1L
  }
  fr
}

# solve(Q): the gradient g and Hessian B of the model, by the forward
# substitution of LinearAlgebra's generic triangular solve
.qd_qsolve <- function(Q) {
  A <- Q$coefs
  y <- Q$rhs
  m <- length(y)
  y[1L] <- y[1L] / A[1L, 1L]
  if (m > 1L) {
    ii <- 2:m
    y[ii] <- y[ii] - A[ii, 1L] * y[1L]
    for (j in 2:m) {
      y[j] <- y[j] / A[j, j]
      if (j < m) {
        ii <- (j + 1L):m
        y[ii] <- y[ii] - A[ii, j] * y[j]
      }
    }
  }
  n <- Q$n
  dp <- Q$dimpiv
  g <- numeric(n)
  B <- matrix(0, n, n)
  k <- 0L
  for (i in seq_len(n)) {
    k <- k + 1L
    g[dp[i]] <- y[k]
    for (j in seq_len(i)) {
      k <- k + 1L
      B[dp[i], dp[j]] <- B[dp[j], dp[i]] <- y[k]
    }
  }
  list(g = g, B = B)
}

.qd_descend <- function(Q, c, b, x0, xbase, scale, skip = FALSE) {
  if (b$qnconverged) return(FALSE)
  Q$bbs <- .qd_bb_fill(Q$bbs, b)
  for (i in seq_len(b$n)) {
    xb <- xbase[i]
    if (isTRUE(min(Q$bbs[i, 1L] - xb, xb - Q$bbs[i, 2L]) > 100 * scale[i])) return(TRUE)
  }
  x <- .qd_position(b, x0)
  thisx <- if (b$splitdim == 0L) 0 else x[b$splitdim]
  if (!skip) {
    leaf <- .qd_find_leaf_at(b, x)
    dx <- (x - xbase) / scale
    if (isTRUE(all(abs(dx) < 100))) {
      .qd_insert(Q, dx, .qd_value(b) - c, b$parent$splitdim, leaf)
    } else {
      sd <- b$parent$splitdim
      if (!(sd %in% Q$dimpiv)) {
        Q$ndims <- Q$ndims + 1L
        Q$dimpiv[Q$ndims] <- sd
      }
    }
  }
  if (b$splitdim != 0L) {
    iskip <- match(thisx, b$xvalues)
    for (i in 1:3) if (!.qd_descend(Q, c, b$children[[i]], x0, xbase, scale, isTRUE(i == iskip))) return(FALSE)
  }
  TRUE
}

.qd_build_quadratic_model <- function(b, x0, scale, thresh) {
  n <- b$n
  Q <- .qd_qmige(n)
  c <- .qd_value(b)
  xbase <- .qd_position(b, x0)
  Q$bbs <- .qd_bb_fill(matrix(NaN, n, 2L), b)
  m <- length(Q$rhs)
  ok <- TRUE
  need_extra <- TRUE
  while (ok && (Q$nzrows < m || need_extra) && !.qd_isroot(b)) {
    if (Q$nzrows == m) need_extra <- isTRUE(min(abs(diag(Q$coefs))) < thresh)
    cindex <- b$parent_cindex
    b <- b$parent
    ok <- all(is.finite(b$fvalues)) && !b$qnconverged
    if (!ok) break
    j <- 1L
    if (j == cindex) j <- j + 1L
    ok <- .qd_descend(Q, c, b$children[[j]], x0, xbase, scale)
    if (!ok) break
    j <- j + 1L
    if (j == cindex) j <- j + 1L
    ok <- .qd_descend(Q, c, b$children[[j]], x0, xbase, scale)
  }
  list(Q = Q, xbase = xbase, c = c, ok = ok)
}

.qd_split_insert <- function(st, Q, rbox, p, rng, xbase, c, xtmp, sd, scale, splits, lower, upper) {
  row <- rng[1L]
  n <- rbox$n
  xcur <- xtmp[sd]
  if (.qd_isroot(p)) {
    xsplit <- splits[[sd]]
    lwr <- lower[sd]; upr <- upper[sd]
  } else {
    bb <- .qd_bb1(p)
    lwr <- bb[1L]; upr <- bb[2L]
    ssd <- scale[sd]
    lo <- max(lwr, 5 * lwr / 6 + xcur / 6, xcur - ssd)
    hi <- min(upr, xcur / 6 + 5 * upr / 6, xcur + ssd)
    xm <- xcur
    if (lo + 1.4901161193847656e-08 * ssd >= xcur) {
      lo <- xcur
      xm <- (xcur + hi) / 2
    } else if (xcur + 1.4901161193847656e-08 * ssd >= hi) {
      xm <- (xcur + lo) / 2
      hi <- xcur
    }
    xsplit <- c(lo, xm, hi)
  }
  if (!isTRUE(xsplit[1L] < xsplit[2L] && xsplit[2L] < xsplit[3L]))
    stop("AssertionError: xsplit[1] < xsplit[2] < xsplit[3]", call. = FALSE)
  j1 <- 1L
  if (xsplit[j1] == xcur) j1 <- j1 + 1L
  j2 <- j1 + 1L
  if (xsplit[j2] == xcur) j2 <- j2 + 1L
  l <- xsplit[j1]
  h <- xsplit[j2]
  Qd <- abs(Q$coefs[row, row])
  if (Qd != 0) {
    xtmp[sd] <- l
    .qd_setrow(Q, (xtmp - xbase) / scale, n, sd)
    .qd_ige(Q, rng, canswap = FALSE)
    if (isTRUE(abs(Q$rowtmp[row]) < Qd)) {
      xtmp[sd] <- h
      .qd_setrow(Q, (xtmp - xbase) / scale, n, sd)
      .qd_ige(Q, rng, canswap = FALSE)
      if (isTRUE(abs(Q$rowtmp[row]) < Qd)) return(list(ins = c(row, -1L), box = rbox))
    }
  }
  val <- .qd_value(rbox)
  fs <- c(val, val, val)
  xtmp[sd] <- l
  fs[j1] <- st$f(xtmp)
  xtmp[sd] <- h
  fs[j2] <- st$f(xtmp)
  .qd_add_children(rbox, sd, xsplit, fs, lwr, upr)
  ins <- c(0L, 0L)
  if (l != xbase[sd]) {
    xtmp[sd] <- l
    ins[1L] <- .qd_insert(Q, (xtmp - xbase) / scale, rbox$fvalues[j1] - c, sd, rbox$children[[j1]], canswap = FALSE)
  }
  if (h != xbase[sd]) {
    xtmp[sd] <- h
    ins[2L] <- .qd_insert(Q, (xtmp - xbase) / scale, rbox$fvalues[j2] - c, sd, rbox$children[[j2]], canswap = FALSE)
  }
  j <- 1L
  if (j == j1) j <- j + 1L
  if (j == j2) j <- j + 1L
  list(ins = ins, box = rbox$children[[j]])
}

.qd_complete_quadratic_model <- function(st, Q, c, b, x0, xbase, scale, splits, lower, upper, thresh) {
  n <- b$n
  rowoffset <- 0L
  for (idim in seq_len(n)) {
    dg <- rowoffset + seq_len(idim + 1L)
    if (isTRUE(all(abs(Q$coefs[cbind(dg, dg)]) >= thresh))) {
      rowoffset <- rowoffset + idim + 1L
      next
    }
    rowb <- row <- rowoffset + idim + 1L
    irow <- idim
    if (is.null(Q$rowbox[[rowb]])) {
      sd <- Q$dimpiv[idim]
      xtmp <- .qd_position(b, x0)
      xcur <- xbase[sd]
      s <- scale[sd]
      bb <- .qd_bb1(.qd_fpws(b, sd))
      l <- max(bb[1L], xcur - s / 2)
      u <- min(bb[2L], xcur + s / 2)
      xm <- xcur
      if (l >= xcur) {
        l <- xcur
        xm <- (xcur + u) / 2
      } else if (xcur >= u) {
        u <- xcur
        xm <- (l + xcur) / 2
      }
      if (abs(u - l) < .qd_epswidth(c(l, u))) return(invisible(Q))
      # QuadDIRECT passes the new middle point as the point of the known
      # value, which gives the box value to the wrong child when the box
      # point is an end of the split; the box point is passed instead
      .qd_split(st, b, xtmp, sd, c(l, xm, u), lower[sd], upper[sd], xcur, .qd_value(b))
      j <- 3L
      j1 <- 1L
      if (b$xvalues[j1] == xcur) {
        j <- j1
        j1 <- j1 + 1L
      }
      j2 <- j1 + 1L
      if (b$xvalues[j2] == xcur) {
        j <- j2
        j2 <- j2 + 1L
      }
      Q$rowbox[[rowb]] <- if (isTRUE(b$fvalues[j1] < b$fvalues[j2])) b$children[[j1]] else b$children[[j2]]
      b <- b$children[[j]]
    }
    while (irow >= 0L) {
      Qd <- abs(Q$coefs[row, row])
      if (isTRUE(Qd < thresh)) {
        rbox <- Q$rowbox[[rowb]]
        if (rbox$splitdim != 0L) stop("AssertionError: isleaf(rbox)", call. = FALSE)
        xtmp <- .qd_position(rbox, x0)
        dx <- xtmp - xbase
        dp <- Q$dimpiv
        if (irow > 0L) {
          sd <- dp[irow]
        } else if (dx[dp[idim]] == 0) {
          sd <- dp[idim]
        } else {
          dxt <- dx / scale
          maxdiag <- 0
          sd <- 0L
          rng <- (rowoffset + 1L):(rowoffset + idim + 1L)
          Q$rowzero[row] <- TRUE
          for (ii in seq_len(idim)) {
            s2 <- dp[ii]
            xold <- dxt[s2]
            dxt[s2] <- if (dxt[s2] == 0) 1 else 2 * dxt[s2]
            .qd_setrow(Q, dxt, n, s2)
            dxt[s2] <- xold
            .qd_ige(Q, rng, canswap = FALSE)
            thisdiag <- abs(Q$rowtmp[row])
            if (isTRUE(thisdiag > maxdiag)) {
              maxdiag <- thisdiag
              sd <- s2
            }
          }
          if (sd == 0L) stop("no direction works", call. = FALSE)
        }
        p <- .qd_fpws(rbox, sd)
        if (!.qd_isroot(p)) {
          bb <- .qd_bb1(p)
          be <- .qd_epswidth(bb)
          if (!isTRUE(bb[2L] - bb[1L] > be) || !isTRUE(bb[2L] - xtmp[sd] > be) || !isTRUE(xtmp[sd] - bb[1L] > be))
            return(invisible(Q))
        }
        r <- .qd_split_insert(st, Q, rbox, p, row:(rowoffset + idim + 1L), xbase, c, xtmp, sd, scale, splits,
                              lower, upper)
        if (Qd == 0 && max(r$ins) != row) stop("AssertionError: maximum(insrows) == row", call. = FALSE)
        if (row == rowb) Q$rowbox[[row]] <- r$box
      }
      rb <- Q$rowbox[[row]]
      if (!is.null(rb) && rb$splitdim == 0L) rowb <- row
      row <- row - 1L
      irow <- irow - 1L
    }
    rowoffset <- rowoffset + idim + 1L
  }
  invisible(Q)
}

# ---------------------------------------------------------------------------
# The quasi Newton step

# lls_bounded!: the box bounded minimum of x' H x / 2 + f' x by coordinate
# descent
.qd_lls_bounded <- function(x, H, f, lower, upper, tol = 1.4901161193847656e-08) {
  n <- length(x)
  niter <- 0L
  while (niter <= 1000L) {
    g <- vapply(seq_len(n), function(i) .jl_dot(H[i, ], x), 0) + f
    xchange <- 0
    for (k in seq_len(n)) {
      xold <- x[k]
      xnew <- .qd_clamp(xold - g[k] / H[k, k], lower[k], upper[k])
      x[k] <- xnew
      dx <- xnew - xold
      xchange <- max(xchange, abs(dx))
      if (!isTRUE(dx == 0)) g <- g + dx * H[, k]
    }
    if (isTRUE(xchange < tol)) break
    niter <- niter + 1L
  }
  x
}

.qd_quasinewton <- function(st, b, B, g, c, scale, x0, lower, upper, itermax = 20L) {
  x <- .qd_position(b, x0)
  fbox <- fmin <- .qd_value(b)
  fail <- function(leaf = b) list(leaf = leaf, fmin = fmin, alpha = 0)
  if (!is.finite(fbox)) return(fail())
  L <- .qd_poschol(B)
  dx <- -.qd_cholsolve(L, g)
  lsc <- (lower - x) / scale
  usc <- (upper - x) / scale
  dx <- .qd_lls_bounded(.qd_clamp(dx, lsc, usc), .qd_llt(L), g, lsc, usc)
  if (isTRUE(all(abs(dx) < sqrt(vapply(scale, .qd_eps, 0))))) {
    b$qnconverged <- TRUE
    return(fail())
  }
  dx <- dx * scale
  alpha <- 1
  inside <- function(y) isTRUE(all(lower <= y & y <= upper))
  iter <- 0L
  while (!inside(x + alpha * dx) && iter < itermax) {
    alpha <- alpha / 2
    iter <- iter + 1L
  }
  if (iter == itermax) return(fail())
  iter <- 0L
  root <- .qd_get_root(b)
  while (iter < itermax) {
    iter <- iter + 1L
    xtarget <- x + alpha * dx
    ftarget <- st$f(xtarget)
    if (isTRUE(ftarget < fbox)) break
    alpha <- alpha / 2
  }
  if (isTRUE(ftarget >= fbox)) return(fail())
  leaf <- .qd_find_leaf_at(root, xtarget)
  if (!is.finite(.qd_value(leaf))) return(fail())
  if (isTRUE(ftarget >= .qd_value(leaf))) return(fail())
  n <- b$n
  targeted <- logical(n)
  # the model; c is reassigned below, as the closure of QuadDIRECT sees
  q <- function(y) .qd_quadform(y, B) / 2 + .jl_dot(g, y) + c
  for (j in seq_len(n)) {
    xleaf <- .qd_position(leaf, x0)
    xtest <- xleaf
    imin <- 0L
    qmin <- Inf
    for (i in seq_len(n)) {
      if (targeted[i]) next
      xtest[i] <- xtarget[i]
      qx <- q(xtest)
      if (isTRUE(qx < qmin)) {
        qmin <- qx
        imin <- i
      }
    }
    if (imin == 0L) break
    xcur <- xleaf[imin]
    xt <- xtarget[imin]
    if (isTRUE(abs(xcur - xt) < sqrt(.qd_eps(scale[imin])))) {
      targeted[imin] <- TRUE
      next
    }
    bb <- .qd_bb(.qd_fpws(leaf, imin), lower[imin], upper[imin])
    abc <- .qd_pick3(xcur, xt, bb)
    a <- abc[1L]; m <- abc[2L]; c <- abc[3L]
    minsep <- 2.220446049250313e-16 * (c - a)
    if (isTRUE(a + minsep >= m) || isTRUE(m + minsep >= c)) {
      targeted[imin] <- TRUE
      next
    }
    .qd_split(st, leaf, xleaf, imin, abc, bb[1L], bb[2L], xleaf[imin], .qd_value(leaf))
    fmin <- min(fmin, leaf$fvalues)
    leaf <- leaf$children[[match(xt, leaf$xvalues)]]
    if (!is.finite(.qd_value(leaf))) return(fail(leaf))
    targeted[imin] <- TRUE
  }
  list(leaf = leaf, fmin = fmin, alpha = alpha)
}

# ---------------------------------------------------------------------------
# analyze and analyze!

.qd_analyze <- function(f, splits, lower, upper, rtol = 1e-3, atol = 0, fvalue = -Inf, maxevals = 2500,
                        nquasinewton = NULL, minwidth = NULL) {
  n <- length(lower)
  if (length(upper) != n || length(splits) != n)
    stop("inconsistent dimensionality among lower, upper, and splits", call. = FALSE)
  for (i in seq_len(n)) {
    s <- splits[[i]]
    if (length(s) != 3L)
      stop("3 values along each dimension required, got ", .qd_jlvec(s), " along dimension ", i, call. = FALSE)
    if (!isTRUE(lower[i] <= s[1L] && s[1L] < s[2L] && s[2L] < s[3L] && s[3L] <= upper[i]))
      stop("splits must be in strictly increasing order and lie between lower and upper, got ", .qd_jlvec(s),
           " along dimension ", i, call. = FALSE)
    if (!all(is.finite(s))) stop("not all splits along dimension ", i, " are finite", call. = FALSE)
  }
  st <- new.env(parent = emptyenv())
  st$nevals <- 0L
  st$f <- function(x) {
    st$nevals <- st$nevals + 1L
    v <- f(x)
    if (is.na(v)) NaN else v
  }
  # init
  x0 <- vapply(splits, function(s) s[2L], 0)
  root <- b <- .qd_box(n = n)
  xtmp <- x0
  fcur <- NaN
  for (i in seq_len(n)) {
    r <- .qd_split(st, b, xtmp, i, splits[[i]], lower[i], upper[i], x0[i], fcur)
    b <- r$box
    xtmp <- r$xtmp
    fcur <- b$parent$fvalues[b$parent_cindex]
  }
  # analyze!, with the evaluations counted from here
  st$nevals <- 0L
  minwidth <- minwidth %||% numeric(n)
  nquasinewton <- nquasinewton %||% (((n + 1) * (n + 2)) %/% 2)
  box <- .qd_minimum(root)
  boxval <- .qd_value(box)
  tol_counter <- 0L
  baseline <- len <- lenold <- length(.qd_leaves(root))
  if (baseline == st$nevals) baseline <- 0L
  used_qn <- nquasinewton > maxevals
  thresh <- 1.4901161193847656e-08
  m <- ((n + 1L) * (n + 2L)) %/% 2L - 1L
  while (isTRUE(boxval > fvalue) && (tol_counter <= n || !used_qn) && len < maxevals) {
    lastval <- boxval
    splitboxes <- .qd_sweep(st, root, x0, splits, lower, upper, TRUE, fvalue, minwidth)
    box <- .qd_minimum(root)
    ne <- baseline + st$nevals
    if (nquasinewton < ne && ne < maxevals && isTRUE(boxval > fvalue) && length(splitboxes)) {
      spbxs <- splitboxes[order(vapply(splitboxes, .qd_value, 0))]
      basin <- Filter(.qd_is_diag_convex, spbxs)
      bad <- logical(length(basin))
      for (ibox in seq_along(basin)) {
        b <- basin[[ibox]]
        if (b$qnconverged || b$splitdim != 0L || bad[ibox]) next
        if (baseline + st$nevals >= maxevals) break
        ns <- .qd_count_splits(b)
        xtmp <- .qd_position(b, x0)
        for (i in seq_len(n)) if (ns[i] == 0L) {
          r <- .qd_split(st, b, xtmp, i, splits[[i]], lower[i], upper[i], xtmp[i], .qd_value(b))
          b <- r$box
          xtmp <- r$xtmp
        }
        scale <- .qd_boxscale(b, splits)
        bq <- .qd_build_quadratic_model(b, x0, scale, thresh)
        Q <- bq$Q
        if (bq$ok && all(Q$dimpiv > 0L)) {
          .qd_complete_quadratic_model(st, Q, bq$c, b, x0, bq$xbase, scale, splits, lower, upper, thresh)
          if (Q$nzrows == m && isTRUE(min(abs(diag(Q$coefs))) >= thresh)) {
            gB <- .qd_qsolve(Q)
            qn <- .qd_quasinewton(st, b, gB$B, gB$g, bq$c, scale, x0, lower, upper)
            used_qn <- used_qn || qn$alpha > 0
            if (isTRUE(qn$fmin < fvalue)) break
            if (qn$alpha > 0) {
              bv <- .qd_value(b)
              if (isTRUE(bv - qn$fmin <= max(atol, rtol * (abs(bv) + abs(qn$fmin)))) &&
                  .qd_issame(.qd_position(b, x0), .qd_position(qn$leaf, x0), scale))
                qn$leaf$qnconverged <- TRUE
              if (qn$alpha == 1) for (j in seq_along(basin)[-seq_len(ibox)]) {
                bj <- basin[[j]]
                dx <- (.qd_position(bj, x0) - bq$xbase) / scale
                pred <- bq$c + .jl_dot(gB$g, dx) + .qd_quadform(dx, gB$B) / 2
                if (isTRUE(abs(pred - .qd_value(bj)) <= 0.1 * abs(pred - qn$fmin))) bad[j] <- TRUE
              }
            }
          }
        }
      }
      box <- .qd_minimum(root)
      boxval <- .qd_value(box)
    }
    len <- baseline + st$nevals
    if (len == lenold) break
    lenold <- len
    if (!isTRUE(boxval <= lastval)) stop("AssertionError: boxval <= lastval", call. = FALSE)
    if (isTRUE(lastval - boxval < atol) || isTRUE(lastval - boxval < rtol * (abs(lastval) + abs(boxval)))) {
      tol_counter <- tol_counter + 1L
    } else {
      tol_counter <- 0L
    }
  }
  list(root = root, x0 = x0)
}
