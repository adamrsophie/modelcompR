# The LP presolve of HiGHS 1.12.0 (presolve/HPresolve.cpp,
# util/HighsLinearSumBounds.cpp, util/HighsSplay.h, util/HighsCDouble.h) and
# its postsolve (presolve/HighsPostsolveStack.cpp), as SciPy 1.18.1 vendors
# them. MIT license notice in R/optimization-scipy-highs-util.R. Indices are
# 1-based and 0 stands for HiGHS's -1; rows and columns keep their original
# indices throughout (shrinkProblem only rebuilds the matrix links).

# HighsCDouble: a value hi + lo kept as c(hi, lo)
.cd_two_sum <- function(a, b) {
  x <- a + b
  z <- x - a
  c(x, (a - (x - z)) + (b - z))
}
.cd_split <- function(a) {
  c <- 134217729 * a
  x <- c - (c - a)
  c(x, a - x)
}
.cd_two_product <- function(a, b) {
  x <- a * b
  sa <- .cd_split(a); sb <- .cd_split(b)
  c(x, sa[2L] * sb[2L] - (((x - sa[1L] * sb[1L]) - sa[2L] * sb[1L]) - sa[1L] * sb[2L]))
}
# CD + double
.cd_add <- function(a, v) {
  r <- .cd_two_sum(a[1L], v)
  c(r[1L], r[2L] + a[2L])
}
# CD += CD (operator+=)
.cd_iadd <- function(a, v) {
  r <- .cd_two_sum(v[1L], a[1L])
  c(r[1L], a[2L] + r[2L] + v[2L])
}
# CD * double
.cd_mul <- function(a, v) {
  r <- .cd_two_product(a[1L], v)
  .cd_iadd(r, c(a[2L] * v, 0))
}
# CD / double
.cd_div <- function(a, v) {
  d <- c(a[1L] / v, a[2L] / v)
  cc <- .cd_mul(d, v)
  cc <- .cd_add(cc, -a[1L]); cc[2L] <- cc[2L] - a[2L]
  cc <- cc / v
  r <- .cd_add(d, -cc[1L]); r[2L] <- r[2L] - cc[2L]
  r
}
.cd_val <- function(a) a[1L] + a[2L]

# top down splay of the row trees (key: column index)
.hs_splay <- function(pr, key, root) {
  if (root == 0L) return(0L)
  L <- pr$ARleft; R <- pr$ARright; K <- pr$Acol
  Nleft <- 0L; Nright <- 0L
  lr <- 0L; rl <- 0L
  repeat {
    if (key < K[root]) {
      left <- L[root]
      if (left == 0L) break
      if (key < K[left]) {
        y <- left
        L[root] <- R[y]; R[y] <- root; root <- y
        if (L[root] == 0L) break
      }
      if (rl == 0L) Nleft <- root else L[rl] <- root
      rl <- root
      root <- L[root]
    } else if (key > K[root]) {
      right <- R[root]
      if (right == 0L) break
      if (key > K[right]) {
        y <- right
        R[root] <- L[y]; L[y] <- root; root <- y
        if (R[root] == 0L) break
      }
      if (lr == 0L) Nright <- root else R[lr] <- root
      lr <- root
      root <- R[root]
    } else break
  }
  if (lr == 0L) Nright <- L[root] else R[lr] <- L[root]
  if (rl == 0L) Nleft <- R[root] else L[rl] <- R[root]
  L[root] <- Nright; R[root] <- Nleft
  pr$ARleft <- L; pr$ARright <- R
  root
}

.hs_splay_link <- function(pr, pos, row) {
  root <- pr$rowroot[row]
  if (root == 0L) {
    pr$ARleft[pos] <- 0L; pr$ARright[pos] <- 0L
  } else {
    root <- .hs_splay(pr, pr$Acol[pos], root)
    if (pr$Acol[pos] < pr$Acol[root]) {
      pr$ARleft[pos] <- pr$ARleft[root]; pr$ARright[pos] <- root; pr$ARleft[root] <- 0L
    } else {
      pr$ARright[pos] <- pr$ARright[root]; pr$ARleft[pos] <- root; pr$ARright[root] <- 0L
    }
  }
  pr$rowroot[row] <- pos
}

# highs_splay_unlink on the subtree whose root is held by holder (a row, or
# a node whose right child holds it)
.hs_splay_unlink_at <- function(pr, node, root, set_root) {
  root <- .hs_splay(pr, pr$Acol[node], root)
  set_root(root)
  if (root != node) {
    .hs_splay_unlink_at(pr, node, pr$ARright[root], function(v) pr$ARright[root] <- v)
    return(invisible())
  }
  if (pr$ARleft[node] == 0L) {
    set_root(pr$ARright[node])
  } else {
    r <- .hs_splay(pr, pr$Acol[node], pr$ARleft[node])
    pr$ARright[r] <- pr$ARright[node]
    set_root(r)
  }
}

.hs_row_preorder <- function(pr, row) {
  out <- integer(); node <- pr$rowroot[row]
  if (node == 0L) return(out)
  stack <- 0L
  repeat {
    out <- c(out, node)
    l <- pr$ARleft[node]; r <- pr$ARright[node]
    if (l != 0L) {
      if (r != 0L) stack <- c(stack, r)
      node <- l
    } else if (r != 0L) {
      node <- r
    } else {
      node <- stack[length(stack)]; stack <- stack[-length(stack)]
      if (node == 0L) break
    }
  }
  out
}
.hs_row_inorder <- function(pr, row) {
  out <- integer()
  walk <- function(node) {
    if (node == 0L) return()
    walk(pr$ARleft[node]); out <<- c(out, node); walk(pr$ARright[node])
  }
  walk(pr$rowroot[row])
  out
}
.hs_col_slots <- function(pr, col) {
  out <- integer(); p <- pr$colhead[col]
  while (p != 0L) { out <- c(out, p); p <- pr$Anext[p] }
  out
}
.hs_find_nonzero <- function(pr, row, col) {
  if (pr$rowroot[row] == 0L) return(0L)
  pr$rowroot[row] <- .hs_splay(pr, col, pr$rowroot[row])
  if (pr$Acol[pr$rowroot[row]] == col) pr$rowroot[row] else 0L
}

# HighsLinearSumBounds: kind "r" sums rows over the column bounds, kind "d"
# sums columns over the row dual bounds
.hs_lsb_new <- function(k) {
  b <- new.env(parent = emptyenv())
  for (w in c("L", "U", "LO", "UO")) {
    b[[paste0(w, "h")]] <- numeric(k); b[[paste0(w, "l")]] <- numeric(k); b[[paste0(w, "n")]] <- integer(k)
  }
  b
}
.hs_lsb_vars <- function(pr, kind) {
  if (kind == "r") list(lo = pr$cl, up = pr$cu, il = pr$implColLower, iu = pr$implColUpper,
                        sl = pr$colLowerSource, su = pr$colUpperSource)
  else list(lo = pr$rdl, up = pr$rdu, il = pr$implRowDualLower, iu = pr$implRowDualUpper,
            sl = pr$rowDualLowerSource, su = pr$rowDualUpperSource)
}
.hs_lsb_impl_up <- function(sum, up, iu, su) if (su == sum) up else min(iu, up)
.hs_lsb_impl_lo <- function(sum, lo, il, sl) if (sl == sum) lo else max(il, lo)
.hs_lsb_update <- function(b, w, sum, dir, bound, coef) {
  if (abs(bound) != Inf) {
    t <- .cd_two_product(bound * dir, coef)
    a <- .cd_iadd(c(b[[paste0(w, "h")]][sum], b[[paste0(w, "l")]][sum]), t)
    b[[paste0(w, "h")]][sum] <- a[1L]; b[[paste0(w, "l")]][sum] <- a[2L]
  } else {
    b[[paste0(w, "n")]][sum] <- b[[paste0(w, "n")]][sum] + dir
  }
}
.hs_lsb_var_up <- function(b, sum, coef, v, dir) .hs_lsb_update(b, if (coef > 0) "UO" else "LO", sum, dir, v, coef)
.hs_lsb_var_lo <- function(b, sum, coef, v, dir) .hs_lsb_update(b, if (coef > 0) "LO" else "UO", sum, dir, v, coef)
.hs_lsb_impl_var_up <- function(b, sum, coef, v, dir) .hs_lsb_update(b, if (coef > 0) "U" else "L", sum, dir, v, coef)
.hs_lsb_impl_var_lo <- function(b, sum, coef, v, dir) .hs_lsb_update(b, if (coef > 0) "L" else "U", sum, dir, v, coef)
.hs_lsb_add <- function(pr, kind, sum, var, coef, dir = 1L) {
  b <- pr$lsb[[kind]]; v <- .hs_lsb_vars(pr, kind)
  .hs_lsb_var_up(b, sum, coef, v$up[var], dir)
  .hs_lsb_var_lo(b, sum, coef, v$lo[var], dir)
  .hs_lsb_impl_var_up(b, sum, coef, .hs_lsb_impl_up(sum, v$up[var], v$iu[var], v$su[var]), dir)
  .hs_lsb_impl_var_lo(b, sum, coef, .hs_lsb_impl_lo(sum, v$lo[var], v$il[var], v$sl[var]), dir)
}
.hs_lsb_updated_var_up <- function(pr, kind, sum, var, coef, old) {
  b <- pr$lsb[[kind]]; v <- .hs_lsb_vars(pr, kind)
  .hs_lsb_var_up(b, sum, coef, old, -1L); .hs_lsb_var_up(b, sum, coef, v$up[var], 1L)
  .hs_lsb_updated_impl_up(pr, kind, sum, var, coef, old, v$iu[var], v$su[var])
}
.hs_lsb_updated_var_lo <- function(pr, kind, sum, var, coef, old) {
  b <- pr$lsb[[kind]]; v <- .hs_lsb_vars(pr, kind)
  .hs_lsb_var_lo(b, sum, coef, old, -1L); .hs_lsb_var_lo(b, sum, coef, v$lo[var], 1L)
  .hs_lsb_updated_impl_lo(pr, kind, sum, var, coef, old, v$il[var], v$sl[var])
}
.hs_lsb_updated_impl_up <- function(pr, kind, sum, var, coef, old_up, old_iu, old_su) {
  b <- pr$lsb[[kind]]; v <- .hs_lsb_vars(pr, kind)
  o <- .hs_lsb_impl_up(sum, old_up, old_iu, old_su)
  nw <- .hs_lsb_impl_up(sum, v$up[var], v$iu[var], v$su[var])
  if (nw == o) return(invisible())
  .hs_lsb_impl_var_up(b, sum, coef, o, -1L); .hs_lsb_impl_var_up(b, sum, coef, nw, 1L)
}
.hs_lsb_updated_impl_lo <- function(pr, kind, sum, var, coef, old_lo, old_il, old_sl) {
  b <- pr$lsb[[kind]]; v <- .hs_lsb_vars(pr, kind)
  o <- .hs_lsb_impl_lo(sum, old_lo, old_il, old_sl)
  nw <- .hs_lsb_impl_lo(sum, v$lo[var], v$il[var], v$sl[var])
  if (nw == o) return(invisible())
  .hs_lsb_impl_var_lo(b, sum, coef, o, -1L); .hs_lsb_impl_var_lo(b, sum, coef, nw, 1L)
}
# getSumLower/Upper(Orig) with a double offset
.hs_lsb_sum <- function(pr, kind, w, sum, offset = 0) {
  b <- pr$lsb[[kind]]
  if (b[[paste0(w, "n")]][sum] != 0L) return(if (w %in% c("L", "LO")) -Inf else Inf)
  .cd_val(.cd_add(c(b[[paste0(w, "h")]][sum], b[[paste0(w, "l")]][sum]), offset))
}
.hs_lsb_ninf <- function(pr, kind, w, sum) pr$lsb[[kind]][[paste0(w, "n")]][sum]
# getResidualSum{Lower,Upper}{,Orig}(sum, var, coef)
.hs_lsb_residual <- function(pr, kind, w, sum, var, coef) {
  b <- pr$lsb[[kind]]; v <- .hs_lsb_vars(pr, kind)
  act <- c(b[[paste0(w, "h")]][sum], b[[paste0(w, "l")]][sum]); ninf <- b[[paste0(w, "n")]][sum]
  lower_side <- w %in% c("L", "LO")
  use_lower <- if (lower_side) coef > 0 else coef < 0
  bound <- if (w %in% c("LO", "UO")) { if (use_lower) v$lo[var] else v$up[var] }
           else if (use_lower) .hs_lsb_impl_lo(sum, v$lo[var], v$il[var], v$sl[var])
           else .hs_lsb_impl_up(sum, v$up[var], v$iu[var], v$su[var])
  if (abs(bound) != Inf) act <- .cd_iadd(act, .cd_two_product(-bound, coef)) else ninf <- ninf - 1L
  if (ninf == 0L) .cd_val(act) else if (lower_side) -Inf else Inf
}

# HPresolve::link and unlink
.hs_link <- function(pr, pos) {
  col <- pr$Acol[pos]; row <- pr$Arow[pos]
  pr$Anext[pos] <- pr$colhead[col]; pr$Aprev[pos] <- 0L
  pr$colhead[col] <- pos
  if (pr$Anext[pos] != 0L) pr$Aprev[pr$Anext[pos]] <- pos
  pr$colsize[col] <- pr$colsize[col] + 1L
  pr$ARleft[pos] <- 0L; pr$ARright[pos] <- 0L
  .hs_splay_link(pr, pos, row)
  .hs_lsb_add(pr, "r", row, col, pr$Avalue[pos])
  .hs_lsb_add(pr, "d", col, row, pr$Avalue[pos])
  pr$rowsize[row] <- pr$rowsize[row] + 1L
}
.hs_unlink <- function(pr, pos) {
  col <- pr$Acol[pos]; row <- pr$Arow[pos]
  nx <- pr$Anext[pos]; pv <- pr$Aprev[pos]
  if (nx != 0L) pr$Aprev[nx] <- pv
  if (pv != 0L) pr$Anext[pv] <- nx else pr$colhead[col] <- nx
  pr$colsize[col] <- pr$colsize[col] - 1L
  if (!pr$colDeleted[col]) {
    if (pr$colsize[col] == 1L) pr$singletonColumns <- c(pr$singletonColumns, col) else .hs_mark_col(pr, col)
    .hs_lsb_add(pr, "d", col, row, pr$Avalue[pos], -1L)
  }
  .hs_splay_unlink_at(pr, pos, pr$rowroot[row], function(v) pr$rowroot[row] <- v)
  pr$rowsize[row] <- pr$rowsize[row] - 1L
  if (!pr$rowDeleted[row]) {
    if (pr$rowsize[row] == 1L) pr$singletonRows <- c(pr$singletonRows, row) else .hs_mark_row(pr, row)
    .hs_lsb_add(pr, "r", row, col, pr$Avalue[pos], -1L)
  }
  .hs_reset_row_dual_from_col(pr, col)
  .hs_reset_col_impl_from_row(pr, row)
  pr$Avalue[pos] <- 0
  pr$freeslots <- c(pr$freeslots, pos)
}
.hs_mark_row <- function(pr, row) {
  if (!pr$changedRowFlag[row]) { pr$changedRowIndices <- c(pr$changedRowIndices, row); pr$changedRowFlag[row] <- TRUE }
}
.hs_mark_col <- function(pr, col) {
  if (!pr$changedColFlag[col]) { pr$changedColIndices <- c(pr$changedColIndices, col); pr$changedColFlag[col] <- TRUE }
}

.hs_add_to_matrix <- function(pr, row, col, val) {
  pos <- .hs_find_nonzero(pr, row, col)
  .hs_mark_row(pr, row); .hs_mark_col(pr, col)
  if (pos == 0L) {
    if (!length(pr$freeslots)) {
      pos <- length(pr$Avalue) + 1L
      pr$Avalue[pos] <- val; pr$Arow[pos] <- row; pr$Acol[pos] <- col
      pr$Anext[pos] <- 0L; pr$Aprev[pos] <- 0L; pr$ARleft[pos] <- 0L; pr$ARright[pos] <- 0L
    } else {
      pos <- pr$freeslots[length(pr$freeslots)]; pr$freeslots <- pr$freeslots[-length(pr$freeslots)]
      pr$Avalue[pos] <- val; pr$Arow[pos] <- row; pr$Acol[pos] <- col; pr$Aprev[pos] <- 0L
    }
    .hs_link(pr, pos)
    .hs_reset_row_dual_from_col(pr, col)
    .hs_reset_col_impl_from_row(pr, row)
  } else {
    s <- pr$Avalue[pos] + val
    if (abs(s) <= 1e-9) {
      .hs_unlink(pr, pos)
    } else {
      .hs_reset_row_dual_from_col(pr, col)
      .hs_reset_col_impl_from_row(pr, row)
      .hs_lsb_add(pr, "r", row, col, pr$Avalue[pos], -1L)
      .hs_lsb_add(pr, "d", col, row, pr$Avalue[pos], -1L)
      pr$Avalue[pos] <- s
      .hs_lsb_add(pr, "r", row, col, s)
      .hs_lsb_add(pr, "d", col, row, s)
    }
  }
}

# implied bound sources
.hs_set_add <- function(s, v) sort(unique(c(s, v)))
.hs_reset_col_impl <- function(pr, col, row) {
  if (!pr$colDeleted[col]) {
    if (pr$colLowerSource[col] != 0L && (row == 0L || pr$colLowerSource[col] == row)) .hs_change_impl_col_lower(pr, col, -Inf, 0L)
    if (pr$colUpperSource[col] != 0L && (row == 0L || pr$colUpperSource[col] == row)) .hs_change_impl_col_upper(pr, col, Inf, 0L)
  } else if (row != 0L && !pr$rowDeleted[row]) {
    pr$colImplSourceByRow[[row]] <- setdiff(pr$colImplSourceByRow[[row]], col)
  }
}
.hs_reset_row_dual_impl <- function(pr, row, col) {
  if (!pr$rowDeleted[row]) {
    if (pr$rowDualLowerSource[row] != 0L && (col == 0L || pr$rowDualLowerSource[row] == col)) .hs_change_impl_row_dual_lower(pr, row, -Inf, 0L)
    if (pr$rowDualUpperSource[row] != 0L && (col == 0L || pr$rowDualUpperSource[row] == col)) .hs_change_impl_row_dual_upper(pr, row, Inf, 0L)
  } else if (col != 0L && !pr$colDeleted[col]) {
    pr$implRowDualSourceByCol[[col]] <- setdiff(pr$implRowDualSourceByCol[[col]], row)
  }
}
.hs_reset_col_impl_from_row <- function(pr, row) {
  for (col in pr$colImplSourceByRow[[row]]) .hs_reset_col_impl(pr, col, row)
}
.hs_reset_row_dual_from_col <- function(pr, col) {
  for (row in pr$implRowDualSourceByCol[[col]]) .hs_reset_row_dual_impl(pr, row, col)
}

.hs_is_lower_implied <- function(pr, col) pr$cl[col] == -Inf || pr$implColLower[col] >= pr$cl[col] - pr$ptol
.hs_is_upper_implied <- function(pr, col) pr$cu[col] == Inf || pr$implColUpper[col] <= pr$cu[col] + pr$ptol
.hs_is_lower_strictly_implied <- function(pr, col, tol) pr$cl[col] == -Inf || pr$implColLower[col] > pr$cl[col] + tol
.hs_is_upper_strictly_implied <- function(pr, col, tol) pr$cu[col] == Inf || pr$implColUpper[col] < pr$cu[col] - tol
.hs_is_implied_free <- function(pr, col) .hs_is_lower_implied(pr, col) && .hs_is_upper_implied(pr, col)
.hs_is_equation <- function(pr, row) pr$rl[row] == pr$ru[row]
.hs_is_ranged <- function(pr, row) pr$rl[row] != -Inf && pr$ru[row] != Inf
.hs_is_dual_implied_free <- function(pr, row) {
  .hs_is_equation(pr, row) || (pr$ru[row] != Inf && pr$implRowDualUpper[row] <= pr$dtol) ||
    (pr$rl[row] != -Inf && pr$implRowDualLower[row] >= -pr$dtol)
}
.hs_is_redundant <- function(pr, row) {
  .hs_lsb_sum(pr, "r", "L", row) >= pr$rl[row] - pr$ptol && .hs_lsb_sum(pr, "r", "U", row) <= pr$ru[row] + pr$ptol
}
.hs_yields_lower <- function(pr, row, val) (val < 0 && pr$ru[row] != Inf) || (val > 0 && pr$rl[row] != -Inf)
.hs_yields_upper <- function(pr, row, val) .hs_yields_lower(pr, row, -val)
.hs_impl_eq_at_lower <- function(pr, row) pr$implRowDualLower[row] > pr$dtol
.hs_impl_eq_at_upper <- function(pr, row) pr$implRowDualUpper[row] < -pr$dtol

.hs_change_col_upper <- function(pr, col, v) {
  old <- pr$cu[col]; pr$cu[col] <- v
  for (p in .hs_col_slots(pr, col)) {
    .hs_lsb_updated_var_up(pr, "r", pr$Arow[p], col, pr$Avalue[p], old)
    .hs_mark_row(pr, pr$Arow[p])
  }
}
.hs_change_col_lower <- function(pr, col, v) {
  old <- pr$cl[col]; pr$cl[col] <- v
  for (p in .hs_col_slots(pr, col)) {
    .hs_lsb_updated_var_lo(pr, "r", pr$Arow[p], col, pr$Avalue[p], old)
    .hs_mark_row(pr, pr$Arow[p])
  }
}
.hs_change_row_dual_upper <- function(pr, row, v) {
  old <- pr$rdu[row]; pr$rdu[row] <- v
  for (p in .hs_row_preorder(pr, row)) {
    .hs_lsb_updated_var_up(pr, "d", pr$Acol[p], row, pr$Avalue[p], old)
    .hs_mark_col(pr, pr$Acol[p])
  }
}
.hs_change_row_dual_lower <- function(pr, row, v) {
  old <- pr$rdl[row]; pr$rdl[row] <- v
  for (p in .hs_row_preorder(pr, row)) {
    .hs_lsb_updated_var_lo(pr, "d", pr$Acol[p], row, pr$Avalue[p], old)
    .hs_mark_col(pr, pr$Acol[p])
  }
}
.hs_change_impl_col_upper <- function(pr, col, v, origin) {
  old <- pr$implColUpper[col]; osrc <- pr$colUpperSource[col]; tol <- pr$ptol
  if (old >= pr$cu[col] - tol && v < pr$cu[col] - tol) .hs_mark_col(pr, col)
  nif <- .hs_is_lower_implied(pr, col) && old > pr$cu[col] + tol && v <= pr$cu[col] + tol
  if (osrc != origin) {
    if (osrc != 0L && osrc != pr$colLowerSource[col]) pr$colImplSourceByRow[[osrc]] <- setdiff(pr$colImplSourceByRow[[osrc]], col)
    if (origin != 0L) pr$colImplSourceByRow[[origin]] <- .hs_set_add(pr$colImplSourceByRow[[origin]], col)
    pr$colUpperSource[col] <- origin
  }
  pr$implColUpper[col] <- v
  if (!nif && min(old, v) >= pr$cu[col]) return(invisible())
  for (p in .hs_col_slots(pr, col)) {
    .hs_lsb_updated_impl_up(pr, "r", pr$Arow[p], col, pr$Avalue[p], pr$cu[col], old, osrc)
    if (nif && .hs_is_dual_implied_free(pr, pr$Arow[p])) .hs_push_subst(pr, pr$Arow[p], col)
    .hs_mark_row(pr, pr$Arow[p])
  }
}
.hs_change_impl_col_lower <- function(pr, col, v, origin) {
  old <- pr$implColLower[col]; osrc <- pr$colLowerSource[col]; tol <- pr$ptol
  if (old <= pr$cl[col] + tol && v > pr$cl[col] + tol) .hs_mark_col(pr, col)
  nif <- .hs_is_upper_implied(pr, col) && old < pr$cl[col] - tol && v >= pr$cl[col] - tol
  if (osrc != origin) {
    if (osrc != 0L && osrc != pr$colUpperSource[col]) pr$colImplSourceByRow[[osrc]] <- setdiff(pr$colImplSourceByRow[[osrc]], col)
    if (origin != 0L) pr$colImplSourceByRow[[origin]] <- .hs_set_add(pr$colImplSourceByRow[[origin]], col)
    pr$colLowerSource[col] <- origin
  }
  pr$implColLower[col] <- v
  if (!nif && max(old, v) <= pr$cl[col]) return(invisible())
  for (p in .hs_col_slots(pr, col)) {
    .hs_lsb_updated_impl_lo(pr, "r", pr$Arow[p], col, pr$Avalue[p], pr$cl[col], old, osrc)
    if (nif && .hs_is_dual_implied_free(pr, pr$Arow[p])) .hs_push_subst(pr, pr$Arow[p], col)
    .hs_mark_row(pr, pr$Arow[p])
  }
}
.hs_change_impl_row_dual_upper <- function(pr, row, v, origin) {
  old <- pr$implRowDualUpper[row]; osrc <- pr$rowDualUpperSource[row]; tol <- pr$dtol
  if (old >= -tol && v < -tol) .hs_mark_row(pr, row)
  ndi <- !.hs_is_dual_implied_free(pr, row) && old > pr$rdu[row] + tol && v <= pr$rdu[row] + tol
  if (osrc != origin) {
    if (osrc != 0L && osrc != pr$rowDualLowerSource[row]) pr$implRowDualSourceByCol[[osrc]] <- setdiff(pr$implRowDualSourceByCol[[osrc]], row)
    if (origin != 0L) pr$implRowDualSourceByCol[[origin]] <- .hs_set_add(pr$implRowDualSourceByCol[[origin]], row)
    pr$rowDualUpperSource[row] <- origin
  }
  pr$implRowDualUpper[row] <- v
  if (!ndi && min(old, v) >= pr$rdu[row]) return(invisible())
  for (p in .hs_row_preorder(pr, row)) {
    .hs_lsb_updated_impl_up(pr, "d", pr$Acol[p], row, pr$Avalue[p], pr$rdu[row], old, osrc)
    .hs_mark_col(pr, pr$Acol[p])
    if (ndi && .hs_is_implied_free(pr, pr$Acol[p])) .hs_push_subst(pr, row, pr$Acol[p])
  }
}
.hs_change_impl_row_dual_lower <- function(pr, row, v, origin) {
  old <- pr$implRowDualLower[row]; osrc <- pr$rowDualLowerSource[row]; tol <- pr$dtol
  if (old <= tol && v > tol) .hs_mark_row(pr, row)
  ndi <- !.hs_is_dual_implied_free(pr, row) && old < pr$rdl[row] - tol && v >= pr$rdl[row] - tol
  if (osrc != origin) {
    if (osrc != 0L && osrc != pr$rowDualUpperSource[row]) pr$implRowDualSourceByCol[[osrc]] <- setdiff(pr$implRowDualSourceByCol[[osrc]], row)
    if (origin != 0L) pr$implRowDualSourceByCol[[origin]] <- .hs_set_add(pr$implRowDualSourceByCol[[origin]], row)
    pr$rowDualLowerSource[row] <- origin
  }
  pr$implRowDualLower[row] <- v
  if (!ndi && max(old, v) <= pr$rdl[row]) return(invisible())
  for (p in .hs_row_preorder(pr, row)) {
    .hs_lsb_updated_impl_lo(pr, "d", pr$Acol[p], row, pr$Avalue[p], pr$rdl[row], old, osrc)
    .hs_mark_col(pr, pr$Acol[p])
    if (ndi && .hs_is_implied_free(pr, pr$Acol[p])) .hs_push_subst(pr, row, pr$Acol[p])
  }
}

# updateColImpliedBounds and updateRowDualImpliedBounds
.hs_check_update_col_impl <- function(pr, row) {
  lo <- if (.hs_impl_eq_at_upper(pr, row)) pr$ru[row] else pr$rl[row]
  up <- if (.hs_impl_eq_at_lower(pr, row)) pr$rl[row] else pr$ru[row]
  ok <- (lo != -Inf && .hs_lsb_ninf(pr, "r", "UO", row) <= 1L) || (up != Inf && .hs_lsb_ninf(pr, "r", "LO", row) <= 1L)
  list(ok = ok, lo = lo, up = up)
}
.hs_update_col_impl_one <- function(pr, row, col, val) {
  ck <- .hs_check_update_col_impl(pr, row)
  if (!ck$ok) return(invisible())
  thr <- 1000 * pr$ptol
  check <- function(bnd, residual, dir) {
    if (dir * residual <= -Inf) return()
    ib <- .cd_val(.cd_div(.cd_add(c(bnd, 0), -residual), val))
    if (abs(ib) * 1e-14 > pr$ptol) return()
    if (dir * val > 0) {
      if (ib < pr$implColUpper[col] - thr) .hs_change_impl_col_upper(pr, col, ib, row)
    } else if (ib > pr$implColLower[col] + thr) .hs_change_impl_col_lower(pr, col, ib, row)
  }
  if (ck$up != Inf) check(ck$up, .hs_lsb_residual(pr, "r", "LO", row, col, val), 1L)
  if (ck$lo != -Inf) check(ck$lo, .hs_lsb_residual(pr, "r", "UO", row, col, val), -1L)
}
.hs_update_col_impl <- function(pr, row) {
  if (!.hs_check_update_col_impl(pr, row)$ok) return(invisible())
  for (p in .hs_row_preorder(pr, row)) .hs_update_col_impl_one(pr, row, pr$Acol[p], pr$Avalue[p])
}
.hs_check_update_row_dual <- function(pr, col) {
  margin <- if (pr$colsize[col] != 1L) pr$ptol else -pr$ptol
  lo <- if (.hs_is_lower_strictly_implied(pr, col, margin)) pr$cost[col] else -Inf
  up <- if (.hs_is_upper_strictly_implied(pr, col, margin)) pr$cost[col] else Inf
  ok <- (lo != -Inf && .hs_lsb_ninf(pr, "d", "UO", col) <= 1L) || (up != Inf && .hs_lsb_ninf(pr, "d", "LO", col) <= 1L)
  list(ok = ok, lo = lo, up = up)
}
.hs_update_row_dual_one <- function(pr, row, col, val) {
  ck <- .hs_check_update_row_dual(pr, col)
  if (!ck$ok) return(invisible())
  thr <- 1000 * pr$dtol
  check <- function(bnd, residual, dir) {
    if (dir * residual <= -Inf) return()
    ib <- .cd_val(.cd_div(.cd_add(c(bnd, 0), -residual), val))
    if (abs(ib) * 1e-14 > pr$dtol) return()
    if (dir * val > 0) {
      if (ib < pr$implRowDualUpper[row] - thr) .hs_change_impl_row_dual_upper(pr, row, ib, col)
    } else if (ib > pr$implRowDualLower[row] + thr) .hs_change_impl_row_dual_lower(pr, row, ib, col)
  }
  if (ck$up != Inf) check(ck$up, .hs_lsb_residual(pr, "d", "LO", col, row, val), 1L)
  if (ck$lo != -Inf) check(ck$lo, .hs_lsb_residual(pr, "d", "UO", col, row, val), -1L)
}
.hs_update_row_dual <- function(pr, col) {
  if (!.hs_check_update_row_dual(pr, col)$ok) return(invisible())
  for (p in .hs_col_slots(pr, col)) .hs_update_row_dual_one(pr, pr$Arow[p], col, pr$Avalue[p])
}

.hs_pre_abort <- function(result) stop(structure(class = c("hs_presolve_abort", "condition"),
                                                 list(message = result, call = NULL, result = result)))

# postsolve stack records (vectors of (index, value) are copied)
.hs_ps_push <- function(pr, rec) pr$stack[[length(pr$stack) + 1L]] <- rec
.hs_col_vec <- function(pr, col) {
  s <- .hs_col_slots(pr, col)
  list(index = pr$Arow[s], value = pr$Avalue[s])
}
.hs_row_vec <- function(pr, row, sorted = FALSE) {
  s <- if (sorted) .hs_row_inorder(pr, row) else .hs_row_preorder(pr, row)
  list(index = pr$Acol[s], value = pr$Avalue[s])
}
.hs_ps_fixed_col <- function(pr, col, fix, type) {
  .hs_ps_push(pr, list(type = "fixedcol", col = col, fix = fix, cost = pr$cost[col], status = type,
                       colv = .hs_col_vec(pr, col)))
}
.hs_ps_redundant_row <- function(pr, row) .hs_ps_push(pr, list(type = "redundantrow", row = row))

.hs_mark_row_deleted <- function(pr, row) {
  if (.hs_is_equation(pr, row)) pr$eqsize[row] <- NA_integer_
  pr$changedRowFlag[row] <- TRUE
  pr$rowDeleted[row] <- TRUE
  pr$numDeletedRows <- pr$numDeletedRows + 1L
}
.hs_mark_col_deleted <- function(pr, col) {
  pr$changedColFlag[col] <- TRUE
  pr$colDeleted[col] <- TRUE
  pr$numDeletedCols <- pr$numDeletedCols + 1L
}
.hs_reinsert_equation <- function(pr, row) {
  if (.hs_is_equation(pr, row) && !is.na(pr$eqsize[row]) && pr$eqsize[row] != pr$rowsize[row])
    pr$eqsize[row] <- pr$rowsize[row]
}
.hs_remove_row <- function(pr, row) {
  .hs_mark_row_deleted(pr, row)
  for (p in .hs_row_inorder(pr, row)) .hs_unlink(pr, p)
}
.hs_shift_rows_for_col <- function(pr, col, fixval) {
  p <- pr$colhead[col]
  while (p != 0L) {
    row <- pr$Arow[p]; v <- pr$Avalue[p]; nx <- pr$Anext[p]
    if (!is.na(fixval)) {
      if (pr$rl[row] != -Inf) pr$rl[row] <- pr$rl[row] - v * fixval
      if (pr$ru[row] != Inf) pr$ru[row] <- pr$ru[row] - v * fixval
    }
    .hs_unlink(pr, p)
    .hs_reinsert_equation(pr, row)
    p <- nx
  }
}
.hs_remove_fixed_col <- function(pr, col) {
  fixval <- pr$cl[col]
  .hs_mark_col_deleted(pr, col)
  .hs_shift_rows_for_col(pr, col, fixval)
  pr$offset <- pr$offset + pr$cost[col] * fixval
  pr$cost[col] <- 0
}
.hs_fix_col_to <- function(pr, col, upper) {
  fixval <- if (upper) pr$cu[col] else pr$cl[col]
  if (abs(fixval) == Inf) .hs_pre_abort("dualinfeasible")
  .hs_ps_fixed_col(pr, col, fixval, if (upper) "upper" else "lower")
  .hs_mark_col_deleted(pr, col)
  .hs_shift_rows_for_col(pr, col, fixval)
  pr$offset <- pr$offset + pr$cost[col] * fixval
  pr$cost[col] <- 0
}
.hs_fix_col_to_zero <- function(pr, col) {
  .hs_ps_fixed_col(pr, col, 0, "zero")
  .hs_mark_col_deleted(pr, col)
  .hs_shift_rows_for_col(pr, col, NA_real_)
  pr$cost[col] <- 0
}

# substitute(substcol, staycol, offset, scale)
.hs_substitute_col <- function(pr, substcol, staycol, offset, scale) {
  p <- pr$colhead[substcol]
  while (p != 0L) {
    row <- pr$Arow[p]; v <- pr$Avalue[p]; nx <- pr$Anext[p]
    .hs_unlink(pr, p)
    if (pr$rl[row] != -Inf) pr$rl[row] <- pr$rl[row] - v * offset
    if (pr$ru[row] != Inf) pr$ru[row] <- pr$ru[row] - v * offset
    .hs_add_to_matrix(pr, row, staycol, scale * v)
    .hs_reinsert_equation(pr, row)
    p <- nx
  }
  if (pr$cost[substcol] != 0) {
    pr$offset <- pr$offset + pr$cost[substcol] * offset
    pr$cost[staycol] <- pr$cost[staycol] + scale * pr$cost[substcol]
    if (abs(pr$cost[staycol]) <= 1e-9) pr$cost[staycol] <- 0
    pr$cost[substcol] <- 0
  }
}
# substitute(row, col, rhs) with the stored row
.hs_substitute_free <- function(pr, row, col, rhs, rowpos) {
  pos <- .hs_find_nonzero(pr, row, col)
  srs <- -1 / pr$Avalue[pos]
  .hs_mark_row_deleted(pr, row); .hs_mark_col_deleted(pr, col)
  p <- pr$colhead[col]
  while (p != 0L) {
    crow <- pr$Arow[p]; cval <- pr$Avalue[p]; nx <- pr$Anext[p]
    if (crow != row) {
      .hs_unlink(pr, p)
      sc <- cval * srs
      if (pr$rl[crow] != -Inf) pr$rl[crow] <- pr$rl[crow] + sc * rhs
      if (pr$ru[crow] != Inf) pr$ru[crow] <- pr$ru[crow] + sc * rhs
      for (q in rowpos) if (pr$Acol[q] != col) .hs_add_to_matrix(pr, crow, pr$Acol[q], sc * pr$Avalue[q])
      .hs_reinsert_equation(pr, crow)
    }
    p <- nx
  }
  if (pr$cost[col] != 0) {
    objscale <- c(pr$cost[col] * srs, 0)
    t <- .cd_mul(objscale, rhs)
    pr$offset <- .cd_val(.cd_add(-t, pr$offset))
    for (q in rowpos) {
      j <- pr$Acol[q]
      pr$cost[j] <- .cd_val(.cd_add(.cd_mul(objscale, pr$Avalue[q]), pr$cost[j]))
      if (abs(pr$cost[j]) <= 1e-9) pr$cost[j] <- 0
    }
    pr$cost[col] <- 0
  }
  for (q in rowpos) .hs_unlink(pr, q)
}

.hs_remove_row_singletons <- function(pr) {
  i <- 1L
  while (i <= length(pr$singletonRows)) {
    row <- pr$singletonRows[i]
    if (!pr$rowDeleted[row] && pr$rowsize[row] <= 1L) .hs_row_presolve(pr, row)
    i <- i + 1L
  }
  pr$singletonRows <- integer()
}

.hs_singleton_row <- function(pr, row) {
  pos <- pr$rowroot[row]
  col <- pr$Acol[pos]; val <- pr$Avalue[pos]; tol <- pr$ptol
  .hs_mark_row_deleted(pr, row)
  .hs_unlink(pr, pos)
  lo <- pr$cl[col]; up <- pr$cu[col]; rl <- pr$rl[row]; ru <- pr$ru[row]
  redundant <- if (val > 0) up * val <= ru + tol && lo * val >= rl - tol
               else lo * val <= ru + tol && up * val >= rl - tol
  if (redundant) {
    .hs_ps_redundant_row(pr, row)
    return(invisible())
  }
  ncu <- Inf; ncl <- -Inf
  if (val > 0) {
    if (ru != Inf) ncu <- ru / val
    if (rl != -Inf) ncl <- rl / val
  } else {
    if (ru != Inf) ncl <- ru / val
    if (rl != -Inf) ncu <- rl / val
  }
  btol <- tol / max(1, abs(val))
  lt <- ncl > lo + btol; ut <- ncu < up - btol
  lb <- if (lt) ncl else lo
  ub <- if (ut) ncu else up
  if (ub <= lb + tol) {
    if (ub < lb - tol) .hs_pre_abort("infeasible")
    maxcol <- if (pr$colsize[col]) max(abs(pr$Avalue[.hs_col_slots(pr, col)])) else 0
    if (ub < lb || (ub > lb && (ub - lb) * max(abs(val), maxcol) <= tol)) {
      if (lt && ut) {
        ub <- 0.5 * (ub + lb); lb <- ub
        lt <- lb > lo; ut <- ub < up
      } else if (lt) {
        lb <- ub; lt <- lb > lo
      } else {
        ub <- lb; ut <- ub < up
      }
    }
  }
  .hs_ps_push(pr, list(type = "singletonrow", row = row, col = col, coef = val, lt = lt, ut = ut))
  if (lt) .hs_change_col_lower(pr, col, lb)
  if (ub == lb) {
    .hs_ps_fixed_col(pr, col, lb, "nonbasic")
    .hs_remove_fixed_col(pr, col)
  } else if (ut) .hs_change_col_upper(pr, col, ub)
  if (!pr$colDeleted[col] && pr$colsize[col] == 0L) .hs_empty_col(pr, col)
}

.hs_doubleton_eq <- function(pr, row, row_type) {
  p1 <- pr$rowroot[row]
  p2 <- if (pr$ARright[p1] != 0L) pr$ARright[p1] else pr$ARleft[p1]
  c1 <- pr$Acol[p1]; c2 <- pr$Acol[p2]
  better1 <- {
    s1 <- pr$colsize[c1]
    if (s1 == 1L) TRUE else {
      s2 <- pr$colsize[c2]
      if (s2 == 1L) FALSE else {
        a1 <- abs(pr$Avalue[p1]); a2 <- abs(pr$Avalue[p2])
        if (s1 != s2 && max(a1, a2) <= 2 * min(a1, a2)) s1 < s2 else a1 > a2
      }
    }
  }
  if (better1) { sc <- c1; st <- c2; scoef <- pr$Avalue[p1]; tcoef <- pr$Avalue[p2] }
  else { sc <- c2; st <- c1; scoef <- pr$Avalue[p2]; tcoef <- pr$Avalue[p1] }
  rhs <- pr$ru[row]
  olo <- pr$cl[st]; oup <- pr$cu[st]; slo <- pr$cl[sc]; sup <- pr$cu[sc]
  imp <- function(b) .cd_val(.cd_div(.cd_add(c(rhs, 0), -(scoef * b)), tcoef))
  if (sign(scoef) != sign(tcoef)) {
    ilo <- if (slo == -Inf) -Inf else imp(slo)
    iup <- if (sup == Inf) Inf else imp(sup)
  } else {
    ilo <- if (sup == Inf) -Inf else imp(sup)
    iup <- if (slo == -Inf) Inf else imp(slo)
  }
  lt <- ilo > olo + pr$ptol
  if (lt) .hs_change_col_lower(pr, st, ilo)
  ut <- iup < oup - pr$ptol
  if (ut) .hs_change_col_upper(pr, st, iup)
  .hs_ps_push(pr, list(type = "doubleton", row = row, colsubst = sc, col = st, coefsubst = scoef, coef = tcoef,
                       rhs = rhs, substcost = pr$cost[sc], lt = lt, ut = ut, rowtype = row_type,
                       colv = .hs_col_vec(pr, sc)))
  .hs_mark_col_deleted(pr, sc)
  .hs_remove_row(pr, row)
  .hs_substitute_col(pr, sc, st, rhs / scoef, -tcoef / scoef)
  .hs_remove_row_singletons(pr)
}

.hs_row_presolve <- function(pr, row) {
  tol <- pr$ptol
  if (.hs_lsb_sum(pr, "r", "L", row) > pr$ru[row] + tol || .hs_lsb_sum(pr, "r", "U", row) < pr$rl[row] - tol)
    .hs_pre_abort("infeasible")
  if (pr$rowsize[row] == 1L) return(.hs_singleton_row(pr, row))
  if (.hs_is_redundant(pr, row)) {
    .hs_ps_redundant_row(pr, row)
    .hs_remove_row(pr, row)
    return(invisible())
  }
  redundant_bounds <- function(col) {
    if (pr$colsize[col] != 1L) return()
    if (pr$cost[col] > 0) {
      if (pr$cl[col] > pr$implColLower[col] - tol) .hs_change_col_lower(pr, col, -Inf)
    } else if (pr$cu[col] < pr$implColUpper[col] + tol) .hs_change_col_upper(pr, col, Inf)
  }
  orl <- pr$rl[row]; oru <- pr$ru[row]
  if (!.hs_is_equation(pr, row)) {
    if (.hs_impl_eq_at_lower(pr, row)) {
      pr$ru[row] <- pr$rl[row]
      .hs_change_row_dual_lower(pr, row, -Inf)
      redundant_bounds(pr$rowDualLowerSource[row])
    } else if (.hs_impl_eq_at_upper(pr, row)) {
      pr$rl[row] <- pr$ru[row]
      .hs_change_row_dual_upper(pr, row, Inf)
      redundant_bounds(pr$rowDualUpperSource[row])
    }
  }
  if (pr$rowsize[row] == 2L && pr$rl[row] == pr$ru[row]) {
    rt <- if (orl == oru) "eq" else if (oru != Inf) "leq" else "geq"
    return(.hs_doubleton_eq(pr, row, rt))
  }
  forcing <- function(dir, side, implied, min_abs, row_type) {
    if (dir * implied > dir * side + tol * min(1, min_abs)) return()
    rv <- .hs_row_vec(pr, row, TRUE)
    nfix <- 0L
    for (k in seq_along(rv$index)) {
      j <- rv$index[k]
      if (dir * rv$value[k] > 0) { if (pr$cu[j] <= pr$implColUpper[j]) nfix <- nfix + 1L }
      else if (pr$cl[j] >= pr$implColLower[j]) nfix <- nfix + 1L
    }
    if (nfix != pr$rowsize[row]) return()
    .hs_ps_push(pr, list(type = "forcingrow", row = row, side = side, rowtype = row_type, rowv = rv))
    .hs_mark_row_deleted(pr, row)
    for (k in seq_along(rv$index)) {
      j <- rv$index[k]
      if (dir * rv$value[k] > 0) {
        .hs_ps_fixed_col(pr, j, pr$cu[j], "upper")
        if (pr$cl[j] < pr$cu[j]) .hs_change_col_lower(pr, j, pr$cu[j])
      } else {
        .hs_ps_fixed_col(pr, j, pr$cl[j], "lower")
        if (pr$cu[j] > pr$cl[j]) .hs_change_col_upper(pr, j, pr$cl[j])
      }
      .hs_remove_fixed_col(pr, j)
    }
    .hs_ps_redundant_row(pr, row)
    .hs_remove_row_singletons(pr)
  }
  min_abs <- min(abs(.hs_row_vec(pr, row, TRUE)$value))
  forcing(1L, pr$rl[row], .hs_lsb_sum(pr, "r", "UO", row), min_abs, "geq")
  if (pr$rowDeleted[row]) return(invisible())
  forcing(-1L, pr$ru[row], .hs_lsb_sum(pr, "r", "LO", row), min_abs, "leq")
  if (pr$rowDeleted[row]) return(invisible())
  .hs_update_col_impl(pr, row)
}

.hs_empty_col <- function(pr, col) {
  if ((pr$cost[col] > 0 && pr$cl[col] == -Inf) || (pr$cost[col] < 0 && pr$cu[col] == Inf)) {
    if (abs(pr$cost[col]) <= pr$dtol) pr$cost[col] <- 0 else .hs_pre_abort("dualinfeasible")
  }
  if (pr$cost[col] > 0) .hs_fix_col_to(pr, col, FALSE)
  else if (pr$cost[col] < 0 || abs(pr$cu[col]) < abs(pr$cl[col])) .hs_fix_col_to(pr, col, TRUE)
  else if (pr$cl[col] != -Inf) .hs_fix_col_to(pr, col, FALSE)
  else .hs_fix_col_to_zero(pr, col)
}

# computeColBounds without a bound column: the worst case bounds
.hs_worst_case_bound <- function(pr, col, lower) {
  out <- if (lower) -Inf else Inf
  for (p in .hs_col_slots(pr, col)) {
    row <- pr$Arow[p]; jv <- pr$Avalue[p]
    hl <- .hs_yields_lower(pr, row, jv); hu <- .hs_yields_upper(pr, row, jv)
    if ((lower && !hl) || (!lower && !hu)) next
    for (dir in c(1L, -1L)) {
      rhs <- if (dir == 1L) pr$ru[row] else pr$rl[row]
      if (dir * rhs == Inf) next
      wants_lower <- dir * jv < 0
      if (wants_lower != lower) next
      if (lower && out == Inf) next
      if (!lower && out == -Inf) next
      res <- .hs_lsb_residual(pr, "r", if (dir < 0) "L" else "U", row, col, jv)
      b <- if (dir < 0) {
        if (res == -Inf) sign(jv) * Inf else NA
      } else if (res == Inf) -sign(jv) * Inf else NA
      if (is.na(b)) b <- .cd_val(.cd_div(.cd_add(c(rhs, 0), -res), jv))
      out <- if (lower) max(out, b) else min(out, b)
    }
  }
  out
}

.hs_detect_dominated_col <- function(pr, col, handle_singletons = TRUE) {
  dt <- pr$dtol
  dual_up <- -.hs_lsb_sum(pr, "d", "L", col, -pr$cost[col])
  dual_lo <- -.hs_lsb_sum(pr, "d", "U", col, -pr$cost[col])
  dominated <- function(dbound, bound, dir) {
    if (dir * dbound <= dt) return()
    if (dir * bound == -Inf) .hs_pre_abort("dualinfeasible")
    .hs_fix_col_to(pr, col, dir < 0)
    if (handle_singletons) .hs_remove_row_singletons(pr)
  }
  weakly <- function(dbound, bound, other, dir) {
    if (dir * dbound < -dt) return()
    if (dir * bound != -Inf) {
      .hs_fix_col_to(pr, col, dir < 0)
      if (handle_singletons) .hs_remove_row_singletons(pr)
      return()
    }
    bod <- if (dir > 0) -.hs_lsb_sum(pr, "d", "UO", col, -pr$cost[col]) else -.hs_lsb_sum(pr, "d", "LO", col, -pr$cost[col])
    if (bod != 0) return()
    .hs_ps_push(pr, list(type = "forcingcol", col = col, cost = pr$cost[col], bound = other, at_inf_upper = dir < 0,
                         colv = .hs_col_vec(pr, col)))
    .hs_mark_col_deleted(pr, col)
    p <- pr$colhead[col]
    while (p != 0L) {
      row <- pr$Arow[p]
      rhs <- if (dir * pr$Avalue[p] > 0) pr$ru[row] else pr$rl[row]
      p <- pr$Anext[p]
      rv <- .hs_row_vec(pr, row)
      keep <- rv$index != col
      .hs_ps_push(pr, list(type = "forcingcolrow", row = row, rhs = rhs,
                           rowv = list(index = rv$index[keep], value = rv$value[keep])))
      .hs_remove_row(pr, row)
    }
  }
  dominated(dual_lo, pr$cl[col], 1L); if (pr$colDeleted[col]) return(invisible())
  dominated(dual_up, pr$cu[col], -1L); if (pr$colDeleted[col]) return(invisible())
  weakly(dual_lo, pr$cl[col], pr$cu[col], 1L); if (pr$colDeleted[col]) return(invisible())
  weakly(dual_up, pr$cu[col], pr$cl[col], -1L)
}

.hs_dual_fixing <- function(pr, col) {
  if (pr$cl[col] == pr$cu[col]) return(invisible())
  down <- 0L; up <- 0L
  if (pr$cost[col] > 0) up <- up + 1L else if (pr$cost[col] < 0) down <- down + 1L
  for (p in .hs_col_slots(pr, col)) {
    if (.hs_yields_lower(pr, pr$Arow[p], pr$Avalue[p])) down <- down + 1L
    if (.hs_yields_upper(pr, pr$Arow[p], pr$Avalue[p])) up <- up + 1L
    if (down > 1L && up > 1L) break
  }
  if (down == 0L || up == 0L) return(.hs_fix_col_to(pr, col, down != 0L))
  tighter <- function(dir, current) {
    if (dir * pr$cost[col] < 0) return(NULL)
    nb <- if (dir > 0) .hs_worst_case_bound(pr, col, TRUE) else -.hs_worst_case_bound(pr, col, FALSE)
    if (nb != -Inf && (nb >= dir * current - pr$ptol || abs(nb) > pr$ptol / 1e-14)) return(NULL)
    nb * dir
  }
  nb <- tighter(1L, pr$cu[col])
  if (!is.null(nb)) {
    nb <- max(nb, pr$cl[col])
    if (nb < pr$cu[col] - pr$ptol && nb == pr$cl[col]) .hs_fix_col_to(pr, col, FALSE)
  } else {
    nb <- tighter(-1L, pr$cl[col])
    if (!is.null(nb)) {
      nb <- min(nb, pr$cu[col])
      if (nb > pr$cl[col] + pr$ptol && nb == pr$cu[col]) .hs_fix_col_to(pr, col, TRUE)
    }
  }
}

.hs_singleton_col <- function(pr, col) {
  pos <- pr$colhead[col]; row <- pr$Arow[pos]; coef <- pr$Avalue[pos]
  if (pr$rowsize[row] == 1L) {
    .hs_singleton_row(pr, row)
    if (!pr$colDeleted[col]) .hs_empty_col(pr, col)
    return(invisible())
  }
  .hs_detect_dominated_col(pr, col, FALSE)
  if (pr$colDeleted[col]) return(invisible())
  .hs_dual_fixing(pr, col)
  if (pr$colDeleted[col]) return(invisible())
  .hs_update_col_impl_one(pr, row, col, coef)
  .hs_update_row_dual_one(pr, row, col, coef)
  if (.hs_is_dual_implied_free(pr, row) && .hs_is_implied_free(pr, col)) {
    rowpos <- .hs_row_inorder(pr, row)
    .hs_substitute_free_col(pr, row, col, rowpos)
  }
}

.hs_substitute_free_col <- function(pr, row, col, rowpos, relax = FALSE) {
  if (.hs_is_equation(pr, row)) { rt <- "eq"; rhs <- pr$ru[row] }
  else if (pr$ru[row] != Inf && pr$implRowDualUpper[row] <= pr$dtol) {
    rt <- "leq"; rhs <- pr$ru[row]
    if (relax) .hs_change_row_dual_upper(pr, row, Inf)
  } else {
    rt <- "geq"; rhs <- pr$rl[row]
    if (relax) .hs_change_row_dual_lower(pr, row, -Inf)
  }
  .hs_ps_push(pr, list(type = "freecol", row = row, col = col, rhs = rhs, cost = pr$cost[col], rowtype = rt,
                       rowv = list(index = pr$Acol[rowpos], value = pr$Avalue[rowpos]), colv = .hs_col_vec(pr, col)))
  .hs_substitute_free(pr, row, col, rhs, rowpos)
}

.hs_col_presolve <- function(pr, col) {
  diff <- pr$cu[col] - pr$cl[col]
  if (diff <= pr$ptol) {
    maxv <- if (pr$colsize[col]) max(abs(pr$Avalue[.hs_col_slots(pr, col)])) else 0
    if (diff <= 1e-9 || maxv * diff <= pr$ptol) {
      if (diff < -pr$ptol) .hs_pre_abort("infeasible")
      if (abs(pr$cl[col]) == Inf) .hs_pre_abort("dualinfeasible")
      .hs_ps_fixed_col(pr, col, pr$cl[col], "nonbasic")
      .hs_remove_fixed_col(pr, col)
      return(invisible())
    }
  }
  if (pr$colsize[col] == 0L) return(.hs_empty_col(pr, col))
  if (pr$colsize[col] == 1L) return(.hs_singleton_col(pr, col))
  .hs_detect_dominated_col(pr, col)
  if (pr$colDeleted[col]) return(invisible())
  .hs_dual_fixing(pr, col)
  if (pr$colDeleted[col]) return(invisible())
  .hs_update_row_dual(pr, col)
}

.hs_presolve_changed_rows <- function(pr) {
  rows <- pr$changedRowIndices; pr$changedRowIndices <- integer()
  for (row in rows) {
    if (pr$rowDeleted[row]) next
    .hs_row_presolve(pr, row)
    pr$changedRowFlag[row] <- pr$rowDeleted[row]
  }
}
.hs_presolve_changed_cols <- function(pr) {
  cols <- pr$changedColIndices; pr$changedColIndices <- integer()
  for (col in cols) {
    if (pr$colDeleted[col]) next
    .hs_col_presolve(pr, col)
    pr$changedColFlag[col] <- pr$colDeleted[col]
  }
}
.hs_presolve_col_singletons <- function(pr) {
  i <- 1L
  while (i <= length(pr$singletonColumns)) {
    col <- pr$singletonColumns[i]
    if (!pr$colDeleted[col]) .hs_col_presolve(pr, col)
    i <- i + 1L
  }
  s <- pr$singletonColumns
  pr$singletonColumns <- s[!pr$colDeleted[s] & pr$colsize[s] <= 1L]
}
# the next equation of the ordered set after (size, row), or the first
.hs_next_equation <- function(pr, after = NULL) {
  rows <- which(!is.na(pr$eqsize))
  if (!length(rows)) return(0L)
  key <- pr$eqsize[rows] * (pr$m + 1) + rows
  if (!is.null(after)) rows <- rows[key > after]
  if (!length(rows)) return(0L)
  key <- pr$eqsize[rows] * (pr$m + 1) + rows
  rows[which.min(key)]
}
.hs_remove_doubleton_equations <- function(pr) {
  row <- .hs_next_equation(pr)
  while (row != 0L) {
    if (pr$rowsize[row] > 2L) return(invisible())
    key <- pr$eqsize[row] * (pr$m + 1) + row
    .hs_row_presolve(pr, row)
    row <- if (pr$rowDeleted[row]) .hs_next_equation(pr) else .hs_next_equation(pr, key)
  }
}

.hs_live_cols <- function(pr) pr$ncol_now - pr$numDeletedCols
.hs_live_rows <- function(pr) pr$nrow_now - pr$numDeletedRows
.hs_store_size <- function(pr) { pr$oldNumCol <- .hs_live_cols(pr); pr$oldNumRow <- .hs_live_rows(pr) }
.hs_size_reduction <- function(pr) {
  r <- 100 * (pr$oldNumRow - .hs_live_rows(pr)) / pr$oldNumRow
  c <- 100 * (pr$oldNumCol - .hs_live_cols(pr)) / pr$oldNumCol
  # std::max(r, c) is (r < c) ? c : r, so a NaN row reduction wins
  if (!is.nan(r) && !is.nan(c) && r < c) c else r
}
.hs_fast_presolve_loop <- function(pr) {
  repeat {
    .hs_store_size(pr)
    .hs_remove_row_singletons(pr)
    .hs_presolve_changed_rows(pr)
    .hs_remove_doubleton_equations(pr)
    .hs_presolve_col_singletons(pr)
    .hs_presolve_changed_cols(pr)
    if (!isTRUE(.hs_size_reduction(pr) > 0.01)) break
  }
}

# shrinkProblem followed by toCSC and okFromCSC: the matrix is relinked in
# column order with the slots of each column in slot order
.hs_shrink <- function(pr) {
  live <- which(pr$Avalue != 0)
  cols <- which(!pr$colDeleted); rows <- which(!pr$rowDeleted)
  pr$ncol_now <- length(cols); pr$nrow_now <- length(rows)
  pr$numDeletedCols <- 0L; pr$numDeletedRows <- 0L
  ord <- live[order(pr$Acol[live], live)]
  pr$Avalue <- pr$Avalue[ord]; pr$Arow <- pr$Arow[ord]; pr$Acol <- pr$Acol[ord]
  k <- length(ord)
  pr$Anext <- pr$Aprev <- pr$ARleft <- pr$ARright <- integer(k)
  pr$freeslots <- integer()
  pr$colhead[] <- 0L; pr$colsize[] <- 0L; pr$rowroot[] <- 0L; pr$rowsize[] <- 0L
  pr$lsb <- list(r = .hs_lsb_new(pr$m), d = .hs_lsb_new(pr$n))
  for (pos in seq_len(k)) .hs_link(pr, pos)
  pr$rcur[rows] <- seq_along(rows) - 1L; pr$ccur[cols] <- seq_along(cols) - 1L
  pr$eqsize <- ifelse(!pr$rowDeleted & pr$rl == pr$ru, pr$rowsize, NA_integer_)
  sb <- pr$subst
  pr$subst <- sb[!pr$rowDeleted[sb[, 1L]] & !pr$colDeleted[sb[, 2L]], , drop = FALSE]
}

.hs_presolve_main <- function(pr) {
  for (row in seq_len(pr$m)) {
    if (pr$rowDeleted[row]) next
    .hs_row_presolve(pr, row)
    pr$changedRowFlag[row] <- FALSE
  }
  for (col in seq_len(pr$n)) {
    if (pr$colDeleted[col]) next
    .hs_col_presolve(pr, col)
    pr$changedColFlag[col] <- FALSE
  }
  parallel_calls <- 0L
  dependent_called <- FALSE
  big <- function() isTRUE(.hs_size_reduction(pr) > 0.05)
  maybe_shrink <- function() {
    if (pr$numDeletedCols >= pr$ncol_now %/% 2L || pr$numDeletedRows >= pr$nrow_now %/% 2L) .hs_shrink(pr)
  }
  # sparsify is off for LPs (lp_presolve_requires_basis_postsolve)
  repeat {
    .hs_fast_presolve_loop(pr)
    .hs_store_size(pr)
    .hs_aggregator(pr)
    if (big()) next
    if (parallel_calls < 5L) {
      maybe_shrink()
      .hs_store_size(pr)
      .hs_parallel_rows_cols(pr)
      parallel_calls <- parallel_calls + 1L
      if (big()) next
    }
    .hs_fast_presolve_loop(pr)
    .hs_fast_presolve_loop(pr)
    if (!dependent_called) {
      maybe_shrink()
      .hs_store_size(pr)
      .hs_dependent_equations(pr)
      dependent_called <- TRUE
      if (big()) next
    }
    break
  }
}

.hs_push_subst <- function(pr, row, col) pr$subst <- rbind(pr$subst, c(row, col))

.hs_max_abs <- function(pr, slots) if (length(slots)) max(abs(pr$Avalue[slots])) else 0

# HighsHashHelpers::hash of the pair (row, col) of 32 bit current indices
.hs_pair_key_hash <- function(r, c) {
  .hs_u64_xor(.hs_pair_hash(1L, r, c), .hs_u64_shift(.hs_pair_hash(0L, r, c), 32L))
}

# HPresolve::aggregator
.hs_aggregator <- function(pr) {
  ok <- function(row, col) !pr$rowDeleted[row] && !pr$colDeleted[col] &&
    .hs_is_implied_free(pr, col) && .hs_is_dual_implied_free(pr, row)
  sb <- pr$subst
  keep <- vapply(seq_len(nrow(sb)), function(k) ok(sb[k, 1L], sb[k, 2L]), TRUE)
  sb <- sb[keep, , drop = FALSE]
  if (nrow(sb) > 1L) {
    rs <- pr$rowsize[sb[, 1L]]; cs <- pr$colsize[sb[, 2L]]
    minlen <- pmin(rs, cs)
    rc <- pr$rcur[sb[, 1L]]; cc <- pr$ccur[sb[, 2L]]
    h <- vapply(seq_len(nrow(sb)), function(k) .hs_pair_key_hash(rc[k], cc[k]), numeric(4L))
    o <- order(minlen != 2L, rs * cs, minlen, h[4L, ], h[3L, ], h[2L, ], h[1L, ], rc, cc)
    sb <- sb[o, , drop = FALSE]
  }
  pr$subst <- sb
  # the loop also visits pairs pushed while it runs; row 0 marks a pair done
  nfail <- 0L
  k <- 0L
  while (k < nrow(pr$subst)) {
    k <- k + 1L
    row <- pr$subst[k, 1L]; col <- pr$subst[k, 2L]
    if (row == 0L) next
    if (!ok(row, col)) { pr$subst[k, 1L] <- 0L; next }
    pos <- .hs_find_nonzero(pr, row, col)
    if (pos == 0L) { pr$subst[k, 1L] <- 0L; next }
    if (pr$rowsize[row] == 2L || pr$colsize[col] == 2L) {
      .hs_substitute_free_col(pr, row, col, .hs_row_inorder(pr, row), TRUE)
      pr$subst[k, 1L] <- 0L
      .hs_remove_row_singletons(pr)
      next
    }
    row_max <- function() .hs_max_abs(pr, .hs_row_preorder(pr, row))
    col_max <- function() .hs_max_abs(pr, .hs_col_slots(pr, col))
    shorter_row <- pr$rowsize[row] < pr$colsize[col]
    maxv <- if (shorter_row) row_max() else col_max()
    if (abs(pr$Avalue[pos]) < maxv * 0.01) {
      maxv <- if (shorter_row) col_max() else row_max()
      if (abs(pr$Avalue[pos]) < maxv * 0.01) { pr$subst[k, 1L] <- 0L; next }
    }
    rowpos <- .hs_row_inorder(pr, row)
    fillin <- -(pr$rowsize[row] + pr$colsize[col] - 1L)
    for (p in .hs_col_slots(pr, col)) {
      other <- pr$Arow[p]
      if (other == row) next
      for (q in rowpos) if (.hs_find_nonzero(pr, other, pr$Acol[q]) == 0L) fillin <- fillin + 1L
      if (fillin > 10L) break
    }
    if (fillin > 10L) {
      nfail <- nfail + 1L
      if (nfail == 3L) break
      next
    }
    nfail <- 0L
    .hs_substitute_free_col(pr, row, col, rowpos, TRUE)
    pr$subst[k, 1L] <- 0L
    .hs_remove_row_singletons(pr)
  }
  pr$subst <- pr$subst[pr$subst[, 1L] != 0L, , drop = FALSE]
}

# HighsHashHelpers::double_hash_code: the binary exponent of val times the
# reciprocal golden ratio in the upper 16 bits, its sign and 15 leading
# mantissa bits in the lower 16 bits
.hs_double_hash_code <- function(val) {
  v <- val * 0x1.3c6ef372fe94fp-1
  e <- floor(log2(abs(v))) + 1
  m <- v / 2^e
  if (abs(m) >= 1) {
    e <- e + 1; m <- m / 2
  } else if (abs(m) < 0.5) {
    e <- e - 1; m <- m * 2
  }
  (e %% 65536) * 65536 + trunc(m * 32768) %% 65536
}

# libc++ unordered_multimap buckets hold the equal key group in order;
# emplace_hint goes before the last element visited
.hs_bucket_group <- function(b, key) if (is.null(b[[key]])) integer() else b[[key]]
.hs_bucket_insert <- function(grp, i) {
  n <- length(grp)
  if (n) c(grp[-n], i, grp[n]) else i
}

# HPresolve::detectParallelRowsAndCols for an LP
.hs_parallel_rows_cols <- function(pr) {
  small <- 1e-9
  live <- which(pr$Avalue != 0)
  row_max <- numeric(pr$m); row_max_i <- integer(pr$m)
  col_max <- numeric(pr$n); col_max_i <- integer(pr$n)
  row_singletons <- integer(pr$m)
  for (i in live) {
    r <- pr$Arow[i]; c <- pr$Acol[i]; v <- pr$Avalue[i]
    if (pr$colsize[c] == 1L) {
      col_max[c] <- v
      row_singletons[r] <- row_singletons[r] + 1L
      next
    }
    a <- abs(v); ar <- abs(row_max[r])
    if (a >= ar - small && (a > ar + small || pr$ccur[c] < row_max_i[r])) {
      row_max[r] <- v; row_max_i[r] <- pr$ccur[c]
    }
    ac <- abs(col_max[c])
    if (a >= ac - small && (a > ac + small || pr$rcur[r] < col_max_i[c])) {
      col_max[c] <- v; col_max_i[c] <- pr$rcur[r]
    }
  }
  row_hash <- lapply(pr$rowsize - row_singletons, function(v) .hs_u64(0, v))
  col_hash <- lapply(pr$colsize, function(v) .hs_u64(0, v))
  for (i in live) {
    r <- pr$Arow[i]; c <- pr$Acol[i]; v <- pr$Avalue[i]
    if (pr$colsize[c] == 1L) {
      col_hash[[c]] <- .hs_u64(0, pr$rcur[r])
    } else {
      row_hash[[r]] <- .hs_sparse_combine(row_hash[[r]], pr$ccur[c], .hs_double_hash_code(v / row_max[r]))
      col_hash[[c]] <- .hs_sparse_combine(col_hash[[c]], pr$rcur[r], .hs_double_hash_code(v / col_max[c]))
    }
  }
  row_key <- vapply(row_hash, paste, "", collapse = ",")
  col_key <- vapply(col_hash, paste, "", collapse = ",")
  dt <- pr$dtol
  buckets <- new.env(parent = emptyenv())
  for (i in seq_len(pr$n)) {
    if (pr$colDeleted[i]) next
    if (pr$colsize[i] == 0L) {
      .hs_col_presolve(pr, i)
      next
    }
    grp <- .hs_bucket_group(buckets, col_key[i])
    del <- 0L; last <- 0L
    for (gi in seq_along(grp)) {
      cand <- grp[gi]; last <- gi
      if (pr$colsize[i] != pr$colsize[cand]) next
      col <- cand; dup <- i
      scale <- col_max[dup] / col_max[col]
      upper_inf <- function() {
        if (scale > 0) .hs_is_upper_strictly_implied(pr, col, pr$ptol) else .hs_is_lower_strictly_implied(pr, col, pr$ptol)
      }
      lower_inf <- function() {
        if (scale > 0) .hs_is_lower_strictly_implied(pr, col, pr$ptol) else .hs_is_upper_strictly_implied(pr, col, pr$ptol)
      }
      dup_upper_inf <- function() .hs_is_upper_strictly_implied(pr, dup, pr$ptol)
      dup_lower_inf <- function() .hs_is_lower_strictly_implied(pr, dup, pr$ptol)
      obj_diff <- .cd_val(.cd_add(.cd_two_product(pr$cost[col], scale), -pr$cost[dup]))
      cl <- pr$cl; cu <- pr$cu
      to_upper <- if (scale > 0) "colupper" else "collower"
      to_lower <- if (scale > 0) "collower" else "colupper"
      case <- "merge"
      if (obj_diff < -dt) {
        if (upper_inf() && cl[dup] != -Inf) case <- "duplower"
        else if (dup_lower_inf() && (scale < 0 || cu[col] != Inf) && (scale > 0 || cl[col] != -Inf)) case <- to_upper
        else next
      } else if (obj_diff > dt) {
        if (lower_inf() && cu[dup] != Inf) case <- "dupupper"
        else if (dup_upper_inf() && (scale < 0 || cl[col] != -Inf) && (scale > 0 || cu[col] != Inf)) case <- to_lower
        else next
      } else {
        if (upper_inf() && cl[dup] != -Inf) case <- "duplower"
        else if (lower_inf() && cu[dup] != Inf) case <- "dupupper"
        else if (dup_upper_inf() && (scale < 0 || cl[col] != -Inf) && (scale > 0 || cu[col] != Inf)) case <- to_lower
        else if (dup_lower_inf() && (scale < 0 || cu[col] != Inf) && (scale > 0 || cl[col] != -Inf)) case <- to_upper
      }
      parallel <- TRUE
      for (p in .hs_col_slots(pr, col)) {
        q <- .hs_find_nonzero(pr, pr$Arow[p], dup)
        parallel <- q != 0L && abs(.cd_val(.cd_add(-.cd_two_product(scale, pr$Avalue[p]), pr$Avalue[q]))) <= small
        if (!parallel) break
      }
      if (!parallel) next
      switch(case,
        duplower = { del <- dup; .hs_fix_col_to(pr, dup, FALSE) },
        dupupper = { del <- dup; .hs_fix_col_to(pr, dup, TRUE) },
        collower = { del <- col; .hs_fix_col_to(pr, col, FALSE) },
        colupper = { del <- col; .hs_fix_col_to(pr, col, TRUE) },
        merge = {
          .hs_ps_push(pr, list(type = "duplicatecol", col = col, dup = dup, scale = scale,
                               cl = cl[col], cu = cu[col], dl = cl[dup], du = cu[dup]))
          .hs_mark_col(pr, col)
          if (scale > 0) {
            ml <- cl[col] + scale * cl[dup]; mu <- cu[col] + scale * cu[dup]
          } else {
            ml <- cl[col] + scale * cu[dup]; mu <- cu[col] + scale * cl[dup]
          }
          .hs_change_col_lower(pr, col, ml)
          .hs_change_col_upper(pr, col, mu)
          .hs_mark_col_deleted(pr, dup)
          .hs_shift_rows_for_col(pr, dup, NA_real_)
          pr$cost[dup] <- 0
          del <- dup
          .hs_reset_col_impl(pr, col, 0L)
        })
      break
    }
    if (del != 0L) {
      if (del != i) buckets[[col_key[i]]] <- grp[-last]
      .hs_remove_row_singletons(pr)
    } else {
      buckets[[col_key[i]]] <- .hs_bucket_insert(grp, i)
    }
  }
  buckets <- new.env(parent = emptyenv())
  for (i in seq_len(pr$m)) {
    if (pr$rowDeleted[i]) next
    if (pr$rowsize[i] <= 1L || (pr$rowsize[i] == 2L && .hs_is_equation(pr, i))) {
      .hs_row_presolve(pr, i)
      next
    }
    # lp_presolve_requires_basis_postsolve skips rows with singleton
    # columns, so only the duplicate row reduction applies
    if (row_singletons[i] != 0L) next
    grp <- .hs_bucket_group(buckets, row_key[i])
    del <- 0L; last <- 0L
    rowpos <- if (length(grp)) .hs_row_inorder(pr, i) else integer()
    for (gi in seq_along(grp)) {
      cand <- grp[gi]; last <- gi
      if (pr$rowsize[i] != pr$rowsize[cand]) next
      scale <- row_max[cand] / row_max[i]
      parallel <- TRUE
      for (p in rowpos) {
        if (pr$colsize[pr$Acol[p]] == 1L) next
        q <- .hs_find_nonzero(pr, cand, pr$Acol[p])
        parallel <- q != 0L && abs(.cd_val(.cd_add(-.cd_two_product(scale, pr$Avalue[p]), pr$Avalue[q]))) <= small
        if (!parallel) break
      }
      if (!parallel) next
      ut <- FALSE; lt <- FALSE
      if (scale > 0) {
        nu <- pr$ru[i] * scale; nl <- pr$rl[i] * scale
      } else {
        nl <- pr$ru[i] * scale; nu <- pr$rl[i] * scale
      }
      if (nu < pr$ru[cand]) {
        if (nu < pr$rl[cand] - pr$ptol) .hs_pre_abort("infeasible")
        if (nu <= pr$rl[cand] + pr$ptol) nu <- pr$rl[cand]
        if (nu < pr$ru[cand]) {
          ut <- TRUE
          if (scale > 0) {
            tmp <- pr$rdl[i] / scale; pr$rdl[i] <- pr$rdl[cand] * scale; pr$rdl[cand] <- tmp
          } else {
            tmp <- pr$rdu[i] / scale; pr$rdu[i] <- pr$rdl[cand] * scale; pr$rdl[cand] <- tmp
          }
          pr$ru[cand] <- nu
        }
      }
      if (nl > pr$rl[cand]) {
        if (nl > pr$ru[cand] + pr$ptol) .hs_pre_abort("infeasible")
        if (nl >= pr$ru[cand] - pr$ptol) nl <- pr$ru[cand]
        if (nl > pr$rl[cand]) {
          lt <- TRUE
          if (scale > 0) {
            tmp <- pr$rdu[i] / scale; pr$rdu[i] <- pr$rdu[cand] * scale; pr$rdu[cand] <- tmp
          } else {
            tmp <- pr$rdl[i] / scale; pr$rdl[i] <- pr$rdu[cand] * scale; pr$rdu[cand] <- tmp
          }
          pr$rl[cand] <- nl
        }
      }
      .hs_reset_row_dual_impl(pr, cand, 0L)
      .hs_ps_push(pr, list(type = "duplicaterow", row = cand, dup = i, scale = scale, ut = ut, lt = lt))
      del <- i
      .hs_mark_row_deleted(pr, i)
      for (p in rowpos) .hs_unlink(pr, p)
      break
    }
    if (del != 0L) {
      if (del != i) buckets[[row_key[i]]] <- grp[-last]
    } else {
      buckets[[row_key[i]]] <- .hs_bucket_insert(grp, i)
    }
  }
}

# HPresolve::removeDependentEquations: HFactor::build on the equations as
# columns (their nonzeros in row order and the right hand side as an extra
# row) leaves the dependent ones without a pivot
.hs_dependent_equations <- function(pr) {
  eqs <- integer()
  row <- .hs_next_equation(pr)
  while (row != 0L) {
    eqs <- c(eqs, row)
    row <- .hs_next_equation(pr, pr$eqsize[row] * (pr$m + 1) + row)
  }
  if (!length(eqs)) return(invisible())
  # the columns of the factor matrix are numbered by their current index
  ncur <- pr$ncol_now
  cols <- lapply(eqs, function(eq) {
    s <- .hs_row_preorder(pr, eq)
    idx <- pr$ccur[pr$Acol[s]] + 1L; val <- pr$Avalue[s]
    if (pr$rl[eq] != 0) {
      idx <- c(idx, ncur + 1L); val <- c(val, pr$rl[eq])
    }
    list(index = idx, value = val)
  })
  f <- .hs_factor_build(ncur + 1L, cols)
  for (k in seq_len(f$rank_deficiency)) {
    eq <- eqs[f$col_with_no_pivot[k]]
    .hs_ps_redundant_row(pr, eq)
    .hs_remove_row(pr, eq)
  }
}

# HPresolve::run on the LP; returns the status ("reduced", "empty",
# "notreduced", "infeasible", "unboundedorinfeasible") and the reduced LP
.hs_presolve <- function(cost, A, row_lower, row_upper, col_lower, col_upper, ptol = 1e-7, dtol = 1e-7) {
  pr <- new.env(parent = emptyenv())
  m <- nrow(A); n <- ncol(A)
  pr$m <- m; pr$n <- n; pr$ncol_now <- n; pr$nrow_now <- m
  pr$ptol <- ptol; pr$dtol <- dtol
  pr$cost <- cost; pr$cl <- col_lower; pr$cu <- col_upper; pr$rl <- row_lower; pr$ru <- row_upper
  pr$offset <- 0
  pr$colLowerSource <- pr$colUpperSource <- integer(n)
  pr$implColLower <- rep(-Inf, n); pr$implColUpper <- rep(Inf, n)
  pr$colImplSourceByRow <- rep(list(integer()), m)
  pr$implRowDualSourceByCol <- rep(list(integer()), n)
  pr$rdl <- ifelse(row_upper == Inf, 0, -Inf); pr$rdu <- ifelse(row_lower == -Inf, 0, Inf)
  pr$implRowDualLower <- rep(-Inf, m); pr$implRowDualUpper <- rep(Inf, m)
  pr$rowDualLowerSource <- pr$rowDualUpperSource <- integer(m)
  nz <- which(A != 0)
  pr$Avalue <- A[nz]; pr$Arow <- as.integer((nz - 1L) %% m + 1L); pr$Acol <- as.integer((nz - 1L) %/% m + 1L)
  k <- length(nz)
  pr$Anext <- pr$Aprev <- pr$ARleft <- pr$ARright <- integer(k)
  pr$freeslots <- integer()
  pr$colhead <- integer(n); pr$colsize <- integer(n); pr$rowroot <- integer(m); pr$rowsize <- integer(m)
  pr$lsb <- list(r = .hs_lsb_new(m), d = .hs_lsb_new(n))
  pr$rowDeleted <- logical(m); pr$colDeleted <- logical(n)
  pr$changedRowFlag <- rep(TRUE, m); pr$changedColFlag <- rep(TRUE, n)
  pr$changedRowIndices <- integer(); pr$changedColIndices <- integer()
  pr$singletonRows <- integer(); pr$singletonColumns <- integer()
  pr$numDeletedRows <- 0L; pr$numDeletedCols <- 0L
  pr$stack <- list()
  pr$rcur <- seq_len(m) - 1L; pr$ccur <- seq_len(n) - 1L
  pr$subst <- matrix(0L, 0L, 2L)
  for (pos in seq_len(k)) .hs_link(pr, pos)
  pr$eqsize <- ifelse(row_lower == row_upper, pr$rowsize, NA_integer_)
  for (row in which(row_lower == row_upper))
    for (p in .hs_row_preorder(pr, row)) if (.hs_is_implied_free(pr, pr$Acol[p])) .hs_push_subst(pr, row, pr$Acol[p])
  status <- tryCatch({
    .hs_presolve_main(pr)
    NULL
  }, hs_presolve_abort = function(e) e$result)
  if (identical(status, "infeasible")) return(list(status = "infeasible", pr = pr))
  if (identical(status, "dualinfeasible")) return(list(status = "unboundedorinfeasible", pr = pr))
  cols <- which(!pr$colDeleted); rows <- which(!pr$rowDeleted)
  st <- if (!length(cols)) "empty" else if (length(pr$stack)) "reduced" else "notreduced"
  Ar <- matrix(0, length(rows), length(cols))
  live <- which(pr$Avalue != 0)
  Ar[cbind(match(pr$Arow[live], rows), match(pr$Acol[live], cols))] <- pr$Avalue[live]
  list(status = st, pr = pr, cols = cols, rows = rows,
       lp = list(cost = pr$cost[cols], A = Ar, row_lower = pr$rl[rows], row_upper = pr$ru[rows],
                 col_lower = pr$cl[cols], col_upper = pr$cu[cols], offset = pr$offset))
}

# HighsPostsolveStack::undo
.hs_postsolve <- function(pre, sol, dtol = 1e-7) {
  pr <- pre$pr
  expand <- function(v, idx, size) {
    out <- c(v, rep(if (is.character(v)) "" else 0, max(0L, size - length(v))))
    for (i in rev(seq_along(idx))) out[idx[i]] <- out[i]
    out
  }
  x <- expand(sol$col_value, pre$cols, pr$n); rv <- expand(sol$row_value, pre$rows, pr$m)
  cd <- expand(sol$col_dual, pre$cols, pr$n); rd <- expand(sol$row_dual, pre$rows, pr$m)
  cs <- expand(sol$col_status, pre$cols, pr$n); rs <- expand(sol$row_status, pre$rows, pr$m)
  row_status_of <- function(dual, rt) if (rt == "eq") { if (dual < 0) "upper" else "lower" } else if (rt == "geq") "lower" else "upper"
  col_status_now <- function(j) {
    st <- if (cd[j] > dtol) "lower" else if (cd[j] < -dtol) "upper" else cs[j]
    cs[j] <<- st
    st
  }
  col_status_row <- function(i) {
    st <- if (rd[i] > dtol) "lower" else if (rd[i] < -dtol) "upper" else rs[i]
    rs[i] <<- st
    st
  }
  for (r in rev(pr$stack)) {
    switch(r$type,
      fixedcol = {
        x[r$col] <- r$fix
        red <- c(r$cost, 0)
        for (k in seq_along(r$colv$index)) red <- .cd_add(red, -(r$colv$value[k] * rd[r$colv$index[k]]))
        cd[r$col] <- .cd_val(red)
        cs[r$col] <- if (r$status == "nonbasic") { if (cd[r$col] >= 0) "lower" else "upper" } else r$status
      },
      redundantrow = { rd[r$row] <- 0; rs[r$row] <- "basic" },
      singletonrow = {
        st <- col_status_now(r$col)
        if ((!r$lt || st != "lower") && (!r$ut || st != "upper")) {
          rs[r$row] <- "basic"; rd[r$row] <- 0
        } else {
          rd[r$row] <- cd[r$col] / r$coef
          cd[r$col] <- 0
          if (st == "lower") rs[r$row] <- if (r$coef > 0) "lower" else "upper"
          else if (st == "upper") rs[r$row] <- if (r$coef > 0) "upper" else "lower"
          cs[r$col] <- "basic"
        }
      },
      forcingrow = {
        basic <- 0L; delta <- 0
        dir <- if (r$rowtype == "leq") 1 else -1
        for (k in seq_along(r$rowv$index)) {
          j <- r$rowv$index[k]; a <- r$rowv$value[k]
          if (dir * (cd[j] - a * delta) * a < 0) { delta <- cd[j] / a; basic <- j }
        }
        if (basic != 0L) {
          rd[r$row] <- rd[r$row] + delta
          for (k in seq_along(r$rowv$index)) {
            j <- r$rowv$index[k]
            cd[j] <- .cd_val(.cd_add(-.cd_mul(c(delta, 0), r$rowv$value[k]), cd[j]))
          }
          cd[basic] <- 0
          rs[r$row] <- if (r$rowtype == "geq") "lower" else "upper"
          cs[basic] <- "basic"
        }
      },
      doubleton = {
        x[r$colsubst] <- .cd_val(.cd_div(.cd_add(-.cd_mul(c(r$coef, 0), x[r$col]), r$rhs), r$coefsubst))
        if (r$row == 0L) next
        st <- col_status_now(r$col)
        rd[r$row] <- 0
        rdual <- c(0, 0)
        for (k in seq_along(r$colv$index)) rdual <- .cd_add(rdual, -(r$colv$value[k] * rd[r$colv$index[k]]))
        rdual <- .cd_div(rdual, r$coefsubst)
        rd[r$row] <- .cd_val(rdual)
        cd[r$colsubst] <- r$substcost
        cd[r$col] <- cd[r$col] + r$substcost * r$coef / r$coefsubst
        if ((r$ut && st == "upper") || (r$lt && st == "lower")) {
          delta <- cd[r$col] / r$coef
          rd[r$row] <- .cd_val(.cd_add(rdual, delta))
          cd[r$col] <- 0
          cd[r$colsubst] <- .cd_val(.cd_add(c(cd[r$colsubst], 0), -(delta * r$coefsubst)))
          cs[r$colsubst] <- if ((sign(r$coef) == sign(r$coefsubst) && cs[r$col] == "upper") ||
                                (sign(r$coef) != sign(r$coefsubst) && cs[r$col] == "lower")) "lower" else "upper"
          cs[r$col] <- "basic"
        } else {
          delta <- cd[r$colsubst] / r$coefsubst
          rd[r$row] <- .cd_val(.cd_add(rdual, delta))
          cd[r$colsubst] <- 0
          cd[r$col] <- .cd_val(.cd_add(c(cd[r$col], 0), -(delta * r$coef)))
          cs[r$colsubst] <- "basic"
        }
        rs[r$row] <- row_status_of(rd[r$row], r$rowtype)
      },
      freecol = {
        coef <- 0; val <- c(0, 0)
        for (k in seq_along(r$rowv$index)) {
          j <- r$rowv$index[k]
          if (j == r$col) coef <- r$rowv$value[k] else val <- .cd_add(val, r$rowv$value[k] * x[j])
        }
        rv[r$row] <- .cd_val(.cd_add(val, coef * x[r$col]))
        x[r$col] <- .cd_val(.cd_div(.cd_add(-val, r$rhs), coef))
        rd[r$row] <- 0
        dv <- c(r$cost, 0)
        for (k in seq_along(r$colv$index)) dv <- .cd_add(dv, -(r$colv$value[k] * rd[r$colv$index[k]]))
        rd[r$row] <- .cd_val(.cd_div(dv, coef))
        cd[r$col] <- 0
        cs[r$col] <- "basic"
        rs[r$row] <- row_status_of(rd[r$row], r$rowtype)
      },
      forcingcol = {
        nb_row <- 0L; nb_status <- ""; val <- r$bound
        dir <- if (r$at_inf_upper) 1 else -1
        for (k in seq_along(r$colv$index)) {
          i <- r$colv$index[k]; a <- r$colv$value[k]
          v <- rv[i] / a
          if (dir * v > dir * val) { nb_row <- i; val <- v; nb_status <- if (dir * a > 0) "lower" else "upper" }
        }
        x[r$col] <- val; cd[r$col] <- 0
        if (nb_row == 0L) cs[r$col] <- if (r$at_inf_upper) "lower" else "upper"
        else { cs[r$col] <- "basic"; rs[nb_row] <- nb_status }
      },
      duplicaterow = {
        if (!r$ut && !r$lt) {
          rd[r$dup] <- 0; rs[r$dup] <- "basic"
        } else {
          st <- col_status_row(r$row)
          tight <- if (st == "upper") r$ut else if (st == "lower") r$lt else NA
          if (isTRUE(tight)) {
            rd[r$dup] <- rd[r$row] / r$scale
            rs[r$dup] <- if (r$scale > 0) "upper" else "lower"
            rd[r$row] <- 0; rs[r$row] <- "basic"
          } else {
            rd[r$dup] <- 0; rs[r$dup] <- "basic"
          }
        }
      },
      duplicatecol = {
        # ponytail: DuplicateColumn::undoFix only runs when the split values
        # break their bounds, which needs a merged value outside its bounds
        merge <- x[r$col]
        cd[r$dup] <- cd[r$col] * r$scale
        st <- cs[r$col]
        if (st == "lower" || st == "upper") {
          x[r$col] <- if (st == "lower") r$cl else r$cu
          at_lower <- (st == "lower") == (r$scale > 0)
          cs[r$dup] <- if (at_lower) "lower" else "upper"
          x[r$dup] <- if (at_lower) r$dl else r$du
        } else if (st == "zero") {
          x[r$col] <- 0; x[r$dup] <- 0; cs[r$dup] <- "zero"
        } else {
          x[r$col] <- if (r$cl != -Inf) r$cl else min(0, r$cu)
          x[r$dup] <- .cd_val(.cd_div(.cd_add(c(merge, 0), -x[r$col]), r$scale))
          recompute <- FALSE
          if (x[r$dup] > r$du) {
            x[r$dup] <- r$du; recompute <- TRUE; cs[r$dup] <- "upper"
          } else if (x[r$dup] < r$dl) {
            x[r$dup] <- r$dl; recompute <- TRUE; cs[r$dup] <- "lower"
          }
          if (recompute) {
            x[r$col] <- merge - r$scale * x[r$dup]
          } else {
            cs[r$dup] <- cs[r$col]; cs[r$col] <- "lower"
          }
        }
      },
      forcingcolrow = {
        val <- c(r$rhs, 0)
        for (k in seq_along(r$rowv$index)) val <- .cd_add(val, -(r$rowv$value[k] * x[r$rowv$index[k]]))
        rv[r$row] <- .cd_val(val); rd[r$row] <- 0; rs[r$row] <- "basic"
      })
  }
  list(col_value = x, row_value = rv, col_dual = cd, row_dual = rd, col_status = cs, row_status = rs)
}

# Highs::run for an LP: presolve, the simplex solve of the reduced LP,
# postsolve and the cleanup solve of the original LP from the postsolved
# basis
.hs_highs_run <- function(cost, A, row_lower, row_upper, col_lower, col_upper, opts = list()) {
  use_presolve <- !isFALSE(opts$presolve)
  if (!use_presolve || nrow(A) == 0L || all(A == 0))
    return(.hs_lp_solve(cost, A, row_lower, row_upper, col_lower, col_upper, opts))
  pre <- .hs_presolve(cost, A, row_lower, row_upper, col_lower, col_upper,
                      opts$ptol %||% 1e-7, opts$dtol %||% 1e-7)
  fail <- function(status) list(status = status, nit = 0L, primal_feasible = FALSE)
  switch(pre$status,
    notreduced = return(.hs_lp_solve(cost, A, row_lower, row_upper, col_lower, col_upper, opts)),
    infeasible = return(fail("Infeasible")),
    unboundedorinfeasible = {
      opts$strategy <- "primal"
      return(.hs_lp_solve(cost, A, row_lower, row_upper, col_lower, col_upper, opts))
    })
  if (pre$status == "empty") {
    red <- list(status = "Optimal", nit = 0L, col_value = numeric(), col_dual = numeric(),
                row_value = numeric(), row_dual = numeric(), col_status = character(), row_status = character())
  } else {
    lp <- pre$lp
    red <- .hs_lp_solve(lp$cost, lp$A, lp$row_lower, lp$row_upper, lp$col_lower, lp$col_upper, opts)
    if (is.null(red$row_status)) red$row_status <- character(length(pre$rows))
    if (!red$status %in% c("Optimal", "Unknown")) return(red)
  }
  post <- .hs_postsolve(pre, red, opts$dtol %||% 1e-7)
  opts$basis <- list(col_status = post$col_status, row_status = post$row_status)
  opts$iteration_start <- red$nit
  .hs_lp_solve(cost, A, row_lower, row_upper, col_lower, col_upper, opts)
}
