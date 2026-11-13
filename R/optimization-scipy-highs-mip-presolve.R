# The MIP rules of the HiGHS 1.12.0 presolve (presolve/HPresolve.cpp: the
# integer branches of the row and column rules, implied integers, coefficient
# strengthening, sparsify, strengthenInequalities, dominatedColumns,
# probing, scaleMIP and the mipsolver hooks) and the primal postsolve of a
# MIP (presolve/HighsPostsolveStack.cpp), as SciPy 1.18.1 vendors them. MIT
# license notice in R/optimization-scipy-highs-util.R. The presolve state
# `pr` is the environment of R/optimization-scipy-highs-presolve.R: indices
# are 1-based and 0 stands for HiGHS's -1. Calls into the clique table, the
# implications and the domain (group A) convert to their 0-based indices.


# ---------------------------------------------------------------------------
# HighsIntegers::mod and modularInverse (util/HighsIntegers.h); the rest of
# HighsIntegers and the HighsCDouble extras live in -mip-util.R (group B1)
# ---------------------------------------------------------------------------

# int64 HighsIntegers::mod
.hmip_mod_int <- function(a, m) {
  r <- a - m * trunc(a / m)
  r + (r < 0) * m
}
# double HighsIntegers::mod: int64_t(std::fmod(a, m)) + (a < 0) * m
.hmip_mod <- function(a, m) trunc(sign(a) * (abs(a) %% abs(m))) + (a < 0) * m
.hmip_modularInverse <- function(a, m) {
  y <- 0; x <- 1
  if (m == 1) return(0)
  a <- .hmip_mod_int(a, m)
  while (a > 1) {
    q <- trunc(a / m); r <- a - q * m
    a <- m; m <- r
    r <- x - q * y; x <- y; y <- r
  }
  x
}

# ---------------------------------------------------------------------------
# HighsLinearSumBounds pieces used only by the MIP rules
# ---------------------------------------------------------------------------

# getImplVarLower/Upper(sum, var)
.hs_lsb_impl_var_lower <- function(pr, sum, col) .hs_lsb_impl_lo(sum, pr$cl[col], pr$implColLower[col], pr$colLowerSource[col])
.hs_lsb_impl_var_upper <- function(pr, sum, col) .hs_lsb_impl_up(sum, pr$cu[col], pr$implColUpper[col], pr$colUpperSource[col])
# sumScaled
.hs_lsb_sum_scaled <- function(pr, kind, sum, scale) {
  b <- pr$lsb[[kind]]
  for (w in c("L", "U", "LO", "UO")) {
    v <- .cd_mul(c(b[[paste0(w, "h")]][sum], b[[paste0(w, "l")]][sum]), scale)
    b[[paste0(w, "h")]][sum] <- v[1L]; b[[paste0(w, "l")]][sum] <- v[2L]
  }
  if (scale < 0) {
    for (p in list(c("L", "U"), c("LO", "UO"))) for (f in c("h", "l", "n")) {
      a1 <- paste0(p[1L], f); a2 <- paste0(p[2L], f)
      t <- b[[a1]][sum]; b[[a1]][sum] <- b[[a2]][sum]; b[[a2]][sum] <- t
    }
  }
}
# getSum{Lower,Upper}{,Orig}(sum, offset) with a HighsCDouble offset
.hs_lsb_sum_cd <- function(pr, kind, w, sum, offset) {
  b <- pr$lsb[[kind]]
  if (b[[paste0(w, "n")]][sum] != 0L) return(if (w %in% c("L", "LO")) -Inf else Inf)
  .cd_val(.cd_add_cd(c(b[[paste0(w, "h")]][sum], b[[paste0(w, "l")]][sum]), offset))
}
# getResidualSum{Lower,Upper}{,Orig}(sum, var, coef, boundVar, boundCoef,
# boundVal) of the row sums
.hs_lsb_residual2 <- function(pr, w, sum, var, coef, bvar, bcoef, bval) {
  b <- pr$lsb$r
  act <- c(b[[paste0(w, "h")]][sum], b[[paste0(w, "l")]][sum]); ninf <- b[[paste0(w, "n")]][sum]
  lower_side <- w %in% c("L", "LO"); orig <- w %in% c("LO", "UO")
  bound_of <- function(v, cf) {
    use_lower <- if (lower_side) cf > 0 else cf < 0
    if (orig) { if (use_lower) pr$cl[v] else pr$cu[v] }
    else if (use_lower) .hs_lsb_impl_var_lower(pr, sum, v) else .hs_lsb_impl_var_upper(pr, sum, v)
  }
  upd <- function(dir, bound, cf) {
    if (abs(bound) != Inf) act <<- .cd_iadd(act, .cd_two_product(dir * bound, cf)) else ninf <<- ninf + dir
  }
  upd(-1, bound_of(var, coef), coef)
  if (bvar != 0L) {
    old <- bound_of(bvar, bcoef)
    if (old != bval) { upd(-1, old, bcoef); upd(1, bval, bcoef) }
  }
  if (ninf == 0L) .cd_val(act) else if (lower_side) -Inf else Inf
}

# HPresolve::computeColBounds; want is a subset of "lower", "upper",
# "wclower", "wcupper"
.hs_compute_col_bounds <- function(pr, col, bcol = 0L, bval = Inf, pattern = 0L, want) {
  lreq <- any(c("lower", "wclower") %in% want); ureq <- any(c("upper", "wcupper") %in% want)
  skip <- function(row, val) {
    hl <- .hs_yields_lower(pr, row, val); hu <- .hs_yields_upper(pr, row, val)
    (!ureq && lreq && !hl) || (!lreq && ureq && !hu) || (ureq && lreq && !hl && !hu)
  }
  nz_row <- integer(); nz_j <- numeric(); nz_k <- numeric()
  store <- function(row, jv, kv) {
    if (pattern != 0L && (jv < 0) != (pattern * kv < 0)) return()
    nz_row <<- c(nz_row, row); nz_j <<- c(nz_j, jv); nz_k <<- c(nz_k, kv)
  }
  if (bcol != 0L) {
    if (pr$colsize[col] < pr$colsize[bcol]) {
      for (p in .hs_col_slots(pr, col)) {
        if (skip(pr$Arow[p], pr$Avalue[p])) next
        q <- .hs_find_nonzero(pr, pr$Arow[p], bcol)
        if (q == 0L) next
        store(pr$Arow[p], pr$Avalue[p], pr$Avalue[q])
      }
    } else {
      for (p in .hs_col_slots(pr, bcol)) {
        q <- .hs_find_nonzero(pr, pr$Arow[p], col)
        if (q == 0L) next
        if (skip(pr$Arow[p], pr$Avalue[q])) next
        store(pr$Arow[p], pr$Avalue[q], pr$Avalue[p])
      }
    }
  } else {
    for (p in .hs_col_slots(pr, col)) {
      if (skip(pr$Arow[p], pr$Avalue[p])) next
      nz_row <- c(nz_row, pr$Arow[p]); nz_j <- c(nz_j, pr$Avalue[p]); nz_k <- c(nz_k, -Inf)
    }
  }
  out <- list(lower = -Inf, upper = Inf, wclower = -Inf, wcupper = Inf)
  bound <- function(k, rhs, dir, wc) {
    jv <- nz_j[k]
    if ((dir > 0 && !wc) || (dir < 0 && wc)) {
      res <- .hs_lsb_residual2(pr, if (wc) "L" else "LO", nz_row[k], col, jv, bcol, nz_k[k], bval)
      if (res == -Inf) return(if (jv < 0) -Inf else Inf)
    } else {
      res <- .hs_lsb_residual2(pr, if (wc) "U" else "UO", nz_row[k], col, jv, bcol, nz_k[k], bval)
      if (res == Inf) return(if (jv < 0) Inf else -Inf)
    }
    .cd_val(.cd_div(.cd_add(c(rhs, 0), -res), jv))
  }
  upd <- function(k, rhs, dir) {
    if (dir * rhs == Inf) return()
    if (dir * nz_j[k] < 0) {
      if ("lower" %in% want) out$lower <<- max(out$lower, bound(k, rhs, dir, FALSE))
      if ("wclower" %in% want && out$wclower != Inf) out$wclower <<- max(out$wclower, bound(k, rhs, dir, TRUE))
    } else {
      if ("upper" %in% want) out$upper <<- min(out$upper, bound(k, rhs, dir, FALSE))
      if ("wcupper" %in% want && out$wcupper != -Inf) out$wcupper <<- min(out$wcupper, bound(k, rhs, dir, TRUE))
    }
  }
  for (k in seq_along(nz_row)) {
    upd(k, pr$ru[nz_row[k]], 1L)
    upd(k, pr$rl[nz_row[k]], -1L)
  }
  out
}

# ---------------------------------------------------------------------------
# Integrality helpers of HPresolve
# ---------------------------------------------------------------------------

.hs_clear_lifting <- function(pr, row) {
  key <- as.character(row)
  if (!is.null(pr$liftingOpportunities[[key]])) pr$liftingOpportunities[[key]] <- list()
}

.hs_row_coefs_integral <- function(pr, row, scale) {
  for (p in .hs_row_preorder(pr, row))
    if (.hs_fractionality(pr$Avalue[p] * scale) > pr$small_matrix_value) return(FALSE)
  TRUE
}

# HPresolve::isImpliedIntegral (an integer column); may tighten row sides
.hs_is_implied_integral <- function(pr, col) {
  dual <- TRUE
  slots <- .hs_col_slots(pr, col)
  for (p in slots) {
    row <- pr$Arow[p]
    if (pr$rowsize[row] < 2L || pr$rowsizeInteger[row] < pr$rowsize[row]) { dual <- FALSE; next }
    rlo <- if (.hs_impl_eq_at_upper(pr, row)) pr$ru[row] else pr$rl[row]
    rup <- if (.hs_impl_eq_at_lower(pr, row)) pr$rl[row] else pr$ru[row]
    if (rup == rlo) {
      dual <- FALSE
      scale <- 1 / pr$Avalue[p]
      if (!.hs_row_coefs_integral(pr, row, scale)) next
      if (.hs_fractionality(rlo * scale) > pr$ptol) .hs_pre_abort("infeasible")
      return(TRUE)
    }
  }
  if (!dual) return(FALSE)
  for (p in slots) {
    row <- pr$Arow[p]; v <- pr$Avalue[p]
    scale <- 1 / v
    if (!.hs_row_coefs_integral(pr, row, scale)) return(FALSE)
    if (pr$ru[row] != Inf) {
      ru <- abs(v) * floor(pr$ru[row] * abs(scale) + pr$ptol)
      if (abs(pr$ru[row] - ru) > pr$small_matrix_value) { pr$ru[row] <- ru; .hs_mark_row(pr, row) }
    }
    if (pr$rl[row] != -Inf) {
      rl <- abs(v) * ceiling(pr$rl[row] * abs(scale) - pr$ptol)
      if (abs(pr$rl[row] - rl) > pr$small_matrix_value) { pr$rl[row] <- rl; .hs_mark_row(pr, row) }
    }
  }
  TRUE
}

# HPresolve::isImpliedInteger (a continuous column)
.hs_is_implied_integer <- function(pr, col) {
  dual <- TRUE
  slots <- .hs_col_slots(pr, col)
  for (p in slots) {
    row <- pr$Arow[p]
    if (pr$rowsize[row] < 2L || pr$rowsizeInteger[row] + pr$rowsizeImplInt[row] < pr$rowsize[row] - 1L) {
      dual <- FALSE; next
    }
    rlo <- if (.hs_impl_eq_at_upper(pr, row)) pr$ru[row] else pr$rl[row]
    rup <- if (.hs_impl_eq_at_lower(pr, row)) pr$rl[row] else pr$ru[row]
    if (rup == rlo) {
      dual <- FALSE
      scale <- 1 / pr$Avalue[p]
      if (.hs_fractionality(rlo * scale) > pr$ptol) next
      if (!.hs_row_coefs_integral(pr, row, scale)) next
      return(TRUE)
    }
  }
  if (!dual) return(FALSE)
  smv <- pr$small_matrix_value
  if ((pr$cl[col] != -Inf && .hs_fractionality(pr$cl[col]) > smv) ||
      (pr$cu[col] != Inf && .hs_fractionality(pr$cu[col]) > smv)) return(FALSE)
  for (p in slots) {
    row <- pr$Arow[p]; scale <- 1 / pr$Avalue[p]
    if (pr$ru[row] != Inf && .hs_fractionality(pr$ru[row] * scale) > pr$ptol) return(FALSE)
    if (pr$rl[row] != -Inf && .hs_fractionality(pr$rl[row] * scale) > pr$ptol) return(FALSE)
    if (!.hs_row_coefs_integral(pr, row, scale)) return(FALSE)
  }
  TRUE
}

# HPresolve::convertImpliedInteger(col, row = -1, skipInputChecks = false)
.hs_convert_implied_integer <- function(pr, col, row = 0L, skip_checks = FALSE) {
  if (pr$colDeleted[col]) return(FALSE)
  if (!skip_checks) {
    if (pr$integrality[col] != 0L) return(FALSE)
    if (!.hs_is_implied_integer(pr, col)) return(FALSE)
  }
  pr$integrality[col] <- 4L
  if (row != 0L) {
    pr$rowsizeImplInt[row] <- pr$rowsizeImplInt[row] + 1L
  } else {
    for (p in .hs_col_slots(pr, col)) pr$rowsizeImplInt[pr$Arow[p]] <- pr$rowsizeImplInt[pr$Arow[p]] + 1L
  }
  .hs_change_col_lower(pr, col, pr$cl[col])
  .hs_change_col_upper(pr, col, pr$cu[col])
  TRUE
}

# HPresolve::transformColumn: x = scale * x' + constant
.hs_transform_column <- function(pr, col, scale, constant) {
  if (!is.null(pr$mipsolver)) .hmip_HighsImplications_columnTransformed(pr$mipsolver$mipdata_$implications,
                                                                       pr$ccur[col], scale, constant)
  .hs_ps_push(pr, list(type = "lineartransform", col = col, scale = scale, constant = constant))
  if (constant != 0) {
    ol <- pr$cl[col]; ou <- pr$cu[col]
    pr$cu[col] <- pr$cu[col] - constant
    pr$cl[col] <- pr$cl[col] - constant
    slots <- .hs_col_slots(pr, col)
    for (p in slots) {
      .hs_lsb_updated_var_lo(pr, "r", pr$Arow[p], col, pr$Avalue[p], ol)
      .hs_lsb_updated_var_up(pr, "r", pr$Arow[p], col, pr$Avalue[p], ou)
    }
    oil <- pr$implColLower[col]; oiu <- pr$implColUpper[col]
    pr$implColLower[col] <- pr$implColLower[col] - constant
    pr$implColUpper[col] <- pr$implColUpper[col] - constant
    for (p in slots) {
      .hs_lsb_updated_impl_lo(pr, "r", pr$Arow[p], col, pr$Avalue[p], pr$cl[col], oil, pr$colLowerSource[col])
      .hs_lsb_updated_impl_up(pr, "r", pr$Arow[p], col, pr$Avalue[p], pr$cu[col], oiu, pr$colUpperSource[col])
    }
  }
  .hs_lsb_sum_scaled(pr, "d", col, scale)
  bs <- 1 / scale
  pr$cl[col] <- pr$cl[col] * bs; pr$cu[col] <- pr$cu[col] * bs
  pr$implColLower[col] <- pr$implColLower[col] * bs; pr$implColUpper[col] <- pr$implColUpper[col] * bs
  if (scale < 0) {
    t <- pr$cl[col]; pr$cl[col] <- pr$cu[col]; pr$cu[col] <- t
    t <- pr$implColLower[col]; pr$implColLower[col] <- pr$implColUpper[col]; pr$implColUpper[col] <- t
    t <- pr$colLowerSource[col]; pr$colLowerSource[col] <- pr$colUpperSource[col]; pr$colUpperSource[col] <- t
  }
  pr$offset <- pr$offset + pr$cost[col] * constant
  pr$cost[col] <- pr$cost[col] * scale
  p <- pr$colhead[col]
  while (p != 0L) {
    v <- pr$Avalue[p]
    pr$Avalue[p] <- v * scale
    row <- pr$Arow[p]; rc <- v * constant
    if (pr$rl[row] != -Inf) pr$rl[row] <- pr$rl[row] - rc
    if (pr$ru[row] != Inf) pr$ru[row] <- pr$ru[row] - rc
    p <- pr$Anext[p]
  }
  if (pr$integrality[col] != 0L) {
    .hs_change_col_lower(pr, col, pr$cl[col])
    .hs_change_col_upper(pr, col, pr$cu[col])
  }
  .hs_mark_col(pr, col)
}

# HPresolve::scaleStoredRow on the positions rowpos
.hs_scale_stored_row <- function(pr, row, rowpos, scale, integral) {
  pr$ru[row] <- pr$ru[row] * scale
  pr$rl[row] <- pr$rl[row] * scale
  pr$implRowDualLower[row] <- pr$implRowDualLower[row] / scale
  pr$implRowDualUpper[row] <- pr$implRowDualUpper[row] / scale
  if (integral) {
    if (pr$ru[row] != Inf) pr$ru[row] <- .hs_round(pr$ru[row])
    if (pr$rl[row] != -Inf) pr$rl[row] <- .hs_round(pr$rl[row])
  }
  for (p in rowpos) {
    pr$Avalue[p] <- pr$Avalue[p] * scale
    if (abs(pr$Avalue[p]) <= pr$small_matrix_value) .hs_unlink(pr, p)
  }
  .hs_lsb_sum_scaled(pr, "r", row, scale)
  if (scale < 0) {
    t <- pr$rdl[row]; pr$rdl[row] <- pr$rdu[row]; pr$rdu[row] <- t
    t <- pr$implRowDualLower[row]; pr$implRowDualLower[row] <- pr$implRowDualUpper[row]; pr$implRowDualUpper[row] <- t
    t <- pr$rowDualLowerSource[row]; pr$rowDualLowerSource[row] <- pr$rowDualUpperSource[row]; pr$rowDualUpperSource[row] <- t
    t <- pr$rl[row]; pr$rl[row] <- pr$ru[row]; pr$ru[row] <- t
  }
}

# a doubleton equation record without a row (HighsPostsolveStack::
# doubletonEquation with row -1 and an empty column vector)
.hs_ps_subst_col <- function(pr, col, stay, scale, offset, lower, upper) {
  .hs_ps_push(pr, list(type = "doubleton", row = 0L, colsubst = col, col = stay, coefsubst = 1, coef = -scale,
                       rhs = offset, substlower = lower, substupper = upper, substcost = 0, lt = FALSE, ut = FALSE,
                       rowtype = "eq", colv = list(index = integer(), value = numeric())))
}
.hs_ps_removed_fixed_col <- function(pr, col, fix) {
  .hs_ps_push(pr, list(type = "fixedcol", col = col, fix = fix, cost = 0, status = "nonbasic",
                       colv = list(index = integer(), value = numeric())))
}

# ---------------------------------------------------------------------------
# rowPresolve: the integer parts
# ---------------------------------------------------------------------------

# the opposite bound test on binary columns (Suhl and Szymanski)
.hs_mip_degree1_tests <- function(pr, row) {
  for (p in .hs_row_preorder(pr, row)) {
    col <- pr$Acol[p]; val <- pr$Avalue[p]
    if (pr$integrality[col] == 0L || pr$cu[col] != pr$cl[col] + 1) next
    offset <- function() .cd_mul(.cd_sub_cd(c(pr$cu[col], 0), c(pr$cl[col], 0)), abs(val))
    test <- function(dir, act, side) {
      if (dir * act >= dir * side - pr$ptol) return()
      if (dir * val > 0) .hs_change_col_lower(pr, col, pr$cl[col] + 1)
      else .hs_change_col_upper(pr, col, pr$cu[col] - 1)
    }
    test(1L, .hs_lsb_sum_cd(pr, "r", "UO", row, -offset()), pr$rl[row])
    test(-1L, .hs_lsb_sum_cd(pr, "r", "LO", row, offset()), pr$ru[row])
  }
}

# the integer block of rowPresolve after the doubleton equations; TRUE when
# rowPresolve returns
.hs_mip_row_presolve <- function(pr, row) {
  smv <- pr$small_matrix_value; tol <- pr$ptol
  rlo <- pr$rl[row]; rup <- pr$ru[row]
  if (rlo == rup) {
    irl <- .hs_lsb_sum(pr, "r", "L", row); iru <- .hs_lsb_sum(pr, "r", "U", row)
    if (irl != -Inf && iru != Inf && abs(irl + iru - 2 * rup) <= smv) {
      bincol <- 0L; bincoef <- abs(iru - rup)
      rowpos <- .hs_row_inorder(pr, row)
      for (p in rowpos) {
        j <- pr$Acol[p]
        if (abs(abs(pr$Avalue[p]) - bincoef) <= smv && pr$integrality[j] == 1L &&
            abs(pr$cu[j] - pr$cl[j] - 1) <= tol) {
          bincol <- j; bincoef <- pr$Avalue[p]
          break
        }
      }
      if (bincol != 0L) {
        for (p in rowpos) {
          col <- pr$Acol[p]
          if (col == bincol) next
          lo <- .hs_lsb_impl_var_lower(pr, row, col); up <- .hs_lsb_impl_var_upper(pr, row, col)
          if (pr$cl[col] == pr$cu[col]) {
            .hs_ps_removed_fixed_col(pr, col, pr$cl[col])
            .hs_remove_fixed_col(pr, col)
            next
          }
          dir <- if ((bincoef < 0) == (pr$Avalue[p] < 0)) 1 else -1
          bound <- if (dir >= 0) up else lo
          scale <- dir * (lo - up)
          offset <- bound - pr$cl[bincol] * scale
          .hs_ps_subst_col(pr, col, bincol, scale, offset, lo, up)
          .hs_substitute_col(pr, col, bincol, offset, scale)
        }
        .hs_remove_row(pr, row)
        .hs_remove_row_singletons(pr)
        return(TRUE)
      }
    }
    if (pr$rowsizeInteger[row] + pr$rowsizeImplInt[row] >= pr$rowsize[row] - 1L) {
      ccol <- 0L; ccoef <- 0; coefs <- numeric()
      rowpos <- .hs_row_inorder(pr, row)
      for (p in rowpos) {
        if (pr$integrality[pr$Acol[p]] == 0L) { ccoef <- pr$Avalue[p]; ccol <- pr$Acol[p]; next }
        coefs <- c(coefs, pr$Avalue[p])
      }
      if (ccoef != 0) {
        is <- .hmip_integralScale(c(coefs, rup), smv, smv)
        if (is != 0 && is <= 1e3) {
          scale <- 1 / abs(ccoef * is)
          if (scale != 1) {
            .hs_transform_column(pr, ccol, scale, 0)
            .hs_convert_implied_integer(pr, ccol, 0L, TRUE)
            if (is != 1) .hs_scale_stored_row(pr, row, rowpos, is, TRUE)
          }
        }
      } else {
        is <- .hmip_integralScale(coefs, smv, smv)
        if (is != 0 && is <= 1e3) {
          rhs <- rup * is
          if (.hs_fractionality(rhs) > tol) .hs_pre_abort("infeasible")
          rhs <- .hs_round(rhs)
          x1c <- 0L; d <- 0
          for (i in seq_along(rowpos)) {
            a <- abs(.hs_round(is * pr$Avalue[rowpos[i]]))
            ng <- if (d == 0) a else .hmip_gcd(a, d)
            if (ng == 1) {
              if (x1c != 0L) { x1c <- 0L; break }
              x1c <- i
            } else d <- ng
          }
          if (x1c != 0L) {
            x1p <- rowpos[x1c]; x1 <- pr$Acol[x1p]
            if (.hs_fractionality(rhs / d) <= pr$epsilon) {
              .hs_transform_column(pr, x1, d, 0)
            } else {
              a1 <- .hmip_mod_int(.hs_round(is * pr$Avalue[x1p]), d)
              b <- .hmip_mod(.hmip_modularInverse(a1, d) * rhs, d)
              zl <- ceiling((pr$cl[x1] - b) / d - tol)
              zu <- floor((pr$cu[x1] - b) / d + tol)
              if (zl == zu) {
                fv <- zl * d + b
                if (fv > pr$cl[x1]) .hs_change_col_lower(pr, x1, fv)
                if (fv < pr$cu[x1]) .hs_change_col_upper(pr, x1, fv)
                .hs_fix_col_to(pr, x1, abs(pr$cl[x1] - fv) > tol)
                rowpos <- rowpos[-x1c]
              } else .hs_transform_column(pr, x1, d, b)
            }
            is <- is / d
          }
          if (is != 1) .hs_scale_stored_row(pr, row, rowpos, is, TRUE)
        }
      }
    }
    return(FALSE)
  }
  rowpos <- .hs_row_inorder(pr, row)
  if (pr$rowsize[row] == pr$rowsizeInteger[row] + pr$rowsizeImplInt[row])
    .hs_mip_integral_row(pr, row, rowpos)
  if (.hs_is_redundant(pr, row)) {
    .hs_ps_redundant_row(pr, row)
    .hs_remove_row(pr, row)
    return(TRUE)
  }
  strengthen <- function(rhs, dir, maxabs) {
    for (p in rowpos) {
      col <- pr$Acol[p]; val <- dir * pr$Avalue[p]
      lo <- .hs_lsb_impl_var_lower(pr, row, col); up <- .hs_lsb_impl_var_upper(pr, row, col)
      if (pr$integrality[col] == 0L) next
      if (val > .cd_val(.cd_add(maxabs, tol))) {
        delta <- .cd_mul(.cd_add(maxabs, -val), dir)
        .hs_add_to_matrix(pr, row, col, .cd_val(delta))
        rhs <- .cd_iadd(rhs, .cd_mul(delta, up))
      } else if (val < .cd_val(.cd_add(-maxabs, -tol))) {
        delta <- .cd_mul(.cd_add(maxabs, val), -dir)
        .hs_add_to_matrix(pr, row, col, .cd_val(delta))
        rhs <- .cd_iadd(rhs, .cd_mul(delta, lo))
      }
    }
    rhs
  }
  if (pr$rl[row] == -Inf) {
    iru <- .hs_lsb_sum(pr, "r", "U", row)
    if (iru != Inf) pr$ru[row] <- .cd_val(strengthen(c(pr$ru[row], 0), 1L, .cd_add(c(iru, 0), -pr$ru[row])))
  }
  if (pr$ru[row] == Inf) {
    irl <- .hs_lsb_sum(pr, "r", "L", row)
    if (irl != -Inf) pr$rl[row] <- .cd_val(strengthen(c(pr$rl[row], 0), -1L, .cd_add(c(pr$rl[row], 0), -irl)))
  }
  FALSE
}

# the integral scaling and Chvatal-Gomory strengthening of an inequality
# whose columns are all integral
.hs_mip_integral_row <- function(pr, row, rowpos) {
  smv <- pr$small_matrix_value; tol <- pr$ptol
  coefs <- pr$Avalue[rowpos]; idx <- pr$Acol[rowpos]
  is <- .hmip_integralScale(coefs, if (pr$rl[row] == -Inf) tol else smv, if (pr$ru[row] == Inf) tol else smv)
  if (is != 0) {
    lfin <- pr$rl[row] != -Inf; rfin <- pr$ru[row] != Inf
    lhs <- c(pr$rl[row], 0); rhs <- c(pr$ru[row], 0)
    if (lfin) lhs <- .cd_mul(lhs, is)
    if (rfin) rhs <- .cd_mul(rhs, is)
    rl_r <- c(-Inf, 0); rr_r <- c(Inf, 0); fl <- c(0, 0); fr <- c(0, 0)
    minl <- 0; minr <- 0; maxval <- 0; lt <- FALSE; rt <- FALSE
    ok <- TRUE
    for (i in seq_along(coefs)) {
      sc <- .cd_mul(c(coefs[i], 0), is)
      ic <- .cd_floor(.cd_add(sc, 0.5))
      cdelta <- .cd_sub_cd(ic, sc); dv <- .cd_val(cdelta)
      coefs[i] <- .cd_val(ic)
      maxval <- max(abs(coefs[i]), maxval)
      ub <- pr$cu[idx[i]]
      if (dv < -smv) {
        if (lfin) {
          if (ub == Inf) { ok <- FALSE; break }
          lhs <- .cd_iadd(lhs, .cd_mul(cdelta, ub))
        }
        minr <- max(-dv, minr)
      } else if (dv > smv) {
        if (rfin) {
          if (ub == Inf) { ok <- FALSE; break }
          rhs <- .cd_iadd(rhs, .cd_mul(cdelta, ub))
        }
        minl <- max(dv, minl)
      }
    }
    if (!ok) return(invisible())
    round_side <- function(v, mint, dir) {
      r <- .cd_mul(.cd_floor(.cd_add(.cd_mul(v, dir), tol)), dir)
      f <- .cd_mul(.cd_sub_cd(v, r), dir)
      list(r = r, f = f, t = .cd_val(f) >= mint - smv)
    }
    if (lfin) { o <- round_side(lhs, minl, -1); rl_r <- o$r; fl <- o$f; lt <- o$t }
    if (rfin) { o <- round_side(rhs, minr, 1); rr_r <- o$r; fr <- o$f; rt <- o$t }
    if (lfin && rfin && .cd_val(rr_r) < .cd_val(rl_r) - 0.5) .hs_pre_abort("infeasible")
    scale_row <- function(scalar, check) {
      if (.cd_val(rl_r) != -Inf) pr$rl[row] <- .cd_val(.cd_div(rl_r, scalar))
      if (.cd_val(rr_r) != Inf) pr$ru[row] <- .cd_val(.cd_div(rr_r, scalar))
      for (i in seq_along(coefs)) {
        delta <- .cd_val(.cd_add(.cd_div(c(coefs[i], 0), scalar), -pr$Avalue[rowpos[i]]))
        if (!check || abs(delta) > smv) .hs_add_to_matrix(pr, row, idx[i], delta)
      }
    }
    both <- lt && rt
    if (both || (lt && pr$ru[row] == Inf) || (rt && pr$rl[row] == -Inf)) {
      if (maxval > 1000 && is > 100) {
        if (both) {
          rl_r <- .cd_div(rl_r, is); rr_r <- .cd_div(rr_r, is)
          if (.cd_val(rr_r) < pr$ru[row] - tol) pr$ru[row] <- .cd_val(rr_r)
          if (.cd_val(rl_r) > pr$rl[row] + tol) pr$rl[row] <- .cd_val(rl_r)
        } else if ((rt && .cd_val(fr) < minr - tol) || (lt && .cd_val(fl) < minl - tol)) {
          scale_row(is, TRUE)
        }
      } else scale_row(1, FALSE)
    }
    return(invisible())
  }
  if (.hs_is_ranged(pr, row)) return(invisible())
  # Chvatal-Gomory strengthening (Achterberg et al., section 3.4)
  dir <- if (pr$ru[row] != Inf) -1 else 1
  rhs <- if (dir < 0) c(-pr$ru[row], 0) else c(pr$rl[row], 0)
  mina <- Inf; maxa <- 0
  for (i in seq_along(coefs)) {
    col <- idx[i]; val <- dir * coefs[i]; av <- abs(val)
    mina <- min(mina, av); maxa <- max(maxa, av)
    if (val < 0 && pr$cu[col] != Inf) rhs <- .cd_sub_cd(rhs, .cd_mul(c(pr$cu[col], 0), val))
    else if (val > 0 && pr$cl[col] != -Inf) rhs <- .cd_sub_cd(rhs, .cd_mul(c(pr$cl[col], 0), val))
    else return(invisible())
  }
  scalars <- 1
  for (t in 1:5) scalars <- c(scalars, t / maxa, t / mina, (2 * t - 1) / (2 * mina))
  scalars <- sort(unique(scalars))
  rc <- numeric(length(coefs)); rr <- c(0, 0)
  round_row <- function(s) {
    accept <- FALSE
    sc <- c(s, 0)
    rr <<- .cd_ceil(.cd_add(.cd_mul_cd(rhs, sc), -tol))
    ratio <- .cd_div_cd(rhs, rr)
    for (i in seq_along(coefs)) {
      ac <- abs(coefs[i])
      rc[i] <<- .cd_val(.cd_ceil(.cd_add(.cd_mul(sc, ac), -1e-14)))
      thr <- .cd_val(.cd_mul(ratio, rc[i]))
      if (ac < thr - smv) return(FALSE)
      accept <- accept || ac > thr + smv
    }
    accept
  }
  found <- FALSE
  for (s in scalars) if (round_row(s)) { found <- TRUE; break }
  if (!found) return(invisible())
  for (i in seq_along(coefs)) {
    col <- idx[i]; val <- dir * coefs[i]
    if (val < 0 && pr$cu[col] != Inf) {
      rc[i] <- -rc[i]
      rr <- .cd_iadd(rr, .cd_mul(c(pr$cu[col], 0), rc[i]))
    } else if (val > 0 && pr$cl[col] != -Inf) {
      rr <- .cd_iadd(rr, .cd_mul(c(pr$cl[col], 0), rc[i]))
    }
    rc[i] <- rc[i] * dir
  }
  rr <- .cd_mul(rr, dir)
  if (dir < 0) pr$ru[row] <- .cd_val(rr) else pr$rl[row] <- .cd_val(rr)
  for (i in seq_along(coefs)) {
    delta <- .cd_val(.cd_add(c(rc[i], 0), -coefs[i]))
    if (abs(delta) > smv) .hs_add_to_matrix(pr, row, idx[i], delta)
  }
}

# ---------------------------------------------------------------------------
# colPresolve and dualFixing: the mipsolver parts
# ---------------------------------------------------------------------------

.hs_mip_col_presolve <- function(pr, col) {
  modify <- function(row, dir, implied, ninf) {
    if (!implied || row == 0L || ninf != 1L || dir * pr$cost[col] < 0 || .hs_is_ranged(pr, row)) return()
    p <- .hs_find_nonzero(pr, row, col)
    if (pr$integrality[col] != 1L ||
        (pr$rowsizeInteger[row] == pr$rowsize[row] && .hs_row_coefs_integral(pr, row, 1 / pr$Avalue[p]))) {
      if (dir * pr$Avalue[p] > 0) .hs_change_impl_row_dual_lower(pr, row, 0, col)
      else .hs_change_impl_row_dual_upper(pr, row, 0, col)
    }
  }
  modify(pr$colLowerSource[col], 1L, .hs_is_lower_implied(pr, col), .hs_lsb_ninf(pr, "d", "UO", col))
  modify(pr$colUpperSource[col], -1L, .hs_is_upper_implied(pr, col), .hs_lsb_ninf(pr, "d", "LO", col))
  .hs_convert_implied_integer(pr, col)
  lo <- pr$cl[col]; up <- pr$cu[col]
  if (pr$integrality[col] != 0L && lo != 0 && (lo != -Inf || up != Inf) && up - lo > 0.5) {
    if (abs(up) > abs(lo)) {
      if (abs(lo) < 1000.5) .hs_transform_column(pr, col, 1, lo)
    } else if (abs(up) < 1000.5) .hs_transform_column(pr, col, -1, up)
  }
}

# substituteCol of dualFixing: a column locked by one row only takes the
# value of a binary column in that row
.hs_mip_dual_fix_subst <- function(pr, col, row, dir, colbound, otherbound) {
  lfin <- pr$rl[row] != -Inf; rfin <- pr$ru[row] != Inf
  for (p in .hs_row_inorder(pr, row)) {
    j <- pr$Acol[p]; v <- pr$Avalue[p]
    if (j == col) next
    if (pr$integrality[j] != 1L || pr$cl[j] != 0 || pr$cu[j] != 1) next
    if ((rfin && .hs_lsb_residual(pr, "r", "UO", row, j, v) > pr$ru[row] + pr$ptol) ||
        (lfin && .hs_lsb_residual(pr, "r", "LO", row, j, v) < pr$rl[row] - pr$ptol)) next
    best <- if (dir > 0) .hs_compute_col_bounds(pr, col, j, pr$cu[j], 0L, "lower")$lower
            else -.hs_compute_col_bounds(pr, col, j, pr$cu[j], 0L, "upper")$upper
    if (best >= dir * colbound - pr$ptol) {
      offset <- otherbound; scale <- colbound - otherbound
      .hs_ps_subst_col(pr, col, j, scale, offset, pr$cl[col], pr$cu[col])
      .hs_mark_col_deleted(pr, col)
      .hs_substitute_col(pr, col, j, offset, scale)
      break
    }
  }
}

# ---------------------------------------------------------------------------
# HPresolve::strengthenInequalities
# ---------------------------------------------------------------------------

.hs_strengthen_inequalities <- function(pr) {
  tol <- pr$ptol
  nstr <- 0L
  for (row in seq_len(pr$m)) {
    if (pr$rowDeleted[row] || pr$rowsize[row] <= 1L) next
    if (.hs_is_ranged(pr, row)) next
    limit <- max(1000L, (pr$ncol_now - pr$numDeletedCols) %/% 20L)
    if (pr$rowsize[row] > limit) next
    cc <- c(0, 0)
    if (pr$rl[row] != -Inf) { maxv <- c(pr$rl[row], 0); scale <- -1 }
    else { maxv <- c(-pr$ru[row], 0); scale <- 1 }
    comp <- integer(); rcost <- numeric(); upper <- numeric(); indices <- integer(); positions <- integer()
    skip <- FALSE
    stack <- pr$rowroot[row]
    while (length(stack)) {
      pos <- stack[length(stack)]; stack <- stack[-length(stack)]
      if (pr$ARright[pos] != 0L) stack <- c(stack, pr$ARright[pos])
      if (pr$ARleft[pos] != 0L) stack <- c(stack, pr$ARleft[pos])
      col <- pr$Acol[pos]
      skip <- pr$cl[col] == -Inf || pr$cu[col] == Inf
      if (skip) break
      w <- pr$Avalue[pos] * scale
      ub <- pr$cu[col] - pr$cl[col]
      if (w > 0) {
        cp <- 1L; maxv <- .cd_iadd(maxv, c(pr$cu[col] * w, 0))
      } else {
        cp <- -1L; maxv <- .cd_iadd(maxv, c(pr$cl[col] * w, 0)); w <- -w
      }
      if (ub <= tol || w <= tol) next
      if (pr$integrality[col] == 0L) { cc <- .cd_iadd(cc, c(w * ub, 0)); next }
      indices <- c(indices, length(rcost) + 1L); positions <- c(positions, pos)
      rcost <- c(rcost, w); comp <- c(comp, cp); upper <- c(upper, ub)
    }
    if (skip) next
    if (.cd_val(maxv) <= tol) {
      .hs_row_presolve(pr, row)
      next
    }
    small <- max(100 * tol, tol * .cd_val(maxv))
    repeat {
      if (.cd_val(.cd_sub_cd(maxv, cc)) <= small || !length(indices)) break
      indices <- indices[order(-rcost[indices], -indices)]
      lambda <- .cd_sub_cd(maxv, cc)
      cover <- integer()
      for (i in rev(seq_along(indices))) {
        ix <- indices[i]; delta <- upper[ix] * rcost[ix]
        if (upper[ix] <= 1000 && rcost[ix] > small && .cd_val(.cd_add(lambda, -delta)) <= small) cover <- c(cover, ix)
        else lambda <- .cd_iadd(lambda, c(-delta, 0))
      }
      if (!length(cover) || .cd_val(lambda) <= small) break
      alpos <- cover[1L]
      for (ix in cover[-1L]) {
        less <- if (rcost[ix] <= 1e-3 || rcost[alpos] <= 1e-3) rcost[ix] > rcost[alpos] else rcost[ix] < rcost[alpos]
        if (less) alpos <- ix
      }
      al <- rcost[alpos]
      coverrhs <- max(ceiling(.cd_val(.cd_add(.cd_div(lambda, al), -tol))), 1)
      slackupper <- c(-coverrhs, 0)
      step <- Inf
      cf <- numeric(length(cover))
      for (i in seq_along(cover)) {
        cf[i] <- ceiling(min(rcost[cover[i]], .cd_val(lambda)) / al - pr$small_matrix_value)
        slackupper <- .cd_iadd(slackupper, c(upper[cover[i]] * cf[i], 0))
        step <- min(step, rcost[cover[i]] / cf[i])
      }
      step <- min(step, .cd_val(.cd_div(maxv, coverrhs)))
      maxv <- .cd_iadd(maxv, c(-(step * coverrhs), 0))
      slackind <- length(rcost) + 1L
      rcost <- c(rcost, step); upper <- c(upper, .cd_val(slackupper))
      for (i in seq_along(cover)) rcost[cover[i]] <- rcost[cover[i]] - step * cf[i]
      indices <- indices[rcost[indices] > tol]
      indices <- c(indices, slackind)
    }
    thr <- .cd_val(.cd_add(maxv, tol))
    indices <- indices[indices <= length(positions) & abs(rcost[indices]) > thr]
    if (!length(indices)) next
    upd <- function(rhs, dir) {
      for (i in indices) {
        cd <- dir * .cd_val(.cd_add(-maxv, rcost[i]))
        col <- pr$Acol[positions[i]]
        if (comp[i] == -1L) {
          rhs <- .cd_iadd(rhs, c(cd * pr$cl[col], 0))
          .hs_add_to_matrix(pr, row, col, cd)
        } else {
          rhs <- .cd_iadd(rhs, c(-(cd * pr$cu[col]), 0))
          .hs_add_to_matrix(pr, row, col, -cd)
        }
      }
      rhs
    }
    if (scale < 0) pr$rl[row] <- .cd_val(upd(c(pr$rl[row], 0), -1))
    else pr$ru[row] <- .cd_val(upd(c(pr$ru[row], 0), 1))
    nstr <- nstr + length(indices)
  }
  if (nstr > 0L) .hs_log_dev(pr, sprintf("Strengthened %d coefficients", nstr))
  nstr
}

.hs_detect_implied_integers <- function(pr) {
  for (col in seq_len(pr$n)) if (!pr$colDeleted[col]) .hs_convert_implied_integer(pr, col)
}

.hs_log_dev <- function(pr, line) pr$log <- c(pr$log, line)

# ---------------------------------------------------------------------------
# HPresolve::dominatedColumns
# ---------------------------------------------------------------------------

# the top five bits of HighsHashHelpers::hash(row) for the 0-based row
.hs_row_sig_bit <- function(row0) .hs_pair_key_hash(row0, 0)[4L] %/% 2048 + 1

.hs_dominated_columns <- function(pr) {
  n <- pr$n; smv <- pr$small_matrix_value; tol <- pr$ptol
  ct <- pr$mipsolver$mipdata_$cliquetable
  sigm <- matrix(FALSE, 32L, n); sigp <- matrix(FALSE, 32L, n)
  bitcache <- integer(pr$m)
  is_binary <- function(i) pr$integrality[i] == 1L && pr$cl[i] == 0 && pr$cu[i] == 1
  for (i in seq_along(pr$Avalue)) {
    if (pr$Avalue[i] == 0) next
    row <- pr$Arow[i]; col <- pr$Acol[i]
    if (bitcache[row] == 0L) bitcache[row] <- as.integer(.hs_row_sig_bit(pr$rcur[row]))
    b <- bitcache[row]
    lf <- pr$rl[row] != -Inf; uf <- pr$ru[row] != Inf
    if (pr$Avalue[i] > 0) { sigm[b, col] <- sigm[b, col] | lf; sigp[b, col] <- sigp[b, col] | uf }
    else { sigm[b, col] <- sigm[b, col] | uf; sigp[b, col] <- sigp[b, col] | lf }
  }
  nfixed <- 0L; nbounds <- 0L
  dom_nz <- function(row, aj, ak) {
    if (.hs_is_ranged(pr, row)) return(abs(aj - ak) <= smv)
    if (pr$ru[row] == Inf) { aj <- -aj; ak <- -ak }
    !(aj > ak + smv)
  }
  check_dom <- function(sj, j, sk, k) {
    if (pr$integrality[j] == 1L && pr$integrality[k] == 0L) return(FALSE)
    jm <- sigm[, j]; jp <- sigp[, j]
    if (sj == -1L) { t <- jp; jp <- jm; jm <- t }
    km <- sigm[, k]; kp <- sigp[, k]
    if (sk == -1L) { t <- kp; kp <- km; km <- t }
    if (any(!jm & km)) return(FALSE)
    if (any(jp & !kp)) return(FALSE)
    if (sj * pr$cost[j] > sk * pr$cost[k] + smv) return(FALSE)
    for (p in .hs_col_slots(pr, j)) {
      row <- pr$Arow[p]; q <- .hs_find_nonzero(pr, row, k)
      if (!dom_nz(row, sj * pr$Avalue[p], sk * (if (q == 0L) 0 else pr$Avalue[q]))) return(FALSE)
    }
    for (p in .hs_col_slots(pr, k)) {
      row <- pr$Arow[p]
      if (.hs_find_nonzero(pr, row, j) != 0L) next
      if (!dom_nz(row, 0, sk * pr$Avalue[p])) return(FALSE)
    }
    TRUE
  }
  fix_col <- function(col, dir) {
    nfixed <<- nfixed + 1L
    .hs_fix_col_to(pr, col, dir > 0)
    .hs_remove_row_singletons(pr)
    .hs_remove_doubleton_equations(pr)
  }
  tighten <- function(col, cbound, at_upper, other, obound, pattern) {
    if (pr$cl[col] == pr$cu[col]) return()
    isint <- pr$integrality[col] != 0L
    if (at_upper) {
      b <- .hs_compute_col_bounds(pr, col, other, obound, pattern, c("lower", "upper", "wcupper"))
      ub <- b$upper; lb <- min(cbound, b$lower)
      if (pr$cost[col] <= 0) {
        wc <- b$wcupper
        if (isint) wc <- floor(wc + tol)
        lb <- max(lb, min(cbound, wc))
      }
    } else {
      b <- .hs_compute_col_bounds(pr, col, other, obound, pattern, c("lower", "upper", "wclower"))
      lb <- b$lower; ub <- max(cbound, b$upper)
      if (pr$cost[col] >= 0) {
        wc <- b$wclower
        if (isint) wc <- ceiling(wc - tol)
        ub <- min(ub, max(cbound, wc))
      }
    }
    if (lb > pr$cl[col] + tol) {
      if (isint) lb <- ceiling(lb - tol)
      if (lb == pr$cu[col]) fix_col(col, 1L)
      else if (isint) { nbounds <<- nbounds + 1L; .hs_change_col_lower(pr, col, lb) }
    }
    if (pr$colDeleted[col]) return()
    if (ub < pr$cu[col] - tol) {
      if (isint) ub <- floor(ub + tol)
      if (ub == pr$cl[col]) fix_col(col, -1L)
      else if (isint) { nbounds <<- nbounds + 1L; .hs_change_col_upper(pr, col, ub) }
    }
  }
  same_type <- function(a, b) {
    (pr$integrality[a] != 0L && pr$integrality[b] != 0L) || (pr$integrality[a] != 1L && pr$integrality[b] != 1L)
  }
  check_cols <- function(col, k, dir, mult, implied, cliques, st) {
    dk <- mult * dir
    domb <- if (dir > 0) pr$cu[col] else pr$cl[col]
    dedb <- if (dk > 0) pr$cl[k] else pr$cu[k]
    domfin <- dir * domb != Inf; dedfin <- dk * dedb != -Inf
    tryfix <- dedfin && (implied || cliques)
    trystr <- (domfin || dedfin) && st
    if (!(tryfix || trystr)) return()
    if (!check_dom(dir, col, dk, k)) return()
    if (tryfix && (implied || .hmip_HighsCliqueTable_haveCommonClique(ct, 2L * pr$ccur[col] + (dir > 0),
                                                                     2L * pr$ccur[k] + (dk > 0)))) {
      fix_col(k, -dk)
    } else if (trystr) {
      if (dedfin) tighten(col, domb, dir > 0, k, dedb, dir * dk)
      if (!pr$colDeleted[col] && domfin) tighten(k, dedb, dk < 0, col, domb, dir * dk)
    }
  }
  check_row <- function(row, col, dir, best, implied, cliques) {
    for (p in .hs_row_inorder(pr, row)) {
      k <- pr$Acol[p]
      if (k == col || pr$colDeleted[k]) next
      ak <- pr$Avalue[p]
      st <- same_type(col, k)
      if (!implied && !cliques && !st) next
      if (dom_nz(row, dir * best, dir * ak)) {
        check_cols(col, k, dir, 1L, implied, cliques, st)
        if (pr$colDeleted[col]) break
      }
      if (!pr$colDeleted[k] && dom_nz(row, dir * best, -dir * ak)) {
        check_cols(col, k, dir, -1L, implied, cliques, st)
        if (pr$colDeleted[col]) break
      }
    }
  }
  ncl <- function(col, val) {
    if (is.null(ct)) return(0L)
    .hmip_HighsCliqueTable_numCliques2(ct, pr$ccur[col], val)
  }
  for (j in seq_len(n)) {
    if (pr$colDeleted[j]) next
    bp <- 0L; bpl <- .Machine$integer.max; ajp <- 0
    bm <- 0L; bml <- .Machine$integer.max; ajm <- 0
    for (p in .hs_col_slots(pr, j)) {
      row <- pr$Arow[p]
      val <- (if (pr$ru[row] != Inf) 1 else -1) * pr$Avalue[p]
      if (val > 0 && pr$rowsize[row] < bpl) { bp <- row; bpl <- pr$rowsize[row]; ajp <- pr$Avalue[p] }
      if (val < 0 && pr$rowsize[row] < bml) { bm <- row; bml <- pr$rowsize[row]; ajm <- pr$Avalue[p] }
    }
    if (bm != 0L) check_row(bm, j, -1L, ajm, .hs_is_lower_implied(pr, j), is_binary(j) && ncl(j, 0L) > 0L)
    if (!pr$colDeleted[j] && bp != 0L)
      check_row(bp, j, 1L, ajp, .hs_is_upper_implied(pr, j), is_binary(j) && ncl(j, 1L) > 0L)
  }
  if (nfixed > 0L || nbounds > 0L)
    .hs_log_dev(pr, sprintf("Fixed %d dominated columns and strengthened %d bounds", nfixed, nbounds))
}

# HPresolve members used by the MIP rules; for an LP every column is
# continuous and pr$mip is FALSE, so the hooks in the LP port stay inert
.hs_mip_fields <- function(pr, integrality, mipsolver = NULL) {
  pr$integrality <- as.integer(integrality)
  pr$rowsizeInteger <- integer(pr$m); pr$rowsizeImplInt <- integer(pr$m)
  pr$mip <- !is.null(mipsolver)
  pr$mipsolver <- mipsolver
  pr$small_matrix_value <- 1e-9
  pr$mip_feastol <- 1e-6
  pr$epsilon <- 1e-9
  pr$liftingOpportunities <- list()
  pr$log <- character()
  invisible(pr)
}

# ---------------------------------------------------------------------------
# The model view. HiGHS renumbers rows and columns in shrinkProblem; the R
# presolve keeps its input indices and records the current numbering in
# pr$cidx / pr$ridx (input index of current column j0 + 1) and pr$ccur /
# pr$rcur (current 0-based index of an input column). The clique table,
# the implications and the domain see the current 0-based indices.
# ---------------------------------------------------------------------------

# .hs_num_nonzeros, .hs_sparsify and .hs_equality_row_addition are B1's (optimization-scipy-highs-presolve.R)

# HPresolve::toCSC over the current numbering (called after a shrink, so every
# live slot has a current row and column)
.hs_mip_to_csc <- function(pr) {
  live <- which(pr$Avalue != 0)
  ord <- live[order(pr$ccur[pr$Acol[live]], live)]
  cnt <- tabulate(pr$ccur[pr$Acol[ord]] + 1L, length(pr$cidx))
  list(format_ = "colwise", start_ = c(0L, cumsum(cnt)), index_ = as.integer(pr$rcur[pr$Arow[ord]]),
       value_ = pr$Avalue[ord])
}

# write the presolve state into the HighsLp the mipsolver points at
# (mipsolver.model_ is the presolve model in HiGHS)
.hs_mip_sync_model <- function(pr, csc = TRUE) {
  m <- pr$mipsolver$model_
  cols <- pr$cidx; rows <- pr$ridx
  m$num_col_ <- length(cols); m$num_row_ <- length(rows)
  m$col_cost_ <- pr$cost[cols]; m$col_lower_ <- pr$cl[cols]; m$col_upper_ <- pr$cu[cols]
  m$integrality_ <- pr$integrality[cols]
  m$row_lower_ <- pr$rl[rows]; m$row_upper_ <- pr$ru[rows]
  m$offset_ <- pr$offset; m$sense_ <- 1L
  if (csc) m$a_matrix_ <- .hs_mip_to_csc(pr)
  invisible(m)
}

# shrinkProblem (with the mipsolver rebuild) followed by toCSC and okFromCSC
.hs_mip_shrink <- function(pr) {
  oldc <- pr$cidx; oldr <- pr$ridx
  .hs_shrink(pr)
  pr$cidx <- which(!pr$colDeleted); pr$ridx <- which(!pr$rowDeleted)
  newc <- ifelse(pr$colDeleted[oldc], -1L, pr$ccur[oldc])
  newr <- ifelse(pr$rowDeleted[oldr], -1L, pr$rcur[oldr])
  ps <- pr$postsolve_stack
  if (!is.null(ps)) .hmip_HighsPostsolveStack_compressIndexMaps(ps, newr, newc)
  ms <- pr$mipsolver
  if (!is.null(ms)) {
    md <- ms$mipdata_
    .hs_mip_sync_model(pr)
    md$rowMatrixSet <- FALSE
    md$objectiveFunction <- .hmip_new_HighsObjectiveFunction(ms)
    md$domain <- .hmip_new_HighsDomain(ms)
    .hmip_HighsCliqueTable_rebuild(md$cliquetable, length(pr$cidx), ps, md$domain, newc, newr)
    .hmip_HighsImplications_rebuild(md$implications, length(pr$cidx), newc, newr)
    opt <- ms$options_mip_
    md$cutpool <- .hmip_new_HighsCutPool(length(pr$cidx), opt$mip_pool_age_limit, opt$mip_pool_soft_limit)
    md$conflictPool <- .hmip_new_HighsConflictPool(5 * opt$mip_pool_age_limit, opt$mip_pool_soft_limit)
  }
  invisible()
}

# toCSC and okFromCSC without a shrink (relinks the matrix in column order)
.hs_mip_relink <- function(pr) {
  .hs_shrink(pr)
  pr$cidx <- which(!pr$colDeleted); pr$ridx <- which(!pr$rowDeleted)
  invisible()
}

# ---------------------------------------------------------------------------
# HPresolve::applyConflictGraphSubstitutions
# ---------------------------------------------------------------------------

.hs_apply_conflict_graph_substitutions <- function(pr) {
  md <- pr$mipsolver$mipdata_
  impl <- md$implications; ct <- md$cliquetable
  # group A keeps both substitution lists as structures of arrays
  for (k in seq_along(impl$subst_substcol)) {
    sc <- pr$cidx[impl$subst_substcol[k] + 1L]; st <- pr$cidx[impl$subst_staycol[k] + 1L]
    if (pr$colDeleted[sc] || pr$colDeleted[st]) next
    scale <- impl$subst_scale[k]; offset <- impl$subst_offset[k]
    pr$probingNumDelCol <- pr$probingNumDelCol + 1L
    .hs_ps_subst_col(pr, sc, st, scale, offset, pr$cl[sc], pr$cu[sc])
    .hs_mark_col_deleted(pr, sc)
    .hs_substitute_col(pr, sc, st, offset, scale)
  }
  impl$subst_substcol <- integer(); impl$subst_staycol <- integer()
  impl$subst_scale <- numeric(); impl$subst_offset <- numeric()
  subs <- .hmip_HighsCliqueTable_getSubstitutions(ct)
  for (k in seq_along(subs$substcol)) {
    replace <- subs$replace[k]
    sc <- pr$cidx[subs$substcol[k] + 1L]; rc <- pr$cidx[.hmip_cv_col(replace) + 1L]
    if (pr$colDeleted[sc] || pr$colDeleted[rc]) next
    if (.hmip_cv_val(replace) == 0L) { scale <- -1; offset <- 1 } else { scale <- 1; offset <- 0 }
    pr$probingNumDelCol <- pr$probingNumDelCol + 1L
    .hs_ps_subst_col(pr, sc, rc, scale, offset, pr$cl[sc], pr$cu[sc])
    .hs_mark_col_deleted(pr, sc)
    .hs_substitute_col(pr, sc, rc, offset, scale)
  }
  ct$subst_substcol <- integer(); ct$subst_replace <- integer()
  invisible()
}

# ---------------------------------------------------------------------------
# HPresolve::runProbing
# ---------------------------------------------------------------------------

# HighsRandom::integer(): 31 bits of pair_hash<0>(state, state >> 32)
.hs_rng_int31 <- function(r) {
  .hs_rng_advance(r)
  v <- .hs_u64_shift(.hs_pair_hash(0L, .hs_lo32(r$state), .hs_hi32(r$state)), 33L)
  v[1L] + 65536 * v[2L]
}

.hs_run_probing <- function(pr) {
  ms <- pr$mipsolver; md <- ms$mipdata_; opt <- ms$options_mip_
  infeasible <- function() .hs_pre_abort("infeasible")
  pr$probingEarlyAbort <- FALSE
  if (pr$numDeletedCols + pr$numDeletedRows != 0L) .hs_mip_shrink(pr) else .hs_mip_relink(pr)
  ct <- md$cliquetable; impl <- md$implications
  if (!is.null(ms)) {
    .hs_mip_sync_model(pr)
    .hmip_HighsCliqueTable_setMaxEntries(ct, .hs_num_nonzeros(pr))
  }
  # kHighsTiny = 1e-14
  huge <- pr$ptol / 1e-14
  for (col in pr$cidx) {
    if (abs(pr$implColLower[col]) <= huge && pr$implColLower[col] > pr$cl[col])
      .hs_change_col_lower(pr, col, pr$implColLower[col])
    if (abs(pr$implColUpper[col]) <= huge && pr$implColUpper[col] < pr$cu[col])
      .hs_change_col_upper(pr, col, pr$implColUpper[col])
  }
  # without a mipsolver (the standalone .hs_mip_presolve) only the bound
  # tightening above runs; the rest needs the domain and the clique table
  if (is.null(ms)) return(invisible())
  oldNumProbed <- pr$numProbed
  .hs_mip_sync_model(pr, csc = FALSE)
  .hmip_HighsMipSolverData_setupDomainPropagation(md)
  domain <- md$domain
  .hmip_HighsDomain_propagate(domain)
  if (.hmip_HighsDomain_infeasible(domain)) infeasible()
  firstCall <- !isTRUE(md$cliquesExtracted)
  md$cliquesExtracted <- TRUE
  if (firstCall) {
    .hmip_HighsCliqueTable_extractCliques(ct, ms)
    if (.hmip_HighsDomain_infeasible(domain)) infeasible()
    if (md$upper_limit != Inf) {
      tmp <- md$upper_limit
      md$upper_limit <- tmp - pr$offset
      .hmip_HighsCliqueTable_extractObjCliques(ct, ms)
      md$upper_limit <- tmp
      if (.hmip_HighsDomain_infeasible(domain)) infeasible()
    }
    .hmip_HighsDomain_propagate(domain)
    if (.hmip_HighsDomain_infeasible(domain)) infeasible()
  }
  .hmip_HighsCliqueTable_cleanupFixed(ct, domain)
  if (.hmip_HighsDomain_infeasible(domain)) infeasible()
  ncol <- length(pr$cidx)
  # the tuples (-min(5000, up * down) / (1 + numProbes), -min(100, up + down),
  # random, col); pdqsort over distinct tuples gives the lexicographic order
  bk1 <- bk2 <- bk3 <- numeric(); bcol <- integer()
  if (!.hmip_HighsCliqueTable_isFull(ct)) {
    rng <- .hs_rng(opt$random_seed %||% 0)
    for (i in seq_len(ncol) - 1L) {
      if (!.hmip_HighsDomain_isBinary(domain, i)) next
      up <- .hmip_HighsCliqueTable_getNumImplications2(ct, i, 1L)
      down <- .hmip_HighsCliqueTable_getNumImplications2(ct, i, 0L)
      bk1 <- c(bk1, -(min(5000, up * down) %/% (1 + pr$numProbes[pr$cidx[i + 1L]])))
      bk2 <- c(bk2, -min(100, up + down)); bk3 <- c(bk3, .hs_rng_int31(rng)); bcol <- c(bcol, i)
    }
  }
  if (!length(bcol)) return(invisible())
  binaries <- bcol[order(bk1, bk2, bk3, bcol)]
  numChangedCols <- 0L
  count_fixed <- function() {
    ch <- .hmip_HighsDomain_getChangedCols(domain)
    while (length(ch) != numChangedCols) {
      numChangedCols <<- numChangedCols + 1L
      if (.hmip_HighsDomain_isFixed(domain, ch[numChangedCols])) pr$probingNumDelCol <- pr$probingNumDelCol + 1L
    }
  }
  count_fixed()
  numCliquesStart <- .hmip_HighsCliqueTable_numCliques(ct)
  numImplicsStart <- .hmip_HighsImplications_getNumImplications(impl)
  numDelStart <- pr$probingNumDelCol
  calcNumDel <- function() pr$probingNumDelCol - numDelStart + length(impl$subst_substcol) +
    length(.hmip_HighsCliqueTable_getSubstitutions(ct)$substcol)
  numDel <- calcNumDel()
  nnz <- .hs_num_nonzeros(pr)
  splayContingent <- ct$numNeighbourhoodQueries + max(if (isTRUE(ms$submip)) 0 else 100000, 10 * nnz)
  numFail <- 0
  maxNumLiftOpps <- max(100000, 10 * length(pr$ridx))
  pr$numLiftOpps <- 0
  lifting <- (opt$mip_lifting_for_probing %||% -1L) != -1L
  if (lifting) {
    ncont <- sum(!pr$colDeleted & pr$integrality == 0L); nall <- sum(!pr$colDeleted)
    lifting <- 100 * ncont >= 2 * nall
  }
  if (lifting) {
    impl$storeLiftingOpportunity <- function(row, col, val, coef) {
      key <- as.character(pr$ridx[row + 1L])
      lo <- pr$liftingOpportunities[[key]]
      if (is.null(lo)) lo <- list(col = integer(), val = integer(), coef = numeric())
      lo$col <- c(lo$col, col); lo$val <- c(lo$val, val); lo$coef <- c(lo$coef, coef)
      pr$liftingOpportunities[[key]] <- lo
      pr$numLiftOpps <- pr$numLiftOpps + 1
    }
  }
  for (i in binaries) {
    if (!is.null(.hmip_HighsCliqueTable_getSubstitution(ct, i)) || !.hmip_HighsDomain_isBinary(domain, i)) next
    tightenLimits <- (pr$numProbed - oldNumProbed) >= 2500
    nrc <- length(pr$ridx) + ncol
    pr$probingEarlyAbort <- if (!tightenLimits) numDel > max(1000, nrc %/% 20) else numDel > min(1000, nrc %/% 20)
    if (pr$probingEarlyAbort) break
    lim <- max(1000000, 2 * .hs_num_nonzeros(pr))
    if (.hmip_HighsCliqueTable_isFull(ct) || .hmip_HighsCliqueTable_numCliques(ct) - numCliquesStart > lim ||
        .hmip_HighsImplications_getNumImplications(impl) - numImplicsStart > lim) break
    if (ct$numNeighbourhoodQueries > splayContingent) break
    if (pr$probingContingent - pr$numProbed < 0) break
    numNewCliques <- -.hmip_HighsCliqueTable_numCliques(ct)
    res <- .hmip_HighsImplications_runProbing(impl, i, 0L)
    if (!isTRUE(res$result)) next
    pr$probingContingent <- pr$probingContingent + res$numReductions
    numNewCliques <- max(numNewCliques + .hmip_HighsCliqueTable_numCliques(ct), 0)
    count_fixed()
    newNumDel <- calcNumDel()
    if (newNumDel > numDel) {
      pr$probingContingent <- pr$probingContingent + numDel
      if (!isTRUE(ms$submip)) {
        splayContingent <- splayContingent + 100 * (newNumDel + numDelStart)
        splayContingent <- splayContingent + 1000 * numNewCliques
      }
      numDel <- newNumDel
      numFail <- 0
    } else if (isTRUE(ms$submip) || numNewCliques == 0) {
      splayContingent <- splayContingent - (if (tightenLimits) 250 else 100) * numFail
      numFail <- numFail + 1
    } else {
      splayContingent <- splayContingent + 1000 * numNewCliques
      numFail <- 0
    }
    pr$numProbed <- pr$numProbed + 1L
    oc <- pr$cidx[i + 1L]
    pr$numProbes[oc] <- pr$numProbes[oc] + 1
    if (pr$numLiftOpps >= maxNumLiftOpps) impl$storeLiftingOpportunity <- NULL
    if (.hmip_HighsDomain_infeasible(domain)) infeasible()
  }
  .hmip_HighsCliqueTable_cleanupFixed(ct, domain)
  if (!firstCall) .hmip_HighsCliqueTable_extractCliques(ct, ms, FALSE)
  .hmip_HighsCliqueTable_runCliqueMerging(ct, domain)
  for (delrow in .hmip_HighsCliqueTable_getDeletedRows(ct)) {
    row <- pr$ridx[delrow + 1L]
    if (!pr$rowDeleted[row]) .hs_remove_row(pr, row)
  }
  ct$deletedrows <- integer()
  ext <- .hmip_HighsCliqueTable_getCliqueExtensions(ct)
  addednnz <- length(ext$first)
  for (k in seq_along(ext$first)) {
    row <- pr$ridx[ext$first[k] + 1L]; v <- ext$second[k]
    if (pr$rowDeleted[row]) { addednnz <- addednnz - 1L; next }
    if (.hmip_cv_val(v) == 0L) {
      pr$rl[row] <- pr$rl[row] - 1; pr$ru[row] <- pr$ru[row] - 1
      val <- -1
    } else val <- 1
    .hs_add_to_matrix(pr, row, pr$cidx[.hmip_cv_col(v) + 1L], val)
  }
  ct$ext_row <- integer(); ct$ext_var <- integer()
  for (i in seq_len(ncol)) {
    col <- pr$cidx[i]
    if (pr$colDeleted[col]) next
    if (pr$cl[col] < domain$col_lower_[i]) .hs_change_col_lower(pr, col, domain$col_lower_[i])
    if (pr$cu[col] > domain$col_upper_[i]) .hs_change_col_upper(pr, col, domain$col_upper_[i])
    if (.hmip_HighsDomain_isFixed(domain, i - 1L)) {
      .hs_ps_removed_fixed_col(pr, col, pr$cl[col])
      .hs_remove_fixed_col(pr, col)
    }
  }
  .hs_apply_conflict_graph_substitutions(pr)
  .hs_log_dev(pr, sprintf("%d probing evaluations: %d deleted rows, %d deleted columns, %d lifted nonzeros",
                          pr$numProbed - oldNumProbed, pr$numDeletedRows, pr$numDeletedCols, addednnz))
  if ((opt$mip_lifting_for_probing %||% -1L) != -1L) {
    if (pr$numDeletedRows == 0L && pr$numDeletedCols == 0L && addednnz == 0L) .hs_lifting_for_probing(pr)
    pr$liftingOpportunities <- list()
    impl$storeLiftingOpportunity <- NULL
  }
  invisible()
}

# HPresolve::liftingForProbing (Achterberg et al. 2019). liftingOpportunities
# is keyed by the input row and holds current 0-based columns in insertion
# order; insert_or_get keeps the first coefficient of a (col, val) pair. The
# rows are visited in insertion order where C++ walks its hash table: only the
# order of the redundantRow records can differ, the lifting itself is sorted.
.hs_lifting_for_probing <- function(pr) {
  ms <- pr$mipsolver
  ct <- ms$mipdata_$cliquetable; domain <- ms$mipdata_$domain
  fillallowed <- ms$options_mip_$mip_lifting_for_probing > 0
  ncol_now <- length(pr$cidx)
  tab_row <- integer(); tab_clique <- list(); tab_score <- numeric(); tab_nfill <- integer()
  bestscoretotal <- -Inf
  numrowsremoved <- 0L
  for (key in names(pr$liftingOpportunities)) {
    row <- as.integer(key)
    if (pr$rowDeleted[row]) next
    dense <- pr$rowsize[row] > max(1000, (ncol_now - pr$numDeletedCols) %/% 20)
    lo <- pr$liftingOpportunities[[key]]
    first <- !duplicated(paste(lo$col, lo$val))
    tcol <- lo$col[first]; tval <- lo$val[first]; tcoef <- lo$coef[first]
    isredundant <- any(paste(tcol, 1L - tval) %in% paste(tcol, tval))
    cv <- integer(); cf <- numeric(); cp <- integer()
    for (k in seq_along(tcol)) {
      col <- pr$cidx[tcol[k] + 1L]
      pos <- .hs_find_nonzero(pr, row, col)
      if (!dense && (fillallowed || pos != 0L) && !pr$colDeleted[col] &&
          !.hmip_HighsDomain_isFixed(domain, tcol[k])) {
        cv <- c(cv, 2L * tcol[k] + tval[k]); cf <- c(cf, tcoef[k]); cp <- c(cp, pos)
      }
    }
    if (isredundant) {
      numrowsremoved <- numrowsremoved + 1L
      .hs_ps_redundant_row(pr, row)
      .hs_remove_row(pr, row)
      next
    }
    if (!length(cv)) next
    # the std::map orders by (col, val), which is the CliqueVar code order
    o <- order(cv); cv <- cv[o]; cf <- cf[o]; cp <- cp[o]
    diffs <- abs(cf - ifelse(cp == 0L, 0, pr$Avalue[pmax(cp, 1L)]))
    bestscore <- -Inf; bestclique <- integer(); bestnfill <- 0L
    for (k in seq_along(cv)) if (diffs[k] > bestscore) {
      bestscore <- diffs[k]; bestnfill <- as.integer(cp[k] == 0L); bestclique <- k
    }
    if (length(cv) > 1L) {
      for (clique in .hmip_HighsCliqueTable_computeMaximalCliques(ct, cv, pr$ptol)) {
        idx <- match(clique, cv)
        score <- c(0, 0)
        for (k in idx) score <- .cd_add(score, diffs[k])
        if (.cd_val(score) > bestscore) {
          bestscore <- .cd_val(score); bestnfill <- sum(cp[idx] == 0L); bestclique <- idx
        }
      }
    }
    tab_row <- c(tab_row, row); tab_clique[[length(tab_row)]] <- list(cv = cv[bestclique], coef = cf[bestclique])
    tab_score <- c(tab_score, bestscore); tab_nfill <- c(tab_nfill, bestnfill)
    bestscoretotal <- max(bestscoretotal, bestscore)
  }
  nelm <- vapply(tab_clique, function(x) length(x$cv), 0L)
  overall <- 0.5 * (tab_score / bestscoretotal) + (1 - 0.5) * (nelm - tab_nfill) / nelm
  nfill <- 0L; nmod <- 0L; numrowsmodified <- 0L
  maxnfill <- max(10L * length(tab_row), .hs_num_nonzeros(pr) %/% 100L)
  for (t in order(-overall, pr$rcur[tab_row])) {
    row <- tab_row[t]; newfill <- tab_nfill[t]
    if (nfill + newfill > maxnfill) break
    nfill <- nfill + newfill; nmod <- nmod + nelm[t] - newfill
    update <- c(0, 0)
    cl <- tab_clique[[t]]
    for (k in seq_along(cl$cv)) {
      .hs_add_to_matrix(pr, row, pr$cidx[cl$cv[k] %/% 2L + 1L], cl$coef[k])
      if (cl$cv[k] %% 2L == 0L) update <- .cd_add(update, cl$coef[k])
    }
    numrowsmodified <- numrowsmodified + 1L
    if (pr$rl[row] != -Inf) pr$rl[row] <- pr$rl[row] + .cd_val(update)
    if (pr$ru[row] != Inf) pr$ru[row] <- pr$ru[row] + .cd_val(update)
  }
  .hs_log_dev(pr, sprintf(paste("Lifting for probing removed %d and modified %d row(s), added %d",
                                "new and modified %d existing nonzero(s)"),
                          numrowsremoved, numrowsmodified, nfill, nmod))
  invisible()
}

# ---------------------------------------------------------------------------
# HPresolve::scaleMIP
# ---------------------------------------------------------------------------

.hs_scale_mip <- function(pr) {
  for (row in pr$ridx) {
    if (pr$rowDeleted[row] || pr$rowsize[row] < 1L ||
        pr$rowsizeInteger[row] + pr$rowsizeImplInt[row] == pr$rowsize[row]) next
    rowpos <- .hs_row_inorder(pr, row)
    mx <- 0
    for (p in rowpos) if (pr$integrality[pr$Acol[p]] == 0L) mx <- max(abs(pr$Avalue[p]), mx)
    scale <- 2^.hmip_round(-log2(mx))
    if (scale == 1) next
    if (pr$ru[row] == Inf) scale <- -scale
    .hs_scale_stored_row(pr, row, rowpos, scale, FALSE)
  }
  for (col in pr$cidx) {
    if (pr$colDeleted[col] || pr$colsize[col] < 1L || pr$integrality[col] != 0L) next
    mx <- 0
    for (p in .hs_col_slots(pr, col)) mx <- max(abs(pr$Avalue[p]), mx)
    scale <- 2^.hmip_round(-log2(mx))
    if (scale == 1) next
    .hs_transform_column(pr, col, scale, 0)
  }
}

# ---------------------------------------------------------------------------
# HPresolve::presolve for a MIP
# ---------------------------------------------------------------------------

.hs_mip_report <- function(pr) {
  .hs_log_dev(pr, sprintf("%d rows, %d cols, %d nonzeros", length(pr$ridx) - pr$numDeletedRows,
                          length(pr$cidx) - pr$numDeletedCols, .hs_num_nonzeros(pr)))
}

.hs_mip_presolve_main <- function(pr) {
  ms <- pr$mipsolver; have_ms <- !is.null(ms)
  ct <- if (have_ms) ms$mipdata_$cliquetable
  if (pr$presolve_on) {
    if (have_ms) .hmip_HighsCliqueTable_setPresolveFlag(ct, TRUE)
    # initialRowAndColPresolve
    for (row in pr$ridx) {
      if (pr$rowDeleted[row]) next
      .hs_row_presolve(pr, row)
      pr$changedRowFlag[row] <- FALSE
    }
    for (col in pr$cidx) {
      if (pr$colDeleted[col]) next
      if (pr$integrality[col] != 0L) {
        .hs_change_col_lower(pr, col, pr$cl[col])
        .hs_change_col_upper(pr, col, pr$cu[col])
      }
      .hs_col_presolve(pr, col)
      pr$changedColFlag[col] <- FALSE
    }
    parallel_calls <- 0L
    try_sparsify <- TRUE
    # the standalone presolve runs the bound part of runProbing as well
    try_probing <- TRUE
    ncliq_before <- -1L
    domcol_after <- FALSE
    last_print <- .Machine$integer.max
    red <- function() .hs_size_reduction(pr)
    big <- function() isTRUE(red() > 0.05)
    maybe_shrink <- function() {
      if (pr$numDeletedCols >= length(pr$cidx) %/% 2L || pr$numDeletedRows >= length(pr$ridx) %/% 2L)
        .hs_mip_shrink(pr)
    }
    nclq <- function() if (have_ms) .hmip_HighsCliqueTable_numCliques(ct) else 0L
    repeat {
      cur <- length(pr$cidx) - pr$numDeletedCols + length(pr$ridx) - pr$numDeletedRows
      if (cur < 0.85 * last_print) { last_print <- cur; .hs_mip_report(pr) }
      .hs_fast_presolve_loop(pr)
      .hs_store_size(pr)
      if (have_ms) .hs_apply_conflict_graph_substitutions(pr)
      .hs_aggregator(pr)
      if (big()) next
      if (try_sparsify) {
        nnz0 <- .hs_num_nonzeros(pr)
        .hs_sparsify(pr)
        nzred <- 100 * (1 - .hs_num_nonzeros(pr) / nnz0)
        if (isTRUE(nzred > 0)) {
          .hs_log_dev(pr, sprintf("Sparsify removed %.1f%% of nonzeros", nzred))
          .hs_fast_presolve_loop(pr)
        }
        try_sparsify <- FALSE
      }
      if (parallel_calls < 5L) {
        maybe_shrink()
        .hs_store_size(pr)
        .hs_parallel_rows_cols_mip(pr)
        parallel_calls <- parallel_calls + 1L
        if (big()) next
      }
      .hs_fast_presolve_loop(pr)
      .hs_strengthen_inequalities(pr)
      .hs_fast_presolve_loop(pr)
      if (ncliq_before == -1L) {
        ncliq_before <- nclq()
        .hs_store_size(pr)
        .hs_dominated_columns(pr)
        if (isTRUE(red() > 0)) .hs_fast_presolve_loop(pr)
        if (big()) next
      }
      if (try_probing) {
        .hs_detect_implied_integers(pr)
        .hs_store_size(pr)
        .hs_run_probing(pr)
        try_probing <- pr$probingContingent > pr$numProbed && (isTRUE(red() > 1) || pr$probingEarlyAbort)
        try_sparsify <- TRUE
        if (big() || try_probing) next
        .hs_fast_presolve_loop(pr)
      }
      if (have_ms && nclq() > ncliq_before && !domcol_after) {
        domcol_after <- TRUE
        .hs_store_size(pr)
        .hs_dominated_columns(pr)
        if (isTRUE(red() > 0)) .hs_fast_presolve_loop(pr)
        if (big()) next
      }
      break
    }
    .hs_mip_report(pr)
  }
  .hs_scale_mip(pr)
}

# HPresolve::run with a mipsolver (or standalone when pr$mipsolver is NULL)
.hs_mip_run <- function(pr) {
  ms <- pr$mipsolver; have_ms <- !is.null(ms)
  status <- tryCatch({
    .hs_mip_presolve_main(pr)
    NULL
  }, hs_presolve_abort = function(e) e$result)
  .hs_mip_flush_stack(pr)
  if (identical(status, "infeasible")) { pr$presolve_status <- "Infeasible"; return("Infeasible") }
  if (identical(status, "dualinfeasible")) {
    pr$presolve_status <- "UnboundedOrInfeasible"; return("UnboundedOrInfeasible")
  }
  .hs_mip_shrink(pr)
  if (have_ms) {
    md <- ms$mipdata_
    .hmip_HighsCliqueTable_setPresolveFlag(md$cliquetable, FALSE)
    .hmip_HighsCliqueTable_setMaxEntries(md$cliquetable, .hs_num_nonzeros(pr))
    .hmip_HighsDomain_addCutpool(md$domain, md$cutpool)
    .hmip_HighsDomain_addConflictPool(md$domain, md$conflictPool)
    if (md$numRestarts != 0L) {
      ps <- pr$postsolve_stack
      numcuts <- 0L
      for (i in rev(seq_along(pr$ridx)) - 1L) {
        if (.hmip_HighsPostsolveStack_getOrigRowIndex(ps, i) < ms$orig_model_$num_row_) break
        numcuts <- numcuts + 1L
        row <- pr$ridx[i + 1L]
        rowpos <- .hs_row_inorder(pr, row)
        integral <- pr$rowsizeInteger[row] + pr$rowsizeImplInt[row] == pr$rowsize[row] &&
          .hs_row_coefs_integral(pr, row, 1)
        .hmip_HighsCutPool_addCut(md$cutpool, ms, pr$ccur[pr$Acol[rowpos]], pr$Avalue[rowpos], length(rowpos),
                                  pr$ru[row], integral, TRUE, FALSE, FALSE)
        .hs_mark_row_deleted(pr, row)
        for (p in rowpos) .hs_unlink(pr, p)
      }
      if (numcuts > 0L) pr$ridx <- pr$ridx[seq_len(length(pr$ridx) - numcuts)]
    }
    .hs_mip_sync_model(pr)
  }
  numred <- if (!is.null(pr$postsolve_stack)) .hmip_HighsPostsolveStack_numReductions(pr$postsolve_stack)
            else length(pr$stack)
  if (!length(pr$cidx)) {
    if (have_ms) {
      if (pr$offset > ms$mipdata_$upper_limit) { pr$presolve_status <- "Infeasible"; return("Infeasible") }
      ms$mipdata_$lower_bound <- 0
    }
    pr$presolve_status <- "ReducedToEmpty"
    feas <- all(!(pr$rl[pr$ridx] > pr$ptol | pr$ru[pr$ridx] < -pr$ptol))
    return(if (feas) "Optimal" else "Infeasible")
  }
  pr$presolve_status <- if (numred > 0) "Reduced" else "NotReduced"
  "Notset"
}

# a standalone MIP presolve of a dense model (no clique table, so probing and
# the conflict graph rules do not run); the presolved model in the same
# shape as .hs_presolve returns
.hs_mip_presolve <- function(cost, A, row_lower, row_upper, col_lower, col_upper, integrality,
                             ptol = 1e-6, dtol = 1e-7, presolve = TRUE) {
  pr <- .hs_mip_new_state(cost, A, row_lower, row_upper, col_lower, col_upper, integrality, ptol, dtol)
  pr$presolve_on <- presolve
  ps <- .hmip_new_HighsPostsolveStack()
  .hmip_HighsPostsolveStack_initializeIndexMaps(ps, nrow(A), ncol(A))
  ps$in2orig_col <- seq_len(ncol(A)); ps$in2orig_row <- seq_len(nrow(A))
  pr$postsolve_stack <- ps
  st <- .hs_mip_run(pr)
  list(status = pr$presolve_status, model_status = st, pr = pr, cols = pr$cidx, rows = pr$ridx,
       lp = .hs_mip_model_list(pr), postsolve_stack = ps)
}

.hs_mip_model_list <- function(pr) {
  csc <- .hs_mip_to_csc(pr)
  list(cost = pr$cost[pr$cidx], col_lower = pr$cl[pr$cidx], col_upper = pr$cu[pr$cidx],
       row_lower = pr$rl[pr$ridx], row_upper = pr$ru[pr$ridx], integrality = pr$integrality[pr$cidx],
       offset = pr$offset, start = csc$start_, index = csc$index_, value = csc$value_)
}

# the HPresolve state for a MIP: the LP presolve state plus the MIP members
.hs_mip_new_state <- function(cost, A, row_lower, row_upper, col_lower, col_upper, integrality,
                              ptol = 1e-6, dtol = 1e-7, mipsolver = NULL) {
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
  pr$cidx <- seq_len(n); pr$ridx <- seq_len(m)
  pr$subst <- matrix(0L, 0L, 2L)
  .hs_mip_fields(pr, integrality, mipsolver)
  pr$mip <- TRUE
  pr$mip_feastol <- ptol
  pr$presolve_on <- TRUE
  pr$probingContingent <- 1000; pr$probingNumDelCol <- 0L; pr$numProbed <- 0L
  pr$numProbes <- numeric(n); pr$probingEarlyAbort <- FALSE
  pr$postsolve_stack <- NULL; pr$flushed <- 0L
  for (pos in seq_len(k)) .hs_link(pr, pos)
  pr$eqsize <- ifelse(row_lower == row_upper, pr$rowsize, NA_integer_)
  for (row in which(row_lower == row_upper))
    for (p in .hs_row_preorder(pr, row)) if (.hs_is_implied_free(pr, pr$Acol[p])) .hs_push_subst(pr, row, pr$Acol[p])
  pr
}

# ---------------------------------------------------------------------------
# HPresolve as a class for the MIP solver (group B1 calls these)
# ---------------------------------------------------------------------------

.hmip_new_HPresolve <- function() {
  self <- new.env(parent = emptyenv())
  self$pr <- NULL
  self$presolve_status_ <- "NotSet"
  self
}

# .hmip_HighsLp_dense(lp) is group B1's (optimization-scipy-highs-mip-solver.R)

# HPresolve::okSetInput(HighsMipSolver&, presolve_reduction_limit)
.hmip_HPresolve_okSetInput <- function(self, mipsolver, presolve_reduction_limit) {
  md <- mipsolver$mipdata_
  if (!identical(mipsolver$model_, md$presolvedModel)) {
    md$presolvedModel <- list2env(as.list(mipsolver$model_), parent = emptyenv())
    mipsolver$model_ <- md$presolvedModel
  } else {
    md$presolvedModel$col_lower_ <- md$domain$col_lower_
    md$presolvedModel$col_upper_ <- md$domain$col_upper_
  }
  lp <- md$presolvedModel
  opt <- mipsolver$options_mip_
  cost <- lp$col_cost_; off <- lp$offset_
  # HPresolve::presolve converts a maximization to a minimization
  if (identical(as.integer(lp$sense_), -1L)) { cost <- -cost; off <- -off }
  integ <- if (length(lp$integrality_)) lp$integrality_ else integer(lp$num_col_)
  pr <- .hs_mip_new_state(cost, .hmip_HighsLp_dense(lp), lp$row_lower_, lp$row_upper_, lp$col_lower_,
                          lp$col_upper_, integ, opt$mip_feasibility_tolerance %||% 1e-6,
                          opt$dual_feasibility_tolerance %||% 1e-7, mipsolver)
  pr$offset <- off
  pr$small_matrix_value <- opt$small_matrix_value %||% 1e-9
  pr$epsilon <- md$epsilon %||% 1e-9
  pr$presolve_on <- !identical(opt$presolve, "off")
  self$pr <- pr
  TRUE
}

# HPresolve::run(HighsPostsolveStack&) -> model status
.hmip_HPresolve_run <- function(self, postsolve_stack) {
  pr <- self$pr
  pr$postsolve_stack <- postsolve_stack
  postsolve_stack$in2orig_col <- postsolve_stack$origColIndex + 1L
  postsolve_stack$in2orig_row <- postsolve_stack$origRowIndex + 1L
  postsolve_stack$pr <- pr
  st <- .hs_mip_run(pr)
  postsolve_stack$pr <- NULL
  self$presolve_status_ <- pr$presolve_status
  st
}

.hmip_HPresolve_getPresolveStatus <- function(self) self$presolve_status_

# ---------------------------------------------------------------------------
# HighsPostsolveStack (the index maps, the record store and the primal undo)
# Records are stored with 1-based original indices; 0 stands for -1.
# ---------------------------------------------------------------------------

.hmip_new_HighsPostsolveStack <- function() {
  s <- new.env(parent = emptyenv())
  s$origNumCol <- 0L; s$origNumRow <- 0L
  s$origColIndex <- integer(); s$origRowIndex <- integer()
  s$linearlyTransformable <- logical()
  s$reductions <- list()
  s$pr <- NULL
  s
}

.hmip_HighsPostsolveStack_initializeIndexMaps <- function(s, numRow, numCol) {
  s$origNumRow <- as.integer(numRow); s$origNumCol <- as.integer(numCol)
  s$origRowIndex <- seq_len(numRow) - 1L; s$origColIndex <- seq_len(numCol) - 1L
  s$linearlyTransformable <- rep(TRUE, numCol)
  invisible(s)
}

.hmip_HighsPostsolveStack_compressIndexMaps <- function(s, newRowIndex, newColIndex) {
  keep <- newRowIndex != -1L
  s$origRowIndex <- s$origRowIndex[seq_along(newRowIndex)][keep][order(newRowIndex[keep])]
  keep <- newColIndex != -1L
  s$origColIndex <- s$origColIndex[seq_along(newColIndex)][keep][order(newColIndex[keep])]
  invisible(s)
}

.hmip_HighsPostsolveStack_getOrigColIndex <- function(s, col) s$origColIndex[col + 1L]
.hmip_HighsPostsolveStack_getOrigRowIndex <- function(s, row) s$origRowIndex[row + 1L]
.hmip_HighsPostsolveStack_getOrigNumCol <- function(s) s$origNumCol
.hmip_HighsPostsolveStack_getOrigNumRow <- function(s) s$origNumRow
.hmip_HighsPostsolveStack_numReductions <- function(s) {
  if (!is.null(s$pr)) .hs_mip_flush_stack(s$pr)
  length(s$reductions)
}
.hmip_HighsPostsolveStack_isColLinearlyTransformable <- function(s, col) {
  if (!is.null(s$pr)) .hs_mip_flush_stack(s$pr)
  s$linearlyTransformable[s$origColIndex[col + 1L] + 1L]
}

.hmip_HighsPostsolveStack_appendCutsToModel <- function(s, numCuts) {
  k <- length(s$origRowIndex)
  if (numCuts > 0L) {
    s$origRowIndex <- c(s$origRowIndex, s$origNumRow + seq_len(numCuts) - 1L)
    s$origNumRow <- s$origNumRow + as.integer(numCuts)
  }
  invisible(s)
}

.hmip_HighsPostsolveStack_removeCutsFromModel <- function(s, numCuts) {
  s$origNumRow <- s$origNumRow - as.integer(numCuts)
  n <- length(s$origRowIndex)
  while (n > 0L && s$origRowIndex[n] >= s$origNumRow) n <- n - 1L
  s$origRowIndex <- s$origRowIndex[seq_len(n)]
  invisible(s)
}

# translate the records the presolve pushed since the last flush from the
# presolve input indices to the original indices
.hs_mip_flush_stack <- function(pr) {
  s <- pr$postsolve_stack
  if (is.null(s)) return(invisible())
  k <- length(pr$stack)
  if (k <= pr$flushed) return(invisible())
  cm <- s$in2orig_col; rm <- s$in2orig_row
  mc <- function(v) cm[v]
  mr <- function(v) ifelse(v == 0L, 0L, rm[pmax(v, 1L)])
  for (i in (pr$flushed + 1L):k) {
    r <- pr$stack[[i]]
    for (f in intersect(c("col", "colsubst"), names(r))) r[[f]] <- mc(r[[f]])
    if (!is.null(r$row)) r$row <- mr(r$row)
    if (!is.null(r$addedeq)) r$addedeq <- mr(r$addedeq)
    if (r$type == "duplicaterow") r$dup <- mr(r$dup)
    else if (!is.null(r$dup)) r$dup <- mc(r$dup)
    if (!is.null(r$colv)) r$colv$index <- mr(r$colv$index)
    if (!is.null(r$rowv)) r$rowv$index <- mc(r$rowv$index)
    if (!is.null(r$addrows)) r$addrows$index <- mr(r$addrows$index)
    if (r$type == "duplicatecol") {
      s$linearlyTransformable[r$col] <- FALSE
      s$linearlyTransformable[r$dup] <- FALSE
    }
    s$reductions[[length(s$reductions) + 1L]] <- r
  }
  pr$flushed <- k
  invisible()
}

# HighsPostsolveStack::undoPrimal: the column values of the original model
# from those of the presolved model (no row values and no duals are given, so
# the row value vector is the zero vector undoIterateBackwards leaves)
.hmip_HighsPostsolveStack_undoPrimal <- function(s, options, col_value) {
  if (is.list(col_value)) col_value <- col_value$col_value
  tol <- options$mip_feasibility_tolerance %||% 1e-6
  ptol <- options$primal_feasibility_tolerance %||% 1e-7
  expand <- function(v, idx, size) {
    out <- c(v, numeric(max(0L, size - length(v))))[seq_len(size)]
    for (i in rev(seq_along(idx))) out[idx[i] + 1L] <- out[i]
    out
  }
  x <- expand(col_value, s$origColIndex, s$origNumCol)
  rv <- numeric(s$origNumRow)
  for (r in rev(s$reductions)) {
    switch(r$type,
      lineartransform = {
        x[r$col] <- x[r$col] * r$scale
        x[r$col] <- x[r$col] + r$constant
      },
      fixedcol = { x[r$col] <- r$fix },
      doubleton = {
        x[r$colsubst] <- .cd_val(.cd_div(.cd_add(-.cd_mul(c(r$coef, 0), x[r$col]), r$rhs), r$coefsubst))
      },
      freecol = {
        coef <- 0; val <- c(0, 0)
        for (k in seq_along(r$rowv$index)) {
          j <- r$rowv$index[k]
          if (j == r$col) coef <- r$rowv$value[k] else val <- .cd_add(val, r$rowv$value[k] * x[j])
        }
        if (r$row <= length(rv)) rv[r$row] <- .cd_val(.cd_add(val, coef * x[r$col]))
        x[r$col] <- .cd_val(.cd_div(.cd_add(-val, r$rhs), coef))
      },
      forcingcol = {
        nb_row <- 0L; val <- r$bound
        dir <- if (r$at_inf_upper) 1 else -1
        for (k in seq_along(r$colv$index)) {
          i <- r$colv$index[k]
          if (i <= length(rv)) {
            v <- rv[i] / r$colv$value[k]
            if (dir * v > dir * val) { nb_row <- i; val <- v }
          }
        }
        if (nb_row != 0L && isTRUE(r$integral)) val <- dir * ceiling(dir * val - tol)
        x[r$col] <- val
      },
      forcingcolrow = {
        val <- c(r$rhs, 0)
        for (k in seq_along(r$rowv$index)) val <- .cd_add(val, -(r$rowv$value[k] * x[r$rowv$index[k]]))
        if (r$row <= length(rv)) rv[r$row] <- .cd_val(val)
      },
      duplicatecol = {
        merge <- x[r$col]
        di <- isTRUE(r$di); ci <- isTRUE(r$ci)
        x[r$col] <- if (r$cl != -Inf) r$cl else min(0, r$cu)
        x[r$dup] <- .cd_val(.cd_div(.cd_add(c(merge, 0), -x[r$col]), r$scale))
        recompute <- FALSE
        if (x[r$dup] > r$du) { x[r$dup] <- r$du; recompute <- TRUE }
        else if (x[r$dup] < r$dl) { x[r$dup] <- r$dl; recompute <- TRUE }
        else if (di && .hs_fractionality(x[r$dup]) > 1e-6) { x[r$dup] <- floor(x[r$dup]); recompute <- TRUE }
        if (recompute) {
          x[r$col] <- merge - r$scale * x[r$dup]
          if (!di && ci) {
            x[r$col] <- ceiling(x[r$col] - 1e-6)
            x[r$dup] <- .cd_val(.cd_div(.cd_add(c(merge, 0), -x[r$col]), r$scale))
          }
        }
        if (.hs_dupcol_split_error(r, x[r$col], x[r$dup], merge, ptol)) {
          v <- .hs_dupcol_undo_fix(r, merge, ptol)
          x[r$col] <- v[1L]; x[r$dup] <- v[2L]
        }
      },
      NULL)
  }
  x
}

# HighsPostsolveStack::getReducedPrimalSolution
.hmip_HighsPostsolveStack_getReducedPrimalSolution <- function(s, origPrimalSolution) {
  x <- origPrimalSolution
  for (r in s$reductions) {
    if (r$type == "duplicatecol") x[r$col] <- x[r$col] + r$scale * x[r$dup]
    else if (r$type == "lineartransform") { x[r$col] <- x[r$col] - r$constant; x[r$col] <- x[r$col] / r$scale }
  }
  x[s$origColIndex + 1L]
}

# ---------------------------------------------------------------------------
# HPresolve::detectParallelRowsAndCols with a mipsolver: weakly implied
# bounds suffice, integer columns restrict the merge scale, rows with
# singleton columns are handled by equalityRowAddition. The hashes and the
# bucket order are those of the LP version (.hs_parallel_rows_cols).
# ---------------------------------------------------------------------------

# fractionality(input, &intval) of util/HighsUtils.h
.hs_frac_int <- function(v) { r <- .hmip_round(v); c(abs(v - r), r) }

# HighsPostsolveStack::DuplicateColumn::okMerge
.hs_dupcol_ok_merge <- function(scale, cl, cu, dl, du, x_int, y_int, tol) {
  x_lo <- if (x_int) ceiling(cl - tol) else cl
  x_up <- if (x_int) floor(cu + tol) else cu
  y_lo <- if (y_int) ceiling(dl - tol) else dl
  y_up <- if (y_int) floor(du + tol) else du
  x_len <- x_up - x_lo; y_len <- y_up - y_lo
  ok <- scale != 0
  a <- abs(scale)
  if (x_int) {
    if (y_int) {
      if (.hmip_fractionality(scale) > tol) ok <- FALSE
      if (a > x_len + 1 + tol) ok <- FALSE
    } else {
      if (y_len == 0) ok <- FALSE else if (a < 1 / y_len) ok <- FALSE
    }
  } else if (y_int) {
    if (a > x_len) ok <- FALSE
  }
  ok
}

.hs_parallel_rows_cols_mip <- function(pr) {
  small <- pr$small_matrix_value
  ptol <- pr$ptol
  live <- which(pr$Avalue != 0)
  row_max <- numeric(pr$m); row_max_i <- integer(pr$m)
  col_max <- numeric(pr$n); col_max_i <- integer(pr$n)
  nsing <- integer(pr$m)
  for (i in live) {
    r <- pr$Arow[i]; c <- pr$Acol[i]; v <- pr$Avalue[i]
    if (pr$colsize[c] == 1L) {
      col_max[c] <- v
      nsing[r] <- nsing[r] + 1L
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
  row_hash <- lapply(pr$rowsize - nsing, function(v) .hs_u64(0, v))
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
  isint <- function(j) pr$integrality[j] == 1L
  dec_sing <- function(j) if (pr$colsize[j] == 1L) {
    r <- pr$Arow[pr$colhead[j]]; nsing[r] <<- nsing[r] - 1L
  }
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
      chk_col <- TRUE; chk_dup <- TRUE
      if (isint(i) && isint(cand)) {
        if (abs(col_max[i]) < abs(col_max[cand])) { col <- i; dup <- cand } else { col <- cand; dup <- i }
        fr <- .hs_frac_int(col_max[dup] / col_max[col])
        if (fr[1L] > small) next
        scale <- fr[2L]
        if (scale != 1) chk_dup <- FALSE
      } else if (isint(i)) {
        col <- i; dup <- cand
        scale <- col_max[dup] / col_max[col]
        chk_col <- FALSE
      } else {
        col <- cand; dup <- i
        scale <- col_max[dup] / col_max[col]
        chk_col <- !isint(cand)
      }
      upper_inf <- function() chk_col && (if (scale > 0) .hs_is_upper_implied(pr, col) else .hs_is_lower_implied(pr, col))
      lower_inf <- function() chk_col && (if (scale > 0) .hs_is_lower_implied(pr, col) else .hs_is_upper_implied(pr, col))
      dup_upper_inf <- function() chk_dup && .hs_is_upper_implied(pr, dup)
      dup_lower_inf <- function() chk_dup && .hs_is_lower_implied(pr, dup)
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
      if (case == "merge" && isint(col)) {
        illegal <- if (!isint(dup)) abs(scale * (cu[dup] - cl[dup])) < 1 - ptol
                   else abs(scale) > cu[col] - cl[col] + 1 + ptol
        if (illegal) next
      }
      parallel <- TRUE
      for (p in .hs_col_slots(pr, col)) {
        q <- .hs_find_nonzero(pr, pr$Arow[p], dup)
        parallel <- q != 0L && abs(.cd_val(.cd_add(-.cd_two_product(scale, pr$Avalue[p]), pr$Avalue[q]))) <= small
        if (!parallel) break
      }
      if (!parallel) next
      switch(case,
        duplower = { del <- dup; dec_sing(dup); .hs_fix_col_to(pr, dup, FALSE) },
        dupupper = { del <- dup; dec_sing(dup); .hs_fix_col_to(pr, dup, TRUE) },
        collower = { del <- col; dec_sing(col); .hs_fix_col_to(pr, col, FALSE) },
        colupper = { del <- col; dec_sing(col); .hs_fix_col_to(pr, col, TRUE) },
        merge = {
          # upstream merges even when okMerge refuses the postsolve record
          # (its break needs debug_report); ported as is
          if (.hs_dupcol_ok_merge(scale, cl[col], cu[col], cl[dup], cu[dup], isint(col), isint(dup),
                                  pr$mip_feastol))
            .hs_ps_push(pr, list(type = "duplicatecol", col = col, dup = dup, scale = scale,
                                 cl = cl[col], cu = cu[col], dl = cl[dup], du = cu[dup],
                                 ci = isint(col), di = isint(dup)))
          rsir <- !isint(dup) && isint(col)
          if (rsir) pr$integrality[col] <- 0L
          .hs_mark_col(pr, col)
          dec_sing(dup)
          if (scale > 0) {
            ml <- cl[col] + scale * cl[dup]; mu <- cu[col] + scale * cu[dup]
          } else {
            ml <- cl[col] + scale * cu[dup]; mu <- cu[col] + scale * cl[dup]
          }
          .hs_change_col_lower(pr, col, ml)
          .hs_change_col_upper(pr, col, mu)
          .hs_mark_col_deleted(pr, dup)
          p <- pr$colhead[dup]
          while (p != 0L) {
            r <- pr$Arow[p]
            if (rsir) pr$rowsizeInteger[r] <- pr$rowsizeInteger[r] - 1L
            nx <- pr$Anext[p]
            .hs_unlink(pr, p)
            .hs_reinsert_equation(pr, r)
            p <- nx
          }
          pr$cost[dup] <- 0
          del <- dup
          .hs_reset_col_impl(pr, col, 0L)
          if (rsir && pr$integrality[dup] == 4L && .hs_is_implied_integer(pr, col))
            .hs_convert_implied_integer(pr, col, 0L, TRUE)
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
    grp <- .hs_bucket_group(buckets, row_key[i])
    ns <- nsing[i]
    del <- 0L; last <- 0L
    rowpos <- if (length(grp)) .hs_row_inorder(pr, i) else integer()
    for (gi in seq_along(grp)) {
      cand <- grp[gi]; last <- gi
      nsc <- nsing[cand]
      if (pr$rowsize[i] - ns != pr$rowsize[cand] - nsc) next
      if (nsc > 1L || ns > 1L) {
        if ((ns != 0L || !.hs_is_equation(pr, i)) && (nsc != 0L || !.hs_is_equation(pr, cand))) next
      } else if (nsc != ns) {
        if (!.hs_is_equation(pr, i) && !.hs_is_equation(pr, cand)) next
      }
      scale <- row_max[cand] / row_max[i]
      parallel <- TRUE
      for (p in rowpos) {
        if (pr$colsize[pr$Acol[p]] == 1L) next
        q <- .hs_find_nonzero(pr, cand, pr$Acol[p])
        parallel <- q != 0L && abs(.cd_val(.cd_add(-.cd_two_product(scale, pr$Avalue[p]), pr$Avalue[q]))) <= small
        if (!parallel) break
      }
      if (!parallel) next
      if (ns == 0L && nsc == 0L) {
        ut <- FALSE; lt <- FALSE
        if (scale > 0) {
          nu <- pr$ru[i] * scale; nl <- pr$rl[i] * scale
        } else {
          nl <- pr$ru[i] * scale; nu <- pr$rl[i] * scale
        }
        if (nu < pr$ru[cand]) {
          if (nu < pr$rl[cand] - ptol) .hs_pre_abort("infeasible")
          if (nu <= pr$rl[cand] + ptol) nu <- pr$rl[cand]
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
          if (nl > pr$ru[cand] + ptol) .hs_pre_abort("infeasible")
          if (nl >= pr$ru[cand] - ptol) nl <- pr$ru[cand]
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
      } else if (.hs_is_equation(pr, i)) {
        .hs_equality_row_addition(pr, i, cand, -scale, rowpos)
        del <- cand
      } else if (.hs_is_equation(pr, cand)) {
        .hs_equality_row_addition(pr, cand, i, -row_max[i] / row_max[cand], .hs_row_preorder(pr, cand))
        del <- i
      }
    }
    if (del != 0L) {
      if (del != i) buckets[[row_key[i]]] <- grp[-last]
    } else {
      buckets[[row_key[i]]] <- .hs_bucket_insert(grp, i)
    }
  }
}

# the mipsolver branch of HPresolve::updateColImpliedBounds: an implied bound
# of an integer column tightens the column bound, and a cut row (one beyond
# the original rows) gives a column bound but no implied bound
.hs_mip_implied_col_upper <- function(pr, row, col, ib, thr) {
  if (pr$integrality[col] != 0L && ib < pr$cu[col] - pr$ptol) .hs_change_col_upper(pr, col, ib)
  if (.hs_mip_is_cut_row(pr, row)) {
    if (ib < pr$cu[col] - thr) .hs_change_col_upper(pr, col, ib)
    ib <- Inf
  }
  ib
}
.hs_mip_implied_col_lower <- function(pr, row, col, ib, thr) {
  if (pr$integrality[col] != 0L && ib > pr$cl[col] + pr$ptol) .hs_change_col_lower(pr, col, ib)
  if (.hs_mip_is_cut_row(pr, row)) {
    if (ib > pr$cl[col] + thr) .hs_change_col_lower(pr, col, ib)
    ib <- -Inf
  }
  ib
}
.hs_mip_is_cut_row <- function(pr, row) {
  ms <- pr$mipsolver; s <- pr$postsolve_stack
  if (is.null(ms) || is.null(s)) return(FALSE)
  s$in2orig_row[row] - 1L >= ms$orig_model_$num_row_
}
