# The basis factorization of HiGHS 1.12.0 (util/HFactor.cpp, util/HFactor.h,
# util/HVectorBase.cpp), as SciPy 1.18.1 vendors it: HFactor::build
# (buildSimple, buildKernel, buildHandleRankDeficiency, buildFinish), the
# solves ftranL, ftranU, btranU and btranL with their sparse and hyper
# sparse paths, and the Forrest-Tomlin update (updateFT, ftranFT, btranFT)
# that HEkk configures. MIT license notice in
# R/optimization-scipy-highs-util.R.
#
# Vectors follow HVector: a list of the dense array, the index list of the
# nonzeros (in the order HiGHS keeps it) and count (-1 when the index list
# is not maintained). Rows, pivots and positions are 1-based.

.hs_tiny <- 1e-14
.hs_zero <- 0x1.dee7a4ad4b81fp-167

# cols holds the basic columns in basic_index order, each a list of the row
# indices (1-based, in matrix storage order), the values and logical = TRUE
# for a slack column. Returns an environment with permute (the pivot row of
# each basic column), the rank deficiency, the columns and rows left
# without a pivot and, for a square basis, the L and U factors.
.hs_factor_build <- function(num_row, cols, pivot_threshold = 0.1, pivot_tolerance = 1e-10) {
  num_basic <- length(cols)
  fa <- new.env(parent = emptyenv())
  fa$num_row <- num_row
  permute <- integer(num_basic)
  mrb <- integer(num_row)
  work <- integer()
  # the pivots in the order they are taken: L columns, U columns
  piv_row <- integer(); piv_val <- numeric()
  l_idx <- list(); l_val <- list(); u_idx <- list(); u_val <- list()
  push <- function(r, v, li = integer(), lv = numeric(), ui = integer(), uv = numeric()) {
    k <- length(piv_row) + 1L
    piv_row[k] <<- r; piv_val[k] <<- v
    l_idx[k] <<- list(li); l_val[k] <<- list(lv); u_idx[k] <<- list(ui); u_val[k] <<- list(uv)
  }
  # the synthetic clock of the build (build_synthetic_tick)
  tick <- 0; bcount <- 0
  # 1. unit columns
  for (c in seq_len(num_basic)) {
    col <- cols[[c]]
    r <- 0L
    if (isTRUE(col$logical)) {
      if (mrb[col$index] >= 0L) r <- col$index else {
        mrb[col$index] <- mrb[col$index] + 1L
        work <- c(work, c); bcount <- bcount + 1
      }
    } else if (length(col$index) == 1L && col$value == 1 && mrb[col$index] >= 0L) {
      r <- col$index
    } else {
      mrb[col$index] <- mrb[col$index] + 1L
      work <- c(work, c); bcount <- bcount + length(col$index)
    }
    if (r > 0L) {
      permute[c] <- r
      push(r, 1)
      mrb[r] <- -num_basic
    }
  }
  tick <- tick + bcount * 60 + (num_row - length(work)) * 80
  t2_search <- 0; t2_l <- 0; t2_u <- 0; t2_p <- length(work)
  # 2. singletons
  while (length(work)) {
    last <- work
    work <- integer()
    for (c in last) {
      idx <- cols[[c]]$index; val <- cols[[c]]$value
      if (isTRUE(cols[[c]]$logical)) val <- 1
      t2_search <- t2_search + length(idx)
      pivot_k <- 0L; found_row_singleton <- FALSE; count <- 0L
      for (k in seq_along(idx)) {
        if (mrb[idx[k]] == 1L) {
          pivot_k <- k; found_row_singleton <- TRUE
          break
        }
        if (mrb[idx[k]] > 1L) {
          pivot_k <- k; count <- count + 1L
        }
      }
      if (found_row_singleton) {
        pm <- 1 / val[pivot_k]
        li <- integer(); lv <- numeric(); ui <- integer(); uv <- numeric()
        for (k in seq_along(idx)[-pivot_k]) {
          i <- idx[k]
          if (mrb[i] > 0L) {
            li <- c(li, i); lv <- c(lv, val[k] * pm)
          } else {
            ui <- c(ui, i); uv <- c(uv, val[k])
          }
          mrb[i] <- mrb[i] - 1L
        }
        mrb[idx[pivot_k]] <- 0L
        permute[c] <- idx[pivot_k]
        push(idx[pivot_k], val[pivot_k], li, lv, ui, uv)
        t2_l <- t2_l + length(li); t2_u <- t2_u + length(ui)
      } else if (count == 1L) {
        mrb[idx[pivot_k]] <- 0L
        permute[c] <- idx[pivot_k]
        push(idx[pivot_k], val[pivot_k], ui = idx[-pivot_k], uv = val[-pivot_k])
        t2_u <- t2_u + length(idx) - 1
      } else {
        work <- c(work, c)
      }
    }
    if (length(last) == length(work)) break
  }
  t2_p <- t2_p - length(work)
  tick <- tick + t2_search * 20 + (t2_p + t2_l + t2_u) * 80
  # 3. the kernel: active column entries, the entries in pivoted rows
  # (colStoreN, read back last first), row lists and count links
  mc_idx <- rep(list(integer()), num_basic); mc_val <- rep(list(numeric()), num_basic)
  mn_idx <- rep(list(integer()), num_basic); mn_val <- rep(list(numeric()), num_basic)
  mc_min <- numeric(num_basic)
  mr <- rep(list(integer()), num_row)
  max_count <- max(num_row, num_basic)
  # links: first[count + 1], next and last per node; -1 is no node and a
  # last of -2 - count marks the head of the list of that count
  cfirst <- rep(-1L, max_count + 2L); cnext <- integer(num_basic); clast <- integer(num_basic)
  rfirst <- rep(-1L, max_count + 2L); rnext <- integer(num_row); rlast <- integer(num_row)
  cadd <- function(j, count) {
    mover <- cfirst[count + 1L]
    clast[j] <<- -2L - count; cnext[j] <<- mover; cfirst[count + 1L] <<- j
    if (mover >= 0L) clast[mover] <<- j
  }
  cdel <- function(j) {
    xl <- clast[j]; xn <- cnext[j]
    if (xl >= 0L) cnext[xl] <<- xn else cfirst[-xl - 2L + 1L] <<- xn
    if (xn >= 0L) clast[xn] <<- xl
  }
  radd <- function(i, count) {
    mover <- rfirst[count + 1L]
    rlast[i] <<- -2L - count; rnext[i] <<- mover; rfirst[count + 1L] <<- i
    if (mover >= 0L) rlast[mover] <<- i
  }
  rdel <- function(i) {
    xl <- rlast[i]; xn <- rnext[i]
    if (xl >= 0L) rnext[xl] <<- xn else rfirst[-xl - 2L + 1L] <<- xn
    if (xn >= 0L) rlast[xn] <<- xl
  }
  col_delete <- function(j, i) {
    k <- match(i, mc_idx[[j]]); n <- length(mc_idx[[j]])
    v <- mc_val[[j]][k]
    mc_idx[[j]][k] <<- mc_idx[[j]][n]; mc_val[[j]][k] <<- mc_val[[j]][n]
    mc_idx[[j]] <<- mc_idx[[j]][-n]; mc_val[[j]] <<- mc_val[[j]][-n]
    v
  }
  row_delete <- function(j, i) {
    k <- match(j, mr[[i]]); n <- length(mr[[i]])
    mr[[i]][k] <<- mr[[i]][n]
    mr[[i]] <<- mr[[i]][-n]
  }
  fix_max <- function(j) mc_min[j] <<- max(0, abs(mc_val[[j]])) * pivot_threshold
  for (i in seq_len(num_row)) if (mrb[i] > 0L) radd(i, mrb[i])
  mr_countx <- 2 * sum(mrb[mrb > 0L]); mc_countx <- 0
  for (c in work) {
    col <- cols[[c]]
    val <- if (isTRUE(col$logical)) 1 else col$value
    mc_countx <- mc_countx + 2 * length(col$index)
    for (k in seq_along(col$index)) {
      i <- col$index[k]
      if (mrb[i] > 0L) {
        mc_idx[[c]] <- c(mc_idx[[c]], i); mc_val[[c]] <- c(mc_val[[c]], val[k])
        mr[[i]] <- c(mr[[i]], c)
      } else {
        mn_idx[[c]] <- c(mn_idx[[c]], i); mn_val[[c]] <- c(mn_val[[c]], val[k])
      }
    }
    fix_max(c)
    cadd(c, length(mc_idx[[c]]))
  }
  nwork <- length(work)
  tick <- tick + (num_row + nwork + mc_countx) * 40 + mr_countx * 20
  fake_search <- 0; fake_fill <- 0; fake_eliminate <- 0
  rank_deficiency <- 0L
  while (nwork > 0L) {
    nwork <- nwork - 1L
    jp <- -1L; ip <- -1L
    search_limit <- min(nwork, 8L); search_count <- 0L
    merit_limit <- 1.0 * num_basic * num_row; merit_pivot <- merit_limit
    found <- FALSE
    if (cfirst[2L] != -1L) {
      jp <- cfirst[2L]; ip <- mc_idx[[jp]][1L]; found <- TRUE
    }
    if (!found && rfirst[2L] != -1L) {
      ip <- rfirst[2L]; jp <- mr[[ip]][1L]; found <- TRUE
    }
    count <- 2L
    while (!found && count <= max_count) {
      if (count <= num_row) {
        j <- cfirst[count + 1L]
        while (j != -1L) {
          mp <- mc_min[j]
          for (k in seq_along(mc_idx[[j]])) {
            if (abs(mc_val[[j]][k]) >= mp) {
              i <- mc_idx[[j]][k]
              rc <- length(mr[[i]])
              ml <- 1.0 * (count - 1L) * (rc - 1L)
              if (merit_pivot > ml) {
                merit_pivot <- ml; jp <- j; ip <- i
                found <- found || rc < count
              }
            }
          }
          if (search_count >= search_limit && merit_pivot < merit_limit) found <- TRUE
          search_count <- search_count + 1L
          if (found) break
          fake_search <- fake_search + count
          j <- cnext[j]
        }
      }
      # the rows of this count are searched even after a column pivot was
      # accepted: the first one can still improve the merit
      if (count <= num_basic) {
        i <- rfirst[count + 1L]
        while (i != -1L) {
          for (j in mr[[i]]) {
            cc <- length(mc_idx[[j]])
            ml <- 1.0 * (count - 1L) * (cc - 1L)
            if (ml < merit_pivot) {
              if (abs(mc_val[[j]][match(i, mc_idx[[j]])]) >= mc_min[j]) {
                merit_pivot <- ml; jp <- j; ip <- i
                found <- found || cc <= count
              }
            }
          }
          if (search_count >= search_limit && merit_pivot < merit_limit) found <- TRUE
          search_count <- search_count + 1L
          if (found) break
          i <- rnext[i]
        }
        fake_search <- fake_search + count
      }
      count <- count + 1L
    }
    if (ip < 0L) {
      rank_deficiency <- nwork + 1L
      break
    }
    pm <- col_delete(jp, ip)
    row_delete(jp, ip)
    cdel(jp)
    rdel(ip)
    if (abs(pm) < pivot_tolerance) {
      if (length(mr[[ip]]) == 0L) {
        cadd(jp, length(mc_idx[[jp]]))
      } else {
        # zeroCol, whose clinkDel repeats the one above on stale links
        for (i in mc_idx[[jp]]) {
          row_delete(jp, i)
          rdel(i)
          radd(i, length(mr[[i]]))
        }
        cdel(jp)
        mc_idx[[jp]] <- integer(); mc_val[[jp]] <- numeric()
        mn_idx[[jp]] <- integer(); mn_val[[jp]] <- numeric()
        radd(ip, length(mr[[ip]]))
      }
      nwork <- nwork + 1L
      next
    }
    permute[jp] <- ip
    mwz_rows <- mc_idx[[jp]]
    mwz <- numeric(num_row); mark <- logical(num_row)
    before <- integer(num_row)
    for (k in seq_along(mwz_rows)) {
      i <- mwz_rows[k]
      mwz[i] <- mc_val[[jp]][k] / pm
      mark[i] <- TRUE
      before[i] <- length(mr[[i]])
      row_delete(jp, i)
    }
    push(ip, pm, mwz_rows, mwz[mwz_rows], rev(mn_idx[[jp]]), rev(mn_val[[jp]]))
    fake_fill <- fake_fill + 2 * length(mwz_rows) + length(mn_idx[[jp]])
    for (ic in mr[[ip]]) {
      my_count <- length(mc_idx[[ic]])
      my_pivot <- col_delete(ic, ip)
      mn_idx[[ic]] <- c(mn_idx[[ic]], ip); mn_val[[ic]] <- c(mn_val[[ic]], my_pivot)
      nfill <- length(mwz_rows); ncancel <- 0L
      for (k in seq_along(mc_idx[[ic]])) {
        i <- mc_idx[[ic]][k]
        if (mark[i]) {
          mark[i] <- FALSE
          nfill <- nfill - 1L
          v <- mc_val[[ic]][k] - my_pivot * mwz[i]
          if (abs(v) < .hs_tiny) {
            v <- 0; ncancel <- ncancel + 1L
          }
          mc_val[[ic]][k] <- v
        }
      }
      if (ncancel > 0L) {
        zero <- mc_val[[ic]] == 0
        for (i in mc_idx[[ic]][zero]) row_delete(ic, i)
        mc_idx[[ic]] <- mc_idx[[ic]][!zero]; mc_val[[ic]] <- mc_val[[ic]][!zero]
      }
      if (nfill > 0L) {
        add <- mwz_rows[mark[mwz_rows]]
        mc_idx[[ic]] <- c(mc_idx[[ic]], add); mc_val[[ic]] <- c(mc_val[[ic]], -my_pivot * mwz[add])
        for (i in add) mr[[i]] <- c(mr[[i]], ic)
      }
      fake_eliminate <- fake_eliminate + length(mwz_rows) + nfill * 2
      mark[mwz_rows] <- TRUE
      fix_max(ic)
      if (my_count != length(mc_idx[[ic]])) {
        cdel(ic)
        cadd(ic, length(mc_idx[[ic]]))
      }
    }
    for (i in mwz_rows) {
      if (before[i] != length(mr[[i]])) {
        rdel(i)
        radd(i, length(mr[[i]]))
      }
    }
  }
  if (rank_deficiency == 0L) tick <- tick + fake_search * 20 + fake_fill * 160 + fake_eliminate * 80
  # buildHandleRankDeficiency
  col_np <- integer(); row_np <- integer()
  rd <- rank_deficiency
  if (rank_deficiency > 0L || num_basic < num_row) {
    pivoted <- integer(num_row)
    for (c in seq_len(num_basic)) if (permute[c] > 0L) pivoted[permute[c]] <- c else col_np <- c(col_np, c)
    if (num_basic < num_row) col_np <- c(col_np, num_basic + seq_len(num_row - num_basic))
    row_np <- which(pivoted == 0L)
    if (num_row < num_basic) row_np <- c(row_np, num_row + seq_len(num_basic - num_row))
    rd <- length(col_np)
    row_rank_deficiency <- rd - max(num_basic - num_row, 0L)
    for (k in seq_len(min(length(col_np), length(row_np)))) {
      if (col_np[k] <= num_basic) permute[col_np[k]] <- row_np[k]
      if (k <= row_rank_deficiency) push(row_np[k], 1)
    }
  }
  fa$build_tick <- tick
  fa$permute <- permute
  fa$rank_deficiency <- rank_deficiency
  fa$col_with_no_pivot <- col_np
  fa$row_with_no_pivot <- row_np
  if (num_basic != num_row) return(fa)
  # buildFinish
  tick <- tick + num_row * 80 + (sum(lengths(l_idx)) + sum(lengths(u_idx))) * 60
  fa$build_tick <- tick
  fa$l_piv <- piv_row; fa$l_idx <- l_idx; fa$l_val <- l_val
  lookup <- integer(num_row); lookup[piv_row] <- seq_len(num_row)
  fa$l_lookup <- lookup
  lr_idx <- rep(list(integer()), num_row); lr_val <- rep(list(numeric()), num_row)
  for (i in seq_len(num_row)) for (k in seq_along(l_idx[[i]])) {
    t <- lookup[l_idx[[i]][k]]
    lr_idx[[t]] <- c(lr_idx[[t]], piv_row[i]); lr_val[[t]] <- c(lr_val[[t]], l_val[[i]][k])
  }
  fa$lr_idx <- lr_idx; fa$lr_val <- lr_val
  fa$u_piv <- piv_row; fa$u_pv <- piv_val; fa$u_idx <- u_idx; fa$u_val <- u_val
  fa$u_lookup <- lookup
  ur_idx <- rep(list(integer()), num_row); ur_val <- rep(list(numeric()), num_row)
  for (i in seq_len(num_row)) for (k in seq_along(u_idx[[i]])) {
    t <- lookup[u_idx[[i]][k]]
    ur_idx[[t]] <- c(ur_idx[[t]], piv_row[i]); ur_val[[t]] <- c(ur_val[[t]], u_val[[i]][k])
  }
  fa$ur_idx <- ur_idx; fa$ur_val <- ur_val
  fa$pf_piv <- integer(); fa$pf_idx <- list(); fa$pf_val <- list()
  fa
}

# HVector helpers
.hs_hv <- function(array, index = which(array != 0)) list(array = array, index = index, count = length(index), tick = 0)

.hs_hv_tight <- function(hv) {
  if (hv$count < 0L) {
    hv$array[abs(hv$array) < .hs_tiny] <- 0
  } else {
    keep <- abs(hv$array[hv$index]) >= .hs_tiny
    hv$array[hv$index[!keep]] <- 0
    hv$index <- hv$index[keep]; hv$count <- length(hv$index)
  }
  hv
}

.hs_hv_reindex <- function(hv) {
  if (hv$count >= 0L && hv$count <= length(hv$array) * 0.1) return(hv)
  hv$index <- which(hv$array != 0); hv$count <- length(hv$index)
  hv
}

# HighsSparseMatrix::collectAj into hv for variable v with multiplier
.hs_hv_collect <- function(ek, hv, v, multiplier) {
  if (v <= ek$n) {
    vals <- ek$A[, v]; rows <- which(vals != 0); vals <- vals[rows]
  } else {
    rows <- v - ek$n; vals <- 1
  }
  for (k in seq_along(rows)) {
    i <- rows[k]
    x0 <- hv$array[i]
    x1 <- x0 + multiplier * vals[k]
    if (x0 == 0) hv$index <- c(hv$index, i)
    hv$array[i] <- if (abs(x1) < .hs_tiny) .hs_zero else x1
  }
  hv$count <- length(hv$index)
  hv
}

# HVectorBase::norm2 in index order
.hs_hv_norm2 <- function(hv) {
  s <- 0
  for (i in hv$index) s <- s + hv$array[i] * hv$array[i]
  s
}

# the sparse solve over pivots `order` with columns (idx, val): each pivot
# value is final when reached and is scattered into its column
.hs_solve_sparse <- function(hv, order, piv, pv, idx, val, m = Inf) {
  a <- hv$array
  out <- integer()
  ft <- 0
  for (p in order) {
    r <- piv[p]
    if (r == 0L) next
    x <- a[r]
    if (abs(x) > .hs_tiny) {
      if (!is.null(pv)) {
        x <- x / pv[p]
        a[r] <- x
      }
      out <- c(out, r)
      ci <- idx[[p]]
      if (p > m) ft <- ft + length(ci)
      if (length(ci)) a[ci] <- a[ci] - x * val[[p]]
    } else a[r] <- 0
  }
  list(array = a, index = out, count = length(out), tick = hv$tick %||% 0, ft = ft)
}

# solveHyper: depth first search from the nonzeros for the pivots reached,
# then the solve in topological order
.hs_solve_hyper <- function(hv, lookup, piv, pv, idx, val, m = Inf) {
  mark <- logical(length(piv))
  list_index <- integer()
  count_pivot <- 0; count_entry <- 0
  for (r0 in hv$index) {
    it <- lookup[r0]
    if (mark[it]) next
    hi <- it; hk <- 1L; stack <- integer()
    mark[hi] <- TRUE
    repeat {
      if (hk <= length(idx[[hi]])) {
        child <- lookup[idx[[hi]][hk]]
        hk <- hk + 1L
        if (!mark[child]) {
          mark[child] <- TRUE
          stack <- c(stack, hi, hk)
          hi <- child; hk <- 1L
          if (hi > m) { count_pivot <- count_pivot + 1; count_entry <- count_entry + length(idx[[hi]]) }
        }
      } else {
        list_index <- c(list_index, hi)
        ns <- length(stack)
        if (!ns) break
        hk <- stack[ns]; hi <- stack[ns - 1L]
        stack <- stack[seq_len(ns - 2L)]
      }
    }
  }
  hv$tick <- (hv$tick %||% 0) + count_pivot * 20 + count_entry * 10
  .hs_solve_sparse(hv, rev(list_index), piv, pv, idx, val)
}

.hs_use_sparse <- function(hv, m, expected, limit) {
  hv$count < 0L || hv$count / m > 0.05 || expected > limit
}

.hs_ftran_l <- function(fa, hv, expected) {
  m <- fa$num_row
  if (.hs_use_sparse(hv, m, expected, 0.15))
    .hs_solve_sparse(hv, seq_len(m), fa$l_piv, NULL, fa$l_idx, fa$l_val)
  else .hs_solve_hyper(hv, fa$l_lookup, fa$l_piv, NULL, fa$l_idx, fa$l_val)
}

.hs_btran_l <- function(fa, hv, expected) {
  m <- fa$num_row
  if (.hs_use_sparse(hv, m, expected, 0.10))
    .hs_solve_sparse(hv, rev(seq_len(m)), fa$l_piv, NULL, fa$lr_idx, fa$lr_val)
  else .hs_solve_hyper(hv, fa$l_lookup, fa$l_piv, NULL, fa$lr_idx, fa$lr_val)
}

# ftranU with ftranFT before it; the packed vector is the spike of the FT
# update
.hs_ftran_u <- function(fa, hv, expected, pack) {
  a <- hv$array; index <- hv$index
  npf <- length(fa$pf_piv); pfx <- sum(lengths(fa$pf_idx))
  tick <- (hv$tick %||% 0) + npf * 20 + pfx * 5
  if (pfx %/% (npf + 1) < 5) tick <- tick + pfx * 5
  for (p in seq_along(fa$pf_piv)) {
    r <- fa$pf_piv[p]
    v0 <- a[r]; v1 <- v0
    pi <- fa$pf_idx[[p]]; pvv <- fa$pf_val[[p]]
    for (k in seq_along(pi)) v1 <- v1 - a[pi[k]] * pvv[k]
    if (v0 != 0 || v1 != 0) {
      if (v0 == 0) index <- c(index, r)
      a[r] <- if (abs(v1) < .hs_tiny) .hs_zero else v1
    }
  }
  hv <- .hs_hv_tight(list(array = a, index = index, count = length(index), tick = tick))
  if (pack) fa$pack_aq <- list(index = hv$index, value = hv$array[hv$index])
  m <- fa$num_row
  if (.hs_use_sparse(hv, m, expected, 0.10)) {
    hv <- .hs_solve_sparse(hv, rev(seq_along(fa$u_piv)), fa$u_piv, fa$u_pv, fa$u_idx, fa$u_val, m)
    hv$tick <- hv$tick + hv$ft * 15 + (length(fa$u_piv) - m) * 10
    hv
  } else .hs_solve_hyper(hv, fa$u_lookup, fa$u_piv, fa$u_pv, fa$u_idx, fa$u_val, m)
}

# btranU with btranFT after it; the packed vector is the row eta of the FT
# update
.hs_btran_u <- function(fa, hv, expected, pack) {
  m <- fa$num_row
  hv <- if (.hs_use_sparse(hv, m, expected, 0.15)) {
    hv <- .hs_solve_sparse(hv, seq_along(fa$u_piv), fa$u_piv, fa$u_pv, fa$ur_idx, fa$ur_val, m)
    hv$tick <- hv$tick + hv$ft * 15 + (length(fa$u_piv) - m) * 10
    hv
  } else .hs_solve_hyper(hv, fa$u_lookup, fa$u_piv, fa$u_pv, fa$ur_idx, fa$ur_val, m)
  hv <- .hs_hv_tight(hv)
  if (pack) fa$pack_ep <- list(index = hv$index, value = hv$array[hv$index])
  a <- hv$array; index <- hv$index
  rt <- 0
  for (p in rev(seq_along(fa$pf_piv))) {
    x <- a[fa$pf_piv[p]]
    if (x != 0) {
      rows <- fa$pf_idx[[p]]
      rt <- rt + length(rows)
      v0 <- a[rows]
      index <- c(index, rows[v0 == 0])
      v1 <- v0 - x * fa$pf_val[[p]]
      a[rows] <- ifelse(abs(v1) < .hs_tiny, .hs_zero, v1)
    }
  }
  .hs_hv_tight(list(array = a, index = index, count = length(index),
                    tick = hv$tick + rt * 15 + length(fa$pf_piv) * 10))
}

# HFactor::ftranCall and btranCall
.hs_factor_ftran <- function(fa, hv, expected, pack = FALSE) {
  use_indices <- hv$count >= 0L
  hv <- .hs_ftran_l(fa, hv, expected)
  hv <- .hs_ftran_u(fa, hv, expected, pack)
  if (use_indices) hv <- .hs_hv_reindex(hv)
  hv
}

.hs_factor_btran <- function(fa, hv, expected, pack = FALSE) {
  use_indices <- hv$count >= 0L
  hv <- .hs_btran_u(fa, hv, expected, pack)
  hv <- .hs_btran_l(fa, hv, expected)
  if (use_indices) hv <- .hs_hv_reindex(hv)
  hv
}

# HFactor::updateFT: column aq (its packed spike and its value alpha in
# row r) replaces the basic variable of row r; ep is the packed row eta
.hs_factor_update_ft <- function(fa, alpha, r) {
  p <- fa$u_lookup[r]
  pivot <- fa$u_pv[p]
  fa$u_piv[p] <- 0L
  for (t in fa$ur_idx[[p]]) {
    q <- fa$u_lookup[t]
    ui <- fa$u_idx[[q]]; n <- length(ui)
    k <- match(r, ui)
    fa$u_idx[[q]][k] <- ui[n]; fa$u_val[[q]][k] <- fa$u_val[[q]][n]
    fa$u_idx[[q]] <- fa$u_idx[[q]][-n]; fa$u_val[[q]] <- fa$u_val[[q]][-n]
  }
  for (t in fa$u_idx[[p]]) {
    q <- fa$u_lookup[t]
    ri <- fa$ur_idx[[q]]; n <- length(ri)
    k <- match(r, ri)
    fa$ur_idx[[q]][k] <- ri[n]; fa$ur_val[[q]][k] <- fa$ur_val[[q]][n]
    fa$ur_idx[[q]] <- fa$ur_idx[[q]][-n]; fa$ur_val[[q]] <- fa$ur_val[[q]][-n]
  }
  aq <- fa$pack_aq
  keep <- aq$index != r
  ni <- aq$index[keep]; nv <- aq$value[keep]
  np <- length(fa$u_piv) + 1L
  for (k in seq_along(ni)) {
    q <- fa$u_lookup[ni[k]]
    fa$ur_idx[[q]] <- c(fa$ur_idx[[q]], r); fa$ur_val[[q]] <- c(fa$ur_val[[q]], nv[k])
  }
  fa$u_idx[[np]] <- ni; fa$u_val[[np]] <- nv
  fa$ur_idx[np] <- list(integer()); fa$ur_val[np] <- list(numeric())
  fa$u_lookup[r] <- np
  fa$u_piv[np] <- r
  fa$u_pv[np] <- pivot * alpha
  ep <- fa$pack_ep
  keep <- ep$index != r
  q <- length(fa$pf_piv) + 1L
  fa$pf_idx[q] <- list(ep$index[keep]); fa$pf_val[q] <- list(-ep$value[keep] * pivot)
  fa$pf_piv[q] <- r
}
