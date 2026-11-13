# BASICLU, the sparse LU factorization with Forrest-Tomlin updates that IPX
# (R/optimization-scipy-highs-ipx.R) uses for its basis matrices, ported from
# HiGHS 1.12.0 (highs/ipm/basiclu and the BasicLu wrapper of
# highs/ipm/ipx/basiclu_wrapper.cc).
#
# Copyright (C) 2016-2019 ERGO-Code; HiGHS is released under the MIT
# License (Copyright (c) 2024 HiGHS).
#
# The object keeps the arrays of struct lu under their C names. Indices are
# the 0-based C indices; an array element a[k] of the C code is a[k + 1L]
# here. The memory sizes Lmem, Umem and Wmem are taken large enough that the
# reallocation requests of the C code never occur; garbage collection keeps
# the order of the lines and of the entries within each line, so the factors
# do not depend on the memory sizes. Multiply adds (a -= b * c, a += b * c)
# use .ipx_madd(), which rounds the product as SciPy's Linux x86-64 build
# of HiGHS does (no FMA contraction). The
# condition estimates of lu_condest(), which only feed debugging output, are
# not computed.

.blu_bigmem <- 1e9

.blu_list_init <- function(L, fl, bl, nelem, nlist) {
  idx <- seq_len(nelem + nlist) - 1L
  L[[fl]] <- idx
  L[[bl]] <- idx
  max(1L, nlist)
}

.blu_list_add <- function(L, fl, bl, elem, list, nelem, minname = NULL) {
  f <- L[[fl]]; b <- L[[bl]]
  temp <- b[nelem + list + 1L]
  b[nelem + list + 1L] <- elem
  b[elem + 1L] <- temp
  f[temp + 1L] <- elem
  f[elem + 1L] <- nelem + list
  L[[fl]] <- f; L[[bl]] <- b
  if (!is.null(minname) && list > 0L && list < L[[minname]]) L[[minname]] <- list
}

.blu_list_remove <- function(L, fl, bl, elem) {
  f <- L[[fl]]; b <- L[[bl]]
  f[b[elem + 1L] + 1L] <- f[elem + 1L]
  b[f[elem + 1L] + 1L] <- b[elem + 1L]
  f[elem + 1L] <- elem
  b[elem + 1L] <- elem
  L[[fl]] <- f; L[[bl]] <- b
}

.blu_list_move <- function(L, fl, bl, elem, list, nelem, minname = NULL) {
  .blu_list_remove(L, fl, bl, elem)
  .blu_list_add(L, fl, bl, elem, list, nelem, minname)
}

.blu_list_swap <- function(L, fl, bl, e1, e2) {
  f <- L[[fl]]; b <- L[[bl]]
  e1next <- f[e1 + 1L]; e2next <- f[e2 + 1L]
  e1prev <- b[e1 + 1L]; e2prev <- b[e2 + 1L]
  if (e1next == e2) {
    f[e2 + 1L] <- e1; b[e1 + 1L] <- e2
    f[e1prev + 1L] <- e2; b[e2 + 1L] <- e1prev
    f[e1 + 1L] <- e2next; b[e2next + 1L] <- e1
  } else if (e2next == e1) {
    f[e1 + 1L] <- e2; b[e2 + 1L] <- e1
    f[e2 + 1L] <- e1next; b[e1next + 1L] <- e2
    f[e2prev + 1L] <- e1; b[e1 + 1L] <- e2prev
  } else {
    f[e2 + 1L] <- e1next; b[e1next + 1L] <- e2
    f[e2prev + 1L] <- e1; b[e1 + 1L] <- e2prev
    f[e1prev + 1L] <- e2; b[e2 + 1L] <- e1prev
    f[e1 + 1L] <- e2next; b[e2next + 1L] <- e1
  }
  L[[fl]] <- f; L[[bl]] <- b
}

# lu_file_empty for the file of nlines lines in Wbegin, Wend, Wflink, Wblink
.blu_file_empty <- function(L, nlines, fmem) {
  L$Wbegin[nlines + 1L] <- 0L
  L$Wend[nlines + 1L] <- fmem
  L$Wbegin[seq_len(nlines)] <- 0L
  L$Wend[seq_len(nlines)] <- 0L
  L$Wflink[seq_len(nlines)] <- seq_len(nlines)
  L$Wblink[seq_len(nlines) + 1L] <- seq_len(nlines) - 1L
  L$Wflink[nlines + 1L] <- 0L
  L$Wblink[1L] <- nlines
}

.blu_file_reappend <- function(L, line, nlines, extra_space) {
  used <- L$Wbegin[nlines + 1L]
  ibeg <- L$Wbegin[line + 1L]; iend <- L$Wend[line + 1L]
  L$Wbegin[line + 1L] <- used
  if (iend > ibeg) {
    k <- ibeg:(iend - 1L)
    L$Windex[used + seq_along(k)] <- L$Windex[k + 1L]
    L$Wvalue[used + seq_along(k)] <- L$Wvalue[k + 1L]
    used <- used + length(k)
  }
  L$Wend[line + 1L] <- used
  L$Wbegin[nlines + 1L] <- used + extra_space
  .blu_list_move(L, "Wflink", "Wblink", line, 0L, nlines)
}

.blu_file_compress <- function(L, nlines, stretch, pad) {
  used <- 0L; extra_space <- 0L; nz <- 0L
  i <- L$Wflink[nlines + 1L]
  while (i < nlines) {
    ibeg <- L$Wbegin[i + 1L]; iend <- L$Wend[i + 1L]
    used <- used + extra_space
    if (used > ibeg) used <- ibeg
    L$Wbegin[i + 1L] <- used
    if (iend > ibeg) {
      k <- ibeg:(iend - 1L)
      vi <- L$Windex[k + 1L]; vx <- L$Wvalue[k + 1L]
      L$Windex[used + seq_along(k)] <- vi
      L$Wvalue[used + seq_along(k)] <- vx
      used <- used + length(k)
    }
    L$Wend[i + 1L] <- used
    extra_space <- as.integer(stretch * (iend - ibeg) + pad)
    nz <- nz + iend - ibeg
    i <- L$Wflink[i + 1L]
  }
  used <- used + extra_space
  if (used > L$Wbegin[nlines + 1L]) used <- L$Wbegin[nlines + 1L]
  L$Wbegin[nlines + 1L] <- used
  nz
}

# basiclu_initialize with the BasicLu wrapper defaults
.blu_new <- function(m) {
  L <- new.env(parent = emptyenv())
  L$m <- m
  L$Lmem <- .blu_bigmem; L$Umem <- .blu_bigmem; L$Wmem <- .blu_bigmem
  L$droptol <- 1e-20; L$abstol <- 1e-14; L$reltol <- 0.1
  L$nzbias <- 1L; L$maxsearch <- 3L; L$pad <- 4L; L$stretch <- 0.3
  L$compress_thres <- 0.5; L$sparse_thres <- 0.05; L$search_rows <- 1L
  L$nfactorize <- 0L; L$nupdate_total <- 0L; L$nforrest_total <- 0L
  L$nsymperm_total <- 0L
  L$iwork0 <- integer(m); L$work0 <- numeric(m); L$work1 <- numeric(m)
  L$col_pivot <- numeric(m); L$row_pivot <- numeric(m)
  L$Lindex <- integer(0); L$Lvalue <- numeric(0)
  L$Uindex <- integer(0); L$Uvalue <- numeric(0)
  L$Windex <- integer(0); L$Wvalue <- numeric(0)
  L$Wbegin <- integer(2L * m + 1L); L$Wend <- integer(2L * m + 1L)
  L$Wflink <- integer(2L * m + 1L); L$Wblink <- integer(2L * m + 1L)
  .blu_reset(L)
  L
}

.blu_reset <- function(L) {
  m <- L$m
  L$nupdate <- -1L; L$nforrest <- 0L
  L$Lnz <- 0L; L$Unz <- 0L; L$Rnz <- 0L
  L$min_pivot <- 0; L$max_pivot <- 0; L$max_eta <- 0
  L$update_cost_numer <- 0; L$update_cost_denom <- 1
  L$onenorm <- 0; L$infnorm <- 0; L$residual_test <- 0
  L$matrix_nz <- 0L; L$rank <- 0L; L$bump_size <- 0L; L$bump_nz <- 0L
  L$nsearch_pivot <- 0L; L$nexpand <- 0L; L$ngarbage <- 0L
  L$factor_flops <- 0; L$pivot_error <- 0
  L$pivot_row <- -1L; L$pivot_col <- -1L
  L$ftran_for_update <- -1L; L$btran_for_update <- -1L
  L$marker <- 0L; L$pivotlen <- 0L; L$rankdef <- 0L
  L$min_colnz <- 1L; L$min_rownz <- 1L
  L$Wend[2L * m + 1L] <- L$Wmem
  L$iwork0 <- integer(m)
  L$work0 <- numeric(m)
}

# basiclu_factorize: B is given by columns, Bcols[[j]] = list(i = 0-based row
# indices, x = values)
.blu_factorize <- function(L, Bcols) {
  .blu_reset(L)
  .blu_singletons(L, Bcols)
  .blu_setup_bump(L, Bcols)
  .blu_factorize_bump(L)
  .blu_build_factors(L)
  L$nupdate <- 0L
  L$ftran_for_update <- -1L
  L$btran_for_update <- -1L
  L$nfactorize <- L$nfactorize + 1L
  .blu_residual_test(L, Bcols)
  factor_cost <- .ipx_madd(0.008, L$factor_flops, .ipx_madd(0.20, L$nsearch_pivot,
    .ipx_madd(0.20, L$bump_nz, .ipx_madd(0.07, L$matrix_nz, 0.04 * L$m))))
  L$update_cost_denom <- factor_cost * 250
  if (L$rank < L$m) "singular" else "ok"
}

.blu_singletons <- function(L, Bcols) {
  m <- L$m
  Bnz <- sum(vapply(Bcols, function(cl) length(cl$i), 1L))
  # B rowwise in Btp, Bti, Btx
  cnt <- integer(m)
  for (j in seq_len(m)) for (i in Bcols[[j]]$i) cnt[i + 1L] <- cnt[i + 1L] + 1L
  Btp <- c(0L, cumsum(cnt))
  nxt <- Btp[seq_len(m)]
  Bti <- integer(Bnz); Btx <- numeric(Bnz)
  for (j in seq_len(m)) {
    cl <- Bcols[[j]]
    for (t in seq_along(cl$i)) {
      i <- cl$i[t]
      put <- nxt[i + 1L]; nxt[i + 1L] <- put + 1L
      Bti[put + 1L] <- j - 1L
      Btx[put + 1L] <- cl$x[t]
    }
  }
  L$pinv <- rep(-1L, m); L$qinv <- rep(-1L, m)
  L$Lbegin_p <- integer(m + 1L); L$Ubegin <- integer(m + 1L)
  L$Lindex <- integer(0); L$Lvalue <- numeric(0)
  L$Uindex <- integer(0); L$Uvalue <- numeric(0)
  rank <- 0L
  rank <- .blu_singleton_cols(L, Bcols, Btp, Bti, Btx, rank)
  rank <- .blu_singleton_rows(L, Bcols, Btp, Bti, Btx, rank)
  L$pinv[L$pinv < 0L] <- -1L
  L$qinv[L$qinv < 0L] <- -1L
  L$matrix_nz <- Bnz
  L$rank <- rank
}

.blu_singleton_cols <- function(L, Bcols, Btp, Bti, Btx, rank) {
  m <- L$m; abstol <- L$abstol
  pinv <- L$pinv; qinv <- L$qinv
  iset <- integer(m); queue <- integer(m)
  rk <- rank
  tail <- 0L
  for (j in 0:(m - 1L)) {
    if (qinv[j + 1L] < 0L) {
      bi <- Bcols[[j + 1L]]$i
      nz <- length(bi)
      i <- 0L
      for (r in bi) i <- bitwXor(i, r)
      iset[j + 1L] <- i
      qinv[j + 1L] <- -nz - 1L
      if (nz == 1L) { queue[tail + 1L] <- j; tail <- tail + 1L }
    }
  }
  put <- L$Ubegin[rank + 1L]
  front <- 0L
  while (front < tail) {
    j <- queue[front + 1L]
    front <- front + 1L
    if (qinv[j + 1L] == -1L) next
    i <- iset[j + 1L]
    pos <- Btp[i + 1L]
    while (Bti[pos + 1L] != j) pos <- pos + 1L
    piv <- Btx[pos + 1L]
    if (piv == 0 || abs(piv) < abstol) next
    qinv[j + 1L] <- rank
    pinv[i + 1L] <- rank
    for (pos in Btp[i + 1L]:(Btp[i + 2L] - 1L)) {
      j2 <- Bti[pos + 1L]
      if (qinv[j2 + 1L] < 0L) {
        L$Uindex[put + 1L] <- j2
        L$Uvalue[put + 1L] <- Btx[pos + 1L]
        put <- put + 1L
        iset[j2 + 1L] <- bitwXor(iset[j2 + 1L], i)
        qinv[j2 + 1L] <- qinv[j2 + 1L] + 1L
        if (qinv[j2 + 1L] == -2L) { queue[tail + 1L] <- j2; tail <- tail + 1L }
      }
    }
    L$Ubegin[rank + 2L] <- put
    L$col_pivot[j + 1L] <- piv
    rank <- rank + 1L
  }
  pos <- L$Lbegin_p[rk + 1L]
  while (rk < rank) {
    L$Lindex[pos + 1L] <- -1L
    pos <- pos + 1L
    L$Lbegin_p[rk + 2L] <- pos
    rk <- rk + 1L
  }
  L$pinv <- pinv; L$qinv <- qinv
  rank
}

.blu_singleton_rows <- function(L, Bcols, Btp, Bti, Btx, rank) {
  m <- L$m; abstol <- L$abstol
  pinv <- L$pinv; qinv <- L$qinv
  iset <- integer(m); queue <- integer(m)
  rk <- rank
  tail <- 0L
  for (i in 0:(m - 1L)) {
    if (pinv[i + 1L] < 0L) {
      nz <- Btp[i + 2L] - Btp[i + 1L]
      j <- 0L
      if (nz > 0L) for (pos in Btp[i + 1L]:(Btp[i + 2L] - 1L)) j <- bitwXor(j, Bti[pos + 1L])
      iset[i + 1L] <- j
      pinv[i + 1L] <- -nz - 1L
      if (nz == 1L) { queue[tail + 1L] <- i; tail <- tail + 1L }
    }
  }
  put <- L$Lbegin_p[rank + 1L]
  front <- 0L
  while (front < tail) {
    i <- queue[front + 1L]
    front <- front + 1L
    if (pinv[i + 1L] == -1L) next
    j <- iset[i + 1L]
    cl <- Bcols[[j + 1L]]
    t <- match(i, cl$i)
    piv <- cl$x[t]
    if (piv == 0 || abs(piv) < abstol) next
    qinv[j + 1L] <- rank
    pinv[i + 1L] <- rank
    for (t in seq_along(cl$i)) {
      i2 <- cl$i[t]
      if (pinv[i2 + 1L] < 0L) {
        L$Lindex[put + 1L] <- i2
        L$Lvalue[put + 1L] <- cl$x[t] / piv
        put <- put + 1L
        iset[i2 + 1L] <- bitwXor(iset[i2 + 1L], j)
        pinv[i2 + 1L] <- pinv[i2 + 1L] + 1L
        if (pinv[i2 + 1L] == -2L) { queue[tail + 1L] <- i2; tail <- tail + 1L }
      }
    }
    L$Lindex[put + 1L] <- -1L
    put <- put + 1L
    L$Lbegin_p[rank + 2L] <- put
    L$col_pivot[j + 1L] <- piv
    rank <- rank + 1L
  }
  pos <- L$Ubegin[rk + 1L]
  while (rk < rank) {
    L$Ubegin[rk + 2L] <- pos
    rk <- rk + 1L
  }
  L$pinv <- pinv; L$qinv <- qinv
  rank
}

.blu_setup_bump <- function(L, Bcols) {
  m <- L$m; rank <- L$rank
  Lnz <- L$Lbegin_p[rank + 1L] - rank
  Unz <- L$Ubegin[rank + 1L]
  abstol <- L$abstol; pad <- L$pad; stretch <- L$stretch
  pinv <- L$pinv; qinv <- L$qinv
  bump_nz <- L$matrix_nz - Lnz - Unz - rank
  L$Windex <- integer(0); L$Wvalue <- numeric(0)
  .blu_file_empty(L, 2L * m, L$Wmem)
  L$colcount_flink <- integer(2L * m + 2L); L$colcount_blink <- integer(2L * m + 2L)
  L$rowcount_flink <- integer(2L * m + 2L); L$rowcount_blink <- integer(2L * m + 2L)
  L$min_colnz <- .blu_list_init(L, "colcount_flink", "colcount_blink", m, m + 2L)
  iwork0 <- integer(m)
  put <- 0L
  for (j in 0:(m - 1L)) {
    if (qinv[j + 1L] >= 0L) next
    cl <- Bcols[[j + 1L]]
    act <- pinv[cl$i + 1L] < 0L
    cnz <- sum(act)
    cmx <- 0
    for (x in cl$x[act]) cmx <- max(cmx, abs(x))
    if (cmx == 0 || cmx < abstol) {
      L$col_pivot[j + 1L] <- 0
      .blu_list_add(L, "colcount_flink", "colcount_blink", j, 0L, m, "min_colnz")
      bump_nz <- bump_nz - cnz
    } else {
      L$col_pivot[j + 1L] <- cmx
      .blu_list_add(L, "colcount_flink", "colcount_blink", j, cnz, m, "min_colnz")
      L$Wbegin[j + 1L] <- put
      ii <- cl$i[act]
      if (cnz > 0L) {
        L$Windex[put + seq_len(cnz)] <- ii
        L$Wvalue[put + seq_len(cnz)] <- cl$x[act]
        iwork0[ii + 1L] <- iwork0[ii + 1L] + 1L
      }
      put <- put + cnz
      L$Wend[j + 1L] <- put
      put <- as.integer(put + stretch * cnz + pad)
      .blu_list_move(L, "Wflink", "Wblink", j, 0L, 2L * m)
    }
  }
  L$min_rownz <- .blu_list_init(L, "rowcount_flink", "rowcount_blink", m, m + 2L)
  for (i in 0:(m - 1L)) {
    if (pinv[i + 1L] >= 0L) next
    rnz <- iwork0[i + 1L]
    iwork0[i + 1L] <- 0L
    .blu_list_add(L, "rowcount_flink", "rowcount_blink", i, rnz, m, "min_rownz")
    L$Wbegin[m + i + 1L] <- put
    L$Wend[m + i + 1L] <- put
    put <- put + rnz
    .blu_list_move(L, "Wflink", "Wblink", m + i, 0L, 2L * m)
    put <- as.integer(put + stretch * rnz + pad)
  }
  for (j in 0:(m - 1L)) {
    b <- L$Wbegin[j + 1L]; e <- L$Wend[j + 1L]
    if (e > b) for (pos in b:(e - 1L)) {
      i <- L$Windex[pos + 1L]
      w <- L$Wend[m + i + 1L]
      L$Windex[w + 1L] <- j
      L$Wend[m + i + 1L] <- w + 1L
    }
  }
  L$Wbegin[2L * m + 1L] <- put
  L$bump_nz <- bump_nz
  L$bump_size <- m - rank
  L$iwork0 <- iwork0
}

.blu_factorize_bump <- function(L) {
  m <- L$m
  while (L$rank + L$rankdef < m) {
    if (L$pivot_col < 0L) .blu_markowitz(L)
    if (L$pivot_row < 0L) {
      .blu_list_remove(L, "colcount_flink", "colcount_blink", L$pivot_col)
      L$pivot_col <- -1L
      L$rankdef <- L$rankdef + 1L
    } else {
      .blu_pivot(L)
      L$pinv[L$pivot_row + 1L] <- L$rank
      L$qinv[L$pivot_col + 1L] <- L$rank
      L$pivot_col <- -1L
      L$pivot_row <- -1L
      L$rank <- L$rank + 1L
    }
  }
}

.blu_markowitz <- function(L) {
  m <- L$m
  Wbegin <- L$Wbegin; Wend <- L$Wend; Windex <- L$Windex; Wvalue <- L$Wvalue
  colmax <- L$col_pivot
  abstol <- L$abstol; reltol <- min(L$reltol, 1); maxsearch <- L$maxsearch
  search_rows <- L$search_rows
  nz_start <- if (search_rows) min(L$min_colnz, L$min_rownz) else L$min_colnz
  pivot_row <- -1L; pivot_col <- -1L
  MC <- as.numeric(m) * m
  nsearch <- 0L; min_colnz <- -1L; min_rownz <- -1L
  ccf <- L$colcount_flink
  done <- FALSE
  if (ccf[m + 1L] != m) {
    pivot_col <- ccf[m + 1L]
    done <- TRUE
  }
  nz <- nz_start
  while (!done && nz <= m) {
    j <- ccf[m + nz + 1L]
    while (!done && j < m) {
      if (min_colnz == -1L) min_colnz <- nz
      cmx <- colmax[j + 1L]
      if (!(cmx == 0 || cmx < abstol)) {
        tol <- max(abstol, reltol * cmx)
        b <- Wbegin[j + 1L]; e <- Wend[j + 1L]
        pos <- b
        while (pos < e) {
          x <- abs(Wvalue[pos + 1L])
          if (!(x == 0 || x < tol)) {
            i <- Windex[pos + 1L]
            nz1 <- nz
            nz2 <- Wend[m + i + 1L] - Wbegin[m + i + 1L]
            mc <- (nz1 - 1) * (nz2 - 1)
            if (mc < MC) {
              MC <- mc
              pivot_row <- i
              pivot_col <- j
              if (search_rows && MC <= (nz1 - 1) * (nz1 - 1)) { done <- TRUE; break }
            }
          }
          pos <- pos + 1L
        }
        if (done) break
        nsearch <- nsearch + 1L
        if (nsearch >= maxsearch) { done <- TRUE; break }
      }
      j <- ccf[j + 1L]
    }
    if (done) break
    if (!search_rows) { nz <- nz + 1L; next }
    i <- L$rowcount_flink[m + nz + 1L]
    while (i < m) {
      if (min_rownz == -1L) min_rownz <- nz
      inext <- L$rowcount_flink[i + 1L]
      cheap <- FALSE; found <- FALSE
      b <- Wbegin[m + i + 1L]; e <- Wend[m + i + 1L]
      pos <- b
      while (pos < e) {
        j <- Windex[pos + 1L]
        nz1 <- nz
        nz2 <- Wend[j + 1L] - Wbegin[j + 1L]
        mc <- (nz1 - 1) * (nz2 - 1)
        pos <- pos + 1L
        if (mc >= MC) next
        cheap <- TRUE
        cmx <- colmax[j + 1L]
        if (cmx == 0 || cmx < abstol) next
        where <- Wbegin[j + 1L]
        while (Windex[where + 1L] != i) where <- where + 1L
        x <- abs(Wvalue[where + 1L])
        if (x >= abstol && x >= reltol * cmx) {
          found <- TRUE
          MC <- mc
          pivot_row <- i
          pivot_col <- j
          if (MC <= nz1 * (nz1 - 1)) { done <- TRUE; break }
        }
      }
      if (done) break
      if (cheap && !found) {
        .blu_list_move(L, "rowcount_flink", "rowcount_blink", i, m + 1L, m)
      } else {
        nsearch <- nsearch + 1L
        if (nsearch >= maxsearch) { done <- TRUE; break }
      }
      i <- inext
    }
    nz <- nz + 1L
  }
  L$pivot_row <- pivot_row
  L$pivot_col <- pivot_col
  L$nsearch_pivot <- L$nsearch_pivot + nsearch
  if (min_colnz >= 0L) L$min_colnz <- min_colnz
  if (min_rownz >= 0L) L$min_rownz <- min_rownz
}

.blu_pivot <- function(L) {
  m <- L$m; rank <- L$rank
  pc <- L$pivot_col; pr <- L$pivot_row
  nz_col <- L$Wend[pc + 1L] - L$Wbegin[pc + 1L]
  nz_row <- L$Wend[m + pr + 1L] - L$Wbegin[m + pr + 1L]
  if (nz_row == 1L) .blu_pivot_singleton_row(L)
  else if (nz_col == 1L) .blu_pivot_singleton_col(L)
  else if (nz_col == 2L) .blu_pivot_doubleton_col(L)
  else .blu_pivot_general(L, small = nz_col - 1L <= 64L)
  b <- L$Ubegin[rank + 1L]; e <- L$Ubegin[rank + 2L]
  if (e > b) for (pos in b:(e - 1L)) {
    j <- L$Uindex[pos + 1L]
    cm <- L$col_pivot[j + 1L]
    if (cm == 0 || cm < L$abstol) .blu_remove_col(L, j)
  }
  L$factor_flops <- L$factor_flops + (nz_col - 1) * (nz_row - 1)
}

.blu_pivot_cleanup <- function(L, pivot, cbeg, rbeg) {
  m <- L$m
  L$col_pivot[L$pivot_col + 1L] <- pivot
  L$Wend[L$pivot_col + 1L] <- cbeg
  L$Wend[m + L$pivot_row + 1L] <- rbeg
  .blu_list_remove(L, "colcount_flink", "colcount_blink", L$pivot_col)
  .blu_list_remove(L, "rowcount_flink", "rowcount_blink", L$pivot_row)
}

# lu_pivot_any and lu_pivot_small: the latter drops cancelled updates
.blu_pivot_general <- function(L, small) {
  m <- L$m; rank <- L$rank
  droptol <- L$droptol; pad <- L$pad; stretch <- L$stretch
  pc <- L$pivot_col; pr <- L$pivot_row
  cbeg <- L$Wbegin[pc + 1L]; cend <- L$Wend[pc + 1L]
  rbeg <- L$Wbegin[m + pr + 1L]; rend <- L$Wend[m + pr + 1L]
  cnz1 <- cend - cbeg - 1L; rnz1 <- rend - rbeg - 1L
  grow <- 0; where <- -1L
  for (pos in cbeg:(cend - 1L)) {
    i <- L$Windex[pos + 1L]
    if (i == pr) where <- pos
    else {
      nz <- L$Wend[m + i + 1L] - L$Wbegin[m + i + 1L]
      grow <- as.integer(grow + nz + rnz1 + stretch * (nz + rnz1) + pad)
    }
  }
  t <- L$Windex[cbeg + 1L]; L$Windex[cbeg + 1L] <- L$Windex[where + 1L]; L$Windex[where + 1L] <- t
  t <- L$Wvalue[cbeg + 1L]; L$Wvalue[cbeg + 1L] <- L$Wvalue[where + 1L]; L$Wvalue[where + 1L] <- t
  pivot <- L$Wvalue[cbeg + 1L]
  where <- -1L
  for (rpos in rbeg:(rend - 1L)) {
    j <- L$Windex[rpos + 1L]
    if (j == pc) where <- rpos
    else {
      nz <- L$Wend[j + 1L] - L$Wbegin[j + 1L]
      grow <- as.integer(grow + nz + cnz1 + stretch * (nz + cnz1) + pad)
    }
  }
  t <- L$Windex[rbeg + 1L]; L$Windex[rbeg + 1L] <- L$Windex[where + 1L]; L$Windex[where + 1L] <- t
  room <- L$Wend[2L * m + 1L] - L$Wbegin[2L * m + 1L]
  if (grow > room) {
    .blu_file_compress(L, 2L * m, stretch, pad)
    cbeg <- L$Wbegin[pc + 1L]; cend <- L$Wend[pc + 1L]
    rbeg <- L$Wbegin[m + pr + 1L]; rend <- L$Wend[m + pr + 1L]
    L$ngarbage <- L$ngarbage + 1L
  }
  Uput <- L$Ubegin[rank + 1L]
  marked <- L$iwork0
  work <- L$work0
  pc_rows <- L$Windex[cbeg + seq_len(cnz1) + 1L]
  marked[pc_rows + 1L] <- seq_len(cnz1)
  wx <- L$Wvalue[cbeg + seq_len(cnz1) + 1L]
  prow_cols <- L$Windex[rbeg + seq_len(rnz1) + 1L]
  cancelled <- matrix(FALSE, max(rnz1, 1L), max(cnz1, 1L))
  for (cn in seq_len(rnz1)) {
    j <- prow_cols[cn]
    cmx <- 0
    where <- -1L
    pos1 <- L$Wbegin[j + 1L]
    put <- pos1
    e <- L$Wend[j + 1L]
    for (pos in pos1:(e - 1L)) {
      i <- L$Windex[pos + 1L]
      position <- marked[i + 1L]
      if (position > 0L) {
        work[position + 1L] <- L$Wvalue[pos + 1L]
      } else {
        if (i == pr) where <- put
        else {
          x <- abs(L$Wvalue[pos + 1L])
          if (x > cmx) cmx <- x
        }
        L$Windex[put + 1L] <- L$Windex[pos + 1L]
        L$Wvalue[put + 1L] <- L$Wvalue[pos + 1L]
        put <- put + 1L
      }
    }
    L$Wend[j + 1L] <- put
    t <- L$Windex[pos1 + 1L]; L$Windex[pos1 + 1L] <- L$Windex[where + 1L]; L$Windex[where + 1L] <- t
    t <- L$Wvalue[pos1 + 1L]; L$Wvalue[pos1 + 1L] <- L$Wvalue[where + 1L]; L$Wvalue[where + 1L] <- t
    xrj <- L$Wvalue[pos1 + 1L]
    room <- L$Wbegin[L$Wflink[j + 1L] + 1L] - put
    if (room < cnz1) {
      nz <- L$Wend[j + 1L] - L$Wbegin[j + 1L]
      room <- as.integer(cnz1 + stretch * (nz + cnz1) + pad)
      .blu_file_reappend(L, j, 2L * m, room)
      put <- L$Wend[j + 1L]
      L$nexpand <- L$nexpand + 1L
    }
    a <- xrj / pivot
    w <- .ipx_madd(-a, wx, work[seq_len(cnz1) + 1L])
    work[seq_len(cnz1) + 1L] <- 0
    if (small) {
      keep <- abs(w) > droptol
      cancelled[cn, ] <- !keep
    } else {
      keep <- rep(TRUE, cnz1)
    }
    nk <- sum(keep)
    if (nk) {
      L$Windex[put + seq_len(nk)] <- pc_rows[keep]
      L$Wvalue[put + seq_len(nk)] <- w[keep]
      cmx <- max(cmx, abs(w[keep]))
      put <- put + nk
    }
    L$Wend[j + 1L] <- put
    if (abs(xrj) > droptol) {
      L$Uindex[Uput + 1L] <- j
      L$Uvalue[Uput + 1L] <- xrj
      Uput <- Uput + 1L
    }
    L$Wbegin[j + 1L] <- L$Wbegin[j + 1L] + 1L
    nz <- L$Wend[j + 1L] - L$Wbegin[j + 1L]
    .blu_list_move(L, "colcount_flink", "colcount_blink", j, nz, m, "min_colnz")
    L$col_pivot[j + 1L] <- cmx
  }
  marked[pc_rows + 1L] <- 0L
  # row file update
  prow_all <- L$Windex[rbeg + seq_len(rnz1 + 1L)]
  marked[prow_all + 1L] <- 1L
  for (pn in seq_len(cnz1)) {
    i <- pc_rows[pn]
    b <- L$Wbegin[m + i + 1L]; e <- L$Wend[m + i + 1L]
    rowj <- L$Windex[b + seq_len(e - b)]
    keep <- rowj[marked[rowj + 1L] == 0L]
    put <- b + length(keep)
    if (length(keep)) L$Windex[b + seq_along(keep)] <- keep
    L$Wend[m + i + 1L] <- put
    room <- L$Wbegin[L$Wflink[m + i + 1L] + 1L] - put
    if (room < rnz1) {
      nz <- L$Wend[m + i + 1L] - L$Wbegin[m + i + 1L]
      room <- as.integer(rnz1 + stretch * (nz + rnz1) + pad)
      .blu_file_reappend(L, m + i, 2L * m, room)
      put <- L$Wend[m + i + 1L]
      L$nexpand <- L$nexpand + 1L
    }
    app <- if (small) prow_cols[!cancelled[seq_len(rnz1), pn]] else prow_cols
    if (length(app)) L$Windex[put + seq_along(app)] <- app
    put <- put + length(app)
    L$Wend[m + i + 1L] <- put
    nz <- put - L$Wbegin[m + i + 1L]
    .blu_list_move(L, "rowcount_flink", "rowcount_blink", i, nz, m, "min_rownz")
  }
  marked[prow_all + 1L] <- 0L
  L$iwork0 <- marked
  L$work0 <- work
  put <- L$Lbegin_p[rank + 1L]
  for (pn in seq_len(cnz1)) {
    x <- L$Wvalue[cbeg + pn + 1L] / pivot
    if (abs(x) > droptol) {
      L$Lindex[put + 1L] <- L$Windex[cbeg + pn + 1L]
      L$Lvalue[put + 1L] <- x
      put <- put + 1L
    }
  }
  L$Lindex[put + 1L] <- -1L
  put <- put + 1L
  L$Lbegin_p[rank + 2L] <- put
  L$Ubegin[rank + 2L] <- Uput
  .blu_pivot_cleanup(L, pivot, cbeg, rbeg)
}

.blu_pivot_singleton_row <- function(L) {
  m <- L$m; rank <- L$rank; droptol <- L$droptol
  pc <- L$pivot_col; pr <- L$pivot_row
  cbeg <- L$Wbegin[pc + 1L]; cend <- L$Wend[pc + 1L]
  rbeg <- L$Wbegin[m + pr + 1L]
  where <- cbeg
  while (L$Windex[where + 1L] != pr) where <- where + 1L
  pivot <- L$Wvalue[where + 1L]
  put <- L$Lbegin_p[rank + 1L]
  for (pos in cbeg:(cend - 1L)) {
    x <- L$Wvalue[pos + 1L] / pivot
    if (pos != where && abs(x) > droptol) {
      L$Lindex[put + 1L] <- L$Windex[pos + 1L]
      L$Lvalue[put + 1L] <- x
      put <- put + 1L
    }
  }
  L$Lindex[put + 1L] <- -1L
  put <- put + 1L
  L$Lbegin_p[rank + 2L] <- put
  L$Ubegin[rank + 2L] <- L$Ubegin[rank + 1L]
  for (pos in cbeg:(cend - 1L)) {
    i <- L$Windex[pos + 1L]
    if (i == pr) next
    w <- L$Wbegin[m + i + 1L]
    while (L$Windex[w + 1L] != pc) w <- w + 1L
    e <- L$Wend[m + i + 1L] - 1L
    L$Wend[m + i + 1L] <- e
    L$Windex[w + 1L] <- L$Windex[e + 1L]
    nz <- e - L$Wbegin[m + i + 1L]
    .blu_list_move(L, "rowcount_flink", "rowcount_blink", i, nz, m, "min_rownz")
  }
  .blu_pivot_cleanup(L, pivot, cbeg, rbeg)
}

.blu_pivot_singleton_col <- function(L) {
  m <- L$m; rank <- L$rank; droptol <- L$droptol
  pc <- L$pivot_col; pr <- L$pivot_row
  cbeg <- L$Wbegin[pc + 1L]
  rbeg <- L$Wbegin[m + pr + 1L]; rend <- L$Wend[m + pr + 1L]
  put <- L$Ubegin[rank + 1L]
  pivot <- L$Wvalue[cbeg + 1L]
  xrj <- 0
  for (rpos in rbeg:(rend - 1L)) {
    j <- L$Windex[rpos + 1L]
    if (j == pc) next
    where <- -1L; cmx <- 0
    for (pos in L$Wbegin[j + 1L]:(L$Wend[j + 1L] - 1L)) {
      if (L$Windex[pos + 1L] == pr) { where <- pos; xrj <- L$Wvalue[pos + 1L] }
      else {
        x <- abs(L$Wvalue[pos + 1L])
        if (x > cmx) cmx <- x
      }
    }
    if (abs(xrj) > droptol) {
      L$Uindex[put + 1L] <- j
      L$Uvalue[put + 1L] <- xrj
      put <- put + 1L
    }
    e <- L$Wend[j + 1L] - 1L
    L$Wend[j + 1L] <- e
    L$Windex[where + 1L] <- L$Windex[e + 1L]
    L$Wvalue[where + 1L] <- L$Wvalue[e + 1L]
    nz <- e - L$Wbegin[j + 1L]
    .blu_list_move(L, "colcount_flink", "colcount_blink", j, nz, m, "min_colnz")
    L$col_pivot[j + 1L] <- cmx
  }
  L$Ubegin[rank + 2L] <- put
  put <- L$Lbegin_p[rank + 1L]
  L$Lindex[put + 1L] <- -1L
  L$Lbegin_p[rank + 2L] <- put + 1L
  .blu_pivot_cleanup(L, pivot, cbeg, rbeg)
}

.blu_pivot_doubleton_col <- function(L) {
  m <- L$m; rank <- L$rank; droptol <- L$droptol
  pad <- L$pad; stretch <- L$stretch
  pc <- L$pivot_col; pr <- L$pivot_row
  cbeg <- L$Wbegin[pc + 1L]
  rbeg <- L$Wbegin[m + pr + 1L]; rend <- L$Wend[m + pr + 1L]
  if (L$Windex[cbeg + 1L] != pr) {
    t <- L$Windex[cbeg + 1L]; L$Windex[cbeg + 1L] <- L$Windex[cbeg + 2L]; L$Windex[cbeg + 2L] <- t
    t <- L$Wvalue[cbeg + 1L]; L$Wvalue[cbeg + 1L] <- L$Wvalue[cbeg + 2L]; L$Wvalue[cbeg + 2L] <- t
  }
  pivot <- L$Wvalue[cbeg + 1L]
  other_row <- L$Windex[cbeg + 2L]
  other_value <- L$Wvalue[cbeg + 2L]
  where <- rbeg
  while (L$Windex[where + 1L] != pc) where <- where + 1L
  t <- L$Windex[rbeg + 1L]; L$Windex[rbeg + 1L] <- L$Windex[where + 1L]; L$Windex[where + 1L] <- t
  nz <- L$Wend[m + other_row + 1L] - L$Wbegin[m + other_row + 1L]
  grow <- as.integer(nz + (rend - rbeg - 1L) + stretch * (nz + rend - rbeg - 1L) + pad)
  room <- L$Wend[2L * m + 1L] - L$Wbegin[2L * m + 1L]
  if (grow > room) {
    .blu_file_compress(L, 2L * m, stretch, pad)
    cbeg <- L$Wbegin[pc + 1L]
    rbeg <- L$Wbegin[m + pr + 1L]; rend <- L$Wend[m + pr + 1L]
    L$ngarbage <- L$ngarbage + 1L
  }
  Uput <- L$Ubegin[rank + 1L]
  put <- rbeg + 1L
  ncancelled <- 0L
  marked <- L$iwork0
  q <- other_value / pivot
  if (rend > rbeg + 1L) for (rpos in (rbeg + 1L):(rend - 1L)) {
    j <- L$Windex[rpos + 1L]
    cmx <- 0
    where_pivot <- -1L; where_other <- -1L
    end <- L$Wend[j + 1L]
    for (pos in L$Wbegin[j + 1L]:(end - 1L)) {
      ix <- L$Windex[pos + 1L]
      if (ix == pr) where_pivot <- pos
      else if (ix == other_row) where_other <- pos
      else {
        x <- abs(L$Wvalue[pos + 1L])
        if (x > cmx) cmx <- x
      }
    }
    xrj <- L$Wvalue[where_pivot + 1L]
    if (abs(xrj) > droptol) {
      L$Uindex[Uput + 1L] <- j
      L$Uvalue[Uput + 1L] <- xrj
      Uput <- Uput + 1L
    }
    if (where_other == -1L) {
      x <- -xrj * q
      xabs <- abs(x)
      if (xabs > droptol) {
        L$Windex[where_pivot + 1L] <- other_row
        L$Wvalue[where_pivot + 1L] <- x
        L$Windex[put + 1L] <- j
        put <- put + 1L
        if (xabs > cmx) cmx <- xabs
      } else {
        end <- L$Wend[j + 1L] - 1L
        L$Wend[j + 1L] <- end
        L$Windex[where_pivot + 1L] <- L$Windex[end + 1L]
        L$Wvalue[where_pivot + 1L] <- L$Wvalue[end + 1L]
        nz <- end - L$Wbegin[j + 1L]
        .blu_list_move(L, "colcount_flink", "colcount_blink", j, nz, m, "min_colnz")
      }
    } else {
      end <- L$Wend[j + 1L] - 1L
      L$Wend[j + 1L] <- end
      L$Windex[where_pivot + 1L] <- L$Windex[end + 1L]
      L$Wvalue[where_pivot + 1L] <- L$Wvalue[end + 1L]
      if (where_other == end) where_other <- where_pivot
      L$Wvalue[where_other + 1L] <- .ipx_madd(-xrj, q, L$Wvalue[where_other + 1L])
      x <- abs(L$Wvalue[where_other + 1L])
      if (x <= droptol) {
        end <- L$Wend[j + 1L] - 1L
        L$Wend[j + 1L] <- end
        L$Windex[where_other + 1L] <- L$Windex[end + 1L]
        L$Wvalue[where_other + 1L] <- L$Wvalue[end + 1L]
        marked[j + 1L] <- 1L
        ncancelled <- ncancelled + 1L
      } else if (x > cmx) cmx <- x
      nz <- L$Wend[j + 1L] - L$Wbegin[j + 1L]
      .blu_list_move(L, "colcount_flink", "colcount_blink", j, nz, m, "min_colnz")
    }
    L$col_pivot[j + 1L] <- cmx
  }
  rend <- put
  L$Ubegin[rank + 2L] <- Uput
  ob <- L$Wbegin[m + other_row + 1L]; oe <- L$Wend[m + other_row + 1L]
  if (ncancelled) {
    marked[pc + 1L] <- 1L
    put <- ob
    for (pos in ob:(oe - 1L)) {
      j <- L$Windex[pos + 1L]
      if (marked[j + 1L]) marked[j + 1L] <- 0L
      else { L$Windex[put + 1L] <- j; put <- put + 1L }
    }
    L$Wend[m + other_row + 1L] <- put
  } else {
    where <- ob
    while (L$Windex[where + 1L] != pc) where <- where + 1L
    end <- oe - 1L
    L$Wend[m + other_row + 1L] <- end
    L$Windex[where + 1L] <- L$Windex[end + 1L]
  }
  L$iwork0 <- marked
  nfill <- rend - (rbeg + 1L)
  room <- L$Wbegin[L$Wflink[m + other_row + 1L] + 1L] - L$Wend[m + other_row + 1L]
  if (nfill > room) {
    nz <- L$Wend[m + other_row + 1L] - L$Wbegin[m + other_row + 1L]
    space <- as.integer(nfill + stretch * (nz + nfill) + pad)
    .blu_file_reappend(L, m + other_row, 2L * m, space)
    L$nexpand <- L$nexpand + 1L
  }
  put <- L$Wend[m + other_row + 1L]
  if (nfill > 0L) {
    L$Windex[put + seq_len(nfill)] <- L$Windex[rbeg + 1L + seq_len(nfill)]
    put <- put + nfill
  }
  L$Wend[m + other_row + 1L] <- put
  nz <- put - L$Wbegin[m + other_row + 1L]
  .blu_list_move(L, "rowcount_flink", "rowcount_blink", other_row, nz, m, "min_rownz")
  put <- L$Lbegin_p[rank + 1L]
  x <- other_value / pivot
  if (abs(x) > droptol) {
    L$Lindex[put + 1L] <- other_row
    L$Lvalue[put + 1L] <- x
    put <- put + 1L
  }
  L$Lindex[put + 1L] <- -1L
  L$Lbegin_p[rank + 2L] <- put + 1L
  .blu_pivot_cleanup(L, pivot, cbeg, rbeg)
}

.blu_remove_col <- function(L, j) {
  m <- L$m
  cbeg <- L$Wbegin[j + 1L]; cend <- L$Wend[j + 1L]
  if (cend > cbeg) for (pos in cbeg:(cend - 1L)) {
    i <- L$Windex[pos + 1L]
    w <- L$Wbegin[m + i + 1L]
    while (L$Windex[w + 1L] != j) w <- w + 1L
    e <- L$Wend[m + i + 1L] - 1L
    L$Wend[m + i + 1L] <- e
    L$Windex[w + 1L] <- L$Windex[e + 1L]
    nz <- e - L$Wbegin[m + i + 1L]
    .blu_list_move(L, "rowcount_flink", "rowcount_blink", i, nz, m, "min_rownz")
  }
  L$col_pivot[j + 1L] <- 0
  L$Wend[j + 1L] <- cbeg
  .blu_list_move(L, "colcount_flink", "colcount_blink", j, 0L, m, "min_colnz")
}

.blu_build_factors <- function(L) {
  m <- L$m; rank <- L$rank; pad <- L$pad; stretch <- L$stretch
  Lnz <- L$Lbegin_p[rank + 1L] - rank
  Unz <- L$Ubegin[rank + 1L]
  pinv <- L$pinv; qinv <- L$qinv
  pivotrow <- integer(m); pivotcol <- integer(m)
  lrank <- rank
  for (i in 0:(m - 1L)) {
    if (pinv[i + 1L] < 0L) { pinv[i + 1L] <- lrank; lrank <- lrank + 1L }
    pivotrow[pinv[i + 1L] + 1L] <- i
  }
  lrank <- rank
  for (j in 0:(m - 1L)) {
    if (qinv[j + 1L] < 0L) { qinv[j + 1L] <- lrank; lrank <- lrank + 1L }
    pivotcol[qinv[j + 1L] + 1L] <- j
  }
  if (rank < m) L$col_pivot[pivotcol[(rank:(m - 1L)) + 1L] + 1L] <- 1
  put <- L$Lbegin_p[rank + 1L]
  k <- rank
  while (k < m) {
    L$Lindex[put + 1L] <- -1L
    put <- put + 1L
    L$Lbegin_p[k + 2L] <- put
    k <- k + 1L
  }
  Lbegin_p <- L$Lbegin_p
  L$Lbegin <- Lbegin_p[pinv + 1L]
  Lindex <- L$Lindex[seq_len(Lnz + m)]
  Lvalue <- L$Lvalue[seq_len(Lnz + m)]
  Lvalue[is.na(Lvalue)] <- 0
  cnt <- tabulate(Lindex[Lindex >= 0L] + 1L, m)
  put <- Lnz + m
  Ltbegin_p <- integer(m); Ltbegin <- integer(m); nxt <- integer(m)
  Lindex2 <- integer(Lnz + m); Lvalue2 <- numeric(Lnz + m)
  for (k in 0:(m - 1L)) {
    i <- pivotrow[k + 1L]
    Ltbegin_p[k + 1L] <- put
    Ltbegin[i + 1L] <- put
    put <- put + cnt[i + 1L]
    Lindex2[put - (Lnz + m) + 1L] <- -1L
    put <- put + 1L
    nxt[i + 1L] <- Ltbegin_p[k + 1L]
  }
  for (k in 0:(m - 1L)) {
    ipivot <- pivotrow[k + 1L]
    get <- Lbegin_p[k + 1L]
    while ((i <- Lindex[get + 1L]) >= 0L) {
      put <- nxt[i + 1L]; nxt[i + 1L] <- put + 1L
      Lindex2[put - (Lnz + m) + 1L] <- ipivot
      Lvalue2[put - (Lnz + m) + 1L] <- Lvalue[get + 1L]
      get <- get + 1L
    }
  }
  L$Lindex <- c(Lindex, Lindex2)
  L$Lvalue <- c(Lvalue, Lvalue2)
  L$Ltbegin <- Ltbegin; L$Ltbegin_p <- Ltbegin_p
  L$Rbegin <- integer(m + 1L); L$Rbegin[1L] <- 2L * (Lnz + m)
  L$eta_row <- integer(m + 1L)
  # U rowwise
  Uindex <- L$Uindex; Uvalue <- L$Uvalue; Ubegin <- L$Ubegin
  L$Wbegin <- integer(m + 1L); L$Wend <- integer(m + 1L)
  L$Wflink <- integer(m + 1L); L$Wblink <- integer(m + 1L)
  .blu_file_empty(L, m, L$Wmem)
  Windex <- integer(0); Wvalue <- numeric(0)
  colcnt <- integer(m)
  put <- 0L
  if (rank == m) {
    for (k in 0:(m - 1L)) {
      jpivot <- pivotcol[k + 1L]
      L$Wbegin[jpivot + 1L] <- put
      nz <- Ubegin[k + 2L] - Ubegin[k + 1L]
      if (nz > 0L) {
        ps <- Ubegin[k + 1L] + seq_len(nz)
        jj <- Uindex[ps]
        Windex[put + seq_len(nz)] <- jj
        Wvalue[put + seq_len(nz)] <- Uvalue[ps]
        colcnt[jj + 1L] <- colcnt[jj + 1L] + 1L
        put <- put + nz
      }
      L$Wend[jpivot + 1L] <- put
      put <- as.integer(put + stretch * nz + pad)
      .blu_list_move(L, "Wflink", "Wblink", jpivot, 0L, m)
    }
  } else {
    Unz <- 0L
    for (k in seq_len(rank) - 1L) {
      jpivot <- pivotcol[k + 1L]
      L$Wbegin[jpivot + 1L] <- put
      nz <- 0L
      b <- Ubegin[k + 1L]; e <- Ubegin[k + 2L]
      if (e > b) for (pos in b:(e - 1L)) {
        j <- Uindex[pos + 1L]
        if (qinv[j + 1L] < rank) {
          Windex[put + 1L] <- j
          Wvalue[put + 1L] <- Uvalue[pos + 1L]
          put <- put + 1L
          colcnt[j + 1L] <- colcnt[j + 1L] + 1L
          nz <- nz + 1L
        }
      }
      L$Wend[jpivot + 1L] <- put
      put <- as.integer(put + stretch * nz + pad)
      .blu_list_move(L, "Wflink", "Wblink", jpivot, 0L, m)
      Unz <- Unz + nz
    }
    for (k in rank:(m - 1L)) {
      jpivot <- pivotcol[k + 1L]
      L$Wbegin[jpivot + 1L] <- put
      L$Wend[jpivot + 1L] <- put
      put <- put + pad
      .blu_list_move(L, "Wflink", "Wblink", jpivot, 0L, m)
    }
  }
  L$Wbegin[m + 1L] <- put
  L$Windex <- Windex; L$Wvalue <- Wvalue
  # U columnwise
  Uindex <- -1L; Uvalue <- 0
  nUbegin <- integer(m + 1L)
  put <- 1L
  for (k in 0:(m - 1L)) {
    j <- pivotcol[k + 1L]; i <- pivotrow[k + 1L]
    nz <- colcnt[j + 1L]
    if (nz == 0L) nUbegin[i + 1L] <- 0L
    else {
      nUbegin[i + 1L] <- put
      put <- put + nz
      Uindex[put + 1L] <- -1L
      Uvalue[put + 1L] <- 0
      put <- put + 1L
    }
    colcnt[j + 1L] <- nUbegin[i + 1L]
  }
  nUbegin[m + 1L] <- put
  for (k in 0:(m - 1L)) {
    jpivot <- pivotcol[k + 1L]; i <- pivotrow[k + 1L]
    b <- L$Wbegin[jpivot + 1L]; e <- L$Wend[jpivot + 1L]
    if (e > b) for (pos in b:(e - 1L)) {
      j <- Windex[pos + 1L]
      put <- colcnt[j + 1L]; colcnt[j + 1L] <- put + 1L
      Uindex[put + 1L] <- i
      Uvalue[put + 1L] <- Wvalue[pos + 1L]
    }
  }
  L$Uindex <- Uindex; L$Uvalue <- Uvalue; L$Ubegin <- nUbegin
  pmap <- integer(m); qmap <- integer(m)
  pmap[pivotcol + 1L] <- pivotrow
  qmap[pivotrow + 1L] <- pivotcol
  L$pmap <- pmap; L$qmap <- qmap
  L$row_pivot <- L$col_pivot[qmap + 1L]
  ap <- abs(L$row_pivot)
  L$max_pivot <- max(ap, 0); L$min_pivot <- min(ap, Inf)
  L$p <- pivotrow
  L$pivotrow <- pivotrow; L$pivotcol <- pivotcol
  L$pivotlen <- m
  L$Lnz <- Lnz; L$Unz <- Unz; L$Rnz <- 0L
  L$pinv <- pinv; L$qinv <- qinv
}

# sequential sum of absolute values (lu_onenorm)
.blu_onenorm <- function(x) { d <- 0; for (v in x) d <- d + abs(v); d }

.blu_residual_test <- function(L, Bcols) {
  m <- L$m; rank <- L$rank; p <- L$p
  pivotcol <- L$pivotcol; pivotrow <- L$pivotrow
  Lindex <- L$Lindex; Lvalue <- L$Lvalue; Uindex <- L$Uindex; Uvalue <- L$Uvalue
  Ubegin <- L$Ubegin; row_pivot <- L$row_pivot
  rhs <- numeric(m); lhs <- numeric(m)
  for (k in 0:(m - 1L)) {
    d <- 0
    pos <- L$Ltbegin_p[k + 1L]
    while ((i <- Lindex[pos + 1L]) >= 0L) { d <- .ipx_madd(lhs[i + 1L], Lvalue[pos + 1L], d); pos <- pos + 1L }
    ip <- p[k + 1L]
    rhs[ip + 1L] <- if (d <= 0) 1 else -1
    lhs[ip + 1L] <- rhs[ip + 1L] - d
  }
  for (k in (m - 1L):0) {
    ip <- pivotrow[k + 1L]
    d <- lhs[ip + 1L] / row_pivot[ip + 1L]
    lhs[ip + 1L] <- d
    pos <- Ubegin[ip + 1L]
    while ((i <- Uindex[pos + 1L]) >= 0L) { lhs[i + 1L] <- .ipx_madd(-d, Uvalue[pos + 1L], lhs[i + 1L]); pos <- pos + 1L }
  }
  for (k in seq_len(rank) - 1L) {
    ip <- pivotrow[k + 1L]; jp <- pivotcol[k + 1L]
    d <- lhs[ip + 1L]
    cl <- Bcols[[jp + 1L]]
    if (length(cl$i)) rhs[cl$i + 1L] <- .ipx_madd(-d, cl$x, rhs[cl$i + 1L])
  }
  if (rank < m) for (k in rank:(m - 1L)) {
    ip <- pivotrow[k + 1L]
    rhs[ip + 1L] <- rhs[ip + 1L] - lhs[ip + 1L]
  }
  norm_ftran <- .blu_onenorm(lhs)
  norm_ftran_res <- .blu_onenorm(rhs)
  for (k in 0:(m - 1L)) {
    ip <- pivotrow[k + 1L]
    d <- 0
    pos <- Ubegin[ip + 1L]
    while ((i <- Uindex[pos + 1L]) >= 0L) { d <- .ipx_madd(lhs[i + 1L], Uvalue[pos + 1L], d); pos <- pos + 1L }
    rhs[ip + 1L] <- if (d <= 0) 1 else -1
    lhs[ip + 1L] <- (rhs[ip + 1L] - d) / row_pivot[ip + 1L]
  }
  for (k in (m - 1L):0) {
    d <- 0
    pos <- L$Lbegin_p[k + 1L]
    while ((i <- Lindex[pos + 1L]) >= 0L) { d <- .ipx_madd(lhs[i + 1L], Lvalue[pos + 1L], d); pos <- pos + 1L }
    lhs[p[k + 1L] + 1L] <- lhs[p[k + 1L] + 1L] - d
  }
  for (k in seq_len(rank) - 1L) {
    ip <- pivotrow[k + 1L]; jp <- pivotcol[k + 1L]
    cl <- Bcols[[jp + 1L]]
    d <- 0
    for (t in seq_along(cl$i)) d <- .ipx_madd(lhs[cl$i[t] + 1L], cl$x[t], d)
    rhs[ip + 1L] <- rhs[ip + 1L] - d
  }
  if (rank < m) for (k in rank:(m - 1L)) {
    ip <- pivotrow[k + 1L]
    rhs[ip + 1L] <- rhs[ip + 1L] - lhs[ip + 1L]
  }
  norm_btran <- .blu_onenorm(lhs)
  norm_btran_res <- .blu_onenorm(rhs)
  # lu_matrix_norm
  rowsum <- numeric(m); onenorm <- 0
  for (k in seq_len(rank) - 1L) {
    cl <- Bcols[[pivotcol[k + 1L] + 1L]]
    colsum <- 0
    for (t in seq_along(cl$i)) {
      colsum <- colsum + abs(cl$x[t])
      rowsum[cl$i[t] + 1L] <- rowsum[cl$i[t] + 1L] + abs(cl$x[t])
    }
    onenorm <- max(onenorm, colsum)
  }
  if (rank < m) for (k in rank:(m - 1L)) {
    ip <- pivotrow[k + 1L]
    rowsum[ip + 1L] <- rowsum[ip + 1L] + 1
    onenorm <- max(onenorm, 1)
  }
  L$onenorm <- onenorm
  L$infnorm <- max(rowsum, 0)
  L$residual_test <- max(norm_ftran_res / .ipx_madd(onenorm, norm_ftran, m),
                         norm_btran_res / .ipx_madd(L$infnorm, norm_btran, m))
}

.blu_garbage_perm <- function(L) {
  m <- L$m
  if (L$pivotlen > m) {
    pc <- L$pivotcol[seq_len(L$pivotlen)]; pr <- L$pivotrow[seq_len(L$pivotlen)]
    keep <- !duplicated(pc, fromLast = TRUE)
    L$pivotcol <- pc[keep]; L$pivotrow <- pr[keep]
    L$pivotlen <- m
  }
}

# lu_solve_dense
.blu_solve_dense <- function(L, rhs, trans) {
  m <- L$m
  .blu_garbage_perm(L)
  pivotcol <- L$pivotcol; pivotrow <- L$pivotrow
  Lindex <- L$Lindex; Lvalue <- L$Lvalue
  work1 <- rhs
  lhs <- numeric(m)
  if (trans == "T") {
    Wbegin <- L$Wbegin; Wend <- L$Wend; Windex <- L$Windex; Wvalue <- L$Wvalue
    for (k in 0:(m - 1L)) {
      jp <- pivotcol[k + 1L]; ip <- pivotrow[k + 1L]
      x <- work1[jp + 1L] / L$col_pivot[jp + 1L]
      b <- Wbegin[jp + 1L]; e <- Wend[jp + 1L]
      if (e > b) {
        ps <- (b + 1L):e
        jj <- Windex[ps] + 1L
        work1[jj] <- .ipx_madd(-x, Wvalue[ps], work1[jj])
      }
      lhs[ip + 1L] <- x
    }
    t <- L$nforrest - 1L
    while (t >= 0L) {
      ip <- L$eta_row[t + 1L]
      x <- lhs[ip + 1L]
      b <- L$Rbegin[t + 1L]; e <- L$Rbegin[t + 2L]
      if (e > b) {
        ps <- (b + 1L):e
        ii <- Lindex[ps] + 1L
        lhs[ii] <- .ipx_madd(-x, Lvalue[ps], lhs[ii])
      }
      t <- t - 1L
    }
    for (k in (m - 1L):0) {
      x <- 0
      pos <- L$Lbegin_p[k + 1L]
      while ((i <- Lindex[pos + 1L]) >= 0L) { x <- .ipx_madd(lhs[i + 1L], Lvalue[pos + 1L], x); pos <- pos + 1L }
      ip <- L$p[k + 1L]
      lhs[ip + 1L] <- lhs[ip + 1L] - x
    }
  } else {
    for (k in 0:(m - 1L)) {
      x <- 0
      pos <- L$Ltbegin_p[k + 1L]
      while ((i <- Lindex[pos + 1L]) >= 0L) { x <- .ipx_madd(work1[i + 1L], Lvalue[pos + 1L], x); pos <- pos + 1L }
      ip <- L$p[k + 1L]
      work1[ip + 1L] <- work1[ip + 1L] - x
    }
    pos <- L$Rbegin[1L]
    for (t in seq_len(L$nforrest) - 1L) {
      ip <- L$eta_row[t + 1L]
      x <- 0
      while (pos < L$Rbegin[t + 2L]) { x <- .ipx_madd(work1[Lindex[pos + 1L] + 1L], Lvalue[pos + 1L], x); pos <- pos + 1L }
      work1[ip + 1L] <- work1[ip + 1L] - x
    }
    Uindex <- L$Uindex; Uvalue <- L$Uvalue; Ubegin <- L$Ubegin
    for (k in (m - 1L):0) {
      jp <- pivotcol[k + 1L]; ip <- pivotrow[k + 1L]
      x <- work1[ip + 1L] / L$row_pivot[ip + 1L]
      pos <- Ubegin[ip + 1L]
      while ((i <- Uindex[pos + 1L]) >= 0L) { work1[i + 1L] <- .ipx_madd(-x, Uvalue[pos + 1L], work1[i + 1L]); pos <- pos + 1L }
      lhs[jp + 1L] <- x
    }
  }
  lhs
}

# lu_dfs: begin/index with end (row file) or -1 terminated lines
.blu_dfs <- function(L, i, begin, end, index, top, xi, marked, M) {
  if (marked[i + 1L] == M) return(list(top = top, xi = xi, marked = marked))
  stack <- i; pstack <- 0L
  head <- 1L
  while (head >= 1L) {
    i <- stack[head]
    if (marked[i + 1L] != M) {
      marked[i + 1L] <- M
      pstack[head] <- begin[i + 1L]
    }
    done <- TRUE
    p <- pstack[head]
    repeat {
      if (is.null(end)) { inext <- index[p + 1L]; if (inext < 0L) break }
      else { if (p >= end[i + 1L]) break; inext <- index[p + 1L] }
      if (marked[inext + 1L] == M) { p <- p + 1L; next }
      pstack[head] <- p + 1L
      head <- head + 1L
      stack[head] <- inext
      done <- FALSE
      break
    }
    if (done) {
      head <- head - 1L
      top <- top - 1L
      xi[top + 1L] <- i
    }
  }
  list(top = top, xi = xi, marked = marked)
}

.blu_solve_symbolic <- function(L, begin, end, index, irhs, M) {
  m <- L$m
  top <- m; xi <- integer(m); marked <- L$iwork0
  for (i in irhs) if (marked[i + 1L] != M) {
    r <- .blu_dfs(L, i, begin, end, index, top, xi, marked, M)
    top <- r$top; xi <- r$xi; marked <- r$marked
  }
  L$iwork0 <- marked
  if (top < m) xi[(top + 1L):m] else integer(0)
}

# lu_solve_triangular on lhs; returns list(lhs, pattern)
.blu_solve_triangular <- function(pattern_symb, begin, end, index, value, pivot,
                                  droptol, lhs) {
  pattern <- integer(0)
  for (ip in pattern_symb) {
    if (lhs[ip + 1L] != 0) {
      if (!is.null(pivot)) lhs[ip + 1L] <- lhs[ip + 1L] / pivot[ip + 1L]
      x <- lhs[ip + 1L]
      if (!is.null(end)) {
        b <- begin[ip + 1L]; e <- end[ip + 1L]
        if (e > b) {
          ps <- (b + 1L):e
          ii <- index[ps] + 1L
          lhs[ii] <- .ipx_madd(-x, value[ps], lhs[ii])
        }
      } else {
        pos <- begin[ip + 1L]
        while ((i <- index[pos + 1L]) >= 0L) { lhs[i + 1L] <- .ipx_madd(-x, value[pos + 1L], lhs[i + 1L]); pos <- pos + 1L }
      }
      if (abs(x) > droptol) pattern <- c(pattern, ip) else lhs[ip + 1L] <- 0
    }
  }
  list(lhs = lhs, pattern = pattern)
}

# lu_solve_for_update; returns the dense solution in column (ftran) or row
# (btran) indexing. irhs, xrhs: 0-based indices and values of the rhs; for
# btran irhs is the column jpivot.
.blu_solve_for_update <- function(L, irhs, xrhs, trans, want = TRUE) {
  m <- L$m; nforrest <- L$nforrest; droptol <- L$droptol
  nz_sparse <- as.integer(L$sparse_thres * m)
  Rflops <- 0
  xlhs <- numeric(m)
  if (trans == "T") {
    jpivot <- irhs
    ipivot <- L$pmap[jpivot + 1L]
    jbegin <- L$Wbegin[jpivot + 1L]; jend <- L$Wend[jpivot + 1L]
    L$marker <- L$marker + 1L; M <- L$marker
    rowpat <- if (jend > jbegin) L$Windex[(jbegin + 1L):jend] else integer(0)
    psymb <- .blu_solve_symbolic(L, L$Wbegin, L$Wend, L$Windex, rowpat, M)
    work <- L$work0
    if (jend > jbegin) work[rowpat + 1L] <- L$Wvalue[(jbegin + 1L):jend]
    tr <- .blu_solve_triangular(psymb, L$Wbegin, L$Wend, L$Windex, L$Wvalue,
                                L$col_pivot, 0, work)
    work <- tr$lhs
    put <- L$Rbegin[nforrest + 1L]
    for (j in psymb) {
      L$Lindex[put + 1L] <- L$pmap[j + 1L]
      L$Lvalue[put + 1L] <- work[j + 1L]
      put <- put + 1L
      work[j + 1L] <- 0
    }
    L$work0 <- work
    L$Rbegin[nforrest + 2L] <- put
    L$eta_row[nforrest + 1L] <- ipivot
    L$btran_for_update <- jpivot
    if (!want) return(invisible(NULL))
    # solution
    L$marker <- L$marker + 1L; M <- L$marker
    marked <- L$iwork0
    pattern <- ipivot
    marked[ipivot + 1L] <- M
    pivot <- L$col_pivot[jpivot + 1L]
    xlhs[ipivot + 1L] <- 1 / pivot
    xdrop <- droptol * abs(pivot)
    b <- L$Rbegin[nforrest + 1L]
    if (put > b) for (pos in b:(put - 1L)) {
      if (abs(L$Lvalue[pos + 1L]) > xdrop) {
        i <- L$Lindex[pos + 1L]
        pattern <- c(pattern, i)
        marked[i + 1L] <- M
        xlhs[i + 1L] <- -L$Lvalue[pos + 1L] / pivot
      }
    }
    t <- nforrest - 1L
    while (t >= 0L) {
      ip <- L$eta_row[t + 1L]
      if (xlhs[ip + 1L] != 0) {
        x <- xlhs[ip + 1L]
        b <- L$Rbegin[t + 1L]; e <- L$Rbegin[t + 2L]
        if (e > b) for (pos in b:(e - 1L)) {
          i <- L$Lindex[pos + 1L]
          if (marked[i + 1L] != M) { marked[i + 1L] <- M; pattern <- c(pattern, i) }
          xlhs[i + 1L] <- .ipx_madd(-x, L$Lvalue[pos + 1L], xlhs[i + 1L])
          Rflops <- Rflops + 1
        }
      }
      t <- t - 1L
    }
    L$iwork0 <- marked
    if (length(pattern) <= nz_sparse) {
      L$marker <- L$marker + 1L; M <- L$marker
      psymb <- .blu_solve_symbolic(L, L$Ltbegin, NULL, L$Lindex, pattern, M)
      tr <- .blu_solve_triangular(psymb, L$Ltbegin, NULL, L$Lindex, L$Lvalue,
                                  NULL, droptol, xlhs)
      xlhs <- tr$lhs
      outpat <- tr$pattern
    } else {
      outpat <- integer(0)
      for (k in (m - 1L):0) {
        ip <- L$p[k + 1L]
        if (xlhs[ip + 1L] != 0) {
          x <- xlhs[ip + 1L]
          pos <- L$Ltbegin_p[k + 1L]
          while ((i <- L$Lindex[pos + 1L]) >= 0L) {
            xlhs[i + 1L] <- .ipx_madd(-x, L$Lvalue[pos + 1L], xlhs[i + 1L])
            pos <- pos + 1L
          }
          if (abs(x) > droptol) outpat <- c(outpat, ip) else xlhs[ip + 1L] <- 0
        }
      }
    }
  } else {
    L$marker <- L$marker + 1L; M <- L$marker
    psymb <- .blu_solve_symbolic(L, L$Lbegin, NULL, L$Lindex, irhs, M)
    work <- L$work0
    work[irhs + 1L] <- xrhs
    tr <- .blu_solve_triangular(psymb, L$Lbegin, NULL, L$Lindex, L$Lvalue, NULL,
                                droptol, work)
    work <- tr$lhs
    pattern <- tr$pattern
    marked <- L$iwork0
    if (length(pattern) < length(psymb)) {
      gone <- setdiff(psymb, pattern)
      marked[gone + 1L] <- marked[gone + 1L] - 1L
    }
    pos <- L$Rbegin[1L]
    for (t in seq_len(nforrest) - 1L) {
      ip <- L$eta_row[t + 1L]
      x <- 0
      while (pos < L$Rbegin[t + 2L]) {
        x <- .ipx_madd(work[L$Lindex[pos + 1L] + 1L], L$Lvalue[pos + 1L], x)
        pos <- pos + 1L
      }
      work[ip + 1L] <- work[ip + 1L] - x
      if (x != 0 && marked[ip + 1L] != M) { marked[ip + 1L] <- M; pattern <- c(pattern, ip) }
    }
    L$iwork0 <- marked
    Rflops <- Rflops + L$Rbegin[nforrest + 1L] - L$Rbegin[1L]
    put <- L$Ubegin[m + 1L]
    for (i in pattern) {
      L$Uindex[put + 1L] <- i
      L$Uvalue[put + 1L] <- work[i + 1L]
      put <- put + 1L
    }
    L$Uindex[put + 1L] <- -1L
    L$Uvalue[put + 1L] <- 0
    L$ftran_for_update <- 0L
    if (!want) {
      L$work0 <- numeric(m)
      L$update_cost_numer <- L$update_cost_numer + Rflops
      return(invisible(NULL))
    }
    if (length(pattern) <= nz_sparse) {
      L$marker <- L$marker + 1L; M <- L$marker
      psymb <- .blu_solve_symbolic(L, L$Ubegin, NULL, L$Uindex, pattern, M)
      tr <- .blu_solve_triangular(psymb, L$Ubegin, NULL, L$Uindex, L$Uvalue,
                                  L$row_pivot, droptol, work)
      outpat <- L$qmap[tr$pattern + 1L]
      xlhs[outpat + 1L] <- tr$lhs[tr$pattern + 1L]
      work <- tr$lhs
      work[tr$pattern + 1L] <- 0
      work[psymb + 1L] <- 0
    } else {
      outpat <- integer(0)
      for (k in (L$pivotlen - 1L):0) {
        ip <- L$pivotrow[k + 1L]; jp <- L$pivotcol[k + 1L]
        if (work[ip + 1L] != 0) {
          x <- work[ip + 1L] / L$row_pivot[ip + 1L]
          work[ip + 1L] <- 0
          pos <- L$Ubegin[ip + 1L]
          while ((i <- L$Uindex[pos + 1L]) >= 0L) {
            work[i + 1L] <- .ipx_madd(-x, L$Uvalue[pos + 1L], work[i + 1L])
            pos <- pos + 1L
          }
          if (abs(x) > droptol) { xlhs[jp + 1L] <- x; outpat <- c(outpat, jp) }
        }
      }
    }
    L$work0 <- numeric(m)
  }
  L$update_cost_numer <- L$update_cost_numer + Rflops
  L$lhs_pattern <- outpat
  xlhs
}

.blu_find <- function(j, index, start, end) {
  if (end >= 0L) {
    while (start < end && index[start + 1L] != j) start <- start + 1L
    start
  } else {
    while (index[start + 1L] != j && index[start + 1L] >= 0L) start <- start + 1L
    if (index[start + 1L] == j) start else end
  }
}

.blu_permute <- function(L, jlist) {
  nswap <- length(jlist) - 1L
  j0 <- jlist[1L]; jn <- jlist[nswap + 1L]
  i0 <- L$pmap[j0 + 1L]; inn <- L$pmap[jn + 1L]
  begin <- L$Wbegin[jn + 1L]; end <- L$Wend[jn + 1L]
  piv <- L$col_pivot[jn + 1L]
  for (n in nswap:1) {
    j <- jlist[n + 1L]; jprev <- jlist[n]
    L$Wbegin[j + 1L] <- L$Wbegin[jprev + 1L]
    L$Wend[j + 1L] <- L$Wend[jprev + 1L]
    .blu_list_swap(L, "Wflink", "Wblink", j, jprev)
    where <- .blu_find(j, L$Windex, L$Wbegin[j + 1L], L$Wend[j + 1L])
    if (n > 1L) {
      L$Windex[where + 1L] <- jprev
      L$col_pivot[j + 1L] <- L$Wvalue[where + 1L]
      L$Wvalue[where + 1L] <- L$col_pivot[jprev + 1L]
    } else {
      L$col_pivot[j + 1L] <- L$Wvalue[where + 1L]
      e <- L$Wend[j + 1L] - 1L
      L$Wend[j + 1L] <- e
      L$Windex[where + 1L] <- L$Windex[e + 1L]
      L$Wvalue[where + 1L] <- L$Wvalue[e + 1L]
    }
    L$min_pivot <- min(L$min_pivot, abs(L$col_pivot[j + 1L]))
    L$max_pivot <- max(L$max_pivot, abs(L$col_pivot[j + 1L]))
  }
  L$Wbegin[j0 + 1L] <- begin
  L$Wend[j0 + 1L] <- end
  where <- .blu_find(j0, L$Windex, L$Wbegin[j0 + 1L], L$Wend[j0 + 1L])
  L$Windex[where + 1L] <- jn
  L$col_pivot[j0 + 1L] <- L$Wvalue[where + 1L]
  L$Wvalue[where + 1L] <- piv
  L$min_pivot <- min(L$min_pivot, abs(L$col_pivot[j0 + 1L]))
  L$max_pivot <- max(L$max_pivot, abs(L$col_pivot[j0 + 1L]))
  begin <- L$Ubegin[i0 + 1L]
  for (n in seq_len(nswap) - 1L) {
    i <- L$pmap[jlist[n + 1L] + 1L]
    inext <- L$pmap[jlist[n + 2L] + 1L]
    L$Ubegin[i + 1L] <- L$Ubegin[inext + 1L]
    where <- .blu_find(i, L$Uindex, L$Ubegin[i + 1L], -1L)
    L$Uindex[where + 1L] <- inext
    L$row_pivot[i + 1L] <- L$Uvalue[where + 1L]
    L$Uvalue[where + 1L] <- L$row_pivot[inext + 1L]
  }
  L$Ubegin[inn + 1L] <- begin
  where <- .blu_find(inn, L$Uindex, L$Ubegin[inn + 1L], -1L)
  L$row_pivot[inn + 1L] <- L$Uvalue[where + 1L]
  end <- where
  while (L$Uindex[end + 1L] >= 0L) end <- end + 1L
  L$Uindex[where + 1L] <- L$Uindex[end]
  L$Uvalue[where + 1L] <- L$Uvalue[end]
  L$Uindex[end] <- -1L
  for (n in nswap:1) {
    j <- jlist[n + 1L]
    i <- L$pmap[jlist[n] + 1L]
    L$pmap[j + 1L] <- i
    L$qmap[i + 1L] <- j
  }
  L$pmap[j0 + 1L] <- inn
  L$qmap[inn + 1L] <- j0
}

.blu_bfs_path <- function(L, j0) {
  m <- L$m
  marked <- L$iwork0
  queue <- j0; front <- 1L; found <- FALSE; j <- j0
  while (front <= length(queue) && !found) {
    j <- queue[front]
    b <- L$Wbegin[j + 1L]; e <- L$Wend[j + 1L]
    if (e > b) for (pos in b:(e - 1L)) {
      k <- L$Windex[pos + 1L]
      if (k == j0) { found <- TRUE; break }
      if (marked[k + 1L] >= 0L) {
        marked[k + 1L] <- -j - 1L
        queue <- c(queue, k)
      }
    }
    front <- front + 1L
  }
  path <- integer(0)
  if (found) {
    while (j != j0) {
      path <- c(j, path)
      j <- -marked[j + 1L] - 1L
    }
    path <- c(j0, path)
  }
  marked[queue + 1L] <- 0L
  L$iwork0 <- marked
  path
}

.blu_compress_packed <- function(L) {
  m <- L$m
  begin <- L$Ubegin; index <- L$Uindex; value <- L$Uvalue
  end <- begin[m + 1L]
  for (i in 0:(m - 1L)) {
    p <- begin[i + 1L]
    if (index[p + 1L] == -1L) begin[i + 1L] <- 0L
    else {
      begin[i + 1L] <- index[p + 1L]
      index[p + 1L] <- -1L - i - 1L
    }
  }
  i <- -1L; put <- 1L
  if (end > 1L) for (get in 1:(end - 1L)) {
    ix <- index[get + 1L]
    if (ix > -1L) {
      index[put + 1L] <- ix; value[put + 1L] <- value[get + 1L]; put <- put + 1L
    } else if (ix < -1L) {
      i <- -1L - ix - 1L
      index[put + 1L] <- begin[i + 1L]
      begin[i + 1L] <- put
      value[put + 1L] <- value[get + 1L]
      put <- put + 1L
    } else if (i >= 0L) {
      i <- -1L
      index[put + 1L] <- -1L
      put <- put + 1L
    }
  }
  begin[m + 1L] <- put
  L$Ubegin <- begin; L$Uindex <- index; L$Uvalue <- value
}

# lu_update; returns "ok" or "singular"
.blu_update <- function(L, xtbl) {
  m <- L$m; nforrest <- L$nforrest; Unz <- L$Unz
  pad <- L$pad; stretch <- L$stretch
  jpivot <- L$btran_for_update
  ipivot <- L$pmap[jpivot + 1L]
  oldpiv <- L$col_pivot[jpivot + 1L]
  spike_diag <- 0; have_diag <- FALSE
  put <- L$Ubegin[m + 1L]
  pos <- put
  while ((i <- L$Uindex[pos + 1L]) >= 0L) {
    if (i != ipivot) {
      L$Uindex[put + 1L] <- i; L$Uvalue[put + 1L] <- L$Uvalue[pos + 1L]; put <- put + 1L
    } else {
      spike_diag <- L$Uvalue[pos + 1L]; have_diag <- TRUE
    }
    pos <- pos + 1L
  }
  if (have_diag) { L$Uindex[put + 1L] <- ipivot; L$Uvalue[put + 1L] <- spike_diag }
  nz_spike <- put - L$Ubegin[m + 1L]
  rb <- L$Rbegin[nforrest + 1L]; re <- L$Rbegin[nforrest + 2L]
  nz_roweta <- re - rb
  L$marker <- L$marker + 1L; M <- L$marker
  marked <- L$iwork0
  work1 <- numeric(m)
  if (re > rb) {
    ii <- L$Lindex[(rb + 1L):re]
    marked[ii + 1L] <- M
    work1[ii + 1L] <- L$Lvalue[(rb + 1L):re]
  }
  L$iwork0 <- marked
  newpiv <- spike_diag
  intersect <- 0L
  sb <- L$Ubegin[m + 1L]
  spike_i <- if (nz_spike) L$Uindex[sb + seq_len(nz_spike)] else integer(0)
  spike_x <- if (nz_spike) L$Uvalue[sb + seq_len(nz_spike)] else numeric(0)
  for (t in seq_len(nz_spike)) {
    i <- spike_i[t]
    if (marked[i + 1L] == M) {
      newpiv <- .ipx_madd(-spike_x[t], work1[i + 1L], newpiv)
      intersect <- intersect + 1L
    }
  }
  if (newpiv == 0 || abs(newpiv) < L$abstol) return("singular")
  piverr <- abs(.ipx_madd(-xtbl, oldpiv, newpiv))
  # remove column jpivot from row file
  nz <- 0L
  pos <- L$Ubegin[ipivot + 1L]
  while ((i <- L$Uindex[pos + 1L]) >= 0L) {
    j <- L$qmap[i + 1L]
    end <- L$Wend[j + 1L]
    L$Wend[j + 1L] <- end - 1L
    where <- .blu_find(jpivot, L$Windex, L$Wbegin[j + 1L], end)
    L$Windex[where + 1L] <- L$Windex[end]
    L$Wvalue[where + 1L] <- L$Wvalue[end]
    nz <- nz + 1L
    pos <- pos + 1L
  }
  Unz <- Unz - nz
  pos <- L$Ubegin[ipivot + 1L]
  while (L$Uindex[pos + 1L] >= 0L) { L$Uindex[pos + 1L] <- -1L; pos <- pos + 1L }
  L$Ubegin[ipivot + 1L] <- L$Ubegin[m + 1L]
  L$Ubegin[m + 1L] <- L$Ubegin[m + 1L] + nz_spike
  L$Uindex[L$Ubegin[m + 1L] + 1L] <- -1L
  L$Uvalue[L$Ubegin[m + 1L] + 1L] <- 0
  L$Ubegin[m + 1L] <- L$Ubegin[m + 1L] + 1L
  pos <- L$Ubegin[ipivot + 1L]
  while ((i <- L$Uindex[pos + 1L]) >= 0L) {
    j <- L$qmap[i + 1L]
    jnext <- L$Wflink[j + 1L]
    if (L$Wend[j + 1L] == L$Wbegin[jnext + 1L]) {
      nzj <- L$Wend[j + 1L] - L$Wbegin[j + 1L]
      room <- as.integer(1 + stretch * (nzj + 1) + pad)
      .blu_file_reappend(L, j, m, room)
    }
    end <- L$Wend[j + 1L]
    L$Wend[j + 1L] <- end + 1L
    L$Windex[end + 1L] <- jpivot
    L$Wvalue[end + 1L] <- L$Uvalue[pos + 1L]
    pos <- pos + 1L
  }
  Unz <- Unz + nz_spike
  L$col_pivot[jpivot + 1L] <- spike_diag
  L$row_pivot[ipivot + 1L] <- spike_diag
  row_reach <- NULL; col_reach <- NULL
  istriangular <- FALSE
  if (have_diag) {
    istriangular <- intersect == 0L
    if (istriangular) {
      L$min_pivot <- min(L$min_pivot, abs(newpiv))
      L$max_pivot <- max(L$max_pivot, abs(newpiv))
      eta_i <- if (re > rb) L$Lindex[(rb + 1L):re] else integer(0)
      row_reach <- c(ipivot, eta_i)
      col_reach <- c(jpivot, L$qmap[eta_i + 1L])
      L$nsymperm_total <- L$nsymperm_total + 1L
    }
  } else {
    path <- .blu_bfs_path(L, jpivot)
    istriangular <- TRUE
    L$marker <- L$marker + 1L; M <- L$marker
    reach <- integer(m); rtop <- m
    np <- length(path)
    t <- 1L
    while (t < np && istriangular) {
      j <- path[t]; jnext <- path[t + 1L]
      where <- .blu_find(jnext, L$Windex, L$Wbegin[j + 1L], L$Wend[j + 1L])
      L$Windex[where + 1L] <- j
      r <- .blu_dfs(L, j, L$Wbegin, L$Wend, L$Windex, rtop, reach, L$iwork0, M)
      rtop <- r$top; reach <- r$xi; L$iwork0 <- r$marked
      reach[rtop + 1L] <- jnext
      L$Windex[where + 1L] <- jnext
      istriangular <- L$iwork0[jnext + 1L] != M
      t <- t + 1L
    }
    if (istriangular) {
      j <- path[np]
      r <- .blu_dfs(L, j, L$Wbegin, L$Wend, L$Windex, rtop, reach, L$iwork0, M)
      rtop <- r$top; reach <- r$xi; L$iwork0 <- r$marked
      reach[rtop + 1L] <- jpivot
      L$iwork0[j + 1L] <- L$iwork0[j + 1L] - 1L
      pos <- L$Ubegin[ipivot + 1L]
      while ((i <- L$Uindex[pos + 1L]) >= 0L) {
        if (L$iwork0[L$qmap[i + 1L] + 1L] == M) istriangular <- FALSE
        pos <- pos + 1L
      }
      L$iwork0[j + 1L] <- L$iwork0[j + 1L] + 1L
    }
    if (istriangular) {
      .blu_permute(L, path)
      Unz <- Unz - 1L
      col_reach <- reach[(rtop + 1L):m]
      row_reach <- L$pmap[col_reach + 1L]
    }
  }
  if (!istriangular) {
    b <- L$Wbegin[jpivot + 1L]; e <- L$Wend[jpivot + 1L]
    if (e > b) for (pos in b:(e - 1L)) {
      j <- L$Windex[pos + 1L]
      where <- -1L
      end <- L$Ubegin[L$pmap[j + 1L] + 1L]
      while ((i <- L$Uindex[end + 1L]) >= 0L) {
        if (i == ipivot) where <- end
        end <- end + 1L
      }
      L$Uindex[where + 1L] <- L$Uindex[end]
      L$Uvalue[where + 1L] <- L$Uvalue[end]
      L$Uindex[end] <- -1L
      Unz <- Unz - 1L
    }
    L$Wend[jpivot + 1L] <- L$Wbegin[jpivot + 1L]
    L$col_pivot[jpivot + 1L] <- newpiv
    L$row_pivot[ipivot + 1L] <- newpiv
    L$min_pivot <- min(L$min_pivot, abs(newpiv))
    L$max_pivot <- max(L$max_pivot, abs(newpiv))
    nz <- 0L
    put <- rb
    max_eta <- 0
    if (re > rb) for (pos in rb:(re - 1L)) {
      if (L$Lvalue[pos + 1L] != 0) {
        max_eta <- max(max_eta, abs(L$Lvalue[pos + 1L]))
        L$Lindex[put + 1L] <- L$Lindex[pos + 1L]
        L$Lvalue[put + 1L] <- L$Lvalue[pos + 1L]
        put <- put + 1L
        nz <- nz + 1L
      }
    }
    L$Rbegin[nforrest + 2L] <- put
    L$Rnz <- L$Rnz + nz
    L$max_eta <- max(L$max_eta, max_eta)
    row_reach <- ipivot; col_reach <- jpivot
    L$nforrest <- L$nforrest + 1L
    L$nforrest_total <- L$nforrest_total + 1L
  }
  nreach <- length(row_reach)
  if (L$pivotlen + nreach > 2L * m) .blu_garbage_perm(L)
  L$pivotrow <- c(L$pivotrow[seq_len(L$pivotlen)], row_reach)
  L$pivotcol <- c(L$pivotcol[seq_len(L$pivotlen)], col_reach)
  L$pivotlen <- L$pivotlen + nreach
  used <- L$Ubegin[m + 1L]
  if (used - Unz - m > L$compress_thres * used) .blu_compress_packed(L)
  used <- L$Wbegin[m + 1L]
  need <- as.integer(Unz + stretch * Unz + m * pad)
  if ((used - need) > L$compress_thres * used) .blu_file_compress(L, m, stretch, pad)
  L$pivot_error <- piverr / (1 + abs(newpiv))
  L$Unz <- Unz
  L$btran_for_update <- -1L
  L$ftran_for_update <- -1L
  L$update_cost_numer <- L$update_cost_numer + nz_roweta
  L$nupdate <- L$nupdate + 1L
  L$nupdate_total <- L$nupdate_total + 1L
  "ok"
}

# basiclu_get_factors after a fresh factorization, with the unit diagonal
# of L removed as BasicLu::_GetFactors does. The factors are given in pivot
# step coordinates as lists of columns list(i = 1-based rows in increasing
# order, x = values); the diagonal of U is the last entry of each column.
# rowperm and colperm are the 1-based pivot rows and columns.
.blu_get_factors <- function(L) {
  m <- L$m
  p <- L$p
  pinvp <- integer(m); pinvp[p + 1L] <- seq_len(m)
  Lcols <- replicate(m, list(i = integer(0), x = numeric(0)), simplify = FALSE)
  for (k in 0:(m - 1L)) {
    pos <- L$Ltbegin_p[k + 1L]
    while ((i <- L$Lindex[pos + 1L]) >= 0L) {
      cl <- pinvp[i + 1L]
      Lcols[[cl]]$i <- c(Lcols[[cl]]$i, k + 1L)
      Lcols[[cl]]$x <- c(Lcols[[cl]]$x, L$Lvalue[pos + 1L])
      pos <- pos + 1L
    }
  }
  pivotcol <- L$pivotcol[seq_len(m)]
  kofj <- integer(m); kofj[pivotcol + 1L] <- seq_len(m)
  Ucols <- replicate(m, list(i = integer(0), x = numeric(0)), simplify = FALSE)
  for (k in 0:(m - 1L)) {
    j <- pivotcol[k + 1L]
    b <- L$Wbegin[j + 1L]; e <- L$Wend[j + 1L]
    if (e > b) for (pos in b:(e - 1L)) {
      cl <- kofj[L$Windex[pos + 1L] + 1L]
      Ucols[[cl]]$i <- c(Ucols[[cl]]$i, k + 1L)
      Ucols[[cl]]$x <- c(Ucols[[cl]]$x, L$Wvalue[pos + 1L])
    }
  }
  for (k in seq_len(m)) {
    Ucols[[k]]$i <- c(Ucols[[k]]$i, k)
    Ucols[[k]]$x <- c(Ucols[[k]]$x, L$col_pivot[pivotcol[k] + 1L])
  }
  list(L = Lcols, U = Ucols, rowperm = L$pivotrow[seq_len(m)] + 1L,
       colperm = pivotcol + 1L)
}

# The LuUpdate contract of IPX (lu_update.h) on BASICLU (BasicLu wrapper)
.ipx_lu_new <- function(m, pivottol) {
  lu <- new.env(parent = emptyenv())
  L <- .blu_new(m)
  L$reltol <- pivottol
  lu$L <- L
  lu$updates <- 0L
  lu$get_pivottol <- function() L$reltol
  lu$set_pivottol <- function(v) L$reltol <- v
  # B: list of columns list(i = 0-based rows, x = values); returns the flag
  # (1 unstable, 2 singular) and the dependent (position, row) pairs
  lu$factorize <- function(Bcols, strict) {
    lu$updates <- 0L
    L$abstol <- if (strict) 1e-3 else 1e-14
    st <- .blu_factorize(L, Bcols)
    flag <- 0L
    if (L$residual_test > 1e-12) flag <- bitwOr(flag, 1L)
    dep <- matrix(integer(0), 0, 2, dimnames = list(NULL, c("p", "i")))
    if (st == "singular") {
      flag <- bitwOr(flag, 2L)
      k <- (L$rank:(m - 1L)) + 1L
      dep <- cbind(p = L$pivotcol[k] + 1L, i = L$pivotrow[k] + 1L)
    }
    list(flag = flag, dependent = dep)
  }
  lu$solve_dense <- function(rhs, trans) .blu_solve_dense(L, rhs, trans)
  # ftran of a column given by 1-based rows and values; btran of position p.
  # lu$pattern is the 1-based pattern of the solution in BASICLU's order.
  # With want = FALSE only the data for the update are computed.
  lu$ftran <- function(rows, vals, want = TRUE) {
    x <- .blu_solve_for_update(L, rows - 1L, vals, "N", want)
    if (want) lu$pattern <- L$lhs_pattern + 1L
    x
  }
  lu$btran <- function(p, want = TRUE) {
    x <- .blu_solve_for_update(L, p - 1L, NULL, "T", want)
    if (want) lu$pattern <- L$lhs_pattern + 1L
    x
  }
  lu$update <- function(pivot) {
    lu$updates <- lu$updates + 1L
    if (.blu_update(L, pivot) == "singular") return(-1L)
    if (L$pivot_error > 1e-8) return(1L)
    0L
  }
  lu$get_factors <- function() .blu_get_factors(L)
  lu$need_fresh <- function() {
    L$nforrest == m || L$update_cost_numer / L$update_cost_denom > 1
  }
  lu
}
