# Sparse LDL factorization for the Clarabel interior point solver: the
# approximate minimum degree ordering of SuiteSparse AMD and the quasidefinite
# LDL factorization of QDLDL.jl.
#
# Ported from AMD 3.3.4 of SuiteSparse, Copyright (c) 1996-2022, Timothy A.
# Davis, Patrick R. Amestoy, and Iain S. Duff, BSD 3-clause license
# (inst/AMD-LICENSE), and from QDLDL.jl (https://github.com/osqp/QDLDL.jl),
# version 0.4.1, Copyright the QDLDL.jl authors, Apache License 2.0
# (inst/CLARABEL-LICENSE).
#
# The AMD functions keep the 0 based index values of the C code; R vectors are
# read and written at value + 1. Integer arithmetic is done in doubles, which
# hold every index exactly.

.AMD_EMPTY <- -1
.amd_flip <- function(i) -i - 2

.amd_valid <- function(n_row, n_col, Ap, Ai) {
  nz <- Ap[n_col + 1L]
  if (Ap[1L] != 0 || nz < 0) return("invalid")
  result <- "ok"
  for (j in seq_len(n_col) - 1L) {
    p1 <- Ap[j + 1L]; p2 <- Ap[j + 2L]
    if (p1 > p2) return("invalid")
    ilast <- .AMD_EMPTY
    p <- p1
    while (p < p2) {
      i <- Ai[p + 1L]
      if (i < 0 || i >= n_row) return("invalid")
      if (i <= ilast) result <- "jumbled"
      ilast <- i
      p <- p + 1
    }
  }
  result
}

# Sorted, duplicate free pattern of the transpose (AMD_preprocess).
.amd_preprocess <- function(n, Ap, Ai) {
  W <- numeric(n); Flag <- rep(.AMD_EMPTY, n)
  for (j in seq_len(n) - 1L) {
    p <- Ap[j + 1L]
    while (p < Ap[j + 2L]) {
      i <- Ai[p + 1L]
      if (Flag[i + 1L] != j) { W[i + 1L] <- W[i + 1L] + 1; Flag[i + 1L] <- j }
      p <- p + 1
    }
  }
  Rp <- c(0, cumsum(W))
  W <- Rp[seq_len(n)]
  Flag[] <- .AMD_EMPTY
  Ri <- numeric(Rp[n + 1L])
  for (j in seq_len(n) - 1L) {
    p <- Ap[j + 1L]
    while (p < Ap[j + 2L]) {
      i <- Ai[p + 1L]
      if (Flag[i + 1L] != j) {
        Ri[W[i + 1L] + 1L] <- j
        W[i + 1L] <- W[i + 1L] + 1
        Flag[i + 1L] <- j
      }
      p <- p + 1
    }
  }
  list(Rp = Rp, Ri = Ri)
}

# Column counts of A + A' without the diagonal (AMD_aat).
.amd_aat <- function(n, Ap, Ai) {
  Len <- numeric(n); Tp <- numeric(n)
  for (k in seq_len(n) - 1L) {
    p <- Ap[k + 1L]; p2 <- Ap[k + 2L]
    while (p < p2) {
      j <- Ai[p + 1L]
      if (j < k) {
        Len[j + 1L] <- Len[j + 1L] + 1; Len[k + 1L] <- Len[k + 1L] + 1
        p <- p + 1
      } else if (j == k) {
        p <- p + 1
        break
      } else break
      pj2 <- Ap[j + 2L]
      pj <- Tp[j + 1L]
      while (pj < pj2) {
        i <- Ai[pj + 1L]
        if (i < k) {
          Len[i + 1L] <- Len[i + 1L] + 1; Len[j + 1L] <- Len[j + 1L] + 1
          pj <- pj + 1
        } else if (i == k) {
          pj <- pj + 1
          break
        } else break
      }
      Tp[j + 1L] <- pj
    }
    Tp[k + 1L] <- p
  }
  for (j in seq_len(n) - 1L) {
    pj <- Tp[j + 1L]
    while (pj < Ap[j + 2L]) {
      i <- Ai[pj + 1L]
      Len[i + 1L] <- Len[i + 1L] + 1; Len[j + 1L] <- Len[j + 1L] + 1
      pj <- pj + 1
    }
  }
  list(Len = Len, nzaat = sum(Len))
}

.amd_post_tree <- function(root, k, Child, Sibling, Order, Stack) {
  head <- 0
  Stack[1L] <- root
  while (head >= 0) {
    i <- Stack[head + 1L]
    if (Child[i + 1L] != .AMD_EMPTY) {
      f <- Child[i + 1L]
      while (f != .AMD_EMPTY) { head <- head + 1; f <- Sibling[f + 1L] }
      h <- head
      f <- Child[i + 1L]
      while (f != .AMD_EMPTY) { Stack[h + 1L] <- f; h <- h - 1; f <- Sibling[f + 1L] }
      Child[i + 1L] <- .AMD_EMPTY
    } else {
      head <- head - 1
      Order[i + 1L] <- k
      k <- k + 1
    }
  }
  list(k = k, Child = Child, Order = Order, Stack = Stack)
}

.amd_postorder <- function(nn, Parent, Nv, Fsize) {
  Child <- rep(.AMD_EMPTY, nn); Sibling <- rep(.AMD_EMPTY, nn)
  for (j in rev(seq_len(nn) - 1L)) {
    if (Nv[j + 1L] > 0) {
      parent <- Parent[j + 1L]
      if (parent != .AMD_EMPTY) {
        Sibling[j + 1L] <- Child[parent + 1L]
        Child[parent + 1L] <- j
      }
    }
  }
  for (i in seq_len(nn) - 1L) {
    if (Nv[i + 1L] > 0 && Child[i + 1L] != .AMD_EMPTY) {
      fprev <- .AMD_EMPTY; maxfrsize <- .AMD_EMPTY; bigfprev <- .AMD_EMPTY; bigf <- .AMD_EMPTY
      f <- Child[i + 1L]
      while (f != .AMD_EMPTY) {
        frsize <- Fsize[f + 1L]
        if (frsize >= maxfrsize) { maxfrsize <- frsize; bigfprev <- fprev; bigf <- f }
        fprev <- f
        f <- Sibling[f + 1L]
      }
      fnext <- Sibling[bigf + 1L]
      if (fnext != .AMD_EMPTY) {
        if (bigfprev == .AMD_EMPTY) Child[i + 1L] <- fnext else Sibling[bigfprev + 1L] <- fnext
        Sibling[bigf + 1L] <- .AMD_EMPTY
        Sibling[fprev + 1L] <- bigf
      }
    }
  }
  Order <- rep(.AMD_EMPTY, nn); Stack <- numeric(nn)
  k <- 0
  for (i in seq_len(nn) - 1L) {
    if (Parent[i + 1L] == .AMD_EMPTY && Nv[i + 1L] > 0) {
      r <- .amd_post_tree(i, k, Child, Sibling, Order, Stack)
      k <- r$k; Child <- r$Child; Order <- r$Order; Stack <- r$Stack
    }
  }
  Order
}

# AMD_1 and AMD_2: build the pattern of A + A' and run the approximate minimum
# degree elimination. Returns the 0 based permutation P (pivot k is row P[k]).
.amd_1_2 <- function(n, Ap, Ai, Len, slen, alpha, aggressive) {
  iwlen <- slen - 6 * n
  Iw <- numeric(iwlen)
  Pe <- numeric(n); Sp <- numeric(n); Tp <- numeric(n)
  pfree <- 0
  for (j in seq_len(n) - 1L) {
    Pe[j + 1L] <- pfree; Sp[j + 1L] <- pfree
    pfree <- pfree + Len[j + 1L]
  }
  for (k in seq_len(n) - 1L) {
    p <- Ap[k + 1L]; p2 <- Ap[k + 2L]
    while (p < p2) {
      j <- Ai[p + 1L]
      if (j < k) {
        Iw[Sp[j + 1L] + 1L] <- k; Sp[j + 1L] <- Sp[j + 1L] + 1
        Iw[Sp[k + 1L] + 1L] <- j; Sp[k + 1L] <- Sp[k + 1L] + 1
        p <- p + 1
      } else if (j == k) {
        p <- p + 1
        break
      } else break
      pj2 <- Ap[j + 2L]
      pj <- Tp[j + 1L]
      while (pj < pj2) {
        i <- Ai[pj + 1L]
        if (i < k) {
          Iw[Sp[i + 1L] + 1L] <- j; Sp[i + 1L] <- Sp[i + 1L] + 1
          Iw[Sp[j + 1L] + 1L] <- i; Sp[j + 1L] <- Sp[j + 1L] + 1
          pj <- pj + 1
        } else if (i == k) {
          pj <- pj + 1
          break
        } else break
      }
      Tp[j + 1L] <- pj
    }
    Tp[k + 1L] <- p
  }
  for (j in seq_len(n) - 1L) {
    pj <- Tp[j + 1L]
    while (pj < Ap[j + 2L]) {
      i <- Ai[pj + 1L]
      Iw[Sp[i + 1L] + 1L] <- j; Sp[i + 1L] <- Sp[i + 1L] + 1
      Iw[Sp[j + 1L] + 1L] <- i; Sp[j + 1L] <- Sp[j + 1L] + 1
      pj <- pj + 1
    }
  }

  # AMD_2
  EMPTY <- .AMD_EMPTY
  FLIP <- .amd_flip
  me <- EMPTY; mindeg <- 0; nel <- 0; lemax <- 0
  dense <- if (alpha < 0) n - 2 else trunc(alpha * sqrt(n))
  dense <- min(n, max(16, dense))
  Last <- rep(EMPTY, n); Head <- rep(EMPTY, n); Next <- rep(EMPTY, n)
  Nv <- rep(1, n); W <- rep(1, n); Elen <- numeric(n); Degree <- Len
  # clear_flag: the counters never come near the integer limit of the C code
  clear_flag <- function(wflg) {
    if (wflg < 2) {
      W[W != 0] <<- 1
      wflg <- 2
    }
    wflg
  }
  wflg <- clear_flag(0)
  ndense <- 0
  for (i in seq_len(n) - 1L) {
    deg <- Degree[i + 1L]
    if (deg == 0) {
      Elen[i + 1L] <- FLIP(1); nel <- nel + 1; Pe[i + 1L] <- EMPTY; W[i + 1L] <- 0
    } else if (deg > dense) {
      ndense <- ndense + 1; Nv[i + 1L] <- 0; Elen[i + 1L] <- EMPTY; nel <- nel + 1; Pe[i + 1L] <- EMPTY
    } else {
      inext <- Head[deg + 1L]
      if (inext != EMPTY) Last[inext + 1L] <- i
      Next[i + 1L] <- inext
      Head[deg + 1L] <- i
    }
  }
  while (nel < n) {
    deg <- mindeg
    while (deg < n) {
      me <- Head[deg + 1L]
      if (me != EMPTY) break
      deg <- deg + 1
    }
    mindeg <- deg
    inext <- Next[me + 1L]
    if (inext != EMPTY) Last[inext + 1L] <- EMPTY
    Head[deg + 1L] <- inext
    elenme <- Elen[me + 1L]
    nvpiv <- Nv[me + 1L]
    nel <- nel + nvpiv
    Nv[me + 1L] <- -nvpiv
    degme <- 0
    if (elenme == 0) {
      pme1 <- Pe[me + 1L]
      pme2 <- pme1 - 1
      p <- pme1
      while (p <= pme1 + Len[me + 1L] - 1) {
        i <- Iw[p + 1L]
        nvi <- Nv[i + 1L]
        if (nvi > 0) {
          degme <- degme + nvi
          Nv[i + 1L] <- -nvi
          pme2 <- pme2 + 1
          Iw[pme2 + 1L] <- i
          ilast <- Last[i + 1L]; inext <- Next[i + 1L]
          if (inext != EMPTY) Last[inext + 1L] <- ilast
          if (ilast != EMPTY) Next[ilast + 1L] <- inext else Head[Degree[i + 1L] + 1L] <- inext
        }
        p <- p + 1
      }
    } else {
      p <- Pe[me + 1L]
      pme1 <- pfree
      slenme <- Len[me + 1L] - elenme
      for (knt1 in seq_len(elenme + 1)) {
        if (knt1 > elenme) {
          e <- me; pj <- p; ln <- slenme
        } else {
          e <- Iw[p + 1L]; p <- p + 1
          pj <- Pe[e + 1L]; ln <- Len[e + 1L]
        }
        knt2 <- 1
        while (knt2 <= ln) {
          i <- Iw[pj + 1L]; pj <- pj + 1
          nvi <- Nv[i + 1L]
          if (nvi > 0) {
            if (pfree >= iwlen) {
              # compress Iw
              Pe[me + 1L] <- p
              Len[me + 1L] <- Len[me + 1L] - knt1
              if (Len[me + 1L] == 0) Pe[me + 1L] <- EMPTY
              Pe[e + 1L] <- pj
              Len[e + 1L] <- ln - knt2
              if (Len[e + 1L] == 0) Pe[e + 1L] <- EMPTY
              for (j in seq_len(n) - 1L) {
                pn <- Pe[j + 1L]
                if (pn >= 0) {
                  Pe[j + 1L] <- Iw[pn + 1L]
                  Iw[pn + 1L] <- FLIP(j)
                }
              }
              psrc <- 0; pdst <- 0; pend <- pme1 - 1
              while (psrc <= pend) {
                j <- FLIP(Iw[psrc + 1L]); psrc <- psrc + 1
                if (j >= 0) {
                  Iw[pdst + 1L] <- Pe[j + 1L]
                  Pe[j + 1L] <- pdst; pdst <- pdst + 1
                  lenj <- Len[j + 1L]
                  knt3 <- 0
                  while (knt3 <= lenj - 2) {
                    Iw[pdst + 1L] <- Iw[psrc + 1L]; pdst <- pdst + 1; psrc <- psrc + 1
                    knt3 <- knt3 + 1
                  }
                }
              }
              p1 <- pdst
              psrc <- pme1
              while (psrc <= pfree - 1) {
                Iw[pdst + 1L] <- Iw[psrc + 1L]; pdst <- pdst + 1
                psrc <- psrc + 1
              }
              pme1 <- p1
              pfree <- pdst
              pj <- Pe[e + 1L]
              p <- Pe[me + 1L]
            }
            degme <- degme + nvi
            Nv[i + 1L] <- -nvi
            Iw[pfree + 1L] <- i; pfree <- pfree + 1
            ilast <- Last[i + 1L]; inext <- Next[i + 1L]
            if (inext != EMPTY) Last[inext + 1L] <- ilast
            if (ilast != EMPTY) Next[ilast + 1L] <- inext else Head[Degree[i + 1L] + 1L] <- inext
          }
          knt2 <- knt2 + 1
        }
        if (e != me) {
          Pe[e + 1L] <- FLIP(me)
          W[e + 1L] <- 0
        }
      }
      pme2 <- pfree - 1
    }
    Degree[me + 1L] <- degme
    Pe[me + 1L] <- pme1
    Len[me + 1L] <- pme2 - pme1 + 1
    Elen[me + 1L] <- FLIP(nvpiv + degme)
    wflg <- clear_flag(wflg)
    pme <- pme1
    while (pme <= pme2) {
      i <- Iw[pme + 1L]
      eln <- Elen[i + 1L]
      if (eln > 0) {
        nvi <- -Nv[i + 1L]
        wnvi <- wflg - nvi
        p <- Pe[i + 1L]
        while (p <= Pe[i + 1L] + eln - 1) {
          e <- Iw[p + 1L]
          we <- W[e + 1L]
          if (we >= wflg) we <- we - nvi else if (we != 0) we <- Degree[e + 1L] + wnvi
          W[e + 1L] <- we
          p <- p + 1
        }
      }
      pme <- pme + 1
    }
    pme <- pme1
    while (pme <= pme2) {
      i <- Iw[pme + 1L]
      p1 <- Pe[i + 1L]
      p2 <- p1 + Elen[i + 1L] - 1
      pn <- p1
      hash <- 0
      deg <- 0
      p <- p1
      while (p <= p2) {
        e <- Iw[p + 1L]
        we <- W[e + 1L]
        if (we != 0) {
          dext <- we - wflg
          if (!aggressive || dext > 0) {
            deg <- deg + dext
            Iw[pn + 1L] <- e; pn <- pn + 1
            hash <- hash + e
          } else {
            Pe[e + 1L] <- FLIP(me)
            W[e + 1L] <- 0
          }
        }
        p <- p + 1
      }
      Elen[i + 1L] <- pn - p1 + 1
      p3 <- pn
      p4 <- p1 + Len[i + 1L]
      p <- p2 + 1
      while (p < p4) {
        j <- Iw[p + 1L]
        nvj <- Nv[j + 1L]
        if (nvj > 0) {
          deg <- deg + nvj
          Iw[pn + 1L] <- j; pn <- pn + 1
          hash <- hash + j
        }
        p <- p + 1
      }
      if (Elen[i + 1L] == 1 && p3 == pn) {
        Pe[i + 1L] <- FLIP(me)
        nvi <- -Nv[i + 1L]
        degme <- degme - nvi
        nvpiv <- nvpiv + nvi
        nel <- nel + nvi
        Nv[i + 1L] <- 0
        Elen[i + 1L] <- EMPTY
      } else {
        Degree[i + 1L] <- min(Degree[i + 1L], deg)
        Iw[pn + 1L] <- Iw[p3 + 1L]
        Iw[p3 + 1L] <- Iw[p1 + 1L]
        Iw[p1 + 1L] <- me
        Len[i + 1L] <- pn - p1 + 1
        hash <- hash %% n
        j <- Head[hash + 1L]
        if (j <= EMPTY) {
          Next[i + 1L] <- FLIP(j)
          Head[hash + 1L] <- FLIP(i)
        } else {
          Next[i + 1L] <- Last[j + 1L]
          Last[j + 1L] <- i
        }
        Last[i + 1L] <- hash
      }
      pme <- pme + 1
    }
    Degree[me + 1L] <- degme
    lemax <- max(lemax, degme)
    wflg <- wflg + lemax
    wflg <- clear_flag(wflg)
    # supervariable detection
    pme <- pme1
    while (pme <= pme2) {
      i <- Iw[pme + 1L]
      if (Nv[i + 1L] < 0) {
        hash <- Last[i + 1L]
        j <- Head[hash + 1L]
        if (j == EMPTY) {
          i <- EMPTY
        } else if (j < EMPTY) {
          i <- FLIP(j)
          Head[hash + 1L] <- EMPTY
        } else {
          i <- Last[j + 1L]
          Last[j + 1L] <- EMPTY
        }
        while (i != EMPTY && Next[i + 1L] != EMPTY) {
          ln <- Len[i + 1L]
          eln <- Elen[i + 1L]
          p <- Pe[i + 1L] + 1
          while (p <= Pe[i + 1L] + ln - 1) {
            W[Iw[p + 1L] + 1L] <- wflg
            p <- p + 1
          }
          jlast <- i
          j <- Next[i + 1L]
          while (j != EMPTY) {
            ok <- (Len[j + 1L] == ln) && (Elen[j + 1L] == eln)
            p <- Pe[j + 1L] + 1
            while (ok && p <= Pe[j + 1L] + ln - 1) {
              if (W[Iw[p + 1L] + 1L] != wflg) ok <- FALSE
              p <- p + 1
            }
            if (ok) {
              Pe[j + 1L] <- FLIP(i)
              Nv[i + 1L] <- Nv[i + 1L] + Nv[j + 1L]
              Nv[j + 1L] <- 0
              Elen[j + 1L] <- EMPTY
              j <- Next[j + 1L]
              Next[jlast + 1L] <- j
            } else {
              jlast <- j
              j <- Next[j + 1L]
            }
          }
          wflg <- wflg + 1
          i <- Next[i + 1L]
        }
      }
      pme <- pme + 1
    }
    # finalize the new element
    p <- pme1
    nleft <- n - nel
    pme <- pme1
    while (pme <= pme2) {
      i <- Iw[pme + 1L]
      nvi <- -Nv[i + 1L]
      if (nvi > 0) {
        Nv[i + 1L] <- nvi
        deg <- Degree[i + 1L] + degme - nvi
        deg <- min(deg, nleft - nvi)
        inext <- Head[deg + 1L]
        if (inext != EMPTY) Last[inext + 1L] <- i
        Next[i + 1L] <- inext
        Last[i + 1L] <- EMPTY
        Head[deg + 1L] <- i
        mindeg <- min(mindeg, deg)
        Degree[i + 1L] <- deg
        Iw[p + 1L] <- i; p <- p + 1
      }
      pme <- pme + 1
    }
    Nv[me + 1L] <- nvpiv
    Len[me + 1L] <- p - pme1
    if (Len[me + 1L] == 0) {
      Pe[me + 1L] <- EMPTY
      W[me + 1L] <- 0
    }
    if (elenme != 0) pfree <- p
  }
  # postordering
  Pe <- FLIP(Pe)
  Elen <- FLIP(Elen)
  for (i in seq_len(n) - 1L) {
    if (Nv[i + 1L] == 0) {
      j <- Pe[i + 1L]
      if (j == EMPTY) next
      while (Nv[j + 1L] == 0) j <- Pe[j + 1L]
      e <- j
      j <- i
      while (Nv[j + 1L] == 0) {
        jnext <- Pe[j + 1L]
        Pe[j + 1L] <- e
        j <- jnext
      }
    }
  }
  W <- .amd_postorder(n, Pe, Nv, Elen)
  Head <- rep(EMPTY, n); Next <- rep(EMPTY, n)
  for (e in seq_len(n) - 1L) {
    k <- W[e + 1L]
    if (k != EMPTY) Head[k + 1L] <- e
  }
  nel <- 0
  for (k in seq_len(n) - 1L) {
    e <- Head[k + 1L]
    if (e == EMPTY) break
    Next[e + 1L] <- nel
    nel <- nel + Nv[e + 1L]
  }
  for (i in seq_len(n) - 1L) {
    if (Nv[i + 1L] == 0) {
      e <- Pe[i + 1L]
      if (e != EMPTY) {
        Next[i + 1L] <- Next[e + 1L]
        Next[e + 1L] <- Next[e + 1L] + 1
      } else {
        Next[i + 1L] <- nel
        nel <- nel + 1
      }
    }
  }
  Last <- numeric(n)
  for (i in seq_len(n) - 1L) Last[Next[i + 1L] + 1L] <- i
  Last
}

# AMD ordering of a square sparse pattern given by 1 based CSC arrays, as
# AMD.jl amd(A, meta): returns a 1 based permutation. `dense` is the AMD_DENSE
# control (default 10), `aggressive` the AMD_AGGRESSIVE control.
.amd_order <- function(n, colptr, rowval, dense = 10, aggressive = TRUE) {
  if (n == 0) return(integer())
  Ap <- colptr - 1; Ai <- rowval - 1
  status <- .amd_valid(n, n, Ap, Ai)
  if (status == "invalid") stop("AMD: invalid input matrix", call. = FALSE)
  if (status == "jumbled") {
    pre <- .amd_preprocess(n, Ap, Ai)
    Ap <- pre$Rp; Ai <- pre$Ri
  }
  aat <- .amd_aat(n, Ap, Ai)
  nzaat <- aat$nzaat
  slen <- nzaat + nzaat %/% 5 + 7 * n
  as.integer(.amd_1_2(n, Ap, Ai, aat$Len, slen, dense, aggressive) + 1)
}

# ---------------------------------------------------------------------------
# QDLDL.jl

# A sparse matrix in compressed sparse column form, as SparseMatrixCSC.
.csc <- function(m, n, colptr, rowval, nzval) list(m = m, n = n, colptr = colptr, rowval = rowval, nzval = nzval)

# Symmetric permutation P A P' of an upper triangular CSC matrix, keeping the
# upper triangle; AtoPAPt maps each entry of A to its place in the result.
.qdldl_permute_symmetric <- function(A, iperm) {
  n <- A$n
  Ar <- A$rowval; Ac <- A$colptr; Av <- A$nzval
  num_entries <- integer(n)
  for (colA in seq_len(n)) {
    colP <- iperm[colA]
    for (ri in seq_len(Ac[colA + 1L] - Ac[colA]) + Ac[colA] - 1L) {
      rowA <- Ar[ri]
      if (rowA <= colA) {
        ci <- max(iperm[rowA], colP)
        num_entries[ci] <- num_entries[ci] + 1L
      }
    }
  }
  Pc <- c(1L, 1L + cumsum(num_entries))
  row_starts <- Pc[seq_len(n)]
  nzA <- length(Av)
  Pr <- integer(nzA); Pv <- numeric(nzA); AtoPAPt <- integer(nzA)
  for (colA in seq_len(n)) {
    colP <- iperm[colA]
    for (ri in seq_len(Ac[colA + 1L] - Ac[colA]) + Ac[colA] - 1L) {
      rowA <- Ar[ri]
      if (rowA <= colA) {
        rowP <- iperm[rowA]
        ci <- max(colP, rowP)
        k <- row_starts[ci]
        Pr[k] <- min(colP, rowP)
        Pv[k] <- Av[ri]
        AtoPAPt[ri] <- k
        row_starts[ci] <- k + 1L
      }
    }
  }
  nz_new <- Pc[n + 1L] - 1L
  list(P = .csc(n, n, Pc, Pr[seq_len(nz_new)], Pv[seq_len(nz_new)]), AtoPAPt = AtoPAPt)
}

.qdldl_etree <- function(n, Ap, Ai) {
  work <- integer(n); Lnz <- integer(n); etree <- rep(-1L, n)
  for (i in seq_len(n)) if (Ap[i] == Ap[i + 1L]) return(NULL)
  for (j in seq_len(n)) {
    work[j] <- j
    for (p in seq_len(Ap[j + 1L] - Ap[j]) + Ap[j] - 1L) {
      i <- Ai[p]
      if (i > j) return(NULL)
      while (work[i] != j) {
        if (etree[i] == -1L) etree[i] <- j
        Lnz[i] <- Lnz[i] + 1L
        work[i] <- j
        i <- etree[i]
      }
    }
  }
  list(etree = etree, Lnz = Lnz)
}

# The numeric (or logical) factorization QDLDL_factor! of workspace `ws`, an
# environment. Returns the number of positive entries of D, or -1.
.qdldl_factor <- function(ws, logical) {
  A <- ws$triuA
  n <- A$n; Ap <- A$colptr; Ai <- A$rowval; Ax <- A$nzval
  Lnz <- ws$Lnz; etree <- ws$etree; Dsigns <- ws$Dsigns
  eps <- ws$regularize_eps; delta <- ws$regularize_delta
  Lp <- ws$Lp; Li <- ws$Li
  Lx <- if (logical) rep(1, length(ws$Lx)) else ws$Lx
  # a logical factorization starts from unit L and D; the numeric one from zero D
  D <- if (logical) rep(1, n) else numeric(n)
  Dinv <- if (logical) rep(1, n) else ws$Dinv
  positive <- 0L
  regcount <- 0L
  yMarkers <- logical(n); yIdx <- integer(n); elimBuffer <- integer(n)
  yVals <- numeric(n)
  Lp[1L] <- 1L
  for (i in seq_len(n)) Lp[i + 1L] <- Lp[i] + Lnz[i]
  LNext <- Lp[seq_len(n)]
  if (!logical) {
    D[1L] <- Ax[1L]
    if (!is.null(Dsigns) && Dsigns[1L] * D[1L] < eps) {
      D[1L] <- delta * Dsigns[1L]
      regcount <- regcount + 1L
    }
    if (D[1L] == 0) return(-1L)
    if (D[1L] > 0) positive <- positive + 1L
    Dinv[1L] <- 1 / D[1L]
  }
  if (n >= 2L) for (k in 2:n) {
    nnzY <- 0L
    for (ii in seq_len(Ap[k + 1L] - Ap[k]) + Ap[k] - 1L) {
      bidx <- Ai[ii]
      if (bidx == k) {
        D[k] <- Ax[ii]
        next
      }
      yVals[bidx] <- Ax[ii]
      nextIdx <- bidx
      if (!yMarkers[nextIdx]) {
        yMarkers[nextIdx] <- TRUE
        elimBuffer[1L] <- nextIdx
        nnzE <- 1L
        nextIdx <- etree[bidx]
        while (nextIdx != -1L && nextIdx < k) {
          if (yMarkers[nextIdx]) break
          yMarkers[nextIdx] <- TRUE
          nnzE <- nnzE + 1L
          elimBuffer[nnzE] <- nextIdx
          nextIdx <- etree[nextIdx]
        }
        while (nnzE != 0L) {
          nnzY <- nnzY + 1L
          yIdx[nnzY] <- elimBuffer[nnzE]
          nnzE <- nnzE - 1L
        }
      }
    }
    for (i in rev(seq_len(nnzY))) {
      cidx <- yIdx[i]
      tmpIdx <- LNext[cidx]
      if (!logical) {
        yc <- yVals[cidx]
        for (j in seq_len(tmpIdx - Lp[cidx]) + Lp[cidx] - 1L) yVals[Li[j]] <- yVals[Li[j]] - Lx[j] * yc
        Lx[tmpIdx] <- yc * Dinv[cidx]
        D[k] <- D[k] - yc * Lx[tmpIdx]
      }
      Li[tmpIdx] <- k
      LNext[cidx] <- LNext[cidx] + 1L
      yVals[cidx] <- 0
      yMarkers[cidx] <- FALSE
    }
    if (!is.null(Dsigns) && Dsigns[k] * D[k] < eps) {
      D[k] <- delta * Dsigns[k]
      regcount <- regcount + 1L
    }
    if (D[k] == 0) return(-1L)
    if (D[k] > 0) positive <- positive + 1L
    Dinv[k] <- 1 / D[k]
  }
  ws$Lp <- Lp; ws$Li <- Li; ws$Lx <- Lx; ws$D <- D; ws$Dinv <- Dinv
  ws$regularize_count <- regcount
  positive
}

# qdldl(A; amd_dense_scale, logical, Dsigns, regularize_eps, regularize_delta)
# for an upper triangular CSC matrix A. Returns an environment.
.qdldl <- function(A, amd_dense_scale = 1, logical = FALSE, Dsigns = NULL, regularize_eps = 1e-12,
                   regularize_delta = 1e-7) {
  perm <- .amd_order(A$n, A$colptr, A$rowval, dense = 10 * amd_dense_scale)
  iperm <- integer(length(perm)); iperm[perm] <- seq_along(perm)
  pa <- .qdldl_permute_symmetric(A, iperm)
  ws <- new.env(parent = emptyenv())
  ws$perm <- perm; ws$iperm <- iperm
  ws$triuA <- pa$P; ws$AtoPAPt <- pa$AtoPAPt
  ws$Dsigns <- if (is.null(Dsigns)) NULL else Dsigns[perm]
  ws$regularize_eps <- regularize_eps; ws$regularize_delta <- regularize_delta
  et <- .qdldl_etree(A$n, pa$P$colptr, pa$P$rowval)
  if (is.null(et)) stop("Input matrix is not upper triangular or has an empty column", call. = FALSE)
  ws$etree <- et$etree; ws$Lnz <- et$Lnz
  sumLnz <- sum(et$Lnz)
  ws$Lp <- integer(A$n + 1L); ws$Li <- integer(sumLnz); ws$Lx <- numeric(sumLnz)
  ws$D <- numeric(A$n); ws$Dinv <- numeric(A$n)
  ws$logical <- logical
  .qdldl_refactor(ws, logical)
  ws
}

.qdldl_refactor <- function(ws, logical = FALSE) {
  ws$logical <- logical
  pos <- .qdldl_factor(ws, logical)
  if (pos < 0L) stop("Zero entry in D (matrix is not quasidefinite)", call. = FALSE)
  ws$positive_inertia <- pos
  invisible(ws)
}

# update_values!: new values for entries `idx` of the unpermuted matrix
.qdldl_update_values <- function(ws, idx, values) {
  ws$triuA$nzval[ws$AtoPAPt[idx]] <- values
  invisible(ws)
}

.qdldl_solve <- function(ws, b) {
  n <- ws$triuA$n
  Lp <- ws$Lp; Li <- ws$Li; Lx <- ws$Lx
  x <- b[ws$perm]
  for (i in seq_len(n)) {
    for (j in seq_len(Lp[i + 1L] - Lp[i]) + Lp[i] - 1L) x[Li[j]] <- x[Li[j]] - Lx[j] * x[i]
  }
  x <- x * ws$Dinv
  for (i in rev(seq_len(n))) {
    for (j in seq_len(Lp[i + 1L] - Lp[i]) + Lp[i] - 1L) x[i] <- x[i] - Lx[j] * x[Li[j]]
  }
  out <- numeric(n)
  out[ws$perm] <- x
  out
}
