# The clique table and implications of the HiGHS 1.12.0 MIP solver
# (mip/HighsCliqueTable.cpp, mip/HighsImplications.cpp) and the hash tree
# they store their incidence lists in (util/HighsHashTree.h), as SciPy
# 1.18.1 vendors them. MIT license notice in
# R/optimization-scipy-highs-util.R. Every stored index is 0-based exactly
# as in C++ (columns, clique ids, CliqueVar index 2 * col + val), and -1
# keeps its C++ meaning.

# ---------------------------------------------------------------------------
# HighsHashTree<HighsInt, V>
#
# A tree is NULL (kEmpty) or a node list. Node types are the C++ Type enum:
# 1 ListLeaf list(t, k, v); 2..5 InnerLeaf of size class t - 1
# list(t, occ, h, k, v) with occ a 64 long logical (bit b at [b + 1]) and h
# the 16 bit hash chunks; 6 BranchNode list(t, occ, child). Entries keep the
# exact C++ order, so .hmip_hht_entries() gives the for_each order.
# ---------------------------------------------------------------------------

.hmip_hht_cache <- new.env(parent = emptyenv())

# HighsHashHelpers::hash(HighsInt key) and its chunks: c6[p + 1] is
# get_hash_chunk(hash, p), c16[p + 1] is get_hash_chunks16(hash, p)
.hmip_hht_hash <- function(key) {
  nm <- as.character(key)
  h <- .hmip_hht_cache[[nm]]
  if (!is.null(h)) return(h)
  a <- if (key < 0) key + 4294967296 else key
  x <- .hs_u64_xor(.hs_pair_hash(1L, a, 0), .hs_u64_shift(.hs_pair_hash(0L, a, 0), 32L))
  bits <- .hs_limbs_bits(x)
  val <- function(lo, hi) sum(bits[(lo + 1L):(hi + 1L)] * 2^(0:(hi - lo)))
  c6 <- vapply(0:8, function(p) val(58L - 6L * p, 63L - 6L * p), 0)
  c16 <- vapply(0:8, function(p) {
    lo <- 48L - 6L * p
    val(lo, min(63L, lo + 15L))
  }, 0)
  h <- list(c6 = as.integer(c6), c16 = as.integer(c16))
  assign(nm, h, envir = .hmip_hht_cache)
  h
}

.hmip_hht_cap <- c(6L, 22L, 38L, 54L)
.hmip_hht_size_class <- function(n) 1L + (n + 9L) %/% 16L
# Occupation::num_set_until(pos): set bits at positions >= pos
.hmip_hht_nsu <- function(occ, pos) sum(occ[(pos + 1L):64L])
.hmip_hht_first <- function(h16) h16 %/% 1024L

.hmip_hht_newleaf <- function(cls) {
  list(t = cls + 1L, occ = logical(64L), h = integer(), k = integer(), v = list())
}

# hashes[i] with the zero sentinel at and after size
.hmip_hht_H <- function(leaf, i) if (i < length(leaf$k)) leaf$h[i + 1L] else 0L

.hmip_hht_find_key <- function(leaf, key, hash, pos) {
  size <- length(leaf$k)
  while (pos != size && leaf$h[pos + 1L] == hash) {
    if (key == leaf$k[pos + 1L]) return(list(found = TRUE, pos = pos))
    pos <- pos + 1L
  }
  list(found = FALSE, pos = pos)
}

# InnerLeaf::insert_entry; returns list(node, inserted)
.hmip_hht_leaf_insert <- function(leaf, hashPos, key, val) {
  hash <- .hmip_hht_hash(key)$c16[hashPos + 1L]
  chunk <- .hmip_hht_first(hash)
  pos <- .hmip_hht_nsu(leaf$occ, chunk)
  size <- length(leaf$k)
  if (leaf$occ[chunk + 1L]) {
    pos <- pos - 1L
    while (.hmip_hht_H(leaf, pos) > hash) pos <- pos + 1L
    fk <- .hmip_hht_find_key(leaf, key, hash, pos)
    if (fk$found) return(list(node = leaf, inserted = FALSE))
    pos <- fk$pos
  } else {
    leaf$occ[chunk + 1L] <- TRUE
    if (pos < size) while (.hmip_hht_H(leaf, pos) > hash) pos <- pos + 1L
  }
  leaf$h <- append(leaf$h, hash, after = pos)
  leaf$k <- append(leaf$k, as.integer(key), after = pos)
  leaf$v <- append(leaf$v, list(val), after = pos)
  list(node = leaf, inserted = TRUE)
}

.hmip_hht_leaf_find <- function(leaf, hashPos, key) {
  hash <- .hmip_hht_hash(key)$c16[hashPos + 1L]
  chunk <- .hmip_hht_first(hash)
  if (!leaf$occ[chunk + 1L]) return(NULL)
  pos <- .hmip_hht_nsu(leaf$occ, chunk) - 1L
  while (.hmip_hht_H(leaf, pos) > hash) pos <- pos + 1L
  fk <- .hmip_hht_find_key(leaf, key, hash, pos)
  if (fk$found) fk$pos else NULL
}

# InnerLeaf::erase_entry; returns the leaf, or NULL when key is absent
.hmip_hht_leaf_erase <- function(leaf, hashPos, key) {
  hash <- .hmip_hht_hash(key)$c16[hashPos + 1L]
  chunk <- .hmip_hht_first(hash)
  if (!leaf$occ[chunk + 1L]) return(NULL)
  startPos <- .hmip_hht_nsu(leaf$occ, chunk) - 1L
  while (.hmip_hht_first(.hmip_hht_H(leaf, startPos)) > chunk) startPos <- startPos + 1L
  pos <- startPos
  while (.hmip_hht_H(leaf, pos) > hash) pos <- pos + 1L
  fk <- .hmip_hht_find_key(leaf, key, hash, pos)
  if (!fk$found) return(NULL)
  pos <- fk$pos
  leaf$h <- leaf$h[-(pos + 1L)]
  leaf$k <- leaf$k[-(pos + 1L)]
  leaf$v <- leaf$v[-(pos + 1L)]
  size <- length(leaf$k)
  if (pos < size) {
    if (.hmip_hht_first(.hmip_hht_H(leaf, startPos)) != chunk) leaf$occ[chunk + 1L] <- !leaf$occ[chunk + 1L]
  } else if (startPos == pos) {
    leaf$occ[chunk + 1L] <- !leaf$occ[chunk + 1L]
  }
  leaf
}

# InnerLeaf::rehash
.hmip_hht_leaf_rehash <- function(leaf, hashPos) {
  size <- length(leaf$k)
  leaf$occ <- logical(64L)
  for (i in seq_len(size)) {
    leaf$h[i] <- .hmip_hht_hash(leaf$k[i])$c16[hashPos + 1L]
    leaf$occ[.hmip_hht_first(leaf$h[i]) + 1L] <- TRUE
  }
  h <- leaf$h; k <- leaf$k; v <- leaf$v
  i <- 0L
  while (i < size) {
    chunk <- .hmip_hht_first(h[i + 1L])
    pos <- .hmip_hht_nsu(leaf$occ, chunk) - 1L
    if (pos > i) {
      t <- h[pos + 1L]; h[pos + 1L] <- h[i + 1L]; h[i + 1L] <- t
      t <- k[pos + 1L]; k[pos + 1L] <- k[i + 1L]; k[i + 1L] <- t
      t <- v[pos + 1L]; v[pos + 1L] <- v[i + 1L]; v[i + 1L] <- t
      next
    }
    while (pos < i && h[pos + 1L] >= h[i + 1L]) pos <- pos + 1L
    if (pos < i) {
      hh <- h[i + 1L]; kk <- k[i + 1L]; vv <- v[i + 1L]
      idx <- (pos + 1L):i
      h[idx + 1L] <- h[idx]; k[idx + 1L] <- k[idx]; v[idx + 1L] <- v[idx]
      h[pos + 1L] <- hh; k[pos + 1L] <- kk; v[pos + 1L] <- vv
    }
    i <- i + 1L
  }
  leaf$h <- h; leaf$k <- k; leaf$v <- v
  leaf
}

.hmip_hht_type <- function(node) if (is.null(node)) 0L else node$t

.hmip_hht_num_estimate <- function(node) {
  t <- .hmip_hht_type(node)
  if (t == 0L) 0L else if (t == 1L) 1L else if (t == 6L) 64L else .hmip_hht_cap[t - 1L]
}
.hmip_hht_num_entries <- function(node) {
  t <- .hmip_hht_type(node)
  if (t == 0L) 0L else if (t == 6L) 64L else length(node$k)
}

# insert_recurse; returns list(node, inserted)
.hmip_hht_insert_recurse <- function(node, hashPos, key, val) {
  hh <- .hmip_hht_hash(key)
  t <- .hmip_hht_type(node)
  if (t == 0L) {
    if (hashPos == 9L) return(list(node = list(t = 1L, k = as.integer(key), v = list(val)), inserted = TRUE))
    return(.hmip_hht_leaf_insert(.hmip_hht_newleaf(1L), hashPos, key, val))
  }
  if (t == 1L) {
    if (any(node$k == key)) return(list(node = node, inserted = FALSE))
    node$k <- c(node$k, as.integer(key)); node$v <- c(node$v, list(val))
    return(list(node = node, inserted = TRUE))
  }
  if (t >= 2L && t <= 4L) {
    if (length(node$k) == .hmip_hht_cap[t - 1L]) {
      if (!is.null(.hmip_hht_leaf_find(node, hashPos, key))) return(list(node = node, inserted = FALSE))
      node$t <- t + 1L
    }
    return(.hmip_hht_leaf_insert(node, hashPos, key, val))
  }
  if (t == 5L) {
    leaf <- node
    if (length(leaf$k) < 54L) return(.hmip_hht_leaf_insert(leaf, hashPos, key, val))
    if (!is.null(.hmip_hht_leaf_find(leaf, hashPos, key))) return(list(node = leaf, inserted = FALSE))
    occ <- leaf$occ
    hashChunk <- hh$c6[hashPos + 1L]
    occ[hashChunk + 1L] <- TRUE
    branchSize <- sum(occ)
    branch <- list(t = 6L, occ = occ, child = vector("list", branchSize))
    lsize <- length(leaf$k)
    if (hashPos + 1L == 9L) {
      addlist <- function(branch, pos, k, v) {
        ch <- branch$child[[pos + 1L]]
        if (is.null(ch)) ch <- list(t = 1L, k = as.integer(k), v = list(v))
        else { ch$k <- c(as.integer(k), ch$k); ch$v <- c(list(v), ch$v) }
        branch$child[pos + 1L] <- list(ch)
        branch
      }
      for (i in seq_len(lsize)) {
        pos <- .hmip_hht_nsu(occ, .hmip_hht_first(leaf$h[i])) - 1L
        branch <- addlist(branch, pos, leaf$k[i], leaf$v[[i]])
      }
      pos <- .hmip_hht_nsu(occ, hashChunk) - 1L
      branch <- addlist(branch, pos, key, val)
      return(list(node = branch, inserted = TRUE))
    }
    if (branchSize > 1L) {
      maxEntriesPerLeaf <- 2L + lsize - branchSize
      if (maxEntriesPerLeaf <= 6L) {
        for (i in seq_len(branchSize)) branch$child[i] <- list(.hmip_hht_newleaf(1L))
        for (i in seq_len(lsize)) {
          pos <- .hmip_hht_nsu(occ, .hmip_hht_first(leaf$h[i])) - 1L
          branch$child[pos + 1L] <- list(.hmip_hht_leaf_insert(branch$child[[pos + 1L]], hashPos + 1L, leaf$k[i], leaf$v[[i]])$node)
        }
        pos <- .hmip_hht_nsu(occ, hashChunk) - 1L
        r <- .hmip_hht_leaf_insert(branch$child[[pos + 1L]], hashPos + 1L, key, val)
        branch$child[pos + 1L] <- list(r$node)
        return(list(node = branch, inserted = r$inserted))
      }
      sizes <- integer(branchSize)
      p0 <- .hmip_hht_nsu(occ, hashChunk)
      sizes[p0] <- sizes[p0] + 1L
      for (i in seq_len(lsize)) {
        pos <- .hmip_hht_nsu(occ, .hmip_hht_first(leaf$h[i]))
        sizes[pos] <- sizes[pos] + 1L
      }
      for (i in seq_len(branchSize)) branch$child[i] <- list(.hmip_hht_newleaf(.hmip_hht_size_class(sizes[i])))
      for (i in seq_len(lsize)) {
        pos <- .hmip_hht_nsu(occ, .hmip_hht_first(leaf$h[i])) - 1L
        branch$child[pos + 1L] <- list(.hmip_hht_leaf_insert(branch$child[[pos + 1L]], hashPos + 1L, leaf$k[i], leaf$v[[i]])$node)
      }
      pos <- .hmip_hht_nsu(occ, hashChunk) - 1L
      r <- .hmip_hht_insert_recurse(branch$child[[pos + 1L]], hashPos + 1L, key, val)
      branch$child[pos + 1L] <- list(r$node)
      return(list(node = branch, inserted = r$inserted))
    }
    leaf <- .hmip_hht_leaf_rehash(leaf, hashPos + 1L)
    r <- .hmip_hht_insert_recurse(leaf, hashPos + 1L, key, val)
    branch$child[1L] <- list(r$node)
    return(list(node = branch, inserted = r$inserted))
  }
  # branch node
  chunk <- hh$c6[hashPos + 1L]
  location <- .hmip_hht_nsu(node$occ, chunk)
  if (node$occ[chunk + 1L]) {
    location <- location - 1L
  } else {
    node$child <- append(node$child, list(NULL), after = location)
    node$occ[chunk + 1L] <- TRUE
  }
  r <- .hmip_hht_insert_recurse(node$child[[location + 1L]], hashPos + 1L, key, val)
  node$child[location + 1L] <- list(r$node)
  list(node = node, inserted = r$inserted)
}

.hmip_hht_merge_into_leaf <- function(leaf, hashPos, m) {
  t <- .hmip_hht_type(m)
  if (t >= 1L && t <= 5L) {
    for (i in seq_along(m$k)) leaf <- .hmip_hht_leaf_insert(leaf, hashPos, m$k[i], m$v[[i]])$node
  }
  leaf
}

# erase_recurse; returns the new node
.hmip_hht_erase_recurse <- function(node, hashPos, key) {
  t <- .hmip_hht_type(node)
  if (t == 0L) return(NULL)
  if (t == 1L) {
    w <- which(node$k == key)
    if (length(w)) { node$k <- node$k[-w[1L]]; node$v <- node$v[-w[1L]] }
    if (!length(node$k)) return(NULL)
    return(node)
  }
  if (t <= 5L) {
    leaf <- .hmip_hht_leaf_erase(node, hashPos, key)
    if (is.null(leaf)) return(node)
    n <- length(leaf$k)
    if (t == 2L && n == 0L) return(NULL)
    if (t > 2L && n == .hmip_hht_cap[t - 2L]) leaf$t <- t - 1L
    return(leaf)
  }
  chunk <- .hmip_hht_hash(key)$c6[hashPos + 1L]
  if (!node$occ[chunk + 1L]) return(node)
  location <- .hmip_hht_nsu(node$occ, chunk) - 1L
  ch <- .hmip_hht_erase_recurse(node$child[[location + 1L]], hashPos + 1L, key)
  node$child[location + 1L] <- list(ch)
  if (!is.null(ch)) return(node)
  node$occ[chunk + 1L] <- FALSE
  # removeChildFromBranchNode
  newNumChild <- sum(node$occ)
  if (newNumChild * 6L <= 54L) {
    childEntries <- 0L
    for (i in 0:newNumChild) {
      childEntries <- childEntries + .hmip_hht_num_estimate(node$child[[i + 1L]])
      if (childEntries > 54L) break
    }
    if (childEntries < 54L) {
      childEntries <- 0L
      for (i in 0:newNumChild) childEntries <- childEntries + .hmip_hht_num_entries(node$child[[i + 1L]])
      if (childEntries < 54L) {
        leaf <- .hmip_hht_newleaf(.hmip_hht_size_class(childEntries))
        for (i in 0:newNumChild) leaf <- .hmip_hht_merge_into_leaf(leaf, hashPos, node$child[[i + 1L]])
        return(leaf)
      }
    }
  }
  node$child <- node$child[-(location + 1L)]
  node
}

.hmip_hht_find_recurse <- function(node, hashPos, key) {
  repeat {
    t <- .hmip_hht_type(node)
    if (t == 0L) return(NULL)
    if (t == 1L) {
      w <- which(node$k == key)
      return(if (length(w)) list(key = key, value = node$v[[w[1L]]]) else NULL)
    }
    if (t <= 5L) {
      pos <- .hmip_hht_leaf_find(node, hashPos, key)
      return(if (is.null(pos)) NULL else list(key = key, value = node$v[[pos + 1L]]))
    }
    chunk <- .hmip_hht_hash(key)$c6[hashPos + 1L]
    if (!node$occ[chunk + 1L]) return(NULL)
    node <- node$child[[.hmip_hht_nsu(node$occ, chunk)]]
    hashPos <- hashPos + 1L
  }
}

.hmip_hht_log2i <- function(mask) max(which(mask)) - 1L

.hmip_hht_common_leaf_leaf <- function(l1, l2) {
  matchMask <- l1$occ & l2$occ
  offset1 <- -1L; offset2 <- -1L
  s1 <- length(l1$k); s2 <- length(l2$k)
  while (any(matchMask)) {
    pos <- .hmip_hht_log2i(matchMask)
    matchMask[pos + 1L] <- FALSE
    i <- .hmip_hht_nsu(l1$occ, pos) + offset1
    while (.hmip_hht_first(.hmip_hht_H(l1, i)) != pos) { i <- i + 1L; offset1 <- offset1 + 1L }
    j <- .hmip_hht_nsu(l2$occ, pos) + offset2
    while (.hmip_hht_first(.hmip_hht_H(l2, j)) != pos) { j <- j + 1L; offset2 <- offset2 + 1L }
    repeat {
      if (l1$h[i + 1L] > l2$h[j + 1L]) {
        i <- i + 1L
        if (i == s1 || .hmip_hht_first(l1$h[i + 1L]) != pos) break
      } else if (l2$h[j + 1L] > l1$h[i + 1L]) {
        j <- j + 1L
        if (j == s2 || .hmip_hht_first(l2$h[j + 1L]) != pos) break
      } else {
        if (l1$k[i + 1L] == l2$k[j + 1L]) return(list(key = l1$k[i + 1L], value = l1$v[[i + 1L]]))
        i <- i + 1L
        if (i == s1 || .hmip_hht_first(l1$h[i + 1L]) != pos) break
        j <- j + 1L
        if (j == s2 || .hmip_hht_first(l2$h[j + 1L]) != pos) break
      }
    }
  }
  NULL
}

.hmip_hht_common_leaf <- function(leaf, n2, hashPos) {
  t2 <- .hmip_hht_type(n2)
  if (t2 >= 2L && t2 <= 5L) return(.hmip_hht_common_leaf_leaf(leaf, n2))
  if (t2 == 6L) {
    matchMask <- n2$occ & leaf$occ
    offset <- -1L
    size <- length(leaf$k)
    while (any(matchMask)) {
      pos <- .hmip_hht_log2i(matchMask)
      matchMask[pos + 1L] <- FALSE
      i <- .hmip_hht_nsu(leaf$occ, pos) + offset
      while (.hmip_hht_first(.hmip_hht_H(leaf, i)) != pos) { i <- i + 1L; offset <- offset + 1L }
      j <- .hmip_hht_nsu(n2$occ, pos) - 1L
      repeat {
        if (!is.null(.hmip_hht_find_recurse(n2$child[[j + 1L]], hashPos + 1L, leaf$k[i + 1L])))
          return(list(key = leaf$k[i + 1L], value = leaf$v[[i + 1L]]))
        i <- i + 1L
        if (!(i < size && .hmip_hht_first(leaf$h[i + 1L]) == pos)) break
      }
    }
  }
  NULL
}

.hmip_hht_common_recurse <- function(n1, n2, hashPos) {
  if (.hmip_hht_type(n1) > .hmip_hht_type(n2)) { t <- n1; n1 <- n2; n2 <- t }
  t1 <- .hmip_hht_type(n1)
  if (t1 == 0L) return(NULL)
  if (t1 == 1L) {
    for (i in seq_along(n1$k)) {
      if (!is.null(.hmip_hht_find_recurse(n2, hashPos, n1$k[i]))) return(list(key = n1$k[i], value = n1$v[[i]]))
    }
    return(NULL)
  }
  if (t1 <= 5L) return(.hmip_hht_common_leaf(n1, n2, hashPos))
  matchMask <- n1$occ & n2$occ
  while (any(matchMask)) {
    pos <- .hmip_hht_log2i(matchMask)
    matchMask[pos + 1L] <- FALSE
    m <- .hmip_hht_common_recurse(n1$child[[.hmip_hht_nsu(n1$occ, pos)]],
                                  n2$child[[.hmip_hht_nsu(n2$occ, pos)]], hashPos + 1L)
    if (!is.null(m)) return(m)
  }
  NULL
}

# public interface (trees are values: insert and erase return the new tree)
.hmip_hht_insert <- function(tree, key, val = NULL) .hmip_hht_insert_recurse(tree, 0L, key, val)
.hmip_hht_erase <- function(tree, key) .hmip_hht_erase_recurse(tree, 0L, key)
.hmip_hht_contains <- function(tree, key) !is.null(.hmip_hht_find_recurse(tree, 0L, key))
.hmip_hht_find <- function(tree, key) .hmip_hht_find_recurse(tree, 0L, key)
.hmip_hht_find_common <- function(t1, t2) .hmip_hht_common_recurse(t1, t2, 0L)
.hmip_hht_empty <- function(tree) is.null(tree)

# replace the value of an existing key in place (insert_or_get followed by
# a write through the returned pointer)
.hmip_hht_set_value <- function(node, key, val, hashPos = 0L) {
  t <- .hmip_hht_type(node)
  if (t == 0L) return(node)
  if (t == 1L) {
    w <- which(node$k == key)
    if (length(w)) node$v[w[1L]] <- list(val)
    return(node)
  }
  if (t <= 5L) {
    pos <- .hmip_hht_leaf_find(node, hashPos, key)
    if (!is.null(pos)) node$v[pos + 1L] <- list(val)
    return(node)
  }
  chunk <- .hmip_hht_hash(key)$c6[hashPos + 1L]
  if (!node$occ[chunk + 1L]) return(node)
  loc <- .hmip_hht_nsu(node$occ, chunk)
  node$child[loc] <- list(.hmip_hht_set_value(node$child[[loc]], key, val, hashPos + 1L))
  node
}

# the for_each order: list(k = keys, v = values)
.hmip_hht_entries <- function(node) {
  t <- .hmip_hht_type(node)
  if (t == 0L) return(list(k = integer(), v = list()))
  if (t <= 5L) return(list(k = node$k, v = node$v))
  ks <- integer(); vs <- list()
  for (ch in node$child) {
    e <- .hmip_hht_entries(ch)
    ks <- c(ks, e$k); vs <- c(vs, e$v)
  }
  list(k = ks, v = vs)
}
.hmip_hht_keys <- function(node) .hmip_hht_entries(node)$k
# apply f to every value in for_each order (for_each with a mutating lambda)
.hmip_hht_map_values <- function(node, f) {
  t <- .hmip_hht_type(node)
  if (t == 0L) return(node)
  if (t <= 5L) {
    for (i in seq_along(node$k)) node$v[i] <- list(f(node$k[i], node$v[[i]]))
    return(node)
  }
  for (i in seq_along(node$child)) node$child[i] <- list(.hmip_hht_map_values(node$child[[i]], f))
  node
}

# ---------------------------------------------------------------------------
# HighsCliqueTable
#
# A CliqueVar is its index 2 * col + val (an integer). The cliques are the
# parallel vectors c_start, c_end, c_origin, c_numZeroFixed, c_equality; the
# substitutions are subst_substcol and subst_replace; the clique extensions
# are ext_row and ext_var. invertedHashList[[v + 1L]] is the
# HighsHashTree<HighsInt, HighsInt> (clique id to entry position) of the
# CliqueVar v, invertedHashListSizeTwo[[v + 1L]] the key only tree of its
# size two cliques. sizeTwoCliques maps the sorted edge to the clique id.
# ---------------------------------------------------------------------------

.hmip_cq_IInf <- .Machine$integer.max
.hmip_cv <- function(col, val) as.integer(2L * col + val)
.hmip_cv_col <- function(v) v %/% 2L
.hmip_cv_val <- function(v) v %% 2L
.hmip_cv_complement <- function(v) bitwXor(v, 1L)
.hmip_cv_weight <- function(v, sol) {
  s <- sol[v %/% 2L + 1L]
  if (v %% 2L == 1L) s else 1.0 - s
}

.hmip_cq_numCol <- function(mipsolver) mipsolver$model_$num_col_
.hmip_cq_numNonzero <- function(mipsolver) {
  m <- mipsolver$model_
  m$a_matrix_$start_[m$num_col_ + 1L]
}
.hmip_HighsDomain_numModelNonzeros <- function(dom) .hmip_cq_numNonzero(dom$mipsolver)

.hmip_new_HighsCliqueTable <- function(ncols) {
  ct <- new.env(parent = emptyenv())
  ct$cliqueentries <- integer()
  ct$invertedHashList <- vector("list", 2L * ncols)
  ct$invertedHashListSizeTwo <- vector("list", 2L * ncols)
  ct$sizeTwoCliques <- new.env(parent = emptyenv())
  ct$sizeTwoCount <- 0L
  ct$free_first <- integer(); ct$free_second <- integer()
  ct$freeslots <- integer()
  ct$c_start <- integer(); ct$c_end <- integer(); ct$c_origin <- integer()
  ct$c_numZeroFixed <- integer(); ct$c_equality <- logical()
  ct$numcliquesvar <- integer(2L * ncols)
  ct$infeasvertexstack <- integer()
  ct$colsubstituted <- integer(ncols)
  ct$subst_substcol <- integer(); ct$subst_replace <- integer()
  ct$deletedrows <- integer()
  ct$ext_row <- integer(); ct$ext_var <- integer()
  ct$iscandidate <- logical()
  ct$colDeleted <- logical(ncols)
  ct$cliquehits <- integer()
  ct$cliquehitinds <- integer()
  ct$randgen <- .hs_rng(0)
  ct$nfixings <- 0L
  ct$numEntries <- 0L
  ct$maxEntries <- .hmip_cq_IInf
  ct$minEntriesForParallelism <- .hmip_cq_IInf
  ct$inPresolve <- FALSE
  ct$numNeighbourhoodQueries <- 0
  ct
}

# HighsHashTable<std::pair<CliqueVar, CliqueVar>, HighsInt> sizeTwoCliques
.hmip_cq_edgekey <- function(v1, v2) {
  if (v1 %/% 2L > v2 %/% 2L) paste0(v2, "_", v1) else paste0(v1, "_", v2)
}
.hmip_cq_s2_find <- function(ct, v1, v2) {
  k <- .hmip_cq_edgekey(v1, v2)
  if (exists(k, envir = ct$sizeTwoCliques, inherits = FALSE)) get(k, envir = ct$sizeTwoCliques) else NULL
}
.hmip_cq_s2_insert <- function(ct, v1, v2, id) {
  k <- .hmip_cq_edgekey(v1, v2)
  if (exists(k, envir = ct$sizeTwoCliques, inherits = FALSE)) return(FALSE)
  assign(k, id, envir = ct$sizeTwoCliques)
  ct$sizeTwoCount <- ct$sizeTwoCount + 1L
  TRUE
}
.hmip_cq_s2_erase <- function(ct, v1, v2) {
  k <- .hmip_cq_edgekey(v1, v2)
  if (exists(k, envir = ct$sizeTwoCliques, inherits = FALSE)) {
    rm(list = k, envir = ct$sizeTwoCliques)
    ct$sizeTwoCount <- ct$sizeTwoCount - 1L
  }
}

.hmip_cq_ihl <- function(ct, v) ct$invertedHashList[[v + 1L]]
.hmip_cq_ihl2 <- function(ct, v) ct$invertedHashListSizeTwo[[v + 1L]]
.hmip_cq_set_ihl <- function(ct, v, tree) ct$invertedHashList[v + 1L] <- list(tree)
.hmip_cq_set_ihl2 <- function(ct, v, tree) ct$invertedHashListSizeTwo[v + 1L] <- list(tree)

.hmip_HighsCliqueTable_unlink <- function(ct, pos, cliqueid) {
  v <- ct$cliqueentries[pos + 1L]
  ct$numcliquesvar[v + 1L] <- ct$numcliquesvar[v + 1L] - 1L
  if (ct$c_end[cliqueid + 1L] - ct$c_start[cliqueid + 1L] == 2L)
    .hmip_cq_set_ihl2(ct, v, .hmip_hht_erase(.hmip_cq_ihl2(ct, v), cliqueid))
  else
    .hmip_cq_set_ihl(ct, v, .hmip_hht_erase(.hmip_cq_ihl(ct, v), cliqueid))
}

.hmip_HighsCliqueTable_link <- function(ct, pos, cliqueid) {
  v <- ct$cliqueentries[pos + 1L]
  ct$numcliquesvar[v + 1L] <- ct$numcliquesvar[v + 1L] + 1L
  if (ct$c_end[cliqueid + 1L] - ct$c_start[cliqueid + 1L] == 2L)
    .hmip_cq_set_ihl2(ct, v, .hmip_hht_insert(.hmip_cq_ihl2(ct, v), cliqueid)$node)
  else
    .hmip_cq_set_ihl(ct, v, .hmip_hht_insert(.hmip_cq_ihl(ct, v), cliqueid, pos)$node)
}

# findCommonCliqueId(numQueries, v1, v2): qenv holds the query counter
# numNeighbourhoodQueries (the table itself or a BronKerboschData)
.hmip_HighsCliqueTable_findCommonCliqueId <- function(ct, v1, v2, qenv = ct) {
  qenv$numNeighbourhoodQueries <- qenv$numNeighbourhoodQueries + 1
  if (!is.null(.hmip_cq_ihl2(ct, v1)) && !is.null(.hmip_cq_ihl2(ct, v2))) {
    id <- .hmip_cq_s2_find(ct, v1, v2)
    if (!is.null(id)) return(id)
  }
  r <- .hmip_hht_find_common(.hmip_cq_ihl(ct, v1), .hmip_cq_ihl(ct, v2))
  if (!is.null(r)) r$key else -1L
}

.hmip_HighsCliqueTable_haveCommonClique <- function(ct, v1, v2, qenv = ct) {
  if (v1 %/% 2L == v2 %/% 2L) return(FALSE)
  .hmip_HighsCliqueTable_findCommonCliqueId(ct, v1, v2, qenv) != -1L
}

# findCommonClique: the CliqueVar indices of the clique entries in stored
# order, NULL when v1 == v2 or there is no common clique
.hmip_HighsCliqueTable_findCommonClique <- function(ct, v1, v2) {
  if (v1 == v2) return(NULL)
  clq <- .hmip_HighsCliqueTable_findCommonCliqueId(ct, v1, v2)
  if (clq == -1L) return(NULL)
  s <- ct$c_start[clq + 1L]; e <- ct$c_end[clq + 1L]
  ct$cliqueentries[(s + 1L):e]
}

.hmip_HighsCliqueTable_resolveSubstitution <- function(ct, v) {
  while (ct$colsubstituted[v %/% 2L + 1L] != 0L) {
    rep <- ct$subst_replace[ct$colsubstituted[v %/% 2L + 1L]]
    v <- if (v %% 2L == 1L) rep else .hmip_cv_complement(rep)
  }
  v
}

# resolveSubstitution(col, val, offset) -> list(col, val, offset)
.hmip_HighsCliqueTable_resolveSubstitution3 <- function(ct, col, val, offset) {
  while (ct$colsubstituted[col + 1L] != 0L) {
    rep <- ct$subst_replace[ct$colsubstituted[col + 1L]]
    if (rep %% 2L == 0L) {
      offset <- offset + val
      val <- -val
    }
    col <- rep %/% 2L
  }
  list(col = col, val = val, offset = offset)
}

.hmip_HighsCliqueTable_getSubstitution <- function(ct, col) {
  k <- ct$colsubstituted[col + 1L]
  if (k == 0L) NULL else list(substcol = ct$subst_substcol[k], replace = ct$subst_replace[k])
}
.hmip_HighsCliqueTable_getSubstitutions <- function(ct) list(substcol = ct$subst_substcol, replace = ct$subst_replace)
.hmip_HighsCliqueTable_getDeletedRows <- function(ct) ct$deletedrows
.hmip_HighsCliqueTable_getCliqueExtensions <- function(ct) list(first = ct$ext_row, second = ct$ext_var)
.hmip_HighsCliqueTable_setPresolveFlag <- function(ct, inPresolve) ct$inPresolve <- as.logical(inPresolve)
.hmip_HighsCliqueTable_getPresolveFlag <- function(ct) ct$inPresolve
.hmip_HighsCliqueTable_getNumEntries <- function(ct) ct$numEntries
.hmip_HighsCliqueTable_setMaxEntries <- function(ct, numNz) ct$maxEntries <- 2000000L + 10L * numNz
.hmip_HighsCliqueTable_setMinEntriesForParallelism <- function(ct, m) ct$minEntriesForParallelism <- m
.hmip_HighsCliqueTable_isFull <- function(ct) ct$numEntries >= ct$maxEntries
.hmip_HighsCliqueTable_getNumFixings <- function(ct) ct$nfixings
.hmip_HighsCliqueTable_numCliques <- function(ct) length(ct$c_start) - length(ct$freeslots)
# cliques.size(), the slot count (separateCliques sizes cliquehits by it)
.hmip_HighsCliqueTable_cliquesSize <- function(ct) length(ct$c_start)
.hmip_HighsCliqueTable_numCliquesVar <- function(ct, v) ct$numcliquesvar[v + 1L]
.hmip_HighsCliqueTable_numCliques2 <- function(ct, col, val) ct$numcliquesvar[2L * col + val + 1L]

# std::remove_if on an integer vector: list(v, n) with the kept prefix
.hmip_cq_remove_if <- function(v, pred) {
  n <- length(v)
  first <- 1L
  while (first <= n && !pred(v[first])) first <- first + 1L
  if (first > n) return(list(v = v, n = n))
  res <- first
  i <- first + 1L
  while (i <= n) {
    if (!pred(v[i])) { v[res] <- v[i]; res <- res + 1L }
    i <- i + 1L
  }
  list(v = v, n = res - 1L)
}

.hmip_cq_free_insert <- function(ct, a, b) {
  f <- ct$free_first; s <- ct$free_second
  if (any(f == a & s == b)) return()
  k <- sum(f < a | (f == a & s < b))
  ct$free_first <- append(f, a, after = k)
  ct$free_second <- append(s, b, after = k)
}

.hmip_HighsCliqueTable_doAddClique <- function(ct, cliquevars, numcliquevars = length(cliquevars),
                                               equality = FALSE, origin = .hmip_cq_IInf) {
  if (!length(ct$freeslots)) {
    cliqueid <- length(ct$c_start)
    ct$c_start <- c(ct$c_start, 0L); ct$c_end <- c(ct$c_end, 0L); ct$c_origin <- c(ct$c_origin, 0L)
    ct$c_numZeroFixed <- c(ct$c_numZeroFixed, 0L); ct$c_equality <- c(ct$c_equality, FALSE)
  } else {
    nf <- length(ct$freeslots)
    cliqueid <- ct$freeslots[nf]
    ct$freeslots <- ct$freeslots[-nf]
  }
  id1 <- cliqueid + 1L
  ct$c_equality[id1] <- equality
  ct$c_origin[id1] <- as.integer(origin)
  f <- ct$free_first; s <- ct$free_second
  it <- if (length(f)) which(f >= numcliquevars)[1L] else NA
  if (is.na(it)) {
    ct$c_start[id1] <- length(ct$cliqueentries)
    ct$c_end[id1] <- ct$c_start[id1] + numcliquevars
    maxEnd <- ct$c_end[id1]
    ct$cliqueentries <- c(ct$cliqueentries, integer(maxEnd - length(ct$cliqueentries)))
  } else {
    fs1 <- f[it]; fs2 <- s[it]
    ct$free_first <- f[-it]; ct$free_second <- s[-it]
    ct$c_start[id1] <- fs2
    ct$c_end[id1] <- fs2 + numcliquevars
    maxEnd <- fs2 + fs1
  }
  ct$c_numZeroFixed[id1] <- 0L
  fixtozero <- FALSE
  start <- ct$c_start[id1]
  k <- start
  for (i in seq_len(numcliquevars)) {
    v <- .hmip_HighsCliqueTable_resolveSubstitution(ct, cliquevars[i])
    if (fixtozero) {
      ct$infeasvertexstack <- c(ct$infeasvertexstack, v)
      next
    }
    vc <- .hmip_cv_complement(v)
    clqHasVCompl <- if (numcliquevars == 2L) .hmip_hht_contains(.hmip_cq_ihl2(ct, vc), cliqueid)
      else .hmip_hht_contains(.hmip_cq_ihl(ct, vc), cliqueid)
    if (clqHasVCompl) {
      fixtozero <- TRUE
      j <- start
      while (j != k) {
        if (ct$cliqueentries[j + 1L] %/% 2L != v %/% 2L)
          ct$infeasvertexstack <- c(ct$infeasvertexstack, ct$cliqueentries[j + 1L])
        .hmip_HighsCliqueTable_unlink(ct, j, cliqueid)
        j <- j + 1L
      }
      k <- start
      next
    }
    if (numcliquevars == 2L) {
      r <- .hmip_hht_insert(.hmip_cq_ihl2(ct, v), cliqueid)
      .hmip_cq_set_ihl2(ct, v, r$node)
    } else {
      r <- .hmip_hht_insert(.hmip_cq_ihl(ct, v), cliqueid, k)
      .hmip_cq_set_ihl(ct, v, r$node)
    }
    if (!r$inserted) {
      ct$infeasvertexstack <- c(ct$infeasvertexstack, v)
      next
    }
    ct$cliqueentries[k + 1L] <- v
    ct$numcliquesvar[v + 1L] <- ct$numcliquesvar[v + 1L] + 1L
    k <- k + 1L
  }
  if (maxEnd > k) {
    if (length(ct$cliqueentries) == maxEnd) ct$cliqueentries <- ct$cliqueentries[seq_len(k)]
    else .hmip_cq_free_insert(ct, maxEnd - k, k)
    if (ct$c_end[id1] > k) {
      sz <- k - start
      if (sz == 0L) {
        ct$c_start[id1] <- -1L; ct$c_end[id1] <- -1L
        ct$freeslots <- c(ct$freeslots, cliqueid)
        return(invisible(NULL))
      } else if (sz == 1L) {
        .hmip_HighsCliqueTable_unlink(ct, start, cliqueid)
        ct$c_start[id1] <- -1L; ct$c_end[id1] <- -1L
        ct$freeslots <- c(ct$freeslots, cliqueid)
        return(invisible(NULL))
      } else if (sz == 2L) {
        .hmip_HighsCliqueTable_unlink(ct, start, cliqueid)
        .hmip_HighsCliqueTable_unlink(ct, start + 1L, cliqueid)
        ct$c_end[id1] <- k
        .hmip_HighsCliqueTable_link(ct, start, cliqueid)
        .hmip_HighsCliqueTable_link(ct, start + 1L, cliqueid)
      } else ct$c_end[id1] <- k
    }
  }
  cliqueLen <- ct$c_end[id1] - ct$c_start[id1]
  ct$numEntries <- ct$numEntries + cliqueLen
  if (cliqueLen == 2L)
    .hmip_cq_s2_insert(ct, ct$cliqueentries[ct$c_start[id1] + 1L], ct$cliqueentries[ct$c_start[id1] + 2L], cliqueid)
  invisible(NULL)
}

.hmip_HighsCliqueTable_removeClique <- function(ct, cliqueid) {
  id1 <- cliqueid + 1L
  o <- ct$c_origin[id1]
  if (o != .hmip_cq_IInf && o != -1L) ct$deletedrows <- c(ct$deletedrows, o)
  start <- ct$c_start[id1]; end <- ct$c_end[id1]
  len <- end - start
  if (len == 2L) .hmip_cq_s2_erase(ct, ct$cliqueentries[start + 1L], ct$cliqueentries[start + 2L])
  i <- start
  while (i != end) { .hmip_HighsCliqueTable_unlink(ct, i, cliqueid); i <- i + 1L }
  ct$freeslots <- c(ct$freeslots, cliqueid)
  .hmip_cq_free_insert(ct, len, start)
  ct$c_start[id1] <- -1L; ct$c_end[id1] <- -1L
  ct$numEntries <- ct$numEntries - len
  invisible(NULL)
}

# queryNeighbourhood: the 0-based positions i of q with a common clique with
# v. The parallel branch of HiGHS gives the same sorted indices and the same
# query count, so it is ported serially.
.hmip_HighsCliqueTable_queryNeighbourhood <- function(ct, v, q, qenv = ct) {
  if (ct$numcliquesvar[v + 1L] == 0L) return(integer())
  out <- integer()
  for (i in seq_along(q)) {
    if (.hmip_HighsCliqueTable_haveCommonClique(ct, v, q[i], qenv)) out <- c(out, i - 1L)
  }
  out
}
# partitionNeighbourhood: list(q, n, neighbourhoodInds)
.hmip_HighsCliqueTable_partitionNeighbourhood <- function(ct, v, q, qenv = ct) {
  inds <- .hmip_HighsCliqueTable_queryNeighbourhood(ct, v, q, qenv)
  for (i in seq_along(inds)) {
    t <- q[i]; q[i] <- q[inds[i] + 1L]; q[inds[i] + 1L] <- t
  }
  list(q = q, n = length(inds), neighbourhoodInds = inds)
}
# shrinkToNeighbourhood: list(q, n, neighbourhoodInds); q[1..n] is the result
.hmip_HighsCliqueTable_shrinkToNeighbourhood <- function(ct, v, q, qenv = ct) {
  inds <- .hmip_HighsCliqueTable_queryNeighbourhood(ct, v, q, qenv)
  for (i in seq_along(inds)) q[i] <- q[inds[i] + 1L]
  list(q = q, n = length(inds), neighbourhoodInds = inds)
}

# BronKerboschData (an environment, also the query counter)
.hmip_new_BronKerboschData <- function(sol) {
  d <- new.env(parent = emptyenv())
  d$sol <- sol
  d$P <- integer(); d$R <- integer(); d$Z <- integer()
  d$cliques <- list()
  d$neighbourhoodInds <- integer()
  d$wR <- 0.0; d$minW <- 1.05; d$feastol <- 1e-6
  d$ncalls <- 0L; d$maxcalls <- 10000L; d$maxcliques <- 100L
  d$maxNeighbourhoodQueries <- Inf
  d$numNeighbourhoodQueries <- 0
  d
}
.hmip_BronKerboschData_stop <- function(d) {
  d$maxcalls == d$ncalls || length(d$cliques) == d$maxcliques ||
    d$numNeighbourhoodQueries > d$maxNeighbourhoodQueries
}

.hmip_HighsCliqueTable_bronKerboschRecurse <- function(ct, data, Plen, X, Xlen = length(X)) {
  X <- X[seq_len(Xlen)]
  sol <- data$sol
  w <- data$wR
  for (i in seq_len(Plen)) w <- w + .hmip_cv_weight(data$P[i], sol)
  if (w < data$minW - data$feastol) return(invisible(NULL))
  if (Plen == 0L && Xlen == 0L) {
    clique <- data$R
    if (data$minW < w - data$feastol) {
      data$maxcliques <- data$maxcliques - length(data$cliques)
      data$cliques <- list()
      data$minW <- w
    }
    data$cliques[[length(data$cliques) + 1L]] <- clique
    return(invisible(NULL))
  }
  data$ncalls <- data$ncalls + 1L
  if (.hmip_BronKerboschData_stop(data)) return(invisible(NULL))
  pivweight <- -1.0
  pivot <- 0L
  for (i in seq_len(Xlen)) {
    wx <- .hmip_cv_weight(X[i], sol)
    if (wx > pivweight) {
      pivweight <- wx; pivot <- X[i]
      if (pivweight >= 1.0 - data$feastol) break
    }
  }
  if (pivweight < 1.0 - data$feastol) {
    for (i in seq_len(Plen)) {
      wp <- .hmip_cv_weight(data$P[i], sol)
      if (wp > pivweight) {
        pivweight <- wp; pivot <- data$P[i]
        if (pivweight >= 1.0 - data$feastol) break
      }
    }
  }
  nb <- .hmip_HighsCliqueTable_queryNeighbourhood(ct, pivot, data$P[seq_len(Plen)], data)
  nb <- c(nb, Plen)
  PminusNu <- integer()
  k <- 0L
  for (i in nb) {
    while (k < i) { PminusNu <- c(PminusNu, data$P[k + 1L]); k <- k + 1L }
    k <- k + 1L
  }
  if (length(PminusNu)) {
    wts <- vapply(PminusNu, .hmip_cv_weight, 0, sol = sol)
    PminusNu <- PminusNu[order(-wts, -PminusNu)]
  }
  localX <- X
  for (v in PminusNu) {
    r <- .hmip_HighsCliqueTable_partitionNeighbourhood(ct, v, data$P[seq_len(Plen)], data)
    data$P[seq_len(Plen)] <- r$q
    newPlen <- r$n
    r <- .hmip_HighsCliqueTable_partitionNeighbourhood(ct, v, localX, data)
    localX <- r$q
    newXlen <- r$n
    data$R <- c(data$R, v)
    wv <- .hmip_cv_weight(v, sol)
    data$wR <- data$wR + wv
    .hmip_HighsCliqueTable_bronKerboschRecurse(ct, data, newPlen, localX[seq_len(newXlen)])
    if (.hmip_BronKerboschData_stop(data)) return(invisible(NULL))
    data$R <- data$R[-length(data$R)]
    data$wR <- data$wR - wv
    w <- w - wv
    if (w < data$minW) return(invisible(NULL))
    vpos <- -1L
    i <- newPlen
    while (i != Plen) {
      if (data$P[i + 1L] == v) { vpos <- i; break }
      i <- i + 1L
    }
    Plen <- Plen - 1L
    t <- data$P[vpos + 1L]; data$P[vpos + 1L] <- data$P[Plen + 1L]; data$P[Plen + 1L] <- t
    localX <- c(localX, v)
  }
  invisible(NULL)
}

.hmip_HighsCliqueTable_computeMaximalCliques <- function(ct, vars, feastol) {
  if (!length(vars)) return(list())
  maxcolindex <- max(vars %/% 2L)
  sol <- numeric(maxcolindex + 1L)
  for (v in vars) sol[v %/% 2L + 1L] <- v %% 2L
  data <- .hmip_new_BronKerboschData(sol)
  data$feastol <- feastol
  for (v in vars) {
    col <- v %/% 2L
    if (ct$colsubstituted[col + 1L] != 0L || ct$colDeleted[col + 1L]) next
    if (ct$numcliquesvar[v + 1L] != 0L) data$P <- c(data$P, v)
  }
  .hmip_HighsCliqueTable_bronKerboschRecurse(ct, data, length(data$P), integer())
  data$cliques
}

# runCliqueSubsumption -> list(nremoved, clique)
.hmip_HighsCliqueTable_runCliqueSubsumption <- function(ct, globaldom, clique) {
  if (length(clique) == 2L) return(list(nremoved = 0L, clique = clique))
  nremoved <- 0L
  redundant <- FALSE
  nc <- length(ct$c_start)
  if (length(ct$cliquehits) < nc) ct$cliquehits <- c(ct$cliquehits, integer(nc - length(ct$cliquehits)))
  clique <- clique[!ct$colDeleted[clique %/% 2L + 1L]]
  hit <- function(cliqueid) {
    if (ct$cliquehits[cliqueid + 1L] == 0L) ct$cliquehitinds <- c(ct$cliquehitinds, cliqueid)
    ct$cliquehits[cliqueid + 1L] <- ct$cliquehits[cliqueid + 1L] + 1L
  }
  for (v in clique) {
    for (id in .hmip_hht_keys(.hmip_cq_ihl(ct, v))) hit(id)
    for (id in .hmip_hht_keys(.hmip_cq_ihl2(ct, v))) hit(id)
  }
  for (cliqueid in ct$cliquehitinds) {
    hits <- ct$cliquehits[cliqueid + 1L]
    ct$cliquehits[cliqueid + 1L] <- 0L
    len <- ct$c_end[cliqueid + 1L] - ct$c_start[cliqueid + 1L] - ct$c_numZeroFixed[cliqueid + 1L]
    if (hits == length(clique)) {
      redundant <- TRUE
    } else if (len == hits) {
      if (ct$c_equality[cliqueid + 1L]) {
        for (v in clique) {
          sizeTwo <- ct$c_end[cliqueid + 1L] - ct$c_start[cliqueid + 1L] == 2L
          vHasClq <- if (sizeTwo) .hmip_hht_contains(.hmip_cq_ihl2(ct, v), cliqueid)
            else .hmip_hht_contains(.hmip_cq_ihl(ct, v), cliqueid)
          if (!vHasClq) ct$infeasvertexstack <- c(ct$infeasvertexstack, v)
        }
      } else {
        nremoved <- nremoved + 1L
        ct$c_origin[cliqueid + 1L] <- .hmip_cq_IInf
        .hmip_HighsCliqueTable_removeClique(ct, cliqueid)
      }
    }
  }
  ct$cliquehitinds <- integer()
  if (redundant) clique <- integer()
  if (length(ct$infeasvertexstack))
    clique <- clique[!vapply(clique, function(v) .hmip_HighsDomain_isFixed(globaldom, v %/% 2L), TRUE)]
  list(nremoved = nremoved, clique = clique)
}

# fix the CliqueVar v to zero as the C++ loops do; TRUE when infeasible
.hmip_cq_fixzero <- function(ct, globaldom, v) {
  col <- v %/% 2L
  wasfixed <- .hmip_HighsDomain_isFixed(globaldom, col)
  .hmip_HighsDomain_fixCol(globaldom, col, as.double(1L - v %% 2L))
  if (globaldom$infeasible_) return(TRUE)
  if (!wasfixed) {
    ct$nfixings <- ct$nfixings + 1L
    ct$infeasvertexstack <- c(ct$infeasvertexstack, v)
  }
  FALSE
}

.hmip_cq_fixzero_process <- function(ct, globaldom, v) {
  wasfixed <- .hmip_HighsDomain_isFixed(globaldom, v %/% 2L)
  .hmip_HighsDomain_fixCol(globaldom, v %/% 2L, as.double(1L - v %% 2L))
  if (!wasfixed) {
    ct$nfixings <- ct$nfixings + 1L
    ct$infeasvertexstack <- c(ct$infeasvertexstack, v)
    .hmip_HighsCliqueTable_processInfeasibleVertices(ct, globaldom)
  }
}

.hmip_HighsCliqueTable_processNewEdge <- function(ct, globaldom, v1, v2) {
  if (v1 %/% 2L == v2 %/% 2L) {
    if (v1 %% 2L == v2 %% 2L) {
      .hmip_cq_fixzero_process(ct, globaldom, v1)
      return(FALSE)
    }
    return(TRUE)
  }
  c1 <- .hmip_cv_complement(v1); c2 <- .hmip_cv_complement(v2)
  if (.hmip_HighsCliqueTable_haveCommonClique(ct, c1, v2)) {
    .hmip_cq_fixzero_process(ct, globaldom, v2)
    return(FALSE)
  } else if (.hmip_HighsCliqueTable_haveCommonClique(ct, c2, v1)) {
    .hmip_cq_fixzero_process(ct, globaldom, v1)
    return(FALSE)
  }
  commonclique <- .hmip_HighsCliqueTable_findCommonCliqueId(ct, c1, c2)
  if (commonclique == -1L) return(FALSE)
  fixothers <- function(a, b) {
    while (commonclique != -1L) {
      start <- ct$c_start[commonclique + 1L]; end <- ct$c_end[commonclique + 1L]
      i <- start
      while (i != end) {
        e <- ct$cliqueentries[i + 1L]
        if (!(e == a || e == b)) {
          if (.hmip_cq_fixzero(ct, globaldom, e)) return(TRUE)
        }
        i <- i + 1L
      }
      .hmip_HighsCliqueTable_removeClique(ct, commonclique)
      commonclique <<- .hmip_HighsCliqueTable_findCommonCliqueId(ct, a, b)
    }
    FALSE
  }
  if (fixothers(c1, c2)) return(TRUE)
  .hmip_HighsCliqueTable_processInfeasibleVertices(ct, globaldom)
  if (globaldom$infeasible_) return(FALSE)
  commonclique <- .hmip_HighsCliqueTable_findCommonCliqueId(ct, v1, v2)
  if (fixothers(v1, v2)) return(TRUE)
  .hmip_HighsCliqueTable_processInfeasibleVertices(ct, globaldom)
  if (.hmip_HighsDomain_isFixed(globaldom, v1 %/% 2L) || .hmip_HighsDomain_isFixed(globaldom, v2 %/% 2L) ||
      globaldom$infeasible_) return(TRUE)
  if (v2 %/% 2L < v1 %/% 2L) {
    if (v1 %% 2L == 1L) v2 <- .hmip_cv_complement(v2)
    substcol <- v1 %/% 2L; replace <- v2
  } else {
    if (v2 %% 2L == 1L) v1 <- .hmip_cv_complement(v1)
    substcol <- v2 %/% 2L; replace <- v1
  }
  ct$subst_substcol <- c(ct$subst_substcol, substcol)
  ct$subst_replace <- c(ct$subst_replace, replace)
  ct$colsubstituted[substcol + 1L] <- length(ct$subst_substcol)
  replacefun <- function(sv, rv) {
    ct$numcliquesvar[rv + 1L] <- ct$numcliquesvar[rv + 1L] + ct$numcliquesvar[sv + 1L]
    ct$numcliquesvar[sv + 1L] <- 0L
    e <- .hmip_hht_entries(.hmip_cq_ihl(ct, sv))
    for (i in seq_along(e$k)) {
      .hmip_cq_set_ihl(ct, rv, .hmip_hht_insert(.hmip_cq_ihl(ct, rv), e$k[i], e$v[[i]])$node)
      ct$cliqueentries[e$v[[i]] + 1L] <- rv
    }
    .hmip_cq_set_ihl(ct, sv, NULL)
    for (cliqueid in .hmip_hht_keys(.hmip_cq_ihl2(ct, sv))) {
      pos <- ct$c_start[cliqueid + 1L]; otherPos <- pos + 1L
      if (ct$cliqueentries[otherPos + 1L] == sv) { t <- pos; pos <- otherPos; otherPos <- t }
      .hmip_cq_set_ihl2(ct, rv, .hmip_hht_insert(.hmip_cq_ihl2(ct, rv), cliqueid)$node)
      ct$cliqueentries[pos + 1L] <- rv
      .hmip_cq_s2_erase(ct, sv, ct$cliqueentries[otherPos + 1L])
      .hmip_cq_s2_insert(ct, rv, ct$cliqueentries[otherPos + 1L], cliqueid)
    }
    .hmip_cq_set_ihl2(ct, sv, NULL)
  }
  replacefun(.hmip_cv(substcol, 1L), replace)
  replacefun(.hmip_cv(substcol, 0L), .hmip_cv_complement(replace))
  TRUE
}

# addClique: returns the cliquevars array as the C++ leaves it (resolved
# substitutions and the remove_if of the fixed entries)
.hmip_HighsCliqueTable_addClique <- function(ct, mipsolver, cliquevars, numcliquevars = length(cliquevars),
                                             equality = FALSE, origin = .hmip_cq_IInf) {
  md <- mipsolver$mipdata_
  globaldom <- md$domain
  for (i in seq_len(numcliquevars)) {
    cliquevars[i] <- .hmip_HighsCliqueTable_resolveSubstitution(ct, cliquevars[i])
    col <- cliquevars[i] %/% 2L
    if (.hmip_HighsDomain_isFixed(globaldom, col)) {
      if (cliquevars[i] %% 2L == globaldom$col_lower_[col + 1L]) {
        for (k in seq_len(numcliquevars)) {
          if (k == i) next
          if (.hmip_cq_fixzero(ct, globaldom, cliquevars[k])) return(invisible(cliquevars))
        }
        .hmip_HighsCliqueTable_processInfeasibleVertices(ct, globaldom)
        return(invisible(cliquevars))
      }
    }
  }
  if (numcliquevars <= 100L) {
    hasNewEdge <- FALSE
    for (i in seq_len(numcliquevars - 1L)) {
      vi <- cliquevars[i]
      if (.hmip_HighsDomain_isFixed(globaldom, vi %/% 2L)) next
      if (ct$numcliquesvar[vi + 1L] == 0L && ct$numcliquesvar[.hmip_cv_complement(vi) + 1L] == 0L) {
        hasNewEdge <- TRUE
        next
      }
      j <- i + 1L
      while (j <= numcliquevars) {
        vj <- cliquevars[j]
        if (.hmip_HighsDomain_isFixed(globaldom, vj %/% 2L)) { j <- j + 1L; next }
        if (.hmip_HighsCliqueTable_haveCommonClique(ct, vi, vj)) { j <- j + 1L; next }
        hasNewEdge <- TRUE
        iscover <- .hmip_HighsCliqueTable_processNewEdge(ct, globaldom, vi, vj)
        if (globaldom$infeasible_) return(invisible(cliquevars))
        nq <- md$nodequeue
        if (!.hmip_HighsNodeQueue_empty(nq)) {
          n1 <- if (vi %% 2L == 1L) .hmip_HighsNodeQueue_getUpNodes(nq, vi %/% 2L) else .hmip_HighsNodeQueue_getDownNodes(nq, vi %/% 2L)
          n2 <- if (vj %% 2L == 1L) .hmip_HighsNodeQueue_getUpNodes(nq, vj %/% 2L) else .hmip_HighsNodeQueue_getDownNodes(nq, vj %/% 2L)
          if (length(n1$first) && length(n2$first)) {
            # lower_bound and upper_bound of (val, kHighsIInf) are both the
            # first node with a bound value above val, so the ranges are
            # always empty in HiGHS 1.12.0
            lb1 <- function(n, val) {
              w <- which(n$first > val | (n$first == val & n$second >= .hmip_cq_IInf))
              if (length(w)) w[1L] else length(n$first) + 1L
            }
            ub1 <- function(n, val) { w <- which(n$first > val); if (length(w)) w[1L] else length(n$first) + 1L }
            it1 <- lb1(n1, vi %% 2L); e1 <- ub1(n1, vi %% 2L)
            it2 <- lb1(n2, vj %% 2L); e2 <- ub1(n2, vj %% 2L)
            if (it1 != e1 && it2 != e2 &&
                (n1$second[it1] <= n2$second[e2 - 1L] || n2$second[it2] <= n1$second[e1 - 1L])) {
              while (it1 != e1 && it2 != e2) {
                if (n1$second[it1] < n2$second[it2]) it1 <- it1 + 1L
                else if (n2$second[it2] < n1$second[it1]) it2 <- it2 + 1L
                else {
                  prunedNode <- n2$second[it2]
                  it1 <- it1 + 1L; it2 <- it2 + 1L
                  md$pruned_treeweight <- .hmip_cd_iaddd(md$pruned_treeweight, .hmip_HighsNodeQueue_pruneNode(nq, prunedNode))
                }
              }
            }
          }
        }
        if (iscover) {
          for (k in seq_len(numcliquevars)) {
            if (k == i || k == j) next
            if (.hmip_cq_fixzero(ct, globaldom, cliquevars[k])) return(invisible(cliquevars))
          }
          .hmip_HighsCliqueTable_processInfeasibleVertices(ct, globaldom)
          return(invisible(cliquevars))
        }
        j <- j + 1L
      }
    }
    if (!hasNewEdge && origin == .hmip_cq_IInf) return(invisible(cliquevars))
  }
  r <- .hmip_cq_remove_if(cliquevars[seq_len(numcliquevars)], function(v) .hmip_HighsDomain_isFixed(globaldom, v %/% 2L))
  cliquevars[seq_len(numcliquevars)] <- r$v
  numcliquevars <- r$n
  if (numcliquevars < 2L) return(invisible(cliquevars))
  .hmip_HighsCliqueTable_doAddClique(ct, cliquevars[seq_len(numcliquevars)], numcliquevars, equality, origin)
  .hmip_HighsCliqueTable_processInfeasibleVertices(ct, globaldom)
  invisible(cliquevars)
}

# the clique loop shared by the extractCliques variants: perm is sorted,
# absval(p) the coefficient size, slack(k) the mincliqueval of position k,
# cv(p) the CliqueVar of position p; stopcheck() ends the loop early
.hmip_cq_cliqueloop <- function(ct, mipsolver, perm, nbin, absval, slack, cv, stopcheck = function() FALSE) {
  gd <- mipsolver$mipdata_$domain
  k <- nbin - 1L
  while (k != 0L && !stopcheck()) {
    mincliqueval <- slack(k)
    # partition_point over perm[0, k)
    cliqueend <- 0L
    while (cliqueend < k && absval(perm[cliqueend + 1L]) > mincliqueval) cliqueend <- cliqueend + 1L
    if (cliqueend == 0L) { k <- k - 1L; next }
    clique <- c(vapply(perm[seq_len(cliqueend)], cv, 0L), cv(perm[k + 1L]))
    if (length(clique) >= 2L) {
      .hmip_HighsCliqueTable_addClique(ct, mipsolver, clique)
      if (gd$infeasible_ || stopcheck()) return(invisible(NULL))
    }
    if (cliqueend == k) return(invisible(NULL))
    k <- k - 1L
  }
  invisible(NULL)
}

# extractCliques(mipsolver, inds, vals, complementation, rhs, nbin, perm,
# clique, feastol) for one transformed row
.hmip_HighsCliqueTable_extractCliquesRow <- function(ct, mipsolver, inds, vals, complementation, rhs, nbin, feastol) {
  md <- mipsolver$mipdata_
  implics <- md$implications
  gd <- md$domain
  n <- length(inds)
  pr <- .hmip_std_partition(seq_len(n) - 1L, function(pos) .hmip_HighsDomain_isBinary(gd, inds[pos + 1L]))
  perm <- pr$v
  nbin <- pr$k
  ntotal <- n
  if (nbin < ntotal) {
    for (i in seq_len(nbin)) {
      pi <- perm[i] + 1L
      bincol <- inds[pi]
      impliedub <- .cd_sub(c(rhs, 0), vals[pi])
      if (.hmip_HighsImplications_tooManyVarBounds(implics)) break
      for (j in (nbin + 1L):ntotal) {
        pj <- perm[j] + 1L
        col <- inds[pj]
        if (.hmip_HighsDomain_isFixed(gd, col)) next
        colub <- .cd_sub(c(gd$col_upper_[col + 1L], 0), gd$col_lower_[col + 1L])
        implcolub <- .cd_div(impliedub, vals[pj])
        if (.hmip_dom_vtype(gd, col) != 0L) implcolub <- c(floor(.cd_val(implcolub) + md$feastol), 0)
        if (.cd_val(implcolub) < .cd_val(.cd_sub(colub, feastol))) {
          if (complementation[pi] == -1L) {
            coef <- .cd_sub_cd(colub, implcolub); constant <- implcolub
          } else {
            coef <- .cd_sub_cd(implcolub, colub); constant <- colub
          }
          if (complementation[pj] == -1L) {
            constant <- .cd_sub(constant, gd$col_upper_[col + 1L])
            .hmip_HighsImplications_addVLB(implics, col, bincol, -.cd_val(coef), -.cd_val(constant))
          } else {
            constant <- .hmip_cd_iaddd(constant, gd$col_lower_[col + 1L])
            .hmip_HighsImplications_addVUB(implics, col, bincol, .cd_val(coef), .cd_val(constant))
          }
        }
      }
    }
  }
  if (nbin <= 1L) return(invisible(NULL))
  b <- perm[seq_len(nbin)]
  perm[seq_len(nbin)] <- b[order(-vals[b + 1L], -b)]
  if (vals[perm[1L] + 1L] + vals[perm[2L] + 1L] <= rhs + feastol) return(invisible(NULL))
  cv <- function(pos) if (complementation[pos + 1L] == -1L) .hmip_cv(inds[pos + 1L], 0L) else .hmip_cv(inds[pos + 1L], 1L)
  # upstream compares vals[0], not vals[perm[0]]
  if (abs(vals[1L] - vals[perm[nbin] + 1L]) <= feastol && rhs < 2 * vals[perm[nbin] + 1L] - feastol) {
    clique <- vapply(perm[seq_len(nbin)], cv, 0L)
    .hmip_HighsCliqueTable_addClique(ct, mipsolver, clique, nbin)
    return(invisible(NULL))
  }
  .hmip_cq_cliqueloop(ct, mipsolver, perm, nbin, function(p) vals[p + 1L],
                      function(k) rhs - vals[perm[k + 1L] + 1L] + feastol, cv)
}

.hmip_HighsCliqueTable_extractCliquesFromCut <- function(ct, mipsolver, inds, vals, len, rhs) {
  if (.hmip_HighsCliqueTable_isFull(ct)) return(invisible(NULL))
  md <- mipsolver$mipdata_
  implics <- md$implications
  gd <- md$domain
  feastol <- md$feastol
  minact <- c(0, 0)
  nbin <- 0L
  for (i in seq_len(len)) {
    col <- inds[i]
    if (.hmip_HighsDomain_isBinary(gd, col)) nbin <- nbin + 1L
    if (vals[i] > 0) {
      if (gd$col_lower_[col + 1L] == -Inf) return(invisible(NULL))
      minact <- .hmip_cd_iaddd(minact, vals[i] * gd$col_lower_[col + 1L])
    } else {
      if (gd$col_upper_[col + 1L] == Inf) return(invisible(NULL))
      minact <- .hmip_cd_iaddd(minact, vals[i] * gd$col_upper_[col + 1L])
    }
  }
  if (.cd_val(.hmip_cd_rsub(rhs, minact)) < 0.0) minact <- c(rhs, 0)
  for (i in seq_len(len)) {
    col <- inds[i]
    if (.hmip_dom_vtype(gd, col) == 0L) next
    boundVal <- .cd_val(.cd_div(.hmip_cd_rsub(rhs, minact), vals[i]))
    if (vals[i] > 0) {
      boundVal <- floor(boundVal + gd$col_lower_[col + 1L] + md$feastol)
      .hmip_HighsDomain_changeBound3(gd, 1L, col, boundVal, .hmip_reason_unspecified())
      if (gd$infeasible_) return(invisible(NULL))
    } else {
      boundVal <- ceiling(boundVal + gd$col_upper_[col + 1L] - md$feastol)
      .hmip_HighsDomain_changeBound3(gd, 0L, col, boundVal, .hmip_reason_unspecified())
      if (gd$infeasible_) return(invisible(NULL))
    }
  }
  if (nbin <= 1L) return(invisible(NULL))
  pr <- .hmip_std_partition(seq_len(len) - 1L, function(pos) .hmip_HighsDomain_isBinary(gd, inds[pos + 1L]))
  perm <- pr$v
  nbin <- pr$k
  if (nbin < len) {
    for (i in seq_len(nbin)) {
      pi <- perm[i] + 1L
      bincol <- inds[pi]
      impliedActivity <- .cd_sub(.hmip_cd_rsub(rhs, minact), abs(vals[pi]))
      for (j in (nbin + 1L):len) {
        pj <- perm[j] + 1L
        col <- inds[pj]
        if (.hmip_HighsDomain_isFixed(gd, col)) next
        if (vals[pj] > 0) {
          implcolub <- .cd_val(.hmip_cd_iaddd(impliedActivity, vals[pj] * gd$col_lower_[col + 1L])) / vals[pj]
          if (.hmip_dom_vtype(gd, col) != 0L) implcolub <- floor(implcolub + md$feastol)
          if (implcolub < gd$col_upper_[col + 1L] - feastol) {
            if (vals[pi] < 0) {
              coef <- gd$col_upper_[col + 1L] - implcolub; constant <- implcolub
            } else {
              if (gd$col_upper_[col + 1L] == Inf) next
              coef <- implcolub - gd$col_upper_[col + 1L]; constant <- gd$col_upper_[col + 1L]
            }
            .hmip_HighsImplications_addVUB(implics, col, bincol, coef, constant)
          }
        } else {
          implcollb <- .cd_val(.hmip_cd_iaddd(impliedActivity, vals[pj] * gd$col_upper_[col + 1L])) / vals[pj]
          if (.hmip_dom_vtype(gd, col) != 0L) implcollb <- ceiling(implcollb - md$feastol)
          if (implcollb > gd$col_lower_[col + 1L] + feastol) {
            if (vals[pi] < 0) {
              coef <- gd$col_lower_[col + 1L] - implcollb; constant <- implcollb
            } else {
              if (gd$col_lower_[col + 1L] == -Inf) next
              coef <- implcollb - gd$col_lower_[col + 1L]; constant <- gd$col_lower_[col + 1L]
            }
            .hmip_HighsImplications_addVLB(implics, col, bincol, coef, constant)
          }
        }
      }
    }
  }
  if (nbin <= 1L) return(invisible(NULL))
  b <- perm[seq_len(nbin)]
  perm[seq_len(nbin)] <- b[order(-abs(vals[b + 1L]), -b)]
  if (abs(vals[perm[1L] + 1L]) + abs(vals[perm[2L] + 1L]) <= .cd_val(.cd_add(.hmip_cd_rsub(rhs, minact), feastol)))
    return(invisible(NULL))
  maxNewEntries <- min(md$numCliqueEntriesAfterPresolve + 100000L + 4L * .hmip_HighsDomain_numModelNonzeros(gd),
                       ct$numEntries + 10L * nbin)
  cv <- function(pos) if (vals[pos + 1L] < 0) .hmip_cv(inds[pos + 1L], 0L) else .hmip_cv(inds[pos + 1L], 1L)
  .hmip_cq_cliqueloop(ct, mipsolver, perm, nbin, function(p) abs(vals[p + 1L]),
                      function(k) .cd_val(.cd_add(.cd_sub(.hmip_cd_rsub(rhs, minact), abs(vals[perm[k + 1L] + 1L])), feastol)),
                      cv, function() ct$numEntries >= maxNewEntries)
}

# extractCliques(mipsolver, transformRows)
.hmip_HighsCliqueTable_extractCliques <- function(ct, mipsolver, transformRows = TRUE) {
  md <- mipsolver$mipdata_; m <- mipsolver$model_
  gd <- md$domain
  for (i in seq_len(m$num_row_) - 1L) {
    start <- md$ARstart_[i + 1L]; end <- md$ARstart_[i + 2L]
    if (.hmip_HighsPostsolveStack_getOrigRowIndex(md$postSolveStack, i) >= mipsolver$orig_model_$num_row_) break
    if (m$row_upper_[i + 1L] == 1.0) {
      issetppc <- TRUE
      clique <- integer()
      j <- start
      while (j != end) {
        col <- md$ARindex_[j + 1L]
        if (gd$col_upper_[col + 1L] == 0.0 && gd$col_lower_[col + 1L] == 0.0) { j <- j + 1L; next }
        if (!.hmip_HighsDomain_isBinary(gd, col)) { issetppc <- FALSE; break }
        if (md$ARvalue_[j + 1L] != 1.0) { issetppc <- FALSE; break }
        clique <- c(clique, .hmip_cv(col, 1L))
        j <- j + 1L
      }
      if (issetppc) {
        equality <- m$row_lower_[i + 1L] == 1.0
        .hmip_HighsCliqueTable_addClique(ct, mipsolver, clique, length(clique), equality, i)
        if (gd$infeasible_) return(invisible(NULL))
        next
      }
    }
    if (!transformRows || .hmip_HighsCliqueTable_isFull(ct)) next
    offset <- 0
    entries <- .hmip_new_HighsHashTable()
    j <- start
    while (j != end) {
      r <- .hmip_HighsCliqueTable_resolveSubstitution3(ct, md$ARindex_[j + 1L], md$ARvalue_[j + 1L], offset)
      offset <- r$offset
      old <- .hmip_HighsHashTable_find(entries, r$col)
      slot <- .hmip_HighsHashTable_ref(entries, r$col, 0)
      .hmip_HighsHashTable_set(entries, slot, (if (is.null(old)) 0 else old) + r$val)
      j <- j + 1L
    }
    ekeys <- .hmip_HighsHashTable_keys(entries)
    evals <- unlist(.hmip_HighsHashTable_values(entries))  # B2 stores values in a list
    for (side in c(1, -1)) {
      if (side == 1 && m$row_upper_[i + 1L] == Inf) next
      if (side == -1 && m$row_lower_[i + 1L] == -Inf) next
      rhs <- if (side == 1) m$row_upper_[i + 1L] - offset else -m$row_lower_[i + 1L] + offset
      inds <- integer(); vals <- numeric(); complementation <- integer()
      freevar <- FALSE
      nbin <- 0L
      for (e in seq_along(ekeys)) {
        col <- ekeys[e]
        val <- if (side == 1) evals[e] else -evals[e]
        if (abs(val) < md$epsilon) next
        if (.hmip_HighsDomain_isBinary(gd, col)) nbin <- nbin + 1L
        if (val < 0) {
          if (gd$col_upper_[col + 1L] == Inf) { freevar <- TRUE; break }
          vals <- c(vals, -val); inds <- c(inds, col); complementation <- c(complementation, -1L)
          rhs <- rhs - val * gd$col_upper_[col + 1L]
        } else {
          if (gd$col_lower_[col + 1L] == -Inf) { freevar <- TRUE; break }
          vals <- c(vals, val); inds <- c(inds, col); complementation <- c(complementation, 1L)
          rhs <- rhs - val * gd$col_lower_[col + 1L]
        }
      }
      if (!freevar && nbin != 0L) {
        .hmip_HighsCliqueTable_extractCliquesRow(ct, mipsolver, inds, vals, complementation, rhs, nbin, md$feastol)
        if (gd$infeasible_) return(invisible(NULL))
      }
    }
  }
  invisible(NULL)
}

.hmip_HighsCliqueTable_extractObjCliques <- function(ct, mipsolver) {
  md <- mipsolver$mipdata_
  nbin <- .hmip_HighsObjectiveFunction_getNumBinariesInObjective(md$objectiveFunction)
  if (nbin <= 1L) return(invisible(NULL))
  gd <- md$domain
  if (.hmip_HighsDomain_getObjectiveLowerBound(gd) == -Inf) return(invisible(NULL))
  pc <- .hmip_HighsDomain_getCutoffConstraint(gd)
  vals <- pc$vals; inds <- pc$inds; len <- pc$len; rhs <- pc$rhs
  pr <- .hmip_std_partition(seq_len(nbin) - 1L, function(pos) vals[pos + 1L] != 0.0 && !.hmip_HighsDomain_isFixed(gd, inds[pos + 1L]))
  perm <- pr$v
  nbin <- pr$k
  if (nbin <= 1L) return(invisible(NULL))
  b <- perm[seq_len(nbin)]
  perm[seq_len(nbin)] <- b[order(-abs(vals[b + 1L]), -b)]
  a <- .hmip_HighsDomain_computeMinActivity(gd, 0L, len, inds, vals)
  minact <- a$activity
  feastol <- md$feastol
  if (abs(vals[perm[1L] + 1L]) + abs(vals[perm[2L] + 1L]) <= .cd_val(.cd_add(.hmip_cd_rsub(rhs, minact), feastol)))
    return(invisible(NULL))
  cv <- function(pos) if (vals[pos + 1L] < 0) .hmip_cv(inds[pos + 1L], 0L) else .hmip_cv(inds[pos + 1L], 1L)
  .hmip_cq_cliqueloop(ct, mipsolver, perm, nbin, function(p) abs(vals[p + 1L]),
                      function(k) .cd_val(.cd_add(.cd_sub(.hmip_cd_rsub(rhs, minact), abs(vals[perm[k + 1L] + 1L])), feastol)),
                      cv)
}

# cliquePartition(clqVars, partitionStart) -> list(clqVars, partitionStart)
.hmip_HighsCliqueTable_cliquePartition <- function(ct, clqVars) {
  clqVars <- .hs_rng_shuffle(ct$randgen, clqVars)
  numClqVars <- length(clqVars)
  extensionEnd <- numClqVars
  partitionStart <- 0L
  for (i in seq_len(numClqVars) - 1L) {
    if (i == extensionEnd) {
      partitionStart <- c(partitionStart, i)
      extensionEnd <- numClqVars
    }
    v <- clqVars[i + 1L]
    extensionStart <- i + 1L
    idx <- seq_len(extensionEnd - extensionStart) + extensionStart
    r <- .hmip_HighsCliqueTable_partitionNeighbourhood(ct, v, clqVars[idx])
    clqVars[idx] <- r$q
    extensionEnd <- r$n + extensionStart
  }
  list(clqVars = clqVars, partitionStart = c(partitionStart, numClqVars))
}

# cliquePartition(objective, clqVars, partitionStart) -> list(clqVars,
# partitionStart); the name B1 calls
.hmip_HighsCliqueTable_cliquePartitionObj <- function(ct, objective, clqVars) {
  clqVars <- .hs_rng_shuffle(ct$randgen, clqVars)
  less <- function(v1, v2) {
    (2 * (v1 %% 2L) - 1) * objective[v1 %/% 2L + 1L] > (2 * (v2 %% 2L) - 1) * objective[v2 %/% 2L + 1L]
  }
  clqVars <- .hmip_pdqsort(clqVars, less, branchless = TRUE)
  numClqVars <- length(clqVars)
  extensionEnd <- numClqVars
  partitionStart <- 0L
  lastSwappedIndex <- 0L
  for (i in seq_len(numClqVars) - 1L) {
    if (i == extensionEnd) {
      partitionStart <- c(partitionStart, i)
      extensionEnd <- numClqVars
      if (lastSwappedIndex >= i) {
        idx <- (i + 1L):(lastSwappedIndex + 1L)
        clqVars[idx] <- .hmip_pdqsort(clqVars[idx], less, branchless = TRUE)
      }
      lastSwappedIndex <- 0L
    }
    v <- clqVars[i + 1L]
    extensionStart <- i + 1L
    idx <- seq_len(extensionEnd - extensionStart) + extensionStart
    r <- .hmip_HighsCliqueTable_partitionNeighbourhood(ct, v, clqVars[idx])
    clqVars[idx] <- r$q
    extensionEnd <- r$n + extensionStart
    if (length(r$neighbourhoodInds))
      lastSwappedIndex <- max(r$neighbourhoodInds[length(r$neighbourhoodInds)] + extensionStart, lastSwappedIndex)
  }
  list(clqVars = clqVars, partitionStart = c(partitionStart, numClqVars))
}

.hmip_HighsCliqueTable_foundCover <- function(ct, globaldom, v1, v2) {
  equality <- FALSE
  commonclique <- .hmip_HighsCliqueTable_findCommonCliqueId(ct, v1, v2)
  if (commonclique != -1L) equality <- TRUE
  while (commonclique != -1L) {
    start <- ct$c_start[commonclique + 1L]; end <- ct$c_end[commonclique + 1L]
    i <- start
    while (i != end) {
      e <- ct$cliqueentries[i + 1L]
      if (!(e == v1 || e == v2)) {
        if (.hmip_cq_fixzero(ct, globaldom, e)) return(equality)
      }
      i <- i + 1L
    }
    .hmip_HighsCliqueTable_removeClique(ct, commonclique)
    commonclique <- .hmip_HighsCliqueTable_findCommonCliqueId(ct, v1, v2)
  }
  .hmip_HighsCliqueTable_processInfeasibleVertices(ct, globaldom)
  equality
}

# the fixings over the cliques of a moved out incidence list; TRUE when the
# domain became infeasible
.hmip_cq_fix_cliques <- function(ct, globaldom, ids, vcol) {
  for (cliqueid in ids) {
    start <- ct$c_start[cliqueid + 1L]; end <- ct$c_end[cliqueid + 1L]
    i <- start
    while (i != end) {
      e <- ct$cliqueentries[i + 1L]
      if (e %/% 2L != vcol) {
        if (.hmip_cq_fixzero(ct, globaldom, e)) return(TRUE)
      }
      i <- i + 1L
    }
    .hmip_HighsCliqueTable_removeClique(ct, cliqueid)
  }
  FALSE
}

.hmip_HighsCliqueTable_processInfeasibleVertices <- function(ct, globaldom) {
  while (length(ct$infeasvertexstack) && !globaldom$infeasible_) {
    n <- length(ct$infeasvertexstack)
    v <- .hmip_cv_complement(ct$infeasvertexstack[n])
    ct$infeasvertexstack <- ct$infeasvertexstack[-n]
    v <- .hmip_HighsCliqueTable_resolveSubstitution(ct, v)
    vcol <- v %/% 2L
    wasfixed <- .hmip_HighsDomain_isFixed(globaldom, vcol)
    .hmip_HighsDomain_fixCol(globaldom, vcol, as.double(v %% 2L))
    if (globaldom$infeasible_) return(invisible(NULL))
    if (!wasfixed) ct$nfixings <- ct$nfixings + 1L
    if (ct$colDeleted[vcol + 1L]) next
    ct$colDeleted[vcol + 1L] <- TRUE
    ids <- .hmip_hht_keys(.hmip_cq_ihl(ct, v))
    ids2 <- .hmip_hht_keys(.hmip_cq_ihl2(ct, v))
    .hmip_cq_set_ihl(ct, v, NULL); .hmip_cq_set_ihl2(ct, v, NULL)
    if (.hmip_cq_fix_cliques(ct, globaldom, ids, vcol)) return(invisible(NULL))
    if (.hmip_cq_fix_cliques(ct, globaldom, ids2, vcol)) return(invisible(NULL))
    vc <- .hmip_cv_complement(v)
    ids <- .hmip_hht_keys(.hmip_cq_ihl(ct, vc))
    ids2 <- .hmip_hht_keys(.hmip_cq_ihl2(ct, vc))
    .hmip_cq_set_ihl(ct, vc, NULL); .hmip_cq_set_ihl2(ct, vc, NULL)
    if (ct$inPresolve) {
      for (cliqueid in ids) {
        ct$c_numZeroFixed[cliqueid + 1L] <- ct$c_numZeroFixed[cliqueid + 1L] + 1L
        if (ct$c_end[cliqueid + 1L] - ct$c_start[cliqueid + 1L] - ct$c_numZeroFixed[cliqueid + 1L] <= 1L)
          .hmip_HighsCliqueTable_removeClique(ct, cliqueid)
      }
      next
    }
    for (cliqueid in ids2) .hmip_HighsCliqueTable_removeClique(ct, cliqueid)
    for (cliqueid in ids) {
      id1 <- cliqueid + 1L
      ct$c_numZeroFixed[id1] <- ct$c_numZeroFixed[id1] + 1L
      len <- ct$c_end[id1] - ct$c_start[id1]
      if (len - ct$c_numZeroFixed[id1] <= 1L) {
        .hmip_HighsCliqueTable_removeClique(ct, cliqueid)
      } else if (ct$c_numZeroFixed[id1] >= max(10L, bitwShiftR(len, 1L))) {
        clq <- ct$cliqueentries[(ct$c_start[id1] + 1L):ct$c_end[id1]]
        .hmip_HighsCliqueTable_removeClique(ct, cliqueid)
        clq <- clq[!vapply(clq, function(x) .hmip_HighsDomain_isFixed(globaldom, x %/% 2L) &&
                             globaldom$col_lower_[x %/% 2L + 1L] == 1L - x %% 2L, TRUE)]
        if (length(clq) > 1L) .hmip_HighsCliqueTable_doAddClique(ct, clq, length(clq))
      }
    }
  }
  .hmip_HighsCliqueTable_propagateAndCleanup(ct, globaldom)
}

.hmip_HighsCliqueTable_propagateAndCleanup <- function(ct, globaldom) {
  start <- length(globaldom$domchgstack_column)
  .hmip_HighsDomain_propagate(globaldom)
  end <- length(globaldom$domchgstack_column)
  while (!globaldom$infeasible_ && start != end) {
    k <- start
    while (k != end) {
      col <- globaldom$domchgstack_column[k + 1L]
      k <- k + 1L
      lo <- globaldom$col_lower_[col + 1L]
      if (lo != globaldom$col_upper_[col + 1L]) next
      if (lo != 1.0 && lo != 0.0) next
      fixval <- as.integer(lo)
      if (ct$numcliquesvar[.hmip_cv(col, 1L - fixval) + 1L] != 0L) {
        .hmip_HighsCliqueTable_vertexInfeasible(ct, globaldom, col, 1L - fixval)
        if (globaldom$infeasible_) return(invisible(NULL))
      }
    }
    start <- length(globaldom$domchgstack_column)
    .hmip_HighsDomain_propagate(globaldom)
    end <- length(globaldom$domchgstack_column)
  }
  invisible(NULL)
}

.hmip_HighsCliqueTable_vertexInfeasible <- function(ct, globaldom, col, val) {
  wasfixed <- .hmip_HighsDomain_isFixed(globaldom, col)
  .hmip_HighsDomain_fixCol(globaldom, col, as.double(1L - val))
  if (globaldom$infeasible_) return(invisible(NULL))
  if (!wasfixed) ct$nfixings <- ct$nfixings + 1L
  ct$infeasvertexstack <- c(ct$infeasvertexstack, .hmip_cv(col, val))
  .hmip_HighsCliqueTable_processInfeasibleVertices(ct, globaldom)
}

.hmip_HighsCliqueTable_addImplications <- function(ct, domain, col, val) {
  v <- .hmip_cv(col, val)
  reason <- .hmip_reason_cliqueTable(col, val)
  while (ct$colsubstituted[v %/% 2L + 1L] != 0L) {
    rep <- ct$subst_replace[ct$colsubstituted[v %/% 2L + 1L]]
    v <- if (v %% 2L == 1L) rep else .hmip_cv_complement(rep)
    if (v %% 2L == 1L) {
      if (domain$col_lower_[v %/% 2L + 1L] == 1.0) next
      .hmip_HighsDomain_changeBound3(domain, 0L, v %/% 2L, 1.0, reason)
      if (domain$infeasible_) return(invisible(NULL))
    } else {
      if (domain$col_upper_[v %/% 2L + 1L] == 0.0) next
      .hmip_HighsDomain_changeBound3(domain, 1L, v %/% 2L, 0.0, reason)
      if (domain$infeasible_) return(invisible(NULL))
    }
  }
  doFixings <- function(cliqueid) {
    start <- ct$c_start[cliqueid + 1L]; end <- ct$c_end[cliqueid + 1L]
    i <- start
    while (i != end) {
      e <- ct$cliqueentries[i + 1L]
      i <- i + 1L
      ecol <- e %/% 2L
      if (ecol == v %/% 2L) next
      if (e %% 2L == 1L) {
        if (domain$col_upper_[ecol + 1L] == 0.0) next
        .hmip_HighsDomain_changeBound3(domain, 1L, ecol, 0.0, reason)
        if (domain$infeasible_) return(TRUE)
      } else {
        if (domain$col_lower_[ecol + 1L] == 1.0) next
        .hmip_HighsDomain_changeBound3(domain, 0L, ecol, 1.0, reason)
        if (domain$infeasible_) return(TRUE)
      }
    }
    FALSE
  }
  for (id in .hmip_hht_keys(.hmip_cq_ihl(ct, v))) if (doFixings(id)) return(invisible(NULL))
  for (id in .hmip_hht_keys(.hmip_cq_ihl2(ct, v))) if (doFixings(id)) return(invisible(NULL))
  invisible(NULL)
}

.hmip_HighsCliqueTable_cleanupFixed <- function(ct, globaldom) {
  numcol <- length(globaldom$col_upper_)
  oldnfixings <- ct$nfixings
  for (i in seq_len(numcol) - 1L) {
    lo <- globaldom$col_lower_[i + 1L]
    if (ct$colDeleted[i + 1L] || lo != globaldom$col_upper_[i + 1L]) next
    if (lo != 1.0 && lo != 0.0) next
    fixval <- as.integer(lo)
    .hmip_HighsCliqueTable_vertexInfeasible(ct, globaldom, i, 1L - fixval)
    if (globaldom$infeasible_) return(invisible(NULL))
  }
  if (ct$nfixings != oldnfixings) .hmip_HighsCliqueTable_propagateAndCleanup(ct, globaldom)
  invisible(NULL)
}

.hmip_HighsCliqueTable_getNumImplications <- function(ct, col) {
  i0 <- .hmip_cv(col, 0L); i1 <- .hmip_cv(col, 1L)
  numimplics <- ct$numcliquesvar[i0 + 1L] + ct$numcliquesvar[i1 + 1L]
  for (v in c(i0, i1)) for (id in .hmip_hht_keys(.hmip_cq_ihl(ct, v))) {
    nimplics <- (ct$c_end[id + 1L] - ct$c_start[id + 1L] - 1L) * (1L + ct$c_equality[id + 1L])
    numimplics <- numimplics + nimplics - 1L
  }
  numimplics
}
.hmip_HighsCliqueTable_getNumImplications2 <- function(ct, col, val) {
  iv <- .hmip_cv(col, as.integer(val))
  numimplics <- ct$numcliquesvar[iv + 1L]
  for (id in .hmip_hht_keys(.hmip_cq_ihl(ct, iv))) {
    nimplics <- (ct$c_end[id + 1L] - ct$c_start[id + 1L] - 1L) * (1L + ct$c_equality[id + 1L])
    numimplics <- numimplics + nimplics - 1L
  }
  numimplics
}

.hmip_cq_ensure_iscandidate <- function(ct) {
  n <- length(ct$invertedHashList)
  if (length(ct$iscandidate) < n) ct$iscandidate <- c(ct$iscandidate, logical(n - length(ct$iscandidate)))
}

# the random extension step of both runCliqueMerging variants
.hmip_cq_random_extend <- function(ct, ext) {
  if (!length(ext)) return(ext)
  ext <- .hs_rng_shuffle(ct$randgen, ext)
  i <- 0L
  while (i < length(ext)) {
    extvar <- ext[i + 1L]
    i <- i + 1L
    tail <- if (i < length(ext)) ext[(i + 1L):length(ext)] else integer()
    r <- .hmip_HighsCliqueTable_shrinkToNeighbourhood(ct, extvar, tail)
    ext <- c(ext[seq_len(i)], r$q[seq_len(r$n)])
  }
  ext
}

.hmip_cq_not_fixed_zero <- function(globaldomain, vs) {
  vs[!vapply(vs, function(v) .hmip_HighsDomain_isFixed(globaldomain, v %/% 2L) &&
               trunc(globaldomain$col_lower_[v %/% 2L + 1L]) == 1L - v %% 2L, TRUE)]
}

# runCliqueMerging(globaldomain, clique, equation) -> the clique
.hmip_HighsCliqueTable_runCliqueMergingClique <- function(ct, globaldomain, clique, equation = FALSE) {
  numcliques <- .hmip_cq_IInf
  .hmip_cq_ensure_iscandidate(ct)
  initialCliqueSize <- length(clique)
  extensionstart <- NA_integer_
  for (i in seq_len(initialCliqueSize)) {
    # upstream tests cliqueentries[i], not clique[i]
    if (.hmip_HighsDomain_isFixed(globaldomain, ct$cliqueentries[i] %/% 2L)) next
    thisNumClqs <- ct$numcliquesvar[clique[i] + 1L]
    if (thisNumClqs < numcliques) {
      numcliques <- thisNumClqs
      extensionstart <- clique[i]
    }
  }
  if (numcliques == .hmip_cq_IInf) return(integer())
  ct$iscandidate[clique[seq_len(initialCliqueSize)] + 1L] <- TRUE
  addCands <- function(cliqueid) {
    i <- ct$c_start[cliqueid + 1L]
    while (i != ct$c_end[cliqueid + 1L]) {
      e <- ct$cliqueentries[i + 1L]
      i <- i + 1L
      if (ct$iscandidate[e + 1L] || .hmip_HighsDomain_isFixed(globaldomain, e %/% 2L)) next
      ct$iscandidate[e + 1L] <- TRUE
      clique <<- c(clique, e)
    }
  }
  for (id in .hmip_hht_keys(.hmip_cq_ihl(ct, extensionstart))) addCands(id)
  for (id in .hmip_hht_keys(.hmip_cq_ihl2(ct, extensionstart))) addCands(id)
  ct$iscandidate[clique + 1L] <- FALSE
  ext <- clique[-seq_len(initialCliqueSize)]
  for (i in seq_len(initialCliqueSize)) {
    if (!length(ext)) break
    if (clique[i] == extensionstart) next
    r <- .hmip_HighsCliqueTable_shrinkToNeighbourhood(ct, clique[i], ext)
    ext <- r$q[seq_len(r$n)]
  }
  ext <- .hmip_cq_random_extend(ct, ext)
  clique <- c(clique[seq_len(initialCliqueSize)], ext)
  if (equation) {
    for (x in ext) .hmip_HighsCliqueTable_vertexInfeasible(ct, globaldomain, x %/% 2L, x %% 2L)
  } else {
    clique <- .hmip_HighsCliqueTable_runCliqueSubsumption(ct, globaldomain, clique)$clique
    if (length(clique)) clique <- .hmip_cq_not_fixed_zero(globaldomain, clique)
  }
  .hmip_HighsCliqueTable_processInfeasibleVertices(ct, globaldomain)
  clique
}

.hmip_HighsCliqueTable_runCliqueMerging <- function(ct, globaldomain) {
  .hmip_cq_ensure_iscandidate(ct)
  nc <- length(ct$c_start)
  if (length(ct$cliquehits) < nc) ct$cliquehits <- c(ct$cliquehits, integer(nc - length(ct$cliquehits)))
  numcliqueslots <- nc
  maxNewEntries <- ct$numEntries + .hmip_HighsDomain_numModelNonzeros(globaldomain)
  haveNonModelCliquesToMerge <- FALSE
  for (k in seq_len(numcliqueslots) - 1L) {
    k1 <- k + 1L
    if (ct$c_start[k1] == -1L) next
    if (!ct$c_equality[k1] && ct$c_origin[k1] == .hmip_cq_IInf) next
    if (ct$c_origin[k1] == -1L) { haveNonModelCliquesToMerge <- TRUE; next }
    numclqvars <- ct$c_end[k1] - ct$c_start[k1]
    if (numclqvars == 0L) next
    clqvars <- ct$cliqueentries[(ct$c_start[k1] + 1L):ct$c_end[k1]]
    extensionstart <- clqvars[1L]
    numcliques <- ct$numcliquesvar[clqvars[1L] + 1L]
    for (i in seq_len(numclqvars - 1L) + 1L) {
      thisNumClqs <- ct$numcliquesvar[clqvars[i] + 1L]
      if (thisNumClqs < numcliques) { numcliques <- thisNumClqs; extensionstart <- clqvars[i] }
    }
    ct$iscandidate[clqvars + 1L] <- TRUE
    extensionvars <- integer()
    addCands <- function(cliqueid) {
      i <- ct$c_start[cliqueid + 1L]
      while (i != ct$c_end[cliqueid + 1L]) {
        e <- ct$cliqueentries[i + 1L]
        i <- i + 1L
        if (ct$iscandidate[e + 1L] || .hmip_HighsDomain_isFixed(globaldomain, e %/% 2L)) next
        ct$iscandidate[e + 1L] <- TRUE
        extensionvars <<- c(extensionvars, e)
      }
    }
    for (id in .hmip_hht_keys(.hmip_cq_ihl(ct, extensionstart))) addCands(id)
    for (id in .hmip_hht_keys(.hmip_cq_ihl2(ct, extensionstart))) addCands(id)
    ct$iscandidate[clqvars + 1L] <- FALSE
    ct$iscandidate[extensionvars + 1L] <- FALSE
    for (x in clqvars) {
      if (!length(extensionvars)) break
      if (x == extensionstart) next
      r <- .hmip_HighsCliqueTable_shrinkToNeighbourhood(ct, x, extensionvars)
      extensionvars <- r$q[seq_len(r$n)]
    }
    extensionvars <- .hmip_cq_random_extend(ct, extensionvars)
    if (ct$c_equality[k1]) {
      for (v in extensionvars) .hmip_HighsCliqueTable_vertexInfeasible(ct, globaldomain, v %/% 2L, v %% 2L)
    } else {
      originrow <- ct$c_origin[k1]
      ct$c_origin[k1] <- .hmip_cq_IInf
      numExtensions <- length(extensionvars)
      tailv <- ct$cliqueentries[(ct$c_start[k1] + 1L):ct$c_end[k1]]
      tailv <- tailv[!ct$colDeleted[tailv %/% 2L + 1L]]
      extensionvars <- c(extensionvars, tailv)
      .hmip_HighsCliqueTable_removeClique(ct, k)
      for (v in extensionvars) {
        for (id in c(.hmip_hht_keys(.hmip_cq_ihl(ct, v)), .hmip_hht_keys(.hmip_cq_ihl2(ct, v)))) {
          if (ct$cliquehits[id + 1L] == 0L) ct$cliquehitinds <- c(ct$cliquehitinds, id)
          ct$cliquehits[id + 1L] <- ct$cliquehits[id + 1L] + 1L
        }
      }
      redundant <- FALSE
      dominatingOrigin <- .hmip_cq_IInf
      for (cliqueid in ct$cliquehitinds) {
        c1 <- cliqueid + 1L
        hits <- ct$cliquehits[c1]
        ct$cliquehits[c1] <- 0L
        if (hits == length(extensionvars)) {
          redundant <- TRUE
          if (ct$c_origin[c1] != .hmip_cq_IInf && ct$c_origin[c1] != -1L) dominatingOrigin <- ct$c_origin[c1]
        } else if (ct$c_end[c1] - ct$c_start[c1] - ct$c_numZeroFixed[c1] == hits) {
          if (ct$c_equality[c1]) {
            for (v in extensionvars) {
              sizeTwo <- ct$c_end[c1] - ct$c_start[c1] == 2L
              vHasClq <- if (sizeTwo) .hmip_hht_contains(.hmip_cq_ihl2(ct, v), cliqueid)
                else .hmip_hht_contains(.hmip_cq_ihl(ct, v), cliqueid)
              if (!vHasClq) ct$infeasvertexstack <- c(ct$infeasvertexstack, v)
            }
          } else .hmip_HighsCliqueTable_removeClique(ct, cliqueid)
        }
      }
      ct$cliquehitinds <- integer()
      if (!redundant) {
        for (i in seq_len(numExtensions)) {
          ct$ext_row <- c(ct$ext_row, originrow); ct$ext_var <- c(ct$ext_var, extensionvars[i])
        }
        extensionvars <- .hmip_cq_not_fixed_zero(globaldomain, extensionvars)
        if (length(extensionvars) > 1L)
          .hmip_HighsCliqueTable_doAddClique(ct, extensionvars, length(extensionvars), FALSE, originrow)
      } else {
        if (dominatingOrigin != .hmip_cq_IInf) ct$deletedrows <- c(ct$deletedrows, originrow)
        else for (i in seq_len(numExtensions)) {
          ct$ext_row <- c(ct$ext_row, originrow); ct$ext_var <- c(ct$ext_var, extensionvars[i])
        }
      }
    }
    .hmip_HighsCliqueTable_processInfeasibleVertices(ct, globaldomain)
    if (ct$numEntries >= maxNewEntries) break
  }
  if (haveNonModelCliquesToMerge) {
    for (k in seq_len(numcliqueslots) - 1L) {
      k1 <- k + 1L
      if (ct$c_start[k1] == -1L) next
      if (ct$c_origin[k1] != -1L) next
      extensionvars <- ct$cliqueentries[(ct$c_start[k1] + 1L):ct$c_end[k1]]
      .hmip_HighsCliqueTable_removeClique(ct, k)
      extensionvars <- .hmip_HighsCliqueTable_runCliqueMergingClique(ct, globaldomain, extensionvars)
      if (length(extensionvars) > 1L) .hmip_HighsCliqueTable_doAddClique(ct, extensionvars, length(extensionvars))
    }
  }
  invisible(NULL)
}

# *this = std::move(newCliqueTable)
.hmip_cq_replace_with <- function(ct, newct) {
  for (nm in ls(ct, all.names = TRUE)) rm(list = nm, envir = ct)
  for (nm in ls(newct, all.names = TRUE)) assign(nm, get(nm, envir = newct), envir = ct)
  invisible(ct)
}

.hmip_HighsCliqueTable_rebuild <- function(ct, ncols, postSolveStack, globaldomain, orig2reducedcol, orig2reducedrow) {
  newct <- .hmip_new_HighsCliqueTable(ncols)
  .hmip_HighsCliqueTable_setPresolveFlag(newct, ct$inPresolve)
  .hmip_HighsCliqueTable_setMinEntriesForParallelism(newct, ct$minEntriesForParallelism)
  for (i in seq_along(ct$c_start)) {
    start <- ct$c_start[i]
    if (start == -1L) next
    idx <- (start + 1L):ct$c_end[i]
    newcols <- vapply(ct$cliqueentries[idx], function(v) {
      col <- orig2reducedcol[v %/% 2L + 1L]
      if (col == -1L || !.hmip_HighsDomain_isBinary(globaldomain, col) ||
          !.hmip_HighsPostsolveStack_isColLinearlyTransformable(postSolveStack, col)) .hmip_cq_IInf else col
    }, 0L)
    keep <- newcols != .hmip_cq_IInf
    numvars <- sum(keep)
    kept <- as.integer(2L * newcols[keep] + ct$cliqueentries[idx][keep] %% 2L)
    ct$cliqueentries[idx[seq_len(numvars)]] <- kept
    if (numvars <= 1L) next
    origin <- if (ct$c_origin[i] != .hmip_cq_IInf) -1L else .hmip_cq_IInf
    .hmip_HighsCliqueTable_doAddClique(newct, kept, numvars, FALSE, origin)
  }
  .hmip_cq_replace_with(ct, newct)
}

.hmip_HighsCliqueTable_buildFrom <- function(ct, origModel, init) {
  ncols <- length(init$colsubstituted)
  newct <- .hmip_new_HighsCliqueTable(ncols)
  .hmip_HighsCliqueTable_setPresolveFlag(newct, ct$inPresolve)
  # upstream passes minEntriesForParallelism to setPresolveFlag
  .hmip_HighsCliqueTable_setPresolveFlag(newct, ct$minEntriesForParallelism != 0L)
  for (i in seq_along(init$c_start)) {
    start <- init$c_start[i]
    if (start == -1L) next
    numvars <- init$c_end[i] - start
    if (numvars - init$c_numZeroFixed[i] <= 1L) next
    clq <- init$cliqueentries[(start + 1L):init$c_end[i]]
    clq <- clq[!(origModel$col_lower_[clq %/% 2L + 1L] != 0.0 | origModel$col_upper_[clq %/% 2L + 1L] != 1.0)]
    if (length(clq) <= 1L) next
    origin <- if (init$c_origin[i] != .hmip_cq_IInf) -1L else .hmip_cq_IInf
    .hmip_HighsCliqueTable_doAddClique(newct, clq, length(clq), FALSE, origin)
  }
  newct$colsubstituted <- init$colsubstituted
  newct$subst_substcol <- init$subst_substcol
  newct$subst_replace <- init$subst_replace
  .hmip_cq_replace_with(ct, newct)
}

# ---------------------------------------------------------------------------
# HighsImplications
#
# vubs[[col + 1L]] and vlbs[[col + 1L]] are HighsHashTree<HighsInt, VarBound>
# trees keyed by the binary column with the VarBound c(coef, constant).
# implics[[2 * col + val + 1L]] is the SoA domain change list of the stored
# implications, implic_computed the computed flags. substitutions are the
# vectors subst_substcol, subst_staycol, subst_scale, subst_offset.
# ---------------------------------------------------------------------------

.hmip_im_vbMin <- function(vb) .cd_val(.cd_add(c(vb[2L], 0), min(vb[1L], 0.0)))
.hmip_im_vbMax <- function(vb) .cd_val(.cd_add(c(vb[2L], 0), max(vb[1L], 0.0)))
.hmip_HighsImplications_calcMaxVarBounds <- function(numcol) 5000000 + 10 * numcol

.hmip_im_emptychg <- function() list(boundval = numeric(), column = integer(), boundtype = integer())

.hmip_HighsImplications_reset <- function(im) {
  numcol <- .hmip_cq_numCol(im$mipsolver)
  im$implics <- vector("list", 2L * numcol)
  im$implic_computed <- logical(2L * numcol)
  im$colsubstituted <- logical(numcol)
  im$numImplications <- 0
  im$vubs <- vector("list", numcol)
  im$vlbs <- vector("list", numcol)
  im$numVarBounds <- 0
  im$maxVarBounds <- .hmip_HighsImplications_calcMaxVarBounds(numcol)
  im$nextCleanupCall <- .hmip_cq_numNonzero(im$mipsolver)
  invisible(NULL)
}

.hmip_new_HighsImplications <- function(mipsolver) {
  im <- new.env(parent = emptyenv())
  im$mipsolver <- mipsolver
  im$subst_substcol <- integer(); im$subst_staycol <- integer()
  im$subst_scale <- numeric(); im$subst_offset <- numeric()
  im$storeLiftingOpportunity <- NULL
  .hmip_HighsImplications_reset(im)
  im
}

.hmip_HighsImplications_getNumImplications <- function(im) as.integer(im$numImplications)
.hmip_HighsImplications_tooManyVarBounds <- function(im) im$numVarBounds >= im$maxVarBounds
.hmip_HighsImplications_implicationsCached <- function(im, col, val) im$implic_computed[2L * col + val + 1L]
.hmip_HighsImplications_getVubs <- function(im, col) im$vubs[[col + 1L]]
.hmip_HighsImplications_getVlbs <- function(im, col) im$vlbs[[col + 1L]]

# getImplications(col, val, infeasible) -> list(implics (SoA), infeasible)
.hmip_HighsImplications_getImplications <- function(im, col, val) {
  loc <- 2L * col + as.integer(val) + 1L
  infeasible <- if (!im$implic_computed[loc]) .hmip_HighsImplications_computeImplications(im, col, val) else FALSE
  imp <- im$implics[[loc]]
  list(implics = if (is.null(imp)) .hmip_im_emptychg() else imp, infeasible = infeasible)
}

.hmip_HighsImplications_computeImplications <- function(im, col, val) {
  ms <- im$mipsolver; md <- ms$mipdata_
  gd <- md$domain
  ct <- md$cliquetable
  val <- as.integer(val)
  .hmip_HighsDomain_propagate(gd)
  if (gd$infeasible_ || .hmip_HighsDomain_isFixed(gd, col)) return(TRUE)
  if (!is.null(im$storeLiftingOpportunity)) .hmip_HighsDomain_setRecordRedundantRows(gd, TRUE)
  changedend <- length(gd$changedcols_)
  stackimplicstart <- length(gd$domchgstack_column) + 1L
  numImplications <- -stackimplicstart
  if (val == 1L) .hmip_HighsDomain_changeBound3(gd, 0L, col, 1)
  else .hmip_HighsDomain_changeBound3(gd, 1L, col, 0)
  storeLifting <- function() {
    f <- im$storeLiftingOpportunity
    if (!is.null(f)) {
      for (row in .hmip_HighsDomain_getRedundantRows(gd))
        f(row, col, if (val == 1L) 1L else 0L, (if (val == 1L) -1 else 1) * .hmip_HighsDomain_getRedundantRowValue(gd, row))
      .hmip_HighsDomain_clearRedundantRows(gd)
      .hmip_HighsDomain_setRecordRedundantRows(gd, FALSE)
    }
  }
  doBacktrack <- function() {
    .hmip_HighsDomain_backtrack(gd)
    .hmip_HighsDomain_clearChangedCols(gd, changedend)
  }
  isInfeasible <- function() {
    if (!gd$infeasible_) return(FALSE)
    storeLifting()
    doBacktrack()
    .hmip_HighsCliqueTable_vertexInfeasible(ct, gd, col, val)
    TRUE
  }
  if (isInfeasible()) return(TRUE)
  .hmip_HighsDomain_propagate(gd)
  if (isInfeasible()) return(TRUE)
  stackimplicend <- length(gd$domchgstack_column)
  numImplications <- numImplications + stackimplicend
  .hmip_HighsPseudocost_addInferenceObservation(md$pseudocost, col, numImplications, val == 1L)
  numEntries <- ct$numEntries
  maxEntries <- 100000 + .hmip_cq_numNonzero(ms)
  sel <- integer()
  i <- stackimplicstart
  while (i < stackimplicend) {
    if (gd$domchgreason_type[i + 1L] == -5L &&
        (bitwShiftR(gd$domchgreason_index[i + 1L], 1L) == col || numEntries >= maxEntries)) { i <- i + 1L; next }
    sel <- c(sel, i)
    i <- i + 1L
  }
  ib <- gd$domchgstack_boundval[sel + 1L]
  ic <- gd$domchgstack_column[sel + 1L]
  it <- gd$domchgstack_boundtype[sel + 1L]
  storeLifting()
  doBacktrack()
  pr <- .hmip_std_partition(seq_along(ic), function(k) !.hmip_HighsDomain_isBinary(gd, ic[k]))
  o <- pr$v
  ib <- ib[o]; ic <- ic[o]; it <- it[o]
  nb <- pr$k
  if (nb > 0L) {
    s <- seq_len(nb)
    oo <- order(ic[s], it[s], ib[s])
    ib[s] <- ib[s][oo]; ic[s] <- ic[s][oo]; it[s] <- it[s][oo]
  }
  cl0 <- .hmip_cv(col, val)
  for (k in seq_len(length(ic) - nb) + nb) {
    cl1 <- if (it[k] == 0L) .hmip_cv(ic[k], 0L) else .hmip_cv(ic[k], 1L)
    res <- .hmip_HighsCliqueTable_addClique(ct, ms, c(cl0, cl1), 2L)
    cl0 <- res[1L]
    if (gd$infeasible_ || .hmip_HighsDomain_isFixed(gd, col)) return(TRUE)
  }
  for (k in seq_len(nb)) {
    c2 <- ic[k]; bv <- ib[k]
    if (it[k] == 0L) {
      lo <- gd$col_lower_[c2 + 1L]
      if (val == 1L) {
        if (lo != -Inf) .hmip_HighsImplications_addVLB(im, c2, col, bv - lo, lo)
      } else .hmip_HighsImplications_addVLB(im, c2, col, lo - bv, bv)
    } else {
      up <- gd$col_upper_[c2 + 1L]
      if (val == 1L) {
        if (up != Inf) .hmip_HighsImplications_addVUB(im, c2, col, bv - up, up)
      } else .hmip_HighsImplications_addVUB(im, c2, col, up - bv, bv)
    }
  }
  loc <- 2L * col + val + 1L
  im$implic_computed[loc] <- TRUE
  if (nb > 0L) {
    s <- seq_len(nb)
    im$implics[[loc]] <- list(boundval = ib[s], column = ic[s], boundtype = it[s])
    im$numImplications <- im$numImplications + nb
  }
  FALSE
}

# getBestVub(col, lpSolution, bestUb) -> list(first = vubCol, second =
# VarBound, bestUb); lpSolution has col_value and col_dual
.hmip_HighsImplications_getBestVub <- function(im, col, lpSolution, bestUb) {
  md <- im$mipsolver$mipdata_; gd <- md$domain; ft <- md$feastol
  best_col <- -1L; best_vb <- c(0.0, Inf)
  minbestUb <- bestUb; bestUbDist <- Inf; bestvubnodes <- 0
  xv <- lpSolution$col_value; xd <- lpSolution$col_dual
  scale <- gd$col_upper_[col + 1L] - gd$col_lower_[col + 1L]
  scale <- if (scale == Inf) 1.0 else 1.0 / scale
  e <- .hmip_hht_entries(im$vubs[[col + 1L]])
  for (i in seq_along(e$k)) {
    vubCol <- e$k[i]; vub <- e$v[[i]]
    if (vub[1L] == Inf) next
    if (.hmip_HighsDomain_isFixed(gd, vubCol)) next
    vubval <- xv[vubCol + 1L] * vub[1L] + vub[2L]
    ubDist <- max(0.0, vubval - xv[col + 1L])
    yDist <- ft + (if (vub[1L] > 0) 1 - xv[vubCol + 1L] else xv[vubCol + 1L])
    norm2 <- 1.0 + vub[1L] * vub[1L]
    if (ubDist * ubDist > yDist * yDist * norm2) next
    ubDist <- ubDist * scale
    if (ubDist <= bestUbDist + ft) {
      minvubval <- .hmip_im_vbMin(vub)
      vubnodes <- if (vub[1L] > 0) .hmip_HighsNodeQueue_numNodesDown(md$nodequeue, vubCol)
        else .hmip_HighsNodeQueue_numNodesUp(md$nodequeue, vubCol)
      better <- (function() {
        if (ubDist < bestUbDist - ft) return(TRUE)
        if (vubnodes > bestvubnodes) return(TRUE)
        if (vubnodes < bestvubnodes) return(FALSE)
        if (minvubval < minbestUb - ft) return(TRUE)
        if (minvubval > minbestUb + ft) return(FALSE)
        xd[vubCol + 1L] / vub[1L] - xd[best_col + 1L] / best_vb[1L] > ft
      })()
      if (isTRUE(better)) {
        bestUb <- vubval; minbestUb <- minvubval
        best_col <- vubCol; best_vb <- vub
        bestvubnodes <- vubnodes; bestUbDist <- ubDist
      }
    }
  }
  list(first = best_col, second = list(coef = best_vb[1L], constant = best_vb[2L]), bestUb = bestUb)
}

.hmip_HighsImplications_getBestVlb <- function(im, col, lpSolution, bestLb) {
  md <- im$mipsolver$mipdata_; gd <- md$domain; ft <- md$feastol
  best_col <- -1L; best_vb <- c(0.0, -Inf)
  maxbestlb <- bestLb; bestLbDist <- Inf; bestvlbnodes <- 0
  xv <- lpSolution$col_value; xd <- lpSolution$col_dual
  scale <- gd$col_upper_[col + 1L] - gd$col_lower_[col + 1L]
  scale <- if (scale == Inf) 1.0 else 1.0 / scale
  e <- .hmip_hht_entries(im$vlbs[[col + 1L]])
  for (i in seq_along(e$k)) {
    vlbCol <- e$k[i]; vlb <- e$v[[i]]
    if (vlb[1L] == -Inf) next
    if (.hmip_HighsDomain_isFixed(gd, vlbCol)) next
    vlbval <- xv[vlbCol + 1L] * vlb[1L] + vlb[2L]
    lbDist <- max(0.0, xv[col + 1L] - vlbval)
    yDist <- ft + (if (vlb[1L] > 0) xv[vlbCol + 1L] else 1 - xv[vlbCol + 1L])
    norm2 <- 1.0 + vlb[1L] * vlb[1L]
    if (lbDist * lbDist > yDist * yDist * norm2) next
    lbDist <- lbDist * scale
    if (lbDist <= bestLbDist + ft) {
      maxvlbval <- .hmip_im_vbMax(vlb)
      vlbnodes <- if (vlb[1L] > 0) .hmip_HighsNodeQueue_numNodesUp(md$nodequeue, vlbCol)
        else .hmip_HighsNodeQueue_numNodesDown(md$nodequeue, vlbCol)
      better <- (function() {
        if (lbDist < bestLbDist - ft) return(TRUE)
        if (vlbnodes > bestvlbnodes) return(TRUE)
        if (vlbnodes < bestvlbnodes) return(FALSE)
        if (maxvlbval > maxbestlb + ft) return(TRUE)
        if (maxvlbval < maxbestlb - ft) return(FALSE)
        xd[vlbCol + 1L] / vlb[1L] - xd[best_col + 1L] / best_vb[1L] < -ft
      })()
      if (isTRUE(better)) {
        bestLb <- vlbval; maxbestlb <- maxvlbval
        best_col <- vlbCol; best_vb <- vlb
        bestvlbnodes <- vlbnodes; bestLbDist <- lbDist
      }
    }
  }
  list(first = best_col, second = list(coef = best_vb[1L], constant = best_vb[2L]), bestLb = bestLb)
}

# runProbing(col, numReductions) -> list(result, numReductions)
.hmip_HighsImplications_runProbing <- function(im, col, numReductions) {
  md <- im$mipsolver$mipdata_; gd <- md$domain; ct <- md$cliquetable
  if (.hmip_HighsDomain_isBinary(gd, col) && !.hmip_HighsImplications_implicationsCached(im, col, 1L) &&
      !.hmip_HighsImplications_implicationsCached(im, col, 0L) &&
      is.null(.hmip_HighsCliqueTable_getSubstitution(ct, col))) {
    out <- function() list(result = TRUE, numReductions = numReductions)
    infeasible <- .hmip_HighsImplications_computeImplications(im, col, 1L)
    if (gd$infeasible_ || infeasible || !is.null(.hmip_HighsCliqueTable_getSubstitution(ct, col))) return(out())
    infeasible <- .hmip_HighsImplications_computeImplications(im, col, 0L)
    if (gd$infeasible_ || infeasible || !is.null(.hmip_HighsCliqueTable_getSubstitution(ct, col))) return(out())
    dn <- .hmip_HighsImplications_getImplications(im, col, 0L)$implics
    up <- .hmip_HighsImplications_getImplications(im, col, 1L)$implics
    nd <- length(dn$column); nu <- length(up$column)
    u <- 0L; d <- 0L
    while (u < nu && d < nd) {
      if (up$column[u + 1L] < dn$column[d + 1L]) u <- u + 1L
      else if (dn$column[d + 1L] < up$column[u + 1L]) d <- d + 1L
      else {
        implcol <- up$column[u + 1L]
        lbDown <- gd$col_lower_[implcol + 1L]; ubDown <- gd$col_upper_[implcol + 1L]
        lbUp <- lbDown; ubUp <- ubDown
        repeat {
          if (dn$boundtype[d + 1L] == 0L) lbDown <- max(lbDown, dn$boundval[d + 1L])
          else ubDown <- min(ubDown, dn$boundval[d + 1L])
          d <- d + 1L
          if (!(d < nd && dn$column[d + 1L] == implcol)) break
        }
        repeat {
          if (up$boundtype[u + 1L] == 0L) lbUp <- max(lbUp, up$boundval[u + 1L])
          else ubUp <- min(ubUp, up$boundval[u + 1L])
          u <- u + 1L
          if (!(u < nu && up$column[u + 1L] == implcol)) break
        }
        if (im$colsubstituted[implcol + 1L] || .hmip_HighsDomain_isFixed(gd, implcol)) next
        if (lbDown == ubDown && lbUp == ubUp && abs(lbDown - lbUp) > md$feastol) {
          im$subst_substcol <- c(im$subst_substcol, implcol)
          im$subst_staycol <- c(im$subst_staycol, col)
          im$subst_offset <- c(im$subst_offset, lbDown)
          im$subst_scale <- c(im$subst_scale, lbUp - lbDown)
          im$colsubstituted[implcol + 1L] <- TRUE
          numReductions <- numReductions + 1L
        } else {
          lb <- min(lbDown, lbUp); ub <- max(ubDown, ubUp)
          if (lb > gd$col_lower_[implcol + 1L]) {
            .hmip_HighsDomain_changeBound3(gd, 0L, implcol, lb, .hmip_reason_unspecified())
            numReductions <- numReductions + 1L
          }
          if (ub < gd$col_upper_[implcol + 1L]) {
            .hmip_HighsDomain_changeBound3(gd, 1L, implcol, ub, .hmip_reason_unspecified())
            numReductions <- numReductions + 1L
          }
        }
      }
    }
    return(out())
  }
  list(result = FALSE, numReductions = numReductions)
}

.hmip_HighsImplications_addVUB <- function(im, col, vubcol, vubcoef, vubconstant) {
  if (.hmip_HighsImplications_tooManyVarBounds(im)) return(invisible(NULL))
  md <- im$mipsolver$mipdata_
  vub <- c(vubcoef, vubconstant)
  minBound <- .hmip_im_vbMin(vub)
  if (minBound >= md$domain$col_upper_[col + 1L] - md$feastol) return(invisible(NULL))
  tree <- im$vubs[[col + 1L]]
  cur <- .hmip_hht_find(tree, vubcol)
  if (is.null(cur)) {
    im$vubs[col + 1L] <- list(.hmip_hht_insert(tree, vubcol, vub)$node)
    im$numVarBounds <- im$numVarBounds + 1
  } else if (minBound < .hmip_im_vbMin(cur$value) - md$feastol) {
    im$vubs[col + 1L] <- list(.hmip_hht_set_value(tree, vubcol, vub))
  }
  invisible(NULL)
}

.hmip_HighsImplications_addVLB <- function(im, col, vlbcol, vlbcoef, vlbconstant) {
  if (.hmip_HighsImplications_tooManyVarBounds(im)) return(invisible(NULL))
  md <- im$mipsolver$mipdata_
  vlb <- c(vlbcoef, vlbconstant)
  maxBound <- .hmip_im_vbMax(vlb)
  if (maxBound <= md$domain$col_lower_[col + 1L] + md$feastol) return(invisible(NULL))
  tree <- im$vlbs[[col + 1L]]
  cur <- .hmip_hht_find(tree, vlbcol)
  if (is.null(cur)) {
    im$vlbs[col + 1L] <- list(.hmip_hht_insert(tree, vlbcol, vlb)$node)
    im$numVarBounds <- im$numVarBounds + 1
  } else if (maxBound > .hmip_im_vbMax(cur$value) + md$feastol) {
    im$vlbs[col + 1L] <- list(.hmip_hht_set_value(tree, vlbcol, vlb))
  }
  invisible(NULL)
}

.hmip_HighsImplications_columnTransformed <- function(im, col, scale, constant) {
  if (scale < 0) {
    t <- im$vubs[[col + 1L]]
    im$vubs[col + 1L] <- list(im$vlbs[[col + 1L]])
    im$vlbs[col + 1L] <- list(t)
  }
  f <- function(k, vb) {
    vb[2L] <- vb[2L] - constant
    vb[2L] <- vb[2L] / scale
    vb[1L] <- vb[1L] / scale
    vb
  }
  im$vlbs[col + 1L] <- list(.hmip_hht_map_values(im$vlbs[[col + 1L]], f))
  im$vubs[col + 1L] <- list(.hmip_hht_map_values(im$vubs[[col + 1L]], f))
  for (i in seq_along(im$subst_substcol)) {
    if (im$subst_substcol[i] == col) {
      im$subst_offset[i] <- im$subst_offset[i] - constant
      im$subst_offset[i] <- im$subst_offset[i] / scale
      im$subst_scale[i] <- im$subst_scale[i] / scale
    }
  }
  invisible(NULL)
}

.hmip_HighsImplications_rebuild <- function(im, ncols, orig2reducedcol, orig2reducedrow) {
  md <- im$mipsolver$mipdata_
  oldvubs <- im$vubs; oldvlbs <- im$vlbs
  im$implics <- vector("list", 2L * ncols)
  im$implic_computed <- logical(2L * ncols)
  im$colsubstituted <- logical(ncols)
  im$subst_substcol <- integer(); im$subst_staycol <- integer()
  im$subst_scale <- numeric(); im$subst_offset <- numeric()
  im$vubs <- vector("list", ncols); im$vlbs <- vector("list", ncols)
  im$numImplications <- 0
  im$numVarBounds <- 0
  im$nextCleanupCall <- .hmip_cq_numNonzero(im$mipsolver)
  pss <- md$postSolveStack
  for (i in seq_along(oldvubs) - 1L) {
    newi <- orig2reducedcol[i + 1L]
    if (newi == -1L || !.hmip_HighsPostsolveStack_isColLinearlyTransformable(pss, newi)) next
    for (side in 1:2) {
      e <- .hmip_hht_entries(if (side == 1L) oldvubs[[i + 1L]] else oldvlbs[[i + 1L]])
      for (j in seq_along(e$k)) {
        nc <- orig2reducedcol[e$k[j] + 1L]
        if (nc == -1L) next
        if (!.hmip_HighsDomain_isBinary(md$domain, nc) || !.hmip_HighsPostsolveStack_isColLinearlyTransformable(pss, nc)) next
        vb <- e$v[[j]]
        if (side == 1L) .hmip_HighsImplications_addVUB(im, newi, nc, vb[1L], vb[2L])
        else .hmip_HighsImplications_addVLB(im, newi, nc, vb[1L], vb[2L])
      }
    }
  }
  invisible(NULL)
}

.hmip_HighsImplications_buildFrom <- function(im, init) {
  gd <- im$mipsolver$mipdata_$domain
  for (i in seq_len(.hmip_cq_numCol(im$mipsolver)) - 1L) {
    e <- .hmip_hht_entries(init$vubs[[i + 1L]])
    for (j in seq_along(e$k)) {
      if (!.hmip_HighsDomain_isBinary(gd, e$k[j])) next
      .hmip_HighsImplications_addVUB(im, i, e$k[j], e$v[[j]][1L], e$v[[j]][2L])
    }
    e <- .hmip_hht_entries(init$vlbs[[i + 1L]])
    for (j in seq_along(e$k)) {
      if (!.hmip_HighsDomain_isBinary(gd, e$k[j])) next
      .hmip_HighsImplications_addVLB(im, i, e$k[j], e$v[[j]][1L], e$v[[j]][2L])
    }
  }
  invisible(NULL)
}

# cleanupVlb / cleanupVub -> list(vb, redundant, infeasible)
.hmip_HighsImplications_cleanupVlb <- function(im, col, vlbCol, vlb, lb, allowBoundChanges = TRUE) {
  vlb <- c(vlb$coef, vlb$constant)
  redundant <- FALSE; infeasible <- FALSE
  if (vlbCol == -1L) return(list(vlb = list(coef = vlb[1L], constant = vlb[2L]), redundant = redundant, infeasible = infeasible))
  md <- im$mipsolver$mipdata_
  maxlb <- .hmip_im_vbMax(vlb); minlb <- .hmip_im_vbMin(vlb)
  if (maxlb <= lb + md$feastol) {
    redundant <- TRUE
  } else if (minlb < lb - md$epsilon) {
    newcoef <- lb - maxlb
    if (vlb[1L] < 0) vlb[1L] <- newcoef
    else { vlb[2L] <- lb; vlb[1L] <- -newcoef }
  } else if (allowBoundChanges && minlb > lb + md$epsilon) {
    .hmip_HighsDomain_changeBound3(md$domain, 0L, col, minlb, .hmip_reason_unspecified())
    infeasible <- md$domain$infeasible_
  }
  list(vlb = list(coef = vlb[1L], constant = vlb[2L]), redundant = redundant, infeasible = infeasible)
}
.hmip_HighsImplications_cleanupVub <- function(im, col, vubCol, vub, ub, allowBoundChanges = TRUE) {
  vub <- c(vub$coef, vub$constant)
  redundant <- FALSE; infeasible <- FALSE
  if (vubCol == -1L) return(list(vub = list(coef = vub[1L], constant = vub[2L]), redundant = redundant, infeasible = infeasible))
  md <- im$mipsolver$mipdata_
  maxub <- .hmip_im_vbMax(vub); minub <- .hmip_im_vbMin(vub)
  if (minub >= ub - md$feastol) {
    redundant <- TRUE
  } else if (maxub > ub + md$epsilon) {
    newcoef <- ub - minub
    if (vub[1L] > 0) vub[1L] <- newcoef
    else { vub[2L] <- ub; vub[1L] <- -newcoef }
  } else if (allowBoundChanges && maxub < ub - md$epsilon) {
    .hmip_HighsDomain_changeBound3(md$domain, 1L, col, maxub, .hmip_reason_unspecified())
    infeasible <- md$domain$infeasible_
  }
  list(vub = list(coef = vub[1L], constant = vub[2L]), redundant = redundant, infeasible = infeasible)
}

.hmip_HighsImplications_cleanupVarbounds <- function(im, col) {
  gd <- im$mipsolver$mipdata_$domain
  ub <- gd$col_upper_[col + 1L]; lb <- gd$col_lower_[col + 1L]
  if (ub == lb) {
    # upstream subtracts the number of columns, not the number of bounds
    im$numVarBounds <- im$numVarBounds - length(im$vlbs)
    im$numVarBounds <- im$numVarBounds - length(im$vubs)
    im$vlbs[col + 1L] <- list(NULL)
    im$vubs[col + 1L] <- list(NULL)
    return(invisible(NULL))
  }
  delVbds <- integer()
  im$vubs[col + 1L] <- list(.hmip_hht_map_values(im$vubs[[col + 1L]], function(k, vb) {
    r <- .hmip_HighsImplications_cleanupVub(im, col, k, list(coef = vb[1L], constant = vb[2L]), ub)
    if (r$redundant) delVbds <<- c(delVbds, k)
    c(r$vub$coef, r$vub$constant)
  }))
  for (k in delVbds) im$vubs[col + 1L] <- list(.hmip_hht_erase(im$vubs[[col + 1L]], k))
  im$numVarBounds <- im$numVarBounds - length(delVbds)
  delVbds <- integer()
  im$vlbs[col + 1L] <- list(.hmip_hht_map_values(im$vlbs[[col + 1L]], function(k, vb) {
    r <- .hmip_HighsImplications_cleanupVlb(im, col, k, list(coef = vb[1L], constant = vb[2L]), lb)
    if (r$redundant) delVbds <<- c(delVbds, k)
    c(r$vlb$coef, r$vlb$constant)
  }))
  for (k in delVbds) im$vlbs[col + 1L] <- list(.hmip_hht_erase(im$vlbs[[col + 1L]], k))
  im$numVarBounds <- im$numVarBounds - length(delVbds)
  invisible(NULL)
}
