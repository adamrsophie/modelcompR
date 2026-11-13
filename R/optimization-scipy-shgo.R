# SciPy's shgo (scipy/optimize/_shgo.py, _shgo_lib/_complex.py and
# _shgo_lib/_vertex.py of SciPy 1.18.1), simplicial homology global
# optimization of Endres, Sandrock and Focke (2018).
#
# The simplicial complex keeps the neighbors of each vertex in Python sets,
# and the order in which shgo walks these sets decides which cells are
# refined first and in which order new vertices enter the vertex cache. The
# port therefore keeps the sets in a model of the hash table of CPython 3.12
# (Objects/setobject.c), with the hash of a vertex being the CPython hash of
# its coordinate tuple (Python/pyhash.c, Objects/tupleobject.c).
#
# Differences from SciPy:
# * only the simplicial sampling is ported: the 'halton', 'sobol' and
#   custom samplers build their complex with the Delaunay triangulation of
#   Qhull;
# * the local minima are sorted with the stable order(); the np.argsort of
#   the NumPy 2 wheel is stable only up to 16 elements, so ties among 17
#   or more local minima may come out in another order;
# * the messages of the disp option go to logging.info, which Python's
#   default logging drops; nothing is printed here either;
# * a warning that the stopping test repeats is signaled each time (the
#   warnings filter of Python shows a repeated text once), and the bounds in
#   the symmetry warning are printed as R prints numbers, not as numpy
#   prints arrays.

# --------------------------------------------------------------------------
# unsigned 64 bit integers as four 16 bit limbs, least significant first

.spsh_carry <- function(r) {
  carry <- 0
  for (k in 1:4) {
    t <- r[k] + carry
    r[k] <- t %% 65536
    carry <- t %/% 65536
  }
  r
}

.spsh_u64_add <- function(a, b) .spsh_carry(a + b)

.spsh_u64_mul <- function(a, b) {
  r <- numeric(4)
  for (i in 1:4) {
    for (j in 1:(5L - i)) r[i + j - 1L] <- r[i + j - 1L] + a[i] * b[j]
    # keep every partial sum below 2^53
    r <- .spsh_carry(r)
  }
  r
}

.spsh_u64_bits <- function(a) as.integer(unlist(lapply(a, function(l) (l %/% 2^(0:15)) %% 2)))

.spsh_bits_u64 <- function(b) vapply(1:4, function(k) sum(b[(16L * k - 15L):(16L * k)] * 2^(0:15)), 0)

.spsh_u64_xor <- function(a, b) as.numeric(bitwXor(as.integer(a), as.integer(b)))

.spsh_u64_rotl <- function(a, k) {
  b <- .spsh_u64_bits(a)
  .spsh_bits_u64(b[((seq_len(64L) - 1L - k) %% 64L) + 1L])
}

# the hash of a Python float (_Py_HashDouble) as 64 bits: for
# |v| = M 2^E with an integer M below 2^53 the hash is M 2^E modulo 2^61 - 1,
# a rotation of the 61 bits of M, negated for negative v; -1 becomes -2
.spsh_float_bits <- function(v) {
  if (v == 0) return(integer(64L))
  a <- abs(v)
  e <- floor(log2(a)) + 1
  m <- a / 2^e
  if (m >= 1) {
    e <- e + 1
    m <- m / 2
  } else if (m < 0.5) {
    e <- e - 1
    m <- m * 2
  }
  M <- m * 2^53
  bits <- as.integer((M %/% 2^(0:60)) %% 2)
  k <- (e - 53) %% 61
  bits <- bits[((seq_len(61L) - 1L - k) %% 61L) + 1L]
  bits <- c(bits, 0L, 0L, 0L)
  if (v < 0) {
    bits <- 1L - bits
    one <- .spsh_bits_u64(bits)
    bits <- .spsh_u64_bits(.spsh_u64_add(one, c(1, 0, 0, 0)))
    if (all(bits == 1L)) bits[1L] <- 0L
  }
  bits
}

.spsh_xx1 <- c(51847, 34283, 31153, 40503)
.spsh_xx2 <- c(60239, 10196, 44605, 49842)
.spsh_xx5 <- c(26565, 5718, 60207, 10196)
.spsh_xxlen <- c(46262, 5731, 60207, 10196)

# the hash of a tuple of floats (tuplehash, the xxHash scheme of CPython 3.8
# and later), as 64 bits
.spsh_tuple_bits <- function(x) {
  acc <- .spsh_xx5
  for (v in x) {
    lane <- .spsh_bits_u64(.spsh_float_bits(v))
    acc <- .spsh_u64_add(acc, .spsh_u64_mul(lane, .spsh_xx2))
    acc <- .spsh_u64_rotl(acc, 31L)
    acc <- .spsh_u64_mul(acc, .spsh_xx1)
  }
  acc <- .spsh_u64_add(acc, .spsh_u64_xor(c(length(x), 0, 0, 0), .spsh_xxlen))
  if (all(acc == 65535)) acc <- c(19412, 23594, 0, 0)
  .spsh_u64_bits(acc)
}

# --------------------------------------------------------------------------
# a set of vertices (keys are vertex numbers; 0 marks an empty slot and -1 a
# deleted one), with the probing, resizing and merging of setobject.c

.spsh_set_new <- function(hb) {
  s <- new.env(parent = emptyenv())
  s$mask <- 7L
  s$fill <- 0L
  s$used <- 0L
  s$keys <- integer(8L)
  s$hb <- hb
  s
}

# bits s .. s + width - 1 of the hash of key k
.spsh_hlow <- function(hb, k, s, width) {
  if (s >= 64L) return(0)
  idx <- s + seq_len(width)
  idx <- idx[idx <= 64L]
  b <- hb(k)
  sum(b[idx] * 2^(idx - s - 1L))
}

.spsh_width <- function(mask) as.integer(round(log2(mask + 1)))

.spsh_insert_clean <- function(keys, mask, k, hb) {
  w <- .spsh_width(mask)
  i <- .spsh_hlow(hb, k, 0L, w)
  sh <- 0L
  repeat {
    if (keys[i + 1L] == 0L) {
      keys[i + 1L] <- k
      return(keys)
    }
    if (i + 9 <= mask) {
      for (j in 1:9) {
        if (keys[i + j + 1L] == 0L) {
          keys[i + j + 1L] <- k
          return(keys)
        }
      }
    }
    sh <- sh + 5L
    i <- (i * 5 + 1 + .spsh_hlow(hb, k, sh, w)) %% (mask + 1)
  }
}

.spsh_resize <- function(s, minused) {
  newsize <- 8
  while (newsize <= minused) newsize <- newsize * 2
  if (newsize == 8 && s$mask == 7L && s$fill == s$used) return(invisible(s))
  old <- s$keys
  live <- old[old > 0L]
  keys <- integer(newsize)
  for (k in live) keys <- .spsh_insert_clean(keys, newsize - 1, k, s$hb)
  s$keys <- keys
  s$mask <- as.integer(newsize - 1)
  s$fill <- s$used
  invisible(s)
}

# the slot of key k, or 0 when absent
.spsh_find <- function(s, k) {
  mask <- s$mask
  w <- .spsh_width(mask)
  keys <- s$keys
  i <- .spsh_hlow(s$hb, k, 0L, w)
  sh <- 0L
  repeat {
    probes <- if (i + 9 <= mask) 9L else 0L
    for (j in 0:probes) {
      kk <- keys[i + j + 1L]
      if (kk == 0L) return(0L)
      if (kk == k) return(as.integer(i + j + 1L))
    }
    sh <- sh + 5L
    i <- (i * 5 + 1 + .spsh_hlow(s$hb, k, sh, w)) %% (mask + 1)
  }
}

.spsh_has <- function(s, k) .spsh_find(s, k) > 0L

# set_add_entry; a new key goes to the last deleted slot met on the probe
# path, else to the empty slot that ends it (checked against the tables of
# CPython 3.12 through ctypes)
.spsh_add <- function(s, k) {
  mask <- s$mask
  w <- .spsh_width(mask)
  keys <- s$keys
  i <- .spsh_hlow(s$hb, k, 0L, w)
  sh <- 0L
  freeslot <- 0L
  repeat {
    probes <- if (i + 9 <= mask) 9L else 0L
    for (j in 0:probes) {
      kk <- keys[i + j + 1L]
      if (kk == 0L) {
        if (freeslot > 0L) {
          s$keys[freeslot] <- k
          s$used <- s$used + 1L
          return(invisible(s))
        }
        s$keys[i + j + 1L] <- k
        s$fill <- s$fill + 1L
        s$used <- s$used + 1L
        if (s$fill * 5 >= mask * 3) .spsh_resize(s, if (s$used > 50000L) s$used * 2 else s$used * 4)
        return(invisible(s))
      }
      if (kk == k) return(invisible(s))
      if (kk == -1L) freeslot <- as.integer(i + j + 1L)
    }
    sh <- sh + 5L
    i <- (i * 5 + 1 + .spsh_hlow(s$hb, k, sh, w)) %% (mask + 1)
  }
}

.spsh_remove <- function(s, k) {
  j <- .spsh_find(s, k)
  if (j > 0L) {
    s$keys[j] <- -1L
    s$used <- s$used - 1L
  }
  invisible(j > 0L)
}

.spsh_items <- function(s) {
  k <- s$keys
  k[k > 0L]
}

# set_merge: so.update(other) for a set other
.spsh_update <- function(s, o) {
  if (identical(s, o) || o$used == 0L) return(invisible(s))
  if ((s$fill + o$used) * 5 >= s$mask * 3) .spsh_resize(s, (s$used + o$used) * 2)
  if (s$fill == 0L && s$mask == o$mask && o$fill == o$used) {
    s$keys <- o$keys
    s$fill <- o$fill
    s$used <- o$used
    return(invisible(s))
  }
  live <- .spsh_items(o)
  if (s$fill == 0L) {
    s$fill <- o$used
    s$used <- o$used
    keys <- s$keys
    for (k in live) keys <- .spsh_insert_clean(keys, s$mask, k, s$hb)
    s$keys <- keys
    return(invisible(s))
  }
  for (k in live) .spsh_add(s, k)
  invisible(s)
}

.spsh_copy <- function(s) {
  r <- .spsh_set_new(s$hb)
  .spsh_update(r, s)
  r
}

# --------------------------------------------------------------------------
# the vertex cache (VertexCacheField) and the simplicial complex (Complex)

# a coordinate tuple as a dictionary key; -0 and 0 are the same key, as in
# Python
.spsh_key <- function(x) paste(sprintf("%a", x + 0), collapse = " ")

.spsh_complex <- function(dim, lb, ub, field, g_cons, workers) {
  HC <- new.env(parent = emptyenv())
  HC$dim <- dim
  HC$lb <- lb
  HC$ub <- ub
  HC$field <- field
  HC$g_cons <- g_cons
  HC$parallel <- !identical(workers, 1L) && !identical(workers, 1)
  HC$key2id <- new.env(hash = TRUE, parent = emptyenv())
  HC$X <- list()
  HC$hb <- list()
  HC$nn <- list()
  HC$f <- numeric()
  HC$feasible <- logical()
  HC$chk <- logical()
  HC$mn <- logical()
  HC$N <- 0L
  HC$nfev <- 0L
  HC$hbf <- function(k) HC$hb[[k]]
  HC$gpool <- .spsh_set_new(HC$hbf)
  HC$fpool <- .spsh_set_new(HC$hbf)
  HC$se <- new.env(hash = TRUE, parent = emptyenv())
  HC
}

.spsh_in_cache <- function(HC, x) !is.null(HC$key2id[[.spsh_key(x)]])

# V[x]: the vertex at x, created (and queued for evaluation) when new
.spsh_V <- function(HC, x) {
  k <- .spsh_key(x)
  id <- HC$key2id[[k]]
  if (!is.null(id)) return(id)
  id <- HC$N + 1L
  HC$N <- id
  assign(k, id, envir = HC$key2id)
  HC$X[[id]] <- x
  HC$hb[[id]] <- .spsh_tuple_bits(x)
  HC$nn[[id]] <- .spsh_set_new(HC$hbf)
  HC$f[id] <- NA_real_
  HC$feasible[id] <- NA
  HC$chk[id] <- TRUE
  HC$mn[id] <- FALSE
  .spsh_add(HC$gpool, id)
  .spsh_add(HC$fpool, id)
  id
}

.spsh_connect <- function(HC, a, b) {
  if (a != b && !.spsh_has(HC$nn[[a]], b)) {
    .spsh_add(HC$nn[[a]], b)
    .spsh_add(HC$nn[[b]], a)
    HC$chk[c(a, b)] <- TRUE
  }
  invisible(NULL)
}

# the Python errors that a refinement generator may raise, and that
# Complex.refine catches
.spsh_raise <- function(cls) stop(structure(class = c(cls, "spsh_pyerror", "error", "condition"),
                                            list(message = cls, call = NULL)))

.spsh_disconnect <- function(HC, a, b) {
  if (.spsh_has(HC$nn[[a]], b)) {
    .spsh_remove(HC$nn[[a]], b)
    if (!.spsh_remove(HC$nn[[b]], a)) .spsh_raise("KeyError")
    HC$chk[c(a, b)] <- TRUE
  }
  invisible(NULL)
}

# Complex.split_edge, memoized on the ordered pair of coordinate tuples
# (functools.cache); a call that raises is not remembered
.spsh_split_edge <- function(HC, x1, x2) {
  key <- paste(.spsh_key(x1), .spsh_key(x2), sep = "|")
  id <- HC$se[[key]]
  if (!is.null(id)) return(id)
  v1 <- .spsh_V(HC, x1)
  v2 <- .spsh_V(HC, x2)
  .spsh_disconnect(HC, v1, v2)
  vct <- (HC$X[[v2]] - HC$X[[v1]]) / 2.0 + HC$X[[v1]]
  vc <- .spsh_V(HC, vct)
  .spsh_connect(HC, vc, v1)
  .spsh_connect(HC, vc, v2)
  assign(key, vc, envir = HC$se)
  vc
}

.spsh_x <- function(HC, v) {
  if (is.null(v)) .spsh_raise("AttributeError")
  HC$X[[v]]
}

# the pairs of itertools.combinations(comb, 2) split in turn
.spsh_split_pairs <- function(HC, comb) {
  m <- length(comb)
  for (p in seq_len(m - 1L)) for (q in (p + 1L):m)
    .spsh_split_edge(HC, .spsh_x(HC, comb[[p]]), .spsh_x(HC, comb[[q]]))
}

# Complex.cyclic_product: the triangulation of the hypercube with its
# centroid; ys() stands at each yield of the generator. (Its IndexError
# branch needs a bounds row with one entry, which SHGO never passes.)
.spsh_cyclic_product <- function(HC, origin, supremum, ys) {
  dim <- length(origin)
  vot <- origin
  vut <- supremum
  .spsh_V(HC, vot)
  ys()
  .spsh_connect(HC, .spsh_V(HC, vut), .spsh_V(HC, vot))
  ys()
  C0x <- list(.spsh_V(HC, vot))
  a_vo <- origin
  a_vo[1L] <- vut[1L]
  a_vo <- .spsh_V(HC, a_vo)
  .spsh_connect(HC, .spsh_V(HC, vot), a_vo)
  ys()
  C1x <- list(a_vo)
  ab_C <- list()
  for (i in seq_len(dim - 1L) - 1L) {
    ip <- i + 2L
    C0x[[ip]] <- integer()
    C1x[[ip]] <- integer()
    cC0x <- C0x[seq_len(i + 1L)]
    cC1x <- C1x[seq_len(i + 1L)]
    for (j in seq_along(cC0x)) {
      VL <- cC0x[[j]]
      VU <- cC1x[[j]]
      for (k in seq_len(min(length(VL), length(VU)))) {
        vl <- VL[k]
        vu <- VU[k]
        a_vl <- HC$X[[vl]]
        a_vu <- HC$X[[vu]]
        a_vl[ip] <- vut[ip]
        a_vu[ip] <- vut[ip]
        a_vl <- .spsh_V(HC, a_vl)
        .spsh_connect(HC, vl, a_vl)
        ys()
        a_vu <- .spsh_V(HC, a_vu)
        .spsh_connect(HC, vu, a_vu)
        .spsh_connect(HC, a_vl, a_vu)
        .spsh_connect(HC, vl, a_vu)
        ab_C[[length(ab_C) + 1L]] <- c(vl, a_vu)
        C0x[[ip]] <- c(C0x[[ip]], vl, vu)
        C1x[[ip]] <- c(C1x[[ip]], a_vl, a_vu)
        C0x[[j]] <- c(C0x[[j]], a_vl)
        C1x[[j]] <- c(C1x[[j]], a_vu)
        ys()
      }
    }
    for (vp in ab_C) {
      b_v <- HC$X[[vp[1L]]]
      ab_v <- HC$X[[vp[2L]]]
      b_v[ip] <- vut[ip]
      ab_v[ip] <- vut[ip]
      b_v <- .spsh_V(HC, b_v)
      ab_v <- .spsh_V(HC, ab_v)
      .spsh_connect(HC, vp[1L], ab_v)
      .spsh_connect(HC, b_v, ab_v)
      ab_C[[length(ab_C) + 1L]] <- c(vp[1L], ab_v)
      ab_C[[length(ab_C) + 1L]] <- c(b_v, ab_v)
    }
  }
  vo <- .spsh_V(HC, vot)
  vs <- .spsh_V(HC, vut)
  .spsh_disconnect(HC, vo, vs)
  vc <- .spsh_split_edge(HC, vot, vut)
  for (v in .spsh_items(HC$nn[[vo]])) .spsh_connect(HC, v, vc)
  ys()
  invisible(NULL)
}

.spsh_tv <- function(o, s) list(o = o, s = s, key = paste(.spsh_key(o), .spsh_key(s), sep = "|"))

# Complex.vpool: the neighbors of origin and supremum inside their box
.spsh_vpool <- function(HC, origin, supremum) {
  vo <- .spsh_V(HC, origin)
  vs <- .spsh_V(HC, supremum)
  bl <- pmin(origin, supremum)
  bu <- pmax(origin, supremum)
  pool <- .spsh_set_new(HC$hbf)
  .spsh_update(pool, HC$nn[[vo]])
  .spsh_update(pool, HC$nn[[vs]])
  for (vn in .spsh_items(.spsh_copy(pool))) {
    x <- HC$X[[vn]]
    if (any(!(bl <= x & x <= bu))) .spsh_remove(pool, vn)
  }
  pool
}

# Complex.refine_local_space: the refinement of the cell from origin to
# supremum; ys() stands at each yield. Python keeps vl, vu, a_vl, a_vu and
# c_vc from one loop to the next, and so does the port.
.spsh_refine_local_space <- function(HC, S, origin, supremum, ys) {
  dim <- length(origin)
  origin_c <- origin
  supremum_c <- supremum
  vl <- NULL
  vu <- NULL
  a_vu <- NULL
  s_origin <- origin
  s_supremum <- supremum
  for (i in seq_len(dim)) {
    if (origin[i] > supremum[i]) {
      s_origin[i] <- supremum[i]
      s_supremum[i] <- origin[i]
    }
  }
  vot <- s_origin
  vut <- s_supremum
  X <- function(v) .spsh_x(HC, v)
  V <- function(x) .spsh_V(HC, x)
  con <- function(a, b) .spsh_connect(HC, a, b)
  se <- function(a, b) .spsh_split_edge(HC, a, b)
  vo <- V(vot)
  vs <- V(vut)
  vco <- se(X(vo), X(vs))
  sup_set <- .spsh_copy(HC$nn[[vco]])
  a_vl <- vot
  a_vl[1L] <- vut[1L]
  if (!.spsh_in_cache(HC, a_vl)) {
    vo <- V(vot)
    vs <- V(vut)
    vco <- se(X(vo), X(vs))
    sup_set <- .spsh_copy(HC$nn[[vco]])
    a_vl <- vot
    a_vl[1L] <- vut[1L]
    a_vl <- V(a_vl)
  } else {
    a_vl <- V(a_vl)
  }
  c_v <- se(X(vo), X(a_vl))
  con(c_v, vco)
  ys()
  Cox <- list(vo)
  Ccx <- list(c_v)
  Cux <- list(a_vl)
  ab_C <- list()
  s_ab_C <- list()
  for (i in seq_len(dim - 1L) - 1L) {
    ip <- i + 2L
    Cox[[ip]] <- integer()
    Ccx[[ip]] <- integer()
    Cux[[ip]] <- integer()
    t_a_vl <- vot
    t_a_vl[ip] <- vut[ip]
    cCox <- Cox[seq_len(i + 1L)]
    cCcx <- Ccx[seq_len(i + 1L)]
    cCux <- Cux[seq_len(i + 1L)]
    ab_Cc <- ab_C
    s_ab_Cc <- s_ab_C
    # the IndexError of Python: the corner of the next direction is missing
    index_error <- !.spsh_in_cache(HC, t_a_vl) || !.spsh_in_cache(HC, vut)
    if (!index_error) {
      for (vectors in s_ab_Cc) {
        bc_vc <- X(vectors[1L])
        b_vl <- X(vectors[2L])
        b_vu <- X(vectors[3L])
        ba_vu <- X(vectors[4L])
        bc_vc[ip] <- vut[ip]
        b_vl[ip] <- vut[ip]
        b_vu[ip] <- vut[ip]
        ba_vu[ip] <- vut[ip]
        bc_vc <- V(bc_vc)
        con(bc_vc, vco)
        ys()
        d_bc_vc <- se(X(vectors[1L]), X(bc_vc))
        con(d_bc_vc, bc_vc)
        con(d_bc_vc, vectors[2L])
        con(d_bc_vc, vectors[3L])
        con(d_bc_vc, vectors[4L])
        ys()
        b_vl <- V(b_vl)
        con(bc_vc, b_vl)
        con(d_bc_vc, b_vl)
        ys()
        b_vu <- V(b_vu)
        con(bc_vc, b_vu)
        con(d_bc_vc, b_vu)
        b_vl_c <- se(X(b_vu), X(b_vl))
        con(bc_vc, b_vl_c)
        ys()
        ba_vu <- V(ba_vu)
        con(bc_vc, ba_vu)
        con(d_bc_vc, ba_vu)
        os_v <- se(X(vectors[2L]), X(ba_vu))
        ss_v <- se(X(b_vl), X(ba_vu))
        b_vu_c <- se(X(b_vu), X(ba_vu))
        con(bc_vc, b_vu_c)
        ys()
        ys()
        ys()
        d_bc_vc <- se(X(vectors[1L]), X(bc_vc))
        con(d_bc_vc, vco)
        ys()
        d_b_vl <- se(X(vectors[2L]), X(b_vl))
        con(d_bc_vc, vco)
        con(d_bc_vc, d_b_vl)
        ys()
        d_b_vu <- se(X(vectors[3L]), X(b_vu))
        con(d_bc_vc, vco)
        con(d_bc_vc, d_b_vu)
        ys()
        d_ba_vu <- se(X(vectors[4L]), X(ba_vu))
        con(d_bc_vc, vco)
        con(d_bc_vc, d_ba_vu)
        ys()
        .spsh_split_pairs(HC, list(vl, vu, a_vu, b_vl, b_vu, ba_vu))
        ab_C[[length(ab_C) + 1L]] <- c(d_bc_vc, vectors[2L], b_vl, a_vu, ba_vu)
        ab_C[[length(ab_C) + 1L]] <- c(d_bc_vc, vl, b_vl, a_vu, ba_vu)
      }
      for (vectors in ab_Cc) {
        bc_vc <- X(vectors[1L])
        b_vl <- X(vectors[2L])
        b_vu <- X(vectors[3L])
        ba_vl <- X(vectors[4L])
        ba_vu <- X(vectors[5L])
        bc_vc[ip] <- vut[ip]
        b_vl[ip] <- vut[ip]
        b_vu[ip] <- vut[ip]
        ba_vl[ip] <- vut[ip]
        ba_vu[ip] <- vut[ip]
        bc_vc <- V(bc_vc)
        con(bc_vc, vco)
        ys()
        d_bc_vc <- se(X(vectors[1L]), X(bc_vc))
        con(d_bc_vc, bc_vc)
        con(d_bc_vc, vectors[2L])
        con(d_bc_vc, vectors[3L])
        con(d_bc_vc, vectors[4L])
        con(d_bc_vc, vectors[5L])
        ys()
        b_vl <- V(b_vl)
        con(bc_vc, b_vl)
        con(d_bc_vc, b_vl)
        ys()
        b_vu <- V(b_vu)
        con(bc_vc, b_vu)
        con(d_bc_vc, b_vu)
        ys()
        ba_vl <- V(ba_vl)
        con(bc_vc, ba_vl)
        con(d_bc_vc, ba_vl)
        se(X(b_vu), X(ba_vl))
        ys()
        ba_vu <- V(ba_vu)
        con(bc_vc, ba_vu)
        con(d_bc_vc, ba_vu)
        os_v <- se(X(vectors[2L]), X(ba_vu))
        ss_v <- se(X(b_vl), X(ba_vu))
        ys()
        ys()
        ys()
        d_bc_vc <- se(X(vectors[1L]), X(bc_vc))
        con(d_bc_vc, vco)
        ys()
        d_b_vl <- se(X(vectors[2L]), X(b_vl))
        con(d_bc_vc, vco)
        con(d_bc_vc, d_b_vl)
        ys()
        d_b_vu <- se(X(vectors[3L]), X(b_vu))
        con(d_bc_vc, vco)
        con(d_bc_vc, d_b_vu)
        ys()
        d_ba_vl <- se(X(vectors[4L]), X(ba_vl))
        con(d_bc_vc, vco)
        con(d_bc_vc, d_ba_vl)
        ys()
        d_ba_vu <- se(X(vectors[5L]), X(ba_vu))
        con(d_bc_vc, vco)
        con(d_bc_vc, d_ba_vu)
        ys()
        c_vc <- vectors[1L]
        vl <- vectors[2L]
        vu <- vectors[3L]
        a_vl <- vectors[4L]
        a_vu <- vectors[5L]
        .spsh_split_pairs(HC, list(vl, vu, a_vl, a_vu, b_vl, b_vu, ba_vl, ba_vu))
        ab_C[[length(ab_C) + 1L]] <- c(bc_vc, b_vl, b_vu, ba_vl, ba_vu)
        ab_C[[length(ab_C) + 1L]] <- c(d_bc_vc, d_b_vl, d_b_vu, d_ba_vl, d_ba_vu)
        ab_C[[length(ab_C) + 1L]] <- c(d_bc_vc, vectors[2L], b_vl, a_vu, ba_vu)
        ab_C[[length(ab_C) + 1L]] <- c(d_bc_vc, vu, b_vu, a_vl, ba_vl)
      }
      for (j in seq_along(cCox)) {
        VL <- cCox[[j]]
        VC <- cCcx[[j]]
        VU <- cCux[[j]]
        for (k in seq_len(min(length(VL), length(VC), length(VU)))) {
          vl <- VL[k]
          vc <- VC[k]
          vu <- VU[k]
          a_vl <- X(vl)
          a_vu <- X(vu)
          a_vl[ip] <- vut[ip]
          a_vu[ip] <- vut[ip]
          a_vl <- V(a_vl)
          a_vu <- V(a_vu)
          c_vc <- se(X(vl), X(a_vu))
          se(X(vl), X(vu))
          con(c_vc, vco)
          con(c_vc, vc)
          con(c_vc, vl)
          con(c_vc, vu)
          con(c_vc, a_vl)
          con(c_vc, a_vu)
          ys()
          c_vl <- se(X(vl), X(a_vl))
          con(c_vl, vco)
          con(c_vc, c_vl)
          ys()
          c_vu <- se(X(vu), X(a_vu))
          con(c_vu, vco)
          con(c_vc, c_vu)
          ys()
          a_vc <- se(X(a_vl), X(a_vu))
          con(a_vc, vco)
          con(a_vc, c_vc)
          ab_C[[length(ab_C) + 1L]] <- c(c_vc, vl, vu, a_vl, a_vu)
          Cox[[ip]] <- c(Cox[[ip]], vl, vc, vu)
          Ccx[[ip]] <- c(Ccx[[ip]], c_vl, c_vc, c_vu)
          Cux[[ip]] <- c(Cux[[ip]], a_vl, a_vc, a_vu)
          Cox[[j]] <- c(Cox[[j]], c_vl, a_vl)
          Ccx[[j]] <- c(Ccx[[j]], c_vc, a_vc)
          Cux[[j]] <- c(Cux[[j]], c_vu, a_vu)
          ys()
        }
      }
    } else {
      for (vectors in ab_Cc) {
        ba_vl <- X(vectors[4L])
        ba_vu <- X(vectors[5L])
        ba_vl[ip] <- vut[ip]
        ba_vu[ip] <- vut[ip]
        ba_vu <- V(ba_vu)
        ys()
        d_bc_vc <- se(X(vectors[2L]), X(ba_vu))
        ys()
        con(d_bc_vc, vectors[2L])
        con(d_bc_vc, vectors[3L])
        con(d_bc_vc, vectors[4L])
        con(d_bc_vc, vectors[5L])
        ys()
        ba_vl <- V(ba_vl)
        ys()
        d_ba_vl <- se(X(vectors[4L]), X(ba_vl))
        d_ba_vu <- se(X(vectors[5L]), X(ba_vu))
        d_ba_vc <- se(X(d_ba_vl), X(d_ba_vu))
        ys()
        ys()
        ys()
        c_vc <- vectors[1L]
        vl <- vectors[2L]
        vu <- vectors[3L]
        a_vl <- vectors[4L]
        a_vu <- vectors[5L]
        .spsh_split_pairs(HC, list(vl, vu, a_vl, a_vu, ba_vl, ba_vu))
      }
      VL <- Cox[[i + 1L]]
      VC <- Ccx[[i + 1L]]
      VU <- Cux[[i + 1L]]
      for (k in seq_len(min(length(VL), length(VC), length(VU)))) {
        vl <- VL[k]
        vc <- VC[k]
        vu <- VU[k]
        a_vu <- X(vu)
        a_vu[ip] <- vut[ip]
        a_vu <- V(a_vu)
        X(a_vl)
        ys()
        c_vc <- se(X(vl), X(a_vu))
        se(X(vl), X(vu))
        con(c_vc, vco)
        con(c_vc, vc)
        con(c_vc, vl)
        con(c_vc, vu)
        con(c_vc, a_vu)
        ys()
        c_vu <- se(X(vu), X(a_vu))
        con(c_vu, vco)
        con(c_vc, c_vu)
        ys()
        Cox[[ip]] <- c(Cox[[ip]], vu)
        Ccx[[ip]] <- c(Ccx[[ip]], c_vu)
        Cux[[ip]] <- c(Cux[[ip]], a_vu)
        s_ab_C[[length(s_ab_C) + 1L]] <- c(c_vc, vl, vu, a_vu)
        ys()
      }
    }
  }
  tv_key <- .spsh_tv(origin_c, supremum_c)$key
  hit <- which(vapply(S$tv, function(t) t$key == tv_key, TRUE))
  if (length(hit)) S$tv <- S$tv[-hit[1L]]
  sups <- .spsh_items(sup_set)
  for (vs in sups) S$tv[[length(S$tv) + 1L]] <- .spsh_tv(X(vco), X(vs))
  c_nn_lists <- lapply(sups, function(vs) .spsh_vpool(HC, X(vco), X(vs)))
  for (t in seq_along(sups)) {
    vcn <- se(X(vco), X(sups[t]))
    for (vnn in .spsh_items(c_nn_lists[[t]])) con(vcn, vnn)
    ys()
  }
  ys()
  invisible(NULL)
}

# --------------------------------------------------------------------------
# evaluation of the pools (VertexCacheField.process_pools)

.spsh_feasible <- function(HC, x) {
  for (g in HC$g_cons) if (any(as.numeric(g(x)) < 0.0)) return(FALSE)
  TRUE
}

.spsh_process_pools <- function(HC) {
  if (!is.null(HC$g_cons)) {
    if (!HC$parallel) {
      for (v in .spsh_items(HC$gpool)) {
        HC$feasible[v] <- .spsh_feasible(HC, HC$X[[v]])
        if (!HC$feasible[v]) HC$f[v] <- Inf
      }
      HC$gpool <- .spsh_set_new(HC$hbf)
    } else {
      # pproc_gpool keeps its pool, so each pass checks every vertex again
      for (v in .spsh_items(HC$gpool)) HC$feasible[v] <- .spsh_feasible(HC, HC$X[[v]])
    }
  }
  for (v in .spsh_items(HC$fpool)) {
    if (!is.null(HC$g_cons) && !HC$feasible[v]) {
      if (HC$parallel) HC$f[v] <- Inf
      next
    }
    fv <- if (HC$parallel) tryCatch(as.numeric(HC$field(HC$X[[v]])), error = function(e) Inf)
          else as.numeric(HC$field(HC$X[[v]]))
    HC$nfev <- HC$nfev + 1L
    HC$f[v] <- if (is.na(fv)) Inf else fv
  }
  HC$fpool <- .spsh_set_new(HC$hbf)
  for (v in seq_len(HC$N)) .spsh_minimiser(HC, v)
  invisible(NULL)
}

# VertexScalarField.minimiser, with its cached answer
.spsh_minimiser <- function(HC, v) {
  if (HC$chk[v]) {
    fv <- HC$f[v]
    HC$mn[v] <- all(vapply(.spsh_items(HC$nn[[v]]), function(u) fv < HC$f[u], TRUE))
    HC$chk[v] <- FALSE
  }
  HC$mn[v]
}

# --------------------------------------------------------------------------
# SHGO

.spsh_lmc_get <- function(S, x) {
  k <- .spsh_key(x)
  if (is.null(S$lmc[[k]])) {
    S$lmc[[k]] <- list(x = x, x_l = NULL, lres = NULL, f_min = NULL)
    S$lmc_keys <- c(S$lmc_keys, k)
  }
  S$lmc[[k]]
}

.spsh_find_lowest_vertex <- function(S) {
  HC <- S$HC
  fl <- Inf
  xl <- NULL
  for (v in seq_len(HC$N)) {
    if (isTRUE(HC$f[v] < fl)) {
      fl <- HC$f[v]
      xl <- HC$X[[v]]
    }
  }
  for (k in S$lmc_keys) {
    m <- S$lmc[[k]]
    if (isTRUE(m$f_min < fl)) {
      fl <- m$f_min
      xl <- m$x_l
    }
  }
  if (fl == Inf) {
    fl <- NULL
    xl <- NULL
  }
  S$f_lowest <- fl
  S$x_lowest <- xl
}

.spsh_stopping_criteria <- function(S) {
  if (!is.null(S$iters) && S$iters_done >= S$iters) S$stop_global <- TRUE
  if (!is.null(S$maxiter) && S$iters_done >= S$maxiter) S$stop_global <- TRUE
  if (!is.null(S$maxfev) && S$fn >= S$maxfev) S$stop_global <- TRUE
  if (!is.null(S$maxev) && S$n_sampled >= S$maxev) S$stop_global <- TRUE
  if (!is.null(S$maxtime) && as.numeric(Sys.time()) - S$init >= S$maxtime) S$stop_global <- TRUE
  if (!is.null(S$f_min_true)) {
    .spsh_find_lowest_vertex(S)
    if (!is.null(S$f_lowest)) {
      if (S$f_min_true == 0.0) {
        if (S$f_lowest <= S$f_tol) S$stop_global <- TRUE
      } else {
        pe <- (S$f_lowest - S$f_min_true) / abs(S$f_min_true)
        if (S$f_lowest <= S$f_min_true) {
          S$stop_global <- TRUE
          if (abs(pe) >= 2 * S$f_tol)
            warning("A much lower value than expected f* = ", .spsq_repr(S$f_min_true), " was found f_lowest = ",
                    .spsq_repr(S$f_lowest), call. = FALSE)
        }
        if (pe <= S$f_tol) S$stop_global <- TRUE
      }
    }
  }
  if (!is.null(S$minhgrd) && S$lmc_size > 0L) {
    hgrd <- S$lmc_size - S$hgr
    S$hgr <- S$lmc_size
    if (hgrd <= S$minhgrd) S$stop_global <- TRUE
  }
  S$stop_global
}

# construct_lcb_simplicial: the box spanned by the neighbors of v_min
.spsh_lcb_simplicial <- function(S, v) {
  HC <- S$HC
  lo <- HC$lb
  hi <- HC$ub
  xv <- HC$X[[v]]
  for (vn in .spsh_items(HC$nn[[v]])) {
    x <- HC$X[[vn]]
    for (i in seq_along(x)) {
      if (x[i] < xv[i] && x[i] > lo[i]) lo[i] <- x[i]
      if (x[i] > xv[i] && x[i] < hi[i]) hi[i] <- x[i]
    }
  }
  list(lb = lo, ub = hi)
}

# SHGO.minimize: the local minimization from the vertex v
.spsh_minimize <- function(S, v) {
  HC <- S$HC
  x_min <- HC$X[[v]]
  m <- .spsh_lmc_get(S, x_min)
  if (!is.null(m$lres)) return(m$lres)
  mk <- S$mk
  g_bounds <- .spsh_lcb_simplicial(S, v)
  opts <- mk$options
  if (!is.null(mk$callback)) opts$callback <- mk$callback
  if (!is.null(mk$hessp)) opts$hessp <- mk$hessp
  lres <- .sp_minimize(S$func, x_min, method = mk$method, jac = mk$jac, hess = mk$hess,
                       bounds = if (S$takes_bounds) g_bounds, constraints = mk$constraints %||% list(),
                       options = opts)
  S$nlfev <- S$nlfev + lres$nfev
  if (!is.null(lres$njev)) S$nljev <- S$nljev + lres$njev
  if (!is.null(lres$nhev)) S$nlhev <- S$nlhev + lres$nhev
  k <- .spsh_key(x_min)
  S$lmc[[k]]$x_l <- lres$x
  S$lmc[[k]]$lres <- lres
  S$lmc[[k]]$f_min <- lres$fun
  S$lmc_size <- S$lmc_size + 1L
  S$xl_maps[[length(S$xl_maps) + 1L]] <- lres$x
  S$f_maps <- c(S$f_maps, lres$fun)
  lres
}

# minimizers(): the vertices below all their neighbors, in cache order,
# skipping the points already found by a local minimization
.spsh_minimizers <- function(S) {
  HC <- S$HC
  pool <- integer()
  for (v in seq_len(HC$N)) {
    x <- HC$X[[v]]
    if (any(vapply(S$xl_maps, function(xl) all(x == xl), TRUE))) next
    if (.spsh_minimiser(HC, v) && !(v %in% pool)) pool <- c(pool, v)
  }
  pool
}

# minimise_pool: the first minimizer of the cache first (X_min is never
# sorted), then always the one farthest from the last local solution
.spsh_minimise_pool <- function(S, X_min, force_iter) {
  HC <- S$HC
  lres <- .spsh_minimize(S, X_min[1L])
  X_min <- X_min[-1L]
  repeat {
    .spsh_stopping_criteria(S)
    if (force_iter) {
      force_iter <- force_iter - 1
      if (force_iter == 0) break
    }
    if (!length(X_min)) break
    d <- vapply(X_min, function(v) {
      s <- 0.0
      x <- HC$X[[v]]
      for (j in seq_along(x)) {
        t <- abs(lres$x[j] - x[j])
        s <- s + t * t
      }
      sqrt(s)
    }, 0.0)
    far <- order(d)[length(d)]
    lres <- .spsh_minimize(S, X_min[far])
    X_min <- X_min[-far]
  }
  invisible(NULL)
}

.spsh_sort_result <- function(S) {
  o <- order(S$f_maps)
  S$res$xl <- do.call(rbind, S$xl_maps[o])
  S$res$funl <- S$f_maps[o]
  S$res$x <- S$xl_maps[[o[1L]]]
  S$res$fun <- S$f_maps[o[1L]]
  S$res$nfev <- S$fn + S$nlfev
}

.spsh_find_minima <- function(S) {
  X_min <- .spsh_minimizers(S)
  if (length(X_min)) {
    .spsh_minimise_pool(S, X_min, S$local_iter)
    .spsh_sort_result(S)
    S$f_lowest <- S$res$fun
    S$x_lowest <- S$res$x
  } else {
    .spsh_find_lowest_vertex(S)
  }
}

# the part of SHGO.iterate after the refinement of the complex
.spsh_end_iteration <- function(S) {
  HC <- S$HC
  if (is.null(S$n)) S$n_sampled <- HC$N else S$n_sampled <- S$n_sampled + S$n
  if (length(S$xl_maps)) {
    # v.star() puts each start vertex of a local minimization into its own
    # neighbor set
    for (k in S$lmc_keys) {
      v <- .spsh_V(HC, S$lmc[[k]]$x)
      .spsh_add(HC$nn[[v]], v)
    }
  }
  .spsh_process_pools(HC)
  S$fn <- HC$nfev
  if (S$minimize_every_iter) .spsh_find_minima(S)
  S$iters_done <- S$iters_done + 1L
  .spsh_stopping_criteria(S)
  if (S$stop_global) stop(structure(class = c("spsh_done", "condition"), list(message = "", call = NULL)))
}

# the iterations of SHGO.iterate_all with the simplicial sampling. The
# generators of Complex (cyclic_product, refine_local_space) run as plain
# functions; at each yield the vertex count is compared with the target of
# Complex.refine, and once it is met the rest of the iteration runs there,
# inside the suspended generator, as it does when Python resumes it later.
.spsh_iterate_simplicial <- function(S) {
  HC <- S$HC
  origin <- HC$lb
  supremum <- HC$ub
  S$tv <- NULL
  S$cp_done <- FALSE
  if (S$symmetry_warn) .spsh_symmetry_warning(HC)
  if (is.null(S$n)) {
    # Complex.refine_all, a whole generator per iteration
    repeat {
      if (is.null(S$tv)) {
        .spsh_cyclic_product(HC, origin, supremum, function() NULL)
        S$tv <- list(.spsh_tv(origin, supremum))
      } else {
        for (vp in S$tv) .spsh_refine_local_space(HC, S, vp$o, vp$s, function() NULL)
      }
      .spsh_end_iteration(S)
    }
  }
  checkpoint <- function() {
    while (HC$N >= S$nt) {
      .spsh_end_iteration(S)
      S$nt <- HC$N + S$n
      if (is.null(S$tv) && S$symmetry_warn) .spsh_symmetry_warning(HC)
    }
  }
  S$nt <- HC$N + S$n
  checkpoint()
  .spsh_cyclic_product(HC, origin, supremum, checkpoint)
  S$tv <- list(.spsh_tv(origin, supremum))
  checkpoint()
  repeat {
    vp <- S$tv[[1L]]
    yielded <- FALSE
    ys <- function() {
      yielded <<- TRUE
      checkpoint()
    }
    failed <- tryCatch({
      .spsh_refine_local_space(HC, S, vp$o, vp$s, ys)
      FALSE
    }, spsh_pyerror = function(e) TRUE)
    if (failed && !yielded) {
      # an error in the first step of a new generator reaches the outer
      # handler of Complex.refine, which goes back to triangulate()
      if (S$symmetry_warn) .spsh_symmetry_warning(HC)
      if (HC$N < S$nt) S$tv[[length(S$tv) + 1L]] <- .spsh_tv(origin, supremum)
      checkpoint()
    }
  }
}

# Complex.triangulate with the symmetry option logs that the bounds of the
# symmetric variables do not match (it compares numpy rows by identity)
.spsh_symmetry_warning <- function(HC) {
  b <- paste0("[", .spsq_repr(HC$lb[1L]), ", ", .spsq_repr(HC$ub[1L]), "]")
  for (i in seq_len(HC$dim - 1L))
    warning("Variable ", i, " was specified as symmetric to variable 0, however, the bounds ", i, " = ", b,
            " and 0 = ", b, " do not match, the mismatch was ignored in the initial triangulation.", call. = FALSE)
}

# old style constraints of shgo (standardize_constraints(..., 'old')); a
# type is kept for the feasibility test when it is a substring of 'ineq',
# as `cons['type'] in 'ineq'` reads, so equality constraints count too
.spsh_g_cons <- function(constraints, dim) {
  if (is.null(constraints)) return(NULL)
  if (!is.null(constraints$type) || !is.null(constraints$fun)) constraints <- list(constraints)
  g <- list()
  for (con in constraints) {
    type <- con$type %||% "nonlinear"
    if (type == "nonlinear") {
      old <- .spsh_new_to_old(con, numeric(dim))
      for (o in old) g[[length(g) + 1L]] <- o$fun
      next
    }
    if (grepl(type, "ineq", fixed = TRUE)) {
      local({
        f <- con$fun
        g[[length(g) + 1L]] <<- function(x) f(x)
      })
    }
  }
  g
}

# standardize_constraints(constraints, x0, meth) of minimize() for COBYLA,
# COBYQA and trust-constr (old_constraint_to_new): an old style
# list(type, fun, jac, args) becomes the NonlinearConstraint 0 <= fun(x) <= 0
# (eq) or 0 <= fun(x) <= Inf (ineq)
.sp_new_constraints <- function(constraints) {
  if (is.null(constraints)) return(list())
  if (!is.null(constraints$type) || !is.null(constraints$fun)) constraints <- list(constraints)
  lapply(seq_along(constraints), function(ic) {
    con <- constraints[[ic]]
    if (is.null(con$type) && !is.null(con$lb)) return(con)
    if (is.null(con$type)) stop("KeyError: 'Constraint ", ic - 1L, " has no type defined.'", call. = FALSE)
    ctype <- tolower(con$type)
    if (ctype %in% c("nonlinear", "linear")) return(con)
    if (!(ctype %in% c("eq", "ineq"))) stop("Unknown constraint type '", con$type, "'.", call. = FALSE)
    if (is.null(con$fun)) stop("Constraint ", ic - 1L, " has no function defined.", call. = FALSE)
    args <- con$args
    fun <- if (is.null(args)) con$fun else function(x) do.call(con$fun, c(list(x), args))
    jac <- "2-point"
    if (!is.null(con$jac)) jac <- if (is.null(args)) con$jac else function(x) do.call(con$jac, c(list(x), args))
    .sp_nonlinear_constraint(fun, 0, if (ctype == "eq") 0 else Inf, jac)
  })
}

# new_constraint_to_old for a NonlinearConstraint: its equality part, then
# its inequality part (both are feasibility tests here)
.spsh_new_to_old <- function(con, x0) {
  fun <- con$fun
  m <- length(as.numeric(fun(x0)))
  lb <- rep_len(as.numeric(con$lb), m)
  ub <- rep_len(as.numeric(con$ub), m)
  i_eq <- lb == ub
  i_bb <- xor(lb != -Inf, i_eq)
  i_ba <- xor(ub != Inf, i_eq)
  out <- list()
  if (any(i_eq)) out[[1L]] <- list(type = "eq", fun = function(x) {
    y <- as.numeric(fun(x))
    y[i_eq] - lb[i_eq]
  })
  if (any(i_bb) || any(i_ba)) out[[length(out) + 1L]] <- list(type = "ineq", fun = function(x) {
    y <- as.numeric(fun(x))
    c(y[i_bb] - lb[i_bb], -(y[i_ba] - ub[i_ba]))
  })
  out
}

# MemoizeJac: a function returning list(value, gradient)
.spsh_memoize_jac <- function(fun) {
  last <- NULL
  val <- NULL
  grad <- NULL
  upd <- function(x) {
    if (is.null(last) || length(last) != length(x) || !all(x == last) || is.null(val) || is.null(grad)) {
      r <- fun(x)
      last <<- x
      val <<- as.numeric(r[[1L]])
      grad <<- as.numeric(r[[2L]])
    }
  }
  list(f = function(x) {
    upd(x)
    val
  }, jac = function(x) {
    upd(x)
    grad
  })
}

.spsh_solver_args <- list(
  "nelder-mead" = character(), powell = character(), cg = "jac", bfgs = "jac",
  "newton-cg" = c("jac", "hess", "hessp"), "l-bfgs-b" = c("jac", "bounds"), tnc = c("jac", "bounds"),
  cobyla = c("constraints", "catol"), cobyqa = c("bounds", "constraints", "feasibility_tol"),
  slsqp = c("jac", "bounds", "constraints"), dogleg = c("jac", "hess"), "trust-ncg" = c("jac", "hess", "hessp"),
  "trust-krylov" = c("jac", "hess", "hessp"), "trust-exact" = c("jac", "hess"),
  "trust-constr" = c("jac", "hess", "hessp", "constraints"))

# scipy.optimize.shgo with the simplicial sampling; the bounds are lb, ub
.sp_shgo <- function(func, lb, ub, constraints = NULL, n = 100L, iters = 1L, callback = NULL,
                     minimizer_kwargs = NULL, options = NULL, sampling_method = "simplicial", workers = 1L) {
  if (is.character(sampling_method) && !(sampling_method %in% c("halton", "sobol", "simplicial")))
    stop("Unknown sampling_method specified. Valid methods: halton, sobol, simplicial", call. = FALSE)
  if (!identical(sampling_method, "simplicial"))
    stop("shgo: the sampling_method '", if (is.character(sampling_method)) sampling_method else "custom",
         "' builds its complex with the Delaunay triangulation of Qhull, which modelcompR does not port; ",
         "use sampling_method = \"simplicial\"", call. = FALSE)
  S <- new.env(parent = emptyenv())
  mk_in <- minimizer_kwargs
  if (isTRUE(options$jac)) {
    if (is.null(mk_in)) mk_in <- list()
    mk_in$jac <- TRUE
    options$jac <- NULL
  }
  fun0 <- func
  if (isTRUE(mk_in$jac)) {
    mj <- .spsh_memoize_jac(func)
    fun0 <- mj$f
    mk_in$jac <- mj$jac
  }
  S$func <- function(x) as.numeric(fun0(x))
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  dim <- length(lb)
  lb[!is.finite(lb)] <- -1e50
  ub[!is.finite(ub)] <- 1e50
  if (any(lb > ub))
    stop("Error: lb > ub in bounds ", paste(ifelse(lb > ub, "True", "False"), collapse = ", "), ".", call. = FALSE)
  g_cons <- .spsh_g_cons(constraints, dim)
  # minimizer_kwargs of the local minimizations
  mk <- list(method = "SLSQP", bounds = list(lb = lb, ub = ub), options = list(), callback = callback)
  if (!is.null(mk_in)) mk <- utils::modifyList(mk, mk_in, keep.null = TRUE) else mk$options <- list(ftol = 1e-12)
  if ((tolower(mk$method) %in% c("slsqp", "cobyla", "cobyqa", "trust-constr") &&
       !is.null(mk_in) && is.null(mk_in$constraints) && !is.null(constraints)) || !is.null(g_cons))
    mk$constraints <- constraints
  S$minimize_every_iter <- TRUE
  S$local_iter <- FALSE
  S$symmetry_warn <- FALSE
  if (!is.null(options)) {
    mk$options <- utils::modifyList(mk$options, options, keep.null = TRUE)
    for (o in c("jac", "hess", "hessp")) {
      if (o %in% names(mk$options)) {
        mk[[o]] <- mk$options[[o]]
        mk$options[[o]] <- NULL
      }
    }
    S$minimize_every_iter <- options$minimize_every_iter %||% TRUE
    S$maxiter <- options$maxiter
    S$maxfev <- options$maxfev
    S$maxev <- options$maxev
    S$init <- as.numeric(Sys.time())
    S$maxtime <- options$maxtime
    if ("f_min" %in% names(options)) {
      S$f_min_true <- options$f_min
      S$f_tol <- options$f_tol %||% 0.0001
    }
    S$minhgrd <- options$minhgrd
    sym <- options$symmetry
    S$symmetry_warn <- length(sym) > 0L && !isFALSE(sym) && !identical(sym, 0) && !identical(sym, 0L) && dim > 1L
    S$local_iter <- options$local_iter %||% FALSE
  }
  meth <- tolower(mk$method)
  if (is.null(.spsh_solver_args[[meth]])) stop("KeyError: '", meth, "'", call. = FALSE)
  keep <- c("fun", "x0", "args", "callback", "options", "method", .spsh_solver_args[[meth]])
  mk <- mk[intersect(names(mk), keep)]
  mk$options <- mk$options[intersect(names(mk$options), c(keep, "ftol"))]
  S$takes_bounds <- "bounds" %in% keep
  S$mk <- mk
  S$stop_global <- FALSE
  S$iters <- iters
  S$iters_done <- 0L
  S$n <- n
  S$n_sampled <- 0
  S$fn <- 0L
  S$hgr <- 0L
  if (is.null(S$n) && is.null(S$iters)) S$n <- 2^dim + 1
  if (is.null(S$iters)) S$iters <- 1L
  if (!is.null(S$n) && S$n == 100) S$n <- 2^dim + 1
  if (!(is.null(S$maxiter) && is.null(S$maxfev) && is.null(S$maxev) && is.null(S$minhgrd) &&
        is.null(S$f_min_true))) S$iters <- NULL
  S$HC <- .spsh_complex(dim, lb, ub, S$func, g_cons, workers)
  S$lmc <- list()
  S$lmc_keys <- character()
  S$lmc_size <- 0L
  S$xl_maps <- list()
  S$f_maps <- numeric()
  S$nlfev <- 0L
  S$nljev <- 0L
  S$nlhev <- 0L
  S$res <- list()
  tryCatch(.spsh_iterate_simplicial(S), spsh_done = function(e) NULL)
  if (!S$minimize_every_iter) .spsh_find_minima(S)
  S$fn <- S$HC$nfev
  res <- list(x = S$res$x, fun = S$res$fun, xl = S$res$xl, funl = S$res$funl, nfev = S$res$nfev %||% 0L,
              nlfev = S$nlfev, nljev = S$nljev, nlhev = S$nlhev, nit = S$iters_done)
  if (!length(S$xl_maps)) {
    .spsh_find_lowest_vertex(S)
    res$fun <- S$f_lowest
    res$x <- S$x_lowest
    res$nfev <- S$fn
    res$tnev <- S$n_sampled
    res$success <- FALSE
    res$message <- paste0("Failed to find a feasible minimizer point. Lowest sampling point = ",
                          if (is.null(S$f_lowest)) "None" else .spsq_repr(S$f_lowest))
  } else {
    res$message <- "Optimization terminated successfully."
    res$success <- TRUE
  }
  res
}

# ScipyShgo: the options of OptimizationSciPy.jl (n = 100, iters = 1,
# simplicial sampling), the constraints lcons <= cons(x) <= ucons as one
# inequality per finite bound
.scipy_run_shgo <- function(f, lb, ub, cons, a, opts, kw) {
  constraints <- NULL
  if (!is.null(cons)) {
    constraints <- list()
    for (i in which(is.finite(a$lcons))) local({
      j <- i
      constraints[[length(constraints) + 1L]] <<- list(type = "ineq", fun = function(x) cons(x)[j] - a$lcons[j])
    })
    for (i in which(is.finite(a$ucons))) local({
      j <- i
      constraints[[length(constraints) + 1L]] <<- list(type = "ineq", fun = function(x) a$ucons[j] - cons(x)[j])
    })
  }
  args <- utils::modifyList(list(n = 100L, iters = 1L, sampling_method = "simplicial"), kw, keep.null = TRUE)
  res <- do.call(.sp_shgo, c(list(function(x) as.numeric(f(x)), lb, ub, constraints = constraints), args))
  if (is.null(res$x)) res$x <- rep(NaN, length(lb))
  if (is.null(res$fun)) res$fun <- NaN
  res$retcode <- if (res$success) "Success" else "Failure"
  res
}
