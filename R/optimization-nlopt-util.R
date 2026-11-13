# The utilities of the NLopt library (src/util of NLopt 2.11.0): the stopping
# criteria (stop.c), the Mersenne Twister and the random numbers derived from
# it (mt19937ar.c), the seconds timer (timer.c), the rescaling helpers
# (rescale.c) and the Bentley McIlroy quicksort (qsort_r.c), which is also
# the BSD qsort_r that nlopt_qsort_r calls on macOS.
#
# NLopt is Copyright (c) 2007-2014 Massachusetts Institute of Technology,
# released under the MIT license:
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions: The
# above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED
# "AS IS", WITHOUT WARRANTY OF ANY KIND.
#
# The Mersenne Twister is Copyright (C) 1997 2002, Makoto Matsumoto and
# Takuji Nishimura, under the BSD license of mt19937ar.c.

# ---------------------------------------------------------------------------
# Unsigned 32 bit arithmetic on doubles

.nl32_and <- function(a, b) {
  bitwAnd(as.integer(a %/% 65536), as.integer(b %/% 65536)) * 65536 +
    bitwAnd(as.integer(a %% 65536), as.integer(b %% 65536))
}
.nl32_xor <- function(a, b) {
  bitwXor(as.integer(a %/% 65536), as.integer(b %/% 65536)) * 65536 +
    bitwXor(as.integer(a %% 65536), as.integer(b %% 65536))
}
.nl32_or <- function(a, b) {
  bitwOr(as.integer(a %/% 65536), as.integer(b %/% 65536)) * 65536 +
    bitwOr(as.integer(a %% 65536), as.integer(b %% 65536))
}
# (a * b) mod 2^32 for a, b below 2^32, exact in doubles by splitting a
.nl32_mul <- function(a, b) {
  ah <- a %/% 65536; al <- a %% 65536
  ((al * b) %% 4294967296 + ((ah * b) %% 65536) * 65536) %% 4294967296
}

# ---------------------------------------------------------------------------
# The Mersenne Twister of NLopt (one state per R session, as NLopt keeps one
# per thread)

.nlopt_mt <- new.env(parent = emptyenv())
.nlopt_mt$mt <- NULL
.nlopt_mt$mti <- 625L
.nlopt_mt$srand_called <- FALSE

.nlopt_init_genrand <- function(s) {
  mt <- numeric(624L)
  mt[1L] <- s %% 4294967296
  for (i in 2L:624L) {
    prev <- mt[i - 1L]
    mt[i] <- (.nl32_mul(.nl32_xor(prev, prev %/% 1073741824), 1812433253) + (i - 1L)) %% 4294967296
  }
  .nlopt_mt$mt <- mt
  .nlopt_mt$mti <- 624L
  invisible(NULL)
}

# the twist of mt19937ar.c genrand_int32 over the 624 words, in the three
# vectorizable blocks of its two loops and the last word
.nlopt_mt_twist <- function(mt) {
  mag <- function(y) ifelse(y %% 2 == 1, 2567483615, 0)
  tw <- function(k, k1, km) {
    y <- .nl32_or(.nl32_and(mt[k], 2147483648), .nl32_and(mt[k1], 2147483647))
    .nl32_xor(.nl32_xor(mt[km], y %/% 2), mag(y))
  }
  k <- 1:227
  mt[k] <- tw(k, k + 1L, k + 397L)
  k <- 228:454
  mt[k] <- tw(k, k + 1L, k - 227L)
  k <- 455:623
  mt[k] <- tw(k, k + 1L, k - 227L)
  mt[624L] <- tw(624L, 1L, 397L)
  mt
}

.nlopt_genrand_int32 <- function() {
  if (.nlopt_mt$mti >= 624L) {
    if (is.null(.nlopt_mt$mt)) .nlopt_init_genrand(5489)
    .nlopt_mt$mt <- .nlopt_mt_twist(.nlopt_mt$mt)
    .nlopt_mt$mti <- 0L
  }
  .nlopt_mt$mti <- .nlopt_mt$mti + 1L
  y <- .nlopt_mt$mt[.nlopt_mt$mti]
  y <- .nl32_xor(y, y %/% 2048)
  y <- .nl32_xor(y, .nl32_and((y * 128) %% 4294967296, 2636928640))
  y <- .nl32_xor(y, .nl32_and((y * 32768) %% 4294967296, 4022730752))
  .nl32_xor(y, y %/% 262144)
}

.nlopt_genrand_res53 <- function() {
  a <- .nlopt_genrand_int32() %/% 32
  b <- .nlopt_genrand_int32() %/% 64
  (a * 67108864 + b) * (1 / 9007199254740992)
}

.nlopt_urand <- function(a, b) a + (b - a) * .nlopt_genrand_res53()
.nlopt_iurand <- function(n) as.integer(.nlopt_genrand_int32() %% n)
.nlopt_nrand <- function(mean, stddev) {
  repeat {
    v1 <- .nlopt_urand(-1, 1)
    v2 <- .nlopt_urand(-1, 1)
    s <- v1 * v1 + v2 * v2
    if (s < 1) break
  }
  if (s == 0) mean else mean + v1 * sqrt(-2 * log(s) / s) * stddev
}

.nlopt_srand <- function(seed) {
  .nlopt_mt$srand_called <- TRUE
  .nlopt_init_genrand(as.numeric(seed) %% 4294967296)
}
.nlopt_srand_time <- function() {
  # nlopt_time_seed: seconds xor microseconds of the clock
  t <- as.numeric(Sys.time())
  .nlopt_srand(.nl32_xor(floor(t) %% 4294967296, floor((t %% 1) * 1e6)))
}
.nlopt_srand_time_default <- function() if (!.nlopt_mt$srand_called) .nlopt_srand_time()

# ---------------------------------------------------------------------------
# Timer and stopping criteria (timer.c, stop.c)

.nlopt_seconds <- function() proc.time()[["elapsed"]]

.nlopt_isinf <- function(x) abs(x) >= .Machine$double.xmax * 0.99 || is.infinite(x)
.nlopt_istiny <- function(x) x == 0 || (is.finite(x) && abs(x) < 0x1.0000000000000p-1022)

.nlopt_vector_norm <- function(vec, w, smin = NULL, smax = NULL) {
  if (!is.null(smin) && !is.null(smax)) vec <- smin + vec * (smax - smin)
  if (!is.null(w)) sum(w * abs(vec)) else sum(abs(vec))
}
.nlopt_diff_norm <- function(x, oldx, w, smin = NULL, smax = NULL) {
  if (!is.null(smin) && !is.null(smax)) {
    x <- smin + x * (smax - smin)
    oldx <- smin + oldx * (smax - smin)
  }
  if (!is.null(w)) sum(w * abs(x - oldx)) else sum(abs(x - oldx))
}

.nlopt_relstop <- function(vold, vnew, reltol, abstol) {
  if (.nlopt_isinf(vold)) return(FALSE)
  abs(vnew - vold) < abstol || abs(vnew - vold) < reltol * (abs(vnew) + abs(vold)) * 0.5 ||
    (reltol > 0 && vnew == vold)
}

.nlopt_stop_ftol <- function(s, f, oldf) .nlopt_relstop(oldf, f, s$ftol_rel, s$ftol_abs)
.nlopt_stop_f <- function(s, f, oldf) f <= s$minf_max || .nlopt_stop_ftol(s, f, oldf)
.nlopt_stop_x <- function(s, x, oldx) {
  if (.nlopt_diff_norm(x, oldx, s$x_weights) < s$xtol_rel * .nlopt_vector_norm(x, s$x_weights)) return(TRUE)
  if (is.null(s$xtol_abs)) return(FALSE)
  all(abs(x - oldx) < s$xtol_abs)
}
.nlopt_stop_dx <- function(s, x, dx) {
  if (.nlopt_vector_norm(dx, s$x_weights) < s$xtol_rel * .nlopt_vector_norm(x, s$x_weights)) return(TRUE)
  if (is.null(s$xtol_abs)) return(FALSE)
  all(abs(dx) < s$xtol_abs)
}
.nlopt_stop_xs <- function(s, xs, oldxs, smin, smax) {
  if (.nlopt_diff_norm(xs, oldxs, s$x_weights, smin, smax) <
      s$xtol_rel * .nlopt_vector_norm(xs, s$x_weights, smin, smax)) return(TRUE)
  if (is.null(s$xtol_abs)) return(FALSE)
  all(abs((smin + xs * (smax - smin)) - (smin + oldxs * (smax - smin))) < s$xtol_abs)
}
.nlopt_stop_evals <- function(s) s$maxeval > 0 && s$nevals >= s$maxeval
.nlopt_stop_time_ <- function(start, maxtime) maxtime > 0 && .nlopt_seconds() - start >= maxtime
.nlopt_stop_time <- function(s) .nlopt_stop_time_(s$start, s$maxtime)
.nlopt_stop_evalstime <- function(s) .nlopt_stop_evals(s) || .nlopt_stop_time(s)
.nlopt_stop_forced <- function(s) isTRUE(s$force_stop$value)

# a stopping environment: the fields of nlopt_stopping, with the evaluation
# count kept in `nevals` (nevals_p) and the force stop flag shared
.nlopt_stopping <- function(n, minf_max = -Inf, ftol_rel = 0, ftol_abs = 0, xtol_rel = 0, xtol_abs = NULL,
                            x_weights = NULL, maxeval = 0L, maxtime = 0, force_stop = NULL) {
  s <- new.env(parent = emptyenv())
  s$n <- n; s$minf_max <- minf_max; s$ftol_rel <- ftol_rel; s$ftol_abs <- ftol_abs; s$xtol_rel <- xtol_rel
  s$xtol_abs <- xtol_abs; s$x_weights <- x_weights; s$nevals <- 0L; s$maxeval <- as.integer(maxeval)
  s$maxtime <- maxtime; s$start <- .nlopt_seconds()
  s$force_stop <- force_stop %||% local({ e <- new.env(parent = emptyenv()); e$value <- FALSE; e })
  s$msg <- NULL
  s
}

.nlopt_count_constraints <- function(cs) sum(vapply(cs, function(c) c$m, 0L))

# ---------------------------------------------------------------------------
# Rescaling (rescale.c)

.nlopt_compute_rescaling <- function(n, dx) {
  s <- rep(1, n)
  if (n == 1L) return(s)
  i <- 2L
  while (i <= n && dx[i] == dx[i - 1L]) i <- i + 1L
  if (i <= n) s[2:n] <- dx[2:n] / dx[1L]
  s
}
.nlopt_rescale <- function(s, x) if (is.null(s)) x else x / s
.nlopt_unscale <- function(s, x) if (is.null(s)) x else x * s
.nlopt_reorder_bounds <- function(lb, ub) {
  sw <- lb > ub
  t <- lb[sw]; lb[sw] <- ub[sw]; ub[sw] <- t
  list(lb = lb, ub = ub)
}

# ---------------------------------------------------------------------------
# qsort_r_fallback of qsort_r.c (Bentley and McIlroy, "Engineering a Sort
# Function") on a vector `a` with a three way comparison `cmp(x, y)`;
# returns the sorted vector. Not stable: ties may come out in another order
# than glibc's merge sort (the qsort_r of NLopt on Linux) gives.

.nlopt_qsort <- function(a, cmp) {
  med3 <- function(i, j, k) {
    if (cmp(a[i], a[j]) < 0) {
      if (cmp(a[j], a[k]) < 0) j else if (cmp(a[i], a[k]) < 0) k else i
    } else {
      if (cmp(a[j], a[k]) > 0) j else if (cmp(a[i], a[k]) < 0) i else k
    }
  }
  swap <- function(i, j) { t <- a[i]; a[i] <<- a[j]; a[j] <<- t }
  vecswap <- function(i, j, n) if (n > 0) for (q in seq_len(n) - 1L) swap(i + q, j + q)
  insertion <- function(lo, n) {
    if (n < 2L) return(invisible())
    for (pm in (lo + 1L):(lo + n - 1L)) {
      pl <- pm
      while (pl > lo && cmp(a[pl - 1L], a[pl]) > 0) { swap(pl, pl - 1L); pl <- pl - 1L }
    }
  }
  rec <- function(lo, n) {
    repeat {
      if (n < 7L) { insertion(lo, n); return(invisible()) }
      pm <- lo + n %/% 2L
      if (n > 7L) {
        pl <- lo; pn <- lo + n - 1L
        if (n > 40L) {
          d <- n %/% 8L
          pl <- med3(pl, pl + d, pl + 2L * d)
          pm <- med3(pm - d, pm, pm + d)
          pn <- med3(pn - 2L * d, pn - d, pn)
        }
        pm <- med3(pl, pm, pn)
      }
      swap(lo, pm)
      pa <- pb <- lo + 1L
      pc <- pd <- lo + n - 1L
      swap_cnt <- FALSE
      repeat {
        while (pb <= pc && (r <- cmp(a[pb], a[lo])) <= 0) {
          if (r == 0) { swap_cnt <- TRUE; swap(pa, pb); pa <- pa + 1L }
          pb <- pb + 1L
        }
        while (pb <= pc && (r <- cmp(a[pc], a[lo])) >= 0) {
          if (r == 0) { swap_cnt <- TRUE; swap(pc, pd); pd <- pd - 1L }
          pc <- pc - 1L
        }
        if (pb > pc) break
        swap(pb, pc); swap_cnt <- TRUE
        pb <- pb + 1L; pc <- pc - 1L
      }
      if (!swap_cnt) { insertion(lo, n); return(invisible()) }
      pn <- lo + n
      d1 <- min(pa - lo, pb - pa)
      vecswap(lo, pb - d1, d1)
      d1 <- min(pd - pc, pn - pd - 1L)
      vecswap(pb, pn - d1, d1)
      d1 <- pb - pa
      d2 <- pd - pc
      if (d1 <= d2) {
        if (d1 > 1L) rec(lo, d1)
        if (d2 > 1L) { lo <- pn - d2; n <- d2 } else return(invisible())
      } else {
        if (d2 > 1L) rec(pn - d2, d2)
        if (d1 > 1L) n <- d1 else return(invisible())
      }
    }
  }
  if (length(a) > 1L) rec(1L, length(a))
  a
}

# indices 1..length(key) sorted by increasing key with .nlopt_qsort
.nlopt_qsort_order <- function(key)
  .nlopt_qsort(seq_along(key), function(i, j) if (key[i] < key[j]) -1 else if (key[i] > key[j]) 1 else 0)
