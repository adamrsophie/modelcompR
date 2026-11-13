# Shared helpers of the HiGHS 1.12.0 MIP solver port: pdqsort
# (extern/pdqsort/pdqsort.h, Copyright (c) 2021 Orson Peters, zlib license),
# the libstdc++ algorithms the MIP code calls (std::partition, the heap
# fallback of pdqsort), HighsHashHelpers::hash for one integer key
# (util/HighsHash.h), HighsIntegers (util/HighsIntegers.h) and the
# HighsCDouble operations R/optimization-scipy-highs-presolve.R does not
# define. MIT license notice of HiGHS in R/optimization-scipy-highs-util.R.
#
# pdqsort is provided 'as-is', without any express or implied warranty. In
# no event will the authors be held liable for any damages arising from the
# use of this software. Permission is granted to anyone to use this software
# for any purpose, including commercial applications, and to alter it and
# redistribute it freely, subject to the following restrictions: 1. The
# origin of this software must not be misrepresented; you must not claim
# that you wrote the original software. If you use this software in a
# product, an acknowledgment in the product documentation would be
# appreciated but is not required. 2. Altered source versions must be
# plainly marked as such, and must not be misrepresented as being the
# original software. 3. This notice may not be removed or altered from any
# source distribution. (This R version is an altered source version.)

# pdqsort(begin, end, comp) on the R vector or list v; less(a, b) is the
# strict weak order. The branchless partition of pdqsort is only chosen for
# std::less on arithmetic types, where every order of equal keys is the
# same, so the branchy partition gives identical results everywhere.
.hmip_pdqsort <- function(v, less, branchless = FALSE) {
  n <- length(v)
  if (n == 0L) return(v)
  is_list <- is.list(v)
  get <- if (is_list) function(i) v[[i]] else function(i) v[i]
  swap <- function(i, j) {
    if (is_list) { t <- v[[i]]; v[[i]] <<- v[[j]]; v[[j]] <<- t }
    else { t <- v[i]; v[i] <<- v[j]; v[j] <<- t }
  }
  put <- if (is_list) function(i, x) v[[i]] <<- x else function(i, x) v[i] <<- x
  # positions are 1-based; `end` is one past the last element
  insertion_sort <- function(begin, end, guarded) {
    if (begin == end) return(invisible())
    cur <- begin + 1L
    while (cur < end) {
      sift <- cur; sift_1 <- cur - 1L
      if (less(get(sift), get(sift_1))) {
        tmp <- get(sift)
        repeat {
          put(sift, get(sift_1)); sift <- sift - 1L
          if (guarded && sift == begin) break
          sift_1 <- sift_1 - 1L
          if (!less(tmp, get(sift_1))) break
        }
        put(sift, tmp)
      }
      cur <- cur + 1L
    }
  }
  partial_insertion_sort <- function(begin, end) {
    if (begin == end) return(TRUE)
    limit <- 0L
    cur <- begin + 1L
    while (cur < end) {
      sift <- cur; sift_1 <- cur - 1L
      if (less(get(sift), get(sift_1))) {
        tmp <- get(sift)
        repeat {
          put(sift, get(sift_1)); sift <- sift - 1L
          if (sift == begin) break
          sift_1 <- sift_1 - 1L
          if (!less(tmp, get(sift_1))) break
        }
        put(sift, tmp)
        limit <- limit + (cur - sift)
      }
      if (limit > 8L) return(FALSE)
      cur <- cur + 1L
    }
    TRUE
  }
  sort2 <- function(a, b) if (less(get(b), get(a))) swap(a, b)
  sort3 <- function(a, b, c) { sort2(a, b); sort2(b, c); sort2(a, b) }
  partition_right <- function(begin, end) {
    pivot <- get(begin)
    first <- begin; last <- end
    repeat { first <- first + 1L; if (!less(get(first), pivot)) break }
    if (first - 1L == begin) {
      while (first < last) { last <- last - 1L; if (less(get(last), pivot)) break }
    } else {
      repeat { last <- last - 1L; if (less(get(last), pivot)) break }
    }
    already <- first >= last
    while (first < last) {
      swap(first, last)
      repeat { first <- first + 1L; if (!less(get(first), pivot)) break }
      repeat { last <- last - 1L; if (less(get(last), pivot)) break }
    }
    pivot_pos <- first - 1L
    put(begin, get(pivot_pos)); put(pivot_pos, pivot)
    c(pivot_pos, already)
  }
  # partition_right_branchless: block partitioning with offset buffers
  partition_right_branchless <- function(begin, end) {
    pivot <- get(begin)
    first <- begin; last <- end
    repeat { first <- first + 1L; if (!less(get(first), pivot)) break }
    if (first - 1L == begin) {
      while (first < last) { last <- last - 1L; if (less(get(last), pivot)) break }
    } else {
      repeat { last <- last - 1L; if (less(get(last), pivot)) break }
    }
    already <- first >= last
    if (!already) {
      swap(first, last)
      first <- first + 1L
      offl <- integer(64L); offr <- integer(64L)
      lbase <- first; rbase <- last
      num_l <- 0L; num_r <- 0L; start_l <- 0L; start_r <- 0L
      while (first < last) {
        num_unknown <- last - first
        left_split <- if (num_l == 0L) (if (num_r == 0L) num_unknown %/% 2L else num_unknown) else 0L
        right_split <- if (num_r == 0L) num_unknown - left_split else 0L
        nl <- if (left_split >= 64L) 64L else left_split
        i <- 0L
        while (i < nl) {
          offl[num_l + 1L] <- i; i <- i + 1L
          if (!less(get(first), pivot)) num_l <- num_l + 1L
          first <- first + 1L
        }
        nr <- if (right_split >= 64L) 64L else right_split
        i <- 0L
        while (i < nr) {
          i <- i + 1L; offr[num_r + 1L] <- i
          last <- last - 1L
          if (less(get(last), pivot)) num_r <- num_r + 1L
        }
        num <- min(num_l, num_r)
        if (num_l == num_r) {
          k <- 0L
          while (k < num) {
            swap(lbase + offl[start_l + k + 1L], rbase - offr[start_r + k + 1L]); k <- k + 1L
          }
        } else if (num > 0L) {
          l <- lbase + offl[start_l + 1L]; r <- rbase - offr[start_r + 1L]
          tmp <- get(l); put(l, get(r))
          k <- 1L
          while (k < num) {
            l <- lbase + offl[start_l + k + 1L]; put(r, get(l))
            r <- rbase - offr[start_r + k + 1L]; put(l, get(r))
            k <- k + 1L
          }
          put(r, tmp)
        }
        num_l <- num_l - num; num_r <- num_r - num
        start_l <- start_l + num; start_r <- start_r + num
        if (num_l == 0L) { start_l <- 0L; lbase <- first }
        if (num_r == 0L) { start_r <- 0L; rbase <- last }
      }
      if (num_l > 0L) {
        while (num_l > 0L) {
          num_l <- num_l - 1L
          last <- last - 1L
          swap(lbase + offl[start_l + num_l + 1L], last)
        }
        first <- last
      }
      if (num_r > 0L) {
        while (num_r > 0L) {
          num_r <- num_r - 1L
          swap(rbase - offr[start_r + num_r + 1L], first)
          first <- first + 1L
        }
        last <- first
      }
    }
    pivot_pos <- first - 1L
    put(begin, get(pivot_pos)); put(pivot_pos, pivot)
    c(pivot_pos, already)
  }
  partition_left <- function(begin, end) {
    pivot <- get(begin)
    first <- begin; last <- end
    repeat { last <- last - 1L; if (!less(pivot, get(last))) break }
    if (last + 1L == end) {
      while (first < last) { first <- first + 1L; if (less(pivot, get(first))) break }
    } else {
      repeat { first <- first + 1L; if (less(pivot, get(first))) break }
    }
    while (first < last) {
      swap(first, last)
      repeat { last <- last - 1L; if (!less(pivot, get(last))) break }
      repeat { first <- first + 1L; if (less(pivot, get(first))) break }
    }
    pivot_pos <- last
    put(begin, get(pivot_pos)); put(pivot_pos, pivot)
    pivot_pos
  }
  # libstdc++ std::make_heap and std::sort_heap with a comparator
  push_heap <- function(first, hole, top, value) {
    parent <- (hole - 1L) %/% 2L
    while (hole > top && less(get(first + parent), value)) {
      put(first + hole, get(first + parent))
      hole <- parent
      parent <- (hole - 1L) %/% 2L
    }
    put(first + hole, value)
  }
  adjust_heap <- function(first, hole, len, value) {
    force(value)
    top <- hole
    second <- hole
    while (second < (len - 1L) %/% 2L) {
      second <- 2L * (second + 1L)
      if (less(get(first + second), get(first + second - 1L))) second <- second - 1L
      put(first + hole, get(first + second))
      hole <- second
    }
    if (bitwAnd(len, 1L) == 0L && second == (len - 2L) %/% 2L) {
      second <- 2L * (second + 1L)
      put(first + hole, get(first + second - 1L))
      hole <- second - 1L
    }
    push_heap(first, hole, top, value)
  }
  heap_sort <- function(begin, end) {
    len <- end - begin
    if (len >= 2L) {
      parent <- (len - 2L) %/% 2L
      repeat {
        adjust_heap(begin, parent, len, get(begin + parent))
        if (parent == 0L) break
        parent <- parent - 1L
      }
    }
    last <- end
    while (last - begin > 1L) {
      last <- last - 1L
      value <- get(last)
      put(last, get(begin))
      adjust_heap(begin, 0L, last - begin, value)
    }
  }
  loop <- function(begin, end, bad_allowed, leftmost) {
    repeat {
      size <- end - begin
      if (size < 24L) {
        insertion_sort(begin, end, leftmost)
        return(invisible())
      }
      s2 <- size %/% 2L
      if (size > 128L) {
        sort3(begin, begin + s2, end - 1L)
        sort3(begin + 1L, begin + (s2 - 1L), end - 2L)
        sort3(begin + 2L, begin + (s2 + 1L), end - 3L)
        sort3(begin + (s2 - 1L), begin + s2, begin + (s2 + 1L))
        swap(begin, begin + s2)
      } else sort3(begin + s2, begin, end - 1L)
      if (!leftmost && !less(get(begin - 1L), get(begin))) {
        begin <- partition_left(begin, end) + 1L
        next
      }
      pr <- if (branchless) partition_right_branchless(begin, end) else partition_right(begin, end)
      pivot_pos <- as.integer(pr[1L]); already <- as.logical(pr[2L])
      l_size <- pivot_pos - begin
      r_size <- end - (pivot_pos + 1L)
      if (l_size < size %/% 8L || r_size < size %/% 8L) {
        bad_allowed <- bad_allowed - 1L
        if (bad_allowed == 0L) {
          heap_sort(begin, end)
          return(invisible())
        }
        if (l_size >= 24L) {
          swap(begin, begin + l_size %/% 4L)
          swap(pivot_pos - 1L, pivot_pos - l_size %/% 4L)
          if (l_size > 128L) {
            swap(begin + 1L, begin + (l_size %/% 4L + 1L))
            swap(begin + 2L, begin + (l_size %/% 4L + 2L))
            swap(pivot_pos - 2L, pivot_pos - (l_size %/% 4L + 1L))
            swap(pivot_pos - 3L, pivot_pos - (l_size %/% 4L + 2L))
          }
        }
        if (r_size >= 24L) {
          swap(pivot_pos + 1L, pivot_pos + (1L + r_size %/% 4L))
          swap(end - 1L, end - r_size %/% 4L)
          if (r_size > 128L) {
            swap(pivot_pos + 2L, pivot_pos + (2L + r_size %/% 4L))
            swap(pivot_pos + 3L, pivot_pos + (3L + r_size %/% 4L))
            swap(end - 2L, end - (1L + r_size %/% 4L))
            swap(end - 3L, end - (2L + r_size %/% 4L))
          }
        }
      } else if (already && partial_insertion_sort(begin, pivot_pos) &&
                 partial_insertion_sort(pivot_pos + 1L, end)) {
        return(invisible())
      }
      loop(begin, pivot_pos, bad_allowed, leftmost)
      begin <- pivot_pos + 1L
      leftmost <- FALSE
    }
  }
  loop(1L, n + 1L, as.integer(floor(log2(n))), TRUE)
  v
}

# libstdc++ std::partition (the bidirectional algorithm random access
# iterators use); k is the offset of the returned iterator
# libstdc++ std::partial_sort(first, first + k, last, comp): __heap_select
# then __sort_heap on the first k elements
.hmip_std_partial_sort <- function(v, k, less) {
  n <- length(v)
  if (k <= 0L || n == 0L) return(v)
  is_list <- is.list(v)
  get <- if (is_list) function(i) v[[i]] else function(i) v[i]
  put <- if (is_list) function(i, x) v[[i]] <<- x else function(i, x) v[i] <<- x
  push_heap <- function(hole, top, value) {
    parent <- (hole - 1L) %/% 2L
    while (hole > top && less(get(1L + parent), value)) {
      put(1L + hole, get(1L + parent)); hole <- parent; parent <- (hole - 1L) %/% 2L
    }
    put(1L + hole, value)
  }
  adjust_heap <- function(hole, len, value) {
    force(value)
    top <- hole; second <- hole
    while (second < (len - 1L) %/% 2L) {
      second <- 2L * (second + 1L)
      if (less(get(1L + second), get(second))) second <- second - 1L
      put(1L + hole, get(1L + second)); hole <- second
    }
    if (bitwAnd(len, 1L) == 0L && second == (len - 2L) %/% 2L) {
      second <- 2L * (second + 1L)
      put(1L + hole, get(second)); hole <- second - 1L
    }
    push_heap(hole, top, value)
  }
  # __make_heap(first, middle)
  if (k >= 2L) {
    parent <- (k - 2L) %/% 2L
    repeat {
      adjust_heap(parent, k, get(1L + parent))
      if (parent == 0L) break
      parent <- parent - 1L
    }
  }
  pos <- k + 1L
  while (pos <= n) {
    if (less(get(pos), get(1L))) {
      value <- get(pos); put(pos, get(1L)); adjust_heap(0L, k, value)
    }
    pos <- pos + 1L
  }
  last <- k
  while (last > 1L) {
    value <- get(last); put(last, get(1L)); last <- last - 1L
    adjust_heap(0L, last, value)
  }
  v
}

.hmip_std_partition <- function(v, pred) {
  first <- 1L; last <- length(v) + 1L
  repeat {
    repeat {
      if (first == last) return(list(v = v, k = first - 1L))
      if (pred(v[[first]])) first <- first + 1L else break
    }
    last <- last - 1L
    repeat {
      if (first == last) return(list(v = v, k = first - 1L))
      if (!pred(v[[last]])) last <- last - 1L else break
    }
    t <- v[[first]]; v[[first]] <- v[[last]]; v[[last]] <- t
    first <- first + 1L
  }
}

# HighsHashHelpers::hash for keys of at most 8 bytes: bytes = (lo32, hi32),
# pair_hash<1>(lo, hi) ^ (pair_hash<0>(lo, hi) >> 32). A 32 bit HighsInt x
# is hash(x); an 8 byte key such as std::pair<HighsInt, HighsInt>(a, r) is
# hash(a, r) (the u64 r * 2^32 + a). Negative values are taken as unsigned.
.hmip_hash <- function(x, hi = 0) {
  if (x < 0) x <- x + 4294967296
  if (hi < 0) hi <- hi + 4294967296
  .hs_u64_xor(.hs_pair_hash(1L, x, hi), .hs_u64_shift(.hs_pair_hash(0L, x, hi), 32L))
}
# comparison of two u64 limb vectors: -1, 0 or 1
.hmip_u64_cmp <- function(a, b) {
  for (i in 4:1) {
    if (a[i] < b[i]) return(-1L)
    if (a[i] > b[i]) return(1L)
  }
  0L
}

# fractionality() of HighsIntegers.h: |x - std::round(x)|
# std::round (half away from zero); x - trunc(x) is exact for every double,
# unlike floor(abs(x) + 0.5) at 0.49999999999999994 and near 2^52
.hmip_round <- function(x) {
  t <- trunc(x)
  up <- which(abs(x - t) >= 0.5)
  t[up] <- t[up] + sign(x[up])
  t
}
.hmip_fractionality <- function(x) abs(x - .hmip_round(x))

# HighsIntegers::nearestInteger: (int64_t)(x + copysign(0.5, x))
.hmip_nearestInteger <- function(x) trunc(x + (if (x < 0 || (x == 0 && 1 / x < 0)) -0.5 else 0.5))

.hmip_gcd <- function(a, b) {
  a <- abs(a); b <- abs(b)
  if (a == 0) return(b)
  if (b == 0) return(a)
  repeat {
    h <- a %% b
    a <- b; b <- h
    if (b == 0) break
  }
  a
}

# HighsCDouble operations beyond those of the presolve port
.cd_sub <- function(a, v) .cd_add(a, -v)
.cd_sub_cd <- function(a, b) { r <- .cd_add(a, -b[1L]); r[2L] <- r[2L] - b[2L]; r }
.cd_add_cd <- function(a, b) { r <- .cd_add(a, b[1L]); r[2L] <- r[2L] + b[2L]; r }
.cd_mul_cd <- function(a, b) .cd_iadd(.cd_mul(a, b[1L]), c(a[1L] * b[2L], 0))
.cd_div_cd <- function(a, b) {
  vd <- b[1L] + b[2L]
  d <- c(a[1L] / vd, a[2L] / vd)
  cc <- .cd_sub_cd(.cd_mul_cd(d, b), a)
  cc <- cc / vd
  .cd_sub_cd(d, cc)
}
.cd_abs <- function(a) if (.cd_val(a) < 0) -a else a
.cd_floor <- function(x) {
  if (.cd_val(.cd_abs(x)) < 1) {
    v <- .cd_val(x)
    return(if (v == 0 || v > 0) c(0, 0) else c(-1, 0))
  }
  fx <- floor(.cd_val(x))
  .cd_two_sum(fx, floor(.cd_val(.cd_sub(x, fx))))
}
.cd_ceil <- function(x) {
  if (.cd_val(.cd_abs(x)) < 1) {
    v <- .cd_val(x)
    return(if (v == 0 || v < 0) c(0, 0) else c(1, 0))
  }
  cx <- ceiling(.cd_val(x))
  .cd_two_sum(cx, ceiling(.cd_val(.cd_sub(x, cx))))
}
.cd_sqrt <- function(v) {
  c0 <- sqrt(v[1L] + v[2L])
  if (c0 == 0) return(c(0, 0))
  res <- .cd_add(.cd_div(v, c0), c0)
  res * 0.5
}

# HighsIntegers::denominator(x, eps, maxdenom)
.hmip_denominator <- function(x, eps, maxdenom) {
  ai <- trunc(x)
  m <- c(ai, 1, 1, 0)
  xi <- c(x, 0)
  fraction <- .cd_sub(xi, ai)
  while (.cd_val(fraction) > eps) {
    xi <- .cd_div_cd(c(1, 0), fraction)
    if (.cd_val(xi) > 2^53) break
    ai <- trunc(.cd_val(xi))
    t <- m[3L] * ai + m[4L]
    if (t > maxdenom) break
    m[4L] <- m[3L]; m[3L] <- t
    t <- m[1L] * ai + m[2L]
    m[2L] <- m[1L]; m[1L] <- t
    fraction <- .cd_sub(xi, ai)
  }
  ai <- trunc((maxdenom - m[4L]) / m[3L])
  m[2L] <- m[2L] + m[1L] * ai
  m[4L] <- m[4L] + m[3L] * ai
  x0 <- m[1L] / m[3L]
  x1 <- m[2L] / m[4L]
  x <- abs(x)
  if (abs(x - x0) < abs(x - x1)) m[3L] else m[4L]
}

# std::frexp exponent
.hmip_frexp_exp <- function(x) {
  if (x == 0 || !is.finite(x)) return(0L)
  e <- floor(log2(abs(x))) + 1
  # correct for rounding of log2 near powers of two
  if (abs(x) / 2^e >= 1) e <- e + 1
  if (abs(x) / 2^e < 0.5) e <- e - 1
  as.integer(e)
}

# HighsIntegers::integralScale(vals, deltadown, deltaup)
.hmip_integralScale <- function(vals, deltadown, deltaup) {
  n <- length(vals)
  if (n == 0L) return(0)
  # std::minmax_element by |value|: first smallest, last largest
  imin <- 1L; imax <- 1L
  if (n >= 2L) {
    # the pairwise scheme of minmax_element gives the first minimum and the
    # last maximum
    a <- abs(vals)
    imin <- which(a == min(a))[1L]
    imax <- max(which(a == max(a)))
  }
  minval <- vals[imin]; maxval <- vals[imax]
  expshift <- 0L
  if (minval < -deltadown || minval > deltaup) expshift <- .hmip_frexp_exp(minval)
  expshift <- max(-expshift, 0L) + 3L
  expMaxVal <- min(.hmip_frexp_exp(maxval), 32L)
  if (expMaxVal + expshift > 32L) expshift <- 32L - expMaxVal
  denom <- 75 * 2^expshift
  startdenom <- denom
  val <- .cd_mul(c(vals[1L], 0), startdenom)
  downval <- .cd_floor(.cd_add(val, deltaup))
  fraction <- .cd_sub_cd(val, downval)
  if (.cd_val(fraction) > deltadown) {
    denom <- denom * .hmip_denominator(.cd_val(fraction), deltaup, 1000)
    val <- .cd_mul(c(vals[1L], 0), denom)
    downval <- .cd_floor(.cd_add(val, deltaup))
    fraction <- .cd_sub_cd(val, downval)
    if (.cd_val(fraction) > deltadown) return(0)
  }
  currgcd <- abs(.cd_val(downval))
  if (n >= 2L) for (i in 2:n) {
    val <- .cd_mul(c(vals[i], 0), denom)
    downval <- .cd_floor(.cd_add(val, deltaup))
    fraction <- .cd_sub_cd(val, downval)
    if (.cd_val(fraction) > deltadown) {
      val <- .cd_mul(c(vals[i], 0), startdenom)
      fraction <- .cd_sub_cd(val, .cd_floor(val))
      denom <- denom * .hmip_denominator(.cd_val(fraction), deltaup, 1000)
      val <- .cd_mul(c(vals[i], 0), denom)
      downval <- .cd_floor(.cd_add(val, deltaup))
      fraction <- .cd_sub_cd(val, downval)
      if (.cd_val(fraction) > deltadown) return(0)
    }
    if (currgcd != 1) {
      currgcd <- .hmip_gcd(currgcd, trunc(.cd_val(downval)))
      if (denom > 4294967295) {
        denom <- denom %/% currgcd
        if (startdenom != 1) startdenom <- startdenom %/% .hmip_gcd(currgcd, startdenom)
        currgcd <- 1
      }
    }
  }
  denom / currgcd
}

# std::nextafter(x, -inf)
.hmip_nextafter_down <- function(x) {
  if (is.na(x) || x == -Inf) return(x)
  if (x == 0) return(-4.9406564584124654e-324)
  if (x == Inf) return(.Machine$double.xmax)
  e <- floor(log2(abs(x)))
  if (2^e > abs(x)) e <- e - 1
  ulp <- if (x > 0 && abs(x) == 2^e) 2^(e - 53) else 2^(e - 52)
  ulp <- max(ulp, 4.9406564584124654e-324)
  x - ulp
}
