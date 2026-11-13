# The random number generator of HiGHS 1.12.0 (util/HighsRandom.h and the
# pair hash of util/HighsHash.h), the version SciPy 1.18.1 vendors in
# subprojects/highs, built with 32 bit HighsInt (highsint64=false).
#
# HiGHS is Copyright (c) 2016-2025 Julian Hall, Ivet Galabova, Leona
# Gottwald, Qi Huangfu, Michael Feldmeier and the HiGHS contributors,
# released under the MIT license:
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the
# following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Unsigned 64 bit integers are vectors of four 16 bit limbs, least
# significant first, held in doubles.

.hs_u64 <- function(hi, lo) c(lo %% 65536, lo %/% 65536, hi %% 65536, hi %/% 65536)
.hs_hi32 <- function(x) x[3L] + 65536 * x[4L]
.hs_lo32 <- function(x) x[1L] + 65536 * x[2L]

.hs_u64_hex <- function(h) {
  d <- strtoi(strsplit(h, "")[[1L]], 16L)
  limb <- function(k) sum(d[(13L - 4L * k):(16L - 4L * k)] * c(4096, 256, 16, 1))
  c(limb(0L), limb(1L), limb(2L), limb(3L))
}

.hs_carry <- function(x) {
  for (i in 1:3) {
    q <- x[i] %/% 65536
    x[i] <- x[i] - q * 65536
    x[i + 1L] <- x[i + 1L] + q
  }
  x[4L] <- x[4L] %% 65536
  x
}

.hs_u64_add <- function(a, b) .hs_carry(a + b)

.hs_u64_mul <- function(a, b) {
  r <- numeric(4L)
  for (i in 1:4) for (j in 1:(5L - i)) r[i + j - 1L] <- r[i + j - 1L] + a[i] * b[j]
  # each sum of up to four 32 bit products fits a double exactly
  .hs_carry(r)
}

.hs_u64_xor <- function(a, b) bitwXor(as.integer(a), as.integer(b)) + 0

# shift right (k > 0) or left (k < 0) by |k| bits, modulo 2^64
.hs_u64_shift <- function(x, k) {
  v <- 0
  bits <- integer(64L)
  for (i in 1:4) bits[(16L * (i - 1L) + 1L):(16L * i)] <- as.integer(intToBits(as.integer(x[i]))[1:16])
  out <- integer(64L)
  if (k > 0) {
    if (k < 64) out[1:(64 - k)] <- bits[(k + 1):64]
  } else {
    k <- -k
    if (k < 64) out[(k + 1):64] <- bits[1:(64 - k)]
  }
  vapply(1:4, function(i) sum(out[(16L * (i - 1L) + 1L):(16L * i)] * 2^(0:15)), 0)
}

.hs_hash_c <- c(
  "c8497d2a400d9551", "80c8963be3e4c2f3", "042d8680e260ae5b", "8a183895eeac1536",
  "a94e9c75f80ad6de", "7e92251dec62835e", "07294165cb671455", "89b0f6212b0a4292",
  "31900011b96bf554", "a44540f8eee2094f", "ce7ffd372e4c64fc", "51c9d471bfe6a10f",
  "758c2a674483826f", "f91a20abe63f8b02", "c2a069024a1fcc6f", "d5bb18b70c5dbd59",
  "d510adac6d1ae289", "571d069b23050a79", "60873b8872933e06", "780481cc19670350",
  "7a48551760216885", "b5d68b918231e6ca", "a7e5571699aa5274", "7b6d309b2cfdcf01",
  "04e77c3d474daeff", "4dbf099fd7247031", "5d70dca901130beb", "9f8b5f0df4182499",
  "293a74c9686092da", "d09bdab6840f52b3", "c05d47f3ab302263", "6b79e62b884b65d6",
  "a581106fc980c34d", "f081b7145ea2293e", "fb27243dd7c3f5ad", "5211bf8860ea667f",
  "9455e65cb2385e7f", "0dfaf6731b449b33", "4ec98b3c6f5e68c7", "007bfd4a42ae936b",
  "65c93061f8674518", "640816f17127c5d1", "6dd4bab17b7c3a74", "34d9268c256fa1ba",
  "0b4d0c6b5b50d7f4", "30aa965bc9fadaff", "c0ac1d0c2771404d", "c5e64509abb76ef2",
  "d606b11990624a36", "0d3f05d242ce2fb7", "469a803cb276fe32", "a4a44d177a3e23f4",
  "b9d9a120dcc1ca03", "2e15af8165234a2e", "10609ba2720573d4", "aa4191b60368d1d5",
  "333dd2300bc57762", "df6ec48f79fb402f", "5ed20fcef1b734fa", "4c94924ec8be21ee",
  "5abe6ad9d131e631", "be10136a522e602d", "53671115c340e779", "9f392fe43e2144da")
.hs_hash_cl <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- lapply(.hs_hash_c, .hs_u64_hex)
    cache
  }
})

# pair_hash<k>(a, b) = (a + c[2k]) * (b + c[2k + 1]) modulo 2^64 for 32 bit
# a and b
.hs_pair_hash <- function(k, a, b) {
  cl <- .hs_hash_cl()
  .hs_u64_mul(.hs_u64_add(.hs_u64(0, a), cl[[2L * k + 1L]]), .hs_u64_add(.hs_u64(0, b), cl[[2L * k + 2L]]))
}

# Arithmetic modulo the Mersenne prime M61 = 2^61 - 1 (HighsHashHelpers) on
# limb vectors: values are reduced into [0, M61), which is what
# multiply_modM61 and the folds after each addition compute.
.hs_limbs_bits <- function(x) {
  unlist(lapply(x, function(l) as.integer(intToBits(as.integer(l))[1:16])))
}
.hs_bits_limbs <- function(b) {
  b <- c(b, integer((-length(b)) %% 16L))
  vapply(seq_len(length(b) %/% 16L), function(i) sum(b[(16L * (i - 1L) + 1L):(16L * i)] * 2^(0:15)), 0)
}
.hs_m61_fold <- function(limbs) {
  repeat {
    b <- .hs_limbs_bits(limbs)
    if (length(b) <= 61L || !any(b[62:length(b)] != 0L)) break
    lo <- .hs_bits_limbs(b[1:61]); hi <- .hs_bits_limbs(b[62:length(b)])
    limbs <- .hs_carry(c(lo, numeric(4L))[1:4] + c(hi, numeric(4L))[1:4])
  }
  limbs <- c(limbs, numeric(4L))[1:4]
  if (all(limbs == c(65535, 65535, 65535, 8191))) numeric(4L) else limbs
}
.hs_m61_add <- function(a, b) .hs_m61_fold(.hs_carry(a + b))
.hs_m61_mul <- function(a, b) {
  r <- numeric(8L)
  for (i in 1:4) for (j in 1:4) r[i + j - 1L] <- r[i + j - 1L] + a[i] * b[j]
  for (i in 1:7) {
    q <- r[i] %/% 65536
    r[i] <- r[i] - q * 65536
    r[i + 1L] <- r[i + 1L] + q
  }
  .hs_m61_fold(r)
}
.hs_m61_modexp <- function(a, e) {
  # the loop of modexp_M61 as written: it squares first and tests the low
  # bit of the remaining exponent, so it is a^e only for some e; the hash
  # values follow it exactly
  res <- a
  while (e != 1) {
    res <- .hs_m61_mul(res, res)
    if (e %% 2 == 1) res <- .hs_m61_mul(res, a)
    e <- e %/% 2
  }
  res
}
# HighsHashHelpers::sparse_combine(hash, index, value)
.hs_sparse_combine <- function(hash, index, value) {
  cl <- .hs_hash_cl()
  a <- .hs_m61_fold_bits(cl[[bitwAnd(index, 63L) + 1L]])
  v <- .hs_m61_fold_bits(.hs_u64(0, value * 2 + 1))
  .hs_m61_add(hash, .hs_m61_mul(v, .hs_m61_modexp(a, index %/% 64L + 1L)))
}
# x & M61: clear the top three bits
.hs_m61_fold_bits <- function(x) {
  x[4L] <- x[4L] %% 8192
  x
}

# the value of the top nbits bits of x
.hs_u64_top <- function(x, nbits) {
  v <- .hs_u64_shift(x, 64L - nbits)
  v[1L] + 65536 * (v[2L] + 65536 * (v[3L] + 65536 * v[4L]))
}

.hs_rng <- function(seed = 0) {
  r <- new.env(parent = emptyenv())
  st <- .hs_u64(0, seed)
  repeat {
    st <- .hs_pair_hash(0L, .hs_lo32(st), .hs_hi32(st))
    st <- .hs_u64_xor(st, .hs_u64_shift(.hs_pair_hash(1L, .hs_hi32(st), seed), 32L))
    if (any(st != 0)) break
  }
  r$state <- st
  r
}

.hs_rng_advance <- function(r) {
  s <- r$state
  s <- .hs_u64_xor(s, .hs_u64_shift(s, 12L))
  s <- .hs_u64_xor(s, .hs_u64_shift(s, -25L))
  s <- .hs_u64_xor(s, .hs_u64_shift(s, 27L))
  r$state <- s
}

# HighsRandom::integer(sup) with 32 bit HighsUInt: drawUniform over the
# outputs 0..7, 9..31 of the pair hash
.hs_rng_integer <- function(r, sup) {
  if (sup <= 1) return(0L)
  nbits <- floor(log2(sup - 1)) + 1
  while (2^(nbits - 1) > sup - 1) nbits <- nbits - 1
  while (2^nbits <= sup - 1) nbits <- nbits + 1
  repeat {
    .hs_rng_advance(r)
    lo <- .hs_lo32(r$state); hi <- .hs_hi32(r$state)
    for (k in c(0:7, 9:31)) {
      val <- .hs_u64_top(.hs_pair_hash(k, lo, hi), nbits)
      if (val < sup) return(as.integer(val))
    }
  }
}

.hs_rng_fraction <- function(r) {
  .hs_rng_advance(r)
  lo <- .hs_lo32(r$state); hi <- .hs_hi32(r$state)
  a <- .hs_u64_shift(.hs_pair_hash(0L, lo, hi), 12L)
  b <- .hs_u64_shift(.hs_pair_hash(1L, lo, hi), 38L)
  x <- .hs_u64_xor(a, b)
  out <- x[1L] + 65536 * (x[2L] + 65536 * (x[3L] + 65536 * x[4L]))
  (1 + out) * 0x1.ffffffffffffep-53
}

.hs_rng_shuffle <- function(r, v) {
  N <- length(v)
  if (N > 1L) for (i in N:2L) {
    pos <- .hs_rng_integer(r, i) + 1L
    tmp <- v[pos]; v[pos] <- v[i]; v[i] <- tmp
  }
  v
}
