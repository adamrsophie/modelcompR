# Fixed noise paths of tools/stochasticdiffeq-probe.jl: the increments are
# z(k) = mod(37 k, 101) / 50 - 1, scaled by sqrt(hg) for W and by `zscale`
# for Z (the Z increments of the second process start at index 7919).
sdeprobe_z <- function(k) ((37 * k) %% 101) / 50 - 1

sdeprobe_grid <- function(m, hg, N, off = 0, mz = m, zscale = sqrt(hg)) {
  W <- matrix(0, N + 1L, m)
  Z <- matrix(0, N + 1L, mz)
  for (k in seq_len(N)) {
    W[k + 1L, ] <- W[k, ] + sqrt(hg) * sdeprobe_z(off + m * (k - 1L) + seq_len(m))
    Z[k + 1L, ] <- Z[k, ] + zscale * sdeprobe_z(off + 7919 + mz * (k - 1L) + seq_len(mz))
  }
  list(t = (0:N) * hg, W = W, Z = Z)
}
