# Writes R/sysdata.rda (the object .STAB_TABLES) from the dump of
# tools/stabilized-dump.jl: the coefficient tables of the stabilized
# Runge-Kutta methods (ROCK2, ROCK4, ESERK4, ESERK5, SERK2) and of the
# stabilized SDE methods (SROCK1, SROCK2, SROCKEM, SROCKC2, SKSROCK,
# TangXiaoSROCK2, KomBurSROCK2). Field names are made ASCII.
#
# Run from the repository root:
#   Rscript tools/stabilized-gen.R tools/stabilized-tables.tsv
args <- commandArgs(TRUE)
infile <- if (length(args)) args[1L] else "tools/stabilized-tables.tsv"
greek <- c("ᵤ" = "u", "ₑ" = "e", "ᵢ" = "i", "σ" = "sigma", "τ" = "tau",
           "α" = "alpha", "n̂" = "nhat", "δ" = "delta", "η" = "eta")
ascii <- function(x) {
  for (k in names(greek)) x <- gsub(k, greek[[k]], x, fixed = TRUE)
  if (grepl("[^ -~]", x)) stop("non ASCII field name left: ", x)
  x
}
.STAB_TABLES <- list()
for (ln in readLines(infile, encoding = "UTF-8")) {
  parts <- strsplit(ln, "\t", fixed = TRUE)[[1L]]
  v <- as.numeric(strsplit(parts[4L], ",", fixed = TRUE)[[1L]])
  dims <- as.integer(strsplit(parts[3L], "x", fixed = TRUE)[[1L]])
  if (length(dims) == 2L) v <- matrix(v, dims[1L], dims[2L])
  .STAB_TABLES[[parts[1L]]][[ascii(enc2utf8(parts[2L]))]] <- v
}
save(.STAB_TABLES, file = "R/sysdata.rda", compress = "xz")
cat(length(.STAB_TABLES), "tables written to R/sysdata.rda\n")
