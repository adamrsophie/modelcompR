# Report function names that mask base R functions.
#
# Run from the repository root: Rscript tools/check_masking.R
# Scans R/*.R for top-level definitions. Exported names (roxygen @export)
# must not mask any function of the packages R attaches by default; internal
# names are reported too, since they shadow those functions inside the
# package namespace. S3 methods (solve.ODEProblem, plot.SciMLSolution, ...)
# are not masks and are skipped.

default_pkgs <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "datasets")
base_names <- unique(unlist(lapply(default_pkgs, getNamespaceExports)))
generics <- base_names[vapply(base_names, function(n) {
  f <- tryCatch(get(n, envir = baseenv()), error = function(e) NULL)
  is.function(f) && (is.primitive(f) || any(grepl("UseMethod|standardGeneric", deparse(f))))
}, TRUE)]

is_s3_method <- function(n) {
  parts <- strsplit(n, ".", fixed = TRUE)[[1L]]
  if (length(parts) < 2L) return(FALSE)
  for (i in seq_len(length(parts) - 1L)) {
    g <- paste(parts[seq_len(i)], collapse = ".")
    if (g %in% c(generics, "print", "format", "plot", "summary", "solve", "as.data.frame",
                 "as.matrix", "length", "names", "diff", "c", "Ops", "Math", "remake"))
      return(TRUE)
  }
  FALSE
}

rows <- list()
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  lines <- readLines(f, warn = FALSE)
  for (i in seq_along(lines)) {
    m <- regmatches(lines[i], regexec("^`?([A-Za-z.][A-Za-z0-9._]*)`?\\s*(<-|=)\\s*function", lines[i]))[[1L]]
    if (!length(m)) next
    n <- m[2L]
    if (!(n %in% base_names) || is_s3_method(n)) next
    # roxygen block directly above
    j <- i - 1L
    exported <- FALSE
    while (j >= 1L && grepl("^#'", lines[j])) {
      if (grepl("@export", lines[j])) exported <- TRUE
      j <- j - 1L
    }
    pkg <- default_pkgs[vapply(default_pkgs, function(p) n %in% getNamespaceExports(p), TRUE)][1L]
    rows[[length(rows) + 1L]] <- data.frame(name = n, masks = paste0(pkg, "::", n),
                                            exported = exported, file = paste0(f, ":", i))
  }
}
out <- if (length(rows)) do.call(rbind, rows) else data.frame()
if (!nrow(out)) {
  cat("No function masks a base R function.\n")
} else {
  print(out, row.names = FALSE)
  if (any(out$exported)) quit(status = 1L)
}
