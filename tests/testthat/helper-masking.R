# Functions of the package namespace whose names are exports of the packages
# R attaches by default, operators such as %||% included. Registered S3
# methods (solve.ODEProblem, `[.SciMLProblem`, ...) are not masks. It reads
# the installed namespace, so it runs on the built package as well as under
# devtools::load_all(); tools/check_masking.R is the source level version.
masking_report <- function(pkg = "modelcompR") {
  default_pkgs <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "datasets")
  ns <- asNamespace(pkg)
  s3 <- getNamespaceInfo(ns, "S3methods")
  methods <- if (length(s3)) s3[, 3L] else character()
  exported <- getNamespaceExports(ns)
  rows <- list()
  for (p in default_pkgs) {
    hits <- intersect(ls(ns, all.names = TRUE), getNamespaceExports(p))
    hits <- hits[vapply(hits, function(n) is.function(get(n, envir = ns, inherits = FALSE)), TRUE)]
    hits <- setdiff(hits, methods)
    if (length(hits))
      rows[[p]] <- data.frame(name = hits, masks = paste0(p, "::", hits), exported = hits %in% exported)
  }
  if (length(rows)) do.call(rbind, unname(rows)) else
    data.frame(name = character(), masks = character(), exported = logical())
}
