# Regenerates R/solvers-tableaus.R. Extract Runge-Kutta tableau constants from the OrdinaryDiffEq.jl sources
# (CompiledFloats branches) into R source. Run from the repository root with
#   Rscript tools/solvers-gen-tableaus.R
# after cloning OrdinaryDiffEq.jl into documentation/.
lib <- "documentation/OrdinaryDiffEq.jl/lib"
out_file <- "R/solvers-tableaus.R"
# the commit the constants come from, recorded in the generated header
commit <- system2("git", c("-C", "documentation/OrdinaryDiffEq.jl", "rev-parse", "HEAD"), stdout = TRUE)
stopifnot(length(commit) == 1L, grepl("^[0-9a-f]{40}$", commit))

extract <- function(file, fname) {
  lines <- readLines(file.path(lib, file), warn = FALSE)
  starts <- grep(paste0("function ", fname, "\\("), lines)
  st <- NA
  for (s in starts) {
    first <- s + which(grepl("convert\\(", lines[(s + 1L):length(lines)]))[1L]
    if (any(grepl("CompiledFloats", lines[s:(first - 1L)]))) { st <- s; break }
  }
  if (is.na(st)) stop("no CompiledFloats method for ", fname)
  en <- st + which(lines[(st + 1L):length(lines)] == "end")[1L]
  body <- lines[st:en]
  body <- body[!grepl("^\\s*#", body)]
  m <- regmatches(body, regexec("^\\s*(\\w+)\\s*=\\s*convert\\((T|T2),\\s*(.+)\\)\\s*$", body))
  m <- m[lengths(m) > 0L]
  # the expressions (numeric literals and exact rationals a // b) are
  # evaluated by Julia, which rounds them correctly; R's parser misrounds
  # some 17 digit literals and has no exact rationals
  exprs <- vapply(m, `[`, "", 4L)
  stopifnot(all(grepl("^[-+0-9.eE/ ()]+$", exprs)))
  names(exprs) <- vapply(m, `[`, "", 2L)
  stopifnot(!anyDuplicated(names(exprs)))
  exprs
}

# repr of convert(Float64, expr) for every expression, in one Julia run
julia_values <- function(exprs) {
  jl <- tempfile(fileext = ".jl")
  writeLines(paste0("println(repr(convert(Float64, ", exprs, ")))"), jl)
  out <- system2("julia", c("--startup-file=no", jl), stdout = TRUE)
  stopifnot(length(out) == length(exprs))
  stats::setNames(out, names(exprs))
}

fmt <- function(v, name) {
  ent <- sprintf("  %s = %s", names(v), ifelse(as.numeric(v) == 0, "0", v))
  paste0(name, " <- c(\n", paste(ent, collapse = ",\n"), "\n)\n")
}

tabs <- list(
  list(".BS3_TAB", "OrdinaryDiffEqLowOrderRK/src/low_order_rk_tableaus.jl", "BS3ConstantCache"),
  list(".DP5_TAB", "OrdinaryDiffEqLowOrderRK/src/low_order_rk_tableaus.jl", "DP5ConstantCacheActual"),
  list(".DP5_DENSE", "OrdinaryDiffEqLowOrderRK/src/low_order_rk_tableaus.jl", "DP5_dense_ds"),
  list(".VERN6_TAB", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern6Tableau"),
  list(".VERN6_EXTRA", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern6ExtraStages"),
  list(".VERN6_INTERP", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern6InterpolationCoefficients"),
  list(".VERN7_TAB", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern7Tableau"),
  list(".VERN7_EXTRA", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern7ExtraStages"),
  list(".VERN7_INTERP", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern7InterpolationCoefficients"),
  list(".VERN8_TAB", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern8Tableau"),
  list(".VERN8_EXTRA", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern8ExtraStages"),
  list(".VERN8_INTERP", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern8InterpolationCoefficients"),
  list(".VERN9_TAB", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern9Tableau"),
  list(".VERN9_EXTRA", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern9ExtraStages"),
  list(".VERN9_INTERP", "OrdinaryDiffEqVerner/src/verner_tableaus.jl", "Vern9InterpolationCoefficients")
)
src <- c(
  "# Runge-Kutta tableau constants, extracted by a script (not typed by hand)",
  "# from the CompiledFloats methods of OrdinaryDiffEq.jl (generator:",
  "#   tools/solvers-gen-tableaus.R):",
  "#   lib/OrdinaryDiffEqLowOrderRK/src/low_order_rk_tableaus.jl",
  "#     BS3ConstantCache, DP5ConstantCacheActual, DP5_dense_ds",
  "#   lib/OrdinaryDiffEqVerner/src/verner_tableaus.jl",
  "#     Vern6Tableau, Vern6ExtraStages, Vern6InterpolationCoefficients,",
  "#     Vern7Tableau, Vern7ExtraStages, Vern7InterpolationCoefficients,",
  "#     Vern8Tableau, Vern8ExtraStages, Vern8InterpolationCoefficients,",
  "#     Vern9Tableau, Vern9ExtraStages, Vern9InterpolationCoefficients",
  paste0("# OrdinaryDiffEq.jl commit ", commit, " (pinned: regenerate"),
  "# from that commit, or update this line with the new one).",
  "# Names are the Julia field names. Do not edit by hand.",
  "")
for (tb in tabs) src <- c(src, fmt(julia_values(extract(tb[[2]], tb[[3]])), tb[[1]]))
writeLines(src, out_file)
cat("wrote", out_file, "; now run python3 tools/exact_literals.py\n")
