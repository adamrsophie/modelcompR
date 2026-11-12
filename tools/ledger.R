# Coverage of the ModelingToolkit public API.
#
# Run from the repository root after roxygen2 has regenerated NAMESPACE:
#   Rscript tools/ledger.R
# Each MTK name in tools/mtk_exports.tsv is matched to an export of modelcompR
# (exact name, without the macro "@" or the mutating "!", or through the
# renames below), or to an entry of tools/not_applicable.tsv. Writes
# tools/ledger_status.tsv and prints the summary.

ledger <- utils::read.delim("tools/mtk_exports.tsv", stringsAsFactors = FALSE)
ns <- readLines("NAMESPACE")
ex <- gsub("[`\"]", "", sub("^export[(](.*)[)]$", "\\1", grep("^export[(]", ns, value = TRUE)))
stubs <- if (file.exists("tools/stubs.tsv"))
  utils::read.delim("tools/stubs.tsv", stringsAsFactors = FALSE) else data.frame(name = character())
na <- if (file.exists("tools/not_applicable.tsv"))
  utils::read.delim("tools/not_applicable.tsv", stringsAsFactors = FALSE) else
  data.frame(name = character(), reason = character())

# MTK name -> R name when they differ.
renames <- c(
  "@variables" = "variables", "@parameters" = "parameters", "@constants" = "constants",
  "@brownians" = "brownians", "@brownian" = "brownians", "@poissonians" = "Poissonian",
  "@discretes" = "discretes", "@independent_variables" = "independent_variable",
  "@named" = "named", "@connector" = "connector", "@component" = "component",
  "@mtkcompile" = "mtkcompile", "@mtkbuild" = "mtkcompile", "@mtkcomplete" = "complete",
  "@derivatives" = "Differential", "@nonamespace" = "nonamespace", "@namespace" = "namespace",
  "t" = "t_nounits", "D" = "D_nounits", "substitute" = "substitute_vars",
  "simplify" = "simplify", "solve" = "solve.ODEProblem", "independent_variable" = "independent_variable",
  "find_solvables!" = "find_solvables", "modified_unknowns!" = "modified_unknowns",
  "reorder_dimension_by_tunables!" = "reorder_dimension_by_tunables",
  "collect_scoped_vars!" = "collect_scoped_vars", "collect_var_to_name!" = "collect_var_to_name",
  "collect_vars!" = "collect_vars", "store_to_mutable_cache!" = "store_to_mutable_cache")

s3 <- sub("^S3method[(]([^,]+),([^)]+)[)]$", "\\1.\\2", grep("^S3method[(]", ns, value = TRUE))
symcomp <- getNamespaceExports("symcompR")
status <- vapply(ledger$name, function(n) {
  if (n %in% na$name) return("not applicable")
  if (n %in% stubs$name) return("stub")
  r <- if (!is.na(renames[n])) renames[[n]] else sub("!$", "", sub("^@", "", n))
  if (r %in% ex || r %in% s3) return(if (identical(r, n)) "ported" else paste0("ported as ", r))
  if (r %in% symcomp) return(paste0("from symcompR (", r, ")"))
  "missing"
}, "")
out <- data.frame(name = ledger$name, source = ledger$source, status = status)
utils::write.table(out, "tools/ledger_status.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
done <- !grepl("^missing", status)
cat(sprintf("%d of %d names accounted for: %d ported, %d stubs that redirect, %d not applicable\n",
            sum(done), length(done), sum(grepl("^ported|symcompR", status)),
            sum(status == "stub"), sum(status == "not applicable")))
if (any(!done)) cat("Missing:", paste(out$name[!done], collapse = " "), "\n")

# The SciML packages the port also covers (tools/upstream_exports.tsv, from
# tools/upstream_exports.jl, and the solver sublibraries in
# tools/upstream_exports_src.tsv, from tools/upstream_exports_src.jl). A name
# counts once, under the first package that exports it;
# tools/upstream_not_applicable.tsv lists the Julia internals with no R
# meaning, each with the R feature that covers it, and
# tools/upstream_aliases.tsv the names that select a Julia library or device
# and run another algorithm of the port (reported as aliases).
if (file.exists("tools/upstream_exports.tsv")) {
  up <- utils::read.delim("tools/upstream_exports.tsv", stringsAsFactors = FALSE)
  if (file.exists("tools/upstream_exports_src.tsv")) {
    src <- utils::read.delim("tools/upstream_exports_src.tsv", stringsAsFactors = FALSE)
    up <- rbind(up, src[!startsWith(src$name, "@reexport "), ])
  }
  up <- up[!duplicated(up$name) & !(up$name %in% ledger$name), ]
  ual <- if (file.exists("tools/upstream_aliases.tsv"))
    utils::read.delim("tools/upstream_aliases.tsv", stringsAsFactors = FALSE) else
    data.frame(name = character(), runs = character())
  una <- if (file.exists("tools/upstream_not_applicable.tsv"))
    utils::read.delim("tools/upstream_not_applicable.tsv", stringsAsFactors = FALSE) else
    data.frame(name = character(), reason = character())
  # Julia names that are R functions of other packages, or need a suffix
  up_renames <- c("solve!" = "solve_integrator", "step!" = "step_integrator",
                  "initialize_dae!" = "initialize_dae", "postamble!" = "postamble", "set_ut!" = "set_ut",
                  "save_discretes!" = "save_discretes", "save_discretes_if_enabled!" = "save_discretes_if_enabled",
                  "save_final_discretes!" = "save_final_discretes", "pop_tstop!" = "pop_tstop",
                  "set_rng!" = "set_rng", "auto_dt_reset!" = "auto_dt_reset",
                  "change_t_via_interpolation!" = "change_t_via_interpolation",
                  "reeval_internals_due_to_modification!" = "reeval_internals_due_to_modification",
                  "addat!" = "addat", "deleteat!" = "deleteat", "resize_non_user_cache!" = "resize_non_user_cache",
                  "addat_non_user_cache!" = "addat_non_user_cache", "deleteat_non_user_cache!" = "deleteat_non_user_cache")
  ustatus <- vapply(up$name, function(n) {
    if (n %in% una$name) return("not applicable")
    r <- if (!is.na(up_renames[n])) up_renames[[n]] else sub("!$", "", sub("^@", "", n))
    if (r %in% ex || r %in% s3) {
      if (n %in% ual$name) return("alias")
      return(if (identical(r, n)) "ported" else paste0("ported as ", r))
    }
    if (r %in% symcomp) return(paste0("from symcompR (", r, ")"))
    "missing"
  }, "")
  uout <- data.frame(package = up$package, name = up$name, status = ustatus)
  utils::write.table(uout, "tools/upstream_status.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  for (pk in unique(uout$package)) {
    st <- uout$status[uout$package == pk]
    cat(sprintf("%-36s %4d names: %4d ported, %3d aliases, %3d not applicable, %4d missing\n", pk, length(st),
                sum(grepl("^ported|symcompR", st)), sum(st == "alias"), sum(st == "not applicable"),
                sum(st == "missing")))
  }
  st <- uout$status
  cat(sprintf("%-36s %4d names: %4d ported, %3d aliases, %3d not applicable, %4d missing\n", "Total", length(st),
              sum(grepl("^ported|symcompR", st)), sum(st == "alias"), sum(st == "not applicable"),
              sum(st == "missing")))
}
