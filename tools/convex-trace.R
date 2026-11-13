# Per iteration comparison with the trace that tools/convex-probe.jl writes:
#   TRACE=convex-trace.tsv Rscript tools/convex-trace.R [case ...]
# Iterate k is the result of a solve with max_iter = k, as on the Julia side.
# For each case it prints the scaled problem data that differ and the first
# iterate whose quantities differ, with each quantity's relative difference.
# Without TRACE it prints the R trace in the same format.
suppressMessages(devtools::load_all(".", quiet = TRUE))
ex <- parse("tests/testthat/test-convex.R")
for (e in ex) if (!(is.call(e) && identical(e[[1]], as.name("test_that")))) eval(e, globalenv())
cases <- commandArgs(TRUE)
if (!length(cases)) cases <- c("pow", "expcon", "geomean", "relentr", "con_explog", "eig_nested", "logdet_con",
                               "logdet", "rootdet_tri_nocd", "rootdet_sq_nocd")
qnames <- c("mu", "sigma", "step_length", "cost_primal", "cost_dual", "res_primal", "res_dual", "gap_abs",
            "gap_rel", "ktratio", "res_primal_inf", "res_dual_inf")
hx <- function(v) paste(vapply(v, function(x) if (is.nan(x)) "NaN" else if (is.infinite(x)) (if (x > 0) "Inf" else "-Inf")
                                   else sprintf("%a", x), ""), collapse = ",")
num <- function(s) if (is.na(s) || !nzchar(s)) numeric() else as.numeric(strsplit(s, ",", fixed = TRUE)[[1]])
rel <- function(a, b) if (length(a) != length(b)) Inf else if (!length(a)) 0 else
  max(ifelse(a == b | (is.nan(a) & is.nan(b)), 0, abs(a - b) / pmax(abs(b), 1e-300)))

with_max_iter <- function(alg, k) {
  oc <- alg$optimizer_constructor
  st <- if (is.function(oc)) oc()$settings else oc$settings
  st$max_iter <- as.integer(k)
  ConvexMOI(structure(list(settings = st), class = "ClarabelOptimizer"))
}
r_iterate <- function(prob, alg, k) {
  cache <- init(prob, with_max_iter(alg, k))
  suppressWarnings(solve_integrator(cache))
  S <- cache$model$solver
  i <- S$info
  list(status = i$status, q = vapply(qnames, function(n) as.numeric(get(n, envir = i)), 0),
       x = S$solution$x, z = S$solution$z, s = S$solution$s,
       data = list(q = S$data$q, b = S$data$b, A = S$data$A$nzval, P = S$data$P$nzval, d = S$data$d, e = S$data$e,
                   c = S$data$c))
}

tracepath <- Sys.getenv("TRACE")
jl <- if (nzchar(tracepath)) strsplit(readLines(tracepath), "\t", fixed = TRUE) else list()
pr <- c(cvx_problems(), cvx_problems_014())
for (nm in cases) {
  base <- sub("_nocd$", "", nm)
  prob <- pr[[base]]
  if (is.null(prob)) { cat(nm, ": no R problem\n"); next }
  alg <- cvx_alg_for(nm)
  rows <- Filter(function(r) r[[1]] == nm, jl)
  if (!length(rows)) {
    for (k in 0:cvx_run(prob, alg)$iters) {
      it <- r_iterate(prob, alg, k)
      cat(nm, k, it$status, hx(it$q), hx(it$x), hx(it$z), hx(it$s), sep = "\t")
      cat("\n")
    }
    next
  }
  first <- TRUE
  for (r in rows) {
    if (r[[2]] == "DATA") {
      it <- r_iterate(prob, alg, 0)
      dd <- mapply(function(a, b) rel(a, num(b)), it$data, r[3:9])
      if (any(dd > 0)) cat(nm, "DATA differs:", paste(names(dd)[dd > 0], signif(dd[dd > 0], 3), collapse = " "), "\n")
      next
    }
    k <- as.integer(r[[2]])
    if (r[[3]] == "ERROR") { cat(nm, k, "julia ERROR", r[[4]], "\n"); next }
    it <- r_iterate(prob, alg, k)
    dq <- mapply(rel, it$q, num(r[[4]]))
    dv <- c(x = rel(it$x, num(r[[5]])), z = rel(it$z, num(r[[6]])), s = rel(it$s, num(r[[7]])))
    d <- c(dq, dv)
    if (!identical(it$status, r[[3]]) || any(d > 0)) {
      if (first) {
        cat(nm, "first difference at k =", k, "status", it$status, "/", r[[3]], "\n")
        o <- d[d > 0]
        cat("  ", paste(names(o), signif(o, 3), collapse = " "), "\n")
        cat("  R    ", hx(it$q), "\n  julia", r[[4]], "\n")
        first <- FALSE
      } else {
        cat(nm, "k =", k, "max rel", signif(max(d), 3), "status", it$status, "/", r[[3]], "\n")
      }
    }
  }
  if (first) cat(nm, ": all", length(rows) - 1, "iterates identical\n")
}
