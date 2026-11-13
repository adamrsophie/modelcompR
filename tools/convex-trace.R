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
if (!length(cases)) cases <- c("pow", "expcon", "geomean", "dualgeomean", "relentr", "con_explog", "eig_nested",
                               "logdet_con", "logdet", "rootdet_tri_nocd", "rootdet_sq_nocd")
qnames <- c("mu", "sigma", "step_length", "cost_primal", "cost_dual", "res_primal", "res_dual", "gap_abs",
            "gap_rel", "ktratio", "res_primal_inf", "res_dual_inf")
hx <- function(v) paste(vapply(v, function(x) if (is.nan(x)) "NaN" else if (is.infinite(x)) (if (x > 0) "Inf" else "-Inf")
                                   else sprintf("%a", x), ""), collapse = ",")
num <- function(s) if (is.na(s) || !nzchar(s)) numeric() else as.numeric(strsplit(s, ",", fixed = TRUE)[[1]])
rel <- function(a, b) if (length(a) != length(b)) Inf else if (!length(a)) 0 else
  max(ifelse(a == b | (is.nan(a) & is.nan(b)), 0, abs(a - b) / pmax(abs(b), 1e-300)))

# Hooks that record the affine step length, the affine step and the
# Mehrotra shift of the last iteration, as the AFF lines of the probe
hook <- new.env()
ns <- asNamespace("modelcompR")
orig_shift <- ns$.cl_combined_ds_shift
orig_len <- ns$.cl_get_step_length
unlockBinding(".cl_combined_ds_shift", ns)
unlockBinding(".cl_get_step_length", ns)
assign(".cl_combined_ds_shift", function(cones, step_z, step_s, sigmamu) {
  out <- orig_shift(cones, step_z, step_s, sigmamu)
  hook$aff <- list(alpha = hook$alpha, sigmamu = sigmamu, step_z = step_z, step_s = step_s, shift = out)
  out
}, ns)
assign(".cl_get_step_length", function(S, steptype, scaling) {
  a <- orig_len(S, steptype, scaling)
  if (steptype == "affine") hook$alpha <- a
  a
}, ns)

# each KKT solve, as the KKTA (affine) and KKTC (combined) lines of the
# probe: right hand side, static regularizer, KKT values, D^-1, the solution
# before refinement, the refinement error norms and candidates, the result.
# The body repeats .cl_kktsolver_solve.
unlockBinding(".cl_kktsolver_solve", ns)
assign(".cl_kktsolver_solve", function(ks) {
  b <- ks$b
  rec <- list(b = b, eps = ks$diagonal_regularizer %||% 0, KKT = ks$KKT$nzval, Dinv = ks$ldl$Dinv)
  x <- ns$.qdldl_solve(ks$ldl, b)
  rec$x0 <- x
  norms <- numeric(); cands <- numeric()
  st <- ks$settings
  ok <- TRUE
  if (st$iterative_refinement_enable) {
    refine_err <- function(xi) ns$.cl_symv(ks$KKT, xi, b, -1)
    normb <- ns$.cl_norminf(b)
    e <- refine_err(x)
    norme <- ns$.cl_norminf(e); norms <- c(norms, norme)
    if (!is.finite(norme)) ok <- FALSE else {
      for (i in seq_len(st$iterative_refinement_max_iter)) {
        if (norme <= st$iterative_refinement_abstol + st$iterative_refinement_reltol * normb) break
        lastnorme <- norme
        dx <- ns$.qdldl_solve(ks$ldl, e)
        dx <- dx + x
        cands <- c(cands, dx)
        e <- refine_err(dx)
        norme <- ns$.cl_norminf(e); norms <- c(norms, norme)
        if (!is.finite(norme)) {
          ok <- FALSE
          break
        }
        ratio <- lastnorme / norme
        if (ratio < st$iterative_refinement_stop_ratio) {
          if (ratio > 1) x <- dx
          break
        }
        x <- dx
      }
    }
  } else ok <- all(is.finite(x))
  ks$x <- x
  rec$norms <- norms; rec$cands <- cands; rec$xout <- x
  hook$kkt <- c(hook$kkt, list(rec))
  ok
}, ns)

with_max_iter <- function(alg, k) {
  oc <- alg$optimizer_constructor
  st <- if (is.function(oc)) oc()$settings else oc$settings
  st$max_iter <- as.integer(k)
  ConvexMOI(structure(list(settings = st), class = "ClarabelOptimizer"))
}
r_iterate <- function(prob, alg, k) {
  hook$aff <- NULL
  hook$kkt <- NULL
  cache <- init(prob, with_max_iter(alg, k))
  suppressWarnings(solve_integrator(cache))
  S <- cache$model$solver
  i <- S$info
  list(status = i$status, q = vapply(qnames, function(n) as.numeric(get(n, envir = i)), 0),
       x = S$solution$x, z = S$solution$z, s = S$solution$s,
       data = list(q = S$data$q, b = S$data$b, A = S$data$A$nzval, P = S$data$P$nzval, d = S$data$d, e = S$data$e,
                   c = S$data$c),
       step = with(S$step_lhs, list(x = x, z = z, s = s, tau = tau, kappa = kappa)),
       cones = lapply(S$cones$cones, function(K) if (is.null(K$H_dual)) NULL else
         list(H_dual = as.numeric(K$H_dual), Hs = as.numeric(K$Hs), grad = K$grad, z = K$z)),
       state = lapply(S$cones$cones, function(K) switch(K$type,
         SecondOrder = list(w = K$w, lambda = K$lambda, eta = K$eta),
         PSDTriangle = list(lambda = K$lambda, R = as.numeric(K$R), Rinv = as.numeric(K$Rinv)), NULL)),
       aff = hook$aff, kkt = tail(hook$kkt, 2L))
}
deep <- function(nm, tag, got, want) {
  d <- mapply(function(a, b) rel(a, num(b)), got, want)
  cat(nm, tag, if (all(d == 0)) "identical" else paste(names(d), signif(d, 3), collapse = " "), "\n")
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
    if (startsWith(r[[2]], "STEP")) {
      deep(nm, r[[2]], r_iterate(prob, alg, as.integer(sub("STEP", "", r[[2]])))$step, r[3:7])
      next
    }
    if (startsWith(r[[2]], "CONE")) {
      kj <- as.integer(strsplit(sub("CONE", "", r[[2]]), ":", fixed = TRUE)[[1]])
      deep(nm, r[[2]], r_iterate(prob, alg, kj[1])$cones[[kj[2]]], r[3:6])
      next
    }
    if (grepl("^(SOC|PSD)", r[[2]])) {
      kj <- as.integer(strsplit(sub("^(SOC|PSD)", "", r[[2]]), ":", fixed = TRUE)[[1]])
      st <- r_iterate(prob, alg, kj[1])$state[[kj[2]]]
      deep(nm, r[[2]], st, r[2 + seq_along(st)])
      next
    }
    if (grepl("^KKT[AC]", r[[2]])) {
      kk <- r_iterate(prob, alg, as.integer(sub("^KKT[AC]", "", r[[2]])))$kkt
      rec <- kk[[if (startsWith(r[[2]], "KKTA")) 1L else 2L]]
      want <- r[3:10]
      want[7] <- gsub(";", ",", want[7], fixed = TRUE)
      deep(nm, r[[2]], rec, want)
      next
    }
    if (startsWith(r[[2]], "AFF")) {
      deep(nm, r[[2]], r_iterate(prob, alg, as.integer(sub("AFF", "", r[[2]])))$aff, r[3:7])
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
        if (nzchar(Sys.getenv("VERBOSE"))) {
          cat("  ", paste(names(d), signif(d, 3), collapse = " "), "\n")
          cat("  R    ", hx(it$q), "\n  julia", r[[4]], "\n")
        }
      }
    }
  }
  if (first) cat(nm, ": all", length(rows) - 1, "iterates identical\n")
}
