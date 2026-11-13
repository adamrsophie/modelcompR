# Benchmarks of DiffEqDevTools.jl: work precision diagrams, shootouts, the
# tags of benchmark entries and their selection (src/benchmark.jl), the
# comparison plot sets (src/autoplot.jl), the WorkPrecision,
# WorkPrecisionSet and Shootout plot recipes (src/plotrecipes.jl) and the
# preset algorithm tags of ext/DiffEqDevToolsOrdinaryDiffEqCoreExt.jl.
#
# Setups are lists with the algorithm `alg` and optionally `abstols`,
# `reltols`, `dts`, `name`, `tags`, `auto_tags`, `prob_choice` and solver
# options. Errors of a WorkPrecision are a named list of vectors with one
# value per tolerance (the StructArray of Julia), with the error names of
# R/devtools-convergence.R. Tags are character vectors (Julia symbols).
# Times are wall clock seconds, so they depend on the machine.

.dv_timeseries_keys <- c("l2", "l_inf", "L2", "L_inf")
.dv_dense_keys <- c("L2", "L_inf")
.dv_weak_keys <- c("weak_final", "weak_l2", "weak_l_inf", "weak_L2", "weak_L_inf")

# SciMLBase.allowedkeywords (src/errors.jl): the setup entries passed on to
# solve
.dv_allowedkeywords <- c(
  "dense", "saveat", "save_idxs", "save_discretes", "tstops", "tspan", "d_discontinuities", "save_everystep",
  "save_on", "save_start", "save_end", "initialize_save", "adaptive", "abstol", "reltol", "dt", "dtmax", "dtmin",
  "force_dtmin", "internalnorm", "controller", "failfactor", "calck", "alias_u0", "maxiters", "maxtime",
  "callback", "isoutofdomain", "domain_checks", "unstable_check", "verbose", "merge_callbacks", "progress",
  "progress_steps", "progress_name", "progress_message", "progress_id", "timeseries_errors", "dense_errors",
  "weak_timeseries_errors", "weak_dense_errors", "wrap", "calculate_error", "initializealg", "alg", "save_noise",
  "delta", "seed", "rng", "alg_hints", "kwargshandle", "trajectories", "batch_size", "sensealg",
  "advance_to_tstop", "stop_at_next_tstop", "u0", "p", "default_set", "second_time", "prob_choice", "alias_jump",
  "alias_noise", "batch", "nlsolve_kwargs", "odesolve_kwargs", "linsolve_kwargs", "ensemblealg", "rng_func",
  "show_trace", "trace_level", "store_trace", "termination_condition", "precondition", "postcondition", "alias",
  "tune_parameters", "optimize_kwargs", "step_limiter", "stage_limiter")

# Entries of b replace those of a (NULL values included).
.dv_merge <- function(a, b) {
  for (n in names(b)) a[n] <- list(b[[n]])
  a
}
.dv_now <- function() as.numeric(Sys.time())
.dv_nan_errors <- function() list(l_inf = NaN, L2 = NaN, final = NaN, l2 = NaN, L_inf = NaN)

# The per tolerance error lists as one list of vectors over the union of the
# names (sorted), NaN where a tolerance has no such error.
.dv_dicts_to_soa <- function(dicts) {
  keys <- sort(unique(unlist(lapply(dicts, names))), method = "radix")
  stats::setNames(lapply(keys, function(k) vapply(dicts, function(d) as.numeric(d[[k]] %||% NaN), 0)), keys)
}

# _default_name: the type name of the algorithm (that of the inner
# algorithm for MethodOfSteps)
.dv_default_name <- function(alg) {
  if (inherits(alg, "DDEAlgorithm") && !is.null(alg$alg)) return(.dv_default_name(alg$alg))
  if (inherits(alg, "CompositeAlgorithm")) return("CompositeAlgorithm")
  class(alg)[1L]
}
.dv_setup_name <- function(setup) setup$name %||% .dv_default_name(setup$alg)
.dv_setup_tags <- function(setup) as.character(unlist(setup$tags))
.dv_setup_options <- function(setup) {
  s <- setup[names(setup) %in% .dv_allowedkeywords]
  s$alg <- NULL
  s
}
# The problem and reference of a setup with prob_choice, or the first of a
# list of problems.
.dv_is_prob_list <- function(prob) is.list(prob) && !inherits(prob, "SciMLProblem")
.dv_choose <- function(x, prob, pc) {
  if (is.null(x)) return(NULL)
  if (!is.null(pc)) return(x[[pc]])
  if (.dv_is_prob_list(prob) && is.list(x) && !inherits(x, c("SciMLSolution", "TestSolution"))) return(x[[1L]])
  x
}

.dv_timed_out <- function(timeout, solve_time, name, abstol) {
  if (is.null(timeout) || solve_time <= timeout) return(FALSE)
  warning(name, " exceeded the ", timeout, "s timeout at abstol=", abstol, " (", signif(solve_time, 3),
          "s); recording this point as NaN.", call. = FALSE)
  TRUE
}

# The fastest of numruns timings (one when it takes longer than seconds),
# after a first run.
.dv_best_time <- function(run, numruns, seconds) {
  timed <- function() {
    t0 <- .dv_now()
    run()
    .dv_now() - t0
  }
  timed()
  b_t <- timed()
  if (b_t > seconds || numruns < 2L) return(b_t)
  min(b_t, vapply(seq_len(numruns - 1L), function(i) timed(), 0))
}

# ---------------------------------------------------------------------------
# Tags

.dv_family_tags <- c("adams", "bdf", "esdirk", "exponential", "explicit_rk", "extrapolation", "firk", "imex",
                     "low_storage_rk", "multistep", "partitioned", "rk", "rkn", "rosenbrock", "sdirk", "ssprk",
                     "stabilized", "stabilized_rk", "symplectic")
.dv_trait_tags <- c("adaptive", "composite", "explicit", "fixed_step", "implicit", "split", "variable_order")
.dv_role_tags <- c("recommended", "reference")
.dv_provider_tags <- c("external", "lsoda", "odeinterface", "ordinarydiffeq", "sundials")
.dv_variant_tags <- c("dae_residual", "dense", "dense_output", "dynamic", "endpoint", "mass_matrix",
                      "matrix_free", "sparse", "static", "timeseries")
.dv_domain_tags <- c("bv", "dae", "dde", "nonstiff", "ode", "pde", "rode", "sde", "stiff")

#' Tags of benchmark entries
#'
#' Entries of a [WorkPrecisionSet()] carry tags (character strings): those
#' of their setup (`tags`) and, unless the setup sets `auto_tags = FALSE`,
#' the tags `auto_tags(alg)` derived from the algorithm: its order
#' (`order_5`; not for variable order methods), `adaptive` or `fixed_step`,
#' and the preset tags of OrdinaryDiffEq algorithms (`explicit` or
#' `implicit`, `dae`, `sde`, `rode`, `rosenbrock`, `exponential`,
#' `partitioned`, `composite`, `adams`, `variable_order`, `multistep`,
#' `firk`, `sdirk`, `esdirk`, `split`).
#'
#' `tag_kind(tag)` classifies tags as `"family"`, `"trait"`, `"role"`,
#' `"provider"`, `"variant"`, `"domain"` or `"unknown"`. `get_tags(set)`
#' lists the tags of each entry and `unique_tags(set)` all tags, sorted.
#' `filter_by_tags(set, ...)` keeps the entries carrying all the given tags,
#' `exclude_by_tags(set, ...)` drops those carrying any of them, and
#' `merge_wp_sets(...)` concatenates sets (keeping the tolerances, problem
#' and error estimate of the first).
#'
#' @param alg An algorithm.
#' @param tag Tags.
#' @param wp_set A `WorkPrecisionSet`.
#' @param ... Tags, or for `merge_wp_sets()` the sets.
#' @return `auto_tags()`, `tag_kind()` and `unique_tags()` return character
#'   vectors, `get_tags()` a list of them, and the others a
#'   `WorkPrecisionSet`.
#' @examples
#' auto_tags(Tsit5())
#' auto_tags(Rodas5P())
#' tag_kind(c("rosenbrock", "adaptive", "reference", "foo"))
#' @export
auto_tags <- function(alg) {
  if (!inherits(alg, c("ODEAlgorithm", "SDEAlgorithm", "RODEAlgorithm", "DDEAlgorithm", "BVPAlgorithm")))
    return(character())
  presets <- .dv_preset_tags(alg)
  tags <- character()
  if (!"variable_order" %in% presets) tags <- c(tags, .dv_order_tag(alg))
  inner <- if (inherits(alg, "DDEAlgorithm") && !is.null(alg$alg)) alg$alg else alg
  # SciMLBase defaults isadaptive to true; the boundary value algorithms of
  # R (which refine their mesh) carry no adaptive field
  adaptive <- if (inherits(inner, "BVPAlgorithm") && is.null(inner$adaptive)) TRUE else isTRUE(isadaptive(inner))
  tags <- c(tags, if (adaptive) "adaptive" else "fixed_step", presets)
  unique(tags)
}

#' @rdname auto_tags
#' @export
tag_kind <- function(tag) {
  vapply(as.character(tag), function(t) {
    if (t %in% .dv_family_tags) return("family")
    if (t %in% .dv_trait_tags || startsWith(t, "order_")) return("trait")
    if (t %in% .dv_role_tags) return("role")
    if (t %in% .dv_provider_tags) return("provider")
    if (t %in% .dv_variant_tags || startsWith(t, "autodiff_")) return("variant")
    if (t %in% .dv_domain_tags) return("domain")
    "unknown"
  }, "", USE.NAMES = FALSE)
}

# order_<n>: the Julia rendering of alg_order (rationals such as 1//2 for
# stochastic methods give order_1_2). Where the order stored by an R
# stochastic algorithm is its weak order, the table gives the Julia alg_order.
.dv_julia_orders <- c(
  COM = "1_1", DRI1 = "1_1", DRI1NM = "1_1", IRI1 = "1_1", NON = "1_1", NON2 = "1_1", PL1WM = "1_1",
  PL1WMA = "1_1", RDI2WM = "1_1", RDI3WM = "1_1", RDI4WM = "1_1", RI1 = "1_1", RI3 = "1_1", RI5 = "1_1",
  RI6 = "1_1", RS1 = "1_1", RS2 = "1_1", SIEA = "1_1", SIEB = "1_1", SMEA = "1_1", SMEB = "1_1",
  W2Ito1 = "1_1", SOSRA = "2_1", SOSRA2 = "2_1", SRA1 = "2_1", SRA2 = "2_1", SRA3 = "2_1")
.dv_order_tag <- function(alg) {
  inner <- if (inherits(alg, "DDEAlgorithm") && !is.null(alg$alg)) alg$alg else alg
  known <- .dv_julia_orders[class(inner)[1L]]
  if (!is.na(known)) return(paste0("order_", known))
  o <- tryCatch(alg_order(inner), error = function(e) NA_real_)
  if (length(o) != 1L || is.na(o)) return(character())
  if (inherits(inner, c("SDEAlgorithm", "RODEAlgorithm"))) {
    q <- which(abs(o * seq_len(12L) - round(o * seq_len(12L))) < 1e-12)[1L]
    if (!is.na(q)) return(paste0("order_", round(o * q), "_", q))
  }
  s <- if (o == round(o)) format(as.integer(o)) else format(o, digits = 15)
  paste0("order_", gsub(".", "_", s, fixed = TRUE))
}

# _preset_algorithm_tags of the OrdinaryDiffEqCore extension: the domain
# (dae and implicit, implicit or explicit, sde, rode) then the family of the
# method (Julia tests the abstract type; R names the methods of each family).
.dv_family_presets <- c(Rodas4 = "rosenbrock", Rodas5P = "rosenbrock", Rosenbrock23 = "rosenbrock",
                        SplitEM = "split", SplitEuler = "exponential")
.dv_preset_tags <- function(alg) {
  if (inherits(alg, "CompositeAlgorithm")) {
    return(unique(c(if (any(vapply(alg$algs, function(a) isTRUE(a$stiff), TRUE))) "implicit" else "explicit",
                    "composite")))
  }
  if (!inherits(alg, c("ODEAlgorithm", "SDEAlgorithm", "RODEAlgorithm"))) return(character())
  domain <- if (isTRUE(alg$dae)) c("dae", "implicit") else if (inherits(alg, "SDEAlgorithm")) "sde" else
    if (inherits(alg, "RODEAlgorithm")) "rode" else if (isTRUE(alg$stiff)) "implicit" else "explicit"
  fam <- .dv_family_presets[class(alg)[1L]]
  c(domain, if (!is.na(fam)) unname(fam))
}

#' @rdname auto_tags
#' @export
get_tags <- function(wp_set) lapply(wp_set$wps, function(wp) wp$tags)

#' @rdname auto_tags
#' @export
unique_tags <- function(wp_set) unique(sort(unlist(get_tags(wp_set)) %||% character(), method = "radix"))

.dv_has_all <- function(wp, tags) all(tags %in% wp$tags)
.dv_has_any <- function(wp, tags) any(tags %in% wp$tags)

#' @rdname auto_tags
#' @export
filter_by_tags <- function(wp_set, ...) {
  tags <- as.character(unlist(list(...)))
  if (!length(tags)) return(wp_set)
  .dv_subset_wps(wp_set, which(vapply(wp_set$wps, .dv_has_all, TRUE, tags = tags)))
}

#' @rdname auto_tags
#' @export
exclude_by_tags <- function(wp_set, ...) {
  tags <- as.character(unlist(list(...)))
  if (!length(tags)) return(wp_set)
  .dv_subset_wps(wp_set, which(!vapply(wp_set$wps, .dv_has_any, TRUE, tags = tags)))
}

#' @rdname auto_tags
#' @export
merge_wp_sets <- function(...) {
  sets <- list(...)
  if (!length(sets)) stop("At least one WorkPrecisionSet is required", call. = FALSE)
  wps <- do.call(c, lapply(sets, function(s) s$wps))
  first <- sets[[1L]]
  .dv_wps_new(wps, length(wps), first$abstols, first$reltols, first$prob,
              do.call(c, lapply(sets, function(s) s$setups)), do.call(c, lapply(sets, function(s) s$names)),
              first$error_estimate, first$numruns, first$active_error_estimates)
}

.dv_subset_wps <- function(wp_set, indices) {
  if (!length(indices)) warning("No entries match the requested tags. Returning an empty WorkPrecisionSet.",
                                call. = FALSE)
  .dv_wps_new(wp_set$wps[indices], length(indices), wp_set$abstols, wp_set$reltols, wp_set$prob,
              if (is.list(wp_set$setups)) wp_set$setups[indices] else wp_set$setups,
              if (length(wp_set$names) > 1L || is.list(wp_set$names)) wp_set$names[indices] else wp_set$names,
              wp_set$error_estimate, wp_set$numruns, wp_set$active_error_estimates)
}

# The entries selected by the tag options of the plot: tags narrows (all of
# them), include_tags adds back the entries with all of those, exclude_tags
# drops the entries with any of those.
.dv_selected_indices <- function(wp_set, tags = NULL, include_tags = NULL, exclude_tags = NULL) {
  idx <- seq_along(wp_set$wps)
  if (!is.null(tags)) idx <- idx[vapply(wp_set$wps[idx], .dv_has_all, TRUE, tags = as.character(tags))]
  if (!is.null(include_tags)) {
    extra <- which(vapply(wp_set$wps, .dv_has_all, TRUE, tags = as.character(include_tags)))
    idx <- sort(union(idx, extra))
  }
  if (!is.null(exclude_tags)) idx <- idx[!vapply(wp_set$wps[idx], .dv_has_any, TRUE, tags = as.character(exclude_tags))]
  idx
}

# ---------------------------------------------------------------------------
# Best of family selection

# (log10 error, log10 time) of the usable tolerances, by increasing error
.dv_wp_points <- function(wp) {
  e <- wp$errors[[wp$error_estimate]]
  tt <- wp$times
  ok <- !is.na(e) & !is.na(tt) & e > 0 & tt > 0
  pts <- cbind(log10(e[ok]), log10(tt[ok]))
  pts[order(pts[, 1L], method = "radix"), , drop = FALSE]
}
.dv_trapezoid <- function(pts) {
  n <- nrow(pts)
  sum(0.5 * (pts[-1L, 2L] + pts[-n, 2L]) * (pts[-1L, 1L] - pts[-n, 1L]))
}

#' Rank benchmark entries
#'
#' `wp_area(wp)` is the area under the work precision curve of a
#' [WorkPrecision()] in log10 coordinates (trapezoid rule over the
#' tolerances with a positive error and time, by increasing error), or `Inf`
#' with fewer than two such points; lower is better.
#'
#' `best_by_tag(set, tag, n)` keeps the `n` best entries carrying `tag`:
#' first those with the most usable points, then the lowest mean log10 time
#' over the error range they cover (the area divided by that range).
#' `best_of_families(set, family_tags, n)` keeps the `n` best of each
#' family, and signals an error when no family has an entry.
#'
#' @param wp A `WorkPrecision`.
#' @param wp_set A `WorkPrecisionSet`.
#' @param tag,family_tags Tags.
#' @param n Number of entries per tag.
#' @param metric The ranking; only `"area"`.
#' @return A number, or a `WorkPrecisionSet`.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
#' tols <- c(1e-3, 1e-4, 1e-5)
#' mk <- function(name, e, t, tags) WorkPrecision(prob, tols, tols, list(final = e), t, NULL, NULL,
#'   name, "final", 3, tags)
#' set <- WorkPrecisionSet(list(mk("A", c(1e-2, 1e-4, 1e-6), c(1, 2, 5) / 1000, "rk"),
#'   mk("B", c(1e-3, 1e-5, 1e-7), c(2, 4, 9) / 1000, "rk")), 2, tols, tols, prob, NULL,
#'   c("A", "B"), "final", NULL)
#' wp_area(set[[1]])
#' best_by_tag(set, "rk")$names
#' @export
wp_area <- function(wp) {
  pts <- .dv_wp_points(wp)
  if (nrow(pts) < 2L) return(Inf)
  .dv_trapezoid(pts)
}

.dv_wp_rank <- function(wp) {
  pts <- .dv_wp_points(wp)
  n <- nrow(pts)
  if (n < 2L) return(c(-n, Inf))
  span <- pts[n, 1L] - pts[1L, 1L]
  if (span == 0) return(c(-n, sum(pts[, 2L]) / n))
  c(-n, .dv_trapezoid(pts) / span)
}

.dv_best_indices <- function(wp_set, tag, n, metric) {
  if (!identical(metric, "area")) stop("Unknown metric: ", metric, ". Supported metrics: area", call. = FALSE)
  cand <- which(vapply(wp_set$wps, function(wp) tag %in% wp$tags, TRUE))
  if (!length(cand)) return(cand)
  r <- do.call(rbind, lapply(wp_set$wps[cand], .dv_wp_rank))
  perm <- order(r[, 1L], r[, 2L], method = "radix")
  cand[perm[seq_len(min(n, length(perm)))]]
}

.dv_best_family_indices <- function(wp_set, family_tags, n, metric)
  sort(unique(unlist(lapply(family_tags, function(tag) .dv_best_indices(wp_set, tag, n, metric)))))

#' @rdname wp_area
#' @export
best_by_tag <- function(wp_set, tag, n = 1L, metric = "area")
  .dv_subset_wps(wp_set, .dv_best_indices(wp_set, as.character(tag), n, metric))

#' @rdname wp_area
#' @export
best_of_families <- function(wp_set, family_tags, n = 1L, metric = "area") {
  family_tags <- as.character(family_tags)
  idx <- .dv_best_family_indices(wp_set, family_tags, n, metric)
  if (!length(idx))
    stop("No entries found for any of the family tags ", paste(family_tags, collapse = ", "), call. = FALSE)
  .dv_subset_wps(wp_set, idx)
}

# ---------------------------------------------------------------------------
# Comparison plot sets

#' Comparison sets of a work precision benchmark
#'
#' `autoplot(set)` splits a [WorkPrecisionSet()] into the sets one would
#' plot to compare its methods, as a named list: `family_<tag>` with the
#' entries of each family present (by default the tags of kind `"family"`
#' of the set), plus the entries tagged with any of `reference_tags`;
#' `best_of_families` with the `best_n` best entries of each family (and
#' the references); and `all`, the whole set. It draws nothing; plot the
#' sets with [plot()].
#'
#' The same method is registered for the `autoplot()` generic of ggplot2,
#' so it also answers `ggplot2::autoplot(set)`, and `autoplot()` passes
#' other objects on to ggplot2 when it is installed.
#'
#' @param object A `WorkPrecisionSet`.
#' @param families Family tags.
#' @param reference_tags Tags of reference entries added to every set.
#' @param best_n Number of entries per family in `best_of_families`.
#' @param ... Unused.
#' @return A named list of `WorkPrecisionSet`s.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
#' tols <- c(1e-3, 1e-4)
#' mk <- function(name, e, tags) WorkPrecision(prob, tols, tols, list(final = e), c(1, 2) / 1000, NULL,
#'   NULL, name, "final", 2, tags)
#' set <- WorkPrecisionSet(list(mk("A", c(1e-2, 1e-4), "rk"), mk("B", c(1e-3, 1e-5), "rosenbrock")), 2,
#'   tols, tols, prob, NULL, c("A", "B"), "final", NULL)
#' names(autoplot(set))
#' @export
autoplot <- function(object, ...) {
  # other objects go to ggplot2; a method autoplot.default here would be
  # found first by the dispatch of ggplot2::autoplot() and loop
  if (!inherits(object, "WorkPrecisionSet") && requireNamespace("ggplot2", quietly = TRUE))
    return(ggplot2::autoplot(object, ...))
  UseMethod("autoplot")
}

#' @rdname autoplot
#' @export
#' @exportS3Method ggplot2::autoplot
autoplot.WorkPrecisionSet <- function(object, families = NULL, reference_tags = NULL, best_n = 2L, ...) {
  families <- if (is.null(families)) {
    u <- unique_tags(object)
    u[tag_kind(u) == "family"]
  } else as.character(families)
  ref <- if (is.null(reference_tags)) integer() else
    which(vapply(object$wps, .dv_has_any, TRUE, tags = as.character(reference_tags)))
  out <- list()
  present <- character()
  for (fam in families) {
    idx <- which(vapply(object$wps, function(wp) fam %in% wp$tags, TRUE))
    if (!length(idx)) next
    present <- c(present, fam)
    out[[paste0("family_", fam)]] <- .dv_subset_wps(object, sort(union(idx, ref)))
  }
  if (length(present)) {
    best <- .dv_best_family_indices(object, present, best_n, "area")
    out[["best_of_families"]] <- .dv_subset_wps(object, sort(union(best, ref)))
  }
  out[["all"]] <- object
  out
}

# ---------------------------------------------------------------------------
# Autodiff variants

#' Setups with every automatic differentiation backend
#'
#' `with_autodiff_variants(setups, ad_backends)` returns the setups tagged
#' `<tag_prefix>_default`, each followed (for algorithms with an `autodiff`
#' option) by one variant per backend: the algorithm with that backend,
#' named `"<name> (<Backend>)"` and tagged `<tag_prefix>_<backend>` (such as
#' `autodiff_finitediff`).
#'
#' @param setups List of setups.
#' @param ad_backends List of backends such as [AutoFiniteDiff()].
#' @param tag_prefix Prefix of the tags.
#' @return A list of setups.
#' @examples
#' v <- with_autodiff_variants(list(list(alg = Rodas5P()), list(alg = Tsit5())),
#'   ad_backends = list(AutoFiniteDiff(), AutoForwardDiff()))
#' vapply(v, function(s) if (is.null(s$name)) "" else s$name, "")
#' @export
with_autodiff_variants <- function(setups, ad_backends, tag_prefix = "autodiff") {
  out <- list()
  for (setup in setups) {
    original <- setup
    original$tags <- c(.dv_setup_tags(setup), paste0(tag_prefix, "_default"))
    out[[length(out) + 1L]] <- original
    ctor <- get0(class(setup$alg)[1L], envir = asNamespace("modelcompR"), mode = "function")
    if (is.null(ctor) || !"autodiff" %in% names(formals(ctor))) next
    for (backend in ad_backends) {
      bname <- .dv_ad_backend_name(backend)
      variant <- setup
      variant$alg <- .with_autodiff(backend, setup$alg)
      variant$name <- paste0(.dv_setup_name(setup), " (", bname, ")")
      variant$tags <- c(.dv_setup_tags(setup), paste0(tag_prefix, "_", tolower(bname)))
      out[[length(out) + 1L]] <- variant
    }
  }
  out
}

.dv_ad_backend_name <- function(backend) {
  name <- class(backend)[1L]
  if (startsWith(name, "Auto")) substring(name, 5L) else name
}

# ---------------------------------------------------------------------------
# Work precision

.dv_wp_new <- function(prob, abstols, reltols, errors, times, dts = NULL, stats = NULL, name = "WP-Alg",
                       error_estimate = "final", N = length(abstols), tags = character()) {
  if (is.list(errors) && length(errors) && is.list(errors[[1L]])) errors <- .dv_dicts_to_soa(errors)
  names(errors) <- vapply(names(errors), .dv_key, "")
  structure(list(prob = prob, abstols = abstols, reltols = reltols, errors = lapply(errors, as.numeric),
                 times = as.numeric(times), dts = dts, stats = stats, name = name,
                 error_estimate = .dv_key(error_estimate), N = as.integer(N), tags = as.character(unlist(tags))),
            class = "WorkPrecision")
}

.dv_combined_tags <- function(alg, tags, include_auto) unique(c(if (isTRUE(include_auto)) auto_tags(alg), tags))

#' Work precision of an algorithm
#'
#' `WorkPrecision(prob, alg, abstols, reltols, dts)` solves the problem
#' with each pair of tolerances (and time step, if `dts` is given) and
#' records the errors (against the analytic solution, or against `appxsol`
#' through [appxtrue()]) and the best of `numruns` timings. A solve that
#' fails or takes longer than `timeout` seconds gets `NaN` errors and time.
#' The errors computed by default depend on `error_estimate`: timeseries
#' errors for `l2`, `l_inf`, `L2`, `L_inf` and dense errors for `L2`,
#' `L_inf`. For a nonlinear problem the error is the `l2` (default) or
#' `l_inf` norm of the residual, or of the difference with `appxsol`.
#'
#' `WorkPrecision(prob, abstols, reltols, errors, times, dts, stats, name,
#' error_estimate, N, tags)` builds one from its fields (with `errors` a
#' named list of vectors or a list of named lists, one per tolerance).
#'
#' `length(wp)` is the number of tolerances and `wp[[i]]` the `i`th time;
#' `plot(wp)` draws the time against the error on log scales.
#'
#' @param prob A problem, or a list of problems selected by the option
#'   `prob_choice`.
#' @param alg The algorithm.
#' @param abstols,reltols Tolerances.
#' @param dts Time steps, or `NULL`.
#' @param dt,p The solver options `dt` and `p`, named in the signature so
#'   that they do not partially match `dts` and `prob`.
#' @param abstol,reltol Ignored: each run uses its own tolerances, as in
#'   Julia (named so that they do not partially match `abstols` and
#'   `reltols`).
#' @param ... Options: `name` (default `"WP-Alg"`), `appxsol`,
#'   `error_estimate` (default `"final"`, `"l2"` for nonlinear problems),
#'   `numruns` (default 20), `seconds` (default 2), `timeout`,
#'   `timeseries_errors`, `dense_errors`, `tags`, `auto_tags` (default
#'   `TRUE`) and solver options; or, with the fields form, the remaining
#'   fields.
#' @return A `WorkPrecision`: a list with `prob`, `abstols`, `reltols`,
#'   `errors`, `times`, `dts`, `stats`, `name`, `error_estimate`, `N` and
#'   `tags`.
#' @examples
#' fn <- ODEFunction(function(u, p, t) -u, analytic = function(u0, p, t) u0 * exp(-t))
#' prob <- ODEProblem(fn, 1, c(0, 1))
#' wp <- WorkPrecision(prob, Tsit5(), c(1e-4, 1e-6), c(1e-4, 1e-6), numruns = 1)
#' wp$errors$final
#' @export
WorkPrecision <- function(prob, alg, abstols, reltols, dts = NULL, ..., dt = NULL, p = NULL, abstol = NULL,
                          reltol = NULL) {
  # dt, p, abstol and reltol are named so that they do not partially match
  # dts, prob, abstols and reltols; abstol and reltol are replaced by the
  # tolerances of each run, as in Julia
  if (missing(alg) || is.numeric(alg)) {
    fields <- c(list(if (!missing(alg)) alg, if (!missing(abstols)) abstols, if (!missing(reltols)) reltols,
                     dts), list(...))
    nm <- c("abstols", "reltols", "errors", "times", "dts", "stats", "name", "error_estimate", "N", "tags")
    named <- names(fields) %||% rep("", length(fields))
    pos <- which(!nzchar(named))
    names(fields)[pos] <- nm[setdiff(seq_along(nm), match(named[nzchar(named)], nm))][seq_along(pos)]
    fields <- fields[!vapply(fields, is.null, TRUE) | names(fields) %in% c("dts", "stats")]
    return(do.call(.dv_wp_new, c(list(prob), fields)))
  }
  .dv_wp_compute(prob, alg, abstols, reltols, dts, .dv_dots(list(...), dt, p))
}

# The work precision functions bind their settings with list2env, so these
# names are local bindings, not globals.
utils::globalVariables(c("appxsol", "appxsol_setup", "dense_errors", "ensemblealg", "error_estimates",
                         "expected_value", "numruns", "numruns_error", "parallel_type", "print_names",
                         "prob_choice", "seconds", "tags", "timeout", "timeseries_errors", "trajectories"))

.dv_wp_compute <- function(prob, alg, abstols, reltols, dts, opts) {
  st <- .dv_settings(opts, list(name = NULL, appxsol = NULL, error_estimate = NULL, numruns = 20L, seconds = 2,
                                timeout = NULL, timeseries_errors = NULL, dense_errors = NULL, tags = character(),
                                auto_tags = TRUE, prob_choice = NULL))
  list2env(st$s, environment())
  kwargs <- st$rest
  kwargs$alg <- NULL
  pc <- prob_choice
  p_ <- if (!is.null(pc)) prob[[pc]] else if (.dv_is_prob_list(prob)) prob[[1L]] else prob
  ap <- .dv_choose(appxsol, prob, pc)
  nonlinear <- inherits(p_, "NonlinearProblem")
  error_estimate <- .dv_key(error_estimate %||% if (nonlinear) "l2" else "final")
  name <- name %||% "WP-Alg"
  N <- length(abstols)
  errs <- vector("list", N)
  times <- numeric(N)
  stats <- vector("list", N)
  ts_err <- timeseries_errors %||% (error_estimate %in% .dv_timeseries_keys)
  dn_err <- dense_errors %||% (error_estimate %in% .dv_dense_keys)
  for (i in seq_len(N)) {
    o <- .dv_merge(kwargs, list(abstol = abstols[i], reltol = reltols[i]))
    if (!is.null(dts)) o$dt <- dts[i]
    if (!nonlinear) o <- .dv_merge(o, list(timeseries_errors = ts_err, dense_errors = dn_err))
    t0 <- .dv_now()
    sol <- .dv_solve(p_, alg, o)
    stats[i] <- list(sol$stats)
    if (.dv_timed_out(timeout, .dv_now() - t0, name, abstols[i])) {
      errs[[i]] <- if (nonlinear) stats::setNames(list(NaN), error_estimate) else .dv_nan_errors()
      times[i] <- NaN
      next
    }
    if (nonlinear) {
      err <- if (is.null(ap)) as.numeric(sol$resid) else
        as.numeric(sol$u) - as.numeric(if (inherits(ap, "SciMLSolution")) ap$u else ap)
      val <- switch(error_estimate, l2 = .dv_norm(err), l_inf = max(abs(err)),
                    stop("Unsupported norm used: ", error_estimate, ".", call. = FALSE))
      errs[[i]] <- stats::setNames(list(val), error_estimate)
      o_t <- o
    } else if (successful_retcode(sol)) {
      e <- if (!is.null(ap)) appxtrue(sol, ap)$errors else sol$errors
      errs[[i]] <- lapply(e, .dv_mean)
      o_t <- .dv_merge(o, list(timeseries_errors = FALSE, dense_errors = FALSE))
    } else {
      errs[[i]] <- .dv_nan_errors()
      times[i] <- NaN
      next
    }
    times[i] <- .dv_best_time(function() .dv_solve(p_, alg, o_t), numruns, seconds)
  }
  .dv_wp_new(prob, abstols, reltols, .dv_dicts_to_soa(errs), times, dts, stats, name, error_estimate, N,
             .dv_combined_tags(alg, tags, auto_tags))
}

#' @export
length.WorkPrecision <- function(x) .subset2(x, "N")

#' @export
`[[.WorkPrecision` <- function(x, i, ...) {
  if (is.numeric(i)) return(.subset2(x, "times")[[i]])
  .subset2(x, i)
}

#' @export
print.WorkPrecision <- function(x, ...) {
  cat("Name: ", x$name, "\n", sep = "")
  cat("Times: ", paste(format(x$times, digits = 4), collapse = " "), "\n", sep = "")
  cat("Errors:\n")
  for (k in names(x$errors)) cat("  ", k, ": ", paste(format(x$errors[[k]], digits = 4), collapse = " "), "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Work precision sets

.dv_wps_new <- function(wps, N, abstols, reltols, prob, setups, names, error_estimate, numruns = NULL,
                        active_error_estimates = NULL) {
  error_estimate <- .dv_key(error_estimate)
  structure(list(wps = wps, N = as.integer(N), abstols = abstols, reltols = reltols, prob = prob, setups = setups,
                 names = names, error_estimate = error_estimate, numruns = numruns,
                 active_error_estimates = vapply(as.character(active_error_estimates %||% error_estimate), .dv_key, "",
                                                 USE.NAMES = FALSE)),
            class = "WorkPrecisionSet")
}

#' Work precision of several algorithms
#'
#' `WorkPrecisionSet(prob, abstols, reltols, setups)` runs [WorkPrecision()]
#' for each setup (a list with the algorithm `alg` and optionally its own
#' `abstols`, `reltols`, `dts`, a `name`, `tags`, `auto_tags`,
#' `prob_choice` and solver options), with the options in `...` (`appxsol`,
#' `error_estimate`, `numruns`, `seconds`, `timeout`, `names`,
#' `print_names` and solver options). With `error_estimates`, a vector of
#' error names, the errors needed by all of them are computed in the same
#' runs; `available_errors(set)` lists them.
#'
#' For an SDE or RODE problem the errors are means over `numruns_error`
#' runs (against the analytic solution, or against a reference solved with
#' `appxsol_setup` on a Brownian grid of step `test_dt`), and the times are
#' the means of `numruns` runs. For an [EnsembleProblem()] the errors are
#' weak errors of `trajectories` paths (with `ensemblealg`): against
#' `expected_value` (`weak_final` or `weak_l2`) or against the ensemble
#' solved with `appxsol_setup`.
#'
#' `WorkPrecisionSet(wps, N, abstols, reltols, prob, setups, names,
#' error_estimate, numruns, active_error_estimates)` builds one from its
#' fields.
#'
#' `length(set)` is the number of entries and `set[[i]]` the `i`th
#' `WorkPrecision`. `plot(set)` draws the time (`y`) against the error
#' (`xvar`, the error estimate by default) of every entry on log scales;
#' `xvar` and `y` may also name tolerances (`"abstols"`, `"reltols"`),
#' `"times"`, `"dts"`, an error or a field of the solver statistics. The
#' options `tags`, `include_tags` and `exclude_tags` select entries,
#' `reference_tags` draws the entries with those tags as thin dashed
#' references, and `view = "dt_convergence"` draws the errors against the
#' time steps with fitted orders.
#'
#' @param prob The problem (a list of problems with `prob_choice`), or for
#'   the fields form the list of `WorkPrecision`s.
#' @param abstols,reltols Tolerances.
#' @param setups List of setups.
#' @param test_dt Time step of the Brownian grid of reference solutions.
#' @param p The solver option `p`, named in the signature so that it does
#'   not partially match `prob`.
#' @param abstol,reltol Ignored: each run uses its own tolerances.
#' @param ... Options, see Details.
#' @param wp_set A `WorkPrecisionSet`.
#' @return A `WorkPrecisionSet`: a list with `wps`, `N`, `abstols`,
#'   `reltols`, `prob`, `setups`, `names`, `error_estimate`, `numruns` and
#'   `active_error_estimates`; `available_errors()` returns error names.
#' @examples
#' fn <- ODEFunction(function(u, p, t) -u, analytic = function(u0, p, t) u0 * exp(-t))
#' prob <- ODEProblem(fn, 1, c(0, 1))
#' setups <- list(list(alg = Euler(), dts = c(1 / 8, 1 / 16)),
#'   list(alg = RK4(), dts = c(1 / 4, 1 / 8)))
#' set <- WorkPrecisionSet(prob, c(1e-3, 1e-4), c(1e-3, 1e-4), setups, numruns = 1,
#'   error_estimates = c("final", "l2"))
#' available_errors(set)
#' set[[2]]$errors$l2
#' @export
WorkPrecisionSet <- function(prob, abstols, reltols, setups, test_dt = NULL, ..., p = NULL, abstol = NULL,
                             reltol = NULL) {
  if (is.list(prob) && !inherits(prob, "SciMLProblem") && all(vapply(prob, inherits, TRUE, "WorkPrecision"))) {
    fields <- c(list(wps = prob, N = abstols, abstols = reltols, reltols = setups, prob = test_dt), list(...))
    pos <- which(!nzchar(names(fields)))
    names(fields)[pos] <- c("setups", "names", "error_estimate", "numruns", "active_error_estimates")[seq_along(pos)]
    return(do.call(.dv_wps_new, fields))
  }
  p1 <- if (.dv_is_prob_list(prob)) prob[[1L]] else prob
  opts <- .dv_dots(list(...), NULL, p)
  if (inherits(prob, "EnsembleProblem")) return(.dv_wps_ensemble(prob, abstols, reltols, setups, test_dt, opts))
  if (inherits(p1, c("SDEProblem", "RODEProblem")) && !.dv_is_prob_list(prob))
    return(.dv_wps_rode(prob, abstols, reltols, setups, test_dt, opts))
  .dv_wps_generic(prob, abstols, reltols, setups, opts)
}

.dv_wps_generic <- function(prob, abstols, reltols, setups, opts) {
  st <- .dv_settings(opts, list(print_names = FALSE, names = NULL, appxsol = NULL, error_estimate = "final",
                                error_estimates = NULL, timeout = NULL))
  list2env(st$s, environment())
  N <- length(setups)
  if (!is.null(names) && length(names) != N) stop("give one name per setup", call. = FALSE)
  names <- names %||% vapply(setups, .dv_setup_name, "")
  error_estimate <- .dv_key(error_estimate)
  active <- if (is.null(error_estimates)) error_estimate else vapply(as.character(error_estimates), .dv_key, "",
                                                                      USE.NAMES = FALSE)
  ts_err <- if (!is.null(error_estimates)) any(active %in% .dv_timeseries_keys)
  dn_err <- if (!is.null(error_estimates)) any(active %in% .dv_dense_keys)
  kwargs <- st$rest
  wps <- lapply(seq_len(N), function(i) {
    s <- setups[[i]]
    if (isTRUE(print_names)) cat(names[i], "\n")
    opts <- list(appxsol = appxsol, error_estimate = error_estimate, timeout = timeout,
                 timeseries_errors = ts_err, dense_errors = dn_err, tags = .dv_setup_tags(s),
                 auto_tags = s$auto_tags %||% TRUE, name = names[i])
    opts <- .dv_merge(.dv_merge(opts, kwargs), .dv_setup_options(s))
    if (!is.null(s$prob_choice)) opts$prob_choice <- s$prob_choice
    .dv_wp_compute(prob, s$alg, s$abstols %||% abstols, s$reltols %||% reltols, s$dts, opts)
  })
  .dv_wps_new(wps, N, abstols, reltols, prob, setups, names, error_estimate, NULL, active)
}

# Solve with a setup at its jth tolerance (dt = 0, the default of setups
# without dts, lets the solver choose the step).
.dv_setup_solve <- function(prob, setup, abstols, reltols, j, kwargs, ts_err, dn_err, extra = list()) {
  ab <- setup$abstols %||% abstols
  rl <- setup$reltols %||% reltols
  dts <- setup$dts %||% numeric(length(ab))
  o <- .dv_merge(.dv_merge(kwargs, .dv_setup_options(setup)),
                 c(list(abstol = ab[j], reltol = rl[j], dt = dts[j], timeseries_errors = ts_err,
                        dense_errors = dn_err), extra))
  o$prob_choice <- NULL
  if (identical(o$dt, 0)) o$dt <- NULL
  o
}

# A NoiseGrid of two Brownian paths with step test_dt, one value per state.
.dv_noise_grid <- function(prob, test_dt) {
  t <- seq(prob$tspan[1L], prob$tspan[2L], by = test_dt)
  m <- length(prob$u0)
  path <- function() Reduce(`+`, lapply(seq_len(length(t) - 1L), function(i) sqrt(test_dt) * stats::rnorm(m)),
                            numeric(m), accumulate = TRUE)
  W <- path()
  NoiseGrid(t, W, path())
}

.dv_with_noise <- function(prob, np) {
  prob$noise_process <- np
  if (inherits(prob, "RODEProblem")) prob$noise <- np
  prob
}

# _calculate_error!: one run of every setup at every tolerance, with errors
# against the analytic solution or a reference on a Brownian grid.
.dv_calculate_error <- function(prob, setups, abstols, reltols, test_dt, appxsol_setup, ts_err, dn_err, kwargs) {
  has_an <- !is.null(.dv_analytic(prob))
  if (!has_an) {
    p_ <- .dv_with_noise(prob, .dv_noise_grid(prob, test_dt))
    s <- as.list(appxsol_setup)
    a <- s$alg
    s$alg <- NULL
    true_sol <- do.call(solve, c(list(p_, a), .dv_merge(kwargs, s)))
  } else {
    p_ <- prob
  }
  .dv_solve(p_, setups[[1L]]$alg, .dv_setup_solve(p_, setups[[1L]], abstols, reltols, 1L, kwargs, FALSE, FALSE))
  lapply(seq_along(abstols), function(j) lapply(seq_along(setups), function(k) {
    sol <- .dv_solve(p_, setups[[k]]$alg, .dv_setup_solve(p_, setups[[k]], abstols, reltols, j, kwargs, ts_err, dn_err))
    if (has_an) sol else appxtrue(sol, true_sol)
  }))
}

.dv_wps_rode <- function(prob, abstols, reltols, setups, test_dt, opts) {
  st <- .dv_settings(opts, list(numruns = 20L, numruns_error = 20L, print_names = FALSE, names = NULL,
                                appxsol_setup = NULL, error_estimate = "final", error_estimates = NULL,
                                parallel_type = "none", timeout = NULL))
  list2env(st$s, environment())
  if (!is.null(names) && length(names) != length(setups)) stop("give one name per setup", call. = FALSE)
  kwargs <- st$rest
  error_estimate <- .dv_key(error_estimate)
  has_an <- !is.null(.dv_analytic(prob))
  ts_err <- has_an && error_estimate %in% .dv_timeseries_keys
  dn_err <- has_an && error_estimate %in% .dv_dense_keys
  weak_ts <- error_estimate %in% c("weak_l2", "weak_l_inf")
  weak_dense <- error_estimate %in% c("weak_L2", "weak_L_inf")
  N <- length(setups)
  M <- length(abstols)
  names <- names %||% vapply(setups, .dv_setup_name, "")
  run <- function(i) .dv_calculate_error(prob, setups, abstols, reltols, test_dt, appxsol_setup, ts_err, dn_err,
                                         kwargs)
  runs <- if (identical(parallel_type, "threads")) {
    parallel::mclapply(seq_len(numruns_error), run, mc.cores = getOption("mc.cores", 2L))
  } else {
    lapply(seq_len(numruns_error), function(i) {
      message("Error calculation: ", i, "/", numruns_error)
      run(i)
    })
  }
  sims <- lapply(seq_len(N), function(k) lapply(seq_len(M), function(j)
    .dv_ensemble_errors(lapply(runs, function(r) r[[j]][[k]]), 0, TRUE, weak_ts, weak_dense)))
  field <- if (error_estimate %in% .dv_weak_keys) "weak_errors" else "error_means"
  times <- matrix(NA_real_, M, N)
  for (k in seq_len(N)) {
    s <- setups[[k]]
    ab <- s$abstols %||% abstols
    sol <- .dv_solve(prob, s$alg, .dv_setup_solve(prob, s, abstols, reltols, 1L, kwargs, FALSE, FALSE))
    x <- if (!length(sol$t)) 0 else round(mean(sol$t) - sum(sol$t) / length(sol$t))
    for (j in seq_len(M)) {
      o <- .dv_setup_solve(prob, s, abstols, reltols, j, kwargs, FALSE, FALSE)
      tt <- numeric(numruns)
      out <- FALSE
      for (i in seq_len(numruns)) {
        t0 <- .dv_now()
        .dv_solve(prob, s$alg, o)
        tt[i] <- .dv_now() - t0
        if (.dv_timed_out(timeout, tt[i], names[k], ab[j])) {
          out <- TRUE
          break
        }
      }
      times[j, k] <- if (out) NaN else .dv_mean(tt) + x
    }
  }
  wps <- lapply(seq_len(N), function(k) {
    s <- setups[[k]]
    # N is the number of tolerances (Julia passes the number of setups)
    ab <- s$abstols %||% abstols
    .dv_wp_new(prob, ab, s$reltols %||% reltols, lapply(sims[[k]], function(sim) .subset2(sim, field)), times[, k],
               s$dts %||% numeric(length(ab)), NULL, names[k], error_estimate, M,
               .dv_combined_tags(s$alg, .dv_setup_tags(s), s$auto_tags %||% TRUE))
  })
  .dv_wps_new(wps, N, abstols, reltols, prob, setups, names, error_estimate, numruns_error,
              error_estimates %||% error_estimate)
}

# The flat values of an ensemble output (a solution or values)
.dv_flat <- function(v) if (inherits(v, "SciMLSolution")) unlist(v$u) else as.numeric(unlist(v))

.dv_wps_ensemble <- function(prob, abstols, reltols, setups, test_dt, opts) {
  st <- .dv_settings(opts, list(numruns = 5L, trajectories = 1000L, print_names = FALSE, names = NULL,
                                appxsol_setup = NULL, expected_value = NULL, error_estimate = "weak_final",
                                error_estimates = NULL, ensemblealg = EnsembleThreads(), timeout = NULL))
  list2env(st$s, environment())
  if (!is.null(names) && length(names) != length(setups)) stop("give one name per setup", call. = FALSE)
  kwargs <- st$rest
  error_estimate <- .dv_key(error_estimate)
  N <- length(setups)
  M <- length(abstols)
  names <- names %||% vapply(setups, .dv_setup_name, "")
  esolve <- function(s, j, extra = list()) {
    o <- .dv_setup_solve(prob, s, abstols, reltols, j, list(), FALSE, FALSE,
                         c(list(trajectories = as.integer(trajectories)), kwargs, extra))
    o$timeseries_errors <- NULL
    o$dense_errors <- NULL
    do.call(solve, c(list(prob, s$alg, ensemblealg), o))
  }
  sols <- lapply(seq_len(N), function(k) {
    out <- lapply(seq_len(M), function(j) esolve(setups[[k]], j))
    message(.dv_default_name(setups[[k]]$alg), " (", k, "/", N, ")")
    out
  })
  if (!error_estimate %in% .dv_weak_keys)
    stop("use RODEProblem instead of EnsembleProblem for strong errors.", call. = FALSE)
  errors <- if (!is.null(expected_value)) {
    lapply(sols, function(sk) vapply(sk, function(sol) {
      d <- lapply(sol$u, function(v) .dv_flat(v) - expected_value)
      switch(error_estimate,
             weak_final = .dv_norm(.dv_vec_mean(d)),
             weak_l2 = abs(.dv_mean(unlist(d))),
             stop("Error estimate ", error_estimate, " is not implemented yet.", call. = FALSE))
    }, 0))
  } else {
    s <- as.list(appxsol_setup)
    a <- s$alg
    s$alg <- NULL
    ref <- do.call(solve, c(list(prob, a, ensemblealg), .dv_merge(kwargs, s),
                            list(trajectories = as.integer(trajectories))))
    lapply(sols, function(sk) vapply(sk, function(sol)
      .dv_norm(.dv_vec_mean(Map(function(a, b) .dv_flat(a) - .dv_flat(b), sol$u, ref$u))), 0))
  }
  times <- matrix(NA_real_, M, N)
  for (k in seq_len(N)) {
    s <- setups[[k]]
    ab <- s$abstols %||% abstols
    esolve(s, 1L)
    for (j in seq_len(M)) {
      tt <- numeric(numruns)
      out <- FALSE
      for (i in seq_len(numruns)) {
        t0 <- .dv_now()
        esolve(s, j)
        tt[i] <- .dv_now() - t0
        if (.dv_timed_out(timeout, tt[i], names[k], ab[j])) {
          out <- TRUE
          break
        }
      }
      times[j, k] <- if (out) NaN else .dv_mean(tt)
    }
  }
  wps <- lapply(seq_len(N), function(k) {
    s <- setups[[k]]
    .dv_wp_new(prob, s$abstols %||% abstols, s$reltols %||% reltols,
               stats::setNames(list(errors[[k]]), error_estimate), times[, k],
               s$dts %||% numeric(length(s$abstols %||% abstols)), NULL, names[k], error_estimate, M,
               .dv_combined_tags(s$alg, .dv_setup_tags(s), s$auto_tags %||% TRUE))
  })
  .dv_wps_new(wps, N, abstols, reltols, prob, setups, names, error_estimate, as.integer(trajectories),
              error_estimates %||% error_estimate)
}

#' @rdname WorkPrecisionSet
#' @export
available_errors <- function(wp_set) wp_set$active_error_estimates

#' @export
length.WorkPrecisionSet <- function(x) .subset2(x, "N")

#' @export
`[[.WorkPrecisionSet` <- function(x, i, ...) {
  if (is.numeric(i)) return(.subset2(x, "wps")[[i]])
  .subset2(x, i)
}

#' @export
print.WorkPrecisionSet <- function(x, ...) {
  cat("WorkPrecisionSet of ", x$N, " wps\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Sample errors

#' Sampling error of SDE benchmarks
#'
#' `get_sample_errors(prob, setup, test_dt, numruns = n, solution_runs = m)`
#' estimates the sampling error of the mean final value over `n` runs of an
#' SDE algorithm: `1.96` times the standard deviation, over `m` repetitions,
#' of the distance between the mean final value and the analytic mean (the
#' mean of `analytic(u0, p, T, W)` over `sample_error_runs` normal `W` with
#' variance `T - t0`, or without an analytic solution the mean over all
#' runs), divided by `sqrt(n)`. With a vector `numruns` it returns one value
#' per entry. Runs without an analytic solution are compared with a
#' reference solved by `appxsol_setup` on a Brownian grid of step `test_dt`.
#'
#' @param prob An SDE or RODE problem.
#' @param setup A setup.
#' @param test_dt Step of the Brownian grid of reference solutions.
#' @param appxsol_setup Setup of the reference.
#' @param numruns Number of runs, or a vector of them.
#' @param error_estimate Unused, as in Julia.
#' @param sample_error_runs Number of samples of the analytic mean.
#' @param solution_runs Number of repetitions.
#' @param parallel_type `"none"` or `"threads"` (forked workers).
#' @param ... Solver options.
#' @param p The solver option `p`, named in the signature so that it does
#'   not partially match `prob` and `parallel_type`.
#' @return A number or a vector.
#' @examples
#' f <- function(u, p, t) 0.5 * u
#' g <- function(u, p, t) 0.2 * u
#' an <- function(u0, p, t, W) u0 * exp(0.48 * t + 0.2 * W)
#' prob <- SDEProblem(SDEFunction(f, g, analytic = an), 1, c(0, 1))
#' set.seed(1)
#' get_sample_errors(prob, list(alg = EM(), dts = 1 / 4), numruns = 3, solution_runs = 2,
#'   sample_error_runs = 100)
#' @export
get_sample_errors <- function(prob, setup, test_dt = NULL, appxsol_setup = NULL, numruns, error_estimate = "final",
                              sample_error_runs = 1e7, solution_runs, parallel_type = "none", ..., p = NULL) {
  kwargs <- .dv_dots(list(...), NULL, p)
  maxnumruns <- max(numruns)
  full <- lapply(seq_len(solution_runs), function(i) {
    message("Solution Run: ", i)
    run <- function(r) .dv_calculate_error(prob, list(setup), 1e-2, 1e-2, test_dt, appxsol_setup, FALSE, FALSE,
                                           kwargs)[[1L]][[1L]]
    if (identical(parallel_type, "threads"))
      parallel::mclapply(seq_len(maxnumruns), run, mc.cores = getOption("mc.cores", 2L)) else
      lapply(seq_len(maxnumruns), run)
  })
  last <- function(s) as.numeric(.dv_last(s$u))
  if (!is.null(.dv_analytic(prob))) {
    dt_ <- prob$tspan[2L] - prob$tspan[1L]
    an <- prob$f$analytic
    analytical_mean_end <- .dv_vec_mean(lapply(seq_len(sample_error_runs), function(i) {
      W <- sqrt(dt_) * stats::rnorm(length(prob$u0))
      an(prob$u0, prob$p, prob$tspan[2L], W)
    }))
  } else {
    analytical_mean_end <- .dv_vec_mean(lapply(full, function(runs) .dv_vec_mean(lapply(runs, last))))
  }
  one <- function(n) {
    ends <- lapply(full, function(runs) .dv_vec_mean(lapply(runs[seq_len(n)], last)))
    1.96 * stats::sd(vapply(ends, function(m) .dv_norm(m - analytical_mean_end), 0)) / sqrt(n)
  }
  if (length(numruns) == 1L) one(numruns) else vapply(numruns, one, 0)
}

# ---------------------------------------------------------------------------
# Shootouts

#' Shootouts between algorithms
#'
#' `Shootout(prob, setups)` solves the problem with each setup and compares
#' the efficiencies \eqn{1 / (error \cdot time)} (the error `error_estimate`
#' against the analytic solution or `appxsol`, the time the best of
#' `numruns` runs). `effratios[i, j]` is the efficiency of `i` relative to
#' `j`, and `winner` the name of the most efficient setup.
#' `ShootoutSet(probs, setups)` runs a shootout for each problem of a list
#' (with the options of `probaux[[i]]` for problem `i`).
#'
#' @param prob A problem.
#' @param probs A list of problems.
#' @param p The solver option `p`, named in the signature so that it does
#'   not partially match `prob` or `probs`.
#' @param setups List of setups.
#' @param ... Options: `appxsol`, `names`, `error_estimate` (default
#'   `"final"`), `numruns` (default 20), `seconds` (default 2) and solver
#'   options; for `ShootoutSet()` also `probaux` and `print_names`.
#' @param x A `Shootout`.
#' @param xlab,ylab Axis labels.
#' @return A `Shootout` (a list with `setups`, `times`, `errors`, `effs`,
#'   `effratios`, `solutions`, `names`, `N`, `bestidx` and `winner`) or a
#'   `ShootoutSet` (with `shootouts`, `probs`, `probaux`, `N` and
#'   `winners`).
#' @examples
#' fn <- ODEFunction(function(u, p, t) -u, analytic = function(u0, p, t) u0 * exp(-t))
#' prob <- ODEProblem(fn, 1, c(0, 1))
#' sh <- Shootout(prob, list(list(alg = Euler()), list(alg = RK4())), dt = 1 / 16, numruns = 1)
#' sh$errors
#' @export
Shootout <- function(prob, setups, ..., p = NULL) .dv_shootout(prob, setups, .dv_dots(list(...), NULL, p))

.dv_shootout <- function(prob, setups, opts) {
  st <- .dv_settings(opts, list(appxsol = NULL, names = NULL, error_estimate = "final", numruns = 20L,
                                seconds = 2))
  list2env(st$s, environment())
  N <- length(setups)
  if (!is.null(names) && length(names) != N) stop("give one name per setup", call. = FALSE)
  error_estimate <- .dv_key(error_estimate)
  ts_err <- error_estimate %in% .dv_timeseries_keys
  dn_err <- error_estimate %in% .dv_dense_keys
  names <- names %||% vapply(setups, function(s) .dv_default_name(s$alg), "")
  kwargs <- st$rest
  errors <- times <- effs <- numeric(N)
  solutions <- vector("list", N)
  for (i in seq_len(N)) {
    s <- setups[[i]]
    pc <- s$prob_choice
    p_ <- if (!is.null(pc)) prob[[pc]] else if (.dv_is_prob_list(prob)) prob[[1L]] else prob
    o <- .dv_merge(.dv_merge(list(timeseries_errors = ts_err, dense_errors = dn_err), kwargs), .dv_setup_options(s))
    o$prob_choice <- NULL
    sol <- .dv_solve(p_, s$alg, o)
    ap <- .dv_choose(appxsol, prob, pc)
    sol <- if (!is.null(ap)) appxtrue(sol, ap) else sol
    errors[i] <- sol$errors[[error_estimate]] %||% NA_real_
    solutions[[i]] <- sol
    o_t <- .dv_merge(kwargs, list(timeseries_errors = FALSE, dense_errors = FALSE))
    times[i] <- .dv_best_time(function() .dv_solve(p_, s$alg, o_t), numruns, seconds)
    effs[i] <- 1 / (errors[i] * times[i])
  }
  bestidx <- which(effs == max(effs))[1L]
  structure(list(setups = setups, times = times, errors = errors, effs = effs, effratios = outer(effs, effs, "/"),
                 solutions = solutions, names = names, N = N, bestidx = bestidx, winner = names[bestidx]),
            class = "Shootout")
}

#' @rdname Shootout
#' @export
ShootoutSet <- function(probs, setups, ..., p = NULL) .dv_shootout_set(probs, setups, .dv_dots(list(...), NULL, p))

.dv_shootout_set <- function(probs, setups, opts) {
  st <- .dv_settings(opts, list(probaux = NULL, names = NULL, print_names = FALSE))
  list2env(st$s, environment())
  N <- length(probs)
  if (!is.null(names) && length(names) != length(setups)) stop("give one name per setup", call. = FALSE)
  names <- names %||% vapply(setups, function(s) .dv_default_name(s$alg), "")
  probaux <- probaux %||% rep(list(list()), N)
  shootouts <- lapply(seq_len(N), function(i) {
    if (isTRUE(print_names)) cat(names[i], "\n")
    .dv_shootout(probs[[i]], setups, .dv_merge(c(list(names = names), st$rest), probaux[[i]]))
  })
  structure(list(shootouts = shootouts, probs = probs, probaux = probaux, N = N,
                 winners = vapply(shootouts, function(s) s$winner, "")),
            class = "ShootoutSet")
}

#' @export
length.Shootout <- function(x) .subset2(x, "N")

#' @export
`[[.Shootout` <- function(x, i, ...) {
  if (is.numeric(i)) return(.subset2(x, "effs")[[i]])
  .subset2(x, i)
}

#' @export
print.Shootout <- function(x, ...) {
  cat("Winner: ", x$winner, "\n", sep = "")
  cat("EffRatios: ", paste(format(x$effratios[x$bestidx, ], digits = 4), collapse = " "), "\n", sep = "")
  invisible(x)
}

#' @export
length.ShootoutSet <- function(x) .subset2(x, "N")

#' @export
`[[.ShootoutSet` <- function(x, i, ...) {
  if (is.numeric(i)) return(.subset2(x, "shootouts")[[i]])
  .subset2(x, i)
}

#' @export
print.ShootoutSet <- function(x, ...) {
  cat("ShootoutSet of ", x$N, " shootouts\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Plots

.dv_key_label <- function(key) {
  switch(key,
    abstols = "Absolute tolerance", reltols = "Relative tolerance", times = "Time (s)", dts = "\u0394t",
    nf = "Number of function evaluations", nf2 = "Number of function 2 evaluations",
    nw = "Number of W matrix evaluations", nsolve = "Number of linear solves",
    njacs = "Number of Jacobians created", nnonliniter = "Number of nonlinear solver iterations",
    nnonlinconvfail = "Number of nonlinear solver convergence failures",
    ncondition = "Number of callback condition calls", naccept = "Number of accepted steps",
    nreject = "Number of rejected steps", maxeig = "Maximum eigenvalue recorded",
    if (key %in% c("final", .dv_timeseries_keys, .dv_weak_keys)) paste0("Error (", key, ")") else key)
}

.dv_wp_value <- function(wp, key) {
  if (key %in% c("abstols", "reltols", "times", "dts")) return(wp[[key]])
  st <- wp$stats
  if (length(st) && is.list(st[[1L]]) && key %in% names(st[[1L]]))
    return(vapply(st, function(s) as.numeric(s[[key]]), 0))
  if (key %in% names(wp$errors)) return(wp$errors[[key]])
  stop("Key ", key, " does not specify a valid property of WorkPrecision", call. = FALSE)
}

#' @rdname WorkPrecision
#' @param x A `WorkPrecision`.
#' @param xlab,ylab,log,type,lwd Plot settings.
#' @export
plot.WorkPrecision <- function(x, xlab = paste0("Error (", x$error_estimate, ")"), ylab = "Time (s)", log = "xy",
                               type = "b", lwd = 3, ...) {
  graphics::plot(x$errors[[x$error_estimate]], x$times, xlab = xlab, ylab = ylab, log = log, type = type,
                 lwd = lwd, main = x$name, ...)
  invisible(x)
}

#' @rdname WorkPrecisionSet
#' @param x A `WorkPrecisionSet`.
#' @param y,xvar Keys of the vertical and horizontal axes.
#' @param view `"benchmark"` or `"dt_convergence"`.
#' @param col Colors of the entries.
#' @param tags,include_tags,exclude_tags,reference_tags Tag selection.
#' @param reference_style Line type, width and transparency of references.
#' @param legend_pos Position of the legend.
#' @export
plot.WorkPrecisionSet <- function(x, y = "times", xvar = x$error_estimate, view = "benchmark", col = NULL,
                                  tags = NULL, include_tags = NULL, exclude_tags = NULL, reference_tags = NULL,
                                  reference_style = list(lty = 2, lwd = 1, alpha = 0.5), legend_pos = "topright",
                                  ...) {
  if (!is.null(tags) || !is.null(include_tags) || !is.null(exclude_tags))
    x <- .dv_subset_wps(x, .dv_selected_indices(x, tags, include_tags, exclude_tags))
  n <- length(x$wps)
  col <- col %||% seq_len(n)
  col <- rep_len(col, n)
  if (identical(view, "benchmark")) {
    xs <- lapply(x$wps, .dv_wp_value, key = .dv_key(xvar))
    ys <- lapply(x$wps, .dv_wp_value, key = .dv_key(y))
    ref <- if (is.null(reference_tags)) integer() else
      which(vapply(x$wps, .dv_has_any, TRUE, tags = as.character(reference_tags)))
    pos <- function(v) v[is.finite(v) & v > 0]
    graphics::plot(range(pos(unlist(xs))), range(pos(unlist(ys))), type = "n", log = "xy",
                   xlab = .dv_key_label(.dv_key(xvar)), ylab = .dv_key_label(.dv_key(y)), ...)
    lty <- rep(1, n); lwd <- rep(3, n)
    if (length(ref) && length(ref) < n) {
      lty[ref] <- reference_style$lty %||% 2
      lwd[ref] <- reference_style$lwd %||% 1
      col[ref] <- grDevices::adjustcolor(col[ref], alpha.f = reference_style$alpha %||% 0.5)
    }
    order_ <- c(ref, setdiff(seq_len(n), ref))
    for (i in order_) graphics::lines(xs[[i]], ys[[i]], type = "b", col = col[i], lty = lty[i], lwd = lwd[i], pch = i)
    graphics::legend(legend_pos, legend = vapply(x$wps, function(wp) wp$name, ""), col = col, lty = lty,
                     pch = seq_len(n), bty = "n")
  } else if (identical(view, "dt_convergence")) {
    idts <- which(vapply(x$setups, function(s) !is.null(s$dts), TRUE))
    if (!length(idts)) stop("Convergence with respect to \u0394t requires runs with fixed time steps", call. = FALSE)
    dts <- lapply(idts, function(i) x$setups[[i]]$dts)
    errs <- lapply(idts, function(i) x$wps[[i]]$errors[[x$error_estimate]])
    fits <- Map(function(d, e) qr.solve(cbind(1, log(d)), log(e)), dts, errs)
    labels <- paste0(x$names[idts], " (\u0394t\u1D56 order p=", vapply(fits, function(f) signif(f[2L], 2), 0), ")")
    graphics::plot(range(unlist(dts)), range(c(unlist(errs), unlist(Map(function(d, f) exp(f[1L]) * d^f[2L], dts, fits)))),
                   type = "n", log = "xy", xlab = "\u0394t", ylab = paste0("Error (", x$error_estimate, ")"), ...)
    for (m in seq_along(idts)) {
      graphics::lines(dts[[m]], exp(fits[[m]][1L]) * dts[[m]]^fits[[m]][2L], lty = 2, lwd = 1, col = col[m])
      graphics::lines(dts[[m]], errs[[m]], type = "b", lwd = 3, col = col[m], pch = m)
    }
    graphics::legend(legend_pos, legend = labels, col = col[seq_along(idts)], lty = 1, pch = seq_along(idts),
                     bty = "n")
  } else {
    stop("view argument `", view, "` not implemented", call. = FALSE)
  }
  invisible(x)
}

#' @rdname Shootout
#' @export
plot.Shootout <- function(x, xlab = "Algorithms", ylab = "Efficiency", ...) {
  graphics::barplot(x$effs, names.arg = x$names, xlab = xlab, ylab = ylab, ...)
  invisible(x)
}

