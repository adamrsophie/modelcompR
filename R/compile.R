# mtkcompile, after ModelingToolkitBase/src/systems/systems.jl.
#
# The entry point validates and normalizes the system, expands connections,
# and then runs either the full structural pipeline of ModelingToolkit
# (alias elimination, index reduction, state selection and tearing; see
# structural-*.R) or, when `structural = FALSE`, the simpler
# ModelingToolkitBase pipeline ported here (order lowering and trivial
# observed equations only).

# The separator used to name lowered derivatives: D(x) becomes xˍt.
DIFF2TERM_SEPARATOR <- "\u02CD"

#' Simplify a system for simulation
#'
#' `mtkcompile(sys)` turns a model into a system ready for building numerical
#' problems: it expands connections, flattens the hierarchy, eliminates
#' aliases and trivially solvable variables, reduces the index of
#' differential algebraic systems, and selects states, leaving the eliminated
#' variables as observed equations. `structural_simplify()` is its former
#' name.
#'
#' @param sys A `System`.
#' @param inputs,outputs,disturbance_inputs Variables to treat as inputs
#'   (turned into parameters) and outputs.
#' @param simplify Simplify the resulting expressions.
#' @param fully_determined Require as many equations as unknowns.
#' @param structural Run the full structural pipeline (`TRUE`) or only the
#'   ModelingToolkitBase lowering (`FALSE`).
#' @param additional_passes A list of functions applied to the result.
#' @param split Unused; kept for compatibility.
#' @param ... Options passed to the structural pipeline.
#' @return A complete, scheduled `System`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)", "y(t)"); parameters("tau = 3")
#' sys <- System(list(D(x) ~ (y - x) / tau, y ~ 2 * x), t, name = "sys")
#' csys <- mtkcompile(sys)
#' unknowns(csys); observed(csys)
#' @export
mtkcompile <- function(sys, inputs = list(), outputs = list(), disturbance_inputs = list(),
                       simplify = FALSE, fully_determined = TRUE, structural = TRUE,
                       additional_passes = list(), split = TRUE, ...) {
  if (!inherits(sys, "System")) stop("mtkcompile() needs a System", call. = FALSE)
  if (isTRUE(.f(sys, "isscheduled")))
    stop("this system was already compiled; call mtkcompile() on the original model", call. = FALSE)
  original <- sys
  inputs <- .flatten_vars(inputs)
  outputs <- .flatten_vars(outputs)
  disturbance_inputs <- .flatten_vars(disturbance_inputs)
  # _mtkcompile calls noise_to_brownians on the unflattened system
  if (!is.null(.f(sys, "noise_eqs")) && length(.f(sys, "systems")) && !length(jumps(sys)) &&
      !length(.collect_hier(sys, function(s) .f(s, "poissonians"))))
    stop("The system must be flattened.", call. = FALSE)
  # extract_top_level_statemachines (structural-statemachines.R)
  sm <- .sm_extract(sys)
  sys <- expand_connections(sm$sys)
  sys <- flatten(sys)
  if (length(sm$machines)) {
    md <- .f(sys, "metadata") %||% list()
    md$mtk_statemachines <- lapply(sm$machines, .sm_prepare)
    sys <- .update_system(sys, metadata = md)
  }
  extra <- list(...)
  compile_fn <- function(sys) {
    full <- .mtk$hooks$structural_compile
    if (structural && !is.null(full) && !is.null(.f(sys, "iv")) && !isTRUE(.f(sys, "is_discrete")))
      do.call(full, c(list(sys, inputs = inputs, outputs = outputs,
                           disturbance_inputs = disturbance_inputs,
                           fully_determined = fully_determined, simplify = simplify), extra))
    else
      .base_compile(sys, inputs = inputs, outputs = outputs, disturbance_inputs = disturbance_inputs,
                    fully_determined = fully_determined)
  }
  if (!is.null(.mtk$hooks$compile_clocked)) compile_fn <- .mtk$hooks$compile_clocked(compile_fn)
  newsys <- .lower_brownians_hook(sys, compile_fn)
  if (simplify) newsys <- .simplify_system(newsys)
  for (pass in additional_passes) newsys <- pass(newsys)
  newsys <- .update_system(newsys, parent = toggle_namespacing(original, FALSE))
  complete(newsys)
}

#' @rdname mtkcompile
#' @export
structural_simplify <- function(sys, ...) mtkcompile(sys, ...)

# Hook for the stochastic module, as _mtkcompile in systems.jl: it lowers
# poissonians to jumps and brownians to noise equations, and calls `compile`
# (the rest of the pipeline) unless the system has jumps.
.lower_brownians_hook <- function(sys, compile) {
  h <- .mtk$hooks$lower_brownians
  stochastic <- length(.f(sys, "brownians")) || length(.f(sys, "poissonians")) ||
    length(.f(sys, "jumps")) || !is.null(.f(sys, "noise_eqs"))
  if (is.null(h) || !stochastic) return(compile(sys))
  h(sys, compile)
}

# Hook for connectors; plain flattening without the connectors module.
#' Expand connections
#'
#' Replaces `connect()` statements by the equations they imply: equality of
#' potential variables and conservation of flow variables. Called by
#' [mtkcompile()].
#'
#' @param sys A `System`.
#' @param ... Options for the connectors module.
#' @return A `System` without connection statements.
#' @export
expand_connections <- function(sys, ...) {
  h <- .mtk$hooks$expand_connections
  if (!is.null(h)) return(h(sys, ...))
  sys
}

.simplify_system <- function(sys) {
  simp <- function(e) Equation(e$lhs, symcompR::simplify(e$rhs))
  .update_system(sys, eqs = lapply(.f(sys, "eqs"), simp), observed = lapply(.f(sys, "observed"), simp))
}

# The lowered name of the order-th derivative of var: x, xˍt, xˍtt, ...
.diff2term <- function(var, iv, order) {
  if (order == 0L) return(var)
  nm <- .name(var)
  ivn <- .name(iv)
  sep <- if (grepl(DIFF2TERM_SEPARATOR, nm, fixed = TRUE)) "" else DIFF2TERM_SEPARATOR
  newname <- paste0(nm, sep, strrep(ivn, order))
  meta <- .meta(var)
  meta$default <- NULL
  meta$kind <- NULL
  meta$guess <- NULL
  out <- .make_var(newname, "unknown", .args(var), meta)
  out
}

#' Lowered names of derivatives
#'
#' `default_toterm(D(x))` returns the variable \ifelse{latex}{\verb{x_t}}{\verb{xˍt}} that `mtkcompile()` uses
#' for a lowered derivative; other expressions are returned unchanged.
#'
#' @param x An expression.
#' @return A `sym` expression.
#' @export
default_toterm <- function(x) {
  x <- as_sym(x)
  # shifted variables lower to x_{t-k} (shift2term, discrete time module)
  if (.is_op(x, "Shift") && !is.null(.mtk$hooks$shift2term)) return(.mtk$hooks$shift2term(x))
  if (!.is_derivative(x)) return(x)
  di <- .deriv_info(x)
  .diff2term(di$var, di$iv, di$order)
}

# Unknown leaves appearing in equations (derivatives give their variable).
.eq_unknowns <- function(eqs) {
  vs <- .collect_vars(c(lapply(eqs, `[[`, "lhs"), lapply(eqs, `[[`, "rhs")))
  out <- list()
  for (v in vs) {
    if (.type(v) == "Derivative") v <- .args(v)[[1L]]
    if (identical(.var_kind(v), "unknown")) out[[length(out) + 1L]] <- v
  }
  .unique_syms(out)
}

# Move algebraic equations var ~ f(...) into observed when f depends only on
# differential variables, parameters and already observed variables.
# Mirrors get_trivial_observed_equations!.
.trivial_observed <- function(diffvars, algeqs, obseqs, all_dvs, iv) {
  dv_keys <- .keys(all_dvs)
  diff_keys <- .keys(diffvars)
  current_obs <- vapply(obseqs, function(e) .key(e$lhs), "")
  incid <- lapply(algeqs, function(e) {
    ks <- .keys(.eq_unknowns(list(Equation(.S0(), e$rhs))))
    setdiff(intersect(ks, dv_keys), c(diff_keys, current_obs))
  })
  keep <- rep(TRUE, length(algeqs))
  new_obs <- character()
  active <- TRUE
  iter <- 0L
  while (active && iter < 100L) {
    active <- FALSE
    iter <- iter + 1L
    for (i in seq_along(algeqs)) {
      if (!keep[i]) next
      lhs <- algeqs[[i]]$lhs
      lk <- .key(lhs)
      if (!(lk %in% dv_keys) || lk %in% new_obs || lk %in% diff_keys) next
      rest <- setdiff(incid[[i]], new_obs)
      if (length(rest)) next
      active <- TRUE
      new_obs <- c(new_obs, lk)
      obseqs[[length(obseqs) + 1L]] <- algeqs[[i]]
      keep[i] <- FALSE
    }
  }
  list(algeqs = algeqs[keep], obseqs = obseqs)
}

.base_compile <- function(sys, inputs = list(), outputs = list(), disturbance_inputs = list(),
                          fully_determined = TRUE) {
  eqs <- Filter(function(e) inherits(e, "Equation"), .f(sys, "eqs"))
  iv <- .f(sys, "iv")
  obs0 <- .f(sys, "observed")
  obs_keys <- vapply(obs0, function(e) .key(e$lhs), "")
  original <- Filter(function(v) !(.key(v) %in% obs_keys), .f(sys, "unknowns"))
  all_dvs <- Filter(function(v) .key(v) %in% .keys(original), .eq_unknowns(eqs))
  io_keys <- .keys(c(inputs, disturbance_inputs))
  unused <- setdiff(.keys(original), .keys(all_dvs))
  all_dvs <- Filter(function(v) !(.key(v) %in% io_keys), all_dvs)
  binds <- .f(sys, "bindings"); ics <- .f(sys, "initial_conditions")
  for (k in unused) {
    binds[[k]] <- NULL
    ics[[k]] <- NULL
  }
  if (fully_determined && length(eqs) > length(all_dvs))
    .stopf("The system is unbalanced. There are %d equations and %d unknowns.", length(eqs), length(all_dvs))
  if (fully_determined && length(eqs) < length(all_dvs))
    .stopf(paste("The system is unbalanced. There are %d equations and %d unknowns.",
                 "This may also be a high-index DAE; compile with structural = TRUE."),
           length(eqs), length(all_dvs))
  new_ps <- .unique_syms(c(.f(sys, "ps"), inputs, disturbance_inputs))
  has_derivs <- any(vapply(eqs, is_diff_equation, TRUE))
  if (!has_derivs && !is.null(h <- .mtk$hooks$compile_shift) && !is.null(r <- h(sys, eqs, all_dvs, obs0, new_ps, inputs, outputs, binds, ics))) return(r)

  if (!has_derivs) {
    tr <- .trivial_observed(list(), eqs, obs0, all_dvs, iv)
    alg <- lapply(tr$algeqs, function(e) Equation(.S0(), e$rhs - e$lhs))
    obs <- .toposort_observed(tr$obseqs)
    ob_keys <- vapply(obs, function(e) .key(e$lhs), "")
    dvs <- Filter(function(v) !(.key(v) %in% ob_keys), all_dvs)
    return(.update_system(sys, eqs = alg, unknowns = dvs, observed = obs, ps = new_ps,
                          inputs = inputs, outputs = outputs, isscheduled = TRUE,
                          bindings = binds, initial_conditions = ics,
                          schedule = list(var_sccs = list(seq_along(alg)))))
  }

  diff_idx <- vapply(eqs, function(e) .is_derivative(e$lhs), TRUE)
  diffeqs <- eqs[diff_idx]
  algeqs <- eqs[!diff_idx]
  # Order lowering on explicit (variable, order) pairs.
  sub_from <- list(); sub_to <- list()
  lowered <- list()
  diffvars <- list()
  for (eq in diffeqs) {
    di <- .deriv_info(eq$lhs)
    if (di$n_specs > 1L) stop("mixed partial derivative in ", .fmt(eq), call. = FALSE)
    var <- di$var
    terms <- lapply(0:(di$order - 1L), function(k) .diff2term(var, iv, k))
    for (k in seq_len(di$order - 1L)) {
      lowered[[length(lowered) + 1L]] <- Equation(.make_derivative(terms[[k]], iv, 1L), terms[[k + 1L]])
      diffvars[[length(diffvars) + 1L]] <- terms[[k]]
      sub_from[[length(sub_from) + 1L]] <- .make_derivative(var, iv, k)
      sub_to[[length(sub_to) + 1L]] <- terms[[k + 1L]]
    }
    last <- terms[[di$order]]
    lowered[[length(lowered) + 1L]] <- Equation(.make_derivative(last, iv, 1L), eq$rhs)
    diffvars[[length(diffvars) + 1L]] <- last
    sub_from[[length(sub_from) + 1L]] <- .make_derivative(var, iv, di$order)
    sub_to[[length(sub_to) + 1L]] <- eq$rhs
    if (di$order > 1L) {
      sub_from[[length(sub_from) + 1L]] <- .make_derivative(last, iv, 1L)
      sub_to[[length(sub_to) + 1L]] <- eq$rhs
    }
  }
  total <- .rule(sub_from, sub_to)
  diffeqs <- lapply(lowered, function(e) Equation(e$lhs, .subst_fixpoint(e$rhs, total)))
  tr <- .trivial_observed(diffvars, algeqs, obs0, all_dvs, iv)
  algeqs <- lapply(tr$algeqs, function(e) Equation(.S0(), .subst_fixpoint(e$rhs - e$lhs, total)))
  obseqs <- lapply(tr$obseqs, function(e) Equation(e$lhs, .subst_fixpoint(e$rhs, total)))
  obs_lhs <- vapply(obseqs, function(e) .key(e$lhs), "")
  alg_vars <- Filter(function(v) !(.key(v) %in% c(.keys(diffvars), obs_lhs, io_keys)), all_dvs)
  new_eqs <- c(diffeqs, algeqs)
  for (e in new_eqs) {
    if (.has_derivative(e$rhs))
      stop("derivative in the right-hand side of ", .fmt(e),
           "; compile with structural = TRUE", call. = FALSE)
  }
  .update_system(sys, eqs = new_eqs, observed = .toposort_observed(obseqs),
                 unknowns = c(diffvars, alg_vars), ps = new_ps, inputs = inputs,
                 outputs = outputs, isscheduled = TRUE, bindings = binds,
                 initial_conditions = ics,
                 schedule = list(var_sccs = list(seq_along(new_eqs)), dummy_sub = total))
}
