# Initialization problems and their use by the problem constructors, after
# ModelingToolkitBase/src/problems/initializationproblem.jl
# (InitializationProblem, get_initialization_problem_type, the messages for
# underdetermined, overdetermined and singular initialization),
# systems/problem_utils.jl (maybe_build_initialization_problem,
# InitializationMetadata, MissingGuessValue), systems/nonlinear/
# initializesystem.jl (remake_initialization_data, late_binding_update_u0_p)
# and ModelingToolkit/src/initialization.jl (analyze_initialization_jacobian).
#
# In ModelingToolkit the initialization problem is stored with the ODE
# function and solved by the integrator. Here it is solved when the problem
# is built or remade: the hook `initialize_problem` fills prob$u0 and the
# solvable parameters in prob$p, and keeps everything needed to run again in
# prob$initialization_data.

# ---------------------------------------------------------------------------
# Compiling the initialization system

# An opaque parameter symbol standing for Initial(x). Inside the compiled
# initialization system Initial(x) must not look like an occurrence of x.
# It is named by the symbolic key of Initial(x), not by its printed form, so
# that a parameter literally named "Initial(x(t))" cannot collide with it.
.init_opaque_symbol <- function(ix) {
  nm <- .key(ix)
  if (!exists(nm, envir = .mtk$vars, inherits = FALSE)) .register_var(nm, "parameter")
  s <- Symbol(nm)
  attr(s, "mtk_meta") <- list(kind = "parameter")
  s
}

# The counterpart of mtkcompile(isys; fully_determined = false): explicit
# equations whose right-hand side is known become observed equations, and so
# do unknowns that appear in a single equation solved for them (such as
# unused derivative variables). What remains is in residual form.
.init_reduce <- function(isys) {
  eqs <- Filter(function(e) inherits(e, "Equation"), .f(isys, "eqs"))
  us <- .f(isys, "unknowns")
  ps <- .f(isys, "ps")
  inodes <- Filter(isinitial, ps)
  isyms <- lapply(inodes, .init_opaque_symbol)
  orule <- if (length(inodes)) .rule(inodes, isyms) else NULL
  sub <- function(e) if (is.null(orule) || !is_sym(e)) e else .subst(e, orule)
  eqs <- lapply(eqs, function(e) Equation(sub(e$lhs), sub(e$rhs)))
  pars <- c(Filter(function(p) !isinitial(p), ps), isyms)
  remap <- function(m) {
    out <- symmap()
    ik <- .keys(inodes)
    for (el in unclass(m)) {
      if (is.null(el$var)) next
      j <- match(.key(el$var), ik)
      v <- if (is.na(j)) el$var else isyms[[j]]
      out[[.key(v)]] <- list(var = v, value = sub(el$value))
    }
    class(out) <- "mtk_symmap"
    out
  }
  uk <- .keys(us)
  n <- length(eqs)
  inc_l <- lapply(eqs, function(e) intersect(.nl_atomic_keys(e$lhs), uk))
  inc_r <- lapply(eqs, function(e) intersect(.nl_atomic_keys(e$rhs), uk))
  lk <- vapply(eqs, function(e) .key(e$lhs), "")
  rk <- vapply(eqs, function(e) .key(e$rhs), "")
  alive <- rep(TRUE, n)
  solved <- character()
  obs <- list()
  # explicit assignments from parameters and already solved variables
  repeat {
    progressed <- FALSE
    for (i in which(alive)) {
      if (lk[i] %in% uk && !(lk[i] %in% solved) && !length(setdiff(inc_r[[i]], solved))) {
        obs[[length(obs) + 1L]] <- eqs[[i]]
        solved <- c(solved, lk[i])
      } else if (rk[i] %in% uk && !(rk[i] %in% solved) && !length(setdiff(inc_l[[i]], solved))) {
        obs[[length(obs) + 1L]] <- Equation(eqs[[i]]$rhs, eqs[[i]]$lhs)
        solved <- c(solved, rk[i])
      } else next
      alive[i] <- FALSE
      progressed <- TRUE
    }
    if (!progressed) break
  }
  # unknowns that occur in one equation only, which is solved for them
  repeat {
    progressed <- FALSE
    live <- which(alive)
    if (!length(live)) break
    occ <- table(unlist(lapply(live, function(i) setdiff(union(inc_l[[i]], inc_r[[i]]), solved))))
    single <- names(occ)[occ == 1L]
    for (i in live) {
      if (lk[i] %in% single && !(lk[i] %in% inc_r[[i]])) {
        obs[[length(obs) + 1L]] <- eqs[[i]]
        solved <- c(solved, lk[i])
      } else if (rk[i] %in% single && !(rk[i] %in% inc_l[[i]])) {
        obs[[length(obs) + 1L]] <- Equation(eqs[[i]]$rhs, eqs[[i]]$lhs)
        solved <- c(solved, rk[i])
      } else next
      alive[i] <- FALSE
      progressed <- TRUE
    }
    if (!progressed) break
  }
  remaining <- us[!(uk %in% solved)]
  resid <- lapply(eqs[alive], function(e) Equation(.S0(), e$rhs - e$lhs))
  out <- System(resid, NULL, unknowns = remaining, ps = pars, observed = obs,
                bindings = remap(.f(isys, "bindings")), initial_conditions = remap(.f(isys, "initial_conditions")),
                guesses = .f(isys, "guesses"), name = nameof(isys), checks = FALSE,
                discover_from_metadata = FALSE)
  md <- .f(isys, "metadata")
  md$initialization$initial_nodes <- inodes
  md$initialization$initial_symbols <- isyms
  out <- .update_system(out, is_initializesystem = TRUE, isscheduled = TRUE, metadata = md,
                        complete = TRUE, namespacing = FALSE)
  .update_system(out, index_cache = .build_index_cache(out))
}

# ---------------------------------------------------------------------------
# Messages

.init_pluralize <- function(n, word) if (n == 1L) word else paste0(word, "s")

.init_overdetermined_message <- function(neqs, nunknown) {
  paste0("Initialization system is overdetermined. ", neqs, " equations for ", nunknown,
         " unknowns. Initialization will default to using least squares.\n\n",
         "Call `analyze_initialization_jacobian(prob)` on the constructed problem to see which ",
         "equations are redundant (and which unknowns, if any, remain underdetermined).\n\n",
         "To suppress this warning, pass `warn_initialize_determined = FALSE`. To turn this warning ",
         "into an error, pass `fully_determined = TRUE`.")
}

.init_underdetermined_message <- function(neqs, nunknown) {
  paste0("Initialization system is underdetermined. ", neqs, " equations for ", nunknown,
         " unknowns. Initialization will default to using least squares.\n\n",
         "Call `analyze_initialization_jacobian(prob)` on the constructed problem to see which ",
         "unknowns are underdetermined (and which equations, if any, are redundant).\n\n",
         "To suppress this warning, pass `warn_initialize_determined = FALSE`. To turn this warning ",
         "into an error, pass `fully_determined = TRUE`.")
}

.init_singular_message <- function(vars) {
  paste0("The initialization system is structurally singular. Guess values may significantly affect ",
         "the initial values of the ODE. The problematic variables are ",
         paste(.nl_names(vars), collapse = ", "), ".\n\n",
         "Note that the identification of problematic variables is a best-effort heuristic.\n\n",
         "To suppress this warning, pass `warn_initialize_determined = FALSE`.")
}

# initialization_structure_error_message: the error raised with
# fully_determined = TRUE.
.init_structure_error <- function(kind, neqs, nunknown, vars = list()) {
  context <- paste0(
    "This is an error in the initialization system, not in the equations being integrated. The ",
    "initialization system solves for the value of every unknown of the model, and of their ",
    "derivatives, at the initial time. It needs as many equations as there are such values to solve ",
    "for. Initial values given for the model become equations of this system.\n\n")
  count <- paste0("The initialization system has ", neqs, " ", .init_pluralize(neqs, "equation"), " for ",
                  nunknown, " ", .init_pluralize(nunknown, "unknown"), ".\n\n")
  notation <- "Note that variables named with a `\u02CDt` suffix are derivatives at the initial time: `x\u02CDt` is `D(x)`."
  msg <- switch(kind,
    under = {
      d <- nunknown - neqs
      paste0("Initialization system is underdetermined.\n\n", context, count, d, " more ",
             .init_pluralize(d, "equation"), if (d == 1L) " is" else " are",
             " needed to determine the initial state.\n\n",
             "Each missing equation is supplied by giving one more initial value, or one more equation ",
             "relating initial values. Guesses are starting values for the initialization solve rather ",
             "than constraints on it, so adding a guess does not supply one.\n\n",
             "Any of the following supplies one missing equation:\n\n",
             "  * Give an initial value in the problem constructor, as in `ODEProblem(sys, list(x ~ 1), tspan)`.\n",
             "  * Give one in the model, via `initial_conditions` or `bindings`.\n",
             "  * Add an equation relating initial values, via the `initialization_eqs` argument of the ",
             "system or of the problem constructor.\n\n",
             "To solve an underdetermined initialization in a least squares sense instead, pass ",
             "`fully_determined = FALSE` to the problem constructor; the guesses then decide which of ",
             "the possible initial states is found.\n\n", notation)
    },
    over = {
      d <- neqs - nunknown
      paste0("Initialization system is overdetermined.\n\n", context, count, "There ",
             if (d == 1L) "is" else "are", " ", d, " ", .init_pluralize(d, "equation"), " too many.\n\n",
             "Remove initial values or equations relating them until as many remain as there are ",
             "unknowns to solve for.\n\n",
             "These come from `initial_conditions`, `bindings` and `initialization_eqs`, and from the ",
             "operating point passed to the problem constructor. To solve an overdetermined ",
             "initialization in a least squares sense instead, pass `fully_determined = FALSE` to the ",
             "problem constructor. That only finds an initial state if the extra equations are ",
             "consistent with the rest.\n\n", notation)
    },
    singular = paste0(
      "Initialization system is structurally singular.\n\n", context, count,
      "No equation determines ", paste(.nl_names(vars), collapse = ", "), ".\n\n",
      "There are as many equations as unknowns, but they do not determine a unique initial state: at ",
      "least one equation is redundant given the others, leaving some unknown undetermined. This ",
      "usually means two of the given initial conditions are related by an equation of the model, so ",
      "that one of them carries no new information. Giving an initial condition for a different ",
      "variable, or for a derivative, in place of the redundant one resolves it.\n\n", notation))
  stop(msg, call. = FALSE)
}

# ---------------------------------------------------------------------------
# Building the problem

.init_default_mgv <- function(mgv) if (is.null(mgv)) MissingGuessValue$HashedRandom() else mgv

# Build the compiled initialization system and its nonlinear problem.
# Returns the data kept in prob$initialization_data.
.init_build <- function(sys, t0, op, guesses = NULL, initialization_eqs = NULL, fully_determined = NULL,
                        warn = TRUE, missing_guess_value = NULL, algebraic_only = FALSE,
                        time_dependent_init = TRUE, jac = FALSE, extra_values = NULL, use_scc = TRUE,
                        use_homotopy_continuation = FALSE) {
  isys <- generate_initializesystem(sys, op = op, initialization_eqs = initialization_eqs %||% list(),
                                    guesses = guesses %||% list(), algebraic_only = algebraic_only,
                                    time_dependent_init = time_dependent_init)
  csys <- .init_reduce(isys)
  us <- unknowns(csys)
  neqs <- length(equations(csys)); nunknown <- length(us)
  unmatched <- list()
  if (neqs && nunknown) {
    m <- .nl_matching(.nl_incidence(.ode_rhs_exprs(csys), us), nunknown)
    unmatched <- us[m$eq_of_var == 0L]
  } else if (nunknown) unmatched <- us
  singular <- neqs == nunknown && length(unmatched) > 0L
  if (isTRUE(fully_determined)) {
    if (neqs < nunknown) .init_structure_error("under", neqs, nunknown)
    if (neqs > nunknown) .init_structure_error("over", neqs, nunknown)
    if (singular) .init_structure_error("singular", neqs, nunknown, unmatched)
  } else if (isTRUE(warn)) {
    if (neqs > nunknown) warning(.init_overdetermined_message(neqs, nunknown), call. = FALSE)
    if (neqs < nunknown) warning(.init_underdetermined_message(neqs, nunknown), call. = FALSE)
    if (singular) warning(.init_singular_message(unmatched), call. = FALSE)
  }
  square <- neqs == nunknown && !singular
  Ab <- if (nunknown > 0L && nunknown <= neqs) .nl_calculate_A_b(csys, throw = FALSE)
  type <- .init_problem_type(csys, neqs, nunknown, square, !is.null(Ab), use_scc, use_homotopy_continuation)
  fn <- NonlinearFunction(csys, jac = jac)
  vals <- .process_op(csys, extra_values %||% list(), require_p = FALSE)
  p <- vals$p
  iv <- .f(sys, "iv")
  t_index <- if (is.null(iv)) NA_integer_ else match(format(iv), names(p))
  if (!is.na(t_index)) p[t_index] <- t0 %||% 0
  if (anyNA(p))
    stop("no value for parameter(s) ", paste(names(p)[is.na(p)], collapse = ", "),
         " of the initialization system; give them in the operating point or as defaults", call. = FALSE)
  u0 <- vals$u0
  if (anyNA(u0)) {
    miss <- which(is.na(u0))
    gv <- .nl_eval_guesses(guesses(csys), us[miss], c(us[-miss], .problem_params(csys)), c(u0[-miss], p),
                           .codegen_observed(csys))
    if (anyNA(gv))
      gv[is.na(gv)] <- .nl_missing_values(us[miss][is.na(gv)], .init_default_mgv(missing_guess_value), TRUE)
    u0[miss] <- gv
  }
  if (type == "HomotopyContinuationProblem") {
    hc <- HomotopyContinuationProblem(csys, as.list(c(u0, p)))
    fn$denominator <- hc$f$denominator
    fn$numerator <- hc$f
  }
  iprob <- .new_problem(type, f = fn, u0 = u0, p = p, sys = csys, kwargs = list())
  iprob$op <- vals$op
  iprob$is_initializeprob <- TRUE
  class(iprob) <- c(type, "AbstractNonlinearProblem", "SciMLProblem")
  if (type == "LinearProblem") {
    # LinearInitializationProblem: A u = b with A and b functions of the
    # parameters; least squares when overdetermined
    pargs <- list(list(name = "p", vars = .problem_params(csys)))
    obs <- .codegen_observed(csys)
    iprob$f$update_A <- .build_closure(Ab$A, pargs, observed = obs, out = "matrix")
    iprob$f$update_b <- .build_closure(Ab$b, pargs, observed = obs, out = "vector")
    iprob$A <- iprob$f$update_A(p)
    iprob$b <- as.numeric(iprob$f$update_b(p))
  } else if (type == "SCCNonlinearProblem") {
    iprob <- .nl_as_scc(iprob)
  } else if (type == "HomotopyProblem") {
    iprob$f$homotopy <- .nl_homotopy_residual(csys)
    iprob$lambdaspan <- c(0, 1)
  } else if (type == "HomotopyContinuationProblem") {
    iprob$degrees <- hc$degrees
  }
  list(initializeprob = iprob, initializesystem = isys, u0_template = u0, p_template = p,
       t_index = t_index, square = square)
}

# get_initialization_problem_type (ModelingToolkit/src/initialization.jl):
# a linear problem for an affine system with at least as many equations as
# unknowns, an SCC problem (or with use_scc = FALSE a homotopy or plain
# nonlinear problem) for a square nonsingular one, least squares otherwise.
# use_homotopy_continuation = TRUE asks for polynomial homotopy continuation
# of a square system instead (an extension of this package).
.init_problem_type <- function(csys, neqs, nunknown, square, affine, use_scc = TRUE,
                               use_homotopy_continuation = FALSE) {
  if (isTRUE(use_homotopy_continuation) && square && nunknown > 0L) return("HomotopyContinuationProblem")
  if (nunknown > 0L && nunknown <= neqs && affine) return("LinearProblem")
  if (square) {
    if (!isFALSE(use_scc) && neqs > 0L) return("SCCNonlinearProblem")
    return(if (.nl_has_call(csys, "homotopy")) "HomotopyProblem" else "NonlinearProblem")
  }
  "NonlinearLeastSquaresProblem"
}

#' Initialization problems
#'
#' `InitializationProblem(sys, t, op)` builds the nonlinear problem whose
#' solution is a consistent set of initial values for `sys` at time `t`:
#' the initialization system of [generate_initializesystem()] is simplified
#' (explicit equations become observed) and wrapped in the problem type
#' ModelingToolkit chooses: a [LinearProblem()] when the equations are
#' affine in the unknowns and there are at least as many equations as
#' unknowns (solved in the least squares sense when overdetermined), an
#' [SCCNonlinearProblem()] when it is square and structurally nonsingular
#' (with `use_scc = FALSE` a [HomotopyProblem()] if it contains [homotopy()]
#' nodes and a [NonlinearProblem()] otherwise), and a
#' [NonlinearLeastSquaresProblem()] for under and overdetermined systems.
#' `use_homotopy_continuation = TRUE` asks for a
#' [HomotopyContinuationProblem()] of a square polynomial system instead,
#' which finds every real root and keeps the one closest to the guesses.
#' Underdetermined, overdetermined and singular systems give a warning, or an
#' error with `fully_determined = TRUE`.
#'
#' The same choices apply to the initialization that [ODEProblem()] and the
#' other problem constructors run: pass `use_scc` and
#' `use_homotopy_continuation` to them.
#'
#' [ODEProblem()] builds and solves this problem itself; this constructor is
#' for inspecting initialization. Solve the result with `solve()`, and index
#' the solution with the variables of `sys`.
#'
#' @param sys A complete `System`.
#' @param t The initial time; `NULL` for time independent systems.
#' @param op Operating point: initial conditions and parameter values.
#' @param guesses Further guesses.
#' @param initialization_eqs Further initialization equations.
#' @param fully_determined `TRUE` turns the warnings about systems that are
#'   not fully determined into errors.
#' @param warn_initialize_determined `FALSE` silences those warnings.
#' @param missing_guess_value Fallback for unknowns without a guess; the
#'   default is `MissingGuessValue$HashedRandom()`.
#' @param algebraic_only Ignore the initialization equations of the system.
#' @param time_dependent_init Solve for the unknowns too (the default for
#'   time dependent systems), or only for solvable parameters.
#' @param jac Generate the symbolic Jacobian.
#' @param use_scc Solve a square system block by block.
#' @param use_homotopy_continuation Solve a square polynomial system by
#'   homotopy continuation.
#' @param ... Solver options stored with the problem.
#' @return A `LinearProblem`, `SCCNonlinearProblem`, `NonlinearProblem`,
#'   `HomotopyProblem`, `HomotopyContinuationProblem` or
#'   `NonlinearLeastSquaresProblem` whose system is the compiled
#'   initialization system.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); y <- Variable("y", t)
#' sys <- mtkcompile(System(list(D(x) ~ -x + y, 0 ~ x^2 + y^2 - 1), t, name = "dae"))
#' iprob <- InitializationProblem(sys, 0, list(x ~ 0.6), guesses = list(y ~ 1))
#' class(iprob)[1]
#' sol <- solve(iprob)
#' sol[y]
#' @export
InitializationProblem <- function(sys, t = NULL, op = list(), guesses = list(), initialization_eqs = list(),
                                  fully_determined = NULL, warn_initialize_determined = TRUE,
                                  missing_guess_value = NULL, algebraic_only = FALSE,
                                  time_dependent_init = !is.null(get_iv(sys)), jac = TRUE, use_scc = TRUE,
                                  use_homotopy_continuation = FALSE, ...) {
  .need_complete(sys, "InitializationProblem")
  b <- .init_build(sys, t %||% 0, op, guesses, initialization_eqs, fully_determined,
                   warn_initialize_determined, missing_guess_value, algebraic_only,
                   time_dependent_init, jac = jac, use_scc = use_scc,
                   use_homotopy_continuation = use_homotopy_continuation)
  prob <- b$initializeprob
  prob$kwargs <- list(...)
  prob$initializesystem <- b$initializesystem
  prob
}

# ---------------------------------------------------------------------------
# The initialize_problem hook

# Solve an initialization problem with the solver of its type and judge it:
# least squares solutions count only when the residual vanishes, since the
# equations are constraints.
.init_solve <- function(iprob, opts = list()) {
  fns <- .nl_problem_fns(iprob)
  o <- utils::modifyList(list(abstol = 1e-10, maxiters = 1000L), opts)
  alg <- o$alg
  o$alg <- NULL
  if (inherits(iprob, "LinearProblem")) {
    A <- iprob$f$update_A(iprob$p)
    if (!is.matrix(A)) A <- matrix(A, ncol = length(iprob$u0))
    u <- .nl_lstsq(A, as.numeric(iprob$f$update_b(iprob$p)))
    res <- .nl_result(u, NULL, "Success", 1L)
  } else if (inherits(iprob, "SCCNonlinearProblem")) {
    res <- .nl_scc_run(iprob, alg, o)
  } else if (inherits(iprob, "HomotopyProblem")) {
    hf <- iprob$f$homotopy
    p <- iprob$p
    at <- function(lam) list(F = function(x) as.numeric(hf$f(x, p, lam)), J = function(x) hf$jac(x, p, lam))
    res <- .nl_sweep(at, iprob$u0, if (inherits(alg, "HomotopySweep")) alg else HomotopySweep(inner = alg), o, 0, 1)
  } else if (inherits(iprob, "HomotopyContinuationProblem")) {
    hp <- iprob
    hp$f <- iprob$f$numerator
    hp$f$denominator <- iprob$f$denominator
    sol <- solve(hp, alg)
    res <- .nl_result(as.numeric(sol$u), NULL, sol$retcode, 0L)
  } else {
    ls <- inherits(iprob, "NonlinearLeastSquaresProblem")
    alg <- alg %||% if (ls) FastShortcutNLLSPolyalg() else .nl_default_alg(iprob$u0, iprob$f$jac)
    res <- .nl_run(alg, fns$F, fns$J, iprob$u0, o, ls = ls)
  }
  res$resid <- fns$F(res$u)
  resid <- res$resid
  scale <- max(1, abs(res$u), abs(iprob$u0))
  res$ok <- all(is.finite(resid)) && max(abs(resid), 0) <= 1e-8 * scale
  res
}

# The conditions of a problem: the signature below and the constants.
.init_digest <- function(op, kw) {
  num <- vapply(unclass(op), function(el) {
    v <- el$value
    if (is_sym(v) && .is_num(v)) v <- as.numeric(v)
    if (is.numeric(v) && length(v) == 1L) v else NA_real_
  }, 0)
  list(.init_signature(op, kw), num)
}

# What changes the equations of the initialization system: the constrained
# keys and whether each value is a constant, an expression or missing.
.init_signature <- function(op, kw) {
  parts <- vapply(names(op), function(k) {
    v <- op[[k]]$value
    paste0(k, "=", if (is.null(v)) "NULL" else if (is_sym(v)) (if (.is_num(v)) "#" else .key(v))
           else if (length(v) == 1L && is.na(v)) "NA" else "#")
  }, "")
  ie <- vapply(.as_equations(kw$initialization_eqs), function(e) paste(.key(e$lhs), .key(e$rhs)), "")
  paste(c(sort(parts), "|", ie, "|", isTRUE(kw$fully_determined)), collapse = ";")
}

.init_supported <- function(prob) {
  !inherits(prob, c("DiscreteProblem", "ImplicitDiscreteProblem", "JumpProblem", "BVProblem",
                    "DAEProblem", "DDEProblem", "SDDEProblem"))
}

.init_hook <- function(prob, vals = NULL, check = TRUE) {
  sys <- prob$sys
  if (!inherits(sys, "System")) return(prob)
  if (is.null(.f(sys, "iv"))) return(.init_static_hook(prob))
  u0 <- prob$u0
  if (!isTRUE(check) || !.init_supported(prob) || is.null(u0) || !is.numeric(u0)) return(prob)
  p <- prob$p
  data <- prob$initialization_data
  t0 <- prob$tspan[1L] %||% 0
  kw <- prob$kwargs %||% list()
  op <- prob$op %||% symmap()
  # nothing changed since the last run: same values, time and conditions
  if (!is.null(data$u0_out) && identical(u0, data$u0_out) && identical(p, data$p_out) &&
      identical(t0, data$t0) && identical(.init_digest(op, kw), data$digest)) return(prob)
  ic <- .f(sys, "index_cache")
  known <- if (is.null(ic)) .keys(c(unknowns(sys), .parameters_of(sys, TRUE)))
           else c(names(ic$unknown_index), names(ic$param_index))
  okeys <- if (is.null(ic)) .keys(observables(sys)) else names(ic$observed_index)
  # values of observed variables and of derivatives, and bindings of
  # observed variables, are conditions only the initialization system can
  # impose
  ukeys <- if (is.null(ic)) .keys(unknowns(sys)) else names(ic$unknown_index)
  extra <- !all(names(op) %in% known) || any(names(.f(sys, "initial_conditions")) %in% okeys) ||
    any(names(.f(sys, "bindings")) %in% c(okeys, ukeys))
  ieqs <- length(.f(sys, "initialization_eqs")) > 0L || length(kw$initialization_eqs) > 0L
  if (!(anyNA(u0) || anyNA(p) || ieqs || extra || !is.null(data$initializeprob))) {
    fn <- prob$f
    M <- if (is.list(fn)) fn$mass_matrix else NULL
    alg_rows <- if (is.matrix(M)) which(rowSums(abs(M)) == 0) else integer()
    if (!length(alg_rows)) return(prob)
    r <- tryCatch(fn$f(u0, p, t0)[alg_rows], error = function(e) NA_real_)
    if (all(is.finite(r)) && max(abs(r)) <= 1e-8 * max(1, abs(u0))) return(prob)
  }
  .init_run(prob, data, t0, kw, op)
}

.init_run <- function(prob, data, t0, kw, op) {
  sys <- prob$sys
  u0 <- prob$u0; p <- prob$p
  sig <- .init_signature(op, kw)
  fresh <- is.null(data$initializeprob) || !identical(data$signature, sig)
  if (fresh) {
    # Numeric values of the unknowns (from the operating point, or given to
    # remake() as a vector) are the constants of the constant conditions,
    # and starting values for the unknowns that are solved for.
    us <- unknowns(sys)
    uk <- .keys(us)
    ics <- .f(sys, "initial_conditions")
    op_b <- op
    gs_b <- .resolve_named(as_symmap(kw$guesses), c(us, .parameters_of(sys, TRUE), observables(sys)))
    for (i in which(!is.na(u0))) {
      el <- op_b[[uk[i]]] %||% ics[[uk[i]]]
      v <- el$value
      if (is_sym(v) && .is_num(v)) v <- as.numeric(v)
      if (is.numeric(v) && length(v) == 1L && !is.na(v)) {
        if (v != u0[[i]]) op_b <- smap_set(op_b, us[[i]], u0[[i]])
      } else if (is.null(el)) gs_b <- smap_set(gs_b, us[[i]], u0[[i]])
    }
    data <- .init_build(sys, t0, op_b, gs_b, kw$initialization_eqs, kw$fully_determined,
                        kw$warn_initialize_determined %||% TRUE, kw$missing_guess_value,
                        use_scc = kw$use_scc %||% TRUE,
                        use_homotopy_continuation = isTRUE(kw$use_homotopy_continuation))
    data$signature <- sig
    data <- .init_maps(data, sys, prob)
  }
  iprob <- data$initializeprob
  ip <- if (fresh || is.null(data$ip_last)) data$p_template else data$ip_last
  iu <- data$u0_template
  # parameters of the system, by name
  common <- intersect(names(ip), names(p))
  common <- common[!is.na(p[common])]
  ip[common] <- p[common]
  if (!is.na(data$t_index)) ip[data$t_index] <- t0
  opnum <- function(k) {
    v <- op[[k]]$value
    if (is_sym(v) && .is_num(v)) v <- as.numeric(v)
    if (is.numeric(v) && length(v) == 1L) v else NA_real_
  }
  consts <- vapply(data$initial_key, opnum, 0, USE.NAMES = FALSE)
  if (!fresh) {
    # constants of the conditions x ~ Initial(x): new values in the
    # operating point, and for unknowns new numeric values given to remake()
    changed <- !identical(u0, data$u0_out)
    for (j in seq_along(data$initial_p)) {
      if (!is.na(consts[j]) && !identical(consts[j], data$op_consts[j])) ip[data$initial_p[j]] <- consts[j]
      ui <- data$initial_u[j]
      if (changed && !is.na(ui) && !is.na(u0[ui])) ip[data$initial_p[j]] <- u0[ui]
    }
    # starting values for the unknowns that are solved for
    if (changed) for (j in seq_along(data$guess_i)) {
      v <- u0[data$guess_u[j]]
      if (!is.na(v)) iu[data$guess_i[j]] <- v
    }
  }
  iprob$p <- ip
  iprob$u0 <- iu
  res <- .init_solve(iprob, kw$initialization_options %||% list())
  unew <- data$umap(res$u, ip)
  names(unew) <- names(u0)
  if (length(data$solvable_names)) p[data$solvable_names] <- data$pmap(res$u, ip)
  if (length(data$sys_initial_names)) p[data$sys_initial_names] <- data$imap(res$u, ip)
  retcode <- if (isTRUE(res$ok)) "Success" else "InitialFailure"
  if (!isTRUE(res$ok))
    warning("Initialization failed (retcode InitialFailure): the initial conditions could not be ",
            "satisfied; the largest residual of the initialization system is ",
            format(max(abs(res$resid), 0), digits = 3), ". The initial conditions may be inconsistent. ",
            "Call `analyze_initialization_jacobian(prob)` to inspect the initialization system.",
            call. = FALSE)
  iprob$u0 <- stats::setNames(res$u, names(iu))
  if (inherits(iprob, "LinearProblem")) {
    iprob$A <- iprob$f$update_A(ip)
    iprob$b <- as.numeric(iprob$f$update_b(ip))
  }
  data$initializeprob <- iprob
  data$retcode <- retcode
  data$resid <- res$resid
  data$u0_out <- unew
  data$p_out <- p
  data$t0 <- t0
  data$ip_last <- ip
  data$op_consts <- consts
  data$digest <- .init_digest(op, kw)
  prob$u0 <- unew
  prob$p <- p
  prob$initialization_data <- data
  # solve() on a problem whose initialization failed returns the retcode
  # InitialFailure instead of integrating; see solve.InitializationFailure()
  cls <- setdiff(class(prob), "InitializationFailure")
  class(prob) <- if (isTRUE(res$ok)) cls else c("InitializationFailure", cls)
  prob
}

# Closures and index maps between the problem and its initialization
# problem.
.init_maps <- function(data, sys, prob) {
  csys <- data$initializeprob$sys
  args <- .codegen_args(csys, with_t = FALSE)
  obs <- .codegen_observed(csys)
  us <- unknowns(sys)
  data$umap <- .build_closure(us, args, observed = obs, out = "vector")
  info <- .f(csys, "metadata")$initialization
  pnames <- names(prob$p)
  ps <- .problem_params(sys)
  solv <- which(.keys(ps) %in% info$solvable)
  data$solvable_names <- pnames[solv]
  data$pmap <- if (length(solv)) .build_closure(ps[solv], args, observed = obs, out = "vector") else NULL
  # Initial(x) parameters of the system take the initial value of x
  sinit <- which(vapply(ps, isinitial, TRUE))
  data$sys_initial_names <- pnames[sinit]
  data$imap <- NULL
  if (length(sinit)) {
    src <- lapply(ps[sinit], .init_source)
    known <- .keys(c(unknowns(csys), observables(csys), .problem_params(csys)))
    obs_sys <- observed(sys)
    rule <- if (length(obs_sys)) .rule(lapply(obs_sys, `[[`, "lhs"), lapply(obs_sys, `[[`, "rhs")) else NULL
    src <- lapply(src, function(s) {
      if (!(.key(s) %in% known) && !is.null(rule)) s <- .subst_fixpoint(s, rule)
      # a source the initialization does not determine keeps the value zero
      if (all(.nl_atomic_keys(s) %in% known)) s else .S0()
    })
    data$imap <- .build_closure(src, args, observed = obs, out = "vector")
  }
  # positions for numeric updates
  ipn <- names(data$p_template)
  ukeys <- .keys(us)
  nodes <- info$initial_nodes; syms <- info$initial_symbols
  initial_p <- integer(); initial_u <- integer(); initial_key <- character()
  for (j in seq_along(nodes)) {
    pi <- match(format(syms[[j]]), ipn)
    if (is.na(pi)) next
    sk <- .key(.init_source(nodes[[j]]))
    initial_p <- c(initial_p, pi)
    initial_u <- c(initial_u, match(sk, ukeys))
    initial_key <- c(initial_key, sk)
  }
  data$initial_p <- initial_p; data$initial_u <- initial_u; data$initial_key <- initial_key
  gi <- match(.keys(unknowns(csys)), ukeys)
  data$guess_i <- which(!is.na(gi)); data$guess_u <- gi[!is.na(gi)]
  data
}

# Time independent problems on remake: unknowns without a value take their
# guesses again, and linear problems refresh A and b.
.init_static_hook <- function(prob) {
  sys <- prob$sys
  if (inherits(prob, "LinearProblem")) {
    if (!is.null(prob$f$update_A)) {
      prob$A <- prob$f$update_A(prob$p)
      prob$b <- as.numeric(prob$f$update_b(prob$p))
    }
    return(prob)
  }
  if (!is.null(prob$p) && anyNA(prob$p) && is.numeric(prob$u0)) {
    pinit <- .init_parameters_hook(sys, list(u0 = prob$u0, p = prob$p, op = prob$op), prob$kwargs %||% list())
    prob$p <- pinit$p
    prob$initialization_data <- pinit$data
  }
  if (is.numeric(prob$u0) && anyNA(prob$u0))
    prob$u0 <- .nl_fill_guesses(sys, prob$u0, prob$p, prob$guesses, prob$missing_guess_value)
  .nl_refresh_bounds(prob)
}

# Time independent initialization: solve for the solvable parameters of a
# nonlinear problem.
.init_parameters_hook <- function(sys, vals, dots = list()) {
  p <- vals$p
  ps <- .problem_params(sys)
  binds <- bindings(sys)
  us <- unknowns(sys)
  # Initial(x) parameters of a time independent system are the starting
  # values of its unknowns
  for (i in which(is.na(p) & vapply(ps, isinitial, TRUE))) {
    j <- match(.key(.init_source(ps[[i]])), .keys(us))
    p[i] <- if (!is.na(j) && !is.na(vals$u0[j])) vals$u0[j] else 0
  }
  if (!anyNA(p)) return(list(p = p, data = NULL))
  miss <- which(is.na(p))
  solvable <- vapply(ps[miss], function(q) {
    .init_is_missing(smap_get(binds, q)) || .init_is_missing(vals$op[[.key(q)]]$value)
  }, TRUE)
  if (any(!solvable))
    stop("no value for parameter(s) ", paste(names(p)[miss][!solvable], collapse = ", "),
         "; give them in the operating point or as defaults", call. = FALSE)
  op <- vals$op %||% symmap()
  for (i in seq_along(us))
    if (!is.na(vals$u0[i]) && is.null(op[[.key(us[[i]])]])) op <- smap_set(op, us[[i]], vals$u0[[i]])
  data <- .init_build(sys, NULL, op, dots$guesses, dots$initialization_eqs, dots$fully_determined,
                      dots$warn_initialize_determined %||% TRUE, dots$missing_guess_value,
                      time_dependent_init = FALSE, use_scc = dots$use_scc %||% TRUE,
                      use_homotopy_continuation = isTRUE(dots$use_homotopy_continuation))
  iprob <- data$initializeprob
  res <- .init_solve(iprob)
  if (!isTRUE(res$ok))
    warning("Initialization of the parameters failed: the largest residual of the initialization system is ",
            format(max(abs(res$resid), 0), digits = 3), ".", call. = FALSE)
  csys <- iprob$sys
  pmap <- .build_closure(ps[miss], .codegen_args(csys, with_t = FALSE), observed = .codegen_observed(csys),
                         out = "vector")
  p[miss] <- pmap(res$u, iprob$p)
  data$retcode <- if (isTRUE(res$ok)) "Success" else "InitialFailure"
  data$resid <- res$resid
  list(p = p, data = data)
}

#' Problems whose initialization failed
#'
#' When the initial conditions of a problem cannot be satisfied, the problem
#' constructor warns and marks the problem with the class
#' `InitializationFailure`. `solve()` on such a problem does not integrate:
#' it returns a solution holding the initial point with the return code
#' `"InitialFailure"`, as the integrators of the SciML ecosystem do. Fix the
#' initial conditions (see [analyze_initialization_jacobian()]) and build or
#' [remake()] the problem again.
#'
#' @param a A problem marked `InitializationFailure`.
#' @param b An algorithm; unused.
#' @param ... Unused.
#' @return A solution with `sol$retcode` equal to `"InitialFailure"`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); y <- Variable("y", t)
#' sys <- mtkcompile(System(list(D(x) ~ -x + y, 0 ~ x^2 + y^2 - 1), t, name = "dae"))
#' prob <- suppressWarnings(ODEProblem(sys, list(x ~ 0.6, y ~ 0.5), c(0, 1)))
#' solve(prob)$retcode
#' @export
solve.InitializationFailure <- function(a, b, ...) {
  .init_failure_solution(a, if (missing(b)) NULL else b)
}

#' @export
print.InitializationFailure <- function(x, ...) {
  cat("Initialization of this problem failed (retcode InitialFailure).\n")
  class(x) <- setdiff(class(x), "InitializationFailure")
  print(x, ...)
}

.init_failure_solution <- function(prob, alg) {
  class(prob) <- setdiff(class(prob), "InitializationFailure")
  if (is.null(prob$tspan))
    return(.nl_new_solution("NonlinearSolution", prob, alg, prob$u0, "InitialFailure"))
  U <- matrix(as.numeric(prob$u0), nrow = 1L, dimnames = list(NULL, names(prob$u0)))
  .new_ode_solution(t = prob$tspan[1L], u = U, prob = prob, alg = alg, retcode = "InitialFailure")
}

# ---------------------------------------------------------------------------
# Diagnostics

#' Rank analysis of the initialization Jacobian
#'
#' `analyze_initialization_jacobian(prob)` diagnoses rank deficiency of the
#' initialization problem of `prob` from the singular value decomposition of
#' its residual Jacobian at the initial guess. It reports the unknowns that
#' span the null space (underdetermined directions, which the guesses decide)
#' and the equations that are locally redundant. The Jacobian is formed by
#' central finite differences.
#'
#' @param prob A problem built from a system (such as an `ODEProblem`), or
#'   an initialization problem from [InitializationProblem()].
#' @param rtol,atol A singular value `s` counts as zero when
#'   `s <= max(atol, rtol * smax)`.
#' @param threshold Only unknowns and equations whose participation exceeds
#'   this value are reported.
#' @param verbose Print the report.
#' @return A list with `jacobian`, `singular_values`, `rank`, `nullity`,
#'   `redundancy`, `underdetermined_unknowns` and `redundant_equations`; the
#'   last two are data frames with a `weight` in `[0, 1]`, the diagonal of
#'   the projector onto the corresponding null space.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); y <- Variable("y", t)
#' sys <- mtkcompile(System(list(D(x) ~ -x + y, 0 ~ x^2 + y^2 - 1), t, name = "dae"))
#' iprob <- InitializationProblem(sys, 0, list(x ~ 0.6), guesses = list(y ~ 1))
#' res <- analyze_initialization_jacobian(iprob, verbose = FALSE)
#' res$rank; res$nullity
#' @export
analyze_initialization_jacobian <- function(prob, rtol = 1e-8, atol = 0, threshold = 1e-3, verbose = TRUE) {
  empty <- list(jacobian = NULL, singular_values = numeric(), rank = 0L, nullity = 0L, redundancy = 0L,
                underdetermined_unknowns = data.frame(unknown = character(), weight = numeric()),
                redundant_equations = data.frame(equation = character(), weight = numeric()))
  iprob <- if (inherits(prob, "AbstractNonlinearProblem") || isTRUE(prob$is_initializeprob)) prob
           else prob$initialization_data$initializeprob
  if (is.null(iprob) && inherits(prob$sys, "System") && !is.null(.f(prob$sys, "iv"))) {
    kw <- prob$kwargs %||% list()
    iprob <- suppressWarnings(InitializationProblem(prob$sys, prob$tspan[1L], prob$op %||% list(),
                                                    guesses = kw$guesses %||% list(),
                                                    initialization_eqs = kw$initialization_eqs %||% list(),
                                                    jac = FALSE))
  }
  if (is.null(iprob)) {
    if (verbose) message("No initialization problem to analyze: the system is fully determined by its initial conditions.")
    return(empty)
  }
  u0 <- as.numeric(iprob$u0)
  if (!length(u0)) {
    if (verbose) message("The initialization problem has no unknowns to solve for.")
    return(empty)
  }
  F <- .nl_problem_fns(iprob)$F
  J <- NULL
  for (i in seq_along(u0)) {
    h <- .Machine$double.eps^(1 / 3) * max(1, abs(u0[i]))
    up <- u0; um <- u0
    up[i] <- u0[i] + h; um[i] <- u0[i] - h
    col <- (F(up) - F(um)) / (2 * h)
    if (is.null(J)) J <- matrix(0, length(col), length(u0))
    J[, i] <- col
  }
  nr <- nrow(J); nc <- ncol(J)
  if (!nr) J <- matrix(0, 0L, nc)
  s <- if (nr) svd(J, nu = nr, nv = nc) else list(d = numeric(), u = matrix(0, 0L, 0L), v = diag(nc))
  S <- s$d
  smax <- if (length(S)) max(S) else 0
  tol <- max(atol, rtol * smax)
  rank <- sum(S > tol)
  sv <- function(i) if (i <= length(S)) S[i] else 0
  col_null <- vapply(seq_len(nc), function(j) sv(j) <= tol, TRUE)
  row_null <- vapply(seq_len(nr), function(i) sv(i) <= tol, TRUE)
  nullity <- sum(col_null); redundancy <- sum(row_null)
  col_w <- if (nullity) rowSums(s$v[, col_null, drop = FALSE]^2) else numeric(nc)
  row_w <- if (redundancy) rowSums(s$u[, row_null, drop = FALSE]^2) else numeric(nr)
  syms <- names(iprob$u0) %||% paste0("u", seq_len(nc))
  eqs <- if (inherits(iprob$sys, "System")) vapply(equations(iprob$sys), .fmt, "") else as.character(seq_len(nr))
  pick <- function(labels, w, what) {
    keep <- which(w > threshold)
    keep <- keep[order(-w[keep])]
    out <- data.frame(labels[keep], w[keep], stringsAsFactors = FALSE)
    names(out) <- c(what, "weight")
    out
  }
  und <- pick(syms, col_w, "unknown")
  red <- pick(eqs, row_w, "equation")
  if (verbose) {
    lines <- c("Initialization Jacobian rank analysis",
               sprintf("  residual Jacobian: %d x %d, rank = %d, nullity = %d, redundancy = %d", nr, nc, rank,
                       nullity, redundancy))
    if (length(S)) {
      k <- min(5L, length(S))
      lines <- c(lines, paste0("  smallest ", k, " singular value(s): ",
                               paste(signif(utils::tail(S, k), 3), collapse = ", ")))
    }
    if (!nullity) lines <- c(lines, "  Full column rank at u0: no underdetermined unknowns.")
    else lines <- c(lines, "  Underdetermined unknowns (participation in [0, 1]):",
                    sprintf("    %8s  %s", format(round(und$weight, 4)), und$unknown))
    if (!redundancy) lines <- c(lines, "  Full row rank at u0: no redundant equations.")
    else lines <- c(lines, "  Redundant equations (participation in [0, 1]):",
                    sprintf("    %8s  %s", format(round(red$weight, 4)), red$equation))
    message(paste(lines, collapse = "\n"))
  }
  list(jacobian = J, singular_values = S, rank = rank, nullity = nullity, redundancy = redundancy,
       underdetermined_unknowns = und, redundant_equations = red)
}

# ---------------------------------------------------------------------------
# prob$ps[Initial(x)]

# The variable whose initial value Initial(x) is, with derivatives lowered.
.init_ps_source <- function(k) {
  src <- .init_source(k)
  if (.is_derivative(src)) default_toterm(src) else src
}

# The constant of the condition x ~ Initial(x), NULL when x has none.
.init_ps_constant <- function(prob, src) {
  k <- .key(src)
  # the constant the last initialization used, which follows numeric remakes
  data <- prob$initialization_data
  j <- match(k, data$initial_key %||% character())
  if (!is.na(j) && !is.null(data$ip_last)) return(unname(data$ip_last[data$initial_p[j]]))
  el <- prob$op[[k]] %||% .f(prob$sys, "initial_conditions")[[k]]
  v <- el$value
  if (is_sym(v) && .is_num(v)) v <- as.numeric(v)
  if (is.numeric(v) && length(v) == 1L && !is.na(v)) v else NULL
}

# prob$ps[k]: Initial(x) that is not a parameter of the problem reads the
# constant of x's initial condition, or the initial value of x.
.init_get_parameter <- function(prob, k) {
  if (!isinitial(k) || is.null(.f(prob$sys, "iv"))) return(NULL)
  src <- .init_ps_source(k)
  data <- prob$initialization_data
  if (!is.null(data$ip_last) && .key(src) %in% (data$initial_key %||% character()))
    return(.init_ps_constant(prob, src))
  # without an initialization problem the values of the unknowns are the
  # initial values
  ui <- match(.key(src), .keys(unknowns(prob$sys)))
  if (!is.na(ui) && is.numeric(prob$u0)) return(unname(prob$u0[ui]))
  .init_ps_constant(prob, src) %||% tryCatch(prob[src], error = function(e) NULL)
}

# prob$ps[k] <- v: a parameter takes the new value and the problem is
# initialized again (in ModelingToolkit this happens when it is solved).
# Initial(x) sets the constant of the condition x ~ Initial(x) when x has a
# constant initial condition and has no effect otherwise, as described in
# the initialization tutorial of ModelingToolkit.
.init_set_parameter <- function(prob, i, value) {
  sys <- prob$sys
  if (is.null(sys)) stop("setting parameters by name needs a problem built from a system", call. = FALSE)
  ks <- if (is.list(i) && !is_sym(i)) i else list(i)
  value <- as.numeric(unlist(value))
  if (length(value) == 1L) value <- rep(value, length(ks))
  if (length(value) != length(ks)) stop("one value per parameter is needed", call. = FALSE)
  pk <- .keys(.problem_params(sys))
  for (j in seq_along(ks)) {
    k <- as_sym(ks[[j]])
    if (.key(k) %in% pk) {
      prob <- remake(prob, p = smap_set(symmap(), k, value[j]))
    } else if (isinitial(k) && !is.null(.f(sys, "iv"))) {
      src <- .init_ps_source(k)
      if (!is.null(.init_ps_constant(prob, src))) prob <- remake(prob, u0 = smap_set(symmap(), src, value[j]))
    } else {
      stop(format(k), " is not a parameter of the problem", call. = FALSE)
    }
  }
  prob
}

# ---------------------------------------------------------------------------
# Registration (see "Extension hooks" in tools/ARCHITECTURE.md)

.mtk$hooks$initialize_problem <- .init_hook
.mtk$hooks$initialize_parameters <- .init_parameters_hook
.mtk$hooks$get_parameter <- .init_get_parameter
.mtk$hooks$set_parameter <- .init_set_parameter
