# The initialization system, after ModelingToolkitBase/src/systems/
# nonlinear/initializesystem.jl (generate_initializesystem,
# generate_initializesystem_timevarying, generate_initializesystem_timeindependent,
# register_dd_var!, initsys_sort_system_bindings!,
# timevaring_initsys_process_op!) and systems/abstractsystem.jl
# (Initial, add_initialization_parameters, isinitial).
#
# The initialization system of a time dependent system is a system of
# nonlinear equations for the initial values of its unknowns, of the
# derivatives of its differential variables (the dummy derivative variables
# xˍt) and of its solvable parameters. A constant initial condition x => 1
# becomes the equation x ~ Initial(x) with the parameter Initial(x) = 1, so
# the constant can change without regenerating the equations.

# Initial(x) in the canonical form of ModelingToolkit: derivatives are
# replaced by their lowered variable (Initial(D(x)) is Initial(xˍt)) and an
# Initial is never wrapped twice.
.init_initial <- function(x) {
  x <- as_sym(x)
  if (isinitial(x)) return(x)
  if (.is_derivative(x)) x <- default_toterm(x)
  Initial(x)
}

# The variable an Initial parameter refers to.
.init_source <- function(ix) .args(ix)[[1L]]

.init_canon <- function(e) {
  .nl_rewrite_call(e, "Initial", function(x) .init_initial(x))
}

# Initial(...) applications referenced by the equations, observed equations,
# initialization equations, bindings and initial conditions of a system.
.init_referenced_initials <- function(sys) {
  if (!exists("Initial", envir = .mtk$vars, inherits = FALSE)) return(list())
  sides <- function(eqs) {
    eqs <- Filter(function(e) inherits(e, "Equation"), eqs)
    c(lapply(eqs, `[[`, "lhs"), lapply(eqs, `[[`, "rhs"))
  }
  vals <- function(m) Filter(is_sym, lapply(unclass(m), `[[`, "value"))
  exprs <- c(sides(.f(sys, "eqs")), sides(.f(sys, "observed")), sides(.f(sys, "initialization_eqs")),
             vals(.f(sys, "bindings")), vals(.f(sys, "initial_conditions")))
  .nl_collect_calls(exprs, "Initial")
}

.init_add_initials <- function(sys, nodes) {
  ps <- .f(sys, "ps")
  have <- .keys(ps)
  ics <- .f(sys, "initial_conditions") %||% symmap()
  binds <- .f(sys, "bindings") %||% symmap()
  for (ix in nodes) {
    k <- .key(ix)
    if (!(k %in% have)) {
      ps[[length(ps) + 1L]] <- ix
      have <- c(have, k)
    }
    # the value of Initial(x) is the initial value of x; it is filled in by
    # initialization when x itself is solved for, hence the missing binding
    if (is.null(ics[[k]])) ics <- smap_set(ics, ix, .init_source(ix))
    if (is.null(binds[[k]])) binds <- smap_set(binds, ix, NA)
  }
  .update_system(sys, ps = ps, initial_conditions = ics, bindings = binds)
}

#' Initial value parameters
#'
#' `add_initialization_parameters(sys)` adds a parameter `Initial(x)` for
#' every unknown, every observed variable, their first derivatives and every
#' solvable parameter of a flat system, as ModelingToolkit does when a system
#' is completed. The value of `Initial(x)` is the initial value of `x`.
#'
#' In this package `complete()` and `mtkcompile()` add only the `Initial`
#' parameters that the model refers to (in equations, bindings such as
#' `q ~ Initial(y)`, or initialization equations), so the parameter vector of
#' a problem stays free of unused entries; call this function to add them
#' all. The initialization system built by [generate_initializesystem()]
#' always has the `Initial` parameters of the constrained variables.
#'
#' @param sys A flat `System` (no subsystems).
#' @param split Unused; kept for compatibility.
#' @return The system with the added parameters.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); k <- Parameter("k")
#' sys <- System(list(D(x) ~ -k * x), t, name = "decay")
#' sys2 <- add_initialization_parameters(sys)
#' parameters(sys2, initial_parameters = TRUE)
#' @export
add_initialization_parameters <- function(sys, split = TRUE) {
  if (length(.f(sys, "systems"))) stop("add_initialization_parameters() needs a flat system", call. = FALSE)
  if (isTRUE(.f(sys, "is_initializesystem"))) return(sys)
  vars <- c(.f(sys, "unknowns"), lapply(.f(sys, "observed"), `[[`, "lhs"))
  iv <- .f(sys, "iv")
  if (!is.null(iv) && !isTRUE(.f(sys, "is_discrete")))
    vars <- c(vars, lapply(Filter(function(v) .type(v) == "Function", vars), .make_derivative, iv = iv, order = 1L))
  binds <- .f(sys, "bindings")
  for (el in unclass(binds))
    if (!is.null(el$var) && !isinitial(el$var) && .init_is_missing(el$value)) vars <- c(vars, list(el$var))
  .init_add_initials(sys, .unique_syms(lapply(vars, .init_initial)))
}

.init_is_missing <- function(v) !is.null(v) && !is_sym(v) && length(v) == 1L && is.na(v)

# The complete hook: make the Initial parameters a model refers to
# parameters of the system, in canonical form.
.init_complete_hook <- function(sys) {
  if (isTRUE(.f(sys, "is_initializesystem"))) return(sys)
  nodes <- .init_referenced_initials(sys)
  if (!length(nodes)) return(sys)
  canon <- lapply(nodes, function(ix) .init_initial(.init_source(ix)))
  if (!identical(.keys(canon), .keys(nodes))) {
    rew <- function(eqs) lapply(eqs, function(e) if (inherits(e, "Equation")) .init_canon(e) else e)
    remap <- function(m) {
      for (k in names(m)) if (is_sym(m[[k]]$value)) m[[k]]$value <- .init_canon(m[[k]]$value)
      m
    }
    sys <- .update_system(sys, eqs = rew(.f(sys, "eqs")), observed = rew(.f(sys, "observed")),
                          initialization_eqs = rew(.f(sys, "initialization_eqs")),
                          bindings = remap(.f(sys, "bindings")),
                          initial_conditions = remap(.f(sys, "initial_conditions")))
  }
  .init_add_initials(sys, .unique_syms(canon))
}

# ---------------------------------------------------------------------------
# generate_initializesystem

# Substitution that leaves operator applications such as Initial(x) intact.
.init_subber <- function(from, to, fixpoint = TRUE) {
  if (!length(from)) return(function(e) e)
  base <- .rule(from, to)
  to_ops <- .nl_op_nodes(to)
  function(e) {
    ops <- c(to_ops, .nl_op_nodes(list(e)))
    rule <- if (length(ops)) .rule_merge(base, .rule(ops, ops)) else base
    if (fixpoint) .subst_fixpoint(e, rule) else .subst(e, rule)
  }
}

#' Generate the initialization system
#'
#' `generate_initializesystem(sys, op = ...)` returns the system of
#' nonlinear equations that determines consistent initial values of a
#' time dependent system. Its unknowns are the unknowns of `sys`, the
#' derivatives of its differential variables (named `xˍt`) and its solvable
#' parameters (those bound to `NA`). Its equations are
#'
#' * the algebraic equations of `sys`, and `xˍt ~ f` for each differential
#'   equation `D(x) ~ f`, with observed variables substituted;
#' * the bindings of variables;
#' * one equation per initial condition in `op` or in the system: a constant
#'   `x ~ 1` becomes `x ~ Initial(x)` with the parameter `Initial(x)` set to
#'   1, an expression is used as it is; conditions on observed variables and
#'   on derivatives (`D(x) ~ 0`) are allowed;
#' * the initialization equations of `sys` and those in `initialization_eqs`.
#'
#' Guesses become the `guesses` of the returned system, and the parameter
#' values and `Initial` values its initial conditions, so it can be passed to
#' [mtkcompile()] with `fully_determined = FALSE` and then to
#' [NonlinearProblem()]. [InitializationProblem()] does all of this.
#'
#' For a time independent system (or `time_dependent_init = FALSE`) the
#' initialization system solves only for the solvable parameters; the
#' unknowns enter through `Initial(x)`.
#'
#' A value of `NULL` in `op` (as the pair `list(x, NULL)`) removes the
#' initial condition of `x`, and `NA` makes a parameter solvable.
#'
#' @param sys A complete `System`.
#' @param op Operating point: initial conditions and parameter values.
#' @param initialization_eqs Further initialization equations.
#' @param guesses Further guesses, as a variable map.
#' @param default_dd_guess Guess for the derivative variables.
#' @param algebraic_only Ignore the initialization equations of the system
#'   (those in `initialization_eqs` are still used).
#' @param time_dependent_init Build the time dependent initialization.
#' @param name Name of the initialization system.
#' @param ... Unused.
#' @return A time independent `System` with `is_initializesystem` set.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); y <- Variable("y", t)
#' sys <- mtkcompile(System(list(D(x) ~ -x + y, 0 ~ x^2 + y^2 - 1), t, name = "dae"))
#' isys <- generate_initializesystem(sys, op = list(x ~ 0.6), guesses = list(y ~ 1))
#' equations(isys)
#' unknowns(isys)
#' @export
generate_initializesystem <- function(sys, op = list(), initialization_eqs = list(), guesses = list(),
                                      default_dd_guess = 0, algebraic_only = FALSE,
                                      time_dependent_init = !is.null(get_iv(sys)), name = nameof(sys), ...) {
  .need_complete(sys, "generate_initializesystem")
  if (isTRUE(time_dependent_init)) {
    if (is.null(.f(sys, "iv")))
      stop("a time dependent initialization requires a time dependent system", call. = FALSE)
    .init_generate_timevarying(sys, op, initialization_eqs, guesses, default_dd_guess, algebraic_only, name)
  } else {
    .init_generate_timeindependent(sys, op, initialization_eqs, guesses, algebraic_only, name)
  }
}

# The operating point of initialization: the initial conditions of the
# system overridden by the user's map, with names resolved and derivative
# keys lowered. Entries whose value is NULL remove a condition.
.init_operating_point <- function(sys, op) {
  allv <- c(unknowns(sys), .parameters_of(sys, TRUE), observables(sys))
  smap_merge(initial_conditions(sys), .resolve_named(as_symmap(op), allv))
}

.init_generate_timevarying <- function(sys, op, initialization_eqs, guesses_in, default_dd_guess,
                                       algebraic_only, name) {
  iv <- .f(sys, "iv")
  us <- unknowns(sys)
  obs <- observed(sys)
  eqs <- equations(sys)
  obs_lhs <- lapply(obs, `[[`, "lhs")
  obs_keys <- .keys(obs_lhs)
  subber <- .init_subber(obs_lhs, lapply(obs, `[[`, "rhs"))
  ps <- .parameters_of(sys, TRUE)
  binds <- bindings(sys)
  allv <- c(us, ps, obs_lhs)
  gs <- smap_merge(guesses(sys), .resolve_named(as_symmap(guesses_in), allv))
  user_keys <- names(.resolve_named(as_symmap(op), allv))
  op <- .init_operating_point(sys, op)

  vars <- us
  var_keys <- .keys(us)
  add_var <- function(v) {
    if (!(.key(v) %in% var_keys)) {
      vars[[length(vars) + 1L]] <<- v
      var_keys <<- c(var_keys, .key(v))
    }
  }
  eqs_ics <- list()
  push <- function(lhs, rhs) {
    if (!identical(.key(lhs), .key(rhs))) eqs_ics[[length(eqs_ics) + 1L]] <<- Equation(lhs, rhs)
  }
  pvals <- symmap()      # values of the parameters of the initialization system
  newbinds <- symmap()
  solvable <- character()

  # Bindings: NA makes a parameter solvable; bindings of variables are
  # equations; bindings of parameters stay bindings.
  for (el in unclass(binds)) {
    k <- el$var; v <- el$value
    if (is.null(k) || isinitial(k)) next
    kk <- .key(k)
    if (.init_is_missing(v)) {
      solvable <- c(solvable, kk)
      next
    }
    if (kk %in% var_keys || kk %in% obs_keys) push(subber(k), subber(as_sym(v)))
    else if (is_sym(v)) newbinds[[kk]] <- el
    else pvals[[kk]] <- el
  }
  # A late binding of NA in the operating point also makes a parameter solvable.
  for (el in unclass(op)) {
    if (!is.null(el$var) && .init_is_missing(el$value) && isparameter(el$var) && !isinitial(el$var))
      solvable <- union(solvable, .key(el$var))
  }
  for (p in ps) if (.key(p) %in% solvable) add_var(p)
  par_keys <- setdiff(.keys(ps), solvable)

  # Dummy derivative variables: xˍt ~ f for each D(x) ~ f. The differential
  # equations are added as they are rather than substituted, so singular
  # systems keep all their unknowns.
  dfrom <- list(); dto <- list()
  dd_keys <- character()
  register_dd <- function(d) {
    ttk <- default_toterm(d)
    dfrom[[length(dfrom) + 1L]] <<- d
    dto[[length(dto) + 1L]] <<- ttk
    # register_dd_var!: the guess is written unless d or ttk has one, also when
    # ttk is an unknown of the system (a state from order lowering)
    if (!smap_has(gs, d) && !smap_has(gs, ttk)) gs <<- smap_set(gs, ttk, default_dd_guess)
    if (!(.key(ttk) %in% c(var_keys, obs_keys))) {
      add_var(ttk)
      dd_keys <<- c(dd_keys, .key(ttk))
    }
    ttk
  }
  residual <- list()
  for (eq in eqs) {
    if (.is_derivative(eq$lhs)) {
      ttk <- register_dd(eq$lhs)
      push(subber(ttk), subber(eq$rhs))
    } else residual[[length(residual) + 1L]] <- eq
  }
  der_subber <- function(e) {
    ds <- Filter(.is_derivative, .collect_vars(list(e)))
    for (d in ds) if (!(.key(d) %in% .keys(dfrom))) register_dd(d)
    if (!length(ds)) return(e)
    .subst(e, .rule(dfrom, dto))
  }
  for (eq in residual) {
    eqs_ics[[length(eqs_ics) + 1L]] <- Equation(subber(der_subber(eq$lhs)), subber(der_subber(eq$rhs)))
  }

  # The operating point: constants become x ~ Initial(x), expressions are
  # equations as they are.
  initials <- symmap()   # Initial(k) => value
  obs_dd <- NULL
  for (el in unclass(op)) {
    k <- el$var; v <- el$value
    if (is.null(k) || is.null(v)) next
    kk <- .key(k)
    if (isinitial(k)) {
      # Initial(x) => x is the default: the value follows the condition on x
      if (is.numeric(v) || .is_num(v)) initials <- smap_set(initials, .init_initial(.init_source(k)), as.numeric(v))
      next
    }
    if (.init_is_missing(v) || identical(kk, .key(iv))) next
    if (kk %in% par_keys) {
      # the binding of a bound parameter wins; it cannot be given a value
      if (is.null(newbinds[[kk]])) pvals[[kk]] <- el
      next
    }
    if (kk %in% var_keys) subk <- k
    else if (kk %in% obs_keys) subk <- subber(k)
    else {
      # the derivative of an observed variable, or of a variable that has
      # no differential equation
      if (is.null(obs_dd)) obs_dd <- .init_derivative_terms(us, obs, iv)
      hit <- obs_dd[[kk]]
      # initial conditions kept by the system for variables it no longer
      # has are ignored; a key given by the user is an error
      if (is.null(hit) && !(kk %in% user_keys)) next
      if (is.null(hit))
        stop("Expected an `Initial` parameter to exist for variable `", format(k), "`, but did not find one. ",
             "This can be because:\n1. The problem constructor was passed an initial condition for `",
             format(k), "` but the variable does not exist in the system. Either ensure that `", format(k),
             "` is a variable/parameter of the system, or change the set of initial conditions provided ",
             "to the problem constructor.\n2. An initial condition was provided for a bound parameter. ",
             "Bound parameters cannot be given initial conditions. Their values are determined solely by ",
             "the binding relation.", call. = FALSE)
      # an initial condition for the derivative of an algebraic variable
      # would need index reduction; it is ignored, as in ModelingToolkit
      if (identical(hit, "ignore")) next
      subk <- subber(der_subber(hit))
    }
    if (is.numeric(v) || .is_num(v)) {
      ik <- .init_initial(k)
      if (!smap_has(initials, ik)) initials <- smap_set(initials, ik, as.numeric(v))
      if (kk %in% var_keys && !smap_has(gs, k)) gs <- smap_set(gs, k, as.numeric(v))
      push(subk, ik)
    } else {
      push(subk, subber(der_subber(as_sym(v))))
    }
  }

  ieqs <- .as_equations(initialization_eqs)
  if (!algebraic_only) ieqs <- c(initialization_equations(sys), ieqs)
  for (eq in ieqs) {
    eq <- .init_canon(eq)
    eqs_ics[[length(eqs_ics) + 1L]] <- Equation(subber(der_subber(eq$lhs)), subber(der_subber(eq$rhs)))
  }

  # Initial parameters: those of the constants above and any referenced by
  # the equations or bindings. The Initial(x) of a variable without a
  # constant condition is the value initialization finds for x, so it is x
  # itself inside this system.
  bind_values <- function() Filter(is_sym, lapply(unclass(newbinds), `[[`, "value"))
  used <- .nl_collect_calls(c(lapply(eqs_ics, `[[`, "lhs"), lapply(eqs_ics, `[[`, "rhs"), bind_values()), "Initial")
  cfrom <- list(); cto <- list()
  for (ix in used) {
    if (smap_has(initials, ix)) next
    src <- .init_source(ix)
    sk <- .key(src)
    v <- op[[sk]]$value
    if (is_sym(v) && .is_num(v)) v <- as.numeric(v)
    if (is.numeric(v) && length(v) == 1L && !is.na(v)) initials <- smap_set(initials, ix, v)
    else if (sk %in% var_keys || sk %in% obs_keys) {
      cfrom[[length(cfrom) + 1L]] <- ix
      cto[[length(cto) + 1L]] <- if (sk %in% var_keys) src else subber(src)
    } else initials <- smap_set(initials, ix, 0)
  }
  if (length(cfrom)) {
    crule <- .rule(cfrom, cto)
    eqs_ics <- lapply(eqs_ics, .subst, rule = crule)
    for (k in names(newbinds)) if (is_sym(newbinds[[k]]$value)) newbinds[[k]]$value <- .subst(newbinds[[k]]$value, crule)
    # a parameter bound to an expression of the unknowns is replaced by it,
    # so the equations show their dependence on those unknowns
    dep <- Filter(function(el) is_sym(el$value) && any(.nl_atomic_keys(el$value) %in% var_keys), unclass(newbinds))
    if (length(dep)) {
      psub <- .init_subber(unname(lapply(dep, `[[`, "var")), unname(lapply(dep, `[[`, "value")))
      eqs_ics <- lapply(eqs_ics, function(e) Equation(psub(e$lhs), psub(e$rhs)))
    }
  }
  init_nodes <- smap_vars(initials)
  pars <- .unique_syms(c(Filter(function(p) .key(p) %in% par_keys && !isinitial(p), ps), list(iv),
                         unname(init_nodes)))
  ics <- smap_merge(pvals, initials)
  # guesses of observed variables cannot be used (they are substituted)
  for (k in names(gs)) if (!(k %in% var_keys)) gs[[k]] <- NULL
  isys <- System(eqs_ics, NULL, unknowns = vars, ps = pars, bindings = newbinds,
                 initial_conditions = ics, guesses = gs, name = name, checks = FALSE,
                 discover_from_metadata = FALSE)
  .update_system(isys, is_initializesystem = TRUE,
                 metadata = list(initialization = list(time_dependent = TRUE, solvable = solvable,
                                                       dd = dd_keys, system_unknowns = .keys(us))))
}

# Keys of the lowered derivatives that may appear in an operating point
# without being unknowns: the derivative of an observed variable maps to the
# time derivative of its expression, and the derivative of an unknown that
# has no differential equation is marked "ignore".
.init_derivative_terms <- function(us, obs, iv) {
  out <- list()
  for (o in obs) {
    if (.type(o$lhs) != "Function") next
    out[[.key(default_toterm(.make_derivative(o$lhs, iv, 1L)))]] <- expand_derivatives(.make_derivative(o$rhs, iv, 1L))
  }
  for (u in us) {
    k <- .key(default_toterm(.make_derivative(u, iv, 1L)))
    if (is.null(out[[k]])) out[[k]] <- "ignore"
  }
  out
}

.init_generate_timeindependent <- function(sys, op, initialization_eqs, guesses_in, algebraic_only, name) {
  us <- unknowns(sys)
  obs_lhs <- observables(sys)
  og_keys <- .keys(c(us, obs_lhs))
  ps <- .parameters_of(sys, TRUE)
  binds <- bindings(sys)
  iv <- .f(sys, "iv")
  allv <- c(us, ps, obs_lhs)
  gs <- smap_merge(guesses(sys), .resolve_named(as_symmap(guesses_in), allv))
  op <- .init_operating_point(sys, op)
  to_initial <- .init_subber(c(us, obs_lhs), lapply(c(us, obs_lhs), .init_initial), fixpoint = FALSE)

  vars <- list(); eqs_ics <- list()
  pvals <- symmap(); newbinds <- symmap(); initials <- symmap()
  solvable <- character()
  for (el in unclass(binds)) {
    k <- el$var; v <- el$value
    if (is.null(k) || isinitial(k)) next
    kk <- .key(k)
    if (.init_is_missing(v)) solvable <- c(solvable, kk)
    else if (kk %in% og_keys) next
    else if (is_sym(v)) newbinds[[kk]] <- el
    else pvals[[kk]] <- el
  }
  for (el in unclass(op))
    if (!is.null(el$var) && .init_is_missing(el$value) && isparameter(el$var) && !isinitial(el$var))
      solvable <- union(solvable, .key(el$var))
  for (p in ps) if (.key(p) %in% solvable) vars[[length(vars) + 1L]] <- p
  par_keys <- setdiff(.keys(ps), solvable)
  valid <- character()
  for (el in unclass(op)) {
    k <- el$var; v <- el$value
    if (is.null(k) || is.null(v) || .init_is_missing(v)) next
    kk <- .key(k)
    if (isinitial(k)) {
      if (is.numeric(v) || .is_num(v)) initials <- smap_set(initials, .init_initial(.init_source(k)), as.numeric(v))
      next
    }
    if (kk %in% og_keys) {
      # unknowns are not solved for here; their values are Initial parameters
      if (is.numeric(v) || .is_num(v)) {
        initials <- smap_set(initials, .init_initial(k), as.numeric(v))
        valid <- c(valid, .key(.init_initial(k)))
      }
      next
    }
    if (kk %in% par_keys) {
      if (is.null(newbinds[[kk]])) pvals[[kk]] <- el
      next
    }
    if (!(kk %in% solvable)) next
    if (is.numeric(v) || .is_num(v)) {
      ik <- .init_initial(k)
      initials <- smap_set(initials, ik, as.numeric(v))
      if (!smap_has(gs, k)) gs <- smap_set(gs, k, as.numeric(v))
      eqs_ics[[length(eqs_ics) + 1L]] <- Equation(k, ik)
    } else eqs_ics[[length(eqs_ics) + 1L]] <- Equation(k, to_initial(as_sym(v)))
  }
  ieqs <- .as_equations(initialization_eqs)
  if (!algebraic_only) ieqs <- c(initialization_equations(sys), ieqs)
  allowed <- c(.keys(ps), if (!is.null(iv)) .key(iv))
  for (eq in ieqs) {
    eq <- .init_canon(eq)
    eq <- Equation(to_initial(eq$lhs), to_initial(eq$rhs))
    used <- .nl_atomic_keys(list(eq$lhs, eq$rhs))
    inits <- .nl_collect_calls(list(eq$lhs, eq$rhs), "Initial")
    bad <- setdiff(used, c(allowed, .keys(inits)))
    if (length(bad))
      stop("Initialization equations of a time independent system may only involve parameters; found ",
           "unknowns in ", .fmt(eq), call. = FALSE)
    # only equations whose Initial parameters have values are usable
    if (all(.keys(inits) %in% names(initials))) eqs_ics[[length(eqs_ics) + 1L]] <- eq
  }
  used <- .nl_collect_calls(c(lapply(eqs_ics, `[[`, "lhs"), lapply(eqs_ics, `[[`, "rhs")), "Initial")
  pars <- .unique_syms(c(Filter(function(p) .key(p) %in% par_keys && !isinitial(p), ps), used))
  ics <- smap_merge(pvals, initials)
  for (k in names(gs)) if (!(k %in% solvable)) gs[[k]] <- NULL
  isys <- System(eqs_ics, NULL, unknowns = vars, ps = pars, bindings = newbinds,
                 initial_conditions = ics, guesses = gs, name = name, checks = FALSE,
                 discover_from_metadata = FALSE)
  .update_system(isys, is_initializesystem = TRUE,
                 metadata = list(initialization = list(time_dependent = FALSE, solvable = solvable,
                                                       dd = character(), system_unknowns = .keys(us))))
}

# ---------------------------------------------------------------------------
# Registration (see "Extension hooks" in tools/ARCHITECTURE.md)

.mtk$hooks$complete <- .init_complete_hook
