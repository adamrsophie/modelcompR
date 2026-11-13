# The System type, after ModelingToolkitBase/src/systems/system.jl and
# abstractsystem.jl. A System is an S3 list; internal code reads its fields
# with .subset2() (or .f()) because `$` is overloaded for namespaced variable
# access, as getproperty is in Julia.

.SYSTEM_FIELDS <- c(
  "tag", "eqs", "noise_eqs", "jumps", "constraints", "costs", "consolidate",
  "unknowns", "ps", "brownians", "poissonians", "iv", "observed", "name",
  "description", "bindings", "initial_conditions", "guesses", "systems",
  "initialization_eqs", "continuous_events", "discrete_events", "connector_type",
  "assertions", "metadata", "is_dde", "tstops", "tspan", "inputs", "outputs",
  "tearing_state", "namespacing", "complete", "index_cache", "parent",
  "initializesystem", "is_initializesystem", "is_discrete", "state_priorities",
  "irreducibles", "maybe_zeros", "isscheduled", "schedule",
  "analytically_integrated", "scopes", "ignored_connections", "preface")

.f <- function(sys, field) .subset2(sys, field)

.new_system <- function(fields) {
  stopifnot(all(names(fields) %in% .SYSTEM_FIELDS))
  sys <- stats::setNames(vector("list", length(.SYSTEM_FIELDS)), .SYSTEM_FIELDS)
  for (n in names(fields)) sys[n] <- list(fields[[n]])
  class(sys) <- "System"
  sys
}

# Copy a system with some fields replaced.
.update_system <- function(sys, ...) {
  upd <- list(...)
  cls <- class(sys)
  sys <- unclass(sys)
  for (n in names(upd)) {
    if (!(n %in% .SYSTEM_FIELDS)) stop("unknown System field ", n, call. = FALSE)
    sys[n] <- list(upd[[n]])
  }
  class(sys) <- cls
  sys
}

#' Build a system of equations
#'
#' `System()` is the single system type of ModelingToolkit 11. A system with an
#' independent variable `iv` describes differential (or differential algebraic,
#' stochastic, jump, delay or discrete) equations; without one it describes a
#' nonlinear system, or an optimization problem when `costs` are given.
#'
#' Unknowns and parameters are inferred from the equations when they are not
#' given. Defaults declared on the variables become initial conditions
#' (numbers) or bindings (expressions and `NA`), and guesses are collected the
#' same way.
#'
#' @param eqs Equations: `Equation` objects, formulas such as `D(x) ~ -x`,
#'   connections, or a list of them.
#' @param iv The independent variable, or `NULL`.
#' @param unknowns,ps Unknowns and parameters, or `NULL` to infer them.
#' @param brownians Brownian variables appearing in `eqs`.
#' @param name The system name, a string. Required.
#' @param systems A list of subsystems.
#' @param observed Observed equations `y ~ expression`.
#' @param bindings,initial_conditions,guesses Variable maps: lists of
#'   formulas `x ~ 1`, of pairs `list(x, 1)`, or named lists.
#' @param defaults Legacy name for `initial_conditions`.
#' @param initialization_eqs Equations that hold only at the initial time.
#' @param continuous_events,discrete_events Events, see
#'   [SymbolicContinuousCallback()] and [SymbolicDiscreteCallback()].
#' @param constraints,costs,consolidate Constraints, costs and the function
#'   combining costs, for optimization and optimal control.
#' @param noise_eqs Noise coefficients of a stochastic system.
#' @param jumps Jumps of a jump system.
#' @param poissonians Poissonian variables.
#' @param connector_type Connector type, set by [connector()].
#' @param assertions Named list mapping conditions to error messages.
#' @param metadata A list of metadata.
#' @param description A string.
#' @param tspan Default time span `c(t0, t1)`.
#' @param tstops Extra integrator stop times.
#' @param inputs,outputs Input and output variables.
#' @param state_priorities,irreducibles State selection hints.
#' @param is_discrete Whether the system is a discrete time system.
#' @param discover_from_metadata Collect defaults and guesses from metadata.
#' @param checks Run consistency checks.
#' @param frame,fixed_point For the multibody `System()` of symcompR (SymPy's
#'   physics.mechanics): a call with no arguments, with a reference frame as
#'   its first argument or with `frame` or `fixed_point` goes to
#'   `symcompR::System()`, so attaching both packages is safe.
#' @return A `System`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)")
#' parameters("tau")
#' sys <- System(list(D(x) ~ (1 - x) / tau), t, name = "sys")
#' sys
#' @export
System <- function(eqs = list(), iv = NULL, unknowns = NULL, ps = NULL, brownians = NULL,
                   name = NULL, systems = list(), observed = list(), bindings = NULL,
                   initial_conditions = NULL, guesses = NULL, defaults = NULL,
                   initialization_eqs = list(), continuous_events = list(),
                   discrete_events = list(), constraints = list(), costs = list(),
                   consolidate = NULL, noise_eqs = NULL, jumps = list(),
                   poissonians = NULL, connector_type = NULL, assertions = list(),
                   metadata = list(), description = "", tspan = NULL, tstops = list(),
                   inputs = list(), outputs = list(), state_priorities = NULL,
                   irreducibles = NULL, is_discrete = FALSE,
                   discover_from_metadata = TRUE, checks = TRUE, frame = NULL, fixed_point = NULL) {
  if ((nargs() == 0L || !missing(frame) || !missing(fixed_point) || inherits(eqs, "physics_obj")) &&
      .symcompr_has("System"))
    return(.symcompr_forward("System", sys.call(), parent.frame()))
  if (is.null(name)) stop("a system needs a name: System(..., name = \"sys\")", call. = FALSE)
  name <- as.character(name)
  if (!is.null(iv) && !is_sym(iv)) {
    if (is.list(iv) && length(iv) >= 1L) iv <- iv[[1L]] else iv <- as_sym(iv)
  }
  eqs <- .as_equations(eqs)
  observed <- .as_equations(observed)
  initialization_eqs <- .as_equations(initialization_eqs)
  if (inherits(systems, "System")) systems <- list(systems)
  systems <- unname(as.list(systems))
  costs <- lapply(if (is_sym(costs) || is.numeric(costs)) list(costs) else costs, as_sym)
  constraints <- lapply(constraints, function(c) if (inherits(c, "Inequality")) c else as_equation(c))
  if (!length(unique(vapply(systems, nameof, ""))) == length(systems))
    stop("subsystem names must be unique in system ", name, call. = FALSE)

  sym_eqs <- Filter(function(e) inherits(e, "Equation"), eqs)
  found <- .collect_vars(c(lapply(sym_eqs, `[[`, "lhs"), lapply(sym_eqs, `[[`, "rhs"),
                           lapply(observed, `[[`, "rhs"), costs,
                           lapply(constraints, `[[`, "lhs"), lapply(constraints, `[[`, "rhs")))
  # jumps and poissonian rates contribute variables (collect_vars! on jumps)
  jx <- .mtk$hooks$jump_exprs
  if (!is.null(jx)) found <- .unique_syms(c(found, .collect_vars(jx(jumps, found))))
  found <- .with_default_vars(found)
  inferred <- .classify_vars(found, iv)
  obs_lhs <- lapply(observed, `[[`, "lhs")
  if (is.null(unknowns)) {
    unknowns <- inferred$unknowns
    if (length(obs_lhs)) unknowns <- Filter(function(u) !.contains_key(obs_lhs, u), unknowns)
  } else {
    unknowns <- .flatten_vars(unknowns)
    if (!is.null(iv)) unknowns <- Filter(function(u) !.is_delay_or_evalat(u, iv), unknowns)
  }
  if (is.null(ps)) {
    ps <- inferred$ps
    # parameters referenced only by events are added by the event processing
  } else {
    ps <- .flatten_vars(ps)
  }
  if (is.null(brownians)) brownians <- inferred$brownians else brownians <- .flatten_vars(brownians)
  if (is.null(poissonians)) poissonians <- inferred$poissonians else poissonians <- .flatten_vars(poissonians)

  initial_conditions <- smap_merge(as_symmap(defaults), as_symmap(initial_conditions))
  bindings <- as_symmap(bindings)
  guesses <- as_symmap(guesses)
  inputs <- .flatten_vars(inputs)
  outputs <- .flatten_vars(outputs)
  sp <- as_symmap(state_priorities)
  irr <- .flatten_vars(irreducibles)
  if (discover_from_metadata) {
    md <- .discover_metadata(c(unknowns, ps, obs_lhs), initial_conditions, bindings, guesses)
    initial_conditions <- md$initial_conditions
    bindings <- md$bindings
    guesses <- md$guesses
    for (u in unknowns) {
      if (isinput(u) && !.contains_key(inputs, u)) inputs[[length(inputs) + 1L]] <- u
      if (isoutput(u) && !.contains_key(outputs, u)) outputs[[length(outputs) + 1L]] <- u
      pr <- .meta_get(u, "state_priority")
      if (!is.null(pr) && !smap_has(sp, u)) sp <- smap_set(sp, u, pr)
      if (isirreducible(u) && !.contains_key(irr, u)) irr[[length(irr) + 1L]] <- u
    }
  }
  for (s in systems) {
    for (v in .f(s, "inputs")) inputs[[length(inputs) + 1L]] <- renamespace(s, v)
    for (v in .f(s, "outputs")) outputs[[length(outputs) + 1L]] <- renamespace(s, v)
  }
  if (is.null(iv) && length(c(continuous_events, discrete_events)))
    stop("time independent systems cannot have events", call. = FALSE)

  sys <- .new_system(list(
    tag = .next_id(), eqs = eqs, noise_eqs = noise_eqs,
    jumps = if (inherits(jumps, "AbstractJump")) list(jumps) else as.list(jumps),
    constraints = constraints, costs = costs, consolidate = consolidate,
    unknowns = unknowns, ps = ps, brownians = brownians, poissonians = poissonians,
    iv = iv, observed = observed, name = name, description = description,
    bindings = bindings, initial_conditions = initial_conditions, guesses = guesses,
    systems = systems, initialization_eqs = initialization_eqs,
    continuous_events = list(), discrete_events = list(),
    connector_type = connector_type, assertions = as.list(assertions),
    metadata = as.list(metadata), is_dde = .check_dde(eqs, iv),
    tstops = as.list(tstops), tspan = tspan, inputs = .unique_syms(inputs),
    outputs = .unique_syms(outputs), namespacing = TRUE, complete = FALSE,
    state_priorities = sp, irreducibles = irr, maybe_zeros = list(),
    isscheduled = FALSE, is_discrete = is_discrete, is_initializesystem = FALSE,
    analytically_integrated = list(), scopes = .collect_scopes(c(sym_eqs, observed))))
  # Events are normalized by the events module when it is loaded; the plain
  # lists are kept otherwise.
  sys <- .set_events(sys, continuous_events, discrete_events)
  # checks: TRUE/FALSE or a bitmask of CheckComponents (2) and CheckUnits (4)
  if (.checks_on(checks, 2L) && !is.null(iv)) .check_system(sys)
  if (.checks_on(checks, 4L) && !is.null(.mtk$hooks$check_units)) .mtk$hooks$check_units(sys)
  sys
}

.checks_on <- function(checks, bit) {
  isTRUE(checks) || (is.numeric(checks) && (checks == 1 || bitwAnd(as.integer(checks), bit) > 0L))
}

.set_events <- function(sys, continuous_events, discrete_events) {
  ce <- .normalize_events(continuous_events, "continuous")
  de <- .normalize_events(discrete_events, "discrete")
  extra <- .event_parameters(c(ce, de), .f(sys, "unknowns"), .f(sys, "ps"), .f(sys, "iv"))
  # defaults of parameters used only by events
  md <- .discover_metadata(extra, .f(sys, "initial_conditions"), .f(sys, "bindings"), .f(sys, "guesses"))
  .update_system(sys, continuous_events = ce, discrete_events = de,
                 ps = .unique_syms(c(.f(sys, "ps"), extra)), initial_conditions = md$initial_conditions,
                 bindings = md$bindings, guesses = md$guesses)
}

# Overridden by the events module; identity defaults.
.normalize_events <- function(events, kind) {
  h <- .mtk$hooks[[paste0("normalize_", kind, "_events")]]
  if (!is.null(h)) return(h(events))
  as.list(events)
}
.event_parameters <- function(events, unknowns, ps, iv) {
  h <- .mtk$hooks$event_parameters
  if (!is.null(h)) return(h(events, unknowns, ps, iv))
  list()
}

.flatten_vars <- function(x) {
  if (is.null(x)) return(list())
  if (is_sym(x)) return(list(x))
  out <- list()
  for (el in x) out <- c(out, .flatten_vars(el))
  out
}

# Split declared variables found in equations into unknowns, parameters and
# brownians. Derivatives contribute their variable as an unknown.
.classify_vars <- function(found, iv) {
  unknowns <- list(); ps <- list(); br <- list(); po <- list()
  ivk <- if (is.null(iv)) "" else .key(iv)
  for (v in found) {
    if (.type(v) == "Derivative") v <- .args(v)[[1L]]
    if (.key(v) == ivk) next
    k <- .var_kind(v)
    if (is.null(k)) next
    if (k == "operator" || k == "iv") next
    if (k == "unknown") {
      if (!is.null(iv) && .is_delay_or_evalat(v, iv)) next
      unknowns[[length(unknowns) + 1L]] <- v
    } else if (k %in% c("parameter", "constant", "discrete")) {
      ps[[length(ps) + 1L]] <- v
    } else if (k == "brownian") {
      br[[length(br) + 1L]] <- v
    } else if (k == "poissonian") {
      po[[length(po) + 1L]] <- v
    }
  }
  list(unknowns = .unique_syms(unknowns), ps = .unique_syms(ps),
       brownians = .unique_syms(br), poissonians = .unique_syms(po))
}

# Variables that appear only in the symbolic defaults of found variables are
# collected too, recursively, as collect_var! does in ModelingToolkit.
.with_default_vars <- function(found) {
  seen <- .keys(found)
  queue <- found
  while (length(queue)) {
    v <- queue[[1L]]
    queue <- queue[-1L]
    if (.type(v) == "Derivative") v <- .args(v)[[1L]]
    d <- .meta_get(v, "default")
    if (!is_sym(d)) next
    for (w in .collect_vars(d)) {
      if (.key(w) %in% seen) next
      seen <- c(seen, .key(w))
      found[[length(found) + 1L]] <- w
      queue[[length(queue) + 1L]] <- w
    }
  }
  found
}

# x(t - 1) or x(0) is a delayed or evaluated form, not an unknown itself.
.is_delay_or_evalat <- function(v, iv) {
  .type(v) == "Function" && length(.args(v)) >= 1L &&
    !any(vapply(.args(v), function(a) identical(.key(a), .key(iv)), TRUE))
}

.check_dde <- function(eqs, iv) {
  if (is.null(iv)) return(FALSE)
  for (e in eqs) {
    if (!inherits(e, "Equation")) next
    for (v in .collect_vars(list(e$lhs, e$rhs))) {
      if (.type(v) == "Function" && identical(.var_kind(v), "unknown") && .is_delay_or_evalat(v, iv)) {
        a <- .args(v)[[1L]]
        if (!.is_num(a)) return(TRUE)
      }
    }
  }
  FALSE
}

# Defaults on variables become initial conditions (numbers) or bindings
# (expressions or NA), as collect_defaults! in MTK. Namespaced variables are
# skipped: their defaults were collected by their own subsystem.
.discover_metadata <- function(vars, ics, binds, guesses) {
  for (v in vars) {
    if (grepl(NAMESPACE_SEPARATOR, .name(v) %||% "", fixed = TRUE)) next
    if (.type(v) == "Derivative") next
    d <- .meta_get(v, "default")
    if (!is.null(d) && !smap_has(ics, v) && !smap_has(binds, v)) {
      if (is.numeric(d) && length(d) == 1L && !is.na(d)) ics <- smap_set(ics, v, d)
      else if (.is_num(d)) ics <- smap_set(ics, v, as.numeric(d))
      else if (is_sym(d) || (is.atomic(d) && length(d) == 1L && is.na(d))) binds <- smap_set(binds, v, d)
      # a nonnumeric value (a function, a string, a list, a matrix) is a
      # constant, an initial condition
      else ics <- smap_set(ics, v, d)
    }
    g <- .meta_get(v, "guess")
    if (!is.null(g) && !smap_has(guesses, v)) guesses <- smap_set(guesses, v, g)
  }
  list(initial_conditions = ics, bindings = binds, guesses = guesses)
}

# Scopes of leaves, read from their metadata attributes.
.collect_scopes <- function(eqs) {
  out <- list()
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    sc <- attr(z, "mtk_meta")$scope
    if (!is.null(sc)) out[[.key(z)]] <<- sc
    for (a in .args(z)) walk(a)
  }
  for (e in eqs) {
    if (!inherits(e, "Equation")) next
    walk(e$lhs)
    walk(e$rhs)
  }
  out
}

.check_system <- function(sys) {
  iv <- .f(sys, "iv")
  # check_variables: an undeclared symbol in a time dependent system is an
  # unknown that does not depend on time, which ModelingToolkit rejects
  for (u in .f(sys, "unknowns")) {
    if (.type(u) == "Symbol" && is.null(attr(u, "mtk_meta")) &&
        !exists(.name(u), envir = .mtk$vars, inherits = FALSE))
      stop("variable ", format(u), " is not a function of the independent variable ", format(iv),
           "; declare parameters with Parameter() or parameters(), and unknowns as ",
           "Variable(\"", format(u), "\", ", format(iv), ")", call. = FALSE)
  }
  for (e in .f(sys, "eqs")) {
    if (!inherits(e, "Equation")) next
    for (d in Filter(.is_derivative, .collect_vars(list(e$lhs, e$rhs)))) {
      di <- .deriv_info(d)
      if (di$n_specs > 1L || !identical(.key(di$iv), .key(iv))) {
        if (!inherits(sys, "PDESystem"))
          stop("equation ", format(e), " differentiates with respect to ", format(di$iv),
               ", not the independent variable ", format(iv), call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Field accessors (get_*), which return the system's own fields.

#' Fields of a system
#'
#' `get_eqs(sys)`, `get_unknowns(sys)` and the other `get_*` functions return
#' a system's own fields, without the contributions of subsystems. The
#' functions without the prefix, such as [equations()], include subsystems.
#'
#' @param sys A `System`.
#' @return The field value.
#' @name system_fields
NULL

#' @rdname system_fields
#' @export
get_eqs <- function(sys) .f(sys, "eqs")
#' @rdname system_fields
#' @export
get_unknowns <- function(sys) .f(sys, "unknowns")
#' @rdname system_fields
#' @export
get_ps <- function(sys) .f(sys, "ps")
#' @rdname system_fields
#' @export
get_iv <- function(sys) .f(sys, "iv")
#' @rdname system_fields
#' @export
get_observed <- function(sys) .f(sys, "observed")
#' @rdname system_fields
#' @export
get_systems <- function(sys) .f(sys, "systems")
#' @rdname system_fields
#' @export
get_name <- function(sys) .f(sys, "name")
#' @rdname system_fields
#' @export
get_description <- function(sys) .f(sys, "description")
#' @rdname system_fields
#' @export
get_bindings <- function(sys) .f(sys, "bindings")
#' @rdname system_fields
#' @export
get_initial_conditions <- function(sys) .f(sys, "initial_conditions")
#' @rdname system_fields
#' @export
get_guesses <- function(sys) .f(sys, "guesses")
#' @rdname system_fields
#' @export
get_initialization_eqs <- function(sys) .f(sys, "initialization_eqs")
#' @rdname system_fields
#' @export
get_continuous_events <- function(sys) .f(sys, "continuous_events")
#' @rdname system_fields
#' @export
get_discrete_events <- function(sys) .f(sys, "discrete_events")
#' @rdname system_fields
#' @export
get_connector_type <- function(sys) .f(sys, "connector_type")
#' @rdname system_fields
#' @export
get_metadata <- function(sys) .f(sys, "metadata")
#' @rdname system_fields
#' @export
get_noise_eqs <- function(sys) .f(sys, "noise_eqs")
#' @rdname system_fields
#' @export
get_jumps <- function(sys) .f(sys, "jumps")
#' @rdname system_fields
#' @export
get_constraints <- function(sys) .f(sys, "constraints")
#' @rdname system_fields
#' @export
get_costs <- function(sys) .f(sys, "costs")
#' @rdname system_fields
#' @export
get_brownians <- function(sys) .f(sys, "brownians")
#' @rdname system_fields
#' @export
get_parent <- function(sys) .f(sys, "parent")
#' @rdname system_fields
#' @export
get_tspan <- function(sys) .f(sys, "tspan")
#' @rdname system_fields
#' @export
get_tstops <- function(sys) .f(sys, "tstops")
#' @rdname system_fields
#' @export
get_schedule <- function(sys) .f(sys, "schedule")
#' @rdname system_fields
#' @export
get_index_cache <- function(sys) .f(sys, "index_cache")
#' @rdname system_fields
#' @export
get_inputs <- function(sys) .f(sys, "inputs")
#' @rdname system_fields
#' @export
get_outputs <- function(sys) .f(sys, "outputs")
#' @rdname system_fields
#' @export
get_assertions <- function(sys) .f(sys, "assertions")
#' @rdname system_fields
#' @export
get_tearing_state <- function(sys) .f(sys, "tearing_state")
#' @rdname system_fields
#' @export
get_initializesystem <- function(sys) .f(sys, "initializesystem")

#' Name of a system
#'
#' @param sys A `System`.
#' @return The name, a string.
#' @export
nameof <- function(sys) .f(sys, "name")

#' Whether a system is complete
#'
#' @param sys A `System`.
#' @return A logical.
#' @export
iscomplete <- function(sys) isTRUE(.f(sys, "complete"))

# a ModelingToolkit system, not the multibody System of symcompR
.is_mtk_system <- function(x) inherits(x, "System") && !inherits(x, "physics_obj")

# Names modelcompR shares with symcompR. A call meant for symcompR is handed
# to it exactly as written: the caller's call, with its argument names, order
# and dots, evaluated again in the caller's frame with symcompR's function.
# Older symcompR versions lack some of the names.
.symcompr_has <- function(name) name %in% getNamespaceExports("symcompR")

.symcompr_forward <- function(name, call, env) {
  call[[1L]] <- call("::", as.name("symcompR"), as.name(name))
  eval(call, env)
}

#' @export
`$.System` <- function(x, name) {
  # the multibody System of symcompR shares the class name
  if (inherits(x, "physics_obj")) return(NextMethod())
  .getproperty(x, name)
}

#' @export
`[[.System` <- function(x, i) {
  if (inherits(x, "physics_obj")) return(NextMethod())
  if (is.character(i) && length(i) == 1L) return(.getproperty(x, i))
  stop("index a System by name", call. = FALSE)
}

#' @export
names.System <- function(x) {
  if (inherits(x, "physics_obj")) return(NextMethod())
  c(vapply(.f(x, "systems"), nameof, ""),
    vapply(c(.f(x, "unknowns"), .f(x, "ps"), lapply(.f(x, "observed"), `[[`, "lhs")),
           function(v) .name(v) %||% format(v), ""))
}

# getproperty: a subsystem, then a variable of this system, then a field.
.getproperty <- function(sys, name) {
  # complete and compiled systems look variables up in their hierarchical
  # parent, as getproperty does through get_parent in ModelingToolkit
  par <- .f(sys, "parent")
  if (!is.null(par) && !isTRUE(.f(sys, "namespacing"))) {
    hit <- tryCatch(.getproperty(toggle_namespacing(par, FALSE), name), error = function(e) NULL)
    if (!is.null(hit)) return(hit)
  }
  ns <- isTRUE(.f(sys, "namespacing"))
  for (s in .f(sys, "systems")) {
    if (identical(nameof(s), name)) {
      if (ns) s <- .rename_system(s, paste0(nameof(sys), NAMESPACE_SEPARATOR, nameof(s)))
      return(s)
    }
  }
  cands <- c(.f(sys, "unknowns"), .f(sys, "ps"), lapply(.f(sys, "observed"), `[[`, "lhs"),
             .f(sys, "brownians"))
  for (v in cands) {
    if (identical(.name(v), name)) return(if (ns) renamespace(sys, v) else v)
  }
  # Analysis points created in this system (getvar in ModelingToolkit).
  for (e in .f(sys, "eqs")) {
    if (inherits(e, "AnalysisPoint") && identical(e$name, name)) return(if (ns) renamespace(sys, e) else e)
  }
  # The independent variable, never namespaced (getvar in ModelingToolkit).
  iv <- .f(sys, "iv")
  if (!is.null(iv) && identical(.name(iv), name)) return(iv)
  # A namespaced variable of a subsystem, written with the separator.
  if (grepl(NAMESPACE_SEPARATOR, name, fixed = TRUE)) {
    parts <- strsplit(name, NAMESPACE_SEPARATOR, fixed = TRUE)[[1L]]
    obj <- sys
    for (p in parts) obj <- .getproperty(obj, p)
    return(obj)
  }
  if (name %in% .SYSTEM_FIELDS) return(.f(sys, name))
  stop("variable ", name, " does not exist in system ", nameof(sys), call. = FALSE)
}

# ---------------------------------------------------------------------------
# Namespacing

.ns_name <- function(prefix, nm) paste0(prefix, NAMESPACE_SEPARATOR, nm)

# Rename one variable leaf (Symbol or applied Function) into `prefix`,
# honoring its scope. Returns the leaf unchanged for global scope.
.rename_leaf <- function(v, prefix, scopes = NULL) {
  kind <- .var_kind(v)
  if (is.null(kind) || kind %in% c("iv", "operator")) return(v)
  sc <- scopes[[.key(v)]] %||% attr(v, "mtk_meta")$scope
  if (!is.null(sc)) {
    if (identical(sc$type, "global")) return(v)
    if (identical(sc$type, "parent")) {
      m <- attr(v, "mtk_meta")
      if (sc$depth <= 1L) m$scope <- NULL else m$scope <- list(type = "parent", depth = sc$depth - 1L)
      attr(v, "mtk_meta") <- m
      return(v)
    }
  }
  newname <- .ns_name(prefix, .name(v))
  ent <- .mtk$vars[[.name(v)]]
  ent <- ent$forms[[.type(v)]] %||% ent
  meta <- attr(v, "mtk_meta")
  if (is.null(meta) && !is.null(ent)) meta <- c(list(kind = ent$kind), ent$meta)
  nent <- .mtk$vars[[newname]]
  if (is.null(nent) || (!is.null(nent$forms) && is.null(nent$forms[[.type(v)]]))) {
    m <- meta
    m$kind <- NULL
    .register_var(newname, kind, m, if (.type(v) == "Function") .args(v) else NULL)
  }
  out <- if (.type(v) == "Function") do.call(Function(newname), .args(v)) else Symbol(newname)
  if (!is.null(meta)) {
    meta$scope <- NULL
    attr(out, "mtk_meta") <- meta
  }
  out
}

#' Namespace a variable or expression
#'
#' `renamespace(sys, x)` prefixes the variables of `x` with the name of
#' `sys`, as they appear when `sys` is a subsystem: `x(t)` in a subsystem
#' named `r1` becomes `r1₊x(t)`. `sys` may also be a name.
#'
#' @param sys A `System` or a name.
#' @param x A variable, expression, equation, or list of them.
#' @return The namespaced object.
#' @export
renamespace <- function(sys, x) {
  prefix <- if (inherits(sys, "System")) nameof(sys) else as.character(sys)
  scopes <- if (inherits(sys, "System")) .f(sys, "scopes") else NULL
  .namespace(x, prefix, scopes)
}

.namespace <- function(x, prefix, scopes = NULL) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "Equation")) return(Equation(.namespace(x$lhs, prefix, scopes), .namespace(x$rhs, prefix, scopes)))
  if (inherits(x, "Inequality")) {
    x$lhs <- .namespace(x$lhs, prefix, scopes)
    x$rhs <- .namespace(x$rhs, prefix, scopes)
    return(x)
  }
  if (inherits(x, "Connection")) return(.namespace_connection(x, prefix))
  if (inherits(x, "System")) return(.rename_system(x, .ns_name(prefix, nameof(x))))
  if (inherits(x, "mtk_event") && !is.null(.mtk$hooks$namespace_event)) return(.mtk$hooks$namespace_event(x, prefix, scopes))
  if (!is_sym(x)) {
    if (is.list(x)) {
      out <- lapply(x, .namespace, prefix = prefix, scopes = scopes)
      if (inherits(x, "AbstractJump")) attributes(out) <- attributes(x)
      return(out)
    }
    return(x)
  }
  if (.type(x) %in% c("Symbol", "Function") && .is_declared(x) &&
      !identical(.var_kind(x), "operator")) {
    if (.type(x) == "Function") {
      # rename the head; arguments such as t stay
      return(.rename_leaf(x, prefix, scopes))
    }
    return(.rename_leaf(x, prefix, scopes))
  }
  found <- .collect_vars(x)
  # derivative nodes are namespaced through their variable
  leaves <- .unique_syms(lapply(found, function(v) if (.type(v) == "Derivative") .args(v)[[1L]] else v))
  leaves <- Filter(function(v) !(.var_kind(v) %in% c("iv", "operator")), leaves)
  if (!length(leaves)) return(x)
  to <- lapply(leaves, .rename_leaf, prefix = prefix, scopes = scopes)
  .subst(x, .rule(leaves, to))
}

# Hook for connections; replaced by the connectors module.
.namespace_connection <- function(x, prefix) {
  h <- .mtk$hooks$namespace_connection
  if (is.null(h)) return(x)
  h(x, prefix)
}

.rename_system <- function(sys, newname) .update_system(sys, name = newname)

#' Namespace the equations of a system
#'
#' @param sys A `System`.
#' @return A list of equations with every variable namespaced by the system
#'   name.
#' @export
namespace_equations <- function(sys) lapply(equations(sys), function(e) renamespace(sys, e))

# ---------------------------------------------------------------------------
# Hierarchical accessors

# Apply f to a system and recursively to its subsystems, namespacing results
# of subsystems into this system.
.collect_hier <- function(sys, f) {
  own <- f(sys)
  subs <- list()
  for (s in .f(sys, "systems")) {
    part <- .collect_hier(s, f)
    subs <- c(subs, lapply(part, function(el) renamespace(s, el)))
  }
  c(own, subs)
}

#' Equations, unknowns and parameters of a system
#'
#' These accessors include the contributions of all subsystems, namespaced
#' by the subsystem names, as the functions of the same names in
#' ModelingToolkit. `*_toplevel` variants return only the system's own
#' entries. For a system returned by [mtkcompile()], `unknowns()` are the
#' states that remain after simplification and `observed()` gives the
#' eliminated variables.
#'
#' @param sys A `System`.
#' @return A list.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)"); parameters("k")
#' sys <- System(list(D(x) ~ -k * x), t, name = "decay")
#' equations(sys); unknowns(sys); parameters(sys)
#' @export
equations <- function(sys) .collect_hier(sys, function(s) .f(s, "eqs"))

#' @rdname equations
#' @export
equations_toplevel <- function(sys) .f(sys, "eqs")

#' @rdname equations
#' @export
unknowns <- function(sys) .unique_syms(.collect_hier(sys, function(s) .f(s, "unknowns")))

#' @rdname equations
#' @export
unknowns_toplevel <- function(sys) .f(sys, "unknowns")

.parameters_of <- function(sys, initial_parameters = FALSE) {
  ps <- .unique_syms(.collect_hier(sys, function(s) .f(s, "ps")))
  if (!initial_parameters) ps <- Filter(function(p) !isinitial(p), ps)
  ps
}

#' @rdname equations
#' @export
parameters_toplevel <- function(sys) .f(sys, "ps")

#' @rdname equations
#' @export
observed <- function(sys) .collect_hier(sys, function(s) .f(s, "observed"))

#' @rdname equations
#' @export
observables <- function(sys) lapply(observed(sys), `[[`, "lhs")

#' @rdname equations
#' @export
independent_variables <- function(sys) {
  iv <- .f(sys, "iv")
  if (is.null(iv)) list() else list(iv)
}

#' @rdname equations
#' @export
initialization_equations <- function(sys) .collect_hier(sys, function(s) .f(s, "initialization_eqs"))

#' @rdname equations
#' @export
full_equations <- function(sys) {
  obs <- observed(sys)
  eqs <- equations(sys)
  if (!length(obs)) return(eqs)
  rule <- .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs"))
  lapply(eqs, function(e) Equation(e$lhs, .subst_fixpoint(e$rhs, rule)))
}

#' @rdname equations
#' @export
costs <- function(sys) .collect_hier(sys, function(s) .f(s, "costs"))

#' @rdname equations
#' @export
cost <- function(sys) {
  cs <- costs(sys)
  cons <- .f(sys, "consolidate")
  if (!is.null(cons)) return(cons(cs))
  if (!length(cs)) return(.S0())
  Reduce(`+`, cs)
}

#' @rdname equations
#' @export
constraints <- function(sys) .collect_hier(sys, function(s) .f(s, "constraints"))

#' @rdname equations
#' @export
jumps <- function(sys) .collect_hier(sys, function(s) .f(s, "jumps"))

#' @rdname equations
#' @export
inputs <- function(sys) .f(sys, "inputs")

#' @rdname equations
#' @export
outputs <- function(sys) .f(sys, "outputs")

#' @rdname equations
#' @export
assertions <- function(sys) .f(sys, "assertions")

.hier_map <- function(sys, field) {
  out <- .f(sys, field) %||% symmap()
  for (s in .f(sys, "systems")) {
    sub <- .hier_map(s, field)
    for (el in unclass(sub)) {
      if (is.null(el$var)) next
      out <- smap_set(out, renamespace(s, el$var),
                      if (is_sym(el$value)) renamespace(s, el$value) else el$value)
    }
  }
  # Own entries override those of subsystems.
  own <- .f(sys, field)
  if (length(own)) for (k in names(own)) out[[k]] <- own[[k]]
  class(out) <- "mtk_symmap"
  out
}

#' Initial conditions, bindings and guesses of a system
#'
#' Hierarchical maps from variables to values. `defaults()` merges initial
#' conditions and bindings, as the older ModelingToolkit API did.
#'
#' @param sys A `System`.
#' @return A variable map.
#' @export
initial_conditions <- function(sys) .hier_map(sys, "initial_conditions")
#' @rdname initial_conditions
#' @export
bindings <- function(sys) .hier_map(sys, "bindings")
#' @rdname initial_conditions
#' @export
guesses <- function(sys) .hier_map(sys, "guesses")
#' @rdname initial_conditions
#' @export
defaults <- function(sys) smap_merge(initial_conditions(sys), bindings(sys))

#' @rdname initial_conditions
#' @param ... Variable maps to add, as formulas `x ~ 1` or pairs.
#' @export
set_defaults <- function(sys, ...) {
  m <- as_symmap(list(...))
  .update_system(sys, initial_conditions = smap_merge(.f(sys, "initial_conditions"), m))
}

#' @rdname initial_conditions
#' @export
bound_parameters <- function(sys) {
  b <- bindings(sys)
  ps <- .parameters_of(sys)
  Filter(function(p) smap_has(b, p), ps)
}

#' @rdname equations
#' @export
diff_equations <- function(sys) Filter(is_diff_equation, equations(sys))
#' @rdname equations
#' @export
alg_equations <- function(sys) Filter(is_alg_equation, equations(sys))
#' @rdname equations
#' @export
has_diff_equations <- function(sys) length(diff_equations(sys)) > 0L
#' @rdname equations
#' @export
has_alg_equations <- function(sys) length(alg_equations(sys)) > 0L
#' @rdname equations
#' @export
get_diff_eqs <- function(sys) Filter(is_diff_equation, .f(sys, "eqs"))
#' @rdname equations
#' @export
get_alg_eqs <- function(sys) Filter(is_alg_equation, .f(sys, "eqs"))
#' @rdname equations
#' @export
has_diff_eqs <- function(sys) length(get_diff_eqs(sys)) > 0L
#' @rdname equations
#' @export
has_alg_eqs <- function(sys) length(get_alg_eqs(sys)) > 0L

# ---------------------------------------------------------------------------
# Building systems from systems

#' Compose, extend, flatten and complete systems
#'
#' `compose(sys, subsystems)` adds subsystems; `extend(sys, base)` merges the
#' equations, variables and maps of `sys` into `base`; `flatten(sys)`
#' collapses the hierarchy into one system whose variables are namespaced;
#' `complete(sys)` marks a system as finished so its variables are accessed
#' without the top level namespace and problems can be built from it.
#'
#' `flatten()`, `compose()` and `extend()` also accept non-system arguments
#' and then behave as symcompR's functions of the same names (`extend()` of a
#' symcompR `DimensionSystem` or `UnitSystem`), so attaching both packages is
#' safe.
#'
#' @param sys A `System`.
#' @param ... Subsystems (for `compose`), or arguments passed on.
#' @param base The system to extend.
#' @param name Name of the result.
#' @param flatten Whether `complete()` flattens the system first.
#' @return A `System`.
#' @export
compose <- function(sys, ..., name = NULL) {
  if (missing(sys) || !.is_mtk_system(sys)) return(.symcompr_forward("compose", sys.call(), parent.frame()))
  if (is.null(name)) name <- nameof(sys)
  subs <- list(...)
  if (length(subs) == 1L && is.list(subs[[1L]]) && !inherits(subs[[1L]], "System")) subs <- subs[[1L]]
  .update_system(sys, systems = c(.f(sys, "systems"), subs), name = as.character(name))
}

#' @rdname compose
#' @export
extend <- function(sys, base, name = nameof(sys), ...) {
  # symcompR's extend() generic (unit and dimension systems) gets the call as
  # written
  if ((missing(sys) || !.is_mtk_system(sys)) && .symcompr_has("extend"))
    return(.symcompr_forward("extend", sys.call(), parent.frame()))
  if (iscomplete(base) || iscomplete(sys)) stop("cannot extend a complete system", call. = FALSE)
  iv <- .f(base, "iv") %||% .f(sys, "iv")
  out <- System(
    eqs = c(.f(base, "eqs"), .f(sys, "eqs")), iv = iv,
    unknowns = .unique_syms(c(.f(base, "unknowns"), .f(sys, "unknowns"))),
    ps = .unique_syms(c(.f(base, "ps"), .f(sys, "ps"))),
    brownians = .unique_syms(c(.f(base, "brownians"), .f(sys, "brownians"))),
    name = name,
    systems = c(.f(base, "systems"), .f(sys, "systems")),
    observed = c(.f(base, "observed"), .f(sys, "observed")),
    bindings = smap_merge(.f(base, "bindings"), .f(sys, "bindings")),
    initial_conditions = smap_merge(.f(base, "initial_conditions"), .f(sys, "initial_conditions")),
    guesses = smap_merge(.f(base, "guesses"), .f(sys, "guesses")),
    initialization_eqs = c(.f(base, "initialization_eqs"), .f(sys, "initialization_eqs")),
    constraints = c(.f(base, "constraints"), .f(sys, "constraints")),
    costs = c(.f(base, "costs"), .f(sys, "costs")),
    jumps = c(.f(base, "jumps"), .f(sys, "jumps")),
    connector_type = .f(sys, "connector_type") %||% .f(base, "connector_type"),
    description = if (nzchar(.f(sys, "description"))) .f(sys, "description") else .f(base, "description"),
    metadata = utils::modifyList(.f(base, "metadata"), .f(sys, "metadata")),
    tspan = .f(sys, "tspan") %||% .f(base, "tspan"),
    discover_from_metadata = FALSE, checks = FALSE)
  .update_system(out,
                 continuous_events = c(.f(base, "continuous_events"), .f(sys, "continuous_events")),
                 discrete_events = c(.f(base, "discrete_events"), .f(sys, "discrete_events")),
                 inputs = .unique_syms(c(.f(base, "inputs"), .f(sys, "inputs"))),
                 outputs = .unique_syms(c(.f(base, "outputs"), .f(sys, "outputs"))),
                 scopes = c(.f(base, "scopes"), .f(sys, "scopes")))
}

#' @rdname compose
#' @export
flatten <- function(sys, ...) {
  if (missing(sys) || !.is_mtk_system(sys)) return(.symcompr_forward("flatten", sys.call(), parent.frame()))
  if (!length(.f(sys, "systems"))) return(sys)
  .update_system(sys,
    eqs = equations(sys),
    unknowns = unknowns(sys),
    ps = .parameters_of(sys, initial_parameters = TRUE),
    brownians = .unique_syms(.collect_hier(sys, function(s) .f(s, "brownians"))),
    poissonians = .unique_syms(.collect_hier(sys, function(s) .f(s, "poissonians"))),
    observed = observed(sys),
    bindings = bindings(sys),
    initial_conditions = initial_conditions(sys),
    guesses = guesses(sys),
    initialization_eqs = initialization_equations(sys),
    continuous_events = .collect_hier(sys, function(s) .f(s, "continuous_events")),
    discrete_events = .collect_hier(sys, function(s) .f(s, "discrete_events")),
    constraints = constraints(sys),
    costs = costs(sys),
    jumps = jumps(sys),
    # Julia's flatten builds a new System without noise_eqs: the noise
    # equations of the root and of the subsystems are dropped
    noise_eqs = NULL,
    assertions = .collect_hier(sys, function(s) .f(s, "assertions")),
    state_priorities = .hier_map(sys, "state_priorities"),
    irreducibles = .unique_syms(.collect_hier(sys, function(s) .f(s, "irreducibles"))),
    systems = list())
}

#' @rdname compose
#' @export
complete <- function(sys, flatten = TRUE) {
  if (flatten) {
    newsys <- expand_connections(sys)
    newsys <- flatten(newsys)
    # keep the hierarchical system for variable access (sys$sub$x)
    if (is.null(.f(newsys, "parent")))
      newsys <- .update_system(newsys, parent = complete(sys, flatten = FALSE))
    sys <- newsys
    h <- .mtk$hooks$complete
    if (!is.null(h)) sys <- h(sys)
  }
  sys <- .update_system(sys, complete = TRUE, namespacing = FALSE)
  .update_system(sys, index_cache = .build_index_cache(sys))
}

#' @rdname compose
#' @param value Whether variables are namespaced on access.
#' @export
toggle_namespacing <- function(sys, value) .update_system(sys, namespacing = isTRUE(value))

#' @rdname compose
#' @export
hierarchy <- function(sys) {
  pr <- function(s, depth) {
    cat(strrep("  ", depth), nameof(s), "\n", sep = "")
    for (c in .f(s, "systems")) pr(c, depth + 1L)
  }
  pr(sys, 0L)
  invisible(sys)
}

#' Rename a system
#'
#' @param sys A `System`.
#' @param name The new name.
#' @return The renamed system.
#' @export
rename_system <- function(sys, name) .rename_system(sys, as.character(name))

#' @export
print.System <- function(x, ...) {
  if (inherits(x, "physics_obj")) return(NextMethod())
  iv <- .f(x, "iv")
  kind <- if (length(.f(x, "costs"))) "Optimization system"
          else if (is.null(iv)) "Nonlinear system"
          else if (!is.null(.f(x, "noise_eqs")) || length(.f(x, "brownians"))) "Stochastic system"
          else if (length(.f(x, "jumps"))) "Jump system"
          else if (isTRUE(.f(x, "is_discrete"))) "Discrete system"
          else "Model"
  eqs <- equations(x)
  cat(kind, " ", nameof(x), if (iscomplete(x)) " (complete)", ":\n", sep = "")
  if (nzchar(.f(x, "description"))) cat(.f(x, "description"), "\n")
  cat("Equations (", length(eqs), "):\n", sep = "")
  show <- utils::head(eqs, 12L)
  for (e in show) cat("  ", .fmt(e), "\n", sep = "")
  if (length(eqs) > 12L) cat("  ...\n")
  u <- unknowns(x)
  cat("Unknowns (", length(u), "): ", paste(vapply(utils::head(u, 10L), format, ""), collapse = ", "),
      if (length(u) > 10L) ", ...", "\n", sep = "")
  p <- .parameters_of(x)
  cat("Parameters (", length(p), "): ", paste(vapply(utils::head(p, 10L), format, ""), collapse = ", "),
      if (length(p) > 10L) ", ...", "\n", sep = "")
  ob <- observed(x)
  if (length(ob)) cat("Observed (", length(ob), "): ",
                      paste(vapply(utils::head(ob, 6L), function(e) format(e$lhs), ""), collapse = ", "),
                      if (length(ob) > 6L) ", ...", "\n", sep = "")
  if (length(.f(x, "systems")))
    cat("Subsystems (", length(.f(x, "systems")), "): ",
        paste(vapply(.f(x, "systems"), nameof, ""), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Legacy system constructors
#'
#' ModelingToolkit 11 replaced `ODESystem`, `NonlinearSystem` and the other
#' system types with [System()]. These wrappers call `System()` and remain
#' for compatibility.
#'
#' @param eqs,iv,unknowns,ps,... Passed to [System()].
#' @return A `System`.
#' @export
ODESystem <- function(eqs, iv, unknowns = NULL, ps = NULL, ...) System(eqs, iv, unknowns, ps, ...)
#' @rdname ODESystem
#' @export
NonlinearSystem <- function(eqs, unknowns = NULL, ps = NULL, ...) System(eqs, NULL, unknowns, ps, ...)
#' @rdname ODESystem
#' @export
SDESystem <- function(eqs, noise_eqs, iv, unknowns = NULL, ps = NULL, ...)
  System(eqs, iv, unknowns, ps, noise_eqs = noise_eqs, ...)
#' @rdname ODESystem
#' @param op Objective expression of an optimization system.
#' @param noise_eqs Noise coefficients of a stochastic system.
#' @param jumps Jumps of a jump system.
#' @export
OptimizationSystem <- function(op, unknowns = NULL, ps = NULL, ...) System(list(), NULL, unknowns, ps, costs = list(op), ...)
#' @rdname ODESystem
#' @export
JumpSystem <- function(jumps, iv, unknowns = NULL, ps = NULL, ...) {
  if (inherits(jumps, c("AbstractJump", "Equation", "formula"))) jumps <- list(jumps)
  is_eq <- vapply(jumps, function(j) inherits(j, c("Equation", "formula")), TRUE)
  System(jumps[is_eq], iv, unknowns, ps, jumps = jumps[!is_eq], ...)
}
#' @rdname ODESystem
#' @export
DiscreteSystem <- function(eqs, iv, unknowns = NULL, ps = NULL, ...) System(eqs, iv, unknowns, ps, is_discrete = TRUE, ...)
#' @rdname ODESystem
#' @export
ImplicitDiscreteSystem <- function(eqs, iv, unknowns = NULL, ps = NULL, ...) System(eqs, iv, unknowns, ps, is_discrete = TRUE, ...)
