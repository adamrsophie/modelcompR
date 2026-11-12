# Ported from StateSelection.jl (https://github.com/JuliaComputing/StateSelection.jl),
# commit 29b3c90. Portions are licensed under the GNU AGPL, version 3; see
# inst/COPYRIGHTS and inst/STATESELECTION-AGPL-LICENSE.
#
# State machines, after
# StateSelection.jl/lib/ModelingToolkitTearing/src/clock_inference/
# state_machines.jl (InitialState, Transition, initial_state, transition,
# activeState, entry, ticksInState, timeInState, namespace_expr) and
# extract_top_level_statemachines and remove_child_equations of
# ModelingToolkit.jl/src/systems/systemstructure.jl.
#
# ModelingToolkit 11 defines the operators, extracts the top level state
# machines before flattening and puts their equations through clock
# inference; simulating them needs an external compilation pass. modelcompR
# supplies one as an extension, following the synchronous state machine
# semantics of the Modelica specification (3.6, chapter 17): every machine is
# a clocked partition and becomes a periodic callback that, at each tick of
# its clock,
#   1. selects the state chosen at the previous tick (the initial state at
#      the first tick),
#   2. fires the immediate transition of highest priority (smallest
#      `priority`) out of it whose condition holds, unless a delayed
#      transition pending from the previous tick has a higher priority,
#   3. resets the variables of the new state (and their past values) to
#      their start values when the fired transition has `reset = TRUE`,
#   4. evaluates the equations of the active state; the variables of the
#      other states keep their values,
#   5. records the delayed transition of highest priority out of the active
#      state whose condition holds, which fires at the next tick, and shifts
#      the past values of all variables of the machine.
# Conditions of immediate transitions see the values before the tick, with
# the counters of the source state already advanced to this tick (at the
# first tick the initial state is being entered, so entry() holds); those of
# delayed transitions see the values after the active state's equations. In
# a transition condition ticksInState() and timeInState() refer to the source
# state. ticksInState() is 1 at the tick a state is entered and timeInState()
# is 0 then; entry() is true at that tick. Hierarchical machines and
# `synchronize = TRUE` are not supported.

# ---------------------------------------------------------------------------
# Operators

.sm_marker <- function(op) .op_node(op, list(0L))

.sm_equation <- function(op, data, rhs = .S0()) {
  structure(list(lhs = .sm_marker(op), rhs = as_sym(rhs), op = op, data = data),
            class = c("StateMachineEquation", "Equation"))
}

.sm_state_name <- function(s) {
  nm <- if (inherits(s, "System")) nameof(s) else as.character(s)
  segs <- strsplit(nm, NAMESPACE_SEPARATOR, fixed = TRUE)[[1L]]
  segs[length(segs)]
}

#' State machines
#'
#' A state machine is a `System` whose equations are all `initial_state()`
#' and `transition()` statements and whose subsystems are its states: systems
#' of clocked discrete equations (with [ShiftIndex()]). `initial_state(s)`
#' names the state active at the first tick and `transition(from, to, cond)`
#' moves from state `from` to state `to` when the condition `cond` holds.
#' With `immediate = TRUE` the transition fires at the tick where `cond`
#' holds, before the equations of `from` are evaluated; otherwise it fires at
#' the next tick. `reset = TRUE` restarts the target state from the start
#' values of its variables. Among the transitions leaving a state the one
#' with the smallest `priority` wins.
#'
#' In the equations of a state, and in transition conditions,
#' `ticksInState()` is the number of ticks since the state was entered
#' (1 at the entry tick), `timeInState()` the time since then, `entry()` is
#' true at the entry tick and `activeState(s)` is true when state `s` is
#' active. With a state argument they refer to that state; in a transition
#' condition, without one, to the source state.
#'
#' ModelingToolkit 11 extracts state machines and infers their clocks, but
#' simulating them needs an external compiler. [mtkcompile()] compiles each
#' machine of a hybrid model into a periodic callback with the semantics of
#' the Modelica specification (chapter 17); this is an extension of
#' modelcompR. Hierarchical machines and `synchronize = TRUE` are not
#' supported.
#'
#' @param s,from,to States (systems), or their names.
#' @param cond The transition condition, a symbolic boolean such as
#'   `ticksInState() >= 3`.
#' @param immediate,reset,synchronize,priority Transition options.
#' @return `initial_state()` and `transition()` return state machine
#'   equations; the other functions return symbolic expressions.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' k <- ShiftIndex(Clock(1))
#' x <- GlobalScope(Variable("x", t, default = 0))
#' up <- System(list(k(x) ~ (k - 1)(x) + 1), t, list(x), list(), name = "up")
#' down <- System(list(k(x) ~ (k - 1)(x) - 1), t, list(x), list(), name = "down")
#' sm <- System(list(initial_state(up), transition(up, down, ticksInState() >= 3),
#'   transition(down, up, ticksInState() >= 3)), t, systems = list(up, down), name = "sm")
#' z <- Variable("z", t, default = 0)
#' model <- System(list(D(z) ~ Hold(x)), t, systems = list(sm), name = "model")
#' sol <- solve(ODEProblem(mtkcompile(model), list(), c(0, 8)))
#' sol(seq(0.5, 7.5, by = 1), idxs = x)
#' @export
initial_state <- function(s) .sm_equation("InitialState", list(s = .sm_state_name(s)))

#' @rdname initial_state
#' @export
transition <- function(from, to, cond, immediate = TRUE, reset = TRUE, synchronize = FALSE, priority = 1L) {
  .sm_equation("Transition", list(from = .sm_state_name(from), to = .sm_state_name(to),
                                  immediate = isTRUE(immediate), reset = isTRUE(reset),
                                  synchronize = isTRUE(synchronize), priority = as.integer(priority)),
               rhs = cond)
}

.sm_op_symbol <- function(op, s) {
  nm <- if (is.null(s)) paste0(op, "()") else paste0(op, "(", .sm_state_name(s), ")")
  if (!exists(nm, envir = .mtk$vars, inherits = FALSE)) .register_var(nm, "operator")
  Symbol(nm)
}

#' @rdname initial_state
#' @export
activeState <- function(s) .sm_op_symbol("activeState", s)

#' @rdname initial_state
#' @export
entry <- function(s = NULL) .sm_op_symbol("entry", s)

#' @rdname initial_state
#' @export
ticksInState <- function(s = NULL) .sm_op_symbol("ticksInState", s)

#' @rdname initial_state
#' @export
timeInState <- function(s = NULL) .sm_op_symbol("timeInState", s)

#' @export
format.StateMachineEquation <- function(x, ...) {
  d <- x$data
  if (x$op == "InitialState") return(paste0("initial_state(", d$s, ")"))
  sprintf("%s \u2192 %s if (%s) [immediate: %d, reset: %d, sync: %d, prio: %d]", d$from, d$to,
          format(x$rhs), d$immediate, d$reset, d$synchronize, d$priority)
}

#' @export
print.StateMachineEquation <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Extraction (extract_top_level_statemachines, remove_child_equations)

.sm_remove_child_equations <- function(sys) {
  .update_system(sys, eqs = list(), systems = lapply(.f(sys, "systems"), .sm_remove_child_equations))
}

# Returns the system without the equations of its state machines and the
# machines, each with the systems between the root and it (`path`).
.sm_extract <- function(sys, path = list()) {
  eqs <- .f(sys, "eqs")
  is_sm <- vapply(eqs, inherits, TRUE, "StateMachineEquation")
  if (length(eqs) && all(is_sm))
    return(list(sys = .sm_remove_child_equations(sys), machines = list(list(sys = sys, path = path))))
  if (any(is_sm))
    stop("Mixing statemachine equations and standard equations in a top-level statemachine is not allowed.",
         call. = FALSE)
  subs <- .f(sys, "systems")
  machines <- list()
  for (i in seq_along(subs)) {
    r <- .sm_extract(subs[[i]], c(path, list(subs[[i]])))
    subs[[i]] <- r$sys
    machines <- c(machines, r$machines)
  }
  if (length(machines)) sys <- .update_system(sys, systems = subs)
  list(sys = sys, machines = machines)
}

# A machine in the namespace of the flattened root: its states' equations,
# its transitions and its initial state.
.sm_prepare <- function(m) {
  smsys <- m$sys
  if (!length(m$path))
    stop("a state machine must be a subsystem of a model with continuous equations", call. = FALSE)
  ancestors <- m$path[-length(m$path)]
  up <- function(x) {
    for (a in rev(ancestors)) x <- renamespace(a, x)
    x
  }
  ns <- function(x) up(renamespace(smsys, x))
  states <- .f(smsys, "systems")
  snames <- vapply(states, nameof, "")
  init <- character(); trans <- list()
  for (eq in .f(smsys, "eqs")) {
    if (eq$op == "InitialState") {
      init <- c(init, eq$data$s)
      next
    }
    d <- eq$data
    if (d$synchronize) stop("transition(synchronize = TRUE) is not supported", call. = FALSE)
    for (s in c(d$from, d$to)) if (!(s %in% snames))
      stop("transition refers to ", s, ", which is not a state of state machine ", nameof(smsys), call. = FALSE)
    trans[[length(trans) + 1L]] <- c(d[c("immediate", "reset", "priority")],
                                     list(from = match(d$from, snames), to = match(d$to, snames), cond = ns(eq$rhs)))
  }
  if (length(init) != 1L || !(init %in% snames))
    stop("state machine ", nameof(smsys), " needs exactly one initial_state() among its states", call. = FALSE)
  for (i in seq_along(snames)) {
    pr <- vapply(Filter(function(tr) tr$from == i, trans), `[[`, 0L, "priority")
    if (anyDuplicated(pr)) stop("the transitions out of state ", snames[i], " need distinct priorities", call. = FALSE)
  }
  state_eqs <- lapply(states, function(s) {
    eqs <- equations(s)
    if (any(vapply(eqs, inherits, TRUE, "StateMachineEquation")))
      stop("hierarchical state machines are not supported", call. = FALSE)
    ns(renamespace(s, eqs))
  })
  list(name = nameof(smsys), states = snames, initial = match(init, snames),
       transitions = trans[order(vapply(trans, `[[`, 0L, "priority"))], state_eqs = state_eqs)
}

# ---------------------------------------------------------------------------
# The compiled machine (an extension; see the header)

# The operator symbols seen by the equations of state `i` (or a condition of
# a transition out of it): their values in the order of .sm_op_values().
.sm_op_vars <- function(states) {
  own <- lapply(c("entry", "ticksInState", "timeInState"), .sm_op_symbol, s = NULL)
  named <- unlist(lapply(states, function(s) lapply(c("entry", "ticksInState", "timeInState", "activeState"),
                                                    .sm_op_symbol, s = s)), recursive = FALSE)
  c(own, named)
}

.sm_op_values <- function(st, i, t) {
  n <- length(st$ticks)
  per <- unlist(lapply(seq_len(n), function(j) c(as.numeric(st$entry[j]), st$ticks[j], t - st$entered[j],
                                                 as.numeric(identical(st$active, j)))))
  c(as.numeric(st$entry[i]), st$ticks[i], t - st$entered[i], per)
}

# Compile a machine partition: `pt` holds the lowered equations of each state
# (pt$low_eqs), the variables each state defines (pt$state_vars), the past
# value chains (pt$hist) and the clock.
.sm_callback <- function(pt, ctx) {
  opv <- .sm_op_vars(pt$states)
  args <- function(z) list(list(name = "u", vars = ctx$us), list(name = "p", vars = ctx$ps),
                           list(name = "sm", vars = opv), list(name = "z", vars = z),
                           list(name = "t", vars = ctx$iv))
  state_fns <- lapply(seq_along(pt$states), function(i) {
    ep <- .explicit_partition(pt$low_eqs[[i]], pt$state_vars[[i]])
    a <- args(ep$newton)
    obs <- c(ep$assigns, ctx$obs)
    list(val = .build_closure(pt$state_vars[[i]], a, observed = obs, out = "vector"),
         res = if (length(ep$residuals)) .build_closure(ep$residuals, a, observed = obs, out = "vector"),
         z_idx = match(.keys(ep$newton), ctx$pk), idx = match(.keys(pt$state_vars[[i]]), ctx$pk))
  })
  cond_fns <- lapply(pt$transitions, function(tr)
    .build_closure(tr$cond, args(list()), observed = ctx$obs, out = "scalar"))
  hist_cur <- list(); hist_old <- list()
  for (h in pt$hist) {
    hist_cur <- c(hist_cur, list(h$var))
    hist_old <- c(hist_old, list(h$vars))
  }
  all_vars <- .unique_syms(c(unlist(pt$state_vars, recursive = FALSE), hist_cur, unlist(hist_old, recursive = FALSE)))
  for (v in all_vars) if (!(.key(v) %in% ctx$pk))
    stop("state machine variable ", format(v), " is not a parameter of the system", call. = FALSE)
  all_idx <- match(.keys(all_vars), ctx$pk)
  hist_idx <- lapply(seq_along(hist_cur), function(i) list(cur = match(.key(hist_cur[[i]]), ctx$pk),
                                                       old = match(.keys(hist_old[[i]]), ctx$pk)))
  n <- length(pt$states)
  # the variables (and past values) a reset of state i restores
  reset_idx <- lapply(seq_len(n), function(i) {
    ks <- .keys(pt$state_vars[[i]])
    hs <- unlist(lapply(pt$hist, function(h) if (.key(h$var) %in% ks) .keys(h$vars)))
    match(c(ks, hs), ctx$pk)
  })
  st <- new.env(parent = emptyenv())
  holds <- function(k, i, u, p, t) isTRUE(as.logical(cond_fns[[k]](u, p, .sm_op_values(st, i, t), numeric(), t)))
  enter <- function(i, t) {
    st$entry <- logical(n)
    st$entry[i] <- TRUE
    st$ticks[i] <- 1
    st$entered[i] <- t
  }
  aff <- function(integ) {
    u <- integ$u; t <- integ$t
    first <- is.null(st$active)
    # the selected state at this tick: entered now at the first tick,
    # otherwise one tick older
    selected <- if (first) pt$initial else st$nxt
    if (first) {
      enter(selected, t)
    } else {
      st$entry <- logical(n)
      st$ticks[selected] <- st$ticks[selected] + 1
    }
    fired <- NULL
    for (k in seq_along(pt$transitions)) {
      tr <- pt$transitions[[k]]
      if (tr$immediate && tr$from == selected && holds(k, selected, u, integ$p, t)) {
        fired <- k
        break
      }
    }
    if (!is.null(st$delayed) && (is.null(fired) || st$delayed < fired)) fired <- st$delayed
    active <- selected
    if (!is.null(fired)) {
      active <- pt$transitions[[fired]]$to
      if (pt$transitions[[fired]]$reset) integ$p[reset_idx[[active]]] <- st$start[match(reset_idx[[active]], all_idx)]
      enter(active, t)
    }
    st$active <- active
    sf <- state_fns[[active]]
    p <- integ$p
    smv <- .sm_op_values(st, active, t)
    z <- numeric()
    if (!is.null(sf$res)) {
      sol <- .newton_solve(function(zz) sf$res(u, p, smv, zz, t), p[sf$z_idx])
      if (!sol$converged) stop("state ", pt$states[active], " could not be solved at t = ", t, call. = FALSE)
      z <- sol$x
    }
    if (length(sf$idx)) integ$p[sf$idx] <- sf$val(u, p, smv, z, t)
    # past values of every variable of the machine move by one tick
    for (h in hist_idx) {
      old <- integ$p[h$old]
      integ$p[h$old] <- c(integ$p[h$cur], old[-length(old)])
    }
    st$delayed <- NULL
    for (k in seq_along(pt$transitions)) {
      tr <- pt$transitions[[k]]
      if (!tr$immediate && tr$from == active && holds(k, active, u, integ$p, t)) {
        st$delayed <- k
        break
      }
    }
    st$nxt <- active
    invisible(NULL)
  }
  clk <- pt$clock
  cb <- .periodic_event_callback(aff, clk$dt, clk$phase)
  own <- cb$initialize
  cb$initialize <- function(integ) {
    st$active <- NULL; st$delayed <- NULL; st$nxt <- NULL
    st$ticks <- numeric(n); st$entered <- numeric(n); st$entry <- logical(n)
    st$start <- integ$p[all_idx]
    if (!is.null(own)) own(integ)
  }
  cb
}

# The machine partitions of a hybrid model: called by .compile_hybrid() with
# the machines prepared by .sm_prepare(), the clock of each machine and the
# independent variable. Returns the partitions and the variables that become
# parameters (with their past values).
.sm_partitions <- function(machines, clock_ids, iv) {
  parts <- list(); vars <- list(); hist <- list()
  for (i in seq_along(machines)) {
    m <- machines[[i]]
    clk <- .clock_of(clock_ids[i])
    if (is.null(clk) || !identical(clk$type, "periodic") || is.null(clk$dt))
      stop("state machine ", m$name, " needs a periodic clock with a sample interval", call. = FALSE)
    eqs <- lapply(unlist(m$state_eqs, recursive = FALSE), function(e) .strip_sample(e))
    eqs <- .substitute_sample_time(eqs, clock_ids[i])
    lhs_vars <- function(es) .unique_syms(Filter(function(v) .is_declared(v) && !isparameter(v),
                                                 lapply(es, `[[`, "lhs")))
    defined <- lhs_vars(eqs)
    low <- .lower_shifts(eqs, defined, iv)
    counts <- lengths(m$state_eqs)
    ends <- cumsum(counts)
    low_eqs <- lapply(seq_along(counts), function(j) low$eqs[seq_len(counts[j]) + ends[j] - counts[j]])
    state_vars <- lapply(low_eqs, lhs_vars)
    parts[[i]] <- list(type = "statemachine", clock = clk, states = m$states, initial = m$initial,
                       transitions = m$transitions, low_eqs = low_eqs, state_vars = state_vars,
                       hist = low$hist)
    vars <- c(vars, defined)
    for (h in low$hist) hist <- c(hist, h$vars)
  }
  list(parts = parts, vars = .unique_syms(vars), hist = .unique_syms(hist))
}
