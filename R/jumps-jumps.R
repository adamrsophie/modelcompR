# Jumps and jump problems, after ModelingToolkitBase/src/problems/jumpproblem.jl
# (JumpProblem, JumpSysMajParamMapper, collect_vars! for jumps,
# modified_unknowns!), systems/codegen.jl (assemble_maj, assemble_crj,
# assemble_vrj, numericrstoich), systems/system.jl (SymbolicMassActionJump)
# and the JumpProcesses.jl types and solvers they build: jumps.jl
# (MassActionJump, ConstantRateJump, VariableRateJump), problem.jl
# (JumpProblem), SSA_stepper.jl (SSAStepper, its step!, solve! and
# savevalues!), variable_rate.jl (VR_FRM) and aggregators/aggregated_api.jl
# (reset_aggregated_jumps!). The aggregators are in R/jumps-aggregators.R.
#
# Every jump stores a symbolic `rate` (for mass action jumps the full
# propensity, rate constant times the falling factorials of the reactants)
# and `affect`, a list of equations `x ~ Pre(x) + ...`. Pure jump problems
# are simulated by SSAStepper with an aggregator. Jumps coupled to ODEs or
# SDEs are handled as in JumpProcesses: the aggregator handles the mass
# action and constant rate jumps through a discrete callback that stops the
# integrator at each jump time, and each variable rate jump integrates its
# rate in an extra state that starts at -Exp(1) and fires the jump when it
# crosses zero (VR_FRM).

# ---------------------------------------------------------------------------
# Jump types

# Stoichiometry as list(list(var = <sym or name>, stoich = <number>), ...).
.jmp_stoich <- function(x) {
  if (is.null(x) || !length(x)) return(list())
  if (is.numeric(x) && !is.null(names(x)))
    return(lapply(names(x), function(n) list(var = n, stoich = unname(x[[n]]))))
  if (inherits(x, c("Equation", "formula"))) x <- list(x)
  lapply(x, function(el) {
    if (inherits(el, "formula")) el <- as_equation(el)
    if (inherits(el, "Equation")) return(list(var = el$lhs, stoich = as.numeric(el$rhs)))
    # an unknown by its index, as the Pair `i => s` of JumpProcesses
    if (is.numeric(el) && is.null(names(el)) && length(el) == 2L) return(list(var = as.integer(el[1L]), stoich = el[2L]))
    if (is.list(el) && length(el) == 2L) return(list(var = as_sym(el[[1L]]), stoich = as.numeric(el[[2L]])))
    stop("cannot interpret a stoichiometry entry; use X ~ 1, list(X, 1) or c(X = 1)", call. = FALSE)
  })
}

.jmp_has_names <- function(st) any(vapply(st, function(s) is.character(s$var) || is.integer(s$var), TRUE))

# Mass action propensity rate * prod_s X_s (X_s - 1) ... (X_s - n_s + 1), as
# evalrxrate in JumpProcesses with pre-scaled rates.
.jmp_mass_action_rate <- function(rate, rs) {
  out <- as_sym(rate)
  for (s in rs) {
    if (.is_num(as_sym(s$var)) && as.numeric(s$var) == 0) next
    for (k in seq_len(s$stoich)) out <- out * (s$var - (k - 1))
  }
  out
}

.jmp_finish_maj <- function(j) {
  if (.jmp_has_names(j$reactant_stoch) || .jmp_has_names(j$net_stoch)) return(j)
  j$rate <- .jmp_mass_action_rate(j$scaled_rates, j$reactant_stoch)
  j$affect <- lapply(j$net_stoch, function(s) Equation(s$var, Pre(s$var) + s$stoich))
  j
}

#' Jumps
#'
#' Jumps of a jump process, as in JumpProcesses.jl, for use in a [System()]
#' (`jumps = list(...)`) or [JumpSystem()]. `ConstantRateJump(rate, affect)`
#' has a rate that changes only when a jump fires; `VariableRateJump(rate,
#' affect)` a rate that may depend on time or on unknowns that change
#' continuously. `affect` is a list of equations giving the new values, such
#' as `list(S ~ Pre(S) - 1, I ~ Pre(I) + 1)`; inside it, `Pre(x)` (or plain
#' `x`) is the value just before the jump.
#'
#' `MassActionJump(rate, reactant_stoch, net_stoch)` is a mass action
#' reaction with propensity `rate * prod(X * (X - 1) * ... )` over the
#' reactants; stoichiometries are lists of formulas `list(S ~ 1, I ~ 1)` (or
#' pairs, or named numeric vectors). As in ModelingToolkit the rate must be
#' pre-scaled (`k / factorial(n)` for `n X -> ...`):
#' `SymbolicMassActionJump()` and `scale_rates = FALSE` do this, while the
#' JumpProcesses default `scale_rates = TRUE` is rejected by systems.
#'
#' @param rate The rate: a symbolic expression, or for problems without a
#'   system an R function `rate(u, p, t)`.
#' @param affect A list of equations, or an R function `affect(integ)`
#'   changing `integ$u`.
#' @param reactant_stoch,net_stoch Reactant and net stoichiometries.
#' @param scale_rates Divide the rate by the factorials of the reactant
#'   stoichiometries.
#' @param save_positions Save before and after the jump (a variable rate
#'   jump coupled to an ODE or SDE); the JumpProcesses default saves after it.
#' @param bounds For [RSSA()] and [RSSACR()] with problems without a system:
#'   a function `bounds(ulow, uhigh, u, p, t)` returning `c(lrate, urate)`,
#'   bounds of the rate while each population stays in `[ulow, uhigh]`. By
#'   default the rate is evaluated at `ulow` and `uhigh`, which assumes it is
#'   monotone in the state.
#' @param lrate,urate,rateinterval Rate bounds of a bounded variable rate jump
#'   (R functions of `(u, p, t)`): for `s` in `[t, t + rateinterval(u, p, t)]`
#'   the rate lies in `[lrate(u, p, t), urate(u, p, t)]`. [Coevolve()]
#'   simulates bounded jumps of problems without a system by thinning, also
#'   in pure jump problems. ModelingToolkit drops these bounds for symbolic
#'   jumps, and so does `JumpProblem()` for systems.
#' @param rootfind,interp_points,abstol,reltol Options of the continuous
#'   callback that fires a variable rate jump coupled to an ODE or SDE;
#'   `rootfind = TRUE` locates the firing time (LeftRootFind). The root
#'   finder of the R driver has its own sampling and tolerances, so
#'   `interp_points`, `abstol` and `reltol` are kept for compatibility only.
#' @return A jump object.
#' @examples
#' t <- t_nounits
#' variables("S(t)", "I(t)", "R(t)"); parameters("beta", "gamma")
#' j1 <- ConstantRateJump(beta * S * I, list(S ~ Pre(S) - 1, I ~ Pre(I) + 1))
#' j2 <- SymbolicMassActionJump(gamma, list(I ~ 1), list(I ~ -1, R ~ 1))
#' js <- JumpSystem(list(j1, j2), t, list(S, I, R), list(beta, gamma), name = "sir")
#' unknowns(js)
#' @name jumps_types
NULL

#' @rdname jumps_types
#' @export
ConstantRateJump <- function(rate, affect, bounds = NULL) {
  if (!is.function(rate)) rate <- as_sym(rate)
  if (!is.function(affect)) affect <- .as_equations(affect)
  structure(list(rate = rate, affect = affect, bounds = bounds), class = c("ConstantRateJump", "AbstractJump"))
}

#' @rdname jumps_types
#' @export
VariableRateJump <- function(rate, affect, lrate = NULL, urate = NULL, rateinterval = NULL,
                             rootfind = TRUE, save_positions = c(FALSE, TRUE), interp_points = 10L,
                             abstol = 1e-12, reltol = 0) {
  if (is.null(urate) != is.null(rateinterval))
    stop("`urate` and `rateinterval` must both be `nothing`, or must both be defined.", call. = FALSE)
  if (!is.null(lrate) && is.null(urate))
    stop("If a lower bound rate, `lrate`, is given than an upper bound rate, `urate`, and rate interval, `rateinterval`, must also be provided.",
         call. = FALSE)
  if (!is.function(rate)) rate <- as_sym(rate)
  if (!is.function(affect)) affect <- .as_equations(affect)
  structure(list(rate = rate, affect = affect, lrate = lrate, urate = urate, rateinterval = rateinterval,
                 rootfind = rootfind, save_positions = save_positions, interp_points = interp_points,
                 abstol = abstol, reltol = reltol),
            class = c("VariableRateJump", "AbstractJump"))
}

#' @rdname jumps_types
#' @export
MassActionJump <- function(rate, reactant_stoch, net_stoch, scale_rates = TRUE) {
  # several reactions at once, as JumpProcesses' MassActionJump(rates,
  # reactant_stoch, net_stoch) with one entry per reaction
  if (is.numeric(rate) && length(rate) > 1L) {
    if (length(reactant_stoch) != length(rate) || length(net_stoch) != length(rate))
      stop("give one reactant and one net stoichiometry per rate", call. = FALSE)
    return(structure(lapply(seq_along(rate), function(k)
      MassActionJump(rate[k], reactant_stoch[[k]], net_stoch[[k]], scale_rates = scale_rates)),
      class = "MassActionJumpSet"))
  }
  rs <- .jmp_stoich(reactant_stoch)
  ns <- .jmp_stoich(net_stoch)
  rate <- as_sym(rate)
  if (scale_rates) {
    fac <- prod(vapply(rs, function(s) factorial(s$stoich), 0))
    if (fac != 1) rate <- rate / fac
  }
  .jmp_finish_maj(structure(list(scaled_rates = rate, reactant_stoch = rs, net_stoch = ns,
                                 scale_rates = scale_rates, rescale_rates_on_update = scale_rates),
                            class = c("MassActionJump", "AbstractJump")))
}

#' @rdname jumps_types
#' @export
SymbolicMassActionJump <- function(rate, reactant_stoch, net_stoch, scale_rates = FALSE) {
  if (scale_rates)
    stop("SymbolicMassActionJump requires pre-scaled rate expressions (scale_rates = FALSE). scale_rates = TRUE is not supported in ModelingToolkitBase.",
         call. = FALSE)
  MassActionJump(rate, reactant_stoch, net_stoch, scale_rates = FALSE)
}

#' @export
print.AbstractJump <- function(x, ...) {
  cat(class(x)[1L], ": rate ", if (is.function(x$rate)) "<function>" else format(x$rate), "\n", sep = "")
  if (is.list(x$affect)) for (e in x$affect) cat("  ", .fmt(e), "\n", sep = "")
  invisible(x)
}

# Expressions whose variables a jump contributes to its system.
.jmp_exprs <- function(j) {
  out <- list()
  if (is_sym(j$rate)) out <- c(out, list(j$rate))
  if (!is.null(j$scaled_rates)) out <- c(out, list(j$scaled_rates))
  if (is.list(j$affect)) for (e in j$affect) out <- c(out, list(e$lhs, e$rhs))
  for (s in c(j$reactant_stoch, j$net_stoch)) if (is_sym(s$var)) out <- c(out, list(s$var))
  out
}

# The jump_exprs hook of System(): variables of jumps and poissonian rates.
.mtk$hooks$jump_exprs <- function(jumps, found) {
  if (inherits(jumps, "AbstractJump")) jumps <- list(jumps)
  out <- list()
  for (j in jumps) {
    if (!inherits(j, "AbstractJump")) next
    if (inherits(j, "MassActionJump") && isTRUE(j$rescale_rates_on_update))
      stop("MassActionJump with rescale_rates_on_update = true is not supported in Systems with jumps or JumpSystems. Rate expressions must be pre-scaled (e.g. k/factorial(n) for n-th order reactions). Use SymbolicMassActionJump or pass scale_rates = false when constructing the MassActionJump.",
           call. = FALSE)
    out <- c(out, .jmp_exprs(j))
  }
  for (v in found) {
    if (identical(.var_kind(v), "poissonian")) {
      r <- .meta_get(v, "rate")
      if (!is.null(r)) out <- c(out, list(as_sym(r)))
    }
  }
  out
}

#' Unknowns changed by a jump
#'
#' `modified_unknowns(jump, sts)` returns the unknowns among `sts` that the
#' jump changes, as `modified_unknowns!` in ModelingToolkit.
#'
#' @param jump A jump.
#' @param sts A list of unknowns.
#' @return A list of variables.
#' @examples
#' t <- t_nounits
#' variables("A(t)", "B(t)")
#' j <- ConstantRateJump(1, list(A ~ Pre(A) + 1))
#' modified_unknowns(j, list(A, B))
#' @export
modified_unknowns <- function(jump, sts) {
  lhs <- if (inherits(jump, "MassActionJump")) lapply(jump$net_stoch, `[[`, "var")
         else lapply(jump$affect, `[[`, "lhs")
  Filter(function(s) .contains_key(lhs, s), sts)
}

# ---------------------------------------------------------------------------
# Assembly: numeric rates and affects

.jmp_unpre <- function(e) {
  pres <- list()
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    if (.is_op(z, "Pre")) {
      pres[[length(pres) + 1L]] <<- z
      return(invisible())
    }
    for (a in .args(z)) walk(a)
  }
  walk(e)
  if (!length(pres)) return(e)
  pres <- .unique_syms(pres)
  .subst(e, .rule(pres, lapply(pres, function(z) .args(z)[[1L]])))
}

# Resolve stoichiometries given by name against the system's unknowns.
.jmp_resolve <- function(j, us) {
  if (!inherits(j, "MassActionJump")) return(j)
  by_name <- stats::setNames(us, vapply(us, function(u) .name(u) %||% format(u), ""))
  fix <- function(st) lapply(st, function(s) {
    if (is.character(s$var)) {
      v <- by_name[[s$var]]
      if (is.null(v)) .stopf("no unknown named %s for a mass action jump", s$var)
      s$var <- v
    }
    s
  })
  j$reactant_stoch <- fix(j$reactant_stoch)
  j$net_stoch <- fix(j$net_stoch)
  .jmp_finish_maj(j)
}


# ---------------------------------------------------------------------------
# Numeric jumps

# Rate and affect closures of a constant or variable rate jump of a system
# (generate_rate_function, generate_affect_function). The state may carry
# the integrated intensities of variable rate jumps after the n unknowns.
.jmp_rate_fn <- function(j, args, obs)
  .build_closure(list(.jmp_unpre(j$rate)), args, observed = obs, out = "scalar")

.jmp_affect_fn <- function(j, ukeys, args, obs) {
  lhs <- lapply(j$affect, `[[`, "lhs")
  idx <- match(.keys(lhs), ukeys)
  if (anyNA(idx))
    .stopf("jump affects may only change unknowns; %s is not an unknown", format(lhs[[which(is.na(idx))[1L]]]))
  fn <- .build_closure(lapply(j$affect, function(e) .jmp_unpre(e$rhs)), args, observed = obs, out = "vector")
  n <- length(ukeys)
  function(integ) {
    x <- integ$u
    if (length(x) != n) x <- x[seq_len(n)]
    integ$u[idx] <- fn(x, integ$p, integ$t)
  }
}

.jmp_vrj_entry <- function(j, rate, affect)
  list(rate = rate, affect = affect, save_positions = j$save_positions %||% c(FALSE, TRUE),
       rootfind = j$rootfind %||% TRUE)

# The jumps of a system in JumpProcesses' order: mass action jumps
# (assemble_maj: numeric stoichiometries, and the rate constants as a function
# of the parameters as JumpSysMajParamMapper), constant rate jumps
# (assemble_crj), then variable rate jumps (assemble_vrj, which drops rate
# bounds). `jumps` keeps the first two kinds for the dependency graphs.
.jmp_numeric <- function(sys, js) {
  us <- unknowns(sys)
  ukeys <- .keys(us)
  js <- lapply(js, .jmp_resolve, us = us)
  majs <- Filter(function(j) inherits(j, "MassActionJump"), js)
  crjs <- Filter(function(j) inherits(j, "ConstantRateJump"), js)
  vrjs <- Filter(function(j) inherits(j, "VariableRateJump"), js)
  args <- .codegen_args(sys)
  obs <- .codegen_observed(sys)
  stoich <- function(st) {
    st <- Filter(function(s) !(.is_num(as_sym(s$var)) && as.numeric(s$var) == 0), st)
    idx <- match(.keys(lapply(st, function(s) as_sym(s$var))), ukeys)
    if (anyNA(idx)) stop("a mass action jump refers to a variable that is not an unknown", call. = FALSE)
    cbind(idx, vapply(st, function(s) as.numeric(s$stoich), 0))
  }
  .jmp_compat(list(
    n = length(us), nmaj = length(majs), ncrj = length(crjs), nbvrj = 0L,
    maj_rs = lapply(majs, function(j) stoich(j$reactant_stoch)),
    maj_ns = lapply(majs, function(j) stoich(j$net_stoch)),
    maj_ratefun = if (length(majs))
      .build_closure(lapply(majs, function(j) as_sym(j$scaled_rates)), args, observed = obs, out = "vector"),
    crj_rate = lapply(crjs, .jmp_rate_fn, args = args, obs = obs),
    crj_affect = lapply(crjs, .jmp_affect_fn, ukeys = ukeys, args = args, obs = obs),
    crj_bounds = rep(list(NULL), length(crjs)),
    bvrj = list(),
    vrj = lapply(vrjs, function(j) .jmp_vrj_entry(j, .jmp_rate_fn(j, args, obs),
                                                   .jmp_affect_fn(j, ukeys, args, obs))),
    jumps = c(majs, crjs)))
}

# Jumps with R function rates and affects, for problems without a system.
.jmp_numeric_raw <- function(js, n) {
  # mass action jumps of a problem without a system refer to unknowns by
  # index: MassActionJump(rate, list(c(i, s), ...), list(c(i, s), ...))
  idx_stoich <- function(st) {
    if (!length(st)) return(matrix(integer(), 0L, 2L))
    if (!all(vapply(st, function(s) is.integer(s$var), TRUE)))
      stop("mass action jumps of a problem without a system refer to unknowns by index, ",
           "for example list(c(1, 1)) for one molecule of the first unknown", call. = FALSE)
    cbind(vapply(st, `[[`, 0L, "var"), vapply(st, function(s) as.numeric(s$stoich), 0))
  }
  majs <- Filter(function(j) inherits(j, "MassActionJump"), js)
  js <- Filter(function(j) !inherits(j, c("MassActionJump", "RegularJump")), js)
  for (j in js) {
    if (!is.function(j$rate) || !is.function(j$affect))
      stop("jumps of a problem without a system need R functions rate(u, p, t) and affect(integ)", call. = FALSE)
  }
  crjs <- Filter(function(j) inherits(j, "ConstantRateJump"), js)
  vrjs <- Filter(function(j) inherits(j, "VariableRateJump"), js)
  bounded <- vapply(vrjs, function(j) is.function(j$urate) && is.function(j$rateinterval), TRUE)
  maj_rates <- vapply(majs, function(j) as.numeric(j$scaled_rates), 0)
  .jmp_compat(list(
    n = n, nmaj = length(majs), ncrj = length(crjs), nbvrj = sum(bounded),
    maj_rs = lapply(majs, function(j) idx_stoich(j$reactant_stoch)),
    maj_ns = lapply(majs, function(j) idx_stoich(j$net_stoch)),
    maj_ratefun = if (length(majs)) function(u, p, t) maj_rates,
    crj_rate = lapply(crjs, `[[`, "rate"), crj_affect = lapply(crjs, `[[`, "affect"),
    crj_bounds = lapply(crjs, function(j) j$bounds),
    bvrj = lapply(vrjs[bounded], function(j) list(rate = j$rate, urate = j$urate, lrate = j$lrate,
                                                  rateinterval = j$rateinterval, affect = j$affect,
                                                  jump = j)),
    vrj = lapply(vrjs[!bounded], function(j) .jmp_vrj_entry(j, j$rate, j$affect)),
    jumps = c(majs, crjs)))
}

# Bounded variable rate jumps go to the aggregator only if it handles them
# (use_vrj_bounds and supports_variablerates); otherwise they are general.
.jmp_partition <- function(x, agg, use_vrj_bounds) {
  if (!x$nbvrj || (use_vrj_bounds && .supports_variablerates(agg))) return(x)
  x$vrj <- c(lapply(x$bvrj, function(b) .jmp_vrj_entry(b$jump, b$rate, b$affect)), x$vrj)
  x$bvrj <- list()
  x$nbvrj <- 0L
  x
}

# rates(u, p, t), all rates of the mass action and constant rate jumps, and
# affect(k, u, p, t), the state after jump k
.jmp_compat <- function(x) {
  x$rates <- function(u, p, t) {
    mr <- if (x$nmaj) as.numeric(x$maj_ratefun(u, p, t)) else numeric()
    c(vapply(seq_len(x$nmaj), function(k) .evalrxrate(u, x$maj_rs[[k]], mr[k]), 0),
      vapply(x$crj_rate, function(f) f(u, p, t), 0))
  }
  x$affect <- function(k, u, p, t) {
    integ <- new.env(parent = emptyenv())
    integ$u <- u; integ$p <- p; integ$t <- t
    if (k <= x$nmaj) {
      ns <- x$maj_ns[[k]]
      integ$u[ns[, 1L]] <- integ$u[ns[, 1L]] + ns[, 2L]
    } else x$crj_affect[[k - x$nmaj]](integ)
    integ$u
  }
  x
}

# The affects of the system's events reset the aggregated jumps
# (process_events with reset_jumps = true).
.jmp_reset_wrap <- function(cbs) {
  if (is.null(cbs)) return(NULL)
  cbs <- .as_callback_set(cbs)
  wrap <- function(f) {
    if (is.null(f)) return(NULL)
    force(f)
    function(integ, ...) {
      f(integ, ...)
      reset_aggregated_jumps(integ)
    }
  }
  fix <- function(cb) {
    cb$affect <- wrap(cb$affect)
    if (!is.null(cb$affect_neg)) cb$affect_neg <- wrap(cb$affect_neg)
    cb$initialize <- wrap(cb$initialize)
    cb
  }
  cbs$continuous <- lapply(cbs$continuous, fix)
  cbs$discrete <- lapply(cbs$discrete, fix)
  cbs
}

# ---------------------------------------------------------------------------
# JumpProblem

#' Jump problems
#'
#' `JumpProblem(sys, op, tspan)` builds a jump problem from a complete system
#' with jumps. Without equations and variable rate jumps the inner problem
#' `jprob$prob` is a `DiscreteProblem`, simulated exactly by [SSAStepper()];
#' otherwise it is an `ODEProblem` (or an `SDEProblem` for systems with
#' noise). As in JumpProcesses, the mass action and constant rate jumps are
#' handled by the `aggregator` (see [jump_aggregators]; `NullAggregator()`
#' chooses one), which stops the ODE or SDE solver at each of their jump
#' times, and every variable rate jump fires when its integrated rate crosses
#' an exponential threshold. Events of the system reset the aggregator after
#' their affects.
#'
#' `JumpProblem(prob, aggregator, jump1, jump2, ...)` couples jumps with R
#' function rates and affects to an existing problem (`prob` a
#' `DiscreteProblem` for pure jump processes).
#'
#' `solve(jprob)` uses [SSAStepper()] for pure jump problems, [Tsit5()] for
#' jumps coupled to ODEs and [SOSRI()] for jumps coupled to SDEs. Options are
#' those of the underlying solver, plus `seed`. `SSAStepper()` takes
#' `saveat`, `save_start`, `save_end`, `tstops` and discrete callbacks
#' (checked after every jump and at every stop), and saves the state after a
#' jump when any of `save_positions` is `TRUE`. For jumps coupled to ODEs
#' and SDEs, `save_positions` controls the saves around each jump as for
#' callbacks in DiffEqBase, also when `saveat` is given.
#'
#' @param sys A complete `System` with jumps, or a problem.
#' @param ... For systems: `op` and `tspan`, then options stored with the
#'   problem (and `bracket_data`, see [BracketData()]). For problems: the
#'   aggregator and the jumps.
#' @param aggregator A jump aggregator.
#' @param save_positions Save the state before and after each jump of the
#'   aggregator; by default both for systems and coupled problems, and after
#'   the jump for `DiscreteProblem`s.
#' @param dependencies For problems without a system: a list with the
#'   dependency graph `dep_graph`, `vartojumps_map` and `jumptovars_map`
#'   (lists of integer vectors) when the aggregator needs them.
#' @param use_vrj_bounds Let an aggregator that supports them ([Coevolve()])
#'   handle bounded variable rate jumps.
#' @param vr_aggregator How general variable rate jumps are simulated, see
#'   [variable_rate_aggregators].
#' @param hopping_constants,spatial_system For spatial problems with [NSM()]
#'   or [DirectCRDirect()], see [spatial_jumps].
#' @return A `JumpProblem`.
#' @examples
#' t <- t_nounits
#' variables("X(t)"); parameters("lam = 10", "mu = 1")
#' birth <- SymbolicMassActionJump(lam, list(), list(X ~ 1))
#' death <- SymbolicMassActionJump(mu, list(X ~ 1), list(X ~ -1))
#' js <- complete(JumpSystem(list(birth, death), t, list(X), list(lam, mu), name = "bd"))
#' jprob <- JumpProblem(js, list(X ~ 0), c(0, 5), aggregator = RSSA())
#' sol <- solve(jprob, SSAStepper(), seed = 1)
#' tail(sol[X], 1)
#' @export
JumpProblem <- function(sys, ..., aggregator = NullAggregator(), save_positions = NULL,
                        dependencies = NULL, use_vrj_bounds = TRUE, hopping_constants = NULL,
                        spatial_system = NULL, vr_aggregator = VR_FRM()) {
  if (!inherits(vr_aggregator, "VariableRateAggregator"))
    stop("`vr_aggregator` must be VR_FRM(), VR_Direct() or VR_DirectFW()", call. = FALSE)
  out <- if (inherits(sys, "SciMLProblem"))
    .jmp_problem_raw(sys, ..., aggregator = aggregator, save_positions = save_positions,
                     dependencies = dependencies, use_vrj_bounds = use_vrj_bounds,
                     hopping_constants = hopping_constants, spatial_system = spatial_system)
  else .jmp_problem_sys(sys, ..., aggregator = aggregator, save_positions = save_positions %||% c(TRUE, TRUE))
  out$vr_aggregator <- vr_aggregator
  out
}

#' Variable rate aggregators
#'
#' How [JumpProblem()] simulates general variable rate jumps coupled to an ODE
#' or SDE (`vr_aggregator`), as in JumpProcesses.jl. `VR_FRM()` (the default)
#' adds one state per jump holding its integrated rate and fires the jump
#' when it reaches an exponential threshold. `VR_Direct()` (and its function
#' wrapper variant `VR_DirectFW()`, the same algorithm in R) integrates the
#' total rate of all variable rate jumps along the solution, with adaptive
#' Gauss-Kronrod quadrature on each step's interpolant, fires at the time
#' its integral reaches one exponential threshold and picks the jump in
#' proportion to the rates there.
#'
#' @return A variable rate aggregator.
#' @examples
#' f <- function(u, p, t) 0 * u
#' vrj <- VariableRateJump(function(u, p, t) 1 + sin(t), function(integ) integ$u[1] <- integ$u[1] + 1)
#' jp <- JumpProblem(ODEProblem(f, 0, c(0, 10)), Direct(), vrj, vr_aggregator = VR_Direct())
#' tail(solve(jp, Tsit5(), seed = 1)$u, 1)
#' @name variable_rate_aggregators
NULL

#' @rdname variable_rate_aggregators
#' @export
VR_FRM <- function() structure(list(name = "VR_FRM"), class = c("VR_FRM", "VariableRateAggregator"))

#' @rdname variable_rate_aggregators
#' @export
VR_Direct <- function() structure(list(name = "VR_Direct"), class = c("VR_Direct", "VariableRateAggregator"))

#' @rdname variable_rate_aggregators
#' @export
VR_DirectFW <- function() structure(list(name = "VR_DirectFW"), class = c("VR_DirectFW", "VR_Direct", "VariableRateAggregator"))

#' @rdname variable_rate_aggregators
#' @export
VariableRateAggregator <- function() VR_FRM()

# VR_Direct (VR_DirectEventCache of variable_rate.jl): one continuous
# callback whose condition is the remaining exponential threshold minus the
# integral of the total rate since the last firing.
.vr_direct_callback <- function(vrjs) {
  st <- new.env(parent = emptyenv())
  rates <- lapply(vrjs, `[[`, "rate")
  total <- function(u, p, t) sum(vapply(rates, function(r) r(u, p, t), 0))
  integral <- function(step, a, b, p) {
    if (b <= a) return(0)
    .sens_quadgk(function(s) total(.step_value(step, s), p, s), a, b, 1e-12, 1e-8)
  }
  cond <- function(u, t, integ) {
    step <- integ$last_step
    if (is.null(step)) return(st$base)
    if (!identical(st$step_t0, step$t0)) {
      # a new step: what the previous one integrated since the base time
      if (!is.null(st$step) && step$t0 > st$base_t)
        st$base <- st$base - integral(st$step, st$base_t, step$t0, integ$p)
      st$base_t <- max(st$base_t, step$t0)
      st$step_t0 <- step$t0
    }
    st$step <- step
    st$base - integral(step, st$base_t, t, integ$p)
  }
  affect <- function(integ) {
    r <- vapply(rates, function(f) f(integ$u, integ$p, integ$t), 0)
    j <- if (length(r) == 1L) 1L else max(.linear_search(r, stats::runif(1L) * sum(r)), 1L)
    vrjs[[j]]$affect(integ)
    st$base <- stats::rexp(1L)
    st$base_t <- integ$t
    st$step <- NULL
    st$step_t0 <- NULL
  }
  init <- function(integ) {
    st$base <- stats::rexp(1L)
    st$base_t <- integ$t
    st$step <- NULL
    st$step_t0 <- NULL
  }
  sp <- Reduce(`|`, lapply(vrjs, function(j) as.logical(j$save_positions %||% c(FALSE, TRUE))))
  # the condition decreases through zero: a downcrossing (affect_neg)
  cb <- ContinuousCallback(cond, affect, affect_neg = affect, initialize = init, rootfind = "left",
                           save_positions = sp)
  cb$force_save <- TRUE
  cb
}

# The aggregator, its dependency maps and the jump problem object.
.jmp_finish <- function(prob, aggregator, jumps, jset, deps, save_positions, sys, kwargs,
                        use_vrj_bounds, bracket_data = NULL) {
  if (!inherits(aggregator, "JumpAggregator")) stop("`aggregator` must be a jump aggregator", call. = FALSE)
  if (aggregator$name == "NullAggregator")
    aggregator <- .select_aggregator(jset$nmaj, jset$ncrj, length(jset$vrj), jset$nbvrj, deps)
  jset <- .jmp_partition(jset, aggregator, use_vrj_bounds)
  if (inherits(prob, "DiscreteProblem") && length(jset$vrj))
    stop("General `VariableRateJump`s require a continuous problem, like an ODE/SDE/DDE/DAE problem. To use a `DiscreteProblem` bounded `VariableRateJump`s must be used. See the JumpProcesses docs.",
         call. = FALSE)
  if (!is.null(deps$dep_graph) && aggregator$name == "Coevolve" &&
      length(deps$dep_graph) != jset$nmaj + jset$ncrj + jset$nbvrj)
    stop("Number of nodes in the dependency graph must be the same as the number of jumps.", call. = FALSE)
  structure(list(prob = prob, aggregator = aggregator, jumps = jumps, jumpset = jset,
                 dependencies = .agg_dependencies(aggregator, deps, jset), save_positions = save_positions,
                 bracket_data = bracket_data, sys = sys, kwargs = kwargs),
            class = c("JumpProblem", "SciMLProblem"))
}

.jmp_problem_sys <- function(sys, op = list(), tspan = NULL, ..., aggregator = NullAggregator(),
                             save_positions = c(TRUE, TRUE), check_initialization = FALSE) {
  .need_complete(sys, "JumpProblem")
  js <- jumps(sys)
  if (!length(js)) stop("A system without jumps cannot be used to construct a JumpProblem", call. = FALSE)
  if (isTRUE(.f(sys, "is_dde"))) stop("JumpProblem does not support systems with delays", call. = FALSE)
  if (length(constraints(sys))) stop("JumpProblem does not support constraints", call. = FALSE)
  if (length(.f(sys, "poissonians"))) stop("call mtkcompile() to turn poissonians into jumps", call. = FALSE)
  if (length(.f(sys, "brownians"))) stop("call mtkcompile() to turn brownians into noise equations", call. = FALSE)
  dots <- list(...)
  if (!is.null(dots$tstops))
    stop("Passing `tstops` directly to `JumpProblem(::System, ...)` is not supported. Define tstops on the `System` via the `tstops` keyword instead.",
         call. = FALSE)
  bracket_data <- dots$bracket_data
  dots$bracket_data <- NULL
  tspan <- .as_tspan(tspan, sys)
  has_vrjs <- any(vapply(js, inherits, TRUE, "VariableRateJump"))
  eqs <- Filter(function(e) inherits(e, "Equation"), equations(sys))
  has_noise <- !is.null(.f(sys, "noise_eqs"))
  vals <- .process_op(sys, op)
  if (anyNA(vals$u0)) .stopf("no initial value for unknown(s) %s", paste(names(vals$u0)[is.na(vals$u0)], collapse = ", "))
  # the system's events, compiled with reset_jumps = true, then the user's callbacks
  h <- .mtk$hooks$problem_callbacks
  syscb <- if (!is.null(h)) .jmp_reset_wrap(h(sys))
  kwargs <- dots
  if (!is.null(syscb)) kwargs$callback <- if (is.null(dots$callback)) syscb else CallbackSet(syscb, dots$callback)
  if (length(.f(sys, "tstops"))) kwargs$sys_tstops <- .f(sys, "tstops")
  drift <- .jmp_drift(sys, eqs)
  if (has_noise || length(eqs) || has_vrjs) {
    fn <- list(f = drift, mass_matrix = NULL, sys = sys)
    if (has_noise) {
      fn <- c(fn, .jmp_noise(sys, eqs))
      class(fn) <- c("SDEFunction", "SciMLFunction")
      inner <- .new_problem("SDEProblem", f = fn, u0 = vals$u0, tspan = tspan, p = vals$p, sys = sys,
                            noise_type = fn$noise_type, kwargs = kwargs)
    } else {
      class(fn) <- c("ODEFunction", "SciMLFunction")
      inner <- .new_problem("ODEProblem", f = fn, u0 = vals$u0, tspan = tspan, p = vals$p, sys = sys,
                            kwargs = kwargs)
    }
  } else {
    fn <- structure(list(f = function(u, p, t) u, sys = sys), class = c("DiscreteFunction", "SciMLFunction"))
    inner <- .new_problem("DiscreteProblem", f = fn, u0 = vals$u0, tspan = tspan, p = vals$p, sys = sys,
                          kwargs = kwargs)
  }
  inner$op <- vals$op
  jset <- .jmp_numeric(sys, js)
  # dependency maps of the mass action and constant rate jumps, as
  # ModelingToolkit builds them (asgraph, variable_dependencies,
  # eqeq_dependencies)
  deps <- list()
  if (length(jset$jumps) && (.needs_depgraph(aggregator) || .needs_vartojumps_map(aggregator) ||
                             aggregator$name == "NullAggregator")) {
    jdeps <- asgraph(sys, eqs = jset$jumps)
    vdeps <- variable_dependencies(sys, eqs = jset$jumps)
    deps$vartojumps_map <- jdeps$badjlist
    deps$jumptovars_map <- vdeps$badjlist
    if (.needs_depgraph(aggregator)) deps$dep_graph <- eqeq_dependencies(jdeps, vdeps)$fadjlist
  }
  .jmp_finish(inner, aggregator, js, jset, deps, save_positions, sys, kwargs, TRUE, bracket_data)
}

# Drift of a jump system, whose equations were not structurally simplified:
# the right-hand side of D(x) ~ f for each unknown x, zero for unknowns that
# change only through jumps.
.jmp_drift <- function(sys, eqs) {
  us <- unknowns(sys)
  rhs <- rep(list(.S0()), length(us))
  ukeys <- .keys(us)
  for (e in eqs) {
    if (!.is_derivative(e$lhs))
      .stopf("jump systems support only differential equations D(x) ~ f; got %s", .fmt(e))
    di <- .deriv_info(e$lhs)
    i <- match(.key(di$var), ukeys)
    if (di$order != 1L || is.na(i)) .stopf("cannot use equation %s in a jump system", .fmt(e))
    rhs[[i]] <- e$rhs
  }
  .build_closure(rhs, .codegen_args(sys), observed = .codegen_observed(sys), out = "vector")
}

# Noise of a jump system in unknown order: the noise rows belong to the
# equations, and unknowns without an equation get zero rows.
.jmp_noise <- function(sys, eqs) {
  ne <- .sde_noise_eqs(sys)
  us <- unknowns(sys)
  n <- length(us)
  row <- vapply(eqs, function(e) match(.key(.deriv_info(e$lhs)$var), .keys(us)), 0L)
  if (is.matrix(ne)) {
    G <- matrix(list(.S0()), n, ncol(ne))
    for (i in seq_along(row)) for (j in seq_len(ncol(ne))) G[[row[i], j]] <- ne[[i, j]]
    type <- .sde_noise_type(ne)
    nb <- ncol(ne)
    if (type == "scalar") G <- lapply(seq_len(n), function(i) G[[i, 1L]])
  } else {
    G <- rep(list(.S0()), n)
    for (i in seq_along(row)) G[[row[i]]] <- ne[[i]]
    type <- "diagonal"
    nb <- n
  }
  g <- .build_closure(G, .codegen_args(sys), observed = .codegen_observed(sys),
                      out = if (is.matrix(G)) "matrix" else "vector")
  list(g = g, noise_type = type, nbrownians = nb)
}

.jmp_problem_raw <- function(prob, ..., aggregator = NullAggregator(), save_positions = NULL,
                             dependencies = NULL, use_vrj_bounds = TRUE, hopping_constants = NULL,
                             spatial_system = NULL) {
  js <- list(...)
  if (length(js) && inherits(js[[1L]], "JumpAggregator")) {
    aggregator <- js[[1L]]
    js <- js[-1L]
  }
  if (length(js) == 1L && is.list(js[[1L]]) && !inherits(js[[1L]], c("AbstractJump", "MassActionJumpSet")))
    js <- js[[1L]]
  if (.is_spatial_agg(aggregator)) {
    js <- unlist(lapply(js, function(j) if (inherits(j, "MassActionJumpSet")) unclass(j) else list(j)),
                 recursive = FALSE)
    sp <- .sp_check_problem(prob, js, hopping_constants, spatial_system)
    return(structure(list(prob = prob, aggregator = aggregator, jumps = js,
                          jumpset = list(n = length(prob$u0), nmaj = 0L, ncrj = 0L, nbvrj = 0L, vrj = list()),
                          spatial = sp, save_positions = save_positions %||% c(FALSE, TRUE),
                          sys = NULL, kwargs = prob$kwargs %||% list()),
                     class = c("JumpProblem", "SciMLProblem")))
  }
  if (!is.null(hopping_constants) || !is.null(spatial_system))
    stop("`hopping_constants` and `spatial_system` need a spatial aggregator, NSM() or DirectCRDirect()",
         call. = FALSE)
  # a MassActionJump with several reactions, as in JumpProcesses
  js <- unlist(lapply(js, function(j) if (inherits(j, "MassActionJumpSet")) unclass(j) else list(j)),
               recursive = FALSE)
  save_positions <- save_positions %||% (if (inherits(prob, "DiscreteProblem")) c(FALSE, TRUE) else c(TRUE, TRUE))
  rjs <- Filter(function(j) inherits(j, "RegularJump"), js)
  if (length(rjs) > 1L) stop("a JumpProblem takes at most one RegularJump", call. = FALSE)
  out <- .jmp_finish(prob, aggregator, js, .jmp_numeric_raw(js, length(prob$u0)), dependencies %||% list(),
                     save_positions, prob$sys, prob$kwargs %||% list(), use_vrj_bounds)
  out$regular_jump <- if (length(rjs)) rjs[[1L]]
  out
}

#' @export
print.JumpProblem <- function(x, ...) {
  cat("JumpProblem with ", length(x$jumps), " jumps and aggregator ", x$aggregator$name,
      "; inner ", class(x$prob)[1L], "\n", sep = "")
  invisible(x)
}

#' @export
`$.JumpProblem` <- function(x, name) {
  if (name == "ps") return(.param_accessor(.subset2(x, "prob")))
  .subset2(x, name)
}

# ---------------------------------------------------------------------------
# Solving

#' @rdname JumpProblem
#' @param a A `JumpProblem`.
#' @param b A solver: [SSAStepper()] for pure jump problems, an ODE or SDE
#'   algorithm otherwise.
#' @export
solve.JumpProblem <- function(a, b = NULL, ...) {
  dots <- list(...)
  opts <- utils::modifyList(a$kwargs %||% list(), dots)
  # the problem's callbacks and those passed to solve are merged (merge_problem_kwargs)
  opts$callback <- CallbackSet(a$kwargs$callback, dots$callback)
  seed <- opts$seed
  opts$seed <- NULL
  inner <- a$prob
  .with_seed(seed, {
    if (inherits(b, "TauLeapingAlgorithm")) return(.leap_solve(a, b, opts))
    if (inherits(a$aggregator, "PureLeaping"))
      stop("a PureLeaping JumpProblem is solved with a tau leaping method such as SimpleTauLeaping()", call. = FALSE)
    if (inherits(inner, "DiscreteProblem")) {
      if (!is.null(b) && !inherits(b, "SSAStepper"))
        stop("pure jump problems are solved with SSAStepper()", call. = FALSE)
      .ssa_solve(a, opts)
    } else {
      if (inherits(b, "SSAStepper")) stop("SSAStepper only supports DiscreteProblems.", call. = FALSE)
      .jmp_coupled_solve(a, b, opts)
    }
  })
}

.jmp_new_agg <- function(jp) {
  if (!is.null(jp$spatial)) return(.sp_agg_new(jp))
  js <- jp$jumpset
  if (js$nmaj + js$ncrj + js$nbvrj == 0L) return(NULL)
  .agg_new(jp$aggregator, js, jp$dependencies, jp$save_positions, jp$bracket_data)
}

# SSAStepper (SSA_stepper.jl): __init, step!, solve! and savevalues!. The
# saved values follow JumpProcesses; the solution also keeps the whole
# piecewise constant path, so sol(t) is exact whatever is saved.
.ssa_solve <- function(jp, opts) {
  prob <- jp$prob
  cbs <- .as_callback_set(opts$callback)
  if (length(cbs$continuous))
    .stopf("SSAStepper does not support continuous callbacks (found %d continuous callback%s). Only DiscreteCallbacks are supported for event detection with SSAStepper. Please use an ODE/SDE solver (e.g., Tsit5()) if you need continuous event detection.",
           length(cbs$continuous), if (length(cbs$continuous) > 1L) "s" else "")
  t0 <- prob$tspan[1L]; tf <- prob$tspan[2L]
  if (tf <= t0)
    stop("The time interval to solve over is non-increasing, i.e. tspan[2] <= tspan[1]. This is not allowed for pure jump problem.",
         call. = FALSE)
  save_everystep <- any(jp$save_positions)
  save_start <- opts$save_start %||% TRUE
  save_end <- opts$save_end %||% TRUE
  saveat <- as.numeric(opts$saveat %||% numeric())
  if (length(saveat) == 1L) saveat <- seq(t0, tf, by = saveat)
  cur_saveat <- if (length(saveat) && saveat[1L] == t0) 2L else 1L
  tstops <- as.numeric(c(unlist(opts$tstops), unlist(opts$sys_tstops)))
  for (cb in cbs$discrete)
    tstops <- c(tstops, if (!is.null(cb$tstops_fn)) cb$tstops_fn(prob$tspan) else cb$tstops)
  tstops <- sort(unique(tstops[tstops > t0]))
  tidx <- 1L

  integ <- new.env(parent = emptyenv())
  integ$u <- as.numeric(prob$u0)
  names(integ$u) <- names(prob$u0)
  integ$p <- prob$p
  integ$t <- t0
  integ$tprev <- t0
  integ$uprev <- integ$u
  integ$prob <- prob
  integ$terminate <- FALSE
  integ$u_modified <- FALSE
  integ$is_ssa <- TRUE
  integ$tstop <- t0
  integ$sol_events <- list()

  ts <- numeric(); us <- list()
  save_state <- function(tt) {
    ts[length(ts) + 1L] <<- tt
    us[[length(us) + 1L]] <<- integ$u
  }
  # savevalues!: the current state, when every step is saved or forced
  savevalues <- function(force = FALSE) {
    if (save_everystep || force) {
      save_state(integ$t)
      return(TRUE)
    }
    FALSE
  }
  # apply_discrete_callback!: returns whether a callback saved
  apply_discrete <- function() {
    saved <- FALSE
    for (cb in cbs$discrete) {
      if (!isTRUE(cb$condition(integ$u, integ$t, integ))) next
      sp <- cb$save_positions %||% c(TRUE, TRUE)
      exact <- savevalues()
      saved <- TRUE
      if (sp[1L] && !exact) savevalues(TRUE)
      integ$u_modified <- TRUE
      cb$affect(integ)
      if (sp[2L]) savevalues(TRUE)
      if (isTRUE(integ$terminate)) break
    }
    saved
  }
  if (save_start) save_state(t0)
  agg <- .jmp_new_agg(jp)
  integ$jump_agg <- agg
  if (!is.null(agg)) {
    .agg_initialize(agg, integ)
    .agg_register(agg, integ)
  } else integ$tstop <- Inf
  .init_callbacks(cbs, integ)
  integ$u_modified <- FALSE
  path_t <- t0; path_u <- list(integ$u)
  # parameters in force from each time on, for observed expressions
  param_t <- t0; param_p <- list(integ$p)
  record_path <- function() {
    if (!identical(integ$u, path_u[[length(path_u)]])) {
      path_t[length(path_t) + 1L] <<- integ$t
      path_u[[length(path_u) + 1L]] <<- integ$u
    }
    if (!identical(integ$p, param_p[[length(param_p)]])) {
      param_t[length(param_t) + 1L] <<- integ$t
      param_p[[length(param_p) + 1L]] <<- integ$p
    }
  }
  njumps <- 0L
  repeat {
    if (length(integ$tstops_added)) {
      add <- integ$tstops_added
      integ$tstops_added <- NULL
      tstops <- sort(unique(c(tstops[seq_along(tstops) >= tidx], add[add > integ$t])))
      tidx <- 1L
    }
    has_tstop <- tidx <= length(tstops) && tstops[tidx] < tf
    has_jump <- integ$t < integ$tstop && integ$tstop < tf
    if (integ$terminate || !(has_jump || has_tstop)) break
    # step!
    integ$tprev <- integ$t
    integ$uprev <- integ$u
    next_jump_time <- if (integ$tstop > integ$t) integ$tstop else Inf
    doaffect <- FALSE
    if (tidx <= length(tstops) && tstops[tidx] < next_jump_time) {
      integ$t <- tstops[tidx]
      tidx <- tidx + 1L
    } else {
      integ$t <- integ$tstop
      doaffect <- TRUE
    }
    while (cur_saveat <= length(saveat) && saveat[cur_saveat] < integ$t) {
      save_state(saveat[cur_saveat])
      cur_saveat <- cur_saveat + 1L
    }
    jump_modified <- TRUE
    first <- TRUE
    while (integ$t == integ$tstop && !integ$terminate) {
      integ$u_modified <- TRUE
      if (doaffect) {
        if (.ssa_affect(agg, integ)) njumps <- njumps + 1L
      }
      jump_modified <- (!first && jump_modified) || integ$u_modified
      first <- FALSE
    }
    saved_in_cb <- if (length(cbs$discrete)) apply_discrete() else FALSE
    if (!saved_in_cb && jump_modified) savevalues()
    record_path()
  }
  if (!integ$terminate) {
    # the advance to the end time is a final step (store_uprev!)
    if (integ$t < tf) {
      integ$tprev <- integ$t
      integ$uprev <- integ$u
    }
    integ$t <- tf
    while (cur_saveat <= length(saveat) && saveat[cur_saveat] < integ$t) {
      save_state(saveat[cur_saveat])
      cur_saveat <- cur_saveat + 1L
    }
    if (length(cbs$discrete)) apply_discrete()
    record_path()
  }
  if (save_end && (!length(ts) || ts[length(ts)] != integ$t)) save_state(integ$t)
  for (cb in cbs$discrete) if (!is.null(cb$finalize)) cb$finalize(integ)
  U <- if (length(us)) do.call(rbind, us) else matrix(numeric(), 0L, length(integ$u))
  colnames(U) <- names(prob$u0)
  # piecewise constant path, right continuous at jump times; a final zero
  # length piece keeps the last time inside a half open interval
  tend <- integ$t
  const <- function(theta, dt, y0, y1, k) y0
  steps <- lapply(seq_along(path_t), function(i) {
    te <- if (i < length(path_t)) path_t[i + 1L] else tend
    list(t0 = path_t[i], dt = max(te - path_t[i], 0), tend = te, y0 = path_u[[i]],
         y1 = path_u[[i]], k = NULL, interp = const)
  })
  steps[[length(steps) + 1L]] <- list(t0 = tend, dt = 0, tend = tend, y0 = integ$u, y1 = integ$u,
                                      k = NULL, interp = const)
  sol <- .new_ode_solution(t = ts, u = U, prob = prob, alg = SSAStepper(),
                           retcode = if (integ$terminate) "Terminated" else "Success",
                           steps = steps, stats = list(njumps = njumps), p = integ$p, dense = TRUE)
  if (length(param_p) > 1L) {
    env <- .sol_env(sol)
    env$param_history <- list(t = param_t, p = do.call(rbind, param_p))
  }
  sol
}

# The aggregator's affect in SSAStepper: one jump (execute_tied_jumps! for
# SSAIntegrator); Coevolve first tests its candidate. Returns whether a jump
# was executed.
.ssa_affect <- function(agg, integ) {
  if (agg$type == "Coevolve" && !.coev_accept(agg, integ)) return(FALSE)
  .agg_execute(agg, integ)
  .agg_generate(agg, integ)
  .agg_register(agg, integ)
  TRUE
}

# Jumps coupled to an ODE or SDE: the aggregator as a discrete callback that
# stops the integrator at each jump time, and VR_FRM for the variable rate
# jumps, one state per jump holding its integrated rate.
.jmp_coupled_solve <- function(jp, alg, opts) {
  prob <- jp$prob
  js <- jp$jumpset
  n <- length(prob$u0)
  keep <- seq_len(n)
  vrjs <- js$vrj
  J <- length(vrjs)
  cbs <- list()
  agg <- .jmp_new_agg(jp)
  if (!is.null(agg)) {
    coev <- agg$type == "Coevolve"
    acb <- DiscreteCallback(
      function(u, t, integ) {
        if (!coev) return(agg$next_jump_time <= t)
        while (agg$next_jump_time <= t && t < agg$end_time) if (.coev_accept(agg, integ)) return(TRUE)
        FALSE
      },
      function(integ) .agg_execute_tied(agg, integ),
      initialize = function(integ) {
        integ$jump_agg <- agg
        .agg_initialize(agg, integ)
        .agg_register(agg, integ)
      },
      save_positions = jp$save_positions)
    acb$force_save <- TRUE
    cbs[[1L]] <- acb
  }
  direct <- inherits(jp$vr_aggregator, "VR_Direct")
  if (direct && J) {
    cbs[[length(cbs) + 1L]] <- .vr_direct_callback(vrjs)
    vrjs <- list()
    J <- 0L
  }
  for (k in seq_len(J)) cbs[[length(cbs) + 1L]] <- local({
    idx <- n + k
    j <- vrjs[[k]]
    cb <- ContinuousCallback(function(u, t, integ) u[idx],
                             function(integ) {
                               j$affect(integ)
                               integ$u[idx] <- -stats::rexp(1L)
                             }, affect_neg = NULL,
                             initialize = function(integ) integ$u[idx] <- -stats::rexp(1L),
                             rootfind = if (isTRUE(j$rootfind)) "left" else "none",
                             save_positions = j$save_positions)
    cb$force_save <- TRUE
    cb
  })
  opts$callback <- do.call(CallbackSet, c(cbs, list(opts$callback)))
  sde <- inherits(prob, "SDEProblem")
  # isadaptive(::JumpProblem, ::ImplicitEM) is false in StochasticDiffEq
  if (inherits(alg, "ImplicitEM")) opts$adaptive <- FALSE
  if (J == 0L) {
    if (sde) return(.as_rode_solution(.sde_integrate(prob, alg %||% .sde_default_alg(prob), opts)))
    alg <- alg %||% Tsit5()
    if (!inherits(alg, "ODEAlgorithm")) stop("jumps coupled to an ODE need an ODE algorithm", call. = FALSE)
    return(.integrate(prob, alg, opts))
  }
  f <- prob$f$f
  rates <- lapply(vrjs, `[[`, "rate")
  faug <- function(u, p, t) {
    x <- u[keep]
    c(f(x, p, t), vapply(rates, function(r) r(x, p, t), 0))
  }
  u0 <- c(prob$u0, stats::setNames(numeric(J), paste0("jump_u", seq_len(J))))
  aug <- prob
  aug$u0 <- u0
  aug$sys <- NULL
  M <- prob$f$mass_matrix
  Maug <- if (is.matrix(M)) {
    out <- diag(n + J)
    out[keep, keep] <- M
    out
  }
  if (sde) {
    g <- prob$f$g
    type <- prob$noise_type %||% prob$f$noise_type %||% "diagonal"
    gaug <- switch(type,
      diagonal = function(u, p, t) c(g(u[keep], p, t), numeric(J)),
      scalar = function(u, p, t) c(g(u[keep], p, t), numeric(J)),
      general = function(u, p, t) rbind(as.matrix(g(u[keep], p, t)), matrix(0, J, prob$f$nbrownians)))
    m <- if (type == "general") prob$f$nbrownians else NULL
    alg <- alg %||% .sde_default_alg(prob)
    aug$f$mass_matrix <- Maug
    sol <- .sde_integrate(aug, alg, opts, f = faug, g = gaug, u0 = u0, type = type, m = m)
  } else {
    alg <- alg %||% Tsit5()
    if (!inherits(alg, "ODEAlgorithm")) stop("jumps coupled to an ODE need an ODE algorithm", call. = FALSE)
    aug$f <- list(f = faug, mass_matrix = Maug)
    sol <- .integrate(aug, alg, opts)
  }
  out <- .jmp_strip(sol, keep, prob)
  if (sde) out <- .as_rode_solution(out)
  out
}

# Drop the integrated intensities from a solution.
.jmp_strip <- function(sol, keep, prob) {
  env <- .sol_env(sol)
  nfull <- ncol(env$U)
  steps <- env$steps
  if (!is.null(steps)) {
    # the dense output stays that of the full state, projected on the
    # unknowns (interpolation caches such as Verner's are not plain vectors)
    steps <- lapply(steps, function(s) {
      f <- s$interp; y0 <- s$y0; y1 <- s$y1
      s$interp <- function(theta, dt, a, b, k) f(theta, dt, y0, y1, k)[keep]
      s$y0 <- y0[keep]; s$y1 <- y1[keep]
      s
    })
  }
  out <- .new_ode_solution(t = env$t, u = env$U[, keep, drop = FALSE], prob = prob, alg = env$alg,
                           retcode = env$retcode, steps = steps, stats = env$stats, p = env$p_final,
                           events = env$events, dense = env$dense)
  oenv <- .sol_env(out)
  oenv$param_history <- env$param_history
  W <- env$W
  # diagonal noise draws one (unused) Wiener process per integrated intensity
  if (!is.null(W) && ncol(W$W) == nfull) W$W <- W$W[, keep, drop = FALSE]
  oenv$W <- W
  out
}
