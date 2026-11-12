# IfLifting, after ModelingToolkit.jl v11.47 (MIT) src/systems/if_lifting.jl
# (CondRewriter, new_cond_sym, the CondRewriter call, rewrite_ifs,
# discontinuities_to_ifelse, VarsUsedInCondition,
# expression_is_time_dependent, generate_condition, generate_affects,
# CONDITION_SIMPLIFIER, IfLifting) and the discontinuity table of
# Symbolics.jl src/discontinuities.jl.
#
# Conditions of if-then-else expressions (symcompR Piecewise, built by
# IfElse()) that depend on time become parameters toggled by continuous
# events, so the integrator steps exactly onto every switch. symcompR turns
# `!(f < 0)` into `f >= 0`, so comparisons are normalized here: a < b is
# `a - b < 0`, a > b is `b - a < 0`, a <= b is not `b - a < 0` and a >= b is
# not `a - b < 0`, which is what Symbolics' discontinuity table produces. A
# condition parameter holds 1 (true) or 0 (false); in conditions it appears
# as `Ne(c, 0)` and its negation as `Eq(c, 0)`, since generated code tests
# conditions with isTRUE(). symcompR simplifies the boolean combinations
# itself; .ifl_ifelse() adds the rules of CONDITION_SIMPLIFIER for
# if-then-else with boolean branches.

.IFL_RELATIONS <- c("StrictLessThan", "StrictGreaterThan", "LessThan", "GreaterThan")

#' Prevent IfLifting of an expression
#'
#' `no_if_lift(e)` is `e`; [IfLifting()] leaves the conditional expressions
#' inside it untouched.
#'
#' @param s An expression.
#' @return A `sym` expression.
#' @examples
#' t <- t_nounits
#' x <- Variable("x", t)
#' no_if_lift(IfElse(x > 1, 1, 0))
#' @export
no_if_lift <- function(s) .NO_IF_LIFT(as_sym(s))

.NO_IF_LIFT <- register_symbolic("no_if_lift", function(s) s, derivatives = list(function(s) 1))

# expression_is_time_dependent: the independent variable, a derivative or a
# variable of time occurs in e.
.ifl_time_dependent <- function(e, iv) {
  ivk <- .key(iv)
  for (v in .collect_vars(e)) {
    if (.key(v) == ivk || .type(v) %in% c("Derivative", "Function")) return(TRUE)
  }
  .key(iv) %in% .keys(free_symbols(e))
}

# ifelse(c, a, b) with the if-then-else rules of CONDITION_SIMPLIFIER.
.ifl_ifelse <- function(c, a, b) {
  if (.type(c) == "BooleanAtom") return(if (identical(.name(c), "true")) a else b)
  if (identical(.key(a), .key(b))) return(a)
  ba <- .type(a) == "BooleanAtom"; bb <- .type(b) == "BooleanAtom"
  if (ba && bb) return(if (identical(.name(a), "true")) c else Not(c))
  IfElse(c, a, b)
}

# The pairs of a Piecewise as (condition, value, rest), rest being the
# Piecewise of the remaining pairs.
.ifl_split_piecewise <- function(e) {
  pr <- .args(e)
  first <- .args(pr[[1L]])
  rest <- if (length(pr) == 1L) as_sym(NaN)
          else if (length(pr) == 2L && .type(.args(pr[[2L]])[[2L]]) == "BooleanAtom" &&
                   identical(.name(.args(pr[[2L]])[[2L]]), "true")) .args(pr[[2L]])[[1L]]
          else do.call(Piecewise, lapply(pr[-1L], function(p) list(.args(p)[[1L]], .args(p)[[2L]])))
  list(cond = first[[2L]], value = first[[1L]], rest = rest)
}

# new_cond_sym: the condition parameter of `f < 0`, created on first use;
# a repeated condition widens its dependency.
.ifl_new_cond <- function(cw, f, dep) {
  k <- .key(f)
  hit <- cw$conds[[k]]
  if (!is.null(hit)) {
    hit$dep <- Or(dep, hit$dep)
    cw$conds[[k]] <- hit
    return(hit$var)
  }
  .mtk$ifl_counter <- (.mtk$ifl_counter %||% 0L) + 1L
  v <- Parameter(paste0("cond#", .mtk$ifl_counter), default = 1,
                 description = paste("IfLifting condition", format(f), "< 0"))
  cw$conds[[k]] <- list(var = v, dep = dep, f = f)
  v
}

# The CondRewriter call: (rewritten condition, when it is true, when it is
# false), with `dep` the condition under which it is evaluated.
.ifl_cond <- function(cw, e, dep) {
  e <- as_sym(e)
  ty <- .type(e)
  same <- function() list(e, e, Not(e))
  if (ty == "BooleanAtom") return(same())
  if (.is_num(e)) return(list(e, as_sym(TRUE), as_sym(TRUE)))
  if (ty %in% c("Symbol", "Function") && .is_declared(e)) return(same())
  if (ty %in% c("Or", "And")) {
    a <- .args(e)
    acc <- .ifl_cond(cw, a[[1L]], dep)
    for (b in a[-1L]) {
      if (ty == "Or") {
        rb <- .ifl_cond(cw, b, And(dep, acc[[3L]]))
        acc <- list(Or(acc[[1L]], rb[[1L]]), Or(acc[[2L]], rb[[2L]]), And(acc[[3L]], rb[[3L]]))
      } else {
        rb <- .ifl_cond(cw, b, And(dep, acc[[2L]]))
        acc <- list(And(acc[[1L]], rb[[1L]]), And(acc[[2L]], rb[[2L]]), Or(acc[[3L]], rb[[3L]]))
      }
    }
    return(acc)
  }
  if (ty == "Piecewise") {
    sp <- .ifl_split_piecewise(e)
    rc <- .ifl_cond(cw, sp$cond, dep)
    ra <- .ifl_cond(cw, sp$value, And(dep, rc[[2L]]))
    rb <- .ifl_cond(cw, sp$rest, And(dep, rc[[3L]]))
    return(list(.ifl_ifelse(rc[[1L]], ra[[1L]], rb[[1L]]),
                Or(And(rc[[2L]], ra[[2L]]), And(rc[[3L]], rb[[2L]])),
                Or(And(rc[[2L]], ra[[3L]]), And(rc[[3L]], rb[[3L]]))))
  }
  if (ty == "Not") {
    r <- .ifl_cond(cw, .args(e)[[1L]], dep)
    return(list(Not(r[[1L]]), r[[3L]], r[[2L]]))
  }
  if (ty %in% .IFL_RELATIONS) {
    if (!.ifl_time_dependent(e, cw$iv)) return(same())
    a <- .args(e)
    f <- switch(ty, StrictLessThan = , GreaterThan = a[[1L]] - a[[2L]], a[[2L]] - a[[1L]])
    cv <- .ifl_new_cond(cw, f, dep)
    yes <- Ne(cv, 0); no <- Eq(cv, 0)
    return(if (ty %in% c("LessThan", "GreaterThan")) list(no, no, yes) else list(yes, yes, no))
  }
  if (ty %in% c("Equality", "Unequality")) return(same())
  if (!.ifl_time_dependent(e, cw$iv)) return(same())
  .stopf(paste("Unsupported expression form in decision variable computation %s. If this is not",
               "meant to be transformed via `IfLifting`, wrap the parent expression in `no_if_lift`."),
         format(e))
}

# rewrite_ifs: replace the conditions of every if-then-else in e. A
# comparison outside a condition is the if-then-else (cond, true, false) of
# Symbolics' discontinuity table.
.ifl_rewrite <- function(cw, e, dep) {
  if (!is_sym(e) || !length(.args(e))) return(e)
  ty <- .type(e)
  if (ty == "Function" && identical(.name(e), "no_if_lift")) return(e)
  if (ty == "Function" && identical(.name(e), "ifl_root")) {
    .ifl_cond(cw, StrictLessThan(.args(e)[[2L]], 0), dep)
    return(.ifl_rewrite(cw, .args(e)[[1L]], dep))
  }
  if (ty == "Piecewise") {
    sp <- .ifl_split_piecewise(e)
    rc <- .ifl_cond(cw, sp$cond, dep)
    return(.ifl_ifelse(rc[[1L]], .ifl_rewrite(cw, sp$value, rc[[2L]]), .ifl_rewrite(cw, sp$rest, rc[[3L]])))
  }
  if (ty %in% .IFL_RELATIONS) return(.ifl_cond(cw, e, dep)[[1L]])
  .rebuild(e, lapply(.args(e), function(a) .ifl_rewrite(cw, a, dep)))
}

# A marker for `value` whose discontinuities are the zero crossings of
# `root`. symcompR merges equal if-then-else branches, so the marker carries
# the root until .ifl_rewrite turns it into a condition parameter and drops
# it.
.IFL_ROOT <- register_symbolic("ifl_root", function(value, root) value,
                               derivatives = list(function(value, root) 1, function(value, root) 0))

# discontinuities_to_ifelse with the table of Symbolics: abs(x) is
# ifelse(x < 0, -x, x), max(x, y) is ifelse(x - y < 0, y, x), min(x, y)
# is ifelse(x - y < 0, x, y) and mod(x, y) (R's `%%`) stops where x is a
# multiple of y. Symbolics registers mod with the root function
# _triangle(x, y), which is never positive (it touches zero at x = y / 2 and
# never crosses it), so ModelingToolkit's event never fires; the root here
# is sin(pi x / y), which changes sign at every multiple of y. Julia's div
# corresponds to R's `%/%`, which symcompR writes as floor(x / y); floor and
# ceiling stop at the integers of their argument. Julia's rem (truncated
# remainder) has no symcompR counterpart.
.ifl_discontinuities <- function(e) {
  if (!is_sym(e) || !length(.args(e))) return(e)
  if (.type(e) == "Function" && identical(.name(e), "no_if_lift")) return(e)
  args <- lapply(.args(e), .ifl_discontinuities)
  nm <- if (.type(e) == "Function") .name(e) else .type(e)
  pairwise <- function(f) Reduce(f, args)
  switch(nm,
    Abs = IfElse(StrictLessThan(args[[1L]], 0), -args[[1L]], args[[1L]]),
    Max = pairwise(function(x, y) IfElse(StrictLessThan(x - y, 0), y, x)),
    Min = pairwise(function(x, y) IfElse(StrictLessThan(x - y, 0), x, y)),
    Mod = .IFL_ROOT(.rebuild(e, args), sin(pi * args[[1L]] / args[[2L]])),
    floor = , ceiling = .IFL_ROOT(.rebuild(e, args), sin(pi * args[[1L]])),
    .rebuild(e, args))
}

# VarsUsedInCondition: variables in the conditions of if-then-else
# expressions, outside no_if_lift.
.ifl_condition_vars <- function(e, acc = list()) {
  if (!is_sym(e) || !length(.args(e))) return(acc)
  if (.type(e) == "Function" && identical(.name(e), "no_if_lift")) return(acc)
  if (.type(e) == "Piecewise") {
    for (p in .args(e)) acc <- c(acc, .collect_vars(.args(p)[[2L]]))
  }
  for (a in .args(e)) acc <- .ifl_condition_vars(a, acc)
  acc
}

# observed_equations_used_by: indices of the observed equations the
# equations depend on, directly or through other observed equations.
.ifl_observed_used <- function(eqs, obs) {
  lhs <- vapply(obs, function(o) .key(o$lhs), "")
  need <- .keys(.collect_vars(unlist(lapply(eqs, function(e) list(e$lhs, e$rhs)), recursive = FALSE)))
  used <- integer()
  repeat {
    new <- setdiff(which(lhs %in% need), used)
    if (!length(new)) break
    used <- c(used, new)
    need <- union(need, .keys(.collect_vars(lapply(obs[new], `[[`, "rhs"))))
  }
  sort(used)
}

#' Turn conditional equations into events
#'
#' `IfLifting` is an additional pass of [mtkcompile()]:
#' `mtkcompile(sys, additional_passes = list(IfLifting))`. Every condition
#' `f < 0` of an if-then-else ([IfElse()]) that depends on time, and every
#' `abs`, `max`, `min`, `%%`, `%/%`, `floor` and `ceiling` (rewritten as
#' if-then-else or stopped at their jumps), becomes a parameter
#' that holds the value of the condition, and a continuous event at the zero
#' crossings of `f` updates it. The integrator then never steps across a
#' discontinuity: it stops on it and continues on the other branch. Wrap an
#' expression in [no_if_lift()] to keep it as it is. This is an experimental
#' pass, as in ModelingToolkit.
#'
#' @param sys A system returned by `mtkcompile()`.
#' @return The system with condition parameters and their events.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); y <- Variable("y", t)
#' sys <- System(list(D(x) ~ abs(y), y ~ sin(t)), t, name = "absy")
#' lifted <- mtkcompile(sys, additional_passes = list(IfLifting))
#' equations(lifted); parameters(lifted)
#' sol <- solve(ODEProblem(lifted, list(x ~ 0), c(0, 5)))
#' sol(pi + 1, idxs = x) # 3 + cos(pi + 1)
#' @export
IfLifting <- function(sys) {
  iv <- .f(sys, "iv")
  if (is.null(iv)) stop("IfLifting is only supported for time-dependent systems.", call. = FALSE)
  cw <- new.env(parent = emptyenv())
  cw$iv <- iv
  cw$conds <- list()
  eqs <- .f(sys, "eqs")
  obs <- .f(sys, "observed")
  used <- .ifl_observed_used(eqs, obs)
  # integer and boolean variables used in conditions are rewritten as
  # conditions themselves when they are observed
  cvars <- .ifl_condition_vars(NULL)
  for (e in c(lapply(eqs, `[[`, "rhs"), lapply(obs[used], `[[`, "rhs"))) cvars <- .ifl_condition_vars(e, cvars)
  cvars <- Filter(function(v) !.is_derivative(v) && !isparameter(v) && .key(v) != .key(iv) &&
                    isTRUE(.meta_get(v, "variable_type", "") %in% c("Bool", "Integer", "logical", "integer")),
                  .unique_syms(cvars))
  cvk <- .keys(cvars)
  eqs <- lapply(eqs, function(e) Equation(e$lhs, .ifl_rewrite(cw, .ifl_discontinuities(e$rhs), as_sym(TRUE))))
  for (i in used) {
    rhs <- .ifl_discontinuities(obs[[i]]$rhs)
    obs[[i]] <- Equation(obs[[i]]$lhs, if (.key(obs[[i]]$lhs) %in% cvk) .ifl_cond(cw, rhs, as_sym(TRUE))[[1L]]
                                      else .ifl_rewrite(cw, rhs, as_sym(TRUE)))
  }
  conds <- unname(cw$conds)
  vars <- lapply(conds, `[[`, "var")
  vk <- .keys(vars)
  # the condition variables each dependency involves
  deps <- lapply(conds, function(cd) which(vk %in% .keys(.collect_vars(cd$dep))))
  callbacks <- list()
  ics <- .f(sys, "initial_conditions")
  for (i in seq_along(conds)) {
    # this condition and those whose dependency involves it, transitively
    upd <- i
    repeat {
      more <- setdiff(which(vapply(deps, function(d) any(d %in% upd), TRUE)), upd)
      if (!length(more)) break
      upd <- c(upd, more)
    }
    names_upd <- vapply(vars[upd], .name, "")
    affect <- function(op) {
      obs_exprs <- lapply(conds[upd], function(cd) if (op == "<") StrictLessThan(cd$f, 0) else LessThan(cd$f, 0))
      ImperativeAffect(function(m, o, ...) o, modified = stats::setNames(vars[upd], names_upd),
                       observed = stats::setNames(obs_exprs, names_upd), skip_checks = TRUE)
    }
    up <- affect("<"); down <- affect("<=")
    dep <- conds[[i]]$dep
    zc <- if (.type(dep) == "BooleanAtom" && identical(.name(dep), "true")) conds[[i]]$f
          else IfElse(dep, conds[[i]]$f, NaN)
    callbacks[[i]] <- SymbolicContinuousCallback(list(Equation(zc, 0)), up, affect_neg = down,
                                                 initialize = up, rootfind = "right")
    ics <- smap_set(ics, vars[[i]], 1)
  }
  .update_system(sys, eqs = eqs, observed = obs, ps = .unique_syms(c(.f(sys, "ps"), vars)),
                 initial_conditions = ics,
                 continuous_events = c(.f(sys, "continuous_events"), callbacks))
}
