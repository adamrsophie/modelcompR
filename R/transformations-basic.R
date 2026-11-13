# System transformations, after
# ModelingToolkitBase/src/systems/diffeqs/basic_transformations.jl:
# liouville_transform, change_of_variables, fractional_to_ordinary,
# linear_fractional_to_ordinary, change_independent_variable,
# add_accumulations, convert_system_indepvar and respecialize. The
# stochastic transformations (stochastic_integral_transform,
# Girsanov_transform, noise_to_brownians) belong to the stochastic module.

# A variable application with new arguments, keeping the metadata.
.recall <- function(v, args) {
  out <- do.call(Function(.name(v)), as.list(args))
  m <- attr(v, "mtk_meta")
  if (is.null(m)) m <- .meta(v)
  attr(out, "mtk_meta") <- m
  out
}

.pair_sides <- function(p) {
  if (inherits(p, "formula")) p <- as_equation(p)
  if (inherits(p, "Equation")) return(list(p$lhs, p$rhs))
  if (is.list(p) && length(p) == 2L) return(list(as_sym(p[[1L]]), as_sym(p[[2L]])))
  stop("a substitution is a formula a ~ b or a pair list(a, b)", call. = FALSE)
}

.as_pairs <- function(x) {
  if (inherits(x, c("formula", "Equation"))) x <- list(x)
  if (is.list(x) && length(x) == 2L && is_sym(x[[1L]]) && !inherits(x[[1L]], "Equation")) x <- list(x)
  lapply(x, .pair_sides)
}

#' Liouville transformation
#'
#' `liouville_transform(sys)` appends the variable `trJ` with
#' `D(trJ) ~ -trJ * tr(J)`, where `J` is the Jacobian of the right-hand side.
#' Started from the density of the initial condition, `trJ` follows the
#' density along each trajectory, which propagates uncertainty in the initial
#' conditions or parameters (Halder, Lee and Bhattacharya, 2013).
#'
#' @param sys An ODE `System` whose equations are `D(x) ~ f`.
#' @param ... Passed to [System()].
#' @return A `System` with the extra unknown `trJ`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)", "y(t)"); parameters("a", "b", "c", "d")
#' sys <- System(list(D(x) ~ a * x - b * x * y, D(y) ~ -d * y + c * x * y), t, name = "lv")
#' equations(liouville_transform(sys))[[3]]
#' @export
liouville_transform <- function(sys, ...) {
  t <- get_iv(sys)
  trJ <- Variable("trJ", t)
  eqs <- Filter(function(e) inherits(e, "Equation"), equations(sys))
  us <- unknowns(sys)
  others <- c(.parameters_of(sys, TRUE), list(t))
  obs <- observed(sys)
  orule <- if (length(obs)) .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs")) else NULL
  tr <- .S0()
  for (e in eqs) {
    if (!.is_derivative(e$lhs)) next
    v <- .deriv_info(e$lhs)$var
    rhs <- if (is.null(orule)) e$rhs else .subst_fixpoint(e$rhs, orule)
    J <- .partials(list(rhs), list(v), c(Filter(function(u) .key(u) != .key(v), us), others))
    tr <- tr + J[[1L, 1L]]
  }
  neweq <- Equation(Differential(t)(trJ), trJ * -tr)
  System(c(eqs, list(neweq)), t, c(us, list(trJ)), .parameters_of(sys), checks = FALSE,
         name = nameof(sys), ...)
}

#' Change of variables
#'
#' `change_of_variables(sys, iv, forward_subs, backward_subs)` rewrites an ODE
#' (or SDE with brownian variables or noise equations) in new variables `z =
#' f(x)`: `forward_subs` gives `f(x) ~ z` and `backward_subs` the inverse `x ~
#' g(z)`. The new equations follow from the chain rule (Itô's lemma for
#' stochastic systems, using the second derivative of each transformation
#' with respect to its own variable). Initial values of the new variables are
#' computed from those of the old ones.
#'
#' @param sys A `System`.
#' @param iv The independent variable.
#' @param forward_subs,backward_subs Lists of formulas `f(x) ~ z` and
#'   `x ~ g(z)` (or pairs).
#' @param simplify Simplify the new equations and compile the result with
#'   [mtkcompile()].
#' @param t0 Time at which initial values are computed, when they depend on
#'   time.
#' @param isSDE Treat the system as stochastic.
#' @return The transformed `System`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); z <- Variable("z", t); a <- Parameter("a")
#' sys <- mtkcompile(System(list(D(x) ~ a * x), t, name = "sys",
#'                          initial_conditions = list(x ~ 1, a ~ -0.5)))
#' new_sys <- change_of_variables(sys, t, list(log(x) ~ z), list(x ~ exp(z)))
#' equations(new_sys)
#' @export
change_of_variables <- function(sys, iv, forward_subs, backward_subs, simplify = TRUE, t0 = NA,
                                isSDE = FALSE) {
  t <- as_sym(iv)
  fwd <- .as_pairs(forward_subs)
  bwd <- .as_pairs(backward_subs)
  old_vars <- lapply(bwd, `[[`, 1L)
  new_vars <- lapply(fwd, `[[`, 2L)
  old_eqs <- Filter(function(e) inherits(e, "Equation"), equations(sys))
  neqs <- .f(sys, "noise_eqs")
  brownvars <- get_brownians(sys)
  if (is.null(neqs) && !length(brownvars)) {
    neqs <- rep(list(list()), length(old_eqs))
  } else if (!is.null(neqs)) {
    isSDE <- TRUE
    if (is.matrix(neqs)) {
      m <- ncol(neqs)
      neqs <- lapply(seq_len(nrow(neqs)), function(i) lapply(seq_len(m), function(j) as_sym(neqs[[i, j]])))
    } else {
      neqs <- lapply(neqs, function(n) list(as_sym(n)))
    }
    brownvars <- lapply(seq_along(neqs[[1L]]), function(i) Brownian(paste0("B_", i)))
  } else {
    isSDE <- TRUE
    zero <- .rule(brownvars, rep(list(.S0()), length(brownvars)))
    neqs <- vector("list", length(old_eqs))
    for (i in seq_along(old_eqs)) {
      rhs <- old_eqs[[i]]$rhs
      neqs[[i]] <- lapply(brownvars, function(B) diff(rhs, B))
      old_eqs[[i]] <- Equation(old_eqs[[i]]$lhs, .subst(rhs, zero))
    }
  }
  D <- Differential(t)
  fwd_rule <- .rule(lapply(fwd, `[[`, 1L), new_vars)
  bwd_rule <- .rule(old_vars, lapply(bwd, `[[`, 2L))
  new_eqs <- list()
  for (k in seq_along(fwd)) {
    f <- fwd[[k]][[1L]]
    ex <- expand_derivatives(D(f))
    first <- diff(f, old_vars[[k]])
    second <- diff(first, old_vars[[k]])
    for (i in seq_along(old_eqs)) {
      eq <- old_eqs[[i]]
      if (!(.key(eq$lhs) %in% .keys(.collect_vars(ex)))) next
      ex <- .subst(ex, .rule(list(eq$lhs), list(eq$rhs)))
      if (isSDE) {
        for (j in seq_along(neqs[[i]])) {
          noise <- neqs[[i]][[j]]
          ex <- ex + noise^2 * second / 2 + noise * first * brownvars[[j]]
        }
      }
    }
    ex <- .subst(.subst(ex, fwd_rule), bwd_rule)
    if (simplify) ex <- symcompR::simplify(symcompR::expand(ex))
    new_eqs[[k]] <- Equation(D(new_vars[[k]]), ex)
  }
  ics <- initial_conditions(sys)
  ic_vars <- smap_vars(ics)
  ic_rule <- .rule(unname(ic_vars), lapply(smap_values(ics), as_sym))
  new_ics <- symmap()
  for (p in fwd) {
    ex <- .subst(p[[1L]], ic_rule)
    if (!is.na(t0)) ex <- .subst(ex, .rule(list(t), list(t0)))
    if (!length(.collect_vars(ex))) ex <- as.numeric(N(ex))
    new_ics <- smap_set(new_ics, p[[2L]], ex)
  }
  for (p in .parameters_of(sys)) if (smap_has(ics, p)) new_ics <- smap_set(new_ics, p, smap_get(ics, p))
  back_eqs <- lapply(bwd, function(p) Equation(p[[1L]], p[[2L]]))
  new_sys <- System(c(new_eqs, back_eqs), t, initial_conditions = new_ics,
                    observed = observed(sys), brownians = if (isSDE) brownvars else NULL,
                    name = "new_sys")
  if (simplify) return(mtkcompile(new_sys))
  new_sys
}

#' Change the independent variable of a system
#'
#' `change_independent_variable(sys, iv)` turns the dependent variable `iv`
#' (for example `u(t)`) into the independent variable of `sys`. Every
#' variable `f(t)` becomes `f(u)`, derivatives follow the chain rule
#' \ifelse{latex}{\verb{d/dt = u_t(u) d/du}}{\verb{d/dt = uˍt(u) d/du}} with the new unknown \ifelse{latex}{\verb{u_t(u)}}{\verb{uˍt(u)}} for `du/dt`, and the
#' old independent variable becomes the unknown `t(u)`. The change is valid
#' when `u` is strictly monotonic in `t`.
#'
#' Extra equations `eqs` relating the variables (such as `t ~ f(u)` for
#' non-autonomous systems) are transformed too; `add_old_diff = TRUE` adds
#' \ifelse{latex}{\verb{D(t(u)) ~ 1 / u_t(u)}}{\verb{D(t(u)) ~ 1 / uˍt(u)}} by the inverse function rule. The change must be
#' made before [mtkcompile()]. Subsystems are transformed recursively; the
#' new variables are global, so `iv` should not be namespaced (use
#' [GlobalScope()] in subsystems).
#'
#' The transformed equations are implicit in the highest derivatives (for
#' example \ifelse{latex}{\verb{u_t^2 y'' + ...}}{\verb{uˍt^2 y'' + ...}}), so compiling them needs the structural
#' pipeline of [mtkcompile()].
#'
#' @param sys A `System` that is not compiled.
#' @param iv The new independent variable, an unknown `u(t)` of `sys`.
#' @param eqs Extra equations.
#' @param add_old_diff Add a differential equation for the old independent
#'   variable.
#' @param simplify Simplify the expanded derivatives.
#' @param fold Unused; kept for compatibility.
#' @return The transformed `System`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); y <- Variable("y", t)
#' M <- System(list(D(D(y)) ~ -9.81, D(x) ~ 1), t, name = "M")
#' M2 <- change_independent_variable(M, x)
#' equations(M2)
#' @export
change_independent_variable <- function(sys, iv, eqs = list(), add_old_diff = FALSE,
                                        simplify = TRUE, fold = FALSE) {
  iv2_of_iv1 <- as_sym(iv)
  iv1 <- get_iv(sys)
  if (isTRUE(.f(sys, "is_dde")))
    .stopf("System %s contains delay differential equations (DDEs). This is currently not supported!", nameof(sys))
  if (isTRUE(.f(sys, "isscheduled")))
    .stopf("System %s is structurally simplified. Change independent variable before structural simplification!", nameof(sys))
  if (is.null(iv1) || .type(iv2_of_iv1) != "Function" || length(.args(iv2_of_iv1)) != 1L ||
      .key(.args(iv2_of_iv1)[[1L]]) != .key(iv1))
    .stopf("Variable %s is not a function of the independent variable %s of the system %s!",
           format(iv2_of_iv1), if (is.null(iv1)) "(none)" else format(iv1), nameof(sys))
  iv1name <- .name(iv1)
  iv2name <- .name(iv2_of_iv1)
  D1 <- Differential(iv1)
  u2 <- .meta_get(iv2_of_iv1, "unit")
  u1 <- .meta_get(iv1, "unit")
  iv2 <- .make_var(iv2name, "iv", NULL, if (is.null(u2)) list() else list(unit = u2))
  iv1_of_iv2 <- GlobalScope(.make_var(iv1name, "unknown", list(iv2),
                                      if (is.null(u2)) list() else list(unit = u1)))
  div2_of_iv1 <- default_toterm(D1(iv2_of_iv1))
  if (!is.null(u2)) div2_of_iv1 <- setmetadata(div2_of_iv1, "unit", paste0("(", u2, ")/(", u1 %||% "1", ")"))
  div2_of_iv1 <- GlobalScope(div2_of_iv1)
  w <- .recall(div2_of_iv1, list(iv2))
  eqs <- .as_equations(eqs)
  if (add_old_diff) eqs <- c(eqs, list(Equation(Differential(iv2)(iv1_of_iv2), 1 / w)))
  sys <- .update_system(sys, eqs = c(get_eqs(sys), eqs),
                        unknowns = c(get_unknowns(sys), list(iv1, div2_of_iv1)))

  # d/dt = w d/du applied n times to g(u)
  chain <- function(g, n) {
    for (k in seq_len(n)) g <- w * diff(g, iv2)
    if (simplify && n > 1L) g <- symcompR::expand(g)
    g
  }
  is_fun_of_iv1 <- function(v) .type(v) == "Function" && length(.args(v)) == 1L &&
    .key(.args(v)[[1L]]) == .key(iv1)
  tr <- function(ex) {
    if (!is_sym(ex)) return(ex)
    ex <- expand_derivatives(ex)
    from <- list(); to <- list()
    for (v in .collect_vars(ex)) {
      if (.type(v) == "Derivative") {
        di <- .deriv_info(v)
        if (di$n_specs != 1L || .key(di$iv) != .key(iv1) || !is_fun_of_iv1(di$var)) next
        rep <- if (.key(di$var) == .key(iv2_of_iv1)) chain(w, di$order - 1L) else
          chain(.recall(di$var, list(iv2)), di$order)
      } else if (is_fun_of_iv1(v)) {
        rep <- if (.key(v) == .key(iv2_of_iv1)) iv2 else .recall(v, list(iv2))
      } else if (.key(v) == .key(iv1)) {
        rep <- iv1_of_iv2
      } else next
      from[[length(from) + 1L]] <- v
      to[[length(to) + 1L]] <- rep
    }
    if (!length(from)) return(ex)
    .subst(ex, .rule(from, to))
  }
  tr_eq <- function(e) {
    if (inherits(e, "Equation")) return(Equation(tr(e$lhs), tr(e$rhs)))
    if (inherits(e, "Connection")) {
      e$systems <- lapply(e$systems, function(s) if (inherits(s, "System")) tr_sys(s) else tr(s))
      return(e)
    }
    e
  }
  tr_map <- function(m) {
    out <- symmap()
    for (el in unclass(m %||% symmap())) {
      if (is.null(el$var)) {
        out[[paste0("@name:", el$name)]] <- el
        next
      }
      out <- smap_set(out, tr(el$var), if (is_sym(el$value)) tr(el$value) else el$value)
    }
    class(out) <- "mtk_symmap"
    out
  }
  tr_sys <- function(s) {
    subs <- lapply(get_systems(s), tr_sys)
    us <- Filter(function(v) .key(v) != .key(iv2), lapply(get_unknowns(s), tr))
    ps <- Filter(function(p) !isinitial(p), lapply(get_ps(s), tr))
    asserts <- lapply(.assertion_pairs(assertions(s)), function(pr) list(tr(pr$cond), pr$msg))
    out <- System(lapply(get_eqs(s), tr_eq), iv2, us, ps, name = nameof(s),
                  observed = lapply(get_observed(s), tr_eq),
                  initialization_eqs = lapply(get_initialization_eqs(s), tr_eq),
                  bindings = tr_map(get_bindings(s)), initial_conditions = tr_map(get_initial_conditions(s)),
                  guesses = tr_map(get_guesses(s)), connector_type = get_connector_type(s),
                  assertions = asserts, description = get_description(s))
    if (length(subs)) out <- compose(out, subs)
    if (iscomplete(s)) out <- complete(out, flatten = !length(subs))
    out
  }
  tr_sys(sys)
}

#' Add accumulation variables
#'
#' `add_accumulations(sys)` adds, for every unknown `x`, a variable
#' `accumulation_x` with `D(accumulation_x) ~ x` and initial value 0. With a
#' list of formulas `acc ~ expr` (or pairs), it adds `D(acc) ~ expr` for each.
#'
#' @param sys A `System`.
#' @param vars Unknowns to accumulate, or a list of formulas `acc ~ expr`.
#' @return The extended `System`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t)
#' sys <- add_accumulations(System(list(D(x) ~ -x), t, name = "s"))
#' equations(sys)
#' @export
add_accumulations <- function(sys, vars = unknowns(sys)) {
  if (inherits(vars, c("formula", "Equation"))) vars <- list(vars)
  pairs <- if (length(vars) && all(vapply(vars, is_sym, TRUE))) {
    lapply(.flatten_vars(vars), function(v) {
      m <- .meta(v)
      m$kind <- NULL
      a <- .make_var(paste0("accumulation_", .name(v)), "unknown",
                     if (.type(v) == "Function") .args(v) else NULL, m)
      list(a, v)
    })
  } else .as_pairs(vars)
  avars <- lapply(pairs, `[[`, 1L)
  ints <- Filter(function(a) .contains_key(unknowns(sys), a), avars)
  if (length(ints)) stop(paste(vapply(ints, format, ""), collapse = ", "), " already exist in the system!",
                         call. = FALSE)
  D <- Differential(get_iv(sys))
  ics <- .f(sys, "initial_conditions")
  for (a in avars) ics <- smap_set(ics, a, 0)
  .update_system(sys, eqs = c(get_eqs(sys), lapply(pairs, function(p) Equation(D(p[[1L]]), p[[2L]]))),
                 unknowns = c(get_unknowns(sys), avars), initial_conditions = ics)
}

#' Substitute the independent variable
#'
#' `convert_system_indepvar(sys, t)` makes `t` the independent variable of
#' `sys` by substitution: unknowns `x(s)` become `x(t)` and time independent
#' unknowns `x` become `x(t)`. Unlike [change_independent_variable()] no
#' chain rule is applied. The system must not be compiled (no observed
#' equations); a nonlinear system becomes a time dependent one.
#'
#' @param sys A `System`.
#' @param t The new independent variable.
#' @param name The name of the result.
#' @return The converted `System`.
#' @examples
#' x <- Variable("x"); a <- Parameter("a")
#' ns <- System(list(0 ~ x^2 - a), name = "ns")
#' tau <- independent_variable("tau")
#' unknowns(convert_system_indepvar(ns, tau))
#' @export
convert_system_indepvar <- function(sys, t, name = nameof(sys)) {
  if (length(observed(sys)))
    stop("convert_system_indepvar() cannot handle reduced model (i.e. observed(sys) is non-empty).",
         call. = FALSE)
  t <- as_sym(t)
  from <- list(); to <- list()
  sts <- unknowns(sys)
  newsts <- lapply(sts, function(s) {
    if (.type(s) == "Function") {
      if (length(.args(s)) != 1L)
        stop("Illegal unknown: ", format(s), ". The unknown can have at most one argument like `x(t)`.",
             call. = FALSE)
      if (.key(.args(s)[[1L]]) == .key(t)) return(s)
      ns <- .recall(s, list(t))
    } else {
      m <- .meta(s)
      m$kind <- NULL
      ns <- .make_var(.name(s), "unknown", list(t), m)
    }
    from[[length(from) + 1L]] <<- s
    to[[length(to) + 1L]] <<- ns
    ns
  })
  old_iv <- get_iv(sys)
  if (!is.null(old_iv)) {
    from[[length(from) + 1L]] <- old_iv
    to[[length(to) + 1L]] <- t
  }
  rule <- .rule(from, to)
  sub <- function(e) if (is_sym(e) || inherits(e, "Equation")) .subst(e, rule) else e
  submap <- function(m) {
    out <- symmap()
    for (el in unclass(m)) out <- smap_set(out, sub(el$var), sub(el$value))
    out
  }
  noise <- .f(sys, "noise_eqs")
  if (!is.null(noise)) noise <- if (is.matrix(noise)) { noise[] <- lapply(noise, sub); noise } else lapply(noise, sub)
  .update_system(sys, eqs = lapply(equations(sys), sub), iv = t, unknowns = newsts,
                 bindings = submap(bindings(sys)), initial_conditions = submap(initial_conditions(sys)),
                 name = as.character(name), noise_eqs = noise,
                 constraints = lapply(get_constraints(sys), function(c) {
                   if (inherits(c, "Inequality")) { c$lhs <- sub(c$lhs); c$rhs <- sub(c$rhs); c } else sub(c)
                 }),
                 costs = lapply(get_costs(sys), sub), systems = list())
}

#' Specialize nonnumeric parameters
#'
#' In ModelingToolkit, `respecialize(sys, mapping)` changes the symbolic type
#' of nonnumeric parameters (such as callable parameters holding
#' interpolants) to the concrete type of their values. R is untyped, so the
#' type change itself does not apply: `respecialize()` checks that each
#' parameter exists and has a nonnumeric value (from `mapping`, a list of
#' parameters or formulas `p ~ value` / pairs `list(p, value)`, or from the
#' initial conditions), sets that value as the parameter's initial condition
#' and records `class(value)[1]` in the `variable_type` metadata. With `all =
#' TRUE` every parameter whose initial value is not a number is respecialized.
#' Only complete systems can be respecialized.
#'
#' @param sys A complete `System`.
#' @param mapping Parameters, or pairs of parameters and values.
#' @param all Respecialize all nonnumeric parameters.
#' @return The respecialized, complete `System`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); f <- Parameter("f")
#' sys <- complete(System(list(D(x) ~ -x), t, ps = list(f), name = "s"))
#' sys2 <- respecialize(sys, list(list(f, splinefun(1:3, c(1, 4, 9)))))
#' getmetadata(parameters(sys2)[[1]], "variable_type")
#' @export
respecialize <- function(sys, mapping = list(), all = length(mapping) == 0L) {
  if (!iscomplete(sys))
    stop("This operation can only be performed on completed systems. Use complete(sys) or mtkcompile(sys).",
         call. = FALSE)
  ics <- .f(sys, "initial_conditions")
  items <- list()
  for (m in mapping) {
    if (is_sym(m)) items[[length(items) + 1L]] <- list(p = m, value = NULL)
    else if (inherits(m, "formula")) {
      # the caller's own formula, evaluated where it was written
      env <- environment(m)
      items[[length(items) + 1L]] <- list(p = eval(m[[2L]], env), value = eval(m[[3L]], env))
    } else items[[length(items) + 1L]] <- list(p = as_sym(m[[1L]]), value = m[[2L]])
  }
  if (all) {
    given <- .keys(lapply(items, `[[`, "p"))
    for (p in get_ps(sys)) {
      v <- smap_get(ics, p)
      if (.key(p) %in% given || is.null(v) || is.numeric(v) || is_sym(v) || is.logical(v)) next
      items[[length(items) + 1L]] <- list(p = p, value = NULL)
    }
  }
  ps <- get_ps(sys)
  for (it in items) {
    k <- it$p
    v <- it$value %||% smap_get(ics, k)
    if (is.null(v)) stop("Parameter ", format(k), " needs an associated value to be respecialized.", call. = FALSE)
    if (is_sym(v)) stop("Parameter ", format(k), " needs an associated value to be respecialized. Found symbolic initial value ",
                        format(v), ".", call. = FALSE)
    if (is.numeric(v) || is.logical(v))
      stop("Numeric types cannot be respecialized; tried to respecialize ", format(k), ".", call. = FALSE)
    idx <- match(.key(k), .keys(ps))
    if (is.na(idx)) stop("Parameter ", format(k), " does not exist in the system.", call. = FALSE)
    new_p <- setmetadata(ps[[idx]], "variable_type", class(v)[1L])
    ps[[idx]] <- new_p
    ics <- smap_set(ics, new_p, v)
  }
  sys <- .update_system(sys, ps = ps, initial_conditions = ics)
  complete(sys, flatten = FALSE)
}

# ---------------------------------------------------------------------------
# Fractional differential equations

# Quadrature of the fractional kernel s^(alpha - 1) / gamma(alpha) as
# sum_i c_i exp(-gamma_i s), valid on [delta, T] to accuracy epsilon
# (fto_helper in basic_transformations.jl).
.fto_quadrature <- function(alpha, epsilon, T) {
  delta <- (gamma(alpha + 1) * epsilon)^(1 / alpha)
  a <- pi / 2 * (1 - (1 - alpha) / ((2 - alpha) * log(1 / epsilon)))
  h <- 2 * pi * a / log(1 + (2 / epsilon * cos(a)^(alpha - 1)))
  x_sub <- (gamma(2 - alpha) * epsilon)^(1 / (1 - alpha))
  x_sup <- -log(gamma(1 - alpha) * epsilon)
  M <- floor(log(x_sub / T) / h)
  N <- ceiling(log(x_sup / delta) / h)
  idx <- M:(N - 1)
  list(c = h * sin(pi * alpha) / pi * exp((1 - alpha) * h * idx), gamma = exp(h * idx),
       M = M, N = N, delta = delta)
}

.fto_var <- function(name, iv) .make_var(name, "unknown", list(iv), list())

#' Fractional differential equations as ordinary ones
#'
#' `fractional_to_ordinary(eqs, variables, alphas, epsilon, T)` approximates
#' the Caputo fractional equations `D^alpha x = eqs` by systems of ordinary
#' equations, with the exponential sum quadrature of the fractional kernel
#' accurate to `epsilon` on `[0, T]` (as ModelingToolkit). Orders above one
#' use repeated integration. `linear_fractional_to_ordinary(degrees, coeffs,
#' rhs, epsilon, T)` does the same for the linear equation
#' `sum(coeffs * D^degrees x) = rhs`. Both return a compiled system. The
#' auxiliary unknowns are named \ifelse{latex}{\verb{z_i}}{\verb{ʐ_i}}; with `matrix = TRUE` each fractional
#' term gets one array unknown \ifelse{latex}{\verb{z_i}}{\verb{ʐ_i}} (a vector, or for orders above one a
#' matrix with one column per integration) whose elements \ifelse{latex}{\verb{z_i[j]}}{\verb{ʐ_i[j]}} follow
#' the same equations, as the array form of ModelingToolkit.
#'
#' @param eqs Right-hand sides (one per variable).
#' @param variables The fractional variables `x(t)`.
#' @param alphas Fractional orders (not integers).
#' @param epsilon Accuracy of the quadrature.
#' @param T Final time.
#' @param initials Initial values (for orders above one, a vector of initial
#'   derivatives per variable).
#' @param additional_eqs Further equations.
#' @param iv The independent variable.
#' @param matrix Use array unknowns for the quadrature states.
#' @return A compiled `System`.
#' @examples
#' t <- t_nounits
#' x <- Variable("x", t)
#' sys <- fractional_to_ordinary(list(1), list(x), 0.5, 1e-2, 1)
#' length(unknowns(sys))
#' msys <- fractional_to_ordinary(list(1), list(x), 0.5, 1e-2, 1, matrix = TRUE)
#' format(unknowns(msys)[[1]])
#' @export
fractional_to_ordinary <- function(eqs, variables, alphas, epsilon, T, initials = 0,
                                   additional_eqs = list(), iv = t_nounits, matrix = FALSE) {
  D <- Differential(iv)
  if (is_sym(eqs) || is.numeric(eqs)) eqs <- list(eqs)
  variables <- .flatten_vars(variables)
  if (!is.list(initials)) initials <- as.list(initials)
  all_eqs <- list()
  all_def <- symmap()
  cnt <- 0L
  add_eq <- function(lhs, rhs) all_eqs[[length(all_eqs) + 1L]] <<- Equation(lhs, rhs)
  add_def <- function(z) all_def <<- smap_set(all_def, z, 0)
  # Julia zips the equations with `initials`, so a scalar `initials` would
  # convert only the first equation; here a scalar is used for every one
  for (j in seq_along(variables)) {
    sub_eq <- as_sym(eqs[[j]])
    alpha <- alphas[[j]]
    init <- initials[[min(j, length(initials))]]
    if (alpha == 1) stop("`fractional_to_ordinary` requires non-integer orders, got alpha = 1", call. = FALSE)
    coeff <- 1
    alpha0 <- alpha
    m <- 1L
    if (alpha > 1) {
      coeff <- 1 / (alpha - 1)
      m <- 2L
      while (alpha - m > 0) {
        coeff <- coeff / (alpha - m)
        m <- m + 1L
      }
      alpha0 <- alpha - m + 1
    }
    q <- .fto_quadrature(alpha0, epsilon, T)
    idx <- q$M:(q$N - 1)
    nq <- length(idx)
    init_poly <- function() {
      out <- .S0()
      for (k in seq_along(init)) out <- out + init[[k]] * iv^(k - 1) / gamma(k)
      out
    }
    if (matrix) {
      # one array unknown per variable: D(z) ~ -diag(gamma) z + [sub_eq, z[, 1:(m-1)] diag(1:(m-1))]
      nm <- paste0("\u0290_", cnt)
      cnt <- cnt + 1L
      if (alpha < 1) {
        z <- Variable(nm, iv, dims = nq)
        rhs <- .S0()
        for (i in seq_len(nq)) {
          add_eq(D(z[[i]]), -q$gamma[i] * z[[i]] + sub_eq)
          add_def(z[[i]])
          rhs <- rhs + q$c[i] * z[[i]]
        }
        rhs <- rhs + as_sym(init)
      } else {
        z <- Variable(nm, iv, dims = c(nq, m))
        for (k in seq_len(m)) for (i in seq_len(nq)) {
          src <- if (k == 1L) sub_eq else (k - 1) * z[[i, k - 1L]]
          add_eq(D(z[[i, k]]), -q$gamma[i] * z[[i, k]] + src)
          add_def(z[[i, k]])
        }
        rhs <- .S0()
        for (i in seq_len(nq)) rhs <- rhs + q$c[i] * z[[i, m]]
        rhs <- coeff * rhs + init_poly()
      }
    } else if (alpha < 1) {
      rhs <- as_sym(init)
      for (i in seq_len(nq)) {
        z <- .fto_var(paste0("\u0290_", cnt), iv)
        cnt <- cnt + 1L
        add_eq(D(z), sub_eq - q$gamma[i] * z)
        add_def(z)
        rhs <- rhs + q$c[i] * z
      }
    } else {
      rhs <- init_poly()
      for (i in seq_len(nq)) {
        cnt <- cnt + 1L
        base <- sub_eq
        z <- NULL
        for (k in seq_len(m)) {
          # Julia names these by quadrature index and order only, which
          # collides between equations; later equations get a suffix
          z <- .fto_var(paste0("\u0290_", idx[i] - q$M, "_", k, if (j > 1L) paste0("_", j)), iv)
          add_eq(D(z), base - q$gamma[i] * z)
          add_def(z)
          base <- k * z
        }
        rhs <- rhs + coeff * q$c[i] * z
      }
    }
    add_eq(variables[[j]], rhs)
  }
  all_eqs <- c(all_eqs, .as_equations(additional_eqs))
  mtkcompile(System(all_eqs, iv, initial_conditions = all_def, name = "sys"))
}

#' @rdname fractional_to_ordinary
#' @param degrees Orders of the derivatives, decreasing, the first the
#'   highest.
#' @param coeffs Coefficients of the derivatives.
#' @param rhs Right-hand side expression.
#' @param symbol Prefix of the generated variable names.
#' @export
linear_fractional_to_ordinary <- function(degrees, coeffs, rhs, epsilon, T, initials = 0,
                                          symbol = "x", iv = t_nounits, matrix = FALSE) {
  D <- Differential(iv)
  previous <- .fto_var(paste0(symbol, "_0"), iv)
  all_eqs <- list()
  all_def <- symmap()
  cnt <- 0L
  initials <- as.numeric(initials)
  for (i in seq_len(ceiling(degrees[1L]))) {
    new_x <- .fto_var(paste0(symbol, "_", i), iv)
    all_eqs[[length(all_eqs) + 1L]] <- Equation(D(previous), new_x)
    all_def <- smap_set(all_def, previous, initials[min(i, length(initials))])
    previous <- new_x
  }
  new_rhs <- -as_sym(rhs)
  for (j in seq_along(degrees)) {
    degree <- degrees[j]
    rounded <- ceiling(degree)
    new_x <- .fto_var(paste0(symbol, "_", rounded), iv)
    if (degree == round(degree)) {
      new_rhs <- new_rhs + coeffs[j] * new_x
    } else {
      q <- .fto_quadrature(rounded - degree, epsilon, T)
      nq <- length(q$c)
      if (matrix) {
        z <- Variable(paste0("\u0290_", cnt), iv, dims = nq)
        cnt <- cnt + 1L
        zs <- lapply(seq_len(nq), function(i) z[[i]])
      } else {
        zs <- lapply(seq_len(nq), function(i) .fto_var(paste0("\u0290_", cnt + i - 1L), iv))
        cnt <- cnt + nq
      }
      s <- .S0()
      for (i in seq_len(nq)) {
        rhs_i <- if (matrix) -q$gamma[i] * zs[[i]] + new_x else new_x - q$gamma[i] * zs[[i]]
        all_eqs[[length(all_eqs) + 1L]] <- Equation(D(zs[[i]]), rhs_i)
        all_def <- smap_set(all_def, zs[[i]], 0)
        s <- s + q$c[i] * zs[[i]]
      }
      new_rhs <- new_rhs + coeffs[j] * s
    }
  }
  all_eqs[[length(all_eqs) + 1L]] <- Equation(.S0(), new_rhs)
  mtkcompile(System(all_eqs, iv, initial_conditions = all_def, name = "sys"))
}
