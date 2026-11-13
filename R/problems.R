# Numerical problems built from systems, after
# ModelingToolkitBase/src/problems/*.jl and systems/problem_utils.jl.
#
# A problem is an S3 list with class c("<Type>Problem", "SciMLProblem"):
#   f       a function object (ODEFunction, NonlinearFunction, ...)
#   u0      named numeric vector of initial values (names are variable names)
#   p       named numeric vector of parameter values
#   tspan   numeric c(t0, t1) for time dependent problems
#   sys     the complete system it was built from (NULL for raw problems)
#   kwargs  solver options stored with the problem (callbacks, tstops, ...)

# Parameters that appear in the numeric parameter vector: bound parameters
# with an expression binding are computed from the others, as observed
# variables are, and are left out.
.problem_params <- function(sys) {
  b <- bindings(sys)
  Filter(function(p) {
    v <- smap_get(b, p)
    is.null(v) || (length(v) == 1L && !is_sym(v) && is.na(v))
  }, .parameters_of(sys, TRUE))
}

# A binding fixes a parameter: an expression of other parameters, or a
# number. NA marks a parameter solved for during initialization instead.
.is_fixed_binding <- function(v) {
  !is.null(v) && (is_sym(v) || (is.numeric(v) && length(v) == 1L && !is.na(v)))
}

# Equations p ~ expr for bound parameters, used as observed in codegen.
.param_binding_eqs <- function(sys) {
  b <- bindings(sys)
  ps <- .parameters_of(sys, TRUE)
  out <- list()
  for (p in ps) {
    v <- smap_get(b, p)
    if (.is_fixed_binding(v)) out[[length(out) + 1L]] <- Equation(p, as_sym(v))
  }
  out
}

.binding_rule <- function(sys) {
  eqs <- .param_binding_eqs(sys)
  .rule(lapply(eqs, `[[`, "lhs"), lapply(eqs, `[[`, "rhs"))
}

# Observed equations including parameter bindings, for generated code.
.codegen_observed <- function(sys) c(.param_binding_eqs(sys), observed(sys))

.codegen_args <- function(sys, with_t = TRUE, du = FALSE) {
  args <- list(list(name = "u", vars = unknowns(sys)),
               list(name = "p", vars = .problem_params(sys)))
  if (with_t) args[[length(args) + 1L]] <- list(name = "t", vars = .f(sys, "iv"))
  if (du) args <- c(list(list(name = "du", vars = lapply(unknowns(sys), function(u)
    .make_derivative(u, .f(sys, "iv"), 1L)))), args)
  args
}

# ---------------------------------------------------------------------------
# Resolving numeric values from operating points, initial conditions and
# bindings.

# Look up names in a symmap whose entries were given by name.
.resolve_named <- function(m, vars) {
  if (!length(m)) return(m)
  out <- symmap()
  by_name <- stats::setNames(vars, vapply(vars, function(v) .name(v) %||% format(v), ""))
  fmt_name <- stats::setNames(vars, vapply(vars, format, ""))
  for (el in unclass(m)) {
    if (!is.null(el$var)) {
      # values given for derivatives refer to their lowered variables
      if (.is_derivative(el$var) || .is_op(el$var, "Shift")) el$var <- default_toterm(el$var)
      out[[.key(el$var)]] <- el
      next
    }
    v <- by_name[[el$name]] %||% fmt_name[[el$name]]
    if (is.null(v)) stop("no variable named ", el$name, " in the system", call. = FALSE)
    out <- smap_set(out, v, el$value)
  }
  class(out) <- "mtk_symmap"
  out
}

# Evaluate a map of values (numbers or expressions of other entries) to
# numbers by fixed point substitution. Returns a named list key -> value
# (NA when unresolved).
.evaluate_varmap <- function(m, extra_rule = NULL, maxiter = 50L) {
  vals <- list()
  pending <- list()
  for (k in names(m)) {
    v <- m[[k]]$value
    if (is_sym(v)) {
      if (.is_num(v)) vals[[k]] <- as.numeric(v) else pending[[k]] <- v
    } else if (is.numeric(v) || is.logical(v)) {
      vals[[k]] <- as.numeric(v)
    } else {
      vals[[k]] <- v
    }
  }
  vars <- lapply(m, `[[`, "var")
  for (it in seq_len(maxiter)) {
    if (!length(pending)) break
    known <- names(vals)[vapply(vals, function(z) is.numeric(z) && length(z) == 1L && !is.na(z), TRUE)]
    rule <- .rule(unname(vars[known]), lapply(known, function(k) as_sym(vals[[k]])))
    if (!is.null(extra_rule)) rule <- .rule_merge(rule, extra_rule)
    progressed <- FALSE
    for (k in names(pending)) {
      e <- .subst(pending[[k]], rule)
      if (.is_num(e)) {
        vals[[k]] <- as.numeric(e)
        pending[[k]] <- NULL
        progressed <- TRUE
      } else if (!length(.collect_vars(e))) {
        num <- tryCatch(as.numeric(N(e)), error = function(err) NA_real_)
        if (!is.na(num)) {
          vals[[k]] <- num
          pending[[k]] <- NULL
          progressed <- TRUE
        } else pending[[k]] <- e
      } else {
        pending[[k]] <- e
      }
    }
    if (!progressed) break
  }
  for (k in names(pending)) vals[[k]] <- NA_real_
  vals
}

# Build the numeric u0 and p vectors for a complete system. Returns a list
# with u0, p (named numeric) and `missing_u` keys of unknowns without value.
.process_op <- function(sys, op, require_p = TRUE) {
  us <- unknowns(sys)
  ps <- .problem_params(sys)
  allv <- c(us, .parameters_of(sys, TRUE), observables(sys))
  op <- .resolve_named(as_symmap(op), allv)
  ics <- initial_conditions(sys)
  binds <- bindings(sys)
  for (k in names(op)) {
    if (!is.null(binds[[k]]) && is_sym(binds[[k]]$value) && isparameter(binds[[k]]$var))
      stop("parameter ", format(binds[[k]]$var), " is bound to ", format(binds[[k]]$value),
           " and cannot be given a value", call. = FALSE)
  }
  merged <- smap_merge(binds, ics, op)
  # observed equations let values of eliminated variables follow from others
  obs <- observed(sys)
  obs_rule <- if (length(obs)) .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs")) else NULL
  for (o in obs) if (!smap_has(merged, o$lhs)) merged <- smap_set(merged, o$lhs, o$rhs)
  vals <- .evaluate_varmap(merged)
  getv <- function(v) {
    x <- vals[[.key(v)]]
    if (is.null(x) || !is.numeric(x)) NA_real_ else x
  }
  u0 <- vapply(us, getv, 0)
  names(u0) <- vapply(us, format, "")
  p <- vapply(ps, getv, 0)
  names(p) <- vapply(ps, format, "")
  if (require_p && anyNA(p)) {
    miss <- names(p)[is.na(p)]
    solvable <- vapply(ps[is.na(p)], function(q) {
      b <- smap_get(binds, q)
      !is.null(b) && !is_sym(b) && is.na(b)
    }, TRUE)
    if (any(!solvable))
      stop("no value for parameter(s) ", paste(miss[!solvable], collapse = ", "),
           "; give them in the operating point or as defaults", call. = FALSE)
  }
  list(u0 = u0, p = p, values = vals, op = op)
}

#' Variable map to numeric vector
#'
#' `varmap_to_vars(varmap, vars)` returns the values of `vars` given in
#' `varmap`, evaluating expressions of other entries.
#'
#' @param varmap A variable map (list of formulas, pairs, or named list).
#' @param vars A list of variables.
#' @param defaults Further values used when `varmap` lacks an entry.
#' @return A named numeric vector.
#' @export
varmap_to_vars <- function(varmap, vars, defaults = NULL) {
  vars <- .flatten_vars(vars)
  m <- smap_merge(.resolve_named(as_symmap(defaults), vars), .resolve_named(as_symmap(varmap), vars))
  vals <- .evaluate_varmap(m)
  out <- vapply(vars, function(v) {
    x <- vals[[.key(v)]]
    if (is.null(x)) NA_real_ else as.numeric(x)
  }, 0)
  stats::setNames(out, vapply(vars, format, ""))
}

.need_complete <- function(sys, what) {
  if (!inherits(sys, "System")) stop(what, " needs a System", call. = FALSE)
  if (!iscomplete(sys))
    stop("a ", what, " can only be built from a complete system; call mtkcompile() or complete() first",
         call. = FALSE)
}

.as_tspan <- function(tspan, sys = NULL) {
  if (is.null(tspan) && !is.null(sys)) tspan <- .f(sys, "tspan")
  if (is.null(tspan)) stop("a time span c(t0, t1) is needed", call. = FALSE)
  tspan <- as.numeric(unlist(tspan))
  if (length(tspan) == 1L) tspan <- c(0, tspan)
  tspan
}

# ---------------------------------------------------------------------------
# ODE functions and problems

#' ODE functions and problems
#'
#' `ODEProblem(sys, op, tspan)` builds an ordinary differential equation
#' problem from a complete system: `op` gives initial values of unknowns and
#' values of parameters (as formulas `x ~ 1`, pairs, or a named list),
#' overriding their defaults. Compiled differential algebraic systems become
#' mass matrix problems. `ODEFunction(sys)` builds the right-hand side,
#' Jacobian and mass matrix without the values.
#'
#' `ODEProblem(f, u0, tspan, p)` with an R function `f(u, p, t)` builds a
#' problem without a system.
#'
#' @param sys A complete `System`, or an R function `f(u, p, t)`.
#' @param op Operating point: values of unknowns and parameters.
#' @param tspan Time span `c(t0, t1)`.
#' @param ... Further values: a legacy parameter map, or solver options
#'   stored with the problem (`callback`, `tstops`, ...). For
#'   `ODEFunction(f, ...)` of an R function: `jac_fn`, `tgrad_fn`,
#'   `mass_matrix`, `analytic` (`function(u0, p, t)`, for
#'   [calculate_solution_errors()]), `paramjac`, `jvp` and `vjp`.
#' @param jac Generate the symbolic Jacobian.
#' @param tgrad Generate the time gradient.
#' @param u0,p Initial values and parameters for raw problems.
#' @param mass_matrix Mass matrix for raw problems.
#' @param check_initialization Check algebraic constraints at t0 and solve
#'   for algebraic unknowns that have no value.
#' @return An `ODEProblem` or `ODEFunction`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)"); parameters("tau")
#' sys <- mtkcompile(System(list(D(x) ~ (1 - x) / tau), t, name = "sys"))
#' prob <- ODEProblem(sys, list(x ~ 0, tau ~ 3), c(0, 10))
#' sol <- solve(prob)
#' sol[x][length(sol$t)]
#' @export
ODEFunction <- function(sys, jac = FALSE, tgrad = FALSE, ...) {
  if (is.function(sys)) {
    dots <- list(...)
    return(structure(list(f = sys, jac = dots$jac_fn %||% dots$jac_fun, tgrad = dots$tgrad_fn,
                          mass_matrix = dots$mass_matrix, sys = NULL, observed = NULL,
                          analytic = dots$analytic, paramjac = dots$paramjac, jvp = dots$jvp, vjp = dots$vjp),
                     class = c("ODEFunction", "SciMLFunction")))
  }
  .need_complete(sys, "ODEFunction")
  obs <- .codegen_observed(sys)
  args <- .codegen_args(sys)
  rhs <- .add_assertions(sys, .ode_rhs_exprs(sys, substitute_observed = FALSE))
  f <- .build_closure(rhs, args, observed = obs, out = "vector")
  J <- NULL
  if (isTRUE(jac)) {
    full <- .ode_rhs_full(sys)
    Jm <- .partials(full, unknowns(sys), c(.parameters_of(sys, TRUE), list(.f(sys, "iv"))))
    J <- .build_closure(Jm, args, out = "matrix")
  }
  tg <- NULL
  if (isTRUE(tgrad)) {
    full <- .ode_rhs_full(sys)
    tm <- .partials(full, list(.f(sys, "iv")), c(unknowns(sys), .parameters_of(sys, TRUE)))
    tg <- .build_closure(lapply(seq_len(nrow(tm)), function(i) tm[[i, 1L]]), args, out = "vector")
  }
  M <- calculate_massmatrix(sys)
  structure(list(f = f, jac = J, tgrad = tg, mass_matrix = M, sys = sys,
                 observed = function(expr) build_explicit_observed_function(sys, expr)),
            class = c("ODEFunction", "SciMLFunction"))
}

# Right-hand sides with observed variables and parameter bindings substituted,
# for differentiation.
.ode_rhs_full <- function(sys) {
  rhs <- .ode_rhs_exprs(sys, substitute_observed = FALSE)
  obs <- .codegen_observed(sys)
  if (!length(obs)) return(rhs)
  rule <- .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs"))
  lapply(rhs, .subst_fixpoint, rule = rule)
}

#' @rdname ODEFunction
#' @export
ODEProblem <- function(sys, op = list(), tspan = NULL, ..., jac = FALSE, tgrad = FALSE,
                       u0 = NULL, p = NULL, mass_matrix = NULL, check_initialization = TRUE) {
  dots <- list(...)
  # ODEProblem(L, u0, tspan) for a SciML operator L: u' = L(u, p, t) u
  if (inherits(sys, "AbstractSciMLOperator")) sys <- .operator_odefunction(sys)
  if (is.function(sys) || inherits(sys, "ODEFunction")) {
    f <- if (inherits(sys, "ODEFunction")) sys else ODEFunction(sys, mass_matrix = mass_matrix)
    # a mass matrix given here applies to an ODEFunction too
    if (!is.null(mass_matrix)) f$mass_matrix <- mass_matrix
    if (is.null(u0)) u0 <- op
    # ODEProblem(f, u0, tspan, p) as in SciMLBase: an unnamed fourth argument is p
    if (is.null(p) && length(dots) && (is.null(names(dots)) || !nzchar(names(dots)[1L]))) {
      p <- dots[[1L]]
      dots <- dots[-1L]
    }
    u0 <- unlist(u0)
    return(.new_problem("ODEProblem", f = f, u0 = stats::setNames(as.numeric(u0), names(u0)),
                        tspan = .as_tspan(tspan), p = if (is.null(p)) numeric() else p,
                        sys = NULL, kwargs = dots))
  }
  .need_complete(sys, "ODEProblem")
  # legacy form ODEProblem(sys, u0map, tspan, pmap)
  legacy_p <- NULL
  if (length(dots) && (is.null(names(dots)) || !nzchar(names(dots)[1L]))) {
    legacy_p <- dots[[1L]]
    dots <- dots[-1L]
  }
  if (!is.null(legacy_p)) op <- c(as.list(op), as.list(legacy_p))
  tspan <- .as_tspan(tspan, sys)
  nlstep <- isTRUE(dots$nlstep)
  dots$nlstep <- NULL
  fn <- ODEFunction(sys, jac = jac, tgrad = tgrad)
  vals <- .process_op(sys, op)
  prob <- .new_problem("ODEProblem", f = fn, u0 = vals$u0, tspan = tspan, p = vals$p,
                       sys = sys, kwargs = .problem_kwargs(sys, dots))
  prob$op <- vals$op
  h <- .mtk$hooks$initialize_problem
  if (!is.null(h)) prob <- h(prob, vals, check = check_initialization)
  else if (check_initialization) prob <- .simple_initialize(prob)
  # the implicit stage system (structural-nlstep.R)
  if (nlstep) prob$f$nlstep_data <- .generate_ode_nlstep_data(sys, prob$u0, prob$p)
  prob
}

.problem_kwargs <- function(sys, dots) {
  h <- .mtk$hooks$problem_callbacks
  cb <- if (!is.null(h)) h(sys) else NULL
  if (!is.null(cb)) dots$callback <- if (is.null(dots$callback)) cb else CallbackSet(cb, dots$callback)
  ts <- .f(sys, "tstops")
  if (length(ts)) dots$sys_tstops <- ts
  dots
}

.new_problem <- function(type, ...) {
  structure(list(...), class = c(type, "SciMLProblem"))
}

# Minimal consistency step used until the initialization module replaces it:
# unknowns without a value are solved for with Newton's method from their
# guesses when the system has algebraic equations, and differential unknowns
# without a value are an error.
.simple_initialize <- function(prob) {
  u0 <- prob$u0
  if (!anyNA(u0)) return(prob)
  sys <- prob$sys
  M <- prob$f$mass_matrix
  alg <- if (is.matrix(M)) which(rowSums(abs(M)) == 0) else integer()
  diffvars <- if (is.matrix(M)) which(colSums(abs(M)) != 0) else seq_along(u0)
  miss <- which(is.na(u0))
  if (any(miss %in% diffvars) || !length(alg))
    stop("no initial value for unknown(s) ", paste(names(u0)[miss], collapse = ", "), call. = FALSE)
  g <- guesses(sys)
  us <- unknowns(sys)
  for (i in miss) {
    gv <- smap_get(g, us[[i]])
    u0[i] <- if (is.numeric(gv)) gv else 0
  }
  f <- prob$f$f
  res <- function(z) {
    u <- u0
    u[miss] <- z
    f(u, prob$p, prob$tspan[1L])[alg]
  }
  sol <- .newton_solve(res, u0[miss])
  if (!sol$converged) stop("could not find consistent initial values for ",
                           paste(names(u0)[miss], collapse = ", "), call. = FALSE)
  u0[miss] <- sol$x
  prob$u0 <- u0
  prob
}

# Damped Newton with a finite difference Jacobian, for small systems.
.newton_solve <- function(F, x0, tol = 1e-10, maxiter = 100L, jac = NULL) {
  x <- as.numeric(x0)
  fx <- F(x)
  for (it in seq_len(maxiter)) {
    if (max(abs(fx)) < tol) return(list(x = x, converged = TRUE, iter = it))
    J <- if (!is.null(jac)) jac(x) else .fd_jacobian(F, x, fx)
    dx <- tryCatch(solve(J, -fx), error = function(e) qr.solve(J, -fx, tol = 1e-14))
    lam <- 1
    nf <- sqrt(sum(fx^2))
    repeat {
      xn <- x + lam * dx
      fn <- F(xn)
      if (all(is.finite(fn)) && sqrt(sum(fn^2)) <= (1 - 1e-4 * lam) * nf) break
      lam <- lam / 2
      if (lam < 1e-8) break
    }
    x <- xn
    fx <- fn
  }
  list(x = x, converged = max(abs(fx)) < sqrt(tol), iter = maxiter)
}

.fd_jacobian <- function(F, x, fx = F(x)) {
  n <- length(x)
  J <- matrix(0, length(fx), n)
  for (j in seq_len(n)) {
    h <- sqrt(.Machine$double.eps) * max(1, abs(x[j]))
    xh <- x
    xh[j] <- x[j] + h
    J[, j] <- (F(xh) - fx) / h
  }
  J
}

#' @export
print.SciMLProblem <- function(x, ...) {
  cat(class(x)[1L], if (!is.null(x$sys)) paste0(" of system ", nameof(x$sys)), "\n", sep = "")
  if (!is.null(x$tspan)) cat("Time span: [", paste(x$tspan, collapse = ", "), "]\n", sep = "")
  if (length(x$u0)) {
    cat("u0:\n")
    print(x$u0)
  }
  if (length(x$p)) {
    cat("p:\n")
    print(x$p)
  }
  invisible(x)
}

#' Rebuild a problem with new values
#'
#' `remake(prob, u0 = , p = , tspan = )` returns a copy of `prob` with new
#' initial values, parameters or time span. `u0` and `p` may be variable
#' maps (only the given entries change) or full numeric vectors.
#'
#' @param prob A problem.
#' @param u0,p,tspan New values.
#' @param ... Further fields to replace.
#' @param check_initialization Run the initialization of the new problem;
#'   with `FALSE` the values are taken as given.
#' @return The new problem.
#' @export
remake <- function(prob, u0 = NULL, p = NULL, tspan = NULL, ..., check_initialization = TRUE) {
  sys <- prob$sys
  if (!is.null(tspan)) prob$tspan <- .as_tspan(tspan)
  if (!is.null(sys) && (!is.null(u0) && !is.numeric(u0) || !is.null(p) && !is.numeric(p))) {
    op <- smap_merge(prob$op %||% symmap(),
                     if (!is.null(u0) && !is.numeric(u0)) .resolve_named(as_symmap(u0), unknowns(sys)),
                     if (!is.null(p) && !is.numeric(p)) .resolve_named(as_symmap(p), .parameters_of(sys, TRUE)))
    vals <- .process_op(sys, op)
    if (inherits(prob, "DAEProblem")) {
      # derivatives D(x) ~ v in the new u0 map replace those entries of du0;
      # the others keep their values, as Julia's coalesce(du0, prob.du0)
      # (R/solvers-dae.R)
      st <- .dae_initial_values(sys, vals)
      # without a new u0 map the current state is kept (prob.u0 in Julia)
      vals$u0 <- if (is.null(u0)) prob$u0 else st$u0
      given <- names(if (!is.null(u0) && !is.numeric(u0)) as_symmap(u0) else symmap())
      iv <- .f(sys, "iv")
      upd <- vapply(unknowns(sys), function(v) {
        d <- .make_derivative(v, iv, 1L)
        any(c(.key(d), .key(default_toterm(d))) %in% given)
      }, TRUE)
      if (length(prob$du0) != length(st$du0)) prob$du0 <- st$du0
      prob$du0[upd] <- st$du0[upd]
    }
    # without initialization, values the map does not determine are kept
    if (!check_initialization && is.numeric(prob$u0) && length(prob$u0) == length(vals$u0))
      vals$u0[is.na(vals$u0)] <- prob$u0[is.na(vals$u0)]
    if (!is.null(u0) && !is.numeric(u0) || is.null(u0)) prob$u0 <- vals$u0
    if (!is.null(p) && !is.numeric(p) || is.null(p)) prob$p <- vals$p
    prob$op <- op
  }
  if (is.numeric(u0)) prob$u0[] <- u0
  if (is.numeric(p)) prob$p[] <- p
  for (n in names(list(...))) prob[[n]] <- list(...)[[n]]
  if (inherits(prob, "DAEProblem") && length(prob$du0) == length(prob$u0)) names(prob$du0) <- names(prob$u0)
  if (!is.null(sys)) {
    h <- .mtk$hooks$initialize_problem
    if (!is.null(h)) prob <- h(prob, NULL, check = check_initialization)
    else if (check_initialization) prob <- .simple_initialize(prob)
  }
  prob
}

# ---------------------------------------------------------------------------
# Symbolic indexing of problems

#' @export
`[.SciMLProblem` <- function(x, i) {
  sys <- x$sys
  if (is.null(sys)) stop("symbolic indexing needs a problem built from a system", call. = FALSE)
  v <- as_sym(i)
  k <- .key(v)
  ui <- match(k, .keys(unknowns(sys)))
  if (!is.na(ui)) return(unname(x$u0[ui]))
  pi <- match(k, .keys(.problem_params(sys)))
  if (!is.na(pi)) return(unname(x$p[pi]))
  fn <- build_explicit_observed_function(sys, v)
  args <- list(x$u0, x$p)
  if (!is.null(.f(sys, "iv"))) args[[3L]] <- x$tspan[1L] %||% 0
  unname(do.call(fn, args))
}

#' @export
`$.SciMLProblem` <- function(x, name) {
  if (name == "ps") return(.param_accessor(x))
  .subset2(x, name)
}

# prob$ps[k] <- v: the replacement on the accessor returns an accessor of the
# updated problem, which this method installs.
#' @export
`$<-.SciMLProblem` <- function(x, name, value) {
  if (identical(name, "ps") && inherits(value, "mtk_param_accessor")) return(.subset2(value, "prob"))
  NextMethod()
}

#' @export
`[<-.mtk_param_accessor` <- function(x, i, value) {
  prob <- .subset2(x, "prob")
  h <- .mtk$hooks$set_parameter
  .param_accessor(if (is.null(h)) remake(prob, p = smap_set(symmap(), as_sym(i), value)) else h(prob, i, value))
}

.param_accessor <- function(prob) {
  structure(list(prob = prob), class = "mtk_param_accessor")
}

#' @export
`[.mtk_param_accessor` <- function(x, i) {
  prob <- .subset2(x, "prob")
  sys <- prob$sys
  v <- as_sym(i)
  pi <- match(.key(v), .keys(.problem_params(sys)))
  if (!is.na(pi)) return(unname(prob$p[pi]))
  b <- smap_get(bindings(sys), v)
  if (.is_fixed_binding(b)) {
    if (!is_sym(b)) return(as.numeric(b))
    rule <- .rule_merge(.rule(.problem_params(sys), lapply(prob$p, as_sym)),
                        .binding_rule(sys))
    return(as.numeric(N(.subst_fixpoint(b, rule))))
  }
  h <- .mtk$hooks$get_parameter
  if (!is.null(h) && !is.null(val <- h(prob, v))) return(val)
  stop(format(v), " is not a parameter of the problem", call. = FALSE)
}
