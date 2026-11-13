# Fully implicit DAEs F(du, u, p, t) = 0 and consistent initialization,
# after:
#   ModelingToolkitBase/src/problems/daeproblem.jl  DAEFunction, DAEProblem
#     (residual rhs - lhs from generate_rhs(sys; implicit_dae = true), du0
#     from the operating point by default_toterm names, differential_vars)
#   ModelingToolkitBase/src/systems/codegen.jl  generate_dae_jacobian
#   OrdinaryDiffEq.jl lib/OrdinaryDiffEqBDF/src/dae_perform_step.jl
#     DImplicitEuler (SPICE divided difference error) and DFBDF
#   OrdinaryDiffEq.jl lib/OrdinaryDiffEqNonlinearSolve/src/initialize_dae.jl
#     BrownFullBasicInit and ShampineCollocationInit for DAEProblem and for
#     mass matrix ODEProblem, check_dae_tolerance
#   OrdinaryDiffEqDefault/src/default_alg.jl  DFBDF is the DAE default

# ---------------------------------------------------------------------------
# Initialization

#' DAE initialization algorithms
#'
#' Pass as `initializealg` to [solve()] to choose how the integrator makes
#' the initial values of a DAE consistent. `BrownFullBasicInit()` keeps
#' the differential variables and solves the algebraic equations for the
#' algebraic variables (and, for fully implicit problems, for the
#' derivatives of the differential variables). `ShampineCollocationInit()`
#' takes a small implicit Euler step to find consistent values of all
#' variables. `CheckInit()` only checks the residual and fails with retcode
#' "InitialFailure" when it is not satisfied; `NoInit()` takes the initial
#' values as they are, without evaluating the residual. `OverrideInit()`
#' solves the problem's initialization problem (its `initialization_data`,
#' see [OverrideInitData()]) for the initial values and parameters, and
#' `DefaultInit()` is OrdinaryDiffEq's default: `OverrideInit()` when the
#' problem has an initialization problem, else `CheckInit()`. After an event
#' that changes the state of a DAE, `OverrideInit()` and `DefaultInit()`
#' make the algebraic values consistent as `BrownFullBasicInit()` does.
#' Without `initializealg` the solvers use `BrownFullBasicInit()`.
#'
#' @param abstol Tolerance of the nonlinear solve; the integrator's `abstol`
#'   when `NULL`.
#' @param initdt Step of `ShampineCollocationInit()`; by default
#'   `min(t0 / 1000, dtmax / 10)` (or `dtmax / 10` when `t0 = 0`).
#' @param reltol,nlsolve Relative tolerance and nonlinear solver of the
#'   initialization solve of `OverrideInit()` (stored; the initialization
#'   problem is solved with its own algorithm).
#' @return An initialization algorithm object.
#' @examples
#' f <- function(du, u, p, t) c(du[1] + u[1], u[2] - 2 * u[1])
#' prob <- DAEProblem(f, c(0, 0), c(1, 5), c(0, 1), differential_vars = c(TRUE, FALSE))
#' sol <- solve(prob, DFBDF(), initializealg = BrownFullBasicInit())
#' sol$u[[1]]
#' @name dae_initialization
NULL

#' @rdname dae_initialization
#' @export
BrownFullBasicInit <- function(abstol = NULL) structure(list(abstol = min(abstol, 1e-10)), class = c("BrownFullBasicInit", "DAEInitializationAlgorithm"))

#' @rdname dae_initialization
#' @export
ShampineCollocationInit <- function(initdt = NULL) structure(list(initdt = initdt), class = c("ShampineCollocationInit", "DAEInitializationAlgorithm"))

#' @rdname dae_initialization
#' @export
CheckInit <- function() structure(list(), class = c("CheckInit", "DAEInitializationAlgorithm"))

#' @rdname dae_initialization
#' @export
NoInit <- function() structure(list(), class = c("NoInit", "DAEInitializationAlgorithm"))

#' @rdname dae_initialization
#' @export
OverrideInit <- function(abstol = NULL, reltol = NULL, nlsolve = NULL)
  structure(list(abstol = abstol, reltol = reltol, nlsolve = nlsolve),
            class = c("OverrideInit", "DAEInitializationAlgorithm"))

#' @rdname dae_initialization
#' @export
DefaultInit <- function() structure(list(), class = c("DefaultInit", "DAEInitializationAlgorithm"))

# DefaultInit of OrdinaryDiffEqCore's _initialize_dae!
.resolve_initializealg <- function(alg, prob) {
  if (!inherits(alg, "DefaultInit")) return(alg)
  if (has_initializeprob(prob) || has_initializeprob(prob$f)) OverrideInit() else CheckInit()
}

# Make (u, du) consistent at t. Returns list(ok, u, du, nf).
.dae_initialize <- function(ctx, t, u, du, alg, dtmax) {
  alg <- alg %||% BrownFullBasicInit()
  # the initialization problem is solved by the driver; reinitializing after
  # an event keeps the differential values
  if (inherits(alg, c("OverrideInit", "DefaultInit"))) alg <- BrownFullBasicInit()
  abstol <- alg$abstol %||% ctx$abstol
  ok <- function(u, du, nf = 0L) list(ok = TRUE, u = u, du = du, nf = nf)
  if (inherits(alg, "NoInit")) return(ok(u, du))
  fail <- list(ok = FALSE, u = u, du = du, nf = 0L)
  if (isTRUE(ctx$dae)) {
    resid <- function(du, u) ctx$f(du, u, ctx$p, t)
  } else {
    alg_rows <- ctx$algebraic
    if (is.null(alg_rows) || !any(alg_rows)) return(ok(u, du))
    resid <- function(du, u) ctx$f(u, ctx$p, t)[alg_rows]
  }
  # the residual check (CheckInit, check_dae_tolerance) is not counted in nf
  r0 <- resid(du, u)
  if (isTRUE(sqrt(mean(r0^2)) / abstol <= 1)) return(ok(u, du))
  if (inherits(alg, "CheckInit")) return(fail)
  if (inherits(alg, "ShampineCollocationInit")) {
    dt <- alg$initdt %||% (if (t != 0) min(abs(t) / 1000, dtmax / 10) else dtmax / 10)
    M <- ctx$mass
    F <- if (isTRUE(ctx$dae)) function(x) ctx$f((x - u) / dt, x, ctx$p, t) else
      function(x) .mass_mul(M, x - u) / dt - ctx$f(x, ctx$p, t)
    sol <- .safe_newton(F, u, tol = min(abstol, 1e-10))
    if (!sol$converged) return(fail)
    return(ok(sol$x, if (isTRUE(ctx$dae)) (sol$x - u) / dt else du))
  }
  # BrownFullBasicInit
  dv <- ctx$differential_vars
  if (isTRUE(ctx$dae)) {
    if (is.null(dv)) stop("differential_vars must be set for DAE initialization", call. = FALSE)
    F <- function(x) ctx$f(ifelse(dv, x, du), ifelse(dv, u, x), ctx$p, t)
    sol <- .safe_newton(F, ifelse(dv, du, u), tol = min(abstol, 1e-10))
    if (!sol$converged) return(fail)
    return(ok(ifelse(dv, u, sol$x), ifelse(dv, sol$x, du)))
  }
  av <- !dv
  if (sum(av) != sum(alg_rows)) return(fail)
  F <- function(x) {
    uu <- u
    uu[av] <- x
    ctx$f(uu, ctx$p, t)[alg_rows]
  }
  sol <- .safe_newton(F, u[av], tol = min(abstol, 1e-10))
  if (!sol$converged) return(fail)
  u[av] <- sol$x
  ok(u, du)
}

# Jacobians of a residual F(du, u, p, t): prob$f$jac(du, u, p, gamma, t)
# (= dF/du + gamma dF/d(du)) when present, else forward differences.
.calc_J_dae <- function(ctx, du, u, t) {
  if (!is.null(ctx$jac)) {
    J <- as.matrix(ctx$jac(du, u, ctx$p, 0, t))
    return(list(J = J, Jdu = as.matrix(ctx$jac(du, u, ctx$p, 1, t)) - J, nf = 0L))
  }
  F0 <- ctx$f(du, u, ctx$p, t)
  J <- .jac_fd(function(x) ctx$f(du, x, ctx$p, t), u, F0)
  Jdu <- .jac_fd(function(x) ctx$f(x, u, ctx$p, t), du, F0)
  list(J = J, Jdu = Jdu, nf = 2L * length(u) + 1L)
}

# ---------------------------------------------------------------------------
# DAE functions and problems

#' DAE functions and problems
#'
#' `DAEProblem(sys, op, tspan)` builds a fully implicit differential
#' algebraic problem `F(du, u, p, t) = 0` from a complete system, as
#' ModelingToolkit's `DAEProblem`: the residual of each equation is
#' `rhs - lhs`, `op` gives initial values of unknowns, parameters and
#' derivatives (`D(x) ~ 0.5`), and `differential_vars` marks the unknowns
#' whose derivatives appear. `DAEFunction(sys)` builds the residual and,
#' with `jac = TRUE`, the Jacobian `dF/du + gamma * dF/d(du)` as a function
#' of `(du, u, p, gamma, t)`.
#'
#' `DAEProblem(f, du0, u0, tspan, p)` with an R function `f(du, u, p, t)`
#' builds a problem without a system. Solve with [DFBDF()] (the default)
#' or [DImplicitEuler()]; the result is a `DAESolution` with the same
#' interface as an ODE solution.
#'
#' @param sys A complete `System`, or an R function `f(du, u, p, t)`.
#' @param op Operating point.
#' @param tspan Time span.
#' @param du0,u0,p Initial derivatives, initial values and parameters of
#'   raw problems.
#' @param differential_vars Logical vector marking differential variables.
#' @param jac Generate the symbolic Jacobian (for a system); for an R
#'   residual function, the Jacobian `jac(du, u, p, gamma, t)` =
#'   `dF/du + gamma * dF/d(du)` (also accepted as `jac_fn`).
#' @param ... Solver options stored with the problem.
#' @return A `DAEProblem` or `DAEFunction`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); y <- Variable("y", t)
#' sys <- mtkcompile(System(list(D(x) ~ -x + y, y ~ 2 * sin(t)), t, name = "dae"))
#' prob <- DAEProblem(sys, list(x ~ 1, D(x) ~ -1), c(0, 1))
#' sol <- solve(prob, DFBDF())
#' sol[x][length(sol$t)]
#' @export
DAEFunction <- function(sys, jac = FALSE, ...) {
  if (is.function(sys)) {
    J <- if (is.function(jac)) jac else list(...)$jac_fn
    return(structure(list(f = sys, jac = J, sys = NULL, observed = NULL, analytic = list(...)$analytic),
                     class = c("DAEFunction", "SciMLFunction")))
  }
  .need_complete(sys, "DAEFunction")
  eqs <- equations(sys)
  rhs <- lapply(eqs, function(e) if (.is_zero(e$lhs)) e$rhs else e$rhs - e$lhs)
  args <- .codegen_args(sys, du = TRUE)
  f <- .build_closure(rhs, args, observed = .codegen_observed(sys), out = "vector")
  J <- if (isTRUE(jac)) generate_dae_jacobian(sys) else NULL
  structure(list(f = f, jac = J, sys = sys,
                 observed = function(expr) build_explicit_observed_function(sys, expr)),
            class = c("DAEFunction", "SciMLFunction"))
}

#' @rdname DAEFunction
#' @export
DAEProblem <- function(sys, op = list(), tspan = NULL, ..., jac = FALSE, du0 = NULL, u0 = NULL,
                       p = NULL, differential_vars = NULL) {
  dots <- list(...)
  if (is.function(sys) || inherits(sys, "DAEFunction")) {
    f <- if (inherits(sys, "DAEFunction")) sys else DAEFunction(sys)
    # positional form DAEProblem(f, du0, u0, tspan, p)
    if (is.null(u0) && is.numeric(tspan)) {
      u0 <- tspan
      tspan <- dots[[1L]]
      dots <- dots[-1L]
      if (is.null(p) && length(dots) && (is.null(names(dots)) || !nzchar(names(dots)[1L]))) {
        p <- dots[[1L]]
        dots <- dots[-1L]
      }
    }
    if (is.null(du0)) du0 <- op
    u0 <- as.numeric(unlist(u0))
    return(.new_problem("DAEProblem", f = f, du0 = as.numeric(unlist(du0)), u0 = u0,
                        tspan = .as_tspan(tspan), p = if (is.null(p)) numeric() else p, sys = NULL,
                        differential_vars = if (is.null(differential_vars)) NULL else as.logical(differential_vars),
                        kwargs = dots))
  }
  .need_complete(sys, "DAEProblem")
  tspan <- .as_tspan(tspan, sys)
  fn <- DAEFunction(sys, jac = jac)
  vals <- .process_op(sys, op)
  st <- .dae_initial_values(sys, vals)
  prob <- .new_problem("DAEProblem", f = fn, du0 = st$du0, u0 = st$u0, tspan = tspan, p = vals$p, sys = sys,
                       differential_vars = st$differential_vars, kwargs = .problem_kwargs(sys, dots))
  prob$op <- vals$op
  prob
}

# u0, du0 and differential_vars of a DAEProblem from the values resolved by
# .process_op: derivatives from their default_toterm names (0 when not
# given), algebraic unknowns without values from their guesses. Shared with
# remake() (R/problems.R).
.dae_initial_values <- function(sys, vals) {
  us <- unknowns(sys)
  iv <- .f(sys, "iv")
  M <- calculate_massmatrix(sys)
  dvars <- if (length(us)) colSums(abs(as.matrix(M))) != 0 else logical()
  du <- vapply(us, function(v) {
    x <- vals$values[[.key(default_toterm(.make_derivative(v, iv, 1L)))]]
    if (is.numeric(x) && length(x) == 1L && !is.na(x)) x else 0
  }, 0)
  names(du) <- names(vals$u0)
  u0 <- vals$u0
  miss <- which(is.na(u0))
  if (length(miss)) {
    if (any(dvars[miss]))
      stop("no initial value for unknown(s) ", paste(names(u0)[miss[dvars[miss]]], collapse = ", "),
           call. = FALSE)
    g <- guesses(sys)
    for (i in miss) {
      gv <- smap_get(g, us[[i]])
      u0[i] <- if (is.numeric(gv)) gv else 0
    }
  }
  list(u0 = u0, du0 = du, differential_vars = dvars)
}

#' @rdname DAEFunction
#' @export
generate_dae_jacobian <- function(sys, ...) {
  .need_complete(sys, "generate_dae_jacobian")
  eqs <- equations(sys)
  us <- unknowns(sys)
  iv <- .f(sys, "iv")
  obs <- .codegen_observed(sys)
  rule <- if (length(obs)) .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs")) else NULL
  rhs <- lapply(eqs, function(e) {
    r <- if (.is_zero(e$lhs)) e$rhs else e$rhs - e$lhs
    if (is.null(rule)) r else .subst_fixpoint(r, rule)
  })
  dus <- lapply(us, function(u) .make_derivative(u, iv, 1L))
  others <- c(.parameters_of(sys, TRUE), list(iv))
  Ju <- .partials(rhs, us, c(dus, others))
  Jdu <- .partials(rhs, dus, c(us, others))
  gam <- Symbol(".gamma_dae")
  n <- length(us)
  W <- .list_matrix(n, n, function(i, j) Ju[[i, j]] + gam * Jdu[[i, j]])
  args <- .codegen_args(sys, du = TRUE)
  args <- c(args[1:3], list(list(name = "gamma", vars = gam)), args[4L])
  .build_closure(W, args, out = "matrix")
}

#' @rdname solve.ODEProblem
#' @export
solve.DAEProblem <- function(a, b, ...) {
  alg <- if (missing(b) || is.null(b)) DFBDF() else b
  # GlobalRichardson wraps a DAE algorithm (R/globaldiffeq.R)
  if (inherits(alg, "ODEAlgorithm") && is.function(alg$solve_override) && isTRUE(alg$dae))
    return(alg$solve_override(a, alg, utils::modifyList(a$kwargs %||% list(), list(...))))
  if (!inherits(alg, "ODEAlgorithm") || !isTRUE(alg$dae))
    stop("a DAEProblem needs a DAE algorithm such as DFBDF() or DImplicitEuler()", call. = FALSE)
  opts <- utils::modifyList(a$kwargs %||% list(), list(...))
  .check_solve_kwargs(opts, names(list(...)))
  if (!is.null(ns <- .null_solution(a, alg, opts))) {
    class(ns) <- c("DAESolution", class(ns))
    return(ns)
  }
  .integrate(a, alg, opts)
}

#' @export
print.DAESolution <- function(x, ...) {
  env <- .sol_env(x)
  cat("DAESolution with retcode ", env$retcode, " (", env$alg$name %||% "", ")\n", sep = "")
  n <- length(env$t)
  cat("Saved time points: ", n, if (n) paste0(" in [", format(env$t[1L]), ", ", format(env$t[n]), "]"), "\n", sep = "")
  if (n) {
    cat("Final state:\n")
    print(env$U[n, ])
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# DAE algorithms

.dimplicit_euler_init <- function(ctx, t, u) {
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver()
  cache$nf0 <- 0L
  # k[1] of initialize!: the derivative after the DAE initialization
  cache$k1 <- ctx$du
  cache
}

.dimplicit_euler_step <- function(ctx, t, u, dt, cache) {
  nl <- .nlsolve(cache$nl, ctx, t, dt, u * 0, u * 0, 1, 1, t + dt, u, u0 = u)
  out <- list(nf = nl$nf, njacs = nl$njacs, nsolve = nl$nsolve, nw = nl$nw, cache = cache)
  if (!nl$ok) return(c(out, list(fail = TRUE, u = u)))
  unew <- u + nl$z
  EEst <- 1
  if (ctx$adaptive && ctx$success_iter > 0L) {
    tprev <- ctx$tprev
    dt1 <- dt * (t + dt - tprev)
    dt2 <- (t - tprev) * (t + dt - tprev)
    r <- 7 / 12 * dt^2
    EEst <- .err_norm(r * abs((unew - u) / dt1 - (u - ctx$uprev2) / dt2), u, unew, ctx$abstol, ctx$reltol)
  }
  # perform_step! shifts k[1] into k[2] and stores (u - uprev) / dt in k[1],
  # also on steps that are then rejected; the generic hermite_interpolant
  # takes them as the derivatives at the start and the end of the step
  kprev <- cache$k1
  cache$k1 <- (unew - u) / dt
  c(out, list(u = unew, EEst = EEst, du = cache$k1, k = list(cache$k1, kprev, ctx$differential_vars)))
}

#' DAE algorithms
#'
#' Solvers for fully implicit [DAEProblem()]s, as in OrdinaryDiffEqBDF.
#' `DFBDF()` is the variable order (1 to 5) fixed leading coefficient BDF
#' for `F(du, u, p, t) = 0`, the default; `DImplicitEuler()` the implicit
#' Euler method with a divided difference error estimate.
#'
#' @param max_order Maximum order of `DFBDF()`.
#' @param qmax,qsteady_min,qsteady_max Largest step growth factor and the
#'   range of factors that keep the step size.
#' @param autodiff How Jacobians are computed when the problem has no
#'   symbolic Jacobian: `AutoForwardDiff()` (central differences, as R has no
#'   dual numbers), `AutoFiniteDiff(fdtype)` or `AutoSparse(ad)`; see
#'   [AutoForwardDiff()]. The default uses forward differences.
#' @return An algorithm object.
#' @examples
#' f <- function(du, u, p, t) c(du[1] + p * u[1], u[2] - u[1]^2)
#' prob <- DAEProblem(f, c(-1, 0), c(1, 1), c(0, 1), 1, differential_vars = c(TRUE, FALSE))
#' sol <- solve(prob, DFBDF(), reltol = 1e-6)
#' sol(1)
#' @name dae_algorithms
NULL

#' @rdname dae_algorithms
#' @export
DFBDF <- function(max_order = 5L, qmax = 10, qsteady_min = 1, qsteady_max = 2, autodiff = NULL) .with_autodiff(autodiff, {
  .new_bdf_alg("DFBDF", max_order, dae = TRUE, qmax = qmax, qsteady_min = qsteady_min,
               qsteady_max = qsteady_max)
})

#' @rdname dae_algorithms
#' @export
DImplicitEuler <- function(autodiff = NULL) .with_autodiff(autodiff, {
  .new_alg("DImplicitEuler", 1L, step = .dimplicit_euler_step, interp = .dabdf2_interp,
           init = .dimplicit_euler_init, stiff = TRUE, dae = TRUE, controller = "I",
           on_reject = .nl_on_reject)
})
