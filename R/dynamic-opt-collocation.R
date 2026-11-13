# Dynamic optimization (optimal control and parameter estimation), after
# ModelingToolkitBase/src/systems/optimal_control_interface.jl
# (process_DynamicOptProblem, process_tspan, constructDefault,
# extract_variable_bounds, add_observed_bounds!, add_cost_function!,
# get_integral_substitution_rules!, process_integral_bounds,
# add_user_constraints!, add_equational_constraints!, DynamicOptSolution),
# the backends in ModelingToolkitBase/ext (MTKInfiniteOptExt.jl:
# JuMPDynamicOptProblem, InfiniteOptDynamicOptProblem,
# add_solve_constraints!, JuMPCollocation, InfiniteOptCollocation;
# MTKCasADiDynamicOptExt.jl and MTKPyomoDynamicOptExt.jl) and
# docs/src/tutorials/dynamic_optimization.md.
#
# ModelingToolkit transcribes the problem for JuMP, InfiniteOpt, CasADi or
# Pyomo and solves it with Ipopt; none of these exists in R. The backend here
# is native direct collocation on a uniform grid of N nodes with a Runge
# Kutta collocation scheme in stage form: the decision variables are the
# states and inputs at the nodes and at the interior stage points of every
# interval, the tunable parameters and, for free final time problems, the
# duration T (the grid is then on [0, 1] and the dynamics are dx/dtau =
# T f(x, u, p, t0 + T tau)). Each interval has the stage equations
# Y_i = x_k + h T sum_j A_ij f(Y_j) and the update x_{k+1} = x_k + h T
# sum_i b_i f(Y_i); a stage at the start or the end of the interval is the
# node itself. The trapezoidal rule (Lobatto IIIA with two stages), Hermite
# Simpson (Lobatto IIIA with three stages, whose midpoint equation is the
# Hermite interpolation condition), backward Euler, Radau IIA, Gauss
# Legendre and any tableau are available. Integral (Lagrange) costs use the
# quadrature weights b of the scheme, and the solution between nodes is the
# collocation polynomial. An explicit tableau of the JuMP and CasADi
# backends has no stage variables, as in MTK's add_solve_constraints!, and
# an implicit one has derivative stage variables K in place of stage states;
# both integrate costs by the trapezoidal rule on the nodes. The nonlinear
# program is solved by the interior
# point method of R/dynamic-opt-nlp.R with exact first derivatives from
# symbolic Jacobians and second derivatives by differences of them.

#' Abstract dynamic optimization types
#'
#' `AbstractDynamicOptProblem` and `AbstractCollocation` are the class names
#' shared by dynamic optimization problems and collocation methods, as the
#' abstract types of ModelingToolkit: `inherits(prob,
#' AbstractDynamicOptProblem)` is `TRUE` for a [DynamicOptProblem()] and
#' `inherits(alg, AbstractCollocation)` for a [JuMPCollocation()].
#'
#' @format Character strings.
#' @name dynamic_opt_types
NULL

#' @rdname dynamic_opt_types
#' @export
AbstractDynamicOptProblem <- "AbstractDynamicOptProblem"

#' @rdname dynamic_opt_types
#' @export
AbstractCollocation <- "AbstractCollocation"

# ---------------------------------------------------------------------------
# Collocation schemes

# Coefficients (increasing powers) of the Legendre polynomial P_n on [-1, 1].
.dopt_legendre <- function(n) {
  p0 <- 1
  if (n == 0L) return(p0)
  p1 <- c(0, 1)
  for (k in seq_len(n - 1L)) {
    xp <- c(0, p1)
    p2 <- ((2 * k + 1) * xp - k * c(p0, rep(0, length(xp) - length(p0)))) / (k + 1)
    p0 <- p1
    p1 <- p2
  }
  p1
}

.dopt_real_roots <- function(cf) {
  while (length(cf) > 1L && cf[length(cf)] == 0) cf <- cf[-length(cf)]
  if (length(cf) < 2L) return(numeric())
  sort(Re(polyroot(cf)))
}

# Collocation nodes on [0, 1]: Gauss Legendre, right Radau (Radau IIA) and
# Lobatto points.
.dopt_nodes <- function(family, s) {
  x <- switch(family,
    gauss = .dopt_real_roots(.dopt_legendre(s)),
    radau = {
      a <- .dopt_legendre(s); b <- c(.dopt_legendre(s - 1L), 0)
      .dopt_real_roots(a - b[seq_along(a)])
    },
    lobatto = {
      dp <- if (s > 2L) { q <- .dopt_legendre(s - 1L); q[-1L] * seq_len(length(q) - 1L) } else numeric()
      c(-1, .dopt_real_roots(dp), 1)
    })
  x <- (x + 1) / 2
  x[abs(x) < 1e-14] <- 0
  x[abs(x - 1) < 1e-14] <- 1
  x
}

# The collocation tableau of nodes c: A[i, j] is the integral of the
# Lagrange basis polynomial l_j from 0 to c_i, b[j] from 0 to 1; `beta` holds
# the coefficients of the integrals of l_j from 0 to theta, which give the
# collocation polynomial between the nodes.
.dopt_collocation_tableau <- function(cc) {
  s <- length(cc)
  V <- outer(cc, 0:(s - 1L), `^`)
  L <- solve(V)                      # column j: coefficients of l_j
  beta <- t(L / seq_len(s))          # integral of theta^(p-1) is theta^p / p
  pw <- function(x) x^seq_len(s)
  A <- t(vapply(cc, function(ci) as.numeric(beta %*% pw(ci)), numeric(s)))
  list(A = matrix(A, s, s), b = as.numeric(beta %*% pw(1)), c = cc, beta = beta)
}

.dopt_scheme_names <- c("trapezoidal", "hermite_simpson", "backward_euler", "RadauIIA5", "RK4", "BackwardEuler",
                        "MidpointEuler", "Backward", "LagrangeRadau(n)", "LagrangeLegendre(n)",
                        "OrthogonalCollocation(n)")

# A collocation scheme from a name or a tableau list(A, b, c).
.dopt_scheme <- function(spec) {
  if (inherits(spec, "dopt_scheme")) return(spec)
  if (is.list(spec)) {
    A <- as.matrix(spec$A); b <- as.numeric(spec$b); cc <- as.numeric(spec$c %||% rowSums(A))
    if (nrow(A) != length(b) || ncol(A) != length(b) || length(cc) != length(b))
      stop("a tableau needs a square matrix A and vectors b and c of the same size", call. = FALSE)
    tab <- list(A = A, b = b, c = cc, beta = NULL)
    if (!anyDuplicated(round(cc, 12))) {
      ct <- .dopt_collocation_tableau(cc)
      if (max(abs(ct$A - A)) < 1e-10 && max(abs(ct$b - b)) < 1e-10) tab$beta <- ct$beta
    }
    return(structure(c(list(name = spec$name %||% "tableau"), tab), class = "dopt_scheme"))
  }
  if (!is.character(spec) || length(spec) != 1L)
    stop("a collocation method is a name or a tableau list(A = , b = , c = )", call. = FALSE)
  m <- regmatches(spec, regexec("^([A-Za-z]+)\\(([0-9]+)\\)$", spec))[[1L]]
  fam <- if (length(m)) m[2L] else spec
  s <- if (length(m)) as.integer(m[3L]) else NA_integer_
  nodes <- switch(fam,
    trapezoidal = .dopt_nodes("lobatto", 2L),
    hermite_simpson = .dopt_nodes("lobatto", 3L),
    backward_euler = , BackwardEuler = , Backward = 1,
    MidpointEuler = 0.5,
    RadauIIA5 = .dopt_nodes("radau", 3L),
    LagrangeRadau = .dopt_nodes("radau", s),
    LagrangeLegendre = .dopt_nodes("gauss", s),
    OrthogonalCollocation = if (s >= 2L) .dopt_nodes("lobatto", s) else NULL,
    RK4 = NULL,
    stop("unknown collocation method `", spec, "`; use one of ", paste(.dopt_scheme_names, collapse = ", "),
         " or a tableau list(A = , b = , c = )", call. = FALSE))
  if (fam == "RK4") {
    A <- matrix(0, 4L, 4L); A[2L, 1L] <- 0.5; A[3L, 2L] <- 0.5; A[4L, 3L] <- 1
    return(structure(list(name = "RK4", A = A, b = c(1, 2, 2, 1) / 6, c = c(0, 0.5, 0.5, 1), beta = NULL),
                     class = "dopt_scheme"))
  }
  if (is.null(nodes)) stop("`", spec, "` needs at least two collocation points", call. = FALSE)
  structure(c(list(name = spec), .dopt_collocation_tableau(nodes)), class = "dopt_scheme")
}

# ---------------------------------------------------------------------------
# Backends

#' Dynamic optimization backends
#'
#' The problem and solver constructors of the JuMP, InfiniteOpt, CasADi and
#' Pyomo backends of ModelingToolkit, with the same signatures. The Julia
#' packages do not exist in R: each constructor builds the native direct
#' collocation problem of [DynamicOptProblem()] and each collocation object
#' selects the transcription the backend uses, which [solve()] then solves
#' with the native interior point method (the `solver` argument, Ipopt in
#' Julia, is kept but not used).
#'
#' * `JuMPCollocation(solver, tableau)` and `CasADiCollocation(solver,
#'   tableau)` transcribe the dynamics with a Runge Kutta tableau, by
#'   default the fifth order Radau IIA tableau of ModelingToolkit's
#'   `constructDefault()` (`"RadauIIA5"`). The tableau is a name (see
#'   below) or a list `list(A = , b = , c = )`, such as an explicit tableau.
#' * `InfiniteOptCollocation(solver, derivative_method)` uses backward
#'   finite differences by default (`"Backward"`), or
#'   `"OrthogonalCollocation(n)"` (Lobatto points; 2 is the trapezoidal rule
#'   and 3 Hermite Simpson).
#' * `PyomoCollocation(solver, derivative_method)` uses
#'   `"LagrangeRadau(n)"` (Radau IIA, the default with `n = 5`),
#'   `"LagrangeLegendre(n)"` (Gauss Legendre), `"BackwardEuler"` or
#'   `"MidpointEuler"`.
#'
#' The names accepted everywhere are `"trapezoidal"`, `"hermite_simpson"`,
#' `"backward_euler"`, `"RadauIIA5"`, `"RK4"` and the ones above. As in
#' ModelingToolkit, JuMP and CasADi hold the inputs constant over a step
#' (every stage uses the input at the start of the step, and `input_sol` is
#' piecewise constant), ignore `nominal_values`, and Pyomo does not accept
#' `initial_trajectory`; with the other backends and with a name or tableau
#' given to `solve()` the inputs are free at every collocation point.
#' With an explicit tableau (such as `"RK4"`), JuMP and CasADi create no
#' stage variables: each step is one equality per state, from the stages
#' evaluated at the state and input of the step's start, and integral costs
#' use the trapezoidal rule on the nodes. With an implicit tableau (such as
#' the default) they have derivative stage variables K, one per stage and
#' state, with the equations `K_i = f(x_k + dt sum_j A_ij K_j)` and the
#' update `x_{k+1} = x_k + dt sum_i b_i K_i`; JuMP declares K at every node
#' (those of the last node are free), CasADi once per step. Their integral
#' costs also use the trapezoidal rule on the nodes. A successful solve has
#' the retcode `Default` (the retcode of `build_solution`), a failed one
#' `ConvergenceFailure`. The JuMP and InfiniteOpt backends
#' fix the initial states (bounds with equal values, which the interior
#' point method removes from the variables), and CasADi adds them as
#' equality constraints.
#' Further arguments of [DynamicOptProblem()], such as `tune_parameters`,
#' `initial_trajectory`, `nominal_values` and `bounds`, are passed through
#' `...`.
#'
#' @param sys A compiled `System` with inputs, costs and constraints.
#' @param op,tspan,dt,steps,guesses,... As in [DynamicOptProblem()].
#' @param solver The optimizer of the Julia backend, such as `"ipopt"`;
#'   unused.
#' @param tableau,derivative_method The transcription, see above.
#' @return A `DynamicOptProblem` of the backend's class, or a collocation
#'   object of class `AbstractCollocation`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t)
#' u <- Variable("u", t, input = TRUE, bounds = c(-1, 1))
#' sys <- System(list(D(x) ~ u), t, costs = list(-EvalAt(1)(x)), name = "s")
#' sys <- mtkcompile(sys, inputs = list(u))
#' prob <- JuMPDynamicOptProblem(sys, list(x ~ 0, u ~ 0), c(0, 1), dt = 0.1)
#' sol <- solve(prob, JuMPCollocation("ipopt"))
#' tail(sol$sol[x], 1)
#' iprob <- InfiniteOptDynamicOptProblem(sys, list(x ~ 0, u ~ 0), c(0, 1), dt = 0.1)
#' tail(solve(iprob, InfiniteOptCollocation("ipopt"))$sol[x], 1)
#' @name dynamic_opt_backends
NULL

.dopt_backend_problem <- function(backend, sys, op, tspan, dt, steps, guesses, ...) {
  dots <- list(...)
  if (backend %in% c("JuMP", "CasADi")) dots$nominal_values <- NULL
  if (backend == "Pyomo" && length(dots$initial_trajectory))
    stop("`initial_trajectory` is not supported by the Pyomo backend.", call. = FALSE)
  prob <- do.call(DynamicOptProblem, c(list(sys, op, tspan, dt = dt, steps = steps, guesses = guesses), dots))
  prob$backend <- backend
  class(prob) <- c(paste0(backend, "DynamicOptProblem"), class(prob))
  prob
}

#' @rdname dynamic_opt_backends
#' @export
JuMPDynamicOptProblem <- function(sys, op, tspan, dt = NULL, steps = NULL, guesses = list(), ...) {
  .dopt_backend_problem("JuMP", sys, op, tspan, dt, steps, guesses, ...)
}
#' @rdname dynamic_opt_backends
#' @export
InfiniteOptDynamicOptProblem <- function(sys, op, tspan, dt = NULL, steps = NULL, guesses = list(), ...) {
  .dopt_backend_problem("InfiniteOpt", sys, op, tspan, dt, steps, guesses, ...)
}
#' @rdname dynamic_opt_backends
#' @export
CasADiDynamicOptProblem <- function(sys, op, tspan, dt = NULL, steps = NULL, guesses = list(), ...) {
  .dopt_backend_problem("CasADi", sys, op, tspan, dt, steps, guesses, ...)
}
#' @rdname dynamic_opt_backends
#' @export
PyomoDynamicOptProblem <- function(sys, op, tspan, dt = NULL, steps = NULL, guesses = list(), ...) {
  .dopt_backend_problem("Pyomo", sys, op, tspan, dt, steps, guesses, ...)
}

.dopt_collocation <- function(backend, solver, spec) {
  scheme <- .dopt_scheme(spec)
  # JuMP and CasADi evaluate every stage of a step with the input at its start
  scheme$hold <- backend %in% c("JuMP", "CasADi")
  # their add_solve_constraints! take an explicit tableau (MTK.is_explicit)
  # without stage variables, and the InfiniteOpt extension behind JuMP and
  # InfiniteOpt fixes the initial states (JuMP's fix) where CasADi adds
  # equality constraints
  A <- scheme$A
  scheme$explicit <- scheme$hold && all(A[upper.tri(A, diag = TRUE)] == 0)
  # an implicit tableau has derivative stage variables K: JuMP declares them
  # as infinite variables on every support (those of the last support appear
  # in no constraint), CasADi once per step
  scheme$kstage <- scheme$hold && !scheme$explicit
  scheme$k_last <- backend == "JuMP"
  # InfiniteOpt transcribes the derivatives of the states as variables at
  # every support, with (D x - T f) / nominal = 0 there and the relations of
  # its derivative method between D x and x
  scheme$dvar <- backend == "InfiniteOpt"
  scheme$fix_u0 <- backend %in% c("JuMP", "InfiniteOpt")
  structure(list(backend = backend, solver = solver, scheme = scheme),
            class = c(paste0(backend, "Collocation"), AbstractCollocation))
}

#' @rdname dynamic_opt_backends
#' @export
JuMPCollocation <- function(solver = NULL, tableau = NULL) .dopt_collocation("JuMP", solver, tableau %||% "RadauIIA5")
#' @rdname dynamic_opt_backends
#' @export
InfiniteOptCollocation <- function(solver = NULL, derivative_method = NULL) {
  .dopt_collocation("InfiniteOpt", solver, derivative_method %||% "Backward")
}
#' @rdname dynamic_opt_backends
#' @export
CasADiCollocation <- function(solver = NULL, tableau = NULL) .dopt_collocation("CasADi", solver, tableau %||% "RadauIIA5")
#' @rdname dynamic_opt_backends
#' @export
PyomoCollocation <- function(solver = NULL, derivative_method = NULL) {
  .dopt_collocation("Pyomo", solver, derivative_method %||% "LagrangeRadau(5)")
}

#' @export
print.AbstractCollocation <- function(x, ...) {
  cat(class(x)[1L], " (", x$scheme$name, ", ", length(x$scheme$b), " stages)\n", sep = "")
  invisible(x)
}

# The scheme a problem is solved with: the collocation object or name given
# to solve(), or the default of the problem's backend.
.dopt_pick_scheme <- function(prob, b) {
  backend <- .subset2(prob, "backend")
  if (is.null(b)) {
    if (is.null(backend)) return(.dopt_scheme("trapezoidal"))
    b <- switch(backend, JuMP = JuMPCollocation(), CasADi = CasADiCollocation(),
                InfiniteOpt = InfiniteOptCollocation(), Pyomo = PyomoCollocation())
  }
  if (inherits(b, AbstractCollocation)) {
    if (!is.null(backend) && !identical(b$backend, backend))
      stop("a ", backend, "DynamicOptProblem is solved with ", backend, "Collocation(), not ", class(b)[1L], "()",
           call. = FALSE)
    return(b$scheme)
  }
  .dopt_scheme(b)
}

# ---------------------------------------------------------------------------
# Problem construction

.dopt_pairs <- function(x) {
  if (is.null(x) || !length(x)) return(list())
  if (inherits(x, "formula")) x <- list(x)
  lapply(x, function(el) {
    if (inherits(el, "formula")) {
      # the caller's own formula, evaluated where it was written
      env <- environment(el)
      return(list(eval(el[[2L]], env), eval(el[[3L]], env)))
    }
    if (inherits(el, "Equation")) return(list(el$lhs, el$rhs))
    list(as_sym(el[[1L]]), el[[2L]])
  })
}

.dopt_num <- function(v, rule) {
  if (is.null(v)) return(NA_real_)
  if (is.numeric(v)) return(as.numeric(v))
  e <- .subst_fixpoint(as_sym(v), rule)
  if (length(.collect_vars(e))) stop("cannot evaluate ", format(v), " to a number", call. = FALSE)
  as.numeric(N(e))
}

# Inline observed equations (and parameter bindings) in an expression,
# including evaluated forms obs(tau) of observed variables.
.dopt_inline_observed <- function(e, obs, iv) {
  if (!length(obs)) return(e)
  by_name <- list()
  for (o in obs) if (!is.null(.name(o$lhs))) by_name[[.name(o$lhs)]] <- o
  for (it in seq_len(50L)) {
    from <- list(); to <- list()
    for (v in .collect_vars(e)) {
      if (.type(v) == "Derivative" || is.null(.name(v))) next
      o <- by_name[[.name(v)]]
      if (is.null(o)) next
      if (.type(v) == "Symbol" || .key(v) == .key(o$lhs)) {
        rep <- o$rhs
      } else if (.type(v) == "Function" && length(.args(v)) == 1L) {
        tau <- .args(v)[[1L]]
        rep <- .subst(EvalAt(tau)(o$rhs), .rule(list(iv), list(tau)))
      } else next
      from[[length(from) + 1L]] <- v
      to[[length(to) + 1L]] <- rep
    }
    if (!length(from)) return(e)
    e <- .subst(e, .rule(from, to))
  }
  e
}

# Integral(g, (t, a, b)) nodes of an expression.
.dopt_integrals <- function(e) {
  out <- list()
  seen <- new.env(parent = emptyenv())
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    k <- .key(z)
    if (!is.null(seen[[k]])) return(invisible())
    assign(k, TRUE, envir = seen)
    if (.type(z) == "Integral") {
      out[[length(out) + 1L]] <<- z
      return(invisible())
    }
    for (a in .args(z)) walk(a)
  }
  walk(e)
  out
}

# The normalized time of a time point: (tau - t0) / (tf - t0), or 1 for the
# free final time.
.dopt_tau_of <- function(a, ctx, what) {
  if (ctx$free && is_sym(a) && .key(a) == .key(ctx$tf)) return(1)
  if (is_sym(a) && !.is_num(a) && length(.collect_vars(a))) a <- .subst_fixpoint(a, ctx$fixed_rule)
  if (is_sym(a) && length(.collect_vars(a))) stop("cannot evaluate ", what, " in a dynamic optimization problem", call. = FALSE)
  if (ctx$free)
    stop("Provided specific time constraint in a free final time problem. This is not supported by the ",
         "collocation solvers at the moment. The offending variable is ", what, ".", call. = FALSE)
  v <- as.numeric(N(as_sym(a)))
  if (v < ctx$t0 - 1e-10 || v > ctx$t0 + ctx$Tfix + 1e-10)
    stop(what, " is evaluated outside of the time span", call. = FALSE)
  (v - ctx$t0) / ctx$Tfix
}

# Compile an expression of the states and inputs (path values x(t), point
# values x(tau) and x(tf), integrals), tunable parameters, time and the
# duration T into a function and its gradient. Point values become
# placeholders E_s and integrals placeholders I_q; both are resolved on the
# collocation grid by .dopt_resolve().
.dopt_compile <- function(e, ctx, integrand = FALSE) {
  e <- .dopt_inline_observed(as_sym(e), ctx$obs, ctx$iv)
  e <- .subst_fixpoint(e, ctx$fixed_rule)
  ints <- list()
  inodes <- .dopt_integrals(e)
  if (length(inodes)) {
    if (integrand) stop("integrals cannot be nested", call. = FALSE)
    iph <- list()
    for (q in seq_along(inodes)) {
      a <- .args(inodes[[q]])
      lim <- .args(a[[2L]])
      if (length(lim) != 3L || .key(lim[[1L]]) != .key(ctx$iv))
        stop("an integral must run over the independent variable: Integral(expr, list(", format(ctx$iv),
             ", lower, upper))", call. = FALSE)
      if (ctx$free) {
        lo <- lim[[2L]]
        if (!(.is_num(as_sym(lo)) && as.numeric(N(as_sym(lo))) == ctx$t0 && is_sym(lim[[3L]]) &&
              .key(lim[[3L]]) == .key(ctx$tf)))
          stop("Free final time problems cannot handle partial timespans.", call. = FALSE)
        span <- c(0, 1)
      } else {
        span <- vapply(lim[2:3], function(x) {
          x <- .subst_fixpoint(as_sym(x), ctx$fixed_rule)
          if (length(.collect_vars(x))) stop("the bounds of an integral must be numbers", call. = FALSE)
          (as.numeric(N(x)) - ctx$t0) / ctx$Tfix
        }, 0)
        if (span[1L] < -1e-10 || span[2L] > 1 + 1e-10) stop("Integral bounds are beyond the timespan.", call. = FALSE)
      }
      iph[[q]] <- Symbol(paste0(".I", q))
      g <- .dopt_compile(a[[1L]], ctx, integrand = TRUE)
      ints[[q]] <- list(ph = iph[[q]], g = g, span = span)
    }
    e <- .subst(e, .rule(inodes, iph))
  }
  names_u <- vapply(ctx$states, .name, "")
  names_w <- vapply(ctx$ctrls, .name, "")
  ivk <- .key(ctx$iv)
  path <- FALSE
  from <- list(); to <- list(); specs <- list()
  for (v in .collect_vars(e)) {
    if (.type(v) == "Derivative") stop("derivatives cannot appear in costs or constraints", call. = FALSE)
    if (.type(v) == "Symbol") {
      if (.key(v) == ivk) path <- TRUE
      else if (isTRUE(.var_kind(v) %in% c("parameter", "constant")) && !.contains_key(ctx$tun, v) &&
               !(ctx$free && .key(v) == .key(ctx$tf)))
        stop("no value for parameter ", format(v), call. = FALSE)
      next
    }
    if (.type(v) != "Function" || length(.args(v)) != 1L) next
    j <- match(.name(v), names_u); q <- match(.name(v), names_w)
    if (is.na(j) && is.na(q)) next
    a <- .args(v)[[1L]]
    if (.key(a) == ivk) {
      path <- TRUE
      next
    }
    if (integrand) stop("an integrand cannot contain values at given times such as ", format(v), call. = FALSE)
    ph <- Symbol(paste0(".E", length(specs) + 1L))
    from[[length(from) + 1L]] <- v
    to[[length(to) + 1L]] <- ph
    specs[[length(specs) + 1L]] <- list(ph = ph, tau = .dopt_tau_of(a, ctx, format(v)),
                                        col = if (!is.na(j)) j else q, kind = if (!is.na(j)) "u" else "w")
  }
  if (length(from)) e <- .subst(e, .rule(from, to))
  if (ctx$free) e <- .subst(e, .rule(list(ctx$tf), list(ctx$t0 + ctx$Tsym)))
  Eph <- lapply(specs, `[[`, "ph")
  Iph <- lapply(ints, `[[`, "ph")
  args <- list(list(name = "u", vars = ctx$states), list(name = "w", vars = ctx$ctrls),
               list(name = "P", vars = ctx$tun), list(name = "t", vars = ctx$iv),
               list(name = "Ts", vars = ctx$Tsym), list(name = "E", vars = Eph), list(name = "I", vars = Iph))
  wrt <- c(ctx$states, ctx$ctrls, ctx$tun, list(ctx$iv, ctx$Tsym), Eph, Iph)
  G <- .partials(list(e), wrt, list())
  list(f = .build_closure(e, args, out = "scalar"),
       g = .build_closure(lapply(seq_along(wrt), function(i) G[[1L, i]]), args, out = "vector"),
       path = path, specs = specs, ints = ints, expr = e)
}

#' Dynamic optimization problems
#'
#' `DynamicOptProblem(sys, op, tspan, dt = )` turns a compiled system with
#' inputs (see [mtkcompile()]'s `inputs`), costs and constraints into an
#' optimal control problem, transcribed by direct collocation. [solve()]
#' returns a [DynamicOptSolution()]. This is the native backend behind the
#' JuMP, InfiniteOpt, CasADi and Pyomo problems of ModelingToolkit (see
#' [dynamic_opt_backends]).
#'
#' Costs are expressions of the states and inputs evaluated at given times,
#' `EvalAt(1)(x)` (written `x(1)` in Julia), or at the final time, and
#' integral (Lagrange) costs `Integral(expr, list(t, t0, tf))` (written
#' `Symbolics.Integral(t in (t0, tf))(expr)` in Julia), which use the
#' quadrature of the collocation scheme; several costs are summed unless the
#' system has a `consolidate` function. Constraints are equations or
#' inequalities: those with `x(t)` hold at every collocation point, those
#' with evaluated variables at their times. Bounds of states, inputs,
#' observed variables and tunable parameters come from their `bounds`
#' metadata (which may refer to parameters) or from `bounds`. A value at a
#' time between collocation points is interpolated linearly.
#'
#' For a free final time, give `tspan = list(0, tf)` with `tf` a parameter
#' and `steps`; the value of `tf` in `op` is the initial guess.
#'
#' `solve(prob, b)` takes the collocation scheme as `b`: a collocation
#' object of [dynamic_opt_backends], or one of the names `"trapezoidal"`
#' (the default), `"hermite_simpson"`, `"backward_euler"`, `"RadauIIA5"`,
#' `"LagrangeRadau(n)"`, `"LagrangeLegendre(n)"`, `"RK4"`, or a tableau
#' `list(A = , b = , c = )`. The solution between the nodes is the
#' collocation polynomial of the scheme.
#'
#' @param sys A compiled `System`.
#' @param op Operating point: initial values of the states (these become
#'   initial constraints), parameter values, and starting values of the
#'   inputs and of a free final time.
#' @param tspan `c(t0, tf)`, or `list(t0, tf)` with `tf` a parameter.
#' @param dt Grid spacing.
#' @param steps Number of grid points (needed for free final time).
#' @param tune_parameters Also optimize the tunable parameters.
#' @param guesses Starting values of states not fixed by `op`.
#' @param initial_trajectory Starting trajectories, formulas `x ~ expr` of
#'   the independent variable and parameters.
#' @param nominal_values Typical magnitudes of states, formulas `x ~ 10`,
#'   used to scale the collocation constraints.
#' @param bounds Bounds overriding metadata, pairs `list(x, c(lo, hi))`.
#' @param ... Unused.
#' @return An object of class `DynamicOptProblem`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t)
#' u <- Variable("u", t, input = TRUE, bounds = c(-1, 1))
#' sys <- System(list(D(x) ~ u), t, costs = list(-EvalAt(1)(x)), name = "s")
#' sys <- mtkcompile(sys, inputs = list(u))
#' prob <- DynamicOptProblem(sys, list(x ~ 0, u ~ 0), c(0, 1), dt = 0.1)
#' sol <- solve(prob)
#' tail(sol$sol[x], 1)
#' # an integral cost with Hermite Simpson collocation
#' isys <- mtkcompile(System(list(D(x) ~ u), t, name = "q",
#'                           costs = list(Integral(x^2 + u^2, list(t, 0, 1)))), inputs = list(u))
#' isol <- solve(DynamicOptProblem(isys, list(x ~ 1, u ~ 0), c(0, 1), dt = 0.1), "hermite_simpson")
#' isol$model$objective
#' @export
DynamicOptProblem <- function(sys, op = list(), tspan, dt = NULL, steps = NULL, tune_parameters = FALSE,
                              guesses = list(), initial_trajectory = list(), nominal_values = list(),
                              bounds = list(), ...) {
  .need_complete(sys, "DynamicOptProblem")
  iv <- get_iv(sys)
  if (is.null(iv)) stop("dynamic optimization needs a time dependent system", call. = FALSE)
  states <- unknowns(sys); ctrls <- inputs(sys)
  n <- length(states); m <- length(ctrls)
  tspan <- as.list(tspan)
  t0 <- tspan[[1L]]; tfin <- tspan[[2L]]
  if (is_sym(t0) && !.is_num(t0))
    stop("Free initial time problems are not currently supported by the collocation solvers.", call. = FALSE)
  t0 <- as.numeric(t0)
  free <- is_sym(tfin) && !.is_num(tfin)
  if (is.null(dt) && is.null(steps))
    stop("Must provide either the dt or the number of intervals to the collocation solvers.", call. = FALSE)
  if (free) {
    if (is.null(steps))
      stop("Free final time problems require specifying the number of steps using the argument `steps`, rather than dt.",
           call. = FALSE)
    if (!is.null(dt)) warning("Specified dt for free final time problem. This will be ignored; dt will be ",
                              "determined by the number of timesteps.", call. = FALSE)
    Np <- as.integer(steps)
    Tfix <- NA_real_
  } else {
    tfin <- as.numeric(tfin)
    if (!is.null(dt)) {
      if (!is.null(steps)) warning("Specified number of steps for problem with concrete tspan. This will be ",
                                   "ignored; number of steps will be determined by dt.", call. = FALSE)
      Np <- as.integer(floor((tfin - t0) / dt + 1e-9)) + 1L
    } else Np <- as.integer(steps)
    Tfix <- tfin - t0
  }
  if (Np < 2L) stop("the collocation grid needs at least two points", call. = FALSE)
  op <- as_symmap(op)
  opv <- .process_op(sys, op, require_p = FALSE)
  seedv <- .process_op(sys, smap_merge(op, as_symmap(guesses)), require_p = FALSE)
  skeys <- .keys(states)
  # initial values are fixed for the states given in op (all states with a
  # value, for DAEs), from op itself so that guesses only seed the solver
  ou0 <- unname(opv$u0)
  u0_idx <- if (has_alg_equations(sys)) seq_len(n) else which(skeys %in% names(opv$op))
  u0_idx <- u0_idx[!is.na(ou0[u0_idx])]
  pp <- .problem_params(sys)
  pk <- .keys(pp)
  pv <- unname(seedv$p)
  ctrl_idx <- match(.keys(ctrls), pk)
  if (anyNA(ctrl_idx))
    stop("the inputs ", paste(vapply(ctrls[is.na(ctrl_idx)], format, ""), collapse = ", "),
         " are not parameters of the compiled system", call. = FALSE)
  tf_idx <- if (free) match(.key(tfin), pk) else integer()
  if (free && is.na(tf_idx)) stop("the final time ", format(tfin), " is not a parameter of the system", call. = FALSE)
  tun_idx <- if (tune_parameters) setdiff(which(vapply(pp, istunable, TRUE)), c(ctrl_idx, tf_idx)) else integer()
  fixed <- setdiff(seq_along(pp), c(ctrl_idx, tun_idx, tf_idx))
  if (anyNA(pv[fixed]))
    stop("no value for parameter(s) ", paste(vapply(pp[fixed][is.na(pv[fixed])], format, ""), collapse = ", "),
         call. = FALSE)
  fixed_rule <- .rule(pp[fixed], lapply(pv[fixed], as_sym))
  bind <- .param_binding_eqs(sys)
  if (length(bind)) fixed_rule <- .rule_merge(fixed_rule, .rule(lapply(bind, `[[`, "lhs"), lapply(bind, `[[`, "rhs")))
  tun <- pp[tun_idx]
  k <- length(tun)
  Tsym <- Symbol(".Tdur")
  ctx <- list(sys = sys, iv = iv, states = states, ctrls = ctrls, tun = tun, n = n, m = m, k = k,
              N = Np, tau = seq(0, 1, length.out = Np), t0 = t0, Tfix = Tfix, free = free,
              tf = if (free) tfin else NULL, Tsym = Tsym, fixed_rule = fixed_rule, obs = observed(sys))
  # dynamics
  rhs <- lapply(.ode_rhs_full(sys), function(r) .subst_fixpoint(r, fixed_rule))
  if (free) rhs <- lapply(rhs, function(r) .subst(r, .rule(list(tfin), list(t0 + Tsym))))
  M <- calculate_massmatrix(sys)
  M <- if (is.matrix(M)) M else diag(n)
  diff_state <- vapply(seq_len(n), function(r) { j <- which(M[r, ] != 0); if (length(j) == 1L) j else NA_integer_ }, 0L)
  dargs <- list(list(name = "u", vars = states), list(name = "w", vars = ctrls),
                list(name = "P", vars = tun), list(name = "t", vars = iv), list(name = "Ts", vars = Tsym))
  Jm <- .partials(rhs, c(states, ctrls, tun, list(iv, Tsym)), list())
  ctx$fF <- .build_closure(rhs, dargs, out = "vector")
  ctx$fJ <- .build_closure(Jm, dargs, out = "matrix")
  ctx$diff_state <- diff_state
  nom <- .resolve_named(as_symmap(nominal_values), states)
  ctx$scale <- vapply(states, function(s) as.numeric(smap_get(nom, s) %||% getnominal(s)), 0)
  # cost and constraints
  ctx$cost <- .dopt_compile(cost(sys), ctx)
  if (ctx$cost$path)
    stop("costs must be evaluated at given times or integrated; use EvalAt(t1)(x), a free final time or ",
         "Integral(expr, list(t, t0, tf))", call. = FALSE)
  eqc <- list(); inc <- list()
  for (c in constraints(sys)) {
    if (inherits(c, "Inequality")) inc[[length(inc) + 1L]] <- .dopt_compile(.ineq_residual(c), ctx)
    else eqc[[length(eqc) + 1L]] <- .dopt_compile(c$lhs - c$rhs, ctx)
  }
  # bounds: user bounds override metadata
  ub_map <- list()
  for (p in .dopt_pairs(bounds)) ub_map[[.key(p[[1L]])]] <- p[[2L]]
  bnd <- function(v) {
    b <- ub_map[[.key(v)]] %||% (if (hasbounds(v)) getbounds(v) else c(-Inf, Inf))
    c(.dopt_num(b[[1L]], fixed_rule), .dopt_num(b[[2L]], fixed_rule))
  }
  for (o in observed(sys)) {
    b <- if (!is.null(ub_map[[.key(o$lhs)]]) || hasbounds(o$lhs)) bnd(o$lhs) else NULL
    if (is.null(b)) next
    if (is.finite(b[1L])) inc[[length(inc) + 1L]] <- .dopt_compile(b[1L] - o$lhs, ctx)
    if (is.finite(b[2L])) inc[[length(inc) + 1L]] <- .dopt_compile(o$lhs - b[2L], ctx)
  }
  ctx$eqc <- eqc; ctx$inc <- inc
  ctx$bx <- matrix(vapply(states, bnd, numeric(2)), ncol = 2L, byrow = TRUE)
  ctx$bw <- matrix(vapply(ctrls, bnd, numeric(2)), ncol = 2L, byrow = TRUE)
  ctx$bP <- matrix(vapply(tun, bnd, numeric(2)), ncol = 2L, byrow = TRUE)
  ctx$bT <- if (free) {
    b <- if (hasbounds(tfin) || !is.null(ub_map[[.key(tfin)]])) bnd(tfin) - t0 else c(0, Inf)
    c(max(b[1L], 1e-8), b[2L])
  }
  # starting values on the nodes
  u0 <- unname(seedv$u0)
  X0 <- matrix(rep(ifelse(is.na(u0), 0, u0), each = Np), Np, n)
  W0 <- matrix(rep(ifelse(is.na(pv[ctrl_idx]), 0, pv[ctrl_idx]), each = Np), Np, m)
  P0 <- ifelse(is.na(pv[tun_idx]), 0, pv[tun_idx])
  T0 <- if (free) (if (is.na(pv[tf_idx])) 1 else pv[tf_idx] - t0)
  tgrid <- if (free) ctx$tau else t0 + ctx$tau * Tfix
  for (p in .dopt_pairs(initial_trajectory)) {
    j <- match(.key(p[[1L]]), skeys)
    if (is.na(j)) next
    ex <- .dopt_inline_observed(as_sym(p[[2L]]), observed(sys), iv)
    rest <- Filter(function(v) .key(v) != .key(iv) && !.contains_key(pp, v), .collect_vars(ex))
    if (length(rest))
      stop("The initial_trajectory for ", format(p[[1L]]), " may only depend on the independent variable ",
           format(iv), " and on parameters, but it also depends on ",
           paste(vapply(rest, format, ""), collapse = ", "), ".", call. = FALSE)
    fn <- generate_trajectory(sys, ex)
    pnum <- pv
    pnum[is.na(pnum)] <- 0
    X0[, j] <- vapply(tgrid, function(tt) fn(pnum, tt), 0)
  }
  ctx$X0 <- X0; ctx$W0 <- W0; ctx$P0 <- P0; ctx$T0 <- T0
  ctx$u0 <- ou0; ctx$u0_idx <- u0_idx
  ctx$pp <- pp; ctx$pv <- pv; ctx$ctrl_idx <- ctrl_idx; ctx$tun_idx <- tun_idx; ctx$tf_idx <- tf_idx
  if (length(cons <- constraints(sys)) && length(cons) + length(u0_idx) > n)
    warning("The control problem is overdetermined. The total number of conditions (# constraints + # fixed ",
            "initial values given by op) exceeds the total number of states.", call. = FALSE)
  ctx <- .dopt_grid(ctx, .dopt_scheme("trapezoidal"))
  structure(list(sys = sys, ctx = ctx, u0 = seedv$u0, p = seedv$p, tspan = tspan),
            class = c("DynamicOptProblem", AbstractDynamicOptProblem, "SciMLProblem"))
}

#' Compile a trajectory expression
#'
#' `generate_trajectory(sys, expr)` returns an R function `f(p, t)` of the
#' parameter vector of a problem built from `sys` and the independent
#' variable, evaluating `expr` (parameter bindings and observed variables
#' that reduce to parameters are inlined). It compiles the
#' `initial_trajectory` guesses of [DynamicOptProblem()].
#'
#' @param sys A complete time dependent `System`.
#' @param expr An expression of the independent variable and parameters.
#' @param ... Unused.
#' @return A function of `p` and `t`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); a <- Parameter("a")
#' sys <- mtkcompile(System(list(D(x) ~ -a * x), t, name = "s"))
#' f <- generate_trajectory(sys, a * t^2)
#' f(2, 3)
#' @export
generate_trajectory <- function(sys, expr, ...) {
  args <- list(list(name = "p", vars = .problem_params(sys)), list(name = "t", vars = get_iv(sys)))
  single <- is_sym(expr) || is.numeric(expr)
  .build_closure(if (single) as_sym(expr) else lapply(expr, as_sym), args,
                 observed = .codegen_observed(sys), out = if (single) "scalar" else "vector")
}

#' @export
print.DynamicOptProblem <- function(x, ...) {
  ctx <- .subset2(x, "ctx")
  cat(class(x)[1L], " of system ", nameof(.subset2(x, "sys")), "\n", sep = "")
  cat("Grid points: ", ctx$N, if (ctx$free) " (free final time)", "\n", sep = "")
  cat("States: ", ctx$n, ", inputs: ", ctx$m, ", tunable parameters: ", ctx$k, "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# The collocation grid

# Evaluation points of a scheme (the nodes and the interior stage points of
# every interval, in time order), the layout of the decision vector
# (states, inputs, tuned parameters, duration), bounds, the starting point,
# and the point values and integrals of the compiled expressions resolved
# on the points.
.dopt_grid <- function(ctx, scheme) {
  N <- ctx$N; n <- ctx$n; m <- ctx$m; k <- ctx$k
  A <- scheme$A; b <- scheme$b; cc <- scheme$c
  s <- length(b)
  # an explicit tableau of JuMP or CasADi: the stages are expressions of the
  # state and input at the start of the step (one equality per step and
  # state), so the nodes are the only points
  explicit <- isTRUE(scheme$explicit) && !anyNA(ctx$diff_state)
  # an implicit tableau of JuMP or CasADi: derivative stage variables K
  # (MTK's add_solve_constraints!), so the nodes are again the only points
  kstage <- isTRUE(scheme$kstage) && !anyNA(ctx$diff_state)
  start <- if (explicit || kstage) rep(TRUE, s) else cc == 0 & apply(A == 0, 1L, all)
  endst <- abs(cc - 1) < 1e-12 & apply(abs(sweep(A, 2L, b)) < 1e-12, 1L, all) & !start
  interior <- !(start | endst)
  ni <- sum(interior)
  node_pt <- (seq_len(N) - 1L) * (ni + 1L) + 1L
  Np <- node_pt[N]
  dtau <- 1 / (N - 1L)
  spt <- matrix(0L, N - 1L, s)
  ptau <- numeric(Np)
  ptau[node_pt] <- ctx$tau
  rank <- cumsum(interior)
  for (i in seq_len(s)) {
    spt[, i] <- if (start[i]) node_pt[-N] else if (endst[i]) node_pt[-1L] else node_pt[-N] + rank[i]
    if (interior[i]) ptau[spt[, i]] <- ctx$tau[-N] + cc[i] * dtau
  }
  # inputs at a node that is no stage point (Gauss and backward Euler
  # schemes) are tied to the input polynomial of the adjacent interval
  used <- unique(as.integer(spt))
  ties <- list()
  for (nd in seq_len(N)) {
    if (node_pt[nd] %in% used) next
    kk <- if (nd > 1L) nd - 1L else 1L
    sp <- unique(spt[kk, ])
    th <- (ptau[sp] - ptau[node_pt[kk]]) / dtau
    th0 <- if (nd > 1L) 1 else 0
    w <- if (anyDuplicated(round(th, 12))) rep(1 / length(th), length(th))
         else vapply(seq_along(th), function(i) prod((th0 - th[-i]) / (th[i] - th[-i])), 0)
    ties[[length(ties) + 1L]] <- list(pt = node_pt[nd], sp = sp, w = w)
  }
  # inputs held over a step (the JuMP and CasADi transcriptions): the inputs
  # of the interior stages equal the input at the start of the interval, and
  # the stage at the end of an interval is evaluated with that input
  hold <- isTRUE(scheme$hold) && m > 0L && !explicit && !kstage
  if (hold) {
    ties <- list()
    for (i in which(interior)) for (kk in seq_len(N - 1L))
      ties[[length(ties) + 1L]] <- list(pt = spt[kk, i], sp = node_pt[kk], w = 1)
  }
  ctx$hold <- hold
  ctx$explicit <- explicit
  ctx$kstage <- kstage
  ctx$stage_end <- endst
  dvar <- isTRUE(scheme$dvar) && !explicit && !kstage
  ctx$dvar <- dvar
  ctx$ties <- if (m && !explicit && !kstage && !dvar) ties else list()
  ctx$scheme <- scheme
  ctx$Npts <- Np; ctx$node_pt <- node_pt; ctx$spt <- spt; ctx$ptau <- ptau; ctx$dtau <- dtau
  ctx$stage_interior <- interior
  # the stage variables K[kk, i, j] of nK nodes follow the duration, node by
  # node and stage by stage; they start at 0 (no start value) and are free
  nK <- if (kstage) (if (isTRUE(scheme$k_last)) N else N - 1L) else 0L
  ctx$nK <- nK
  ctx$iT <- Np * (n + m) + k + 1L
  ctx$kbase <- Np * (n + m) + k + as.integer(ctx$free)
  # the derivative variables of InfiniteOpt: D[p, q] for point p and the
  # q-th differential state, column by column
  ctx$drows <- which(!is.na(ctx$diff_state))
  nD <- if (dvar) Np * length(ctx$drows) else 0L
  nx <- nK * s * n + nD
  ctx$nz <- ctx$kbase + nx
  ctx$var_key <- c(rep(seq_len(Np), n), rep(seq_len(Np), m), rep(Inf, k + as.integer(ctx$free)),
                   rep(node_pt[seq_len(nK)], each = s * n), rep(seq_len(Np), length.out = nD))
  ctx$lower <- c(rep(ctx$bx[, 1L], each = Np), rep(ctx$bw[, 1L], each = Np), ctx$bP[, 1L], ctx$bT[1L],
                 rep(-Inf, nx))
  ctx$upper <- c(rep(ctx$bx[, 2L], each = Np), rep(ctx$bw[, 2L], each = Np), ctx$bP[, 2L], ctx$bT[2L],
                 rep(Inf, nx))
  if (dvar) {
    # the relations of the derivative method on each finite element
    # (derivative_expr_data of InfiniteOpt): sum_k C[k, j] D(q_k) - x(q_j) +
    # x(lb) = 0 for the points q after the element's start lb, with
    # C = M1t \ M2t from the support offsets (Backward: C = dt)
    sc_t <- if (ctx$free) 1 else ctx$Tfix
    ep <- order(cc)[!start[order(cc)]]
    off <- cc[ep] * dtau * sc_t
    nq <- length(off)
    M1t <- outer(seq_len(nq), seq_len(nq), function(kq, j) kq * off[j]^(kq - 1L))
    M2t <- outer(seq_len(nq), seq_len(nq), function(kq, j) off[j]^kq)
    ctx$dv_pts <- ep
    ctx$dv_C <- solve(M1t, M2t)
  }
  spread <- function(M) {
    if (!ncol(M)) return(numeric())
    as.numeric(vapply(seq_len(ncol(M)), function(j) stats::approx(ctx$tau, M[, j], xout = ptau, rule = 2)$y,
                      numeric(Np)))
  }
  z0 <- c(spread(ctx$X0), spread(ctx$W0), ctx$P0, ctx$T0, numeric(nx))
  # fixed initial states are bounds with equal values, which Ipopt turns
  # into parameters (fixed_variable_treatment = make_parameter)
  ctx$u0_fixed <- isTRUE(scheme$fix_u0)
  if (ctx$u0_fixed && length(ctx$u0_idx)) {
    iz <- (ctx$u0_idx - 1L) * Np + node_pt[1L]
    ctx$lower[iz] <- ctx$u0[ctx$u0_idx]
    ctx$upper[iz] <- ctx$u0[ctx$u0_idx]
  }
  ctx$z0 <- pmin(pmax(z0, ctx$lower), ctx$upper)
  # point values are interpolated between the collocation points of a
  # collocation method, between the nodes otherwise
  ip <- if (is.null(scheme$beta)) node_pt else seq_len(Np)
  ctx$interp_pts <- ip[order(ptau[ip])]
  ctx$cost <- .dopt_resolve(ctx$cost, ctx)
  ctx$eqc <- lapply(ctx$eqc, .dopt_resolve, ctx = ctx)
  ctx$inc <- lapply(ctx$inc, .dopt_resolve, ctx = ctx)
  ctx
}

.dopt_resolve <- function(it, ctx) {
  ip <- ctx$interp_pts
  tt <- ctx$ptau[ip]
  for (s in seq_along(it$specs)) {
    tau <- it$specs[[s]]$tau
    hit <- which(abs(tt - tau) <= 1e-10)
    if (length(hit)) {
      it$specs[[s]]$pts <- ip[hit[1L]]
      it$specs[[s]]$wts <- 1
    } else {
      r <- max(1L, min(findInterval(tau, tt), length(tt) - 1L))
      th <- (tau - tt[r]) / (tt[r + 1L] - tt[r])
      it$specs[[s]]$pts <- ip[c(r, r + 1L)]
      it$specs[[s]]$wts <- c(1 - th, th)
    }
  }
  for (q in seq_along(it$ints)) {
    span <- it$ints[[q]]$span
    kl <- which(abs(ctx$tau - span[1L]) <= 1e-9); kh <- which(abs(ctx$tau - span[2L]) <= 1e-9)
    if (!length(kl) || !length(kh))
      stop("the bounds of an integral must be times of the collocation grid; adjust dt", call. = FALSE)
    ks <- seq_len(ctx$N - 1L)
    ks <- ks[ks >= kl[1L] & ks < kh[1L]]
    if (isTRUE(ctx$dvar)) {
      # InfiniteOpt's UniTrapezoid on all supports in the span, internal
      # collocation points included (_trapezoid_coeff)
      ord <- order(ctx$ptau)
      pts <- ord[ctx$ptau[ord] >= span[1L] - 1e-9 & ctx$ptau[ord] <= span[2L] + 1e-9]
      tq <- ctx$ptau[pts]
      L <- length(tq)
      wq <- c((tq[2L] - tq[1L]) / 2, if (L > 2L) (tq[3:L] - tq[1:(L - 2L)]) / 2, (tq[L] - tq[L - 1L]) / 2)
    } else if (isTRUE(ctx$explicit) || isTRUE(ctx$kstage)) {
      # the trapezoidal rule on the supports (InfiniteOpt's default measure
      # and the lowered_integral of CasADi)
      pts <- c(ctx$node_pt[ks], ctx$node_pt[ks + 1L])
      wq <- rep(ctx$dtau / 2, 2L * length(ks))
    } else {
      pts <- as.integer(ctx$spt[ks, , drop = FALSE])
      wq <- rep(ctx$dtau * ctx$scheme$b, each = length(ks))
    }
    agg <- rowsum(wq, pts)
    it$ints[[q]]$qpts <- as.integer(rownames(agg))
    it$ints[[q]]$qw <- agg[, 1L]
  }
  it
}

# ---------------------------------------------------------------------------
# Evaluation of the nonlinear program

# Columns of the decision vector for the variables of a point (states,
# inputs, tuned parameters, duration).
.dopt_point_cols <- function(ctx, pt) {
  Np <- ctx$Npts
  c((seq_len(ctx$n) - 1L) * Np + pt, Np * ctx$n + (seq_len(ctx$m) - 1L) * Np + pt,
    Np * (ctx$n + ctx$m) + seq_len(ctx$k), if (ctx$free) ctx$iT)
}

# The gradient over (states, inputs, parameters, t, T) of a function at a
# point, as columns and values of the decision vector.
.dopt_local_to_z <- function(ctx, pt, g) {
  nmk <- ctx$n + ctx$m + ctx$k
  cols <- .dopt_point_cols(ctx, pt)
  vals <- g[seq_len(nmk)]
  if (ctx$free) vals <- c(vals, g[nmk + 1L] * ctx$ptau[pt] + g[nmk + 2L])
  list(cols = cols, vals = vals)
}

# Hessian of a scalar function of the variables of point pt by forward
# differences of its gradient grad(x, u, P, Tv) over (x, u, P, T); the
# inputs may be taken from another point ptu (inputs held over a step).
.dopt_point_hess <- function(ctx, pt, X, W, P, Tv, grad, ptu = pt) {
  n <- ctx$n; m <- ctx$m; k <- ctx$k; Np <- ctx$Npts
  v <- c(X[pt, ], W[ptu, ], P, if (ctx$free) Tv)
  d <- length(v)
  gfun <- function(v) grad(v[seq_len(n)], v[n + seq_len(m)], v[n + m + seq_len(k)], if (ctx$free) v[d] else Tv)
  H <- .dopt_fd_hess(gfun, v)
  cols <- c((seq_len(n) - 1L) * Np + pt, Np * n + (seq_len(m) - 1L) * Np + ptu,
            Np * (n + m) + seq_len(k), if (ctx$free) ctx$iT)
  nzh <- which(H != 0, arr.ind = TRUE)
  list(i = cols[nzh[, 1L]], j = cols[nzh[, 2L]], v = H[nzh])
}

# The symmetric Hessian of a function with gradient g by forward differences.
.dopt_fd_hess <- function(g, v) {
  d <- length(v)
  H <- matrix(0, d, d)
  if (!d) return(H)
  g0 <- g(v)
  for (q in seq_len(d)) {
    h <- 1e-7 * max(1, abs(v[q]))
    vq <- v
    vq[q] <- v[q] + h
    H[, q] <- (g(vq) - g0) / h
  }
  (H + t(H)) / 2
}

# Values, gradients (triplets with local rows) and Hessian data of a
# compiled expression: one row per point for path expressions, one row
# otherwise.
.dopt_expr <- function(it, ctx, X, W, P, Tv, tt, level) {
  n <- ctx$n; m <- ctx$m; k <- ctx$k; Np <- ctx$Npts; nmk <- n + m + k
  E <- vapply(it$specs, function(s) {
    V <- if (s$kind == "u") X[, s$col] else W[, s$col]
    sum(V[s$pts] * s$wts)
  }, 0)
  Ival <- numeric(length(it$ints))
  Igrad <- vector("list", length(it$ints))
  for (q in seq_along(it$ints)) {
    iq <- it$ints[[q]]
    g <- iq$g
    gv <- vapply(iq$qpts, function(pt) g$f(X[pt, ], W[pt, ], P, tt[pt], Tv, numeric(), numeric()), 0)
    Ival[q] <- Tv * sum(iq$qw * gv)
    if (level >= 1L) {
      cols <- list(); vals <- list()
      for (r in seq_along(iq$qpts)) {
        pt <- iq$qpts[r]
        gr <- g$g(X[pt, ], W[pt, ], P, tt[pt], Tv, numeric(), numeric())
        lz <- .dopt_local_to_z(ctx, pt, Tv * iq$qw[r] * gr)
        if (ctx$free) lz$vals[length(lz$vals)] <- lz$vals[length(lz$vals)] + iq$qw[r] * gv[r]
        cols[[r]] <- lz$cols; vals[[r]] <- lz$vals
      }
      Igrad[[q]] <- list(cols = unlist(cols), vals = unlist(vals))
    }
  }
  pts <- if (it$path) seq_len(Np) else 1L
  v <- vapply(pts, function(pt) it$f(X[pt, ], W[pt, ], P, tt[pt], Tv, E, Ival), 0)
  out <- list(v = v, E = E, I = Ival, Igrad = Igrad, pts = pts)
  # rows are placed after the last point they depend on
  dep <- c(unlist(lapply(it$specs, `[[`, "pts")), unlist(lapply(it$ints, `[[`, "qpts")))
  out$key <- if (it$path) pmax(pts, max(c(0, dep))) else max(c(0, dep))
  if (level < 1L) return(out)
  ti <- list(); tj <- list(); tv <- list()
  ns <- length(it$specs)
  for (r in seq_along(pts)) {
    pt <- pts[r]
    gr <- it$g(X[pt, ], W[pt, ], P, tt[pt], Tv, E, Ival)
    if (it$path) {
      lz <- .dopt_local_to_z(ctx, pt, gr)
      cols <- lz$cols; vals <- lz$vals
    } else {
      cols <- c(Np * (n + m) + seq_len(k), if (ctx$free) ctx$iT)
      vals <- c(gr[n + m + seq_len(k)], if (ctx$free) gr[nmk + 2L])
    }
    for (s in seq_len(ns)) {
      sp <- it$specs[[s]]
      base <- if (sp$kind == "u") (sp$col - 1L) * Np else Np * n + (sp$col - 1L) * Np
      cols <- c(cols, base + sp$pts)
      vals <- c(vals, gr[nmk + 2L + s] * sp$wts)
    }
    for (q in seq_along(it$ints)) {
      cols <- c(cols, Igrad[[q]]$cols)
      vals <- c(vals, gr[nmk + 2L + ns + q] * Igrad[[q]]$vals)
    }
    ti[[r]] <- rep(r, length(cols)); tj[[r]] <- cols; tv[[r]] <- vals
  }
  out$J <- list(i = unlist(ti), j = unlist(tj), v = unlist(tv))
  out
}

# Hessian triplets of sum(w * e) for a compiled expression with weights w
# per row. A row depends on the variables of its point (path expressions),
# the tuned parameters, the duration and the point values E and integrals
# I; second derivatives in all of them are mapped to the decision vector by
# the chain rule, and the curvature of the integrands is added at their
# quadrature points.
.dopt_expr_hess <- function(it, ex, w, ctx, X, W, P, Tv, tt) {
  n <- ctx$n; m <- ctx$m; k <- ctx$k; Np <- ctx$Npts; nmk <- n + m + k
  ns <- length(it$specs); nq <- length(it$ints)
  free <- ctx$free
  out <- list()
  spec_maps <- lapply(it$specs, function(sp) {
    base <- if (sp$kind == "u") (sp$col - 1L) * Np else Np * n + (sp$col - 1L) * Np
    list(cols = base + sp$pts, vals = sp$wts)
  })
  for (r in which(w != 0)) {
    wr <- w[r]
    pt <- ex$pts[r]
    if (it$path && !ns && !nq) {
      # a function of one point only
      out[[length(out) + 1L]] <- .dopt_point_hess(ctx, pt, X, W, P, Tv, function(x, u, Pp, Tt) {
        tp <- ctx$t0 + ctx$ptau[pt] * Tt
        gr <- it$g(x, u, Pp, tp, Tt, ex$E, ex$I)
        c(wr * gr[seq_len(nmk)], if (free) wr * (gr[nmk + 1L] * ctx$ptau[pt] + gr[nmk + 2L]))
      })
      next
    }
    loc <- if (it$path) n + m else 0L
    a0 <- c(if (it$path) c(X[pt, ], W[pt, ]), P, if (free) Tv, ex$E, ex$I)
    ga <- function(a) {
      x <- if (it$path) a[seq_len(n)] else X[pt, ]
      u <- if (it$path) a[n + seq_len(m)] else W[pt, ]
      Pp <- a[loc + seq_len(k)]
      Tt <- if (free) a[loc + k + 1L] else Tv
      off <- loc + k + as.integer(free)
      tp <- ctx$t0 + ctx$ptau[pt] * Tt
      gr <- it$g(x, u, Pp, tp, Tt, a[off + seq_len(ns)], a[off + ns + seq_len(nq)])
      wr * c(if (it$path) gr[seq_len(n + m)], gr[n + m + seq_len(k)],
             if (free) gr[nmk + 1L] * ctx$ptau[pt] + gr[nmk + 2L], gr[nmk + 2L + seq_len(ns + nq)])
    }
    g0 <- ga(a0)
    Ha <- .dopt_fd_hess(ga, a0)
    maps <- c(if (it$path) c(lapply(seq_len(n), function(j) list(cols = (j - 1L) * Np + pt, vals = 1)),
                             lapply(seq_len(m), function(q) list(cols = Np * n + (q - 1L) * Np + pt, vals = 1))),
              lapply(seq_len(k), function(l) list(cols = Np * (n + m) + l, vals = 1)),
              if (free) list(list(cols = ctx$iT, vals = 1)), spec_maps, ex$Igrad)
    nzh <- which(Ha != 0, arr.ind = TRUE)
    for (e in seq_len(nrow(nzh))) {
      ma <- maps[[nzh[e, 1L]]]; mb <- maps[[nzh[e, 2L]]]
      out[[length(out) + 1L]] <- list(i = rep(ma$cols, times = length(mb$cols)),
                                      j = rep(mb$cols, each = length(ma$cols)),
                                      v = Ha[nzh[e, 1L], nzh[e, 2L]] * rep(ma$vals, times = length(mb$vals)) *
                                        rep(mb$vals, each = length(ma$vals)))
    }
    # curvature of the integrals: sum over the quadrature points of T g
    for (q in seq_len(nq)) {
      coef <- g0[length(g0) - nq + q]
      if (coef == 0) next
      iq <- it$ints[[q]]
      for (rr in seq_along(iq$qpts)) {
        qp <- iq$qpts[rr]
        om <- coef * iq$qw[rr]
        out[[length(out) + 1L]] <- .dopt_point_hess(ctx, qp, X, W, P, Tv, function(x, u, Pp, Tt) {
          tp <- ctx$t0 + ctx$ptau[qp] * Tt
          gr <- iq$g$g(x, u, Pp, tp, Tt, numeric(), numeric())
          g <- om * Tt * gr[seq_len(nmk)]
          if (free) {
            gv <- iq$g$f(x, u, Pp, tp, Tt, numeric(), numeric())
            g <- c(g, om * (gv + Tt * (gr[nmk + 1L] * ctx$ptau[qp] + gr[nmk + 2L])))
          }
          g
        })
      }
    }
  }
  out
}

# Values and derivatives of the nonlinear program at z: the cost f, its
# gradient gf, the equality constraints cE (collocation equations,
# algebraic equations, initial values, user equations) and inequalities cI
# with Jacobians as triplets (and as matrices with dense = TRUE), the keys
# placing the rows in time, and with level 2 the Hessian of the Lagrangian
# obj f + lamE' cE + lamI' cI.
.dopt_eval <- function(ctx, z, level = 1L, dense = TRUE, lamE = NULL, lamI = NULL, obj = 1) {
  Np <- ctx$Npts; N <- ctx$N; n <- ctx$n; m <- ctx$m; k <- ctx$k; nz <- ctx$nz
  nc <- n + m + k + 2L
  X <- matrix(z[seq_len(Np * n)], Np, n)
  W <- matrix(z[Np * n + seq_len(Np * m)], Np, m)
  P <- z[Np * (n + m) + seq_len(k)]
  Tv <- if (ctx$free) z[ctx$iT] else ctx$Tfix
  tt <- ctx$t0 + ctx$ptau * Tv
  FF <- matrix(0, Np, n)
  for (i in seq_len(Np)) FF[i, ] <- ctx$fF(X[i, ], W[i, ], P, tt[i], Tv)
  if (level >= 1L) {
    JJ <- array(0, c(Np, n, nc))
    for (i in seq_len(Np)) JJ[i, , ] <- ctx$fJ(X[i, ], W[i, ], P, tt[i], Tv)
  }
  hold <- isTRUE(ctx$hold)
  if (hold) {
    # the end stage of each interval with the input at its start
    hx <- ctx$node_pt[-1L]; hu <- ctx$node_pt[-N]
    FH <- matrix(0, N - 1L, n)
    for (i in seq_len(N - 1L)) FH[i, ] <- ctx$fF(X[hx[i], ], W[hu[i], ], P, tt[hx[i]], Tv)
    if (level >= 1L) {
      JH <- array(0, c(N - 1L, n, nc))
      for (i in seq_len(N - 1L)) JH[i, , ] <- ctx$fJ(X[hx[i], ], W[hu[i], ], P, tt[hx[i]], Tv)
    }
  }
  ti <- list(); tj <- list(); tv <- list()
  add <- function(rows, cols, vals) {
    ti[[length(ti) + 1L]] <<- rows
    tj[[length(tj) + 1L]] <<- cols
    tv[[length(tv) + 1L]] <<- vals
  }
  # d f_r at `pts` (scaled by coef) with respect to the decision variables;
  # the inputs are those of `upts` and the derivatives rows `idx` of JA
  add_node_jac <- function(rows, pts, r, coef, upts = pts, JA = JJ, idx = pts) {
    for (q in seq_len(nc)) {
      vals <- coef * JA[idx, r, q]
      if (q <= n) cols <- (q - 1L) * Np + pts
      else if (q <= n + m) cols <- Np * n + (q - n - 1L) * Np + upts
      else if (q <= n + m + k) cols <- rep(Np * (n + m) + q - n - m, length(pts))
      else if (ctx$free) {
        cols <- rep(ctx$iT, length(pts))
        if (q == n + m + k + 1L) vals <- vals * ctx$ptau[pts]
      } else next
      add(rows, cols, vals)
    }
  }
  A <- ctx$scheme$A; b <- ctx$scheme$b
  spt <- ctx$spt
  nk <- ctx$node_pt[-N]; nk1 <- ctx$node_pt[-1L]
  h <- ctx$dtau
  cE <- numeric(); keyE <- numeric()
  dyn <- list()                      # (rows, r, stage coefficients) for the Hessian
  row <- 0L
  expl <- list()
  if (isTRUE(ctx$explicit)) {
    # U_k + dt sum_i alpha_i K_i - U_{k+1}, as MTK's add_solve_constraints!
    for (kk in seq_len(N - 1L)) {
      st <- .dopt_explicit_step(ctx, kk, X[nk[kk], ], W[nk[kk], ], P, Tv, level)
      rows <- row + seq_len(n)
      cE <- c(cE, st$val - X[nk1[kk], ])
      keyE <- c(keyE, rep(nk1[kk], n))
      if (level >= 1L) {
        cols <- .dopt_point_cols(ctx, nk[kk])
        add(rep(rows, length(cols)), rep(cols, each = n), as.numeric(st$D))
        add(rows, (seq_len(n) - 1L) * Np + nk1[kk], rep(-1, n))
      }
      expl[[kk]] <- rows
      row <- row + n
    }
  }
  kst <- list()
  Kv <- NULL
  if (isTRUE(ctx$kstage)) {
    # MTK's add_solve_constraints! for an implicit tableau: per step, the
    # stage equations K_i - T f(U_k + dt sum_l A_il K_l, V_k, P, t_k + c_i dt)
    # (one row per state) and U_k + dt sum_i alpha_i K_i - U_{k+1}
    s <- length(b)
    Kv <- array(z[ctx$kbase + seq_len(ctx$nK * s * n)], c(n, s, ctx$nK))
    for (kk in seq_len(N - 1L)) {
      Kk <- Kv[, , kk, drop = FALSE]
      dim(Kk) <- c(n, s)
      st <- .dopt_kstage_step(ctx, kk, X[nk[kk], ], W[nk[kk], ], P, Tv, Kk, level)
      for (i in seq_len(s)) {
        rows <- row + seq_len(n)
        cE <- c(cE, Kk[, i] - st$val[, i])
        keyE <- c(keyE, rep(nk1[kk], n))
        if (level >= 1L) {
          D <- -st$D[[i]]
          own <- cbind(seq_len(n), st$np + (i - 1L) * n + seq_len(n))
          D[own] <- D[own] + 1
          add(rep(rows, length(st$cols)), rep(st$cols, each = n), as.numeric(D))
        }
        kst[[length(kst) + 1L]] <- list(rows = rows, kk = kk, i = i)
        row <- row + n
      }
      rows <- row + seq_len(n)
      cE <- c(cE, X[nk[kk], ] + .dopt_support_dt(ctx) * as.numeric(Kk %*% b) - X[nk1[kk], ])
      keyE <- c(keyE, rep(nk1[kk], n))
      if (level >= 1L) {
        add(rows, (seq_len(n) - 1L) * Np + nk[kk], rep(1, n))
        add(rows, (seq_len(n) - 1L) * Np + nk1[kk], rep(-1, n))
        for (i in seq_len(s))
          add(rows, ctx$kbase + ((kk - 1L) * s + i - 1L) * n + seq_len(n), rep(.dopt_support_dt(ctx) * b[i], n))
      }
      row <- row + n
    }
  }
  dvr <- list()
  if (isTRUE(ctx$dvar)) {
    # InfiniteOpt: (D x - T f) / nominal at every support, then the
    # relations of the derivative method on every element
    Tsc <- if (ctx$free) Tv else 1
    dcol <- function(q, pts) ctx$kbase + (q - 1L) * Np + pts
    for (q in seq_along(ctx$drows)) {
      r <- ctx$drows[q]; j <- ctx$diff_state[r]; sc <- ctx$scale[j]
      Dq <- z[dcol(q, seq_len(Np))]
      rows <- row + seq_len(Np)
      cE <- c(cE, (Dq - Tsc * FF[, r]) / sc)
      keyE <- c(keyE, seq_len(Np))
      if (level >= 1L) {
        add(rows, dcol(q, seq_len(Np)), rep(1 / sc, Np))
        add_node_jac(rows, seq_len(Np), r, -Tsc / sc)
        if (ctx$free) add(rows, rep(ctx$iT, Np), -FF[, r] / sc)
      }
      dvr[[length(dvr) + 1L]] <- list(rows = rows, r = r, sc = sc)
      row <- row + Np
    }
    ep <- ctx$dv_pts; C <- ctx$dv_C
    for (kk in seq_len(N - 1L)) {
      qp <- spt[kk, ep]
      for (q in seq_along(ctx$drows)) {
        j <- ctx$diff_state[ctx$drows[q]]
        Dq <- z[dcol(q, qp)]
        for (e in seq_along(ep)) {
          row <- row + 1L
          cE <- c(cE, sum(C[, e] * Dq) - X[qp[e], j] + X[nk[kk], j])
          keyE <- c(keyE, qp[e])
          if (level >= 1L) {
            add(rep(row, length(qp)), dcol(q, qp), C[, e])
            add(c(row, row), (j - 1L) * Np + c(qp[e], nk[kk]), c(-1, 1))
          }
        }
      }
    }
  }
  for (r in if (isTRUE(ctx$explicit) || isTRUE(ctx$kstage) || isTRUE(ctx$dvar)) integer()
       else which(!is.na(ctx$diff_state))) {
    j <- ctx$diff_state[r]
    sc <- ctx$scale[j]
    eqs <- c(lapply(which(ctx$stage_interior), function(i) list(target = spt[, i], coef = A[i, ])),
             list(list(target = nk1, coef = b)))
    for (e in eqs) {
      rows <- row + seq_len(N - 1L)
      fsum <- numeric(N - 1L)
      held <- hold & ctx$stage_end
      for (l in which(e$coef != 0))
        fsum <- fsum + e$coef[l] * (if (held[l]) FH[, r] else FF[spt[, l], r])
      # stage form scaled by the step and the nominal value, as the
      # derivative constraints of the InfiniteOpt backend
      cE <- c(cE, ((X[e$target, j] - X[nk, j]) / h - Tv * fsum) / sc)
      keyE <- c(keyE, nk1)
      if (level >= 1L) {
        add(rows, (j - 1L) * Np + e$target, rep(1 / (h * sc), N - 1L))
        add(rows, (j - 1L) * Np + nk, rep(-1 / (h * sc), N - 1L))
        for (l in which(e$coef != 0)) {
          if (held[l]) add_node_jac(rows, hx, r, -Tv * e$coef[l] / sc, upts = hu, JA = JH, idx = seq_len(N - 1L))
          else add_node_jac(rows, spt[, l], r, -Tv * e$coef[l] / sc)
        }
        if (ctx$free) add(rows, rep(ctx$iT, N - 1L), -fsum / sc)
      }
      dyn[[length(dyn) + 1L]] <- list(rows = rows, r = r, coef = e$coef / sc)
      row <- row + N - 1L
    }
  }
  alg <- list()
  for (r in which(is.na(ctx$diff_state))) {
    rows <- row + seq_len(Np)
    cE <- c(cE, FF[, r])
    keyE <- c(keyE, seq_len(Np))
    if (level >= 1L) add_node_jac(rows, seq_len(Np), r, 1)
    alg[[length(alg) + 1L]] <- list(rows = rows, r = r)
    row <- row + Np
  }
  for (tie in ctx$ties) {
    for (q in seq_len(m)) {
      row <- row + 1L
      base <- Np * n + (q - 1L) * Np
      cE <- c(cE, W[tie$pt, q] - sum(tie$w * W[tie$sp, q]))
      keyE <- c(keyE, max(tie$pt, tie$sp))
      if (level >= 1L) add(c(row, rep(row, length(tie$sp))), base + c(tie$pt, tie$sp), c(1, -tie$w))
    }
  }
  for (j in if (isTRUE(ctx$u0_fixed)) integer() else ctx$u0_idx) {
    row <- row + 1L
    cE <- c(cE, X[ctx$node_pt[1L], j] - ctx$u0[j])
    keyE <- c(keyE, 1)
    if (level >= 1L) add(row, (j - 1L) * Np + ctx$node_pt[1L], 1)
  }
  exE <- list()
  for (it in ctx$eqc) {
    ex <- .dopt_expr(it, ctx, X, W, P, Tv, tt, level)
    ex$row0 <- row
    if (level >= 1L) add(row + ex$J$i, ex$J$j, ex$J$v)
    cE <- c(cE, ex$v); keyE <- c(keyE, ex$key)
    row <- row + length(ex$v)
    exE[[length(exE) + 1L]] <- ex
  }
  nE <- row
  JE <- list(i = unlist(ti), j = unlist(tj), v = unlist(tv))
  cI <- numeric(); keyI <- numeric()
  JI <- list(i = integer(), j = integer(), v = numeric())
  exI <- list()
  rowi <- 0L
  for (it in ctx$inc) {
    ex <- .dopt_expr(it, ctx, X, W, P, Tv, tt, level)
    ex$row0 <- rowi
    if (level >= 1L) JI <- list(i = c(JI$i, rowi + ex$J$i), j = c(JI$j, ex$J$j), v = c(JI$v, ex$J$v))
    cI <- c(cI, ex$v); keyI <- c(keyI, ex$key)
    rowi <- rowi + length(ex$v)
    exI[[length(exI) + 1L]] <- ex
  }
  cf <- .dopt_expr(ctx$cost, ctx, X, W, P, Tv, tt, level)
  out <- list(f = cf$v, cE = cE, cI = cI, keyE = keyE, keyI = keyI, X = X, W = W, P = P, T = Tv, t = tt, FF = FF,
              K = Kv)
  if (level < 1L) return(out)
  out$gf <- .nlp_accum(cf$J$j, cf$J$v, nz)
  out$JE <- JE; out$JI <- JI
  if (dense) {
    D <- matrix(0, nE, nz)
    if (length(JE$i)) { s <- rowsum(JE$v, (JE$j - 1) * nE + JE$i); D[as.numeric(rownames(s))] <- s[, 1L] }
    out$JE <- D
    D <- matrix(0, length(cI), nz)
    if (length(JI$i)) { s <- rowsum(JI$v, (JI$j - 1) * length(cI) + JI$i); D[as.numeric(rownames(s))] <- s[, 1L] }
    out$JI <- D
  }
  if (level < 2L) return(out)
  # Hessian of the Lagrangian: the dynamics point by point, with weights
  # alpha (algebraic rows) and beta (collocation rows, times T)
  lamE <- lamE %||% numeric(nE); lamI <- lamI %||% numeric(length(cI))
  al <- matrix(0, Np, n); be <- matrix(0, Np, n); bh <- matrix(0, N - 1L, n)
  for (d in dyn) {
    lr <- lamE[d$rows]
    for (l in which(d$coef != 0)) {
      if (hold && ctx$stage_end[l]) bh[, d$r] <- bh[, d$r] - lr * d$coef[l]
      else be[spt[, l], d$r] <- be[spt[, l], d$r] - lr * d$coef[l]
    }
  }
  for (a in alg) al[, a$r] <- al[, a$r] + lamE[a$rows]
  for (d in dvr) {
    if (ctx$free) be[, d$r] <- be[, d$r] - lamE[d$rows] / d$sc
    else al[, d$r] <- al[, d$r] - lamE[d$rows] / d$sc
  }
  H <- list()
  nmk <- n + m + k
  for (pt in seq_len(Np)) {
    alp <- al[pt, ]; bep <- be[pt, ]
    if (all(alp == 0) && all(bep == 0)) next
    H[[length(H) + 1L]] <- .dopt_point_hess(ctx, pt, X, W, P, Tv, function(x, u, Pp, Tt) {
      tp <- ctx$t0 + ctx$ptau[pt] * Tt
      J <- ctx$fJ(x, u, Pp, tp, Tt)
      wv <- alp + Tt * bep
      g <- as.numeric(crossprod(J[, seq_len(nmk), drop = FALSE], wv))
      if (ctx$free) {
        Fv <- ctx$fF(x, u, Pp, tp, Tt)
        g <- c(g, sum(bep * Fv) + sum(wv * (J[, nmk + 1L] * ctx$ptau[pt] + J[, nmk + 2L])))
      }
      g
    })
  }
  if (hold) for (i in which(rowSums(bh != 0) > 0)) {
    bep <- bh[i, ]
    ptx <- hx[i]
    H[[length(H) + 1L]] <- .dopt_point_hess(ctx, ptx, X, W, P, Tv, function(x, u, Pp, Tt) {
      tp <- ctx$t0 + ctx$ptau[ptx] * Tt
      J <- ctx$fJ(x, u, Pp, tp, Tt)
      g <- Tt * as.numeric(crossprod(J[, seq_len(nmk), drop = FALSE], bep))
      if (ctx$free) {
        Fv <- ctx$fF(x, u, Pp, tp, Tt)
        g <- c(g, sum(bep * Fv) + Tt * sum(bep * (J[, nmk + 1L] * ctx$ptau[ptx] + J[, nmk + 2L])))
      }
      g
    }, ptu = hu[i])
  }
  for (kk in seq_along(expl)) {
    lr <- lamE[expl[[kk]]]
    if (all(lr == 0)) next
    H[[length(H) + 1L]] <- .dopt_point_hess(ctx, nk[kk], X, W, P, Tv, function(x, u, Pp, Tt) {
      as.numeric(crossprod(.dopt_explicit_step(ctx, kk, x, u, Pp, Tt, 1L)$D, lr))
    })
  }
  for (e in kst) {
    lr <- lamE[e$rows]
    if (all(lr == 0)) next
    Kk <- Kv[, , e$kk, drop = FALSE]
    dim(Kk) <- dim(Kv)[1:2]
    H[[length(H) + 1L]] <- .dopt_kstage_hess(ctx, e$kk, e$i, lr, X, W, P, Tv, Kk)
  }
  H <- c(H, .dopt_expr_hess(ctx$cost, cf, obj, ctx, X, W, P, Tv, tt))
  for (q in seq_along(ctx$eqc)) {
    ex <- exE[[q]]
    H <- c(H, .dopt_expr_hess(ctx$eqc[[q]], ex, lamE[ex$row0 + seq_along(ex$v)], ctx, X, W, P, Tv, tt))
  }
  for (q in seq_along(ctx$inc)) {
    ex <- exI[[q]]
    H <- c(H, .dopt_expr_hess(ctx$inc[[q]], ex, lamI[ex$row0 + seq_along(ex$v)], ctx, X, W, P, Tv, tt))
  }
  out$H <- list(i = unlist(lapply(H, `[[`, "i")), j = unlist(lapply(H, `[[`, "j")),
                v = unlist(lapply(H, `[[`, "v")))
  out
}

# One step of an explicit tableau from node kk: the value U_k + dt sum_i
# alpha_i K_i with K_i = T f(U_k + dt sum_{j<i} A_ij K_j, V_k, P,
# t_k + c_i dt) and, for level >= 1, its Jacobian over the point variables
# (state, input, tuned parameters and, for a free final time, T).
.dopt_explicit_step <- function(ctx, kk, x, u, P, Tv, level = 1L) {
  n <- ctx$n; m <- ctx$m; k <- ctx$k
  A <- ctx$scheme$A; b <- ctx$scheme$b; cc <- ctx$scheme$c
  h <- ctx$dtau
  rj <- match(seq_len(n), ctx$diff_state)
  nv <- n + m + k + as.integer(ctx$free)
  tau0 <- ctx$tau[kk]
  K <- vector("list", length(b)); dK <- vector("list", length(b))
  val <- x
  D <- if (level >= 1L) cbind(diag(n), matrix(0, n, nv - n)) else NULL
  for (i in seq_along(b)) {
    Y <- x
    dY <- if (level >= 1L) cbind(diag(n), matrix(0, n, nv - n)) else NULL
    for (j in seq_len(i - 1L)) if (A[i, j] != 0) {
      Y <- Y + h * A[i, j] * K[[j]]
      if (level >= 1L) dY <- dY + h * A[i, j] * dK[[j]]
    }
    ti <- tau0 + cc[i] * h
    tp <- ctx$t0 + ti * Tv
    Fv <- ctx$fF(Y, u, P, tp, Tv)[rj]
    K[[i]] <- Tv * Fv
    if (level >= 1L) {
      J <- ctx$fJ(Y, u, P, tp, Tv)[rj, , drop = FALSE]
      dF <- J[, seq_len(n), drop = FALSE] %*% dY
      if (m + k) dF[, n + seq_len(m + k)] <- dF[, n + seq_len(m + k)] + J[, n + seq_len(m + k)]
      dK[[i]] <- Tv * dF
      if (ctx$free) dK[[i]][, nv] <- dK[[i]][, nv] + Fv + Tv * (J[, n + m + k + 1L] * ti + J[, n + m + k + 2L])
    }
    if (b[i] != 0) {
      val <- val + h * b[i] * K[[i]]
      if (level >= 1L) D <- D + h * b[i] * dK[[i]]
    }
  }
  list(val = val, D = D)
}

# The stages of an implicit tableau with derivative stage variables Kk (n x
# s) at step kk: val[, i] = T f(U_k + dt sum_l A_il K_l, V_k, P, t_k + c_i dt)
# and, for level >= 1, its Jacobian D[[i]] over the columns `cols` (the
# point variables of node kk, then the np + seq_len(s * n) stage variables
# of the step).
.dopt_kstage_step <- function(ctx, kk, x, u, P, Tv, Kk, level = 1L) {
  n <- ctx$n; m <- ctx$m; k <- ctx$k
  A <- ctx$scheme$A; cc <- ctx$scheme$c
  s <- length(cc)
  h <- ctx$dtau
  # K is in the units of the supports (the times, or [0, 1] scaled by the
  # free duration): K = Ts f with dt the support spacing
  Tsc <- if (ctx$free) Tv else 1
  hs <- .dopt_support_dt(ctx)
  rj <- match(seq_len(n), ctx$diff_state)
  cols <- c(.dopt_point_cols(ctx, ctx$node_pt[kk]), ctx$kbase + (kk - 1L) * s * n + seq_len(s * n))
  np <- n + m + k + as.integer(ctx$free)
  val <- matrix(0, n, s)
  D <- vector("list", s)
  for (i in seq_len(s)) {
    Y <- x + as.numeric(Kk %*% A[i, ]) * hs
    ti <- ctx$tau[kk] + cc[i] * h
    tp <- ctx$t0 + ti * Tv
    Fv <- ctx$fF(Y, u, P, tp, Tv)[rj]
    val[, i] <- Tsc * Fv
    if (level >= 1L) {
      J <- ctx$fJ(Y, u, P, tp, Tv)[rj, , drop = FALSE]
      Jx <- Tsc * J[, seq_len(n), drop = FALSE]
      Di <- matrix(0, n, np + s * n)
      Di[, seq_len(n)] <- Jx
      if (m + k) Di[, n + seq_len(m + k)] <- Tsc * J[, n + seq_len(m + k)]
      if (ctx$free) Di[, np] <- Fv + Tv * (J[, n + m + k + 1L] * ti + J[, n + m + k + 2L])
      for (l in which(A[i, ] != 0)) Di[, np + (l - 1L) * n + seq_len(n)] <- Jx * (hs * A[i, l])
      D[[i]] <- Di
    }
  }
  list(val = val, D = D, cols = cols, np = np)
}

# The spacing of the supports of the time parameter (InfiniteOpt and CasADi
# work in the problem's time, or on [0, 1] for a free final time).
.dopt_support_dt <- function(ctx) ctx$dtau * (if (ctx$free) 1 else ctx$Tfix)

# Hessian triplets of -sum_j lr_j K-stage value (stage i of step kk) by
# forward differences of its gradient.
.dopt_kstage_hess <- function(ctx, kk, i, lr, X, W, P, Tv, Kk) {
  n <- ctx$n; m <- ctx$m; k <- ctx$k
  s <- ncol(Kk)
  pt <- ctx$node_pt[kk]
  v <- c(X[pt, ], W[pt, ], P, if (ctx$free) Tv, as.numeric(Kk))
  np <- n + m + k + as.integer(ctx$free)
  gfun <- function(v) {
    st <- .dopt_kstage_step(ctx, kk, v[seq_len(n)], v[n + seq_len(m)], v[n + m + seq_len(k)],
                            if (ctx$free) v[np] else Tv, matrix(v[np + seq_len(s * n)], n, s), 1L)
    -as.numeric(crossprod(st$D[[i]], lr))
  }
  H <- .dopt_fd_hess(gfun, v)
  cols <- c(.dopt_point_cols(ctx, pt), ctx$kbase + (kk - 1L) * s * n + seq_len(s * n))
  nzh <- which(H != 0, arr.ind = TRUE)
  list(i = cols[nzh[, 1L]], j = cols[nzh[, 2L]], v = H[nzh])
}

# ---------------------------------------------------------------------------
# Solving

# The interior point method of R/dynamic-opt-nlp.R on the collocation
# problem.
.dopt_solve_ip <- function(ctx, tol, maxiter, verbose) {
  keys <- .dopt_eval(ctx, ctx$z0, level = 0L)
  evalf <- function(z, level, lamE = NULL, lamI = NULL, obj = 1) .dopt_eval(ctx, z, level, dense = FALSE, lamE, lamI, obj)
  r <- .nlp_interior_point(evalf, ctx$z0, ctx$lower, ctx$upper, ctx$var_key, keys$keyE, keys$keyI, tol = tol,
                           maxiter = maxiter, verbose = verbose)
  # successful_solve of ModelingToolkit: Ipopt's Solve_Succeeded,
  # Solved_To_Acceptable_Level and Feasible_Point_Found (LOCALLY_SOLVED and
  # ALMOST_LOCALLY_SOLVED with a feasible point); otherwise ConvergenceFailure
  status <- if (r$status %in% c("Success", "Acceptable", "FeasiblePointFound")) "Success" else "ConvergenceFailure"
  list(z = r$z, v = .dopt_eval(ctx, r$z, level = 0L), lambda = r$lamE, nu = r$lamI, violation = r$violation,
       stationarity = NA_real_, iterations = r$iterations, status = status)
}

# The augmented Lagrangian of the optimization module (AugLag() in
# optimization-solvers.R), after the same feasibility phase.
.dopt_shared_auglag <- function(evalf, ctx, tol, maxiter, inner_maxit) {
  cache <- new.env(parent = emptyenv())
  get_v <- function(z) {
    if (!identical(cache$z, z)) {
      cache$z <- z
      cache$v <- evalf(z)
    }
    cache$v
  }
  z0 <- .dopt_restore(get_v, ctx$z0, get_v(ctx$z0), ctx$lower, ctx$upper, tol, inner_maxit)
  v0 <- get_v(z0)
  nE <- length(v0$cE); nI <- length(v0$cI)
  r <- .opt_auglag(f = function(z) get_v(z)$f, g = function(z) get_v(z)$gf,
                   cons = function(z) { v <- get_v(z); c(v$cE, v$cI) },
                   cons_j = function(z) { v <- get_v(z); rbind(v$JE, v$JI) },
                   lcons = c(rep(0, nE), rep(-Inf, nI)), ucons = rep(0, nE + nI), x0 = z0,
                   lb = ctx$lower, ub = ctx$upper, alg = AugLag(ctol = tol, max_outer = maxiter),
                   maxiters = inner_maxit)
  v <- get_v(r$u)
  list(z = r$u, v = v, lambda = r$lambda, nu = r$mu, violation = r$viol, stationarity = NA_real_,
       iterations = r$outer, status = r$retcode)
}

# Feasibility phase: from an infeasible start (such as a too short final
# time), minimize the squared constraint violation first, so the augmented
# Lagrangian does not settle in an infeasible local minimum.
.dopt_restore <- function(get_v, z, v, lower, upper, tol, inner_maxit) {
  if (max(c(0, abs(v$cE), pmax(v$cI, 0))) <= tol) return(z)
  Ff <- function(z) {
    v <- get_v(z)
    (sum(v$cE^2) + sum(pmax(v$cI, 0)^2)) / 2
  }
  Fg <- function(z) {
    v <- get_v(z)
    as.numeric(crossprod(v$JE, v$cE) + crossprod(v$JI, pmax(v$cI, 0)))
  }
  stats::optim(z, Ff, Fg, method = "L-BFGS-B", lower = lower, upper = upper,
               control = list(maxit = inner_maxit, factr = 10, pgtol = 0))$par
}

# Augmented Lagrangian method for min f(z) s.t. cE(z) = 0, cI(z) <= 0 and
# lower <= z <= upper, with L-BFGS-B for the subproblems. Stops when the
# constraints hold to `tol` and the projected gradient of the Lagrangian is
# small, or after `maxiter` outer iterations.
.dopt_auglag <- function(evalf, z, lower, upper, tol = 1e-6, maxiter = 30L, inner_maxit = 2000L,
                         verbose = FALSE) {
  cache <- new.env(parent = emptyenv())
  get_v <- function(z) {
    if (!identical(cache$z, z)) {
      cache$z <- z
      cache$v <- evalf(z)
    }
    cache$v
  }
  v <- get_v(z)
  z <- .dopt_restore(get_v, z, v, lower, upper, tol, inner_maxit)
  v <- get_v(z)
  lam <- numeric(length(v$cE)); nu <- numeric(length(v$cI))
  mu <- 10
  viol_prev <- Inf
  f_prev <- Inf
  status <- "MaxIters"
  iters <- 0L
  kkt <- Inf
  for (it in seq_len(maxiter)) {
    iters <- it
    Lf <- function(z) {
      v <- get_v(z)
      s <- pmax(0, nu + mu * v$cI)
      v$f + sum(lam * v$cE) + mu / 2 * sum(v$cE^2) + (sum(s^2) - sum(nu^2)) / (2 * mu)
    }
    Lg <- function(z) {
      v <- get_v(z)
      s <- pmax(0, nu + mu * v$cI)
      as.numeric(v$gf + crossprod(v$JE, lam + mu * v$cE) + crossprod(v$JI, s))
    }
    opt <- stats::optim(z, Lf, Lg, method = "L-BFGS-B", lower = lower, upper = upper,
                        control = list(maxit = inner_maxit, factr = 1e5, pgtol = 0))
    z <- opt$par
    v <- get_v(z)
    lam <- lam + mu * v$cE
    nu <- pmax(0, nu + mu * v$cI)
    viol <- max(c(0, abs(v$cE), pmax(v$cI, 0)))
    g <- as.numeric(v$gf + crossprod(v$JE, lam) + crossprod(v$JI, nu))
    kkt <- max(c(0, abs(z - pmin(pmax(z - g, lower), upper))))
    if (verbose) message(sprintf("iteration %d: objective %.8g, violation %.3g, stationarity %.3g, mu %.3g, evaluations %d",
                                 it, v$f, viol, kkt, mu, opt$counts[[1L]]))
    fchange <- abs(v$f - f_prev) <= sqrt(tol) * (1 + abs(v$f))
    if (viol <= tol && (kkt <= 1e-3 * (1 + max(abs(v$gf))) || fchange)) {
      status <- "Success"
      break
    }
    if (viol > tol && viol > 0.25 * viol_prev) mu <- min(mu * 10, 1e10)
    viol_prev <- viol
    f_prev <- v$f
  }
  if (status != "Success" && viol > sqrt(tol)) status <- "Infeasible"
  list(z = z, v = v, lambda = lam, nu = nu, violation = viol, stationarity = kkt, iterations = iters,
       status = status)
}


# The solution between the nodes: the collocation polynomial of each
# interval, x(t_k + theta h) = x_k + h sum_l beta_l(theta) f_l, for the
# differential states, Lagrange interpolation of the values at the stage
# points of the interval for the algebraic states and the inputs, and cubic
# Hermite interpolation with the node derivatives (linear interpolation for
# the inputs) for tableaus that are not collocation methods.
.dopt_steps <- function(ctx, v, values, deriv = TRUE) {
  N <- ctx$N; s <- length(ctx$scheme$b)
  beta <- ctx$scheme$beta
  Tv <- v$T
  h <- ctx$dtau * Tv
  diffj <- ctx$diff_state[!is.na(ctx$diff_state)]
  rows_of <- which(!is.na(ctx$diff_state))
  steps <- vector("list", N - 1L)
  for (k in seq_len(N - 1L)) {
    a <- ctx$node_pt[k]; b <- ctx$node_pt[k + 1L]
    pts <- if (!deriv && isTRUE(ctx$hold)) a else unique(ctx$spt[k, ])
    th <- (ctx$ptau[pts] - ctx$ptau[a]) / ctx$dtau
    kk <- list(V = values[pts, , drop = FALSE], th = th, beta = beta, diff = integer(), K = NULL)
    if (deriv && length(diffj)) {
      kk$diff <- diffj
      if (!is.null(beta) && isTRUE(ctx$kstage)) {
        # the derivative stage variables (T f at the stages)
        Kk <- v$K[, , k, drop = FALSE]
        dim(Kk) <- dim(v$K)[1:2]
        kk$K <- t(Kk)[, diffj, drop = FALSE] / (if (ctx$free) Tv else 1)
      } else if (!is.null(beta)) {
        kk$K <- v$FF[ctx$spt[k, ], rows_of, drop = FALSE]
      } else {
        kk$f0 <- v$FF[a, rows_of]; kk$f1 <- v$FF[b, rows_of]
      }
    }
    steps[[k]] <- list(t0 = v$t[a], dt = h, y0 = values[a, ], y1 = values[b, ], k = kk,
                       interp = .dopt_interp)
  }
  steps
}

.dopt_interp <- function(theta, dt, y0, y1, k) {
  th <- k$th
  if (length(unique(round(th, 12))) == length(th)) {
    # Lagrange interpolation through the values at the points of the interval
    lw <- vapply(seq_along(th), function(i) prod((theta - th[-i]) / (th[i] - th[-i])), 0)
    out <- as.numeric(crossprod(k$V, lw))
  } else {
    out <- y0 + theta * (y1 - y0)
  }
  if (length(k$diff)) {
    j <- k$diff
    if (!is.null(k$K)) {
      bt <- as.numeric(k$beta %*% theta^seq_len(ncol(k$beta)))
      out[j] <- y0[j] + dt * as.numeric(crossprod(k$K, bt))
    } else {
      h00 <- 2 * theta^3 - 3 * theta^2 + 1; h10 <- theta^3 - 2 * theta^2 + theta
      h01 <- -2 * theta^3 + 3 * theta^2; h11 <- theta^3 - theta^2
      out[j] <- h00 * y0[j] + h10 * dt * k$f0 + h01 * y1[j] + h11 * dt * k$f1
    }
  }
  names(out) <- names(y0)
  out
}

#' @rdname DynamicOptProblem
#' @param a A `DynamicOptProblem`.
#' @param b The collocation scheme: a collocation object such as
#'   [JuMPCollocation()], a name such as `"hermite_simpson"`, or a tableau;
#'   by default the scheme of the problem's backend, and the trapezoidal
#'   rule for a plain `DynamicOptProblem`.
#' @param tol Tolerance of the optimality conditions (the constraint
#'   violation for the augmented Lagrangian methods).
#' @param maxiter Maximum number of iterations.
#' @param inner_maxit Maximum number of L-BFGS-B iterations per subproblem
#'   of the augmented Lagrangian methods.
#' @param verbose Print the progress of the optimizer.
#' @param optimizer `"interior_point"` is a primal dual interior point
#'   method in the manner of Ipopt; `"optimization"` uses the augmented
#'   Lagrangian of the optimization module (`AugLag()`) and `"native"` an
#'   augmented Lagrangian of this file. `"auto"` (the default) is the
#'   interior point method, with Ipopt's restoration phase; a solve that
#'   does not converge has the retcode `ConvergenceFailure`, as with the
#'   JuMP backend of ModelingToolkit.
#' @export
solve.DynamicOptProblem <- function(a, b, tol = NULL, maxiter = NULL, inner_maxit = 2000L, verbose = FALSE,
                                    optimizer = c("auto", "interior_point", "native", "optimization"), ...) {
  optimizer <- match.arg(optimizer)
  scheme <- .dopt_pick_scheme(a, if (missing(b)) NULL else b)
  ctx <- .dopt_grid(.subset2(a, "ctx"), scheme)
  sys <- .subset2(a, "sys")
  ip <- optimizer %in% c("auto", "interior_point")
  tol <- tol %||% if (ip) 1e-8 else 1e-6
  maxiter <- maxiter %||% if (ip) 3000L else 30L
  evalf <- function(z) .dopt_eval(ctx, z)
  res <- if (ip) .dopt_solve_ip(ctx, tol, maxiter, verbose)
  else if (optimizer == "optimization") .dopt_shared_auglag(evalf, ctx, tol, maxiter, inner_maxit)
  else .dopt_auglag(evalf, ctx$z0, ctx$lower, ctx$upper, tol = tol, maxiter = maxiter,
                    inner_maxit = inner_maxit, verbose = verbose)
  v <- .dopt_eval(ctx, res$z, level = 0L)
  # ModelingToolkit's solve builds the solutions with build_solution (the
  # retcode Default) and sets ConvergenceFailure when successful_solve fails
  retcode <- if (identical(res$status, "Success")) "Default" else "ConvergenceFailure"
  nodes <- ctx$node_pt
  times <- v$t[nodes]
  pv <- ctx$pv
  pv[ctx$tun_idx] <- v$P
  if (ctx$free) pv[ctx$tf_idx] <- ctx$t0 + v$T
  pv[ctx$ctrl_idx] <- if (ctx$m) v$W[ctx$Npts, ] else numeric()
  pv[is.na(pv)] <- 0
  names(pv) <- vapply(ctx$pp, format, "")
  Xall <- v$X
  colnames(Xall) <- vapply(ctx$states, format, "")
  X <- Xall[nodes, , drop = FALSE]
  prob <- .new_problem("ODEProblem", f = ODEFunction(sys), u0 = X[1L, ], tspan = range(times), p = pv,
                       sys = sys, kwargs = list())
  alg <- list(name = paste0("DirectCollocation (", scheme$name, ")"))
  sol <- .new_ode_solution(t = times, u = X, prob = prob, alg = alg, retcode = retcode, p = pv,
                           steps = .dopt_steps(ctx, v, Xall), dense = TRUE)
  if (ctx$m) {
    Pm <- matrix(pv, ctx$N, length(pv), byrow = TRUE)
    Pm[, ctx$ctrl_idx] <- v$W[nodes, , drop = FALSE]
    env <- .sol_env(sol)
    env$param_history <- list(t = times, p = Pm)
  }
  input_sol <- NULL
  if (ctx$m) {
    Wall <- v$W
    colnames(Wall) <- vapply(ctx$ctrls, format, "")
    iprob <- .new_problem("ODEProblem", f = NULL, u0 = Wall[1L, ], tspan = range(times), p = numeric(),
                          sys = NULL, kwargs = list())
    input_sol <- .new_ode_solution(t = times, u = Wall[nodes, , drop = FALSE], prob = iprob, alg = alg,
                                   retcode = retcode, steps = .dopt_steps(ctx, v, Wall, deriv = FALSE),
                                   dense = TRUE)
  }
  model <- list(objective = v$f, z = res$z, constraint_violation = res$violation,
                stationarity = res$stationarity, iterations = res$iterations,
                multipliers = list(equality = res$lambda, inequality = res$nu),
                parameters = stats::setNames(v$P, vapply(ctx$tun, format, "")),
                final_time = ctx$t0 + v$T, method = scheme$name, backend = .subset2(a, "backend"),
                status = res$status)
  DynamicOptSolution(model, sol, input_sol)
}

#' Solutions of dynamic optimization problems
#'
#' `DynamicOptSolution(model, sol, input_sol)` holds the result of
#' [solve()] on a [DynamicOptProblem()]: `model` describes the solved
#' nonlinear program (objective value, constraint violation, iterations,
#' optimized parameters and final time), `sol` is an `ODESolution` of the
#' states on the collocation grid (index it with variables, `sol$sol[x]`; its
#' `retcode` reports convergence), and `input_sol` an `ODESolution` of the
#' inputs (or `NULL` without inputs).
#'
#' @param model A list describing the optimization result.
#' @param sol An `ODESolution` of the states.
#' @param input_sol An `ODESolution` of the inputs, or `NULL`.
#' @return An object of class `DynamicOptSolution`.
#' @export
DynamicOptSolution <- function(model, sol, input_sol = NULL) {
  structure(list(model = model, sol = sol, input_sol = input_sol), class = "DynamicOptSolution")
}

#' @export
print.DynamicOptSolution <- function(x, ...) {
  cat("retcode:", x$sol$retcode, "\n\n")
  cat("Optimal control solution (direct collocation, ", x$model$method, ")\n", sep = "")
  cat("Objective: ", format(x$model$objective), "\n", sep = "")
  cat("Constraint violation: ", format(x$model$constraint_violation), "\n", sep = "")
  cat("Iterations: ", x$model$iterations, "\n", sep = "")
  if (length(x$model$parameters)) {
    cat("Tuned parameters:\n")
    print(x$model$parameters)
  }
  cat("\nQuery the trajectory of the states with sol$sol and of the inputs with sol$input_sol.\n")
  invisible(x)
}
