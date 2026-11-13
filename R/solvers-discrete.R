# Implicit discrete problems, semilinear ODE problems and generate_W, after:
#   ModelingToolkitBase/src/problems/implicitdiscreteproblem.jl
#     ImplicitDiscreteFunction (generate_rhs with implicit_dae and
#     override_discrete: residual rhs - lhs, algebraic equations shifted one
#     step forward), ImplicitDiscreteProblem
#   OrdinaryDiffEq.jl lib/SimpleImplicitDiscreteSolve and
#     lib/ImplicitDiscreteSolve/src/solve.jl  IDSolve: Newton solve of
#     f(u_next, u, p, t + dt) = 0 per step, the initial state made
#     consistent by one solve at t0 (_initialize_dae!)
#   ModelingToolkit.jl src/problems/semilinearodeproblem.jl and
#     src/systems/codegen.jl  calculate_semiquadratic_form,
#     SemilinearODEFunction, SemilinearODEProblem (a SplitODEProblem)
#   ModelingToolkitBase/src/systems/codegen.jl  generate_W
# DiscreteProblem, DiscreteFunction and FunctionMap are in
# R/clocks-discrete.R.

# ---------------------------------------------------------------------------
# Implicit discrete problems

#' Implicit discrete functions and problems
#'
#' `ImplicitDiscreteProblem(sys, op, tspan)` builds the implicit map
#' `f(u_next, u, p, t) = 0` of a compiled discrete system with algebraic
#' equations, as ModelingToolkit's `ImplicitDiscreteProblem`: equations
#' `Shift(x, 1) ~ rhs` give `rhs - x_next` and algebraic equations `0 ~ g`
#' are evaluated at the next step. `ImplicitDiscreteFunction(sys)` builds
#' `f`. `ImplicitDiscreteProblem(f, u0, tspan, p)` takes an R function
#' `f(u_next, u, p, t)`. Solve with `IDSolve()`, which steps by `dt`
#' (default 1) and solves each step with Newton's method; the initial state
#' is made consistent by one solve at the initial time.
#'
#' @param sys A compiled discrete `System`, or an R function.
#' @param op Values of (past) unknowns and parameters.
#' @param tspan Time span.
#' @param u0,p Initial state and parameters of raw problems.
#' @param ... Solver options stored with the problem.
#' @return An `ImplicitDiscreteProblem`, `ImplicitDiscreteFunction`, or
#'   algorithm object.
#' @examples
#' # u_next solves u_next + u_next^3 = u
#' f <- function(un, u, p, t) un + un^3 - u
#' prob <- ImplicitDiscreteProblem(f, 2, c(0, 3))
#' sol <- solve(prob, IDSolve())
#' sol$u
#' @export
ImplicitDiscreteFunction <- function(sys, ...) {
  if (is.function(sys))
    return(structure(list(f = sys, sys = NULL), class = c("ImplicitDiscreteFunction", "SciMLFunction")))
  .need_complete(sys, "ImplicitDiscreteFunction")
  iv <- .f(sys, "iv")
  if (is.null(iv)) stop("an ImplicitDiscreteFunction needs a time dependent system", call. = FALSE)
  us <- unknowns(sys)
  nexts <- lapply(us, .make_shift, steps = 1L)
  obs <- .codegen_observed(sys)
  orule <- if (length(obs)) .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs")) else NULL
  srule <- .rule(us, nexts)
  res <- lapply(equations(sys), function(e) {
    if (.is_zero(e$lhs)) {
      r <- if (is.null(orule)) e$rhs else .subst_fixpoint(e$rhs, orule)
      .subst(r, srule)
    } else e$rhs - e$lhs
  })
  args <- c(list(list(name = "u_next", vars = nexts)), .codegen_args(sys))
  f <- .build_closure(res, args, observed = obs, out = "vector")
  structure(list(f = f, sys = sys, observed = function(expr) build_explicit_observed_function(sys, expr)),
            class = c("ImplicitDiscreteFunction", "SciMLFunction"))
}

#' @rdname ImplicitDiscreteFunction
#' @export
ImplicitDiscreteProblem <- function(sys, op = list(), tspan = NULL, ..., u0 = NULL, p = NULL) {
  dots <- list(...)
  if (is.function(sys) || inherits(sys, "ImplicitDiscreteFunction")) {
    f <- if (inherits(sys, "ImplicitDiscreteFunction")) sys else ImplicitDiscreteFunction(sys)
    if (is.null(u0)) u0 <- op
    # ImplicitDiscreteProblem(f, u0, tspan, p) as in SciMLBase: an unnamed
    # fourth argument is p
    if (is.null(p) && length(dots) && (is.null(names(dots)) || !nzchar(names(dots)[1L]))) {
      p <- dots[[1L]]
      dots <- dots[-1L]
    }
    return(.new_problem("ImplicitDiscreteProblem", f = f, u0 = as.numeric(unlist(u0)),
                        tspan = .as_tspan(tspan), p = if (is.null(p)) numeric() else p, sys = NULL,
                        kwargs = dots))
  }
  .need_complete(sys, "ImplicitDiscreteProblem")
  tspan <- .as_tspan(tspan, sys)
  fn <- ImplicitDiscreteFunction(sys)
  vals <- .process_op(sys, .discrete_op(sys, op))
  u0 <- vals$u0
  # unknowns without a value start from their guesses (or 0) and are found
  # by the initial solve
  g <- guesses(sys)
  us <- unknowns(sys)
  for (i in which(is.na(u0))) {
    gv <- smap_get(g, us[[i]])
    u0[i] <- if (is.numeric(gv)) gv else 0
  }
  # the initialization of the system: history unknowns x_(t-i) keep their
  # values and the algebraic equations at the initial step give the others
  # (the OverrideInit that ImplicitDiscreteSolve uses for MTK problems)
  alg <- Filter(function(e) .is_zero(e$lhs), equations(sys))
  free <- which(vapply(us, function(u) is.null(.meta_get(u, "unshifted")), TRUE))
  initialized <- FALSE
  if (length(alg) && length(alg) == length(free)) {
    g <- .build_closure(lapply(alg, `[[`, "rhs"), .codegen_args(sys), observed = .codegen_observed(sys),
                        out = "vector")
    F <- function(x) {
      u <- u0
      u[free] <- x
      g(u, vals$p, tspan[1L])
    }
    sol <- .safe_newton(F, u0[free], tol = 1e-12)
    if (!sol$converged) stop("could not solve the algebraic equations at the initial step", call. = FALSE)
    u0[free] <- sol$x
    initialized <- TRUE
  }
  prob <- .new_problem("ImplicitDiscreteProblem", f = fn, u0 = u0, tspan = tspan, p = vals$p, sys = sys,
                       kwargs = .problem_kwargs(sys, dots))
  prob$op <- vals$op
  prob$initialized <- initialized
  prob
}

# .newton_solve (R/problems.R) that reports a singular Jacobian as failure.
.safe_newton <- function(F, x0, ...) {
  tryCatch(.newton_solve(F, x0, ...), error = function(e) list(x = x0, converged = FALSE))
}

#' @rdname ImplicitDiscreteFunction
#' @export
IDSolve <- function() structure(list(name = "IDSolve"), class = c("IDSolve", "DiscreteAlgorithm"))

#' @rdname solve.ODEProblem
#' @export
solve.ImplicitDiscreteProblem <- function(a, b, ...) {
  opts <- utils::modifyList(a$kwargs %||% list(), list(...))
  # SimpleIDSolve (R/solvers-simpleid.R)
  if (!missing(b) && inherits(b, "SimpleIDSolve")) return(.sids_solve(a, opts))
  if (!missing(b) && !is.null(b) && !inherits(b, "IDSolve"))
    stop("an ImplicitDiscreteProblem is solved with IDSolve() or SimpleIDSolve()", call. = FALSE)
  t0 <- a$tspan[1L]; tf <- a$tspan[2L]
  dt <- abs(opts$dt %||% 1)
  tol <- opts$abstol %||% 1e-10
  n <- floor(abs(tf - t0) / dt + 1e-9)
  ts <- t0 + sign(tf - t0) * (0:n) * dt
  f <- a$f$f
  p <- a$p
  nf <- 0L
  step <- function(u, t) {
    F <- function(x) {
      nf <<- nf + 1L
      as.numeric(f(x, u, p, t))
    }
    .safe_newton(F, u, tol = tol)
  }
  out <- function(ts, us, retcode) {
    U <- do.call(rbind, us)
    colnames(U) <- names(a$u0)
    .new_ode_solution(t = ts, u = U, prob = a, alg = IDSolve(), retcode = retcode, stats = list(nf = nf))
  }
  u0 <- as.numeric(a$u0)
  names(u0) <- names(a$u0)
  # the initial state is made consistent by one solve at t0, unless the
  # problem's own initialization did it (_initialize_dae!, OverrideInit)
  if (isTRUE(a$initialized)) {
    u <- u0
  } else {
    ini <- step(u0, t0)
    if (!ini$converged) return(out(t0, list(u0), "InitialFailure"))
    u <- stats::setNames(ini$x, names(u0))
  }
  us <- list(u)
  for (i in seq_along(ts)[-1L]) {
    s <- step(u, ts[i])
    if (!s$converged) {
      warning("the nonlinear solve failed at t = ", ts[i], call. = FALSE)
      return(out(ts[seq_len(i - 1L)], us, "ConvergenceFailure"))
    }
    u <- stats::setNames(s$x, names(u0))
    us[[length(us) + 1L]] <- u
  }
  out(ts, us, "Success")
}

# ---------------------------------------------------------------------------
# Semilinear ODE problems

#' Semiquadratic form of a system
#'
#' `calculate_semiquadratic_form(sys)` splits the right-hand side of each
#' equation of a compiled ODE system into `A x + x' B_i x + C_i`: `A` holds
#' the coefficients of the linear terms, `B[[i]]` the upper triangular
#' coefficients of the quadratic terms of equation `i`, and `C` the
#' remaining terms. Coefficients may contain parameters. Entries are `NULL`
#' when the corresponding part is zero.
#'
#' @param sys A complete ODE `System`.
#' @param ... Unused.
#' @return A list with `A` (a list matrix), `B` (a list of list matrices)
#'   and `C` (a list of expressions).
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); y <- Variable("y", t); k <- Parameter("k")
#' sys <- mtkcompile(System(list(D(x) ~ -k * x + x * y, D(y) ~ sin(x) - y), t, name = "s"))
#' calculate_semiquadratic_form(sys)$A
#' @export
calculate_semiquadratic_form <- function(sys, ...) {
  .need_complete(sys, "calculate_semiquadratic_form")
  rhs <- .ode_rhs_exprs(sys)
  dvs <- unknowns(sys)
  n <- length(dvs)
  ph <- lapply(seq_len(n), function(i) Symbol(paste0(".sq", i)))
  low <- .rule(dvs, ph)
  back <- .rule(ph, dvs)
  phk <- .keys(ph)
  zero <- .S0()
  A <- .list_matrix(length(rhs), n, function(i, j) zero)
  B <- lapply(seq_along(rhs), function(i) .list_matrix(n, n, function(i, j) zero))
  C <- rep(list(zero), length(rhs))
  free_of <- function(e) !length(intersect(.keys(free_symbols(e)), phk))
  for (i in seq_along(rhs)) {
    e <- expand(.subst(rhs[[i]], low))
    terms <- if (.type(e) == "Add") .args(e) else list(e)
    for (term in terms) {
      js <- match(intersect(.keys(free_symbols(term)), phk), phk)
      placed <- FALSE
      if (length(js) == 1L) {
        d1 <- diff(term, ph[[js]])
        if (free_of(d1)) {
          A[[i, js]] <- A[[i, js]] + .subst(d1, back)
          placed <- TRUE
        } else {
          d2 <- diff(d1, ph[[js]])
          if (free_of(d2) && .is_zero(expand(term - d2 / 2 * ph[[js]]^2))) {
            B[[i]][[js, js]] <- B[[i]][[js, js]] + .subst(d2 / 2, back)
            placed <- TRUE
          }
        }
      } else if (length(js) == 2L) {
        js <- sort(js)
        d <- diff(diff(term, ph[[js[1L]]]), ph[[js[2L]]])
        if (free_of(d) && .is_zero(expand(term - d * ph[[js[1L]]] * ph[[js[2L]]]))) {
          B[[i]][[js[1L], js[2L]]] <- B[[i]][[js[1L], js[2L]]] + .subst(d, back)
          placed <- TRUE
        }
      }
      if (!placed) C[[i]] <- C[[i]] + .subst(term, back)
    }
  }
  allzero <- function(m) all(vapply(m, .is_zero, TRUE))
  B <- lapply(B, function(b) if (allzero(b)) NULL else b)
  list(A = if (allzero(A)) NULL else A,
       B = if (all(vapply(B, is.null, TRUE))) NULL else B,
       C = if (all(vapply(C, .is_zero, TRUE))) NULL else C)
}

#' Semilinear ODE functions and problems
#'
#' `SemilinearODEProblem(sys, op, tspan)` builds a split ODE problem
#' `du/dt = f1(u) + f2(u)` from the semiquadratic form of the system
#' ([calculate_semiquadratic_form()]): by default the linear part `A u` is
#' the stiff part `f1` and the quadratic and nonlinear parts form `f2`, as
#' ModelingToolkit's `SemilinearODEProblem`. The problem is an
#' `ODEProblem` too (with `f = f1 + f2`), so every ODE solver applies.
#'
#' @param sys A complete ODE `System`.
#' @param op Operating point.
#' @param tspan Time span.
#' @param stiff_linear,stiff_quadratic,stiff_nonlinear Which parts go into
#'   the stiff function `f1`.
#' @param jac Generate the Jacobian of the full right-hand side.
#' @param ... Passed to [ODEProblem()].
#' @return A `SemilinearODEProblem` (class `SplitODEProblem`, `ODEProblem`)
#'   or `SemilinearODEFunction`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); y <- Variable("y", t)
#' sys <- mtkcompile(System(list(D(x) ~ -2 * x + x * y, D(y) ~ -y), t, name = "s"))
#' prob <- SemilinearODEProblem(sys, list(x ~ 1, y ~ 1), c(0, 1))
#' prob$f$f1(c(1, 1), prob$p, 0)
#' sol <- solve(prob, Rodas5P())
#' @export
SemilinearODEFunction <- function(sys, stiff_linear = TRUE, stiff_quadratic = FALSE,
                                  stiff_nonlinear = FALSE, jac = FALSE, ...) {
  .need_complete(sys, "SemilinearODEFunction")
  sq <- calculate_semiquadratic_form(sys)
  A <- sq$A; B <- sq$B; C <- sq$C
  if (is.null(A) && is.null(B))
    stop("cannot generate split form for the system: it has no linear or quadratic part", call. = FALSE)
  if ((stiff_linear || is.null(A)) && (stiff_quadratic || is.null(B)) && stiff_nonlinear)
    stop("A, B and C cannot all be stiff", call. = FALSE)
  if ((!stiff_linear || is.null(A)) && (!stiff_quadratic || is.null(B)) && !stiff_nonlinear)
    stop("A, B and C cannot all be nonstiff", call. = FALSE)
  dvs <- unknowns(sys)
  neq <- length(equations(sys))
  zero <- .S0()
  part <- function(i, which) {
    out <- zero
    if (which == "A" && !is.null(A))
      for (j in seq_along(dvs)) out <- out + A[[i, j]] * dvs[[j]]
    if (which == "B" && !is.null(B) && !is.null(B[[i]]))
      for (j in seq_along(dvs)) for (k in j:length(dvs)) out <- out + B[[i]][[j, k]] * dvs[[j]] * dvs[[k]]
    if (which == "C" && !is.null(C)) out <- C[[i]]
    out
  }
  stiff <- c(A = stiff_linear, B = stiff_quadratic, C = stiff_nonlinear)
  build <- function(sel) {
    exprs <- lapply(seq_len(neq), function(i) {
      e <- zero
      for (w in names(stiff)[sel]) e <- e + part(i, w)
      e
    })
    .build_closure(exprs, .codegen_args(sys), observed = .codegen_observed(sys), out = "vector")
  }
  fn <- ODEFunction(sys, jac = jac)
  fn$f1 <- build(stiff)
  fn$f2 <- build(!stiff)
  fn$semiquadratic_form <- sq
  class(fn) <- c("SemilinearODEFunction", "SplitFunction", class(fn))
  fn
}

#' @rdname SemilinearODEFunction
#' @export
SemilinearODEProblem <- function(sys, op = list(), tspan = NULL, stiff_linear = TRUE, stiff_quadratic = FALSE,
                                 stiff_nonlinear = FALSE, jac = FALSE, ...) {
  .need_complete(sys, "SemilinearODEProblem")
  prob <- ODEProblem(sys, op, tspan, ..., jac = jac)
  sf <- SemilinearODEFunction(sys, stiff_linear = stiff_linear, stiff_quadratic = stiff_quadratic,
                              stiff_nonlinear = stiff_nonlinear, jac = jac)
  f1 <- sf$f1; f2 <- sf$f2
  sf$f <- function(u, p, t) f1(u, p, t) + f2(u, p, t)
  prob$f <- sf
  class(prob) <- c("SemilinearODEProblem", "SplitODEProblem", class(prob))
  prob
}

# ---------------------------------------------------------------------------
# generate_W

#' Generate the W matrix function
#'
#' `generate_W(sys)` returns an R function `(u, p, gamma, t)` computing
#' `gamma * M + J`, with `M` the mass matrix and `J` the Jacobian of the
#' right-hand side, the matrix the implicit solvers factorize (as
#' ModelingToolkit's `generate_W`).
#'
#' @param sys A complete ODE `System`.
#' @param ... Unused.
#' @return An R function.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); k <- Parameter("k")
#' sys <- mtkcompile(System(list(D(x) ~ -k * x^2), t, name = "s"))
#' W <- generate_W(sys)
#' W(2, c(k = 3), -1, 0)
#' @export
generate_W <- function(sys, ...) {
  .need_complete(sys, "generate_W")
  J <- .symmatrix_to_list_matrix(calculate_jacobian(sys))
  M <- as.matrix(calculate_massmatrix(sys))
  gam <- Symbol(".gamma_W")
  n <- nrow(M)
  W <- .list_matrix(n, n, function(i, j) gam * M[i, j] + J[[i, j]])
  args <- .codegen_args(sys)
  .build_closure(W, c(args[1:2], list(list(name = "gamma", vars = gam)), args[3L]), out = "matrix")
}
