# Boundary value problems, after ModelingToolkitBase/src/problems/bvproblem.jl
# (BVProblem), systems/codegen.jl (generate_boundary_conditions) and the
# solvers of BoundaryValueDiffEq.jl (Shooting, MultipleShooting from
# BoundaryValueDiffEqShooting; MIRK4, MIRK5, MIRK6 with defect control from
# BoundaryValueDiffEqMIRK, see the MIRK section below).
#
# Boundary conditions are the system's constraints, equations between
# evaluated variables such as EvalAt(0)(x) ~ 1 (x(0) ~ 1 in Julia), plus a
# hard initial condition x(t0) = value for every unknown given a value in the
# operating point; guesses only seed the solver.

#' Boundary value problems
#'
#' `BVProblem(sys, op, tspan, guesses = )` builds a boundary value problem
#' from a compiled system with constraints such as `EvalAt(0)(x) ~ 1` and
#' `EvalAt(1)(x) ~ 2`. Unknowns given values in `op` are fixed at the initial
#' time; `guesses` (and the system's guesses) start the iteration.
#'
#' `BVProblem(f, bc, u0, tspan, p)` builds one from an R function `f(u, p,
#' t)` and a residual `bc(sol, p, t)`, where `sol(s)` returns the state at
#' time `s`.
#'
#' Solve with `solve(bvp, MIRK4(), dt = 0.01)` (collocation with mesh
#' refinement from an initial step `dt`; also `MIRK5()`, `MIRK6()`),
#' `Shooting(Tsit5())` or `MultipleShooting(10, Tsit5())`. The result is a
#' `BVPSolution`, an ODE solution with dense output.
#'
#' @param sys A compiled `System`, or the function `f`.
#' @param ... For systems: `op`, `tspan`. For functions: `bc`, `u0`,
#'   `tspan`, `p`.
#' @param guesses Initial guesses for unknowns without fixed values.
#' @return A `BVProblem`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("y(t)")
#' sys <- mtkcompile(System(list(D(D(y)) ~ -y), t, name = "bvp",
#'   constraints = list(EvalAt(0)(y) ~ 0, EvalAt(pi / 2)(y) ~ 1)))
#' bvp <- BVProblem(sys, list(), c(0, pi / 2), guesses = list(y ~ 0, D(y) ~ 1))
#' sol <- solve(bvp, MIRK4(), dt = 0.05)
#' sol(pi / 4)[1]
#' @export
BVProblem <- function(sys, ..., guesses = NULL) {
  if (is.function(sys) || inherits(sys, "ODEFunction")) return(.bvp_problem_raw(sys, ...))
  .bvp_problem_sys(sys, ..., guesses = guesses)
}

.bvp_problem_raw <- function(f, bc, u0, tspan, p = numeric(), ...) {
  fn <- if (inherits(f, "ODEFunction")) f else ODEFunction(f)
  tspan <- .as_tspan(tspan)
  .new_problem("BVProblem", f = fn, bc = function(S, p) bc(S, p, tspan), u0 = as.numeric(unlist(u0)),
               tspan = tspan, p = p, sys = NULL, kwargs = list(...))
}

.bvp_problem_sys <- function(sys, op = list(), tspan = NULL, ..., guesses = NULL, callback = NULL) {
  .need_complete(sys, "BVProblem")
  if (!is.null(callback)) stop("BVP solvers do not support callbacks.", call. = FALSE)
  if (is_dde(sys)) stop("BVProblem does not support systems with delays", call. = FALSE)
  if (length(jumps(sys))) stop("BVProblem does not support systems with jumps", call. = FALSE)
  if (!is.null(.f(sys, "noise_eqs")) || length(.f(sys, "brownians")))
    stop("BVProblem does not support systems with noise", call. = FALSE)
  if (length(continuous_events(sys)) || length(discrete_events(sys)))
    stop("BVP solvers do not support events.", call. = FALSE)
  tspan <- .as_tspan(tspan, sys)
  dvs <- unknowns(sys)
  allv <- c(dvs, .parameters_of(sys, TRUE), observables(sys))
  opm <- .resolve_named(as_symmap(op), allv)
  gm <- .resolve_named(as_symmap(guesses), allv)
  is_full_ic <- has_alg_equations(sys) && !length(costs(sys))
  vals <- .process_op(sys, if (is_full_ic) opm else smap_merge(gm, opm))
  u0 <- vals$u0
  sg <- guesses(sys)
  for (i in which(is.na(u0))) {
    g <- smap_get(sg, dvs[[i]])
    u0[i] <- if (is.numeric(g) && length(g) == 1L) g else 0
  }
  u0_idxs <- if (is_full_ic) seq_along(dvs) else which(.keys(dvs) %in% names(opm))
  fn <- ODEFunction(sys)
  M <- fn$mass_matrix
  if (is.matrix(M) && !isTRUE(all.equal(M, diag(nrow(M)))))
    stop("boundary value problems with algebraic equations are not supported", call. = FALSE)
  bc <- .bvp_bc_fn(sys, u0, u0_idxs, tspan[1L])
  nbc <- length(u0_idxs) + length(constraints(sys))
  if (nbc > length(dvs))
    warning("The BVProblem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by op) exceeds the total number of states. The BVP solvers will default to doing a nonlinear least-squares optimization.",
            call. = FALSE)
  prob <- .new_problem("BVProblem", f = fn, bc = bc, u0 = u0, tspan = tspan, p = vals$p, sys = sys,
                       u0_idxs = u0_idxs, kwargs = list(...))
  prob$op <- vals$op
  prob
}

# The boundary residual bc(S, p): x(t0) - u0 for fixed unknowns, then the
# constraints lhs - rhs. After generate_boundary_conditions.
.bvp_bc_fn <- function(sys, u0, u0_idxs, t0) {
  cons <- constraints(sys)
  for (c in cons) if (!inherits(c, "Equation")) stop("BVP constraints must be equations", call. = FALSE)
  exprs <- lapply(cons, function(c) c$lhs - c$rhs)
  pargs <- list(list(name = "p", vars = .problem_params(sys)))
  cg <- if (length(exprs)) .evalat_codegen(exprs, sys, pargs, .param_binding_eqs(sys)) else NULL
  fixed <- u0[u0_idxs]
  function(S, p) {
    r0 <- if (length(u0_idxs)) S(t0)[u0_idxs] - fixed else numeric()
    if (is.null(cg)) return(r0)
    vals <- if (cg$n) .evalat_values(cg$tfun(p), cg$idx, S) else numeric()
    c(r0, cg$fn(p, vals))
  }
}

#' @export
print.BVProblem <- function(x, ...) {
  cat("BVProblem", if (!is.null(x$sys)) paste0(" of system ", nameof(x$sys)), "\n", sep = "")
  cat("Time span: [", paste(x$tspan, collapse = ", "), "]\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Solvers

#' Boundary value problem solvers
#'
#' `Shooting(ode_alg)` solves for the initial state with Newton's method,
#' integrating the ODE with `ode_alg` for each residual. `MultipleShooting(
#' nshoots, ode_alg)` shoots from `nshoots` nodes and adds continuity
#' conditions. `MIRK4()`, `MIRK5()` and `MIRK6()` are the mono-implicit
#' Runge-Kutta collocation methods of orders 4, 5 and 6 of
#' BoundaryValueDiffEq.jl with their continuous extensions and defect
#' control: starting from a uniform mesh of step `dt`, the mesh is refined
#' until the relative defect of the continuous solution is below `abstol`.
#' Options of `solve()`: `dt` (initial mesh step of the MIRK methods,
#' required), `abstol` (defect tolerance, 1e-6; for shooting the
#' nonlinear tolerance), `adaptive = FALSE` (keep the initial mesh),
#' `nlsolve_kwargs = list(abstol = )` (Newton tolerance of the MIRK
#' methods, 1e-10), `maxiters`, and for shooting `reltol` and `abstol` of
#' the inner ODE solves through `odesolve_kwargs`.
#'
#' @param ode_alg An ODE algorithm.
#' @param nshoots Number of shooting intervals.
#' @return An algorithm object.
#' @examples
#' Shooting(Tsit5())
#' # y'' = 100 y with boundary layers: the mesh is refined from 4 intervals
#' bvp <- BVProblem(function(u, p, t) c(u[2], 100 * u[1]),
#'                  function(sol, p, t) c(sol(0)[1] - 1, sol(1)[1] - 1), c(1, 0), c(0, 1))
#' sol <- solve(bvp, MIRK4(), dt = 0.25)
#' length(sol$t)
#' sol(0.5)[1] - 1 / cosh(5)
#' @name bvp_algorithms
NULL

#' @rdname bvp_algorithms
#' @export
Shooting <- function(ode_alg = Tsit5()) structure(list(name = "Shooting", ode_alg = ode_alg),
                                                  class = c("Shooting", "BVPAlgorithm"))
#' @rdname bvp_algorithms
#' @export
MultipleShooting <- function(nshoots = 10L, ode_alg = Tsit5()) {
  structure(list(name = "MultipleShooting", nshoots = as.integer(nshoots), ode_alg = ode_alg),
            class = c("MultipleShooting", "BVPAlgorithm"))
}
.new_mirk <- function(name, defect_threshold, max_num_subintervals) {
  structure(list(name = name, defect_threshold = defect_threshold,
                 max_num_subintervals = as.integer(max_num_subintervals)),
            class = c(name, "MIRK", "BVPAlgorithm"))
}

#' @rdname bvp_algorithms
#' @param defect_threshold Largest acceptable relative defect; above it the
#'   mesh is halved and the iteration restarts.
#' @param max_num_subintervals Largest number of mesh intervals.
#' @export
MIRK4 <- function(defect_threshold = 0.1, max_num_subintervals = 3000L)
  .new_mirk("MIRK4", defect_threshold, max_num_subintervals)

#' @rdname bvp_algorithms
#' @export
MIRK5 <- function(defect_threshold = 0.1, max_num_subintervals = 3000L)
  .new_mirk("MIRK5", defect_threshold, max_num_subintervals)

#' @rdname bvp_algorithms
#' @export
MIRK6 <- function(defect_threshold = 0.1, max_num_subintervals = 3000L)
  .new_mirk("MIRK6", defect_threshold, max_num_subintervals)

#' @export
print.BVPAlgorithm <- function(x, ...) {
  cat(x$name, "()\n", sep = "")
  invisible(x)
}

# Forward difference Jacobian with a step suited to residuals that hold
# integration error.
.bvp_fd_jac <- function(F, x, fx, rel = 1e-6) {
  J <- matrix(0, length(fx), length(x))
  for (j in seq_along(x)) {
    h <- rel * max(1, abs(x[j]))
    xh <- x
    xh[j] <- x[j] + h
    J[, j] <- (F(xh) - fx) / h
  }
  J
}

.bvp_newton <- function(F, x0, tol, maxiters, jac) {
  res <- .newton_solve(F, x0, tol = tol, maxiter = maxiters, jac = jac)
  if (!res$converged) warning("the boundary value solver did not converge", call. = FALSE)
  res
}

#' @rdname bvp_algorithms
#' @param a A `BVProblem`.
#' @param b A BVP algorithm; `MIRK4()` by default.
#' @param ... Solver options.
#' @export
solve.BVProblem <- function(a, b = NULL, ...) {
  alg <- b %||% MIRK4()
  opts <- utils::modifyList(a$kwargs %||% list(), list(...))
  nbc <- length(a$bc(function(s) a$u0, a$p))
  if (nbc < length(a$u0))
    .stopf("the boundary value problem is underdetermined: %d conditions for %d states", nbc, length(a$u0))
  sol <- if (inherits(alg, "MIRK")) .bvp_mirk(a, alg, opts)
         else if (inherits(alg, "MultipleShooting")) .bvp_multiple_shooting(a, alg, opts)
         else if (inherits(alg, "Shooting")) .bvp_shooting(a, alg, opts)
         else stop("not a BVP algorithm: ", class(alg)[1L], call. = FALSE)
  env <- .sol_env(sol)
  env$prob <- a
  env$alg <- alg
  class(sol) <- c("BVPSolution", class(sol))
  sol
}

.bvp_ode_solve <- function(a, alg, u0, tspan, opts) {
  ok <- utils::modifyList(list(reltol = 1e-10, abstol = 1e-12), opts$odesolve_kwargs %||% list())
  prob <- .new_problem("ODEProblem", f = a$f, u0 = stats::setNames(u0, names(a$u0)), tspan = tspan,
                       p = a$p, sys = a$sys, kwargs = list())
  .integrate(prob, alg$ode_alg, ok)
}

.bvp_shooting <- function(a, alg, opts) {
  F <- function(x) {
    s <- .bvp_ode_solve(a, alg, x, a$tspan, opts)
    a$bc(function(tt) as.numeric(s(tt)), a$p)
  }
  res <- .bvp_newton(F, a$u0, tol = opts$abstol %||% 1e-8, maxiters = opts$maxiters %||% 50L,
                     jac = function(x) .bvp_fd_jac(F, x, F(x)))
  sol <- .bvp_ode_solve(a, alg, res$x, a$tspan, opts)
  if (!res$converged) {
    env <- .sol_env(sol)
    env$retcode <- "Failure"
  }
  sol
}

.bvp_multiple_shooting <- function(a, alg, opts) {
  n <- length(a$u0)
  N <- alg$nshoots
  nodes <- seq(a$tspan[1L], a$tspan[2L], length.out = N + 1L)
  seg <- function(k, xk) .bvp_ode_solve(a, alg, xk, nodes[c(k, k + 1L)], opts)
  block <- function(X, k) X[(k - 1L) * n + seq_len(n)]
  S_of <- function(segs) function(s) {
    k <- min(N, max(1L, findInterval(s, nodes, rightmost.closed = TRUE)))
    as.numeric(segs[[k]](s))
  }
  resid <- function(segs, X) {
    cont <- unlist(lapply(seq_len(N - 1L), function(k) as.numeric(segs[[k]](nodes[k + 1L])) - block(X, k + 1L)))
    c(cont, a$bc(S_of(segs), a$p))
  }
  all_segs <- function(X) lapply(seq_len(N), function(k) seg(k, block(X, k)))
  F <- function(X) resid(all_segs(X), X)
  jac <- function(X) {
    segs <- all_segs(X)
    r0 <- resid(segs, X)
    J <- matrix(0, length(r0), length(X))
    for (k in seq_len(N)) for (i in seq_len(n)) {
      j <- (k - 1L) * n + i
      h <- 1e-6 * max(1, abs(X[j]))
      X2 <- X
      X2[j] <- X[j] + h
      s2 <- segs
      s2[[k]] <- seg(k, block(X2, k))
      J[, j] <- (resid(s2, X2) - r0) / h
    }
    J
  }
  # initial nodes from one integration of the initial guess
  s0 <- .bvp_ode_solve(a, alg, a$u0, a$tspan, opts)
  X0 <- unlist(lapply(seq_len(N), function(k) as.numeric(s0(nodes[k]))))
  if (!all(is.finite(X0))) X0 <- rep(a$u0, N)
  res <- .bvp_newton(F, X0, tol = opts$abstol %||% 1e-8, maxiters = opts$maxiters %||% 50L, jac = jac)
  segs <- all_segs(res$x)
  # join the segments into one dense solution
  steps <- unlist(lapply(segs, function(s) .sol_env(s)$steps), recursive = FALSE)
  ts <- c(unlist(lapply(segs, function(s) utils::head(s$t, -1L))), a$tspan[2L])
  U <- rbind(do.call(rbind, lapply(segs, function(s) utils::head(as.matrix(s), -1L))),
             as.numeric(segs[[N]](a$tspan[2L])))
  colnames(U) <- names(a$u0)
  .new_ode_solution(t = ts, u = U, prob = a, alg = alg, retcode = if (res$converged) "Success" else "Failure",
                    steps = steps, stats = list(), dense = TRUE)
}

# ---------------------------------------------------------------------------
# MIRK collocation with defect control (BoundaryValueDiffEqMIRK: mirk.jl
# __perform_mirk_iteration, collocation.jl Phi!, adaptivity.jl
# error_estimate!(DefectControl), mesh_selector!, redistribute!, half_mesh!,
# interpolation.jl interp_setup!, interp_weights, EvalSol;
# mirk_tableaus.jl constructMIRK4/5/6; BoundaryValueDiffEq.jl commit b662553).
#
# On a mesh t_1 < ... < t_(N+1) with h_i = t_(i+1) - t_i the stages of
# interval i are K_r = f((1 - v_r) y_i + v_r y_(i+1) + h_i sum_s x_rs K_s,
# t_i + c_r h_i) and the residual is y_(i+1) - y_i - h_i sum_r b_r K_r. The
# continuous extension adds s_star - stage stages (c_star, v_star, x_star)
# and evaluates y_i + h_i sum_j w_j(tau) K_j. Boundary conditions see that
# interpolant. The defect (z' - f(z)) / (|f(z)| + 1) at tau_star and
# 1 - tau_star of each interval drives the mesh: above defect_threshold the
# mesh is halved and the guess restarts from zero, above abstol the mesh is
# redistributed to equidistribute the defect (or halved when it is already
# even) and the old interpolant gives the new guess.

.MIRK_TABS <- list(
  MIRK4 = list(order = 4L, stage = 3L, c = c(0, 1, 1 / 2), v = c(0, 1, 1 / 2), b = c(1 / 6, 1 / 6, 2 / 3),
               x = rbind(c(0, 0, 0), c(0, 0, 0), c(1 / 8, -1 / 8, 0)),
               s_star = 4L, c_star = 3 / 4, v_star = 27 / 32,
               x_star = rbind(c(3 / 64, -9 / 64, 0, 0)), tau_star = 0.226),
  MIRK5 = list(order = 5L, stage = 4L, c = c(0, 1, 3 / 4, 3 / 10), v = c(0, 1, 27 / 32, 837 / 1250),
               b = c(5 / 54, 1 / 14, 32 / 81, 250 / 567),
               x = rbind(c(0, 0, 0, 0), c(0, 0, 0, 0), c(3 / 64, -9 / 64, 0, 0),
                         c(21 / 1000, 63 / 5000, -252 / 625, 0)),
               s_star = 6L, c_star = c(4 / 5, 13 / 23), v_star = c(4 / 5, 13 / 23),
               x_star = rbind(c(14 / 1125, -74 / 875, -128 / 3375, 104 / 945, 0, 0),
                              c(1 / 2, 4508233 / 1958887, 48720832 / 2518569, -27646420 / 17629983,
                                -11517095 / 559682, 0)),
               tau_star = 0.3),
  MIRK6 = list(order = 6L, stage = 5L, c = c(0, 1, 1 / 4, 3 / 4, 1 / 2), v = c(0, 1, 5 / 32, 27 / 32, 1 / 2),
               b = c(7 / 90, 7 / 90, 16 / 45, 16 / 45, 2 / 15),
               x = rbind(c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0), c(9 / 64, -3 / 64, 0, 0, 0),
                         c(3 / 64, -9 / 64, 0, 0, 0), c(-5 / 24, 5 / 24, 2 / 3, -2 / 3, 0)),
               s_star = 9L, c_star = c(7 / 16, 3 / 8, 9 / 16, 1 / 8), v_star = c(7 / 16, 3 / 8, 9 / 16, 1 / 8),
               x_star = rbind(c(1547 / 32768, -1225 / 32768, 749 / 4096, -287 / 2048, -861 / 16384, 0, 0, 0, 0),
                              c(83 / 1536, -13 / 384, 283 / 1536, -167 / 1536, -49 / 512, 0, 0, 0, 0),
                              c(1225 / 32768, -1547 / 32768, 287 / 2048, -749 / 4096, 861 / 16384, 0, 0, 0, 0),
                              c(233 / 3456, -19 / 1152, 0, 0, 0, -5 / 72, 7 / 72, -17 / 216, 0)),
               tau_star = 0.7156))

# interp_weights: the weights w(tau) of the continuous extension and their
# derivatives, as polynomial coefficients in tau (constant term first).
.mirk_poly <- function(co, tau) as.numeric(co %*% tau^(seq_len(ncol(co)) - 1L))
.MIRK_W <- list(
  MIRK5 = list(
    w = rbind(c(0, 1, -83910 / 22464, 143041 / 22464, -113808 / 22464, 33256 / 22464),
              c(0, 0, -2418 / 3360, 12303 / 3360, -19512 / 3360, 10904 / 3360),
              -8 / 81 * c(0, 0, -78, 209, -204, 8),
              -25 / 1134 * c(0, 0, -390, 1045, -1020, 328),
              -25 / 5184 * c(0, 0, 390, 255, -1680, 2072),
              279841 / 168480 * c(0, 0, -6, 21, -24, 8)),
    wp = rbind(c(1, -13985 / 1872, 143041 / 7488, -2371 / 117, 20785 / 2808),
               c(0, -403 / 280, 12303 / 1120, -813 / 35, 1363 / 84),
               c(0, 416 / 27, -1672 / 27, 2176 / 27, -320 / 81),
               c(0, 3250 / 189, -26125 / 378, 17000 / 189, -20500 / 567),
               c(0, -1625 / 432, -2125 / 576, 875 / 27, -32375 / 648),
               c(0, -279841 / 14040, 1958887 / 18720, -279841 / 1755, 279841 / 4212))),
  MIRK6 = list(
    w = rbind(c(0, 1, -28607 / 7434, -166210 / 33453, 334780 / 11151, -1911296 / 55755, 406528 / 33453),
              c(0, 0, 777 / 590, -2534158 / 234171, 2088580 / 78057, -10479104 / 390285, 11328512 / 1170855),
              c(0, 0, -1008 / 59, 222176 / 1593, -180032 / 531, 876544 / 2655, -180224 / 1593),
              c(0, 0, -1008 / 59, 222176 / 1593, -180032 / 531, 876544 / 2655, -180224 / 1593),
              c(0, 0, -378 / 59, 27772 / 531, -22504 / 177, 109568 / 885, -22528 / 531),
              c(0, 0, -95232 / 413, 62384128 / 33453, -49429504 / 11151, 46759936 / 11151, -46661632 / 33453),
              c(0, 0, 896 / 5, -4352 / 3, 3456, -16384 / 5, 16384 / 15),
              c(0, 0, 50176 / 531, -179554304 / 234171, 143363072 / 78057, -136675328 / 78057,
                137363456 / 234171),
              c(0, 0, 0, 16384 / 441, -16384 / 147, 16384 / 147, -16384 / 441)),
    wp = rbind(c(1, -28607 / 3717, -166210 / 11151, 1339120 / 11151, -1911296 / 11151, 813056 / 11151),
               c(0, 777 / 295, -2534158 / 78057, 8354320 / 78057, -10479104 / 78057, 22657024 / 390285),
               c(0, -2016 / 59, 222176 / 531, -720128 / 531, 876544 / 531, -360448 / 531),
               c(0, -2016 / 59, 222176 / 531, -720128 / 531, 876544 / 531, -360448 / 531),
               c(0, -756 / 59, 27772 / 177, -90016 / 177, 109568 / 177, -45056 / 177),
               c(0, -190464 / 413, 62384128 / 11151, -197718016 / 11151, 233799680 / 11151, -93323264 / 11151),
               c(0, 1792 / 5, -4352, 13824, -16384, 32768 / 5),
               c(0, 100352 / 531, -179554304 / 78057, 573452288 / 78057, -683376640 / 78057, 274726912 / 78057),
               c(0, 0, 16384 / 147, -65536 / 147, 81920 / 147, -32768 / 147))))

.mirk_weights <- function(name, tau) {
  if (name == "MIRK4") {
    t2 <- tau * tau; tm1 <- tau - 1; t4m3 <- 4 * tau - 3; t2m1 <- 2 * tau - 1
    return(list(w = c(-tau * (2 * tau - 3) * (2 * t2 - 3 * tau + 2) / 6, t2 * (12 * t2 - 20 * tau + 9) / 6,
                      2 * t2 * (6 * t2 - 14 * tau + 9) / 3, -16 * t2 * tm1 * tm1 / 3),
                wp = c(-tm1 * t4m3 * t2m1 / 3, tau * t2m1 * t4m3, 4 * tau * t4m3 * tm1, -32 * tau * t2m1 * tm1 / 3)))
  }
  W <- .MIRK_W[[name]]
  list(w = .mirk_poly(W$w, tau), wp = .mirk_poly(W$wp, tau))
}

# Stages of every interval: a list of n x N matrices, the discrete stages
# then (extra = TRUE) the stages of the continuous extension. Stages at an
# end of the interval (v = c = 0 or 1 with no coupling) reuse f at the
# nodes, `Fn`.
.mirk_stages <- function(tab, f, p, mesh, h, Ym, Fn = NULL, extra = FALSE, idx = NULL) {
  N <- length(h)
  idx <- idx %||% seq_len(N)
  n <- nrow(Ym)
  Y0 <- Ym[, idx, drop = FALSE]
  Y1 <- Ym[, idx + 1L, drop = FALSE]
  H <- matrix(h[idx], n, length(idx), byrow = TRUE)
  t0 <- mesh[idx]
  feval <- function(G, tt) matrix(vapply(seq_along(tt), function(i) as.numeric(f(G[, i], p, tt[i])), numeric(n)),
                                  n, length(tt))
  node_f <- function(j) if (!is.null(Fn)) Fn[, j, drop = FALSE] else feval(Ym[, j, drop = FALSE], mesh[j])
  K <- vector("list", tab$stage + if (extra) tab$s_star - tab$stage else 0L)
  for (r in seq_len(tab$stage)) {
    if (all(tab$x[r, ] == 0) && tab$v[r] == tab$c[r] && tab$c[r] %in% c(0, 1)) {
      K[[r]] <- node_f(idx + as.integer(tab$c[r]))
      next
    }
    G <- (1 - tab$v[r]) * Y0 + tab$v[r] * Y1
    for (s in seq_len(r - 1L)) if (tab$x[r, s] != 0) G <- G + H * (tab$x[r, s] * K[[s]])
    K[[r]] <- feval(G, t0 + tab$c[r] * h[idx])
  }
  if (extra) for (r in seq_len(tab$s_star - tab$stage)) {
    G <- (1 - tab$v_star[r]) * Y0 + tab$v_star[r] * Y1
    xs <- tab$x_star[r, ]
    for (j in seq_len(tab$stage + r - 1L)) if (xs[j] != 0) G <- G + H * (xs[j] * K[[j]])
    K[[tab$stage + r]] <- feval(G, t0 + tab$c_star[r] * h[idx])
  }
  K
}

# Phi!: the collocation residuals, an n x N matrix.
.mirk_residual <- function(tab, Ym, h, K) {
  N <- length(h)
  S <- 0
  for (r in seq_len(tab$stage)) S <- S + tab$b[r] * K[[r]]
  Ym[, -1L, drop = FALSE] - Ym[, -(N + 1L), drop = FALSE] - matrix(h, nrow(Ym), N, byrow = TRUE) * S
}

# The continuous extension on interval i from its stages (an n x s_star
# matrix): the value at s, and its derivative with deriv = TRUE.
.mirk_eval <- function(name, Kall, y0, h, tau, deriv = FALSE) {
  w <- .mirk_weights(name, tau)
  if (deriv) as.numeric(Kall %*% w$wp) else y0 + h * as.numeric(Kall %*% w$w)
}

# interval(): t in (t_i, t_(i+1)] is in interval i, clamped to 1..N.
.mirk_interval <- function(mesh, s) {
  N <- length(mesh) - 1L
  min(N, max(1L, sum(mesh < s)))
}

# EvalSol: the state function handed to the boundary conditions. Stages
# of the intervals it touches are computed on demand.
.mirk_evalsol <- function(tab, name, f, p, mesh, h, Ym, touched = NULL) {
  cache <- list()
  N <- length(h)
  function(s) {
    if (s == mesh[1L]) return(Ym[, 1L])
    if (s == mesh[N + 1L]) return(Ym[, N + 1L])
    i <- .mirk_interval(mesh, s)
    if (!is.null(touched)) touched$i <- c(touched$i, i)
    key <- as.character(i)
    Kall <- cache[[key]]
    if (is.null(Kall)) {
      Kall <- do.call(cbind, .mirk_stages(tab, f, p, mesh, h, Ym, extra = TRUE, idx = i))
      cache[[key]] <<- Kall
    }
    .mirk_eval(name, Kall, Ym[, i], h[i], (s - mesh[i]) / h[i])
  }
}

# error_estimate!(DefectControl): per interval the larger of the relative
# defects at tau_star and 1 - tau_star (max norm per interval).
.mirk_defect <- function(tab, name, f, p, mesh, h, Ym, Kall) {
  N <- length(h)
  n <- nrow(Ym)
  est <- function(tau) {
    w <- .mirk_weights(name, tau)
    Z <- Ym[, -(N + 1L), drop = FALSE]
    Zp <- 0
    for (j in seq_along(Kall)) {
      Z <- Z + matrix(h, n, N, byrow = TRUE) * (w$w[j] * Kall[[j]])
      Zp <- Zp + w$wp[j] * Kall[[j]]
    }
    tt <- mesh[-(N + 1L)] + tau * h
    Fz <- matrix(vapply(seq_len(N), function(i) as.numeric(f(Z[, i], p, tt[i])), numeric(n)), n, N)
    apply(abs((Zp - Fz) / (abs(Fz) + 1)), 2L, max)
  }
  pmax(est(tab$tau_star), est(1 - tab$tau_star))
}

# half_mesh!
.mirk_half_mesh <- function(mesh) {
  N <- length(mesh) - 1L
  out <- numeric(2L * N + 1L)
  out[seq(1L, 2L * N + 1L, by = 2L)] <- mesh
  out[seq(2L, 2L * N, by = 2L)] <- (mesh[-1L] + mesh[-(N + 1L)]) / 2
  out
}

# redistribute!: Nsub intervals with equal integrals of the density s_hat.
.mirk_redistribute <- function(mesh, h, shat, Nsub) {
  N <- length(h)
  zeta <- sum(shat * h) / Nsub
  out <- numeric(Nsub + 1L)
  out[1L] <- mesh[1L]
  t <- mesh[1L]
  k <- 1L
  i <- 0L
  integral <- 0
  while (k <= N) {
    nxt <- integral + shat[k] * (mesh[k + 1L] - t)
    if (nxt > zeta) {
      if (i + 2L > Nsub) break
      out[i + 2L] <- (zeta - integral) / shat[k] + t
      t <- out[i + 2L]
      i <- i + 1L
      integral <- 0
    } else {
      integral <- nxt
      t <- mesh[k + 1L]
      k <- k + 1L
    }
  }
  out[Nsub + 1L] <- mesh[N + 1L]
  out
}

# mesh_selector!(DefectControl): the new mesh, or NULL when it would exceed
# max_num_subintervals.
.mirk_mesh_selector <- function(mesh, h, defects, order, abstol, max_sub) {
  Nn <- length(mesh)
  n <- Nn - 1L
  shat <- (defects / abstol)^(1 / (order + 1))
  r1 <- max(shat); r2 <- sum(shat); r3 <- r2 / n
  n_predict <- round(1.3 * r2 + 1)
  if (abs(n_predict - n) < 0.1 * n) n_predict <- round(n + 0.1 * n)
  if (r1 <= r3) {
    if (2L * n > max_sub) return(NULL)
    return(.mirk_half_mesh(mesh))
  }
  Nsub <- as.integer(min(max(n_predict, Nn %/% 2L), 4L * n))
  if (Nsub > max_sub) return(NULL)
  .mirk_redistribute(mesh, h, shat / h, Nsub)
}

.bvp_mirk <- function(a, alg, opts) {
  dt <- opts$dt
  if (is.null(dt) || !is.numeric(dt) || dt <= 0)
    stop(alg$name, "() needs a positive mesh step: solve(bvp, ", alg$name, "(), dt = ...)", call. = FALSE)
  tab <- .MIRK_TABS[[alg$name]]
  name <- alg$name
  t0 <- a$tspan[1L]; tf <- a$tspan[2L]
  mesh <- seq(t0, tf, length.out = ceiling((tf - t0) / dt - 1e-10) + 1L)
  n <- length(a$u0)
  f <- a$f$f
  p <- a$p
  abstol <- opts$abstol %||% 1e-6
  adaptive <- opts$adaptive %||% TRUE
  # the Newton tolerance on the residual (Julia passes abstol); kept tight
  # so that the discrete solution is solved well below the defect tolerance
  ntol <- opts$nlsolve_kwargs$abstol %||% min(abstol, 1e-10)
  maxiters <- opts$maxiters %||% 50L
  Y <- rep(a$u0, length(mesh))
  iterations <- 0L
  repeat {
    N <- length(mesh) - 1L
    h <- diff(mesh)
    fn_nodes <- function(Ym) matrix(vapply(seq_len(N + 1L), function(i) as.numeric(f(Ym[, i], p, mesh[i])), numeric(n)),
                                    n, N + 1L)
    colloc <- function(Ym, Fn) .mirk_residual(tab, Ym, h, .mirk_stages(tab, f, p, mesh, h, Ym, Fn = Fn))
    bcres <- function(Ym, touched = NULL) a$bc(.mirk_evalsol(tab, name, f, p, mesh, h, Ym, touched), p)
    F <- function(Y) {
      Ym <- matrix(Y, n, N + 1L)
      c(as.numeric(colloc(Ym, fn_nodes(Ym))), bcres(Ym))
    }
    jac <- function(Y) {
      Ym <- matrix(Y, n, N + 1L)
      Fn <- fn_nodes(Ym)
      C0 <- colloc(Ym, Fn)
      touched <- new.env(parent = emptyenv())
      touched$i <- integer()
      b0 <- bcres(Ym, touched)
      J <- matrix(0, N * n + length(b0), (N + 1L) * n)
      hstep <- 1e-7 * pmax(abs(Ym), 1)
      # collocation rows: interval i depends on nodes i and i + 1 only, so
      # perturbing every other node at once recovers all columns
      for (parity in 0:1) {
        nodes <- which((seq_len(N + 1L) - 1L) %% 2L == parity)
        for (comp in seq_len(n)) {
          Y2 <- Ym
          Y2[comp, nodes] <- Ym[comp, nodes] + hstep[comp, nodes]
          F2 <- Fn
          F2[, nodes] <- fn_nodes(Y2)[, nodes, drop = FALSE]
          D <- colloc(Y2, F2) - C0
          for (i in seq_len(N)) {
            node <- if ((i - 1L) %% 2L == parity) i else i + 1L
            J[(i - 1L) * n + seq_len(n), (node - 1L) * n + comp] <- D[, i] / hstep[comp, node]
          }
        }
      }
      # boundary rows: the nodes of the intervals the conditions evaluate
      # (the end values are exact nodes)
      bnodes <- sort(unique(c(1L, N + 1L, touched$i, touched$i + 1L)))
      for (node in bnodes) for (comp in seq_len(n)) {
        Y2 <- Ym
        Y2[comp, node] <- Ym[comp, node] + hstep[comp, node]
        J[N * n + seq_along(b0), (node - 1L) * n + comp] <- (bcres(Y2) - b0) / hstep[comp, node]
      }
      J
    }
    res <- .newton_solve(F, Y, tol = ntol, maxiter = maxiters, jac = jac)
    iterations <- iterations + res$iter
    Y <- res$x
    ok <- isTRUE(res$converged) && all(is.finite(Y))
    error_norm <- 2 * abstol
    if (!adaptive) break
    Ym <- matrix(Y, n, N + 1L)
    Kall <- NULL
    if (ok) {
      Kall <- .mirk_stages(tab, f, p, mesh, h, Ym, extra = TRUE)
      defects <- .mirk_defect(tab, name, f, p, mesh, h, Ym, Kall)
      error_norm <- max(defects)
      ok <- is.finite(error_norm) && error_norm <= alg$defect_threshold
    }
    if (ok) {
      if (error_norm <= abstol) break
      new_mesh <- .mirk_mesh_selector(mesh, h, defects, tab$order, abstol, alg$max_num_subintervals)
      if (is.null(new_mesh)) {
        ok <- FALSE
        break
      }
      # the new guess from the interpolant on the old mesh
      Y <- as.numeric(vapply(new_mesh, function(s) {
        i <- .mirk_interval(mesh, s)
        .mirk_eval(name, do.call(cbind, lapply(Kall, function(K) K[, i])), Ym[, i], h[i], (s - mesh[i]) / h[i])
      }, numeric(n)))
      mesh <- new_mesh
    } else {
      # the nonlinear solve failed or the defect is too large: halve the mesh
      # and restart from a zero guess
      if (2L * N > alg$max_num_subintervals) break
      mesh <- .mirk_half_mesh(mesh)
      Y <- rep(0, n * length(mesh))
    }
  }
  N <- length(mesh) - 1L
  h <- diff(mesh)
  Ym <- matrix(Y, n, N + 1L)
  Kall <- .mirk_stages(tab, f, p, mesh, h, Ym, extra = TRUE)
  if (!ok) warning("the boundary value solver did not converge", call. = FALSE)
  steps <- lapply(seq_len(N), function(i)
    list(t0 = mesh[i], dt = h[i], y0 = Ym[, i], y1 = Ym[, i + 1L],
         k = list(name = name, K = do.call(cbind, lapply(Kall, function(K) K[, i]))), interp = .mirk_interp))
  U <- t(Ym)
  colnames(U) <- names(a$u0)
  .new_ode_solution(t = mesh, u = U, prob = a, alg = alg, retcode = if (ok) "Success" else "Failure",
                    steps = steps, stats = list(iterations = iterations, nmesh = N + 1L), dense = TRUE)
}

.mirk_interp <- function(theta, dt, y0, y1, k) .mirk_eval(k$name, k$K, y0, dt, theta)
