# Linearization, after ModelingToolkit/src/linearization.jl
# (LinearizationOpPoint, linearization_function, _hold_opened_signals,
# LinearizationFunction, LinearizationProblem, solve(::LinearizationProblem),
# linearize_symbolic, linearize, similarity_transform, reorder_unknowns) and
# src/systems/analysis_points.jl (the analysis point method of
# linearization_function).
#
# Jacobians are symbolic (the right-hand sides with observed variables and
# parameter bindings substituted, differentiated with .partials) and compiled
# to R closures, where ModelingToolkit uses ForwardDiff.

#' Operating point from a solution
#'
#' `LinearizationOpPoint(sol, t, op)` uses the values of the differential
#' unknowns and parameters of the ODE solution `sol` at time `t` as the
#' operating point of [linearize()]; entries of `op` override them and may
#' be expressions, evaluated from the solution. With a vector `t`,
#' `linearize()` returns one linearization per time point. Initial
#' conditions of the system are ignored.
#'
#' @param sol An ODE solution.
#' @param t A time or vector of times.
#' @param op Further operating point values.
#' @return An object of class `LinearizationOpPoint`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); u <- Variable("u", t, input = TRUE)
#' sys <- System(list(D(x) ~ -x^3 + u), t, name = "s")
#' prob <- ODEProblem(mtkcompile(sys, inputs = list(u)), list(x ~ 1, u ~ 0), c(0, 1))
#' sol <- solve(prob)
#' lin <- linearize(sys, list(u), list(x), op = LinearizationOpPoint(sol, 1))
#' lin$A
#' @export
LinearizationOpPoint <- function(sol, t, op = list()) {
  structure(list(sol = sol, t = t, op = op), class = "LinearizationOpPoint")
}

# _build_op_from_solution for one time point: differential unknowns and
# parameters of the solution, then `op` (expressions evaluated from sol).
.build_op_from_solution <- function(opp, ti) {
  sol <- opp$sol
  prob <- sol$prob
  ssys <- prob$sys
  if (is.null(ssys)) stop("LinearizationOpPoint needs a solution of a problem built from a system", call. = FALSE)
  out <- symmap()
  us <- unknowns(ssys)
  uval <- as.numeric(sol(ti))
  for (i in seq_along(equations(ssys))) {
    if (.is_derivative(equations(ssys)[[i]]$lhs)) out <- smap_set(out, us[[i]], uval[i])
  }
  pp <- .problem_params(ssys)
  pt <- .sol_p_at(sol, ti)
  for (i in seq_along(pp)) out <- smap_set(out, pp[[i]], unname(pt[i]))
  extra <- .resolve_named(as_symmap(opp$op), c(us, .parameters_of(ssys, TRUE), observables(ssys)))
  for (el in unclass(extra)) {
    out <- smap_set(out, el$var, .sol_value_at(sol, el$value, ti))
  }
  out
}

.sol_value_at <- function(sol, expr, ti) {
  if (!is_sym(expr) || .is_num(expr)) return(if (is_sym(expr)) as.numeric(expr) else expr)
  ssys <- sol$prob$sys
  fn <- build_explicit_observed_function(ssys, expr)
  fn(as.numeric(sol(ti)), .sol_p_at(sol, ti), ti)
}

# Parameter values in force at time ti, after any event that changed them.
.sol_p_at <- function(sol, ti) {
  p <- .sol_params_at(.sol_env(sol), ti)
  if (is.matrix(p)) p[1L, ] else p
}

#' Linearization functions
#'
#' `linearization_function(sys, inputs, outputs)` compiles `sys` with the
#' given inputs and outputs and returns a list with `lin_fun`, a function
#' of `(u, p, t)` returning the Jacobian blocks `f_x`, `f_z`, `g_x`, `g_z`,
#' `f_u`, `g_u`, `h_x`, `h_z`, `h_u` of
#' \deqn{x' = f(x, z, u), \quad 0 = g(x, z, u), \quad y = h(x, z, u)}
#' at a consistent point (algebraic unknowns are solved for first), and the
#' simplified system `sys`. `u` and `p` are numeric vectors ordered as
#' `unknowns(sys)` and the problem parameters of `sys`; inputs have become
#' parameters. [linearize()] is the higher level interface.
#'
#' When `inputs` are analysis points (or their names) they are transformed
#' with [linearization_ap_transform()] first.
#'
#' @param sys A `System` (not compiled).
#' @param inputs,outputs Lists of variables, or analysis points.
#' @param op Operating point: values of unknowns, parameters and inputs.
#' @param t Time of the linearization.
#' @param simplify Passed to [mtkcompile()].
#' @param initialize Solve for consistent algebraic unknowns.
#' @param guesses Initial guesses for algebraic unknowns.
#' @param p Further parameter values.
#' @param zero_dummy_der Set the values of unknowns added by compilation
#'   (dummy derivatives) to zero.
#' @param ignore_system_initial_conditions Ignore the initial conditions of
#'   the system and use only `op`.
#' @param loop_openings Analysis points whose connections are removed.
#' @param loop_opening_params,loop_opening_drivers Opened signals and the
#'   variables that drove them (set by the analysis functions).
#' @param system_modifier A function applied to the system before
#'   compilation.
#' @param verbose Warn about an empty operating point.
#' @param ... Passed to [mtkcompile()].
#' @return A list with `lin_fun` (class `LinearizationFunction`) and `sys`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); u <- Variable("u", t, input = TRUE)
#' sys <- System(list(D(x) ~ -2 * x + u), t, name = "s")
#' lf <- linearization_function(sys, list(u), list(x), op = list(x ~ 0, u ~ 0))
#' blocks <- lf$lin_fun(c(0), c(0), 0)
#' blocks$f_x
#' @export
linearization_function <- function(sys, inputs, outputs, op = list(), t = 0, simplify = FALSE,
                                   initialize = TRUE, guesses = list(), p = list(),
                                   zero_dummy_der = FALSE, ignore_system_initial_conditions = FALSE,
                                   loop_openings = list(), loop_opening_params = list(),
                                   loop_opening_drivers = list(), system_modifier = identity,
                                   verbose = TRUE, ...) {
  if (inherits(op, "LinearizationOpPoint")) {
    t <- op$t[1L]
    op <- .build_op_from_solution(op, t)
    ignore_system_initial_conditions <- TRUE
  }
  if (!length(op) && isTRUE(verbose))
    warning("An empty operating point was passed to `linearization_function`. An operating point ",
            "containing the variables that will be changed in `linearize` should be provided. ",
            "Disable this warning with verbose = FALSE.", call. = FALSE)
  if (.is_ap_like(inputs)) {
    tr <- linearization_ap_transform(sys, inputs, outputs, loop_openings)
    sys <- tr$sys
    inputs <- tr$inputs
    outputs <- tr$outputs
    loop_opening_params <- c(loop_opening_params, tr$loop_opening_params)
    loop_opening_drivers <- c(loop_opening_drivers, tr$loop_opening_drivers)
  } else if (length(loop_openings)) {
    lo <- .handle_loop_openings(sys, loop_openings)
    sys <- lo$sys
    loop_opening_params <- c(loop_opening_params, lo$params)
    loop_opening_drivers <- c(loop_opening_drivers, lo$drivers)
  }
  inputs <- .flatten_vars(inputs)
  outputs <- .flatten_vars(outputs)
  sys <- system_modifier(sys)
  .check_io_found(sys, inputs, "inputs")
  .check_io_found(sys, outputs, "outputs")
  ssys <- mtkcompile(sys, inputs = inputs, outputs = outputs, simplify = simplify, ...)
  if (ignore_system_initial_conditions) {
    ics <- .f(ssys, "initial_conditions") %||% symmap()
    keep <- vapply(unclass(ics), function(el) isTRUE(.meta_get(el$var, "analysis_variable")), TRUE)
    ics <- unclass(ics)[keep]
    class(ics) <- "mtk_symmap"
    ssys <- .update_system(ssys, initial_conditions = ics)
  }
  allv <- c(unknowns(ssys), .parameters_of(ssys, TRUE), observables(ssys))
  op <- .resolve_named(as_symmap(op), allv)
  hold <- .hold_opened_signals(op, loop_opening_params, loop_opening_drivers, initial_conditions(ssys))
  op <- hold$op
  if (zero_dummy_der) {
    dummy <- Filter(function(v) !(.key(v) %in% .keys(unknowns(sys))), unknowns(ssys))
    ics <- .f(ssys, "initial_conditions") %||% symmap()
    for (v in dummy) ics <- smap_set(ics, v, 0)
    ssys <- .update_system(ssys, initial_conditions = ics)
  }
  st <- .lin_state(ssys, inputs, outputs, hold$held, hold$drivers)
  st$initialize <- isTRUE(initialize)
  st$guesses <- smap_merge(guesses(ssys), .resolve_named(as_symmap(guesses), allv))
  vals <- .process_op(ssys, smap_merge(op, .resolve_named(as_symmap(p), allv)), require_p = FALSE)
  st$u0 <- vals$u0
  st$p0 <- .lin_fill_params(st, vals$p, vals$u0)
  st$t0 <- t
  lin_fun <- structure(function(u, p, t) .linfun_eval(st, u, p, t),
                       class = c("LinearizationFunction", "function"))
  attr(lin_fun, "state") <- st
  list(lin_fun = lin_fun, sys = ssys)
}

# IONotFoundError: inputs and outputs must be variables of the flattened
# system, named without the root system's name.
.check_io_found <- function(sys, vars, what) {
  if (!length(vars)) return(invisible())
  known <- .keys(c(unknowns(sys), .parameters_of(sys, TRUE), observables(sys)))
  miss <- Filter(function(v) !all(.keys(Filter(function(z) !(.var_kind(z) %in% "iv"), .collect_vars(v))) %in% known), vars)
  if (!length(miss)) return(invisible())
  nms <- vapply(miss, format, "")
  hint <- if (any(startsWith(nms, paste0(nameof(sys), NAMESPACE_SEPARATOR))))
    paste0(" Some are namespaced with the name of the system `", nameof(sys), "`; take them from ",
           "toggle_namespacing(sys, FALSE) or complete(sys) instead of sys.") else ""
  stop("The following ", what, " provided to `mtkcompile` were not found in the system: ",
       paste(nms, collapse = ", "), ".", hint, call. = FALSE)
}

#' @export
print.LinearizationFunction <- function(x, ...) {
  st <- attr(x, "state")
  cat("LinearizationFunction of system ", nameof(st$sys), " with ", length(st$inputs), " input(s) and ",
      length(st$outputs), " output(s)\n", sep = "")
  invisible(x)
}

# _hold_opened_signals: opened signals without a value in op or an initial
# condition are solved for with the equation opened ~ driver.
.hold_opened_signals <- function(op, params, drivers, ics) {
  held <- list(); hd <- list()
  for (i in seq_along(params)) {
    q <- params[[i]]
    if (smap_has(op, q) || smap_has(ics, q)) next
    op <- smap_set(op, q, NA)
    held[[length(held) + 1L]] <- q
    hd[[length(hd) + 1L]] <- drivers[[i]]
  }
  list(op = op, held = held, drivers = hd)
}

# Everything lin_fun needs: index sets and compiled Jacobian closures.
.lin_state <- function(ssys, inputs, outputs, held, drivers) {
  st <- new.env(parent = emptyenv())
  st$sys <- ssys
  eqs <- equations(ssys)
  isd <- vapply(eqs, function(e) .is_derivative(e$lhs), TRUE)
  st$diff_idxs <- which(isd)
  st$alge_idxs <- which(!isd)
  us <- unknowns(ssys)
  pp <- .problem_params(ssys)
  st$params <- pp
  st$inputs <- inputs
  st$outputs <- outputs
  st$input_idx <- match(.keys(inputs), .keys(pp))
  if (anyNA(st$input_idx))
    stop("input(s) ", paste(vapply(inputs[is.na(st$input_idx)], format, ""), collapse = ", "),
         " are not parameters of the simplified system", call. = FALSE)
  st$held <- held
  st$held_idx <- match(.keys(held), .keys(pp))
  args <- .codegen_args(ssys)
  obs <- .codegen_observed(ssys)
  st$F <- .build_closure(.ode_rhs_exprs(ssys, substitute_observed = FALSE), args, observed = obs, out = "vector")
  st$driver_fn <- .build_closure(drivers, args, observed = obs, out = "vector")
  rule <- if (length(obs)) .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs")) else NULL
  full <- .ode_rhs_exprs(ssys, substitute_observed = TRUE)
  hfull <- lapply(outputs, function(o) if (is.null(rule)) o else .subst_fixpoint(o, rule))
  st$h <- .build_closure(hfull, args, out = "vector")
  allp <- .parameters_of(ssys, TRUE)
  iv <- if (is.null(.f(ssys, "iv"))) list() else list(.f(ssys, "iv"))
  others_u <- c(allp, iv)
  others_in <- c(us, Filter(function(q) !(.key(q) %in% .keys(inputs)), allp), iv)
  st$J_xz <- .lin_mat_fn(.partials(full, us, others_u), args)
  st$J_u <- .lin_mat_fn(.partials(full, inputs, others_in), args)
  st$H_xz <- .lin_mat_fn(.partials(hfull, us, others_u), args)
  st$H_u <- .lin_mat_fn(.partials(hfull, inputs, others_in), args)
  st
}

.lin_mat_fn <- function(J, args, observed = list()) {
  nr <- nrow(J); nc <- ncol(J)
  if (nr == 0L || nc == 0L) return(function(u, p, t) matrix(0, nr, nc))
  .build_closure(J, args, observed = observed, out = "matrix")
}

# Parameters without a value: inputs default to zero, Initial(x) to the value
# of x, held signals are solved for; anything else is an error.
.lin_fill_params <- function(st, p, u0) {
  pp <- st$params
  for (i in which(is.na(p))) {
    q <- pp[[i]]
    if (i %in% st$input_idx) p[i] <- 0
    else if (i %in% st$held_idx) p[i] <- smap_get(st$guesses, q) %||% 0
    else if (isinitial(q)) {
      j <- match(.key(.args(q)[[1L]]), .keys(unknowns(st$sys)))
      p[i] <- if (is.na(j) || is.na(u0[j])) 0 else u0[j]
    } else stop("no value for parameter ", format(q), "; give it in the operating point", call. = FALSE)
  }
  p
}

# Consistent algebraic unknowns and held signals, with the differential
# unknowns fixed (OverrideInit with algebraic_only = true).
.lin_initialize <- function(st, u, p, t) {
  us <- unknowns(st$sys)
  for (i in which(is.na(u))) {
    g <- smap_get(st$guesses, us[[i]])
    u[i] <- if (is.numeric(g)) g else 0
  }
  a <- st$alge_idxs
  h <- st$held_idx
  if (!st$initialize || (!length(a) && !length(h))) return(list(u = u, p = p))
  res <- function(w) {
    u[a] <- w[seq_along(a)]
    p[h] <- w[length(a) + seq_along(h)]
    c(st$F(u, p, t)[a], if (length(h)) p[h] - st$driver_fn(u, p, t))
  }
  w0 <- c(u[a], p[h])
  r0 <- res(w0)
  if (all(is.finite(r0)) && max(abs(r0), 0) < 1e-10) return(list(u = u, p = p))
  sol <- .newton_solve(res, w0)
  if (!sol$converged)
    stop("Initialization failed: no consistent values for the algebraic unknowns at the operating point", call. = FALSE)
  u[a] <- sol$x[seq_along(a)]
  p[h] <- sol$x[length(a) + seq_along(h)]
  list(u = u, p = p)
}

.linfun_eval <- function(st, u, p, t) {
  if (length(u) != length(unknowns(st$sys)))
    stop("Number of unknown variables (", length(unknowns(st$sys)),
         ") does not match the number of input unknowns (", length(u), ")", call. = FALSE)
  u <- as.numeric(u); p <- as.numeric(p)
  init <- .lin_initialize(st, u, p, t)
  u <- init$u; p <- init$p
  fg_xz <- st$J_xz(u, p, t)
  fg_u <- st$J_u(u, p, t)
  h_xz <- st$H_xz(u, p, t)
  h_u <- st$H_u(u, p, t)
  d <- st$diff_idxs; a <- st$alge_idxs
  list(f_x = fg_xz[d, d, drop = FALSE], f_z = fg_xz[d, a, drop = FALSE],
       g_x = fg_xz[a, d, drop = FALSE], g_z = fg_xz[a, a, drop = FALSE],
       f_u = fg_u[d, , drop = FALSE], g_u = fg_u[a, , drop = FALSE],
       h_x = h_xz[, d, drop = FALSE], h_z = h_xz[, a, drop = FALSE], h_u = h_u,
       x = stats::setNames(u, vapply(unknowns(st$sys), format, "")),
       p = stats::setNames(p, vapply(st$params, format, "")), t = t)
}

# State space matrices from the Jacobian blocks (solve(::LinearizationProblem)).
# With eliminate = TRUE the algebraic unknowns are eliminated through
# z = -g_z^-1 (g_x x + g_u u); otherwise the matrices of ModelingToolkit are
# returned, whose state includes z.
.lin_matrices <- function(b, inputs, allow_input_derivatives = FALSE, eliminate = TRUE) {
  nx <- nrow(b$f_x); nu <- ncol(b$f_u); nz <- ncol(b$f_z); ny <- nrow(b$h_x)
  if (nz == 0L) return(list(A = b$f_x, B = b$f_u, C = b$h_x, D = b$h_u))
  gz <- tryCatch(qr(b$g_z), error = function(e) NULL)
  if (is.null(gz) || nrow(b$g_z) != nz || gz$rank < nz)
    stop("g_z not invertible, this indicates that the DAE is of index > 1.", call. = FALSE)
  lsolve <- function(B) if (ncol(B)) solve(b$g_z, B) else matrix(0, nz, 0L)
  gzgx <- -lsolve(b$g_x)
  gzgu <- -lsolve(b$g_u)
  if (eliminate) {
    return(list(A = b$f_x + b$f_z %*% gzgx, B = b$f_u + b$f_z %*% gzgu,
                C = b$h_x + b$h_z %*% gzgx, D = b$h_u + b$h_z %*% gzgu))
  }
  A <- rbind(cbind(b$f_x, b$f_z), cbind(gzgx %*% b$f_x, gzgx %*% b$f_z))
  B <- rbind(b$f_u, gzgx %*% b$f_u)
  C <- cbind(b$h_x, b$h_z)
  D <- b$h_u
  if (any(gzgu != 0)) {
    if (!allow_input_derivatives) {
      der <- which(colSums(gzgu != 0) > 0)
      stop("Input derivatives appeared in expressions (-g_z\\g_u != 0), the following inputs appeared ",
           "differentiated: [", paste(vapply(inputs[der], format, ""), collapse = ", "), "]. Call with ",
           "allow_input_derivatives = TRUE to allow this and have the returned B matrix be of double ",
           "width (", 2L * nu, "), where the last ", nu, " inputs are the derivatives of the first ", nu,
           " inputs.", call. = FALSE)
    }
    B <- cbind(B, rbind(matrix(0, nx, nu), gzgu))
    D <- cbind(D, matrix(0, ny, nu))
  }
  list(A = A, B = B, C = C, D = D)
}

#' Linearize a system
#'
#' `linearize(sys, inputs, outputs)` linearizes `sys` around the operating
#' point `op` at time `t` and returns the matrices of
#' \deqn{x' = A x + B u, \quad y = C x + D u.}
#' `inputs` and `outputs` are lists of variables or analysis points. The
#' state `x` is all unknowns of the simplified system `sys` (returned in
#' the result), as in ModelingToolkit's `linearize`: for a system with
#' algebraic unknowns `z` the state is `(x, z)` and
#' \eqn{A = [f_x, f_z; G f_x, G f_z]} with \eqn{G = -g_z^{-1} g_x}; input
#' derivatives then need `allow_input_derivatives = TRUE`. With
#' `eliminate_algebraic = TRUE` the algebraic unknowns are eliminated
#' through \eqn{z = -g_z^{-1}(g_x x + g_u u)}, so that
#' \eqn{A = f_x - f_z g_z^{-1} g_x} and the state is the differential
#' unknowns only.
#'
#' `linearize(ssys, lin_fun)` linearizes with a function from
#' [linearization_function()], which avoids compiling again for new
#' operating points.
#'
#' @param sys A `System`, or the simplified system for the second form.
#' @param inputs,outputs Lists of variables or analysis points; for the
#'   second form, `inputs` is the `LinearizationFunction`.
#' @param op Operating point, a variable map or a
#'   [LinearizationOpPoint()].
#' @param t Time of the linearization.
#' @param allow_input_derivatives Allow derivatives of inputs, see above.
#' @param eliminate_algebraic Eliminate algebraic unknowns from the state.
#' @param zero_dummy_der Set dummy derivatives to zero at the operating
#'   point.
#' @param p Further parameter values.
#' @param ... Passed to [linearization_function()].
#' @return A list with matrices `A`, `B`, `C`, `D`, the simplified system
#'   `sys`, the names of the `states`, and `extras`, the operating point
#'   `list(x, p, t)` after initialization. With a vector of times in a
#'   `LinearizationOpPoint`, a list of such results.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); v <- Variable("v", t)
#' F <- Variable("F", t, input = TRUE); y <- Variable("y", t, output = TRUE)
#' m <- Parameter("m", default = 2); k <- Parameter("k", default = 4)
#' c <- Parameter("c", default = 1)
#' msd <- System(list(D(x) ~ v, D(v) ~ (F - k * x - c * v) / m, y ~ x), t, name = "msd")
#' lin <- linearize(msd, list(F), list(y), op = list(x ~ 0, v ~ 0, F ~ 0))
#' lin$A; lin$B; lin$C; lin$D
#' @export
linearize <- function(sys, inputs, outputs = NULL, op = list(), t = 0, allow_input_derivatives = FALSE,
                      eliminate_algebraic = FALSE, zero_dummy_der = FALSE, p = list(), ...) {
  multi <- inherits(op, "LinearizationOpPoint") && length(op$t) > 1L
  if (inherits(inputs, "LinearizationFunction")) {
    if (!multi)
      return(.linearize_with(inputs, op = op, t = t, allow_input_derivatives = allow_input_derivatives,
                             eliminate_algebraic = eliminate_algebraic, p = p))
    lf <- list(lin_fun = inputs)
  } else if (multi) {
    lf <- linearization_function(sys, inputs, outputs, op = LinearizationOpPoint(op$sol, op$t[1L], op$op),
                                 zero_dummy_der = zero_dummy_der, p = p, ...)
  }
  if (multi) {
    return(lapply(op$t, function(ti) {
      .linearize_with(lf$lin_fun, op = .build_op_from_solution(op, ti), t = ti,
                      allow_input_derivatives = allow_input_derivatives,
                      eliminate_algebraic = eliminate_algebraic, p = p)
    }))
  }
  lf <- linearization_function(sys, inputs, outputs, op = op, t = t, zero_dummy_der = zero_dummy_der, p = p, ...)
  if (inherits(op, "LinearizationOpPoint")) {
    t <- op$t
    op <- .build_op_from_solution(op, t)
  }
  .linearize_with(lf$lin_fun, op = op, t = t, allow_input_derivatives = allow_input_derivatives,
                  eliminate_algebraic = eliminate_algebraic, p = p)
}

# linearize(ssys, lin_fun; op, t): apply op to the stored operating point
# and solve the linearization problem.
.linearize_with <- function(lin_fun, op = list(), t = 0, allow_input_derivatives = FALSE,
                            eliminate_algebraic = FALSE, p = list()) {
  if (inherits(op, "LinearizationOpPoint")) {
    t <- op$t[1L]
    op <- .build_op_from_solution(op, t)
  }
  prob <- LinearizationProblem(lin_fun, t = t)
  prob <- .lin_apply_op(prob, op, check_held = TRUE)
  prob <- .lin_apply_op(prob, p)
  solve(prob, allow_input_derivatives = allow_input_derivatives, eliminate_algebraic = eliminate_algebraic)
}

.lin_apply_op <- function(prob, op, check_held = FALSE) {
  if (!length(op)) return(prob)
  st <- attr(prob$f, "state")
  us <- unknowns(st$sys)
  pp <- st$params
  op <- .resolve_named(as_symmap(op), c(us, .parameters_of(st$sys, TRUE), observables(st$sys)))
  if (check_held && length(st$held)) {
    given <- Filter(function(q) {
      v <- smap_get(op, q)
      !is.null(v) && !(length(v) == 1L && !is_sym(v) && is.na(v))
    }, st$held)
    if (length(given))
      stop("The operating point provides values for the signal(s) ",
           paste(vapply(given, format, ""), collapse = ", "), ", which are disconnected by loop ",
           "openings and held at the value of the signals that drive them. Pass the values in the ",
           "operating point given to `linearization_function` (or to the analysis function that ",
           "builds it) instead.", call. = FALSE)
  }
  cur <- symmap()
  for (i in seq_along(us)) cur <- smap_set(cur, us[[i]], prob$u0[i])
  for (i in seq_along(pp)) cur <- smap_set(cur, pp[[i]], prob$p[i])
  given <- Filter(function(el) !(length(el$value) == 1L && !is_sym(el$value) && is.na(el$value)), unclass(op))
  vals <- .evaluate_varmap(smap_merge(cur, structure(given, class = "mtk_symmap")))
  for (el in given) {
    k <- .key(el$var)
    v <- vals[[k]]
    if (is.null(v) || !is.numeric(v) || is.na(v)) next
    i <- match(k, .keys(us))
    if (!is.na(i)) prob$u0[i] <- v
    j <- match(k, .keys(pp))
    if (!is.na(j)) prob$p[j] <- v
  }
  prob
}

#' Linearization problems
#'
#' `LinearizationProblem(sys, inputs, outputs, t = 0, ...)` builds a problem
#' that is linearized by `solve(prob)` at its current operating point.
#' Values of unknowns and parameters are changed with `prob[x] <- value` and
#' the time with `prob$t <- value` between calls, which avoids compiling the
#' system again. `solve()` returns a list with the matrices `A`, `B`, `C`,
#' `D` and `extras`.
#'
#' @param sys A `System`, or a `LinearizationFunction`.
#' @param inputs,outputs Lists of variables.
#' @param t Time of the linearization.
#' @param ... Passed to [linearization_function()].
#' @return A `LinearizationProblem`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); u <- Variable("u", t, input = TRUE)
#' a <- Parameter("a", default = 1)
#' sys <- System(list(D(x) ~ -a * x^2 + u), t, name = "s")
#' prob <- LinearizationProblem(sys, list(u), list(x), op = list(x ~ 1, u ~ 0))
#' solve(prob)$A
#' prob[x] <- 2
#' solve(prob)$A
#' @export
LinearizationProblem <- function(sys, inputs = NULL, outputs = NULL, t = 0, ...) {
  lin_fun <- if (inherits(sys, "LinearizationFunction")) sys else
    linearization_function(sys, inputs, outputs, t = t, ...)$lin_fun
  st <- attr(lin_fun, "state")
  structure(list(f = lin_fun, t = t, u0 = st$u0, p = st$p0), class = "LinearizationProblem")
}

#' @export
print.LinearizationProblem <- function(x, ...) {
  cat("LinearizationProblem at time ", x$t, " of system ", nameof(attr(x$f, "state")$sys), "\n", sep = "")
  invisible(x)
}

#' @export
`[.LinearizationProblem` <- function(x, i) {
  st <- attr(x$f, "state")
  k <- .key(as_sym(i))
  j <- match(k, .keys(unknowns(st$sys)))
  if (!is.na(j)) return(unname(x$u0[j]))
  j <- match(k, .keys(st$params))
  if (!is.na(j)) return(unname(x$p[j]))
  stop(format(i), " is not an unknown or parameter of the linearization problem", call. = FALSE)
}

#' @export
`[<-.LinearizationProblem` <- function(x, i, value) {
  st <- attr(x$f, "state")
  k <- .key(as_sym(i))
  j <- match(k, .keys(unknowns(st$sys)))
  if (!is.na(j)) {
    x$u0[j] <- value
    return(x)
  }
  j <- match(k, .keys(st$params))
  if (!is.na(j)) {
    x$p[j] <- value
    return(x)
  }
  stop(format(i), " is not an unknown or parameter of the linearization problem", call. = FALSE)
}

#' @rdname LinearizationProblem
#' @param a A `LinearizationProblem`.
#' @param b Unused.
#' @param allow_input_derivatives,eliminate_algebraic See [linearize()].
#' @export
solve.LinearizationProblem <- function(a, b, allow_input_derivatives = FALSE, eliminate_algebraic = FALSE, ...) {
  st <- attr(a$f, "state")
  blocks <- a$f(a$u0, a$p, a$t)
  mats <- .lin_matrices(blocks, st$inputs, allow_input_derivatives, eliminate_algebraic)
  us <- vapply(unknowns(st$sys), format, "")
  states <- if (eliminate_algebraic) us[st$diff_idxs] else us[c(st$diff_idxs, st$alge_idxs)]
  c(mats, list(sys = st$sys, states = states, extras = list(x = blocks$x, p = blocks$p, t = blocks$t)))
}

#' Symbolic linearization
#'
#' `linearize_symbolic(sys, inputs, outputs)` returns the symbolic matrices
#' `A`, `B`, `C`, `D` (symcompR matrices) of the linearization of `sys`
#' together with the Jacobian blocks `f_x`, ..., `h_u`, and the simplified
#' system. Algebraic unknowns are eliminated as in [linearize()] unless
#' `eliminate_algebraic = FALSE`.
#'
#' @inheritParams linearize
#' @param simplify Simplify the resulting expressions.
#' @return A list with the matrices, the blocks and `sys`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); u <- Variable("u", t, input = TRUE)
#' y <- Variable("y", t, output = TRUE); k <- Parameter("k")
#' sys <- System(list(D(x) ~ -k * x + u, y ~ x), t, name = "s")
#' linearize_symbolic(sys, list(u), list(y))$A
#' @export
linearize_symbolic <- function(sys, inputs, outputs, simplify = FALSE, allow_input_derivatives = FALSE,
                               eliminate_algebraic = FALSE, ...) {
  inputs <- .flatten_vars(inputs); outputs <- .flatten_vars(outputs)
  ssys <- mtkcompile(sys, inputs = inputs, outputs = outputs, simplify = simplify, ...)
  eqs <- equations(ssys)
  isd <- vapply(eqs, function(e) .is_derivative(e$lhs), TRUE)
  d <- which(isd); a <- which(!isd)
  us <- unknowns(ssys)
  obs <- .codegen_observed(ssys)
  rule <- if (length(obs)) .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs")) else NULL
  full <- .ode_rhs_exprs(ssys, substitute_observed = TRUE)
  hfull <- lapply(outputs, function(o) if (is.null(rule)) o else .subst_fixpoint(o, rule))
  allp <- .parameters_of(ssys, TRUE)
  iv <- list(.f(ssys, "iv"))
  others_in <- c(us, Filter(function(q) !(.key(q) %in% .keys(inputs)), allp), iv)
  fg_xz <- .partials(full, us, c(allp, iv))
  fg_u <- .partials(full, inputs, others_in)
  h_xz <- .partials(hfull, us, c(allp, iv))
  h_u <- .partials(hfull, inputs, others_in)
  sub <- function(M, r, c) .sym_mat(M[r, c, drop = FALSE])
  bl <- list(f_x = sub(fg_xz, d, d), f_z = sub(fg_xz, d, a), g_x = sub(fg_xz, a, d), g_z = sub(fg_xz, a, a),
             f_u = sub(fg_u, d, seq_along(inputs)), g_u = sub(fg_u, a, seq_along(inputs)),
             h_x = sub(h_xz, seq_along(outputs), d), h_z = sub(h_xz, seq_along(outputs), a),
             h_u = sub(h_u, seq_along(outputs), seq_along(inputs)))
  if (!length(a)) {
    mats <- list(A = bl$f_x, B = bl$f_u, C = bl$h_x, D = bl$h_u)
  } else {
    gzgx <- .sym_neg(.sym_solve(bl$g_z, bl$g_x))
    gzgu <- .sym_neg(.sym_solve(bl$g_z, bl$g_u))
    if (eliminate_algebraic) {
      mats <- list(A = .sym_add(bl$f_x, .sym_mul(bl$f_z, gzgx)), B = .sym_add(bl$f_u, .sym_mul(bl$f_z, gzgu)),
                   C = .sym_add(bl$h_x, .sym_mul(bl$h_z, gzgx)), D = .sym_add(bl$h_u, .sym_mul(bl$h_z, gzgu)))
    } else {
      if (!.sym_iszero(gzgu) && !allow_input_derivatives)
        stop("Input derivatives appeared in expressions (-g_z\\g_u != 0); call linearize_symbolic with ",
             "allow_input_derivatives = TRUE.", call. = FALSE)
      mats <- list(A = .sym_rbind(.sym_cbind(bl$f_x, bl$f_z), .sym_cbind(.sym_mul(gzgx, bl$f_x), .sym_mul(gzgx, bl$f_z))),
                   B = .sym_rbind(bl$f_u, .sym_mul(gzgx, bl$f_u)), C = .sym_cbind(bl$h_x, bl$h_z), D = bl$h_u)
      if (!.sym_iszero(gzgu)) {
        mats$B <- .sym_cbind(mats$B, .sym_rbind(.sym_zeros(length(d), length(inputs)), gzgu))
        mats$D <- .sym_cbind(mats$D, .sym_zeros(length(outputs), length(inputs)))
      }
    }
  }
  if (simplify) mats <- lapply(mats, function(M) { M[] <- lapply(M, symcompR::simplify); M })
  c(mats, bl, list(sys = ssys))
}

# Symbolic matrices are kept as list matrices of sym entries.
.sym_mat <- function(M) M
.sym_neg <- function(A) { A[] <- lapply(A, function(e) -e); A }
.sym_zeros <- function(r, c) .list_matrix(r, c, function(i, j) .S0())
.sym_add <- function(A, B) .list_matrix(nrow(A), ncol(A), function(i, j) A[[i, j]] + B[[i, j]])
.sym_mul <- function(A, B) {
  .list_matrix(nrow(A), ncol(B), function(i, j) {
    if (!ncol(A)) return(.S0())
    Reduce(`+`, lapply(seq_len(ncol(A)), function(k) A[[i, k]] * B[[k, j]]))
  })
}
.sym_solve <- function(A, B) {
  if (!ncol(B)) return(.sym_zeros(nrow(A), 0L))
  X <- solve(.as_symmatrix(A), .as_symmatrix(B))
  .list_matrix(nrow(A), ncol(B), function(i, j) X[i, j])
}
.sym_cbind <- function(A, B) {
  out <- matrix(list(), nrow(A), ncol(A) + ncol(B))
  if (ncol(A)) out[, seq_len(ncol(A))] <- A
  if (ncol(B)) out[, ncol(A) + seq_len(ncol(B))] <- B
  out
}
.sym_rbind <- function(A, B) {
  out <- matrix(list(), nrow(A) + nrow(B), ncol(A))
  if (nrow(A)) out[seq_len(nrow(A)), ] <- A
  if (nrow(B)) out[nrow(A) + seq_len(nrow(B)), ] <- B
  out
}
.sym_iszero <- function(M) all(vapply(M, .is_zero, TRUE))

#' Transform linear state space models
#'
#' `similarity_transform(sys, T)` changes the state coordinates of a linear
#' model `list(A, B, C, D)` (as returned by [linearize()]) with `T x2 = x`:
#' \eqn{A_2 = T^{-1} A T}, \eqn{B_2 = T^{-1} B}, \eqn{C_2 = C T}. With
#' `unitary = TRUE` the transpose of `T` is used instead of its inverse.
#' `reorder_unknowns(sys, old, new)` permutes the states from the order
#' `old` (the `states` of a [linearize()] result) to `new`.
#'
#' @param sys A list with matrices `A`, `B`, `C`, `D`.
#' @param T The transformation matrix.
#' @param unitary Whether `T` is unitary.
#' @param old,new The current and desired orders of the states: lists of
#'   variables or character vectors. For a result of [linearize()] the
#'   current order is its `states`.
#' @return A list with matrices `A`, `B`, `C`, `D`.
#' @examples
#' lsys <- list(A = matrix(c(-1, 0, 0, -2), 2), B = matrix(c(1, 0), 2),
#'   C = matrix(c(0, 1), 1), D = matrix(0, 1, 1))
#' similarity_transform(lsys, matrix(c(0, 1, 1, 0), 2), unitary = TRUE)$A
#' reorder_unknowns(lsys, c("a", "b"), c("b", "a"))$A
#' @export
similarity_transform <- function(sys, T, unitary = FALSE) {
  if (unitary) {
    A <- t(T) %*% sys$A %*% T
    B <- t(T) %*% sys$B
  } else {
    A <- solve(T, sys$A %*% T)
    B <- solve(T, sys$B)
  }
  list(A = A, B = B, C = sys$C %*% T, D = sys$D)
}

#' @rdname similarity_transform
#' @export
reorder_unknowns <- function(sys, old, new) {
  nm <- function(x) if (is.character(x)) x else vapply(.flatten_vars(x), format, "")
  old <- nm(old); new <- nm(new)
  nx <- length(old)
  if (length(new) != nx) stop("old and new must have the same length", call. = FALSE)
  perm <- match(new, old)
  if (anyNA(perm)) stop("new contains states not in old", call. = FALSE)
  if (!is.unsorted(perm)) return(sys)
  P <- matrix(0, nx, nx)
  for (i in seq_len(nx)) P[perm[i], i] <- 1
  similarity_transform(sys, P, unitary = TRUE)
}
