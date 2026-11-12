# The Modelica homotopy operator, the limited operator and polynomial
# homotopy continuation, after ModelingToolkitBase/src/systems/
# homotopy_operator.jl (homotopy, has_homotopy, strip_homotopy, HomotopyCtx,
# homotopy_enabled, lower_homotopy, generate_homotopy_residual),
# problems/homotopyproblem.jl (HomotopyProblem), systems/limited_operator.jl
# (limited, limitnew, limitold, lower_limited, the postcondition corrector),
# systems/nonlinear/homotopy_continuation.jl (HomotopyContinuationProblem and
# the polynomial checks), NonlinearSolveBase/src/homotopy_sweep.jl
# (HomotopySweep) and NonlinearSolveHomotopyContinuation (the solve
# interface). The path tracker for polynomial systems is a native total
# degree homotopy, since HomotopyContinuation.jl has no R counterpart.

# ---------------------------------------------------------------------------
# homotopy(actual, simplified)

.nl_register_operator <- function(name) {
  if (!exists(name, envir = .mtk$funs, inherits = FALSE))
    register_symbolic(name, function(actual, other) actual,
                      derivatives = list(function(actual, other) 1, function(actual, other) 0))
  Function(name)
}

#' The Modelica homotopy operator
#'
#' `homotopy(actual, simplified)` declares that `simplified` is an easy to
#' solve approximation of `actual`. Wherever the operator is evaluated
#' numerically it is `actual`. A [HomotopyProblem()] built from a system
#' that contains the operator replaces every node by the blend
#' `(1 - lambda) * simplified + lambda * actual` and is solved by sweeping
#' `lambda` from 0 to 1, which walks from the easy solution to the true one.
#'
#' `strip_homotopy(sys)` replaces every node by its `actual` branch and
#' records `FALSE` under the metadata key `HomotopyCtx`, which
#' `homotopy_enabled(sys)` reads; in ModelingToolkit this is what
#' `mtkcompile(sys; homotopy = false)` does.
#'
#' @param actual The real expression.
#' @param simplified An approximation whose solution is easy to reach.
#' @param sys A `System`.
#' @return A `sym` expression; a `System` for `strip_homotopy()`; a logical
#'   for `homotopy_enabled()`.
#' @examples
#' y <- Variable("y")
#' sys <- complete(System(list(0 ~ homotopy(atan(y - 3), y)), name = "h"))
#' prob <- HomotopyProblem(sys, list(y ~ 12))
#' solve(prob)[y]
#' homotopy_enabled(sys)
#' equations(strip_homotopy(sys))
#' @export
homotopy <- function(actual, simplified) {
  if (is.numeric(actual) && is.numeric(simplified)) return(actual)
  do.call(.nl_register_operator("homotopy"), list(as_sym(actual), as_sym(simplified)))
}

#' @rdname homotopy
#' @format `HomotopyCtx` is the metadata key, a string.
#' @export
HomotopyCtx <- "HomotopyCtx"

#' @rdname homotopy
#' @export
homotopy_enabled <- function(sys) !identical(.f(sys, "metadata")[[HomotopyCtx]], FALSE)

# Rebuild an expression bottom up, replacing every application of the
# function `fname` by replace(args...).
.nl_rewrite_call <- function(e, fname, replace) {
  cache <- new.env(parent = emptyenv())
  walk <- function(z) {
    if (!is_sym(z)) return(z)
    a <- .args(z)
    if (!length(a)) return(z)
    k <- .key(z)
    hit <- cache[[k]]
    if (!is.null(hit)) return(hit)
    args <- lapply(a, walk)
    out <- if (.type(z) == "Function" && identical(.name(z), fname)) do.call(replace, args)
           else if (identical(.keys(args), .keys(a))) z else .rebuild(z, args)
    assign(k, out, envir = cache)
    out
  }
  if (inherits(e, "Equation")) return(Equation(walk(e$lhs), walk(e$rhs)))
  if (is_sym(e)) walk(e) else lapply(e, .nl_rewrite_call, fname = fname, replace = replace)
}

.nl_strip_map <- function(m, rew) {
  for (k in names(m)) if (is_sym(m[[k]]$value)) m[[k]]$value <- rew(m[[k]]$value)
  m
}

#' @rdname homotopy
#' @export
strip_homotopy <- function(sys) {
  rew <- function(x) .nl_rewrite_call(x, "homotopy", function(a, s) a)
  md <- .f(sys, "metadata")
  md[[HomotopyCtx]] <- FALSE
  out <- .update_system(sys,
    eqs = lapply(.f(sys, "eqs"), function(e) if (inherits(e, "Equation")) rew(e) else e),
    observed = lapply(.f(sys, "observed"), rew),
    initialization_eqs = lapply(.f(sys, "initialization_eqs"), rew),
    costs = lapply(.f(sys, "costs"), rew),
    constraints = lapply(.f(sys, "constraints"), function(cn) if (inherits(cn, "Equation")) rew(cn) else cn),
    bindings = .nl_strip_map(.f(sys, "bindings"), rew),
    initial_conditions = .nl_strip_map(.f(sys, "initial_conditions"), rew),
    guesses = .nl_strip_map(.f(sys, "guesses"), rew),
    systems = lapply(.f(sys, "systems"), strip_homotopy),
    metadata = md)
  if (!is.null(.f(out, "noise_eqs")) && is.list(.f(out, "noise_eqs")))
    out <- .update_system(out, noise_eqs = lapply(.f(out, "noise_eqs"), rew))
  out
}

# The continuation parameter is a fixed reserved symbol, never a parameter
# of the system (HOMOTOPY_LAMBDA).
.nl_homotopy_lambda <- function() Symbol("__homotopy_\u03BB\u2098\u209C\u2096")

# generate_homotopy_residual: the residual f(u, p, lambda) with every
# homotopy(actual, simplified) node replaced by the blend, and its Jacobian.
.nl_homotopy_residual <- function(sys) {
  lam <- .nl_homotopy_lambda()
  resid <- .nl_rewrite_call(.ode_rhs_exprs(sys), "homotopy", function(a, s) lam * a + (1 - lam) * s)
  args <- c(.codegen_args(sys, with_t = FALSE), list(list(name = "lambda", vars = lam)))
  list(f = .build_closure(resid, args, out = "vector"),
       jac = .build_closure(.partials(resid, unknowns(sys), c(.nl_hold(sys), list(lam))), args, out = "matrix"))
}

#' Homotopy problems
#'
#' `HomotopyProblem(sys, op)` builds a continuation problem from a system
#' whose equations contain [homotopy()] nodes: the residual is `f(u, p,
#' lambda)` with every node replaced by `(1 - lambda) * simplified + lambda *
#' actual`. `solve(prob, HomotopySweep())` solves it at `lambdaspan[1]` and
#' follows the solution to `lambdaspan[2]`, predicting each step by secant
#' extrapolation and correcting with a nonlinear solver; a failed step halves
#' the increment. `AbstractNonlinearProblem(sys, op)` returns a
#' `HomotopyProblem` when the system contains `homotopy` nodes and a
#' [NonlinearProblem()] otherwise.
#'
#' @param sys A complete `System`.
#' @param op Operating point: values of unknowns and parameters.
#' @param lambdaspan The range of the continuation parameter.
#' @param ... Solver options stored with the problem.
#' @param guesses,missing_guess_value As in [NonlinearProblem()].
#' @param inner The nonlinear algorithm used at each step.
#' @param initial_step_factor First increment as a fraction of the span.
#' @param min_dlambda Smallest increment before giving up.
#' @param max_step_factor Largest increment as a fraction of the span.
#' @param expand_factor,expand_threshold Growth of the increment after
#'   `expand_threshold` accepted steps in a row.
#' @param maxsteps Maximum number of steps.
#' @return A `HomotopyProblem`, or the algorithm object.
#' @examples
#' y <- Variable("y")
#' sys <- complete(System(list(0 ~ homotopy(atan(y - 3), y)), name = "h"))
#' prob <- HomotopyProblem(sys, list(y ~ 12))
#' sol <- solve(prob, HomotopySweep())
#' sol$u; sol$retcode
#' class(AbstractNonlinearProblem(sys, list(y ~ 12)))[1]
#' @export
HomotopyProblem <- function(sys, op = list(), lambdaspan = c(0, 1), ..., guesses = list(),
                            missing_guess_value = NULL) {
  if (inherits(sys, "System") && !is.null(.f(sys, "iv")) && iscomplete(sys)) sys <- .nl_from_time_dependent(sys)
  .nl_check_system(sys, "HomotopyProblem")
  if (!.nl_has_call(sys, "homotopy"))
    stop("`HomotopyProblem` requires a system containing `homotopy(actual, simplified)` nodes; none ",
         "were found. Use `NonlinearProblem` for systems without `homotopy`.", call. = FALSE)
  if (length(equations(sys)) != length(unknowns(sys)))
    .stopf("`HomotopyProblem` requires as many equations as unknowns; the system has %d equations and %d unknowns.",
           length(equations(sys)), length(unknowns(sys)))
  lam <- .nl_homotopy_lambda()
  for (v in c(unknowns(sys), .parameters_of(sys, TRUE), observables(sys)))
    if (identical(.name(v), .name(lam)))
      stop("the homotopy continuation sentinel `", .name(lam), "` collides with a symbol of the system; ",
           "this name is reserved for internal use.", call. = FALSE)
  hfn <- .nl_homotopy_residual(sys)
  vals <- .nl_values(sys, op, guesses, missing_guess_value)
  fn <- structure(list(f = hfn$f, jac = hfn$jac, sys = sys,
                       observed = function(expr) build_explicit_observed_function(sys, expr)),
                  class = c("NonlinearFunction", "SciMLFunction"))
  prob <- .new_problem("HomotopyProblem", f = fn, u0 = vals$u0, p = vals$p, sys = sys, kwargs = list(...))
  prob$op <- vals$op
  prob$lambdaspan <- as.numeric(lambdaspan)
  prob$guesses <- guesses
  prob$missing_guess_value <- missing_guess_value
  class(prob) <- c("HomotopyProblem", "AbstractNonlinearProblem", "SciMLProblem")
  prob
}

#' @rdname HomotopyProblem
#' @export
AbstractNonlinearProblem <- function(sys, op = list(), ...) {
  s2 <- if (inherits(sys, "System") && !is.null(.f(sys, "iv")) && iscomplete(sys)) .nl_from_time_dependent(sys) else sys
  if (inherits(s2, "System") && .nl_has_call(s2, "homotopy")) HomotopyProblem(s2, op, ...)
  else NonlinearProblem(s2, op, ...)
}

#' @rdname HomotopyProblem
#' @export
HomotopySweep <- function(inner = NULL, initial_step_factor = 0.1, min_dlambda = sqrt(.Machine$double.eps),
                          max_step_factor = 1, expand_factor = 2, expand_threshold = 2L, maxsteps = 10000L) {
  .nl_alg("HomotopySweep", inner = inner, initial_step_factor = initial_step_factor,
          min_dlambda = min_dlambda, max_step_factor = max_step_factor, expand_factor = expand_factor,
          expand_threshold = as.integer(expand_threshold), maxsteps = as.integer(maxsteps))
}

#' @rdname solve_nonlinear
#' @export
solve.HomotopyProblem <- function(a, b, ...) {
  alg <- if (missing(b) || is.null(b)) HomotopySweep() else b
  if (!inherits(alg, "HomotopySweep")) alg <- HomotopySweep(inner = alg)
  opts <- .nl_opts(a, list(...))
  # an initialization problem keeps the continuation residual under f$homotopy
  hf <- a$f$homotopy %||% a$f
  f <- hf$f; jac <- hf$jac; p <- a$p
  at <- function(lam) list(F = function(x) as.numeric(f(x, p, lam)), J = function(x) jac(x, p, lam))
  res <- .nl_sweep(at, a$u0, alg, opts, a$lambdaspan[1L], a$lambdaspan[2L])
  resid <- as.numeric(f(res$u, p, a$lambdaspan[2L]))
  .nl_new_solution("NonlinearSolution", a, alg, res$u, res$retcode, resid = resid,
                   stats = list(nsteps = res$nsteps))
}

# Natural parameter continuation (HomotopySweep): solve at l0 and follow the
# solution to l1 with a secant predictor; a failed step halves the
# increment. `at(lambda)` returns the residual F and Jacobian J there.
.nl_sweep <- function(at, u0, alg, opts, l0, l1) {
  # the continuation algorithms of nonlinear-homotopy-extras.R carry their own driver
  if (is.function(alg$sweep)) return(alg$sweep(at, u0, opts, l0, l1))
  inner <- alg$inner %||% .nl_default_alg(u0, TRUE)
  width <- l1 - l0
  fn <- at(l0)
  res <- .nl_run(inner, fn$F, fn$J, u0, opts, ls = FALSE)
  nsteps <- 0L
  retcode <- res$retcode
  u <- res$u
  if (identical(retcode, "Success") && width != 0) {
    lam <- l0
    dl <- alg$initial_step_factor * width
    uprev <- NULL; lprev <- NULL
    streak <- 0L
    while (abs(l1 - lam) > 0 && nsteps < alg$maxsteps) {
      if (abs(dl) > abs(l1 - lam)) dl <- l1 - lam
      ln <- lam + dl
      # secant predictor along the solution path
      guess <- if (is.null(uprev)) u else u + (u - uprev) * (ln - lam) / (lam - lprev)
      fn <- at(ln)
      res <- .nl_run(inner, fn$F, fn$J, guess, opts, ls = FALSE)
      nsteps <- nsteps + 1L
      if (identical(res$retcode, "Success")) {
        uprev <- u; lprev <- lam
        u <- res$u; lam <- ln
        streak <- streak + 1L
        if (streak >= alg$expand_threshold) {
          dl <- sign(width) * min(abs(dl) * alg$expand_factor, alg$max_step_factor * abs(width))
          streak <- 0L
        }
      } else {
        dl <- dl / 2
        streak <- 0L
        if (abs(dl) < alg$min_dlambda) {
          retcode <- "ConvergenceFailure"
          break
        }
      }
    }
    if (identical(retcode, "Success") && lam != l1) retcode <- "MaxIters"
  }
  list(u = u, retcode = retcode, nsteps = nsteps, iters = nsteps)
}

# ---------------------------------------------------------------------------
# limited(actual, limiter)

.nl_sentinel <- function(name) {
  e <- Symbol(name)
  attr(e, "mtk_meta") <- list(kind = "parameter", default = 0, tunable = FALSE, scope = list(type = "global"))
  e
}

#' Limited quantities of a nonlinear solve
#'
#' `limited(actual, limiter)` declares that the Newton updates of the
#' expression `actual` must be corrected at every iteration, in the manner of
#' SPICE junction voltage limiting. `limiter` is an expression in the
#' placeholders `limitnew` (the proposed value) and `limitold` (the previously
#' accepted value) and the parameters of the system; it must return
#' `limitnew` when `limitnew == limitold`.
#'
#' A [NonlinearProblem()] built from a system with `limited` nodes follows
#' the predictor and corrector formulation of Aadithya, Keiter and Mei: each
#' node becomes an auxiliary unknown `limited_k` with the equation
#' `limited_k ~ actual`, and the limiter is applied to that unknown after
#' every iterate the solver accepts. In ModelingToolkit the lowering happens
#' in `mtkcompile`; here it happens when the problem is built. Wherever the
#' operator is evaluated numerically, as in time dependent problems, it is
#' `actual`.
#'
#' @param actual The expression being limited.
#' @param limiter The correction rule in terms of `limitnew` and `limitold`.
#' @return A `sym` expression.
#' @format `limitnew` and `limitold` are `sym` symbols.
#' @examples
#' v <- Variable("v")
#' sys <- complete(System(list(0 ~ exp(limited(v, Min(limitnew, limitold + 1))) - 5), name = "lim"))
#' prob <- NonlinearProblem(sys, list(v ~ 0))
#' sol <- solve(prob, NewtonRaphson())
#' sol[v]
#' @export
limited <- function(actual, limiter) {
  if (is.numeric(actual) && is.numeric(limiter)) return(actual)
  .nl_register_sentinels()
  do.call(.nl_register_operator("limited"), list(as_sym(actual), as_sym(limiter)))
}

#' @rdname limited
#' @export
limitnew <- .nl_sentinel("__limitnew_\u2098\u209C\u2096")

#' @rdname limited
#' @export
limitold <- .nl_sentinel("__limitold_\u2098\u209C\u2096")

.nl_register_sentinels <- function() {
  for (s in list(limitnew, limitold))
    if (!exists(.name(s), envir = .mtk$vars, inherits = FALSE))
      .register_var(.name(s), "parameter", list(default = 0, tunable = FALSE))
  invisible(NULL)
}

.nl_collect_calls <- function(exprs, fname) {
  seen <- new.env(parent = emptyenv())
  out <- list()
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    k <- .key(z)
    if (!is.null(seen[[k]])) return(invisible())
    assign(k, TRUE, envir = seen)
    if (.type(z) == "Function" && identical(.name(z), fname)) {
      out[[length(out) + 1L]] <<- z
      return(invisible())
    }
    for (a in .args(z)) walk(a)
  }
  for (z in exprs) walk(z)
  out
}

# lower_limited: the predictor and corrector augmentation of a time
# independent system. Returns the augmented system and the postcondition
# function(xnew, xold, p).
.lower_limited <- function(sys) {
  eqs <- equations(sys)
  nodes <- .nl_collect_calls(c(lapply(eqs, `[[`, "lhs"), lapply(eqs, `[[`, "rhs")), "limited")
  if (!length(nodes)) return(list(sys = sys, postcondition = NULL))
  sent <- .keys(list(limitnew, limitold))
  new_vars <- list(); new_eqs <- list(); limiters <- list()
  ics <- .f(sys, "initial_conditions"); gs <- .f(sys, "guesses")
  for (k in seq_along(nodes)) {
    a <- .args(nodes[[k]])
    if (.nl_exprs_have_call(a, "limited"))
      stop("`limited` operators may not be nested; found ", format(nodes[[k]]), ".", call. = FALSE)
    v <- .make_var(paste0("#limited_", k), "unknown", NULL, list(irreducible = TRUE))
    new_vars[[k]] <- v
    new_eqs[[k]] <- Equation(v, a[[1L]])
    limiters[[k]] <- a[[2L]]
    ics <- smap_set(ics, v, a[[1L]])
    gs <- smap_set(gs, v, a[[1L]])
  }
  rule <- .rule(nodes, new_vars)
  lowered <- lapply(eqs, .subst, rule = rule)
  drop <- function(vs) Filter(function(p) !(.key(p) %in% sent), vs)
  nsys <- .update_system(sys, eqs = c(lowered, new_eqs), unknowns = c(drop(.f(sys, "unknowns")), new_vars),
                         ps = drop(.f(sys, "ps")), initial_conditions = ics, guesses = gs,
                         irreducibles = c(.f(sys, "irreducibles"), new_vars))
  nsys <- .update_system(nsys, index_cache = .build_index_cache(nsys))
  idx <- match(.keys(new_vars), .keys(unknowns(nsys)))
  largs <- list(list(name = "p", vars = .problem_params(nsys)), list(name = "lnew", vars = limitnew),
                list(name = "lold", vars = limitold))
  fns <- lapply(limiters, function(l) .build_closure(l, largs, observed = .param_binding_eqs(nsys), out = "scalar"))
  post <- function(xnew, xold, p) {
    for (k in seq_along(fns)) xnew[idx[k]] <- fns[[k]](p, xnew[idx[k]], xold[idx[k]])
    xnew
  }
  list(sys = nsys, postcondition = post)
}

# ---------------------------------------------------------------------------
# Polynomial homotopy continuation

#' Homotopy continuation for polynomial systems
#'
#' `HomotopyContinuationProblem(sys, op)` builds a problem for a nonlinear
#' system whose equations are polynomials (or rational functions) in the
#' unknowns. `solve(prob, HomotopyContinuationJL())` tracks solution paths
#' from the roots of the total degree start system `x_i^d_i = 1` to the
#' roots of the system, in complex arithmetic, with the homotopy
#' `(1 - s) * gamma * G(x) + s * F(x)`. The solution holds the real root
#' closest to the initial values in `sol$u`; all roots found are in
#' `sol$original` (`real_solutions`, `solutions`, and the number of tracked
#' and failed paths). Roots at which a denominator vanishes are discarded.
#'
#' The number of paths is the product of the degrees of the equations, so
#' this is meant for small systems. Paths that diverge or end at a singular
#' root may be lost; `sol$original$n_failed` counts them, so the list of
#' roots is not guaranteed to be complete. ModelingToolkit delegates the
#' tracking to HomotopyContinuation.jl; the tracker here is native R.
#'
#' @param sys A complete time independent `System` with polynomial
#'   equations.
#' @param op Operating point: values of unknowns (used to choose the
#'   reported root) and parameters.
#' @param ... Solver options stored with the problem.
#' @param gamma The complex constant of the homotopy; a fixed generic value
#'   by default so results are reproducible.
#' @param max_paths The largest number of paths to track.
#' @return A `HomotopyContinuationProblem`, or the algorithm object.
#' @examples
#' x <- Variable("x"); y <- Variable("y")
#' sys <- complete(System(list(0 ~ x^2 + y^2 - 5, 0 ~ x * y - 2), name = "poly"))
#' prob <- HomotopyContinuationProblem(sys, list(x ~ 1, y ~ 1))
#' sol <- solve(prob)
#' sol$u
#' nrow(sol$original$real_solutions)
#' @export
HomotopyContinuationProblem <- function(sys, op = list(), ...) {
  .nl_check_system(sys, "HomotopyContinuationProblem")
  us <- unknowns(sys)
  n <- length(us)
  if (length(equations(sys)) != n)
    .stopf("`HomotopyContinuationProblem` requires as many equations as unknowns; found %d equations and %d unknowns.",
           length(equations(sys)), n)
  resid <- .nl_rewrite_call(.ode_rhs_exprs(sys), "homotopy", function(a, s) a)
  ph <- lapply(seq_len(n), function(i) Symbol(paste0(".hc", i)))
  low <- .nl_subst(resid, us, ph)
  nums <- vector("list", n); dens <- vector("list", n); degs <- integer(n)
  for (i in seq_len(n)) {
    fr <- fraction(together(low[[i]]))
    # A denominator free of the unknowns is a constant scale; keep the
    # equation as written so that, for example, exact rationals from floating
    # point coefficients do not inflate the coefficients the path tracker
    # works with.
    if (!any(.keys(free_symbols(fr[[2L]])) %in% .keys(ph))) fr <- list(low[[i]], as_sym(1))
    ok <- tryCatch(isTRUE(do.call(is_polynomial, c(list(fr[[1L]]), ph))) &&
                     isTRUE(do.call(is_polynomial, c(list(fr[[2L]]), ph))), error = function(e) FALSE)
    if (!ok)
      stop("Cannot convert the system to polynomial form: equation ", .fmt(equations(sys)[[i]]),
           " is not a polynomial in the unknowns. Allowed polynomial operations are `*, /, +, -, ^` ",
           "with integer exponents.", call. = FALSE)
    d <- do.call(total_degree, c(list(expand(fr[[1L]])), ph))
    if (!is.finite(d) || d < 1)
      stop("equation ", .fmt(equations(sys)[[i]]), " does not depend on the unknowns", call. = FALSE)
    nums[[i]] <- fr[[1L]]; dens[[i]] <- fr[[2L]]; degs[i] <- as.integer(d)
  }
  back <- .rule(ph, us)
  nums <- lapply(nums, .subst, rule = back)
  dens <- lapply(dens, .subst, rule = back)
  args <- .codegen_args(sys, with_t = FALSE)
  fn <- structure(list(f = .build_closure(nums, args, out = "vector", cse = FALSE),
                       jac = .build_closure(.partials(nums, us, .nl_hold(sys)), args, out = "matrix", cse = FALSE),
                       denominator = .build_closure(dens, args, out = "vector", cse = FALSE),
                       sys = sys, observed = function(expr) build_explicit_observed_function(sys, expr)),
                  class = c("HomotopyNonlinearFunction", "NonlinearFunction", "SciMLFunction"))
  vals <- .process_op(sys, op)
  u0 <- vals$u0
  if (anyNA(u0)) {
    g <- .nl_eval_guesses(guesses(sys), us[is.na(u0)], .problem_params(sys), vals$p)
    g[is.na(g)] <- 0
    u0[is.na(u0)] <- g
  }
  prob <- .new_problem("HomotopyContinuationProblem", f = fn, u0 = u0, p = vals$p, sys = sys,
                       kwargs = list(...))
  prob$op <- vals$op
  prob$degrees <- degs
  class(prob) <- c("HomotopyContinuationProblem", "AbstractNonlinearProblem", "SciMLProblem")
  prob
}

#' @rdname HomotopyContinuationProblem
#' @export
HomotopyContinuationJL <- function(gamma = complex(modulus = 1, argument = 0.6180339887 * pi), max_paths = 4096L) {
  .nl_alg("HomotopyContinuationJL", gamma = gamma, max_paths = as.integer(max_paths))
}

.hc_solve <- function(A, b) tryCatch(solve(A, b), error = function(e) NULL)

# Track one path of H(x, s) = (1 - s) gamma G(x) + s F(x) from s = 0 to 1:
# Euler predictor, Newton corrector, step halving on failure.
.hc_track <- function(F, JF, degs, x, gamma) {
  n <- length(x)
  G <- function(z) z^degs - 1
  JG <- function(z) diag(degs * z^(degs - 1), n)
  H <- function(z, s) (1 - s) * gamma * G(z) + s * F(z)
  Hx <- function(z, s) (1 - s) * gamma * JG(z) + s * JF(z)
  s <- 0; ds <- 0.02; streak <- 0L
  status <- "ok"
  for (step in seq_len(20000L)) {
    if (s >= 1) break
    ds <- min(ds, 1 - s)
    dxds <- .hc_solve(Hx(x, s), -(F(x) - gamma * G(x)))
    accepted <- FALSE
    if (!is.null(dxds) && all(is.finite(dxds))) {
      xc <- x + ds * dxds
      sn <- s + ds
      for (it in 1:4) {
        dx <- .hc_solve(Hx(xc, sn), -H(xc, sn))
        if (is.null(dx) || !all(is.finite(dx))) break
        xc <- xc + dx
        if (max(Mod(dx)) <= 1e-9 * max(1, max(Mod(xc)))) {
          accepted <- TRUE
          break
        }
      }
    }
    if (accepted) {
      x <- xc; s <- sn
      streak <- streak + 1L
      if (streak >= 3L) {
        ds <- min(2 * ds, 0.05)
        streak <- 0L
      }
      if (max(Mod(x)) > 1e8) {
        status <- "diverged"
        break
      }
    } else {
      ds <- ds / 2
      streak <- 0L
      if (ds < 1e-12) {
        status <- if (s > 0.99) "singular" else "failed"
        break
      }
    }
  }
  if (s < 1 && status == "ok") status <- "failed"
  if (status %in% c("ok", "singular")) {
    # endgame: Newton on the target system
    for (it in 1:50) {
      fx <- F(x)
      if (max(Mod(fx)) <= 1e-13 * max(1, max(Mod(x))^max(degs))) break
      Jx <- JF(x)
      dx <- .hc_solve(Jx, -fx)
      if (is.null(dx) || !all(is.finite(dx))) {
        dx <- tryCatch(qr.solve(Jx, -fx, tol = 1e-14), error = function(e) NULL)
        if (is.null(dx) || !all(is.finite(dx))) break
      }
      x <- x + dx
      if (max(Mod(dx)) <= 1e-15 * max(1, max(Mod(x)))) break
    }
    fx <- F(x)
    if (!all(is.finite(fx)) || max(Mod(fx)) > 1e-7 * max(1, max(Mod(x))^max(degs))) status <- "failed"
  }
  list(x = x, status = status)
}

#' @rdname solve_nonlinear
#' @export
solve.HomotopyContinuationProblem <- function(a, b, ...) {
  alg <- if (missing(b) || is.null(b)) HomotopyContinuationJL() else b
  p <- a$p
  n <- length(a$u0)
  degs <- a$degrees
  npaths <- prod(degs)
  if (npaths > alg$max_paths)
    .stopf("the total degree homotopy needs %s paths, more than max_paths = %d", format(npaths), alg$max_paths)
  hf <- a$f$numerator %||% a$f
  f <- hf$f; jac <- hf$jac; den <- a$f$denominator %||% hf$denominator
  F <- function(z) as.complex(f(z, p))
  JF <- function(z) matrix(as.complex(jac(z, p)), n, n)
  starts <- as.matrix(expand.grid(lapply(degs, function(d) seq_len(d) - 1L)))
  sols <- list(); nfail <- 0L; statuses <- character(nrow(starts))
  for (r in seq_len(nrow(starts))) {
    x0 <- exp(2i * pi * starts[r, ] / degs)
    tr <- .hc_track(F, JF, degs, x0, alg$gamma)
    statuses[r] <- tr$status
    if (tr$status == "failed" || tr$status == "diverged") {
      nfail <- nfail + (tr$status == "failed")
      next
    }
    dv <- as.complex(den(tr$x, p))
    if (any(Mod(dv) <= 1e-10)) {
      statuses[r] <- "infeasible"
      next
    }
    dup <- FALSE
    for (s in sols) if (max(Mod(s - tr$x)) <= 1e-6 * max(1, max(Mod(s)))) {
      dup <- TRUE
      break
    }
    if (!dup) sols[[length(sols) + 1L]] <- tr$x
  }
  nms <- names(a$u0)
  allsol <- if (length(sols)) do.call(rbind, sols) else matrix(complex(), 0L, n)
  # a root is real when every component is, and its real part is still a root
  maxdeg <- max(degs)
  isreal <- if (length(sols)) apply(allsol, 1L, function(z) {
    all(abs(Im(z)) <= 1e-8 * pmax(1, Mod(z))) &&
      max(Mod(F(Re(z)))) <= 1e-6 * max(1, max(abs(Re(z)))^maxdeg)
  }) else logical()
  realsol <- Re(allsol[isreal, , drop = FALSE])
  colnames(allsol) <- nms; colnames(realsol) <- nms
  original <- list(solutions = allsol, real_solutions = realsol, n_paths = nrow(starts), n_failed = nfail,
                   path_status = statuses)
  if (!nrow(realsol)) {
    retcode <- if (any(statuses == "infeasible")) "Infeasible" else "ConvergenceFailure"
    return(.nl_new_solution("NonlinearSolution", a, alg, a$u0, retcode, resid = Re(F(a$u0)),
                            extra = list(original = original)))
  }
  d <- apply(realsol, 1L, function(z) sum((z - a$u0)^2))
  u <- realsol[which.min(d), ]
  .nl_new_solution("NonlinearSolution", a, alg, u, "Success", resid = Re(F(u)) / Re(as.complex(den(u, p))),
                   stats = list(npaths = nrow(starts)), extra = list(original = original))
}
