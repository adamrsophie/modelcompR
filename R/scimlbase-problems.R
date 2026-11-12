# Problem and function types of SciMLBase.jl that restate an ODE, SDE, DDE
# or BVP in structured form (src/problems/ode_problems.jl: SplitODEProblem,
# DynamicalODEProblem, SecondOrderODEProblem, IncrementingODEProblem;
# sde_problems.jl: SplitSDEProblem, DynamicalSDEProblem; dde_problems.jl:
# DynamicalDDEProblem, SecondOrderDDEProblem; bvp_problems.jl:
# TwoPointBVProblem, SecondOrderBVProblem, TwoPointSecondOrderBVProblem;
# scimlfunctions.jl: the matching function types), the analytical problem
# (analytical_problems.jl) and integral problems (integral_problems.jl,
# IntegralFunction, BatchIntegralFunction, IntegralSolution).
#
# Each structured problem is an ordinary problem of the package with the
# parts kept: a SplitODEProblem is an ODEProblem with f = f1 + f2 and
# `prob$f$f1`, `prob$f$f2`; a DynamicalODEProblem integrates c(v, u) with
# `prob$f$f1` (dv) and `prob$f$f2` (du). Every solver therefore applies, and
# the methods that use the structure (IMEX, symplectic and Nystrom methods)
# read the parts.

#' Structured differential equation problems
#'
#' `SplitODEProblem(f1, f2, u0, tspan, p)` is `du/dt = f1(u, p, t) + f2(u,
#' p, t)`, with `f1` the stiff (implicit) part for IMEX methods.
#' `DynamicalODEProblem(f1, f2, v0, u0, tspan, p)` is `dv/dt = f1(v, u, p,
#' t)`, `du/dt = f2(v, u, p, t)` on the state `c(v, u)` (velocities first),
#' and `SecondOrderODEProblem(f, du0, u0, tspan, p)` the second order system
#' `d2u/dt2 = f(du, u, p, t)` as such a problem with `f2 = v`.
#' `IncrementingODEProblem(f, u0, tspan, p)` takes `f(du, u, p, t, alpha,
#' beta) = alpha F(u, p, t) + beta du`, the form of low storage methods.
#' `SplitSDEProblem(f1, f2, g, u0, tspan, p)` and `DynamicalSDEProblem(f1,
#' f2, g, v0, u0, tspan, p)` (noise `g(v, u, p, t)` on the velocities) are
#' the stochastic versions, `DynamicalDDEProblem(f1, f2, v0, u0, h, tspan,
#' p)` (`f1(v, u, h, p, t)`, with `h(p, t)` the history of `c(v, u)`) and
#' `SecondOrderDDEProblem(f, du0, u0, h, tspan, p)` the delay versions.
#' `TwoPointBVProblem(f, bc, u0, tspan, p)` has separate conditions `bc =
#' list(bca, bcb)` with `bca(ua, p)` at the start and `bcb(ub, p)` at the
#' end; `SecondOrderBVProblem(f, bc, u0, tspan, p)` is `d2u/dt2 = f(du, u, p,
#' t)` with `bc(du, u, p, t)` (`du(s)` and `u(s)` give the solution), and
#' `TwoPointSecondOrderBVProblem()` the second order problem with `bc =
#' list(bca, bcb)` of `(dua, ua, p)` and `(dub, ub, p)`.
#'
#' The function types `SplitFunction()`, `DynamicalODEFunction()`,
#' `IncrementingODEFunction()`, `SplitSDEFunction()`,
#' `DynamicalSDEFunction()`, `DynamicalDDEFunction()`, `SDDEFunction()`,
#' `BVPFunction()`, `TwoPointBVPFunction()`, `DynamicalBVPFunction()` and
#' `TwoPointDynamicalBVPFunction()` hold the same parts and may be passed in
#' place of the functions.
#'
#' @param f,f1,f2,g Right-hand sides and diffusion, see Details.
#' @param bc Boundary conditions.
#' @param u0,v0,du0 Initial states and velocities (or guesses for BVPs).
#' @param h History function.
#' @param tspan,p Time span and parameters.
#' @param bcresid_prototype,twopoint Options of the BVP functions.
#' @param ... Problem options.
#' @return A problem or a function object.
#' @examples
#' prob <- SecondOrderODEProblem(function(du, u, p, t) -u, 0, 1, c(0, pi))
#' sol <- solve(prob, Tsit5())
#' round(sol(pi), 4)
#' sp <- SplitODEProblem(function(u, p, t) -2 * u, function(u, p, t) u, 1, c(0, 1))
#' solve(sp, Tsit5())(1)
#' @name structured_problems
NULL

.with_class <- function(x, cls) {
  class(x) <- unique(c(cls, class(x)))
  x
}

#' @rdname structured_problems
#' @export
SplitFunction <- function(f1, f2, ...)
  structure(list(f1 = f1, f2 = f2, f = function(u, p, t) f1(u, p, t) + f2(u, p, t), ...),
            class = c("SplitFunction", "ODEFunction", "SciMLFunction"))

#' @rdname structured_problems
#' @export
SplitODEProblem <- function(f1, f2 = NULL, u0, tspan, p = numeric(), ...) {
  fn <- if (inherits(f1, "SplitFunction")) f1 else SplitFunction(f1, f2)
  if (inherits(f1, "SplitFunction") && !missing(f2) && !is.null(f2)) {
    # SplitODEProblem(SplitFunction, u0, tspan, p)
    return(SplitODEProblem(f1, NULL, f2, u0, if (missing(tspan)) numeric() else tspan, ...))
  }
  prob <- ODEProblem(fn$f, u0, tspan, p, ...)
  for (nm in c("f1", "f2")) prob$f[[nm]] <- fn[[nm]]
  class(prob$f) <- unique(c("SplitFunction", class(prob$f)))
  .with_class(prob, "SplitODEProblem")
}

#' @rdname structured_problems
#' @export
DynamicalODEFunction <- function(f1, f2, ...)
  structure(list(f1 = f1, f2 = f2, ...), class = c("DynamicalODEFunction", "ODEFunction", "SciMLFunction"))

.dyn_split <- function(z, nv) list(v = z[seq_len(nv)], u = z[nv + seq_len(length(z) - nv)])

#' @rdname structured_problems
#' @export
DynamicalODEProblem <- function(f1, f2, v0, u0, tspan, p = numeric(), ...) {
  fn <- if (inherits(f1, "DynamicalODEFunction")) {
    # DynamicalODEProblem(DynamicalODEFunction, v0, u0, tspan, p)
    args <- list(f2, v0, u0, tspan)
    if (!missing(p)) args <- c(args, list(p))
    return(do.call(DynamicalODEProblem, c(list(f1$f1, f1$f2), args, list(...))))
  } else DynamicalODEFunction(f1, f2)
  v0 <- as.numeric(unlist(v0)); u0 <- as.numeric(unlist(u0))
  nv <- length(v0)
  f <- function(z, p, t) {
    s <- .dyn_split(z, nv)
    c(f1(s$v, s$u, p, t), f2(s$v, s$u, p, t))
  }
  prob <- ODEProblem(f, c(v0, u0), tspan, p, ...)
  prob$f$f1 <- fn$f1; prob$f$f2 <- fn$f2
  prob$nv <- nv
  .with_class(prob, "DynamicalODEProblem")
}

#' @rdname structured_problems
#' @export
SecondOrderODEProblem <- function(f, du0, u0, tspan, p = numeric(), ...) {
  prob <- DynamicalODEProblem(f, function(v, u, p, t) v, du0, u0, tspan, p, ...)
  prob$second_order <- TRUE
  .with_class(prob, "SecondOrderODEProblem")
}

#' @rdname structured_problems
#' @export
IncrementingODEFunction <- function(f, ...) {
  F <- if (length(formals(f)) >= 6L) function(u, p, t) f(0 * u, u, p, t, 1, 0) else f
  structure(list(f = F, incrementing = f, ...), class = c("IncrementingODEFunction", "ODEFunction", "SciMLFunction"))
}

#' @rdname structured_problems
#' @export
IncrementingODEProblem <- function(f, u0, tspan, p = numeric(), ...) {
  fn <- if (inherits(f, "IncrementingODEFunction")) f else IncrementingODEFunction(f)
  prob <- ODEProblem(fn$f, u0, tspan, p, ...)
  prob$f$incrementing <- fn$incrementing
  .with_class(prob, "IncrementingODEProblem")
}

#' @rdname structured_problems
#' @export
SplitSDEFunction <- function(f1, f2, g, ...)
  structure(list(f1 = f1, f2 = f2, g = g, f = function(u, p, t) f1(u, p, t) + f2(u, p, t), ...),
            class = c("SplitSDEFunction", "SDEFunction", "SciMLFunction"))

#' @rdname structured_problems
#' @export
SplitSDEProblem <- function(f1, f2, g, u0, tspan, p = numeric(), ...) {
  fn <- SplitSDEFunction(f1, f2, g)
  prob <- SDEProblem(fn$f, g, u0, tspan, p, ...)
  prob$f$f1 <- f1; prob$f$f2 <- f2
  .with_class(prob, "SplitSDEProblem")
}

#' @rdname structured_problems
#' @export
DynamicalSDEFunction <- function(f1, f2, g, ...)
  structure(list(f1 = f1, f2 = f2, g = g, ...), class = c("DynamicalSDEFunction", "SDEFunction", "SciMLFunction"))

#' @rdname structured_problems
#' @export
DynamicalSDEProblem <- function(f1, f2, g, v0, u0, tspan, p = numeric(), ...) {
  v0 <- as.numeric(unlist(v0)); u0 <- as.numeric(unlist(u0))
  nv <- length(v0); nu <- length(u0)
  f <- function(z, p, t) {
    s <- .dyn_split(z, nv)
    c(f1(s$v, s$u, p, t), f2(s$v, s$u, p, t))
  }
  gg <- function(z, p, t) {
    s <- .dyn_split(z, nv)
    c(rep_len(as.numeric(g(s$v, s$u, p, t)), nv), numeric(nu))
  }
  prob <- SDEProblem(f, gg, c(v0, u0), tspan, p, ...)
  prob$f$f1 <- f1; prob$f$f2 <- f2; prob$f$g_dyn <- g
  prob$nv <- nv
  .with_class(prob, "DynamicalSDEProblem")
}

#' @rdname structured_problems
#' @export
DynamicalDDEFunction <- function(f1, f2, ...)
  structure(list(f1 = f1, f2 = f2, ...), class = c("DynamicalDDEFunction", "DDEFunction", "SciMLFunction"))

#' @rdname structured_problems
#' @export
DynamicalDDEProblem <- function(f1, f2, v0, u0, h, tspan, p = numeric(), ...) {
  v0 <- as.numeric(unlist(v0)); u0 <- as.numeric(unlist(u0))
  nv <- length(v0)
  f <- function(z, hist, p, t) {
    s <- .dyn_split(z, nv)
    c(f1(s$v, s$u, hist, p, t), f2(s$v, s$u, hist, p, t))
  }
  prob <- DDEProblem(f, c(v0, u0), h, tspan, p, ...)
  prob$f$f1 <- f1; prob$f$f2 <- f2
  prob$nv <- nv
  .with_class(prob, "DynamicalDDEProblem")
}

#' @rdname structured_problems
#' @export
SecondOrderDDEProblem <- function(f, du0, u0, h, tspan, p = numeric(), ...) {
  prob <- DynamicalDDEProblem(f, function(v, u, hist, p, t) v, du0, u0, h, tspan, p, ...)
  .with_class(prob, "SecondOrderDDEProblem")
}

#' @rdname structured_problems
#' @export
SDDEFunction <- function(f, g, ...)
  structure(list(f = f, g = g, mass_matrix = NULL, sys = NULL, ...), class = c("SDDEFunction", "DDEFunction", "SciMLFunction"))

#' @rdname structured_problems
#' @export
BVPFunction <- function(f, bc, bcresid_prototype = NULL, twopoint = FALSE, ...)
  structure(list(f = f, bc = bc, bcresid_prototype = bcresid_prototype, twopoint = twopoint, ...),
            class = c("BVPFunction", "SciMLFunction"))

#' @rdname structured_problems
#' @export
TwoPointBVPFunction <- function(f, bc, ...) BVPFunction(f, bc, twopoint = TRUE, ...)

#' @rdname structured_problems
#' @export
DynamicalBVPFunction <- function(f, bc, ...)
  structure(list(f = f, bc = bc, twopoint = FALSE, ...), class = c("DynamicalBVPFunction", "BVPFunction", "SciMLFunction"))

#' @rdname structured_problems
#' @export
TwoPointDynamicalBVPFunction <- function(f, bc, ...)
  structure(list(f = f, bc = bc, twopoint = TRUE, ...),
            class = c("TwoPointDynamicalBVPFunction", "DynamicalBVPFunction", "BVPFunction", "SciMLFunction"))

.twopoint_bc <- function(bc) {
  if (!is.list(bc) || length(bc) != 2L) stop("a two point BVP takes bc = list(bca, bcb)", call. = FALSE)
  force(bc)
  function(sol, p, t) c(bc[[1L]](sol(t[1L]), p), bc[[2L]](sol(t[2L]), p))
}

#' @rdname structured_problems
#' @export
TwoPointBVProblem <- function(f, bc, u0, tspan, p = numeric(), ...) {
  if (inherits(f, "BVPFunction")) return(TwoPointBVProblem(f$f, f$bc, bc, u0, if (missing(tspan)) numeric() else tspan, ...))
  prob <- BVProblem(f, .twopoint_bc(bc), u0, tspan, p, ...)
  prob$bca <- bc[[1L]]; prob$bcb <- bc[[2L]]
  .with_class(prob, "TwoPointBVProblem")
}

# A second order BVP on c(du, u): the residual sees du(s) and u(s).
.second_order_bvp <- function(f, bcfun, u0, tspan, p, ...) {
  u0 <- as.numeric(unlist(u0))
  n <- length(u0)
  du0 <- numeric(n)
  rhs <- function(z, p, t) c(f(z[seq_len(n)], z[n + seq_len(n)], p, t), z[seq_len(n)])
  BVProblem(rhs, bcfun, c(du0, u0), tspan, p, ...)
}

#' @rdname structured_problems
#' @export
SecondOrderBVProblem <- function(f, bc, u0, tspan, p = numeric(), ...) {
  n <- length(unlist(u0))
  bcfun <- function(sol, p, t) bc(function(s) sol(s)[seq_len(n)], function(s) sol(s)[n + seq_len(n)], p, t)
  .with_class(.second_order_bvp(f, bcfun, u0, tspan, p, ...), "SecondOrderBVProblem")
}

#' @rdname structured_problems
#' @export
TwoPointSecondOrderBVProblem <- function(f, bc, u0, tspan, p = numeric(), ...) {
  n <- length(unlist(u0))
  bcfun <- function(sol, p, t) {
    a <- sol(t[1L]); b <- sol(t[2L])
    c(bc[[1L]](a[seq_len(n)], a[n + seq_len(n)], p), bc[[2L]](b[seq_len(n)], b[n + seq_len(n)], p))
  }
  .with_class(.second_order_bvp(f, bcfun, u0, tspan, p, ...), c("TwoPointSecondOrderBVProblem", "SecondOrderBVProblem"))
}

# ---------------------------------------------------------------------------
# Analytical problems

#' Analytical problems
#'
#' `AnalyticalProblem(f, u0, tspan, p)` is a problem whose solution is known:
#' `f(u0, p, t)` is the state at `t`. `solve(prob, dt = , saveat = )`
#' evaluates it on the requested times (every `dt`, `saveat`, or the two
#' ends of `tspan`) and returns an [ODESolution()] interpolated linearly.
#'
#' @param f Solution `f(u0, p, t)`.
#' @param u0,tspan,p Initial state, time span and parameters.
#' @param ... Options passed to `solve()`.
#' @return A problem.
#' @examples
#' prob <- AnalyticalProblem(function(u0, p, t) u0 * exp(-t), 1, c(0, 1))
#' solve(prob, dt = 0.5)$u
#' @name analytical_problems
NULL

#' @rdname analytical_problems
#' @export
AnalyticalProblem <- function(f, u0, tspan, p = numeric(), ...)
  structure(list(f = f, u0 = as.numeric(unlist(u0)), tspan = .as_tspan(tspan), p = p, kwargs = list(...)),
            class = c("AnalyticalProblem", "SciMLProblem"))

#' @export
solve.AnalyticalProblem <- function(a, b = NULL, ...) {
  opts <- utils::modifyList(a$kwargs, list(...))
  ts <- if (!is.null(opts$saveat) && length(opts$saveat) > 1L) as.numeric(opts$saveat)
        else if (!is.null(opts$dt)) unique(c(seq(a$tspan[1L], a$tspan[2L], by = opts$dt), a$tspan[2L]))
        else if (!is.null(opts$saveat)) unique(c(seq(a$tspan[1L], a$tspan[2L], by = opts$saveat), a$tspan[2L]))
        else a$tspan
  U <- do.call(rbind, lapply(ts, function(s) as.numeric(a$f(a$u0, a$p, s))))
  pr <- ODEProblem(function(u, p, t) 0 * u, a$u0, a$tspan, a$p)
  ODESolution(ts, U, pr, alg = list(name = "AnalyticalSolution"))
}

# ---------------------------------------------------------------------------
# Integral problems

#' Integral problems
#'
#' `IntegralProblem(f, domain, p)` asks for the integral of `f(u, p)` over
#' `domain = c(lb, ub)` (vectors `lb` and `ub` give a box; infinite bounds
#' are allowed). `IntegralFunction(f)` is the integrand,
#' `BatchIntegralFunction(f)` an integrand evaluated at many points at once
#' (`f(U, p)` with one point per column of `U`, returning one value or
#' column per point). `SampledIntegralProblem(y, x, dim)` integrates sampled
#' values `y` over the points `x` (along dimension `dim` of an array `y`).
#'
#' `solve(prob)` integrates: one dimension adaptively with Gauss-Kronrod
#' quadrature ([stats::integrate()]), several dimensions by nested
#' quadrature; vector valued integrands are integrated per component.
#' `solve(sampled, alg)` uses the trapezoidal rule (`TrapezoidalRule()`, the
#' default) or Simpson's rule (`SimpsonsRule()`). The result is an
#' `IntegralSolution` with `u` (the integral), `resid` (the error estimate)
#' and `retcode`. Options: `abstol`, `reltol` (default `1e-8`), `maxiters`.
#'
#' @param f Integrand.
#' @param domain Bounds `c(lb, ub)` or `list(lb, ub)`.
#' @param p Parameters.
#' @param integrand_prototype,max_batch Options of the Julia signature.
#' @param y,x,dim Sampled values, points and dimension.
#' @param u,resid,retcode Solution fields.
#' @param prob,alg A problem and an algorithm.
#' @param ... Options passed to `solve()`.
#' @return A problem, a function, an algorithm or a solution.
#' @examples
#' solve(IntegralProblem(function(x, p) exp(-x^2), c(-Inf, Inf)))$u
#' solve(IntegralProblem(function(x, p) x[1] * x[2], list(c(0, 0), c(1, 2))))$u
#' x <- seq(0, pi, length.out = 101)
#' solve(SampledIntegralProblem(sin(x), x), SimpsonsRule())$u
#' @name integral_problems
NULL

#' @rdname integral_problems
#' @export
IntegralFunction <- function(f, integrand_prototype = NULL)
  structure(list(f = f, integrand_prototype = integrand_prototype, batch = FALSE),
            class = c("IntegralFunction", "SciMLFunction"))

#' @rdname integral_problems
#' @export
BatchIntegralFunction <- function(f, integrand_prototype = NULL, max_batch = 1000L)
  structure(list(f = f, integrand_prototype = integrand_prototype, batch = TRUE, max_batch = max_batch),
            class = c("BatchIntegralFunction", "IntegralFunction", "SciMLFunction"))

#' @rdname integral_problems
#' @export
IntegralProblem <- function(f, domain, p = numeric(), ...) {
  fn <- if (inherits(f, "IntegralFunction")) f else IntegralFunction(f)
  if (is.list(domain)) { lb <- as.numeric(domain[[1L]]); ub <- as.numeric(domain[[2L]]) }
  else if (length(domain) == 2L) { lb <- domain[1L]; ub <- domain[2L] }
  else stop("give domain = c(lb, ub) or list(lb, ub)", call. = FALSE)
  structure(list(f = fn, domain = list(lb, ub), lb = lb, ub = ub, p = p, kwargs = list(...)),
            class = c("IntegralProblem", "SciMLProblem"))
}

#' @rdname integral_problems
#' @export
SampledIntegralProblem <- function(y, x, dim = NULL, ...) {
  dim <- dim %||% (if (is.null(base::dim(y))) 1L else length(base::dim(y)))
  structure(list(y = y, x = as.numeric(x), dim = as.integer(dim), kwargs = list(...)),
            class = c("SampledIntegralProblem", "SciMLProblem"))
}

#' @rdname integral_problems
#' @export
IntegralSolution <- function(u, resid = NA_real_, retcode = "Success", prob = NULL, alg = NULL)
  structure(list(u = u, resid = resid, retcode = retcode, prob = prob, alg = alg), class = "IntegralSolution")

#' @export
print.IntegralSolution <- function(x, ...) {
  cat("IntegralSolution with retcode ", x$retcode, "\nu: ", sep = "")
  cat(format(x$u), "\n")
  invisible(x)
}

#' @rdname integral_problems
#' @export
TrapezoidalRule <- function() structure(list(name = "TrapezoidalRule"), class = c("TrapezoidalRule", "SampledIntegralAlgorithm"))

#' @rdname integral_problems
#' @export
SimpsonsRule <- function() structure(list(name = "SimpsonsRule"), class = c("SimpsonsRule", "SampledIntegralAlgorithm"))

# value of the integrand at one point, as a vector
.int_eval <- function(fn, x, p) {
  if (isTRUE(fn$batch)) as.numeric(fn$f(matrix(x, ncol = 1L), p)) else as.numeric(fn$f(x, p))
}

#' @export
solve.IntegralProblem <- function(a, b = NULL, ...) {
  opts <- utils::modifyList(a$kwargs, list(...))
  reltol <- opts$reltol %||% 1e-8
  abstol <- opts$abstol %||% 1e-8
  maxit <- opts$maxiters %||% 1000L
  fn <- a$f; p <- a$p; lb <- a$lb; ub <- a$ub
  d <- length(lb)
  nout <- length(.int_eval(fn, ifelse(is.finite(lb), lb, ifelse(is.finite(ub), ub, 0)) +
                             ifelse(is.finite(ub - lb), (ub - lb) / 3, 0), p))
  ok <- TRUE
  err_total <- numeric(nout)
  # nested adaptive quadrature over the last coordinates
  integ_dim <- function(k, fixed, comp) {
    g <- function(xs) vapply(xs, function(xk) {
      x <- c(fixed, xk)
      if (k == d) .int_eval(fn, x, p)[comp] else integ_dim(k + 1L, x, comp)
    }, 0)
    r <- tryCatch(stats::integrate(g, lb[k], ub[k], rel.tol = reltol, abs.tol = abstol, subdivisions = maxit),
                  error = function(e) { ok <<- FALSE; list(value = NA_real_, abs.error = NA_real_) })
    if (k == 1L) err_total[comp] <<- r$abs.error
    r$value
  }
  vals <- vapply(seq_len(nout), function(comp) integ_dim(1L, numeric(), comp), 0)
  IntegralSolution(vals, err_total, if (ok && all(is.finite(vals))) "Success" else "Failure", a,
                   b %||% list(name = if (d == 1L) "QuadGKJL" else "NestedQuadrature"))
}

#' @export
solve.SampledIntegralProblem <- function(a, b = TrapezoidalRule(), ...) {
  x <- a$x; y <- a$y
  rule <- function(v) {
    h <- diff(x)
    if (inherits(b, "SimpsonsRule")) {
      n <- length(x)
      if (n < 3L) return(sum(h * (v[-1L] + v[-n]) / 2))
      # composite Simpson on nonuniform points, pairwise intervals
      s <- 0; i <- 1L
      while (i + 2L <= n) {
        h0 <- h[i]; h1 <- h[i + 1L]
        s <- s + (h0 + h1) / 6 * ((2 - h1 / h0) * v[i] + (h0 + h1)^2 / (h0 * h1) * v[i + 1L] + (2 - h0 / h1) * v[i + 2L])
        i <- i + 2L
      }
      if (i < n) s <- s + h[i] * (v[i] + v[i + 1L]) / 2
      return(s)
    }
    sum(h * (v[-1L] + v[-length(v)]) / 2)
  }
  u <- if (is.null(dim(y))) rule(as.numeric(y)) else apply(y, setdiff(seq_along(dim(y)), a$dim), rule)
  IntegralSolution(u, NA_real_, "Success", a, b)
}
