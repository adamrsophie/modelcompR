# SimpleIDSolve, after SimpleImplicitDiscreteSolve.jl v2.0.4
# (src/SimpleImplicitDiscreteSolve.jl): every step of an
# ImplicitDiscreteProblem f(u_next, u, p, t) = 0 is solved with
# SimpleNonlinearSolve's SimpleNewtonRaphson (src/raphson.jl, the iteration
# .simple_newton_raphson of R/optimization-simple.R) on the times
# t0:dt:tf of a Julia float range (Base.twiceprecision.jl).

# Base.rat: the continued fraction approximation n / d of x with
# |n|, |d| at most maxintfloat(Float32) = 2^24 (0 / 0 when there is none).
.sids_rat <- function(x) {
  y <- x
  a <- 1; d <- 1; b <- 0; cc <- 0
  m <- 16777216
  while (abs(y) <= m) {
    f <- trunc(y)
    y <- y - f
    an <- f * a + cc; cc <- a; a <- an
    bn <- f * b + d; d <- b; b <- bn
    if (max(abs(a), abs(b)) > m) return(c(cc, d))
    if (a / b == x) break
    y <- 1 / y
  }
  c(a, b)
}

# The elements of the Julia range start:step:stop of Float64 values: when
# start, step and stop are ratios of small integers the elements are
# (start_n + k step_n) / den (correctly rounded, as the TwicePrecision range),
# otherwise start + k step.
.sids_float_range <- function(start, step, stop) {
  if (step == 0) stop("range step cannot be zero", call. = FALSE)
  gcd <- function(x, y) if (y == 0) abs(x) else gcd(y, x %% y)
  st <- .sids_rat(step)
  if (st[2L] != 0 && st[1L] / st[2L] == step) {
    s0 <- .sids_rat(start)
    s1 <- .sids_rat(stop)
    if (s0[2L] != 0 && s1[2L] != 0 && s0[1L] / s0[2L] == start && s1[1L] / s1[2L] == stop) {
      den <- s0[2L] / gcd(s0[2L], st[2L]) * st[2L]
      m <- 2^53
      if (den != 0 && abs(start * den) <= m && abs(step * den) <= m) {
        start_n <- round(start * den)
        step_n <- round(step * den)
        len <- max(0, (den * s1[1L] - s1[2L] * start_n + step_n * s1[2L]) %/% (step_n * s1[2L]))
        between <- function(a, x, b) (a <= x && x <= b) || (b <= x && x <= a)
        if (between(start, start + (len - 1) * step, stop + step / 2) && !between(start, start + len * step, stop))
          return((start_n + (seq_len(len) - 1) * step_n) / den)
      }
    }
  }
  lf <- (stop - start) / step
  len <- if (lf < 0) 0 else if (lf == 0) 1 else {
    l <- round(lf) + 1
    stop2 <- start + (l - 1) * step
    l - ((start < stop && stop < stop2) + (start > stop && stop > stop2))
  }
  start + (seq_len(len) - 1) * step
}

# The Jacobian of f(x, u, p, t) by x: the problem's `jac` when it has one;
# otherwise the function is traced with symbolic vectors (as
# modelingtoolkitize) and differentiated exactly, which stands in for the
# ForwardDiff Jacobian of Julia; central differences when tracing fails.
.sids_jacobian <- function(prob) {
  jac <- prob$f$jac
  f <- prob$f$f
  if (is.function(jac)) return(function(x, u, p, t) as.matrix(jac(x, u, p, t)))
  n <- length(prob$u0)
  p <- prob$p
  fd <- function(x, u, p, t) .jac_central(function(z) as.numeric(f(z, u, p, t)), x)
  if (!is.numeric(p)) return(fd)
  xs <- lapply(seq_len(n), function(i) Symbol(paste0(".sidx", i)))
  us <- lapply(seq_len(n), function(i) Symbol(paste0(".sidu", i)))
  ps <- lapply(seq_along(p), function(i) Symbol(paste0(".sidp", i)))
  ts <- Symbol(".sidt")
  res <- tryCatch(.trace(f, list(.symvec(xs), .symvec(us, names(prob$u0)), .symvec(ps, names(p)), .symvec(list(ts)))),
                  error = function(e) NULL)
  if (is.null(res) || length(res) != n) return(fd)
  J <- .partials(res, xs, c(us, ps, list(ts)))
  args <- list(list(name = "x", vars = xs), list(name = "u", vars = us), list(name = "p", vars = ps),
               list(name = "t", vars = list(ts)))
  cj <- tryCatch(.build_closure(J, args, out = "matrix"), error = function(e) NULL)
  if (is.null(cj)) return(fd)
  function(x, u, p, t) matrix(as.numeric(cj(x, u, p, t)), n, n)
}

#' Simple implicit discrete solver
#'
#' `SimpleIDSolve()` solves an [ImplicitDiscreteProblem()] `f(u_next, u, p,
#' t) = 0` as SimpleImplicitDiscreteSolve.jl: the state is first made
#' consistent by solving `f(u, u0, p, t0) = 0` from `u0`, then each step
#' from `t0` to `tf` by `dt` (default 1, the times of the Julia range
#' `t0:dt:tf`) solves `f(u_next, u_prev, p, t_next) = 0` with the simple
#' Newton Raphson iteration of SimpleNonlinearSolve (absolute tolerance
#' `3e-13` on the largest residual, 1000 iterations). The Jacobian is the
#' problem's `jac(u_next, u, p, t)` when given, the exact derivative of the
#' traced function otherwise (central differences when the function cannot
#' be traced).
#'
#' The solution has `retcode` `"Default"` (as Julia's), `"InitialFailure"`
#' when the initial solve fails (the solution then holds `u0` at `t0`) or
#' `"ConvergenceFailure"` when a step fails (the solution ends with the last
#' iterate at the time of the failed step). Options of `solve()`: `dt`,
#' `save_everystep`, `save_start`, `save_end`; `adaptive` and `dense` must
#' be `FALSE`.
#'
#' @return An algorithm object.
#' @examples
#' f <- function(un, u, p, t) un + un * un * un - u
#' prob <- ImplicitDiscreteProblem(f, 2, c(0, 3))
#' sol <- solve(prob, SimpleIDSolve())
#' unlist(sol$u)
#' @export
SimpleIDSolve <- function() structure(list(name = "SimpleIDSolve"), class = c("SimpleIDSolve", "DiscreteAlgorithm"))

.sids_solve <- function(prob, opts) {
  dt <- opts$dt %||% 1
  if (isTRUE(opts$adaptive)) stop("SimpleIDSolve has fixed steps: adaptive must be FALSE", call. = FALSE)
  if (isTRUE(opts$dense)) stop("SimpleIDSolve has no dense output: dense must be FALSE", call. = FALSE)
  save_everystep <- opts$save_everystep %||% TRUE
  save_start <- opts$save_start %||% TRUE
  save_end <- opts$save_end %||% TRUE
  f <- prob$f$f
  p <- prob$p
  J <- .sids_jacobian(prob)
  newton <- function(uprev, t) {
    .simple_newton_raphson(function(x) as.numeric(f(x, uprev, p, t)), function(x) J(x, uprev, p, t),
                           uprev, maxiters = 1000L, abstol = 3e-13)
  }
  nm <- names(prob$u0)
  u0 <- stats::setNames(as.numeric(prob$u0), nm)
  t0 <- prob$tspan[1L]
  build <- function(ts, us, retcode) {
    U <- do.call(rbind, lapply(us, function(u) stats::setNames(as.numeric(u), nm)))
    if (is.null(U)) U <- matrix(numeric(), 0L, length(u0), dimnames = list(NULL, nm))
    .new_ode_solution(t = ts, u = U, prob = prob, alg = SimpleIDSolve(), retcode = retcode)
  }
  ini <- newton(u0, t0)
  if (ini$retcode != "Success") return(build(t0, list(u0), "InitialFailure"))
  ts <- .sids_float_range(t0, dt, prob$tspan[2L])
  u <- ini$u
  keep_t <- if (save_start) ts[1L] else numeric()
  keep_u <- if (save_start) list(u) else list()
  for (i in seq_along(ts)[-1L]) {
    s <- newton(u, ts[i])
    u <- s$u
    if (save_everystep || s$retcode != "Success") {
      keep_t <- c(keep_t, ts[i])
      keep_u[[length(keep_u) + 1L]] <- u
    }
    if (s$retcode != "Success") return(build(keep_t, keep_u, "ConvergenceFailure"))
  }
  if (!save_everystep && save_end && length(ts) > 1L) {
    keep_t <- c(keep_t, ts[length(ts)])
    keep_u[[length(keep_u) + 1L]] <- u
  }
  build(keep_t, keep_u, "Default")
}
