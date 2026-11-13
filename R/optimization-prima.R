# OptimizationPRIMA: the derivative free solvers of PRIMA (UOBYQA, NEWUOA,
# BOBYQA, LINCOA and COBYLA) as OptimizationProblem algorithms, after
# lib/OptimizationPRIMA/src/OptimizationPRIMA.jl. COBYLA runs the native
# port of R/optimization-prima-cobyla.R in its Fortran variant.

.prima_alg <- function(name, traits) {
  structure(list(name = name, traits = utils::modifyList(list(callback = TRUE), traits), solve = .prima_solve),
            class = c(name, "PRIMASolvers", "OptimizationAlgorithm"))
}

#' PRIMA optimizers
#'
#' The derivative free solvers of Powell as modernized by PRIMA
#' (OptimizationPRIMA.jl), for [solve()] on an [OptimizationProblem()].
#'
#' * `UOBYQA()`: unconstrained, quadratic models by interpolation.
#' * `NEWUOA()`: unconstrained, quadratic models with minimum Frobenius norm
#'   updates.
#' * `BOBYQA()`: bound constraints.
#' * `LINCOA()`: bounds and linear constraints.
#' * `COBYLA()`: bounds, linear and nonlinear constraints, with linear
#'   models on a simplex. The constraints `lcons <= cons(u) <= ucons` of the
#'   problem are split, as in the Julia wrapper, into linear equalities,
#'   linear inequalities and nonlinear ones by their Hessians at a vector of
#'   ones (numerically when the problem has no constraint Hessian).
#'
#' `maxiters` sets the largest number of function evaluations (`maxfun`);
#' further arguments of [solve()] such as `rhobeg`, `rhoend`, `ftarget`,
#' `npt` (not for `UOBYQA()`, whose `npt` is `(n + 1)(n + 2)/2`), `scale` or
#' `iprint` go to the solver. `iprint` is a level of PRIMA.jl's `Message`:
#' 0 (`"MSG_NONE"`), 1 (`"MSG_EXIT"`, the exit message), 2 (`"MSG_RHO"`,
#' also each reduction of the trust region bound) or 3 (`"MSG_FEVL"`, also
#' each function evaluation); the messages of `UOBYQA()` and `NEWUOA()` are
#' those of libprima, printed to the console. Invalid `maxfun`, `npt`,
#' `ftarget` or `rhoend` values are revised with a warning, as in libprima.
#' `abstol`, `reltol` and `maxtime` are not supported, as in the Julia
#' wrapper.
#'
#' `UOBYQA()`, `NEWUOA()` and `COBYLA()` are ported; `BOBYQA()` and
#' `LINCOA()` are not yet ported and raise an error on [solve()].
#'
#' @return An algorithm object for [solve()].
#' @examples
#' prob <- OptimizationProblem(function(u, p) (u[1] - 1) * (u[1] - 1) + (u[2] + 2) * (u[2] + 2),
#'                             c(0, 0), numeric(), cons = function(u, p) u[1] + u[2],
#'                             lcons = -2, ucons = 0)
#' solve(prob, COBYLA())$u
#' rosen <- function(u, p) (1 - u[1]) * (1 - u[1]) + 100 * (u[2] - u[1] * u[1]) * (u[2] - u[1] * u[1])
#' solve(OptimizationProblem(rosen, c(-1.2, 1), numeric()), UOBYQA())$u
#' solve(OptimizationProblem(rosen, c(-1.2, 1), numeric()), NEWUOA(), npt = 4)$u
#' @name prima_optimizers
NULL

#' @rdname prima_optimizers
#' @export
UOBYQA <- function() .prima_alg("UOBYQA", list())

#' @rdname prima_optimizers
#' @export
NEWUOA <- function() .prima_alg("NEWUOA", list())

#' @rdname prima_optimizers
#' @export
BOBYQA <- function() .prima_alg("BOBYQA", list(bounds = TRUE))

#' @rdname prima_optimizers
#' @export
LINCOA <- function() .prima_alg("LINCOA", list(bounds = TRUE, cons = TRUE))

#' @rdname prima_optimizers
#' @export
COBYLA <- function() .prima_alg("COBYLA", list(bounds = TRUE, cons = TRUE))

# split the constraints into linear equalities, linear inequalities and
# nonlinear ones: by the Hessians at a vector of ones as the Julia wrapper, or
# by comparing each constraint with its linear model at two further points
.prima_split_cons <- function(cons, cons_j, cons_h, n, m) {
  x1 <- rep(1, n)
  c1 <- cons(x1)
  J <- if (is.function(cons_j)) cons_j(x1) else {
    h <- 1e-6
    vapply(seq_len(n), function(k) {
      e <- replace(numeric(n), k, h)
      (cons(x1 + e) - cons(x1 - e)) / (2 * h)
    }, numeric(m))
  }
  J <- matrix(J, m, n)
  linear <- if (is.function(cons_h)) {
    H <- cons_h(x1)
    vapply(seq_len(m), function(i) all(H[[i]] == 0), TRUE)
  } else {
    pts <- list(x1 + seq_len(n) / (n + 1), x1 - 2 * rev(seq_len(n)) / (n + 1))
    ok <- rep(TRUE, m)
    for (z in pts) {
      cz <- cons(z)
      model <- c1 + as.numeric(J %*% (z - x1))
      ok <- ok & abs(cz - model) <= 1e-6 * (1 + abs(cz))
    }
    ok
  }
  list(J = J, c1 = c1, linear = linear)
}

# the constraints lcons <= cons(u) <= ucons of the problem as linear
# equalities Aeq u = beq, linear inequalities Aineq u <= bineq and the
# indices of the nonlinear ones with a finite lower (nll) or upper (nlu) bound
.prima_cons_parts <- function(a, n) {
  fn <- a$f
  p <- a$p
  m <- length(a$lcons %||% numeric())
  cons <- function(x) as.numeric(fn$cons(x, p))
  cons_j <- if (is.function(fn$cons_j)) function(x) fn$cons_j(x, p)
  cons_h <- if (is.function(fn$cons_h)) function(x) fn$cons_h(x, p)
  sp <- .prima_split_cons(cons, cons_j, cons_h, n, m)
  lcons <- rep_len(as.numeric(a$lcons), m)
  ucons <- rep_len(as.numeric(a$ucons), m)
  # the Julia wrapper passes the bounds lcons and ucons as the right hand
  # sides and the raw constraint values as nonlinear inequalities; here the
  # constant part of each linear constraint and both bounds are taken into
  # account (see the porting notes)
  off <- sp$c1 - as.numeric(sp$J %*% rep(1, n))
  eq <- which(sp$linear & lcons == ucons)
  ineq <- which(sp$linear & lcons != ucons)
  il <- ineq[is.finite(lcons[ineq])]
  iu <- ineq[is.finite(ucons[ineq])]
  nl <- which(!sp$linear)
  list(cons = cons, lcons = lcons, ucons = ucons, nl = nl,
       nll = nl[is.finite(lcons[nl])], nlu = nl[is.finite(ucons[nl])],
       Aeq = if (length(eq)) sp$J[eq, , drop = FALSE], beq = if (length(eq)) lcons[eq] - off[eq],
       Aineq = if (length(il) + length(iu)) rbind(-sp$J[il, , drop = FALSE], sp$J[iu, , drop = FALSE]),
       bineq = if (length(il) + length(iu)) c(-(lcons[il] - off[il]), ucons[iu] - off[iu]))
}

# the checks of rhobeg and rhoend in PRIMA.jl (_check_rho)
.prima_check_rho <- function(rhobeg, rhoend) {
  if (!(is.finite(rhobeg) && rhobeg > 0)) stop("`rhobeg` must be finite and positive", call. = FALSE)
  if (!(is.finite(rhoend) && rhoend >= 0)) stop("`rhoend` must be finite and non-negative", call. = FALSE)
  if (rhoend > rhobeg) stop("`rhoend` must be less or equal `rhobeg`", call. = FALSE)
}

# the SciML return code of a PRIMA exit flag (sciml_prima_retcode)
.prima_retcode <- function(status) {
  switch(status, FTARGET_ACHIEVED = , SMALL_TR_RADIUS = "Success",
         MAXFUN_REACHED = , MAXTR_REACHED = "MaxIters",
         NO_SPACE_BETWEEN_BOUNDS = , ZERO_LINEAR_CONSTRAINT = , INVALID_INPUT = "InitialFailure",
         NAN_INF_X = , NAN_INF_F = , NAN_INF_MODEL = "Unstable",
         DAMAGING_ROUNDING = "ConvergenceFailure", "Failure")
}

.prima_solve <- function(a, alg, opts) {
  cls <- class(alg)[1L]
  if (!is.null(opts$maxtime) || !is.null(opts$abstol) || !is.null(opts$reltol))
    stop("maxtime, abstol and reltol kwargs not supported in ", cls, "()", call. = FALSE)
  if (identical(cls, "COBYLA")) return(.prima_solve_cobyla(a, alg, opts))
  if (!cls %in% c("UOBYQA", "NEWUOA"))
    stop(cls, "() of OptimizationPRIMA is not yet ported to native R in modelcompR", call. = FALSE)
  fn <- a$f
  p <- a$p
  x0 <- as.numeric(a$u0)
  n <- length(x0)
  kw <- opts[setdiff(names(opts), c("maxiters", "callback", "progress", "verbose", "seed"))]
  known <- c("rhobeg", "rhoend", "ftarget", "iprint", "scale", if (cls != "UOBYQA") "npt")
  unknown <- setdiff(names(kw), known)
  if (length(unknown)) stop("unsupported keyword argument `", unknown[1L], "` of ", cls, "()", call. = FALSE)
  has_cons <- is.function(fn$cons) && length(a$lcons %||% numeric()) > 0L
  if (has_cons)
    stop(IncompatibleOptimizerError(paste0("The algorithm ", cls, " does not support constraints. Either ",
                                           "remove the `cons` function passed to `OptimizationFunction` or use ",
                                           "a different algorithm.")))
  bounded <- !is.null(a$lb) || !is.null(a$ub)
  if (bounded)
    stop(IncompatibleOptimizerError(paste0("The algorithm ", cls, " does not support box constraints. ",
                                           "Either remove the `lb` or `ub` bounds passed to ",
                                           "`OptimizationProblem` or use a different algorithm.")))
  args <- list(rhobeg = 1, rhoend = NULL, ftarget = -Inf, maxfun = 500L * n, npt = 2L * n + 1L, iprint = 0L,
               scale = NULL)
  if (!is.null(opts$maxiters)) args$maxfun <- .opt_maxiters(opts$maxiters)
  args <- utils::modifyList(args, kw)
  if (is.null(args$rhoend)) args$rhoend <- 1e-6 * args$rhobeg
  if (cls != "UOBYQA" && !(n + 2 <= args$npt && args$npt <= (n + 1) * (n + 2) / 2))
    stop("`n + 2 <= npt <= (n + 1)*(n + 2)/2` must hold", call. = FALSE)
  .prima_check_rho(args$rhobeg, args$rhoend)
  # iprint is a PRIMA.Message of PRIMA.jl: MSG_NONE, MSG_EXIT, MSG_RHO or MSG_FEVL
  msgs <- c("MSG_NONE", "MSG_EXIT", "MSG_RHO", "MSG_FEVL")
  if (is.character(args$iprint)) args$iprint <- match(args$iprint, msgs) - 1L
  if (length(args$iprint) != 1L || !isTRUE(args$iprint %in% 0:3))
    stop("iprint must be 0 to 3 or one of ", paste(msgs, collapse = ", "), call. = FALSE)
  # the scaling of the variables of PRIMA.jl: the solver works on x / scale
  scl <- if (is.null(args$scale)) rep(1, n) else as.numeric(args$scale)
  if (length(scl) != n) stop("vector of scaling factors have incompatible number of elements", call. = FALSE)
  if (!all(is.finite(scl) & scl > 0)) stop("scaling factor(s) must be finite and positive", call. = FALSE)
  scaled <- !is.null(args$scale)
  f <- if (scaled) function(x) as.numeric(fn$f(scl * x, p))[1L] else function(x) as.numeric(fn$f(x, p))[1L]
  xs <- if (scaled) x0 / scl else x0
  common <- list(fun = f, x0 = xs, rhobeg = args$rhobeg, rhoend = args$rhoend, ftarget = args$ftarget,
                 maxfun = as.integer(args$maxfun), iprint = as.integer(args$iprint))
  res <- if (cls == "UOBYQA") do.call(.uo_uobyqa, common) else
    do.call(.nw_newuoa, c(common, list(npt = as.integer(args$npt))))
  status <- names(.pr_info)[match(res$info, .pr_info)]
  original <- list(fx = res$f, nf = res$nf, status = status, cstrv = res$cstrv %||% 0)
  .opt_solution(a, alg, if (scaled) res$x * scl else res$x, res$f, .prima_retcode(status), original)
}

.prima_solve_cobyla <- function(a, alg, opts) {
  fn <- a$f
  p <- a$p
  x0 <- as.numeric(a$u0)
  n <- length(x0)
  kw <- opts[setdiff(names(opts), c("maxiters", "callback", "progress", "verbose", "seed"))]
  f <- function(x) as.numeric(fn$f(x, p))[1L]
  m <- length(a$lcons %||% numeric())
  if (!is.function(fn$cons) || !m) stop("COBYLA requires constraints", call. = FALSE)
  lc <- .prima_cons_parts(a, n)
  nl <- lc$nl
  nll <- lc$nll
  nlu <- lc$nlu
  calcfc <- function(x) {
    cv <- if (length(nl)) lc$cons(x)
    list(f = f(x), c = c(lc$lcons[nll] - cv[nll], cv[nlu] - lc$ucons[nlu]))
  }
  args <- list(rhobeg = 1, rhoend = NULL, ftarget = -Inf, maxfun = 500L * n)
  if (!is.null(opts$maxiters)) args$maxfun <- as.integer(opts$maxiters)
  args <- utils::modifyList(args, kw)
  if (is.null(args$rhoend)) args$rhoend <- 1e-6 * args$rhobeg
  .prima_check_rho(args$rhobeg, args$rhoend)
  r0 <- calcfc(x0)
  res <- .pr_cobyla(calcfc, length(r0$c), x0, Aineq = lc$Aineq, bineq = lc$bineq, Aeq = lc$Aeq, beq = lc$beq,
                    xl = a$lb, xu = a$ub, f0 = r0$f, nlconstr0 = r0$c, rhobeg = args$rhobeg, rhoend = args$rhoend,
                    ftarget = args$ftarget, maxfun = args$maxfun, fortran = TRUE)
  status <- names(.pr_info)[match(res$info, .pr_info)]
  original <- list(fx = res$f, nf = res$nf, status = status, cstrv = res$cstrv, constr = res$constr)
  .opt_solution(a, alg, res$x, res$f, .prima_retcode(status), original)
}
