# Adjoint sensitivities of fully implicit DAEs, after SciMLSensitivity.jl
# v7.119 src/dae_adjoint.jl (DAEAdjointProblem, DAEAdjointSensitivityFunction,
# DAEAdjointResidual, dae_adjoint_structure), the DAE branches of
# src/adjoint_common.jl (ReverseLossCallback: the consistent jumps of the
# discrete cost) and src/sensitivity_interface.jl (dG/du0 = w(t0) and the
# parameter corrections of the jumps).
#
# For F(du, u, p, t) = 0 the adjoint of Cao, Li, Petzold and Serban (SIAM
# J. Sci. Comput. 24(3), 2003) is written for z = (w, lambda), w = F_du'
# lambda, as the semi explicit system
#   w_j' = [F_u' lambda]_j - g_u[j]          (j differential)
#   0    = w_j                               (j algebraic)
#   0    = w_j - [F_du' lambda]_j            (j differential)
#   0    = [F_u' lambda]_j - g_u[j]          (j algebraic)
#   grad' = F_p' lambda - g_p                (InterpolatingAdjoint)
# solved backward from z(T) = 0, as a singular mass matrix ODEProblem for
# ODE algorithms or as a DAEProblem for DAE algorithms. dG/du0 = w(t0).
#
# Julia differentiates the residual by automatic differentiation and reads
# du(t) from the derivative of the interpolant, sol(t, Val{1}). Here a
# problem built from a System gets the Jacobians F_u, F_du and F_p
# symbolically and one built from an R function uses its Jacobian
# jac(du, u, p, gamma, t) and central differences; du(t) is the central
# difference of the step interpolant.

# The forward trajectory with its time derivative.
.sens_dae_trajectory <- function(sol) {
  tr <- .sens_trajectory(sol)
  steps <- .sol_env(sol)$steps
  starts <- vapply(steps, `[[`, 0, "t0")
  sgn <- if (length(starts) > 1L && starts[length(starts)] < starts[1L]) -1 else 1
  tr$dy <- function(t) {
    st <- steps[[max(1L, min(length(steps), findInterval(sgn * t, sgn * starts)))]]
    h <- .Machine$double.eps^(1 / 3) * abs(st$dt)
    if (h == 0) return(0 * st$y0)
    (.step_value(st, t + h) - .step_value(st, t - h)) / (2 * h)
  }
  tr
}

# Symbolic Jacobians of the residual of a DAE system by u, du and p.
.sens_dae_symbolic <- function(sys) {
  us <- unknowns(sys)
  iv <- .f(sys, "iv")
  obs <- .codegen_observed(sys)
  rule <- if (length(obs)) .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs"))
  rhs <- lapply(equations(sys), function(e) {
    r <- if (.is_zero(e$lhs)) e$rhs else e$rhs - e$lhs
    if (is.null(rule)) r else .subst_fixpoint(r, rule)
  })
  dus <- lapply(us, function(u) .make_derivative(u, iv, 1L))
  ps <- .problem_params(sys)
  others <- c(.parameters_of(sys, TRUE), list(iv))
  args <- .codegen_args(sys, du = TRUE)
  list(Ju = .build_closure(.partials(rhs, us, c(dus, others)), args, out = "matrix"),
       Jdu = .build_closure(.partials(rhs, dus, c(us, others)), args, out = "matrix"),
       Jp = if (length(ps)) .build_closure(.partials(rhs, ps, c(us, dus, list(iv))), args, out = "matrix"))
}

# The derivative bundle of a DAE problem along its solution: the residual
# Jacobians Ju, Jdu, Jp at (du, u, p, t) and the vector Jacobian product
# vjp(l, u, p, t) at the trajectory's du(t).
.sens_dae_fns <- function(prob, sol, sensealg = NULL) {
  n <- length(prob$u0)
  np <- length(prob$p)
  F <- prob$f$f
  tr <- .sens_dae_trajectory(sol)
  sys <- prob$sys
  out <- list(n = n, np = np, tr = tr)
  if (!is.null(sys) && inherits(sys, "System") && !identical(sensealg$autojacvec, FALSE)) {
    d <- .sens_dae_symbolic(sys)
    out$Ju <- d$Ju
    out$Jdu <- d$Jdu
    out$Jp <- function(du, u, p, t) if (np) d$Jp(du, u, p, t) else matrix(0, n, 0L)
  } else {
    jac <- prob$f$jac
    if (is.function(jac)) {
      out$Ju <- function(du, u, p, t) as.matrix(jac(du, u, p, 0, t))
      out$Jdu <- function(du, u, p, t) as.matrix(jac(du, u, p, 1, t)) - as.matrix(jac(du, u, p, 0, t))
    } else {
      out$Ju <- function(du, u, p, t) .sens_fd_jac(function(z) F(du, z, p, t), u)
      out$Jdu <- function(du, u, p, t) .sens_fd_jac(function(z) F(z, u, p, t), du)
    }
    out$Jp <- function(du, u, p, t) if (!np) matrix(0, n, 0L) else
      .sens_fd_jac(function(q) F(du, u, .sens_pvec(p, q), t), as.numeric(p))
  }
  out$vjp <- function(l, u, p, t) {
    du <- tr$dy(t)
    list(u = as.numeric(crossprod(out$Ju(du, u, p, t), l)), du = as.numeric(crossprod(out$Jdu(du, u, p, t), l)),
         p = as.numeric(crossprod(out$Jp(du, u, p, t), l)))
  }
  out
}

# LU without a zero pivot (issuccess(lu(A, check = false))).
.sens_lu_ok <- function(A) !inherits(tryCatch(solve(A, tol = 0), error = identity), "error")

.SENS_DAE_HIGH_INDEX <- paste(
  "The algebraic subsystem of the DAE could not be classified as index-1 (dF_alg/du_alg nonsingular) or",
  "Hessenberg index-2 (dF_alg/du_alg identically zero with dF_alg/du_diff * (dF_diff/d(du)_diff)^-1 *",
  "dF_diff/du_alg nonsingular) at the terminal time of the forward solution. Mixed index-1/index-2 systems",
  "and DAEs of index 3 or higher are not currently supported by the `DAEProblem` adjoint. Consider",
  "index-reducing the system (e.g. via `mtkcompile`) before differentiating it.")

# dae_adjoint_structure: algebraic equations are the zero rows of F_du,
# algebraic variables its zero columns (or !differential_vars).
.sens_dae_structure <- function(fns, p, t, differential_vars) {
  n <- fns$n
  y <- fns$tr$at(t)
  dy <- fns$tr$dy(t)
  Jdu <- fns$Jdu(dy, y, p, t)
  alg_eq <- which(rowSums(Jdu != 0) == 0)
  algv <- if (is.null(differential_vars)) which(colSums(Jdu != 0) == 0) else which(!differential_vars)
  diff_eq <- setdiff(seq_len(n), alg_eq)
  diffv <- setdiff(seq_len(n), algv)
  if (length(alg_eq) != length(algv))
    stop("The number of algebraic equations (rows of dF/d(du) that are identically zero) does not match the ",
         "number of algebraic variables (from `differential_vars`, or columns of dF/d(du) that are identically ",
         "zero). Adjoint sensitivity analysis for `DAEProblem` requires a square algebraic subsystem. If ",
         "`differential_vars` was not provided to the `DAEProblem`, provide it explicitly so that the ",
         "algebraic variables can be identified structurally.", call. = FALSE)
  index2 <- FALSE
  if (length(algv)) {
    Ju <- fns$Ju(dy, y, p, t)
    dhda <- Ju[alg_eq, algv, drop = FALSE]
    if (!.sens_lu_ok(dhda)) {
      if (any(dhda != 0)) stop(.SENS_DAE_HIGH_INDEX, call. = FALSE)
      Fdd <- Jdu[diff_eq, diffv, drop = FALSE]
      if (!.sens_lu_ok(Fdd)) stop(.SENS_DAE_HIGH_INDEX, call. = FALSE)
      N <- t(solve(Fdd, Ju[diff_eq, algv, drop = FALSE]))
      if (!.sens_lu_ok(N %*% t(Ju[alg_eq, diffv, drop = FALSE]))) stop(.SENS_DAE_HIGH_INDEX, call. = FALSE)
      index2 <- TRUE
    }
  }
  list(diff = diffv, alg = algv, diff_eq = diff_eq, alg_eq = alg_eq, index2 = index2)
}

#' Adjoint problem of a fully implicit DAE
#'
#' `DAEAdjointProblem(sol, sensealg, alg, ...)` builds the backward adjoint
#' problem of the solution of a [DAEProblem()] `F(du, u, p, t) = 0`, the
#' augmented adjoint of Cao, Li, Petzold and Serban (2003) for `z = (w,
#' lambda)` with `w = F_du' lambda` (and the parameter gradient for
#' [InterpolatingAdjoint()]):
#'
#' `w_j' = (F_u' lambda - g_u)_j` and `0 = w_j - (F_du' lambda)_j` for the
#' differential variables, `0 = w_j` and `0 = (F_u' lambda - g_u)_j` for the
#' algebraic ones, solved backward from `z(T) = 0`. With an ODE algorithm
#' (a mass matrix method such as [Rodas5P()] or [FBDF()]) it is a singular
#' mass matrix `ODEProblem`, with a DAE algorithm ([DFBDF()], the default
#' when the forward solve used one) a `DAEProblem`. Discrete costs jump `w`
#' by their gradient, corrected so that the jump is consistent when the cost
#' involves algebraic variables (index 1) or the DAE is Hessenberg index 2.
#' Then `dG/du0 = w(t0)` and `dG/dp = int (g_p - lambda' F_p) dt` (the
#' state of [InterpolatingAdjoint()], or a quadrature for
#' [QuadratureAdjoint()], [GaussAdjoint()] and [GaussKronrodAdjoint()]).
#'
#' [adjoint_sensitivities()] uses it for DAE solutions. The forward
#' solution must be dense; callbacks are not supported. The Jacobians of
#' the residual come from the System the problem was built from
#' (symbolically) or from the problem's Jacobian and central differences;
#' `du(t)` is the derivative of the forward interpolant (by central
#' differences of the step interpolant).
#'
#' @inheritParams adjoint_sensitivities
#' @param sol The dense solution of a `DAEProblem`.
#' @param sensealg [InterpolatingAdjoint()] (without checkpointing),
#'   [QuadratureAdjoint()], [GaussAdjoint()], [GaussKronrodAdjoint()] or
#'   [SundialsAdjoint()] (solved as `InterpolatingAdjoint`).
#' @param alg The algorithm of the backward solve; the forward algorithm when
#'   `NULL`.
#' @param callback Must be `NULL`: callbacks are not supported.
#' @return An `ODEProblem` or `DAEProblem` with elements `alg` and `finish`
#'   (a function of the backward solution returning `list(du0, dp)`).
#' @examples
#' f <- function(du, u, p, t) c(du[1] + p[1] * u[1] - u[2], u[2] - p[2] * u[1])
#' prob <- DAEProblem(f, du0 = c(-0.5, 0), u0 = c(1, 0.5), tspan = c(0, 1), p = c(1, 0.5),
#'                    differential_vars = c(TRUE, FALSE))
#' sol <- solve(prob, DFBDF(), abstol = 1e-10, reltol = 1e-10)
#' # G = u1(1) = exp(-(p1 - p2)): dG/dp = (-1, 1) exp(-0.5)
#' adjoint_sensitivities(sol, DFBDF(), t = 1, dgdu_discrete = function(u, p, t, i) c(1, 0),
#'                       abstol = 1e-10, reltol = 1e-10)$dp / exp(-0.5)
#' @export
DAEAdjointProblem <- function(sol, sensealg = InterpolatingAdjoint(), alg = NULL, t = NULL,
                              dgdu_discrete = NULL, dgdp_discrete = NULL, dgdu_continuous = NULL,
                              dgdp_continuous = NULL, g = NULL, callback = NULL, no_start = FALSE,
                              abstol = 1e-6, reltol = 1e-3) {
  env <- .sol_env(sol)
  prob <- env$prob
  if (!inherits(prob, "DAEProblem")) stop("DAEAdjointProblem() needs the solution of a DAEProblem", call. = FALSE)
  if (is.null(dgdu_discrete) && is.null(dgdu_continuous) && is.null(g))
    stop("Either `dgdu_discrete`, `dgdu_continuous`, or `g` must be specified.", call. = FALSE)
  if (!is.null(t) && is.null(dgdu_discrete) && is.null(dgdp_discrete))
    stop("It looks like you're using the direct `adjoint_sensitivities` interface with a discrete cost ",
         "function but no specified `dgdu_discrete` or `dgdp_discrete`.", call. = FALSE)
  sundials <- inherits(sensealg, "SundialsAdjoint")
  method <- if (sundials) "InterpolatingAdjoint" else class(sensealg)[1L]
  if (!method %in% c("InterpolatingAdjoint", "QuadratureAdjoint", "GaussAdjoint", "GaussKronrodAdjoint"))
    stop("`", method, "` is not currently supported for adjoint sensitivity analysis of fully implicit ",
         "`DAEProblem`s. Supported methods are `InterpolatingAdjoint`, `QuadratureAdjoint`, `GaussAdjoint`, ",
         "and `GaussKronrodAdjoint`. Alternatively, express the DAE in mass-matrix form as an `ODEProblem`, ",
         "for which all continuous adjoint methods are supported.", call. = FALSE)
  if (!sundials && isTRUE(sensealg$checkpointing))
    stop("Checkpointing is not currently supported for DAEProblem adjoints. Use a dense forward solution with ",
         "`InterpolatingAdjoint(checkpointing = false)`.", call. = FALSE)
  cbs <- .as_callback_set(callback)
  fcb <- .as_callback_set(prob$kwargs$callback)
  if (length(cbs$continuous) + length(cbs$discrete) + length(fcb$continuous) + length(fcb$discrete))
    stop("Callbacks are not currently supported for DAEProblem adjoints.", call. = FALSE)
  withgrad <- method == "InterpolatingAdjoint"
  alg <- alg %||% env$alg %||% DFBDF()
  fns <- .sens_dae_fns(prob, sol, sensealg)
  tr <- fns$tr
  n <- fns$n
  np <- fns$np
  p <- prob$p
  t0 <- prob$tspan[1L]
  tend <- if (identical(env$retcode, "Terminated")) env$t[length(env$t)] else prob$tspan[2L]
  st <- .sens_dae_structure(fns, p, tend, prob$differential_vars)
  if (!withgrad && st$index2)
    stop("`", method, "` does not support Hessenberg index-2 DAEs: the index-2 adjoint multiplier is ",
         "impulsive at discrete cost times and its parameter-gradient contribution is concentrated at the ",
         "quadrature interval endpoints, where it is not captured by a quadrature over the adjoint ",
         "interpolant. Use `InterpolatingAdjoint`, whose state-augmented gradient accumulation integrates ",
         "it correctly.", call. = FALSE)
  tol <- function(x) 1e-12 * max(1, abs(x))
  dir <- if (tend < t0) -1 else 1
  inside <- function(x) dir * (x - t0) > tol(t0) & dir * (tend - x) > tol(tend)

  # costs
  cont <- !is.null(g) || !is.null(dgdu_continuous) || !is.null(dgdp_continuous)
  gu_c <- dgdu_continuous %||% if (!is.null(g)) function(u, p, tt) .sens_fd_grad(function(z) g(z, p, tt), u)
  gp_c <- dgdp_continuous %||% if (!is.null(g) && is.null(dgdu_continuous) && np)
    function(u, p, tt) .sens_fd_grad(function(q) g(u, .sens_pvec(p, q), tt), as.numeric(p))
  ts <- if (is.null(t)) numeric() else as.numeric(t)
  ord <- order(ts)
  cost <- list(t = ts[ord], i = ord)

  # z = (w, lambda, grad)
  iw <- seq_len(n)
  il <- n + seq_len(n)
  ig <- if (withgrad) 2L * n + seq_len(np) else integer()
  N <- 2L * n + length(ig)
  dv <- st$diff
  av <- st$alg
  acc <- new.env(parent = emptyenv())
  acc$mu <- numeric(np)    # quadrature parameter gradient (Gauss methods) and discrete dgdp
  acc$corr <- numeric(np)  # corrections of the consistent jumps

  rhs <- function(z, p, tt) {
    y <- tr$at(tt)
    dy <- tr$dy(tt)
    l <- z[il]
    w <- z[iw]
    vu <- as.numeric(crossprod(fns$Ju(dy, y, p, tt), l))
    vdu <- as.numeric(crossprod(fns$Jdu(dy, y, p, tt), l))
    dgrad <- if (withgrad && np) as.numeric(crossprod(fns$Jp(dy, y, p, tt), l))
    if (cont) {
      if (!is.null(gu_c)) vu <- vu - as.numeric(gu_c(y, p, tt))
      if (withgrad && np && !is.null(gp_c)) dgrad <- dgrad - as.numeric(gp_c(y, p, tt))
    }
    ow <- vu
    ol <- vdu
    ol[av] <- vu[av]
    ow[av] <- w[av]
    ol[dv] <- w[dv] - vdu[dv]
    c(ow, ol, dgrad)
  }
  rhs <- structure(rhs, class = c("DAEAdjointSensitivityFunction", SensitivityFunction, "function"))
  # the adjoint is linear in z
  jac <- function(z, p, tt) {
    y <- tr$at(tt)
    dy <- tr$dy(tt)
    Ju <- fns$Ju(dy, y, p, tt)
    A <- matrix(0, N, N)
    A[iw, il] <- t(Ju)
    A[iw[av], ] <- 0
    A[cbind(iw[av], iw[av])] <- 1
    A[il, il] <- -t(fns$Jdu(dy, y, p, tt))
    A[il[av], il] <- t(Ju)[av, , drop = FALSE]
    A[cbind(il[dv], iw[dv])] <- 1
    if (length(ig)) A[ig, il] <- t(fns$Jp(dy, y, p, tt))
    A
  }
  mmdiag <- numeric(N)
  mmdiag[dv] <- 1
  mmdiag[ig] <- 1

  # one discrete cost sample: the jump of w, consistent with the algebraic
  # constraints of the adjoint (ReverseLossCallback)
  sample_jump <- function(z, p, tt, i) {
    y <- tr$at(tt)
    gu <- if (!is.null(dgdu_discrete)) as.numeric(dgdu_discrete(y, p, tt, i)) else numeric(n)
    gp <- if (!is.null(dgdp_discrete) && np) as.numeric(dgdp_discrete(y, p, tt, i)) else numeric(np)
    if (length(av)) {
      dy <- tr$dy(tt)
      Ju <- fns$Ju(dy, y, p, tt)
      dhdd <- Ju[st$alg_eq, dv, drop = FALSE]
      dla <- if (st$index2) {
        if (any(gu[av] != 0))
          stop("Discrete cost contributions on the algebraic variables of an index-2 DAE are not supported: ",
               "they correspond to derivatives of delta distributions in the adjoint. Reformulate the cost ",
               "in terms of the differential variables.", call. = FALSE)
        Jdu <- fns$Jdu(dy, y, p, tt)
        Nm <- t(solve(Jdu[st$diff_eq, dv, drop = FALSE], Ju[st$diff_eq, av, drop = FALSE]))
        -solve(Nm %*% t(dhdd), Nm %*% gu[dv])
      } else -solve(t(Ju[st$alg_eq, av, drop = FALSE]), gu[av])
      dld <- as.numeric(crossprod(dhdd, dla)) + gu[dv]
      if (np) acc$corr <- acc$corr + as.numeric(crossprod(fns$Jp(dy, y, p, tt)[st$alg_eq, , drop = FALSE], dla))
    } else dld <- gu[dv]
    z[iw[dv]] <- z[iw[dv]] + dld
    if (np && any(gp != 0)) {
      if (withgrad) z[ig] <- z[ig] + gp else acc$mu <- acc$mu + gp
    }
    z
  }
  point <- function(tt, z, p) {
    for (k in rev(which(abs(cost$t - tt) <= tol(tt)))) z <- sample_jump(z, p, tt, cost$i[k])
    z
  }

  # the parameter gradient integrand of the quadrature methods, g_p -
  # lambda' F_p
  mu_integrand <- function(l, p, s) {
    v <- -fns$vjp(l, tr$at(s), p, s)$p
    if (cont && !is.null(gp_c)) v <- v + as.numeric(gp_c(tr$at(s), p, s))
    v
  }
  z0 <- point(tend, numeric(N), p)
  inner <- sort(unique(cost$t))
  inner <- inner[inside(inner)]
  cb_rev <- if (length(inner)) {
    cb <- DiscreteCallback(function(u, tt, integ) .near_any(inner, tt), function(integ) {
      integ$u <- point(integ$t, integ$u, integ$p)
    }, save_positions = c(FALSE, FALSE))
    cb$tstops <- inner
    cb
  }
  gauss_cb <- if (method %in% c("GaussAdjoint", "GaussKronrodAdjoint") && np) {
    DiscreteCallback(function(u, tt, integ) {
      s <- integ$last_step
      a <- s$t0
      b <- s$tend %||% (s$t0 + s$dt)
      integrand <- function(x) mu_integrand(.step_value(s, x)[il], integ$p, x)
      acc$mu <- acc$mu + if (method == "GaussAdjoint") .sens_gauss7(integrand, b, a) else
        .sens_quadgk(integrand, b, a, abstol, reltol)
      FALSE
    }, function(integ) NULL)
  }
  quadrature <- method == "QuadratureAdjoint"
  kw <- list(callback = CallbackSet(gauss_cb, cb_rev), dense = quadrature, save_everystep = quadrature,
             initializealg = if (st$index2) NoInit() else BrownFullBasicInit())
  if (isTRUE(alg$dae)) {
    resid <- structure(function(dz, z, p, tt) mmdiag * dz - rhs(z, p, tt),
                       class = c("DAEAdjointResidual", SensitivityFunction, "function"))
    rjac <- function(dz, z, p, gamma, tt) gamma * diag(mmdiag, N) - jac(z, p, tt)
    # automatic initial step selection does not run backward for DAEs: seed
    # the step with the last forward step
    ft <- env$t
    kw$dt <- if (length(ft) > 1L) sign(t0 - tend) * abs(ft[length(ft)] - ft[length(ft) - 1L]) else
      (t0 - tend) / 1000
    aprob <- do.call(DAEProblem, c(list(DAEFunction(resid, jac = rjac), du0 = ifelse(mmdiag == 1, rhs(z0, p, tend), 0),
                                        u0 = z0, tspan = c(tend, t0), p = p, differential_vars = mmdiag == 1), kw))
  } else {
    aprob <- .new_problem("ODEProblem", f = ODEFunction(rhs, mass_matrix = diag(mmdiag, N), jac_fn = jac),
                          u0 = z0, tspan = c(tend, t0), p = p, sys = NULL, kwargs = kw)
  }
  aprob$alg <- alg
  aprob$finish <- function(adj) {
    aenv <- .sol_env(adj)
    zf <- aenv$U[nrow(aenv$U), ]
    if (!no_start) zf <- point(t0, zf, p)
    du0 <- zf[iw]
    dp <- if (!np) numeric() else if (withgrad) zf[ig] + acc$corr else if (method != "QuadratureAdjoint")
      acc$mu + acc$corr else {
        cuts <- sort(unique(c(t0, cost$t[inside(cost$t)], tend)))
        lam_at <- function(x) .interp_states(aenv, x)[1L, il]
        mu <- acc$mu + acc$corr
        for (k in seq_len(length(cuts) - 1L))
          mu <- mu + .sens_quadgk(function(x) mu_integrand(lam_at(x), p, x), cuts[k], cuts[k + 1L], abstol, reltol)
        mu
      }
    names(du0) <- names(prob$u0)
    if (np) names(dp) <- names(prob$p)
    list(du0 = du0, dp = dp)
  }
  aprob
}
