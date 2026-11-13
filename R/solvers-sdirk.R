# SDIRK, ESDIRK and IMEX additive Runge-Kutta methods of
# OrdinaryDiffEqSDIRK. Every method runs through the out of place step of
# generic_imex_perform_step.jl (_perform_step_oop! for the :standard
# tableau flavor) on the tableaus of R/solvers-sdirk-tableaus.R:
#   lib/OrdinaryDiffEqSDIRK/src/generic_imex_perform_step.jl  step, FSAL, EEst
#   lib/OrdinaryDiffEqSDIRK/src/algorithms.jl                  keywords
#   lib/OrdinaryDiffEqSDIRK/src/alg_utils.jl                   orders, traits
# Each implicit stage solves dt f(tmp + gamma z, p, t + c dt) = M z with
# the Newton iteration of R/solvers-implicit.R (.nlsolve, method "dirk"),
# gamma being the last diagonal entry of the tableau. An IMEX method on a
# SplitODEProblem treats f1 implicitly (with the Jacobian of the split
# function) and f2 explicitly (nf2 counts the f2 calls).

# z / d elementwise, 0 where d is 0 (_mmdiv); d is the diagonal of the mass
# matrix of a method with an explicit first stage, or NULL.
.sdirk_mmdiv <- function(x, d) if (is.null(d)) x else ifelse(d == 0, 0, x / d)

# Sum of w[j] * Z[, j] over j in idx, added left to right to `acc` (the
# unrolled stage sums of the out of place step, zero weights included).
.sdirk_acc <- function(acc, w, Z, idx) {
  for (j in idx) acc <- acc + w[j] * Z[, j]
  acc
}

.sdirk_init <- function(spec) function(ctx, t, u) {
  tab <- spec$tab
  M <- ctx$mass
  if (tab$diag_mm && !is.null(M) && any(M[row(M) != col(M)] != 0))
    stop(spec$name, " only works with diagonal mass matrices. Please choose a solver suitable for your problem (e.g. Rodas5P)",
         call. = FALSE)
  cache <- new.env(parent = emptyenv())
  cache$nl <- .new_nlsolver(spec$nlsolve)
  cache$spec <- spec
  cache$fsal_val <- ctx$f(u, ctx$p, t)
  # integrator.reeval_fsal: .sdirk_on_event evaluates `fsal` again after an
  # event (the active binding records it); after a step limiter the driver
  # sets reeval_fsal and the next step evaluates it (reset_fsal! of
  # update_fsal!, which the last step never reaches)
  cache$fsal_ext <- FALSE
  cache$defer_fsal <- TRUE
  cache$reeval_fsal <- FALSE
  makeActiveBinding("fsal", function(v) {
    if (missing(v)) return(cache$fsal_val)
    cache$fsal_val <- v
    cache$fsal_ext <- TRUE
  }, cache)
  cache$reset_pending <- FALSE
  cache$last_stepfail <- FALSE
  cache$k_acc <- NULL
  cache
}

# After an event the Newton state is kept and the FSAL derivative is
# evaluated at the new state (reset_fsal!).
.sdirk_on_event <- function(ctx, t, u, cache) {
  cache$fsal <- ctx$f(u, ctx$p, t)
  list(cache = cache, nf = 1L)
}

# The step the MaxOrder, VariableOrder and CutoffOrder predictors
# extrapolate: the state at its start and the derivatives at both ends of
# the last accepted step. OrdinaryDiffEqSDIRK reads integrator.uprev2, which
# only methods with alg_extrapolates update (it stays u0 here), and
# integrator.k, which a rejected step overwrites; R uses the accepted step
# (see PORTING.md).
.sdirk_history <- function(ctx, cache) list(y0 = ctx$uprev2, k = cache$k_acc)

# The stage guess of those predictors: the Hermite interpolant of the step
# `h` (.sdirk_history) in power form, extrapolated to t + c dt. `ci` is the
# stage abscissa, `tmp` the explicit part of the stage and `gam` the
# diagonal of the tableau.
.sdirk_hermite_guess <- function(h, ctx, pred, ci, t, dt, u, tmp, gam) {
  k1 <- h$k[[1L]]
  k2 <- h$k[[2L]]
  y0 <- h$y0
  tprev <- ctx$tprev
  th <- (.fma(ci, dt, t) - tprev) / (t - tprev)
  dtp <- t - tprev
  q <- switch(pred, MaxOrder = 3L, CutoffOrder = if (ci <= 1 / 2) 3L else 1L,
              if (th <= 3 / 2) 3L else if (th <= 5 / 2) 2L else 1L)
  # @muladd puts the bare summands first and the product last
  up <- y0
  if (q >= 2L) up <- up + (th * th) * (3 * (u - y0) - (2 * dtp) * k1 - dtp * k2)
  if (q >= 3L) up <- up + (th * th * th) * ((-2 * (u - y0) + dtp * k1) + dtp * k2)
  up <- up + (th * dtp) * k1
  (up - tmp) * (1 / gam)
}

.sdirk_step <- function(ctx, t, u, dt, cache) {
  spec <- cache$spec
  tab <- spec$tab
  s <- tab$s
  A <- tab$A
  cc <- tab$c
  p <- ctx$p
  nl <- cache$nl
  pred <- spec$predictor
  is_imex <- tab$split && inherits(ctx$fobj, "SplitFunction")
  cnt <- list(nf = 0L, nf2 = 0L, njacs = 0L, nsolve = 0L, nw = 0L)
  eig <- NULL
  # update_fsal! of the accepted step before: f is evaluated again after a
  # step limiter and, for a method restricted to diagonal mass matrices, when
  # it is not adaptive, unless an event already did
  reeval <- isTRUE(ctx$derivative_discontinuity) || cache$fsal_ext || cache$reeval_fsal
  if ((cache$reset_pending || cache$reeval_fsal) && !cache$fsal_ext) {
    cache$fsal_val <- ctx$f(u, p, t)
    cnt$nf <- 1L
  }
  cache$reset_pending <- FALSE
  cache$reeval_fsal <- FALSE
  cache$fsal_ext <- FALSE
  fsalfirst <- cache$fsal_val
  # the nonlinear solver sees f1 of an IMEX split and the whole f otherwise;
  # the Jacobian of a split function is that of f1
  cx <- ctx
  if (is_imex) {
    cx$f <- ctx$f1
    cx$jac <- ctx$jac1
    cx$tgrad <- NULL
  } else if (inherits(ctx$fobj, "SplitFunction")) cx$jac <- ctx$jac1 %||% ctx$jac
  gam <- A[s, s]
  first <- TRUE
  fail <- function() {
    cache$last_stepfail <- TRUE
    c(cnt, list(fail = TRUE, u = u, cache = cache, eigen_est = eig))
  }
  stage <- function(z, tmp, ci) {
    r <- .nlsolve(nl, cx, t, dt, z, tmp, gam, 1, .fma(ci, dt, t), u, method = "dirk",
                  firststage = first, eta_first = tab$adaptive, pred_c = ci)
    first <<- FALSE
    for (nm in c("nf", "njacs", "nsolve", "nw")) cnt[[nm]] <<- cnt[[nm]] + r[[nm]]
    if (!is.null(r$eigen_est)) eig <<- r$eigen_est
    r
  }
  hermite <- pred %in% c("MaxOrder", "VariableOrder", "CutoffOrder")
  if (hermite) hh <- .sdirk_history(ctx, cache)
  hermite_zero <- !hermite || ctx$success_iter == 0L || reeval || cache$last_stepfail || is.null(hh$k)
  M <- ctx$mass
  mmd <- if (!is.null(M) && tab$efs) diag(M) else NULL
  Z <- matrix(0, length(u), s)
  K <- if (is_imex) matrix(0, length(u), s)
  tmp <- u
  # stage 1
  if (tab$efs) {
    if (is_imex) {
      Z[, 1L] <- dt * .sdirk_mmdiv(ctx$f1(u, p, t), mmd)
      K[, 1L] <- dt * .sdirk_mmdiv(fsalfirst, mmd) - Z[, 1L]
    } else Z[, 1L] <- dt * .sdirk_mmdiv(fsalfirst, mmd)
  } else {
    z0 <- if (tab$s1x && pred == "Linear") dt * fsalfirst else u * 0
    r <- stage(z0, u, cc[1L])
    if (!r$ok) return(fail())
    Z[, 1L] <- r$z
    if (is_imex) {
      K[, 1L] <- dt * ctx$f2(u + gam * Z[, 1L], p, t)
      cnt$nf2 <- cnt$nf2 + 1L
    }
  }
  # stages 2 to s
  for (i in seq_len(s)[-1L]) {
    tmp <- .sdirk_acc(u, A[i, ], Z, seq_len(i - 1L))
    if (is_imex) tmp <- .sdirk_acc(tmp, tab$AE[i, ], K, seq_len(i - 1L))
    zg <- if (tab$efs && is_imex) {
      Z[, 1L]
    } else if (pred == "Trivial") {
      u * 0
    } else if (pred == "CopyPrev") {
      Z[, i - 1L]
    } else if (pred == "StageExtrap") {
      if (i == 2L) Z[, 1L] else
        Z[, i - 1L] + (Z[, i - 1L] - Z[, i - 2L]) * ((cc[i] - cc[i - 1L]) / (cc[i - 1L] - cc[i - 2L]))
    } else if (hermite) {
      if (hermite_zero) u * 0 else .sdirk_hermite_guess(hh, ctx, pred, cc[i], t, dt, u, tmp, gam)
    } else if (length(tab$csg) && tab$csg[i] != 0) {
      tab$csg[i]
    } else if (length(tab$alpha) && any(tab$alpha[[i]] != 0)) {
      .sdirk_acc(tab$alpha[[i]][1L] * Z[, 1L], tab$alpha[[i]], Z, seq_len(i - 1L)[-1L])
    } else u * 0
    r <- stage(zg, tmp, cc[i])
    if (!r$ok) return(fail())
    Z[, i] <- r$z
    if (s > i && is_imex) {
      K[, i] <- dt * ctx$f2(tmp + gam * Z[, i], p, .fma(tab$ce[i], dt, t))
      cnt$nf2 <- cnt$nf2 + 1L
    }
  }
  # the new state
  if (is_imex) {
    K[, s] <- dt * ctx$f2(tmp + gam * Z[, s], p, .fma(tab$ce[s], dt, t))
    cnt$nf2 <- cnt$nf2 + 1L
    unew <- .sdirk_acc(.sdirk_acc(u, tab$bi, Z, seq_len(s)), tab$be, K, seq_len(s))
  } else if (tab$sa) {
    unew <- tmp + gam * Z[, s]
  } else unew <- .sdirk_acc(u, tab$bi, Z, seq_len(s))
  # the embedded error estimate, smoothed by W = J - M / (gamma dt)
  EEst <- NULL
  if (ctx$adaptive && length(tab$btilde)) {
    est <- .sdirk_acc(tab$btilde[1L] * Z[, 1L], tab$btilde, Z, seq_len(s)[-1L])
    if (is_imex && length(tab$ebtilde)) est <- .sdirk_acc(est, tab$ebtilde, K, seq_len(s))
    if (spec$smooth_est) {
      est <- .wsolve(nl$W, est)
      cnt$nsolve <- cnt$nsolve + 1L
    }
    EEst <- .err_norm(est, u, unew, ctx$abstol, ctx$reltol)
  }
  # the derivative at the end of the step
  if (is_imex) {
    fnew <- ctx$f(unew, p, t + dt)
  } else if (tab$efl) {
    fnew <- ctx$f(unew, p, .fma(tab$fl_c, dt, t))
    cnt$nf <- cnt$nf + 1L
  } else {
    fnew <- Z[, s] / dt
    if (!is.null(mmd)) fnew <- mmd * fnew
  }
  cache$k_last <- list(fsalfirst, fnew)
  if (!ctx$adaptive || isTRUE(EEst <= 1)) {
    cache$fsal_val <- fnew
    cache$k_acc <- cache$k_last
    cache$last_stepfail <- FALSE
    cache$reset_pending <- tab$diag_mm && !ctx$adaptive
  }
  c(cnt, list(u = unew, EEst = EEst, k = list(.mm_deriv(ctx, fsalfirst), .mm_deriv(ctx, fnew)),
              cache = cache, eigen_est = eig))
}

# Predictor or the deprecated `extrapolant` (_resolve_predictor).
.sdirk_predictor <- function(predictor, extrapolant) {
  if (!is.null(extrapolant)) {
    .Deprecated(msg = "The `extrapolant` argument is deprecated; use `predictor` (a value of `Predictor`).")
    if (extrapolant %in% unlist(Predictor)) return(extrapolant)
    predictor <- switch(extrapolant, constant = Predictor$Trivial, linear = Predictor$Linear,
                        interpolant = Predictor$MaxOrder, predictor)
  }
  if (!(length(predictor) == 1L && predictor %in% unlist(Predictor)))
    stop("`predictor` must be a value of `Predictor`, for example Predictor$Linear", call. = FALSE)
  predictor
}

# `smooth_est` is the flag the method honors (_esdirk_smooth_est): FALSE for
# the methods without one.
.sdirk_alg <- function(name, autodiff = NULL, smooth_est = FALSE, predictor = NULL, extrapolant = NULL,
                       step_limiter = NULL, tab = .SDIRK_TABS[[name]], nlsolve = NLNewton()) {
  if (!is.null(step_limiter) && !is.function(step_limiter))
    stop("step_limiter must be a function(u, integ, p, t)", call. = FALSE)
  tab$A <- matrix(tab$Ai, tab$s, tab$s, byrow = TRUE)
  if (length(tab$Ae)) tab$AE <- matrix(tab$Ae, tab$s, tab$s, byrow = TRUE)
  spec <- list(name = name, tab = tab, predictor = .sdirk_predictor(predictor %||% tab$predictor, extrapolant),
               smooth_est = isTRUE(smooth_est), nlsolve = .check_nlsolve(nlsolve))
  alg <- .new_alg(name, tab$order, adaptive = tab$adaptive, step = .sdirk_step, interp = .hermite_dae_interp,
                  init = .sdirk_init(spec), on_event = .sdirk_on_event, stiff = TRUE, mass_matrix = TRUE,
                  qsteady_max = if (tab$adaptive) 6 / 5, on_reject = .nl_on_reject, step_limiter = step_limiter,
                  nlsolve = spec$nlsolve)
  .with_autodiff(autodiff, alg)
}
