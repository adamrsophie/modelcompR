# Extrapolation methods of OrdinaryDiffEqExtrapolation 2.6.4: AitkenNeville
# (explicit Euler with the Aitken and Neville scheme on the Romberg
# sequence), ExtrapolationMidpointDeuflhard and
# ExtrapolationMidpointHairerWanner (explicit midpoint rule),
# ImplicitDeuflhardExtrapolation and ImplicitHairerWannerExtrapolation
# (linearly implicit midpoint rule), ImplicitEulerExtrapolation (linearly
# implicit Euler with the Richardson table) and
# ImplicitEulerBarycentricExtrapolation (linearly implicit Euler with
# barycentric weights). The steps follow the out of place ConstantCache
# versions of extrapolation_perform_step.jl, the caches
# extrapolation_caches.jl and the order and step size control
# controllers.jl (ExtrapolationController); the Float64 coefficient tables
# are in R/solvers-extrapolation-tables.R.
#
# The order and step size control runs through the "Extrapolation"
# controller kind of the driver (R/solvers-ode.R): the hooks
# extrap_stepsize (stepsize_controller!, which returns 0 here),
# extrap_accept (step_accept_controller!) and extrap_reject
# (step_reject_controller!). The internal controller
# (stepsize_controller_internal!) runs inside the step and fills the
# scaling factors Q of the cache, which is an environment kept across
# rejected steps, as the mutable ConstantCache of Julia. AitkenNeville uses
# the PIController in Julia (a FIXME upstream) and calls it inside the step
# for every candidate order; its PI state (q11, errold, beta1) lives in the
# cache and its hooks are the PI controller of OrdinaryDiffEqCore.
#
# Arithmetic as in Julia: the Euler substeps of AitkenNeville are fused
# (@muladd under @..), the midpoint and implicit substeps, the Richardson
# tables and the barycentric sums are not; sum() of the weighted internal
# discretizations is a left fold; the error estimate is calculate_residuals
# with ODE_DEFAULT_NORM (.ns_err of R/solvers-nordsieck.R); W \ b is
# LinearAlgebra's backslash on W = dt * J - I (.ojl_backslash).
#
# Reproduced upstream behavior (it changes the order and step size path,
# never an accepted solution):
#   * AitkenNeville keeps the PIController: the candidate loop leaves errold
#     at the value of the last candidate and selects the candidate with the
#     least work 2^(i - 1) / dtpropose, a signed quantity, so a backward
#     integration selects the most expensive order.
#   * The Hairer and Wanner order selection compares work[n_curr - 1] and
#     work[n_curr] (the orders n_curr - 2 and n_curr - 1).
#   * ImplicitEulerExtrapolation never resets its scaling factors Q and
#     replaces them by |dt| / (|dt| / Q) at every accepted step.
#   * The divergence monitor of ImplicitDeuflhardExtrapolation compares
#     every step of the internal discretization with the first one.
#   * The "orders changed" warning of the constructors compares a tuple with
#     itself and never fires; the orders are adjusted silently.
#
# threading: R runs the internal discretizations one after another; the
# values equal those of the threaded code of Julia. Where the threaded code
# counts differently (the evaluations of ExtrapolationMidpointDeuflhard and
# ExtrapolationMidpointHairerWanner, the evaluations, solves and second
# Jacobian of ImplicitEulerExtrapolation) the counts follow it. The implicit
# column methods run the columns in the chunks of Julia's threaded loop, where
# a divergence ends only its own chunk, so a failed step counts the same work.
# The threaded divergence monitor of ImplicitDeuflhardExtrapolation with the
# Romberg sequence compares the first component only, as in Julia.

# ---------------------------------------------------------------------------
# Helpers

# eltype(uprev)(scalar) * sum(broadcast(*, T[idx], w)): a left fold of the
# weighted vectors, then the scaling factor.
.extrap_combine <- function(Tl, idx, w, s) {
  acc <- Tl[[idx[1L]]] * w[1L]
  for (k in seq_along(idx)[-1L]) acc <- acc + Tl[[idx[k]]] * w[k]
  s * acc
}

# W = h * J - I (the identity mass matrix of the extrapolation methods)
.extrap_W <- function(J, h) {
  W <- h * J
  diag(W) <- diag(W) - 1
  W
}

.extrap_solve <- function(W, b) {
  if (nrow(W) <= 16L) .ojl_backslash(W, b) else as.numeric(solve(W, b))
}

# timedepentdtmin of DiffEqBase
.extrap_dtmin <- function(t, dtmin) abs(max(.ojl_ulp(t), dtmin))

# internalnorm(reltol, t)
.extrap_tol <- function(reltol) if (length(reltol) == 1L) abs(reltol) else .ns_norm(reltol)

# The sequence keyword: a string or a symbol; anything else warns and
# becomes "harmonic", as the constructors of algorithms.jl.
.extrap_sequence <- function(sequence, name) {
  s <- if (is.symbol(sequence) || is.character(sequence)) as.character(sequence) else ""
  if (length(s) != 1L || !(s %in% c("harmonic", "romberg", "bulirsch"))) {
    warning("The `sequence` given to the `", name, "` algorithm is not valid: it must match ",
            "`:harmonic`, `:romberg` or `:bulirsch`. Thus it has been changed from :",
            paste(format(sequence), collapse = " "), " to :harmonic", call. = FALSE)
    s <- "harmonic"
  }
  s
}

.extrap_sequence_factor <- function(sequence_factor) {
  sequence_factor <- as.integer(sequence_factor)
  if (sequence_factor %% 2L != 0L) {
    warning("A non-even number cannot be used as sequence factor. Thus it has been changed from ",
            sequence_factor, " to 2", call. = FALSE)
    sequence_factor <- 2L
  }
  sequence_factor
}

# generate_sequence of ImplicitEulerExtrapolation
.extrap_iee_sequence <- function(sequence, max_order) {
  n <- 0:max_order
  switch(sequence,
         harmonic = as.integer(n + 1L),
         romberg = as.integer(2^n),
         bulirsch = as.integer(ifelse(n == 0L, 1, ifelse(n %% 2L == 1L, 2^((n + 1L) %/% 2L), 3 * 2^(n %/% 2L - 1L)))))
}

# ---------------------------------------------------------------------------
# Internal discretizations

# The internal discretization of order i (column i + 1 of T) of the midpoint
# and barycentric families: explicit midpoint ("midpoint"), linearly
# implicit midpoint ("imid", Deuflhard; "ihw", Hairer and Wanner with the
# final smoothing step) and linearly implicit Euler with smoothing ("bary").
# check: the divergence monitor of the first two discretizations.
.extrap_column <- function(ex, ctx, t, uprev, dt, fsal, J, i, check) {
  j_int <- ex$factor * ex$seq[i + 1L]
  dt_int <- dt / j_int
  base <- ex$base
  if (base == "midpoint") {
    u2 <- uprev
    u1 <- u2 + dt_int * fsal
    Ti <- u1
    for (j in seq_len(j_int - 1L) + 1L) {
      Ti <- u2 + 2 * dt_int * ctx$f(u1, ctx$p, t + (j - 1L) * dt_int)
      u2 <- u1
      u1 <- Ti
    }
    return(list(T = Ti, nf = j_int - 1L, nw = 0L, nsolve = 0L, fail = FALSE))
  }
  W <- .extrap_W(J, dt_int)
  u2 <- uprev
  u1 <- u2 + .extrap_solve(W, -(dt_int * fsal))
  diff1 <- u1 - u2
  jmax <- if (base == "imid") j_int else j_int + 1L
  nf <- 0L
  Ti <- u1
  for (j in seq_len(jmax - 1L) + 1L) {
    fj <- ctx$f(u1, ctx$p, t + (j - 1L) * dt_int)
    nf <- nf + 1L
    if (base == "bary") {
      Ti <- u1 + .extrap_solve(W, -(dt_int * fj))
      if (j == jmax) Ti <- 0.25 * (Ti + 2 * u1 + u2)
    } else {
      Ti <- 2 * u1 - u2 + 2 * .extrap_solve(W, -(dt_int * fj - (u1 - u2)))
      if (base == "ihw" && j == jmax) Ti <- 0.5 * (Ti + u2)
    }
    u2 <- u1
    u1 <- Ti
    if (check) {
      diff2 <- u1 - u2
      bad <- if (isTRUE(ex$scalar_check)) abs(diff1[1L]) < abs(0.5 * (diff2[1L] - diff1[1L])) else
        .ns_norm(diff1) < .ns_norm(0.5 * (diff2 - diff1))
      if (isTRUE(bad)) return(list(T = NULL, nf = nf, nw = 1L, nsolve = 0L, fail = TRUE))
    }
    if (base != "imid") diff1 <- u1 - u2
  }
  list(T = Ti, nf = nf, nw = 1L, nsolve = 0L, fail = FALSE)
}

# The linearly implicit Euler discretization of ImplicitEulerExtrapolation
# with m substeps (the first column of row i + 1 of the Richardson table).
.extrap_iee_column <- function(ex, ctx, t, uprev, dt, fsal, J, i, check) {
  m <- ex$seq[i + 1L]
  dt_temp <- dt / m
  W <- .extrap_W(J, dt_temp)
  k_copy <- fsal
  u_tmp <- uprev
  diff1 <- numeric(length(uprev))
  nf <- 0L
  ns <- 0L
  for (j in seq_len(m)) {
    k <- .extrap_solve(W, -(dt_temp * k_copy))
    ns <- ns + 1L
    u_tmp2 <- u_tmp
    u_tmp <- u_tmp + k
    if (check && j >= 2L) {
      diff2 <- u_tmp - u_tmp2
      diff2 <- 0.5 * (diff2 - diff1)
      if (isTRUE(.ns_norm(diff1) < .ns_norm(diff2))) return(list(T = NULL, nf = nf, nw = 1L, nsolve = ns, fail = TRUE))
    }
    diff1 <- u_tmp - u_tmp2
    k_copy <- ctx$f(u_tmp, ctx$p, t + j * dt_temp)
    nf <- nf + 1L
  }
  list(T = u_tmp, nf = nf, nw = 1L, nsolve = ns, fail = FALSE)
}

# Row i of the Richardson table of ImplicitEulerExtrapolation from its first
# column and row i - 1.
.extrap_iee_row <- function(Tm, i, sq) {
  for (j in seq_len(i - 1L) + 1L) {
    r <- sq[i] / sq[i - j + 1L]
    Tm[[i]][[j]] <- (r * Tm[[i]][[j - 1L]] - Tm[[i - 1L]][[j - 1L]]) / (r - 1)
  }
  Tm
}

# ---------------------------------------------------------------------------
# ExtrapolationController

# stepsize_controller_internal!: the scaling factor of order n.
.extrap_qint <- function(ex, EEst, n, co) {
  if (isTRUE(EEst == 0)) return(1 / co$qmax)
  beta1 <- ex$beta1(n)
  g <- .fastpower(ex$gfac, beta1)
  max(1 / co$qmax, min(1 / co$qmin, .fastpower(EEst, beta1) / g))
}

# stepsize_predictor! of the Deuflhard family: the scaling factor of order
# n_new from the error estimate of the current order.
.extrap_predict <- function(st, n_new, co) {
  EEst <- st$EEst
  if (isTRUE(EEst == 0)) return(1 / co$qmax)
  mo <- st$ex$min_order
  s_curr <- st$stage[st$n_curr - mo + 1L]
  s_new <- st$stage[n_new - mo + 1L]
  beta1 <- 1 / (2 * st$n_curr + 1)
  g <- .fastpower(0.25, beta1)
  qtmp <- EEst * .fastpower(.fastpower(st$tolr, 1.0 - s_curr / s_new), beta1) / g
  max(1 / co$qmax, min(1 / co$qmin, qtmp))
}

.extrap_stepsize <- function(cache, EEst, co, first) 0

# step_accept_controller! of the Deuflhard family (dt is |integrator.dt|)
.extrap_deuf_accept <- function(cache, q, dt, co) {
  st <- cache
  mo <- st$ex$min_order
  n_curr <- st$n_curr
  s <- st$stage
  dtmin <- .extrap_dtmin(st$t_new, st$dtmin)
  dtmax <- abs(st$dtmax)
  m <- n_curr - mo + 1L
  dt_new <- pmax(dtmin, pmin(dtmax, abs(dt / st$Q[seq_len(m)])))
  work <- s[seq_len(m)] / dt_new
  n_new <- which.min(work) + mo - 1L
  if (n_new == n_curr && n_curr < min(st$ex$max_order, st$n_old + 1L)) {
    st$Q[m + 1L] <- .extrap_predict(st, n_new + 1L, co)
    dt_new[m + 1L] <- max(dtmin, min(dtmax, abs(dt / st$Q[m + 1L])))
    if (work[m] > s[m + 1L] / dt_new[m + 1L]) n_new <- n_new + 1L
  }
  st$n_curr <- n_new
  dt_new[n_new - mo + 1L]
}

# step_reject_controller! of the Deuflhard family
.extrap_deuf_reject <- function(cache, EEst, dt, co) {
  st <- cache
  mo <- st$ex$min_order
  if (st$n_curr < st$n_old) st$Q[st$n_old - mo + 1L] <- .extrap_predict(st, st$n_old, co)
  st$n_curr <- st$n_old
  dt_red <- dt / st$Q[st$n_old - mo + 1L]
  max(.extrap_dtmin(st$t_old, st$dtmin), min(abs(st$dtmax), abs(dt_red)))
}

# step_accept_controller! of the Hairer and Wanner family
.extrap_hw_accept <- function(cache, q, dt, co) {
  st <- cache
  mo <- st$ex$min_order
  maxo <- st$ex$max_order
  n_curr <- st$n_curr
  n_old <- st$n_old
  s <- st$stage
  sigma <- 0.9
  lo <- min(n_old, n_curr) - 1L
  hi <- max(n_curr, n_old) + 1L
  dt_new <- numeric(length(st$Q))
  absdt <- abs(dt)
  Qd <- absdt / st$Q
  dt_new[lo:hi] <- Qd[lo:hi]
  st$Q <- absdt / Qd
  dtmin <- .extrap_dtmin(st$t_new, st$dtmin)
  work <- numeric(length(st$Q))
  work[lo:hi] <- s[lo:hi] / dt_new[lo:hi]
  n_new <- n_old
  if (n_curr == mo) {
    n_new <- mo + 1L
  } else if (n_curr <= n_old) {
    if (work[n_curr - 1L] < sigma * work[n_curr]) {
      n_new <- max(n_curr - 1L, n_old - 1L, mo + 1L)
    } else if (work[n_curr] < sigma * work[n_curr - 1L]) {
      n_new <- min(n_curr + 1L, maxo - 1L)
    } else {
      n_new <- n_curr
    }
  } else {
    if (work[n_old] < sigma * work[n_old + 1L]) n_new <- max(n_old - 1L, mo + 1L)
    if (work[n_curr + 1L] < sigma * work[n_new + 1L]) n_new <- min(n_new + 1L, maxo - 1L)
  }
  st$n_curr <- n_new
  if (n_new == n_curr + 1L) {
    dt_new[n_new + 1L] <- s[n_curr + 2L] / s[n_curr + 1L] * dt_new[n_curr + 1L]
    dt_new[n_new + 1L] <- max(dtmin, min(abs(st$dtmax), abs(dt_new[n_new + 1L])))
  }
  dt_new[n_new + 1L]
}

# step_reject_controller! of the Hairer and Wanner family
.extrap_hw_reject <- function(cache, EEst, dt, co) {
  st <- cache
  n_red <- st$n_old
  if (st$n_curr == st$n_old - 1L) n_red <- max(st$ex$min_order + 1L, st$n_old - 1L)
  st$n_curr <- n_red
  dt_red <- dt / st$Q[n_red + 1L]
  max(.extrap_dtmin(st$t_old, st$dtmin), min(abs(st$dtmax), abs(dt_red)))
}

# ---------------------------------------------------------------------------
# Steps of the midpoint, barycentric and implicit Euler families

.extrap_init <- function(ctx, t, u, ex) {
  st <- new.env(parent = emptyenv())
  if (ex$base != "ieuler" && ex$max_order > 15L)
    stop("max_order > 15 not allowed for Float32 or Float64 with this algorithm. That's a bad idea.",
         call. = FALSE)
  if (ex$base == "ieuler") {
    ex$seq <- .extrap_iee_sequence(ex$sequence, ex$max_order)
  } else {
    tab <- .EXTRAP_TABLES[[if (ex$base == "bary") "barycentric" else "midpoint"]][[ex$sequence]]
    ex$seq <- tab$seq
    ex$tab <- tab
  }
  S <- cumsum(as.numeric(ex$seq))
  mo <- ex$min_order
  maxo <- ex$max_order
  st$stage <- switch(ex$stage,
    deuf = { n <- seq_len(maxo - mo + 1L); ex$factor * S[mo + n] - mo - n + 3 - ex$factor },
    ideuf = { n <- seq_len(maxo - mo + 1L); 4 * S[mo + n] - mo - n - 1 + length(u) },
    hw = { n <- seq_len(maxo + 1L); ex$factor * S[n] - n + 3 - ex$factor },
    ihw = { n <- seq_len(maxo + 1L); 8 * S[n] - n + 3 },
    iee = { n <- seq_len(maxo + 1L); 2 * S[n] - n + 7 },
    bary = { n <- seq_len(maxo + 1L); 2 * ex$factor * S[n] - n + 7 })
  st$ex <- ex
  st$fsal <- ctx$f(u, ctx$p, t)
  st$nf0 <- 0L
  st$n_curr <- ex$init_order
  st$n_old <- ex$init_order
  st$Q <- numeric(if (ex$family == "deuf") maxo - mo + 1L else maxo + 1L)
  st$EEst <- 1
  st$dtmin <- ctx$dtmin %||% 0
  st$dtmax <- ctx$dtmax %||% Inf
  st$tolr <- .extrap_tol(ctx$reltol)
  st
}

.extrap_step <- function(ctx, t, u, dt, cache) {
  st <- cache
  ex <- st$ex
  co <- ctx$co
  adaptive <- isTRUE(ctx$adaptive)
  uprev <- u
  fsal <- st$fsal
  mo <- ex$min_order
  hw <- ex$family == "hw"
  iee <- ex$base == "ieuler"
  n_curr <- st$n_curr
  if (!iee) st$Q[] <- 0
  if (adaptive) {
    if (hw) {
      if (!(mo + 1L <= n_curr && n_curr <= ex$max_order - 1L))
        stop("Something went wrong while setting up the order window: ", n_curr, " not in [",
             mo + 1L, ",", ex$max_order - 1L, "]. Please report this error", call. = FALSE)
      win_min <- n_curr - 1L
      win_max <- n_curr + 1L
    } else {
      win_min <- max(mo, n_curr - 1L)
      win_max <- min(ex$max_order, n_curr + 1L)
    }
    st$n_old <- n_curr
    n_curr <- win_min
  }
  st$t_old <- t
  st$t_new <- t + dt
  nf <- 0L; njacs <- 0L; nw <- 0L; nsolve <- 0L
  J <- NULL
  if (ex$implicit) {
    jj <- .calc_J(ctx, uprev, t)
    J <- jj$J
    nf <- nf + jj$nf
    njacs <- 1L
    if (iee && ex$threaded) {
      # the threaded branch evaluates the Jacobian a second time
      nf <- nf + jj$nf
      njacs <- 2L
    }
  }
  column <- if (iee) .extrap_iee_column else .extrap_column
  Tl <- vector("list", ex$max_order + 1L)
  cf <- 0L; cw <- 0L; cs <- 0L
  # the work split of Julia's threaded loop: a divergence ends only the
  # current chunk (the return leaves one iteration), the others still run;
  # this decides the counts of a failed step
  chunks <- if (!ex$threaded || iee) list(0:n_curr) else if (ex$sequence == "romberg")
    list(seq_len(n_curr) - 1L, n_curr) else
    lapply(0:(n_curr %/% 2L), function(i) unique(c(i, n_curr - i)))
  failed <- FALSE
  for (ch in chunks) for (i in ch) {
    r <- column(ex, ctx, t, uprev, dt, fsal, J, i, check = i <= 1L)
    cf <- cf + r$nf; cw <- cw + r$nw; cs <- cs + r$nsolve
    if (r$fail) { failed <- TRUE; break }
    Tl[[i + 1L]] <- r$T
  }
  if (failed) {
    # force_stepfail: dt / failfactor, not counted as a rejection
    if (iee && ex$threaded) { cf <- 0L; cw <- 0L; cs <- 0L }
    return(list(u = uprev, fail = TRUE, nreject = 0L, cache = st, k = list(fsal, fsal),
                nf = nf + cf, njacs = njacs, nw = nw + cw, nsolve = nsolve + cs))
  }
  if (ex$threaded && ex$base == "midpoint") {
    cf <- as.integer(st$stage[if (hw) n_curr + 1L else n_curr - mo + 1L] - 1)
  } else if (ex$threaded && iee) {
    cf <- as.integer(sum(ex$seq[seq_len(n_curr + 1L)]) - 1L)
    cs <- cf
    cw <- n_curr + 1L
  }
  nf <- nf + cf; nw <- nw + cw; nsolve <- nsolve + cs
  if (iee) {
    Tm <- lapply(Tl, function(v) list(v))
    for (i in seq_len(n_curr + 1L)) Tm <- .extrap_iee_row(Tm, i, ex$seq)
    approx <- function(i) list(u = Tm[[i + 1L]][[i + 1L]], ut = Tm[[i + 1L]][[i]])
  } else {
    tab <- ex$tab
    approx <- function(i) list(u = .extrap_combine(Tl, seq_len(i + 1L), tab$w[[i + 1L]], tab$s[i + 1L]),
                               ut = .extrap_combine(Tl, seq_len(i) + 1L, tab$w2[[i]], tab$s2[i]))
  }
  qidx <- function(i) if (hw) i + 1L else i - mo + 1L
  EEst <- st$EEst
  if (adaptive) {
    for (i in (if (hw) win_min - 1L else mo):n_curr) {
      a <- approx(i)
      EEst <- .ns_err(a$u - a$ut, uprev, a$u, ctx)
      st$n_curr <- i
      st$Q[qidx(i)] <- .extrap_qint(ex, EEst, i, co)
    }
    while (n_curr <= win_max) {
      if (isTRUE(EEst <= 1)) break
      pass <- if (hw) {
        n_curr < mo + 1L ||
          EEst <= (if (n_curr + 2L <= win_max + 1L) prod(as.numeric(ex$seq[(n_curr + 2L):(win_max + 1L)])) else 1)
      } else {
        # the convergence monitor; Deuflhard's tolerance ||reltol||, the
        # implicit method takes ||utilde - u|| / EEst
        tol <- if (ex$implicit) .ns_norm(a$ut - a$u) / EEst else st$tolr
        EEst <= tol^(st$stage[n_curr - mo + 1L] / st$stage[win_max - mo + 1L] - 1)
      }
      if (!isTRUE(pass)) break
      n_curr <- n_curr + 1L
      st$n_curr <- n_curr
      r <- column(ex, ctx, t, uprev, dt, fsal, J, n_curr, check = FALSE)
      nf <- nf + r$nf; nw <- nw + r$nw; nsolve <- nsolve + r$nsolve
      Tl[[n_curr + 1L]] <- r$T
      if (iee) {
        Tm[[n_curr + 1L]] <- list(r$T)
        Tm <- .extrap_iee_row(Tm, n_curr + 1L, ex$seq)
      }
      a <- approx(n_curr)
      EEst <- .ns_err(a$u - a$ut, uprev, a$u, ctx)
      st$Q[qidx(n_curr)] <- .extrap_qint(ex, EEst, n_curr, co)
    }
    st$EEst <- EEst
    unew <- a$u
  } else {
    unew <- if (iee) Tm[[n_curr + 1L]][[n_curr + 1L]] else
      .extrap_combine(Tl, seq_len(n_curr + 1L), ex$tab$w[[n_curr + 1L]], ex$tab$s[n_curr + 1L])
  }
  fnew <- ctx$f(unew, ctx$p, t + dt)
  if (ex$count_final) nf <- nf + 1L
  if (!adaptive || isTRUE(EEst <= 1)) st$fsal <- fnew
  list(u = unew, EEst = if (adaptive) EEst, k = list(fsal, fnew), cache = st,
       nf = nf, njacs = njacs, nw = nw, nsolve = nsolve)
}

# ---------------------------------------------------------------------------
# AitkenNeville

.aitken_init <- function(ctx, t, u, ex) {
  st <- new.env(parent = emptyenv())
  st$ex <- ex
  st$fsal <- ctx$f(u, ctx$p, t)
  st$nf0 <- 1L
  st$step_no <- 1L
  st$cur_order <- max(ex$init_order, ex$min_order)
  # PIControllerCache: q11, errold (qoldinit) and the integrator's EEst; the
  # controller's beta1, which the step changes, is set on the first step
  st$q11 <- 1
  st$errold <- 1e-4
  st$EEst <- 1
  st$beta1 <- NULL
  st
}

# stepsize_controller! of the PIController
.aitken_q <- function(cache, EEst, co, first) {
  st <- cache
  if (is.null(st$beta1)) st$beta1 <- co$beta1
  qmax <- if (first) co$qmax_first else co$qmax
  if (isTRUE(EEst == 0)) return(1 / qmax)
  q11 <- .fastpower(EEst, st$beta1)
  st$q11 <- q11
  # Julia writes (q11 / fastpower(errold, beta2)) / gamma with the second
  # division under @fastmath, which LLVM reassociates to q11 / (p * gamma)
  min(1 / co$qmin, max(1 / qmax, q11 / (.fastpower(st$errold, co$beta2) * co$gamma)))
}

# step_accept_controller! of the PIController (dt signed in the step, |dt|
# from the driver)
.aitken_accept_dt <- function(st, q, EEst, dt, co) {
  if (co$qsteady_min <= q && q <= co$qsteady_max) q <- 1
  st$errold <- max(EEst, 1e-4)
  dt / q
}

.aitken_accept <- function(cache, q, dt, co) .aitken_accept_dt(cache, q, cache$EEst, dt, co)

.aitken_reject <- function(cache, EEst, dt, co) dt / min(1 / co$qmin, cache$q11 / co$gamma)

.aitken_step <- function(ctx, t, u, dt, cache) {
  st <- cache
  ex <- st$ex
  co <- ctx$co
  if (is.null(st$beta1)) st$beta1 <- co$beta1
  uprev <- u
  fsal <- st$fsal
  maxo <- min(ex$max_order, st$cur_order + 1L)
  Tm <- vector("list", maxo)
  nf <- 0L
  for (i in seq_len(maxo)) {
    m <- 2^(i - 1L)
    dt_temp <- dt / m
    v <- .fma(dt_temp, fsal, uprev)
    k <- ctx$f(v, ctx$p, t + dt_temp)
    nf <- nf + 1L
    for (j in seq_len(m - 1) + 1) {
      v <- .fma(dt_temp, k, v)
      k <- ctx$f(v, ctx$p, t + j * dt_temp)
      nf <- nf + 1L
    }
    Tm[[i]] <- list(v)
  }
  # Richardson extrapolation (the Aitken and Neville scheme)
  tmp <- 1
  for (j in seq_len(maxo - 1L) + 1L) {
    tmp <- tmp * 2
    for (i in j:maxo) Tm[[i]][[j]] <- (tmp * Tm[[i]][[j - 1L]] - Tm[[i - 1L]][[j - 1L]]) / (tmp - 1)
  }
  disco <- ctx$disco_checkpoint
  disco_used <- FALSE
  if (isTRUE(ctx$adaptive)) {
    minimum_work <- Inf
    range_start <- if (st$step_no == 1L) 2L else max(2L, st$cur_order - 1L)
    first <- ctx$success_iter == 0L
    for (i in seq_len(max(0L, maxo - range_start + 1L)) + range_start - 1L) {
      A <- 2^(i - 1L)
      EEst <- .ns_err(Tm[[i]][[i]] - Tm[[i]][[i - 1L]], uprev, Tm[[i]][[i]], ctx)
      beta1 <- st$beta1
      qold <- st$q11
      st$beta1 <- 1 / (i + 1)
      dtpropose <- .aitken_accept_dt(st, .aitken_q(st, EEst, co, first), EEst, dt, co)
      if (!is.null(disco)) {
        # handle_disco_accept! of the first candidate uses up the pending
        # discontinuity step (t is still the start of the step)
        dtpropose <- min((disco - t) / 4, dtpropose)
        disco <- NULL
        disco_used <- TRUE
      }
      st$beta1 <- beta1
      st$q11 <- qold
      work <- A / dtpropose
      if (isTRUE(work < minimum_work)) {
        st$beta1 <- 1 / (i + 1)
        st$cur_order <- i
        minimum_work <- work
        st$EEst <- EEst
      }
    }
  }
  st$step_no <- st$step_no + 1L
  unew <- Tm[[st$cur_order]][[st$cur_order]]
  fnew <- ctx$f(unew, ctx$p, t + dt)
  nf <- nf + 1L
  if (!isTRUE(ctx$adaptive) || isTRUE(st$EEst <= 1)) st$fsal <- fnew
  list(u = unew, EEst = if (isTRUE(ctx$adaptive)) st$EEst, k = list(fsal, fnew), cache = st, nf = nf,
       disco_used = disco_used)
}

# ---------------------------------------------------------------------------
# Constructors

.extrap_alg <- function(name, order, ex, implicit, autodiff = NULL) {
  deuf <- ex$family == "deuf"
  alg <- .new_alg(name, order, step = .extrap_step, init = function(ctx, t, u) .extrap_init(ctx, t, u, ex),
                  stiff = implicit, on_event = .implicit_on_event, controller = "Extrapolation",
                  extrap_stepsize = .extrap_stepsize,
                  extrap_accept = if (deuf) .extrap_deuf_accept else .extrap_hw_accept,
                  extrap_reject = if (deuf) .extrap_deuf_reject else .extrap_hw_reject,
                  min_order = ex$min_order, init_order = ex$init_order, max_order = ex$max_order,
                  sequence = ex$sequence, threading = ex$threading)
  # the implicit midpoint methods subdivide with the fixed factor 4
  if (ex$base %in% c("midpoint", "bary")) alg$sequence_factor <- ex$factor
  if (implicit) .with_autodiff(autodiff, alg) else alg
}

.extrap_spec <- function(family, base, stage, min_order, init_order, max_order, sequence, factor,
                         threading, count_final, beta1, gfac, implicit) {
  list(family = family, base = base, stage = stage, min_order = as.integer(min_order),
       init_order = as.integer(init_order), max_order = as.integer(max_order), sequence = sequence,
       factor = as.integer(factor), threading = threading, threaded = !isFALSE(threading),
       count_final = count_final, beta1 = beta1, gfac = gfac, implicit = implicit)
}

#' Extrapolation methods
#'
#' The extrapolation methods of OrdinaryDiffEqExtrapolation.jl, adaptive in
#' both the step size and the order. Each step computes internal
#' discretizations with the step divided by the terms of a subdividing
#' sequence and extrapolates them to a high order; the difference to the
#' extrapolation of the next lower order estimates the error.
#'
#' * `AitkenNeville()`: explicit Euler on the Romberg sequence with the
#'   Aitken and Neville scheme (as in Julia, it keeps the PI step size
#'   controller of OrdinaryDiffEq).
#' * `ExtrapolationMidpointDeuflhard()` and
#'   `ExtrapolationMidpointHairerWanner()`: explicit midpoint rule with the
#'   order and step size control of Deuflhard and of Hairer and Wanner.
#' * `ImplicitDeuflhardExtrapolation()` and
#'   `ImplicitHairerWannerExtrapolation()`: the linearly implicit midpoint
#'   rule for stiff problems.
#' * `ImplicitEulerExtrapolation()`: linearly implicit Euler with the
#'   Richardson table.
#' * `ImplicitEulerBarycentricExtrapolation()`: linearly implicit Euler
#'   combined with barycentric weights.
#'
#' The implicit methods evaluate one Jacobian per step (the `jac_fn` of the
#' [ODEFunction()] if given, else as chosen by `autodiff`) and solve with
#' `W = h J - I`; they
#' do not support mass matrices, as in Julia. The orders are adjusted to
#' their valid ranges as in Julia (for example the Hairer and Wanner methods
#' need `min_order + 1 <= init_order <= max_order - 1`). `max_order` above
#' 15 is an error for the methods with Float64 coefficient tables (all but
#' `AitkenNeville()` and `ImplicitEulerExtrapolation()`).
#'
#' @param max_order,min_order,init_order The largest, smallest and initial
#'   extrapolation order.
#' @param threading Accepted for compatibility: R computes the internal
#'   discretizations one after another (the result equals Julia's threaded
#'   one; the evaluation counts follow Julia's threaded bookkeeping where it
#'   differs).
#' @param sequence The subdividing sequence: `"harmonic"`, `"romberg"` or
#'   `"bulirsch"` (a string or a symbol); other values warn and use
#'   `"harmonic"`.
#' @param sequence_factor The even factor applied to the subdividing
#'   sequence (odd values warn and use 2).
#' @param autodiff As for [ImplicitEuler()]: how the Jacobian is formed when
#'   the problem has none.
#' @param concrete_jac,linsolve Accepted for compatibility and ignored: the
#'   out of place methods of Julia solve with `W \ b`.
#' @return An ODE algorithm for [solve.ODEProblem()].
#' @examples
#' prob <- ODEProblem(function(u, p, t) c(u[2], -u[1]), c(1, 0), c(0, 1))
#' sol <- solve(prob, ExtrapolationMidpointDeuflhard(), reltol = 1e-10, abstol = 1e-10)
#' sol$u[[length(sol$u)]] - c(cos(1), -sin(1))
#' stiff <- ODEProblem(function(u, p, t) -50 * (u - t * t), 1, c(0, 1))
#' solve(stiff, ImplicitEulerExtrapolation(), reltol = 1e-6)$stats$naccept
#' @name extrapolation_algorithms
NULL

#' @rdname extrapolation_algorithms
#' @export
AitkenNeville <- function(max_order = 10L, min_order = 1L, init_order = 5L, threading = FALSE) {
  ex <- list(min_order = as.integer(min_order), init_order = as.integer(init_order),
             max_order = as.integer(max_order), threading = threading)
  .new_alg("AitkenNeville", max(ex$init_order, ex$min_order), step = .aitken_step,
           init = function(ctx, t, u) .aitken_init(ctx, t, u, ex), on_event = .implicit_on_event,
           controller = "Extrapolation", extrap_stepsize = .aitken_q, extrap_accept = .aitken_accept,
           extrap_reject = .aitken_reject,
           # PIController(alg): beta1 = 7 / (10 alg_order), beta2 = 2 / (5 alg_order)
           # with alg_order = init_order
           beta1 = 7 / (10 * ex$init_order), beta2 = 2 / (5 * ex$init_order), pi_controller = TRUE,
           min_order = ex$min_order, init_order = ex$init_order, max_order = ex$max_order,
           threading = threading)
}

#' @rdname extrapolation_algorithms
#' @export
ExtrapolationMidpointDeuflhard <- function(min_order = 1L, init_order = 5L, max_order = 10L,
                                           sequence = "harmonic", threading = TRUE, sequence_factor = 2L) {
  min_order <- max(1L, min_order)
  init_order <- max(min_order, init_order)
  max_order <- max(init_order, max_order)
  sequence_factor <- .extrap_sequence_factor(sequence_factor)
  sequence <- .extrap_sequence(sequence, "ExtrapolationMidpointDeuflhard")
  ex <- .extrap_spec("deuf", "midpoint", "deuf", min_order, init_order, max_order, sequence, sequence_factor,
                     threading, FALSE, function(n) 1 / (2 * n + 1), 0.25, FALSE)
  .extrap_alg("ExtrapolationMidpointDeuflhard", 2L * (ex$init_order + 1L), ex, FALSE)
}

#' @rdname extrapolation_algorithms
#' @export
ExtrapolationMidpointHairerWanner <- function(min_order = 2L, init_order = 5L, max_order = 10L,
                                              sequence = "harmonic", threading = TRUE, sequence_factor = 2L) {
  min_order <- max(2L, min_order)
  init_order <- max(min_order + 1L, init_order)
  max_order <- max(init_order + 1L, max_order)
  sequence_factor <- .extrap_sequence_factor(sequence_factor)
  sequence <- .extrap_sequence(sequence, "ExtrapolationMidpointHairerWanner")
  ex <- .extrap_spec("hw", "midpoint", "hw", min_order, init_order, max_order, sequence, sequence_factor,
                     threading, TRUE, function(n) 1 / (2 * n + 1), 0.65, FALSE)
  .extrap_alg("ExtrapolationMidpointHairerWanner", 2L * (ex$init_order + 1L), ex, FALSE)
}

#' @rdname extrapolation_algorithms
#' @export
ImplicitDeuflhardExtrapolation <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL,
                                           min_order = 1L, init_order = 5L, max_order = 10L,
                                           sequence = "harmonic", threading = FALSE) {
  min_order <- max(1L, min_order)
  init_order <- max(min_order, init_order)
  max_order <- max(init_order, max_order)
  sequence <- .extrap_sequence(sequence, "ImplicitDeuflhardExtrapolation")
  ex <- .extrap_spec("deuf", "imid", "ideuf", min_order, init_order, max_order, sequence, 4L,
                     threading, TRUE, function(n) 1 / (2 * n + 1), 0.25, TRUE)
  # the threaded branch with the Romberg sequence checks the first component
  ex$scalar_check <- ex$threaded && sequence == "romberg"
  .extrap_alg("ImplicitDeuflhardExtrapolation", 2L * (ex$init_order + 1L), ex, TRUE, autodiff)
}

#' @rdname extrapolation_algorithms
#' @export
ImplicitHairerWannerExtrapolation <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL,
                                              min_order = 2L, init_order = 5L, max_order = 10L,
                                              sequence = "harmonic", threading = FALSE) {
  min_order <- max(2L, min_order)
  init_order <- max(min_order + 1L, init_order)
  max_order <- max(init_order + 1L, max_order)
  sequence <- .extrap_sequence(sequence, "ImplicitHairerWannerExtrapolation")
  ex <- .extrap_spec("hw", "ihw", "ihw", min_order, init_order, max_order, sequence, 4L,
                     threading, FALSE, function(n) 1 / (2 * n + 1), 0.65, TRUE)
  .extrap_alg("ImplicitHairerWannerExtrapolation", 2L * (ex$init_order + 1L), ex, TRUE, autodiff)
}

#' @rdname extrapolation_algorithms
#' @export
ImplicitEulerExtrapolation <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL,
                                       max_order = 12L, min_order = 3L, init_order = 5L,
                                       threading = FALSE, sequence = "harmonic") {
  min_order <- max(3L, min_order)
  init_order <- max(min_order + 1L, init_order)
  max_order <- max(init_order + 1L, max_order)
  sequence <- .extrap_sequence(sequence, "ImplicitEulerExtrapolation")
  ex <- .extrap_spec("hw", "ieuler", "iee", min_order, init_order, max_order, sequence, 1L,
                     threading, TRUE, function(n) 1 / n, 0.65, TRUE)
  .extrap_alg("ImplicitEulerExtrapolation", ex$init_order + 1L, ex, TRUE, autodiff)
}

#' @rdname extrapolation_algorithms
#' @export
ImplicitEulerBarycentricExtrapolation <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL,
                                                  min_order = 3L, init_order = 5L, max_order = 12L,
                                                  sequence = "harmonic", threading = FALSE,
                                                  sequence_factor = 2L) {
  min_order <- max(3L, min_order)
  init_order <- max(min_order + 1L, init_order)
  max_order <- max(init_order + 1L, max_order)
  sequence <- .extrap_sequence(sequence, "ImplicitEulerBarycentricExtrapolation")
  ex <- .extrap_spec("hw", "bary", "bary", min_order, init_order, max_order, sequence,
                     as.integer(sequence_factor), threading, FALSE, function(n) 1 / (n - 1), 0.65, TRUE)
  .extrap_alg("ImplicitEulerBarycentricExtrapolation", ex$init_order, ex, TRUE, autodiff)
}
