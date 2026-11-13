# Stabilized explicit SDE methods of StochasticDiffEqROCK
# (perform_step/SROCK_perform_step.jl, SROCK_utils.jl): Chebyshev type
# stages for the drift, as many as the spectral radius of the drift
# Jacobian requires, followed by a finishing procedure for the noise. The
# coefficient tables are in R/sysdata.rda (tools/stabilized-dump.jl).

#' Stabilized explicit SDE methods
#'
#' Stabilized Runge-Kutta methods for stiff SDEs, named as in
#' StochasticDiffEqROCK. Their drift stages follow Chebyshev type
#' recurrences whose number grows with the spectral radius of the drift
#' Jacobian, estimated every step by a nonlinear power iteration (or given
#' by `eigen_est`), so stiff drifts are integrated without linear solves.
#' All are fixed step methods (pass `dt`) and accept diagonal, scalar and
#' general noise.
#'
#' * `SROCK1()`: the first order S-ROCK method of Abdulle and Cirilli, Ito
#'   (strong order 1 for one dimensional, diagonal or commutative noise) or
#'   Stratonovich.
#' * `SROCK2()`: the weak order 2 S-ROCK2 method of Abdulle, Vilmart and
#'   Zygalakis (Ito).
#' * `SROCKEM()`: the ROCK stabilized Euler-Maruyama method; with
#'   `strong_order_1 = TRUE` a derivative free Milstein correction is added.
#' * `SKSROCK()`: the second kind S-ROCK method of Abdulle, Almuslimani and
#'   Vilmart for stiff Ito problems; `post_processing = TRUE` applies their
#'   post processing on the last step.
#' * `TangXiaoSROCK2()`: the weak order 2 methods of Tang and Xiao, in five
#'   versions (`version_num`) with different stability domains.
#' * `KomBurSROCK2()`: the weak order 2 method of Komori and Burrage for
#'   Stratonovich problems.
#' * `SROCKC2()`: the weak order 2 method of Komori and Burrage for Ito
#'   problems.
#'
#' The methods with weak order 2 (`SROCK2`, `TangXiaoSROCK2`,
#' `KomBurSROCK2`) draw two point random variables for general noise and
#' `TangXiaoSROCK2` draws two random signs every step, with R's random
#' number generator (seed it, or pass `seed` to [solve()]).
#'
#' @param interpretation `"Ito"` or `"Stratonovich"` (`SROCK1()`).
#' @param strong_order_1 Whether `SROCKEM()` adds the Milstein correction.
#' @param post_processing Whether `SKSROCK()` post processes the last step.
#' @param version_num The version (1 to 5) of `TangXiaoSROCK2()`.
#' @param eigen_est A bound on the spectral radius of the drift Jacobian: a
#'   number, or a function `eigen_est(u, p, t)` evaluated every step. By
#'   default it is estimated by power iteration.
#' @return An SDE algorithm.
#' @examples
#' n <- 10
#' lap <- function(u, p, t) (c(0, u[-n]) - 2 * u + c(u[-1], 0)) * (n + 1)^2
#' prob <- SDEProblem(lap, function(u, p, t) 0.1 * u, rep(1, n), c(0, 0.5))
#' sol <- solve(prob, SROCK2(), dt = 0.05, seed = 1)
#' tail(sol$u, 1)
#' sol2 <- solve(prob, SKSROCK(eigen_est = 4 * (n + 1)^2), dt = 0.05, seed = 1)
#' @name sde_stabilized_methods
NULL

# Uniform draws for the random signs and two point variables; one place to
# draw them so they can be replayed in tests.
.srock_unif <- function(n) stats::runif(n)

# init_chi!: Rademacher variables
.srock_chi <- function(m) ifelse(.srock_unif(m) < 0.5, -1, 1)

# The state that persists between steps: the power iteration's last
# eigenvector and estimate. A new cache (after a callback changed the state)
# marks a discontinuity.
.srock_init <- function(ctx, t, u) {
  st <- new.env(parent = emptyenv())
  st$zprev <- u
  st$eig <- 0
  st$disc <- TRUE
  list(st = st)
}

# maxeig! of StochasticDiffEqROCK (every step, from f(uprev)) or the user's
# bound; returns the number of drift evaluations.
.srock_eig <- function(ctx, st, t, u, a) {
  e <- a$eigen_est
  nf <- 0L
  if (!is.null(e)) {
    st$eig <- .as_num(if (is.function(e)) e(u, ctx$p, t) else e)
  } else {
    nf <- 1L + .stab_maxeig(ctx, st, t, u, ctx$f(u, ctx$p, t), FALSE, 50L)
  }
  st$disc <- FALSE
  nf
}

# choose_deg! of the second order methods: the first tabulated degree at
# least mdeg and the start of its recurrence coefficients.
.srock_deg <- function(ms, mdeg) {
  start <- 1L
  for (i in seq_along(ms)) {
    if (ms[i] >= mdeg) return(list(deg_index = i, mdeg = as.integer(ms[i]), start = start))
    start <- start + as.integer(ms[i]) * 2L - 1L
  }
  stop("no tabulated degree for ", mdeg, " stages", call. = FALSE)
}

# choose_deg! of SROCK1 and SROCKEM: the damping of the largest tabulated
# degree not above mdeg.
.srock_eta <- function(tab, mdeg) {
  eta <- 1
  for (i in seq_along(tab$ms)) {
    if (tab$ms[i] <= mdeg) eta <- tab$meta[i] else break
  }
  eta
}

# The Chebyshev parameters of a damped first order recurrence.
.srock_omega <- function(eta, mdeg) {
  w0 <- 1 + (eta / (mdeg^2))
  sw <- sqrt((w0^2) - 1)
  cosh_inv <- log(w0 + sw)
  w1 <- (sw * cosh(mdeg * cosh_inv)) / (mdeg * sinh(mdeg * cosh_inv))
  list(w0 = w0, w1 = w1, cosh_inv = cosh_inv)
}

.srock_g <- function(g, u, p, t, type) {
  G <- g(u, p, t)
  if (type == "general") as.matrix(G) else .as_num(G)
}

.srock_mv <- function(G, v) .as_num(G %*% v)

.srock_alg <- function(name, order, step, interpretation = "Ito", ...) {
  alg <- .new_sde_alg(name, order, FALSE, NULL, interpretation = interpretation, ...)
  alg$init <- .srock_init
  alg$step <- function(ctx, g, u, t, dt, dW, type, delta, cache, dZ = NULL, tspan = NULL, ...)
    step(ctx, g, u, t, dt, dW, type, cache, alg, tspan)
  alg
}

.srock_result <- function(u, cache, nf) list(u = u, err = NULL, k = NULL, cache = cache, nf = nf)

# The damped Chebyshev stages of SROCK1 and SROCKEM from u: returns the last
# two stages, their times and the function evaluations. `hook(i, unew, uim1,
# tim1)` adds the noise of SROCK1 inside the recurrence.
.srock_cheb <- function(f, p, t, u, dt, mdeg, w0, w1, hook = NULL) {
  uim2 <- u
  k <- f(u, p, t)
  Tim2 <- 1
  Tim1 <- w0
  tim1 <- t + dt * (w1 / w0)
  ti <- tim1
  tim2 <- t
  uim1 <- u + (dt * w1 / w0) * k
  unew <- uim1
  for (i in seq.int(2L, length.out = mdeg - 1L)) {
    Ti <- 2 * w0 * Tim1 - Tim2
    mu <- 2 * w1 * (Tim1 / Ti)
    nu <- 2 * w0 * (Tim1 / Ti)
    kap <- -Tim2 / Ti
    k <- f(uim1, p, tim1)
    unew <- dt * mu * k + nu * uim1 + kap * uim2
    if (!is.null(hook)) unew <- hook(i, unew, uim1, tim1)
    ti <- mu * dt + nu * tim1 + kap * tim2
    if (i < mdeg) {
      uim2 <- uim1
      uim1 <- unew
      tim2 <- tim1
      tim1 <- ti
      Tim2 <- Tim1
      Tim1 <- Ti
    }
  }
  list(u = unew, t = ti, nf = mdeg)
}

.srock1_step <- function(ctx, g, u, t, dt, dW, type, cache, a, tspan) {
  st <- cache$st
  p <- ctx$p
  nf <- .srock_eig(ctx, st, t, u, a)
  # The fewest stages whose fitted stability length s^2 (131.97 - 0.45 s) / 197,
  # twice the mean square stability domain (about 0.33 s^2), covers
  # 2 |dt| rho, by binary search. StochasticDiffEqROCK compares the fit with
  # the stage estimate floor(sqrt(2 |dt| rho) + 1) instead, which takes too
  # few stages for stiff problems (see PORTING.md).
  target <- 2 * abs(dt) * abs(st$eig)
  mn <- 3L
  mx <- 200L
  while (mx - mn > 1L) {
    mid <- floor((mn + mx) * 0.5)
    if (mid^2 * (131.97 / 197 - 0.45 * mid / 197) > target) mx <- mid else mn <- mid
  }
  mdeg <- as.integer(if (target > mn^2 * (131.97 / 197 - 0.45 * mn / 197)) mx else mn)
  om <- .srock_omega(.srock_eta(.STAB_TABLES$SROCK1, mdeg), mdeg)
  strat <- a$interpretation == "Stratonovich"
  if (strat) {
    alpha <- cosh(mdeg * om$cosh_inv) / (2 * om$w0 * cosh((mdeg - 1) * om$cosh_inv))
    gam <- 1 / (2 * alpha)
    beta <- -gam
  }
  diag <- type == "diagonal"
  gm2 <- NULL
  hook <- function(i, unew, uim1, tim1) {
    if (strat && i > mdeg - 2L) {
      if (i == mdeg - 1L) {
        gm2 <<- .srock_g(g, uim1, p, tim1, type)
        if (diag) unew + alpha * gm2 * dW else unew + alpha * .gdw(gm2, dW, type)
      } else {
        gm1 <- .srock_g(g, uim1, p, tim1, type)
        if (diag) unew + (beta * gm2 + gam * gm1) * dW
        else unew + (beta * .gdw(gm2, dW, type) + gam * .gdw(gm1, dW, type))
      }
    } else if (!strat && i == mdeg) {
      gm2 <- .srock_g(g, uim1, p, tim1, type)
      if (type == "general") return(unew + .gdw(gm2, dW, type))
      sq <- sqrt(abs(dt))
      gm1 <- .srock_g(g, uim1 + sq * gm2, p, tim1, type)
      unew + (gm2 * dW + (1 / (2 * sq)) * (gm1 - gm2) * (dW^2 - abs(dt)))
    } else unew
  }
  r <- .srock_cheb(ctx$f, p, t, u, dt, mdeg, om$w0, om$w1, hook)
  .srock_result(r$u, cache, nf + r$nf)
}

.srockem_step <- function(ctx, g, u, t, dt, dW, type, cache, a, tspan) {
  st <- cache$st
  p <- ctx$p
  nf <- .srock_eig(ctx, st, t, u, a)
  mdeg <- floor(sqrt(abs(dt) * abs(st$eig) / (if (a$strong_order_1) 0.19 else 0.33)) + 1)
  mdeg <- as.integer(max(3, min(mdeg, 200)))
  om <- .srock_omega(.srock_eta(.STAB_TABLES$SROCKEM, mdeg), mdeg)
  r <- .srock_cheb(ctx$f, p, t, u, dt, mdeg, om$w0, om$w1)
  unew <- r$u
  Gs <- .srock_g(g, unew, p, r$t, type)
  unew <- unew + .gdw(Gs, dW, type)
  if (a$strong_order_1) {
    if (type != "general") {
      uim2 <- 0.5 * Gs * (dW^2 - abs(dt))
      Gs <- .srock_g(g, unew + uim2, p, r$t, type)
      uim1 <- 0.5 * Gs
      Gs <- .srock_g(g, unew - uim2, p, r$t, type)
      uim1 <- uim1 - 0.5 * Gs
      unew <- unew + uim1
    } else {
      m <- length(dW)
      for (i in seq_len(m)) {
        uim2 <- .srock_mv(Gs, 0.5 * (dW * dW[i] - (seq_len(m) == i) * abs(dt)))
        uim1 <- 0.5 * .srock_g(g, unew + uim2, p, r$t, type)[, i]
        uim1 <- uim1 - 0.5 * .srock_g(g, unew - uim2, p, r$t, type)[, i]
        unew <- unew + uim1
      }
    }
  }
  .srock_result(unew, cache, nf + r$nf)
}

.sksrock_step <- function(ctx, g, u, t, dt, dW, type, cache, a, tspan) {
  st <- cache$st
  f <- ctx$f
  p <- ctx$p
  tab <- .STAB_TABLES$SKSROCK
  nf <- .srock_eig(ctx, st, t, u, a)
  eta <- 0.05
  mdeg <- floor(sqrt((abs(dt) * abs(st$eig) + 1.5) / (2 - eta * 4 / 3)) + 1)
  mdeg <- as.integer(max(3, min(mdeg, 200)))
  om <- .srock_omega(eta, mdeg)
  w0 <- om$w0
  w1 <- om$w1
  mu <- w1 / w0
  nu <- mdeg * w1 / 2
  kap <- mdeg * w1 / w0
  Tim2 <- 1
  Tim1 <- w0
  tim1 <- t + dt * (w1 / w0)
  ti <- tim1
  tim2 <- t
  # stage 1
  uu <- .gdw(.srock_g(g, u, p, t, type), dW, type)
  if (a$post_processing) {
    al <- tab$malpha[mdeg - 1L]
    uim1 <- u + nu * uu
    uim2 <- f(u, p, t)
    uim1 <- f(uim1, p, t)
    uim1 <- u + (mu * dt) * uim1 + kap * uu + al * dt * (uim1 - 2 * uim2)
    uim2 <- f(u - nu * uu, p, t)
    uim1 <- uim1 + (al * dt) * uim2
    nf <- nf + 3L
  } else {
    uim2 <- f(u + nu * uu, p, t)
    uim1 <- u + (mu * dt) * uim2 + kap * uu
    nf <- nf + 1L
  }
  uim2 <- u
  for (i in seq.int(2L, length.out = mdeg - 1L)) {
    Ti <- 2 * w0 * Tim1 - Tim2
    mu <- 2 * w1 * (Tim1 / Ti)
    nu <- 2 * w0 * (Tim1 / Ti)
    kap <- -Tim2 / Ti
    uu <- f(uim1, p, tim1)
    uu <- dt * mu * uu + nu * uim1 + kap * uim2
    ti <- mu * dt + nu * tim1 + kap * tim2
    if (i < mdeg) {
      uim2 <- uim1
      uim1 <- uu
      tim2 <- tim1
      tim1 <- ti
      Tim2 <- Tim1
      Tim1 <- Ti
    }
  }
  nf <- nf + mdeg - 1L
  if (a$post_processing && (t + dt >= tspan[2L])) {
    Gs <- .srock_g(g, uu, p, ti, type)
    uim1 <- if (type == "general") Gs[, 1L] else Gs
    winc <- .srock_unif(1L) * 6
    if (winc < 1) uu <- uu - (sqrt(3 * dt) * tab$mc[mdeg - 1L]) * uim1
    else if (winc < 2) uu <- uu + (sqrt(3 * dt) * tab$mc[mdeg - 1L]) * uim1
  }
  .srock_result(uu, cache, nf)
}

# Stages 1 to s - 2 of the second order methods (SROCK2, SROCKC2), scaled
# by alpha: returns the last two stages, their times and the coefficient
# pointer.
.srock2_stages <- function(f, p, t, u, dt, alpha, recf, start, mdeg) {
  mu <- recf[start]
  ti <- t + alpha * dt * mu
  tim1 <- ti
  tim2 <- t
  uim2 <- u
  ui <- f(u, p, t)
  uim1 <- u + alpha * dt * mu * ui
  for (i in seq.int(2L, length.out = mdeg - 1L)) {
    mu <- recf[start + 2L * (i - 2L) + 1L]
    kap <- recf[start + 2L * (i - 2L) + 2L]
    nu <- 1 + kap
    ui <- f(uim1, p, tim1)
    ui <- alpha * dt * mu * ui + nu * uim1 - kap * uim2
    uim2 <- uim1
    uim1 <- ui
    ti <- alpha * dt * mu + nu * tim1 - kap * tim2
    tim2 <- tim1
    tim1 <- ti
  }
  list(uim1 = uim1, uim2 = uim2, tim1 = tim1, tim2 = tim2, nf = mdeg)
}

.srock2_step <- function(ctx, g, u, t, dt, dW, type, cache, a, tspan) {
  st <- cache$st
  f <- ctx$f
  p <- ctx$p
  tab <- .STAB_TABLES$SROCK2
  chi <- if (type == "general") .srock_chi(length(dW))
  nf <- .srock_eig(ctx, st, t, u, a)
  mdeg <- floor(sqrt((2 * abs(dt) * abs(st$eig) + 1.5) / 0.811) + 1)
  d <- .srock_deg(tab$ms, max(3, min(mdeg, 200)) - 2)
  mdeg <- d$mdeg
  di <- d$deg_index
  alpha <- tab$malpha[di]
  sg <- tab$msigma[di]
  sigma <- (1 - alpha) * 0.5 + alpha * sg
  tau <- 0.5 * ((1 - alpha)^2) + 2 * alpha * (1 - alpha) * sg + (alpha^2) * (sg * (sg + tab$mtau[di]))
  sqrt_dt <- sqrt(abs(dt))
  s <- .srock2_stages(f, p, t, u, dt, alpha, tab$recf, d$start, mdeg)
  uim1 <- s$uim1
  uim2 <- s$uim2
  tim1 <- s$tim1
  tim2 <- s$tim2
  # stage s - 1
  r2 <- tab$recf2[(di - 1L) * 4L + 1:4]
  mu <- r2[1L]
  kap <- r2[2L]
  nu <- 1 + kap
  ui <- f(uim1, p, tim1)
  tx <- tim1 + 2 * dt * tau
  ux <- uim1 + 2 * dt * tau * ui
  unew <- uim1 + (2 * sigma - 0.5) * dt * ui
  ui <- alpha * dt * mu * ui + nu * uim1 - kap * uim2
  ti <- alpha * dt * mu + nu * tim1 - kap * tim2
  tim2 <- tim1
  tim1 <- ti
  uim2 <- uim1
  uim1 <- ui
  # stage s
  mu <- r2[3L]
  kap <- r2[4L]
  nu <- 1 + kap
  ui <- f(uim1, p, tim1)
  ui <- alpha * dt * mu * ui + nu * uim1 - kap * uim2
  ti <- alpha * dt * mu + nu * tim1 - kap * tim2
  nf <- nf + s$nf + 2L
  if (type != "general") {
    Gs <- .srock_g(g, uim1, p, tim1, type)
    unew <- unew + Gs * dW
    Gs <- .srock_g(g, ui, p, ti, type)
    ux <- ux + Gs * dW
    ux <- f(ux, p, tx)
    unew <- unew + (0.5 * dt) * ux
    ux <- 0.5 * Gs * (dW^2 - abs(dt))
    unew <- unew + 0.5 * .srock_g(g, ui + ux, p, ti, type)
    unew <- unew - 0.5 * .srock_g(g, ui - ux, p, ti, type)
    ux <- Gs * sqrt_dt
    unew <- unew + 0.25 * dW * (.srock_g(g, ui + ux, p, ti, type) - 2 * Gs)
    unew <- unew + 0.25 * dW * .srock_g(g, ui - ux, p, ti, type)
  } else {
    m <- length(dW)
    idx <- seq_len(m)
    Gs <- .srock_g(g, uim1, p, tim1, type)
    unew <- unew + .srock_mv(Gs, dW)
    Gs <- .srock_g(g, ui, p, ti, type)
    ux <- ux + .srock_mv(Gs, dW)
    ux <- f(ux, p, tx)
    unew <- unew + 0.5 * dt * ux
    for (i in idx) {
      wr <- 0.5 * (dW * dW[i] - (idx == i) * abs(dt) - (idx > i) * abs(dt) * chi + (idx < i) * abs(dt) * chi[i])
      ux <- .srock_mv(Gs, wr)
      unew <- unew + 0.5 * .srock_g(g, ui + ux, p, ti, type)[, i]
      unew <- unew - 0.5 * .srock_g(g, ui - ux, p, ti, type)[, i]
    }
    ux <- sqrt_dt * .srock_mv(Gs, chi)
    unew <- unew + (0.25 * .srock_mv(.srock_g(g, ui + ux, p, ti, type), dW) - 0.5 * .srock_mv(Gs, dW))
    unew <- unew + 0.25 * .srock_mv(.srock_g(g, ui - ux, p, ti, type), dW)
  }
  nf <- nf + 1L
  .srock_result(unew, cache, nf)
}

.tangxiao_step <- function(ctx, g, u, t, dt, dW, type, cache, a, tspan) {
  st <- cache$st
  f <- ctx$f
  p <- ctx$p
  tab <- .STAB_TABLES$TangXiaoSROCK2
  v <- a$version_num
  nhat <- tab$mnhat[v]
  nf <- .srock_eig(ctx, st, t, u, a)
  mdeg <- floor(sqrt((abs(dt) * abs(st$eig) + 1.5) / (if (v <= 2) 0.811 else 0.611)) + 1)
  d <- .srock_deg(tab$ms, max(4, min(mdeg, 200)) - 2)
  mdeg <- d$mdeg
  di <- d$deg_index
  start <- d$start
  start_mcs <- 1L + sum(tab$mnhat[seq_len(v - 1L)])
  alpha <- tab$malpha[v]
  sg <- tab$msigma[di]
  sigma <- (1 - alpha) * 0.5 + alpha * sg
  tau <- 0.5 * ((1 - alpha)^2) + 2 * alpha * (1 - alpha) * sg + (alpha^2) * (sg * (sg + tab$mtau[di]))
  eta1 <- if (.srock_unif(1L) < 0.5) -1 else 1
  eta2 <- if (.srock_unif(1L) < 0.5) -1 else 1
  sqrt_dt <- sqrt(abs(dt))
  n <- length(u)
  U1 <- numeric(n)
  U2 <- numeric(n)
  th1 <- 0
  th2 <- 0
  ti <- tim1 <- tim2 <- tx <- t
  ui <- uim1 <- ux <- uim2 <- u
  unew <- u
  for (i in 0:(mdeg + 1L)) {
    if (i == 1L) {
      mu <- tab$recf[start]
      ti <- tim1 <- t + alpha * dt * mu
      ui <- f(u, p, t)
      uim1 <- u + alpha * dt * mu * ui
    } else if (i > 1L && i <= mdeg) {
      mu <- tab$recf[start + 2L * (i - 2L) + 1L]
      kap <- tab$recf[start + 2L * (i - 2L) + 2L]
      nu <- 1 + kap
      ui <- f(uim1, p, tim1)
      ui <- alpha * dt * mu * ui + nu * uim1 - kap * uim2
      ti <- alpha * dt * mu + nu * tim1 - kap * tim2
    } else if (i == mdeg + 1L) {
      mu <- tab$recf2[(di - 1L) * 4L + 1L]
      kap <- tab$recf2[(di - 1L) * 4L + 2L]
      nu <- 1 + kap
      ui <- f(uim1, p, tim1)
      tx <- tim1 + 2 * tau * dt
      ux <- uim1 + (2 * tau * dt) * ui
      unew <- uim1 + (2 * sigma - 0.5) * dt * ui
      ui <- alpha * dt * mu * ui + nu * uim1 - kap * uim2
      ti <- alpha * dt * mu + nu * tim1 - kap * tim2
    }
    j <- i - mdeg - 1L + nhat
    if (j > 0) {
      j <- j + start_mcs - 1L
      if (i == 0L) {
        uc <- u
        tc <- t
      } else if (i == 1L) {
        uc <- uim1
        tc <- tim1
      } else {
        uc <- ui
        tc <- ti
      }
      U1 <- U1 + tab$c1[j] * uc
      th1 <- th1 + tab$c1[j] * tc
      U2 <- U2 + tab$c2[j] * uc
      th2 <- th2 + tab$c2[j] * tc
    }
    if (i > 1L && i < mdeg + 1L) {
      uim2 <- uim1
      uim1 <- ui
      tim2 <- tim1
      tim1 <- ti
    }
  }
  nf <- nf + mdeg + 1L
  h <- 0.5 * eta1 * sqrt_dt
  if (type == "scalar") {
    Gs <- .srock_g(g, U1, p, th1, type)
    ux <- f(ux + Gs * dW, p, tx)
    unew <- unew + ((0.5 * dt) * ux + Gs * ((dW^2 - abs(dt)) / (eta1 * sqrt_dt) - dW))
    U1 <- U1 - (eta1 * sqrt_dt / 2) * Gs
    U2 <- U2 + (eta1 * sqrt_dt / 2) * Gs
    unew <- unew + .srock_g(g, U2, p, th2, type) * dW
    unew <- unew + .srock_g(g, U1, p, th1, type) * (dW - (dW^2 - abs(dt)) / (eta1 * sqrt_dt))
  } else if (type == "diagonal") {
    Gs <- .srock_g(g, U1, p, th1, type)
    uim1 <- Gs * dW
    U1 <- U1 - h * Gs
    U2 <- U2 + h * Gs
    ux <- f(ux + uim1, p, tx)
    unew <- unew + 0.5 * dt * ux
    unew <- unew + Gs * ((dW^2 - abs(dt)) / (eta1 * sqrt_dt) - dW)
    unew <- unew + .srock_g(g, U2, p, th2, type) * dW
    unew <- unew - .srock_g(g, U1, p, th1, type) * ((dW^2 - abs(dt)) / (eta1 * sqrt_dt) - dW)
  } else {
    m <- length(dW)
    Gs <- .srock_g(g, U1, p, th1, type)
    uim1 <- Gs[, 1L] * dW[1L]
    for (i in seq_len(m)[-1L]) uim1 <- uim1 + Gs[, i] * dW[i]
    ux <- f(ux + uim1, p, tx)
    unew <- unew + ((0.5 * dt) * ux - uim1)
    for (i in seq_len(m)) {
      uim1 <- U1 - h * Gs[, i]
      Gs1 <- .srock_g(g, uim1, p, th1, type)
      unew <- unew + (Gs1[, i] * dW[i] + (Gs[, i] - Gs1[, i]) * ((dW[i]^2 - abs(dt)) / (eta1 * sqrt_dt)))
    }
    # as in StochasticDiffEqROCK, the argument accumulates over i from the
    # last one of the loop above
    for (i in seq_len(m)) {
      for (k in seq_len(m)) {
        wik <- if (i > k) 0.5 * (1 + eta2) * dW[k] else if (i < k) 0.5 * (1 - eta2) * dW[k] else 0.5 * (eta1 * sqrt_dt)
        uim1 <- uim1 + Gs[, k] * wik
      }
      unew <- unew + .srock_g(g, uim1, p, th2, type)[, i] * dW[i]
    }
  }
  nf <- nf + 1L
  .srock_result(unew, cache, nf)
}

.komburk_step <- function(ctx, g, u, t, dt, dW, type, cache, a, tspan) {
  st <- cache$st
  f <- ctx$f
  p <- ctx$p
  tab <- .STAB_TABLES$KomBurSROCK2
  nf <- .srock_eig(ctx, st, t, u, a)
  mdeg <- floor(sqrt((2 * abs(dt) * abs(st$eig) + 1.5) / 0.811) + 1)
  d <- .srock_deg(tab$ms, max(6, min(mdeg, 200)) - 2)
  # mdeg is s in the paper
  mdeg <- d$mdeg + 2L
  start <- d$start
  di <- d$deg_index
  sigma <- tab$msigma[di]
  tau <- tab$mtau[di]
  sqrt_dt <- sqrt(abs(dt))
  gen <- type == "general"
  m <- length(dW)
  chi <- if (gen) .srock_chi(m)
  recf <- tab$recf
  tim2 <- t
  uim2 <- u
  mu <- recf[start]
  ti <- tim1 <- t + dt * mu
  uu <- f(u, p, t)
  uim1 <- u + dt * mu * uu
  for (i in seq.int(2L, length.out = mdeg - 5L)) {
    mu <- recf[start + 2L * (i - 2L) + 1L]
    th <- recf[start + 2L * (i - 2L) + 2L]
    uu <- f(uim1, p, tim1)
    uu <- dt * mu * uu + (1 + th) * uim1 - th * uim2
    ti <- dt * mu + (1 + th) * tim1 - th * tim2
    uim2 <- uim1
    uim1 <- uu
    tim2 <- tim1
    tim1 <- ti
  }
  mus3 <- dt * recf[start + 2L * (mdeg - 5L) + 1L]
  ths3 <- recf[start + 2L * (mdeg - 5L) + 2L]
  mus2 <- dt * recf[start + 2L * (mdeg - 4L) + 1L]
  ths2 <- recf[start + 2L * (mdeg - 4L) + 2L]
  del <- tab$mdelta[(di - 1L) * 8L + 1:8]
  d1 <- dt * del[1L]
  d2 <- dt * del[2L]
  d3 <- dt * del[3L]
  d4 <- del[4L]
  d5 <- del[5L]
  d6 <- del[6L]
  d7 <- del[7L]
  d8 <- del[8L]
  C1 <- (1 + ths2) * mus3
  uu <- uim1 + ths3 * (1 + ths2) * (uim1 - uim2)
  uim1 <- uim1 + ths3 * (uim1 - uim2)
  uim2 <- uu
  ttmp <- tim1 + ths3 * (1 + ths2) * (tim1 - tim2)
  tim1 <- tim1 + ths3 * (tim1 - tim2)
  tim2 <- ttmp
  nf <- nf + mdeg - 4L + 4L
  if (!gen) {
    # stage s - 3
    ys3 <- f(uim1, p, tim1)
    Xs3 <- .srock_g(g, uim1 + mus3 * ys3, p, tim1 + mus3, type)
    uu <- uu + (C1 * ys3 + 1 / 8 * dW * Xs3)
    # stage s - 2
    ys2 <- f(uim1 + mus3 * ys3 + d4 * dW * Xs3, p, tim1 + mus3)
    Xs2 <- .srock_g(g, uim2 + C1 * ys3 + 2 / 3 * dW * Xs3, p, tim2 + C1, type)
    uu <- uu + (mus2 * ys2 + 3 / 8 * dW * Xs2)
    # stage s - 1
    ys1 <- f(uim2 + mus2 * ys2 + C1 * ys3 + d5 * dW * Xs3 + d4 * dW * Xs2, p, tim2 + mus2 + C1)
    Xs1 <- .srock_g(g, uim1 + mus3 * ys3 + d1 * ys2 + 1 / 12 * dW * Xs3 + 1 / 4 * dW * Xs2, p,
                    tim1 + mus3 + d1, type)
    uu <- uu + ((sigma - tau) * dt * ys1 + 3 / 8 * dW * Xs1)
    # stage s
    utmp <- f(uim2 + C1 * ys3 + mus2 * ys2 + sigma * dt * ys1 + d6 * dW * Xs3 + d7 * dW * Xs2 + d8 * dW * Xs1,
              p, tim2 + C1 + mus2 + sigma * dt)
    uu <- uu + (sigma + tau) * dt * utmp
    Xs1 <- .srock_g(g, uim1 + mus3 * ys3 + d2 * ys2 + d3 * ys1 - 5 / 4 * dW * Xs3 + 1 / 4 * dW * Xs2 +
                      2 * dW * Xs1, p, tim1 + mus3 + d2 + d3, type)
    uu <- uu + 1 / 8 * dW * Xs1
  } else {
    # stage s - 3
    ys3 <- f(uim1, p, tim1)
    Xs3 <- .srock_g(g, uim1 + mus3 * ys3, p, tim1 + mus3, type)
    SX3 <- .srock_mv(Xs3, dW)
    uu <- uu + (C1 * ys3 + 1 / 8 * SX3)
    # stage s - 2
    ys2 <- f(uim1 + mus3 * ys3 + d4 * SX3, p, tim1 + mus3)
    Xs2 <- Xs3 * 0
    for (i in seq_len(m))
      Xs2[, i] <- .srock_g(g, uim2 + C1 * ys3 + 2 / 3 * (Xs3[, i] * dW[i]), p, tim2 + C1, type)[, i]
    SX2 <- .srock_mv(Xs2, dW)
    uu <- uu + (mus2 * ys2 + 3 / 8 * SX2)
    # stage s - 1
    ys1 <- f(uim2 + mus2 * ys2 + C1 * ys3 + d5 * SX3 + d4 * SX2, p, tim2 + mus2 + C1)
    Xs1 <- Xs2 * 0
    for (i in seq_len(m))
      Xs1[, i] <- .srock_g(g, uim1 + mus3 * ys3 + d1 * ys2 + 1 / 4 * SX3 + 3 / 4 * SX2 -
                             1 / 6 * (Xs3[, i] * dW[i]) - 1 / 2 * (Xs2[, i] * dW[i]),
                           p, tim1 + mus3 + d1, type)[, i]
    SX1 <- .srock_mv(Xs1, dW)
    uu <- uu + ((sigma - tau) * dt * ys1 + 3 / 8 * SX1)
    # stage s
    utmp <- f(uim2 + C1 * ys3 + mus2 * ys2 + sigma * dt * ys1 + d6 * SX3 + d7 * SX2 + d8 * SX1,
              p, tim2 + C1 + mus2 + sigma * dt)
    uu <- uu + (sigma + tau) * dt * utmp
    S <- numeric(length(u))
    for (i in seq_len(m))
      S <- S + .srock_g(g, uim1 + mus3 * ys3 + d2 * ys2 + d3 * ys1 + 1 / 4 * SX3 + 3 / 4 * SX2 -
                          3 / 2 * (Xs3[, i] * dW[i]) - 1 / 2 * (Xs2[, i] * dW[i]) + 2 * (Xs1[, i] * dW[i]),
                        p, tim1 + mus3 + d2 + d3, type)[, i] * dW[i]
    uu <- uu + 1 / 8 * S
    for (i in seq_len(m)) {
      S <- numeric(length(u))
      for (k in seq_len(m)) {
        if (k != i) S <- S + (if (i > k) -chi[i] * dW[k] else chi[k] * dW[i]) * Xs3[, k]
      }
      ttmp <- tim1 + mus3
      Gm <- .srock_g(g, uim1 + mus3 * ys3 - 1 / 4 * S, p, ttmp, type)
      Gp <- .srock_g(g, uim1 + mus3 * ys3 + 1 / 4 * S, p, ttmp, type)
      # as in StochasticDiffEqROCK, the factor m - 1 compounds over i
      sqrt_dt <- sqrt_dt * (m - 1)
      uu <- uu + sqrt_dt * (Gp[, i] - Gm[, i])
    }
  }
  .srock_result(uu, cache, nf)
}

.srockc2_step <- function(ctx, g, u, t, dt, dW, type, cache, a, tspan) {
  st <- cache$st
  f <- ctx$f
  p <- ctx$p
  tab <- .STAB_TABLES$SROCKC2
  nf <- .srock_eig(ctx, st, t, u, a)
  mdeg <- floor(sqrt((2 * abs(dt) * abs(st$eig) + 1.5) / 0.811) + 1)
  d <- .srock_deg(tab$ms, max(3, min(mdeg, 200)) - 2)
  mdeg <- d$mdeg
  di <- d$deg_index
  sigma <- tab$msigma[di]
  tau <- tab$mtau[di]
  sqrt_dt <- sqrt(dt)
  s <- .srock2_stages(f, p, t, u, dt, 1, tab$recf, d$start, mdeg)
  uim1 <- s$uim1
  tim1 <- s$tim1
  # the two stage finishing procedure
  ui <- f(uim1, p, tim1)
  uim2 <- uim1 + dt * sigma * ui
  unew <- uim1 + dt * (sigma - tau) * ui
  tim2 <- tim1 + dt * sigma
  uim2 <- f(uim2, p, tim2)
  unew <- unew + dt * (sigma + tau) * uim2
  nf <- nf + s$nf + 2L
  if (type != "general") {
    Gs <- .srock_g(g, uim1, p, tim1, type)
    unew <- unew + Gs * dW
    w <- (dW^2 - dt) / (2 * sqrt_dt)
    unew <- unew + 0.5 * sqrt_dt * .srock_g(g, uim1 + w * Gs, p, tim1, type)
    unew <- unew - 0.5 * sqrt_dt * .srock_g(g, uim1 - w * Gs, p, tim1, type)
  } else {
    Gs <- .srock_g(g, uim1, p, tim1, type)
    uim2 <- .srock_mv(Gs, dW)
    unew <- unew + uim2
    for (i in seq_along(dW)) {
      ui <- Gs[, i]
      tmp <- uim1 - 0.5 * sqrt_dt * ui + 0.5 * (dW[i] / sqrt_dt) * uim2
      unew <- unew + 0.5 * sqrt_dt * .srock_g(g, tmp, p, tim1, type)[, i]
      tmp <- uim1 + 0.5 * sqrt_dt * ui - 0.5 * (dW[i] / sqrt_dt) * uim2
      unew <- unew - 0.5 * sqrt_dt * .srock_g(g, tmp, p, tim1, type)[, i]
    }
  }
  .srock_result(unew, cache, nf)
}

#' @rdname sde_stabilized_methods
#' @export
SROCK1 <- function(interpretation = c("Ito", "Stratonovich"), eigen_est = NULL) {
  interpretation <- match.arg(interpretation)
  .srock_alg("SROCK1", 0.5, .srock1_step, interpretation = interpretation, eigen_est = eigen_est)
}

#' @rdname sde_stabilized_methods
#' @export
SROCK2 <- function(eigen_est = NULL) .srock_alg("SROCK2", 1, .srock2_step, eigen_est = eigen_est)

#' @rdname sde_stabilized_methods
#' @export
SROCKEM <- function(strong_order_1 = TRUE, eigen_est = NULL)
  .srock_alg("SROCKEM", if (strong_order_1) 1 else 0.5, .srockem_step, strong_order_1 = strong_order_1,
             eigen_est = eigen_est)

#' @rdname sde_stabilized_methods
#' @export
SKSROCK <- function(post_processing = FALSE, eigen_est = NULL)
  .srock_alg("SKSROCK", 0.5, .sksrock_step, post_processing = post_processing, eigen_est = eigen_est)

#' @rdname sde_stabilized_methods
#' @export
TangXiaoSROCK2 <- function(version_num = 5, eigen_est = NULL) {
  if (!(version_num %in% 1:5)) stop("version_num must be 1, 2, 3, 4 or 5", call. = FALSE)
  .srock_alg("TangXiaoSROCK2", 1, .tangxiao_step, version_num = as.integer(version_num), eigen_est = eigen_est)
}

#' @rdname sde_stabilized_methods
#' @export
KomBurSROCK2 <- function(eigen_est = NULL)
  .srock_alg("KomBurSROCK2", 1, .komburk_step, interpretation = "Stratonovich", eigen_est = eigen_est)

#' @rdname sde_stabilized_methods
#' @export
SROCKC2 <- function(eigen_est = NULL) .srock_alg("SROCKC2", 1, .srockc2_step, eigen_est = eigen_est)
