# Noise processes, after DiffEqNoiseProcess.jl: types.jl (NoiseProcess,
# SimpleNoiseProcess, NoiseWrapper, NoiseFunction, NoiseGrid,
# NoiseApproximation, VirtualBrownianTree, BoxWedgeTail, NoiseTransport),
# wiener.jl, geometric_bm.jl, ornstein_uhlenbeck.jl, compoundpoisson.jl,
# bridges.jl, correlated_noisefunc.jl, pCN.jl, rswm.jl, solve.jl
# (NoiseProblem) and noise_interfaces/ (accept_step!, reject_step!,
# calculate_step!, setup_next_step!, save_noise!, interpolation).
#
# A noise process is a function of time (`W(t)` gives its value, bridging
# between saved points and stepping past the last one) whose state lives in
# an environment: `W$t`, `W$W` (the saved path), `W$curt`, `W$curW`, `W$dt`,
# `W$dW` (the pending step). The increment distribution is
# `dist(dW, W, dt, u, p, t, rng)` and the bridge `bridge(dW, W, W0, Wh, q, h,
# u, p, t, rng)`, as in Julia: `dist` returns an increment over `dt` from the
# current value, `bridge` the increment over the fraction `q` of a step of
# length `h` that starts at `W0` and has the increment `Wh`. R draws from its
# own generator, so `rng` is `NULL`.

.np_env <- function(x) environment(x)$e

#' @export
`$.AbstractNoiseProcess` <- function(x, name) get0(name, envir = .np_env(x), inherits = FALSE)

#' @export
`$<-.AbstractNoiseProcess` <- function(x, name, value) {
  assign(name, value, envir = .np_env(x))
  x
}

#' @export
`[[.AbstractNoiseProcess` <- function(x, i, ...) get0(i, envir = .np_env(x), inherits = FALSE)

#' @export
print.AbstractNoiseProcess <- function(x, ...) {
  e <- .np_env(x)
  cat(class(x)[1L], " at t = ", format(e$curt), " with ", length(e$t), " saved points\n", sep = "")
  invisible(x)
}

# A noise object: the function `W(t, u, p)` gives its value, its closure
# holds the state environment `e`; `e$kind` is "process" (random increments
# with rejection memory) or "path" (a deterministic function of t once
# constructed, such as a grid or a function).
.np_object <- function(e, cls) {
  f <- function(t, u = NULL, p = NULL) .np_value(e, t, u, p)
  obj <- structure(f, class = unique(c(cls, "AbstractNoiseProcess", "function")))
  e$self <- obj
  e$cls <- class(obj)
  obj
}

# `e` holds the state; `.np_value` dispatches on `e$kind`.
.np_value <- function(e, t, u = NULL, p = NULL) {
  switch(e$kind,
    process = .np_interpolate(e, t, u, p),
    path = e$path(t, u, p),
    stop("unknown noise kind", call. = FALSE))
}

.np_call_dist <- function(e, dW, dt, u, p, t) .as_num(e$dist(dW, e$self, dt, u, p, t, NULL))
.np_call_bridge <- function(e, dW, W0, Wh, q, h, u, p, t)
  .as_num(e$bridge(dW, e$self, W0, Wh, q, h, u, p, t, NULL))

#' Noise processes
#'
#' The noise processes of DiffEqNoiseProcess.jl, for the `noise` argument of
#' [SDEProblem()] and `RODEProblem()`, or on their own with
#' `NoiseProblem()`. A process is a function: `W(t)` returns its value at
#' `t`, bridged between saved points (the bridged point is then saved) or
#' simulated past the last one; `W$t` and `W$W` hold the saved path.
#'
#' `WienerProcess(t0, W0, Z0)` is Brownian motion started at `W0` (a vector
#' gives independent components); `Z0` adds a second, independent process
#' (`W$Z`), used by solvers that need two Wiener processes. `RealWienerProcess()`,
#' `SimpleWienerProcess()` (no rejection memory: adaptive rejections are an
#' error) and `CorrelatedWienerProcess(Gamma, t0, W0)` (covariance `Gamma` per
#' unit time) are variants. `GeometricBrownianMotionProcess(mu, sigma, t0,
#' W0)` and `OrnsteinUhlenbeckProcess(Theta, mu, sigma, t0, W0)` are exact
#' samples of those processes; `CompoundPoissonProcess(rate, t0, W0)` counts
#' Poisson events with intensity `rate(u, p, t)` (with `computerates = FALSE`
#' the rate is fixed at `rate`). `BrownianBridge(t0, tend, W0, Wend)`,
#' `GeometricBrownianBridge()`, `OrnsteinUhlenbeckBridge()` and
#' `CompoundPoissonBridge()` are the same processes conditioned to reach
#' `Wend` at `tend`.
#'
#' `NoiseProcess(t0, W0, Z0, dist, bridge)` builds a process from its
#' increment distribution `dist(dW, W, dt, u, p, t, rng)` and its bridge
#' `bridge(dW, W, W0, Wh, q, h, u, p, t, rng)` (both return increments; `rng`
#' is `NULL` as R uses its own generator). `rswm = RSWM()` chooses the
#' rejection sampling with memory (`adaptivealg = "RSwM3"` keeps the rejected
#' future increments, `"RSwM0"` discards them, for state dependent noise).
#' `SimpleNoiseProcess()` has no rejection memory.
#'
#' `NoiseGrid(t, W, Z)` is a fixed path interpolated linearly;
#' `NoiseFunction(t0, W, Z)` a path given by `W(u, p, t)`;
#' `NoiseWrapper(source)` replays a saved process (a noise process, the
#' Brownian path `sol$W` of an SDE solution, or the solution itself), bridged
#' between its points, optionally reversed in time or restricted to the
#' components `indx`; `NoiseApproximation(source1, source2)` takes the path
#' from the solution of a differential equation integrator (see [init()]);
#' `NoiseTransport(t0, W, RV, rv)` is `W(u, p, t, rv)` with a random variable
#' `rv` drawn from `RV()` at every reset; `VirtualBrownianTree(t0, W0, tend =
#' , Wend = )` is Brownian motion evaluated by bisection from stored seeds, so
#' any value can be recomputed without storing the path; `BoxWedgeTail(t0,
#' W0)` is two dimensional Brownian motion whose steps also sample the Levy
#' area `W$A` from its exact joint density with the increment (Gaines and
#' Lyons). `pCN(source, rho)` is the preconditioned Crank-Nicolson proposal
#' `rho W + sqrt(1 - rho^2) W'` with an independent Brownian path `W'`.
#'
#' `accept_step()`, `reject_step()`, `calculate_step()`, `setup_next_step()`
#' and `save_noise()` step a process by hand, as an SDE integrator does:
#' `calculate_step(W, dt)` draws the pending increment `W$dW` over `dt`,
#' `accept_step(W, dt)` adds it and draws the next one, `reject_step(W,
#' dtnew)` replaces it by a bridged increment over `dtnew`, keeping the rest
#' for later steps. `resize_noise()`, `addat_noise()` and `deleteat_noise()`
#' change the number of components (as callbacks that resize the state do).
#'
#' @param t0,tend Initial and final times.
#' @param W0,Z0,Wend,Zend Initial and final values.
#' @param mu,sigma,Theta Process parameters.
#' @param Gamma Covariance matrix per unit time.
#' @param rate Event intensity: a function `rate(u, p, t)` or a number.
#' @param computerates Recompute the rate at every step.
#' @param dist,bridge Increment distribution and bridge, see Details.
#' @param rswm Rejection sampling with memory, from `RSWM()`.
#' @param discard_length Shortest stored future increment.
#' @param adaptivealg `"RSwM3"`, `"RSwM2"`, `"RSwM1"` or `"RSwM0"`.
#' @param save_everystep Save every step of the path.
#' @param reset,reseed Restart from `t0` (and draw a new path) when a solve
#'   uses the process.
#' @param continuous Whether the path is continuous (for interpolation).
#' @param covariance Covariance of a correlated process.
#' @param t,W,Z,RV,rv Path data, see Details.
#' @param source,source1,source2 Noise sources, see Details.
#' @param reverse,indx Reverse time; components of the source to use.
#' @param rho Correlation of the pCN proposal.
#' @param tree_depth,search_depth,atol Virtual Brownian tree resolution.
#' @param rtol,nr,na,nz,box_grouping,sqeezing Box wedge tail sampler
#'   settings; R samples the area by inverting its conditional distribution
#'   to `rtol`.
#' @param noise A noise process.
#' @param tspan Time span.
#' @param seed Seed of the random number generator.
#' @param dt,dtnew Step sizes.
#' @param u,p State and parameters for state dependent noise.
#' @param setup_next Draw the next increment after accepting.
#' @param i,n Components to add or delete; new number of components.
#' @param ... Options.
#' @return A noise process, a problem or the process stepped in place.
#' @examples
#' W <- WienerProcess(0, 0)
#' prob <- NoiseProblem(W, c(0, 1))
#' sol <- solve(prob, dt = 0.1, seed = 1)
#' sol(0.55)
#' ou <- OrnsteinUhlenbeckProcess(1, 0, 0.3, 0, 1)
#' f <- function(u, p, t) -u
#' g <- function(u, p, t) 0.1
#' s <- solve(SDEProblem(f, g, 1, c(0, 1), noise = ou), EM(), dt = 0.01, seed = 2)
#' s$retcode
#' @name noise_processes
NULL

# ---------------------------------------------------------------------------
# Rejection sampling with memory

#' @rdname noise_processes
#' @export
RSWM <- function(discard_length = 1e-15, adaptivealg = c("RSwM3", "RSwM2", "RSwM1", "RSwM0"))
  structure(list(discard_length = discard_length, adaptivealg = match.arg(adaptivealg)), class = "RSWM")


# ---------------------------------------------------------------------------
# The general process

.np_process <- function(t0, W0, Z0, dist, bridge, rswm = RSWM(), save_everystep = TRUE, covariance = NULL,
                        reset = TRUE, reseed = TRUE, continuous = TRUE, simple = FALSE, cls = "NoiseProcess",
                        extra = list()) {
  e <- new.env(parent = emptyenv())
  e$kind <- "process"
  e$t0 <- t0; e$W0 <- .as_num(W0); e$Z0 <- if (!is.null(Z0)) .as_num(Z0)
  e$dist <- dist; e$bridge <- bridge
  e$rswm <- rswm; e$save_everystep <- save_everystep; e$covariance <- covariance
  e$reset <- reset; e$reseed <- reseed; e$continuous <- continuous; e$simple <- simple
  for (nm in names(extra)) assign(nm, extra[[nm]], envir = e)
  .np_reinit(e, t0)
  .np_object(e, c(cls, "NoiseProcess"))
}

# Back to t0 with an empty path and empty stacks (reinit! of the process).
.np_reinit <- function(e, t0 = e$t0, dt = 0) {
  e$curt <- t0
  e$curW <- e$W0
  e$curZ <- e$Z0
  e$t <- t0
  e$W <- list(e$W0)
  e$Z <- if (!is.null(e$Z0)) list(e$Z0)
  e$dt <- dt
  e$dW <- 0 * e$W0
  e$dZ <- if (!is.null(e$Z0)) 0 * e$Z0
  e$S1 <- list(); e$S2 <- list()
  e$iter <- 0L
  if (!is.null(e$on_reinit)) e$on_reinit(e)
  invisible(e)
}

#' @rdname noise_processes
#' @export
NoiseProcess <- function(t0, W0, Z0 = NULL, dist, bridge = NULL, rswm = RSWM(), save_everystep = TRUE,
                         covariance = NULL, reset = TRUE, reseed = TRUE, continuous = TRUE) {
  .np_process(t0, W0, Z0, dist, bridge, rswm, save_everystep, covariance, reset, reseed, continuous)
}

#' @rdname noise_processes
#' @export
SimpleNoiseProcess <- function(t0, W0, Z0 = NULL, dist, bridge = NULL, save_everystep = TRUE,
                               covariance = NULL, reset = TRUE, reseed = TRUE) {
  .np_process(t0, W0, Z0, dist, bridge, RSWM(), save_everystep, covariance, reset, reseed, TRUE,
              simple = TRUE, cls = "SimpleNoiseProcess")
}

# WHITE_NOISE_DIST and WHITE_NOISE_BRIDGE (wiener.jl)
.white_dist <- function(dW, W, dt, u, p, t, rng) sqrt(abs(dt)) * stats::rnorm(length(dW))
.white_bridge <- function(dW, W, W0, Wh, q, h, u, p, t, rng)
  sqrt((1 - q) * q * abs(h)) * stats::rnorm(length(dW)) + q * Wh

#' @rdname noise_processes
#' @export
WienerProcess <- function(t0, W0, Z0 = NULL, ...) {
  # symcompR's WienerProcess(sym) is the stochastic process of sympy.stats,
  # named by a string or a symbol
  if ((missing(t0) || is.character(t0) || is_sym(t0)) && .symcompr_has("WienerProcess"))
    return(.symcompr_forward("WienerProcess", sys.call(), parent.frame()))
  .np_process(t0, W0, Z0, .white_dist, .white_bridge, ..., cls = "WienerProcess")
}

#' @rdname noise_processes
#' @export
RealWienerProcess <- function(t0, W0, Z0 = NULL, ...)
  .np_process(t0, W0, Z0, .white_dist, .white_bridge, ..., cls = "RealWienerProcess")

#' @rdname noise_processes
#' @export
SimpleWienerProcess <- function(t0, W0, Z0 = NULL, ...)
  .np_process(t0, W0, Z0, .white_dist, .white_bridge, ..., simple = TRUE, cls = "SimpleWienerProcess")

#' @rdname noise_processes
#' @export
CorrelatedWienerProcess <- function(Gamma, t0, W0, Z0 = NULL, ...) {
  s <- svd(as.matrix(Gamma))
  A <- s$u %*% diag(sqrt(s$d), length(s$d))
  dist <- function(dW, W, dt, u, p, t, rng) .as_num(A %*% (sqrt(abs(dt)) * stats::rnorm(ncol(A))))
  # the bridge of A B for a standard Brownian motion B (DiffEqNoiseProcess
  # has none and forbids adaptive rejections)
  bridge <- function(dW, W, W0, Wh, q, h, u, p, t, rng)
    .as_num(A %*% (sqrt((1 - q) * q * abs(h)) * stats::rnorm(ncol(A)))) + q * Wh
  .np_process(t0, W0, Z0, dist, bridge, ..., covariance = Gamma, cls = "CorrelatedWienerProcess")
}

#' @rdname noise_processes
#' @export
GeometricBrownianMotionProcess <- function(mu, sigma, t0, W0, Z0 = NULL, ...) {
  dist <- function(dW, W, dt, u, p, t, rng) {
    cur <- W$curW
    cur * expm1((mu - sigma^2 / 2) * dt + sigma * sqrt(dt) * stats::rnorm(length(dW)))
  }
  bridge <- function(dW, W, W0, Wh, q, h, u, p, t, rng) {
    m <- log(abs(W0)) + q * (log(abs(W0 + Wh)) - log(abs(W0)))
    exp(m + sqrt((1 - q) * q * abs(h)) * sigma * stats::rnorm(length(dW))) - W0
  }
  .np_process(t0, W0, Z0, dist, bridge, ..., cls = "GeometricBrownianMotionProcess",
              extra = list(mu = mu, sigma = sigma))
}

#' @rdname noise_processes
#' @export
OrnsteinUhlenbeckProcess <- function(Theta, mu, sigma, t0, W0, Z0 = NULL, ...) {
  dist <- function(dW, W, dt, u, p, t, rng) {
    cur <- W$curW
    mu + (cur - mu) * exp(-Theta * dt) + stats::rnorm(length(dW)) * sigma * sqrt(-expm1(-2 * Theta * dt) / (2 * Theta)) - cur
  }
  bridge <- function(dW, W, W0, Wh, q, h, u, p, t, rng) {
    s <- function(x) sinh(Theta * x)
    var <- sigma^2 * s(h * (1 - q)) * s(q * h) / (Theta * s(h))
    (W0 - mu) * (s(h * (1 - q)) / s(h) - 1) + (Wh + W0 - mu) * s(q * h) / s(h) + sqrt(var) * stats::rnorm(length(dW))
  }
  .np_process(t0, W0, Z0, dist, bridge, ..., cls = "OrnsteinUhlenbeckProcess",
              extra = list(Theta = Theta, mu = mu, sigma = sigma))
}

#' @rdname noise_processes
#' @export
CompoundPoissonProcess <- function(rate, t0, W0, computerates = TRUE, rswm = RSWM(adaptivealg = "RSwM0"), ...) {
  cur <- new.env(parent = emptyenv())
  cur$rate <- if (is.function(rate)) NULL else rate
  dist <- function(dW, W, dt, u, p, t, rng) {
    if (computerates && is.function(rate)) cur$rate <- rate(u, p, t)
    r <- rep_len(.as_num(cur$rate), length(dW))
    .as_num(stats::rpois(length(dW), dt * r))
  }
  bridge <- function(dW, W, W0, Wh, q, h, u, p, t, rng)
    .as_num(stats::rbinom(length(Wh), as.integer(round(Wh)), q))
  .np_process(t0, W0, NULL, dist, bridge, rswm = rswm, ..., continuous = FALSE, cls = "CompoundPoissonProcess",
              extra = list(rate = rate, computerates = computerates))
}

# A bridge is the process with the increment to the end already on the stack.
.np_with_end <- function(W, tend, Wend, Zend = NULL) {
  e <- .np_env(W)
  e$bridge_end <- list(h = tend - e$t0, dW = .as_num(Wend) - e$W0,
                       dZ = if (!is.null(Zend)) .as_num(Zend) - e$Z0)
  e$on_reinit <- function(e) e$S1 <- list(e$bridge_end)
  e$S1 <- list(e$bridge_end)
  W
}

#' @rdname noise_processes
#' @export
BrownianBridge <- function(t0, tend, W0, Wend, Z0 = NULL, Zend = NULL, ...)
  .np_with_end(.np_process(t0, W0, Z0, .white_dist, .white_bridge, ..., cls = "BrownianBridge"), tend, Wend, Zend)

#' @rdname noise_processes
#' @export
GeometricBrownianBridge <- function(mu, sigma, t0, tend, W0, Wend, Z0 = NULL, Zend = NULL, ...) {
  W <- GeometricBrownianMotionProcess(mu, sigma, t0, W0, Z0, ...)
  class(W)[1L] <- "GeometricBrownianBridge"
  .np_with_end(W, tend, Wend, Zend)
}

#' @rdname noise_processes
#' @export
OrnsteinUhlenbeckBridge <- function(Theta, mu, sigma, t0, tend, W0, Wend, Z0 = NULL, ...) {
  W <- OrnsteinUhlenbeckProcess(Theta, mu, sigma, t0, W0, Z0, ...)
  class(W)[1L] <- "OrnsteinUhlenbeckBridge"
  .np_with_end(W, tend, Wend)
}

#' @rdname noise_processes
#' @export
CompoundPoissonBridge <- function(rate, t0, tend, W0, Wend, rswm = RSWM(), ...) {
  W <- CompoundPoissonProcess(rate, t0, W0, rswm = rswm, ...)
  class(W)[1L] <- "CompoundPoissonBridge"
  .np_with_end(W, tend, Wend)
}

# ---------------------------------------------------------------------------
# Stepping a process (noise_process_interface.jl)

# The pending increment over dt from the current value: future increments on
# the stack S1 first (RSwM), then fresh draws.
.np_draw <- function(e, dt, u, p) {
  need <- dt
  dW <- 0 * e$curW
  dZ <- if (!is.null(e$curZ)) 0 * e$curZ
  tol <- 1e-12 * max(abs(dt), 0x1.56e1fc2f8f359p-997)
  while (abs(need) > tol && length(e$S1)) {
    s <- e$S1[[1L]]
    if (abs(s$h) <= abs(need) + tol) {
      dW <- dW + s$dW
      if (!is.null(dZ)) dZ <- dZ + (s$dZ %||% 0)
      need <- need - s$h
      e$S1 <- e$S1[-1L]
    } else {
      q <- need / s$h
      W0 <- e$curW + dW
      w1 <- .np_call_bridge(e, s$dW, W0, s$dW, q, s$h, u, p, e$curt)
      z1 <- if (!is.null(dZ)) .np_call_bridge(e, s$dZ, e$curZ + dZ, s$dZ, q, s$h, u, p, e$curt)
      e$S1[[1L]] <- list(h = s$h - need, dW = s$dW - w1, dZ = if (!is.null(dZ)) s$dZ - z1)
      dW <- dW + w1
      if (!is.null(dZ)) dZ <- dZ + z1
      need <- 0
    }
  }
  if (abs(need) > tol) {
    old <- e$curW
    e$curW <- old + dW
    dW <- dW + .np_call_dist(e, e$curW, need, u, p, e$curt)
    e$curW <- old
    if (!is.null(dZ)) {
      oz <- e$curZ
      e$curZ <- oz + dZ
      dZ <- dZ + .np_call_dist(e, e$curZ, need, u, p, e$curt)
      e$curZ <- oz
    }
  }
  list(dW = dW, dZ = dZ)
}

#' @rdname noise_processes
#' @export
calculate_step <- function(W, dt, u = NULL, p = NULL) {
  e <- .np_env(W)
  if (e$kind == "path") {
    e$dW <- .as_num(e$path(e$curt + dt, u, p)) - e$curW
    e$dt <- dt
    return(invisible(W))
  }
  d <- .np_draw(e, dt, u, p)
  e$dW <- d$dW; e$dZ <- d$dZ; e$dt <- dt
  invisible(W)
}

#' @rdname noise_processes
#' @export
setup_next_step <- function(W, u = NULL, p = NULL) calculate_step(W, .np_env(W)$dt, u, p)

#' @rdname noise_processes
#' @export
accept_step <- function(W, dt, u = NULL, p = NULL, setup_next = TRUE) {
  e <- .np_env(W)
  e$curt <- e$curt + e$dt
  e$iter <- e$iter + 1L
  e$curW <- e$curW + e$dW
  if (!is.null(e$curZ)) e$curZ <- e$curZ + e$dZ
  if (!is.null(e$on_accept)) e$on_accept(e)
  if (e$save_everystep) save_noise(W)
  e$dt <- dt
  if (setup_next) setup_next_step(W, u, p)
  invisible(W)
}

#' @rdname noise_processes
#' @export
reject_step <- function(W, dtnew, u = NULL, p = NULL) {
  e <- .np_env(W)
  if (isTRUE(e$simple)) stop(class(W)[1L], " cannot be used with adaptivity rejections", call. = FALSE)
  if (isTRUE(e$no_reject)) stop(class(W)[1L], " cannot be used with adaptivity rejections", call. = FALSE)
  if (e$kind == "path") return(calculate_step(W, dtnew, u, p))
  q <- dtnew / e$dt
  w1 <- .np_call_bridge(e, e$dW, e$curW, e$dW, q, e$dt, u, p, e$curt)
  z1 <- if (!is.null(e$dZ)) .np_call_bridge(e, e$dZ, e$curZ, e$dZ, q, e$dt, u, p, e$curt)
  cut <- e$dt - dtnew
  if (e$rswm$adaptivealg != "RSwM0" && abs(cut) > e$rswm$discard_length)
    e$S1 <- c(list(list(h = cut, dW = e$dW - w1, dZ = if (!is.null(z1)) e$dZ - z1)), e$S1)
  e$dW <- w1; e$dZ <- z1; e$dt <- dtnew
  invisible(W)
}

#' @rdname noise_processes
#' @export
save_noise <- function(W) {
  e <- .np_env(W)
  if (abs(e$t[length(e$t)] - e$curt) > 0 || length(e$t) == 0L) {
    e$t <- c(e$t, e$curt)
    e$W[[length(e$W) + 1L]] <- e$curW
    if (!is.null(e$Z)) e$Z[[length(e$Z) + 1L]] <- e$curZ
  }
  invisible(W)
}

# Value at t: past the last saved point the process steps there, inside the
# saved path it is bridged between the neighbors and the point is saved.
.np_interpolate <- function(e, t, u = NULL, p = NULL) {
  ts <- e$t
  n <- length(ts)
  tol <- 100 * .Machine$double.eps * max(1, abs(t))
  hit <- which(abs(ts - t) <= tol)
  if (length(hit)) return(e$W[[hit[1L]]])
  if (t > ts[n]) {
    old <- e$curW
    e$curW <- e$W[[n]]
    val <- e$curW + .np_call_dist(e, e$curW, t - ts[n], u, p, ts[n])
    e$curW <- old
    e$t <- c(ts, t)
    e$W[[n + 1L]] <- val
    if (!is.null(e$Z)) {
      oz <- e$curZ
      e$curZ <- e$Z[[n]]
      e$Z[[n + 1L]] <- e$curZ + .np_call_dist(e, e$curZ, t - ts[n], u, p, ts[n])
      e$curZ <- oz
    }
    return(val)
  }
  if (t < ts[1L]) stop("cannot interpolate the noise before its first time point", call. = FALSE)
  i <- findInterval(t, ts)
  W0 <- e$W[[i]]; Wh <- e$W[[i + 1L]] - W0
  h <- ts[i + 1L] - ts[i]; q <- (t - ts[i]) / h
  if (is.null(e$bridge)) stop("this noise process has no bridge", call. = FALSE)
  val <- W0 + .np_call_bridge(e, Wh, W0, Wh, q, h, u, p, ts[i])
  e$t <- append(ts, t, after = i)
  e$W <- append(e$W, list(val), after = i)
  if (!is.null(e$Z)) {
    Z0 <- e$Z[[i]]; Zh <- e$Z[[i + 1L]] - Z0
    e$Z <- append(e$Z, list(Z0 + .np_call_bridge(e, Zh, Z0, Zh, q, h, u, p, ts[i])), after = i)
  }
  val
}

#' @rdname noise_processes
#' @export
resize_noise <- function(W, n) {
  e <- .np_env(W)
  fit <- function(v, fill) if (length(v) >= n) v[seq_len(n)] else c(v, rep(fill, n - length(v)))
  e$curW <- fit(e$curW, 0); e$dW <- fit(e$dW, 0); e$W0 <- fit(e$W0, 0)
  if (!is.null(e$curZ)) { e$curZ <- fit(e$curZ, 0); e$dZ <- fit(e$dZ, 0); e$Z0 <- fit(e$Z0, 0) }
  e$W <- lapply(e$W, fit, fill = 0)
  if (!is.null(e$Z)) e$Z <- lapply(e$Z, fit, fill = 0)
  e$S1 <- lapply(e$S1, function(s) list(h = s$h, dW = fit(s$dW, NA), dZ = if (!is.null(s$dZ)) fit(s$dZ, NA)))
  # new components of stored future increments get fresh draws
  e$S1 <- lapply(e$S1, function(s) {
    na <- is.na(s$dW)
    if (any(na)) s$dW[na] <- sqrt(abs(s$h)) * stats::rnorm(sum(na))
    if (!is.null(s$dZ) && any(is.na(s$dZ))) s$dZ[is.na(s$dZ)] <- sqrt(abs(s$h)) * stats::rnorm(sum(is.na(s$dZ)))
    s
  })
  invisible(W)
}

#' @rdname noise_processes
#' @export
addat_noise <- function(W, i) {
  e <- .np_env(W)
  add <- function(v) append(v, rep(0, length(i)), after = min(i) - 1L)
  e$curW <- add(e$curW); e$dW <- add(e$dW); e$W0 <- add(e$W0)
  if (!is.null(e$curZ)) { e$curZ <- add(e$curZ); e$dZ <- add(e$dZ); e$Z0 <- add(e$Z0) }
  e$W <- lapply(e$W, add)
  if (!is.null(e$Z)) e$Z <- lapply(e$Z, add)
  e$S1 <- lapply(e$S1, function(s) {
    s$dW <- append(s$dW, sqrt(abs(s$h)) * stats::rnorm(length(i)), after = min(i) - 1L)
    if (!is.null(s$dZ)) s$dZ <- append(s$dZ, sqrt(abs(s$h)) * stats::rnorm(length(i)), after = min(i) - 1L)
    s
  })
  invisible(W)
}

#' @rdname noise_processes
#' @export
deleteat_noise <- function(W, i) {
  e <- .np_env(W)
  del <- function(v) v[-i]
  e$curW <- del(e$curW); e$dW <- del(e$dW); e$W0 <- del(e$W0)
  if (!is.null(e$curZ)) { e$curZ <- del(e$curZ); e$dZ <- del(e$dZ); e$Z0 <- del(e$Z0) }
  e$W <- lapply(e$W, del)
  if (!is.null(e$Z)) e$Z <- lapply(e$Z, del)
  e$S1 <- lapply(e$S1, function(s) { s$dW <- del(s$dW); if (!is.null(s$dZ)) s$dZ <- del(s$dZ); s })
  invisible(W)
}

#' @rdname noise_processes
#' @export
fill_new_noise_caches <- function(W, i) invisible(W)

# ---------------------------------------------------------------------------
# Path noises: grids, functions, wrappers, approximations, transports, trees

.np_path <- function(t0, W0, path, Z0 = NULL, reset = TRUE, cls, extra = list(), no_reject = FALSE) {
  e <- new.env(parent = emptyenv())
  e$kind <- "path"
  e$t0 <- t0; e$W0 <- .as_num(W0); e$Z0 <- if (!is.null(Z0)) .as_num(Z0)
  e$path <- path
  e$reset <- reset; e$reseed <- FALSE; e$save_everystep <- TRUE
  e$rswm <- RSWM(); e$no_reject <- no_reject
  for (nm in names(extra)) assign(nm, extra[[nm]], envir = e)
  .np_reinit(e, t0)
  .np_object(e, cls)
}

#' @rdname noise_processes
#' @export
NoiseGrid <- function(t, W, Z = NULL, reset = TRUE) {
  Wm <- if (is.list(W)) do.call(rbind, lapply(W, as.numeric)) else as.matrix(W)
  if (!is.list(W) && is.null(dim(W))) Wm <- matrix(W, ncol = 1L)
  Zm <- if (!is.null(Z)) { if (is.list(Z)) do.call(rbind, lapply(Z, as.numeric)) else as.matrix(Z) }
  if (nrow(Wm) != length(t)) stop("NoiseGrid needs one value of W per time point", call. = FALSE)
  # a grid for a backward solve runs from the last time to the first; the
  # interpolation below needs increasing times
  rev_grid <- length(t) > 1L && t[1L] > t[length(t)]
  ta <- if (rev_grid) rev(t) else t
  ascend <- function(M) if (rev_grid && !is.null(M)) M[rev(seq_len(nrow(M))), , drop = FALSE] else M
  Wa <- ascend(Wm)
  Za <- ascend(Zm)
  interp <- function(M, s) {
    if (s > max(t) + 10 * .Machine$double.eps * max(1, abs(s)))
      stop("Solution interpolation cannot extrapolate past the final timepoint. Build a longer NoiseGrid to cover the integration.",
           call. = FALSE)
    if (s < min(t) - 10 * .Machine$double.eps * max(1, abs(s)))
      stop("Solution interpolation cannot extrapolate before the first timepoint. Build a longer NoiseGrid to cover the integration.",
           call. = FALSE)
    apply(M, 2L, function(col) stats::approx(ta, col, xout = s, rule = 2, ties = "ordered")$y)
  }
  path <- function(s, u = NULL, p = NULL) interp(Wa, s)
  .np_path(t[1L], Wm[1L, ], path, if (!is.null(Zm)) Zm[1L, ], reset, "NoiseGrid",
           extra = list(grid_t = t, grid_W = Wm, grid_Z = Zm, zpath = if (!is.null(Zm)) function(s, u, p) interp(Za, s)))
}

#' @rdname noise_processes
#' @export
NoiseFunction <- function(t0, W, Z = NULL, reset = TRUE) {
  path <- function(s, u = NULL, p = NULL) .as_num(W(u, p, s))
  .np_path(t0, W(NULL, NULL, t0), path, if (!is.null(Z)) Z(NULL, NULL, t0), reset, "NoiseFunction",
           extra = list(zpath = if (!is.null(Z)) function(s, u, p) .as_num(Z(u, p, s))))
}

#' @rdname noise_processes
#' @export
NoiseTransport <- function(t0, W, RV, rv = NULL, Z = NULL, reset = TRUE, reseed = TRUE) {
  holder <- new.env(parent = emptyenv())
  holder$rv <- rv %||% RV()
  path <- function(s, u = NULL, p = NULL) .as_num(W(u, p, s, holder$rv))
  obj <- .np_path(t0, W(NULL, NULL, t0, holder$rv), path, if (!is.null(Z)) Z(NULL, NULL, t0, holder$rv),
                  reset, "NoiseTransport",
                  extra = list(RV = RV, holder = holder,
                               zpath = if (!is.null(Z)) function(s, u, p) .as_num(Z(u, p, s, holder$rv))))
  e <- .np_env(obj)
  e$reseed <- reseed
  e$on_reinit <- function(e) if (isTRUE(e$reseed) && !is.null(e$RV)) e$holder$rv <- e$RV()
  obj
}

# The saved path of a source: a noise process, sol$W, or an SDE solution.
.np_source_path <- function(source) {
  if (inherits(source, "AbstractNoiseProcess")) {
    e <- .np_env(source)
    return(list(t = e$t, W = do.call(rbind, e$W), Z = if (!is.null(e$Z)) do.call(rbind, e$Z), proc = source))
  }
  if (is.function(source) && !is.null(source$W)) source <- source$W
  if (is.list(source) && !is.null(source$t) && !is.null(source$W))
    return(list(t = source$t, W = as.matrix(source$W), Z = source$Z, proc = NULL))
  stop("NoiseWrapper needs a noise process or the Brownian path of a solution", call. = FALSE)
}

#' @rdname noise_processes
#' @export
NoiseWrapper <- function(source, reset = TRUE, reverse = FALSE, indx = NULL) {
  sp <- .np_source_path(source)
  ts <- sp$t; Wm <- sp$W
  if (!is.null(indx)) Wm <- Wm[, indx, drop = FALSE]
  if (reverse) { ts <- rev(ts); Wm <- Wm[rev(seq_len(nrow(Wm))), , drop = FALSE] }
  store <- new.env(parent = emptyenv())
  store$t <- ts; store$W <- Wm
  bridge <- if (!is.null(sp$proc)) .np_env(sp$proc)$bridge else .white_bridge
  path <- function(s, u = NULL, p = NULL) {
    tt <- store$t
    tol <- 100 * .Machine$double.eps * max(1, abs(s))
    hit <- which(abs(tt - s) <= tol)
    if (length(hit)) return(store$W[hit[1L], ])
    up <- tt[1L] <= tt[length(tt)]
    o <- if (up) order(tt) else order(-tt)
    tt <- tt[o]; M <- store$W[o, , drop = FALSE]
    inside <- if (up) s > tt[1L] && s < tt[length(tt)] else s < tt[1L] && s > tt[length(tt)]
    if (!inside) stop("NoiseWrapper cannot extrapolate past the saved noise", call. = FALSE)
    i <- if (up) findInterval(s, tt) else findInterval(-s, -tt)
    W0 <- M[i, ]; Wh <- M[i + 1L, ] - W0; h <- tt[i + 1L] - tt[i]; q <- (s - tt[i]) / h
    val <- W0 + .as_num(bridge(Wh, NULL, W0, Wh, q, h, u, p, tt[i], NULL))
    store$t <- c(store$t, s); store$W <- rbind(store$W, val)
    val
  }
  .np_path(ts[1L], Wm[1L, ], path, NULL, reset, "NoiseWrapper",
           extra = list(source = source, store = store, reverse = reverse))
}

#' @rdname noise_processes
#' @export
NoiseApproximation <- function(source1, source2 = NULL, reset = TRUE) {
  value_of <- function(integ, s) {
    while (integ$t < s && !isTRUE(.integ_I(integ)$done)) step_integrator(integ)
    .as_num(integ$interp(s))
  }
  path <- function(s, u = NULL, p = NULL) value_of(source1, s)
  .np_path(source1$t, .as_num(source1$u), path,
           if (!is.null(source2)) .as_num(source2$u), reset, "NoiseApproximation",
           extra = list(source1 = source1, source2 = source2,
                        zpath = if (!is.null(source2)) function(s, u, p) value_of(source2, s)))
}

#' @rdname noise_processes
#' @export
VirtualBrownianTree <- function(t0, W0, Z0 = NULL, dist = .white_dist, bridge = NULL, tend = NULL, Wend = NULL,
                                Zend = NULL, atol = 1e-10, tree_depth = 4L, search_depth = NULL, seed = NULL,
                                reset = TRUE) {
  tend <- tend %||% (t0 + 1)
  W0 <- .as_num(W0)
  seed <- seed %||% sample.int(.Machine$integer.max, 1L)
  if (is.null(Wend)) Wend <- .with_seed(seed, W0 + dist(W0, NULL, tend - t0, NULL, NULL, t0, NULL))
  search_depth <- search_depth %||% max(tree_depth, ceiling(log2(abs(tend - t0) / atol)))
  # each node of the bisection draws its midpoint from its own seed, so a
  # value is recomputed identically from the root
  value <- function(s, a, b, Wa, Wb, node_seed, depth) {
    if (abs(s - a) <= atol) return(Wa)
    if (abs(s - b) <= atol) return(Wb)
    if (depth >= search_depth) return(Wa + (s - a) / (b - a) * (Wb - Wa))
    m <- (a + b) / 2
    seeds <- .with_seed(node_seed, sample.int(.Machine$integer.max, 3L))
    h <- b - a
    Wm <- .with_seed(seeds[1L], Wa + sqrt(h / 4) * stats::rnorm(length(Wa)) + (Wb - Wa) / 2)
    if (s < m) value(s, a, m, Wa, Wm, seeds[2L], depth + 1L) else value(s, m, b, Wm, Wb, seeds[3L], depth + 1L)
  }
  path <- function(s, u = NULL, p = NULL) {
    if (s > tend + atol || s < t0 - atol) stop("VirtualBrownianTree is defined on [t0, tend]", call. = FALSE)
    value(s, t0, tend, W0, .as_num(Wend), seed, 0L)
  }
  .np_path(t0, W0, path, Z0, reset, "VirtualBrownianTree",
           extra = list(tend = tend, Wend = Wend, seed = seed, tree_depth = tree_depth, search_depth = search_depth))
}

# ---------------------------------------------------------------------------
# Box wedge tail: 2d Brownian increments with their Levy area

# P(|a| <= x | r) for the area a of a unit time 2d Brownian motion whose
# increment has radius r: the joint density of Gaines and Lyons is
# (r / pi) int_0^Inf (x / sinh x) exp(-r^2 x / (2 tanh x)) cos(a x) dx.
.bwt_area_cdf <- function(a, r, rtol) {
  phi <- function(x) ifelse(x < 1e-8, 1, (x / sinh(x)) * exp(-r^2 / 2 * (x / tanh(x) - 1)))
  (2 / pi) * stats::integrate(function(x) phi(x) * ifelse(x < 1e-12, a, sin(a * x) / x), 0, Inf,
                              rel.tol = rtol, subdivisions = 1000L)$value
}

.bwt_sample <- function(dt, rtol) {
  r <- sqrt(-2 * log(stats::runif(1)))
  U <- stats::runif(1)
  hi <- 1
  while (.bwt_area_cdf(hi, r, rtol) < U) hi <- 2 * hi
  a <- stats::uniroot(function(x) .bwt_area_cdf(x, r, rtol) - U, c(0, hi), tol = rtol * max(1, hi))$root
  theta <- stats::runif(1, 0, 2 * pi)
  list(dW = sqrt(abs(dt)) * r * c(cos(theta), sin(theta)), A = sample(c(-1, 1), 1L) * a * dt)
}

#' @rdname noise_processes
#' @export
BoxWedgeTail <- function(t0, W0, Z0 = NULL, dist = .white_dist, bridge = .white_bridge, rtol = 1e-8, nr = 4L, na = 4L,
                         nz = 10L, box_grouping = c("MinEntropy", "Columns", "none"), sqeezing = TRUE,
                         save_everystep = TRUE, reset = TRUE, reseed = TRUE) {
  if (length(W0) != 2L) stop("The BoxWedgeTail algorithms can only be used for 2d Brownian processes.", call. = FALSE)
  box_grouping <- match.arg(box_grouping)
  rt <- max(rtol, 1e-10)
  bwt_dist <- function(dW, W, dt, u, p, t, rng) {
    s <- .bwt_sample(dt, rt)
    assign("curA", s$A, envir = .np_env(W))
    s$dW
  }
  W <- .np_process(t0, W0, Z0, bwt_dist, bridge, save_everystep = save_everystep, reset = reset, reseed = reseed,
                   cls = "BoxWedgeTail", extra = list(rtol = rtol, box_grouping = box_grouping, curA = 0, A = 0,
                                                      no_reject = TRUE))
  e <- .np_env(W)
  e$on_reinit <- function(e) { e$curA <- 0; e$A <- 0 }
  e$on_accept <- function(e) e$A <- c(e$A, e$curA)
  W
}

# ---------------------------------------------------------------------------
# pCN

#' @rdname noise_processes
#' @export
pCN <- function(source, rho, reset = TRUE, reverse = FALSE, indx = NULL) {
  sp <- .np_source_path(source)
  ts <- sp$t
  innov <- function(M) {
    inc <- sqrt(diff(ts)) * matrix(stats::rnorm((length(ts) - 1L) * ncol(M)), ncol = ncol(M))
    rbind(0 * M[1L, ], apply(inc, 2L, cumsum))
  }
  Wn <- rho * sp$W + sqrt(1 - rho^2) * innov(sp$W)
  if (inherits(source, "NoiseGrid")) return(NoiseGrid(ts, Wn, sp$Z, reset = reset))
  NoiseWrapper(list(t = ts, W = Wn), reset = reset, reverse = reverse, indx = indx)
}

# ---------------------------------------------------------------------------
# NoiseProblem

#' @rdname noise_processes
#' @export
NoiseProblem <- function(noise, tspan, seed = 0L, ...)
  structure(list(noise = noise, tspan = .as_num(tspan), seed = seed, kwargs = list(...)),
            class = c("NoiseProblem", "SciMLProblem"))

.np_copy <- function(W) {
  e <- .np_env(W)
  e2 <- list2env(as.list(e, all.names = TRUE), parent = emptyenv())
  .np_object(e2, class(W))
}

# The process a solve runs on (the noise setup of StochasticDiffEq's
# __init): a copy of the problem's process, reinitialized at t only when it
# is not there. Upstream the problem's process is aliased, so a solve leaves
# it at the end of the span and the next solve reinitializes it (and a
# NoiseTransport draws a new realization); `solved` on the original stands
# for that, since here every solve runs on a copy. The first solve keeps a
# NoiseTransport's supplied realization.
.np_start <- function(W, t) {
  e0 <- .np_env(W)
  proc <- .np_copy(W)
  e <- .np_env(proc)
  if (isTRUE(e$reset) && (isTRUE(e0$solved) || e$curt != t)) .np_reinit(e, t)
  e0$solved <- TRUE
  proc
}

#' @export
solve.NoiseProblem <- function(a, b = NULL, ...) {
  opts <- utils::modifyList(a$kwargs, list(...))
  dt <- opts$dt
  if (is.null(dt) || dt == 0) stop("dt must be provided to simulate a noise process. Please pass dt=...", call. = FALSE)
  seed <- opts$seed %||% (if (a$seed != 0) a$seed)
  .with_seed(seed, {
    W <- .np_copy(a$noise)
    e <- .np_env(W)
    if (isTRUE(e$reset)) .np_reinit(e, a$tspan[1L], dt)
    tend <- a$tspan[2L]
    if (e$kind == "path") {
      while (e$curt < tend - 100 * .Machine$double.eps * max(1, abs(tend))) {
        h <- min(dt, tend - e$curt)
        calculate_step(W, h)
        accept_step(W, h, setup_next = FALSE)
      }
    } else {
      e$dt <- dt
      setup_next_step(W)
      while (e$curt < tend) {
        endtol <- 100 * .Machine$double.eps * max(abs(tend), abs(e$curt))
        if (abs(tend - (e$curt + e$dt)) <= endtol) {
          accept_step(W, dt)
          e$curt <- tend
          e$t[length(e$t)] <- tend
        } else if (e$curt + e$dt > tend + endtol) {
          calculate_step(W, tend - e$curt)
          accept_step(W, tend - e$curt)
          e$curt <- tend
          e$t[length(e$t)] <- tend
        } else accept_step(W, dt)
      }
    }
    W
  })
}

# ---------------------------------------------------------------------------
# Two and three point random variables of the weak methods
# (StochasticDiffEqCore/src/weak_utils.jl)

#' Random variables of the weak SDE methods
#'
#' The discrete increments that replace the Brownian increment in the weak
#' methods, computed from a Gaussian sample so they keep its moments:
#' `calc_twopoint_random(sqdt, dW)` is `sqdt` where `dW > 0` and `-sqdt`
#' otherwise; `calc_threepoint_random(sq3dt, quantile, dW_scaled)` is `0`
#' when `|dW_scaled| <= -quantile` and `+-sq3dt` by the sign beyond it, with
#' `quantile` the standard normal 1/6 quantile, so the outcomes have the
#' probabilities 1/6, 2/3, 1/6 of a three point variable.
#'
#' @param sqdt,sq3dt `sqrt(dt)` and `sqrt(3 dt)`.
#' @param dW,dW_scaled Gaussian increments (`dW_scaled = dW / sqrt(dt)`).
#' @param quantile The standard normal quantile `qnorm(1/6)`.
#' @return A numeric vector.
#' @examples
#' calc_twopoint_random(0.1, c(-0.3, 0.2))
#' calc_threepoint_random(sqrt(0.03), qnorm(1 / 6), c(-2, 0.1, 1.5))
#' @name weak_random_variables
NULL

#' @rdname weak_random_variables
#' @export
calc_twopoint_random <- function(sqdt, dW) ifelse(dW > 0, sqdt, -sqdt)

#' @rdname weak_random_variables
#' @export
calc_threepoint_random <- function(sq3dt, quantile, dW_scaled)
  ifelse(abs(dW_scaled) > -quantile, ifelse(dW_scaled < quantile, -sq3dt, sq3dt), 0)
