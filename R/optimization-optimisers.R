# The rules of Optimisers.jl (src/rules.jl, v0.4.11) and the solver loop of
# OptimizationOptimisers (lib/OptimizationOptimisers/src/OptimizationOptimisers.jl,
# v0.3.26): gradient based update rules that run for a fixed number of
# iterations and return the best point evaluated.

.optimisers_rule <- function(name, ...) {
  r <- .opt_alg(name, ...)
  class(r) <- c(name, "AbstractRule", "OptimizationAlgorithm")
  r$solver <- .optimisers_solve
  r$traits <- list(grad = TRUE)
  r
}

#' @export
print.AbstractRule <- function(x, ...) {
  if (inherits(x, "OptimiserChain")) {
    cat("OptimiserChain(", sep = "")
    for (i in seq_along(x$opts)) {
      if (i > 1L) cat(", ")
      cat(x$opts[[i]]$name, "()", sep = "")
    }
    cat(")\n")
    return(invisible(x))
  }
  pars <- x[setdiff(names(x), c("name", "solver", "traits"))]
  vals <- vapply(pars, function(v) paste(format(v, digits = 6), collapse = ", "), "")
  cat(x$name, "(", paste(names(vals), vals, sep = " = ", collapse = ", "), ")\n", sep = "")
  invisible(x)
}

#' Optimisers.jl update rules
#'
#' The update rules of Optimisers.jl, used through OptimizationOptimisers:
#' `solve(prob, Adam(0.05), maxiters = 500)`. Each iteration evaluates the
#' gradient and the objective at the current point and moves it by the
#' rule's step; `maxiters` (or `epochs`) is required. With `save_best = TRUE`
#' (the default) the solution is the best point evaluated, and the last
#' iteration only evaluates. An iteration whose gradient is not finite does
#' not move the point.
#'
#' * `Descent(eta)`: `x - eta * g`. The default `eta` is the single
#'   precision 0.1 of Julia, 0.10000000149011612.
#' * `Momentum(eta, rho)`: velocity `v = rho * v + eta * g`, step `v`.
#' * `Nesterov(eta, rho)`: Nesterov momentum.
#' * `RMSProp(eta, rho, epsilon, centred)`: scaling by the root mean square
#'   of the gradients (when centered: of their variance).
#' * `Rprop(eta, ell, gamma)`: steps by the sign of the gradient with
#'   adaptive sizes, multiplied by `ell` and kept within `gamma`.
#' * `Adam(eta, beta, epsilon)`, `AdaMax`, `NAdam`, `RAdam`, `OAdam`,
#'   `AMSGrad`, `AdaBelief`: the Adam family.
#' * `AdamW(eta, beta, lambda, epsilon, couple)`: Adam with weight decay
#'   (coupled with the learning rate as in pytorch, or decoupled).
#' * `Lion(eta, beta)`: steps by the sign of an interpolated momentum.
#' * `AdaGrad(eta, epsilon)`, `AdaDelta(rho, epsilon)`.
#' * `WeightDecay(lambda)` adds `lambda * x` to the gradient and
#'   `SignDecay(lambda)` adds `lambda * sign(x)`; `ClipGrad(delta)` clamps
#'   every component to `[-delta, delta]` and `ClipNorm(omega, p)` rescales
#'   a gradient whose `p` norm exceeds `omega`. These are meant to be
#'   composed with `OptimiserChain(...)`, which applies its rules in order.
#' * `AccumGrad(n)` in a chain accumulates the mean of `n` gradients before
#'   the later rules see it; in between the point does not move.
#'
#' `ADAM`, `ADAMW`, `NADAM`, `RADAM`, `OADAM`, `ADAGrad` and `ADADelta` are
#' the deprecated Optimisers.jl spellings of the same rules. `Optimisers` is
#' the module: a list of the rules and of `setup(rule, x)` (the initial
#' state), `update(state, x, dx)` (returns `list(state, x)`) and
#' `apply(rule, state, x, dx)` (returns `list(state, dx)`, the step to
#' subtract).
#'
#' @param eta Learning rate.
#' @param rho Momentum or decay factor.
#' @param beta Decay of the first and second moment estimates.
#' @param epsilon Constant against division by zero.
#' @param centred,centered Use the centered variant of RMSProp (`centred` is
#'   the name of the Julia keyword).
#' @param ell Decrease and increase factors of Rprop.
#' @param gamma Smallest and largest step of Rprop.
#' @param lambda Weight decay strength.
#' @param couple Couple the weight decay of AdamW with the learning rate.
#' @param delta Clipping bound of each component.
#' @param omega Largest norm of the gradient.
#' @param p Order of the norm.
#' @param throw Signal an error when the norm is not finite.
#' @param n Number of gradients to accumulate.
#' @param ... Rules of the chain, or the arguments of the rule.
#' @return A rule (an optimization algorithm).
#' @examples
#' f <- function(u, p) (1 - u[1])^2 + 100 * (u[2] - u[1]^2)^2
#' g <- function(u, p) c(-2 * (1 - u[1]) - 400 * u[1] * (u[2] - u[1]^2), 200 * (u[2] - u[1]^2))
#' prob <- OptimizationProblem(f, c(-1.2, 1), grad = g)
#' solve(prob, Adam(0.05), maxiters = 300)$u
#' solve(prob, OptimiserChain(ClipNorm(1), Descent(1e-3)), maxiters = 100)$objective
#' @name optimisers_rules
#' @export
Descent <- function(eta = 0.100000001490116119384765625) .optimisers_rule("Descent", eta = eta)
#' @rdname optimisers_rules
#' @export
Momentum <- function(eta = 0.01, rho = 0.9) .optimisers_rule("Momentum", eta = eta, rho = rho)
#' @rdname optimisers_rules
#' @export
Nesterov <- function(eta = 0.001, rho = 0.9) .optimisers_rule("Nesterov", eta = eta, rho = rho)
#' @rdname optimisers_rules
#' @export
RMSProp <- function(eta = 0.001, rho = 0.9, epsilon = 1e-8, centred = FALSE, centered = FALSE) {
  if (eta < 0) stop("the learning rate cannot be negative", call. = FALSE)
  .optimisers_rule("RMSProp", eta = eta, rho = rho, epsilon = epsilon, centred = isTRUE(centred) || isTRUE(centered))
}
#' @rdname optimisers_rules
#' @export
Rprop <- function(eta = 1e-3, ell = c(0.5, 1.2), gamma = c(1e-6, 50)) .optimisers_rule("Rprop", eta = eta, ell = ell,
                                                                                      gamma = gamma)
#' @rdname optimisers_rules
#' @export
Adam <- function(eta = 0.001, beta = c(0.9, 0.999), epsilon = 1e-8)
  .optimisers_rule("Adam", eta = eta, beta = beta, epsilon = epsilon)
#' @rdname optimisers_rules
#' @export
Lion <- function(eta = 0.001, beta = c(0.9, 0.999)) .optimisers_rule("Lion", eta = eta, beta = beta)
#' @rdname optimisers_rules
#' @export
RAdam <- function(eta = 0.001, beta = c(0.9, 0.999), epsilon = 1e-8)
  .optimisers_rule("RAdam", eta = eta, beta = beta, epsilon = epsilon)
#' @rdname optimisers_rules
#' @export
AdaMax <- function(eta = 0.001, beta = c(0.9, 0.999), epsilon = 1e-8)
  .optimisers_rule("AdaMax", eta = eta, beta = beta, epsilon = epsilon)
#' @rdname optimisers_rules
#' @export
OAdam <- function(eta = 0.001, beta = c(0.5, 0.9), epsilon = 1e-8)
  .optimisers_rule("OAdam", eta = eta, beta = beta, epsilon = epsilon)
#' @rdname optimisers_rules
#' @export
AdaGrad <- function(eta = 0.1, epsilon = 1e-8) .optimisers_rule("AdaGrad", eta = eta, epsilon = epsilon)
#' @rdname optimisers_rules
#' @export
AdaDelta <- function(rho = 0.9, epsilon = 1e-8) .optimisers_rule("AdaDelta", rho = rho, epsilon = epsilon)
#' @rdname optimisers_rules
#' @export
AMSGrad <- function(eta = 0.001, beta = c(0.9, 0.999), epsilon = 1e-8)
  .optimisers_rule("AMSGrad", eta = eta, beta = beta, epsilon = epsilon)
#' @rdname optimisers_rules
#' @export
NAdam <- function(eta = 0.001, beta = c(0.9, 0.999), epsilon = 1e-8)
  .optimisers_rule("NAdam", eta = eta, beta = beta, epsilon = epsilon)
#' @rdname optimisers_rules
#' @export
AdamW <- function(eta = 0.001, beta = c(0.9, 0.999), lambda = 0, epsilon = 1e-8, couple = TRUE) {
  if (eta < 0) stop("the learning rate cannot be negative", call. = FALSE)
  .optimisers_rule("AdamW", eta = eta, beta = beta, lambda = lambda, epsilon = epsilon, couple = isTRUE(couple))
}
#' @rdname optimisers_rules
#' @export
AdaBelief <- function(eta = 0.001, beta = c(0.9, 0.999), epsilon = 1e-16)
  .optimisers_rule("AdaBelief", eta = eta, beta = beta, epsilon = epsilon)
#' @rdname optimisers_rules
#' @export
WeightDecay <- function(lambda = 5e-4) .optimisers_rule("WeightDecay", lambda = lambda)
#' @rdname optimisers_rules
#' @export
SignDecay <- function(lambda = 1e-3) .optimisers_rule("SignDecay", lambda = lambda)
#' @rdname optimisers_rules
#' @export
ClipGrad <- function(delta = 10) .optimisers_rule("ClipGrad", delta = delta)
#' @rdname optimisers_rules
#' @export
ClipNorm <- function(omega = 10, p = 2, throw = TRUE) .optimisers_rule("ClipNorm", omega = omega, p = p,
                                                                       throw = isTRUE(throw))
#' @rdname optimisers_rules
#' @export
OptimiserChain <- function(...) {
  opts <- list(...)
  if (!all(vapply(opts, inherits, TRUE, "AbstractRule"))) stop("OptimiserChain() takes Optimisers rules", call. = FALSE)
  .optimisers_rule("OptimiserChain", opts = opts)
}
#' @rdname optimisers_rules
#' @export
AccumGrad <- function(n) {
  if (!(n > 0)) stop("AccumGrad must accumulate at least one gradient", call. = FALSE)
  r <- .optimisers_rule("AccumGrad")
  r$n <- as.integer(n) # set after: `n =` would partially match the `name` argument
  r
}

#' @rdname optimisers_rules
#' @export
ADAM <- function(...) Adam(...)
#' @rdname optimisers_rules
#' @export
ADAMW <- function(...) AdamW(...)
#' @rdname optimisers_rules
#' @export
NADAM <- function(...) NAdam(...)
#' @rdname optimisers_rules
#' @export
RADAM <- function(...) RAdam(...)
#' @rdname optimisers_rules
#' @export
OADAM <- function(...) OAdam(...)
#' @rdname optimisers_rules
#' @export
ADAGrad <- function(...) AdaGrad(...)
#' @rdname optimisers_rules
#' @export
ADADelta <- function(...) AdaDelta(...)

# init(rule, x) of Optimisers.jl
.optimisers_init <- function(o, x) {
  z <- numeric(length(x))
  switch(o$name,
    Descent = , WeightDecay = , SignDecay = , ClipGrad = , ClipNorm = NULL,
    Momentum = , Nesterov = , Lion = z,
    RMSProp = list(quad = z, lin = z),
    Rprop = list(g = z, eta = rep(o$eta, length(x))),
    Adam = , AdaMax = , NAdam = , AdamW = , AdaBelief = list(m = z, v = z, bt = o$beta),
    RAdam = list(m = z, v = z, bt = o$beta, t = 1L),
    OAdam = list(m = z, v = z, bt = o$beta, term = z),
    AdaGrad = rep(o$epsilon, length(x)),
    AdaDelta = list(acc = z, dacc = z),
    AMSGrad = list(m = rep(o$epsilon, length(x)), v = rep(o$epsilon, length(x)), vhat = rep(o$epsilon, length(x))),
    OptimiserChain = lapply(o$opts, .optimisers_init, x = x),
    AccumGrad = list(acc = z, counter = 1L),
    stop("unknown rule ", o$name, call. = FALSE))
}

# norm(dx, p) as Optimisers' ClipNorm computes it
.optimisers_norm <- function(dx, p) {
  if (p == 2) sqrt(sum(dx * dx))
  else if (p == 1) sum(abs(dx))
  else if (p == Inf) max(abs(dx))
  else if (p == 0) as.numeric(sum(dx != 0))
  else if (p == -Inf) min(abs(dx))
  else sum(abs(dx)^p)^(1 / p)
}

# apply!(rule, state, x, dx) of Optimisers.jl: list(state, step); a NULL
# step (AccumGrad between updates) leaves x unchanged.
.optimisers_apply <- function(o, s, x, dx) {
  b1 <- o$beta[1L]; b2 <- o$beta[2L]
  switch(o$name,
    Descent = list(s, dx * o$eta),
    Momentum = {
      s <- o$rho * s + o$eta * dx
      list(s, s)
    },
    Nesterov = {
      newdx <- -o$rho^2 * s + (1 + o$rho) * o$eta * dx
      list(o$rho * s - o$eta * dx, newdx)
    },
    RMSProp = {
      s$quad <- o$rho * s$quad + (1 - o$rho) * (dx * dx)
      if (o$centred) {
        s$lin <- o$rho * s$lin + (1 - o$rho) * dx
        list(s, dx * o$eta / (sqrt(s$quad - s$lin * s$lin) + o$epsilon))
      } else list(s, dx * o$eta / (sqrt(s$quad) + o$epsilon))
    },
    Rprop = {
      gd <- s$g * dx
      eta <- ifelse(gd > 0, pmin(s$eta * o$ell[2L], o$gamma[2L]),
                    ifelse(gd < 0, pmax(s$eta * o$ell[1L], o$gamma[1L]), s$eta))
      g <- ifelse(gd < 0, 0, dx)
      list(list(g = g, eta = eta), eta * sign(g))
    },
    Adam = {
      m <- b1 * s$m + (1 - b1) * dx
      v <- b2 * s$v + (1 - b2) * (dx * dx)
      step <- m / (1 - s$bt[1L]) / (sqrt(v / (1 - s$bt[2L])) + o$epsilon) * o$eta
      list(list(m = m, v = v, bt = s$bt * o$beta), step)
    },
    Lion = {
      s <- b2 * dx + (1 - b2) * s
      list(s, o$eta * sign((b2 - b1) * dx + b1 * s))
    },
    RAdam = {
      rinf <- 2 / (1 - b2) - 1
      m <- b1 * s$m + (1 - b1) * dx
      v <- b2 * s$v + (1 - b2) * (dx * dx)
      rho <- rinf - 2 * s$t * s$bt[2L] / (1 - s$bt[2L])
      step <- if (rho > 4) {
        r <- sqrt((rho - 4) * (rho - 2) * rinf / ((rinf - 4) * (rinf - 2) * rho))
        m / (1 - s$bt[1L]) / (sqrt(v / (1 - s$bt[2L])) + o$epsilon) * o$eta * r
      } else m / (1 - s$bt[1L]) * o$eta
      list(list(m = m, v = v, bt = s$bt * o$beta, t = s$t + 1L), step)
    },
    AdaMax = {
      m <- b1 * s$m + (1 - b1) * dx
      u <- pmax(b2 * s$v, abs(dx))
      list(list(m = m, v = u, bt = s$bt * o$beta), (o$eta / (1 - s$bt[1L])) * m / (u + o$epsilon))
    },
    OAdam = {
      m <- b1 * s$m + (1 - b1) * dx
      v <- b2 * s$v + (1 - b2) * (dx * dx)
      term <- o$eta * m / (1 - s$bt[1L]) / (sqrt(v / (1 - s$bt[2L])) + o$epsilon)
      list(list(m = m, v = v, bt = s$bt * o$beta, term = term), 2 * term - s$term)
    },
    AdaGrad = {
      acc <- s + dx * dx
      list(acc, dx * o$eta / (sqrt(acc) + o$epsilon))
    },
    AdaDelta = {
      acc <- o$rho * s$acc + (1 - o$rho) * (dx * dx)
      step <- dx * sqrt(s$dacc + o$epsilon) / sqrt(acc + o$epsilon)
      list(list(acc = acc, dacc = o$rho * s$dacc + (1 - o$rho) * (step * step)), step)
    },
    AMSGrad = {
      m <- b1 * s$m + (1 - b1) * dx
      v <- b2 * s$v + (1 - b2) * (dx * dx)
      vhat <- pmax(s$vhat, v)
      list(list(m = m, v = v, vhat = vhat), o$eta * m / (sqrt(vhat) + o$epsilon))
    },
    NAdam = {
      m <- b1 * s$m + (1 - b1) * dx
      v <- b2 * s$v + (1 - b2) * (dx * dx)
      step <- (b1 * m / (1 - b1 * s$bt[1L]) + (1 - b1) * dx / (1 - s$bt[1L])) /
        (sqrt(v * b2 / (1 - s$bt[2L])) + o$epsilon) * o$eta
      list(list(m = m, v = v, bt = s$bt * o$beta), step)
    },
    AdamW = {
      m <- b1 * s$m + (1 - b1) * dx
      v <- b2 * s$v + (1 - b2) * (dx * dx)
      step <- m / (1 - s$bt[1L]) / (sqrt(v / (1 - s$bt[2L])) + o$epsilon)
      step <- if (o$couple) o$eta * (step + o$lambda * x) else o$eta * step + o$lambda * x
      list(list(m = m, v = v, bt = s$bt * o$beta), step)
    },
    AdaBelief = {
      m <- b1 * s$m + (1 - b1) * dx
      d <- dx - m
      st <- b2 * s$v + (1 - b2) * (d * d) + o$epsilon
      list(list(m = m, v = st, bt = s$bt * o$beta), o$eta * m / (1 - s$bt[1L]) / (sqrt(st / (1 - s$bt[2L])) + o$epsilon))
    },
    WeightDecay = list(s, dx + o$lambda * x),
    SignDecay = list(s, dx + o$lambda * sign(x)),
    ClipGrad = list(s, pmin(pmax(dx, -o$delta), o$delta)),
    ClipNorm = {
      nrm <- .optimisers_norm(dx, o$p)
      if (o$throw && !is.finite(nrm)) stop("gradient has ", o$p, " norm ", nrm, call. = FALSE)
      list(s, dx * min(o$omega / nrm, 1))
    },
    OptimiserChain = {
      for (i in seq_along(o$opts)) {
        if (is.null(dx)) break
        r <- .optimisers_apply(o$opts[[i]], s[[i]], x, dx)
        if (!is.null(r[[1L]])) s[[i]] <- r[[1L]]
        dx <- r[[2L]]
      }
      list(s, dx)
    },
    AccumGrad = {
      acc <- if (s$counter == 1L) dx / o$n else s$acc + dx / o$n
      if (s$counter == o$n) list(list(acc = acc, counter = 1L), acc)
      else list(list(acc = acc, counter = s$counter + 1L), NULL)
    },
    stop("unknown rule ", o$name, call. = FALSE))
}

# Optimisers.update(state, x, dx): list(state, x)
.optimisers_update <- function(rule, state, x, dx) {
  r <- .optimisers_apply(rule, state, x, dx)
  list(r[[1L]], if (is.null(r[[2L]])) x else x - r[[2L]])
}

#' @rdname optimisers_rules
#' @export
Optimisers <- list(
  Descent = Descent, Momentum = Momentum, Nesterov = Nesterov, RMSProp = RMSProp, Rprop = Rprop, Adam = Adam,
  Lion = Lion, RAdam = RAdam, AdaMax = AdaMax, OAdam = OAdam, AdaGrad = AdaGrad, AdaDelta = AdaDelta,
  AMSGrad = AMSGrad, NAdam = NAdam, AdamW = AdamW, AdaBelief = AdaBelief, WeightDecay = WeightDecay,
  SignDecay = SignDecay, ClipGrad = ClipGrad, ClipNorm = ClipNorm, OptimiserChain = OptimiserChain,
  AccumGrad = AccumGrad, ADAM = ADAM, ADAMW = ADAMW, NADAM = NADAM, RADAM = RADAM, OADAM = OADAM,
  ADAGrad = ADAGrad, ADADelta = ADADelta,
  setup = function(rule, x) .optimisers_init(rule, x),
  update = function(state, x, dx, rule) .optimisers_update(rule, state, x, dx),
  apply = function(rule, state, x, dx) .optimisers_apply(rule, state, x, dx))

# The objective and gradient of a problem as plain functions of the
# unknowns as a vector; the functions of the problem receive them in the
# shape of u0 (a matrix u0 gives matrices, as in Julia)
.opt_fns <- function(prob) {
  fn <- prob$f
  p <- prob$p
  sh <- function(x) .opt_shape(x, prob)
  f <- function(x) as.numeric(fn$f(sh(x), p))[1L]
  g <- if (is.function(fn$grad)) function(x) as.numeric(fn$grad(sh(x), p)) else function(x) .opt_fd_gradient(f, x)
  list(f = f, g = g, p = p, fn = fn)
}

# the unknowns in the shape of u0
.opt_shape <- function(x, prob) {
  d <- dim(prob$u0)
  if (!is.null(d) && length(x) == prod(d)) dim(x) <- d
  x
}

.opt_maxiters <- function(m) {
  if (is.null(m)) return(NULL)
  if (m <= 0) stop("The number of maxiters has to be a non-negative and non-zero number.", call. = FALSE)
  as.integer(round(m))
}

# __solve of OptimizationOptimisers
.optimisers_solve <- function(prob, rule, opts) {
  fs <- .opt_fns(prob)
  if (is.function(prob$f$cons)) stop(rule$name, "() does not handle constraints", call. = FALSE)
  ndata <- 1L
  epochs <- opts$epochs; maxiters <- opts$maxiters
  if (is.null(maxiters) && is.null(epochs))
    stop("The number of iterations must be specified with either the epochs or maxiters kwarg. ",
         "Where maxiters = epochs * length(data).", call. = FALSE)
  if (!is.null(maxiters) && !is.null(epochs) && maxiters != epochs * ndata)
    stop("Both maxiters and epochs were passed but maxiters != epochs * length(data).", call. = FALSE)
  if (is.null(maxiters)) maxiters <- epochs * ndata
  if (is.null(epochs)) epochs <- maxiters / ndata
  epochs <- .opt_maxiters(epochs); maxiters <- .opt_maxiters(maxiters)
  save_best <- opts$save_best %||% TRUE
  theta <- as.numeric(prob$u0)
  x <- Inf
  min_err <- Inf; min_theta <- theta
  state <- .optimisers_init(rule, theta)
  iterations <- 0L; fevals <- 0L; gevals <- 0L
  t0 <- proc.time()[["elapsed"]]
  for (epoch in seq_len(epochs)) {
    G <- fs$g(theta)
    x <- fs$f(theta)
    iterations <- iterations + 1L
    fevals <- fevals + 2L
    gevals <- gevals + 1L
    if (save_best) {
      if (x < min_err) {
        min_err <- x
        min_theta <- theta
      }
      if (iterations == ndata * epochs) {
        x <- min_err
        theta <- min_theta
        G <- fs$g(theta)
        break
      }
    }
    # a non finite gradient skips the update
    if (all(is.finite(G))) {
      r <- .optimisers_update(rule, state, theta, G)
      state <- r[[1L]]; theta <- r[[2L]]
    } else if (!identical(opts$verbose, FALSE)) {
      warning("Skipping parameter update due to NaN or Inf in gradients at iteration ", iterations, call. = FALSE)
    }
  }
  if (iterations == 0L) {
    x <- fs$f(as.numeric(prob$u0))
    fevals <- fevals + 1L
  }
  sol <- .opt_solution(prob, rule, theta, x, "Default", list(state = state))
  .opt_set_fields(sol, stats = OptimizationStats(iterations = iterations, time = proc.time()[["elapsed"]] - t0,
                                                fevals = fevals, gevals = gevals))
}
