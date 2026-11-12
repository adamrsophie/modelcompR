# OptimizationNOMAD: # NOMADOpt() runs a mesh adaptive direct search (MADS,
# Audet and Dennis 2006) in native R, standing in for the NOMAD 4 library
# that NOMAD.jl calls: a quadratic model search, orthogonal 2n poll
# directions from Householder reflections, opportunistic polling with a
# speculative search along the last successful direction, a cache of
# evaluated points, the extreme barrier for constraints and a first phase
# that minimizes the constraint violation from an infeasible start.

#' Mesh adaptive direct search
#'
#' `NOMADOpt()` (OptimizationNOMAD.jl) minimizes a black box objective with
#' the mesh adaptive direct search method of NOMAD, subject to the bounds
#' of the problem and to `lcons <= cons(u) <= ucons`. Each constraint enters
#' as `max(lcons - cons(u), cons(u) - ucons) <= 0`. The native method uses
#' a quadratic model search after each failed poll (quadratic regression models of the
#' objective and the constraints on the cached points near the iterate,
#' minimized in a trust box) and 2n orthogonal poll directions, doubles the
#' frame after a success, halves it after a failure and stops when the
#' frame reaches the precision of the iterate, after `maxiters` (NOMAD's
#' `max_bb_eval`, default 20000) evaluations or after `maxtime` seconds.
#' NOMAD's Nelder Mead search is not part of the native method, so iterates
#' and counts differ from NOMAD while the solutions agree.
#'
#' Options of `solve()`: `cons_method` (`"ExtremeBarrierMethod"`, the
#' default, or `"ProgressiveBarrierMethod"`, which the native method treats
#' as the extreme barrier after its first phase), and the NOMAD options
#' `max_bb_eval`, `opportunistic_eval`, `use_cache` and the display options,
#' which have no effect. `abstol` and `reltol` are not used, as in
#' OptimizationNOMAD.jl.
#'
#' @return An algorithm object for [solve()].
#' @examples
#' prob <- OptimizationProblem(function(u, p) (u[1] - 1) * (u[1] - 1) + (u[2] + 2) * (u[2] + 2),
#'                             c(0, 0), numeric())
#' solve(prob, NOMADOpt())$u
#' @export
NOMADOpt <- function() {
  structure(list(name = "NOMADOpt", traits = list(callback = FALSE, bounds = TRUE, cons = TRUE),
                 solve = .nomad_solve),
            class = c("NOMADOpt", "OptimizationAlgorithm"))
}

.nomad_options <- c("display_degree", "display_all_eval", "display_infeasible", "display_unsuccessful",
                    "max_bb_eval", "opportunistic_eval", "use_cache", "lh_search", "speculative_search", "nm_search",
                    "cons_method", "maxtime")

.nomad_solve <- function(a, alg, opts) {
  fn <- a$f
  p <- a$p
  x0 <- as.numeric(a$u0)
  n <- length(x0)
  extra <- setdiff(names(opts), c("maxiters", "abstol", "reltol", "callback", "progress", "verbose", .nomad_options))
  if (length(extra)) stop("type NomadOptions has no field ", extra[1L], call. = FALSE)
  method <- opts$cons_method %||% "ExtremeBarrierMethod"
  if (!method %in% c("ExtremeBarrierMethod", "ProgressiveBarrierMethod"))
    stop("cons_method must be ExtremeBarrierMethod or ProgressiveBarrierMethod", call. = FALSE)
  max_eval <- as.integer(opts$maxiters %||% opts$max_bb_eval %||% 20000L)
  deadline <- if (is.null(opts$maxtime)) Inf else proc.time()[["elapsed"]] + opts$maxtime
  lb <- if (is.null(a$lb)) rep(-Inf, n) else rep_len(as.numeric(a$lb), n)
  ub <- if (is.null(a$ub)) rep(Inf, n) else rep_len(as.numeric(a$ub), n)
  m <- length(a$ucons %||% numeric())
  cons <- if (is.function(fn$cons) && m) function(x) as.numeric(fn$cons(x, p))
  cache <- new.env(hash = TRUE)
  pts <- list()
  neval <- 0L
  bb <- function(x) {
    key <- paste(sprintf("%a", x), collapse = ",")
    hit <- cache[[key]]
    if (!is.null(hit)) return(hit)
    neval <<- neval + 1L
    f <- as.numeric(fn$f(x, p))[1L]
    cv <- if (is.null(cons)) numeric() else cons(x)
    h <- sum(pmax(0, a$lcons - cv, cv - a$ucons)^2)
    out <- c(if (is.finite(f)) f else Inf, h)
    if (is.finite(f) && all(is.finite(cv))) pts[[length(pts) + 1L]] <<- c(x, f, cv)
    assign(key, out, envir = cache)
    out
  }
  x0 <- pmin(pmax(x0, lb), ub)
  frame <- ifelse(is.finite(lb) & is.finite(ub), (ub - lb) / 10, ifelse(x0 != 0, abs(x0) / 10, 1))
  frame[frame <= 0] <- 1
  primes <- sqrt(c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53))
  dirs <- function(k) {
    v <- ((k * primes[(seq_len(n) - 1L) %% 16L + 1L] + seq_len(n) / (n + 1)) %% 1) - 0.5
    H <- diag(n) - 2 * tcrossprod(v) / sum(v * v)
    H <- sweep(H, 2L, apply(abs(H), 2L, max), "/")
    cbind(H, -H)
  }
  # phase one minimizes the violation h; then the objective under the barrier
  score <- function(o, phase) if (phase == 1L) o[2L] else if (o[2L] > 0) Inf else o[1L]
  x <- x0
  ox <- bb(x)
  phase <- if (ox[2L] > 0) 1L else 2L
  last <- NULL
  k <- 0L
  failed <- FALSE
  repeat {
    if (all(frame < pmax(1e-13, 1e-13 * abs(x)))) { stop_reason <- "precision"; break }
    if (neval >= max_eval) { stop_reason <- "max_bb_eval"; break }
    if (proc.time()[["elapsed"]] > deadline) { stop_reason <- "max_time"; break }
    k <- k + 1L
    cur <- score(ox, phase)
    trials <- if (is.null(last)) list() else list(pmin(pmax(x + last, lb), ub))
    if (phase == 2L && failed) {
      ms <- .nomad_model_search(pts, x, 2 * max(frame), n, a, lb, ub)
      if (!is.null(ms)) trials <- c(list(structure(ms, model = TRUE)), trials)
    }
    D <- dirs(k)
    for (j in seq_len(ncol(D))) trials[[length(trials) + 1L]] <- pmin(pmax(x + frame * D[, j], lb), ub)
    success <- FALSE
    for (t in trials) {
      if (identical(t, x)) next
      if (neval >= max_eval) break
      by_model <- isTRUE(attr(t, "model"))
      t <- as.numeric(t)
      ot <- bb(t)
      if (score(ot, phase) < cur) {
        last <- if (!by_model) t - x
        x <- t; ox <- ot
        success <- TRUE
        break
      }
    }
    failed <- !success
    if (success && by_model) {
      failed <- TRUE
    } else if (success) {
      frame <- frame * 2
      if (phase == 1L && ox[2L] == 0) phase <- 2L
    } else {
      frame <- frame / 2
      last <- NULL
    }
  }
  feasible <- ox[2L] == 0
  retcode <- if (feasible) "Default" else "Infeasible"
  .opt_solution(a, alg, if (feasible) x else rep(NaN, n), if (feasible) ox[1L] else NaN, retcode,
                list(x_best_feas = if (feasible) x, bbo_best_feas = if (feasible) ox[1L],
                     x_best_inf = if (!feasible) x, h_best_inf = if (!feasible) ox[2L], nb_evals = neval,
                     iterations = k, stop_reason = stop_reason))
}

# quadratic model search: regression models of the objective and of the
# constraints on the cached points within `r` of x, minimized over the box
# of half width r by COBYLA, or by L-BFGS-B without constraints
.nomad_model_search <- function(pts, x, r, n, a, lb, ub) {
  q <- (n + 1L) * (n + 2L) / 2L
  if (length(pts) < q || !is.finite(r) || r <= 0) return(NULL)
  P <- do.call(rbind, utils::tail(pts, 10L * q))
  Z <- sweep(P[, seq_len(n), drop = FALSE], 2L, x) / r
  near <- apply(abs(Z), 1L, max) <= 1
  if (sum(near) < q) return(NULL)
  Z <- Z[near, , drop = FALSE]
  Y <- P[near, -seq_len(n), drop = FALSE]
  ij <- which(upper.tri(diag(n), diag = TRUE), arr.ind = TRUE)
  basis <- function(z) c(1, z, z[ij[, 1L]] * z[ij[, 2L]])
  B <- t(apply(Z, 1L, basis))
  qrB <- qr(B)
  if (qrB$rank < q) return(NULL)
  coef <- qr.coef(qrB, Y)
  model <- function(z) as.numeric(crossprod(coef, basis(z)))
  m <- ncol(Y) - 1L
  lo <- pmax(-1, (lb - x) / r); hi <- pmin(1, (ub - x) / r)
  z0 <- pmin(pmax(numeric(n), lo), hi)
  sol <- if (m) {
    lc <- a$lcons; uc <- a$ucons
    keep <- is.finite(c(lc, uc))
    calcfc <- function(z) { v <- model(z); list(f = v[1L], c = c(lc - v[-1L], v[-1L] - uc)[keep]) }
    tryCatch({
      r0 <- calcfc(z0)
      rc <- .pr_cobyla(calcfc, length(r0$c), z0, xl = lo, xu = hi, f0 = r0$f, nlconstr0 = r0$c, rhobeg = 0.5,
                      rhoend = 1e-6, maxfun = 30L * n)
      list(u = rc$x)
    }, error = function(e) NULL)
  } else {
    tryCatch(list(u = stats::optim(z0, function(w) model(w)[1L], method = "L-BFGS-B", lower = lo, upper = hi)$par),
             error = function(e) NULL)
  }
  if (is.null(sol) || !all(is.finite(sol$u))) return(NULL)
  pmin(pmax(x + r * as.numeric(sol$u), lb), ub)
}
