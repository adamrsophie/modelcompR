# Composite SDE algorithms with automatic stiffness switching, after
# StochasticDiffEqCore (src/algorithms.jl: StochasticCompositeAlgorithm;
# src/composite_algs.jl: AutoSwitch, AutoSwitchCache, is_stiff,
# AutoAlgSwitch; src/perform_step/composite.jl; src/alg_utils.jl:
# alg_order, isadaptive and alg_stability_size of composites) and
# StochasticDiffEq.jl (src/StochasticDiffEq.jl: AutoSOSRI2, AutoSOSRA2),
# with the generic choose_algorithm! of OrdinaryDiffEqCore
# (perform_step/composite_perform_step.jl).
#
# A composite runs as one method of the core driver, so it has one step
# size controller, as the PIController of a StochasticDiffEq composite (its
# gains come from the largest order of the members). Before every step
# attempt the driver calls the runner's pre_step hook with the step size the
# integrator holds at that point, which is where loopheader! calls
# choose_algorithm!: the choice may switch the member and scale the step
# size. Each member keeps its cache for the whole solve; the stiffness
# estimate eigen_est is set by the four stage SRI step when the first member
# is SOSRI2 and by the three stage SRA step when it is SOSRA2, as in
# sri.jl and sra.jl.

#' Composite SDE algorithms and automatic stiffness switching
#'
#' `StochasticCompositeAlgorithm(algs, choice_function)` steps with one of
#' the SDE algorithms in the list `algs`, chosen before every step by
#' `choice_function`. The choice is either an [AutoSwitch()] object (the
#' switch between its nonstiff and stiff methods) or a function of one
#' argument returning the index of the method; the argument is a list with
#' the `integrator`, its `t`, `u`, `p`, the step count `iter`, the step size
#' `dt` and the stiffness estimate `eigen_est`.
#'
#' `AutoAlgSwitch(nonstiffalg, stiffalg, ...)` is the composite that starts
#' with `nonstiffalg` (or `stiffalg` with `stiffalgfirst = TRUE`) and
#' switches to `stiffalg` after more than `maxstiffstep` consecutive steps
#' whose stiffness ratio `abs(eigen_est * dt / alg_stability_size(nonstiffalg))`
#' exceeds `nonstifftol`, and back after more than `maxnonstiffstep`
#' consecutive steps below `stifftol`; the step size is multiplied by
#' `dtfac` when switching to the stiff method and divided by it when
#' switching back. `switch_max` is the number of nonstiff verdicts after
#' which StochasticDiffEq checks the solution for errors again. The
#' stiffness estimate comes from the steps of `SOSRI2()` and `SOSRA2()` when
#' they are the first method: `AutoSOSRI2(alg, ...)` and `AutoSOSRA2(alg,
#' ...)` pair them with a stiff method `alg`. A method without a stability
#' region size counts as size 0, so any nonzero estimate is stiff for it.
#'
#' The solution records the index of the method of each saved point in
#' `sol$alg_choice`.
#'
#' @param algs A list of SDE algorithms.
#' @param choice_function An [AutoSwitch()] object or a function, see
#'   Details.
#' @param nonstiffalg,stiffalg SDE algorithms.
#' @param alg A stiff SDE algorithm.
#' @param maxstiffstep,maxnonstiffstep Consecutive verdicts before switching.
#' @param nonstifftol,stifftol Stiffness ratio tolerances.
#' @param dtfac Step size factor applied at a switch.
#' @param stiffalgfirst Start with the stiff method.
#' @param switch_max Nonstiff verdicts after which errors are checked again.
#' @param ... Options of `AutoAlgSwitch()`.
#' @return An SDE algorithm of class `StochasticCompositeAlgorithm`.
#' @examples
#' f <- function(u, p, t) -50 * u
#' g <- function(u, p, t) 0.1 * u
#' prob <- SDEProblem(f, g, 1, c(0, 0.5))
#' sol <- solve(prob, AutoSOSRI2(LambaEM()), dt = 0.01, seed = 1)
#' table(sol$alg_choice)
#' @name sde_composite
NULL

.sde_autoswitch_spec <- function(nonstiffalg, stiffalg, maxstiffstep = 10L, maxnonstiffstep = 3L,
                                 nonstifftol = 9 / 10, stifftol = 9 / 10, dtfac = 2, stiffalgfirst = FALSE,
                                 switch_max = 5L)
  structure(list(nonstiffalg = nonstiffalg, stiffalg = stiffalg, maxstiffstep = maxstiffstep,
                 maxnonstiffstep = maxnonstiffstep, nonstifftol = nonstifftol, stifftol = stifftol, dtfac = dtfac,
                 stiffalgfirst = stiffalgfirst, switch_max = switch_max), class = "SDEAutoSwitch")

# An AutoSwitch() object of the ODE solvers carries its settings in the
# closure of its choice function.
.sde_spec_from_autoswitch <- function(x) {
  env <- environment(x$choose)
  .sde_autoswitch_spec(x$algs[[1L]], x$algs[[2L]], env$maxstiffstep, env$maxnonstiffstep, env$nonstifftol,
                       env$stifftol, env$dtfac, env$stiffalgfirst, env$switch_max %||% 5L)
}

#' @rdname sde_composite
#' @export
StochasticCompositeAlgorithm <- function(algs, choice_function) {
  if (!is.list(algs) || inherits(algs, "SDEAlgorithm")) algs <- list(algs)
  for (a in algs) if (!inherits(a, "SDEAlgorithm")) stop("not an SDE algorithm: ", class(a)[1L], call. = FALSE)
  spec <- if (inherits(choice_function, "SDEAutoSwitch")) choice_function else
    if (inherits(choice_function, "CompositeAlgorithm")) .sde_spec_from_autoswitch(choice_function) else NULL
  if (is.null(spec) && !is.function(choice_function))
    stop("choice_function is an AutoSwitch() object or a function", call. = FALSE)
  zl <- Filter(Negate(is.null), lapply(algs, `[[`, "z_length"))
  structure(list(name = "StochasticCompositeAlgorithm", algs = algs, choice_function = choice_function,
                 spec = spec, order = max(vapply(algs, function(a) .as_num(a$order), 0)),
                 adaptive = all(vapply(algs, function(a) isTRUE(a$adaptive), TRUE)),
                 needs_dZ = any(vapply(algs, function(a) isTRUE(a$needs_dZ), TRUE)),
                 noise = unique(unlist(lapply(algs, `[[`, "noise"))), delta = 1,
                 interpretation = algs[[1L]]$interpretation,
                 z_length = if (length(zl)) function(m, dt) max(vapply(zl, function(z) z(m, dt), 0))),
            class = c("StochasticCompositeAlgorithm", "CompositeAlgorithm", "SDEAlgorithm"))
}

#' @export
print.StochasticCompositeAlgorithm <- function(x, ...) {
  cat("StochasticCompositeAlgorithm(", paste(vapply(x$algs, `[[`, "", "name"), collapse = ", "), ")\n", sep = "")
  invisible(x)
}

#' @rdname sde_composite
#' @export
AutoAlgSwitch <- function(nonstiffalg, stiffalg, maxstiffstep = 10L, maxnonstiffstep = 3L, nonstifftol = 9 / 10,
                          stifftol = 9 / 10, dtfac = 2, stiffalgfirst = FALSE, switch_max = 5L) {
  spec <- .sde_autoswitch_spec(nonstiffalg, stiffalg, maxstiffstep, maxnonstiffstep, nonstifftol, stifftol,
                               dtfac, stiffalgfirst, switch_max)
  StochasticCompositeAlgorithm(list(nonstiffalg, stiffalg), spec)
}

#' @rdname sde_composite
#' @export
AutoSOSRI2 <- function(alg, ...) AutoAlgSwitch(SOSRI2(), alg, ...)

#' @rdname sde_composite
#' @export
AutoSOSRA2 <- function(alg, ...) AutoAlgSwitch(SOSRA2(), alg, ...)

# (AS::AutoSwitchCache)(integrator) with is_stiff: the method index and the
# step size after the choice.
.sde_autoswitch_choose <- function(spec, AS, iter, eigen_est, dt) {
  if (iter == 0L) return(list(current = if (spec$stiffalgfirst) 2L else 1L, dt = dt))
  stiffness <- abs(eigen_est * dt / (spec$nonstiffalg$stability_size %||% 0))
  tol <- if (AS$is_stiffalg) spec$stifftol else spec$nonstifftol
  stiff <- isTRUE(stiffness > tol)
  AS$successive_switches <- if (stiff) 0L else AS$successive_switches + 1L
  AS$do_error_check <- AS$successive_switches > spec$switch_max || !stiff
  AS$count <- if (stiff) { if (AS$count < 0L) 1L else AS$count + 1L } else
    { if (AS$count > 0L) -1L else AS$count - 1L }
  if (!AS$is_stiffalg && AS$count > spec$maxstiffstep) {
    dt <- dt * spec$dtfac
    AS$is_stiffalg <- TRUE
  } else if (AS$is_stiffalg && AS$count < -spec$maxnonstiffstep) {
    dt <- dt / spec$dtfac
    AS$is_stiffalg <- FALSE
  }
  list(current = if (AS$is_stiffalg) 2L else 1L, dt = dt)
}

# The core driver method of a composite; run_step(a, ctx, t, u, dt, cache)
# takes the noise increments and steps with member a (see .sde_setup).
.sde_composite_runner <- function(alg, adaptive, run_step, mm) {
  algs <- alg$algs
  first <- algs[[1L]]$name
  spec <- alg$spec
  st <- new.env(parent = emptyenv())
  st$started <- FALSE
  member_init <- function(i, ctx, t, u) (algs[[i]]$init %||% function(ctx, t, u) list())(ctx, t, u)
  choose <- function(iter, eigen_est, dt, integ) {
    if (!is.null(spec)) return(.sde_autoswitch_choose(spec, st$AS, iter, eigen_est, dt))
    i <- as.integer(alg$choice_function(list(integrator = integ, t = integ$t, u = integ$u, p = integ$p,
                                             iter = iter, dt = dt, eigen_est = eigen_est)))
    if (length(i) != 1L || is.na(i) || i < 1L || i > length(algs))
      stop("the choice function must return the index of a method", call. = FALSE)
    list(current = i, dt = dt)
  }
  init <- function(ctx, t, u) {
    if (!st$started) {
      st$started <- TRUE
      st$ctx <- ctx
      st$AS <- new.env(parent = emptyenv())
      st$AS$count <- 0L
      st$AS$successive_switches <- 0L
      st$AS$is_stiffalg <- isTRUE(spec$stiffalgfirst)
      st$AS$do_error_check <- TRUE
      st$caches <- vector("list", length(algs))
      st$current <- choose(0L, 1, NA_real_, list(t = t, u = u, p = ctx$p))$current
      st$first <- st$current
      st$log <- list()
    }
    # initialize! of the current member (also after the state was modified)
    st$caches[[st$current]] <- member_init(st$current, ctx, t, u)
    list()
  }
  step <- function(ctx, t, u, dt, cache) {
    i <- st$current
    res <- run_step(algs[[i]], ctx, t, u, dt, st$caches[[i]])
    if (!is.null(res$cache)) st$caches[[i]] <- res$cache
    st$log[[length(st$log) + 1L]] <- c(t, t + dt, i)
    res$eigen_est <- if (first == "SOSRI2") res$eigen_sri else if (first == "SOSRA2") res$eigen_sra
    res$cache <- cache
    res
  }
  runner <- .new_alg(alg$name, alg$order, adaptive = adaptive, step = step, interp = .sde_linear_interp,
                     init = init, mass_matrix = mm,
                     on_reject = function(cache) {
                       i <- st$current
                       if (!is.null(algs[[i]]$step)) .nl_on_reject(st$caches[[i]])
                     })
  # choose_algorithm! in loopheader!: dt_seen is the step size of the
  # integrator, dt the controller's proposal the driver keeps otherwise
  runner$pre_step <- function(cache, iter, dt_seen, eigen_est, dt, integ) {
    ch <- choose(iter, eigen_est, dt_seen, integ)
    if (ch$current == st$current) return(if (ch$dt == dt_seen) dt else ch$dt)
    new <- ch$current
    if (is.null(st$caches[[new]]) || isTRUE(algs[[new]]$init_on_switch))
      st$caches[[new]] <- member_init(new, st$ctx, integ$t, integ$u)
    st$current <- new
    ch$dt
  }
  runner$alg_choice <- function(ts) {
    if (!length(st$log)) return(rep(st$first, length(ts)))
    lg <- do.call(rbind, st$log)
    vapply(seq_along(ts), function(k) {
      if (k == 1L) return(as.integer(st$first))
      tol <- 1e-12 * max(1, abs(ts[k]))
      idx <- which(lg[, 1L] < ts[k] - tol & lg[, 2L] >= ts[k] - tol)
      as.integer(if (length(idx)) lg[max(idx), 3L] else lg[nrow(lg), 3L])
    }, 0L)
  }
  runner
}
