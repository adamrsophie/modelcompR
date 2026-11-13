# Automatic stiffness switching and the default ODE algorithm, after
# OrdinaryDiffEq.jl:
#   lib/OrdinaryDiffEqCore/src/composite_algs.jl  AutoSwitch, is_stiff,
#     AutoAlgSwitch (maxstiffstep 10, maxnonstiffstep 3, tolerances 9/10,
#     dtfac 2)
#   lib/OrdinaryDiffEqDefault/src/default_alg.jl  DefaultODEAlgorithm,
#     default_autoswitch, nonstiffchoice, stiffchoice
#   lib/OrdinaryDiffEqTsit5/src/alg_utils.jl, OrdinaryDiffEqVerner/src/
#     alg_utils.jl  alg_stability_size (Tsit5 3.5068, Vern7 4.64)
#
# A composite algorithm has class c("CompositeAlgorithm", "ODEAlgorithm")
# with `algs`, `choice_state()` (a fresh AutoSwitchCache environment) and
# `choose(AS, info)`, which returns the index of the method for the next
# step (and on later calls list(current, dt)). The driver calls it before
# every step with the latest stiffness estimate eigen_est and dt.

.composite_alg <- function(name, algs, choice_state, choose) {
  structure(list(name = name, algs = algs, order = algs[[1L]]$order, adaptive = TRUE,
                 mass_matrix = all(vapply(algs, function(a) isTRUE(a$mass_matrix), TRUE)) ||
                   name %in% c("DefaultODEAlgorithm", "DefaultImplicitODEAlgorithm"),
                 stiff = FALSE, choice_state = choice_state, choose = choose),
            class = c(name, "CompositeAlgorithm", "ODEAlgorithm"))
}

.autoswitch_state <- function(stiffalgfirst = FALSE) {
  AS <- new.env(parent = emptyenv())
  AS$count <- 0L
  AS$is_stiffalg <- stiffalgfirst
  AS$current <- 0L
  AS
}

# is_stiff: the stiffness ratio |eigen_est * dt / stability size| against
# tol; NaN counts as stiff, and Inf is not stiff against tol = Inf
.is_stiff <- function(eigen_est, dt, size, tol) {
  !isTRUE(abs(eigen_est * dt / size) <= tol)
}

.count_update <- function(AS, stiff) {
  AS$count <- if (stiff) { if (AS$count < 0L) 1L else AS$count + 1L } else
    { if (AS$count > 0L) -1L else AS$count - 1L }
}

#' Automatic switching between a nonstiff and a stiff method
#'
#' `AutoSwitch(nonstiffalg, stiffalg)` builds a composite algorithm that
#' starts with the nonstiff method and switches to the stiff one after
#' `maxstiffstep` consecutive steps flagged as stiff (the stiffness estimate
#' times the step size exceeds the nonstiff method's stability region), and
#' back after `maxnonstiffstep` steps flagged as nonstiff, as
#' `AutoAlgSwitch` in OrdinaryDiffEq.jl. `AutoTsit5(alg)` and `AutoVern7(alg)`
#' pair `Tsit5()` or `Vern7()` with a stiff method.
#'
#' @param nonstiffalg,stiffalg Algorithms.
#' @param alg A stiff algorithm such as [Rosenbrock23()].
#' @param maxstiffstep,maxnonstiffstep Consecutive steps before switching.
#' @param nonstifftol,stifftol Stiffness thresholds.
#' @param dtfac Step size factor applied when switching.
#' @param stiffalgfirst Start with the stiff method.
#' @return A composite algorithm.
#' @examples
#' f <- function(u, p, t) c(-0.04 * u[1] + 1e4 * u[2] * u[3],
#'                          0.04 * u[1] - 1e4 * u[2] * u[3] - 3e7 * u[2]^2,
#'                          3e7 * u[2]^2)
#' prob <- ODEProblem(f, c(1, 0, 0), c(0, 1))
#' sol <- solve(prob, AutoTsit5(Rosenbrock23()))
#' table(sol$alg_choice)
#' @export
AutoSwitch <- function(nonstiffalg, stiffalg, maxstiffstep = 10L, maxnonstiffstep = 3L,
                       nonstifftol = 9 / 10, stifftol = 9 / 10, dtfac = 2, stiffalgfirst = FALSE) {
  size <- nonstiffalg$stability_size %||% 1
  choose <- function(AS, info) {
    if (isTRUE(info$init)) {
      AS$current <- if (stiffalgfirst) 2L else 1L
      return(AS$current)
    }
    tol <- if (AS$is_stiffalg) stifftol else nonstifftol
    .count_update(AS, .is_stiff(info$eigen_est, info$dt, size, tol))
    dt <- info$dt
    if (!AS$is_stiffalg && AS$count > maxstiffstep) {
      dt <- dt * dtfac
      AS$is_stiffalg <- TRUE
    } else if (AS$is_stiffalg && AS$count < -maxnonstiffstep) {
      dt <- dt / dtfac
      AS$is_stiffalg <- FALSE
    }
    AS$current <- if (AS$is_stiffalg) 2L else 1L
    list(current = AS$current, dt = dt)
  }
  .composite_alg(paste0("AutoSwitch(", nonstiffalg$name, ", ", stiffalg$name, ")"),
                 list(nonstiffalg, stiffalg), function() .autoswitch_state(stiffalgfirst), choose)
}

#' @rdname AutoSwitch
#' @param ... Options for [AutoSwitch()].
#' @export
AutoTsit5 <- function(alg, ...) AutoSwitch(Tsit5(), alg, ...)

#' @rdname AutoSwitch
#' @export
AutoVern7 <- function(alg, ...) AutoSwitch(Vern7(), alg, ...)

#' @rdname AutoSwitch
#' @export
AutoVern6 <- function(alg, ...) AutoSwitch(Vern6(), alg, ...)

#' @rdname AutoSwitch
#' @export
AutoVern8 <- function(alg, ...) AutoSwitch(Vern8(), alg, ...)

#' @rdname AutoSwitch
#' @export
AutoVern9 <- function(alg, ...) AutoSwitch(Vern9(), alg, ...)

# DefaultSolverChoice indices: Tsit5 = 1, Vern7 = 2, Rosenbrock23 = 3,
# Rodas5P = 4, FBDF = 5 (KrylovFBDF, choice 6, is FBDF with a dense solver
# here).
.DEFAULT_LOW_TOL <- 1e-6
.DEFAULT_SMALLSIZE <- 50L
.DEFAULT_MEDIUMSIZE <- 500L

.nonstiffchoice <- function(reltol) if (reltol < .DEFAULT_LOW_TOL) 2L else 1L

.stiffchoice <- function(reltol, len, mass) {
  if (len > .DEFAULT_SMALLSIZE) return(5L)
  if (reltol < .DEFAULT_LOW_TOL || !.is_identity_mass(mass)) 4L else 3L
}

#' The default ODE algorithm
#'
#' `DefaultODEAlgorithm()` is what [solve()] uses when no algorithm is
#' given, as in OrdinaryDiffEqDefault. With an identity mass matrix it
#' starts with `Tsit5()` (`Vern7()` when `reltol < 1e-6`) and switches to a
#' stiff method when stiffness is detected: `Rosenbrock23()` (`Rodas5P()`
#' when `reltol < 1e-6`) up to 50 states and `FBDF()` above, switching
#' back when the problem becomes nonstiff. With a non identity mass matrix
#' it uses `Rodas5P()` up to 50 states and `FBDF()` above, without
#' switching. `sol$alg_choice` records the method of every step (1 Tsit5,
#' 2 Vern7, 3 Rosenbrock23, 4 Rodas5P, 5 FBDF).
#'
#' `DefaultImplicitODEAlgorithm()` starts with the stiff method and
#' compares the stiffness estimate with `stol` while stiff and `ntol` while
#' nonstiff; with the defaults 0 and `Inf` it keeps the stiff method.
#'
#' @param stiffalgfirst Start with the stiff method.
#' @return A composite algorithm.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -u, 1, c(0, 1))
#' sol <- solve(prob)
#' unique(sol$alg_choice)
#' sol <- solve(ODEProblem(function(u, p, t) -1000 * u, 1, c(0, 1)), DefaultImplicitODEAlgorithm())
#' unique(sol$alg_choice)
#' @export
DefaultODEAlgorithm <- function(stiffalgfirst = FALSE) .default_alg("DefaultODEAlgorithm", stiffalgfirst, 9 / 10, 9 / 10)

#' @rdname DefaultODEAlgorithm
#' @param stol,ntol Stiffness thresholds while the stiff and the nonstiff
#'   method runs.
#' @export
DefaultImplicitODEAlgorithm <- function(stol = 0, ntol = Inf) .default_alg("DefaultImplicitODEAlgorithm", TRUE, stol, ntol)

.default_alg <- function(name, stiffalgfirst, stol, ntol) {
  algs <- list(Tsit5(), Vern7(), Rosenbrock23(), Rodas5P(), FBDF())
  sizes <- c(algs[[1L]]$stability_size, algs[[2L]]$stability_size)
  maxstiffstep <- 10L; maxnonstiffstep <- 3L; dtfac <- 2
  choose <- function(AS, info) {
    identity_mass <- .is_identity_mass(info$mass)
    if (isTRUE(info$init)) {
      AS$current <- if (stiffalgfirst || !identity_mass) .stiffchoice(info$reltol, info$n, info$mass) else
        .nonstiffchoice(info$reltol)
      AS$is_stiffalg <- AS$current > 2L
      return(AS$current)
    }
    dt <- info$dt
    stiff <- .is_stiff(info$eigen_est, dt, sizes[.nonstiffchoice(info$reltol)], if (AS$is_stiffalg) stol else ntol)
    .count_update(AS, stiff)
    if (!identity_mass) {
      # a mass matrix problem keeps its stiff method
    } else if (!AS$is_stiffalg && AS$count > maxstiffstep) {
      dt <- dt * dtfac
      AS$is_stiffalg <- TRUE
      AS$current <- .stiffchoice(info$reltol, info$n, info$mass)
    } else if (AS$is_stiffalg && AS$count < -maxnonstiffstep) {
      dt <- dt / dtfac
      AS$is_stiffalg <- FALSE
      AS$current <- .nonstiffchoice(info$reltol)
    }
    list(current = AS$current, dt = dt)
  }
  .composite_alg(name, algs, function() .autoswitch_state(stiffalgfirst), choose)
}

#' @export
print.CompositeAlgorithm <- function(x, ...) {
  cat(x$name, if (inherits(x, c("DefaultODEAlgorithm", "DefaultImplicitODEAlgorithm"))) "()", "\n", sep = "")
  invisible(x)
}
