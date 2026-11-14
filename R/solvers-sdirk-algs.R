# Constructors of the SDIRK, ESDIRK and IMEX methods of OrdinaryDiffEqSDIRK
# (algorithms.jl); the step is in R/solvers-sdirk.R.

#' SDIRK, ESDIRK and IMEX Runge-Kutta methods
#'
#' Singly diagonally implicit Runge-Kutta methods of OrdinaryDiffEqSDIRK for
#' stiff problems, mass matrix DAEs of index 1 and split problems. Each
#' implicit stage is solved by the modified Newton iteration of
#' [ImplicitEuler()], which reuses the Jacobian across stages and steps; the
#' Jacobian is `prob$f$jac` when the problem has one and finite differences
#' otherwise.
#'
#' Adaptive methods with an embedded error estimate:
#' * `TRBDF2()`: the second order A, B, L and S stable TR-BDF2 method, with
#'   a smoothed error estimate; good for loose tolerances.
#' * `SDIRK2()`: a second order A, B and L stable SDIRK method.
#' * `Kvaerno3()`, `Kvaerno4()`, `Kvaerno5()`: A and L stable, stiffly
#'   accurate ESDIRK methods of orders 3, 4 and 5.
#' * `KenCarp3()`, `KenCarp4()`, `KenCarp47()`, `KenCarp5()`,
#'   `KenCarp58()`: A and L stable, stiffly accurate ESDIRK methods of orders
#'   3 to 5 (`KenCarp47()` with seven, `KenCarp58()` with eight stages) that
#'   are additive (IMEX) methods on a [SplitODEProblem()].
#' * `Cash4()`, `Hairer4()`, `Hairer42()`: A and L stable fourth order SDIRK
#'   methods. `Cash4(embedding = 3)` uses the third order embedded method of
#'   Cash, any other value the second embedding.
#' * `ESDIRK325L2SA()`, `ESDIRK436L2SA2()`, `ESDIRK437L2SA()`,
#'   `ESDIRK547L2SA2()`, `ESDIRK54I8L2SA()`, `ESDIRK659L2SA()`: the optimized
#'   ESDIRK methods of Kennedy and Carpenter, stiffly accurate with an L
#'   stable embedding, of orders 3 to 6 (the first digit).
#'
#' Fixed step methods (give `dt`):
#' * `ImplicitMidpoint()`: the second order symmetric and symplectic
#'   implicit midpoint rule.
#' * `SSPSDIRK2()`: a second order strong stability preserving SDIRK method
#'   (SSP coefficient 4).
#' * `SFSDIRK4()` to `SFSDIRK8()`: fourth order SDIRK methods with four to
#'   eight stages.
#' * `CFNLIRK3()`: a third order IMEX method.
#' * `ARS222()`, `ARS232()`, `ARS343()`, `ARS443()`: the IMEX methods of
#'   Ascher, Ruuth and Spiteri (ARS(s, sigma, p)).
#' * `IMEXSSP222()`, `IMEXSSP2322()`, `IMEXSSP3332()`, `IMEXSSP3433()`: the
#'   IMEX strong stability preserving methods of Pareschi and Russo, with an
#'   implicit first stage.
#' * `BHR553()`: the five stage third order IMEX method BHR(5,5,3)* of
#'   Boscarino and Russo.
#'
#' `SDIRK22()` is the [Trapezoid()] method, as in OrdinaryDiffEqSDIRK, where
#' the constructor returns a `Trapezoid` (whose step does not use the stage
#' predictor).
#'
#' The IMEX methods (`KenCarp*`, `ARS*`, `IMEXSSP*`, `BHR553()` and
#' `CFNLIRK3()`) treat `f1` of a [SplitODEProblem()] implicitly, with the
#' Jacobian given to the split function, and `f2` explicitly; on other
#' problems they integrate `f` implicitly. The methods with an explicit
#' first stage (the ESDIRK, `TRBDF2()`, `ARS*`, `BHR553()` and `CFNLIRK3()`)
#' accept only diagonal mass matrices. Dense output is the cubic Hermite
#' interpolant.
#'
#' @param autodiff How Jacobians are computed when the problem has no
#'   symbolic Jacobian; see [ImplicitEuler()].
#' @param smooth_est Multiply the embedded error estimate by the inverse of
#'   the Newton matrix (a smoothed estimate for very stiff problems).
#' @param predictor The initial guess of the Newton iteration of each stage,
#'   a value of [Predictor]: `Trivial` (zero), `Linear` (the linear
#'   extrapolation of the derivative for an implicit first stage, the
#'   tableau's stage guesses afterwards), `Tableau` (the tableau's stage
#'   guesses), `CopyPrev` (the previous stage), `StageExtrap` (the linear
#'   extrapolation of the two previous stages), or `MaxOrder`,
#'   `VariableOrder` and `CutoffOrder` (the dense output of the last step,
#'   of order 3, of an order falling with the distance of the extrapolation,
#'   or of order 3 only for stages in the first half of the step).
#' @param extrapolant Deprecated: `"constant"`, `"linear"` or
#'   `"interpolant"`, the older names of `Predictor$Trivial`,
#'   `Predictor$Linear` and `Predictor$MaxOrder`.
#' @param step_limiter A function `step_limiter(u, integ, p, t)` returning
#'   the state after each accepted step, or `NULL`.
#' @param embedding The embedded method of `Cash4()`, 3 or 2.
#' @param nlsolve The nonlinear solver of the implicit stages: [NLNewton()]
#'   (with its `kappa`, `max_iter`, `new_W_dt_cutoff`, `always_new`,
#'   `check_div`, `relax` and `predictor`; `fast_convergence_cutoff` is
#'   stored, as OrdinaryDiffEq stores it without using it), the fixed point
#'   iteration [NLFunctional()] or its Anderson accelerated form
#'   [NLAnderson()]. A `predictor(uprev, p, t, dt)` returns a guess of the
#'   stage value at the stage time `t`; it replaces the guess of `predictor`
#'   at every stage.
#' @return An algorithm object.
#' @examples
#' prob <- ODEProblem(function(u, p, t) -50 * (u - cos(t)), 0, c(0, 1))
#' sol <- solve(prob, KenCarp4())
#' sol(1)
#' solve(prob, TRBDF2(), reltol = 1e-6)$stats$naccept
#' solve(prob, ImplicitMidpoint(), dt = 0.05)$u[[21]]
#' sp <- SplitODEProblem(function(u, p, t) -50 * u, function(u, p, t) 50 * cos(t), 0, c(0, 1))
#' solve(sp, ARS343(), dt = 0.05)$u[[21]]
#' @name sdirk_algorithms
NULL

#' @rdname sdirk_algorithms
#' @export
ImplicitMidpoint <- function(autodiff = NULL, predictor = Predictor$Linear, extrapolant = NULL, step_limiter = NULL,
                             nlsolve = NLNewton())
  .sdirk_alg("ImplicitMidpoint", autodiff, FALSE, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
TRBDF2 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                   step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("TRBDF2", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
SDIRK2 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                   step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("SDIRK2", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
SDIRK22 <- function(autodiff = NULL, predictor = Predictor$Linear, extrapolant = NULL, step_limiter = NULL,
                    nlsolve = NLNewton()) {
  .sdirk_predictor(predictor, extrapolant)
  if (!is.null(step_limiter) && !is.function(step_limiter))
    stop("step_limiter must be a function(u, integ, p, t)", call. = FALSE)
  alg <- Trapezoid(autodiff = autodiff, nlsolve = nlsolve)
  alg$step_limiter <- step_limiter
  alg
}

#' @rdname sdirk_algorithms
#' @export
SSPSDIRK2 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Trivial, extrapolant = NULL,
                      nlsolve = NLNewton())
  .sdirk_alg("SSPSDIRK2", autodiff, FALSE, predictor, extrapolant, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
Kvaerno3 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                     step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("Kvaerno3", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
KenCarp3 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                     step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("KenCarp3", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
CFNLIRK3 <- function(autodiff = NULL, predictor = Predictor$Linear, extrapolant = NULL, nlsolve = NLNewton())
  .sdirk_alg("CFNLIRK3", autodiff, FALSE, predictor, extrapolant, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
Cash4 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                  embedding = 3L, nlsolve = NLNewton())
  .sdirk_alg("Cash4", autodiff, smooth_est, predictor, extrapolant,
             tab = .SDIRK_TABS[[if (embedding == 3) "Cash4" else "Cash4_e2"]], nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
SFSDIRK4 <- function(autodiff = NULL, predictor = Predictor$Linear, extrapolant = NULL, nlsolve = NLNewton())
  .sdirk_alg("SFSDIRK4", autodiff, FALSE, predictor, extrapolant, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
SFSDIRK5 <- function(autodiff = NULL, predictor = Predictor$Linear, extrapolant = NULL, nlsolve = NLNewton())
  .sdirk_alg("SFSDIRK5", autodiff, FALSE, predictor, extrapolant, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
SFSDIRK6 <- function(autodiff = NULL, predictor = Predictor$Linear, extrapolant = NULL, nlsolve = NLNewton())
  .sdirk_alg("SFSDIRK6", autodiff, FALSE, predictor, extrapolant, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
SFSDIRK7 <- function(autodiff = NULL, predictor = Predictor$Linear, extrapolant = NULL, nlsolve = NLNewton())
  .sdirk_alg("SFSDIRK7", autodiff, FALSE, predictor, extrapolant, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
SFSDIRK8 <- function(autodiff = NULL, predictor = Predictor$Linear, extrapolant = NULL, nlsolve = NLNewton())
  .sdirk_alg("SFSDIRK8", autodiff, FALSE, predictor, extrapolant, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
Hairer4 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                    nlsolve = NLNewton())
  .sdirk_alg("Hairer4", autodiff, smooth_est, predictor, extrapolant, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
Hairer42 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                     nlsolve = NLNewton())
  .sdirk_alg("Hairer42", autodiff, smooth_est, predictor, extrapolant, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
Kvaerno4 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                     step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("Kvaerno4", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
Kvaerno5 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                     step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("Kvaerno5", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
KenCarp4 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                     step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("KenCarp4", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
KenCarp47 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                      step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("KenCarp47", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
KenCarp5 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                     step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("KenCarp5", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
KenCarp58 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                      step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("KenCarp58", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
ESDIRK54I8L2SA <- function(autodiff = NULL, smooth_est = FALSE, predictor = Predictor$StageExtrap,
                           extrapolant = NULL, step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("ESDIRK54I8L2SA", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
ESDIRK436L2SA2 <- function(autodiff = NULL, smooth_est = FALSE, predictor = Predictor$StageExtrap,
                           extrapolant = NULL, step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("ESDIRK436L2SA2", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
ESDIRK437L2SA <- function(autodiff = NULL, smooth_est = FALSE, predictor = Predictor$StageExtrap,
                          extrapolant = NULL, step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("ESDIRK437L2SA", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
ESDIRK547L2SA2 <- function(autodiff = NULL, smooth_est = FALSE, predictor = Predictor$StageExtrap,
                           extrapolant = NULL, step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("ESDIRK547L2SA2", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
ESDIRK659L2SA <- function(autodiff = NULL, smooth_est = FALSE, predictor = Predictor$StageExtrap,
                          extrapolant = NULL, step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("ESDIRK659L2SA", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
ESDIRK325L2SA <- function(autodiff = NULL, smooth_est = FALSE, predictor = Predictor$StageExtrap,
                          extrapolant = NULL, step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("ESDIRK325L2SA", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
ARS222 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                   step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("ARS222", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
ARS232 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                   step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("ARS232", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
ARS343 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                   step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("ARS343", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
ARS443 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                   step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("ARS443", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
IMEXSSP222 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                       step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("IMEXSSP222", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
IMEXSSP2322 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                        step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("IMEXSSP2322", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
IMEXSSP3332 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                        step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("IMEXSSP3332", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
IMEXSSP3433 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                        step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("IMEXSSP3433", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)

#' @rdname sdirk_algorithms
#' @export
BHR553 <- function(autodiff = NULL, smooth_est = TRUE, predictor = Predictor$Linear, extrapolant = NULL,
                   step_limiter = NULL, nlsolve = NLNewton())
  .sdirk_alg("BHR553", autodiff, smooth_est, predictor, extrapolant, step_limiter, nlsolve = nlsolve)
