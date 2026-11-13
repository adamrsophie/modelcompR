# The Rosenbrock and Rosenbrock-W methods of OrdinaryDiffEqRosenbrock
# without a step limiter (algorithms.jl, the second @eval loop, and
# RosenbrockW6S4OS): RodasTableau methods run by the unified loop of
# R/solvers-rosenbrock.R. None has an H matrix, so the dense output is the
# Hermite interpolant of f at both ends of the step; all are FSAL in
# OrdinaryDiffEqCore's sense (reset_fsal! after callbacks).

#' Classical Rosenbrock and Rosenbrock-W algorithms
#'
#' Rosenbrock methods of OrdinaryDiffEq.jl with an embedded error estimate
#' and Hermite dense output: `ROS2()` (order 2, L-stable), `ROS2PR()` and the
#' W-method `ROS2S()` (order 2, stiffly accurate), `ROS3()`, `ROS3PR()`,
#' `Scholz4_7()`, `ROS3PRL()` and `ROS3PRL2()` (order 3, Rang), the
#' Rosenbrock-W methods `ROS34PW1a()`, `ROS34PW1b()`, `ROS34PW2()`,
#' `ROS34PW3()` and `ROS34PRw()` of Rang and Angermann, the Rosenbrock-Krylov
#' method `ROK4a()`, the fourth order methods `RosShamp4()` (Shampine),
#' `Veldd4()`, `Velds4()` (van Veldhuizen), `GRK4T()`, `GRK4A()` (Kaps and
#' Rentrop) and `Ros4LStab()` (Hairer and Wanner), and the fixed step fourth
#' order W-method `RosenbrockW6S4OS()`, which needs `dt`. All accept mass
#' matrices; the W-methods keep the Jacobian along accepted steps (see
#' [Rosenbrock23()]).
#'
#' @inheritParams rosenbrock_algorithms
#' @return An algorithm object.
#' @examples
#' prob <- ODEProblem(function(u, p, t) c(-50 * (u[1] - u[2]), u[1] - u[2]), c(1, 0), c(0, 1))
#' solve(prob, ROS34PW2())$stats$njacs
#' solve(prob, GRK4A())(1)
#' solve(prob, RosenbrockW6S4OS(), dt = 0.1)(1)
#' @name rosenbrock_w_algorithms
NULL

#' @rdname rosenbrock_w_algorithms
#' @export
ROS2 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                 jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS2", 2L, "ROS2RodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
ROS2PR <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                   jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS2PR", 2L, "ROS2PRRodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
ROS2S <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                  jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS2S", 2L, "ROS2SRodasTableau", wmethod = TRUE)

#' @rdname rosenbrock_w_algorithms
#' @export
ROS3 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                 jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS3", 3L, "ROS3RodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
ROS3PR <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                   jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS3PR", 3L, "ROS3PRRodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
Scholz4_7 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                      jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Scholz4_7", 3L, "Scholz4_7RodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
ROS34PW1a <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                      jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS34PW1a", 3L, "ROS34PW1aRodasTableau", wmethod = TRUE)

#' @rdname rosenbrock_w_algorithms
#' @export
ROS34PW1b <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                      jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS34PW1b", 3L, "ROS34PW1bRodasTableau", wmethod = TRUE)

#' @rdname rosenbrock_w_algorithms
#' @export
ROS34PW2 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                     jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS34PW2", 3L, "ROS34PW2RodasTableau", wmethod = TRUE)

#' @rdname rosenbrock_w_algorithms
#' @export
ROS34PW3 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                     jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS34PW3", 4L, "ROS34PW3RodasTableau", wmethod = TRUE)

#' @rdname rosenbrock_w_algorithms
#' @export
ROS34PRw <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                     jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS34PRw", 3L, "ROS34PRwRodasTableau", wmethod = TRUE)

#' @rdname rosenbrock_w_algorithms
#' @export
ROS3PRL <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                    jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS3PRL", 3L, "ROS3PRLRodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
ROS3PRL2 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                     jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROS3PRL2", 3L, "ROS3PRL2RodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
ROK4a <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                  jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "ROK4a", 4L, "ROK4aRodasTableau", wmethod = TRUE)

#' @rdname rosenbrock_w_algorithms
#' @export
RosShamp4 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                      jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "RosShamp4", 4L, "RosShamp4RodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
Veldd4 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                   jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Veldd4", 4L, "Veldd4RodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
Velds4 <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                   jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Velds4", 4L, "Velds4RodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
GRK4T <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                  jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "GRK4T", 4L, "GRK4TRodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
GRK4A <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                  jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "GRK4A", 4L, "GRK4ARodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
Ros4LStab <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                      jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "Ros4LStab", 4L, "Ros4LStabRodasTableau")

#' @rdname rosenbrock_w_algorithms
#' @export
RosenbrockW6S4OS <- function(autodiff = NULL, concrete_jac = NULL, linsolve = NULL, max_jac_age = 20L,
                             jac_reuse_gamma_tol = 0.03)
  .rosen_from(environment(), "RosenbrockW6S4OS", 4L, "RosenbrockW6S4OSRodasTableau", wmethod = TRUE, adaptive = FALSE)
