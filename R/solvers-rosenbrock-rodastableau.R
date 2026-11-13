# The exported tableau constructors of OrdinaryDiffEqRosenbrockTableaus 2.4.2
# (rosenbrock_tableaus.jl): RodasTableau and one constructor per method,
# returning the Float64 tableau of R/solvers-rosenbrock-tableaus.R as a list
# with Julia's field names.

#' Rosenbrock tableaus
#'
#' The coefficient tableaus of the Rosenbrock methods, as exported by
#' OrdinaryDiffEqRosenbrockTableaus.jl. `RodasTableau()` builds a tableau
#' from its fields; the other functions return the tableau of one method
#' (for example `Rodas4Tableau()` for [Rodas4()]) as the Julia constructor
#' builds it for `Float64` coefficients. A tableau is a list of class
#' `"RodasTableau"` with the stage matrices `A` and `C`, the diagonal
#' coefficient `gamma`, the abscissae `c`, the time derivative coefficients
#' `d`, the dense output matrix `H` (no rows for methods without a stiff
#' aware interpolant), the solution weights `b` and the embedded error
#' weights `btilde` (`NULL` for fixed step methods).
#'
#' @param A,C Stage coefficient matrices (lower triangular).
#' @param gamma Diagonal coefficient.
#' @param c Abscissae.
#' @param d Time derivative coefficients.
#' @param H Dense output matrix, one row per interpolation coefficient.
#' @param b Solution weights.
#' @param btilde Embedded error estimate weights, or `NULL`.
#' @param T,T2 The element types of the coefficients and of `gamma` and `c`;
#'   only `"Float64"` is available.
#' @return A list of class `"RodasTableau"`.
#' @examples
#' tab <- Rodas4Tableau()
#' tab$gamma
#' dim(tab$A)
#' is.null(RosenbrockW6S4OSRodasTableau()$btilde)
#' @name rosenbrock_tableaus
NULL

.rosen_check_type <- function(x, arg) {
  if (!(is.character(x) && length(x) == 1L && x %in% c("Float64", "double", "numeric")))
    stop("`", arg, "` must be \"Float64\": the tableaus are available in double precision only", call. = FALSE)
  invisible(x)
}

.rosen_tableau <- function(name, T, T2) {
  .rosen_check_type(T, "T")
  .rosen_check_type(T2, "T2")
  structure(.ROSEN_TABLEAUS[[name]], class = "RodasTableau")
}

#' @rdname rosenbrock_tableaus
#' @export
RodasTableau <- function(A, C, gamma, c, d, H, b, btilde)
  structure(list(A = as.matrix(A), C = as.matrix(C), gamma = as.numeric(gamma), c = as.numeric(c),
                 d = as.numeric(d), H = as.matrix(H), b = as.numeric(b),
                 btilde = if (is.null(btilde)) NULL else as.numeric(btilde)),
            class = "RodasTableau")

#' @rdname rosenbrock_tableaus
#' @export
Rodas5Tableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Rodas5Tableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Rodas4Tableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Rodas4Tableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Rodas42Tableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Rodas42Tableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Rodas4PTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Rodas4PTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Rodas4P2Tableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Rodas4P2Tableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Rodas4PWTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Rodas4PWTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS3PRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS3PRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Rodas3RodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Rodas3RodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Rodas3dRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Rodas3dRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Rodas3PRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Rodas3PRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
RosShamp4RodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("RosShamp4RodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Veldd4RodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Veldd4RodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Velds4RodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Velds4RodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
GRK4TRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("GRK4TRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
GRK4ARodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("GRK4ARodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Ros4LStabRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Ros4LStabRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS2RodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS2RodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS2PRRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS2PRRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS2SRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS2SRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS3RodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS3RodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS3PRRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS3PRRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
Scholz4_7RodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("Scholz4_7RodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS34PW1aRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS34PW1aRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS34PW1bRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS34PW1bRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS34PW2RodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS34PW2RodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS34PRwRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS34PRwRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS3PRLRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS3PRLRodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
ROS3PRL2RodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("ROS3PRL2RodasTableau", T, T2)

#' @rdname rosenbrock_tableaus
#' @export
RosenbrockW6S4OSRodasTableau <- function(T = "Float64", T2 = "Float64") .rosen_tableau("RosenbrockW6S4OSRodasTableau", T, T2)
