# MathOptInterface vector sets, status codes and the Clarabel MOI wrapper with
# the bridges ConvexOptimization reaches (bounds, Nonpositives, NormOneCone,
# NormInfinityCone, RotatedSecondOrderCone, PositiveSemidefiniteConeTriangle).
#
# Ported from MathOptInterface.jl (https://github.com/jump-dev/MathOptInterface.jl),
# version 1.54.0, Copyright (c) 2017: Miles Lubin and contributors, Copyright
# (c) 2017: Google Inc., MIT License (inst/CONVEX-LICENSES), and from the MOI
# wrapper of Clarabel.jl 0.11.1 (https://github.com/oxfordcontrol/Clarabel.jl),
# Copyright 2022 University of Oxford Control Group, Apache License 2.0
# (inst/CLARABEL-LICENSE).

.moi_set <- function(name, dimension, ...)
  structure(list(name = name, dimension = as.integer(dimension), ...), class = c(paste0("MOI_", name), "MOI_AbstractVectorSet"))

.moi_check_dim <- function(name, d, min = 0L) {
  if (!is.numeric(d) || length(d) != 1L || d != round(d) || d < min)
    stop("Dimension of ", name, " must be >= ", min, ", not ", format(d), ".", call. = FALSE)
  as.integer(d)
}

#' MathOptInterface sets and status codes
#'
#' `MOI` collects the vector sets of MathOptInterface that a [ConeConstraint()]
#' names, `MOI$dimension(set)`, and the status codes of a solve. The sets are
#' `Zeros(d)`, `Nonnegatives(d)`, `Nonpositives(d)`, `SecondOrderCone(d)`
#' (`(t, x)` with `norm(x) <= t`), `RotatedSecondOrderCone(d)` (`(t, u, x)`
#' with `2 t u >= sum(x^2)`, `t, u >= 0`), `NormOneCone(d)`,
#' `NormInfinityCone(d)`, `ExponentialCone()` (`(x, y, z)` with
#' `y exp(x / y) <= z`, `y > 0`), `DualExponentialCone()`, `PowerCone(a)`
#' (`x^a y^(1 - a) >= abs(z)`), `DualPowerCone(a)`,
#' `PositiveSemidefiniteConeTriangle(n)` (the upper triangle of a symmetric
#' `n` by `n` matrix, column by column), `ScaledPositiveSemidefiniteConeTriangle(n)`
#' (the same with the off diagonal entries scaled by `sqrt(2)`) and Clarabel's
#' `GenPowerCone(a, dim2)` (`prod(x^a) >= norm(z)` for `x` of length
#' `length(a)` and `z` of length `dim2`, `sum(a) == 1`).
#' `TerminationStatusCode` and `ResultStatusCode` list the status names that
#' `sol$original$termination_status`, `primal_status` and `dual_status` take.
#'
#' @format A list of set constructors, `dimension`, and the status code vectors.
#' @examples
#' MOI$dimension(MOI$SecondOrderCone(3))
#' MOI$dimension(MOI$PositiveSemidefiniteConeTriangle(3))
#' MOI$Zeros(2)
#' @export
MOI <- list(
  Zeros = function(dimension) .moi_set("Zeros", .moi_check_dim("Zeros", dimension)),
  Nonnegatives = function(dimension) .moi_set("Nonnegatives", .moi_check_dim("Nonnegatives", dimension)),
  Nonpositives = function(dimension) .moi_set("Nonpositives", .moi_check_dim("Nonpositives", dimension)),
  SecondOrderCone = function(dimension) .moi_set("SecondOrderCone", .moi_check_dim("SecondOrderCone", dimension)),
  RotatedSecondOrderCone = function(dimension)
    .moi_set("RotatedSecondOrderCone", .moi_check_dim("RotatedSecondOrderCone", dimension)),
  NormOneCone = function(dimension) .moi_set("NormOneCone", .moi_check_dim("NormOneCone", dimension, 1L)),
  NormInfinityCone = function(dimension) .moi_set("NormInfinityCone", .moi_check_dim("NormInfinityCone", dimension, 1L)),
  ExponentialCone = function() .moi_set("ExponentialCone", 3L),
  DualExponentialCone = function() .moi_set("DualExponentialCone", 3L),
  PowerCone = function(exponent) .moi_set("PowerCone", 3L, exponent = as.numeric(exponent)),
  DualPowerCone = function(exponent) .moi_set("DualPowerCone", 3L, exponent = as.numeric(exponent)),
  PositiveSemidefiniteConeTriangle = function(side_dimension) {
    n <- .moi_check_dim("PositiveSemidefiniteConeTriangle", side_dimension)
    .moi_set("PositiveSemidefiniteConeTriangle", (n * (n + 1L)) %/% 2L, side_dimension = n)
  },
  ScaledPositiveSemidefiniteConeTriangle = function(side_dimension) {
    n <- .moi_check_dim("ScaledPositiveSemidefiniteConeTriangle", side_dimension)
    .moi_set("ScaledPositiveSemidefiniteConeTriangle", (n * (n + 1L)) %/% 2L, side_dimension = n)
  },
  GenPowerCone = function(alpha, dim2) {
    alpha <- as.numeric(alpha)
    # Clarabel.jl: the exponents are positive and sum to one
    if (!all(alpha > 0) || abs(sum(alpha) - 1) > .Machine$double.eps * length(alpha) * 8)
      stop("GenPowerCone requires positive exponents that sum to one", call. = FALSE)
    .moi_set("GenPowerCone", length(alpha) + as.integer(dim2), alpha = alpha, dim2 = as.integer(dim2))
  },
  dimension = function(set) {
    if (!inherits(set, "MOI_AbstractVectorSet")) stop("`set` is not a MathOptInterface vector set", call. = FALSE)
    set$dimension
  },
  TerminationStatusCode = c("OPTIMIZE_NOT_CALLED", "OPTIMAL", "INFEASIBLE", "DUAL_INFEASIBLE", "LOCALLY_SOLVED",
                            "LOCALLY_INFEASIBLE", "INFEASIBLE_OR_UNBOUNDED", "ALMOST_OPTIMAL", "ALMOST_INFEASIBLE",
                            "ALMOST_DUAL_INFEASIBLE", "ALMOST_LOCALLY_SOLVED", "ITERATION_LIMIT", "TIME_LIMIT",
                            "NODE_LIMIT", "SOLUTION_LIMIT", "MEMORY_LIMIT", "OBJECTIVE_LIMIT", "NORM_LIMIT",
                            "OTHER_LIMIT", "SLOW_PROGRESS", "NUMERICAL_ERROR", "INVALID_MODEL", "INVALID_OPTION",
                            "INTERRUPTED", "OTHER_ERROR"),
  ResultStatusCode = c("NO_SOLUTION", "FEASIBLE_POINT", "NEARLY_FEASIBLE_POINT", "INFEASIBLE_POINT",
                       "INFEASIBILITY_CERTIFICATE", "NEARLY_INFEASIBILITY_CERTIFICATE", "REDUCTION_CERTIFICATE",
                       "NEARLY_REDUCTION_CERTIFICATE", "UNKNOWN_RESULT_STATUS", "OTHER_RESULT_STATUS")
)

.moi_format_set <- function(set) {
  extra <- switch(set$name,
    ExponentialCone = , DualExponentialCone = "",
    PowerCone = , DualPowerCone = format(set$exponent),
    PositiveSemidefiniteConeTriangle = , ScaledPositiveSemidefiniteConeTriangle = set$side_dimension,
    GenPowerCone = paste0(deparse(set$alpha), ", ", set$dim2),
    set$dimension)
  paste0("MOI$", set$name, "(", extra, ")")
}

#' @export
format.MOI_AbstractVectorSet <- function(x, ...) .moi_format_set(x)

#' @export
print.MOI_AbstractVectorSet <- function(x, ...) {
  cat(.moi_format_set(x), "\n", sep = "")
  invisible(x)
}

#' The Clarabel conic solver
#'
#' `ClarabelOptimizer(...)` is the interior point solver of Clarabel.jl as
#' an optimizer for [ConvexMOI()]: it solves
#' `min x'Px/2 + q'x` subject to `Ax + s = b`, `s` in a product of zero,
#' nonnegative, second order, exponential, power, generalized power and
#' positive semidefinite cones, with the homogeneous embedding, Ruiz
#' equilibration, a regularized LDL factorization (AMD ordering) with
#' iterative refinement and Mehrotra correction. `...` sets solver options
#' by their Clarabel names: `max_iter` (200), `time_limit` (`Inf`),
#' `tol_gap_abs`, `tol_gap_rel`, `tol_feas` (all `1e-8`), `tol_infeas_abs`,
#' `tol_infeas_rel` (`1e-8`), `tol_ktratio` (`1e-6`), the `reduced_tol_*`
#' tolerances, `max_step_fraction` (0.99), `equilibrate_enable` (`TRUE`),
#' `equilibrate_max_iter` (10), `equilibrate_min_scaling` (`1e-4`),
#' `equilibrate_max_scaling` (`1e4`), `linesearch_backtrack_step` (0.8),
#' `min_switch_step_length` (0.1), `min_terminate_step_length` (`1e-4`),
#' `static_regularization_enable`, `static_regularization_constant`
#' (`1e-8`), `dynamic_regularization_enable`, `dynamic_regularization_eps`
#' (`1e-13`), `dynamic_regularization_delta` (`2e-7`),
#' `iterative_refinement_enable`, the `iterative_refinement_*` settings and
#' `presolve_enable`.
#'
#' @param ... Solver options.
#' @return A `ClarabelOptimizer` object.
#' @examples
#' ClarabelOptimizer(max_iter = 50)$settings$max_iter
#' @export
ClarabelOptimizer <- function(...) {
  st <- .cl_default_settings()
  opts <- list(...)
  bad <- setdiff(names(opts), names(st))
  if (length(bad)) stop("Unrecognized option(s): ", paste(bad, collapse = ", "), call. = FALSE)
  st[names(opts)] <- opts
  structure(list(settings = st), class = "ClarabelOptimizer")
}

.moi_optimizer <- function(oc) {
  if (is.function(oc)) oc <- oc()
  if (!inherits(oc, "ClarabelOptimizer"))
    stop("`optimizer_constructor` must be ClarabelOptimizer or ClarabelOptimizer(...)", call. = FALSE)
  oc
}

# Termination, primal and dual status of each Clarabel status (Clarabel.jl
# MOI wrapper)
.CL_MOI_STATUS <- list(
  SOLVED = c("OPTIMAL", "FEASIBLE_POINT", "FEASIBLE_POINT"),
  MAX_ITERATIONS = c("ITERATION_LIMIT", "OTHER_RESULT_STATUS", "OTHER_RESULT_STATUS"),
  MAX_TIME = c("TIME_LIMIT", "OTHER_RESULT_STATUS", "OTHER_RESULT_STATUS"),
  PRIMAL_INFEASIBLE = c("INFEASIBLE", "INFEASIBLE_POINT", "INFEASIBILITY_CERTIFICATE"),
  DUAL_INFEASIBLE = c("DUAL_INFEASIBLE", "INFEASIBILITY_CERTIFICATE", "INFEASIBLE_POINT"),
  ALMOST_SOLVED = c("ALMOST_OPTIMAL", "NEARLY_FEASIBLE_POINT", "NEARLY_FEASIBLE_POINT"),
  ALMOST_PRIMAL_INFEASIBLE = c("ALMOST_INFEASIBLE", "INFEASIBLE_POINT", "NEARLY_INFEASIBILITY_CERTIFICATE"),
  ALMOST_DUAL_INFEASIBLE = c("ALMOST_DUAL_INFEASIBLE", "NEARLY_INFEASIBILITY_CERTIFICATE", "INFEASIBLE_POINT"),
  NUMERICAL_ERROR = c("NUMERICAL_ERROR", "OTHER_RESULT_STATUS", "OTHER_RESULT_STATUS"),
  INSUFFICIENT_PROGRESS = c("SLOW_PROGRESS", "OTHER_RESULT_STATUS", "OTHER_RESULT_STATUS"))

# ---------------------------------------------------------------------------
# The model: variables, constraints in insertion order, the bridges, and the
# copy into Clarabel's (P, q, A, b, cones)

.moi_model_new <- function(optimizer) {
  M <- new.env(parent = emptyenv())
  M$optimizer <- optimizer
  .moi_empty(M)
  M
}

.moi_empty <- function(M) {
  M$nvar <- 0L
  M$inner <- list()     # constraints in the solver's native sets
  M$outer <- list()     # constraints as added, with their bridge
  M$counts <- list()
  M$nbridge <- 0L
  M$obj_terms <- numeric(); M$obj_const <- 0; M$sense <- "FEASIBILITY_SENSE"
  M$termination_status <- "OPTIMIZE_NOT_CALLED"
  M$primal_status <- "NO_SOLUTION"; M$dual_status <- "NO_SOLUTION"
  M$result_count <- 0L
  M$solution <- NULL
  invisible(M)
}

.moi_add_variables <- function(M, k) {
  idx <- M$nvar + seq_len(k)
  M$nvar <- M$nvar + as.integer(k)
  idx
}

# F is a dense coefficient matrix over the variables vars, g the constants:
# the function F x[vars] + g in a native set
.moi_add_inner <- function(M, type, F, vars, g, extra = list()) {
  M$inner[[length(M$inner) + 1L]] <- c(list(type = type, F = F, vars = vars, g = g), extra)
  length(M$inner)
}

.moi_index <- function(M, key, bridged) {
  if (bridged) {
    M$nbridge <- M$nbridge + 1L
    return(list(type = key, value = -M$nbridge))
  }
  M$counts[[key]] <- (M$counts[[key]] %||% 0L) + 1L
  list(type = key, value = M$counts[[key]])
}

# Add F x[vars] + g in set; returns the constraint index of the outer model
.moi_add_constraint <- function(M, F, vars, g, set) {
  F <- matrix(F, nrow = length(g))
  d <- length(g)
  key <- paste0("VectorAffineFunction{Float64}-in-", set$name)
  rec <- list(set = set, inner = integer(), bridge = "none")
  native <- c(Zeros = "Zeros", Nonnegatives = "Nonnegatives", SecondOrderCone = "SecondOrderCone",
              ExponentialCone = "ExponentialCone", PowerCone = "PowerCone", GenPowerCone = "GenPowerCone",
              ScaledPositiveSemidefiniteConeTriangle = "ScaledPSD")
  switch(set$name,
    Zeros = , Nonnegatives = , SecondOrderCone = , ExponentialCone = , PowerCone = , GenPowerCone = ,
    ScaledPositiveSemidefiniteConeTriangle = {
      rec$inner <- .moi_add_inner(M, native[[set$name]], F, vars, g,
                                  list(alpha = set$exponent %||% set$alpha, dim2 = set$dim2, side = set$side_dimension))
    },
    Nonpositives = {
      rec$bridge <- "Nonpositives"
      rec$inner <- .moi_add_inner(M, "Nonnegatives", -F, vars, -g)
    },
    NormOneCone = {
      # (t, x) in NormOne(d): t - sum(y) >= 0, y - x >= 0, y + x >= 0
      rec$bridge <- "NormOneCone"
      k <- d - 1L
      y <- .moi_add_variables(M, k)
      Fx <- F[-1L, , drop = FALSE]; gx <- g[-1L]
      Fn <- rbind(cbind(F[1L, , drop = FALSE], matrix(-1, 1L, k)), cbind(-Fx, diag(1, k)), cbind(Fx, diag(1, k)))
      rec$inner <- .moi_add_inner(M, "Nonnegatives", Fn, c(vars, y), c(g[1L], -gx, gx))
    },
    NormInfinityCone = {
      # (t, x) in NormInf(d): t - x >= 0, t + x >= 0
      rec$bridge <- "NormInfinityCone"
      k <- d - 1L
      Ft <- F[rep(1L, k), , drop = FALSE]; Fx <- F[-1L, , drop = FALSE]
      rec$inner <- .moi_add_inner(M, "Nonnegatives", rbind(Ft - Fx, Ft + Fx), vars,
                                  c(rep(g[1L], k) - g[-1L], rep(g[1L], k) + g[-1L]))
    },
    RotatedSecondOrderCone = {
      # (t, u, x) to ((t + u) / sqrt(2), (t - u) / sqrt(2), x) in SOC
      rec$bridge <- "RotatedSecondOrderCone"
      s2 <- 1 / sqrt(2)
      Tm <- diag(1, d); Tm[1:2, 1:2] <- matrix(c(s2, s2, s2, -s2), 2L)
      rec$inner <- .moi_add_inner(M, "SecondOrderCone", Tm %*% F, vars, as.numeric(Tm %*% g))
    },
    PositiveSemidefiniteConeTriangle = {
      rec$bridge <- "PositiveSemidefiniteConeTriangle"
      sc <- .moi_psd_scaling(set$side_dimension)
      rec$inner <- .moi_add_inner(M, "ScaledPSD", F * sc, vars, g * sc, list(side = set$side_dimension))
    },
    stop("Constraints of type MathOptInterface.VectorAffineFunction{Float64}-in-", .moi_format_set(set),
         " are not supported by the solver, and there are no bridges that can reformulate it into supported ",
         "constraints.", call. = FALSE))
  M$outer[[length(M$outer) + 1L]] <- rec
  idx <- .moi_index(M, key, rec$bridge != "none")
  idx$outer <- length(M$outer)
  idx
}

# sqrt(2) on the off diagonal entries of an upper triangle stored column by column
.moi_psd_scaling <- function(n) {
  sc <- numeric((n * (n + 1L)) %/% 2L)
  k <- 1L
  for (j in seq_len(n)) for (i in seq_len(j)) {
    sc[k] <- if (i == j) 1 else sqrt(2)
    k <- k + 1L
  }
  sc
}

# Variable bounds: x >= lb is x - lb in Nonnegatives, x <= ub is -x + ub
.moi_add_bound <- function(M, var, value, lower) {
  rec <- list(set = NULL, bridge = "bound",
              inner = .moi_add_inner(M, "Nonnegatives", matrix(if (lower) 1 else -1, 1L, 1L), var,
                                     if (lower) -value else value))
  M$outer[[length(M$outer) + 1L]] <- rec
  .moi_index(M, if (lower) "VariableIndex-in-GreaterThan{Float64}" else "VariableIndex-in-LessThan{Float64}", TRUE)
}

.moi_set_objective <- function(M, coefs, constant, sense) {
  M$obj_terms <- coefs
  M$obj_const <- constant
  M$sense <- sense
}

.MOI_ROW_ORDER <- c("Zeros", "Nonnegatives", "SecondOrderCone", "ExponentialCone", "ScaledPSD", "PowerCone",
                    "GenPowerCone")

.moi_optimize <- function(M) {
  n <- M$nvar
  types <- vapply(M$inner, function(c) c$type, "")
  ord <- unlist(lapply(.MOI_ROW_ORDER, function(t) which(types == t)))
  I <- integer(); J <- integer(); V <- numeric(); b <- numeric()
  rows <- vector("list", length(M$inner))
  cones <- list()
  row0 <- 0L
  for (k in ord) {
    con <- M$inner[[k]]
    d <- length(con$g)
    nzr <- which(con$F != 0, arr.ind = TRUE)
    if (length(nzr)) {
      I <- c(I, row0 + nzr[, 1L]); J <- c(J, con$vars[nzr[, 2L]]); V <- c(V, -con$F[nzr])
    }
    b <- c(b, con$g)
    rows[[k]] <- row0 + seq_len(d)
    row0 <- row0 + d
    spec <- switch(con$type,
      Zeros = .cl_spec("ZeroConeT", d), Nonnegatives = .cl_spec("NonnegativeConeT", d),
      SecondOrderCone = .cl_spec("SecondOrderConeT", d), ExponentialCone = .cl_spec("ExponentialConeT"),
      ScaledPSD = .cl_spec("PSDTriangleConeT", con$side), PowerCone = .cl_spec("PowerConeT", alpha = con$alpha),
      GenPowerCone = .cl_spec("GenPowerConeT", alpha = con$alpha, dim2 = con$dim2))
    last <- length(cones)
    # the wrapper merges consecutive zero and nonnegative cones
    if (last && cones[[last]]$type == spec$type && spec$type %in% c("ZeroConeT", "NonnegativeConeT"))
      cones[[last]]$dim <- cones[[last]]$dim + d
    else cones[[last + 1L]] <- spec
  }
  m <- row0
  A <- .cl_sparse(I, J, V, m, n)
  q <- numeric(n)
  q[seq_along(M$obj_terms)] <- M$obj_terms
  const <- M$obj_const
  if (M$sense == "MAX_SENSE") {
    q <- -q
    const <- -const
  }
  P <- .csc(n, n, rep(1L, n + 1L), integer(), numeric())
  S <- .cl_setup(P, q, A, b, cones, M$optimizer$settings)
  sol <- .cl_solve(S)
  M$solver <- S
  M$cones <- cones
  M$rows <- rows
  M$A <- A; M$b <- b; M$q <- q
  M$solution <- sol
  M$objconstant <- const
  st <- .CL_MOI_STATUS[[sol$status]]
  M$termination_status <- st[1L]; M$primal_status <- st[2L]; M$dual_status <- st[3L]
  M$raw_status <- .CL_STATUS[[sol$status]]
  M$result_count <- 1L
  M$iterations <- sol$iterations
  M$solve_time <- sol$solve_time
  # as the Clarabel wrapper: the info cost, which is finite also on an infeasible exit
  rawobj <- S$info$cost_primal + const
  M$objective_value <- if (M$sense == "MAX_SENSE") -rawobj else rawobj
  invisible(M)
}

# ConstraintDual of an outer constraint, mapped back through its bridge
.moi_constraint_dual <- function(M, idx) {
  rec <- M$outer[[idx$outer]]
  z <- M$solution$z[M$rows[[rec$inner]]]
  switch(rec$bridge,
    none = z,
    Nonpositives = -z,
    NormOneCone = {
      k <- (length(z) - 1L) %/% 2L
      c(z[1L], z[k + 1L + seq_len(k)] - z[1L + seq_len(k)])
    },
    NormInfinityCone = {
      k <- length(z) %/% 2L
      c(.cl_ssum(z), z[k + seq_len(k)] - z[seq_len(k)])
    },
    RotatedSecondOrderCone = {
      s2 <- 1 / sqrt(2)
      c((z[1L] + z[2L]) * s2, (z[1L] - z[2L]) * s2, z[-(1:2)])
    },
    PositiveSemidefiniteConeTriangle = z / .moi_psd_scaling(rec$set$side_dimension),
    z)
}
