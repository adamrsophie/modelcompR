# MathOptInterface vector sets, status codes and the Clarabel MOI wrapper with
# the constraint bridges a ConeConstraint reaches: Nonpositives, NormOne,
# NormInfinity, RSOCtoSOC, SetDotScaling, GeoMeanToPower, DualGeoMean,
# RelativeEntropy, NormToPower, NormSpectral, NormNuclear, LogDet, RootDet,
# Square, SplitHyperRectangle, Vectorize and LessToGreater.
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

.moi_check_odd <- function(name, d) {
  if (!is.numeric(d) || length(d) != 1L || d != round(d) || d < 3 || d %% 2 == 0)
    stop("Dimension of ", name, " must be an odd integer >= 3, not ", format(d), ".", call. = FALSE)
  as.integer(d)
}

.moi_matrix_norm_set <- function(name, row_dim, column_dim) {
  if (row_dim < 0 || column_dim < 0) {
    if (name == "NormSpectralCone")
      stop("Dimensions of NormSpectralCone must be >= 0, not (", row_dim, ", ", column_dim, ").", call. = FALSE)
    stop("Dimension of NormNuclearCone must be >= 0, not (", row_dim, ", ", column_dim, ").", call. = FALSE)
  }
  r <- as.integer(row_dim); k <- as.integer(column_dim)
  .moi_set(name, 1L + r * k, row_dim = r, column_dim = k)
}

.moi_det_set <- function(name, side_dimension, offset, triangle) {
  if (!is.numeric(side_dimension) || side_dimension < 0)
    stop("Side dimension of ", name, " must be >= 0, not ", format(side_dimension), ".", call. = FALSE)
  n <- as.integer(side_dimension)
  .moi_set(name, offset + if (triangle) (n * (n + 1L)) %/% 2L else n * n, side_dimension = n)
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
#' `length(a)` and `z` of length `dim2`, `sum(a) == 1`). Sets Clarabel
#' lacks are reformulated by the bridges of MathOptInterface:
#' `NormCone(p, d)` (`(t, x)` with `norm(x, p) <= t`),
#' `GeometricMeanCone(d)` (`t <= prod(x)^(1 / (d - 1))`, `x >= 0`),
#' `DualGeometricMeanCone(d)`, `RelativeEntropyCone(d)` (`(u, v, w)` with
#' `sum(w * log(w / v)) <= u`), `NormSpectralCone(m, n)` and
#' `NormNuclearCone(m, n)` (`t` and an `m` by `n` matrix stored column by
#' column), `PositiveSemidefiniteConeSquare(n)`, `LogDetConeTriangle(n)` and
#' `LogDetConeSquare(n)` (`(t, u, X)` with `u logdet(X / u) >= t`),
#' `RootDetConeTriangle(n)`, `RootDetConeSquare(n)` (`(t, X)` with
#' `det(X)^(1 / n) >= t`) and `HyperRectangle(lower, upper)`.
#' `DualExponentialCone`, `DualPowerCone` and `DualRelativeEntropyCone` have
#' no reformulation for this solver, as in MathOptInterface.
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
  SecondOrderCone = function(dimension) .moi_set("SecondOrderCone", .moi_check_dim("SecondOrderCone", dimension, 1L)),
  RotatedSecondOrderCone = function(dimension)
    .moi_set("RotatedSecondOrderCone", .moi_check_dim("RotatedSecondOrderCone", dimension, 2L)),
  NormCone = function(p, dimension) {
    d <- .moi_check_dim("NormCone", dimension, 1L)
    if (!(p >= 1)) stop("ArgumentError: The p- argument to `NormCone` must be `>= 1`", call. = FALSE)
    .moi_set("NormCone", d, p = as.numeric(p))
  },
  GeometricMeanCone = function(dimension) .moi_set("GeometricMeanCone", .moi_check_dim("GeometricMeanCone", dimension, 2L)),
  DualGeometricMeanCone = function(dimension)
    .moi_set("DualGeometricMeanCone", .moi_check_dim("DualGeometricMeanCone", dimension, 2L)),
  RelativeEntropyCone = function(dimension) .moi_set("RelativeEntropyCone", .moi_check_odd("RelativeEntropyCone", dimension)),
  DualRelativeEntropyCone = function(dimension)
    .moi_set("DualRelativeEntropyCone", .moi_check_odd("RelativeEntropyCone", dimension)),
  NormSpectralCone = function(row_dim, column_dim) .moi_matrix_norm_set("NormSpectralCone", row_dim, column_dim),
  NormNuclearCone = function(row_dim, column_dim) .moi_matrix_norm_set("NormNuclearCone", row_dim, column_dim),
  PositiveSemidefiniteConeSquare = function(side_dimension) {
    n <- .moi_check_dim("PositiveSemidefiniteConeSquare", side_dimension)
    .moi_set("PositiveSemidefiniteConeSquare", n * n, side_dimension = n)
  },
  LogDetConeTriangle = function(side_dimension) .moi_det_set("LogDetConeTriangle", side_dimension, 2L, TRUE),
  LogDetConeSquare = function(side_dimension) .moi_det_set("LogDetConeSquare", side_dimension, 2L, FALSE),
  RootDetConeTriangle = function(side_dimension) .moi_det_set("RootDetConeTriangle", side_dimension, 1L, TRUE),
  RootDetConeSquare = function(side_dimension) .moi_det_set("RootDetConeSquare", side_dimension, 1L, FALSE),
  HyperRectangle = function(lower, upper) {
    if (length(lower) != length(upper))
      stop("DimensionMismatch: length of lower (=", length(lower), ") and upper (=", length(upper),
           ") bounds must match.", call. = FALSE)
    .moi_set("HyperRectangle", length(lower), lower = as.numeric(lower), upper = as.numeric(upper))
  },
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

# Julia's show of a Float64: the shortest digits that round trip, plain
# between 1e-4 and 1e6 and in scientific notation outside
.jl_show_float <- function(x) {
  if (is.nan(x)) return("NaN")
  if (is.infinite(x)) return(if (x > 0) "Inf" else "-Inf")
  if (x == 0) return(if (1 / x < 0) "-0.0" else "0.0")
  for (p in 1:17) {
    txt <- sprintf("%.*e", p - 1L, x)
    if (as.numeric(txt) == x) break
  }
  e <- as.integer(sub("^.*e", "", txt))
  d <- sub("0+$", "", gsub("[-.]", "", sub("e.*$", "", txt)))
  if (!nzchar(d)) d <- "0"
  sg <- if (x < 0) "-" else ""
  if (e >= -4 && e <= 5) {
    if (e < 0) return(paste0(sg, "0.", strrep("0", -e - 1L), d))
    n <- nchar(d)
    ip <- if (n > e + 1L) substr(d, 1L, e + 1L) else paste0(d, strrep("0", e + 1L - n))
    fp <- if (n > e + 1L) substr(d, e + 2L, n) else "0"
    return(paste0(sg, ip, ".", fp))
  }
  paste0(sg, substr(d, 1L, 1L), ".", if (nchar(d) > 1L) substr(d, 2L, nchar(d)) else "0", "e", e)
}
.jl_show_vec <- function(x) paste0("[", paste(vapply(x, .jl_show_float, ""), collapse = ", "), "]")

# a set as Julia shows it, so the messages read as in ConvexOptimization.jl
.moi_format_set <- function(set) {
  extra <- switch(set$name,
    ExponentialCone = , DualExponentialCone = "",
    PowerCone = , DualPowerCone = .jl_show_float(set$exponent),
    NormCone = paste0(.jl_show_float(set$p), ", ", set$dimension),
    NormSpectralCone = , NormNuclearCone = paste0(set$row_dim, ", ", set$column_dim),
    HyperRectangle = paste0(.jl_show_vec(set$lower), ", ", .jl_show_vec(set$upper)),
    PositiveSemidefiniteConeTriangle = , ScaledPositiveSemidefiniteConeTriangle = , PositiveSemidefiniteConeSquare = ,
    LogDetConeTriangle = , LogDetConeSquare = , RootDetConeTriangle = , RootDetConeSquare = set$side_dimension,
    GenPowerCone = paste0(.jl_show_vec(set$alpha), ", ", set$dim2),
    set$dimension)
  typ <- if (set$name %in% c("PowerCone", "DualPowerCone", "GenPowerCone", "HyperRectangle")) "{Float64}" else ""
  paste0("MathOptInterface.", set$name, typ, "(", extra, ")")
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
  M$nodes <- character(); M$lessthan <- NULL; M$active <- character()
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

# An affine function F x[vars] + g: F is dense over the variables vars
.moi_aff <- function(F, vars, g) list(F = matrix(as.numeric(F), length(g), length(vars)), vars = as.integer(vars),
                                      g = as.numeric(g))
.moi_aff_vars <- function(v) .moi_aff(diag(1, length(v)), v, numeric(length(v)))
.moi_aff_rows <- function(a, i) .moi_aff(a$F[i, , drop = FALSE], a$vars, a$g[i])
.moi_aff_scale <- function(a, s) .moi_aff(a$F * s, a$vars, a$g * s)
.moi_aff_widen <- function(a, vars) {
  F <- matrix(0, length(a$g), length(vars))
  F[, match(a$vars, vars)] <- a$F
  F
}
.moi_aff_vcat <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  vars <- unique(unlist(lapply(parts, function(p) p$vars)))
  .moi_aff(do.call(rbind, lapply(parts, .moi_aff_widen, vars = vars)), vars,
           unlist(lapply(parts, function(p) p$g)))
}
# a + s b for two functions of the same output dimension
.moi_aff_add <- function(a, b, s = 1) {
  vars <- unique(c(a$vars, b$vars))
  .moi_aff(.moi_aff_widen(a, vars) + s * .moi_aff_widen(b, vars), vars, a$g + s * b$g)
}
.moi_aff_sum <- function(v) .moi_aff(matrix(1, 1L, length(v)), v, 0)

# the position of entry (i, j) in the upper triangle stored column by column
.moi_trimap <- function(i, j) if (i < j) .moi_trimap(j, i) else j + ((i - 1L) * i) %/% 2L

# a constraint in a set Clarabel takes as is; a VectorOfVariables function
# is listed before every VectorAffineFunction when the model is copied
.moi_add_inner <- function(M, type, a, vov = FALSE, extra = list()) {
  M$inner[[length(M$inner) + 1L]] <- c(list(type = type, F = a$F, vars = a$vars, g = a$g, vov = vov), extra)
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
  key <- paste0("VectorAffineFunction{Float64}-in-", set$name)
  node <- .moi_bridge_vector(M, .moi_aff(F, vars, g), set)
  M$outer[[length(M$outer) + 1L]] <- node
  idx <- .moi_index(M, key, node$bridge != "none")
  idx$outer <- length(M$outer)
  idx
}

.moi_native <- function(M, a, set, vov = FALSE) {
  type <- switch(set$name, ScaledPositiveSemidefiniteConeTriangle = "ScaledPSD", set$name)
  list(set = set, bridge = "none",
       inner = .moi_add_inner(M, type, a, vov, list(alpha = set$exponent %||% set$alpha, dim2 = set$dim2,
                                                    side = set$side_dimension)))
}

# active bridges are named as print_active_bridges names them, less "MOIB."
.moi_bridged <- function(M, set, bridge, kids, ...) {
  M$active <- union(M$active, paste0("Constraint.", bridge))
  list(set = set, bridge = bridge, kids = kids, ...)
}

# MathOptInterface's bridge graph picks, for each set Clarabel lacks, the
# chain of least total cost (one per bridge, nothing for a native set),
# ties going to the bridge registered first. The costs for this solver:
# Nonpositives, PSDTriangle, RSOC, NormOne, NormInfinity, GeometricMean (by
# GeoMeanToPowerBridge; GeoMeantoRelEntrBridge costs 3, GeoMeanBridge 4) and
# HyperRectangle 1; NormCone, RelativeEntropy, NormSpectral and the scalar
# GreaterThan and EqualTo 1 or 2; NormNuclear, RootDetTriangle and
# PSDSquare 3; LogDetTriangle 4 and the square forms one more than their
# triangle plus an equality. DualExponentialCone and DualPowerCone have no
# bridge, so DualRelativeEntropyCone (which bridges to the former) has none
# either.
.moi_bridge_vector <- function(M, a, set) {
  d <- length(a$g)
  switch(set$name,
    Zeros = , Nonnegatives = , SecondOrderCone = , ExponentialCone = , PowerCone = , GenPowerCone = ,
    ScaledPositiveSemidefiniteConeTriangle = .moi_native(M, a, set),
    Nonpositives = {
      M$nodes <- union(M$nodes, "Nonpositives")
      .moi_bridged(M, set, "NonposToNonnegBridge",
                   list(.moi_native(M, .moi_aff_scale(a, -1), MOI$Nonnegatives(d))))
    },
    NormOneCone = {
      # (t, x): t - sum(y) >= 0, y - x >= 0, y + x >= 0
      k <- d - 1L
      y <- .moi_aff_vars(.moi_add_variables(M, k))
      t <- .moi_aff_rows(a, 1L); x <- .moi_aff_rows(a, -1L)
      f <- .moi_aff_vcat(.moi_aff_add(t, .moi_aff_sum(y$vars), -1), .moi_aff_add(y, x, -1), .moi_aff_add(y, x))
      .moi_bridged(M, set, "NormOneBridge", list(.moi_native(M, f, MOI$Nonnegatives(2L * k + 1L))))
    },
    NormInfinityCone = {
      # (t, x): t - x >= 0, t + x >= 0
      k <- d - 1L
      t <- .moi_aff_rows(a, rep(1L, k)); x <- .moi_aff_rows(a, -1L)
      .moi_bridged(M, set, "NormInfinityBridge",
                   list(.moi_native(M, .moi_aff_vcat(.moi_aff_add(t, x, -1), .moi_aff_add(t, x)), MOI$Nonnegatives(2L * k))))
    },
    RotatedSecondOrderCone = {
      # (t, u, x) to ((t + u) / sqrt(2), (t - u) / sqrt(2), x) in SOC
      s2 <- 1 / sqrt(2)
      Tm <- diag(1, d); Tm[1:2, 1:2] <- matrix(c(s2, s2, s2, -s2), 2L)
      .moi_bridged(M, set, "RSOCtoSOCBridge",
                   list(.moi_native(M, .moi_aff(Tm %*% a$F, a$vars, as.numeric(Tm %*% a$g)), MOI$SecondOrderCone(d))))
    },
    PositiveSemidefiniteConeTriangle = {
      sc <- .moi_psd_scaling(set$side_dimension)
      .moi_bridged(M, set, "SetDotScalingBridge",
                   list(.moi_native(M, .moi_aff(a$F * sc, a$vars, a$g * sc),
                                    MOI$ScaledPositiveSemidefiniteConeTriangle(set$side_dimension))))
    },
    GeometricMeanCone = .moi_geomean_to_power(M, a, set),
    DualGeometricMeanCone = {
      f1 <- .moi_aff_rows(a, 1L)
      nn <- .moi_native(M, .moi_aff_scale(f1, -1), MOI$Nonnegatives(1L))
      gm <- .moi_bridge_vector(M, .moi_aff_vcat(.moi_aff_scale(f1, 1 / -(d - 1)), .moi_aff_rows(a, -1L)),
                               MOI$GeometricMeanCone(d))
      .moi_bridged(M, set, "DualGeoMeanBridge", list(nn, gm))
    },
    RelativeEntropyCone = {
      v <- (d - 1L) %/% 2L
      y <- .moi_add_variables(M, v)
      ge <- .moi_bridge_scalar(M, .moi_aff_add(.moi_aff_rows(a, 1L), .moi_aff_sum(y), -1), "GreaterThan")
      ex <- lapply(seq_len(v), function(i)
        .moi_native(M, .moi_aff_vcat(.moi_aff_scale(.moi_aff_vars(y[i]), -1), .moi_aff_rows(a, 1L + i + v),
                                     .moi_aff_rows(a, 1L + i)), MOI$ExponentialCone()))
      .moi_bridged(M, set, "RelativeEntropyBridge", c(list(ge), ex), v = v)
    },
    NormCone = {
      # sum(r) = t with (r_i, t, x_i) in PowerCone(1 / p)
      r <- .moi_add_variables(M, d - 1L)
      t <- .moi_aff_rows(a, 1L)
      pw <- lapply(seq_along(r), function(i)
        .moi_native(M, .moi_aff_vcat(.moi_aff_vars(r[i]), t, .moi_aff_rows(a, i + 1L)), MOI$PowerCone(1 / set$p)))
      eq <- .moi_bridge_scalar(M, .moi_aff_add(.moi_aff_sum(r), t, -1), "EqualTo")
      .moi_bridged(M, set, "NormToPowerBridge", c(pw, list(eq)))
    },
    NormSpectralCone = .moi_norm_spectral(M, a, set),
    NormNuclearCone = .moi_norm_nuclear(M, a, set),
    LogDetConeTriangle = .moi_logdet(M, a, set),
    RootDetConeTriangle = .moi_rootdet(M, a, set),
    PositiveSemidefiniteConeSquare = , LogDetConeSquare = , RootDetConeSquare = .moi_square(M, a, set),
    HyperRectangle = .moi_hyperrectangle(M, a, set),
    stop("Constraints of type MathOptInterface.VectorAffineFunction{Float64}-in-", .moi_format_set(set),
         " are not supported by the solver, and there are no bridges that can reformulate it into supported ",
         "constraints.", call. = FALSE))
}

# A scalar f in GreaterThan(0), LessThan(0) or EqualTo(0), the constant of
# f already moved by normalize_and_add_constraint. LessThan has two chains
# of cost 2: LessToGreaterBridge, registered first, wins unless the
# Nonpositives node is in the graph and the GreaterThan node is not when the
# LessThan node is created, in which case Bellman-Ford settles
# VectorizeBridge one sweep earlier.
.moi_bridge_scalar <- function(M, a, type) {
  switch(type,
    GreaterThan = {
      M$nodes <- union(M$nodes, "GreaterThan")
      .moi_bridged(M, type, "VectorizeBridge", list(.moi_native(M, a, MOI$Nonnegatives(1L))))
    },
    EqualTo = .moi_bridged(M, type, "VectorizeBridge", list(.moi_native(M, a, MOI$Zeros(1L)))),
    LessThan = {
      if (is.null(M$lessthan)) {
        M$lessthan <- if ("Nonpositives" %in% M$nodes && !("GreaterThan" %in% M$nodes)) "VectorizeBridge"
          else "LessToGreaterBridge"
        M$nodes <- union(M$nodes, c("GreaterThan", "Nonpositives"))
      }
      if (M$lessthan == "VectorizeBridge")
        .moi_bridged(M, type, "VectorizeBridge", list(.moi_bridge_vector(M, a, MOI$Nonpositives(1L))))
      else .moi_bridged(M, type, "LessToGreaterBridge", list(.moi_bridge_scalar(M, .moi_aff_scale(a, -1), "GreaterThan")))
    })
}

.moi_geomean_to_power <- function(M, a, set) {
  d <- length(a$g)
  fi <- function(i) .moi_aff_rows(a, i)
  if (d <= 3L) {
    pw <- .moi_native(M, .moi_aff_rows(a, if (d == 2L) c(2L, 2L, 1L) else c(2L, 3L, 1L)), MOI$PowerCone(1 / 2))
    return(.moi_bridged(M, set, "GeoMeanToPowerBridge", list(pw), d = d))
  }
  ts <- .moi_add_variables(M, d - 3L)
  nn <- .moi_native(M, .moi_aff_vars(ts), MOI$Nonnegatives(d - 3L), vov = TRUE)
  pws <- list()
  z <- fi(1L); x <- fi(2L); n <- d - 1L
  for (i in seq_len(d - 3L)) {
    y <- .moi_aff_vars(ts[i])
    pws[[i]] <- .moi_native(M, .moi_aff_vcat(x, y, z), MOI$PowerCone(1 / n))
    n <- n - 1L
    x <- fi(2L + i); z <- y
  }
  pws[[d - 2L]] <- .moi_native(M, .moi_aff_vcat(fi(d - 1L), fi(d), .moi_aff_vars(ts[d - 3L])), MOI$PowerCone(1 / 2))
  .moi_bridged(M, set, "GeoMeanToPowerBridge", pws, d = d, nn = nn)
}

# [t I, X'; X, t I] in PSDTriangle(m + n) for the m by n matrix X
.moi_norm_spectral <- function(M, a, set) {
  m <- set$row_dim; n <- set$column_dim; side <- m + n
  F <- matrix(0, (side * (side + 1L)) %/% 2L, length(a$vars)); g <- numeric(nrow(F))
  for (i in seq_len(side)) {
    k <- .moi_trimap(i, i)
    F[k, ] <- F[k, ] + a$F[1L, ]; g[k] <- g[k] + a$g[1L]
  }
  r <- 2L
  for (j in seq_len(n)) for (i in n + seq_len(m)) {
    k <- .moi_trimap(i, j)
    F[k, ] <- F[k, ] + a$F[r, ]; g[k] <- g[k] + a$g[r]
    r <- r + 1L
  }
  psd <- .moi_bridge_vector(M, .moi_aff(F, a$vars, g), MOI$PositiveSemidefiniteConeTriangle(side))
  .moi_bridged(M, set, "NormSpectralBridge", list(psd))
}

# t >= (tr U + tr V) / 2 with [U, X'; X, V] PSD
.moi_norm_nuclear <- function(M, a, set) {
  m <- set$row_dim; n <- set$column_dim
  U <- .moi_add_variables(M, (n * (n + 1L)) %/% 2L)
  V <- .moi_add_variables(M, (m * (m + 1L)) %/% 2L)
  dg <- c(U[vapply(seq_len(n), function(i) .moi_trimap(i, i), 0L)], V[vapply(seq_len(m), function(i) .moi_trimap(i, i), 0L)])
  ge <- .moi_bridge_scalar(M, .moi_aff_add(.moi_aff_rows(a, 1L), .moi_aff(matrix(1, 1L, length(dg)) / 2, dg, 0), -1),
                           "GreaterThan")
  parts <- list(.moi_aff_vars(U))
  for (i in seq_len(m)) {
    parts[[length(parts) + 1L]] <- .moi_aff_rows(a, seq.int(1L + i, 1L + m * n, by = m))
    parts[[length(parts) + 1L]] <- .moi_aff_vars(V[vapply(seq_len(i), function(j) .moi_trimap(i, j), 0L)])
  }
  psd <- .moi_bridge_vector(M, do.call(.moi_aff_vcat, parts), MOI$PositiveSemidefiniteConeTriangle(m + n))
  .moi_bridged(M, set, "NormNuclearBridge", list(ge, psd), m = m, n = n)
}

# [X, D'; D, diag(D)] in PSDTriangle(2 d) for a lower triangular D of new
# variables; returns the eigenvalue variables diag(D)
.moi_extract_eigenvalues <- function(M, a, d, offset) {
  n <- (d * (d + 1L)) %/% 2L
  X <- .moi_aff_rows(a, offset + seq_len(n))
  Delta <- .moi_add_variables(M, n)
  N <- (2L * d * (2L * d + 1L)) %/% 2L
  zv <- rep(NA_integer_, N - n)
  for (j in seq_len(d)) {
    for (i in j:d) zv[.moi_trimap(i, d + j) - n] <- Delta[.moi_trimap(i, j)]
    zv[.moi_trimap(d + j, d + j) - n] <- Delta[.moi_trimap(j, j)]
  }
  Fz <- matrix(0, N - n, length(Delta))
  Fz[cbind(which(!is.na(zv)), match(zv[!is.na(zv)], Delta))] <- 1
  Y <- .moi_aff_vcat(X, .moi_aff(Fz, Delta, numeric(N - n)))
  sd <- .moi_bridge_vector(M, Y, MOI$PositiveSemidefiniteConeTriangle(2L * d))
  list(sd = sd, D = Delta[vapply(seq_len(d), function(i) .moi_trimap(i, i), 0L)], n = n)
}

.moi_logdet <- function(M, a, set) {
  d <- set$side_dimension
  ev <- .moi_extract_eigenvalues(M, a, d, 2L)
  l <- .moi_add_variables(M, d)
  u <- .moi_aff_rows(a, 2L)
  ex <- lapply(seq_len(d), function(i)
    .moi_native(M, .moi_aff_vcat(.moi_aff_vars(l[i]), u, .moi_aff_vars(ev$D[i])), MOI$ExponentialCone()))
  tl <- .moi_bridge_scalar(M, .moi_aff_add(.moi_aff_rows(a, 1L), .moi_aff_sum(l), -1), "LessThan")
  .moi_bridged(M, set, "LogDetBridge", c(list(ev$sd, tl), ex), n = ev$n)
}

.moi_rootdet <- function(M, a, set) {
  d <- set$side_dimension
  ev <- .moi_extract_eigenvalues(M, a, d, 1L)
  gm <- .moi_bridge_vector(M, .moi_aff_vcat(.moi_aff_rows(a, 1L), .moi_aff_vars(ev$D)), MOI$GeometricMeanCone(d + 1L))
  .moi_bridged(M, set, "RootDetBridge", list(ev$sd, gm), n = ev$n)
}

# the upper triangle of a square matrix set, with f_ij = f_ji added as an
# equality wherever the two entries differ
.moi_square <- function(M, a, set) {
  dim <- set$side_dimension
  offset <- switch(set$name, LogDetConeSquare = 2L, RootDetConeSquare = 1L, 0L)
  idx <- seq_len(offset)
  sym <- list()
  for (j in seq_len(dim)) for (i in seq_len(j)) {
    idx <- c(idx, offset + i + (j - 1L) * dim)
    if (i == j) next
    diff <- .moi_aff_add(.moi_aff_rows(a, offset + i + (j - 1L) * dim), .moi_aff_rows(a, offset + j + (i - 1L) * dim), -1)
    big <- max(abs(c(diff$F, diff$g)), 0)
    if (big <= 1e-10) next
    if (big <= 1e-8)
      warning("The entries (", i, ", ", j, ") and (", j, ", ", i, ") of the matrix are almost identical, but a ",
              "constraint has been added to ensure their equality because the largest difference between the ",
              "coefficients is smaller than 1e-8 but larger than 1e-10. This usually means that there is a ",
              "modeling error in your formulation.", call. = FALSE)
    sym[[length(sym) + 1L]] <- list(ij = c(i, j), node = .moi_bridge_scalar(M, diff, "EqualTo"))
  }
  tset <- switch(set$name, LogDetConeSquare = MOI$LogDetConeTriangle(dim), RootDetConeSquare = MOI$RootDetConeTriangle(dim),
                 MOI$PositiveSemidefiniteConeTriangle(dim))
  tri <- .moi_bridge_vector(M, .moi_aff_rows(a, idx), tset)
  # print_active_bridges lists the bridge graph, which holds the EqualTo
  # chain even when every pair is symmetric and no equality is added
  M$active <- union(M$active, "Constraint.VectorizeBridge")
  .moi_bridged(M, set, "SquareBridge", list(tri), sym = sym, dim = dim, offset = offset)
}

.moi_hyperrectangle <- function(M, a, set) {
  N <- length(a$g)
  lo <- lapply(seq_len(N), function(i) if (is.finite(set$lower[i]))
    .moi_aff_add(.moi_aff_rows(a, i), .moi_aff(matrix(0, 1L, 0L), integer(), set$lower[i]), -1))
  up <- lapply(seq_len(N), function(i) if (is.finite(set$upper[i]))
    .moi_aff_add(.moi_aff(matrix(0, 1L, 0L), integer(), set$upper[i]), .moi_aff_rows(a, i), -1))
  rows <- c(lo, up)
  rows <- rows[!vapply(rows, is.null, TRUE)]
  if (!length(rows)) return(.moi_bridged(M, set, "SplitHyperRectangleBridge", list()))
  .moi_bridged(M, set, "SplitHyperRectangleBridge",
               list(.moi_native(M, do.call(.moi_aff_vcat, rows), MOI$Nonnegatives(length(rows)))))
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

# Variable bounds: x >= lb is x - lb in Nonnegatives, x <= ub is -x + ub.
# The bridge graph creates the GreaterThan and LessThan nodes on the way.
# The bound is a constrained variable: GreaterThan through
# Variable.VectorizeBridge, LessThan through it and then
# Variable.NonposToNonnegBridge
.moi_add_bound <- function(M, var, value, lower) {
  M$active <- union(M$active, c("Variable.VectorizeBridge", if (!lower) "Variable.NonposToNonnegBridge"))
  M$nodes <- union(M$nodes, "GreaterThan")
  if (!lower && is.null(M$lessthan)) {
    M$lessthan <- "LessToGreaterBridge"
    M$nodes <- union(M$nodes, "Nonpositives")
  }
  node <- list(set = NULL, bridge = "bound",
               kids = list(.moi_native(M, .moi_aff(if (lower) 1 else -1, var, if (lower) -value else value),
                                       MOI$Nonnegatives(1L))))
  M$outer[[length(M$outer) + 1L]] <- node
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
  vov <- vapply(M$inner, function(c) c$vov, TRUE)
  # the cache lists VectorOfVariables constraints before VectorAffineFunction
  # ones, each by set in the order of MOI.Utilities.Model
  ord <- unlist(lapply(c(TRUE, FALSE), function(v)
    unlist(lapply(.MOI_ROW_ORDER, function(t) which(types == t & vov == v)))))
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

# ConstraintDual of a constraint, mapped back through its bridges
.moi_node_dual <- function(M, node) {
  kd <- function(k) .moi_node_dual(M, node$kids[[k]])
  seqsum <- function(v) Reduce(`+`, v, 0)
  switch(node$bridge,
    none = M$solution$z[M$rows[[node$inner]]],
    NonposToNonnegBridge = , LessToGreaterBridge = -kd(1L),
    VectorizeBridge = kd(1L),
    NormOneBridge = {
      z <- kd(1L)
      k <- (length(z) - 1L) %/% 2L
      c(z[1L], z[k + 1L + seq_len(k)] - z[1L + seq_len(k)])
    },
    NormInfinityBridge = {
      z <- kd(1L)
      k <- length(z) %/% 2L
      c(.cl_ssum(z), z[k + seq_len(k)] - z[seq_len(k)])
    },
    RSOCtoSOCBridge = {
      z <- kd(1L)
      s2 <- 1 / sqrt(2)
      c((z[1L] + z[2L]) * s2, (z[1L] - z[2L]) * s2, z[-(1:2)])
    },
    SetDotScalingBridge = kd(1L) / .moi_psd_scaling(node$set$side_dimension),
    GeoMeanToPowerBridge = {
      z <- kd(1L)
      if (node$d == 2L) return(c(z[3L], z[1L] + z[2L]))
      if (node$d == 3L) return(z[c(3L, 1L, 2L)])
      np <- length(node$kids)
      g <- z[c(3L, 1L)]
      for (i in seq_len(np - 2L) + 1L) g <- c(g, kd(i)[1L])
      c(g, kd(np)[1:2])
    },
    DualGeoMeanBridge = {
      z <- kd(2L)
      z[1L] <- z[1L] / -(length(z) - 1)
      z
    },
    RelativeEntropyBridge = {
      v <- node$v
      out <- numeric(1L + 2L * v)
      out[1L] <- kd(1L)
      for (i in seq_len(v)) {
        e <- kd(1L + i)
        out[1L + i] <- e[3L]; out[1L + i + v] <- e[2L]
      }
      out
    },
    NormSpectralBridge = , NormNuclearBridge = {
      st <- node$set
      m <- st$row_dim; n <- st$column_dim
      spectral <- node$bridge == "NormSpectralBridge"
      z <- kd(if (spectral) 1L else 2L)
      ix <- unlist(lapply(seq_len(n), function(j) vapply(n + seq_len(m), function(i) .moi_trimap(i, j), 0L)))
      t <- if (spectral) seqsum(z[vapply(seq_len(m + n), function(i) .moi_trimap(i, i), 0L)]) else kd(1L)
      c(t, 2 * z[ix])
    },
    LogDetBridge = {
      d <- node$set$side_dimension
      c(kd(2L), seqsum(vapply(seq_len(d), function(i) kd(2L + i)[2L], 0)), kd(1L)[seq_len(node$n)])
    },
    RootDetBridge = c(kd(2L)[1L], kd(1L)[seq_len(node$n)]),
    SquareBridge = {
      tri <- kd(1L)
      dim <- node$dim; offset <- node$offset
      out <- numeric(offset + dim^2)
      out[seq_len(offset)] <- tri[seq_len(offset)]
      k <- offset; s <- 1L
      for (j in seq_len(dim)) for (i in seq_len(j)) {
        k <- k + 1L
        up <- offset + i + (j - 1L) * dim; lo <- offset + j + (i - 1L) * dim
        if (i == j) out[up] <- tri[k]
        else if (s <= length(node$sym) && identical(node$sym[[s]]$ij, c(i, j))) {
          piv <- .moi_node_dual(M, node$sym[[s]]$node)
          out[up] <- 2 * tri[k] + piv
          out[lo] <- -piv
          s <- s + 1L
        } else out[up] <- out[lo] <- tri[k]
      }
      out
    },
    SplitHyperRectangleBridge = {
      st <- node$set
      out <- numeric(length(st$lower))
      if (!length(node$kids)) return(out)
      v <- kd(1L)
      r <- 0L
      for (i in seq_along(st$lower)) if (is.finite(st$lower[i])) { r <- r + 1L; out[i] <- out[i] + v[r] }
      for (i in seq_along(st$upper)) if (is.finite(st$upper[i])) { r <- r + 1L; out[i] <- out[i] - v[r] }
      out
    },
    # MathOptInterface implements no dual for this bridge
    NormToPowerBridge = stop("GetAttributeNotAllowed{MathOptInterface.ConstraintDual}: Getting attribute ",
      "MathOptInterface.ConstraintDual(1) cannot be performed because: Bridge of type `NormToPowerBridge` does not ",
      "support accessing the attribute `MathOptInterface.ConstraintDual(1)`. If you encountered this error ",
      "unexpectedly, it probably means your model has been reformulated using the bridge, and you are attempting ",
      "to query an attribute that we haven't implemented yet for this bridge.", call. = FALSE),
    stop("no dual for bridge ", node$bridge, call. = FALSE))
}

.moi_constraint_dual <- function(M, idx) .moi_node_dual(M, M$outer[[idx$outer]])
