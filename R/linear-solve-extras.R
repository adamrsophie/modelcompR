# Parts of LinearSolve.jl reexported by NonlinearSolve: the WarmStart and
# NonstructuralZeros options (src/iterative_wrappers.jl, src/common.jl),
# PartitionedSolversAlgorithm (src/extension_algs.jl), update_gamma! (src/lhl.jl),
# the supernodal panel solves (src/SupernodalLU/solve.jl), slate_isavailable
# (src/slate.jl) and show_algorithm_choices (src/preferences.jl).

#' Linear solver options and helpers
#'
#' `WarmStart` lists the initial guesses of the Krylov methods for a
#' `LinearCache` solved repeatedly, passed as `KrylovJL_GMRES(warm_start =
#' WarmStart$Hegedus)`: `Auto` and `None` start from zero, `Previous` from the
#' previous solution and `Hegedus` from the previous solution rescaled to
#' minimize the initial residual (used only when it halves the residual). The
#' stopping threshold stays the one of a cold start.
#'
#' `NonstructuralZeros` tells `OperatorAssumptions()` how stored zeros of a
#' sparse matrix behave over a sequence of solves: `None` keeps the matrix,
#' `Present` and `Persistent` drop the stored zeros of each matrix, and `Auto`
#' (the default) drops them when at least a tenth of the stored entries are
#' zero. The direct solvers of R factorize sparse matrices as dense ones, so
#' the option is recorded in the assumptions without changing the solution.
#'
#' `PartitionedSolversAlgorithm(solver)` is the algorithm of
#' PartitionedSolvers.jl: the default solver is a direct LU solve, and
#' `solver = "cg"` (or another Krylov method name) an iterative one.
#'
#' `update_gamma(cache, gamma)` sets the `gamma` of the `WOperator` of a
#' `LinearCache` (the system matrix `J - M / gamma`) and invalidates the
#' cached factorization.
#'
#' `supernodal_panel_solve(W, B, np, operation, algorithm)` applies a
#' triangular solve with the factored diagonal block `W[1:np, 1:np]` of a
#' supernodal LU panel to `B` and returns the result: `"factor_right_upper"`
#' (`B / U11`), `"factor_lower"` and `"lower"` (`L11 \ B`, unit lower) and
#' `"upper"` (`U11 \ B`). `supernodal_panel_solve_backend()` runs one backend
#' (`"kernel"`, column oriented substitution; `"blas"` and
#' `"triangularsolve"`, the triangular solves of R); `"auto"` picks the kernel
#' for blocks up to 256 wide or a single column, as LinearSolve.
#'
#' `slate_isavailable()` is `TRUE`: `SLATEFactorization()` runs as an LU
#' factorization in R. `show_algorithm_choices()` prints the default
#' algorithm for dense square matrices of several sizes.
#'
#' @param solver A Krylov method name, or `NULL` for the direct solver.
#' @param ... Options of the solver.
#' @param cache A `LinearCache` from `init()`.
#' @param gamma The new shift.
#' @param W A panel whose leading `np` by `np` block holds the factors.
#' @param B The right-hand side panel.
#' @param np The width of the diagonal block.
#' @param operation The triangular operation.
#' @param algorithm The backend.
#' @param libpath Ignored; the path of a SLATE library.
#' @return `WarmStart` and `NonstructuralZeros` are lists of names;
#'   `PartitionedSolversAlgorithm()` a linear algorithm; `update_gamma()` the
#'   cache; the panel solves the updated `B`; `slate_isavailable()` a logical;
#'   `show_algorithm_choices()` is called for its output.
#' @examples
#' A <- diag(4) + 0.1; b <- 1:4
#' cache <- init(LinearProblem(A, b), KrylovJL_GMRES(warm_start = WarmStart$Hegedus))
#' solve_integrator(cache)$u
#' solve(LinearProblem(A, b), PartitionedSolversAlgorithm("cg"))$u
#' W <- matrix(c(2, 0.5, 1, 3), 2)
#' supernodal_panel_solve(W, matrix(1:4, 2), 2, "upper")
#' OperatorAssumptions(TRUE, nonstructural_zeros = NonstructuralZeros$Present)$nonstructural_zeros
#' @name linear_solve_extras
NULL

#' @rdname linear_solve_extras
#' @export
WarmStart <- list(Auto = "Auto", None = "None", Previous = "Previous", Hegedus = "Hegedus")

#' @rdname linear_solve_extras
#' @export
NonstructuralZeros <- list(Auto = "Auto", None = "None", Persistent = "Persistent", Present = "Present")

#' @rdname linear_solve_extras
#' @export
PartitionedSolversAlgorithm <- function(solver = NULL, ...) {
  if (is.null(solver)) return(.lin_alg("PartitionedSolversAlgorithm", "lu", solver = NULL, ...))
  .lin_alg("PartitionedSolversAlgorithm", "krylov", krylov = match.arg(solver, c("cg", "gmres", "minres", "bicgstab")),
           solver = solver, ...)
}

# The operator handed to the solver: SciMLOperators as matrices.
.lin_prepare_A <- function(A, opts, cache) {
  if (inherits(A, "AbstractSciMLOperator")) A <- concretize(A)
  A
}

# Warm start of the Krylov methods (_krylov_warm_start!): the initial guess
# and tolerances keeping the cold start threshold.
.lin_warm_start <- function(mode, A, b, u, opts) {
  if (is.null(mode) || !mode %in% c("Previous", "Hegedus") || !is.null(opts$u0)) return(opts)
  if (is.null(u) || length(u) != (if (is.function(A)) length(b) else ncol(A))) return(opts)
  unorm <- sqrt(sum(u * u))
  if (unorm == 0 || !is.finite(unorm)) return(opts)
  mv <- if (is.function(A)) A else function(x) as.numeric(A %*% x)
  if (mode == "Hegedus") {
    Au <- mv(u)
    d <- sum(Au * Au)
    if (d == 0 || !is.finite(d)) return(opts)
    Aub <- sum(Au * b)
    if (!is.finite(Aub)) return(opts)
    bnorm <- sqrt(sum(b * b))
    if (bnorm == 0 || !is.finite(bnorm)) return(opts)
    if (abs(Aub) < sqrt(1 - 0.5^2) * sqrt(d) * bnorm) return(opts)
    opts$u0 <- (Aub / d) * u
  } else opts$u0 <- u
  # the Krylov methods of R measure reltol against |b|, so the stopping
  # threshold is the one of a cold start
  opts
}

#' @rdname linear_solve_extras
#' @export
update_gamma <- function(cache, gamma) {
  A <- cache$A
  if (!inherits(A, "WOperator"))
    stop("update_gamma! needs the cache's `A` to be a `WOperator`, got ", class(A)[1L], call. = FALSE)
  cache$A <- WOperator(A$mass_matrix, gamma, A$J)
  cache$fact <- NULL; cache$fact_A <- NULL
  invisible(cache)
}

.panel_check <- function(W, np, operation) {
  if (!operation %in% c("factor_right_upper", "factor_lower", "lower", "upper"))
    stop("unknown supernodal panel operation: ", operation, call. = FALSE)
  np <- as.integer(np)
  if (np < 0 || np > min(dim(W))) stop("np must be at most the size of W", call. = FALSE)
  np
}

#' @rdname linear_solve_extras
#' @export
supernodal_panel_solve <- function(W, B, np, operation, algorithm = "auto") {
  W <- as.matrix(W); B <- as.matrix(B)
  np <- .panel_check(W, np, operation)
  if (algorithm == "auto") {
    algorithm <- if (operation %in% c("factor_right_upper", "factor_lower")) "triangularsolve"
                 else if (np <= 256L || ncol(B) == 1L) "kernel"
                 else if (np > 1792L) "blas" else "triangularsolve"
  }
  supernodal_panel_solve_backend(algorithm, W, B, np, operation)
}

#' @rdname linear_solve_extras
#' @export
supernodal_panel_solve_backend <- function(algorithm, W, B, np, operation) {
  W <- as.matrix(W); B <- as.matrix(B)
  np <- .panel_check(W, np, operation)
  algorithm <- match.arg(algorithm, c("kernel", "blas", "triangularsolve"))
  if (operation %in% c("factor_right_upper", "factor_lower") || algorithm == "triangularsolve") {
    if (operation == "factor_right_upper") {
      U <- W[seq_len(np), seq_len(np), drop = FALSE]
      U[lower.tri(U)] <- 0
      return(t(forwardsolve(t(U), t(B))))
    }
    if (operation == "factor_lower") operation <- "lower"
    algorithm <- "blas"
  }
  rows <- seq_len(np)
  if (algorithm == "blas") {
    D <- W[rows, rows, drop = FALSE]
    if (operation == "lower") {
      D[upper.tri(D)] <- 0; diag(D) <- 1
      B[rows, ] <- forwardsolve(D, B[rows, , drop = FALSE])
    } else B[rows, ] <- backsolve(D, B[rows, , drop = FALSE])
    return(B)
  }
  # the column oriented kernels
  for (r in seq_len(ncol(B))) {
    x <- B[, r]
    if (operation == "lower") {
      for (j in rows) {
        xj <- x[j]
        if (xj == 0 || j == np) next
        i <- (j + 1L):np
        x[i] <- x[i] - W[i, j] * xj
      }
    } else {
      for (j in rev(rows)) {
        xj <- x[j] / W[j, j]
        x[j] <- xj
        if (xj == 0 || j == 1L) next
        i <- seq_len(j - 1L)
        x[i] <- x[i] - W[i, j] * xj
      }
    }
    B[, r] <- x
  }
  B
}

#' @rdname linear_solve_extras
#' @export
slate_isavailable <- function(libpath = NULL) TRUE

#' @rdname linear_solve_extras
#' @export
show_algorithm_choices <- function() {
  line <- function(ch, n) cat(strrep(ch, n), "\n", sep = "")
  line("=", 60)
  cat("LinearSolve.jl Algorithm Choice Analysis\n")
  line("=", 60)
  cat("Current Preferences:\n")
  line("-", 60)
  cat("No autotune preferences currently set.\n")
  cat("\nDefault Algorithm Choices:\n")
  line("-", 80)
  cat("Size       Category    Float64            ComplexF64\n")
  line("-", 80)
  for (tc in list(list(8, "tiny"), list(50, "small"), list(200, "medium"), list(500, "large"), list(1500, "big"))) {
    n <- tc[[1L]]
    A <- diag(n)
    alg <- .lin_default(A, OperatorAssumptions(TRUE))$name
    cat(formatC(sprintf("%dx%d", n, n), width = 10), " ", formatC(tc[[2L]], width = -11), " ",
        formatC(alg, width = -18), " ", formatC(alg, width = -18), "\n", sep = "")
  }
  cat("\nSystem Information:\n")
  line("-", 60)
  cat("LAPACK: ", La_library(), "\n", sep = "")
  cat("BLAS: ", extSoftVersion()[["BLAS"]], "\n", sep = "")
  line("=", 60)
  invisible(NULL)
}
