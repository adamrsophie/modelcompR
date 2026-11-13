# Geodesic disciplined convex programming (gDCP): the geodesic curvature of
# SymbolicAnalysis.jl on the manifold of symmetric positive definite matrices
# and on the Lorentz model, with the matrix atoms it certifies.
#
# Ported from SymbolicAnalysis.jl 0.5.8 (src/gdcp/gdcp_rules.jl, spd.jl,
# lorentz.jl, canon.jl and fold.jl;
# https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl), Copyright (c) 2024
# Vaibhav Dixit and contributors, MIT License (inst/CONVEX-LICENSES).
#
# A matrix variable is a symcompR MatrixSymbol (MatrixSymbol("X", n, n)); a
# Lorentz point is an n x 1 MatrixSymbol. An atom applied to one returns an
# uninterpreted function node named as the Julia atom. Constant matrix and
# vector arguments are stored as their entries (by columns), function
# arguments (affine_map(conjugation, ...)) and the manifold of distance as
# symbols named as they print; on numbers every atom computes its value.

# ---------------------------------------------------------------------------
# Manifolds

#' @rdname manifolds
#' @export
SymmetricPositiveDefinite <- function(n)
  structure(list(n = .mf_check_dim(n)), class = c("SymmetricPositiveDefinite", "AbstractManifold"))
#' @rdname manifolds
#' @export
Lorentz <- function(n) structure(list(n = .mf_check_dim(n)), class = c("Lorentz", "AbstractManifold"))

# ---------------------------------------------------------------------------
# Node structure

# symbols that stand for a function or manifold argument, not a variable
.sa_state$tags <- character()
.sg_tag <- function(name) {
  .sa_state$tags <- union(.sa_state$tags, name)
  Symbol(name)
}

.sg_leaf_types <- c("Symbol", "MatrixSymbol")
.sg_is_var <- function(e) is_sym(e) && .type(e) %in% .sg_leaf_types &&
  !(.type(e) == "Symbol" && .name(e) %in% .sa_state$tags)
# a constant holds no variable (Julia's isconstarg)
.sg_has_var <- function(e) {
  if (!is_sym(e)) return(FALSE)
  if (.type(e) %in% .sg_leaf_types) return(.sg_is_var(e))
  any(vapply(.args(e), .sg_has_var, TRUE))
}
.sg_const <- function(e) !.sg_has_var(e)
# issym(e) || iscall(e) of Julia: anything that is not a constant
.sg_iscall <- function(e) .sg_has_var(e) && !(.type(e) %in% .sg_leaf_types)

# the value of a constant factor: a number, or 1 for an identity matrix
.sg_constval <- function(e) {
  if (!is_sym(e)) return(as.numeric(e))
  if (.cvx_is_num(e)) return(.cvx_num(e))
  if (.type(e) == "Identity") return(1)
  if (.type(e) %in% c("Mul", "MatMul")) return(prod(vapply(.args(e), .sg_constval, 0)))
  NA_real_
}

.sg_matrix_types <- c("MatrixSymbol", "MatMul", "MatAdd", "Inverse", "Transpose", "Adjoint", "MatPow", "Identity",
                      "HadamardProduct")
.sg_matrix_fns <- c("conjugation", "affine_map", "hadamard_product", "scalar_mat", "adjoint", "inv", "diag")
# whether an expression is array valued
.sg_is_mat <- function(e) {
  if (!is_sym(e)) return(FALSE)
  t <- .type(e)
  if (t %in% .sg_matrix_types) return(TRUE)
  if (t == "Pow") return(.sg_is_mat(.args(e)[[1L]]))
  if (t == "Function") {
    nm <- .name(e)
    if (nm %in% .sg_matrix_fns) return(TRUE)
    if (nm %in% c("log", "exp", "sqrt")) return(.sg_is_mat(.args(e)[[1L]]))
  }
  FALSE
}

# the Julia operation of a node
.sg_op <- function(e) {
  if (!is_sym(e)) return("")
  t <- .type(e)
  switch(t,
    Add = , MatAdd = "+",
    Mul = , MatMul = "*",
    # sqrt(X) is a Pow, or a MatPow for a matrix symbol in newer symcompR
    Pow = , MatPow = {
      i <- .args(e)[[2L]]
      if (.cvx_is_num(i) && .cvx_num(i) == 0.5) "sqrt" else "^"
    },
    Inverse = "inv", Transpose = , Adjoint = "adjoint", Trace = "tr",
    Function = switch(.name(e), Abs = "abs", Max = "max", Min = "min", .name(e)),
    t)
}

# ---------------------------------------------------------------------------
# Rule table: add_gdcprule

.sg_rules <- new.env(parent = emptyenv())

.SG_CURV <- c("GConvex", "GConcave", "GLinear", "GUnknownCurvature")
.SG_MONO <- c("GIncreasing", "GDecreasing", "GAnyMono")

#' Register a geodesic DCP rule
#'
#' `add_gdcprule(f, manifold, sign, curvature, monotonicity)` registers the
#' geodesic curvature of an atom for [analyze()] on a manifold, as
#' SymbolicAnalysis.jl. The atom is a function node named `f` (build it with
#' `Function(f)`). As upstream the table is keyed by the atom alone: the
#' manifold is recorded but not consulted, so a rule applies on either
#' supported manifold, and a later registration replaces an earlier one.
#'
#' @param f The atom name.
#' @param manifold The manifold constructor or name the rule is stated for,
#'   such as `"Lorentz"`.
#' @param sign `"Positive"`, `"Negative"` or `"AnySign"`.
#' @param curvature `"GConvex"`, `"GConcave"`, `"GLinear"` or
#'   `"GUnknownCurvature"`.
#' @param monotonicity `"GIncreasing"`, `"GDecreasing"` or `"GAnyMono"`, one
#'   per argument or one for all.
#' @return The rule, invisibly.
#' @examples
#' add_gdcprule("first_coordinate", "Lorentz", "Positive", "GConvex", "GIncreasing")
#' p <- MatrixSymbol("p", 3, 1)
#' analyze(Function("first_coordinate")(p), Lorentz(3))$gcurvature
#' @export
add_gdcprule <- function(f, manifold, sign, curvature, monotonicity) {
  if (is.function(manifold)) manifold <- deparse(substitute(manifold))
  rule <- list(manifold = manifold, sign = sign, gcurvature = curvature, gmonotonicity = monotonicity)
  assign(f, rule, envir = .sg_rules)
  invisible(rule)
}

.sg_rule <- function(f) if (nzchar(f) && exists(f, envir = .sg_rules, inherits = FALSE)) get(f, envir = .sg_rules) else NULL

local({
  spd <- "SymmetricPositiveDefinite"; lor <- "Lorentz"
  add_gdcprule("logdet", spd, "AnySign", "GLinear", "GIncreasing")
  add_gdcprule("conjugation", spd, "Positive", "GConvex", "GIncreasing")
  add_gdcprule("tr", spd, "Positive", "GConvex", "GIncreasing")
  add_gdcprule("sum", spd, "Positive", "GConvex", "GIncreasing")
  add_gdcprule("adjoint", spd, "Positive", "GLinear", "GIncreasing")
  add_gdcprule("scalar_mat", spd, "Positive", "GConvex", "GIncreasing")
  add_gdcprule("diag", spd, "Positive", "GConvex", "GIncreasing")
  add_gdcprule("sdivergence", spd, "Positive", "GConvex", "GIncreasing")
  add_gdcprule("distance", spd, "Positive", "GConvex", "GAnyMono")
  add_gdcprule("quad_form", spd, "Positive", "GConvex", "GIncreasing")
  add_gdcprule("eigmax", spd, "Positive", "GConvex", "GIncreasing")
  add_gdcprule("log_quad_form", spd, "AnySign", "GConvex", "GIncreasing")
  add_gdcprule("inv", spd, "Positive", "GConvex", "GDecreasing")
  add_gdcprule("eigsummax", spd, "Positive", "GConvex", "GIncreasing")
  add_gdcprule("schatten_norm", spd, "Positive", "GConvex", "GIncreasing")
  add_gdcprule("sum_log_eigmax", spd, "AnySign", "GConvex", "GIncreasing")
  add_gdcprule("affine_map", spd, "Positive", "GConvex", "GIncreasing")
  add_gdcprule("hadamard_product", spd, "Positive", "GConvex", "GIncreasing")
  # lorentz.jl registers distance again, over the SPD rule
  add_gdcprule("distance", lor, "Positive", "GConvex", "GAnyMono")
  add_gdcprule("lorentz_log_barrier", lor, "AnySign", "GConvex", "GIncreasing")
  add_gdcprule("lorentz_homogeneous_quadratic", lor, "AnySign", "GConvex", "GAnyMono")
  add_gdcprule("lorentz_homogeneous_diagonal", lor, "AnySign", "GConvex", "GAnyMono")
  # these two carry the Euclidean AnyMono upstream
  add_gdcprule("lorentz_nonhomogeneous_quadratic", lor, "AnySign", "GConvex", "AnyMono")
  add_gdcprule("lorentz_least_squares", lor, "Positive", "GConvex", "AnyMono")
})

# ---------------------------------------------------------------------------
# Euclidean rules of the matrix atoms (atoms.jl); NULL when there is none

.sg_euclid_rule <- function(e) {
  args <- .args(e)
  one <- function(s, c, m) list(s, c, m, lapply(args, list))
  switch(.sg_op(e),
    tr = , sum = , diag = , adjoint = one("AnySign", "Affine", "Increasing"),
    inv = one("AnySign", "Convex", "Decreasing"),
    sqrt = one("Positive", "Concave", "Increasing"),
    `^` = one("AnySign", "UnknownCurvature", "AnyMono"),
    eigmax = , eigsummax = one("AnySign", "Convex", "AnyMono"),
    eigmin = , eigsummin = one("AnySign", "Concave", "AnyMono"),
    opnorm = one("Positive", "Convex", "AnyMono"),
    # LinearAlgebra.norm of a matrix: the norm of its entries
    norm = one("Positive", "Convex", .sa_incr_if_pos),
    logdet = one("AnySign", "Concave", "AnyMono"),
    # a quad_form whose matrix is not a constant
    quad_form = one("AnySign", "UnknownCurvature", "AnyMono"),
    NULL)
}

# whether a node is one of the matrix forms above rather than a traced atom
.sg_matrix_node <- function(e) {
  t <- .type(e)
  if (t %in% c("Inverse", "Transpose", "Adjoint", "Trace", "MatPow")) return(TRUE)
  if (t == "Pow") return(.sg_is_mat(.args(e)[[1L]]))
  t == "Function" && length(.args(e)) > 0L && (.sg_is_mat(.args(e)[[1L]]) ||
    (.name(e) %in% c("tr", "sum", "diag", "eigsummax", "eigsummin") && .sg_has_var(.args(e)[[1L]])))
}

# ---------------------------------------------------------------------------
# Geodesic curvature (find_gcurvature, node_gcurvature and the fold)

.sg_flip <- function(g) switch(g, GConvex = "GConcave", GConcave = "GConvex", g)

.sg_gcurv <- function(e) {
  if (!is_sym(e) || .sg_const(e)) return("GLinear")
  if (.type(e) %in% .sg_leaf_types) return("GLinear")
  k <- .key(e)
  hit <- .sa_state$gmemo[[k]]
  if (!is.null(hit)) return(hit)
  op <- .sg_op(e)
  g <- if (op == "*") .sg_mul_gcurv(.args(e)) else if (op == "+") .sg_add_gcurv(.args(e)) else .sg_find_gcurv(e)
  assign(k, g, envir = .sa_state$gmemo)
  g
}

.sg_mul_gcurv <- function(args) {
  nonconst <- NULL; count <- 0L; prod <- 1
  for (a in args) {
    if (.sg_has_var(a)) {
      count <- count + 1L
      nonconst <- a
      if (count > 1L) {
        warning("DGCP does not support multiple non-constant arguments in multiplication", call. = FALSE)
        return("GUnknownCurvature")
      }
    } else prod <- prod * .sg_constval(a)
  }
  if (is.null(nonconst)) return("GLinear")
  g <- .sg_gcurv(nonconst)
  if (!is.na(prod) && prod < 0) .sg_flip(g) else g
}

.sg_add_gcurv <- function(args) {
  cvx <- FALSE; ccv <- FALSE
  for (a in args) {
    g <- .sg_gcurv(a)
    if (g == "GLinear") next
    if (g == "GConvex") {
      cvx <- TRUE
      if (ccv) return("GUnknownCurvature")
    } else if (g == "GConcave") {
      ccv <- TRUE
      if (cvx) return("GUnknownCurvature")
    } else return("GUnknownCurvature")
  }
  if (cvx) "GConvex" else if (ccv) "GConcave" else "GLinear"
}

# Julia indexes the monotonicity of a concave rule directly: a rule built by
# add_dcprule holds a tuple (an index past its end is a BoundsError), one
# built by power_rule or the norm rule a bare value (a MethodError)
.sg_mono_at <- function(mono, i, op) {
  if (is.function(mono) || op %in% c("^", "norm"))
    stop("MethodError: no method matching getindex(::", if (is.function(mono)) "typeof(increasing_if_positive)"
         else "SymbolicAnalysis.Monotonicity", ", ::Int64) The function `getindex` exists, but no method is defined ",
         "for this combination of argument types.", call. = FALSE)
  if (i > length(mono))
    stop("BoundsError: attempt to access Tuple{SymbolicAnalysis.Monotonicity} at index [", i, "]", call. = FALSE)
  mono[[i]]
}

.sg_mono_prop <- function(mono, i, arg) {
  m <- if (is.function(mono)) mono else mono[[min(i, length(mono))]]
  if (is.function(m)) m(list(arg)) else m
}

.sg_find_gcurv <- function(e) {
  f <- .sg_op(e)
  args <- .args(e)
  calls <- vapply(args, .sg_iscall, TRUE)
  known <- FALSE
  fc <- "GUnknownCurvature"; fm <- "GAnyMono"
  rule <- .sg_rule(f)
  a1 <- if (length(args)) args[[1L]] else NULL
  if (!is.null(rule) && !any(calls)) {
    fc <- rule$gcurvature; fm <- rule$gmonotonicity; known <- TRUE
  } else if (f == "logdet") {
    if (.sg_op(a1) %in% c("conjugation", "diag", "+", "affine_map", "hadamard_product")) return("GConvex")
  } else if (f == "log" && .sg_iscall(a1) && .sg_op(a1) %in% c("tr", "quad_form")) {
    return("GConvex")
  } else if (f %in% c("schatten_norm", "eigsummax") && .sg_op(a1) == "log") {
    return("GConvex")
  } else if (f == "affine_map") {
    # upstream compares the stored function argument with tr, conjugation
    # and diag, but the term holds it wrapped as a constant node, so the
    # test never holds; its branch for sum_log_eigmax (hasdcprule of that
    # wrapped argument) is likewise never taken
    return("GUnknownCurvature")
  } else if (!is.null(rule) && any(calls)) {
    for (i in which(calls)) {
      o <- .sg_op(args[[i]])
      if (o == "inv") {
        # the flip compares the tuple of monotonicities with one value, so
        # it never matches and gives GAnyMono
        fc <- rule$gcurvature; fm <- "GAnyMono"; known <- TRUE
      } else if (o %in% c("broadcast", "affine_map") ||
                 (o == "+" && sum(vapply(.args(args[[i]]), .sg_has_var, TRUE)) <= 1L)) {
        fc <- rule$gcurvature; fm <- rule$gmonotonicity; known <- TRUE
      }
    }
  } else if (f == "*") {
    a <- .sg_constval(a1)
    if (!is.na(a) && a > 0) return(.sg_gcurv(args[[2L]]))
    if (!is.na(a) && a < 0) return(.sg_flip(.sg_gcurv(args[[2L]])))
    warning("Disciplined Programming does not support multiple non-constant arguments in multiplication",
            call. = FALSE)
    return("GUnknownCurvature")
  }
  if (!known) {
    r <- .sa_rule(e)
    if (is.null(r) || !any(vapply(args, function(a) .sg_gcurv(a) %in% c("GConvex", "GConcave"), TRUE)))
      return("GUnknownCurvature")
    fc <- r[[2L]]; fm <- r[[3L]]
  }
  if (fc %in% c("Convex", "Affine")) {
    for (i in seq_along(args)) {
      ac <- .sg_gcurv(args[[i]])
      m <- .sg_mono_prop(fm, i, args[[i]])
      ok <- if (ac == "GConvex") m == "Increasing" else if (ac == "GConcave") m == "Decreasing" else ac == "GLinear"
      if (!ok) return("GUnknownCurvature")
    }
    return("GConvex")
  }
  if (fc == "Concave") {
    for (i in seq_along(args)) {
      ac <- .sg_gcurv(args[[i]])
      m <- .sg_mono_at(fm, i, f)
      ok <- if (ac == "GConcave") m == "Increasing" else if (ac == "GConvex") m == "Decreasing" else ac == "GLinear"
      if (!ok) return("GUnknownCurvature")
    }
    return("GConcave")
  }
  if (fc %in% .SG_CURV) return(fc)
  "GUnknownCurvature"
}

# The canonical forms of canon.jl: a' Y a becomes quad_form(a, Y). Upstream
# lists a second rule for B' X B to conjugation(X, B) with the same pattern,
# so the first rule always wins.
.sg_canonize <- function(e) {
  if (!is_sym(e)) return(e)
  a <- .args(e)
  if (length(a) && !(.type(e) %in% .sg_leaf_types)) {
    na <- lapply(a, .sg_canonize)
    if (!identical(.keys(na), .keys(a))) e <- .rebuild(e, na)
  }
  if (.type(e) == "MatMul") {
    a <- .args(e)
    if (length(a) == 3L && .type(a[[1L]]) %in% c("Transpose", "Adjoint") &&
        identical(.key(.args(a[[1L]])[[1L]]), .key(a[[3L]])))
      return(Function("quad_form")(a[[3L]], a[[2L]]))
  }
  e
}

# ---------------------------------------------------------------------------
# Atoms

# a node of atom `name`: symbolic arguments as they are, a function or a
# manifold as a symbol, a numeric array as its entries
.sg_node <- function(name, ...) {
  parts <- list(...)
  args <- list()
  for (x in parts) {
    if (is_sym(x)) args <- c(args, list(x))
    else if (inherits(x, "AbstractManifold")) args <- c(args, list(.sg_tag(paste0(class(x)[1L], "(", x$n, ")"))))
    else if (is.character(x)) args <- c(args, list(.sg_tag(x)))
    else args <- c(args, as.list(as.numeric(unlist(x))))
  }
  do.call(Function(name), args)
}
.sg_symbolic <- function(...) any(vapply(list(...), is_sym, TRUE))

.sg_isposdef <- function(B) {
  if (!isSymmetric(unname(B), tol = 0)) return(FALSE)
  !inherits(tryCatch(chol(B), error = function(e) e), "error")
}
# Symmetric(B): the upper triangle
.sg_symmetric <- function(B) {
  B[lower.tri(B)] <- t(B)[lower.tri(B)]
  B
}
.sg_eig_sym <- function(B) sort(eigen(.sg_symmetric(B), symmetric = TRUE, only.values = TRUE)$values)
.sg_check_psd <- function(B) {
  if (!.sg_isposdef(B) || !(.sg_eig_sym(B)[1L] >= 0))
    stop("DomainError: B must be positive semi-definite.", call. = FALSE)
}

#' Geodesically convex matrix atoms
#'
#' The atoms of the geodesic analysis of SymbolicAnalysis.jl ([analyze()]
#' with a manifold). On numbers they compute their value; on a matrix
#' variable (a symcompR `MatrixSymbol`) they return the atom for analysis.
#' On the manifold of symmetric positive definite matrices: `tr(X)` (the
#' trace), `conjugation(X, B)` (\eqn{B^T X B}), `scalar_mat(X, k)`
#' (\eqn{tr(X) I_k}), `sdivergence(X, Y)` (the symmetric Stein divergence
#' \eqn{logdet((X + Y)/2) - logdet(X Y)/2}), `distance(M, p, q)` (the
#' Riemannian distance on `M`), `log_quad_form(y, X)` (\eqn{log(y^T X y)}, or
#' the log of the sum over a list of vectors), `eigsummax(m, k)` and
#' `eigsummin(m, k)` (the sum of the `k` largest or smallest eigenvalues),
#' `schatten_norm(X, p)` (the `p`-norm of the eigenvalues),
#' `sum_log_eigmax(f, X, k)` (the sum of `f` of the logs of the `k` largest
#' eigenvalues; `f` defaults to the identity), `affine_map(f, X, B, Y)`
#' (\eqn{B} plus `f` of `X`: `conjugation` with a matrix or a list of
#' matrices `Y`, `tr` or `diag`, or `hadamard_product` with `Y`) and
#' `hadamard_product(X, B)` (the entrywise product). On the Lorentz model:
#' `lorentz_log_barrier(p)`, `lorentz_homogeneous_quadratic(A, p)`,
#' `lorentz_homogeneous_diagonal(a, p)`,
#' `lorentz_nonhomogeneous_quadratic(A, b, c, p)`,
#' `lorentz_least_squares(X, y, p)` and `lorentz_transform(O, p)`, which
#' check the conditions for geodesic convexity on numbers. [logdet()],
#' [eigmax()], [eigmin()], [opnorm()] and [quad_form()] also take a matrix
#' variable, and the inverse is symcompR's `Inverse(X)`. The sum and the
#' diagonal of a matrix variable are `Function("sum")(X)` and
#' `Function("diag")(X)`.
#'
#' @param X,Y,m,B,A,O Numeric matrices or matrix variables.
#' @param M A manifold: [SymmetricPositiveDefinite()] or [Lorentz()]. For
#'   two points of symcompR geometry, `distance(p, q)` is symcompR's.
#' @param p,q Points: matrices on the positive definite manifold, vectors
#'   on the Lorentz model.
#' @param k An integer.
#' @param y A vector, or a list of vectors.
#' @param f A function name: `"conjugation"`, `"tr"`, `"diag"` or
#'   `"hadamard_product"` for `affine_map`; an elementwise function name for
#'   `sum_log_eigmax`.
#' @param a,b Numeric vectors.
#' @param c A number.
#' @param ... The arguments of `sum_log_eigmax`: `(f, X, k)` or `(X, k)`.
#' @return A number or matrix, or the atom.
#' @examples
#' A <- matrix(c(2, 0.5, 0.5, 1), 2)
#' conjugation(diag(2), A)
#' sdivergence(A, diag(2))
#' X <- MatrixSymbol("X", 2, 2)
#' analyze(logdet(conjugation(X, A)) - logdet(X), SymmetricPositiveDefinite(2))$gcurvature
#' @name gdcp_atoms
NULL

#' @rdname gdcp_atoms
#' @export
tr <- function(X) if (is_sym(X)) .sg_node("tr", X) else sum(diag(as.matrix(X)))

#' @rdname gdcp_atoms
#' @export
conjugation <- function(X, B) {
  if (.sg_symbolic(X, B)) return(.sg_node("conjugation", X, B))
  t(B) %*% X %*% B
}

#' @rdname gdcp_atoms
#' @export
scalar_mat <- function(X, k = nrow(X)) {
  if (is_sym(X)) return(.sg_node("scalar_mat", X, k))
  tr(X) * diag(k)
}

.sg_logdet_num <- function(M) {
  d <- determinant(M, logarithm = TRUE)
  if (d$sign < 0) stop("DomainError: logdet of a matrix with a negative determinant", call. = FALSE)
  as.numeric(d$modulus)
}

#' @rdname gdcp_atoms
#' @export
sdivergence <- function(X, Y) {
  if (.sg_symbolic(X, Y)) return(.sg_node("sdivergence", X, Y))
  .sg_logdet_num((X + Y) / 2) - 1 / 2 * .sg_logdet_num(X %*% Y)
}

#' @rdname gdcp_atoms
#' @export
distance <- function(M, p, q) {
  if (!inherits(M, "AbstractManifold")) return(symcompR::distance(M, p))
  if (.sg_symbolic(p, q)) return(.sg_node("distance", M, p, q))
  if (inherits(M, "SymmetricPositiveDefinite")) {
    # the norm of the logarithmic map of the affine invariant metric
    e <- eigen(.sg_symmetric(p), symmetric = TRUE)
    pinv <- e$vectors %*% diag(1 / sqrt(e$values), length(e$values)) %*% t(e$vectors)
    s <- eigen(.sg_symmetric(pinv %*% q %*% pinv), symmetric = TRUE, only.values = TRUE)$values
    return(sqrt(sum(log(s)^2)))
  }
  if (inherits(M, "Lorentz")) {
    # the Minkowski norm of q - p, the logarithmic map of the flat metric
    w <- as.numeric(q) - as.numeric(p)
    n <- length(w)
    ip <- -w[n] * w[n] + sum(w[-n] * w[-n])
    if (ip < 0) stop("DomainError with ", .jl_show_float(ip), ": sqrt was called with a negative real argument",
                     call. = FALSE)
    return(sqrt(ip))
  }
  stop("`distance` is defined on SymmetricPositiveDefinite and Lorentz manifolds", call. = FALSE)
}

#' @rdname gdcp_atoms
#' @export
log_quad_form <- function(y, X) {
  if (.sg_symbolic(X) || (!is.list(y) && .sg_symbolic(y))) return(.sg_node("log_quad_form", y, X))
  if (is.list(y)) return(log(sum(vapply(y, function(v) sum(v * (X %*% v)), 0))))
  log(sum(y * (X %*% y)))
}

.sg_check_k <- function(m, k) if (k < 1 || k > nrow(m)) stop("DomainError with ", k,
                                                            ": k must be between 1 and size(m, 1)", call. = FALSE)

#' @rdname gdcp_atoms
#' @export
eigsummax <- function(m, k) {
  if (is_sym(m)) return(.sg_node("eigsummax", m, k))
  .sg_check_k(m, k)
  ev <- .sg_eig_sym(m)
  sum(ev[(nrow(m) - k + 1L):nrow(m)])
}

#' @rdname gdcp_atoms
#' @export
eigsummin <- function(m, k) {
  if (is_sym(m)) return(.sg_node("eigsummin", m, k))
  .sg_check_k(m, k)
  sum(.sg_eig_sym(m)[seq_len(k)])
}

#' @rdname gdcp_atoms
#' @export
schatten_norm <- function(X, p = 2L) {
  if (is_sym(X)) return(.sg_node("schatten_norm", X, p))
  ev <- eigen(X, only.values = TRUE)$values
  if (is.infinite(p)) max(Mod(ev)) else sum(Mod(ev)^p)^(1 / p)
}

#' @rdname gdcp_atoms
#' @export
sum_log_eigmax <- function(...) {
  a <- list(...)
  if (length(a) == 3L) { f <- a[[1L]]; X <- a[[2L]]; k <- a[[3L]] } else { f <- NULL; X <- a[[1L]]; k <- a[[2L]] }
  if (is_sym(X)) return(if (is.null(f)) .sg_node("sum_log_eigmax", X, k) else .sg_node("sum_log_eigmax", f, X, k))
  n <- nrow(X)
  ev <- .sg_eig_sym(X)[(n - k + 1L):n]
  if (is.null(f)) sum(log(ev)) else sum(match.fun(f)(log(ev)))
}

#' @rdname gdcp_atoms
#' @export
affine_map <- function(f, X, B, Y) {
  if (!is.character(f) || length(f) != 1L)
    stop("`f` must be one of \"conjugation\", \"tr\", \"diag\" or \"hadamard_product\"", call. = FALSE)
  if (is_sym(X)) return(if (missing(Y)) .sg_node("affine_map", f, X, B) else .sg_node("affine_map", f, X, B, Y))
  if (f == "hadamard_product") {
    # affine_map(hadamard_product, X, Y, B): the offset is the last argument
    .sg_check_psd(Y)
    return(Y + hadamard_product(X, B))
  }
  .sg_check_psd(B)
  switch(f,
    conjugation = if (is.list(Y)) B + Reduce(`+`, lapply(Y, function(y) conjugation(X, y))) else B + conjugation(X, Y),
    tr = stop("MethodError: no method matching +(::Matrix{Float64}, ::Float64)", call. = FALSE),
    diag = stop("DimensionMismatch: dimensions must match: a has dims (Base.OneTo(", nrow(B), "), Base.OneTo(",
                ncol(B), ")), b has dims (Base.OneTo(", nrow(X), "),), mismatch at 2", call. = FALSE),
    stop("`f` must be one of \"conjugation\", \"tr\", \"diag\" or \"hadamard_product\"", call. = FALSE))
}

#' @rdname gdcp_atoms
#' @export
hadamard_product <- function(X, B) {
  if (is_sym(X) && .type(X) == "MatrixSymbol" && !is_sym(B)) return(.sg_node("hadamard_product", X, B))
  if (is_sym(X) || is_sym(B) || inherits(X, "symmatrix") || inherits(B, "symmatrix"))
    return(symcompR::hadamard_product(X, B))
  if ((!.sg_isposdef(B) || !(.sg_eig_sym(B)[1L] >= 0)) && !any(apply(B, 1L, prod) == 0))
    stop("DomainError: B must be positive semi-definite and have no zero rows.", call. = FALSE)
  B * X
}

#' @rdname gdcp_atoms
#' @export
lorentz_log_barrier <- function(p) {
  if (is_sym(p)) return(.sg_node("lorentz_log_barrier", p))
  -log(-1 + p[length(p)])
}

#' @rdname gdcp_atoms
#' @export
lorentz_homogeneous_quadratic <- function(A, p) {
  if (is_sym(p)) return(.sg_node("lorentz_homogeneous_quadratic", A, p))
  d <- nrow(A) - 1L
  a_vec <- A[seq_len(d), d + 1L]
  sigma <- A[d + 1L, d + 1L]
  lambda_min <- min(eigen(A[seq_len(d), seq_len(d), drop = FALSE], only.values = TRUE)$values)
  c1 <- sqrt(sum(a_vec^2)) <= 1e-10 && sigma >= -lambda_min
  c2 <- sigma + lambda_min > 2 * sqrt(sum(a_vec * a_vec))
  if (!(c1 || c2)) stop("ArgumentError: Matrix A does not satisfy geodesic convexity conditions", call. = FALSE)
  sum(p * (A %*% p))
}

#' @rdname gdcp_atoms
#' @export
lorentz_homogeneous_diagonal <- function(a, p) {
  if (is_sym(p)) return(.sg_node("lorentz_homogeneous_diagonal", a, p))
  if (length(a) != length(p)) stop("DimensionMismatch: Vectors must have same length", call. = FALSE)
  if (min(a[-length(a)]) + a[length(a)] < 0)
    stop("ArgumentError: For geodesic convexity, min(a[1:end-1]) + a[end] \u2265 0 is required", call. = FALSE)
  sum(a * p^2)
}

#' @rdname gdcp_atoms
#' @export
lorentz_nonhomogeneous_quadratic <- function(A, b, c, p) {
  if (is_sym(p)) return(.sg_node("lorentz_nonhomogeneous_quadratic", A, b, c, p))
  n <- length(b)
  if (!(sum(b[-n]^2) <= b[n]^2 && b[n] >= 0))
    stop("ArgumentError: Vector b must be in the Lorentz cone for geodesic convexity", call. = FALSE)
  lorentz_homogeneous_quadratic(A, p) + sum(b * p) + c
}

#' @rdname gdcp_atoms
#' @export
lorentz_least_squares <- function(X, y, p) {
  if (is_sym(p)) return(.sg_node("lorentz_least_squares", X, y, p))
  lorentz_nonhomogeneous_quadratic(t(X) %*% X, as.numeric(-2 * t(X) %*% y), sum(y * y), p)
}

#' @rdname gdcp_atoms
#' @export
lorentz_transform <- function(O, p) {
  if (is_sym(p)) return(.sg_node("lorentz_transform", O, p))
  d <- length(p) - 1L
  J <- diag(c(rep(1, d), -1))
  G <- t(O) %*% J %*% O
  if (!(sqrt(sum((G - J)^2)) <= 1e-10 * max(sqrt(sum(G^2)), sqrt(sum(J^2)))))
    stop("ArgumentError: Matrix is not in the Lorentz group", call. = FALSE)
  if ((O %*% c(rep(0, d), 1))[d + 1L] <= 0)
    stop("ArgumentError: Matrix does not preserve the positive time direction", call. = FALSE)
  as.numeric(O %*% p)
}
