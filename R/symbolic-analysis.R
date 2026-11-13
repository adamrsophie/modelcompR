# Disciplined convex programming analysis of symbolic expressions (the
# curvature and sign of SymbolicAnalysis.jl) and the linear expansion of
# Symbolics.jl, with the matrix products and norm atoms that trace a convex
# problem.
#
# Ported from SymbolicAnalysis.jl 0.3.9
# (https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl), Copyright (c) 2024
# Vaibhav Dixit and contributors, MIT License (inst/CONVEX-LICENSES), and
# from linear_expansion of Symbolics.jl
# (https://github.com/JuliaSymbolics/Symbolics.jl), Copyright (c) 2021:
# Shashi Gowda, Yingbo Ma, Chris Rackauckas, Julia Computing, MIT License
# (inst/COPYRIGHTS).

# ---------------------------------------------------------------------------
# Tracing: the symbolic vectors of modelingtoolkitize (mtk_symvec,
# R/modelingtoolkitize-trace.R), with matrix products and norm atoms

# Matrix products and norm seen by a traced function
.cvx_matmul <- function(x, y) {
  if (!inherits(x, "mtk_symvec") && !inherits(y, "mtk_symvec")) return(base::`%*%`(x, y))
  dot <- function(cf, es) {
    out <- as_sym(0)
    for (j in seq_along(es)) if (cf[j] != 0) out <- out + cf[j] * es[[j]]
    out
  }
  if (inherits(y, "mtk_symvec")) {
    X <- if (is.matrix(x)) x else matrix(as.numeric(x), nrow = 1L)
    es <- .symvec_elems(y)
    if (ncol(X) != length(es)) stop("non-conformable arguments", call. = FALSE)
    return(.symvec(lapply(seq_len(nrow(X)), function(i) dot(X[i, ], es))))
  }
  Y <- if (is.matrix(y)) y else matrix(as.numeric(y), ncol = 1L)
  es <- .symvec_elems(x)
  if (nrow(Y) != length(es)) stop("non-conformable arguments", call. = FALSE)
  .symvec(lapply(seq_len(ncol(Y)), function(j) dot(Y[, j], es)))
}

# norm(w, type) of a symbolic vector is the atom norm(w..., order); the types
# of base::norm map to the vector norms they compute
.cvx_norm <- function(x, type = "O", ...) {
  if (!inherits(x, "mtk_symvec") && !is_sym(x)) return(base::norm(x, type, ...))
  order <- if (is.character(type)) {
    switch(toupper(substr(type, 1L, 1L)), O = , "1" = as_sym(1), I = , M = Symbol(".norm_Inf"), F = , E = , "2" = as_sym(2),
           stop("argument type[1]='", type, "' must be one of 'M','1','O','I','F' or 'E'", call. = FALSE))
  } else if (inherits(type, "mtk_symvec") || is_sym(type)) .symvec_elems(type)[[1L]]
  else if (is.numeric(type) && length(type) == 1L) {
    if (is.infinite(type) && type > 0) Symbol(".norm_Inf") else as_sym(type)
  } else stop("invalid norm type", call. = FALSE)
  .symvec(list(do.call(Function("norm"), c(.symvec_elems(x), list(order)))))
}

.cvx_traceable <- function(fun) {
  if (!is.function(fun) || is.primitive(fun)) return(fun)
  env <- new.env(parent = environment(fun))
  env$norm <- .cvx_norm
  env$`%*%` <- .cvx_matmul
  environment(fun) <- env
  fun
}

# ---------------------------------------------------------------------------
# Expression helpers

.cvx_is_num <- function(e) is_sym(e) && .type(e) %in% c("Integer", "Rational", "Float")
.cvx_num <- function(e) as.numeric(e)
.cvx_free <- function(e) if (.cvx_is_num(e)) list() else free_symbols(e)
.cvx_depends <- function(e, keys) any(.keys(.cvx_free(e)) %in% keys)
.cvx_is_norm <- function(e) .type(e) == "Function" && identical(.name(e), "norm")
.cvx_fname <- function(e) if (.type(e) == "Function") .name(e) else .type(e)

# ---------------------------------------------------------------------------
# Curvature and sign (SymbolicAnalysis.jl rules)

.SA_SIGN <- c("Positive", "Negative", "AnySign")
.SA_CURV <- c("Convex", "Concave", "Affine", "UnknownCurvature")

.sa_num_sign <- function(v) if (v < 0) "Negative" else "Positive"

.sa_incr_if_pos <- function(arg) {
  s <- .sa_sign(arg)
  if (s == "AnySign") "AnyMono" else if (s == "Positive") "Increasing" else "Decreasing"
}

# The DCP rule of a call: list(sign, curvature, monotonicity (vector, the last
# repeats; a function of the argument for increasing_if_positive))
.sa_rule <- function(e) {
  args <- .args(e)
  f <- .cvx_fname(e)
  if (.type(e) == "Pow") {
    i <- args[[2L]]
    if (!.cvx_is_num(i)) return(NULL)
    i <- .cvx_num(i)
    return(if (i == 1) list("AnySign", "Affine", "Increasing")
      else if (i == round(i) && i %% 2 == 0) list("Positive", "Convex", .sa_incr_if_pos)
      else if (i == round(i)) list("Positive", "Convex", "Increasing")
      else if (i >= 1) list("Positive", "Convex", "Increasing")
      else if (i > 0) list("Positive", "Concave", "Increasing")
      else list("Positive", "Convex", "Increasing"))
  }
  switch(f,
    exp = list("Positive", "Convex", "Increasing"),
    log = list("AnySign", "Concave", "Increasing"),
    Abs = , abs = list("Positive", "Convex", .sa_incr_if_pos),
    Max = list("AnySign", "Convex", "Increasing"),
    Min = list("AnySign", "Concave", "Increasing"),
    log1p = list("Negative", "Concave", "Increasing"),
    norm = {
      p <- args[[length(args)]]
      pv <- if (identical(.key(p), .key(Symbol(".norm_Inf")))) Inf else if (.cvx_is_num(p)) .cvx_num(p) else NA
      if (is.na(pv)) list("Positive", "UnknownCurvature", "AnyMono")
      else if (pv >= 1) list("Positive", "Convex", .sa_incr_if_pos)
      else if (pv > 0) list("Positive", "UnknownCurvature", "Increasing")
      else list("Positive", "UnknownCurvature", "AnyMono")
    },
    NULL)
}

.sa_sign <- function(e) {
  if (!is_sym(e)) return(.sa_num_sign(e))
  if (.cvx_is_num(e)) return(.sa_num_sign(.cvx_num(e)))
  t <- .type(e)
  if (t == "Mul") {
    neg <- 0L
    for (a in .args(e)) {
      s <- .sa_sign(a)
      if (s == "AnySign") return("AnySign")
      if (s == "Negative") neg <- neg + 1L
    }
    return(if (neg %% 2L == 1L) "Negative" else "Positive")
  }
  if (t == "Add") {
    all_neg <- TRUE; all_pos <- TRUE
    for (a in .args(e)) {
      s <- .sa_sign(a)
      if (s == "AnySign") return("AnySign")
      if (s == "Negative") all_pos <- FALSE else all_neg <- FALSE
    }
    return(if (all_neg) "Negative" else if (all_pos) "Positive" else "AnySign")
  }
  r <- if (t == "Symbol") NULL else .sa_rule(e)
  if (is.null(r)) "AnySign" else r[[1L]]
}

.sa_composes_as <- function(target, mono, args) {
  flipped <- if (target == "Convex") "Concave" else "Convex"
  for (i in seq_along(args)) {
    ac <- .sa_curvature(args[[i]])
    if (ac == "Affine") next
    m <- if (is.function(mono)) mono(args[[i]]) else mono[[min(i, length(mono))]]
    if (is.function(m)) m <- m(args[[i]])
    ok <- if (ac == target) m == "Increasing" else if (ac == flipped) m == "Decreasing" else FALSE
    if (!ok) return(FALSE)
  }
  TRUE
}

.sa_curvature <- function(e) {
  if (!is_sym(e) || .cvx_is_num(e) || .type(e) == "Symbol") return("Affine")
  t <- .type(e)
  args <- .args(e)
  if (t == "Mul") {
    nonconst <- NULL; count <- 0L; prod <- 1
    for (a in args) {
      if (!.cvx_is_num(a)) {
        count <- count + 1L
        nonconst <- a
        if (count > 1L) {
          warning("DCP does not support multiple non-constant arguments in multiplication", call. = FALSE)
          return("UnknownCurvature")
        }
      } else prod <- prod * .cvx_num(a)
    }
    if (is.null(nonconst)) return("Affine")
    cv <- .sa_curvature(nonconst)
    if (.sa_num_sign(prod) == "Negative") return(switch(cv, Convex = "Concave", Concave = "Convex", cv))
    return(cv)
  }
  if (t == "Add") {
    has_cvx <- FALSE; has_ccv <- FALSE
    for (a in args) {
      cv <- .sa_curvature(a)
      if (cv == "Affine") next
      if (cv == "Convex") {
        has_cvx <- TRUE
        if (has_ccv) return("UnknownCurvature")
      } else if (cv == "Concave") {
        has_ccv <- TRUE
        if (has_cvx) return("UnknownCurvature")
      } else return("UnknownCurvature")
    }
    return(if (has_cvx) "Convex" else if (has_ccv) "Concave" else "Affine")
  }
  r <- .sa_rule(e)
  if (is.null(r)) return("UnknownCurvature")
  # the exponent of a power is not an argument the rule composes over
  cargs <- if (t == "Pow") args[1L] else if (.cvx_is_norm(e)) args[-length(args)] else args
  fc <- r[[2L]]; mono <- r[[3L]]
  if (fc == "Convex") return(if (.sa_composes_as("Convex", mono, cargs)) "Convex" else "UnknownCurvature")
  if (fc == "Concave") return(if (.sa_composes_as("Concave", mono, cargs)) "Concave" else "UnknownCurvature")
  if (fc == "Affine") {
    cvx <- .sa_composes_as("Convex", mono, cargs)
    ccv <- .sa_composes_as("Concave", mono, cargs)
    return(if (cvx && ccv) "Affine" else if (cvx) "Convex" else if (ccv) "Concave" else "UnknownCurvature")
  }
  "UnknownCurvature"
}

#' Disciplined convex programming analysis
#'
#' `analyze(ex)` reports the curvature (`"Convex"`, `"Concave"`, `"Affine"`
#' or `"UnknownCurvature"`) and the sign (`"Positive"`, `"Negative"` or
#' `"AnySign"`) of a symbolic expression by the composition rules of
#' disciplined convex programming, as SymbolicAnalysis.jl: sums keep a
#' common curvature, a product may have one non constant factor (more emit
#' a warning and give `"UnknownCurvature"`), and an atom (`exp`, `log`,
#' `abs`, `Max`, `Min`, powers, `norm`) composes with arguments of matching
#' curvature and monotonicity.
#'
#' @param ex A symbolic expression (or a number).
#' @return An `AnalysisResult`: a list with `curvature`, `sign` and
#'   `gcurvature` (`NULL`; the geodesic analysis on manifolds is not yet
#'   ported).
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' analyze(exp(x) - log(y))$curvature
#' analyze((x - 1)^2 + 2 * abs(y))$curvature
#' analyze(-x^2)$sign
#' @export
analyze <- function(ex) {
  if (inherits(ex, "mtk_symvec")) ex <- .symvec_elems(ex)[[1L]]
  structure(list(curvature = .sa_curvature(ex), sign = .sa_sign(ex), gcurvature = NULL), class = "AnalysisResult")
}

#' @export
print.AnalysisResult <- function(x, ...) {
  cat("AnalysisResult: curvature = ", x$curvature, ", sign = ", x$sign, "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Linear expansion: es == A %*% vars + b (Symbolics.linear_expansion). As
# upstream it is greedy: a product is expanded in one factor (here the factor
# with the earliest variable) and the other factors, which may still hold
# variables, become the coefficient.

.cvx_linexp1 <- function(e, keys) {
  nv <- length(keys)
  zero <- function() rep(list(as_sym(0)), nv)
  if (!is_sym(e)) e <- as_sym(e)
  if (!.cvx_depends(e, keys)) return(list(a = zero(), b = e, ok = TRUE))
  t <- .type(e)
  if (t == "Symbol") {
    a <- zero()
    a[[match(.key(e), keys)]] <- as_sym(1)
    return(list(a = a, b = as_sym(0), ok = TRUE))
  }
  if (t == "Add") {
    a <- zero(); b <- as_sym(0)
    for (arg in .args(e)) {
      r <- .cvx_linexp1(arg, keys)
      if (!r$ok) return(r)
      for (j in seq_len(nv)) a[[j]] <- a[[j]] + r$a[[j]]
      b <- b + r$b
    }
    return(list(a = a, b = b, ok = TRUE))
  }
  if (t == "Mul") {
    args <- .args(e)
    first <- vapply(args, function(f) {
      w <- which(keys %in% .keys(.cvx_free(f)))
      if (length(w)) min(w) else NA_integer_
    }, 0L)
    k <- which(first == min(first, na.rm = TRUE))[1L]
    r <- .cvx_linexp1(args[[k]], keys)
    if (!r$ok) return(r)
    rest <- as_sym(1)
    for (j in seq_along(args)[-k]) rest <- rest * args[[j]]
    return(list(a = lapply(r$a, function(x) x * rest), b = r$b * rest, ok = TRUE))
  }
  if (t == "Pow" && .cvx_is_num(.args(e)[[2L]]) && .cvx_num(.args(e)[[2L]]) == 1)
    return(.cvx_linexp1(.args(e)[[1L]], keys))
  list(a = zero(), b = as_sym(0), ok = FALSE)
}

.cvx_linear_expansion <- function(es, vars) {
  keys <- .keys(vars)
  es <- .symvec_elems(es)
  A <- matrix(list(), length(es), length(vars))
  b <- vector("list", length(es))
  for (i in seq_along(es)) {
    r <- .cvx_linexp1(es[[i]], keys)
    if (!r$ok) return(list(A = A, b = b, islin = FALSE))
    for (j in seq_along(vars)) A[[i, j]] <- r$a[[j]]
    b[[i]] <- r$b
  }
  list(A = A, b = b, islin = TRUE)
}
