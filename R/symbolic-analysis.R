# Disciplined convex programming analysis of symbolic expressions (the
# curvature and sign of SymbolicAnalysis.jl) and the linear expansion of
# Symbolics.jl, with the matrix products and norm atoms that trace a convex
# problem.
#
# Ported from SymbolicAnalysis.jl 0.5.8
# (https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl), Copyright (c) 2024
# Vaibhav Dixit and contributors, MIT License (inst/CONVEX-LICENSES), and
# from linear_expansion of Symbolics.jl
# (https://github.com/JuliaSymbolics/Symbolics.jl), Copyright (c) 2021:
# Shashi Gowda, Yingbo Ma, Chris Rackauckas, Julia Computing, MIT License
# (inst/COPYRIGHTS).

# ---------------------------------------------------------------------------
# Tracing: the symbolic vectors of modelingtoolkitize (mtk_symvec,
# R/modelingtoolkitize-trace.R), with matrix products and the atoms of a
# convex problem. Julia traces array expressions (`sum(abs.(w))`,
# `maximum(w)`, `u' * P * u`) to single terms; while a convex problem is
# traced (`.cvx_state$active`) the symbolic vector methods keep the same
# structure as one uninterpreted function node per atom, with its arguments
# in the order written: sum_abs, maximum_abs, sum_sq, maximum, minimum, max,
# min, abs, self_prod (w' w), quad_prod (v' P v, P appended by columns) and
# rdiv (a division by an expression in the variables).

.cvx_state <- new.env(parent = emptyenv())
.cvx_state$active <- FALSE
.cvx_state$optkeys <- character()

.cvx_fn <- function(name, args) do.call(Function(name), unname(args))
.cvx_sv1 <- function(e) .symvec(list(e))
.cvx_same <- function(a, b) length(a) == length(b) && all(vapply(seq_along(a), function(i)
  identical(.key(as_sym(a[[i]])), .key(as_sym(b[[i]]))), TRUE))

.cvx_trace_abs <- function(x) {
  el <- .symvec_elems(x)
  out <- .symvec(lapply(el, function(z) if (.cvx_is_num(z)) as_sym(abs(.cvx_num(z))) else .cvx_fn("abs", list(z))),
                 names(x))
  attr(out, "cvx_map") <- "abs"
  attr(out, "cvx_arg") <- el
  out
}

.cvx_trace_ops <- function(gen, e1, e2) {
  if (gen == "^" && inherits(e1, "mtk_symvec") && is.numeric(e2) && length(e2) == 1L && e2 == 2) {
    el <- .symvec_elems(e1)
    out <- .symvec(lapply(el, function(z) z^2), names(e1))
    attr(out, "cvx_map") <- "sq"
    attr(out, "cvx_arg") <- el
    return(out)
  }
  if (gen == "/") {
    b <- .symvec_elems(e2)
    if (!any(vapply(b, .cvx_depends, TRUE, keys = .cvx_state$optkeys))) return(NULL)
    a <- .symvec_elems(e1)
    n <- max(length(a), length(b))
    return(.symvec(lapply(seq_len(n), function(i) {
      num <- a[[(i - 1L) %% length(a) + 1L]]; den <- b[[(i - 1L) %% length(b) + 1L]]
      if (.cvx_depends(den, .cvx_state$optkeys)) .cvx_fn("rdiv", list(num, den)) else num / den
    })))
  }
  NULL
}

.cvx_flatten_fn <- function(el, name) {
  out <- list()
  for (z in el) {
    if (is_sym(z) && .type(z) == "Function" && identical(.name(z), name)) out <- c(out, .cvx_flatten_fn(.args(z), name))
    else out[[length(out) + 1L]] <- z
  }
  out
}

.cvx_trace_summary <- function(gen, args) {
  if (length(args) == 1L && !is.null(m <- attr(args[[1L]], "cvx_map"))) {
    w <- attr(args[[1L]], "cvx_arg")
    if (gen == "sum") return(.cvx_sv1(.cvx_fn(if (m == "abs") "sum_abs" else "sum_sq", w)))
    if (m == "abs") {
      if (gen == "min")
        stop("`minimum(abs.(w))` is neither convex nor concave and cannot be lowered. Route to a general ",
             "OptimizationProblem/NLP solver.", call. = FALSE)
      return(.cvx_sv1(.cvx_fn("maximum_abs", w)))
    }
  }
  el <- .symvec_elems(args)
  if (gen == "sum") return(.cvx_sv1(if (length(el)) Reduce(`+`, el) else as_sym(0)))
  if (all(vapply(el, .cvx_is_num, TRUE))) {
    f <- if (gen == "max") base::max else base::min
    return(.cvx_sv1(as_sym(f(vapply(el, .cvx_num, 0)))))
  }
  if (length(el) == 1L) return(.cvx_sv1(el[[1L]]))
  # one vector argument is Julia's maximum(w) / minimum(w); several arguments
  # its scalar max(a, b, ...), which the lowering flattens
  if (length(args) == 1L) return(.cvx_sv1(.cvx_fn(if (gen == "max") "maximum" else "minimum", el)))
  .cvx_sv1(.cvx_fn(gen, .cvx_flatten_fn(el, gen)))
}

#' @export
t.mtk_symvec <- function(x) structure(list(v = .symvec_elems(x), mid = list()), class = "mtk_symadj")

# v' M as the scalar elements of a row vector
.cvx_adj_row <- function(x) {
  v <- x$v
  if (!length(x$mid)) return(v)
  M <- Reduce(`%*%`, x$mid)
  lapply(seq_len(ncol(M)), function(j) {
    out <- as_sym(0)
    for (i in seq_along(v)) if (M[i, j] != 0) out <- out + M[i, j] * v[[i]]
    out
  })
}

.cvx_dot <- function(a, b) {
  out <- as_sym(0)
  for (i in seq_along(a)) out <- out + a[[i]] * b[[i]]
  out
}

# Matrix products and norm seen by a traced function
.cvx_matmul <- function(x, y) {
  if (inherits(x, "mtk_symadj")) {
    if (inherits(y, "mtk_symvec")) {
      w <- .symvec_elems(y)
      if (.cvx_same(x$v, w)) {
        if (!length(x$mid)) return(.cvx_sv1(.cvx_fn("self_prod", w)))
        P <- Reduce(`%*%`, x$mid)
        return(.cvx_sv1(.cvx_fn("quad_prod", c(w, as.list(as.numeric(P))))))
      }
      r <- .cvx_adj_row(x)
      if (length(r) != length(w)) stop("non-conformable arguments", call. = FALSE)
      return(.cvx_sv1(.cvx_dot(r, w)))
    }
    if (is.numeric(y) && is.null(dim(y))) {
      r <- .cvx_adj_row(x)
      if (length(r) != length(y)) stop("non-conformable arguments", call. = FALSE)
      return(.cvx_sv1(.cvx_dot(r, as.list(y))))
    }
    if (!is.numeric(y))
      stop("`u' * P * u` / `quad_form(u, P)` is lowered only for a constant numeric matrix `P`; got `",
           paste(format(y), collapse = ", "), "`. A `P` built from `p` moves the cone matrix with theta and ",
           "cannot be canonicalized: pass the quadratic form differently, or use ",
           "`solve(remake(prob, p = ...), alg)`.", call. = FALSE)
    x$mid[[length(x$mid) + 1L]] <- as.matrix(y)
    return(x)
  }
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
# of base::norm map to the vector norms they compute. Of a symbolic matrix,
# type "2" is the spectral norm (opnorm) and "F" the norm of its entries.
.cvx_norm <- function(x, type = "O", ...) {
  if (is.list(x) && !is.null(dim(x)) && !inherits(x, "mtk_symvec")) {
    t1 <- toupper(substr(as.character(type), 1L, 1L))
    if (t1 == "2") return(opnorm(x))
    if (t1 %in% c("F", "E")) return(.cvx_norm(.symvec(as.list(x)), "2"))
    stop("the norm type '", type, "' of a symbolic matrix is not a supported atom; use \"2\" or \"F\"", call. = FALSE)
  }
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
# Matrix and quadratic atoms: numeric on numbers, symbolic nodes on traced
# values

# A symbolic matrix: a list matrix (matrix(c(u[1], 0, 0, u[2]), 2)) of
# expressions; NULL for a numeric matrix
.cvx_symmat <- function(X, name) {
  if (inherits(X, "mtk_symvec") || is_sym(X)) stop("`", name, "` requires a matrix argument.", call. = FALSE)
  if (is.numeric(X)) {
    if (is.null(dim(X)) || length(dim(X)) != 2L) stop("`", name, "` requires a matrix argument.", call. = FALSE)
    return(NULL)
  }
  if (!is.list(X) || length(dim(X)) != 2L) stop("`", name, "` requires a matrix argument.", call. = FALSE)
  el <- lapply(as.vector(X), function(z) if (inherits(z, "mtk_symvec")) .symvec_elems(z)[[1L]] else as_sym(z))
  if (all(vapply(el, .cvx_is_num, TRUE)))
    return(list(numeric = matrix(vapply(el, .cvx_num, 0), nrow(X), ncol(X))))
  list(el = el, nrow = nrow(X), ncol = ncol(X))
}

.cvx_matrix_atom <- function(X, name, numfun) {
  s <- .cvx_symmat(X, name)
  if (is.null(s)) return(numfun(X))
  if (!is.null(s$numeric)) return(numfun(s$numeric))
  if (!(s$nrow > 0 && s$ncol > 0)) stop("`", name, "` requires a nonempty matrix.", call. = FALSE)
  .cvx_sv1(.cvx_fn(name, c(s$el, list(s$nrow))))
}

#' Matrix and quadratic atoms of disciplined convex programs
#'
#' The atoms that [ConvexMOI()] lowers to cones, as in ConvexOptimization.jl
#' and SymbolicAnalysis.jl: `eigmax(X)` and `eigmin(X)` (the extreme
#' eigenvalues of a symmetric matrix, lowered to a positive semidefinite
#' cone), `opnorm(X)` (the spectral norm, `norm(X, "2")`), `logdet(X)` (the
#' log determinant, a log determinant cone) and `quad_form(x, P)`
#' (`x' P x` for a constant positive semidefinite `P`, a rotated second order
#' cone). On numbers they compute the value; on the traced variables of a
#' convex problem they return the atom. A symbolic matrix is a list matrix of
#' expressions, such as `matrix(c(u[1], 0, 0, 1 - u[1]), 2)`.
#'
#' @param X A numeric or symbolic matrix.
#' @param x A numeric or symbolic vector.
#' @param P A numeric matrix.
#' @return A number, or the traced atom.
#' @examples
#' eigmax(matrix(c(2, 1, 1, 2), 2))
#' quad_form(c(1, 2), diag(2))
#' prob <- ConvexOptimizationProblem(function(u, p) eigmax(matrix(c(u[1], 0, 0, 1 - u[1]), 2)), 0)
#' solve(prob, ConvexMOI())$u
#' @export
eigmax <- function(X) .cvx_matrix_atom(X, "eigmax", function(M) {
  if (nrow(M) != ncol(M)) stop("`eigmax` requires a square symmetric matrix; got size (", nrow(M), ", ", ncol(M), ").",
                               call. = FALSE)
  max(eigen(M, symmetric = isSymmetric(M), only.values = TRUE)$values)
})

#' @rdname eigmax
#' @export
eigmin <- function(X) .cvx_matrix_atom(X, "eigmin", function(M) {
  if (nrow(M) != ncol(M)) stop("`eigmin` requires a square symmetric matrix; got size (", nrow(M), ", ", ncol(M), ").",
                               call. = FALSE)
  min(eigen(M, symmetric = isSymmetric(M), only.values = TRUE)$values)
})

#' @rdname eigmax
#' @export
opnorm <- function(X) .cvx_matrix_atom(X, "opnorm", function(M) if (length(M)) max(svd(M, 0L, 0L)$d) else 0)

#' @rdname eigmax
#' @export
logdet <- function(X) .cvx_matrix_atom(X, "logdet", function(M) {
  d <- determinant(M, logarithm = TRUE)
  if (d$sign < 0) stop("DomainError: logdet of a matrix with a negative determinant", call. = FALSE)
  as.numeric(d$modulus)
})

#' @rdname eigmax
#' @export
quad_form <- function(x, P) {
  if (!inherits(x, "mtk_symvec") && !is_sym(x)) {
    if (length(x) != nrow(P)) stop("DimensionMismatch: x and P must have same length", call. = FALSE)
    return(sum(x * (P %*% x)))
  }
  w <- .symvec_elems(x)
  if (!is.numeric(P))
    stop("`u' * P * u` / `quad_form(u, P)` is lowered only for a constant numeric matrix `P`; got `",
         paste(format(P), collapse = ", "), "`. A `P` built from `p` moves the cone matrix with theta and cannot be ",
         "canonicalized: pass the quadratic form differently, or use `solve(remake(prob, p = ...), alg)`.",
         call. = FALSE)
  P <- as.matrix(P)
  .cvx_sv1(.cvx_fn("quad_form", c(w, as.list(as.numeric(P)))))
}

# the vector and the matrix of a quad_form / quad_prod node
.cvx_quad_parts <- function(e) {
  a <- .args(e)
  L <- length(a)
  n <- as.integer(round((-1 + sqrt(1 + 4 * L)) / 2))
  list(v = a[seq_len(n)], P = matrix(vapply(a[n + seq_len(n * n)], .cvx_num, 0), n, n))
}

# the entries, rows and columns of a matrix atom node
.cvx_mat_parts <- function(e) {
  a <- .args(e)
  m <- as.integer(.cvx_num(a[[length(a)]]))
  el <- a[-length(a)]
  list(el = el, nrow = m, ncol = length(el) %/% m)
}

# ---------------------------------------------------------------------------
# Expression helpers

.cvx_is_num <- function(e) is_sym(e) && .type(e) %in% c("Integer", "Rational", "Float")
.cvx_num <- function(e) as.numeric(e)
.cvx_free <- function(e) if (.cvx_is_num(e)) list() else free_symbols(e)
.cvx_depends <- function(e, keys) is_sym(e) && !.cvx_is_num(e) && any(.keys(.cvx_free(e)) %in% keys)
.cvx_is_norm <- function(e) .type(e) == "Function" && identical(.name(e), "norm")
.cvx_fname <- function(e) if (.type(e) == "Function") .name(e) else .type(e)

# ---------------------------------------------------------------------------
# Curvature and sign (SymbolicAnalysis.jl rules)

.SA_SIGN <- c("Positive", "Negative", "AnySign")
.SA_CURV <- c("Convex", "Concave", "Affine", "UnknownCurvature")

.sa_num_sign <- function(v) if (v < 0) "Negative" else "Positive"

# the sign and curvature of an array argument (getsign / getcurvature of an
# AbstractArray)
.sa_vsign <- function(es) {
  s <- vapply(es, .sa_sign, "")
  if (length(s) && all(s == "Negative")) "Negative" else if (all(s == "Positive")) "Positive" else "AnySign"
}

.sa_vcurv <- function(es) {
  has_cvx <- FALSE; has_ccv <- FALSE
  for (e in es) {
    cv <- .sa_curvature(e)
    if (cv == "Affine") next
    if (cv == "Convex") {
      if (has_ccv) return("UnknownCurvature")
      has_cvx <- TRUE
    } else if (cv == "Concave") {
      if (has_cvx) return("UnknownCurvature")
      has_ccv <- TRUE
    } else return("UnknownCurvature")
  }
  if (has_cvx) "Convex" else if (has_ccv) "Concave" else "Affine"
}

.sa_incr_if_pos <- function(group) {
  s <- .sa_vsign(group)
  if (s == "AnySign") "AnyMono" else if (s == "Positive") "Increasing" else "Decreasing"
}

# The curvature of f(group_1, ...) for an atom of curvature fc and
# monotonicity mono (a vector whose last entry repeats, or a function of the
# argument group such as increasing_if_positive)
.sa_compose <- function(fc, mono, groups) {
  composes <- function(target) {
    flipped <- if (target == "Convex") "Concave" else "Convex"
    for (i in seq_along(groups)) {
      ac <- .sa_vcurv(groups[[i]])
      if (ac == "Affine") next
      m <- if (is.function(mono)) mono else mono[[min(i, length(mono))]]
      if (is.function(m)) m <- m(groups[[i]])
      ok <- if (ac == target) m == "Increasing" else if (ac == flipped) m == "Decreasing" else FALSE
      if (!ok) return(FALSE)
    }
    TRUE
  }
  if (fc == "Convex") return(if (composes("Convex")) "Convex" else "UnknownCurvature")
  if (fc == "Concave") return(if (composes("Concave")) "Concave" else "UnknownCurvature")
  if (fc == "Affine") {
    cvx <- composes("Convex"); ccv <- composes("Concave")
    return(if (cvx && ccv) "Affine" else if (cvx) "Convex" else if (ccv) "Concave" else "UnknownCurvature")
  }
  "UnknownCurvature"
}

.sa_is_posdef <- function(P) {
  if (!isSymmetric(unname(P), tol = 0)) return(FALSE)
  ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
  all(ev > 0)
}

# The DCP rule of a call: list(sign, curvature, monotonicity, argument
# groups), or NULL for a call without one
.sa_rule <- function(e) {
  args <- .args(e)
  t <- .type(e)
  one <- function(s, c, m) list(s, c, m, lapply(args, list))
  if (t == "Pow") {
    i <- args[[2L]]
    if (!.cvx_is_num(i)) {
      # power_rule for a non-constant exponent: a positive constant base
      b <- args[[1L]]
      if (.cvx_is_num(b) && .cvx_num(b) > 0) {
        bv <- .cvx_num(b)
        return(one("Positive", "Convex", c("AnyMono", if (bv > 1) "Increasing" else if (bv < 1) "Decreasing" else "AnyMono")))
      }
      return(one("AnySign", "UnknownCurvature", "AnyMono"))
    }
    i <- .cvx_num(i)
    return(if (i == 1) one("AnySign", "Affine", "Increasing")
      else if (i == round(i) && i %% 2 == 0) one("Positive", "Convex", .sa_incr_if_pos)
      else if (i == round(i)) one("Positive", "Convex", "Increasing")
      else if (i >= 1) one("Positive", "Convex", "Increasing")
      else if (i > 0) one("Positive", "Concave", "Increasing")
      else one("Positive", "Convex", "Increasing"))
  }
  if (t != "Function") return(NULL)
  vec1 <- function(s, c, m) list(s, c, m, list(args))
  mat <- function(s, c, m) list(s, c, m, list(args[-length(args)]))
  switch(.name(e),
    exp = , exp2 = , exp10 = one("Positive", "Convex", "Increasing"),
    expm1 = one("AnySign", "Convex", "Increasing"),
    log = , log2 = , log10 = , log1p = one("AnySign", "Concave", "Increasing"),
    Abs = , abs = , abs2 = , cosh = , hypot = one("Positive", "Convex", .sa_incr_if_pos),
    Max = , max = one("AnySign", "Convex", "Increasing"),
    Min = , min = one("AnySign", "Concave", "Increasing"),
    maximum = vec1("AnySign", "Convex", "Increasing"),
    minimum = vec1("AnySign", "Concave", "Increasing"),
    # sum and maximum over the elementwise abs or square of w (a convex
    # broadcast of sign Positive): the reducer keeps its own sign AnySign
    sum_abs = , maximum_abs = , sum_sq = {
      c <- .sa_compose("Convex", .sa_incr_if_pos, list(args))
      list("AnySign", c, "Increasing", list())
    },
    # w' w: a product of two non constant factors
    self_prod = {
      warning("DCP does not support multiple non-constant arguments in multiplication", call. = FALSE)
      list("AnySign", "UnknownCurvature", "AnyMono", list())
    },
    quad_form = , quad_prod = {
      q <- .cvx_quad_parts(e)
      if (.sa_is_posdef(q$P))
        list("Positive", "Convex", list(if (all(q$P >= 0)) .sa_incr_if_pos else "AnyMono", "Increasing"), list(q$v))
      else list("AnySign", "UnknownCurvature", "AnyMono", list(q$v))
    },
    eigmax = mat("AnySign", "Convex", "AnyMono"),
    eigmin = mat("AnySign", "Concave", "AnyMono"),
    opnorm = mat("Positive", "Convex", "AnyMono"),
    logdet = mat("AnySign", "Concave", "AnyMono"),
    norm = {
      p <- args[[length(args)]]
      w <- list(args[-length(args)])
      pv <- if (identical(.key(p), .key(Symbol(".norm_Inf")))) Inf else if (.cvx_is_num(p)) .cvx_num(p) else NA
      if (is.na(pv)) list("Positive", "UnknownCurvature", "AnyMono", w)
      else if (pv >= 1) list("Positive", "Convex", .sa_incr_if_pos, w)
      else if (pv > 0) list("Positive", if (.sa_vsign(w[[1L]]) == "Positive") "Concave" else "UnknownCurvature",
                            "Increasing", w)
      else list("Positive", "UnknownCurvature", "AnyMono", w)
    },
    # the rule of `/` for a denominator in the variables
    rdiv = {
      num <- args[[1L]]; den <- args[[2L]]
      if (.cvx_is_num(num) && .cvx_num(num) != 0 && .sa_sign(den) == "Positive") {
        if (.cvx_num(num) > 0) one("Positive", "Convex", c("AnyMono", "Decreasing"))
        else one("Negative", "Concave", c("AnyMono", "Increasing"))
      } else one("AnySign", "UnknownCurvature", "AnyMono")
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
    if (cv == "Affine") return("Affine")
    if (.sa_num_sign(prod) == "Negative") return(switch(cv, Convex = "Concave", Concave = "Convex", cv))
    return(cv)
  }
  if (t == "Add") return(.sa_vcurv(args))
  r <- .sa_rule(e)
  if (is.null(r)) return("UnknownCurvature")
  if (!length(r[[4L]])) return(r[[2L]])
  .sa_compose(r[[2L]], r[[3L]], r[[4L]])
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
