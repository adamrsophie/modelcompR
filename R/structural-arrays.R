# Arithmetic on array variables, after the array support of Symbolics and
# the scalarization of array equations in
# StateSelection.jl/lib/ModelingToolkitTearing/src/tearingstate.jl
# (scalarize_tearing_state_eqs!). An array variable (`Variable(..., dims =)`
# or `variables("x(t)[1:3]")`) is a list of scalar variables x[1](t), ...;
# these methods make elementwise operations, `%*%` and `sum` return arrays of
# scalar expressions, and an equation between arrays is a list of scalar
# equations, so the structural pipeline only ever sees scalars.

.arr_elems <- function(x) if (inherits(x, "mtk_array")) unclass(x) else x

.arr_get <- function(x, i) {
  if (is_sym(x)) return(x)
  if (is.list(x)) return(x[[if (length(x) == 1L) 1L else i]])
  x[[if (length(x) == 1L) 1L else i]]
}

.arr_new <- function(vals, like = NULL) {
  .arr_install_ops_hook()
  d <- if (inherits(like, "mtk_array")) dim(like) else NULL
  out <- structure(vals, class = c("mtk_array", "list"))
  if (!is.null(d) && prod(d) == length(vals)) dim(out) <- d
  out
}

.arr_ops <- function(gen, e1, e2) {
  f <- get(gen, envir = baseenv(), mode = "function")
  if (missing(e2)) return(.arr_new(lapply(unclass(e1), f), e1))
  a <- .arr_elems(e1); b <- .arr_elems(e2)
  la <- if (is_sym(a)) 1L else length(a)
  lb <- if (is_sym(b)) 1L else length(b)
  n <- max(la, lb)
  if (la != lb && min(la, lb) != 1L)
    stop("arrays of lengths ", la, " and ", lb, " cannot be combined elementwise", call. = FALSE)
  .arr_new(lapply(seq_len(n), function(i) f(.arr_get(a, i), .arr_get(b, i))),
           if (inherits(e1, "mtk_array") && la == n) e1 else e2)
}

#' @export
Ops.mtk_array <- function(e1, e2) if (missing(e2)) .arr_ops(.Generic, e1) else .arr_ops(.Generic, e1, e2)

# symcompR's Ops.sym always wins the dispatch of sym * array
# (chooseOpsMethod.sym); its dispatch hook hands arrays back to .arr_ops.
# Installed when the first array variable is created.
.arr_install_ops_hook <- function() {
  hooks <- get0(".ops_hooks", envir = asNamespace("symcompR"), inherits = FALSE)
  if (!is.environment(hooks) || isTRUE(attr(hooks$dispatch, "mtk_array"))) return(invisible())
  old <- hooks$dispatch
  hooks$dispatch <- structure(function(gen, e1, e2) {
    if (inherits(e1, "mtk_array") || inherits(e2, "mtk_array")) return(.arr_ops(gen, e1, e2))
    if (is.null(old)) NULL else old(gen, e1, e2)
  }, mtk_array = TRUE)
  invisible()
}

#' @export
Math.mtk_array <- function(x, ...) {
  f <- get(.Generic, envir = baseenv(), mode = "function")
  .arr_new(lapply(unclass(x), function(z) f(z, ...)), x)
}

#' @export
Summary.mtk_array <- function(..., na.rm = FALSE) {
  if (!(.Generic %in% c("sum", "prod"))) stop(.Generic, " is not defined for array variables", call. = FALSE)
  els <- unlist(lapply(list(...), function(z) if (is.list(z)) unclass(z) else as.list(z)), recursive = FALSE)
  Reduce(if (.Generic == "sum") `+` else `*`, els)
}

# Rows by columns view of a matrix, a vector (a column, or a row when
# `row`) or an array variable.
.arr_matrix <- function(x, row = FALSE) {
  d <- dim(x)
  v <- if (is.list(x)) unclass(x) else as.list(x)
  if (is.null(d)) d <- if (row) c(1L, length(v)) else c(length(v), 1L)
  list(v = v, nr = d[1L], nc = d[2L])
}

# %*% with the conventions of base R for vectors: a vector on the left is a
# row when that conforms, and one on the right is a row when the left
# operand has a single column.
# The matrixOps group generic exists from R 4.4.0; the method is registered
# under its name with a function that roxygen does not take for an S3 method.
#' @rawNamespace if (getRversion() >= "4.4.0") S3method(matrixOps, mtk_array, .arr_matrixops)
#' @noRd
.arr_matrixops <- function(x, y) {
  if (.Generic != "%*%") stop(.Generic, " is not defined for array variables", call. = FALSE)
  b <- .arr_matrix(y)
  a <- .arr_matrix(x, row = is.null(dim(x)) && length(x) == b$nr)
  if (a$nc != b$nr && is.null(dim(y)) && a$nc == 1L) b <- .arr_matrix(y, row = TRUE)
  if (a$nc != b$nr) stop("non-conformable arguments", call. = FALSE)
  out <- vector("list", a$nr * b$nc)
  for (i in seq_len(a$nr)) for (j in seq_len(b$nc)) {
    terms <- lapply(seq_len(a$nc), function(k) a$v[[i + (k - 1L) * a$nr]] * b$v[[k + (j - 1L) * b$nr]])
    out[[i + (j - 1L) * a$nr]] <- as_sym(Reduce(`+`, terms))
  }
  res <- structure(out, class = c("mtk_array", "list"))
  if (b$nc > 1L || a$nr == 1L) dim(res) <- c(a$nr, b$nc)
  res
}
