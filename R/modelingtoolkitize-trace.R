# modelingtoolkitize, after ModelingToolkitBase/src/modelingtoolkitize/
# common.jl (construct_vars, define_params, construct_params,
# lhs_from_mass_matrix, trace_rhs, defaults_from_u0_p), odeproblem.jl,
# nonlinearproblem.jl, optimizationproblem.jl and sdeproblem.jl.
#
# Julia traces the numerical function with symbolic `Num` arrays. Here the
# state, parameter and time arguments are "symbolic vectors" (class
# mtk_symvec): lists of symcompR expressions whose indexing, arithmetic, math
# and summary methods return symbolic vectors, so a function written for
# numeric vectors (u[1], p["k"], sum(u), c(...)) computes the symbolic
# right-hand side. As for numeric vectors, u[1] and u[[1]] are vectors of
# length one. Control flow on symbolic values cannot be traced, as in Julia.

.symvec <- function(xs, nms = NULL) {
  xs <- lapply(xs, function(z) if (inherits(z, "mtk_symvec")) unclass(z)[[1L]] else as_sym(z))
  structure(xs, names = nms, class = "mtk_symvec")
}

# Elements of a traced value as a list of sym.
.symvec_elems <- function(z) {
  if (is.null(z)) return(list())
  if (inherits(z, "mtk_symvec")) return(unclass(z))
  if (is_sym(z)) return(list(z))
  if (is.numeric(z) || is.logical(z)) return(lapply(as.numeric(z), as_sym))
  if (is.list(z)) return(do.call(c, c(list(list()), lapply(unname(z), .symvec_elems))))
  stop("cannot trace a value of class ", class(z)[1L],
       "; return a vector built with c() or list() from u, p and t", call. = FALSE)
}

#' @export
`[.mtk_symvec` <- function(x, i, ...) {
  if (missing(i)) return(x)
  out <- unclass(x)[i]
  if (any(vapply(out, is.null, TRUE))) stop("index out of bounds in a traced vector", call. = FALSE)
  structure(out, class = "mtk_symvec")
}

#' @export
`[[.mtk_symvec` <- function(x, i, ...) {
  el <- unclass(x)[[i]]
  .symvec(list(el), if (is.character(i)) i else NULL)
}

#' @export
`$.mtk_symvec` <- function(x, name) x[[name]]

#' @export
`[<-.mtk_symvec` <- function(x, i, value) {
  y <- unclass(x)
  y[i] <- .symvec_elems(value)
  .symvec(y, names(y))
}

#' @export
`[[<-.mtk_symvec` <- function(x, i, value) {
  y <- unclass(x)
  v <- .symvec_elems(value)
  if (length(v) != 1L) stop("replacement has length ", length(v), call. = FALSE)
  y[[i]] <- v[[1L]]
  .symvec(y, names(y))
}

#' @export
c.mtk_symvec <- function(...) {
  parts <- list(...)
  nms <- character()
  out <- list()
  argn <- names(parts) %||% rep("", length(parts))
  for (k in seq_along(parts)) {
    el <- .symvec_elems(parts[[k]])
    pn <- names(parts[[k]]) %||% rep("", length(el))
    if (length(pn) != length(el)) pn <- rep("", length(el))
    # names as base::c(): a = x gives "a", a = c(x, y) "a1", "a2", a = c(b = x) "a.b"
    if (nzchar(argn[k])) {
      pn <- if (length(el) == 1L && !nzchar(pn)) argn[k]
            else ifelse(nzchar(pn), paste0(argn[k], ".", pn), paste0(argn[k], seq_along(el)))
    }
    out <- c(out, el)
    nms <- c(nms, pn)
  }
  .symvec(out, if (any(nzchar(nms))) nms else NULL)
}

#' @export
as.list.mtk_symvec <- function(x, ...) {
  out <- lapply(seq_along(x), function(i) .symvec(list(unclass(x)[[i]])))
  names(out) <- names(x)
  out
}

#' @export
Ops.mtk_symvec <- function(e1, e2) {
  # the convex tracer keeps elementwise squares and divisions by variables as
  # atoms (R/symbolic-analysis.R)
  if (.cvx_state$active && !missing(e2)) {
    r <- .cvx_trace_ops(.Generic, e1, e2)
    if (!is.null(r)) return(r)
  }
  op <- get(.Generic, envir = baseenv())
  a <- .symvec_elems(e1)
  if (missing(e2)) return(.symvec(lapply(a, function(z) op(z)), names(e1)))
  b <- .symvec_elems(e2)
  na <- length(a); nb <- length(b)
  if (!na || !nb) return(.symvec(list()))
  n <- max(na, nb)
  res <- lapply(seq_len(n), function(i) {
    r <- op(a[[(i - 1L) %% na + 1L]], b[[(i - 1L) %% nb + 1L]])
    if (is.logical(r)) as_sym(r) else r
  })
  nms <- if (na >= nb && inherits(e1, "mtk_symvec")) names(e1) else names(e2)
  .symvec(res, if (length(nms) == n) nms else NULL)
}

#' @export
chooseOpsMethod.mtk_symvec <- function(x, y, mx, my, cl, reverse) TRUE

#' @export
Math.mtk_symvec <- function(x, ...) {
  if (.cvx_state$active && .Generic == "abs") return(.cvx_trace_abs(x))
  f <- get(.Generic, envir = baseenv())
  .symvec(lapply(unclass(x), function(z) f(z, ...)), names(x))
}

#' @export
Summary.mtk_symvec <- function(..., na.rm = FALSE) {
  if (.cvx_state$active && .Generic %in% c("sum", "max", "min")) return(.cvx_trace_summary(.Generic, list(...)))
  el <- .symvec_elems(list(...))
  r <- switch(.Generic,
    sum = if (length(el)) Reduce(`+`, el) else .S0(),
    prod = if (length(el)) Reduce(`*`, el) else .S1(),
    max = do.call(symcompR::Max, el),
    min = do.call(symcompR::Min, el),
    stop(.Generic, "() cannot be traced symbolically", call. = FALSE))
  .symvec(list(r))
}

#' @export
format.mtk_symvec <- function(x, ...) vapply(unclass(x), format, "")

#' @export
print.mtk_symvec <- function(x, ...) {
  cat("Traced vector:\n")
  print(stats::setNames(format.mtk_symvec(x), names(x)), quote = FALSE)
  invisible(x)
}

# ---------------------------------------------------------------------------

.SUBSCRIPTS <- c("\u2080", "\u2081", "\u2082", "\u2083", "\u2084", "\u2085", "\u2086", "\u2087",
                 "\u2088", "\u2089")
.subscript <- function(i) paste(.SUBSCRIPTS[as.integer(strsplit(as.character(i), "")[[1L]]) + 1L], collapse = "")

.mtkize_fun <- function(f) if (is.function(f)) f else f$f

# construct_vars: names from u_names, the problem's system, the names of u0,
# or x₁, x₂, ...
.mtkize_vars <- function(prob, t, u_names) {
  u0 <- prob$u0
  n <- length(u0)
  if (!n) return(list())
  sys <- prob$sys
  bounds <- NULL
  if (!is.null(u_names)) {
    if (length(u_names) != n)
      stop("Number of unknowns (", n, ") does not match number of names (", length(u_names), ").", call. = FALSE)
    nms <- as.character(u_names)
  } else if (!is.null(sys)) {
    us <- unknowns(sys)
    nms <- vapply(us, .name, "")
    bounds <- lapply(us, function(u) if (hasbounds(u)) getbounds(u) else NULL)
  } else if (!is.null(names(u0)) && all(nzchar(names(u0)))) {
    nms <- names(u0)
  } else {
    nms <- paste0("x", vapply(seq_len(n), .subscript, ""))
  }
  lapply(seq_len(n), function(i) {
    Variable(nms[i], t, bounds = bounds[[i]])
  })
}

# define_params: the traced parameter object (a symbolic vector, or a named
# list of them for list parameters) and the flat list of parameters.
.mtkize_params <- function(prob, p_names) {
  p <- prob$p
  if (is.null(p) || !length(p)) return(list(obj = .symvec(list()), vars = list(), values = numeric()))
  sys <- prob$sys
  if (is.null(p_names) && !is.null(sys)) p_names <- vapply(.problem_params(sys), .name, "")
  scalar <- is.numeric(p) && length(p) == 1L && is.null(names(p))
  nms <- if (!is.null(p_names)) {
    if (length(p_names) != length(p))
      stop("Number of parameters (", length(p), ") does not match number of names (", length(p_names), ").",
           call. = FALSE)
    as.character(p_names)
  } else if (!is.null(names(p)) && all(nzchar(names(p)))) {
    names(p)
  } else if (scalar) {
    "\u03B1"
  } else {
    paste0("\u03B1", vapply(seq_along(p), .subscript, ""))
  }
  if (is.list(p)) {
    bad <- !vapply(p, function(v) is.numeric(v) && length(v) == 1L, TRUE)
    if (any(bad)) stop("modelingtoolkitize() supports numeric parameters only; element(s) ",
                       paste(which(bad), collapse = ", "), " are not numbers", call. = FALSE)
  }
  vars <- lapply(nms, Parameter)
  obj <- if (is.list(p)) stats::setNames(lapply(vars, function(v) .symvec(list(v))), names(p))
         else .symvec(vars, names(p))
  list(obj = obj, vars = vars, values = as.numeric(unlist(p)))
}

.trace <- function(f, args) {
  r <- tryCatch(do.call(f, args), error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("mtk_symvec", msg, fixed = TRUE))
      stop("tracing mixed a symbolic expression with an element of u or p taken by iterating over ",
           "them (for loops, lapply); index with u[i] instead", call. = FALSE)
    stop("modelingtoolkitize() could not trace the function: ", msg, call. = FALSE)
  })
  out <- .symvec_elems(r)
  out
}

# lhs_from_mass_matrix
.lhs_from_mass_matrix <- function(D, M, vars) {
  n <- length(vars)
  if (is.null(M) || (is.matrix(M) && all(M == diag(n)))) return(lapply(vars, D))
  M <- as.matrix(M)
  lapply(seq_len(nrow(M)), function(i) {
    nz <- which(M[i, ] != 0)
    if (!length(nz)) return(.S0())
    if (length(nz) == 1L && M[i, nz] == 1) return(D(vars[[nz]]))
    stop("Non-permutation mass matrix is not supported.", call. = FALSE)
  })
}

.mtkize_defaults <- function(vars, u0, pr) {
  m <- symmap()
  for (i in seq_along(vars)) m <- smap_set(m, vars[[i]], unname(u0[i]))
  for (i in seq_along(pr$vars)) m <- smap_set(m, pr$vars[[i]], pr$values[i])
  m
}

.mtkize_name <- function(prefix, dots) dots$name %||% paste0(prefix, .next_id())

#' Build a symbolic system from a numerical problem
#'
#' `modelingtoolkitize(prob)` traces the function of a numerical problem
#' with symbolic arguments and returns the equivalent `System`, so the
#' symbolic tools (Jacobians, simplification, index reduction) apply to code
#' written for the numerical solvers. Methods exist for `ODEProblem` (`f(u,
#' p, t)`), `NonlinearProblem` and `NonlinearLeastSquaresProblem` (`f(u,
#' p)`), `OptimizationProblem` (objective `f(u, p)`, constraints, bounds) and
#' `SDEProblem` (drift and diffusion).
#'
#' The function is called with symbolic vectors for `u` and `p` (and the
#' time `t`): indexing (`u[1]`, `p["k"]`, `p$k`), arithmetic, math functions,
#' `sum()`, `prod()` and `c()` work as for numeric vectors, and the function
#' must return its result (a vector built with `c()` or a list).
#' Branching on symbolic values cannot be traced.
#'
#' Unknowns are named from `u_names`, the system the problem was built from,
#' the names of `u0`, or \ifelse{latex}{\verb{x_1, x_2, ...}}{\verb{x₁, x₂, ...}}; parameters from `p_names`, the
#' system, the names of `p`, or \ifelse{latex}{\verb{alpha_1, alpha_2, ...}}{\verb{α₁, α₂, ...}} (\ifelse{latex}{\verb{alpha}}{\verb{α}} for a single number). The
#' values of `u0` and `p` become initial conditions of the system.
#'
#' @param prob A problem.
#' @param u_names,p_names Names of the unknowns and parameters.
#' @param ... Passed to [System()], for example `name`.
#' @return A `System`.
#' @examples
#' lotka <- function(u, p, t) {
#'   x <- u[1]; y <- u[2]
#'   c(p[1] * x - p[2] * x * y, -p[3] * y + p[4] * x * y)
#' }
#' prob <- ODEProblem(lotka, c(1, 1), c(0, 10), c(1.5, 1, 3, 1))
#' sys <- modelingtoolkitize(prob, name = "lv")
#' equations(sys)
#' @export
modelingtoolkitize <- function(prob, ...) UseMethod("modelingtoolkitize")

#' @rdname modelingtoolkitize
#' @param return_symbolic_u0_p Internal: also return the traced variables.
#' @export
modelingtoolkitize.ODEProblem <- function(prob, u_names = NULL, p_names = NULL, ...,
                                          return_symbolic_u0_p = FALSE) {
  dots <- list(...)
  t <- t_nounits
  vars <- .mtkize_vars(prob, t, u_names)
  pr <- .mtkize_params(prob, p_names)
  lhs <- .lhs_from_mass_matrix(D_nounits, prob$f$mass_matrix, vars)
  rhs <- .trace(.mtkize_fun(prob$f), list(.symvec(vars, names(prob$u0)), pr$obj, .symvec(list(t))))
  if (length(rhs) != length(vars))
    stop("the traced function returned ", length(rhs), " values for ", length(vars), " unknowns", call. = FALSE)
  eqs <- lapply(seq_along(vars), function(i) Equation(lhs[[i]], rhs[[i]]))
  dots$name <- .mtkize_name("MTKizedODE", dots)
  sys <- do.call(System, c(list(eqs, t, vars, pr$vars,
                                initial_conditions = .mtkize_defaults(vars, prob$u0, pr)), dots))
  if (return_symbolic_u0_p) return(list(sys = sys, vars = vars, params = pr))
  sys
}

#' @rdname modelingtoolkitize
#' @export
modelingtoolkitize.AbstractNonlinearProblem <- function(prob, u_names = NULL, p_names = NULL, ...) {
  dots <- list(...)
  vars <- .mtkize_vars(prob, NULL, u_names)
  pr <- .mtkize_params(prob, p_names)
  rhs <- .trace(.mtkize_fun(prob$f), list(.symvec(vars, names(prob$u0)), pr$obj))
  eqs <- lapply(Filter(function(r) !.is_zero(r), rhs), function(r) Equation(.S0(), r))
  dots$name <- .mtkize_name("MTKizedNonlin", dots)
  do.call(System, c(list(eqs, NULL, vars, pr$vars,
                         initial_conditions = .mtkize_defaults(vars, prob$u0, pr)), dots))
}

#' @rdname modelingtoolkitize
#' @export
modelingtoolkitize.NonlinearProblem <- function(prob, u_names = NULL, p_names = NULL, ...)
  modelingtoolkitize.AbstractNonlinearProblem(prob, u_names, p_names, ...)

#' @rdname modelingtoolkitize
#' @export
modelingtoolkitize.OptimizationProblem <- function(prob, u_names = NULL, p_names = NULL, ...) {
  dots <- list(...)
  vars <- .mtkize_vars(prob, NULL, u_names)
  lb <- prob$lb; ub <- prob$ub
  if (!is.null(ub)) {
    vars <- lapply(seq_along(vars), function(i) {
      lo <- lb[i]; hi <- ub[i]
      if ((lo == 0 && hi == 0) || (is.infinite(lo) && lo < 0 && is.infinite(hi) && hi > 0)) vars[[i]]
      else setmetadata(vars[[i]], "bounds", c(lo, hi))
    })
  }
  pr <- .mtkize_params(prob, p_names)
  u <- .symvec(vars, names(prob$u0))
  fobj <- .mtkize_fun(prob$f)
  objective <- .trace(fobj, list(u, pr$obj))
  if (length(objective) != 1L) stop("the objective must return one value", call. = FALSE)
  consf <- if (is.list(prob$f)) prob$f$cons else prob$cons
  cons <- list()
  if (!is.null(consf)) {
    lhs <- .trace(consf, list(u, pr$obj))
    lc <- prob$lcons; uc <- prob$ucons
    if ((is.null(lc) || all(is.infinite(lc))) && (is.null(uc) || all(is.infinite(uc))))
      stop("Constraints passed have no proper bounds defined. Ensure you pass equal bounds (the scalar ",
           "that the constraint should evaluate to) for equality constraints or pass the lower and upper ",
           "bounds for inequality constraints.", call. = FALSE)
    for (i in seq_along(lhs)) {
      if (!is.null(lc) && !is.infinite(lc[i])) {
        cons[[length(cons) + 1L]] <- if (!is.null(uc) && lc[i] == uc[i]) Equation(lhs[[i]], uc[i])
                                     else Inequality(lc[i], lhs[[i]], "<=")
      }
    }
    for (i in seq_along(lhs)) {
      if (!is.null(uc) && !is.infinite(uc[i]) && (is.null(lc) || lc[i] != uc[i]))
        cons[[length(cons) + 1L]] <- Inequality(lhs[[i]], uc[i], "<=")
    }
  }
  dots$name <- .mtkize_name("MTKizedOpt", dots)
  do.call(System, c(list(list(), NULL, vars, pr$vars, costs = objective, constraints = cons,
                         initial_conditions = .mtkize_defaults(vars, prob$u0, pr)), dots))
}

#' @rdname modelingtoolkitize
#' @export
modelingtoolkitize.SDEProblem <- function(prob, u_names = NULL, p_names = NULL, ...) {
  dots <- list(...)
  dots$name <- .mtkize_name("MTKizedSDE", dots)
  ode <- prob
  class(ode) <- c("ODEProblem", "SciMLProblem")
  res <- do.call(modelingtoolkitize.ODEProblem,
                 c(list(ode, u_names = u_names, p_names = p_names, return_symbolic_u0_p = TRUE), dots))
  sys <- res$sys
  args <- list(.symvec(res$vars, names(prob$u0)), res$params$obj, .symvec(list(t_nounits)))
  neqs <- .trace(prob$f$g, args)
  n <- length(res$vars)
  type <- prob$f$noise_type %||% prob$noise_type %||% "diagonal"
  if (type == "diagonal") {
    if (length(neqs) != n) stop("diagonal noise needs one diffusion term per unknown", call. = FALSE)
  } else {
    neqs <- matrix(neqs, nrow = n)
  }
  .update_system(sys, noise_eqs = neqs)
}
