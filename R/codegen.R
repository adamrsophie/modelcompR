# Code generation: symbolic expressions to R closures, after
# ModelingToolkitBase/src/systems/codegen.jl and Symbolics' build_function.
#
# Variables are first replaced by placeholder symbols (".u1", ".p2", ...),
# common subexpressions are hoisted with symcompR::cse(), and each expression
# tree is translated into an R call. Placeholders map back to vector slots
# such as u[1] when the call is built.

#' Register a numeric R function for use in symbolic expressions
#'
#' `register_symbolic("f", fn)` returns a symbolic function `f` whose
#' applications `f(x)` stay unevaluated in expressions and call `fn` in
#' generated code, like `@register_symbolic` in ModelingToolkit. Derivatives
#' are taken by central differences unless `derivatives` gives a function for
#' each argument.
#'
#' @param name Function name.
#' @param fn An R function of numeric arguments.
#' @param derivatives Optional list of R functions, the partial derivatives
#'   with respect to each argument, taking the same arguments as `fn`.
#' @return A symbolic function.
#' @examples
#' sq <- register_symbolic("sq", function(x) x^2)
#' x <- Variable("x", t_nounits)
#' sq(x)
#' @export
register_symbolic <- function(name, fn, derivatives = NULL) {
  stopifnot(is.function(fn))
  assign(name, list(fn = fn, derivatives = derivatives), envir = .mtk$funs)
  Function(name)
}

.is_registered_fun <- function(e) .type(e) == "Function" && exists(.name(e), envir = .mtk$funs, inherits = FALSE)

# Numeric partial derivative of a registered function.
.fd_partial <- function(fn, i, args) {
  h <- 1e-6 * max(1, abs(args[[i]]))
  a1 <- args; a2 <- args
  a1[[i]] <- args[[i]] + h
  a2[[i]] <- args[[i]] - h
  (do.call(fn, a1) - do.call(fn, a2)) / (2 * h)
}

.R_MATH <- c(sin = "sin", cos = "cos", tan = "tan", asin = "asin", acos = "acos", atan = "atan",
             sinh = "sinh", cosh = "cosh", tanh = "tanh", asinh = "asinh", acosh = "acosh",
             atanh = "atanh", exp = "exp", Abs = "abs", sign = "sign", floor = "floor",
             ceiling = "ceiling", gamma = "gamma", loggamma = "lgamma", factorial = "factorial",
             beta = "beta", atan2 = "atan2", re = "Re", im = "Im", sqrt = "sqrt")

.cg_error <- function(e) stop("cannot generate code for ", format(e), call. = FALSE)

# Translate a lowered expression into an R call. `names` maps placeholder
# symbol names to R language objects.
.to_rcall <- function(e, names) {
  ty <- .type(e)
  switch(ty,
    Integer = , Rational = as.double(as.numeric(e)),
    Float = as.numeric(e),
    Symbol = , Dummy = {
      nm <- .name(e)
      hit <- names[[nm]]
      if (!is.null(hit)) return(hit)
      stop("symbol ", nm, " is not an argument of the generated function", call. = FALSE)
    },
    Constant = switch(EXPR = .name(e), pi = quote(pi), E = exp(1), oo = Inf, `-oo` = -Inf, nan = NaN,
                      zoo = Inf, I = 1i, as.numeric(N(e))),
    Add = Reduce(function(a, b) call("+", a, b), lapply(.args(e), .to_rcall, names = names)),
    Mul = {
      a <- .args(e)
      if (.is_num(a[[1L]]) && as.numeric(a[[1L]]) == -1 && length(a) >= 2L)
        return(call("-", Reduce(function(x, y) call("*", x, y), lapply(a[-1L], .to_rcall, names = names))))
      Reduce(function(x, y) call("*", x, y), lapply(a, .to_rcall, names = names))
    },
    Pow = {
      b <- .args(e)[[1L]]; x <- .args(e)[[2L]]
      if (.is_num(x)) {
        xv <- as.numeric(x)
        if (xv == 0.5) return(call("sqrt", .to_rcall(b, names)))
        if (xv == -1) return(call("/", 1, .to_rcall(b, names)))
        if (xv == 2) return(call("^", .to_rcall(b, names), 2))
      }
      if (.type(b) == "Constant" && identical(.name(b), "E")) return(call("exp", .to_rcall(x, names)))
      call("^", .to_rcall(b, names), .to_rcall(x, names))
    },
    Function = .fun_rcall(e, names),
    Piecewise = {
      pairs <- .args(e)
      out <- NaN
      for (pr in rev(pairs)) {
        val <- .to_rcall(.args(pr)[[1L]], names)
        cond <- .args(pr)[[2L]]
        if (.type(cond) == "BooleanAtom" && identical(.name(cond), "true")) out <- val
        else out <- call("if", call("isTRUE", .to_rcall(cond, names)), val, out)
      }
      out
    },
    BooleanAtom = identical(.name(e), "true"),
    Equality = call("==", .to_rcall(.args(e)[[1L]], names), .to_rcall(.args(e)[[2L]], names)),
    Unequality = call("!=", .to_rcall(.args(e)[[1L]], names), .to_rcall(.args(e)[[2L]], names)),
    StrictLessThan = call("<", .to_rcall(.args(e)[[1L]], names), .to_rcall(.args(e)[[2L]], names)),
    LessThan = call("<=", .to_rcall(.args(e)[[1L]], names), .to_rcall(.args(e)[[2L]], names)),
    StrictGreaterThan = call(">", .to_rcall(.args(e)[[1L]], names), .to_rcall(.args(e)[[2L]], names)),
    GreaterThan = call(">=", .to_rcall(.args(e)[[1L]], names), .to_rcall(.args(e)[[2L]], names)),
    And = Reduce(function(a, b) call("&&", a, b), lapply(.args(e), .to_rcall, names = names)),
    Or = Reduce(function(a, b) call("||", a, b), lapply(.args(e), .to_rcall, names = names)),
    Not = call("!", .to_rcall(.args(e)[[1L]], names)),
    Derivative = .registered_derivative_rcall(e, names),
    Subs = .subs_rcall(e, names),
    .cg_error(e))
}

.fun_rcall <- function(e, names) {
  nm <- .name(e)
  args <- lapply(.args(e), .to_rcall, names = names)
  if (exists(nm, envir = .mtk$funs, inherits = FALSE)) {
    fn <- get(nm, envir = .mtk$funs)$fn
    return(as.call(c(list(fn), args)))
  }
  r <- .R_MATH[nm]
  if (!is.na(r)) return(as.call(c(list(as.name(r)), args)))
  switch(nm,
    log = if (length(args) == 1L) call("log", args[[1L]]) else call("log", args[[1L]], args[[2L]]),
    Max = as.call(c(list(as.name("max")), args)),
    Min = as.call(c(list(as.name("min")), args)),
    sec = call("/", 1, call("cos", args[[1L]])),
    csc = call("/", 1, call("sin", args[[1L]])),
    cot = call("/", 1, call("tan", args[[1L]])),
    Heaviside = call("+", call("as.numeric", call(">", args[[1L]], 0)),
                     call("*", if (length(args) > 1L) args[[2L]] else 0.5, call("as.numeric", call("==", args[[1L]], 0)))),
    erf = call("-", call("*", 2, call("pnorm", call("*", args[[1L]], sqrt(2)))), 1),
    erfc = call("*", 2, call("pnorm", call("*", call("-", args[[1L]]), sqrt(2)))),
    Mod = call("%%", args[[1L]], args[[2L]]),
    DiracDelta = 0,
    sinc = call("ifelse", call("==", args[[1L]], 0), 1, call("/", call("sin", args[[1L]]), args[[1L]])),
    .cg_error(e))
}

# Derivative(f(x), x) of a registered function f, evaluated numerically or
# with the user's derivative.
.registered_derivative_rcall <- function(e, names) {
  inner <- .args(e)[[1L]]
  if (.type(inner) != "Function" || !exists(.name(inner), envir = .mtk$funs, inherits = FALSE))
    .cg_error(e)
  reg <- get(.name(inner), envir = .mtk$funs)
  out <- NULL
  fargs <- .args(inner)
  for (spec in .args(e)[-1L]) {
    v <- .args(spec)[[1L]]
    n <- as.integer(as.numeric(.args(spec)[[2L]]))
    i <- which(vapply(fargs, function(a) identical(.key(a), .key(v)), TRUE))
    if (length(i) != 1L || n != 1L || length(.args(e)) > 2L) .cg_error(e)
    argcalls <- lapply(fargs, .to_rcall, names = names)
    if (!is.null(reg$derivatives) && length(reg$derivatives) >= i) {
      out <- as.call(c(list(reg$derivatives[[i]]), argcalls))
    } else {
      out <- call(".fd_partial", reg$fn, i, as.call(c(list(as.name("list")), argcalls)))
    }
  }
  out
}

# Subs(Derivative(f(xi), xi), xi, point): the derivative evaluated at point.
.subs_rcall <- function(e, names) {
  a <- .args(e)
  expr <- a[[1L]]
  vars <- if (.type(a[[2L]]) == "Tuple") .args(a[[2L]]) else list(a[[2L]])
  pts <- if (.type(a[[3L]]) == "Tuple") .args(a[[3L]]) else list(a[[3L]])
  inner_names <- names
  bind <- list()
  for (i in seq_along(vars)) {
    nm <- paste0(".s", .next_id())
    inner_names[[.name(vars[[i]])]] <- as.name(nm)
    bind[[nm]] <- .to_rcall(pts[[i]], names)
  }
  body <- .to_rcall(expr, inner_names)
  as.call(c(list(as.name("local")), list(as.call(c(list(as.name("{")),
    lapply(names(bind), function(n) call("<-", as.name(n), bind[[n]])), list(body))))))
}

# ---------------------------------------------------------------------------
# Lowering: variables to placeholder symbols

# Build the lowering for argument groups. `args` is a list of argument
# specifications, each list(name = "u", vars = list(<sym>...)) for a vector
# argument or list(name = "t", vars = <sym>) for a scalar. Returns the
# substitution rule (sym -> placeholder) and the placeholder name map.
.lowering <- function(args) {
  from <- list(); to <- list(); names <- list()
  for (a in args) {
    if (is.null(a$vars)) next
    if (is_sym(a$vars)) {
      ph <- paste0(".", a$name)
      from[[length(from) + 1L]] <- a$vars
      to[[length(to) + 1L]] <- Symbol(ph)
      names[[ph]] <- as.name(a$name)
      next
    }
    for (i in seq_along(a$vars)) {
      ph <- paste0(".", a$name, i)
      from[[length(from) + 1L]] <- a$vars[[i]]
      to[[length(to) + 1L]] <- Symbol(ph)
      names[[ph]] <- call("[", as.name(a$name), i)
    }
  }
  list(rule = .rule(from, to), names = names)
}

# Observed equations needed by `exprs`, in dependency order.
.needed_observed <- function(exprs, obs) {
  if (!length(obs)) return(list())
  lhs_keys <- vapply(obs, function(o) .key(o$lhs), "")
  need <- rep(FALSE, length(obs))
  frontier <- .keys(.collect_vars(exprs))
  repeat {
    hit <- which(lhs_keys %in% frontier & !need)
    if (!length(hit)) break
    need[hit] <- TRUE
    frontier <- .keys(.collect_vars(lapply(obs[hit], `[[`, "rhs")))
  }
  .toposort_observed(obs[need])
}

# Sort observed equations so each lhs is defined before it is used.
.toposort_observed <- function(obs) {
  if (length(obs) <= 1L) return(obs)
  lhs <- vapply(obs, function(o) .key(o$lhs), "")
  deps <- lapply(obs, function(o) intersect(.keys(.collect_vars(o$rhs)), lhs))
  done <- character()
  out <- list()
  remaining <- seq_along(obs)
  while (length(remaining)) {
    ready <- remaining[vapply(remaining, function(i) all(deps[[i]] %in% c(done, lhs[i])), TRUE)]
    if (!length(ready)) stop("observed equations are cyclic", call. = FALSE)
    for (i in ready) {
      out[[length(out) + 1L]] <- obs[[i]]
      done <- c(done, lhs[i])
    }
    remaining <- setdiff(remaining, ready)
  }
  out
}

# c() for generated functions: the names go from plain numeric results only, so
# symbolic (modelingtoolkitize) and dual number inputs keep their own c methods
.cg_c <- function(...) {
  r <- c(...)
  if (is.numeric(r) && !is.object(r)) names(r) <- NULL
  r
}

#' Build an R function from symbolic expressions
#'
#' `build_function(exprs, ...)` returns an R function whose arguments are the
#' groups in `...`: each group is a list of variables (passed as one numeric
#' vector) or a single variable (passed as a number). Observed equations, if
#' given, are computed first, so `exprs` may refer to observed variables.
#' This plays the role of Symbolics' `build_function`.
#'
#' @param exprs An expression, a list of expressions, or a list matrix.
#' @param ... Argument groups. Named arguments name the function arguments;
#'   unnamed ones are called `u`, `p` and `t` in order, or `a1`, `a2`, ...
#' @param observed Observed equations to compute first.
#' @param cse Use common subexpression elimination.
#' @param expression Return the function's R code instead of the function.
#' @return An R function, or a language object when `expression = TRUE`.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y"); a <- Symbol("a")
#' f <- build_function(list(x * a, sin(y)), list(x, y), a)
#' f(c(1, 2), 3)
#' @export
build_function <- function(exprs, ..., observed = list(), cse = TRUE, expression = FALSE) {
  groups <- list(...)
  default_names <- c("u", "p", "t")
  nms <- names(groups) %||% rep("", length(groups))
  for (i in seq_along(groups)) {
    if (!nzchar(nms[i])) nms[i] <- if (length(groups) <= 3L) default_names[i] else paste0("a", i)
  }
  args <- lapply(seq_along(groups), function(i) {
    g <- groups[[i]]
    list(name = nms[i], vars = if (is.null(g)) NULL else if (is_sym(g)) g else .flatten_vars(g))
  })
  .build_closure(exprs, args, observed = observed, cse = cse, expression = expression)
}

# Core builder used by every generator. `args` as in .lowering(); `exprs` a
# sym, a list of syms (vector output) or a list matrix (matrix output).
.build_closure <- function(exprs, args, observed = list(), cse = TRUE, expression = FALSE,
                           out = c("auto", "vector", "matrix", "scalar"), sparse_matrix = FALSE) {
  out <- match.arg(out)
  is_mat <- is.matrix(exprs) || inherits(exprs, "symmatrix")
  if (inherits(exprs, "symmatrix")) exprs <- .symmatrix_to_list_matrix(exprs)
  dims <- if (is_mat) dim(exprs) else NULL
  scalar <- is_sym(exprs) || is.numeric(exprs)
  flat <- if (scalar) list(as_sym(exprs)) else lapply(unclass(exprs), as_sym)
  if (out == "auto") out <- if (is_mat) "matrix" else if (scalar) "scalar" else "vector"

  low <- .lowering(args)
  names <- low$names
  obs <- .needed_observed(flat, observed)
  obs_lines <- list()
  rule <- low$rule
  if (length(obs)) {
    for (i in seq_along(obs)) {
      ph <- paste0(".o", i)
      rhs <- .subst(obs[[i]]$rhs, rule)
      obs_lines[[i]] <- list(ph = ph, rhs = rhs)
      rule <- .rule_merge(rule, .rule(list(obs[[i]]$lhs), list(Symbol(ph))))
      names[[ph]] <- as.name(ph)
    }
  }
  lowered <- lapply(flat, function(e) .subst(e, rule))

  stmts <- list()
  for (ol in obs_lines) stmts[[length(stmts) + 1L]] <- call("<-", as.name(ol$ph), .to_rcall(ol$rhs, names))
  nz <- !vapply(lowered, .is_zero, TRUE)
  if (cse && sum(nz) > 0L && !any(vapply(lowered[nz], .has_piecewise, TRUE))) {
    cs <- symcompR::cse(lowered[nz], prefix = ".c")
    for (r in cs$replacements) {
      ph <- .name(r[[1L]])
      names[[ph]] <- as.name(ph)
      stmts[[length(stmts) + 1L]] <- call("<-", as.name(ph), .to_rcall(r[[2L]], names))
    }
    lowered[nz] <- cs$reduced
  }
  calls <- lapply(lowered, .to_rcall, names = names)
  result <- switch(out,
    # .cg_c drops names: arithmetic on a named u would otherwise label the
    # entries with the names of whatever operand came first
    scalar = if (is.numeric(calls[[1L]])) calls[[1L]] else call(".cg_c", calls[[1L]]),
    vector = as.call(c(list(as.name(".cg_c")), calls)),
    matrix = if (sparse_matrix) {
      # the structural nonzeros as a Matrix::sparseMatrix, as the SparseMatrixCSC
      # that Symbolics builds for sparse = true
      k <- which(nz)
      as.call(list(quote(Matrix::sparseMatrix), i = as.integer((k - 1L) %% dims[1L] + 1L),
           j = as.integer((k - 1L) %/% dims[1L] + 1L),
           x = as.call(c(list(as.name(".cg_c")), calls[k])),
           dims = as.integer(dims)))
    } else {
      body <- list(call("<-", as.name(".J"), call("matrix", 0, dims[1L], dims[2L])))
      for (k in which(nz)) {
        i <- ((k - 1L) %% dims[1L]) + 1L
        j <- ((k - 1L) %/% dims[1L]) + 1L
        body[[length(body) + 1L]] <- call("<-", call("[", as.name(".J"), i, j), calls[[k]])
      }
      body[[length(body) + 1L]] <- as.name(".J")
      as.call(c(list(as.name("{")), body))
    })
  if (out == "vector" && !length(calls)) result <- quote(numeric(0))
  body <- as.call(c(list(as.name("{")), stmts, list(result)))
  fargs <- stats::setNames(rep(list(quote(expr = )), length(args)), vapply(args, `[[`, "", "name"))
  fn <- as.function(c(fargs, body))
  environment(fn) <- asNamespace("modelcompR")
  if (expression) return(body)
  fn
}

.has_piecewise <- function(e) {
  if (.type(e) == "Piecewise") return(TRUE)
  for (a in .args(e)) if (.has_piecewise(a)) return(TRUE)
  FALSE
}

.symmatrix_to_list_matrix <- function(M) {
  d <- dim(M)
  out <- matrix(list(), d[1L], d[2L])
  for (i in seq_len(d[1L])) for (j in seq_len(d[2L])) out[[i, j]] <- M[i, j]
  out
}

# ---------------------------------------------------------------------------
# Symbolic derivatives of systems

.list_matrix <- function(rows, cols, f) {
  out <- matrix(list(), rows, cols)
  for (i in seq_len(rows)) for (j in seq_len(cols)) out[[i, j]] <- f(i, j)
  out
}

# Differentiate expressions holding the given variables independent: each
# variable is replaced by a placeholder symbol first, so d/dt of t*x(t) is x(t)
# rather than the total derivative.
.partials <- function(exprs, wrt, others = list()) {
  hold <- c(wrt, others)
  ph <- lapply(seq_along(hold), function(i) Symbol(paste0(".h", i)))
  rule <- .rule(hold, ph)
  back <- .rule(ph, hold)
  low <- lapply(exprs, function(e) .subst(e, rule))
  wph <- ph[seq_along(wrt)]
  wkeys <- .keys(wph)
  # differentiate each expression only by the variables it contains
  present <- lapply(low, function(e) intersect(.keys(free_symbols(e)), wkeys))
  zero <- .S0()
  .list_matrix(length(exprs), length(wrt), function(i, j) {
    if (!(wkeys[j] %in% present[[i]])) return(zero)
    .subst(diff(low[[i]], wph[[j]]), back)
  })
}

.as_symmatrix <- function(lm) {
  rows <- lapply(seq_len(nrow(lm)), function(i) lapply(seq_len(ncol(lm)), function(j) lm[[i, j]]))
  if (!length(rows)) return(lm)
  Matrix(rows)
}

# Right-hand sides of a compiled ODE system, observed variables and bound
# parameters substituted: for D(x) ~ f equations, f; for algebraic 0 ~ g
# equations, g.
.ode_rhs_exprs <- function(sys, substitute_observed = TRUE) {
  eqs <- equations(sys)
  obs <- if (iscomplete(sys)) .codegen_observed(sys) else observed(sys)
  rule <- if (substitute_observed && length(obs))
    .rule(lapply(obs, `[[`, "lhs"), lapply(obs, `[[`, "rhs")) else NULL
  lapply(eqs, function(e) {
    r <- if (.is_derivative(e$lhs)) e$rhs else e$rhs - e$lhs
    if (!is.null(rule)) .subst_fixpoint(r, rule) else r
  })
}

#' Symbolic Jacobians and related derivatives
#'
#' `calculate_jacobian(sys)` is the Jacobian of the right-hand side of a
#' compiled system with respect to its unknowns, `calculate_tgrad(sys)` the
#' partial derivative with respect to time, `calculate_paramjac(sys)` the
#' Jacobian with respect to the parameters, `calculate_massmatrix(sys)` the
#' mass matrix, and `calculate_hessian(sys)` the Hessians of the
#' (scalar) cost of an optimization system or of each equation of a
#' nonlinear system. The `generate_*` functions return R functions of
#' `(u, p, t)` (or `(u, p)` for time independent systems).
#'
#' @param sys A complete `System`.
#' @param sparse For `generate_jacobian()`, return a function that builds a
#'   `Matrix::sparseMatrix` of the structural nonzeros. Symbolic matrices
#'   are always stored dense, so `calculate_jacobian()` ignores it.
#' @param simplify Simplify each entry.
#' @param ... Unused.
#' @return A symbolic matrix, or an R function.
#' @export
calculate_jacobian <- function(sys, sparse = FALSE, simplify = FALSE, ...) {
  rhs <- .ode_rhs_exprs(sys)
  J <- .partials(rhs, unknowns(sys), c(.parameters_of(sys, TRUE), independent_variables(sys)))
  if (simplify) J[] <- lapply(J, symcompR::simplify)
  .as_symmatrix(J)
}

#' @rdname calculate_jacobian
#' @export
calculate_tgrad <- function(sys, simplify = FALSE, ...) {
  iv <- .f(sys, "iv")
  if (is.null(iv)) stop("a time independent system has no time gradient", call. = FALSE)
  rhs <- .ode_rhs_exprs(sys)
  J <- .partials(rhs, list(iv), c(unknowns(sys), .parameters_of(sys, TRUE)))
  out <- lapply(seq_len(nrow(J)), function(i) J[[i, 1L]])
  if (simplify) out <- lapply(out, symcompR::simplify)
  out
}

#' @rdname calculate_jacobian
#' @export
calculate_paramjac <- function(sys, simplify = FALSE, ...) {
  rhs <- .ode_rhs_exprs(sys)
  J <- .partials(rhs, .problem_params(sys), c(unknowns(sys), independent_variables(sys)))
  if (simplify) J[] <- lapply(J, symcompR::simplify)
  .as_symmatrix(J)
}

#' @rdname calculate_jacobian
#' @export
calculate_massmatrix <- function(sys, ...) {
  eqs <- equations(sys)
  us <- unknowns(sys)
  n <- length(us)
  M <- diag(0, n)
  for (i in seq_along(eqs)) {
    lhs <- eqs[[i]]$lhs
    if (.is_derivative(lhs)) {
      v <- .args(lhs)[[1L]]
      j <- match(.key(v), .keys(us))
      if (is.na(j)) stop("derivative of a non-unknown in ", .fmt(eqs[[i]]), call. = FALSE)
      M[i, j] <- 1
    }
  }
  if (n && all(M == diag(n))) return(diag(n))
  M
}

#' @rdname calculate_jacobian
#' @export
calculate_hessian <- function(sys, ...) {
  vars <- unknowns(sys)
  others <- c(.parameters_of(sys, TRUE), independent_variables(sys))
  if (length(costs(sys))) {
    g <- .partials(list(cost(sys)), vars, others)
    grad <- lapply(seq_along(vars), function(j) g[[1L, j]])
    return(.as_symmatrix(.partials(grad, vars, others)))
  }
  rhs <- .ode_rhs_exprs(sys)
  lapply(rhs, function(r) {
    g <- .partials(list(r), vars, others)
    grad <- lapply(seq_along(vars), function(j) g[[1L, j]])
    .as_symmatrix(.partials(grad, vars, others))
  })
}

# Arguments of generated functions: the unknowns, the problem parameter
# vector (bound parameters excluded, see .problem_params) and time.
.sys_args <- function(sys, with_t = TRUE) .codegen_args(sys, with_t = with_t)

# Assertions of the system are added to the last right-hand side, so a failing
# assertion makes it NaN (hook of the debugging module; see debugging-*.R).
.add_assertions <- function(sys, rhs) {
  h <- .mtk$hooks$rhs_assertions
  if (is.null(h) || !length(rhs) || !length(assertions(sys))) return(rhs)
  h(sys, rhs)
}

#' @rdname calculate_jacobian
#' @export
generate_rhs <- function(sys, ...) {
  rhs <- .add_assertions(sys, .ode_rhs_exprs(sys, substitute_observed = FALSE))
  .build_closure(rhs, .sys_args(sys), observed = .codegen_observed(sys), out = "vector", ...)
}

#' @rdname calculate_jacobian
#' @export
generate_jacobian <- function(sys, sparse = FALSE, ...) {
  J <- calculate_jacobian(sys)
  .build_closure(.symmatrix_to_list_matrix(J), .sys_args(sys), out = "matrix", sparse_matrix = sparse, ...)
}

#' @rdname calculate_jacobian
#' @export
generate_tgrad <- function(sys, ...) {
  .build_closure(calculate_tgrad(sys), .sys_args(sys), out = "vector", ...)
}

#' @rdname calculate_jacobian
#' @export
generate_paramjac <- function(sys, ...) {
  J <- calculate_paramjac(sys)
  .build_closure(.symmatrix_to_list_matrix(J), .sys_args(sys), out = "matrix", ...)
}

#' @rdname calculate_jacobian
#' @param exprs Expressions of the system's variables.
#' @export
build_explicit_observed_function <- function(sys, exprs, ...) {
  single <- is_sym(exprs) || is.numeric(exprs)
  exprs <- if (single) list(as_sym(exprs)) else lapply(.flatten_vars(exprs), as_sym)
  fn <- .build_closure(exprs, .sys_args(sys, !is.null(.f(sys, "iv"))), observed = .codegen_observed(sys),
                       out = if (single) "scalar" else "vector", ...)
  fn
}

#' @rdname calculate_jacobian
#' @param expr An expression to compile.
#' @export
generate_custom_function <- function(sys, expr, ...) build_explicit_observed_function(sys, expr, ...)

# ---------------------------------------------------------------------------
# Index cache: positions of unknowns and parameters in the numeric vectors.

.build_index_cache <- function(sys) {
  us <- unknowns(sys)
  # positions in prob$p: bound parameters are not stored (see .problem_params)
  ps <- .problem_params(sys)
  list(unknown_index = stats::setNames(seq_along(us), .keys(us)),
       param_index = stats::setNames(seq_along(ps), .keys(ps)),
       observed_index = stats::setNames(seq_along(observed(sys)),
                                        vapply(observed(sys), function(o) .key(o$lhs), "")))
}

#' Symbolic if-then-else
#'
#' `IfElse(cond, a, b)` is `a` when `cond` holds and `b` otherwise, as
#' `ifelse` in ModelingToolkit models. It is a symcompR `Piecewise`.
#'
#' @param cond A symbolic condition, such as `x > 0`.
#' @param a,b Expressions.
#' @return A `sym` expression.
#' @export
IfElse <- function(cond, a, b) {
  if (isTRUE(cond)) return(as_sym(a))
  if (identical(cond, FALSE)) return(as_sym(b))
  Piecewise(list(as_sym(a), cond), list(as_sym(b), TRUE))
}
