# Internal helpers shared by every module: the variable registry, symbolic
# maps keyed by symcompR keys, and small expression utilities.
#
# symcompR expressions are immutable lists of class "sym" whose `key` field is
# a canonical string for the whole tree. Every map in modelcompR is keyed by
# that string, so two structurally equal expressions share one entry.

.mtk <- new.env(parent = emptyenv())
.mtk$vars <- new.env(parent = emptyenv())   # name -> list(kind, meta)
.mtk$funs <- new.env(parent = emptyenv())   # registered symbolic functions
.mtk$counter <- 0L
.mtk$hooks <- new.env(parent = emptyenv())  # extension points between modules

# The namespace separator used by ModelingToolkit (U+208A, subscript plus).
NAMESPACE_SEPARATOR <- "\u208A"

`%||%` <- function(a, b) if (is.null(a)) b else a

.key <- function(e) .subset2(e, "key")
.type <- function(e) .subset2(e, "type")
.args <- function(e) .subset2(e, "args")
.name <- function(e) .subset2(e, "name")

.next_id <- function() {
  .mtk$counter <- .mtk$counter + 1L
  .mtk$counter
}

# ---------------------------------------------------------------------------
# Symbolic maps: named lists whose names are keys and whose elements are
# list(var = <sym>, value = <anything>). They keep insertion order, which the
# problem constructors rely on for deterministic output.

symmap <- function() structure(list(), class = "mtk_symmap")

smap_set <- function(m, var, value) {
  var <- as_sym(var)
  m[[.key(var)]] <- list(var = var, value = value)
  if (!inherits(m, "mtk_symmap")) class(m) <- "mtk_symmap"
  m
}

smap_get <- function(m, var, default = NULL) {
  hit <- m[[.key(as_sym(var))]]
  if (is.null(hit)) default else hit$value
}

smap_has <- function(m, var) !is.null(m[[.key(as_sym(var))]])

smap_del <- function(m, var) {
  m[[.key(as_sym(var))]] <- NULL
  m
}

smap_vars <- function(m) lapply(unclass(m), `[[`, "var")
smap_values <- function(m) lapply(unclass(m), `[[`, "value")

# Merge maps; later maps win.
smap_merge <- function(...) {
  out <- symmap()
  for (m in list(...)) {
    if (is.null(m)) next
    m <- as_symmap(m)
    for (k in names(m)) out[[k]] <- m[[k]]
  }
  class(out) <- "mtk_symmap"
  out
}

# Coerce user input into a symmap. Accepted forms: a symmap; a list of
# Equation objects or formulas (`x ~ 1`); a list of pairs list(var, value);
# a named list or named numeric vector whose names are variable names (looked
# up later by name, see .resolve_named()).
as_symmap <- function(x) {
  if (is.null(x)) return(symmap())
  if (inherits(x, "mtk_symmap")) return(x)
  out <- symmap()
  if (is.numeric(x) && !is.null(names(x))) x <- as.list(x)
  if (inherits(x, "Equation") || inherits(x, "formula")) x <- list(x)
  if (!is.list(x)) stop("cannot interpret ", class(x)[1], " as a variable map", call. = FALSE)
  nms <- names(x)
  for (i in seq_along(x)) {
    el <- x[[i]]
    if (inherits(el, "formula") && length(el) == 3L) {
      # evaluate both sides here so `p ~ NA` (a solvable parameter) works
      env <- environment(el) %||% parent.frame()
      lhs <- eval(el[[2L]], env)
      rhs <- eval(el[[3L]], env)
      if (inherits(lhs, "mtk_array")) {
        # an array variable takes one value per element
        lv <- unclass(lhs); rv <- if (inherits(rhs, "mtk_array")) unclass(rhs) else rhs
        for (j in seq_along(lv)) out <- smap_set(out, lv[[j]], .unwrap_value(rv[[if (length(rv) == 1L) 1L else j]]))
        next
      }
      if (is.logical(rhs) && length(rhs) == 1L && is.na(rhs)) {
        out <- smap_set(out, lhs, NA_real_)
        next
      }
      el <- Equation(lhs, rhs)
    }
    if (inherits(el, "Equation")) {
      out <- smap_set(out, el$lhs, .unwrap_value(el$rhs))
    } else if (is.list(el) && !is_sym(el) && length(el) == 2L && is_sym(el[[1L]])) {
      out <- smap_set(out, el[[1L]], el[[2L]])
    } else if (!is.null(nms) && nzchar(nms[i])) {
      out[[paste0("@name:", nms[i])]] <- list(var = NULL, name = nms[i], value = el)
    } else {
      stop("cannot interpret element ", i, " of a variable map", call. = FALSE)
    }
  }
  class(out) <- "mtk_symmap"
  out
}

# Numbers wrapped as sym become plain doubles where that is lossless.
.unwrap_value <- function(v) {
  if (is_sym(v) && .type(v) %in% c("Integer", "Rational", "Float")) return(as.numeric(v))
  v
}

#' @export
print.mtk_symmap <- function(x, ...) {
  if (!length(x)) {
    cat("<empty map>\n")
    return(invisible(x))
  }
  for (el in unclass(x)) {
    lhs <- if (is.null(el$var)) el$name else format(el$var)
    val <- el$value
    rhs <- if (is_sym(val)) format(val) else if (length(val) == 1L && is.na(val)) "missing" else format(val)
    cat(lhs, " => ", paste(rhs, collapse = " "), "\n", sep = "")
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# Expression utilities

.is_num <- function(e) is_sym(e) && .type(e) %in% c("Integer", "Rational", "Float")
.is_zero <- function(e) .is_num(e) && as.numeric(e) == 0
.S0 <- function() as_sym(0)
.S1 <- function() as_sym(1)

# A variable application x(t) or x(t, s): a Function node whose head is a
# registered variable (not a mathematical function).
.is_varcall <- function(e) {
  if (!is_sym(e) || .type(e) != "Function") return(FALSE)
  k <- .var_kind(e)
  !is.null(k) && !(k %in% c("operator", "iv"))
}

.var_kind <- function(e) {
  if (!is_sym(e)) return(NULL)
  ty <- .type(e)
  if (ty == "Derivative") return(.var_kind(.args(e)[[1L]]))
  if (!(ty %in% c("Symbol", "Function"))) return(NULL)
  k <- attr(e, "mtk_meta")$kind
  if (!is.null(k)) return(k)
  ent <- .mtk$vars[[.name(e)]]
  if (is.null(ent)) return(.undeclared_kind(e))
  # a name can be both a symbol and a function (t and t(u)); see .register_var
  (ent$forms[[ty]] %||% ent)$kind %||% .undeclared_kind(e)
}

# Symbols and function applications made with symcompR directly, without
# Variable() or Parameter(), are unknowns, as symbols without metadata are in
# ModelingToolkit (collect_var!). Internal placeholders (names starting with
# "." or "_"), registered symbolic functions and the mathematical functions
# of symcompR and base R are not variables.
.undeclared_kind <- function(e) {
  nm <- .name(e)
  if (is.null(nm) || !nzchar(nm) || startsWith(nm, ".") || startsWith(nm, "_")) return(NULL)
  if (.type(e) == "Function") {
    if (exists(nm, envir = .mtk$funs, inherits = FALSE)) return(NULL)
    if (nm %in% .math_function_names()) return(NULL)
  }
  "unknown"
}

# Heads of mathematical functions: those symcompR evaluates (the generated
# .MATH_HEADS, see tools/gen_math_heads.R) and those code generation maps to
# R functions. Any other undeclared function head is a variable, so a model
# may use S(t), I(t) or N(t) even though symcompR and base R export
# functions with those names.
.math_function_names <- function() {
  if (is.null(.mtk$math_names))
    .mtk$math_names <- unique(c(.MATH_HEADS, names(.R_MATH), "log", "Max", "Min", "sec", "csc",
                                "cot", "Heaviside", "erf", "erfc", "Mod", "DiracDelta", "sinc"))
  .mtk$math_names
}

# Any declared symbolic variable: unknown, parameter, constant, brownian or
# independent variable, as a leaf node (Symbol or applied Function).
.is_declared <- function(e) !is.null(.var_kind(e)) && .type(e) != "Derivative"

.is_derivative <- function(e) is_sym(e) && .type(e) == "Derivative"

# For Derivative(x(t), (t, n)) return list(var = x(t), iv = t, order = n).
# Mixed partials return the first spec only; callers that care check n_specs.
.deriv_info <- function(e) {
  a <- .args(e)
  spec <- .args(a[[2L]])
  list(var = a[[1L]], iv = spec[[1L]], order = as.integer(as.numeric(spec[[2L]])),
       n_specs = length(a) - 1L)
}

.make_derivative <- function(var, iv, order = 1L) {
  if (order <= 0L) return(var)
  Derivative(var, iv, as.integer(order))
}

# All declared variable leaves and derivative nodes in an expression, unique
# by key, in first-seen order. Derivative nodes are returned as nodes and their
# inner variable is not returned separately unless it appears elsewhere.
.collect_vars <- function(e, acc = new.env(parent = emptyenv()), order = NULL) {
  out <- list()
  seen <- acc
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    k <- .key(z)
    if (!is.null(seen[[k]])) return(invisible())
    ty <- .type(z)
    if (ty == "Derivative" && .is_declared(.args(z)[[1L]])) {
      assign(k, TRUE, envir = seen)
      out[[length(out) + 1L]] <<- z
      return(invisible())
    }
    if (ty %in% c("Symbol", "Function") && .is_declared(z) &&
        !identical(.var_kind(z), "operator")) {
      assign(k, TRUE, envir = seen)
      out[[length(out) + 1L]] <<- z
      if (ty == "Function") for (a in .args(z)) walk(a)
      return(invisible())
    }
    assign(k, TRUE, envir = seen)
    for (a in .args(z)) walk(a)
  }
  if (is.list(e) && !is_sym(e)) for (z in e) walk(z) else walk(e)
  out
}

.unique_syms <- function(xs) {
  if (!length(xs)) return(list())
  ks <- vapply(xs, .key, "")
  xs[!duplicated(ks)]
}

.keys <- function(xs) if (length(xs)) vapply(xs, .key, "") else character()

.contains_key <- function(xs, e) .key(e) %in% .keys(xs)

# Top-down replacement by key. Replaces the largest matching subtree first, so
# Derivative(x(t), t) can map to one symbol while x(t) maps to another.
# `rule` is a named list: names are keys, values are replacement sym objects.
.replace <- function(e, rule) {
  if (!length(rule)) return(e)
  env <- attr(rule, "env") %||% list2env(unclass(rule), parent = emptyenv())
  walk <- function(x) {
    hit <- env[[.subset2(x, "key")]]
    if (!is.null(hit)) return(hit)
    a <- .subset2(x, "args")
    if (!length(a)) return(x)
    args <- lapply(a, walk)
    if (identical(vapply(args, .key, ""), vapply(a, .key, ""))) return(x)
    .rebuild(x, args)
  }
  walk(e)
}

# xreplace needs the left-hand expressions, so the rule carries them in an
# attribute; build rules with .rule() rather than by hand.
.rule <- function(from, to) {
  out <- list()
  lhs <- list()
  for (i in seq_along(from)) {
    k <- .key(from[[i]])
    out[[k]] <- as_sym(to[[i]])
    lhs[[k]] <- from[[i]]
  }
  attr(out, "lhs") <- lhs
  attr(out, "env") <- list2env(out, parent = emptyenv())
  out
}

.rule_lhs <- function(rule, k) attr(rule, "lhs")[[k]]

.rule_merge <- function(a, b) {
  out <- a
  lhs <- attr(a, "lhs") %||% list()
  for (k in names(b)) {
    out[[k]] <- b[[k]]
    lhs[[k]] <- attr(b, "lhs")[[k]]
  }
  attr(out, "lhs") <- lhs
  attr(out, "env") <- list2env(unclass(out)[names(out)], parent = emptyenv())
  out
}

.subst <- function(e, rule) {
  if (!length(rule)) return(e)
  if (inherits(e, "Equation")) return(Equation(.subst(e$lhs, rule), .subst(e$rhs, rule)))
  if (is.list(e) && !is_sym(e)) return(lapply(e, .subst, rule = rule))
  .replace(e, rule)
}

# Fixed point substitution: apply until no more change (bounded).
.subst_fixpoint <- function(e, rule, maxiter = 100L) {
  for (i in seq_len(maxiter)) {
    e2 <- .subst(e, rule)
    if (identical(.key(e2), .key(e))) return(e2)
    e <- e2
  }
  e
}

.fmt <- function(e) {
  if (is.null(e)) return("NULL")
  if (is_sym(e)) return(format(e))
  if (inherits(e, "Equation")) return(paste(format(e$lhs), "~", format(e$rhs)))
  paste(format(e), collapse = " ")
}

.stopf <- function(...) stop(sprintf(...), call. = FALSE)
.warnf <- function(...) warning(sprintf(...), call. = FALSE)
