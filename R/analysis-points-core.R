# Analysis points and their transformations, after
# ModelingToolkitBase/src/systems/analysis_points.jl (AnalysisPoint,
# remove_analysis_points, ap_var, modify_nested_subsystem, Break, GetInput,
# PerturbOutput, AddVariable, ComplementarySensitivityTransform,
# LoopTransferTransform, canonicalize_ap, open_loop,
# generate_control_function) and ModelingToolkit/src/systems/analysis_points.jl
# (handle_loop_openings, get_linear_analysis_function, get_sensitivity,
# get_comp_sensitivity, get_looptransfer, linearization_ap_transform).
#
# An analysis point is an object of class c("AnalysisPoint", "Connection")
# placed directly in an equation list; `sys$name` returns it namespaced.

#' Analysis points
#'
#' An analysis point marks a causal connection between an output and one or
#' more inputs for linear analysis. It is usually created by
#' `connect(out, "name", in)` and referred to as `sys$name` on the
#' unsimplified model; `AnalysisPoint("name")` creates an empty point with a
#' name only. `ap_var()` returns the variable of one end of a point: the
#' variable `u` of a connector, its only unknown, or a causal variable
#' itself. `expand_connections()` turns analysis points into ordinary
#' connections.
#'
#' @param input The output connector (or variable) feeding the point, or a
#'   name when it is the only argument.
#' @param name The name of the point, a string.
#' @param outputs A list of input connectors (or variables) it feeds.
#' @param verbose Warn when the causality looks reversed.
#' @param x A connector system or a variable.
#' @return An `AnalysisPoint`; `ap_var()` returns a variable.
#' @examples
#' t <- t_nounits
#' u_out <- Variable("y", t, output = TRUE)
#' u_in <- Variable("u", t, input = TRUE)
#' ap <- connect(u_out, "plant_input", u_in)
#' ap
#' ap_var(ap$input)
#' @export
AnalysisPoint <- function(input = NULL, name = NULL, outputs = NULL, verbose = TRUE) {
  if (is.character(input) && is.null(name)) {
    name <- input
    input <- NULL
  }
  if (is.null(name)) name <- ""
  if (!is.null(outputs) && (inherits(outputs, "System") || is_sym(outputs))) outputs <- list(outputs)
  ap <- structure(list(input = input, name = as.character(name), outputs = outputs),
                  class = c("AnalysisPoint", "Connection"))
  if (isTRUE(verbose)) {
    if (!is.null(input) && !isoutput(ap_var(input))) .ap_warning(1L, name, TRUE)
    for (i in seq_along(outputs)) if (!isinput(ap_var(outputs[[i]]))) .ap_warning(2L + i, name, FALSE)
  }
  ap
}

.ap_warning <- function(arg, name, should_be_output) {
  .warnf(paste("The %d-th argument to analysis point %s was not a %s. This is supported in order",
               "to handle inverse models, but may not be what you intended. Silence this message",
               "with connect(out, name, in, verbose = FALSE)."),
         arg, name, if (should_be_output) "output" else "input")
}

#' @export
format.AnalysisPoint <- function(x, ...) {
  if (is.null(x$input)) return(paste0("AnalysisPoint(", x$name, ")"))
  outs <- vapply(x$outputs, function(o) format(ap_var(o)), "")
  paste0("AnalysisPoint(", x$name, " from ", format(ap_var(x$input)), " to ",
         if (length(outs) == 1L) outs else paste0("[", paste(outs, collapse = ", "), "]"), ")")
}

#' @export
print.AnalysisPoint <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' @rdname AnalysisPoint
#' @export
ap_var <- function(x) {
  if (inherits(x, "System")) {
    u <- tryCatch(.getproperty(x, "u"), error = function(e) NULL)
    if (is_sym(u)) return(u)
    us <- unknowns(x)
    if (length(us) == 1L) return(renamespace(x, us[[1L]]))
    stop("Could not determine the analysis-point variable in system ", nameof(x),
         ". To use an analysis point, apply it to a connection between causal blocks which ",
         "have a variable named `u` or a single unknown of the same size.", call. = FALSE)
  }
  as_sym(x)
}

.to_connection <- function(ap) {
  if (is.null(ap$input)) stop("Empty analysis point ", ap$name, call. = FALSE)
  Connection(c(list(ap$input), ap$outputs))
}

#' Remove analysis points
#'
#' Replaces every analysis point of a system and its subsystems by the
#' equivalent connection.
#'
#' @param sys A `System`.
#' @return The `System`.
#' @examples
#' t <- t_nounits
#' y <- Variable("y", t, output = TRUE); u <- Variable("u", t, input = TRUE)
#' src <- System(list(y ~ 1), t, name = "src")
#' snk <- System(list(), t, list(u), list(), name = "snk")
#' sys <- System(list(connect(src$y, "ap", snk$u)), t, systems = list(src, snk), name = "s")
#' get_eqs(remove_analysis_points(sys))
#' @keywords internal
#' @export
remove_analysis_points <- function(sys) {
  eqs <- lapply(.f(sys, "eqs"), function(e) if (inherits(e, "AnalysisPoint")) .to_connection(e) else e)
  .update_system(sys, eqs = eqs, systems = lapply(.f(sys, "systems"), remove_analysis_points))
}

# modify_nested_subsystem: apply fn to the subsystem of root at `hierarchy`
# (which starts with the name of root); fn returns list(sys, vars) and the
# variables are namespaced up to, but excluding, root.
.modify_nested <- function(fn, root, hierarchy) {
  if (!length(hierarchy)) return(fn(root))
  if (nameof(root) != hierarchy[1L])
    stop("Invalid analysis point name `", paste(hierarchy, collapse = NAMESPACE_SEPARATOR),
         "`. The name must include the name of the root system `", nameof(root), "`. This ",
         "typically happens when using an analysis point obtained from a complete system to ",
         "linearize a system that is not complete.", call. = FALSE)
  hierarchy <- hierarchy[-1L]
  helper <- function(sys, i) {
    if (i > length(hierarchy)) {
      res <- fn(sys)
    } else {
      subs <- .f(sys, "systems")
      idx <- match(hierarchy[i], vapply(subs, nameof, ""))
      if (is.na(idx))
        stop("System ", paste(c(nameof(root), hierarchy[seq_len(i - 1L)]), collapse = "."),
             " does not have a subsystem named ", hierarchy[i], ".", call. = FALSE)
      res <- helper(subs[[idx]], i + 1L)
      subs[[idx]] <- res$sys
      res$sys <- .update_system(sys, systems = subs)
    }
    if (i != 1L) res$vars <- lapply(res$vars, function(v) renamespace(res$sys, v))
    res
  }
  helper(root, 1L)
}

.ap_hierarchy <- function(ap) utils::head(.namespace_hierarchy(ap$name), -1L)

.ap_index <- function(sys, ap) {
  nm <- utils::tail(.namespace_hierarchy(ap$name), 1L)
  for (i in seq_along(.f(sys, "eqs"))) {
    e <- .f(sys, "eqs")[[i]]
    if (inherits(e, "AnalysisPoint") && identical(e$name, nm)) return(i)
  }
  stop("Analysis point ", ap$name, " not found in system ", nameof(sys), ".", call. = FALSE)
}

# get_analysis_variable: a new variable name(t), or d_name(t) for a
# perturbation, marked as an analysis variable; its default is zero.
.analysis_variable <- function(name, iv, perturb = TRUE) {
  if (perturb) name <- paste0("d_", name)
  .make_var(name, "unknown", iv, list(analysis_variable = TRUE))
}

.with_ap_ignored <- function(sys, ap) {
  .update_system(sys, ignored_connections = c(.f(sys, "ignored_connections"), list(.to_connection(ap))))
}

.unknown_to_parameter <- function(sys, var) {
  us <- .f(sys, "unknowns")
  k <- .key(var)
  if (!(k %in% .keys(us))) stop("Variable ", format(var), " is not an unknown of the system.", call. = FALSE)
  .update_system(sys, unknowns = us[.keys(us) != k], ps = c(.f(sys, "ps"), list(var)))
}

# Primitive transformations. Each returns list(sys, vars).

.ap_break <- function(ap, sys, outputs_to_params = FALSE) {
  .modify_nested(function(bsys) {
    idx <- .ap_index(bsys, ap)
    eqs <- .f(bsys, "eqs")
    apo <- eqs[[idx]]
    bsys <- .with_ap_ignored(.update_system(bsys, eqs = eqs[-idx]), apo)
    out_vars <- list()
    for (out in apo$outputs) {
      if (inherits(out, "System")) {
        if (outputs_to_params) {
          res <- .modify_nested(function(osys) {
            ov <- ap_var(toggle_namespacing(osys, FALSE))
            list(sys = .unknown_to_parameter(osys, ov), vars = list(ov))
          }, bsys, .namespace_hierarchy(nameof(renamespace(bsys, out))))
          bsys <- res$sys
          ov <- res$vars[[1L]]
        } else ov <- ap_var(out)
      } else {
        ov <- ap_var(out)
        if (outputs_to_params) {
          path <- c(nameof(bsys), .namespace_hierarchy(.name(ov)))
          bsys <- .modify_nested(function(asys) {
            v <- .getproperty(toggle_namespacing(asys, FALSE), utils::tail(path, 1L))
            list(sys = .unknown_to_parameter(asys, v), vars = list())
          }, bsys, utils::head(path, -1L))$sys
        }
      }
      out_vars[[length(out_vars) + 1L]] <- ov
    }
    list(sys = bsys, vars = list(out_vars))
  }, sys, .ap_hierarchy(ap))
}

.ap_get_input <- function(ap, sys) {
  .modify_nested(function(asys) {
    apo <- .f(asys, "eqs")[[.ap_index(asys, ap)]]
    list(sys = asys, vars = list(ap_var(apo$input)))
  }, sys, .ap_hierarchy(ap))
}

.ap_perturb_output <- function(ap, sys, with_output = FALSE) {
  iv <- .f(sys, "iv")
  .modify_nested(function(asys) {
    idx <- .ap_index(asys, ap)
    eqs <- .f(asys, "eqs")
    apo <- eqs[[idx]]
    eqs <- eqs[-idx]
    ivar <- ap_var(apo$input)
    new_var <- .analysis_variable(apo$name, iv)
    for (o in apo$outputs) eqs[[length(eqs) + 1L]] <- Equation(ap_var(o), ivar + new_var)
    unks <- c(.f(asys, "unknowns"), list(new_var))
    ics <- smap_set(.f(asys, "initial_conditions") %||% symmap(), new_var, 0)
    vars <- list(new_var)
    if (with_output) {
      out_var <- .analysis_variable(apo$name, iv, perturb = FALSE)
      eqs[[length(eqs) + 1L]] <- Equation(out_var, ivar + new_var)
      unks <- c(unks, list(out_var))
      vars <- c(vars, list(out_var))
    }
    asys <- .update_system(asys, eqs = eqs, unknowns = unks, initial_conditions = ics)
    list(sys = .with_ap_ignored(asys, apo), vars = vars)
  }, sys, .ap_hierarchy(ap))
}

.ap_add_variable <- function(ap, sys, name = NULL) {
  iv <- .f(sys, "iv")
  .modify_nested(function(asys) {
    apo <- .f(asys, "eqs")[[.ap_index(asys, ap)]]
    v <- .analysis_variable(name %||% apo$name, iv, perturb = FALSE)
    list(sys = .update_system(asys, unknowns = c(.f(asys, "unknowns"), list(v))), vars = list(v))
  }, sys, .ap_hierarchy(ap))
}

# Derived transformations; each returns list(sys, vars = list(du, u)).
.sensitivity_transform <- function(ap, sys) .ap_perturb_output(ap, sys, TRUE)

.comp_sensitivity_transform <- function(ap, sys) {
  r <- .ap_get_input(ap, sys)
  u <- r$vars[[1L]]
  nm <- paste0(utils::tail(.namespace_hierarchy(ap$name), 1L), "_comp_sens_du")
  r <- .ap_add_variable(ap, r$sys, nm)
  du <- r$vars[[1L]]
  r <- .ap_perturb_output(ap, r$sys)
  ddu <- r$vars[[1L]]
  # PerturbOutput adds input + _du ~ output, the complementary sensitivity
  # wants output + du ~ input, so du ~ -_du.
  sys <- .update_system(r$sys, eqs = c(.f(r$sys, "eqs"), list(Equation(du, -ddu))),
                        initial_conditions = smap_set(.f(r$sys, "initial_conditions") %||% symmap(), du, -ddu))
  list(sys = sys, vars = list(du, u))
}

.looptransfer_transform <- function(ap, sys) {
  r <- .ap_get_input(ap, sys)
  u <- r$vars[[1L]]
  r <- .ap_break(ap, r$sys)
  list(sys = r$sys, vars = list(r$vars[[1L]], u))
}

# canonicalize_ap: a flat list of analysis points namespaced by sys; names
# (strings) become analysis points, anything else is kept.
.canonicalize_ap <- function(sys, ap) {
  if (is.null(ap)) return(list())
  if (inherits(ap, "AnalysisPoint")) {
    return(list(if (isTRUE(.f(sys, "namespacing"))) ap else renamespace(sys, ap)))
  }
  if (is.character(ap)) return(lapply(ap, function(a) AnalysisPoint(name = .ns_name(nameof(sys), a), verbose = FALSE)))
  if (is_sym(ap)) return(list(ap))
  if (is.list(ap)) return(do.call(c, c(list(list()), lapply(ap, .canonicalize_ap, sys = sys))))
  list(ap)
}

.is_ap_like <- function(x) {
  if (inherits(x, "AnalysisPoint") || is.character(x)) return(TRUE)
  is.list(x) && !is_sym(x) && length(x) > 0L && all(vapply(x, .is_ap_like, TRUE))
}

#' Open a loop at an analysis point
#'
#' `open_loop(sys, ap)` breaks the connection of the analysis point `ap` of
#' the unsimplified model `sys`, as for a loop transfer function. The
#' signal that fed the inputs of the point becomes the new variable `du`
#' (the inputs are now driven by it), and `u` is the variable that drove the
#' point.
#'
#' @param sys A `System` that has not been completed or compiled.
#' @param ap An analysis point `sys$name`, or its name.
#' @param system_modifier A function applied to the modified system.
#' @return A list with the modified system `sys`, the list `du` of new input
#'   variables and the driving variable `u`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); u <- Variable("u", t, input = TRUE)
#' y <- Variable("y", t, output = TRUE); e <- Variable("e", t, output = TRUE)
#' sys <- System(list(D(x) ~ -x + u, y ~ x, e ~ -y, connect(e, "loop", u)),
#'   t, name = "fb")
#' ol <- open_loop(sys, sys$loop)
#' ol$du; ol$u
#' @export
open_loop <- function(sys, ap, system_modifier = identity) {
  ap <- .canonicalize_ap(sys, ap)
  if (length(ap) != 1L) stop("open_loop() takes one analysis point", call. = FALSE)
  r <- .looptransfer_transform(ap[[1L]], sys)
  list(sys = system_modifier(r$sys), du = r$vars[[1L]], u = r$vars[[2L]])
}

# handle_loop_openings: break each point, turning its outputs into solvable
# parameters (bindings NA) held at the value of their driver.
.handle_loop_openings <- function(sys, aps) {
  params <- list(); drivers <- list()
  for (ap in .canonicalize_ap(sys, aps)) {
    r <- .ap_get_input(ap, sys)
    driver <- r$vars[[1L]]
    r <- .ap_break(ap, r$sys, TRUE)
    sys <- r$sys
    dv <- r$vars[[1L]]
    params <- c(params, dv)
    drivers <- c(drivers, rep(list(driver), length(dv)))
  }
  if (length(params)) {
    b <- .f(sys, "bindings") %||% symmap()
    for (p in params) b <- smap_set(b, p, NA)
    sys <- .update_system(sys, bindings = b)
  }
  list(sys = sys, params = params, drivers = drivers)
}

.get_linear_analysis_function <- function(sys, .tf, aps, system_modifier = identity,
                                          loop_openings = list(), ...) {
  lo <- .handle_loop_openings(sys, loop_openings)
  sys <- lo$sys
  dus <- list(); us <- list()
  for (ap in .canonicalize_ap(sys, aps)) {
    r <- .tf(ap, sys)
    sys <- r$sys
    du <- .flatten_vars(r$vars[[1L]]); u <- .flatten_vars(r$vars[[2L]])
    dus <- c(dus, du); us <- c(us, u)
    if (identical(.tf, .looptransfer_transform)) {
      for (d in du) {
        lo$params <- c(lo$params, list(d))
        lo$drivers <- c(lo$drivers, u[1L])
      }
    }
  }
  linearization_function(system_modifier(sys), dus, us, loop_opening_params = lo$params,
                         loop_opening_drivers = lo$drivers, ...)
}

#' Sensitivity functions at analysis points
#'
#' `get_sensitivity(sys, ap)` linearizes the sensitivity function at the
#' analysis point `ap`: the transfer from a perturbation added to the input
#' of the point to the signal at its outputs. `get_comp_sensitivity()` gives
#' the complementary sensitivity function and `get_looptransfer()` the loop
#' transfer function, which includes the negative feedback of the loop.
#' The `*_function` variants return the linearization function and the
#' simplified system instead, see [linearization_function()].
#'
#' @param sys The unsimplified, uncompleted model.
#' @param ap An analysis point `sys$name` (or `sys$inner$name`), its name,
#'   or a list of them.
#' @param loop_openings Analysis points whose connections are removed
#'   before linearization; the opened signals keep the value of their
#'   drivers unless `op` gives one.
#' @param system_modifier A function applied to the transformed system.
#' @param allow_input_derivatives,op,t,eliminate_algebraic Passed to
#'   [linearize()].
#' @param ... Passed to [linearization_function()].
#' @return For `get_*`, the result of [linearize()]: a list with matrices
#'   `A`, `B`, `C`, `D`, the simplified system `sys` and `extras`. For
#'   `get_*_function`, a list with `lin_fun` and `sys`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); u <- Variable("u", t, input = TRUE)
#' y <- Variable("y", t, output = TRUE); e <- Variable("e", t, output = TRUE)
#' sys <- System(list(D(x) ~ -x + u, y ~ x, e ~ -y, connect(e, "loop", u)),
#'   t, name = "fb")
#' S <- get_sensitivity(sys, sys$loop, verbose = FALSE)
#' S$A; S$D
#' @export
get_sensitivity <- function(sys, ap, loop_openings = list(), system_modifier = identity,
                            allow_input_derivatives = TRUE, op = list(), t = 0,
                            eliminate_algebraic = FALSE, ...) {
  lf <- get_sensitivity_function(sys, ap, loop_openings = loop_openings,
                                 system_modifier = system_modifier, op = op, t = t, ...)
  linearize(lf$sys, lf$lin_fun, op = op, t = t, allow_input_derivatives = allow_input_derivatives,
            eliminate_algebraic = eliminate_algebraic)
}

#' @rdname get_sensitivity
#' @export
get_comp_sensitivity <- function(sys, ap, loop_openings = list(), system_modifier = identity,
                                 allow_input_derivatives = TRUE, op = list(), t = 0,
                                 eliminate_algebraic = FALSE, ...) {
  lf <- get_comp_sensitivity_function(sys, ap, loop_openings = loop_openings,
                                      system_modifier = system_modifier, op = op, t = t, ...)
  linearize(lf$sys, lf$lin_fun, op = op, t = t, allow_input_derivatives = allow_input_derivatives,
            eliminate_algebraic = eliminate_algebraic)
}

#' @rdname get_sensitivity
#' @export
get_looptransfer <- function(sys, ap, loop_openings = list(), system_modifier = identity,
                             allow_input_derivatives = TRUE, op = list(), t = 0,
                             eliminate_algebraic = FALSE, ...) {
  lf <- get_looptransfer_function(sys, ap, loop_openings = loop_openings,
                                  system_modifier = system_modifier, op = op, t = t, ...)
  linearize(lf$sys, lf$lin_fun, op = op, t = t, allow_input_derivatives = allow_input_derivatives,
            eliminate_algebraic = eliminate_algebraic)
}

#' @rdname get_sensitivity
#' @export
get_sensitivity_function <- function(sys, ap, ...) .get_linear_analysis_function(sys, .sensitivity_transform, ap, ...)

#' @rdname get_sensitivity
#' @export
get_comp_sensitivity_function <- function(sys, ap, ...) .get_linear_analysis_function(sys, .comp_sensitivity_transform, ap, ...)

#' @rdname get_sensitivity
#' @export
get_looptransfer_function <- function(sys, ap, ...) .get_linear_analysis_function(sys, .looptransfer_transform, ap, ...)

#' Prepare a model with analysis points for linearization
#'
#' `linearization_ap_transform(sys, inputs, outputs, loop_openings)` applies
#' the analysis point transformations used by [linearize()] when its inputs
#' or outputs are analysis points: input points get a perturbation variable
#' added to their signal (or are broken when also listed in
#' `loop_openings`), output points get a variable equal to their signal.
#'
#' @param sys The unsimplified, uncompleted model.
#' @param inputs,outputs Analysis points (or names); outputs may also be
#'   variables.
#' @param loop_openings Analysis points to break.
#' @return A list with the transformed `sys`, `inputs`, `outputs`,
#'   `loop_opening_params` and `loop_opening_drivers`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' x <- Variable("x", t); u <- Variable("u", t, input = TRUE)
#' y <- Variable("y", t, output = TRUE); e <- Variable("e", t, output = TRUE)
#' sys <- System(list(D(x) ~ -x + u, y ~ x, e ~ -y, connect(e, "loop", u)),
#'   t, name = "fb")
#' tr <- linearization_ap_transform(sys, sys$loop, list(y), list())
#' tr$inputs
#' @export
linearization_ap_transform <- function(sys, inputs, outputs, loop_openings = list()) {
  lo_names <- vapply(.canonicalize_ap(sys, loop_openings), function(a) a$name, "")
  inputs <- .canonicalize_ap(sys, inputs)
  outputs <- .canonicalize_ap(sys, outputs)
  in_names <- vapply(inputs, function(a) a$name, "")
  for (o in outputs) {
    if (inherits(o, "AnalysisPoint") && o$name %in% in_names)
      stop("The analysis point ", o$name, " is used both as an input and as an output of the ",
           "linearization. Use get_sensitivity() or get_comp_sensitivity() instead.", call. = FALSE)
  }
  input_vars <- list(); opened <- list(); opened_drivers <- list()
  for (inp in inputs) {
    if (inp$name %in% lo_names) {
      lo_names <- setdiff(lo_names, inp$name)
      r <- .ap_get_input(inp, sys)
      driver <- r$vars[[1L]]
      r <- .ap_break(inp, r$sys)
      sys <- r$sys
      iv <- r$vars[[1L]]
      input_vars <- c(input_vars, iv)
      opened <- c(opened, iv)
      opened_drivers <- c(opened_drivers, rep(list(driver), length(iv)))
    } else {
      r <- .ap_perturb_output(inp, sys)
      sys <- r$sys
      input_vars <- c(input_vars, r$vars[1L])
    }
  }
  output_vars <- list()
  for (o in outputs) {
    if (inherits(o, "AnalysisPoint")) {
      r <- .ap_add_variable(o, sys)
      ov <- r$vars[[1L]]
      r <- .ap_get_input(o, r$sys)
      sys <- .update_system(r$sys, eqs = c(.f(r$sys, "eqs"), list(Equation(ov, r$vars[[1L]]))))
    } else ov <- o
    output_vars <- c(output_vars, list(ov))
  }
  lo <- .handle_loop_openings(sys, lapply(lo_names, function(n) AnalysisPoint(name = n, verbose = FALSE)))
  list(sys = lo$sys, inputs = input_vars, outputs = output_vars,
       loop_opening_params = c(lo$params, opened), loop_opening_drivers = c(lo$drivers, opened_drivers))
}

