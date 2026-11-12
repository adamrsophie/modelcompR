# Ported from StateSelection.jl (https://github.com/JuliaComputing/StateSelection.jl),
# commit 29b3c90. Portions are licensed under the GNU AGPL, version 3; see
# inst/COPYRIGHTS and inst/STATESELECTION-AGPL-LICENSE.
#
# Clock inference, after
# StateSelection.jl/lib/ModelingToolkitTearing/src/clock_inference/
# clock_inference.jl (ClockVertex, ClockInference, InferEquationClosure,
# InferVariableClosure, infer_clocks!, split_system), clock.jl
# (InferredClock, get_time_domain, input_timedomain and output_timedomain of
# Differential and Shift) and the operator time domains of Sample, Hold and
# SampleTime in ModelingToolkit.jl/src/systems/clock_inference.jl.
#
# The inference graph is a hypergraph over variables (every unknown and every
# application of a time varying operator, Derivative, Shift, Sample, Hold and
# SampleTime, as the fullvars of a TearingState), equations, initialization
# equations, clocks, inferred clocks, operator arguments (Expression
# vertices) and the integer sequence of purely discrete systems. Each
# equation, and each argument of an operator, is a hyperedge; the clock
# partitions are its connected components, found here by union find. The
# result gives the clock of every equation and unknown: 0 for the continuous
# clock, otherwise the id of the clock in .mtk$clocks. split_system is done by
# .compile_hybrid() in clocks-hybrid.R, which compiles the continuous
# partition and turns each clocked one into a periodic callback.
#
# Time domains are written as strings: "CONT" (ContinuousClock), "CLK:<id>"
# (a concrete clock), "INTSEQ" (IntegerSequence), "ID:<i>" (InferredDiscrete
# with relative index i) and "INF" (Inferred, no information).

.CI_OPS <- c("Sample", "Hold", "Shift", "SampleTime")

.ci_is_timevarying <- function(e) .is_derivative(e) || (.type(e) == "Function" && .is_op(e) && .name(e) %in% .CI_OPS)

# A concrete clock domain from a clock id (0 is no clock).
.ci_clock_domain <- function(cid) {
  if (is.null(cid) || cid == 0L) return(NULL)
  clk <- .clock_of(cid)
  if (identical(clk$type, "integer")) return("INTSEQ")
  if (identical(clk$type, "continuous")) return("CONT")
  paste0("CLK:", cid)
}

# The time domain declared on a variable (VariableTimeDomain metadata).
.ci_meta_domain <- function(v) {
  td <- attr(v, "mtk_meta")$timedomain
  if (!inherits(td, "TimeDomain")) return(NULL)
  if (identical(td$type, "integer")) return("INTSEQ")
  if (identical(td$type, "continuous")) return("CONT")
  paste0("CLK:", .clock_id(td))
}

.ci_clock_arg <- function(e, pos) {
  a <- .args(e)
  if (length(a) >= pos) as.integer(as.numeric(a[[pos]])) else 0L
}

# get_time_domain: the output domain of an operator application, else the
# declared domain of a variable.
.ci_get_time_domain <- function(e) {
  if (.ci_is_timevarying(e)) return(.ci_output_timedomain(e))
  .ci_meta_domain(e)
}

# input_timedomain of an operator application, one domain per argument.
.ci_input_timedomain <- function(e) {
  if (.is_derivative(e)) return(list("CONT"))
  switch(.name(e),
    Sample = list("CONT"),
    Hold = list(.ci_meta_domain(.args(e)[[1L]]) %||% "ID:0"),
    Shift = list(.ci_clock_domain(.ci_clock_arg(e, 3L)) %||% .ci_meta_domain(.args(e)[[1L]]) %||% "ID:0"),
    SampleTime = list())
}

.ci_output_timedomain <- function(e) {
  if (.is_derivative(e)) return("CONT")
  switch(.name(e),
    Sample = .ci_clock_domain(.ci_clock_arg(e, 2L)) %||% "ID:0",
    Hold = "CONT",
    Shift = "ID:0",
    SampleTime = .ci_clock_domain(.ci_clock_arg(e, 1L)) %||% "ID:1")
}

# The arguments of an operator application that carry time domains: the
# clock ids of Shift, Sample and SampleTime are not arguments.
.ci_op_args <- function(e) {
  if (.is_derivative(e)) return(list(.args(e)[[1L]]))
  switch(.name(e), SampleTime = list(), list(.args(e)[[1L]]))
}

# Variables of an expression as search_variables! sees them: unknown leaves
# and applications of time varying operators, which are atomic.
.ci_search_vars <- function(e, unknown_keys) {
  out <- list()
  seen <- new.env(parent = emptyenv())
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    k <- .key(z)
    if (!is.null(seen[[k]])) return(invisible())
    assign(k, TRUE, envir = seen)
    if (.ci_is_timevarying(z) || (k %in% unknown_keys)) {
      out[[length(out) + 1L]] <<- z
      return(invisible())
    }
    if (.type(z) == "Function" && .is_op(z)) return(invisible())
    for (a in .args(z)) walk(a)
  }
  walk(e)
  out
}

.ci_vertex_of_domain <- function(d) {
  if (is.null(d) || identical(d, "INF")) return(NULL)
  if (startsWith(d, "ID:")) return(NULL)
  d
}

# infer_clocks!: returns the domain of every equation (0 continuous, else a
# clock id), of every initialization equation, and a function giving the
# domain of an unknown.
.ci_infer_clocks <- function(eqs, unknown_keys, init_eqs = list()) {
  parent <- new.env(parent = emptyenv())
  uf_find <- function(a) {
    p <- parent[[a]]
    if (is.null(p)) {
      assign(a, a, envir = parent)
      return(a)
    }
    while (!identical(p, a)) {
      a <- p
      p <- parent[[a]]
    }
    a
  }
  add_edge <- function(vs) {
    vs <- unique(vs)
    if (!length(vs)) return(invisible())
    r <- uf_find(vs[1L])
    for (v in vs[-1L]) {
      rv <- uf_find(v)
      if (!identical(rv, r)) assign(rv, r, envir = parent)
    }
    invisible()
  }
  must <- character()
  var_terms <- list()
  vvert <- function(z) {
    k <- paste0("v:", .key(z))
    if (is.null(var_terms[[k]])) var_terms[[k]] <<- z
    k
  }
  # InferVariableClosure: adds to the parent hyperedge, returns it
  infer_var <- function(var, hyper) {
    hyper <- c(hyper, vvert(var))
    d <- .ci_get_time_domain(var)
    hyper <- c(hyper, .ci_vertex_of_domain(d))
    if (!.ci_is_timevarying(var)) return(hyper)
    if (.is_op(var, "Shift")) hyper <- c(hyper, .ci_vertex_of_domain(.ci_input_timedomain(var)[[1L]]))
    tdoms <- .ci_input_timedomain(var)
    args <- .ci_op_args(var)
    relative <- list()
    for (i in seq_along(args)) {
      arg <- args[[i]]; dom <- tdoms[[i]]
      ah <- character()
      for (v in .ci_search_vars(arg, unknown_keys)) ah <- infer_var(v, ah)
      if (!.is_num(arg)) ah <- c(ah, paste0("e:", .key(arg)))
      if (startsWith(dom, "ID:")) {
        relative[[dom]] <- c(relative[[dom]], ah)
        must <<- c(must, ah)
      } else if (!identical(dom, "INF")) {
        ah <- c(ah, dom)
      }
      add_edge(ah)
    }
    out <- .ci_output_timedomain(var)
    if (startsWith(out, "ID:")) {
      if (!is.null(relative[[out]])) {
        hyper <- c(hyper, relative[[out]])
        relative[[out]] <- NULL
      }
      must <<- c(must, hyper)
    } else if (!identical(out, "INF")) {
      hyper <- c(hyper, out)
    }
    for (r in relative) add_edge(r)
    hyper
  }
  infer_eq <- function(node, eq) {
    hyper <- node
    for (side in list(eq$lhs, eq$rhs)) for (v in .ci_search_vars(side, unknown_keys)) hyper <- infer_var(v, hyper)
    add_edge(hyper)
  }
  for (i in seq_along(eqs)) infer_eq(paste0("eq:", i), eqs[[i]])
  for (i in seq_along(init_eqs)) infer_eq(paste0("ieq:", i), init_eqs[[i]])
  # all shifted variables have the same clock as the unshifted variant
  for (k in names(var_terms)) {
    z <- var_terms[[k]]
    if (!.is_op(z, "Shift")) next
    base <- .args(z)[[1L]]
    td <- .ci_vertex_of_domain(.ci_input_timedomain(z)[[1L]])
    add_edge(c(k, vvert(base), td))
  }
  nodes <- ls(parent, all.names = TRUE)
  roots <- vapply(nodes, uf_find, "")
  dom <- list()
  for (r in unique(roots)) {
    members <- nodes[roots == r]
    clks <- unique(grep("^(CLK:|CONT$)", members, value = TRUE))
    if (!length(clks)) {
      if ("INTSEQ" %in% members) {
        clks <- "CONT"
      } else if (any(members %in% must)) {
        stop("Found a clock partition that must be discrete (due to the presence of an ",
             "`InferredDiscrete`) but does not have any associated clock (and would otherwise ",
             "then default to being on the continuous clock). This likely means that the ",
             "partition was not assigned a valid discrete clock and the model is incorrect.",
             call. = FALSE)
      } else {
        clks <- "CONT"
      }
    }
    if (length(clks) > 1L) {
      vs <- vapply(var_terms[intersect(names(var_terms), members)], format, "")
      stop("Found clock partition with multiple associated clocks. Involved variables: ",
           paste(vs, collapse = ", "), ". Involved clocks: ",
           paste(vapply(clks, function(c) if (c == "CONT") "ContinuousClock()" else
             format(.clock_of(as.integer(sub("CLK:", "", c)))), ""), collapse = ", "), ".", call. = FALSE)
    }
    if (clks == "CONT" && any(members %in% must))
      stop("Found a clock partition that must be discrete (due to the presence of an ",
           "`InferredDiscrete`) but is associated with a continuous clock. This is likely ",
           "a modeling error.", call. = FALSE)
    d <- if (clks == "CONT") 0L else as.integer(sub("CLK:", "", clks))
    for (m in members) dom[[m]] <- d
  }
  list(eq = vapply(seq_along(eqs), function(i) dom[[paste0("eq:", i)]] %||% 0L, 0L),
       init_eq = vapply(seq_along(init_eqs), function(i) dom[[paste0("ieq:", i)]] %||% 0L, 0L),
       var = function(v) dom[[paste0("v:", .key(v))]] %||% 0L)
}
