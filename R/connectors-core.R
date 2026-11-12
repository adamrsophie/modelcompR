# Connectors and connection expansion, after
# ModelingToolkitBase/src/systems/connectors.jl (connect, domain_connect,
# connector_type, instream, generate_connection_set, expand_connections,
# expand_instream, get_domain_bindings, warn_unbalanced_connectors) and
# connectiongraph.jl (ConnectionVertex, the connection hypergraph,
# remove_negative_connections!, connectionsets).
#
# A connection is an S3 object of class "Connection" placed directly in the
# equation list of a system; analysis points (R/analysis-points-core.R)
# inherit from it. `expand_connections()` replaces them by equations before
# flattening.

#' Connector variable types
#'
#' `Flow` and `Stream` are the values of the `connect` metadata of connector
#' variables: `Variable("i", t, connect = Flow)`. In a connection, flow
#' variables sum to zero (outside connectors with a negative sign) and all
#' other variables, of type `"Equality"` (the default), are set equal.
#' Stream variables carry quantities transported by a flow and are read with
#' [instream()]. ModelingToolkit's `Equality` is the string `"Equality"`
#' here, since the name is taken by symcompR.
#'
#' @format Character strings.
#' @examples
#' i <- Variable("i", t_nounits, connect = Flow)
#' getconnect(i)
#' @export
Flow <- "Flow"

#' @rdname Flow
#' @export
Stream <- "Stream"

.connect_type <- function(v) .meta_get(v, "connect") %||% "Equality"

#' Connection objects
#'
#' `Connection(systems)` is the object returned by [connect()]: a list of
#' connected connector systems or causal variables. It is placed directly
#' in the equation list of a system.
#'
#' @param systems A list of connector systems or variables.
#' @param domain Whether this is a domain only connection, see
#'   [domain_connect()].
#' @return An object of class `Connection`.
#' @examples
#' t <- t_nounits
#' pin <- function(name) connector(System(list(), t, list(Variable("v", t),
#'   Variable("i", t, connect = Flow)), list(), name = name))
#' Connection(list(pin("a"), pin("b")))
#' @export
Connection <- function(systems = list(), domain = FALSE) {
  structure(list(systems = systems, domain = isTRUE(domain)), class = "Connection")
}

#' @export
format.Connection <- function(x, ...) {
  nm <- vapply(x$systems, function(s) {
    n <- if (inherits(s, "System")) nameof(s) else format(s)
    gsub(NAMESPACE_SEPARATOR, ".", n, fixed = TRUE)
  }, "")
  paste0(if (isTRUE(x$domain)) "domain_connect(" else "connect(", paste(nm, collapse = ", "), ")")
}

#' @export
print.Connection <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' Connect connectors
#'
#' `connect(a, b, ...)` connects two or more connector systems (built with
#' [connector()]) or causal variables, and is placed directly in the
#' equation list of a system. Connected potential (`"Equality"`) variables
#' are set equal and [Flow] variables sum to zero. Causal variables must be
#' marked `input` or `output`; the first is normally the output and the
#' rest are inputs.
#'
#' `connect(out, "name", in, ...)` instead creates an [AnalysisPoint] called
#' `name` between the output `out` and the inputs, for linear analysis.
#' `domain_connect()` connects domains only: it shares the parameters of a
#' domain source without generating equations.
#'
#' @param ... Connectors or causal variables; for an analysis point, the
#'   output, the name (a string or an [AnalysisPoint]) and the inputs.
#' @param verbose For analysis points, warn when the causality looks
#'   reversed.
#' @return A `Connection` (or `AnalysisPoint`).
#' @examples
#' t <- t_nounits
#' Pin <- connector(function(name) System(list(), t,
#'   list(Variable("v", t), Variable("i", t, connect = Flow)), list(), name = name))
#' a <- Pin(name = "a"); b <- Pin(name = "b")
#' connect(a, b)
#' @export
connect <- function(..., verbose = TRUE) {
  args <- list(...)
  if (length(args) < 2L) stop("connect() needs at least two arguments", call. = FALSE)
  second <- args[[2L]]
  if (is.character(second) || inherits(second, "AnalysisPoint")) {
    name <- if (is.character(second)) second else second$name
    if (length(args) < 3L) stop("an analysis point needs an input and at least one output", call. = FALSE)
    outs <- args[-c(1L, 2L)]
    if (!inherits(args[[1L]], "System")) {
      allvars <- lapply(c(args[1L], outs), as_sym)
      .validate_causal_connection(allvars)
      return(AnalysisPoint(allvars[[1L]], name, allvars[-1L], verbose = verbose))
    }
    return(AnalysisPoint(args[[1L]], name, outs, verbose = verbose))
  }
  if (all(vapply(args, inherits, TRUE, "System"))) {
    nms <- vapply(args, nameof, "")
    if (anyDuplicated(nms)) stop("connect takes distinct systems!", call. = FALSE)
    return(Connection(args))
  }
  allvars <- lapply(args, as_sym)
  .validate_causal_connection(allvars)
  Connection(allvars)
}

#' @rdname connect
#' @export
domain_connect <- function(...) {
  args <- list(...)
  if (length(args) < 2L || !all(vapply(args, inherits, TRUE, "System")))
    stop("domain_connect() needs at least two connector systems", call. = FALSE)
  if (anyDuplicated(vapply(args, nameof, ""))) stop("connect takes distinct systems!", call. = FALSE)
  Connection(args, domain = TRUE)
}

.validate_causal_connection <- function(allvars) {
  for (v in allvars) {
    if (!identical(.var_kind(v), "unknown"))
      stop("expected ", format(v), " to be a variable", call. = FALSE)
  }
  if (anyDuplicated(.keys(allvars)))
    stop("expected all connection variables to be unique", call. = FALSE)
  bad <- Filter(function(v) !(isinput(v) || isoutput(v)), allvars)
  if (length(bad))
    stop("Only causal variables can be used in a `connect` statement. Each variable must be ",
         "either an input or an output (metadata input = TRUE or output = TRUE). ",
         "The following variables were found to be non-causal: ",
         paste(vapply(bad, format, ""), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

#' Stream connector inflow
#'
#' `instream(x)` is the value of the stream variable `x` flowing into the
#' connector from the connection, used in models with [Stream] variables.
#' It is replaced when connections are expanded.
#'
#' @param x A stream variable.
#' @return A `sym` expression.
#' @examples
#' h <- Variable("h", t_nounits, connect = Stream)
#' instream(h)
#' @export
instream <- function(x) .op_apply("instream", as_sym(x))

# ---------------------------------------------------------------------------
# Connector systems

#' Connector type of a system
#'
#' `connector_type(sys)` classifies a connector: `"DomainConnector"` when it
#' has a single unknown which is a flow variable, `"StreamConnector"` when it
#' has stream variables, and `"RegularConnector"` otherwise.
#' `isconnector(sys)` tests whether a system was marked as a connector.
#'
#' @param sys A `System`.
#' @return A string, or a logical.
#' @examples
#' t <- t_nounits
#' p <- System(list(), t, list(Variable("v", t), Variable("i", t, connect = Flow)),
#'   list(), name = "p")
#' connector_type(p)
#' @export
connector_type <- function(sys) {
  cnt <- .connector_variable_counts(sys)
  if (cnt[["stream"]] > 0L && cnt[["flow"]] > 1L)
    stop("There are multiple flow variables in the stream connector ", nameof(sys), "!", call. = FALSE)
  if (cnt[["flow"]] == 1L && length(.f(sys, "unknowns")) == 1L) return("DomainConnector")
  if (cnt[["stream"]] > 0L) "StreamConnector" else "RegularConnector"
}

#' @rdname connector_type
#' @export
isconnector <- function(sys) inherits(sys, "System") && !is.null(.f(sys, "connector_type"))

.is_domain_connector <- function(s) isconnector(s) && identical(.f(s, "connector_type"), "DomainConnector")

.connector_variable_counts <- function(sys) {
  n_stream <- 0L; n_flow <- 0L; n_regular <- 0L
  for (s in .f(sys, "unknowns")) {
    ty <- .connect_type(s)
    if (ty == "Stream") n_stream <- n_stream + 1L
    else if (ty == "Flow") n_flow <- n_flow + 1L
    else if (!(isinput(s) || isoutput(s))) n_regular <- n_regular + 1L
  }
  c(stream = n_stream, flow = n_flow, regular = n_regular)
}

# Warn once per connector type with different numbers of flow and regular
# variables (warn_unbalanced_connectors).
.warn_unbalanced_connectors <- function(sys) {
  seen <- character()
  stack <- list(sys)
  while (length(stack)) {
    s <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL
    stack <- c(stack, .f(s, "systems"))
    if (!isconnector(s) || .is_domain_connector(s)) next
    if (isTRUE(.f(s, "metadata")$is_frame)) next
    key <- .f(s, "metadata")$gui_metadata %||% nameof(s)
    if (key %in% seen) next
    cnt <- .connector_variable_counts(s)
    if (cnt[["flow"]] == cnt[["regular"]]) next
    seen <- c(seen, key)
    .warnf(paste("Connector %s contains %d flow variables, yet %d regular (non-flow, non-stream,",
                 "non-input, non-output) variables. This could lead to imbalanced model that are",
                 "difficult to debug. Consider marking some of the regular variables as",
                 "input/output variables."), key, cnt[["flow"]], cnt[["regular"]])
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# The connection hypergraph (connectiongraph.jl)

# Names separated by "." or the namespace separator.
.namespace_hierarchy <- function(name) strsplit(name, paste0("[.", NAMESPACE_SEPARATOR, "]"))[[1L]]

# A vertex: the path of names from the root (excluding the root) to a
# variable or connector system, whether it is an outside connector, and its
# connection type ("Equality", "Flow", "Stream", "InputVar", "OutputVar").
.cvertex <- function(namespace, var = NULL, isouter = TRUE, type = NULL) {
  if (is.null(var)) return(list(name = namespace, isouter = TRUE, type = "Flow"))
  nm <- if (inherits(var, "System")) nameof(var) else .name(var)
  if (is.null(type)) type <- if (inherits(var, "System")) "Equality" else .connect_type(var)
  list(name = c(namespace, .namespace_hierarchy(nm)), isouter = isouter, type = type)
}

.cv_key <- function(v) paste(paste(v$name, collapse = NAMESPACE_SEPARATOR), v$type, v$isouter, sep = "|")

.hg_new <- function() {
  g <- new.env(parent = emptyenv())
  g$labels <- new.env(parent = emptyenv())
  g$invmap <- list()
  g$edges <- list()
  g$inputs <- integer()
  g
}

.hg_add_vertex <- function(g, v) {
  k <- .cv_key(v)
  j <- g$labels[[k]]
  if (!is.null(j)) return(j)
  j <- length(g$invmap) + 1L
  g$invmap[[j]] <- v
  assign(k, j, envir = g$labels)
  j
}

# Add a hyperedge; the first vertex is recorded as the input, which only the
# negative graph uses.
.hg_add_edge <- function(g, verts) {
  ids <- vapply(verts, function(v) .hg_add_vertex(g, v), 0L)
  g$edges[[length(g$edges) + 1L]] <- unique(ids)
  g$inputs <- c(g$inputs, ids[1L])
  invisible(length(g$edges))
}

# Connected components by union find, ordered by the first vertex of each
# set, as Graphs.connected_components(::HyperGraph).
.hg_components <- function(g) {
  n <- length(g$invmap)
  if (!n) return(list())
  parent <- seq_len(n)
  root_of <- function(i) {
    while (parent[i] != i) i <- parent[i]
    i
  }
  for (e in g$edges) {
    if (length(e) < 2L) next
    r <- root_of(e[1L])
    for (k in e[-1L]) {
      rk <- root_of(k)
      if (rk != r) parent[rk] <- r
    }
  }
  roots <- vapply(seq_len(n), root_of, 0L)
  lapply(unique(roots), function(r) g$invmap[roots == r])
}

# remove_negative_connections!: for each negative (causal) hyperedge whose
# vertices are all in `g`, remove its non-input vertices from every hyperedge
# of `g` that holds its input and at least one other of its vertices.
.hg_remove_negative <- function(g, neg) {
  if (!length(neg$edges)) return(invisible())
  edgesets <- g$edges
  rm <- vector("list", length(edgesets))
  for (i in seq_along(neg$edges)) {
    neg_verts <- neg$invmap[neg$edges[[i]]]
    jdx <- vapply(neg_verts, function(v) g$labels[[.cv_key(v)]] %||% 0L, 0L)
    if (any(jdx == 0L)) next
    input_j <- g$labels[[.cv_key(neg$invmap[[neg$inputs[i]]])]]
    others <- setdiff(jdx, input_j)
    for (ej in seq_along(edgesets)) {
      es <- edgesets[[ej]]
      if (!(input_j %in% es) || !any(others %in% es)) next
      rm[[ej]] <- c(rm[[ej]], others)
    }
  }
  for (ej in seq_along(rm)) if (length(rm[[ej]])) g$edges[[ej]] <- setdiff(g$edges[[ej]], rm[[ej]])
  invisible()
}

# ---------------------------------------------------------------------------
# Building connection sets

# Names of the connectors that are direct subsystems of `sys`.
.generate_isouter <- function(sys) {
  outer <- character()
  for (s in .f(sys, "systems")) if (isconnector(s)) outer <- c(outer, nameof(s))
  function(x) {
    if (is.character(x)) return(x %in% outer)
    if (!isconnector(x)) stop(nameof(x), " is not a connector!", call. = FALSE)
    .namespace_hierarchy(nameof(x))[1L] %in% outer
  }
}

.connection_error <- function(ss) {
  stop("Different types of connectors are in one connection statement: <",
       paste(vapply(ss, nameof, ""), collapse = ", "), ">", call. = FALSE)
}

.cset_state <- function() list(cg = .hg_new(), dg = .hg_new(), negative = FALSE)

.handle_connection <- function(state, conn, namespace, isouter) {
  if (isTRUE(conn$domain)) {
    edge <- lapply(conn$systems, function(s) .cvertex(c(namespace, .namespace_hierarchy(nameof(s)))))
    .hg_add_edge(state$dg, edge)
    return(invisible())
  }
  conn_items <- conn$systems
  if (length(conn_items) && !inherits(conn_items[[1L]], "System"))
    return(.csets_from_vars(state, namespace, conn_items, isouter))
  .csets_from_systems(state, namespace, conn_items, isouter)
}

.csets_from_vars <- function(state, namespace, vars, isouter) {
  edge <- lapply(vars, function(v) {
    ns <- .namespace_hierarchy(.name(v))
    i <- isinput(v); o <- isoutput(v)
    names <- paste(vapply(vars, format, ""), collapse = ", ")
    if (i && o) stop("Variable ", format(v), " in connection `connect(", names, ")` is both input and output.", call. = FALSE)
    if (!i && !o) stop("Variable ", format(v), " in connection `connect(", names, ")` is neither input nor output.", call. = FALSE)
    list(name = c(namespace, ns), isouter = length(ns) == 1L || isouter(ns[1L]),
         type = if (i) "InputVar" else "OutputVar")
  })
  .hg_add_edge(state$cg, edge)
  # removed analysis points also remove the Equality variant of the edge
  if (state$negative) .hg_add_edge(state$cg, lapply(edge, function(v) { v$type <- "Equality"; v }))
  invisible()
}

.csets_from_systems <- function(state, namespace, systems, isouter) {
  regular <- list()
  domain_system <- NULL
  for (s in systems) {
    if (.is_domain_connector(s)) {
      if (!is.null(domain_system))
        stop("connect(", paste(vapply(systems, nameof, ""), collapse = ","),
             ") contains multiple source domain connectors. There can only be one!", call. = FALSE)
      domain_system <- s
    } else regular[[length(regular) + 1L]] <- s
  }
  if (!length(regular)) stop("a connection needs a connector that is not a domain source", call. = FALSE)
  if (!is.null(domain_system)) {
    dv <- unknowns(domain_system)[[1L]]
    dns <- c(namespace, .namespace_hierarchy(nameof(domain_system)))
    edge <- list(.cvertex(dns, dv, FALSE))
    dedge <- list(.cvertex(dns))
    for (s in regular) {
      sns <- c(namespace, .namespace_hierarchy(nameof(s)))
      so <- isouter(s)
      for (v in unknowns(s)) {
        if (!(.connect_type(v) == "Flow" && identical(.key(v), .key(dv)))) next
        edge[[length(edge) + 1L]] <- .cvertex(sns, v, so)
        dedge[[length(dedge) + 1L]] <- .cvertex(sns)
      }
    }
    .hg_add_edge(state$cg, edge)
    .hg_add_edge(state$dg, dedge)
    return(invisible())
  }
  # ponytail: multibody frames (IsFrame orientation variables) are not handled;
  # add when a frame library exists.
  dvs1 <- unknowns(regular[[1L]])
  keys1 <- .keys(dvs1)
  var_sets <- rep(list(list()), length(dvs1))
  dedge <- list()
  for (i in seq_along(regular)) {
    s <- regular[[i]]
    uv <- unknowns(s)
    if (i != 1L && (length(uv) != length(dvs1) || !all(.keys(uv) %in% keys1))) .connection_error(regular)
    sns <- c(namespace, .namespace_hierarchy(nameof(s)))
    so <- isouter(s)
    for (j in seq_along(uv)) var_sets[[j]][[length(var_sets[[j]]) + 1L]] <- .cvertex(sns, uv[[j]], so)
    dedge[[length(dedge) + 1L]] <- .cvertex(sns)
  }
  for (vs in var_sets) {
    if (length(unique(vapply(vs, `[[`, "", "type"))) > 1L) .connection_error(regular)
    .hg_add_edge(state$cg, vs)
  }
  .hg_add_edge(state$dg, dedge)
  invisible()
}

# Preorder walk (generate_connection_set!): collects connection sets, removes
# connection equations and records ignored connections in the negative state.
.gen_cset_rec <- function(state, neg, sys, namespace) {
  if (isTRUE(.f(sys, "namespacing"))) namespace <- c(namespace, nameof(sys))
  subsys <- .f(sys, "systems")
  isouter <- .generate_isouter(sys)
  eqs <- list()
  for (eq in .f(sys, "eqs")) {
    if (inherits(eq, "Connection")) .handle_connection(state, eq, namespace, isouter)
    else eqs[[length(eqs) + 1L]] <- eq
  }
  for (conn in .f(sys, "ignored_connections")) .handle_connection(neg, conn, namespace, isouter)
  # every flow variable of a connector gets at least its inside singleton set
  for (s in subsys) {
    if (!isconnector(s) || .is_domain_connector(s)) next
    for (v in unknowns(s)) {
      if (.connect_type(v) != "Flow") next
      .hg_add_edge(state$cg, list(.cvertex(c(namespace, nameof(s)), v, FALSE)))
    }
  }
  new_systems <- lapply(subsys, function(s) .gen_cset_rec(state, neg, s, namespace))
  .update_system(sys, systems = new_systems, eqs = eqs)
}

.generate_connection_set <- function(sys) {
  state <- .cset_state()
  neg <- .cset_state()
  neg$negative <- TRUE
  sys <- toggle_namespacing(sys, FALSE)
  sys <- .gen_cset_rec(state, neg, sys, character())
  .hg_remove_negative(state$cg, neg$cg)
  .hg_remove_negative(state$dg, neg$dg)
  list(sys = sys, csets = .hg_components(state$cg), domain_csets = .hg_components(state$dg))
}

# The variable or subsystem a vertex refers to (iterative_getproperty).
.variable_from_vertex <- function(sys, cvert) {
  obj <- toggle_namespacing(sys, FALSE)
  for (nm in cvert$name) obj <- .getproperty(obj, nm)
  obj
}

# ---------------------------------------------------------------------------
# Connection equations

.generate_connection_equations <- function(sys, csets) {
  eqs <- list()
  stream_csets <- list()
  add <- function(lhs, rhs) eqs[[length(eqs) + 1L]] <<- Equation(lhs, rhs)
  for (cset in csets) {
    vtype <- cset[[1L]]$type
    if (vtype %in% c("InputVar", "OutputVar")) {
      if (length(cset) <= 1L) next
      inner_output <- NULL; outer_input <- NULL
      for (cv in cset) {
        if (cv$isouter && cv$type == "InputVar") {
          if (!is.null(outer_input)) stop("Found two outer input connectors in the same connection set.", call. = FALSE)
          outer_input <- cv
        } else if (!cv$isouter && cv$type == "OutputVar") {
          if (!is.null(inner_output)) stop("Found two inner output connectors in the same connection set.", call. = FALSE)
          inner_output <- cv
        }
      }
      root <- inner_output %||% outer_input
      if (is.null(root))
        stop("Found a causal connection set with no source. For a causal connection set to be valid ",
             "it must have one outer input or one inner output. Problematic connection set: ",
             paste(vapply(cset, function(v) paste(v$name, collapse = "."), ""), collapse = ", "), call. = FALSE)
      root_var <- .variable_from_vertex(sys, root)
      rk <- .cv_key(root)
      for (cv in cset) if (.cv_key(cv) != rk) add(.variable_from_vertex(sys, cv), root_var)
    } else if (vtype == "Stream") {
      stream_csets[[length(stream_csets) + 1L]] <- cset
    } else if (vtype == "Flow") {
      terms <- lapply(cset, function(cv) {
        v <- .variable_from_vertex(sys, cv)
        if (cv$isouter) -v else v
      })
      add(.S0(), Reduce(`+`, terms))
    } else {
      vars <- lapply(cset, function(cv) .variable_from_vertex(sys, cv))
      outer_input <- 0L; inner_output <- 0L
      all_io <- TRUE
      for (i in seq_along(vars)) {
        is_i <- isinput(vars[[i]]); is_o <- isoutput(vars[[i]])
        all_io <- all_io && (is_i || is_o)
        if (!all_io) break
        if (cset[[i]]$isouter && is_i && outer_input == 0L) outer_input <- i
        else if (!cset[[i]]$isouter && is_o && inner_output == 0L) inner_output <- i
      }
      if (all_io && xor(outer_input != 0L, inner_output != 0L)) {
        r <- if (outer_input == 0L) inner_output else outer_input
        for (i in seq_along(vars)) if (i != r) add(vars[[i]], vars[[r]])
      } else {
        for (i in seq_along(vars)[-1L]) add(vars[[1L]], vars[[i]])
      }
    }
  }
  list(eqs = eqs, stream_csets = stream_csets)
}

# The flow variable of the connector holding the stream variable of `cv`.
.get_flowvar <- function(sys, cv) {
  par <- .variable_from_vertex(sys, list(name = utils::head(cv$name, -1L)))
  for (v in unknowns(par)) if (.connect_type(v) == "Flow") return(renamespace(par, v))
  stop("There is no flow variable in system `", nameof(par), "`", call. = FALSE)
}

.positivemax_num <- function(m, si) {
  eps <- 1e-4
  alpha <- ifelse(si > eps, 1, ifelse(si > 0, (si / eps)^2 * (3 - 2 * si / eps), 0))
  alpha * pmax(m, 0) + (1 - alpha) * eps
}

# Runtime of the mixing rule of instream (instream_rt). As in the Julia
# source, the outer flows enter positivemax with a negative sign.
.instream_rt_num <- function(inner_n, outer_n, ...) {
  v <- c(...)
  inner_n <- as.integer(inner_n); outer_n <- as.integer(outer_n)
  inner_f <- v[seq_len(inner_n)]
  inner_s <- v[inner_n + seq_len(inner_n)]
  off <- 2L * inner_n
  outer_f <- v[off + seq_len(outer_n)]
  outer_s <- v[off + outer_n + seq_len(outer_n)]
  si <- sum(pmax(-inner_f, 0)) + sum(pmax(outer_f, 0))
  wi <- .positivemax_num(-inner_f, si)
  wo <- .positivemax_num(-outer_f, si)
  (sum(wi * inner_s) + sum(wo * outer_s)) / (sum(wi) + sum(wo))
}

.stream_funs <- function() {
  list(positivemax = register_symbolic("_positivemax", .positivemax_num),
       instream_rt = register_symbolic("instream_rt", .instream_rt_num))
}

.collect_instream <- function(eqs) {
  out <- list()
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    if (.is_op(z, "instream")) out[[.key(z)]] <<- z
    for (a in .args(z)) walk(a)
  }
  for (e in eqs) if (inherits(e, "Equation")) { walk(e$lhs); walk(e$rhs) }
  out
}

# expand_instream: equations and instream substitutions for stream sets,
# after the Modelica specification 3.6, section 15.2.
.expand_instream <- function(csets, sys, tol = 1e-8) {
  from <- list(); to <- list()
  eqs <- list()
  sub <- function(a, b) { from[[length(from) + 1L]] <<- instream(a); to[[length(to) + 1L]] <<- b }
  if (!length(csets)) return(list(eqs = eqs, rule = .rule(from, to)))
  fns <- .stream_funs()
  for (cset in csets) {
    is_out <- vapply(cset, `[[`, TRUE, "isouter")
    n_outer <- sum(is_out); n_inner <- length(cset) - n_outer
    sv <- lapply(cset, function(cv) .variable_from_vertex(sys, cv))
    if (n_inner == 1L && n_outer == 0L) {
      sub(sv[[1L]], sv[[1L]])
    } else if (n_inner == 2L && n_outer == 0L) {
      sub(sv[[1L]], sv[[2L]])
      sub(sv[[2L]], sv[[1L]])
    } else if (n_inner == 1L && n_outer == 1L) {
      inner <- sv[[which(!is_out)]]; outer <- sv[[which(is_out)]]
      sub(inner, instream(outer))
      eqs[[length(eqs) + 1L]] <- Equation(outer, inner)
    } else if (n_inner == 0L && n_outer == 2L) {
      eqs[[length(eqs) + 1L]] <- Equation(sv[[1L]], instream(sv[[2L]]))
      eqs[[length(eqs) + 1L]] <- Equation(sv[[2L]], instream(sv[[1L]]))
    } else {
      fv <- lapply(cset, function(cv) .get_flowvar(sys, cv))
      inner_s <- sv[!is_out]; outer_s <- sv[is_out]
      inner_f <- fv[!is_out]; outer_f <- fv[is_out]
      for (k in seq_along(inner_s)) {
        args <- c(list(as_sym(n_inner - 1L), as_sym(n_outer)), inner_f[-k], inner_s[-k], outer_f, outer_s)
        sub(inner_s[[k]], do.call(fns$instream_rt, args))
      }
      for (q in seq_len(n_outer)) {
        # As in the Julia source, positivemax(m, sq; tol) evaluates
        # _positivemax(m, tol); sq itself is unused.
        pm <- function(m) fns$positivemax(m, as_sym(tol))
        num <- c(lapply(seq_along(inner_f), function(i) pm(-inner_f[[i]]) * inner_s[[i]]),
                 lapply(seq_along(outer_f)[-q], function(i) pm(outer_f[[i]]) * instream(outer_s[[i]])))
        den <- c(lapply(inner_f, function(f) pm(-f)), lapply(outer_f[-q], pm))
        eqs[[length(eqs) + 1L]] <- Equation(outer_s[[q]], Reduce(`+`, num) / Reduce(`+`, den))
      }
    }
  }
  list(eqs = eqs, rule = .rule(from, to))
}

# get_domain_bindings: parameters of connectors in a domain network are bound
# to the parameters of the same name of its domain source.
.get_domain_bindings <- function(sys, domain_csets) {
  binds <- symmap()
  for (cset in domain_csets) {
    systems <- lapply(cset, function(cv) .variable_from_vertex(sys, cv))
    idx <- which(vapply(systems, .is_domain_connector, TRUE))
    if (!length(idx)) next
    dsys <- systems[[idx[1L]]]
    dbinds <- bindings(dsys); dics <- initial_conditions(dsys)
    for (j in seq_along(systems)) {
      if (j == idx[1L]) next
      csys <- systems[[j]]
      if (.is_domain_connector(csys))
        stop("Domain sources ", nameof(dsys), " and ", nameof(csys), " are connected!", call. = FALSE)
      for (par in .parameters_of(csys)) {
        if (!smap_has(dbinds, par) && !smap_has(dics, par)) next
        binds <- smap_set(binds, renamespace(csys, par), renamespace(dsys, par))
      }
    }
  }
  binds
}

# The expand_connections hook.
.expand_connections <- function(sys, tol = 1e-10, ...) {
  if (!.has_connections(sys)) return(sys)
  .warn_unbalanced_connectors(sys)
  ns <- .f(sys, "namespacing")
  sys <- remove_analysis_points(sys)
  gc <- .generate_connection_set(sys)
  sys <- gc$sys
  ce <- .generate_connection_equations(sys, gc$csets)
  eqs <- equations(sys)
  st <- .expand_instream(ce$stream_csets, sys, tol = tol)
  eqs <- c(eqs, ce$eqs, st$eqs)
  if (length(st$rule)) {
    it <- max(length(st$rule), 10L)
    eqs <- lapply(eqs, function(e) Equation(.subst_fixpoint(e$lhs, st$rule, it), .subst_fixpoint(e$rhs, st$rule, it)))
  }
  dbinds <- .get_domain_bindings(sys, gc$domain_csets)
  out <- flatten(sys)
  binds <- .f(out, "bindings") %||% symmap()
  common <- intersect(names(dbinds), names(binds))
  if (length(common)) stop("Cannot merge without overriding: common key ", binds[[common[1L]]]$var, call. = FALSE)
  .update_system(out, eqs = eqs, bindings = smap_merge(binds, dbinds), namespacing = ns,
                 ignored_connections = NULL)
}

# Whether any system in the hierarchy holds a connection, an analysis point,
# an ignored connection or a connector (whose flows default to zero).
.has_connections <- function(sys) {
  if (length(.f(sys, "ignored_connections"))) return(TRUE)
  for (e in .f(sys, "eqs")) if (inherits(e, "Connection")) return(TRUE)
  for (s in .f(sys, "systems")) if (isconnector(s) || .has_connections(s)) return(TRUE)
  FALSE
}

.mtk$hooks$expand_connections <- .expand_connections

# Namespacing of connections and analysis points inside equation lists.
.mtk$hooks$namespace_connection <- function(x, prefix) {
  ns1 <- function(s) if (inherits(s, "System")) .rename_system(s, .ns_name(prefix, nameof(s))) else .namespace(s, prefix)
  if (inherits(x, "AnalysisPoint")) {
    x$input <- if (is.null(x$input)) NULL else ns1(x$input)
    x$outputs <- if (is.null(x$outputs)) NULL else lapply(x$outputs, ns1)
    x$name <- .ns_name(prefix, x$name)
    return(x)
  }
  x$systems <- lapply(x$systems, ns1)
  x
}
