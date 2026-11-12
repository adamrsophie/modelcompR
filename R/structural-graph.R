# Bipartite graphs, matchings and strongly connected components, after
# BipartiteGraphs.jl (MIT): src/bipartite_graph.jl (BipartiteGraph, Matching),
# src/maximal_matching.jl (construct_augmenting_path!, maximal_matching),
# src/dicmobigraph.jl (DiCMOBiGraph neighbor rules) and
# src/condensation_graphs.jl (MatchedCondensationGraph); find_var_sccs of
# ModelingToolkit.jl v10 (MIT) src/structural_transformation/utils.jl; and the
# dependency graphs of ModelingToolkitBase/src/systems/dependency_graphs.jl
# (equation_dependencies, asgraph, variable_dependencies, asdigraph,
# eqeq_dependencies, varvar_dependencies).
#
# A BipartiteGraph is a mutable environment: `fadjlist[[e]]` are the sorted
# destination vertices (variables) of source vertex (equation) e, and
# `badjlist[[v]]` the sorted sources of destination v. A matching is an
# integer vector indexed by destination: the matched source, 0L for unassigned
# and -1L for a selected state (StateSelection.SelectedState).

.UNASSIGNED <- 0L
.SELECTED_STATE <- -1L

#' Bipartite graphs of equations and variables
#'
#' `BipartiteGraph(nsrcs, ndsts)` builds an empty bipartite graph with `nsrcs`
#' source vertices (equations) and `ndsts` destination vertices (variables),
#' as the type of the same name in ModelingToolkit. `fadjlist` gives the
#' neighbors of each source instead. The result is a mutable object with
#' fields `fadjlist` (variables of each equation) and `badjlist` (equations of
#' each variable). [asgraph()] builds one from a system.
#'
#' @param nsrcs,ndsts Numbers of source and destination vertices.
#' @param fadjlist A list of integer vectors: the destinations of each source.
#' @return An object of class `BipartiteGraph`.
#' @examples
#' g <- BipartiteGraph(fadjlist = list(c(1L, 2L), 2L), ndsts = 3L)
#' g$fadjlist
#' g$badjlist
#' @export
BipartiteGraph <- function(nsrcs = 0L, ndsts = 0L, fadjlist = NULL) {
  g <- new.env(parent = emptyenv())
  if (!is.null(fadjlist)) {
    nsrcs <- length(fadjlist)
    ndsts <- max(as.integer(ndsts), 0L, unlist(fadjlist))
  }
  g$fadjlist <- rep(list(integer()), nsrcs)
  g$badjlist <- rep(list(integer()), ndsts)
  class(g) <- "BipartiteGraph"
  for (e in seq_along(fadjlist)) .bg_set_neighbors(g, e, fadjlist[[e]])
  g
}

#' @export
print.BipartiteGraph <- function(x, ...) {
  cat("BipartiteGraph with (", length(x$fadjlist), ", ", length(x$badjlist), ") (sources, destinations) and ",
      sum(lengths(x$fadjlist)), " edges\n", sep = "")
  for (e in seq_along(x$fadjlist)) cat("  ", e, ": ", paste(x$fadjlist[[e]], collapse = " "), "\n", sep = "")
  invisible(x)
}

.bg_nsrcs <- function(g) length(g$fadjlist)
.bg_ndsts <- function(g) length(g$badjlist)

.bg_copy <- function(g) {
  if (is.null(g)) return(NULL)
  h <- new.env(parent = emptyenv())
  h$fadjlist <- g$fadjlist
  h$badjlist <- g$badjlist
  class(h) <- "BipartiteGraph"
  h
}

.bg_has_edge <- function(g, e, v) v %in% g$fadjlist[[e]]

.bg_add_edge <- function(g, e, v) {
  if (v %in% g$fadjlist[[e]]) return(invisible(FALSE))
  g$fadjlist[[e]] <- sort.int(c(g$fadjlist[[e]], v))
  g$badjlist[[v]] <- sort.int(c(g$badjlist[[v]], e))
  invisible(TRUE)
}

.bg_rem_edge <- function(g, e, v) {
  g$fadjlist[[e]] <- setdiff(g$fadjlist[[e]], v)
  g$badjlist[[v]] <- setdiff(g$badjlist[[v]], e)
  invisible(TRUE)
}

# Replace the variables of equation e (set_neighbors! on the graph).
.bg_set_neighbors <- function(g, e, vars) {
  vars <- sort.int(unique(as.integer(vars)))
  old <- g$fadjlist[[e]]
  for (v in setdiff(old, vars)) g$badjlist[[v]] <- setdiff(g$badjlist[[v]], e)
  for (v in setdiff(vars, old)) g$badjlist[[v]] <- sort.int(c(g$badjlist[[v]], e))
  g$fadjlist[[e]] <- vars
  invisible(g)
}

# Replace the equations of variable v (set_neighbors! on invview(graph)).
.bg_set_dneighbors <- function(g, v, eqs) {
  eqs <- sort.int(unique(as.integer(eqs)))
  old <- g$badjlist[[v]]
  for (e in setdiff(old, eqs)) g$fadjlist[[e]] <- setdiff(g$fadjlist[[e]], v)
  for (e in setdiff(eqs, old)) g$fadjlist[[e]] <- sort.int(c(g$fadjlist[[e]], v))
  g$badjlist[[v]] <- eqs
  invisible(g)
}

.bg_add_src <- function(g) {
  g$fadjlist[[length(g$fadjlist) + 1L]] <- integer()
  length(g$fadjlist)
}

.bg_add_dst <- function(g) {
  g$badjlist[[length(g$badjlist) + 1L]] <- integer()
  length(g$badjlist)
}

# Inverse of a matching: source -> destination, 0L when unmatched.
.inv_matching <- function(m, nsrcs) {
  inv <- integer(nsrcs)
  hit <- which(m > 0L)
  inv[m[hit]] <- hit
  inv
}

# construct_augmenting_path!: `st` is an environment holding the matching
# `m` (indexed by the far side), the far side colors `dcolor` and optionally
# the near side colors `scolor`. `adj` lists the far side neighbors of each
# near side vertex, so the same code matches equations to variables
# (adj = fadjlist) and variables to equations (adj = badjlist). `ok` is the
# logical filter on the far side.
.augment <- function(st, adj, vsrc, ok) {
  if (!is.null(st$scolor)) st$scolor[vsrc] <- TRUE
  nb <- adj[[vsrc]]
  for (vdst in nb) {
    if (ok[vdst] && st$m[vdst] == .UNASSIGNED) {
      st$m[vdst] <- vsrc
      return(TRUE)
    }
  }
  for (vdst in nb) {
    if (!ok[vdst] || st$dcolor[vdst]) next
    st$dcolor[vdst] <- TRUE
    msrc <- st$m[vdst]
    if (msrc > 0L && .augment(st, adj, msrc, ok)) {
      st$m[vdst] <- vsrc
      return(TRUE)
    }
  }
  FALSE
}

# maximal_matching(g; srcfilter, dstfilter): a matching of variables to
# equations of length max(nsrcs, ndsts).
.maximal_matching <- function(g, srcfilter = NULL, dstfilter = NULL) {
  ns <- .bg_nsrcs(g); nd <- .bg_ndsts(g)
  st <- new.env(parent = emptyenv())
  st$m <- integer(max(ns, nd))
  ok <- if (is.null(dstfilter)) rep(TRUE, nd) else dstfilter
  for (e in seq_len(ns)) {
    if (!is.null(srcfilter) && !srcfilter[e]) next
    st$dcolor <- logical(nd)
    .augment(st, g$fadjlist, e, ok)
  }
  st$m
}

# Strongly connected components of a directed graph given by a successor
# function, with Tarjan's algorithm. Components are returned in the order
# they complete, which is reverse topological order.
# ponytail: recursive depth first search, limited by R's expression depth
# (chains of a few thousand vertices); make it iterative if that is reached.
.tarjan_scc <- function(n, succ) {
  index <- integer(n); low <- integer(n); onstack <- logical(n)
  stack <- integer(n); sp <- 0L; counter <- 0L
  out <- list()
  visit <- function(v) {
    counter <<- counter + 1L
    index[v] <<- counter; low[v] <<- counter
    sp <<- sp + 1L; stack[sp] <<- v; onstack[v] <<- TRUE
    for (w in succ(v)) {
      if (index[w] == 0L) {
        visit(w)
        if (low[w] < low[v]) low[v] <<- low[w]
      } else if (onstack[w] && index[w] < low[v]) {
        low[v] <<- index[w]
      }
    }
    if (low[v] == index[v]) {
      comp <- integer()
      repeat {
        w <- stack[sp]; sp <<- sp - 1L
        onstack[w] <<- FALSE
        comp <- c(comp, w)
        if (w == v) break
      }
      out[[length(out) + 1L]] <<- comp
    }
  }
  for (v in seq_len(n)) if (index[v] == 0L) visit(v)
  out
}

# Topological order of a directed acyclic graph, exactly as
# Graphs.topological_sort_by_dfs (the order of ties follows its stack).
.toposort_dfs <- function(n, succ) {
  color <- integer(n)   # 0 new, 1 on the current path, 2 done
  verts <- integer()
  for (v in seq_len(n)) {
    if (color[v] != 0L) next
    stack <- v
    while (length(stack)) {
      u <- stack[length(stack)]
      if (color[u] == 0L) {
        color[u] <- 1L
        for (w in succ(u)) {
          if (color[w] == 1L) stop("the graph has a cycle", call. = FALSE)
          if (color[w] == 0L) stack <- c(stack, w)
        }
      } else {
        stack <- stack[-length(stack)]
        if (color[u] == 1L) {
          color[u] <- 2L
          verts <- c(u, verts)
        }
      }
    }
  }
  verts
}

# Successors of variable v in DiCMOBiGraph{true}(g, m): the variables matched
# to the equations v appears in. An edge v -> w says w is computed from v.
.cmo_var_succ <- function(g, m) {
  function(v) {
    es <- g$badjlist[[v]]
    if (!length(es)) return(integer())
    inv <- attr(m, "inv")
    w <- inv[es]
    w[w > 0L & w != v]
  }
}

# find_var_sccs(g, assign): strongly connected components of the variables
# under the orientation given by the matching, in dependency order (each
# component depends only on earlier ones), each sorted.
.find_var_sccs <- function(g, m) {
  nd <- .bg_ndsts(g)
  m <- m[seq_len(nd)]
  m[m < 0L] <- 0L
  attr(m, "inv") <- .inv_matching(m, .bg_nsrcs(g))
  succ <- .cmo_var_succ(g, m)
  sccs <- .tarjan_scc(nd, succ)
  if (!length(sccs)) return(list())
  comp <- integer(nd)
  for (i in seq_along(sccs)) comp[sccs[[i]]] <- i
  csucc <- function(cc) {
    out <- integer()
    for (v in sccs[[cc]]) out <- c(out, comp[succ(v)])
    out[out != cc]
  }
  lapply(sccs[.toposort_dfs(length(sccs), csucc)], sort.int)
}

# ---------------------------------------------------------------------------
# Dependency graphs (ModelingToolkitBase dependency_graphs.jl)

.dep_vars <- function(target, variables) {
  vs <- .collect_vars(target)
  vs <- lapply(vs, function(v) if (.type(v) == "Derivative") .args(v)[[1L]] else v)
  keys <- unique(.keys(vs))
  variables[.keys(variables) %in% keys]
}

#' Dependency graphs of a system
#'
#' `equation_dependencies(sys)` lists, for each equation, the variables its
#' right-hand side depends on. `variable_dependencies(sys)` is the bipartite
#' graph from variables to the equations that modify them (the variables on
#' the left-hand sides). `asgraph()` turns the equation dependencies into a
#' [BipartiteGraph()] from equations to variables, `asdigraph()` turns a
#' bipartite graph into a directed graph on the joint vertex set, and
#' `eqeq_dependencies()` and `varvar_dependencies()` give the directed graphs
#' of equations (or variables) that depend on each other. These are the
#' functions of the same names in ModelingToolkit.
#'
#' @param sys A `System`.
#' @param variables The variables to track; defaults to the unknowns.
#' @param eqs The equations (or jumps); defaults to the system's equations.
#' @param eqdeps,vardeps Bipartite graphs from `asgraph()` and
#'   `variable_dependencies()`. `asgraph()` also accepts the list returned by
#'   `equation_dependencies()` as `sys`, with `variables` the variable list.
#' @param g A `BipartiteGraph`.
#' @param equationsfirst Number the equations before the variables.
#' @return A list of variable lists, a `BipartiteGraph`, or a directed graph:
#'   a list with `fadjlist` (successors) and `badjlist` (predecessors).
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)", "y(t)"); parameters("k")
#' sys <- System(list(D(x) ~ -k * y, D(y) ~ x), t, name = "sys")
#' equation_dependencies(sys)
#' eg <- asgraph(sys)
#' vg <- variable_dependencies(sys)
#' eqeq_dependencies(eg, vg)$fadjlist
#' varvar_dependencies(eg, vg)$fadjlist
#' asdigraph(eg, sys)$fadjlist
#' @export
equation_dependencies <- function(sys, variables = unknowns(sys), eqs = equations(sys)) {
  lapply(eqs, function(eq) {
    target <- if (inherits(eq, "Equation")) eq$rhs else if (is_sym(eq)) eq else eq$rate %||% eq
    .dep_vars(target, variables)
  })
}

#' @rdname equation_dependencies
#' @export
asgraph <- function(sys, variables = unknowns(sys), eqs = equations(sys)) {
  eqdeps <- if (inherits(sys, "System")) equation_dependencies(sys, variables, eqs) else sys
  vk <- .keys(variables)
  BipartiteGraph(fadjlist = lapply(eqdeps, function(dep) match(.keys(dep), vk)),
                 ndsts = length(variables))
}

#' @rdname equation_dependencies
#' @export
variable_dependencies <- function(sys, variables = unknowns(sys), eqs = equations(sys)) {
  modified <- lapply(eqs, function(eq) {
    target <- if (inherits(eq, "Equation")) list(eq$lhs)
              else lapply(.as_equations(eq$affect %||% list()), `[[`, "lhs")
    match(.keys(.dep_vars(target, variables)), .keys(variables))
  })
  # sources are the variables here, destinations the equations
  g <- BipartiteGraph(length(variables), length(eqs))
  for (e in seq_along(modified)) for (v in modified[[e]]) .bg_add_edge(g, v, e)
  g
}

.digraph <- function(fadjlist, badjlist) {
  structure(list(fadjlist = fadjlist, badjlist = badjlist, ne = sum(lengths(fadjlist))),
            class = "mtk_digraph")
}

#' @rdname equation_dependencies
#' @export
asdigraph <- function(g, sys, variables = unknowns(sys), equationsfirst = TRUE, eqs = equations(sys)) {
  neqs <- length(eqs); nvars <- length(variables)
  offset <- if (equationsfirst) neqs else nvars
  fadj <- c(lapply(g$fadjlist, function(x) x + offset), rep(list(integer()), length(g$badjlist)))
  badj <- c(rep(list(integer()), length(g$fadjlist)), g$badjlist)
  .digraph(fadj, badj)
}

#' @rdname equation_dependencies
#' @export
eqeq_dependencies <- function(eqdeps, vardeps) {
  n <- length(eqdeps$fadjlist)
  fadj <- rep(list(integer()), n)
  for (eqidx in seq_along(vardeps$badjlist)) {
    for (sidx in vardeps$badjlist[[eqidx]])
      fadj[[eqidx]] <- union(fadj[[eqidx]], eqdeps$badjlist[[sidx]])
  }
  fadj <- lapply(fadj, sort.int)
  badj <- rep(list(integer()), n)
  for (i in seq_len(n)) for (j in fadj[[i]]) badj[[j]] <- c(badj[[j]], i)
  .digraph(fadj, badj)
}

#' @rdname equation_dependencies
#' @export
varvar_dependencies <- function(eqdeps, vardeps) eqeq_dependencies(vardeps, eqdeps)
