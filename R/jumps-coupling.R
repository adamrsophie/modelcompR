# Coupled jump problems and jump utilities, after JumpProcesses.jl:
# coupling.jl (SplitCoupledJumpProblem, cat_problems, build_split_jumps),
# jumps.jl (JumpSet, RateBounds, get_num_majumps), jump_dependencies.jl
# (JumpDependencies), aggregators/aggregators.jl (needs_depgraph,
# needs_vartojumps_map) and isolated_copy.jl.

#' Coupled jump problems and jump utilities
#'
#' `SplitCoupledJumpProblem(prob, prob_control, aggregator, coupling_map)`
#' simulates two jump problems with common random numbers, the split
#' coupling of Anderson and Higham: each pair `c(i, j)` of `coupling_map`
#' couples jump `i` of `prob` with jump `j` of `prob_control` (they fire
#' together at the minimum of their rates), the other jumps fire alone. The
#' state is `c(u, u_control)`. Both problems are built from R functions with
#' [ConstantRateJump()]s; their inner problems are discrete, ODE or SDE
#' problems.
#'
#' `JumpSet(...)` groups jumps; `JumpProblem()` takes it like a list.
#' `RateBounds(lrate, urate, rateinterval)` is a bound of a rate over an
#' interval of states or time (the return value of the `bounds` function of
#' a [ConstantRateJump()] for [RSSA()]). `JumpDependencies()` holds the
#' dependency maps passed to `JumpProblem(dependencies = )`.
#' `get_num_majumps()` counts mass action reactions; `needs_depgraph()` and
#' `needs_vartojumps_map()` tell which maps an aggregator uses.
#' `isolated_copy(jprob)` returns a problem whose solves share no state with
#' `jprob`; R problems are values, so this is `jprob` itself.
#'
#' @param prob,prob_control Jump problems.
#' @param aggregator A jump aggregator.
#' @param coupling_map List of index pairs `c(i, j)`.
#' @param ... Jumps, or for `SplitCoupledJumpProblem()` options of
#'   [JumpProblem()].
#' @param lrate,urate Lower and upper rate bounds.
#' @param rateinterval Time over which the bounds hold.
#' @param dep_graph,vartojumps_map,jumptovars_map Dependency maps.
#' @param jumps Mass action jumps.
#' @param agg An aggregator.
#' @param jprob A jump problem.
#' @return A jump problem, a jump set, bounds, maps, a count or a flag.
#' @examples
#' birth <- ConstantRateJump(function(u, p, t) 10, function(integ) integ$u[1] <- integ$u[1] + 1)
#' death <- ConstantRateJump(function(u, p, t) u[1], function(integ) integ$u[1] <- integ$u[1] - 1)
#' death2 <- ConstantRateJump(function(u, p, t) 1.1 * u[1], function(integ) integ$u[1] <- integ$u[1] - 1)
#' j1 <- JumpProblem(DiscreteProblem(0, c(0, 5)), Direct(), birth, death)
#' j2 <- JumpProblem(DiscreteProblem(0, c(0, 5)), Direct(), birth, death2)
#' cp <- SplitCoupledJumpProblem(j1, j2, Direct(), list(c(1, 1), c(2, 2)))
#' tail(solve(cp, SSAStepper(), seed = 1)$u, 1)
#' @name jump_coupling
NULL

#' @rdname jump_coupling
#' @export
SplitCoupledJumpProblem <- function(prob, prob_control, aggregator, coupling_map, ...) {
  if (!inherits(prob, "JumpProblem") || !inherits(prob_control, "JumpProblem"))
    stop("SplitCoupledJumpProblem couples two JumpProblems", call. = FALSE)
  js <- prob$jumpset; jc <- prob_control$jumpset
  if (js$nmaj || jc$nmaj || length(js$vrj) || length(jc$vrj))
    stop("SplitCoupledJumpProblem couples problems built from ConstantRateJumps", call. = FALSE)
  l <- length(prob$prob$u0)
  if (length(prob_control$prob$u0) != l) stop("the coupled problems need states of the same length", call. = FALSE)
  first <- seq_len(l); second <- l + seq_len(l)
  flip <- function(integ) integ$u[c(first, second)] <- integ$u[c(second, first)]
  on_control <- function(aff) function(integ) {
    flip(integ)
    aff(integ)
    flip(integ)
  }
  coupling_map <- lapply(coupling_map, as.integer)
  r <- js$crj_rate; a <- js$crj_affect
  rc <- jc$crj_rate; ac <- jc$crj_affect
  jumps <- list()
  add <- function(rate, affect) jumps[[length(jumps) + 1L]] <<- ConstantRateJump(rate, affect)
  for (k in setdiff(seq_along(r), vapply(coupling_map, `[`, 0L, 1L)))
    local({ rr <- r[[k]]; add(function(u, p, t) rr(u[first], p, t), a[[k]]) })
  for (k in setdiff(seq_along(rc), vapply(coupling_map, `[`, 0L, 2L)))
    local({ rr <- rc[[k]]; add(function(u, p, t) rr(u[second], p, t), on_control(ac[[k]])) })
  for (cm in coupling_map) local({
    rate <- r[[cm[1L]]]; ratec <- rc[[cm[2L]]]
    aff <- a[[cm[1L]]]; affc <- ac[[cm[2L]]]
    m <- function(u, p, t) min(rate(u[first], p, t), ratec(u[second], p, t))
    add(m, function(integ) {
      aff(integ)
      on_control(affc)(integ)
    })
    add(function(u, p, t) rate(u[first], p, t) - m(u, p, t), aff)
    # the control jump alone (build_split_jumps applies the first problem's
    # affect here, an evident slip; the control's affect is meant)
    add(function(u, p, t) ratec(u[second], p, t) - m(u, p, t), on_control(affc))
  })
  inner <- .cat_problems(prob$prob, prob_control$prob)
  do.call(JumpProblem, c(list(inner, aggregator), jumps, list(...)))
}

# cat_problems: the state c(u, u_control) and the two right-hand sides
.cat_problems <- function(p1, p2) {
  u0 <- c(as.numeric(p1$u0), as.numeric(p2$u0))
  l <- length(p1$u0)
  first <- seq_len(l); second <- l + seq_len(l)
  disc <- function(p) inherits(p, "DiscreteProblem")
  fz <- function(p) if (disc(p)) function(u, pp, t) 0 * u else p$f$f
  if (disc(p1) && disc(p2)) return(DiscreteProblem(u0, p1$tspan, p1$p))
  f1 <- fz(p1); f2 <- fz(p2)
  f <- function(u, p, t) c(f1(u[first], p, t), f2(u[second], p, t))
  if (inherits(p1, "SDEProblem") || inherits(p2, "SDEProblem")) {
    g1 <- if (inherits(p1, "SDEProblem")) p1$f$g else function(u, p, t) 0 * u
    g2 <- if (inherits(p2, "SDEProblem")) p2$f$g else function(u, p, t) 0 * u
    return(SDEProblem(f, function(u, p, t) c(g1(u[first], p, t), g2(u[second], p, t)), u0, p1$tspan, p1$p))
  }
  ODEProblem(f, u0, p1$tspan, p1$p)
}

#' @rdname jump_coupling
#' @export
JumpSet <- function(...) {
  js <- unlist(lapply(list(...), function(j) {
    if (is.null(j)) list()
    else if (inherits(j, c("AbstractJump", "MassActionJumpSet"))) list(j)
    else as.list(j)
  }), recursive = FALSE)
  structure(js, class = "JumpSet")
}

#' @rdname jump_coupling
#' @export
RateBounds <- function(lrate = NULL, urate = NULL, rateinterval = Inf) {
  if (is.null(lrate) && is.null(urate)) stop("`RateBounds` requires at least one of `lrate` or `urate`.", call. = FALSE)
  structure(list(lrate = lrate %||% 0, urate = urate %||% Inf, rateinterval = rateinterval), class = "RateBounds")
}

#' @rdname jump_coupling
#' @export
JumpDependencies <- function(dep_graph = NULL, vartojumps_map = NULL, jumptovars_map = NULL)
  list(dep_graph = dep_graph, vartojumps_map = vartojumps_map, jumptovars_map = jumptovars_map)

#' @rdname jump_coupling
#' @export
get_num_majumps <- function(jumps) {
  if (is.null(jumps)) return(0L)
  if (inherits(jumps, "MassActionJumpSet")) return(length(jumps))
  if (inherits(jumps, "SpatialMassActionJump")) return(length(jumps$rs))
  if (inherits(jumps, "MassActionJump")) return(1L)
  0L
}

#' @rdname jump_coupling
#' @export
needs_depgraph <- function(agg) .needs_depgraph(agg)

#' @rdname jump_coupling
#' @export
needs_vartojumps_map <- function(agg) .needs_vartojumps_map(agg)

#' @rdname jump_coupling
#' @export
isolated_copy <- function(jprob) jprob
