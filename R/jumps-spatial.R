# Spatial jump processes, after JumpProcesses.jl src/spatial: topology.jl
# (CartesianGrid, CartesianGridRej, num_sites, neighbors, outdegree,
# nth_nbr), hop_rates.jl (the HopRates forms D_s, D_si, D_sij, D_s L_ij and
# D_si L_ij), spatial_massaction_jump.jl (SpatialMassActionJump,
# rate_at_site, evalrxrate), reaction_rates.jl (RxRates), utils.jl
# (SpatialJump, sample_jump_direct), nsm.jl (NSM), directcrdirect.jl
# (DirectCRDirect) and flatten.jl.
#
# The state is a matrix with one row per species and one column per site;
# the SSA driver works on it in column major order. For every species s and
# site i the hopping constants are kept as the vector w[[s, i]] of the
# per molecule rates to the neighbors of i (in ascending order), so the hop
# rate of s out of i is u[s, i] sum(w[[s, i]]) and the target is drawn in
# proportion to w[[s, i]]. NSM keeps a next firing time per site;
# DirectCRDirect draws the site in proportion to its total rate (the same
# law as its composition rejection table).

#' Spatial jump processes
#'
#' Reactions at the sites of a graph or a Cartesian grid and hops of the
#' species between neighboring sites, as in JumpProcesses.jl.
#' `CartesianGrid(dims)` (or `CartesianGridRej(dims)`, or
#' `CartesianGridRej(dimension, linear_size)`) is a grid with one, two or
#' three dimensions whose sites are numbered in column major order; a graph
#' is given as a list of neighbor index vectors or as a symmetric adjacency
#' matrix. `num_sites()`, `neighbors()` and `outdegree()` query either.
#'
#' `SpatialMassActionJump(uniform_rates, spatial_rates, reactant_stoch,
#' net_stoch)` holds mass action reactions whose rates are the same at all
#' sites (`uniform_rates`, a vector) or depend on the site (`spatial_rates`,
#' a matrix with one row per reaction and one column per site); uniform
#' reactions come first. The stoichiometry refers to species by index, as for
#' [MassActionJump()]; a plain `MassActionJump` has uniform rates.
#'
#' `JumpProblem(prob, NSM(), jumps, hopping_constants = , spatial_system = )`
#' with `prob = DiscreteProblem(u0, tspan, p)` and `u0` a species by sites
#' matrix builds the spatial problem; [SSAStepper()] solves it with the next
#' subvolume method `NSM()` or with `DirectCRDirect()`. The solution holds
#' the state in column major order; `matrix(sol$u[[k]], nspecies)` is the
#' species by sites matrix. `hopping_constants` is one of
#' * a vector `D` (species `s` hops to each neighbor with rate `D[s]`),
#' * a species by sites matrix `D` (rate `D[s, i]` from site `i`),
#' * a species by sites matrix of lists, each element the vector of rates to
#'   the neighbors of the site (`D[[s, i]][j]`),
#' * `list(D, L)` with `D` a vector or a species by sites matrix and `L` a
#'   list with, for each site, the vector of factors for its neighbors (rate
#'   `D[s] L[[i]][j]` or `D[s, i] L[[i]][j]`).
#'
#' `flatten_spatial(jumps, prob, spatial_system, hopping_constants)` returns
#' the equivalent non spatial problem: a `DiscreteProblem` with one unknown
#' per species and site and mass action jumps for the reactions at each site
#' and for each hop, for any aggregator (JumpProcesses' `flatten`).
#'
#' @param dims Side lengths of the grid.
#' @param linear_size With `linear_size`, `dims` is the number of dimensions of a
#'   hypercube grid.
#' @param spatial_system A grid or a graph.
#' @param site A site index.
#' @param uniform_rates,spatial_rates Reaction rate constants.
#' @param reactant_stoch,net_stoch Stoichiometry per reaction.
#' @param scale_rates Divide the rates by the factorials of the reactant
#'   stoichiometry.
#' @param useiszero Treat a single reactant entry with species 0 as no
#'   reactant.
#' @param jumps A `SpatialMassActionJump` or `MassActionJump`.
#' @param prob A `DiscreteProblem` with a species by sites initial matrix.
#' @param hopping_constants Hopping constants, see Details.
#' @return A grid, a jump, an aggregator, a number or neighbor vector, or (for
#'   `flatten_spatial()`) a list with `prob` and `jumps`.
#' @examples
#' grid <- CartesianGrid(c(3, 3))
#' neighbors(grid, 1)
#' u0 <- matrix(0, 1, 9); u0[1, 5] <- 100
#' jp <- JumpProblem(DiscreteProblem(u0, c(0, 1)), NSM(),
#'                   SpatialMassActionJump(0.1, NULL, list(list(c(1, 1))), list(list(c(1, -1)))),
#'                   hopping_constants = 1, spatial_system = grid)
#' sol <- solve(jp, SSAStepper(), seed = 1)
#' matrix(tail(sol$u, 1)[[1]], 1)
#' @name spatial_jumps
NULL

#' @rdname spatial_jumps
#' @export
CartesianGridRej <- function(dims, linear_size = NULL) {
  if (!is.null(linear_size)) dims <- rep(as.integer(linear_size), as.integer(dims))
  dims <- as.integer(dims)
  if (!length(dims) || length(dims) > 3L) stop("Cartesian grids have one, two or three dimensions", call. = FALSE)
  n <- prod(dims)
  # offsets in the order of potential_offsets (ascending linear index)
  strides <- cumprod(c(1L, dims))[seq_along(dims)]
  nbrs <- lapply(seq_len(n), function(site) {
    ci <- as.integer(arrayInd(site, dims))
    out <- integer()
    for (d in rev(seq_along(dims))) if (ci[d] > 1L) out <- c(out, site - strides[d])
    for (d in seq_along(dims)) if (ci[d] < dims[d]) out <- c(out, site + strides[d])
    sort(out)
  })
  structure(list(dims = dims, nbrs = nbrs), class = c("CartesianGridRej", "SpatialSystem"))
}

#' @rdname spatial_jumps
#' @export
CartesianGrid <- function(dims) CartesianGridRej(dims)

#' @export
print.CartesianGridRej <- function(x, ...) {
  cat("A Cartesian grid with dimensions (", paste(x$dims, collapse = ", "), ")\n", sep = "")
  invisible(x)
}

# The neighbor lists of a spatial system: a grid, a list of neighbor
# vectors or a symmetric adjacency matrix.
.sp_nbrs <- function(spatial_system) {
  if (inherits(spatial_system, "CartesianGridRej")) return(spatial_system$nbrs)
  if (is.matrix(spatial_system)) {
    if (!isSymmetric(unname(spatial_system))) stop("an adjacency matrix must be symmetric", call. = FALSE)
    return(lapply(seq_len(nrow(spatial_system)), function(i) which(spatial_system[i, ] != 0 & seq_len(nrow(spatial_system)) != i)))
  }
  if (is.list(spatial_system)) return(lapply(spatial_system, function(v) sort(as.integer(v))))
  stop("a spatial system is a CartesianGrid, a list of neighbor vectors or an adjacency matrix", call. = FALSE)
}

#' @rdname spatial_jumps
#' @export
num_sites <- function(spatial_system) length(.sp_nbrs(spatial_system))

#' @rdname spatial_jumps
#' @export
neighbors <- function(spatial_system, site) .sp_nbrs(spatial_system)[[site]]

#' @rdname spatial_jumps
#' @export
outdegree <- function(spatial_system, site) length(.sp_nbrs(spatial_system)[[site]])

#' @rdname spatial_jumps
#' @export
SpatialMassActionJump <- function(uniform_rates = NULL, spatial_rates = NULL, reactant_stoch, net_stoch,
                                  scale_rates = TRUE, useiszero = TRUE) {
  if (inherits(uniform_rates, "MassActionJumpSet") || inherits(uniform_rates, "MassActionJump")) {
    js <- if (inherits(uniform_rates, "MassActionJump")) list(uniform_rates) else unclass(uniform_rates)
    rs <- lapply(js, function(j) lapply(j$reactant_stoch, function(s) c(s$var, s$stoich)))
    ns <- lapply(js, function(j) lapply(j$net_stoch, function(s) c(s$var, s$stoich)))
    return(SpatialMassActionJump(vapply(js, function(j) as.numeric(j$scaled_rates), 0), NULL, rs, ns,
                                 scale_rates = FALSE, useiszero = useiszero))
  }
  if (is.matrix(uniform_rates) && is.null(spatial_rates)) {
    spatial_rates <- uniform_rates
    uniform_rates <- NULL
  }
  nu <- length(uniform_rates)
  ns_ <- if (is.null(spatial_rates)) 0L else nrow(spatial_rates)
  if (length(reactant_stoch) != nu + ns_ || length(net_stoch) != nu + ns_)
    stop("give one reactant and one net stoichiometry per reaction (uniform reactions first)", call. = FALSE)
  as_mat <- function(st) {
    if (!length(st)) return(matrix(integer(), 0L, 2L))
    m <- do.call(rbind, lapply(st, function(v) as.numeric(v)))
    if (useiszero && nrow(m) == 1L && m[1L, 1L] == 0) return(matrix(integer(), 0L, 2L))
    m
  }
  rs <- lapply(reactant_stoch, as_mat)
  ns <- lapply(net_stoch, as_mat)
  if (scale_rates) {
    fac <- vapply(rs, function(m) prod(factorial(m[, 2L])), 0)
    if (nu) uniform_rates <- uniform_rates / fac[seq_len(nu)]
    if (ns_) spatial_rates <- spatial_rates / fac[nu + seq_len(ns_)]
  }
  structure(list(uniform_rates = uniform_rates, spatial_rates = spatial_rates, rs = rs, ns = ns),
            class = c("SpatialMassActionJump", "AbstractJump"))
}

#' @rdname spatial_jumps
#' @export
NSM <- function() .new_aggregator("NSM")

#' @rdname spatial_jumps
#' @export
DirectCRDirect <- function() .new_aggregator("DirectCRDirect")

.is_spatial_agg <- function(agg) agg$name %in% c("NSM", "DirectCRDirect")

# rate_at_site
.sp_rate_const <- function(smaj, rx, site) {
  nu <- length(smaj$uniform_rates)
  if (rx <= nu) smaj$uniform_rates[rx] else smaj$spatial_rates[rx - nu, site]
}

# The per neighbor hopping constants w[[s, i]] (HopRates).
.sp_hop_weights <- function(hc, nbrs, nspec) {
  nsites <- length(nbrs)
  deg <- lengths(nbrs)
  w <- matrix(list(), nspec, nsites)
  site_factor <- NULL
  if (is.list(hc) && !is.matrix(hc) && length(hc) == 2L && is.list(hc[[2L]])) {
    site_factor <- hc[[2L]]
    hc <- hc[[1L]]
    if (length(site_factor) != nsites) stop("give one vector of site hopping constants per site", call. = FALSE)
  }
  for (s in seq_len(nspec)) for (i in seq_len(nsites)) {
    v <- if (is.matrix(hc) && is.list(hc)) as.numeric(hc[[s, i]])
         else if (is.matrix(hc)) rep(hc[s, i], deg[i])
         else rep(as.numeric(hc)[if (length(hc) == 1L) 1L else s], deg[i])
    if (!is.null(site_factor)) v <- v * as.numeric(site_factor[[i]])
    if (length(v) != deg[i])
      .stopf("site %d has %d neighbors but %d hopping constants were given for species %d", i, deg[i], length(v), s)
    w[[s, i]] <- v
  }
  w
}

.sp_check_problem <- function(prob, jumps, hopping_constants, spatial_system) {
  if (!inherits(prob, "DiscreteProblem")) stop("spatial jump problems need a DiscreteProblem", call. = FALSE)
  if (is.null(spatial_system) || is.null(hopping_constants))
    stop("spatial jump problems need `hopping_constants` and `spatial_system`", call. = FALSE)
  u0 <- prob$u0
  nbrs <- .sp_nbrs(spatial_system)
  nspec <- if (is.matrix(u0)) nrow(u0) else 1L
  if (length(u0) != nspec * length(nbrs))
    stop("u0 must be a species by sites matrix with one column per site", call. = FALSE)
  smaj <- if (!length(jumps)) SpatialMassActionJump(NULL, NULL, list(), list())
          else if (inherits(jumps[[1L]], "SpatialMassActionJump")) jumps[[1L]]
          else SpatialMassActionJump(structure(jumps, class = "MassActionJumpSet"))
  if (length(Filter(function(j) !inherits(j, c("SpatialMassActionJump", "MassActionJump")), jumps)))
    stop("spatial jump problems take mass action reactions only", call. = FALSE)
  list(nbrs = nbrs, nspec = nspec, nsites = length(nbrs), smaj = smaj,
       w = .sp_hop_weights(hopping_constants, nbrs, nspec))
}

# The spatial aggregator: rates of reactions and hops at each site.
.sp_agg_new <- function(jp) {
  a <- new.env(parent = emptyenv())
  sp <- jp$spatial
  a$type <- jp$aggregator$name
  a$spatial <- sp
  a$save_positions <- jp$save_positions
  a$nrx <- length(sp$smaj$rs)
  a$next_jump_time <- Inf
  a$end_time <- Inf
  a
}

.sp_u <- function(a, integ) matrix(integ$u[seq_len(a$spatial$nspec * a$spatial$nsites)], a$spatial$nspec)

# rates of all reactions and hops at one site
.sp_site_rates <- function(a, u, site) {
  sp <- a$spatial
  rx <- vapply(seq_len(a$nrx), function(r) .evalrxrate(u[, site], sp$smaj$rs[[r]], .sp_rate_const(sp$smaj, r, site)), 0)
  hop <- vapply(seq_len(sp$nspec), function(s) u[s, site] * sum(sp$w[[s, site]]), 0)
  a$rx[, site] <- rx
  a$hop[, site] <- hop
  a$site_rate[site] <- sum(rx) + sum(hop)
}

.sp_new_time <- function(a, site, t) {
  r <- a$site_rate[site]
  a$pq[site] <- if (r > 0) t + stats::rexp(1L) / r else Inf
}

.sp_initialize <- function(a, integ) {
  sp <- a$spatial
  a$end_time <- integ$prob$tspan[2L]
  u <- .sp_u(a, integ)
  a$rx <- matrix(0, a$nrx, sp$nsites)
  a$hop <- matrix(0, sp$nspec, sp$nsites)
  a$site_rate <- numeric(sp$nsites)
  for (i in seq_len(sp$nsites)) .sp_site_rates(a, u, i)
  if (a$type == "NSM") {
    a$pq <- numeric(sp$nsites)
    for (i in seq_len(sp$nsites)) .sp_new_time(a, i, integ$t)
  }
  .sp_generate(a, integ)
}

# sample_jump_direct at a site: a reaction or a hop and its target
.sp_sample_at_site <- function(a, site) {
  sp <- a$spatial
  rxs <- sum(a$rx[, site])
  if (stats::runif(1L) * a$site_rate[site] < rxs) {
    r <- .linear_search(a$rx[, site], stats::runif(1L) * rxs)
    return(list(src = site, rx = max(r, 1L), dst = site))
  }
  hs <- a$hop[, site]
  s <- max(.linear_search(hs, stats::runif(1L) * sum(hs)), 1L)
  wv <- sp$w[[s, site]]
  k <- max(.linear_search(wv, stats::runif(1L) * sum(wv)), 1L)
  list(src = site, species = s, dst = sp$nbrs[[site]][k])
}

.sp_generate <- function(a, integ) {
  t <- integ$t
  if (a$type == "NSM") {
    site <- which.min(a$pq)
    a$next_jump_time <- if (length(site)) a$pq[site] else Inf
  } else {
    total <- sum(a$site_rate)
    a$next_jump_time <- if (total > 0) t + stats::rexp(1L) / total else Inf
    site <- if (total > 0) max(.linear_search(a$site_rate, stats::runif(1L) * total), 1L) else integer()
  }
  if (a$next_jump_time >= a$end_time || !length(site)) return(invisible())
  a$next_jump <- .sp_sample_at_site(a, site)
  invisible()
}

.sp_execute <- function(a, integ) {
  sp <- a$spatial
  j <- a$next_jump
  n <- sp$nspec
  idx <- function(s, i) (i - 1L) * n + s
  if (!is.null(j$rx)) {
    ns <- sp$smaj$ns[[j$rx]]
    for (r in seq_len(nrow(ns))) {
      k <- idx(ns[r, 1L], j$src)
      integ$u[k] <- integ$u[k] + ns[r, 2L]
    }
    sites <- j$src
  } else {
    integ$u[idx(j$species, j$src)] <- integ$u[idx(j$species, j$src)] - 1
    integ$u[idx(j$species, j$dst)] <- integ$u[idx(j$species, j$dst)] + 1
    sites <- c(j$src, j$dst)
  }
  u <- .sp_u(a, integ)
  for (i in sites) {
    .sp_site_rates(a, u, i)
    if (a$type == "NSM") .sp_new_time(a, i, integ$t)
  }
  invisible()
}

#' @rdname spatial_jumps
#' @export
flatten_spatial <- function(jumps, prob, spatial_system, hopping_constants) {
  if (inherits(jumps, "AbstractJump")) jumps <- list(jumps)
  sp <- .sp_check_problem(prob, jumps, hopping_constants, spatial_system)
  n <- sp$nspec
  idx <- function(s, i) (i - 1L) * n + s
  rates <- numeric(); rs <- list(); ns <- list()
  add <- function(k, r, nn) {
    rates[length(rates) + 1L] <<- k
    rs[[length(rs) + 1L]] <<- r
    ns[[length(ns) + 1L]] <<- nn
  }
  # hops first, then the reactions at each site (flatten.jl)
  for (src in seq_len(sp$nsites)) for (s in seq_len(n)) {
    for (k in seq_along(sp$nbrs[[src]])) {
      dst <- sp$nbrs[[src]][k]
      add(sp$w[[s, src]][k], list(c(idx(s, src), 1)), list(c(idx(s, src), -1), c(idx(s, dst), 1)))
    }
  }
  for (site in seq_len(sp$nsites)) for (r in seq_along(sp$smaj$rs)) {
    m <- sp$smaj$rs[[r]]; nm <- sp$smaj$ns[[r]]
    add(.sp_rate_const(sp$smaj, r, site),
        lapply(seq_len(nrow(m)), function(q) c(idx(m[q, 1L], site), m[q, 2L])),
        lapply(seq_len(nrow(nm)), function(q) c(idx(nm[q, 1L], site), nm[q, 2L])))
  }
  fprob <- DiscreteProblem(as.numeric(prob$u0), prob$tspan, prob$p)
  jumps <- if (length(rates) > 1L) MassActionJump(rates, rs, ns, scale_rates = FALSE)
           else if (length(rates)) MassActionJump(rates, rs[[1L]], ns[[1L]], scale_rates = FALSE)
  list(prob = fprob, jumps = jumps)
}
