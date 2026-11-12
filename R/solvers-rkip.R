# Runge-Kutta in the interaction picture, after OrdinaryDiffEqRKIP.jl 2.1
# (algorithms.jl, rkip_cache.jl, rkip_utils.jl, rkip_perform_step.jl), for
# semilinear problems u' = A u + f(u, t) given as SplitODEProblem(L, f, ...)
# with a constant SciML operator L. The exponentials exp(+-A c_i dt) of the
# stages are cached for a geometric grid of step sizes, and the proposed
# steps are moved onto that grid. The cache lives in the algorithm object
# (an environment), so solving again with the same RKIP() reuses it.

# The Verner6 tableau of OrdinaryDiffEqExplicitTableaus (Verner's "most
# efficient" 6(5) pair), as the Float64 values used upstream.
.rkip_verner6 <- function() {
  A <- matrix(0, 8L, 8L)
  A[2L, 1L] <- 0.14285714285714285
  A[3L, 1:2] <- c(0.04938271604938271, 0.1728395061728395)
  A[4L, 1:3] <- c(0.21209912536443148, -0.5510204081632653, 0.7674927113702624)
  A[5L, 1:4] <- c(0.2895622895622896, -0.42424242424242425, 0.2937062937062937, 0.5076405076405076)
  A[6L, 1:5] <- c(-0.011851917613636364, 0.2684659090909091, 0.10772747760052448, 0.1602190777972028,
                  0.225439453125)
  A[7L, 1:6] <- c(1.584056712962963, -2.15625, -1.3347480028195489, 4.356905864197531,
                  -3.8984933035714286, 2.4485287292304836)
  A[8L, 1:7] <- c(0.07896675084175084, -0.23863636363636365, 0.6312080189547294, 0.17677523094189762,
                  -0.6515625, 1.0032488628979856, 0)
  ExplicitRKTableau(A, c(0, 0.14285714285714285, 0.2222222222222222, 0.42857142857142855,
                         0.6666666666666666, 0.75, 1, 1),
                    c(0.07314814814814814, 0, 0.28460092539039905, 0.19951329534662868, 0,
                      0.3659519168291098, 0.01, 0.06678571428571428),
                    6, c(0.07064814814814815, 0, 0.298830971659919, 0.1687312440645774, 0.06075,
                         0.3210396361273554, 0, 0.08), adaptiveorder = 5)
}

# Base.logrange(a, b, n): a geometric progression with exact endpoints
.rkip_logrange <- function(a, b, n) {
  n <- as.integer(n)
  if (n == 1L) return(a)
  i <- seq_len(n)
  x <- exp((n - i) * (log(a) / (n - 1L)) + (i - 1L) * (log(b) / (n - 1L)))
  x[1L] <- a
  x[n] <- b
  x
}

# exp(L, h): the exponential of h times the operator (a matrix or a number)
.rkip_exp <- function(L, h) {
  if (inherits(L, "ScalarOperator")) return(exp(h * L$value))
  if (inherits(L, "DiagonalOperator")) return(diag(exp(h * L$diag), length(L$diag)))
  .eu_exp_base(h * .sop_matrix(L))
}

#' Runge-Kutta in the interaction picture
#'
#' `RKIP()` solves a semilinear problem `u' = A u + f(u, t)`, given as
#' `SplitODEProblem(L, f, u0, tspan)` with a constant SciML operator `L`
#' (such as `MatrixOperator(A)` or `DiagonalOperator(d)`), with an explicit
#' Runge-Kutta method (by default Verner's 6(5) pair) in the interaction
#' picture `v = exp(-A t) u`, as in OrdinaryDiffEqRKIP. The exponentials
#' `exp(A c_i dt)` are cached for `nb_of_cache_step` step sizes in geometric
#' progression from `dtmin` to `dtmax`; the proposed steps are moved onto
#' these sizes (and clamped to them when `clamp_lower_dt` or
#' `clamp_higher_dt`). The cache is kept in the algorithm object, so solving
#' several problems with the same linear part and the same `RKIP()` object
#' reuses it.
#'
#' @param dtmin,dtmax The smallest and largest cached step sizes.
#' @param nb_of_cache_step The number of cached step sizes.
#' @param tableau An [ExplicitRKTableau()] with an embedded error estimate;
#'   by default Verner's 6(5) tableau.
#' @param clamp_lower_dt,clamp_higher_dt Whether to clamp the proposed step
#'   to the smallest or largest cached step.
#' @param use_ldiv Whether to apply `exp(-A h)` by solving with `exp(A h)`
#'   instead of caching it.
#' @return An algorithm object.
#' @examples
#' L <- MatrixOperator(matrix(c(-2, 1, 0, -3), 2))
#' prob <- SplitODEProblem(L, function(u, p, t) -u * u * u, c(1, 0.5), c(0, 1))
#' sol <- solve(prob, RKIP())
#' sol$u[[length(sol$u)]]
#' @export
RKIP <- function(dtmin = 1e-3, dtmax = 1, nb_of_cache_step = 100, tableau = NULL,
                 clamp_lower_dt = FALSE, clamp_higher_dt = TRUE, use_ldiv = FALSE) {
  tab <- tableau %||% .rkip_verner6()
  if (is.null(tab$alphaEEst)) stop("RKIP() needs a tableau with an embedded error estimate", call. = FALSE)
  grid <- .rkip_logrange(dtmin, dtmax, nb_of_cache_step)
  ng <- length(grid)
  s <- tab$stages
  cc <- tab$c
  c_unique <- unique(c(cc, 1))
  c_mapping <- match(c(cc, 1), c_unique)
  # the cache of the algorithm (RKIPCache): kept across solves
  st <- new.env(parent = emptyenv())
  st$last_step <- 0
  st$cached <- logical(ng)
  st$pos <- st$neg <- vector("list", length(c_unique))
  st$cache_pos <- st$cache_neg <- vector("list", ng)

  # cache_exp_op_for_this_step!
  prepare <- function(L, dt) {
    if (dt != 0 && !.eu_isapprox(dt, st$last_step)) {
      dt_abs <- abs(dt)
      action <- "single_use"
      step_index <- min(max(findInterval(dt_abs, grid), 1L), ng)
      if (.eu_isapprox(grid[step_index], dt_abs))
        action <- if (st$cached[step_index]) "use_cached" else "cache"
      if (action == "cache") st$cache_pos[[step_index]] <- st$cache_neg[[step_index]] <- vector("list", length(c_unique))
      for (k in seq_along(c_unique)) {
        h <- abs(dt * c_unique[k])
        if (action == "use_cached") {
          st$pos[[k]] <- st$cache_pos[[step_index]][[k]]
          st$neg[[k]] <- st$cache_neg[[step_index]][[k]]
        } else {
          st$pos[[k]] <- .rkip_exp(L, h)
          if (!use_ldiv) st$neg[[k]] <- .rkip_exp(L, -h)
          if (action == "cache") {
            st$cache_pos[[step_index]][[k]] <- st$pos[[k]]
            st$cache_neg[[step_index]][[k]] <- st$neg[[k]]
          }
        }
      }
      st$cached[step_index] <- st$cached[step_index] || action == "cache"
    }
    st$last_step <- dt
  }
  # expmv_rkip_mip: exp(+-A c_i dt) v for stage i (the last index is c = 1)
  expmv <- function(v, positive, i) {
    k <- c_mapping[i]
    if (positive) return(.exprk_mv(st$pos[[k]], v))
    if (use_ldiv) return(if (is.matrix(st$pos[[k]])) as.numeric(solve(st$pos[[k]], v)) else v / st$pos[[k]])
    .exprk_mv(st$neg[[k]], v)
  }
  operator <- function(ctx) {
    if (!inherits(ctx$fobj, "SplitFunction") || !inherits(ctx$fobj$f1, "AbstractSciMLOperator"))
      stop("RKIP() is only implemented for semilinear problems: a SplitODEProblem whose ",
           "first function is a SciML operator", call. = FALSE)
    ctx$fobj$f1
  }
  init <- function(ctx, t, u) {
    L <- operator(ctx)
    f2 <- .as_rhs(ctx$fobj$f2)
    u <- as.numeric(u)
    list(fsal = .eu_mul(L, u) + f2(u, ctx$p, t), nf20 = 1L)
  }
  A <- tab$A
  al <- tab$alpha
  ae <- tab$alphaEEst
  step <- function(ctx, t, u, dt, cache) {
    L <- operator(ctx)
    f2 <- .as_rhs(ctx$fobj$f2)
    p <- ctx$p
    uprev <- as.numeric(u)
    prepare(L, dt)
    unew <- uprev
    kk <- vector("list", s)
    utilde <- NULL
    for (i in seq_len(s)) {
      ki <- uprev
      for (j in seq_len(i - 1L)) ki <- ki + (A[i, j] * dt) * kk[[j]]
      t_ <- t + cc[i] * dt
      # to the true coordinates, the nonlinear part, back to the interaction picture
      ki <- expmv(ki, TRUE, i)
      ki <- f2(ki, p, t_)
      ki <- expmv(ki, FALSE, i)
      kk[[i]] <- ki
      unew <- unew + (al[i] * dt) * ki
      if (ctx$adaptive)
        utilde <- if (i == 1L) dt * (al[i] - ae[i]) * ki else utilde + (dt * (al[i] - ae[i])) * ki
    }
    unew <- expmv(unew, TRUE, s + 1L)
    fl <- .eu_mul(L, unew) + f2(unew, p, t + dt)
    res <- list(u = unew, k = list(cache$fsal, fl), cache = list(fsal = fl),
                nf = 2L * s + 1L, nf2 = s + 1L)
    if (ctx$adaptive) res$err <- expmv(utilde, TRUE, s + 1L)
    res
  }
  # dtnew_modification: move the proposed step onto the grid
  dtnew_mod <- function(dtnew, eigen_est = NULL) {
    if (grid[1L] > dtnew && clamp_lower_dt) return(grid[1L])
    if (grid[ng] < dtnew && clamp_higher_dt) return(grid[ng])
    i <- sum(grid < dtnew) + 1L
    # upstream indexes past the grid (a BoundsError) for a step above dtmax
    # when clamp_higher_dt = FALSE; the step is then kept
    if (i > ng) dtnew else grid[i]
  }
  .new_alg("RKIP", tab$order, adaptive = TRUE, step = step, init = init,
           adaptive_order = tab$adaptiveorder, dtnew_mod = dtnew_mod, tableau = tab,
           dt_for_exp_caching = grid, clamp_lower_dt = clamp_lower_dt,
           clamp_higher_dt = clamp_higher_dt, use_ldiv = use_ldiv, cache = st)
}
