# Rosenbrock methods with mass matrices, after OrdinaryDiffEq.jl:
#   lib/OrdinaryDiffEqRosenbrock/src/rosenbrock_perform_step.jl
#     perform_step! for Rosenbrock23ConstantCache and the unified Rodas loop
#     (RosenbrockCombinedConstantCache)
#   lib/OrdinaryDiffEqRosenbrock/src/rosenbrock_interpolants.jl
#     Shampine's interpolant of Rosenbrock23 and the H matrix dense output
#   lib/OrdinaryDiffEqRosenbrock/src/rosenbrock_tableaus.jl
#     Rosenbrock23Tableau, RODAS5PA/C/c/d/H, Rodas5PTableau
#   lib/OrdinaryDiffEqRosenbrockTableaus/src/rosenbrock_tableaus.jl
#     RODAS4A/C/c/d/H, Rodas4Tableau
#   lib/OrdinaryDiffEqDifferentiation/src/derivative_utils.jl
#     calc_rosenbrock_differentiation (W = J - M / (dt gamma)), do_newJW
#
# The Jacobian is prob$f$jac when present, else forward differences; the
# time derivative is prob$f$tgrad when present, else a forward difference.
# Both are reused when the previous attempt failed the error test (the
# non Newton branch of do_newJW). Singular mass matrices (index 1 DAEs) are
# supported: W = J - M / (dt gamma) stays invertible.

# Rosenbrock23Tableau
.ROS23 <- list(c32 = 6 + sqrt(2), d = 1 / (2 + sqrt(2)))

# RodasTableau: A, C, gamma, c, d, H; b is the last row of A followed by 1
# and btilde the last unit vector (stiffly accurate embedded method).
.rodas_tab <- function(A, C, gamma, c, d, H) {
  s <- nrow(A)
  list(A = A, C = C, gamma = gamma, c = c, d = d, H = H,
       b = c(A[s, seq_len(s - 1L)], 1), btilde = c(numeric(s - 1L), 1))
}

.RODAS4 <- .rodas_tab(
  A = rbind(
    c(0, 0, 0, 0, 0, 0),
    c(1.544, 0, 0, 0, 0, 0),
    c(0.9466785280815826, 0.2557011698983284, 0, 0, 0, 0),
    c(3.314825187068521, 2.896124015972201, 0x1.ff4dfe3212151p-1, 0, 0, 0),
    c(1.221224509226641, 6.019134481288629, 12.53708332932087, -0.687886036105895, 0, 0),
    c(1.221224509226641, 6.019134481288629, 12.53708332932087, -0.687886036105895, 1, 0)),
  C = rbind(
    c(0, 0, 0, 0, 0),
    c(-5.6688, 0, 0, 0, 0),
    c(-2.430093356833875, -0.2063599157091915, 0, 0, 0),
    c(-0.1073529058151375, -0x1.3306a769f29b4p+3, -20.47028614809616, 0, 0),
    c(7.496443313967647, -10.24680431464352, -33.99990352819905, 11.7089089320616, 0),
    c(8.083246795921522, -7.981132988064893, -31.52159432874371, 16.31930543123136, -6.058818238834054)),
  gamma = 0.25,
  c = c(0, 0.386, 0.21, 0.63, 1, 1),
  d = c(0.25, -0.1043, 0.1035, -0.0362, 0, 0),
  H = rbind(
    c(10.12623508344586, -7.487995877610167, -34.80091861555747, -7.992771707568823, 1.025137723295662, 0),
    c(-0.6762803392801253, 6.087714651680015, 16.43084320892478, 24.76722511418386, -6.594389125716872, 0)))

.RODAS5P <- .rodas_tab(
  A = rbind(
    c(0, 0, 0, 0, 0, 0, 0, 0),
    c(3.0, 0, 0, 0, 0, 0, 0, 0),
    c(2.849394379747939, 0x1.d56caff9bd456p-2, 0, 0, 0, 0, 0, 0),
    c(-6.954028509809101, 2.489845061869568, -10.358996098473584, 0, 0, 0, 0, 0),
    c(2.8029986275628964, 0.5072464736228206, -0.3988312541770524, -0.04721187230404641, 0, 0, 0, 0),
    c(-7.502846399306121, 2.561846144803919, -11.627539656261098, -0.18268767659942256, 0x1.eec4504b32a44p-6, 0, 0, 0),
    c(-7.502846399306121, 2.561846144803919, -11.627539656261098, -0.18268767659942256, 0x1.eec4504b32a44p-6, 1, 0, 0),
    c(-7.502846399306121, 2.561846144803919, -11.627539656261098, -0.18268767659942256, 0x1.eec4504b32a44p-6, 1, 1, 0)),
  C = rbind(
    c(0, 0, 0, 0, 0, 0, 0),
    c(-0x1.c4f6adfeb3d76p+3, 0, 0, 0, 0, 0, 0),
    c(-17.97296035885952, -2.859693295451294, 0, 0, 0, 0, 0),
    c(147.12150275711716, -1.41221402718213, 71.68940251302358, 0, 0, 0, 0),
    c(165.43517024871676, -0.4592823456491126, 42.90938336958603, -5.961986721573306, 0, 0, 0),
    c(24.854864614690072, -3.0009227002832186, 47.4931110020768, 5.5814197821558125, -0.6610691825249471, 0, 0),
    c(30.91273214028599, -3.1208243349937974, 77.79954646070892, 34.28646028294783, -19.097331116725623, -0x1.c1683717071e4p+4, 0),
    c(37.80277123390563, -3.2571969029072276, 112.26918849496327, 66.9347231244047, -40.06618937091002, -54.66780262877968, -9.48861652309627)),
  gamma = 0x1.b20c5235b5100p-3,
  c = c(0, 0.6358126895828704, 0.4095798393397535, 0.9769306725060716, 0.4288403609558664, 1, 1, 1),
  d = c(0x1.b20c5235b5100p-3, -0x1.b20c5235b5100p-2, -0.3384627126235924, 0x1.cdfd3bc9c89dep+0,
        2.325825639765069, 0, 0, 0),
  H = rbind(
    c(0x1.9f2e3b2083ef4p+4, -2.5579724845846235, 0x1.4de1d0d80ff5ap+3, -2.3679251022685204, 0.524948541321073,
      1.1241088310450404, 0.4272876194431874, -0x1.604d2e4930bdep-3),
    c(-9.91568850695171, -0.9689944594115154, 3.0438037242978453, -24.495224566215796, 20.176138334709044,
      15.98066361424651, -6.789040303419874, -6.710236069923372),
    c(11.419903575922262, 2.8879645146136994, 72.92137995996029, 80.12511834622643, -52.072871366152654,
      -59.78993625266729, -0.15582684282751913, 4.883087185713722)))

.ros_init <- function(ctx, t, u) {
  cache <- new.env(parent = emptyenv())
  cache$fsal <- ctx$f(u, ctx$p, t)
  cache$J <- NULL
  cache$dT <- NULL
  cache
}

# J and dT for this step: fresh unless the previous attempt failed the
# error test (do_newJW, non Newton branch: J is kept across rejections, but
# not for composite algorithms).
.ros_derivatives <- function(ctx, t, u, fu, cache) {
  reuse <- !is.null(cache$J) && !isTRUE(ctx$composite) && ctx$iter > 1L && ctx$EEst_prev > 1 &&
    isTRUE(cache$J_t == t)
  cnt <- list(nf = 0L, njacs = 0L)
  if (!reuse) {
    jj <- .calc_J(ctx, u, t, fu)
    tg <- .calc_tgrad(ctx, u, t, fu)
    cache$J <- jj$J
    cache$dT <- tg$dT
    cache$J_t <- t
    cnt <- list(nf = jj$nf + tg$nf, njacs = 1L)
  }
  cnt$eigen_est <- if (isTRUE(ctx$composite)) max(rowSums(abs(cache$J))) else NULL
  cnt
}

.ros_W <- function(ctx, J, dtgamma) {
  .factorize(if (is.null(ctx$mass)) J - diag(1 / dtgamma, nrow(J)) else J - ctx$mass / dtgamma)
}

.rosenbrock23_step <- function(ctx, t, u, dt, cache) {
  f <- ctx$f; p <- ctx$p; M <- ctx$mass
  c32 <- .ROS23$c32; d <- .ROS23$d
  dtgamma <- dt * d
  neginvdtgamma <- -1 / dtgamma
  fsalfirst <- cache$fsal
  der <- .ros_derivatives(ctx, t, u, fsalfirst, cache)
  W <- .ros_W(ctx, cache$J, dtgamma)
  out <- list(nf = der$nf, njacs = der$njacs, nw = 1L, nsolve = 0L, eigen_est = der$eigen_est, cache = cache)
  # a singular W rejects the step (EEst = 2 in Julia) and fails a fixed step
  if (is.null(W)) return(c(out, list(u = u, fail = TRUE)))
  dT <- cache$dT
  k1 <- .wsolve(W, fsalfirst + dtgamma * dT) * neginvdtgamma
  f1 <- f(u + dt / 2 * k1, p, t + dt / 2)
  k2 <- .wsolve(W, f1 - .mass_mul(M, k1)) * neginvdtgamma + k1
  unew <- u + dt * k2
  out$nf <- out$nf + 1L
  out$nsolve <- 2L
  EEst <- 0
  fsallast <- NULL
  if (ctx$adaptive) {
    fsallast <- f(unew, p, t + dt)
    out$nf <- out$nf + 1L
    lt <- if (is.null(M)) fsallast - c32 * (k2 - f1) - 2 * (k1 - fsalfirst) + dt * dT else
      fsallast - .mass_mul(M, c32 * k2 + 2 * k1) + c32 * f1 + 2 * fsalfirst + dt * dT
    k3 <- .wsolve(W, lt) * neginvdtgamma
    out$nsolve <- 3L
    utilde <- .mass_mul(M, dt / 6 * (k1 - 2 * k2 + k3))
    EEst <- .err_norm(utilde, u, unew, ctx$abstol, ctx$reltol)
    if (!is.null(M)) EEst <- EEst + .alg_penalty(ctx, fsallast)
  }
  if (isTRUE(EEst <= 1)) {
    if (is.null(fsallast)) {
      fsallast <- f(unew, p, t + dt)
      out$nf <- out$nf + 1L
    }
    cache$fsal <- fsallast
  }
  c(out, list(u = unew, EEst = EEst, k = list(k1, k2)))
}

# Shampine's interpolant (MATLAB ode23s), rosenbrock2332pre0.
.rosenbrock23_interp <- function(theta, dt, y0, y1, k) {
  d <- .ROS23$d
  c1 <- theta * (1 - theta) / (1 - 2 * d)
  c2 <- theta * (theta - 2 * d) / (1 - 2 * d)
  y0 + dt * (c1 * k[[1L]] + c2 * k[[2L]])
}

# The unified Rodas loop of rosenbrock_perform_step.jl.
.rodas_step <- function(tab) {
  force(tab)
  function(ctx, t, u, dt, cache) {
    f <- ctx$f; p <- ctx$p; M <- ctx$mass
    A <- tab$A; C <- tab$C; cc <- tab$c; d <- tab$d
    s <- nrow(A)
    invdt <- 1 / dt
    dtgamma <- dt * tab$gamma
    du <- cache$fsal
    der <- .ros_derivatives(ctx, t, u, du, cache)
    W <- .ros_W(ctx, cache$J, dtgamma)
    out <- list(nf = der$nf, njacs = der$njacs, nw = 1L, nsolve = 0L, eigen_est = der$eigen_est,
                cache = cache)
    if (is.null(W)) return(c(out, list(u = u, fail = TRUE)))
    dT <- cache$dT
    ks <- matrix(0, length(u), s)
    ks[, 1L] <- .wsolve(W, -(du + (dt * d[1L]) * dT))
    for (stage in 2:s) {
      ui <- u + as.numeric(ks[, seq_len(stage - 1L), drop = FALSE] %*% A[stage, seq_len(stage - 1L)])
      # stages that repeat the previous abscissa and row reuse its f value
      if (!(stage > 2L && cc[stage] == cc[stage - 1L] &&
            all(A[stage, seq_len(stage - 1L)] == A[stage - 1L, seq_len(stage - 1L)]))) {
        du <- f(ui, p, t + cc[stage] * dt)
        out$nf <- out$nf + 1L
      }
      lt <- as.numeric(ks[, seq_len(stage - 1L), drop = FALSE] %*% (C[stage, seq_len(stage - 1L)] * invdt))
      lt <- du + (dt * d[stage]) * dT + .mass_mul(M, lt)
      ks[, stage] <- .wsolve(W, -lt)
    }
    out$nsolve <- s
    unew <- u + as.numeric(ks %*% tab$b)
    EEst <- if (ctx$adaptive) .err_norm(as.numeric(ks %*% tab$btilde), u, unew, ctx$abstol, ctx$reltol) else 0
    kd <- lapply(seq_len(nrow(tab$H)), function(j) as.numeric(ks %*% tab$H[j, ]))
    if (isTRUE(EEst <= 1)) {
      cache$fsal <- f(unew, p, t + dt)
      out$nf <- out$nf + 1L
    }
    c(out, list(u = unew, EEst = EEst, k = kd))
  }
}

# The H matrix dense output (interp_order 2 for two rows, 3 for three).
.rodas_interp <- function(theta, dt, y0, y1, k) {
  t1 <- 1 - theta
  if (length(k) == 2L) return(t1 * y0 + theta * (y1 + t1 * (k[[1L]] + theta * k[[2L]])))
  t1 * y0 + theta * (y1 + t1 * (k[[1L]] + theta * (k[[2L]] + theta * k[[3L]])))
}

#' Rosenbrock algorithms
#'
#' Linearly implicit Rosenbrock methods for stiff problems and for
#' differential algebraic equations in mass matrix form (index 1, singular
#' mass matrices as produced by [mtkcompile()]), as in OrdinaryDiffEq.jl.
#' `Rosenbrock23()` is the order 2(3) method of MATLAB's ode23s with
#' Shampine's interpolant, `Rodas4()` an order 4 and `Rodas5P()` an order 5
#' stiffly accurate method with their stiff aware dense output. Each step
#' uses the Jacobian `prob$f$jac` (from `ODEProblem(sys, ..., jac = TRUE)`)
#' when present and finite differences otherwise, and `prob$f$tgrad` for
#' non autonomous problems when present.
#'
#' @param autodiff How Jacobians are computed when the problem has no
#'   symbolic Jacobian: `AutoForwardDiff()` (central differences, as R has no
#'   dual numbers), `AutoFiniteDiff(fdtype)` or `AutoSparse(ad)`; see
#'   [AutoForwardDiff()]. The default uses forward differences.
#' @return An algorithm object.
#' @examples
#' vdp <- function(u, p, t) c(u[2], p * (1 - u[1]^2) * u[2] - u[1])
#' prob <- ODEProblem(vdp, c(2, 0), c(0, 1), p = 1000)
#' sol <- solve(prob, Rodas5P())
#' sol$retcode
#' @name rosenbrock_algorithms
NULL

#' @rdname rosenbrock_algorithms
#' @export
Rosenbrock23 <- function(autodiff = NULL) .with_autodiff(autodiff, {
  .new_alg("Rosenbrock23", 2L, step = .rosenbrock23_step, interp = .rosenbrock23_interp,
           init = .ros_init, stiff = TRUE, mass_matrix = TRUE, qsteady_max = 6 / 5)
})

#' @rdname rosenbrock_algorithms
#' @export
Rodas4 <- function(autodiff = NULL) .with_autodiff(autodiff, {
  .new_alg("Rodas4", 4L, step = .rodas_step(.RODAS4), interp = .rodas_interp,
           init = .ros_init, stiff = TRUE, mass_matrix = TRUE, qsteady_max = 6 / 5)
})

#' @rdname rosenbrock_algorithms
#' @export
Rodas5P <- function(autodiff = NULL) .with_autodiff(autodiff, {
  .new_alg("Rodas5P", 5L, step = .rodas_step(.RODAS5P), interp = .rodas_interp,
           init = .ros_init, stiff = TRUE, mass_matrix = TRUE, qsteady_max = 6 / 5)
})
