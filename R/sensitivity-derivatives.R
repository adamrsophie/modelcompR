# Derivatives used by the sensitivity methods, after SciMLSensitivity.jl
# src/derivative_wrappers.jl (vecjacobian!, jacobian!, gradient!,
# accumulate_cost!), src/adjoint_common.jl (adjointdiffcache: the mass
# matrix structure) and DiffEqBase's UJacobianWrapper/ParamJacobianWrapper.
#
# Julia obtains Jacobian vector and vector Jacobian products by automatic
# differentiation of the user function. Here a problem built from a System
# gets them symbolically: the Jacobians df/du and df/dp of the compiled
# right-hand side, and the gradient of the scalar lambda' f with respect to
# the unknowns and the parameters, all from .partials() and compiled with
# .build_closure(). A problem built from a plain R function uses its
# Jacobian `jac_fn` when given and central finite differences otherwise.

# Central difference step: the cube root of the machine epsilon balances
# truncation and rounding errors.
.sens_h <- function(x) .Machine$double.eps^(1 / 3) * max(1, abs(x))

# Jacobian of a vector function F at x by central differences.
.sens_fd_jac <- function(F, x) {
  x <- as.numeric(x)
  cols <- lapply(seq_along(x), function(j) {
    h <- .sens_h(x[j])
    xp <- x; xm <- x
    xp[j] <- x[j] + h
    xm[j] <- x[j] - h
    (as.numeric(F(xp)) - as.numeric(F(xm))) / (2 * h)
  })
  if (!length(cols)) return(matrix(0, length(F(x)), 0L))
  do.call(cbind, cols)
}

.sens_fd_grad <- function(F, x) as.numeric(.sens_fd_jac(function(z) F(z), x))

# The derivative bundle of an ODE problem: f, the Jacobians and the vector
# Jacobian product vjp(l, u, p, t) = list(u = J' l, p = fp' l).
.sens_fns <- function(prob, sensealg = NULL) {
  f <- prob$f$f
  sys <- prob$sys
  M <- prob$f$mass_matrix
  if (!is.matrix(M) || .is_identity_mass(M)) M <- NULL
  n <- length(prob$u0)
  np <- length(prob$p)
  symbolic <- !is.null(sys) && !is.null(prob$f$sys) && !identical(sensealg$autojacvec, FALSE) &&
    inherits(sys, "System")
  out <- list(f = function(u, p, t) as.numeric(f(u, p, t)), n = n, np = np, mass = M,
              symbolic = symbolic)
  if (symbolic) {
    d <- .sens_symbolic(sys)
    out$jac <- d$jac
    out$pjac <- function(u, p, t) if (np) d$pjac(u, p, t) else matrix(0, n, 0L)
    out$vjp <- function(l, u, p, t) {
      v <- unname(d$vjp(l, u, p, t))
      list(u = v[seq_len(n)], p = v[n + seq_len(np)])
    }
    return(out)
  }
  ujac <- prob$f$jac
  out$jac <- if (is.function(ujac)) function(u, p, t) as.matrix(ujac(u, p, t)) else
    function(u, p, t) .sens_fd_jac(function(z) f(z, p, t), u)
  out$pjac <- function(u, p, t) {
    if (!np) return(matrix(0, n, 0L))
    J <- .sens_fd_jac(function(q) {
      qq <- p
      qq[] <- q
      f(u, qq, t)
    }, p)
    J
  }
  out$vjp <- function(l, u, p, t) list(u = as.numeric(crossprod(out$jac(u, p, t), l)),
                                       p = as.numeric(crossprod(out$pjac(u, p, t), l)))
  out
}

# Symbolic Jacobians and vector Jacobian product of a compiled system,
# cached by the keys of its right-hand sides, unknowns and parameters.
.sens_symbolic <- function(sys) {
  rhs <- .ode_rhs_exprs(sys)
  us <- unknowns(sys)
  ps <- .problem_params(sys)
  iv <- .f(sys, "iv")
  key <- paste(c(.keys(rhs), "|", .keys(us), "|", .keys(ps)), collapse = ";")
  if (is.null(.mtk$sens_cache)) .mtk$sens_cache <- new.env(parent = emptyenv())
  hit <- .mtk$sens_cache[[key]]
  if (!is.null(hit)) return(hit)
  others <- c(.parameters_of(sys, TRUE), if (!is.null(iv)) list(iv))
  args <- .codegen_args(sys, with_t = !is.null(iv))
  Ju <- .partials(rhs, us, others)
  jac <- .build_closure(Ju, args, out = "matrix")
  pjac <- if (length(ps)) .build_closure(.partials(rhs, ps, c(us, if (!is.null(iv)) list(iv))), args,
                                         out = "matrix")
  # the gradient of lambda' f by the unknowns and the parameters
  lam <- lapply(seq_along(rhs), function(i) Symbol(paste0("_sl", i)))
  L <- Reduce(`+`, Map(function(l, r) l * r, lam, rhs), .S0())
  wk <- .keys(c(us, ps))
  G <- .partials(list(L), c(us, ps), c(lam, Filter(function(v) !(.key(v) %in% wk), others)))
  gexprs <- lapply(seq_len(ncol(G)), function(j) G[[1L, j]])
  vargs <- c(list(list(name = "l", vars = lam)), args)
  vjp <- .build_closure(gexprs, vargs, out = "vector")
  if (is.null(iv)) {
    jac0 <- jac; pjac0 <- pjac; vjp0 <- vjp
    jac <- function(u, p, t) jac0(u, p)
    pjac <- if (!is.null(pjac0)) function(u, p, t) pjac0(u, p)
    vjp <- function(l, u, p, t) vjp0(l, u, p)
  }
  out <- list(jac = jac, pjac = pjac, vjp = vjp)
  assign(key, out, envir = .mtk$sens_cache)
  out
}

# ---------------------------------------------------------------------------
# Semi-explicit DAE structure of a constant mass matrix (adjointdiffcache):
# differential variables are the columns of M with a nonzero entry, the
# algebraic equations the rows without one; the differential block must be
# nonsingular.

.sens_mass_structure <- function(M, n) {
  if (is.null(M)) return(list(diff = seq_len(n), alg = integer(), Mdd = diag(1, n), dae = FALSE))
  dv <- which(colSums(abs(M)) != 0)
  av <- setdiff(seq_len(n), dv)
  ae <- which(rowSums(abs(M)) == 0)
  de <- setdiff(seq_len(n), ae)
  if (length(ae) != length(av))
    stop("the mass matrix must have as many zero rows as zero columns (a semi-explicit DAE)", call. = FALSE)
  Mdd <- M[de, dv, drop = FALSE]
  if (length(dv) && abs(det(Mdd)) < 1e-300)
    stop("The submatrix corresponding to the differential variables of the mass matrix must be nonsingular!",
         call. = FALSE)
  list(diff = dv, alg = av, diff_eq = de, alg_eq = ae, Mdd = Mdd, dae = length(av) > 0L)
}

# ---------------------------------------------------------------------------
# Forward trajectory: the state y(t) of the forward solution for the
# adjoint. The adjoint integrates each segment between events backward, so
# lookups are restricted to the current segment and the side of a
# discontinuity is chosen by the segment, not by comparing times.

.sens_trajectory <- function(sol, event_t = numeric(), checkpoint = NULL) {
  env <- .sol_env(sol)
  tr <- new.env(parent = emptyenv())
  tr$event_t <- sort(event_t)
  tr$seg <- length(event_t)
  if (is.null(checkpoint)) {
    steps <- env$steps
    if (is.null(steps) || !length(steps))
      stop("the adjoint needs a dense forward solution: solve without `saveat`, with `dense = TRUE`, ",
           "or with `sensealg = ...`, or use a checkpointing sensitivity algorithm", call. = FALSE)
    starts <- vapply(steps, `[[`, 0, "t0")
    by_seg <- split(seq_along(steps), factor(findInterval(starts, tr$event_t), levels = 0:length(event_t)))
    # a solution over a decreasing time span has decreasing step starts
    sgn <- if (length(starts) > 1L && starts[length(starts)] < starts[1L]) -1 else 1
    tr$at <- function(t) {
      idx <- by_seg[[tr$seg + 1L]]
      if (!length(idx)) idx <- seq_along(steps)
      k <- findInterval(sgn * t, sgn * starts[idx])
      i <- idx[max(1L, min(length(idx), k))]
      .step_value(steps[[i]], t)
    }
    return(tr)
  }
  # checkpointing: re-solve the interval between saved points containing t
  ck <- checkpoint
  tr$cache <- NULL
  tr$at <- function(t) {
    k <- max(1L, min(length(ck$times) - 1L, findInterval(t, ck$times, rightmost.closed = TRUE)))
    if (is.null(tr$cache) || tr$cache$k != k) {
      pr <- ck$prob
      pr$u0 <- ck$states[[k]]
      pr$tspan <- ck$times[c(k, k + 1L)]
      pr$kwargs$callback <- NULL
      pr$kwargs$saveat <- NULL
      s <- .integrate(pr, ck$alg, c(ck$opts, list(dense = TRUE, save_everystep = TRUE)))
      tr$cache <- list(k = k, steps = .sol_env(s)$steps)
      tr$cache$starts <- vapply(tr$cache$steps, `[[`, 0, "t0")
    }
    st <- tr$cache
    i <- max(1L, min(length(st$steps), findInterval(t, st$starts)))
    .step_value(st$steps[[i]], t)
  }
  tr
}

# ---------------------------------------------------------------------------
# Quadrature

# Gauss-Kronrod (7, 15) nodes and weights on [-1, 1].
.GK15 <- local({
  xk <- c(0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
          0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
          0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
          0.207784955007898467600689403773245, 0)
  wk <- c(0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
          0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
          0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
          0.204432940075298892414161999234649, 0.209482141084727828012999174891714)
  wg <- c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
          0.381830050505118944950369775488975, 0.417959183673469387755102040816327)
  list(x = c(-xk[1:7], xk[8], rev(xk[1:7])), wk = c(wk[1:7], wk[8], rev(wk[1:7])),
       wg = c(0, wg[1], 0, wg[2], 0, wg[3], 0, wg[4], 0, wg[3], 0, wg[2], 0, wg[1], 0))
})

# 7 point Gauss-Legendre rule on [-1, 1].
.GL7 <- list(x = c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788,
                   -0.405845151377397166906606412076961, 0, 0.405845151377397166906606412076961,
                   0.741531185599394439863864773280788, 0.949107912342758524526189684047851),
             w = c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                   0.381830050505118944950369775488975, 0.417959183673469387755102040816327,
                   0.381830050505118944950369775488975, 0.279705391489276667901467771423780,
                   0.129484966168869693270611432679082))

.sens_gk15 <- function(fn, a, b) {
  c0 <- (a + b) / 2; h <- (b - a) / 2
  vals <- lapply(c0 + h * .GK15$x, fn)
  K <- h * Reduce(`+`, Map(`*`, vals, .GK15$wk))
  G <- h * Reduce(`+`, Map(`*`, vals, .GK15$wg))
  list(value = K, err = max(abs(K - G)))
}

# Adaptive Gauss-Kronrod quadrature of a vector valued function, as QuadGK.
.sens_quadgk <- function(fn, a, b, abstol = 1e-6, reltol = 1e-3, maxdepth = 30L) {
  if (a == b) return(0 * fn(a))
  segs <- list(c(.sens_gk15(fn, a, b), list(a = a, b = b, depth = 0L)))
  repeat {
    total <- Reduce(`+`, lapply(segs, `[[`, "value"))
    err <- sum(vapply(segs, `[[`, 0, "err"))
    if (err <= max(abstol, reltol * max(abs(total)))) return(total)
    i <- which.max(vapply(segs, `[[`, 0, "err"))
    s <- segs[[i]]
    if (s$depth >= maxdepth) return(total)
    m <- (s$a + s$b) / 2
    segs[[i]] <- c(.sens_gk15(fn, s$a, m), list(a = s$a, b = m, depth = s$depth + 1L))
    segs[[length(segs) + 1L]] <- c(.sens_gk15(fn, m, s$b), list(a = m, b = s$b, depth = s$depth + 1L))
  }
}

.sens_gauss7 <- function(fn, a, b) {
  c0 <- (a + b) / 2; h <- (b - a) / 2
  h * Reduce(`+`, Map(function(x, w) w * fn(c0 + h * x), .GL7$x, .GL7$w))
}

# ---------------------------------------------------------------------------
# Initialization: the map from the problem parameters to the initial values
# (and to the parameters the initialization solves for), differentiated by
# central differences over remake(), which reruns the initialization of the
# problem. Problems without initialization data have u0 independent of p.

.sens_init_jacobian <- function(prob) {
  n <- length(prob$u0); np <- length(prob$p)
  none <- list(du0 = matrix(0, n, np), dp = diag(1, np), active = FALSE)
  if (is.null(prob$sys) || is.null(prob$initialization_data) || !np) return(none)
  # ponytail: finite differences through the initialization solve; an
  # implicit function theorem pullback of the initialization system would
  # avoid the 2 np re-initializations.
  map <- function(q) {
    pr <- suppressWarnings(remake(prob, p = q))
    c(as.numeric(pr$u0), as.numeric(pr$p))
  }
  J <- .sens_fd_jac(map, as.numeric(prob$p))
  du0 <- J[seq_len(n), , drop = FALSE]
  dp <- J[n + seq_len(np), , drop = FALSE]
  list(du0 = du0, dp = dp, active = TRUE)
}
