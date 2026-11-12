# scipy.optimize.root with the quasi Newton methods of
# scipy/optimize/_nonlin.py: nonlin_solve with the Jacobian approximations
# BroydenFirst, BroydenSecond, Anderson, LinearMixing, DiagBroyden and
# ExcitingMixing, the armijo and wolfe line searches and the termination
# test of TerminationCondition.

.nl_norm <- function(v) sqrt(sum(v * v))
.nl_maxnorm <- function(v) max(abs(v))

# LowRankMatrix: alpha I + sum_k c_k d_k', collapsed to a dense matrix when
# the rank exceeds the dimension
.nl_lowrank <- function(alpha, n) {
  e <- new.env(parent = emptyenv())
  e$alpha <- alpha; e$cs <- list(); e$ds <- list(); e$collapsed <- NULL
  mv <- function(v, alpha, cs, ds) {
    w <- alpha * v
    for (k in seq_along(cs)) w <- w + sum(ds[[k]] * v) * cs[[k]]
    w
  }
  sv <- function(v, alpha, cs, ds) {
    if (!length(cs)) return(v / alpha)
    m <- length(cs)
    A <- alpha * diag(m)
    for (i in seq_len(m)) for (j in seq_len(m)) A[i, j] <- A[i, j] + sum(ds[[i]] * cs[[j]])
    q <- vapply(ds, function(d) sum(d * v), 0) / alpha
    q <- solve(A, q)
    w <- v / alpha
    for (k in seq_len(m)) w <- w - q[k] * cs[[k]]
    w
  }
  e$matvec <- function(v) if (!is.null(e$collapsed)) as.numeric(e$collapsed %*% v) else mv(v, e$alpha, e$cs, e$ds)
  e$rmatvec <- function(v) if (!is.null(e$collapsed)) as.numeric(crossprod(e$collapsed, v)) else
    mv(v, e$alpha, e$ds, e$cs)
  e$solve <- function(v) if (!is.null(e$collapsed)) as.numeric(solve(e$collapsed, v)) else sv(v, e$alpha, e$cs, e$ds)
  e$append <- function(c, d) {
    if (!is.null(e$collapsed)) { e$collapsed <- e$collapsed + outer(c, d); return(invisible()) }
    e$cs[[length(e$cs) + 1L]] <- c
    e$ds[[length(e$ds) + 1L]] <- d
    if (length(e$cs) > length(c)) {
      G <- e$alpha * diag(n)
      for (k in seq_along(e$cs)) G <- G + outer(e$cs[[k]], e$ds[[k]])
      e$collapsed <- G; e$cs <- NULL; e$ds <- NULL; e$alpha <- NULL
    }
  }
  e$restart_reduce <- function(rank) {
    if (is.null(e$collapsed) && length(e$cs) > rank) { e$cs <- list(); e$ds <- list() }
  }
  e$simple_reduce <- function(rank) {
    if (!is.null(e$collapsed)) return(invisible())
    while (length(e$cs) > rank) { e$cs[[1L]] <- NULL; e$ds[[1L]] <- NULL }
  }
  e$svd_reduce <- function(max_rank, to_retain = NULL) {
    if (!is.null(e$collapsed)) return(invisible())
    p <- max_rank
    q <- to_retain %||% (p - 2)
    if (length(e$cs)) p <- min(p, length(e$cs[[1L]]))
    q <- max(0, min(q, p - 1))
    if (length(e$cs) < p) return(invisible())
    C <- do.call(cbind, e$cs); D <- do.call(cbind, e$ds)
    qd <- qr(D)
    R <- qr.R(qd)[, order(qd$pivot), drop = FALSE]
    D <- qr.Q(qd)
    C <- C %*% t(R)
    s <- svd(C)
    WH <- t(s$v)
    C <- C %*% solve(WH)
    D <- D %*% t(WH)
    e$cs <- lapply(seq_len(q), function(k) C[, k]); e$ds <- lapply(seq_len(q), function(k) D[, k])
  }
  e
}

# the Jacobian approximations; each has setup(x, F), solve(f) and update(x, F)
.scipy_nl_jacobian <- function(method, opts) {
  allowed <- switch(method, broyden1 = , broyden2 = c("alpha", "reduction_method", "max_rank"),
                    anderson = c("alpha", "w0", "M"), linearmixing = , diagbroyden = "alpha",
                    excitingmixing = c("alpha", "alphamax"))
  bad <- setdiff(names(opts), allowed)
  if (length(bad)) stop("__init__() got an unexpected keyword argument '", bad[1L], "'", call. = FALSE)
  J <- new.env(parent = emptyenv())
  J$alpha <- opts$alpha
  J$setup <- function(x, F) {
    J$last_f <- F; J$last_x <- x
    if (is.null(J$alpha)) {
      nf <- .nl_norm(F)
      J$alpha <- if (nf != 0) 0.5 * max(.nl_norm(x), 1) / nf else 1
    }
    J$init(x, F)
  }
  J$update <- function(x, f) {
    df <- f - J$last_f; dx <- x - J$last_x
    J$upd(x, f, dx, df, .nl_norm(dx), .nl_norm(df))
    J$last_f <- f; J$last_x <- x
  }
  J$init <- function(x, F) NULL
  if (method %in% c("broyden1", "broyden2")) {
    max_rank <- opts$max_rank %||% Inf
    rm <- opts$reduction_method %||% "restart"
    params <- if (is.character(rm) && length(rm) == 1L) list() else as.list(rm)[-1L]
    rm <- as.character(rm[[1L]])
    if (!rm %in% c("svd", "simple", "restart")) stop("Unknown rank reduction method '", rm, "'", call. = FALSE)
    params <- c(list(max_rank - 1), params)
    J$reduce <- function() do.call(J$Gm[[paste0(rm, "_reduce")]], params)
    J$init <- function(x, F) J$Gm <- .nl_lowrank(-J$alpha, length(F))
    J$solve <- function(f) {
      r <- J$Gm$matvec(f)
      if (!all(is.finite(r))) { J$setup(J$last_x, J$last_f); r <- J$Gm$matvec(f) }
      r
    }
    J$upd <- if (method == "broyden1") function(x, f, dx, df, dx_norm, df_norm) {
      J$reduce()
      v <- J$Gm$rmatvec(dx)
      c <- dx - J$Gm$matvec(df)
      J$Gm$append(c, v / sum(df * v))
    } else function(x, f, dx, df, dx_norm, df_norm) {
      J$reduce()
      c <- dx - J$Gm$matvec(df)
      J$Gm$append(c, df / df_norm^2)
    }
  } else if (method == "anderson") {
    J$M <- opts$M %||% 5L; J$w0 <- opts$w0 %||% 0.01
    J$dx <- list(); J$df <- list()
    J$solve <- function(f) {
      dx <- -J$alpha * f
      n <- length(J$dx)
      if (!n) return(dx)
      df_f <- vapply(J$df, function(d) sum(d * f), 0)
      gamma <- tryCatch(solve(J$a, df_f), error = function(e) NULL)
      if (is.null(gamma)) { J$dx <- list(); J$df <- list(); return(dx) }
      for (m in seq_len(n)) dx <- dx + gamma[m] * (J$dx[[m]] + J$alpha * J$df[[m]])
      dx
    }
    J$upd <- function(x, f, dx, df, dx_norm, df_norm) {
      if (J$M == 0) return(invisible())
      J$dx[[length(J$dx) + 1L]] <- dx; J$df[[length(J$df) + 1L]] <- df
      while (length(J$dx) > J$M) { J$dx[[1L]] <- NULL; J$df[[1L]] <- NULL }
      n <- length(J$dx)
      a <- matrix(0, n, n)
      for (i in seq_len(n)) for (j in i:n) a[i, j] <- (1 + if (i == j) J$w0^2 else 0) * sum(J$df[[i]] * J$df[[j]])
      a[lower.tri(a)] <- t(a)[lower.tri(a)]
      J$a <- a
    }
  } else if (method == "linearmixing") {
    J$solve <- function(f) -f * J$alpha
    J$upd <- function(...) NULL
  } else if (method == "diagbroyden") {
    J$init <- function(x, F) J$d <- rep(1 / J$alpha, length(F))
    J$solve <- function(f) -f / J$d
    J$upd <- function(x, f, dx, df, dx_norm, df_norm) J$d <- J$d - (df + J$d * dx) * dx / dx_norm^2
  } else if (method == "excitingmixing") {
    J$alphamax <- opts$alphamax %||% 1
    J$init <- function(x, F) J$beta <- rep(J$alpha, length(F))
    J$solve <- function(f) -f * J$beta
    J$upd <- function(x, f, dx, df, dx_norm, df_norm) {
      incr <- f * J$last_f > 0
      J$beta[incr] <- J$beta[incr] + J$alpha
      J$beta[!incr] <- J$alpha
      J$beta <- pmin(pmax(J$beta, 0), J$alphamax)
    }
  }
  J
}

# scalar_search_armijo of scipy/optimize/_linesearch.py
.sp_scalar_search_armijo <- function(phi, phi0, derphi0, c1 = 1e-4, alpha0 = 1, amin = 0) {
  phi_a0 <- phi(alpha0)
  if (isTRUE(phi_a0 <= phi0 + c1 * alpha0 * derphi0)) return(list(alpha = alpha0, phi = phi_a0))
  alpha1 <- -(derphi0) * alpha0^2 / 2 / (phi_a0 - phi0 - derphi0 * alpha0)
  phi_a1 <- phi(alpha1)
  if (isTRUE(phi_a1 <= phi0 + c1 * alpha1 * derphi0)) return(list(alpha = alpha1, phi = phi_a1))
  # comparisons with NaN are false, as in Python
  while (isTRUE(alpha1 > amin)) {
    factor <- alpha0^2 * alpha1^2 * (alpha1 - alpha0)
    a <- (alpha0^2 * (phi_a1 - phi0 - derphi0 * alpha1) - alpha1^2 * (phi_a0 - phi0 - derphi0 * alpha0)) / factor
    b <- (-alpha0^3 * (phi_a1 - phi0 - derphi0 * alpha1) + alpha1^3 * (phi_a0 - phi0 - derphi0 * alpha0)) / factor
    alpha2 <- (-b + sqrt(abs(b^2 - 3 * a * derphi0))) / (3 * a)
    phi_a2 <- phi(alpha2)
    if (isTRUE(phi_a2 <= phi0 + c1 * alpha2 * derphi0)) return(list(alpha = alpha2, phi = phi_a2))
    if (isTRUE((alpha1 - alpha2) > alpha1 / 2) || isTRUE((1 - alpha2 / alpha1) < 0.96)) alpha2 <- alpha1 / 2
    alpha0 <- alpha1; alpha1 <- alpha2
    phi_a0 <- phi_a1; phi_a1 <- phi_a2
  }
  list(alpha = NULL, phi = phi_a1)
}

.nl_line_search <- function(func, x, Fx, dx, search_type, rdiff = 1e-8, smin = 1e-2) {
  tmp <- list(s = 0, Fx = Fx, phi = .nl_norm(Fx)^2)
  s_norm <- .nl_norm(x) / .nl_norm(dx)
  phi <- function(s, store = TRUE) {
    if (s == tmp$s) return(tmp$phi)
    v <- func(x + s * dx)
    p <- if (all(is.finite(v))) .nl_norm(v)^2 else Inf
    if (store) tmp <<- list(s = s, Fx = v, phi = p)
    p
  }
  derphi <- function(s) {
    ds <- (abs(s) + s_norm + 1) * rdiff
    (phi(s + ds, store = FALSE) - phi(s)) / ds
  }
  if (search_type == "wolfe") {
    # arguments are forced here: phi() changes tmp
    phi0 <- tmp$phi
    derphi0 <- derphi(0)
    s <- .sp_dcsrch(phi, derphi, 1e-4, 0.9, 1e-2, smin, 50, 1, phi0, derphi0)$stp
  } else {
    phi0 <- tmp$phi
    s <- .sp_scalar_search_armijo(phi, phi0, -phi0, amin = smin)$alpha
  }
  if (is.null(s)) s <- 1
  xn <- x + s * dx
  Fn <- if (s == tmp$s) tmp$Fx else func(xn)
  list(s = s, x = xn, Fx = Fn, Fx_norm = .nl_norm(Fn))
}

# nonlin_solve(full_output = TRUE, raise_exception = FALSE) as
# _root_nonlin_solve calls it
.sp_nonlin_solve <- function(func, x0, method, o) {
  known <- c("nit", "disp", "maxiter", "ftol", "fatol", "xtol", "xatol", "tol_norm", "line_search", "jac_options")
  bad <- setdiff(names(o), known)
  if (length(bad)) warning("Unknown solver options: ", paste(bad, collapse = ", "), call. = FALSE)
  f_tol <- o$fatol %||% .Machine$double.eps^(1 / 3)
  f_rtol <- o$ftol %||% Inf
  x_tol <- o$xatol %||% Inf
  x_rtol <- o$xtol %||% Inf
  iter <- o$nit
  tol_norm <- o$tol_norm %||% .nl_maxnorm
  ls <- if ("line_search" %in% names(o)) o$line_search else "armijo"
  if (isTRUE(ls)) ls <- "armijo"
  if (isFALSE(ls)) ls <- NULL
  if (!is.null(ls) && !ls %in% c("armijo", "wolfe")) stop("Invalid line search", call. = FALSE)
  x <- as.numeric(x0)
  dx <- rep(Inf, length(x))
  Fx <- func(x)
  Fx_norm <- .nl_norm(Fx)
  J <- .scipy_nl_jacobian(method, o$jac_options %||% list())
  J$setup(x, Fx)
  maxiter <- o$maxiter %||% (if (!is.null(iter)) iter + 1L else 100L * (length(x) + 1L))
  gamma <- 0.9; eta_max <- 0.9999; eta_treshold <- 0.1; eta <- 1e-3
  iteration <- 0L
  f0_norm <- NULL
  status <- 0L
  for (k in seq_len(maxiter)) {
    iteration <- iteration + 1L
    f_norm <- tol_norm(Fx)
    if (is.null(f0_norm)) f0_norm <- f_norm
    status <- if (f_norm == 0) 1L else if (!is.null(iter)) 2L * (iteration > iter) else
      as.integer(f_norm <= f_tol && f_norm / f_rtol <= f0_norm && tol_norm(dx) <= x_tol &&
                 tol_norm(dx) / x_rtol <= tol_norm(x))
    if (status) break
    dx <- -J$solve(Fx)
    if (.nl_norm(dx) == 0)
      stop("Jacobian inversion yielded zero vector. This indicates a bug in the Jacobian approximation.",
           call. = FALSE)
    if (!is.null(ls)) {
      r <- .nl_line_search(func, x, Fx, dx, ls)
      x <- r$x; Fx <- r$Fx; Fx_norm_new <- r$Fx_norm
    } else {
      x <- x + dx; Fx <- func(x); Fx_norm_new <- .nl_norm(Fx)
    }
    J$update(x, Fx)
    eta_A <- gamma * Fx_norm_new^2 / Fx_norm^2
    eta <- if (gamma * eta^2 < eta_treshold) min(eta_max, eta_A) else min(eta_max, max(eta_A, gamma * eta^2))
    Fx_norm <- Fx_norm_new
  }
  if (!status) status <- 2L
  list(x = x, method = method, nit = iteration, fun = Fx, status = status, success = status == 1L,
       message = if (status == 1L) "A solution was found at the specified tolerance." else
         "The maximum number of iterations allowed has been reached.")
}
