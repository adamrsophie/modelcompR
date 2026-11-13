# IPX, the interior point solver of HiGHS 1.12.0 (highs/ipm/ipx), with its
# crossover to a basic solution, and the HiGHS glue around it
# (highs/ipm/IpxWrapper.cpp and the IPX conversions of
# lp_data/HighsSolution.cpp), as SciPy 1.18.1 runs it for
# linprog(method = "highs-ipm"). The entry point is .highs_ipx_solve(lp,
# options) on an LP in HiGHS form (row_lower <= A x <= row_upper,
# col_lower <= x <= col_upper); it returns what solveLpIpx() leaves in the
# HiGHS solver object. The options include run_centring (with
# max_centring_steps and centring_ratio_tolerance) and the settings of the
# MIP solver's analytic centre (zero costs, ipm_iteration_limit = 200,
# run_crossover = "off"). The simplex cleanup that HiGHS' solveLp() runs
# after an "Unknown" or "UnboundedOrInfeasible" IPX status, the KKT check of
# Highs::lpKktCheck(), presolve and the linprog interface live with the
# simplex code.
#
# Index conventions: C++ indices j = 0..n+m-1 become j = 1..n+m; a basis
# position p = 0..m-1 becomes p = 1..m. Matrices are dense (the LP sizes
# linprog meets in this package are small); the basis matrix is factorized
# by BASICLU (R/optimization-scipy-highs-basiclu.R) behind the LuUpdate
# contract of lu_update.h (.ipx_lu_new).

# IPX status and error codes (ipx_status.h)
.ipx_st <- list(not_run = 0L, solved = 1000L, stopped = 1005L, no_model = 1006L,
                optimal = 1L, imprecise = 2L, primal_infeas = 3L, dual_infeas = 4L,
                user_interrupt = 5L, time_limit = 6L, iter_limit = 7L,
                no_progress = 8L, failed = 9L, debug = 10L)
.ipx_err <- list(invalid_dimension = 103L, invalid_matrix = 104L, invalid_vector = 105L,
                 cr_iter_limit = 201L, cr_matrix_not_posdef = 202L,
                 cr_precond_not_posdef = 203L, cr_no_progress = 204L,
                 cr_inf_or_nan = 205L, basis_singular = 301L,
                 basis_too_ill_conditioned = 306L, time_interrupt = 999L)

.ipx_infnorm <- function(x) if (length(x)) max(abs(x)) else 0

# Floating point order of the reference build. The reference is SciPy's
# Linux x86-64 build of HiGHS: GCC without FMA instructions and without
# reassociation, so every a * b + c rounds the product, and every reduction
# adds in loop order. .ipx_madd() marks the multiply adds of the C code.
.ipx_madd <- function(a, b, c) a * b + c

# ipx::Dot (utils.cc) and loops compiled like it: products added in order in
# double precision (sum() accumulates in long double, which is 80 bits on x86)
.ipx_dot <- function(x, y) Reduce(`+`, x * y, 0)

# DotColumn for every column of the dense matrix A, in row order
.ipx_coldots <- function(A, v) {
  d <- numeric(ncol(A))
  for (i in seq_len(nrow(A))) d <- d + A[i, ] * v[i]
  d
}

# ScatterColumn of every column j with weight d[j] onto lhs, in column order
.ipx_scatter <- function(A, d, lhs) {
  for (j in seq_len(ncol(A))) if (d[j] != 0) lhs <- .ipx_madd(A[, j], d[j], lhs)
  lhs
}

# FindMaxAbs (utils.cc): first index of the maximum absolute entry
.ipx_find_max_abs <- function(x) {
  a <- abs(x)
  if (!length(a) || max(a) <= 0) return(1L)
  which.max(a)
}

# Sortperm (utils.cc): order of (value, index) pairs
.ipx_sortperm <- function(v, reverse = FALSE) {
  idx <- seq_along(v)
  if (reverse) order(v, idx, decreasing = TRUE) else order(v, idx)
}

# exponent of frexp() for x >= 0 (frexp(0) gives 0)
.ipx_frexp_exp <- function(x) {
  e <- numeric(length(x))
  pos <- x > 0
  if (any(pos)) {
    xp <- x[pos]
    ep <- floor(log2(xp)) + 1
    ep <- ep - (xp < 2^(ep - 1))
    ep <- ep + (xp >= 2^ep)
    e[pos] <- ep
  }
  e
}

# EquilibrationFactor (model.cc) with C integer division
.ipx_eq_factor <- function(expmin, expmax, e) {
  f <- rep(1, length(e))
  lo <- e < expmin
  hi <- e > expmax
  f[lo] <- 2^((expmin - e[lo] + 1) %/% 2)
  f[hi] <- 2^(-((e[hi] - expmax + 1) %/% 2))
  f
}

# Model::EquilibrateMatrix
.ipx_equilibrate <- function(A) {
  expmin <- 0; expmax <- 3
  nz <- A != 0
  e <- .ipx_frexp_exp(abs(A[nz]))
  if (!any(e < expmin | e > expmax)) return(list(A = A, colscale = NULL, rowscale = NULL))
  cs <- rep(1, ncol(A)); rs <- rep(1, nrow(A))
  for (round in 1:10) {
    aa <- abs(A)
    rowmax <- apply(aa, 1, max)
    colmax <- apply(aa, 2, max)
    fr <- .ipx_eq_factor(expmin, expmax, .ipx_frexp_exp(rowmax))
    fc <- .ipx_eq_factor(expmin, expmax, .ipx_frexp_exp(colmax))
    rs <- rs * fr
    cs <- cs * fc
    if (all(fr == 1) && all(fc == 1)) break
    A <- A * outer(fr, fc)
  }
  list(A = A, colscale = cs, rowscale = rs)
}

# ---------------------------------------------------------------------------
# Model (model.cc): scaling and the choice between the primal and the dual
# computational form AI x = b, lb <= x <= ub with AI = [A I].

.ipx_model_load <- function(A, rhs, ctype, offset, obj, lbu, ubu, par) {
  nc <- nrow(A); nv <- ncol(A)
  if (nc < 0 || nv <= 0) return(list(errflag = .ipx_err$invalid_dimension))
  bad_vec <- any(!is.finite(rhs)) || any(!is.finite(obj)) ||
    any(!is.finite(lbu) & lbu != -Inf) || any(!is.finite(ubu) & ubu != Inf) ||
    any(lbu > ubu) || any(!ctype %in% c("=", "<", ">"))
  if (bad_vec) return(list(errflag = .ipx_err$invalid_vector))
  if (any(!is.finite(A))) return(list(errflag = .ipx_err$invalid_matrix))
  M <- list(errflag = 0L, num_constr = nc, num_var = nv, ctype = ctype, offset = offset)
  M$boxed <- which(is.finite(lbu) & is.finite(ubu))
  M$norm_obj <- .ipx_infnorm(obj)
  M$norm_rhs <- max(.ipx_infnorm(rhs), abs(lbu[is.finite(lbu)]), abs(ubu[is.finite(ubu)]), 0)
  # ScaleModel: flip variables with only an upper bound, then equilibrate
  flipped <- which(is.finite(ubu) & is.infinite(lbu))
  slb <- lbu; sub <- ubu; sobj <- obj
  if (length(flipped)) {
    slb[flipped] <- -ubu[flipped]
    sub[flipped] <- Inf
    A[, flipped] <- -A[, flipped]
    sobj[flipped] <- -sobj[flipped]
  }
  M$flipped <- flipped
  eq <- if (par$scale >= 1) .ipx_equilibrate(A) else list(A = A)
  A <- eq$A
  M$colscale <- eq$colscale
  M$rowscale <- eq$rowscale
  srhs <- rhs
  if (!is.null(M$colscale)) {
    sobj <- sobj * M$colscale
    slb <- slb / M$colscale
    sub <- sub / M$colscale
  }
  if (!is.null(M$rowscale)) srhs <- srhs * M$rowscale
  M$As <- A; M$sobj <- sobj; M$srhs <- srhs; M$slb <- slb; M$sub <- sub
  dualize <- par$dualize
  if (dualize == -1) dualize <- nc > 2 * nv
  if (dualize == -2) dualize <- FALSE
  if (dualize) {
    m <- nv; n <- nc + length(M$boxed)
    nb <- length(M$boxed)
    Eb <- matrix(0, nv, nb)
    Eb[cbind(M$boxed, seq_len(nb))] <- -1
    AI <- cbind(t(A), Eb, diag(1, nv))
    b <- sobj
    cc <- c(-srhs, sub[M$boxed], ifelse(is.finite(slb), -slb, 0))
    lb <- ub <- numeric(n + m)
    lb[seq_len(nc)] <- ifelse(ctype == ">", 0, -Inf)
    ub[seq_len(nc)] <- ifelse(ctype == "<", 0, Inf)
    if (nb) { lb[nc + seq_len(nb)] <- 0; ub[nc + seq_len(nb)] <- Inf }
    lb[n + seq_len(m)] <- 0
    ub[n + seq_len(m)] <- ifelse(is.finite(slb), Inf, 0)
  } else {
    m <- nc; n <- nv
    AI <- cbind(A, diag(1, nc))
    b <- srhs
    cc <- c(sobj, numeric(nc))
    lb <- c(slb, ifelse(ctype == ">", -Inf, 0))
    ub <- c(sub, ifelse(ctype == "<", Inf, 0))
  }
  dimnames(AI) <- NULL
  M$dualized <- as.logical(dualize)
  M$m <- m; M$n <- n; M$AI <- AI; M$b <- b; M$c <- cc; M$lb <- lb; M$ub <- ub
  M$norm_c <- .ipx_infnorm(cc)
  M$norm_bounds <- max(.ipx_infnorm(b), abs(lb[is.finite(lb)]), abs(ub[is.finite(ub)]), 0)
  M
}

# DualizeBackInteriorSolution and ScaleBackInteriorSolution
.ipx_dualize_back_interior <- function(M, x, xl, xu, y, zl, zu) {
  n <- M$n; m <- M$m; nc <- M$num_constr; nv <- M$num_var; ct <- M$ctype
  if (M$dualized) {
    nb <- length(M$boxed)
    xx <- -y
    yy <- numeric(nc)
    yy[ct == "="] <- x[which(ct == "=")]
    yy[ct == "<"] <- -xu[which(ct == "<")]
    yy[ct == ">"] <- xl[which(ct == ">")]
    zl_u <- xl[n + seq_len(nv)]
    zl_u[!is.finite(M$slb)] <- 0
    zu_u <- numeric(nv)
    zu_u[M$boxed] <- xl[nc + seq_len(nb)]
    xl_u <- rep(Inf, nv)
    fin <- is.finite(M$slb)
    xl_u[fin] <- zl[n + which(fin)]
    xu_u <- rep(Inf, nv)
    xu_u[M$boxed] <- zl[nc + seq_len(nb)]
    slack <- numeric(nc)
    slack[ct == "<"] <- zu[which(ct == "<")]
    slack[ct == ">"] <- -zl[which(ct == ">")]
  } else {
    xx <- x[seq_len(nv)]
    yy <- y
    yy[ct == "<"] <- -zl[n + which(ct == "<")]
    yy[ct == ">"] <- zu[n + which(ct == ">")]
    zl_u <- zl[seq_len(nv)]; zu_u <- zu[seq_len(nv)]
    xl_u <- xl[seq_len(nv)]; xu_u <- xu[seq_len(nv)]
    slack <- numeric(nc)
    slack[ct == "<"] <- xl[n + which(ct == "<")]
    slack[ct == ">"] <- -xu[n + which(ct == ">")]
  }
  list(x = xx, xl = xl_u, xu = xu_u, slack = slack, y = yy, zl = zl_u, zu = zu_u)
}

.ipx_scale_back_interior <- function(M, s) {
  if (!is.null(M$colscale)) {
    cs <- M$colscale
    s$x <- s$x * cs; s$xl <- s$xl * cs; s$xu <- s$xu * cs
    s$zl <- s$zl / cs; s$zu <- s$zu / cs
  }
  if (!is.null(M$rowscale)) {
    s$y <- s$y * M$rowscale
    s$slack <- s$slack / M$rowscale
  }
  for (j in M$flipped) {
    s$x[j] <- -s$x[j]
    s$xu[j] <- s$xl[j]; s$xl[j] <- Inf
    s$zu[j] <- s$zl[j]; s$zl[j] <- 0
  }
  s
}

# Model::EvaluateInteriorSolution
.ipx_evaluate_interior <- function(M, it) {
  s <- .ipx_dualize_back_interior(M, it$x, it$xl, it$xu, it$y, it$zl, it$zu)
  fl <- is.finite(M$slb); fu <- is.finite(M$sub)
  rl <- numeric(M$num_var); ru <- numeric(M$num_var)
  rl[fl] <- M$slb[fl] - s$x[fl] + s$xl[fl]
  ru[fu] <- M$sub[fu] - s$x[fu] - s$xu[fu]
  # MultiplyWithScaledMatrix: by the columns of AI, which are the rows of
  # the user matrix when the model is dualized
  if (M$dualized) {
    At <- t(M$As)
    rb <- -.ipx_coldots(At, s$x)
    rc <- .ipx_scatter(At, -s$y, numeric(M$num_var))
  } else {
    rb <- .ipx_scatter(M$As, -s$x, numeric(nrow(M$As)))
    rc <- -.ipx_coldots(M$As, s$y)
  }
  rb <- (rb - s$slack) + M$srhs
  rc <- (rc - (s$zl - s$zu)) + M$sobj
  if (!is.null(M$colscale)) { rc <- rc / M$colscale; rl <- rl * M$colscale; ru <- ru * M$colscale }
  if (!is.null(M$rowscale)) rb <- rb / M$rowscale
  for (j in M$flipped) { rc[j] <- -rc[j]; ru[j] <- -rl[j]; rl[j] <- 0 }
  pres <- max(.ipx_infnorm(rb), .ipx_infnorm(rl), .ipx_infnorm(ru))
  dres <- .ipx_infnorm(rc)
  pobj <- M$offset + .ipx_dot(M$sobj, s$x)
  dobj <- M$offset + .ipx_dot(M$srhs, s$y)
  for (j in seq_len(M$num_var)) {
    if (fl[j]) dobj <- .ipx_madd(M$slb[j], s$zl[j], dobj)
    if (fu[j]) dobj <- .ipx_madd(-M$sub[j], s$zu[j], dobj)
  }
  list(abs_presidual = pres, abs_dresidual = dres,
       rel_presidual = pres / (1 + M$norm_rhs), rel_dresidual = dres / (1 + M$norm_obj),
       pobjval = pobj, dobjval = dobj,
       rel_objgap = (pobj - dobj) / .ipx_madd(0.5, abs(pobj + dobj), 1))
}

# DualizeBackBasicSolution, DualizeBackBasis, CorrectScaledBasicSolution
.ipx_basic_back <- function(M, x, y, z, status) {
  n <- M$n; nc <- M$num_constr; nv <- M$num_var
  if (M$dualized) {
    nb <- length(M$boxed)
    xx <- -y
    slack <- -z[seq_len(nc)]
    yy <- x[seq_len(nc)]
    zz <- x[n + seq_len(nv)]
    zz[M$boxed] <- zz[M$boxed] - x[nc + seq_len(nb)]
    cbasis <- ifelse(status[seq_len(nc)] == 0L, -1L, 0L)
    vbasis <- ifelse(status[n + seq_len(nv)] == 0L,
                     ifelse(is.finite(M$slb), -1L, -3L), 0L)
    kb <- status[nc + seq_len(nb)] == 0L
    vbasis[M$boxed[kb]] <- -2L
  } else {
    xx <- x[seq_len(nv)]
    slack <- x[n + seq_len(nc)]
    yy <- y
    zz <- z[seq_len(nv)]
    cbasis <- ifelse(status[n + seq_len(nc)] == 0L, 0L, -1L)
    vbasis <- status[seq_len(nv)]
  }
  # CorrectScaledBasicSolution
  xx[vbasis == -1L] <- M$slb[vbasis == -1L]
  xx[vbasis == -2L] <- M$sub[vbasis == -2L]
  zz[vbasis == 0L] <- 0
  slack[cbasis == -1L] <- 0
  yy[cbasis == 0L] <- 0
  list(x = xx, slack = slack, y = yy, z = zz, cbasis = cbasis, vbasis = vbasis)
}

.ipx_scale_back_basic <- function(M, s) {
  if (!is.null(M$colscale)) { s$x <- s$x * M$colscale; s$z <- s$z / M$colscale }
  if (!is.null(M$rowscale)) { s$y <- s$y * M$rowscale; s$slack <- s$slack / M$rowscale }
  if (length(M$flipped)) {
    s$x[M$flipped] <- -s$x[M$flipped]
    s$z[M$flipped] <- -s$z[M$flipped]
  }
  s
}

# Model::EvaluateBasicSolution
.ipx_evaluate_basic <- function(M, x, y, z, status) {
  s <- .ipx_basic_back(M, x, y, z, status)
  xinf <- numeric(M$num_var); zinf <- numeric(M$num_var)
  sinf <- numeric(M$num_constr); yinf <- numeric(M$num_constr)
  lo <- s$x < M$slb; xinf[lo] <- s$x[lo] - M$slb[lo]
  hi <- s$x > M$sub; xinf[hi] <- s$x[hi] - M$sub[hi]
  k <- s$vbasis != -1L & s$z > 0; zinf[k] <- s$z[k]
  k <- s$vbasis != -2L & s$z < 0; zinf[k] <- s$z[k]
  ct <- M$ctype
  k <- ct == "<" & s$slack < 0; sinf[k] <- s$slack[k]
  k <- ct == "<" & s$y > 0; yinf[k] <- s$y[k]
  k <- ct == ">" & s$slack > 0; sinf[k] <- s$slack[k]
  k <- ct == ">" & s$y < 0; yinf[k] <- s$y[k]
  inf <- .ipx_scale_back_basic(M, list(x = xinf, slack = sinf, y = yinf, z = zinf))
  list(primal_infeas = max(.ipx_infnorm(inf$x), .ipx_infnorm(inf$slack)),
       dual_infeas = max(.ipx_infnorm(inf$z), .ipx_infnorm(inf$y)))
}

# ---------------------------------------------------------------------------
# Iterate (iterate.cc). State codes: 1 barrier lb, 2 barrier ub, 3 barrier
# boxed, 4 free, 5 fixed, 6 implied lb, 7 implied ub, 8 implied eq.

.ipx_iterate_new <- function(M) {
  it <- new.env(parent = emptyenv())
  N <- M$n + M$m
  fl <- is.finite(M$lb); fu <- is.finite(M$ub)
  it$M <- M
  it$st <- ifelse(fl & fu, 3L, ifelse(fl, 1L, ifelse(fu, 2L, 4L)))
  it$x <- numeric(N); it$y <- numeric(M$m)
  it$xl <- ifelse(fl, 1, Inf); it$xu <- ifelse(fu, 1, Inf)
  it$zl <- ifelse(fl, 1, 0); it$zu <- ifelse(fu, 1, 0)
  it$bounds_measure <- 1 + M$norm_bounds
  it$costs_measure <- 1 + M$norm_c
  it$feasibility_tol <- 1e-6; it$optimality_tol <- 1e-8; it$start_crossover_tol <- -1
  it$postprocessed <- FALSE; it$evaluated <- FALSE
  it
}

.ipx_has_lb <- function(it) it$st == 1L | it$st == 3L
.ipx_has_ub <- function(it) it$st == 2L | it$st == 3L
.ipx_is_fixed <- function(it) it$st == 5L
.ipx_is_freest <- function(it) it$st == 4L | it$st >= 6L
.ipx_is_implied <- function(it) it$st >= 6L
.ipx_is_barrier <- function(it) it$st <= 3L

.ipx_it_initialize <- function(it, x, xl, xu, y, zl, zu) {
  M <- it$M
  it$x <- x; it$xl <- xl; it$xu <- xu; it$y <- y; it$zl <- zl; it$zu <- zu
  fl <- is.finite(M$lb); fu <- is.finite(M$ub)
  it$st <- ifelse(fl & fu, 3L, ifelse(fl, 1L, ifelse(fu, 2L, 4L)))
  it$evaluated <- FALSE; it$postprocessed <- FALSE
}

.ipx_it_update <- function(it, sp, dx, dxl, dxu, sd, dy, dzl, dzu) {
  if (!is.null(dx)) { k <- !.ipx_is_fixed(it); it$x[k] <- .ipx_madd(sp, dx[k], it$x[k]) }
  hl <- .ipx_has_lb(it); hu <- .ipx_has_ub(it)
  if (!is.null(dxl)) it$xl[hl] <- pmax(.ipx_madd(sp, dxl[hl], it$xl[hl]), 0x1.4484bfeebc2a0p-100)
  if (!is.null(dxu)) it$xu[hu] <- pmax(.ipx_madd(sp, dxu[hu], it$xu[hu]), 0x1.4484bfeebc2a0p-100)
  if (!is.null(dy)) it$y <- .ipx_madd(sd, dy, it$y)
  if (!is.null(dzl)) it$zl[hl] <- pmax(.ipx_madd(sd, dzl[hl], it$zl[hl]), 0x1.4484bfeebc2a0p-100)
  if (!is.null(dzu)) it$zu[hu] <- pmax(.ipx_madd(sd, dzu[hu], it$zu[hu]), 0x1.4484bfeebc2a0p-100)
  it$evaluated <- FALSE
}

.ipx_it_make_fixed <- function(it, j, value = NULL) {
  if (!is.null(value)) it$x[j] <- value
  it$xl[j] <- 0; it$xu[j] <- 0; it$zl[j] <- 0; it$zu[j] <- 0
  it$st[j] <- 5L; it$evaluated <- FALSE
}
.ipx_it_make_implied <- function(it, j, code) {
  it$xl[j] <- Inf; it$xu[j] <- Inf
  if (code == 8L) { it$zl[j] <- 0; it$zu[j] <- 0 }
  it$st[j] <- code; it$evaluated <- FALSE
}

.ipx_it_scaling <- function(it) {
  d <- 1 / sqrt(it$zl / it$xl + it$zu / it$xu)
  d[.ipx_is_fixed(it)] <- 0
  d[.ipx_is_freest(it)] <- Inf
  d
}

.ipx_it_evaluate <- function(it) {
  if (it$evaluated) return(invisible())
  M <- it$M; AI <- M$AI
  fixed <- .ipx_is_fixed(it); hl <- .ipx_has_lb(it); hu <- .ipx_has_ub(it)
  impl <- .ipx_is_implied(it)
  x <- it$x; y <- it$y; zl <- it$zl; zu <- it$zu
  aty <- .ipx_coldots(AI, y)
  rb <- .ipx_scatter(AI, -x, M$b)
  rc <- (M$c - zl + zu) - aty
  if (!it$postprocessed) rc[fixed] <- 0
  rl <- numeric(length(x)); ru <- numeric(length(x))
  rl[hl] <- M$lb[hl] - x[hl] + it$xl[hl]
  ru[hu] <- M$ub[hu] - x[hu] - it$xu[hu]
  it$rb <- rb; it$rc <- rc; it$rl <- rl; it$ru <- ru
  it$presidual <- max(.ipx_infnorm(rb), .ipx_infnorm(rl), .ipx_infnorm(ru))
  it$dresidual <- .ipx_infnorm(rc)
  lb <- M$lb; ub <- M$ub
  it$offset <- 0
  if (it$postprocessed) {
    it$pobjective <- M$offset + .ipx_dot(M$c, x)
    d <- M$offset + .ipx_dot(M$b, y)
    for (j in seq_along(x)) {
      if (is.finite(lb[j])) d <- .ipx_madd(lb[j], zl[j], d)
      if (is.finite(ub[j])) d <- .ipx_madd(-ub[j], zu[j], d)
    }
    it$dobjective <- d
  } else {
    p <- M$offset; off <- 0
    for (j in seq_along(x)) {
      if (!fixed[j]) p <- .ipx_madd(M$c[j], x[j], p) else off <- .ipx_madd(M$c[j], x[j], off)
      if (impl[j]) {
        p <- .ipx_madd(-(zl[j] - zu[j]), x[j], p)
        off <- .ipx_madd(zl[j] - zu[j], x[j], off)
      }
    }
    d <- M$offset + .ipx_dot(M$b, y)
    for (j in seq_along(x)) {
      if (hl[j]) d <- .ipx_madd(lb[j], zl[j], d)
      if (hu[j]) d <- .ipx_madd(-ub[j], zu[j], d)
      if (fixed[j]) d <- .ipx_madd(-x[j], aty[j], d)
    }
    it$pobjective <- p; it$offset <- off; it$dobjective <- d
  }
  cp <- 0
  for (j in which(hl)) cp <- .ipx_madd(it$xl[j], zl[j], cp)
  for (j in which(hu)) cp <- .ipx_madd(it$xu[j], zu[j], cp)
  nf <- sum(hl) + sum(hu)
  it$complementarity <- cp
  it$mu <- if (nf > 0) cp / nf else 0
  it$evaluated <- TRUE
  invisible()
}

.ipx_it_get <- function(it, name) { .ipx_it_evaluate(it); get(name, envir = it) }
.ipx_it_pobj_pp <- function(it) { .ipx_it_evaluate(it); it$pobjective + it$offset }
.ipx_it_dobj_pp <- function(it) { .ipx_it_evaluate(it); it$dobjective + it$offset }

.ipx_it_feasible <- function(it) {
  .ipx_it_evaluate(it)
  it$presidual <= it$feasibility_tol * it$bounds_measure &&
    it$dresidual <= it$feasibility_tol * it$costs_measure
}
.ipx_it_optimal <- function(it) {
  p <- .ipx_it_pobj_pp(it); d <- .ipx_it_dobj_pp(it)
  abs(p - d) <= it$optimality_tol * (1 + abs(0.5 * (p + d)))
}

# Iterate::ResidualsFromDropping
.ipx_it_drop_residuals <- function(it) {
  M <- it$M; st <- it$st
  N <- length(st)
  xdrop <- numeric(N); zdrop <- numeric(N)
  x <- it$x; xl <- it$xl; xu <- it$xu; zl <- it$zl; zu <- it$zu
  for (j in seq_len(N)) {
    s <- st[j]
    if (s == 3L) s <- if (zl[j] / xl[j] >= zu[j] / xu[j]) 1L else 2L
    if (s == 1L) {
      if (zl[j] >= xl[j]) xdrop[j] <- x[j] - M$lb[j] else zdrop[j] <- zl[j] - zu[j]
    } else if (s == 2L) {
      if (zu[j] >= xu[j]) xdrop[j] <- x[j] - M$ub[j] else zdrop[j] <- zl[j] - zu[j]
    }
  }
  amax <- apply(abs(M$AI), 2, max)
  list(pres = max(abs(xdrop) * amax, 0), dres = max(abs(zdrop), 0))
}

.ipx_it_term_crit <- function(it) {
  if (.ipx_it_feasible(it) && .ipx_it_optimal(it)) {
    if (it$start_crossover_tol <= 0) return(TRUE)
    r <- .ipx_it_drop_residuals(it)
    if (r$pres <= it$start_crossover_tol * (1 + it$M$norm_bounds) &&
        r$dres <= it$start_crossover_tol * (1 + it$M$norm_c)) return(TRUE)
  }
  FALSE
}

# Iterate::Postprocess
.ipx_it_postprocess <- function(it) {
  M <- it$M
  aty <- .ipx_coldots(M$AI, it$y)
  for (j in which(.ipx_is_fixed(it))) {
    it$xl[j] <- it$x[j] - M$lb[j]
    it$xu[j] <- M$ub[j] - it$x[j]
    if (M$lb[j] == M$ub[j]) {
      z <- M$c[j] - aty[j]
      if (z >= 0) it$zl[j] <- z else it$zu[j] <- -z
    }
  }
  for (j in which(.ipx_is_implied(it))) {
    z <- M$c[j] - aty[j]
    s <- it$st[j]
    if (s == 8L) {
      if (z >= 0) { it$zl[j] <- z; it$zu[j] <- 0 } else { it$zl[j] <- 0; it$zu[j] <- -z }
      it$x[j] <- M$lb[j]
    } else if (s == 6L) {
      it$zl[j] <- z; it$zu[j] <- 0; it$x[j] <- M$lb[j]
    } else {
      it$zl[j] <- 0; it$zu[j] <- -z; it$x[j] <- M$ub[j]
    }
    it$xl[j] <- it$x[j] - M$lb[j]
    it$xu[j] <- M$ub[j] - it$x[j]
  }
  it$postprocessed <- TRUE; it$evaluated <- FALSE
}

# Iterate::DropToComplementarity
.ipx_it_drop_to_complementarity <- function(it) {
  M <- it$M; lb <- M$lb; ub <- M$ub
  N <- length(it$x)
  x <- numeric(N); z <- numeric(N)
  for (j in seq_len(N)) {
    xlj <- it$xl[j]; xuj <- it$xu[j]; zlj <- it$zl[j]; zuj <- it$zu[j]
    xj <- min(max(it$x[j], lb[j]), ub[j])
    if (lb[j] == ub[j]) {
      x[j] <- lb[j]; z[j] <- zlj - zuj
    } else if (is.finite(lb[j]) && is.finite(ub[j])) {
      if (zlj * xuj >= zuj * xlj) {
        if (zlj >= xlj) { x[j] <- lb[j]; z[j] <- max(0, zlj - zuj) } else { x[j] <- xj; z[j] <- 0 }
      } else {
        if (zuj >= xuj) { x[j] <- ub[j]; z[j] <- min(0, zlj - zuj) } else { x[j] <- xj; z[j] <- 0 }
      }
    } else if (is.finite(lb[j])) {
      if (zlj >= xlj) { x[j] <- lb[j]; z[j] <- max(0, zlj - zuj) } else { x[j] <- xj; z[j] <- 0 }
    } else if (is.finite(ub[j])) {
      if (zuj >= xuj) { x[j] <- ub[j]; z[j] <- min(0, zlj - zuj) } else { x[j] <- xj; z[j] <- 0 }
    } else {
      x[j] <- xj; z[j] <- 0
    }
  }
  list(x = x, y = it$y, z = z)
}

# ---------------------------------------------------------------------------
# Conjugate residuals (conjugate_residuals.cc). Capply(v) returns
# list(lhs, dot = v'lhs); Papply likewise for the preconditioner.

.ipx_cr <- function(Capply, rhs, tol, resscale, maxiter, Papply = NULL) {
  m <- length(rhs)
  if (maxiter < 0) maxiter <- m + 100
  lhs <- numeric(m)
  residual <- rhs
  errflag <- 0L; iter <- 0L
  resnorm_of <- function(r) if (is.null(resscale)) .ipx_infnorm(r) else .ipx_infnorm(resscale * r)
  if (is.null(Papply)) {
    a <- Capply(residual); Cresidual <- a$lhs; cdot <- a$dot
    step <- residual; Cstep <- Cresidual
    repeat {
      if (resnorm_of(residual) <= tol) break
      if (iter == maxiter) { errflag <- .ipx_err$cr_iter_limit; break }
      if (cdot <= 0) { errflag <- .ipx_err$cr_matrix_not_posdef; break }
      alpha <- cdot / .ipx_dot(Cstep, Cstep)
      if (!is.finite(alpha)) { errflag <- .ipx_err$cr_inf_or_nan; break }
      lhs <- lhs + alpha * step
      residual <- residual - alpha * Cstep
      a <- Capply(residual); Cresidual <- a$lhs
      beta <- a$dot / cdot
      step <- residual + beta * step
      Cstep <- Cresidual + beta * Cstep
      cdot <- a$dot
      iter <- iter + 1L
    }
  } else {
    p <- Papply(residual); sresidual <- p$lhs; resnorm_pre <- p$dot
    a <- Capply(sresidual); Csres <- a$lhs; cdot <- a$dot
    step <- sresidual; Cstep <- Csres
    repeat {
      if (resnorm_of(residual) <= tol) break
      if (iter == maxiter) { errflag <- .ipx_err$cr_iter_limit; break }
      if (cdot <= 0) { errflag <- .ipx_err$cr_matrix_not_posdef; break }
      pc <- Papply(Cstep)
      if (pc$dot <= 0) { errflag <- .ipx_err$cr_precond_not_posdef; break }
      alpha <- cdot / pc$dot
      if (!is.finite(alpha)) { errflag <- .ipx_err$cr_inf_or_nan; break }
      lhs <- lhs + alpha * step
      residual <- residual - alpha * Cstep
      sresidual <- sresidual - alpha * pc$lhs
      a <- Capply(sresidual); Csres <- a$lhs
      beta <- a$dot / cdot
      step <- sresidual + beta * step
      Cstep <- Csres + beta * Cstep
      cdot <- a$dot
      iter <- iter + 1L
      if (iter %% 5L == 0L) {
        p <- Papply(residual); sresidual <- p$lhs
        if (p$dot >= resnorm_pre) { errflag <- .ipx_err$cr_no_progress; break }
        resnorm_pre <- p$dot
      }
    }
  }
  list(lhs = lhs, iter = iter, errflag = errflag)
}

# KKTSolverDiag (kkt_solver_diag.cc) with NormalMatrix and DiagonalPrecond:
# normal equations AI W AI' preconditioned by their diagonal.
.ipx_kkt_diag <- function(M, info, maxiter = -1L) {
  k <- new.env(parent = emptyenv())
  n <- M$n; m <- M$m
  A <- M$AI[, seq_len(n), drop = FALSE]
  k$maxiter <- maxiter
  k$iter <- 0L
  k$basis <- NULL
  k$factorize <- function(it) {
    k$iter <- 0L
    if (!is.null(it)) {
      g <- it$zl / it$xl + it$zu / it$xu
      regval <- .ipx_it_get(it, "mu")
      gp <- g[g != 0]
      if (length(gp)) regval <- min(regval, gp)
      W <- 1 / g
      W[is.infinite(W)] <- 1 / regval
    } else {
      W <- rep(1, n + m)
    }
    k$W <- W
    k$resscale <- 1 / sqrt(W[n + seq_len(m)])
    dg <- W[n + seq_len(m)]
    for (j in seq_len(n)) dg <- .ipx_madd(A[, j] * W[j], A[, j], dg)
    k$diagonal <- dg
  }
  k$solve <- function(a, b, tol) {
    W <- k$W
    rhs <- .ipx_scatter(M$AI, W * a, -b)
    Wn <- W[seq_len(n)]; Ws <- W[n + seq_len(m)]
    Capply <- function(v) {
      lhs <- .ipx_scatter(A, .ipx_coldots(A, v) * Wn, v * Ws)
      list(lhs = lhs, dot = .ipx_dot(v, lhs))
    }
    Papply <- function(v) { lhs <- v / k$diagonal; list(lhs = lhs, dot = .ipx_dot(lhs, v)) }
    cr <- .ipx_cr(Capply, rhs, tol, k$resscale, k$maxiter, Papply)
    info$errflag <- cr$errflag
    info$kktiter1 <- info$kktiter1 + cr$iter
    k$iter <- k$iter + cr$iter
    y <- cr$lhs
    xs <- Wn * (a[seq_len(n)] - .ipx_coldots(A, y))
    x <- c(xs, .ipx_scatter(A, -xs, b))
    list(x = x, y = y)
  }
  k
}

# ---------------------------------------------------------------------------
# Basis (basis.cc). map[j]: -1 nonbasic, -2 nonbasic fixed, 1..m basic at
# that position, m+1..2m basic free (position map-m).

.ipx_basis_new <- function(M, par, info) {
  b <- new.env(parent = emptyenv())
  b$M <- M; b$par <- par; b$info <- info
  b$lu <- .ipx_lu_new(M$m, par$lu_pivottol)
  .ipx_basis_slack(b)
  b
}

.ipx_b_pos <- function(b, j) { v <- b$map[j]; if (v > b$M$m) v - b$M$m else if (v >= 1) v else 0L }
.ipx_b_isbasic <- function(b, j) b$map[j] >= 1L
.ipx_b_isnonbasic <- function(b, j) b$map[j] < 0L

.ipx_basis_slack <- function(b) {
  n <- b$M$n; m <- b$M$m
  b$basis <- n + seq_len(m)
  b$map <- c(rep(-1L, n), seq_len(m))
  .ipx_basis_factorize(b)
}

# basis matrix by columns for BASICLU; empty slots give empty columns
.ipx_basis_matrix <- function(b) {
  AI <- b$M$AI
  lapply(b$basis, function(j) {
    if (j <= 0) return(list(i = integer(0), x = numeric(0)))
    col <- AI[, j]
    nz <- which(col != 0)
    list(i = nz - 1L, x = col[nz])
  })
}

.ipx_basis_adapt <- function(b, dep) {
  n <- b$M$n
  for (t in seq_len(nrow(dep))) {
    p <- dep[t, 1]; i <- dep[t, 2]
    jb <- b$basis[p]; jn <- n + i
    b$basis[p] <- jn
    b$map[jn] <- p
    if (jb > 0) b$map[jb] <- -1L
  }
  nrow(dep)
}

.ipx_basis_factorize <- function(b) {
  err <- 0L
  repeat {
    f <- b$lu$factorize(.ipx_basis_matrix(b), FALSE)
    if (bitwAnd(f$flag, 2L)) {
      .ipx_basis_adapt(b, f$dependent)
      err <- .ipx_err$basis_singular
      break
    }
    if (bitwAnd(f$flag, 1L) && .ipx_basis_tighten(b)) next
    break
  }
  b$fresh <- TRUE
  err
}

.ipx_basis_crash_factorize <- function(b) {
  f <- b$lu$factorize(.ipx_basis_matrix(b), TRUE)
  nd <- 0L
  if (bitwAnd(f$flag, 2L)) nd <- .ipx_basis_adapt(b, f$dependent)
  b$fresh <- TRUE
  nd
}

.ipx_basis_tighten <- function(b) {
  tol <- b$lu$get_pivottol()
  if (tol <= 0.05) b$lu$set_pivottol(0.1)
  else if (tol <= 0.25) b$lu$set_pivottol(0.3)
  else if (tol <= 0.5) b$lu$set_pivottol(0.9)
  else return(FALSE)
  TRUE
}

.ipx_basis_solve <- function(b, rhs, trans = "N") b$lu$solve_dense(rhs, trans)

# SolveForUpdate: ftran of column j if nonbasic, btran of its position if
# basic. With want = FALSE (the overload without lhs) only the data for the
# update are computed. The pattern of the solution, in BASICLU's order, is
# left in b$lu$pattern.
.ipx_basis_sfu <- function(b, j, want = TRUE) {
  p <- .ipx_b_pos(b, j)
  if (p > 0L) return(b$lu$btran(p, want))
  col <- b$M$AI[, j]
  nz <- which(col != 0)
  b$lu$ftran(nz, col[nz], want)
}

# The order in which for_each_nonzero visits an IndexedVector of dimension
# dim: its pattern if known and sparse (IndexedVector::sparse), else all
# indices.
.ipx_iv_order <- function(pattern, dim) {
  if (!is.null(pattern) && length(pattern) <= 0.1 * dim) pattern else seq_len(dim)
}

# First index in visiting order ord at which v attains a maximum larger
# than thr (the strict > updates of the upstream loops), 0 if none
.ipx_first_max <- function(v, ord, thr = 0) {
  if (!length(ord)) return(0L)
  vo <- v[ord]
  k <- which.max(vo)
  if (length(k) && vo[k] > thr) ord[k] else 0L
}

# Basis::TableauRow. Returns the btran solution and its pattern and the
# tableau row with its pattern (NULL when computed densely).
.ipx_basis_tableau_row <- function(b, jb, ignore_fixed = FALSE) {
  M <- b$M; AI <- M$AI; n <- M$n; m <- M$m
  btran <- .ipx_basis_sfu(b, jb)
  bpat <- b$lu$pattern
  elig <- b$map == -1L | (b$map == -2L & !ignore_fixed)
  is_sparse <- length(bpat) <= 0.1 * m
  if (is_sparse) {
    nz <- 0L
    for (i in bpat) nz <- nz + sum(AI[i, ] != 0)
    nz <- nz %/% 2L
    if (nz > 0.1 * n) is_sparse <- FALSE
  }
  row <- numeric(n + m)
  if (is_sparse) {
    rpat <- integer(0)
    marked <- logical(n + m)
    for (i in bpat) {
      temp <- btran[i]
      for (j in which(AI[i, ] != 0)) {
        if (elig[j] && !marked[j]) { marked[j] <- TRUE; rpat <- c(rpat, j) }
        if (marked[j]) row[j] <- .ipx_madd(temp, AI[i, j], row[j])
      }
    }
  } else {
    rpat <- NULL
    k <- which(elig)
    row[k] <- .ipx_coldots(AI[, k, drop = FALSE], btran)
  }
  list(btran = btran, btran_pattern = bpat, row = row, row_pattern = rpat)
}

# ExchangeIfStable: returns list(err, exchanged)
.ipx_basis_exchange <- function(b, jb, jn, tbl, sys) {
  if (sys > 0) .ipx_basis_sfu(b, jn, FALSE)
  if (sys < 0) .ipx_basis_sfu(b, jb, FALSE)
  err <- b$lu$update(tbl)
  if (err != 0L) {
    if (b$fresh && !.ipx_basis_tighten(b))
      return(list(err = .ipx_err$basis_too_ill_conditioned, exchanged = FALSE))
    return(list(err = .ipx_basis_factorize(b), exchanged = FALSE))
  }
  ib <- .ipx_b_pos(b, jb)
  b$basis[ib] <- jn
  b$map[jn] <- ib
  b$map[jb] <- -1L
  b$fresh <- FALSE
  if (b$lu$need_fresh()) return(list(err = .ipx_basis_factorize(b), exchanged = TRUE))
  list(err = 0L, exchanged = TRUE)
}

.ipx_basis_crash_exchange <- function(b, jb, jn, tbl) {
  ib <- .ipx_b_pos(b, jb)
  b$basis[ib] <- jn
  b$map[jn] <- ib
  b$map[jb] <- -1L
  b$fresh <- FALSE
  err <- b$lu$update(tbl)
  if (err != 0L || b$lu$need_fresh()) .ipx_basis_crash_factorize(b)
}

# ComputeBasicSolution
.ipx_basis_basic_solution <- function(b, x, z) {
  M <- b$M; AI <- M$AI
  nb <- b$map < 0L
  y <- .ipx_scatter(AI, ifelse(nb, -x, 0), M$b)
  x[b$basis] <- .ipx_basis_solve(b, y, "N")
  y <- .ipx_basis_solve(b, M$c[b$basis] - z[b$basis], "T")
  z[nb] <- M$c[nb] - .ipx_coldots(AI[, nb, drop = FALSE], y)
  list(x = x, y = y, z = z)
}

# guess_basis.cc
.ipx_guess_basis <- function(M, w) {
  m <- M$m; n <- M$n; AI <- M$AI
  N <- n + m
  colpat <- lapply(seq_len(N), function(j) which(AI[, j] != 0))
  basis <- integer()
  rownumber <- rep(-1L, m)
  active <- rep(TRUE, N)
  # ProcessFreeColumns: incomplete left-looking LU of the free columns
  Li <- list(); Lx <- list()
  for (j in which(w == Inf)) {
    marked <- rep(FALSE, m)
    fin <- integer()
    for (i0 in colpat[[j]]) {
      if (marked[i0]) next
      istack <- i0; pstack <- 0L
      while (length(istack)) {
        h <- length(istack)
        i <- istack[h]
        kcol <- rownumber[i]
        if (!marked[i]) { marked[i] <- TRUE; pstack[h] <- 1L }
        ch <- if (kcol >= 1L) Li[[kcol]] else integer()
        done <- TRUE
        p <- pstack[h]
        while (p <= length(ch)) {
          inext <- ch[p]
          if (marked[inext]) { p <- p + 1L; next }
          pstack[h] <- p + 1L
          istack <- c(istack, inext); pstack <- c(pstack, 0L)
          done <- FALSE
          break
        }
        if (done) {
          istack <- istack[-h]; pstack <- pstack[-h]
          fin <- c(fin, i)
        }
      }
    }
    pattern <- rev(fin)
    lhs <- numeric(m)
    lhs[colpat[[j]]] <- AI[colpat[[j]], j]
    lhsmax <- 0; imax <- 0L
    for (i in pattern) {
      temp <- lhs[i]; kcol <- rownumber[i]
      if (temp != 0) {
        if (kcol >= 1L) {
          lhs[Li[[kcol]]] <- lhs[Li[[kcol]]] - Lx[[kcol]] * temp
        } else if (abs(temp) > lhsmax) {
          lhsmax <- abs(temp); imax <- i
        }
      }
    }
    pivot <- if (imax > 0L) lhs[imax] else 0
    if (abs(pivot) > 1e-3) {
      basis <- c(basis, j)
      rownumber[imax] <- length(basis)
      rr <- colpat[[j]]
      rr <- rr[rownumber[rr] < 0L & lhs[rr] != 0]
      Li[[length(basis)]] <- rr
      Lx[[length(basis)]] <- lhs[rr] / pivot
    }
    active[j] <- FALSE
  }
  # ProcessSingletons
  ncol_nz <- vapply(colpat, length, 1L)
  for (i in seq_len(m)) {
    if (rownumber[i] >= 0L) next
    rowmax <- 0; maxs <- 0; js <- 0L
    for (j in which(AI[i, ] != 0)) {
      if (!active[j]) next
      a <- abs(AI[i, j]) * w[j]
      rowmax <- max(rowmax, a)
      if (a > maxs && ncol_nz[j] == 1L) { maxs <- a; js <- j }
    }
    if (maxs > 0 && maxs >= 0.5 * rowmax) {
      basis <- c(basis, js)
      rownumber[i] <- length(basis)
      active[js] <- FALSE
    }
  }
  # ProcessRemaining: matching in decreasing order of weight
  colperm <- .ipx_sortperm(w, reverse = TRUE)
  jmatch <- ifelse(rownumber >= 0L, -2L, -1L)
  marked <- rep(-1L, N)
  cheap <- rep(1L, N)
  num_failed <- 0L
  augment <- function(jstart) {
    jstack <- jstart; istack <- integer(); pstack <- integer()
    found <- FALSE
    head <- 1L
    while (head >= 1L) {
      j <- jstack[head]; pj <- colpat[[j]]
      if (marked[j] != jstart) {
        marked[j] <<- jstart
        p <- cheap[j]; i <- 0L
        while (p <= length(pj) && !found) { i <- pj[p]; found <- jmatch[i] == -1L; p <- p + 1L }
        cheap[j] <<- p
        if (found) { istack[head] <- i; break }
        pstack[head] <- 1L
      }
      p <- pstack[head]; adv <- FALSE
      while (p <= length(pj)) {
        i <- pj[p]
        if (jmatch[i] < -1L || marked[jmatch[i]] == jstart) { p <- p + 1L; next }
        pstack[head] <- p + 1L
        istack[head] <- i
        head <- head + 1L
        jstack[head] <- jmatch[i]
        adv <- TRUE
        break
      }
      if (!adv) head <- head - 1L
    }
    if (found) for (h in head:1L) jmatch[istack[h]] <<- jstack[h]
    found
  }
  for (j in colperm) {
    if (!active[j]) next
    if (w[j] == 0) break
    if (augment(j)) basis <- c(basis, j) else num_failed <- num_failed + 1L
    if (num_failed >= 10L * (m - length(basis))) break
  }
  rownumber[jmatch >= 1L] <- m
  c(basis, n + which(rownumber < 0L))
}

# InverseSearch: rook search for a large entry of inverse(B)
.ipx_inverse_search <- function(b) {
  m <- b$M$m
  work <- 1 / seq_len(m)
  inverse_max <- 0
  repeat {
    work <- .ipx_basis_solve(b, work, "N")
    if (any(!is.finite(work))) break
    pmax <- .ipx_find_max_abs(work)
    work <- numeric(m); work[pmax] <- 1
    work <- .ipx_basis_solve(b, work, "T")
    if (any(!is.finite(work))) break
    imax <- .ipx_find_max_abs(work)
    entry <- work[imax]
    if (abs(entry) <= 2 * inverse_max) return(c(pmax, imax, entry))
    inverse_max <- abs(entry)
    work <- numeric(m); work[imax] <- 1
  }
  c(0, 0, Inf)
}

.ipx_basis_repair <- function(b) {
  n <- b$M$n
  b$info$basis_repairs <- 0L
  repeat {
    e <- .ipx_inverse_search(b)
    pmax <- e[1]; imax <- e[2]; pivot <- e[3]
    if (pmax < 1 || imax < 1 || !is.finite(pivot)) { b$info$basis_repairs <- -1L; break }
    if (abs(pivot) < 1e5) break
    jb <- b$basis[pmax]; jn <- n + imax
    if (!.ipx_b_isnonbasic(b, jn)) { b$info$basis_repairs <- -2L; break }
    if (b$info$basis_repairs >= 200L) { b$info$basis_repairs <- -3L; break }
    .ipx_basis_sfu(b, jb, FALSE)
    .ipx_basis_sfu(b, jn, FALSE)
    .ipx_basis_crash_exchange(b, jb, jn, pivot)
    b$info$basis_repairs <- b$info$basis_repairs + 1L
  }
}

.ipx_pivot_free_in <- function(b, w) {
  M <- b$M; n <- M$n; m <- M$m; info <- b$info
  dtol <- max(0, b$par$dependency_tol)
  info$errflag <- 0L; info$dependent_cols <- 0L
  remaining <- which(is.infinite(w) & b$map < 0L)
  while (length(remaining)) {
    jn <- remaining[length(remaining)]
    ftran <- .ipx_basis_sfu(b, jn)
    ord <- .ipx_iv_order(b$lu$pattern, m)
    af <- abs(ftran)
    nonfree <- !is.infinite(w[b$basis])
    pmax <- .ipx_first_max(af, ord)
    fmax <- if (pmax > 0L) af[pmax] else 0
    afn <- ifelse(nonfree, af, 0)
    pmax_nf <- .ipx_first_max(afn, ord)
    fmax_nf <- if (pmax_nf > 0L) afn[pmax_nf] else 0
    if (fmax > 4 && fmax_nf < 1) {
      jb <- b$basis[pmax]
      ex <- .ipx_basis_exchange(b, jb, jn, ftran[pmax], -1L)
      info$errflag <- ex$err
      if (ex$err != 0L) return(invisible())
      if (!ex$exchanged) next
      remaining[length(remaining)] <- jb
      info$updates_start <- info$updates_start + 1L
    } else if (fmax_nf <= dtol) {
      if (!info$cols_inconsistent) {
        delta <- M$c[jn]
        for (p in ord) if (!nonfree[p]) delta <- .ipx_madd(-M$c[b$basis[p]], ftran[p], delta)
        if (abs(delta) > dtol) info$cols_inconsistent <- TRUE
      }
      info$dependent_cols <- info$dependent_cols + 1L
      remaining <- remaining[-length(remaining)]
    } else {
      jb <- b$basis[pmax_nf]
      ex <- .ipx_basis_exchange(b, jb, jn, ftran[pmax_nf], -1L)
      info$errflag <- ex$err
      if (ex$err != 0L) return(invisible())
      if (!ex$exchanged) next
      remaining <- remaining[-length(remaining)]
      info$updates_start <- info$updates_start + 1L
    }
  }
}

.ipx_pivot_fixed_out <- function(b, w) {
  M <- b$M; n <- M$n; m <- M$m; info <- b$info
  dtol <- max(0, b$par$dependency_tol)
  info$errflag <- 0L; info$dependent_rows <- 0L
  fx <- M$lb == M$ub & M$lb != 0
  bmf <- .ipx_scatter(M$AI, ifelse(fx, -M$lb, 0), M$b)
  remaining <- n + which(w[n + seq_len(m)] == 0 & b$map[n + seq_len(m)] >= 1L)
  elig <- seq_len(n + m) > n | w != 0
  while (length(remaining)) {
    jb <- remaining[length(remaining)]
    tr <- .ipx_basis_tableau_row(b, jb)
    row <- tr$row
    ord <- .ipx_iv_order(tr$row_pattern, n + m)
    ar <- ifelse(elig, abs(row), 0)
    jmax <- .ipx_first_max(ar, ord)
    rmax <- if (jmax > 0L) ar[jmax] else 0
    arn <- ifelse(w != 0, ar, 0)
    rmax_nf <- max(arn[ord], 0)
    if (rmax > 4 && rmax_nf < 1) {
      ex <- .ipx_basis_exchange(b, jb, jmax, row[jmax], 1L)
      info$errflag <- ex$err
      if (ex$err != 0L) return(invisible())
      if (!ex$exchanged) next
      remaining[length(remaining)] <- jmax
      info$updates_start <- info$updates_start + 1L
    } else if (rmax_nf <= dtol) {
      if (!info$rows_inconsistent) {
        # Dot(btran, b_minus_fixed_columns): over the sparse pattern, else
        # over all entries
        bo <- .ipx_iv_order(tr$btran_pattern, m)
        if (length(bo) < m || m == 0L) {
          delta <- 0
          for (i in bo) delta <- .ipx_madd(bmf[i], tr$btran[i], delta)
        } else {
          delta <- .ipx_dot(bmf, tr$btran)
        }
        if (abs(delta) > dtol) info$rows_inconsistent <- TRUE
      }
      info$dependent_rows <- info$dependent_rows + 1L
      remaining <- remaining[-length(remaining)]
    } else {
      ra <- abs(row)
      sc <- ra * w
      okk <- ra >= 0.1 * rmax_nf & !is.na(sc)
      sc[!okk] <- 0
      js <- .ipx_first_max(sc, ord)
      ex <- .ipx_basis_exchange(b, jb, js, row[js], 1L)
      info$errflag <- ex$err
      if (ex$err != 0L) return(invisible())
      if (!ex$exchanged) next
      remaining <- remaining[-length(remaining)]
      info$updates_start <- info$updates_start + 1L
    }
  }
}

.ipx_basis_from_weights <- function(b, w) {
  info <- b$info
  info$errflag <- 0L; info$dependent_rows <- 0L; info$dependent_cols <- 0L
  if (b$par$crash_basis) {
    b$basis <- .ipx_guess_basis(b$M, w)
    b$map <- rep(-1L, b$M$n + b$M$m)
    b$map[b$basis] <- seq_along(b$basis)
    .ipx_basis_crash_factorize(b)
    .ipx_basis_repair(b)
    if (info$basis_repairs < 0L) .ipx_basis_slack(b)
  } else {
    .ipx_basis_slack(b)
  }
  .ipx_pivot_free_in(b, w)
  if (info$errflag) return(invisible())
  .ipx_pivot_fixed_out(b, w)
}

# starting_basis.cc
.ipx_starting_basis <- function(it, b, info) {
  M <- it$M; n <- M$n; m <- M$m
  info$errflag <- 0L
  w <- .ipx_it_scaling(it)
  w[M$lb == M$ub] <- 0
  .ipx_basis_from_weights(b, w)
  if (info$errflag) return(invisible())
  for (j in which(w == 0 | is.infinite(w))) {
    if (b$map[j] >= 1L) {
      if (b$map[j] <= m) b$map[j] <- b$map[j] + m
    } else {
      b$map[j] <- -2L
    }
  }
  for (j in which(M$lb == M$ub & b$map == -2L)) .ipx_it_make_fixed(it, j, M$lb[j])
  # PostprocessDependencies
  N <- n + m
  dx <- numeric(N); dy <- numeric(m)
  dep_cols <- integer(); dep_rows <- integer()
  if (info$dependent_cols > 0L) {
    dxb <- numeric(m)
    for (j in seq_len(n)) {
      if (is.infinite(M$lb[j]) && is.infinite(M$ub[j]) && b$map[j] < 0L) {
        dx[j] <- -it$x[j]
        dxb <- dxb + M$AI[, j] * it$x[j]
        dep_cols <- c(dep_cols, j)
      }
    }
    dxb <- .ipx_basis_solve(b, dxb, "N")
    dx[b$basis] <- dxb
  }
  if (info$dependent_rows > 0L) {
    for (p in seq_len(m)) {
      j <- b$basis[p]
      if (j > n && M$lb[j] == M$ub[j]) {
        dy[p] <- -it$y[j - n]
        dep_rows <- c(dep_rows, j - n)
      }
    }
    dy <- .ipx_basis_solve(b, dy, "T")
    dy[dep_rows] <- -it$y[dep_rows]
  }
  .ipx_it_update(it, 1, dx, NULL, NULL, 1, dy, NULL, NULL)
  for (j in dep_cols) .ipx_it_make_fixed(it, j, 0)
  for (i in dep_rows) .ipx_it_make_implied(it, n + i, 8L)
}

# ---------------------------------------------------------------------------
# Maxvolume (maxvolume.cc), the heuristic variant with row slices

.ipx_maxvolume <- function(b, colscale) {
  M <- b$M; m <- M$m; n <- M$n; AI <- M$AI
  par <- b$par
  volumetol <- max(par$volume_tol, 1)
  maxskip <- par$maxskip_updates
  updates <- 0L; errflag <- 0L
  num_slices <- min(5L + max(m %/% par$rows_per_slice, 0L), m)
  inv <- numeric(m)
  pb <- b$map[b$basis] <= m
  inv[pb] <- 1 / colscale[b$basis[pb]]
  cs <- numeric(n + m)
  nbm <- b$map == -1L
  cs[nbm] <- colscale[nbm]
  perm <- .ipx_sortperm(inv)
  for (s in seq_len(num_slices) - 1L) {
    used <- logical(m)
    used[perm] <- (seq_len(m) - 1L) %% num_slices == s
    work <- .ipx_basis_solve(b, ifelse(used, inv, 0), "T")
    cw <- numeric(n + m)
    k <- cs != 0
    cw[k] <- .ipx_coldots(AI[, k, drop = FALSE], work) * cs[k]
    candidates <- integer()
    skipped <- 0L
    repeat {
      if (!length(candidates)) {
        aw <- abs(cw); jmax <- 1L; jmax2 <- 1L; wmax <- 0; wmax2 <- 0
        for (j in seq_along(aw)) {
          if (aw[j] > wmax) { wmax2 <- wmax; wmax <- aw[j]; jmax2 <- jmax; jmax <- j }
          else if (aw[j] > wmax2) { wmax2 <- aw[j]; jmax2 <- j }
        }
        candidates <- c(jmax2, jmax)
      }
      jn <- candidates[length(candidates)]
      if (cw[jn] == 0) break
      lhs <- .ipx_basis_sfu(b, jn)
      ord <- .ipx_iv_order(b$lu$pattern, m)
      # ScaleFtran
      sc <- lhs * cs[jn] * inv
      v <- ifelse(abs(lhs) > 1e-7, abs(sc), 0)
      pmax <- .ipx_first_max(v, ord)
      if (pmax == 0L) pmax <- 1L
      lhs[ord] <- sc[ord]
      vmax <- abs(lhs[pmax])
      if (vmax <= volumetol) {
        cw[jn] <- 0; cs[jn] <- 0
        candidates <- candidates[-length(candidates)]
        skipped <- skipped + 1L
        if (skipped > maxskip && maxskip >= 0) break
        next
      }
      wrec <- 0
      for (p in ord) if (used[p]) wrec <- wrec + lhs[p]
      jb <- b$basis[pmax]
      tr <- .ipx_basis_tableau_row(b, jb, TRUE)
      pivot <- tr$row[jn]
      ex <- .ipx_basis_exchange(b, jb, jn, pivot, 0L)
      errflag <- ex$err
      if (errflag != 0L) break
      if (!ex$exchanged) next
      updates <- updates + 1L
      dn <- cs[jn]; dbinv <- inv[pmax]
      cs[jb] <- 1 / inv[pmax]
      inv[pmax] <- 1 / cs[jn]
      cs[jn] <- 0
      alpha <- (used[pmax] - wrec) / (dn * pivot)
      ro <- .ipx_iv_order(tr$row_pattern, n + m)
      cw[ro] <- .ipx_madd(alpha * tr$row[ro], cs[ro], cw[ro])
      cw[jb] <- used[pmax] + alpha / dbinv
      cw[jn] <- 0
      candidates <- integer()
    }
    if (errflag != 0L) break
  }
  list(errflag = errflag, updates = updates)
}

# ---------------------------------------------------------------------------
# KKTSolverBasis (kkt_solver_basis.cc) with SplittedNormalMatrix

.ipx_kkt_basis <- function(M, b, par, info) {
  k <- new.env(parent = emptyenv())
  n <- M$n; m <- M$m; AI <- M$AI
  k$maxiter <- -1L
  k$iter <- 0L
  k$basis_changes <- 0L
  k$basis <- b
  drop_primal <- function(it) {
    xl <- it$xl; xu <- it$xu; zl <- it$zl; zu <- it$zu
    cand <- integer()
    for (p in seq_len(m)) {
      jb <- b$basis[p]
      if (b$map[jb] > m) next
      if (xl[jb] <= xu[jb]) { xj <- xl[jb]; zj <- zl[jb] } else { xj <- xu[jb]; zj <- zu[jb] }
      if (xj < 0.01 * zj && xj <= par$ipm_drop_primal) cand <- c(cand, jb)
    }
    if (!length(cand)) return(invisible())
    inv <- 1 / k$colscale[b$basis]
    while (length(cand)) {
      jb <- cand[length(cand)]
      p <- .ipx_b_pos(b, jb)
      s <- inv[p]
      tr <- .ipx_basis_tableau_row(b, jb, TRUE)
      ord <- .ipx_iv_order(tr$row_pattern, n + m)
      piv <- abs(tr$row)
      v <- ifelse(piv > 1e-7, piv * k$colscale * s, 0)
      v[is.na(v)] <- 0
      jmax <- .ipx_first_max(v, ord, 2)
      if (jmax > 0L) {
        ex <- .ipx_basis_exchange(b, jb, jmax, tr$row[jmax], 1L)
        info$errflag <- ex$err
        if (ex$err != 0L) return(invisible())
        if (!ex$exchanged) next
        inv[p] <- 1 / k$colscale[jmax]
        info$updates_ipm <- info$updates_ipm + 1L
        k$basis_changes <- k$basis_changes + 1L
      } else {
        if (zl[jb] / xl[jb] > zu[jb] / xu[jb]) .ipx_it_make_implied(it, jb, 6L)
        else .ipx_it_make_implied(it, jb, 7L)
        b$map[jb] <- b$map[jb] + m
        inv[p] <- 0
        k$colscale[jb] <- Inf
        info$primal_dropped <- info$primal_dropped + 1L
      }
      cand <- cand[-length(cand)]
    }
  }
  drop_dual <- function(it) {
    xl <- it$xl; xu <- it$xu; zl <- it$zl; zu <- it$zu
    cand <- integer()
    for (jn in which(b$map == -1L)) {
      if (zl[jn] >= zu[jn]) { xj <- xl[jn]; zj <- zl[jn] } else { xj <- xu[jn]; zj <- zu[jn] }
      if (zj < 0.01 * xj && zj <= par$ipm_drop_dual) cand <- c(cand, jn)
    }
    if (!length(cand)) return(invisible())
    inv <- 1 / k$colscale[b$basis]
    while (length(cand)) {
      jn <- cand[length(cand)]
      s <- k$colscale[jn]
      ftran <- .ipx_basis_sfu(b, jn)
      ord <- .ipx_iv_order(b$lu$pattern, m)
      piv <- abs(ftran)
      v <- ifelse(piv > 1e-7, piv * inv * s, 0)
      pmax <- .ipx_first_max(v, ord, 2)
      if (pmax > 0L) {
        jb <- b$basis[pmax]
        ex <- .ipx_basis_exchange(b, jb, jn, ftran[pmax], -1L)
        info$errflag <- ex$err
        if (ex$err != 0L) return(invisible())
        if (!ex$exchanged) next
        inv[pmax] <- 1 / k$colscale[jn]
        info$updates_ipm <- info$updates_ipm + 1L
        k$basis_changes <- k$basis_changes + 1L
      } else {
        .ipx_it_make_fixed(it, jn)
        b$map[jn] <- -2L
        k$colscale[jn] <- 0
        info$dual_dropped <- info$dual_dropped + 1L
      }
      cand <- cand[-length(cand)]
    }
  }
  # SplittedNormalMatrix::Prepare. The factors stay in pivot step
  # coordinates; N keeps the row order of AI, which is the order in which
  # the upstream loops visit the entries of its row permuted copy.
  prepare <- function() {
    f <- b$lu$get_factors()
    cs <- k$colscale
    Uc <- f$U
    for (q in seq_len(m)) {
      j <- b$basis[f$colperm[q]]
      if (b$map[j] <= m) Uc[[q]]$x <- Uc[[q]]$x * cs[j]
    }
    k$Lc <- f$L; k$Uc <- Uc
    k$colperm <- f$colperm
    rpi <- integer(m); rpi[f$rowperm] <- seq_len(m)
    k$rowperm_inv <- rpi
    nb <- which(b$map == -1L)
    k$nonbasic <- nb
    k$N <- AI[, nb, drop = FALSE] * rep(cs[nb], each = m)
    k$free_pos <- which(b$map[b$basis[f$colperm]] > m)
  }
  # SplittedNormalMatrix::_Apply: BackwardSolve, AddNormalProduct,
  # ForwardSolve (sparse_utils.cc) and the dot product with rhs
  apply_normal <- function(v) {
    x <- v
    Uc <- k$Uc; Lc <- k$Lc
    for (i in seq_len(m)) {
      cl <- Uc[[i]]; ne <- length(cl$i)
      d <- 0
      for (t in seq_len(ne - 1L)) d <- .ipx_madd(x[cl$i[t]], cl$x[t], d)
      x[i] <- (x[i] - d) / cl$x[ne]
    }
    for (i in rev(seq_len(m))) {
      cl <- Lc[[i]]
      d <- 0
      for (t in seq_along(cl$i)) d <- .ipx_madd(x[cl$i[t]], cl$x[t], d)
      x[i] <- x[i] - d
    }
    rpi <- k$rowperm_inv
    lo <- .ipx_scatter(k$N, .ipx_coldots(k$N, x[rpi]), numeric(m))
    x <- numeric(m); x[rpi] <- lo
    for (j in seq_len(m)) {
      temp <- x[j]
      if (temp != 0) {
        cl <- Lc[[j]]
        if (length(cl$i)) x[cl$i] <- .ipx_madd(-cl$x, temp, x[cl$i])
      }
    }
    for (j in rev(seq_len(m))) {
      cl <- Uc[[j]]; ne <- length(cl$i)
      x[j] <- x[j] / cl$x[ne]
      temp <- x[j]
      if (temp != 0 && ne > 1L) {
        r <- cl$i[-ne]
        x[r] <- .ipx_madd(-cl$x[-ne], temp, x[r])
      }
    }
    lhs <- x + v
    lhs[k$free_pos] <- 0
    list(lhs = lhs, dot = .ipx_dot(v, lhs))
  }
  k$factorize <- function(it) {
    info$errflag <- 0L
    k$iter <- 0L
    k$basis_changes <- 0L
    k$colscale <- .ipx_it_scaling(it)
    .ipx_it_evaluate(it)
    if (it$pobjective >= it$dobjective) {
      drop_primal(it)
      if (info$errflag) return(invisible())
      drop_dual(it)
      if (info$errflag) return(invisible())
    }
    mv <- .ipx_maxvolume(b, k$colscale)
    info$updates_ipm <- info$updates_ipm + mv$updates
    k$basis_changes <- k$basis_changes + mv$updates
    info$errflag <- mv$errflag
    if (info$errflag) return(invisible())
    if (!b$fresh) {
      info$errflag <- .ipx_basis_factorize(b)
      if (info$errflag) return(invisible())
    }
    prepare()
  }
  k$solve <- function(a, bb, tol) {
    info$errflag <- 0L
    cs <- k$colscale
    bas <- b$basis
    isfree <- b$map[bas] > m
    work <- numeric(m)
    nf <- sum(isfree)
    if (nf > 0) {
      work[isfree] <- a[bas[isfree]]
      work <- .ipx_basis_solve(b, work, "T")
    }
    nb <- k$nonbasic
    d2 <- cs[nb] * cs[nb]
    ANb <- AI[, nb, drop = FALSE]
    alpha <- if (nf > 0) (a[nb] - .ipx_coldots(ANb, work)) * d2 else d2 * a[nb]
    rhs <- .ipx_scatter(ANb, alpha, numeric(m))
    rhs <- .ipx_basis_solve(b, rhs, "N")
    work <- .ipx_basis_solve(b, bb, "N")
    d <- cs[bas]
    rhs <- ifelse(isfree, 0, .ipx_madd(a[bas], d, (rhs - work) / d))
    cr <- .ipx_cr(apply_normal, rhs[k$colperm], tol, NULL, k$maxiter)
    info$errflag <- cr$errflag
    info$kktiter2 <- info$kktiter2 + cr$iter
    k$iter <- k$iter + cr$iter
    y <- numeric(m)
    y[k$colperm] <- cr$lhs
    y <- ifelse(isfree, a[bas], y / d)
    y <- .ipx_basis_solve(b, y, "T")
    x <- numeric(n + m)
    x[nb] <- (a[nb] - .ipx_coldots(ANb, y)) * d2
    work <- .ipx_scatter(ANb, -x[nb], bb)
    x[bas] <- .ipx_basis_solve(b, work, "N")
    list(x = x, y = y)
  }
  k
}

# ---------------------------------------------------------------------------
# IPM (ipm.cc)

.ipx_step_to_boundary <- function(x, dx) {
  damp <- 1 - .Machine$double.eps
  alpha <- 1; iblock <- 0L
  for (i in which(dx < 0)) {
    if (.ipx_madd(alpha, dx[i], x[i]) < 0) {
      alpha <- -(x[i] * damp) / dx[i]
      iblock <- i
    }
  }
  list(alpha = alpha, block = iblock)
}

.ipx_ipm_new <- function(par, info) {
  ipm <- new.env(parent = emptyenv())
  ipm$par <- par; ipm$info <- info
  ipm$best_complementarity <- 0
  ipm$maxiter <- -1L
  ipm$step_primal <- 0; ipm$step_dual <- 0
  ipm
}

.ipx_ipm_log <- function(ipm) {
  it <- ipm$it; info <- ipm$info
  .ipx_it_evaluate(it)
  p <- .ipx_it_pobj_pp(it); d <- .ipx_it_dobj_pp(it)
  row <- c(iter = info$iter, pobj = p, dobj = d,
           pinf = it$presidual / it$bounds_measure,
           dinf = it$dresidual / it$costs_measure,
           gap = abs(p - d) / .ipx_madd(0.5, abs(p + d), 1),
           kktiter = ipm$kkt$iter, step_primal = ipm$step_primal,
           step_dual = ipm$step_dual)
  info$trace[[length(info$trace) + 1L]] <- row
  if (isTRUE(ipm$par$log))
    cat(sprintf(" %3d%s  %15.8e  %15.8e  %9.2e  %9.2e  %8.2e  %4.2f %4.2f  %7d %7d  %7d %7d\n",
                info$iter, if (.ipx_it_feasible(it) && .ipx_it_optimal(it)) "*" else " ",
                p, d, row[["pinf"]], row[["dinf"]], row[["gap"]], ipm$step_primal,
                ipm$step_dual, as.integer(ipm$kkt$basis_changes %||% 0L), ipm$kkt$iter,
                info$dual_dropped, info$primal_dropped))
}

.ipx_ipm_solve_newton <- function(ipm, rb, rc, rl, ru, sl, su) {
  it <- ipm$it; M <- it$M; AI <- M$AI
  xl <- it$xl; xu <- it$xu; zl <- it$zl; zu <- it$zu
  hl <- .ipx_has_lb(it); hu <- .ipx_has_ub(it)
  fixed <- .ipx_is_fixed(it); fr <- .ipx_is_freest(it); bar <- .ipx_is_barrier(it)
  rhs1 <- -rc
  rhs1[hl] <- rhs1[hl] + .ipx_madd(zl[hl], rl[hl], sl[hl]) / xl[hl]
  rhs1[hu] <- rhs1[hu] - .ipx_madd(-zu[hu], ru[hu], su[hu]) / xu[hu]
  rhs1[fixed] <- 0
  tol <- ipm$par$kkt_tol * sqrt(.ipx_it_get(it, "mu"))
  s <- ipm$kkt$solve(rhs1, rb, tol)
  if (ipm$info$errflag) {
    # The C code solves into the components of the reused Step object and
    # returns before recovering the bound steps: dx and dy (not negated)
    # are new, the bound steps are those of the last step taken. Only the
    # centring steps go on with such a step.
    N <- length(rhs1); ls <- ipm$last_step
    if (is.null(ls)) ls <- list(xl = numeric(N), xu = numeric(N), zl = numeric(N), zu = numeric(N))
    return(list(x = s$x, xl = ls$xl, xu = ls$xu, y = s$y, zl = ls$zl, zu = ls$zu, partial = TRUE))
  }
  dx <- s$x; dy <- -s$y
  N <- length(dx)
  dxl <- numeric(N); dzl <- numeric(N); dxu <- numeric(N); dzu <- numeric(N)
  dxl[bar] <- dx[bar] - rl[bar]
  dzl[bar] <- .ipx_madd(-zl[bar], dxl[bar], sl[bar]) / xl[bar]
  dxu[bar] <- ru[bar] - dx[bar]
  dzu[bar] <- .ipx_madd(-zu[bar], dxu[bar], su[bar]) / xu[bar]
  atdy <- .ipx_coldots(AI, dy)
  fl <- is.finite(xl); fu <- is.finite(xu)
  for (j in which(bar)) {
    if (fl[j] && fu[j]) {
      if (zl[j] * xu[j] >= zu[j] * xl[j]) dzl[j] <- rc[j] + dzu[j] - atdy[j]
      else dzu[j] <- -rc[j] + dzl[j] + atdy[j]
    } else if (fl[j]) {
      dzl[j] <- rc[j] + dzu[j] - atdy[j]
    } else {
      dzu[j] <- -rc[j] + dzl[j] + atdy[j]
    }
  }
  list(x = dx, xl = dxl, xu = dxu, y = dy, zl = dzl, zu = dzu)
}

.ipx_ipm_predictor <- function(ipm) {
  it <- ipm$it
  .ipx_it_evaluate(it)
  N <- length(it$x)
  hl <- .ipx_has_lb(it); hu <- .ipx_has_ub(it)
  sl <- numeric(N); su <- numeric(N)
  sl[hl] <- -it$xl[hl] * it$zl[hl]
  su[hu] <- -it$xu[hu] * it$zu[hu]
  .ipx_ipm_solve_newton(ipm, it$rb, it$rc, it$rl, it$ru, sl, su)
}

.ipx_ipm_corrector <- function(ipm, st) {
  it <- ipm$it
  .ipx_it_evaluate(it)
  xl <- it$xl; xu <- it$xu; zl <- it$zl; zu <- it$zu
  mu <- it$mu
  maxp <- min(.ipx_step_to_boundary(xl, st$xl)$alpha, .ipx_step_to_boundary(xu, st$xu)$alpha)
  maxd <- min(.ipx_step_to_boundary(zl, st$zl)$alpha, .ipx_step_to_boundary(zu, st$zu)$alpha)
  hl <- .ipx_has_lb(it); hu <- .ipx_has_ub(it)
  muaff <- .ipx_mu_after_step(xl, xu, zl, zu, st, maxp, maxd, hl, hu)
  muaff <- muaff / (sum(hl) + sum(hu))
  ratio <- muaff / mu
  sigma <- ratio * ratio * ratio
  N <- length(xl)
  sl <- numeric(N); su <- numeric(N)
  sm <- sigma * mu
  sl[hl] <- .ipx_madd(-st$xl[hl], st$zl[hl], .ipx_madd(-xl[hl], zl[hl], sm))
  su[hu] <- .ipx_madd(-st$xu[hu], st$zu[hu], .ipx_madd(-xu[hu], zu[hu], sm))
  .ipx_ipm_solve_newton(ipm, it$rb, it$rc, it$rl, it$ru, sl, su)
}

# sum over j of the complementarity products after the step, lower bound
# term before upper bound term for each j
.ipx_mu_after_step <- function(xl, xu, zl, zu, st, maxp, maxd, hl, hu) {
  pl <- .ipx_madd(maxp, st$xl, xl); ql <- .ipx_madd(maxd, st$zl, zl)
  pu <- .ipx_madd(maxp, st$xu, xu); qu <- .ipx_madd(maxd, st$zu, zu)
  mu <- 0
  for (j in which(hl | hu)) {
    if (hl[j]) mu <- .ipx_madd(pl[j], ql[j], mu)
    if (hu[j]) mu <- .ipx_madd(pu[j], qu[j], mu)
  }
  mu
}

.ipx_ipm_step_sizes <- function(ipm, st, centring = FALSE) {
  it <- ipm$it
  xl <- it$xl; xu <- it$xu; zl <- it$zl; zu <- it$zu
  gammaf <- 0.9; gammaa <- 1 / (1 - gammaf)
  sxl <- .ipx_step_to_boundary(xl, st$xl); sxu <- .ipx_step_to_boundary(xu, st$xu)
  szl <- .ipx_step_to_boundary(zl, st$zl); szu <- .ipx_step_to_boundary(zu, st$zu)
  maxp <- min(sxl$alpha, sxu$alpha); maxd <- min(szl$alpha, szu$alpha)
  hl <- .ipx_has_lb(it); hu <- .ipx_has_ub(it)
  mufull <- .ipx_mu_after_step(xl, xu, zl, zu, st, maxp, maxd, hl, hu)
  mufull <- mufull / (sum(hl) + sum(hu)) / gammaa
  alphap <- 1; alphad <- 1
  if (maxp < 1) {
    if (sxl$alpha <= sxu$alpha) {
      j <- sxl$block
      alphap <- (xl[j] - mufull / .ipx_madd(maxd, st$zl[j], zl[j])) / (-st$xl[j])
    } else {
      j <- sxu$block
      alphap <- (xu[j] - mufull / .ipx_madd(maxd, st$zu[j], zu[j])) / (-st$xu[j])
    }
    alphap <- min(max(alphap, gammaf * maxp), 1)
  }
  if (maxd < 1) {
    if (szl$alpha <= szu$alpha) {
      j <- szl$block
      alphad <- (zl[j] - mufull / .ipx_madd(maxp, st$xl[j], xl[j])) / (-st$zl[j])
    } else {
      j <- szu$block
      alphad <- (zu[j] - mufull / .ipx_madd(maxp, st$xu[j], xu[j])) / (-st$zu[j])
    }
    alphad <- min(max(alphad, gammaf * maxd), 1)
  }
  ipm$step_primal <- min(alphap, 1 - 1e-6)
  ipm$step_dual <- min(alphad, 1 - 1e-6)
  if (centring) {
    # centring steps are scaled by centringAlphaScaling
    ipm$step_primal <- alphap * ipm$par$centring_alpha_scaling
    ipm$step_dual <- alphad * ipm$par$centring_alpha_scaling
  }
}

.ipx_ipm_make_step <- function(ipm, st, centring = FALSE) {
  .ipx_ipm_step_sizes(ipm, st, centring)
  it <- ipm$it
  .ipx_it_update(it, ipm$step_primal, st$x, st$xl, st$xu, ipm$step_dual, st$y, st$zl, st$zu)
  ipm$last_step <- st
  if (centring) return(invisible())
  if (min(ipm$step_primal, ipm$step_dual) < 0.05) ipm$num_bad_iter <- ipm$num_bad_iter + 1L
  else ipm$num_bad_iter <- 0L
  ipm$best_complementarity <- min(ipm$best_complementarity, .ipx_it_get(it, "complementarity"))
}

.ipx_ipm_status_from_err <- function(info) {
  if (info$errflag == .ipx_err$time_interrupt) {
    info$errflag <- 0L; info$status_ipm <- .ipx_st$time_limit
  } else {
    info$status_ipm <- .ipx_st$failed
  }
}

# IPM::StartingPoint
.ipx_ipm_starting_point <- function(ipm, kkt, it) {
  ipm$kkt <- kkt; ipm$it <- it
  info <- ipm$info; M <- it$M; n <- M$n; m <- M$m; AI <- M$AI
  lb <- M$lb; ub <- M$ub; cc <- M$c
  N <- n + m
  kkt$factorize(NULL)
  x <- pmin(pmax(numeric(N), lb), ub)
  rb <- .ipx_scatter(AI, -x, M$b)
  s <- kkt$solve(numeric(N), rb, 0.1 * .ipx_infnorm(rb))
  if (info$errflag) return(.ipx_ipm_status_from_err(info))
  x <- x + s$x
  y <- s$y
  xl <- x - lb; xu <- ub - x
  xinfeas <- max(0, -xl, -xu)
  xshift1 <- .ipx_madd(1.5, xinfeas, 1)
  xl <- xl + xshift1; xu <- xu + xshift1
  fl <- is.finite(lb); fu <- is.finite(ub)
  cnorm <- sqrt(.ipx_dot(cc, cc))
  if (cnorm == 0) {
    # zero objective: y stays the solution of the first KKT solve
    zl <- ifelse(fl, 1, 0); zu <- ifelse(fu, 1, 0)
  } else {
    s <- kkt$solve(cc, numeric(m), 0.1 * .ipx_infnorm(cc))
    if (info$errflag) return(.ipx_ipm_status_from_err(info))
    y <- s$y
    zl <- cc - .ipx_coldots(AI, y)
    if (sqrt(.ipx_dot(zl, zl)) < 0.05 * cnorm) {
      zl <- zl + 0.05 * cc
      y <- y * (1 - 0.05)
    }
    zval <- zl
    zl <- numeric(N); zu <- numeric(N)
    bx <- fl & fu
    zl[bx] <- 0.5 * zval[bx]; zu[bx] <- -0.5 * zval[bx]
    zl[fl & !fu] <- zval[fl & !fu]
    zu[fu & !fl] <- -zval[fu & !fl]
    zinfeas <- max(0, -zl, -zu)
    zshift1 <- .ipx_madd(1.5, zinfeas, 1)
    zl[fl] <- zl[fl] + zshift1
    zu[fu] <- zu[fu] + zshift1
  }
  xsum <- 1; zsum <- 1; mu <- 1
  for (j in which(fl | fu)) {
    if (fl[j]) { xsum <- xsum + xl[j]; zsum <- zsum + zl[j]; mu <- .ipx_madd(xl[j], zl[j], mu) }
    if (fu[j]) { xsum <- xsum + xu[j]; zsum <- zsum + zu[j]; mu <- .ipx_madd(xu[j], zu[j], mu) }
  }
  xl <- xl + 0.5 * mu / zsum; xu <- xu + 0.5 * mu / zsum
  zl[fl] <- zl[fl] + 0.5 * mu / xsum
  zu[fu] <- zu[fu] + 0.5 * mu / xsum
  .ipx_it_initialize(it, x, xl, xu, y, zl, zu)
  ipm$best_complementarity <- .ipx_it_get(it, "complementarity")
  .ipx_ipm_log(ipm)
  info$status_ipm <- .ipx_st$not_run
}

# IPM::Driver
.ipx_ipm_driver <- function(ipm, kkt, it) {
  ipm$kkt <- kkt; ipm$it <- it
  info <- ipm$info
  ipm$num_bad_iter <- 0L
  ipm$last_step <- NULL
  repeat {
    if (.ipx_it_term_crit(it)) { info$status_ipm <- .ipx_st$optimal; break }
    if (ipm$num_bad_iter >= 5L ||
        (ipm$best_complementarity > 0 &&
         .ipx_it_get(it, "complementarity") > 1e6 * ipm$best_complementarity)) {
      p <- .ipx_it_pobj_pp(it); d <- .ipx_it_dobj_pp(it)
      if (d > max(10 * abs(p), 1)) {
        info$status_ipm <- if (it$M$dualized) .ipx_st$dual_infeas else .ipx_st$primal_infeas
      } else if (p < -max(10 * abs(d), 1)) {
        info$status_ipm <- if (it$M$dualized) .ipx_st$primal_infeas else .ipx_st$dual_infeas
      } else {
        info$status_ipm <- .ipx_st$no_progress
      }
      break
    }
    if (info$iter >= ipm$maxiter) { info$status_ipm <- .ipx_st$iter_limit; break }
    if (is.finite(ipm$par$time_limit) && ipm$par$time_limit >= 0 &&
        proc.time()[[3]] - ipm$par$t0 > ipm$par$time_limit) {
      info$errflag <- .ipx_err$time_interrupt; break
    }
    kkt$factorize(it)
    if (info$errflag) break
    st <- .ipx_ipm_predictor(ipm)
    if (info$errflag) break
    st <- .ipx_ipm_corrector(ipm, st)
    if (info$errflag) break
    .ipx_ipm_make_step(ipm, st)
    info$iter <- info$iter + 1L
    .ipx_ipm_log(ipm)
  }
  if (info$errflag) .ipx_ipm_status_from_err(info)
  if (ipm$par$run_centring && info$status_ipm == .ipx_st$optimal && !info$centring_tried)
    .ipx_ipm_centring(ipm)
}

# IPM::AssessCentrality: the ratio max(xj * zj) / min(xj * zj) over the
# barrier terms (lower bounds first), with mu inside the range, and the
# number of products outside [0.1 * mu, 10 * mu]
.ipx_assess_centrality <- function(it, xl, xu, zl, zu, mu) {
  hl <- .ipx_has_lb(it); hu <- .ipx_has_ub(it)
  p <- c((xl * zl)[hl], (xu * zu)[hu])
  gamma <- 0.1
  bad <- sum(p < gamma * mu | p > mu / gamma)
  list(ratio = max(p, 0, mu) / min(p, Inf, mu), bad = bad)
}

# the centring steps at the end of IPM::Driver (run_centring): Newton steps
# towards the central point of the frozen mu, accepted while the centrality
# improves
.ipx_ipm_centring <- function(ipm) {
  it <- ipm$it; info <- ipm$info; par <- ipm$par
  mu_frozen <- .ipx_it_get(it, "mu")
  ac <- .ipx_assess_centrality(it, it$xl, it$xu, it$zl, it$zu, mu_frozen)
  info$centring_success <- FALSE
  if (ac$ratio < par$centring_ratio_tolerance) {
    info$centring_success <- TRUE
  } else {
    for (ii in seq_len(par$max_centring_steps)) {
      # IPM::Centring: sigma = 1
      .ipx_it_evaluate(it)
      N <- length(it$x)
      hl <- .ipx_has_lb(it); hu <- .ipx_has_ub(it)
      sl <- numeric(N); su <- numeric(N)
      sl[hl] <- -it$xl[hl] * it$zl[hl] + mu_frozen
      su[hu] <- -it$xu[hu] * it$zu[hu] + mu_frozen
      st <- .ipx_ipm_solve_newton(ipm, it$rb, it$rc, it$rl, it$ru, sl, su)
      # IPM::EvaluateCentringStep
      .ipx_ipm_step_sizes(ipm, st, TRUE)
      xl <- it$xl; xu <- it$xu; zl <- it$zl; zu <- it$zu
      xl[hl] <- xl[hl] + ipm$step_primal * st$xl[hl]
      xu[hu] <- xu[hu] + ipm$step_primal * st$xu[hu]
      zl[hl] <- zl[hl] + ipm$step_dual * st$zl[hl]
      zu[hu] <- zu[hu] + ipm$step_dual * st$zu[hu]
      mu_temp <- 0
      for (j in which(hl | hu)) {
        if (hl[j]) mu_temp <- mu_temp + xl[j] * zl[j]
        if (hu[j]) mu_temp <- mu_temp + xu[j] * zu[j]
      }
      mu_temp <- mu_temp / (sum(hl) + sum(hu))
      tc <- .ipx_assess_centrality(it, xl, xu, zl, zu, mu_temp)
      if (!(tc$ratio < par$centring_ratio_reduction * ac$ratio || tc$bad < ac$bad)) break
      .ipx_ipm_make_step(ipm, st, TRUE)
      info$iter <- info$iter + 1L
      .ipx_ipm_log(ipm)
      ac <- .ipx_assess_centrality(it, it$xl, it$xu, it$zl, it$zu, .ipx_it_get(it, "mu"))
      if (ac$ratio < par$centring_ratio_tolerance) {
        info$centring_success <- TRUE
        break
      }
    }
  }
  info$centring_tried <- TRUE
}

# ---------------------------------------------------------------------------
# Crossover (crossover.cc)

.ipx_push_dual <- function(b, y, z, vars, x, par, info) {
  M <- b$M; lb <- M$lb; ub <- M$ub
  feastol <- if (M$dualized) par$pfeasibility_tol else par$dfeasibility_tol
  sr1 <- x != ub; sr2 <- x != lb
  pivots <- 0L
  nxt <- 1L
  while (nxt <= length(vars)) {
    jb <- vars[nxt]
    if (z[jb] == 0) { nxt <- nxt + 1L; next }
    tr <- .ipx_basis_tableau_row(b, jb)
    row <- tr$row
    ord <- .ipx_iv_order(tr$row_pattern, length(row))
    step <- z[jb]
    # DualRatioTest
    jn <- 0L
    for (j in ord[abs(row[ord]) > 1e-5]) {
      piv <- row[j]
      if (sr1[j] && .ipx_madd(-step, piv, z[j]) < -feastol) { step <- (z[j] + feastol) / piv; jn <- j }
      if (sr2[j] && .ipx_madd(-step, piv, z[j]) > feastol) { step <- (z[j] - feastol) / piv; jn <- j }
    }
    if (jn > 0L) {
      jn <- 0L; maxp <- 1e-5
      for (j in ord) {
        piv <- row[j]
        if (abs(piv) > maxp && abs(z[j] / piv) <= abs(step)) {
          if (sr1[j] && step * piv > 0) { jn <- j; maxp <- abs(piv) }
          if (sr2[j] && step * piv < 0) { jn <- j; maxp <- abs(piv) }
        }
      }
    }
    if (jn > 0L) {
      ex <- .ipx_basis_exchange(b, jb, jn, row[jn], 1L)
      info$errflag <- ex$err
      if (ex$err != 0L) break
      if (!ex$exchanged) next
      pivots <- pivots + 1L
      step <- z[jn] / row[jn]
    }
    if (step != 0) {
      bo <- .ipx_iv_order(tr$btran_pattern, length(y))
      y[bo] <- .ipx_madd(step, tr$btran[bo], y[bo])
      z[ord] <- .ipx_madd(-step, row[ord], z[ord])
      k <- ord[sr1[ord]]; z[k] <- pmax(z[k], 0)
      k <- ord[sr2[ord]]; z[k] <- pmin(z[k], 0)
      z[jb] <- z[jb] - step
    }
    if (jn > 0L) z[jn] <- 0
    nxt <- nxt + 1L
  }
  info$status_crossover <- if (info$errflag != 0L) .ipx_st$failed else .ipx_st$optimal
  list(y = y, z = z, pivots = pivots)
}

.ipx_push_primal <- function(b, x, vars, par, info) {
  M <- b$M; lb <- M$lb; ub <- M$ub
  feastol <- if (M$dualized) par$dfeasibility_tol else par$pfeasibility_tol
  xb <- x[b$basis]; lbb <- lb[b$basis]; ubb <- ub[b$basis]
  pivots <- 0L
  nxt <- 1L
  while (nxt <= length(vars)) {
    jn <- vars[nxt]
    if (x[jn] == lb[jn] || x[jn] == ub[jn] ||
        (x[jn] == 0 && is.infinite(lb[jn]) && is.infinite(ub[jn]))) { nxt <- nxt + 1L; next }
    move_to <- 0
    if (is.finite(lb[jn]) && is.finite(ub[jn]))
      move_to <- if (x[jn] - lb[jn] <= ub[jn] - x[jn]) lb[jn] else ub[jn]
    else if (is.finite(lb[jn])) move_to <- lb[jn]
    else if (is.finite(ub[jn])) move_to <- ub[jn]
    step <- x[jn] - move_to
    ftran <- .ipx_basis_sfu(b, jn)
    ord <- .ipx_iv_order(b$lu$pattern, length(ftran))
    # PrimalRatioTest
    pblock <- 0L; at_lb <- TRUE
    for (p in ord[abs(ftran[ord]) > 1e-5]) {
      piv <- ftran[p]
      if (.ipx_madd(step, piv, xb[p]) < lbb[p] - feastol) {
        step <- (lbb[p] - xb[p] - feastol) / piv; pblock <- p; at_lb <- TRUE
      }
      if (.ipx_madd(step, piv, xb[p]) > ubb[p] + feastol) {
        step <- (ubb[p] - xb[p] + feastol) / piv; pblock <- p; at_lb <- FALSE
      }
    }
    if (pblock > 0L) {
      pblock <- 0L; maxp <- 1e-5
      for (p in ord) {
        piv <- ftran[p]
        if (abs(piv) > maxp) {
          if (step * piv < 0 && abs((lbb[p] - xb[p]) / piv) <= abs(step)) {
            pblock <- p; at_lb <- TRUE; maxp <- abs(piv)
          }
          if (step * piv > 0 && abs((ubb[p] - xb[p]) / piv) <= abs(step)) {
            pblock <- p; at_lb <- FALSE; maxp <- abs(piv)
          }
        }
      }
    }
    jb <- if (pblock > 0L) b$basis[pblock] else 0L
    if (pblock > 0L) {
      ex <- .ipx_basis_exchange(b, jb, jn, ftran[pblock], -1L)
      info$errflag <- ex$err
      if (ex$err != 0L) break
      if (!ex$exchanged) next
      pivots <- pivots + 1L
      step <- if (at_lb) (lbb[pblock] - xb[pblock]) / ftran[pblock]
              else (ubb[pblock] - xb[pblock]) / ftran[pblock]
    }
    if (step != 0) {
      xb[ord] <- pmin(pmax(.ipx_madd(step, ftran[ord], xb[ord]), lbb[ord]), ubb[ord])
      x[jn] <- x[jn] - step
    }
    if (pblock > 0L) {
      x[jb] <- if (at_lb) lbb[pblock] else ubb[pblock]
      xb[pblock] <- x[jn]; lbb[pblock] <- lb[jn]; ubb[pblock] <- ub[jn]
    } else {
      x[jn] <- move_to
    }
    nxt <- nxt + 1L
  }
  x[b$basis] <- xb
  info$status_crossover <- if (info$errflag != 0L) .ipx_st$failed else .ipx_st$optimal
  list(x = x, pivots = pivots)
}

.ipx_crossover <- function(b, x, y, z, weights, par, info) {
  M <- b$M; lb <- M$lb; ub <- M$ub
  perm <- .ipx_sortperm(weights)
  dvars <- perm[b$map[perm] >= 1L & z[perm] != 0]
  pd <- .ipx_push_dual(b, y, z, dvars, x, par, info)
  y <- pd$y; z <- pd$z
  info$updates_crossover <- pd$pivots
  if (info$status_crossover != .ipx_st$optimal) return(list(x = x, y = y, z = z))
  rp <- rev(perm)
  pvars <- rp[b$map[rp] < 0L & x[rp] != lb[rp] & x[rp] != ub[rp] &
                !(is.infinite(lb[rp]) & is.infinite(ub[rp]) & x[rp] == 0)]
  pp <- .ipx_push_primal(b, x, pvars, par, info)
  info$updates_crossover <- pd$pivots + pp$pivots
  list(x = pp$x, y = y, z = z)
}

# ---------------------------------------------------------------------------
# LpSolver (lp_solver.cc)

.ipx_lp_solve <- function(M, par) {
  info <- new.env(parent = emptyenv())
  for (f in c("iter", "errflag", "kktiter1", "kktiter2", "updates_start", "updates_ipm",
              "updates_crossover", "dependent_rows", "dependent_cols", "basis_repairs",
              "primal_dropped", "dual_dropped", "status_ipm", "status_crossover"))
    assign(f, 0L, envir = info)
  info$rows_inconsistent <- FALSE; info$cols_inconsistent <- FALSE
  info$trace <- list()
  # InteriorPointSolve
  it <- .ipx_iterate_new(M)
  it$feasibility_tol <- par$ipm_feasibility_tol
  it$optimality_tol <- par$ipm_optimality_tol
  if (par$run_crossover != 0) it$start_crossover_tol <- par$start_crossover_tol
  ipm <- .ipx_ipm_new(par, info)
  b <- NULL
  info$centring_tried <- FALSE
  info$centring_success <- FALSE
  run_ipm <- function() {
    .ipx_ipm_starting_point(ipm, .ipx_kkt_diag(M, info), it)
    if (info$status_ipm != .ipx_st$not_run) return(invisible())
    # RunInitialIPM
    kkt <- .ipx_kkt_diag(M, info, if (par$switchiter < 0) min(500L, 10L + M$m %/% 20L) else -1L)
    ipm$maxiter <- if (par$switchiter < 0) par$ipm_maxiter else min(par$switchiter, par$ipm_maxiter)
    .ipx_ipm_driver(ipm, kkt, it)
    s <- info$status_ipm
    if (s == .ipx_st$optimal || s == .ipx_st$no_progress) info$status_ipm <- .ipx_st$not_run
    if (s == .ipx_st$failed) { info$status_ipm <- .ipx_st$not_run; info$errflag <- 0L }
    if (s == .ipx_st$iter_limit && info$iter < par$ipm_maxiter) info$status_ipm <- .ipx_st$not_run
    if (info$status_ipm != .ipx_st$not_run) return(invisible())
    # BuildStartingBasis
    b <<- .ipx_basis_new(M, par, info)
    .ipx_starting_basis(it, b, info)
    if (info$errflag) {
      if (info$errflag == .ipx_err$time_interrupt) {
        info$errflag <- 0L; info$status_ipm <- .ipx_st$time_limit
      } else info$status_ipm <- .ipx_st$failed
      return(invisible())
    }
    if (M$dualized) {
      t <- info$dependent_rows; info$dependent_rows <- info$dependent_cols; info$dependent_cols <- t
      t <- info$rows_inconsistent; info$rows_inconsistent <- info$cols_inconsistent
      info$cols_inconsistent <- t
    }
    if (info$rows_inconsistent) { info$status_ipm <- .ipx_st$primal_infeas; return(invisible()) }
    if (info$cols_inconsistent) { info$status_ipm <- .ipx_st$dual_infeas; return(invisible()) }
    if (info$centring_tried) return(invisible())
    # RunMainIPM
    ipm$maxiter <- par$ipm_maxiter
    .ipx_ipm_driver(ipm, .ipx_kkt_basis(M, b, par, info), it)
  }
  run_ipm()
  .ipx_it_postprocess(it)
  ev <- .ipx_evaluate_interior(M, it)
  if (info$status_ipm == .ipx_st$optimal &&
      (abs(ev$rel_objgap) > par$ipm_optimality_tol ||
       ev$rel_presidual > par$ipm_feasibility_tol || ev$rel_dresidual > par$ipm_feasibility_tol))
    info$status_ipm <- .ipx_st$imprecise
  if (info$centring_tried)
    info$status_ipm <- if (info$centring_success) .ipx_st$optimal else .ipx_st$imprecise
  run_co <- (info$status_ipm == .ipx_st$optimal && par$run_crossover == 1) ||
    (info$status_ipm == .ipx_st$imprecise && par$run_crossover != 0)
  basic <- NULL
  if (run_co) {
    st <- .ipx_it_drop_to_complementarity(it)
    w <- .ipx_it_scaling(it)
    co <- .ipx_crossover(b, st$x, st$y, st$z, w, par, info)
    if (info$status_crossover == .ipx_st$optimal) {
      bs <- .ipx_basis_basic_solution(b, co$x, co$z)
      N <- M$n + M$m
      status <- integer(N)
      for (j in seq_len(N)) {
        if (b$map[j] >= 1L) status[j] <- 0L
        else if (M$lb[j] == M$ub[j]) status[j] <- if (bs$z[j] >= 0) -1L else -2L
        else if (bs$x[j] == M$lb[j]) status[j] <- -1L
        else if (bs$x[j] == M$ub[j]) status[j] <- -2L
        else status[j] <- -3L
      }
      basic <- list(x = bs$x, y = bs$y, z = bs$z, status = status)
      eb <- .ipx_evaluate_basic(M, bs$x, bs$y, bs$z, status)
      if (eb$primal_infeas > par$pfeasibility_tol || eb$dual_infeas > par$dfeasibility_tol)
        info$status_crossover <- .ipx_st$imprecise
    }
  }
  if (info$status_ipm %in% c(.ipx_st$primal_infeas, .ipx_st$dual_infeas)) {
    info$status <- .ipx_st$solved
  } else {
    ms <- if (run_co) info$status_crossover else info$status_ipm
    info$status <- if (ms %in% c(.ipx_st$optimal, .ipx_st$imprecise)) .ipx_st$solved else .ipx_st$stopped
  }
  list(info = info, it = it, basic = basic)
}

# ---------------------------------------------------------------------------
# HiGHS glue (IpxWrapper.cpp, HighsSolution.cpp)

# fillInIpxData: boxed rows get an explicit slack column, free rows are
# dropped. Upstream indexes the slack entries by the original row and reads
# past the end of free_rows; both are fixed here (see PORTING.md).
.ipx_fill_data <- function(lp) {
  A <- as.matrix(lp$A)
  rl <- lp$row_lower; ru <- lp$row_upper
  nr <- length(rl)
  boxed_rows <- which(rl < ru & is.finite(rl) & is.finite(ru))
  free_rows <- which(is.infinite(rl) & rl < 0 & is.infinite(ru) & ru > 0)
  keep <- setdiff(seq_len(nr), free_rows)
  ctype <- character(length(keep)); rhs <- numeric(length(keep))
  for (t in seq_along(keep)) {
    r <- keep[t]
    if (is.finite(rl[r]) && !is.finite(ru[r])) { rhs[t] <- rl[r]; ctype[t] <- ">" }
    else if (!is.finite(rl[r]) && is.finite(ru[r])) { rhs[t] <- ru[r]; ctype[t] <- "<" }
    else if (rl[r] == ru[r]) { rhs[t] <- ru[r]; ctype[t] <- "=" }
    else { rhs[t] <- 0; ctype[t] <- "=" }
  }
  nb <- length(boxed_rows)
  S <- matrix(0, length(keep), nb)
  if (nb) S[cbind(match(boxed_rows, keep), seq_len(nb))] <- -1
  sense <- if (is.null(lp$sense)) 1 else lp$sense
  list(A = cbind(A[keep, , drop = FALSE], S), rhs = rhs, ctype = ctype,
       obj = c(sense * lp$c, numeric(nb)),
       lb = c(lp$col_lower, rl[boxed_rows]), ub = c(lp$col_upper, ru[boxed_rows]),
       offset = sense * (if (is.null(lp$offset)) 0 else lp$offset),
       boxed_rows = boxed_rows, free_rows = free_rows, keep = keep, sense = sense)
}

.ipx_par <- function(options) {
  o <- function(name, default) if (is.null(options[[name]])) default else options[[name]]
  pft <- o("primal_feasibility_tolerance", 1e-7)
  dft <- o("dual_feasibility_tolerance", 1e-7)
  rc <- o("run_crossover", "on")
  centring <- isTRUE(o("run_centring", FALSE))
  if (centring) rc <- "off"
  par <- list(ipm_feasibility_tol = min(pft, dft),
              ipm_optimality_tol = o("ipm_optimality_tolerance", 1e-8),
              start_crossover_tol = o("start_crossover_tolerance", 1e-8),
              ipm_maxiter = o("ipm_iteration_limit", .Machine$integer.max),
              time_limit = o("time_limit", Inf),
              run_crossover = switch(rc, on = 1, off = 0, choose = -1),
              dualize = switch(as.character(o("ipx_dualize_strategy", 2L)),
                               "-1" = 0, "1" = 1, "2" = -1, "3" = -2, -1),
              scale = 1, ipm_drop_primal = 1e-9, ipm_drop_dual = 1e-9, kkt_tol = 0.3,
              crash_basis = 1, dependency_tol = 1e-6, volume_tol = 2,
              rows_per_slice = 10000L, maxskip_updates = 10L, lu_pivottol = 0.0625,
              pfeasibility_tol = 1e-7, dfeasibility_tol = 1e-7, switchiter = -1L,
              update_heuristic = 1L, run_centring = centring,
              max_centring_steps = o("max_centring_steps", 5L),
              centring_ratio_tolerance = o("centring_ratio_tolerance", 100),
              centring_ratio_reduction = 1.5, centring_alpha_scaling = 0.5, log = isTRUE(options$log), t0 = proc.time()[[3]])
  kkt <- o("kkt_tolerance", 1e-7)
  if (kkt != 1e-7) {
    par$ipm_feasibility_tol <- kkt
    par$ipm_optimality_tol <- 0.1 * kkt
    par$start_crossover_tol <- 0.1 * kkt
  }
  if (par$run_crossover == 0) par$start_crossover_tol <- -1
  par
}

.highs_ipx_solve <- function(lp, options = list()) {
  par <- .ipx_par(options)
  d <- .ipx_fill_data(lp)
  nc <- length(lp$c); nr <- length(lp$row_lower)
  out <- list(status = "kSolveError", x = NULL, row_activity = NULL, col_dual = NULL,
              row_dual = NULL, col_basis = NULL, row_basis = NULL,
              ipm_iterations = 0L, crossover_iterations = 0L, objective = NULL)
  M <- .ipx_model_load(d$A, d$rhs, d$ctype, d$offset, d$obj, d$lb, d$ub, par)
  if (M$errflag != 0L) return(out)
  r <- .ipx_lp_solve(M, par)
  info <- r$info
  out$ipm_iterations <- info$iter
  out$crossover_iterations <- info$updates_crossover
  out$log <- do.call(rbind, info$trace)
  st <- .ipx_st
  sense <- d$sense
  A <- as.matrix(lp$A)
  interior <- function() {
    s <- .ipx_scale_back_interior(M, .ipx_dualize_back_interior(
      M, r$it$x, r$it$xl, r$it$xu, r$it$y, r$it$zl, r$it$zu))
    x <- s$x[seq_len(nc)]
    cd <- (s$zl - s$zu)[seq_len(nc)]
    rv <- .ipx_scatter(A, x, numeric(nr)); rd <- numeric(nr)
    for (t in seq_along(d$keep)) {
      row <- d$keep[t]
      k <- match(row, d$boxed_rows)
      if (!is.na(k)) {
        rv[row] <- s$x[nc + k]; rd[row] <- s$zl[nc + k] - s$zu[nc + k]
      } else {
        rv[row] <- d$rhs[t] - s$slack[t]; rd[row] <- s$y[t]
      }
    }
    rd[d$free_rows] <- 0
    out$x <<- x; out$col_dual <<- sense * cd
    out$row_activity <<- rv; out$row_dual <<- sense * rd
  }
  fill_objective <- function() {
    # HighsLp::objectiveValue: accumulated in order from the offset
    v <- if (is.null(lp$offset)) 0 else lp$offset
    for (j in seq_along(lp$c)) v <- .ipx_madd(lp$c[j], out$x[j], v)
    out$objective <<- v
  }
  bad <- c(st$failed, st$debug)
  if (info$status == st$stopped) {
    interior(); fill_objective()
    if (info$status_crossover %in% bad || info$status_ipm %in% bad) {
      out$status <- "kSolveError"
    } else if (info$status_crossover == st$time_limit) {
      out$status <- "kTimeLimit"
    } else if (info$status_ipm == st$time_limit) {
      out$status <- "kTimeLimit"
    } else if (info$status_ipm == st$iter_limit) {
      out$status <- "kIterationLimit"
    } else if (info$status_ipm == st$no_progress) {
      out$status <- "kUnknown"
    } else {
      out$status <- "kSolveError"
    }
    return(out)
  }
  if (info$status_crossover %in% bad) { out$status <- "kSolveError"; return(out) }
  if (info$status_ipm == st$primal_infeas || info$status_ipm == st$dual_infeas) {
    interior(); fill_objective()
    out$status <- if (info$status_ipm == st$primal_infeas) "kInfeasible" else "kUnboundedOrInfeasible"
    return(out)
  }
  imprecise <- info$status_crossover == st$imprecise || info$status_ipm == st$imprecise
  if (!is.null(r$basic)) {
    bs <- r$basic
    s <- .ipx_scale_back_basic(M, .ipx_basic_back(M, bs$x, bs$y, bs$z, bs$status))
    s$vbasis[M$flipped][s$vbasis[M$flipped] == -1L] <- -2L
    names_st <- c("0" = "kBasic", "-1" = "kLower", "-2" = "kUpper", "-3" = "kZero")
    cb <- unname(names_st[as.character(s$vbasis[seq_len(nc)])])
    x <- s$x[seq_len(nc)]
    cd <- ifelse(s$vbasis[seq_len(nc)] == 0L, 0, s$z[seq_len(nc)])
    rv <- .ipx_scatter(A, x, numeric(nr)); rd <- numeric(nr); rb <- rep("kBasic", nr)
    for (t in seq_along(d$keep)) {
      row <- d$keep[t]
      k <- match(row, d$boxed_rows)
      if (!is.na(k)) {
        rv[row] <- s$x[nc + k]
        sv <- s$vbasis[nc + k]
        if (s$cbasis[t] == 0L || sv == 0L) { rb[row] <- "kBasic"; rd[row] <- 0 }
        else if (sv == -1L) { rb[row] <- "kLower"; rd[row] <- s$z[nc + k] }
        else { rb[row] <- "kUpper"; rd[row] <- s$z[nc + k] }
      } else if (s$cbasis[t] == 0L) {
        rv[row] <- d$rhs[t] - s$slack[t]; rd[row] <- 0; rb[row] <- "kBasic"
      } else {
        rv[row] <- d$rhs[t] - s$slack[t]; rd[row] <- s$y[t]
        rb[row] <- switch(d$ctype[t], ">" = "kLower", "<" = "kUpper",
                          "=" = if (rd[row] >= 0) "kLower" else "kUpper")
      }
    }
    out$x <- x; out$col_dual <- sense * cd
    out$row_activity <- rv; out$row_dual <- sense * rd
    out$col_basis <- cb; out$row_basis <- rb
    fill_objective()
  } else {
    interior(); fill_objective()
  }
  out$status <- if (imprecise) "kUnknown" else "kOptimal"
  out
}
