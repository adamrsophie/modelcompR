# Stochastic differential equation systems and problems, after
# ModelingToolkitBase/src/systems/systems.jl (_mtkcompile,
# _brownians_to_noise_eqs, extract_brownians_to_noise_eqs,
# _poissonians_to_jumps, extract_poissonians_to_jumps, simplify_sde_system,
# __num_isdiag_noise), problems/sdeproblem.jl (SDEFunction, SDEProblem,
# calculate_noise_and_rate_prototype), systems/codegen.jl
# (generate_diffusion_function) and systems/diffeqs/basic_transformations.jl
# (stochastic_integral_transform, Girsanov_transform, noise_to_brownians).
#
# Noise equations are stored as in Julia: a list of expressions (diagonal
# noise, one independent Wiener process per equation) or a list matrix with
# one column per Wiener process (an n x 1 matrix is scalar noise).

# ---------------------------------------------------------------------------
# Lowering inside mtkcompile

# Coefficients of `vars` in `e`, which must be affine in them jointly.
# Returns list(coeffs = <list>, resid = <sym>), as Symbolics.linear_expansion.
.sde_linear_split <- function(e, vars, what, eq) {
  keys <- .keys(vars)
  coeffs <- lapply(vars, function(v) diff(e, v))
  for (i in seq_along(vars)) {
    if (length(intersect(.keys(.collect_vars(coeffs[[i]])), keys)))
      .stopf("Expected %s variables to appear linearly in equations. %s %s appears non-linearly in equation %s.",
             what, what, format(vars[[i]]), .fmt(eq))
  }
  resid <- .subst(e, .rule(vars, rep(list(.S0()), length(vars))))
  list(coeffs = coeffs, resid = resid)
}

# Brownian coefficients of each equation: (new equations, noise equations).
# Diagonal noise is detected when every row and every column of the noise
# matrix has at most one nonzero entry. Julia checks rows only, which would
# make two equations driven by one brownian independent; see the report.
.sde_brownians_to_noise_eqs <- function(eqs, bvars) {
  n <- length(eqs); m <- length(bvars)
  g <- matrix(list(.S0()), n, m)
  new_eqs <- eqs
  for (i in seq_len(n)) {
    sp <- .sde_linear_split(eqs[[i]]$rhs, bvars, "brownian", eqs[[i]])
    for (j in seq_len(m)) g[[i, j]] <- sp$coeffs[[j]]
    new_eqs[[i]] <- Equation(eqs[[i]]$lhs, sp$resid)
  }
  list(eqs = new_eqs, noise_eqs = .sde_noise_shape(g))
}

.sde_noise_shape <- function(g) {
  if (ncol(g) == 1L) return(g)
  nz <- matrix(!vapply(g, .is_zero, TRUE), nrow(g), ncol(g))
  if (all(rowSums(nz) <= 1L) && all(colSums(nz) <= 1L) && ncol(g) == nrow(g)) {
    return(lapply(seq_len(nrow(g)), function(i) {
      j <- which(nz[i, ])
      if (length(j)) g[[i, j]] else .S0()
    }))
  }
  g
}

# Noise equations as a list (diagonal) or list matrix (general or scalar).
.sde_noise_eqs <- function(sys) {
  ne <- .f(sys, "noise_eqs")
  if (is.null(ne)) return(NULL)
  if (inherits(ne, "symmatrix")) ne <- .symmatrix_to_list_matrix(ne)
  if (is.matrix(ne)) {
    out <- matrix(lapply(ne, as_sym), nrow(ne), ncol(ne))
    return(out)
  }
  lapply(.flatten_vars(lapply(ne, as_sym)), as_sym)
}

.sde_noise_type <- function(ne) {
  if (!is.matrix(ne)) return("diagonal")
  if (ncol(ne) == 1L) "scalar" else "general"
}

# Linear coefficient of each poissonian becomes a jump; pure jump equations
# D(x) ~ 0 are removed. After _poissonians_to_jumps.
.sde_poissonians_to_jumps <- function(sys, save_positions = c(FALSE, TRUE)) {
  pvars <- .f(sys, "poissonians")
  iv <- .f(sys, "iv")
  eqs <- .f(sys, "eqs")
  us <- .f(sys, "unknowns")
  new_jumps <- list()
  remove <- integer()
  extra_ps <- list()
  for (dN in pvars) {
    rate <- .meta_get(dN, "rate")
    if (is.null(rate)) next
    rate <- as_sym(rate)
    # a poissonian of a subsystem keeps its rate in the subsystem's names
    nm <- .name(dN)
    if (grepl(NAMESPACE_SEPARATOR, nm, fixed = TRUE)) {
      parts <- strsplit(nm, NAMESPACE_SEPARATOR, fixed = TRUE)[[1L]]
      for (pre in rev(parts[-length(parts)])) rate <- .namespace(rate, pre)
    }
    affects <- list()
    for (i in seq_along(eqs)) {
      eq <- eqs[[i]]
      if (!inherits(eq, "Equation") || !.is_derivative(eq$lhs)) next
      di <- .deriv_info(eq$lhs)
      if (di$order != 1L)
        .stopf("Higher-order derivative equation %s found in system with poissonians. Poissonians are only supported in first-order differential equations.", .fmt(eq))
      sp <- .sde_linear_split(eq$rhs, list(dN), "poissonian", eq)
      coeff <- sp$coeffs[[1L]]
      if (.is_zero(coeff)) next
      affects[[length(affects) + 1L]] <- Equation(di$var, Pre(di$var) + Pre(coeff))
      if (.is_zero(sp$resid)) remove <- c(remove, i) else eqs[[i]] <- Equation(eq$lhs, sp$resid)
    }
    if (!length(affects)) next
    rvars <- .collect_vars(rate)
    for (v in rvars) if (isparameter(v)) extra_ps[[length(extra_ps) + 1L]] <- v
    variable <- any(vapply(rvars, function(v) .contains_key(c(us, list(iv)), v), TRUE))
    new_jumps[[length(new_jumps) + 1L]] <- if (variable)
      VariableRateJump(rate, affects, save_positions = save_positions) else ConstantRateJump(rate, affects)
  }
  if (length(remove)) eqs <- eqs[-remove]
  .update_system(sys, eqs = eqs, jumps = c(new_jumps, .f(sys, "jumps")), poissonians = list(),
                 ps = .unique_syms(c(.f(sys, "ps"), extra_ps)))
}

# The lower_brownians hook of mtkcompile, after _mtkcompile and
# simplify_sde_system. `compile` runs the remaining (structural) pipeline.
.sde_lower <- function(sys, compile) {
  if (length(.f(sys, "poissonians"))) sys <- .sde_poissonians_to_jumps(sys)
  if (length(.f(sys, "jumps"))) {
    # jump systems keep every unknown: no structural simplification
    if (length(.f(sys, "brownians"))) {
      eqs <- Filter(function(e) inherits(e, "Equation"), .f(sys, "eqs"))
      low <- .sde_brownians_to_noise_eqs(eqs, .f(sys, "brownians"))
      sys <- .update_system(sys, eqs = low$eqs, noise_eqs = low$noise_eqs, brownians = list())
    }
    return(.update_system(sys, isscheduled = TRUE))
  }
  if (!is.null(.f(sys, "noise_eqs"))) sys <- noise_to_brownians(sys, names = "\u03B1\u2098\u209C\u2096")
  bvars <- .f(sys, "brownians")
  if (!length(bvars)) return(compile(sys))
  sys <- .update_system(sys, brownians = list())
  out <- compile(sys)
  eqs <- .f(out, "eqs")
  low <- .sde_brownians_to_noise_eqs(eqs, bvars)
  for (o in .f(out, "observed")) {
    if (length(intersect(.keys(.collect_vars(o$rhs)), .keys(bvars))))
      .stopf("observed equation %s depends on a brownian variable", .fmt(o))
  }
  sched <- .f(out, "schedule") %||% list()
  sched$var_sccs <- list(seq_along(low$eqs))
  .update_system(out, eqs = low$eqs, noise_eqs = low$noise_eqs, schedule = sched)
}
.mtk$hooks$lower_brownians <- .sde_lower

#' Convert noise equations to brownian variables
#'
#' `noise_to_brownians(sys)` returns a system equivalent to `sys` whose noise,
#' given as noise equations (`noise_eqs`), is written with brownian variables
#' added to the equations instead: `D(x) ~ f` with noise `g` becomes
#' `D(x) ~ f + g * α_1`.
#'
#' @param sys A flattened `System` with noise equations.
#' @param names A name prefix (brownians `<names>_1`, `<names>_2`, ...) or
#'   one name per brownian.
#' @return A `System` with brownians and no noise equations.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)"); parameters("a", "b")
#' sde <- SDESystem(list(D(x) ~ a * x), list(b * x), t, name = "gbm")
#' equations(noise_to_brownians(sde))
#' @export
noise_to_brownians <- function(sys, names = "\u03B1") {
  ne <- .sde_noise_eqs(sys)
  if (is.null(ne)) stop("Expected a system with `noise_eqs`.", call. = FALSE)
  if (length(.f(sys, "systems"))) stop("The system must be flattened.", call. = FALSE)
  nb <- if (is.matrix(ne)) ncol(ne) else length(ne)
  if (length(names) == 1L && nb != 1L) names <- paste0(names, "_", seq_len(nb))
  else if (length(names) == 1L) names <- paste0(names, "_1")
  if (length(names) != nb)
    .stopf("The system has %d brownian variables. Received %d names for the brownian variables.",
           nb, length(names))
  bvars <- lapply(names, Brownian)
  eqs <- .f(sys, "eqs")
  if (length(eqs) != (if (is.matrix(ne)) nrow(ne) else length(ne)))
    stop("the number of noise equations does not match the number of equations", call. = FALSE)
  new_eqs <- lapply(seq_along(eqs), function(i) {
    term <- if (is.matrix(ne)) Reduce(`+`, lapply(seq_len(nb), function(j) ne[[i, j]] * bvars[[j]]))
            else ne[[i]] * bvars[[i]]
    Equation(eqs[[i]]$lhs, eqs[[i]]$rhs + term)
  })
  noise_ps <- Filter(isparameter, .collect_vars(if (is.matrix(ne)) as.list(ne) else ne))
  .update_system(sys, eqs = new_eqs, brownians = bvars, noise_eqs = NULL,
                 ps = .unique_syms(c(.f(sys, "ps"), noise_ps)))
}

# ---------------------------------------------------------------------------
# Transformations

#' Convert between Itô and Stratonovich SDEs
#'
#' `stochastic_integral_transform(sys, -1/2)` converts an Itô system to the
#' equivalent Stratonovich system and `1/2` converts back, adding
#' `correction_factor * sum_k (J g_k) g_k` to the drift, where `g_k` are the
#' noise columns and `J` their Jacobian. For diagonal noise only the
#' derivative of each noise term with respect to its own state enters, which
#' is the correct correction; ModelingToolkit uses the full Jacobian here,
#' which agrees whenever each noise term depends only on its own state.
#'
#' @param sys A flattened `System` with noise equations.
#' @param correction_factor `-1/2` (Itô to Stratonovich) or `1/2`.
#' @return The transformed `System`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)"); parameters("a", "b")
#' sde <- SDESystem(list(D(x) ~ a * x), list(b * x), t, name = "gbm")
#' equations(stochastic_integral_transform(sde, -1/2))
#' @export
stochastic_integral_transform <- function(sys, correction_factor) {
  if (length(.f(sys, "systems"))) stop("The system must be flattened.", call. = FALSE)
  ne <- .sde_noise_eqs(sys)
  if (is.null(ne))
    stop("`stochastic_integral_transform` expects a system with noise_eqs. If your noise is specified using brownian variables, consider calling `mtkcompile`.",
         call. = FALSE)
  eqs <- .f(sys, "eqs")
  dvs <- .f(sys, "unknowns")
  others <- c(.f(sys, "ps"), list(.f(sys, "iv")))
  n <- length(dvs)
  corr <- rep(list(.S0()), n)
  if (!is.matrix(ne)) {
    J <- .partials(ne, dvs, others)
    for (i in seq_len(n)) corr[[i]] <- J[[i, i]] * ne[[i]]
  } else {
    for (k in seq_len(ncol(ne))) {
      col <- lapply(seq_len(n), function(i) ne[[i, k]])
      J <- .partials(col, dvs, others)
      for (i in seq_len(n)) for (j in seq_len(n)) corr[[i]] <- corr[[i]] + J[[i, j]] * col[[j]]
    }
  }
  cf <- as_sym(correction_factor)
  deqs <- lapply(seq_along(eqs), function(i)
    Equation(eqs[[i]]$lhs, eqs[[i]]$rhs + cf * symcompR::simplify(corr[[i]])))
  .update_system(sys, eqs = deqs)
}

#' Girsanov transformation of an SDE system
#'
#' Measure transformation that reduces the variance of Monte Carlo estimators
#' of `E[g(X_t)]`: the returned system has an extra state `θ(t)` with
#' initial value `θ0` and an observed `weight(t) ~ θ(t) / θ0`, so that
#' `E[g(X_t) weight_t]` estimates the same quantity with smaller variance for a
#' good choice of `u`. Reference: Kloeden, Platen and Schurz (2012),
#' Numerical solution of SDE through computer experiments.
#'
#' @param sys A `System` with noise equations.
#' @param u A scalar symbolic function of the unknowns defining the
#'   adjustable parameters of the transformation.
#' @param theta0 Initial value of `θ` (`θ0` in ModelingToolkit).
#' @return The transformed `System`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)"); parameters("a", "b")
#' sde <- SDESystem(list(D(x) ~ a * x), list(b * x), t, name = "gbm")
#' demod <- Girsanov_transform(sde, x, theta0 = 0.1)
#' equations(demod)
#' @export
Girsanov_transform <- function(sys, u, theta0 = 1) {
  iv <- .f(sys, "iv")
  Dt <- Differential(iv)
  theta <- Variable("\u03B8", iv)
  weight <- Variable("weight", iv)
  dvs <- .f(sys, "unknowns")
  others <- c(.f(sys, "ps"), list(iv))
  u <- as_sym(u)
  gr <- .partials(list(u), dvs, others)
  grad <- lapply(seq_along(dvs), function(j) gr[[1L, j]])
  ne <- .sde_noise_eqs(sys)
  n <- length(dvs)
  if (!is.matrix(ne)) {
    d <- lapply(seq_len(n), function(i) symcompR::simplify(-(ne[[i]] * grad[[i]]) / u))
    drift <- lapply(seq_len(n), function(i) ne[[i]] * d[[i]])
    G <- matrix(list(.S0()), n, n)
    for (i in seq_len(n)) G[[i, i]] <- ne[[i]]
  } else {
    m <- ncol(ne)
    d <- lapply(seq_len(m), function(k)
      symcompR::simplify(-Reduce(`+`, lapply(seq_len(n), function(i) ne[[i, k]] * grad[[i]])) / u))
    drift <- lapply(seq_len(n), function(i) Reduce(`+`, lapply(seq_len(m), function(k) ne[[i, k]] * d[[k]])))
    G <- ne
  }
  eqs <- .f(sys, "eqs")
  deqs <- lapply(seq_along(eqs), function(i) Equation(eqs[[i]]$lhs, eqs[[i]]$rhs - drift[[i]]))
  deqs[[length(deqs) + 1L]] <- Equation(Dt(theta), 0)
  newrow <- matrix(lapply(d, function(dk) theta * dk), 1L, length(d))
  G <- rbind(G, newrow)
  ics <- smap_set(.f(sys, "initial_conditions"), theta, theta0)
  out <- .update_system(sys, eqs = deqs, noise_eqs = G, unknowns = c(dvs, list(theta)),
                        initial_conditions = ics,
                        observed = c(list(Equation(weight, theta / theta0)), .f(sys, "observed")))
  par <- .f(sys, "parent")
  if (!is.null(par))
    out <- .update_system(out, parent = .update_system(par, unknowns = c(.f(par, "unknowns"), list(theta, weight))))
  if (iscomplete(sys)) out <- .update_system(out, index_cache = .build_index_cache(out))
  out
}

# ---------------------------------------------------------------------------
# Code generation

#' Generate the diffusion function of an SDE system
#'
#' Returns an R function `g(u, p, t)` evaluating the noise equations of a
#' compiled stochastic system: a vector for diagonal and scalar noise, a
#' matrix (one column per Wiener process) for general noise. For systems with
#' delays the function is `g(u, h, p, t)` with `h(p, t)` the history.
#'
#' @param sys A complete `System` with noise equations.
#' @param ... Unused.
#' @return An R function.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)"); parameters("a", "b")
#' B <- Brownian("B")
#' sys <- mtkcompile(System(list(D(x) ~ a * x + b * x * B), t, name = "gbm"))
#' g <- generate_diffusion_function(sys)
#' g(2, c(1, 0.5), 0)
#' @export
generate_diffusion_function <- function(sys, ...) {
  ne <- .sde_noise_eqs(sys)
  if (is.null(ne)) stop("the system has no noise equations", call. = FALSE)
  if (is.matrix(ne) && ncol(ne) == 1L) ne <- lapply(seq_len(nrow(ne)), function(i) ne[[i, 1L]])
  out <- if (is.matrix(ne)) "matrix" else "vector"
  if (isTRUE(.f(sys, "is_dde"))) return(.dde_build(ne, sys, out = out))
  .build_closure(ne, .codegen_args(sys), observed = .codegen_observed(sys), out = out)
}

# ---------------------------------------------------------------------------
# SDEFunction and SDEProblem

#' SDE functions and problems
#'
#' `SDEProblem(sys, op, tspan)` builds a stochastic differential equation
#' problem from a compiled system with brownian variables (or noise
#' equations); `op` gives initial values and parameters as for
#' [ODEProblem()]. `SDEFunction(sys)` builds the drift `f(u, p, t)` and
#' diffusion `g(u, p, t)`.
#'
#' `SDEProblem(f, g, u0, tspan, p)` with R functions builds a problem without
#' a system. Diagonal noise is the default; give `noise_rate_prototype` (a
#' matrix with one column per Wiener process) for general noise, or
#' `noise = "scalar"` for one Wiener process shared by all equations.
#' `SDEFunction(f, g, mass_matrix = M)` gives a mass matrix (singular for
#' algebraic equations), solved with the implicit methods
#' ([ImplicitEM()], [ImplicitRKMil()], ...).
#'
#' @param sys A compiled `System`, or the drift function `f(u, p, t)`.
#' @param ... For systems: `op`, `tspan` and solver options stored with the
#'   problem (`callback`, `tstops`, ...). For functions: `g`, `u0`, `tspan`,
#'   `p`, then options.
#' @param noise_rate_prototype,noise Noise shape for raw problems.
#' @param jac,tgrad Generate the drift Jacobian and time gradient.
#' @return An `SDEProblem` or `SDEFunction`.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("x(t)"); parameters("mu = 1.5", "sigma = 0.3")
#' B <- Brownian("B")
#' sys <- mtkcompile(System(list(D(x) ~ mu * x + sigma * x * B), t, name = "gbm"))
#' prob <- SDEProblem(sys, list(x ~ 1), c(0, 1))
#' sol <- solve(prob, EM(), dt = 0.01, seed = 1)
#' tail(sol[x], 1)
#' @export
SDEFunction <- function(sys, ..., jac = FALSE, tgrad = FALSE) {
  if (is.function(sys)) {
    dots <- list(...)
    return(structure(list(f = sys, g = dots[[1L]], jac = NULL, tgrad = NULL, mass_matrix = dots$mass_matrix,
                          sys = NULL, observed = NULL),
                     class = c("SDEFunction", "SciMLFunction")))
  }
  .need_complete(sys, "SDEFunction")
  .sde_check_noise(sys, "SDEFunction")
  fn <- ODEFunction(sys, jac = jac, tgrad = tgrad)
  fn$g <- generate_diffusion_function(sys)
  ne <- .sde_noise_eqs(sys)
  fn$noise_type <- .sde_noise_type(ne)
  fn$nbrownians <- if (is.matrix(ne)) ncol(ne) else length(ne)
  class(fn) <- c("SDEFunction", "SciMLFunction")
  fn
}

.sde_check_noise <- function(sys, what) {
  if (length(.f(sys, "brownians")))
    .stopf("%s needs a compiled system: call mtkcompile() to turn brownians into noise equations", what)
  if (is.null(.f(sys, "noise_eqs")))
    .stopf("A system without noise cannot be used to construct a %s; use ODEProblem()", what)
  if (length(jumps(sys))) .stopf("A system with jumps cannot be used to construct a %s; use JumpProblem()", what)
  if (isTRUE(.f(sys, "is_dde")) && what %in% c("SDEProblem", "SDEFunction"))
    .stopf("A system with delays cannot be used to construct a %s; use SDDEProblem()", what)
}

#' @rdname SDEFunction
#' @export
SDEProblem <- function(sys, ..., noise_rate_prototype = NULL, noise = NULL) {
  # SDEProblem(SDEFunction(f, g), u0, tspan, p) takes g from the function
  if (inherits(sys, "SDEFunction") && !(length(list(...)) && is.function(..1)))
    return(.sde_problem_raw(sys, sys$g, ..., noise_rate_prototype = noise_rate_prototype, noise = noise))
  if (is.function(sys) || inherits(sys, "SDEFunction"))
    return(.sde_problem_raw(sys, ..., noise_rate_prototype = noise_rate_prototype, noise = noise))
  .sde_problem_sys(sys, ...)
}

.sde_problem_raw <- function(f, g = NULL, u0, tspan, p = numeric(), ..., noise_rate_prototype = NULL,
                             noise = NULL) {
  fn <- if (inherits(f, "SDEFunction")) f else SDEFunction(f, g)
  type <- if (!is.null(noise_rate_prototype)) "general" else if (identical(noise, "scalar")) "scalar" else "diagonal"
  u0 <- as.numeric(unlist(u0))
  fn$noise_type <- type
  fn$nbrownians <- switch(type, general = ncol(as.matrix(noise_rate_prototype)), scalar = 1L,
                          diagonal = length(u0))
  .new_problem("SDEProblem", f = fn, u0 = u0, tspan = .as_tspan(tspan), p = p, sys = NULL,
               noise_type = type, kwargs = list(...))
}

.sde_problem_sys <- function(sys, op = list(), tspan = NULL, ..., check_initialization = TRUE,
                             jac = FALSE, tgrad = FALSE) {
  .need_complete(sys, "SDEProblem")
  .sde_check_noise(sys, "SDEProblem")
  prob <- ODEProblem(sys, op, tspan, ..., jac = jac, tgrad = tgrad,
                     check_initialization = check_initialization)
  fn <- SDEFunction(sys, jac = jac, tgrad = tgrad)
  for (nm in setdiff(names(prob$f), names(fn))) fn[[nm]] <- prob$f[[nm]]
  prob$f <- fn
  prob$noise_type <- fn$noise_type
  class(prob) <- c(if (inherits(prob, "InitializationFailure")) "InitializationFailure",
                   "SDEProblem", "SciMLProblem")
  prob
}

# ---------------------------------------------------------------------------
# Poissonian declarations, the @poissonians macro

#' Declare poissonian variables
#'
#' `poissonians(dN = lambda)` declares a poissonian `dN`, the increment of a
#' Poisson counting process with rate `lambda`, and assigns it in the calling
#' environment, as `@poissonians dN(λ)`. Written linearly in differential
#' equations, `D(x) ~ f + h * dN`, it becomes a jump `x ~ Pre(x) + Pre(h)`
#' at rate `lambda` during [mtkcompile()]: a constant rate jump when `lambda`
#' depends only on parameters, a variable rate jump otherwise.
#'
#' @param ... Named rates (numbers or expressions), or strings
#'   `"dN(lambda)"` evaluated in the calling environment. Given a system,
#'   `poissonians(sys)` returns its poissonian variables instead.
#' @param .env Environment to assign into.
#' @return Invisibly, a named list of the declared poissonians.
#' @examples
#' t <- t_nounits; D <- D_nounits
#' variables("X(t)"); parameters("lam = 2")
#' poissonians(dN = lam)
#' sys <- mtkcompile(System(list(D(X) ~ dN), t, name = "counter"))
#' length(jumps(sys))
#' @export
poissonians <- function(..., .env = parent.frame()) {
  dots <- list(...)
  if (length(dots) && inherits(dots[[1L]], "System"))
    return(.unique_syms(.collect_hier(dots[[1L]], function(s) .f(s, "poissonians"))))
  nms <- names(dots) %||% rep("", length(dots))
  out <- list()
  for (i in seq_along(dots)) {
    if (nzchar(nms[i])) {
      nm <- nms[i]; rate <- dots[[i]]
    } else {
      # the declaration string is the caller's own code, evaluated in the
      # caller's environment as the @poissonians macro would
      ex <- parse(text = dots[[i]], keep.source = FALSE)[[1L]]
      if (!is.call(ex)) .stopf("poissonians requires a rate expression: use dN(rate)")
      nm <- as.character(ex[[1L]]); rate <- eval(ex[[2L]], .env)
    }
    v <- Poissonian(nm, rate)
    assign(nm, v, envir = .env)
    out[[nm]] <- v
  }
  invisible(out)
}
