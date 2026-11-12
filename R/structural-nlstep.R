# The stage system of implicit solvers, after ModelingToolkit.jl v11.47 (MIT)
# src/systems/solver_nlprob.jl (generate_ODENLStepData, inner_nlsystem,
# assemble_nlstep_data, generate_nlprobmap, ODE_GAMMA, ODE_C).
#
# An implicit stage of an ODE solver with mass matrix M solves
#   0 = outer_tmp + g1 f(g2 z + inner_tmp, p, c) - g3 M z
# for z, with the solver's coefficients g1, g2, g3, the stage time c and the
# buffers outer_tmp and inner_tmp as parameters. `ODEProblem(sys, ...,
# nlstep = TRUE)` compiles this system with mtkcompile(), so that tearing
# reduces it, and stores the resulting nonlinear problem with the ODE
# function as `f$nlstep_data`, as ODENLStepData in SciMLBase. The native
# implicit solvers of modelcompR do not use it; it is available to code
# that steps a problem itself.

.nlstep_param <- function(name) Parameter(name)

# inner_nlsystem: the stage system with the unknowns of sys as stage
# unknowns and the solver coefficients as parameters.
.inner_nlsystem <- function(sys, mm, compile = TRUE) {
  dvs <- unknowns(sys)
  eqs <- full_equations(sys)
  N <- length(dvs)
  if (length(eqs) != N) stop("the stage system needs as many equations as unknowns", call. = FALSE)
  g <- lapply(1:3, function(i) .nlstep_param(paste0("\u03B3", c("\u2081", "\u2082", "\u2083")[i], "\u2098\u209C\u2096")))
  cc <- .nlstep_param("c\u2098\u209C\u2096")
  outer <- lapply(seq_len(N), function(i) .nlstep_param(paste0("outer_tmp\u2098\u209C\u2096[", i, "]")))
  inner <- lapply(seq_len(N), function(i) .nlstep_param(paste0("inner_tmp\u2098\u209C\u2096[", i, "]")))
  rule <- .rule(c(dvs, list(.f(sys, "iv"))), c(lapply(seq_len(N), function(i) g[[2L]] * dvs[[i]] + inner[[i]]), list(cc)))
  rhs <- lapply(eqs, function(e) .subst(e$rhs, rule))
  new_eqs <- lapply(seq_len(N), function(i) {
    mz <- Reduce(`+`, lapply(seq_len(N), function(j) mm[i, j] * dvs[[j]]))
    Equation(.S0(), outer[[i]] + g[[1L]] * rhs[[i]] - g[[3L]] * mz)
  })
  ps <- c(parameters(sys), g, list(cc), inner, outer)
  nlsys <- System(new_eqs, NULL, dvs, ps, name = "nlsys", checks = FALSE)
  nlsys <- if (compile) mtkcompile(nlsys) else complete(nlsys)
  list(sys = nlsys, gamma = g, c = cc, outer = outer, inner = inner)
}

# generate_ODENLStepData and assemble_nlstep_data: the stage problem with
# g1 = g2 = g3 = 1, c = 0 and zero buffers, the positions of its unknowns
# among those of sys, setters for the solver coefficients and buffers, and
# the map from a stage solution to the unknowns of sys.
.generate_ode_nlstep_data <- function(sys, u0, p, mm = NULL, compile = TRUE) {
  if (is.null(mm)) mm <- calculate_massmatrix(sys)
  N <- length(unknowns(sys))
  if (!is.matrix(mm)) mm <- diag(N)
  inn <- .inner_nlsystem(sys, mm, compile)
  nlsys <- inn$sys
  op <- list()
  for (gm in inn$gamma) op[[length(op) + 1L]] <- Equation(gm, 1)
  op[[length(op) + 1L]] <- Equation(inn$c, 0)
  for (v in c(inn$outer, inn$inner)) op[[length(op) + 1L]] <- Equation(v, 0)
  vals <- c(as.list(u0), as.list(p))
  sk <- c(.keys(unknowns(sys)), .keys(.problem_params(sys)))
  for (v in c(unknowns(nlsys), .problem_params(nlsys))) {
    i <- match(.key(v), sk)
    if (!is.na(i) && !.contains_key(c(inn$gamma, list(inn$c), inn$outer, inn$inner), v))
      op[[length(op) + 1L]] <- Equation(v, vals[[i]])
  }
  nlprob <- NonlinearProblem(nlsys, op)
  pk <- .keys(.problem_params(nlsys))
  obs <- build_explicit_observed_function(nlsys, unknowns(sys))
  setter <- function(vars) {
    idx <- match(.keys(vars), pk)
    function(prob, values) {
      prob$p[idx[!is.na(idx)]] <- values[!is.na(idx)]
      prob
    }
  }
  structure(list(nlprob = nlprob, subsetidxs = match(.keys(unknowns(nlsys)), .keys(unknowns(sys))),
                 set_gamma_c = setter(c(inn$gamma, list(inn$c))), set_outer_tmp = setter(inn$outer),
                 set_inner_tmp = setter(inn$inner),
                 nlprobmap = function(u, p = nlprob$p) obs(u, p)),
            class = "ODENLStepData")
}
