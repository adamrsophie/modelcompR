# modelcompR 0.1.0

* First version: a native R port of ModelingToolkit.jl 11.47 on top of
  symcompR. All 327 public names of ModelingToolkit and ModelingToolkitBase
  are accounted for (319 ported, 8 not applicable).
* The SciML packages it builds on are ported natively as well:
  OrdinaryDiffEq with all of its solver sublibraries, StochasticDiffEq,
  JumpProcesses, DelayDiffEq, NonlinearSolve, LinearSolve, SciMLOperators,
  SciMLBase, DiffEqCallbacks, DiffEqDevTools, GlobalDiffEq,
  ImplicitDiscreteSolve, SciMLSensitivity, GlobalSensitivity,
  QuasiMonteCarlo, and Optimization with its solver wrappers (Optim, NLopt,
  Ipopt, MadNLP, HiGHS through SciPy, BlackBoxOptim, CMA-ES, Manopt, NOMAD,
  PRIMA, Metaheuristics and the others). Every one of their 1833 public
  names is accounted for (1522 ported, 72 aliases of ported algorithms, 239
  without an R meaning); `tools/upstream_status.tsv` records the status of
  each name.
* Where an upstream method is biased or wrong, the port follows the
  correct method and documents the difference: for example the randomized
  quasi Monte Carlo designs and corrected estimators of `gsa()`, the Ito
  correction of the SDE adjoints, and complex states in the Radau and
  Gauss-Legendre methods.
* Beyond the reference packages, the sensitivity methods differentiate
  events in mass matrix problems and SDEs, `DAEProblem`s with events,
  losses that read `sol$t`, interpolate `sol(t)` or read parameters before
  an event changed them, and decreasing time spans with events and
  checkpointing.
* Numerical contracts are checked against their definitions: mass
  matrices in every representation, observed and parameter chains of any
  length, dual numbers that follow base R, linear solves of sparse and low
  rank operators, residual based retcodes of the eigensolvers, integer
  domains and the objective sense of the optimizers, failure accounting
  of ensembles, and refusal of gradients of failed solves. Regression
  coefficients and partial correlations of `gsa()` are invariant to
  shifts, and the DeltaMoment interval uses the bootstrap spread.
* `arithmetic_mode("native")` runs the solvers with R arithmetic and LAPACK
  instead of the bit exact replay of Julia's floating point operations.
* Sparse Jacobians (`jac` returning a sparse matrix, or `AutoSparse()` with
  a `sparsity` pattern) stay sparse through the implicit and Rosenbrock
  methods; `ODEFunction(f, jac = J)` now uses J.
* Inline linear and nonlinear solves follow the branch of each integrator
  and have implicit derivatives.
* New options: `loss_grad` and `loss_rows` of `solve_gradient()`,
  `save_idxs_storage` of `solve()`, and `retries`, `rng_streams`,
  `checkpoint` and `run_id` of ensembles; `cache_info()` reports the caches.
* The derivative of a loss is taken in every saved value unless
  `loss_rows` or `loss_grad` says otherwise; it is never guessed from probes
  of the loss. Gradients see the `save_idxs` view of the loss, and a failed
  adjoint or replay stage raises a `sens_stage_failure` condition.
* Solutions and gradients carry an `execution` record: the algorithm, the
  arithmetic, the Jacobian and mass routes, tolerances, storage, the saved
  view, the derivative route of the loss and the registrations in force.
* Sparse structures: symmetric and unit triangular `Matrix` patterns are
  mirrored before coloring, `jac_prototype` marks every stored entry
  structural, the Krylov methods keep `Matrix` inputs and sparse
  preconditioners sparse, constant sparse masses keep sparse products and
  Newton matrices, and `generate_jacobian(sparse = TRUE)` builds only the
  structural nonzeros.
* Masses given as SciML operators that change with the state or time are
  supported (nonsingular ones), and non-finite masses are rejected.
* Krylov methods work at any scale of the right hand side and report
  Success only for the requested residual of the original system; low rank
  updates fall back to the assembled matrix when the base is singular;
  Arnoldi uses thick restarts and Rayleigh Ritz extraction and reports
  Success only with every requested pair converged.
* Regression sensitivity standardizes inputs and outputs, ranks ties by
  their average, and gives NA for constant inputs; Shapley effects keep
  their sampling design apart from the random numbers of the model.
* Ensemble checkpoints are versioned and written atomically, record the
  configuration and an identity of the run, and resume with the stored
  random streams.
* Inline nonlinear blocks commit their branch only at accepted states, and
  sensitivities, observed values and the integrator interface replay it;
  problems keep the function registrations they were built with.
* `native` arithmetic keeps complex states in the implicit methods.
