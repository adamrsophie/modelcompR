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
* New options: `loss_grad` of `solve_gradient()`, `save_idxs_storage` of
  `solve()`, and `retries`, `rng_streams` and `checkpoint` of ensembles;
  `cache_info()` reports the caches.
