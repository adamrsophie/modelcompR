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
