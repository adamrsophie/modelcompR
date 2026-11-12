# Porting status

Reference: ModelingToolkit.jl 11.47.0 (commit `d6e93c5`), with
ModelingToolkitBase and SciCompDSL from the same repository, StateSelection.jl
(commit `29b3c90`, including ModelingToolkitTearing), and the SciML solver
packages listed in `inst/COPYRIGHTS`. Numerical reference values in
`tests/testthat/fixtures/julia/` come from ModelingToolkit 11.46.1 run by
`tools/julia/reference.jl`.

## API coverage

All 327 public names of ModelingToolkit and ModelingToolkitBase (their
`export` and `@public` lists, in `tools/mtk_exports.tsv`) are accounted for:

* 318 are ported, among them the JuMP, InfiniteOpt, CasADi and Pyomo
  dynamic optimization backends, whose constructors build the native
  collocation problem with the backend's transcription (`tools/stubs.tsv`
  is empty);
* 9 have no R meaning (`tools/not_applicable.tsv`: the Julia mutable cache
  functions, `CompilerOptions`, `SciMLLogging`, `StructuralTransformations`,
  the precompilation workloads and `MTKParameters`).

This counts names, not behavior: a ported name can support fewer cases than
its Julia counterpart. The behaviors not yet ported are listed below.
`Rscript tools/ledger.R` recomputes the table in `tools/ledger_status.tsv`.

Julia macros are functions in R:

| Julia | R |
|-------|---|
| `@variables`, `@parameters`, `@constants`, `@brownians`, `@discretes`, `@poissonians` | `variables()`, `parameters()`, `constants()`, `brownians()`, `discretes()`, `poissonians()`, or `Variable()`, `Parameter()`, ... |
| `@independent_variables` | `independent_variable()` |
| `@named` | the `name` argument |
| `@component`, `@mtkmodel` | R functions returning a `System`; `component()`, `mtkmodel()` |
| `@connector` | `connector()` |
| `@mtkcompile`, `@mtkbuild`, `@mtkcomplete` | `mtkcompile()`, `complete()` |
| `t`, `D` | `t_nounits`, `D_nounits` |
| `substitute` | `substitute_vars()` (base R's `substitute` is not masked) |
| `find_solvables!`, `modified_unknowns!`, `collect_vars!`, ... | the same names without `!` |

No modelcompR function masks a base R function (`tools/check_masking.R`,
`tests/testthat/test-namespace.R`). Two names are shared with symcompR,
`flatten()` and `compose()`; both pass anything that is not a `System` to
symcompR.

## Modules

| Area | Status | Notes |
| --- | --- | --- |
| Variables, parameters, metadata | Ported | Metadata travels as an attribute of the symcompR object. Numeric defaults become initial conditions, symbolic and `NA` defaults become bindings, as in MTK 11. Scopes `LocalScope`, `ParentScope`, `GlobalScope`. Undeclared symcompR symbols and function applications are unknowns, as symbols without metadata are in MTK. |
| `System` and hierarchy | Ported | One `System` type with the v11 fields; `compose`, `extend`, `flatten`, `complete`; namespacing with `₊`; `sys$sub$x`. |
| `mtkcompile` | Ported | Connection expansion, alias elimination (Bareiss), Pantelides index reduction, dummy derivative state selection, tearing, reassembly, observed sorting, analytic integration; the ModelingToolkitBase pipeline with `structural = FALSE`. Unknowns match Julia for the pendulum, the RC circuit and the reference DAEs. Clock inference over clocked partitions (the hypergraph of ModelingToolkitTearing, partitions as Julia's for the SampledData tutorial), array variables (scalarized elementwise; `+`, `*`, `%*%`, `sum` and math functions on arrays), `IfLifting` and `no_if_lift`, inline solves of linear blocks (`inline_linear_sccs`, `analytical_linear_scc_limit`) and, as an extension, of nonlinear blocks (`inline_nonlinear_sccs`), the implicit stage system of `solver_nlprob.jl` (`ODEProblem(..., nlstep = TRUE)`, not used by the native solvers), state machines (`initial_state`, `transition`, `activeState`, `entry`, `ticksInState`, `timeInState`; extracted as in MTK 11, simulated as an extension with Modelica's semantics, without hierarchy or `synchronize`). `eliminate_mm_zeros` defaults to `FALSE`; with `TRUE` events update the integration constants. |
| Code generation | Ported | R closures from symcompR expressions with CSE; Jacobians, time gradients, parameter Jacobians, mass matrices, Hessians; `build_function`, `register_symbolic`. |
| ODE solvers | Ported | Tsit5, BS3, DP5, Vern6, Vern7, Vern8, Vern9 (and `AutoVern6`, `AutoVern8`, `AutoVern9`), the `autodiff` choices `AutoForwardDiff`, `AutoFiniteDiff`, `AutoSparse` for the implicit methods, Euler, RK4, Rosenbrock23, Rodas4, Rodas5P, ImplicitEuler, Trapezoid, FBDF, the default algorithm with stiffness switching (`AutoTsit5`, `AutoVern7`), mass matrix DAEs, the `step_limiter` option. Tableaus are generated from OrdinaryDiffEq.jl (commit `4460ff8`, recorded in `R/solvers-tableaus.R`). FBDF has stability limit detection, the optional time filters and the controller options of OrdinaryDiffEqBDF 2.4.14; its accepted steps on Robertson and Van der Pol (mu = 1000) are within 5% of Julia's. |
| DAE problems | Ported | `DAEProblem` (also through `remake`), `DFBDF`, `DImplicitEuler`, Brown and Shampine initialization. |
| Events and callbacks | Ported | Symbolic continuous and discrete events with `Pre`, `ImperativeAffect`, `AssignmentAffect`, numerical `ContinuousCallback`, `VectorContinuousCallback`, `DiscreteCallback`, `PresetTimeCallback`, `PeriodicCallback`. Parameter history for observed expressions. |
| Discrete time | Ported | `Shift`, `ShiftIndex`, `Sample`, `Hold`, `Clock`, `DiscreteProblem`, `FunctionMap`, `ImplicitDiscreteProblem`, `IDSolve`. Hybrid sampled data systems are an extension; MTK 11 needs an external pass for them. |
| Initialization | Ported | `Initial` parameters, initialization systems, guesses, solvable parameters, `InitializationProblem`, consistency checks and MTK's messages for under and overdetermined cases. Initialization systems become the problem types of `get_initialization_problem_type`: `LinearProblem` for affine systems, `SCCNonlinearProblem` (`use_scc = TRUE`, the default; linear, homotopy and nonlinear blocks), `HomotopyProblem` for `homotopy()` nodes, `NonlinearLeastSquaresProblem` for under and overdetermined systems, and `HomotopyContinuationProblem` with `use_homotopy_continuation = TRUE` (an extension). `prob$ps[Initial(x)] <- v` and `prob$ps[p] <- v` reinitialize the problem; `remake(check_initialization = FALSE)` does not. |
| Nonlinear | Ported | `NonlinearProblem`, least squares, steady state, interval problems, SCC problems, `LinearProblem`, homotopy continuation; NewtonRaphson, TrustRegion (dogleg with the Simple, NLsolve, Fan and Bastin radius updates), LevenbergMarquardt, GaussNewton, Broyden, Klement, the polyalgorithms of NonlinearSolve in their order (with `start_index` for small problems), ITP, Bisection, Brent, Falsi. NonlinearSolve's default More subproblem of `TrustRegion` is not ported. |
| Optimization | Ported | `OptimizationProblem` with constraints and bounds, symbolic gradients and Hessians; `optim`, `nlminb` and an augmented Lagrangian; multiobjective problems by `NSGA2()` and `WeightedSumScan()` (Pareto set and front), integer and binary unknowns by `BranchAndBound()` over the continuous relaxation (globally optimal for convex relaxations). |
| Dynamic optimization | Ported | Native direct collocation (`DynamicOptProblem`) with Runge Kutta schemes in stage form: trapezoidal, Hermite Simpson, backward Euler, Radau IIA (ModelingToolkit's default tableau), Gauss Legendre, any tableau; integral (Lagrange) costs with the scheme's quadrature; the collocation polynomial between nodes; an interior point solver in the manner of Ipopt (filter line search, second order corrections, watchdog; no restoration phase, the augmented Lagrangian is the fallback). The JuMP, InfiniteOpt, CasADi and Pyomo problem and collocation constructors build it with their default transcriptions; JuMP and CasADi hold the inputs constant over a step, as in ModelingToolkit. |
| SDEs | Ported | Brownians lowered to noise equations, `SDEProblem`, EM, EulerHeun, LambaEM, LambaEulerHeun, RKMil, RKMilCommute, SRIW1, SOSRI, the additive noise methods SRA1, SRA2, SRA3, SOSRA, SOSRA2 and SKenCarp, the drift implicit ImplicitEM (STrapezoid, SImplicitMidpoint), ImplicitEulerHeun and ImplicitRKMil, the split step ISSEM and ISSEulerHeun, StochasticDiffEq's default algorithm with `alg_hints`, SDEs with algebraic equations (singular mass matrices, with StochasticDiffEq's `alg_mass_matrix_compatible` rules), ensembles, `stochastic_integral_transform`, `Girsanov_transform`, `noise_to_brownians`. Events inside a step are located on the linear interpolant, as StochasticDiffEq does. Noise equations of subsystems are dropped by `flatten`, as in ModelingToolkit. |
| Jumps | Ported | `MassActionJump`, `ConstantRateJump`, `VariableRateJump` (with rate bounds), `SymbolicMassActionJump`, poissonians; the aggregators `Direct`, `DirectFW`, `SortingDirect`, `RDirect`, `FRM`, `FRMFW`, `NRM`, `CCNRM`, `DirectCR`, `RSSA`, `RSSACR` and `Coevolve` (all of JumpProcesses' `JUMP_AGGREGATORS`) with dependency graphs built as ModelingToolkit builds them, and `NullAggregator` (JumpProcesses' `select_aggregator`); `SSAStepper` with discrete callbacks, tstops and the system's events (`reset_aggregated_jumps`); jumps coupled to ODEs and SDEs through the aggregator's callback and VR_FRM, with `save_positions`. Index based `MassActionJump`s with several reactions for problems without a system, `DiscreteProblem(u0, tspan, p)`. Spatial jumps: `SpatialMassActionJump`, `CartesianGrid`, `CartesianGridRej`, graphs as neighbor lists or adjacency matrices, `num_sites`, `neighbors`, `outdegree`, all hopping constant forms of `HopRates`, the aggregators `NSM` and `DirectCRDirect`, and `flatten_spatial()` (JumpProcesses' `flatten`). The variable rate aggregators `VR_FRM`, `VR_Direct`, `VR_DirectFW` (`vr_aggregator`), `SplitCoupledJumpProblem`, `JumpSet`, `RateBounds`, `JumpDependencies`, `get_num_majumps`, `needs_depgraph`, `needs_vartojumps_map`, `isolated_copy`. Tau leaping: `RegularJump`, `PureLeaping`, `SimpleTauLeaping`, `SimpleExplicitTauLeaping`, `SimpleImplicitTauLeaping`, `SimpleTrapezoidalLeaping`, `SimpleAdaptiveTauLeaping` (JumpProcesses) and `TauLeaping`, `CaoTauLeaping`, `ImplicitTauLeaping`, `ThetaTrapezoidalTauLeaping` (StochasticDiffEq). |
| Integrator interface (SciMLBase) | Ported | `init()`, `step_integrator()` (`step!`), `solve_integrator()` (`solve!`), `reinit()`, `set_u()`, `set_t()`, `set_proposed_dt()`, `get_proposed_dt()`, `set_abstol()`, `set_reltol()`, `add_saveat()`, `savevalues()`, `derivative_discontinuity()`, `get_du()`, `add_tstop()`, `terminate()`, `u_modified()`, for ODE, DAE and SDE problems; `solve()` is `init()` followed by `solve_integrator()`. |
| Callback library (DiffEqCallbacks) | Ported | All 22 public names: `SavingCallback` with `SavedValues`, `FunctionCallingCallback`, `LinearizingSavingCallback` with `IndependentlyLinearizedSolution`, `PositiveDomain`, `GeneralDomain`, `ManifoldProjection`, `AutoAbstol`, `StepsizeLimiter`, `TerminateSteadyState`, `ProbIntsUncertainty`, `AdaptiveProbIntsUncertainty`, `IterativeCallback`, `PeriodicCallback`, `PresetTimeCallback`, `IntegratingCallback`, `IntegratingSumCallback`, `IntegratingGKCallback`, `IntegratingGKSumCallback` with `IntegrandValues` and `IntegrandValuesSum`. |
| Ensembles (SciMLBase) | Ported | `EnsembleProblem`, `EnsembleSolution`, `EnsembleSummary`, the algorithms `EnsembleSerial`, `EnsembleThreads` (forked workers), `EnsembleDistributed` and `EnsembleSplitThreads` (socket clusters), `EnsembleGPUKernel` and `EnsembleGPUArray` (CPU workers), and the functions of `EnsembleAnalysis` (`timestep_*`, `timepoint_*`, `timeseries_steps_*`, `timeseries_point_*`, `componentwise_*`). |
| Delays and BVPs | Ported | `DDEProblem` with `MethodOfSteps` (discontinuities with orders propagated by constant lags, neutral problems, state dependent lags located after rejected steps), `SDDEProblem`, `BVProblem` with shooting, multiple shooting and MIRK4, MIRK5, MIRK6 with defect control and mesh refinement. |
| Connectors and components | Ported | `connector`, `connect`, `Flow`, `Stream`, `instream`, `domain_connect`, `component`, `mtkmodel`. |
| Linear analysis | Ported | `linearize` (Julia's form by default; `eliminate_algebraic = TRUE` reduces to differential states), analysis points, sensitivity and loop transfer functions, `linearize_symbolic`, `similarity_transform`. |
| Transformations | Ported | `liouville_transform`, `change_independent_variable`, `change_of_variables`, `add_accumulations`, `fractional_to_ordinary` and `linear_fractional_to_ordinary` (scalar and array forms), `convert_system_indepvar`, `modelingtoolkitize`. |
| PDE, units, debugging | Ported | `PDESystem` container (no discretization, as in MTK), unit metadata and validation, `debug_system` and assertions. |
| Sensitivity analysis (SciMLSensitivity.jl) | Ported | Forward sensitivities, the interpolating, backsolve, quadrature and Gauss adjoints with checkpointing, events, mass matrix DAEs, steady state adjoints, LSS, NILSS and NILSAS, second order and pathwise SDE sensitivities, `solve_gradient()`; see the SciMLSensitivity section. |

## SciMLSensitivity

Reference: SciMLSensitivity.jl 7.119 (`documentation/SciMLSensitivity.jl`),
ported in `R/sensitivity-*.R` and tested in `tests/testthat/test-sensitivity.R`
against analytic derivatives and central finite differences of solves.
Derivatives of the right-hand side come from the System: the Jacobians
`df/du`, `df/dp` and the vector Jacobian product, the gradient of
`lambda' f` with respect to the unknowns and parameters, are built with
`.partials()` and compiled to R. Problems built from R functions use their
`jac_fn` or central finite differences. Event affects and conditions are
differentiated by central differences on a stand-in integrator.

| Julia | R |
| --- | --- |
| `Zygote.gradient((u0, p) -> loss(solve(remake(prob; u0, p), alg; sensealg)), u0, p)` | `solve_gradient(loss, prob, alg, ..., sensealg = )`; the loss is an R function of the solution, its cotangent with respect to the saved values is taken by central differences |
| `solve(prob, alg; sensealg = ...)` | the same; for ODE problems it records the events for the adjoint and keeps a dense solution |
| `ODEForwardSensitivityProblem`, `extract_local_sensitivities` | the same names; `ForwardDiffSensitivity()` takes the seed matrices `du0`, `dp` and differentiates events |
| `adjoint_sensitivities`, `ODEAdjointProblem` | the same names, with `dgdu_discrete`, `dgdp_discrete`, `dgdu_continuous`, `dgdp_continuous`, `g`, `checkpoints`, `no_start` |
| `ForwardSensitivity`, `ForwardDiffSensitivity` | forward sensitivity equations of the augmented system (both); u0 dependence on p through initialization is included |
| `InterpolatingAdjoint(checkpointing)`, `BacksolveAdjoint(checkpointing)`, `QuadratureAdjoint`, `GaussAdjoint`, `GaussKronrodAdjoint` | ported |
| `SteadyStateAdjoint` | ported for `NonlinearProblem` and `SteadyStateProblem` (implicit function theorem) |
| callback tracking (`track_callbacks`, `setup_reverse_callbacks`) | continuous, vector continuous, discrete and preset time events that change states or parameters, moving event times, terminating events, samples saved at event times |
| mass matrix DAE adjoints | semi-explicit index 1 constant mass matrices, with the jumps of `ReverseLossCallback` |
| `ForwardLSS`, `AdjointLSS`, `TimeDilation`, `CosWindowing`, `Cos2Windowing`, `ForwardLSSProblem`, `AdjointLSSProblem`, `NILSS`, `NILSSProblem`, `NILSAS`, `NILSASProblem`, `shadow_forward`, `shadow_adjoint` | ported; the LSS Schur complement is solved by block tridiagonal elimination |
| `second_order_sensitivities`, `second_order_sensitivity_product`, `ForwardDiffOverAdjoint` | central differences of the adjoint gradient (finite differences over adjoint) |
| SDE sensitivities | pathwise forward sensitivities of `SDEProblem` (fixed step method, same seed) |

Automatic differentiation and hardware backends have no R counterpart:
`autojacvec = ZygoteVJP()`, `EnzymeVJP()`, `ReverseDiffVJP()`,
`TrackerVJP()`, `MooncakeVJP()` and `ReactantVJP()` map to `autojacvec =
TRUE` (symbolic vector Jacobian products from the System, the default) or
`FALSE` (finite differences); the discretize then optimize adjoints
`ZygoteAdjoint`, `ReverseDiffAdjoint`, `TrackerAdjoint`, `EnzymeAdjoint`
and `MooncakeAdjoint` (automatic differentiation through the solver) map
to the continuous methods above, `ForwardDiffSensitivity()` for the forward
case; GPU arrays and Reactant compilation do not apply.

Not ported:

* Fully implicit `DAEProblem` adjoints (`DAEAdjointProblem`,
  dae_adjoint.jl): the R `DAEProblem` solvers give no dense interpolant of
  `du`; mass matrix DAEs are supported.
* SDE and RODE adjoints (`SDEAdjointProblem`, `RODEAdjointProblem`,
  `StochasticTransformedFunction`): they need the replayed Brownian path and
  the Ito to Stratonovich correction; SDEs get pathwise forward
  sensitivities.
* `OptimizationAdjoint`, `UnconstrainedOptimizationAdjoint`
  (optimization_adjoint.jl): they need the KKT system and the active set of
  an optimization solution, which the R optimizers do not report.
* `SundialsAdjoint`: CVODES and IDAS are C libraries.
* `dgdt_discrete` (losses that depend on `sol$t` itself), losses that read
  parameters in force before an event changed them (an error), and losses
  that interpolate `sol(t)` (an error; give `saveat`).
* Events in mass matrix problems, checkpointing together with events,
  `QuadratureAdjoint` with events whose parameter changes depend on the
  state, and decreasing time spans together with events, checkpointing or
  `QuadratureAdjoint` (all raise errors). Decreasing time spans work for
  the forward methods and the other adjoints.
* `ODELocalSensitivityProblem` (deprecated alias in Julia),
  `SensitivityAlg` (deprecated), the internals `SensitivityFunction`,
  `AdjointSensitivityIntegrand` and the VJP traits.

## Not yet ported

* `TrustRegion` with More's subproblem (`TrustRegionSubproblem.More`, the
  default of NonlinearSolve); the dogleg subproblem is used.
* The restoration phase of Ipopt in the dynamic optimization solver; a hard
  start (such as the cart pole from rest) needs an `initial_trajectory`.
* The public names of StochasticDiffEq, NonlinearSolve, SciMLBase and
  Optimization marked missing in `tools/upstream_status.tsv` (see
  `Rscript tools/ledger.R`), and the solver families of the OrdinaryDiffEq
  sublibraries beyond those listed in the ODE solvers row.

## Known limitations

* Variable kinds are found by name for leaves that lost their metadata, so a
  plain symcompR `Symbol("a")` takes the kind of an earlier
  `Parameter("a")` declared in the same session. Declare variables with
  `Variable()` and `Parameter()`, or use distinct names.
* Inside initialization systems `Initial(x)` is represented by a parameter
  named by the symbolic key of `Initial(x)`.
* Hierarchical state machines and `transition(synchronize = TRUE)` are not
  supported; ModelingToolkit 11 itself only extracts state machines.
* With `inline_nonlinear_sccs = TRUE` the warm start of each nested Newton
  solve belongs to the compiled system and is shared by all problems built
  from it; blocks need a locally unique solution.
* Undeclared function heads are classified with a snapshot of symcompR's
  mathematical functions (`R/math-heads.R`, regenerated by
  `tools/gen_math_heads.R`); a function added to symcompR later is treated
  as a variable until the snapshot is regenerated.

## Deliberate differences

* `step!` and `solve!` are `step_integrator()` and `solve_integrator()`
  (`stats::step` exists, and `!` is not part of R names). The integrator
  of `init()` is an environment; `integ$interp(t)` interpolates inside the
  last step (Julia's `integrator(t)`).
* `SplitCoupledJumpProblem`: the jump in which only the control problem
  fires applies the control problem's affect; JumpProcesses'
  `build_split_jumps` applies the first problem's affect there.
* The parallel ensembles give each trajectory its own L'Ecuyer-CMRG random
  stream; `EnsembleSerial()` draws all trajectories from one stream in
  turn, as before.

* `IfLifting` stops at the jumps of `%%` (Julia's `mod`) with the root
  function `sin(pi x / y)`. Symbolics registers `mod`, `rem` and `div` with
  `_triangle(x, y)`, which is never positive, so ModelingToolkit's event
  for them never fires. `%/%` (Julia's `div`), `floor` and `ceiling` stop at
  the integers of their argument; Julia's `rem` has no symcompR counterpart.
* Solutions are R closures: `sol(t)` interpolates, `sol[x]` and `sol[expr]`
  index symbolically, `sol$t` and `sol$u` hold the saved values.
* At an event time the first of the two saved rows uses the parameters
  before the event; Julia uses the parameters after it for both.
* Falling below the minimum step returns `"DtLessThanMin"`.
* Diagonal noise requires one entry per row and per column of the noise
  matrix.
* `stats$nf` counts the stages the Verner methods compute only for
  interpolation (OrdinaryDiffEq's `_ode_addsteps!` calls `f` without
  counting them).
* Implicit methods reuse the Jacobian across factorizations as Julia does
  for in place problems; for out of place problems Julia recomputes it with
  every new factorization.
* FBDF applies its `step_limiter` (Julia's FBDF stores `step_limiter!` but
  never calls it); the `step_limiter` option of `solve()` is applied after
  each accepted step by every method, as in OrdinaryDiffEqCore.
* Dependent delays: after a rejected step DelayDiffEq evaluates the
  current interpolant with Theta measured from `tprev`, the start of the
  previous step, so its search can miss a discontinuity; the port searches
  the attempted step on `[t, t + dt]`.
* MIRK methods solve the collocation equations to a residual of 1e-10
  (`nlsolve_kwargs`); BoundaryValueDiffEq passes `abstol` (the defect
  tolerance) to the nonlinear solver. MIRK5's continuous extension, as
  tabulated in BoundaryValueDiffEq, has a defect of order 4 rather than 5,
  so MIRK5 refines more than MIRK4 on some problems.
* Only the jump callbacks save around their events as DiffEqBase does when
  `saveat` is given (forced saves); other callbacks keep the behavior
  above. With `save_everystep` every callback honors `save_positions`.
* An inconsistent initial state of an SDE with a singular mass matrix and
  without a system is made consistent (StochasticDiffEq's CheckInit raises
  an error). `ISSEM`, the default for such SDEs with non diagonal noise,
  keeps the algebraic equations; StochasticDiffEq's ISSEM gives wrong
  results with mass matrices. `SKenCarp` accepts mass matrices (in
  StochasticDiffEq it errors on them).
* `CaoTauLeaping` chooses each step with Cao's formula (the net changes of
  the jumps are read from `c`); StochasticDiffEq's CaoTauLeaping leaves the
  step size unset. A rejected `TauLeaping` step thins its Poisson counts
  binomially for the shorter step; the remainder of the rejected interval
  is not kept after the step is accepted.
* SortingDirect learns its search order per solve; JumpProcesses keeps it in
  the problem across solves (same law).
