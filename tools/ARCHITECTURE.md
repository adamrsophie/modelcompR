# modelcompR architecture and contributor contract

modelcompR is a native R port of ModelingToolkit.jl v11.47 (MIT) built on the
pure R computer algebra system symcompR (SymPy port,
https://github.com/adamrsophie/symcompR). No compiled code; base R plus
symcompR, stats, graphics, utils. The goal is a complete port: every name in
`tools/mtk_exports.tsv` is ported, renamed, or explicitly not applicable.

Reference sources (read only, under documentation/, gitignored):
ModelingToolkit.jl (v11.47, lib/ModelingToolkitBase, lib/SciCompDSL, docs, test),
ModelingToolkit-v10 (v10.32.1), StateSelection.jl (newest, AGPL-3 since 2025-12-01;
lib/ModelingToolkitTearing), BipartiteGraphs.jl, OrdinaryDiffEq.jl, SciMLBase.jl,
NonlinearSolve.jl, ModelingToolkitStandardLibrary.jl, SciMLSensitivity.jl.
A Julia reference environment is at tools/julia (`julia --project=tools/julia`).

## Licensing rule

Files translated from StateSelection.jl (including lib/ModelingToolkitTearing)
start with this header and are AGPL-3:

    # Ported from StateSelection.jl (https://github.com/JuliaComputing/StateSelection.jl),
    # commit 29b3c90. Portions are licensed under the GNU AGPL, version 3; see
    # inst/COPYRIGHTS and inst/STATESELECTION-AGPL-LICENSE.

All other files are GPL-3 and may port MIT code from the other repositories.

## Symbolic layer (symcompR)

* Expressions: immutable S3 lists of class "sym" with fields type, args, name,
  value, extra, key. `key` is a canonical string; equality and every map in
  modelcompR use it (`.key(e)`). Read fields with `.type()`, `.args()`,
  `.name()` (helpers in R/aaa-utils.R), never `$<-`.
* `x(t)` is `Function("x")(t)`: node type "Function", name "x", args [t].
  Parameters and time independent variables are `Symbol`s.
* `Derivative(x(t), t)` stays unevaluated; nested derivatives merge into
  `Derivative(x(t), (t, n))`. Use `.deriv_info(d)` -> list(var, iv, order).
  Build with `.make_derivative(var, iv, order)`.
* `xreplace` is top down: put derivative nodes in the same rule as the
  variables when lowering. Use `.subst(e, .rule(from, to))`,
  `.subst_fixpoint()`, never `subs()` for variable replacement.
* `diff(e, v)` works for v a Symbol, an applied function x(t), or a
  Derivative node. For Jacobians use `.partials(exprs, wrt, others)`, which
  holds every variable independent (placeholder symbols) and differentiates
  only incident variables.
* symcompR registers S3 methods on base generics (`solve.sym`, `diff.sym`,
  `solve.list`); modelcompR does the same (`solve.ODEProblem`, ...). Never
  mask a base or symcompR function unless it delegates for foreign
  arguments (see `flatten()` in R/system.R).

## Variables and metadata (R/variables.R)

* `Variable(name, iv, default=, guess=, bounds=, ...)`, `Parameter()`,
  `Constant()`, `Brownian()`, `Poissonian()`, `Discrete()`;
  string forms `variables("x(t) = 1")`, `parameters("k = 0.5")` assign into
  the caller. `t_nounits`, `D_nounits` are the default time and derivative.
* Kinds live in the registry `.mtk$vars[[name]]$kind`: "unknown",
  "parameter", "constant", "discrete", "brownian", "poissonian", "iv",
  "operator". `.var_kind(e)`, `.is_declared(e)`, `.is_varcall(e)`.
* Metadata is an attribute `attr(e, "mtk_meta")` that survives symcompR
  arithmetic; `.meta(e)`, `.meta_get(e, field)`; the registry is a fallback.
* Operators `Pre(x)`, `Initial(x)` are symcompR function nodes whose head is
  registered with kind "operator" (`.op_apply`, `.is_op(e, "Pre")`).
  `.collect_vars()` walks through operators and returns variable leaves and
  derivative nodes.

## Equations (R/equations.R)

* `Equation(lhs, rhs)`; formulas `D(x) ~ -k*x` are converted by
  `as_equation()` / `.as_equations()` (evaluated in the formula environment).
* `Inequality(lhs, rhs, "<=")`, infix `%<~%`, `%>~%`.
* `Differential(t)`, `expand_derivatives()`, `Pre`, `Initial`, `EvalAt`.

## System (R/system.R)

* S3 list of class "System" with the v11 fields in `.SYSTEM_FIELDS`. Read
  fields with `.f(sys, "eqs")` (or `.subset2`), never `sys$eqs`: `$` is the
  namespaced variable access (`sys$x`, `sys$sub$x`), as getproperty in Julia.
  Copy with `.update_system(sys, field = value)`.
* Variable maps (`bindings`, `initial_conditions`, `guesses`) are "symmaps":
  named lists keyed by symcompR key with entries list(var, value); helpers
  `symmap()`, `smap_set/get/has/del/merge`, `as_symmap()`.
  Numeric defaults -> initial_conditions; symbolic or NA defaults -> bindings.
* Hierarchical accessors include subsystems, namespaced with "₊"
  (NAMESPACE_SEPARATOR): `equations()`, `unknowns()`, `parameters(sys)`
  (= `.parameters_of(sys)`), `observed()`, `bindings()`,
  `initial_conditions()`, `guesses()`, ... `renamespace(sys, x)`.
* `flatten()`, `complete()`, `compose()`, `extend()`.

## Compilation (R/compile.R)

`mtkcompile(sys)`:
1. `expand_connections(sys)` (hook `expand_connections`)
2. `flatten(sys)`; hook `lower_brownians` for stochastic systems
3. hook `structural_compile(sys, inputs, outputs, disturbance_inputs,
   fully_determined, simplify, ...)` when present and the system has an
   independent variable; otherwise `.base_compile()` (MTKBase pipeline:
   order lowering with `default_toterm` names xˍt, trivial observed).
4. sets `parent`, `isscheduled = TRUE`, then `complete()`.
A compiled ODE/DAE system has equations `D(x) ~ f` (differential) and
`0 ~ g` (algebraic); `unknowns()` in the same order as the equations' states;
`observed()` sorted so every lhs is defined before use.

## Code generation (R/codegen.R)

`.build_closure(exprs, args, observed, out = "vector"|"matrix"|"scalar")`
with `args = list(list(name = "u", vars = <list>), list(name = "p", vars =
<list>), list(name = "t", vars = <sym>))` returns an R closure; observed
equations become local assignments; CSE via symcompR. Public:
`build_function()`, `generate_rhs()`, `calculate_jacobian()`, ...
`register_symbolic(name, fn)` for user functions.
Problem code uses `.codegen_args(sys)` and `.codegen_observed(sys)`
(observed plus bound parameter equations); the parameter vector is
`.problem_params(sys)` (parameters minus those bound to expressions).

## Problems, solutions and solvers

* Problems (R/problems.R): `.new_problem("ODEProblem", f = <ODEFunction>,
  u0 = named numeric, tspan, p = named numeric, sys, kwargs)`; class
  c("ODEProblem", "SciMLProblem"). `.process_op(sys, op)` resolves values.
  Hook `initialize_problem(prob, vals, check)` (initialization module)
  replaces `.simple_initialize()`. Hook `problem_callbacks(sys)` returns
  callbacks from symbolic events.
* `solve(prob, alg, ...)` is an S3 method on base::solve:
  `solve.ODEProblem(a, b, ...)`. Algorithms are objects of class
  c("<Name>", "ODEAlgorithm") created by `.new_alg()`; see the header of
  R/solvers-ode.R for the step/interp/init contract. The driver
  `.integrate(prob, alg, opts)` does step control, saving, tstops and
  callbacks.
* Solutions (R/solution.R) are closures of class c("ODESolution",
  "SciMLSolution", "function"): `sol(t)`, `sol$t`, `sol$u`, `sol$retcode`,
  `sol[x]`, `sol[expr]`, `sol[list(x, y)]`. Build with `.new_solution(cls,
  data)` where data has t, U (matrix, rows = times), prob, sys, retcode.
* Numerical callbacks (R/callbacks.R): `ContinuousCallback`,
  `VectorContinuousCallback`, `DiscreteCallback`, `PresetTimeCallback`,
  `PeriodicCallback`, `CallbackSet`; affects modify the integrator
  environment (`integ$u`, `integ$p`, `terminate(integ)`).

## Extension hooks

Modules register functions in `.mtk$hooks` from their own file at load time:

    .mtk$hooks$problem_callbacks <- function(sys) ...

Top-level assignments into `.mtk$hooks` run when the package namespace is
built, because R/aaa-utils.R (which creates `.mtk`) is collated first
(alphabetical). Name files so they sort after aaa-utils.R.

Hook signatures that differ from the plain form: `lower_brownians(sys, compile)`
runs for brownians, poissonians, jumps and noise_eqs and receives the rest of
the compile pipeline as `compile`; `jump_exprs(jumps)` returns expressions
whose variables System() collects. `linearize()` defaults to Julia's form
(algebraic unknowns kept in the state); `eliminate_algebraic = TRUE` removes them.

Hooks in use: expand_connections, namespace_connection, namespace_event, lower_brownians, jump_exprs,
compile_shift, compile_clocked, shift2term, check_units, rhs_assertions, initialize_parameters, complete,
structural_compile, initialize_problem, problem_callbacks,
normalize_continuous_events, normalize_discrete_events, event_parameters.

## Rules for every module

* No function name may mask a base R function: no top-level definition, exported
  or internal, may reuse a name exported by base, stats, utils, graphics,
  grDevices, methods or datasets. Use an S3 method on the base generic when the
  meaning is the same (solve, plot, print, summary, simulate, ...), otherwise a
  different name (`substitute_vars`, not `substitute`). Run
  `Rscript tools/check_masking.R` before you finish; it must print
  "No function masks a base R function".

* Port faithfully from the v11.47 sources (cite the Julia file and function
  in a comment at the top of each R file). Match MTK names; R-incompatible
  macros become functions (`@mtkmodel` -> function factories, `@named` ->
  `name =` argument).
* Write roxygen2 docs with `@export` and a runnable `@examples` for every
  exported function. Do NOT run roxygen2::roxygenise() or
  devtools::document(); the lead regenerates NAMESPACE. Test with
  `devtools::load_all()` (exports everything) and testthat.
* One test file per module: tests/testthat/test-<module>.R, with analytic or
  Julia reference values. Keep each test file under 60 s.
* Writing style in code comments and docs: no em or en dashes and no spaced
  hyphen as punctuation; American spelling.
* Do not edit another module's files. Small changes to core files
  (aaa-utils.R, variables.R, equations.R, system.R, codegen.R, problems.R,
  solution.R, solvers-ode.R, callbacks.R, compile.R) are allowed when
  unavoidable; keep them minimal and list them in your final report.
* Never commit, push, or change git configuration.
