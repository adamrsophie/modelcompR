# modelcompR: Symbolic Numeric Modeling in R

<!-- badges: start -->
[![R-CMD-check](https://github.com/adamrsophie/modelcompR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/adamrsophie/modelcompR/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

***modelcompR*** is a native R port of the Julia library
[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) (version
11.47). Models are written as symbolic equations: variables, parameters,
derivatives and hierarchical components. `mtkcompile()` simplifies them
structurally, and modelcompR generates R functions for the numerical solvers.
The symbolic layer is [symcompR](https://github.com/adamrsophie/symcompR), a
native R port of SymPy. Everything runs in R: there is no Julia, no Python and
no compiled code.

## Installation

```r
# install.packages("remotes")
remotes::install_github("adamrsophie/symcompR")
remotes::install_github("adamrsophie/modelcompR")
```

## Quick start

The first order lag from the ModelingToolkit tutorial:

```r
library(modelcompR)
t <- t_nounits
D <- D_nounits

variables("x(t)")
parameters("tau")

sys <- System(list(D(x) ~ (1 - x) / tau), t, name = "lag")
sys <- mtkcompile(sys)

prob <- ODEProblem(sys, list(x ~ 0, tau ~ 3), c(0, 10))
sol <- solve(prob)

sol[x]                 # time series of x
sol(5)                 # dense interpolation at t = 5
sol[(1 - x) / tau]     # any expression of the model
plot(sol)
```

Equations are R formulas inside `list()`, so they read as in Julia:
`D(x) ~ -k * x`. `Equation(lhs, rhs)` builds the same object without a formula.
Variables are declared with `variables("x(t) = 1")` and
`parameters("k = 0.5")`, which mirror `@variables` and `@parameters`, or one
at a time with `Variable()` and `Parameter()`.

A second order equation with an observed variable:

```r
variables("y(t)", "v(t)")
parameters("k = 2")
osc <- mtkcompile(System(list(D(D(y)) ~ -k * y, v ~ 2 * D(y)), t, name = "osc"))
unknowns(osc)   # y(t), yˍt(t): the derivative is lowered as in ModelingToolkit
observed(osc)   # v(t) ~ 2*yˍt(t)
sol <- solve(ODEProblem(osc, list(y ~ 1, D(y) ~ 0), c(0, 5)))
sol[v]
```

## Relation to ModelingToolkit.jl

modelcompR follows ModelingToolkit 11: one `System` type, `mtkcompile()`,
bindings and initial conditions, the problem constructors, symbolic indexing of
problems and solutions, and the same names wherever R allows them. Where Julia
syntax has no R equivalent:

| Julia | R |
|-------|---|
| `@variables x(t) = 1` | `variables("x(t) = 1")` or `x <- Variable("x", t, default = 1)` |
| `@parameters k` | `parameters("k")` or `k <- Parameter("k")` |
| `D(x) ~ -k * x` | the formula `D(x) ~ -k * x`, or `Equation(D(x), -k * x)` |
| `@named sys = System(...)` | `sys <- System(..., name = "sys")` |
| `@mtkmodel`, `@component` | R functions that return a `System` |
| `sys.x`, `sys.sub.x` | `sys$x`, `sys$sub$x` |
| `sol[x]`, `sol(t)` | `sol[x]`, `sol(t)` |
| `ModelingToolkit.t_nounits` | `t_nounits` (the name `t` would mask `base::t`) |

`solve()` is a method of the base R generic, the same pattern symcompR uses,
so `solve()` keeps working on matrices and symcompR expressions. The numerical
solvers are pure R ports of the corresponding OrdinaryDiffEq.jl,
NonlinearSolve.jl, StochasticDiffEq.jl and JumpProcesses.jl methods and keep
their names (`Tsit5()`, `Rodas5P()`, `FBDF()`, `NewtonRaphson()`, ...).

## What is covered

All 327 public names of ModelingToolkit and ModelingToolkitBase are accounted
for: 310 are ported, 8 are the Julia dynamic optimization backends (JuMP,
InfiniteOpt, CasADi, Pyomo), which point to the native collocation backend,
and 9 are Julia internals with no R meaning. Some ported functions cover
fewer cases than their Julia counterparts; [PORTING.md](PORTING.md) lists the
module by module status, what is not yet ported and the deliberate
differences.

| Area | Highlights |
|------|------------|
| Modeling | `System`, variables and parameters with metadata, bindings and initial conditions, scopes, `compose`, `extend`, connectors with flow and stream variables, analysis points |
| Compilation | `mtkcompile`: connection expansion, alias elimination, Pantelides index reduction, dummy derivative state selection, tearing, observed equations |
| Problems | ODE (with mass matrix DAEs), DAE, nonlinear, steady state, interval, linear, optimization, SDE, jump, delay, discrete, implicit discrete, boundary value, initialization, linearization, dynamic optimization |
| Solvers | Tsit5, Vern7, Vern9, DP5, BS3, Rosenbrock23, Rodas4, Rodas5P, FBDF, DFBDF, automatic stiffness switching; NewtonRaphson, TrustRegion, LevenbergMarquardt, homotopy continuation; EM, SRIW1, SOSRI; Direct SSA; MethodOfSteps; MIRK4 and shooting |
| Events | Continuous and discrete symbolic events with `Pre`, imperative affects, clocks, `Sample`, `Hold`, `Shift` |
| Analysis | `linearize`, sensitivity and loop transfer functions, `modelingtoolkitize`, variable and independent variable changes |

## Validation

The tests compare modelcompR with ModelingToolkit itself: `tools/julia/reference.jl`
writes reference trajectories and compiled system structures (unknowns,
observed variables and parameters) to `tests/testthat/fixtures/julia/`, and the
R tests check the RC circuit, the Cartesian pendulum, the Robertson problem, an
index 1 DAE, the Lotka Volterra model and a nonlinear system against them.
Other tests use analytic solutions and convergence orders.

## Compatibility with symcompR

Expressions in modelcompR are symcompR expressions, so every symcompR function
(`simplify`, `expand`, `latex`, `diff`, `subs`, `lambdify`, ...) works on
equations, Jacobians and observed expressions. Symbols and functions created
with symcompR directly can be used in models; as in ModelingToolkit, undeclared
symbols are unknowns. Attaching both packages masks nothing: `flatten()` and
`compose()` hand any argument that is not a system to symcompR, and no
modelcompR function masks a base R function.

## License and credits

modelcompR is distributed under the GNU General Public License, version 3. It
ports algorithms from ModelingToolkit.jl, OrdinaryDiffEq.jl and the other
SciML packages (MIT License), and from StateSelection.jl, whose recent code is
licensed under the GNU Affero General Public License, version 3; the R files
derived from it carry that license. See `inst/COPYRIGHTS` for details and the
retained license notices.
