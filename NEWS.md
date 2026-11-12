# modelcompR 0.1.0

* First version: a native R port of ModelingToolkit.jl 11.47 on top of
  symcompR. All 327 public names of ModelingToolkit and ModelingToolkitBase
  are accounted for (318 ported, 9 not applicable). The port of the SciML
  packages it builds on (OrdinaryDiffEq, StochasticDiffEq, JumpProcesses,
  NonlinearSolve, LinearSolve, Optimization, DelayDiffEq, DiffEqCallbacks,
  SciMLBase, SciMLSensitivity) is in progress; `tools/upstream_status.tsv`
  records the status of every name.
