# Public names of the SciML packages modelcompR ports, for tools/ledger.R.
#
# Run from the repository root:
#   julia --project=tools/julia tools/upstream_exports.jl
# Writes tools/upstream_exports.tsv (package, name). ModelingToolkit itself
# is listed in tools/mtk_exports.tsv.
using ModelingToolkit, JumpProcesses, StochasticDiffEq, OrdinaryDiffEq, NonlinearSolve, Optimization, DelayDiffEq,
      ModelingToolkitStandardLibrary
# indirect dependencies, loaded through the packages that use them
const SciMLBase = OrdinaryDiffEq.SciMLBase
const DiffEqCallbacks = Base.loaded_modules[Base.PkgId(Base.UUID("459566f4-90b8-5000-8ac3-15dfb0a30def"), "DiffEqCallbacks")]
const DiffEqNoiseProcess = StochasticDiffEq.DiffEqNoiseProcess
const LinearSolve = NonlinearSolve.LinearSolve
const PKGS = [JumpProcesses, StochasticDiffEq, OrdinaryDiffEq, NonlinearSolve, Optimization, DelayDiffEq,
              SciMLBase, DiffEqCallbacks, DiffEqNoiseProcess, LinearSolve, ModelingToolkitStandardLibrary]
function public_names(m::Module)
    ns = Base.names(m; all = false)
    filter(n -> n != nameof(m) && !startswith(string(n), "#"), ns)
end
open(joinpath(@__DIR__, "upstream_exports.tsv"), "w") do io
    println(io, "package\tname")
    for m in PKGS
        for n in sort(string.(public_names(m)))
            println(io, string(nameof(m)), "\t", n)
        end
    end
end
