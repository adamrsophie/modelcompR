# Coefficients and reference runs of the low storage Runge-Kutta methods of
# OrdinaryDiffEqLowStorageRK (the out of place ConstantCache forms of
# low_storage_rk_perform_step.jl), for R/solvers-lowstorage.R and
# tests/testthat/test-lowstorage.R.
#
# Run from the repository root:
#   julia --project=tools/julia tools/lowstorage-probe.jl CACHES.tsv RUNS.tsv
# CACHES.tsv: method, cache type, adaptive, order, then field=value pairs
# (hex floats, vectors comma separated; CACHES.tsv - skips them). RUNS.tsv: method, case, steps, nf,
# naccept, nreject, final u, times (hex).
using OrdinaryDiffEqLowStorageRK, Printf
const OC = OrdinaryDiffEqLowStorageRK.OrdinaryDiffEqCore

hx(x::Real) = @sprintf("%a", Float64(x))
hx(x::AbstractVector) = join(hx.(x), ",")
hx(x::Tuple) = join(hx.(collect(x)), ",")

algs = sort([n for n in names(OrdinaryDiffEqLowStorageRK)
             if isdefined(OrdinaryDiffEqLowStorageRK, n) &&
                (T = getfield(OrdinaryDiffEqLowStorageRK, n); T isa Type && T <: OC.OrdinaryDiffEqAlgorithm)])

# polynomial right hand sides: the elementary functions of Julia and of the C
# library may differ in the last bit
f(u, p, t) = [u[2], -u[1] + 0.3 * (1 - t / 3) * u[1]^2]
g(u, p, t) = [-20 * (u[1] - t * (2 - t)), 0.5 * u[1] * u[2]]
P1 = ODEProblem(f, [1.0, 0.5], (0.0, 3.0))
P2 = ODEProblem(g, [0.0, 1.0], (0.0, 2.0))

ARGS[1] != "-" && open(ARGS[1], "w") do io
    for n in algs
        alg = getfield(OrdinaryDiffEqLowStorageRK, n)()
        integ = init(P1, alg; dt = 0.01, adaptive = false)
        c = integ.cache
        fields = String[]
        for fn in fieldnames(typeof(c))
            v = getfield(c, fn)
            if v isa Real || (v isa Union{AbstractVector, Tuple} && all(x -> x isa Real, v))
                push!(fields, string(fn, "=", hx(v)))
            end
        end
        ctrl = OC.isadaptive(alg) ? OC.default_controller(Float64, alg) : nothing
        beta = ctrl isa OC.PIDController ? hx(collect(ctrl.beta)) : ""
        println(io, join(vcat([string(n), string(nameof(typeof(c))), string(OC.isadaptive(alg)),
                               string(OC.alg_order(alg)), string(OC.isfsal(alg)), "pid=" * beta], fields), "\t"))
    end
end

function row(io, n, case, sol)
    s = sol.stats
    println(io, join([string(n), case, string(length(sol.t)), string(s.nf), string(s.naccept),
                      string(s.nreject), hx(sol.u[end]), hx(sol.t)], "\t"))
end

open(ARGS[2], "w") do io
    for n in algs
        A = getfield(OrdinaryDiffEqLowStorageRK, n)
        alg = A()
        row(io, n, "fixed", solve(P1, alg; dt = 0.05, adaptive = false))
        row(io, n, "fixed_saveat", solve(P1, alg; dt = 0.05, adaptive = false, saveat = 0.37))
        if OC.isadaptive(alg)
            row(io, n, "adaptive", solve(P1, alg; reltol = 1e-6, abstol = 1e-8))
            row(io, n, "reject", solve(P2, alg; reltol = 1e-4, abstol = 1e-6, dt = 0.5))
            row(io, n, "saveat", solve(P1, alg; reltol = 1e-6, abstol = 1e-8, saveat = 0.37))
        end
    end
end
