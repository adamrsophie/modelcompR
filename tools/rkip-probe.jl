# Reference solutions of the OrdinaryDiffEqRKIP method ported in
# R/solvers-rkip.R, for tests/testthat/test-rkip.R, with the Verner6 tableau
# and the logrange step grids.
#
# Run from the repository root:
#   julia --project=tools/julia tools/rkip-probe.jl OUT.txt
# Writes one line per value: a name, then the dimensions and the entries in
# column major order, all separated by spaces. For each solve the saved
# states (one column per time), the times and the statistics
# (nf, nf2, naccept, nreject) are written.
#
# RKIP needs exp(A, t) to return an operator of the same type, which a
# MatrixOperator does not do, so the probe defines a small operator whose
# exp(A, t) is the operator of LinearAlgebra.exp(t * A) (what the R port
# computes for a MatrixOperator).
using OrdinaryDiffEqRKIP, SciMLOperators, SciMLBase, LinearAlgebra
using OrdinaryDiffEqRKIP: Verner6

const OUT = open(length(ARGS) > 0 ? ARGS[1] : "rkip-probe.txt", "w")
emit(name, x::Number) = println(OUT, name, " 1 1 ", repr(Float64(x)))
emit(name, x::AbstractVector) = println(OUT, name, " ", length(x), " 1 ", join(repr.(Float64.(x)), " "))
emit(name, x::AbstractMatrix) = println(OUT, name, " ", size(x, 1), " ", size(x, 2), " ", join(repr.(Float64.(vec(x))), " "))
function emitsol(name, sol)
    emit("u_$name", reduce(hcat, [vec(u) for u in sol.u]))
    emit("t_$name", sol.t)
    return emit("s_$name", [sol.stats.nf, sol.stats.nf2, sol.stats.naccept, sol.stats.nreject])
end
function tryemit(name, thunk)
    return try
        emitsol(name, thunk())
    catch err
        println(OUT, "err_$name 1 1 1.0")
        println(stderr, name, ": ", sprint(showerror, err))
    end
end

struct ExpOp{T} <: SciMLOperators.AbstractSciMLOperator{T}
    A::Matrix{T}
end
Base.size(L::ExpOp) = size(L.A)
Base.size(L::ExpOp, i) = size(L.A, i)
LinearAlgebra.exp(L::ExpOp, t) = ExpOp(exp(t * L.A))
(L::ExpOp)(v, u, p, t) = L.A * v
(L::ExpOp)(w, v, u, p, t) = mul!(w, L.A, v)
SciMLOperators.isconstant(::ExpOp) = true
SciMLOperators.islinear(::ExpOp) = true
SciMLOperators.update_coefficients(L::ExpOp, u, p, t) = L
SciMLOperators.update_coefficients!(L::ExpOp, u, p, t) = nothing
Base.convert(::Type{AbstractMatrix}, L::ExpOp) = L.A
Base.:*(L::ExpOp, v::AbstractVector) = L.A * v

tab = Verner6(Float64)
emit("tab_A", tab.A)
emit("tab_c", tab.c)
emit("tab_alpha", tab.α)
emit("tab_alphaEEst", tab.αEEst)
emit("tab_orders", [tab.stages, tab.order, tab.adaptiveorder])
emit("grid_default", collect(logrange(1.0e-3, 1.0, 100)))
emit("grid_b", collect(logrange(1.0e-2, 10.0, 100)))
emit("grid_c", collect(logrange(1.0e-4, 0.05, 37)))

A4 = [
    -2.0 0.5 0.1 0.0;
    0.3 -1.5 0.4 0.2;
    0.0 0.2 -1.0 0.5;
    0.1 0.0 0.3 -8.0
]
# time independent nonlinear part (the stage times are formed under @fastmath)
f2(u, p, t) = [-u[1] * u[1] * u[1] + 0.2 * u[2] * u[3], 0.5 * u[1] * u[2], -0.2 * u[3] * u[3] + 0.3 * u[4], 0.1 * u[4] * u[1]]
u0 = [1.0, 0.5, -0.5, 0.2]
mk(tspan) = SplitODEProblem(SplitFunction{false}(ExpOp(A4), f2), u0, tspan)
prob = mk((0.0, 2.0))
emit("isinplace", isinplace(prob) ? 1.0 : 0.0)

tryemit("default", () -> solve(prob, RKIP()))
tryemit("tight", () -> solve(prob, RKIP(); reltol = 1.0e-8, abstol = 1.0e-10))
tryemit("gridb", () -> solve(prob, RKIP(1.0e-2, 10.0); reltol = 1.0e-6, abstol = 1.0e-8))
tryemit("clamplow", () -> solve(prob, RKIP(clamp_lower_dt = true); reltol = 1.0e-9, abstol = 1.0e-11))
tryemit("noclamphigh", () -> solve(mk((0.0, 10.0)), RKIP(1.0e-4, 0.05, nb_of_cache_step = 37, clamp_higher_dt = false)))
tryemit("clamphigh", () -> solve(mk((0.0, 10.0)), RKIP(1.0e-4, 0.05, nb_of_cache_step = 37)))
tryemit("fixed", () -> solve(prob, RKIP(); dt = 0.1, adaptive = false))
tryemit("fixedoff", () -> solve(prob, RKIP(); dt = 0.15, adaptive = false))
# one algorithm object for two solves: the second reuses the cached exponentials
alg = RKIP()
tryemit("reuse1", () -> solve(prob, alg; reltol = 1.0e-5, abstol = 1.0e-7))
tryemit("reuse2", () -> solve(mk((0.0, 3.0)), alg; reltol = 1.0e-6, abstol = 1.0e-8))
close(OUT)
