# Reference solutions of the OrdinaryDiffEqLinear methods ported in
# R/solvers-linear-lie.R, for tests/testthat/test-linear-lie.R.
#
# Run from the repository root:
#   julia --project=tools/julia tools/linear-probe.jl OUT.txt
# Writes one line per value: a name, then the dimensions and the entries in
# column major order, all separated by spaces. For each solve the saved
# states (one column per time), the times and the statistics
# (nf, naccept, nreject) are written.
using OrdinaryDiffEqLinear, SciMLOperators, LinearAlgebra

const OUT = open(length(ARGS) > 0 ? ARGS[1] : "linear-probe.txt", "w")
emit(name, x::Number) = println(OUT, name, " 1 1 ", repr(Float64(x)))
emit(name, x::AbstractVector) = println(OUT, name, " ", length(x), " 1 ", join(repr.(Float64.(x)), " "))
emit(name, x::AbstractMatrix) = println(OUT, name, " ", size(x, 1), " ", size(x, 2), " ", join(repr.(Float64.(vec(x))), " "))
function emitsol(name, sol)
    emit("u_$name", reduce(hcat, [vec(u) for u in sol.u]))
    emit("t_$name", sol.t)
    return emit("s_$name", [sol.stats.nf, sol.stats.naccept, sol.stats.nreject])
end

# A(u, t) for a 3 x 3 state and time dependent problem
function upd1!(A, u, p, t)
    A[1, 1] = -0.5; A[1, 2] = 1.0 + 0.2 * t; A[1, 3] = 0.1 * u[1] * u[2]
    A[2, 1] = -1.0; A[2, 2] = -0.3; A[2, 3] = 0.5 * t
    A[3, 1] = 0.2 * u[3] * u[3]; A[3, 2] = -0.4; A[3, 3] = -0.1 - 0.3 * t * t
    return nothing
end
L1 = MatrixOperator(zeros(3, 3); update_func! = upd1!)
upd1!(L1.A, [1.0, 0.5, -0.3], nothing, 0.0)
prob1 = ODEProblem(L1, [1.0, 0.5, -0.3], (0.0, 1.0))
emit("isinplace1", isinplace(prob1) ? 1.0 : 0.0)

# 6 x 6 nonsymmetric and symmetric problems for the Krylov variants
base6 = [((3i + 5j) % 7 - 3) / 10 - (i == j ? 1.0 : 0.0) for i in 1:6, j in 1:6]
sym6 = (base6 + base6') / 2
function upd2!(A, u, p, t)
    A .= base6
    for i in 1:6
        A[i, i] -= 0.2 * t
    end
    A[1, 6] += 0.1 * u[1] * u[1]
    return nothing
end
function upd3!(A, u, p, t)
    A .= sym6
    for i in 1:6
        A[i, i] -= 0.2 * t
    end
    A[1, 6] += 0.05 * u[1] * u[1]
    A[6, 1] += 0.05 * u[1] * u[1]
    return nothing
end
u06 = [1.0 - 0.1 * (i - 1) for i in 1:6]
L2 = MatrixOperator(copy(base6); update_func! = upd2!)
L3 = MatrixOperator(copy(sym6); update_func! = upd3!)
prob2 = ODEProblem(L2, u06, (0.0, 1.0))
prob3 = ODEProblem(L3, u06, (0.0, 1.0))

algs = [
    ("MagnusMidpoint", MagnusMidpoint), ("MagnusLeapfrog", MagnusLeapfrog), ("LieEuler", LieEuler),
    ("MagnusGauss4", MagnusGauss4), ("MagnusNC6", MagnusNC6), ("MagnusGL6", MagnusGL6),
    ("MagnusGL8", MagnusGL8), ("MagnusNC8", MagnusNC8), ("MagnusGL4", MagnusGL4),
    ("RKMK2", RKMK2), ("RKMK4", RKMK4), ("LieRK4", LieRK4), ("CG2", CG2), ("CG3", CG3), ("CG4a", CG4a),
]
krylov_algs = ["MagnusMidpoint", "MagnusLeapfrog", "LieEuler", "MagnusGauss4", "MagnusNC6", "MagnusGL6",
    "MagnusGL8", "MagnusNC8", "MagnusGL4", "RKMK2", "RKMK4", "LieRK4"]
for (nm, A) in algs
    emitsol("p1_$nm", solve(prob1, A(); dt = 0.1))
    emitsol("p1dt_$nm", solve(prob1, A(); dt = 0.25))
end
for (nm, A) in algs
    nm in krylov_algs || continue
    emitsol("p1k_$nm", solve(prob1, A(krylov = true, m = 2, iop = 1); dt = 0.1))
    emitsol("p2k_$nm", solve(prob2, A(krylov = true, m = 4, iop = 2); dt = 0.1))
    emitsol("p2kd_$nm", solve(prob2, A(krylov = true); dt = 0.2))
    emitsol("p3k_$nm", solve(prob3, A(krylov = true, m = 4); dt = 0.1))
    emitsol("p2_$nm", solve(prob2, A(); dt = 0.1))
end

# adaptive MagnusAdapt4
prob1l = ODEProblem(L1, [1.0, 0.5, -0.3], (0.0, 3.0))
emitsol("a1_default", solve(prob1l, MagnusAdapt4()))
emitsol("a1_tight", solve(prob1l, MagnusAdapt4(); reltol = 1.0e-7, abstol = 1.0e-9))
emitsol("a2_default", solve(ODEProblem(L2, u06, (0.0, 2.0)), MagnusAdapt4(); reltol = 1.0e-5, abstol = 1.0e-7))
emitsol("a1_fixed", solve(prob1, MagnusAdapt4(); dt = 0.1, adaptive = false))

# LinearExponential on a constant operator
LC = MatrixOperator(copy(base6))
probc = ODEProblem(LC, u06, (0.0, 2.0))
emit("isinplacec", isinplace(probc) ? 1.0 : 0.0)
emitsol("le_nodt", solve(probc, LinearExponential()))
emitsol("le_off", solve(probc, LinearExponential(); dt = 0.5))
emitsol("le_simple", solve(probc, LinearExponential(krylov = :simple); dt = 0.5))
emitsol("le_simple4", solve(probc, LinearExponential(krylov = :simple, m = 4, iop = 2); dt = 0.5))
emitsol("le_adapt", solve(probc, LinearExponential(krylov = :adaptive); dt = 0.5))
emitsol("le_adapt3", solve(probc, LinearExponential(krylov = :adaptive, m = 3); dt = 0.5, reltol = 1.0e-8))
emitsol("le_adapt3n", solve(probc, LinearExponential(krylov = :adaptive, m = 3, iop = 1)))
LS = MatrixOperator(copy(sym6))
emitsol("le_sym", solve(ODEProblem(LS, u06, (0.0, 2.0)), LinearExponential(krylov = :simple, m = 4); dt = 0.5))

# CayleyEuler on a matrix state
function updc!(A, U, p, t)
    A[1, 1] = 0.0; A[1, 2] = 1.0 + 0.1 * t; A[1, 3] = 0.2 * U[1, 1] * U[2, 2]
    A[2, 1] = -1.0 - 0.1 * t; A[2, 2] = 0.0; A[2, 3] = 0.3
    A[3, 1] = -0.2 * U[1, 1] * U[2, 2]; A[3, 2] = -0.3; A[3, 3] = -0.05 * U[3, 3] * U[3, 3]
    return nothing
end
U0 = [1.0 0.2 0.0; 0.1 0.9 0.3; 0.0 -0.2 1.1]
LCa = MatrixOperator(zeros(3, 3); update_func! = updc!)
updc!(LCa.A, U0, nothing, 0.0)
probca = ODEProblem(LCa, copy(U0), (0.0, 1.0))
emit("isinplaceca", isinplace(probca) ? 1.0 : 0.0)
emitsol("cay", solve(probca, CayleyEuler(); dt = 0.1))
f2!(du, u, p, t) = (du .= 0.1 .* u .* u; nothing)
probcs = SplitODEProblem(LCa, f2!, copy(U0), (0.0, 1.0))
emit("isinplacecs", isinplace(probcs) ? 1.0 : 0.0)
emitsol("cays", solve(probcs, CayleyEuler(); dt = 0.25))
close(OUT)
