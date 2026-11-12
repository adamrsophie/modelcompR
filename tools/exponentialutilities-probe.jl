# Reference values of the ExponentialUtilities.jl kernels ported in
# R/expo-utilities.R (matrix exponentials, phi functions, Arnoldi and
# Lanczos, expv, phiv, phiv_timestep), for tests/testthat/test-expo-utilities.R.
#
# Run from the repository root:
#   julia --project=tools/julia tools/exponentialutilities-probe.jl OUT.txt
# Writes one line per value: a name, then the dimensions and the entries in
# column major order, all separated by spaces.
using OrdinaryDiffEqExponentialRK, LinearAlgebra
# ExponentialUtilities is reached through OrdinaryDiffEqExponentialRK, which
# loads it (it is not a direct dependency of tools/julia)
const ExponentialUtilities = OrdinaryDiffEqExponentialRK.ExponentialUtilities
using .ExponentialUtilities

const OUT = open(length(ARGS) > 0 ? ARGS[1] : "expo-probe.txt", "w")
emit(name, x::Number) = println(OUT, name, " 1 1 ", repr(Float64(x)))
emit(name, x::AbstractVector) = println(OUT, name, " ", length(x), " 1 ", join(repr.(Float64.(x)), " "))
emit(name, x::AbstractMatrix) = println(OUT, name, " ", size(x, 1), " ", size(x, 2), " ", join(repr.(Float64.(vec(x))), " "))

A3 = [-1.0 0.3 0.1; 0.2 -2.0 0.5; 0.05 0.4 -0.7]
S3 = [-2.0 0.5 0.1; 0.5 -1.0 0.3; 0.1 0.3 -0.5]
B3 = [1.0 1.0e5 0.0; 1.0e-5 1.0 2.0e3; 0.0 3.0e-4 -2.0]   # badly scaled
T3 = [1.0 2.0 3.0; 0.0 0.5 1.0; 0.0 0.0 -1.0]            # triangular (permutations)
P4 = [0.0 1.0 0.0 0.0; 2.0 0.0 0.5 0.0; 0.0 0.0 0.0 0.0; 1.0 0.0 0.3 -1.0]
mats = Dict("A3" => A3, "S3" => S3, "B3" => B3, "T3" => T3, "P4" => P4)
scales = [0.005, 0.1, 0.3, 0.6, 1.5, 10.0, 100.0]

for (nm, M) in sort(collect(mats); by = first)
    for (i, c) in enumerate(scales)
        X = c * M
        emit("gen_$(nm)_$i", exponential!(copy(X), ExpMethodGeneric()))
        emit("hb_$(nm)_$i", exponential!(copy(X), ExpMethodHigham2005Base()))
        emit("base_$(nm)_$i", exp(X))
    end
    for k in 1:4, (i, c) in enumerate([0.1, 1.0, 8.0])
        P = phi(c * M, k)
        for j in 0:k
            emit("phim_$(nm)_$(k)_$(i)_$j", P[j + 1])
        end
    end
    v = collect(1.0:size(M, 1))
    for k in 1:3
        emit("phivd_$(nm)_$k", ExponentialUtilities.phiv_dense(0.7 * M, v, k))
    end
end
for z in [-30.0, -2.5, -0.01, 0.0, 0.3, 4.0], k in 1:4
    emit("phis_$(z)_$k", phi(z, k))
end
D = Diagonal([-1.0, -0.2, 0.5])
for k in 1:3
    P = phi(0.4 * D, k)
    for j in 0:k
        emit("phidiag_$(k)_$j", Matrix(P[j + 1]))
    end
end

# Krylov
n = 15
A15 = [((i * 7 + j * 13) % 11 - 5) / 10 + (i == j ? -2.0 : 0.0) for i in 1:n, j in 1:n]
S15 = (A15 + A15') / 2
b15 = [1.0 + 0.1 * i for i in 1:n]
for (nm, M) in [("A15", A15), ("S15", S15)]
    for (m, iop) in [(5, 0), (8, 2), (15, 0)]
        Ks = arnoldi(M, b15; m = m, iop = iop)
        emit("arnV_$(nm)_$(m)_$iop", ExponentialUtilities.getV(Ks))
        emit("arnH_$(nm)_$(m)_$iop", ExponentialUtilities.getH(Ks))
        emit("arnm_$(nm)_$(m)_$iop", Ks.m)
        emit("expv_$(nm)_$(m)_$iop", expv(0.3, M, b15; m = m, iop = iop))
        w, err = phiv!(zeros(n, 4), 0.3, Ks, 3; errest = true)
        emit("phiv_$(nm)_$(m)_$iop", w)
        emit("phiverr_$(nm)_$(m)_$iop", err)
        wc = phiv!(zeros(n, 3), 0.3, Ks, 2; correct = true)
        emit("phivc_$(nm)_$(m)_$iop", wc)
    end
end
# happy breakdown: an invariant subspace
Kb = arnoldi(A3, [1.0, 0.0, 0.0]; m = 3)
emit("arnm_break", Kb.m)
emit("expv_break", expv(0.5, A3, [1.0, 0.0, 0.0]; m = 3))

# adaptive time stepping
B = hcat(b15, 0.5 * b15, [Float64(i % 3) for i in 1:n])
for (nm, M) in [("A15", A15), ("S15", S15)]
    for adaptive in (false, true), (m, iop) in [(5, 0), (8, 2)]
        tag = "$(nm)_$(adaptive)_$(m)_$iop"
        emit("pts_opn_$tag", phiv_timestep(2.0, M, B; m = m, iop = iop, adaptive = adaptive, tol = 1.0e-8, opnorm = opnorm))
        emit("pts_arn_$tag", phiv_timestep(2.0, M, B; m = m, iop = iop, adaptive = adaptive, tol = 1.0e-8))
        emit("pts_tau_$tag", phiv_timestep([0.5, 1.2, 2.0], M, B; m = m, iop = iop, adaptive = adaptive, tol = 1.0e-6,
                                           opnorm = opnorm, tau = 0.4))
        emit("ets_$tag", expv_timestep(1.5, M, b15; m = m, iop = iop, adaptive = adaptive, tol = 1.0e-9))
    end
end
close(OUT)
