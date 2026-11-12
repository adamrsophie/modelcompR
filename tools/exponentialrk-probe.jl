# Reference solutions of the OrdinaryDiffEqExponentialRK methods ported in
# R/solvers-exponential.R, for tests/testthat/test-exponential-rk.R. The
# problems are out of place, so the ConstantCache steps run.
#
# Run from the repository root:
#   julia --project=tools/julia tools/exponentialrk-probe.jl OUT.txt
# Writes one line per value: a name, then the dimensions and the entries in
# column major order, all separated by spaces. For each solve the saved
# states (one column per time), the times and the statistics
# (nf, nf2, njacs, naccept, nreject) are written.
using OrdinaryDiffEqExponentialRK, SciMLOperators, LinearAlgebra
using OrdinaryDiffEqDifferentiation: AutoFiniteDiff

const OUT = open(length(ARGS) > 0 ? ARGS[1] : "exponentialrk-probe.txt", "w")
emit(name, x::Number) = println(OUT, name, " 1 1 ", repr(Float64(x)))
emit(name, x::AbstractVector) = println(OUT, name, " ", length(x), " 1 ", join(repr.(Float64.(x)), " "))
emit(name, x::AbstractMatrix) = println(OUT, name, " ", size(x, 1), " ", size(x, 2), " ", join(repr.(Float64.(vec(x))), " "))
function emitsol(name, sol)
    emit("u_$name", reduce(hcat, [vec(u) for u in sol.u]))
    emit("t_$name", sol.t)
    return emit("s_$name", [sol.stats.nf, sol.stats.nf2, sol.stats.njacs, sol.stats.naccept, sol.stats.nreject])
end
function tryemit(name, thunk)
    return try
        emitsol(name, thunk())
    catch err
        println(OUT, "err_$name 1 1 1.0")
        println(stderr, name, ": ", sprint(showerror, err))
    end
end

# the literal coefficients of EPIRK5P1 and EPIRK5P2
emit(
    "consts", [
        0.35129592695058193092, 0.84405472011657126298, 0.71111095364366870359,
        0.62378111953371494809, 2.37739153404186675585, 1.78975267532362337976, 9.35854650579261718128,
        0.46629408528088195806, 0.8821791265336386514, 0.92074916488140031449, 0.79791561832664517267,
        2.80620373289331259751, 2.52806310256246280783, -0.128486782657005566142,
        -0.161028033172809183351, 5.26726331616909606251,
    ]
)

A4 = [
    -2.0 0.5 0.1 0.0;
    0.3 -1.5 0.4 0.2;
    0.0 0.2 -1.0 0.5;
    0.1 0.0 0.3 -2.5
]
S4 = (A4 + A4') / 2
f2(u, p, t) = [
    -u[1] * u[1] * u[1] + 0.1 * t, 0.5 * u[1] * u[2], -0.2 * u[3] * u[3] + 0.3 * u[4],
    0.1 * u[4] * u[1] - 0.05 * t * t,
]
jf2(u, p, t) = [
    -3 * u[1] * u[1] 0.0 0.0 0.0;
    0.5 * u[2] 0.5 * u[1] 0.0 0.0;
    0.0 0.0 -0.4 * u[3] 0.3;
    0.1 * u[4] 0.0 0.0 0.1 * u[1]
]
u0 = [1.0, 0.5, -0.5, 0.2]
sp1 = SplitODEProblem(MatrixOperator(A4), f2, u0, (0.0, 1.0))
sp2 = SplitODEProblem(MatrixOperator(S4), f2, u0, (0.0, 1.0))
sp3 = SplitODEProblem(DiagonalOperator([-1.0, -2.0, -0.5, -3.0]), f2, u0, (0.0, 1.0))
emit("isinplace_sp1", isinplace(sp1) ? 1.0 : 0.0)
fn(u, p, t) = A4 * u + f2(u, p, t)
jn(u, p, t) = A4 + jf2(u, p, t)
np1 = ODEProblem(ODEFunction(fn; jac = jn), u0, (0.0, 1.0))
emit("isinplace_np1", isinplace(np1) ? 1.0 : 0.0)

classical = [
    ("LawsonEuler", LawsonEuler), ("NorsettEuler", NorsettEuler), ("ETDRK2", ETDRK2),
    ("ETDRK3", ETDRK3), ("ETDRK4", ETDRK4), ("HochOst4", HochOst4), ("Friedli", Friedli),
]
for (nm, A) in classical
    tryemit("sp1_$nm", () -> solve(sp1, A(); dt = 0.1))
    tryemit("sp3_$nm", () -> solve(sp3, A(); dt = 0.1))
    tryemit("sp1k_$nm", () -> solve(sp1, A(krylov = true, m = 3, iop = 1); dt = 0.1))
    tryemit("sp1kd_$nm", () -> solve(sp1, A(krylov = true); dt = 0.25))
    tryemit("sp2k_$nm", () -> solve(sp2, A(krylov = true, m = 3); dt = 0.1))
    tryemit("np1k_$nm", () -> solve(np1, A(krylov = true, m = 3); dt = 0.1))
end
tryemit("sp1_ETD2", () -> solve(sp1, ETD2(); dt = 0.1))
tryemit("sp3_ETD2", () -> solve(sp3, ETD2(); dt = 0.25))
# a dt that does not divide the interval (the cached methods cannot change dt)
tryemit("sp1_LawsonEuler_dt03", () -> solve(sp1, LawsonEuler(); dt = 0.3))
tryemit("sp1_ETDRK4_dt03", () -> solve(sp1, ETDRK4(krylov = true); dt = 0.3))

epirk = [
    ("Exp4", Exp4), ("EPIRK4s3A", EPIRK4s3A), ("EPIRK4s3B", EPIRK4s3B), ("EPIRK5s3", EPIRK5s3),
    ("EXPRB53s3", EXPRB53s3), ("EPIRK5P1", EPIRK5P1), ("EPIRK5P2", EPIRK5P2),
]
for (nm, A) in epirk
    tryemit("np1_$nm", () -> solve(np1, A(); dt = 0.1))
    tryemit("np1t_$nm", () -> solve(np1, A(m = 3, iop = 1); dt = 0.1, reltol = 1.0e-8))
    tryemit("np1n_$nm", () -> solve(np1, A(adaptive_krylov = false, m = 3); dt = 0.25))
    tryemit("sp1e_$nm", () -> solve(sp1, A(); dt = 0.1))
end
for (nm, A) in [("Exprb32", Exprb32), ("Exprb43", Exprb43)]
    tryemit("np1a_$nm", () -> solve(np1, A()))
    tryemit("np1at_$nm", () -> solve(np1, A(m = 3); reltol = 1.0e-7, abstol = 1.0e-9))
    tryemit("np1af_$nm", () -> solve(np1, A(); dt = 0.1, adaptive = false))
    tryemit("sp1a_$nm", () -> solve(sp1, A()))
end

# no Jacobian: forward finite differences, the default of the R port
npfd = ODEProblem(fn, u0, (0.0, 1.0))
fd = AutoFiniteDiff()
tryemit("npfdk_ETDRK4", () -> solve(npfd, ETDRK4(krylov = true, m = 3, autodiff = fd); dt = 0.1))
tryemit("npfdk_LawsonEuler", () -> solve(npfd, LawsonEuler(krylov = true, autodiff = fd); dt = 0.25))
for (nm, A) in [("EPIRK4s3A", EPIRK4s3A), ("EPIRK5P2", EPIRK5P2)]
    tryemit("npfd_$nm", () -> solve(npfd, A(autodiff = fd); dt = 0.1))
    tryemit("sp1fd_$nm", () -> solve(sp1, A(autodiff = fd); dt = 0.1))
end
for (nm, A) in [("Exprb32", Exprb32), ("Exprb43", Exprb43)]
    tryemit("npfda_$nm", () -> solve(npfd, A(autodiff = fd)))
    tryemit("sp1fda_$nm", () -> solve(sp1, A(autodiff = fd)))
end
# a split operator with an update function: out of place, the steps use the
# operator as constructed and only the right-hand side updates it
upd(A, u, p, t) = A4 + t * [0.1 0.0 0.0 0.0; 0.0 -0.2 0.0 0.0; 0.0 0.0 0.3 0.0; 0.0 0.0 0.0 -0.1]
spu = SplitODEProblem(MatrixOperator(copy(A4); update_func = upd), f2, u0, (0.0, 1.0))
tryemit("spu_ETDRK4", () -> solve(spu, ETDRK4(); dt = 0.1))
tryemit("spuk_ETDRK4", () -> solve(spu, ETDRK4(krylov = true, m = 3); dt = 0.1))
tryemit("spuk_LawsonEuler", () -> solve(spu, LawsonEuler(krylov = true); dt = 0.25))
tryemit("spu_EPIRK4s3A", () -> solve(spu, EPIRK4s3A(autodiff = fd); dt = 0.1))
close(OUT)
