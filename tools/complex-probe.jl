# Reference behavior of upstream for complex and array states and for the
# interpolant derivatives, for tests/testthat/test-complex-states.R:
#  A. the Landau-Zener problem i psi' = H(t) psi of
#     OrdinaryDiffEq test/InterfaceI/complex_tests.jl (out of place) with the
#     explicit methods and the implicit methods with AutoFiniteDiff;
#  B. a scalar complex ODE and the complex SDE of StochasticDiffEq
#     test/complex_tests.jl on a fixed noise grid;
#  C. a matrix state (u0 a 2 x 2 matrix);
#  D. sol(t, Val{n}) of RadauIIA5, NordsieckBDF and Tsit5 and the linear
#     interpolation of a saveat solution outside its saved range;
#  E. NewtonRaphson on a complex u0.
# Every case runs in try/catch; an error prints its message.
# Run: julia --project=tools/julia --startup-file=no tools/complex-probe.jl
using OrdinaryDiffEq, StochasticDiffEq, DiffEqNoiseProcess, NonlinearSolve, LinearAlgebra
using OrdinaryDiffEqBDF, OrdinaryDiffEqFIRK, OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqTsit5
using OrdinaryDiffEqVerner, OrdinaryDiffEqSDIRK, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqDifferentiation: AutoFiniteDiff
import Pkg

for (k, v) in sort(collect(Pkg.dependencies()); by = x -> x[2].name)
    v.name in ("OrdinaryDiffEq", "OrdinaryDiffEqCore", "StochasticDiffEq", "NonlinearSolve",
        "OrdinaryDiffEqBDF", "OrdinaryDiffEqFIRK", "SciMLBase") &&
        println("version ", v.name, " ", v.version)
end

flat(x) = collect(Iterators.flatten(x))
cplx(x) = [real.(flat(x)); imag.(flat(x))]
out(name, key, x) = println(name, " ", key, " ", join(repr.(Float64.(x)), " "))
outs(name, key, s) = println(name, " ", key, " ", s)
function go(name, thunk)
    try
        sol = thunk()
        outs(name, "retcode", string(sol.retcode))
        out(name, "t", sol.t)
        out(name, "uend", cplx(sol.u[end]))
        st = sol.stats
        out(name, "stats", [st.nf, st.naccept, st.nreject])
        return sol
    catch e
        outs(name, "error", replace(first(sprint(showerror, e), 300), '\n' => ' '))
        return nothing
    end
end
function try_out(name, key, thunk)
    try
        out(name, key, thunk())
    catch e
        outs(name, key * "_error", replace(first(sprint(showerror, e), 300), '\n' => ' '))
    end
end

# A. Landau-Zener
H(t) = -im * [t 1; 1 -t]
fun(psi, p, t) = H(t) * psi
const T = 0.1
for alg in (Midpoint(), RK4(), DP5(), Tsit5(), Vern7(),
        ImplicitEuler(autodiff = AutoFiniteDiff()), Trapezoid(autodiff = AutoFiniteDiff()),
        Kvaerno3(autodiff = AutoFiniteDiff()), Rosenbrock23(autodiff = AutoFiniteDiff()))
    name = "LZ_" * string(nameof(typeof(alg)))
    sol = go(name, () -> solve(ODEProblem(fun, [1.0 + 0.0im; 0.0], (-T, T)), alg))
    sol === nothing && continue
    try_out(name, "solT", () -> cplx(sol(T)))
    try_out(name, "sol0", () -> cplx(sol(0.0)))
    try_out(name, "dsol0", () -> cplx(sol(0.0, Val{1})))
end

# B. scalar complex ODE (u' = (-1 + 2i) u) and the complex SDE on a grid
go("SC_Tsit5", () -> solve(ODEProblem((u, p, t) -> (-1 + 2im) * u, 1.0 + 1.0im, (0.0, 1.0)), Tsit5()))
go("SC_Tsit5_tol", () -> solve(ODEProblem((u, p, t) -> (-1 + 2im) * u, 1.0 + 1.0im, (0.0, 1.0)), Tsit5();
    abstol = 1e-10, reltol = 1e-10))
tg = collect(0.0:0.01:1.0)
Wg = [0.0; cumsum(0.1 .* sin.(1:100))]
sf(u, p, t) = 1.01u
try
    for alg in (SRIW1(), SOSRI(), SRA3(), SOSRA(), EM())
        name = "SDE_" * string(nameof(typeof(alg)))
        prob = SDEProblem(sf, sf, 1.0 + 1.0im, (0.0, 1.0); noise = NoiseGrid(tg, Wg))
        go(name, () -> solve(prob, alg; dt = 0.01, adaptive = false))
    end
catch e
    outs("SDE", "error", replace(first(sprint(showerror, e), 300), '\n' => ' '))
end

# C. matrix state
A2 = [-1.0 0.5; 0.0 -2.0]
go("MAT_Tsit5", () -> solve(ODEProblem((u, p, t) -> A2 * u, [1.0 2.0; 3.0 4.0], (0.0, 1.0)), Tsit5()))

# D. interpolant derivatives
rf(u, p, t) = [-0.5 * u[1] + sin(t), u[1] - 2.0 * u[2]]
for alg in (RadauIIA5(), NordsieckBDF(), Tsit5())
    name = "D_" * string(nameof(typeof(alg)))
    sol = go(name, () -> solve(ODEProblem(rf, [1.0, 0.5], (0.0, 2.0)), alg; abstol = 1e-8, reltol = 1e-8))
    sol === nothing && continue
    for n in 0:4
        try_out(name, "d$n", () -> sol(0.73, Val{n}))
    end
end
sol = go("LIN_saveat", () -> solve(ODEProblem(rf, [1.0, 0.5], (0.0, 2.0)), Tsit5(); saveat = 0.5))
if sol !== nothing
    for tq in (-0.25, 0.3, 2.4)
        try_out("LIN_saveat", "v$tq", () -> sol(tq))
        try_out("LIN_saveat", "dv$tq", () -> sol(tq, Val{1}))
    end
end

# E. NewtonRaphson on a complex u0
nf(u, p) = u .^ 2 .+ 1
for u0 in ([1.0 + 1.0im], [0.5 - 2.0im, 1.0 + 0.5im])
    name = "NR_" * string(length(u0))
    try
        s = solve(NonlinearProblem(nf, u0), NewtonRaphson(autodiff = AutoFiniteDiff()))
        outs(name, "retcode", string(s.retcode))
        out(name, "u", cplx(s.u))
        out(name, "stats", [s.stats.nf, s.stats.nsteps])
    catch e
        outs(name, "error", replace(first(sprint(showerror, e), 300), '\n' => ' '))
    end
end
