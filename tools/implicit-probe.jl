# Out of place ImplicitEuler, Trapezoid and FBDF on a nonlinear heat equation;
# the tab separated output feeds tools/implicit-gen-test.R
# (tests/testthat/test-implicit-julia.R). Run from the repository root:
#   julia --project=tools/julia tools/implicit-probe.jl > implicit.tsv
using OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF, OrdinaryDiffEqCore, SciMLBase
const N = 8
heat(u, p, t) = [p[1] * ((i == 1 ? zero(u[1]) : u[i-1]) - 2u[i] + (i == N ? zero(u[1]) : u[i+1])) * (N + 1)^2 - 0.5 * u[i]^3 + 0.1 * t for i in 1:N]
function hjac(u, p, t)
    J = zeros(N, N)
    for i in 1:N
        J[i, i] = -2 * p[1] * (N + 1)^2 - 1.5 * u[i]^2
        i > 1 && (J[i, i-1] = p[1] * (N + 1)^2)
        i < N && (J[i, i+1] = p[1] * (N + 1)^2)
    end
    J
end
u0 = [i * (N + 1 - i) / N^2 for i in 1:N]
prob = ODEProblem(ODEFunction(heat; jac = hjac), u0, (0.0, 0.3), [1.0])
cb = DiscreteCallback((u, t, integ) -> t == 0.1, integ -> (integ.u = integ.u * 0.5))
for (nm, a) in [("ImplicitEuler", ImplicitEuler()), ("Trapezoid", Trapezoid()), ("FBDF", FBDF())], sc in ("plain", "tstop", "cb")
    s = sc == "plain" ? solve(prob, a; reltol = 1e-4, abstol = 1e-6) :
        sc == "tstop" ? solve(prob, a; reltol = 1e-4, abstol = 1e-6, tstops = [0.1]) :
        solve(prob, a; reltol = 1e-4, abstol = 1e-6, tstops = [0.1], callback = cb)
    st = s.stats
    println(nm, "\t", sc, "\t", length(s.t), "\t", st.nf, "\t", st.njacs, "\t", st.nw, "\t", st.nsolve, "\t", st.nreject, "\t", join(repr.(s.u[end]), ","))
end
