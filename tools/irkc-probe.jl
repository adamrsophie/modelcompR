# Reference values for tests/testthat/test-stabilized-irk.R (IRKC of
# OrdinaryDiffEqStabilizedIRK on a split nonlinear heat equation, and
# DImplicitEuler of OrdinaryDiffEqBDF on an index 1 DAE with an analytic
# Jacobian). Run from the repository root:
#   julia --project=tools/julia tools/irkc-probe.jl > irkc.tsv
#   REF=irkc.tsv OUT=tests/testthat/test-stabilized-irk.R Rscript tools/irkc-gen-test.R
using OrdinaryDiffEqStabilizedIRK, OrdinaryDiffEqBDF, SciMLBase
const N = 8
f1(u, p, t) = [p[1] * ((i == 1 ? zero(u[1]) : u[i-1]) - 2u[i] + (i == N ? zero(u[1]) : u[i+1])) * (N + 1)^2 - 0.5 * u[i] * u[i] * u[i] for i in 1:N]
f2(u, p, t) = [p[2] * u[i] * (1 - u[i]) + 0.1 * t for i in 1:N]
function j1(u, p, t)
    J = zeros(N, N)
    for i in 1:N
        J[i, i] = -2 * p[1] * (N + 1)^2 - 1.5 * u[i]^2
        i > 1 && (J[i, i-1] = p[1] * (N + 1)^2)
        i < N && (J[i, i+1] = p[1] * (N + 1)^2)
    end
    J
end
u0 = [i * (N + 1 - i) / N^2 for i in 1:N]
cb = DiscreteCallback((u, t, integ) -> t == 0.1, integ -> (integ.u = integ.u * 0.5))
for (pp, rt) in (([1.0, 5.0], 1e-4), ([0.2, 40.0], 1e-3))
    prob = SplitODEProblem(SplitFunction(f1, f2; jac = j1), u0, (0.0, 0.3), pp)
    for sc in ("plain", "tstop", "cb", "eig", "fixed")
        s = sc == "plain" ? solve(prob, IRKC(); reltol = rt, abstol = 1e-6) :
            sc == "tstop" ? solve(prob, IRKC(); reltol = rt, abstol = 1e-6, tstops = [0.1]) :
            sc == "cb" ? solve(prob, IRKC(); reltol = rt, abstol = 1e-6, tstops = [0.1], callback = cb) :
            sc == "eig" ? solve(prob, IRKC(eigen_est = integ -> (integ.eigen_est = 30.0)); reltol = rt, abstol = 1e-6) :
            solve(prob, IRKC(); dt = 1 / 32, adaptive = false)
        st = s.stats
        println("IRKC\t", pp[2], "\t", rt, "\t", sc, "\t", length(s.t), "\t", st.nf, "\t", st.nf2, "\t", st.njacs, "\t", st.nw, "\t", st.nsolve, "\t", st.nreject, "\t", join(repr.(s.u[end]), ","))
    end
end
F(du, u, p, t) = [du[1] - u[2], u[2] * u[2] * u[2] + u[2] + u[1] * u[1] * u[1] + u[1]]
J(du, u, p, gamma, t) = [gamma -1.0; 3u[1]^2+1 3u[2]^2+1]
dprob = DAEProblem(DAEFunction(F; jac = J), [-1.0, -1.0], [1.0, -1.0], (0.0, 2.0); differential_vars = [true, false])
for (rt, at) in ((1e-4, 1e-6), (1e-6, 1e-8))
    s = solve(dprob, DImplicitEuler(); reltol = rt, abstol = at, initializealg = SciMLBase.NoInit())
    st = s.stats
    println("DImplicitEuler\t", rt, "\t", at, "\tplain\t", length(s.t), "\t", st.nf, "\t0\t", st.njacs, "\t", st.nw, "\t", st.nsolve, "\t", st.nreject, "\t", join(repr.(s.u[end]), ","))
end
