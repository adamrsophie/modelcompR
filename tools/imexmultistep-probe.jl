# Reference values for tests/testthat/test-imex-multistep.R: CNAB2 and CNLF2
# of OrdinaryDiffEqIMEXMultistep on an out of place split problem
# du/dt = f1(u) + f2(u) (f1 implicit, f2 explicit), with forward difference
# or analytic Jacobians of f1, fixed steps, and an event that halves the
# state at t = 1/2 (the next step restarts with the startup formula).
#   julia --project=tools/julia tools/imexmultistep-probe.jl
using OrdinaryDiffEqIMEXMultistep, SciMLBase
using OrdinaryDiffEqDifferentiation: AutoFiniteDiff

f1(u, p, t) = [-p[1] * u[1] + u[2] * u[2] + 0.5 * t, -p[1] * (u[2] - u[1]) - 0.5 * t * u[2]]
f2(u, p, t) = [u[2] - 0.1 * u[1] * u[1] * u[1] + t, -u[1] + 0.2 * u[2] * t]
j1(u, p, t) = [-p[1] 2*u[2]; p[1] -p[1]-0.5*t]
u0 = [1.0, 0.5]
fd = AutoFiniteDiff()
cb = DiscreteCallback((u, t, integ) -> t == 0.5, integ -> (integ.u = integ.u * 0.5))

function show(name, case, s)
    st = s.stats
    println(name, "\t", case, "\t", length(s.t), "\t", st.nf, "\t", st.nf2, "\t", st.njacs, "\t", st.nw, "\t",
        st.nsolve, "\t", st.naccept, "\t", st.nreject, "\t", join(repr.(s.u[end]), ","))
end

for p in ([20.0], [200.0])
    P = SplitODEProblem(f1, f2, u0, (0.0, 1.0), p)
    PJ = SplitODEProblem(SplitFunction(f1, f2; jac = j1), u0, (0.0, 1.0), p)
    for (nm, alg) in (("CNAB2", CNAB2(autodiff = fd)), ("CNLF2", CNLF2(autodiff = fd)))
        show(nm, "fd16_" * string(p[1]), solve(P, alg; dt = 1 / 16))
        show(nm, "fd40_" * string(p[1]), solve(P, alg; dt = 1 / 40))
        show(nm, "jac16_" * string(p[1]), solve(PJ, alg; dt = 1 / 16))
        show(nm, "cb16_" * string(p[1]), solve(P, alg; dt = 1 / 16, callback = cb, tstops = [0.5]))
    end
end
