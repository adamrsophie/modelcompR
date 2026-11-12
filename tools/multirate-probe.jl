# Reference values for tests/testthat/test-multirate.R: the explicit methods
# of OrdinaryDiffEqMultirate (MREEF, MRAB, MRIGARKERK22a/22b/33a/45a, MIS)
# and the implicit ones of OrdinaryDiffEqMultirateImplicit (MREIL,
# MRIGARKIRK21a, MRIGARKESDIRK34a, MRIGARKESDIRK46a) on out of place split
# problems du/dt = f1(u) + f2(u) (f1 fast, f2 slow), with forward difference
# Jacobians. OrdinaryDiffEqMultirateImplicit (1.0.1) is not registered: add it
# to a copy of the tools/julia environment with
#   Pkg.develop(path = "<OrdinaryDiffEq.jl>/lib/OrdinaryDiffEqMultirateImplicit")
# (with the [sources] entries of its Project.toml removed) and run
#   julia --project=<that environment> tools/multirate-probe.jl
using OrdinaryDiffEqMultirate, OrdinaryDiffEqMultirateImplicit, SciMLBase
using OrdinaryDiffEqDifferentiation: AutoFiniteDiff
using LinearAlgebra

# a fast part with a mildly nonlinear coupling and time dependence, a slow
# nonlinear part
f1(u, p, t) = [-p[1] * u[1] + u[2] * u[2] + 0.5 * t, -p[1] * (u[2] - u[1]) - 0.5 * t * u[2]]
f2(u, p, t) = [u[2] - 0.1 * u[1] * u[1] * u[1] + t, -u[1] + 0.2 * u[2] * t]
# a stiff slow part for the implicit MRI-GARK methods
g1(u, p, t) = [-2 * u[1] + u[2] + 0.25 * t, u[1] - 3 * u[2]]
g2(u, p, t) = [-p[1] * (u[1] - u[2] * u[2]) + t, -p[1] * u[2] + 0.5 * u[1] * u[1]]
u0 = [1.0, 0.5]
prob(p) = SplitODEProblem(f1, f2, u0, (0.0, 1.0), p)
probs(p) = SplitODEProblem(g1, g2, u0, (0.0, 1.0), p)
fd = AutoFiniteDiff()

function show(name, case, s)
    st = s.stats
    println(name, "\t", case, "\t", length(s.t), "\t", st.nf, "\t", st.nf2, "\t", st.njacs, "\t", st.nw, "\t",
        st.nsolve, "\t", st.naccept, "\t", st.nreject, "\t", join(repr.(s.u[end]), ","))
end

# tableau values the R port hard codes
tab = OrdinaryDiffEqMultirate.MIS2Tableau(Float64)
println("MIS2\td\t", join(repr.(tab.d), ","))
println("MIS2\tc\t", join(repr.(tab.c), ","))
println("MIS2\tctilde\t", join(repr.(tab.ctilde), ","))
println("MIS2\tbhat\t", join(repr.(tab.bhat), ","))
t34 = OrdinaryDiffEqMultirateImplicit.MRIGARKESDIRK34aTableau(Float64)
println("ESDIRK34a\tW0\t", join(repr.(vec(t34.W0)), ","))

for p in ([20.0], [60.0])
    P = prob(p)
    for (nm, alg) in (("MREEF", MREEF()), ("MREEF_r3", MREEF(m = 3, order = 3, seq = :romberg)),
            ("MREEF_o2", MREEF(m = 6, order = 2)),
            ("MRAB", MRAB()), ("MRAB_k1", MRAB(k = 1, m = 5)), ("MRAB_k3", MRAB(k = 3, m = 6)),
            ("MRAB_k4", MRAB(k = 4, m = 8)), ("MRAB_k5", MRAB(k = 5, m = 8)),
            ("MRIGARKERK22a", MRIGARKERK22a(m = 4)), ("MRIGARKERK22b", MRIGARKERK22b(m = 4)),
            ("MRIGARKERK33a", MRIGARKERK33a(m = 3)), ("MRIGARKERK45a", MRIGARKERK45a(m = 2)),
            ("MIS", MIS(m = 6)), ("MREIL", MREIL(autodiff = fd)),
            ("MREIL_r3", MREIL(m = 2, order = 3, seq = :romberg, autodiff = fd)))
        # MRAB with k = 1 has no error estimate (EEst = 1): the step never
        # grows, so the adaptive solve is cut at 50 steps
        kw = nm == "MRAB_k1" ? (maxiters = 50,) : NamedTuple()
        show(nm, "adaptive_" * string(p[1]), solve(P, alg; reltol = 1e-6, abstol = 1e-8, kw...))
        show(nm, "fixed_" * string(p[1]), solve(P, alg; dt = 1 / 16, adaptive = false))
    end
end
for p in ([30.0], [400.0])
    P = probs(p)
    for (nm, alg) in (("MRIGARKIRK21a", MRIGARKIRK21a(m = 3, autodiff = fd)),
            ("MRIGARKESDIRK34a", MRIGARKESDIRK34a(m = 2, autodiff = fd)),
            ("MRIGARKESDIRK46a", MRIGARKESDIRK46a(m = 2, autodiff = fd)))
        show(nm, "adaptive_" * string(p[1]), solve(P, alg; reltol = 1e-6, abstol = 1e-8))
        show(nm, "fixed_" * string(p[1]), solve(P, alg; dt = 1 / 16, adaptive = false))
    end
end

# an event that halves the state at t = 1/2, and the analytic Jacobian of f1
# for MREIL
cb = DiscreteCallback((u, t, integ) -> t == 0.5, integ -> (integ.u = integ.u * 0.5))
j1(u, p, t) = [-p[1] 2*u[2]; p[1] -p[1]-0.5*t]
P = prob([20.0])
for (nm, alg) in (("MREEF", MREEF()), ("MRAB_k3", MRAB(k = 3, m = 6)), ("MRIGARKERK33a", MRIGARKERK33a(m = 3)),
        ("MIS", MIS(m = 6)), ("MREIL", MREIL(autodiff = fd)))
    show(nm, "cbfixed_20.0", solve(P, alg; dt = 1 / 16, adaptive = false, callback = cb, tstops = [0.5]))
    show(nm, "cbadaptive_20.0", solve(P, alg; reltol = 1e-6, abstol = 1e-8, callback = cb, tstops = [0.5]))
end
PS = probs([30.0])
alg = MRIGARKESDIRK46a(m = 2, autodiff = fd)
show("MRIGARKESDIRK46a", "cbfixed_30.0", solve(PS, alg; dt = 1 / 16, adaptive = false, callback = cb, tstops = [0.5]))
show("MRIGARKESDIRK46a", "cbadaptive_30.0", solve(PS, alg; reltol = 1e-6, abstol = 1e-8, callback = cb, tstops = [0.5]))
PJ = SplitODEProblem(SplitFunction(f1, f2; jac = j1), u0, (0.0, 1.0), [20.0])
show("MREIL", "jacfixed_20.0", solve(PJ, MREIL(autodiff = fd); dt = 1 / 16, adaptive = false))
show("MREIL", "jacadaptive_20.0", solve(PJ, MREIL(autodiff = fd); reltol = 1e-6, abstol = 1e-8))

# loose tolerances: short adaptive solves of the methods whose solves above
# take thousands of steps
for (nm, alg, P, pv) in (("MREEF_o2", MREEF(m = 6, order = 2), prob, 20.0), ("MRAB", MRAB(), prob, 20.0),
        ("MRIGARKERK22a", MRIGARKERK22a(m = 4), prob, 20.0), ("MRIGARKERK22b", MRIGARKERK22b(m = 4), prob, 20.0),
        ("MIS", MIS(m = 6), prob, 20.0), ("MRIGARKIRK21a", MRIGARKIRK21a(m = 3, autodiff = fd), probs, 30.0),
        ("MRIGARKESDIRK34a", MRIGARKESDIRK34a(m = 2, autodiff = fd), probs, 30.0))
    show(nm, "loose_" * string(pv), solve(P([pv]), alg; reltol = 1e-3, abstol = 1e-5))
end
# MREIL with a diagonal mass matrix: W = J - M / h_fast
PM = SplitODEProblem(SplitFunction(f1, f2; mass_matrix = Diagonal([1.0, 2.0])), u0, (0.0, 1.0), [20.0])
show("MREIL", "massfixed_20.0", solve(PM, MREIL(autodiff = fd); dt = 1 / 16, adaptive = false))
show("MREIL", "massadaptive_20.0", solve(PM, MREIL(autodiff = fd); reltol = 1e-6, abstol = 1e-8))
