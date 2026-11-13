# Reference values for tests/testthat/test-bdf.R: the out of place methods of
# OrdinaryDiffEqBDF (ABDF2, QNDF1, QNDF2, QNDF, QBDF1, QBDF2, QBDF, MEBDF2,
# NordsieckBDF, the SBDF family, DABDF2 and DNordsieckBDF) on a nonlinear heat
# equation with an analytic Jacobian, its split form and an index 1 DAE.
# Each line: method, case, number of saved points, nf, nf2, njacs, nw, nsolve,
# naccept, nreject and the final state (repr). Lines starting with "dense"
# give the dense output at t = 0.05, 0.123 and 0.29 (DAE: 0.3, 0.77, 1.51;
# Robertson: 1e-3, 0.5, 4e4); "ctrl" lines the controller settings; "initdt"
# lines the automatic first step, which the "dt0" cases pass explicitly.
# Run from the repository root:
#   julia --project=tools/julia --startup-file=no tools/bdf-probe.jl > bdf.tsv
using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, SciMLBase
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
f1(u, p, t) = [p[1] * ((i == 1 ? zero(u[1]) : u[i-1]) - 2u[i] + (i == N ? zero(u[1]) : u[i+1])) * (N + 1)^2 - 0.5 * u[i] * u[i] * u[i] for i in 1:N]
f2(u, p, t) = [p[2] * u[i] * (1 - u[i]) + 0.1 * t for i in 1:N]
u0 = [i * (N + 1 - i) / N^2 for i in 1:N]
prob = ODEProblem(ODEFunction(heat; jac = hjac), u0, (0.0, 0.3), [1.0])
sprob = SplitODEProblem(SplitFunction(f1, f2; jac = hjac), u0, (0.0, 0.3), [1.0, 5.0])
cb = DiscreteCallback((u, t, integ) -> t == 0.1, integ -> (integ.u = integ.u * 0.5))

pr(nm, sc, s) = println(nm, "\t", sc, "\t", length(s.t), "\t", s.stats.nf, "\t", s.stats.nf2, "\t", s.stats.njacs, "\t",
    s.stats.nw, "\t", s.stats.nsolve, "\t", s.stats.naccept, "\t", s.stats.nreject, "\t", join(repr.(s.u[end]), ","))

for (nm, a) in [("ABDF2", ABDF2()), ("QNDF1", QNDF1()), ("QNDF2", QNDF2()), ("QNDF", QNDF()),
                ("QBDF1", QBDF1()), ("QBDF2", QBDF2()), ("QBDF", QBDF()), ("NordsieckBDF", NordsieckBDF()),
                ("MEBDF2", MEBDF2()), ("SBDF2", SBDF2()), ("DABDF2", DABDF2()), ("DNordsieckBDF", DNordsieckBDF())]
    p = nm in ("SBDF2",) ? sprob : prob
    if nm in ("DABDF2", "DNordsieckBDF")
        continue
    end
    it = init(p, a; dt = 0.01, adaptive = nm in ("MEBDF2", "SBDF2") ? false : true)
    println("ctrl\t", nm, "\t", OrdinaryDiffEqCore.get_qmin(it), "\t", OrdinaryDiffEqCore.get_qmax(it), "\t",
        OrdinaryDiffEqCore.get_gamma(it), "\t", OrdinaryDiffEqCore.get_qsteady_min(it), "\t",
        OrdinaryDiffEqCore.get_qsteady_max(it), "\t", OrdinaryDiffEqCore.get_failfactor(it), "\t",
        typeof(it.controller_cache).name.name, "\t", OrdinaryDiffEqCore.isadaptive(a))
end

tols = Dict("lo" => (1e-4, 1e-6), "hi" => (1e-7, 1e-9))
adaptive_algs = [("ABDF2", () -> ABDF2()), ("QNDF1", () -> QNDF1()), ("QNDF2", () -> QNDF2()), ("QNDF", () -> QNDF()),
                 ("QBDF1", () -> QBDF1()), ("QBDF2", () -> QBDF2()), ("QBDF", () -> QBDF()),
                 ("NordsieckBDF", () -> NordsieckBDF()), ("NordsieckBDFstald", () -> NordsieckBDF(stald = true)),
                 ("QNDF3", () -> QNDF(max_order = Val(3))), ("NordsieckBDF2", () -> NordsieckBDF(max_order = Val(2)))]
for (nm, mk) in adaptive_algs, sc in ("plain", "tstop", "cb", "hi")
    rt, at = sc == "hi" ? tols["hi"] : tols["lo"]
    s = sc == "tstop" ? solve(prob, mk(); reltol = rt, abstol = at, tstops = [0.1]) :
        sc == "cb" ? solve(prob, mk(); reltol = rt, abstol = at, tstops = [0.1], callback = cb) :
        solve(prob, mk(); reltol = rt, abstol = at)
    pr(nm, sc, s)
    if sc == "plain"
        println("dense\t", nm, "\t", join(repr.(vcat(s(0.05), s(0.123), s(0.29))), ","))
    end
end

# The automatic first step of the heat problem, and runs that pass it as dt
# (QNDF chooses with exact comparisons of dt, which a one ulp change of the
# first step can flip)
for (sc, rt, at) in (("lo", 1e-4, 1e-6), ("hi", 1e-7, 1e-9))
    println("initdt\t", sc, "\t", repr(init(prob, QNDF(); reltol = rt, abstol = at).dt))
end
dt_lo = init(prob, QNDF(); reltol = 1e-4, abstol = 1e-6).dt
dt_hi = init(prob, QNDF(); reltol = 1e-7, abstol = 1e-9).dt
for (nm, mk) in [("QNDF", () -> QNDF()), ("QBDF", () -> QBDF()), ("QNDF3", () -> QNDF(max_order = Val(3)))],
        sc in ("plain_dt0", "tstop_dt0", "hi_dt0")
    rt, at, d = sc == "hi_dt0" ? (1e-7, 1e-9, dt_hi) : (1e-4, 1e-6, dt_lo)
    s = sc == "tstop_dt0" ? solve(prob, mk(); reltol = rt, abstol = at, dt = d, tstops = [0.1]) :
        solve(prob, mk(); reltol = rt, abstol = at, dt = d)
    pr(nm, sc, s)
end

# Robertson's problem: QNDF meets a Newton failure on its first steps
rober(u, p, t) = [-0.04 * u[1] + 1e4 * u[2] * u[3], 0.04 * u[1] - 1e4 * u[2] * u[3] - 3e7 * u[2] * u[2], 3e7 * u[2] * u[2]]
rjac(u, p, t) = [-0.04 1e4*u[3] 1e4*u[2]; 0.04 -1e4*u[3]-6e7*u[2] -1e4*u[2]; 0.0 6e7*u[2] 0.0]
rprob = ODEProblem(ODEFunction(rober; jac = rjac), [1.0, 0.0, 0.0], (0.0, 1e5))
rsol = solve(rprob, QNDF(); reltol = 1e-3, abstol = 1e-6)
pr("QNDF", "rober", rsol)
println("dense\tQNDFrober\t", join(repr.(vcat(rsol(1e-3), rsol(0.5), rsol(4e4))), ","))
println("nonlinconvfail\tQNDFrober\t", rsol.stats.nnonlinconvfail)

# fixed step runs
for (nm, a) in [("MEBDF2", MEBDF2()), ("ABDF2", ABDF2()), ("QNDF", QNDF()), ("NordsieckBDF", NordsieckBDF()),
                ("QNDF2", QNDF2()), ("QNDF1", QNDF1())]
    s = solve(prob, a; dt = 1 / 64, adaptive = false)
    pr(nm, "fixed", s)
    nm == "MEBDF2" && println("dense\t", nm, "\t", join(repr.(vcat(s(0.05), s(0.123), s(0.29))), ","))
end
for (nm, a) in [("IMEXEuler", IMEXEuler()), ("IMEXEulerARK", IMEXEulerARK()), ("SBDF2", SBDF2()), ("SBDF3", SBDF3()),
                ("SBDF4", SBDF4()), ("SBDF3ark", SBDF(3; ark = true)), ("SBDF4ark", SBDF4(ark = true))]
    for sc in ("fixed", "tstop", "cb")
        s = sc == "fixed" ? solve(sprob, a; dt = 1 / 64) :
            sc == "tstop" ? solve(sprob, a; dt = 1 / 64, tstops = [0.1]) :
            solve(sprob, a; dt = 1 / 64, tstops = [0.1], callback = cb)
        pr(nm, sc, s)
    end
end

# DAE: u1' = u2, 0 = u2^3 + u2 + u1^3 + u1
F(du, u, p, t) = [du[1] - u[2], u[2] * u[2] * u[2] + u[2] + u[1] * u[1] * u[1] + u[1]]
J(du, u, p, gamma, t) = [gamma -1.0; 3u[1]^2+1 3u[2]^2+1]
dprob = DAEProblem(DAEFunction(F; jac = J), [-1.0, -1.0], [1.0, -1.0], (0.0, 2.0); differential_vars = [true, false])
for (nm, a) in [("DABDF2", DABDF2()), ("DNordsieckBDF", DNordsieckBDF())], (sc, rt, at) in (("lo", 1e-4, 1e-6), ("hi", 1e-6, 1e-8))
    s = solve(dprob, a; reltol = rt, abstol = at, initializealg = SciMLBase.NoInit())
    pr(nm, sc, s)
    sc == "lo" && println("dense\t", nm, "\t", join(repr.(vcat(s(0.3), s(0.77), s(1.51))), ","))
end
