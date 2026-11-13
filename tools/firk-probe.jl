# Reference runs of OrdinaryDiffEqFIRK (RadauIIA3, RadauIIA5, RadauIIA9,
# AdaptiveRadau, GaussLegendre), OrdinaryDiffEqPDIRK (PDIRK44),
# OrdinaryDiffEqQPRK (QPRK98) and OrdinaryDiffEqPRK (KuttaPRK2p5) on out of
# place problems: the tab separated lines (method, case, saved points, nf,
# njacs, nw, nsolve, nreject, naccept, final state) are the numbers of
# tests/testthat/test-firk.R. Run from the repository root:
#   julia --project=tools/julia tools/firk-probe.jl [group...]
# with groups among prk, pdirk, radau, gauss, adaptive, extra, defects,
# dense (default: all).
using OrdinaryDiffEqFIRK, OrdinaryDiffEqPDIRK, OrdinaryDiffEqQPRK, OrdinaryDiffEqPRK
using OrdinaryDiffEqCore, SciMLBase

# Lotka-Volterra (nonstiff, polynomial)
lv(u, p, t) = [p[1] * u[1] - u[1] * u[2], -p[2] * u[2] + u[1] * u[2]]
lvprob = ODEProblem(lv, [1.0, 1.0], (0.0, 2.0), [1.5, 3.0])

# nonlinear heat equation with an analytic Jacobian (stiff)
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
heatprob = ODEProblem(ODEFunction(heat; jac = hjac), [i * (N + 1 - i) / N^2 for i in 1:N], (0.0, 0.3), [1.0])

# Van der Pol with mu = 10 and its Jacobian (stiff, polynomial)
vdp(u, p, t) = [u[2], p[1] * ((1 - u[1] * u[1]) * u[2]) - u[1]]
vdpjac(u, p, t) = [0.0 1.0; -2 * p[1] * u[1] * u[2]-1 p[1] * (1 - u[1] * u[1])]
vdpprob = ODEProblem(ODEFunction(vdp; jac = vdpjac), [2.0, 0.0], (0.0, 6.0), [10.0])

cbheat = DiscreteCallback((u, t, integ) -> t == 0.1, integ -> (integ.u = integ.u * 0.5))
cblv = DiscreteCallback((u, t, integ) -> t == 1.0, integ -> (integ.u = integ.u * 0.5))

function report(nm, sc, s)
    st = s.stats
    println(nm, "\t", sc, "\t", length(s.t), "\t", st.nf, "\t", st.njacs, "\t", st.nw, "\t", st.nsolve, "\t",
        st.nreject, "\t", st.naccept, "\t", join(repr.(s.u[end]), ","))
end

want(g) = isempty(ARGS) || "all" in ARGS || g in ARGS
if want("prk")
    report("KuttaPRK2p5", "fixed", solve(lvprob, KuttaPRK2p5(), dt = 0.05))
    report("KuttaPRK2p5", "cb", solve(lvprob, KuttaPRK2p5(), dt = 0.05, tstops = [1.0], callback = cblv))
    report("QPRK98", "plain", solve(lvprob, QPRK98(), reltol = 1e-8, abstol = 1e-10))
    report("QPRK98", "tstop", solve(lvprob, QPRK98(), reltol = 1e-8, abstol = 1e-10, tstops = [1.0]))
    report("QPRK98", "cb", solve(lvprob, QPRK98(), reltol = 1e-8, abstol = 1e-10, tstops = [1.0], callback = cblv))
    report("QPRK98", "fixed", solve(lvprob, QPRK98(), dt = 0.25, adaptive = false))
end
if want("pdirk")
    for th in (true, false)
        nm = th ? "PDIRK44" : "PDIRK44_serial"
        report(nm, "heat", solve(heatprob, PDIRK44(threading = th), dt = 0.01))
        report(nm, "heatcb", solve(heatprob, PDIRK44(threading = th), dt = 0.01, tstops = [0.1], callback = cbheat))
        report(nm, "vdp", solve(vdpprob, PDIRK44(threading = th), dt = 0.01))
    end
end
const FIRK = [("RadauIIA3", RadauIIA3), ("RadauIIA5", RadauIIA5), ("RadauIIA9", RadauIIA9)]
if want("radau")
    for (nm, A) in FIRK
        report(nm, "plain", solve(heatprob, A(), reltol = 1e-4, abstol = 1e-6))
        report(nm, "tstop", solve(heatprob, A(), reltol = 1e-4, abstol = 1e-6, tstops = [0.1]))
        report(nm, "cb", solve(heatprob, A(), reltol = 1e-4, abstol = 1e-6, tstops = [0.1], callback = cbheat))
        report(nm, "vdp", solve(vdpprob, A(), reltol = 1e-6, abstol = 1e-8))
        report(nm, "fixed", solve(vdpprob, A(), dt = 0.05, adaptive = false))
    end
end
if want("gauss")
    for s in (2, 3, 5)
        nm = "GaussLegendre$s"
        report(nm, "plain", solve(heatprob, GaussLegendre(num_stages = s), reltol = 1e-4, abstol = 1e-6))
        report(nm, "cb", solve(heatprob, GaussLegendre(num_stages = s), reltol = 1e-4, abstol = 1e-6, tstops = [0.1], callback = cbheat))
        report(nm, "vdp", solve(vdpprob, GaussLegendre(num_stages = s), reltol = 1e-6, abstol = 1e-8))
    end
    report("GaussLegendre1", "fixed", solve(vdpprob, GaussLegendre(num_stages = 1), dt = 0.05, adaptive = false))
    report("GaussLegendre3", "fixed", solve(vdpprob, GaussLegendre(num_stages = 3), dt = 0.05, adaptive = false))
end
if want("adaptive")
    report("AdaptiveRadau", "plain", solve(heatprob, AdaptiveRadau(), reltol = 1e-4, abstol = 1e-6))
    report("AdaptiveRadau", "cb", solve(heatprob, AdaptiveRadau(), reltol = 1e-4, abstol = 1e-6, tstops = [0.1], callback = cbheat))
    report("AdaptiveRadau", "vdp", solve(vdpprob, AdaptiveRadau(), reltol = 1e-6, abstol = 1e-8))
    report("AdaptiveRadau", "vdp_tight", solve(vdpprob, AdaptiveRadau(), reltol = 1e-10, abstol = 1e-12))
    report("AdaptiveRadau_9_13", "vdp", solve(vdpprob, AdaptiveRadau(min_order = 9, max_order = 13), reltol = 1e-8, abstol = 1e-10))
    report("AdaptiveRadau_5_9", "heat", solve(heatprob, AdaptiveRadau(min_order = 5, max_order = 9), reltol = 1e-8, abstol = 1e-10))
end
# Robertson's problem in mass matrix form (index one DAE) with its Jacobian
rober(u, p, t) = [-p[1] * u[1] + p[2] * u[2] * u[3], p[1] * u[1] - p[2] * u[2] * u[3] - p[3] * (u[2] * u[2]),
    u[1] + u[2] + u[3] - 1]
roberjac(u, p, t) = [-p[1] p[2]*u[3] p[2]*u[2]; p[1] -p[2]*u[3]-2*p[3]*u[2] -p[2]*u[2]; 1.0 1.0 1.0]
roberprob = ODEProblem(ODEFunction(rober; mass_matrix = [1.0 0 0; 0 1.0 0; 0 0 0], jac = roberjac),
    [1.0, 0.0, 0.0], (0.0, 10.0), [0.04, 1.0e4, 3.0e7])
if want("extra")
    # tableaus Julia computes at run time (9 Radau stages, 6 Gauss-Legendre stages)
    report("AdaptiveRadau_13_17", "vdp", solve(vdpprob, AdaptiveRadau(min_order = 13, max_order = 17), reltol = 1e-10, abstol = 1e-12))
    report("GaussLegendre6", "vdp", solve(vdpprob, GaussLegendre(num_stages = 6), reltol = 1e-8, abstol = 1e-10))
    report("RadauIIA5_constant", "heat", solve(heatprob, RadauIIA5(extrapolant = :constant), reltol = 1e-4, abstol = 1e-6))
    report("RadauIIA5_kappa", "vdp", solve(vdpprob, RadauIIA5(κ = 0.1, smooth_est = false), reltol = 1e-6, abstol = 1e-8))
    report("RadauIIA5", "rober", solve(roberprob, RadauIIA5(), reltol = 1e-6, abstol = 1e-8))
    report("RadauIIA9", "rober", solve(roberprob, RadauIIA9(), reltol = 1e-6, abstol = 1e-8))
    report("AdaptiveRadau", "rober", solve(roberprob, AdaptiveRadau(), reltol = 1e-6, abstol = 1e-8))
    # with the default initial step the first error estimate of RadauIIA9 is
    # rounding (about 1e-18), so it is given here
    report("RadauIIA9", "rober_dt", solve(roberprob, RadauIIA9(), dt = 5e-4, reltol = 1e-6, abstol = 1e-8))
end
# Defects of the out of place steps that the R port fixes (see
# tests/testthat/test-firk.R): RadauIIA3 has no Newton convergence check,
# AdaptiveRadau computes mass_matrix * w across the stages. Also the six
# stage Gauss-Legendre tableau Julia builds at run time.
cub = ODEProblem(ODEFunction((u, p, t) -> -p[1] .* (u .* u .* u); jac = (u, p, t) -> reshape([-3 * p[1] * (u[1] * u[1])], 1, 1)),
    [1.0], (0.0, 2.0), [1000.0])
if want("defects")
    report("RadauIIA3_loose", "rober", solve(roberprob, RadauIIA3(), reltol = 1e-3, abstol = 1e-6))
    report("RadauIIA3_cubic", "fixed", solve(cub, RadauIIA3(), dt = 0.25, adaptive = false))
    s = solve(roberprob, AdaptiveRadau(), reltol = 1e-6, abstol = 1e-8)
    println("AdaptiveRadau\trober_retcode\t", s.retcode, "\t", repr(s.t[end]))
    tab = OrdinaryDiffEqFIRK.GaussLegendreTableau(Float64, Float64, 6)
    for nm in (:A, :b, :c, :e)
        println("GaussLegendreTableau6\t", nm, "\t", join(repr.(vec(getfield(tab, nm))), ","))
    end
end
# Dense output between the steps: method, "dense", case, time, state
function dense(nm, sc, s, ts)
    for t in ts
        println(nm, "\tdense\t", sc, "\t", repr(t), "\t", join(repr.(s(t)), ","))
    end
end
if want("dense")
    th = (0.0137, 0.1234, 0.2718)
    tv = (0.37, 2.5, 5.55)
    for (nm, A) in FIRK
        dense(nm, "heat", solve(heatprob, A(), reltol = 1e-4, abstol = 1e-6), th)
        dense(nm, "vdp", solve(vdpprob, A(), reltol = 1e-6, abstol = 1e-8), tv)
    end
    dense("RadauIIA5_constant", "heat", solve(heatprob, RadauIIA5(extrapolant = :constant), reltol = 1e-4, abstol = 1e-6), th)
    dense("GaussLegendre3", "heat", solve(heatprob, GaussLegendre(num_stages = 3), reltol = 1e-4, abstol = 1e-6), th)
    dense("GaussLegendre3", "vdp", solve(vdpprob, GaussLegendre(num_stages = 3), reltol = 1e-6, abstol = 1e-8), tv)
    dense("AdaptiveRadau", "heat", solve(heatprob, AdaptiveRadau(), reltol = 1e-4, abstol = 1e-6), th)
    dense("AdaptiveRadau", "vdp", solve(vdpprob, AdaptiveRadau(), reltol = 1e-6, abstol = 1e-8), tv)
    dense("AdaptiveRadau_9_13", "vdp", solve(vdpprob, AdaptiveRadau(min_order = 9, max_order = 13), reltol = 1e-8, abstol = 1e-10), tv)
    dense("PDIRK44", "heat", solve(heatprob, PDIRK44(), dt = 0.01), th)
    dense("QPRK98", "lv", solve(lvprob, QPRK98(), reltol = 1e-8, abstol = 1e-10), (0.37, 1.23))
    dense("KuttaPRK2p5", "lv", solve(lvprob, KuttaPRK2p5(), dt = 0.05), (0.37, 1.23))
end
