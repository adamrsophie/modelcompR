# Reference values for tests/testthat/test-newmark.R: OrdinaryDiffEqNewmark
# on second order test problems, printed as R list entries (out of place).
#
#   julia --project=tools/julia tools/newmark-probe.jl
#
# Problems (state ArrayPartition((v, u)), R state c(v, u)):
#   L1  u'' = -4 u,               v0 = (0.5, -0.3), u0 = (1, 0.2), t in (0, 1)
#   P1  u'' = -u - u^3,           v0 = (0.5, -0.3), u0 = (1, 0.2), t in (0, 1)
#   P5  u'' = -u + t,             v0 = (0.5, -0.3), u0 = (1, 0.2), t in (0, 1)
using OrdinaryDiffEq, OrdinaryDiffEqNewmark
const C = OrdinaryDiffEq.OrdinaryDiffEqCore

const PROBS = Dict(
    "L1" => DynamicalODEProblem((v, u, p, t) -> -4 .* u, (v, u, p, t) -> v, [0.5, -0.3], [1.0, 0.2], (0.0, 1.0)),
    "P1" => DynamicalODEProblem((v, u, p, t) -> -u .- u .* u .* u, (v, u, p, t) -> v, [0.5, -0.3], [1.0, 0.2], (0.0, 1.0)),
    "P5" => DynamicalODEProblem((v, u, p, t) -> -u .+ t, (v, u, p, t) -> v, [0.5, -0.3], [1.0, 0.2], (0.0, 1.0)),
)

rv(x) = "c(" * join(repr.(collect(x)), ", ") * ")"
function emit(name, case, sol; dt = nothing, tol = nothing, tq = nothing)
    s = "  list(alg = \"$name\", case = \"$case\""
    dt === nothing || (s *= ", dt = $(repr(dt))")
    tol === nothing || (s *= ", tol = $(repr(tol))")
    st = sol.stats
    s *= ", nf = $(st.nf), nf2 = $(st.nf2), naccept = $(st.naccept), nreject = $(st.nreject), nt = $(length(sol.t)),\n" *
        "       t2 = $(repr(sol.t[2])), u = $(rv(sol.u[end]))"
    tq === nothing || (s *= ",\n       tq = $(repr(tq)), uq = $(rv(sol(tq)))")
    println(s * "),")
end

const ALGS = [
    ("NewmarkBeta()", NewmarkBeta()),
    ("NewmarkBeta(0.3, 0.6)", NewmarkBeta(0.3, 0.6)),
    ("GeneralizedAlpha(rho_inf = 0.5)", GeneralizedAlpha(rho_inf = 0.5)),
    ("GeneralizedAlpha(alpha_hht = -0.1)", GeneralizedAlpha(alpha_hht = -0.1)),
    ("GeneralizedAlpha(-0.2, 0.2, 0.36, 0.9)", GeneralizedAlpha(-0.2, 0.2, 0.36, 0.9)),
]
for (nm, alg) in ALGS
    println("# $nm order $(C.alg_order(alg)) fields ", alg isa GeneralizedAlpha ?
        rv((alg.αm, alg.αf, alg.β, alg.γ)) : rv((alg.β, alg.γ)))
end
println("# fixed dt")
for (nm, alg) in ALGS, case in ("L1", "P1", "P5")
    emit(nm, case, solve(PROBS[case], alg, dt = 0.1, adaptive = false); dt = 0.1, tq = 0.55)
end
println("# adaptive")
for (nm, alg) in ALGS
    emit(nm, "L1", solve(PROBS["L1"], alg); tq = 0.55)
    emit(nm, "P1", solve(PROBS["P1"], alg); tq = 0.55)
    emit(nm, "P1", solve(PROBS["P1"], alg, reltol = 1e-6, abstol = 1e-6); tol = 1e-6)
end
# the linear acceleration method: blind error estimate
try
    solve(PROBS["P1"], NewmarkBeta(1 / 6, 1 / 2))
catch e
    println("# NewmarkBeta(1/6, 1/2) adaptive: ", typeof(e))
end
emit("NewmarkBeta(1/6, 1/2)", "P1", solve(PROBS["P1"], NewmarkBeta(1 / 6, 1 / 2), dt = 0.1, adaptive = false); dt = 0.1)
# saveat during a fixed step solve against the dense solution
let s = solve(PROBS["P1"], NewmarkBeta(), dt = 0.1, adaptive = false, saveat = [0.15, 0.55]),
        d = solve(PROBS["P1"], NewmarkBeta(), dt = 0.1, adaptive = false)
    println("# saveat t = $(rv(s.t))")
    println("#   saveat u(0.15) = $(rv(s.u[1])), u(0.55) = $(rv(s.u[2]))")
    println("#   dense  u(0.15) = $(rv(d(0.15))), u(0.55) = $(rv(d(0.55)))")
end
# start at rest: the first steps leave the acceleration unchanged to rounding,
# so the Zienkiewicz-Xie estimate is 0 (taken as 1) until the Newton
# correction moves it
let p0 = SecondOrderODEProblem((v, u, p, t) -> -u, [0.0], [1.0], (0.0, 1.0))
    for (nm, alg) in (("NewmarkBeta()", NewmarkBeta()), ("GeneralizedAlpha(rho_inf = 0.8)", GeneralizedAlpha(rho_inf = 0.8)))
        s = solve(p0, alg, reltol = 1e-6)
        println("# at rest $nm: $(s.retcode) nf = $(s.stats.nf), naccept = $(s.stats.naccept), ",
                "nreject = $(s.stats.nreject), t2 = $(repr(s.t[2])), t3 = $(repr(s.t[3]))")
    end
end
