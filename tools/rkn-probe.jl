# Reference values for tests/testthat/test-rkn.R and test-symplectic.R:
# OrdinaryDiffEqRKN and OrdinaryDiffEqSymplecticRK on second order test
# problems, printed as R list entries (out of place functions).
#
#   julia --project=tools/julia tools/rkn-probe.jl
#
# Problems (state ArrayPartition((v, u)), R state c(v, u)):
#   P1  u'' = -u - u^3,           v0 = (0.5, -0.3), u0 = (1, 0.2), t in (0, 1)
#   P2  u'' = -u / |u|^3 (Kepler), v0 = (0, 1.2),    u0 = (1, 0),   t in (0, 3)
#   P3  u'' = -u - u^3 - 0.1 u',   v0 = (0.5, -0.3), u0 = (1, 0.2), t in (0, 1)
#   P4  P1 with u' = v / 2
#   P5  u'' = -u + t,             v0 = (0.5, -0.3), u0 = (1, 0.2), t in (0, 1)
using OrdinaryDiffEq, OrdinaryDiffEqRKN, OrdinaryDiffEqSymplecticRK

kep(v, u, p, t) = (r = sqrt(u[1] * u[1] + u[2] * u[2]); -u ./ (r * r * r))
const PROBS = Dict(
    "P1" => DynamicalODEProblem((v, u, p, t) -> -u .- u .* u .* u, (v, u, p, t) -> v, [0.5, -0.3], [1.0, 0.2], (0.0, 1.0)),
    "P2" => DynamicalODEProblem(kep, (v, u, p, t) -> v, [0.0, 1.2], [1.0, 0.0], (0.0, 3.0)),
    "P3" => DynamicalODEProblem((v, u, p, t) -> -u .- u .* u .* u .- 0.1 .* v, (v, u, p, t) -> v, [0.5, -0.3], [1.0, 0.2], (0.0, 1.0)),
    "P4" => DynamicalODEProblem((v, u, p, t) -> -u .- u .* u .* u, (v, u, p, t) -> 0.5 .* v, [0.5, -0.3], [1.0, 0.2], (0.0, 1.0)),
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
fixed(name, alg, case; dt = 0.1, tq = 0.55) =
    emit(name, case, solve(PROBS[case], alg, dt = dt, adaptive = false); dt = dt, tq = tq)
adapt(name, alg, case; tol = 1e-8, tq = nothing) =
    emit(name, case, solve(PROBS[case], alg, reltol = tol, abstol = tol); tol = tol, tq = tq)

const SYMP = (:SymplecticEuler, :VelocityVerlet, :VerletLeapfrog, :LeapfrogDriftKickDrift, :PseudoVerletLeapfrog,
    :McAte2, :Ruth3, :McAte3, :CandyRoz4, :McAte4, :CalvoSanz4, :McAte42, :McAte5, :Yoshida6, :KahanLi6,
    :McAte8, :KahanLi8, :SofSpa10)
println("# symplectic, fixed dt")
for nm in SYMP, case in ("P1", "P5")
    fixed(String(nm), getfield(OrdinaryDiffEqSymplecticRK, nm)(), case)
end
for nm in (:SymplecticEuler, :LeapfrogDriftKickDrift, :McAte3, :KahanLi6, :SofSpa10), case in ("P3", "P4")
    fixed(String(nm), getfield(OrdinaryDiffEqSymplecticRK, nm)(), case)
end
# VerletLeapfrog: linear interpolation of the saved steps by default
let sol = solve(PROBS["P1"], VerletLeapfrog(), dt = 0.1)
    println("# VerletLeapfrog dense = $(sol.dense), sol(0.55) = $(rv(sol(0.55)))")
end
let sol = solve(PROBS["P1"], McAte5(), dt = 0.1, saveat = [0.15, 0.55])
    println("# McAte5 saveat: t = $(rv(sol.t)), u = $(rv(sol.u[1])), $(rv(sol.u[2]))")
end

const VI = (:Nystrom4VelocityIndependent, :Nystrom5VelocityIndependent, :DPRKN4, :DPRKN5, :DPRKN6, :DPRKN6FM,
    :DPRKN8, :DPRKN12, :ERKN4, :ERKN5, :ERKN7)
const VD = (:Nystrom4, :FineRKN4, :FineRKN5, :RKN4)
const ADAPT = (:DPRKN4, :DPRKN5, :DPRKN6, :DPRKN6FM, :DPRKN8, :DPRKN12, :ERKN4, :ERKN5, :ERKN7, :FineRKN4, :FineRKN5)
println("# RKN, fixed dt")
for nm in (VI..., VD...), case in ("P1", "P5")
    fixed(String(nm), getfield(OrdinaryDiffEqRKN, nm)(), case)
end
for nm in VD
    fixed(String(nm), getfield(OrdinaryDiffEqRKN, nm)(), "P3")
end
println("# RKN, adaptive")
for nm in ADAPT
    alg = getfield(OrdinaryDiffEqRKN, nm)()
    adapt(String(nm), alg, "P1"; tq = 0.55)
    adapt(String(nm), alg, "P2"; tq = 1.7)
    adapt(String(nm), alg, "P1"; tol = 1e-4)
end
for nm in (:FineRKN4, :FineRKN5)
    adapt(String(nm), getfield(OrdinaryDiffEqRKN, nm)(), "P3"; tq = 0.55)
end
let sol = solve(PROBS["P1"], DPRKN6())
    emit("DPRKN6", "P1", sol; tq = 0.55)
end
let sol = solve(PROBS["P2"], DPRKN6(), reltol = 1e-8, abstol = 1e-8, saveat = [0.5, 1.7, 2.9])
    println("# DPRKN6 P2 saveat: t = $(rv(sol.t)), u2 = $(rv(sol.u[2]))")
end
println("# IRKN (upstream code, IRKN3 in place: its out of place step fails)")
f1!(dv, v, u, p, t) = (dv .= -u .- u .* u .* u)
f5!(dv, v, u, p, t) = (dv .= -u .+ t)
f2!(du, v, u, p, t) = (du .= v)
const PIIP = Dict(
    "P1" => DynamicalODEProblem(f1!, f2!, [0.5, -0.3], [1.0, 0.2], (0.0, 1.0)),
    "P5" => DynamicalODEProblem(f5!, f2!, [0.5, -0.3], [1.0, 0.2], (0.0, 1.0)),
)
for case in ("P1", "P5")
    emit("IRKN3", case, solve(PIIP[case], IRKN3(), dt = 0.1, adaptive = false); dt = 0.1, tq = 0.55)
    fixed("IRKN4", IRKN4(), case)
end
