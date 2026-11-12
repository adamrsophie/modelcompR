# Reference values of the discontinuity detection of OrdinaryDiffEqCore
# (disco.jl, handle_disco_accept! / handle_disco_reject!) for
# tests/testthat/test-disco.R. Each line: case, dd, nt, nf, naccept,
# nreject, the saved times and the final state.
using OrdinaryDiffEqTsit5, OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqVerner,
      OrdinaryDiffEqLowStorageRK, OrdinaryDiffEqStabilizedRK, OrdinaryDiffEqCore, SciMLBase
import NonlinearSolve
const C = OrdinaryDiffEqCore
f1(u, p, t) = [u[2] < 0.5 ? u[1] : -u[1] * u[1], 1.0]
f2(u, p, t) = [u[2] < 0.5 ? u[1] : (u[2] < 0.7 ? -u[1] * u[1] : 2 * u[1]), 1.0]
cbA = ContinuousCallback((u, t, integ) -> u[2] - 0.5, integ -> nothing)
cbB = ContinuousCallback((u, t, integ) -> u[2] - 0.52, integ -> nothing)
cbN = ContinuousCallback((u, t, integ) -> u[2] - 0.5, integ -> nothing; maybe_discontinuity = false)
vcb = VectorContinuousCallback((out, u, t, integ) -> (out[1] = u[2] - 0.5; out[2] = u[2] - 0.7; nothing),
                               (integ, idx) -> nothing, 2)
ctrl(alg::RDPK3Sp35, dd) = C.PIDController(Float64, alg; beta = (0.64, -0.31, 0.04), discontinuity_detection = dd)
ctrl(alg::RKC, dd) = C.PredictiveController(alg; discontinuity_detection = dd)
ctrl(alg, dd) = C.PIController(alg; discontinuity_detection = dd)
cases = [
    ("Tsit5", f1, Tsit5(), cbA), ("BS3", f1, BS3(), cbA), ("Vern7", f1, Vern7(), cbA),
    ("RDPK3Sp35", f1, RDPK3Sp35(), cbA), ("RKC", f1, RKC(), cbA),
    ("Tsit5_vector", f2, Tsit5(), vcb), ("Tsit5_BA", f1, Tsit5(), CallbackSet(cbB, cbA)),
    ("Tsit5_nodisco", f1, Tsit5(), cbN),
]
for (name, f, alg, cb) in cases, dd in (false, true)
    prob = ODEProblem(f, [1.0, 0.0], (0.0, 1.0))
    s = solve(prob, alg; callback = cb, controller = ctrl(alg, dd), reltol = 1e-6, abstol = 1e-8)
    println(name, "\t", dd, "\t", length(s.t), "\t", s.stats.nf, "\t", s.stats.naccept, "\t", s.stats.nreject,
            "\t", join(repr.(s.t), ","), "\t", join(repr.(s.u[end]), ","))
end
# OrdinaryDiffEq brackets the second callback on [0, theta] without checking
# the sign at theta: with cbA before cbB the later root of cbB gives theta 0
prob = ODEProblem(f1, [1.0, 0.0], (0.0, 1.0))
s = solve(prob, Tsit5(); callback = CallbackSet(cbA, cbB), controller = ctrl(Tsit5(), true), reltol = 1e-6, abstol = 1e-8)
println("Tsit5_AB_upstream\ttrue\t", length(s.t), "\t", s.stats.nf, "\t", s.stats.naccept, "\t", s.stats.nreject,
        "\t", join(repr.(s.t), ","), "\t", join(repr.(s.u[end]), ","))
