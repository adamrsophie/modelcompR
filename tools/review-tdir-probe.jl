# Reference behavior of upstream for four review findings, for
# tests/testthat/test-sde-tdir.R, test-exponential-rk.R and test-rosenbrock.R:
#  A. Exp4, EPIRK* and EXPRB53s3 integrating backward in time, and
#     phiv_timestep with negative times;
#  B. ROS3P, ROS3PR, ROS34PW1a and ROS34PW1b (Rodas3 as a control) on linear
#     problems at abstol = reltol = 1e-8;
#  C. adaptive RKMil and RKMilGeneral at the default tolerances (RKMilCommute
#     and SRIW1 as controls), with random noise and on a fixed grid;
#  D. SDE problems backward in time (test/tdir_tests.jl) on a fixed reversed
#     NoiseGrid for every algorithm, and zero length spans of SDE and DDE
#     problems.
# Every case runs in try/catch; an error prints its message.
# Run: julia --project=tools/julia --startup-file=no tools/review-tdir-probe.jl
using StochasticDiffEq, DiffEqNoiseProcess, DelayDiffEq, LinearAlgebra, Random
using OrdinaryDiffEqExponentialRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqTsit5
using OrdinaryDiffEqDifferentiation: AutoFiniteDiff
import Pkg
const ExponentialUtilities = OrdinaryDiffEqExponentialRK.ExponentialUtilities

for (k, v) in sort(collect(Pkg.dependencies()); by = x -> x[2].name)
    v.name in ("StochasticDiffEq", "StochasticDiffEqCore", "StochasticDiffEqLowOrder",
        "StochasticDiffEqMilstein", "DiffEqNoiseProcess", "DelayDiffEq", "OrdinaryDiffEqCore",
        "OrdinaryDiffEqExponentialRK", "ExponentialUtilities", "OrdinaryDiffEqRosenbrock") &&
        println("version ", v.name, " ", v.version)
end

out(name, key, x) = println(name, " ", key, " ", join(repr.(Float64.(collect(Iterators.flatten(x)))), " "))
outs(name, key, s) = println(name, " ", key, " ", s)
function go(name, thunk)
    try
        sol = thunk()
        outs(name, "retcode", string(sol.retcode))
        out(name, "t", sol.t)
        out(name, "uend", sol.u[end])
        st = sol.stats
        out(name, "stats", [st.nf, st.naccept, st.nreject])
        return sol
    catch e
        outs(name, "error", replace(first(sprint(showerror, e), 300), '\n' => ' '))
        return nothing
    end
end

# ---- A. exponential integrators backward in time ---------------------------
A2 = [-1.0 0.3; 0.2 -2.0]
fe(u, p, t) = A2 * u + [0.1 * u[1]^2, -0.05 * t]
je(u, p, t) = A2 + [0.2*u[1] 0.0; 0.0 0.0]
pe_fwd = ODEProblem(ODEFunction(fe; jac = je), [1.0, 0.5], (0.0, 1.0))
pe_bwd = ODEProblem(ODEFunction(fe; jac = je), [1.0, 0.5], (1.0, 0.0))
pl_bwd = ODEProblem(ODEFunction((u, p, t) -> -u; jac = (u, p, t) -> fill(-1.0, 1, 1)), [exp(-1.0)], (1.0, 0.0))
for (nm, A) in [("Exp4", Exp4), ("EPIRK4s3A", EPIRK4s3A), ("EPIRK4s3B", EPIRK4s3B),
        ("EPIRK5s3", EPIRK5s3), ("EXPRB53s3", EXPRB53s3), ("EPIRK5P1", EPIRK5P1), ("EPIRK5P2", EPIRK5P2)]
    go("expo_fwd_$nm", () -> solve(pe_fwd, A(); dt = 0.1))
    for rep in 1:2
        go("expo_bwd$(rep)_$nm", () -> solve(pe_bwd, A(); dt = -0.1))
        go("expo_lin$(rep)_$nm", () -> solve(pl_bwd, A(); dt = -0.01))
    end
    go("expo_bwdpos_$nm", () -> solve(pe_bwd, A(); dt = 0.1))
    go("expo_bwdnoadapt_$nm", () -> solve(pe_bwd, A(adaptive_krylov = false); dt = -0.1))
end
# phiv_timestep(ts, A, B): B has the columns b_0, b_1 of the phi functions
Bm = [1.0 0.2; 0.5 -0.1]
for (nm, ts) in [("neg1", [-0.5, -1.0]), ("neg2", [-0.5, -1.0]), ("pos", [0.5, 1.0])]
    try
        out("phiv_timestep_$nm", "U", ExponentialUtilities.phiv_timestep(ts, A2, Bm))
    catch e
        outs("phiv_timestep_$nm", "error", replace(first(sprint(showerror, e), 300), '\n' => ' '))
    end
end

# ---- B. Rosenbrock methods on linear problems -------------------------------
pr_dec = ODEProblem(ODEFunction((u, p, t) -> -u; jac = (u, p, t) -> fill(-1.0, 1, 1)), [1.0], (0.0, 1.0))
pr_osc = ODEProblem(ODEFunction((u, p, t) -> [u[2], -u[1]]; jac = (u, p, t) -> [0.0 1.0; -1.0 0.0]),
    [1.0, 0.0], (0.0, 1.0))
pr_sq = ODEProblem(ODEFunction((u, p, t) -> -u .^ 2; jac = (u, p, t) -> fill(-2.0 * u[1], 1, 1)), [1.0], (0.0, 1.0))
for (nm, A) in [("ROS3P", ROS3P), ("ROS3PR", ROS3PR), ("ROS34PW1a", ROS34PW1a),
        ("ROS34PW1b", ROS34PW1b), ("Rodas3", Rodas3), ("ROS3", ROS3)]
    s = go("ros_dec_$nm", () -> solve(pr_dec, A(); abstol = 1e-8, reltol = 1e-8))
    s === nothing || out("ros_dec_$nm", "err", [abs(s.u[end][1] - exp(-1.0))])
    s = go("ros_osc_$nm", () -> solve(pr_osc, A(); abstol = 1e-8, reltol = 1e-8))
    s === nothing || out("ros_osc_$nm", "err", [maximum(abs.(s.u[end] .- [cos(1.0), -sin(1.0)]))])
    s = go("ros_sq_$nm", () -> solve(pr_sq, A(); abstol = 1e-8, reltol = 1e-8))
    s === nothing || out("ros_sq_$nm", "err", [abs(s.u[end][1] - 0.5)])
end

# ---- C. adaptive RKMil --------------------------------------------------------
z(k) = mod(37 * k, 101) / 50 - 1
fm(u, p, t) = -u
gm(u, p, t) = 0.1 * u
for (nm, A) in [("RKMil", RKMil), ("RKMilGeneral", RKMilGeneral), ("RKMilCommute", RKMilCommute), ("SRIW1", SRIW1)]
    for seed in 1:5
        go("rkmil_seed$(seed)_$nm", () -> solve(SDEProblem(fm, gm, 1.0, (0.0, 1.0)), A(); seed = seed, maxiters = 2000))
    end
    go("rkmil_loose_$nm", () -> solve(SDEProblem(fm, gm, 1.0, (0.0, 1.0)), A(); seed = 7,
        abstol = 0.1, reltol = 0.1, maxiters = 10^5))
    hg = 1 / 4096
    tg = collect(0:4096) .* hg
    Wg = [0.0]; Zg = [0.0]
    for k in 1:4096
        push!(Wg, Wg[end] + sqrt(hg) * z(k))
        push!(Zg, Zg[end] + sqrt(hg) * z(k + 7919))
    end
    go("rkmil_grid_$nm", () -> solve(SDEProblem(fm, gm, 1.0, (0.0, 1.0); noise = NoiseGrid(tg, Wg, Zg)), A();
        maxiters = 2000))
end

# ---- D. backward in time ------------------------------------------------------
fd(u, p, t) = -0.5 .* u .+ 0.1 * t
gd(u, p, t) = 0.3 .* u
hb = 1 / 64
tb = collect(1.0 .- (0:64) .* hb)
Wb = [zeros(2)]; Zb = [zeros(2)]
for k in 1:64
    push!(Wb, Wb[end] .+ sqrt(hb) .* [z(2k - 1), z(2k)])
    push!(Zb, Zb[end] .+ sqrt(hb) .* [z(7919 + 2k - 1), z(7919 + 2k)])
end
pbwd() = SDEProblem(fd, gd, [1.0, 0.5], (1.0, 0.0); noise = NoiseGrid(tb, Wb, Zb))
algnames = ["EM", "EulerHeun", "LambaEM", "LambaEulerHeun", "RKMil", "RKMilCommute", "RKMilGeneral",
    "SimplifiedEM", "SRI", "SRIW1", "SRIW2", "SOSRI", "SOSRI2", "SRA", "SRA1", "SRA2", "SRA3", "SOSRA",
    "SOSRA2", "ImplicitEM", "ImplicitEulerHeun", "ImplicitRKMil", "STrapezoid", "SImplicitMidpoint",
    "ISSEM", "ISSEulerHeun", "SKenCarp", "SROCK1", "SROCK2", "SROCKEM", "SROCKC2", "SKSROCK",
    "TangXiaoSROCK2", "KomBurSROCK2", "WangLi3SMil_A", "WangLi3SMil_B", "WangLi3SMil_C",
    "WangLi3SMil_D", "WangLi3SMil_E", "WangLi3SMil_F", "DRI1", "DRI1NM", "RI1", "RI3", "RI5", "RI6",
    "RDI1WM", "RDI2WM", "RDI3WM", "RDI4WM", "RS1", "RS2", "PL1WM", "PL1WMA", "NON", "NON2", "COM",
    "SIEA", "SIEB", "SMEA", "SMEB", "W2Ito1", "IRI1"]
for nm in algnames
    alg = try
        getfield(StochasticDiffEq, Symbol(nm))()
    catch e
        outs("bwd_$nm", "error", "constructor: " * replace(first(sprint(showerror, e), 200), '\n' => ' '))
        continue
    end
    go("bwd_$nm", () -> solve(pbwd(), alg; dt = -hb, adaptive = false))
end
# the step sign follows the span (a positive dt on a reversed span), events
go("bwd_EM_posdt", () -> solve(pbwd(), EM(); dt = hb))
go("bwd_EM_tstops", () -> solve(pbwd(), EM(); dt = -1 / 16, tstops = [0.33]))
go("bwd_SRIW1_tstops", () -> solve(pbwd(), SRIW1(); dt = -1 / 16, adaptive = false, tstops = [0.33]))
# adaptive solves backward on the grid
for nm in ["SRIW1", "SOSRI", "SRA1", "LambaEM", "RKMil", "RKMilCommute", "ImplicitEM", "EM"]
    go("bwdad_$nm", () -> solve(pbwd(), getfield(StochasticDiffEq, Symbol(nm))(); dt = -hb))
end
# the reversed noise of tdir_tests.jl, on a fixed forward grid
tf = collect(0:100) .* 0.01
Wf = [0.0]; for k in 1:100; push!(Wf, Wf[end] + 0.1 * z(k + 11)); end
pf = SDEProblem((u, p, t) -> u, (u, p, t) -> u, 0.5, (0.0, 1.0); noise = NoiseGrid(tf, Wf))
sf = go("rev_fwd_EulerHeun", () -> solve(pf, EulerHeun(); dt = 0.01))
if sf !== nothing
    out("rev_fwd_EulerHeun", "u", sf.u)
    pr = SDEProblem((u, p, t) -> u, (u, p, t) -> u, sf.u[end], (1.0, 0.0); noise = NoiseGrid(reverse(tf), reverse(Wf)))
    sr = go("rev_bwd_EulerHeun", () -> solve(pr, EulerHeun(); dt = 0.01))
    sr === nothing || out("rev_bwd_EulerHeun", "u", sr.u)
end
# random noise backward: the structure of tdir_tests.jl
Random.seed!(100)
go("bwd_random_EM", () -> solve(SDEProblem(fd, gd, [1.0, 0.5], (1.0, 0.0)), EM(); dt = -1 / 16, tstops = [0.33]))
go("bwd_random_SRIW1", () -> solve(SDEProblem(fd, gd, [1.0, 0.5], (1.0, 0.0)), SRIW1(); tstops = [0.33]))
go("bwd_random_scalar_EM", () -> solve(SDEProblem((u, p, t) -> 1.01u, (u, p, t) -> 0.87u, 0.5, (1.0, 0.0)), EM(); dt = -1 / 16))

# ---- zero length spans ---------------------------------------------------------
for (nm, A) in [("EM", EM), ("SRIW1", SRIW1), ("SOSRI", SOSRI), ("LambaEM", LambaEM)]
    s = go("zero_sde_$nm", () -> solve(SDEProblem(fd, gd, [1.0, 0.5], (0.0, 0.0)), A(); dt = 0.1))
    try; s === nothing || out("zero_sde_$nm", "W", s.W === nothing ? [NaN] : s.W.t); catch; end
end
go("zero_sde_default", () -> solve(SDEProblem(fd, gd, [1.0, 0.5], (0.0, 0.0))))
go("zero_sde_nodt", () -> solve(SDEProblem(fd, gd, [1.0, 0.5], (0.0, 0.0)), SRIW1()))
go("zero_sde_init", () -> (ig = init(SDEProblem(fd, gd, [1.0, 0.5], (0.0, 0.0)), EM(); dt = 0.1); solve!(ig)))
hdde(p, t) = [1.0]
fdde(u, h, p, t) = -h(p, t - 1.0)
go("zero_dde", () -> solve(DDEProblem(fdde, [1.0], hdde, (0.0, 0.0); constant_lags = [1.0]), MethodOfSteps(Tsit5())))
go("zero_dde_dt", () -> solve(DDEProblem(fdde, [1.0], hdde, (0.0, 0.0); constant_lags = [1.0]), MethodOfSteps(Tsit5()); dt = 0.1))
go("zero_ode", () -> solve(ODEProblem((u, p, t) -> -u, [1.0], (0.0, 0.0)), Tsit5()))
