# Reference runs of the out of place Rosenbrock methods of OrdinaryDiffEqRosenbrock
# (tests/testthat/test-rosenbrock-family.R). Prints one tab separated line per
# run: method, case, number of saved points, nf, njacs, nw, nsolve, naccept,
# nreject, the final state and dense output values. Run from the repository root:
#   julia --project=tools/julia tools/rosenbrock-probe.jl > rosenbrock.tsv
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqCore, SciMLBase, LinearAlgebra
const AutoFiniteDiff = OrdinaryDiffEqRosenbrock.AutoFiniteDiff

# Robertson
rober(u, p, t) = [-p[1] * u[1] + p[2] * u[2] * u[3], p[1] * u[1] - p[2] * u[2] * u[3] - p[3] * u[2] * u[2], p[3] * u[2] * u[2]]
rjac(u, p, t) = [-p[1] p[2]*u[3] p[2]*u[2]; p[1] -p[2]*u[3]-2*p[3]*u[2] -p[2]*u[2]; 0.0 2*p[3]*u[2] 0.0]
# Robertson as a DAE (singular mass matrix)
roberdae(u, p, t) = [-p[1] * u[1] + p[2] * u[2] * u[3], p[1] * u[1] - p[2] * u[2] * u[3] - p[3] * u[2] * u[2], u[1] + u[2] + u[3] - 1]
rdjac(u, p, t) = [-p[1] p[2]*u[3] p[2]*u[2]; p[1] -p[2]*u[3]-2*p[3]*u[2] -p[2]*u[2]; 1.0 1.0 1.0]
# a mildly stiff non autonomous polynomial system
pf(u, p, t) = [-50 * (u[1] - u[2] * u[2]) + t, -u[2] + 0.5 * u[1] * u[3], -2 * u[3] + u[1] * u[1] * u[1] - 0.1 * t * t]
pjac(u, p, t) = [-50.0 100*u[2] 0.0; 0.5*u[3] -1.0 0.5*u[1]; 3*u[1]*u[1] 0.0 -2.0]
ptgrad(u, p, t) = [1.0, 0.0, -0.2 * t]
# a semi explicit index 1 DAE: two differential variables, one algebraic
df(u, p, t) = [-u[1] + u[3], u[1] - 2 * u[2] + 0.1 * t, u[3] - u[2] * u[2] - 0.5 * u[1]]
djac(u, p, t) = [-1.0 0.0 1.0; 1.0 -2.0 0.0; -0.5 -2*u[2] 1.0]
dtgrad(u, p, t) = [0.0, 0.1, 0.0]
const M3 = [1.0 0 0; 0 1 0; 0 0 0]

prob_rob = ODEProblem(ODEFunction(rober; jac = rjac), [1.0, 0.0, 0.0], (0.0, 10.0), [0.04, 1e4, 3e7])
prob_robdae = ODEProblem(ODEFunction(roberdae; jac = rdjac, mass_matrix = M3), [1.0, 0.0, 0.0], (0.0, 10.0), [0.04, 1e4, 3e7])
prob_p = ODEProblem(ODEFunction(pf; jac = pjac, tgrad = ptgrad), [1.0, 0.5, 0.2], (0.0, 2.0))
prob_pfd = ODEProblem(pf, [1.0, 0.5, 0.2], (0.0, 2.0))
prob_d = ODEProblem(ODEFunction(df; jac = djac, tgrad = dtgrad, mass_matrix = M3), [1.0, 0.5, 0.75], (0.0, 2.0))
prob_dfd = ODEProblem(ODEFunction(df; mass_matrix = M3), [1.0, 0.5, 0.75], (0.0, 2.0))
cb = DiscreteCallback((u, t, integ) -> t == 0.5, integ -> (integ.u = integ.u * 0.5))

const METHODS = [:Rosenbrock23, :Rosenbrock32, :ROS3P, :Rodas3, :Rodas3d, :Rodas23W, :Rodas3P, :Rodas4,
    :Rodas42, :Rodas4P, :Rodas4P2, :Rodas4PW, :Rodas5, :Rodas5P, :Rodas5Pe, :Rodas5Pr, :Rodas6P,
    :ROS2, :ROS2PR, :ROS2S, :ROS3, :ROS3PR, :Scholz4_7, :ROS34PW1a, :ROS34PW1b, :ROS34PW2, :ROS34PW3,
    :ROS34PRw, :ROS3PRL, :ROS3PRL2, :ROK4a, :RosShamp4, :Veldd4, :Velds4, :GRK4T, :GRK4A, :Ros4LStab, :Tsit5DA]
mk(nm; kw...) = getfield(OrdinaryDiffEqRosenbrock, nm)(; kw...)

fmt(v) = join(repr.(v), ",")
function out(nm, case, s, ts = Float64[])
    st = s.stats
    interp = join([fmt(s(tt)) for tt in ts], ";")
    println(nm, "\t", case, "\t", s.retcode, "\t", length(s.t), "\t", st.nf, "\t", st.njacs, "\t", st.nw, "\t",
        st.nsolve, "\t", st.naccept, "\t", st.nreject, "\t", fmt(s.u[end]), "\t", interp)
    flush(stdout)
end
function run(nm, case, f)
    try
        out(nm, case, f()...)
    catch e
        println(nm, "\t", case, "\tERROR\t", replace(sprint(showerror, e), r"\s+" => " ")[1:min(end, 200)])
        flush(stdout)
    end
end

const TS = [0.123, 0.777, 1.5]
for nm in METHODS
    tol = (reltol = 1e-5, abstol = 1e-7)
    nm == :Tsit5DA || run(nm, "rob", () -> (solve(prob_rob, mk(nm); reltol = 1e-6, abstol = 1e-8), [0.5, 7.0]))
    run(nm, "plain", () -> (solve(prob_p, mk(nm); tol...), TS))
    run(nm, "tstop", () -> (solve(prob_p, mk(nm); tstops = [0.37], tol...), TS))
    run(nm, "cb", () -> (solve(prob_p, mk(nm); tstops = [0.5], callback = cb, tol...), [0.3, 0.7]))
    run(nm, "nodense", () -> (solve(prob_p, mk(nm); save_everystep = false, tol...),))
    run(nm, "saveat", () -> (solve(prob_p, mk(nm); saveat = 0.25, tol...),))
    run(nm, "fd", () -> (solve(prob_pfd, mk(nm; autodiff = AutoFiniteDiff()); tol...), TS))
    run(nm, "dae", () -> (solve(prob_d, mk(nm); tol...), TS))
    run(nm, "daefd", () -> (solve(prob_dfd, mk(nm; autodiff = AutoFiniteDiff()); tol...), TS))
    nm == :Tsit5DA || run(nm, "robdae", () -> (solve(prob_robdae, mk(nm); reltol = 1e-6, abstol = 1e-8), [0.5, 7.0]))
    run(nm, "fixed", () -> (solve(prob_p, mk(nm); dt = 0.05, adaptive = false), TS))
    run(nm, "fixeddae", () -> (solve(prob_d, mk(nm); dt = 0.05, adaptive = false), TS))
end
run(:RosenbrockW6S4OS, "fixed", () -> (solve(prob_p, RosenbrockW6S4OS(); dt = 0.05), TS))
run(:RosenbrockW6S4OS, "fixeddae", () -> (solve(prob_d, RosenbrockW6S4OS(); dt = 0.05), TS))
run(:RosenbrockW6S4OS, "fixedfd", () -> (solve(prob_pfd, RosenbrockW6S4OS(autodiff = AutoFiniteDiff()); dt = 0.05), TS))
# Jacobian reuse knobs of the W-methods
run(:ROS34PW2_age3, "rob", () -> (solve(prob_rob, ROS34PW2(max_jac_age = 3); reltol = 1e-6, abstol = 1e-8),))
run(:ROS34PW2_tol, "rob", () -> (solve(prob_rob, ROS34PW2(jac_reuse_gamma_tol = 0.5); reltol = 1e-6, abstol = 1e-8),))
run(:Rosenbrock23_age20, "rob", () -> (solve(prob_rob, Rosenbrock23(max_jac_age = 20); reltol = 1e-6, abstol = 1e-8),))
# a step limiter (FSAL and non FSAL methods)
lim(u, integ, p, t) = (u .= max.(u, 0.0); nothing)
run(:ROS3P_lim, "plain", () -> (solve(prob_p, ROS3P(step_limiter! = lim); reltol = 1e-5, abstol = 1e-7),))
run(:Rodas4_lim, "plain", () -> (solve(prob_p, Rodas4(step_limiter! = lim); reltol = 1e-5, abstol = 1e-7),))
run(:Rosenbrock23_lim, "plain", () -> (solve(prob_p, Rosenbrock23(step_limiter! = lim); reltol = 1e-5, abstol = 1e-7),))
# the tableau constructors
for c in (:Rodas4Tableau, :ROS3PRodasTableau, :RosenbrockW6S4OSRodasTableau)
    tab = getfield(OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockTableaus, c)(Float64, Float64)
    println(c, "\ttab\t", fmt(vec(tab.A)), "\t", repr(tab.gamma), "\t", tab.btilde === nothing ? "nothing" : fmt(tab.btilde))
end
const TOL = (reltol = 1e-5, abstol = 1e-7)
# the default autodiff (AutoForwardDiff) without a Jacobian function; R
# differentiates by finite differences, so only the step counts agree
for nm in (:Rodas5P, :Rosenbrock23, :ROS34PW2)
    run(nm, "ad", () -> (solve(prob_pfd, mk(nm); TOL...), TS))
    run(nm, "daead", () -> (solve(prob_dfd, mk(nm); TOL...), TS))
end
# backward integration: diffdir makes the finite difference in t one sided
# toward t0
prob_bjac = ODEProblem(ODEFunction(df; jac = djac, tgrad = dtgrad), [1.0, 0.5, 0.75], (2.0, 0.0))
prob_bfd = ODEProblem(df, [1.0, 0.5, 0.75], (2.0, 0.0))
prob_bdaefd = ODEProblem(ODEFunction(df; mass_matrix = M3), [1.0, 0.5, 0.75], (2.0, 0.0))
for nm in (:Rodas5P, :Rosenbrock23, :ROS34PW2, :Rodas4)
    run(nm, "bjac", () -> (solve(prob_bjac, mk(nm); TOL...), reverse(TS)))
    run(nm, "bfd", () -> (solve(prob_bfd, mk(nm; autodiff = AutoFiniteDiff()); TOL...), reverse(TS)))
    run(nm, "bdaefd", () -> (solve(prob_bdaefd, mk(nm; autodiff = AutoFiniteDiff()); TOL...), reverse(TS)))
    run(nm, "bfixedfd", () -> (solve(prob_bfd, mk(nm; autodiff = AutoFiniteDiff()); dt = 0.05, adaptive = false),
        reverse(TS)))
end
# HybridExplicitImplicitRK ignores its tableau (tabtype = Tsit5DATableau)
for (lab, tab) in (("Hybrid_Tsit5DATableau", OrdinaryDiffEqRosenbrock.Tsit5DATableau(Float64, Float64)),
        ("Hybrid_Rodas4Tableau", OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockTableaus.Rodas4Tableau(Float64, Float64)))
    run(lab, "dae", () -> (solve(prob_dfd, HybridExplicitImplicitRK(tab; order = 5); TOL...), TS))
end
# lu(W; check = false) \ b (OpenBLAS) and ODE_DEFAULT_NORM on random data,
# printed with %a: the LU and NORM data of the test file
using Printf, Random
hx(v) = join([@sprintf("%a", x) for x in vec(v)], ",")
for (seed, sizes, reps) in ((7, (2, 3, 4), 150), (7, (5, 6, 7, 8, 10, 12, 16, 20), 60))
    Random.seed!(seed)
    for n in sizes, r in 1:reps
        W = randn(n, n) .* exp.(3 .* randn(n, n))
        b = randn(n) .* exp.(2 .* randn(n))
        println("LU\t", n, "\t", hx(W), "\t", hx(b), "\t", hx(lu(W; check = false) \ b))
    end
end
Random.seed!(11)
for n in vcat(1:40, [47, 64, 100]), r in 1:30
    u = randn(n) .* exp.(2 .* randn(n))
    println("NORM\t", n, "\t", hx(u), "\t", @sprintf("%a", OrdinaryDiffEqCore.DiffEqBase.ODE_DEFAULT_NORM(u, 0.0)))
end
