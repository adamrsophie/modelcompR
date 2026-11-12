# Reference values of DiffEqDevTools.jl 3.6.3 for the tests of the R port
# (tests/testthat/test-devtools-*.R). Run with
#   julia --project=tools/julia --startup-file=no tools/DiffEqDevTools-probe.jl
# Every line is "label value" with full precision values (repr). Cases where
# DiffEqDevTools or SciMLBase throw print the error, as the proof of the
# deliberate differences listed in the tests.
using DiffEqDevTools, OrdinaryDiffEq
using OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqTsit5, OrdinaryDiffEqVerner, OrdinaryDiffEqSSPRK,
      OrdinaryDiffEqSDIRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqBDF, OrdinaryDiffEqStabilizedRK,
      OrdinaryDiffEqExplicitRK
using OrdinaryDiffEqExplicitRK: constructTsit5ExplicitRK, ExplicitRK
using StochasticDiffEq, NonlinearSolve, DelayDiffEq
using DiffEqNoiseProcess: NoiseGrid
using SciMLBase, Statistics, LinearAlgebra
const DiffEqBase = DiffEqDevTools.DiffEqBase
const StructArray = DiffEqDevTools.StructArray
pr(name, x) = println(name, " ", repr(x))
err(f) = try
    f()
catch e
    sprint(showerror, e)
end

# ---------------------------------------------------------------------------
# Tableaus
tab = constructTsit5ExplicitRK()
A4 = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0.0]
rk4 = DiffEqBase.ExplicitRKTableau(A4, [0, 0.5, 0.5, 1.0], [1 / 6, 1 / 3, 1 / 3, 1 / 6], 4)
pr("res_none", residual_order_condition(tab, 3))
pr("res4", residual_order_condition(tab, 4))
pr("res_sum", [residual_order_condition(tab, i, +, abs) for i in 1:6])
pr("res_sum_emb", [residual_order_condition(tab, i, +, abs; embedded = true) for i in 1:5])
pr("res_rk4", [residual_order_condition(rk4, i, +, abs) for i in 1:5])
pr("res_rk4_sq", [residual_order_condition(rk4, i, +) for i in 1:5])
pr("check_tsit5", err(() -> check_tableau(tab)))
pr("check_rk4", check_tableau(rk4))
# embedded weights with sum(bhat c) = 1/2 but sum(bhat) = 1/2: Julia checks
# only the conditions of order adaptiveorder (all(f, 2) iterates over 2)
bad = DiffEqBase.ExplicitRKTableau(A4, [0, 0.5, 0.5, 1.0], [1 / 6, 1 / 3, 1 / 3, 1 / 6], 4;
    adaptiveorder = 2, αEEst = [0.0, 0.0, 0.0, 0.5])
pr("check_embedded", err(() -> check_tableau(bad)))
pr("res_bad_emb", [residual_order_condition(bad, i, +, abs; embedded = true) for i in 1:2])
for z in (-1.0, -2.5, 0.3)
    pr("sr $z", [stability_region(z, tab), stability_region(z, tab; embedded = true), stability_region(z, rk4)])
end
for z in (-1.0 + 2.0im, 0.5im, -3.0 + 0.1im)
    pr("src $z", [stability_region(z, tab), stability_region(z, rk4)])
end
pr("sreg", [stability_region(tab), stability_region(rk4)])
pr("isi_tab", err(() -> imaginary_stability_interval(tab)))
pr("isi_rk4", imaginary_stability_interval(rk4))
for alg in (Euler(), RK4(), Tsit5(), SSPRK33(), Midpoint(), Heun(), BS3(), DP5())
    pr("alg $(nameof(typeof(alg)))", [stability_region(-1.5, alg), stability_region(alg)])
end
for alg in (SSPRK33(), RK4(), Tsit5())
    pr("isialg $(nameof(typeof(alg)))", imaginary_stability_interval(alg))
    pr("srcalg $(nameof(typeof(alg)))", stability_region(-0.5 + 1.0im, alg))
end
pr("ie", stability_region(-2.0, ImplicitEuler()))
pr("ie_complex", err(() -> stability_region(-1.0 + 1.0im, ImplicitEuler())))
for alg in (Euler(), RK4(), SSPRK33(), Heun(), BS3(), ExplicitRK(tableau = tab))
    A, b, c = deduce_Butcher_tableau(alg)
    pr("deduce $(nameof(typeof(alg)))", (A, b, c))
end

# ---------------------------------------------------------------------------
# Convergence of ODEs: u1' = -u1^2, u2' = u1 u2, u1 = 1/(1+t), u2 = 2(1+t)
f(u, p, t) = [-u[1] * u[1], u[1] * u[2]]
an(u0, p, t) = [1 / (1 + t), 2 * (1 + t)]
prob = ODEProblem(ODEFunction(f; analytic = an), [1.0, 2.0], (0.0, 1.0))
dts = 1 ./ 2 .^ (3:6)
for alg in (Euler(), RK4(), Tsit5())
    sim = test_convergence(dts, prob, alg; dense_errors = true)
    n = nameof(typeof(alg))
    for k in (:final, :l2, :l∞, :L2, :L∞)
        pr("conv $n $k", sim.errors[k])
        pr("Oest $n $k", sim.𝒪est[k])
    end
end
pr("keys_default", sort(collect(keys(test_convergence(dts, prob, Euler()).errors))))
pr("keys_noevery", sort(collect(keys(test_convergence(dts, prob, Euler(); save_everystep = false).errors))))
sol = solve(prob, Euler(); dt = 1 / 16)
ref = solve(prob, RK4(); dt = 1 / 256, adaptive = false)
e1 = appxtrue(sol, ref)
for k in (:final, :l2, :l∞, :L2, :L∞)
    pr("appx $k", e1.errors[k])
end
pr("appxTS", appxtrue(sol, TestSolution(ref)).errors)
pr("appxTSnointerp", appxtrue(sol, TestSolution(ref.t, ref.u)).errors)
pr("appxTS3", appxtrue(sol, TestSolution(ref.t, ref.u, ref)).errors)
ref16 = solve(prob, RK4(); dt = 1 / 16, adaptive = false)
pr("appxTSsame", err(() -> appxtrue(sol, TestSolution(ref16.t, ref16.u)).errors))
psc = ODEProblem((u, p, t) -> -u * u, 1.0, (0.0, 1.0))
ssc = solve(psc, Euler(); dt = 1 / 16)
rsc = solve(psc, RK4(); dt = 1 / 16, adaptive = false)
pr("appxTSsame_scalar", appxtrue(ssc, TestSolution(rsc.t, rsc.u)).errors)
rsc2 = solve(psc, RK4(); dt = 1 / 16, adaptive = false, dense = false)
pr("appx_nondense_scalar", appxtrue(ssc, rsc2).errors)
prob2 = ODEProblem(f, [1.0, 2.0], (0.0, 1.0))
sim = analyticless_test_convergence(dts, prob2, Heun(), Dict(:alg => RK4(), :dt => 1 / 512, :adaptive => false))
for k in (:final, :l2, :l∞, :L2, :L∞)
    pr("aless $k", sim.errors[k])
    pr("alessO $k", sim.𝒪est[k])
end
probs = [remake(prob, tspan = (0.0, T)) for T in (0.5, 1.0)]
pr("probs final", test_convergence(probs, [0.5, 1.0], RK4(); dt = 1 / 8, adaptive = false).errors[:final])

# ---------------------------------------------------------------------------
# Random problems on a fixed Brownian path: u1 = 1 + t^2 + W1/2, u2 = 2 + t + W2/4
grid(h) = collect(0:h:1)
pathW(tg) = (W = [[0.3t - 0.7t * t + 0.01 * (k % 3), -0.2t + 0.5t * t * t - 0.02 * (k % 5)]
                  for (k, t) in enumerate(tg)]; W[1] = [0.0, 0.0]; W)
fs(u, p, t) = [2t, 1.0]
gs(u, p, t) = [0.5, 0.25]
ans_(u0, p, t, W) = [u0[1] + t * t + 0.5W[1], u0[2] + t + 0.25W[2]]
ans1(u0, p, t, W) = u0 + t * t + 0.5W
sdeprob(h) = SDEProblem(SDEFunction(fs, gs; analytic = ans_), [1.0, 2.0], (0.0, 1.0); noise = NoiseGrid(grid(h), pathW(grid(h))))
sdeprob1(h) = SDEProblem(SDEFunction((u, p, t) -> 2t, (u, p, t) -> 0.5; analytic = ans1), 1.0, (0.0, 1.0);
    noise = NoiseGrid(grid(h), [w[1] for w in pathW(grid(h))]))
sdts = [1 / 8, 1 / 16, 1 / 32]
sprob = sdeprob(1 / 64)
# SciMLBase takes sol.W.u[i], the ith grid value, for the ith step: compare
# with the errors at W(t[i])
for dt in sdts
    s = solve(sprob, EM(); dt = dt)
    pr("sde upstream $dt", s.errors)
    ua = [ans_(s.prob.u0, nothing, t, s.W(t)[1]) for t in s.t]
    d = reduce(vcat, s.u .- ua)
    pr("sde at W(t) $dt", Dict(:final => mean(abs.(s.u[end] - ua[end])), :l2 => sqrt(mean(d .^ 2)), :l∞ => maximum(abs.(d))))
end
sprob1 = sdeprob1(1 / 64)
sim = test_convergence(sdts, sprob1, EM(), EnsembleSerial(); trajectories = 2, dense_errors = true)
pr("sde1 L2", sim.errors[:L2]); pr("sde1 L∞", sim.errors[:L∞])
sim = test_convergence(sdts, sprob1, EM(), EnsembleSerial(); trajectories = 2, weak_dense_errors = true)
pr("sde weak_L2", sim.errors[:weak_L2]); pr("sde weak_L∞", sim.errors[:weak_L∞])
pr("sde vec dense", err(() -> solve(sprob, EM(); dt = 1 / 8, dense_errors = true).errors))
sim = test_convergence(sdts, EnsembleProblem(sprob; output_func = (sol, i) -> (sol.u[end][1], false)), EM(),
    EnsembleSerial(); trajectories = 2, expected_value = 2.0)
pr("ev weak_final", sim.errors[:weak_final]); pr("ev O", sim.𝒪est[:weak_final])
sim = test_convergence(sdts, EnsembleProblem(sprob; output_func = (sol, i) -> (sol.u[end], false)), EM(),
    EnsembleSerial(); trajectories = 2, expected_value = [2.0, 3.0])
pr("ev2 weak_l2", sim.errors[:weak_l2])
sim = test_convergence(sdts, sdeprob(1 / 8 / 4), EM(), EnsembleSerial(); trajectories = 2)
pr("sde_red keys", sort(collect(keys(sim.errors)), by = string))
pr("sde aligned", solve(sdeprob(1 / 8), EM(); dt = 1 / 8).errors)
pr("sde1 aligned", solve(sdeprob1(1 / 8), EM(); dt = 1 / 8, dense_errors = true).errors)
# DAE: du1 + u1^2 = 0, u1 u2 - 2 = 0
fd(du, u, p, t) = [du[1] + u[1] * u[1], u[1] * u[2] - 2.0]
and_(du0, u0, p, t) = [1 / (1 + t), 2 * (1 + t)]
dprob = DAEProblem(DAEFunction(fd; analytic = and_), [-1.0, 2.0], [1.0, 2.0], (0.0, 1.0); differential_vars = [true, false])
sim = test_convergence(sdts, dprob, DImplicitEuler(); dense_errors = true)
for k in (:final, :l2, :l∞, :L2)
    pr("dae $k", sim.errors[k])
end
# DDE: u'(t) = -u(t - 1), h = 1
fdd(u, h, p, t) = -h(p, t - 1)
hh(p, t) = 1.0
andd(u0, h, p, t) = t <= 1 ? 1 - t : -(2 * (t - 1) - (t * t - 1) / 2)
ddprob = DDEProblem(DDEFunction(fdd; analytic = andd), 1.0, hh, (0.0, 2.0); constant_lags = [1.0])
for dt in sdts
    pr("dde $dt", solve(ddprob, MethodOfSteps(Euler()); dt = dt).errors)
end

# ---------------------------------------------------------------------------
# Work precision
abstols = [1e-3, 1e-4, 1e-5]; reltols = [1e-3, 1e-4, 1e-5]
setups = [Dict(:alg => Euler(), :dts => [1 / 8, 1 / 16, 1 / 32]),
          Dict(:alg => RK4(), :dts => [1 / 4, 1 / 8, 1 / 16], :adaptive => false)]
wps = WorkPrecisionSet(prob, abstols, reltols, setups; numruns = 2, error_estimates = [:final, :l2, :L2])
for i in 1:2
    pr("wps names $i", wps.names[i])
    for k in (:final, :l2, :L2)
        pr("wps $i $k", getproperty(wps[i].errors, k))
    end
    pr("wps tags $i", wps[i].tags)
end
pr("avail", available_errors(wps))
wp = WorkPrecision(prob, RK4(), abstols, reltols, [1 / 4, 1 / 8, 1 / 16]; numruns = 2, error_estimate = :l∞,
    name = "R", adaptive = false)
pr("wp linf", wp.errors.l∞); pr("wp props", propertynames(wp.errors))
wp = WorkPrecision(prob, Tsit5(), abstols, reltols; numruns = 2)
pr("wp tsit5 final", wp.errors.final)
pr("wp tsit5 naccept", [s.naccept for s in wp.stats])
wp = WorkPrecision(prob2, Heun(), abstols, reltols, [1 / 8, 1 / 16, 1 / 32]; numruns = 2, appxsol = ref, error_estimate = :L2)
pr("wp appx L2", wp.errors.L2)
sh = Shootout(prob, [Dict(:alg => Euler()), Dict(:alg => RK4())]; dt = 1 / 16, numruns = 2, adaptive = false)
pr("shoot errors", sh.errors); pr("shoot names", sh.names)
pcsetups = [Dict(:alg => RK4(), :dts => [1 / 4, 1 / 8, 1 / 16], :adaptive => false),
            Dict(:alg => RK4(), :dts => [1 / 4, 1 / 8, 1 / 16], :adaptive => false, :prob_choice => 2)]
pprobs = [prob, remake(prob, tspan = (0.0, 0.5))]
pr("pc noappx", err(() -> WorkPrecisionSet(pprobs, abstols, reltols, pcsetups; numruns = 2).N))
refs = [solve(p, Tsit5(); abstol = 1e-12, reltol = 1e-12) for p in pprobs]
w = WorkPrecisionSet(pprobs, abstols, reltols, pcsetups; numruns = 2, appxsol = refs)
pr("pc 1", w[1].errors.final); pr("pc 2", w[2].errors.final)
dw = WorkPrecisionSet(ddprob, [1e-3, 1e-4], [1e-3, 1e-4], [Dict(:alg => MethodOfSteps(Euler()), :dts => [1 / 8, 1 / 16])];
    numruns = 2, error_estimate = :l2)
pr("dde names", dw.names); pr("dde l2", dw[1].errors.l2); pr("dde tags", dw[1].tags)
# nonlinear: x^3 + x = 2, y^2 = 3
fnl(u, p) = [u[1] * u[1] * u[1] + u[1] - p[1], u[2] * u[2] - p[2]]
nprob = NonlinearProblem(fnl, [1.0, 1.0], [2.0, 3.0])
nwp = WorkPrecisionSet(nprob, [1e-6, 1e-10], [1e-6, 1e-10], [Dict(:alg => NewtonRaphson())]; numruns = 2, error_estimate = :l2)
pr("nl l2", nwp[1].errors.l2)
pr("nl bad", err(() -> WorkPrecision(nprob, NewtonRaphson(), [1e-6], [1e-6]; numruns = 2, error_estimate = :final)))
# SDE and RODE sets on a NoiseGrid of step 1/16 with dt = 1/16
sprob16 = sdeprob(1 / 16)
ssetups = [Dict(:alg => EM(), :dts => [1 / 16, 1 / 16, 1 / 16]), Dict(:alg => EulerHeun(), :dts => [1 / 16, 1 / 16, 1 / 16])]
for ee in (:final, :l2, :weak_final)
    rw = WorkPrecisionSet(sprob16, abstols, reltols, ssetups; numruns = 2, numruns_error = 2, error_estimate = ee)
    for k in 1:2
        pr("rode $ee $k", getproperty(rw[k].errors, ee))
        pr("rode $ee $k N", (rw[k].N, length(getproperty(rw[k].errors, ee))))
        pr("rode $ee $k props", propertynames(rw[k].errors))
    end
    pr("rode $ee tags", rw[1].tags)
end
eprob = EnsembleProblem(sprob16; output_func = (sol, i) -> (sol.u[end][1], false))
eprob2 = EnsembleProblem(sprob16; output_func = (sol, i) -> (sol.u[end], false))
eset(p, ee; kw...) = WorkPrecisionSet(p, abstols, reltols, ssetups; numruns = 2, trajectories = 2,
    ensemblealg = EnsembleSerial(), error_estimate = ee, kw...)
ew = eset(eprob, :weak_final; expected_value = 2.0)
pr("ens weak_final", (ew[1].errors.weak_final, ew[2].errors.weak_final, ew[1].N))
pr("ens weak_l2 scalar", err(() -> eset(eprob, :weak_l2; expected_value = 2.0)[1].errors.weak_l2))
ew = eset(eprob2, :weak_l2; expected_value = [2.0, 3.0])
pr("ens vec weak_l2", (ew[1].errors.weak_l2, ew[2].errors.weak_l2))
pr("ens vec weak_final", err(() -> WorkPrecisionSet(eprob2, abstols, reltols, ssetups; numruns = 2, trajectories = 3,
    ensemblealg = EnsembleSerial(), expected_value = [2.0, 3.0], error_estimate = :weak_final)[1].errors.weak_final))
ew = eset(eprob, :weak_final; appxsol_setup = Dict(:alg => EM(), :dt => 1 / 16))
pr("ens appx", ew[2].errors.weak_final)
pr("ens strong", err(() -> eset(eprob, :final; expected_value = 2.0)))
pr("ens notimpl", err(() -> eset(eprob, :weak_l∞; expected_value = 2.0)))
# the BVP WorkPrecision reads errsol, which is assigned only when appxsol is
# given (benchmark.jl, WorkPrecision for AbstractBVProblem): the same scoping
function bvp_like(cur)
    local out
    for i in 1:1
        if cur !== nothing
            errsol = cur
            out = errsol
        else
            out = errsol
        end
    end
    return out
end
pr("bvp scoping", err(() -> bvp_like(nothing)))

# ---------------------------------------------------------------------------
# Tags and selection. The strings are fixed constructor calls written in this
# file (eval runs nothing else), so each label names the call.
for nm in ("Euler()", "RK4()", "Tsit5()", "BS3()", "Vern7()", "ImplicitEuler()", "Rodas4()", "Rodas5P()", "Rosenbrock23()",
           "FBDF()", "DFBDF()", "DImplicitEuler()", "Trapezoid()", "SSPRK33()", "SSPRK104()", "ROCK2()", "EM()", "EulerHeun()",
           "SRIW1()", "SOSRI()", "RKMil()", "SRA1()", "SOSRA()", "DRI1()", "RS1()", "ImplicitEM()", "SKenCarp()", "SplitEM()",
           "RandomEM()", "RandomTaylor15()", "AutoTsit5(Rosenbrock23())", "MethodOfSteps(Tsit5())", "NewtonRaphson()")
    pr("tags $nm", err(() -> auto_tags(eval(Meta.parse(nm)))))
end
# Boundary value algorithms: BoundaryValueDiffEq is not in tools/julia, so
# this block runs only in an environment that has it.
if Base.find_package("BoundaryValueDiffEq") !== nothing
    @eval using BoundaryValueDiffEq
    for nm in ("MIRK4()", "MIRK6()", "Shooting(Tsit5())", "MultipleShooting(4, Tsit5())")
        pr("tags $nm", Base.invokelatest(auto_tags, eval(Meta.parse(nm))))
    end
end
for t in (:rk, :sdirk, :adaptive, :order_5, :reference, :sundials, :autodiff_x, :dense, :stiff, :foo)
    pr("kind $t", tag_kind(t))
end
mk(name, errs, times, tags) = WorkPrecision(prob, abstols, reltols, StructArray((final = errs,)), times, nothing, nothing,
    name, :final, 3, tags)
set = WorkPrecisionSet([mk("A", [1e-2, 1e-4, 1e-6], [1e-3, 2e-3, 5e-3], [:rk, :explicit]),
                        mk("B", [1e-3, 1e-5, 1e-7], [2e-3, 4e-3, 9e-3], [:rk, :explicit]),
                        mk("C", [1e-2, NaN, 1e-6], [3e-3, 5e-3, 7e-3], [:rosenbrock, :implicit]),
                        mk("D", [2e-2, 3e-4, 1e-5], [1e-3, 2e-3, 3e-3], [:rosenbrock, :implicit, :reference])],
    4, abstols, reltols, prob, [Dict(:alg => Euler()) for _ in 1:4], ["A", "B", "C", "D"], :final, nothing)
pr("area", [wp_area(w) for w in set.wps])
pr("best rk", best_by_tag(set, :rk).names)
pr("best rk2", best_by_tag(set, :rk; n = 2).names)
pr("best ros", best_by_tag(set, :rosenbrock).names)
pr("best fam", best_of_families(set, [:rk, :rosenbrock]).names)
pr("best none", err(() -> best_of_families(set, [:sdirk])))
pr("best metric", err(() -> best_by_tag(set, :rk; metric = :time)))
pr("filter", filter_by_tags(set, :rk, :explicit).names)
pr("filter2", filter_by_tags(set, :implicit).names)
pr("exclude", exclude_by_tags(set, :reference).names)
pr("exclude2", exclude_by_tags(set, :rk, :reference).names)
pr("unique", unique_tags(set))
pr("get", get_tags(set))
ap = autoplot(set; reference_tags = :reference, best_n = 1)
for k in sort(collect(keys(ap)))
    pr("autoplot $k", ap[k].names)
end
ap = autoplot(set)
for k in sort(collect(keys(ap)))
    pr("autoplot0 $k", ap[k].names)
end
m = merge_wp_sets(set, best_by_tag(set, :rk))
pr("merge", (m.names, m.N))
pr("empty", filter_by_tags(set, :nothing).N)
v = with_autodiff_variants([Dict(:alg => Rodas5P()), Dict(:alg => Tsit5(), :tags => [:nonstiff])];
    ad_backends = [AutoFiniteDiff(), AutoForwardDiff()])
for d in v
    pr("adv", (get(d, :name, nothing), d[:tags], nameof(typeof(d[:alg]))))
end
