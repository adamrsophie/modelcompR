# Reference runs of OrdinaryDiffEqTaylorSeries (ExplicitTaylor2,
# ExplicitTaylor, ExplicitTaylorAdaptiveOrder) on out of place problems: the
# tab separated lines (method, case, saved points, nf, nreject, naccept,
# final state) are the numbers of tests/testthat/test-taylor.R. Run from the
# repository root:
#   julia --project=tools/julia tools/taylor-probe.jl [steps] [dense]
# With the manifest of tools/julia, TaylorDiff 0.3.5 fails to precompile:
# ChainRules 1.74.13 defines the tanh rule with a branch on
# abs(real(x)) < 1, which TaylorDiff's rule generation evaluates on a
# Symbolics.Num. The numbers were produced in a copy of the environment with
# ChainRules developed from a local copy whose tanh rule is
# `@scalar_rule tanh(x) 1 - Omega ^ 2` (Omega the result), which none of these
# problems use.
using OrdinaryDiffEqTaylorSeries, OrdinaryDiffEqCore, SciMLBase

# Lotka-Volterra and a forced oscillator with a rational forcing term
lv(u, p, t) = [p[1] * u[1] - u[1] * u[2], -p[2] * u[2] + u[1] * u[2]]
lvprob = ODEProblem(lv, [1.0, 1.0], (0.0, 2.0), [1.5, 3.0])
osc(u, p, t) = [u[2], -u[1] + p[1] * t * t / (1 + t)]
oscprob = ODEProblem(osc, [1.0, 0.0], (0.0, 3.0), [0.5])
cblv = DiscreteCallback((u, t, integ) -> t == 1.0, integ -> (integ.u = integ.u * 0.5))

function report(nm, sc, s)
    st = s.stats
    println(nm, "\t", sc, "\t", length(s.t), "\t", st.nf, "\t", st.nreject, "\t", st.naccept, "\t",
        join(repr.(s.u[end]), ","))
end

want(g) = isempty(ARGS) || "all" in ARGS || g in ARGS
if want("steps")
report("ExplicitTaylor2", "lv", solve(lvprob, ExplicitTaylor2(), dt = 0.01))
report("ExplicitTaylor2", "osc", solve(oscprob, ExplicitTaylor2(), dt = 0.01))
report("ExplicitTaylor2", "cb", solve(lvprob, ExplicitTaylor2(), dt = 0.01, tstops = [1.0], callback = cblv))
for P in (1, 2, 4, 8)
    report("ExplicitTaylorOrder$P", "lv", solve(lvprob, ExplicitTaylor(order = Val(P)), reltol = 1e-8, abstol = 1e-10))
    report("ExplicitTaylorOrder$P", "osc", solve(oscprob, ExplicitTaylor(order = Val(P)), reltol = 1e-8, abstol = 1e-10))
end
report("ExplicitTaylorOrder4", "cb", solve(lvprob, ExplicitTaylor(order = Val(4)), reltol = 1e-8, abstol = 1e-10, tstops = [1.0], callback = cblv))
report("ExplicitTaylorOrder4", "fixed", solve(lvprob, ExplicitTaylor(order = Val(4)), dt = 0.05, adaptive = false))
report("ExplicitTaylorAdaptiveOrder", "lv", solve(lvprob, ExplicitTaylorAdaptiveOrder(), reltol = 1e-8, abstol = 1e-10))
report("ExplicitTaylorAdaptiveOrder", "osc", solve(oscprob, ExplicitTaylorAdaptiveOrder(), reltol = 1e-8, abstol = 1e-10))
report("ExplicitTaylorAdaptiveOrder", "cb", solve(lvprob, ExplicitTaylorAdaptiveOrder(), reltol = 1e-8, abstol = 1e-10, tstops = [1.0], callback = cblv))
report("ExplicitTaylorAdaptiveOrder_4_9", "lv", solve(lvprob, ExplicitTaylorAdaptiveOrder(min_order = Val(4), max_order = Val(9)), reltol = 1e-12, abstol = 1e-12))
report("ExplicitTaylorAdaptiveOrder_2_5", "osc", solve(oscprob, ExplicitTaylorAdaptiveOrder(min_order = Val(2), max_order = Val(5)), reltol = 1e-6, abstol = 1e-8))
end
# Dense output between the steps: method, "dense", case, time, state
function dense(nm, sc, s, ts)
    for t in ts
        println(nm, "\tdense\t", sc, "\t", repr(t), "\t", join(repr.(s(t)), ","))
    end
end
if want("dense")
    dense("ExplicitTaylor2", "lv", solve(lvprob, ExplicitTaylor2(), dt = 0.01), (0.373, 1.234))
    dense("ExplicitTaylorOrder4", "lv", solve(lvprob, ExplicitTaylor(order = Val(4)), reltol = 1e-8, abstol = 1e-10), (0.373, 1.234))
    dense("ExplicitTaylorOrder4", "osc", solve(oscprob, ExplicitTaylor(order = Val(4)), reltol = 1e-8, abstol = 1e-10), (0.373, 2.345))
    dense("ExplicitTaylorAdaptiveOrder", "lv", solve(lvprob, ExplicitTaylorAdaptiveOrder(), reltol = 1e-8, abstol = 1e-10), (0.373, 1.234))
    dense("ExplicitTaylorAdaptiveOrder", "osc", solve(oscprob, ExplicitTaylorAdaptiveOrder(), reltol = 1e-8, abstol = 1e-10), (0.373, 2.345))
end
