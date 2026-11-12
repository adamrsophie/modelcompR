# Reference values for tests/testthat/test-optimization-metaheuristics.R:
# the deterministic parts of Metaheuristics.jl 3.6 (operators with their
# random numbers replayed, sorting, crowding, hypervolume contributions,
# SPEA2, NSGA3 and MOEAD_DE helpers, MCCGA bits and Hooke and Jeeves, ECA,
# CGSA, SA and RDEx helpers) and the strings of the OptimizationMetaheuristics
# keyword errors. Run with
#   julia --project=<env with Metaheuristics 3.6> tools/metaheuristics-probe.jl
import Metaheuristics
using Metaheuristics: ECA, εDE
using Random
const MH = Metaheuristics
pr(name, x) = println(name, " = ", repr(x))

pr("refdirs_3_4", MH.gen_ref_dirs(3, 4))
pr("refdirs_2_5", MH.gen_ref_dirs(2, 5))

pr("vs1", MH.violationsSum([1.0, -2.0, 0.5], [1.0e-5, -0.3]; ε = 1.0e-4))
pr("vs2", MH.violationsSum(Float64[], [0.2, -0.1]))
pr("vs3", MH.violationsSum([0.3], Float64[]))

# SBX and polynomial mutation, with the random numbers they draw replayed
b = MH.BoxConstrainedSpace(lb = [-1.0, -1, -1, -1], ub = [1.0, 1, 1, 1])
v1 = [0.1, 0.5, -0.3, 0.9]; v2 = [0.4, 0.2, 0.8, 0.9]
for (k, (η, p)) in enumerate(((15.0, 1.0), (2.0, 0.9)))
    c = MH.SBX_crossover(v1, v2, b, η, p, MersenneTwister(7 + k))
    rng = MersenneTwister(7 + k)
    pr("sbx$(k)_u", rand(rng, 4)); pr("sbx$(k)_R", rand(rng, 4)); pr("sbx$(k)_swap", rand(rng, Bool, 4))
    pr("sbx$(k)_c1", c[1]); pr("sbx$(k)_c2", c[2])
end
w = MH.polynomial_mutation!(copy(v1), b, 20.0, 0.7, MersenneTwister(11))
rng = MersenneTwister(11)
u = rand(rng, 4)
pr("pm_u", u); pr("pm_R", rand(rng, count(<(0.7), u))); pr("pm_v", w)

# ECA, CGSA, SA, RDEx and MOEAD_DE helpers
pr("fmass", MH.fitnessToMass([3.0, -1.0, 2.5, 0.0]))
pr("fmass0", MH.fitnessToMass([0.0, 0.0]))
U = [MH.create_child([0.1, 0.2], 1.0), MH.create_child([-0.5, 0.3], 2.0), MH.create_child([0.7, -0.9], 0.5)]
pr("center", MH.center(U, [0.3, 1.2, 0.5]))
pr("masscalc", MH.massCalculation([3.0, 1.0, 2.0, 5.0], :minimize))
pr("gconst", MH.Gconstant(5, 50))
pr("chaos", [MH.chaos(i, 7, 50, 3.0) for i in 1:10])
pr("newsol", MH.newSol([0.3, -0.7, 0.0], 10.0^2.5))
pr("wmean", MH.RDExModule.AdaptationUtils.weighted_mean([0.5, 0.9, 0.2, 0.7], [1.0, 3.0, 0.5, 2.0], 3, zeros(4)))
pr("moead_g", MH.g([0.3, 0.7], [0.4, 0.6], [0.1, 0.0]))
mp = MH.MOEAD_DE(MH.gen_ref_dirs(2, 7)).parameters
MH.initialize_closest_weight_vectors!(mp)
pr("moead_B", mp.B)

# MCCGA
pr("bits1", MH.bits(0.1))
pr("bits2", MH.bits([0.1, -2.5]))
pr("floats", MH.floats(MH.bits([0.1, -2.5, 3.75e-3])))
pr("hj", MH.hj(x -> (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2, [-1.2, 1.0]))

# Multiobjective: nondominated sorting, crowding distance, hypervolume,
# SPEA2 and NSGA3
F2 = [[0.1, 0.9], [0.4, 0.5], [0.2, 0.8], [0.6, 0.3], [0.5, 0.6], [0.9, 0.05], [0.3, 0.95], [0.8, 0.2]]
sv2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
mkpop(F, x0 = 0.0) = [MH.create_child([x0 + i, 0.0], (F[i], [0.0], [0.0])) for i in eachindex(F)]
pop = mkpop(F2)
pr("nds2", MH.non_dominated_sort(pop))
F3 = [[0.1, 0.9, 0.5], [0.4, 0.5, 0.6], [0.2, 0.8, 0.1], [0.6, 0.3, 0.4], [0.5, 0.6, 0.2], [0.9, 0.05, 0.7],
    [0.3, 0.7, 0.3], [0.8, 0.2, 0.9]]
pop3 = mkpop(F3)
pr("nds3", MH.non_dominated_sort(pop3))
popc = [MH.create_child([1.0 * i, 0.0], (F2[i], [i % 3 == 0 ? 0.5 : -1.0], [0.0])) for i in eachindex(F2)]
pr("ndsc", MH.non_dominated_sort(popc))
front2 = pop[[1, 2, 3, 4, 6, 8]]
cd = MH.compute_crowding_distance(front2)
pr("crowd2_order", [s.x[1] for s in front2]); pr("crowd2", cd)
front3 = copy(pop3)
cd = MH.compute_crowding_distance(front3)
pr("crowd3_order", [s.x[1] for s in front3]); pr("crowd3", cd)

P2 = [0.1 0.2 0.4 0.6 0.8 0.9; 0.9 0.8 0.5 0.3 0.2 0.05]
pr("hv2", MH.hv_contrib_2d(P2))
P3 = reduce(hcat, F3[[1, 2, 3, 4, 6, 7, 8]])
pr("hv3_ref", [maximum(P3[m, :]) + 0.1 * (maximum(P3[m, :]) - minimum(P3[m, :])) for m in 1:3])
pr("hv3", MH.compute_contribution(mkpop(F3[[1, 2, 3, 4, 6, 7, 8]])))

pr("spea2_fit", MH.compute_fitness(pop))
pr("spea2_trunc", MH.truncation(pop, 3))
pr("spea2_fit3", MH.compute_fitness(pop3))

pn = mkpop(F3)
MH.fast_non_dominated_sort!(pn)
pr("nsga3_order", [s.x[1] for s in pn])
Fn = MH.normalize(pn)
pr("nsga3_norm", Fn)
refs = MH.gen_ref_dirs(3, 3)
nich = zeros(Int, length(pn)); freq = zeros(Int, length(refs)); dist = fill(Inf, length(pn))
MH.associate!(nich, freq, dist, Fn, refs, 4)
pr("nsga3_nich", nich); pr("nsga3_freq", freq); pr("nsga3_dist", dist)

pr("igd", MH.PerformanceIndicators.igd([0.1 0.9; 0.5 0.5; 0.95 0.1], [0.0 1.0; 0.25 0.75; 0.5 0.5; 0.75 0.25; 1.0 0.0]))

# The OptimizationMetaheuristics keyword errors name typeof(opt) and its
# supertype, as printed in a session that loaded the package
pr("type_eca", string(typeof(ECA())))
pr("type_eps", string(typeof(εDE())))
pr("super", string(typeof(ECA()).super))
pr("option_fields", propertynames(MH.Options()))
pr("information_fields", propertynames(MH.Information()))
