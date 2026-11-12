# Reference values for the BBO_borg_moea port (R/optimization-bbo.R).
# BlackBoxOptim 0.6.13, Julia 1.13. Run with
#   julia --project=<env with BlackBoxOptim> tools/bbo-probe.jl
# The run of Borg MOEA uses Julia's random numbers, which R cannot draw;
# the deterministic parts (epsilon box indexing, comparisons, the archive)
# are printed for exact comparison, and each genetic operator is printed
# with the random numbers it drew (replayed by the R tests).
using BlackBoxOptim, Random, LinearAlgebra, Printf
const B = BlackBoxOptim

# BBO_HEX=1 prints hexadecimal floats, which R reads exactly
f17(x) = get(ENV, "BBO_HEX", "") == "1" ? @sprintf("%a", x) : @sprintf("%.17g", x)
vec17(v) = join(map(f17, v), " ")

println("## epsilon box index")
for (u, eps, minimizing) in [
        ((0.123, 4.56), [0.1, 0.1], true),
        ((-0.35, 2.0), [0.1, 0.1], true),
        ((0.3, 0.7), [0.1, 0.1], true),
        ((NaN, 1.5), [0.1, 0.1], true),
        ((0.123, -4.56), [0.1, 0.5], false),
        ((1.0e-17, 3.0), [0.25, 0.1], true),
    ]
    f = B.IndexedTupleFitness(u, sum(u), eps, Val{minimizing})
    println("IDX ", vec17(u), " | ", vec17(eps), " | ", minimizing, " | ", f17(f.agg), " ",
        join(string.(f.index), " "), " ", f17(f.dist))
end

println("## hat compare")
fs = B.EpsBoxDominanceFitnessScheme{2, Float64}(0.1)
mk(u) = convert(B.IndexedTupleFitness{2, Float64}, u, fs)
for (u, v) in [((0.12, 0.33), (0.25, 0.41)), ((0.12, 0.45), (0.25, 0.41)), ((0.12, 0.33), (0.13, 0.34)),
               ((0.11, 0.31), (0.19, 0.39)), ((0.15, 0.35), (0.15, 0.35)), ((0.5, 0.1), (0.1, 0.5))]
    println("HAT ", vec17(u), " | ", vec17(v), " | ", join(B.hat_compare(mk(u), mk(v), fs), " "))
end

println("## archive")
arch = B.EpsBoxArchive(fs, max_size = 100)
x0 = [0.0, 0.0]
for k in 1:60
    u = (abs(sin(1.7 * k)) * 2, abs(cos(0.9 * k)) * (1 + k / 30))
    x = [k / 10, -k / 7]
    B.add_candidate!(arch, mk(u), x, mod(k, 7) + 1, k)
    if k == 30
        B.notify!(arch, :restart)
    end
    b = arch.best_front_elem
    println("ARCH ", vec17(u), " | ", k, " ", length(arch), " ", arch.num_candidates, " ", arch.last_progress, " ",
        B.noprogress_streak(arch, since_restart = true), " ", vec17(b.params), " ", f17(b.fitness.agg))
end
for el in B.pareto_frontier(arch)
    println("FRONT ", vec17(el.fitness.orig), " ", el.tag, " ", el.num_fevals, " ", el.n_restarts, " ", vec17(el.params))
end
tc = B.tagcounts(arch, 0.9)
println("TAGS ", join([string(k, "=", f17(tc[k])) for k in sort(collect(keys(tc)))], " "))

println("## operators")
pop = [0.1 0.7 -0.3 0.45 0.9 -0.8;
       0.5 -0.2 0.6 0.05 -0.4 0.3;
       -0.6 0.2 0.8 -0.1 0.35 0.0;
       0.25 0.9 -0.7 0.6 -0.15 0.55]
ss = B.RectSearchSpace([(-1.0, 1.0) for _ in 1:4])
for seed in 1:3
    # DiffEvoRandBin1
    Random.seed!(seed)
    t = fill(NaN, 4)
    B.apply!(B.DiffEvoRandBin1(B.DE_DefaultOptions), t, 0, pop, [2, 5, 3])
    Random.seed!(seed)
    m = rand(1:4)
    d = Float64[m]
    for i in 1:4
        i == m || push!(d, rand())
    end
    println("OP de ", seed, " | ", vec17(t), " | ", vec17(d))
    # SimulatedBinaryCrossover
    Random.seed!(seed)
    ts = [fill(NaN, 4), fill(NaN, 4)]
    B.apply!(B.SimulatedBinaryCrossover(B.SBX_DefaultOptions), ts, [0, 0], pop, [4, 1])
    Random.seed!(seed)
    d = Float64[]
    for i in 1:4
        r = rand()
        push!(d, r)
        if r <= 0.2
            push!(d, rand())
            push!(d, rand(Bool))
        end
    end
    println("OP sbx ", seed, " | ", vec17(vcat(ts...)), " | ", vec17(d))
    # SimplexCrossover{3}
    Random.seed!(seed)
    t = fill(NaN, 4)
    B.apply!(B.SimplexCrossover{3}(B.SPX_DefaultOptions), t, 0, pop, [6, 2, 4])
    Random.seed!(seed)
    println("OP spx ", seed, " | ", vec17(t), " | ", vec17([rand(), rand()]))
    # ParentCentricCrossover{2}, {3}
    for (np, pix) in [(2, [3, 5]), (3, [1, 4, 6])]
        Random.seed!(seed)
        t = fill(NaN, 4)
        B.apply!(B.ParentCentricCrossover{np}(B.PCX_DefaultOptions), t, 0, pop, pix)
        Random.seed!(seed)
        d = vcat(randn(np - 1), randn())
        println("OP pcx", np, " ", seed, " | ", vec17(t), " | ", vec17(d))
    end
    # UnimodalNormalDistributionCrossover{2}, {3}
    for (np, pix) in [(2, [2, 6]), (3, [5, 3, 1])]
        Random.seed!(seed)
        t = fill(NaN, 4)
        B.apply!(B.UnimodalNormalDistributionCrossover{np}(B.UNDX_DefaultOptions), t, 0, pop, pix)
        Random.seed!(seed)
        d = [randn() for _ in 1:min(4, np)]
        println("OP undx", np, " ", seed, " | ", vec17(t), " | ", vec17(d))
    end
    # MutationClock over PolynomialMutation and UniformMutation
    for (nm, inner) in [("pm", B.PolynomialMutation(ss, B.PM_DefaultOptions)), ("um", B.UniformMutation(ss))]
        Random.seed!(seed)
        mc = B.MutationClock(inner, 1 / 4)
        t = copy(pop[:, 2])
        B.apply!(mc, t, 1)
        B.apply!(mc, t, 1)
        Random.seed!(seed)
        d = Float64[]
        nv() = (r = rand(); push!(d, r); ceil(Int, -log(r) * 4))
        i = 1 + nv()
        for rep in 1:2
            while i <= 4
                push!(d, rand())
                i += nv()
            end
            i -= 4
        end
        println("OP clock", nm, " ", seed, " | ", vec17(t), " ", mc.nextVarToMutate, " | ", vec17(d))
    end
    # RandomBound
    Random.seed!(seed)
    t = [1.3, -0.2, -1.7, 0.95]
    B.apply!(B.RandomBound(ss), t, pop, 2)
    Random.seed!(seed)
    println("OP bound ", seed, " | ", vec17(t), " | ", vec17([rand(), rand()]))
end

println("## runs")
sphere_rastrigin(x) = (sum(x .^ 2), sum(x .^ 2 .- 10 .* cos.(2pi .* x) .+ 10))
rosen_ackley(x) = ((1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2,
    -20.0 * exp(-0.2 * sqrt(0.5 * (x[1]^2 + x[2]^2))) - exp(0.5 * (cos(2pi * x[1]) + cos(2pi * x[2]))) + exp(1) + 20.0)
function zdt1(x)
    f1 = x[1]
    g = 1 + 9 * sum(x[2:end]) / (length(x) - 1)
    return (f1, g * (1 - sqrt(f1 / g)))
end
# the arguments of OptimizationBBO's __solve for BBO_borg_moea()
function run(f, maxsteps; cb = nothing)
    args = Dict{Symbol, Any}(:Method => :borg_moea, :SearchRange => [(0.0, 2.0), (0.0, 2.0)],
        :NumDimensions => 2, :FitnessScheme => ParetoFitnessScheme{2}(is_minimizing = true), :TraceMode => :silent)
    maxsteps === nothing || (args[:MaxSteps] = maxsteps)
    if cb !== nothing
        args[:CallbackFunction] = cb
        args[:CallbackInterval] = 0.0
    end
    setup = bbsetup(f; args...)
    return bboptimize(setup, [0.25, 0.25])
end
for (nm, f) in [("sphere_rastrigin", sphere_rastrigin), ("rosen_ackley", rosen_ackley), ("zdt1", zdt1)]
    for seed in 1:5
        Random.seed!(seed)
        r = run(f, nothing)
        println("RUN ", nm, " ", seed, " | ", vec17(best_fitness(r)), " | ", vec17(best_candidate(r)), " | ",
            r.iterations, " ", r.f_calls, " ", length(pareto_frontier(r)), " | ", r.stop_reason)
    end
end
Random.seed!(1)
r = run(sphere_rastrigin, 100)
println("SHORT ", r.iterations, " ", r.f_calls, " | ", r.stop_reason)
n = Ref(0)
Random.seed!(1)
r = run(sphere_rastrigin, nothing, cb = ctrl -> (n[] += 1; n[] == 25 && B.shutdown_optimizer!(ctrl)))
println("SHUTDOWN ", n[], " ", r.iterations, " ", r.f_calls, " | ", r.stop_reason)
