# Deterministic pieces of CMAEvolutionStrategy.jl (0.2.6) and GCMAES.jl
# (0.1.34) used by tests/testthat/test-optimization-cma.R: the strategy
# parameters, the recombination weights and the bound transformation.
using CMAEvolutionStrategy, GCMAES

const C = CMAEvolutionStrategy
for (n, popsize) in ((2, 6), (10, 10), (2, 20))
    p = C.Parameters(zeros(n), 0.1; lower = fill(-1.0, n), upper = fill(0.8, n), popsize = popsize)
    println("cma n=$n lambda=$(p.λ)")
    println("  weights = ", repr(p.weights.weights))
    println("  mueff = ", repr(p.weights.μeff))
    println("  sigma c d e chin = ", repr([p.sigma.c, p.sigma.d, p.sigma.e, p.sigma.χn]))
    println("  cov c c1 cmu e = ", repr([p.cov.c, p.cov.c1, p.cov.cμ, p.cov.e]))
    println("  mean = ", repr(p.mean))
end
b = C.BoxConstraints([-1.0, -1.0], [0.8, 0.8])
for x in ([0.0, 0.0], [-1.03, 0.79], [0.83, -1.2], [5.0, -7.3], [0.76, -0.97])
    println("transform ", repr(x), " = ", repr(C.transform(b, x)))
end
for y in ([0.0, 0.0], [-0.99, 0.79], [0.8, -1.0])
    println("backtransform ", repr(y), " = ", repr(C.backtransform(b, y)))
end
b2 = C.BoxConstraints([-10.0, 2.0], [Inf, 3.0])
println("transform2 = ", repr(C.transform(b2, [-12.0, 3.4])))

for (n, s0) in ((2, 0.2), (10, 0.2), (1, 0.2))
    o = GCMAES.CMAESOpt(x -> sum(abs2, x), x -> 2x, zeros(n), s0, fill(-10.0, n), fill(10.0, n))
    println("gcmaes n=$n lambda=$(o.λ) mu=$(o.μ)")
    println("  w = ", repr(o.w))
    println("  mueff cc cs c1 cmu ds chin = ", repr([o.μeff, o.cc, o.cσ, o.c1, o.cμ, o.dσ, o.χₙ]))
end
