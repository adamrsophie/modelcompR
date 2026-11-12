# Reference iterates of the Optimisers.jl rules driven by OptimizationOptimisers
# (tests/testthat/test-optimization-optimisers.R). The analytic gradient is given
# to both sides. Run with an environment holding Optimization,
# OptimizationOptimisers and Optimisers (OptimizationOptimisers 0.3.26,
# Optimisers 0.4.11):
#   julia --project=<env> tools/optimisers-probe.jl
using Optimization, OptimizationOptimisers, Optimisers, SciMLBase, ForwardDiff, ADTypes

quad(x, p) = 1.5 * x[1] * x[1] + x[1] * x[2] + x[2] * x[2] - x[1] + x[2]
quadg!(G, x, p) = (G[1] = 3 * x[1] + x[2] - 1; G[2] = x[1] + 2 * x[2] + 1; nothing)
rosen(x, p) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
rosg!(G, x, p) = (G[1] = -2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2); G[2] = 200 * (x[2] - x[1]^2); nothing)

probs = [("quad", OptimizationProblem(OptimizationFunction(quad, AutoForwardDiff(); grad = quadg!), [0.7, -0.4], [1.0])),
         ("rosen", OptimizationProblem(OptimizationFunction(rosen, AutoForwardDiff(); grad = rosg!), [-1.2, 1.0], [1.0]))]

rules = [
    ("Descent", () -> Descent()), ("Descent05", () -> Descent(0.05)),
    ("Momentum", () -> Momentum()), ("Nesterov", () -> Nesterov(0.01)),
    ("RMSProp", () -> RMSProp(0.01)), ("RMSPropc", () -> RMSProp(0.01; centred = true)),
    ("Rprop", () -> Rprop(0.01)), ("Adam", () -> Adam(0.05)), ("Lion", () -> Lion(0.01)),
    ("RAdam", () -> RAdam(0.05)), ("AdaMax", () -> AdaMax(0.05)), ("OAdam", () -> OAdam(0.05)),
    ("AdaGrad", () -> AdaGrad(0.1)), ("AdaDelta", () -> AdaDelta(0.9)), ("AMSGrad", () -> AMSGrad(0.05)),
    ("NAdam", () -> NAdam(0.05)), ("AdamW", () -> AdamW(0.05, (0.9, 0.999), 0.1)),
    ("AdamWd", () -> AdamW(0.05, (0.9, 0.999), 0.1; couple = false)), ("AdaBelief", () -> AdaBelief(0.05)),
    ("ChainWD", () -> OptimiserChain(WeightDecay(0.01), Adam(0.05))),
    ("ChainSD", () -> OptimiserChain(SignDecay(0.01), Descent(0.01))),
    ("ChainCG", () -> OptimiserChain(ClipGrad(0.5), Descent(0.05))),
    ("ChainCN", () -> OptimiserChain(ClipNorm(1.0), Momentum(0.02))),
    ("ChainCN1", () -> OptimiserChain(ClipNorm(1.0, 1), Descent(0.02))),
    ("ChainAG", () -> OptimiserChain(AccumGrad(3), Descent(0.02))),
]

for (pn, prob) in probs, (rn, mk) in rules, n in (1, 2, 3, 7, 40)
    sol = solve(prob, mk(); maxiters = n)
    println(pn, "\t", rn, "\t", n, "\t", sol.retcode, "\t", repr(sol.u[1]), "\t", repr(sol.u[2]), "\t", repr(sol.objective))
end
# save_best = false keeps the last iterate
for (pn, prob) in probs
    sol = solve(prob, Adam(0.05); maxiters = 40, save_best = false)
    println(pn, "\tAdamNoBest\t40\t", sol.retcode, "\t", repr(sol.u[1]), "\t", repr(sol.u[2]), "\t", repr(sol.objective))
end
