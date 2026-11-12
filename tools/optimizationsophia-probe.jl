# Reference results of OptimizationSophia (tests/testthat/test-optimization-polyalgorithms.R).
# Run with an environment holding Optimization and OptimizationSophia (1.4.2):
#   julia --project=<env> tools/optimizationsophia-probe.jl
# With k = 1 the Hessian estimate is never made (the index of the single data
# batch is 1 and 1 % 1 is 0), so the iterates are deterministic: clipped sign
# steps with weight decay. Columns: problem, settings, maxiters, retcode, u,
# objective (hexadecimal).
using Optimization, OptimizationSophia, ADTypes, ForwardDiff, Printf

rosen(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
function rosen_grad!(g, x, p)
    g[1] = -2.0 * (p[1] - x[1]) - 4.0 * p[2] * (x[2] - x[1]^2) * x[1]
    g[2] = 2.0 * p[2] * (x[2] - x[1]^2)
    return nothing
end

h(v) = join([@sprintf("%a", a) for a in v], " ")
# the AD backend binds the parameters (Sophia calls the objective as f(u));
# the gradient is the analytic one
rf = OptimizationFunction(rosen, AutoForwardDiff(); grad = rosen_grad!)
prob = OptimizationProblem(rf, [-1.2, 1.0], [1.0, 100.0])
for (sn, alg) in (("k1", Sophia(k = 1)), ("k1eta", Sophia(η = 0.01, k = 1, λ = 0.0, ρ = 0.5))), n in (1, 2, 5, 50, 0)
    sol = n == 0 ? solve(prob, alg) : solve(prob, alg; maxiters = n)
    # the solution carries no statistics (sol.stats is nothing)
    println("rosen\t", sn, "\t", n, "\t", sol.retcode, "\t", h(sol.u), "\t", @sprintf("%a", sol.objective))
end
