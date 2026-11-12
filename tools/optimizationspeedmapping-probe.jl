# Reference results of OptimizationSpeedMapping and SpeedMapping's ACX method
# (tests/testthat/test-optimization-speedmapping.R). Run with an environment
# holding Optimization, OptimizationSpeedMapping (0.2.9) and SpeedMapping
# (0.4.1):
#   julia --project=<env> tools/optimizationspeedmapping-probe.jl
# Columns: problem, maxiters (0: the default), retcode, u, objective
# (hexadecimal), maps, f_calls, iterations, status, residual norm and last
# learning rate (hexadecimal).
using Optimization, OptimizationSpeedMapping, Printf

rosen(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
function rosen_grad!(g, x, p)
    g[1] = -2.0 * (p[1] - x[1]) - 4.0 * p[2] * (x[2] - x[1]^2) * x[1]
    g[2] = 2.0 * p[2] * (x[2] - x[1]^2)
    return nothing
end
quad3(x, p) = x[1]^2 + 10.0 * (x[2] - 1.0)^2 + 100.0 * (x[3] + 0.5)^2
function quad3_grad!(g, x, p)
    g[1] = 2.0 * x[1]
    g[2] = 20.0 * (x[2] - 1.0)
    g[3] = 200.0 * (x[3] + 0.5)
    return nothing
end

h(v) = join([@sprintf("%a", a) for a in v], " ")
rf = OptimizationFunction(rosen; grad = rosen_grad!)
qf = OptimizationFunction(quad3; grad = quad3_grad!)
probs = [("rosen", OptimizationProblem(rf, [-1.2, 1.0], [1.0, 100.0])),
         ("rosenbox", OptimizationProblem(rf, [0.1, 0.2], [1.0, 100.0]; lb = [-0.5, -0.5], ub = [0.5, 0.8])),
         ("rosenlow", OptimizationProblem(rf, [-1.2, 1.0], [1.0, 100.0]; lb = [-2.0, 1.5], ub = [2.0, 3.0])),
         ("quad3", OptimizationProblem(qf, [3.0, -2.0, 1.0], nothing)),
         ("quad3opt", OptimizationProblem(qf, [0.0, 1.0, -0.5], nothing))]
for (pn, prob) in probs, n in (5, 10, 20, 50, 0)
    sol = n == 0 ? solve(prob, SpeedMappingOpt()) : solve(prob, SpeedMappingOpt(); maxiters = n)
    o = sol.original
    println(pn, "\t", n, "\t", sol.retcode, "\t", h(sol.u), "\t", @sprintf("%a", sol.objective), "\t",
            o.maps, "\t", o.f_calls, "\t", o.iterations, "\t", o.status, "\t", @sprintf("%a", o.residual_norm),
            "\t", @sprintf("%a", o.last_learning_rate))
end
