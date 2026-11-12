# Reference results of OptimizationPolyalgorithms
# (tests/testthat/test-optimization-polyalgorithms.R). Run with an environment
# holding Optimization and OptimizationPolyalgorithms (0.3.12):
#   julia --project=<env> tools/optimizationpolyalgorithms-probe.jl
# Columns: problem, maxiters (0: the default), retcode, u, objective
# (hexadecimal), iterations, fevals and gevals.
using Optimization, OptimizationPolyalgorithms, ADTypes, ForwardDiff, Printf

rosen(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
function rosen_grad!(g, x, p)
    g[1] = -2.0 * (p[1] - x[1]) - 4.0 * p[2] * (x[2] - x[1]^2) * x[1]
    g[2] = 2.0 * p[2] * (x[2] - x[1]^2)
    return nothing
end

h(v) = join([@sprintf("%a", a) for a in v], " ")
# the AD backend binds the parameters (OptimizationOptimisers calls the
# objective as f(u)); the gradient is the analytic one
rf = OptimizationFunction(rosen, AutoForwardDiff(); grad = rosen_grad!)
# PolyOpt does not allow bounds (SciMLBase.allowsbounds is false), so its
# bounded branch is never reached
probs = [("rosen", OptimizationProblem(rf, [-1.2, 1.0], [1.0, 100.0]))]
for (pn, prob) in probs, n in (5, 20, 100, 0)
    sol = n == 0 ? solve(prob, PolyOpt()) : solve(prob, PolyOpt(); maxiters = n)
    st = sol.stats
    println(pn, "\t", n, "\t", sol.retcode, "\t", h(sol.u), "\t", @sprintf("%a", sol.objective), "\t",
            st.iterations, "\t", st.fevals, "\t", st.gevals)
end
