# Reference results of SimpleOptimization (tests/testthat/test-optimization-simple.R).
# SimpleOptimization differentiates the objective with ForwardDiff; the R
# tests give the gradient and Hessian formulas that repeat the ForwardDiff
# arithmetic, and the first block below checks that they agree bit for bit.
# Run with an environment holding Optimization, SimpleOptimization (2.0.2),
# SimpleNonlinearSolve, LineSearch and ForwardDiff:
#   julia --project=<env> tools/simpleoptimization-probe.jl
using Optimization, SimpleOptimization, ForwardDiff, ADTypes, Random, Printf, LinearAlgebra

quad(x, p) = 1.5 * x[1] * x[1] + x[1] * x[2] + x[2] * x[2] - x[1] + x[2]
rosen(x, p) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
quadg(x) = [3 * x[1] + x[2] - 1, x[1] + 2 * x[2] + 1]
function roseng(x)
    t2 = x[2] - x[1] * x[1]
    return [-(2 * (1 - x[1])) + 100 * ((2 * t2) * (-(2 * x[1]))), 100 * (2 * t2)]
end
function rosenh(x)
    t2 = x[2] - x[1] * x[1]
    h11 = 2 + 100 * (8 * (x[1] * x[1]) - 4 * t2)
    h12 = 100 * (-(4 * x[1]))
    return [h11 h12; h12 200.0]
end

Random.seed!(11)
bad = zeros(Int, 3)
for k in 1:20000
    x = randn(2) .* 2.0 .^ rand(-3:3)
    bad[1] += ForwardDiff.gradient(z -> quad(z, nothing), x) != quadg(x)
    bad[2] += ForwardDiff.gradient(z -> rosen(z, nothing), x) != roseng(x)
    bad[3] += ForwardDiff.jacobian(y -> ForwardDiff.gradient(z -> rosen(z, nothing), y), x) != rosenh(x)
end
println("# derivative formula mismatches (quad gradient, rosen gradient, rosen Hessian): ", bad)

h(v) = join([@sprintf("%a", a) for a in v], " ")
function show_sol(pn, mn, n, sol)
    println(pn, "\t", mn, "\t", n, "\t", sol.retcode, "\t", h(vec(sol.u)), "\t", @sprintf("%a", sol.objective), "\t",
            sol.stats.iterations, "\t", sol.stats.fevals)
end

qf = OptimizationFunction(quad, AutoForwardDiff())
rf = OptimizationFunction(rosen, AutoForwardDiff())
probs = [("quad", OptimizationProblem(qf, [0.7, -0.4], [1.0])),
         ("rosen", OptimizationProblem(rf, [-1.2, 1.0], [1.0])),
         ("rosenbox", OptimizationProblem(rf, [0.1, 0.2], [1.0]; lb = [-0.5, -0.5], ub = [0.5, 0.8]))]
methods = [("SimpleLBFGS", () -> SimpleLBFGS()),
           ("SimpleLBFGS3", () -> SimpleLBFGS(threshold = 3)),
           ("SimpleBFGS", () -> SimpleBFGS()),
           ("SimpleNewton", () -> SimpleNewton()),
           ("SimpleGradientDescent", () -> SimpleGradientDescent()),
           ("SimpleGradientDescent3", () -> SimpleGradientDescent(eta = 0.001)),
           ("SimpleSOAP", () -> SimpleSOAP())]
for (pn, prob) in probs[1:2], (mn, mk) in methods, n in (1, 2, 3, 5, 10, 25, 0)
    sol = n == 0 ? solve(prob, mk()) : solve(prob, mk(); maxiters = n)
    show_sol(pn, mn, n, sol)
end
for (mn, mk) in methods[1:2], n in (1, 2, 3, 5, 10, 0)
    sol = n == 0 ? solve(probs[3][2], mk()) : solve(probs[3][2], mk(); maxiters = n)
    show_sol("rosenbox", mn, n, sol)
end
for n in (1, 2, 3, 5, 10, 0)
    sol = n == 0 ? solve(probs[1][2], SimpleLBFGS(); abstol = 1e-3) :
          solve(probs[1][2], SimpleLBFGS(); maxiters = n, abstol = 1e-3)
    show_sol("quad", "SimpleLBFGSabstol", n, sol)
end

# the matrix path of SimpleSOAP. With max_dim below both dimensions it keeps
# the identity bases (exact comparison); with one basis (max_dim = 2) and with
# both, the eigen decompositions and QR factorizations of LAPACK enter.
mf(x, p) = sum(abs2, x .- p)
P = reshape(collect(1.0:6.0), 3, 2) ./ 4
mprob = OptimizationProblem(OptimizationFunction(mf, AutoForwardDiff()), zeros(3, 2), P)
for (mn, mk) in (("SOAPid", () -> SimpleSOAP(eta = 0.05, max_dim = 1)),
                 ("SOAPright", () -> SimpleSOAP(eta = 0.05, max_dim = 2, freq = 3)),
                 ("SOAPboth", () -> SimpleSOAP(eta = 0.05, weight_decay = 0.0))), n in (1, 2, 5, 25, 400)
    show_sol("matrix", mn, n, solve(mprob, mk(); maxiters = n))
end
