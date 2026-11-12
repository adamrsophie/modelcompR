# Reference results of OptimizationODE (tests/testthat/test-optimization-ode.R).
# Run with an environment holding Optimization, OptimizationODE (0.1.12) and
# OrdinaryDiffEqBDF (for the last section):
#   julia --project=<env> tools/optimizationode-probe.jl
# Columns: problem, method, maxiters (0: the default), retcode, u, objective
# (hexadecimal), iterations, fevals, gevals and the final time (hexadecimal).
using Optimization, OptimizationODE, Printf

quad(x, p) = (x[1] - p[1])^2 + (x[2] - p[2])^2
function quad_grad!(g, x, p)
    g[1] = 2.0 * (x[1] - p[1])
    g[2] = 2.0 * (x[2] - p[2])
    return nothing
end
rosen(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
function rosen_grad!(g, x, p)
    g[1] = -2.0 * (p[1] - x[1]) - 4.0 * p[2] * (x[2] - x[1]^2) * x[1]
    g[2] = 2.0 * p[2] * (x[2] - x[1]^2)
    return nothing
end
cobj(x, p) = x[1]^2 + x[2]^2
function cobj_grad!(g, x, p)
    g[1] = 2.0 * x[1]
    g[2] = 2.0 * x[2]
    return nothing
end
ccons(res, x, p) = (res[1] = x[1] + x[2] - 1.0; nothing)
function ccons_j!(J, x, p)
    J[1, 1] = 1.0
    J[1, 2] = 1.0
    return nothing
end

h(v) = join([@sprintf("%a", a) for a in v], " ")
function show_sol(pn, mn, n, sol)
    st = sol.stats
    println(pn, "\t", mn, "\t", n, "\t", sol.retcode, "\t", h(sol.u), "\t", @sprintf("%a", sol.objective), "\t",
            st.iterations, "\t", st.fevals, "\t", st.gevals, "\t", @sprintf("%a", st.time))
end

probs = [("quad", OptimizationProblem(OptimizationFunction(quad; grad = quad_grad!), [2.0, 2.0], [1.0, 1.0])),
         ("rosen", OptimizationProblem(OptimizationFunction(rosen; grad = rosen_grad!), [1.5, 2.0], [1.0, 100.0]))]
methods = [("ODEGradientDescent", ODEGradientDescent), ("RKChebyshevDescent", RKChebyshevDescent),
           ("RKAccelerated", RKAccelerated), ("HighOrderDescent", HighOrderDescent)]
for (pn, prob) in probs, (mn, mk) in methods, n in (1, 2, 5, 50, 0)
    pn == "rosen" && mn == "ODEGradientDescent" && n == 0 && continue
    sol = n == 0 ? solve(prob, mk(); dt = 0.001, maxiters = 1_000_000) : solve(prob, mk(); dt = 0.001, maxiters = n)
    show_sol(pn, mn, n, sol)
end
# the gradient descent on the Rosenbrock function with a larger step
for n in (1, 5, 50, 0)
    sol = n == 0 ? solve(probs[2][2], ODEGradientDescent(); dt = 0.002, maxiters = 1_000_000) :
          solve(probs[2][2], ODEGradientDescent(); dt = 0.002, maxiters = n)
    show_sol("rosen", "ODEGradientDescent2", n, sol)
end
# adaptive methods without an initial dt
for (mn, mk) in methods[2:4]
    show_sol("quad", mn * "nodt", 0, solve(probs[1][2], mk()))
end

# the mass matrix DAE: the last unknown is replaced by the constraint
cf = OptimizationFunction(cobj; grad = cobj_grad!, cons = ccons, cons_j = ccons_j!)
cprob = OptimizationProblem(cf, [1.0, 1.0], [1.0]; lcons = [-10.0], ucons = [10.0])
for n in (5, 50, 0)
    sol = n == 0 ? solve(cprob, DAEMassMatrix(); dt = 0.01, maxiters = 1_000_000) :
          solve(cprob, DAEMassMatrix(); dt = 0.01, maxiters = n)
    println("cons\tDAEMassMatrix\t", n, "\t", sol.retcode, "\t", h(sol.u), "\t", @sprintf("%a", sol.objective))
end

# DAEOptimizer with a fully implicit DAE solver (solve_dae_implicit). Its
# residual du - grad f climbs the objective (DFBDF ends near -4.86e8, 4.86e8
# with the objective 4.7e17); R uses the residual of the gradient flow,
# du + grad f, as the mass matrix form does. The rows below integrate that
# corrected DAE directly with the settings of solve_dae_implicit. They need
# OrdinaryDiffEqBDF in the environment as well.
using OrdinaryDiffEqBDF
println("# upstream DFBDF: ", solve(cprob, DAEOptimizer(DFBDF()); dt = 0.01).u)
function cresid!(res, du, u, p, t)
    cobj_grad!(res, u, p)
    @. res = du + res
    ccons(view(res, 2:2), u, p)
    nothing
end
dprob = OptimizationODE.SciMLBase.DAEProblem(cresid!, zeros(2), [1.0, 1.0], (0.0, 10.0), [1.0])
for (an, alg) in (("DImplicitEuler", DImplicitEuler()), ("DFBDF", DFBDF())), n in (5, 50, 0)
    kw = n == 0 ? (;) : (; maxiters = n)
    sol = solve(dprob, alg; dt = 0.01, initializealg = OptimizationODE.DiffEqBase.ShampineCollocationInit(), kw...)
    println("consimpl\t", an, "\t", n, "\t", sol.retcode, "\t", h(sol.u[end]))
end
