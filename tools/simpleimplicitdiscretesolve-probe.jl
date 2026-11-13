# Reference values for SimpleIDSolve (SimpleImplicitDiscreteSolve 2.0.4,
# SimpleNonlinearSolve's SimpleNewtonRaphson with a ForwardDiff Jacobian).
# Run with: julia --project=tools/julia tools/simpleimplicitdiscretesolve-probe.jl
using SimpleImplicitDiscreteSolve

show_sol(name, sol) = begin
    println("# ", name)
    println("retcode = ", sol.retcode)
    println("t = ", repr(collect(sol.t)))
    println("u = ", repr(reduce(vcat, sol.u)))
end

# u_next solves u_next + u_next^3 = u
cubic(un, u, p, t) = [un[1] + un[1] * un[1] * un[1] - u[1]]
prob = ImplicitDiscreteProblem(cubic, [2.0], (0.0, 3.0))
show_sol("cubic dt 1", solve(prob, SimpleIDSolve()))

# implicit Euler of Lotka-Volterra with a time dependent forcing, dt = 0.1
function lv(un, u, p, t)
    h = 0.1
    return [
        un[1] - u[1] - h * (p[1] * un[1] - p[2] * un[1] * un[2]),
        un[2] - u[2] - h * (-p[3] * un[2] + p[4] * un[1] * un[2] + t * t),
    ]
end
prob = ImplicitDiscreteProblem(lv, [1.0, 1.0], (0.0, 1.0), [1.5, 1.0, 3.0, 1.0])
show_sol("lv dt 0.1", solve(prob, SimpleIDSolve(); dt = 0.1))
show_sol("lv dt 0.1 save_everystep false", solve(prob, SimpleIDSolve(); dt = 0.1, save_everystep = false))
show_sol("lv dt 0.1 save_everystep false save_start false",
    solve(prob, SimpleIDSolve(); dt = 0.1, save_everystep = false, save_start = false))
try
    solve(prob, SimpleIDSolve(); dt = 0.1, save_start = false)
    println("# save_start false: no error")
catch e
    println("# save_start false: ", nameof(typeof(e)))
end

# the initial solve f(u, u0, p, t0) = 0 replaces u0
alg2(un, u, p, t) = [un[1] - 0.5 * u[1], un[2] - un[1] * un[1]]
prob = ImplicitDiscreteProblem(alg2, [1.0, 5.0], (0.0, 2.0))
show_sol("init replaces u0", solve(prob, SimpleIDSolve()))

# no real root at the initial time
noroot(un, u, p, t) = [un[1] * un[1] + 1.0]
prob = ImplicitDiscreteProblem(noroot, [1.3], (0.0, 2.0))
show_sol("initial failure", solve(prob, SimpleIDSolve()))

# a root until u drops below 1
shrink(un, u, p, t) = [un[1] * un[1] - (u[1] - 1.0)]
prob = ImplicitDiscreteProblem(shrink, [5.0], (0.0, 5.0))
show_sol("convergence failure", solve(prob, SimpleIDSolve()))
