# Reference values of the GlobalDiffEq.jl port (R/globaldiffeq.R,
# tests/testthat/test-globaldiffeq.R).
#   julia --project=tools/julia --startup-file=no tools/globaldiffeq-probe.jl
# Out of place vector problems with polynomial right-hand sides, so libm does
# not enter the comparison.
using GlobalDiffEq, OrdinaryDiffEqTsit5, OrdinaryDiffEqLowOrderRK

lv(u, p, t) = [p[1] * u[1] - p[2] * u[1] * u[2], -p[3] * u[2] + p[4] * u[1] * u[2]]
prob = ODEProblem(lv, [1.0, 1.0], (0.0, 5.0), [1.5, 1.0, 3.0, 1.0])
tq = [0.3, 1.7, 4.1]

show(x) = println(repr(x))

for alg in (GLEE23(), GLEE24(), GLEE35(), MM5GEE())
    name = string(nameof(typeof(alg)))
    sol = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-8)
    println("# ", name, " adaptive")
    println("n = ", length(sol.t), " nf = ", sol.stats.nf, " naccept = ", sol.stats.naccept,
            " nreject = ", sol.stats.nreject)
    show(sol.t[2]); show(sol.u[end]); show(global_error_estimate(sol)[end])
    show(global_error_estimate(sol, 3))
    show(reduce(vcat, [sol(t) for t in tq]))
    solf = solve(prob, alg; dt = 0.125, adaptive = false)
    println("# ", name, " fixed dt = 0.125")
    println("n = ", length(solf.t), " nf = ", solf.stats.nf)
    show(solf.u[end]); show(solf.global_error[end])
    sols = solve(prob, alg; reltol = 1.0e-5, abstol = 1.0e-7, saveat = [0.5, 2.5])
    println("# ", name, " saveat")
    show(sols.t); show(reduce(vcat, sols.u)); show(reduce(vcat, sols.global_error))
end

# scalar decay with a length one vector state
dec(u, p, t) = [-p[1] * u[1] * u[1]]
probd = ODEProblem(dec, [2.0], (0.0, 3.0), [0.7])
for alg in (GLEE24(), MM5GEE())
    sol = solve(probd, alg; reltol = 1.0e-8, abstol = 1.0e-10)
    println("# decay ", nameof(typeof(alg)))
    println("n = ", length(sol.t), " nf = ", sol.stats.nf, " nreject = ", sol.stats.nreject)
    show(sol.u[end]); show(sol.global_error[end])
end

# GlobalRichardson: the last inner solve
for (inner, dt) in ((Tsit5(), 0.5), (Euler(), 0.25), (RK4(), 0.25))
    sol = solve(prob, GlobalRichardson(inner); dt = dt, reltol = 1.0e-6, abstol = 1.0e-8)
    println("# GlobalRichardson ", nameof(typeof(inner)), " dt = ", dt)
    println("n = ", length(sol.t), " nf = ", sol.stats.nf)
    show(sol.t[2]); show(sol.u[end])
end
sol = solve(prob, GlobalRichardson(Tsit5()); dt = 0.5)
println("# GlobalRichardson Tsit5 default tolerances")
println("n = ", length(sol.t)); show(sol.u[end])
