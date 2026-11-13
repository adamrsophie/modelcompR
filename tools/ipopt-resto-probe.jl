# Ipopt iteration logs for the restoration phase of the interior point
# solver (R/dynamic-opt-nlp.R): small NLPs whose solves enter the
# restoration phase, and the cart pole of ModelingToolkit's tests solved
# from rest without an initial trajectory, with the JuMP backend and Ipopt
# at print_level 5. Usage (packages: JuMP Ipopt InfiniteOpt):
#   julia --project=tools/julia tools/ipopt-resto-probe.jl > ipopt-resto.txt
import JuMP
import Ipopt

println("== versions")
for p in ("Ipopt", "JuMP")
    try
        println(p, " ", pkgversion(Base.require(Main, Symbol(p))))
    catch
        println(p, " ?")
    end
end
flush(stdout)

# Small NLPs that need the restoration phase
function small(name, build)
    println("== ", name)
    flush(stdout)
    try
        m = JuMP.Model(Ipopt.Optimizer)
        JuMP.set_attribute(m, "print_level", 5)
        xs = build(m)
        JuMP.optimize!(m)
        flush(stdout)
        println("status ", JuMP.termination_status(m), " ", JuMP.raw_status(m))
        println("x ", repr(JuMP.value.(xs)))
        println("objective ", repr(JuMP.objective_value(m)))
        println("iterations ", JuMP.barrier_iterations(m))
    catch err
        println("error ", sprint(showerror, err))
    end
    flush(stdout)
end
# infeasible start for x^2 + y^2 = 1, y = x^2 + 2 (locally infeasible)
small("infeasible_circle", m -> begin
    JuMP.@variable(m, x, start = 0.0); JuMP.@variable(m, y, start = 0.0)
    JuMP.@objective(m, Min, x^2 + y^2)
    JuMP.@constraint(m, x^2 + y^2 == 1); JuMP.@constraint(m, y == x^2 + 2)
    [x, y]
end)
# Waechter and Biegler's example (2006) from x1 = -2
small("wachter_biegler", m -> begin
    JuMP.@variable(m, x1, start = -2.0); JuMP.@variable(m, x2 >= 0, start = 1.0)
    JuMP.@variable(m, x3 >= 0, start = 1.0)
    JuMP.@objective(m, Min, x1)
    JuMP.@constraint(m, x1^2 - x2 - 1 == 0); JuMP.@constraint(m, x1 - x3 - 0.5 == 0)
    [x1, x2, x3]
end)
# HS071 from a start far outside
small("hs071_far", m -> begin
    JuMP.@variable(m, 1 <= x[1:4] <= 5)
    start = [5.0, 1.0, 1.0, 5.0]
    for i in 1:4; JuMP.set_start_value(x[i], start[i]); end
    JuMP.@objective(m, Min, x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3])
    JuMP.@constraint(m, x[1] * x[2] * x[3] * x[4] >= 25)
    JuMP.@constraint(m, sum(x .^ 2) == 40)
    x
end)

# The cart pole of ModelingToolkit's tests (dynamic_optimization.jl)
try
    @eval using ModelingToolkit
    @eval import InfiniteOpt
    @eval import OrdinaryDiffEqExplicitTableaus as ExplicitTableaus
    @eval begin
        using ModelingToolkit: t_nounits as t, D_nounits as D
        const MTK = ModelingToolkit
        MTK.@parameters g l mₚ mₖ
        MTK.@variables x(..) θ(..) u(t) [input = true, bounds = (-10, 10)]
        s = sin(θ(t)); c = cos(θ(t))
        H = [mₖ+mₚ mₚ*l*c; mₚ*l*c mₚ*l^2]
        C = [0 -mₚ*D(θ(t))*l*s; 0 0]
        qd = [D(x(t)), D(θ(t))]
        G = [0, mₚ * g * l * s]
        B = [1, 0]
        tf = 5
        rhss = -H \ Vector(C * qd + G - B * u)
        eqs = [D(D(x(t))) ~ rhss[1], D(D(θ(t))) ~ rhss[2]]
        cons = [θ(tf) ~ π, x(tf) ~ 0, D(θ(tf)) ~ 0, D(x(tf)) ~ 0]
        costs = [MTK.Symbolics.Integral(t in (0, tf))(u^2)]
        MTK.@named cartpole = MTK.System(eqs, t; costs, constraints = cons)
        cartpole = MTK.mtkcompile(cartpole; inputs = [u])
        var_order = [θ(t), x(t), D(θ(t)), D(x(t))]
        u0map = [D(x(t)) => 0.0, D(θ(t)) => 0.0, θ(t) => 0.0, x(t) => 0.0]
        pmap = [mₖ => 1.0, mₚ => 0.2, l => 0.5, g => 9.81, u => 0]
        for (name, dt) in (("rk4_dt0.04", 0.04), ("rk4_dt0.25", 0.25))
            println("== cartpole ", name)
            flush(stdout)
            try
                jprob = MTK.JuMPDynamicOptProblem(cartpole, [u0map; pmap], (0, tf); dt = dt)
                jsol = MTK.solve(jprob, MTK.JuMPCollocation(Ipopt.Optimizer, ExplicitTableaus.RK4()); verbose = true)
                flush(stdout)
                println("retcode ", jsol.sol.retcode)
                println("final ", repr(Float64.(jsol.sol[var_order][end])))
            catch err
                println("error ", sprint(showerror, err))
            end
            flush(stdout)
        end
    end
catch err
    println("error ", sprint(showerror, err))
end
