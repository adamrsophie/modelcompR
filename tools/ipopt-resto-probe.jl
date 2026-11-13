# Ipopt iteration logs for the restoration phase of the interior point
# solver (R/dynamic-opt-nlp.R): the cart pole of ModelingToolkit's tests
# solved from rest without an initial trajectory, with the JuMP backend and
# Ipopt at print_level 5, plus small NLPs whose solves enter the
# restoration phase. Usage (packages: JuMP Ipopt InfiniteOpt):
#   julia --project=tools/julia tools/ipopt-resto-probe.jl > ipopt-resto.txt
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
import InfiniteOpt
using JuMP, Ipopt
import OrdinaryDiffEqExplicitTableaus as ExplicitTableaus
const Symbolics = ModelingToolkit.Symbolics

println("== versions")
for p in ("Ipopt", "Ipopt_jll", "JuMP", "InfiniteOpt", "ModelingToolkit")
    try
        println(p, " ", pkgversion(Base.require(Main, Symbol(p))))
    catch
        println(p, " ?")
    end
end

@parameters g l mₚ mₖ
@variables x(..) θ(..) u(t) [input = true, bounds = (-10, 10)]
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
costs = [Symbolics.Integral(t in (0, tf))(u^2)]
@named cartpole = System(eqs, t; costs, constraints = cons)
cartpole = mtkcompile(cartpole; inputs = [u])
var_order = [θ(t), x(t), D(θ(t)), D(x(t))]
u0map = [D(x(t)) => 0.0, D(θ(t)) => 0.0, θ(t) => 0.0, x(t) => 0.0]
pmap = [mₖ => 1.0, mₚ => 0.2, l => 0.5, g => 9.81, u => 0]

for (name, dt, tab) in (("rk4_dt0.04", 0.04, ExplicitTableaus.RK4()),
                        ("rk4_dt0.25", 0.25, ExplicitTableaus.RK4()))
    println("== cartpole ", name)
    jprob = JuMPDynamicOptProblem(cartpole, [u0map; pmap], (0, tf); dt = dt)
    flush(stdout)
    jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, tab); verbose = true)
    flush(stdout)
    println("retcode ", jsol.sol.retcode)
    println("final ", repr(Float64.(jsol.sol[var_order][end])))
    println("u ", repr(Float64.(jsol.input_sol.u isa AbstractVector ? reduce(vcat, jsol.input_sol.u) : jsol.input_sol.u)))
end

# Small NLPs that need the restoration phase
function small(name, build)
    println("== ", name)
    m = Model(Ipopt.Optimizer)
    set_attribute(m, "print_level", 5)
    xs = build(m)
    optimize!(m)
    println("status ", termination_status(m), " ", raw_status(m))
    println("x ", repr(value.(xs)))
    println("objective ", repr(objective_value(m)))
    println("iterations ", barrier_iterations(m))
end
# infeasible start for x^2 + y^2 = 1, y = x^2 + 2 (locally infeasible)
small("infeasible_circle", m -> begin
    @variable(m, x, start = 0.0); @variable(m, y, start = 0.0)
    @objective(m, Min, x^2 + y^2)
    @constraint(m, x^2 + y^2 == 1); @constraint(m, y == x^2 + 2)
    [x, y]
end)
# feasible but a bad start: Waechter and Biegler's example (2006), x1 = 1 + ...
small("wachter_biegler", m -> begin
    @variable(m, x1, start = -2.0); @variable(m, x2 >= 0, start = 1.0); @variable(m, x3 >= 0, start = 1.0)
    @objective(m, Min, x1)
    @constraint(m, x1^2 - x2 - 1 == 0); @constraint(m, x1 - x3 - 0.5 == 0)
    [x1, x2, x3]
end)
# HS071 from a start far outside
small("hs071_far", m -> begin
    @variable(m, 1 <= x[1:4] <= 5)
    start = [5.0, 1.0, 1.0, 5.0]
    for i in 1:4; set_start_value(x[i], start[i]); end
    @objective(m, Min, x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3])
    @constraint(m, x[1] * x[2] * x[3] * x[4] >= 25)
    @constraint(m, sum(x .^ 2) == 40)
    x
end)
