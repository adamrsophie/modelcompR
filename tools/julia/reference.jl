# Reference values from ModelingToolkit.jl for the modelcompR test suite.
#
# Run from the repository root:
#   julia --startup-file=no --project=tools/julia tools/julia/reference.jl
# Writes CSV files to tests/testthat/fixtures/julia/. Ordinary R tests read
# these files and never need Julia.

using ModelingToolkit, OrdinaryDiffEq, NonlinearSolve
using ModelingToolkit: t_nounits as t, D_nounits as D
using DelimitedFiles

const OUT = joinpath(@__DIR__, "..", "..", "tests", "testthat", "fixtures", "julia")
mkpath(OUT)

function write_series(name, sol, vars, ts)
    data = hcat(ts, [sol(ts; idxs = v) for v in vars]...)
    header = ["t"; [string(v) for v in vars]]
    open(joinpath(OUT, name * ".csv"), "w") do io
        println(io, join(header, ","))
        writedlm(io, data, ',')
    end
end

function write_names(name, sys)
    open(joinpath(OUT, name * "_structure.csv"), "w") do io
        println(io, "kind,name")
        for u in unknowns(sys)
            println(io, "unknown,", string(u))
        end
        for eq in observed(sys)
            println(io, "observed,", string(eq.lhs))
        end
        for p in parameters(sys)
            println(io, "parameter,", string(p))
        end
    end
end

# 1. First order lag
@variables x(t)
@parameters τ
@mtkcompile lag = System([D(x) ~ (1 - x) / τ], t)
sol = solve(ODEProblem(lag, [x => 0.0, τ => 3.0], (0.0, 10.0)), Tsit5(); abstol = 1e-10, reltol = 1e-10)
write_series("lag", sol, [x], collect(0.0:1.0:10.0))

# 2. Second order oscillator with an observed variable
@variables y(t) v(t)
@parameters k = 2.0
@mtkcompile osc = System([D(D(y)) ~ -k * y, v ~ 2 * D(y)], t)
write_names("osc", osc)
sol = solve(ODEProblem(osc, [y => 1.0, D(y) => 0.0], (0.0, 5.0)), Tsit5(); abstol = 1e-10, reltol = 1e-10)
write_series("osc", sol, [y, v], collect(0.0:0.5:5.0))

# 3. Lotka Volterra
@variables prey(t) pred(t)
@parameters α = 1.5 β = 1.0 γ = 3.0 δ = 1.0
@mtkcompile lv = System([D(prey) ~ α * prey - β * prey * pred, D(pred) ~ -γ * pred + δ * prey * pred], t)
sol = solve(ODEProblem(lv, [prey => 1.0, pred => 1.0], (0.0, 10.0)), Tsit5(); abstol = 1e-10, reltol = 1e-10)
write_series("lotka_volterra", sol, [prey, pred], collect(0.0:0.5:10.0))

# 4. Robertson kinetics (stiff)
@variables y1(t) y2(t) y3(t)
@parameters k1 = 0.04 k2 = 3e7 k3 = 1e4
@mtkcompile rober = System([D(y1) ~ -k1 * y1 + k3 * y2 * y3,
                            D(y2) ~ k1 * y1 - k2 * y2^2 - k3 * y2 * y3,
                            D(y3) ~ k2 * y2^2], t)
sol = solve(ODEProblem(rober, [y1 => 1.0, y2 => 0.0, y3 => 0.0], (0.0, 1e5)), Rodas5P(); abstol = 1e-10, reltol = 1e-10)
write_series("robertson", sol, [y1, y2, y3], [0.0, 1.0, 10.0, 100.0, 1000.0, 1e4, 1e5])

# 5. Index 1 DAE kept as a mass matrix problem
@variables a(t) b(t)
@parameters c = 2.0
@mtkcompile dae1 = System([D(a) ~ -a + b, 0 ~ b^3 + b - c * a], t)
write_names("dae1", dae1)
sol = solve(ODEProblem(dae1, [a => 1.0], (0.0, 3.0); guesses = [b => 0.5]), Rodas5P(); abstol = 1e-10, reltol = 1e-10)
write_series("dae1", sol, [a, b], collect(0.0:0.5:3.0))

# 6. Cartesian pendulum (index 3)
@variables px(t) py(t) λ(t)
@parameters g = 9.81 L = 1.0
peqs = [D(D(px)) ~ λ * px, D(D(py)) ~ λ * py - g, px^2 + py^2 ~ L^2]
@mtkcompile pend = System(peqs, t)
write_names("pendulum", pend)
prob = ODEProblem(pend, [px => 1.0, D(px) => 0.0], (0.0, 2.0); guesses = [λ => 0.0, py => 0.0])
sol = solve(prob, Rodas5P(); abstol = 1e-10, reltol = 1e-10)
write_series("pendulum", sol, [px, py], collect(0.0:0.25:2.0))

# 7. RC circuit from the acausal components tutorial, in the functional
# component style of ModelingToolkit 11
@connector function Pin(; name)
    vars = @variables v(t) i(t) [connect = Flow]
    System(Equation[], t, vars, []; name)
end
@component function Ground(; name)
    @named g = Pin()
    System([g.v ~ 0], t, [], []; systems = [g], name)
end
@component function OnePort(; name)
    @named p = Pin()
    @named n = Pin()
    vars = @variables v(t) i(t)
    System([v ~ p.v - n.v, 0 ~ p.i + n.i, i ~ p.i], t, vars, []; systems = [p, n], name)
end
@component function Resistor(; name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R = R
    extend(System([v ~ i * R], t, [], ps; name), oneport)
end
@component function Capacitor(; name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C = C
    extend(System([D(v) ~ i / C], t, [], ps; name), oneport)
end
@component function ConstantVoltage(; name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V = V
    extend(System([V ~ v], t, [], ps; name), oneport)
end
@named resistor = Resistor(R = 1.0)
@named capacitor = Capacitor(C = 1.0)
@named source = ConstantVoltage(V = 1.0)
@named ground = Ground()
rc_eqs = [connect(source.p, resistor.p), connect(resistor.n, capacitor.p),
          connect(capacitor.n, source.n), connect(capacitor.n, ground.g)]
@named rc_model = System(rc_eqs, t; systems = [resistor, capacitor, source, ground])
rc = mtkcompile(rc_model)
write_names("rc", rc)
sol = solve(ODEProblem(rc, [rc.capacitor.v => 0.0], (0.0, 10.0)), Tsit5(); abstol = 1e-10, reltol = 1e-10)
write_series("rc", sol, [rc.capacitor.v, rc.resistor.i], collect(0.0:1.0:10.0))

# 8. Nonlinear system
@variables nx ny nz
@parameters σ = 10.0 ρ = 26.0 βn = 8 / 3
@mtkcompile ns = System([0 ~ σ * (ny - nx), 0 ~ nx * (ρ - nz) - ny, 0 ~ nx * ny - βn * nz])
write_names("nonlinear", ns)
nsol = solve(NonlinearProblem(ns, [nx => 1.0, ny => 0.0, nz => 0.0]), NewtonRaphson())
open(joinpath(OUT, "nonlinear.csv"), "w") do io
    println(io, "name,value")
    for v in [nx, ny, nz]
        println(io, string(v), ",", nsol[v])
    end
end

println("reference values written to ", OUT)
