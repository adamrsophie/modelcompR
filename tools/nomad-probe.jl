# Reference values for NOMADOpt() (R/optimization-nomad.R): runs NOMAD.jl
# as OptimizationNOMAD.jl calls it on polynomial problems and prints R
# assignments. NOMAD_jll has no library for macOS on arm64 (NOMAD.solve
# fails with "libnomadInterface not defined"), so the tests compare with the
# exact minimizers; run this on a supported platform for NOMAD's values.
# Usage:
#   julia --project=<env with NOMAD> tools/nomad-probe.jl > refs.R
using NOMAD

rv(v) = "c(" * join(repr.(Float64.(v)), ", ") * ")"

function run(name, n, f, x0; cons = nothing, lcons = nothing, ucons = nothing, lb = nothing, ub = nothing,
             method = "EB")
    bounds = (;)
    lb === nothing || (bounds = (; bounds..., lower_bound = lb))
    ub === nothing || (bounds = (; bounds..., upper_bound = ub))
    if cons === nothing
        bb = x -> (true, true, [f(x)])
        p = NomadProblem(n, 1, ["OBJ"], bb; bounds...)
    else
        bb = function (x)
            c = cons(x)
            c = max.(lcons .- c, c .- ucons)
            (true, true, vcat(f(x), c))
        end
        p = NomadProblem(n, 1 + length(ucons), vcat("OBJ", fill(method, length(ucons))), bb; bounds...)
    end
    p.options.display_degree = 0
    r = NOMAD.solve(p, x0)
    println("ref\$", name, " <- list(x = ", rv(r.x_best_feas), ", f = ", repr(first(r.bbo_best_feas)), ")")
end

quad(x) = (x[1] - 1) * (x[1] - 1) + 2 * (x[2] + 0.5) * (x[2] + 0.5) + x[1] * x[2]
run("quad", 2, quad, [3.0, 2.0])
run("quad_b", 2, quad, [3.0, 2.0]; lb = [1.5, -2.0], ub = [4.0, 4.0])
run("quad_c", 2, quad, [0.0, 0.0]; cons = x -> [x[1] + x[2]], lcons = [1.0], ucons = [Inf])
run("quad_pb", 2, quad, [0.0, 0.0]; cons = x -> [x[1] + x[2]], lcons = [1.0], ucons = [Inf], method = "PB")
