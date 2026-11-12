# Reference values for tests/testthat/test-optimization-quaddirect.R: the
# evaluations of QuadDIRECT.analyze (master of timholy/QuadDIRECT.jl) as
# called by OptimizationQuadDIRECT, in order, and the minimum.
# The test functions use only +, -, * and / (and squares, which Julia and R
# both compute as x * x) so that their values agree in R.
# Run with an environment holding QuadDIRECT:
#   julia --project=<env> tools/quaddirect-probe.jl
using QuadDIRECT

function run_case(name, f, splits, lower, upper; kwargs...)
    evals = Tuple{Vector{Float64},Float64}[]
    lf = x -> (v = f(x); push!(evals, (copy(x), v)); v)
    root, x0 = QuadDIRECT.analyze(lf, splits, lower, upper; kwargs...)
    box = minimum(root)
    println("CASE ", name)
    for (x, v) in evals
        println("E ", join(repr.(x), " "), " ", repr(v))
    end
    println("U ", join(repr.(QuadDIRECT.position(box, x0)), " "))
    println("O ", repr(QuadDIRECT.value(box)))
    println("N ", length(evals))
end

rosen(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
s2 = ([-0.5, 0.0, 0.5], [-0.5, 0.0, 0.5])
run_case("rosen2", rosen, s2, [-1.0, -1.0], [0.8, 0.8])
run_case("rosen_maxevals", rosen, s2, [-1.0, -1.0], [0.8, 0.8]; maxevals = 100, rtol = 1e-6)
canyon(x) = 0.1 * (x[1] + x[2])^2 + 10 * (x[1] - x[2])^2
run_case("canyon", canyon, ([-11.0, -10.0, -9.0], [-7.0, -6.0, -5.0]), [-Inf, -Inf], [Inf, Inf]; atol = 0.01)
quad3(x) = (x[1] - 0.3)^2 + 2 * (x[2] + 0.2)^2 + 3 * (x[3] - 0.1)^2 + 0.5 * x[1] * x[2] - 0.3 * x[2] * x[3] +
           0.4 * x[1] * x[2] * x[3]
run_case("quad3", quad3, ([-0.5, 0.0, 0.5], [-0.5, 0.0, 0.5], [-0.5, 0.0, 0.5]), fill(-1.0, 3), fill(1.0, 3);
         maxevals = 400)
bound2(x) = (x[1] - 2)^2 + (x[2] + 0.5)^2 + x[1] * x[2]
run_case("bound2", bound2, s2, [-1.0, -1.0], [1.0, 1.0])
oned(x) = (x[1] - 0.7)^2 * (x[1] + 0.3)^2 + 0.1 * x[1]
run_case("oned", oned, ([-1.0, 0.0, 1.0],), [-2.0], [2.0])
camel(x) = (4 - 2.1 * x[1]^2 + x[1]^2 * x[1]^2 / 3) * x[1]^2 + x[1] * x[2] + (-4 + 4 * x[2]^2) * x[2]^2
run_case("camel", camel, ([-2.0, 0.0, 2.0], [-1.0, 0.0, 1.0]), [-3.0, -2.0], [3.0, 2.0])
# A run through the branch of complete_quadratic_model! where the box point
# is an end of the new split. QuadDIRECT gives the box value to the middle
# point there; the R port passes the box point instead, so its evaluations
# differ from the 79th on (the minimum found is the same).
const Hedge = [2.2033 -0.8346 -0.1717; -0.8346 1.6376 -0.1638; -0.1717 -0.1638 1.393]
const cedge = [0.5042, 1.2136, -0.5949]
function edge(x)
    d = x - cedge; s = 0.0
    for i = 1:3, j = 1:3
        s = s + d[i] * (Hedge[i, j] * d[j])
    end
    s
end
run_case("edge", edge, ([-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]), fill(-1.0, 3), fill(1.0, 3);
         maxevals = 300)
