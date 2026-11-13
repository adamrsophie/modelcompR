# Reference values for tests/testthat/test-convex.R (ConvexOptimization.jl
# 0.1.0 with Clarabel.jl 0.11.1). It builds its own environment with the
# pinned versions (SymbolicAnalysis 0.3 follows from ConvexOptimization's
# compat). Run from the repository root:
#   julia tools/convex-probe.jl > convex.tsv
#   REF=convex.tsv OUT=tests/testthat/test-convex.R Rscript tools/convex-gen-test.R
# One line per case: name, retcode, Clarabel iterations, objective, u, the
# duals (one comma separated group per ConeConstraint, groups separated by
# ';'), the cones Clarabel received after presolve, n, m, the lowered atom
# sets in order and the AMD permutation of the KKT matrix. Floats are printed
# with %a so that they are read back exactly.
using Pkg
Pkg.activate(; temp = true)
Pkg.add([
    PackageSpec(name = "ConvexOptimization", version = "0.1.0"),
    PackageSpec(name = "Clarabel", version = "0.11.1"),
    PackageSpec(name = "MathOptInterface", version = "1.54.0"),
    PackageSpec(name = "SciMLBase"),
]; io = devnull)
using ConvexOptimization, SciMLBase, Clarabel, LinearAlgebra, Printf
using SciMLBase: ConvexOptimizationProblem
import MathOptInterface as MOI

const ALG = ConvexMOI(Clarabel.Optimizer)
hx(x) = isnan(x) ? "NaN" : isinf(x) ? (x > 0 ? "Inf" : "-Inf") : @sprintf("%a", x)
hv(v) = join(hx.(v), ",")
conestr(c::Clarabel.ZeroConeT) = "Zero:$(c.dim)"
conestr(c::Clarabel.NonnegativeConeT) = "Nonnegative:$(c.dim)"
conestr(c::Clarabel.SecondOrderConeT) = "SecondOrder:$(c.dim)"
conestr(c::Clarabel.ExponentialConeT) = "Exponential"
conestr(c::Clarabel.PowerConeT) = "Power"
conestr(c::Clarabel.GenPowerConeT) = "GenPower:$(length(getfield(c, Symbol("\u03b1")))):$(c.dim2)"
conestr(c::Clarabel.PSDTriangleConeT) = "PSDTriangle:$(c.dim)"
setstr(s) = string(nameof(typeof(s)), ":", MOI.dimension(s))

function report(name, cache, sol)
    solver = MOI.get(cache.model, MOI.RawSolver())
    perm = solver.kktsystem.kktsolver.ldlsolver.factors.perm
    println(join([
        name, string(sol.retcode), string(solver.info.iterations), hx(sol.objective), hv(sol.u),
        join(hv.(sol.dual), ";"), join(conestr.(solver.data.cones), ","),
        string(solver.data.n), string(solver.data.m),
        join(setstr.(cache.dpp.atomsets), ","), join(perm, ","),
    ], "\t"))
end
function probe(name, prob)
    cache = SciMLBase.init(prob, ALG)
    sol = SciMLBase.solve!(cache)
    report(name, cache, sol)
end
function sweep(name, prob, thetas)
    cache = SciMLBase.init(prob, ALG)
    for (k, th) in enumerate(thetas)
        cache = SciMLBase.reinit!(cache; p = th)
        sol = SciMLBase.solve!(cache)
        report("$(name)_$k", cache, sol)
    end
end

cons = [
    ConeConstraint((u, p) -> [u[1] + u[2] - 1.0], MOI.Zeros(1)),
    ConeConstraint((u, p) -> [u[1], u[2]], MOI.Nonnegatives(2)),
]
c1 = [ConeConstraint((u, p) -> [u[1] + u[2] - 1.0], MOI.Zeros(1))]
A = reshape([1.0, 1.0], 2, 1)
b = [0.0, 2.0]
probe("lp", ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + 2 * u[2]), [0.5, 0.5]; constraints = cons))
probe("soc", ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[3]), [0.0, 0.0, 0.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] - 3.0, u[2] - 4.0], MOI.Zeros(2)),
    ConeConstraint((u, p) -> [u[3], u[1], u[2]], MOI.SecondOrderCone(3))]))
probe("norm2", ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(A * u - b, 2)), [0.0]))
probe("norm2_cons", ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(u .- 1.0, 2)), [0.5, 0.5]; constraints = cons))
probe("norm1", ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(u .- 1.0, 1)), [0.5, 0.5]; constraints = c1))
probe("norminf", ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(u .- 1.0, Inf)), [0.5, 0.5]; constraints = c1))
probe("exp", ConvexOptimizationProblem(OptimizationFunction((u, p) -> exp(u[1]) + u[2]), [0.0, 0.0];
    constraints = [ConeConstraint((u, p) -> [u[1] - 1.0, u[2]], MOI.Zeros(2))]))
probe("log", ConvexOptimizationProblem(OptimizationFunction((u, p) -> -log(u[1])), [1.0];
    constraints = [ConeConstraint((u, p) -> [u[1] - 2.0], MOI.Zeros(1))]))
# cones and bridges beyond the upstream tests
probe("rsoc", ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + u[2]), [1.0, 1.0, 0.0]; constraints = [
    ConeConstraint((u, p) -> [u[3] - 2.0], MOI.Zeros(1)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3]], MOI.RotatedSecondOrderCone(3))]))
probe("nonpos", ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1] - u[2]), [0.0, 0.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] + 2 * u[2] - 4.0, 3 * u[1] + u[2] - 6.0], MOI.Nonpositives(2)),
    ConeConstraint((u, p) -> [u[1], u[2]], MOI.Nonnegatives(2))]))
probe("soc6", ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[6]), zeros(6); constraints = [
    ConeConstraint((u, p) -> [u[1] - 1.0, u[2] - 2.0, u[3] - 3.0, u[4] - 4.0, u[5] - 5.0], MOI.Zeros(5)),
    ConeConstraint((u, p) -> [u[6], u[1], u[2], u[3], u[4], u[5]], MOI.SecondOrderCone(6))]))
probe("pow", ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[3]), [1.0, 1.0, 0.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] - 1.0, u[2] - 4.0], MOI.Zeros(2)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3]], MOI.PowerCone(0.3))]))
probe("genpow", ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[3] - u[4]), [1.0, 1.0, 0.0, 0.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] - 1.0, u[2] - 2.0], MOI.Zeros(2)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3], u[4]], Clarabel.MOI.GenPowerCone([0.4, 0.6], 2))]))
probe("psd", ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + u[3]), [1.0, 0.0, 1.0]; constraints = [
    ConeConstraint((u, p) -> [u[2] - 1.0], MOI.Zeros(1)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3]], MOI.PositiveSemidefiniteConeTriangle(2))]))
probe("psd3", ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + u[3] + u[6]), [1.0, 0.0, 1.0, 0.0, 0.0, 1.0]; constraints = [
    ConeConstraint((u, p) -> [u[2] - 1.0, u[4] - 0.5, u[5] + 0.25], MOI.Zeros(3)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3], u[4], u[5], u[6]], MOI.PositiveSemidefiniteConeTriangle(3))]))
probe("expcon", ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[3]), [0.0, 1.0, 1.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] - 0.5, u[2] - 2.0], MOI.Zeros(2)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3]], MOI.ExponentialCone())]))
probe("lpbox", ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] - 2 * u[2]), [0.0, 0.0];
    constraints = c1, lb = [-1.0, -3.0], ub = [2.0, 1.5]))
probe("infeas", ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), [0.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] - 1.0], MOI.Nonnegatives(1)), ConeConstraint((u, p) -> [-u[1]], MOI.Nonnegatives(1))]))

# parametric sweeps of the upstream tests
lp_cons = [
    ConeConstraint((u, p) -> [u[1] + u[2] - 1.0 - p[2]], MOI.Zeros(1)),
    ConeConstraint((u, p) -> [u[1], u[2]], MOI.Nonnegatives(2)),
]
sweep("lp_p", ConvexOptimizationProblem(OptimizationFunction((u, p) -> p[1] * u[1] + u[2]), [0.5, 0.5], [0.5, 0.0];
    constraints = lp_cons), ([0.5, 0.0], [2.0, 1.0], [-1.0, 0.5], [0.0, 3.0], [0.5, 0.0]))
Al = [1.0 0.5; 0.2 1.0; 0.7 0.7]
bl = [1.0, 2.0, 0.5]
sweep("lasso", ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(Al * u - bl, 2) + p[1] * norm(u, 1)),
    [0.0, 0.0], [0.0]; lb = [-10.0, -10.0], ub = [10.0, 10.0]), ([0.0], [0.25], [0.5], [2.0], [10.0]))
sweep("proj", ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(u .- p, 2)), [0.0, 0.0], [1.0, 1.0];
    constraints = c1), ([1.0, 1.0], [3.0, -1.0], [0.0, 0.0]))
sweep("cross", ConvexOptimizationProblem(OptimizationFunction((u, p) -> p[1] * u[1] + u[2]), [0.5, 0.5], [0.0];
    constraints = c1, lb = [-10.0, -10.0], ub = [10.0, 10.0]), ([0.0], [3.0], [-5.0], [0.0]))
sweep("maxsense", ConvexOptimizationProblem(OptimizationFunction((u, p) -> p[1] * u[1] + u[2] - norm(u .- 1.0, 2)),
    [0.0, 0.0], [1.0]; sense = SciMLBase.MaxSense, lb = [-5.0, -5.0], ub = [5.0, 5.0]), ([1.0], [-2.0], [0.0]))
sweep("layout", ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(u .- p[1:2], 2) + p[3] * u[1]),
    [0.0, 0.0], [1.0, 1.0, 0.0]; constraints = [
        ConeConstraint((u, p) -> [u[1] + u[2] - p[3]], MOI.Zeros(1)),
        ConeConstraint((u, p) -> [u[1], u[2]], MOI.Nonnegatives(2))],
    lb = [-10.0, -10.0], ub = [10.0, 10.0]), ([1.0, 1.0, 0.0], [3.0, -1.0, 2.0], [0.0, 0.0, -4.0]))
sweep("failed", ConvexOptimizationProblem(OptimizationFunction((u, p) -> p[1] * u[1] + u[2]), [0.0, 0.0], [1.0];
    constraints = [ConeConstraint((u, p) -> [u[1], u[2]], MOI.Nonnegatives(2))]), ([1.0], [-1.0], [2.0]))
