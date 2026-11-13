# Reference values for tests/testthat/test-convex.R (ConvexOptimization.jl
# 0.1.4 with Clarabel.jl 0.11.1, MathOptInterface.jl 1.54.0 and
# SymbolicAnalysis.jl 0.5.8, the newest releases). It builds its own
# environment with these versions; the resolved versions are printed to
# stderr. Run from the repository root:
#   julia tools/convex-probe.jl > convex.tsv
#   REF=convex.tsv OUT=tests/testthat/test-convex.R Rscript tools/convex-gen-test.R
# One line per case: name, retcode, Clarabel iterations, objective, u, the
# duals (one comma separated group per ConeConstraint, groups separated by
# ';'), the cones Clarabel received after presolve and chordal
# decomposition, n, m, the lowered atom sets in order, the AMD permutation of
# the KKT matrix, whether a chordal decomposition was applied and the active
# bridges. A case that throws prints name, ERROR and the message. Floats are
# printed with %a so that they are read back exactly.
using Pkg
Pkg.activate(; temp = true)
Pkg.add([
    PackageSpec(name = "ConvexOptimization", version = "0.1.4"),
    PackageSpec(name = "SymbolicAnalysis", version = "0.5.8"),
    PackageSpec(name = "Clarabel", version = "0.11.1"),
    PackageSpec(name = "MathOptInterface", version = "1.54.0"),
    PackageSpec(name = "SciMLBase"),
]; io = devnull)
Pkg.status(; io = stderr)
using ConvexOptimization, SciMLBase, Clarabel, LinearAlgebra, Printf, SparseArrays
using SciMLBase: ConvexOptimizationProblem
import MathOptInterface as MOI
import SymbolicAnalysis

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

# MOI.RawSolver through a CachingOptimizer tries to map the solver object
# through the index map and throws, so unwrap the layers by hand
rawsolver(m::Clarabel.Optimizer) = m.solver
rawsolver(m::MOI.Bridges.LazyBridgeOptimizer) = rawsolver(m.model)
rawsolver(m::MOI.Utilities.CachingOptimizer) = rawsolver(m.optimizer)
rawsolver(m) = MOI.get(m, MOI.RawSolver())
# a case is passed as a function building the problem, so that an error in
# the constructor is reported as that case's ERROR line
build(prob) = prob isa Function ? prob() : prob
function report(name, cache, sol)
    solver = rawsolver(cache.model)
    perm = solver.kktsystem.kktsolver.ldlsolver.factors.perm
    println(join([
        name, string(sol.retcode), string(solver.info.iterations), hx(sol.objective), hv(sol.u),
        join(hv.(sol.dual), ";"), join(conestr.(solver.data.cones), ","),
        string(solver.data.n), string(solver.data.m),
        join(setstr.(cache.dpp.atomsets), ","), join(perm, ","),
        isnothing(solver.data.chordal_info) ? "none" : "chordal",
        bridges_str(cache.model),
    ], "\t"))
end
function bridges_str(model)
    io = IOBuffer()
    try
        MOI.Bridges.print_active_bridges(io, model)
    catch e
        print(io, "unavailable: ", sprint(showerror, e))
    end
    return join(filter(!isempty, strip.(split(String(take!(io)), "\n"))), " | ")
end
errline(name, e) = println(join([name, "ERROR", replace(sprint(showerror, e), "\n" => " ", "\t" => " ")], "\t"))
function probe(name, prob; alg = ALG)
    try
        cache = SciMLBase.init(build(prob), alg)
        sol = SciMLBase.solve!(cache)
        report(name, cache, sol)
        name in TRACE_CASES && trace(name, prob, alg, Int(rawsolver(cache.model).info.iterations))
    catch e
        e isa InterruptException && rethrow()
        errline(name, e)
    end
end

# Per iteration trace of the cases whose iterates drift from the R port. The
# solver loop has no callback, so iterate k is the result of a fresh solve
# with max_iter = k; Clarabel is deterministic, so this is the k-th iterate of
# the full solve. Written to the file named by TRACE (default
# convex-trace.tsv): one DATA line per case with the scaled problem data, then
# per k the status, mu, sigma, step length, primal and dual cost, primal and
# dual residuals, absolute and relative gap, ktratio, the infeasibility
# residuals and the unscaled x, z and s.
const TRACE_CASES = Set(["pow", "expcon", "geomean", "relentr", "con_explog", "eig_nested", "logdet_con", "logdet",
    "rootdet_tri_nocd", "rootdet_sq_nocd"])
const TRACE_PATH = get(ENV, "TRACE", "convex-trace.tsv")
isfile(TRACE_PATH) && rm(TRACE_PATH)
function with_max_iter(alg, k)
    c = alg.optimizer_constructor
    c isa MOI.OptimizerWithAttributes ||
        return ConvexMOI(MOI.OptimizerWithAttributes(c, "max_iter" => k))
    return ConvexMOI(MOI.OptimizerWithAttributes(c.optimizer_constructor, c.params...,
        MOI.RawOptimizerAttribute("max_iter") => k))
end
nz(M) = M isa SparseMatrixCSC ? nonzeros(M) : nonzeros(parent(M))
function trace(name, prob, alg, iters)
    open(TRACE_PATH, "a") do io
        for k in 0:iters
            try
                cache = SciMLBase.init(build(prob), with_max_iter(alg, k))
                # the solver state is read even if the wrapper rejects an early exit
                try
                    SciMLBase.solve!(cache)
                catch e
                    e isa InterruptException && rethrow()
                end
                s = rawsolver(cache.model)
                d = s.data
                eq = d.equilibration
                k == 0 && println(io, join([name, "DATA", hv(d.q), hv(d.b), hv(nz(d.A)), hv(nz(d.P)),
                    hv(eq.d), hv(eq.e), hx(eq.c[])], "\t"))
                i = s.info
                println(io, join([name, string(k), string(i.status),
                    hv([i.mu, i.sigma, i.step_length, i.cost_primal, i.cost_dual, i.res_primal, i.res_dual,
                        i.gap_abs, i.gap_rel, i.ktratio, i.res_primal_inf, i.res_dual_inf]),
                    hv(s.solution.x), hv(s.solution.z), hv(s.solution.s)], "\t"))
                # the last scaling of the nonsymmetric cones and the last
                # step direction, to find where the iterates first part
                if 1 <= k <= 4
                    st = s.step_lhs
                    println(io, join([name, "STEP" * string(k), hv(st.x), hv(st.z), hv(st.s), hx(st.τ), hx(st.κ)], "\t"))
                    for (j, K) in enumerate(s.cones)
                        hasproperty(K, :H_dual) || continue
                        println(io, join([name, "CONE" * string(k) * ":" * string(j), hv(K.H_dual), hv(K.Hs),
                            hv(K.grad), hv(K.z)], "\t"))
                    end
                end
            catch e
                e isa InterruptException && rethrow()
                println(io, join([name, string(k), "ERROR", replace(sprint(showerror, e), "\n" => " ", "\t" => " ")], "\t"))
            end
        end
    end
end
function sweep(name, prob, thetas; alg = ALG)
    cache = try
        SciMLBase.init(build(prob), alg)
    catch e
        e isa InterruptException && rethrow()
        return errline(name, e)
    end
    for (k, th) in enumerate(thetas)
        try
            cache = SciMLBase.reinit!(cache; p = th)
            sol = SciMLBase.solve!(cache)
            report("$(name)_$k", cache, sol)
        catch e
            e isa InterruptException && rethrow()
            errline("$(name)_$k", e)
        end
    end
end

cons = [
    ConeConstraint((u, p) -> [u[1] + u[2] - 1.0], MOI.Zeros(1)),
    ConeConstraint((u, p) -> [u[1], u[2]], MOI.Nonnegatives(2)),
]
c1 = [ConeConstraint((u, p) -> [u[1] + u[2] - 1.0], MOI.Zeros(1))]
A = reshape([1.0, 1.0], 2, 1)
b = [0.0, 2.0]
probe("lp", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + 2 * u[2]), [0.5, 0.5]; constraints = cons))
probe("soc", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[3]), [0.0, 0.0, 0.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] - 3.0, u[2] - 4.0], MOI.Zeros(2)),
    ConeConstraint((u, p) -> [u[3], u[1], u[2]], MOI.SecondOrderCone(3))]))
probe("norm2", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(A * u - b, 2)), [0.0]))
probe("norm2_cons", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(u .- 1.0, 2)), [0.5, 0.5]; constraints = cons))
probe("norm1", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(u .- 1.0, 1)), [0.5, 0.5]; constraints = c1))
probe("norminf", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(u .- 1.0, Inf)), [0.5, 0.5]; constraints = c1))
probe("exp", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> exp(u[1]) + u[2]), [0.0, 0.0];
    constraints = [ConeConstraint((u, p) -> [u[1] - 1.0, u[2]], MOI.Zeros(2))]))
probe("log", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -log(u[1])), [1.0];
    constraints = [ConeConstraint((u, p) -> [u[1] - 2.0], MOI.Zeros(1))]))
# cones and bridges beyond the upstream tests
probe("rsoc", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + u[2]), [1.0, 1.0, 0.0]; constraints = [
    ConeConstraint((u, p) -> [u[3] - 2.0], MOI.Zeros(1)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3]], MOI.RotatedSecondOrderCone(3))]))
probe("nonpos", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1] - u[2]), [0.0, 0.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] + 2 * u[2] - 4.0, 3 * u[1] + u[2] - 6.0], MOI.Nonpositives(2)),
    ConeConstraint((u, p) -> [u[1], u[2]], MOI.Nonnegatives(2))]))
probe("soc6", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[6]), zeros(6); constraints = [
    ConeConstraint((u, p) -> [u[1] - 1.0, u[2] - 2.0, u[3] - 3.0, u[4] - 4.0, u[5] - 5.0], MOI.Zeros(5)),
    ConeConstraint((u, p) -> [u[6], u[1], u[2], u[3], u[4], u[5]], MOI.SecondOrderCone(6))]))
probe("pow", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[3]), [1.0, 1.0, 0.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] - 1.0, u[2] - 4.0], MOI.Zeros(2)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3]], MOI.PowerCone(0.3))]))
probe("genpow", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[3] - u[4]), [1.0, 1.0, 0.0, 0.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] - 1.0, u[2] - 2.0], MOI.Zeros(2)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3], u[4]], Clarabel.MOI.GenPowerCone([0.4, 0.6], 2))]))
probe("psd", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + u[3]), [1.0, 0.0, 1.0]; constraints = [
    ConeConstraint((u, p) -> [u[2] - 1.0], MOI.Zeros(1)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3]], MOI.PositiveSemidefiniteConeTriangle(2))]))
probe("psd3", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + u[3] + u[6]), [1.0, 0.0, 1.0, 0.0, 0.0, 1.0]; constraints = [
    ConeConstraint((u, p) -> [u[2] - 1.0, u[4] - 0.5, u[5] + 0.25], MOI.Zeros(3)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3], u[4], u[5], u[6]], MOI.PositiveSemidefiniteConeTriangle(3))]))
probe("expcon", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[3]), [0.0, 1.0, 1.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] - 0.5, u[2] - 2.0], MOI.Zeros(2)),
    ConeConstraint((u, p) -> [u[1], u[2], u[3]], MOI.ExponentialCone())]))
probe("lpbox", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] - 2 * u[2]), [0.0, 0.0];
    constraints = c1, lb = [-1.0, -3.0], ub = [2.0, 1.5]))
probe("infeas", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), [0.0]; constraints = [
    ConeConstraint((u, p) -> [u[1] - 1.0], MOI.Nonnegatives(1)), ConeConstraint((u, p) -> [-u[1]], MOI.Nonnegatives(1))]))

# parametric sweeps of the upstream tests
lp_cons = [
    ConeConstraint((u, p) -> [u[1] + u[2] - 1.0 - p[2]], MOI.Zeros(1)),
    ConeConstraint((u, p) -> [u[1], u[2]], MOI.Nonnegatives(2)),
]
sweep("lp_p", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> p[1] * u[1] + u[2]), [0.5, 0.5], [0.5, 0.0];
    constraints = lp_cons), ([0.5, 0.0], [2.0, 1.0], [-1.0, 0.5], [0.0, 3.0], [0.5, 0.0]))
Al = [1.0 0.5; 0.2 1.0; 0.7 0.7]
bl = [1.0, 2.0, 0.5]
sweep("lasso", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(Al * u - bl, 2) + p[1] * norm(u, 1)),
    [0.0, 0.0], [0.0]; lb = [-10.0, -10.0], ub = [10.0, 10.0]), ([0.0], [0.25], [0.5], [2.0], [10.0]))
sweep("proj", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(u .- p, 2)), [0.0, 0.0], [1.0, 1.0];
    constraints = c1), ([1.0, 1.0], [3.0, -1.0], [0.0, 0.0]))
sweep("cross", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> p[1] * u[1] + u[2]), [0.5, 0.5], [0.0];
    constraints = c1, lb = [-10.0, -10.0], ub = [10.0, 10.0]), ([0.0], [3.0], [-5.0], [0.0]))
sweep("maxsense", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> p[1] * u[1] + u[2] - norm(u .- 1.0, 2)),
    [0.0, 0.0], [1.0]; sense = SciMLBase.MaxSense, lb = [-5.0, -5.0], ub = [5.0, 5.0]), ([1.0], [-2.0], [0.0]))
sweep("layout", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(u .- p[1:2], 2) + p[3] * u[1]),
    [0.0, 0.0], [1.0, 1.0, 0.0]; constraints = [
        ConeConstraint((u, p) -> [u[1] + u[2] - p[3]], MOI.Zeros(1)),
        ConeConstraint((u, p) -> [u[1], u[2]], MOI.Nonnegatives(2))],
    lb = [-10.0, -10.0], ub = [10.0, 10.0]), ([1.0, 1.0, 0.0], [3.0, -1.0, 2.0], [0.0, 0.0, -4.0]))
sweep("failed", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> p[1] * u[1] + u[2]), [0.0, 0.0], [1.0];
    constraints = [ConeConstraint((u, p) -> [u[1], u[2]], MOI.Nonnegatives(2))]), ([1.0], [-1.0], [2.0]))

# ConvexOptimization 0.1.4: the new atoms
const ALG_TIGHT = ConvexMOI(MOI.OptimizerWithAttributes(Clarabel.Optimizer, "tol_gap_abs" => 1.0e-10,
    "tol_gap_rel" => 1.0e-10, "tol_feas" => 1.0e-10, "tol_ktratio" => 1.0e-10))
QA = Float64[1 0 2; 0 1 1; 1 1 0; 2 0 1]
QB = Float64[1, 2, 3, 4]
P2 = [2.0 0.5; 0.5 1.0]
Xd(u) = [u[1] 0.0; 0.0 1.0 - u[1]]
probe("abs", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> abs(u[1] - 2.0) + abs(u[2] + 1.0)), [0.0, 0.0]))
probe("lad", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> sum(abs.(Al * u - bl))), [0.0, 0.0]))
probe("cheb", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> maximum(abs.(Al * u - bl))), [0.0, 0.0]))
probe("max4", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> max(u[1] - 1.0, 1.0 - u[1], u[2] + 0.5, -u[2])), [0.0, 0.0]))
probe("minimum_max", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> minimum(u)), [0.5, 0.5];
    constraints = c1, sense = SciMLBase.MaxSense))
probe("minimum_min", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> minimum(u)), [0.5, 0.5]; constraints = c1))
probe("minabs", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> minimum(abs.(u))), [0.5, 0.5]))
probe("sumsq", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> sum((QA * u - QB) .^ 2)), zeros(3)))
probe("normsq", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> norm(QA * u - QB)^2), zeros(3)))
probe("square", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> (u[1] - 1.0)^2 + (u[2] + 2.0)^2 + u[1]), [0.0, 0.0]))
probe("quadform", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u' * P2 * u + u[1] - 2 * u[2]), [0.0, 0.0]))
probe("quad_form_fn", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> SymbolicAnalysis.quad_form(u, P2) + u[1]), [0.0, 0.0]))
probe("selfprod", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> (QA * u - QB)' * (QA * u - QB)), zeros(3)))
probe("quad_nonpsd", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u' * [1.0 0.0; 0.0 -1.0] * u), [0.0, 0.0]))
probe("eigmax", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> eigmax(Xd(u))), [0.0]); alg = ALG_TIGHT)
probe("eigmin", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> eigmin(Xd(u))), [0.0]; sense = SciMLBase.MaxSense); alg = ALG_TIGHT)
probe("opnorm", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> opnorm([u[1] - 1.0 0.0 0.0; 0.0 u[2] + 2.0 0.0])), [0.0, 0.0]); alg = ALG_TIGHT)
probe("logdet", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> logdet([u[1] 0.0; 0.0 u[2]])), [1.0, 1.0];
    sense = SciMLBase.MaxSense, constraints = [ConeConstraint((u, p) -> [2.0 - u[1] - u[2]], MOI.Nonnegatives(1))]); alg = ALG_TIGHT)
probe("eig_asym", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> eigmax([u[1] u[1]; 0.0 1.0])), [0.0]))
probe("eig_scale", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> p[1] * eigmax(Xd(u))), [0.0], [1.0]))
probe("eig_nested", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> exp(eigmax(Xd(u)))), [0.0]); alg = ALG_TIGHT)
probe("eig_max", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> max(eigmax(Xd(u)), 1.0)), [0.0]); alg = ALG_TIGHT)
sweep("eig_pdata", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> eigmax(Xd(u)) + eigmax([p[1] 0.0; 0.0 -p[1]])),
    [0.0], [-3.0]), ([-3.0], [1.0], [2.0]); alg = ALG_TIGHT)
sweep("eig_param", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> eigmax([u[1] + p[1] 0.0; 0.0 1.0 - u[1] + p[2]])),
    [0.0], [0.0, 0.0]), ([0.0, 0.0], [1.0, -0.5], [-2.0, 3.0]); alg = ALG_TIGHT)
# nested atoms and parameter-only data
probe("nest_exp", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> exp(norm(u .- [1.0, 2.0]))), [0.0, 0.0]))
probe("nest_max", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> max(norm(u), 1.0)), [0.4, 0.6]; constraints = c1))
probe("nest_sumabs_sq", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> sum(abs.(u))^2), [0.4, 0.6]; constraints = c1))
probe("nest_abs_norm", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> abs(norm(u .- 1.0))), [0.0, 0.0]))
probe("nest_box", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> max(u[1], abs(u[2]))), [0.3, -0.4];
    lb = [-1.0, -1.0], ub = [1.0, 1.0]))
sweep("nest_p", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> exp(norm(u .- p))), [0.0, 0.0], [1.0, 2.0]),
    ([1.0, 2.0], [-1.0, 0.5], [0.0, 0.0]))
sweep("pdata_abs", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> abs(u[1] + abs(p[1]) - 1) + u[1]^2), [0.0], [2.0]),
    ([2.0], [-0.25], [0.0], [-1.0]); alg = ALG_TIGHT)
sweep("pdata_exp", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> abs(u[1] + exp(p[1]) - 2) + u[1]^2), [0.0], [0.0]),
    ([0.0], [log(2.25)], [log(1.5)]); alg = ALG_TIGHT)
sweep("pdata_max", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> max(p[1], p[2]) * 1 + norm(u .- p)), [0.0, 0.0], [1.0, -2.0]),
    ([1.0, -2.0], [-2.0, 3.0], [0.0, 0.0]); alg = ALG_TIGHT)
probe("scalar_p", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), [0.0], 1.0; lb = [0.0], ub = [10.0]))
# atoms in constraints
probe("con_norm", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1] - u[2]), [0.0, 0.0];
    constraints = [ConeConstraint((u, p) -> [norm(u) - 1.0], MOI.Nonpositives(1))]))
probe("con_l1", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1] - 2 * u[2]), [0.0, 0.0];
    constraints = [ConeConstraint((u, p) -> [sum(abs.(u)) - 1.0], MOI.Nonpositives(1))]))
probe("con_explog", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1] - u[2]), [0.0, 0.0];
    constraints = [ConeConstraint((u, p) -> [exp(u[1]) - 2.0, u[2]^2 - 4.0], MOI.Nonpositives(2)),
        ConeConstraint((u, p) -> [log(u[1] + 2.0)], MOI.Nonnegatives(1))]))
probe("con_quad", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1]), [0.0, 0.0];
    constraints = [ConeConstraint((u, p) -> [u' * P2 * u - 1.0], MOI.Nonpositives(1))]))
probe("con_minimum", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + u[2]), [1.0, 1.0];
    constraints = [ConeConstraint((u, p) -> [minimum(u) - 0.5], MOI.Nonnegatives(1))]))
probe("con_two", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1] - u[2]), [0.0, 0.0];
    constraints = [ConeConstraint((u, p) -> [norm(u) + abs(u[1] - u[2]) - 2.0], MOI.Nonpositives(1))]))
probe("con_wrong", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), [0.0, 0.0];
    constraints = [ConeConstraint((u, p) -> [norm(u) - 1.0], MOI.Nonnegatives(1))]))
probe("con_eq", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), [0.0, 0.0];
    constraints = [ConeConstraint((u, p) -> [norm(u) - 1.0], MOI.Zeros(1))]))
probe("con_sign", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), [0.0, 0.0], [1.0];
    constraints = [ConeConstraint((u, p) -> [p[1] - norm(u)], MOI.Nonpositives(1))]))
sweep("con_p", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + u[2]), [0.0, 0.0], [1.0, 1.0, 1.0];
    constraints = [ConeConstraint((u, p) -> [norm(u .- p[1:2]) - p[3]], MOI.Nonpositives(1))]),
    ([1.0, 1.0, 1.0], [0.0, 2.0, 0.5], [-1.0, 0.0, 2.0]))
probe("logdet_con", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + u[2]), [1.0, 1.0];
    constraints = [ConeConstraint((u, p) -> [logdet([u[1] 0.0; 0.0 u[2]])], MOI.Nonnegatives(1))]);
    alg = ConvexMOI(MOI.OptimizerWithAttributes(Clarabel.Optimizer, "tol_gap_abs" => 1.0e-12, "tol_gap_rel" => 1.0e-12,
        "tol_feas" => 1.0e-12, "tol_ktratio" => 1.0e-12)))

# MathOptInterface sets reached through bridges
cz(vals) = ConeConstraint((u, p) -> [u[i] - v for (i, v) in vals], MOI.Zeros(length(vals)))
probe("geomean", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1]), [1.0, 1.0, 1.0];
    constraints = [cz([(2, 2.0), (3, 8.0)]), ConeConstraint((u, p) -> [u[1], u[2], u[3]], MOI.GeometricMeanCone(3))]))
probe("geomean4", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1]), ones(4);
    constraints = [cz([(2, 1.0), (3, 2.0), (4, 4.0)]), ConeConstraint((u, p) -> [u[1], u[2], u[3], u[4]], MOI.GeometricMeanCone(4))]))
probe("relentr", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), ones(5);
    constraints = [cz([(2, 1.0), (3, 2.0), (4, 2.0), (5, 1.0)]), ConeConstraint((u, p) -> [u[i] for i in 1:5], MOI.RelativeEntropyCone(5))]))
probe("normcone2", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), ones(3);
    constraints = [cz([(2, 3.0), (3, 4.0)]), ConeConstraint((u, p) -> [u[i] for i in 1:3], MOI.NormCone(2.0, 3))]))
probe("normcone3", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), ones(3);
    constraints = [cz([(2, 3.0), (3, 4.0)]), ConeConstraint((u, p) -> [u[i] for i in 1:3], MOI.NormCone(3.0, 3))]))
probe("normone_set", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), ones(3);
    constraints = [cz([(2, 3.0), (3, -4.0)]), ConeConstraint((u, p) -> [u[i] for i in 1:3], MOI.NormOneCone(3))]))
probe("norminf_set", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), ones(3);
    constraints = [cz([(2, 3.0), (3, -4.0)]), ConeConstraint((u, p) -> [u[i] for i in 1:3], MOI.NormInfinityCone(3))]))
Mf = [2.0, 0.5, 1.0, 0.5, 1.0, 3.0]   # a 2x3 matrix by columns
# Clarabel's PSD completion of the dual indexes an OrderedSet and throws a
# MethodError whenever a cone is decomposed, so each such case also runs
# without the completion to give the decomposed solve a reference
const ALG_NOCD = ConvexMOI(MOI.OptimizerWithAttributes(Clarabel.Optimizer,
    "chordal_decomposition_complete_dual" => false))
probe_both(name, prob) = (probe(name, prob); probe(name * "_nocd", prob; alg = ALG_NOCD))
probe_both("normspec", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), ones(7);
    constraints = [cz([(i + 1, Mf[i]) for i in 1:6]), ConeConstraint((u, p) -> [u[i] for i in 1:7], MOI.NormSpectralCone(2, 3))]))
probe("normnuc", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), ones(7);
    constraints = [cz([(i + 1, Mf[i]) for i in 1:6]), ConeConstraint((u, p) -> [u[i] for i in 1:7], MOI.NormNuclearCone(2, 3))]))
probe("logdet_tri", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1]), ones(5);
    constraints = [cz([(2, 1.0), (3, 2.0), (4, 0.5), (5, 1.0)]), ConeConstraint((u, p) -> [u[i] for i in 1:5], MOI.LogDetConeTriangle(2))]))
probe_both("logdet_sq", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1]), ones(6);
    constraints = [cz([(2, 1.0), (3, 2.0), (4, 0.5), (5, 0.5), (6, 1.0)]), ConeConstraint((u, p) -> [u[i] for i in 1:6], MOI.LogDetConeSquare(2))]))
probe_both("rootdet_tri", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1]), ones(4);
    constraints = [cz([(2, 2.0), (3, 0.5), (4, 1.0)]), ConeConstraint((u, p) -> [u[i] for i in 1:4], MOI.RootDetConeTriangle(2))]))
probe_both("rootdet_sq", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1]), ones(5);
    constraints = [cz([(2, 2.0), (3, 0.5), (4, 0.5), (5, 1.0)]), ConeConstraint((u, p) -> [u[i] for i in 1:5], MOI.RootDetConeSquare(2))]))
probe("psdsq", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1] + u[4]), [1.0, 0.0, 0.0, 1.0];
    constraints = [cz([(2, 1.0), (3, 1.0)]), ConeConstraint((u, p) -> [u[i] for i in 1:4], MOI.PositiveSemidefiniteConeSquare(2))]))
probe("hyperrect", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> -u[1] + u[2]), [0.0, 0.0];
    constraints = [ConeConstraint((u, p) -> [u[1], u[2]], MOI.HyperRectangle([0.0, -1.0], [1.0, 2.0]))]))
probe("dualexp", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[3]), ones(3);
    constraints = [ConeConstraint((u, p) -> [u[i] for i in 1:3], MOI.DualExponentialCone())]))
probe("dualpow", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[3]), ones(3);
    constraints = [ConeConstraint((u, p) -> [u[i] for i in 1:3], MOI.DualPowerCone(0.3))]))
probe("dualrelentr", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), ones(3);
    constraints = [ConeConstraint((u, p) -> [u[i] for i in 1:3], MOI.DualRelativeEntropyCone(3))]))
probe("dualgeomean", () -> ConvexOptimizationProblem(OptimizationFunction((u, p) -> u[1]), ones(3);
    constraints = [cz([(2, 1.0), (3, 2.0)]), ConeConstraint((u, p) -> [u[i] for i in 1:3], MOI.DualGeometricMeanCone(3))]))

# sparse PSD cones: Clarabel's chordal decomposition (PSDSquare keeps it on)
function bandsq(u, n)
    X = Matrix{Any}(undef, n, n)
    for j in 1:n, i in 1:n
        X[i, j] = i == j ? u[i] : abs(i - j) == 1 ? u[n + min(i, j)] : 0.0
    end
    return vec(X)
end
function arrowsq(u, n)
    X = Matrix{Any}(undef, n, n)
    for j in 1:n, i in 1:n
        X[i, j] = i == j ? u[i] : (i == 1 || j == 1) ? u[n + max(i, j) - 1] : 0.0
    end
    return vec(X)
end
# band of half width w; the entries of diagonal k follow those of diagonals
# 1..k-1. Width 3 gives 4-cliques overlapping in 3 (merge weight
# 64 + 64 - 125 >= 0) and width 4 gives 5-cliques overlapping in 4, so the
# clique graph merges run, with ties between equal weights
function wbandsq(u, n, w)
    X = Matrix{Any}(undef, n, n)
    for j in 1:n, i in 1:n
        k = abs(i - j)
        X[i, j] = k == 0 ? u[i] : k <= w ? u[n + sum(n - t for t in 1:(k - 1); init = 0) + min(i, j)] : 0.0
    end
    return vec(X)
end
nband(n, w) = sum(n - t for t in 1:w)
for (nm, f, n, nv) in (("band6", bandsq, 6, 5), ("arrow6", arrowsq, 6, 5), ("band9", bandsq, 9, 8),
        ("wband8", (u, n) -> wbandsq(u, n, 3), 8, nband(8, 3)), ("wband12", (u, n) -> wbandsq(u, n, 4), 12, nband(12, 4)))
    prob = ConvexOptimizationProblem(OptimizationFunction((u, p) -> sum(u[i] for i in 1:n)), [ones(n); zeros(nv)];
        constraints = [ConeConstraint((u, p) -> [u[n + k] - 1.0 for k in 1:nv], MOI.Zeros(nv)),
            ConeConstraint((u, p) -> f(u, n), MOI.PositiveSemidefiniteConeSquare(n))])
    probe("chordal_" * nm, prob)
    for (mm, merge) in (("pc", :parent_child), ("nm", :none))
        probe("chordal_$(nm)_$mm", () -> prob; alg = ConvexMOI(MOI.OptimizerWithAttributes(Clarabel.Optimizer,
            "chordal_decomposition_merge_method" => merge)))
    end
    probe("chordal_$(nm)_std", () -> prob; alg = ConvexMOI(MOI.OptimizerWithAttributes(Clarabel.Optimizer,
        "chordal_decomposition_compact" => false)))
    probe("chordal_$(nm)_nocd", () -> prob; alg = ConvexMOI(MOI.OptimizerWithAttributes(Clarabel.Optimizer,
        "chordal_decomposition_complete_dual" => false)))
    probe("chordal_$(nm)_nm_nocd", () -> prob; alg = ConvexMOI(MOI.OptimizerWithAttributes(Clarabel.Optimizer,
        "chordal_decomposition_merge_method" => :none, "chordal_decomposition_complete_dual" => false)))
    probe("chordal_$(nm)_std_nocd", () -> prob; alg = ConvexMOI(MOI.OptimizerWithAttributes(Clarabel.Optimizer,
        "chordal_decomposition_compact" => false, "chordal_decomposition_complete_dual" => false)))
end
