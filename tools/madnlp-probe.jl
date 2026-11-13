# Iteration logs of MadNLP for MadNLPOptimizer (R/optimization-ipopt.R):
# HS071 from the default and a far start, the infeasible circle and
# Waechter and Biegler's example (restoration phase), and the options of
# MadNLPOptimizer that change the iterates, through
# OptimizationMadNLP.jl as SciML calls it (print_level INFO, exact
# derivatives from ForwardDiff). Usage (packages: Optimization
# OptimizationMadNLP MadNLP ADTypes ForwardDiff):
#   julia --project=tools/julia tools/madnlp-probe.jl > madnlp-probe.txt
import Optimization
import OptimizationMadNLP
import MadNLP
import ADTypes
import ForwardDiff

println("== versions")
for p in ("Optimization", "OptimizationMadNLP", "MadNLP")
    try
        println(p, " ", pkgversion(Base.require(Main, Symbol(p))))
    catch
        println(p, " ?")
    end
end
flush(stdout)

hs071_f(x, p) = x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
hs071_c(res, x, p) = (res .= [x[1] * x[2] * x[3] * x[4], sum(abs2, x)])
circle_f(x, p) = x[1]^2 + x[2]^2
circle_c(res, x, p) = (res .= [x[1]^2 + x[2]^2, x[2] - x[1]^2])
wb_f(x, p) = x[1]
wb_c(res, x, p) = (res .= [x[1]^2 - x[2] - 1, x[1] - x[3] - 0.5])

function run(name, f, c, x0, lcons, ucons, lb, ub, opt)
    println("== ", name)
    flush(stdout)
    try
        of = Optimization.OptimizationFunction(f, ADTypes.AutoForwardDiff(); cons = c)
        prob = Optimization.OptimizationProblem(of, x0, nothing; lcons = lcons, ucons = ucons, lb = lb, ub = ub)
        sol = Optimization.solve(prob, opt; verbose = verbose)
        flush(stdout)
        println("retcode ", sol.retcode)
        println("status ", sol.original.status)
        println("u ", repr(sol.u))
        println("objective ", repr(sol.objective))
        println("iterations ", sol.original.counters.k)
        println("multipliers ", repr(sol.original.multipliers))
    catch err
        println("error ", sprint(showerror, err))
    end
    flush(stdout)
end

verbose = MadNLP.INFO
hs = (hs071_f, hs071_c, [1.0, 5.0, 5.0, 1.0], [25.0, 40.0], [Inf, 40.0], fill(1.0, 4), fill(5.0, 4))
hs_far = (hs071_f, hs071_c, [5.0, 1.0, 1.0, 5.0], [25.0, 40.0], [Inf, 40.0], fill(1.0, 4), fill(5.0, 4))
circle = (circle_f, circle_c, [0.0, 0.0], [1.0, 2.0], [1.0, 2.0], [-Inf, -Inf], [Inf, Inf])
wb = (wb_f, wb_c, [-2.0, 1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [-Inf, 0.0, 0.0], [Inf, Inf, Inf])

# OptimizationMadNLP takes the print level from solve's verbose keyword
mk(; kw...) = OptimizationMadNLP.MadNLPOptimizer(; kw...)
ao(; kw...) = mk(additional_options = Dict{Symbol, Any}(kw...))

# the same with no constraints (bounds only)
function run_nocons(name, f, x0, lb, ub, opt)
    println("== ", name)
    flush(stdout)
    try
        of = Optimization.OptimizationFunction(f, ADTypes.AutoForwardDiff())
        prob = Optimization.OptimizationProblem(of, x0, nothing; lb = lb, ub = ub)
        sol = Optimization.solve(prob, opt; verbose = verbose)
        flush(stdout)
        println("retcode ", sol.retcode)
        println("status ", sol.original.status)
        println("u ", repr(sol.u))
        println("objective ", repr(sol.objective))
        println("iterations ", sol.original.counters.k)
    catch err
        println("error ", sprint(showerror, err))
    end
    flush(stdout)
end
rosen_f(x, p) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
eq_f(x, p) = (x[1] - 1)^2 + (x[2] - 2)^2
eq_c(res, x, p) = (res .= [x[1]^2 + x[2]^2])
eq = (eq_f, eq_c, [0.5, 0.5], [1.0], [1.0], [-Inf, -Inf], [Inf, Inf])
hs_fixed = (hs071_f, hs071_c, [1.0, 5.0, 5.0, 1.0], [25.0, 40.0], [Inf, 40.0], [1.0, 5.0, 1.0, 1.0], [5.0, 5.0, 5.0, 5.0])

run("hs071_default", hs..., mk())
run("hs071_far", hs_far..., mk())
run("infeasible_circle", circle..., mk())
run("wachter_biegler", wb..., mk())
run("hs071_lbfgs", hs..., mk(hessian_approximation = MadNLP.CompactLBFGS))
run("hs071_quality_function", hs..., mk(barrier = MadNLP.QualityFunctionUpdate()))
run("hs071_loqo", hs..., mk(barrier = MadNLP.LOQOUpdate()))
run("hs071_no_scaling", hs..., mk(nlp_scaling = false))
run("hs071_scaling_max_gradient", hs..., mk(nlp_scaling_max_gradient = 1.0))
run("hs071_acceptable", hs..., mk(acceptable_tol = 1e-2, acceptable_iter = 2))
run("hs071_lapack", hs..., mk(linear_solver = MadNLP.LapackCPUSolver))
run("hs071_sparse_kkt", hs..., mk(kkt_system = MadNLP.SparseKKTSystem))
run("hs071_condensed_kkt", hs..., mk(kkt_system = MadNLP.SparseCondensedKKTSystem))
run("hs071_far_lbfgs", hs_far..., mk(hessian_approximation = MadNLP.CompactLBFGS))
run("hs071_far_quality_function", hs_far..., mk(barrier = MadNLP.QualityFunctionUpdate()))
run("hs071_far_loqo", hs_far..., mk(barrier = MadNLP.LOQOUpdate()))
run("hs071_monotone_mu_init", hs..., mk(barrier = MadNLP.MonotoneUpdate(mu_init = 0.5)))
run("hs071_fixed", hs_fixed..., mk())
run("hs071_fixed_relax", hs_fixed..., ao(fixed_variable_treatment = MadNLP.RelaxBound))
run("hs071_relax_equality", hs..., ao(equality_treatment = MadNLP.RelaxEquality))
run("hs071_inertia_free", hs..., ao(inertia_correction_method = MadNLP.InertiaFree))
run("hs071_inertia_ignore", hs..., ao(inertia_correction_method = MadNLP.InertiaIgnore))
run("hs071_umfpack", hs..., mk(linear_solver = MadNLP.UmfpackSolver))
run("hs071_bound_push", hs..., ao(bound_push = 0.1, bound_fac = 0.2))
run("hs071_dual_zero", hs..., ao(dual_initialization_method = MadNLP.DualInitializeSetZero))
run("hs071_far_dual_zero", hs_far..., ao(dual_initialization_method = MadNLP.DualInitializeSetZero))
run("hs071_kappa_d", hs..., ao(kappa_d = 1e-2, bound_relax_factor = 0.0))
run("equality_only", eq..., mk())
run("equality_only_lbfgs", eq..., mk(hessian_approximation = MadNLP.CompactLBFGS))
run("infeasible_circle_lbfgs", circle..., mk(hessian_approximation = MadNLP.CompactLBFGS))
run("wachter_biegler_lbfgs", wb..., mk(hessian_approximation = MadNLP.CompactLBFGS))
run("wachter_biegler_inertia_free", wb..., ao(inertia_correction_method = MadNLP.InertiaFree))
run_nocons("rosenbrock_bounds", rosen_f, [-1.2, 1.0], [-1.5, -0.5], [0.5, 2.0], mk())
run_nocons("rosenbrock_bounds_lbfgs", rosen_f, [-1.2, 1.0], [-1.5, -0.5], [0.5, 2.0],
    mk(hessian_approximation = MadNLP.CompactLBFGS))


# CHOLMODSolver (Cholesky and LDL' factorizations), the dense KKT systems
# with LAPACK, the BFGS and DampedBFGS approximations (which need a dense
# KKT system) and the condensed KKT system on more problems
cholmod(; kw...) = mk(; linear_solver = MadNLP.CHOLMODSolver, kw...)
condensed(; kw...) = mk(; kkt_system = MadNLP.SparseCondensedKKTSystem, kw...)
dense(; kw...) = mk(; kkt_system = MadNLP.DenseKKTSystem, linear_solver = MadNLP.LapackCPUSolver, kw...)
dense_condensed(; kw...) = mk(; kkt_system = MadNLP.DenseCondensedKKTSystem, linear_solver = MadNLP.LapackCPUSolver,
    kw...)
run("hs071_cholmod_condensed", hs..., cholmod(kkt_system = MadNLP.SparseCondensedKKTSystem))
run("hs071_cholmod_ldl_condensed", hs..., cholmod(kkt_system = MadNLP.SparseCondensedKKTSystem,
    additional_options = Dict{Symbol, Any}(:cholmod_algorithm => MadNLP.LDL)))
run("hs071_far_cholmod_condensed", hs_far..., cholmod(kkt_system = MadNLP.SparseCondensedKKTSystem))
run("hs071_cholmod", hs..., cholmod())
run("hs071_far_condensed", hs_far..., condensed())
run("hs071_fixed_condensed", hs_fixed..., condensed())
run("infeasible_circle_condensed", circle..., condensed())
run("wachter_biegler_condensed", wb..., condensed())
run("equality_only_condensed", eq..., condensed())
run_nocons("rosenbrock_bounds_condensed", rosen_f, [-1.2, 1.0], [-1.5, -0.5], [0.5, 2.0], condensed())
run("hs071_dense", hs..., dense())
run("hs071_far_dense", hs_far..., dense())
run("hs071_dense_condensed", hs..., dense_condensed())
run("hs071_far_dense_condensed", hs_far..., dense_condensed())
run("wachter_biegler_dense", wb..., dense())
run("hs071_bfgs", hs..., dense(hessian_approximation = MadNLP.BFGS))
run("hs071_damped_bfgs", hs..., dense(hessian_approximation = MadNLP.DampedBFGS))
run("hs071_far_bfgs", hs_far..., dense(hessian_approximation = MadNLP.BFGS))
run("hs071_far_damped_bfgs", hs_far..., dense(hessian_approximation = MadNLP.DampedBFGS))
run("hs071_bfgs_dense_condensed", hs..., dense_condensed(hessian_approximation = MadNLP.BFGS))
run("hs071_damped_bfgs_dense_condensed", hs..., dense_condensed(hessian_approximation = MadNLP.DampedBFGS))
run("equality_only_bfgs", eq..., dense(hessian_approximation = MadNLP.BFGS))
run("wachter_biegler_bfgs", wb..., dense(hessian_approximation = MadNLP.BFGS))
run_nocons("rosenbrock_bounds_bfgs", rosen_f, [-1.2, 1.0], [-1.5, -0.5], [0.5, 2.0],
    dense(hessian_approximation = MadNLP.BFGS))
run_nocons("rosenbrock_bounds_damped_bfgs", rosen_f, [-1.2, 1.0], [-1.5, -0.5], [0.5, 2.0],
    dense(hessian_approximation = MadNLP.DampedBFGS))
