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
        sol = Optimization.solve(prob, opt)
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

verbose = MadNLP.INFO
hs = (hs071_f, hs071_c, [1.0, 5.0, 5.0, 1.0], [25.0, 40.0], [Inf, 40.0], fill(1.0, 4), fill(5.0, 4))
hs_far = (hs071_f, hs071_c, [5.0, 1.0, 1.0, 5.0], [25.0, 40.0], [Inf, 40.0], fill(1.0, 4), fill(5.0, 4))
circle = (circle_f, circle_c, [0.0, 0.0], [1.0, 2.0], [1.0, 2.0], [-Inf, -Inf], [Inf, Inf])
wb = (wb_f, wb_c, [-2.0, 1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [-Inf, 0.0, 0.0], [Inf, Inf, Inf])

mk(; kw...) = OptimizationMadNLP.MadNLPOptimizer(; print_level = verbose, kw...)

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
