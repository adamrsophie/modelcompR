# Reference values for the options of IpoptOptimizer (R/optimization-ipopt.R):
# HS071 solved by Ipopt.jl with each option that changes the iterates, at
# print_level 5 so that the iteration tables can be compared line by line.
# Usage (packages: Ipopt):
#   julia --project=tools/julia tools/ipopt-options-probe.jl > ipopt-options.txt
using Ipopt

function hs071(opts; x0 = [1.0, 5.0, 5.0, 1.0])
    n, m = 4, 2
    f(x) = x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
    g!(x, gv) = (gv[1] = x[1] * x[2] * x[3] * x[4]; gv[2] = sum(x .^ 2); nothing)
    function gradf!(x, gv)
        gv[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3]); gv[2] = x[1] * x[4]
        gv[3] = x[1] * x[4] + 1; gv[4] = x[1] * (x[1] + x[2] + x[3]); nothing
    end
    function jac!(x, rows, cols, v)
        if v === nothing
            k = 0
            for i in 1:2, j in 1:4; k += 1; rows[k] = i; cols[k] = j; end
        else
            v[1] = x[2] * x[3] * x[4]; v[2] = x[1] * x[3] * x[4]; v[3] = x[1] * x[2] * x[4]; v[4] = x[1] * x[2] * x[3]
            v[5] = 2x[1]; v[6] = 2x[2]; v[7] = 2x[3]; v[8] = 2x[4]
        end
        nothing
    end
    function hess!(x, rows, cols, of, l, v)
        if v === nothing
            k = 0
            for i in 1:4, j in 1:i; k += 1; rows[k] = i; cols[k] = j; end
        else
            H = zeros(4, 4)
            H[1, 1] = of * 2x[4]; H[2, 1] = of * x[4]; H[3, 1] = of * x[4]; H[4, 1] = of * (2x[1] + x[2] + x[3])
            H[4, 2] = of * x[1]; H[4, 3] = of * x[1]
            H[2, 1] += l[1] * x[3] * x[4]; H[3, 1] += l[1] * x[2] * x[4]; H[4, 1] += l[1] * x[2] * x[3]
            H[3, 2] += l[1] * x[1] * x[4]; H[4, 2] += l[1] * x[1] * x[3]; H[4, 3] += l[1] * x[1] * x[2]
            for i in 1:4; H[i, i] += 2l[2]; end
            k = 0
            for i in 1:4, j in 1:i; k += 1; v[k] = H[i, j]; end
        end
        nothing
    end
    prob = Ipopt.CreateIpoptProblem(n, fill(1.0, 4), fill(5.0, 4), m, [25.0, 40.0], [2.0e19, 40.0], 8, 10,
                                    f, g!, gradf!, jac!, hess!)
    Ipopt.AddIpoptIntOption(prob, "print_level", 5)
    for (k, v) in opts
        v isa Int ? Ipopt.AddIpoptIntOption(prob, k, v) :
        v isa Float64 ? Ipopt.AddIpoptNumOption(prob, k, v) : Ipopt.AddIpoptStrOption(prob, k, v)
    end
    prob.x .= x0
    status = Ipopt.IpoptSolve(prob)
    return prob, status
end

cases = [
    ("default", Pair[]),
    ("far_start", Pair[]),
    ("lbfgs", ["hessian_approximation" => "limited-memory"]),
    ("lsr1", ["hessian_approximation" => "limited-memory", "limited_memory_update_type" => "sr1"]),
    ("lbfgs_hist3", ["hessian_approximation" => "limited-memory", "limited_memory_max_history" => 3]),
    ("warm_start", ["warm_start_init_point" => "yes"]),
    ("mu_init", ["mu_init" => 1e-3]),
    ("scaling_none", ["nlp_scaling_method" => "none"]),
    ("scaling_max_grad", ["nlp_scaling_max_gradient" => 1.0]),
    ("accept_every", ["accept_every_trial_step" => "yes"]),
    ("expect_infeasible", ["expect_infeasible_problem" => "yes"]),
    ("honor_bounds", ["honor_original_bounds" => "yes"]),
    ("adaptive_qf", ["mu_strategy" => "adaptive"]),
    ("adaptive_probing", ["mu_strategy" => "adaptive", "mu_oracle" => "probing"]),
    ("adaptive_loqo", ["mu_strategy" => "adaptive", "mu_oracle" => "loqo"]),
    ("adaptive_kkt_error", ["mu_strategy" => "adaptive", "adaptive_mu_globalization" => "kkt-error"]),
    ("adaptive_never", ["mu_strategy" => "adaptive", "adaptive_mu_globalization" => "never-monotone-mode"]),
    ("penalty_ls", ["line_search_method" => "penalty"]),
    ("cg_penalty_ls", ["line_search_method" => "cg-penalty"]),
    ("acceptable", ["acceptable_tol" => 1e-2, "acceptable_iter" => 2, "tol" => 1e-14]),
    ("constr_viol_tol", ["constr_viol_tol" => 1e-12, "tol" => 1e-3]),
    ("dual_inf_tol", ["dual_inf_tol" => 1e-12, "tol" => 1e-3]),
    ("compl_inf_tol", ["compl_inf_tol" => 1e-12, "tol" => 1e-3]),
    ("ma27", ["linear_solver" => "ma27"]),
    ("spral", ["linear_solver" => "spral"]),
    ("equilibration", ["nlp_scaling_method" => "equilibration-based"]),
    ("mc19", ["linear_system_scaling" => "mc19"]),
]
for (name, opts) in cases
    println("== ", name)
    flush(stdout)
    x0 = name == "far_start" ? [5.0, 1.0, 1.0, 5.0] : [1.0, 5.0, 5.0, 1.0]
    try
        prob, status = hs071(opts; x0 = x0)
        flush(stdout)
        println("status ", status)
        println("x ", repr(prob.x))
        println("objective ", repr(prob.obj_val))
        println("mult_g ", repr(prob.mult_g))
    catch err
        println("error ", err)
    end
    flush(stdout)
end
