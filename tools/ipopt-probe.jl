# Reference values for IpoptOptimizer() (R/optimization-ipopt.R): solves
# HS071 with Ipopt.jl and prints
# R assignments with full precision. Usage:
#   julia --project=<env with Ipopt> tools/ipopt-probe.jl > refs.R
using Ipopt

rv(v) = "c(" * join([repr(Float64(t)) for t in v], ", ") * ")"

function solve_ipopt(n, xl, xu, m, gl, gu, f, g!, gradf!, jac!, hess!, x0)
    jrows = [i for i in 1:m for j in 1:n]
    jcols = [j for i in 1:m for j in 1:n]
    hrows = [i for i in 1:n for j in 1:i]
    hcols = [j for i in 1:n for j in 1:i]
    eval_jac_g(x, rows, cols, values) = values === nothing ? (rows .= jrows; cols .= jcols; nothing) : (jac!(values, x); nothing)
    eval_h(x, rows, cols, of, lam, values) = values === nothing ? (rows .= hrows; cols .= hcols; nothing) : (hess!(values, x, of, lam); nothing)
    prob = Ipopt.CreateIpoptProblem(n, xl, xu, m, gl, gu, m * n, n * (n + 1) ÷ 2, f, (x, gv) -> (g!(gv, x); nothing),
                                    (x, gv) -> (gradf!(gv, x); nothing), eval_jac_g, eval_h)
    Ipopt.AddIpoptIntOption(prob, "print_level", 0)
    prob.x .= x0
    status = Ipopt.IpoptSolve(prob)
    return prob, status
end

# HS071
f(x) = x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
g!(gv, x) = (gv[1] = x[1] * x[2] * x[3] * x[4]; gv[2] = sum(x .^ 2))
gradf!(gv, x) = (gv[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3]); gv[2] = x[1] * x[4]; gv[3] = x[1] * x[4] + 1;
                 gv[4] = x[1] * (x[1] + x[2] + x[3]))
jac!(v, x) = (v[1] = x[2] * x[3] * x[4]; v[2] = x[1] * x[3] * x[4]; v[3] = x[1] * x[2] * x[4]; v[4] = x[1] * x[2] * x[3];
              v[5] = 2x[1]; v[6] = 2x[2]; v[7] = 2x[3]; v[8] = 2x[4])
function hess!(v, x, of, l)
    H = zeros(4, 4)
    H[1, 1] = of * 2x[4]; H[2, 1] = of * x[4]; H[3, 1] = of * x[4]; H[4, 1] = of * (2x[1] + x[2] + x[3])
    H[4, 2] = of * x[1]; H[4, 3] = of * x[1]
    H[2, 1] += l[1] * x[3] * x[4]; H[3, 1] += l[1] * x[2] * x[4]; H[4, 1] += l[1] * x[2] * x[3]
    H[3, 2] += l[1] * x[1] * x[4]; H[4, 2] += l[1] * x[1] * x[3]; H[4, 3] += l[1] * x[1] * x[2]
    for i in 1:4; H[i, i] += 2l[2]; end
    k = 0
    for i in 1:4, j in 1:i; k += 1; v[k] = H[i, j]; end
end
prob, status = solve_ipopt(4, fill(1.0, 4), fill(5.0, 4), 2, [25.0, 40.0], [2.0e19, 40.0], f, g!, gradf!, jac!, hess!,
                           [1.0, 5.0, 5.0, 1.0])
println("ref\$hs071 <- list(x = $(rv(prob.x)), f = $(repr(prob.obj_val)), status = \"$(status)\", mult_g = $(rv(prob.mult_g)))")
