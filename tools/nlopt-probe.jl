# Reference values for the NLopt interface (R/optimization-nlopt.R): runs
# NLopt.jl and prints R assignments with full precision. Usage:
#   julia --project=<env with NLopt> tools/nlopt-probe.jl > refs.R
using NLopt

rv(v) = "c(" * join([repr(Float64(t)) for t in v], ", ") * ")"
emit(name, minf, minx, ret, opt) = println("ref\$$name <- list(x = $(rv(minx)), f = $(repr(minf)), ret = \"$(ret)\", ",
                                           "numevals = $(opt.numevals)L)")

function poly3(x, g)
    a = x[1] - 1; b = x[2] - 2; c = x[3] + 1
    if length(g) > 0
        g[1] = 2a + 2a * b * b + 0.5c; g[2] = 4b + 2a * a * b; g[3] = 6c + 0.5a
    end
    a * a + 2b * b + 3c * c + a * a * b * b + 0.5a * c
end
function rosen(x, g)
    a = 1 - x[1]; b = x[2] - x[1] * x[1]
    if length(g) > 0
        g[1] = -2a - 400x[1] * b; g[2] = 200b
    end
    a * a + 100b * b
end

opt = Opt(:LN_NELDERMEAD, 3); opt.min_objective = poly3; opt.xtol_rel = 1e-10; opt.maxeval = 2000
(minf, minx, ret) = optimize(opt, [0.0, 0.0, 0.0]); emit("nm_poly3", minf, minx, ret, opt)
opt = Opt(:LD_LBFGS, 2); opt.min_objective = rosen; opt.ftol_rel = 1e-14
(minf, minx, ret) = optimize(opt, [-1.2, 1.0]); emit("lbfgs_rosen", minf, minx, ret, opt)
opt = Opt(:LD_LBFGS, 2); opt.min_objective = rosen; opt.lower_bounds = [-2.0, -2.0]; opt.upper_bounds = [0.8, 2.0]
opt.ftol_rel = 1e-14
(minf, minx, ret) = optimize(opt, [-1.2, 1.0]); emit("lbfgs_rosen_bounds", minf, minx, ret, opt)
opt = Opt(:LN_COBYLA, 3); opt.min_objective = poly3; opt.xtol_rel = 1e-10; opt.maxeval = 5000
inequality_constraint!(opt, (x, g) -> x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 4, 1e-8)
equality_constraint!(opt, (x, g) -> x[1] + x[2] + x[3] - 1, 1e-8)
(minf, minx, ret) = optimize(opt, [0.0, 0.0, 0.0]); emit("cobyla_poly3", minf, minx, ret, opt)
