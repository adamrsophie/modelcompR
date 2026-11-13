# Reference values for the NLopt interface (R/optimization-nlopt.R and the
# native ports of the NLopt algorithms): runs NLopt.jl (the NLopt C library)
# and prints R assignments with full precision. Every run is seeded with
# NLopt.srand(1234) so the stochastic algorithms are reproducible. Usage:
#   julia --project=<env with NLopt> tools/nlopt-probe.jl > refs.R
using NLopt

rv(v) = "c(" * join([repr(Float64(t)) for t in v], ", ") * ")"
emit(name, minf, minx, ret, opt) = println("ref\$$name <- list(x = $(rv(minx)), f = $(repr(minf)), ret = \"$(ret)\", ",
                                           "numevals = $(opt.numevals)L)")
emit_err(name, e) = println("ref\$$name <- list(error = $(repr(sprint(showerror, e))))")

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
function ball(x, g)
    if length(g) > 0
        g .= 2 .* x
    end
    sum(abs2, x) - 4
end
function plane(x, g)
    if length(g) > 0
        g .= 1.0
    end
    sum(x) - 1
end

# the earlier references, kept
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

# stopval and maxeval stops
opt = Opt(:LN_NELDERMEAD, 3); opt.min_objective = poly3; opt.stopval = 1.0
(minf, minx, ret) = optimize(opt, [0.0, 0.0, 0.0]); emit("nm_poly3_stopval", minf, minx, ret, opt)
opt = Opt(:LN_NELDERMEAD, 3); opt.max_objective = (x, g) -> -poly3(x, g); opt.stopval = -1.0
(minf, minx, ret) = optimize(opt, [0.0, 0.0, 0.0]); emit("nm_poly3_max_stopval", minf, minx, ret, opt)
opt = Opt(:LN_NELDERMEAD, 3); opt.min_objective = poly3; opt.maxeval = 25
(minf, minx, ret) = optimize(opt, [0.0, 0.0, 0.0]); emit("nm_poly3_maxeval", minf, minx, ret, opt)

# every algorithm on three problems: poly3 in [-5, 5]^3 from 0, the
# Rosenbrock function in [-2, 2]^2 from (-1.2, 1), and poly3 with the
# constraints |x|^2 <= 4 (and sum(x) = 1 where equalities are allowed)
const ALL = [:GN_DIRECT, :GN_DIRECT_L, :GN_DIRECT_L_RAND, :GN_DIRECT_NOSCAL, :GN_DIRECT_L_NOSCAL,
    :GN_DIRECT_L_RAND_NOSCAL, :GN_ORIG_DIRECT, :GN_ORIG_DIRECT_L, :GD_STOGO, :GD_STOGO_RAND,
    :LD_LBFGS_NOCEDAL, :LD_LBFGS, :LN_PRAXIS, :LD_VAR1, :LD_VAR2, :LD_TNEWTON, :LD_TNEWTON_RESTART,
    :LD_TNEWTON_PRECOND, :LD_TNEWTON_PRECOND_RESTART, :GN_CRS2_LM, :GN_MLSL, :GD_MLSL, :GN_MLSL_LDS,
    :GD_MLSL_LDS, :LD_MMA, :LN_COBYLA, :LN_NEWUOA, :LN_NEWUOA_BOUND, :LN_NELDERMEAD, :LN_SBPLX,
    :LN_AUGLAG, :LD_AUGLAG, :LN_AUGLAG_EQ, :LD_AUGLAG_EQ, :LN_BOBYQA, :GN_ISRES, :AUGLAG, :AUGLAG_EQ,
    :G_MLSL, :G_MLSL_LDS, :LD_SLSQP, :LD_CCSAQ, :GN_ESCH, :GN_AGS]
const EQ = Set([:LN_COBYLA, :GN_ISRES, :LN_AUGLAG, :LD_AUGLAG, :LN_AUGLAG_EQ, :LD_AUGLAG_EQ, :AUGLAG, :AUGLAG_EQ,
    :LD_SLSQP])
const INEQ = union(EQ, Set([:GN_ORIG_DIRECT, :GN_ORIG_DIRECT_L, :LD_MMA, :LD_CCSAQ, :GN_AGS]))
isglobal(a) = startswith(string(a), "G")
needslocal(a) = a in (:AUGLAG, :AUGLAG_EQ, :G_MLSL, :G_MLSL_LDS)

function setup(a, n, f, lb, ub)
    opt = Opt(a, n)
    opt.min_objective = f
    opt.lower_bounds = lb; opt.upper_bounds = ub
    opt.xtol_rel = 1e-8
    opt.maxeval = isglobal(a) ? 3000 : 5000
    if needslocal(a)
        lo = Opt(:LD_LBFGS, n); lo.xtol_rel = 1e-8; lo.maxeval = 1000
        opt.local_optimizer = lo
    end
    opt
end

function run(name, a, n, f, lb, ub, x0; cons = :none)
    NLopt.srand(1234)
    try
        opt = setup(a, n, f, lb, ub)
        if cons != :none
            inequality_constraint!(opt, ball, 1e-8)
            cons == :eq && equality_constraint!(opt, plane, 1e-8)
        end
        (minf, minx, ret) = optimize(opt, copy(x0))
        emit(name, minf, minx, ret, opt)
    catch e
        emit_err(name, e)
    end
end

for a in ALL
    s = lowercase(string(a))
    run("$(s)_poly3", a, 3, poly3, fill(-5.0, 3), fill(5.0, 3), [0.0, 0.0, 0.0])
    run("$(s)_rosen", a, 2, rosen, [-2.0, -2.0], [2.0, 2.0], [-1.2, 1.0])
    if a in INEQ
        run("$(s)_cons", a, 3, poly3, fill(-5.0, 3), fill(5.0, 3), [0.0, 0.0, 0.0];
            cons = a in EQ ? :eq : :ineq)
    end
end

# the Mersenne Twister of NLopt (util/mt19937ar.c) through the seeded
# uniform numbers of the stochastic algorithms is checked by their runs; a
# fixed variable (lb == ub) exercises the dimension elimination
run("ln_neldermead_fixed", :LN_NELDERMEAD, 3, poly3, [-5.0, 2.0, -5.0], [5.0, 2.0, 5.0], [0.0, 2.0, 0.0])
run("gn_crs2_lm_fixed", :GN_CRS2_LM, 3, poly3, [-5.0, 2.0, -5.0], [5.0, 2.0, 5.0], [0.0, 2.0, 0.0])
