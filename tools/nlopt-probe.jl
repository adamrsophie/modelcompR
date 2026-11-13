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

# ---------------------------------------------------------------------------
# Bit level references of the Luksan (LD_LBFGS, LD_VAR1, LD_VAR2, the
# LD_TNEWTON family), COBYLA, NEWUOA, BOBYQA and SLSQP ports: the minimum,
# the point and every evaluated objective value in hexadecimal (%a, which R
# reads back as a numeric literal). The objective trace goes to the file
# named by TRACE (one line per case: the name, then the values, tab
# separated), the summary to standard output as R assignments of `hex`.
using Printf
hx(t) = isnan(t) ? "NaN" : isinf(t) ? (t > 0 ? "Inf" : "-Inf") : @sprintf("%a", t)
hv(v) = "c(" * join(["\"$(hx(Float64(t)))\"" for t in v], ", ") * ")"
const TRACEIO = open(get(ENV, "TRACE", "nlopt-trace.tsv"), "w")
println("# NLopt ", NLopt.version())

function quart5(x, g)
    s = sum(x)
    if length(g) > 0
        for i in 1:5
            g[i] = 2 * i * (x[i] - 0.1 * i) + 0.4 * s * s * s
        end
    end
    v = 0.1 * s * s * s * s
    for i in 1:5
        d = x[i] - 0.1 * i
        v += i * d * d
    end
    v
end

function runhex(name, a, n, f, lb, ub, x0; cons = :none, xtol_rel = 0.0, ftol_rel = 0.0, maxeval = 0,
                params = Dict{String,Float64}(), vector_storage = 0)
    tr = Float64[]
    fw = (x, g) -> (v = f(x, g); push!(tr, v); v)
    try
        opt = Opt(a, n)
        opt.min_objective = fw
        lb === nothing || (opt.lower_bounds = lb)
        ub === nothing || (opt.upper_bounds = ub)
        opt.xtol_rel = xtol_rel; opt.ftol_rel = ftol_rel; opt.maxeval = maxeval
        vector_storage > 0 && (opt.vector_storage = vector_storage)
        for (k, v) in params
            opt.params[k] = v
        end
        if cons != :none
            inequality_constraint!(opt, ball, 1e-8)
            cons == :eq && equality_constraint!(opt, plane, 1e-8)
        end
        (minf, minx, ret) = optimize(opt, copy(x0))
        println("hex\$$name <- list(x = $(hv(minx)), f = \"$(hx(minf))\", ret = \"$(ret)\", numevals = $(opt.numevals)L)")
    catch e
        println("hex\$$name <- list(error = $(repr(sprint(showerror, e))))")
    end
    println(TRACEIO, name, "\t", join(hx.(tr), "\t"))
end

const BOX3 = (fill(-5.0, 3), fill(5.0, 3))
for a in (:LD_LBFGS, :LD_VAR1, :LD_VAR2, :LD_TNEWTON, :LD_TNEWTON_RESTART, :LD_TNEWTON_PRECOND,
          :LD_TNEWTON_PRECOND_RESTART, :LN_COBYLA, :LN_NEWUOA, :LN_NEWUOA_BOUND, :LN_BOBYQA, :LD_SLSQP)
    s = lowercase(string(a))
    nd = a in (:LN_NEWUOA,)
    runhex("$(s)_poly3", a, 3, poly3, nothing, nothing, [0.0, 0.0, 0.0]; xtol_rel = 1e-10, maxeval = 2000)
    runhex("$(s)_poly3_box", a, 3, poly3, BOX3..., [0.0, 0.0, 0.0]; xtol_rel = 1e-8, maxeval = 2000)
    runhex("$(s)_rosen", a, 2, rosen, nothing, nothing, [-1.2, 1.0]; ftol_rel = 1e-14, maxeval = 3000)
    runhex("$(s)_rosen_box", a, 2, rosen, [-2.0, -2.0], [0.8, 2.0], [-1.2, 1.0]; ftol_rel = 1e-14, maxeval = 3000)
    runhex("$(s)_rosen_max", a, 2, rosen, [-2.0, -2.0], [2.0, 2.0], [-1.2, 1.0]; maxeval = 25)
    runhex("$(s)_quart5", a, 5, quart5, fill(-3.0, 5), fill(3.0, 5), fill(1.0, 5); xtol_rel = 1e-9, maxeval = 4000)
    if a in (:LN_COBYLA, :LD_SLSQP)
        runhex("$(s)_cons", a, 3, poly3, BOX3..., [0.0, 0.0, 0.0]; cons = :eq, xtol_rel = 1e-10, maxeval = 5000)
        runhex("$(s)_cons_ineq", a, 3, poly3, BOX3..., [0.0, 0.0, 0.0]; cons = :ineq, xtol_rel = 1e-10,
               maxeval = 5000)
    end
end
# the Luksan parameters: a small vector storage and a gradient tolerance
runhex("ld_lbfgs_rosen_m3", :LD_LBFGS, 2, rosen, nothing, nothing, [-1.2, 1.0]; ftol_rel = 1e-14, maxeval = 3000,
       vector_storage = 3)
runhex("ld_var2_quart5_tolg", :LD_VAR2, 5, quart5, nothing, nothing, fill(1.0, 5); maxeval = 3000,
       params = Dict("tolg" => 1e-4))
close(TRACEIO)

# ---------------------------------------------------------------------------
# Bit level references of the DIRECT family (GN_DIRECT*, GN_ORIG_DIRECT*),
# LN_PRAXIS, StoGO (GD_STOGO, GD_STOGO_RAND) and AGS (GN_AGS) ports: every
# run is seeded with NLopt.srand(1234); the summary goes to standard output
# as R assignments of `ghex`, the objective trace to the file named by
# TRACE_GLOBAL (one line per case: the name, then the values, tab separated).
const GTRACEIO = open(get(ENV, "TRACE_GLOBAL", "probe-out/nlopt-trace-global.tsv"), "w")

function runghex(name, a, n, f, lb, ub, x0; cons = :none, xtol_rel = 0.0, ftol_rel = 0.0, maxeval = 0,
                 params = Dict{String,Float64}())
    tr = Float64[]
    fw = (x, g) -> (v = f(x, g); push!(tr, v); v)
    NLopt.srand(1234)
    try
        opt = Opt(a, n)
        opt.min_objective = fw
        lb === nothing || (opt.lower_bounds = lb)
        ub === nothing || (opt.upper_bounds = ub)
        opt.xtol_rel = xtol_rel; opt.ftol_rel = ftol_rel; opt.maxeval = maxeval
        for (k, v) in params
            opt.params[k] = v
        end
        cons == :ineq && inequality_constraint!(opt, ball, 1e-8)
        (minf, minx, ret) = optimize(opt, copy(x0))
        println("ghex\$$name <- list(x = $(hv(minx)), f = \"$(hx(minf))\", ret = \"$(ret)\", numevals = $(opt.numevals)L)")
    catch e
        println("ghex\$$name <- list(error = $(repr(sprint(showerror, e))))")
    end
    println(GTRACEIO, name, "\t", join(hx.(tr), "\t"))
end

const GBOX3 = (fill(-5.0, 3), fill(5.0, 3))
const GBOX2 = ([-2.0, -2.0], [2.0, 2.0])
for a in (:GN_DIRECT, :GN_DIRECT_L, :GN_DIRECT_L_RAND, :GN_DIRECT_NOSCAL, :GN_DIRECT_L_NOSCAL,
          :GN_DIRECT_L_RAND_NOSCAL, :GN_ORIG_DIRECT, :GN_ORIG_DIRECT_L)
    s = lowercase(string(a))
    runghex("$(s)_poly3", a, 3, poly3, GBOX3..., [0.0, 0.0, 0.0]; xtol_rel = 1e-8, maxeval = 600)
    runghex("$(s)_rosen", a, 2, rosen, GBOX2..., [-1.2, 1.0]; xtol_rel = 1e-8, maxeval = 500)
    runghex("$(s)_quart5", a, 5, quart5, fill(-3.0, 5), fill(2.0, 5), fill(1.0, 5); maxeval = 800)
    runghex("$(s)_xtol", a, 2, rosen, GBOX2..., [-1.2, 1.0]; xtol_rel = 1e-3, maxeval = 20000)
end
runghex("gn_direct_magic", :GN_DIRECT, 2, rosen, GBOX2..., [-1.2, 1.0]; maxeval = 400,
        params = Dict("magic_eps" => 1e-4))
runghex("gn_direct_l_fixed", :GN_DIRECT_L, 3, poly3, [-5.0, 2.0, -5.0], [5.0, 2.0, 5.0], [0.0, 2.0, 0.0];
        maxeval = 400)
runghex("gn_direct_l_ftol", :GN_DIRECT_L, 2, rosen, GBOX2..., [-1.2, 1.0]; ftol_rel = 1e-3, maxeval = 5000)
for a in (:GN_ORIG_DIRECT, :GN_ORIG_DIRECT_L)
    s = lowercase(string(a))
    runghex("$(s)_cons", a, 3, poly3, GBOX3..., [0.0, 0.0, 0.0]; cons = :ineq, maxeval = 600)
end
runghex("gn_orig_direct_magic", :GN_ORIG_DIRECT, 2, rosen, GBOX2..., [-1.2, 1.0]; maxeval = 400,
        params = Dict("magic_eps" => 1e-4, "magic_eps_abs" => 1e-6))

runghex("ln_praxis_poly3", :LN_PRAXIS, 3, poly3, nothing, nothing, [0.0, 0.0, 0.0]; xtol_rel = 1e-10, maxeval = 2000)
runghex("ln_praxis_poly3_box", :LN_PRAXIS, 3, poly3, GBOX3..., [0.0, 0.0, 0.0]; xtol_rel = 1e-8, maxeval = 2000)
runghex("ln_praxis_poly3_half", :LN_PRAXIS, 3, poly3, fill(-1.0, 3), fill(Inf, 3), [0.0, 0.0, 0.0];
        xtol_rel = 1e-8, maxeval = 2000)
runghex("ln_praxis_poly3_upper", :LN_PRAXIS, 3, poly3, fill(-Inf, 3), fill(1.5, 3), [0.0, 0.0, 0.0];
        xtol_rel = 1e-8, maxeval = 2000)
runghex("ln_praxis_rosen", :LN_PRAXIS, 2, rosen, nothing, nothing, [-1.2, 1.0]; ftol_rel = 1e-14, maxeval = 3000)
runghex("ln_praxis_rosen_box", :LN_PRAXIS, 2, rosen, [-2.0, -2.0], [0.8, 2.0], [-1.2, 1.0]; ftol_rel = 1e-14,
        maxeval = 3000)
runghex("ln_praxis_rosen_max", :LN_PRAXIS, 2, rosen, nothing, nothing, [-1.2, 1.0]; maxeval = 25)
runghex("ln_praxis_quart5", :LN_PRAXIS, 5, quart5, fill(-3.0, 5), fill(3.0, 5), fill(1.0, 5); xtol_rel = 1e-9,
        maxeval = 4000)
runghex("ln_praxis_t0", :LN_PRAXIS, 2, rosen, nothing, nothing, [-1.2, 1.0]; maxeval = 3000,
        params = Dict("t0_tol" => 1e-6))

for a in (:GD_STOGO, :GD_STOGO_RAND)
    s = lowercase(string(a))
    runghex("$(s)_poly3", a, 3, poly3, GBOX3..., [0.0, 0.0, 0.0]; xtol_rel = 1e-8, maxeval = 600)
    runghex("$(s)_rosen", a, 2, rosen, GBOX2..., [-1.2, 1.0]; xtol_rel = 1e-8, maxeval = 600)
    runghex("$(s)_quart5", a, 5, quart5, fill(-3.0, 5), fill(2.0, 5), fill(1.0, 5); maxeval = 800)
end

runghex("gn_ags_poly3", :GN_AGS, 3, poly3, GBOX3..., [0.0, 0.0, 0.0]; maxeval = 600)
runghex("gn_ags_rosen", :GN_AGS, 2, rosen, GBOX2..., [-1.2, 1.0]; maxeval = 500)
runghex("gn_ags_rosen_xtol", :GN_AGS, 2, rosen, GBOX2..., [-1.2, 1.0]; xtol_rel = 1e-3, maxeval = 5000)
runghex("gn_ags_cons", :GN_AGS, 3, poly3, GBOX3..., [0.0, 0.0, 0.0]; cons = :ineq, maxeval = 600)
close(GTRACEIO)
