# Reference values for the PRIMA ports (R/optimization-prima*.R): runs
# PRIMA.jl (libprima) and prints R assignments with full precision. Usage:
#   julia --project=<env with PRIMA> tools/prima-probe.jl > refs.R
# The UOBYQA, NEWUOA, BOBYQA and LINCOA values are printed as hexadecimal
# literals (exact in R) for tests/testthat/test-optimization-prima.R.
using PRIMA, Printf

rv(v) = "c(" * join([repr(Float64(t)) for t in v], ", ") * ")"
function emit(name, x, info)
    println("ref\$$name <- list(x = $(rv(x)), fx = $(repr(info.fx)), nf = $(info.nf)L, status = \"$(info.status)\", ",
            "cstrv = $(repr(info.cstrv)))")
end

poly3(x) = (a = x[1] - 1; b = x[2] - 2; c = x[3] + 1; a * a + 2 * b * b + 3 * c * c + a * a * b * b + 0.5 * a * c)
rosen(x) = (a = 1 - x[1]; b = x[2] - x[1] * x[1]; a * a + 100 * b * b)

x, info = PRIMA.cobyla(poly3, [0.0, 0.0, 0.0]; nonlinear_ineq = x -> [x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 4])
emit("cobyla_poly3_ball", x, info)
x, info = PRIMA.cobyla(rosen, [-1.2, 1.0]; linear_ineq = ([1.0 1.0], [1.5]),
                       nonlinear_ineq = x -> [x[1] * x[1] + x[2] * x[2] - 1.5])
emit("cobyla_rosen_lin_disk", x, info)
x, info = PRIMA.cobyla(poly3, [0.0, 0.0, 0.0]; linear_eq = ([1.0 1.0 1.0], [1.0]), xl = [-1.0, -1.0, -2.0],
                       xu = [2.0, 1.5, 2.0], nonlinear_ineq = x -> [x[1] * x[3] - 0.5])
emit("cobyla_poly3_eq_bounds", x, info)
x, info = PRIMA.cobyla(rosen, [-1.2, 1.0]; nonlinear_ineq = x -> [-x[1]], rhobeg = 0.5, rhoend = 1e-8, maxfun = 400)
emit("cobyla_rosen_half", x, info)

hx(t) = @sprintf("%a", Float64(t))
hv(v) = "c(" * join(hx.(v), ", ") * ")"
function emith(name, x, info)
    println("ref\$$name <- list(x = $(hv(x)), fx = $(hx(info.fx)), nf = $(info.nf)L, status = $(Int(info.status))L)")
end
# objectives written as left to right loops, the same operations as in R
function chained(x)
    s = 0.0
    for i in 1:length(x)-1
        a = 1 - x[i]; b = x[i+1] - x[i] * x[i]
        s = s + (a * a + 4 * (b * b))
    end
    s
end

# NEWUOA
emith("newuoa_rosen", PRIMA.newuoa(rosen, [-1.2, 1.0])...)
emith("newuoa_poly3", PRIMA.newuoa(poly3, [0.0, 0.0, 0.0])...)
emith("newuoa_rosen_npt4", PRIMA.newuoa(rosen, [-1.2, 1.0]; npt = 4)...)
emith("newuoa_rosen_npt6", PRIMA.newuoa(rosen, [-1.2, 1.0]; npt = 6)...)
emith("newuoa_poly3_rho", PRIMA.newuoa(poly3, [0.0, 0.0, 0.0]; rhobeg = 0.5, rhoend = 1e-8)...)
emith("newuoa_chained5", PRIMA.newuoa(chained, fill(-1.0, 5))...)
emith("newuoa_chained5_npt15", PRIMA.newuoa(chained, fill(-1.0, 5); npt = 15)...)
emith("newuoa_chained5_npt7", PRIMA.newuoa(chained, fill(-1.0, 5); npt = 7)...)
emith("newuoa_rosen_maxfun", PRIMA.newuoa(rosen, [-1.2, 1.0]; maxfun = 30)...)
emith("newuoa_rosen_ftarget", PRIMA.newuoa(rosen, [-1.2, 1.0]; ftarget = 1e-3)...)
emith("newuoa_chained10_npt30", PRIMA.newuoa(chained, fill(-1.0, 10); npt = 30)...)
emith("newuoa_rosen_scale", PRIMA.newuoa(rosen, [-1.2, 1.0]; scale = [2.0, 0.5])...)
# edge cases: maxfun below n + 3 (npt revised), rhoend = 0, rhoend = rhobeg,
# ftarget met at x0, NaN and Inf values of the objective
nanrosen(x) = x[1] + x[2] > 1.5 ? NaN : rosen(x)
infrosen(x) = x[1] > 0.9 ? Inf : rosen(x)
emith("newuoa_rosen_maxfun3", PRIMA.newuoa(rosen, [-1.2, 1.0]; maxfun = 3)...)
emith("newuoa_rosen_rhoend0", PRIMA.newuoa(rosen, [-1.2, 1.0]; rhoend = 0.0)...)
emith("newuoa_rosen_rhoend1", PRIMA.newuoa(rosen, [-1.2, 1.0]; rhoend = 1.0)...)
emith("newuoa_rosen_ftarget0", PRIMA.newuoa(rosen, [-1.2, 1.0]; ftarget = 100.0)...)
emith("newuoa_nanrosen", PRIMA.newuoa(nanrosen, [-1.2, 1.0])...)
emith("newuoa_infrosen", PRIMA.newuoa(infrosen, [-1.2, 1.0])...)

# UOBYQA (npt is (n + 1)(n + 2)/2 and cannot be set)
one(x) = (x[1] - 2) * (x[1] - 2) + 0.1 * (x[1] * x[1] * x[1] * x[1])
emith("uobyqa_rosen", PRIMA.uobyqa(rosen, [-1.2, 1.0])...)
emith("uobyqa_poly3", PRIMA.uobyqa(poly3, [0.0, 0.0, 0.0])...)
emith("uobyqa_poly3_rho", PRIMA.uobyqa(poly3, [0.0, 0.0, 0.0]; rhobeg = 0.5, rhoend = 1e-8)...)
emith("uobyqa_chained5", PRIMA.uobyqa(chained, fill(-1.0, 5))...)
emith("uobyqa_chained7", PRIMA.uobyqa(chained, fill(-1.0, 7))...)
emith("uobyqa_one", PRIMA.uobyqa(one, [0.0])...)
# a gradient above 1e8 at the start (the trust region step scales the model)
emith("uobyqa_one_far", PRIMA.uobyqa(one, [1000.0])...)
emith("uobyqa_rosen_maxfun", PRIMA.uobyqa(rosen, [-1.2, 1.0]; maxfun = 30)...)
emith("uobyqa_rosen_ftarget", PRIMA.uobyqa(rosen, [-1.2, 1.0]; ftarget = 1e-3)...)
emith("uobyqa_rosen_scale", PRIMA.uobyqa(rosen, [-1.2, 1.0]; scale = [2.0, 0.5])...)
# edge cases: maxfun below (n + 1)(n + 2)/2 + 1, rhoend = 0, rhoend = rhobeg,
# ftarget met at x0, NaN and Inf values of the objective
emith("uobyqa_rosen_maxfun3", PRIMA.uobyqa(rosen, [-1.2, 1.0]; maxfun = 3)...)
emith("uobyqa_rosen_rhoend0", PRIMA.uobyqa(rosen, [-1.2, 1.0]; rhoend = 0.0)...)
emith("uobyqa_rosen_rhoend1", PRIMA.uobyqa(rosen, [-1.2, 1.0]; rhoend = 1.0)...)
emith("uobyqa_rosen_ftarget0", PRIMA.uobyqa(rosen, [-1.2, 1.0]; ftarget = 100.0)...)
emith("uobyqa_nanrosen", PRIMA.uobyqa(nanrosen, [-1.2, 1.0])...)
emith("uobyqa_infrosen", PRIMA.uobyqa(infrosen, [-1.2, 1.0])...)

# BOBYQA and LINCOA (hexadecimal literals, with the constraint violation)
function emithc(name, x, info)
    println("ref\$$name <- list(x = $(hv(x)), fx = $(hx(info.fx)), nf = $(info.nf)L, status = $(Int(info.status))L, ",
            "cstrv = $(hx(info.cstrv)))")
end
hm(A) = "matrix(" * hv(vec(Float64.(A))) * ", $(size(A, 1)), $(size(A, 2)))"

# BOBYQA: x0 is revised to the bounds or rhobeg inside them (honour_x0 is
# false since PRIMA.jl passes rhobeg), rhobeg above minimum(xu - xl)/2 is
# revised to minimum(xu - xl)/4
emithc("bobyqa_rosen", PRIMA.bobyqa(rosen, [-1.2, 1.0]; xl = [-2.0, -2.0], xu = [2.0, 2.0])...)
emithc("bobyqa_rosen_active", PRIMA.bobyqa(rosen, [-1.2, 1.0]; xl = [-2.0, -2.0], xu = [0.5, 2.0])...)
emithc("bobyqa_rosen_free", PRIMA.bobyqa(rosen, [-1.2, 1.0])...)
emithc("bobyqa_rosen_inside", PRIMA.bobyqa(rosen, [-0.5, 1.5]; xl = [-2.0, -2.0], xu = [0.8, 4.0])...)
emithc("bobyqa_rosen_outside", PRIMA.bobyqa(rosen, [5.0, -5.0]; xl = [-2.0, -2.0], xu = [2.0, 2.0])...)
emithc("bobyqa_poly3", PRIMA.bobyqa(poly3, [0.0, 0.0, 0.0]; xl = [-1.0, -1.0, -3.0], xu = [0.5, 3.0, 3.0])...)
emithc("bobyqa_poly3_rho", PRIMA.bobyqa(poly3, [0.0, 0.0, 0.0]; xl = [-3.0, -3.0, -3.0], xu = [3.0, 3.0, 3.0],
                                       rhobeg = 0.5, rhoend = 1e-8)...)
emithc("bobyqa_rosen_bigrho", PRIMA.bobyqa(rosen, [-1.2, 1.0]; xl = [-2.0, -2.0], xu = [2.0, 2.0], rhobeg = 3.0)...)
emithc("bobyqa_chained5", PRIMA.bobyqa(chained, fill(-1.0, 5); xl = fill(-2.0, 5), xu = fill(0.8, 5))...)
emithc("bobyqa_chained5_npt7", PRIMA.bobyqa(chained, fill(-1.0, 5); xl = fill(-2.0, 5), xu = fill(2.0, 5), npt = 7)...)
emithc("bobyqa_chained5_npt21", PRIMA.bobyqa(chained, fill(-1.0, 5); xl = fill(-2.0, 5), xu = fill(2.0, 5),
                                            npt = 21)...)
emithc("bobyqa_chained10_npt30", PRIMA.bobyqa(chained, fill(-1.0, 10); xl = fill(-3.0, 10), xu = fill(0.9, 10),
                                             npt = 30)...)
emithc("bobyqa_one", PRIMA.bobyqa(one, [0.0]; xl = [-5.0], xu = [5.0])...)
emithc("bobyqa_one_far", PRIMA.bobyqa(one, [1000.0]; xl = [-2000.0], xu = [2000.0])...)
emithc("bobyqa_rosen_maxfun", PRIMA.bobyqa(rosen, [-1.2, 1.0]; xl = [-2.0, -2.0], xu = [2.0, 2.0], maxfun = 30)...)
emithc("bobyqa_rosen_ftarget", PRIMA.bobyqa(rosen, [-1.2, 1.0]; xl = [-2.0, -2.0], xu = [2.0, 2.0], ftarget = 1e-3)...)
emithc("bobyqa_rosen_scale", PRIMA.bobyqa(rosen, [-1.2, 1.0]; xl = [-2.0, -2.0], xu = [2.0, 2.0],
                                         scale = [2.0, 0.5])...)
emithc("bobyqa_rosen_maxfun3", PRIMA.bobyqa(rosen, [-1.2, 1.0]; xl = [-2.0, -2.0], xu = [2.0, 2.0], maxfun = 3)...)
emithc("bobyqa_rosen_rhoend0", PRIMA.bobyqa(rosen, [-1.2, 1.0]; xl = [-2.0, -2.0], xu = [2.0, 2.0], rhoend = 0.0)...)
emithc("bobyqa_nanrosen", PRIMA.bobyqa(nanrosen, [-1.2, 1.0]; xl = [-2.0, -2.0], xu = [2.0, 2.0])...)
emithc("bobyqa_infrosen", PRIMA.bobyqa(infrosen, [-1.2, 1.0]; xl = [-2.0, -2.0], xu = [2.0, 2.0])...)

# LINCOA: linear_ineq = (A, b) is A x <= b, linear_eq = (A, b) is A x = b
emithc("lincoa_rosen", PRIMA.lincoa(rosen, [-1.2, 1.0])...)
emithc("lincoa_rosen_ineq", PRIMA.lincoa(rosen, [-1.2, 1.0]; linear_ineq = ([1.0 1.0], [1.0]))...)
emithc("lincoa_rosen_eq_bounds", PRIMA.lincoa(rosen, [0.0, 0.0]; linear_eq = ([1.0 -1.0], [0.0]), xl = [-2.0, -2.0],
                                             xu = [0.5, 2.0])...)
emithc("lincoa_rosen_bounds", PRIMA.lincoa(rosen, [-1.2, 1.0]; xl = [-2.0, -2.0], xu = [0.5, 2.0])...)
emithc("lincoa_poly3_ineq", PRIMA.lincoa(poly3, [0.0, 0.0, 0.0]; linear_ineq = ([1.0 1.0 1.0; -1.0 0.0 0.0], [1.5, 0.0]))...)
# an infeasible x0 (the right hand sides are relaxed) and a zero row
emithc("lincoa_poly3_infeasible", PRIMA.lincoa(poly3, [1.0, 1.0, 1.0]; linear_ineq = ([1.0 1.0 1.0], [1.5]))...)
emithc("lincoa_rosen_zerorow", PRIMA.lincoa(rosen, [-1.2, 1.0]; linear_ineq = ([0.0 0.0; 1.0 1.0], [1.0, 1.0]))...)
emithc("lincoa_chained5", PRIMA.lincoa(chained, fill(-1.0, 5); linear_ineq = (ones(1, 5), [3.0]))...)
emithc("lincoa_chained5_npt7", PRIMA.lincoa(chained, fill(-1.0, 5); linear_ineq = (ones(1, 5), [3.0]), npt = 7)...)
emithc("lincoa_chained5_npt21", PRIMA.lincoa(chained, fill(-1.0, 5); linear_ineq = (ones(1, 5), [3.0]), npt = 21)...)
emithc("lincoa_chained10_npt30", PRIMA.lincoa(chained, fill(-1.0, 10);
                                             linear_ineq = ([ones(1, 10); -ones(1, 5) zeros(1, 5)], [6.0, 0.0]),
                                             npt = 30)...)
emithc("lincoa_rosen_rho", PRIMA.lincoa(rosen, [-1.2, 1.0]; linear_ineq = ([1.0 1.0], [1.0]), rhobeg = 0.5,
                                       rhoend = 1e-8)...)
emithc("lincoa_rosen_maxfun", PRIMA.lincoa(rosen, [-1.2, 1.0]; linear_ineq = ([1.0 1.0], [1.0]), maxfun = 30)...)
emithc("lincoa_rosen_ftarget", PRIMA.lincoa(rosen, [-1.2, 1.0]; linear_ineq = ([1.0 1.0], [1.0]), ftarget = 0.2)...)
emithc("lincoa_rosen_scale", PRIMA.lincoa(rosen, [-1.2, 1.0]; linear_ineq = ([1.0 1.0], [1.0]), scale = [2.0, 0.5])...)
emithc("lincoa_nanrosen", PRIMA.lincoa(nanrosen, [-1.2, 1.0]; linear_ineq = ([1.0 1.0], [1.0]))...)
emithc("lincoa_infrosen", PRIMA.lincoa(infrosen, [-1.2, 1.0]; linear_ineq = ([1.0 1.0], [1.0]))...)
emithc("lincoa_one", PRIMA.lincoa(one, [0.0]; linear_ineq = ([1.0;;], [1.0]))...)

# seeded sweeps of BOBYQA and LINCOA: the inputs are printed with the result,
# so that the tests need no random number generator
using Random
function sweepf(x, w, c, q)
    s = 0.0
    for i in eachindex(x)
        t = x[i] - c[i]
        s = s + w[i] * (t * t)
    end
    for i in 1:length(x)-1
        u = x[i] * x[i+1] - 0.5
        s = s + q[i] * (u * u)
    end
    s
end
rng = MersenneTwister(20261112)
for k in 1:24
    n = rand(rng, 1:5)
    c = randn(rng, n); w = 0.5 .+ 2 .* rand(rng, n); q = rand(rng, n - 1)
    xl = randn(rng, n) .- 1.0; xu = xl .+ 0.3 .+ 3 .* rand(rng, n)
    x0 = rand(rng) < 0.2 ? 3 .* randn(rng, n) : xl .+ (xu .- xl) .* rand(rng, n)
    npt = rand(rng, unique([n + 2, 2 * n + 1, div((n + 1) * (n + 2), 2)]))
    rhobeg = rand(rng) < 0.5 ? 1.0 : 0.2 * minimum(xu .- xl)
    x, info = PRIMA.bobyqa(z -> sweepf(z, w, c, q), x0; xl, xu, npt, rhobeg)
    println("ref\$sweepbq_$k <- list(w = $(hv(w)), c = $(hv(c)), q = $(hv(q)), x0 = $(hv(x0)), xl = $(hv(xl)), ",
            "xu = $(hv(xu)), npt = $(npt)L, rhobeg = $(hx(rhobeg)), x = $(hv(x)), fx = $(hx(info.fx)), ",
            "nf = $(info.nf)L, status = $(Int(info.status))L)")
end
for k in 1:24
    n = rand(rng, 2:5)
    c = randn(rng, n); w = 0.5 .+ 2 .* rand(rng, n); q = rand(rng, n - 1)
    x0 = randn(rng, n)
    mi = rand(rng, 0:3); me = rand(rng) < 0.3 ? 1 : 0
    Ai = randn(rng, mi, n); bi = Ai * x0 .+ (rand(rng, mi) .- 0.2)
    Ae = randn(rng, me, n); be = Ae * x0 .+ 0.1 .* randn(rng, me)
    bnd = rand(rng) < 0.5
    xl = bnd ? x0 .- 0.5 .- rand(rng, n) : fill(-Inf, n)
    xu = bnd ? x0 .+ 0.5 .+ rand(rng, n) : fill(Inf, n)
    npt = rand(rng, unique([n + 2, 2 * n + 1, div((n + 1) * (n + 2), 2)]))
    kw = (; npt, xl, xu)
    mi > 0 && (kw = (; kw..., linear_ineq = (Ai, bi)))
    me > 0 && (kw = (; kw..., linear_eq = (Ae, be)))
    x, info = PRIMA.lincoa(z -> sweepf(z, w, c, q), x0; kw...)
    println("ref\$sweeplc_$k <- list(w = $(hv(w)), c = $(hv(c)), q = $(hv(q)), x0 = $(hv(x0)), ",
            "Aineq = $(mi > 0 ? hm(Ai) : "NULL"), bineq = $(mi > 0 ? hv(bi) : "NULL"), ",
            "Aeq = $(me > 0 ? hm(Ae) : "NULL"), beq = $(me > 0 ? hv(be) : "NULL"), ",
            "xl = $(bnd ? hv(xl) : "NULL"), xu = $(bnd ? hv(xu) : "NULL"), npt = $(npt)L, x = $(hv(x)), ",
            "fx = $(hx(info.fx)), nf = $(info.nf)L, status = $(Int(info.status))L, cstrv = $(hx(info.cstrv)))")
end
