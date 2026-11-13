# Reference values for the PRIMA ports (R/optimization-prima*.R): runs
# PRIMA.jl (libprima) and prints R assignments with full precision. Usage:
#   julia --project=<env with PRIMA> tools/prima-probe.jl > refs.R
# The UOBYQA and NEWUOA values are printed as hexadecimal
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
