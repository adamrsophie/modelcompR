# Reference values for the PRIMA ports (R/optimization-prima*.R): runs
# PRIMA.jl (libprima) and prints R assignments with full precision. Usage:
#   julia --project=<env with PRIMA> tools/prima-probe.jl > refs.R
using PRIMA

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
