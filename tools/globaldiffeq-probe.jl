# Reference values of the GlobalDiffEq.jl port (R/globaldiffeq.R and the
# companion and adjoint files below, tests/testthat/test-globaldiffeq*.R).
#   julia --project=tools/julia --startup-file=no tools/globaldiffeq-probe.jl
# Out of place vector problems with polynomial right-hand sides, so libm does
# not enter the comparison (the adjoint section at the end also uses exp and
# sin). The adjoint section needs SciMLSensitivity in the environment.
using GlobalDiffEq, OrdinaryDiffEqTsit5, OrdinaryDiffEqLowOrderRK

lv(u, p, t) = [p[1] * u[1] - p[2] * u[1] * u[2], -p[3] * u[2] + p[4] * u[1] * u[2]]
prob = ODEProblem(lv, [1.0, 1.0], (0.0, 5.0), [1.5, 1.0, 3.0, 1.0])
tq = [0.3, 1.7, 4.1]

show(x) = println(repr(x))

for alg in (GLEE23(), GLEE24(), GLEE35(), MM5GEE())
    name = string(nameof(typeof(alg)))
    sol = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-8)
    println("# ", name, " adaptive")
    println("n = ", length(sol.t), " nf = ", sol.stats.nf, " naccept = ", sol.stats.naccept,
            " nreject = ", sol.stats.nreject)
    show(sol.t[2]); show(sol.u[end]); show(global_error_estimate(sol)[end])
    show(global_error_estimate(sol, 3))
    show(reduce(vcat, [sol(t) for t in tq]))
    solf = solve(prob, alg; dt = 0.125, adaptive = false)
    println("# ", name, " fixed dt = 0.125")
    println("n = ", length(solf.t), " nf = ", solf.stats.nf)
    show(solf.u[end]); show(solf.global_error[end])
    sols = solve(prob, alg; reltol = 1.0e-5, abstol = 1.0e-7, saveat = [0.5, 2.5])
    println("# ", name, " saveat")
    show(sols.t); show(reduce(vcat, sols.u)); show(reduce(vcat, sols.global_error))
end

# scalar decay with a length one vector state
dec(u, p, t) = [-p[1] * u[1] * u[1]]
probd = ODEProblem(dec, [2.0], (0.0, 3.0), [0.7])
for alg in (GLEE24(), MM5GEE())
    sol = solve(probd, alg; reltol = 1.0e-8, abstol = 1.0e-10)
    println("# decay ", nameof(typeof(alg)))
    println("n = ", length(sol.t), " nf = ", sol.stats.nf, " nreject = ", sol.stats.nreject)
    show(sol.u[end]); show(sol.global_error[end])
end

# GlobalRichardson: the last inner solve
for (inner, dt) in ((Tsit5(), 0.5), (Euler(), 0.25), (RK4(), 0.25))
    sol = solve(prob, GlobalRichardson(inner); dt = dt, reltol = 1.0e-6, abstol = 1.0e-8)
    println("# GlobalRichardson ", nameof(typeof(inner)), " dt = ", dt)
    println("n = ", length(sol.t), " nf = ", sol.stats.nf)
    show(sol.t[2]); show(sol.u[end])
end
sol = solve(prob, GlobalRichardson(Tsit5()); dt = 0.5)
println("# GlobalRichardson Tsit5 default tolerances")
println("n = ", length(sol.t)); show(sol.u[end])

# ---------------------------------------------------------------------------
# Companion global error estimation (R/globaldiffeq-companion.R,
# tests/testthat/test-globaldiffeq-companion.R). One line per value,
# `name = repr(value)`; failures print `name error: message`.
using LinearAlgebra
using OrdinaryDiffEqVerner, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF
import SciMLBase

pr(name, x) = println(name, " = ", repr(x))
function tryp(f, name)
    try
        pr(name, f())
    catch e
        println(name, " error: ", replace(sprint(showerror, e), "\n" => " "))
    end
end

lvc(u, p, t) = [1.5 * u[1] - u[1] * u[2], -3.0 * u[2] + u[1] * u[2]]
probc = ODEProblem(lvc, [1.0, 1.0], (0.0, 2.0))
refc = solve(probc, Vern9(); abstol = 1.0e-14, reltol = 1.0e-14)
pr("gee_ref_uT", refc.u[end])

# estimates for every equation and mode
for eq in (DefectCorrection, ErrorTransport), mode in (InterpolatingMode, SimultaneousMode), tol in (1.0e-4, 1.0e-6)
    tryp("gee_$(eq)_$(mode)_$(tol)") do
        global_error_estimate(probc, GlobalErrorEstimation(Tsit5(); equation = eq, mode = mode);
            abstol = tol, reltol = tol)
    end
end
tryp("gee_tsit5_1e-4_uT") do
    solve(probc, Tsit5(); abstol = 1.0e-4, reltol = 1.0e-4).u[end]
end

# the companion with its own tolerances and solver
tryp("gee_companion_vern7") do
    global_error_estimate(probc, GlobalErrorEstimation(Tsit5(); companion_alg = Vern7(),
        companion_abstol = 1.0e-12, companion_reltol = 1.0e-10); abstol = 1.0e-5, reltol = 1.0e-5)
end

# other dense solvers (DefectCorrection, InterpolatingMode, companion_alg = alg)
for (name, alg) in (("DP5", DP5()), ("Vern7", Vern7()), ("BS3", BS3()), ("Rosenbrock23", Rosenbrock23()),
                    ("Rodas5P", Rodas5P()), ("ImplicitEuler", ImplicitEuler()), ("FBDF", FBDF()))
    tryp("gee_alg_$(name)") do
        global_error_estimate(probc, GlobalErrorEstimation(alg); abstol = 1.0e-5, reltol = 1.0e-5)
    end
    tryp("gee_alg_$(name)_tsit5companion") do
        global_error_estimate(probc, GlobalErrorEstimation(alg; companion_alg = Tsit5());
            abstol = 1.0e-5, reltol = 1.0e-5)
    end
    tryp("gee_alg_$(name)_uT") do
        solve(probc, alg; abstol = 1.0e-5, reltol = 1.0e-5).u[end]
    end
    # the dense first derivative: interior times, a step node from the left
    # and the right, and the end point from the right
    dsol = try
        solve(probc, alg; abstol = 1.0e-6, reltol = 1.0e-6)
    catch e
        println("dense_$(name) error: ", replace(sprint(showerror, e), "\n" => " "))
        nothing
    end
    dsol === nothing && continue
    pr("dense_$(name)_n", length(dsol.t))
    tryp("dense_$(name)_d") do
        reduce(vcat, [dsol(t, Val{1}) for t in (0.3, 1.1)])
    end
    tryp("dense_$(name)_node_left") do
        dsol(dsol.t[3], Val{1})
    end
    tryp("dense_$(name)_node_right") do
        dsol(dsol.t[3], Val{1}; continuity = :right)
    end
    tryp("dense_$(name)_T_right") do
        dsol(2.0, Val{1}; continuity = :right)
    end
    # the solve of gee_alg_<name>: its statistics and the dense first
    # derivative at the midpoint of every step (the companion integrates it)
    tryp("dense5_$(name)") do
        esol = solve(probc, alg; abstol = 1.0e-5, reltol = 1.0e-5)
        (length(esol.t), esol.stats.naccept, esol.stats.nreject, esol.stats.nf,
            reduce(vcat, [esol((esol.t[i] + esol.t[i + 1]) / 2, Val{1}) for i in 1:(length(esol.t) - 1)]))
    end
end
# ImplicitEuler at the tolerance of the R test, with an explicit companion
tryp("gee_alg_ImplicitEuler_tsit5companion_1e-3") do
    global_error_estimate(probc, GlobalErrorEstimation(ImplicitEuler(); companion_alg = Tsit5());
        abstol = 1.0e-3, reltol = 1.0e-3)
end
tryp("gee_alg_ImplicitEuler_1e-3") do
    global_error_estimate(probc, GlobalErrorEstimation(ImplicitEuler()); abstol = 1.0e-3, reltol = 1.0e-3)
end
tryp("gee_alg_ImplicitEuler_1e-3_uT") do
    solve(probc, ImplicitEuler(); abstol = 1.0e-3, reltol = 1.0e-3).u[end]
end
# the ImplicitEuler self companion of gee_alg_ImplicitEuler step by step:
# the companion solve's statistics, its first times and how many of its
# right hand side evaluations fall off the forward step nodes (R: 829
# accepted steps, no rejects, 36 of 4144 evaluations off the nodes,
# eps(T) = [3.4e-15, -6.6e-16])
tryp("gee_ie_trace") do
    fsol = solve(probc, ImplicitEuler(); abstol = 1.0e-5, reltol = 1.0e-5, dense = true)
    rhs = GlobalDiffEq._defect_correction_rhs(probc, fsol, true)
    evals = Float64[]
    rec = (e, p, t) -> (push!(evals, t); rhs(e, p, t))
    opts = GlobalErrorEstimation(ImplicitEuler()).options
    csol = solve(ODEProblem{false}(rec, zeros(2), (0.0, 2.0)), ImplicitEuler();
        abstol = opts.companion_abstol, reltol = opts.companion_reltol, tstops = fsol.t)
    (length(fsol.t), fsol.t[2:5], csol.stats.naccept, csol.stats.nreject, csol.stats.nf,
        csol.t[1:6], count(t -> !(t in fsol.t), evals), length(evals), csol.u[end])
end
# a length one state, y' = 2 y on [0, 1]
tryp("gee_lin") do
    global_error_estimate(ODEProblem((u, p, t) -> p .* u, [1.0], (0.0, 1.0), 2.0), GlobalErrorEstimation(Tsit5());
        abstol = 1.0e-5, reltol = 1.0e-5)
end

# global error control
for eq in (DefectCorrection, ErrorTransport), mode in (InterpolatingMode, SimultaneousMode)
    tryp("geectl_$(eq)_$(mode)") do
        csol = solve(probc, GlobalErrorEstimation(Tsit5(); equation = eq, mode = mode, gtol = 1.0e-7);
            abstol = 1.0e-3, reltol = 1.0e-3)
        (length(csol.t), csol.stats.nf, csol.retcode, csol.u[end])
    end
end

# validation
tryp("gee_nogtol") do
    solve(probc, GlobalErrorEstimation(Tsit5()))
end
tryp("gee_callback") do
    cb = DiscreteCallback((u, t, integrator) -> false, integrator -> nothing)
    global_error_estimate(ODEProblem(lvc, [1.0, 1.0], (0.0, 2.0); callback = cb), GlobalErrorEstimation(Tsit5()))
end
tryp("gee_backwards") do
    global_error_estimate(ODEProblem(lvc, [1.0, 1.0], (2.0, 0.0)), GlobalErrorEstimation(Tsit5()))
end
tryp("gee_mass") do
    f = ODEFunction(lvc; mass_matrix = [1.0 0.0; 0.0 2.0])
    global_error_estimate(ODEProblem(f, [1.0, 1.0], (0.0, 2.0)), GlobalErrorEstimation(Rodas5P()))
end
tryp("gee_tolerance") do
    global_error_estimate(probc, GlobalErrorEstimation(Tsit5()); abstol = -1.0)
end
tryp("gee_companion_tolerance") do
    GlobalErrorEstimation(Tsit5(); companion_abstol = 0.0, companion_reltol = 0.0)
end
tryp("gee_gtol") do
    GlobalErrorEstimation(Tsit5(); gtol = -1.0)
end
tryp("gee_maxiters_fail") do
    solve(probc, GlobalErrorEstimation(Tsit5(); gtol = 1.0e-9, maxiters = 1); abstol = 1.0e-3, reltol = 1.0e-3)
end
tryp(() -> SciMLBase.has_global_error(GlobalErrorEstimation(Tsit5())), "has_global_error_gee")
tryp(() -> SciMLBase.has_global_error(GLEE24()), "has_global_error_glee")

# ---------------------------------------------------------------------------
# Adjoint global error estimation (R/globaldiffeq-adjoint.R,
# tests/testthat/test-globaldiffeq-adjoint.R). Needs SciMLSensitivity (the
# GlobalDiffEqSciMLSensitivityExt extension). The estimates below have one
# state or as many samples as states, so they do not depend on the draws.
using SciMLSensitivity, Random

linear!(du, u, p, t) = (du[1] = p * u[1]; nothing)
linear(u, p, t) = [p * u[1]]
lprob = ODEProblem(linear!, [1.0], (0.0, 2.0), 2.0)
tryp("ga_lin_est") do
    adjoint_error_estimate(lprob, GlobalAdjoint(Tsit5(); samples = 1, rng = Xoshiro(123));
        abstol = 1.0e-4, reltol = 1.0e-4)
end
tryp("ga_lin_est_oop") do
    adjoint_error_estimate(ODEProblem(linear, [1.0], (0.0, 2.0), 2.0),
        GlobalAdjoint(Tsit5(); samples = 1, rng = Xoshiro(123)); abstol = 1.0e-4, reltol = 1.0e-4)
end
tryp("ga_lin_est_quadrature") do
    adjoint_error_estimate(lprob, GlobalAdjoint(Tsit5(); samples = 1, sensealg = QuadratureAdjoint());
        abstol = 1.0e-4, reltol = 1.0e-4)
end
tryp("ga_lin_actual") do
    abs(solve(lprob, Tsit5(); abstol = 1.0e-4, reltol = 1.0e-4).u[end][1] - exp(4.0))
end
tryp("ga_lin_ctl") do
    sol = solve(lprob, GlobalAdjoint(Tsit5(); gtol = 1.0e-6, samples = 1, rng = Xoshiro(321));
        abstol = 1.0e-3, reltol = 1.0e-3)
    (length(sol.t), sol.stats.nf, sol.u[end][1], sol.global_error[end], count(!iszero, sol.global_error))
end
for (tag, tspan, gtol, tol) in (("", (0.0, 2.0), 1.0e-5, 1.0e-3), ("_cheap", (0.0, 1.0), 1.0e-3, 1.0e-3))
    p2 = ODEProblem(linear!, [1.0], tspan, 2.0)
    tryp("ga_traj$(tag)") do
        sol = solve(p2, GlobalAdjoint(Tsit5(); gtol, scope = TrajectoryError, samples = 1, rng = Xoshiro(11));
            abstol = tol, reltol = tol)
        (sol.t, sol.global_error, [u[1] for u in sol.u])
    end
    tryp("ga_grid$(tag)") do
        sol = solve(p2, GlobalAdjoint(Tsit5(); gtol, control = StepGridRefinement, samples = 1, rng = Xoshiro(22));
            abstol = tol, reltol = tol)
        (length(sol.t), sol.retcode, sol.t[2], sol.u[end][1], sol.global_error[end])
    end
end

function lvi!(du, u, p, t)
    du[1] = 1.5 * u[1] - u[1] * u[2]
    du[2] = -3.0 * u[2] + u[1] * u[2]
    return nothing
end
lvprob = ODEProblem(lvi!, [1.0, 1.0], (0.0, 2.0))
tryp("ga_lv_est") do
    adjoint_error_estimate(lvprob, GlobalAdjoint(Tsit5(); samples = 2, rng = Xoshiro(5));
        abstol = 1.0e-5, reltol = 1.0e-5)
end
tryp("ga_lv_actual") do
    norm(solve(lvprob, Tsit5(); abstol = 1.0e-5, reltol = 1.0e-5).u[end] - refc.u[end])
end
tryp("ga_lv_ctl") do
    sol = solve(lvprob, GlobalAdjoint(Tsit5(); gtol = 1.0e-6, samples = 2, rng = Xoshiro(6));
        abstol = 1.0e-3, reltol = 1.0e-3)
    (length(sol.t), sol.u[end], sol.global_error[end])
end

tryp("ga_callback") do
    cb = DiscreteCallback((u, t, integrator) -> false, integrator -> nothing)
    adjoint_error_estimate(ODEProblem(linear!, [1.0], (0.0, 2.0), 2.0; callback = cb),
        GlobalAdjoint(Tsit5(); samples = 1, rng = Xoshiro(789)); abstol = 1.0e-4, reltol = 1.0e-4)
end
forced(u, p, t) = [-u[1] + p(t)]
tryp("ga_forced") do
    adjoint_error_estimate(ODEProblem(forced, [0.0], (0.0, 1.0), t -> 0.1 * sin(t)),
        GlobalAdjoint(Tsit5(); samples = 1, rng = Xoshiro(456)); abstol = 1.0e-4, reltol = 1.0e-4)
end
tryp("ga_nogtol") do
    solve(lprob, GlobalAdjoint(Tsit5()))
end
tryp("ga_samples") do
    GlobalAdjoint(Tsit5(); samples = 0)
end
tryp("ga_gtol") do
    GlobalAdjoint(Tsit5(); gtol = -1.0)
end
tryp(() -> SciMLBase.has_global_error(GlobalAdjoint(Tsit5())), "has_global_error_ga")
# the directions of a three state problem from Xoshiro(123) (R draws with rng)
tryp(() -> randn(Xoshiro(123), 3), "ga_randn_123")

# ---------------------------------------------------------------------------
# The ImplicitEuler self companion of gee_alg_ImplicitEuler: the initial step
# size pieces (initdt.jl: f0, f1 at t0 + 1e-6), the initial dt of the
# integrator, the first companion steps and states and the first evaluations.
tryp("gee_ie_trace2") do
    fsol = solve(probc, ImplicitEuler(); abstol = 1.0e-5, reltol = 1.0e-5, dense = true)
    rhs = GlobalDiffEq._defect_correction_rhs(probc, fsol, true)
    evals = Float64[]
    rec = (e, p, t) -> (push!(evals, t); rhs(e, p, t))
    opts = GlobalErrorEstimation(ImplicitEuler()).options
    cprob = ODEProblem{false}(rec, zeros(2), (0.0, 2.0))
    integ = SciMLBase.init(cprob, ImplicitEuler(); abstol = opts.companion_abstol,
        reltol = opts.companion_reltol, tstops = fsol.t)
    dt0 = integ.dt
    empty!(evals)
    csol = solve(cprob, ImplicitEuler(); abstol = opts.companion_abstol, reltol = opts.companion_reltol,
        tstops = fsol.t)
    (rhs(zeros(2), nothing, 0.0), rhs(zeros(2), nothing, 1.0e-6), dt0, csol.t[1:12],
        reduce(vcat, csol.u[2:12]), evals[1:40], opts.companion_abstol, opts.companion_reltol,
        csol.stats.naccept, csol.stats.nreject, csol.stats.nnonlinconvfail)
end

# ---------------------------------------------------------------------------
# Dense output of the methods without an interpolant of their own
# (R/solvers-erk.R .ERK_DENSE, tests/testthat/test-dense-output.R): fixed
# steps of 0.1 on a nonautonomous problem, the values and first derivatives
# at points in the first and later steps, a saveat solve per kind of dense
# output and the companion estimate of GlobalDiffEq.
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK, OrdinaryDiffEqRKN, OrdinaryDiffEqNewmark
lvt(u, p, t) = [1.5 * u[1] - u[1] * u[2], -3.0 * u[2] + u[1] * u[2] + 0.1 * t]
probt = ODEProblem(lvt, [1.0, 1.0], (0.0, 1.0))
dtq = (0.03, 0.47, 0.71)
dense_algs = (
    "Alshina2" => Alshina2(), "Alshina3" => Alshina3(), "Alshina6" => Alshina6(),
    "CarpenterKennedy2N54" => CarpenterKennedy2N54(), "DGLDDRK73_C" => DGLDDRK73_C(),
    "DGLDDRK84_C" => DGLDDRK84_C(), "DGLDDRK84_F" => DGLDDRK84_F(), "NDBLSRK124" => NDBLSRK124(),
    "NDBLSRK134" => NDBLSRK134(), "NDBLSRK144" => NDBLSRK144(), "ORK256" => ORK256(),
    "pRRK22" => pRRK22(kappa = 0.5), "pRRK33" => pRRK33(kappa = 0.5), "pRRK54" => pRRK54(kappa = 0.5),
    "PSRK3p5q4" => PSRK3p5q4(), "PSRK3p6q5" => PSRK3p6q5(), "PSRK4p7q6" => PSRK4p7q6(),
    "RDPK3Sp35" => RDPK3Sp35(), "RDPK3Sp49" => RDPK3Sp49(), "RDPK3Sp510" => RDPK3Sp510(),
    "RK46NL" => RK46NL(), "RKO65" => RKO65(), "RKV76IIa" => RKV76IIa(), "SHLDDRK64" => SHLDDRK64(),
    "SSPRK104" => SSPRK104(), "SSPRK22" => SSPRK22(), "SSPRK33" => SSPRK33(), "SSPRK43" => SSPRK43(),
    "SSPRK432" => SSPRK432(), "SSPRK53" => SSPRK53(), "SSPRK53_2N1" => SSPRK53_2N1(),
    "SSPRK53_2N2" => SSPRK53_2N2(), "SSPRK53_H" => SSPRK53_H(), "SSPRK54" => SSPRK54(),
    "SSPRK63" => SSPRK63(), "SSPRK73" => SSPRK73(), "SSPRK83" => SSPRK83(), "SSPRK932" => SSPRK932(),
)
for (name, alg) in dense_algs
    tryp("dense41_$(name)") do
        s = solve(probt, alg; dt = 0.1, adaptive = false)
        (s.stats.nf, s.u[end], reduce(vcat, [s(t) for t in dtq]), reduce(vcat, [s(t, Val{1}) for t in dtq]))
    end
end
for (name, alg) in ("SSPRK33" => SSPRK33(), "SSPRK53" => SSPRK53(), "Alshina3" => Alshina3(),
        "RKO65" => RKO65(), "pRRK33" => pRRK33(kappa = 0.5), "RK46NL" => RK46NL())
    tryp("dense41_saveat_$(name)") do
        s = solve(probt, alg; dt = 0.1, adaptive = false, saveat = [0.25, 0.65])
        (s.stats.nf, s.t, reduce(vcat, s.u))
    end
end
for (name, alg) in ("SSPRK33" => SSPRK33(), "SSPRK53" => SSPRK53(), "Alshina2" => Alshina2(),
        "RKO65" => RKO65(), "pRRK33" => pRRK33(kappa = 0.5))
    tryp("dense41_gee_$(name)") do
        global_error_estimate(probc, GlobalErrorEstimation(alg; companion_alg = Tsit5()); dt = 0.01)
    end
end
for (name, alg) in ("SSPRK43" => SSPRK43(), "RDPK3Sp35" => RDPK3Sp35(), "SSPRK932" => SSPRK932())
    tryp("dense41_gee_$(name)") do
        global_error_estimate(probc, GlobalErrorEstimation(alg); abstol = 1.0e-5, reltol = 1.0e-5)
    end
end
# DImplicitEuler on the index 1 DAE of tools/irkc-probe.jl: k = [du(t1), du(t0)]
Fd(du, u, p, t) = [du[1] - u[2], u[2] * u[2] * u[2] + u[2] + u[1] * u[1] * u[1] + u[1]]
Jd(du, u, p, gamma, t) = [gamma -1.0; 3u[1]^2+1 3u[2]^2+1]
dprobd = DAEProblem(DAEFunction(Fd; jac = Jd), [-1.0, -1.0], [1.0, -1.0], (0.0, 2.0); differential_vars = [true, false])
tryp("dense41_DImplicitEuler") do
    s = solve(dprobd, DImplicitEuler(); reltol = 1.0e-4, abstol = 1.0e-6, initializealg = SciMLBase.NoInit())
    tt = (s.t[2] / 2, 0.3, 1.1)
    (length(s.t), s.stats.nreject, reduce(vcat, [s(t) for t in tt]), reduce(vcat, [s(t, Val{1}) for t in tt]))
end
# DPRKN6 (Val{0} of its own, Val{1} the generic Hermite on its stages) and
# NewmarkBeta on u'' = -u - u^3 (state (v, u))
p1 = DynamicalODEProblem((v, u, p, t) -> -u .- u .* u .* u, (v, u, p, t) -> v, [0.5, -0.3], [1.0, 0.2], (0.0, 1.0))
tryp("dense41_DPRKN6") do
    s = solve(p1, DPRKN6())
    (length(s.t), reduce(vcat, [collect(s(t)) for t in dtq]), reduce(vcat, [collect(s(t, Val{1})) for t in dtq]))
end
tryp("dense41_NewmarkBeta") do
    s = solve(p1, NewmarkBeta(); dt = 0.1, adaptive = false)
    (length(s.t), reduce(vcat, [collect(s(t)) for t in dtq]), reduce(vcat, [collect(s(t, Val{1})) for t in dtq]))
end

# The stiff forward solvers with a Tsit5 companion (the forward interpolant
# alone) and the stats of their self companion solves
for (name, alg) in ("Rosenbrock23" => Rosenbrock23(), "Rodas5P" => Rodas5P(), "FBDF" => FBDF())
    tryp("dense41_geetsit_$(name)") do
        global_error_estimate(probc, GlobalErrorEstimation(alg; companion_alg = Tsit5());
            abstol = 1.0e-5, reltol = 1.0e-5)
    end
    tryp("dense41_geeself_$(name)") do
        fsol = solve(probc, alg; abstol = 1.0e-5, reltol = 1.0e-5, dense = true)
        rhs = GlobalDiffEq._defect_correction_rhs(probc, fsol, true)
        opts = GlobalErrorEstimation(alg).options
        csol = solve(ODEProblem{false}(rhs, zeros(2), (0.0, 2.0)), alg; abstol = opts.companion_abstol,
            reltol = opts.companion_reltol, tstops = fsol.t)
        (csol.stats.naccept, csol.stats.nreject, csol.stats.nf, csol.stats.nnonlinconvfail, csol.u[end])
    end
end
# Val{1} of methods accepted through the shared registrations: DABDF2 on the
# DAE above, PDIRK44 and CFRLDDRK64 on probt
using OrdinaryDiffEqPDIRK
tryp("dense41_DABDF2") do
    s = solve(dprobd, DABDF2(); reltol = 1.0e-4, abstol = 1.0e-6, initializealg = SciMLBase.NoInit())
    tt = (s.t[2] / 2, 0.3, 1.1)
    (length(s.t), reduce(vcat, [s(t) for t in tt]), reduce(vcat, [s(t, Val{1}) for t in tt]))
end
for (name, alg) in ("PDIRK44" => PDIRK44(), "CFRLDDRK64" => CFRLDDRK64())
    tryp("dense41_$(name)") do
        s = solve(probt, alg; dt = 0.1, adaptive = false)
        (reduce(vcat, [s(t) for t in dtq]), reduce(vcat, [s(t, Val{1}) for t in dtq]))
    end
end
# How @muladd nests a sum of three products
tryp("dense41_muladd") do
    MM = Base.loaded_modules[Base.PkgId(Base.UUID("46d2c3a1-f734-5fdb-9937-b9b9aeba4221"), "MuladdMacro")]
    string(MM.to_muladd(:(a * b + c * d + e * f)))
end
