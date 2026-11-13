# Reference values for tests/testthat/test-sensitivity-misc.R,
# test-sensitivity-dae.R and test-sensitivity-sde.R (SciMLSensitivity
# 7.119.14).
#
# Beyond tools/julia/Project.toml this script needs SciMLSensitivity
# (7.119.14), Zygote (0.7.16), QuadGK (2.11.3) and ForwardDiff (1.4.7): run
# it in a copy of that environment with these four added, e.g.
#   julia --project=<env> tools/scimlsensitivity-probe.jl
using SciMLSensitivity, OrdinaryDiffEq, StochasticDiffEq, Zygote, QuadGK, ForwardDiff, Printf
using Optimization, OptimizationOptimJL, DiffEqNoiseProcess, SciMLBase
using OrdinaryDiffEqBDF: DFBDF, DImplicitEuler, FBDF

pr(name, x) = println(name, " = ", repr(x))

# ---------------------------------------------------------------------------
# supports_functor_params, SundialsAdjoint, SensitivityAlg
for a in (GaussAdjoint(), GaussKronrodAdjoint(), QuadratureAdjoint(), InterpolatingAdjoint(),
          BacksolveAdjoint(), ForwardSensitivity(), ForwardDiffSensitivity(), SteadyStateAdjoint())
    println(nameof(typeof(a)), " supports_functor_params ", supports_functor_params(a))
end
sa = SundialsAdjoint()
println("SundialsAdjoint steps ", sa.steps, " interp ", sa.interp, " quad_error_control ", sa.quad_error_control)
try
    SundialsAdjoint(interp = :linear)
catch e
    println("SundialsAdjoint error: ", e.msg)
end
SensitivityAlg()

# ---------------------------------------------------------------------------
# ODEForwardSensitivityFunction: the right-hand side of the forward
# sensitivity problem of Lotka Volterra
function lv!(du, u, p, t)
    du[1] = p[1] * u[1] - p[2] * u[1] * u[2]
    du[2] = -p[3] * u[2] + p[4] * u[1] * u[2] + t * p[1]
    return nothing
end
p = [1.5, 1.0, 3.0, 1.0]
u0 = [1.0, 1.0]
fprob = ODEForwardSensitivityProblem(lv!, u0, (0.0, 2.0), p)
z = [1.2, 0.7, 0.1, -0.2, 0.3, 0.4, -0.5, 0.6, 0.7, -0.8]
dz = similar(z)
fprob.f(dz, z, p, 0.3)
pr("forward_rhs", dz)
println("forward f isa SensitivityFunction ", fprob.f isa SciMLSensitivity.SensitivityFunction)
println("adjoint f isa SensitivityFunction ", ODEAdjointProblem(solve(ODEProblem((u, p, t) -> -p[1] .* u,
    [1.0], (0.0, 1.0), [2.0]), Tsit5()), InterpolatingAdjoint(autojacvec = ZygoteVJP()), Tsit5(), [1.0],
    (out, u, p, t, i) -> (out .= 1.0)).f.f isa SciMLSensitivity.SensitivityFunction)

# ---------------------------------------------------------------------------
# AdjointSensitivityIntegrand: G = sum(u(2)) + int_0^2 u1^2 / 2 dt
lvo(u, p, t) = [p[1] * u[1] - p[2] * u[1] * u[2], -p[3] * u[2] + p[4] * u[1] * u[2] + t * p[1]]
oprob = ODEProblem(lvo, u0, (0.0, 2.0), p)
sol = solve(oprob, Tsit5(); abstol = 1e-12, reltol = 1e-12)
dgdu_d(out, u, p, t, i) = (out .= 1.0)
dgdu_c(out, u, p, t) = (out[1] = u[1]; out[2] = 0.0)
dgdp_c(out, u, p, t) = (out .= 0.0)
adj = ODEAdjointProblem(sol, QuadratureAdjoint(autojacvec = ZygoteVJP()), Tsit5(), [2.0], dgdu_d, nothing, dgdu_c, dgdp_c)
adj_sol = solve(adj, Tsit5(); abstol = 1e-12, reltol = 1e-12)
integrand = AdjointSensitivityIntegrand(sol, adj_sol, QuadratureAdjoint(autojacvec = ZygoteVJP()))
for tq in (0.0, 0.5, 1.3, 2.0)
    pr("integrand($tq)", vec(integrand(tq)))
end
dp, _ = quadgk(integrand, 0.0, 2.0; atol = 1e-12, rtol = 1e-12)
pr("integrand_quadgk", vec(dp))
dgdp_q(out, u, p, t) = (out .= [u[1] * u[2], 0.0, 0.0, t])
integrand2 = AdjointSensitivityIntegrand(sol, adj_sol, QuadratureAdjoint(autojacvec = ZygoteVJP()), dgdp_q)
pr("integrand_dgdp(0.5)", vec(integrand2(0.5)))

# ---------------------------------------------------------------------------
# StochasticTransformedFunction: drift minus the Ito correction J' vec(g)
fs(u, p, t) = [p[1] * u[1] - u[2], p[2] * u[2] + u[1] * t]
gs(u, p, t) = [p[3] * u[1] * u[2], p[4] * sin(u[2])]
sprob = SDEProblem(fs, gs, [1.0, 0.5], (0.0, 1.0), [1.01, -0.5, 0.3, 0.7])
ssol = solve(sprob, EM(); dt = 0.1, seed = 1)
tf = SciMLSensitivity.StochasticTransformedFunction(ssol, sprob.f, sprob.g)
pr("transformed_diagonal", tf([1.3, -0.4], sprob.p, 0.25))
gm(u, p, t) = [p[3] * u[1] u[2]; p[4] * u[1] * u[2] cos(u[1])]
mprob = SDEProblem(fs, gm, [1.0, 0.5], (0.0, 1.0), [1.01, -0.5, 0.3, 0.7]; noise_rate_prototype = zeros(2, 2))
msol = solve(mprob, EM(); dt = 0.1, seed = 1)
tm = SciMLSensitivity.StochasticTransformedFunction(msol, mprob.f, mprob.g)
pr("transformed_general", tm([1.3, -0.4], mprob.p, 0.25))
tc = SciMLSensitivity.StochasticTransformedFunction(ssol, sprob.f, sprob.g, (u, p, t) -> [u[1], 2u[2]])
pr("transformed_analytic", tc([1.3, -0.4], sprob.p, 0.25))

# ---------------------------------------------------------------------------
# UnconstrainedOptimizationAdjoint and OptimizationAdjoint
# unconstrained: UnconstrainedOptimizationAdjoint through the solve
fu(x, p) = (x[1] - p[1])^2 + p[2] * (x[2] - x[1]^2)^2 + p[3] * x[1]
optf = OptimizationFunction(fu, AutoForwardDiff())
pu = [1.0, 3.0, 0.5]
probu = OptimizationProblem(optf, [0.0, 0.0], pu)
solu = solve(probu, BFGS(); g_tol = 1e-14)
pr("unc_u", solu.u)
gr = Zygote.gradient(p -> sum(solve(remake(probu, p = p), BFGS(); g_tol = 1e-14,
    sensealg = UnconstrainedOptimizationAdjoint()).u .* [1.0, 2.0]), pu)
pr("unc_dp", gr[1])
gr2 = Zygote.gradient(p -> sum(solve(remake(probu, p = p), BFGS(); g_tol = 1e-14,
    sensealg = OptimizationAdjoint()).u .* [1.0, 2.0]), pu)
pr("unc_dp_kkt", gr2[1])

# constrained: an equality, an active and an inactive inequality, an active bound
fc(x, p) = (x[1] - p[1])^2 + (x[2] - p[2])^2 + p[4] * (x[3] - p[3])^2 + p[6] * x[1] * x[2]
cons(res, x, p) = (res[1] = x[1] + 2 * x[2] + x[3] - p[5]; res[2] = x[1]^2 + x[2]^2 - 1.0;
                   res[3] = x[3] - x[1] - 10.0)
pc = [2.0, 1.5, -1.0, 1.0, 1.5, 0.3]
optc = OptimizationFunction(fc, AutoForwardDiff(); cons = cons)
probc = OptimizationProblem(optc, [0.1, 0.1, 0.5], pc; lcons = [0.0, -Inf, -Inf], ucons = [0.0, 0.0, 0.0],
    lb = [-Inf, -Inf, 0.0], ub = [Inf, Inf, Inf])
solc = solve(probc, IPNewton(); g_tol = 1e-14, x_tol = 1e-14, f_tol = 1e-14, maxiters = 2000)
pr("con_u", solc.u)
cv = zeros(3); cons(cv, solc.u, pc)
pr("con_cons", cv)
dpc = adjoint_sensitivities(solc, nothing; sensealg = OptimizationAdjoint(),
    dgdu = (out, u, p, t, i) -> (out .= [1.0, 2.0, -1.0]))
pr("con_dp", dpc)

# ---------------------------------------------------------------------------
# DAEAdjointProblem: index 1, the cost on the algebraic variable as well
function dae1!(res, du, u, p, t)
    res[1] = -p[1] * u[1] + p[2] * u[1] * u[2] - du[1]
    res[2] = p[3] * u[1] - u[2] * u[3] - du[2]
    res[3] = u[3] - u[1]^2 - p[4] * u[2]
    return nothing
end
p1 = [0.8, 0.6, 1.2, 0.4]
u01 = [1.0, 0.5, 1.2]
du01 = [-0.5, 0.6, 0.0]
prob1 = DAEProblem(dae1!, du01, u01, (0.0, 1.5), p1, differential_vars = [true, true, false])
sol1 = solve(prob1, DFBDF(), abstol = 1e-10, reltol = 1e-10)
pr("dae1_nsteps", length(sol1.t))
pr("dae1_uT", sol1.u[end])
ts1 = [0.5, 1.0, 1.5]
sq(out, u, p, t, i) = (out .= u)
function G1(p, u0d)
    u3 = u0d[1]^2 + p[4] * u0d[2]
    u0 = [u0d[1], u0d[2], u3]
    du0 = [-p[1] * u0[1] + p[2] * u0[1] * u0[2], p[3] * u0[1] - u0[2] * u0[3], 0.0]
    _prob = DAEProblem(dae1!, promote_type(eltype(p), eltype(u0d)).(du0), promote_type(eltype(p), eltype(u0d)).(u0), (0.0, 1.5), promote_type(eltype(p), eltype(u0d)).(p), differential_vars = [true, true, false])
    _sol = solve(_prob, DFBDF(), abstol = 1e-12, reltol = 1e-12, saveat = ts1)
    sum(sum(abs2, u) / 2 for u in _sol.u if true) - (_sol.t[1] == 0.0 ? sum(abs2, _sol.u[1]) / 2 : 0.0)
end
pr("dae1_fd_p", ForwardDiff.gradient(p -> G1(p, u01[1:2]), p1))
pr("dae1_fd_u", ForwardDiff.gradient(u -> G1(p1, u), u01[1:2]))
for (name, alg, sa) in (("interp_fbdf", FBDF(), InterpolatingAdjoint(autojacvec = ReverseDiffVJP())),
                        ("interp_dfbdf", DFBDF(), InterpolatingAdjoint(autojacvec = ReverseDiffVJP())),
                        ("quad_fbdf", FBDF(), QuadratureAdjoint(autojacvec = ReverseDiffVJP())),
                        ("gauss_fbdf", FBDF(), GaussAdjoint(autojacvec = ReverseDiffVJP())),
                        ("gk_fbdf", FBDF(), GaussKronrodAdjoint(autojacvec = ReverseDiffVJP())))
    du0, dp = adjoint_sensitivities(sol1, alg; t = ts1, dgdu_discrete = sq, sensealg = sa,
                                    abstol = 1e-10, reltol = 1e-10)
    pr("dae1_$(name)_du0", vec(du0))
    pr("dae1_$(name)_dp", vec(dp))
end
# continuous cost int u1^2/2 + u3 dt
dgc(out, u, p, t) = (out .= [u[1], 0.0, 1.0])
du0, dp = adjoint_sensitivities(sol1, FBDF(); dgdu_continuous = dgc,
                                sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP()),
                                abstol = 1e-10, reltol = 1e-10)
pr("dae1_cont_du0", vec(du0))
pr("dae1_cont_dp", vec(dp))
function G1c(p)
    u0 = [1.0, 0.5, 1.0 + p[4] * 0.5]
    du0 = [-p[1] + p[2] * 0.5, p[3] - 0.5 * u0[3], 0.0]
    _prob = DAEProblem(dae1!, eltype(p).(du0), eltype(p).(u0), (0.0, 1.5), p, differential_vars = [true, true, false])
    _sol = solve(_prob, DFBDF(), abstol = 1e-12, reltol = 1e-12)
    # trapezoid on a fine grid is too coarse; use the dense output
    tt = range(0.0, 1.5, length = 3001)
    v = [(_sol(x)[1]^2 / 2 + _sol(x)[3]) for x in tt]
    h = step(tt)
    h / 3 * (v[1] + v[end] + 4 * sum(v[2:2:end-1]) + 2 * sum(v[3:2:end-2]))
end
pr("dae1_cont_fd_p", ForwardDiff.gradient(G1c, p1))

# DAEAdjointProblem: Hessenberg index 2 (Core3/dae_adjoint.jl of the tests)
function hess2!(res, du, u, p, t)
    res[1] = p[1] * u[2] + u[3] - du[1]
    res[2] = -p[2] * u[1] - du[2]
    res[3] = u[1] - sin(p[3] * t)
    return nothing
end
p2 = [1.0, 2.0, 3.0]
prob2 = DAEProblem(hess2!, [p2[3], 0.0, 0.0], [0.0, 1.0, p2[3] - p2[1]], (0.0, 1.0), p2,
                   differential_vars = [true, true, false])
sol2 = solve(prob2, DImplicitEuler(), adaptive = false, dt = 1e-3, initializealg = SciMLBase.NoInit())
pr("hess2_uT", sol2.u[end])
ts2 = collect(0.2:0.2:1.0)
dg2(out, u, p, t, i) = (out .= [u[1], u[2], 0.0])
du0, dp = adjoint_sensitivities(sol2, FBDF(); t = ts2, dgdu_discrete = dg2,
                                sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP()),
                                abstol = 1e-8, reltol = 1e-8)
pr("hess2_du0", vec(du0))
pr("hess2_dp", vec(dp))
for sa in (QuadratureAdjoint(autojacvec = ReverseDiffVJP()), BacksolveAdjoint(autojacvec = ReverseDiffVJP()))
    try
        adjoint_sensitivities(sol2, FBDF(); t = ts2, dgdu_discrete = dg2, sensealg = sa)
    catch e
        println("hess2 error ", nameof(typeof(sa)), ": ", e.msg)
    end
end
try
    adjoint_sensitivities(sol1, FBDF(); t = ts1, dgdu_discrete = sq,
                          sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP(), checkpointing = true))
catch e
    println("checkpoint error: ", e.msg)
end

# ---------------------------------------------------------------------------
# SDEAdjointProblem and RODEAdjointProblem on a fixed noise path
# a fixed noise path on the grid 0:0.05:1
tg = collect(0.0:0.05:1.0)
w1(t) = 0.3 * sin(3.1 * t) + 0.25 * sin(7.3 * t) - 0.1 * sin(13.0 * t)
w2(t) = 0.4 * sin(5.2 * t) - 0.2 * sin(11.1 * t)
w3(t) = 0.35 * sin(2.3 * t) + 0.15 * sin(17.0 * t)
grid2() = NoiseGrid(tg, [[w1(t), w2(t)] for t in tg])
grid1() = NoiseGrid(tg, [w1(t) for t in tg])
grid3() = NoiseGrid(tg, [[w1(t), w2(t), w3(t)] for t in tg])

function f!(du, u, p, t)
    du[1] = p[1] * u[1] - u[1] * u[2]
    du[2] = -p[2] * u[2] + u[1] * u[2] + 0.1 * t
    return nothing
end
function gdiag!(du, u, p, t)
    du[1] = p[3] * u[1]
    du[2] = p[4] * u[2]^2
    return nothing
end
function gmix!(du, u, p, t)
    du[1] = p[3] * u[1] * u[2]
    du[2] = p[4] * u[1]
    return nothing
end
function ggen!(du, u, p, t)
    du[1, 1] = p[3] * u[1]
    du[1, 2] = 0.2 * u[2]
    du[1, 3] = p[4] * u[1] * u[2]
    du[2, 1] = 0.1
    du[2, 2] = p[3] * u[2]^2
    du[2, 3] = p[4] * cos(u[1])
    return nothing
end
p = [1.2, 0.9, 0.3, 0.25]
u0 = [1.0, 0.8]
ts = [0.5, 1.0]
dg(out, u, p, t, i) = (out .= u)

function run(name, prob, alg)
    sol = solve(prob, alg; dt = 0.05)
    pr("$(name)_uT", sol.u[end])
    for (sname, sa) in (("back", BacksolveAdjoint(autojacvec = ReverseDiffVJP())),
                        ("interp", InterpolatingAdjoint(autojacvec = ReverseDiffVJP())))
        du0, dp = adjoint_sensitivities(sol, EulerHeun(); t = ts, dgdu_discrete = dg, sensealg = sa, dt = 0.05)
        pr("$(name)_$(sname)_du0", vec(du0))
        pr("$(name)_$(sname)_dp", vec(dp))
    end
end
run("diag_strat", SDEProblem(f!, gdiag!, u0, (0.0, 1.0), p, noise = grid2()), EulerHeun())
run("diag_ito", SDEProblem(f!, gdiag!, u0, (0.0, 1.0), p, noise = grid2()), EM())
run("scalar", SDEProblem(f!, gdiag!, u0, (0.0, 1.0), p, noise = grid1()), EulerHeun())
run("general", SDEProblem(f!, ggen!, u0, (0.0, 1.0), p, noise = grid3(), noise_rate_prototype = zeros(2, 3)),
    EulerHeun())
# diagonal noise whose terms mix the states
solm = solve(SDEProblem(f!, gmix!, u0, (0.0, 1.0), p, noise = grid2()), EulerHeun(); dt = 0.05)
pr("mix_uT", solm.u[end])
for (sname, sa) in (("back", BacksolveAdjoint(autojacvec = ReverseDiffVJP(), noisemixing = true)),
                    ("interp", InterpolatingAdjoint(autojacvec = ReverseDiffVJP(), noisemixing = true)),
                    ("back_nomix", BacksolveAdjoint(autojacvec = ReverseDiffVJP())))
    du0, dp = adjoint_sensitivities(solm, EulerHeun(); t = ts, dgdu_discrete = dg, sensealg = sa, dt = 0.05)
    pr("mix_$(sname)_du0", vec(du0))
    pr("mix_$(sname)_dp", vec(dp))
end
# an analytic Ito correction
solc = solve(SDEProblem(f!, gdiag!, u0, (0.0, 1.0), p, noise = grid2()), EM(); dt = 0.05)
du0, dp = adjoint_sensitivities(solc, EulerHeun(); t = ts, dgdu_discrete = dg, dt = 0.05,
    sensealg = BacksolveAdjoint(autojacvec = ReverseDiffVJP()),
    corfunc_analytical = (du, u, p, t) -> (du[1] = p[3]^2 * u[1]; du[2] = 2 * p[4]^2 * u[2]^3))
pr("diag_ito_cor_du0", vec(du0))
pr("diag_ito_cor_dp", vec(dp))

# RODE. rode_interp samples the cost at t = 0.5 with the stale state of the
# last drift evaluation (u(0.55) for RandomEM); R samples u(0.5), so the test
# holds the R values and quotes the Julia ones in a comment.
function frode!(du, u, p, t, W)
    du[1] = p[1] * u[1] * sin(W[1]) - p[2] * u[2]
    du[2] = p[3] * u[2] * W[2] + p[4] * u[1]
    return nothing
end
rprob = RODEProblem(frode!, u0, (0.0, 1.0), p, noise = grid2())
rsol = solve(rprob, RandomEM(); dt = 0.05)
pr("rode_uT", rsol.u[end])
for (sname, sa) in (("back", BacksolveAdjoint(autojacvec = ReverseDiffVJP())),
                    ("interp", InterpolatingAdjoint(autojacvec = ReverseDiffVJP())))
    du0, dp = adjoint_sensitivities(rsol, RandomEM(); t = ts, dgdu_discrete = dg, sensealg = sa, dt = 0.05)
    pr("rode_$(sname)_du0", vec(du0))
    pr("rode_$(sname)_dp", vec(dp))
end
# the cotangent fixed at the forward states: the adjoint of a solve saved at
# ts with the loss sum(abs2, u) / 2 (solve_gradient with saveat)
du0, dp = adjoint_sensitivities(rsol, RandomEM(); t = ts, dt = 0.05,
    sensealg = BacksolveAdjoint(autojacvec = ReverseDiffVJP()), dgdu_discrete = (out, u, p, t, i) -> (out .= rsol(t)))
pr("rode_backfix_du0", vec(du0))
pr("rode_backfix_dp", vec(dp))

# a discrete cost at every saved time, t0 included; with the cotangent fixed
# at the saved states it is the adjoint of a solve, the reference of
# solve_gradient(function(sol) sum(as.matrix(sol)^2), ...)
sas = (("back", BacksolveAdjoint(autojacvec = ReverseDiffVJP())),
       ("interp", InterpolatingAdjoint(autojacvec = ReverseDiffVJP())))
solg = solve(SDEProblem(f!, gdiag!, u0, (0.0, 1.0), p, noise = grid2()), EulerHeun(); dt = 0.05)
# the saved time near 0.5 (R: 0.49999999999999994): a cost at t = 0.5 lies
# above it and is applied before the BacksolveAdjoint checkpoint reset
pr("sde_t11", solg.t[11])
for (sname, sa) in sas
    du0, dp = adjoint_sensitivities(solg, EulerHeun(); t = solg.t, dt = 0.05, sensealg = sa,
        dgdu_discrete = (out, u, p, t, i) -> (out .= 2 .* u))
    pr("grid_$(sname)_du0", vec(du0))
    pr("grid_$(sname)_dp", vec(dp))
    du0, dp = adjoint_sensitivities(solg, EulerHeun(); t = solg.t, dt = 0.05, sensealg = sa,
        dgdu_discrete = (out, u, p, t, i) -> (out .= 2 .* solg.u[i]))
    pr("gridfix_$(sname)_du0", vec(du0))
    pr("gridfix_$(sname)_dp", vec(dp))
end
# continuous costs: Julia also subtracts their gradient in the adjoint
# diffusion; R does not (the test holds the R values)
for (sname, sa) in sas
    du0, dp = adjoint_sensitivities(solg, EulerHeun(); dt = 0.05, sensealg = sa,
        dgdu_continuous = (out, u, p, t) -> (out .= u))
    pr("cont_$(sname)_du0", vec(du0))
    pr("cont_$(sname)_dp", vec(dp))
    du0, dp = adjoint_sensitivities(solg, EulerHeun(); dt = 0.05, sensealg = sa,
        dgdu_continuous = (out, u, p, t) -> (out .= u),
        dgdp_continuous = (out, u, p, t) -> (out .= 0; out[1] = u[1] * u[2]))
    pr("contp_$(sname)_du0", vec(du0))
    pr("contp_$(sname)_dp", vec(dp))
end
# InterpolatingAdjoint solved with EM: the cost reads the stale state of the
# drift evaluation at the start of the step (R reads u(0.5))
sole = solve(SDEProblem(f!, gdiag!, u0, (0.0, 1.0), p, noise = grid2()), EM(); dt = 0.05)
du0, dp = adjoint_sensitivities(sole, EM(); t = ts, dt = 0.05,
    sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP()), dgdu_discrete = dg)
pr("em_interp_du0", vec(du0))
pr("em_interp_dp", vec(dp))
