# Reference values of tests/testthat/test-nonlinear-polyalg.R: the forcing
# terms of EisenstatWalkerForcing2 in a Newton Krylov solve and the entries of
# TraceAll (NonlinearSolve 4.32.0, NonlinearSolveFirstOrder 2.10.0).
#
#   julia --project=<env> tools/nonlinearsolve-forcing-probe.jl
#
# Each step prints the forcing term used by the step (eta), the residual norm
# before the step (rnorm_prev after the post step hook) and after it (rnorm).
# The R test replays the residual norms through the forcing cache and expects
# the same forcing terms. The field eta is written with an escape to keep the
# file ASCII.
using NonlinearSolve
const LS = NonlinearSolve.NonlinearSolveFirstOrder.LinearSolve
const ETA = Symbol("\u03b7")

function trace_forcing(name, prob, forcing)
    alg = NewtonRaphson(; linsolve = LS.KrylovJL_GMRES(), forcing)
    cache = init(prob, alg; abstol = 1.0e-10)
    println("INIT ", name, " ", repr(cache.forcing_cache.rnorm))
    k = 0
    while !cache.force_stop && k < 30
        step!(cache)
        fc = cache.forcing_cache
        println("STEP ", name, " ", k, " ", repr(getfield(fc, ETA)), " ", repr(fc.rnorm_prev), " ", repr(fc.rnorm))
        k += 1
    end
end

circ(u, p) = [u[1] * u[1] + u[2] * u[2] - 4, u[1] * u[2] - 1]
broyden(u, p) = [(3 - 2 * u[i]) * u[i] - (i > 1 ? u[i - 1] : 0.0) - 2 * (i < length(u) ? u[i + 1] : 0.0) + 1
                 for i in eachindex(u)]

trace_forcing("circ", NonlinearProblem(circ, [3.0, 0.1]), EisenstatWalkerForcing2())
trace_forcing("broyden", NonlinearProblem(broyden, fill(-1.0, 6)), EisenstatWalkerForcing2())
# positional fields: eta_0, eta_max, gamma, alpha, safeguard, safeguard_threshold
trace_forcing("circ_nosg", NonlinearProblem(circ, [3.0, 0.1]), EisenstatWalkerForcing2(0.3, 0.6, 0.8, 1.5, false, 0.1))
trace_forcing("circ_sg", NonlinearProblem(circ, [3.0, 0.1]), EisenstatWalkerForcing2(0.8, 0.95, 0.9, 1.2, true, 0.05))

# TraceAll: the stored state of each trace entry (the problem has no
# Jacobian, so Julia differentiates with ForwardDiff and R with differences)
g(u, p) = [u[1] * u[1] - 2, u[2] - 1]
gprob = NonlinearProblem(g, [1.0, 0.0]; lb = [0.0, 0.0], ub = [2.0, 2.0])
for (name, alg) in (("Dogbox", Dogbox()), ("TRR", TrustRegionReflective()), ("BTR", BoundedTrustRegion()))
    sol = solve(gprob, alg; store_trace = Val(true), trace_level = TraceAll())
    for e in sol.trace.history
        s = e.storage
        println("TA ", name, " ", e.iteration, " ", repr(e.fnorm), " ", repr(e.stepnorm), " ", repr(e.condJ), " u=",
                repr(s.u), " fu=", repr(s.fu), " du=", repr(getfield(s, Symbol("\u03b4u"))), " J=", repr(s.J))
    end
end
