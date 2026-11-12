# Reference runs of OrdinaryDiffEqAdamsBashforthMoulton (v2.2.1, out of place
# ConstantCache methods) for tests/testthat/test-abm.R. Run from the
# repository root:
#   julia --project=tools/julia tools/abm-probe.jl > abm.tsv
# Lines "S" hold alg, case, saved points, nf, naccept, nreject and the final
# state; lines "T" trace every accepted step (t, nf, nreject, cache order).
# The cases "<name>fix" rerun the variable step methods with the corrections
# of R/solvers-abm.R:
#   VCAB3, VCAB4, VCAB5, VCABM3, VCABM4, VCABM5: a rejected starting step
#     (Bogacki-Shampine or Runge-Kutta) leaves the step counter, the step
#     size history and the divided differences as they were; the installed
#     methods keep them, so the retried step builds on the rejected one.
#   VCABM: at order 1 the corrector starts from uprev instead of the
#     predicted value (the installed method adds the predictor's Euler term
#     twice); the derivative at the end of the step is evaluated also with
#     adaptive = false (the installed method leaves it at zero, so every
#     later step starts from a zero derivative); a derivative discontinuity
#     (a callback that changes u) restarts the method at order 1 (the
#     installed method keeps the divided differences from before it).
using OrdinaryDiffEqAdamsBashforthMoulton, OrdinaryDiffEqCore, SciMLBase
const ABM = OrdinaryDiffEqAdamsBashforthMoulton

lv(u, p, t) = [p[1] * u[1] - u[1] * u[2], -p[2] * u[2] + u[1] * u[2] + 0.1 * t]
const prob = ODEProblem(lv, [1.0, 1.0], (0.0, 2.0), [1.5, 3.0])
halve!(integ) = (integ.u = integ.u * 0.5)
const cb = DiscreteCallback((u, t, integ) -> t == 1.0, halve!)

order_of(integ) = hasproperty(integ.cache, :order) ? integ.cache.order : 0

function kwargs_for(case, adaptive)
    if adaptive
        tol = case == "loose" ? (reltol = 1e-2, abstol = 1e-2) :
              case == "mid" ? (reltol = 1e-4, abstol = 1e-4) : (reltol = 1e-6, abstol = 1e-8)
        case == "tstop" && return (; tol..., tstops = [0.73])
        case == "cb" && return (; tol..., tstops = [1.0], callback = cb)
        return tol
    end
    case == "tstop" && return (dt = 0.05, adaptive = false, tstops = [0.73])
    case == "cb" && return (dt = 0.05, adaptive = false, tstops = [1.0], callback = cb)
    return (dt = 0.05, adaptive = false)
end

function run(name, alg, case, adaptive)
    kw = kwargs_for(case, adaptive)
    integ = init(prob, alg; kw...)
    for _ in integ
        println("T\t", name, "\t", case, "\t", repr(integ.t), "\t", integ.stats.nf, "\t",
                integ.stats.nreject, "\t", order_of(integ))
    end
    s = solve(prob, alg; kw...)
    st = s.stats
    println("S\t", name, "\t", case, "\t", length(s.t), "\t", st.nf, "\t", st.naccept, "\t", st.nreject,
            "\t", join(repr.(s.u[end]), ","))
end

for (nm, a) in [("AB3", AB3()), ("AB4", AB4()), ("AB5", AB5()), ("ABM32", ABM32()), ("ABM43", ABM43()),
                ("ABM54", ABM54())], case in ("plain", "tstop", "cb")
    run(nm, a, case, false)
end
for (nm, a) in [("VCAB3", VCAB3()), ("VCAB4", VCAB4()), ("VCAB5", VCAB5()), ("VCABM3", VCABM3()),
                ("VCABM4", VCABM4()), ("VCABM5", VCABM5()), ("VCABM", VCABM())],
    case in ("plain", "mid", "loose", "tstop", "cb")
    run(nm, a, case, true)
end
run("VCABM", VCABM(), "fixed", false)

# The corrected order 1 step of VCABM (see the header).
@eval ABM @muladd function perform_step!(integrator, cache::VCABMConstantCache, repeat_step = false)
    @inbounds begin
        (; t, dt, uprev, u, f, p) = integrator
        (; dts, g, ϕ_n, ϕ_np1, ϕstar_n, ϕstar_nm1, order, max_order) = cache
        k1 = integrator.fsalfirst
        step = integrator.iter
        if integrator.derivative_discontinuity
            cache.order = order = 1
        end
        k = order
        tmp = dts[13]
        for i in 12:-1:1
            dts[i + 1] = dts[i]
        end
        dts[1] = dt
        ϕ_and_ϕstar!(cache, k1, k)
        g_coefs!(cache, k + 1)
        u = muladd(g[1], ϕstar_n[1], uprev)
        for i in 2:(k - 1)
            u = muladd(g[i], ϕstar_n[i], u)
        end
        du_np1 = f(u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        ϕ_np1!(cache, du_np1, k + 1)
        u = muladd(g[k], ϕ_np1[k], k == 1 ? uprev : u)
        if !integrator.opts.adaptive
            integrator.fsallast = f(u, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
        if integrator.opts.adaptive
            utilde = (g[k + 1] - g[k]) * ϕ_np1[k + 1]
            atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
            if OrdinaryDiffEqCore.get_EEst(integrator) > one(OrdinaryDiffEqCore.get_EEst(integrator))
                for i in 1:12
                    dts[i] = dts[i + 1]
                end
                dts[13] = tmp
                return nothing
            end
            integrator.fsallast = f(u, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            if step <= 4 || order < 3
                cache.order = min(order + 1, 3)
            else
                utildem2 = (g[k - 1] - g[k - 2]) * ϕ_np1[k - 1]
                utildem1 = (g[k] - g[k - 1]) * ϕ_np1[k]
                expand_ϕ_and_ϕstar!(cache, k + 1)
                ϕ_np1!(cache, integrator.fsallast, k + 2)
                utildep1 = dt * γstar[(k + 1) + 1] * ϕ_np1[k + 2]
                atmpm2 = calculate_residuals(utildem2, uprev, u, integrator.opts.abstol,
                    integrator.opts.reltol, integrator.opts.internalnorm, t)
                atmpm1 = calculate_residuals(utildem1, uprev, u, integrator.opts.abstol,
                    integrator.opts.reltol, integrator.opts.internalnorm, t)
                atmpp1 = calculate_residuals(utildep1, uprev, u, integrator.opts.abstol,
                    integrator.opts.reltol, integrator.opts.internalnorm, t)
                errm2 = integrator.opts.internalnorm(atmpm2, t)
                errm1 = integrator.opts.internalnorm(atmpm1, t)
                errp1 = integrator.opts.internalnorm(atmpp1, t)
                if max(errm2, errm1) <= OrdinaryDiffEqCore.get_EEst(integrator)
                    cache.order = order - 1
                elseif errp1 < OrdinaryDiffEqCore.get_EEst(integrator)
                    cache.order = min(order + 1, max_order)
                    OrdinaryDiffEqCore.set_EEst!(integrator, one(OrdinaryDiffEqCore.get_EEst(integrator)))
                end
            end
        end
        cache.ϕstar_nm1, cache.ϕstar_n = ϕstar_n, ϕstar_nm1
        integrator.k[1] = integrator.fsalfirst
        integrator.k[2] = integrator.fsallast
        integrator.u = u
        return nothing
    end
end
# A rejected starting step of the fixed order methods leaves no trace.
for C in (:VCAB3ConstantCache, :VCAB4ConstantCache, :VCAB5ConstantCache, :VCABM3ConstantCache,
          :VCABM4ConstantCache, :VCABM5ConstantCache)
    @eval ABM function perform_step!(integrator, cache::$C)
        s = integrator.derivative_discontinuity ? 1 : cache.step
        dts = copy(cache.dts)
        nm1 = copy(cache.ϕstar_nm1)
        n = copy(cache.ϕstar_n)
        perform_step!(integrator, cache, false)
        if s < cache.order && integrator.opts.adaptive &&
           OrdinaryDiffEqCore.get_EEst(integrator) > 1
            cache.step = s
            cache.dts .= dts
            cache.ϕstar_nm1 = nm1
            cache.ϕstar_n = n
        end
        return nothing
    end
end
for (nm, a) in [("VCAB3", VCAB3()), ("VCAB4", VCAB4()), ("VCAB5", VCAB5()), ("VCABM3", VCABM3()),
                ("VCABM4", VCABM4()), ("VCABM5", VCABM5()), ("VCABM", VCABM())],
    case in ("plain", "mid", "loose", "tstop", "cb")
    Base.invokelatest(run, nm * "fix", a, case, true)
end
Base.invokelatest(run, "VCABMfix", VCABM(), "fixed", false)
let s = Base.invokelatest(solve, ODEProblem((u, p, t) -> [1.0], [0.0], (0.0, 1.0)), VCABM(); dt = 0.1)
    println("U1\t", join(repr.(s.t), ","), "\t", join(repr.(first.(s.u)), ","))
end
