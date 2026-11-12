# Reference runs of OrdinaryDiffEqNordsieck (v2.3.2, out of place
# ConstantCache methods) for tests/testthat/test-nordsieck.R. Run from the
# repository root:
#   julia --project=tools/julia tools/nordsieck-probe.jl > nordsieck.tsv
# The problem is a Lotka-Volterra system with a time dependent term on
# (0, 2); the case "long" integrates to t = 10, where the Adams order also
# decreases, and the case "stiff" takes p[2] = 20 with tolerance 0.1, where
# functional iterations fail. Lines "S" hold alg, case, saved points, nf, naccept, nreject,
# the number of failed functional iterations (force_stepfail, not counted in
# nreject), the final state and the dense output at t = 0.005, 0.3 and 1.5;
# lines "T" trace every accepted step (t, nf, nreject, order).
# The runs named "...fix" use the corrections of R/solvers-nordsieck.R:
# AN5
#   a step rejected by its error estimate rewinds the Nordsieck vector and
#     restores the step size history without correcting it (the installed
#     method corrects the predicted vector and keeps it, so the retried step
#     starts from a vector at the rejected time);
#   a failed functional iteration also restores the step size history;
#   the starting step clears the corrector coefficients as a fresh cache
#     does (after a restart the installed method applies the coefficients of
#     the step before the callback);
#   the dense output of the starting step uses the derivative at its end
#     (the installed method passes the second Tsit5 stage).
# JVODE
#   an Adams order decrease measures the past steps in units of the step the
#     Nordsieck vector is scaled to (dts[2]; the installed method divides by
#     the new step dts[1]);
#   a restart after a callback also resets L and n_wait as the fresh cache
#     does.
using OrdinaryDiffEqNordsieck, OrdinaryDiffEqCore, SciMLBase
const NS = OrdinaryDiffEqNordsieck

lv(u, p, t) = [p[1] * u[1] - u[1] * u[2], -p[2] * u[2] + u[1] * u[2] + 0.1 * t]
const prob = ODEProblem(lv, [1.0, 1.0], (0.0, 2.0), [1.5, 3.0])
const prob10 = ODEProblem(lv, [1.0, 1.0], (0.0, 10.0), [1.5, 3.0])
const probstiff = ODEProblem(lv, [1.0, 1.0], (0.0, 2.0), [1.5, 20.0])
halve!(integ) = (integ.u = integ.u * 0.5)
const cb = DiscreteCallback((u, t, integ) -> t == 1.0, halve!)
const FAILS = Ref(0)

# count the failed functional iterations; redefined after every patch of
# perform_step! (whose default argument also defines the two argument method)
function count_fails()
    for C in (:AN5ConstantCache, :JVODEConstantCache)
        @eval NS function perform_step!(integrator, cache::$C)
            perform_step!(integrator, cache, false)
            integrator.force_stepfail && (Main.FAILS[] += 1)
            return nothing
        end
    end
end
count_fails()

function kwargs_for(case)
    tol = case == "stiff" ? (reltol = 1e-1, abstol = 1e-1) :
          case == "loose" ? (reltol = 1e-2, abstol = 1e-2) :
          case in ("mid", "long") ? (reltol = 1e-4, abstol = 1e-4) : (reltol = 1e-6, abstol = 1e-8)
    case == "tstop" && return (; tol..., tstops = [0.73])
    case == "cb" && return (; tol..., tstops = [1.0], callback = cb)
    case == "fixed" && return (dt = 0.05, adaptive = false)
    case == "fixedcb" && return (dt = 0.05, adaptive = false, tstops = [0.73, 1.0], callback = cb)
    return tol
end

function run(name, alg, case)
    kw = kwargs_for(case)
    pr = case == "long" ? prob10 : case == "stiff" ? probstiff : prob
    integ = init(pr, alg; kw...)
    for _ in integ
        println("T\t", name, "\t", case, "\t", repr(integ.t), "\t", integ.stats.nf, "\t",
                integ.stats.nreject, "\t", integ.cache.order)
    end
    FAILS[] = 0
    s = solve(pr, alg; kw...)
    st = s.stats
    dense = vcat(s(0.005), s(0.3), s(1.5))
    println("S\t", name, "\t", case, "\t", length(s.t), "\t", st.nf, "\t", st.naccept, "\t", st.nreject,
            "\t", FAILS[], "\t", join(repr.(s.u[end]), ","), "\t", join(repr.(dense), ","))
end

const CASES = ("plain", "mid", "loose", "tstop", "cb", "fixed", "fixedcb", "long", "stiff")
for (nm, a) in [("AN5", AN5()), ("JVODE_Adams", JVODE_Adams()), ("JVODE_BDF", JVODE_BDF())], case in CASES
    run(nm, a, case)
end

# lines "F": a functional iteration that fails without adaptivity ends the
# solve (retcode, saved points, nf, naccept, nreject, final state)
for (nm, a) in [("AN5", AN5()), ("JVODE_Adams", JVODE_Adams()), ("JVODE_BDF", JVODE_BDF())]
    s = solve(probstiff, a; dt = 0.2, adaptive = false)
    println("F\t", nm, "\t", s.retcode, "\t", length(s.t), "\t", s.stats.nf, "\t", s.stats.naccept, "\t",
            s.stats.nreject, "\t", join(repr.(s.u[end]), ","))
end

@eval NS @muladd function perform_step!(integrator, cache::AN5ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, differential_vars) = integrator
    (; z, l, m, c_LTE, dts, tsit5tab) = cache
    if integrator.derivative_discontinuity
        cache.order = 1
    end
    if cache.order == 1
        perform_step!(integrator, tsit5tab, repeat_step)
        cache.order = 4
        z[1] = integrator.uprev
        z[2] = integrator.k[1] * dt
        z[3] = ode_interpolant(
            zero(t), dt, nothing, nothing, integrator.k, tsit5tab, nothing,
            Val{2}, differential_vars
        ) * dt^2 / 2
        z[4] = ode_interpolant(
            zero(t), dt, nothing, nothing, integrator.k, tsit5tab, nothing,
            Val{3}, differential_vars
        ) * dt^3 / 6
        z[5] = ode_interpolant(
            zero(t), dt, nothing, nothing, integrator.k, tsit5tab, nothing,
            Val{4}, differential_vars
        ) * dt^4 / 24
        z[6] = zero(cache.z[6])
        fill!(dts, dt)
        perform_predict!(cache)
        cache.Δ = integrator.u - integrator.uprev
        fill!(l, zero(eltype(l)))
        update_nordsieck_vector!(cache)
        integrator.k[2] = integrator.fsallast
        if integrator.opts.adaptive && get_EEst(integrator) >= one(get_EEst(integrator))
            cache.order = 1
        end
    else
        tmp = dts[6]
        for i in 5:-1:1
            dts[i + 1] = dts[i]
        end
        dts[1] = dt
        dt != dts[2] && nordsieck_rescale!(cache)
        integrator.k[1] = z[2] / dt
        perform_predict!(cache)
        calc_coeff!(cache)
        isucceed = nlsolve_functional!(integrator, cache)
        if !isucceed
            integrator.force_stepfail = true
            nordsieck_restore_failed_step!(cache, tmp)
            return nothing
        end
        if integrator.opts.adaptive
            atmp = calculate_residuals(
                cache.Δ, uprev, integrator.u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t
            )
            set_EEst!(integrator, integrator.opts.internalnorm(atmp, t) * cache.c_LTE)
            if get_EEst(integrator) > one(get_EEst(integrator))
                nordsieck_restore_failed_step!(cache, tmp)
                return nothing
            end
        end
        cache.order = 5
        update_nordsieck_vector!(cache)
        integrator.k[2] = cache.z[2] / dt
    end
    return nothing
end
count_fails()
for case in CASES
    Base.invokelatest(run, "AN5fix", AN5(), case)
end

# JVODE corrections: an Adams order decrease measures the past steps in units
# of the step the Nordsieck vector is scaled to (dts[2]; the installed method
# divides by the new step dts[1]), and a restart after a callback also resets
# L and n_wait as the fresh cache does.
@eval NS function nordsieck_adjust_order!(cache::T, dorder) where {T}
    (; order, dts) = cache
    if dorder == 1
        cache.z[order + 2] = zero(cache.z[order + 2])
    else
        cache.l .= 0
        cache.l[2] = 1
        dt = dts[2]
        hsum = zero(eltype(cache.dts))
        for j in 1:(order - 2)
            hsum += cache.dts[j + 1]
            ξ = hsum / dt
            for i in (j + 1):-1:1
                cache.l[i + 1] = cache.l[i + 1] * ξ + cache.l[i]
            end
        end
        for j in 2:(order - 1)
            cache.l[j + 1] = order * cache.l[j] / j
        end
        for j in 3:order
            cache.z[j] = muladd.(-cache.l[j], cache.z[order + 1], cache.z[j])
        end
    end
    return nothing
end
@eval NS @muladd function perform_step!(integrator, cache::JVODEConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, differential_vars) = integrator
    (; z, l, m, c_LTE, dts, tsit5tab) = cache
    if integrator.derivative_discontinuity || integrator.iter == 1
        cache.order = 1
        cache.nextorder = 1
        cache.L = 2
        cache.n_wait = 2
        z[1] = integrator.uprev
        z[2] = f(uprev, p, t) * dt
        increment_nf!(integrator.stats, 1)
        dts[1] = dt
    end
    tmp = dts[13]
    for i in 12:-1:1
        dts[i + 1] = dts[i]
    end
    dts[1] = dt
    nordsieck_needs_adjust(integrator, cache) && nordsieck_adjust!(integrator, cache)
    integrator.k[1] = z[2] / dt
    perform_predict!(cache)
    is_bdf(integrator.alg) ? calc_coeff_bdf!(cache) : calc_coeff!(cache)
    isucceed = nlsolve_functional!(integrator, cache)
    if !isucceed
        integrator.force_stepfail = true
        nordsieck_restore_failed_step!(cache, tmp)
        return nothing
    end
    if integrator.opts.adaptive
        atmp = calculate_residuals(
            cache.Δ, uprev, integrator.u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        set_EEst!(integrator, integrator.opts.internalnorm(atmp, t) * cache.c_LTE)
    end
    nordsieck_complete_step!(integrator, cache, tmp)
    integrator.k[2] = cache.z[2] / dt
    return nothing
end
count_fails()
for (nm, a) in [("JVODE_Adamsfix", JVODE_Adams()), ("JVODE_BDFfix", JVODE_BDF())], case in CASES
    Base.invokelatest(run, nm, a, case)
end
