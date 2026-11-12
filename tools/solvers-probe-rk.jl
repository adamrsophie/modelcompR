# Butcher tableaus of the explicit Runge-Kutta methods of the OrdinaryDiffEq
# sublibraries, recovered from one step of the Julia solvers, for
# tools/solvers-gen-rk.R.
#
# Run from the repository root:
#   julia --project=tools/julia tools/solvers-probe-rk.jl OUT.tsv
# The right-hand side returns a new unit vector for every new point it is
# evaluated at, so the stage inputs of one step from u0 = 0 with dt = 1 are
# the rows of A, their times c, and the step result b. The embedded error
# estimate is read through internalnorm (abstol = 1, reltol = 0), one
# component at a time. Writes one line per method: name, then key=value
# fields (vectors comma separated, A row major).
using OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqHighOrderRK, OrdinaryDiffEqSSPRK,
      OrdinaryDiffEqLowStorageRK, OrdinaryDiffEqFeagin, OrdinaryDiffEqVerner, OrdinaryDiffEqTsit5,
      OrdinaryDiffEqSIMDRK
const OrdinaryDiffEqCore = OrdinaryDiffEqLowOrderRK.OrdinaryDiffEqCore
const SciMLBase = OrdinaryDiffEqLowOrderRK.SciMLBase
const N = 40
function probe(alg)
    pts = Dict{Vector{Float64}, Int}()
    calls = Tuple{Float64, Vector{Float64}}[]
    function f!(du, u, p, t)
        key = vcat(t, round.(u; digits = 14))
        i = get!(pts, key) do
            push!(calls, (t, copy(u)))
            length(calls)
        end
        fill!(du, 0.0)
        du[i] = 1.0
        nothing
    end
    prob = ODEProblem(f!, zeros(N), (0.0, 1.0))
    integ = init(prob, alg; dt = 1.0, adaptive = false, dense = false, save_everystep = false)
    step!(integ)
    u1 = copy(integ.u)
    s = length(calls)
    A = zeros(s, s); c = zeros(s)
    for (j, (t, u)) in enumerate(calls)
        c[j] = t
        A[j, 1:s] = u[1:s]
    end
    b = u1[1:s]
    btilde = Float64[]
    linear_est = true
    if SciMLBase.isadaptive(alg)
        # EEst = internalnorm(atmp) with atmp = dt * sum(btilde_i k_i) when
        # abstol = 1, reltol = 0; a shift keeps the norm positive
        est(w) = begin
            empty!(pts); empty!(calls)
            ig = init(prob, alg; dt = 1.0, adaptive = true, abstol = 1.0, reltol = 0.0, dtmin = 1.0,
                      force_dtmin = true, dense = false, save_everystep = false,
                      internalnorm = (u, t) -> u isa AbstractArray ? sum(w[1:length(u)] .* u) + 10.0 : abs(u))
            step!(ig)
            OrdinaryDiffEqCore.get_EEst(ig) - 10.0
        end
        btilde = [est(Float64.(1:N .== k)) for k in 1:s]
        w = collect(range(0.3, 1.7; length = N))
        linear_est = isapprox(est(w), sum(w[1:s] .* btilde); atol = 1e-10, rtol = 1e-8)
    end
    (A = A, b = b, c = c, btilde = btilde, s = s, linear_est = linear_est)
end
num(x) = repr(Float64(x))
vec_str(v) = join(num.(v), ",")
function meta(alg)
    d = Dict{String, String}()
    d["order"] = string(OrdinaryDiffEqCore.alg_order(alg))
    d["adaptive"] = string(SciMLBase.isadaptive(alg))
    if SciMLBase.isadaptive(alg)
        d["adaptive_order"] = string(try OrdinaryDiffEqCore.alg_adaptive_order(alg) catch; OrdinaryDiffEqCore.alg_order(alg) - 1 end)
    end
    d["fsal"] = string(OrdinaryDiffEqCore.isfsal(alg))
    for (k, fn) in (("qmax", OrdinaryDiffEqCore.qmax_default), ("qmin", OrdinaryDiffEqCore.qmin_default),
                    ("gamma", OrdinaryDiffEqCore.gamma_default), ("beta2", OrdinaryDiffEqCore.beta2_default))
        v = try fn(alg) catch; nothing end
        v === nothing || (d[k] = num(v))
    end
    b2 = try OrdinaryDiffEqCore.beta2_default(alg) catch; nothing end
    if b2 !== nothing
        b1 = try OrdinaryDiffEqCore.beta1_default(alg, b2) catch; nothing end
        b1 === nothing || (d["beta1"] = num(b1))
    end
    ssp = try OrdinaryDiffEqCore.ssp_coefficient(alg) catch; nothing end
    ssp === nothing || (d["ssp"] = num(ssp))
    cc = try OrdinaryDiffEqCore.default_controller(Float64, alg) catch; nothing end
    cc === nothing || (d["controller"] = string(nameof(typeof(cc))))
    if cc !== nothing && occursin("PID", string(typeof(cc)))
        d["pid_beta"] = vec_str(collect(cc.beta))
        d["pid_accept_safety"] = num(cc.accept_safety)
    end
    d
end
names_ = split(read(ARGS[2], String))
open(ARGS[1], "w") do io
    for nm in names_
        try
            alg = getfield(Main, Symbol(nm))()
            r = probe(alg)
            m = meta(alg)
            fields = ["s=$(r.s)", "A=" * vec_str(vec(permutedims(r.A))), "b=" * vec_str(r.b), "c=" * vec_str(r.c)]
            isempty(r.btilde) || push!(fields, "btilde=" * vec_str(r.btilde))
            r.linear_est || push!(fields, "nonlinear_est=true")
            for (k, v) in sort(collect(m)); push!(fields, "$k=$v"); end
            println(io, nm, "\t", join(fields, "\t"))
            flush(io)
            println(stderr, "OK ", nm)
            flush(stderr)
        catch e
            println(stderr, "FAIL ", nm, ": ", first(sprint(showerror, e), 300))
            flush(stderr)
        end
    end
end
