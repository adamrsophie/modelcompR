# Reference values of StochasticDiffEq for tests/testthat/test-sde-composite.R,
# test-sde-baoab-iif.R and test-sde-milgeneral.R: composite algorithms with
# stiffness switching, BAOAB, the IIF methods and RKMilGeneral on fixed noise
# paths (NoiseGrid). The increments are fixed numbers z(k) = mod(37k, 101) /
# 50 - 1 scaled by sqrt(grid step) for W and Z; for nondiagonal noise the
# RKMilGeneral Z increments are the z(k) themselves (the Fourier coefficients).
# Run: julia --project=tools/julia --startup-file=no tools/stochasticdiffeq-probe.jl
using StochasticDiffEq, DiffEqNoiseProcess, LinearAlgebra
import OrdinaryDiffEqCore
const StochasticDiffEqIIF = StochasticDiffEq.StochasticDiffEqIIF
const SciMLBase = OrdinaryDiffEqCore.SciMLBase

z(k) = mod(37 * k, 101) / 50 - 1
function grid(m, hg, N; off = 0, mz = m, zscale = sqrt(hg))
    W = [zeros(m)]; Z = [zeros(mz)]
    for k in 1:N
        push!(W, W[end] .+ sqrt(hg) .* [z(off + m * (k - 1) + j) for j in 1:m])
        push!(Z, Z[end] .+ zscale .* [z(off + 7919 + mz * (k - 1) + j) for j in 1:mz])
    end
    return W, Z
end
rle(v) = (r = Tuple{Int, Int}[]; for x in v; (!isempty(r) && r[end][1] == x) ? (r[end] = (x, r[end][2] + 1)) : push!(r, (x, 1)); end; r)
flat(u) = collect(Iterators.flatten(u))

# ---- composite algorithms --------------------------------------------------
hg = 1 / 2048; N = 2150
ts = collect(0:N) .* hg
W, Z = grid(2, hg, N)
fstiff(u, p, t) = -(1 + 8000 * t * t * (1 - t) * (1 - t)) .* u
gmul(u, p, t) = 0.1 .* u
gadd(u, p, t) = [0.1, 0.2]
function comp(name, alg; g = gmul)
    prob = SDEProblem(fstiff, g, [1.0, 0.5], (0.0, 1.0), noise = NoiseGrid(ts, W, Z))
    sol = solve(prob, alg; dt = 0.01)
    println(name, " n ", length(sol.t), " acc ", sol.stats.naccept, " rej ", sol.stats.nreject,
            " choice ", repr(sol.alg_choice === nothing ? nothing : rle(sol.alg_choice)))
    println(name, " u ", repr([sol.u[100]; sol.u[300]; sol.u[end]]), " t ", repr([sol.t[100], sol.t[300]]))
end
comp("sosri2", SOSRI2())
comp("auto_dtfac1", AutoAlgSwitch(SOSRI2(), SOSRI(); maxstiffstep = 2, maxnonstiffstep = 1, dtfac = 1))
# LambaEM inside a composite: is_split_step(integrator.alg) reads the
# composite, so the step is the one of LambaEM(split = false)
comp("auto_lamba", AutoSOSRI2(LambaEM(); maxstiffstep = 1, maxnonstiffstep = 1, dtfac = 1))
comp("auto_sra", AutoSOSRA2(SOSRA(); maxstiffstep = 2, maxnonstiffstep = 1, nonstifftol = 0.2, stifftol = 0.2, dtfac = 1);
     g = gadd)
comp("choice_fn", StochasticCompositeAlgorithm((SOSRI2(), SOSRI()), integrator -> integrator.t > 0.5 ? 2 : 1))
comp("auto_dtfac2_upstream", AutoAlgSwitch(SOSRI2(), SOSRI(); maxstiffstep = 2, maxnonstiffstep = 1))

# ---- BAOAB -------------------------------------------------------------------
hb = 1 / 64
tb = collect(0:80) .* hb
Wb, _ = grid(2, hb, 80; off = 11)
f1b(v, u, p, t) = -u .- 0.1 .* u .^ 3
f2b(v, u, p, t) = v
gb(u, p, t) = 0.3 .+ 0.1 .* u
for (name, alg) in (("baoab", BAOAB(gamma = 0.7)), ("baoab_noscale", BAOAB(gamma = 0.7, scale_noise = false)))
    prob = DynamicalSDEProblem(f1b, f2b, gb, [0.5, -0.2], [1.0, 0.3], (0.0, 1.0), noise = NoiseGrid(tb, Wb))
    sol = solve(prob, alg; dt = 1 / 16)
    println(name, " n ", length(sol.t), " u ", repr(flat([collect(sol.u[i]) for i in (2, 9, 17)])))
end

# ---- IIF methods ---------------------------------------------------------------
hi = 1 / 64
ti = collect(0:80) .* hi
Wi1 = [0.0]; for k in 1:80; push!(Wi1, Wi1[end] + sqrt(hi) * z(k + 23)); end
fs1(u, p, t) = -1.3
fs2(u, p, t) = 0.4 * u + 0.1 * t
gs(u, p, t) = 0.3 * u
for (name, alg) in (("iif1m_scalar", IIF1M()), ("iif2m_scalar", IIF2M()), ("iif1mil_scalar", IIF1Mil()))
    prob = SDEProblem(SplitSDEFunction(fs1, fs2, gs), 0.8, (0.0, 1.0), noise = NoiseGrid(ti, Wi1))
    try
        sol = solve(prob, alg; dt = 1 / 16)
        println(name, " u ", repr(sol.u))
    catch e
        println(name, " error ", sprint(showerror, e))
    end
end
const Ai = [-2.0 1.0; 1.0 -2.0]
const Bi = [0.87 0.0; 0.0 0.87]
fv1(du, u, p, t) = Ai
fv2(du, u, p, t) = (du .= 1.01 .* u .+ 0.05 .* u .^ 2)
function gv(du, u, p, t)
    mul!(@view(du[:, 1]), Bi, u)
    du[:, 2] .= 0.2 .* u .^ 2
    return nothing
end
Wi, _ = grid(2, hi, 80; off = 5)
function iifvec(name, alg)
    prob = SDEProblem(SplitSDEFunction(fv1, fv2, gv), [0.3, 0.6], (0.0, 1.0), noise = NoiseGrid(ti, Wi),
                      noise_rate_prototype = zeros(2, 2))
    try
        sol = solve(prob, alg; dt = 1 / 16)
        println(name, " u ", repr(flat(sol.u)))
    catch e
        println(name, " error ", typeof(e))
    end
end
iifvec("iif1m_vec", IIF1M())
iifvec("iif2m_vec", IIF2M())
iifvec("iif1mil_vec_upstream", IIF1Mil())
# the in place IIF1Mil step with f1(u, p, t) for the old f1(t, u, du)
function OrdinaryDiffEqCore.perform_step!(integrator, cache::StochasticDiffEqIIF.IIF1MilCache)
    (; rtmp1, rtmp2, rtmp3, tmp, noise_tmp) = cache
    (; uhold, rhs, nl_rhs) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    alg = StochasticDiffEq.StochasticDiffEqCore.unwrap_alg(integrator, true)
    dW = W.dW
    sqdt = integrator.sqdt
    g = integrator.f.g
    A = integrator.f.f1(rtmp1, uprev, p, t)
    M = exp(A * dt)
    integrator.f.g(rtmp2, uprev, p, t)
    rtmp2 = M * rtmp2
    (; gtmp, gtmp2) = cache
    mul!(rtmp3, rtmp2, W.dW)
    I = zeros(length(dW), length(dW))
    mil_correction = zeros(length(dW))
    for i in 1:length(dW), j in 1:length(dW)
        I[j, i] = 0.5 * dW[i] * dW[j]
        j == i && (I[i, i] -= 0.5 * dt)
    end
    for j in 1:length(uprev)
        Kj = uprev .+ sqdt * rtmp2[:, j]
        g(gtmp, Kj, p, t)
        mul!(gtmp2, M, gtmp)
        Dgj = (gtmp2 - rtmp2) / sqdt
        mil_correction .+= Dgj * I[:, j]
    end
    rtmp3 .+= mil_correction
    mul!(tmp, M, uprev)
    tmp .+= rtmp3
    rhs.t = t
    rhs.dt = dt
    rhs.tmp = tmp
    rhs.sizeu = size(u)
    nlres = alg.nlsolve(nl_rhs, uhold)
    copyto!(uhold, nlres)
end
iifvec("iif1mil_vec_fixed", IIF1Mil())

# ---- RKMilGeneral ------------------------------------------------------------
# the method declares a second process; for diagonal noise its values are not
# read, so the grids below carry an unused Z
hr = 1 / 64
tr = collect(0:80) .* hr
Wd, Zd = grid(2, hr, 80; off = 3)
fr(u, p, t) = -u .+ 0.2 .* u .^ 2
gdiag(u, p, t) = 0.3 .* u .+ 0.1 .* u .^ 2
for (name, alg) in (("rkmg_diag_ito", RKMilGeneral()), ("rkmg_diag_strat",
        RKMilGeneral(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich)))
    prob = SDEProblem(fr, gdiag, [0.8, 0.4], (0.0, 1.0), noise = NoiseGrid(tr, Wd, Zd))
    sol = solve(prob, alg; dt = 1 / 16, adaptive = false)
    println(name, " u ", repr(flat(sol.u)))
end
hf = 1 / 1024
tf = collect(0:1100) .* hf
Wf, Zf = grid(2, hf, 1100; off = 9)
prob = SDEProblem(fr, gdiag, [0.8, 0.4], (0.0, 1.0), noise = NoiseGrid(tf, Wf, Zf))
sol = solve(prob, RKMilGeneral(); dt = 1 / 16, abstol = 1e-3, reltol = 1e-3)
println("rkmg_diag_adaptive n ", length(sol.t), " acc ", sol.stats.naccept, " rej ", sol.stats.nreject,
        " u ", repr([sol.u[5]; sol.u[end]]), " t ", repr(sol.t[1:6]))

function gnd(u, p, t)
    return [0.2 * u[1] 0.1 * u[2]; 0.05 * u[1] * u[2] 0.3 * u[2]]
end
# Z increments: the z(k) themselves, 2 m n + (m^2 + m) / 2 = 15 numbers for p = 3
Wn, Zn = grid(2, hr, 80; off = 13, mz = 15, zscale = 1.0)
for (name, alg, dt) in (("rkmg_nd_ito_coef", RKMilGeneral(p = 3), 1 / 64),
        ("rkmg_nd_strat_coef", RKMilGeneral(p = 3, interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich), 1 / 64),
        ("rkmg_nd_ito_grid", RKMilGeneral(p = 3), 1 / 16))
    prob = SDEProblem(fr, gnd, [0.8, 0.4], (0.0, 1.0), noise = NoiseGrid(tr, Wn, Zn), noise_rate_prototype = zeros(2, 2))
    sol = solve(prob, alg; dt = dt, adaptive = false)
    idx = dt == 1 / 64 ? (2, 17, 33, 65) : (2, 5, 9, 17)
    println(name, " n ", length(sol.t), " u ", repr(flat([sol.u[i] for i in idx])))
end
Wa, Za = grid(2, hf, 1100; off = 17, mz = 15, zscale = 1.0)
prob = SDEProblem(fr, gnd, [0.8, 0.4], (0.0, 1.0), noise = NoiseGrid(tf, Wa, Za), noise_rate_prototype = zeros(2, 2))
try
    sol = solve(prob, RKMilGeneral(p = 3); dt = 1 / 16, abstol = 1e-3, reltol = 1e-3)
    println("rkmg_nd_adaptive n ", length(sol.t), " acc ", sol.stats.naccept, " rej ", sol.stats.nreject,
            " u ", repr([sol.u[5]; sol.u[end]]), " t ", repr(sol.t[1:6]))
catch e
    println("rkmg_nd_adaptive error ", sprint(showerror, e))
end

# ---- composite with the noise following an increased step size -------------
# upstream shrinks the noise to the step size after the algorithm choice but
# does not grow it, so the step after a switch to the stiff method uses the
# increment of the old step size; the fix lets it follow both ways
function OrdinaryDiffEqCore.shrink_noise_to_integrator_dt!(integrator)
    W = OrdinaryDiffEqCore._get_W(integrator)
    isnothing(W) && return nothing
    if abs(integrator.dt) != abs(W.dt)
        OrdinaryDiffEqCore.reject_noise!(W, integrator.dt, integrator.u, integrator.p)
        OrdinaryDiffEqCore.reject_noise!(OrdinaryDiffEqCore._get_P(integrator), integrator.dt, integrator.u, integrator.p)
        integrator.sqdt = integrator.tdir * sqrt(abs(integrator.dt))
    end
    return nothing
end
comp("auto_dtfac2_fixed", AutoAlgSwitch(SOSRI2(), SOSRI(); maxstiffstep = 2, maxnonstiffstep = 1))
comp("auto_first_fixed", AutoAlgSwitch(SOSRI2(), SOSRI(); maxstiffstep = 2, maxnonstiffstep = 1, stiffalgfirst = true))
