# Reference values for tests/testthat/test-nonlinear-wrappers.R.
# Run with an environment holding NonlinearSolve 4.32.0 and the wrapped
# packages (MINPACK, NLsolve, LeastSquaresOptim,
# FastLevenbergMarquardt, NLSolvers, SpeedMapping, FixedPointAcceleration,
# SIAMFANLEquations):
#   julia --project=<env> tools/nonlinearsolve-wrappers-probe.jl
# Each line: RES problem algorithm retcode u resid stats, where stats is
# [nf, njacs, nfactors, nsolve, nsteps] or nothing.
using NonlinearSolve, LinearAlgebra
import MINPACK, NLsolve, LeastSquaresOptim, FastLevenbergMarquardt
import NLSolvers, SpeedMapping, FixedPointAcceleration, SIAMFANLEquations

circf(u, p) = [u[1]^2 + u[2]^2 - 4, u[1] * u[2] - 1]
circj(u, p) = [2u[1] 2u[2]; u[2] u[1]]
rosenf(u, p) = [10 * (u[2] - u[1]^2), 1 - u[1]]
rosenj(u, p) = [-20u[1] 10.0; -1.0 0.0]
fixpf(u, p) = [0.5 * cos(u[2]) - u[1] + 0.1, 0.5 * sin(u[1]) - u[2] + 0.2]
fixpj(u, p) = [-1.0 -0.5 * sin(u[2]); 0.5 * cos(u[1]) -1.0]
burdenf(u, p) = [3u[1] - cos(u[2] * u[3]) - 0.5,
                 u[1]^2 - 81 * (u[2] + 0.1)^2 + sin(u[3]) + 1.06,
                 exp(-u[1] * u[2]) + 20u[3] + (10pi - 3) / 3]
burdenj(u, p) = [3.0 u[3] * sin(u[2] * u[3]) u[2] * sin(u[2] * u[3]);
                 2u[1] -162 * (u[2] + 0.1) cos(u[3]);
                 -u[2] * exp(-u[1] * u[2]) -u[1] * exp(-u[1] * u[2]) 20.0]
const TT = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
const YY = [1.0, 1.6, 2.8, 4.4, 7.6, 12.0]
expf(u, p) = u[1] .* exp.(u[2] .* TT) .- YY
expj(u, p) = hcat(exp.(u[2] .* TT), u[1] .* TT .* exp.(u[2] .* TT))

nlp(f, j, u0) = NonlinearProblem(NonlinearFunction(f; jac = j), u0)
nlp0(f, u0) = NonlinearProblem(NonlinearFunction(f), u0)
nlls(f, j, u0) = NonlinearLeastSquaresProblem(
    NonlinearFunction(f; jac = j, resid_prototype = zeros(6)), u0)
nlls0(f, u0) = NonlinearLeastSquaresProblem(
    NonlinearFunction(f; resid_prototype = zeros(6)), u0)

problems = (
    ("circ", nlp(circf, circj, [3.0, 0.1])),
    ("rosen", nlp(rosenf, rosenj, [-1.2, 1.0])),
    ("fixp", nlp(fixpf, fixpj, [0.0, 0.0])),
    ("burden", nlp(burdenf, burdenj, [0.1, 0.1, -0.1])),
)
problems0 = (
    ("circ0", nlp0(circf, [3.0, 0.1])),
    ("rosen0", nlp0(rosenf, [-1.2, 1.0])),
    ("burden0", nlp0(burdenf, [0.1, 0.1, -0.1])),
)
lsproblems = (
    ("expfit", nlls(expf, expj, [1.0, 0.5])),
    ("expfit0", nlls0(expf, [1.0, 0.5])),
)

# The constructor accepts only `factorupdate = :nielson` while
# FastLevenbergMarquardt tests `:nielsen`, so the Nielsen update is never used
# through NonlinearSolve; this variant builds the struct directly to reach it.
const NIELSEN = FastLevenbergMarquardtJL(nothing, :cholesky, 1.0e-6, 13.0, 3.0, :nielsen,
    1.0e-12, 1.0e16, 1.0e-28, 1.0e32, :FastLevenbergMarquardtJL)
# LeastSquaresOptimJL with its default `autodiff = :central` fails on a problem
# with a Jacobian (select_jacobian_autodiff has no Symbol method); the `J`
# variants pass `autodiff = nothing` so the problem's Jacobian is used.
algs = (
    ("NLsolveTR", NLsolveJL(), (;)),
    ("NLsolveTRNoScale", NLsolveJL(; autoscale = false), (;)),
    ("NLsolveNewton", NLsolveJL(; method = :newton), (;)),
    ("NLsolveNewtonBT", NLsolveJL(; method = :newton, linesearch = NLsolve.LineSearches.BackTracking()), (;)),
    ("NLsolveAnderson", NLsolveJL(; method = :anderson), (;)),
    ("NLsolveAnderson2", NLsolveJL(; method = :anderson, m = 2, beta = 0.5), (;)),
    ("NLsolveBroyden", NLsolveJL(; method = :broyden), (;)),
    ("NLSolvers", NLSolversJL(NLSolvers.LineSearch(NLSolvers.Newton(), NLSolvers.Static(1))), (;)),
    ("LSOlm", LeastSquaresOptimJL(:lm), (;)),
    ("LSOdogleg", LeastSquaresOptimJL(:dogleg), (;)),
    ("LSOlmChol", LeastSquaresOptimJL(:lm; linsolve = :cholesky), (;)),
    ("LSOlmJ", LeastSquaresOptimJL(:lm; autodiff = nothing), (;)),
    ("LSOdoglegJ", LeastSquaresOptimJL(:dogleg; autodiff = nothing), (;)),
    ("LSOlmCholJ", LeastSquaresOptimJL(:lm; linsolve = :cholesky, autodiff = nothing), (;)),
    ("FastLM", FastLevenbergMarquardtJL(), (;)),
    ("FastLMqr", FastLevenbergMarquardtJL(:qr), (;)),
    ("FastLMnielson", FastLevenbergMarquardtJL(; factorupdate = :nielson), (;)),
    ("FastLMnielsen", NIELSEN, (;)),
    ("SIAMnewton", SIAMFANLEquationsJL(), (;)),
    ("SIAMptc", SIAMFANLEquationsJL(; method = :pseudotransient), (;)),
    ("SIAMnewtonGMRES", SIAMFANLEquationsJL(; linsolve = :gmres), (;)),
    ("SIAMnewtonBiCGSTAB", SIAMFANLEquationsJL(; linsolve = :bicgstab), (;)),
    # NonlinearSolve does not pass `delta` to ptcsoli, which then uses 1e-6
    ("SIAMptcGMRES", SIAMFANLEquationsJL(; method = :pseudotransient, linsolve = :gmres, delta = 1.0e-6), (;)),
    ("SIAMptcBiCGSTAB", SIAMFANLEquationsJL(; method = :pseudotransient, linsolve = :bicgstab, delta = 1.0e-6), (;)),
    ("SIAMaa0", SIAMFANLEquationsJL(; method = :anderson), (;)),
    ("SIAMaa2", SIAMFANLEquationsJL(; method = :anderson, m = 2), (;)),
    ("SpeedMapping", SpeedMappingJL(), (;)),
    ("FPAAnderson", FixedPointAccelerationJL(), (;)),
    ("FPAMPE", FixedPointAccelerationJL(; algorithm = :MPE), (;)),
    ("FPARRE", FixedPointAccelerationJL(; algorithm = :RRE), (;)),
    ("FPAVEA", FixedPointAccelerationJL(; algorithm = :VEA), (;)),
    ("FPASEA", FixedPointAccelerationJL(; algorithm = :SEA), (;)),
    ("FPASimple", FixedPointAccelerationJL(; algorithm = :Simple), (;)),
    ("FPAAitken", FixedPointAccelerationJL(; algorithm = :Aitken), (;)),
    ("FPANewton", FixedPointAccelerationJL(; algorithm = :Newton), (;)),
)
lsalgs = (
    ("LSOlm", LeastSquaresOptimJL(:lm), (;)),
    ("LSOdogleg", LeastSquaresOptimJL(:dogleg), (;)),
    ("LSOlmChol", LeastSquaresOptimJL(:lm; linsolve = :cholesky), (;)),
    ("LSOlmJ", LeastSquaresOptimJL(:lm; autodiff = nothing), (;)),
    ("LSOdoglegJ", LeastSquaresOptimJL(:dogleg; autodiff = nothing), (;)),
    ("FastLM", FastLevenbergMarquardtJL(), (;)),
    ("FastLMqr", FastLevenbergMarquardtJL(:qr), (;)),
    ("FastLMnielson", FastLevenbergMarquardtJL(; factorupdate = :nielson), (;)),
    ("FastLMnielsen", NIELSEN, (;)),
)
st(s) = s === nothing ? "nothing" : repr([s.nf, s.njacs, s.nfactors, s.nsolve, s.nsteps])
function run(pn, prob, an, alg, kw)
    sol = try
        solve(prob, alg; kw...)
    catch err
        println("ERR $pn $an ", replace(sprint(showerror, err)[1:min(end, 150)], "\n" => " "))
        return
    end
    println("RES $pn $an ", sol.retcode, " ", repr(collect(sol.u)), " ", repr(collect(sol.resid)), " ", st(sol.stats))
end
for (pn, prob) in problems, (an, alg, kw) in algs
    run(pn, prob, an, alg, kw)
end
for (pn, prob) in lsproblems, (an, alg, kw) in lsalgs
    run(pn, prob, an, alg, kw)
end
# scalar secant
run("cubic", NonlinearProblem((u, p) -> u^3 - 2u - 5, 2.0), "SIAMsecant", SIAMFANLEquationsJL(; method = :secant), (;))
run("cubicfar", NonlinearProblem((u, p) -> u^3 - 2u - 5, 10.0), "SIAMsecant", SIAMFANLEquationsJL(; method = :secant), (;))

# CMINPACK. MINPACK.jl passes closures to @cfunction, which aarch64 does not
# support, so the CMINPACK extension cannot run on this platform. The same
# cminpack routines are called here with the arguments of the MINPACK.jl
# wrappers (hybrd1, hybrj with diag = 1, mode = 2, factor = 100 and
# maxfev = typemax(Cint); lmdif1; lmder with ftol = xtol = tol, gtol = 0,
# mode = 1, factor = 100 and maxfev = 200(n + 1)) and the AlgoTrace counting.
#
# The results are not usable as a reference: the cminpack_jll 1.3.11 binary
# calls BLAS dnrm2_ from enorm without a prototype, so enorm converts the
# integer register to a double (`scvtf d0, w0` in its disassembly) and every
# norm is garbage. The ENORM line below shows it: enorm([5.01, -0.7]) should
# be 5.0586658320153940. The R port is checked instead against the MINPACK
# translation of Eigen (tools/nonlinearsolve-cminpack-reference.cpp).
const LIB = MINPACK.libcminpack
println("ENORM ", repr(ccall((:enorm, LIB), Cdouble, (Cint, Ptr{Cdouble}), 2, [5.01, -0.7])))
const CBF = Ref{Any}(nothing)
const CBJ = Ref{Any}(nothing)
const CNT = zeros(Int, 3)          # f_calls, g_calls, prev_iflag
const MAXIT = Ref(1000)
const TRACING = Ref(false)
function tpush(iflag)
    TRACING[] || return
    if iflag == 2
        CNT[3] == 1 && (CNT[2] += 1)
    elseif iflag == 1
        CNT[1] += 1
    end
    CNT[3] = iflag
end
ret() = CNT[1] > MAXIT[] ? Cint(-1) : Cint(0)
function cb_hybrd(p::Ptr{Cvoid}, n::Cint, x::Ptr{Cdouble}, fv::Ptr{Cdouble}, iflag::Cint)
    iflag < 0 && return Cint(0)
    unsafe_wrap(Array, fv, n) .= CBF[](copy(unsafe_wrap(Array, x, n)))
    tpush(iflag)
    return ret()
end
function cb_hybrj(p::Ptr{Cvoid}, n::Cint, x::Ptr{Cdouble}, fv::Ptr{Cdouble}, fj::Ptr{Cdouble}, ld::Cint, iflag::Cint)
    iflag < 0 && return Cint(0)
    xv = copy(unsafe_wrap(Array, x, n))
    iflag == 1 && (unsafe_wrap(Array, fv, n) .= CBF[](xv))
    iflag == 2 && (reshape(unsafe_wrap(Array, fj, ld * n), Int(ld), Int(n)) .= CBJ[](xv))
    tpush(iflag)
    return ret()
end
function cb_lmdif(p::Ptr{Cvoid}, m::Cint, n::Cint, x::Ptr{Cdouble}, fv::Ptr{Cdouble}, iflag::Cint)
    iflag < 0 && return Cint(0)
    unsafe_wrap(Array, fv, m) .= CBF[](copy(unsafe_wrap(Array, x, n)))
    tpush(iflag)
    return ret()
end
function cb_lmder(p::Ptr{Cvoid}, m::Cint, n::Cint, x::Ptr{Cdouble}, fv::Ptr{Cdouble}, fj::Ptr{Cdouble}, ld::Cint, iflag::Cint)
    iflag < 0 && return Cint(0)
    xv = copy(unsafe_wrap(Array, x, n))
    iflag == 1 && (unsafe_wrap(Array, fv, m) .= CBF[](xv))
    iflag == 2 && (reshape(unsafe_wrap(Array, fj, ld * n), Int(ld), Int(n)) .= CBJ[](xv))
    tpush(iflag)
    return ret()
end
function cminpack(f, j, x0, m, method; tol = eps()^(4 // 5), tracing = false, maxit = 1000)
    CBF[] = f; CBJ[] = j; fill!(CNT, 0); MAXIT[] = maxit; TRACING[] = tracing
    n = length(x0); x = copy(x0); fvec = zeros(m)
    if method == :hybr && j === nothing
        lwa = ceil(Int, (n * (3n + 13)) / 2)
        info = ccall((:hybrd1, LIB), Cint,
            (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cint),
            @cfunction(cb_hybrd, Cint, (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint)),
            C_NULL, n, x, fvec, tol, ones(lwa), lwa)
    elseif method == :hybr
        lr = ceil(Cint, n * (n + 1) / 2)
        info = ccall((:hybrj, LIB), Cint,
            (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Cint,
             Ptr{Cdouble}, Cint, Cdouble, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint,
             Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
            @cfunction(cb_hybrj, Cint, (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint)),
            C_NULL, n, x, fvec, zeros(n, n), n, tol, typemax(Cint), ones(n), 2, 100.0, 0, Cint[0], Cint[0],
            zeros(lr), lr, zeros(n), zeros(n), zeros(n), zeros(n), zeros(n))
    elseif j === nothing
        lwa = m * n + 5n + m
        info = ccall((:lmdif1, LIB), Cint,
            (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cdouble}, Cint),
            @cfunction(cb_lmdif, Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint)),
            C_NULL, m, n, x, fvec, tol, zeros(Cint, n), zeros(lwa), lwa)
    else
        info = ccall((:lmder, LIB), Cint,
            (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
             Cint, Cdouble, Cdouble, Cdouble, Cint, Ptr{Cdouble}, Cint, Cdouble,
             Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble},
             Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
            @cfunction(cb_lmder, Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint)),
            C_NULL, m, n, x, fvec, zeros(m, n), m, tol, tol, 0.0, (n + 1) * 200, zeros(n), 1, 100.0,
            0, Cint[0], Cint[0], zeros(Cint, n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(m))
    end
    rc = maximum(abs, fvec) <= tol ? "Success" : "Failure"
    return rc, x, fvec, [CNT[1], CNT[2], CNT[2], CNT[2], -1], info
end
cm(f, j) = (f === nothing ? nothing : (u -> f(u, nothing)), j === nothing ? nothing : (u -> j(u, nothing)))
mpcases = (
    ("circ", circf, circj, [3.0, 0.1], 2), ("rosen", rosenf, rosenj, [-1.2, 1.0], 2),
    ("fixp", fixpf, fixpj, [0.0, 0.0], 2), ("burden", burdenf, burdenj, [0.1, 0.1, -0.1], 3),
)
for (pn, f, j, u0, m) in mpcases, (suffix, jj) in (("", j), ("0", nothing))
    for (an, method, tracing) in (("CMINPACK", :hybr, false), ("CMINPACKTrace", :hybr, true), ("CMINPACKlm", :lm, true))
        ff, fj = cm(f, jj)
        rc, x, fv, st, info = cminpack(ff, fj, u0, m, method; tracing)
        println("RES $(pn)$(suffix) $an $rc ", repr(x), " ", repr(fv), " ", repr(st))
    end
end
for (pn, jj) in (("expfit", expj), ("expfit0", nothing))
    ff, fj = cm(expf, jj)
    rc, x, fv, st, info = cminpack(ff, fj, [1.0, 0.5], 6, :lm; tracing = true)
    println("RES $pn CMINPACK $rc ", repr(x), " ", repr(fv), " ", repr(st))
end
