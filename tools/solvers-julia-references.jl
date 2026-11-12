# Julia references for the solver tests (tests/testthat/test-solvers.R,
# test-delay.R). Run from the repository root:
#   julia --startup-file=no --project=tools/julia tools/solvers-julia-references.jl
# with the tools/julia environment (OrdinaryDiffEqBDF 2.4.14, DelayDiffEq
# 6.4.6, OrdinaryDiffEqCore 4.19.1; the R port follows OrdinaryDiffEq.jl commit
# 4460ff8). Writes
#   tests/testthat/fixtures/julia/solver_stats.csv   step statistics of FBDF,
#     DFBDF and MethodOfSteps(Tsit5()) (in place functions, as MTK builds them)
#   tests/testthat/fixtures/julia/stald_detect.csv   verdicts of
#     OrdinaryDiffEqBDF._stald_detect on random scaled derivative data, plus
#     rows recorded from an FBDF run of an oscillator (u1' = -u1 + 100 u2,
#     u2' = -100 u1 - u2, u3' = -0.05 u3 + 0.01 u1, reltol 1e-6, [0, 100]) where
#     the detection fires; those rows are kept from the existing file and
#     checked against Julia here.
using OrdinaryDiffEq, DelayDiffEq
const OBDF = only(m for m in values(Base.loaded_modules) if nameof(m) == :OrdinaryDiffEqBDF)
const DFBDF = OBDF.DFBDF
function rober!(du, u, p, t)
    du[1] = -0.04u[1] + 1e4 * u[2] * u[3]; du[2] = 0.04u[1] - 1e4 * u[2] * u[3] - 3e7 * u[2]^2; du[3] = 3e7 * u[2]^2
end
function vdp!(du, u, p, t)
    du[1] = u[2]; du[2] = p * (1 - u[1]^2) * u[2] - u[1]
end
function roberdae!(r, du, u, p, t)
    r[1] = -0.04u[1] + 1e4 * u[2] * u[3] - du[1]
    r[2] = 0.04u[1] - 3e7 * u[2]^2 - 1e4 * u[2] * u[3] - du[2]
    r[3] = u[1] + u[2] + u[3] - 1
end
function vdpdae!(r, du, u, p, t)
    r[1] = u[2] - du[1]; r[2] = p * (1 - u[1]^2) * u[2] - u[1] - du[2]
end
rows = String[]
add(name, sol) = push!(rows, join([name, string(sol.retcode), sol.stats.naccept, sol.stats.nreject, sol.stats.njacs, sol.stats.nw], ","))
rp = ODEProblem(rober!, [1.0, 0, 0], (0.0, 1e5))
add("rober_FBDF_1e-6", solve(rp, FBDF(), reltol = 1e-6, abstol = 1e-10))
add("rober_FBDF_tf_1e-6", solve(rp, FBDF(time_filter = true), reltol = 1e-6, abstol = 1e-10))
add("rober_FBDF_qmax2_1e-6", solve(rp, FBDF(qmax = 2 // 1), reltol = 1e-6, abstol = 1e-10))
vp = ODEProblem(vdp!, [2.0, 0], (0.0, 1000.0), 1000.0)
add("vdp_FBDF_1e-5", solve(vp, FBDF(), reltol = 1e-5, abstol = 1e-7))
add("vdp_FBDF_tf_1e-5", solve(vp, FBDF(time_filter = true), reltol = 1e-5, abstol = 1e-7))
dp = DAEProblem(roberdae!, [-0.04, 0.04, 0.0], [1.0, 0, 0], (0.0, 1e5), differential_vars = [true, true, false])
add("roberdae_DFBDF_1e-6", solve(dp, DFBDF(), reltol = 1e-6, abstol = 1e-10))
vdp_d = DAEProblem(vdpdae!, [0.0, -2.0], [2.0, 0], (0.0, 1000.0), 1000.0, differential_vars = [true, true])
add("vdpdae_DFBDF_1e-5", solve(vdp_d, DFBDF(), reltol = 1e-5, abstol = 1e-7))
h1(p, t; idxs = nothing) = [1.0]
f1(u, h, p, t) = [-h(p, t - 1)[1]]
add("dde_lag1_Tsit5", solve(DDEProblem(f1, [1.0], h1, (0.0, 8.0); constant_lags = [1.0]), MethodOfSteps(Tsit5()),
                            reltol = 1e-8, abstol = 1e-10))
f2(u, h, p, t) = [-h(p, t - 1)[1] - 0.5 * h(p, t - 1.5)[1]]
add("dde_lags2_Tsit5", solve(DDEProblem(f2, [1.0], h1, (0.0, 6.0); constant_lags = [1.0, 1.5]), MethodOfSteps(Tsit5()),
                             reltol = 1e-8, abstol = 1e-10))
f3(u, h, p, t) = [-h(p, t - 1 - u[1])[1]]
add("dde_sd_Tsit5", solve(DDEProblem(f3, [1.0], h1, (0.0, 2.0); dependent_lags = [(u, p, t) -> 1 + u[1]]),
                          MethodOfSteps(Tsit5()), reltol = 1e-9, abstol = 1e-11))
open("tests/testthat/fixtures/julia/solver_stats.csv", "w") do io
    println(io, "case,retcode,naccept,nreject,njacs,nw")
    foreach(r -> println(io, r), rows)
end

# STALD
using Random, Printf
const STALD_FILE = "tests/testthat/fixtures/julia/stald_detect.csv"
# the recorded rows are the ones that order a reduction (flags 4 to 6)
recorded = isfile(STALD_FILE) ? filter(l -> parse(Int, split(l, ",")[end]) > 3, readlines(STALD_FILE)[2:end]) :
           String[]
Random.seed!(20261111)
srows = String[]
seen = Dict{Int, Int}()
for trial in 1:200000
    q = rand(3:5)
    r = 0.85 + 0.3 * rand()
    st = OBDF.StabilityLimitDetectionState(Float64)
    mode = rand(1:3)
    for k in 1:3, i in 1:5
        a = 10.0^(-3 * rand())
        b = 10.0^(-3 * rand())
        st.ssdat[i, k] = mode == 1 ? (1 + 0.001randn()) * b * r^(2 * (5 - i)) :
                         mode == 2 ? a + b * r^(2 * (5 - i)) : a * (1 + 0.5 * rand())
    end
    if mode == 2
        for k in 1:3
            a = 10.0^(-3 * rand()); b = 10.0^(-2 * rand())
            for i in 1:5
                st.ssdat[i, k] = a + b * r^(2 * (5 - i)) * (1 + 1e-4 * randn())
            end
        end
    end
    flag = OBDF._stald_detect(st, q)
    c = get(seen, flag, 0)
    if c < 6
        seen[flag] = c + 1
        push!(srows, join([string(q); [@sprintf("%.17g", x) for x in vec(st.ssdat)]; string(flag)], ","))
    end
end
for l in recorded
    v = parse.(Float64, split(l, ","))
    st = OBDF.StabilityLimitDetectionState(Float64)
    st.ssdat .= reshape(v[2:16], 5, 3)
    @assert OBDF._stald_detect(st, Int(v[1])) == Int(v[17])
    push!(srows, l)
end
open(STALD_FILE, "w") do io
    println(io, "q," * join(["s$(i)_$(k)" for k in 1:3 for i in 1:5], ",") * ",flag")
    foreach(r -> println(io, r), srows)
end
