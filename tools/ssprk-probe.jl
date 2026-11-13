# Reference runs of the SSP methods of OrdinaryDiffEqSSPRK (the out of place
# ConstantCache forms of ssprk_perform_step.jl), for R/solvers-ssprk.R,
# R/solvers-erk-extra.R (SSPRKMSVS32, SSPRKMSVS43) and
# tests/testthat/test-ssprk.R (written by tools/ssprk-gen-test.R).
#
# Run from the repository root:
#   julia --project=tools/julia tools/ssprk-probe.jl > runs.tsv
# One line per run: method, case, saved points, nf, naccept, nreject, all
# saved states (hex, row by row), all saved times (hex) and the dense output
# at t = 0.03, 0.47, 0.71 (hex; empty for saveat runs).
using OrdinaryDiffEqSSPRK, Printf

hx(x::Real) = @sprintf("%a", Float64(x))
hx(x::AbstractVector) = join(hx.(x), ",")

# polynomial right hand sides: the elementary functions of Julia and of the C
# library may differ in the last bit
lv(u, p, t) = [1.5 * u[1] - u[1] * u[2], -3.0 * u[2] + u[1] * u[2] + 0.1 * t]
g(u, p, t) = [-20 * (u[1] - t * (2 - t)), 0.5 * u[1] * u[2]]
P1 = ODEProblem(lv, [1.0, 1.0], (0.0, 1.0))
P2 = ODEProblem(g, [0.0, 1.0], (0.0, 2.0))
tq = (0.03, 0.47, 0.71)

algs = ("SSPRK22" => SSPRK22, "SSPRK33" => SSPRK33, "SSPRK53" => SSPRK53, "SSPRK53_2N1" => SSPRK53_2N1,
        "SSPRK53_2N2" => SSPRK53_2N2, "SSPRK53_H" => SSPRK53_H, "SSPRK63" => SSPRK63, "SSPRK73" => SSPRK73,
        "SSPRK83" => SSPRK83, "SSPRK43" => SSPRK43, "SSPRK432" => SSPRK432, "SSPRKMSVS32" => SSPRKMSVS32,
        "SSPRKMSVS43" => SSPRKMSVS43, "SSPRK932" => SSPRK932, "SSPRK54" => SSPRK54, "SSPRK104" => SSPRK104,
        "KYKSSPRK42" => KYKSSPRK42, "KYK2014DGSSPRK_3S2" => KYK2014DGSSPRK_3S2)

function row(n, case, sol; dense = true)
    s = sol.stats
    d = dense ? hx(reduce(vcat, [sol(t) for t in tq])) : ""
    println(join([n, case, string(length(sol.t)), string(s.nf), string(s.naccept), string(s.nreject),
                  hx(reduce(vcat, sol.u)), hx(sol.t), d], "\t"))
end

for (n, A) in algs
    tryrun(f) = try f() catch e; println(stderr, n, ": ", e) end
    tryrun(() -> row(n, "fixed", solve(P1, A(); dt = 0.1, adaptive = false)))
    tryrun(() -> row(n, "fixed_saveat", solve(P1, A(); dt = 0.1, adaptive = false, saveat = [0.25, 0.65]); dense = false))
    if n in ("SSPRK43", "SSPRK432", "SSPRK932")
        tryrun(() -> row(n, "adaptive", solve(P1, A(); reltol = 1e-6, abstol = 1e-8)))
        tryrun(() -> row(n, "reject", solve(P2, A(); reltol = 1e-4, abstol = 1e-6, dt = 0.5)))
        tryrun(() -> row(n, "saveat", solve(P1, A(); reltol = 1e-6, abstol = 1e-8, saveat = 0.37); dense = false))
    end
end
