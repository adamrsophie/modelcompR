# The derivatives OptimizationBase builds from ForwardDiff for
# AutoForwardDiff() (R/forwarddiff.R): the gradient, the constraint
# Jacobian, the Hessian (DifferentiationInterface's SecondOrder(AutoForwardDiff,
# AutoForwardDiff), a pushforward of the ForwardDiff gradient per seed), the
# constraint Hessians and the Hessian of the Lagrangian (matrix form, which
# OptimizationMadNLP calls, and vector form, which OptimizationIpopt calls)
# at a few points, printed with all digits. Usage (packages: Optimization
# OptimizationBase ADTypes ForwardDiff DifferentiationInterface):
#   julia --project=tools/julia tools/forwarddiff-probe.jl > forwarddiff-probe.txt
import Optimization
import OptimizationBase
import ADTypes
import ForwardDiff
import DifferentiationInterface

println("== versions")
for p in ("Optimization", "OptimizationBase", "ForwardDiff", "DifferentiationInterface")
    try
        println(p, " ", pkgversion(Base.require(Main, Symbol(p))))
    catch
        println(p, " ?")
    end
end
flush(stdout)

hs071_f(x, p) = x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
hs071_c(res, x, p) = (res .= [x[1] * x[2] * x[3] * x[4], sum(abs2, x)])
circle_f(x, p) = x[1]^2 + x[2]^2
circle_c(res, x, p) = (res .= [x[1]^2 + x[2]^2, x[2] - x[1]^2])
wb_f(x, p) = x[1]
wb_c(res, x, p) = (res .= [x[1]^2 - x[2] - 1, x[1] - x[3] - 0.5])
rosen_f(x, p) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
# the elementary functions of R/forwarddiff.R (the R side writes the same
# expressions in the same order)
fun_f(x, p) = exp(x[1]) * sin(x[2]) + log(x[3]) / sqrt(x[1]) + tan(x[2]) * atan(x[3]) -
    tanh(x[1] * x[2]) + x[3]^2.5 + abs(x[1] - x[3]) + x[2]^3 / (1 + x[1]^2) + cos(x[1] * x[3]) +
    x[1]^x[2] + 2.0^x[3] + inv(x[2]) + asin(x[1] / 4) + acos(x[2] / 4) + sinh(x[3] / 3) +
    cosh(x[1] / 3) + log10(x[2]) + log2(x[3]) + exp2(x[1]) + cbrt(x[2]) + x[1]^4
fun_c(res, x, p) = (res .= [expm1(x[1]) * log1p(x[2]), asinh(x[3]) + x[1] / x[2],
    acosh(x[2] + 1) - atanh(x[1] / 4) * x[3]])

fmt(v) = repr(Float64.(v))

function run(name, f, c, m, points, lams)
    for (k, x) in enumerate(points)
        println("== ", name, " point ", k)
        flush(stdout)
        try
            of = Optimization.OptimizationFunction(f, ADTypes.AutoForwardDiff(); cons = c)
            fi = OptimizationBase.instantiate_function(of, OptimizationBase.ReInitCache(x, nothing),
                ADTypes.AutoForwardDiff(), m; g = true, h = true, cons_j = true, cons_h = true,
                lag_h = m > 0)
            n = length(x)
            G = zeros(n); fi.grad(G, x)
            println("x ", fmt(x))
            println("grad ", fmt(G))
            H = zeros(n, n); fi.hess(H, x)
            println("hess ", fmt(vec(H)))
            if m > 0
                J = zeros(m, n); fi.cons_j(J, x)
                println("cons_j ", fmt(vec(J)))
                Hs = [zeros(n, n) for _ in 1:m]; fi.cons_h(Hs, x)
                for i in 1:m
                    println("cons_h", i, " ", fmt(vec(Hs[i])))
                end
                for (j, (σ, λ)) in enumerate(lams)
                    L = zeros(n, n); fi.lag_h(L, x, σ, λ)
                    println("lag_h", j, " ", fmt(vec(L)))
                    h = zeros(div(n * (n + 1), 2)); fi.lag_h(h, x, σ, λ)
                    println("lag_hv", j, " ", fmt(h))
                end
            end
        catch err
            println("error ", sprint(showerror, err))
        end
        flush(stdout)
    end
end

run("hs071", hs071_f, hs071_c, 2, ([1.0, 5.0, 5.0, 1.0], [1.3, 4.1, 3.7, 1.9]),
    ((1.0, [0.5, -0.3]), (0.0, [0.7, 0.0]), (0.25, [1.5, 2.5])))
run("circle", circle_f, circle_c, 2, ([0.1, 0.2], [0.7, 1.3]), ((1.0, [0.5, -0.3]), (0.0, [0.0, 2.0])))
run("wachter_biegler", wb_f, wb_c, 2, ([-2.0, 1.0, 1.0],), ((1.0, [0.5, -0.3]), (0.0, [1.0, 1.0])))
run("rosenbrock", rosen_f, nothing, 0, ([-1.2, 1.0], [0.3, 0.7]), ())
run("functions", fun_f, fun_c, 3, ([0.7, 1.3, 2.1], [1.1, 0.4, 0.9]),
    ((1.0, [0.5, -0.3, 0.2]), (0.0, [0.7, 0.0, -1.1]), (0.3, [1.5, 2.5, -0.5])))
# 13 unknowns: two chunks of 7 (ForwardDiff's pickchunksize)
long_f(x, p) = sum(i * x[i]^2 * x[mod1(i + 1, 13)] for i in 1:13) + exp(x[1] * x[13])
long_c(res, x, p) = (res .= [sum(x) * x[2], x[7]^3 - x[12] * x[3]])
run("long", long_f, long_c, 2, ([0.1 * i for i in 1:13],), ((1.0, [0.5, -0.3]), (0.0, [0.2, 0.0])))
