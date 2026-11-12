# Reference values for the polyalgorithms, multistart, Jacobian reuse,
# tracing, Sobol points, pickchunksize, sparsity detection and colorings
# reexported by NonlinearSolve (NonlinearSolve 4.32.0, NonlinearSolveBase
# 2.55.1, NonlinearSolveFirstOrder 2.10.0, Sobol 1.5.0,
# SparseConnectivityTracer 1.2.4, SparseMatrixColorings 0.4.29).
# Run with an environment holding these packages.
using NonlinearSolve, LinearAlgebra, Sobol, SparseArrays, ADTypes
using SparseConnectivityTracer, SparseMatrixColorings

fmt(v) = "[" * join(repr.(v), ", ") * "]"
st(s) = (s.stats.nf, s.stats.njacs, s.stats.nfactors, s.stats.nsolve, s.stats.nsteps)
report(tag, s) = println("CASE ", tag, " ", s.retcode, " ", fmt(s.u), " ", join(st(s), " "))

rosen(u, p) = [10 * (u[2] - u[1] * u[1]), 1 - u[1]]
rosenj(u, p) = [-20*u[1] 10.0; -1.0 0.0]
s5 = sqrt(5.0); s10 = sqrt(10.0)
powell(u, p) = [u[1] + 10*u[2], s5*(u[3] - u[4]), (u[2] - 2*u[3])*(u[2] - 2*u[3]), s10*(u[1] - u[4])*(u[1] - u[4])]
powellj(u, p) = [1.0 10.0 0.0 0.0;
                 0.0 0.0 s5 -s5;
                 0.0 2*(u[2]-2*u[3]) -4*(u[2]-2*u[3]) 0.0;
                 2*s10*(u[1]-u[4]) 0.0 0.0 -2*s10*(u[1]-u[4])]
tt = [1.0, 2.0, 4.0, 8.0, 16.0]; yy = [0.9, 1.5, 2.1, 2.6, 2.75]
mm(u, p) = [u[1]*t/(u[2] + t) - y for (t, y) in zip(tt, yy)]
mmj(u, p) = reduce(vcat, [[t/(u[2]+t) -u[1]*t/((u[2]+t)*(u[2]+t))] for t in tt])
circ(u, p) = [u[1]*u[1] + u[2]*u[2] - 4, u[1]*u[2] - 1]
circj(u, p) = [2*u[1] 2*u[2]; u[2] u[1]]
# Broyden tridiagonal: (3 - 2x_i) x_i - x_(i-1) - 2 x_(i+1) + 1
function broyden(u, p)
    n = length(u)
    [(3 - 2*u[i])*u[i] - (i > 1 ? u[i-1] : 0.0) - 2*(i < n ? u[i+1] : 0.0) + 1 for i in 1:n]
end
function broydenj(u, p)
    n = length(u)
    J = zeros(n, n)
    for i in 1:n
        J[i, i] = 3 - 4*u[i]
        i > 1 && (J[i, i-1] = -1.0)
        i < n && (J[i, i+1] = -2.0)
    end
    J
end
# a root outside the box: the bounded solve stalls on the bound
outside(u, p) = [u[1] - 3, u[2] - 0.5]
outsidej(u, p) = [1.0 0.0; 0.0 1.0]
# a root inside the box that a start on the bound misses: x^2 - 1 with x in [-2, 2]
# from u0 = 0 (zero gradient)
sq(u, p) = [u[1]*u[1] - 1, u[2] - 1]
sqj(u, p) = [2*u[1] 0.0; 0.0 1.0]

nlp(f, j, u0, lb, ub) = NonlinearProblem(NonlinearFunction(f; jac = j), u0; lb = lb, ub = ub)
nlls(f, j, u0, lb, ub) = NonlinearLeastSquaresProblem(NonlinearFunction(f; jac = j, resid_prototype = zeros(5)), u0; lb = lb, ub = ub)

cases = [
    ("rosen_in", nlp(rosen, rosenj, [-1.2, 1.0], [-2.0, -2.0], [2.0, 2.0])),
    ("rosen_bd", nlp(rosen, rosenj, [-1.2, 1.0], [-2.0, -2.0], [0.5, 2.0])),
    ("powell", nlp(powell, powellj, [3.0, -1.0, 0.0, 1.0], fill(-5.0, 4), fill(5.0, 4))),
    ("circ", nlp(circ, circj, [3.0, 0.5], [0.0, 0.0], [3.0, 3.0])),
    ("mm_ls", nlls(mm, mmj, [1.0, 1.0], [0.0, 0.0], [2.5, 10.0])),
    ("outside", nlp(outside, outsidej, [1.0, 1.0], [0.0, 0.0], [2.0, 2.0])),
    ("sq", nlp(sq, sqj, [0.0, 0.0], [-2.0, -2.0], [2.0, 2.0])),
    ("rosen_half", nlp(rosen, rosenj, [-1.2, 1.0], [-2.0, -Inf], [2.0, Inf])),
]

for (name, prob) in cases
    report("$name FSBP", solve(prob, FastShortcutBoundedPolyalg()))
    report("$name FSBP_post", solve(prob, FastShortcutBoundedPolyalg(; must_support_postcondition = true)))
    report("$name Default", solve(prob))
    report("$name Sobol4", solve(prob, SobolMultistart(; nstarts = 4)))
    report("$name SobolAll", solve(prob, SobolMultistart(; nstarts = 6, early_exit = false)))
    report("$name SobolNoRest", solve(prob, SobolMultistart(; nstarts = 3, restoration = false)))
    report("$name Poly_DT", solve(prob, NonlinearSolvePolyAlgorithm((Dogbox(), TrustRegionReflective()))))
    report("$name Poly_BLMTR", solve(prob, NonlinearSolvePolyAlgorithm((BoundedLevenbergMarquardt(), BoundedTrustRegion()); start_index = 2)))
end

# Jacobian reuse of BoundedTrustRegion
for n in (8, 16, 20)
    u0 = fill(-1.0, n)
    prob = nlp(broyden, broydenj, u0, fill(-2.0, n), fill(2.0, n))
    report("broyden$n BTR", solve(prob, BoundedTrustRegion()))
    prob2 = NonlinearProblem(NonlinearFunction(broyden; jac = broydenj), u0)
    for (tag, reuse) in (("auto", nothing), ("off", false), ("on", true), ("age3", JacobianReuse(; max_age = 3)),
                         ("ratio", JacobianReuse(; max_age = 10, max_residual_ratio = 0.5)))
        alg = GeneralizedFirstOrderAlgorithm(; trustregion = NonlinearSolveFirstOrder.BoundedTrustRegionScheme(),
            descent = Dogleg(), jacobian_reuse = reuse, name = :BoundedTrustRegion)
        report("broyden$n GFOA_$tag", solve(prob2, alg))
    end
end

# traces
let prob = cases[1][2]
    for (tag, alg) in (("Dogbox", Dogbox()), ("BTR", BoundedTrustRegion()), ("TRR", TrustRegionReflective()))
        for (lt, lvl) in (("min", TraceMinimal()), ("cond", TraceWithJacobianConditionNumber()), ("min2", TraceMinimal(2)))
            sol = solve(prob, alg; store_trace = Val(true), trace_level = lvl)
            for e in sol.trace.history
                println("TRACE $tag $lt ", e.iteration, " ", repr(e.fnorm), " ", repr(e.stepnorm), " ",
                        e.condJ === nothing ? "nothing" : repr(e.condJ))
            end
        end
    end
    sol = solve(cases[5][2], Dogbox(); store_trace = Val(true), trace_level = TraceMinimal())
    for e in sol.trace.history
        println("TRACE mmDogbox min ", e.iteration, " ", repr(e.fnorm), " ", repr(e.stepnorm), " ", e.norm_type)
    end
    solve(prob, Dogbox(); show_trace = Val(true), trace_level = TraceWithJacobianConditionNumber(3))
end

# Sobol points
for d in (1, 2, 5)
    s = SobolSeq(d)
    for k in 1:12
        println("SOBOL $d $k ", fmt(next!(s)))
    end
end
let s = SobolSeq([-1.0, 0.0, 2.0], [1.0, 10.0, 3.0])
    for k in 1:6
        println("SOBOLS $k ", fmt(next!(s)))
    end
end
let s = SobolSeq(40)
    x = zeros(40)
    for k in 1:40
        next!(s, x)
    end
    println("SOBOL40 ", fmt(x))
end

# pickchunksize
println("CHUNK ", join([NonlinearSolve.pickchunksize(n) for n in (1, 5, 12, 13, 24, 25, 30, 100, 1000)], " "))
println("CHUNKV ", NonlinearSolve.pickchunksize(zeros(30)))

# sparsity
spf1(x) = [x[1] * x[2], sin(x[3]), x[1] + x[4] * x[4], abs(x[2]) + floor(x[3]), x[4] / x[2]]
spf2(x) = [sum(x[1:3]), prod(x[2:4]), max(x[1], x[4]), x[3]^2, exp(x[1]) * x[4]]
hsf1(x) = x[1] * x[2] + exp(x[3]) + x[4] / x[1] + x[2]^3
hsf2(x) = sum(abs2, x) + x[1] * x[3] * x[5] + sin(x[2] + x[4])
function prpat(tag, S)
    I, J, _ = findnz(sparse(S))
    println("SPAR $tag ", size(S, 1), " ", size(S, 2), " ", join(I, ","), " ", join(J, ","))
end
prpat("j1", jacobian_sparsity(spf1, rand(4), TracerSparsityDetector()))
prpat("j2", jacobian_sparsity(spf2, rand(4), TracerSparsityDetector()))
prpat("h1", hessian_sparsity(hsf1, rand(4), TracerSparsityDetector()))
prpat("h2", hessian_sparsity(hsf2, rand(5), TracerSparsityDetector()))

# colorings on deterministic patterns
function pat(m, n, a, b, c)
    S = spzeros(Bool, m, n)
    for i in 1:m, j in 1:n
        if (a * i + b * j) % c == 0 || i == j
            S[i, j] = true
        end
    end
    S
end
for (m, n, a, b, c) in ((6, 6, 3, 5, 7), (8, 5, 2, 3, 5), (5, 9, 4, 1, 6), (10, 10, 1, 1, 4))
    S = pat(m, n, a, b, c)
    println("COLOR col $m $n $a $b $c ", join(ADTypes.column_coloring(S, GreedyColoringAlgorithm()), " "))
    println("COLOR row $m $n $a $b $c ", join(ADTypes.row_coloring(S, GreedyColoringAlgorithm()), " "))
    if m == n
        Ssym = S .| S'
        println("COLOR sym $m $n $a $b $c ", join(ADTypes.symmetric_coloring(Ssym, GreedyColoringAlgorithm()), " "))
    end
end
let n = 12
    S = spzeros(Bool, n, n)
    for i in 1:n
        S[i, i] = true
        i < n && (S[i, i+1] = S[i+1, i] = true)
        i + 3 <= n && (S[i, i+3] = S[i+3, i] = true)
    end
    println("COLOR sym banded12 ", join(ADTypes.symmetric_coloring(S, GreedyColoringAlgorithm()), " "))
    A = spzeros(Bool, n, n)
    for i in 1:n
        A[i, i] = true; A[1, i] = A[i, 1] = true
    end
    println("COLOR sym arrow12 ", join(ADTypes.symmetric_coloring(A, GreedyColoringAlgorithm()), " "))
end
