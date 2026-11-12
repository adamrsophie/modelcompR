# Reference values for tests/testthat/test-nonlinear-homotopy-extras.R.
# Run with an environment holding NonlinearSolve 4.32.0 (NonlinearSolveBase
# 2.55.1, SimpleNonlinearSolve with SimpleHomotopySweep):
#   julia --project=<env> tools/nonlinearsolve-homotopy-probe.jl
# Every residual call records its lambda; the line LAM lists the distinct
# consecutive lambda values (the continuation path), RES the retcode and u.
using NonlinearSolve
const SHS = NonlinearSolve.SimpleNonlinearSolve.SimpleHomotopySweep

const LAMS = Float64[]
rec(l) = l isa Float64 && (isempty(LAMS) || LAMS[end] != l) && push!(LAMS, Float64(l))

atanf(u, p, l) = (rec(l); [l * atan(u[1] - 3) + (1 - l) * u[1]])
atanj(u, p, l) = reshape([l / (1 + (u[1] - 3)^2) + (1 - l)], 1, 1)
scurf(u, p, l) = (rec(l); [u[1]^3 - 3 * u[1] - (8 * l - 4)])
scurj(u, p, l) = reshape([3 * u[1]^2 - 3], 1, 1)
circf(u, p, l) = (rec(l); [l * (u[1]^2 + u[2]^2 - 4) + (1 - l) * (u[1] - 2),
                           l * (u[1] * u[2] - 1) + (1 - l) * (u[2] - 0.5)])
circj(u, p, l) = [l * 2 * u[1] + (1 - l) l * 2 * u[2]; l * u[2] l * u[1] + (1 - l)]

problems = (
    ("atan", atanf, atanj, [12.0]),
    ("scurve", scurf, scurj, [-2.2]),
    ("circ", circf, circj, [2.0, 0.5]),
)
algs = (
    ("SimpleSweep", SHS()),
    ("SimpleSweepConst", SHS(; predictor = :constant)),
    ("SimpleSweepFixed", SHS(; adaptive = false, nsteps = 4)),
    ("KantNR", KantorovichHomotopy(; inner = NewtonRaphson())),
    ("KantNRSecant", KantorovichHomotopy(; inner = NewtonRaphson(), predictor = :secant)),
    ("Kant", KantorovichHomotopy()),
    ("ArcLength", ArcLengthContinuation()),
    ("ArcLengthTangent", ArcLengthContinuation(; predictor = :tangent)),
    ("ArcLengthNR", ArcLengthContinuation(; inner = NewtonRaphson())),
    ("Poly", HomotopyPolyAlgorithm()),
    ("PolyNoHandoff", HomotopyPolyAlgorithm(; warm_handoff = false)),
    ("FSHP", FastShortcutHomotopyPolyalg()),
)
for (pn, f, j, u0) in problems
    nf = NonlinearFunction(f; jac = j, lambda_extended = true)
    prob = HomotopyProblem(nf, u0)
    for (an, alg) in algs
        empty!(LAMS)
        sol = try
            solve(prob, alg)
        catch err
            println("ERR $pn $an ", sprint(showerror, err)[1:min(end, 200)])
            continue
        end
        println("RES $pn $an ", sol.retcode, " ", repr(sol.u))
        println("LAM $pn $an ", repr(LAMS))
    end
end
