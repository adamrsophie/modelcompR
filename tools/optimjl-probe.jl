# Reference results of the Optim.jl methods driven by OptimizationOptimJL
# (tests/testthat/test-optimization-optimjl.R). The analytic gradient and
# Hessian are given to both sides. Run with an environment holding
# Optimization, OptimizationOptimJL (0.4.23), Optim (2.3.2) and LineSearches
# (7.8.2):
#   julia --project=<env> tools/optimjl-probe.jl
using Optimization, OptimizationOptimJL, Optim, LineSearches, SciMLBase, ForwardDiff, ADTypes

quad(x, p) = 1.5 * x[1] * x[1] + x[1] * x[2] + x[2] * x[2] - x[1] + x[2]
quadg!(G, x, p) = (G[1] = 3 * x[1] + x[2] - 1; G[2] = x[1] + 2 * x[2] + 1; nothing)
quadh!(H, x, p) = (H[1, 1] = 3.0; H[1, 2] = 1.0; H[2, 1] = 1.0; H[2, 2] = 2.0; nothing)
rosen(x, p) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
rosg!(G, x, p) = (G[1] = -2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2); G[2] = 200 * (x[2] - x[1]^2); nothing)
function rosh!(H, x, p)
    H[1, 1] = 2 - 400 * (x[2] - x[1]^2) + 800 * x[1]^2
    H[1, 2] = -400 * x[1]
    H[2, 1] = -400 * x[1]
    H[2, 2] = 200.0
    nothing
end
qf = OptimizationFunction(quad, AutoForwardDiff(); grad = quadg!, hess = quadh!)
rf = OptimizationFunction(rosen, AutoForwardDiff(); grad = rosg!, hess = rosh!)
probs = [("quad", OptimizationProblem(qf, [0.7, -0.4], [1.0])),
         ("rosen", OptimizationProblem(rf, [-1.2, 1.0], [1.0])),
         ("rosenbox", OptimizationProblem(rf, [0.1, 0.2], [1.0]; lb = [-0.5, -0.5], ub = [0.5, 0.8]))]

methods = [
    ("GradientDescent", () -> GradientDescent()),
    ("ConjugateGradient", () -> ConjugateGradient()),
    ("AcceleratedGradientDescent", () -> AcceleratedGradientDescent()),
    ("MomentumGradientDescent", () -> MomentumGradientDescent()),
    ("BFGS", () -> Optim.BFGS()),
    ("BFGSstep", () -> Optim.BFGS(initial_stepnorm = 0.01)),
    ("LBFGS", () -> Optim.LBFGS()),
    ("LBFGS3", () -> Optim.LBFGS(m = 3)),
    ("Newton", () -> Newton()),
    ("NewtonTrustRegion", () -> NewtonTrustRegion()),
    ("NGMRES", () -> NGMRES()),
    ("OACCEL", () -> OACCEL()),
    ("GDback", () -> GradientDescent(linesearch = LineSearches.BackTracking())),
    ("BFGSback2", () -> Optim.BFGS(linesearch = LineSearches.BackTracking(order = 2))),
    ("CGquad", () -> ConjugateGradient(alphaguess = LineSearches.InitialQuadratic())),
    ("MGDcc", () -> MomentumGradientDescent(alphaguess = LineSearches.InitialConstantChange())),
    ("GDstatic", () -> GradientDescent(alphaguess = LineSearches.InitialStatic(alpha = 0.001), linesearch = LineSearches.Static())),
    ("NGMRESlbfgs", () -> NGMRES(nlprecon = Optim.LBFGS())),
    ("BFGSmt", () -> Optim.BFGS(linesearch = LineSearches.MoreThuente())),
    ("GDmt", () -> GradientDescent(linesearch = LineSearches.MoreThuente())),
    ("CGsw", () -> ConjugateGradient(linesearch = LineSearches.StrongWolfe())),
    ("LBFGSsw", () -> Optim.LBFGS(linesearch = LineSearches.StrongWolfe())),
]
boxmethods = [
    ("FminboxLBFGS", () -> Fminbox()),
    ("FminboxGD", () -> Fminbox(GradientDescent())),
    ("FminboxCG", () -> Fminbox(ConjugateGradient())),
    ("FminboxBFGS", () -> Fminbox(Optim.BFGS())),
    ("GradientDescent", () -> GradientDescent()),
]

function show_sol(pn, mn, n, sol)
    println(pn, "\t", mn, "\t", n, "\t", sol.retcode, "\t", repr(sol.u[1]), "\t", repr(sol.u[2]), "\t",
            repr(sol.objective), "\t", sol.stats.iterations, "\t", sol.stats.fevals, "\t", sol.stats.gevals,
            "\t", sol.stats.hevals)
end

for (pn, prob) in probs[1:2], (mn, mk) in methods, n in (1, 2, 3, 5, 10, 25, 0)
    sol = n == 0 ? solve(prob, mk()) : solve(prob, mk(); maxiters = n)
    show_sol(pn, mn, n, sol)
end
for (mn, mk) in boxmethods, n in (1, 2, 3, 0)
    sol = n == 0 ? solve(probs[3][2], mk()) : solve(probs[3][2], mk(); maxiters = n)
    show_sol("rosenbox", mn, n, sol)
end

# derivative free and interior point methods
for (pn, prob) in probs[1:2], n in (1, 2, 3, 5, 10, 25, 0)
    sol = n == 0 ? solve(prob, Optim.NelderMead()) : solve(prob, Optim.NelderMead(); maxiters = n)
    show_sol(pn, "NelderMead", n, sol)
end
for n in (1, 2, 3, 0)
    sol = n == 0 ? solve(probs[3][2], Fminbox(Optim.NelderMead())) :
          solve(probs[3][2], Fminbox(Optim.NelderMead()); maxiters = n)
    show_sol("rosenbox", "FminboxNM", n, sol)
end
c1!(r, x, p) = (r[1] = x[1]^2 + x[2]^2; nothing)
c1j!(J, x, p) = (J[1, 1] = 2 * x[1]; J[1, 2] = 2 * x[2]; nothing)
c1h!(H, x, p) = (H[1][1, 1] = 2.0; H[1][1, 2] = 0.0; H[1][2, 1] = 0.0; H[1][2, 2] = 2.0; nothing)
eqc!(r, x, p) = (r[1] = x[1] + x[2]; nothing)
eqj!(J, x, p) = (J[1, 1] = 1.0; J[1, 2] = 1.0; nothing)
eqh!(H, x, p) = (H[1] .= 0.0; nothing)
rfc = OptimizationFunction(rosen, AutoForwardDiff(); grad = rosg!, hess = rosh!, cons = c1!, cons_j = c1j!,
                           cons_h = c1h!)
qfe = OptimizationFunction(quad, AutoForwardDiff(); grad = quadg!, hess = quadh!, cons = eqc!, cons_j = eqj!,
                           cons_h = eqh!)
cprobs = [("rosenbox", probs[3][2]),
          ("rosencons", OptimizationProblem(rfc, [0.1, 0.2], [1.0]; lcons = [-Inf], ucons = [0.8])),
          ("quadeq", OptimizationProblem(qfe, [0.7, -0.4], [1.0]; lcons = [0.5], ucons = [0.5])),
          ("rosenboth", OptimizationProblem(rfc, [0.1, 0.2], [1.0]; lb = [-0.5, -0.5], ub = [0.5, 0.8],
                                            lcons = [0.01], ucons = [0.5]))]
for (pn, prob) in cprobs, (mn, mk) in (("IPNewton", () -> IPNewton()), ("IPNewtonmu", () -> IPNewton(μ0 = 0.1))),
    n in (1, 2, 3, 5, 10, 0)
    sol = n == 0 ? solve(prob, mk()) : solve(prob, mk(); maxiters = n)
    show_sol(pn, mn, n, sol)
end

# Optim's Adam, AdaMax and KrylovTrustRegion (whose Hessian vector product is
# the Hessian times the vector), and the first order methods on the Sphere
# and Stiefel manifolds. The Stiefel problem f(X) = tr(X'AX) + tr(B'X) on 3x2
# matrices is written with plain loops on both sides
using LinearAlgebra
quadhv!(Hv, x, v, p) = (H = zeros(2, 2); quadh!(H, x, p); mul!(Hv, H, v); nothing)
roshv!(Hv, x, v, p) = (H = zeros(2, 2); rosh!(H, x, p); mul!(Hv, H, v); nothing)
qfv = OptimizationFunction(quad, AutoForwardDiff(); grad = quadg!, hess = quadh!, hv = quadhv!)
rfv = OptimizationFunction(rosen, AutoForwardDiff(); grad = rosg!, hess = rosh!, hv = roshv!)
hvprobs = [("quad", OptimizationProblem(qfv, [0.7, -0.4], [1.0])),
           ("rosen", OptimizationProblem(rfv, [-1.2, 1.0], [1.0]))]
extmethods = [
    ("Adam", () -> Optim.Adam()),
    ("Adam2", () -> Optim.Adam(alpha = 0.01)),
    ("Adamsched", () -> Optim.Adam(alpha = t -> 0.01 / sqrt(t))),
    ("AdaMax", () -> Optim.AdaMax()),
    ("KrylovTrustRegion", () -> Optim.KrylovTrustRegion()),
]
for (pn, prob) in hvprobs, (mn, mk) in extmethods, n in (1, 2, 3, 5, 10, 25, 0)
    sol = n == 0 ? solve(prob, mk()) : solve(prob, mk(); maxiters = n)
    show_sol(pn, mn, n, sol)
end
S = Optim.Sphere()
spheremethods = [
    ("GDsphere", () -> GradientDescent(manifold = S)),
    ("CGsphere", () -> ConjugateGradient(manifold = S)),
    ("AGDsphere", () -> AcceleratedGradientDescent(manifold = S)),
    ("MGDsphere", () -> MomentumGradientDescent(manifold = S)),
    ("BFGSsphere", () -> Optim.BFGS(manifold = S)),
    ("LBFGSsphere", () -> Optim.LBFGS(manifold = S)),
    ("NGMRESsphere", () -> NGMRES(manifold = S)),
    ("OACCELsphere", () -> OACCEL(manifold = S)),
    ("Adamsphere", () -> Optim.Adam(0.01, 0.9, 0.999, 1e-8, S)),
]
for (mn, mk) in spheremethods, n in (1, 2, 3, 5, 10, 25, 0)
    sol = n == 0 ? solve(probs[2][2], mk()) : solve(probs[2][2], mk(); maxiters = n)
    show_sol("rosen", mn, n, sol)
end
const SA = [2.0 1.0 0.0; 1.0 3.0 1.0; 0.0 1.0 4.0]
const SB = [1.0 0.0; 0.0 1.0; 1.0 1.0]
function stief(X, p)
    s = 0.0
    for j in 1:2, i in 1:3
        ax = 0.0
        for k in 1:3
            ax = ax + SA[i, k] * X[k, j]
        end
        s = s + X[i, j] * ax + SB[i, j] * X[i, j]
    end
    return s
end
function stiefg!(G, X, p)
    for j in 1:2, i in 1:3
        ax = 0.0
        for k in 1:3
            ax = ax + SA[i, k] * X[k, j]
        end
        G[i, j] = 2 * ax + SB[i, j]
    end
    nothing
end
sf = OptimizationFunction(stief, AutoForwardDiff(); grad = stiefg!)
X0 = [1.0 0.0; 0.0 1.0; 0.5 0.5]
sprob = OptimizationProblem(sf, X0, [1.0])
for (mn, mk) in (("GDstiefel", () -> GradientDescent(manifold = Optim.Stiefel())),
                 ("LBFGSstiefel", () -> Optim.LBFGS(manifold = Optim.Stiefel())),
                 ("CGstiefelchol", () -> ConjugateGradient(manifold = Optim.Stiefel(:CholQR)))),
    n in (1, 2, 3, 5, 10, 25, 0)
    sol = n == 0 ? solve(sprob, mk()) : solve(sprob, mk(); maxiters = n)
    println("stiefel\t", mn, "\t", n, "\t", sol.retcode, "\t", join(repr.(vec(sol.u)), " "), "\t",
            repr(sol.objective), "\t", sol.stats.iterations, "\t", sol.stats.fevals, "\t", sol.stats.gevals,
            "\t", sol.stats.hevals)
end
# Optim's PowerManifold and ProductManifold on f(x) = x'Px + c'x with a 4x4
# matrix written with plain loops
const PA = [3.0 1.0 0.0 0.5; 1.0 2.0 0.5 0.0; 0.0 0.5 4.0 1.0; 0.5 0.0 1.0 1.0]
const PC = [1.0, -2.0, 0.5, 1.5]
function pmf(x, p)
    s = 0.0
    for i in 1:4
        ax = 0.0
        for k in 1:4
            ax = ax + PA[i, k] * x[k]
        end
        s = s + x[i] * ax + PC[i] * x[i]
    end
    return s
end
function pmg!(G, x, p)
    for i in 1:4
        ax = 0.0
        for k in 1:4
            ax = ax + PA[i, k] * x[k]
        end
        G[i] = 2 * ax + PC[i]
    end
    nothing
end
pmfun = OptimizationFunction(pmf, AutoForwardDiff(); grad = pmg!)
pmprob = OptimizationProblem(pmfun, [0.6, 0.8, -0.5, 0.2], [1.0])
PW = Optim.PowerManifold(Optim.Sphere(), (2,), (2,))
PR = Optim.ProductManifold(Optim.Sphere(), Optim.Flat(), (2,), (2,))
for (mn, mk) in (("GDpower", () -> GradientDescent(manifold = PW)),
                 ("LBFGSpower", () -> Optim.LBFGS(manifold = PW)),
                 ("CGpower", () -> ConjugateGradient(manifold = PW)),
                 ("GDproduct", () -> GradientDescent(manifold = PR)),
                 ("LBFGSproduct", () -> Optim.LBFGS(manifold = PR)),
                 ("CGproduct", () -> ConjugateGradient(manifold = PR))),
    n in (1, 2, 3, 5, 10, 25, 0)
    sol = n == 0 ? solve(pmprob, mk()) : solve(pmprob, mk(); maxiters = n)
    println("manifolds\t", mn, "\t", n, "\t", sol.retcode, "\t", join(repr.(vec(sol.u)), " "), "\t",
            repr(sol.objective), "\t", sol.stats.iterations, "\t", sol.stats.fevals, "\t", sol.stats.gevals,
            "\t", sol.stats.hevals)
end
