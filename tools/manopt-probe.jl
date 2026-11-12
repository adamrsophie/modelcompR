# Reference values for tests/testthat/test-optimization-manopt.R.
# Manopt 0.6.12, Manifolds 0.11.38, ManifoldsBase 2.6.5.
# Run: julia --project=<env with Manopt and Manifolds> --startup-file=no tools/manopt-probe.jl
using Manopt, Manifolds, ManifoldsBase, LinearAlgebra

# a cubic on the sphere S^2, written as scalar expressions so that R evaluates
# it with the same floating point operations
f(M, p) = 2 * p[1]^2 + p[2]^2 + 3 * p[3]^2 + p[1] * p[2] + 0.2 * p[1] * p[3] + 0.6 * p[2] * p[3] + p[1] * p[2] * p[3]
egrad(p) = [4 * p[1] + p[2] + 0.2 * p[3] + p[2] * p[3],
            2 * p[2] + p[1] + 0.6 * p[3] + p[1] * p[3],
            6 * p[3] + 0.2 * p[1] + 0.6 * p[2] + p[1] * p[2]]
ehess(p, v) = [4 * v[1] + (1 + p[3]) * v[2] + (0.2 + p[2]) * v[3],
               (1 + p[3]) * v[1] + 2 * v[2] + (0.6 + p[1]) * v[3],
               (0.2 + p[2]) * v[1] + (0.6 + p[1]) * v[2] + 6 * v[3]]
M = Sphere(2)
grad_f(M, p) = project(M, p, egrad(p))
function Hess_f(M, p, X)
    Y = project(M, p, ehess(p, X))
    G = egrad(p)
    return Y + (-dot(p, G - project(M, p, G))) * X
end
p0 = [0.6, -0.48, 0.64]

pr(name, x) = println(name, " = ", repr(x))
iters(st) = Manopt.get_count(st, :Iterations)

st = gradient_descent(M, f, grad_f, copy(p0); stopping_criterion = StopAfterIteration(25), return_state = true)
pr("gd_exp", get_solver_result(st))
st = gradient_descent(M, f, grad_f, copy(p0); stopping_criterion = StopAfterIteration(25),
                      retraction_method = ProjectionRetraction(), return_state = true)
pr("gd_proj", get_solver_result(st))
st = conjugate_gradient_descent(M, f, grad_f, copy(p0); stopping_criterion = StopAfterIteration(25), return_state = true)
pr("cg", get_solver_result(st))
st = quasi_Newton(M, f, grad_f, copy(p0); stopping_criterion = StopAfterIteration(12), return_state = true)
pr("qn", get_solver_result(st))
st = trust_regions(M, f, grad_f, Hess_f, copy(p0); stopping_criterion = StopAfterIteration(10), return_state = true)
pr("tr", get_solver_result(st))
st = adaptive_regularization_with_cubics(M, f, grad_f, Hess_f, copy(p0); stopping_criterion = StopAfterIteration(8),
                                         return_state = true)
pr("arc", get_solver_result(st))
st = Frank_Wolfe_method(M, f, grad_f, copy(p0); stopping_criterion = StopAfterIteration(6), return_state = true)
pr("fw", get_solver_result(st))
simplex = [[0.6, -0.48, 0.64], [0.0, 0.6, 0.8], [-0.48, 0.64, 0.6]]
st = NelderMead(M, f, NelderMeadSimplex(simplex); stopping_criterion = StopAfterIteration(30), return_state = true)
pr("nm", get_solver_result(st))
# full default runs: the minimizer and the iterations
st = gradient_descent(M, f, grad_f, copy(p0); return_state = true)
pr("gd_default", (get_solver_result(st), Manopt.get_stopping_criterion(st).criteria[1].at_iteration,
                  Manopt.get_stopping_criterion(st).criteria[2].at_iteration))
st = quasi_Newton(M, f, grad_f, copy(p0); return_state = true)
pr("qn_default", (get_solver_result(st), Manopt.get_stopping_criterion(st).criteria[2].at_iteration))
st = trust_regions(M, f, grad_f, Hess_f, copy(p0); return_state = true)
pr("tr_default", (get_solver_result(st), Manopt.get_stopping_criterion(st).criteria[2].at_iteration))

# Rosenbrock on Euclidean(2), as the OptimizationManopt tests (x0 = 0, p = [1, 100])
R2 = Euclidean(2)
pp = [1.0, 100.0]
rf(M, x) = (pp[1] - x[1])^2 + pp[2] * (x[2] - x[1]^2)^2
rg(M, x) = [-2.0 * (pp[1] - x[1]) - 4.0 * pp[2] * (x[2] - x[1]^2) * x[1], 2.0 * pp[2] * (x[2] - x[1]^2)]
rh(M, x, v) = [(2 - 4 * pp[2] * (x[2] - 3 * x[1]^2)) * v[1] + (-4 * pp[2] * x[1]) * v[2],
               (-4 * pp[2] * x[1]) * v[1] + (2 * pp[2]) * v[2]]
x0 = [0.0, 0.0]
st = gradient_descent(R2, rf, rg, copy(x0); stopping_criterion = StopAfterIteration(50), return_state = true)
pr("rosen_gd50", get_solver_result(st))
st = quasi_Newton(R2, rf, rg, copy(x0); stopping_criterion = StopAfterIteration(20), return_state = true)
pr("rosen_qn20", get_solver_result(st))
# the solver defaults as OptimizationManopt calls them without keywords
st = adaptive_regularization_with_cubics(R2, rf, rg, rh, copy(x0); return_state = true)
x = get_solver_result(st)
pr("rosen_arc_default", (x, rf(R2, x), [c.at_iteration for c in Manopt.get_stopping_criterion(st).criteria]))
st = trust_regions(R2, rf, rg, rh, copy(x0); return_state = true)
x = get_solver_result(st)
pr("rosen_tr_default", (x, rf(R2, x), [c.at_iteration for c in Manopt.get_stopping_criterion(st).criteria]))
st = quasi_Newton(R2, rf, rg, copy(x0); stopping_criterion = StopWhenAny(StopAfterIteration(30), StopWhenGradientNormLess(1e-6)),
                  return_state = true)
x = get_solver_result(st)
pr("rosen_qn_maxiters30", (x, rf(R2, x), [c.at_iteration for c in Manopt.get_stopping_criterion(st).criteria]))
st = gradient_descent(R2, rf, rg, copy(x0); return_state = true)
x = get_solver_result(st)
pr("rosen_gd_default", (x, rf(R2, x)))
