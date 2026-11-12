# Reference values for tests/testthat/test-amf.R: AMF of OrdinaryDiffEqAMF
# around the Rosenbrock methods on a 4 dimensional problem
#   du/dt = (A1 + A2) u - u^3 + (0.1 t, 0, 0, 0)
# whose Jacobian splits into the constant upper triangular A1 and the lower
# triangular A2 - 3 diag(u^2): W = -(I - gamma J1)(I - gamma J2) / gamma
# (split), the exact W = -(I - gamma J) / gamma (jac only) and the same
# product from prebuilt factors (amf_factors), with an event that halves the
# state at t = 1/2; then a three way split of
#   du/dt = (S + T) u - u^3 + (0.1 t, 0, 0, 0)
# into the symmetric S (Bunch-Kaufman factor), the tridiagonal T (LU of a
# tridiagonal factor) and the diagonal -3 diag(u^2). The time derivative is
# a forward difference. Julia's AMF needs an in place f! (the out of place
# Rosenbrock caches do not accept an operator W).
#   julia --project=tools/julia tools/amf-probe.jl
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqAMF, SciMLOperators, LinearAlgebra
using OrdinaryDiffEqDifferentiation: AutoFiniteDiff

A1 = [-3.0 1.0 0.5 0.0; 0.0 -2.0 0.25 1.0; 0.0 0.0 -4.0 0.5; 0.0 0.0 0.0 -1.0]
A2 = [0.0 0.0 0.0 0.0; 0.5 0.0 0.0 0.0; 0.25 1.0 0.0 0.0; 0.0 0.5 0.25 0.0]
function f!(du, u, p, t)
    mul!(du, A1 + A2, u)
    @. du -= u * u * u
    du[1] += 0.1 * t
    return nothing
end
upd2!(J, u, p, t) = (J .= A2; for i in 1:4; J[i, i] -= 3 * u[i] * u[i]; end; nothing)
u0 = [1.0, 0.5, -0.25, 0.75]
fd = AutoFiniteDiff()

function show(name, case, s)
    st = s.stats
    println(name, "\t", case, "\t", length(s.t), "\t", st.nf, "\t", st.njacs, "\t", st.nw, "\t", st.nsolve, "\t",
        st.naccept, "\t", st.nreject, "\t", join(repr.(s.u[end]), ","))
end

function fun(form)
    J1 = MatrixOperator(copy(A1))
    J2 = MatrixOperator(copy(A2); update_func! = upd2!)
    Jop = cache_operator(J1 + J2, zeros(4))
    if form == "split"
        build_amf_function(f!; jac = Jop, split = (J1, J2))
    elseif form == "jac"
        build_amf_function(f!; jac = Jop)
    else
        g() = ScalarOperator(1.0; update_func = (old, u, p, t; gamma = 1.0) -> gamma, accepted_kwargs = Val((:gamma,)))
        build_amf_function(f!; jac = Jop, amf_factors = (IdentityOperator(4) - g() * J1, IdentityOperator(4) - g() * J2))
    end
end

for form in ("split", "jac", "factors")
    prob = ODEProblem(fun(form), u0, (0.0, 2.0))
    for (nm, alg) in (("Rosenbrock23", Rosenbrock23), ("Rodas4", Rodas4), ("Rodas5P", Rodas5P))
        show(nm, form * "_adaptive", solve(prob, AMF(alg; autodiff = fd); reltol = 1e-6, abstol = 1e-8))
        show(nm, form * "_fixed", solve(prob, AMF(alg; autodiff = fd); dt = 1 / 32, adaptive = false))
    end
end

cb = DiscreteCallback((u, t, integ) -> t == 0.5, integ -> (integ.u .*= 0.5))
prob = ODEProblem(fun("split"), u0, (0.0, 2.0))
for (nm, alg) in (("Rosenbrock23", Rosenbrock23), ("Rodas4", Rodas4))
    show(nm, "cb_adaptive", solve(prob, AMF(alg; autodiff = fd); reltol = 1e-6, abstol = 1e-8, callback = cb,
        tstops = [0.5]))
    show(nm, "cb_fixed", solve(prob, AMF(alg; autodiff = fd); dt = 1 / 32, adaptive = false, callback = cb,
        tstops = [0.5]))
end

S = [-4.0 1.0 0.5 0.25; 1.0 -3.0 0.5 0.0; 0.5 0.5 -5.0 1.0; 0.25 0.0 1.0 -2.0]
T = [-1.0 0.5 0.0 0.0; 0.25 -1.0 0.5 0.0; 0.0 0.25 -1.0 0.5; 0.0 0.0 0.25 -1.0]
function g!(du, u, p, t)
    mul!(du, S + T, u)
    @. du -= u * u * u
    du[1] += 0.1 * t
    return nothing
end
updD!(J, u, p, t) = (fill!(J, 0.0); for i in 1:4; J[i, i] = -3 * u[i] * u[i]; end; nothing)
JS = MatrixOperator(copy(S))
JT = MatrixOperator(copy(T))
JD = MatrixOperator(zeros(4, 4); update_func! = updD!)
prob3 = ODEProblem(build_amf_function(g!; jac = cache_operator(JS + JT + JD, zeros(4)), split = (JS, JT, JD)),
    u0, (0.0, 2.0))
for (nm, alg) in (("Rosenbrock23", Rosenbrock23), ("Rodas5P", Rodas5P))
    show(nm, "sym_adaptive", solve(prob3, AMF(alg; autodiff = fd); reltol = 1e-6, abstol = 1e-8))
    show(nm, "sym_fixed", solve(prob3, AMF(alg; autodiff = fd); dt = 1 / 32, adaptive = false))
end
