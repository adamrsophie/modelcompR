# Reference values for tests/testthat/test-optimization-nlpmodels.R.
# Run in an environment with NLPModels (0.21):
#   julia --project=<env with NLPModels> tools/nlpmodels-probe.jl
using NLPModels, LinearAlgebra

show_meta(m) = println(repr((m.ifix, m.ilow, m.iupp, m.irng, m.ifree, m.iinf, m.jfix, m.jlow, m.jupp,
    m.jrng, m.jfree, m.jinf, m.nnzj, m.nnzh, m.nln, m.nlin, m.nnln, m.minimize)))

# metadata analysis
m1 = NLPModelMeta(4; lvar = [0.0, -Inf, 1.0, -2.0], uvar = [0.0, 2.0, Inf, 3.0], ncon = 3,
    lcon = [0.0, -Inf, 1.0], ucon = [0.0, 1.0, Inf], lin = [2])
show_meta(m1)
m2 = NLPModelMeta(m1; minimize = false, name = "copy")
show_meta(m2)
println(repr((has_bounds(m1), bound_constrained(m1), unconstrained(m1), linearly_constrained(m1),
    equality_constrained(m1), inequality_constrained(m1), has_equalities(m1), has_inequalities(m1))))
m3 = NLPModelMeta(2)
println(repr((has_bounds(m3), unconstrained(m3), m3.nnzj, m3.nnzh, m3.jac_available)))
# build_nlpmodel_meta of a problem with bounds and two constraints
m4 = NLPModelMeta(3; ncon = 2, nnzj = 6, nnzh = 6, x0 = [1.0, 2.0, 3.0], y0 = zeros(2),
    lvar = [-1.0, -1.0, -1.0], uvar = [1.0, 1.0, 1.0], lcon = [0.0, -Inf], ucon = [0.0, 1.0], minimize = false)
show_meta(m4)
ns = NLSMeta(5, 3; lin = [1, 4])
println(repr((ns.nnzj, ns.nnzh, ns.nln, ns.nnln, ns.nlin)))

# a model with the counter conventions of ADNLPModels
mutable struct Toy <: AbstractNLPModel{Float64, Vector{Float64}}
    meta::NLPModelMeta{Float64, Vector{Float64}}
    counters::Counters
end
NLPModels.obj(nlp::Toy, x::AbstractVector) = (increment!(nlp, :neval_obj); (1 - x[1])^2 + 10 * (x[2] - x[1]^2)^2)
function NLPModels.grad!(nlp::Toy, x::AbstractVector, g::AbstractVector)
    increment!(nlp, :neval_grad)
    g[1] = -2 * (1 - x[1]) - 40 * x[1] * (x[2] - x[1]^2)
    g[2] = 20 * (x[2] - x[1]^2)
    g
end
NLPModels.cons_nln!(nlp::Toy, x::AbstractVector, c::AbstractVector) = (increment!(nlp, :neval_cons_nln); c[1] = x[1]^2 + x[2]^2; c)
NLPModels.jac_nln_structure!(nlp::Toy, r::AbstractVector, c::AbstractVector) = (r .= [1, 1]; c .= [1, 2]; (r, c))
NLPModels.jac_nln_coord!(nlp::Toy, x::AbstractVector, v::AbstractVector) = (increment!(nlp, :neval_jac_nln); v .= [2x[1], 2x[2]]; v)
NLPModels.jprod_nln!(nlp::Toy, x::AbstractVector, v::AbstractVector, Jv::AbstractVector) = (increment!(nlp, :neval_jprod_nln); Jv[1] = 2x[1] * v[1] + 2x[2] * v[2]; Jv)
NLPModels.jtprod_nln!(nlp::Toy, x::AbstractVector, v::AbstractVector, Jtv::AbstractVector) = (increment!(nlp, :neval_jtprod_nln); Jtv .= [2x[1], 2x[2]] .* v[1]; Jtv)
toy = Toy(NLPModelMeta(2; ncon = 1, lcon = [-Inf], ucon = [1.0], x0 = [0.5, 0.5]), Counters())
x = [0.3, -0.4]
println(repr((obj(toy, x), grad(toy, x), cons(toy, x), jac_coord(toy, x), jprod(toy, x, [1.0, 2.0]),
    jtprod(toy, x, [3.0]))))
c = toy.counters
println(repr((c.neval_obj, c.neval_grad, c.neval_cons, c.neval_cons_nln, c.neval_jac, c.neval_jac_nln,
    c.neval_jprod, c.neval_jprod_nln, c.neval_jtprod, c.neval_jtprod_nln, sum_counters(toy))))

# a least squares model: counters of obj and grad (the residual of
# Rosenbrock in least squares form)
mutable struct ToyNLS <: AbstractNLSModel{Float64, Vector{Float64}}
    meta::NLPModelMeta{Float64, Vector{Float64}}
    nls_meta::NLSMeta{Float64, Vector{Float64}}
    counters::NLSCounters
end
NLPModels.residual!(nls::ToyNLS, x::AbstractVector, F::AbstractVector) = (increment!(nls, :neval_residual); F .= [1 - x[1], 10 * (x[2] - x[1]^2)]; F)
J(x) = [-1.0 0.0; -20x[1] 10.0]
NLPModels.jtprod_residual!(nls::ToyNLS, x::AbstractVector, v::AbstractVector, Jtv::AbstractVector) = (increment!(nls, :neval_jtprod_residual); Jtv .= J(x)' * v; Jtv)
nls = ToyNLS(NLPModelMeta(2; x0 = [-1.2, 1.0]), NLSMeta(2, 2), NLSCounters())
println(repr((obj(nls, x), grad(nls, x))))
c = nls.counters
println(repr((c.neval_obj, c.neval_grad, c.neval_residual, c.neval_jtprod_residual, sum_counters(nls))))
