# Reference values for TrustRegion of NonlinearSolveFirstOrder (trust_region.jl:
# More's subproblem and radius update by default, the dogleg subproblem and
# the other radius updates on request), with analytic Jacobians. Prints R
# assignments. Usage:
#   julia --project=<env with NonlinearSolve> tools/trustregion-probe.jl > refs.R
using NonlinearSolve, LinearAlgebra

rv(v) = "c(" * join([repr(Float64(t)) for t in v], ", ") * ")"
function emit(name, sol)
    println("tr\$$name <- list(u = $(rv(sol.u)), resid = $(rv(sol.resid)), retcode = \"$(sol.retcode)\", ",
            "nsteps = $(sol.stats.nsteps)L, nf = $(sol.stats.nf)L)")
end

rosen(u, p) = [10 * (u[2] - u[1] * u[1]), 1 - u[1]]
rosenj(u, p) = [-20*u[1] 10.0; -1.0 0.0]
circ(u, p) = [u[1] * u[1] + u[2] * u[2] - 2, u[1] - u[2]]
circj(u, p) = [2*u[1] 2*u[2]; 1.0 -1.0]
lsq(u, p) = [u[1] - 1, u[1] + 1, u[2] - 3]
lsqj(u, p) = [1.0 0.0; 1.0 0.0; 0.0 1.0]
blsq(u, p) = [u[1] - 2 * u[2] - 1, u[1] - u[2]]
blsqj(u, p) = [1.0 -2.0; 1.0 -1.0]

pr = NonlinearProblem(NonlinearFunction(rosen; jac = rosenj), [-1.2, 1.0])
pc = NonlinearProblem(NonlinearFunction(circ; jac = circj), [1.0, 2.0])
pl = NonlinearLeastSquaresProblem(NonlinearFunction(lsq; jac = lsqj, resid_prototype = zeros(3)), [5.0, 5.0])
pb = NonlinearLeastSquaresProblem(NonlinearFunction(blsq; jac = blsqj, resid_prototype = zeros(2)), [0.0, 0.0];
                                  lb = [0.0, 0.0], ub = [Inf, Inf])

emit("rosen_default", solve(pr, TrustRegion()))
emit("circ_default", solve(pc, TrustRegion()))
emit("lsq_default", solve(pl, TrustRegion()))
emit("blsq_default", solve(pb, TrustRegion()))
emit("rosen_dogleg", solve(pr, TrustRegion(; subproblem = TrustRegionSubproblem.Dogleg)))
emit("circ_dogleg", solve(pc, TrustRegionDogleg()))
emit("rosen_robust", solve(pr, TrustRegionRobust()))
for (nm, s) in (("Simple", RadiusUpdateSchemes.Simple), ("NLsolve", RadiusUpdateSchemes.NLsolve),
                ("NocedalWright", RadiusUpdateSchemes.NocedalWright), ("Hei", RadiusUpdateSchemes.Hei),
                ("Yuan", RadiusUpdateSchemes.Yuan), ("More", RadiusUpdateSchemes.More),
                ("Bastin", RadiusUpdateSchemes.Bastin), ("Fan", RadiusUpdateSchemes.Fan))
    emit("rosen_$nm", solve(pr, TrustRegion(; radius_update_scheme = s)))
    emit("circ_$nm", solve(pc, TrustRegion(; radius_update_scheme = s)))
    emit("rosen_dogleg_$nm", solve(pr, TrustRegion(; subproblem = TrustRegionSubproblem.Dogleg,
                                                   radius_update_scheme = s)))
end
