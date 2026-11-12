# Reference values for the native bounded solvers of NonlinearSolveFirstOrder
# (BoundedGaussNewton, BoundedLevenbergMarquardt, TrustRegionReflective, Dogbox,
# BoundedTrustRegion), with analytic Jacobians.
using NonlinearSolve, LinearAlgebra

println("LineSearch isa Module: ", NonlinearSolve.LineSearch isa Module)

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
# a square system whose root sits on the boundary region interior
circ(u, p) = [u[1]*u[1] + u[2]*u[2] - 4, u[1]*u[2] - 1]
circj(u, p) = [2*u[1] 2*u[2]; u[2] u[1]]

cases = [
  ("rosen_in", NonlinearProblem, rosen, rosenj, [-1.2, 1.0], [-2.0, -2.0], [2.0, 2.0]),
  ("rosen_bd", NonlinearProblem, rosen, rosenj, [-1.2, 1.0], [-2.0, -2.0], [0.5, 2.0]),
  ("powell", NonlinearProblem, powell, powellj, [3.0, -1.0, 0.0, 1.0], [-5.0, -5.0, -5.0, -5.0], [5.0, 5.0, 5.0, 5.0]),
  ("powell_bd", NonlinearProblem, powell, powellj, [3.0, -1.0, 0.0, 1.0], [0.5, -5.0, -5.0, -5.0], [5.0, 5.0, 5.0, 5.0]),
  ("mm_ls", NonlinearLeastSquaresProblem, mm, mmj, [1.0, 1.0], [0.0, 0.0], [2.5, 10.0]),
  ("mm_ls_free", NonlinearLeastSquaresProblem, mm, mmj, [1.0, 1.0], [0.0, 0.0], [10.0, 10.0]),
  ("circ", NonlinearProblem, circ, circj, [3.0, 0.5], [0.0, 0.0], [3.0, 3.0]),
]
algs = [
  ("BGN_ls", () -> BoundedGaussNewton()),
  ("BGN_tr", () -> BoundedGaussNewton(globalization = :trustregion)),
  ("BGN_trls", () -> BoundedGaussNewton(globalization = :trustregion_linesearch)),
  ("BLM", () -> BoundedLevenbergMarquardt()),
  ("TRR", () -> TrustRegionReflective()),
  ("Dogbox", () -> Dogbox()),
  ("BTR", () -> BoundedTrustRegion()),
]
for (cname, P, f, j, u0, lb, ub) in cases, (aname, mk) in algs, maxiters in (1000, 3)
  prob = P(NonlinearFunction(f; jac = j, resid_prototype = P === NonlinearLeastSquaresProblem ? zeros(5) : nothing), u0; lb = lb, ub = ub)
  sol = try
    solve(prob, mk(); maxiters = maxiters)
  catch e
    println("CASE ", cname, " ", aname, " ", maxiters, " ERROR ", sprint(showerror, e)[1:min(end, 200)])
    continue
  end
  s = sol.stats
  println("CASE ", cname, " ", aname, " ", maxiters, " ", sol.retcode, " ", repr(sol.u), " ", s.nf, " ", s.njacs, " ", s.nfactors, " ", s.nsolve, " ", s.nsteps)
end

# the leaves of the existing R port, to document whether they are iterate exact
for (aname, alg) in (("NewtonRaphson", NewtonRaphson()), ("TrustRegion", TrustRegion()), ("LevenbergMarquardt", LevenbergMarquardt()))
  prob = NonlinearProblem(NonlinearFunction(powell; jac = powellj), [3.0, -1.0, 0.0, 1.0])
  sol = solve(prob, alg)
  s = sol.stats
  println("LEAF ", aname, " ", sol.retcode, " ", repr(sol.u), " ", s.nf, " ", s.njacs, " ", s.nfactors, " ", s.nsolve, " ", s.nsteps)
end
