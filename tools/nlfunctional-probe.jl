# f call counts of the FSAL low order explicit RK methods, and the implicit
# methods with the functional, Anderson and Newton (always_new, relax)
# nonlinear solvers on a small linear system
using OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve: NLFunctional, NLAnderson, NLNewton
g(u, p, t) = -u
prob1 = ODEProblem(g, [1.0], (0.0, 1.0))
for (name, alg) in (("Midpoint", Midpoint()), ("Heun", Heun()), ("Ralston", Ralston()))
    for (mode, kw) in (("adaptive", (;)), ("fixed", (dt = 0.1, adaptive = false)))
        sol = solve(prob1, alg; kw...)
        println(join([name, mode, sol.stats.nf, sol.stats.naccept, sol.stats.nreject, repr(sol.u[end][1])], "\t"))
    end
end
f(u, p, t) = [-2 * u[1] + u[2], u[1] - 3 * u[2] + sin(t)]
prob = ODEProblem(f, [1.0, 0.5], (0.0, 2.0))
nlss = (("fun", NLFunctional()), ("and", NLAnderson()), ("anew", NLNewton(always_new = true)),
        ("relax", NLNewton(relax = 0.3)), ("newton", NLNewton()))
for (name, A) in (("ImplicitEuler", ImplicitEuler), ("Trapezoid", Trapezoid), ("TRBDF2", TRBDF2),
                  ("KenCarp4", KenCarp4), ("SDIRK2", SDIRK2))
    for (k, nls) in nlss
        s = try
            sol = solve(prob, A(nlsolve = nls))
            st = sol.stats
            join([sol.retcode, st.nf, st.naccept, st.nreject, st.njacs, st.nsolve, st.nw,
                  repr(sol.u[end][1]), repr(sol.u[end][2])], "\t")
        catch e
            "ERROR\t" * sprint(showerror, e)[1:min(end, 200)]
        end
        println(name, "\t", k, "\t", s)
    end
end
# always_new with the in place form: the out of place form writes
# uprev + gamma z into integrator.u, which aliases uprev
f!(du, u, p, t) = (du[1] = -2 * u[1] + u[2]; du[2] = u[1] - 3 * u[2] + sin(t); nothing)
probip = ODEProblem(f!, [1.0, 0.5], (0.0, 2.0))
for (name, A) in (("ImplicitEuler", ImplicitEuler), ("Trapezoid", Trapezoid), ("TRBDF2", TRBDF2),
                  ("KenCarp4", KenCarp4), ("SDIRK2", SDIRK2))
    sol = solve(probip, A(nlsolve = NLNewton(always_new = true)))
    st = sol.stats
    println(join([name, "anew_iip", sol.retcode, st.nf, st.naccept, st.nreject, st.njacs, st.nsolve, st.nw,
                  repr(sol.u[end][1]), repr(sol.u[end][2])], "\t"))
end
