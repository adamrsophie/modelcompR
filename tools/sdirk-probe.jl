# Reference values for tests/testthat/test-sdirk.R: the SDIRK, ESDIRK and
# IMEX methods of OrdinaryDiffEqSDIRK, out of place, on a nonlinear heat
# equation (plain and split into a stiff diffusion part f1 and a reaction
# part f2) and on the Robertson problem in mass matrix form, all with
# analytic Jacobians. Each line holds the label, the case, the number of
# saved points, nf, nf2, njacs, nw, nsolve, naccept, nreject, the final
# state and the dense output at t = 0.123. Run from the repository root:
#   julia --project=tools/julia tools/sdirk-probe.jl > sdirk.tsv
using OrdinaryDiffEqSDIRK, OrdinaryDiffEqCore, SciMLBase
const P = OrdinaryDiffEqSDIRK.Predictor
const N = 8
heat(u, p, t) = [p[1] * ((i == 1 ? zero(u[1]) : u[i-1]) - 2u[i] + (i == N ? zero(u[1]) : u[i+1])) * (N + 1)^2 - 0.5 * u[i]^3 + 0.1 * t for i in 1:N]
f1(u, p, t) = [p[1] * ((i == 1 ? zero(u[1]) : u[i-1]) - 2u[i] + (i == N ? zero(u[1]) : u[i+1])) * (N + 1)^2 - 0.5 * u[i] * u[i] * u[i] for i in 1:N]
f2(u, p, t) = [p[2] * u[i] * (1 - u[i]) + 0.1 * t for i in 1:N]
function hjac(u, p, t)
    J = zeros(N, N)
    for i in 1:N
        J[i, i] = -2 * p[1] * (N + 1)^2 - 1.5 * u[i]^2
        i > 1 && (J[i, i-1] = p[1] * (N + 1)^2)
        i < N && (J[i, i+1] = p[1] * (N + 1)^2)
    end
    J
end
rober(u, p, t) = [-0.04 * u[1] + 1.0e4 * u[2] * u[3], 0.04 * u[1] - 1.0e4 * u[2] * u[3] - 3.0e7 * u[2] * u[2], u[1] + u[2] + u[3] - 1]
rjac(u, p, t) = [-0.04 1.0e4*u[3] 1.0e4*u[2]; 0.04 -1.0e4*u[3]-6.0e7*u[2] -1.0e4*u[2]; 1.0 1.0 1.0]
u0 = [i * (N + 1 - i) / N^2 for i in 1:N]
prob = ODEProblem(ODEFunction(heat; jac = hjac), u0, (0.0, 0.3), [1.0])
sprob = SplitODEProblem(SplitFunction(f1, f2; jac = hjac), u0, (0.0, 0.3), [1.0, 5.0])
mprob = ODEProblem(ODEFunction(rober; jac = rjac, mass_matrix = [1.0 0 0; 0 1.0 0; 0 0 0]), [1.0, 0.0, 0.0], (0.0, 0.5))
cb = DiscreteCallback((u, t, integ) -> t == 0.1, integ -> (integ.u = integ.u * 0.5))
lim(u, integ, p, t) = nothing

function probe(label, case, pr, a)
    s = try
        case == "plain" ? solve(pr, a; reltol = 1e-4, abstol = 1e-6) :
        case == "loose" ? solve(pr, a; reltol = 1e-2, abstol = 1e-4) :
        case == "tstop" ? solve(pr, a; reltol = 1e-4, abstol = 1e-6, tstops = [0.1]) :
        case == "cb" ? solve(pr, a; reltol = 1e-4, abstol = 1e-6, tstops = [0.1], callback = cb) :
        case == "fixed" ? solve(pr, a; dt = 1 / 32, adaptive = false) :
        case == "fixedcb" ? solve(pr, a; dt = 1 / 32, adaptive = false, tstops = [0.1], callback = cb) :
        case == "mm" ? solve(pr, a; reltol = 1e-4, abstol = 1e-6, initializealg = SciMLBase.NoInit()) :
        case == "mmfixed" ? solve(pr, a; dt = 1 / 64, adaptive = false, initializealg = SciMLBase.NoInit()) :
        error("case")
    catch e
        println(label, "\t", case, "\tERROR\t", replace(sprint(showerror, e), r"\s+" => " ")[1:min(end, 200)])
        return
    end
    st = s.stats
    println(label, "\t", case, "\t", length(s.t), "\t", st.nf, "\t", st.nf2, "\t", st.njacs, "\t", st.nw, "\t",
        st.nsolve, "\t", st.naccept, "\t", st.nreject, "\t", join(repr.(s.u[end]), ","), "\t", join(repr.(s(0.123)), ","))
end

adaptive_pure = [("TRBDF2", TRBDF2()), ("SDIRK2", SDIRK2()), ("Kvaerno3", Kvaerno3()), ("Kvaerno4", Kvaerno4()),
    ("Kvaerno5", Kvaerno5()), ("KenCarp3", KenCarp3()), ("KenCarp4", KenCarp4()), ("KenCarp47", KenCarp47()),
    ("KenCarp5", KenCarp5()), ("KenCarp58", KenCarp58()), ("Cash4", Cash4()), ("Cash4_e2", Cash4(embedding = 2)),
    ("Hairer4", Hairer4()), ("Hairer42", Hairer42()), ("ESDIRK54I8L2SA", ESDIRK54I8L2SA()),
    ("ESDIRK436L2SA2", ESDIRK436L2SA2()), ("ESDIRK437L2SA", ESDIRK437L2SA()), ("ESDIRK547L2SA2", ESDIRK547L2SA2()),
    ("ESDIRK659L2SA", ESDIRK659L2SA()), ("ESDIRK325L2SA", ESDIRK325L2SA())]
fixed_pure = [("ImplicitMidpoint", ImplicitMidpoint()), ("SSPSDIRK2", SSPSDIRK2()), ("SFSDIRK4", SFSDIRK4()),
    ("SFSDIRK5", SFSDIRK5()), ("SFSDIRK6", SFSDIRK6()), ("SFSDIRK7", SFSDIRK7()), ("SFSDIRK8", SFSDIRK8())]
adaptive_split = [("KenCarp3", KenCarp3()), ("KenCarp4", KenCarp4()), ("KenCarp47", KenCarp47()),
    ("KenCarp5", KenCarp5()), ("KenCarp58", KenCarp58())]
fixed_split = [("ARS222", ARS222()), ("ARS232", ARS232()), ("ARS343", ARS343()), ("ARS443", ARS443()),
    ("IMEXSSP222", IMEXSSP222()), ("IMEXSSP2322", IMEXSSP2322()), ("IMEXSSP3332", IMEXSSP3332()),
    ("IMEXSSP3433", IMEXSSP3433()), ("BHR553", BHR553()), ("CFNLIRK3", CFNLIRK3())]

for (n, a) in adaptive_pure, c in ("plain", "loose", "tstop", "cb")
    probe(n, c, prob, a)
end
for (n, a) in fixed_pure, c in ("fixed", "fixedcb")
    probe(n, c, prob, a)
end
for (n, a) in [("TRBDF2", TRBDF2()), ("SDIRK2", SDIRK2()), ("Kvaerno4", Kvaerno4()), ("KenCarp4", KenCarp4()), ("Cash4", Cash4())]
    probe(n, "fixed", prob, a)
end
for (n, a) in adaptive_split, c in ("plain", "loose", "tstop", "cb")
    probe("split_" * n, c, sprob, a)
end
for (n, a) in fixed_split, c in ("fixed", "fixedcb")
    probe("split_" * n, c, sprob, a)
end
probe("split_KenCarp4", "fixed", sprob, KenCarp4())
probe("split_KenCarp4", "fixedcb", sprob, KenCarp4())
# split methods on a plain problem, plain methods on a split problem
for (n, a) in fixed_split
    probe(n, "fixed", prob, a)
end
probe("split_Kvaerno4", "plain", sprob, Kvaerno4())
probe("split_TRBDF2", "plain", sprob, TRBDF2())
probe("split_SFSDIRK4", "fixed", sprob, SFSDIRK4())
# options
probe("TRBDF2_nosmooth", "plain", prob, TRBDF2(smooth_est = false))
probe("KenCarp4_nosmooth", "plain", prob, KenCarp4(smooth_est = false))
probe("ESDIRK436L2SA2_smooth", "plain", prob, ESDIRK436L2SA2(smooth_est = true))
probe("TRBDF2_limiter", "plain", prob, TRBDF2(step_limiter! = lim))
probe("split_KenCarp4_limiter", "plain", sprob, KenCarp4(step_limiter! = lim))
probe("split_ARS343_limiter", "fixed", sprob, ARS343(step_limiter! = lim))
# stage predictors
preds = [("Trivial", P.Trivial), ("Linear", P.Linear), ("CopyPrev", P.CopyPrev), ("StageExtrap", P.StageExtrap),
    ("Tableau", P.Tableau), ("MaxOrder", P.MaxOrder), ("VariableOrder", P.VariableOrder), ("CutoffOrder", P.CutoffOrder)]
for (pn, pv) in preds
    probe("TRBDF2_" * pn, "plain", prob, TRBDF2(predictor = pv))
    probe("Kvaerno4_" * pn, "plain", prob, Kvaerno4(predictor = pv))
    probe("SDIRK2_" * pn, "plain", prob, SDIRK2(predictor = pv))
    probe("Cash4_" * pn, "plain", prob, Cash4(predictor = pv))
    probe("ESDIRK436L2SA2_" * pn, "plain", prob, ESDIRK436L2SA2(predictor = pv))
    probe("SSPSDIRK2_" * pn, "fixed", prob, SSPSDIRK2(predictor = pv))
    probe("SFSDIRK5_" * pn, "fixed", prob, SFSDIRK5(predictor = pv))
    probe("Kvaerno4_" * pn, "fixed", prob, Kvaerno4(predictor = pv))
    probe("split_KenCarp4_" * pn, "plain", sprob, KenCarp4(predictor = pv))
    probe("split_IMEXSSP3433_" * pn, "fixed", sprob, IMEXSSP3433(predictor = pv))
end
# mass matrix (index 1 Robertson problem)
for (n, a) in [("TRBDF2", TRBDF2()), ("Kvaerno3", Kvaerno3()), ("KenCarp4", KenCarp4()), ("ESDIRK436L2SA2", ESDIRK436L2SA2()),
    ("SDIRK2", SDIRK2()), ("Cash4", Cash4()), ("Hairer4", Hairer4())]
    probe("mm_" * n, "mm", mprob, a)
end
for (n, a) in [("ImplicitMidpoint", ImplicitMidpoint()), ("SFSDIRK4", SFSDIRK4()), ("ARS343", ARS343()), ("TRBDF2", TRBDF2()), ("IMEXSSP3433", IMEXSSP3433())]
    probe("mm_" * n, "mmfixed", mprob, a)
end
