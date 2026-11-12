# Reference values for the SciPy ports (R/optimization-scipy*.R): runs
# SciPy 1.18.1 with the options OptimizationSciPy.jl passes and prints R
# assignments with full precision. Usage:
#   python3 tools/scipy-probe.py [section ...] > refs.R
# (a virtual environment with scipy==1.18.1 is enough).
import sys
import numpy as np
import scipy
import scipy.optimize as so

assert scipy.__version__ == "1.18.1", scipy.__version__


def rv(v):
    a = np.atleast_1d(np.asarray(v, dtype=float))
    return "c(" + ", ".join(repr(float(t)) for t in a) + ")"


def emit(name, res, fields=("x", "fun", "nit", "nfev", "njev", "nhev", "status")):
    parts = []
    for f in fields:
        if f in res and res[f] is not None:
            val = res[f]
            if isinstance(val, (bool, np.bool_)):
                parts.append(f"{f} = {'TRUE' if val else 'FALSE'}")
            elif f in ("nit", "nfev", "njev", "nhev", "status"):
                parts.append(f"{f} = {int(val)}L")
            else:
                parts.append(f"{f} = {rv(val)}")
    print(f"ref${name} <- list({', '.join(parts)})")


def rosen(x):
    a = 1.0 - x[0]
    b = x[1] - x[0] * x[0]
    return a * a + 100.0 * b * b


def rosen_grad(x):
    b = x[1] - x[0] * x[0]
    return np.array([-2.0 * (1.0 - x[0]) - 400.0 * x[0] * b, 200.0 * b])


def rosen_hess(x):
    return np.array([[2.0 - 400.0 * x[1] + 1200.0 * x[0] * x[0], -400.0 * x[0]],
                     [-400.0 * x[0], 200.0]])


def quartic(x):
    # x^4 - 3 x^3 + 2, minimum at 9/4
    return x * x * x * x - 3.0 * x * x * x + 2.0


def poly3(x):
    # a smooth polynomial in three variables with minimum at (1, 2, -1)
    a = x[0] - 1.0
    b = x[1] - 2.0
    c = x[2] + 1.0
    return a * a + 2.0 * b * b + 3.0 * c * c + a * a * b * b + 0.5 * a * c


def poly3_grad(x):
    a = x[0] - 1.0
    b = x[1] - 2.0
    c = x[2] + 1.0
    return np.array([2.0 * a + 2.0 * a * b * b + 0.5 * c, 4.0 * b + 2.0 * a * a * b, 6.0 * c + 0.5 * a])


def poly3_hess(x):
    a = x[0] - 1.0
    b = x[1] - 2.0
    return np.array([[2.0 + 2.0 * b * b, 4.0 * a * b, 0.5], [4.0 * a * b, 4.0 + 2.0 * a * a, 0.0],
                     [0.5, 0.0, 6.0]])


sections = set(sys.argv[1:]) or {"minimize", "scalar"}
x0 = np.array([-1.2, 1.0])

if "minimize" in sections:
    for m in ["Nelder-Mead", "Powell", "CG", "BFGS"]:
        jac = rosen_grad if m in ("CG", "BFGS") else None
        emit(m.replace("-", "") + "_rosen", so.minimize(rosen, x0, method=m, jac=jac, options={}))
    emit("NelderMead_poly3", so.minimize(poly3, [0.0, 0.0, 0.0], method="Nelder-Mead"))
    emit("Powell_poly3", so.minimize(poly3, [0.0, 0.0, 0.0], method="Powell"))
    emit("BFGS_poly3_fd", so.minimize(poly3, [0.0, 0.0, 0.0], method="BFGS"))
    emit("CG_poly3_fd", so.minimize(poly3, [0.0, 0.0, 0.0], method="CG"))
    emit("NewtonCG_rosen", so.minimize(rosen, x0, method="Newton-CG", jac=rosen_grad, hess=rosen_hess))
    b = so.Bounds([-2.0, 0.5], [0.8, 2.0])
    emit("NelderMead_rosen_bounds", so.minimize(rosen, x0, method="Nelder-Mead", bounds=b))
    emit("Powell_rosen_bounds", so.minimize(rosen, x0, method="Powell", bounds=b))
    emit("BFGS_rosen_maxiter", so.minimize(rosen, x0, method="BFGS", jac=rosen_grad, options={"maxiter": 5}))

if "scalar" in sections:
    emit("brent_quartic", so.minimize_scalar(quartic, method="brent"), ("x", "fun", "nit", "nfev"))
    emit("golden_quartic", so.minimize_scalar(quartic, method="golden"), ("x", "fun", "nit", "nfev"))
    emit("bounded_quartic", so.minimize_scalar(quartic, bounds=(-1.0, 4.0), method="bounded"),
         ("x", "fun", "nit", "nfev", "status"))
    r = so.brute(rosen, [(-2.0, 2.0), (-1.0, 3.0)], Ns=20, full_output=True, finish=so.fmin, workers=1)
    print(f"ref$brute_rosen <- list(x = {rv(r[0])}, fun = {rv(r[1])})")

if "trust" in sections:
    tro = {"gtol": 1e-10, "maxiter": 50000}
    emit("dogleg_rosen", so.minimize(rosen, x0, method="dogleg", jac=rosen_grad, hess=rosen_hess, options=tro))
    emit("trustncg_rosen", so.minimize(rosen, x0, method="trust-ncg", jac=rosen_grad, hess=rosen_hess, options=tro))
    emit("trustexact_rosen", so.minimize(rosen, x0, method="trust-exact", jac=rosen_grad, hess=rosen_hess,
                                         options=tro))
    emit("trustncg_rosen_bfgs", so.minimize(rosen, x0, method="trust-ncg", jac=rosen_grad, hess=so.BFGS(),
                                            options=tro))
    emit("trustexact_poly3", so.minimize(poly3, [3.0, -1.0, 0.5], method="trust-exact", jac=poly3_grad,
                                         hess=poly3_hess, options=tro))
    emit("dogleg_poly3", so.minimize(poly3, [3.0, -1.0, 0.5], method="dogleg", jac=poly3_grad,
                                     hess=poly3_hess, options=tro))
    emit("trustkrylov_rosen", so.minimize(rosen, x0, method="trust-krylov", jac=rosen_grad, hess=rosen_hess,
                                          options=tro))
    emit("trustkrylov_poly3", so.minimize(poly3, [3.0, -1.0, 0.5], method="trust-krylov", jac=poly3_grad,
                                          hess=poly3_hess, options=tro))

if "lbfgsb" in sections:
    emit("lbfgsb_rosen", so.minimize(rosen, x0, method="L-BFGS-B", jac=rosen_grad))
    emit("lbfgsb_rosen_fd", so.minimize(rosen, x0, method="L-BFGS-B"))
    emit("lbfgsb_rosen_m2", so.minimize(rosen, x0, method="L-BFGS-B", jac=rosen_grad, options={"maxcor": 2}))
    emit("lbfgsb_rosen_bounds", so.minimize(rosen, x0, method="L-BFGS-B", jac=rosen_grad,
                                            bounds=so.Bounds([-2.0, 1.5], [0.8, 3.0])))
    emit("lbfgsb_rosen_lower", so.minimize(rosen, [2.0, 3.0], method="L-BFGS-B", jac=rosen_grad,
                                           bounds=so.Bounds([1.2, -np.inf], [np.inf, np.inf])))
    emit("lbfgsb_poly3_bounds_fd", so.minimize(poly3, [0.0, 0.0, 0.0], method="L-BFGS-B",
                                               bounds=so.Bounds([-1.0, -1.0, -0.5], [0.5, 1.5, 1.0])))
    emit("lbfgsb_poly3_upper", so.minimize(poly3, [0.0, 0.0, 0.0], method="L-BFGS-B", jac=poly3_grad,
                                           bounds=so.Bounds([-np.inf, -np.inf, -np.inf], [np.inf, 1.0, np.inf]),
                                           options={"maxcor": 1}))
    emit("lbfgsb_rosen_maxiter", so.minimize(rosen, x0, method="L-BFGS-B", jac=rosen_grad, options={"maxiter": 7}))

if "cobyla" in sections:
    cf = ("x", "fun", "nfev", "status", "maxcv")
    emit("cobyla_rosen", so.minimize(rosen, x0, method="COBYLA"), cf)
    emit("cobyla_poly3", so.minimize(poly3, [0.0, 0.0, 0.0], method="COBYLA"), cf)
    emit("cobyla_poly3_bounds", so.minimize(poly3, [0.0, 0.0, 0.0], method="COBYLA",
                                            bounds=so.Bounds([-1.0, -1.0, -0.5], [0.5, 1.5, 1.0])), cf)
    disk = so.NonlinearConstraint(lambda x: np.array([x[0] * x[0] + x[1] * x[1]]), -np.inf, 1.5)
    emit("cobyla_rosen_disk", so.minimize(rosen, x0, method="COBYLA", constraints=[disk]), cf)
    curve = so.NonlinearConstraint(lambda x: np.array([x[0] + x[1] * x[1], x[0] * x[2]]), [1.0, -0.5], [1.0, 0.5])
    emit("cobyla_poly3_eq", so.minimize(poly3, [0.0, 0.0, 0.0], method="COBYLA", constraints=[curve]), cf)
    emit("cobyla_poly3_fixed", so.minimize(poly3, [0.0, 0.0, 0.0], method="COBYLA",
                                           bounds=so.Bounds([-1.0, 0.25, -0.5], [0.5, 0.25, 1.0])), cf)
    emit("cobyla_rosen_tol", so.minimize(rosen, x0, method="COBYLA", options={"tol": 1e-8, "maxiter": 3000}), cf)

if "zeros" in sections:
    def cub(x):
        return x * x * x - 2.0 * x - 5.0

    def cub_p(x):
        return 3.0 * x * x - 2.0

    def cub_pp(x):
        return 6.0 * x

    def remit(name, r):
        print(f"ref${name} <- list(root = {float(r.root)!r}, iterations = {int(r.iterations)}L, "
              f"function_calls = {int(r.function_calls)}L, converged = {'TRUE' if r.converged else 'FALSE'})")

    for m in ["brentq", "brenth", "bisect", "ridder", "toms748"]:
        remit(m + "_cub", so.root_scalar(cub, method=m, bracket=(2.0, 3.0)))
        remit(m + "_cos", so.root_scalar(lambda x: np.cos(x) - x, method=m, bracket=(0.0, 1.0)))
        remit(m + "_wide", so.root_scalar(cub, method=m, bracket=(-1.0, 4.0)))
    remit("newton_cub", so.root_scalar(cub, method="newton", x0=1.5, fprime=cub_p))
    remit("newton_fd_cub", so.root_scalar(cub, method="newton", x0=1.5))
    remit("secant_cub", so.root_scalar(cub, method="secant", x0=1.5))
    remit("secant2_cub", so.root_scalar(cub, method="secant", x0=1.5, x1=2.5))
    remit("halley_cub", so.root_scalar(cub, method="halley", x0=1.5, fprime=cub_p, fprime2=cub_pp))
    remit("brentq_limited", so.root_scalar(cub, method="brentq", bracket=(2.0, 3.0), xtol=1e-3, maxiter=3))

if "root" in sections:
    # scipy.optimize.root (MINPACK hybr and lm) and least_squares(method="lm")
    def sys2(x):
        d = x[0] - x[1]
        return np.array([x[0] + 0.5 * d * d * d - 1.0, -0.5 * d * d * d + x[1]])

    def sys2_jac(x):
        d = x[0] - x[1]
        return np.array([[1 + 1.5 * d * d, -1.5 * d * d], [-1.5 * d * d, 1 + 1.5 * d * d]])

    def sys3(x):
        return np.array([x[0] * x[0] + x[1] * x[1] - 4.0, x[0] * x[1] - 1.0, x[2] * x[2] * x[2] + x[0] - 2.0])

    def res3(x):
        return np.array([10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0], 0.5 * (x[0] * x[1] - 1.0)])

    for name, fun, jac, x0 in [("sys2", sys2, None, [0.0, 0.0]), ("sys2j", sys2, sys2_jac, [0.0, 0.0]),
                               ("sys3", sys3, None, [1.0, 1.0, 1.0]), ("sys3far", sys3, None, [3.0, -2.0, 5.0])]:
        for method in ("hybr", "lm"):
            emit(f"{name}_{method}", so.root(fun, x0, jac=jac, method=method),
                 ("x", "fun", "nfev", "njev", "status"))
    for name, fun, x0 in [("res3", res3, [-1.2, 1.0]), ("sys3", sys3, [3.0, -2.0, 5.0])]:
        r = so.least_squares(fun, x0, method="lm", jac="2-point")
        emit(f"lsq_{name}_lm", r, ("x", "cost", "fun", "nfev", "status", "optimality"))

if "lsq" in sections:
    # least_squares(method = "trf" / "dogbox") as OptimizationSciPy.jl calls
    # it (jac "2-point"), and root(method = "df-sane")
    def res3(x):
        return np.array([10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0], 0.5 * (x[0] * x[1] - 1.0)])

    def sys3(x):
        return np.array([x[0] * x[0] + x[1] * x[1] - 4.0, x[0] * x[1] - 1.0, x[2] * x[2] * x[2] + x[0] - 2.0])

    tt = np.linspace(0.0, 2.0, 9)
    yy = 1.0 + 0.5 * tt - 0.3 * tt * tt + np.array([0.1, -0.2, 0.05, 3.0, -0.1, 0.0, 0.2, -2.5, 0.1])

    def fit(x):
        return x[0] + x[1] * tt + x[2] * tt * tt * x[0] - yy

    fields = ("x", "cost", "fun", "nfev", "njev", "status", "optimality")
    for method in ("trf", "dogbox"):
        emit(f"{method}_res3", so.least_squares(res3, [-1.2, 1.0], method=method, jac="2-point"), fields)
        emit(f"{method}_res3_b", so.least_squares(res3, [-1.2, 1.0], method=method, jac="2-point",
                                                  bounds=([-2.0, -2.0], [0.8, 2.0])), fields)
        emit(f"{method}_sys3", so.least_squares(sys3, [3.0, -2.0, 5.0], method=method, jac="2-point"), fields)
        for loss in ("soft_l1", "huber", "cauchy", "arctan"):
            emit(f"{method}_fit_{loss}", so.least_squares(fit, [1.0, 0.0, 0.0], method=method, jac="2-point",
                                                          loss=loss), fields)
        emit(f"{method}_fit_b", so.least_squares(fit, [1.0, 0.7, 0.0], method=method, jac="2-point", loss="soft_l1",
                                                 bounds=([0.0, 0.6, -1.0], [2.0, 2.0, 1.0])), fields)
    for ls in ("cruz", "cheng"):
        emit(f"dfsane_sys3_{ls}", so.root(sys3, [3.0, -2.0, 5.0], method="df-sane", options={"line_search": ls}),
             ("x", "fun", "nfev", "nit"))

if "lp" in sections:
    # linprog and milp (HiGHS) as OptimizationSciPy.jl calls them, on
    # problems with a unique optimum
    lpf = ("x", "fun", "slack", "con", "status", "nit")

    def emit_lp(name, r):
        emit(name, r, lpf)
        d = {}
        for k in ("ineqlin", "eqlin", "lower", "upper"):
            if hasattr(r, k) and getattr(r, k) is not None and getattr(r, k).marginals is not None:
                d[k] = getattr(r, k).marginals
        print(f"ref${name}_marg <- list(" + ", ".join(f"{k} = c({', '.join(repr(float(v)) for v in np.atleast_1d(m))})"
                                                       for k, m in d.items()) + ")")

    c1 = [-1.0, -2.0, 0.5]
    A1 = [[1.0, 1.0, 1.0], [2.0, 0.5, 0.0], [0.0, 1.0, -1.0]]
    b1 = [4.0, 5.0, 1.5]
    Ae = [[1.0, -1.0, 2.0]]
    be = [0.5]
    bnd = [(0.0, 3.0), (-1.0, None), (None, 2.0)]
    for meth in ("highs", "highs-ds", "highs-ipm"):
        emit_lp(f"lp1_{meth}", so.linprog(c1, A_ub=A1, b_ub=b1, A_eq=Ae, b_eq=be, bounds=bnd, method=meth))
    emit_lp("lp_free", so.linprog([1.0, 1.0], A_ub=[[-1.0, -2.0], [-3.0, -1.0]], b_ub=[-4.0, -6.0]))
    emit_lp("lp_infeas", so.linprog([1.0, 1.0], A_ub=[[1.0, 1.0]], b_ub=[-1.0], bounds=[(0.0, None), (0.0, None)]))
    emit_lp("lp_unb", so.linprog([-1.0, 0.0], A_ub=[[1.0, -1.0]], b_ub=[1.0], bounds=[(0.0, None), (0.0, None)]))
    mf = ("x", "fun", "status", "mip_node_count", "mip_dual_bound", "mip_gap")
    cm = np.array([-5.0, -4.0, -3.0, -1.0])
    Am = np.array([[2.0, 3.0, 1.0, 1.0], [4.0, 1.0, 2.0, 0.5], [3.0, 4.0, 2.0, 1.0]])
    emit("milp1", so.milp(c=cm, integrality=[1, 1, 1, 0], bounds=so.Bounds([0.0] * 4, [10.0, 10.0, 10.0, 2.5]),
                          constraints=so.LinearConstraint(Am, [-np.inf] * 3, [5.0, 11.0, 8.0])), mf)
    emit("milp2", so.milp(c=np.array([-1.0, -1.0]), integrality=[1, 1], bounds=so.Bounds([0.0, 0.0], [np.inf, np.inf]),
                          constraints=so.LinearConstraint(np.array([[2.0, 2.0], [1.0, -1.0]]), [-np.inf, 0.0], [7.0, 0.0])), mf)
    emit("milp_infeas", so.milp(c=np.array([1.0]), integrality=[1], bounds=so.Bounds([0.2], [0.8])), mf)

if "nonlin" in sections:
    # root() with the quasi Newton methods of _nonlin.py
    def fq(x):
        return np.array([x[0] + 0.1 * x[0] * x[0] * x[0] - 0.2 * x[1] - 1.0,
                         x[1] + 0.1 * x[1] * x[1] * x[1] + 0.2 * x[0] - 0.5,
                         x[2] + 0.05 * x[0] * x[2] - 0.3])

    nf = ("x", "fun", "nit", "status")
    for m in ("broyden1", "broyden2", "anderson", "linearmixing", "diagbroyden", "excitingmixing"):
        emit(f"nl_{m}", so.root(fq, [0.0, 0.0, 0.0], method=m), nf)
    emit("nl_b1_wolfe", so.root(fq, [0.0, 0.0, 0.0], method="broyden1", options={"line_search": "wolfe"}), nf)
    emit("nl_b1_none", so.root(fq, [0.0, 0.0, 0.0], method="broyden1", options={"line_search": None}), nf)
    emit("nl_b1_svd", so.root(fq, [0.0, 0.0, 0.0], method="broyden1",
                              options={"jac_options": {"reduction_method": "svd", "max_rank": 2}}), nf)
    emit("nl_b2_simple", so.root(fq, [0.0, 0.0, 0.0], method="broyden2",
                                 options={"jac_options": {"reduction_method": "simple", "max_rank": 2}}), nf)
    emit("nl_and_tol", so.root(fq, [0.0, 0.0, 0.0], method="anderson", tol=1e-10,
                               options={"jac_options": {"M": 3, "alpha": 0.8}}), nf)
    emit("nl_lm_nit", so.root(fq, [0.0, 0.0, 0.0], method="linearmixing", options={"nit": 4}), nf)

if "lp2" in sections:
    # active lower and upper bounds, a free variable, the legacy simplex and
    # semicontinuous milp variables
    def emit_marg(name, r):
        emit(name, r, ("x", "fun", "status"))
        print(f"ref${name}_marg <- list(" + ", ".join(
            f"{k} = c({', '.join(repr(float(v)) for v in np.atleast_1d(getattr(r, k).marginals))})"
            for k in ("ineqlin", "eqlin", "lower", "upper")) + ")")

    c2 = [1.0, -1.0, 2.0, 0.5]
    A2 = [[1.0, 1.0, 1.0, 1.0], [1.0, -2.0, 0.0, 1.0]]
    b2 = [5.0, 2.0]
    Ae2 = [[0.0, 1.0, 1.0, -1.0]]
    be2 = [1.0]
    bnd2 = [(1.0, 4.0), (-2.0, 3.0), (None, None), (0.0, 2.0)]
    emit_marg("lp2", so.linprog(c2, A_ub=A2, b_ub=b2, A_eq=Ae2, b_eq=be2, bounds=bnd2))
    r = so.linprog([1.0, 1.0], A_ub=[[-1.0, -2.0], [-3.0, -1.0]], b_ub=[-4.0, -6.0], method="simplex")
    emit("lp_simplex", r, ("x", "fun", "status", "nit"))
    print(f"ref$lp_simplex_msg <- {r.message!r}")
    emit("milp_semi", so.milp(c=np.array([-1.0, -2.0]), integrality=[2, 3], bounds=so.Bounds([2.0, 1.5], [5.0, 4.0]),
                              constraints=so.LinearConstraint(np.array([[1.0, 1.0], [1.0, 3.0]]), [-np.inf, -np.inf],
                                                              [5.5, 10.0])), ("x", "fun", "status"))
