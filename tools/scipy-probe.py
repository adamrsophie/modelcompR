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

if "krylov" in sections:
    # root(method = "krylov") (KrylovJacobian with lgmres) and least_squares
    # with tr_solver = "lsmr"
    def fq(x):
        return np.array([x[0] + 0.1 * x[0] * x[0] * x[0] - 0.2 * x[1] - 1.0,
                         x[1] + 0.1 * x[1] * x[1] * x[1] + 0.2 * x[0] - 0.5,
                         x[2] + 0.05 * x[0] * x[2] - 0.3])

    def sys2(x):
        d = x[0] - x[1]
        return np.array([x[0] + 0.5 * d * d * d - 1.0, -0.5 * d * d * d + x[1]])

    def res3(x):
        return np.array([10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0], 0.5 * (x[0] * x[1] - 1.0)])

    def sys3(x):
        return np.array([x[0] * x[0] + x[1] * x[1] - 4.0, x[0] * x[1] - 1.0, x[2] * x[2] * x[2] + x[0] - 2.0])

    tt = np.linspace(0.0, 2.0, 9)
    yy = 1.0 + 0.5 * tt - 0.3 * tt * tt + np.array([0.1, -0.2, 0.05, 3.0, -0.1, 0.0, 0.2, -2.5, 0.1])

    def fit(x):
        return x[0] + x[1] * tt + x[2] * tt * tt * x[0] - yy

    nf = ("x", "fun", "nit", "status")
    emit("kr_fq", so.root(fq, [0.0, 0.0, 0.0], method="krylov"), nf)
    emit("kr_fq_none", so.root(fq, [2.0, -1.0, 3.0], method="krylov", options={"line_search": None}), nf)
    emit("kr_sys2", so.root(sys2, [0.3, 0.1], method="krylov"), nf)
    emit("kr_fq_k2", so.root(fq, [0.0, 0.0, 0.0], method="krylov",
                             options={"jac_options": {"outer_k": 2, "inner_inner_m": 2}}), nf)
    fields = ("x", "cost", "nfev", "njev", "status", "optimality")
    for method in ("trf", "dogbox"):
        emit(f"lsmr_{method}_res3", so.least_squares(res3, [-1.2, 1.0], method=method, jac="2-point",
                                                     tr_solver="lsmr"), fields)
        emit(f"lsmr_{method}_res3_b", so.least_squares(res3, [-1.2, 1.0], method=method, jac="2-point",
                                                       tr_solver="lsmr", bounds=([-2.0, -2.0], [0.8, 2.0])), fields)
        emit(f"lsmr_{method}_sys3", so.least_squares(sys3, [3.0, -2.0, 5.0], method=method, jac="2-point",
                                                     tr_solver="lsmr"), fields)
        emit(f"lsmr_{method}_fit_b", so.least_squares(fit, [1.0, 0.7, 0.0], method=method, jac="2-point",
                                                      loss="soft_l1", tr_solver="lsmr",
                                                      bounds=([0.0, 0.6, -1.0], [2.0, 2.0, 1.0])), fields)
    emit("lsmr_trf_res3_noreg", so.least_squares(res3, [-1.2, 1.0], jac="2-point", tr_solver="lsmr",
                                                 tr_options={"regularize": False}), fields)
    emit("lsmr_trf_fit_damp", so.least_squares(fit, [1.0, 0.0, 0.0], jac="2-point", tr_solver="lsmr",
                                               tr_options={"damp": 0.1}), fields)

if "trlib" in sections:
    # trust-krylov (trlib) beyond the two cases of the trust section
    def chain(x):
        f = 0.0
        for i in range(len(x) - 1):
            a = 1.0 - x[i]
            b = x[i + 1] - x[i] * x[i]
            f += a * a + 100.0 * b * b
        return f

    def chain_grad(x):
        g = np.zeros(len(x))
        for i in range(len(x) - 1):
            b = x[i + 1] - x[i] * x[i]
            g[i] += -2.0 * (1.0 - x[i]) - 400.0 * x[i] * b
            g[i + 1] += 200.0 * b
        return g

    def chain_hess(x):
        n = len(x)
        H = np.zeros((n, n))
        for i in range(n - 1):
            H[i, i] += 2.0 - 400.0 * x[i + 1] + 1200.0 * x[i] * x[i]
            H[i, i + 1] += -400.0 * x[i]
            H[i + 1, i] += -400.0 * x[i]
            H[i + 1, i + 1] += 200.0
        return H

    def well(x):
        # a double well in each coordinate with a weak chain coupling
        f = 0.1 * x[0]
        for i in range(len(x)):
            u = x[i] * x[i] - 1.0
            f += u * u
        for i in range(len(x) - 1):
            f += 0.5 * x[i] * x[i + 1]
        return f

    def well_grad(x):
        g = 4.0 * x * (x * x - 1.0)
        g[0] += 0.1
        g[:-1] += 0.5 * x[1:]
        g[1:] += 0.5 * x[:-1]
        return g

    def well_hess(x):
        n = len(x)
        H = np.diag(12.0 * x * x - 4.0)
        for i in range(n - 1):
            H[i, i + 1] = 0.5
            H[i + 1, i] = 0.5
        return H

    tro = {"gtol": 1e-10, "maxiter": 50000}
    trx = {"gtol": 1e-10, "maxiter": 50000, "inexact": False}
    c0 = np.array([-1.2, 1.0, -1.2, 1.0, -1.2, 1.0])
    w0 = np.array([0.1, -0.1, 0.1, -0.1, 0.1])
    emit("tk_rosen_exact", so.minimize(rosen, x0, method="trust-krylov", jac=rosen_grad, hess=rosen_hess,
                                       options=trx))
    emit("tk_poly3_exact", so.minimize(poly3, [3.0, -1.0, 0.5], method="trust-krylov", jac=poly3_grad,
                                       hess=poly3_hess, options=trx))
    emit("tk_chain", so.minimize(chain, c0, method="trust-krylov", jac=chain_grad, hess=chain_hess, options=tro))
    emit("tk_chain_exact", so.minimize(chain, c0, method="trust-krylov", jac=chain_grad, hess=chain_hess,
                                       options=trx))
    emit("tk_well", so.minimize(well, w0, method="trust-krylov", jac=well_grad, hess=well_hess, options=tro))
    emit("tk_well_exact", so.minimize(well, w0, method="trust-krylov", jac=well_grad, hess=well_hess,
                                      options=trx))
    emit("tk_rosen_default", so.minimize(rosen, x0, method="trust-krylov", jac=rosen_grad, hess=rosen_hess))

if "slsqp" in sections:
    # SLSQP: one case per branch of slsqp.c (dgelsy, lsi with nnls, lsei with
    # equalities, the augmented problem, singular and dependent constraints)
    import warnings
    warnings.simplefilter("ignore")
    lo3 = [0.0, -np.inf, -0.5]
    hi3 = [0.8, 1.5, np.inf]

    def ceq_lin(x):
        return np.array([x[0] + x[1] + x[2] - 1.0])

    def ceq_lin_jac(x):
        return np.array([[1.0, 1.0, 1.0]])

    def ceq_sph(x):
        return np.array([x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 2.0])

    def cin_diff(x):
        return np.array([x[0] - x[2] - 0.2])

    def c_two(x):
        return np.array([x[0] * x[0] + x[1] * x[1] - 2.0, x[0] - x[1]])

    def c_two_jac(x):
        return np.array([[2.0 * x[0], 2.0 * x[1]], [1.0, -1.0]])

    def q2(x):
        a = x[0] - 2.0
        b = x[1] - 1.0
        return a * a + b * b

    def q2_grad(x):
        return np.array([2.0 * (x[0] - 2.0), 2.0 * (x[1] - 1.0)])

    sm = "SLSQP"
    emit("sq_rosen", so.minimize(rosen, x0, method=sm, jac=rosen_grad))
    emit("sq_rosen_fd", so.minimize(rosen, x0, method=sm))
    emit("sq_poly3_fd", so.minimize(poly3, [3.0, -1.0, 0.5], method=sm))
    emit("sq_rosen_iter", so.minimize(rosen, x0, method=sm, jac=rosen_grad, options={"maxiter": 5}))
    emit("sq_bounds", so.minimize(poly3, [0.5, 1.0, 0.0], method=sm, jac=poly3_grad, bounds=list(zip(lo3, hi3))))
    emit("sq_eq", so.minimize(poly3, [3.0, -1.0, 0.5], method=sm, jac=poly3_grad,
                              constraints={"type": "eq", "fun": ceq_lin, "jac": ceq_lin_jac}))
    emit("sq_mixed", so.minimize(poly3, [1.0, 0.5, 0.0], method=sm, jac=poly3_grad,
                                 bounds=[(-np.inf, np.inf), (-1.0, 0.9), (-np.inf, np.inf)],
                                 constraints=[{"type": "eq", "fun": ceq_sph}, {"type": "ineq", "fun": cin_diff}]))
    emit("sq_dep", so.minimize(poly3, [3.0, -1.0, 0.5], method=sm, jac=poly3_grad,
                               constraints=[{"type": "eq", "fun": ceq_lin, "jac": ceq_lin_jac},
                                            {"type": "eq", "fun": lambda x: 2.0 * ceq_lin(x),
                                             "jac": lambda x: 2.0 * ceq_lin_jac(x)}]))
    emit("sq_square", so.minimize(q2, [0.0, 0.0], method=sm, jac=q2_grad,
                                  constraints={"type": "eq", "fun": c_two, "jac": c_two_jac}))
    # the linearization at the start is inconsistent with the bounds
    emit("sq_incons", so.minimize(lambda x: (x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 1.0) * (x[1] - 1.0), [0.5, 0.5],
                                  method=sm, jac=lambda x: np.array([2.0 * (x[0] - 0.5), 2.0 * (x[1] - 1.0)]),
                                  bounds=[(0.0, 1.0), (0.0, 2.5)],
                                  constraints={"type": "ineq", "fun": lambda x: x[0] * x[0] + x[1] * x[1] - 4.0,
                                               "jac": lambda x: np.array([[2.0 * x[0], 2.0 * x[1]]])}))
    nlc = so.NonlinearConstraint(lambda x: np.array([x[0] * x[1] + x[2], x[0] - x[1]]), [-np.inf, -0.25],
                                 [1.0, -0.25])
    emit("sq_new_fixed", so.minimize(poly3, [0.5, 1.0, -0.5], method=sm, jac=poly3_grad,
                                     bounds=[(0.0, 2.0), (-np.inf, np.inf), (-0.5, -0.5)], constraints=nlc))
    emit("sq_new_two", so.minimize(poly3, [0.5, 1.0, 0.0], method=sm, jac=poly3_grad,
                                   constraints=so.NonlinearConstraint(lambda x: np.array([x[0] * x[1] + x[2]]),
                                                                      -0.5, 0.75)))
    # every variable fixed by its bounds, at a point that violates nlc
    r = so.minimize(poly3, [0.0, 0.0, 0.0], method=sm, bounds=[(1.0, 1.0), (2.0, 2.0), (0.5, 0.5)], constraints=nlc)
    emit("sq_allfixed", r, ("x", "fun", "nit", "nfev", "njev", "status", "success"))
    print("# sq_allfixed message:", r.message)

if "tnc" in sections:
    # TNC (tnc.c): free and bounded problems, the finite difference gradient,
    # fixed variables (always removed for TNC), the options and the limits
    import warnings
    warnings.simplefilter("ignore")
    tm = "TNC"

    def chain(x):
        f = 0.0
        for i in range(len(x) - 1):
            a = 1.0 - x[i]
            b = x[i + 1] - x[i] * x[i]
            f += a * a + 100.0 * b * b
        return f

    def chain_grad(x):
        g = np.zeros(len(x))
        for i in range(len(x) - 1):
            b = x[i + 1] - x[i] * x[i]
            g[i] += -2.0 * (1.0 - x[i]) - 400.0 * x[i] * b
            g[i + 1] += 200.0 * b
        return g

    c0 = [-1.2, 1.0, -1.2, 1.0, -1.2, 1.0]
    emit("tn_rosen", so.minimize(rosen, x0, method=tm, jac=rosen_grad))
    emit("tn_rosen_fd", so.minimize(rosen, x0, method=tm))
    emit("tn_poly3_fd", so.minimize(poly3, [3.0, -1.0, 0.5], method=tm))
    emit("tn_bounds", so.minimize(poly3, [0.5, 1.0, 0.0], method=tm, jac=poly3_grad,
                                  bounds=[(0.0, 0.8), (-np.inf, 1.5), (-0.5, np.inf)]))
    emit("tn_bounds_fd", so.minimize(poly3, [0.5, 1.0, 0.0], method=tm,
                                     bounds=[(0.0, 0.8), (-np.inf, 1.5), (-0.5, np.inf)]))
    emit("tn_fixed", so.minimize(poly3, [0.5, 1.0, 0.0], method=tm, jac=poly3_grad,
                                 bounds=[(0.0, 2.0), (1.25, 1.25), (-np.inf, np.inf)]))
    emit("tn_chain", so.minimize(chain, c0, method=tm, jac=chain_grad,
                                 bounds=[(-1.0, 0.9)] * 3 + [(-np.inf, np.inf)] * 3))
    emit("tn_chain_free", so.minimize(chain, c0, method=tm, jac=chain_grad, options={"maxfun": 2000}))
    # nine variables: two blocks of the vectorized dot products in msolve
    c9 = [-1.2 if i % 2 == 0 else 1.0 for i in range(9)]
    emit("tn_chain9", so.minimize(chain, c9, method=tm, jac=chain_grad,
                                  bounds=[(-2.0, 1.5)] * 2 + [(0.5, 1.5)] + [(-2.0, 1.5)] * 6, options={"maxfun": 3000}))
    emit("tn_maxfun", so.minimize(rosen, x0, method=tm, jac=rosen_grad, options={"maxfun": 20}))
    emit("tn_scale", so.minimize(rosen, x0, method=tm, jac=rosen_grad,
                                 options={"scale": np.array([2.0, 0.5]), "offset": np.array([0.1, -0.2])}))
    emit("tn_steepest", so.minimize(rosen, x0, method=tm, jac=rosen_grad, options={"maxCGit": 0, "maxfun": 300}))
    emit("tn_opts", so.minimize(rosen, x0, method=tm, jac=rosen_grad,
                                options={"eta": 0.5, "stepmx": 1.0, "accuracy": 1e-10, "ftol": 1e-12, "xtol": 1e-12,
                                         "gtol": 1e-12, "rescale": 0.0, "minfev": -1.0}))
    emit("tn_out", so.minimize(poly3, [3.0, 5.0, -4.0], method=tm, jac=poly3_grad,
                               bounds=[(0.0, 2.0), (0.0, 2.5), (-2.0, 0.0)]))

if "direct" in sections:
    # DIRECT (DIRect.c, DIRsubrout.c): the options of OptimizationSciPy.jl
    # (eps 1e-4, maxiter 1000, locally biased, vol_tol 1e-16, len_tol 1e-6),
    # the original algorithm, the limits and the stopping rules
    jl = {"eps": 1e-4, "maxiter": 1000, "locally_biased": True, "vol_tol": 1e-16, "len_tol": 1e-6}

    def drun(name, fun, bounds, **kw):
        opts = dict(jl)
        opts.update(kw)
        r = so.direct(fun, bounds, **opts)
        emit(name, r)
        print(f"ref${name}$message <- {r.message!r}")

    b2 = [(-2.0, 2.0), (-1.0, 3.0)]
    b3 = [(-1.0, 2.5), (0.0, 3.0), (-2.0, 1.0)]
    drun("dr_rosen", rosen, b2)
    drun("dr_poly3", poly3, b3)
    drun("dr_orig", rosen, b2, locally_biased=False)
    drun("dr_orig3", poly3, b3, locally_biased=False, maxfun=600)
    drun("dr_maxfun", poly3, b3, maxfun=200)
    drun("dr_maxiter", rosen, b2, maxiter=12)
    drun("dr_vol", poly3, b3, vol_tol=1e-8)
    drun("dr_fmin", rosen, b2, f_min=0.0, f_min_rtol=1e-3)
    drun("dr_eps", poly3, b3, eps=-1e-3, len_tol=1e-4)
    drun("dr_quartic", lambda x: quartic(x[0]), [(-1.0, 4.0)])

    # The same runs with the log of the C code switched on (disp, which
    # direct() always passes as False). DIRsubrout.c guards its `goto L60`
    # for a rectangle off the convex hull (helpgreater > helplower) with the
    # `if (logfile)` of a commented out fprintf, so only with a log are such
    # rectangles dropped as DIRECT prescribes. The R port drops them: these
    # are its references. The log goes to /dev/null.
    import os
    import ctypes
    import scipy.optimize._direct_py as dpy
    plain = dpy._direct
    libc = ctypes.CDLL(None)

    def logged(*a):
        a = list(a)
        a[4] = True
        sys.stdout.flush()
        saved = os.dup(1)
        null = os.open(os.devnull, os.O_WRONLY)
        os.dup2(null, 1)
        try:
            return plain(*a)
        finally:
            libc.fflush(None)
            os.dup2(saved, 1)
            os.close(null)
            os.close(saved)

    dpy._direct = logged
    drun("drf_rosen", rosen, b2)
    drun("drf_poly3", poly3, b3)
    drun("drf_orig", rosen, b2, locally_biased=False)
    drun("drf_orig3", poly3, b3, locally_biased=False, maxfun=600)
    drun("drf_maxfun", poly3, b3, maxfun=200)
    drun("drf_maxiter", rosen, b2, maxiter=12)
    drun("drf_vol", poly3, b3, vol_tol=1e-8)
    drun("drf_fmin", rosen, b2, f_min=0.0, f_min_rtol=1e-3)
    drun("drf_eps", poly3, b3, eps=-1e-3, len_tol=1e-4)
    drun("drf_quartic", lambda x: quartic(x[0]), [(-1.0, 4.0)])
    dpy._direct = plain

if "shgo" in sections:
    # shgo with simplicial sampling and the options of OptimizationSciPy.jl
    # (n 100, iters 1); the local minimizer is SLSQP with forward differences
    # unless minimizer_kwargs names another method. The double wells
    # (x0^2 - 1)^2 + (x1^2 - 1)^2 (plus a tilt) have four local minima.
    import multiprocessing

    def wells(x):
        a = x[0] * x[0] - 1.0
        b = x[1] * x[1] - 1.0
        return a * a + b * b + 0.1 * x[0] + 0.2 * x[1]

    def sym(x):
        a = x[0] * x[0] - 1.0
        b = x[1] * x[1] - 1.0
        return a * a + b * b

    def wells3(x):
        a = x[0] * x[0] - 1.0
        b = x[1] * x[1] - 1.0
        c = x[2] - 0.5
        return a * a + b * b + c * c * c * c + 0.1 * x[0] + 0.3 * x[1] * x[2]

    def quart1(x):
        return x[0] * x[0] * x[0] * x[0] - 3.0 * x[0] * x[0] * x[0] + 2.0 * x[0] * x[0] + 0.5 * x[0]

    def shrun(name, fun, bounds, **kw):
        a = {"n": 100, "iters": 1, "sampling_method": "simplicial"}
        a.update(kw)
        r = so.shgo(fun, bounds, **a)
        parts = [f"x = {rv(r.x) if r.x is not None else 'NULL'}",
                 f"fun = {repr(float(r.fun)) if r.fun is not None else 'NULL'}"]
        for f in ("nfev", "nlfev", "nljev", "nit"):
            parts.append(f"{f} = {int(getattr(r, f))}L")
        parts.append(f"success = {'TRUE' if r.success else 'FALSE'}")
        if "funl" in r:
            parts.append(f"funl = {rv(r.funl)}")
        print(f"ref${name} <- list({', '.join(parts)})")
        print(f"ref${name}$message <- {r.message!r}")

    w2 = [(-2.0, 2.0), (-2.0, 2.0)]
    w3 = [(-2.0, 2.0), (-2.0, 2.0), (-1.0, 2.0)]
    rb = [(-2.0, 2.0), (-1.0, 3.0)]
    cons_w = [{"type": "ineq", "fun": lambda x: 1.5 - x[0] - x[1]}]
    cons_r = [{"type": "ineq", "fun": lambda x: 2.0 - x[0] * x[0] - x[1] * x[1]}]
    shrun("sh_rosen", rosen, rb)
    shrun("sh_poly3", poly3, [(-1.0, 2.5), (0.0, 3.0), (-2.0, 1.0)])
    shrun("sh_wells", wells, w2)
    shrun("sh_sym4", sym, w2, iters=4)
    shrun("sh_quart5", quart1, [(-1.0, 3.0)], iters=5)
    shrun("sh_wells3", wells3, w3, n=40, iters=3)
    shrun("sh_rosen_none", rosen, rb, n=None, iters=3)
    shrun("sh_cons", wells, w2, constraints=cons_w, iters=2)
    shrun("sh_cons_r", rosen, rb, constraints=cons_r, iters=3)
    shrun("sh_eq", wells, w2, constraints=[{"type": "eq", "fun": lambda x: x[0] - x[1]}])
    shrun("sh_infeasible", wells, w2, constraints=[{"type": "ineq", "fun": lambda x: x[0] + x[1] - 10.0}])
    shrun("sh_maxfev", wells, w2, options={"maxfev": 40})
    shrun("sh_maxiter", wells, w2, options={"maxiter": 3})
    shrun("sh_maxev", wells, w2, options={"maxev": 30})
    shrun("sh_fmin", wells, w2, options={"f_min": -0.3, "f_tol": 1e-3})
    shrun("sh_minhgrd", wells, w2, n=10, options={"minhgrd": 0})
    shrun("sh_local", wells, w2, iters=2, options={"local_iter": 1})
    shrun("sh_every", wells, w2, iters=3, options={"minimize_every_iter": False})
    shrun("sh_lbfgsb", wells, w2, iters=2, minimizer_kwargs={"method": "L-BFGS-B"})
    shrun("sh_symmetry", sym, w2, options={"symmetry": True})
    calls = []
    shrun("sh_callback", wells, w2, iters=2, callback=lambda x: calls.append(np.copy(x)))
    print(f"ref$sh_callback$ncalls <- {len(calls)}L")
    print(f"ref$sh_callback$last <- {rv(calls[-1])}")
    # the parallel evaluation of the sampling points (the map of a pool of
    # two forked processes)
    with multiprocessing.get_context("fork").Pool(2) as pool:
        shrun("sh_workers", wells, w2, iters=2, workers=pool.map)

if "trustconstr" in sections:
    # trust-constr with the options of OptimizationSciPy.jl (initial_tr_radius
    # 1, finite_diff_rel_step 1e-8, gtol 1e-10, maxiter 50000) and the BFGS or
    # SR1 updates; the constraints are NonlinearConstraint objects with the
    # finite_diff_rel_step of the wrapper. Equality constraints alone run the
    # SQP method, inequalities and bounds the interior point method.
    import io
    import contextlib
    import warnings

    tc_opts = {"initial_tr_radius": 1.0, "verbose": 0, "finite_diff_rel_step": 1e-8, "gtol": 1e-10,
               "maxiter": 50000}
    tc_x0 = np.array([-1.2, 1.0])
    tc_box = so.Bounds([-2.0, -2.0], [2.0, 0.8])

    def circ(x):
        return np.array([x[0] * x[0] + x[1] * x[1]])

    def circ_jac(x):
        return np.array([[2.0 * x[0], 2.0 * x[1]]])

    def two(x):
        return np.array([x[0] * x[0] + x[1] * x[1], x[0] - x[1]])

    def two_jac(x):
        return np.array([[2.0 * x[0], 2.0 * x[1]], [1.0, -1.0]])

    def dup(x):
        return np.array([x[0] + x[1], x[0] + x[1]])

    def dup_jac(x):
        return np.array([[1.0, 1.0], [1.0, 1.0]])

    def f3(x):
        a = x[0] - 1.0
        b = x[1] - 2.0
        c = x[2] - 3.0
        return a * a + b * b + c * c + x[0] * x[1] * x[2]

    def f3_grad(x):
        return np.array([2.0 * (x[0] - 1.0) + x[1] * x[2], 2.0 * (x[1] - 2.0) + x[0] * x[2],
                         2.0 * (x[2] - 3.0) + x[0] * x[1]])

    def g3(x):
        return np.array([x[0] + x[1] + x[2], x[0] * x[1]])

    def g3_jac(x):
        return np.array([[1.0, 1.0, 1.0], [x[1], x[0], 0.0]])

    def tnlc(fun, lb, ub, jac="2-point", **kw):
        return so.NonlinearConstraint(fun, lb, ub, jac=jac, finite_diff_rel_step=1e-8, **kw)

    def tcrun(name, fun, x0, options=None, **kw):
        o = dict(tc_opts)
        o.update(options or {})
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r = so.minimize(fun, x0, method="trust-constr", options=o, **kw)
        parts = [f"x = {rv(r.x)}", f"fun = {repr(float(r.fun))}"]
        for f in ("nit", "nfev", "njev", "nhev", "cg_niter", "status"):
            parts.append(f"{f} = {int(getattr(r, f))}L")
        cn = ", ".join(str(int(t)) + "L" for t in r.constr_nfev)
        cj = ", ".join(str(int(t)) + "L" for t in r.constr_njev)
        parts.append(f"constr_nfev = c({cn})")
        parts.append(f"constr_njev = c({cj})")
        parts.append(f"nwarn = {len(w)}L")
        print(f"ref${name} <- list({', '.join(parts)})")
        print(f"ref${name}$message <- {r.message!r}")

    tcrun("tc_none", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS())
    tcrun("tc_hess", rosen, tc_x0, jac=rosen_grad, hess=rosen_hess)
    tcrun("tc_sr1", rosen, tc_x0, jac=rosen_grad, hess=so.SR1())
    tcrun("tc_fd", rosen, tc_x0, hess=so.BFGS())
    tcrun("tc_maxiter", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(), options={"maxiter": 5})
    tcrun("tc_bounds", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(), bounds=tc_box)
    tcrun("tc_ineq", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(),
          constraints=[tnlc(circ, -np.inf, 1.5, circ_jac)])
    tcrun("tc_ineq_fd", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(), constraints=[tnlc(circ, -np.inf, 1.5)])
    tcrun("tc_greater", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(),
          constraints=[tnlc(circ, 0.5, np.inf, circ_jac)])
    tcrun("tc_eq", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(), constraints=[tnlc(circ, 1.5, 1.5, circ_jac)])
    tcrun("tc_both", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(), bounds=tc_box,
          constraints=[tnlc(circ, -np.inf, 1.5, circ_jac)])
    tcrun("tc_mixed", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(),
          constraints=[tnlc(two, [0.5, 0.0], [1.5, 0.0], two_jac)])
    tcrun("tc_mixed2", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(),
          constraints=[tnlc(two, [0.5, -np.inf], [1.5, 0.3], two_jac)])
    tcrun("tc_keep", rosen, np.array([0.5, 0.5]), jac=rosen_grad, hess=so.BFGS(),
          constraints=[tnlc(circ, -np.inf, 1.5, circ_jac, keep_feasible=True)])
    tcrun("tc_svd", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(),
          constraints=[tnlc(circ, -np.inf, 1.5, circ_jac)], options={"factorization_method": "SVDFactorization"})
    tcrun("tc_singular", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(),
          constraints=[tnlc(dup, 1.0, 1.0, dup_jac)])
    tcrun("tc_linear", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(),
          constraints=[so.LinearConstraint([[1.0, 1.0]], -np.inf, 1.0)])
    tcrun("tc_three", f3, np.array([0.5, 0.5, 0.5]), jac=f3_grad, hess=so.BFGS(),
          bounds=so.Bounds([-1.0, -1.0, -1.0], [2.0, 2.0, 2.0]),
          constraints=[tnlc(g3, [1.0, -np.inf], [1.0, 0.5], g3_jac)])
    tcrun("tc_three_eq", f3, np.array([0.5, 0.5, 0.5]), jac=f3_grad, hess=so.BFGS(),
          constraints=[tnlc(g3, [1.0, 0.2], [1.0, 0.2], g3_jac)])
    # finite difference Hessians (linear operators of gradient or Jacobian
    # differences), a Hessian vector product and a constraint Hessian
    def rosen_hessp(x, p):
        h11 = 2.0 - 400.0 * x[1] + 1200.0 * x[0] * x[0]
        h12 = -400.0 * x[0]
        return np.array([h11 * p[0] + h12 * p[1], h12 * p[0] + 200.0 * p[1]])

    def circ_hess(x, v):
        return np.array([[2.0 * v[0], 0.0], [0.0, 2.0 * v[0]]])

    tcrun("tc_hess_fd", rosen, tc_x0, jac=rosen_grad, hess="2-point")
    tcrun("tc_hess_fd3", rosen, tc_x0, jac=rosen_grad, hess="3-point")
    tcrun("tc_hessp", rosen, tc_x0, jac=rosen_grad, hessp=rosen_hessp)
    tcrun("tc_con_hess", rosen, tc_x0, jac=rosen_grad, hess=rosen_hess,
          constraints=[tnlc(circ, -np.inf, 1.5, circ_jac, hess=circ_hess)])
    tcrun("tc_con_hess_fd", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(),
          constraints=[tnlc(circ, -np.inf, 1.5, circ_jac, hess="2-point")])

    # a callback that raises StopIteration at the third iteration
    def tc_stop(intermediate_result):
        if intermediate_result.nit == 3:
            raise StopIteration

    tcrun("tc_stop", rosen, tc_x0, jac=rosen_grad, hess=so.BFGS(), callback=tc_stop)
    # the report of verbose = 2 and 3 (the last line holds the run time)
    for v in (2, 3):
        for nm, kw in (("sqp", {}), ("ip", {"constraints": [tnlc(circ, -np.inf, 1.5, circ_jac)]})):
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                so.minimize(rosen, tc_x0, method="trust-constr", jac=rosen_grad, hess=so.BFGS(),
                            options={"initial_tr_radius": 1.0, "verbose": v, "maxiter": 3}, **kw)
            lines = buf.getvalue().splitlines()[:-1]
            print(f"ref$tc_report_{nm}{v} <- c({', '.join(repr(t) for t in lines)})")

    # the errors of SciPy
    def tc_ftol():
        so.minimize(rosen, tc_x0, method="trust-constr", jac=rosen_grad, options={"ftol": 1e-8})

    def tc_factorization():
        so.minimize(rosen, tc_x0, method="trust-constr", jac=rosen_grad, hess=so.BFGS(),
                    constraints=[tnlc(circ, 1.5, 1.5, circ_jac)], options={"factorization_method": "LU"})

    def tc_fd_both():
        so.minimize(rosen, tc_x0, method="trust-constr", hess="2-point")

    for nm, call in (("ftol", tc_ftol), ("factorization", tc_factorization), ("fd_both", tc_fd_both)):
        try:
            call()
            print(f"ref$tc_err_{nm} <- NULL")
        except Exception as e:
            print(f"ref$tc_err_{nm} <- {str(e)!r}")

if "cobyqa" in sections:
    # COBYQA through minimize(method="COBYQA") with the defaults of
    # _minimize_cobyqa. The functions are objects without __name__, so the
    # report of disp=True names the objective "fun" and leaves out the
    # constraint values, as the R port does.
    import io
    import contextlib
    import warnings

    class CqF:
        def __init__(self, f):
            self.f = f

        def __call__(self, x):
            return self.f(x)

    def cq_quad(x):
        return (x[0] - 1.0) * (x[0] - 1.0) + 2.0 * (x[1] + 0.5) * (x[1] + 0.5) + x[0] * x[1]

    def cq_rosen(x):
        return (1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0])

    def cq_poly3(x):
        return x[0] * x[0] * x[0] * x[0] + x[1] * x[1] + x[2] * x[2] * x[2] * x[2] - x[0] * x[2] + x[1] - 2.0 * x[2]

    def cq_circ(x):
        return np.array([x[0] * x[0] + x[1] * x[1]])

    def cq_two(x):
        return np.array([x[0] * x[0] + x[1] * x[1], x[0] - x[1]])

    cnl = so.NonlinearConstraint
    clc = so.LinearConstraint
    cq_cases = {
        "cq_quad": (cq_quad, [0.5, 0.5], {}),
        "cq_rosen": (cq_rosen, [-1.2, 1.0], {}),
        "cq_poly3": (cq_poly3, [0.5, 0.5, 0.5], {}),
        "cq_bounds": (cq_rosen, [-1.2, 1.0], {"bounds": so.Bounds([-2.0, -0.5], [0.8, 2.0])}),
        "cq_scale": (cq_quad, [0.5, 0.5], {"bounds": so.Bounds([-2.0, -3.0], [2.0, 1.0]), "options": {"scale": True}}),
        "cq_fixed": (cq_poly3, [0.5, 0.5, 0.5], {"bounds": so.Bounds([-2.0, 0.3, -2.0], [2.0, 0.3, 2.0])}),
        "cq_linear": (cq_quad, [0.5, 0.5], {"constraints": [clc([[1.0, 1.0]], -np.inf, 0.2)]}),
        "cq_lineq": (cq_poly3, [0.5, 0.5, 0.5], {"constraints": [clc([[1.0, 1.0, 1.0]], 1.0, 1.0)]}),
        "cq_ineq": (cq_rosen, [-1.2, 1.0], {"constraints": [cnl(CqF(cq_circ), -np.inf, 1.5)]}),
        "cq_greater": (cq_quad, [0.5, 0.5], {"constraints": [cnl(CqF(cq_circ), 2.0, np.inf)]}),
        "cq_eq": (cq_quad, [0.5, 0.5], {"constraints": [cnl(CqF(cq_circ), 1.0, 1.0)]}),
        "cq_mixed": (cq_rosen, [-1.2, 1.0], {"bounds": so.Bounds([-2.0, -2.0], [2.0, 0.9]),
                                             "constraints": [cnl(CqF(cq_two), [0.5, -np.inf], [1.5, 0.3])]}),
        "cq_maxfev": (cq_rosen, [-1.2, 1.0], {"options": {"maxfev": 20}}),
        "cq_maxiter": (cq_rosen, [-1.2, 1.0], {"options": {"maxiter": 5}}),
        "cq_target": (cq_quad, [0.5, 0.5], {"options": {"f_target": 0.5}}),
        "cq_radius": (cq_quad, [0.5, 0.5], {"options": {"initial_tr_radius": 0.1, "final_tr_radius": 1e-8}}),
        "cq_smallbox": (cq_quad, [0.5, 0.5], {"bounds": so.Bounds([0.0, -0.6], [0.4, 0.6])}),
        "cq_initfev": (cq_quad, [0.5, 0.5], {"options": {"maxfev": 3}}),
        # Problem.violation evaluates the constraints at the reduced point:
        # these two runs are wrong in SciPy (see PORTING.md)
        "cq_scale_nl": (cq_quad, [0.5, 0.5], {"bounds": so.Bounds([-2.0, -3.0], [2.0, 1.0]),
                                              "constraints": [cnl(CqF(cq_circ), -np.inf, 1.0)],
                                              "options": {"scale": True}}),
        "cq_fixed_nl": (cq_poly3, [0.5, 0.5, 0.5], {"bounds": so.Bounds([-2.0, 0.3, -2.0], [2.0, 0.3, 2.0]),
                                                    "constraints": [cnl(CqF(cq_circ), -np.inf, 0.5)]}),
    }

    def cqrun(name, fun, x0, kw, **extra):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r = so.minimize(CqF(fun), x0, method="COBYQA", **kw, **extra)
        parts = [f"x = {rv(r.x)}", f"fun = {repr(float(r.fun))}", f"maxcv = {repr(float(r.maxcv))}",
                 f"nfev = {int(r.nfev)}L", f"nit = {int(r.nit)}L", f"status = {int(r.status)}L",
                 f"success = {'TRUE' if r.success else 'FALSE'}", f"nwarn = {len(w)}L"]
        print(f"ref${name} <- list({', '.join(parts)})")
        print(f"ref${name}$message <- {r.message!r}")

    for nm, (fun, x0, kw) in cq_cases.items():
        cqrun(nm, fun, x0, kw)

    # callbacks: the best point after each evaluation, StopIteration ends the run
    cq_seen = []

    def cq_cb_x(x):
        cq_seen.append(np.array(x))
        if len(cq_seen) == 9:
            raise StopIteration

    def cq_cb_ir(intermediate_result):
        if intermediate_result.fun < -0.7:
            raise StopIteration

    cqrun("cq_cb_x", cq_quad, [0.5, 0.5], {}, callback=cq_cb_x)
    print(f"ref$cq_cb_x$seen <- {rv(np.concatenate(cq_seen))}")
    cqrun("cq_cb_ir", cq_quad, [0.5, 0.5], {}, callback=cq_cb_ir)

    # the report of disp=True up to the first trust region step
    def cq_sq7(x):
        d = x - 0.3
        s = d[0] * d[0]
        for i in range(1, 7):
            s = s + d[i] * d[i]
        return s

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        so.minimize(CqF(cq_quad), [0.5, 0.5], method="COBYQA", options={"disp": True, "maxiter": 12},
                    constraints=[cnl(CqF(cq_circ), -np.inf, 1.0)])
    lines = buf.getvalue().splitlines()[:18]
    print(f"ref$cq_report <- c({', '.join(repr(t) for t in lines)})")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        so.minimize(CqF(cq_sq7), np.arange(7.0) / 10 - 2e-5, method="COBYQA",
                    options={"disp": True, "maxiter": 1})
    lines = buf.getvalue().splitlines()[:10]
    print(f"ref$cq_report7 <- c({', '.join(repr(t) for t in lines)})")

    # the errors and warnings
    def cq_err(nm, call):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                call()
                print(f"ref$cq_err_{nm} <- NULL")
            except Exception as e:
                print(f"ref$cq_err_{nm} <- {str(e)!r}")
        print(f"ref$cq_warn_{nm} <- c({', '.join(repr(str(t.message)) for t in w)})")

    def cq_go(**kw):
        return lambda: so.minimize(CqF(cq_quad), [0.5, 0.5], method="COBYQA", **kw)

    cq_err("radius", cq_go(options={"initial_tr_radius": 0.0}))
    cq_err("final", cq_go(options={"final_tr_radius": -1.0}))
    cq_err("order", cq_go(options={"initial_tr_radius": 1e-3, "final_tr_radius": 1e-2}))
    cq_err("maxfev", cq_go(options={"maxfev": 0}))
    cq_err("maxiter", cq_go(options={"maxiter": 0}))
    cq_err("allfixed", cq_go(bounds=so.Bounds([1.0, 2.0], [1.0, 2.0])))
    cq_err("unknown", cq_go(options={"maxiter": 2, "foo": 1, "bar": 2}))
