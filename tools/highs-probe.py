# Reference values for the HiGHS port behind ScipyLinprog and ScipyMilp
# (R/optimization-scipy-highs*.R): SciPy 1.18.1 (HiGHS 1.12.0) linprog with
# the methods highs, highs-ds and highs-ipm and milp, on the problems of
# tests/testthat/test-optimization-scipy.R and on LPs and MILPs that
# exercise the dual revised simplex, presolve, the interior point method
# and the branch and bound. Each LP also runs in a layer matrix (presolve,
# scaling, cost perturbation and edge weights switched on one at a time)
# so that a difference in the iteration count points at one layer. Prints
# every number in hex (float.hex) and the HiGHS log of each case.
# Usage (scipy==1.18.1):
#   python3 tools/highs-probe.py > highs-probe.txt 2>&1
import os
import sys
import tempfile
import warnings

import numpy as np
import scipy
import scipy.optimize as so

warnings.simplefilter("ignore")
print("== versions")
print("scipy", scipy.__version__, "numpy", np.__version__)
LOGDIR = tempfile.mkdtemp()


def hx(v):
    if v is None:
        return "NULL"
    a = np.atleast_1d(np.asarray(v, dtype=float))
    return "c(" + ", ".join(float(t).hex() for t in a) + ")"


def show(name, res, logfile):
    print("==", name)
    for k in ("x", "fun", "status", "message", "nit", "crossover_nit", "success",
              "mip_node_count", "mip_dual_bound", "mip_gap"):
        if k in res:
            v = res[k]
            print(k, v if isinstance(v, (str, bool, int, np.integer)) or v is None else hx(v))
    for k in ("ineqlin", "eqlin", "lower", "upper"):
        if k in res and res[k] is not None:
            print(k, "marginals", hx(res[k].marginals), "residual", hx(res[k].residual))
    sys.stdout.flush()
    if os.path.exists(logfile):
        with open(logfile) as fh:
            for line in fh:
                print("log|", line.rstrip())
        os.remove(logfile)
    sys.stdout.flush()


def logopts(name, dev):
    f = os.path.join(LOGDIR, name + ".log")
    return f, {"output_flag": True, "log_to_console": False, "log_file": f, "log_dev_level": dev}


def lp(name, c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=(0, None), method="highs", dev=0, **opts):
    f, lo = logopts(name, dev)
    lo["time_limit"] = 60.0
    lo.update(opts)
    try:
        res = so.linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method=method, options=lo)
        show(name, res, f)
    except Exception as err:
        print("==", name)
        print("error", repr(err))
        sys.stdout.flush()


def milp(name, c, integrality, lb, ub, A=None, lbc=None, ubc=None, dev=0, **opts):
    f, lo = logopts(name, dev)
    lo["time_limit"] = 60.0
    lo.update(opts)
    cons = () if A is None else so.LinearConstraint(np.atleast_2d(A), lbc, ubc)
    try:
        res = so.milp(c, integrality=integrality, bounds=so.Bounds(lb, ub), constraints=cons, options=lo)
        show(name, res, f)
    except Exception as err:
        print("==", name)
        print("error", repr(err))
        sys.stdout.flush()


def raw(name, c, A, lhs, rhs, lb, ub, integrality=None, dev=0, **opts):
    # HiGHS through the bindings SciPy ships, for options that
    # scipy.optimize does not pass on (the cost perturbation multiplier,
    # mip_allow_restart)
    from scipy.optimize._highspy import _core as _h
    f, lo = logopts(name, dev)
    lo["time_limit"] = 60.0
    lo.update(opts)
    print("==", name)
    try:
        A = np.atleast_2d(np.asarray(A, dtype=float)) if A is not None and len(A) else np.zeros((0, len(c)))
        h = _h._Highs()
        for k, val in lo.items():
            if k == "presolve" and isinstance(val, bool):
                val = "on" if val else "off"
            h.setOptionValue(k, val)
        lpm = _h.HighsLp()
        lpm.num_col_ = len(c)
        lpm.num_row_ = A.shape[0]
        lpm.col_cost_ = np.asarray(c, dtype=float)
        lpm.col_lower_ = np.asarray(lb, dtype=float)
        lpm.col_upper_ = np.asarray(ub, dtype=float)
        lpm.row_lower_ = np.asarray(lhs, dtype=float)
        lpm.row_upper_ = np.asarray(rhs, dtype=float)
        lpm.a_matrix_.format_ = _h.MatrixFormat.kColwise
        lpm.a_matrix_.num_col_ = len(c)
        lpm.a_matrix_.num_row_ = A.shape[0]
        start = [0]; index = []; value = []
        for j in range(len(c)):
            for i in range(A.shape[0]):
                if A[i, j] != 0:
                    index.append(i); value.append(A[i, j])
            start.append(len(index))
        lpm.a_matrix_.start_ = np.asarray(start, dtype=np.int32)
        lpm.a_matrix_.index_ = np.asarray(index, dtype=np.int32)
        lpm.a_matrix_.value_ = np.asarray(value, dtype=float)
        if integrality is not None:
            lpm.integrality_ = [_h.HighsVarType(int(t)) for t in integrality]
        h.passModel(lpm)
        h.run()
        info = h.getInfo()
        print("status", h.modelStatusToString(h.getModelStatus()))
        print("nit", info.simplex_iteration_count, "ipm_nit", info.ipm_iteration_count,
              "mip_node_count", info.mip_node_count)
        print("fun", hx(info.objective_function_value))
        sol = h.getSolution()
        if sol.value_valid:
            print("x", hx(sol.col_value))
            print("row_value", hx(sol.row_value))
        if sol.dual_valid:
            print("col_dual", hx(sol.col_dual))
            print("row_dual", hx(sol.row_dual))
        bas = h.getBasis()
        if bas.valid:
            print("col_status", [str(t).split(".")[-1] for t in bas.col_status])
            print("row_status", [str(t).split(".")[-1] for t in bas.row_status])
    except Exception as err:
        print("error", repr(err))
    sys.stdout.flush()
    if os.path.exists(f):
        with open(f) as fh:
            for line in fh:
                print("log|", line.rstrip())
    sys.stdout.flush()


def raw_lp(name, c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=(0, None), dev=0, **opts):
    n = len(c)
    bl = bounds if isinstance(bounds, list) else [bounds] * n
    lb = [-inf if b[0] is None else b[0] for b in bl]
    ub = [inf if b[1] is None else b[1] for b in bl]
    rows = []; lhs = []; rhs = []
    if A_ub is not None:
        rows += [list(r) for r in np.atleast_2d(A_ub)]; lhs += [-inf] * len(b_ub); rhs += list(b_ub)
    if A_eq is not None:
        rows += [list(r) for r in np.atleast_2d(A_eq)]; lhs += list(b_eq); rhs += list(b_eq)
    raw(name, c, rows, lhs, rhs, lb, ub, dev=dev, **opts)


inf = np.inf
# the layer matrix: bare (presolve, scaling, perturbation off, Dantzig
# pricing), then one layer on at a time, then all on with presolve off, then
# the defaults
BARE = {"presolve": False, "simplex_scale_strategy": 0, "dual_simplex_cost_perturbation_multiplier": 0.0,
        "simplex_dual_edge_weight_strategy": 0}
LAYERS = [
    ("bare", dict(BARE)),
    ("scale", dict(BARE, simplex_scale_strategy=2)),
    ("perturb", dict(BARE, dual_simplex_cost_perturbation_multiplier=1.0)),
    ("dse", {k: v for k, v in BARE.items() if k != "simplex_dual_edge_weight_strategy"}),
    ("nopresolve", {"presolve": False}),
    ("default", {}),
]

LPS = []
# the problems of tests/testthat/test-optimization-scipy.R
LPS.append(("t_lp1", dict(c=[-1.0, -2.0, 0.5], A_ub=[[1.0, 1.0, 1.0], [2.0, 0.5, 0.0], [0.0, 1.0, -1.0]],
                          b_ub=[4.0, 5.0, 1.5], A_eq=[[1.0, -1.0, 2.0]], b_eq=[0.5],
                          bounds=[(0, 3), (-1, None), (None, 2)])))
LPS.append(("t_lp2", dict(c=[1.0, 1.0], A_ub=[[-1.0, -2.0], [-3.0, -1.0]], b_ub=[-4.0, -6.0])))
LPS.append(("t_lpmax", dict(c=[-1.0, -1.0], A_ub=[[1.0, 2.0], [3.0, 1.0]], b_ub=[4.0, 6.0])))
LPS.append(("t_lpinf", dict(c=[1.0, 1.0], A_ub=[[1.0, 1.0]], b_ub=[-1.0])))
LPS.append(("t_lpunb", dict(c=[-1.0, 0.0], A_ub=[[1.0, -1.0]], b_ub=[1.0])))
LPS.append(("t_lp3", dict(c=[1.0, -1.0, 2.0, 0.5], A_ub=[[1.0, 1.0, 1.0, 1.0], [1.0, -2.0, 0.0, 1.0]], b_ub=[5.0, 2.0],
                          A_eq=[[0.0, 1.0, 1.0, -1.0]], b_eq=[1.0], bounds=[(1, 4), (-2, 3), (None, None), (0, 2)])))
rng = np.random.default_rng(20261113)
A1 = rng.uniform(-1, 1, (8, 12)); b1 = A1 @ rng.uniform(0, 1, 12) + 0.5
Ae1 = rng.uniform(-1, 1, (3, 12)); be1 = Ae1 @ rng.uniform(0, 1, 12)
c1 = rng.uniform(-1, 1, 12)
LPS.append(("lp_dense", dict(c=c1, A_ub=A1, b_ub=b1, A_eq=Ae1, b_eq=be1, bounds=[(0, 2)] * 12)))
n = 5
ckm = -np.array([2.0 ** (n - 1 - j) for j in range(n)])
Akm = np.zeros((n, n)); bkm = np.array([5.0 ** (i + 1) for i in range(n)])
for i in range(n):
    for j in range(i):
        Akm[i, j] = 2.0 ** (i - j + 1)
    Akm[i, i] = 1.0
LPS.append(("lp_klee_minty", dict(c=ckm, A_ub=Akm, b_ub=bkm)))
LPS.append(("lp_free_fixed", dict(c=[1.0, -1.0, 2.0, 0.5], A_ub=[[1.0, 1.0, 0.0, 1.0], [-1.0, 2.0, 1.0, 0.0]],
                                  b_ub=[4.0, 3.0], A_eq=[[1.0, -1.0, 1.0, 0.0], [2.0, -2.0, 2.0, 0.0]], b_eq=[1.0, 2.0],
                                  bounds=[(None, None), (-1, None), (0, 0), (0, 3)])))
supply = np.array([20.0, 30.0, 25.0]); demand = np.array([10.0, 25.0, 15.0, 25.0])
cost = np.array([[8.0, 6.0, 10.0, 9.0], [9.0, 12.0, 13.0, 7.0], [14.0, 9.0, 16.0, 5.0]])
At = np.zeros((7, 12))
for i in range(3):
    At[i, 4 * i:4 * i + 4] = 1.0
for j in range(4):
    At[3 + j, j::4] = 1.0
LPS.append(("lp_transport", dict(c=cost.ravel(), A_ub=At[:3], b_ub=supply, A_eq=At[3:], b_eq=demand)))
# a larger random LP with free and boxed variables and mixed rows
A2 = rng.uniform(-1, 1, (15, 20)); A2[rng.uniform(0, 1, (15, 20)) < 0.5] = 0.0
x2 = rng.uniform(-1, 2, 20)
b2 = A2 @ x2 + rng.uniform(0, 1, 15)
c2 = rng.uniform(-1, 1, 20)
bnd2 = [(-1, 3)] * 10 + [(None, None)] * 3 + [(0, None)] * 4 + [(None, 2)] * 3
LPS.append(("lp_sparse", dict(c=c2, A_ub=A2[:12], b_ub=b2[:12], A_eq=A2[12:], b_eq=A2[12:] @ x2, bounds=bnd2)))
LPS.append(("lp_infeas_eq", dict(c=[1.0, 2.0, 3.0], A_eq=[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], b_eq=[1.0, 2.0])))
LPS.append(("lp_unb_free", dict(c=[1.0, -1.0, 0.0], A_ub=[[1.0, 1.0, 1.0]], b_ub=[3.0],
                                bounds=[(None, None), (None, None), (0, 1)])))

# the LP data, so that the R port runs the same problems
for name, kw in LPS:
    print("== data", name)
    for k, val in kw.items():
        if k == "bounds":
            bl = val if isinstance(val, list) else [val] * len(kw["c"])
            print("lb", hx([-inf if b[0] is None else b[0] for b in bl]))
            print("ub", hx([inf if b[1] is None else b[1] for b in bl]))
        else:
            a = np.atleast_2d(np.asarray(val, dtype=float))
            print(k, a.shape[0], a.shape[1], hx(a.ravel()))
sys.stdout.flush()

for name, kw in LPS:
    for lname, opts in LAYERS:
        dev = 3 if lname in ("bare", "nopresolve", "default") else 0
        if lname in ("nopresolve", "default"):
            lp(f"{name}_ds_{lname}", method="highs-ds", dev=dev, **kw, **opts)
        else:
            # the dual simplex method as highs-ds selects it
            raw_lp(f"{name}_ds_{lname}", dev=dev, solver="simplex", simplex_strategy=1, **kw, **opts)
    raw_lp(f"{name}_raw_ds", solver="simplex", simplex_strategy=1, **kw)
    raw_lp(f"{name}_raw_ds_nopresolve", solver="simplex", simplex_strategy=1, presolve=False, **kw)
    lp(f"{name}_highs", method="highs", **kw)
    lp(f"{name}_ipm", method="highs-ipm", dev=3, **kw)
    lp(f"{name}_ipm_nopresolve", method="highs-ipm", dev=3, presolve=False, **kw)
    lp(f"{name}_ipm_nocrossover", method="highs-ipm", dev=0, run_crossover="off", **kw)

# a random LP with parallel rows where the dual simplex method without
# presolve ends on a taboo basis in the port
fh = lambda v: [float.fromhex(t) if t not in ("inf", "-inf") else float(t) for t in v]
TABOO = dict(c=fh(["-0x1.999999999999ap-1", "-0x1p-1", "-0x1.b333333333333p+0", "-0x1.2666666666666p+1",
                   "-0x1.3333333333333p-2", "-0x1.999999999999ap-2", "0x1p-1", "-0x1.70a3d70a3d70bp-1"]),
             A=[fh(["0x0p+0", "0x0p+0", "0x0p+0", "-0x1.999999999999ap+0", "0x0p+0", "0x0p+0",
                    "-0x1.3333333333333p-2", "0x0p+0"]),
                fh(["0x0p+0", "-0x1.b333333333333p+0", "0x0p+0", "0x0p+0", "0x0p+0", "-0x1.3333333333333p-1",
                    "0x0p+0", "0x0p+0"])],
             lhs=fh(["-0x1p-1", "0x1.999999999999ap-2"]), rhs=fh(["0x1.f333333333333p+1", "0x1.999999999999ap-2"]),
             lb=fh(["-0x1.3333333333333p-2", "-0x1.3333333333333p-1", "-0x1.e666666666666p+0", "-0x1.ccccccccccccdp-1",
                    "-0x1.6666666666666p-1", "-0x1.4cccccccccccdp+0", "-0x1.0cccccccccccdp+1", "-0x1.ccccccccccccdp-1"]),
             ub=fh(["inf", "inf", "0x1.3333333333333p+0", "inf", "0x1.999999999999ap+1", "inf",
                    "-0x1.0cccccccccccdp+1", "inf"]))
for pre in (False, True):
    raw("lp_taboo_raw_ds" + ("" if pre else "_nopresolve"), TABOO["c"], TABOO["A"], TABOO["lhs"], TABOO["rhs"],
        TABOO["lb"], TABOO["ub"], dev=3, solver="simplex", simplex_strategy=1, presolve=pre)

# small iteration limits: where the stand in stopped with MaxIters
for it in (0, 1, 2, 3, 5):
    lp(f"t_lp1_maxiter{it}", method="highs", maxiter=it, **LPS[0][1])
    lp(f"t_lp1_maxiter{it}_nopresolve", method="highs", maxiter=it, presolve=False, **LPS[0][1])
    lp(f"lp_dense_maxiter{it}", method="highs-ds", maxiter=it, **LPS[6][1])
    lp(f"lp_dense_maxiter{it}_nopresolve", method="highs-ds", maxiter=it, presolve=False, **LPS[6][1])
    lp(f"lp_dense_ipm_maxiter{it}", method="highs-ipm", maxiter=it, **LPS[6][1])

# MILPs: the problems of the tests, a knapsack, a set cover, an infeasible
# one, a gap limit, each with the defaults, presolve off, and presolve,
# heuristics, symmetry detection and restarts off
MOFF = {"presolve": False, "mip_heuristic_effort": 0.0, "mip_detect_symmetry": False, "mip_allow_restart": False,
        "mip_heuristic_run_feasibility_jump": False, "mip_heuristic_run_rins": False, "mip_heuristic_run_rens": False,
        "mip_heuristic_run_root_reduced_cost": False, "mip_heuristic_run_zi_round": False,
        "mip_heuristic_run_shifting": False}
MILPS = []
MILPS.append(("t_milp1", dict(c=np.array([-5.0, -4.0, -3.0, -1.0]), integrality=np.array([1, 1, 1, 0]),
                              lb=np.zeros(4), ub=np.array([10.0, 10.0, 10.0, 2.5]),
                              A=[[2.0, 3.0, 1.0, 1.0], [4.0, 1.0, 2.0, 0.5], [3.0, 4.0, 2.0, 1.0]],
                              lbc=[-inf] * 3, ubc=[5.0, 11.0, 8.0])))
MILPS.append(("t_milp2", dict(c=np.array([-1.0, -1.0]), integrality=np.array([1, 1]), lb=[0.0, 0.0], ub=[inf, inf],
                              A=[[2.0, 2.0], [1.0, -1.0]], lbc=[-inf, 0.0], ubc=[7.0, 0.0])))
MILPS.append(("t_milp3", dict(c=np.array([1.0]), integrality=np.array([1]), lb=[0.2], ub=[0.8])))
MILPS.append(("t_milp_semi", dict(c=np.array([-1.0, -2.0]), integrality=np.array([2, 3]), lb=[2.0, 1.5], ub=[5.0, 4.0],
                                  A=[[1.0, 1.0], [1.0, 3.0]], lbc=[-inf, -inf], ubc=[5.5, 10.0])))
w = np.array([12.0, 2.0, 1.0, 1.0, 4.0, 7.0, 9.0, 3.0]); v = np.array([4.0, 2.0, 1.0, 2.0, 10.0, 7.0, 8.0, 3.0])
MILPS.append(("milp_knapsack", dict(c=-v, integrality=np.ones(8), lb=np.zeros(8), ub=np.ones(8), A=[w],
                                    lbc=[-inf], ubc=[15.0])))
cover = (rng.uniform(0, 1, (10, 12)) < 0.3).astype(float); cover[:, 0] = 1.0
ccov = rng.uniform(1, 3, 12)
MILPS.append(("milp_cover", dict(c=ccov, integrality=np.ones(12), lb=np.zeros(12), ub=np.ones(12), A=cover,
                                 lbc=np.ones(10), ubc=np.full(10, inf))))
MILPS.append(("milp_general", dict(c=np.array([-5.0, -4.0, -3.0]), integrality=np.array([1, 1, 0]),
                                   lb=np.zeros(3), ub=np.array([10.0, 10.0, 2.5]),
                                   A=[[2.0, 3.0, 1.0], [4.0, 1.0, 2.0], [3.0, 4.0, 2.0]], lbc=[-inf] * 3,
                                   ubc=[5.0, 11.0, 8.0])))
MILPS.append(("milp_infeasible", dict(c=np.array([1.0, 1.0]), integrality=np.ones(2), lb=np.zeros(2),
                                      ub=np.full(2, 10.0), A=[[2.0, 2.0]], lbc=[3.0], ubc=[3.0])))
Ai = rng.integers(-3, 6, (6, 8)).astype(float)
MILPS.append(("milp_int8", dict(c=-rng.integers(1, 9, 8).astype(float), integrality=np.ones(8), lb=np.zeros(8),
                                ub=np.full(8, 5.0), A=Ai, lbc=np.full(6, -inf), ubc=rng.integers(8, 20, 6).astype(float))))
for name, kw in MILPS:
    print("== data", name)
    for k, val in kw.items():
        if val is not None:
            a = np.atleast_2d(np.asarray(val, dtype=float))
            print(k, a.shape[0], a.shape[1], hx(a.ravel()))
sys.stdout.flush()
for name, kw in MILPS:
    milp(f"{name}_default", dev=1, **kw)
    milp(f"{name}_nopresolve", dev=1, presolve=False, **kw)
    raw(f"{name}_bare", kw["c"], kw.get("A"), kw.get("lbc", []), kw.get("ubc", []), kw["lb"], kw["ub"],
        integrality=kw["integrality"], dev=1, **MOFF)
milp("milp_knapsack_gap", mip_rel_gap=0.5, **MILPS[4][1])
milp("milp_knapsack_nodes1", node_limit=1, presolve=False, **MILPS[4][1])
