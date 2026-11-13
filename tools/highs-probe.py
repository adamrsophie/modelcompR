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
    lo.update(opts)
    try:
        res = so.linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method=method, options=lo)
        show(name, res, f)
    except Exception as err:
        print("==", name)
        print("error", repr(err))


def milp(name, c, integrality, lb, ub, A=None, lbc=None, ubc=None, dev=0, **opts):
    f, lo = logopts(name, dev)
    lo.update(opts)
    cons = () if A is None else so.LinearConstraint(np.atleast_2d(A), lbc, ubc)
    try:
        res = so.milp(c, integrality=integrality, bounds=so.Bounds(lb, ub), constraints=cons, options=lo)
        show(name, res, f)
    except Exception as err:
        print("==", name)
        print("error", repr(err))


inf = np.inf
# the layer matrix: bare (presolve, scaling, perturbation off, Dantzig
# pricing), then one layer on at a time, then all on with presolve off, then
# the defaults
BARE = {"presolve": False, "simplex_scale_strategy": 0, "dual_simplex_cost_perturbation_multiplier": 0.0,
        "simplex_dual_edge_weight_strategy": "dantzig"}
LAYERS = [
    ("bare", dict(BARE)),
    ("scale", dict(BARE, simplex_scale_strategy=1)),
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

for name, kw in LPS:
    for lname, opts in LAYERS:
        lp(f"{name}_ds_{lname}", method="highs-ds", dev=3 if lname in ("bare", "nopresolve", "default") else 0,
           **kw, **opts)
    lp(f"{name}_highs", method="highs", **kw)
    lp(f"{name}_ipm", method="highs-ipm", dev=3, **kw)
    lp(f"{name}_ipm_nopresolve", method="highs-ipm", dev=3, presolve=False, **kw)
    lp(f"{name}_ipm_nocrossover", method="highs-ipm", dev=0, run_crossover="off", **kw)

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
    milp(f"{name}_default", dev=1, **kw)
    milp(f"{name}_nopresolve", dev=1, presolve=False, **kw)
    milp(f"{name}_bare", dev=1, **kw, **MOFF)
milp("milp_knapsack_gap", mip_rel_gap=0.5, **MILPS[4][1])
milp("milp_knapsack_nodes1", node_limit=1, presolve=False, **MILPS[4][1])
