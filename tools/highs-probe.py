# Reference values for HiGHS behind ScipyLinprog and ScipyMilp
# (R/optimization-scipy-linprog.R): SciPy's linprog and milp with the HiGHS
# methods on LPs that exercise the dual revised simplex (degenerate,
# unbounded, infeasible, free and boxed variables, equality rows), presolve
# on and off, and MILPs whose branch and bound needs several nodes. Prints
# the solution, status, message, iteration and node counts and the
# marginals with full precision. Usage (scipy==1.18.1):
#   python3 tools/highs-probe.py > highs-probe.txt
import numpy as np
import scipy
import scipy.optimize as so

print("== versions")
print("scipy", scipy.__version__, "numpy", np.__version__)


def rv(v):
    if v is None:
        return "NULL"
    a = np.atleast_1d(np.asarray(v, dtype=float))
    return "c(" + ", ".join(repr(float(t)) for t in a) + ")"


def show(name, res):
    print("==", name)
    for k in ("x", "fun", "status", "message", "nit", "success", "mip_node_count", "mip_dual_bound", "mip_gap"):
        if k in res:
            v = res[k]
            print(k, v if isinstance(v, (str, bool, int, np.integer)) or v is None else rv(v))
    for k in ("ineqlin", "eqlin", "lower", "upper"):
        if k in res and res[k] is not None:
            print(k, "marginals", rv(res[k].marginals), "residual", rv(res[k].residual))


def case(name, f):
    try:
        show(name, f())
    except Exception as err:
        print("==", name)
        print("error", repr(err))


rng = np.random.default_rng(20261113)

# a dense LP with equality rows and boxed variables
A1 = rng.uniform(-1, 1, (8, 12)); b1 = A1 @ rng.uniform(0, 1, 12) + 0.5
Ae1 = rng.uniform(-1, 1, (3, 12)); be1 = Ae1 @ rng.uniform(0, 1, 12)
c1 = rng.uniform(-1, 1, 12)
for meth in ("highs", "highs-ds", "highs-ipm"):
    for pre in (True, False):
        case(f"lp_dense_{meth}_presolve{pre}",
             lambda meth=meth, pre=pre: so.linprog(c1, A_ub=A1, b_ub=b1, A_eq=Ae1, b_eq=be1,
                                                   bounds=[(0, 2)] * 12, method=meth,
                                                   options={"presolve": pre}))
# a degenerate LP (Klee-Minty cube in 5 dimensions)
n = 5
ckm = -np.array([2.0 ** (n - 1 - j) for j in range(n)])
Akm = np.zeros((n, n)); bkm = np.array([5.0 ** (i + 1) for i in range(n)])
for i in range(n):
    for j in range(i):
        Akm[i, j] = 2.0 ** (i - j + 1)
    Akm[i, i] = 1.0
for pre in (True, False):
    case(f"lp_klee_minty_presolve{pre}",
         lambda pre=pre: so.linprog(ckm, A_ub=Akm, b_ub=bkm, method="highs-ds", options={"presolve": pre}))
# free variables, a redundant equality and a fixed variable
case("lp_free_fixed", lambda: so.linprog([1.0, -1.0, 2.0, 0.5], A_ub=[[1.0, 1.0, 0.0, 1.0], [-1.0, 2.0, 1.0, 0.0]],
                                         b_ub=[4.0, 3.0], A_eq=[[1.0, -1.0, 1.0, 0.0], [2.0, -2.0, 2.0, 0.0]],
                                         b_eq=[1.0, 2.0], bounds=[(None, None), (-1, None), (0, 0), (0, 3)],
                                         method="highs"))
case("lp_infeasible", lambda: so.linprog([1.0, 1.0], A_ub=[[1.0, 1.0], [-1.0, -1.0]], b_ub=[1.0, -2.0],
                                         method="highs", options={"presolve": False}))
case("lp_unbounded", lambda: so.linprog([-1.0, -1.0], A_ub=[[1.0, -1.0]], b_ub=[1.0],
                                        method="highs", options={"presolve": False}))
case("lp_maxiter", lambda: so.linprog(c1, A_ub=A1, b_ub=b1, A_eq=Ae1, b_eq=be1, bounds=[(0, 2)] * 12,
                                      method="highs-ds", options={"maxiter": 2, "presolve": False}))
# a transportation problem (sparse, many ties)
supply = np.array([20.0, 30.0, 25.0]); demand = np.array([10.0, 25.0, 15.0, 25.0])
cost = np.array([[8.0, 6.0, 10.0, 9.0], [9.0, 12.0, 13.0, 7.0], [14.0, 9.0, 16.0, 5.0]])
At = np.zeros((7, 12))
for i in range(3):
    At[i, 4 * i:4 * i + 4] = 1.0
for j in range(4):
    At[3 + j, j::4] = 1.0
for pre in (True, False):
    case(f"lp_transport_presolve{pre}",
         lambda pre=pre: so.linprog(cost.ravel(), A_ub=At[:3], b_ub=supply, A_eq=At[3:], b_eq=demand,
                                    method="highs", options={"presolve": pre}))

# MILPs: a knapsack, a small set cover, an integer program with a
# fractional LP relaxation at the root, an infeasible one, a time limit
w = np.array([12.0, 2.0, 1.0, 1.0, 4.0, 7.0, 9.0, 3.0]); v = np.array([4.0, 2.0, 1.0, 2.0, 10.0, 7.0, 8.0, 3.0])
for pre in (True, False):
    case(f"milp_knapsack_presolve{pre}",
         lambda pre=pre: so.milp(-v, integrality=np.ones(8), bounds=so.Bounds(0, 1),
                                 constraints=so.LinearConstraint(w[None, :], -np.inf, 15.0),
                                 options={"presolve": pre}))
cover = (rng.uniform(0, 1, (10, 12)) < 0.3).astype(float); cover[:, 0] = 1.0
case("milp_cover", lambda: so.milp(rng.uniform(1, 3, 12), integrality=np.ones(12), bounds=so.Bounds(0, 1),
                                   constraints=so.LinearConstraint(cover, 1.0, np.inf)))
case("milp_general", lambda: so.milp(np.array([-5.0, -4.0, -3.0]), integrality=np.array([1, 1, 0]),
                                     bounds=so.Bounds([0, 0, 0], [10, 10, 2.5]),
                                     constraints=so.LinearConstraint([[2.0, 3.0, 1.0], [4.0, 1.0, 2.0],
                                                                      [3.0, 4.0, 2.0]], -np.inf,
                                                                     [5.0, 11.0, 8.0])))
case("milp_infeasible", lambda: so.milp(np.array([1.0, 1.0]), integrality=np.ones(2), bounds=so.Bounds(0, 10),
                                        constraints=so.LinearConstraint([[2.0, 2.0]], 3.0, 3.0)))
case("milp_gap", lambda: so.milp(-v, integrality=np.ones(8), bounds=so.Bounds(0, 1),
                                 constraints=so.LinearConstraint(w[None, :], -np.inf, 15.0),
                                 options={"mip_rel_gap": 0.5}))
