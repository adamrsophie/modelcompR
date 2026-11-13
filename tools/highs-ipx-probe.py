# Reference values for IPX, the interior point solver of HiGHS 1.12.0, with
# its crossover (R/optimization-scipy-highs-ipx.R): SciPy's
# linprog(method = "highs-ipm") with presolve off, so that IPX sees the LP as
# given, on LPs with equality and inequality rows, boxed, free, fixed and
# upper bounded variables, badly scaled entries (equilibration), more rows
# than twice the columns (IPX dualizes), redundant rows, infeasible,
# unbounded and degenerate problems and iteration limits. Prints the status,
# the solution, the objective, the iteration counts and the marginals in hex
# and the HiGHS log with the IPX debug lines (log_dev_level = 1 sets the IPX
# debug level 2). Usage (scipy==1.18.1):
#   python3 tools/highs-ipx-probe.py > highs-ipx-probe.txt
import ctypes
import sys
import warnings

import numpy as np
import scipy
import scipy.optimize as so

warnings.simplefilter("ignore")
libc = ctypes.CDLL(None)


def flush():
    sys.stdout.flush()
    libc.fflush(None)


print("== versions")
print("scipy", scipy.__version__, "numpy", np.__version__)


def hx(v):
    if v is None:
        return "NULL"
    a = np.atleast_1d(np.asarray(v, dtype=float))
    return "c(" + ", ".join(float(t).hex() for t in a) + ")"


def case(name, c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=(0, None), **opts):
    options = {"presolve": False, "disp": True, "log_dev_level": 1, "time_limit": 60.0}
    options.update(opts)
    print("== case", name)
    print("-- log")
    flush()
    try:
        res = so.linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                         method="highs-ipm", options=options)
    except Exception as err:
        flush()
        print("-- result")
        print("error", repr(err))
        return
    flush()
    print("-- result")
    print("status", res.status)
    print("message", res.message)
    print("nit", res.nit)
    print("crossover_nit", res.get("crossover_nit"))
    print("x", hx(res.x))
    print("fun", hx(res.fun))
    for k in ("ineqlin", "eqlin", "lower", "upper"):
        if k in res and res[k] is not None:
            print(k, "marginals", hx(res[k].marginals), "residual", hx(res[k].residual))
    flush()


rng = np.random.default_rng(20261113)

# equality form: min c'x, A x = b, x >= 0
Aeq = np.array([[1.0, 1.0, 1.0, 0.0, 0.0], [1.0, -1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 0.0, 1.0]])
case("equality", [-1.0, -2.0, 0.0, 0.0, 0.0], A_eq=Aeq, b_eq=[4.0, 1.0, 3.0])
# inequality rows, 3 rows and 2 columns (3 > 2*2 is false: primal)
case("inequality", [-1.0, -2.0], A_ub=[[1.0, 1.0], [-1.0, 1.0], [1.0, -2.0]], b_ub=[4.0, 2.0, 1.0])
# 5 rows and 2 columns: IPX dualizes
case("dualized", [-1.0, -1.0], A_ub=[[1.0, 2.0], [2.0, 1.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
     b_ub=[8.0, 8.0, 3.5, 3.5, 5.0])
# dualized with boxed and free columns
case("dualized_boxed_free", [1.0, -1.0], A_ub=[[1.0, 1.0], [-1.0, 1.0], [1.0, -1.0], [-1.0, -1.0], [2.0, 1.0]],
     b_ub=[3.0, 1.0, 2.0, 1.0, 5.0], bounds=[(None, None), (-1.0, 2.0)])
# boxed variables
case("boxed", [-3.0, -1.0, -2.0], A_ub=[[1.0, 1.0, 3.0], [2.0, 2.0, 5.0], [4.0, 1.0, 2.0]],
     b_ub=[30.0, 24.0, 36.0], bounds=[(0, 5), (1, 4), (0, 6)])
# free variables
case("free", [1.0, 1.0, 2.0], A_ub=[[-1.0, 1.0, 0.0], [1.0, 1.0, 1.0]], b_ub=[1.0, 6.0],
     A_eq=[[1.0, -1.0, 1.0]], b_eq=[1.0], bounds=[(None, None), (None, None), (0, None)])
# a variable with only an upper bound and a fixed variable
case("upper_fixed", [1.0, -1.0, 1.0], A_ub=[[1.0, 1.0, 1.0], [1.0, -1.0, 0.0]], b_ub=[4.0, 1.0],
     A_eq=[[0.0, 1.0, 1.0]], b_eq=[2.0], bounds=[(None, 3.0), (0, None), (0.5, 0.5)])
# entries outside [0.5, 16): equilibration
case("scaled", [-1.0, -1e-3, -50.0], A_ub=[[1e3, 2.0, 1e-2], [3e-3, 4e2, 5.0], [1.0, 1.0, 1e4]],
     b_ub=[1e4, 2e2, 3e5])
# a redundant equality row
case("redundant", [1.0, 2.0, 3.0], A_eq=[[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [1.0, -1.0, 0.0]],
     b_eq=[3.0, 6.0, 0.5])
case("infeasible", [1.0, 1.0], A_ub=[[1.0, 1.0], [-1.0, -1.0]], b_ub=[1.0, -2.0])
case("infeasible_eq", [1.0, 1.0, 1.0], A_eq=[[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 2.0, 1.0]],
     b_eq=[1.0, 1.0, 3.0])
case("unbounded", [-1.0, -1.0], A_ub=[[1.0, -1.0]], b_ub=[1.0])
case("unbounded_free", [-1.0, 0.0, 1.0], A_eq=[[1.0, 1.0, 1.0]], b_eq=[1.0],
     bounds=[(None, None), (None, None), (0, None)])
# degenerate: Klee-Minty cube in 5 dimensions and a degenerate vertex
n = 5
ckm = -np.array([2.0 ** (n - 1 - j) for j in range(n)])
Akm = np.zeros((n, n)); bkm = np.array([5.0 ** (i + 1) for i in range(n)])
for i in range(n):
    for j in range(i):
        Akm[i, j] = 2.0 ** (i - j + 1)
    Akm[i, i] = 1.0
case("klee_minty", ckm, A_ub=Akm, b_ub=bkm)
case("degenerate", [-1.0, -1.0, 0.0], A_ub=[[1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 2.0]],
     b_ub=[1.0, 1.0, 2.0, 2.0])
# a dense LP with equality rows and boxed variables
A1 = rng.uniform(-1, 1, (8, 12)); b1 = A1 @ rng.uniform(0, 1, 12) + 0.5
Ae1 = rng.uniform(-1, 1, (3, 12)); be1 = Ae1 @ rng.uniform(0, 1, 12)
c1 = rng.uniform(-1, 1, 12)
case("dense", c1, A_ub=A1, b_ub=b1, A_eq=Ae1, b_eq=be1, bounds=[(0, 2)] * 12)
# a transportation problem (sparse, many ties)
supply = np.array([20.0, 30.0, 25.0]); demand = np.array([10.0, 25.0, 15.0, 25.0])
cost = np.array([[8.0, 6.0, 10.0, 9.0], [9.0, 12.0, 13.0, 7.0], [14.0, 9.0, 16.0, 5.0]])
At = np.zeros((7, 12))
for i in range(3):
    At[i, 4 * i:4 * i + 4] = 1.0
for j in range(4):
    At[3 + j, j::4] = 1.0
case("transport", cost.ravel(), A_ub=At[:3], b_ub=supply, A_eq=At[3:], b_eq=demand)
# a larger sparse LP: the initial IPM hands over to the basis preconditioned
# IPM once the diagonal preconditioned CR method needs more than 10 + m/20
# iterations
m2, n2 = 40, 70
A2 = rng.uniform(-1, 1, (m2, n2)) * (rng.uniform(0, 1, (m2, n2)) < 0.25)
x2 = rng.uniform(0, 1, n2)
b2 = A2 @ x2 + rng.uniform(0, 1, m2)
c2 = rng.uniform(-1, 1, n2)
case("sparse40", c2, A_ub=A2, b_ub=b2, bounds=[(0, 3)] * n2)
Ae3 = rng.uniform(-1, 1, (25, 60)) * (rng.uniform(0, 1, (25, 60)) < 0.3)
be3 = Ae3 @ rng.uniform(0, 2, 60)
c3 = rng.uniform(0, 1, 60)
case("sparse_eq25", c3, A_eq=Ae3, b_eq=be3)
# iteration limits
for k in (0, 1, 3):
    case(f"dense_maxiter{k}", c1, A_ub=A1, b_ub=b1, A_eq=Ae1, b_eq=be1, bounds=[(0, 2)] * 12, maxiter=k)
for k in (5, 12):
    case(f"sparse40_maxiter{k}", c2, A_ub=A2, b_ub=b2, bounds=[(0, 3)] * n2, maxiter=k)
# the matrices of the random cases, for the R tests
print("== data")
for nm, v in (("A1", A1), ("b1", b1), ("Ae1", Ae1), ("be1", be1), ("c1", c1), ("A2", A2), ("b2", b2),
              ("c2", c2), ("Ae3", Ae3), ("be3", be3), ("c3", c3)):
    print(nm, hx(np.ravel(v, order="F")))
flush()
