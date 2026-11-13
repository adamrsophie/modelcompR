# Reference values for the analytic centre modes of IPX in HiGHS 1.12.0
# (R/optimization-scipy-highs-ipx.R): the run_centring option (centring
# steps after the IPM, no crossover) and the analytic centre computation of
# the MIP solver (HighsMipSolverData::startAnalyticCenterComputation: zero
# costs, presolve off, ipm_iteration_limit = 200, run_crossover = off).
# SciPy's linprog does not pass these options, so the probe uses highspy,
# built from the same HiGHS release. The LPs are those of
# tools/highs-ipx-probe.py (same random draws). A first case repeats a
# SciPy probe case with the default options, to check that both builds
# agree. Prints the IPX iteration log, the model status, the iteration
# count and the primal values in hex.
# Usage (highspy==1.12.0, numpy):
#   python3 tools/highs-ipx-centring-probe.py > highs-ipx-centring-probe.txt
import ctypes
import sys

import highspy
import numpy as np

libc = ctypes.CDLL(None)


def flush():
    sys.stdout.flush()
    libc.fflush(None)


print("== versions")
print("highspy", highspy.Highs().version(), "numpy", np.__version__)


def hx(v):
    a = np.atleast_1d(np.asarray(v, dtype=float))
    return "c(" + ", ".join(float(t).hex() for t in a) + ")"


def lp_of(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=(0, None)):
    c = np.asarray(c, dtype=float)
    n = len(c)
    Aub = np.zeros((0, n)) if A_ub is None else np.atleast_2d(np.asarray(A_ub, dtype=float))
    Aeq = np.zeros((0, n)) if A_eq is None else np.atleast_2d(np.asarray(A_eq, dtype=float))
    A = np.vstack([Aub, Aeq])
    rl = np.concatenate([np.full(Aub.shape[0], -np.inf), np.asarray(b_eq if b_eq is not None else [], dtype=float)])
    ru = np.concatenate([np.asarray(b_ub if b_ub is not None else [], dtype=float),
                         np.asarray(b_eq if b_eq is not None else [], dtype=float)])
    if isinstance(bounds, tuple):
        bounds = [bounds] * n
    lo = np.array([-np.inf if b[0] is None else b[0] for b in bounds], dtype=float)
    up = np.array([np.inf if b[1] is None else b[1] for b in bounds], dtype=float)
    return c, A, rl, ru, lo, up


def case(name, lp, zero_cost=False, **opts):
    c, A, rl, ru, lo, up = lp
    if zero_cost:
        c = np.zeros_like(c)
    print("== case", name)
    print("-- log")
    flush()
    h = highspy.Highs()
    h.setOptionValue("output_flag", True)
    h.setOptionValue("log_dev_level", 1)
    h.setOptionValue("presolve", "off")
    h.setOptionValue("solver", "ipm")
    h.setOptionValue("time_limit", 60.0)
    for k, v in opts.items():
        h.setOptionValue(k, v)
    m, n = A.shape
    model = highspy.HighsLp()
    model.num_col_ = n
    model.num_row_ = m
    model.col_cost_ = c
    model.col_lower_ = lo
    model.col_upper_ = up
    model.row_lower_ = rl
    model.row_upper_ = ru
    start, index, value = [0], [], []
    for j in range(n):
        for i in range(m):
            if A[i, j] != 0.0:
                index.append(i)
                value.append(A[i, j])
        start.append(len(index))
    model.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    model.a_matrix_.start_ = start
    model.a_matrix_.index_ = index
    model.a_matrix_.value_ = value
    h.passModel(model)
    h.run()
    flush()
    info = h.getInfo()
    sol = h.getSolution()
    print("-- result")
    print("status", h.modelStatusToString(h.getModelStatus()))
    print("ipm_iteration_count", info.ipm_iteration_count)
    print("crossover_iteration_count", info.crossover_iteration_count)
    if sol.value_valid:
        print("x", hx(sol.col_value))
        print("row_value", hx(sol.row_value))
    if sol.dual_valid:
        print("col_dual", hx(sol.col_dual))
        print("row_dual", hx(sol.row_dual))
    print("objective", float(info.objective_function_value).hex())
    flush()


rng = np.random.default_rng(20261113)
Aeq = np.array([[1.0, 1.0, 1.0, 0.0, 0.0], [1.0, -1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 0.0, 1.0]])
equality = lp_of([-1.0, -2.0, 0.0, 0.0, 0.0], A_eq=Aeq, b_eq=[4.0, 1.0, 3.0])
boxed = lp_of([-3.0, -1.0, -2.0], A_ub=[[1.0, 1.0, 3.0], [2.0, 2.0, 5.0], [4.0, 1.0, 2.0]],
              b_ub=[30.0, 24.0, 36.0], bounds=[(0, 5), (1, 4), (0, 6)])
degenerate = lp_of([-1.0, -1.0, 0.0], A_ub=[[1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 2.0]],
                   b_ub=[1.0, 1.0, 2.0, 2.0])
A1 = rng.uniform(-1, 1, (8, 12)); b1 = A1 @ rng.uniform(0, 1, 12) + 0.5
Ae1 = rng.uniform(-1, 1, (3, 12)); be1 = Ae1 @ rng.uniform(0, 1, 12)
c1 = rng.uniform(-1, 1, 12)
dense = lp_of(c1, A_ub=A1, b_ub=b1, A_eq=Ae1, b_eq=be1, bounds=[(0, 2)] * 12)
supply = np.array([20.0, 30.0, 25.0]); demand = np.array([10.0, 25.0, 15.0, 25.0])
cost = np.array([[8.0, 6.0, 10.0, 9.0], [9.0, 12.0, 13.0, 7.0], [14.0, 9.0, 16.0, 5.0]])
At = np.zeros((7, 12))
for i in range(3):
    At[i, 4 * i:4 * i + 4] = 1.0
for j in range(4):
    At[3 + j, j::4] = 1.0
transport = lp_of(cost.ravel(), A_ub=At[:3], b_ub=supply, A_eq=At[3:], b_eq=demand)
m2, n2 = 40, 70
A2 = rng.uniform(-1, 1, (m2, n2)) * (rng.uniform(0, 1, (m2, n2)) < 0.25)
x2 = rng.uniform(0, 1, n2)
b2 = A2 @ x2 + rng.uniform(0, 1, m2)
c2 = rng.uniform(-1, 1, n2)
sparse40 = lp_of(c2, A_ub=A2, b_ub=b2, bounds=[(0, 3)] * n2)
Ae3 = rng.uniform(-1, 1, (25, 60)) * (rng.uniform(0, 1, (25, 60)) < 0.3)
be3 = Ae3 @ rng.uniform(0, 2, 60)
c3 = rng.uniform(0, 1, 60)
sparse_eq25 = lp_of(c3, A_eq=Ae3, b_eq=be3)

# the same build as SciPy's?
case("check_dense", dense)
lps = (("equality", equality), ("boxed", boxed), ("degenerate", degenerate), ("dense", dense),
       ("transport", transport), ("sparse40", sparse40), ("sparse_eq25", sparse_eq25))
for nm, lp in lps:
    case("centring_" + nm, lp, run_centring=True)
for nm, lp in lps:
    case("mipac_" + nm, lp, zero_cost=True, ipm_iteration_limit=200, run_crossover="off")
for nm, lp in lps:
    case("acentring_" + nm, lp, zero_cost=True, run_centring=True)
case("centring_steps1_dense", dense, run_centring=True, max_centring_steps=1)
case("centring_tol1_dense", dense, run_centring=True, centring_ratio_tolerance=1.0)
