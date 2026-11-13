# Step level reference traces of the feasibility jump heuristic of HiGHS
# 1.12.0 (mip/feasibilityjump.hh, vendored by SciPy 1.18.1) for the R port in
# R/optimization-scipy-highs-mip-heur.R.
#
# Usage: python3 tools/fj-probe.py <highs-probe output>
# Reads the "== data <name>" blocks of the probe output, builds the input
# HighsMipSolverData::feasibilityJump hands the solver for each MILP with
# presolve off (no incumbent, so variables start at a finite bound or 0),
# compiles tools/fj-harness.cpp and prints, per case, the input and the trace.
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
FEASTOL = 1e-6   # mip_feasibility_tolerance
EPSILON = 1e-9   # small_matrix_value
SEED = 0         # random_seed


def parse(path):
    cases, cur = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("== "):
            parts = line.split()
            cur = parts[2] if len(parts) == 3 and parts[1] == "data" else None
            if cur is not None:
                cases[cur] = {}
            continue
        if cur is None or not line:
            continue
        try:
            key, nr, nc, rest = line.split(" ", 3)
            vals = [float.fromhex(t) for t in rest[2:-1].split(", ")]
            cases[cur][key] = (int(nr), int(nc), vals)
        except ValueError:
            # the LP blocks have other layouts; only the MILP blocks matter here
            continue
    return cases


def fj_input(d):
    import math
    c = d["c"][2]
    n = len(c)
    integ = d["integrality"][2]
    lb, ub = d["lb"][2], d["ub"][2]
    lines = [f"{SEED} {EPSILON.hex()} {FEASTOL.hex()}"]
    vars_, init = [], []
    for j in range(n):
        lo, up = lb[j], ub[j]
        isint = integ[j] != 0
        if isint:
            lo = math.ceil(lo - FEASTOL) if math.isfinite(lo) else lo
            up = math.floor(up + FEASTOL) if math.isfinite(up) else up
            lo, up = float(lo), float(up)
        vars_.append(f"{int(isint)} {lo.hex() if math.isfinite(lo) else lo} {up.hex() if math.isfinite(up) else up} "
                     f"{float(c[j]).hex()}")
        init.append(lo if math.isfinite(lo) else (up if math.isfinite(up) else 0.0))
    cons, nnz = [], 0
    if "A" in d:
        m, nc, a = d["A"]
        lbc, ubc = d["lbc"][2], d["ubc"][2]
        for i in range(m):
            idx = [j for j in range(nc) if a[i * nc + j] != 0]
            val = [a[i * nc + j] for j in idx]
            nnz += len(idx)
            body = f"{len(idx)} " + " ".join(map(str, idx)) + " " + " ".join(v.hex() for v in val)
            if math.isfinite(lbc[i]):
                cons.append(f"2 {lbc[i].hex()} {body}")
            if math.isfinite(ubc[i]):
                cons.append(f"1 {ubc[i].hex()} {body}")
    lines[0] += f" {nnz}"
    lines.append(str(n))
    lines += vars_
    lines.append(str(len(cons)))
    lines += cons
    lines.append(" ".join(float(v).hex() for v in init))
    return "\n".join(lines) + "\n"


def main():
    cases = parse(sys.argv[1])
    exe = os.path.join(tempfile.mkdtemp(), "fj")
    subprocess.check_call(["g++", "-O2", "-ffp-contract=off", "-std=c++17", "-I" + os.path.join(HERE, "fj"),
                           os.path.join(HERE, "fj-harness.cpp"), "-o", exe])
    for name, d in cases.items():
        if "integrality" not in d:
            continue
        inp = fj_input(d)
        print("== fj", name)
        for line in inp.splitlines():
            print("in|", line)
        out = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True).stdout
        sys.stdout.write(out)
        sys.stdout.flush()


main()
