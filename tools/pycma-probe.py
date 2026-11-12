# Deterministic pieces of pycma (the Python package cma, version 4.5.0)
# used by tests/testthat/test-optimization-cma.R: the strategy parameters,
# the recombination weights, the initial step size scaling with bounds and
# the bound transformation. Run with the cma package importable, for
# example `python3 -I tools/pycma-probe.py <directory containing cma>`.
import sys
if len(sys.argv) > 1:
    sys.path.insert(0, sys.argv[1])
import cma
from cma.utilities.math import Mh

for n, popsize, bounds in ((2, None, [[-1, -1], [0.8, 0.8]]), (10, None, [None, None]), (1, None, [-5, 5]),
                           (2, 20, [[-1, -1], [0.8, 0.8]])):
    opts = {'bounds': bounds, 'verbose': -9}
    if popsize:
        opts['popsize'] = popsize
    es = cma.CMAEvolutionStrategy(n * [0.0], 1, opts)
    sp = es.sp
    print("n=%d popsize=%d mu=%d lam_mirr=%d" % (n, sp.popsize, sp.weights.mu, sp.lam_mirr))
    print("  weights =", repr([float(w) for w in sp.weights]))
    print("  mueff cc c1 cmu =", repr([float(sp.weights.mueff), float(sp.cc), float(sp.c1), float(sp.cmu)]))
    print("  cs damps chiN =", repr([float(es.adapt_sigma.cs), float(es.adapt_sigma.damps), float(Mh.chiN(n))]))
    print("  maxiter tolstagnation =", repr([float(es.opts['maxiter']), float(es.opts['tolstagnation'])]))
    print("  scaling =", repr([float(s) for s in (es.sigma_vec.scaling * cma.utilities.utils.np.ones(n))]))
    print("  C diag =", repr([float(c) for c in es.sm.C.diagonal()]))

tf = cma.transformations.BoxConstraintsLinQuadTransformation([[-1, 0.8], [-1, 0.8]])
for x in ([0.0, 0.0], [-1.03, 0.79], [0.83, -1.2], [5.0, -7.3], [0.76, -0.97]):
    print("transform", repr(x), "=", repr([float(v) for v in tf(x)]))
tf2 = cma.transformations.BoxConstraintsLinQuadTransformation([[-10, None], [2, 3]])
print("transform2 =", repr([float(v) for v in tf2([-12.0, 3.4])]))
