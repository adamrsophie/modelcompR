/* Reference for the MINPACK port (R/optimization-scipy-minpack.R): runs
 * LMDIF of SciPy 1.18.1 (scipy/optimize/src/minpack.c) compiled without
 * floating point contraction on the system of the ScipyRoot("lm") test from
 * (1, 1, 1), where the Jacobian is singular, and prints every evaluated point
 * and the final result. SciPy's arm64 wheels contract to fused multiply adds
 * and take another path there. Build and run with
 *   clang -O2 -ffp-contract=off -I<scipy>/scipy/optimize/src \
 *     tools/minpack-nofma-check.c <scipy>/scipy/optimize/src/minpack.c -lm
 * MINPACK counts 25 evaluations; scipy.optimize.root adds the 2 calls of its
 * shape checks, so the port reports nfev 27. */
#include <stdio.h>
#include "minpack.h"

static int sys3(int *m, int *n, double *x, double *f, int *iflag) {
    if (*iflag == 1 || *iflag == 2) printf("%.17g %.17g %.17g\n", x[0], x[1], x[2]);
    f[0] = x[0] * x[0] + x[1] * x[1] - 4.0;
    f[1] = x[0] * x[1] - 1.0;
    f[2] = x[2] * x[2] * x[2] + x[0] - 2.0;
    return 0;
}

int main(void) {
    double x[3] = {1.0, 1.0, 1.0}, fvec[3], diag[3], fjac[9], qtf[3], wa1[3], wa2[3], wa3[3], wa4[3];
    int ipvt[3], info, nfev;
    LMDIF(sys3, 3, 3, x, fvec, 1.49012e-8, 1.49012e-8, 0.0, 800, 0.0, diag, 1, 100.0, 0, &info, &nfev, fjac, 3, ipvt,
          qtf, wa1, wa2, wa3, wa4);
    printf("info %d nfev %d x %.17g %.17g %.17g\n", info, nfev, x[0], x[1], x[2]);
    return 0;
}
