// Independent MINPACK reference: Eigen's port (unsupported/NonLinearOptimization)
// of hybrj and lmder, with fdjac1/fdjac2 forward differences for the
// derivative free variants. Mirrors the MINPACK.jl calls used by
// NonlinearSolve's CMINPACK (hybrj: mode 2, diag 1, factor 100; hybrd1:
// maxfev 200(n+1); lmder/lmdif1: mode 1, factor 100, gtol 0, maxfev 200(n+1)).
// Build and run (Eigen 5 headers):
//   clang++ -std=c++17 -O0 -ffp-contract=off -I/opt/homebrew/include/eigen3 \
//     nonlinearsolve-cminpack-reference.cpp -o cminpack-reference && ./cminpack-reference
#include <cstdio>
#include <cmath>
#include <vector>
#include <functional>
#include <limits>
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
using namespace Eigen;
typedef std::function<VectorXd(const VectorXd&)> Fun;
typedef std::function<MatrixXd(const VectorXd&)> Jac;

static int cf, cg, cprev;
static bool tracing;
static void tpush(int iflag) {
  if (!tracing) return;
  if (iflag == 2) { if (cprev == 1) cg++; } else cf++;
  cprev = iflag;
}
struct Base {
  Fun f; Jac j; int m, n;
  std::vector<std::pair<VectorXd, VectorXd>> cache;
  VectorXd evalf(const VectorXd& x, int iflag) {
    VectorXd v = f(x); tpush(iflag);
    if (iflag == 1) cache.push_back({x, v});
    return v;
  }
  MatrixXd evalj(const VectorXd& x) {
    if (j) { tpush(2); return j(x); }
    // fdjac1 / fdjac2 (dense): h = sqrt(eps) |x_j|, or sqrt(eps) when x_j = 0
    VectorXd fx;
    for (int k = (int)cache.size() - 1; k >= 0; --k) if (cache[k].first == x) { fx = cache[k].second; break; }
    double eps = std::sqrt(std::numeric_limits<double>::epsilon());
    MatrixXd J(m, n); VectorXd y = x;
    for (int c = 0; c < n; ++c) {
      double t = y[c], h = eps * std::fabs(t);
      if (h == 0.) h = eps;
      y[c] = t + h;
      VectorXd w = evalf(y, 2);
      y[c] = t;
      J.col(c) = (w - fx) / h;
    }
    return J;
  }
};
struct HFun : Base {
  int operator()(const VectorXd& x, VectorXd& fv) { fv = evalf(x, 1); return 0; }
  int df(const VectorXd& x, MatrixXd& fj) { fj = evalj(x); return 0; }
};
struct LFun : Base {
  typedef double Scalar;
  enum { InputsAtCompileTime = Dynamic, ValuesAtCompileTime = Dynamic };
  typedef VectorXd InputType; typedef VectorXd ValueType; typedef MatrixXd JacobianType;
  int inputs() const { return n; }
  int values() const { return m; }
  int operator()(const VectorXd& x, VectorXd& fv) { fv = evalf(x, 1); return 0; }
  int df(const VectorXd& x, MatrixXd& fj) { fj = evalj(x); return 0; }
};
static void pv(const VectorXd& v) {
  printf("[");
  for (int i = 0; i < v.size(); ++i) printf("%s%.17g", i ? ", " : "", v[i]);
  printf("]");
}
static void run(const char* pn, const char* an, Fun f, Jac j, VectorXd x0, int m, bool lm, bool trace) {
  const double tol = std::pow(std::numeric_limits<double>::epsilon(), 0.8);
  cf = cg = cprev = 0; tracing = trace;
  int n = x0.size(); VectorXd x = x0, fv; int info;
  if (!lm) {
    HFun fun; fun.f = f; fun.j = j; fun.m = m; fun.n = n;
    HybridNonLinearSolver<HFun> s(fun);
    s.useExternalScaling = true; s.diag.setConstant(n, 1.);
    s.parameters.xtol = tol; s.parameters.factor = 100.;
    s.parameters.maxfev = j ? std::numeric_limits<int>::max() : 200 * (n + 1);
    info = s.solve(x);
    fv = s.fvec;
  } else {
    LFun fun; fun.f = f; fun.j = j; fun.m = m; fun.n = n;
    LevenbergMarquardt<LFun> s(fun);
    s.parameters.ftol = tol; s.parameters.xtol = tol; s.parameters.gtol = 0.;
    s.parameters.factor = 100.; s.parameters.maxfev = 200 * (n + 1);
    info = s.minimize(x);
    fv = s.fvec;
  }
  const char* rc = fv.cwiseAbs().maxCoeff() <= tol ? "Success" : "Failure";
  printf("RES %s %s %s ", pn, an, rc); pv(x); printf(" "); pv(fv);
  printf(" [%d, %d, %d, %d, -1] info %d\n", cf, cg, cg, cg, info);
}
int main() {
  Fun circf = [](const VectorXd& u) { VectorXd r(2); r << u[0]*u[0] + u[1]*u[1] - 4, u[0]*u[1] - 1; return r; };
  Jac circj = [](const VectorXd& u) { MatrixXd J(2, 2); J << 2*u[0], 2*u[1], u[1], u[0]; return J; };
  Fun rosenf = [](const VectorXd& u) { VectorXd r(2); r << 10*(u[1] - u[0]*u[0]), 1 - u[0]; return r; };
  Jac rosenj = [](const VectorXd& u) { MatrixXd J(2, 2); J << -20*u[0], 10, -1, 0; return J; };
  Fun fixpf = [](const VectorXd& u) { VectorXd r(2); r << 0.5*std::cos(u[1]) - u[0] + 0.1, 0.5*std::sin(u[0]) - u[1] + 0.2; return r; };
  Jac fixpj = [](const VectorXd& u) { MatrixXd J(2, 2); J << -1, -0.5*std::sin(u[1]), 0.5*std::cos(u[0]), -1; return J; };
  Fun burdenf = [](const VectorXd& u) { VectorXd r(3);
    r << 3*u[0] - std::cos(u[1]*u[2]) - 0.5, u[0]*u[0] - 81*(u[1] + 0.1)*(u[1] + 0.1) + std::sin(u[2]) + 1.06,
         std::exp(-u[0]*u[1]) + 20*u[2] + (10*M_PI - 3)/3; return r; };
  Jac burdenj = [](const VectorXd& u) { MatrixXd J(3, 3);
    J << 3, u[2]*std::sin(u[1]*u[2]), u[1]*std::sin(u[1]*u[2]),
         2*u[0], -162*(u[1] + 0.1), std::cos(u[2]),
         -u[1]*std::exp(-u[0]*u[1]), -u[0]*std::exp(-u[0]*u[1]), 20; return J; };
  double TT[6] = {0, 0.5, 1, 1.5, 2, 2.5}, YY[6] = {1, 1.6, 2.8, 4.4, 7.6, 12};
  Fun expf = [&](const VectorXd& u) { VectorXd r(6); for (int i = 0; i < 6; ++i) r[i] = u[0]*std::exp(u[1]*TT[i]) - YY[i]; return r; };
  Jac expj = [&](const VectorXd& u) { MatrixXd J(6, 2);
    for (int i = 0; i < 6; ++i) { J(i, 0) = std::exp(u[1]*TT[i]); J(i, 1) = u[0]*TT[i]*std::exp(u[1]*TT[i]); } return J; };
  struct P { const char* name; Fun f; Jac j; VectorXd u0; int m; };
  VectorXd a(2), b(2), c(2), d(3), e(2);
  a << 3, 0.1; b << -1.2, 1; c << 0, 0; d << 0.1, 0.1, -0.1; e << 1, 0.5;
  std::vector<P> ps = {{"circ", circf, circj, a, 2}, {"rosen", rosenf, rosenj, b, 2},
                       {"fixp", fixpf, fixpj, c, 2}, {"burden", burdenf, burdenj, d, 3}};
  for (auto& p : ps) for (int z = 0; z < 2; ++z) {
    char nm[32]; snprintf(nm, sizeof nm, "%s%s", p.name, z ? "0" : "");
    Jac jj = z ? Jac() : p.j;
    run(nm, "CMINPACK", p.f, jj, p.u0, p.m, false, false);
    run(nm, "CMINPACKTrace", p.f, jj, p.u0, p.m, false, true);
    run(nm, "CMINPACKlm", p.f, jj, p.u0, p.m, true, true);
  }
  run("expfit", "CMINPACK", expf, expj, e, 6, true, true);
  run("expfit0", "CMINPACK", expf, Jac(), e, 6, true, true);
  return 0;
}
