// Step level trace of the feasibility jump heuristic that HiGHS 1.12.0
// vendors (mip/feasibilityjump.hh, copied with its stub includes under
// tools/fj), driven as HighsMipSolverData::feasibilityJump drives it.
// Build: g++ -O2 -ffp-contract=off -std=c++17 -Itools/fj tools/fj-harness.cpp
//
// Input on stdin (whitespace separated, doubles in any strtod format,
// hex floats included):
//   seed equalityTolerance violationTolerance nnz
//   nvars, then per variable: type(0 continuous, 1 integer) lb ub obj
//   ncons, then per constraint: sense(1 Lte, 2 Gte) rhs k idx_1..idx_k val_1..val_k
//   nvars initial values
// Output: the first 24 outputs of std::mt19937(seed) and of
// std::uniform_real_distribution<double>(0, 1) over a fresh mt19937(seed),
// then one line per step of the solve loop (replicated below with the
// private members made visible), the solutions the callback receives, and
// the result of the unmodified solve() for comparison.
#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <random>
#include <string>
#include <vector>
// make every member visible (class members before the first access
// specifier are private, hence class -> struct; enum class -> enum struct)
#include "util/HighsRandom.h"
#define private public
#define class struct
#include "feasibilityjump.hh"
#undef class
#undef private

using namespace external_feasibilityjump;

static double rd() {
  char buf[128];
  if (scanf("%127s", buf) != 1) exit(2);
  return strtod(buf, nullptr);
}

struct Input {
  int seed;
  double eqtol, violtol;
  long nnz;
  std::vector<int> type;
  std::vector<double> lb, ub, obj, init;
  std::vector<int> sense;
  std::vector<double> rhs;
  std::vector<std::vector<int>> idx;
  std::vector<std::vector<double>> val;
};

static void build(const Input& in, FeasibilityJumpSolver& s) {
  for (size_t j = 0; j < in.lb.size(); ++j)
    s.addVar(in.type[j] ? VarType::Integer : VarType::Continuous, in.lb[j], in.ub[j], in.obj[j]);
  for (size_t i = 0; i < in.rhs.size(); ++i) {
    std::vector<int> ix = in.idx[i];
    std::vector<double> vx = in.val[i];
    s.addConstraint(in.sense[i] == 1 ? RowType::Lte : RowType::Gte, in.rhs[i], (int)ix.size(), ix.data(),
                    vx.data(), 0);
  }
}

int main() {
  Input in;
  in.seed = (int)rd();
  in.eqtol = rd();
  in.violtol = rd();
  in.nnz = (long)rd();
  int n = (int)rd();
  for (int j = 0; j < n; ++j) {
    in.type.push_back((int)rd());
    in.lb.push_back(rd());
    in.ub.push_back(rd());
    in.obj.push_back(rd());
  }
  int m = (int)rd();
  for (int i = 0; i < m; ++i) {
    in.sense.push_back((int)rd());
    in.rhs.push_back(rd());
    int k = (int)rd();
    std::vector<int> ix(k);
    std::vector<double> vx(k);
    for (int t = 0; t < k; ++t) ix[t] = (int)rd();
    for (int t = 0; t < k; ++t) vx[t] = rd();
    in.idx.push_back(ix);
    in.val.push_back(vx);
  }
  for (int j = 0; j < n; ++j) in.init.push_back(rd());

  {
    std::mt19937 r(in.seed);
    printf("mt");
    for (int t = 0; t < 24; ++t) printf(" %lu", (unsigned long)r());
    printf("\n");
    std::mt19937 r2(in.seed);
    printf("unif");
    for (int t = 0; t < 24; ++t) printf(" %a", std::uniform_real_distribution<double>(0., 1.)(r2));
    printf("\n");
    // HighsRandom (util/HighsRandom.h) for the primal heuristics
    HighsRandom h1(in.seed), h2(in.seed), h3(in.seed), h4(in.seed);
    printf("hrand_closed");
    for (int t = 0; t < 8; ++t) printf(" %a", h1.closedFraction());
    printf("\nhrand_real");
    for (int t = 0; t < 8; ++t) printf(" %a", h2.real(0.1, 0.9));
    printf("\nhrand_bit");
    for (int t = 0; t < 8; ++t) printf(" %d", (int)h3.bit());
    printf("\nhrand_mixed");
    for (int t = 0; t < 4; ++t) printf(" %a %d", h4.real(0.6, 0.6), (int)h4.integer(7));
    printf("\n");
  }

  const size_t kMaxTotalEffort = (size_t)in.nnz << 10;
  const size_t kMaxEffortSinceLastImprovement = (size_t)in.nnz << 8;
  std::vector<double> col_value;
  double solobj = NAN;
  bool found = false;
  auto cb = [&](FJStatus status) -> CallbackControlFlow {
    printf("cb effort %zu since %zu sol %d", status.totalEffort, status.effortSinceLastImprovement,
           status.solution != nullptr);
    if (status.solution != nullptr) {
      found = true;
      col_value = std::vector<double>(status.solution, status.solution + status.numVars);
      solobj = status.solutionObjectiveValue;
      printf(" obj %a", solobj);
    }
    printf("\n");
    if (status.effortSinceLastImprovement > kMaxEffortSinceLastImprovement ||
        status.totalEffort > kMaxTotalEffort)
      return CallbackControlFlow::Terminate;
    return CallbackControlFlow::Continue;
  };

  HighsLogOptions lo;
  {
    // the loop of FeasibilityJumpSolver::solve without the logging, traced
    FeasibilityJumpSolver s(lo, in.seed, in.eqtol, in.violtol);
    build(in, s);
    std::vector<double> init = in.init;
    std::function<CallbackControlFlow(FJStatus)> f = cb;
    s.init(init.data());
    printf("init effort %zu viol %zu good %zu obj %a\n", s.totalEffort, s.problem.violatedConstraints.size(),
           s.goodVarsSet.size(), s.problem.incumbentObjective);
    for (int step = 0; step < INT_MAX; step += 1) {
      if (s.user_terminate(f, nullptr)) break;
      if (s.problem.violatedConstraints.size() < s.bestViolationScore) {
        s.effortAtLastImprovement = s.totalEffort;
        s.bestViolationScore = s.problem.violatedConstraints.size();
      }
      if (s.problem.violatedConstraints.empty() && s.problem.incumbentObjective < s.bestObjective) {
        double rel = (s.bestObjective - s.problem.incumbentObjective) /
                     std::max(1.0, std::fabs(s.problem.incumbentObjective));
        if (rel > s.kMinRelativeObjectiveImprovement) {
          s.effortAtLastImprovement = s.totalEffort;
          s.bestObjective = s.problem.incumbentObjective;
          if (s.user_terminate(f, s.problem.incumbentAssignment.data())) break;
        }
      }
      if (s.problem.vars.size() == 0) break;
      uint32_t var = s.selectVariable();
      s.doVariableMove(var);
      if (step < 3000 || step % 500 == 0)
        printf("s %d %u %a %zu %zu %zu %zu %a\n", step, var, s.problem.incumbentAssignment[var], s.totalEffort,
               s.problem.violatedConstraints.size(), s.goodVarsSet.size(), s.nBumps,
               s.problem.incumbentObjective);
    }
    printf("end effort %zu bumps %zu found %d", s.totalEffort, s.nBumps, (int)found);
    if (found) {
      printf(" obj %a x", solobj);
      for (double v : col_value) printf(" %a", v);
    }
    printf("\n");
  }
  {
    // the unmodified solve() for comparison with the traced loop above
    found = false;
    FeasibilityJumpSolver s(lo, in.seed, in.eqtol, in.violtol);
    build(in, s);
    std::vector<double> init = in.init;
    printf("solve\n");
    s.solve(init.data(), cb);
    printf("solve end effort %zu found %d", s.totalEffort, (int)found);
    if (found) {
      printf(" obj %a x", solobj);
      for (double v : col_value) printf(" %a", v);
    }
    printf("\n");
  }
  return 0;
}
