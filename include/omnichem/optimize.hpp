#pragma once
// Unconstrained minimisation: limited-memory BFGS with a backtracking
// (Armijo) line search and a per-coordinate step cap. Used by the UFF force
// field and the distance-geometry refinement.

#include <functional>
#include <vector>

namespace omnichem {

// Returns f(x) and writes the gradient into g (same size as x).
using Objective = std::function<double(const std::vector<double>& x, std::vector<double>& g)>;

struct LbfgsOptions {
  int max_iterations = 2000;
  int memory = 10;
  double grad_tol = 1e-4;   // converged when max |g_i| < grad_tol
  double f_rel_tol = 1e-14;  // or when the relative decrease of f stalls
  double max_step = 0.3;    // largest change of any coordinate per iteration
};

struct LbfgsResult {
  bool converged = false;
  int iterations = 0;
  int evaluations = 0;
  double f = 0.0;
  double grad_max = 0.0;
};

LbfgsResult minimize_lbfgs(const Objective& f, std::vector<double>& x, const LbfgsOptions& options = {});

}  // namespace omnichem
