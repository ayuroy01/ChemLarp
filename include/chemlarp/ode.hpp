#pragma once
// Stiff ODE integration: variable-order (1-5), variable-step BDF.
//
// Formulation: at each step the corrector satisfies p'(t_{n+1}) = f(t_{n+1}, y_{n+1}),
// where p interpolates y_{n+1}, y_n, ..., y_{n+1-k} at their actual (unequal)
// times. Its coefficients are derivatives of the Lagrange basis, so variable
// steps need no Nordsieck rescaling. The corrector is solved by modified Newton
// with (c0 I - J) factored once per step attempt.
//
// Error: Milne's device, LTE ~ ||y_corr - y_pred|| / (k+2) in a weighted RMS
// norm (w_i = rtol |y_i| + atol_i). Order is raised or lowered by comparing the step
// sizes predicted for k-1, k, k+1. Output at requested times uses the same
// interpolating polynomial (dense output), so outputs do not shorten steps.

#include <functional>
#include <vector>

#include "chemlarp/diagnostics.hpp"

namespace chemlarp {

using OdeRhs = std::function<void(double t, const std::vector<double>& y, std::vector<double>& dydt)>;
// Row-major n x n, J[i*n + j] = d f_i / d y_j.
using OdeJacobian = std::function<void(double t, const std::vector<double>& y, std::vector<double>& J)>;
using OdeObserver = std::function<void(double t, const std::vector<double>& y)>;

struct BdfOptions {
  double rtol = 1e-6;
  std::vector<double> atol{1e-12};  // one value for all components, or one per component
  double h_initial = 0.0;           // 0 = automatic
  double h_max = 0.0;               // 0 = unlimited
  std::size_t max_steps = 500000;
  int max_order = 5;
  // Linear solver for the Newton iteration matrix: 0 = automatic (sparse LU when n >= 60 and at most a fifth
  // of the Jacobian is nonzero), 1 = dense LU, 2 = sparse LU (see sparse.hpp).
  int linear_solver = 0;
};

struct BdfStats {
  std::size_t steps = 0;
  std::size_t rejected = 0;
  std::size_t newton_failures = 0;
  std::size_t rhs_evals = 0;
  std::size_t jacobian_evals = 0;
  int max_order_used = 0;
  bool sparse = false;          // the last factorisation used the sparse solver
  std::size_t factor_nnz = 0;   // nonzeros of L + U in the last sparse factorisation
  double last_h = 0.0;
};

// Integrates from t0 to t_end. `y` holds y(t0) on entry and y(t_end) on success.
// `on_output` is called for every entry of `output_times` within (t0, t_end]
// (sorted ascending). `on_step` (optional) is called after every accepted step.
// If `jacobian` is null, a forward-difference Jacobian is used.
bool integrate_bdf(const OdeRhs& rhs, const OdeJacobian* jacobian, double t0, std::vector<double>& y, double t_end,
                   const std::vector<double>& output_times, const OdeObserver& on_output, const OdeObserver* on_step,
                   const BdfOptions& options, BdfStats& stats, Diagnostics& diag);

// LU factorisation with partial pivoting, reusable for several right-hand sides.
class LuFactor {
 public:
  bool factor(std::vector<double> a, std::size_t n);  // row-major copy
  void solve(std::vector<double>& b) const;           // in place
 private:
  std::size_t n_ = 0;
  std::vector<double> a_;
  std::vector<std::size_t> piv_;
};

}  // namespace chemlarp
