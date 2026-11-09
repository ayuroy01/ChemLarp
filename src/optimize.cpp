#include "omnichem/optimize.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <vector>

namespace omnichem {

namespace {
double dot(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}
double max_abs(const std::vector<double>& a) {
  double m = 0.0;
  for (double v : a) m = std::max(m, std::fabs(v));
  return m;
}
}  // namespace

LbfgsResult minimize_lbfgs(const Objective& func, std::vector<double>& x, const LbfgsOptions& opt) {
  LbfgsResult res;
  const std::size_t n = x.size();
  std::vector<double> g(n), d(n), x_new(n), g_new(n);
  double f = func(x, g);
  ++res.evaluations;
  std::deque<std::vector<double>> S, Y;
  std::deque<double> rho;

  for (res.iterations = 0; res.iterations < opt.max_iterations; ++res.iterations) {
    res.grad_max = max_abs(g);
    if (res.grad_max < opt.grad_tol || !std::isfinite(f)) {
      res.converged = std::isfinite(f);
      break;
    }
    // Two-loop recursion: d = -H g.
    d = g;
    std::vector<double> alpha(S.size());
    for (std::size_t k = S.size(); k-- > 0;) {
      alpha[k] = rho[k] * dot(S[k], d);
      for (std::size_t i = 0; i < n; ++i) d[i] -= alpha[k] * Y[k][i];
    }
    if (!S.empty()) {
      const double gamma = dot(S.back(), Y.back()) / dot(Y.back(), Y.back());
      for (double& v : d) v *= gamma;
    }
    for (std::size_t k = 0; k < S.size(); ++k) {
      const double beta = rho[k] * dot(Y[k], d);
      for (std::size_t i = 0; i < n; ++i) d[i] += (alpha[k] - beta) * S[k][i];
    }
    for (double& v : d) v = -v;
    double slope = dot(d, g);
    if (!(slope < 0.0)) {  // not a descent direction: restart with steepest descent
      S.clear(); Y.clear(); rho.clear();
      for (std::size_t i = 0; i < n; ++i) d[i] = -g[i];
      slope = dot(d, g);
    }
    // Cap the largest coordinate change.
    const double dmax = max_abs(d);
    if (dmax > opt.max_step) {
      const double s = opt.max_step / dmax;
      for (double& v : d) v *= s;
      slope *= s;
    }
    // Backtracking Armijo line search.
    double step = 1.0, f_new = f;
    bool accepted = false;
    for (int bt = 0; bt < 40; ++bt) {
      for (std::size_t i = 0; i < n; ++i) x_new[i] = x[i] + step * d[i];
      f_new = func(x_new, g_new);
      ++res.evaluations;
      if (std::isfinite(f_new) && f_new <= f + 1e-4 * step * slope) { accepted = true; break; }
      step *= 0.5;
    }
    if (!accepted) {
      if (S.empty()) break;  // even steepest descent cannot decrease f: at a (numerical) minimum
      S.clear(); Y.clear(); rho.clear();
      continue;
    }
    std::vector<double> s(n), yv(n);
    for (std::size_t i = 0; i < n; ++i) {
      s[i] = x_new[i] - x[i];
      yv[i] = g_new[i] - g[i];
    }
    const double sy = dot(s, yv);
    if (sy > 1e-12 * std::sqrt(dot(s, s) * dot(yv, yv))) {
      S.push_back(std::move(s));
      Y.push_back(std::move(yv));
      rho.push_back(1.0 / sy);
      if (static_cast<int>(S.size()) > opt.memory) { S.pop_front(); Y.pop_front(); rho.pop_front(); }
    }
    const double df = f - f_new;
    x.swap(x_new);
    g.swap(g_new);
    f = f_new;
    if (df >= 0 && df <= opt.f_rel_tol * std::max(1.0, std::fabs(f)) && max_abs(g) < 100 * opt.grad_tol) {
      res.converged = true;
      ++res.iterations;
      break;
    }
  }
  res.f = f;
  res.grad_max = max_abs(g);
  if (res.grad_max < opt.grad_tol) res.converged = true;
  return res;
}

}  // namespace omnichem
