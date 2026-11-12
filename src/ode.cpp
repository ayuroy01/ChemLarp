#include "chemlarp/ode.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <deque>
#include <string>

namespace chemlarp {

bool LuFactor::factor(std::vector<double> a, std::size_t n) {
  n_ = n;
  a_ = std::move(a);
  piv_.resize(n);
  double scale = 0.0;
  for (double v : a_) scale = std::max(scale, std::fabs(v));
  if (scale == 0.0 || !std::isfinite(scale)) return false;
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t p = k;
    for (std::size_t r = k + 1; r < n; ++r)
      if (std::fabs(a_[r * n + k]) > std::fabs(a_[p * n + k])) p = r;
    piv_[k] = p;
    if (std::fabs(a_[p * n + k]) < 1e-300 + 1e-15 * scale) return false;
    if (p != k)
      for (std::size_t c = 0; c < n; ++c) std::swap(a_[k * n + c], a_[p * n + c]);
    const double inv = 1.0 / a_[k * n + k];
    for (std::size_t r = k + 1; r < n; ++r) {
      double& l = a_[r * n + k];
      if (l == 0.0) continue;
      l *= inv;
      for (std::size_t c = k + 1; c < n; ++c) a_[r * n + c] -= l * a_[k * n + c];
    }
  }
  return true;
}

void LuFactor::solve(std::vector<double>& b) const {
  const std::size_t n = n_;
  // factor() swaps whole rows (multipliers included), so L is stored in the final
  // row order: apply every interchange to b first, then forward-substitute.
  for (std::size_t k = 0; k < n; ++k) std::swap(b[k], b[piv_[k]]);
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t r = k + 1; r < n; ++r) b[r] -= a_[r * n + k] * b[k];
  for (std::size_t i = n; i-- > 0;) {
    double s = b[i];
    for (std::size_t c = i + 1; c < n; ++c) s -= a_[i * n + c] * b[c];
    b[i] = s / a_[i * n + i];
  }
}

namespace {

constexpr std::size_t kHistory = 7;

double wrms(const std::vector<double>& v, const std::vector<double>& w) {
  double s = 0.0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    const double r = v[i] / w[i];
    s += r * r;
  }
  return v.empty() ? 0.0 : std::sqrt(s / static_cast<double>(v.size()));
}

// d/dt of the Lagrange basis polynomials at nodes[0].
void derivative_coefficients(const std::vector<double>& nodes, std::vector<double>& c) {
  const std::size_t m = nodes.size();
  c.assign(m, 0.0);
  for (std::size_t i = 1; i < m; ++i) c[0] += 1.0 / (nodes[0] - nodes[i]);
  for (std::size_t j = 1; j < m; ++j) {
    double num = 1.0, den = 1.0;
    for (std::size_t i = 0; i < m; ++i) {
      if (i == j) continue;
      den *= nodes[j] - nodes[i];
      if (i != 0) num *= nodes[0] - nodes[i];
    }
    c[j] = num / den;
  }
}

struct Point {
  double t;
  std::vector<double> y;
};

// Value at t of the polynomial through `pts` (all used).
void interpolate(const std::vector<const Point*>& pts, double t, std::vector<double>& out) {
  const std::size_t n = pts.front()->y.size();
  out.assign(n, 0.0);
  for (std::size_t j = 0; j < pts.size(); ++j) {
    double L = 1.0;
    for (std::size_t i = 0; i < pts.size(); ++i)
      if (i != j) L *= (t - pts[i]->t) / (pts[j]->t - pts[i]->t);
    for (std::size_t c = 0; c < n; ++c) out[c] += L * pts[j]->y[c];
  }
}

}  // namespace

bool integrate_bdf(const OdeRhs& rhs, const OdeJacobian* jacobian, double t0, std::vector<double>& y, double t_end,
                   const std::vector<double>& output_times, const OdeObserver& on_output, const OdeObserver* on_step,
                   const BdfOptions& opt, BdfStats& stats, Diagnostics& diag) {
  const std::size_t n = y.size();
  if (n == 0 || !(t_end > t0) || !(opt.rtol > 0.0) || opt.atol.empty() || (opt.atol.size() != 1 && opt.atol.size() != n)) {
    diag.error(DiagCode::InvalidReactorSetup, "integrator needs a non-empty state, t_end > t0, rtol > 0 and 1 or n atol values");
    return false;
  }
  std::vector<double> atol(n);
  for (std::size_t i = 0; i < n; ++i) atol[i] = opt.atol.size() == 1 ? opt.atol[0] : opt.atol[i];
  const int max_order = std::clamp(opt.max_order, 1, 5);

  auto f = [&](double t, const std::vector<double>& yy, std::vector<double>& out) {
    ++stats.rhs_evals;
    out.resize(n);
    rhs(t, yy, out);
  };
  auto weights = [&](const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& w) {
    w.resize(n);
    for (std::size_t i = 0; i < n; ++i) w[i] = opt.rtol * std::max(std::fabs(a[i]), std::fabs(b[i])) + atol[i];
  };

  std::deque<Point> hist;
  hist.push_back({t0, y});
  std::vector<double> f0, w, tmp;
  f(t0, y, f0);
  for (double v : f0)
    if (!std::isfinite(v)) {
      diag.error(DiagCode::InvalidReactorSetup, "right-hand side is not finite at the initial state");
      return false;
    }

  double h = opt.h_initial;
  if (!(h > 0.0)) {
    weights(y, y, w);
    const double d0 = wrms(y, w), d1 = wrms(f0, w);
    h = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 * (t_end - t0) : 0.01 * d0 / d1;
  }
  h = std::min(h, t_end - t0);
  if (opt.h_max > 0.0) h = std::min(h, opt.h_max);

  std::size_t next_out = 0;
  while (next_out < output_times.size() && output_times[next_out] <= t0) {
    if (output_times[next_out] == t0) on_output(t0, y);
    ++next_out;
  }

  int k = 1, steps_at_order = 0, consecutive_fail = 0;
  double t = t0;
  std::vector<double> c, y_pred, y_new, fv, G, J(n * n), M(n * n), fp;
  std::vector<double> nodes;
  LuFactor lu;

  while (t < t_end) {
    if (stats.steps >= opt.max_steps) {
      diag.error(DiagCode::TooManySteps, "integration stopped after " + std::to_string(stats.steps) +
                                             " steps at t = " + std::to_string(t));
      return false;
    }
    if (t + h >= t_end || t_end - (t + h) < 1e-12 * std::fabs(t_end)) h = t_end - t;
    if (h <= 4.0 * DBL_EPSILON * std::max(std::fabs(t), 1.0)) {
      diag.error(DiagCode::StepSizeUnderflow, "step size underflow at t = " + std::to_string(t));
      return false;
    }
    const double t_new = t + h;
    const int ke = std::min<int>(k, static_cast<int>(hist.size()));

    nodes.assign(1, t_new);
    for (int j = 0; j < ke; ++j) nodes.push_back(hist[hist.size() - 1 - static_cast<std::size_t>(j)].t);
    derivative_coefficients(nodes, c);

    // Predictor: extrapolate through up to ke+1 past points (explicit Euler on the very first step).
    const Point& last = hist.back();
    if (hist.size() == 1) {
      y_pred = last.y;
      for (std::size_t i = 0; i < n; ++i) y_pred[i] += h * f0[i];
    } else {
      std::vector<const Point*> pts;
      const std::size_t use = std::min<std::size_t>(static_cast<std::size_t>(ke) + 1, hist.size());
      for (std::size_t j = 0; j < use; ++j) pts.push_back(&hist[hist.size() - 1 - j]);
      interpolate(pts, t_new, y_pred);
    }

    // Jacobian at the prediction; iteration matrix c0 I - J.
    ++stats.jacobian_evals;
    if (jacobian) {
      (*jacobian)(t_new, y_pred, J);
    } else {
      f(t_new, y_pred, fp);
      double ymax = 0.0;
      for (double v : y_pred) ymax = std::max(ymax, std::fabs(v));
      for (std::size_t j = 0; j < n; ++j) {
        tmp = y_pred;
        // Scale-aware perturbation: tiny components get a step tied to the state's
        // magnitude, otherwise rounding noise in f is amplified by 1/d.
        const double d = std::sqrt(DBL_EPSILON) * std::max({std::fabs(y_pred[j]), 1e-6 * ymax, 1e3 * atol[j]});
        tmp[j] += d;
        f(t_new, tmp, fv);
        for (std::size_t i = 0; i < n; ++i) J[i * n + j] = (fv[i] - fp[i]) / d;
      }
    }
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j) M[i * n + j] = (i == j ? c[0] : 0.0) - J[i * n + j];
    bool converged = lu.factor(M, n);

    std::vector<double> s(n, 0.0);  // history part of the corrector
    for (int j = 1; j <= ke; ++j) {
      const auto& yj = hist[hist.size() - static_cast<std::size_t>(j)].y;
      for (std::size_t i = 0; i < n; ++i) s[i] += c[static_cast<std::size_t>(j)] * yj[i];
    }
    y_new = y_pred;
    if (converged) {
      converged = false;
      double prev = 0.0;
      for (int it = 0; it < 7; ++it) {
        f(t_new, y_new, fv);
        G.resize(n);
        for (std::size_t i = 0; i < n; ++i) G[i] = -(c[0] * y_new[i] + s[i] - fv[i]);
        lu.solve(G);
        bool finite = true;
        for (std::size_t i = 0; i < n; ++i) {
          y_new[i] += G[i];
          finite = finite && std::isfinite(y_new[i]);
        }
        if (!finite) break;
        weights(y_new, last.y, w);
        const double nrm = wrms(G, w);
        if (nrm < 1e-2) { converged = true; break; }
        if (it > 0 && nrm > 2.0 * prev) break;
        prev = nrm;
      }
    }
    if (!converged) {
      ++stats.newton_failures;
      h *= 0.25;
      if (++consecutive_fail > 3) { k = 1; steps_at_order = 0; }
      continue;
    }

    weights(y_new, last.y, w);
    tmp.resize(n);
    for (std::size_t i = 0; i < n; ++i) tmp[i] = y_new[i] - y_pred[i];
    const double err = wrms(tmp, w) / (hist.size() == 1 ? 2.0 : static_cast<double>(ke + 2));
    if (!(err <= 1.0)) {
      ++stats.rejected;
      const double factor = std::isfinite(err) ? std::max(0.2, 0.9 * std::pow(err, -1.0 / (ke + 1))) : 0.2;
      h *= factor;
      if (++consecutive_fail >= 3) { k = std::max(1, k - 1); steps_at_order = 0; }
      continue;
    }

    // Accepted. Error estimates at neighbouring orders, from past points only.
    auto order_error = [&](int q) {
      std::vector<const Point*> pts;
      for (int j = 0; j <= q; ++j) pts.push_back(&hist[hist.size() - 1 - static_cast<std::size_t>(j)]);
      std::vector<double> p;
      interpolate(pts, t_new, p);
      for (std::size_t i = 0; i < n; ++i) p[i] = y_new[i] - p[i];
      return wrms(p, w) / static_cast<double>(q + 2);
    };
    auto step_for = [&](double e, int q) { return h * std::min(5.0, 0.9 * std::pow(std::max(e, 1e-10), -1.0 / (q + 1))); };
    double h_next = step_for(err, ke);
    int k_next = ke;
    if (steps_at_order >= ke + 1 && hist.size() >= 2) {
      if (ke > 1) {
        const double hl = step_for(order_error(ke - 1), ke - 1);
        if (hl / 1.3 > h_next) { h_next = hl; k_next = ke - 1; }
      }
      if (ke < max_order && hist.size() >= static_cast<std::size_t>(ke) + 2) {
        const double hh = step_for(order_error(ke + 1), ke + 1);
        if (hh / 1.4 > h_next) { h_next = hh; k_next = ke + 1; }
      }
    }

    // Dense output through the corrector polynomial.
    if (next_out < output_times.size() && output_times[next_out] <= t_new) {
      Point np{t_new, y_new};
      std::vector<const Point*> pts{&np};
      for (int j = 0; j < ke; ++j) pts.push_back(&hist[hist.size() - 1 - static_cast<std::size_t>(j)]);
      std::vector<double> yo;
      while (next_out < output_times.size() && output_times[next_out] <= t_new) {
        if (output_times[next_out] == t_new) yo = y_new;
        else interpolate(pts, output_times[next_out], yo);
        on_output(output_times[next_out], yo);
        ++next_out;
      }
    }

    hist.push_back({t_new, y_new});
    if (hist.size() > kHistory) hist.pop_front();
    t = t_new;
    ++stats.steps;
    stats.last_h = h;
    stats.max_order_used = std::max(stats.max_order_used, ke);
    consecutive_fail = 0;
    if (on_step) (*on_step)(t, y_new);
    steps_at_order = (k_next == ke) ? steps_at_order + 1 : 0;
    k = k_next;
    h = std::max(h_next, 0.2 * h);
    if (opt.h_max > 0.0) h = std::min(h, opt.h_max);
  }
  y = hist.back().y;
  return true;
}

}  // namespace chemlarp
