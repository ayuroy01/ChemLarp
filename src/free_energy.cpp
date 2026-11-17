// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
#include "chemlarp/free_energy.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <random>
#include <sstream>

#include "chemlarp/json.hpp"
#include "chemlarp/parallel.hpp"

namespace chemlarp {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

double logsumexp(const double* a, std::size_t n) {
  double m = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < n; ++i) m = std::max(m, a[i]);
  if (!std::isfinite(m)) return m;
  double s = 0.0;
  for (std::size_t i = 0; i < n; ++i) s += std::exp(a[i] - m);
  return m + std::log(s);
}

double mean_of(const std::vector<double>& a) {
  double s = 0.0;
  for (double x : a) s += x;
  return a.empty() ? 0.0 : s / static_cast<double>(a.size());
}

// Cyclic Jacobi eigen-decomposition of a symmetric n × n matrix (row-major); eigenvectors in the columns of V.
void jacobi(std::vector<double> a, std::size_t n, std::vector<double>& w, std::vector<double>& V) {
  V.assign(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) V[i * n + i] = 1.0;
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = i + 1; j < n; ++j) off += a[i * n + j] * a[i * n + j];
    if (off < 1e-30) break;
    for (std::size_t p = 0; p < n; ++p)
      for (std::size_t q = p + 1; q < n; ++q) {
        const double apq = a[p * n + q];
        if (std::fabs(apq) < 1e-300) continue;
        const double theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
        const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
        for (std::size_t k = 0; k < n; ++k) {
          const double akp = a[k * n + p], akq = a[k * n + q];
          a[k * n + p] = c * akp - s * akq;
          a[k * n + q] = s * akp + c * akq;
        }
        for (std::size_t k = 0; k < n; ++k) {
          const double apk = a[p * n + k], aqk = a[q * n + k];
          a[p * n + k] = c * apk - s * aqk;
          a[q * n + k] = s * apk + c * aqk;
        }
        for (std::size_t k = 0; k < n; ++k) {
          const double vkp = V[k * n + p], vkq = V[k * n + q];
          V[k * n + p] = c * vkp - s * vkq;
          V[k * n + q] = s * vkp + c * vkq;
        }
      }
  }
  w.resize(n);
  for (std::size_t i = 0; i < n; ++i) w[i] = a[i * n + i];
}

// Solves A x = b in place (n × n, row-major) by Gaussian elimination with partial pivoting.
bool gauss_solve(std::vector<double> A, std::vector<double>& b, std::size_t n) {
  for (std::size_t c = 0; c < n; ++c) {
    std::size_t piv = c;
    for (std::size_t r = c + 1; r < n; ++r)
      if (std::fabs(A[r * n + c]) > std::fabs(A[piv * n + c])) piv = r;
    if (std::fabs(A[piv * n + c]) < 1e-300) return false;
    if (piv != c) {
      for (std::size_t k = 0; k < n; ++k) std::swap(A[c * n + k], A[piv * n + k]);
      std::swap(b[c], b[piv]);
    }
    for (std::size_t r = c + 1; r < n; ++r) {
      const double m = A[r * n + c] / A[c * n + c];
      if (m == 0.0) continue;
      for (std::size_t k = c; k < n; ++k) A[r * n + k] -= m * A[c * n + k];
      b[r] -= m * b[c];
    }
  }
  for (std::size_t c = n; c-- > 0;) {
    double s = b[c];
    for (std::size_t k = c + 1; k < n; ++k) s -= A[c * n + k] * b[k];
    b[c] = s / A[c * n + c];
  }
  return true;
}

// Reproducible deviates (std::normal_distribution differs between standard libraries).
class Rng {
 public:
  explicit Rng(std::uint64_t seed) : g_(seed) {}
  double uniform() { return (static_cast<double>(g_() >> 11) + 0.5) * 0x1.0p-53; }
  double normal() {
    if (have_) {
      have_ = false;
      return spare_;
    }
    const double u = uniform(), v = uniform(), r = std::sqrt(-2.0 * std::log(u));
    spare_ = r * std::sin(2.0 * std::numbers::pi * v);
    have_ = true;
    return r * std::cos(2.0 * std::numbers::pi * v);
  }

 private:
  std::mt19937_64 g_;
  bool have_ = false;
  double spare_ = 0.0;
};

}  // namespace

// ---- estimators -----------------------------------------------------------------------------------------

FreeEnergyEstimate exp_average(const std::vector<double>& w) {
  FreeEnergyEstimate r;
  const std::size_t n = w.size();
  if (n == 0) return {kNaN, kNaN};
  std::vector<double> neg(n);
  for (std::size_t i = 0; i < n; ++i) neg[i] = -w[i];
  const double lse = logsumexp(neg.data(), n);
  r.df = -(lse - std::log(static_cast<double>(n)));
  // Delta method: var(ΔF) = var(e^{−w}) / (n ⟨e^{−w}⟩²), computed relative to the mean to avoid overflow.
  double s2 = 0.0;
  const double logmean = lse - std::log(static_cast<double>(n));
  for (double x : neg) {
    const double e = std::exp(x - logmean);
    s2 += (e - 1.0) * (e - 1.0);
  }
  r.error = n > 1 ? std::sqrt(s2 / static_cast<double>(n - 1) / static_cast<double>(n)) : kNaN;
  return r;
}

FreeEnergyEstimate bar(const std::vector<double>& wf, const std::vector<double>& wr) {
  const double nf = static_cast<double>(wf.size()), nr = static_cast<double>(wr.size());
  if (wf.empty() || wr.empty()) return {kNaN, kNaN};
  const double M = std::log(nf / nr);
  auto fermi = [](double x) { return x > 0 ? std::exp(-x) / (1.0 + std::exp(-x)) : 1.0 / (1.0 + std::exp(x)); };
  // g(ΔF) = Σ_F f(M + w − ΔF) − Σ_R f(−M + w + ΔF) increases with ΔF; its root is the BAR estimate.
  auto g = [&](double df) {
    double a = 0.0, b = 0.0;
    for (double w : wf) a += fermi(M + w - df);
    for (double w : wr) b += fermi(-M + w + df);
    return a - b;
  };
  double lo = exp_average(wf).df, hi = -exp_average(wr).df;
  if (!std::isfinite(lo)) lo = -1e3;
  if (!std::isfinite(hi)) hi = 1e3;
  if (lo > hi) std::swap(lo, hi);
  lo -= 1.0, hi += 1.0;
  for (int k = 0; k < 60 && g(lo) > 0; ++k) lo -= (hi - lo);
  for (int k = 0; k < 60 && g(hi) < 0; ++k) hi += (hi - lo);
  for (int it = 0; it < 200 && hi - lo > 1e-12 * std::max(1.0, std::fabs(lo)); ++it) {
    const double mid = 0.5 * (lo + hi);
    (g(mid) < 0 ? lo : hi) = mid;
  }
  FreeEnergyEstimate r;
  r.df = 0.5 * (lo + hi);
  // Asymptotic variance (Shirts et al. 2003).
  double fF = 0.0, fF2 = 0.0, fR = 0.0, fR2 = 0.0;
  for (double w : wf) {
    const double x = fermi(M + w - r.df);
    fF += x, fF2 += x * x;
  }
  for (double w : wr) {
    const double x = fermi(-M + w + r.df);
    fR += x, fR2 += x * x;
  }
  fF /= nf, fF2 /= nf, fR /= nr, fR2 /= nr;
  const double var = (fF2 / (fF * fF) - 1.0) / nf + (fR2 / (fR * fR) - 1.0) / nr;
  r.error = var > 0 ? std::sqrt(var) : 0.0;
  return r;
}

MbarResult mbar(const std::vector<std::vector<double>>& u, const std::vector<std::size_t>& nk, double tol, int max_it) {
  MbarResult r;
  const std::size_t K = u.size();
  if (K == 0 || nk.size() != K) return r;
  const std::size_t N = u[0].size();
  std::size_t total = 0;
  for (std::size_t k = 0; k < K; ++k) {
    if (u[k].size() != N) return r;
    total += nk[k];
  }
  if (total != N || N == 0) return r;
  std::vector<double> lognk(K);
  for (std::size_t k = 0; k < K; ++k) lognk[k] = nk[k] ? std::log(static_cast<double>(nk[k])) : -std::numeric_limits<double>::infinity();
  std::vector<double> f(K, 0.0), logden(N), tmp(std::max(K, N));
  auto denominators = [&] {
    for (std::size_t n = 0; n < N; ++n) {
      for (std::size_t k = 0; k < K; ++k) tmp[k] = lognk[k] + f[k] - u[k][n];
      logden[n] = logsumexp(tmp.data(), K);
    }
  };
  // f_i = −ln Σ_n exp(−u_in − logden_n)
  auto update = [&](std::vector<double>& out) {
    for (std::size_t i = 0; i < K; ++i) {
      for (std::size_t n = 0; n < N; ++n) tmp[n] = -u[i][n] - logden[n];
      out[i] = -logsumexp(tmp.data(), N);
    }
    const double f0 = out[0];
    for (double& x : out) x -= f0;
  };
  std::vector<double> fn(K);
  // Start from exponential averaging on the samples of state 0 (when it has any); the iteration corrects
  // a poor start.
  if (nk[0]) {
    for (std::size_t i = 1; i < K; ++i) {
      std::vector<double> w(nk[0]);
      for (std::size_t n = 0; n < nk[0]; ++n) w[n] = u[i][n] - u[0][n];
      const double e = exp_average(w).df;
      if (std::isfinite(e)) f[i] = e;
    }
  }
  int it = 0;
  bool done = false;
  for (; it < max_it; ++it) {
    denominators();
    update(fn);
    double dmax = 0.0;
    for (std::size_t i = 0; i < K; ++i) dmax = std::max(dmax, std::fabs(fn[i] - f[i]));
    f = fn;
    if (dmax < 1e-6) break;
  }
  // Newton on the convex objective over the sampled states (first sampled state fixed).
  std::vector<std::size_t> idx;
  for (std::size_t k = 0; k < K; ++k)
    if (nk[k]) idx.push_back(k);
  std::vector<double> Wcol(N * K);
  auto weights = [&] {
    denominators();
    for (std::size_t n = 0; n < N; ++n)
      for (std::size_t k = 0; k < K; ++k) Wcol[n * K + k] = std::exp(f[k] - u[k][n] - logden[n]);
  };
  for (int nt = 0; nt < 50 && idx.size() > 1; ++nt, ++it) {
    weights();
    const std::size_t m = idx.size() - 1;
    std::vector<double> g(m), H(m * m, 0.0);
    double gmax = 0.0;
    for (std::size_t a = 0; a < m; ++a) {
      const std::size_t i = idx[a + 1];
      double s = 0.0;
      for (std::size_t n = 0; n < N; ++n) s += Wcol[n * K + i];
      g[a] = static_cast<double>(nk[i]) * (s - 1.0);
      gmax = std::max(gmax, std::fabs(s - 1.0));
      for (std::size_t b = 0; b < m; ++b) {
        const std::size_t j = idx[b + 1];
        double sij = 0.0;
        for (std::size_t n = 0; n < N; ++n) sij += Wcol[n * K + i] * Wcol[n * K + j];
        H[a * m + b] = -static_cast<double>(nk[i]) * static_cast<double>(nk[j]) * sij + (a == b ? static_cast<double>(nk[i]) * s : 0.0);
      }
    }
    if (gmax < tol) {
      done = true;
      break;
    }
    std::vector<double> step = g;
    if (!gauss_solve(H, step, m)) break;
    // Damped step: halve until the gradient norm drops (the objective is convex, so this terminates).
    const std::vector<double> f_old = f;
    double t = 1.0;
    for (int ls = 0; ls < 30; ++ls, t *= 0.5) {
      f = f_old;
      for (std::size_t a = 0; a < m; ++a) f[idx[a + 1]] -= t * step[a];
      weights();
      double gm = 0.0;
      for (std::size_t a = 0; a < m; ++a) {
        double s = 0.0;
        for (std::size_t n = 0; n < N; ++n) s += Wcol[n * K + idx[a + 1]];
        gm = std::max(gm, std::fabs(s - 1.0));
      }
      if (gm < gmax) break;
    }
  }
  if (idx.size() <= 1) done = true;
  // Final free energies for every state (sampled or not), with f_0 = 0.
  denominators();
  update(f);
  weights();
  r.f = f;
  r.converged = done;
  r.iterations = it;
  // Overlap and covariance.
  std::vector<double> WtW(K * K, 0.0);
  for (std::size_t n = 0; n < N; ++n)
    for (std::size_t i = 0; i < K; ++i) {
      const double wi = Wcol[n * K + i];
      if (wi == 0.0) continue;
      for (std::size_t j = 0; j < K; ++j) WtW[i * K + j] += wi * Wcol[n * K + j];
    }
  r.overlap.assign(K, std::vector<double>(K, 0.0));
  for (std::size_t i = 0; i < K; ++i)
    for (std::size_t j = 0; j < K; ++j) r.overlap[i][j] = WtW[i * K + j] * static_cast<double>(nk[j]);
  // Θ over the sampled states: C = N^½ WᵀW N^½ = V Λ Vᵀ; Θ = N^-½ V diag(λ / (1 − λ)) Vᵀ N^-½, the λ = 1
  // eigenvector (the normalisation, a null direction of the free energies) dropped.
  const std::size_t S = idx.size();
  std::vector<double> C(S * S), lam, V;
  for (std::size_t a = 0; a < S; ++a)
    for (std::size_t b = 0; b < S; ++b)
      C[a * S + b] = std::sqrt(static_cast<double>(nk[idx[a]]) * static_cast<double>(nk[idx[b]])) * WtW[idx[a] * K + idx[b]];
  jacobi(C, S, lam, V);
  std::vector<double> Th(S * S, 0.0);
  for (std::size_t e = 0; e < S; ++e) {
    const double d = 1.0 - lam[e];
    if (d < 1e-10) continue;
    const double s = lam[e] / d;
    for (std::size_t a = 0; a < S; ++a)
      for (std::size_t b = 0; b < S; ++b) Th[a * S + b] += V[a * S + e] * s * V[b * S + e];
  }
  for (std::size_t a = 0; a < S; ++a)
    for (std::size_t b = 0; b < S; ++b) Th[a * S + b] /= std::sqrt(static_cast<double>(nk[idx[a]]) * static_cast<double>(nk[idx[b]]));
  std::vector<int> pos(K, -1);
  for (std::size_t a = 0; a < S; ++a) pos[idx[a]] = static_cast<int>(a);
  r.df.assign(K, std::vector<double>(K, 0.0));
  r.ddf.assign(K, std::vector<double>(K, 0.0));
  for (std::size_t i = 0; i < K; ++i)
    for (std::size_t j = 0; j < K; ++j) {
      r.df[i][j] = f[j] - f[i];
      if (i == j) continue;
      if (pos[i] < 0 || pos[j] < 0) {
        r.ddf[i][j] = kNaN;
        continue;
      }
      const std::size_t a = static_cast<std::size_t>(pos[i]), b = static_cast<std::size_t>(pos[j]);
      const double v = Th[a * S + a] + Th[b * S + b] - 2.0 * Th[a * S + b];
      r.ddf[i][j] = v > 0 ? std::sqrt(v) : 0.0;
    }
  return r;
}

FreeEnergyEstimate thermodynamic_integration(const std::vector<double>& lam, const std::vector<double>& mean,
                                             const std::vector<double>& err) {
  FreeEnergyEstimate r;
  if (lam.size() < 2 || mean.size() != lam.size() || err.size() != lam.size()) return {kNaN, kNaN};
  double var = 0.0;
  for (std::size_t k = 0; k < lam.size(); ++k) {
    const double wl = k > 0 ? 0.5 * (lam[k] - lam[k - 1]) : 0.0;
    const double wr = k + 1 < lam.size() ? 0.5 * (lam[k + 1] - lam[k]) : 0.0;
    const double w = wl + wr;
    r.df += w * mean[k];
    var += w * w * err[k] * err[k];
  }
  r.error = std::sqrt(var);
  return r;
}

double statistical_inefficiency(const std::vector<double>& a) {
  const std::size_t n = a.size();
  if (n < 4) return 1.0;
  const double mu = mean_of(a);
  double c0 = 0.0;
  for (double x : a) c0 += (x - mu) * (x - mu);
  c0 /= static_cast<double>(n);
  if (!(c0 > 0)) return 1.0;
  double g = 1.0;
  // The sum is stopped at the first non-positive autocorrelation, but never before t = 3 (as in pymbar).
  for (std::size_t t = 1; t < n - 1; ++t) {
    double c = 0.0;
    for (std::size_t i = 0; i + t < n; ++i) c += (a[i] - mu) * (a[i + t] - mu);
    c /= static_cast<double>(n - t) * c0;
    if (c <= 0.0 && t > 3) break;
    g += 2.0 * c * (1.0 - static_cast<double>(t) / static_cast<double>(n));
  }
  return std::max(1.0, g);
}

std::vector<std::size_t> subsample_indices(std::size_t n, double g) {
  const std::size_t stride = static_cast<std::size_t>(std::ceil(std::max(1.0, g)));
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < n; i += stride) out.push_back(i);
  return out;
}

// ---- alchemical Lennard-Jones particle --------------------------------------------------------------------

namespace {

struct LjSystem {
  std::size_t N = 0;
  double L = 0, rc2 = 0, ushift = 0, alpha = 0.5;
  std::vector<double> x, v, f;

  void lattice(int cells) {
    const double a = L / cells;
    const double basis[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    std::size_t p = 0;
    for (int i = 0; i < cells; ++i)
      for (int j = 0; j < cells; ++j)
        for (int k = 0; k < cells; ++k)
          for (const auto& b : basis) {
            if (p >= N) return;
            x[3 * p] = (i + b[0] + 0.25) * a;
            x[3 * p + 1] = (j + b[1] + 0.25) * a;
            x[3 * p + 2] = (k + b[2] + 0.25) * a;
            ++p;
          }
  }
  double min_image(double d) const { return d - L * std::round(d / L); }
  double r2_between(std::size_t i, std::size_t j) const {
    double s = 0.0;
    for (int k = 0; k < 3; ++k) {
      const double d = min_image(x[3 * i + k] - x[3 * j + k]);
      s += d * d;
    }
    return s;
  }
};

// Soft-core pair energy and its λ derivative at r² (truncated and shifted at rc).
inline void softcore(double r2, double lam, double alpha, double rc2, double& u, double& dudl) {
  auto raw = [&](double rr2, double& e, double& de) {
    const double r6 = rr2 * rr2 * rr2;
    const double D = std::max(alpha * (1.0 - lam) + r6, 1e-300);  // no 0/0 at a zero distance
    const double iD = 1.0 / D, iD2 = iD * iD;
    e = 4.0 * lam * (iD2 - iD);
    de = 4.0 * (iD2 - iD) + 4.0 * lam * alpha * (2.0 * iD2 * iD - iD2);
  };
  double e, de, ec, dec;
  raw(r2, e, de);
  raw(rc2, ec, dec);
  u = e - ec;
  dudl = de - dec;
}

}  // namespace

std::optional<LjSolvationResult> lj_solvation(const LjSolvationOptions& o, Diagnostics& diag) {
  if (o.cells < 2 || o.cells > 6 || !(o.density > 0 && o.density < 1.2) || !(o.temperature > 0) || !(o.cutoff > 1.0) ||
      !(o.dt > 0 && o.dt < 0.02) || !(o.gamma > 0) || o.lambdas.size() < 2 || o.production < 10 * o.sample_every ||
      o.sample_every == 0) {
    diag.error(DiagCode::InvalidInput,
               "LJ solvation: 2-6 cells, 0 < density < 1.2, T > 0, cutoff > 1, 0 < dt < 0.02, gamma > 0, at least two "
               "lambda values, and at least 10 samples per window");
    return std::nullopt;
  }
  std::vector<double> lam = o.lambdas;
  if (!std::is_sorted(lam.begin(), lam.end()) || lam.front() != 0.0 || lam.back() != 1.0) {
    diag.error(DiagCode::InvalidInput, "lambda values must ascend from 0 to 1");
    return std::nullopt;
  }
  const std::size_t N = static_cast<std::size_t>(4 * o.cells * o.cells * o.cells);
  const double L = std::cbrt(static_cast<double>(N) / o.density);
  if (o.cutoff > 0.5 * L) {
    diag.error(DiagCode::InvalidInput, "the cutoff must not exceed half the box length (use more cells)");
    return std::nullopt;
  }
  const std::size_t K = lam.size();
  const double T = o.temperature, beta = 1.0 / T;
  const double rc2 = o.cutoff * o.cutoff, ir6c = 1.0 / (rc2 * rc2 * rc2), ushift = 4.0 * ir6c * (ir6c - 1.0);
  LjSolvationResult res;
  res.particles = N;
  res.box = L;

  // Per window: the solute energy at every λ (reduced) and ∂u/∂λ at its own λ, per stored sample.
  struct WindowData {
    std::vector<std::vector<double>> u;  // [state][sample]
    std::vector<double> dudl;
    bool ok = true;
  };
  std::vector<WindowData> data(K);

  auto simulate = [&](bool& ok, std::size_t natoms, bool solute, double lw, std::uint64_t seed,
                      const std::function<void(const LjSystem&)>& sample) {
    LjSystem s;
    s.N = natoms;
    s.L = L;
    s.rc2 = rc2;
    s.alpha = o.alpha;
    s.x.assign(3 * natoms, 0.0);
    s.v.assign(3 * natoms, 0.0);
    s.f.assign(3 * natoms, 0.0);
    s.lattice(o.cells);
    Rng rng(seed);
    for (double& q : s.v) q = std::sqrt(T) * rng.normal();
    double U = 0.0;
    auto forces = [&] {
      std::fill(s.f.begin(), s.f.end(), 0.0);
      U = 0.0;
      for (std::size_t i = 0; i < natoms; ++i)
        for (std::size_t j = i + 1; j < natoms; ++j) {
          double d[3], r2 = 0.0;
          for (int k = 0; k < 3; ++k) {
            d[k] = s.min_image(s.x[3 * i + k] - s.x[3 * j + k]);
            r2 += d[k] * d[k];
          }
          if (r2 >= rc2) continue;
          double ff;
          if (solute && i == 0) {
            // Soft-core: F·r̂ / r = 24 λ r⁴ (2/D³ − 1/D²).
            const double r6 = r2 * r2 * r2, D = std::max(o.alpha * (1.0 - lw) + r6, 1e-300), iD = 1.0 / D;
            double e, de;
            softcore(r2, lw, o.alpha, rc2, e, de);
            U += e;
            ff = 24.0 * lw * r2 * r2 * (2.0 * iD * iD * iD - iD * iD);
          } else {
            const double ir2 = 1.0 / r2, ir6 = ir2 * ir2 * ir2;
            U += 4.0 * ir6 * (ir6 - 1.0) - ushift;
            ff = 48.0 * ir6 * (ir6 - 0.5) * ir2;
          }
          for (int k = 0; k < 3; ++k) {
            s.f[3 * i + k] += ff * d[k];
            s.f[3 * j + k] -= ff * d[k];
          }
        }
    };
    forces();
    const double dt = o.dt, c1 = std::exp(-o.gamma * dt), c2 = std::sqrt((1.0 - c1 * c1) * T);
    const std::size_t total = o.equilibration + o.production;
    for (std::size_t step = 1; step <= total; ++step) {
      // BAOAB
      for (std::size_t k = 0; k < 3 * natoms; ++k) s.v[k] += 0.5 * dt * s.f[k], s.x[k] += 0.5 * dt * s.v[k];
      for (double& q : s.v) q = c1 * q + c2 * rng.normal();
      for (std::size_t k = 0; k < 3 * natoms; ++k) s.x[k] += 0.5 * dt * s.v[k];
      forces();
      for (std::size_t k = 0; k < 3 * natoms; ++k) s.v[k] += 0.5 * dt * s.f[k];
      if (!std::isfinite(U)) {
        ok = false;
        return;
      }
      if (step > o.equilibration && (step - o.equilibration) % o.sample_every == 0) sample(s);
    }
  };

  parallel_for(K, [&](std::size_t w) {
    WindowData& D = data[w];
    D.u.assign(K, {});
    simulate(D.ok, N, true, lam[w], o.seed * 7919ULL + w * 104729ULL + 1, [&](const LjSystem& s) {
      std::vector<double> e(K, 0.0);
      double dl = 0.0;
      for (std::size_t j = 1; j < N; ++j) {
        const double r2 = s.r2_between(0, j);
        if (r2 >= rc2) continue;
        for (std::size_t k = 0; k < K; ++k) {
          double uu, de;
          softcore(r2, lam[k], o.alpha, rc2, uu, de);
          e[k] += uu;
          if (k == w) dl += de;
        }
      }
      for (std::size_t k = 0; k < K; ++k) D.u[k].push_back(beta * e[k]);
      D.dudl.push_back(beta * dl);
    });
  });
  for (std::size_t w = 0; w < K; ++w)
    if (!data[w].ok) {
      diag.error(DiagCode::NotConverged, "the dynamics became unstable at lambda " + std::to_string(lam[w]) + "; use a shorter time step");
      return std::nullopt;
    }

  // Subsample each window by the statistical inefficiency of its energy differences to the neighbours.
  std::vector<std::vector<std::size_t>> keep(K);
  for (std::size_t w = 0; w < K; ++w) {
    const WindowData& D = data[w];
    const std::size_t n = D.dudl.size();
    double g = statistical_inefficiency(D.dudl);
    for (std::size_t nb : {w > 0 ? w - 1 : w, w + 1 < K ? w + 1 : w}) {
      if (nb == w) continue;
      std::vector<double> du(n);
      for (std::size_t i = 0; i < n; ++i) du[i] = D.u[nb][i] - D.u[w][i];
      g = std::max(g, statistical_inefficiency(du));
    }
    keep[w] = subsample_indices(n, g);
    LjSolvationWindow win;
    win.lambda = lam[w];
    win.inefficiency = g;
    win.samples = keep[w].size();
    std::vector<double> sub;
    for (std::size_t i : keep[w]) sub.push_back(D.dudl[i]);
    win.dudl = mean_of(sub);
    double var = 0.0;
    for (double q : sub) var += (q - win.dudl) * (q - win.dudl);
    win.dudl_error = sub.size() > 1 ? std::sqrt(var / static_cast<double>(sub.size() - 1) / static_cast<double>(sub.size())) : 0.0;
    res.windows.push_back(win);
  }

  // MBAR on the subsampled samples.
  std::vector<std::size_t> nk(K);
  std::vector<std::vector<double>> u(K);
  for (std::size_t w = 0; w < K; ++w) {
    nk[w] = keep[w].size();
    for (std::size_t k = 0; k < K; ++k)
      for (std::size_t i : keep[w]) u[k].push_back(data[w].u[k][i]);
  }
  res.mbar_detail = mbar(u, nk);
  if (o.keep_matrix) res.u_kn = u, res.n_k = nk;
  res.mbar_converged = res.mbar_detail.converged;
  res.mbar = {res.mbar_detail.df[0][K - 1], res.mbar_detail.ddf[0][K - 1]};
  res.min_neighbour_overlap = 1.0;
  for (std::size_t k = 0; k + 1 < K; ++k)
    res.min_neighbour_overlap = std::min(res.min_neighbour_overlap, res.mbar_detail.overlap[k][k + 1]);

  // BAR summed over neighbouring windows; EXP from λ = 0 straight to 1 (a single-step estimate: poor overlap
  // makes it biased, which is why it is shown).
  double bsum = 0.0, bvar = 0.0;
  for (std::size_t k = 0; k + 1 < K; ++k) {
    std::vector<double> wf, wr;
    for (std::size_t i : keep[k]) wf.push_back(data[k].u[k + 1][i] - data[k].u[k][i]);
    for (std::size_t i : keep[k + 1]) wr.push_back(data[k + 1].u[k][i] - data[k + 1].u[k + 1][i]);
    const FreeEnergyEstimate b = bar(wf, wr);
    bsum += b.df;
    bvar += b.error * b.error;
  }
  res.bar = {bsum, std::sqrt(bvar)};
  {
    std::vector<double> w0;
    for (std::size_t i : keep[0]) w0.push_back(data[0].u[K - 1][i] - data[0].u[0][i]);
    res.exp_forward = exp_average(w0);
  }
  std::vector<double> means, errs;
  for (const auto& win : res.windows) means.push_back(win.dudl), errs.push_back(win.dudl_error);
  res.ti = thermodynamic_integration(lam, means, errs);

  // Widom insertion in the pure fluid of N − 1 particles, same box.
  if (o.insertions > 0) {
    std::vector<double> frame_means;
    Rng ins(o.seed * 31337ULL + 99);
    bool widom_ok = true;
    simulate(widom_ok, N - 1, false, 1.0, o.seed * 7919ULL + 424242ULL, [&](const LjSystem& s) {
      double acc = 0.0;
      for (std::size_t t = 0; t < o.insertions; ++t) {
        const double p[3] = {ins.uniform() * L, ins.uniform() * L, ins.uniform() * L};
        double dU = 0.0;
        for (std::size_t j = 0; j < s.N; ++j) {
          double r2 = 0.0;
          for (int k = 0; k < 3; ++k) {
            const double d = s.min_image(p[k] - s.x[3 * j + k]);
            r2 += d * d;
          }
          if (r2 >= rc2) continue;
          const double ir2 = 1.0 / r2, ir6 = ir2 * ir2 * ir2;
          dU += 4.0 * ir6 * (ir6 - 1.0) - ushift;
        }
        acc += std::exp(-beta * dU);
      }
      frame_means.push_back(acc / static_cast<double>(o.insertions));
    });
    if (!widom_ok || frame_means.size() < 2) {
      diag.warn(DiagCode::NotConverged, "the Widom run did not finish; no insertion estimate");
      res.widom = {kNaN, kNaN};
    } else {
      // Block average of the per-frame Boltzmann factors, then −ln of the mean (error by the delta method).
      const double g = statistical_inefficiency(frame_means);
      const double m = mean_of(frame_means);
      double var = 0.0;
      for (double q : frame_means) var += (q - m) * (q - m);
      var /= static_cast<double>(frame_means.size() - 1);
      const double neff = static_cast<double>(frame_means.size()) / g;
      res.widom = {-std::log(m), std::sqrt(var / neff) / m};
    }
  }
  if (res.min_neighbour_overlap < 0.03)
    diag.warn(DiagCode::NotConverged, "neighbouring lambda windows overlap poorly (minimum " + std::to_string(res.min_neighbour_overlap) +
                                          "); add windows where the overlap is small");
  return res;
}

// ---- JSON ------------------------------------------------------------------------------------------------

namespace {

std::string jn(double v) {
  if (!std::isfinite(v)) return "null";
  char b[32];
  std::snprintf(b, sizeof b, "%.10g", v);
  return b;
}

std::string est(const FreeEnergyEstimate& e) { return "{\"df\":" + jn(e.df) + ",\"error\":" + jn(e.error) + "}"; }

std::string matrix(const std::vector<std::vector<double>>& m) {
  std::string s = "[";
  for (std::size_t i = 0; i < m.size(); ++i) {
    s += i ? ",[" : "[";
    for (std::size_t j = 0; j < m[i].size(); ++j) s += (j ? "," : "") + jn(m[i][j]);
    s += "]";
  }
  return s + "]";
}

}  // namespace

std::string mbar_json(std::string_view text) {
  Diagnostics d;
  std::istringstream in{std::string(text)};
  std::string line;
  std::vector<std::vector<double>> rows;
  std::size_t lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    const auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#') continue;
    std::istringstream ls(line);
    std::vector<double> r;
    std::string tok;
    while (ls >> tok) {
      char* end = nullptr;
      const double v = std::strtod(tok.c_str(), &end);
      if (end == tok.c_str() || *end != '\0' || !std::isfinite(v)) {
        d.error(DiagCode::InvalidInput, "line " + std::to_string(lineno) + ": '" + tok + "' is not a number");
        return "{\"ok\":false,\"diagnostics\":" + diagnostics_to_json(d) + "}";
      }
      r.push_back(v);
    }
    rows.push_back(r);
  }
  if (rows.size() < 3) {
    d.error(DiagCode::InvalidInput, "need a line of sample counts and at least two states");
    return "{\"ok\":false,\"diagnostics\":" + diagnostics_to_json(d) + "}";
  }
  const std::size_t K = rows[0].size();
  std::vector<std::size_t> nk;
  std::size_t total = 0;
  for (double x : rows[0]) {
    if (x < 0 || x != std::floor(x)) {
      d.error(DiagCode::InvalidInput, "the first line must list a non-negative integer sample count per state");
      return "{\"ok\":false,\"diagnostics\":" + diagnostics_to_json(d) + "}";
    }
    nk.push_back(static_cast<std::size_t>(x));
    total += nk.back();
  }
  std::vector<std::vector<double>> u(rows.begin() + 1, rows.end());
  if (u.size() != K) {
    d.error(DiagCode::InvalidInput, std::to_string(K) + " sample counts but " + std::to_string(u.size()) + " rows of reduced potentials");
    return "{\"ok\":false,\"diagnostics\":" + diagnostics_to_json(d) + "}";
  }
  for (std::size_t k = 0; k < K; ++k)
    if (u[k].size() != total) {
      d.error(DiagCode::InvalidInput, "row " + std::to_string(k + 1) + " has " + std::to_string(u[k].size()) +
                                          " values; the sample counts add up to " + std::to_string(total));
      return "{\"ok\":false,\"diagnostics\":" + diagnostics_to_json(d) + "}";
    }
  const MbarResult r = mbar(u, nk);
  if (!r.converged) d.warn(DiagCode::NotConverged, "MBAR did not converge to the requested tolerance");
  double minov = 1.0;
  for (std::size_t k = 0; k + 1 < K; ++k) minov = std::min(minov, r.overlap[k][k + 1]);
  if (minov < 0.03) d.warn(DiagCode::NotConverged, "neighbouring states overlap poorly (minimum " + std::to_string(minov) + ")");
  std::string s = "{\"ok\":true,\"converged\":" + std::string(r.converged ? "true" : "false") + ",\"iterations\":" +
                  std::to_string(r.iterations) + ",\"f\":[";
  for (std::size_t k = 0; k < K; ++k) s += (k ? "," : "") + jn(r.f[k]);
  s += "],\"df\":" + matrix(r.df) + ",\"ddf\":" + matrix(r.ddf) + ",\"overlap\":" + matrix(r.overlap) +
       ",\"min_neighbour_overlap\":" + jn(minov);
  return s + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

std::string lj_solvation_json(const LjSolvationOptions& o) {
  Diagnostics d;
  const auto r = lj_solvation(o, d);
  std::string s = std::string("{\"ok\":") + (r ? "true" : "false");
  if (r) {
    const double T = o.temperature;
    auto both = [&](const FreeEnergyEstimate& e) {
      return "{\"kT\":" + est(e) + ",\"epsilon\":" + est({e.df * T, e.error * T}) + "}";
    };
    s += ",\"particles\":" + std::to_string(r->particles) + ",\"box\":" + jn(r->box) + ",\"temperature\":" + jn(T) +
         ",\"density\":" + jn(o.density) + ",\"mu_ex\":{\"mbar\":" + both(r->mbar) + ",\"bar\":" + both(r->bar) +
         ",\"ti\":" + both(r->ti) + ",\"exp\":" + both(r->exp_forward) + ",\"widom\":" + both(r->widom) +
         "},\"mbar_converged\":" + (r->mbar_converged ? "true" : "false") +
         ",\"min_neighbour_overlap\":" + jn(r->min_neighbour_overlap) + ",\"windows\":[";
    for (std::size_t k = 0; k < r->windows.size(); ++k) {
      const auto& w = r->windows[k];
      s += std::string(k ? "," : "") + "{\"lambda\":" + jn(w.lambda) + ",\"samples\":" + std::to_string(w.samples) +
           ",\"inefficiency\":" + jn(w.inefficiency) + ",\"dudl\":" + jn(w.dudl) + ",\"dudl_error\":" + jn(w.dudl_error) + "}";
    }
    s += "],\"overlap\":" + matrix(r->mbar_detail.overlap);
    if (o.keep_matrix) {
      s += ",\"n_k\":[";
      for (std::size_t k = 0; k < r->n_k.size(); ++k) s += (k ? "," : "") + std::to_string(r->n_k[k]);
      s += "],\"u_kn\":" + matrix(r->u_kn);
    }
  }
  return s + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

}  // namespace chemlarp
