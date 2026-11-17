// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
// Free-energy estimators against exact results, and the alchemical Lennard-Jones particle.
//
// Harmonic oscillators u_k(x) = K_k (x - x0_k)^2 / 2 (kT units) have f_k = ln(K_k / 2 pi) / 2 exactly.
// scripts/validate_fep.py checks the same estimators against pymbar and Widom insertion.
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "chemlarp/free_energy.hpp"
#include "chemlarp/parallel.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

struct Oscillators {
  std::vector<double> K, x0;
  std::vector<std::vector<double>> xs;  // samples per state
  double u(std::size_t k, double x) const { return 0.5 * K[k] * (x - x0[k]) * (x - x0[k]); }
  std::vector<std::vector<double>> matrix() const {
    std::vector<std::vector<double>> m(K.size());
    for (std::size_t k = 0; k < K.size(); ++k)
      for (const auto& s : xs)
        for (double x : s) m[k].push_back(u(k, x));
    return m;
  }
};

Oscillators draw(std::mt19937_64& g, const std::vector<double>& K, const std::vector<double>& x0, std::vector<std::size_t> n) {
  std::normal_distribution<double> nd(0.0, 1.0);
  Oscillators o{K, x0, {}};
  o.xs.resize(K.size());
  for (std::size_t k = 0; k < K.size(); ++k)
    for (std::size_t i = 0; i < n[k]; ++i) o.xs[k].push_back(x0[k] + nd(g) / std::sqrt(K[k]));
  return o;
}

}  // namespace

TEST(free_energy_estimators_match_exact_oscillators) {
  const std::vector<double> K = {1, 2, 4, 8, 16}, x0 = {0, 0.5, 1.0, 1.5, 2.0};
  std::mt19937_64 g(3);
  const Oscillators o = draw(g, K, x0, {2000, 2000, 2000, 2000, 2000});
  const MbarResult r = mbar(o.matrix(), {2000, 2000, 2000, 2000, 2000});
  CHECK(r.converged);
  for (std::size_t k = 1; k < K.size(); ++k) {
    const double exact = 0.5 * std::log(K[k] / K[0]);
    CHECK(std::fabs(r.df[0][k] - exact) < 4.0 * r.ddf[0][k]);
    CHECK(r.ddf[0][k] > 0.005 && r.ddf[0][k] < 0.1);
    CHECK(std::fabs(r.df[0][k] + r.df[k][0]) < 1e-12);
  }
  // Overlap rows sum to 1.
  for (const auto& row : r.overlap) {
    double s = 0;
    for (double x : row) s += x;
    CHECK(std::fabs(s - 1.0) < 1e-9);
  }
  // BAR and EXP between the first two states.
  std::vector<double> wf, wr;
  for (double x : o.xs[0]) wf.push_back(o.u(1, x) - o.u(0, x));
  for (double x : o.xs[1]) wr.push_back(o.u(0, x) - o.u(1, x));
  const double exact01 = 0.5 * std::log(2.0);
  const FreeEnergyEstimate b = bar(wf, wr), e = exp_average(wf);
  CHECK(std::fabs(b.df - exact01) < 4.0 * b.error);
  CHECK(std::fabs(e.df - exact01) < 4.0 * e.error);
}

TEST(free_energy_bar_equals_two_state_mbar) {
  // With two states MBAR reduces to BAR: the same estimate, and the same asymptotic error.
  std::mt19937_64 g(5);
  const Oscillators o = draw(g, {1.0, 3.0}, {0.0, 0.7}, {700, 400});
  const MbarResult r = mbar(o.matrix(), {700, 400});
  std::vector<double> wf, wr;
  for (double x : o.xs[0]) wf.push_back(o.u(1, x) - o.u(0, x));
  for (double x : o.xs[1]) wr.push_back(o.u(0, x) - o.u(1, x));
  const FreeEnergyEstimate b = bar(wf, wr);
  CHECK(std::fabs(r.df[0][1] - b.df) < 1e-8);
  CHECK(std::fabs(r.ddf[0][1] - b.error) < 1e-6);
}

TEST(free_energy_mbar_edge_cases) {
  // Identical states: zero difference with zero uncertainty.
  std::mt19937_64 g(9);
  Oscillators o = draw(g, {2.0, 2.0}, {0.0, 0.0}, {300, 300});
  MbarResult r = mbar(o.matrix(), {300, 300});
  CHECK(std::fabs(r.df[0][1]) < 1e-10);
  CHECK(r.ddf[0][1] < 1e-6);
  // A state that is evaluated but not sampled: its free energy, no uncertainty.
  o = draw(g, {1.0, 2.0, 4.0}, {0.0, 0.3, 0.6}, {800, 0, 800});
  r = mbar(o.matrix(), {800, 0, 800});
  CHECK(std::fabs(r.df[0][1] - 0.5 * std::log(2.0)) < 0.05);
  CHECK(std::isnan(r.ddf[0][1]));
  CHECK(!std::isnan(r.ddf[0][2]));
  // Inconsistent input: nothing.
  CHECK(mbar({{1.0, 2.0}, {1.0}}, {1, 1}).f.empty());
  CHECK(mbar({{1.0, 2.0}, {1.0, 2.0}}, {1, 2}).f.empty());
}

TEST(free_energy_ti_and_inefficiency) {
  // Trapezoidal rule: exact for a linear integrand, second order otherwise.
  CHECK(std::fabs(thermodynamic_integration({0, 0.3, 1}, {1, 1.6, 3}, {0, 0, 0}).df - 2.0) < 1e-12);
  auto err = [](int m) {
    std::vector<double> l, f, e;
    for (int k = 0; k <= m; ++k) {
      const double x = double(k) / m;
      l.push_back(x), f.push_back(7.5 / (1 + 15 * x)), e.push_back(0.1);
    }
    return std::fabs(thermodynamic_integration(l, f, e).df - 0.5 * std::log(16.0));
  };
  const double e20 = err(20), e40 = err(40);
  CHECK(e40 < e20 / 3.5 && e40 > e20 / 4.5);
  CHECK(std::fabs(thermodynamic_integration({0, 0.5, 1}, {0, 0, 0}, {1, 1, 1}).error - std::sqrt(0.0625 + 0.25 + 0.0625)) < 1e-12);
  // AR(1) with phi = 0.8: g = (1 + phi) / (1 - phi) = 9; the average over 40 series is within 10 %.
  std::mt19937_64 g(17);
  std::normal_distribution<double> nd(0.0, 1.0);
  double sum = 0.0;
  for (int s = 0; s < 40; ++s) {
    std::vector<double> a(8000);
    double y = 0.0;
    for (auto& v : a) y = 0.8 * y + nd(g), v = y;
    sum += statistical_inefficiency(a);
  }
  CHECK(std::fabs(sum / 40 - 9.0) < 0.9);
  std::vector<double> white(4000);
  for (auto& v : white) v = nd(g);
  CHECK(statistical_inefficiency(white) < 1.3);
  CHECK(subsample_indices(10, 2.2) == (std::vector<std::size_t>{0, 3, 6, 9}));
}

TEST(free_energy_lj_solvation_consistent_and_deterministic) {
  LjSolvationOptions o;
  o.lambdas = {0.0, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.7, 1.0};
  o.density = 0.5;
  o.temperature = 1.5;
  o.equilibration = 1000;
  o.production = 6000;
  o.insertions = 50;
  std::string first;
  for (int threads : {1, 4}) {
    set_thread_count(threads);
    Diagnostics d;
    const auto r = lj_solvation(o, d);
    CHECK(r.has_value());
    if (!r) continue;
    CHECK(r->mbar_converged);
    // The estimators agree within their errors (generous: short runs).
    const double s = std::hypot(r->mbar.error, r->ti.error) + 0.05;
    CHECK(std::fabs(r->mbar.df - r->bar.df) < 4 * s);
    CHECK(std::fabs(r->mbar.df - r->ti.df) < 4 * s);
    CHECK(std::fabs(r->mbar.df - r->widom.df) < 4 * (r->mbar.error + r->widom.error) + 0.1);
    CHECK(r->mbar.df < 0.0);  // attractive at this density and temperature
    const std::string j = lj_solvation_json(o);
    if (first.empty()) first = j;
    CHECK(j == first);  // same bits on 1 and 4 threads
  }
  set_thread_count(4);
  // Invalid options are refused.
  LjSolvationOptions bad;
  bad.lambdas = {0.0, 0.5};
  Diagnostics d;
  CHECK(!lj_solvation(bad, d).has_value());
}

TEST(free_energy_mbar_text_input) {
  const std::string ok = "# two states\n2 2\n0 0 1 1\n1 1 0 0\n";
  CHECK(mbar_json(ok).rfind("{\"ok\":true", 0) == 0);
  CHECK(mbar_json("2 2\n0 0 1\n1 1 0 0\n").find("has 3 values") != std::string::npos);
  CHECK(mbar_json("2 x\n0 0 1 1\n1 1 0 0\n").find("not a number") != std::string::npos);
  CHECK(mbar_json("2 2 2\n0 0 1 1\n1 1 0 0\n").find("rows of reduced potentials") != std::string::npos);
}
