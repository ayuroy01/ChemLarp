// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
// Ewald / PME and explicit-water simulation.
//
// PME is checked against the direct reciprocal-space sum (and that against finite differences), the full Ewald sum
// against the rock-salt Madelung constant 1.747565 (Kittel), the FFT against a direct DFT. The water tests are
// short smoke runs; scripts/validate_hydration.py has the long runs against published TIP3P results.
#include <cmath>
#include <complex>
#include <numbers>
#include <random>
#include <string>
#include <vector>

#include "chemlarp/ewald.hpp"
#include "chemlarp/solvation.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

TEST(fft_matches_direct_dft) {
  std::mt19937_64 g(3);
  std::uniform_real_distribution<double> u(-1, 1);
  for (std::size_t n : {1u, 2u, 8u, 12u, 30u, 45u, 64u, 75u}) {
    CHECK_EQ(fft_size(n), n);
    std::vector<std::complex<double>> a(n);
    for (auto& z : a) z = {u(g), u(g)};
    for (int sign : {-1, 1}) {
      auto b = a;
      fft_inplace(b, sign);
      double err = 0;
      for (std::size_t k = 0; k < n; ++k) {
        std::complex<double> s = 0;
        for (std::size_t j = 0; j < n; ++j) s += a[j] * std::polar(1.0, sign * 2.0 * std::numbers::pi * double(j * k % n) / double(n));
        err = std::max(err, std::abs(s - b[k]));
      }
      CHECK(err < 1e-12 * double(n));
    }
  }
  CHECK_EQ(fft_size(7), 8u);
  CHECK_EQ(fft_size(49), 50u);
  CHECK_EQ(fft_size(97), 100u);
}

TEST(pme_matches_exact_ewald) {
  std::mt19937_64 g(1);
  std::uniform_real_distribution<double> u(0, 1);
  const std::array<double, 3> L = {20.0, 22.0, 24.0};
  const int N = 60;
  std::vector<double> x(3 * N), q(N);
  for (int i = 0; i < N; ++i) {
    for (int d = 0; d < 3; ++d) x[3 * i + d] = u(g) * L[d];
    q[i] = (i % 2 ? 1 : -1) * (0.3 + u(g));
  }
  double sum = 0;
  for (double v : q) sum += v;
  for (double& v : q) v -= sum / N;
  const double beta = ewald_beta(9.0, 1e-5);
  CHECK_NEAR(std::erfc(beta * 9.0), 1e-5, 1e-9);
  std::vector<double> fe(3 * N, 0.0);
  const double Ee = ewald_reciprocal_exact(L, beta, x, q, 14, &fe);
  for (int c : {0, 7, 100, 151}) {
    auto xp = x, xm = x;
    const double h = 1e-4;
    xp[c] += h;
    xm[c] -= h;
    const double fd = -(ewald_reciprocal_exact(L, beta, xp, q, 14, nullptr) - ewald_reciprocal_exact(L, beta, xm, q, 14, nullptr)) / (2 * h);
    CHECK(std::fabs(fd - fe[c]) < 1e-5);
  }
  Pme pme(L, beta, 6, 1.0);
  std::vector<double> fp(3 * N, 0.0);
  const double Ep = pme.energy(x, q, &fp);
  CHECK(std::fabs(Ep - Ee) < 1e-3 * std::fabs(Ee) + 1e-3);
  double fmax = 0;
  for (int i = 0; i < 3 * N; ++i) fmax = std::max(fmax, std::fabs(fp[i] - fe[i]));
  CHECK(fmax < 0.02);
  // A finer grid converges further.
  Pme fine(L, beta, 6, 0.6);
  CHECK(std::fabs(fine.energy(x, q, nullptr) - Ee) < std::max(1e-6, std::fabs(Ep - Ee)));
}

TEST(ewald_gives_rock_salt_madelung_constant) {
  const std::array<double, 3> C = {2, 2, 2};
  std::vector<double> xn, qn;
  for (int a = 0; a < 2; ++a)
    for (int b = 0; b < 2; ++b)
      for (int c = 0; c < 2; ++c) {
        xn.insert(xn.end(), {double(a), double(b), double(c)});
        qn.push_back(((a + b + c) % 2) ? 1.0 : -1.0);
      }
  const double bt = 2.5;
  double Er = 0;
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j)
      for (int nx = -4; nx <= 4; ++nx)
        for (int ny = -4; ny <= 4; ++ny)
          for (int nz = -4; nz <= 4; ++nz) {
            if (i == j && !nx && !ny && !nz) continue;
            const double dx = xn[3 * i] - xn[3 * j] + 2 * nx, dy = xn[3 * i + 1] - xn[3 * j + 1] + 2 * ny,
                         dz = xn[3 * i + 2] - xn[3 * j + 2] + 2 * nz, r = std::sqrt(dx * dx + dy * dy + dz * dz);
            Er += 0.5 * kCoulomb * qn[i] * qn[j] * std::erfc(bt * r) / r;
          }
  const double exact = -(Er + ewald_reciprocal_exact(C, bt, xn, qn, 10, nullptr) + ewald_self(bt, qn)) / (4 * kCoulomb);
  CHECK_NEAR(exact, 1.747565, 1e-6);
  Pme pme(C, bt, 6, 0.05);
  CHECK_NEAR(-(Er + pme.energy(xn, qn, nullptr) + ewald_self(bt, qn)) / (4 * kCoulomb), 1.747565, 1e-4);
}

TEST(tip3p_water_box_short_run) {
  SolvationOptions o;
  o.waters = 216;
  o.equilibration_ps = 1.0;
  Diagnostics d;
  const auto r = simulate_water(o, 2.0, d);
  CHECK(r.has_value());
  if (!r) return;
  CHECK_EQ(r->waters, 216);
  CHECK(r->density > 0.9 && r->density < 1.06);
  CHECK(r->potential_per_molecule > -10.5 && r->potential_per_molecule < -8.5);
  CHECK(r->temperature > 270 && r->temperature < 340);
  double peak = 0, at = 0;
  for (const auto& [rr, gg] : r->g_oo)
    if (gg > peak) peak = gg, at = rr;
  CHECK(at > 2.6 && at < 2.95);  // first O–O shell
  CHECK(peak > 2.0);
}

TEST(hydration_rejects_bad_input) {
  SolvationOptions o;
  o.waters = 300;
  Diagnostics d1, d2, d3;
  CHECK(!hydration_free_energy("[NH4+]", o, d1).has_value());
  CHECK(d1.has_errors());
  CHECK(!hydration_free_energy("not a smiles((", o, d2).has_value());
  o.waters = 100;
  CHECK(!hydration_free_energy("CO", o, d3).has_value());
  const std::string j = hydration_json("[Na+]", o);
  CHECK(j.rfind("{\"ok\":false", 0) == 0);
}

TEST(hydration_pipeline_runs) {
  // Three windows and a fraction of a picosecond: checks the plumbing (windows, MBAR/BAR/TI, JSON), not the value.
  SolvationOptions o;
  o.waters = 300;
  o.equilibration_ps = 0.3;
  o.production_ps = 1.0;
  o.sample_every_ps = 0.1;
  o.states = {{1.0, 1.0}, {0.0, 1.0}, {0.0, 0.0}};
  Diagnostics d;
  const auto r = hydration_free_energy("C", o, d);
  CHECK(r.has_value());
  if (!r) return;
  CHECK_EQ(r->windows.size(), 3u);
  CHECK_EQ(r->solute_atoms, 5);
  CHECK(std::isfinite(r->dg_mbar.df) && std::isfinite(r->dg_bar.df) && std::isfinite(r->dg_ti.df));
  CHECK_NEAR(r->dg_elec + r->dg_vdw, r->dg_mbar.df, 1e-6);
  CHECK(std::fabs(r->dg_mbar.df) < 20.0);
}
