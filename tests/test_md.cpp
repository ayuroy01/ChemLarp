// Molecular dynamics, thermostats, normal modes and the Lennard-Jones fluid.
//
// Independent references: the Lennard-Jones values come from a Metropolis Monte Carlo simulation of the
// same model (scripts/validate_md.py); equipartition and the harmonic limit are exact results of classical
// statistical mechanics; the time-reversibility and dt² error scaling are properties of velocity Verlet.
#include <cmath>
#include <string>
#include <vector>

#include "chemlarp/chemlarp_c.h"
#include "chemlarp/md.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

MoleculeMdResult mol(const char* smiles, MdOptions o, bool modes = false) {
  Diagnostics d;
  auto r = molecule_md(smiles, o, modes, d);
  CHECK(r.has_value());
  return *r;
}

constexpr double kB = 0.0019872043;

}  // namespace

TEST(md_nve_conserves_energy_and_scales_as_dt_squared) {
  MdOptions o;
  o.thermostat = Thermostat::None;
  o.temperature_K = 300.0;
  double sd[2];
  int k = 0;
  for (double dt : {0.4, 0.2}) {
    o.dt_fs = dt;
    o.steps = std::size_t(1000.0 / dt);
    o.sample_every = 5;
    const auto r = mol("CCO", o);
    CHECK(std::fabs(r.md.drift_kcal_mol_ps) < 0.01);
    sd[k++] = r.md.fluctuation;
  }
  CHECK(sd[0] / sd[1] > 3.5 && sd[0] / sd[1] < 4.5);
}

TEST(md_velocity_verlet_is_time_reversible) {
  // Forward, reverse the velocities, run the same number of steps: back to the start (Morse diatomic).
  const PotentialFunction morse = [](const std::vector<double>& x, std::vector<double>& g) {
    g.assign(6, 0.0);
    double d[3], r2s = 0.0;
    for (int k = 0; k < 3; ++k) d[k] = x[3 + k] - x[k], r2s += d[k] * d[k];
    const double rr = std::sqrt(r2s), e = std::exp(-2.0 * (rr - 1.2));
    const double U = 100.0 * (1 - e) * (1 - e), dUdr = 100.0 * 2 * (1 - e) * 2.0 * e;
    for (int k = 0; k < 3; ++k) g[3 + k] = dUdr * d[k] / rr, g[k] = -g[3 + k];
    return U;
  };
  const std::vector<double> m = {12.0, 16.0}, xs = {0, 0, 0, 1.0, 0.2, 0.1};
  MdOptions mo;
  mo.thermostat = Thermostat::None;
  mo.steps = 2000;
  mo.dt_fs = 0.5;
  mo.sample_every = 2000;
  mo.temperature_K = 500;
  Diagnostics d3;
  auto f = run_md(m, xs, morse, mo, d3);
  CHECK(f.has_value());
  std::vector<double> vb = f->v;
  for (double& q : vb) q = -q;
  auto b = run_md(m, f->x, morse, mo, d3, &vb);
  CHECK(b.has_value());
  for (std::size_t i = 0; i < xs.size(); ++i) CHECK_NEAR(b->x[i], xs[i], 1e-8);
}

TEST(md_runs_are_reproducible_for_a_seed) {
  MdOptions o;
  o.steps = 400;
  o.seed = 11;
  const auto a = mol("CCN", o), b = mol("CCN", o);
  for (std::size_t i = 0; i < a.md.x.size(); ++i) CHECK_EQ(a.md.x[i], b.md.x[i]);
  o.seed = 12;
  const auto c = mol("CCN", o);
  bool differ = false;
  for (std::size_t i = 0; i < a.md.x.size(); ++i) differ = differ || a.md.x[i] != c.md.x[i];
  CHECK(differ);
}

TEST(md_normal_modes_and_harmonic_spectrum) {
  MdOptions o;
  o.thermostat = Thermostat::None;
  o.temperature_K = 1.0;
  o.dt_fs = 0.25;
  o.steps = 40000;
  o.sample_every = 4;
  o.keep_velocities = true;
  const auto r = mol("O", o, true);
  CHECK(r.modes.has_value());
  const auto& f = r.modes->frequencies_cm1;
  CHECK_EQ(f.size(), std::size_t{3});
  CHECK_NEAR(f[0], 1757.4, 0.5);  // UFF bend and stretches (regression values)
  CHECK_NEAR(f[2], 3759.6, 0.5);
  CHECK_EQ(r.md.degrees_of_freedom, std::size_t{3});
  // In the harmonic limit each mode is a peak of the velocity spectrum.
  const auto spec = velocity_spectrum(r.md, r.masses, 1.0, 4000.0, 4001);
  for (double nu : f) {
    std::size_t best = 0;
    for (std::size_t i = std::size_t(nu) - 15; i <= std::size_t(nu) + 15; ++i)
      if (spec[i] > spec[best] || best == 0) best = i;
    CHECK(std::fabs(double(best) - nu) < 6.0);
    CHECK(spec[best] > 0.01);
  }
  // A linear molecule has 3N − 5 modes with a degenerate bend.
  const auto co2 = mol("O=C=O", [] { MdOptions q; q.steps = 10; return q; }(), true);
  CHECK(co2.modes.has_value());
  CHECK_EQ(co2.modes->frequencies_cm1.size(), std::size_t{4});
  CHECK_NEAR(co2.modes->frequencies_cm1[0], co2.modes->frequencies_cm1[1], 0.5);
  CHECK_EQ(co2.md.degrees_of_freedom, std::size_t{4});
}

TEST(md_thermostats_give_equipartition) {
  // Classical harmonic limit: <U − Umin> = (3N − 6) kB T / 2, and the mean kinetic temperature is T.
  for (Thermostat th : {Thermostat::Bussi, Thermostat::Langevin}) {
    MdOptions o;
    o.thermostat = th;
    o.temperature_K = 10.0;
    o.tau_fs = 50.0;
    o.steps = 120000;
    o.sample_every = 20;
    const auto r = mol("O", o);
    const double want = 0.5 * 3 * kB * 10.0;
    CHECK_NEAR(r.md.mean_potential - r.minimized_energy, want, 0.06 * want);
    CHECK_NEAR(r.md.mean_temperature_K, 10.0, 0.4);
  }
}

TEST(md_lennard_jones_against_monte_carlo) {
  // Monte Carlo (scripts/validate_md.py), N = 108, rc = 2.5 shifted, rho 0.8, T 1.0:
  // U/N = −4.702 ± 0.005, P = 1.604 ± 0.027.
  LjOptions o;
  o.equilibration = 2000;
  o.production = 15000;
  Diagnostics d;
  auto r = lennard_jones_fluid(o, d);
  CHECK(r.has_value());
  CHECK_EQ(r->particles, std::size_t{108});
  CHECK_NEAR(r->energy_per_particle, -4.702, 0.04);
  CHECK_NEAR(r->pressure, 1.604, 0.15);
  CHECK_NEAR(r->temperature, 1.0, 0.03);
  // NVE: the energy is conserved.
  o.thermostat = Thermostat::None;
  o.equilibration = 0;
  o.production = 5000;
  o.dt = 0.002;
  auto n = lennard_jones_fluid(o, d);
  CHECK(n.has_value());
  CHECK(std::fabs(n->drift) < 1e-3);
}

TEST(md_input_errors) {
  Diagnostics d;
  const PotentialFunction zero = [](const std::vector<double>& x, std::vector<double>& g) {
    g.assign(x.size(), 0.0);
    return 0.0;
  };
  MdOptions o;
  CHECK(!run_md({1.0}, {0, 0}, zero, o, d).has_value());  // wrong coordinate count
  o.dt_fs = 50;
  CHECK(!run_md({1.0}, {0, 0, 0}, zero, o, d).has_value());
  o.dt_fs = 0.5;
  CHECK(!run_md({-1.0}, {0, 0, 0}, zero, o, d).has_value());
  LjOptions lj;
  lj.cells = 1;  // box too small for the cutoff
  CHECK(!lennard_jones_fluid(lj, d).has_value());
  char* bad = cl_md_json("C1CC", 300, 10, 0.5, "bussi", 100, 1, 1, 0, "uff");
  CHECK(std::string(bad).find("\"ok\":false") != std::string::npos);
  cl_free_string(bad);
  char* good = cl_md_json("CO", 300, 200, 0.5, "langevin", 100, 10, 1, 1, "mmff94");
  const std::string s(good);
  cl_free_string(good);
  CHECK(s.rfind("{\"ok\":true", 0) == 0);
  CHECK(s.find("\"positions\":[[") != std::string::npos);
  CHECK(s.find("\"modes\":[") != std::string::npos);
}
