#include <cmath>
#include <string>
#include <utility>

#include "chemlarp/gibbs.hpp"
#include "chemlarp/kinetics.hpp"
#include "chemlarp/ode.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

Mechanism mech(const std::string& text) {
  Diagnostics d;
  auto m = parse_mechanism(text, d);
  if (!m) throw cltest::Failure("mechanism failed to parse: " + d.items().front().message);
  return *m;
}

Diagnostics mech_errors(const std::string& text) {
  Diagnostics d;
  parse_mechanism(text, d);
  return d;
}

// Constant-volume reactor whose volume is exactly 1 m^3 (P = N R T), so moles == mol/m^3.
ReactorSetup unit_volume(const Mechanism& m, double T, std::initializer_list<std::pair<const char*, double>> feed) {
  ReactorSetup su;
  su.model = ReactorModel::IsothermalConstantVolume;
  su.temperature_K = T;
  su.moles.assign(m.species_count(), 0.0);
  double N = 0;
  for (auto& [name, v] : feed) {
    su.moles[static_cast<std::size_t>(m.find(name))] = v;
    N += v;
  }
  su.pressure_Pa = N * kGasConstant * T;
  return su;
}

double sample_value(const ReactorResult& r, double t, std::size_t sp) {
  for (const auto& s : r.samples)
    if (std::fabs(s.t - t) <= 1e-12 * std::max(1.0, t)) return s.n[sp];
  throw cltest::Failure("no sample at t = " + std::to_string(t));
}

}  // namespace

TEST(lu_factor_with_late_pivots) {
  // Needs row interchanges at several elimination steps (regression: swaps must be
  // applied to b before forward substitution).
  const std::size_t n = 4;
  std::vector<double> A = {1e-3, 2, 3, 4,  //
                           5, 6e-3, 7, 8,  //
                           9, 10, 1e-3, 12,  //
                           13, 14, 15, 1e-3};
  const std::vector<double> x = {1, -2, 3, -4};
  std::vector<double> b(n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) b[i] += A[i * n + j] * x[j];
  LuFactor lu;
  CHECK(lu.factor(A, n));
  lu.solve(b);
  for (std::size_t i = 0; i < n; ++i) CHECK_NEAR(b[i], x[i], 1e-10);
}

// Prothero-Robinson: y' = L (y - cos t) - sin t, y(0) = 1, exact y = cos t, L = -1e6 (very stiff).
TEST(bdf_prothero_robinson_stiff) {
  const double L = -1e6;
  OdeRhs f = [&](double t, const std::vector<double>& y, std::vector<double>& dy) {
    dy[0] = L * (y[0] - std::cos(t)) - std::sin(t);
  };
  OdeJacobian jac = [&](double, const std::vector<double>&, std::vector<double>& J) { J.assign(1, L); };
  std::vector<double> y{1.0};
  std::vector<double> outs = {1.0, 2.0, 5.0, 10.0};
  double worst = 0;
  OdeObserver obs = [&](double t, const std::vector<double>& yy) { worst = std::max(worst, std::fabs(yy[0] - std::cos(t))); };
  BdfOptions o;
  o.rtol = 1e-8;
  o.atol = {1e-10};
  BdfStats st;
  Diagnostics d;
  CHECK(integrate_bdf(f, &jac, 0.0, y, 10.0, outs, obs, nullptr, o, st, d));
  CHECK(worst < 1e-6);
  CHECK(st.steps < 3000);  // an explicit method would need ~1e7 steps for stability
}

TEST(kinetics_first_order_analytic) {
  auto m = mech("species\n A CCO\n B COC\nend\nreactions J/mol m3\n A => B  2.5 0 0\nend\n");
  auto su = unit_volume(m, 300, {{"A", 1.0}});
  su.output_times = {0.1, 0.5, 1.0, 2.0};
  su.t_end = 2.0;
  su.rtol = 1e-9;
  Diagnostics d;
  auto r = run_reactor(m, su, d);
  CHECK(r.ok);
  for (double t : su.output_times) {
    CHECK_NEAR(sample_value(r, t, 0) / std::exp(-2.5 * t), 1.0, 1e-6);
    CHECK_NEAR(sample_value(r, t, 0) + sample_value(r, t, 1), 1.0, 1e-10);
  }
}

TEST(kinetics_second_order_stoichiometry) {
  // 2 C2H4 => C4H8 with rate k [E]^2: d[E]/dt = -2 k [E]^2  =>  [E] = E0 / (1 + 2 k E0 t)
  auto m = mech("species\n E C=C\n B C=CCC\nend\nreactions J/mol m3\n 2 E => B  0.3 0 0\nend\n");
  auto su = unit_volume(m, 300, {{"E", 2.0}});
  su.output_times = {0.5, 1.0, 4.0};
  su.t_end = 4.0;
  su.rtol = 1e-9;
  Diagnostics d;
  auto r = run_reactor(m, su, d);
  CHECK(r.ok);
  for (double t : su.output_times) CHECK_NEAR(sample_value(r, t, 0) / (2.0 / (1 + 2 * 0.3 * 2.0 * t)), 1.0, 1e-6);
}

// Robertson (1966): the classic stiff kinetics benchmark, with propanol isomers as
// A, B, C so the "isomerisations" are element-balanced. Reference values at t = 40
// from Hairer & Wanner, "Solving ODEs II".
TEST(kinetics_robertson_benchmark) {
  auto m = mech(
      "species\n A CCCO\n B CC(C)O\n C CCOC\nend\n"
      "reactions J/mol m3\n"
      "  A => B      0.04  0 0\n"
      "  2 B => B + C  3e7 0 0\n"
      "  B + C => A + C  1e4 0 0\n"
      "end\n");
  auto su = unit_volume(m, 300, {{"A", 1.0}});
  su.output_times = {40.0};
  su.t_end = 4e5;
  su.rtol = 1e-9;
  Diagnostics d;
  auto r = run_reactor(m, su, d);
  CHECK(r.ok);
  CHECK_NEAR(sample_value(r, 40.0, 0), 0.7158270687, 1e-6);
  CHECK_NEAR(sample_value(r, 40.0, 1) / 9.185534764e-6, 1.0, 1e-4);
  CHECK_NEAR(sample_value(r, 40.0, 2), 0.2841637457, 1e-6);
  double sum = 0;
  for (double v : r.final_state.n) sum += v;
  CHECK_NEAR(sum, 1.0, 1e-9);
  CHECK(r.stats.steps < 2000);  // 4e5 s of a problem with stiffness ratio ~1e11
}

TEST(kinetics_jacobian_matches_finite_difference) {
  Diagnostics d;
  auto m = load_mechanism(cltest::data_path("mechanisms/h2o2.ocm"), d);
  CHECK(m.has_value());
  const std::size_t S = m->species_count();
  std::vector<double> C = {20, 10, 0.3, 0.05, 0.08, 4, 0.002, 1e-4, 50};
  std::vector<double> J(S * S), wp(S), wm(S);
  production_rates_jacobian(*m, 1800.0, C.data(), J.data());
  for (std::size_t j = 0; j < S; ++j) {
    auto cp = C, cm = C;
    const double h = 1e-3 * std::max(C[j], 1.0);  // rates are low-degree polynomials in C: large h is exact enough
    cp[j] += h;
    cm[j] -= h;
    production_rates(*m, 1800.0, cp.data(), wp.data());
    production_rates(*m, 1800.0, cm.data(), wm.data());
    for (std::size_t i = 0; i < S; ++i) {
      const double fd = (wp[i] - wm[i]) / (2 * h);
      if (std::fabs(fd - J[i * S + j]) > 1e-6 * (std::fabs(fd) + std::fabs(J[i * S + j])) + 1e-6)
        throw cltest::Failure("J[" + m->names[i] + "][" + m->names[j] + "] analytic " + std::to_string(J[i * S + j]) +
                              " vs finite difference " + std::to_string(fd));
    }
  }
}

// Detailed balance: a reversible mechanism must relax to the Gibbs-minimiser equilibrium.
TEST(kinetics_relaxes_to_gibbs_equilibrium) {
  auto m = mech(
      "species\n CO [C-]#[O+]\n H2O O\n CO2 O=C=O\n H2 [H][H]\nend\n"
      "reactions kJ/mol cm3\n CO + H2O <=> CO2 + H2   1e12 0 80\nend\n");
  ReactorSetup su;
  su.model = ReactorModel::IsothermalConstantVolume;
  su.temperature_K = 1000;
  su.pressure_Pa = 1e5;
  su.moles = {1.0, 1.0, 0.0, 0.0};
  su.t_end = 10.0;
  Diagnostics d;
  auto r = run_reactor(m, su, d);
  CHECK(r.ok);
  EquilibriumInput in;
  for (std::uint32_t i = 0; i < 4; ++i) {
    in.composition.push_back(m.state.species[i].composition);
    in.g_RT.push_back(m.thermo[i]->g_RT(1000));
  }
  in.n0 = su.moles;
  in.pressure_Pa = 1e5;
  auto eq = minimize_gibbs_ideal_gas(in, d);
  CHECK(eq.converged);
  for (std::size_t i = 0; i < 4; ++i) CHECK_NEAR(r.final_state.n[i], eq.n[i], 1e-6);
}

TEST(kinetics_h2o2_adiabatic_ignition_invariants) {
  Diagnostics d;
  auto m = load_mechanism(cltest::data_path("mechanisms/h2o2.ocm"), d);
  CHECK(m.has_value());
  ReactorSetup su;
  su.model = ReactorModel::AdiabaticConstantVolume;
  su.temperature_K = 1100;
  su.pressure_Pa = 1e5;
  su.moles.assign(m->species_count(), 0.0);
  su.moles[static_cast<std::size_t>(m->find("H2"))] = 2;
  su.moles[static_cast<std::size_t>(m->find("O2"))] = 1;
  su.moles[static_cast<std::size_t>(m->find("AR"))] = 7;
  su.t_end = 0.05;
  su.output_times = {0.05};
  auto r = run_reactor(*m, su, d);
  CHECK(r.ok);
  CHECK(r.ignition_time > 1e-6 && r.ignition_time < 1e-2);
  const auto& f = r.final_state;
  CHECK(f.T > 2500.0);
  CHECK(r.stats.steps < 5000);

  // Elements conserved.
  for (std::size_t z : {1u, 8u, 18u}) {
    double a = 0, b = 0;
    for (std::size_t i = 0; i < m->species_count(); ++i) {
      a += m->state.species[i].composition[z] * su.moles[i];
      b += m->state.species[i].composition[z] * f.n[i];
    }
    CHECK_NEAR(b / a, 1.0, 1e-9);
  }
  // Internal energy conserved (closed, rigid, adiabatic).
  auto U = [&](const std::vector<double>& n, double T) {
    double u = 0;
    for (std::size_t i = 0; i < n.size(); ++i) u += n[i] * (m->thermo[i]->h_RT(T) - 1.0) * kGasConstant * T;
    return u;
  };
  double scale = 0;
  for (std::size_t i = 0; i < f.n.size(); ++i) scale += f.n[i] * m->thermo[i]->cp_R(f.T) * kGasConstant * f.T;
  CHECK(std::fabs(U(f.n, f.T) - U(su.moles, su.temperature_K)) < 1e-6 * scale);

  // The burnt state is the Gibbs equilibrium at its own (T, P).
  EquilibriumInput in;
  for (std::uint32_t i = 0; i < m->species_count(); ++i) {
    in.composition.push_back(m->state.species[i].composition);
    in.g_RT.push_back(m->thermo[i]->g_RT(f.T));
  }
  in.n0 = su.moles;
  in.pressure_Pa = f.P;
  auto eq = minimize_gibbs_ideal_gas(in, d);
  CHECK(eq.converged);
  for (std::size_t i = 0; i < f.n.size(); ++i)
    if (eq.n[i] > 1e-4) CHECK_NEAR(f.n[i] / eq.n[i], 1.0, 1e-4);
}

TEST(kinetics_constant_pressure_expands) {
  // H2O2 => 2 OH at fixed T and P: moles double, so volume doubles at completion.
  auto m = mech("species\n P OO\n OH [OH]\nend\nreactions J/mol m3\n P => 2 OH  50 0 0\nend\n");
  ReactorSetup su;
  su.model = ReactorModel::IsothermalConstantPressure;
  su.temperature_K = 800;
  su.pressure_Pa = 2e5;
  su.moles = {1.0, 0.0};
  su.t_end = 2.0;
  Diagnostics d;
  auto r = run_reactor(m, su, d);
  CHECK(r.ok);
  CHECK_NEAR(r.final_state.P, 2e5, 1e-6);
  CHECK_NEAR(r.final_state.n[1], 2.0, 1e-6);
  CHECK_NEAR(r.final_state.V / r.samples.front().V, 2.0, 1e-6);
}

TEST(mechanism_error_vectors) {
  const std::string sp = "species\n H2 [H][H]\n O2 O=O\n H2O O\nend\n";
  CHECK(mech_errors(sp + "reactions\n H2 + O2 => H2O 1 0 0\nend\n").contains(DiagCode::ElementImbalance));
  CHECK(mech_errors(sp + "reactions\n H2 + X => H2O 1 0 0\nend\n").contains(DiagCode::UnknownSpecies));
  CHECK(mech_errors(sp + "reactions\n 2H2 + O2 => 2H2O -1 0 0\nend\n").contains(DiagCode::InvalidRateParameter));
  CHECK(mech_errors(sp + "reactions\n 2H2 + O2 => 2H2O 1 0\nend\n").contains(DiagCode::MechanismSyntax));
  CHECK(mech_errors(sp + "reactions\n H2 + O2 (+M) => H2O 1 0 0\nend\n").contains(DiagCode::MechanismSyntax));  // (+M) on one side
  CHECK(mech_errors(sp + "reactions\n H2 + M => H2 + O2 1 0 0\nend\n").contains(DiagCode::MechanismSyntax));
  CHECK(mech_errors("species\n A O\n B [OH2]\nend\n").contains(DiagCode::DuplicateSpecies));  // same structure
  CHECK(mech_errors("species\n A O\n A CCO\nend\n").contains(DiagCode::DuplicateSpecies));
  CHECK(mech_errors("species\n A C(C)(C)(C)(C)C\nend\n").contains(DiagCode::ValenceExceeded));
  // Reversible reactions need thermo for K_c; silanes have neither database data nor Joback groups.
  CHECK(mech_errors("species\n A C[Si](C)(C)C\n B CC[SiH](C)C\nend\nreactions\n A <=> B 1 0 0\nend\n")
            .contains(DiagCode::NoThermoData));
  CHECK(!mech_errors("species\n A C[Si](C)(C)C\n B CC[SiH](C)C\nend\nreactions\n A => B 1 0 0\nend\n").has_errors());
  CHECK(mech_errors(sp + "reactions\n 2H2 + O2 => 2H2O 1 0 0\n").contains(DiagCode::MechanismSyntax));  // no 'end'
  // A clean mechanism with CHEMKIN-style coefficients parses.
  CHECK(!mech_errors(sp + "reactions\n 2H2+O2=>2H2O 1e10 0 0\nend\n").has_errors());
}

TEST(reactor_setup_error_vectors) {
  auto m = mech("species\n A CCO\n B COC\nend\nreactions\n A => B 1 0 0\nend\n");
  ReactorSetup su;
  su.moles = {1.0};  // wrong length
  Diagnostics d;
  run_reactor(m, su, d);
  CHECK(d.contains(DiagCode::InvalidReactorSetup));
  su.moles = {-1.0, 0.0};
  Diagnostics d2;
  run_reactor(m, su, d2);
  CHECK(d2.contains(DiagCode::NegativeAmount));
  su.moles = {1.0, 0.0};
  su.model = ReactorModel::AdiabaticConstantVolume;  // Joback species have no heat capacity
  Diagnostics d3;
  run_reactor(m, su, d3);
  CHECK(d3.contains(DiagCode::MissingHeatCapacity));
}
