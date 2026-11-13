#include <cmath>
#include <string>
#include <vector>

#include "chemlarp/eos.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {
CriticalProperties comp(const char* smiles) {
  auto c = critical_properties(smiles);
  if (!c) throw cltest::Failure(std::string("no critical data for ") + smiles);
  return *c;
}
double psat(const char* smiles, double T, CubicEos e = CubicEos::PengRobinson) {
  Diagnostics d;
  auto p = saturation_pressure(comp(smiles), e, T, d);
  if (!p) throw cltest::Failure(std::string("no saturation pressure for ") + smiles);
  return *p;
}
}  // namespace

// Experimental vapour pressures: NIST Chemistry WebBook (fluid properties).
TEST(eos_vapour_pressures_against_experiment) {
  CHECK_NEAR(psat("CCC", 300.0) / 0.99757e6, 1.0, 0.015);        // propane
  CHECK_NEAR(psat("CCCC", 300.0) / 0.25806e6, 1.0, 0.015);       // n-butane
  CHECK_NEAR(psat("O=C=O", 280.0) / 4.1607e6, 1.0, 0.015);       // carbon dioxide
  CHECK_NEAR(psat("C", 150.0) / 1.0400e6, 1.0, 0.015);           // methane
  CHECK_NEAR(psat("CCCCCCC", 371.53) / 0.101325e6, 1.0, 0.015);  // n-heptane at its normal boiling point
  CHECK_NEAR(psat("O", 373.124) / 0.101325e6, 1.0, 0.06);        // water: PR is known to be ~5% low
  CHECK_NEAR(psat("CCC", 300.0, CubicEos::SoaveRedlichKwong) / 0.99757e6, 1.0, 0.02);
  CHECK_NEAR(psat("CCCC", 300.0, CubicEos::SoaveRedlichKwong) / 0.25806e6, 1.0, 0.02);
}

TEST(eos_ideal_gas_limit) {
  Diagnostics d;
  auto m = mixture_from_smiles({"C", "CCC"}, CubicEos::PengRobinson, d);
  CHECK(m.has_value());
  const auto p = phase_properties(*m, 400.0, 1.0, {0.4, 0.6}, Root::Vapor);
  CHECK_NEAR(p.Z, 1.0, 1e-6);
  for (double l : p.ln_phi) CHECK_NEAR(l, 0.0, 1e-6);
}

// ln φ_i must be the partial molar derivative of n·G^R/RT (thermodynamic consistency).
TEST(eos_fugacity_is_partial_molar_derivative) {
  Diagnostics d;
  for (auto e : {CubicEos::PengRobinson, CubicEos::SoaveRedlichKwong}) {
    auto m = mixture_from_smiles({"C", "CCC", "CCCCC", "O=C=O"}, e, d);
    CHECK(m.has_value());
    const double T = 320.0, P = 2.0e6;
    const std::vector<double> n = {0.3, 0.25, 0.35, 0.1};
    for (Root r : {Root::Liquid, Root::Vapor}) {
      auto nG = [&](const std::vector<double>& nn) {
        double N = 0;
        for (double v : nn) N += v;
        std::vector<double> x(nn.size());
        for (std::size_t i = 0; i < nn.size(); ++i) x[i] = nn[i] / N;
        return N * phase_properties(*m, T, P, x, r).gibbs_residual_RT;
      };
      const auto p = phase_properties(*m, T, P, n, r);
      for (std::size_t i = 0; i < n.size(); ++i) {
        auto up = n, dn = n;
        up[i] += 1e-6;
        dn[i] -= 1e-6;
        CHECK_NEAR((nG(up) - nG(dn)) / 2e-6, p.ln_phi[i], 1e-6);
      }
    }
  }
}

TEST(flash_two_phase_is_consistent) {
  Diagnostics d;
  auto m = mixture_from_smiles({"C", "CCC", "CCCCC"}, CubicEos::PengRobinson, d);
  CHECK(m.has_value());
  const std::vector<double> z = {0.5, 0.3, 0.2};
  const auto f = flash_tp(*m, 300.0, 3.0e6, z, d);
  CHECK(f.converged);
  CHECK(f.region == PhaseRegion::TwoPhase);
  CHECK(f.vapor_fraction > 0.0 && f.vapor_fraction < 1.0);
  double sx = 0, sy = 0;
  for (std::size_t i = 0; i < 3; ++i) {
    sx += f.x[i];
    sy += f.y[i];
    CHECK_NEAR(f.vapor_fraction * f.y[i] + (1 - f.vapor_fraction) * f.x[i], z[i], 1e-10);  // material balance
    CHECK_NEAR(std::log(f.x[i]) + f.liquid.ln_phi[i], std::log(f.y[i]) + f.vapor.ln_phi[i], 1e-9);  // equal fugacity
  }
  CHECK_NEAR(sx, 1.0, 1e-12);
  CHECK_NEAR(sy, 1.0, 1e-12);
  CHECK(f.y[0] > z[0] && f.x[2] > z[2]);  // light component enriched in the vapour, heavy in the liquid
  // Feed scaling does not matter.
  const auto g = flash_tp(*m, 300.0, 3.0e6, {1.5, 0.9, 0.6}, d);
  CHECK_NEAR(g.vapor_fraction, f.vapor_fraction, 1e-10);
}

TEST(flash_single_phase_regions) {
  Diagnostics d;
  auto m = mixture_from_smiles({"C", "CCC", "CCCCC"}, CubicEos::PengRobinson, d);
  CHECK(flash_tp(*m, 300.0, 0.1e6, {0.5, 0.3, 0.2}, d).region == PhaseRegion::Vapor);
  const auto dense = flash_tp(*m, 300.0, 30e6, {0.5, 0.3, 0.2}, d);
  CHECK(dense.region != PhaseRegion::TwoPhase);
  // Pure propane just below / above its vapour pressure.
  auto pr = mixture_from_smiles({"CCC"}, CubicEos::PengRobinson, d);
  const double ps = psat("CCC", 300.0);
  CHECK(flash_tp(*pr, 300.0, 0.95 * ps, {1.0}, d).region == PhaseRegion::Vapor);
  CHECK(flash_tp(*pr, 300.0, 1.05 * ps, {1.0}, d).region == PhaseRegion::Liquid);
  // Hot, compressed CO2 is supercritical.
  auto co2 = mixture_from_smiles({"O=C=O"}, CubicEos::PengRobinson, d);
  CHECK(flash_tp(*co2, 350.0, 10e6, {1.0}, d).region == PhaseRegion::Supercritical);
}

TEST(eos_lookup_and_error_vectors) {
  CHECK(critical_properties("OCC")->name == "ethanol");  // any spelling of the structure
  CHECK(critical_properties("CCCCCCCCCCCCCCCC")->name == "hexadecane");  // from the property database
  CHECK(!critical_properties("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC").has_value());
  Diagnostics d;
  CHECK(!mixture_from_smiles({"C", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"}, CubicEos::PengRobinson, d));
  CHECK(d.contains(DiagCode::NoCriticalData));
  Diagnostics d2;
  CHECK(!saturation_pressure(comp("C"), CubicEos::PengRobinson, 200.0, d2));  // above Tc = 190.56 K
  CHECK(d2.contains(DiagCode::SupercriticalTemperature));
  Diagnostics d3;
  auto m = mixture_from_smiles({"C"}, CubicEos::PengRobinson, d3);
  CHECK(!flash_tp(*m, 300, 1e5, {-1.0}, d3).ok);
  CHECK(d3.contains(DiagCode::NegativeAmount));
}
