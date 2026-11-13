// Pure-component property database and correlations.
//
// Reference values come from sources other than the ones the table was built from:
//  * ideal-gas Cp: NIST-JANAF Thermochemical Tables (1998);
//  * enthalpy of vaporization and saturated liquid volume: VDI Heat Atlas PPDS correlations;
//  * vapour pressure: Wagner equation of McGarry, Ind. Eng. Chem. Process Des. Dev. 22 (1983) 313.
// The tolerances reflect the accuracy of the correlations (see "Pure-component properties: validation" in
// docs/ARCHITECTURE.md), not numerical precision.
#include <cmath>

#include "chemlarp/eos.hpp"
#include "chemlarp/properties.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {
PureComponent pc(const char* smiles) {
  auto c = pure_component(smiles);
  CHECK(c.has_value());
  return *c;
}
double rel(double a, double b) { return a / b - 1.0; }
}  // namespace

TEST(property_lookup) {
  CHECK(pure_component_count() >= 350);
  CHECK(pc("OCC").name == "ethanol");                // any spelling of the structure
  CHECK(pc("C1=CC=CC=C1").cas == "71-43-2");
  CHECK(!pure_component("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC").has_value());
  CHECK(!pure_component("not a smiles").has_value());
  const PureComponent b = pc("c1ccccc1");
  CHECK_NEAR(b.Tc_K, 562.05, 0.01);
  CHECK_NEAR(b.Pc_Pa, 4.895e6, 1.0);
  CHECK_NEAR(b.omega, 0.210, 0.005);
  CHECK_NEAR(b.Tb_K, 353.23, 0.01);
  CHECK_EQ(b.critical_source, 'I');
}

TEST(property_table_consistency_checks) {
  // The IUPAC table has the rows of cyclohexane and methylcyclopentane exchanged; CRC is used instead.
  CHECK_NEAR(pc("C1CCCCC1").Tc_K, 553.4, 0.5);
  CHECK_EQ(pc("C1CCCCC1").critical_source, 'C');
  CHECK_NEAR(pc("CC1CCCC1").Tc_K, 532.8, 0.5);
  // Mis-scaled Cp polynomial (3-methylpentane) dropped; acetylene's kept against JANAF (44.0 J/mol/K).
  Diagnostics d;
  CHECK(!ideal_gas_cp(pc("CCC(C)CC"), 300.0, d));
  CHECK(d.contains(DiagCode::NoPropertyData));
  Diagnostics e;
  CHECK_NEAR(*ideal_gas_cp(pc("C#C"), 298.15, e), 44.04, 0.5);
  // omega is defined from the chosen Tc and Pc. 2-methyl-1-pentanol has no vapour-pressure equation at
  // 0.7 Tc and the PSRK value (0.726) belongs with PSRK's Tc = 582 K, not 604.4 K: no omega, rather
  // than an inconsistent one.
  CHECK(std::isnan(pc("CCCC(C)CO").omega));
  CHECK(!critical_properties("CCCC(C)CO").has_value());
}

TEST(ideal_gas_cp_against_janaf) {
  struct Row { const char* smiles; double T, cp; };
  const Row rows[] = {{"C", 500.0, 46.342}, {"C", 1000.0, 71.795}, {"O", 500.0, 35.226},
                      {"O", 1000.0, 41.268}, {"O=C=O", 500.0, 44.627}, {"O=C=O", 1000.0, 54.308}};
  for (const Row& r : rows) {
    Diagnostics d;
    auto cp = ideal_gas_cp(pc(r.smiles), r.T, d);
    CHECK(cp.has_value());
    CHECK(std::fabs(rel(*cp, r.cp)) < 0.015);  // the polynomials are within 1.3% here, worst at 1000 K
    CHECK(d.empty());
  }
}

TEST(ideal_gas_enthalpy_and_entropy_integrate_cp) {
  const PureComponent b = pc("c1ccccc1");
  Diagnostics d;
  CHECK_NEAR(*ideal_gas_enthalpy(b, 298.15, d), 0.0, 1e-9);
  CHECK_NEAR(*ideal_gas_entropy(b, 298.15, 1e5, d), b.S_gas_JmolK, 1e-9);
  for (double T : {350.0, 600.0, 900.0}) {
    const double h = 1e-3;
    const double dH = (*ideal_gas_enthalpy(b, T + h, d) - *ideal_gas_enthalpy(b, T - h, d)) / (2 * h);
    const double dS = (*ideal_gas_entropy(b, T + h, 1e5, d) - *ideal_gas_entropy(b, T - h, 1e5, d)) / (2 * h);
    const double cp = *ideal_gas_cp(b, T, d);
    CHECK_NEAR(dH, cp, 1e-5 * cp);
    CHECK_NEAR(dS * T, cp, 1e-5 * cp);
  }
  // S(P) = S(1 bar) − R ln(P / 1 bar)
  CHECK_NEAR(*ideal_gas_entropy(b, 400.0, 1e6, d) - *ideal_gas_entropy(b, 400.0, 1e5, d), -8.314462618 * std::log(10.0), 1e-9);
  CHECK(!d.has_errors());
}

TEST(enthalpy_of_vaporization_watson) {
  const PureComponent b = pc("c1ccccc1");
  Diagnostics d;
  CHECK_NEAR(*enthalpy_of_vaporization(b, b.Tb_K, d), b.Hvap_Tb_Jmol, 1e-9);
  CHECK(std::fabs(rel(*enthalpy_of_vaporization(b, 400.0, d), 27862.5)) < 0.02);  // VDI PPDS
  CHECK(*enthalpy_of_vaporization(b, b.Tc_K - 1e-3, d) < 0.01 * b.Hvap_Tb_Jmol);    // vanishes at Tc
  CHECK(!d.has_errors());
  Diagnostics e;
  CHECK(!enthalpy_of_vaporization(b, 600.0, e));
  CHECK(e.contains(DiagCode::SupercriticalTemperature));
  Diagnostics f;
  CHECK(!enthalpy_of_vaporization(pc("C#C"), 250.0, f));  // no tabulated ΔHvap at Tb
  CHECK(f.contains(DiagCode::NoPropertyData));
}

TEST(saturated_liquid_volume) {
  // VDI PPDS density where tabulated: benzene against Perry's DIPPR-105 (independent), water against IAPWS-95.
  Diagnostics d;
  auto v = saturated_liquid_volume(pc("c1ccccc1"), 300.0, d);
  CHECK(v.has_value());
  CHECK(std::fabs(rel(*v, 8.9667e-5)) < 0.005);
  CHECK(std::fabs(rel(18.015268e-3 / *saturated_liquid_volume(pc("O"), 300.0, d), 996.513)) < 0.003);
  CHECK(std::fabs(rel(18.015268e-3 / *saturated_liquid_volume(pc("O"), 450.0, d), 890.341)) < 0.003);
  CHECK(d.empty());
  // Rackett with Spencer–Danner Z_RA where VDI has no equation (2,4-dimethylpentane): no warning; with Zc
  // when no Z_RA is tabulated either (isopropylamine): a warning.
  Diagnostics r;
  CHECK(saturated_liquid_volume(pc("CC(C)CC(C)C"), 290.0, r).has_value());
  CHECK(r.empty());
  CHECK(saturated_liquid_volume(pc("CC(C)N"), 290.0, r).has_value());
  CHECK(r.contains(DiagCode::PropertyExtrapolated));
  // Below the melting point: extrapolation is reported.
  Diagnostics e;
  CHECK(saturated_liquid_volume(pc("c1ccccc1"), 250.0, e).has_value());
  CHECK(e.contains(DiagCode::PropertyExtrapolated));
}

TEST(property_errors_and_warnings) {
  const PureComponent b = pc("c1ccccc1");
  Diagnostics d;
  CHECK(!ideal_gas_cp(b, -1.0, d));
  CHECK(!ideal_gas_entropy(b, 300.0, 0.0, d));
  CHECK_EQ(d.count(DiagCode::InvalidMixtureInput), std::size_t{2});
  Diagnostics e;
  CHECK(ideal_gas_cp(b, 1500.0, e).has_value());  // fitted to 1000 K
  CHECK(e.contains(DiagCode::PropertyExtrapolated));
  CHECK(!e.has_errors());
}

TEST(eos_uses_property_database) {
  // n-Decane is not in the EoS's own table: constants from the database, vapour pressure against
  // McGarry's Wagner equation at 500 K (3.277 bar).
  auto c = critical_properties("CCCCCCCCCC");
  CHECK(c.has_value());
  Diagnostics d;
  auto mix = mixture_from_smiles({"CCCCCCCCCC"}, CubicEos::PengRobinson, d);
  CHECK(mix.has_value());
  auto p = saturation_pressure(mix->components[0], CubicEos::PengRobinson, 500.0, d);
  CHECK(p.has_value());
  CHECK(std::fabs(rel(*p, 327693.6)) < 0.01);
}

// Transport properties: Perry's (DIPPR) correlations against the independent VDI Heat Atlas (PPDS)
// correlations, and for water against REFPROP's IAPWS-based values. Our evaluation of the Perry
// coefficients agrees with thermo 0.6.1's to rounding over ~2,700 points (see ARCHITECTURE.md).
TEST(transport_against_reference_correlations) {
  struct Row { const char* smiles; int prop; double T, ref, tol; };
  // prop: 0 liquid viscosity, 1 liquid conductivity, 2 surface tension, 3 vapour viscosity, 4 vapour conductivity
  const Row rows[] = {
      {"c1ccccc1", 0, 300.0, 5.8613e-4, 0.02}, {"c1ccccc1", 1, 300.0, 0.14279, 0.02},
      {"c1ccccc1", 2, 300.0, 0.02795, 0.02},   {"c1ccccc1", 3, 400.0, 1.0121e-5, 0.02},
      {"c1ccccc1", 4, 400.0, 0.019287, 0.02},  {"CCO", 0, 300.0, 1.0264e-3, 0.02},
      {"CCO", 1, 300.0, 0.16762, 0.02},        {"CCO", 2, 300.0, 0.021882, 0.02},
      {"CCO", 3, 400.0, 1.1775e-5, 0.02},      {"CCO", 4, 400.0, 0.025564, 0.02},
      // Water against REFPROP (IAPWS): Perry's liquid viscosity is 2.7% high at 300 K.
      {"O", 0, 300.0, 8.5365e-4, 0.03},        {"O", 1, 300.0, 0.6095, 0.01},
      {"O", 2, 300.0, 0.071686, 0.01},         {"O", 3, 500.0, 1.7299e-5, 0.02},
      {"O", 4, 500.0, 0.036033, 0.03},
  };
  for (const Row& r : rows) {
    const PureComponent c = pc(r.smiles);
    Diagnostics d;
    std::optional<double> v = r.prop == 0   ? liquid_viscosity(c, r.T, d)
                              : r.prop == 1 ? liquid_thermal_conductivity(c, r.T, d)
                              : r.prop == 2 ? surface_tension(c, r.T, d)
                              : r.prop == 3 ? gas_viscosity(c, r.T, d)
                                            : gas_thermal_conductivity(c, r.T, d);
    CHECK(v.has_value());
    CHECK(std::fabs(rel(*v, r.ref)) < r.tol);
    CHECK(d.empty());
  }
}

TEST(transport_limits_and_errors) {
  const PureComponent b = pc("c1ccccc1");
  Diagnostics d;
  CHECK(!liquid_viscosity(b, 600.0, d));  // above Tc
  CHECK(!surface_tension(b, 600.0, d));
  CHECK_EQ(d.count(DiagCode::SupercriticalTemperature), std::size_t{2});
  CHECK(gas_viscosity(b, 600.0, d).has_value());  // the vapour has no such limit
  // Surface tension falls towards zero at Tc; viscosity of a liquid falls with temperature.
  CHECK(*surface_tension(b, 560.0, d) < 0.05 * *surface_tension(b, 300.0, d));
  CHECK(*liquid_viscosity(b, 350.0, d) < *liquid_viscosity(b, 300.0, d));
  CHECK(*gas_viscosity(b, 500.0, d) > *gas_viscosity(b, 400.0, d));
  Diagnostics e;
  CHECK(gas_viscosity(b, 2000.0, e).has_value());  // fitted to 1000 K
  CHECK(e.contains(DiagCode::PropertyExtrapolated));
  // A missing correlation is reported, not estimated: tert-butanol has no surface-tension correlation,
  // 2-ethylphenol no transport data at all.
  Diagnostics f;
  CHECK(liquid_viscosity(pc("CC(C)(C)O"), 320.0, f).has_value());
  CHECK(!surface_tension(pc("CC(C)(C)O"), 320.0, f));
  CHECK(!liquid_viscosity(pc("CCc1ccccc1O"), 320.0, f));
  CHECK_EQ(f.count(DiagCode::NoPropertyData), std::size_t{2});
}

TEST(liquid_heat_capacity) {
  // Water against IAPWS-95 (75.315 J/mol/K at 300 K, 75.564 at 350 K), and the enthalpy integral.
  const PureComponent w = pc("O");
  Diagnostics d;
  CHECK_NEAR(*liquid_cp(w, 300.0, d), 75.315, 0.004 * 75.315);
  CHECK_NEAR(*liquid_cp(w, 350.0, d), 75.564, 0.004 * 75.564);
  const double h = 1e-3;
  CHECK_NEAR((*liquid_enthalpy_change(w, 300.0, 330.0 + h, d) - *liquid_enthalpy_change(w, 300.0, 330.0 - h, d)) / (2 * h),
             *liquid_cp(w, 330.0, d), 1e-6);
  CHECK_NEAR(*liquid_enthalpy_change(w, 330.0, 300.0, d), -*liquid_enthalpy_change(w, 300.0, 330.0, d), 1e-9);
  CHECK(!d.has_errors());
  // Rows that agree with neither Poling's nor CRC's Cp(298) are dropped (acetic anhydride).
  Diagnostics e;
  CHECK(!liquid_cp(pc("CC(=O)OC(C)=O"), 300.0, e));
  CHECK(e.contains(DiagCode::NoPropertyData));
}
