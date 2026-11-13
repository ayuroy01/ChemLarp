// Activity-coefficient models and modified-Raoult VLE.
//
// Reference values:
//  * ln γ for NRTL and UNIFAC from thermo 0.6.1 (Caleb Bell, MIT licence), an independent implementation;
//    the scratch comparison over 300 random NRTL and 400 random UNIFAC systems agreed to 5e-14.
//    thermo's NRTL uses α = c + d T; ChemLarp uses α = c + d (T − 273.15 K) (Aspen convention), so the
//    reference was computed with c' = c − 273.15 d.
//  * Normal boiling points: CRC Handbook of Chemistry and Physics.
//  * Azeotropes at 101.325 kPa: CRC Handbook (44th ed.) / Lange's Handbook (10th ed.), via the compilation in
//    the Wikipedia "Azeotrope tables" article (composition as weight %). Tolerances are documented in
//    docs/ARCHITECTURE.md and reflect UNIFAC's predictive accuracy, not the numerical tolerance.
#include <cmath>
#include <string>
#include <vector>

#include "chemlarp/activity.hpp"
#include "chemlarp/chemlarp_c.h"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

ActivitySystem unifac(const std::vector<std::string>& smiles) {
  Diagnostics d;
  auto s = unifac_system_from_smiles(smiles, d);
  CHECK(s.has_value());
  return *s;
}

NrtlParameters nrtl_square(std::size_t n) {
  NrtlParameters p;
  p.n = n;
  p.a.assign(n * n, 0.0);
  p.b.assign(n * n, 0.0);
  p.alpha_c.assign(n * n, 0.3);
  return p;
}

// Liquid mole fraction of B from weight % of B.
double mole_fraction(double wB, double Ma, double Mb) { return (wB / Mb) / (wB / Mb + (100.0 - wB) / Ma); }

struct Azeotrope {
  double T_K = 0.0, xB = 0.0;
  bool found = false;
};

// Locate a binary azeotrope at P by the sign change of y_A − x_A along the bubble-point curve.
Azeotrope find_azeotrope(const ActivitySystem& sys, double P) {
  auto diff = [&](double xA, double* T) {
    Diagnostics d;
    auto r = bubble_temperature(sys, P, {xA, 1.0 - xA}, d);
    CHECK(r.has_value());
    if (T) *T = r->T_K;
    return r->y[0] - xA;
  };
  Azeotrope az;
  double prev_x = 0.01, prev = diff(prev_x, nullptr);
  for (int i = 2; i < 100; ++i) {
    const double x = i / 100.0, f = diff(x, nullptr);
    if (prev * f <= 0.0) {
      double lo = prev_x, hi = x, flo = prev, T = 0.0;
      for (int k = 0; k < 50; ++k) {
        const double m = 0.5 * (lo + hi), fm = diff(m, &T);
        if (flo * fm <= 0.0) hi = m; else { lo = m; flo = fm; }
      }
      az.found = true;
      az.xB = 1.0 - 0.5 * (lo + hi);
      az.T_K = T;
      return az;
    }
    prev_x = x;
    prev = f;
  }
  return az;
}

}  // namespace

TEST(nrtl_matches_reference_implementation) {
  Diagnostics d;
  NrtlParameters p = nrtl_square(2);
  p.a = {0, 0.5, -0.3, 0};
  p.b = {0, -120.0, 650.0, 0};
  p.alpha_c = {0, 0.3, 0.3, 0};
  p.alpha_d = {0, 0.001, 0.001, 0};
  auto g = nrtl_ln_gamma(p, 343.15, {0.25, 0.75}, d);
  CHECK(g.has_value());
  CHECK_NEAR((*g)[0], 0.7039497258148246, 1e-12);
  CHECK_NEAR((*g)[1], 0.1332494307590076, 1e-12);

  NrtlParameters q;
  q.n = 3;
  q.a = {0, 1.1, -0.4, 0.2, 0, 0.9, 0.7, -0.6, 0};
  q.b = {0, 100, -50, 200, 0, 30, -80, 150, 0};
  q.e = {0, 0.01, 0, 0, 0, 0.02, 0.03, 0, 0};
  q.f = {0, 1e-4, 0, 0, 0, -2e-4, 0, 3e-4, 0};
  q.alpha_c = {0, 0.2, 0.3, 0.2, 0, 0.47, 0.3, 0.47, 0};
  auto h = nrtl_ln_gamma(q, 300.0, {0.2, 0.3, 0.5}, d);
  CHECK(h.has_value());
  CHECK_NEAR((*h)[0], 0.3114055742151451, 1e-12);
  CHECK_NEAR((*h)[1], 0.5489818517196942, 1e-12);
  CHECK_NEAR((*h)[2], -0.035399147985800225, 1e-12);
  CHECK(!d.has_errors());
}

TEST(nrtl_limits) {
  Diagnostics d;
  // τ = 0 everywhere: ideal solution.
  auto g = nrtl_ln_gamma(nrtl_square(3), 320.0, {0.2, 0.5, 0.3}, d);
  CHECK(g.has_value());
  for (double v : *g) CHECK_NEAR(v, 0.0, 1e-15);
  // Pure component: γ = 1; infinitely dilute component: ln γ∞ = τ_21 + τ_12 exp(−α τ_12).
  NrtlParameters p = nrtl_square(2);
  p.a = {0, 0.8, 1.5, 0};
  auto pure = nrtl_ln_gamma(p, 300.0, {1.0, 0.0}, d);
  CHECK(pure.has_value());
  CHECK_NEAR((*pure)[0], 0.0, 1e-15);
  CHECK_NEAR((*pure)[1], 0.8 + 1.5 * std::exp(-0.3 * 1.5), 1e-14);
  // Mole fractions are normalised.
  auto a = nrtl_ln_gamma(p, 300.0, {0.2, 0.6}, d);
  auto b = nrtl_ln_gamma(p, 300.0, {0.25, 0.75}, d);
  CHECK_NEAR((*a)[0], (*b)[0], 1e-15);
  CHECK(!d.has_errors());
}

TEST(unifac_matches_reference_implementation) {
  Diagnostics d;
  auto g = unifac_ln_gamma({{{1, 1}, {2, 1}, {14, 1}}, {{16, 1}}}, 323.15, {0.3, 0.7}, d);
  CHECK(g.has_value());
  CHECK_NEAR((*g)[0], 0.49854592321978186, 1e-12);
  CHECK_NEAR((*g)[1], 0.20649363596105091, 1e-12);
  // acetone + chloroform + methanol
  auto h = unifac_ln_gamma({{{1, 1}, {18, 1}}, {{50, 1}}, {{15, 1}}}, 330.0, {0.4, 0.35, 0.25}, d);
  CHECK(h.has_value());
  CHECK_NEAR((*h)[0], -0.17307596699774236, 1e-12);
  CHECK_NEAR((*h)[1], -0.05471861375714562, 1e-12);
  CHECK_NEAR((*h)[2], 0.5418851764461461, 1e-12);
  // toluene infinitely dilute in hexane
  auto t = unifac_ln_gamma({{{9, 5}, {11, 1}}, {{1, 2}, {2, 4}}}, 298.15, {0.0, 1.0}, d);
  CHECK(t.has_value());
  CHECK_NEAR((*t)[0], 0.48540061456612177, 1e-12);
  CHECK_NEAR((*t)[1], 0.0, 1e-14);
  CHECK(!d.has_errors());
}

TEST(unifac_limits_and_gibbs_duhem) {
  Diagnostics d;
  const std::vector<UnifacGroups> etoh_water = {{{1, 1}, {2, 1}, {14, 1}}, {{16, 1}}};
  // Pure components have γ = 1; identical components are an ideal mixture.
  for (double x1 : {0.0, 1.0}) {
    auto g = unifac_ln_gamma(etoh_water, 330.0, {x1, 1.0 - x1}, d);
    CHECK(g.has_value());
    CHECK_NEAR(x1 == 1.0 ? (*g)[0] : (*g)[1], 0.0, 1e-14);
  }
  auto same = unifac_ln_gamma({{{1, 1}, {2, 1}, {14, 1}}, {{1, 1}, {2, 1}, {14, 1}}}, 330.0, {0.4, 0.6}, d);
  CHECK(same.has_value());
  CHECK_NEAR((*same)[0], 0.0, 1e-14);
  CHECK_NEAR((*same)[1], 0.0, 1e-14);
  // Gibbs–Duhem at constant T: x1 d ln γ1/dx1 + x2 d ln γ2/dx1 = 0.
  for (double x1 : {0.1, 0.35, 0.6, 0.9}) {
    const double h = 1e-6;
    auto p = unifac_ln_gamma(etoh_water, 330.0, {x1 + h, 1.0 - x1 - h}, d);
    auto m = unifac_ln_gamma(etoh_water, 330.0, {x1 - h, 1.0 - x1 + h}, d);
    const double d1 = ((*p)[0] - (*m)[0]) / (2 * h), d2 = ((*p)[1] - (*m)[1]) / (2 * h);
    CHECK_NEAR(x1 * d1 + (1.0 - x1) * d2, 0.0, 1e-7);
  }
  CHECK(!d.has_errors());
}

TEST(activity_invalid_input_is_reported) {
  const std::vector<UnifacGroups> etoh_water = {{{1, 1}, {2, 1}, {14, 1}}, {{16, 1}}};
  {
    Diagnostics d;
    CHECK(!unifac_ln_gamma(etoh_water, 330.0, {0.5}, d));
    CHECK(d.contains(DiagCode::InvalidMixtureInput));
  }
  {
    Diagnostics d;
    CHECK(!unifac_ln_gamma(etoh_water, 330.0, {-0.1, 1.1}, d));
    CHECK(!unifac_ln_gamma(etoh_water, 330.0, {0.0, 0.0}, d));
    CHECK(!unifac_ln_gamma(etoh_water, -5.0, {0.5, 0.5}, d));
    CHECK(!unifac_ln_gamma({{{999, 1}}, {{16, 1}}}, 330.0, {0.5, 0.5}, d));
    CHECK(!nrtl_ln_gamma(NrtlParameters{2, {}, {}, {}, {}, {}, {}}, 300.0, {0.5, 0.5}, d));
    CHECK_EQ(d.count(DiagCode::InvalidMixtureInput), std::size_t{5});
  }
  {
    // Water + iodomethane: no published a_mn between main groups 7 (H2O) and 32 (I); never taken as zero.
    Diagnostics d;
    CHECK(!unifac_ln_gamma({{{16, 1}}, {{1, 1}, {63, 1}}}, 300.0, {0.5, 0.5}, d));
    CHECK(d.contains(DiagCode::NoInteractionParameter));
  }
}

TEST(unifac_groups_from_structure) {
  auto groups = [](const char* smiles) {
    Diagnostics d;
    auto g = unifac_groups_from_smiles(smiles, d);
    return g ? *g : UnifacGroups{};
  };
  CHECK(groups("CCO") == (UnifacGroups{{1, 1}, {2, 1}, {14, 1}}));
  CHECK(groups("O") == (UnifacGroups{{16, 1}}));
  CHECK(groups("Cc1ccccc1") == (UnifacGroups{{9, 5}, {11, 1}}));
  CHECK(groups("CC(C)=O") == (UnifacGroups{{1, 1}, {18, 1}}));
  CHECK(groups("CCOC(C)=O") == (UnifacGroups{{1, 1}, {2, 1}, {21, 1}}));  // CH3COO + CH2 + CH3 (DDBST)
  CHECK(groups("C1CCOC1") == (UnifacGroups{{2, 3}, {27, 1}}));
  Diagnostics d;
  CHECK(!unifac_groups_from_smiles("C[Si](C)(C)C", d));
  CHECK(d.contains(DiagCode::NoUnifacGroups));
}

TEST(antoine_normal_boiling_points) {
  // CRC normal boiling points (°C); the Poling et al. Antoine fits reproduce them within 0.5 K.
  const std::vector<std::pair<const char*, double>> tb = {
      {"O", 100.0},         {"CCO", 78.4},       {"CO", 64.7},           {"CC(C)=O", 56.5},
      {"ClC(Cl)Cl", 61.1},  {"c1ccccc1", 80.2},  {"C1CCCCC1", 80.7},     {"CCCCCC", 68.9},
      {"Cc1ccccc1", 110.8}, {"CCCO", 97.3},      {"ClC(Cl)(Cl)Cl", 76.8}};
  for (const auto& [smiles, t] : tb) {
    Diagnostics d;
    auto r = bubble_temperature(unifac({smiles}), 101325.0, {1.0}, d);
    CHECK(r.has_value());
    CHECK_NEAR(r->T_K - 273.15, t, 0.5);
    CHECK(!r->extrapolated);
  }
}

TEST(vle_bubble_dew_consistency) {
  const ActivitySystem sys = unifac({"CCO", "O"});
  Diagnostics d;
  auto bub = bubble_temperature(sys, 101325.0, {0.3, 0.7}, d);
  CHECK(bub.has_value());
  // The dew point of the bubble vapour is the bubble point, with the original liquid.
  auto dew = dew_temperature(sys, 101325.0, bub->y, d);
  CHECK(dew.has_value());
  CHECK_NEAR(dew->T_K, bub->T_K, 1e-6);
  CHECK_NEAR(dew->x[0], 0.3, 1e-6);
  // Bubble pressure at the bubble temperature returns the pressure; dew P ≤ bubble P at fixed T and z.
  auto bp = bubble_pressure(sys, bub->T_K, {0.3, 0.7}, d);
  CHECK(bp.has_value());
  CHECK_NEAR(bp->P_Pa, 101325.0, 1e-3);
  auto dp = dew_pressure(sys, bub->T_K, {0.3, 0.7}, d);
  CHECK(dp.has_value());
  CHECK(dp->P_Pa < bp->P_Pa);
  double sy = 0.0;
  for (double v : bub->y) sy += v;
  CHECK_NEAR(sy, 1.0, 1e-9);
  CHECK(!d.has_errors());
}

TEST(vle_reproduces_experimental_azeotropes) {
  // {A, B, M_A, M_B, T_az (°C), wt % B}; UNIFAC predictions agree within 0.6 K and 0.05 in x_B
  // (see ARCHITECTURE.md for the full table).
  struct Row { const char *a, *b; double Ma, Mb, T, wB; };
  const std::vector<Row> rows = {
      {"O", "CCO", 18.015, 46.069, 78.1, 95.5},                 // minimum-boiling
      {"CO", "c1ccccc1", 32.042, 78.114, 58.3, 60.5},
      {"CC(C)=O", "ClC(Cl)Cl", 58.080, 119.37, 64.7, 80.0},     // maximum-boiling
      {"c1ccccc1", "C1CCCCC1", 78.114, 84.162, 77.8, 45.0},
      {"O", "CC(C)O", 18.015, 60.096, 80.4, 87.7},
  };
  for (const Row& r : rows) {
    const Azeotrope az = find_azeotrope(unifac({r.a, r.b}), 101325.0);
    CHECK(az.found);
    CHECK_NEAR(az.T_K - 273.15, r.T, 0.6);
    CHECK_NEAR(az.xB, mole_fraction(r.wB, r.Ma, r.Mb), 0.05);
  }
  // Methanol + water has no azeotrope (CRC / Lange).
  CHECK(!find_azeotrope(unifac({"CO", "O"}), 101325.0).found);
}

TEST(vle_extrapolated_vapour_pressure_is_flagged) {
  Diagnostics d;
  auto a = antoine_coefficients("O");
  CHECK(a.has_value());
  auto r = bubble_pressure(unifac({"O"}), a->Tmax_K + 20.0, {1.0}, d);
  CHECK(r.has_value());
  CHECK(r->extrapolated);
  CHECK(d.contains(DiagCode::NoVaporPressureData));
  CHECK(!d.has_errors());
  // A compound with UNIFAC groups but no Antoine coefficients is refused, not guessed.
  Diagnostics e;
  CHECK(!unifac_system_from_smiles({"CCO", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"}, e));
  CHECK(e.contains(DiagCode::NoVaporPressureData));
}

TEST(vle_c_api) {
  auto call = [](const char* lines, const char* mode, double v) {
    char* p = cl_vle_text_json(lines, mode, v);
    CHECK(p != nullptr);
    std::string s(p);
    cl_free_string(p);
    return s;
  };
  const std::string ok = call("CCO 0.3\nO 0.7", "bubble-T", 101325.0);
  CHECK(ok.rfind("{\"ok\":true,\"T_K\":354.8", 0) == 0);
  CHECK(ok.find("\"name\":\"water\"") != std::string::npos);
  CHECK(call("CCO 1", "triple-point", 300.0).find("\"ok\":false") != std::string::npos);
  CHECK(call("CCO 0.5\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 0.5", "bubble-P", 300.0).find("NoVaporPressureData") != std::string::npos);
}

TEST(vle_ideal_nrtl_is_raoults_law) {
  // NRTL with τ = 0 is an ideal solution: bubble P = Σ x_i P_i^sat, dew 1/P = Σ y_i / P_i^sat.
  ActivitySystem sys;
  sys.model = ActivityModel::Nrtl;
  sys.nrtl = nrtl_square(2);
  sys.antoine = {*antoine_coefficients("c1ccccc1"), *antoine_coefficients("Cc1ccccc1")};
  const double T = 360.0, p1 = sys.antoine[0].pressure(T), p2 = sys.antoine[1].pressure(T);
  Diagnostics d;
  auto b = bubble_pressure(sys, T, {0.4, 0.6}, d);
  CHECK(b.has_value());
  CHECK_NEAR(b->P_Pa / (0.4 * p1 + 0.6 * p2), 1.0, 1e-12);
  CHECK_NEAR(b->y[0], 0.4 * p1 / b->P_Pa, 1e-12);
  auto w = dew_pressure(sys, T, {0.4, 0.6}, d);
  CHECK(w.has_value());
  CHECK_NEAR(w->P_Pa * (0.4 / p1 + 0.6 / p2), 1.0, 1e-9);
  CHECK(!d.has_errors());
}
