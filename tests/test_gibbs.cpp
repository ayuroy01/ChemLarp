#include <cmath>
#include <utility>

#include "omnichem/gibbs.hpp"
#include "omnichem/smiles.hpp"
#include "omnichem/thermo.hpp"
#include "test_framework.hpp"

using namespace omnichem;

namespace {

struct System {
  ChemicalSystemState st;
  Diagnostics d;
  std::vector<std::int32_t> idx;
  System(double T, double P_bar, std::initializer_list<std::pair<const char*, double>> species) {
    st.conditions = {T, P_bar * 1e5};
    for (auto& [smi, n] : species) {
      const std::int32_t sp = parse_smiles(st, smi, d);
      CHECK(sp >= 0);
      st.moles[static_cast<std::size_t>(sp)] = n;
      idx.push_back(sp);
    }
  }
  double n(std::size_t i) const { return st.moles[static_cast<std::size_t>(idx[i])]; }
};

double element_total(const ChemicalSystemState& st, std::size_t z) {
  double s = 0;
  for (std::size_t i = 0; i < st.species.size(); ++i) s += st.species[i].composition[z] * st.moles[i];
  return s;
}

}  // namespace

TEST(solve_dense_basic) {
  std::vector<double> m = {2, 1, 1, 3};
  std::vector<double> r = {3, 5};
  CHECK(solve_dense(m, r, 2));
  CHECK_NEAR(r[0], 0.8, 1e-14);
  CHECK_NEAR(r[1], 1.4, 1e-14);
  std::vector<double> sing = {1, 2, 2, 4};
  std::vector<double> rs = {1, 2};
  CHECK(!solve_dense(sing, rs, 2));
}

// CO + H2O <=> CO2 + H2 : equal moles of gas on both sides, so x-ratio = K exactly.
TEST(water_gas_shift_matches_equilibrium_constant) {
  System s(1000.0, 1.0, {{"[C-]#[O+]", 1.0}, {"O", 1.0}, {"O=C=O", 0.0}, {"[H][H]", 0.0}});
  auto r = equilibrate(s.st, s.d);
  CHECK(r.converged);
  CHECK(!s.d.has_errors());
  const auto& g = s.st.g_RT;
  const double K = std::exp(-(g[2] + g[3] - g[0] - g[1]));
  CHECK(K > 1.0 && K < 2.0);  // literature: K(1000 K) ≈ 1.4
  const double Q = (s.n(2) * s.n(3)) / (s.n(0) * s.n(1));
  CHECK_NEAR(Q / K, 1.0, 1e-8);
  CHECK(r.max_relative_element_error < 1e-10);
}

TEST(methane_combustion_conserves_elements) {
  System s(2000.0, 1.0,
           {{"C", 1.0}, {"O=O", 2.0}, {"O=C=O", 0.0}, {"O", 0.0}, {"[C-]#[O+]", 0.0}, {"[H][H]", 0.0},
            {"[OH]", 0.0}, {"[H]", 0.0}, {"[O]", 0.0}});
  auto r = equilibrate(s.st, s.d);
  CHECK(r.converged);
  CHECK_NEAR(element_total(s.st, 6), 1.0, 1e-10);
  CHECK_NEAR(element_total(s.st, 1), 4.0, 1e-10);
  CHECK_NEAR(element_total(s.st, 8), 4.0, 1e-10);
  CHECK(s.n(0) < 1e-12);     // methane is fully consumed
  CHECK(s.n(2) > 0.9);       // mostly CO2
  CHECK(s.n(3) > 1.9);       // mostly H2O
  CHECK(s.n(4) > 1e-4);      // some CO from dissociation at 2000 K
}

TEST(high_temperature_shifts_dissociation) {
  auto co2_frac = [](double T) {
    System s(T, 1.0, {{"O=C=O", 1.0}, {"[C-]#[O+]", 0.0}, {"O=O", 0.0}});
    auto r = equilibrate(s.st, s.d);
    CHECK(r.converged);
    return s.n(0);
  };
  const double low = co2_frac(1500.0), high = co2_frac(3000.0);
  CHECK(low > 0.99);
  CHECK(high < low);
}

TEST(pressure_le_chatelier) {
  // 2 CO2 <=> 2 CO + O2 increases moles: higher pressure suppresses dissociation.
  auto co = [](double P) {
    System s(3000.0, P, {{"O=C=O", 1.0}, {"[C-]#[O+]", 0.0}, {"O=O", 0.0}});
    CHECK(equilibrate(s.st, s.d).converged);
    return s.n(1);
  };
  CHECK(co(10.0) < co(1.0));
}

TEST(rank_deficient_element_matrix) {
  // Only water present: H row = 2 x O row. The solver must drop the dependent row.
  System s(500.0, 1.0, {{"O", 2.0}});
  auto r = equilibrate(s.st, s.d);
  CHECK(r.converged);
  CHECK_EQ(r.elements.size(), 1u);
  CHECK_NEAR(s.n(0), 2.0, 1e-10);
}

TEST(species_without_available_elements_are_zero) {
  System s(1000.0, 1.0, {{"[H][H]", 1.0}, {"N#N", 0.0}, {"O", 0.0}});
  auto r = equilibrate(s.st, s.d);
  CHECK(r.converged);
  CHECK_EQ(s.n(1), 0.0);
  CHECK_EQ(s.n(2), 0.0);
  CHECK_NEAR(s.n(0), 1.0, 1e-12);
}

TEST(equilibrium_error_vectors) {
  {
    System s(1000.0, 1.0, {{"[NH4+]", 1.0}});
    equilibrate(s.st, s.d);
    CHECK(s.d.contains(DiagCode::ChargedSpeciesUnsupported));
  }
  {
    System s(1000.0, 1.0, {{"O", -1.0}});
    equilibrate(s.st, s.d);
    CHECK(s.d.contains(DiagCode::NegativeAmount));
  }
  {
    System s(1000.0, 1.0, {{"O", 0.0}});
    equilibrate(s.st, s.d);
    CHECK(s.d.contains(DiagCode::EmptySystem));
  }
  {
    System s(1000.0, 1.0, {{"C[Si](C)(C)C", 1.0}});
    equilibrate(s.st, s.d);
    CHECK(s.d.contains(DiagCode::NoThermoData));
  }
}
