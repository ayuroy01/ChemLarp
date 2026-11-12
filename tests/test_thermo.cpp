#include "chemlarp/joback.hpp"
#include "chemlarp/reaction.hpp"
#include "chemlarp/smiles.hpp"
#include "chemlarp/thermo.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

// ΔfH°(298) in kJ/mol and S°(298) in J/mol/K from the built-in database.
void check_species(const char* smiles, double dHf_kJ, double S_J, double tolH, double tolS) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smiles, d);
  CHECK(sp >= 0);
  auto t = resolve_thermo(st, static_cast<std::uint32_t>(sp), kTref, d);
  CHECK(t.has_value());
  CHECK(t->source == ThermoSource::Nasa7Database);
  CHECK_NEAR(t->h_RT(kTref) * kGasConstant * kTref / 1000.0, dHf_kJ, tolH);
  CHECK_NEAR(t->s_R(kTref) * kGasConstant, S_J, tolS);
}

double joback_hf(const char* smiles, bool* ok = nullptr) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smiles, d);
  CHECK(sp >= 0);
  auto r = joback_estimate(st, static_cast<std::uint32_t>(sp), d);
  if (ok) *ok = r.ok;
  return r.dHf_kJ;
}

}  // namespace

// Reference values: NIST-JANAF / CODATA.
TEST(nasa7_matches_reference_at_298) {
  check_species("O", -241.83, 188.84, 0.5, 0.5);         // H2O
  check_species("O=C=O", -393.51, 213.79, 0.5, 0.5);     // CO2
  check_species("C", -74.87, 186.25, 0.6, 0.5);          // CH4
  check_species("[C-]#[O+]", -110.53, 197.66, 0.5, 0.5);  // CO
  check_species("[H][H]", 0.0, 130.68, 0.1, 0.3);
  check_species("O=O", 0.0, 205.15, 0.1, 0.3);
  check_species("N#N", 0.0, 191.61, 0.1, 0.3);
  check_species("[H]", 217.998, 114.72, 0.5, 0.5);
  check_species("[O]", 249.18, 161.06, 0.5, 0.5);
}

TEST(nasa7_continuous_at_tmid) {
  auto t = builtin_thermo("H2O");
  CHECK(t.has_value());
  const double lo = 999.999999, hi = 1000.0;
  CHECK_NEAR(t->nasa.h_RT(lo), t->nasa.h_RT(hi), 1e-4);
  CHECK_NEAR(t->nasa.s_R(lo), t->nasa.s_R(hi), 1e-4);
}

TEST(joback_against_experiment) {
  CHECK_NEAR(joback_hf("CCO"), -234.8, 5.0);      // ethanol
  CHECK_NEAR(joback_hf("CCC"), -104.7, 3.0);      // propane
  CHECK_NEAR(joback_hf("CC(C)=O"), -217.1, 3.0);  // acetone
  CHECK_NEAR(joback_hf("c1ccccc1"), 82.9, 5.0);   // benzene
  CHECK_NEAR(joback_hf("CC(=O)O"), -432.8, 15.0); // acetic acid
}

TEST(joback_reports_unassignable_atoms) {
  bool ok = true;
  joback_hf("C", &ok);  // methane has no Joback group
  CHECK(!ok);
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, "C[Si](C)(C)C", d);
  CHECK(sp >= 0);
  auto r = joback_estimate(st, static_cast<std::uint32_t>(sp), d);
  CHECK(!r.ok);
  CHECK(d.contains(DiagCode::NoJobackGroup));
}

TEST(joback_fallback_is_flagged) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, "CCO", d);
  auto t = resolve_thermo(st, static_cast<std::uint32_t>(sp), 298.15, d);
  CHECK(t.has_value());
  CHECK(t->source == ThermoSource::JobackEstimate);
  CHECK(d.contains(DiagCode::ThermoApproximation));
  // Ethanol S°(gas, 298) is 281.6 J/mol/K; Joback-derived value should be in range.
  CHECK_NEAR(t->s298_J, 281.6, 15.0);
}

// ---- stoichiometry error vectors -------------------------------------------

TEST(reaction_balanced_combustion) {
  ChemicalSystemState st;
  Diagnostics d;
  auto rxn = parse_reaction(st, "1 C + 2 O=O -> 1 O=C=O + 2 O", d);
  CHECK(rxn.has_value());
  CHECK(check_balance(st, *rxn, d));
  CHECK(!d.has_errors());
  auto th = reaction_thermo(st, *rxn, kTref, d);
  CHECK(th.has_value());
  CHECK_NEAR(th->dH_J / 1000.0, -802.3, 1.5);  // CH4 + 2 O2 -> CO2 + 2 H2O(g)
}

TEST(reaction_imbalance_vector) {
  ChemicalSystemState st;
  Diagnostics d;
  auto rxn = parse_reaction(st, "C + O=O -> O=C=O + O", d);
  CHECK(rxn.has_value());
  CHECK(!check_balance(st, *rxn, d));
  CHECK_EQ(d.count(DiagCode::ElementImbalance), 2u);
  bool sawH = false, sawO = false;
  for (const auto& x : d.items()) {
    if (x.code != DiagCode::ElementImbalance) continue;
    if (x.element == 1) { sawH = true; CHECK_NEAR(x.delta, -2.0, 1e-12); }
    if (x.element == 8) { sawO = true; CHECK_NEAR(x.delta, 1.0, 1e-12); }
  }
  CHECK(sawH && sawO);
}

TEST(reaction_charge_imbalance) {
  ChemicalSystemState st;
  Diagnostics d;
  auto rxn = parse_reaction(st, "[NH4+] -> N", d);
  CHECK(rxn.has_value());
  CHECK(!check_balance(st, *rxn, d));
  CHECK(d.contains(DiagCode::ChargeImbalance));
  CHECK(d.contains(DiagCode::ElementImbalance));
  ChemicalSystemState st2;
  Diagnostics d2;
  auto ok = parse_reaction(st2, "[NH4+] -> N + [H+]", d2);
  CHECK(ok.has_value());
  CHECK(check_balance(st2, *ok, d2));
}

TEST(reaction_malformed_inputs) {
  auto errs = [](const char* eq) {
    ChemicalSystemState st;
    Diagnostics d;
    parse_reaction(st, eq, d);
    return d;
  };
  CHECK(errs("0 C -> C").contains(DiagCode::InvalidCoefficient));
  CHECK(errs("C + O=O").contains(DiagCode::MalformedReaction));
  CHECK(errs("-> C").contains(DiagCode::EmptyReactionSide));
  CHECK(errs("C O -> C").contains(DiagCode::MalformedReaction));
  CHECK(errs("C + -> C").contains(DiagCode::MalformedReaction));
  auto bad = errs("C(C)(C)(C)(C)C -> C");  // valence error inside a term, position is equation-relative
  CHECK(bad.contains(DiagCode::ValenceExceeded));
  auto pos = errs("2 O=O -> 2 O(C)(C)C");
  CHECK(pos.contains(DiagCode::ValenceExceeded));
  CHECK_EQ(pos.items()[0].position, 11u);
}
