#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../src/yaml_lite.hpp"
#include "omnichem/kinetics.hpp"
#include "omnichem/thermo.hpp"
#include "test_framework.hpp"

using namespace omnichem;

namespace {

Mechanism mech(const std::string& text) {
  Diagnostics d;
  auto m = parse_mechanism(text, d);
  if (!m) throw octest::Failure("mechanism failed to parse: " + d.items().front().message);
  return *m;
}

Diagnostics errors_of(const std::string& text) {
  Diagnostics d;
  parse_mechanism(text, d);
  return d;
}

std::string read(const std::string& relative) {
  std::ifstream in(octest::data_path(relative), std::ios::binary);
  if (!in) throw octest::Failure("cannot read " + relative);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

const std::string kHO2 = "species\n H [H]\n O2 O=O\n HO2 [O]O\n AR [Ar]\nend\n";

}  // namespace

TEST(falloff_lindemann_limits) {
  // k_inf = 2, k_0 = 4 (SI, m3/J): Pr = 2 [M].
  auto m = mech(kHO2 + "reactions J/mol m3\n H + O2 (+M) <=> HO2 (+M)  2 0 0\n  LOW / 4 0 0 /\nend\n");
  const auto& r = m.reactions.at(0);
  CHECK(r.kind == RateKind::FallOff);
  CHECK_NEAR(falloff_rate_constant(r, 1000, 1e-9) / (4 * 1e-9), 1.0, 1e-8);  // low-pressure limit k0 [M]
  CHECK_NEAR(falloff_rate_constant(r, 1000, 1e9) / 2.0, 1.0, 1e-8);          // high-pressure limit k_inf
  CHECK_NEAR(falloff_rate_constant(r, 1000, 0.5), 1.0, 1e-12);               // Pr = 1: k_inf / 2
}

TEST(falloff_troe_reference_values) {
  // a = 0.5, T3 = 1e30, T1 = 1e-30 => Fcent = 0.5 at any T. Reference k/k_inf from an
  // independent evaluation of the Troe formula (Gilbert, Luther & Troe 1983).
  auto m = mech(kHO2 + "reactions J/mol m3\n H + O2 (+M) <=> HO2 (+M)  2 0 0\n  LOW / 4 0 0 /\n  TROE / 0.5 1e30 1e-30 /\nend\n");
  const auto& r = m.reactions.at(0);
  CHECK_NEAR(r.troe->F_cent(1234.0), 0.5, 1e-15);
  const std::vector<std::pair<double, double>> ref = {
      {1.0, 0.2549688514647912}, {0.01, 0.008040574224104261}, {100.0, 0.8663404417220719}};
  for (auto [Pr, k_over_kinf] : ref) CHECK_NEAR(falloff_rate_constant(r, 1000, Pr / 2.0) / 2.0, k_over_kinf, 1e-13);
  // dk/dM analytic vs central difference
  double dk = 0;
  falloff_rate_constant(r, 1000, 0.7, &dk);
  const double h = 1e-6;
  const double fd = (falloff_rate_constant(r, 1000, 0.7 + h) - falloff_rate_constant(r, 1000, 0.7 - h)) / (2 * h);
  CHECK_NEAR(dk / fd, 1.0, 1e-7);
}

TEST(chemically_activated_limits) {
  // HIGH: main line is k_0 = 4, HIGH is k_inf = 2.  k = k0 / (1 + Pr) F.
  auto m = mech(kHO2 + "reactions J/mol m3\n H + O2 (+M) <=> HO2 (+M)  4 0 0\n  HIGH / 2 0 0 /\nend\n");
  const auto& r = m.reactions.at(0);
  CHECK(r.kind == RateKind::ChemicallyActivated);
  CHECK_NEAR(falloff_rate_constant(r, 1000, 1e-9) / 4.0, 1.0, 1e-8);                 // low P: k0
  CHECK_NEAR(falloff_rate_constant(r, 1000, 1e6) * (4.0 * 1e6 / 2.0) / 4.0, 1.0, 1e-5);  // high P: k0 / Pr
}

TEST(specific_collider_uses_only_that_species) {
  auto m = mech(kHO2 + "reactions J/mol m3\n H + O2 (+AR) <=> HO2 (+AR)  2 0 0\n  LOW / 4 0 0 /\nend\n");
  const auto& r = m.reactions.at(0);
  CHECK_EQ(r.efficiency[static_cast<std::size_t>(m.find("AR"))], 1.0);
  CHECK_EQ(r.efficiency[static_cast<std::size_t>(m.find("O2"))], 0.0);
}

TEST(plog_interpolation_and_clamping) {
  auto m = mech(kHO2 + "reactions J/mol m3\n H + O2 <=> HO2  1 0 0\n"
                       "  PLOG / 1.0  10 0 0 /\n  PLOG / 10.0 1000 0 0 /\n"
                       "  PLOG / 10.0 1000 0 0 /\nend\n");  // two entries at one pressure are summed
  const auto& r = m.reactions.at(0);
  CHECK(r.kind == RateKind::PressureLog);
  CHECK_EQ(r.plog.size(), 2u);
  const double atm = 101325.0;
  CHECK_NEAR(plog_rate_constant(r, 1000, 1 * atm), 10.0, 1e-9);
  CHECK_NEAR(plog_rate_constant(r, 1000, 10 * atm), 2000.0, 1e-7);
  CHECK_NEAR(plog_rate_constant(r, 1000, std::sqrt(10.0) * atm), std::sqrt(10.0 * 2000.0), 1e-7);  // log-log midpoint
  CHECK_NEAR(plog_rate_constant(r, 1000, 0.01 * atm), 10.0, 1e-9);   // clamped below
  CHECK_NEAR(plog_rate_constant(r, 1000, 100 * atm), 2000.0, 1e-7);  // clamped above
}

TEST(explicit_reverse_rate) {
  auto m = mech("species\n A CCO\n B COC\nend\nreactions J/mol m3\n A <=> B  3 0 0\n  REV / 7 0 0 /\nend\n");
  std::vector<double> C = {2.0, 5.0}, q(1);
  rates_of_progress(m, 300, C.data(), q.data());
  CHECK_NEAR(q[0], 3 * 2.0 - 7 * 5.0, 1e-12);
}

TEST(jacobian_exact_for_all_rate_kinds) {
  Diagnostics d;
  auto m = load_mechanism(octest::data_path("mechanisms/h2o2.ocm"), d);
  CHECK(m.has_value());
  auto plog = mech(kHO2 + "reactions\n H + O2 <=> HO2  1e12 0 0\n  PLOG / 0.1 1e11 0 0 /\n  PLOG / 10 4e13 0.5 300 /\nend\n");
  for (const Mechanism* mm : {&*m, &plog}) {
    const std::size_t S = mm->species_count();
    std::vector<double> C(S);
    for (std::size_t i = 0; i < S; ++i) C[i] = 0.5 + 3.0 * static_cast<double>((i * 7) % 5);
    std::vector<double> J(S * S), wp(S), wm(S);
    production_rates_jacobian(*mm, 1500.0, C.data(), J.data());
    for (std::size_t j = 0; j < S; ++j) {
      auto cp = C, cm = C;
      const double h = 1e-4 * C[j];
      cp[j] += h;
      cm[j] -= h;
      production_rates(*mm, 1500.0, cp.data(), wp.data());
      production_rates(*mm, 1500.0, cm.data(), wm.data());
      for (std::size_t i = 0; i < S; ++i) {
        const double fd = (wp[i] - wm[i]) / (2 * h), a = J[i * S + j];
        if (std::fabs(fd - a) > 1e-6 * (std::fabs(fd) + std::fabs(a)) + 1e-8)
          throw octest::Failure("J[" + mm->names[i] + "][" + mm->names[j] + "] analytic " + std::to_string(a) + " fd " +
                                std::to_string(fd));
      }
    }
  }
}

// The same mechanism written as .ocm, CHEMKIN-II (+therm.dat, with and without a
// SMILES map) and Cantera YAML must import to identical kinetics.
TEST(three_formats_give_identical_rates) {
  Diagnostics d;
  auto ocm = load_mechanism(octest::data_path("mechanisms/h2o2.ocm"), d);
  auto ck = import_chemkin(read("mechanisms/h2o2/chem.inp"), read("mechanisms/h2o2/therm.dat"), "", d);
  auto ck_smi = import_chemkin(read("mechanisms/h2o2/chem.inp"), read("mechanisms/h2o2/therm.dat"),
                               read("mechanisms/h2o2/species.smi"), d);
  auto yml = load_mechanism(octest::data_path("mechanisms/h2o2.yaml"), d);
  if (!ocm || !ck || !ck_smi || !yml) throw octest::Failure("import failed: " + d.items().front().message);
  CHECK(!d.has_errors());
  CHECK(detect_mechanism_format(read("mechanisms/h2o2/chem.inp")) == MechanismFormat::Chemkin);
  CHECK(detect_mechanism_format(read("mechanisms/h2o2.yaml")) == MechanismFormat::CanteraYaml);
  CHECK(detect_mechanism_format(read("mechanisms/h2o2.ocm")) == MechanismFormat::OmniChem);
  const std::size_t S = ocm->species_count();
  CHECK_EQ(ocm->reactions.size(), 23u);
  for (const Mechanism* other : {&*ck, &*ck_smi, &*yml}) {
    CHECK_EQ(other->species_count(), S);
    CHECK_EQ(other->reactions.size(), ocm->reactions.size());
    for (std::size_t i = 0; i < S; ++i) CHECK_EQ(other->names[i], ocm->names[i]);
    for (double T : {600.0, 1200.0, 2400.0}) {
      std::vector<double> C(S), w0(S), w1(S);
      for (std::size_t i = 0; i < S; ++i) C[i] = 1.0 + static_cast<double>(i);
      production_rates(*ocm, T, C.data(), w0.data());
      production_rates(*other, T, C.data(), w1.data());
      for (std::size_t i = 0; i < S; ++i) CHECK_NEAR(w1[i], w0[i], 1e-10 * (std::fabs(w0[i]) + 1.0));
    }
  }
  // Thermo cards parsed from fixed columns reproduce the built-in GRI-Mech data.
  const auto h2o = builtin_thermo("H2O");
  const auto& card = *ck->thermo[static_cast<std::size_t>(ck->find("H2O"))];
  for (double T : {300.0, 999.0, 1001.0, 3000.0}) CHECK_NEAR(card.g_RT(T), h2o->g_RT(T), 1e-12);
  // Without a SMILES map, CHEMKIN species are composition-only; with it they have structures.
  CHECK(ck->state.species[0].atom_count == 0 && ck->state.species[0].formula == "H2");
  CHECK(ck_smi->state.species[0].canonical_smiles == "[H][H]");
}

TEST(yaml_reader_subset) {
  const std::string text =
      "# comment\n"
      "description: |-\n  multi\n  line # not a comment\n"
      "units: {length: cm, quantity: mol}\n"
      "list:\n- a\n- 'b: c'\n- {x: 1, y: [1, 2,\n    3]}\n"
      "nested:\n  key: \"quoted # value\"\n  seq:\n  - name: s1\n    v: 2\n  - name: s2\n";
  yaml::Node root;
  yaml::ParseError err;
  CHECK(yaml::parse(text, root, err));
  CHECK_EQ(root.get("units")->get("length")->scalar, std::string("cm"));
  CHECK(root.get("description")->scalar.find("line # not a comment") != std::string::npos);
  const auto& list = *root.get("list");
  CHECK_EQ(list.seq.size(), 3u);
  CHECK_EQ(list.seq[1].scalar, std::string("b: c"));
  CHECK_EQ(list.seq[2].get("y")->seq.size(), 3u);
  CHECK_EQ(root.get("nested")->get("key")->scalar, std::string("quoted # value"));
  CHECK_EQ(root.get("nested")->get("seq")->seq[1].get("name")->scalar, std::string("s2"));
  yaml::Node bad;
  CHECK(!yaml::parse("a: [1, 2\nb: 3\n", bad, err));
}

TEST(falloff_and_import_error_vectors) {
  const std::string r = "reactions J/mol m3\n";
  CHECK(errors_of(kHO2 + r + " H + O2 (+M) <=> HO2 (+M) 2 0 0\nend\n").contains(DiagCode::MechanismSyntax));  // no LOW
  CHECK(errors_of(kHO2 + r + " H + O2 (+M) <=> HO2 2 0 0\n  LOW/4 0 0/\nend\n").contains(DiagCode::MechanismSyntax));  // one side
  CHECK(errors_of(kHO2 + r + " H + O2 + M <=> HO2 + M 2 0 0\n  LOW/4 0 0/\nend\n").contains(DiagCode::MechanismSyntax));
  CHECK(errors_of(kHO2 + r + " H + O2 (+M) <=> HO2 (+M) 2 0 0\n  LOW/4 0 0/ HIGH/1 0 0/\nend\n").contains(DiagCode::MechanismSyntax));
  CHECK(errors_of(kHO2 + r + " H + O2 (+M) <=> HO2 (+M) 2 0 0\n  LOW/4 0 0/\n  TROE/0.5 1/\nend\n").contains(DiagCode::InvalidRateParameter));
  CHECK(errors_of(kHO2 + r + " H + O2 (+M) <=> HO2 (+M) 2 0 0\n  LOW/4 0 0/\n  SRI/1 2 3/\nend\n").contains(DiagCode::UnsupportedReactionType));
  CHECK(errors_of(kHO2 + r + " H + O2 (+M) <=> HO2 (+M) 2 0 0\n  LOW/4 0 0/\n  REV/1 0 0/\nend\n").contains(DiagCode::UnsupportedReactionType));
  CHECK(errors_of(kHO2 + r + " H + O2 => HO2 2 0 0\n  REV/1 0 0/\nend\n").contains(DiagCode::MechanismSyntax));
  CHECK(errors_of(kHO2 + r + " H + O2 + M <=> HO2 + M 2 0 0\n  PLOG/1 1 0 0/\nend\n").contains(DiagCode::MechanismSyntax));
  CHECK(errors_of(kHO2 + r + " H + O2 <=> HO2 2 0 0\n  H2/2/\nend\n").contains(DiagCode::UnknownSpecies));
  CHECK(errors_of(kHO2 + r + " H + O2 <=> HO2 -2 0 0\nend\n").contains(DiagCode::InvalidRateParameter));
  CHECK(errors_of(kHO2 + "reactions molecules\n H + O2 <=> HO2 2 0 0\nend\n").contains(DiagCode::UnsupportedReactionType));
  CHECK(errors_of(kHO2 + r + " LOW/4 0 0/\nend\n").contains(DiagCode::MechanismSyntax));  // aux before any reaction
  auto dup = errors_of(kHO2 + r + " H + O2 <=> HO2 2 0 0\n H + O2 <=> HO2 3 0 0\nend\n");
  CHECK(!dup.has_errors() && dup.contains(DiagCode::DuplicateReaction));
  auto ok_dup = errors_of(kHO2 + r + " H + O2 <=> HO2 2 0 0\n DUP\n H + O2 <=> HO2 -1 0 0\n DUP\nend\n");
  CHECK(!ok_dup.has_errors() && !ok_dup.contains(DiagCode::DuplicateReaction));

  // CHEMKIN: a species with neither thermo nor SMILES has no composition.
  Diagnostics d1;
  CHECK(!import_chemkin("ELEMENTS H END\nSPECIES H2 X END\nREACTIONS\nEND\n", read("mechanisms/h2o2/therm.dat"), "", d1));
  CHECK(d1.contains(DiagCode::UnknownSpecies));
  // A SMILES map entry that disagrees with the thermo card is reported.
  Diagnostics d2;
  CHECK(!import_chemkin(read("mechanisms/h2o2/chem.inp"), read("mechanisms/h2o2/therm.dat"), "H2O2 O\n", d2));
  CHECK(d2.contains(DiagCode::CompositionMismatch));
  // Malformed thermo card columns.
  Diagnostics d3;
  std::string therm = read("mechanisms/h2o2/therm.dat");
  therm.replace(therm.find("E+00"), 4, "XX+0");
  CHECK(!import_chemkin(read("mechanisms/h2o2/chem.inp"), therm, "", d3));
  CHECK(d3.contains(DiagCode::MechanismSyntax));
  // YAML syntax errors carry the line number.
  Diagnostics d4;
  CHECK(!import_cantera_yaml("species:\n- name: A\n  composition: {H: 2\n", "", d4));
  CHECK(d4.items().front().message.find("YAML line") != std::string::npos);
}
