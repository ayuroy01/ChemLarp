#include <set>
#include <string>
#include <vector>

#include "omnichem/canonical.hpp"
#include "omnichem/cip.hpp"
#include "omnichem/smiles.hpp"
#include "omnichem/thermo.hpp"
#include "test_framework.hpp"

using namespace omnichem;

namespace {

std::string canon(const std::string& smi) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smi, d);
  if (sp < 0) throw octest::Failure("failed to parse " + smi + ": " + d.items().front().message);
  return st.species[static_cast<std::size_t>(sp)].canonical_smiles;
}

std::vector<CipLabel> centre_labels(const std::string& smi) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smi, d);
  CHECK(sp >= 0);
  std::vector<CipLabel> out;
  for (const auto& c : cip_centres(st, static_cast<std::uint32_t>(sp))) out.push_back(c.label);
  return out;
}

CipLabel bond_label(const std::string& smi) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smi, d);
  CHECK(sp >= 0);
  auto v = cip_double_bonds(st, static_cast<std::uint32_t>(sp));
  CHECK_EQ(v.size(), 1u);
  return v[0].label;
}

void expect_same(std::initializer_list<const char*> spellings) {
  const std::string ref = canon(*spellings.begin());
  for (const char* s : spellings)
    if (canon(s) != ref) throw octest::Failure(std::string(s) + " -> " + canon(s) + " but expected " + ref);
}

// A corpus covering aromatics, heteroatoms, charges, isotopes, rings, fused
// systems and every stereo feature.
const std::vector<std::string> kCorpus = {
    "CCO",
    "c1ccccc1",
    "Cc1ccc(O)cc1",
    "c1ccncc1",
    "c1cc[nH]c1",
    "Cn1cccc1",
    "c1ccoc1",
    "c1ccsc1",
    "O=c1cccc[nH]1",
    "c1ccc2ccccc2c1",
    "c1ccc2cccc2cc1",
    "c1ccc2[nH]ccc2c1",
    "c1ccc(cc1)-c1ccccc1",
    "C1CC2CCC1CC2",
    "CC(=O)Oc1ccccc1C(=O)O",
    "[O-][N+](=O)c1ccccc1",
    "C[N+](C)(C)C",
    "[NH4+].[Cl-]",
    "[2H]C([2H])([2H])O",
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "C[C@@H](C(=O)O)N",
    "CC[C@@H](C)O",
    "F[C@H](Cl)Br",
    "C[C@H]1CCCC[C@@H]1C",
    "C[C@@H](O)[C@H](C)O",
    "C[S@](=O)CC",
    "F/C=C/F",
    "F/C=C\\F",
    "C/C=C/C=C/C",
    "C/C=C\\C=C/C",
    "CC/C(C)=C(\\C)CO",
    "OC(=O)/C=C/c1ccccc1",
    "C/C=C/[C@H](O)C",
};

}  // namespace

TEST(canonical_equivalent_spellings) {
  expect_same({"CCO", "OCC", "C(O)C", "[CH3][CH2][OH]"});
  expect_same({"c1ccccc1", "C1=CC=CC=C1", "C=1C=CC=CC=1"});
  expect_same({"Cc1ccccc1", "CC1=CC=CC=C1", "c1ccc(C)cc1", "C1=CC(C)=CC=C1"});
  expect_same({"c1ccncc1", "C1=CC=NC=C1", "n1ccccc1"});
  expect_same({"c1cc[nH]c1", "C1=CNC=C1", "[nH]1cccc1"});
  expect_same({"c1ccc2ccccc2c1", "C1=CC=C2C=CC=CC2=C1", "C1=CC2=CC=CC=C2C=C1"});
  expect_same({"C", "[CH4]", "[H]C([H])([H])[H]"});
  expect_same({"[O-][N+](=O)c1ccccc1", "c1ccccc1[N+]([O-])=O"});
  expect_same({"[Na+].[Cl-]", "[Cl-].[Na+]"});
}

TEST(canonical_distinguishes_isomers) {
  std::set<std::string> forms = {canon("CCO"), canon("COC"), canon("CCC"), canon("CC=C"), canon("C1CC1")};
  CHECK_EQ(forms.size(), 5u);
  CHECK(canon("[13CH4]") != canon("C"));
  CHECK(canon("F/C=C/F") != canon("F/C=C\\F"));
  CHECK(canon("F[C@H](Cl)Br") != canon("F[C@@H](Cl)Br"));
  CHECK(canon("C[C@@H](O)[C@H](C)O") != canon("C[C@@H](O)[C@@H](C)O"));  // meso vs chiral 2,3-butanediol
}

TEST(canonical_is_idempotent) {
  for (const auto& s : kCorpus) {
    const std::string c = canon(s);
    if (canon(c) != c) throw octest::Failure(s + ": canonical " + c + " re-canonicalises to " + canon(c));
  }
}

// The strongest check: many random atom orders, in aromatic and Kekulé form,
// must parse back to the same canonical string (this also proves the writer's
// stereo parity and '/' '\' placement are correct for every output order).
TEST(canonical_random_round_trip) {
  for (const auto& s : kCorpus) {
    ChemicalSystemState st;
    Diagnostics d;
    const std::int32_t sp = parse_smiles(st, s, d);
    CHECK(sp >= 0);
    const std::string ref = st.species[static_cast<std::size_t>(sp)].canonical_smiles;
    for (std::uint64_t seed = 1; seed <= 40; ++seed) {
      SmilesWriteOptions o;
      o.canonical = false;
      o.random_seed = seed;
      o.kekule = (seed % 2) == 0;
      const std::string rnd = write_smiles(st, static_cast<std::uint32_t>(sp), o);
      std::string back;
      try {
        back = canon(rnd);
      } catch (const octest::Failure& f) {
        throw octest::Failure(s + ": random form " + rnd + " does not parse: " + f.what());
      }
      if (back != ref)
        throw octest::Failure(s + ": random form " + rnd + " canonicalises to " + back + ", expected " + ref);
    }
  }
}

TEST(aromaticity_perception) {
  ChemicalSystemState st;
  Diagnostics d;
  const auto benz = parse_smiles(st, "C1=CC=CC=C1", d);
  for (std::uint32_t a = 0; a < 6; ++a) CHECK(st.atoms.flags[a] & kAtomAromatic);
  const auto chd = parse_smiles(st, "C1=CCC=CC1", d);  // 1,4-cyclohexadiene: not aromatic
  const auto& s = st.species[static_cast<std::size_t>(chd)];
  for (std::uint32_t a = s.atom_begin; a < s.atom_begin + s.atom_count; ++a) CHECK(!(st.atoms.flags[a] & kAtomAromatic));
  const auto cot = parse_smiles(st, "C1=CC=CC=CC=C1", d);  // cyclooctatetraene, 8 pi: not aromatic
  const auto& c8 = st.species[static_cast<std::size_t>(cot)];
  CHECK(!(st.atoms.flags[c8.atom_begin] & kAtomAromatic));
  CHECK(benz == 0);
  CHECK(canon("C1=CC=CC=CC=C1").find('c') == std::string::npos);  // written in Kekulé form
}

TEST(cip_reference_assignments) {
  // L-alanine, (S); written as in PubChem and as an alternative spelling.
  CHECK(centre_labels("C[C@@H](C(=O)O)N") == std::vector<CipLabel>{CipLabel::S});
  CHECK(centre_labels("N[C@@H](C)C(=O)O") == std::vector<CipLabel>{CipLabel::S});
  CHECK(centre_labels("N[C@H](C)C(=O)O") == std::vector<CipLabel>{CipLabel::R});  // D-alanine
  // (R)-butan-2-ol: ethyl vs methyl is decided in the second sphere.
  CHECK(centre_labels("CC[C@@H](C)O") == std::vector<CipLabel>{CipLabel::R});
  CHECK(centre_labels("F[C@H](Cl)Br") == std::vector<CipLabel>{CipLabel::R});
  CHECK(centre_labels("[C@@H](F)(Cl)Br") == std::vector<CipLabel>{CipLabel::R});
  // Isotopes break ties by Rule 2: [2H] outranks H. refs [2H, H, C, O] -> [H, O, C, 2H] is even, so @ = R.
  CHECK(centre_labels("[2H][C@H](C)O") == std::vector<CipLabel>{CipLabel::R});
  CHECK(centre_labels("[2H][C@@H](C)O") == std::vector<CipLabel>{CipLabel::S});
  CHECK(centre_labels("[2H][C@@H](C)O") != centre_labels("[H][C@@H](C)O"));  // plain H: no centre at all
}

TEST(cip_double_bonds) {
  CHECK(bond_label("F/C=C/F") == CipLabel::E);
  CHECK(bond_label("F/C=C\\F") == CipLabel::Z);
  CHECK(bond_label("C/C=C/C") == CipLabel::E);
  CHECK(bond_label("C(\\F)=C/F") == CipLabel::E);
  CHECK(bond_label("F/C(Cl)=C/F") == CipLabel::Z);   // Cl outranks F on the first carbon
  CHECK(bond_label("OC(=O)/C=C/c1ccccc1") == CipLabel::E);  // trans-cinnamic acid
}

TEST(cip_ties_report_and_drop_stereo) {
  ChemicalSystemState st;
  Diagnostics d;
  const auto sp = parse_smiles(st, "C[C@H](C)O", d);  // two methyls: not a stereocentre
  CHECK(sp >= 0);
  CHECK(d.contains(DiagCode::CipUndetermined));
  CHECK_EQ(st.species[static_cast<std::size_t>(sp)].canonical_smiles, std::string("CC(C)O"));
}

TEST(stereo_error_vectors) {
  ChemicalSystemState st;
  Diagnostics d;
  CHECK(parse_smiles(st, "F/C(\\F)=C/F", d) < 0);  // both fluorines marked onto the same side
  CHECK(d.contains(DiagCode::ConflictingBondStereo));
  Diagnostics d2;
  CHECK(parse_smiles(st, "C1/C=C/CC1", d2) >= 0);  // cis/trans in a 5-ring is implied
  CHECK(d2.contains(DiagCode::StereoIgnored));
  Diagnostics d3;
  CHECK(parse_smiles(st, "C[C@](C)=C", d3) >= 0);  // 3 neighbours, no lone pair element
  CHECK(st.tetrahedral.empty() || d3.contains(DiagCode::StereoIgnored) || d3.contains(DiagCode::CipUndetermined));
}

TEST(thermo_lookup_uses_canonical_structure) {
  // Different spellings of water and CO2 hit the same database entry.
  for (const char* smi : {"O", "[OH2]", "[H]O[H]", "O=C=O", "C(=O)=O"}) {
    ChemicalSystemState st;
    Diagnostics d;
    const auto sp = parse_smiles(st, smi, d);
    auto t = resolve_thermo(st, static_cast<std::uint32_t>(sp), 298.15, d);
    CHECK(t.has_value());
    CHECK(t->source == ThermoSource::Nasa7Database);
    CHECK(!d.contains(DiagCode::ThermoFormulaMatch));
  }
  // A different valence-bond spelling of CO falls back to formula, with a warning.
  ChemicalSystemState st;
  Diagnostics d;
  const auto sp = parse_smiles(st, "[C]=O", d);
  auto t = resolve_thermo(st, static_cast<std::uint32_t>(sp), 298.15, d);
  CHECK(t.has_value());
  CHECK(d.contains(DiagCode::ThermoFormulaMatch));
}
