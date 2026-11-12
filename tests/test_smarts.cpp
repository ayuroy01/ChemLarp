// SMARTS matching and Morgan fingerprints. Every expected count below is
// RDKit 2026.09's answer for the same SMILES and SMARTS (unique matches;
// GetSparseCountFingerprint totals and distinct identifiers).
#include <string>
#include <vector>

#include "omnichem/fingerprint.hpp"
#include "omnichem/rings.hpp"
#include "omnichem/smarts.hpp"
#include "omnichem/smiles.hpp"
#include "test_framework.hpp"

using namespace omnichem;

namespace {

std::size_t count(const char* smiles, const char* smarts, bool chiral = false, bool unique = true) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smiles, d);
  if (sp < 0) throw octest::Failure(std::string("bad SMILES ") + smiles);
  const SmartsQuery q = parse_smarts(smarts, d);
  if (!q.valid()) throw octest::Failure(std::string("bad SMARTS ") + smarts + ": " + d.items().front().message);
  SmartsMatchOptions o;
  o.use_chirality = chiral;
  o.unique = unique;
  o.max_matches = 0;
  return match_smarts(st, static_cast<std::uint32_t>(sp), q, o).size();
}

struct Case {
  const char* smiles;
  const char* smarts;
  std::size_t expected;
};

const Case kCases[] = {
    {"CC(=O)Oc1ccccc1C(=O)O", "[CX3](=O)[OX2H1]", 1}, {"CC(=O)Oc1ccccc1C(=O)O", "c", 6},
    {"CC(=O)Oc1ccccc1C(=O)O", "[c;H1]", 4},           {"CC(=O)Oc1ccccc1C(=O)O", "[#6]~[#8]", 5},
    {"CC(=O)Oc1ccccc1C(=O)O", "[R]", 6},              {"CC(=O)Oc1ccccc1C(=O)O", "[!R]", 7},
    {"CC(=O)Oc1ccccc1C(=O)O", "[D3]", 4},             {"CC(=O)Oc1ccccc1C(=O)O", "[X4]", 1},
    {"CC(=O)Oc1ccccc1C(=O)O", "[v4]", 9},             {"CC(=O)Oc1ccccc1C(=O)O", "[h1]", 5},
    {"CC(=O)Oc1ccccc1C(=O)O", "[$(C=O)]", 2},         {"CC(=O)Oc1ccccc1C(=O)O", "[$([OH]C=O)]", 1},
    {"CC(=O)Oc1ccccc1C(=O)O", "C=,#O", 2},            {"CC(=O)Oc1ccccc1C(=O)O", "[C,N;!R]", 3},
    {"CC(=O)Oc1ccccc1C(=O)O", "c:c", 6},              {"CC(=O)Oc1ccccc1C(=O)O", "[+0]", 13},
    {"CC(=O)Oc1ccccc1C(=O)O", "*", 13},
    // bicyclo[2.2.2]octane: three relevant six-rings, so bridgeheads are in 3 rings
    {"C1CC2CCC1CC2", "[R3]", 2}, {"C1CC2CCC1CC2", "[R2]", 6}, {"C1CC2CCC1CC2", "[r6]", 8},
    {"C1CC2CCC1CC2", "[x3]", 2}, {"C1CC2CCC1CC2", "[x2]", 6}, {"C1CC2CCC1CC2", "C1CCCCC1", 3},
    {"C1CC2CCC1CC2", "C@C", 9},
    {"C[N+](C)(C)CC(=O)[O-]", "[N+]", 1}, {"C[N+](C)(C)CC(=O)[O-]", "[O-]", 1},
    {"C[N+](C)(C)CC(=O)[O-]", "[#7X4+]", 1}, {"C[N+](C)(C)CC(=O)[O-]", "[OX1-]", 1},
    {"Nc1ccc(cc1)S(=O)(=O)N", "[NX3;H2]", 2}, {"Nc1ccc(cc1)S(=O)(=O)N", "[S;X4;v6](=O)(=O)", 1},
    {"Nc1ccc(cc1)S(=O)(=O)N", "[NH2]c", 1}, {"Nc1ccc(cc1)S(=O)(=O)N", "[nH]", 0},
    {"c1ccc2[nH]ccc2c1", "[nH]", 1}, {"c1ccc2[nH]ccc2c1", "[r5]", 5}, {"c1ccc2[nH]ccc2c1", "[R2]", 2},
    {"c1ccc2[nH]ccc2c1", "a1aaaa1", 1}, {"c1ccc2[nH]ccc2c1", "[c;R2]", 2},
    {"[2H]C([2H])([2H])O", "[2H]", 3}, {"[2H]C([2H])([2H])O", "[CH0]", 0}, {"[2H]C([2H])([2H])O", "[C;H3]", 1},
    {"[2H]C([2H])([2H])O", "[OH]", 1},
};

std::pair<std::size_t, std::size_t> envs(const char* smiles, int radius) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smiles, d);
  MorganOptions o;
  o.radius = radius;
  std::size_t total = 0;
  const auto e = morgan_environments(st, static_cast<std::uint32_t>(sp), o);
  for (const auto& x : e) total += x.second;
  return {total, e.size()};
}

}  // namespace

TEST(smarts_counts_match_rdkit) {
  for (const Case& c : kCases) {
    const std::size_t got = count(c.smiles, c.smarts);
    if (got != c.expected)
      throw octest::Failure(std::string(c.smiles) + " / " + c.smarts + ": " + std::to_string(got) + " matches, expected " +
                            std::to_string(c.expected));
  }
  CHECK_EQ(count("c1ccccc1", "c1ccccc1", false, false), std::size_t{12});  // every automorphism
  CHECK_EQ(count("c1ccccc1", "c1ccccc1"), std::size_t{1});
}

TEST(smarts_chirality) {
  // L-alanine written two ways; the query is L.
  CHECK_EQ(count("C[C@H](N)C(=O)O", "N[C@@H](C)C(=O)O", true), std::size_t{1});
  CHECK_EQ(count("N[C@@H](C)C(=O)O", "N[C@@H](C)C(=O)O", true), std::size_t{1});
  CHECK_EQ(count("C[C@@H](N)C(=O)O", "N[C@@H](C)C(=O)O", true), std::size_t{0});
  CHECK_EQ(count("C[C@@H](N)C(=O)O", "N[C@@H](C)C(=O)O", false), std::size_t{1});  // ignored by default
  CHECK_EQ(count("CC(N)C(=O)O", "N[C@@H](C)C(=O)O", true), std::size_t{0});         // no stereo in the target
  // A query that leaves out the hydrogen: [C@@](N)(C)C(=O)O names three neighbours.
  CHECK_EQ(count("C[C@H](N)C(=O)O", "N[C@@](C)C(=O)O", true), std::size_t{1});
}

TEST(smarts_errors) {
  for (const char* bad : {"", "[C", "C(", "C)", "C1CC", "[Q]", "C==", "[#]", "$(C", "C>C", "(C).(C)"}) {
    Diagnostics d;
    CHECK(!parse_smarts(bad, d).valid());
    CHECK(d.has_errors());
  }
  Diagnostics d;
  parse_smarts("C>>C", d);
  CHECK(d.contains(DiagCode::SmartsUnsupported));
}

TEST(relevant_rings) {
  ChemicalSystemState st;
  Diagnostics d;
  // bicyclo[2.2.2]octane: SSSR has 2 rings, relevant cycles 3; norbornane: 2 and 2.
  auto sp = static_cast<std::uint32_t>(parse_smiles(st, "C1CC2CCC1CC2", d));
  CHECK_EQ(find_sssr(st, sp).size(), std::size_t{2});
  CHECK_EQ(find_relevant_rings(st, sp).size(), std::size_t{3});
  sp = static_cast<std::uint32_t>(parse_smiles(st, "C1CC2CCC1C2", d));
  CHECK_EQ(find_relevant_rings(st, sp).size(), std::size_t{2});
  sp = static_cast<std::uint32_t>(parse_smiles(st, "c1ccc2ccccc2c1", d));  // naphthalene: the 10-ring is a sum
  CHECK_EQ(find_relevant_rings(st, sp).size(), std::size_t{2});
  sp = static_cast<std::uint32_t>(parse_smiles(st, "C12C3C4C1C5C2C3C45", d));  // cubane: six 4-rings
  CHECK_EQ(find_relevant_rings(st, sp).size(), std::size_t{6});
}

TEST(morgan_environments_match_rdkit) {
  struct E { const char* smiles; int radius; std::size_t total, distinct; };
  const E cases[] = {
      {"c1ccccc1", 1, 12, 2}, {"c1ccccc1", 2, 18, 3}, {"c1ccccc1", 3, 19, 4},
      {"CC(=O)Oc1ccccc1C(=O)O", 1, 26, 17}, {"CC(=O)Oc1ccccc1C(=O)O", 2, 35, 25}, {"CC(=O)Oc1ccccc1C(=O)O", 3, 42, 32},
      {"CCO", 2, 6, 6}, {"C1CC2CCC1CC2", 1, 16, 4}, {"C1CC2CCC1CC2", 2, 24, 6}, {"C1CC2CCC1CC2", 3, 25, 7},
      {"CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 2, 37, 25}, {"CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 3, 46, 34},
  };
  for (const E& e : cases) {
    const auto [total, distinct] = envs(e.smiles, e.radius);
    if (total != e.total || distinct != e.distinct)
      throw octest::Failure(std::string(e.smiles) + " r" + std::to_string(e.radius) + ": " + std::to_string(total) + "/" +
                            std::to_string(distinct) + ", expected " + std::to_string(e.total) + "/" +
                            std::to_string(e.distinct));
  }
}

TEST(morgan_similarity) {
  ChemicalSystemState st;
  Diagnostics d;
  const auto a = static_cast<std::uint32_t>(parse_smiles(st, "CC(=O)Oc1ccccc1C(=O)O", d));
  const auto b = static_cast<std::uint32_t>(parse_smiles(st, "OC(=O)c1ccccc1O", d));
  const auto a2 = static_cast<std::uint32_t>(parse_smiles(st, "OC(=O)c1ccccc1OC(C)=O", d));  // aspirin again
  const Fingerprint fa = morgan_fingerprint(st, a), fb = morgan_fingerprint(st, b), fa2 = morgan_fingerprint(st, a2);
  CHECK_EQ(fa.bits, 2048);
  CHECK_EQ(tanimoto(fa, fa2), 1.0);  // spelling-independent
  CHECK_NEAR(tanimoto(fa, fb), tanimoto(fb, fa), 1e-15);
  CHECK(tanimoto(fa, fb) > 0.3 && tanimoto(fa, fb) < 0.7);
  CHECK(dice(fa, fb) >= tanimoto(fa, fb));
  CHECK(fa.count() <= 25);  // at most one bit per distinct environment
  MorganOptions small;
  small.bits = 64;
  bool threw = false;
  try { tanimoto(fa, morgan_fingerprint(st, b, small)); } catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw);
}
