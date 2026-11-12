// CIP rules beyond constitution (1b, 3, 4a-c, 5). Every expected descriptor
// here was cross-checked against RDKit's rdCIPLabeler (the "centres" algorithm
// of Hanson et al., J. Chem. Inf. Model. 58:1755, 2018).
#include <string>
#include <vector>

#include "chemlarp/cip.hpp"
#include "chemlarp/smiles.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

// "a<atom><label>" for centres and "b<bond begin>-<end><label>" for double bonds, in index order.
std::string labels(const std::string& smi) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smi, d);
  if (sp < 0) throw cltest::Failure("parse failed: " + smi);
  std::string out;
  for (const auto& c : cip_centres(st, static_cast<std::uint32_t>(sp)))
    if (c.label != CipLabel::None) out += "a" + std::to_string(c.index) + to_string(c.label) + " ";
  for (const auto& c : cip_double_bonds(st, static_cast<std::uint32_t>(sp)))
    if (c.label != CipLabel::None)
      out += "b" + std::to_string(st.bonds.begin[c.index]) + "-" + std::to_string(st.bonds.end[c.index]) + to_string(c.label) + " ";
  if (!out.empty()) out.pop_back();
  return out;
}

void expect(const std::string& smi, const std::string& want) {
  const std::string got = labels(smi);
  if (got != want) throw cltest::Failure(smi + ": got [" + got + "], expected [" + want + "]");
}

}  // namespace

TEST(cip_pseudoasymmetric_centres) {
  // pentane-2,3,4-triol: C3 is r/s only when its two branches are mirror images
  expect("C[C@H](O)[C@H](O)[C@@H](C)O", "a1S a3s a5R");
  expect("C[C@H](O)[C@@H](O)[C@H](C)O", "a1S a5S");  // (S,S): C3 is not a stereocentre
  expect("C[C@H](O)[C@H](O)[C@H](C)O", "a1S a5S");
  // ring cis/trans isomers with no chiral centre
  expect("C1[C@H](CC[C@H](C)C1)C", "a1r a4r");
  expect("C1[C@@H](C[C@@H](C[C@@H]1Cl)Cl)Cl", "a1s a3s a5s");
  // pseudoasymmetric phosphorus
  expect("CC[C@@H](C)O[P@@](C)(=O)O[C@@H](C)CC", "a2R a5s a9S");
}

TEST(cip_rule_3_z_precedes_e) {
  expect("C/C=C/[C@H](O)/C=C\\C", "a3S b1-2E b5-6Z");
  expect("C/C=C/[C@@H](O)/C=C\\C", "a3R b1-2E b5-6Z");
}

TEST(cip_rule_4b_like_precedes_unlike) {
  expect("C[C@H](Cl)[C@H](F)[C@@H](F)[C@H](C)Cl", "a1S a3R a5R a7S");
  expect("C[C@@H](O)C(=C[C@H](C)O)[C@H](C)O", "a1R a5S a8S");
  expect("OC(=O)[C@H](O)[C@@H](O)C(=O)O", "a3R a5R");  // L-tartaric acid
  expect("OC(=O)[C@H](O)[C@H](O)C(=O)O", "a3R a5S");   // meso
}

TEST(cip_aromatic_duplicates_and_ring_closures) {
  // Kekulé-independent ranking of pyridyl against phenyl and against each other
  expect("O[C@@H](c1ccccc1)c1ccccn1", "a1S");
  expect("C[C@H]([C@@H](O)[C@@H](c1ccccn1)C)c1ncccc1", "a1S a2r a4R");
  // spiro centre whose rings differ only by the far stereocentres
  expect("O[C@H]1CC[C@@]2(CC[C@H](O)CC2)CC1", "a1S a4R a7S");
  // sulfoxide lone pair; exocyclic double bond with stereogenic ring atoms
  expect("C[S@](=O)CC", "a1S");
  expect("F/C=C1/C[C@@H](C)C[C@H](C)C1", "a4S a7S");
}

TEST(smiles_direction_mark_at_ring_closure) {
  // the mark before a closing digit reads from the closing atom
  expect("C1=C/c2cccc(n2)/C=C/c2ccc/1cc2", "b0-1E b8-9E");
  expect("C1=C\\c2cccc(n2)/C=C/c2ccc/1cc2", "b0-1Z b8-9E");
  expect("F/C=C/1CCCCC1C", "b1-2E");  // mark at the opening digit
  expect("F/C=C(/C1C)CCCC1", "b1-2E");
  expect("F/C=C1/CCCCC1C", "b1-2Z");
}
