#include <cstdint>

#include "omnichem/smiles.hpp"
#include "test_framework.hpp"

using namespace omnichem;

namespace {

struct Parsed {
  ChemicalSystemState st;
  Diagnostics d;
  std::int32_t sp = -1;
  const Species& s() const { return st.species.at(static_cast<std::size_t>(sp)); }
};

Parsed parse(const char* smi) {
  Parsed p;
  p.sp = parse_smiles(p.st, smi, p.d);
  return p;
}

int count_order(const Parsed& p, int order) {
  int n = 0;
  for (std::size_t b = 0; b < p.st.bonds.size(); ++b) n += p.st.bonds.order[b] == order;
  return n;
}

}  // namespace

TEST(smiles_ethanol_formula_and_implicit_h) {
  auto p = parse("CCO");
  CHECK(p.sp == 0);
  CHECK_EQ(p.s().formula, std::string("C2H6O"));
  CHECK_EQ(p.st.atoms.implicit_h[0], 3);
  CHECK_EQ(p.st.atoms.implicit_h[1], 2);
  CHECK_EQ(p.st.atoms.implicit_h[2], 1);
  CHECK_NEAR(p.s().molar_mass, 46.069, 0.01);
  CHECK(!(p.st.atoms.flags[0] & kAtomInRing));
}

TEST(smiles_benzene_kekulized) {
  auto p = parse("c1ccccc1");
  CHECK(p.sp == 0);
  CHECK_EQ(p.s().formula, std::string("C6H6"));
  CHECK_EQ(count_order(p, 2), 3);
  CHECK_EQ(count_order(p, 1), 3);
  for (std::size_t i = 0; i < 6; ++i) CHECK(p.st.atoms.flags[i] & kAtomInRing);
}

TEST(smiles_heteroaromatics) {
  CHECK_EQ(parse("c1ccncc1").s().formula, std::string("C5H5N"));     // pyridine
  CHECK_EQ(parse("c1cc[nH]c1").s().formula, std::string("C4H5N"));   // pyrrole
  CHECK_EQ(parse("c1ccsc1").s().formula, std::string("C4H4S"));      // thiophene
  CHECK_EQ(parse("c1ccoc1").s().formula, std::string("C4H4O"));      // furan
  CHECK_EQ(parse("c1ccc2ccccc2c1").s().formula, std::string("C10H8"));  // naphthalene
  CHECK_EQ(parse("Cn1cccc1").s().formula, std::string("C5H7N"));     // N-methylpyrrole
  CHECK_EQ(parse("c1cc[n+](C)cc1").s().formula, std::string("C6H8N+"));  // N-methylpyridinium
}

TEST(smiles_nitrobenzene_charge_separated) {
  auto p = parse("c1ccccc1[N+](=O)[O-]");
  CHECK(p.sp == 0);
  CHECK_EQ(p.s().formula, std::string("C6H5NO2"));
  CHECK_EQ(p.s().net_charge, 0);
}

TEST(smiles_bracket_atoms) {
  auto p = parse("[NH4+]");
  CHECK_EQ(p.s().formula, std::string("H4N+"));
  CHECK_EQ(p.s().net_charge, 1);
  CHECK_EQ(parse("[H][H]").s().formula, std::string("H2"));
  CHECK_EQ(parse("[C-]#[O+]").s().formula, std::string("CO"));
  CHECK_EQ(parse("[13CH4]").st.atoms.isotope[0], 13);
  CHECK_EQ(parse("[Na+].[Cl-]").s().formula, std::string("ClNa"));
  CHECK_EQ(parse("[Fe+3]").s().net_charge, 3);
}

TEST(smiles_chirality_and_directional_bonds) {
  auto p = parse("C[C@H](O)CC");
  CHECK(p.sp == 0);
  CHECK(p.st.atoms.chirality[1] == static_cast<std::int8_t>(Chirality::CounterClockwise));
  auto q = parse("C[C@@H](O)CC");
  CHECK(q.st.atoms.chirality[1] == static_cast<std::int8_t>(Chirality::Clockwise));
  auto e = parse("F/C=C/F");
  CHECK(e.sp == 0);
  CHECK(e.st.bonds.flags[0] & kBondUp);
}

TEST(smiles_ring_closures_and_percent) {
  CHECK_EQ(parse("C1CCCCC1").s().formula, std::string("C6H12"));
  CHECK_EQ(parse("C%10CCCCC%10").s().formula, std::string("C6H12"));
  CHECK_EQ(parse("C=1CCCCC1").s().formula, std::string("C6H10"));
}

// ---- error vectors ---------------------------------------------------------

TEST(error_valence_exceeded_carbon) {
  auto p = parse("C(C)(C)(C)(C)C");
  CHECK(p.sp < 0);
  CHECK(p.d.contains(DiagCode::ValenceExceeded));
  CHECK_EQ(p.d.items()[0].position, 0u);
  CHECK_EQ(p.d.items()[0].atom, 0);
}

TEST(error_valence_reports_every_bad_atom) {
  auto p = parse("C(C)(C)(C)(C)C.O(C)(C)C");
  CHECK(p.sp < 0);
  CHECK_EQ(p.d.count(DiagCode::ValenceExceeded), 2u);
}

TEST(error_valence_charge_aware) {
  CHECK(parse("C#O").d.contains(DiagCode::ValenceExceeded));  // O with valence 3 needs O+
  CHECK(parse("[C-]#[O+]").sp == 0);
  CHECK(parse("[F-]C").d.contains(DiagCode::ValenceExceeded));  // F- is isoelectronic with Ne
  CHECK(parse("C[N+](C)(C)C").sp == 0);                           // N+ behaves like C
  CHECK(parse("C[N+](C)(C)(C)C").d.contains(DiagCode::ValenceExceeded));
}

TEST(error_kekulization_failure) {
  auto p = parse("c1ccnc1");  // pyrrole nitrogen missing its H
  CHECK(p.sp < 0);
  CHECK(p.d.contains(DiagCode::KekulizationFailed));
  CHECK(parse("c1cccc1").d.contains(DiagCode::KekulizationFailed));
  auto both = parse("C(C)(C)(C)(C)C.c1ccnc1");  // independent errors are both reported
  CHECK(both.d.contains(DiagCode::KekulizationFailed));
  CHECK(both.d.contains(DiagCode::ValenceExceeded));
}

TEST(error_aromatic_outside_ring) {
  CHECK(parse("c").d.contains(DiagCode::InvalidAromaticAtom));
  CHECK(parse("Cc").d.contains(DiagCode::InvalidAromaticAtom));
}

TEST(error_syntax_positions) {
  auto a = parse("CC)C");
  CHECK(a.d.contains(DiagCode::UnmatchedCloseParen));
  CHECK_EQ(a.d.items()[0].position, 2u);
  CHECK(parse("C(C").d.contains(DiagCode::UnmatchedOpenParen));
  CHECK(parse("C1CC").d.contains(DiagCode::UnclosedRing));
  CHECK(parse("C=").d.contains(DiagCode::BondWithoutAtom));
  CHECK(parse("=C").d.contains(DiagCode::BondWithoutAtom));
  CHECK(parse("C11").d.contains(DiagCode::RingSelfBond));
  CHECK(parse("C1C1").d.contains(DiagCode::DuplicateBond));
  CHECK(parse("C=1CC#1").d.contains(DiagCode::RingBondConflict));
  CHECK(parse("C()C").d.contains(DiagCode::EmptyBranch));
  CHECK(parse("").d.contains(DiagCode::EmptyInput));
  CHECK(parse("C C").d.contains(DiagCode::UnexpectedCharacter));
  CHECK(parse("[CH4").d.contains(DiagCode::UnclosedBracketAtom));
  CHECK(parse("[C@TH1](F)(Cl)(Br)I").d.contains(DiagCode::UnsupportedFeature));
}

TEST(error_unknown_elements) {
  CHECK(parse("[Xx]").d.contains(DiagCode::UnknownElement));
  auto na = parse("NaCl");
  CHECK(na.d.contains(DiagCode::UnknownElement));
  CHECK(na.d.items()[0].message.find("[Na]") != std::string::npos);
  CHECK(parse("CHC").d.contains(DiagCode::UnknownElement));
}

TEST(error_collects_multiple_syntax_errors) {
  auto p = parse("C(C))C1[Qq]");
  CHECK(p.d.contains(DiagCode::UnmatchedCloseParen));
  CHECK(p.d.contains(DiagCode::UnknownElement));
  CHECK(p.d.contains(DiagCode::UnclosedRing));
  CHECK(p.d.error_count() >= 3);
}

TEST(parse_is_transactional) {
  ChemicalSystemState st;
  Diagnostics d;
  CHECK(parse_smiles(st, "CCO", d) == 0);
  const auto atoms = st.atoms.size(), bonds = st.bonds.size();
  CHECK(parse_smiles(st, "C(C)(C)(C)(C)C", d) == -1);
  CHECK_EQ(st.atoms.size(), atoms);
  CHECK_EQ(st.bonds.size(), bonds);
  CHECK_EQ(st.species.size(), 1u);
  CHECK(parse_smiles(st, "O", d) == 1);
  CHECK_EQ(st.species[1].atom_begin, 3u);
  CHECK_EQ(st.adjacency.offsets.size(), st.atoms.size() + 1);
}

TEST(soa_columns_are_simd_aligned) {
  auto p = parse("CCCCCCCCCC");
  CHECK(reinterpret_cast<std::uintptr_t>(p.st.atoms.x.data()) % 64 == 0);
  CHECK(reinterpret_cast<std::uintptr_t>(p.st.atoms.element.data()) % 64 == 0);
  CHECK(reinterpret_cast<std::uintptr_t>(p.st.bonds.order.data()) % 64 == 0);
}
