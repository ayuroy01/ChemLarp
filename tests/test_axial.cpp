// Axial chirality: allenes from SMILES and 3D coordinates, hindered biaryls from 3D coordinates.
//
// The Ra/Sa descriptors are checked against the helicity definition (IUPAC 2013, P-93.5.1.1): viewed along
// the axis, the turn from the higher-ranked near ligand to the higher-ranked far ligand by the smaller angle
// is clockwise for P and counterclockwise for M, and Ra = M, Sa = P. A positive torsion angle near-ligand –
// near end – far end – far-ligand is a clockwise turn, so it means P (Sa). This is computed from the 3D
// structure here, independently of the permutation used by cip_label().
#include <cmath>
#include <string>

#include "chemlarp/canonical.hpp"
#include "chemlarp/cip.hpp"
#include "chemlarp/conformer.hpp"
#include "chemlarp/molio.hpp"
#include "chemlarp/smiles.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

std::string axis_label(const ChemicalSystemState& st, std::int32_t sp) {
  std::string s;
  for (const auto& a : cip_axes(st, static_cast<std::uint32_t>(sp))) s += (s.empty() ? "" : " ") + std::string(to_string(a.label));
  return s;
}

std::string canon(const std::string& smiles) {
  ChemicalSystemState st;
  Diagnostics d;
  const auto sp = parse_smiles(st, smiles, d);
  return sp < 0 ? std::string("?") : st.species[static_cast<std::size_t>(sp)].canonical_smiles;
}

// z -> -z on the atom lines of a V2000 molfile.
std::string mirror(const std::string& sdf) {
  std::string out;
  std::size_t pos = 0;
  int line = 0, natoms = 0;
  while (pos < sdf.size()) {
    std::size_t e = sdf.find('\n', pos);
    if (e == std::string::npos) e = sdf.size();
    std::string l = sdf.substr(pos, e - pos);
    if (line == 3) natoms = std::stoi(l.substr(0, 3));
    if (line >= 4 && line < 4 + natoms) {
      char b[16];
      std::snprintf(b, sizeof b, "%10.4f", -std::stod(l.substr(20, 10)));
      l.replace(20, 10, b);
    }
    out += l + "\n";
    pos = e + 1;
    ++line;
  }
  return out;
}

std::int32_t read_first(ChemicalSystemState& st, const std::string& text) {
  Diagnostics d;
  const auto recs = read_molecules(st, text, d);
  return recs.empty() ? -1 : recs[0].species;
}

}  // namespace

TEST(axial_allene_smiles_labels_and_canonical) {
  ChemicalSystemState st;
  Diagnostics d;
  const auto a = parse_smiles(st, "CC=[C@]=CC", d);
  CHECK(a >= 0);
  CHECK_EQ(st.axial.size(), std::size_t{1});
  CHECK_EQ(axis_label(st, a), std::string("Ra"));
  CHECK_EQ(st.species[static_cast<std::size_t>(a)].canonical_smiles, std::string("CC=[C@]=CC"));
  // Same molecule written from the other end and with the implicit H first: same canonical string.
  CHECK_EQ(canon("C(C)=[C@@]=CC"), canon("CC=[C@]=CC"));
  CHECK_EQ(canon("FC(C)=[C@]=C(Cl)O"), canon("OC(Cl)=[C@]=C(C)F"));  // reversed: an even permutation
  // Enantiomers differ.
  CHECK(canon("CC=[C@@]=CC") != canon("CC=[C@]=CC"));
  // Longer cumulenes with an even number of double bonds are axial too; canonical strings round-trip.
  for (const char* s : {"CC=C=[C@]=C=CC", "OC(Cl)=[C@]=C(C)F", "C[C@H](O)C=[C@]=C[C@@H](C)O", "C1CCCC=[C@]=CCC1"}) {
    const std::string c1 = canon(s);
    CHECK(c1.find('@') != std::string::npos);
    CHECK_EQ(canon(c1), c1);
  }
}

TEST(axial_not_stereogenic_or_invalid) {
  // Two equal ligands on one end: no descriptor, nothing written, a diagnostic.
  ChemicalSystemState st;
  Diagnostics d;
  const auto a = parse_smiles(st, "CC(C)=[C@]=CC", d);
  CHECK(a >= 0);
  CHECK_EQ(axis_label(st, a), std::string(""));
  CHECK(st.species[static_cast<std::size_t>(a)].canonical_smiles.find('@') == std::string::npos);
  CHECK(!d.empty());
  // A butatriene (three double bonds) is planar: '@' on an inner atom is not axial and is ignored.
  ChemicalSystemState st2;
  Diagnostics d2;
  const auto b = parse_smiles(st2, "CC=[C@]=C=CC", d2);
  CHECK(b >= 0);
  CHECK(st2.axial.empty());
  CHECK(!d2.empty());
}

TEST(axial_allene_3d_round_trip_and_helicity) {
  for (const char* smiles : {"OC(Cl)=[C@]=C(C)F", "OC(Cl)=[C@@]=C(C)F"}) {
    ChemicalSystemState st;
    Diagnostics d;
    EmbedResult r;
    const auto sp = embed_smiles(st, smiles, d, &r);
    CHECK(sp >= 0);
    CHECK_EQ(r.stereo_centres_ok, r.stereo_centres);
    const std::string label = axis_label(st, sp);
    // The end groups are perpendicular (UFF cumulene term), and the descriptor follows the helicity:
    // Cl outranks O on the near end, F outranks C on the far end.
    const auto& s = st.species[static_cast<std::size_t>(sp)];
    std::uint32_t cl = 0, f = 0;
    for (std::uint32_t i = s.atom_begin; i < s.atom_begin + s.atom_count; ++i) {
      if (st.atoms.element[i] == 17) cl = i;
      if (st.atoms.element[i] == 9) f = i;
    }
    const AxialStereo& ax = st.axial.back();
    const std::uint32_t near_end = st.adjacency.neighbors(cl)[0], far_end = st.adjacency.neighbors(f)[0];
    const double tors = dihedral_degrees(st, cl, near_end, far_end, f);
    CHECK(std::fabs(std::fabs(tors) - 90.0) < 15.0);
    CHECK_EQ(label, std::string(tors > 0 ? "Sa" : "Ra"));
    CHECK(ax.centre != kNoAtom);
    // Through an SDF file and back, and its mirror image.
    const std::string sdf = write_sdf(st, static_cast<std::uint32_t>(sp));
    ChemicalSystemState s2, s3;
    const auto q = read_first(s2, sdf), m = read_first(s3, mirror(sdf));
    CHECK(q >= 0 && m >= 0);
    CHECK_EQ(axis_label(s2, q), label);
    CHECK_EQ(axis_label(s3, m), std::string(label == "Ra" ? "Sa" : "Ra"));
    CHECK_EQ(s2.species[static_cast<std::size_t>(q)].canonical_smiles, canon(smiles));
  }
}

TEST(axial_biaryl_from_3d) {
  // 1,1'-bi-2-naphthol (BINOL, four ortho positions substituted or fused) is atropisomeric; 2,2'-dimethyl-
  // biphenyl (two) and biphenyl (none) are not.
  ChemicalSystemState st;
  Diagnostics d;
  const auto sp = embed_smiles(st, "Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12", d);
  CHECK(sp >= 0);
  const std::string sdf = write_sdf(st, static_cast<std::uint32_t>(sp));
  ChemicalSystemState s2, s3;
  const auto q = read_first(s2, sdf), m = read_first(s3, mirror(sdf));
  CHECK(q >= 0 && m >= 0);
  CHECK_EQ(s2.axial.size(), std::size_t{1});
  const std::string l = axis_label(s2, q);
  CHECK(l == "Ra" || l == "Sa");
  CHECK_EQ(axis_label(s3, m), std::string(l == "Ra" ? "Sa" : "Ra"));
  CHECK(s2.axial[0].centre == kNoAtom);
  // Helicity: the higher-ranked ortho ligand of each ring is the C-OH carbon (O beats the fused carbon).
  auto ortho_oh = [&](std::uint32_t end) {
    for (std::uint32_t n : s2.adjacency.neighbors(end))
      for (std::uint32_t x : s2.adjacency.neighbors(n))
        if (s2.atoms.element[x] == 8) return n;
    return end;
  };
  const AxialStereo& ax = s2.axial[0];
  const double tors = dihedral_degrees(s2, ortho_oh(ax.begin), ax.begin, ax.end, ortho_oh(ax.end));
  CHECK_EQ(l, std::string(tors > 0 ? "Sa" : "Ra"));
  for (const char* s : {"Cc1ccccc1-c1ccccc1C", "c1ccccc1-c1ccccc1"}) {
    ChemicalSystemState a;
    Diagnostics da;
    const auto e = embed_smiles(a, s, da);
    ChemicalSystemState b;
    CHECK(read_first(b, write_sdf(a, static_cast<std::uint32_t>(e))) >= 0);
    CHECK(b.axial.empty());
  }
}
