// Planar (cyclophane, metallocene) and helical (helicene) chirality from 3D coordinates.
//
// Reference: three crystal structures of known absolute configuration from the Crystallography Open
// Database (data/planar_chirality_cod.sdf, Flack parameters near 0): a bis[(Rp)-4-bromo[2.2]paracyclophane]
// ether and two (Sp) ferrocenes (1,2- and 1,3-disubstituted). Mirror images must give the opposite labels.
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "chemlarp/cip.hpp"
#include "chemlarp/conformer.hpp"
#include "chemlarp/molio.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

std::string labels(const ChemicalSystemState& st, std::int32_t sp) {
  std::string s;
  for (const auto& c : cip_planes(st, static_cast<std::uint32_t>(sp))) s += (s.empty() ? "" : " ") + std::string(to_string(c.label));
  return s;
}

std::string swap_labels(std::string s) {
  for (auto& [a, b] : {std::pair<std::string, std::string>{"Rp", "Sp"}, {"P", "M"}}) {
    std::string out;
    for (std::size_t i = 0; i < s.size();) {
      if (s.compare(i, a.size(), a) == 0) out += b, i += a.size();
      else if (s.compare(i, b.size(), b) == 0) out += a, i += b.size();
      else out += s[i++];
    }
    s = out;
  }
  return s;
}

// z -> -z on the atom lines of every V2000 record.
std::string mirror(const std::string& sdf) {
  std::istringstream in(sdf);
  std::string l, out;
  int line = 0, natoms = 0;
  while (std::getline(in, l)) {
    if (line == 3) natoms = std::stoi(l.substr(0, 3));
    if (line >= 4 && line < 4 + natoms) {
      char b[16];
      std::snprintf(b, sizeof b, "%10.4f", -std::stod(l.substr(20, 10)));
      l.replace(20, 10, b);
    }
    out += l + "\n";
    ++line;
    if (l == "$$$$") line = 0, natoms = 0;
  }
  return out;
}

std::string embedded_sdf(const std::string& smiles, std::uint64_t seed = 42) {
  ChemicalSystemState st;
  Diagnostics d;
  EmbedOptions o;
  o.seed = seed;
  const auto sp = embed_smiles(st, smiles, d, nullptr, o);
  return sp < 0 ? std::string() : write_sdf(st, static_cast<std::uint32_t>(sp));
}

}  // namespace

TEST(planar_crystal_structures_match_published_descriptors) {
  std::ifstream f(std::string(CHEMLARP_TEST_DATA_DIR) + "/planar_chirality_cod.sdf");
  std::stringstream ss;
  ss << f.rdbuf();
  const std::string text = ss.str();
  CHECK(!text.empty());
  for (int mirrored = 0; mirrored < 2; ++mirrored) {
    ChemicalSystemState st;
    Diagnostics d;
    const auto recs = read_molecules(st, mirrored ? mirror(text) : text, d);
    CHECK_EQ(recs.size(), std::size_t{3});
    if (recs.size() != 3) continue;
    // Both decks of both paracyclophane units carry a substituent; all four are Rp.
    const std::string want[3] = {"Rp Rp Rp Rp", "Sp", "Sp"};
    for (int k = 0; k < 3; ++k) {
      CHECK(recs[k].species >= 0);
      const std::string got = labels(st, recs[k].species);
      CHECK_EQ(got, mirrored ? swap_labels(want[k]) : want[k]);
    }
    // The ferrocenes: hapto bonds cut to the ionic form; only the substituted ring is a chirality plane.
    CHECK(d.contains(DiagCode::HaptoBondsDisconnected));
    CHECK(st.species[static_cast<std::size_t>(recs[1].species)].canonical_smiles.find("[Fe+2]") != std::string::npos);
    for (int k = 1; k < 3; ++k) {
      int n = 0;
      for (const auto& u : st.planar)
        if (st.atoms.species[u.atoms[0]] == static_cast<std::uint32_t>(recs[k].species)) n += u.kind == PlanarKind::Metallocene;
      CHECK_EQ(n, 1);
    }
  }
}

TEST(planar_cyclophanes_from_embedding) {
  // 4-Bromo[2.2]paracyclophane: one chirality plane (the brominated deck); the parent hydrocarbon has none.
  const std::string sdf = embedded_sdf("Brc1cc2ccc1CCc1ccc(cc1)CC2");
  CHECK(!sdf.empty());
  ChemicalSystemState a, b, c;
  Diagnostics d;
  const auto ra = read_molecules(a, sdf, d), rb = read_molecules(b, mirror(sdf), d);
  CHECK(!ra.empty() && !rb.empty());
  const std::string l = labels(a, ra[0].species);
  CHECK(l == "Rp" || l == "Sp");
  CHECK_EQ(labels(b, rb[0].species), swap_labels(l));
  CHECK(a.planar.size() == 1 && a.planar[0].kind == PlanarKind::Cyclophane && a.planar[0].pilots.size() == 1);
  const auto rc = read_molecules(c, embedded_sdf("c1cc2ccc1CCc1ccc(cc1)CC2"), d);
  CHECK(!rc.empty() && c.planar.empty());
  // Planar stereo has no SMILES notation: the canonical string is the same for both enantiomers.
  CHECK_EQ(a.species[0].canonical_smiles, b.species[0].canonical_smiles);
}

TEST(planar_helicenes) {
  // [6]- and [5]helicene are helical; [4]helicene (configurationally labile), chrysene and picene (zigzag)
  // and planar coronene are not taken.
  for (auto [s, rings] : {std::pair<const char*, std::uint32_t>{"c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c12", 6},
                          {"c1ccc2c(c1)ccc1ccc3ccc4ccccc4c3c12", 5}}) {
    const std::string sdf = embedded_sdf(s);
    ChemicalSystemState a, b;
    Diagnostics d;
    const auto ra = read_molecules(a, sdf, d), rb = read_molecules(b, mirror(sdf), d);
    CHECK(!ra.empty() && !rb.empty());
    CHECK_EQ(a.planar.size(), std::size_t{1});
    if (a.planar.size() != 1) continue;
    const PlanarStereo& u = a.planar[0];
    CHECK(u.kind == PlanarKind::Helicene);
    CHECK_EQ(u.rings, rings);
    CHECK_EQ(u.atoms.size(), std::size_t{u.rings + 1});
    // P: right-handed, so the torsion angle in the middle of the inner rim is positive.
    const std::size_t m = u.atoms.size() / 2 - 2;
    const double t = dihedral_degrees(a, u.atoms[m], u.atoms[m + 1], u.atoms[m + 2], u.atoms[m + 3]);
    const std::string l = labels(a, ra[0].species);
    CHECK_EQ(l, std::string(t > 0 ? "P" : "M"));
    CHECK_EQ(labels(b, rb[0].species), swap_labels(l));
  }
  for (const char* s : {"c1ccc2c(c1)ccc1ccc3ccccc3c12", "c1ccc2c(c1)ccc1c2ccc2c1ccc1ccccc12",
                        "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61"}) {
    ChemicalSystemState a;
    Diagnostics d;
    const auto ra = read_molecules(a, embedded_sdf(s), d);
    CHECK(!ra.empty());
    CHECK(a.planar.empty());
  }
}
