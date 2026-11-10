#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "omnichem/depict.hpp"
#include "omnichem/smiles.hpp"
#include "test_framework.hpp"

using namespace omnichem;

namespace {

struct Parsed {
  ChemicalSystemState st;
  Diagnostics d;
  std::uint32_t sp = 0;
};

Parsed parse(const std::string& smi) {
  Parsed p;
  const auto sp = parse_smiles(p.st, smi, p.d);
  if (sp < 0) throw octest::Failure("parse failed: " + smi);
  p.sp = static_cast<std::uint32_t>(sp);
  return p;
}

const std::vector<std::string> kDrugs = {
    "CC(=O)Oc1ccccc1C(=O)O",                                                   // aspirin
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",                                              // caffeine
    "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",                                         // (S)-ibuprofen
    "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O",            // testosterone
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",                                  // glucose
    "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)Cc3ccccc3)C(=O)O)C",              // penicillin G
    "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",  // atorvastatin
    "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",  // cholesterol
    "c1ccc2ccccc2c1",                                                          // naphthalene
    "C1CCC2(CC1)CCCC2",                                                        // spiro[4.5]decane
    "CCCCCCCCCCCCCCCC(=O)O",                                                   // palmitic acid
};

// Chirality read back from the drawing alone: wedge = toward the viewer.
int drawn_parity(const Depiction& dp, std::uint32_t centre_idx, const std::array<std::int64_t, 4>& refs) {
  std::array<std::array<double, 3>, 4> v{};
  int missing = -1;
  for (int i = 0; i < 4; ++i) {
    const auto r = refs[static_cast<std::size_t>(i)];
    if (r < 0) {
      missing = i;
      continue;
    }
    const auto& c = dp.atoms[centre_idx];
    const auto& a = dp.atoms[static_cast<std::size_t>(r)];
    double dx = a.x - c.x, dy = -(a.y - c.y);  // SVG y points down
    const double n = std::hypot(dx, dy);
    dx /= n;
    dy /= n;
    double z = 0.0;
    for (const auto& b : dp.bonds)
      if (b.a == centre_idx && b.b == static_cast<std::uint32_t>(r))
        z = b.stereo == BondStereo::Wedge ? 0.8 : b.stereo == BondStereo::Hash ? -0.8 : 0.0;
    v[static_cast<std::size_t>(i)] = {dx, dy, z};
  }
  if (missing >= 0) {
    std::array<double, 3> s{0, 0, 0};
    for (int i = 0; i < 4; ++i)
      if (i != missing)
        for (int k = 0; k < 3; ++k) s[static_cast<std::size_t>(k)] -= v[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
    v[static_cast<std::size_t>(missing)] = s;
  }
  auto sub = [](auto x, auto y) { return std::array<double, 3>{x[0] - y[0], x[1] - y[1], x[2] - y[2]}; };
  const auto b = sub(v[1], v[0]), c = sub(v[2], v[0]), d = sub(v[3], v[0]);
  const double vol = b[0] * (c[1] * d[2] - c[2] * d[1]) + b[1] * (c[2] * d[0] - c[0] * d[2]) + b[2] * (c[0] * d[1] - c[1] * d[0]);
  return vol < 0 ? 1 : 2;  // '@' : '@@'
}

}  // namespace

TEST(depict_bond_lengths_are_uniform) {
  for (const auto& smi : kDrugs) {
    Parsed p = parse(smi);
    const auto xy = layout_2d(p.st, p.sp);
    const Species& s = p.st.species[p.sp];
    for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
      const auto& u = xy[p.st.bonds.begin[b] - s.atom_begin];
      const auto& v = xy[p.st.bonds.end[b] - s.atom_begin];
      const double len = std::hypot(u[0] - v[0], u[1] - v[1]);
      if (std::fabs(len - 1.0) > 0.02) throw octest::Failure(smi + ": bond length " + std::to_string(len));
    }
  }
}

TEST(depict_has_no_atom_overlaps) {
  for (const auto& smi : kDrugs) {
    Parsed p = parse(smi);
    const auto xy = layout_2d(p.st, p.sp);
    const Species& s = p.st.species[p.sp];
    for (std::uint32_t i = 0; i < s.atom_count; ++i)
      for (std::uint32_t j = i + 1; j < s.atom_count; ++j) {
        const double d = std::hypot(xy[i][0] - xy[j][0], xy[i][1] - xy[j][1]);
        if (d < 0.6) throw octest::Failure(smi + ": atoms " + std::to_string(i + 1) + " and " + std::to_string(j + 1) + " overlap");
      }
  }
}

TEST(depict_rings_are_regular_polygons) {
  Parsed p = parse("c1ccccc1");
  const auto xy = layout_2d(p.st, p.sp);
  double cx = 0, cy = 0;
  for (const auto& a : xy) {
    cx += a[0] / 6.0;
    cy += a[1] / 6.0;
  }
  for (const auto& a : xy) CHECK_NEAR(std::hypot(a[0] - cx, a[1] - cy), 1.0, 1e-9);
}

TEST(depict_wedges_encode_every_stereocentre) {
  for (const auto& smi : kDrugs) {
    Parsed p = parse(smi);
    Diagnostics d;
    const Depiction dp = depict(p.st, p.sp, d);
    CHECK(dp.ok);
    CHECK_EQ(dp.stereo_centres_ok, dp.stereo_centres);
    const std::uint32_t begin = p.st.species[p.sp].atom_begin;
    for (const auto& t : p.st.tetrahedral) {
      if (p.st.atoms.species[t.atom] != p.sp) continue;
      const std::uint32_t c = t.atom - begin;
      std::array<std::int64_t, 4> refs{};
      for (int i = 0; i < 4; ++i) {
        const auto r = t.refs[static_cast<std::size_t>(i)];
        if (r == kImplicitHRef) {
          refs[static_cast<std::size_t>(i)] = -1;
          for (const auto& b : dp.bonds)  // an explicitly drawn hydrogen, if any
            if (b.a == c && dp.atoms[b.b].atom < 0) refs[static_cast<std::size_t>(i)] = b.b;
        } else {
          refs[static_cast<std::size_t>(i)] = r - begin;
        }
      }
      const int expect = t.parity == Chirality::CounterClockwise ? 1 : 2;
      if (drawn_parity(dp, c, refs) != expect) throw octest::Failure(smi + ": wrong wedge at atom " + std::to_string(c + 1));
    }
  }
}

TEST(depict_keeps_double_bond_geometry) {
  for (const char* smi : {"F/C=C/F", "F/C=C\\F", "C/C=C/C=C/C", "c1ccccc1/C=C\\c1ccccc1", "CC/C(C)=C(/CC)C"}) {
    Parsed p = parse(smi);
    Diagnostics d;
    const Depiction dp = depict(p.st, p.sp, d);
    CHECK(dp.stereo_bonds > 0);
    if (dp.stereo_bonds_ok != dp.stereo_bonds) throw octest::Failure(std::string(smi) + ": E/Z not kept");
  }
}

TEST(depict_ring_fusion_hydrogens_are_drawn) {
  Parsed p = parse("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O");
  Diagnostics d;
  const Depiction dp = depict(p.st, p.sp, d);
  std::size_t drawn_h = 0, wedges = 0;
  for (const auto& a : dp.atoms) drawn_h += a.atom < 0;
  for (const auto& b : dp.bonds) wedges += b.stereo != BondStereo::None;
  CHECK_EQ(drawn_h, 3u);  // 8β-H, 9α-H, 14α-H
  CHECK_EQ(wedges, 6u);
}

TEST(depict_labels_and_svg) {
  Diagnostics d;
  const Depiction eth = depict_smiles("CCO", d);
  CHECK(eth.ok);
  CHECK(eth.svg.rfind("<svg", 0) == 0);
  // the hydrogen sits on the side away from the bond: OH or HO
  const bool oh = eth.atoms[2].label == "OH";
  CHECK(oh || eth.atoms[2].label == "HO");
  CHECK(eth.svg.find(oh ? ">O</tspan><tspan>H</tspan>" : ">H</tspan><tspan>O</tspan>") != std::string::npos);
  CHECK((eth.atoms[2].x > eth.atoms[1].x) == oh);
  const Depiction nh4 = depict_smiles("[NH4+]", d);
  CHECK_EQ(nh4.atoms[0].label, std::string("NH4+"));
  CHECK(nh4.svg.find("&#8722;") == std::string::npos);
  const Depiction water = depict_smiles("O", d);
  CHECK_EQ(water.atoms[0].label, std::string("H2O"));
  const Depiction cl = depict_smiles("[Cl-]", d);
  CHECK(cl.svg.find("&#8722;") != std::string::npos);
  const Depiction amp = depict_smiles("C[N+](C)(C)C", d);  // title escaping keeps the markup valid
  CHECK(amp.svg.find("<title>") != std::string::npos);
}

TEST(depict_separates_fragments) {
  Diagnostics d;
  const Depiction dp = depict_smiles("[Na+].[Cl-]", d);
  CHECK_EQ(dp.atoms.size(), 2u);
  CHECK(std::fabs(dp.atoms[0].x - dp.atoms[1].x) > 30.0);
}

TEST(depict_molfile) {
  Diagnostics d;
  const Depiction dp = depict_smiles("N[C@@H](C)C(=O)O", d);
  CHECK(dp.molfile.find("V2000") != std::string::npos);
  CHECK(dp.molfile.find("2D") != std::string::npos);
  CHECK(dp.molfile.find("M  END") != std::string::npos);
  // counts line: 6 atoms, 5 bonds, chiral flag set
  CHECK(dp.molfile.find("  6  5  0  0  1") != std::string::npos);
  std::size_t stereo = 0;
  for (const auto& b : dp.bonds) stereo += b.stereo != BondStereo::None;
  CHECK_EQ(stereo, 1u);
}
