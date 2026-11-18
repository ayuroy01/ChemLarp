// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
// Docking: the scoring function against hand-computed values, the pose gradient against finite differences,
// symmetry-corrected RMSD, and a redock into a synthetic pocket.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "chemlarp/docking.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

// A file under data/, valid from any working directory.
std::string read_data_file(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

// --- the scoring function ---------------------------------------------------------------------

// Two hydrophobic carbons 4.0 A apart: surface distance d = 4.0 - 1.9 - 1.9 = 0.2, so gauss1 =
// exp(-(0.2/0.5)^2), gauss2 = exp(-((0.2-3)/2)^2), no repulsion (d > 0), the hydrophobic term on its
// plateau (d < 0.5) and no hydrogen bond.
TEST(docking_scoring_terms_of_a_carbon_pair) {
  const ScoringTerms t = pair_terms(XsType::C_H, XsType::C_H, 4.0);
  CHECK_NEAR(t.gauss1, std::exp(-0.16), 1e-12);
  CHECK_NEAR(t.gauss1, 0.8521437890, 1e-9);
  CHECK_NEAR(t.gauss2, std::exp(-1.96), 1e-12);
  CHECK_NEAR(t.gauss2, 0.1408584209, 1e-9);
  CHECK_EQ(t.repulsion, 0.0);
  CHECK_EQ(t.hydrophobic, 1.0);
  CHECK_EQ(t.hbond, 0.0);
  // -0.035579*0.85214379 - 0.005156*0.14085842 - 0.035069*1
  CHECK_NEAR(pair_energy(XsType::C_H, XsType::C_H, 4.0), -0.0661136899, 1e-9);
}

// A donor nitrogen and an acceptor oxygen 2.9 A apart: d = 2.9 - 1.8 - 1.7 = -0.6, so the pair is inside
// each other's radii (repulsion 0.36) and the hydrogen-bond term is on its ramp, 0.6/0.7 of the way up.
TEST(docking_scoring_terms_of_a_hydrogen_bond) {
  const ScoringTerms t = pair_terms(XsType::N_D, XsType::O_A, 2.9);
  CHECK_NEAR(t.gauss1, std::exp(-1.44), 1e-12);
  CHECK_NEAR(t.gauss2, std::exp(-3.24), 1e-12);
  CHECK_NEAR(t.repulsion, 0.36, 1e-12);
  CHECK_EQ(t.hydrophobic, 0.0);
  CHECK_NEAR(t.hbond, 6.0 / 7.0, 1e-12);
  CHECK_NEAR(pair_energy(XsType::N_D, XsType::O_A, 2.9), -0.2096625246, 1e-9);
  // The same pair with the nitrogen unable to donate: only the hydrogen-bond term goes.
  const double no_hbond = pair_energy(XsType::N_P, XsType::O_A, 2.9);
  CHECK_NEAR(no_hbond - pair_energy(XsType::N_D, XsType::O_A, 2.9), 0.587439 * 6.0 / 7.0, 1e-9);
  // Closer than -0.7 the term saturates at 1.
  CHECK_EQ(pair_terms(XsType::N_D, XsType::O_A, 2.5).hbond, 1.0);
  // Beyond touching it is gone.
  CHECK_EQ(pair_terms(XsType::N_D, XsType::O_A, 3.6).hbond, 0.0);
}

TEST(docking_scoring_cutoff_and_boundaries) {
  CHECK_EQ(pair_energy(XsType::C_H, XsType::C_H, 8.0), 0.0);
  CHECK_EQ(pair_energy(XsType::C_H, XsType::C_H, 50.0), 0.0);
  CHECK(pair_energy(XsType::C_H, XsType::C_H, 7.99) != 0.0);
  // An unknown element scores nothing, and so does a distance that is not a number.
  CHECK_EQ(pair_energy(XsType::None, XsType::C_H, 4.0), 0.0);
  CHECK_EQ(pair_energy(XsType::C_H, XsType::C_H, std::nan("")), 0.0);
  // Overlapping atoms are pure repulsion: d = -3.8, so 0.840245 * 14.44 dominates.
  CHECK(pair_energy(XsType::C_H, XsType::C_H, 0.0) > 11.0);
}

TEST(docking_pair_energy_derivative_matches_finite_difference) {
  const XsType pairs[][2] = {{XsType::C_H, XsType::C_H},
                             {XsType::N_D, XsType::O_A},
                             {XsType::C_H, XsType::O_DA},
                             {XsType::Met_D, XsType::O_A},
                             {XsType::Cl_H, XsType::C_H}};
  for (const auto& p : pairs)
    for (double r = 0.5; r < 7.9; r += 0.07) {
      double analytic = 0;
      pair_energy(p[0], p[1], r, {}, &analytic);
      const double h = 1e-6;
      const double numeric = (pair_energy(p[0], p[1], r + h) - pair_energy(p[0], p[1], r - h)) / (2 * h);
      // The hydrophobic and hydrogen-bond ramps have corners; skip a window around each.
      const double d = r - xs_radius(p[0]) - xs_radius(p[1]);
      const bool corner = std::fabs(d - 0.5) < 1e-4 || std::fabs(d - 1.5) < 1e-4 || std::fabs(d + 0.7) < 1e-4 ||
                          std::fabs(d) < 1e-4;
      if (corner) continue;
      CHECK_NEAR(analytic, numeric, 1e-5 * std::max(1.0, std::fabs(numeric)));
    }
}

TEST(docking_atom_types_and_radii) {
  // Carbon is hydrophobic unless it is bonded to a nitrogen or an oxygen.
  CHECK(xs_type(6, false, 3, 4) == XsType::C_H);
  CHECK(xs_type(6, true, 2, 4) == XsType::C_P);
  CHECK(xs_hydrophobic(XsType::C_H));
  CHECK(!xs_hydrophobic(XsType::C_P));
  // Nitrogen: a lone pair left over (two bonds or fewer) makes an acceptor, a hydrogen makes a donor.
  CHECK(xs_type(7, false, 0, 2) == XsType::N_A);   // pyridine
  CHECK(xs_type(7, false, 1, 3) == XsType::N_D);   // amide, amidine, pyrrole
  CHECK(xs_type(7, false, 3, 4) == XsType::N_D);   // ammonium
  CHECK(xs_type(7, false, 0, 3) == XsType::N_P);   // nitro, tertiary amine, proline
  CHECK(xs_type(7, false, 1, 2) == XsType::N_DA);  // imine
  CHECK(xs_donor(XsType::N_D) && !xs_acceptor(XsType::N_D));
  CHECK(!xs_donor(XsType::N_A) && xs_acceptor(XsType::N_A));
  CHECK(!xs_donor(XsType::N_P) && !xs_acceptor(XsType::N_P));
  // Oxygen always accepts; a hydrogen makes it donate as well.
  CHECK(xs_type(8, false, 0, 1) == XsType::O_A);
  CHECK(xs_type(8, false, 1, 2) == XsType::O_DA);
  CHECK(xs_acceptor(XsType::O_A) && xs_acceptor(XsType::O_DA));
  CHECK(xs_donor(XsType::O_DA) && !xs_donor(XsType::O_A));
  // Halogens are hydrophobic, metals donate, hydrogen and the elements the function does not know score
  // nothing at all.
  CHECK(xs_type(17, false, 0, 1) == XsType::Cl_H && xs_hydrophobic(XsType::Cl_H));
  CHECK(xs_type(30, false, 0, 0) == XsType::Met_D && xs_donor(XsType::Met_D));
  CHECK(xs_type(1, false, 0, 1) == XsType::None);
  CHECK(xs_type(5, false, 0, 3) == XsType::None);   // boron
  CHECK(xs_type(14, false, 0, 4) == XsType::None);  // silicon
  // X-Score radii.
  CHECK_NEAR(xs_radius(XsType::C_H), 1.9, 1e-12);
  CHECK_NEAR(xs_radius(XsType::N_A), 1.8, 1e-12);
  CHECK_NEAR(xs_radius(XsType::O_A), 1.7, 1e-12);
  CHECK_NEAR(xs_radius(XsType::S_P), 2.0, 1e-12);
  CHECK_NEAR(xs_radius(XsType::P_P), 2.1, 1e-12);
  CHECK_NEAR(xs_radius(XsType::F_H), 1.5, 1e-12);
  CHECK_NEAR(xs_radius(XsType::I_H), 2.2, 1e-12);
  CHECK_NEAR(xs_radius(XsType::Met_D), 1.2, 1e-12);
  CHECK_EQ(xs_radius(XsType::None), 0.0);
}

// --- the ligand's torsion tree ----------------------------------------------------------------

std::size_t torsions_of(const std::string& smiles) {
  Diagnostics diag;
  const auto lig = ligand_from_smiles(smiles, diag);
  if (!lig) throw cltest::Failure("cannot build the ligand " + smiles);
  return lig->rotatable_bonds();
}

TEST(docking_rotatable_bonds) {
  CHECK_EQ(torsions_of("CC"), 0u);                   // both ends have one heavy neighbour
  CHECK_EQ(torsions_of("CCCC"), 1u);                 // only the middle bond
  CHECK_EQ(torsions_of("CCCCC"), 2u);
  CHECK_EQ(torsions_of("Cc1ccccc1"), 0u);            // a methyl turns nothing a heavy atom can see
  CHECK_EQ(torsions_of("c1ccccc1-c1ccccc1"), 1u);    // biphenyl
  CHECK_EQ(torsions_of("C1CCCCC1"), 0u);             // ring bonds are fixed by the embedding
  CHECK_EQ(torsions_of("C1CCCCC1C1CCCCC1"), 1u);     // the bond between two rings turns
  CHECK_EQ(torsions_of("CC(=O)NC"), 0u);             // amide C-N, and both methyls are terminal
  CHECK_EQ(torsions_of("CCC(=O)NCC"), 2u);           // the two C-C bonds next to the amide, not the C-N
  CHECK_EQ(torsions_of("CC#CC"), 0u);                // a triple bond is not single
  CHECK_EQ(torsions_of("N=C(N)c1ccccc1"), 1u);       // benzamidine: only the ring-to-amidine bond
  CHECK_EQ(torsions_of("OCCCCO"), 3u);               // butane-1,4-diol: the three C-C bonds
}

TEST(docking_fragments_and_torsion_atoms) {
  Diagnostics diag;
  const auto lig = ligand_from_smiles("OCCCCO", diag);
  CHECK(lig.has_value());
  CHECK_EQ(lig->rotatable_bonds(), 3u);
  CHECK_EQ(lig->heavy_atoms(), 6u);
  // Four rigid fragments, cut at the three rotatable bonds.
  const std::uint32_t fragments = 1 + *std::max_element(lig->fragment.begin(), lig->fragment.end());
  CHECK_EQ(fragments, 4u);
  // Every torsion turns a proper, non-empty subset of the molecule, and the atom it turns about is in it
  // while the atom it turns around is not.
  for (const LigandTorsion& t : lig->torsions) {
    CHECK(!t.moves.empty());
    CHECK(t.moves.size() < lig->coordinates.size());
    CHECK(std::find(t.moves.begin(), t.moves.end(), t.end) != t.moves.end());
    CHECK(std::find(t.moves.begin(), t.moves.end(), t.begin) == t.moves.end());
  }
}

// --- symmetry-corrected RMSD ------------------------------------------------------------------

PointCloud pose_of(const DockLigand& lig) { return lig.coordinates; }

std::array<double, 3> heavy_centroid(const DockLigand& lig, const PointCloud& pose) {
  std::array<double, 3> c{0, 0, 0};
  for (std::uint32_t i : lig.heavy)
    for (int k = 0; k < 3; ++k) c[k] += pose[i][k];
  for (int k = 0; k < 3; ++k) c[k] /= static_cast<double>(lig.heavy.size());
  return c;
}

// Plain RMSD in the atoms' own order: what the number would be without the isomorphism search.
double rmsd_by_index(const DockLigand& lig, const PointCloud& a, const PointCloud& b) {
  double s = 0;
  for (std::uint32_t i : lig.heavy)
    for (int k = 0; k < 3; ++k) s += (a[i][k] - b[i][k]) * (a[i][k] - b[i][k]);
  return std::sqrt(s / static_cast<double>(lig.heavy.size()));
}

TEST(docking_symmetry_rmsd_sees_through_a_ring_flip) {
  Diagnostics diag;
  const auto lig = ligand_from_smiles("c1ccccc1", diag);
  CHECK(lig.has_value());
  CHECK_EQ(lig->heavy_atoms(), 6u);
  const PointCloud a = pose_of(*lig);
  CHECK_NEAR(symmetry_rmsd(*lig, a, *lig, a), 0.0, 1e-9);

  // Turning benzene by 60 degrees about its own axis maps the ring onto itself: every carbon lands where
  // its neighbour was. By atom order that is a whole C-C bond length of displacement; as a molecule it is
  // the same pose, which is what the isomorphism search has to see.
  const std::array<double, 3> c = heavy_centroid(*lig, a);
  const std::uint32_t p0 = lig->heavy[0], p1 = lig->heavy[1], p2 = lig->heavy[2];
  const std::array<double, 3> u{a[p1][0] - a[p0][0], a[p1][1] - a[p0][1], a[p1][2] - a[p0][2]};
  const std::array<double, 3> v{a[p2][0] - a[p0][0], a[p2][1] - a[p0][1], a[p2][2] - a[p0][2]};
  std::array<double, 3> n{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
  const double nl = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  for (int k = 0; k < 3; ++k) n[k] /= nl;

  PointCloud b = a;
  const double angle = 3.14159265358979323846 / 3.0, cs = std::cos(angle), sn = std::sin(angle);
  for (auto& q : b) {
    const std::array<double, 3> d{q[0] - c[0], q[1] - c[1], q[2] - c[2]};
    const std::array<double, 3> cr{n[1] * d[2] - n[2] * d[1], n[2] * d[0] - n[0] * d[2], n[0] * d[1] - n[1] * d[0]};
    const double dn = n[0] * d[0] + n[1] * d[1] + n[2] * d[2];
    for (int k = 0; k < 3; ++k) q[k] = c[k] + cs * d[k] + sn * cr[k] + (1 - cs) * dn * n[k];
  }
  CHECK(rmsd_by_index(*lig, a, b) > 1.2);  // about one C-C bond
  std::size_t matchings = 0;
  CHECK_NEAR(symmetry_rmsd(*lig, a, *lig, b, &matchings), 0.0, 0.02);
  CHECK(matchings >= 1);  // the search prunes on the best matching so far, so not all twelve are completed

  // A shift of the whole molecule is a real displacement, symmetry or not.
  PointCloud shifted = a;
  for (auto& q : shifted) q[0] += 0.75;
  CHECK_NEAR(symmetry_rmsd(*lig, a, *lig, shifted), 0.75, 1e-6);
}

// A molfile of a ligand with its atoms written in the opposite order: what a crystal ligand from the PDB
// looks like next to the same molecule built from its SMILES.
std::string reversed_molfile(const DockLigand& lig) {
  const std::size_t n = lig.coordinates.size();
  char buf[128];
  std::string out = "reversed\n  ChemLarp          3D\n\n";
  std::snprintf(buf, sizeof buf, "%3zu%3zu  0  0  0  0  0  0  0  0999 V2000\n", n, lig.bonds.size());
  out += buf;
  for (std::size_t k = 0; k < n; ++k) {
    const std::size_t i = n - 1 - k;
    const char* sym = lig.element[i] == 6 ? "C" : (lig.element[i] == 8 ? "O" : "H");
    std::snprintf(buf, sizeof buf, "%10.4f%10.4f%10.4f %-3.3s 0  0  0  0  0  0  0  0  0  0  0  0\n",
                  lig.coordinates[i][0], lig.coordinates[i][1], lig.coordinates[i][2], sym);
    out += buf;
  }
  for (const auto& b : lig.bonds) {
    std::snprintf(buf, sizeof buf, "%3zu%3zu  1  0\n", n - b.first, n - b.second);
    out += buf;
  }
  return out + "M  END\n$$$$\n";
}

TEST(docking_symmetry_rmsd_matches_two_atom_orders) {
  Diagnostics diag;
  const auto lig = ligand_from_smiles("OCCCCO", diag);
  CHECK(lig.has_value());
  const auto other = ligand_from_molfile(reversed_molfile(*lig), diag, 1);
  CHECK(other.has_value());
  CHECK_EQ(other->heavy_atoms(), lig->heavy_atoms());
  // The same coordinates in the opposite order: the RMSD must be zero, which only the isomorphism search
  // can see. The diol is symmetric end to end, so the identity is not the only matching either. The
  // tolerance is the four decimal places a molfile writes coordinates to.
  CHECK_NEAR(symmetry_rmsd(*lig, pose_of(*lig), *other, pose_of(*other)), 0.0, 1e-3);
}

TEST(docking_symmetry_rmsd_rejects_a_different_molecule) {
  Diagnostics diag;
  const auto benzene = ligand_from_smiles("c1ccccc1", diag);
  const auto toluene = ligand_from_smiles("Cc1ccccc1", diag);
  const auto pyridine = ligand_from_smiles("c1ccncc1", diag);
  CHECK(benzene && toluene && pyridine);
  CHECK(symmetry_rmsd(*benzene, pose_of(*benzene), *toluene, pose_of(*toluene)) < 0);   // different sizes
  CHECK(symmetry_rmsd(*benzene, pose_of(*benzene), *pyridine, pose_of(*pyridine)) < 0);  // different elements
}

// --- a synthetic receptor ---------------------------------------------------------------------

ProteinAtom make_atom(const char* name, const char* element, double x, double y, double z) {
  ProteinAtom a;
  a.name = name;
  a.element = element;
  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

// Two parallel sheets of carbon 8 A apart: a flat hydrophobic slot whose best pose for an aromatic ring is
// the ring lying in the middle plane, which gives the redock below something unambiguous to find.
Protein hydrophobic_slot() {
  Protein p;
  ProteinResidue r;
  r.name = "PKT";
  r.chain = "A";
  r.number = 1;
  r.hetero = true;
  int n = 0;
  for (int sign = -1; sign <= 1; sign += 2)
    for (int ix = -2; ix <= 2; ++ix)
      for (int iy = -2; iy <= 2; ++iy)
        r.atoms.push_back(make_atom(("C" + std::to_string(++n)).c_str(), "C", 2.0 * ix, 2.0 * iy, 4.0 * sign));
  p.residues.push_back(std::move(r));
  return p;
}

ReceptorOptions bare_receptor_options() {
  ReceptorOptions o;
  o.protonate = false;  // the hydrogens are placed by hand here, not predicted
  o.repair = false;
  return o;
}

TEST(docking_receptor_typing_from_geometry) {
  Protein p;
  ProteinResidue r;
  r.name = "TST";
  r.chain = "A";
  r.number = 1;
  r.hetero = true;
  // A nitrogen with two heavy neighbours and a hydrogen: a donor that cannot accept, like a backbone N-H.
  r.atoms.push_back(make_atom("N1", "N", 0, 0, 0));
  r.atoms.push_back(make_atom("H1", "H", 0, 0, 1.01));
  r.atoms.push_back(make_atom("CA", "C", 1.45, 0, 0));
  r.atoms.push_back(make_atom("CB", "C", -1.45, 0, 0));
  // A nitrogen with two heavy neighbours and no hydrogen: an acceptor, like an imidazole N.
  r.atoms.push_back(make_atom("N2", "N", 10, 0, 0));
  r.atoms.push_back(make_atom("CC", "C", 11.45, 0, 0));
  r.atoms.push_back(make_atom("CD", "C", 8.55, 0, 0));
  // A carbonyl: the carbon turns polar, the oxygen accepts.
  r.atoms.push_back(make_atom("CE", "C", 20, 0, 0));
  r.atoms.push_back(make_atom("O1", "O", 21.23, 0, 0));
  // A hydroxyl oxygen: donor and acceptor.
  r.atoms.push_back(make_atom("O2", "O", 25, 0, 0));
  r.atoms.push_back(make_atom("H2", "H", 25.96, 0, 0));
  // A lone carbon stays hydrophobic, and a zinc ion counts as a donor.
  r.atoms.push_back(make_atom("CF", "C", 30, 0, 0));
  r.atoms.push_back(make_atom("ZN", "Zn", 40, 0, 0));
  p.residues.push_back(std::move(r));

  Diagnostics diag;
  const DockReceptor rec = prepare_receptor(p, bare_receptor_options(), {}, diag);
  auto type_of = [&](const std::string& name) {
    for (const ReceptorAtom& a : rec.atoms)
      if (a.name == name) return a.type;
    return XsType::None;
  };
  CHECK(type_of("N1") == XsType::N_D);
  CHECK(type_of("N2") == XsType::N_A);
  CHECK(type_of("CE") == XsType::C_P);
  CHECK(type_of("O1") == XsType::O_A);
  CHECK(type_of("O2") == XsType::O_DA);
  CHECK(type_of("CF") == XsType::C_H);
  CHECK(type_of("ZN") == XsType::Met_D);
  CHECK(type_of("H1") == XsType::None);  // hydrogens are not scored
  CHECK_EQ(rec.atoms.size(), 11u);       // the twelve heavy atoms, less nothing; the two hydrogens are gone
  CHECK(rec.kept_hetero.size() == 1 && rec.kept_hetero[0] == "A 1 TST");
}

TEST(docking_receptor_drops_waters_and_the_reference_site) {
  Protein p = hydrophobic_slot();
  ProteinResidue w;
  w.name = "HOH";
  w.chain = "A";
  w.number = 2;
  w.hetero = true;
  w.atoms.push_back(make_atom("O", "O", 0, 0, 0));
  p.residues.push_back(w);
  ProteinResidue other;  // a hetero group sitting in the site the reference ligand occupies
  other.name = "LIG";
  other.chain = "A";
  other.number = 3;
  other.hetero = true;
  other.atoms.push_back(make_atom("C1", "C", 0.5, 0.5, 0.0));
  p.residues.push_back(other);

  Diagnostics diag;
  ReceptorOptions dry = bare_receptor_options();
  dry.remove_waters = true;
  const DockReceptor with_reference = prepare_receptor(p, dry, {{0.0, 0.0, 0.0}}, diag);
  CHECK_EQ(with_reference.waters_removed, 1u);
  CHECK_EQ(with_reference.hetero_removed, 1u);
  CHECK_EQ(with_reference.atoms.size(), 50u);  // the slot alone
  // Without a reference there is nothing to clear, so the other group stays; the water still goes.
  const DockReceptor without = prepare_receptor(p, dry, {}, diag);
  CHECK_EQ(without.waters_removed, 1u);
  CHECK_EQ(without.hetero_removed, 0u);
  CHECK_EQ(without.atoms.size(), 51u);
  // By default the waters stay, because a water a structure resolves is part of the site — but one the
  // reference ligand sits on top of is a water the ligand displaces, and it goes with the rest of the site.
  const DockReceptor wet = prepare_receptor(p, bare_receptor_options(), {{0.0, 0.0, 0.0}}, diag);
  CHECK_EQ(wet.waters_removed, 0u);
  CHECK_EQ(wet.hetero_removed, 2u);  // the water on the reference ligand, and the other group
  CHECK_EQ(wet.atoms.size(), 50u);
  // A water away from the ligand stays.
  Protein far = hydrophobic_slot();
  ProteinResidue w2;
  w2.name = "HOH";
  w2.chain = "A";
  w2.number = 9;
  w2.hetero = true;
  w2.atoms.push_back(make_atom("O", "O", 0, 0, 20));
  far.residues.push_back(w2);
  const DockReceptor kept = prepare_receptor(far, bare_receptor_options(), {{0.0, 0.0, 0.0}}, diag);
  CHECK_EQ(kept.waters_removed, 0u);
  CHECK_EQ(kept.hetero_removed, 0u);
  CHECK_EQ(kept.atoms.size(), 51u);  // the slot and the distant water's oxygen
}

// A residue written on HETATM records but bonded into the chain (a modified residue such as
// selenomethionine or S-hydroxycysteine) is part of the protein and must survive the clearing of the
// ligand's site, even though a ligand can carry backbone atom names of its own.
TEST(docking_receptor_keeps_a_modified_residue_in_the_chain) {
  Protein p;
  auto residue = [](const char* name, int number, bool hetero, double x) {
    ProteinResidue r;
    r.name = name;
    r.chain = "A";
    r.number = number;
    r.hetero = hetero;
    r.atoms.push_back(make_atom("N", "N", x, 0, 0));
    r.atoms.push_back(make_atom("CA", "C", x + 1.46, 0, 0));
    r.atoms.push_back(make_atom("C", "C", x + 2.43, 0, 0));
    return r;
  };
  // Three residues in a row, 3.4 A apart, so each C is 1.33 A from the next N: a peptide chain. The middle
  // one is on HETATM records, like a modified residue.
  p.residues.push_back(residue("ALA", 1, false, 0.0));
  p.residues.push_back(residue("CSO", 2, true, 3.76));
  p.residues.push_back(residue("ALA", 3, false, 7.52));
  // A real ligand, not bonded to anything, which does carry backbone atom names.
  ProteinResidue ligand = residue("MTX", 4, true, 40.0);
  p.residues.push_back(ligand);

  Diagnostics diag;
  // A reference ligand on top of both hetero groups: only the one that is not in the chain may go.
  const DockReceptor rec =
      prepare_receptor(p, bare_receptor_options(), {{3.76, 0.0, 0.0}, {40.0, 0.0, 0.0}}, diag);
  CHECK_EQ(rec.hetero_removed, 1u);
  CHECK_EQ(rec.atoms.size(), 9u);  // three residues of three atoms; the ligand is gone
  CHECK(rec.kept_hetero.empty());  // the modified residue is chain, not a kept hetero group
}

// --- the pose gradient ------------------------------------------------------------------------

TEST(docking_pose_gradient_matches_finite_differences) {
  Diagnostics diag;
  const Protein p = hydrophobic_slot();
  const DockReceptor rec = prepare_receptor(p, bare_receptor_options(), {}, diag);
  const auto lig = ligand_from_smiles("OCCCCO", diag);
  CHECK(lig.has_value());
  const std::size_t n = pose_variable_count(*lig);
  CHECK_EQ(n, 9u);  // three translations, three for the rotation vector, three torsions

  // Three poses off every special point: away from the origin, with a rotation vector of real length, and
  // with torsions nowhere near zero. The rotation's exponential map has to be differentiated exactly, not
  // only at the identity.
  for (int trial = 0; trial < 3; ++trial) {
    const double s = 0.4 * (trial + 1);
    std::vector<double> v(n, 0.0);
    v[0] = 0.4 * s;
    v[1] = -0.3 * s;
    v[2] = 0.2 * s;
    v[3] = 0.5 * s;
    v[4] = -0.8 * s;
    v[5] = 0.3 * s;
    for (std::size_t k = 6; k < n; ++k) v[k] = s * (0.4 + 0.3 * static_cast<double>(k));
    std::vector<double> g;
    pose_variable_energy(rec, *lig, v, {}, &g);
    CHECK_EQ(g.size(), n);
    for (std::size_t k = 0; k < n; ++k) {
      const double h = 1e-6;
      std::vector<double> a = v, b = v;
      a[k] += h;
      b[k] -= h;
      const double numeric =
          (pose_variable_energy(rec, *lig, a) - pose_variable_energy(rec, *lig, b)) / (2 * h);
      CHECK_NEAR(g[k], numeric, 1e-4 * std::max(1.0, std::fabs(numeric)));
    }
  }
}

TEST(docking_pose_variables_reproduce_the_input_pose) {
  Diagnostics diag;
  const Protein p = hydrophobic_slot();
  const DockReceptor rec = prepare_receptor(p, bare_receptor_options(), {}, diag);
  const auto lig = ligand_from_smiles("OCCCCO", diag);
  CHECK(lig.has_value());
  // The variables are relative to the heavy-atom centroid, so putting the centroid back where it was, with
  // no rotation and no torsion, has to give the coordinates the ligand came with.
  const std::array<double, 3> c = heavy_centroid(*lig, pose_of(*lig));
  std::vector<double> v(pose_variable_count(*lig), 0.0);
  v[0] = c[0];
  v[1] = c[1];
  v[2] = c[2];
  PointCloud back;
  const double e = pose_variable_energy(rec, *lig, v, {}, nullptr, &back);
  CHECK_EQ(back.size(), lig->coordinates.size());
  for (std::size_t i = 0; i < back.size(); ++i)
    for (int k = 0; k < 3; ++k) CHECK_NEAR(back[i][k], lig->coordinates[i][k], 1e-9);
  // And the energy is the one score_pose gives for those coordinates (the ligand's own pairs included).
  const auto scored = score_pose(rec, *lig, lig->coordinates);
  CHECK(scored.has_value());
  CHECK_NEAR(e, scored->intermolecular + scored->internal, 1e-9);
  // A wrong number of variables is rejected rather than read past the end.
  std::vector<double> short_vector(3, 0.0);
  CHECK_EQ(pose_variable_energy(rec, *lig, short_vector), 0.0);
}

// --- rescoring and refining a pose ------------------------------------------------------------

TEST(docking_rescores_a_pose_without_moving_it) {
  Diagnostics diag;
  const Protein p = hydrophobic_slot();
  const DockReceptor rec = prepare_receptor(p, bare_receptor_options(), {}, diag);
  const auto lig = ligand_from_smiles("c1ccccc1", diag);
  CHECK(lig.has_value());

  // The ring flat in the middle of the slot: every carbon is about 4 A from the sheets, which is the
  // hydrophobic plateau, so the score has to be clearly negative.
  PointCloud centred = lig->coordinates;
  const std::array<double, 3> c = heavy_centroid(*lig, centred);
  for (auto& q : centred)
    for (int k = 0; k < 3; ++k) q[k] -= c[k];
  const auto in_slot = score_pose(rec, *lig, centred);
  CHECK(in_slot.has_value());
  CHECK(in_slot->intermolecular < -3.0);
  CHECK(in_slot->terms.hydrophobic > 10.0);  // many hydrophobic contacts
  CHECK_EQ(in_slot->terms.hbond, 0.0);       // nothing here donates or accepts
  CHECK_EQ(in_slot->internal, 0.0);          // benzene is one rigid fragment
  CHECK_NEAR(in_slot->affinity, in_slot->intermolecular, 1e-12);  // no rotatable bonds, so no N_rot factor

  // Far away there is nothing to interact with.
  PointCloud away = centred;
  for (auto& q : away) q[0] += 60.0;
  const auto outside = score_pose(rec, *lig, away);
  CHECK(outside.has_value());
  CHECK_EQ(outside->intermolecular, 0.0);
  // Jammed against a sheet it is repulsive.
  PointCloud clash = centred;
  for (auto& q : clash) q[2] += 4.0;
  const auto overlapping = score_pose(rec, *lig, clash);
  CHECK(overlapping.has_value());
  CHECK(overlapping->terms.repulsion > 1.0);
  CHECK(overlapping->intermolecular > in_slot->intermolecular);
  // A pose of the wrong size is refused.
  CHECK(!score_pose(rec, *lig, PointCloud{{0, 0, 0}}).has_value());
}

TEST(docking_refining_a_pose_keeps_the_molecule_rigid_and_lowers_the_score) {
  Diagnostics diag;
  const Protein p = hydrophobic_slot();
  const DockReceptor rec = prepare_receptor(p, bare_receptor_options(), {}, diag);
  auto lig = ligand_from_smiles("OCCCCO", diag);
  CHECK(lig.has_value());
  // Put the ligand in the slot, tilted, and refine it where it stands.
  const std::array<double, 3> c = heavy_centroid(*lig, lig->coordinates);
  for (auto& q : lig->coordinates)
    for (int k = 0; k < 3; ++k) q[k] -= c[k];
  const auto before = score_pose(rec, *lig, lig->coordinates);
  CHECK(before.has_value());

  DockOptions o;
  o.box = DockBox{0, 0, 0, 16, 16, 16};
  o.refine_only = true;
  o.modes = 1;
  const auto result = dock(rec, *lig, o, diag);
  CHECK(result.has_value());
  CHECK_EQ(result->poses.size(), 1u);
  CHECK_EQ(result->steps, 0);
  const DockedPose& pose = result->poses.front();
  CHECK(pose.intermolecular + pose.internal <= before->intermolecular + before->internal + 1e-9);

  // Only the position, the orientation and the torsions may change, so every bond length is untouched and
  // so is every distance inside a rigid fragment.
  for (const auto& [a, b] : lig->bonds) {
    auto d = [](const std::array<double, 3>& u, const std::array<double, 3>& v) {
      return std::sqrt((u[0] - v[0]) * (u[0] - v[0]) + (u[1] - v[1]) * (u[1] - v[1]) + (u[2] - v[2]) * (u[2] - v[2]));
    };
    CHECK_NEAR(d(pose.coordinates[a], pose.coordinates[b]), d(lig->coordinates[a], lig->coordinates[b]), 1e-7);
  }
}

// --- a redock ---------------------------------------------------------------------------------

TEST(docking_finds_the_pocket_in_a_synthetic_receptor) {
  Diagnostics diag;
  const Protein p = hydrophobic_slot();
  const DockReceptor rec = prepare_receptor(p, bare_receptor_options(), {}, diag);
  CHECK_EQ(rec.atoms.size(), 50u);
  const auto lig = ligand_from_smiles("c1ccccc1", diag);
  CHECK(lig.has_value());

  DockOptions o;
  o.box = DockBox{0, 0, 0, 14, 14, 14};
  o.runs = 4;
  o.steps = 30;
  o.modes = 3;
  o.seed = 7;
  const auto result = dock(rec, *lig, o, diag);
  CHECK(result.has_value());
  CHECK(!result->poses.empty());
  CHECK_EQ(result->rotatable_bonds, 0u);
  CHECK_EQ(result->heavy_atoms, 6u);
  CHECK(result->evaluations > 0);

  // The slot's best pose is the ring lying in the middle plane: the centroid near the origin and the ring's
  // normal along z. Nothing in the run knows that; it comes out of the search.
  const DockedPose& best = result->poses.front();
  CHECK(best.affinity < -3.0);
  const std::array<double, 3> c = heavy_centroid(*lig, best.coordinates);
  CHECK(std::fabs(c[0]) < 2.0);
  CHECK(std::fabs(c[1]) < 2.0);
  CHECK(std::fabs(c[2]) < 1.0);
  const std::uint32_t a0 = lig->heavy[0], a1 = lig->heavy[1], a2 = lig->heavy[2];
  const std::array<double, 3> u{best.coordinates[a1][0] - best.coordinates[a0][0],
                               best.coordinates[a1][1] - best.coordinates[a0][1],
                               best.coordinates[a1][2] - best.coordinates[a0][2]};
  const std::array<double, 3> v{best.coordinates[a2][0] - best.coordinates[a0][0],
                               best.coordinates[a2][1] - best.coordinates[a0][1],
                               best.coordinates[a2][2] - best.coordinates[a0][2]};
  std::array<double, 3> nrm{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
  const double nl = std::sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
  CHECK(std::fabs(nrm[2] / nl) > 0.85);

  // The poses reported are distinct and ordered.
  for (std::size_t k = 1; k < result->poses.size(); ++k)
    CHECK(symmetry_rmsd(*lig, result->poses[k].coordinates, *lig, result->poses[0].coordinates) >= o.min_rmsd ||
          result->poses[k].affinity >= result->poses[0].affinity);

  // The same seed gives the same answer, whatever the thread count sees.
  const auto again = dock(rec, *lig, o, diag);
  CHECK(again.has_value());
  CHECK_EQ(again->poses.size(), result->poses.size());
  CHECK_NEAR(again->poses.front().affinity, best.affinity, 1e-12);

  // The poses come out as an SDF with the affinity attached.
  const std::string sdf = write_poses_sdf(*lig, result->poses);
  CHECK(sdf.find("V2000") != std::string::npos);
  CHECK(sdf.find("> <affinity_kcal_mol>") != std::string::npos);
  CHECK_EQ(static_cast<std::size_t>(std::count(sdf.begin(), sdf.end(), '$') / 4), result->poses.size());
}

// A box far larger than any site would ask for gigabytes of grid maps at the default spacing. The spacing
// widens instead, so the job stays inside memory and the result says what it was built at.
TEST(docking_coarsens_the_grid_for_an_enormous_box) {
  Diagnostics diag;
  const Protein p = hydrophobic_slot();
  const DockReceptor rec = prepare_receptor(p, bare_receptor_options(), {}, diag);
  const auto lig = ligand_from_smiles("c1ccccc1", diag);
  CHECK(lig.has_value());
  DockOptions o;
  o.box = DockBox{0, 0, 0, 40, 40, 40};  // 44 A of grid a side: 118^3 points, well inside the limit
  o.runs = 1;
  o.steps = 2;
  o.modes = 1;
  const auto small = dock(rec, *lig, o, diag);
  CHECK(small.has_value());
  CHECK_NEAR(small->grid_spacing, o.grid_spacing, 1e-12);

  o.box = DockBox{0, 0, 0, 180, 180, 180};  // 490 a side at 0.375 A: 1.2e8 points, about a gigabyte a map
  const auto huge = dock(rec, *lig, o, diag);
  CHECK(huge.has_value());
  CHECK(huge->grid_spacing > o.grid_spacing);
  CHECK(!huge->poses.empty());
}

TEST(docking_rejects_impossible_requests) {
  Diagnostics diag;
  const Protein p = hydrophobic_slot();
  const DockReceptor rec = prepare_receptor(p, bare_receptor_options(), {}, diag);
  const auto lig = ligand_from_smiles("c1ccccc1", diag);
  CHECK(lig.has_value());
  DockOptions o;
  o.runs = 1;
  o.steps = 2;
  // No box.
  CHECK(!dock(rec, *lig, o, diag).has_value());
  o.box = DockBox{0, 0, 0, 10, 10, 10};
  // No runs, no poses.
  DockOptions zero_runs = o;
  zero_runs.runs = 0;
  CHECK(!dock(rec, *lig, zero_runs, diag).has_value());
  DockOptions zero_modes = o;
  zero_modes.modes = 0;
  CHECK(!dock(rec, *lig, zero_modes, diag).has_value());
  // An empty receptor.
  CHECK(!dock(DockReceptor{}, *lig, o, diag).has_value());
  // A ligand that is not a molecule.
  Diagnostics bad;
  CHECK(!ligand_from_smiles("C(C(C", bad).has_value());
  CHECK(bad.has_errors());
  Diagnostics flat;
  CHECK(!ligand_from_molfile("not a molfile at all\n\n\n\n", flat).has_value());
}

// The whole job on a real complex, with the one structure the repository carries (data/dock, PDB entry
// 3PTB: bovine beta-trypsin with benzamidine). The deposited ligand gives only its chemical identity — its
// canonical SMILES — and a conformation is generated from scratch, so the crystal geometry is the reference
// for the RMSD and never an input to the search. docs/ARCHITECTURE.md reports the other eleven complexes of
// the set, which are not in the repository.
TEST(docking_redocks_benzamidine_into_trypsin) {
  const std::string pdb = read_data_file(cltest::data_path("dock/3PTB.pdb"));
  const std::string sdf = read_data_file(cltest::data_path("dock/3PTB_BEN.sdf"));
  CHECK(pdb.size() > 100000);
  CHECK(sdf.size() > 500);

  Diagnostics diag;
  const std::optional<Protein> protein = read_protein(pdb, diag);
  CHECK(protein.has_value());
  const auto crystal = ligand_from_molfile(sdf, diag);
  CHECK(crystal.has_value());
  CHECK_EQ(crystal->smiles, std::string("N=C(N)c1ccccc1"));
  CHECK_EQ(crystal->heavy_atoms(), 9u);

  PointCloud site;
  for (std::uint32_t i : crystal->heavy) site.push_back(crystal->coordinates[i]);
  const DockReceptor rec = prepare_receptor(*protein, ReceptorOptions{}, site, diag);
  // The ligand's own copy is cleared from the receptor; the calcium ion and the waters away from the site
  // stay, since a resolved water is part of the structure.
  CHECK_EQ(rec.hetero_removed, 1u);
  CHECK_EQ(rec.waters_removed, 0u);
  CHECK(rec.atoms.size() > 1600);
  CHECK(rec.donors > 100 && rec.acceptors > 100 && rec.hydrophobic > 100);

  // The molecule rebuilt from its SMILES, which is what gets docked.
  const auto ligand = ligand_from_smiles(crystal->smiles, diag);
  CHECK(ligand.has_value());
  CHECK_EQ(ligand->heavy_atoms(), 9u);
  CHECK_EQ(ligand->rotatable_bonds(), 1u);  // only the ring-to-amidine bond

  DockOptions o;
  o.box = box_around(site, 4.0);
  o.runs = 4;
  o.seed = 42;
  const auto result = dock(rec, *ligand, o, diag, &*crystal);
  CHECK(result.has_value());
  CHECK(!result->poses.empty());
  CHECK(result->has_reference);

  const DockedPose& best = result->poses.front();
  CHECK(best.rmsd >= 0);
  CHECK(best.rmsd < 2.0);          // the crystal pose, reproduced
  CHECK(best.affinity < -4.0);     // benzamidine in trypsin scores about -6 kcal/mol
  CHECK(best.affinity > -12.0);
  // The deposited pose, minimised where it stands, stays where it is and scores in the same range: that is
  // the yardstick that tells a search failure from a scoring failure.
  CHECK(result->reference_rmsd >= 0 && result->reference_rmsd < 1.0);
  CHECK(result->reference_affinity < -4.0);
  CHECK(!diag.has_errors());
}

TEST(docking_json_reports_a_run_and_refuses_bad_input) {
  // A receptor as a PDB file, so the JSON entry point is exercised end to end.
  std::string pdb;
  char buf[96];
  int serial = 0;
  for (int sign = -1; sign <= 1; sign += 2)
    for (int ix = -2; ix <= 2; ++ix)
      for (int iy = -2; iy <= 2; ++iy) {
        ++serial;
        std::snprintf(buf, sizeof buf, "HETATM%5d  C%-2d PKT A   1    %8.3f%8.3f%8.3f  1.00  0.00           C\n",
                      serial, serial % 100, 2.0 * ix, 2.0 * iy, 4.0 * sign);
        pdb += buf;
      }
  pdb += "END\n";

  DockJsonOptions o;
  o.receptor = ReceptorOptions{};
  o.receptor.protonate = false;
  o.dock.box = DockBox{0, 0, 0, 14, 14, 14};
  o.dock.runs = 2;
  o.dock.steps = 10;
  o.dock.modes = 2;
  o.dock.seed = 3;
  const std::string json = dock_json(pdb, "c1ccccc1", o);
  CHECK(json.rfind("{\"ok\":true", 0) == 0);
  CHECK(json.find("\"smiles\":\"c1ccccc1\"") != std::string::npos);
  CHECK(json.find("\"heavy_atoms\":6") != std::string::npos);
  CHECK(json.find("\"rotatable_bonds\":0") != std::string::npos);
  CHECK(json.find("\"affinity\":-") != std::string::npos);
  CHECK(json.find("\"sdf\":") != std::string::npos);

  // No box and no reference: there is nothing to say where to search.
  DockJsonOptions no_box = o;
  no_box.dock.box = DockBox{};
  const std::string refused = dock_json(pdb, "c1ccccc1", no_box);
  CHECK(refused.rfind("{\"ok\":false", 0) == 0);
  CHECK(refused.find("InvalidInput") != std::string::npos);
  // A reference residue that is not in the file.
  DockJsonOptions missing = no_box;
  missing.reference = "XYZ";
  CHECK(dock_json(pdb, "c1ccccc1", missing).find("InvalidInput") != std::string::npos);
  // A receptor that is not a structure file, and a ligand that is not a molecule.
  CHECK(dock_json("hello", "c1ccccc1", o).rfind("{\"ok\":false", 0) == 0);
  CHECK(dock_json(pdb, "C(C(C", o).rfind("{\"ok\":false", 0) == 0);
}

}  // namespace
