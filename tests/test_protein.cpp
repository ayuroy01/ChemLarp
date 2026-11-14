// Proteins: PDB and mmCIF reading, disulfides, missing atoms and loops, pKa, protonation, hydrogen-bond
// network, hydrogen placement, sequence properties.
//
// data/1crn.pdb is crambin (PDB 1CRN, X-ray, no hydrogens); data/1l2y_model1.pdb is the first model of the
// Trp-cage NMR ensemble (PDB 1L2Y) with its deposited hydrogens. Reference pI and net charges are from
// Biopython's IsoelectricPoint (same Bjellqvist pKa set); see "Proteins: validation" in docs/ARCHITECTURE.md.
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <tuple>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include "chemlarp/chemlarp_c.h"
#include "chemlarp/protein.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

std::string slurp(const std::string& name) {
  std::ifstream f(cltest::data_path(name));
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

double distance(const ProteinAtom& a, const ProteinAtom& b) {
  return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

}  // namespace

TEST(protein_reads_crambin) {
  Diagnostics d;
  auto p = read_protein_pdb(slurp("1crn.pdb"), d);
  CHECK(p.has_value());
  const auto seqs = chain_sequences(*p);
  CHECK_EQ(seqs.size(), std::size_t{1});
  CHECK_EQ(seqs[0].second, std::string("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"));
  CHECK_EQ(p->atom_count(), std::size_t{327});
  CHECK(p->chain_breaks.empty());
  CHECK(p->residues.front().n_terminal && p->residues.back().c_terminal);
  // Disulfides Cys3–Cys40, Cys4–Cys32, Cys16–Cys26.
  CHECK_EQ(p->disulfides.size(), std::size_t{3});
  std::map<int, int> ss;
  for (const auto& [a, b] : p->disulfides) ss[p->residues[a].number] = p->residues[b].number;
  CHECK_EQ(ss[3], 40);
  CHECK_EQ(ss[4], 32);
  CHECK_EQ(ss[16], 26);
}

TEST(protein_protonation_states) {
  for (const auto& [pH, net, tym] : {std::tuple{2.0, 3, 0}, std::tuple{7.0, 0, 0}, std::tuple{12.0, -3, 2}}) {
    Diagnostics d;
    auto p = read_protein_pdb(slurp("1crn.pdb"), d);
    CHECK(p.has_value());
    ProtonationOptions o;
    o.pH = pH;
    const auto rep = protonate(*p, o, d);
    CHECK_EQ(rep.net_charge, net);
    int n_tym = 0, n_cyx = 0;
    for (const auto& r : p->residues) n_tym += r.state == "TYM", n_cyx += r.state == "CYX";
    CHECK_EQ(n_tym, tym);
    CHECK_EQ(n_cyx, 6);
    CHECK(rep.missing_atoms.empty());
    // Every disulfide cysteine stays without HG; the C-terminal residue gains OXT.
    for (const auto& r : p->residues)
      if (r.state == "CYX") CHECK(!r.atom("HG"));
    CHECK(p->residues.back().atom("OXT") != nullptr);
    CHECK_EQ(p->residues.back().atom("HXT") != nullptr, pH < 3.2);
  }
}

TEST(protein_hydrogens_match_nmr_structure) {
  // Rebuild all hydrogens of the Trp-cage NMR model and compare with the deposited ones by atom name.
  const std::string text = slurp("1l2y_model1.pdb");
  Diagnostics d;
  auto ref = read_protein_pdb(text, d);
  auto p = read_protein_pdb(text, d);
  CHECK(ref.has_value() && p.has_value());
  const auto rep = protonate(*p, ProtonationOptions{}, d);
  CHECK_EQ(rep.hydrogens_removed, std::size_t{150});
  CHECK_EQ(rep.hydrogens_added, std::size_t{150});
  double sum2 = 0.0;
  int n = 0;
  for (std::size_t i = 0; i < ref->residues.size(); ++i) {
    const auto& rr = ref->residues[i];
    const auto& pr = p->residues[i];
    for (const auto& h : rr.atoms) {
      if (h.element != "H") continue;
      // Hydroxyl rotors are placed by a hydrogen-bond search, not compared.
      if ((rr.name == "SER" && h.name == "HG") || (rr.name == "TYR" && h.name == "HH")) continue;
      // NH3+ groups rotate in the hydrogen-bond network as well.
      if ((rr.name == "LYS" && h.name.rfind("HZ", 0) == 0) || (rr.n_terminal && (h.name == "H1" || h.name == "H2" || h.name == "H3"))) continue;
      // Methyl and NH3+ hydrogens: compare with the nearest of the group (rotation, not naming).
      double best = 1e9;
      const std::string stem = h.name.substr(0, h.name.size() - 1);
      const bool group = std::isdigit(static_cast<unsigned char>(h.name.back())) && pr.atom(stem + "1") && pr.atom(stem + "3");
      for (const auto& q : pr.atoms) {
        if (q.element != "H") continue;
        if (group ? q.name.substr(0, q.name.size() - 1) == h.name.substr(0, h.name.size() - 1) : q.name == h.name)
          best = std::min(best, distance(h, q));
      }
      CHECK(best < 0.3);
      sum2 += best * best;
      ++n;
    }
  }
  CHECK(n > 135);
  CHECK(std::sqrt(sum2 / n) < 0.12);
}

TEST(protein_write_read_round_trip) {
  Diagnostics d;
  auto p = read_protein_pdb(slurp("1crn.pdb"), d);
  CHECK(p.has_value());
  protonate(*p, ProtonationOptions{}, d);
  const std::string out = write_protein_pdb(*p);
  auto q = read_protein_pdb(out, d);
  CHECK(q.has_value());
  CHECK_EQ(q->atom_count(), p->atom_count());
  CHECK_EQ(q->disulfides.size(), std::size_t{3});
  // Re-protonating the written structure reproduces it.
  const auto rep = protonate(*q, ProtonationOptions{}, d);
  CHECK_EQ(rep.hydrogens_removed, rep.hydrogens_added);
  for (std::size_t i = 0; i < p->residues.size(); ++i)
    for (const auto& a : p->residues[i].atoms) {
      const ProteinAtom* b = q->residues[i].atom(a.name);
      CHECK(b != nullptr);
      if (b) CHECK(distance(a, *b) < 2e-3);
    }
}

TEST(protein_sequence_properties) {
  // Biopython IsoelectricPoint / ProtParam references.
  const std::string all20 = "ACDEFGHIKLMNPQRSTVWY";
  CHECK_NEAR(isoelectric_point(all20), 6.7846, 2e-4);
  CHECK_NEAR(sequence_net_charge(all20, 7.0), -0.12490, 1e-4);
  CHECK_NEAR(sequence_net_charge(all20, 4.0), 2.51810, 1e-4);
  CHECK_NEAR(isoelectric_point("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"), 5.7252, 2e-4);
  CHECK_NEAR(isoelectric_point("EKKK"), 9.7053, 2e-4);
  // Gly-Gly is C4H8N2O3: 132.119 g/mol with standard atomic weights.
  CHECK_NEAR(sequence_molecular_weight("GG"), 4 * 12.011 + 8 * 1.008 + 2 * 14.007 + 3 * 15.999, 1e-9);
  CHECK_EQ(sequence_molecular_weight(""), 0.0);
  // PROPKA model set: one Lys peptide at pH 7 is almost +1 (N-term + Lys − C-term).
  CHECK_NEAR(sequence_net_charge("GKG", 7.0, PkaSet::PropkaModel), 1.0, 0.1);
}

TEST(protein_errors_and_c_api) {
  Diagnostics d;
  CHECK(!read_protein_pdb("HEADER nothing here\n", d).has_value());
  CHECK(d.has_errors());
  char* s = cl_protein_json(slurp("1crn.pdb").c_str(), 7.0);
  const std::string j(s);
  cl_free_string(s);
  CHECK(j.rfind("{\"ok\":true", 0) == 0);
  CHECK(j.find("\"disulfides\":[[\"A\",3,\"A\",40]") != std::string::npos);
  CHECK(j.find("\"hydrogens_added\":315") != std::string::npos);
  char* bad = cl_protein_json("ATOM", -1.0);
  CHECK(std::string(bad).find("\"ok\":false") != std::string::npos);
  cl_free_string(bad);
}

namespace {

bool backbone_name(const std::string& n) { return n == "N" || n == "CA" || n == "C" || n == "O" || n == "OXT"; }

double dihedral(const ProteinAtom& a, const ProteinAtom& b, const ProteinAtom& c, const ProteinAtom& d) {
  const double b0[3] = {a.x - b.x, a.y - b.y, a.z - b.z}, b1r[3] = {c.x - b.x, c.y - b.y, c.z - b.z}, b2[3] = {d.x - c.x, d.y - c.y, d.z - c.z};
  const double n1 = std::sqrt(b1r[0] * b1r[0] + b1r[1] * b1r[1] + b1r[2] * b1r[2]);
  const double b1[3] = {b1r[0] / n1, b1r[1] / n1, b1r[2] / n1};
  auto dot = [](const double* u, const double* v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
  double v[3], w[3];
  for (int k = 0; k < 3; ++k) v[k] = b0[k] - dot(b0, b1) * b1[k], w[k] = b2[k] - dot(b2, b1) * b1[k];
  const double x[3] = {b1[1] * v[2] - b1[2] * v[1], b1[2] * v[0] - b1[0] * v[2], b1[0] * v[1] - b1[1] * v[0]};
  return std::atan2(dot(x, w), dot(v, w)) * 180.0 / 3.14159265358979323846;
}

double angle_difference(double a, double b) {
  double d = std::fmod(a - b, 360.0);
  if (d > 180.0) d -= 360.0;
  if (d < -180.0) d += 360.0;
  return std::fabs(d);
}

// His 15 of hen lysozyme (PDB 1AKI) as a lone residue, with one extra HETATM line.
std::string histidine_with(const std::string& extra) {
  return "ATOM      1  N   HIS A  15      39.681  23.247   5.463  1.00 17.84           N\n"
         "ATOM      2  CA  HIS A  15      39.075  24.591   5.526  1.00 15.99           C\n"
         "ATOM      3  C   HIS A  15      38.310  24.861   6.814  1.00 17.70           C\n"
         "ATOM      4  O   HIS A  15      37.608  25.875   6.952  1.00 22.68           O\n"
         "ATOM      5  CB  HIS A  15      38.223  24.918   4.285  1.00 19.19           C\n"
         "ATOM      6  CG  HIS A  15      39.085  25.388   3.171  1.00 20.14           C\n"
         "ATOM      7  ND1 HIS A  15      39.457  24.635   2.091  1.00 24.53           N\n"
         "ATOM      8  CD2 HIS A  15      39.739  26.570   2.993  1.00 24.44           C\n"
         "ATOM      9  CE1 HIS A  15      40.211  25.312   1.258  1.00 20.68           C\n"
         "ATOM     10  NE2 HIS A  15      40.524  26.419   1.889  1.00 27.34           N\n" +
         extra + "END\n";
}

// A HETATM line at distance d from `at`, along the bisector pointing away from its two ring neighbours.
std::string probe(const std::string& text, const char* at, const char* n1, const char* n2, double d, const char* element) {
  Diagnostics diag;
  auto p = read_protein_pdb(text, diag);
  const auto& r = p->residues[0];
  const ProteinAtom *a = r.atom(at), *b = r.atom(n1), *c = r.atom(n2);
  double u[3] = {a->x - 0.5 * (b->x + c->x), a->y - 0.5 * (b->y + c->y), a->z - 0.5 * (b->z + c->z)};
  const double n = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  char line[128];
  const bool metal = std::string(element) == "ZN";
  std::snprintf(line, sizeof line, "HETATM   11 %-4s %3s A 201    %8.3f%8.3f%8.3f  1.00 10.00          %2s\n", metal ? "ZN" : " O1", metal ? " ZN" : "ACT",
                a->x + d * u[0] / n, a->y + d * u[1] / n, a->z + d * u[2] / n, element);
  return line;
}

}  // namespace

TEST(protein_mmcif_round_trip_and_format) {
  Diagnostics d;
  auto p = read_protein_pdb(slurp("1crn.pdb"), d);
  CHECK(p.has_value());
  CHECK(fits_pdb_format(*p));
  const std::string cif = write_protein_cif(*p);
  auto q = read_protein(cif, d);  // detected as mmCIF
  CHECK(q.has_value());
  CHECK_EQ(q->atom_count(), p->atom_count());
  CHECK_EQ(q->residues.size(), p->residues.size());
  CHECK_EQ(q->disulfides.size(), std::size_t{3});
  CHECK(q->chain_breaks.empty());
  for (std::size_t i = 0; i < p->residues.size(); ++i) {
    CHECK_EQ(q->residues[i].name, p->residues[i].name);
    CHECK_EQ(q->residues[i].number, p->residues[i].number);
    for (const auto& a : p->residues[i].atoms) {
      const ProteinAtom* b = q->residues[i].atom(a.name);
      CHECK(b != nullptr);
      if (b) CHECK(distance(a, *b) < 1e-3);
    }
  }
  // Quoted values, a multi-character chain name, and residues listed without coordinates.
  const std::string small =
      "data_T\n_struct.title 'two residues; one missing'\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
      "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n"
      "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n"
      "_atom_site.B_iso_or_equiv\n_atom_site.pdbx_PDB_model_num\n"
      "ATOM 1 N N . GLY AB1 1 ? 0.000 0.000 0.000 1.00 5.0 1\nATOM 2 C CA . GLY AB1 1 ? 1.458 0.000 0.000 1.00 5.0 1\n"
      "ATOM 3 C C . GLY AB1 1 ? 2.009 1.420 0.000 1.00 5.0 1\nATOM 4 O O . GLY AB1 1 ? 1.251 2.390 0.000 1.00 5.0 1\n"
      "HETATM 5 O \"O5'\" . LIG AB1 9 ? 9.000 9.000 9.000 1.00 5.0 1\n#\n"
      "loop_\n_pdbx_poly_seq_scheme.asym_id\n_pdbx_poly_seq_scheme.mon_id\n_pdbx_poly_seq_scheme.pdb_seq_num\n"
      "_pdbx_poly_seq_scheme.auth_seq_num\n_pdbx_poly_seq_scheme.pdb_strand_id\n_pdbx_poly_seq_scheme.pdb_ins_code\n"
      "A GLY 1 1 AB1 .\nA ALA 2 ? AB1 .\n#\n";
  auto s = read_protein(small, d);
  CHECK(s.has_value());
  CHECK_EQ(s->residues.size(), std::size_t{2});
  CHECK_EQ(s->residues[0].chain, std::string("AB1"));
  CHECK(s->residues[1].atom("O5'") != nullptr);
  CHECK_EQ(s->missing_residues.size(), std::size_t{1});
  CHECK_EQ(s->missing_residues[0].name, std::string("ALA"));
  CHECK(!fits_pdb_format(*s));
  CHECK(write_protein_cif(*s).find("\"O5'\"") != std::string::npos);
}

TEST(protein_rebuilds_missing_side_chains) {
  Diagnostics d;
  auto ref = read_protein_pdb(slurp("1crn.pdb"), d);
  CHECK(ref.has_value());
  // One side chain at a time: Tyr 29, Arg 10, Ile 33, Leu 18, Phe 13 (buried or well packed in crambin).
  for (int number : {29, 10, 33, 18, 13}) {
    Protein p = *ref;
    for (std::size_t i = 0; i < p.residues.size(); ++i) {
      auto& r = p.residues[i];
      if (r.number != number) continue;
      const std::size_t full = r.atoms.size();
      r.atoms.erase(std::remove_if(r.atoms.begin(), r.atoms.end(), [](const ProteinAtom& a) { return !backbone_name(a.name); }), r.atoms.end());
      const auto rep = repair(p, RepairOptions{}, d);
      CHECK_EQ(p.residues[i].atoms.size(), full);
      CHECK_EQ(rep.atoms_added, full - 4);
      CHECK_EQ(rep.completed_residues.size(), std::size_t{1});
      const auto& o = ref->residues[i];
      const auto& n = p.residues[i];
      CHECK(distance(*o.atom("CB"), *n.atom("CB")) < 0.25);  // rms 0.09 Å over crambin
      const char* g = o.atom("CG") ? "CG" : "CG1";
      CHECK(angle_difference(dihedral(*o.atom("N"), *o.atom("CA"), *o.atom("CB"), *o.atom(g)),
                             dihedral(*n.atom("N"), *n.atom("CA"), *n.atom("CB"), *n.atom(g))) < 40.0);
      for (const auto& a : n.atoms) CHECK(a.modelled == !backbone_name(a.name));
    }
  }
  // All side chains removed: most chi1 angles are recovered (26 of 32 = 81% when this test was written).
  Protein p = *ref;
  for (auto& r : p.residues)
    r.atoms.erase(std::remove_if(r.atoms.begin(), r.atoms.end(), [](const ProteinAtom& a) { return !backbone_name(a.name); }), r.atoms.end());
  repair(p, RepairOptions{}, d);
  CHECK_EQ(p.atom_count(), ref->atom_count());
  int n = 0, ok = 0;
  for (std::size_t i = 0; i < p.residues.size(); ++i) {
    const auto& o = ref->residues[i];
    const auto& r = p.residues[i];
    const char* g = nullptr;
    for (const char* c : {"CG", "CG1", "OG", "OG1", "SG"})
      if (o.atom(c)) g = c;
    if (!g || o.name == "PRO") continue;
    ++n;
    ok += angle_difference(dihedral(*o.atom("N"), *o.atom("CA"), *o.atom("CB"), *o.atom(g)),
                           dihedral(*r.atom("N"), *r.atom("CA"), *r.atom("CB"), *r.atom(g))) < 40.0;
  }
  CHECK_EQ(n, 32);
  CHECK(ok >= 22);
  // Nothing to do on a complete structure.
  Protein c = *ref;
  CHECK_EQ(repair(c, RepairOptions{}, d).atoms_added, std::size_t{0});
}

TEST(protein_builds_missing_loop) {
  Diagnostics d;
  auto ref = read_protein_pdb(slurp("1crn.pdb"), d);
  CHECK(ref.has_value());
  // Residues 20–23 (GTPE) removed and listed as missing.
  Protein p = *ref;
  std::size_t first = 0;
  for (std::size_t i = 0; i < p.residues.size(); ++i)
    if (p.residues[i].number == 20) first = i;
  for (std::size_t k = 0; k < 4; ++k) {
    const auto& r = p.residues[first + k];
    p.missing_residues.push_back({r.name, r.chain, r.number, r.icode});
  }
  p.residues.erase(p.residues.begin() + long(first), p.residues.begin() + long(first) + 4);
  for (auto& [a, b] : p.disulfides) a -= a > first ? 4 : 0, b -= b > first ? 4 : 0;  // indices after the removal
  const auto rep = repair(p, RepairOptions{}, d);
  CHECK_EQ(rep.built_residues.size(), std::size_t{4});
  CHECK(rep.not_built.empty());
  CHECK(p.missing_residues.empty());
  CHECK(p.chain_breaks.empty());
  CHECK_EQ(p.residues.size(), ref->residues.size());
  CHECK_EQ(p.atom_count(), ref->atom_count());
  CHECK_EQ(p.disulfides.size(), std::size_t{3});
  double sum2 = 0.0;
  int n = 0;
  for (std::size_t k = 0; k < 4; ++k) {
    const auto& r = p.residues[first + k];
    CHECK(r.predicted);
    CHECK_EQ(r.name, ref->residues[first + k].name);
    CHECK_EQ(r.number, 20 + int(k));
    for (const char* a : {"N", "CA", "C"}) sum2 += std::pow(distance(*r.atom(a), *ref->residues[first + k].atom(a)), 2), ++n;
  }
  // Peptide bonds along the rebuilt stretch, including both junctions.
  for (std::size_t i = first; i <= first + 4; ++i) CHECK_NEAR(distance(*p.residues[i - 1].atom("C"), *p.residues[i].atom("N")), 1.33, 0.06);
  CHECK(std::sqrt(sum2 / n) < 2.5);  // 4-residue loops: median 0.7 Å over the validation set
  // The disulfides still join the same residues after the insertion.
  for (const auto& [a, b] : p.disulfides) CHECK(distance(*p.residues[a].atom("SG"), *p.residues[b].atom("SG")) < 2.5);
  // Without loop building the residues stay missing.
  Protein q = *ref;
  q.missing_residues.push_back({"GLY", "A", 20, ' '});
  q.residues.erase(q.residues.begin() + long(first));
  RepairOptions no;
  no.loops = false;
  CHECK(repair(q, no, d).built_residues.empty());
  CHECK_EQ(q.missing_residues.size(), std::size_t{1});
}

TEST(protein_structure_pka) {
  Diagnostics d;
  auto p = read_protein_pdb(slurp("1crn.pdb"), d);
  CHECK(p.has_value());
  const auto sites = predict_pka(*p);
  // N-terminus, Arg 10, Arg 17, Glu 23, Tyr 29, Asp 43, Tyr 44, C-terminus; the six cysteines are in disulfides.
  CHECK_EQ(sites.size(), std::size_t{8});
  std::map<std::string, double> v;
  for (const auto& s : sites) {
    v[s.group + std::to_string(p->residues[s.residue].number)] = s.pka;
    CHECK_NEAR(s.pka, s.reference + s.desolvation + s.hydrogen_bonds + s.charges, 1e-12);
    CHECK(s.buried >= 0.0 && s.buried <= 1.0);
  }
  CHECK_NEAR(v["ARG10"], 12.5, 1e-12);  // no model for Arg, Tyr, Cys: reference values
  CHECK_NEAR(v["TYR29"], 10.0, 1e-12);
  CHECK_NEAR(v["GLU23"], 4.066, 0.01);  // regression values of the fitted model
  CHECK_NEAR(v["ASP43"], 3.700, 0.01);
  CHECK_NEAR(v["C-46"], 2.858, 0.01);
  CHECK_NEAR(v["N+1"], 7.806, 0.01);
  // A positive charge next to a carboxylate lowers its pKa; reference values on request.
  ProtonationOptions o;
  o.pH = 3.0;
  o.structure_pka = false;
  auto q = *p;
  const auto rep = protonate(q, o, d);
  CHECK(rep.pka.empty());
  for (const auto& r : q.residues)
    if (r.name == "ASP") CHECK_EQ(r.state, std::string("ASH"));  // 3.0 < 3.8
}

TEST(protein_histidine_tautomer_and_metal) {
  const std::string base = histidine_with("");
  struct Case {
    const char *at, *n1, *n2, *element;
    double d;
    const char* expect;
  };
  // An acceptor oxygen in line with a ring nitrogen puts the hydrogen on that nitrogen; a zinc ion on a
  // nitrogen leaves it bare and protonates the other one.
  for (const Case& c : {Case{"ND1", "CG", "CE1", " O", 2.8, "HID"}, Case{"NE2", "CD2", "CE1", " O", 2.8, "HIE"},
                        Case{"NE2", "CD2", "CE1", "ZN", 2.1, "HID"}, Case{"ND1", "CG", "CE1", "ZN", 2.1, "HIE"}}) {
    Diagnostics d;
    auto p = read_protein_pdb(histidine_with(probe(base, c.at, c.n1, c.n2, c.d, c.element)), d);
    CHECK(p.has_value());
    ProtonationOptions o;
    o.pH = 9.0;
    o.flips = false;
    const auto rep = protonate(*p, o, d);
    CHECK_EQ(p->residues[0].state, std::string(c.expect));
    CHECK_EQ(p->residues[0].atom("HD1") != nullptr, std::string(c.expect) == "HID");
    CHECK_EQ(p->residues[0].atom("HE2") != nullptr, std::string(c.expect) == "HIE");
    CHECK(p->residues[0].atom("HD2") && p->residues[0].atom("HE1"));
    CHECK(rep.undetermined.empty());
  }
  // With flips allowed, an acceptor facing CD2 turns the ring so that a nitrogen donates to it.
  {
    Diagnostics d;
    auto p = read_protein_pdb(histidine_with(probe(base, "CD2", "CG", "NE2", 2.8, " O")), d);
    const ProteinAtom old_cd2 = *p->residues[0].atom("CD2");
    ProtonationOptions o;
    o.pH = 9.0;
    const auto rep = protonate(*p, o, d);
    CHECK_EQ(rep.flipped.size(), std::size_t{1});
    CHECK(distance(*p->residues[0].atom("ND1"), old_cd2) < 1e-9);
    CHECK_EQ(p->residues[0].state, std::string("HID"));
  }
  // Alone, nothing decides: HIE by default, and reported; charged at low pH.
  {
    Diagnostics d;
    auto p = read_protein_pdb(base, d);
    ProtonationOptions o;
    o.pH = 9.0;
    auto rep = protonate(*p, o, d);
    CHECK_EQ(p->residues[0].state, std::string("HIE"));
    CHECK_EQ(rep.undetermined.size(), std::size_t{1});
    o.pH = 4.0;
    rep = protonate(*p, o, d);
    CHECK_EQ(p->residues[0].state, std::string("HIP"));
    CHECK(p->residues[0].atom("HD1") && p->residues[0].atom("HE2"));
  }
}

TEST(protein_amide_flips_restored) {
  // Trp-cage: Gln 5 donates a hydrogen bond from NE2. Turning the amide by 180° in the input, the network
  // turns it back. For Asn 1 the evidence (between 0.3 and 0.5 of a hydrogen bond) is below what a flip
  // needs: it stays as read and is reported as undetermined.
  const std::string text = slurp("1l2y_model1.pdb");
  Diagnostics d;
  auto ref = read_protein_pdb(text, d);
  auto p = read_protein_pdb(text, d);
  auto swap_atoms = [](ProteinResidue& r, const char* a, const char* b) {
    ProteinAtom *x = r.atom(a), *y = r.atom(b);
    std::swap(x->x, y->x), std::swap(x->y, y->y), std::swap(x->z, y->z);
  };
  for (auto& r : p->residues) {
    if (r.name == "ASN") swap_atoms(r, "OD1", "ND2");
    if (r.name == "GLN") swap_atoms(r, "OE1", "NE2");
  }
  const auto rep = protonate(*p, ProtonationOptions{}, d);
  CHECK_EQ(rep.flipped.size(), std::size_t{1});
  CHECK_EQ(rep.flipped[0], std::string("A 5 GLN"));
  CHECK_EQ(rep.undetermined.size(), std::size_t{1});
  CHECK(rep.undetermined[0].rfind("A 1 ASN", 0) == 0);
  for (std::size_t i = 0; i < p->residues.size(); ++i) {
    const auto& r = p->residues[i];
    if (r.name == "GLN") CHECK(distance(*r.atom("NE2"), *ref->residues[i].atom("NE2")) < 1e-9);
  }
  // As deposited, nothing is flipped; with flips disabled the turned input stays as it is.
  auto q = read_protein_pdb(text, d);
  CHECK(protonate(*q, ProtonationOptions{}, d).flipped.empty());
  ProtonationOptions no;
  no.flips = false;
  for (auto& r : q->residues)
    if (r.name == "GLN") swap_atoms(r, "OE1", "NE2");
  CHECK(protonate(*q, no, d).flipped.empty());
}

TEST(protein_json_options) {
  const std::string text = slurp("1crn.pdb");
  char* s = cl_protein_json_opts(text.c_str(), 7.0, 0);
  const std::string j(s);
  cl_free_string(s);
  CHECK(j.rfind("{\"ok\":true", 0) == 0);
  CHECK(j.find("\"format\":\"pdb\"") != std::string::npos);
  CHECK(j.find("\"group\":\"GLU\"") != std::string::npos);
  CHECK(j.find("\"titration_structure\":[[0,") != std::string::npos);
  CHECK(j.find("\"repair\":{\"atoms_added\":0") != std::string::npos);
  CHECK(j.find("\"undetermined\":[\"A 12 ASN: amide orientation\"") != std::string::npos);
  // mmCIF input through the same entry point gives the same result.
  Diagnostics d;
  auto p = read_protein_pdb(text, d);
  char* c = cl_protein_json_opts(write_protein_cif(*p).c_str(), 7.0, 0);
  const std::string jc(c);
  cl_free_string(c);
  CHECK(jc.find("\"hydrogens_added\":315") != std::string::npos);
  CHECK(jc.find("\"disulfides\":[[\"A\",3,\"A\",40]") != std::string::npos);
  char* bad = cl_protein_json_opts(text.c_str(), 15.0, 0);
  CHECK(std::string(bad).find("\"ok\":false") != std::string::npos);
  cl_free_string(bad);
}
