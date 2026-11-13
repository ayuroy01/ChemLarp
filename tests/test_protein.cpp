// Proteins: PDB reading, disulfides, protonation, hydrogen placement, sequence properties.
//
// data/1crn.pdb is crambin (PDB 1CRN, X-ray, no hydrogens); data/1l2y_model1.pdb is the first model of the
// Trp-cage NMR ensemble (PDB 1L2Y) with its deposited hydrogens. Reference pI and net charges are from
// Biopython's IsoelectricPoint (same Bjellqvist pKa set); see "Proteins: validation" in docs/ARCHITECTURE.md.
#include <cctype>
#include <cmath>
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
  CHECK(n > 140);
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
