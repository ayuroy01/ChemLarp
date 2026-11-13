#include "chemlarp/protein.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <numbers>
#include <set>
#include <sstream>

#include "chemlarp/elements.hpp"
#include "chemlarp/json.hpp"

namespace chemlarp {

const ProteinAtom* ProteinResidue::atom(std::string_view n) const {
  for (const auto& a : atoms)
    if (a.name == n) return &a;
  return nullptr;
}
ProteinAtom* ProteinResidue::atom(std::string_view n) {
  for (auto& a : atoms)
    if (a.name == n) return &a;
  return nullptr;
}
std::size_t Protein::atom_count() const {
  std::size_t n = 0;
  for (const auto& r : residues) n += r.atoms.size();
  return n;
}

namespace {

// ---------------------------------------------------------------------------------------------------------
// Templates: heavy-atom bonds inside each standard residue and the hydrogens on each parent atom.

struct HDef {
  const char* parent;
  std::vector<const char*> names;
  char geom;          // 'T' tetrahedral, 'P' planar sp2, 'R' rotatable OH/SH, 'A' aromatic OH (in plane)
  const char* prev;   // T with two heavy neighbours: the neighbour toward the backbone (prochirality);
                      // P with one heavy neighbour: the reference atom whose side is "cis"
  bool first_cis;     // P with two H: whether the first name is cis to `prev`
};

struct ResidueDef {
  char code;
  std::vector<std::pair<const char*, const char*>> bonds;  // side chain (CA-CB onward)
  std::vector<HDef> h;
  std::array<int, 5> formula;  // C H N O S of the residue in a chain (free amino acid minus water)
};

const std::map<std::string, ResidueDef>& residue_defs() {
  static const std::map<std::string, ResidueDef> d = {
      {"ALA", {'A', {{"CA", "CB"}}, {{"CB", {"HB1", "HB2", "HB3"}, 'T', nullptr, false}}, {3, 5, 1, 1, 0}}},
      {"ARG", {'R', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "NE"}, {"NE", "CZ"}, {"CZ", "NH1"}, {"CZ", "NH2"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"CG", {"HG2", "HG3"}, 'T', "CB", false},
                {"CD", {"HD2", "HD3"}, 'T', "CG", false}, {"NE", {"HE"}, 'P', nullptr, false},
                {"NH1", {"HH11", "HH12"}, 'P', "NE", true}, {"NH2", {"HH21", "HH22"}, 'P', "NE", true}},
               {6, 12, 4, 1, 0}}},
      {"ASN", {'N', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "OD1"}, {"CG", "ND2"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"ND2", {"HD21", "HD22"}, 'P', "OD1", false}}, {4, 6, 2, 2, 0}}},
      {"ASP", {'D', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "OD1"}, {"CG", "OD2"}}, {{"CB", {"HB2", "HB3"}, 'T', "CA", false}},
               {4, 5, 1, 3, 0}}},
      {"CYS", {'C', {{"CA", "CB"}, {"CB", "SG"}}, {{"CB", {"HB2", "HB3"}, 'T', "CA", false}}, {3, 5, 1, 1, 1}}},
      {"GLN", {'Q', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "OE1"}, {"CD", "NE2"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"CG", {"HG2", "HG3"}, 'T', "CB", false},
                {"NE2", {"HE21", "HE22"}, 'P', "OE1", false}},
               {5, 8, 2, 2, 0}}},
      {"GLU", {'E', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "OE1"}, {"CD", "OE2"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"CG", {"HG2", "HG3"}, 'T', "CB", false}}, {5, 7, 1, 3, 0}}},
      {"GLY", {'G', {}, {}, {2, 3, 1, 1, 0}}},
      {"HIS", {'H', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "ND1"}, {"CG", "CD2"}, {"ND1", "CE1"}, {"CE1", "NE2"}, {"NE2", "CD2"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"CD2", {"HD2"}, 'P', nullptr, false},
                {"CE1", {"HE1"}, 'P', nullptr, false}},
               {6, 7, 3, 1, 0}}},
      {"ILE", {'I', {{"CA", "CB"}, {"CB", "CG1"}, {"CB", "CG2"}, {"CG1", "CD1"}},
               {{"CB", {"HB"}, 'T', nullptr, false}, {"CG1", {"HG12", "HG13"}, 'T', "CB", false},
                {"CG2", {"HG21", "HG22", "HG23"}, 'T', nullptr, false}, {"CD1", {"HD11", "HD12", "HD13"}, 'T', nullptr, false}},
               {6, 11, 1, 1, 0}}},
      {"LEU", {'L', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CG", "CD2"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"CG", {"HG"}, 'T', nullptr, false},
                {"CD1", {"HD11", "HD12", "HD13"}, 'T', nullptr, false}, {"CD2", {"HD21", "HD22", "HD23"}, 'T', nullptr, false}},
               {6, 11, 1, 1, 0}}},
      {"LYS", {'K', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "CE"}, {"CE", "NZ"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"CG", {"HG2", "HG3"}, 'T', "CB", false},
                {"CD", {"HD2", "HD3"}, 'T', "CG", false}, {"CE", {"HE2", "HE3"}, 'T', "CD", false}},
               {6, 12, 2, 1, 0}}},
      {"MET", {'M', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "SD"}, {"SD", "CE"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"CG", {"HG2", "HG3"}, 'T', "CB", false},
                {"CE", {"HE1", "HE2", "HE3"}, 'T', nullptr, false}},
               {5, 9, 1, 1, 1}}},
      {"PHE", {'F', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CG", "CD2"}, {"CD1", "CE1"}, {"CD2", "CE2"}, {"CE1", "CZ"}, {"CE2", "CZ"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"CD1", {"HD1"}, 'P', nullptr, false}, {"CD2", {"HD2"}, 'P', nullptr, false},
                {"CE1", {"HE1"}, 'P', nullptr, false}, {"CE2", {"HE2"}, 'P', nullptr, false}, {"CZ", {"HZ"}, 'P', nullptr, false}},
               {9, 9, 1, 1, 0}}},
      {"PRO", {'P', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "N"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"CG", {"HG2", "HG3"}, 'T', "CB", false},
                {"CD", {"HD2", "HD3"}, 'T', "CG", false}},
               {5, 7, 1, 1, 0}}},
      {"SER", {'S', {{"CA", "CB"}, {"CB", "OG"}}, {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"OG", {"HG"}, 'R', nullptr, false}},
               {3, 5, 1, 2, 0}}},
      {"THR", {'T', {{"CA", "CB"}, {"CB", "OG1"}, {"CB", "CG2"}},
               {{"CB", {"HB"}, 'T', nullptr, false}, {"OG1", {"HG1"}, 'R', nullptr, false},
                {"CG2", {"HG21", "HG22", "HG23"}, 'T', nullptr, false}},
               {4, 7, 1, 2, 0}}},
      {"TRP", {'W', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CG", "CD2"}, {"CD1", "NE1"}, {"NE1", "CE2"}, {"CE2", "CD2"},
                     {"CE2", "CZ2"}, {"CD2", "CE3"}, {"CE3", "CZ3"}, {"CZ3", "CH2"}, {"CH2", "CZ2"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"CD1", {"HD1"}, 'P', nullptr, false}, {"NE1", {"HE1"}, 'P', nullptr, false},
                {"CE3", {"HE3"}, 'P', nullptr, false}, {"CZ2", {"HZ2"}, 'P', nullptr, false}, {"CZ3", {"HZ3"}, 'P', nullptr, false},
                {"CH2", {"HH2"}, 'P', nullptr, false}},
               {11, 10, 2, 1, 0}}},
      {"TYR", {'Y', {{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CG", "CD2"}, {"CD1", "CE1"}, {"CD2", "CE2"}, {"CE1", "CZ"}, {"CE2", "CZ"}, {"CZ", "OH"}},
               {{"CB", {"HB2", "HB3"}, 'T', "CA", false}, {"CD1", {"HD1"}, 'P', nullptr, false}, {"CD2", {"HD2"}, 'P', nullptr, false},
                {"CE1", {"HE1"}, 'P', nullptr, false}, {"CE2", {"HE2"}, 'P', nullptr, false}},
               {9, 9, 1, 2, 0}}},
      {"VAL", {'V', {{"CA", "CB"}, {"CB", "CG1"}, {"CB", "CG2"}},
               {{"CB", {"HB"}, 'T', nullptr, false}, {"CG1", {"HG11", "HG12", "HG13"}, 'T', nullptr, false},
                {"CG2", {"HG21", "HG22", "HG23"}, 'T', nullptr, false}},
               {5, 9, 1, 1, 0}}},
  };
  return d;
}

// Common alternative residue names (protonation-state names of force fields, selenomethionine excluded).
std::string canonical_name(const std::string& n) {
  static const std::map<std::string, std::string> alias = {
      {"HID", "HIS"}, {"HIE", "HIS"}, {"HIP", "HIS"}, {"HSD", "HIS"}, {"HSE", "HIS"}, {"HSP", "HIS"}, {"HSE", "HIS"},
      {"ASH", "ASP"}, {"GLH", "GLU"}, {"LYN", "LYS"}, {"CYX", "CYS"}, {"CYM", "CYS"}, {"TYM", "TYR"}, {"ARN", "ARG"}};
  auto it = alias.find(n);
  return it != alias.end() ? it->second : n;
}

using Vec = std::array<double, 3>;
Vec pos(const ProteinAtom& a) { return {a.x, a.y, a.z}; }
Vec sub(Vec a, Vec b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec add(Vec a, Vec b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec mul(Vec a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(Vec a, Vec b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec cross(Vec a, Vec b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
double norm(Vec a) { return std::sqrt(dot(a, a)); }
Vec unit(Vec a) {
  const double n = norm(a);
  return n > 0 ? mul(a, 1.0 / n) : Vec{1, 0, 0};
}
double dist(const ProteinAtom& a, const ProteinAtom& b) { return norm(sub(pos(a), pos(b))); }

// Position of atom D given A-B-C, bond C-D, angle B-C-D and dihedral A-B-C-D (degrees).
Vec place(Vec A, Vec B, Vec C, double bond, double angle_deg, double dihedral_deg) {
  const double ang = angle_deg * std::numbers::pi / 180.0, dih = dihedral_deg * std::numbers::pi / 180.0;
  const Vec bc = unit(sub(C, B));
  const Vec n = unit(cross(sub(B, A), bc));
  const Vec m = cross(n, bc);
  const Vec d2 = {-bond * std::cos(ang), bond * std::sin(ang) * std::cos(dih), bond * std::sin(ang) * std::sin(dih)};
  return add(C, add(add(mul(bc, d2[0]), mul(m, d2[1])), mul(n, d2[2])));
}

std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(' ');
  if (a == std::string::npos) return "";
  return s.substr(a, s.find_last_not_of(' ') - a + 1);
}

std::string residue_label(const ProteinResidue& r) {
  return (r.chain.empty() ? std::string("_") : r.chain) + " " + std::to_string(r.number) +
         (r.icode != ' ' ? std::string(1, r.icode) : "") + " " + r.name;
}

// Model pKa values (PROPKA 3).
struct Pka {
  double asp = 3.80, glu = 4.50, his = 6.50, cys = 9.00, tyr = 10.00, lys = 10.50, arg = 12.50, nterm = 8.00, cterm = 3.20;
};

}  // namespace

std::optional<Protein> read_protein_pdb(std::string_view text, Diagnostics& diag) {
  Protein p;
  std::istringstream in{std::string(text)};
  std::map<std::string, std::size_t> index;  // residue key -> index
  struct Pending {
    std::size_t res;
    ProteinAtom atom;
    char alt;
  };
  std::map<std::pair<std::size_t, std::string>, std::pair<ProteinAtom, double>> chosen;  // (res, name) -> atom, occupancy
  std::vector<std::pair<std::size_t, std::string>> order;
  std::vector<std::array<std::string, 4>> ssbond;
  bool in_model = false, model_done = false;
  for (std::string line; std::getline(in, line);) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const std::string rec = line.substr(0, 6);
    if (rec == "TITLE " && p.title.size() < 300) p.title += (p.title.empty() ? "" : " ") + trim(line.size() > 10 ? line.substr(10) : "");
    if (rec == "MODEL ") {
      if (in_model || model_done) model_done = true;
      in_model = true;
    }
    if (rec == "ENDMDL") {
      model_done = true;
      in_model = false;
    }
    if (rec == "SSBOND" && line.size() >= 35)
      ssbond.push_back({line.substr(15, 1), trim(line.substr(17, 5)), line.substr(29, 1), trim(line.substr(31, 5))});
    if ((rec != "ATOM  " && rec != "HETATM") || model_done) continue;
    if (line.size() < 54) {
      diag.warn(DiagCode::MalformedRecord, "short ATOM record skipped: " + line);
      continue;
    }
    ProteinAtom a;
    a.name = trim(line.substr(12, 4));
    const char alt = line[16];
    const std::string resname = trim(line.substr(17, 3));
    const std::string chain = trim(line.substr(21, 1));
    const int num = std::atoi(line.substr(22, 4).c_str());
    const char icode = line.size() > 26 ? line[26] : ' ';
    a.x = std::atof(line.substr(30, 8).c_str());
    a.y = std::atof(line.substr(38, 8).c_str());
    a.z = std::atof(line.substr(46, 8).c_str());
    a.occupancy = line.size() >= 60 ? std::atof(line.substr(54, 6).c_str()) : 1.0;
    a.bfactor = line.size() >= 66 ? std::atof(line.substr(60, 6).c_str()) : 0.0;
    std::string el = line.size() >= 78 ? trim(line.substr(76, 2)) : "";
    if (el.empty()) {  // from the atom name: first letter(s) that are letters
      std::string nm = trim(line.substr(12, 2));
      el = nm.empty() ? "" : std::string(1, nm[0] == ' ' || std::isdigit(static_cast<unsigned char>(nm[0])) ? nm[1] : nm[0]);
    }
    if (el.size() == 2) el[1] = char(std::tolower(static_cast<unsigned char>(el[1])));
    a.element = el;
    if (line.size() >= 80) {
      const std::string c = trim(line.substr(78, 2));
      if (c.size() == 2 && std::isdigit(static_cast<unsigned char>(c[0]))) a.charge = (c[1] == '-' ? -1 : 1) * (c[0] - '0');
    }
    const std::string key = chain + "|" + std::to_string(num) + "|" + std::string(1, icode) + "|" + resname;
    auto it = index.find(key);
    if (it == index.end()) {
      ProteinResidue r;
      r.name = resname;
      r.chain = chain;
      r.number = num;
      r.icode = icode;
      r.hetero = rec == "HETATM" && !residue_defs().count(canonical_name(resname));
      p.residues.push_back(r);
      it = index.emplace(key, p.residues.size() - 1).first;
    }
    const auto k = std::make_pair(it->second, a.name);
    auto c = chosen.find(k);
    if (c == chosen.end()) {
      chosen.emplace(k, std::make_pair(a, a.occupancy));
      order.push_back(k);
    } else if (alt != ' ' && a.occupancy > c->second.second) {
      c->second = {a, a.occupancy};
    }
  }
  for (const auto& k : order) p.residues[k.first].atoms.push_back(chosen[k].first);
  if (p.residues.empty()) {
    diag.error(DiagCode::UnknownFileFormat, "no ATOM or HETATM records");
    return std::nullopt;
  }
  // Disulfides: SSBOND records and S–S distance.
  std::set<std::pair<std::size_t, std::size_t>> ss;
  std::vector<std::size_t> cys;
  for (std::size_t i = 0; i < p.residues.size(); ++i)
    if (canonical_name(p.residues[i].name) == "CYS" && p.residues[i].atom("SG")) cys.push_back(i);
  for (std::size_t a = 0; a < cys.size(); ++a)
    for (std::size_t b = a + 1; b < cys.size(); ++b)
      if (dist(*p.residues[cys[a]].atom("SG"), *p.residues[cys[b]].atom("SG")) < 2.5) ss.insert({cys[a], cys[b]});
  for (const auto& s : ssbond) {
    std::size_t ia = SIZE_MAX, ib = SIZE_MAX;
    for (std::size_t i : cys) {
      const auto& r = p.residues[i];
      if ((r.chain.empty() ? " " : r.chain) == (s[0] == "" ? " " : s[0]) && std::to_string(r.number) == s[1]) ia = i;
      if ((r.chain.empty() ? " " : r.chain) == (s[2] == "" ? " " : s[2]) && std::to_string(r.number) == s[3]) ib = i;
    }
    if (ia != SIZE_MAX && ib != SIZE_MAX) ss.insert({std::min(ia, ib), std::max(ia, ib)});
  }
  p.disulfides.assign(ss.begin(), ss.end());
  // Termini and chain breaks among amino-acid residues of each chain.
  for (std::size_t i = 0; i < p.residues.size(); ++i) {
    auto& r = p.residues[i];
    if (!residue_defs().count(canonical_name(r.name))) continue;
    std::size_t prev = SIZE_MAX, next = SIZE_MAX;
    for (std::size_t j = i; j-- > 0;)
      if (p.residues[j].chain == r.chain && residue_defs().count(canonical_name(p.residues[j].name))) { prev = j; break; }
    for (std::size_t j = i + 1; j < p.residues.size(); ++j)
      if (p.residues[j].chain == r.chain && residue_defs().count(canonical_name(p.residues[j].name))) { next = j; break; }
    const ProteinAtom* N = r.atom("N");
    const ProteinAtom* C = r.atom("C");
    bool linked_prev = false, linked_next = false;
    if (prev != SIZE_MAX && N && p.residues[prev].atom("C")) linked_prev = dist(*p.residues[prev].atom("C"), *N) < 2.0;
    if (next != SIZE_MAX && C && p.residues[next].atom("N")) linked_next = dist(*C, *p.residues[next].atom("N")) < 2.0;
    r.n_terminal = !linked_prev;
    r.c_terminal = !linked_next;
    if (prev != SIZE_MAX && !linked_prev) p.chain_breaks.push_back(i);
  }
  return p;
}

std::vector<std::pair<std::string, std::string>> chain_sequences(const Protein& p) {
  std::vector<std::pair<std::string, std::string>> out;
  for (const auto& r : p.residues) {
    auto it = residue_defs().find(canonical_name(r.name));
    if (it == residue_defs().end()) continue;
    if (out.empty() || out.back().first != r.chain) out.push_back({r.chain, ""});
    out.back().second += it->second.code;
  }
  return out;
}

ProtonationReport protonate(Protein& p, const ProtonationOptions& opt, Diagnostics& diag) {
  ProtonationReport rep;
  const Pka K;
  const double pH = opt.pH;
  // Remove hydrogens.
  for (auto& r : p.residues) {
    if (r.hetero) continue;
    const std::size_t before = r.atoms.size();
    r.atoms.erase(std::remove_if(r.atoms.begin(), r.atoms.end(), [](const ProteinAtom& a) { return a.element == "H" || a.element == "D"; }),
                  r.atoms.end());
    rep.hydrogens_removed += before - r.atoms.size();
  }
  std::set<std::size_t> in_ss;
  for (const auto& [a, b] : p.disulfides) in_ss.insert(a), in_ss.insert(b);
  // Heavy atoms of everything, for hydrogen-bond acceptors and metals.
  struct Site {
    Vec x;
    std::string el;
    std::size_t res;
    std::string name;
  };
  std::vector<Site> heavy;
  for (std::size_t i = 0; i < p.residues.size(); ++i)
    for (const auto& a : p.residues[i].atoms) heavy.push_back({pos(a), a.element, i, a.name});
  auto is_metal = [](const std::string& e) {
    return e == "Zn" || e == "Fe" || e == "Cu" || e == "Mg" || e == "Ca" || e == "Mn" || e == "Co" || e == "Ni" || e == "Na" || e == "K";
  };
  // Acceptor score of a polar hydrogen position: distance to the nearest O/N of another residue (or of
  // this residue's backbone carbonyl), ignoring the parent's own neighbours.
  auto acceptor_distance = [&](Vec h, std::size_t res, const std::string& parent) {
    double best = 1e9;
    for (const auto& s : heavy) {
      if (s.el != "O" && s.el != "N") continue;
      if (s.res == res && (s.name == parent || s.name != "O")) continue;
      best = std::min(best, norm(sub(s.x, h)));
    }
    return best;
  };
  const auto& defs = residue_defs();
  for (std::size_t ri = 0; ri < p.residues.size(); ++ri) {
    auto& r = p.residues[ri];
    const std::string base = canonical_name(r.name);
    auto dit = defs.find(base);
    if (dit == defs.end()) {
      if (!r.hetero || (r.name != "HOH" && r.name != "WAT"))
        if (std::find(rep.unknown_residues.begin(), rep.unknown_residues.end(), r.name) == rep.unknown_residues.end() && !r.hetero)
          rep.unknown_residues.push_back(r.name);
      continue;
    }
    const ResidueDef& def = dit->second;
    // Missing heavy atoms.
    std::vector<std::string> need = {"N", "CA", "C", "O"};
    for (const auto& [a, b] : def.bonds) need.push_back(a), need.push_back(b);
    std::sort(need.begin(), need.end());
    need.erase(std::unique(need.begin(), need.end()), need.end());
    for (const auto& n : need)
      if (!r.atom(n)) rep.missing_atoms.push_back(residue_label(r) + " " + n);
    // Terminal OXT.
    if (r.c_terminal && opt.add_oxt && !r.atom("OXT") && r.atom("C") && r.atom("O") && r.atom("CA")) {
      const Vec c = pos(*r.atom("C"));
      const Vec dir = unit(add(unit(sub(pos(*r.atom("CA")), c)), unit(sub(pos(*r.atom("O")), c))));
      ProteinAtom o;
      o.name = "OXT";
      o.element = "O";
      const Vec x = sub(c, mul(dir, 1.25));
      o.x = x[0], o.y = x[1], o.z = x[2];
      o.added = true;
      r.atoms.push_back(o);
    }
    // Protonation state.
    int charge = 0;
    std::string state = base;
    bool his_d = false, his_e = false;
    if (base == "ASP") {
      state = pH < K.asp ? "ASH" : "ASP";
      charge = pH < K.asp ? 0 : -1;
    } else if (base == "GLU") {
      state = pH < K.glu ? "GLH" : "GLU";
      charge = pH < K.glu ? 0 : -1;
    } else if (base == "LYS") {
      state = pH < K.lys ? "LYS" : "LYN";
      charge = pH < K.lys ? 1 : 0;
    } else if (base == "ARG") {
      state = pH < K.arg ? "ARG" : "ARN";
      charge = pH < K.arg ? 1 : 0;
    } else if (base == "TYR") {
      state = pH < K.tyr ? "TYR" : "TYM";
      charge = pH < K.tyr ? 0 : -1;
    } else if (base == "CYS") {
      if (in_ss.count(ri)) state = "CYX";
      else {
        state = pH < K.cys ? "CYS" : "CYM";
        charge = pH < K.cys ? 0 : -1;
      }
    } else if (base == "HIS") {
      const ProteinAtom* nd = r.atom("ND1");
      const ProteinAtom* ne = r.atom("NE2");
      bool nd_metal = false, ne_metal = false, nd_hb = false, ne_hb = false;
      for (const auto& s : heavy) {
        if (s.res == ri) continue;
        for (int k = 0; k < 2; ++k) {
          const ProteinAtom* n = k ? ne : nd;
          if (!n) continue;
          const double d = norm(sub(s.x, pos(*n)));
          if (is_metal(s.el) && d < 2.6) (k ? ne_metal : nd_metal) = true;
          if ((s.el == "O" || s.el == "N") && d < 3.2) (k ? ne_hb : nd_hb) = true;
        }
      }
      if (pH < K.his && !nd_metal && !ne_metal) {
        state = "HIP";
        his_d = his_e = true;
        charge = 1;
      } else if (ne_metal || (nd_hb && !ne_hb && !nd_metal)) {
        state = "HID";
        his_d = true;
      } else {
        state = "HIE";
        his_e = true;
      }
    }
    // Termini.
    const bool nh3 = r.n_terminal && pH < K.nterm;
    const bool cooh = r.c_terminal && pH < K.cterm;
    if (r.n_terminal) charge += nh3 ? 1 : 0;
    if (r.c_terminal) charge += cooh ? 0 : -1;
    r.state = state;
    r.charge = charge;
    rep.net_charge += charge;
    // Neighbours of an atom (inside the residue, plus the peptide links).
    auto neighbours = [&](const std::string& name) {
      std::vector<Vec> out;
      std::vector<std::string> names;
      auto addn = [&](const ProteinAtom* a, const std::string& n) {
        if (a) out.push_back(pos(*a)), names.push_back(n);
      };
      std::vector<std::pair<std::string, std::string>> bonds = {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"C", "OXT"}};
      bonds.insert(bonds.end(), def.bonds.begin(), def.bonds.end());
      for (const auto& [a, b] : bonds) {
        if (name == a) addn(r.atom(b), b);
        if (name == b) addn(r.atom(a), a);
      }
      if (name == "N" && !r.n_terminal)
        for (std::size_t j = ri; j-- > 0;)
          if (p.residues[j].chain == r.chain && defs.count(canonical_name(p.residues[j].name))) {
            addn(p.residues[j].atom("C"), "C-");
            break;
          }
      return std::make_pair(out, names);
    };
    auto push = [&](const std::string& name, Vec x) {
      ProteinAtom h;
      h.name = name;
      h.element = "H";
      h.x = x[0], h.y = x[1], h.z = x[2];
      h.added = true;
      r.atoms.push_back(h);
      ++rep.hydrogens_added;
    };
    auto bond_len = [](const std::string& parent_el, bool aromatic) {
      if (parent_el == "N") return 1.01;
      if (parent_el == "O") return 0.96;
      if (parent_el == "S") return 1.34;
      return aromatic ? 1.08 : 1.09;
    };
    // Hydrogen definitions for this residue and state.
    std::vector<HDef> hs;
    if (base == "GLY") hs.push_back({"CA", {"HA2", "HA3"}, 'T', "N", false});
    else hs.push_back({"CA", {"HA"}, 'T', nullptr, false});
    if (r.n_terminal) {
      if (base == "PRO") hs.push_back({"N", nh3 ? std::vector<const char*>{"H2", "H3"} : std::vector<const char*>{"H2"}, 'T', nullptr, false});
      else hs.push_back({"N", nh3 ? std::vector<const char*>{"H1", "H2", "H3"} : std::vector<const char*>{"H1", "H2"}, 'T', nullptr, false});
    } else if (base != "PRO") {
      hs.push_back({"N", {"H"}, 'P', nullptr, false});
    }
    for (const auto& h : def.h) {
      HDef x = h;
      if (base == "LYS" && state == "LYN" && std::string(h.parent) == "NZ") continue;
      hs.push_back(x);
    }
    if (base == "LYS") hs.push_back({"NZ", state == "LYS" ? std::vector<const char*>{"HZ1", "HZ2", "HZ3"} : std::vector<const char*>{"HZ1", "HZ2"}, 'T', nullptr, false});
    if (base == "ARG" && state == "ARN") {
      for (auto& x : hs)
        if (std::string(x.parent) == "NH2") x.names = {"HH21"};
    }
    if (base == "SER") {}
    if (base == "CYS" && state == "CYS") hs.push_back({"SG", {"HG"}, 'R', nullptr, false});
    if (base == "TYR" && state == "TYR") hs.push_back({"OH", {"HH"}, 'A', nullptr, false});
    if (base == "HIS") {
      if (his_d) hs.push_back({"ND1", {"HD1"}, 'P', nullptr, false});
      if (his_e) hs.push_back({"NE2", {"HE2"}, 'P', nullptr, false});
    }
    if (state == "ASH") hs.push_back({"OD2", {"HD2"}, 'C', "OD1", false});
    if (state == "GLH") hs.push_back({"OE2", {"HE2"}, 'C', "OE1", false});
    if (cooh && r.atom("OXT")) hs.push_back({"OXT", {"HXT"}, 'C', "O", false});
    // Place them.
    for (const auto& h : hs) {
      const ProteinAtom* P = r.atom(h.parent);
      if (!P) continue;
      const Vec c = pos(*P);
      auto [nb, nbn] = neighbours(h.parent);
      if (nb.empty()) continue;
      const bool aromatic = h.geom == 'P' && P->element == "C";
      const double d = bond_len(P->element, aromatic);
      const std::size_t n = h.names.size();
      if (h.geom == 'T') {
        if (nb.size() >= 3 && n == 1) {
          Vec s{0, 0, 0};
          for (const auto& v : nb) s = add(s, unit(sub(v, c)));
          push(h.names[0], add(c, mul(unit(mul(s, -1.0)), d)));
        } else if (nb.size() == 2 && n == 2) {
          // Prochirality: H2 has a negative sign of ((prev − C) × (next − C)) · (H − C).
          Vec prev = nb[0], next = nb[1];
          if (h.prev && nbn[1] == h.prev) std::swap(prev, next);
          const Vec u1 = unit(sub(prev, c)), u2 = unit(sub(next, c));
          const Vec b = unit(mul(add(u1, u2), -1.0));
          const Vec nn = unit(cross(u1, u2));
          const double a = 54.735 * std::numbers::pi / 180.0;
          const Vec ha = add(c, mul(add(mul(b, std::cos(a)), mul(nn, std::sin(a))), d));
          const Vec hb = add(c, mul(add(mul(b, std::cos(a)), mul(nn, -std::sin(a))), d));
          const bool a_is_2 = dot(cross(sub(prev, c), sub(next, c)), sub(ha, c)) < 0;
          push(h.names[0], a_is_2 ? ha : hb);
          push(h.names[1], a_is_2 ? hb : ha);
        } else if (nb.size() == 2 && n == 1) {
          Vec s = add(unit(sub(nb[0], c)), unit(sub(nb[1], c)));
          push(h.names[0], add(c, mul(unit(mul(s, -1.0)), d)));
        } else if (nb.size() == 1) {
          // Methyl / NH3+: staggered against a reference atom bonded to the neighbour.
          const std::string nbname = nbn[0];
          auto [nb2, nb2n] = neighbours(nbname == "C-" ? "C" : nbname);
          Vec ref = c;
          for (std::size_t k = 0; k < nb2.size(); ++k)
            if (nb2n[k] != h.parent) {
              ref = nb2[k];
              break;
            }
          if (nbname == "CA" && std::string(h.parent) == "N") ref = r.atom("C") ? pos(*r.atom("C")) : ref;
          const double dihs[3] = {180.0, 60.0, -60.0};
          for (std::size_t k = 0; k < n && k < 3; ++k) push(h.names[k], place(ref, nb[0], c, d, 109.47, dihs[k]));
        }
      } else if (h.geom == 'P') {
        if (nb.size() >= 2 && n == 1) {
          Vec s = add(unit(sub(nb[0], c)), unit(sub(nb[1], c)));
          push(h.names[0], add(c, mul(unit(mul(s, -1.0)), d)));
        } else if (nb.size() == 1) {
          // sp2 NH2 (or NH): in the plane of the neighbour X and the reference atom R bonded to X.
          const Vec X = nb[0];
          auto [nb2, nb2n] = neighbours(nbn[0]);
          Vec R = c;
          for (std::size_t k = 0; k < nb2.size(); ++k)
            if ((h.prev && nb2n[k] == h.prev) || (!h.prev && nb2n[k] != h.parent)) {
              R = nb2[k];
              break;
            }
          // cis to R: dihedral R-X-P-H = 0; trans: 180.
          const Vec hc = place(R, X, c, d, 120.0, 0.0), ht = place(R, X, c, d, 120.0, 180.0);
          if (n == 2) {
            push(h.names[0], h.first_cis ? hc : ht);
            push(h.names[1], h.first_cis ? ht : hc);
          } else {
            push(h.names[0], ht);
          }
        }
      } else if (h.geom == 'R' || h.geom == 'A' || h.geom == 'C') {
        // One H on O or S with a single heavy neighbour X; dihedral about X–P chosen below.
        const Vec X = nb[0];
        auto [nb2, nb2n] = neighbours(nbn[0]);
        Vec R = c;
        for (std::size_t k = 0; k < nb2.size(); ++k)
          if ((h.prev && nb2n[k] == h.prev) || (!h.prev && nb2n[k] != h.parent)) {
            R = nb2[k];
            break;
          }
        const double angle = P->element == "S" ? 96.0 : (h.geom == 'C' ? 108.0 : 109.47);
        std::vector<double> cand;
        if (h.geom == 'C') cand = {0.0};               // carboxylic acid: syn
        else if (h.geom == 'A') cand = {0.0, 180.0};    // phenol: in the ring plane
        else cand = {180.0, 60.0, -60.0};               // alcohol, thiol: staggered
        double best_d = 1e9;
        Vec best = place(R, X, c, d, angle, cand[0]);
        for (double t : cand) {
          const Vec x = place(R, X, c, d, angle, t);
          const double a = acceptor_distance(x, ri, h.parent);
          if (a < best_d - 1e-9 && a < 2.6) best_d = a, best = x;
        }
        push(h.names[0], best);
      }
    }
  }
  if (!rep.missing_atoms.empty())
    diag.warn(DiagCode::MalformedRecord, std::to_string(rep.missing_atoms.size()) +
                                             " heavy atoms of standard residues are missing (e.g. " + rep.missing_atoms[0] +
                                             "); their hydrogens are not placed");
  if (!rep.unknown_residues.empty()) {
    std::string l;
    for (const auto& u : rep.unknown_residues) l += (l.empty() ? "" : ", ") + u;
    diag.warn(DiagCode::FileContentIgnored, "no template for residue(s) " + l + ": left without hydrogens");
  }
  return rep;
}

std::string write_protein_pdb(const Protein& p) {
  std::string out;
  char b[128];
  if (!p.title.empty()) out += "REMARK   1 " + p.title.substr(0, 68) + "\n";
  out += "REMARK   2 WRITTEN BY CHEMLARP\n";
  int serial = 0;
  for (std::size_t k = 0; k < p.disulfides.size(); ++k) {
    const auto& a = p.residues[p.disulfides[k].first];
    const auto& c = p.residues[p.disulfides[k].second];
    std::snprintf(b, sizeof b, "SSBOND %3zu CYS %1s %4d%c   CYS %1s %4d%c\n", k + 1, a.chain.empty() ? " " : a.chain.c_str(),
                  a.number, a.icode, c.chain.empty() ? " " : c.chain.c_str(), c.number, c.icode);
    out += b;
  }
  std::string chain = "\x01";
  for (std::size_t i = 0; i < p.residues.size(); ++i) {
    const auto& r = p.residues[i];
    if (chain != "\x01" && r.chain != chain) {
      out += "TER\n";
    }
    chain = r.chain;
    const std::string rn = r.state.empty() ? r.name : r.name;
    for (const auto& a : r.atoms) {
      ++serial;
      // Atom names: four characters, one-letter elements start in column 14.
      std::string nm = a.name;
      if (nm.size() < 4 && a.element.size() == 1) nm = " " + nm;
      char chg[3] = "  ";
      if (a.charge) std::snprintf(chg, sizeof chg, "%d%c", std::abs(a.charge), a.charge > 0 ? '+' : '-');
      std::string el = a.element;
      for (char& ch : el) ch = char(std::toupper(static_cast<unsigned char>(ch)));
      std::snprintf(b, sizeof b, "%-6s%5d %-4.4s %3.3s %1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
                    r.hetero ? "HETATM" : "ATOM", serial % 100000, nm.c_str(), rn.c_str(), r.chain.empty() ? " " : r.chain.c_str(),
                    r.number, r.icode, a.x, a.y, a.z, a.occupancy, a.bfactor, el.c_str(), chg);
      out += b;
    }
  }
  out += "TER\nEND\n";
  return out;
}

namespace {

struct SeqCounts {
  std::map<char, int> n;
  char first = 0, last = 0;
};
SeqCounts counts(std::string_view s) {
  SeqCounts c;
  for (char ch : s) {
    const char u = char(std::toupper(static_cast<unsigned char>(ch)));
    if (u < 'A' || u > 'Z') continue;
    if (!c.first) c.first = u;
    c.last = u;
    ++c.n[u];
  }
  return c;
}

}  // namespace

double sequence_net_charge(std::string_view sequence, double pH, PkaSet set) {
  const SeqCounts c = counts(sequence);
  if (!c.first) return 0.0;
  auto pos = [&](double pk) { return 1.0 / (1.0 + std::pow(10.0, pH - pk)); };
  auto neg = [&](double pk) { return 1.0 / (1.0 + std::pow(10.0, pk - pH)); };
  auto cnt = [&](char k) {
    auto it = c.n.find(k);
    return it == c.n.end() ? 0 : it->second;
  };
  double q = 0.0;
  if (set == PkaSet::Bjellqvist) {
    // ExPASy / Bjellqvist: terminal pKa depend on the terminal residue.
    static const std::map<char, double> nterm = {{'A', 7.59}, {'M', 7.0}, {'S', 6.93}, {'P', 8.36}, {'T', 6.82}, {'V', 7.44}, {'E', 7.7}};
    static const std::map<char, double> cterm = {{'D', 4.55}, {'E', 4.75}};
    const double pn = nterm.count(c.first) ? nterm.at(c.first) : 7.5;
    const double pc = cterm.count(c.last) ? cterm.at(c.last) : 3.55;
    q += pos(pn) + cnt('K') * pos(10.0) + cnt('R') * pos(12.0) + cnt('H') * pos(5.98);
    q -= neg(pc) + cnt('D') * neg(4.05) + cnt('E') * neg(4.45) + cnt('C') * neg(9.0) + cnt('Y') * neg(10.0);
  } else {
    const Pka K;
    q += pos(K.nterm) + cnt('K') * pos(K.lys) + cnt('R') * pos(K.arg) + cnt('H') * pos(K.his);
    q -= neg(K.cterm) + cnt('D') * neg(K.asp) + cnt('E') * neg(K.glu) + cnt('C') * neg(K.cys) + cnt('Y') * neg(K.tyr);
  }
  return q;
}

double isoelectric_point(std::string_view sequence, PkaSet set) {
  double lo = 0.0, hi = 14.0;
  for (int i = 0; i < 100 && hi - lo > 1e-6; ++i) {
    const double m = 0.5 * (lo + hi);
    (sequence_net_charge(sequence, m, set) > 0 ? lo : hi) = m;
  }
  return 0.5 * (lo + hi);
}

double sequence_molecular_weight(std::string_view sequence) {
  const double mC = find_element_by_z(6)->mass, mH = find_element_by_z(1)->mass, mN = find_element_by_z(7)->mass,
               mO = find_element_by_z(8)->mass, mS = find_element_by_z(16)->mass;
  double w = 2 * mH + mO;  // one water
  bool any = false;
  for (char ch : sequence) {
    const char u = char(std::toupper(static_cast<unsigned char>(ch)));
    for (const auto& [name, d] : residue_defs())
      if (d.code == u) {
        const auto& f = d.formula;
        w += f[0] * mC + f[1] * mH + f[2] * mN + f[3] * mO + f[4] * mS;
        any = true;
      }
  }
  return any ? w : 0.0;
}

std::string protein_json(std::string_view pdb_text, double pH) {
  Diagnostics d;
  auto p = read_protein_pdb(pdb_text, d);
  std::string o = std::string("{\"ok\":") + (p ? "true" : "false");
  auto jn = [](double v) {
    char b[40];
    std::snprintf(b, sizeof b, "%.10g", v);
    return std::isfinite(v) ? std::string(b) : std::string("null");
  };
  if (p) {
    ProtonationOptions opt;
    opt.pH = pH;
    const auto rep = protonate(*p, opt, d);
    o += ",\"title\":" + json_escape(p->title) + ",\"chains\":[";
    std::string all;
    const auto seqs = chain_sequences(*p);
    for (std::size_t i = 0; i < seqs.size(); ++i) {
      std::size_t n = 0;
      for (const auto& r : p->residues) n += r.chain == seqs[i].first;
      o += (i ? "," : "") + std::string("{\"id\":") + json_escape(seqs[i].first) + ",\"sequence\":" + json_escape(seqs[i].second) +
           ",\"residues\":" + std::to_string(n) + ",\"pI\":" + jn(isoelectric_point(seqs[i].second)) +
           ",\"molecular_weight\":" + jn(sequence_molecular_weight(seqs[i].second)) + "}";
      all += seqs[i].second;
    }
    o += "],\"residues\":[";
    bool first = true;
    static const std::set<std::string> titratable = {"ASP", "GLU", "HIS", "LYS", "ARG", "CYS", "TYR"};
    for (const auto& r : p->residues) {
      if (!titratable.count(canonical_name(r.name)) && !r.n_terminal && !r.c_terminal) continue;
      if (!residue_defs().count(canonical_name(r.name))) continue;
      o += (first ? "" : ",") + std::string("{\"chain\":") + json_escape(r.chain) + ",\"number\":" + std::to_string(r.number) +
           ",\"icode\":" + json_escape(std::string(r.icode == ' ' ? "" : std::string(1, r.icode))) + ",\"name\":" +
           json_escape(r.name) + ",\"state\":" + json_escape(r.state) + ",\"charge\":" + std::to_string(r.charge) +
           ",\"n_terminal\":" + (r.n_terminal ? "true" : "false") + ",\"c_terminal\":" + (r.c_terminal ? "true" : "false") + "}";
      first = false;
    }
    o += "],\"disulfides\":[";
    for (std::size_t k = 0; k < p->disulfides.size(); ++k) {
      const auto& a = p->residues[p->disulfides[k].first];
      const auto& b = p->residues[p->disulfides[k].second];
      o += (k ? ",[" : "[") + json_escape(a.chain) + "," + std::to_string(a.number) + "," + json_escape(b.chain) + "," +
           std::to_string(b.number) + "]";
    }
    o += "],\"chain_breaks\":[";
    for (std::size_t k = 0; k < p->chain_breaks.size(); ++k) {
      const auto& r = p->residues[p->chain_breaks[k]];
      o += (k ? "," : "") + json_escape(r.chain + " " + std::to_string(r.number));
    }
    o += "],\"missing_atoms\":[";
    for (std::size_t k = 0; k < rep.missing_atoms.size(); ++k) o += (k ? "," : "") + json_escape(rep.missing_atoms[k]);
    o += "],\"unknown_residues\":[";
    for (std::size_t k = 0; k < rep.unknown_residues.size(); ++k) o += (k ? "," : "") + json_escape(rep.unknown_residues[k]);
    o += "],\"atoms\":" + std::to_string(p->atom_count()) + ",\"hydrogens_added\":" + std::to_string(rep.hydrogens_added) +
         ",\"hydrogens_removed\":" + std::to_string(rep.hydrogens_removed) + ",\"net_charge\":" + std::to_string(rep.net_charge) +
         ",\"pI\":" + jn(all.empty() ? 0.0 : isoelectric_point(all)) + ",\"titration\":[";
    for (int k = 0; k <= 56; ++k) {
      const double ph = k * 0.25;
      o += (k ? ",[" : "[") + jn(ph) + "," + jn(sequence_net_charge(all, ph)) + "]";
    }
    o += "],\"pdb\":" + json_escape(write_protein_pdb(*p));
  }
  return o + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

}  // namespace chemlarp
