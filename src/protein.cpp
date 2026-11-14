#include "chemlarp/protein.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <numbers>
#include <set>
#include <sstream>

#include "chemlarp/elements.hpp"
#include "chemlarp/json.hpp"
#include "protein_internal.hpp"

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

using namespace detail;


namespace {

// One coordinate record of either file format.
struct RawAtom {
  ProteinAtom atom;
  std::string resname, chain;
  int number = 0;
  char icode = ' ';
  char alt = ' ';
  bool hetatm = false;
};

struct RawFile {
  std::string title;
  std::vector<RawAtom> atoms;                    // first model only
  std::vector<std::array<std::string, 4>> ssbond;  // chain, number, chain, number
  std::vector<MissingResidue> missing;
};

void set_element(ProteinAtom& a, std::string el, const std::string& name_field) {
  if (el.empty()) {  // from the atom name: the first letter
    for (char c : name_field)
      if (std::isalpha(static_cast<unsigned char>(c))) {
        el = std::string(1, c);
        break;
      }
  }
  if (!el.empty()) el[0] = char(std::toupper(static_cast<unsigned char>(el[0])));
  if (el.size() == 2) el[1] = char(std::tolower(static_cast<unsigned char>(el[1])));
  a.element = el;
}

// Residues, alternate locations, disulfides, termini and chain breaks from the raw records.
std::optional<Protein> assemble(RawFile& raw, Diagnostics& diag) {
  Protein p;
  p.title = raw.title;
  p.missing_residues = std::move(raw.missing);
  std::map<std::string, std::size_t> index;  // residue key -> index
  std::map<std::pair<std::size_t, std::string>, std::pair<ProteinAtom, double>> chosen;  // (res, name) -> atom, occupancy
  std::vector<std::pair<std::size_t, std::string>> order;
  for (auto& ra : raw.atoms) {
    const std::string key = ra.chain + "|" + std::to_string(ra.number) + "|" + std::string(1, ra.icode) + "|" + ra.resname;
    auto it = index.find(key);
    if (it == index.end()) {
      ProteinResidue r;
      r.name = ra.resname;
      r.chain = ra.chain;
      r.number = ra.number;
      r.icode = ra.icode;
      r.hetero = ra.hetatm && !residue_defs().count(canonical_name(ra.resname));
      p.residues.push_back(r);
      it = index.emplace(key, p.residues.size() - 1).first;
    }
    const auto k = std::make_pair(it->second, ra.atom.name);
    auto c = chosen.find(k);
    if (c == chosen.end()) {
      chosen.emplace(k, std::make_pair(ra.atom, ra.atom.occupancy));
      order.push_back(k);
    } else if (ra.alt != ' ' && ra.atom.occupancy > c->second.second) {
      c->second = {ra.atom, ra.atom.occupancy};
    }
  }
  for (const auto& k : order) p.residues[k.first].atoms.push_back(chosen[k].first);
  if (p.residues.empty()) {
    diag.error(DiagCode::UnknownFileFormat, "no atom records");
    return std::nullopt;
  }
  // Disulfides: SSBOND records and S–S distance.
  std::set<std::pair<std::size_t, std::size_t>> ss;
  std::vector<std::size_t> cys;
  for (std::size_t i = 0; i < p.residues.size(); ++i)
    if (canonical_name(p.residues[i].name) == "CYS" && p.residues[i].atom("SG")) cys.push_back(i);
  // Cell list on SG positions (large structures have thousands of cysteines).
  {
    std::map<std::array<int, 3>, std::vector<std::size_t>> cells;
    auto cell = [](const ProteinAtom& a) { return std::array<int, 3>{int(std::floor(a.x / 3.0)), int(std::floor(a.y / 3.0)), int(std::floor(a.z / 3.0))}; };
    for (std::size_t i : cys) cells[cell(*p.residues[i].atom("SG"))].push_back(i);
    for (std::size_t i : cys) {
      const ProteinAtom& s = *p.residues[i].atom("SG");
      const auto c = cell(s);
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
          for (int dz = -1; dz <= 1; ++dz) {
            auto it = cells.find({c[0] + dx, c[1] + dy, c[2] + dz});
            if (it == cells.end()) continue;
            for (std::size_t j : it->second)
              if (j > i && dist(s, *p.residues[j].atom("SG")) < 2.5) ss.insert({i, j});
          }
    }
  }
  for (const auto& s : raw.ssbond) {
    std::size_t ia = SIZE_MAX, ib = SIZE_MAX;
    for (std::size_t i : cys) {
      const auto& r = p.residues[i];
      if (r.chain == s[0] && std::to_string(r.number) == s[1]) ia = i;
      if (r.chain == s[2] && std::to_string(r.number) == s[3]) ib = i;
    }
    if (ia != SIZE_MAX && ib != SIZE_MAX && ia != ib) ss.insert({std::min(ia, ib), std::max(ia, ib)});
  }
  p.disulfides.assign(ss.begin(), ss.end());
  // Termini and chain breaks among amino-acid residues of each chain.
  std::map<std::string, std::size_t> last_in_chain;  // previous amino-acid residue of the chain
  std::vector<std::size_t> prev_of(p.residues.size(), SIZE_MAX);
  for (std::size_t i = 0; i < p.residues.size(); ++i) {
    if (!residue_defs().count(canonical_name(p.residues[i].name))) continue;
    auto it = last_in_chain.find(p.residues[i].chain);
    if (it != last_in_chain.end()) prev_of[i] = it->second;
    last_in_chain[p.residues[i].chain] = i;
  }
  for (std::size_t i = 0; i < p.residues.size(); ++i) {
    auto& r = p.residues[i];
    if (!residue_defs().count(canonical_name(r.name))) continue;
    r.n_terminal = true;
    r.c_terminal = true;
  }
  for (std::size_t i = 0; i < p.residues.size(); ++i) {
    const std::size_t prev = prev_of[i];
    if (prev == SIZE_MAX) continue;
    const ProteinAtom* N = p.residues[i].atom("N");
    const ProteinAtom* C = p.residues[prev].atom("C");
    if (N && C && dist(*C, *N) < 2.0) {
      p.residues[i].n_terminal = false;
      p.residues[prev].c_terminal = false;
    } else {
      p.chain_breaks.push_back(i);
    }
  }
  return p;
}

// ---- mmCIF ----------------------------------------------------------------------------------------------

struct CifTable {
  std::vector<std::string> columns;
  std::vector<std::string_view> cells;  // row-major
  std::size_t rows() const { return columns.empty() ? 0 : cells.size() / columns.size(); }
  int column(std::string_view name) const {
    for (std::size_t i = 0; i < columns.size(); ++i)
      if (columns[i] == name) return int(i);
    return -1;
  }
};

// The categories named in `wanted` (e.g. "_atom_site") of the first data block. Values keep their text;
// "." and "?" (not given) become empty.
std::map<std::string, CifTable> parse_cif(std::string_view t, const std::set<std::string>& wanted) {
  std::map<std::string, CifTable> out;
  std::size_t i = 0;
  const std::size_t n = t.size();
  bool quoted = false;
  auto next = [&](std::string_view& tok) -> bool {
    for (;;) {
      while (i < n && (t[i] == ' ' || t[i] == '\t' || t[i] == '\r' || t[i] == '\n')) ++i;
      if (i >= n) return false;
      if (t[i] == '#') {
        while (i < n && t[i] != '\n') ++i;
        continue;
      }
      break;
    }
    quoted = false;
    if (t[i] == ';' && (i == 0 || t[i - 1] == '\n')) {  // text field up to a line starting with ';'
      const std::size_t b = i + 1;
      std::size_t e = t.find("\n;", b);
      if (e == std::string_view::npos) e = n;
      tok = t.substr(b, e - b);
      i = std::min(n, e + 2);
      quoted = true;
      return true;
    }
    if (t[i] == '\'' || t[i] == '"') {
      const char q = t[i];
      const std::size_t b = ++i;
      while (i < n && !(t[i] == q && (i + 1 >= n || t[i + 1] == ' ' || t[i + 1] == '\t' || t[i + 1] == '\r' || t[i + 1] == '\n')) && t[i] != '\n') ++i;
      tok = t.substr(b, i - b);
      if (i < n) ++i;
      quoted = true;
      return true;
    }
    const std::size_t b = i;
    while (i < n && t[i] != ' ' && t[i] != '\t' && t[i] != '\r' && t[i] != '\n') ++i;
    tok = t.substr(b, i - b);
    return true;
  };
  auto value = [&](std::string_view v) { return !quoted && (v == "." || v == "?") ? std::string_view{} : v; };
  auto category = [](std::string_view name) { return std::string(name.substr(0, name.find('.'))); };
  auto item = [](std::string_view name) {
    const auto d = name.find('.');
    return std::string(d == std::string_view::npos ? std::string_view{} : name.substr(d + 1));
  };
  std::string_view tok;
  bool have = next(tok);
  int blocks = 0;
  while (have) {
    if (!quoted && tok.substr(0, 5) == "data_") {
      if (++blocks > 1) break;
      have = next(tok);
    } else if (!quoted && tok == "loop_") {
      std::vector<std::string_view> names;
      while ((have = next(tok)) && !quoted && !tok.empty() && tok[0] == '_') names.push_back(tok);
      CifTable* tab = nullptr;
      if (!names.empty() && wanted.count(category(names[0]))) {
        tab = &out[category(names[0])];
        for (auto nm : names) tab->columns.push_back(item(nm));
      }
      while (have && (quoted || !(tok.empty() || tok[0] == '_' || tok == "loop_" || tok.substr(0, 5) == "data_"))) {
        if (tab) tab->cells.push_back(value(tok));
        have = next(tok);
      }
      if (tab && !tab->columns.empty()) tab->cells.resize(tab->cells.size() - tab->cells.size() % tab->columns.size());
    } else if (!quoted && !tok.empty() && tok[0] == '_') {
      const std::string_view name = tok;
      have = next(tok);
      if (!have) break;
      if (wanted.count(category(name))) {
        auto& tab = out[category(name)];
        tab.columns.push_back(item(name));
        tab.cells.push_back(value(tok));
      }
      have = next(tok);
    } else {
      have = next(tok);
    }
  }
  return out;
}

int to_int(std::string_view s, int fallback = 0) {
  if (s.empty()) return fallback;
  return std::atoi(std::string(s).c_str());
}
double to_double(std::string_view s, double fallback = 0.0) {
  if (s.empty()) return fallback;
  return std::atof(std::string(s).c_str());
}

bool looks_like_cif(std::string_view text) {
  const auto a = text.find_first_not_of(" \t\r\n");
  if (a == std::string_view::npos) return false;
  if (text.substr(a, 5) == "data_") return true;
  return text.find("_atom_site.") != std::string_view::npos && text.find("\nATOM  ") == std::string_view::npos && text.substr(0, 6) != "ATOM  ";
}

}  // namespace

std::optional<Protein> read_protein_pdb(std::string_view text, Diagnostics& diag) {
  RawFile raw;
  std::istringstream in{std::string(text)};
  bool in_model = false, model_done = false;
  for (std::string line; std::getline(in, line);) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const std::string rec = line.substr(0, 6);
    if (rec == "TITLE " && raw.title.size() < 300) raw.title += (raw.title.empty() ? "" : " ") + trim(line.size() > 10 ? line.substr(10) : "");
    if (rec == "MODEL ") {
      if (in_model || model_done) model_done = true;
      in_model = true;
    }
    if (rec == "ENDMDL") {
      model_done = true;
      in_model = false;
    }
    if (rec == "SSBOND" && line.size() >= 35)
      raw.ssbond.push_back({trim(line.substr(15, 1)), trim(line.substr(17, 5)), trim(line.substr(29, 1)), trim(line.substr(31, 5))});
    if (rec == "REMARK" && line.size() >= 26 && line.substr(6, 4) == " 465") {
      // "REMARK 465   M RES C SSSEQI": optional model, residue name, chain, number, insertion code.
      const std::string res = trim(line.substr(15, 3));
      const std::string num = trim(line.substr(21, 5));
      const std::string model = trim(line.substr(10, 4));
      const bool numeric = !num.empty() && (std::isdigit(static_cast<unsigned char>(num.back())) != 0) &&
                           num.find_first_not_of("-0123456789") == std::string::npos;
      if (res.size() == 3 && numeric && (model.empty() || model == "1") && std::isupper(static_cast<unsigned char>(res[0]))) {
        MissingResidue m;
        m.name = res;
        m.chain = trim(line.substr(19, 1));
        m.number = std::atoi(num.c_str());
        m.icode = line.size() > 26 ? line[26] : ' ';
        raw.missing.push_back(m);
      }
    }
    if ((rec != "ATOM  " && rec != "HETATM") || model_done) continue;
    if (line.size() < 54) {
      diag.warn(DiagCode::MalformedRecord, "short ATOM record skipped: " + line);
      continue;
    }
    RawAtom ra;
    ProteinAtom& a = ra.atom;
    a.name = trim(line.substr(12, 4));
    ra.alt = line[16];
    ra.resname = trim(line.substr(17, 3));
    ra.chain = trim(line.substr(21, 1));
    ra.number = std::atoi(line.substr(22, 4).c_str());
    ra.icode = line.size() > 26 ? line[26] : ' ';
    ra.hetatm = rec == "HETATM";
    a.x = std::atof(line.substr(30, 8).c_str());
    a.y = std::atof(line.substr(38, 8).c_str());
    a.z = std::atof(line.substr(46, 8).c_str());
    a.occupancy = line.size() >= 60 ? std::atof(line.substr(54, 6).c_str()) : 1.0;
    a.bfactor = line.size() >= 66 ? std::atof(line.substr(60, 6).c_str()) : 0.0;
    // Without an element column, a name starting in column 13 is a two-letter element only for HETATM.
    set_element(a, line.size() >= 78 ? trim(line.substr(76, 2)) : "", line.substr(12, 4));
    if (line.size() >= 80) {
      const std::string c = trim(line.substr(78, 2));
      if (c.size() == 2 && std::isdigit(static_cast<unsigned char>(c[0]))) a.charge = (c[1] == '-' ? -1 : 1) * (c[0] - '0');
    }
    raw.atoms.push_back(std::move(ra));
  }
  return assemble(raw, diag);
}

std::optional<Protein> read_protein_cif(std::string_view text, Diagnostics& diag) {
  const auto tabs = parse_cif(text, {"_atom_site", "_struct", "_pdbx_poly_seq_scheme", "_struct_conn"});
  RawFile raw;
  if (auto it = tabs.find("_struct"); it != tabs.end()) {
    const int c = it->second.column("title");
    if (c >= 0 && it->second.rows() > 0) {
      raw.title = std::string(it->second.cells[std::size_t(c)]);
      for (char& ch : raw.title)
        if (ch == '\n' || ch == '\r') ch = ' ';
      raw.title = trim(raw.title);
    }
  }
  auto it = tabs.find("_atom_site");
  if (it == tabs.end() || it->second.rows() == 0) {
    diag.error(DiagCode::UnknownFileFormat, "no _atom_site records");
    return std::nullopt;
  }
  const CifTable& t = it->second;
  auto col = [&](const char* a, const char* b = nullptr) {
    int c = t.column(a);
    if (c < 0 && b) c = t.column(b);
    return c;
  };
  const int c_group = col("group_PDB"), c_el = col("type_symbol"), c_name = col("auth_atom_id", "label_atom_id"),
            c_name2 = col("label_atom_id"), c_alt = col("label_alt_id"), c_res = col("auth_comp_id", "label_comp_id"),
            c_chain = col("auth_asym_id", "label_asym_id"), c_seq = col("auth_seq_id", "label_seq_id"),
            c_ins = col("pdbx_PDB_ins_code"), c_x = col("Cartn_x"), c_y = col("Cartn_y"), c_z = col("Cartn_z"),
            c_occ = col("occupancy"), c_b = col("B_iso_or_equiv"), c_q = col("pdbx_formal_charge"),
            c_model = col("pdbx_PDB_model_num");
  if (c_name < 0 || c_res < 0 || c_x < 0 || c_y < 0 || c_z < 0 || c_seq < 0) {
    diag.error(DiagCode::MalformedRecord, "_atom_site lacks atom names, residue names, numbers or coordinates");
    return std::nullopt;
  }
  const std::size_t nc = t.columns.size();
  std::string_view first_model;
  raw.atoms.reserve(t.rows());
  for (std::size_t r = 0; r < t.rows(); ++r) {
    const std::string_view* row = &t.cells[r * nc];
    auto get = [&](int c) { return c < 0 ? std::string_view{} : row[c]; };
    if (c_model >= 0) {
      if (r == 0) first_model = get(c_model);
      else if (get(c_model) != first_model) break;
    }
    RawAtom ra;
    ProteinAtom& a = ra.atom;
    a.name = std::string(get(c_name2).empty() ? get(c_name) : get(c_name2));
    ra.alt = get(c_alt).empty() ? ' ' : get(c_alt)[0];
    ra.resname = std::string(get(c_res));
    ra.chain = std::string(get(c_chain));
    ra.number = to_int(get(c_seq));
    ra.icode = get(c_ins).empty() ? ' ' : get(c_ins)[0];
    ra.hetatm = get(c_group) == "HETATM";
    a.x = to_double(get(c_x));
    a.y = to_double(get(c_y));
    a.z = to_double(get(c_z));
    a.occupancy = to_double(get(c_occ), 1.0);
    a.bfactor = to_double(get(c_b));
    a.charge = to_int(get(c_q));
    set_element(a, std::string(get(c_el)), a.name);
    raw.atoms.push_back(std::move(ra));
  }
  if (auto sc = tabs.find("_struct_conn"); sc != tabs.end()) {
    const CifTable& s = sc->second;
    const int ty = s.column("conn_type_id"), a1 = s.column("ptnr1_auth_asym_id"), n1 = s.column("ptnr1_auth_seq_id"),
              a2 = s.column("ptnr2_auth_asym_id"), n2 = s.column("ptnr2_auth_seq_id");
    if (ty >= 0 && a1 >= 0 && n1 >= 0 && a2 >= 0 && n2 >= 0)
      for (std::size_t r = 0; r < s.rows(); ++r) {
        const std::string_view* row = &s.cells[r * s.columns.size()];
        if (row[ty] == "disulf") raw.ssbond.push_back({std::string(row[a1]), std::string(row[n1]), std::string(row[a2]), std::string(row[n2])});
      }
  }
  if (auto ps = tabs.find("_pdbx_poly_seq_scheme"); ps != tabs.end()) {
    // A residue of the sequence with no auth_seq_num has no coordinates.
    const CifTable& s = ps->second;
    const int mon = s.column("mon_id"), pdbn = s.column("pdb_seq_num"), auth = s.column("auth_seq_num"),
              strand = s.column("pdb_strand_id"), ins = s.column("pdb_ins_code");
    if (mon >= 0 && pdbn >= 0 && auth >= 0 && strand >= 0)
      for (std::size_t r = 0; r < s.rows(); ++r) {
        const std::string_view* row = &s.cells[r * s.columns.size()];
        if (!row[auth].empty()) continue;
        MissingResidue m;
        m.name = std::string(row[mon]);
        m.chain = std::string(row[strand]);
        m.number = to_int(row[pdbn]);
        m.icode = ins >= 0 && !row[ins].empty() ? row[ins][0] : ' ';
        raw.missing.push_back(m);
      }
  }
  return assemble(raw, diag);
}

std::optional<Protein> read_protein(std::string_view text, Diagnostics& diag) {
  return looks_like_cif(text) ? read_protein_cif(text, diag) : read_protein_pdb(text, diag);
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

std::string write_protein_pdb(const Protein& p) {
  std::string out;
  char b[128];
  if (!p.title.empty()) out += "REMARK   1 " + p.title.substr(0, 68) + "\n";
  out += "REMARK   2 WRITTEN BY CHEMLARP\n";
  {
    std::size_t n = 0;
    for (const auto& r : p.residues) n += r.predicted;
    if (n) out += "REMARK   3 " + std::to_string(n) + " RESIDUES WITHOUT EXPERIMENTAL COORDINATES WERE BUILT (OCCUPANCY 0.00)\n";
  }
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
                    r.number, r.icode, a.x, a.y, a.z, r.predicted ? 0.0 : a.occupancy, a.bfactor, el.c_str(), chg);
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

bool fits_pdb_format(const Protein& p) {
  if (p.atom_count() > 99999) return false;
  for (const auto& r : p.residues)
    if (r.chain.size() > 1 || r.number > 9999 || r.number < -999) return false;
  return true;
}

std::string write_protein_cif(const Protein& p) {
  auto q = [](const std::string& s) {
    if (s.empty()) return std::string("?");
    if (s.find('\'') != std::string::npos) return "\"" + s + "\"";
    if (s.find('"') != std::string::npos || s.find(' ') != std::string::npos || s[0] == '_' || s[0] == '#') return "'" + s + "'";
    return s;
  };
  std::string out = "data_chemlarp\n#\n";
  if (!p.title.empty()) {
    std::string t = p.title;
    for (char& c : t)
      if (c == '\n' || c == '\r') c = ' ';
    out += "_struct.title\n;" + t + "\n;\n#\n";
  }
  out += "loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n"
         "_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n"
         "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n"
         "_atom_site.B_iso_or_equiv\n_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n"
         "_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\n";
  char b[256];
  std::size_t serial = 0;
  for (const auto& r : p.residues) {
    const std::string chain = q(r.chain.empty() ? std::string("A") : r.chain), comp = q(r.name);
    const std::string ins = r.icode == ' ' ? std::string("?") : std::string(1, r.icode);
    for (const auto& a : r.atoms) {
      std::string el = a.element;
      for (char& ch : el) ch = char(std::toupper(static_cast<unsigned char>(ch)));
      const std::string name = q(a.name);
      std::snprintf(b, sizeof b, "%s %zu %s %s . %s %s 1 %d %s %.3f %.3f %.3f %.2f %.2f %s %d %s %s %s 1\n", r.hetero ? "HETATM" : "ATOM",
                    ++serial, el.empty() ? "X" : el.c_str(), name.c_str(), comp.c_str(), chain.c_str(), r.number, ins.c_str(), a.x, a.y, a.z,
                    a.occupancy, a.bfactor, a.charge ? std::to_string(a.charge).c_str() : "?", r.number, comp.c_str(), chain.c_str(),
                    name.c_str());
      out += b;
    }
  }
  return out + "#\n";
}

std::string protein_json(std::string_view pdb_text, double pH) { return protein_json(pdb_text, ProteinJsonOptions{pH}); }

std::string protein_json(std::string_view text, const ProteinJsonOptions& jo) {
  Diagnostics d;
  auto p = read_protein(text, d);
  std::string o = std::string("{\"ok\":") + (p ? "true" : "false");
  auto jn = [](double v) {
    char b[40];
    std::snprintf(b, sizeof b, "%.10g", v);
    return std::isfinite(v) ? std::string(b) : std::string("null");
  };
  auto list = [](const std::vector<std::string>& v) {
    std::string s = "[";
    for (std::size_t k = 0; k < v.size(); ++k) s += (k ? "," : "") + json_escape(v[k]);
    return s + "]";
  };
  if (p) {
    RepairReport fix;
    const std::size_t listed_missing = p->missing_residues.size();
    if (jo.repair) {
      RepairOptions ro;
      ro.loops = jo.loops;
      ro.termini = jo.termini;
      fix = repair(*p, ro, d);
    }
    ProtonationOptions opt;
    opt.pH = jo.pH;
    opt.flips = jo.flips;
    opt.structure_pka = jo.structure_pka;
    const auto rep = protonate(*p, opt, d);
    const auto& pka = rep.pka;
    o += ",\"title\":" + json_escape(p->title) + ",\"chains\":[";
    std::string all;
    const auto seqs = chain_sequences(*p);
    std::map<std::string, std::size_t> chain_count;
    for (const auto& r : p->residues) ++chain_count[r.chain];
    for (std::size_t i = 0; i < seqs.size(); ++i) {
      o += (i ? "," : "") + std::string("{\"id\":") + json_escape(seqs[i].first) + ",\"sequence\":" + json_escape(seqs[i].second) +
           ",\"residues\":" + std::to_string(chain_count[seqs[i].first]) + ",\"pI\":" + jn(isoelectric_point(seqs[i].second)) +
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
           ",\"n_terminal\":" + (r.n_terminal ? "true" : "false") + ",\"c_terminal\":" + (r.c_terminal ? "true" : "false") +
           ",\"predicted\":" + (r.predicted ? "true" : "false") + "}";
      first = false;
    }
    o += "],\"pka\":[";
    for (std::size_t k = 0; k < pka.size(); ++k) {
      const auto& s = pka[k];
      const auto& r = p->residues[s.residue];
      o += (k ? "," : "") + std::string("{\"chain\":") + json_escape(r.chain) + ",\"number\":" + std::to_string(r.number) +
           ",\"name\":" + json_escape(r.name) + ",\"group\":" + json_escape(s.group) + ",\"reference\":" + jn(s.reference) +
           ",\"pKa\":" + jn(s.pka) + ",\"desolvation\":" + jn(s.desolvation) + ",\"hydrogen_bonds\":" + jn(s.hydrogen_bonds) +
           ",\"charges\":" + jn(s.charges) + ",\"buried\":" + jn(s.buried) + ",\"features\":[";
      for (int i = 0; i < kPkaFeatures; ++i) o += (i ? "," : "") + jn(s.features[i]);
      o += "]}";
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
    o += "],\"missing_atoms\":" + list(rep.missing_atoms) + ",\"unknown_residues\":" + list(rep.unknown_residues) +
         ",\"flipped\":" + list(rep.flipped) + ",\"undetermined\":" + list(rep.undetermined) +
         ",\"repair\":{\"atoms_added\":" + std::to_string(fix.atoms_added) + ",\"completed_residues\":" + list(fix.completed_residues) +
         ",\"built_residues\":" + list(fix.built_residues) + ",\"not_built\":" + list(fix.not_built) +
         ",\"missing_residues_listed\":" + std::to_string(listed_missing) + ",\"missing_residues_left\":" +
         std::to_string(p->missing_residues.size()) + "}";
    o += ",\"atoms\":" + std::to_string(p->atom_count()) + ",\"hydrogens_added\":" + std::to_string(rep.hydrogens_added) +
         ",\"hydrogens_removed\":" + std::to_string(rep.hydrogens_removed) + ",\"net_charge\":" + std::to_string(rep.net_charge) +
         ",\"pI\":" + jn(all.empty() ? 0.0 : isoelectric_point(all)) + ",\"titration\":[";
    for (int k = 0; k <= 56; ++k) {
      const double ph = k * 0.25;
      o += (k ? ",[" : "[") + jn(ph) + "," + jn(sequence_net_charge(all, ph)) + "]";
    }
    // Net charge against pH from the structure-based pKa values (Henderson–Hasselbalch for each group).
    o += "],\"titration_structure\":[";
    for (int k = 0; k <= 56; ++k) {
      const double ph = k * 0.25;
      double q = 0.0;
      for (const auto& s : pka) {
        const bool acid = s.group == "ASP" || s.group == "GLU" || s.group == "CYS" || s.group == "TYR" || s.group == "C-";
        q += acid ? -1.0 / (1.0 + std::pow(10.0, s.pka - ph)) : 1.0 / (1.0 + std::pow(10.0, ph - s.pka));
      }
      o += (k ? ",[" : "[") + jn(ph) + "," + jn(q) + "]";
    }
    const bool pdb_ok = fits_pdb_format(*p);
    o += "],\"format\":" + std::string(pdb_ok ? "\"pdb\"" : "\"cif\"") + ",\"pdb\":" + json_escape(pdb_ok ? write_protein_pdb(*p) : std::string()) +
         ",\"cif\":" + json_escape(pdb_ok ? std::string() : write_protein_cif(*p));
  }
  return o + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

}  // namespace chemlarp
