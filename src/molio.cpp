#include "omnichem/molio.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "omnichem/canonical.hpp"
#include "omnichem/cip.hpp"
#include "omnichem/elements.hpp"
#include "omnichem/smiles.hpp"

namespace omnichem {
namespace {

// ---- raw record ------------------------------------------------------------------

enum : std::uint8_t {
  kStereoNone = 0,
  kStereoWedge = 1,       // V2000 1, V3000 CFG=1, CML W
  kStereoEitherDb = 3,    // V2000 3 on a double bond: cis or trans unknown
  kStereoEither = 4,      // V2000 4, V3000 CFG=2 on a single bond: wavy
  kStereoHash = 6,        // V2000 6, V3000 CFG=3, CML H
};

struct RAtom {
  std::uint8_t z = 0;
  int charge = 0;
  std::uint16_t isotope = 0;
  int radical = 0;   // MDL code: 1 singlet, 2 doublet, 3 triplet
  int hcount = -1;   // total hydrogens when the file states it (CML hydrogenCount)
  int valence = -1;  // total valence when the file states it (V2000 vvv, V3000 VAL)
  double x = 0, y = 0, z3 = 0;
  std::size_t line = 0;
};

struct RBond {
  std::uint32_t a = 0, b = 0;  // a is the narrow end of a wedge
  std::uint8_t order = 1;      // 1, 2, 3, or 4 = aromatic
  std::uint8_t stereo = kStereoNone;
  std::size_t line = 0;
};

struct RawMol {
  std::vector<RAtom> atoms;
  std::vector<RBond> bonds;
  bool implicit_h = true;  // fill valences with implicit hydrogens (MDL convention)
  bool failed = false;     // the reader already reported an error
};

std::string at_line(std::size_t line) { return "line " + std::to_string(line) + ": "; }

bool organic(std::uint8_t z) {
  switch (z) {
    case 5: case 6: case 7: case 8: case 9: case 15: case 16: case 17: case 35: case 53: return true;
    default: return false;
  }
}

bool aromatic_organic(std::uint8_t z) {
  return z == 5 || z == 6 || z == 7 || z == 8 || z == 15 || z == 16;
}

int radical_electrons(int code) { return code == 2 ? 1 : (code == 1 || code == 3) ? 2 : 0; }

struct Vec3 {
  double x = 0, y = 0, z = 0;
};
Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 cross(Vec3 a, Vec3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
double norm(Vec3 a) { return std::sqrt(dot(a, a)); }
Vec3 unit(Vec3 a) {
  const double n = norm(a);
  return n > 1e-12 ? Vec3{a.x / n, a.y / n, a.z / n} : Vec3{};
}

void ring_label(std::string& out, int label) {
  if (label < 10) {
    out += static_cast<char>('0' + label);
  } else {
    out += '%';
    out += static_cast<char>('0' + label / 10);
    out += static_cast<char>('0' + label % 10);
  }
}

// Smallest ring through state bond `b` (0 if it is not in a ring).
int ring_size_through(const ChemicalSystemState& st, std::uint32_t b) {
  if (!(st.bonds.flags[b] & kBondInRing)) return 0;
  const std::uint32_t from = st.bonds.begin[b], to = st.bonds.end[b];
  std::unordered_map<std::uint32_t, int> dist{{from, 0}};
  std::vector<std::uint32_t> queue{from};
  for (std::size_t head = 0; head < queue.size(); ++head) {
    const std::uint32_t v = queue[head];
    const int d = dist[v];
    if (d >= 8) break;
    const auto nb = st.adjacency.neighbors(v);
    const auto bs = st.adjacency.bonds_of(v);
    for (std::size_t k = 0; k < nb.size(); ++k) {
      if (bs[k] == b || dist.count(nb[k])) continue;
      if (nb[k] == to) return d + 2;
      dist[nb[k]] = d + 1;
      queue.push_back(nb[k]);
    }
  }
  return 99;
}

// ---- raw record -> species --------------------------------------------------------

class Builder {
 public:
  Builder(ChemicalSystemState& st, const RawMol& m, bool keep_h, Diagnostics& diag)
      : st_(st), m_(m), keep_h_(keep_h), diag_(diag) {}

  std::int32_t build(int dims);

 private:
  bool spell(std::string& smiles);
  void add_stereo(std::uint32_t sp, int dims);
  std::uint32_t global(std::uint32_t raw) const { return base_ + index_[raw]; }

  ChemicalSystemState& st_;
  const RawMol& m_;
  bool keep_h_;
  Diagnostics& diag_;
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> adj_;  // raw atom -> (neighbour, bond)
  std::vector<char> folded_;
  std::vector<std::uint32_t> index_;  // raw atom -> species-local atom
  std::vector<int> implicit_;         // implicit hydrogens added to each raw atom
  std::uint32_t base_ = 0;
};

bool Builder::spell(std::string& smiles) {
  const std::size_t n = m_.atoms.size();
  adj_.assign(n, {});
  for (std::uint32_t k = 0; k < m_.bonds.size(); ++k) {
    adj_[m_.bonds[k].a].push_back({m_.bonds[k].b, k});
    adj_[m_.bonds[k].b].push_back({m_.bonds[k].a, k});
  }
  folded_.assign(n, 0);
  if (!keep_h_) {
    for (std::size_t i = 0; i < n; ++i) {
      const RAtom& a = m_.atoms[i];
      if (a.z != 1 || a.isotope || a.charge || a.radical || adj_[i].size() != 1) continue;
      const RBond& b = m_.bonds[adj_[i][0].second];
      if (b.order == 1 && m_.atoms[adj_[i][0].first].z != 1) folded_[i] = 1;
    }
  }
  index_.assign(n, 0);
  implicit_.assign(n, 0);
  std::uint32_t kept = 0;
  for (std::size_t i = 0; i < n; ++i)
    if (!folded_[i]) index_[i] = kept++;

  std::vector<int> label_of(m_.bonds.size(), -1);
  std::array<bool, 100> used{};
  bool first = true;
  for (std::uint32_t i = 0; i < n; ++i) {
    if (folded_[i]) continue;
    const RAtom& a = m_.atoms[i];
    int sum = 0, fh = 0, kept_h = 0;
    bool arom = false;
    for (auto [nb, k] : adj_[i]) {
      const RBond& b = m_.bonds[k];
      if (folded_[nb]) ++fh;
      else if (m_.atoms[nb].z == 1) ++kept_h;
      if (b.order == 4) {
        arom = true;
        sum += 1;
      } else {
        sum += b.order;
      }
    }
    if (!first) smiles += '.';
    first = false;

    const std::string_view sym = element_symbol(a.z);
    bool bare = false;
    int h = fh;
    if (a.hcount >= 0) {
      h = std::max(0, a.hcount - kept_h);
    } else if (a.valence >= 0) {
      h = fh + std::max(0, a.valence - sum);
    } else if (arom) {
      bare = aromatic_organic(a.z) && !a.charge && !a.isotope && !a.radical && fh == 0;
    } else if (organic(a.z) && m_.implicit_h) {
      bare = !a.charge && !a.isotope && !a.radical && fh == 0;
      const int used_valence = sum + radical_electrons(a.radical);
      for (std::uint8_t v : default_valences(a.z, a.charge)) {
        if (v >= used_valence) { h = fh + (v - used_valence); break; }
      }
    }
    implicit_[i] = h - fh;
    std::string atom;
    if (arom && aromatic_organic(a.z)) {
      atom = std::string(sym);
      atom[0] = static_cast<char>(atom[0] - 'A' + 'a');
    } else {
      atom = std::string(sym);
    }
    if (bare) {
      smiles += atom;
    } else {
      smiles += '[';
      if (a.isotope) smiles += std::to_string(a.isotope);
      smiles += atom;
      if (h > 0) smiles += 'H';
      if (h > 1) smiles += std::to_string(h);
      if (a.charge) {
        smiles += a.charge > 0 ? '+' : '-';
        if (std::abs(a.charge) > 1) smiles += std::to_string(std::abs(a.charge));
      }
      smiles += ']';
    }
    for (auto [nb, k] : adj_[i]) {
      if (folded_[nb]) continue;
      if (label_of[k] >= 0) {  // close
        ring_label(smiles, label_of[k]);
        used[static_cast<std::size_t>(label_of[k])] = false;
        continue;
      }
      int label = 0;
      while (label < 100 && used[static_cast<std::size_t>(label)]) ++label;
      if (label == 100) {
        diag_.error(DiagCode::UnsupportedFileFeature,
                    at_line(a.line) + "more than 100 bonds span this point of the atom order; renumber the atoms");
        return false;
      }
      used[static_cast<std::size_t>(label)] = true;
      label_of[k] = label;
      static constexpr char kSym[] = {'-', '-', '=', '#', ':'};
      smiles += kSym[std::min<int>(m_.bonds[k].order, 4)];
      ring_label(smiles, label);
    }
  }
  return true;
}

std::int32_t Builder::build(int dims) {
  std::string smiles;
  if (!spell(smiles)) return -1;
  Diagnostics local;
  const std::int32_t sp = parse_smiles(st_, smiles, local);
  // Map parser diagnostics (positions in the generated SMILES) back to file lines.
  std::vector<std::uint32_t> raw_of;
  for (std::uint32_t i = 0; i < m_.atoms.size(); ++i)
    if (!folded_[i]) raw_of.push_back(i);
  const std::uint32_t base = sp >= 0 ? st_.species[static_cast<std::size_t>(sp)].atom_begin : 0;
  for (const Diagnostic& d : local.items()) {
    if (d.code == DiagCode::StereoIgnored || d.code == DiagCode::CipUndetermined) continue;
    std::string where;
    std::int64_t local_atom = d.atom;
    if (sp >= 0 && local_atom >= static_cast<std::int64_t>(base)) local_atom -= base;
    if (local_atom >= 0 && static_cast<std::size_t>(local_atom) < raw_of.size()) {
      const std::uint32_t r = raw_of[static_cast<std::size_t>(local_atom)];
      where = at_line(m_.atoms[r].line) + "atom " + std::to_string(r + 1) + ": ";
    } else if (!m_.atoms.empty()) {
      where = at_line(m_.atoms[0].line);
    }
    Diagnostic& out = d.severity == Severity::Error ? diag_.error(d.code, where + d.message)
                                                     : diag_.warn(d.code, where + d.message);
    if (sp >= 0 && d.atom >= 0) out.on_atom(d.atom);
    if (sp >= 0) out.on_species(sp);
  }
  if (sp < 0) return -1;

  base_ = base;
  for (std::uint32_t i = 0; i < m_.atoms.size(); ++i) {
    if (folded_[i]) continue;
    const std::uint32_t g = global(i);
    st_.atoms.x[g] = m_.atoms[i].x;
    st_.atoms.y[g] = m_.atoms[i].y;
    st_.atoms.z[g] = m_.atoms[i].z3;
  }
  if (dims > 0) add_stereo(static_cast<std::uint32_t>(sp), dims);
  Species& s = st_.species[static_cast<std::size_t>(sp)];
  Diagnostics quiet;
  s.canonical_smiles = canonical_smiles(st_, static_cast<std::uint32_t>(sp), quiet);
  s.smiles = s.canonical_smiles;
  ++st_.revision;
  return sp;
}

void Builder::add_stereo(std::uint32_t sp, int dims) {
  const std::size_t tetra_begin = st_.tetrahedral.size();
  const std::size_t db_begin = st_.double_bonds.size();
  auto pos = [&](std::uint32_t r) { return Vec3{m_.atoms[r].x, m_.atoms[r].y, m_.atoms[r].z3}; };

  // Tetrahedral centres.
  for (std::uint32_t i = 0; i < m_.atoms.size(); ++i) {
    if (folded_[i]) continue;
    const RAtom& a = m_.atoms[i];
    const int ih = implicit_[i];
    int hydrogens = ih;
    for (auto [nb, k] : adj_[i]) hydrogens += folded_[nb] ? 1 : 0;
    const int total = static_cast<int>(adj_[i].size()) + ih;
    if (hydrogens > 1) continue;
    bool candidate = false;
    if (total == 4) {
      switch (a.z) {
        case 5: case 6: case 14: case 15: case 16: case 32: case 33: case 34: case 50: candidate = true; break;
        case 7: candidate = a.charge == 1; break;
        default: break;
      }
    } else if (total == 3) {
      candidate = a.z == 15 || a.z == 16 || a.z == 33 || a.z == 34;
    }
    if (!candidate) continue;
    for (auto [nb, k] : adj_[i]) candidate = candidate && m_.bonds[k].order != 4;
    if (!candidate) continue;

    const Vec3 c = pos(i);
    std::vector<std::uint32_t> refs;
    std::vector<Vec3> pts;
    bool wedged = false, wavy = false;
    for (auto [nb, k] : adj_[i]) {
      const RBond& b = m_.bonds[k];
      Vec3 d = pos(nb) - c;
      if (dims == 2) {
        d.z = 0;
        d = unit(d);
        if (b.a == i && (b.stereo == kStereoWedge || b.stereo == kStereoHash)) {
          d.z = b.stereo == kStereoWedge ? 1.0 : -1.0;
          wedged = true;
        }
        if (b.a == i && b.stereo == kStereoEither) wavy = true;
      }
      if (folded_[nb] && dims == 3) continue;  // placed below, opposite the heavy neighbours
      refs.push_back(folded_[nb] ? kImplicitHRef : global(nb));
      pts.push_back(unit(d));
    }
    if (dims == 3 && hydrogens == 1 && ih == 0) {
      // A written hydrogen keeps its place in the neighbour order; its coordinates
      // are not used, so a poorly placed H cannot flip the centre.
      std::size_t slot = 0;
      for (auto [nb, k] : adj_[i]) {
        if (folded_[nb]) break;
        ++slot;
      }
      Vec3 sum{};
      for (const Vec3& p : pts) sum = sum + p;
      if (norm(sum) < 1e-3) continue;
      refs.insert(refs.begin() + static_cast<std::ptrdiff_t>(slot), kImplicitHRef);
      pts.insert(pts.begin() + static_cast<std::ptrdiff_t>(slot), unit(Vec3{-sum.x, -sum.y, -sum.z}));
    }
    if (dims == 2 && (!wedged || wavy)) continue;
    if (pts.size() < 4) {
      Vec3 sum{};
      for (const Vec3& p : pts) sum = sum + p;
      Vec3 missing = unit(Vec3{-sum.x, -sum.y, -sum.z});
      if (norm(sum) < 1e-3) continue;  // planar: no side to put the hydrogen or lone pair on
      refs.push_back(ih == 1 ? kImplicitHRef : kLonePairRef);
      pts.push_back(missing);
    }
    if (pts.size() != 4) continue;
    const double v = dot(pts[1] - pts[0], cross(pts[2] - pts[0], pts[3] - pts[0]));
    if (std::fabs(v) < 0.05) continue;
    TetrahedralStereo t;
    t.atom = global(i);
    for (std::size_t k = 0; k < 4; ++k) t.refs[k] = refs[k];
    t.parity = v < 0 ? Chirality::CounterClockwise : Chirality::Clockwise;
    st_.tetrahedral.push_back(t);
  }

  // Double bonds.
  for (const RBond& b : m_.bonds) {
    if (b.order != 2 || b.stereo == kStereoEitherDb || folded_[b.a] || folded_[b.b]) continue;
    const std::uint32_t ga = global(b.a), gb = global(b.b);
    std::uint32_t bond = UINT32_MAX;
    {
      const auto nb = st_.adjacency.neighbors(ga);
      const auto bs = st_.adjacency.bonds_of(ga);
      for (std::size_t k = 0; k < nb.size(); ++k)
        if (nb[k] == gb) bond = bs[k];
    }
    if (bond == UINT32_MAX || (st_.bonds.flags[bond] & kBondAromatic)) continue;
    const int ring = ring_size_through(st_, bond);
    if (ring > 0 && ring < 8) continue;
    bool skip = false;
    std::array<std::int64_t, 2> ref{-1, -1};
    const std::array<std::uint32_t, 2> ends{b.a, b.b};
    for (int side = 0; side < 2; ++side) {
      const std::uint32_t e = ends[static_cast<std::size_t>(side)], other = ends[static_cast<std::size_t>(1 - side)];
      int count = 0;
      for (auto [nb, k] : adj_[e]) {
        if (nb == other) continue;
        if (m_.bonds[k].order >= 2) skip = true;  // cumulated (allene, azide, ketenimine): no E/Z
        ++count;
        if (m_.bonds[k].stereo == kStereoEither) skip = true;
        if (ref[static_cast<std::size_t>(side)] < 0 && !folded_[nb]) ref[static_cast<std::size_t>(side)] = nb;
      }
      if (count == 0 || count > 2) skip = true;
    }
    if (skip || ref[0] < 0 || ref[1] < 0) continue;
    const std::uint32_t ra = static_cast<std::uint32_t>(ref[0]), rb = static_cast<std::uint32_t>(ref[1]);
    const Vec3 axis = pos(b.b) - pos(b.a);
    const Vec3 n1 = cross(pos(b.a) - pos(ra), axis), n2 = cross(axis, pos(rb) - pos(b.b));
    const double l = norm(n1) * norm(n2);
    if (l < 1e-8) continue;
    const double cosang = dot(n1, n2) / l;
    // Plane normals of (ra, a, b) and (a, b, rb): parallel when ra and rb are cis.
    if (std::fabs(cosang) < 0.17) continue;
    const bool cis = cosang > 0;
    DoubleBondStereo d;
    d.bond = bond;
    const bool forward = st_.bonds.begin[bond] == ga;
    d.ref_begin = forward ? global(ra) : global(rb);
    d.ref_end = forward ? global(rb) : global(ra);
    d.cis = cis;
    st_.double_bonds.push_back(d);
  }

  // Keep only what the CIP rules find stereogenic. Removing one element can
  // change whether another (pseudo-asymmetric) one is, so repeat until stable.
  for (int round = 0; round < 4; ++round) {
    bool changed = false;
    std::vector<TetrahedralStereo> keep_t(st_.tetrahedral.begin(), st_.tetrahedral.begin() + static_cast<std::ptrdiff_t>(tetra_begin));
    for (std::size_t k = tetra_begin; k < st_.tetrahedral.size(); ++k) {
      if (cip_label(st_, st_.tetrahedral[k]) != CipLabel::None) keep_t.push_back(st_.tetrahedral[k]);
      else changed = true;
    }
    std::vector<DoubleBondStereo> keep_d(st_.double_bonds.begin(), st_.double_bonds.begin() + static_cast<std::ptrdiff_t>(db_begin));
    for (std::size_t k = db_begin; k < st_.double_bonds.size(); ++k) {
      if (cip_label(st_, st_.double_bonds[k]) != CipLabel::None) keep_d.push_back(st_.double_bonds[k]);
      else changed = true;
    }
    st_.tetrahedral = std::move(keep_t);
    st_.double_bonds = std::move(keep_d);
    if (!changed) break;
  }
  for (std::size_t k = tetra_begin; k < st_.tetrahedral.size(); ++k)
    st_.atoms.chirality[st_.tetrahedral[k].atom] = static_cast<std::int8_t>(st_.tetrahedral[k].parity);
  (void)sp;
}

std::int32_t build_species(ChemicalSystemState& st, const RawMol& m, bool keep_h, Diagnostics& diag) {
  if (m.failed) return -1;
  if (m.atoms.empty()) {
    diag.error(DiagCode::MalformedRecord, "record has no atoms");
    return -1;
  }
  int dims = 0;
  for (const RAtom& a : m.atoms) {
    if (std::fabs(a.z3) > 1e-4) { dims = 3; break; }
    if (std::fabs(a.x) > 1e-4 || std::fabs(a.y) > 1e-4) dims = 2;
  }
  Builder b(st, m, keep_h, diag);
  return b.build(dims);
}

int record_dimensions(const RawMol& m) {
  int dims = 0;
  for (const RAtom& a : m.atoms) {
    if (std::fabs(a.z3) > 1e-4) return 3;
    if (std::fabs(a.x) > 1e-4 || std::fabs(a.y) > 1e-4) dims = 2;
  }
  return dims;
}

// ---- text helpers -----------------------------------------------------------------

std::vector<std::string_view> split_lines(std::string_view text) {
  std::vector<std::string_view> lines;
  std::size_t start = 0;
  while (start <= text.size()) {
    std::size_t end = text.find('\n', start);
    if (end == std::string_view::npos) end = text.size();
    std::string_view l = text.substr(start, end - start);
    if (!l.empty() && l.back() == '\r') l.remove_suffix(1);
    lines.push_back(l);
    if (end == text.size()) break;
    start = end + 1;
  }
  if (!lines.empty() && lines.back().empty()) lines.pop_back();
  return lines;
}

std::string_view trim(std::string_view s) {
  while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.remove_prefix(1);
  while (!s.empty() && (s.back() == ' ' || s.back() == '\t')) s.remove_suffix(1);
  return s;
}

std::string_view field(std::string_view line, std::size_t start, std::size_t len) {
  if (start >= line.size()) return {};
  return trim(line.substr(start, len));
}

bool to_int(std::string_view s, int& out) {
  s = trim(s);
  if (!s.empty() && s.front() == '+') s.remove_prefix(1);
  if (s.empty()) return false;
  const auto r = std::from_chars(s.data(), s.data() + s.size(), out);
  return r.ec == std::errc() && r.ptr == s.data() + s.size();
}

bool to_double(std::string_view s, double& out) {
  s = trim(s);
  if (s.empty()) return false;
  std::string tmp(s);
  char* end = nullptr;
  out = std::strtod(tmp.c_str(), &end);
  return end == tmp.c_str() + tmp.size();
}

std::vector<std::string_view> words(std::string_view s) {
  std::vector<std::string_view> out;
  std::size_t i = 0;
  while (i < s.size()) {
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
    std::size_t j = i;
    while (j < s.size() && s[j] != ' ' && s[j] != '\t') ++j;
    if (j > i) out.push_back(s.substr(i, j - i));
    i = j;
  }
  return out;
}

bool starts_with(std::string_view s, std::string_view p) { return s.substr(0, p.size()) == p; }

// Element symbol as written in a molfile; D and T are hydrogen isotopes.
bool mdl_element(std::string_view sym, RAtom& a) {
  if (sym == "D") { a.z = 1; a.isotope = 2; return true; }
  if (sym == "T") { a.z = 1; a.isotope = 3; return true; }
  const ElementInfo* e = find_element_by_symbol(sym);
  if (!e) return false;
  a.z = e->z;
  return true;
}

// ---- MDL molfile / SDF -----------------------------------------------------------

struct SdfRecord {
  RawMol mol;
  std::string title;
  std::vector<std::pair<std::string, std::string>> properties;
  std::size_t line = 0;
};

class MolfileReader {
 public:
  MolfileReader(const std::vector<std::string_view>& lines, Diagnostics& diag) : L_(lines), diag_(diag) {}

  // Reads the record starting at line i; leaves i after its "$$$$" (or at the end).
  SdfRecord next(std::size_t& i);

 private:
  void fail(std::size_t line, DiagCode code, const std::string& msg) {
    diag_.error(code, at_line(line + 1) + msg);
    rec_.mol.failed = true;
  }
  void v2000(std::size_t& i, int na, int nb);
  void v3000(std::size_t& i);
  void data_items(std::size_t& i);
  bool bond_type(int type, std::size_t line, RBond& b);

  const std::vector<std::string_view>& L_;
  Diagnostics& diag_;
  SdfRecord rec_;
  bool sgroup_warned_ = false;
};

bool MolfileReader::bond_type(int type, std::size_t line, RBond& b) {
  if (type >= 1 && type <= 4) {
    b.order = static_cast<std::uint8_t>(type);
    return true;
  }
  if (type == 9 || type == 10) {
    diag_.warn(DiagCode::FileContentIgnored, at_line(line + 1) + (type == 9 ? "coordination" : "hydrogen") +
                                                 std::string(" bond left out"));
    return false;
  }
  fail(line, DiagCode::UnsupportedFileFeature, "query bond type " + std::to_string(type) + " is not a structure");
  return false;
}

SdfRecord MolfileReader::next(std::size_t& i) {
  rec_ = SdfRecord{};
  rec_.line = i + 1;
  rec_.title = std::string(trim(L_[i]));
  const std::size_t counts = i + 3;
  if (counts >= L_.size()) {
    fail(i, DiagCode::MalformedRecord, "molfile header is shorter than 4 lines");
    i = L_.size();
    return std::move(rec_);
  }
  const std::string_view c = L_[counts];
  i = counts + 1;
  if (c.find("V3000") != std::string_view::npos) {
    v3000(i);
  } else {
    int na = 0, nb = 0;
    if (!to_int(field(c, 0, 3), na) || !to_int(field(c, 3, 3), nb)) {
      const auto w = words(c);
      if (w.size() < 2 || !to_int(w[0], na) || !to_int(w[1], nb)) {
        fail(counts, DiagCode::MalformedRecord, "counts line does not start with the atom and bond counts");
        na = nb = -1;
      }
    }
    if (na >= 0) v2000(i, na, nb);
  }
  data_items(i);
  return std::move(rec_);
}

void MolfileReader::v2000(std::size_t& i, int na, int nb) {
  RawMol& m = rec_.mol;
  if (i + static_cast<std::size_t>(na) + static_cast<std::size_t>(nb) > L_.size()) {
    fail(i - 1, DiagCode::MalformedRecord, "the file ends inside the atom or bond block");
    i = L_.size();
    return;
  }
  for (int k = 0; k < na; ++k, ++i) {
    const std::string_view l = L_[i];
    RAtom a;
    a.line = i + 1;
    std::string_view sym = field(l, 31, 3);
    bool ok = to_double(field(l, 0, 10), a.x) && to_double(field(l, 10, 10), a.y) && to_double(field(l, 20, 10), a.z3);
    std::vector<std::string_view> w;
    if (!ok) {  // not column-aligned: fall back to whitespace fields
      w = words(l);
      ok = w.size() >= 4 && to_double(w[0], a.x) && to_double(w[1], a.y) && to_double(w[2], a.z3);
      if (ok) sym = w[3];
    }
    if (!ok) {
      fail(i, DiagCode::MalformedRecord, "atom line needs x, y, z and an element symbol");
      m.atoms.push_back(a);
      continue;
    }
    if (!mdl_element(sym, a)) {
      fail(i, DiagCode::UnsupportedFileFeature, "atom symbol '" + std::string(sym) + "' is not an element (query atoms are not read)");
    }
    auto num = [&](std::size_t col, std::size_t wi) {
      int v = 0;
      if (w.empty()) to_int(field(l, col, 3), v);
      else if (w.size() > wi) to_int(w[wi], v);
      return v;
    };
    int massdiff = 0;
    if (w.empty()) to_int(field(l, 34, 2), massdiff);
    else if (w.size() > 4) to_int(w[4], massdiff);
    const int ccode = num(36, 5);
    static constexpr int kCharge[] = {0, 3, 2, 1, 0, -1, -2, -3};
    if (ccode >= 1 && ccode <= 7) a.charge = kCharge[ccode];
    if (ccode == 4) a.radical = 2;
    const int vcode = num(48, 9);
    if (vcode == 15) a.valence = 0;
    else if (vcode > 0) a.valence = vcode;
    if (massdiff && a.z) {
      if (const ElementInfo* e = find_element_by_z(a.z))
        a.isotope = static_cast<std::uint16_t>(std::lround(e->mass) + massdiff);
    }
    m.atoms.push_back(a);
  }
  for (int k = 0; k < nb; ++k, ++i) {
    const std::string_view l = L_[i];
    int a = 0, b = 0, type = 0, stereo = 0;
    bool ok = to_int(field(l, 0, 3), a) && to_int(field(l, 3, 3), b) && to_int(field(l, 6, 3), type);
    if (ok) {
      to_int(field(l, 9, 3), stereo);
    } else {
      const auto w = words(l);
      ok = w.size() >= 3 && to_int(w[0], a) && to_int(w[1], b) && to_int(w[2], type);
      if (ok && w.size() > 3) to_int(w[3], stereo);
    }
    if (!ok || a < 1 || b < 1 || a > na || b > na || a == b) {
      fail(i, DiagCode::MalformedRecord, "bond line needs two different atom numbers between 1 and " + std::to_string(na) + " and a type");
      continue;
    }
    RBond bd;
    bd.a = static_cast<std::uint32_t>(a - 1);
    bd.b = static_cast<std::uint32_t>(b - 1);
    bd.line = i + 1;
    if (!bond_type(type, i, bd)) continue;
    if (bd.order == 1 && (stereo == 1 || stereo == 4 || stereo == 6)) bd.stereo = static_cast<std::uint8_t>(stereo);
    if (bd.order == 2 && stereo == 3) bd.stereo = kStereoEitherDb;
    m.bonds.push_back(bd);
  }
  // Properties block.
  bool chg_reset = false, iso_reset = false;
  for (; i < L_.size(); ++i) {
    const std::string_view l = L_[i];
    if (starts_with(l, "M  END") || starts_with(l, "$$$$") || starts_with(l, ">")) {
      if (starts_with(l, "M  END")) ++i;
      return;
    }
    if (starts_with(l, "A  ") || starts_with(l, "G  ")) { ++i; continue; }
    if (starts_with(l, "S  SKP")) {
      int n = 0;
      to_int(field(l, 6, 3), n);
      i += static_cast<std::size_t>(std::max(0, n));
      continue;
    }
    const std::string_view tag = l.substr(0, std::min<std::size_t>(6, l.size()));
    if (tag == "M  CHG" || tag == "M  RAD" || tag == "M  ISO") {
      const auto w = words(l.substr(6));
      int n = 0;
      if (w.empty() || !to_int(w[0], n) || w.size() < 1 + 2 * static_cast<std::size_t>(n)) {
        fail(i, DiagCode::MalformedRecord, std::string(tag) + " line has fewer entries than it announces");
        continue;
      }
      if (tag != "M  ISO" && !chg_reset) {  // first CHG/RAD replaces the atom block's charges and radicals
        for (RAtom& a : m.atoms) { a.charge = 0; a.radical = 0; }
        chg_reset = true;
      }
      if (tag == "M  ISO" && !iso_reset) {
        for (RAtom& a : m.atoms) if (a.z != 1 || a.isotope > 3) a.isotope = 0;
        iso_reset = true;
      }
      for (int k = 0; k < n; ++k) {
        int at = 0, v = 0;
        if (!to_int(w[1 + 2 * static_cast<std::size_t>(k)], at) || !to_int(w[2 + 2 * static_cast<std::size_t>(k)], v) ||
            at < 1 || at > static_cast<int>(m.atoms.size())) {
          fail(i, DiagCode::MalformedRecord, std::string(tag) + " refers to an atom that does not exist");
          continue;
        }
        RAtom& a = m.atoms[static_cast<std::size_t>(at - 1)];
        if (tag == "M  CHG") a.charge = v;
        else if (tag == "M  RAD") a.radical = v;
        else a.isotope = static_cast<std::uint16_t>(std::max(0, v));
      }
    } else if (tag == "M  ALS" || tag == "M  RGP" || tag == "M  LIN" || tag == "M  SUB" || tag == "M  UNS" ||
               tag == "M  RBC") {
      fail(i, DiagCode::UnsupportedFileFeature, "query property '" + std::string(trim(tag)) + "' is not a structure");
    } else if (starts_with(l, "M  S") && !sgroup_warned_) {
      diag_.warn(DiagCode::FileContentIgnored, at_line(i + 1) + "S-groups (abbreviations, polymers, data) are not read");
      sgroup_warned_ = true;
    }
  }
}

// Splits a V3000 line into fields, keeping "(...)", "[...]" and quoted strings whole.
std::vector<std::string> v30_fields(std::string_view s) {
  std::vector<std::string> out;
  std::size_t i = 0;
  while (i < s.size()) {
    while (i < s.size() && s[i] == ' ') ++i;
    if (i >= s.size()) break;
    std::string f;
    int depth = 0;
    bool quote = false;
    while (i < s.size() && (depth > 0 || quote || s[i] != ' ')) {
      const char ch = s[i];
      if (ch == '"') quote = !quote;
      else if (!quote && (ch == '(' || ch == '[')) ++depth;
      else if (!quote && (ch == ')' || ch == ']')) --depth;
      f += ch;
      ++i;
    }
    out.push_back(std::move(f));
  }
  return out;
}

void MolfileReader::v3000(std::size_t& i) {
  RawMol& m = rec_.mol;
  // Join continuation lines ('-' at the end) into logical lines.
  std::vector<std::pair<std::string, std::size_t>> logical;
  for (; i < L_.size(); ++i) {
    const std::string_view l = L_[i];
    if (starts_with(l, "M  END")) { ++i; break; }
    if (starts_with(l, "$$$$")) break;
    if (!starts_with(l, "M  V30")) continue;
    std::string body(trim(l.substr(6)));
    const std::size_t first = i;
    while (!body.empty() && body.back() == '-' && i + 1 < L_.size() && starts_with(L_[i + 1], "M  V30")) {
      body.pop_back();
      ++i;
      body += std::string(L_[i].substr(std::min<std::size_t>(7, L_[i].size())));
    }
    logical.push_back({body, first});
  }
  std::unordered_map<long, std::uint32_t> atom_index;
  std::string block;
  bool ignored_warned = false;
  for (const auto& [body, line] : logical) {
    const auto f = v30_fields(body);
    if (f.empty()) continue;
    if (f[0] == "BEGIN" && f.size() > 1) {
      block = f[1];
      if (block != "CTAB" && block != "ATOM" && block != "BOND" && !ignored_warned) {
        diag_.warn(DiagCode::FileContentIgnored, at_line(line + 1) + block + " block is not read");
        ignored_warned = true;
      }
      if (block == "RGROUP" || block == "TEMPLATE") fail(line, DiagCode::UnsupportedFileFeature, block + " blocks are not read");
      continue;
    }
    if (f[0] == "END") { block = "CTAB"; continue; }
    auto kv = [&](const std::string& key, int& out) {
      for (std::size_t k = 0; k < f.size(); ++k)
        if (f[k].size() > key.size() && f[k].compare(0, key.size() + 1, key + "=") == 0)
          return to_int(std::string_view(f[k]).substr(key.size() + 1), out);
      return false;
    };
    if (block == "ATOM") {
      RAtom a;
      a.line = line + 1;
      int idx = 0;
      if (f.size() < 5 || !to_int(f[0], idx) || !to_double(f[2], a.x) || !to_double(f[3], a.y) || !to_double(f[4], a.z3)) {
        fail(line, DiagCode::MalformedRecord, "V3000 atom needs index, type, x, y, z");
        continue;
      }
      std::string sym = f[1];
      if (sym.size() > 1 && sym.front() == '"') sym = sym.substr(1, sym.size() - 2);
      if (!mdl_element(sym, a))
        fail(line, DiagCode::UnsupportedFileFeature, "atom type '" + sym + "' is not an element (query atoms are not read)");
      int v = 0;
      if (kv("CHG", v)) a.charge = v;
      if (kv("RAD", v)) a.radical = v;
      if (kv("MASS", v)) a.isotope = static_cast<std::uint16_t>(std::max(0, v));
      if (kv("VAL", v)) a.valence = v < 0 ? 0 : v;
      atom_index[idx] = static_cast<std::uint32_t>(m.atoms.size());
      m.atoms.push_back(a);
    } else if (block == "BOND") {
      int type = 0, a1 = 0, a2 = 0, idx = 0;
      if (f.size() < 4 || !to_int(f[0], idx) || !to_int(f[1], type) || !to_int(f[2], a1) || !to_int(f[3], a2) ||
          !atom_index.count(a1) || !atom_index.count(a2) || a1 == a2) {
        fail(line, DiagCode::MalformedRecord, "V3000 bond needs index, type and two different existing atoms");
        continue;
      }
      RBond b;
      b.a = atom_index[a1];
      b.b = atom_index[a2];
      b.line = line + 1;
      if (!bond_type(type, line, b)) continue;
      int cfg = 0;
      if (kv("CFG", cfg)) {
        if (b.order == 1 && cfg == 1) b.stereo = kStereoWedge;
        if (b.order == 1 && cfg == 3) b.stereo = kStereoHash;
        if (b.order == 1 && cfg == 2) b.stereo = kStereoEither;
        if (b.order == 2 && cfg == 2) b.stereo = kStereoEitherDb;
      }
      m.bonds.push_back(b);
    }
  }
}

void MolfileReader::data_items(std::size_t& i) {
  for (; i < L_.size(); ++i) {
    const std::string_view l = L_[i];
    if (starts_with(l, "$$$$")) { ++i; return; }
    if (!starts_with(l, ">")) continue;
    std::string name;
    const std::size_t lt = l.find('<'), gt = l.rfind('>');
    if (lt != std::string_view::npos && gt != std::string_view::npos && gt > lt) name = std::string(l.substr(lt + 1, gt - lt - 1));
    std::string value;
    while (i + 1 < L_.size() && !trim(L_[i + 1]).empty() && !starts_with(L_[i + 1], "$$$$")) {
      ++i;
      if (!value.empty()) value += '\n';
      value += std::string(L_[i]);
    }
    rec_.properties.push_back({name, value});
  }
}

// ---- PDB ---------------------------------------------------------------------------

// Single-bond covalent radii, A (Cordero et al., Dalton Trans. 2008, 2832).
double covalent_radius(std::uint8_t z) {
  switch (z) {
    case 1: return 0.31;
    case 5: return 0.84;
    case 6: return 0.76;
    case 7: return 0.71;
    case 8: return 0.66;
    case 9: return 0.57;
    case 14: return 1.11;
    case 15: return 1.07;
    case 16: return 1.05;
    case 17: return 1.02;
    case 33: return 1.19;
    case 34: return 1.20;
    case 35: return 1.20;
    case 53: return 1.39;
    default: return 0.0;  // metals and the rest: no distance bonds
  }
}

// Heavy-atom bonds of the standard amino acids by PDB atom name (neutral
// forms, Kekulé rings); '=' marks a double bond. The backbone N-CA, CA-C, C=O,
// C-OXT and CA-CB are added to every residue except where noted.
struct ResidueTemplate {
  bool has_cb = true;
  std::vector<std::string> bonds;
};

const std::map<std::string, ResidueTemplate>& residue_templates() {
  static const std::map<std::string, ResidueTemplate> t = {
      {"ALA", {true, {}}},
      {"ARG", {true, {"CB-CG", "CG-CD", "CD-NE", "NE-CZ", "CZ-NH1", "CZ=NH2"}}},
      {"ASN", {true, {"CB-CG", "CG=OD1", "CG-ND2"}}},
      {"ASP", {true, {"CB-CG", "CG=OD1", "CG-OD2"}}},
      {"CYS", {true, {"CB-SG"}}},
      {"GLN", {true, {"CB-CG", "CG-CD", "CD=OE1", "CD-NE2"}}},
      {"GLU", {true, {"CB-CG", "CG-CD", "CD=OE1", "CD-OE2"}}},
      {"GLY", {false, {}}},
      {"HIS", {true, {"CB-CG", "CG-ND1", "CG=CD2", "ND1=CE1", "CE1-NE2", "NE2-CD2"}}},
      {"ILE", {true, {"CB-CG1", "CB-CG2", "CG1-CD1"}}},
      {"LEU", {true, {"CB-CG", "CG-CD1", "CG-CD2"}}},
      {"LYS", {true, {"CB-CG", "CG-CD", "CD-CE", "CE-NZ"}}},
      {"MET", {true, {"CB-CG", "CG-SD", "SD-CE"}}},
      {"MSE", {true, {"CB-CG", "CG-SE", "SE-CE"}}},
      {"PHE", {true, {"CB-CG", "CG=CD1", "CG-CD2", "CD1-CE1", "CD2=CE2", "CE1=CZ", "CE2-CZ"}}},
      {"PRO", {true, {"CB-CG", "CG-CD", "CD-N"}}},
      {"SER", {true, {"CB-OG"}}},
      {"THR", {true, {"CB-OG1", "CB-CG2"}}},
      {"TRP", {true, {"CB-CG", "CG=CD1", "CG-CD2", "CD1-NE1", "NE1-CE2", "CE2-CD2", "CE2=CZ2", "CZ2-CH2", "CH2=CZ3",
                      "CZ3-CE3", "CE3=CD2"}}},
      {"TYR", {true, {"CB-CG", "CG=CD1", "CG-CD2", "CD1-CE1", "CD2=CE2", "CE1=CZ", "CE2-CZ", "CZ-OH"}}},
      {"VAL", {true, {"CB-CG1", "CB-CG2"}}},
  };
  return t;
}

std::string canonical_residue(const std::string& r) {
  if (r == "HID" || r == "HIE" || r == "HIP" || r == "HSD" || r == "HSE" || r == "HSP") return "HIS";
  if (r == "CYX" || r == "CYM") return "CYS";
  if (r == "ASH") return "ASP";
  if (r == "GLH") return "GLU";
  if (r == "LYN") return "LYS";
  return r;
}

struct PdbAtom {
  std::string name, resname, key;  // key: chain + resSeq + iCode
  int serial = 0;
  bool hetatm = false;
};

RawMol read_pdb(const std::vector<std::string_view>& L, std::string& title, Diagnostics& diag) {
  RawMol m;
  std::vector<PdbAtom> info;
  std::unordered_map<int, std::uint32_t> by_serial;
  std::unordered_map<std::string, char> altloc_of;  // residue key + atom name -> chosen altLoc
  std::map<std::pair<std::uint32_t, std::uint32_t>, int> conect;
  std::vector<std::size_t> links;  // LINK lines, resolved once all atoms are read
  bool more_models = false, model_done = false, altloc_warned = false;
  int models = 0;
  for (std::size_t i = 0; i < L.size(); ++i) {
    const std::string_view l = L[i];
    const std::string_view rec = l.substr(0, std::min<std::size_t>(6, l.size()));
    if (rec == "MODEL ") {
      if (++models > 1) more_models = true;
      continue;
    }
    if (rec == "ENDMDL") { model_done = true; continue; }
    if (trim(l) == "END") break;
    if (title.empty() && (rec == "HEADER" || rec == "COMPND" || rec == "TITLE ")) {
      title = std::string(trim(l.substr(std::min<std::size_t>(10, l.size()))));
      continue;
    }
    if (rec == "CONECT") {
      int a = 0;
      if (!to_int(field(l, 6, 5), a)) continue;
      for (std::size_t c = 11; c < 31; c += 5) {
        int b = 0;
        if (!to_int(field(l, c, 5), b)) continue;
        if (!by_serial.count(a) || !by_serial.count(b)) continue;
        const std::uint32_t x = by_serial[a], y = by_serial[b];
        if (x == y) continue;
        ++conect[{std::min(x, y), std::max(x, y)} ];
      }
      continue;
    }
    if (rec == "LINK  ") {
      links.push_back(i);
      continue;
    }
    if (rec != "ATOM  " && rec != "HETATM") continue;
    if (model_done) continue;
    RAtom a;
    a.line = i + 1;
    if (!to_double(field(l, 30, 8), a.x) || !to_double(field(l, 38, 8), a.y) || !to_double(field(l, 46, 8), a.z3)) {
      diag.error(DiagCode::MalformedRecord, at_line(i + 1) + "atom record needs x, y, z in columns 31-54");
      m.failed = true;
      continue;
    }
    PdbAtom p;
    p.hetatm = rec == "HETATM";
    to_int(field(l, 6, 5), p.serial);
    const std::string raw_name = l.size() > 12 ? std::string(l.substr(12, std::min<std::size_t>(4, l.size() - 12))) : "";
    p.name = std::string(trim(raw_name));
    p.resname = std::string(field(l, 17, 3));
    p.key = std::string(l.size() > 21 ? l.substr(21, 1) : " ") + std::string(field(l, 22, 4)) +
            std::string(l.size() > 26 ? l.substr(26, 1) : " ");
    const char alt = l.size() > 16 ? l[16] : ' ';
    if (alt != ' ') {
      auto [it, inserted] = altloc_of.insert({p.key + "|" + p.name, alt});
      if (!inserted && it->second != alt) {
        if (!altloc_warned)
          diag.warn(DiagCode::FileContentIgnored, at_line(i + 1) + "alternate locations after the first are not read");
        altloc_warned = true;
        continue;
      }
    }
    std::string sym(field(l, 76, 2));
    if (sym.empty()) {  // derive from the atom name: right-justified one-letter elements start in column 14
      std::string nm = raw_name;
      while (nm.size() < 4) nm += ' ';
      if (nm[0] == ' ' || (nm[0] >= '0' && nm[0] <= '9')) sym = std::string(1, nm[1]);
      else sym = p.hetatm ? nm.substr(0, 2) : std::string(1, nm[0]);
      if (!find_element_by_symbol(sym) && sym.size() == 2) sym = sym.substr(0, 1);
    }
    if (sym.size() == 2) sym[1] = static_cast<char>(std::tolower(static_cast<unsigned char>(sym[1])));
    if (!mdl_element(sym, a)) {
      diag.error(DiagCode::UnknownElement, at_line(i + 1) + "unknown element '" + sym + "'");
      m.failed = true;
      continue;
    }
    const std::string_view chg = field(l, 78, 2);
    if (chg.size() == 2 && chg[0] >= '0' && chg[0] <= '9') a.charge = (chg[1] == '-' ? -1 : 1) * (chg[0] - '0');
    by_serial[p.serial] = static_cast<std::uint32_t>(m.atoms.size());
    m.atoms.push_back(a);
    info.push_back(p);
  }
  if (more_models) diag.warn(DiagCode::FileContentIgnored, "only the first MODEL is read");

  const std::size_t n = m.atoms.size();
  bool any_h = false;
  for (const RAtom& a : m.atoms) any_h = any_h || a.z == 1;
  m.implicit_h = !any_h;

  // Bonds: CONECT first (a repeated pair is a double bond), then distances.
  std::vector<char> in_conect(n, 0);
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint8_t> bonds;
  {
    std::unordered_map<std::string, std::uint32_t> by_name;  // residue key + "|" + atom name
    for (std::uint32_t k = 0; k < n; ++k) by_name.insert({info[k].key + "|" + info[k].name, k});
    for (std::size_t li : links) {
      const std::string_view l = L[li];
      auto atom_at = [&](std::size_t name_col, std::size_t chain_col) -> std::int64_t {
        const std::string key = std::string(l.size() > chain_col ? l.substr(chain_col, 1) : " ") +
                                std::string(field(l, chain_col + 1, 4)) +
                                std::string(l.size() > chain_col + 5 ? l.substr(chain_col + 5, 1) : " ");
        auto it = by_name.find(key + "|" + std::string(field(l, name_col, 4)));
        return it == by_name.end() ? -1 : static_cast<std::int64_t>(it->second);
      };
      const std::int64_t x = atom_at(12, 21), y = atom_at(42, 51);
      if (x < 0 || y < 0 || x == y) continue;
      const auto a = static_cast<std::uint32_t>(std::min(x, y)), b = static_cast<std::uint32_t>(std::max(x, y));
      if (!conect.count({a, b})) conect[{a, b}] = 2;  // as if listed from both ends: a single bond
    }
  }
  bool metal_warned = false;
  for (const auto& [pair, count] : conect) {
    // Each bond is normally listed from both ends, so halve symmetric counts.
    const auto x = pair.first, y = pair.second;
    in_conect[x] = in_conect[y] = 1;
    if (covalent_radius(m.atoms[x].z) == 0 || covalent_radius(m.atoms[y].z) == 0) {
      if (!metal_warned)
        diag.warn(DiagCode::FileContentIgnored, at_line(m.atoms[x].line) + "CONECT bonds to metals are left out");
      metal_warned = true;
      continue;
    }
    bonds[pair] = static_cast<std::uint8_t>(std::clamp(count >= 2 ? (count + 1) / 2 : 1, 1, 3));
  }
  {
    const double cell = 2.6;
    std::unordered_map<long long, std::vector<std::uint32_t>> grid;
    auto key = [&](long long gx, long long gy, long long gz) { return (gx * 73856093LL) ^ (gy * 19349663LL) ^ (gz * 83492791LL); };
    auto gc = [&](double v) { return static_cast<long long>(std::floor(v / cell)); };
    // Distance bonds inside a residue, plus the standard links between residues;
    // other links between residues come from CONECT and LINK only.
    auto linkable = [&](std::uint32_t x, std::uint32_t y) {
      const PdbAtom& p = info[x];
      const PdbAtom& q = info[y];
      if (p.hetatm && q.hetatm && in_conect[x] && in_conect[y]) return false;  // CONECT is authoritative
      if (p.key == q.key) {
        // Heavy-atom bonds of standard residues come from the template, not from distances.
        const bool heavy = m.atoms[x].z != 1 && m.atoms[y].z != 1;
        return !(heavy && residue_templates().count(canonical_residue(p.resname)));
      }
      auto pair = [&](const char* u, const char* v) { return (p.name == u && q.name == v) || (p.name == v && q.name == u); };
      return pair("C", "N") || pair("O3'", "P") || pair("O3*", "P") || pair("SG", "SG");
    };
    for (std::uint32_t k = 0; k < n; ++k)
      if (covalent_radius(m.atoms[k].z) > 0)
        grid[key(gc(m.atoms[k].x), gc(m.atoms[k].y), gc(m.atoms[k].z3))].push_back(k);
    for (std::uint32_t k = 0; k < n; ++k) {
      if (covalent_radius(m.atoms[k].z) == 0) continue;
      const RAtom& a = m.atoms[k];
      for (long long dx = -1; dx <= 1; ++dx)
        for (long long dy = -1; dy <= 1; ++dy)
          for (long long dz = -1; dz <= 1; ++dz) {
            auto it = grid.find(key(gc(a.x) + dx, gc(a.y) + dy, gc(a.z3) + dz));
            if (it == grid.end()) continue;
            for (std::uint32_t j : it->second) {
              if (j <= k) continue;
              const RAtom& b = m.atoms[j];
              if ((a.z == 1 && b.z == 1) || bonds.count({k, j}) || !linkable(k, j)) continue;
              const double d = std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z3 - b.z3) * (a.z3 - b.z3));
              if (d > 0.4 && d < covalent_radius(a.z) + covalent_radius(b.z) + 0.45) bonds[{k, j}] = 1;
            }
          }
    }
  }
  // A hydrogen keeps only its shortest bond.
  {
    std::vector<std::pair<double, std::pair<std::uint32_t, std::uint32_t>>> best(n, {1e9, {0, 0}});
    for (const auto& [pair, order] : bonds) {
      for (const std::uint32_t h : {pair.first, pair.second}) {
        if (m.atoms[h].z != 1) continue;
        const RAtom& a = m.atoms[pair.first];
        const RAtom& b = m.atoms[pair.second];
        const double d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z3 - b.z3) * (a.z3 - b.z3);
        if (d < best[h].first) best[h] = {d, pair};
      }
    }
    for (auto it = bonds.begin(); it != bonds.end();) {
      bool drop = false;
      for (const std::uint32_t h : {it->first.first, it->first.second})
        if (m.atoms[h].z == 1 && best[h].second != it->first) drop = true;
      it = drop ? bonds.erase(it) : std::next(it);
    }
  }
  // Residue templates.
  std::map<std::string, std::unordered_map<std::string, std::uint32_t>> residues;  // key -> name -> atom
  std::map<std::string, std::string> resname_of;
  for (std::uint32_t k = 0; k < n; ++k) {
    residues[info[k].key][info[k].name] = k;
    resname_of[info[k].key] = info[k].resname;
  }
  std::vector<std::string> guessed;
  const auto& templates = residue_templates();
  for (const auto& [key, atoms] : residues) {
    const std::string res = canonical_residue(resname_of[key]);
    auto t = templates.find(res);
    auto add = [&](const std::string& x, const std::string& y, std::uint8_t order) {
      auto ix = atoms.find(x), iy = atoms.find(y);
      if (ix == atoms.end() || iy == atoms.end()) return;
      bonds[{std::min(ix->second, iy->second), std::max(ix->second, iy->second)}] = order;
    };
    if (t != templates.end()) {
      add("N", "CA", 1);
      add("CA", "C", 1);
      add("C", "O", 2);
      add("C", "OXT", 1);
      if (t->second.has_cb) add("CA", "CB", 1);
      for (const std::string& b : t->second.bonds) {
        const std::size_t sep = b.find_first_of("-=");
        add(b.substr(0, sep), b.substr(sep + 1), b[sep] == '=' ? 2 : 1);
      }
      continue;
    }
    std::size_t heavy = 0;
    for (const auto& [name, k] : atoms) heavy += m.atoms[k].z != 1;
    bool has_conect_order = false;
    for (const auto& [name, k] : atoms)
      for (const auto& [pair, count] : conect)
        if ((pair.first == k || pair.second == k) && count > 2) has_conect_order = true;
    if (heavy > 1 && res != "HOH" && res != "WAT" && !has_conect_order &&
        std::find(guessed.begin(), guessed.end(), res) == guessed.end())
      guessed.push_back(res);
  }
  if (!guessed.empty()) {
    std::string list;
    for (const auto& g : guessed) list += (list.empty() ? "" : ", ") + g;
    diag.warn(DiagCode::BondOrdersGuessed, "PDB gives no bond orders for residue(s) " + list +
                                               "; their bonds are read as single (supply CONECT orders or a molfile)");
  }
  for (const auto& [pair, order] : bonds) {
    RBond b;
    b.a = pair.first;
    b.b = pair.second;
    b.order = order;
    m.bonds.push_back(b);
  }
  // With explicit hydrogens, the protonation is given: charge 4-valent N and
  // lone carboxylate / phenolate-like O in standard residues accordingly.
  if (any_h) {
    std::vector<int> valence(n, 0);
    for (const RBond& b : m.bonds) { valence[b.a] += b.order; valence[b.b] += b.order; }
    for (std::uint32_t k = 0; k < n; ++k) {
      if (m.atoms[k].charge || !templates.count(canonical_residue(info[k].resname))) continue;
      if (m.atoms[k].z == 7 && valence[k] == 4) m.atoms[k].charge = 1;
      if (m.atoms[k].z == 8 && valence[k] == 1) m.atoms[k].charge = -1;
    }
  }
  return m;
}

// ---- CML ---------------------------------------------------------------------------

std::string xml_unescape(std::string_view s) {
  std::string out;
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (s[i] != '&') { out += s[i]; continue; }
    const std::size_t semi = s.find(';', i);
    if (semi == std::string_view::npos) { out += s[i]; continue; }
    const std::string_view ent = s.substr(i + 1, semi - i - 1);
    if (ent == "amp") out += '&';
    else if (ent == "lt") out += '<';
    else if (ent == "gt") out += '>';
    else if (ent == "quot") out += '"';
    else if (ent == "apos") out += '\'';
    else if (!ent.empty() && ent[0] == '#') {
      const long code = ent.size() > 1 && (ent[1] == 'x' || ent[1] == 'X') ? std::strtol(std::string(ent.substr(2)).c_str(), nullptr, 16)
                                                                          : std::strtol(std::string(ent.substr(1)).c_str(), nullptr, 10);
      if (code > 0 && code < 128) out += static_cast<char>(code);
      else out += '?';
    } else {
      out += std::string(s.substr(i, semi - i + 1));
    }
    i = semi;
  }
  return out;
}

struct CmlMolecule {
  RawMol mol;
  std::string title, id;
  std::size_t line = 0;
  std::unordered_map<std::string, std::uint32_t> ids;
};

std::vector<CmlMolecule> read_cml(std::string_view text, Diagnostics& diag) {
  std::vector<CmlMolecule> done;
  std::vector<CmlMolecule> open;  // stack of <molecule> elements
  std::vector<std::string> path;  // element names
  std::size_t line = 1, pos = 0;
  std::int64_t stereo_bond = -1;
  std::string text_buf;
  auto advance = [&](std::size_t to) {
    for (std::size_t k = pos; k < to && k < text.size(); ++k) line += text[k] == '\n';
    pos = to;
  };
  auto local = [](std::string_view name) {
    const std::size_t colon = name.find(':');
    return std::string(colon == std::string_view::npos ? name : name.substr(colon + 1));
  };
  auto add_atom = [&](CmlMolecule& cm, const std::map<std::string, std::string>& at, std::size_t ln) {
    RAtom a;
    a.line = ln;
    auto get = [&](const char* k) -> const std::string* {
      auto it = at.find(k);
      return it == at.end() ? nullptr : &it->second;
    };
    const std::string* el = get("elementType");
    if (!el || !mdl_element(*el, a)) {
      diag.error(DiagCode::UnsupportedFileFeature, at_line(ln) + "atom has no element or an unknown elementType '" +
                                                       (el ? *el : std::string()) + "'");
      cm.mol.failed = true;
    }
    double v = 0;
    if (const auto* s = get("x3"); s && to_double(*s, v)) {
      a.x = v;
      if (const auto* y = get("y3")) to_double(*y, a.y);
      if (const auto* z = get("z3")) to_double(*z, a.z3);
    } else if (const auto* s2 = get("x2"); s2 && to_double(*s2, v)) {
      a.x = v;
      if (const auto* y = get("y2")) to_double(*y, a.y);
    }
    int iv = 0;
    if (const auto* s = get("formalCharge"); s && to_int(*s, iv)) a.charge = iv;
    if (const auto* s = get("hydrogenCount"); s && to_int(*s, iv)) a.hcount = iv;
    if (const auto* s = get("isotopeNumber"); s && to_int(*s, iv)) a.isotope = static_cast<std::uint16_t>(iv);
    else if (const auto* s3 = get("isotope"); s3 && to_int(*s3, iv)) a.isotope = static_cast<std::uint16_t>(iv);
    if (const auto* s = get("spinMultiplicity"); s && to_int(*s, iv)) a.radical = iv == 2 ? 2 : iv == 3 ? 3 : 0;
    const std::string id = get("id") ? *get("id") : "a" + std::to_string(cm.mol.atoms.size() + 1);
    cm.ids[id] = static_cast<std::uint32_t>(cm.mol.atoms.size());
    cm.mol.atoms.push_back(a);
  };
  auto add_bond = [&](CmlMolecule& cm, const std::string& r1, const std::string& r2, const std::string& order, std::size_t ln) {
    auto a = cm.ids.find(r1), b = cm.ids.find(r2);
    if (a == cm.ids.end() || b == cm.ids.end() || a->second == b->second) {
      diag.error(DiagCode::MalformedRecord, at_line(ln) + "bond refers to atoms '" + r1 + "', '" + r2 + "' that are not defined");
      cm.mol.failed = true;
      return;
    }
    RBond bd;
    bd.a = a->second;
    bd.b = b->second;
    bd.line = ln;
    if (order == "1" || order == "S" || order.empty()) bd.order = 1;
    else if (order == "2" || order == "D") bd.order = 2;
    else if (order == "3" || order == "T") bd.order = 3;
    else if (order == "A" || order == "4") bd.order = 4;
    else {
      diag.error(DiagCode::UnsupportedFileFeature, at_line(ln) + "bond order '" + order + "' is not read");
      cm.mol.failed = true;
      return;
    }
    cm.mol.bonds.push_back(bd);
  };

  while (pos < text.size()) {
    const std::size_t lt = text.find('<', pos);
    if (lt == std::string_view::npos) break;
    text_buf += std::string(text.substr(pos, lt - pos));
    advance(lt);
    if (text.substr(lt, 4) == "<!--") {
      const std::size_t e = text.find("-->", lt);
      advance(e == std::string_view::npos ? text.size() : e + 3);
      continue;
    }
    if (text.substr(lt, 9) == "<![CDATA[") {
      const std::size_t e = text.find("]]>", lt);
      text_buf += std::string(text.substr(lt + 9, (e == std::string_view::npos ? text.size() : e) - lt - 9));
      advance(e == std::string_view::npos ? text.size() : e + 3);
      continue;
    }
    if (text.substr(lt, 2) == "<?" || text.substr(lt, 2) == "<!") {
      const std::size_t e = text.find('>', lt);
      advance(e == std::string_view::npos ? text.size() : e + 1);
      continue;
    }
    // Find the end of the tag, skipping '>' inside quoted attribute values.
    std::size_t e = lt + 1;
    char q = 0;
    for (; e < text.size(); ++e) {
      if (q) { if (text[e] == q) q = 0; }
      else if (text[e] == '"' || text[e] == '\'') q = text[e];
      else if (text[e] == '>') break;
    }
    if (e >= text.size()) {
      diag.error(DiagCode::MalformedRecord, at_line(line) + "unterminated XML tag");
      break;
    }
    std::string_view tag = text.substr(lt + 1, e - lt - 1);
    const std::size_t tag_line = line;
    advance(e + 1);
    const bool closing = !tag.empty() && tag[0] == '/';
    const bool self_closing = !tag.empty() && tag.back() == '/';
    if (closing) tag.remove_prefix(1);
    if (self_closing) tag.remove_suffix(1);
    std::size_t k = 0;
    while (k < tag.size() && tag[k] != ' ' && tag[k] != '\t' && tag[k] != '\n' && tag[k] != '\r') ++k;
    const std::string name = local(tag.substr(0, k));
    std::map<std::string, std::string> at;
    while (k < tag.size()) {
      while (k < tag.size() && (tag[k] == ' ' || tag[k] == '\t' || tag[k] == '\n' || tag[k] == '\r')) ++k;
      std::size_t eq = tag.find('=', k);
      if (eq == std::string_view::npos) break;
      const std::string key = local(trim(tag.substr(k, eq - k)));
      std::size_t v = eq + 1;
      while (v < tag.size() && tag[v] == ' ') ++v;
      if (v >= tag.size()) break;
      const char quote = tag[v];
      std::size_t ve = tag.find(quote, v + 1);
      if (ve == std::string_view::npos) ve = tag.size();
      at[key] = xml_unescape(tag.substr(v + 1, ve - v - 1));
      k = ve + 1;
    }

    if (closing) {
      const std::string body = xml_unescape(trim(text_buf));
      text_buf.clear();
      if (name == "bondStereo" && stereo_bond >= 0 && !open.empty()) {
        RBond& b = open.back().mol.bonds[static_cast<std::size_t>(stereo_bond)];
        if (b.order == 1 && body == "W") b.stereo = kStereoWedge;
        if (b.order == 1 && body == "H") b.stereo = kStereoHash;
        stereo_bond = -1;
      } else if (name == "name" && !open.empty() && path.size() >= 2 && path[path.size() - 2] == "molecule" &&
                 open.back().title.empty()) {
        open.back().title = body;
      } else if (name == "molecule" && !open.empty()) {
        if (open.back().title.empty()) open.back().title = open.back().id;  // title, then <name>, then id
        if (!open.back().mol.atoms.empty() || open.back().mol.failed) done.push_back(std::move(open.back()));
        open.pop_back();
      }
      if (!path.empty()) path.pop_back();
      continue;
    }
    text_buf.clear();
    if (name == "molecule") {
      CmlMolecule cm;
      cm.line = tag_line;
      if (at.count("title")) cm.title = at["title"];
      if (at.count("id")) cm.id = at["id"];
      open.push_back(std::move(cm));
    } else if (!open.empty() && name == "atom") {
      add_atom(open.back(), at, tag_line);
    } else if (!open.empty() && name == "atomArray" && at.count("atomID")) {
      // Array form: parallel space-separated attribute lists.
      std::map<std::string, std::vector<std::string_view>> cols;
      std::map<std::string, std::string> keep = at;
      for (const auto& [key, value] : keep) cols[key] = words(value);
      const std::size_t na = cols["atomID"].size();
      static const std::map<std::string, std::string> rename = {{"atomID", "id"}};
      for (std::size_t a = 0; a < na; ++a) {
        std::map<std::string, std::string> one;
        for (const auto& [key, list] : cols)
          if (a < list.size()) one[rename.count(key) ? rename.at(key) : key] = std::string(list[a]);
        add_atom(open.back(), one, tag_line);
      }
    } else if (!open.empty() && name == "bond") {
      const auto refs = words(at["atomRefs2"]);
      if (refs.size() != 2) {
        diag.error(DiagCode::MalformedRecord, at_line(tag_line) + "bond needs atomRefs2 with two atom ids");
        open.back().mol.failed = true;
      } else {
        add_bond(open.back(), std::string(refs[0]), std::string(refs[1]), at["order"], tag_line);
        stereo_bond = static_cast<std::int64_t>(open.back().mol.bonds.size()) - 1;
      }
    } else if (!open.empty() && name == "bondArray" && at.count("atomRef1")) {
      const auto r1 = words(at["atomRef1"]), r2 = words(at["atomRef2"]), ord = words(at["order"]);
      for (std::size_t b = 0; b < r1.size() && b < r2.size(); ++b)
        add_bond(open.back(), std::string(r1[b]), std::string(r2[b]), b < ord.size() ? std::string(ord[b]) : "1", tag_line);
    }
    if (!self_closing) path.push_back(name);
    else if (name == "molecule" && !open.empty()) {
      open.pop_back();
    }
  }
  return done;
}

}  // namespace

// ---- public ------------------------------------------------------------------------

int mdl_valence_field(std::uint8_t z, int charge, int bond_order_sum, int hydrogens) noexcept {
  int fill = 0;
  if (organic(z)) {
    for (std::uint8_t v : default_valences(z, charge)) {
      if (v >= bond_order_sum) { fill = v - bond_order_sum; break; }
    }
  }
  // Readers disagree on implicit hydrogens for elements outside B C N O P S F Cl Br I,
  // so those always get an explicit valence.
  if (organic(z) && fill == hydrogens) return 0;
  const int total = bond_order_sum + hydrogens;
  return total == 0 ? 15 : total;
}

const char* to_string(MolFormat f) noexcept {
  switch (f) {
    case MolFormat::Auto: return "auto";
    case MolFormat::Molfile: return "molfile";
    case MolFormat::Pdb: return "pdb";
    case MolFormat::Cml: return "cml";
  }
  return "auto";
}

MolFormat detect_format(std::string_view text, std::string_view filename) {
  std::string ext;
  if (const std::size_t dot = filename.rfind('.'); dot != std::string_view::npos)
    for (char c : filename.substr(dot + 1)) ext += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  if (ext == "mol" || ext == "sdf" || ext == "sd" || ext == "mdl") return MolFormat::Molfile;
  if (ext == "pdb" || ext == "ent") return MolFormat::Pdb;
  if (ext == "cml") return MolFormat::Cml;
  std::string_view t = text;
  while (!t.empty() && (t.front() == ' ' || t.front() == '\n' || t.front() == '\r' || t.front() == '\t' ||
                        static_cast<unsigned char>(t.front()) == 0xEF || static_cast<unsigned char>(t.front()) == 0xBB ||
                        static_cast<unsigned char>(t.front()) == 0xBF))
    t.remove_prefix(1);
  if (!t.empty() && t.front() == '<') return MolFormat::Cml;
  const std::string_view head = text.substr(0, 1 << 16);
  if (head.find("V2000") != std::string_view::npos || head.find("V3000") != std::string_view::npos ||
      head.find("M  END") != std::string_view::npos)
    return MolFormat::Molfile;
  const auto lines = split_lines(head);
  for (const auto& l : lines)
    if (starts_with(l, "ATOM  ") || starts_with(l, "HETATM")) return MolFormat::Pdb;
  return MolFormat::Auto;
}

std::vector<MolRecord> read_molecules(ChemicalSystemState& state, std::string_view text, Diagnostics& diag,
                                      const ReadOptions& options) {
  std::vector<MolRecord> out;
  const MolFormat format = options.format == MolFormat::Auto ? detect_format(text) : options.format;
  auto finish = [&](RawMol& m, MolRecord r) {
    r.dimensions = record_dimensions(m);
    if (!m.failed && m.atoms.empty()) diag.error(DiagCode::MalformedRecord, at_line(r.line) + "record has no atoms");
    else r.species = build_species(state, m, options.keep_hydrogens, diag);
    out.push_back(std::move(r));
    return options.max_records == 0 || out.size() < options.max_records;
  };
  switch (format) {
    case MolFormat::Auto:
      diag.error(DiagCode::UnknownFileFormat, "not a recognised structure file (MOL/SDF, PDB or CML)");
      return out;
    case MolFormat::Molfile: {
      const auto lines = split_lines(text);
      MolfileReader reader(lines, diag);
      std::size_t i = 0;
      while (i < lines.size()) {
        std::size_t j = i;
        while (j < lines.size() && trim(lines[j]).empty()) ++j;
        if (j == lines.size()) break;
        SdfRecord rec = reader.next(i);
        MolRecord r;
        r.title = rec.title;
        r.properties = std::move(rec.properties);
        r.line = rec.line;
        if (!finish(rec.mol, std::move(r))) break;
      }
      break;
    }
    case MolFormat::Pdb: {
      const auto lines = split_lines(text);
      MolRecord r;
      r.line = 1;
      RawMol m = read_pdb(lines, r.title, diag);
      finish(m, std::move(r));
      break;
    }
    case MolFormat::Cml: {
      for (CmlMolecule& cm : read_cml(text, diag)) {
        MolRecord r;
        r.title = cm.title;
        r.line = cm.line;
        if (!finish(cm.mol, std::move(r))) break;
      }
      if (out.empty()) diag.error(DiagCode::MalformedRecord, "no <molecule> with atoms found");
      break;
    }
  }
  return out;
}

}  // namespace omnichem
