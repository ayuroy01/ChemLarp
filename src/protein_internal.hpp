#pragma once
// Internal to the protein sources: residue templates, vector helpers, a cell grid.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <numbers>
#include <string>
#include <unordered_map>
#include <vector>

#include "chemlarp/protein.hpp"

namespace chemlarp::detail {

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

inline const std::map<std::string, ResidueDef>& residue_defs() {
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
inline std::string canonical_name(const std::string& n) {
  static const std::map<std::string, std::string> alias = {
      {"HID", "HIS"}, {"HIE", "HIS"}, {"HIP", "HIS"}, {"HSD", "HIS"}, {"HSE", "HIS"}, {"HSP", "HIS"}, {"HSE", "HIS"},
      {"ASH", "ASP"}, {"GLH", "GLU"}, {"LYN", "LYS"}, {"CYX", "CYS"}, {"CYM", "CYS"}, {"TYM", "TYR"}, {"ARN", "ARG"}};
  auto it = alias.find(n);
  return it != alias.end() ? it->second : n;
}

using Vec = std::array<double, 3>;
inline Vec pos(const ProteinAtom& a) { return {a.x, a.y, a.z}; }
inline Vec sub(Vec a, Vec b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline Vec add(Vec a, Vec b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline Vec mul(Vec a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
inline double dot(Vec a, Vec b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline Vec cross(Vec a, Vec b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
inline double norm(Vec a) { return std::sqrt(dot(a, a)); }
inline Vec unit(Vec a) {
  const double n = norm(a);
  return n > 0 ? mul(a, 1.0 / n) : Vec{1, 0, 0};
}
inline double dist(const ProteinAtom& a, const ProteinAtom& b) { return norm(sub(pos(a), pos(b))); }

// Position of atom D given A-B-C, bond C-D, angle B-C-D and dihedral A-B-C-D (degrees).
inline Vec place(Vec A, Vec B, Vec C, double bond, double angle_deg, double dihedral_deg) {
  const double ang = angle_deg * std::numbers::pi / 180.0, dih = dihedral_deg * std::numbers::pi / 180.0;
  const Vec bc = unit(sub(C, B));
  const Vec n = unit(cross(sub(B, A), bc));
  const Vec m = cross(n, bc);
  const Vec d2 = {-bond * std::cos(ang), bond * std::sin(ang) * std::cos(dih), bond * std::sin(ang) * std::sin(dih)};
  return add(C, add(add(mul(bc, d2[0]), mul(m, d2[1])), mul(n, d2[2])));
}

inline std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(' ');
  if (a == std::string::npos) return "";
  return s.substr(a, s.find_last_not_of(' ') - a + 1);
}

inline std::string residue_label(const ProteinResidue& r) {
  return (r.chain.empty() ? std::string("_") : r.chain) + " " + std::to_string(r.number) +
         (r.icode != ' ' ? std::string(1, r.icode) : "") + " " + r.name;
}

// Model pKa values (PROPKA 3).
struct Pka {
  double asp = 3.80, glu = 4.50, his = 6.50, cys = 9.00, tyr = 10.00, lys = 10.50, arg = 12.50, nterm = 8.00, cterm = 3.20;
};


inline Vec pos3(double x, double y, double z) { return {x, y, z}; }
inline void set_pos(ProteinAtom& a, Vec x) { a.x = x[0], a.y = x[1], a.z = x[2]; }
inline double dihedral_deg(Vec a, Vec b, Vec c, Vec d) {
  const Vec b1 = unit(sub(c, b));
  const Vec b0 = sub(a, b), b2 = sub(d, c);
  const Vec v = sub(b0, mul(b1, dot(b0, b1))), w = sub(b2, mul(b1, dot(b2, b1)));
  return std::atan2(dot(cross(b1, v), w), dot(v, w)) * 180.0 / std::numbers::pi;
}

// Points in cubic cells for neighbour queries; ids are whatever the caller stores.
class Grid {
 public:
  explicit Grid(double cell = 4.0) : cell_(cell) {}
  void add(Vec x, std::uint32_t id) { cells_[key(x)].push_back({x, id}); }
  // Calls f(id, position, distance squared) for every stored point within r of x.
  template <class F>
  void query(Vec x, double r, F&& f) const {
    const int n = int(std::ceil(r / cell_));
    const std::array<int, 3> c = {int(std::floor(x[0] / cell_)), int(std::floor(x[1] / cell_)), int(std::floor(x[2] / cell_))};
    const double r2 = r * r;
    for (int dx = -n; dx <= n; ++dx)
      for (int dy = -n; dy <= n; ++dy)
        for (int dz = -n; dz <= n; ++dz) {
          auto it = cells_.find(pack(c[0] + dx, c[1] + dy, c[2] + dz));
          if (it == cells_.end()) continue;
          for (const auto& e : it->second) {
            const Vec d = sub(e.x, x);
            const double d2 = dot(d, d);
            if (d2 <= r2) f(e.id, e.x, d2);
          }
        }
  }

 private:
  struct Entry {
    Vec x;
    std::uint32_t id;
  };
  static std::uint64_t pack(int a, int b, int c) {
    return (std::uint64_t(std::uint32_t(a) & 0x1FFFFF) << 42) | (std::uint64_t(std::uint32_t(b) & 0x1FFFFF) << 21) |
           std::uint64_t(std::uint32_t(c) & 0x1FFFFF);
  }
  std::uint64_t key(Vec x) const { return pack(int(std::floor(x[0] / cell_)), int(std::floor(x[1] / cell_)), int(std::floor(x[2] / cell_))); }
  double cell_;
  std::unordered_map<std::uint64_t, std::vector<Entry>> cells_;
};

}  // namespace chemlarp::detail
