#include "chemlarp/uff.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "chemlarp/elements.hpp"

namespace chemlarp {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;

// Rappé et al. (1992), Table 1.  label, r1, theta0, x1, D1, zeta, Z1, Vi, Uj, Xi
constexpr UffParameters kTable[] = {
    {"H_", 0.354, 180.0, 2.886, 0.044, 12.000, 0.712, 0.000, 0.00, 4.528},
    {"B_3", 0.838, 109.47, 4.083, 0.180, 12.052, 1.755, 0.000, 2.00, 5.110},
    {"B_2", 0.828, 120.0, 4.083, 0.180, 12.052, 1.755, 0.000, 2.00, 5.110},
    {"C_3", 0.757, 109.47, 3.851, 0.105, 12.730, 1.912, 2.119, 2.00, 5.343},
    {"C_R", 0.729, 120.0, 3.851, 0.105, 12.730, 1.912, 0.000, 2.00, 5.343},
    {"C_2", 0.732, 120.0, 3.851, 0.105, 12.730, 1.912, 0.000, 2.00, 5.343},
    {"C_1", 0.706, 180.0, 3.851, 0.105, 12.730, 1.912, 0.000, 2.00, 5.343},
    {"N_3", 0.700, 106.7, 3.660, 0.069, 13.407, 2.544, 0.450, 2.00, 6.899},
    {"N_R", 0.699, 120.0, 3.660, 0.069, 13.407, 2.544, 0.000, 2.00, 6.899},
    {"N_2", 0.685, 111.2, 3.660, 0.069, 13.407, 2.544, 0.000, 2.00, 6.899},
    {"N_1", 0.656, 180.0, 3.660, 0.069, 13.407, 2.544, 0.000, 2.00, 6.899},
    {"O_3", 0.658, 104.51, 3.500, 0.060, 14.085, 2.300, 0.018, 2.00, 8.741},
    {"O_R", 0.680, 110.0, 3.500, 0.060, 14.085, 2.300, 0.000, 2.00, 8.741},
    {"O_2", 0.634, 120.0, 3.500, 0.060, 14.085, 2.300, 0.000, 2.00, 8.741},
    {"O_1", 0.639, 180.0, 3.500, 0.060, 14.085, 2.300, 0.000, 2.00, 8.741},
    {"F_", 0.668, 180.0, 3.364, 0.050, 14.762, 1.735, 0.000, 2.00, 10.874},
    {"Si3", 1.117, 109.47, 4.295, 0.402, 12.175, 2.323, 1.225, 1.25, 4.168},
    {"P_3+3", 1.101, 93.8, 4.147, 0.305, 13.072, 2.863, 2.400, 1.25, 5.463},
    {"P_3+5", 1.056, 109.47, 4.147, 0.305, 13.072, 2.863, 2.400, 1.25, 5.463},
    {"S_3+2", 1.064, 92.1, 4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928},
    {"S_3+4", 1.049, 103.2, 4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928},
    {"S_3+6", 1.027, 109.47, 4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928},
    {"S_R", 1.077, 92.2, 4.035, 0.274, 13.969, 2.703, 0.000, 1.25, 6.928},
    {"S_2", 0.854, 120.0, 4.035, 0.274, 13.969, 2.703, 0.000, 1.25, 6.928},
    {"Cl", 1.044, 180.0, 3.947, 0.227, 14.866, 2.348, 0.000, 1.25, 8.564},
    {"Br", 1.192, 180.0, 4.189, 0.251, 15.000, 2.519, 0.000, 0.70, 7.790},
    {"I_", 1.382, 180.0, 4.500, 0.339, 15.000, 2.650, 0.000, 0.20, 6.822},
};
constexpr int kTableSize = static_cast<int>(sizeof(kTable) / sizeof(kTable[0]));

int type_of(const char* label) {
  for (int i = 0; i < kTableSize; ++i)
    if (std::string(kTable[i].label) == label) return i;
  return -1;
}

enum class Hyb { Terminal, Sp, Sp2, Sp3 };

Hyb hybridisation(int t) {
  const std::string l = kTable[t].label;
  if (l.size() < 3 || (l[1] != '_' && l != "Si3")) return Hyb::Terminal;  // H_, F_, Cl, Br, I_
  if (l == "Si3") return Hyb::Sp3;
  switch (l[2]) {
    case '1': return Hyb::Sp;
    case '2': case 'R': return Hyb::Sp2;
    case '3': return Hyb::Sp3;
    default: return Hyb::Terminal;
  }
}

bool group16(std::uint8_t z) { return z == 8 || z == 16 || z == 34 || z == 52; }

struct V3 {
  double x, y, z;
};
inline V3 operator-(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V3 operator+(V3 a, V3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline V3 operator*(double s, V3 a) { return {s * a.x, s * a.y, s * a.z}; }
inline double dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline V3 cross(V3 a, V3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline double norm(V3 a) { return std::sqrt(dot(a, a)); }
inline V3 at(const std::vector<double>& x, std::uint32_t i) { return {x[3 * i], x[3 * i + 1], x[3 * i + 2]}; }
inline void add(std::vector<double>& g, std::uint32_t i, V3 v) {
  g[3 * i] += v.x;
  g[3 * i + 1] += v.y;
  g[3 * i + 2] += v.z;
}

bool is_amide_cn(const ChemicalSystemState& st, std::uint32_t c, std::uint32_t n) {
  if (st.atoms.element[c] != 6 || st.atoms.element[n] != 7) return false;
  if ((st.atoms.flags[c] & kAtomAromatic) || (st.atoms.flags[n] & kAtomAromatic)) return false;
  auto nb = st.adjacency.neighbors(c);
  auto bd = st.adjacency.bonds_of(c);
  for (std::size_t k = 0; k < nb.size(); ++k)
    if (st.bonds.order[bd[k]] == 2 && (st.atoms.element[nb[k]] == 8 || st.atoms.element[nb[k]] == 16)) return true;
  return false;
}

}  // namespace

const UffParameters& uff_parameters(int type) { return kTable[type]; }
int uff_type_by_label(const std::string& label) { return type_of(label.c_str()); }

double uff_bond_order(const ChemicalSystemState& st, std::uint32_t b) {
  const std::uint32_t u = st.bonds.begin[b], v = st.bonds.end[b];
  if (st.bonds.flags[b] & kBondAromatic) return 1.5;
  const int o = st.bonds.order[b];
  if (o == 1 && (is_amide_cn(st, u, v) || is_amide_cn(st, v, u))) return 1.41;
  return o;
}

double uff_rest_length(int ti, int tj, double bo) {
  const auto& a = kTable[ti];
  const auto& b = kTable[tj];
  const double rbo = -0.1332 * (a.r1 + b.r1) * std::log(bo);
  const double sx = std::sqrt(a.Xi) - std::sqrt(b.Xi);
  const double ren = a.r1 * b.r1 * sx * sx / (a.Xi * a.r1 + b.Xi * b.r1);
  return a.r1 + b.r1 + rbo - ren;
}

std::vector<int> assign_uff_types(const ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag) {
  const Species& s = st.species[sp];
  std::vector<int> types(s.atom_count, -1);
  for (std::uint32_t g = s.atom_begin; g < s.atom_begin + s.atom_count; ++g) {
    const std::uint8_t z = st.atoms.element[g];
    const bool arom = st.atoms.flags[g] & kAtomAromatic;
    auto nb = st.adjacency.neighbors(g);
    auto bd = st.adjacency.bonds_of(g);
    int doubles = 0, triples = 0, valence = st.atoms.implicit_h[g];
    for (std::size_t k = 0; k < nb.size(); ++k) {
      const int o = st.bonds.order[bd[k]];
      valence += o;
      if (st.bonds.flags[bd[k]] & kBondAromatic) continue;
      doubles += o == 2;
      triples += o == 3;
    }
    const std::size_t degree = nb.size() + st.atoms.implicit_h[g];
    const char* label = nullptr;
    switch (z) {
      case 1: label = "H_"; break;
      case 5: label = degree >= 4 ? "B_3" : "B_2"; break;
      case 6:
        if (arom) label = "C_R";
        else if (triples > 0 || doubles >= 2) label = "C_1";
        else if (doubles == 1 || degree <= 3) label = "C_2";
        else label = "C_3";
        break;
      case 7: {
        bool amide = false;
        for (std::size_t k = 0; k < nb.size() && !amide; ++k)
          amide = st.bonds.order[bd[k]] == 1 && is_amide_cn(st, nb[k], g);
        if (arom) label = "N_R";
        else if (triples > 0 || doubles >= 2) label = "N_1";
        else if (doubles == 1) label = "N_2";
        else if (amide && degree == 3) label = "N_R";
        else label = "N_3";
        break;
      }
      case 8:
        if (arom) label = "O_R";
        else if (triples > 0) label = "O_1";
        else if (doubles > 0) label = "O_2";
        else label = "O_3";
        break;
      case 9: label = "F_"; break;
      case 14: label = "Si3"; break;
      case 15: label = (valence >= 5 || degree >= 4) ? "P_3+5" : "P_3+3"; break;
      case 16:
        if (arom) label = "S_R";
        else if (valence >= 5) label = "S_3+6";
        else if (valence == 4) label = "S_3+4";
        else if (doubles >= 1) label = "S_2";
        else label = "S_3+2";
        break;
      case 17: label = "Cl"; break;
      case 35: label = "Br"; break;
      case 53: label = "I_"; break;
      default: break;
    }
    if (!label) {
      diag.error(DiagCode::NoForceFieldType, "no UFF atom type for " + std::string(element_symbol(z)) + " (atom " +
                                                 std::to_string(g - s.atom_begin) + ")")
          .on_atom(g)
          .on_species(sp);
      continue;
    }
    types[g - s.atom_begin] = type_of(label);
  }
  return types;
}

std::optional<UffForceField> UffForceField::build(const ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag) {
  const Species& s = st.species[sp];
  const std::uint32_t base = s.atom_begin;
  for (std::uint32_t g = base; g < base + s.atom_count; ++g) {
    if (st.atoms.implicit_h[g] > 0) {
      diag.error(DiagCode::ImplicitHydrogens,
                 "species " + s.smiles + " has implicit hydrogens; call add_explicit_hydrogens first")
          .on_species(sp);
      return std::nullopt;
    }
  }
  UffForceField ff;
  ff.n_ = s.atom_count;
  ff.types_ = assign_uff_types(st, sp, diag);
  for (int t : ff.types_)
    if (t < 0) return std::nullopt;
  const auto& T = ff.types_;
  auto local = [&](std::uint32_t g) { return g - base; };

  // Bonds and a local bond-order lookup.
  std::vector<std::vector<std::pair<std::uint32_t, double>>> adj(ff.n_);  // (neighbour, bond order)
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    const std::uint32_t i = local(st.bonds.begin[b]), j = local(st.bonds.end[b]);
    const double bo = uff_bond_order(st, b);
    const double r0 = uff_rest_length(T[i], T[j], bo);
    const double k = 664.12 * kTable[T[i]].Z1 * kTable[T[j]].Z1 / (r0 * r0 * r0);
    ff.bonds_.push_back({i, j, k, r0});
    adj[i].push_back({j, bo});
    adj[j].push_back({i, bo});
  }
  auto bond_order = [&](std::uint32_t a, std::uint32_t b) {
    for (auto [n, bo] : adj[a])
      if (n == b) return bo;
    return 1.0;
  };

  // Angles.
  std::set<std::pair<std::uint32_t, std::uint32_t>> pairs12_13;
  for (const auto& b : ff.bonds_) pairs12_13.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
  for (std::uint32_t j = 0; j < ff.n_; ++j) {
    const double th0 = kTable[T[j]].theta0 * kDeg;
    const double c = std::cos(th0);
    for (std::size_t a = 0; a < adj[j].size(); ++a) {
      for (std::size_t b = a + 1; b < adj[j].size(); ++b) {
        const std::uint32_t i = adj[j][a].first, k = adj[j][b].first;
        pairs12_13.insert({std::min(i, k), std::max(i, k)});
        const double rij = uff_rest_length(T[i], T[j], adj[j][a].second);
        const double rjk = uff_rest_length(T[j], T[k], adj[j][b].second);
        const double rik = std::sqrt(rij * rij + rjk * rjk - 2 * rij * rjk * c);
        const double K = 664.12 * kTable[T[i]].Z1 * kTable[T[k]].Z1 / std::pow(rik, 5) *
                         (3 * rij * rjk * (1 - c * c) - rik * rik * c);
        Angle ang{i, j, k, K, 0, 0, 0, false};
        if (kTable[T[j]].theta0 > 179.0) {
          ang.linear = true;  // E = K (1 + cos θ)
        } else {
          const double sn = std::sin(th0);
          ang.c2 = 1.0 / (4.0 * sn * sn);
          ang.c1 = -4.0 * ang.c2 * c;
          ang.c0 = ang.c2 * (2.0 * c * c + 1.0);
        }
        ff.angles_.push_back(ang);
      }
    }
  }

  // Torsions about each bond j-k.
  for (const auto& bd : ff.bonds_) {
    const std::uint32_t j = bd.i, k = bd.j;
    const Hyb hj = hybridisation(T[j]), hk = hybridisation(T[k]);
    if (hj == Hyb::Terminal || hk == Hyb::Terminal || hj == Hyb::Sp || hk == Hyb::Sp) continue;
    if (adj[j].size() < 2 || adj[k].size() < 2) continue;
    const double bo = bond_order(j, k);
    const std::uint8_t zj = st.atoms.element[base + j], zk = st.atoms.element[base + k];
    double V = 0, phi0 = 0;
    int n = 0;
    auto sp2_sp2 = [&](std::uint32_t a, std::uint32_t b) {
      return 5.0 * std::sqrt(kTable[T[a]].Uj * kTable[T[b]].Uj) * (1.0 + 4.18 * std::log(bo));
    };
    if (hj == Hyb::Sp3 && hk == Hyb::Sp3) {
      if (group16(zj) && group16(zk)) {
        const double vj = zj == 8 ? 2.0 : 6.8, vk = zk == 8 ? 2.0 : 6.8;
        V = std::sqrt(vj * vk);
        n = 2;
        phi0 = 90;
      } else {
        V = std::sqrt(kTable[T[j]].Vi * kTable[T[k]].Vi);
        n = 3;
        phi0 = 180;
      }
    } else if (hj == Hyb::Sp2 && hk == Hyb::Sp2) {
      V = sp2_sp2(j, k);
      n = 2;
      phi0 = 180;
    } else {
      const std::uint32_t s3 = hj == Hyb::Sp3 ? j : k, s2 = hj == Hyb::Sp3 ? k : j;
      const std::uint8_t z3 = st.atoms.element[base + s3];
      bool s2_conjugated = false;
      for (auto [nbr, o] : adj[s2])
        if (nbr != s3 && hybridisation(T[nbr]) == Hyb::Sp2) s2_conjugated = true;
      if (group16(z3)) {
        V = sp2_sp2(s3, s2);
        n = 2;
        phi0 = 90;
      } else if (s2_conjugated) {
        V = 2.0;
        n = 3;
        phi0 = 180;
      } else {
        V = 1.0;
        n = 6;
        phi0 = 0;
      }
    }
    if (V == 0.0) continue;
    const double share = static_cast<double>((adj[j].size() - 1) * (adj[k].size() - 1));
    for (auto [i, o1] : adj[j]) {
      if (i == k) continue;
      for (auto [l, o2] : adj[k]) {
        if (l == j || l == i) continue;
        ff.torsions_.push_back({i, j, k, l, V / share, std::cos(n * phi0 * kDeg), n});
      }
    }
  }

  // Inversion at sp2 C and N with three neighbours.
  for (std::uint32_t c = 0; c < ff.n_; ++c) {
    if (adj[c].size() != 3) continue;
    const std::string lab = kTable[T[c]].label;
    if (lab != "C_2" && lab != "C_R" && lab != "N_2" && lab != "N_R") continue;
    double K = 6.0;
    if (lab[0] == 'C')
      for (auto [nbr, o] : adj[c])
        if (std::string(kTable[T[nbr]].label) == "O_2") K = 50.0;
    const std::uint32_t a = adj[c][0].first, b = adj[c][1].first, d = adj[c][2].first;
    // Each neighbour takes a turn as the out-of-plane atom; K is shared between the three.
    ff.inversions_.push_back({c, a, b, d, K / 3.0, 1.0, -1.0, 0.0});
    ff.inversions_.push_back({c, b, d, a, K / 3.0, 1.0, -1.0, 0.0});
    ff.inversions_.push_back({c, d, a, b, K / 3.0, 1.0, -1.0, 0.0});
  }

  // van der Waals beyond 1-3.
  for (std::uint32_t i = 0; i < ff.n_; ++i) {
    for (std::uint32_t j = i + 1; j < ff.n_; ++j) {
      if (pairs12_13.count({i, j})) continue;
      const double D = std::sqrt(kTable[T[i]].D1 * kTable[T[j]].D1);
      const double x = std::sqrt(kTable[T[i]].x1 * kTable[T[j]].x1);
      ff.vdw_.push_back({i, j, D, x});
    }
  }
  return ff;
}

double UffForceField::eval(const std::vector<double>& X, std::vector<double>* G, UffEnergyBreakdown* parts) const {
  UffEnergyBreakdown e;
  if (G) G->assign(3 * n_, 0.0);

  for (const auto& b : bonds_) {
    const V3 d = at(X, b.i) - at(X, b.j);
    const double r = norm(d);
    const double dr = r - b.r0;
    e.bond += 0.5 * b.k * dr * dr;
    if (G && r > 0) {
      const V3 f = (b.k * dr / r) * d;
      add(*G, b.i, f);
      add(*G, b.j, -1.0 * f);
    }
  }

  for (const auto& a : angles_) {
    const V3 u = at(X, a.i) - at(X, a.j), v = at(X, a.k) - at(X, a.j);
    const double lu = norm(u), lv = norm(v);
    if (lu == 0 || lv == 0) continue;
    double c = dot(u, v) / (lu * lv);
    c = std::clamp(c, -1.0, 1.0);
    double dE_dc;
    if (a.linear) {
      e.angle += a.K * (1.0 + c);
      dE_dc = a.K;
    } else {
      e.angle += a.K * (a.c0 + a.c1 * c + a.c2 * (2 * c * c - 1));
      dE_dc = a.K * (a.c1 + 4 * a.c2 * c);
    }
    if (G) {
      const V3 gi = (dE_dc / lu) * ((1.0 / lv) * v - (c / lu) * u);
      const V3 gk = (dE_dc / lv) * ((1.0 / lu) * u - (c / lv) * v);
      add(*G, a.i, gi);
      add(*G, a.k, gk);
      add(*G, a.j, -1.0 * (gi + gk));
    }
  }

  for (const auto& t : torsions_) {
    // Blondel & Karplus, J. Comput. Chem. 17 (1996) 1132: F = r1-r2, G = r2-r3, H = r4-r3,
    // A = F×G, B = H×G, cos φ = A·B/|A||B|, sin φ = (B×A)·G/(|A||B||G|).
    const V3 p1 = at(X, t.i), p2 = at(X, t.j), p3 = at(X, t.k), p4 = at(X, t.l);
    const V3 F = p1 - p2, Gv = p2 - p3, H = p4 - p3;
    const V3 A = cross(F, Gv), Bv = cross(H, Gv);
    const double A2 = dot(A, A), B2 = dot(Bv, Bv), lg = norm(Gv);
    if (A2 < 1e-20 || B2 < 1e-20 || lg == 0) continue;
    const double phi = std::atan2(dot(cross(Bv, A), Gv) / lg, dot(A, Bv));
    e.torsion += 0.5 * t.V * (1.0 - t.cos_nphi0 * std::cos(t.n * phi));
    if (G) {
      const double dE = 0.5 * t.V * t.cos_nphi0 * t.n * std::sin(t.n * phi);
      const double fg = dot(F, Gv), hg = dot(H, Gv);
      const V3 d1 = (-lg / A2) * A;
      const V3 d4 = (lg / B2) * Bv;
      const V3 d2 = (lg / A2) * A + (fg / (A2 * lg)) * A - (hg / (B2 * lg)) * Bv;
      const V3 d3 = (-lg / B2) * Bv - (fg / (A2 * lg)) * A + (hg / (B2 * lg)) * Bv;
      add(*G, t.i, dE * d1);
      add(*G, t.j, dE * d2);
      add(*G, t.k, dE * d3);
      add(*G, t.l, dE * d4);
    }
  }

  for (const auto& w : inversions_) {
    // ω: angle between bond center->l and the plane (center, j, k); s = sin ω.
    const V3 c = at(X, w.center);
    const V3 rj = at(X, w.j) - c, rk = at(X, w.k) - c, rl = at(X, w.l) - c;
    const V3 nrm = cross(rj, rk);
    const double gn = norm(nrm), hl = norm(rl);
    if (gn < 1e-12 || hl < 1e-12) continue;
    const double f = dot(nrm, rl);
    double s = std::clamp(f / (gn * hl), -1.0, 1.0);
    const double cw = std::sqrt(std::max(0.0, 1.0 - s * s));
    e.inversion += w.K * (w.c0 + w.c1 * cw + w.c2 * (2 * cw * cw - 1));
    if (G) {
      // dE/ds with cos ω = sqrt(1 - s²); clamp near |s| = 1 where cos ω -> 0.
      const double dE_dcw = w.K * (w.c1 + 4 * w.c2 * cw);
      const double dcw_ds = cw > 1e-8 ? -s / cw : 0.0;
      const double dE_ds = dE_dcw * dcw_ds;
      const V3 df_rj = cross(rk, rl), df_rk = cross(rl, rj), df_rl = nrm;
      const V3 dg_rj = (1.0 / gn) * cross(rk, nrm), dg_rk = (1.0 / gn) * cross(nrm, rj);
      const V3 dh_rl = (1.0 / hl) * rl;
      const double inv = 1.0 / (gn * hl);
      const V3 ds_rj = inv * df_rj - (f / (gn * gn * hl)) * dg_rj;
      const V3 ds_rk = inv * df_rk - (f / (gn * gn * hl)) * dg_rk;
      const V3 ds_rl = inv * df_rl - (f / (gn * hl * hl)) * dh_rl;
      add(*G, w.j, dE_ds * ds_rj);
      add(*G, w.k, dE_ds * ds_rk);
      add(*G, w.l, dE_ds * ds_rl);
      add(*G, w.center, -dE_ds * (ds_rj + ds_rk + ds_rl));
    }
  }

  for (const auto& p : vdw_) {
    const V3 d = at(X, p.i) - at(X, p.j);
    const double r2 = std::max(dot(d, d), 1e-6);
    const double q2 = p.x * p.x / r2;
    const double q6 = q2 * q2 * q2;
    e.vdw += p.D * (q6 * q6 - 2.0 * q6);
    if (G) {
      // dE/dr * (1/r) = D (-12 q12 + 12 q6) / r²
      const double f = p.D * (-12.0 * q6 * q6 + 12.0 * q6) / r2;
      add(*G, p.i, f * d);
      add(*G, p.j, -f * d);
    }
  }
  if (parts) *parts = e;
  return e.total();
}

double UffForceField::energy(const std::vector<double>& xyz) const { return eval(xyz, nullptr, nullptr); }

double UffForceField::energy_gradient(const std::vector<double>& xyz, std::vector<double>& grad) const {
  return eval(xyz, &grad, nullptr);
}

UffEnergyBreakdown UffForceField::breakdown(const std::vector<double>& xyz) const {
  UffEnergyBreakdown e;
  eval(xyz, nullptr, &e);
  return e;
}

UffMinimizeResult uff_minimize(ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag, const LbfgsOptions& options) {
  UffMinimizeResult res;
  auto ff = UffForceField::build(st, sp, diag);
  if (!ff) return res;
  const Species& s = st.species[sp];
  std::vector<double> x(3 * s.atom_count);
  for (std::uint32_t i = 0; i < s.atom_count; ++i) {
    x[3 * i] = st.atoms.x[s.atom_begin + i];
    x[3 * i + 1] = st.atoms.y[s.atom_begin + i];
    x[3 * i + 2] = st.atoms.z[s.atom_begin + i];
  }
  res.energy_initial = ff->energy(x);
  auto r = minimize_lbfgs([&](const std::vector<double>& xx, std::vector<double>& g) { return ff->energy_gradient(xx, g); },
                          x, options);
  for (std::uint32_t i = 0; i < s.atom_count; ++i) {
    st.atoms.x[s.atom_begin + i] = x[3 * i];
    st.atoms.y[s.atom_begin + i] = x[3 * i + 1];
    st.atoms.z[s.atom_begin + i] = x[3 * i + 2];
  }
  ++st.revision;
  res.ok = true;
  res.converged = r.converged;
  res.iterations = r.iterations;
  res.energy_final = r.f;
  res.grad_max = r.grad_max;
  if (!r.converged)
    diag.warn(DiagCode::MinimizationNotConverged, "UFF minimisation stopped after " + std::to_string(r.iterations) +
                                                      " iterations (max |gradient| " + std::to_string(r.grad_max) + ")")
        .on_species(sp);
  return res;
}

}  // namespace chemlarp
