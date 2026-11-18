// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
// Protein-ligand docking: see include/chemlarp/docking.hpp for the model and its references.

#include "chemlarp/docking.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <string>
#include <utility>

#include "chemlarp/conformer.hpp"
#include "chemlarp/elements.hpp"
#include "chemlarp/json.hpp"
#include "chemlarp/molio.hpp"
#include "chemlarp/parallel.hpp"
#include "chemlarp/smiles.hpp"

namespace chemlarp {
namespace {

constexpr double kPi = 3.14159265358979323846;
// Weight of the quadratic penalty that keeps a ligand atom inside the search box, kcal/mol per A^2. The
// repulsion term's weight is 0.84 per A^2, so this is of the same order: enough to confine, smooth enough to
// minimise through.
constexpr double kBoxPenalty = 1.0;
// How far the grid reaches beyond the search box, so that an atom in the penalty region still reads a
// receptor energy instead of a clamped one.
constexpr double kGridMargin = 2.0;
// Cap on the graph matchings the symmetry-corrected RMSD enumerates, so that a highly symmetric molecule
// cannot turn the search into a combinatorial explosion. The best matching found so far prunes the rest.
constexpr std::size_t kMaxMatchings = 200000;
// While the Monte Carlo hunts, repulsion is squashed smoothly towards these ceilings (kcal/mol) instead of
// standing up as a wall: e -> e·c/(c + e) for e > 0. A random pose that lands inside the protein can then
// slide out, and a ligand can reach a buried site through a narrow mouth, instead of being thrown straight
// back by a gradient of hundreds. The receptor ceiling applies to a ligand atom's whole interaction with
// the receptor and the ligand ceiling to each of the ligand's own pairs, which is how the published search
// applies them. Every pose that is kept is then minimised again on the function as it is, which is also the
// one reported, so the smoothing changes where the search looks and never what it reports.
constexpr double kHuntCapIntermolecular = 1.5;
constexpr double kHuntCapInternal = 10.0;
// Most points the grid maps may hold per atom type, which is what bounds the memory a job asks for: at eight
// bytes a point this is 32 MB a map, and a box of 55 A a side at the default 0.375 A spacing fits inside it.
// A larger box gets a coarser spacing rather than a larger allocation.
constexpr std::size_t kMaxGridPoints = 4000000;

// Smooth saturation of a repulsive energy at `cap`, and the factor its derivative picks up.
double squash(double e, double cap, double* derivative_factor = nullptr) {
  if (!(cap > 0.0) || !(e > 0.0)) {
    if (derivative_factor) *derivative_factor = 1.0;
    return e;
  }
  const double f = cap / (cap + e);
  if (derivative_factor) *derivative_factor = f * f;
  return e * f;
}

// ---- small vector helpers -----------------------------------------------------------------------

using Vec3 = std::array<double, 3>;

Vec3 operator-(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 operator+(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 operator*(double s, const Vec3& a) { return {s * a[0], s * a[1], s * a[2]}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double length(const Vec3& a) { return std::sqrt(dot(a, a)); }
Vec3 normalized(const Vec3& a) {
  const double n = length(a);
  return n > 1e-12 ? (1.0 / n) * a : Vec3{1, 0, 0};
}

// Rotation of `v` about the unit axis `a` by `angle` (right-handed), Rodrigues' formula.
Vec3 rotate_about(const Vec3& v, const Vec3& a, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  return c * v + s * cross(a, v) + (1.0 - c) * dot(a, v) * a;
}

// ---- deterministic random numbers ---------------------------------------------------------------

// splitmix64; each Monte Carlo run has its own generator seeded from the run index, so the result does not
// depend on how the runs are spread over threads.
struct Rng {
  std::uint64_t s = 0;
  explicit Rng(std::uint64_t seed) : s(seed + 0x9E3779B97F4A7C15ull) {}
  std::uint64_t next() {
    std::uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }
  double uniform() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }  // [0, 1)
  double uniform(double lo, double hi) { return lo + (hi - lo) * uniform(); }
  std::size_t index(std::size_t n) { return n ? static_cast<std::size_t>(next() % n) : 0; }
  double normal() {  // Box-Muller
    const double u = std::max(uniform(), 1e-16), v = uniform();
    return std::sqrt(-2.0 * std::log(u)) * std::cos(2.0 * kPi * v);
  }
  Vec3 in_unit_ball() {
    for (int k = 0; k < 100; ++k) {
      const Vec3 p{uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)};
      if (dot(p, p) <= 1.0) return p;
    }
    return {0, 0, 0};
  }
  // Uniform over the sphere of rotations (Shoemake): a random unit quaternion.
  std::array<double, 4> orientation() {
    const double u1 = uniform(), u2 = uniform(), u3 = uniform();
    const double r1 = std::sqrt(1.0 - u1), r2 = std::sqrt(u1);
    return {r1 * std::sin(2 * kPi * u2), r1 * std::cos(2 * kPi * u2), r2 * std::sin(2 * kPi * u3),
            r2 * std::cos(2 * kPi * u3)};
  }
};

// ---- quaternions and the exponential map --------------------------------------------------------

using Quat = std::array<double, 4>;  // (w, x, y, z)

Quat quat_identity() { return {1, 0, 0, 0}; }

Quat quat_multiply(const Quat& a, const Quat& b) {
  return {a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3], a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
          a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1], a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]};
}

Quat quat_normalize(const Quat& q) {
  const double n = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  if (!(n > 1e-12)) return quat_identity();
  return {q[0] / n, q[1] / n, q[2] / n, q[3] / n};
}

// The rotation vector r (axis * angle) as a quaternion.
Quat quat_from_rotation_vector(const Vec3& r) {
  const double t = length(r);
  if (t < 1e-9) return quat_normalize({1, 0.5 * r[0], 0.5 * r[1], 0.5 * r[2]});
  const double s = std::sin(0.5 * t) / t;
  return {std::cos(0.5 * t), s * r[0], s * r[1], s * r[2]};
}

using Mat3 = std::array<double, 9>;  // row-major

Mat3 matrix_from_quat(const Quat& q) {
  const double w = q[0], x = q[1], y = q[2], z = q[3];
  return {1 - 2 * (y * y + z * z), 2 * (x * y - w * z),     2 * (x * z + w * y),
          2 * (x * y + w * z),     1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
          2 * (x * z - w * y),     2 * (y * z + w * x),     1 - 2 * (x * x + y * y)};
}

Vec3 mat_apply(const Mat3& m, const Vec3& v) {
  return {m[0] * v[0] + m[1] * v[1] + m[2] * v[2], m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
          m[6] * v[0] + m[7] * v[1] + m[8] * v[2]};
}

Vec3 mat_apply_t(const Mat3& m, const Vec3& v) {
  return {m[0] * v[0] + m[3] * v[1] + m[6] * v[2], m[1] * v[0] + m[4] * v[1] + m[7] * v[2],
          m[2] * v[0] + m[5] * v[1] + m[8] * v[2]};
}

// Transpose of the right Jacobian of SO(3) applied to w:  J_r(r)^T w = w + a (r x w) + b (r x (r x w)),
// a = (1 - cos t) / t^2, b = (t - sin t) / t^3 (series near 0). With it, the gradient of an energy with
// respect to the rotation vector r is J_r(r)^T sum_i ( v_i x (R(r)^T g_i) ), which is exact for every r,
// not only at r = 0.
Vec3 right_jacobian_transpose(const Vec3& r, const Vec3& w) {
  const double t2 = dot(r, r), t = std::sqrt(t2);
  double a, b;
  if (t < 1e-5) {
    a = 0.5 - t2 / 24.0;
    b = 1.0 / 6.0 - t2 / 120.0;
  } else {
    a = (1.0 - std::cos(t)) / t2;
    b = (t - std::sin(t)) / (t2 * t);
  }
  const Vec3 rxw = cross(r, w);
  return w + a * rxw + b * cross(r, rxw);
}

// ---- cell list ----------------------------------------------------------------------------------

// A uniform grid of cells over a set of points, for "every point within R" queries. Sizes are computed in
// double and clamped before any conversion to int, so absurd coordinates (a fuzzed PDB file) cannot overflow.
class CellList {
 public:
  void build(const std::vector<Vec3>& points, double cell) {
    points_ = &points;
    cell_ = std::max(cell, 1e-3);
    lo_ = {0, 0, 0};
    n_ = {1, 1, 1};
    start_.assign(2, 0);
    items_.clear();
    if (points.empty()) return;
    Vec3 hi = points[0];
    lo_ = points[0];
    for (const Vec3& p : points)
      for (int k = 0; k < 3; ++k) {
        lo_[k] = std::min(lo_[k], p[k]);
        hi[k] = std::max(hi[k], p[k]);
      }
    for (int k = 0; k < 3; ++k) {
      const double span = hi[k] - lo_[k];
      const double cells = std::isfinite(span) ? std::floor(span / cell_) + 1.0 : 1.0;
      n_[k] = static_cast<int>(std::min(std::max(cells, 1.0), 512.0));
    }
    const std::size_t total = static_cast<std::size_t>(n_[0]) * n_[1] * n_[2];
    std::vector<std::uint32_t> count(total + 1, 0);
    std::vector<std::uint32_t> cell_of(points.size(), 0);
    for (std::size_t i = 0; i < points.size(); ++i) {
      const std::uint32_t c = cell_index(points[i]);
      cell_of[i] = c;
      ++count[c + 1];
    }
    for (std::size_t c = 0; c < total; ++c) count[c + 1] += count[c];
    start_ = count;
    items_.resize(points.size());
    std::vector<std::uint32_t> fill(start_.begin(), start_.end() - 1);
    for (std::size_t i = 0; i < points.size(); ++i) items_[fill[cell_of[i]]++] = static_cast<std::uint32_t>(i);
  }

  // Calls body(index, r2) for every stored point within `radius` of p.
  template <class F>
  void for_each_near(const Vec3& p, double radius, F&& body) const {
    if (!points_ || points_->empty()) return;
    int lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
      const double a = std::floor((p[k] - radius - lo_[k]) / cell_);
      const double b = std::floor((p[k] + radius - lo_[k]) / cell_);
      if (!(std::isfinite(a) && std::isfinite(b))) return;
      lo[k] = static_cast<int>(std::min(std::max(a, 0.0), static_cast<double>(n_[k] - 1)));
      hi[k] = static_cast<int>(std::min(std::max(b, 0.0), static_cast<double>(n_[k] - 1)));
      if (b < 0.0 || a > static_cast<double>(n_[k] - 1)) return;
    }
    const double r2max = radius * radius;
    for (int ix = lo[0]; ix <= hi[0]; ++ix)
      for (int iy = lo[1]; iy <= hi[1]; ++iy)
        for (int iz = lo[2]; iz <= hi[2]; ++iz) {
          const std::size_t c = (static_cast<std::size_t>(ix) * n_[1] + iy) * n_[2] + iz;
          for (std::uint32_t k = start_[c]; k < start_[c + 1]; ++k) {
            const std::uint32_t i = items_[k];
            const Vec3 d = (*points_)[i] - p;
            const double r2 = dot(d, d);
            if (r2 <= r2max) body(i, r2);
          }
        }
  }

 private:
  std::uint32_t cell_index(const Vec3& p) const {
    std::size_t idx[3];
    for (int k = 0; k < 3; ++k) {
      const double f = std::floor((p[k] - lo_[k]) / cell_);
      const double c = std::isfinite(f) ? std::min(std::max(f, 0.0), static_cast<double>(n_[k] - 1)) : 0.0;
      idx[k] = static_cast<std::size_t>(c);
    }
    return static_cast<std::uint32_t>((idx[0] * n_[1] + idx[1]) * n_[2] + idx[2]);
  }

  const std::vector<Vec3>* points_ = nullptr;
  double cell_ = 8.0;
  Vec3 lo_{0, 0, 0};
  std::array<int, 3> n_{1, 1, 1};
  std::vector<std::uint32_t> start_;
  std::vector<std::uint32_t> items_;
};

// ---- the scored elements ------------------------------------------------------------------------

// Metal ions that appear in protein structures and are typed as donors (Vina does the same for these).
bool is_scored_metal(std::uint8_t z) {
  switch (z) {
    case 11:  // Na
    case 12:  // Mg
    case 19:  // K
    case 20:  // Ca
    case 25:  // Mn
    case 26:  // Fe
    case 27:  // Co
    case 28:  // Ni
    case 29:  // Cu
    case 30:  // Zn
    case 48:  // Cd
    case 80:  // Hg
      return true;
    default:
      return false;
  }
}

}  // namespace

// ---- atom types and the scoring function --------------------------------------------------------

const char* to_string(XsType t) noexcept {
  switch (t) {
    case XsType::C_H: return "C_H";
    case XsType::C_P: return "C_P";
    case XsType::N_P: return "N_P";
    case XsType::N_D: return "N_D";
    case XsType::N_A: return "N_A";
    case XsType::N_DA: return "N_DA";
    case XsType::O_A: return "O_A";
    case XsType::O_DA: return "O_DA";
    case XsType::S_P: return "S_P";
    case XsType::P_P: return "P_P";
    case XsType::F_H: return "F_H";
    case XsType::Cl_H: return "Cl_H";
    case XsType::Br_H: return "Br_H";
    case XsType::I_H: return "I_H";
    case XsType::Met_D: return "Met_D";
    case XsType::None: break;
  }
  return "none";
}

double xs_radius(XsType t) noexcept {
  switch (t) {
    case XsType::C_H:
    case XsType::C_P: return 1.9;
    case XsType::N_P:
    case XsType::N_D:
    case XsType::N_A:
    case XsType::N_DA: return 1.8;
    case XsType::O_A:
    case XsType::O_DA: return 1.7;
    case XsType::S_P: return 2.0;
    case XsType::P_P: return 2.1;
    case XsType::F_H: return 1.5;
    case XsType::Cl_H: return 1.8;
    case XsType::Br_H: return 2.0;
    case XsType::I_H: return 2.2;
    case XsType::Met_D: return 1.2;
    case XsType::None: break;
  }
  return 0.0;
}

bool xs_hydrophobic(XsType t) noexcept {
  return t == XsType::C_H || t == XsType::F_H || t == XsType::Cl_H || t == XsType::Br_H || t == XsType::I_H;
}

bool xs_donor(XsType t) noexcept {
  return t == XsType::N_D || t == XsType::N_DA || t == XsType::O_DA || t == XsType::Met_D;
}

bool xs_acceptor(XsType t) noexcept {
  return t == XsType::N_A || t == XsType::N_DA || t == XsType::O_A || t == XsType::O_DA;
}

XsType xs_type(std::uint8_t z, bool heteroatom_neighbour, int hydrogens, int connections) noexcept {
  switch (z) {
    case 6: return heteroatom_neighbour ? XsType::C_P : XsType::C_H;
    case 7: {
      // A nitrogen accepts only when it still has a lone pair free, which the number of bonds stands in for:
      // an amide, ammonium, amidine or pyrrole nitrogen (three or more bonds) donates but does not accept,
      // while a pyridine, imidazole or nitrile nitrogen (two bonds or fewer) accepts. This is the rule the
      // AutoDock atom typer uses, crude for a tertiary amine but the same as the reference program's.
      const bool donor = hydrogens > 0, acceptor = connections <= 2;
      if (donor) return acceptor ? XsType::N_DA : XsType::N_D;
      return acceptor ? XsType::N_A : XsType::N_P;
    }
    case 8: return hydrogens > 0 ? XsType::O_DA : XsType::O_A;
    case 16: return XsType::S_P;
    case 15: return XsType::P_P;
    case 9: return XsType::F_H;
    case 17: return XsType::Cl_H;
    case 35: return XsType::Br_H;
    case 53: return XsType::I_H;
    default: return is_scored_metal(z) ? XsType::Met_D : XsType::None;
  }
}

double ScoringTerms::weighted(const ScoringWeights& w) const noexcept {
  return w.gauss1 * gauss1 + w.gauss2 * gauss2 + w.repulsion * repulsion + w.hydrophobic * hydrophobic +
         w.hbond * hbond;
}

ScoringTerms pair_terms(XsType a, XsType b, double r, const ScoringWeights& w, double* derivative) {
  ScoringTerms t;
  if (derivative) *derivative = 0.0;
  if (a == XsType::None || b == XsType::None) return t;
  if (!(r >= 0.0) || r >= w.cutoff) return t;  // a NaN distance scores nothing rather than poisoning the sum
  const double d = r - xs_radius(a) - xs_radius(b);

  // gauss1: exp(-(d / 0.5)^2)
  const double g1 = std::exp(-4.0 * d * d);
  double dg1 = -8.0 * d * g1;
  // gauss2: exp(-((d - 3) / 2)^2)
  const double e2 = (d - 3.0) * 0.5;
  const double g2 = std::exp(-e2 * e2);
  double dg2 = -0.5 * (d - 3.0) * g2;
  // repulsion: d^2 for d < 0
  const double rep = d < 0.0 ? d * d : 0.0;
  double drep = d < 0.0 ? 2.0 * d : 0.0;

  double hyd = 0.0, dhyd = 0.0;
  if (xs_hydrophobic(a) && xs_hydrophobic(b)) {
    if (d < 0.5) {
      hyd = 1.0;
    } else if (d < 1.5) {
      hyd = 1.5 - d;
      dhyd = -1.0;
    }
  }
  double hb = 0.0, dhb = 0.0;
  if ((xs_donor(a) && xs_acceptor(b)) || (xs_acceptor(a) && xs_donor(b))) {
    if (d < -0.7) {
      hb = 1.0;
    } else if (d < 0.0) {
      hb = -d / 0.7;
      dhb = -1.0 / 0.7;
    }
  }
  t.gauss1 = g1;
  t.gauss2 = g2;
  t.repulsion = rep;
  t.hydrophobic = hyd;
  t.hbond = hb;
  if (derivative)
    *derivative = w.gauss1 * dg1 + w.gauss2 * dg2 + w.repulsion * drep + w.hydrophobic * dhyd + w.hbond * dhb;
  return t;
}

double pair_energy(XsType a, XsType b, double r, const ScoringWeights& w, double* derivative) {
  return pair_terms(a, b, r, w, derivative).weighted(w);
}

// ---- receptor -----------------------------------------------------------------------------------

namespace {

bool is_water(const std::string& name) { return name == "HOH" || name == "DOD" || name == "WAT" || name == "H2O"; }

// Is residue `i` joined to the residue before or after it by a peptide bond? Modified residues
// (selenomethionine, S-hydroxycysteine, phosphoserine) are deposited as HETATM, and without this one
// sitting next to the ligand would be deleted along with the ligand's site-mates, leaving a hole in the
// protein. Having backbone atom names is not enough to tell: a ligand can have them too — methotrexate's
// glutamate carries N, CA and C — so the test is the bond itself.
bool bonded_into_chain(const Protein& p, std::size_t i) {
  const ProteinResidue& r = p.residues[i];
  const ProteinAtom* n = r.atom("N");
  const ProteinAtom* c = r.atom("C");
  if (!n && !c) return false;
  auto close = [](const ProteinAtom* a, const ProteinAtom* b) {
    if (!a || !b) return false;
    const double dx = a->x - b->x, dy = a->y - b->y, dz = a->z - b->z;
    return dx * dx + dy * dy + dz * dz < 1.5 * 1.5;  // a peptide bond is 1.33 A
  };
  if (i > 0 && p.residues[i - 1].chain == r.chain && close(n, p.residues[i - 1].atom("C"))) return true;
  const std::size_t j = i + 1;
  return j < p.residues.size() && p.residues[j].chain == r.chain && close(c, p.residues[j].atom("N"));
}

std::string residue_label(const ProteinResidue& r) {
  std::string s = r.chain.empty() ? std::string("_") : r.chain;
  s += " " + std::to_string(r.number);
  if (r.icode != ' ') s += r.icode;
  return s + " " + r.name;
}

std::uint8_t element_of(const ProteinAtom& a) {
  if (!a.element.empty())
    if (const ElementInfo* e = find_element_by_symbol(a.element)) return e->z;
  // Fall back to the first letter of the atom name (a file with a blank element column that the reader
  // could not resolve either).
  for (char c : a.name)
    if (std::isalpha(static_cast<unsigned char>(c))) {
      const std::string s(1, static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
      if (const ElementInfo* e = find_element_by_symbol(s)) return e->z;
      break;
    }
  return 0;
}

}  // namespace

PointCloud protein_ligand_atoms(const Protein& p, std::string_view residue_name) {
  PointCloud out;
  for (const ProteinResidue& r : p.residues) {
    if (r.name != residue_name) continue;
    for (const ProteinAtom& a : r.atoms) {
      const std::uint8_t z = element_of(a);
      if (z == 0 || z == 1) continue;
      out.push_back({a.x, a.y, a.z});
    }
  }
  return out;
}

DockBox box_around(const PointCloud& points, double margin) {
  DockBox b;
  if (points.empty()) return b;
  Vec3 lo{points[0][0], points[0][1], points[0][2]}, hi = lo;
  for (const auto& p : points)
    for (int k = 0; k < 3; ++k) {
      lo[k] = std::min(lo[k], p[k]);
      hi[k] = std::max(hi[k], p[k]);
    }
  b.cx = 0.5 * (lo[0] + hi[0]);
  b.cy = 0.5 * (lo[1] + hi[1]);
  b.cz = 0.5 * (lo[2] + hi[2]);
  b.sx = hi[0] - lo[0] + 2 * margin;
  b.sy = hi[1] - lo[1] + 2 * margin;
  b.sz = hi[2] - lo[2] + 2 * margin;
  return b;
}

DockReceptor prepare_receptor(const Protein& protein, const ReceptorOptions& options, const PointCloud& reference,
                              Diagnostics& diag) {
  DockReceptor out;
  Protein p = protein;
  if (options.repair) {
    RepairOptions ro;
    ro.loops = false;
    ro.termini = false;
    repair(p, ro, diag);
  }
  if (options.protonate) {
    ProtonationOptions po;
    po.pH = options.pH;
    protonate(p, po, diag);
  }

  // Which residues to keep. Waters go; a hetero group with an atom close to the reference ligand is either
  // that ligand or something sharing its site, so it goes too. Everything else (haem, metals, cofactors,
  // the chains) stays.
  std::vector<Vec3> ref;
  ref.reserve(reference.size());
  for (const auto& q : reference) ref.push_back({q[0], q[1], q[2]});
  CellList ref_cells;
  ref_cells.build(ref, std::max(options.ligand_overlap, 1.0));

  std::vector<char> keep(p.residues.size(), 1);
  for (std::size_t i = 0; i < p.residues.size(); ++i) {
    const ProteinResidue& r = p.residues[i];
    if (options.remove_waters && is_water(r.name)) {
      keep[i] = 0;
      ++out.waters_removed;
      continue;
    }
    if (!r.hetero || ref.empty() || bonded_into_chain(p, i)) continue;
    bool overlaps = false;
    for (const ProteinAtom& a : r.atoms) {
      if (overlaps) break;
      const std::uint8_t z = element_of(a);
      if (z == 0 || z == 1) continue;
      ref_cells.for_each_near({a.x, a.y, a.z}, options.ligand_overlap, [&](std::uint32_t, double) { overlaps = true; });
    }
    if (overlaps) {
      keep[i] = 0;
      ++out.hetero_removed;
    }
  }
  for (std::size_t i = 0; i < p.residues.size(); ++i)
    if (keep[i] && p.residues[i].hetero && !is_water(p.residues[i].name) && !bonded_into_chain(p, i))
      out.kept_hetero.push_back(residue_label(p.residues[i]));

  // Collect the atoms that survive, heavy atoms and hydrogens apart.
  struct Raw {
    Vec3 pos;
    std::uint8_t z;
    std::size_t residue;
    const ProteinAtom* atom;
  };
  std::vector<Raw> heavy;
  std::vector<Vec3> hydrogens;
  for (std::size_t i = 0; i < p.residues.size(); ++i) {
    if (!keep[i]) continue;
    for (const ProteinAtom& a : p.residues[i].atoms) {
      if (!(std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z))) continue;
      const std::uint8_t z = element_of(a);
      if (z == 0) continue;
      if (z == 1) {
        hydrogens.push_back({a.x, a.y, a.z});
      } else {
        heavy.push_back({{a.x, a.y, a.z}, z, i, &a});
      }
    }
  }

  // Typing needs connectivity, which a PDB file does not carry: a hydrogen within 1.45 A makes a donor
  // (N-H 1.01, O-H 0.96, S-H 1.34), a nitrogen or oxygen within 1.75 A makes a carbon polar
  // (C-N 1.47, C-O 1.43, aromatic C-N 1.34), and a nitrogen's heavy neighbours within 1.8 A decide whether
  // it still has a lone pair to accept with.
  std::vector<Vec3> heavy_pos(heavy.size());
  for (std::size_t i = 0; i < heavy.size(); ++i) heavy_pos[i] = heavy[i].pos;
  CellList h_cells, heavy_cells;
  h_cells.build(hydrogens, 2.0);
  heavy_cells.build(heavy_pos, 2.0);

  out.atoms.reserve(heavy.size());
  for (std::size_t i = 0; i < heavy.size(); ++i) {
    int hcount = 0;
    if (heavy[i].z == 7 || heavy[i].z == 8 || heavy[i].z == 16)
      h_cells.for_each_near(heavy[i].pos, 1.45, [&](std::uint32_t, double) { ++hcount; });
    bool hetero_neighbour = false;
    int heavy_neighbours = 0;
    if (heavy[i].z == 6)
      heavy_cells.for_each_near(heavy[i].pos, 1.75, [&](std::uint32_t j, double r2) {
        if (j != i && r2 > 0.16 && (heavy[j].z == 7 || heavy[j].z == 8)) hetero_neighbour = true;
      });
    if (heavy[i].z == 7)
      heavy_cells.for_each_near(heavy[i].pos, 1.8, [&](std::uint32_t j, double r2) {
        if (j != i && r2 > 0.16) ++heavy_neighbours;
      });
    const XsType t = xs_type(heavy[i].z, hetero_neighbour, hcount, heavy_neighbours + hcount);
    if (t == XsType::None) continue;
    ReceptorAtom ra;
    ra.x = heavy[i].pos[0];
    ra.y = heavy[i].pos[1];
    ra.z = heavy[i].pos[2];
    ra.type = t;
    ra.residue = residue_label(p.residues[heavy[i].residue]);
    ra.name = heavy[i].atom->name;
    out.atoms.push_back(std::move(ra));
    out.donors += xs_donor(t);
    out.acceptors += xs_acceptor(t);
    out.hydrophobic += xs_hydrophobic(t);
  }
  if (out.atoms.empty()) diag.error(DiagCode::InvalidInput, "the receptor has no atoms the scoring function knows");
  return out;
}

// ---- ligand -------------------------------------------------------------------------------------

namespace {

// A bond is in a ring exactly when it is not a bridge of the molecular graph. Iterative Tarjan so that a
// long chain cannot overflow the stack.
std::vector<char> ring_bonds(const std::vector<std::pair<std::uint32_t, std::uint32_t>>& bonds, std::size_t atoms) {
  std::vector<char> in_ring(bonds.size(), 0);
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> adj(atoms);  // (neighbour, bond)
  for (std::uint32_t b = 0; b < bonds.size(); ++b) {
    adj[bonds[b].first].push_back({bonds[b].second, b});
    adj[bonds[b].second].push_back({bonds[b].first, b});
  }
  std::vector<int> disc(atoms, -1), low(atoms, 0);
  std::vector<std::uint32_t> parent_bond(atoms, 0xFFFFFFFFu);
  std::vector<std::size_t> next(atoms, 0);
  int timer = 0;
  std::vector<std::uint32_t> stack;
  for (std::uint32_t root = 0; root < atoms; ++root) {
    if (disc[root] >= 0) continue;
    disc[root] = low[root] = timer++;
    stack.push_back(root);
    while (!stack.empty()) {
      const std::uint32_t v = stack.back();
      if (next[v] < adj[v].size()) {
        const auto [w, b] = adj[v][next[v]++];
        if (b == parent_bond[v]) continue;
        if (disc[w] >= 0) {
          low[v] = std::min(low[v], disc[w]);
        } else {
          parent_bond[w] = b;
          disc[w] = low[w] = timer++;
          stack.push_back(w);
        }
      } else {
        stack.pop_back();
        if (!stack.empty()) {
          const std::uint32_t u = stack.back();
          low[u] = std::min(low[u], low[v]);
          if (low[v] > disc[u]) in_ring[parent_bond[v]] = 0;  // bridge
          else in_ring[parent_bond[v]] = 1;
        }
      }
    }
  }
  return in_ring;
}

// The largest ligand docking will take, in heavy atoms. A drug-like ligand has 20 to 60 and a peptide ligand
// a hundred or two; past this the caller means a protein, not a ligand, and generating a conformation for it
// (distance geometry smooths a full distance matrix) would not finish.
constexpr std::uint32_t kMaxLigandHeavyAtoms = 300;

bool ligand_size_ok(const ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag) {
  const Species& s = state.species[species];
  std::uint32_t n = 0;
  for (std::uint32_t i = s.atom_begin; i < s.atom_begin + s.atom_count; ++i) n += state.atoms.element[i] != 1;
  if (n <= kMaxLigandHeavyAtoms) return true;
  diag.error(DiagCode::InvalidInput, "the ligand has " + std::to_string(n) +
                                         " heavy atoms; docking takes at most " +
                                         std::to_string(kMaxLigandHeavyAtoms))
      .on_species(species);
  return false;
}

// Builds the type, heavy-atom list, bonds, fragments and torsion tree of a ligand already in `state`.
bool finish_ligand(DockLigand& lig, Diagnostics& diag) {
  if (!ligand_size_ok(lig.state, lig.species, diag)) return false;
  const ChemicalSystemState& st = lig.state;
  const Species& s = st.species[lig.species];
  const std::uint32_t base = s.atom_begin, n = s.atom_count;
  if (n == 0) {
    diag.error(DiagCode::InvalidInput, "the ligand has no atoms");
    return false;
  }
  lig.smiles = s.canonical_smiles;
  lig.element.assign(n, 0);
  lig.coordinates.assign(n, Vec3{0, 0, 0});
  for (std::uint32_t i = 0; i < n; ++i) {
    lig.element[i] = st.atoms.element[base + i];
    lig.coordinates[i] = {st.atoms.x[base + i], st.atoms.y[base + i], st.atoms.z[base + i]};
    if (!(std::isfinite(lig.coordinates[i][0]) && std::isfinite(lig.coordinates[i][1]) &&
          std::isfinite(lig.coordinates[i][2]))) {
      diag.error(DiagCode::InvalidInput, "the ligand has a coordinate that is not a number");
      return false;
    }
  }
  lig.bonds.clear();
  std::vector<std::uint8_t> order;
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    lig.bonds.push_back({st.bonds.begin[b] - base, st.bonds.end[b] - base});
    order.push_back(st.bonds.order[b]);
  }

  std::vector<std::vector<std::uint32_t>> adj(n);
  for (const auto& [a, b] : lig.bonds) {
    adj[a].push_back(b);
    adj[b].push_back(a);
  }
  // Atom types: hydrogens carried explicitly plus the implicit count of the record.
  lig.type.assign(n, XsType::None);
  for (std::uint32_t i = 0; i < n; ++i) {
    const std::uint8_t z = lig.element[i];
    if (z == 1) continue;
    int hcount = st.atoms.implicit_h[base + i];
    int heavy_neighbours = 0;
    bool hetero = false;
    for (std::uint32_t j : adj[i]) {
      if (lig.element[j] == 1) ++hcount;
      else ++heavy_neighbours;
      if (lig.element[j] == 7 || lig.element[j] == 8) hetero = true;
    }
    lig.type[i] = xs_type(z, hetero, hcount, heavy_neighbours + hcount);
    if (lig.type[i] != XsType::None) lig.heavy.push_back(i);
  }
  if (lig.heavy.empty()) {
    diag.error(DiagCode::InvalidInput, "the ligand has no atoms the scoring function knows");
    return false;
  }

  // Rotatable bonds: single, acyclic, not an amide C-N, at least two heavy neighbours at each end.
  const std::vector<char> in_ring = ring_bonds(lig.bonds, n);
  auto heavy_degree = [&](std::uint32_t i) {
    std::size_t k = 0;
    for (std::uint32_t j : adj[i]) k += lig.element[j] != 1;
    return k;
  };
  auto is_amide_carbon = [&](std::uint32_t c) {
    if (lig.element[c] != 6) return false;
    for (std::size_t k = 0; k < lig.bonds.size(); ++k) {
      if (order[k] < 2) continue;
      const std::uint32_t a = lig.bonds[k].first, b = lig.bonds[k].second;
      const std::uint32_t other = a == c ? b : (b == c ? a : 0xFFFFFFFFu);
      if (other == 0xFFFFFFFFu) continue;
      if (lig.element[other] == 8 || lig.element[other] == 16) return true;
    }
    return false;
  };
  std::vector<std::size_t> rotatable;
  for (std::size_t k = 0; k < lig.bonds.size(); ++k) {
    const std::uint32_t a = lig.bonds[k].first, b = lig.bonds[k].second;
    if (order[k] != 1 || in_ring[k]) continue;
    if (lig.element[a] == 1 || lig.element[b] == 1) continue;
    if (heavy_degree(a) < 2 || heavy_degree(b) < 2) continue;
    const bool amide = (lig.element[a] == 7 && is_amide_carbon(b)) || (lig.element[b] == 7 && is_amide_carbon(a));
    if (amide) continue;
    rotatable.push_back(k);
  }

  // Rigid fragments: the connected components left when the rotatable bonds are cut.
  std::vector<char> cut(lig.bonds.size(), 0);
  for (std::size_t k : rotatable) cut[k] = 1;
  lig.fragment.assign(n, 0xFFFFFFFFu);
  std::uint32_t fragments = 0;
  std::vector<std::vector<std::pair<std::uint32_t, std::size_t>>> adj_bond(n);
  for (std::size_t k = 0; k < lig.bonds.size(); ++k) {
    adj_bond[lig.bonds[k].first].push_back({lig.bonds[k].second, k});
    adj_bond[lig.bonds[k].second].push_back({lig.bonds[k].first, k});
  }
  for (std::uint32_t i = 0; i < n; ++i) {
    if (lig.fragment[i] != 0xFFFFFFFFu) continue;
    std::deque<std::uint32_t> q{i};
    lig.fragment[i] = fragments;
    while (!q.empty()) {
      const std::uint32_t v = q.front();
      q.pop_front();
      for (const auto& [w, k] : adj_bond[v])
        if (!cut[k] && lig.fragment[w] == 0xFFFFFFFFu) {
          lig.fragment[w] = fragments;
          q.push_back(w);
        }
    }
    ++fragments;
  }

  // Torsion tree: root at the fragment with the most heavy atoms, then breadth-first over the cut bonds, so
  // a torsion's own axis is already in place when it is applied.
  std::vector<std::size_t> fragment_heavy(fragments, 0);
  for (std::uint32_t i = 0; i < n; ++i)
    if (lig.element[i] != 1) ++fragment_heavy[lig.fragment[i]];
  std::uint32_t root = 0;
  for (std::uint32_t f = 1; f < fragments; ++f)
    if (fragment_heavy[f] > fragment_heavy[root]) root = f;

  std::vector<char> seen(fragments, 0);
  seen[root] = 1;
  std::deque<std::uint32_t> fq{root};
  std::vector<std::pair<std::uint32_t, std::uint32_t>> order_bonds;  // (begin atom, end atom) in tree order
  while (!fq.empty()) {
    const std::uint32_t f = fq.front();
    fq.pop_front();
    for (std::size_t k : rotatable) {
      const std::uint32_t a = lig.bonds[k].first, b = lig.bonds[k].second;
      std::uint32_t begin = 0, end = 0;
      if (lig.fragment[a] == f && !seen[lig.fragment[b]]) { begin = a; end = b; }
      else if (lig.fragment[b] == f && !seen[lig.fragment[a]]) { begin = b; end = a; }
      else continue;
      seen[lig.fragment[end]] = 1;
      fq.push_back(lig.fragment[end]);
      order_bonds.push_back({begin, end});
    }
  }
  // `moves` of a torsion: every atom reachable from `end` without crossing the torsion's own bond.
  for (const auto& [begin, end] : order_bonds) {
    LigandTorsion t;
    t.begin = begin;
    t.end = end;
    std::vector<char> visited(n, 0);
    visited[begin] = 1;
    std::deque<std::uint32_t> q{end};
    visited[end] = 1;
    while (!q.empty()) {
      const std::uint32_t v = q.front();
      q.pop_front();
      t.moves.push_back(v);
      for (std::uint32_t w : adj[v])
        if (!visited[w]) {
          visited[w] = 1;
          q.push_back(w);
        }
    }
    std::sort(t.moves.begin(), t.moves.end());
    lig.torsions.push_back(std::move(t));
  }
  return true;
}

}  // namespace

std::optional<DockLigand> ligand_from_smiles(std::string_view smiles, Diagnostics& diag, std::uint64_t seed) {
  DockLigand lig;
  // The size is checked before the conformation is generated: distance-geometry embedding smooths a full
  // distance matrix, so a structure of thousands of atoms would take hours. A ligand is not that large.
  const std::int32_t parsed = parse_smiles(lig.state, smiles, diag);
  if (parsed < 0) return std::nullopt;
  if (!ligand_size_ok(lig.state, static_cast<std::uint32_t>(parsed), diag)) return std::nullopt;
  const std::int32_t with_h =
      add_explicit_hydrogens(lig.state, static_cast<std::uint32_t>(parsed), diag);
  if (with_h < 0) return std::nullopt;
  EmbedOptions eo;
  eo.seed = seed;
  const EmbedResult r = embed_species(lig.state, static_cast<std::uint32_t>(with_h), diag, eo);
  if (!r.ok) return std::nullopt;
  lig.species = static_cast<std::uint32_t>(with_h);
  if (!finish_ligand(lig, diag)) return std::nullopt;
  return lig;
}

std::optional<DockLigand> ligand_from_molfile(std::string_view text, Diagnostics& diag, std::uint64_t seed) {
  (void)seed;
  DockLigand lig;
  ReadOptions ro;
  ro.keep_hydrogens = true;
  ro.max_records = 1;
  const std::vector<MolRecord> records = read_molecules(lig.state, text, diag, ro);
  if (records.empty() || records[0].species < 0) {
    if (!diag.has_errors()) diag.error(DiagCode::MalformedRecord, "the ligand record could not be read");
    return std::nullopt;
  }
  if (records[0].dimensions != 3) {
    diag.error(DiagCode::InvalidInput, "the ligand record has no 3D coordinates");
    return std::nullopt;
  }
  lig.species = static_cast<std::uint32_t>(records[0].species);
  if (!finish_ligand(lig, diag)) return std::nullopt;
  return lig;
}

// ---- grid maps ----------------------------------------------------------------------------------

namespace {

constexpr int kTypeCount = static_cast<int>(XsType::Met_D) + 1;

struct GridMaps {
  DockBox box;
  double spacing = 0.375;
  Vec3 origin{0, 0, 0};
  std::array<int, 3> n{1, 1, 1};
  std::array<int, kTypeCount> slot{};  // XsType -> map index, -1 when the ligand has no such atom
  std::vector<std::vector<double>> maps;

  std::size_t at(int ix, int iy, int iz) const {
    return (static_cast<std::size_t>(ix) * n[1] + iy) * n[2] + iz;
  }

  // Trilinear interpolation, with the gradient of the interpolated surface. Outside the grid the nearest
  // face value is used and the gradient is zero in that direction: the box penalty does the pushing there.
  double value(int map, const Vec3& p, Vec3& gradient) const {
    const std::vector<double>& m = maps[static_cast<std::size_t>(map)];
    int i0[3];
    double f[3];
    gradient = {0, 0, 0};
    bool clamped[3] = {false, false, false};
    for (int k = 0; k < 3; ++k) {
      const double u = (p[k] - origin[k]) / spacing;
      const double hi = static_cast<double>(n[k] - 1);
      if (!(u > 0.0)) {
        i0[k] = 0;
        f[k] = 0.0;
        clamped[k] = true;
      } else if (u >= hi) {
        i0[k] = std::max(n[k] - 2, 0);
        f[k] = n[k] >= 2 ? 1.0 : 0.0;
        clamped[k] = true;
      } else {
        i0[k] = static_cast<int>(u);
        if (i0[k] > n[k] - 2) i0[k] = std::max(n[k] - 2, 0);
        f[k] = u - i0[k];
      }
    }
    const int i1[3] = {std::min(i0[0] + 1, n[0] - 1), std::min(i0[1] + 1, n[1] - 1), std::min(i0[2] + 1, n[2] - 1)};
    const double c000 = m[at(i0[0], i0[1], i0[2])], c001 = m[at(i0[0], i0[1], i1[2])];
    const double c010 = m[at(i0[0], i1[1], i0[2])], c011 = m[at(i0[0], i1[1], i1[2])];
    const double c100 = m[at(i1[0], i0[1], i0[2])], c101 = m[at(i1[0], i0[1], i1[2])];
    const double c110 = m[at(i1[0], i1[1], i0[2])], c111 = m[at(i1[0], i1[1], i1[2])];
    const double x = f[0], y = f[1], z = f[2];
    const double c00 = c000 + (c100 - c000) * x, c01 = c001 + (c101 - c001) * x;
    const double c10 = c010 + (c110 - c010) * x, c11 = c011 + (c111 - c011) * x;
    const double c0 = c00 + (c10 - c00) * y, c1 = c01 + (c11 - c01) * y;
    if (!clamped[0]) {
      const double d00 = c100 - c000, d01 = c101 - c001, d10 = c110 - c010, d11 = c111 - c011;
      const double d0 = d00 + (d10 - d00) * y, d1 = d01 + (d11 - d01) * y;
      gradient[0] = (d0 + (d1 - d0) * z) / spacing;
    }
    if (!clamped[1]) {
      const double d0 = c10 - c00, d1 = c11 - c01;
      gradient[1] = (d0 + (d1 - d0) * z) / spacing;
    }
    if (!clamped[2]) gradient[2] = (c1 - c0) / spacing;
    return c0 + (c1 - c0) * z;
  }
};

// One map per atom type the ligand actually has: the receptor's contribution to an atom of that type at every
// grid point. The grid reaches kGridMargin beyond the search box so that an atom pushed out by the penalty
// still reads a real energy.
GridMaps build_grid(const DockReceptor& receptor, const DockLigand& ligand, const DockOptions& options,
                    const CellList& cells, const DockLigand* also) {
  GridMaps g;
  g.box = options.box;
  g.spacing = std::max(options.grid_spacing, 0.05);
  g.slot.fill(-1);
  int used = 0;
  // Every type either molecule has. `also` is the reference pose, which is minimised on the same grid: it
  // can carry a type the docked ligand does not, since one comes from a SMILES with explicit hydrogens and
  // the other from a file with implicit ones, and a nitrogen's type turns on how many it has.
  auto want = [&](const DockLigand& l) {
    for (std::uint32_t i : l.heavy) {
      const int t = static_cast<int>(l.type[i]);
      if (g.slot[t] < 0) g.slot[t] = used++;
    }
  };
  want(ligand);
  if (also) want(*also);
  const Vec3 size{options.box.sx + 2 * kGridMargin, options.box.sy + 2 * kGridMargin,
                  options.box.sz + 2 * kGridMargin};
  g.origin = {options.box.cx - 0.5 * size[0], options.box.cy - 0.5 * size[1], options.box.cz - 0.5 * size[2]};
  // A box of ordinary size at 0.375 A is under a million points, a few megabytes a map. A much larger box
  // would ask for gigabytes, so the spacing is widened until the maps fit instead: a coarse map of a huge
  // box is a poor search, but it is a search, and it cannot exhaust the machine.
  auto axis_points = [&](double edge) { return std::floor(edge / g.spacing) + 2.0; };
  while (axis_points(size[0]) * axis_points(size[1]) * axis_points(size[2]) > static_cast<double>(kMaxGridPoints))
    g.spacing *= 1.5;
  for (int k = 0; k < 3; ++k)
    g.n[k] = static_cast<int>(std::min(std::max(axis_points(size[k]), 2.0), 4096.0));
  const std::size_t total = static_cast<std::size_t>(g.n[0]) * g.n[1] * g.n[2];
  g.maps.assign(static_cast<std::size_t>(used), std::vector<double>(total, 0.0));

  // Which type each map holds.
  std::vector<XsType> map_type(static_cast<std::size_t>(used), XsType::None);
  for (int t = 0; t < kTypeCount; ++t)
    if (g.slot[t] >= 0) map_type[static_cast<std::size_t>(g.slot[t])] = static_cast<XsType>(t);

  const ScoringWeights w = options.weights;
  // One index per x slice: each writes only its own part of every map.
  parallel_for(static_cast<std::size_t>(g.n[0]), [&](std::size_t ix) {
    for (int iy = 0; iy < g.n[1]; ++iy)
      for (int iz = 0; iz < g.n[2]; ++iz) {
        const Vec3 p{g.origin[0] + static_cast<double>(ix) * g.spacing, g.origin[1] + iy * g.spacing,
                     g.origin[2] + iz * g.spacing};
        const std::size_t cell = g.at(static_cast<int>(ix), iy, iz);
        cells.for_each_near(p, w.cutoff, [&](std::uint32_t j, double r2) {
          const double r = std::sqrt(r2);
          const XsType rt = receptor.atoms[j].type;
          for (std::size_t m = 0; m < map_type.size(); ++m)
            g.maps[m][cell] += pair_energy(map_type[m], rt, r, w);
        });
      }
  });
  return g;
}

// ---- the pose model -----------------------------------------------------------------------------

// Everything a score evaluation needs: the ligand in its body frame, the intramolecular pairs, the receptor
// grid, and the box the pose must stay in.
struct Model {
  const DockLigand* lig = nullptr;
  // Either a grid of precomputed receptor terms (the search) or the receptor itself with a cell list over
  // it (an exact sum, for scoring and for checking the gradient).
  const GridMaps* grid = nullptr;
  const DockReceptor* receptor = nullptr;
  const CellList* cells = nullptr;
  ScoringWeights w;
  ScoringWeights intra_w;  // `w` with the attraction inside the ligand scaled by internal_attraction
  // > 0 while hunting: squashes repulsion, per ligand atom for the receptor and per pair inside the ligand.
  double inter_cap = 0.0;
  double intra_cap = 0.0;
  std::vector<Vec3> base;        // input coordinates minus the rotation centre
  Vec3 centre{0, 0, 0};
  std::vector<std::pair<std::uint32_t, std::uint32_t>> intra;  // heavy pairs whose distance can change
  std::size_t torsions = 0;
  mutable std::size_t evaluations = 0;

  std::size_t variables() const { return 6 + torsions; }
};

// `w` with every attractive term inside the ligand scaled by internal_attraction (repulsion untouched).
ScoringWeights internal_weights(const ScoringWeights& w) {
  ScoringWeights i = w;
  const double s = std::max(w.internal_attraction, 0.0);
  i.gauss1 *= s;
  i.gauss2 *= s;
  i.hydrophobic *= s;
  i.hbond *= s;
  return i;
}

Model build_model(const DockLigand& lig, const GridMaps* grid, const ScoringWeights& w) {
  Model m;
  m.lig = &lig;
  m.grid = grid;
  m.w = w;
  m.intra_w = internal_weights(w);
  m.torsions = lig.torsions.size();
  Vec3 c{0, 0, 0};
  for (std::uint32_t i : lig.heavy) c = c + lig.coordinates[i];
  m.centre = (1.0 / static_cast<double>(lig.heavy.size())) * c;
  m.base.resize(lig.coordinates.size());
  for (std::size_t i = 0; i < lig.coordinates.size(); ++i) m.base[i] = lig.coordinates[i] - m.centre;

  // Intramolecular pairs: more than three bonds apart and in different rigid fragments, so that the terms
  // only see distances the torsions can actually change.
  const std::size_t n = lig.coordinates.size();
  std::vector<std::vector<std::uint32_t>> adj(n);
  for (const auto& [a, b] : lig.bonds) {
    adj[a].push_back(b);
    adj[b].push_back(a);
  }
  std::vector<int> dist(n);
  for (std::size_t k = 0; k < lig.heavy.size(); ++k) {
    const std::uint32_t i = lig.heavy[k];
    std::fill(dist.begin(), dist.end(), -1);
    dist[i] = 0;
    std::deque<std::uint32_t> q{i};
    while (!q.empty()) {
      const std::uint32_t v = q.front();
      q.pop_front();
      if (dist[v] >= 4) continue;  // further pairs are all included anyway
      for (std::uint32_t u : adj[v])
        if (dist[u] < 0) {
          dist[u] = dist[v] + 1;
          q.push_back(u);
        }
    }
    for (std::size_t l = k + 1; l < lig.heavy.size(); ++l) {
      const std::uint32_t j = lig.heavy[l];
      if (lig.fragment[i] == lig.fragment[j]) continue;
      if (dist[j] >= 0 && dist[j] < 4) continue;
      m.intra.push_back({i, j});
    }
  }
  return m;
}

// A pose: translation, orientation and torsion angles.
struct Pose {
  Vec3 translation{0, 0, 0};
  Quat orientation = quat_identity();
  std::vector<double> torsions;
};

// Ligand coordinates of a pose, with the rotation vector `r` applied on top of the stored orientation.
void pose_coordinates(const Model& m, const Pose& pose, const Vec3& r, std::vector<Vec3>& out, Mat3& rot,
                      std::vector<Vec3>& body) {
  const DockLigand& lig = *m.lig;
  body = m.base;
  for (std::size_t j = 0; j < lig.torsions.size(); ++j) {
    const LigandTorsion& t = lig.torsions[j];
    const Vec3 axis = normalized(body[t.end] - body[t.begin]);
    const Vec3 pivot = body[t.end];
    const double angle = pose.torsions[j];
    for (std::uint32_t i : t.moves) body[i] = pivot + rotate_about(body[i] - pivot, axis, angle);
  }
  const Mat3 q = matrix_from_quat(pose.orientation);
  for (Vec3& v : body) v = mat_apply(q, v);
  rot = matrix_from_quat(quat_from_rotation_vector(r));
  out.resize(body.size());
  for (std::size_t i = 0; i < body.size(); ++i) out[i] = pose.translation + mat_apply(rot, body[i]);
}

// Energy of a pose and its Cartesian gradient on every ligand atom.
struct Energy {
  double intermolecular = 0, internal = 0, penalty = 0;
  double total() const { return intermolecular + internal + penalty; }
};

Energy evaluate(const Model& m, const std::vector<Vec3>& x, std::vector<Vec3>* gradient) {
  const DockLigand& lig = *m.lig;
  Energy e;
  if (gradient) gradient->assign(x.size(), Vec3{0, 0, 0});
  ++m.evaluations;

  if (m.grid) {
    // Receptor: the precomputed map of the atom's type.
    const GridMaps& g = *m.grid;
    for (std::uint32_t i : lig.heavy) {
      const int slot = g.slot[static_cast<int>(lig.type[i])];
      if (slot < 0) continue;  // no map for this type: it was not in the set the grid was built for
      Vec3 gr{0, 0, 0};
      double squashed = 1.0;
      e.intermolecular += squash(g.value(slot, x[i], gr), m.inter_cap, &squashed);
      if (gradient) (*gradient)[i] = (*gradient)[i] + squashed * gr;
    }
    // The box the pose must stay in.
    const double lo[3] = {g.box.cx - 0.5 * g.box.sx, g.box.cy - 0.5 * g.box.sy, g.box.cz - 0.5 * g.box.sz};
    const double hi[3] = {g.box.cx + 0.5 * g.box.sx, g.box.cy + 0.5 * g.box.sy, g.box.cz + 0.5 * g.box.sz};
    for (std::uint32_t i : lig.heavy)
      for (int k = 0; k < 3; ++k) {
        double d = 0;
        if (x[i][k] < lo[k]) d = x[i][k] - lo[k];
        else if (x[i][k] > hi[k]) d = x[i][k] - hi[k];
        if (d == 0.0) continue;
        e.penalty += kBoxPenalty * d * d;
        if (gradient) (*gradient)[i][k] += 2.0 * kBoxPenalty * d;
      }
  } else if (m.receptor && m.cells) {
    // The same sum, exactly: every receptor atom within the cutoff. No box, so no penalty.
    for (std::uint32_t i : lig.heavy)
      m.cells->for_each_near(x[i], m.w.cutoff, [&](std::uint32_t j, double r2) {
        const double r = std::sqrt(r2);
        double de = 0;
        e.intermolecular += pair_energy(lig.type[i], m.receptor->atoms[j].type, r, m.w, gradient ? &de : nullptr);
        if (gradient && r > 1e-9) {
          const Vec3 d{m.receptor->atoms[j].x - x[i][0], m.receptor->atoms[j].y - x[i][1],
                       m.receptor->atoms[j].z - x[i][2]};
          (*gradient)[i] = (*gradient)[i] - (de / r) * d;
        }
      });
  }
  // The ligand's own pairs.
  for (const auto& [i, j] : m.intra) {
    const Vec3 d = x[j] - x[i];
    const double r2 = dot(d, d);
    if (r2 >= m.w.cutoff * m.w.cutoff) continue;
    const double r = std::sqrt(r2);
    double de = 0, squashed = 1.0;
    e.internal += squash(pair_energy(lig.type[i], lig.type[j], r, m.intra_w, gradient ? &de : nullptr),
                         m.intra_cap, &squashed);
    de *= squashed;
    if (gradient && r > 1e-9) {
      const Vec3 f = (de / r) * d;
      (*gradient)[i] = (*gradient)[i] - f;
      (*gradient)[j] = (*gradient)[j] + f;
    }
  }
  return e;
}

// The objective handed to minimize_lbfgs: x = [translation(3), rotation vector(3), torsions].
class PoseObjective {
 public:
  PoseObjective(const Model& m, const Pose& base) : m_(m), base_(base) {}

  double operator()(const std::vector<double>& v, std::vector<double>& grad) {
    Pose p = base_;
    p.translation = {v[0], v[1], v[2]};
    const Vec3 r{v[3], v[4], v[5]};
    for (std::size_t j = 0; j < p.torsions.size(); ++j) p.torsions[j] = v[6 + j];
    Mat3 rot{};
    pose_coordinates(m_, p, r, world_, rot, body_);
    const Energy e = evaluate(m_, world_, &cart_);
    grad.assign(v.size(), 0.0);
    Vec3 force{0, 0, 0}, torque{0, 0, 0};
    for (std::uint32_t i : m_.lig->heavy) {
      force = force + cart_[i];
      torque = torque + cross(body_[i], mat_apply_t(rot, cart_[i]));
    }
    grad[0] = force[0];
    grad[1] = force[1];
    grad[2] = force[2];
    const Vec3 gr = right_jacobian_transpose(r, torque);
    grad[3] = gr[0];
    grad[4] = gr[1];
    grad[5] = gr[2];
    for (std::size_t j = 0; j < m_.lig->torsions.size(); ++j) {
      const LigandTorsion& t = m_.lig->torsions[j];
      const Vec3 axis = normalized(world_[t.end] - world_[t.begin]);
      const Vec3 pivot = world_[t.end];
      double d = 0;
      for (std::uint32_t i : t.moves)
        if (m_.lig->type[i] != XsType::None) d += dot(cart_[i], cross(axis, world_[i] - pivot));
      grad[6 + j] = d;
    }
    last_ = e;
    return e.total();
  }

  const Energy& last() const { return last_; }
  const std::vector<Vec3>& world() const { return world_; }

 private:
  const Model& m_;
  Pose base_;
  std::vector<Vec3> world_, body_, cart_;
  Energy last_;
};

// Local minimisation. The rotation vector starts at zero and is folded into the pose's quaternion
// afterwards, so the next minimisation starts from a fresh, singularity-free parameterisation.
Energy minimize_pose(const Model& m, Pose& pose, const LbfgsOptions& options, std::vector<Vec3>& coordinates) {
  std::vector<double> v(m.variables(), 0.0);
  v[0] = pose.translation[0];
  v[1] = pose.translation[1];
  v[2] = pose.translation[2];
  for (std::size_t j = 0; j < pose.torsions.size(); ++j) v[6 + j] = pose.torsions[j];
  PoseObjective obj(m, pose);
  Objective f = [&obj](const std::vector<double>& x, std::vector<double>& g) { return obj(x, g); };
  minimize_lbfgs(f, v, options);
  // One last evaluation at the returned point, so energy and coordinates match it exactly.
  std::vector<double> g;
  obj(v, g);
  pose.translation = {v[0], v[1], v[2]};
  pose.orientation = quat_normalize(quat_multiply(quat_from_rotation_vector({v[3], v[4], v[5]}), pose.orientation));
  for (std::size_t j = 0; j < pose.torsions.size(); ++j) pose.torsions[j] = v[6 + j];
  coordinates = obj.world();
  return obj.last();
}

// RMSD between two poses of the same molecule, over the symmetry group given as permutations of the
// heavy-atom list (the identity is always the first one), so that two poses related by a flip of a phenyl
// ring or a carboxylate count as the same pose and are not both reported.
double heavy_rmsd(const DockLigand& lig, const std::vector<std::vector<std::uint32_t>>& symmetry,
                  const std::vector<Vec3>& a, const std::vector<Vec3>& b) {
  const std::size_t n = lig.heavy.size();
  double best = -1;
  for (const std::vector<std::uint32_t>& perm : symmetry) {
    double s = 0;
    for (std::size_t k = 0; k < n; ++k) {
      const Vec3 d = a[lig.heavy[k]] - b[lig.heavy[perm[k]]];
      s += dot(d, d);
      if (best >= 0 && s >= best) break;
    }
    if (best < 0 || s < best) best = s;
  }
  return std::sqrt(std::max(best, 0.0) / static_cast<double>(n));
}

struct Candidate {
  double rank = 0;  // total energy over the rotatable-bond factor: what the poses are ordered by
  double intermolecular = 0, internal = 0;
  std::vector<Vec3> coordinates;
};

}  // namespace

// ---- symmetry-corrected RMSD --------------------------------------------------------------------

namespace {

// Refined invariants (element and the multiset of neighbours' invariants) over the heavy-atom graph: atoms
// that cannot map onto each other get different values, which prunes the search below.
std::vector<std::uint64_t> heavy_invariants(const std::vector<std::uint8_t>& element,
                                           const std::vector<std::vector<std::uint32_t>>& adj) {
  const std::size_t n = element.size();
  std::vector<std::uint64_t> inv(n);
  for (std::size_t i = 0; i < n; ++i) inv[i] = element[i] * 131u + adj[i].size();
  std::vector<std::uint64_t> next(n);
  for (int round = 0; round < 8; ++round) {
    for (std::size_t i = 0; i < n; ++i) {
      std::vector<std::uint64_t> ns;
      ns.reserve(adj[i].size());
      for (std::uint32_t j : adj[i]) ns.push_back(inv[j]);
      std::sort(ns.begin(), ns.end());
      std::uint64_t h = inv[i] * 0x9E3779B97F4A7C15ull;
      for (std::uint64_t v : ns) h = (h ^ v) * 0x100000001B3ull;
      next[i] = h;
    }
    inv = next;
  }
  return inv;
}

struct HeavyGraph {
  std::vector<std::uint32_t> atoms;  // ligand-local indices, in graph order
  std::vector<std::uint8_t> element;
  std::vector<std::vector<std::uint32_t>> adj;  // indices into `atoms`
  std::vector<std::vector<char>> connected;
  std::vector<std::uint64_t> invariant;
};

HeavyGraph heavy_graph(const DockLigand& lig) {
  HeavyGraph g;
  g.atoms = lig.heavy;
  std::vector<std::uint32_t> local(lig.coordinates.size(), 0xFFFFFFFFu);
  for (std::size_t k = 0; k < g.atoms.size(); ++k) local[g.atoms[k]] = static_cast<std::uint32_t>(k);
  const std::size_t n = g.atoms.size();
  g.element.resize(n);
  for (std::size_t k = 0; k < n; ++k) g.element[k] = lig.element[g.atoms[k]];
  g.adj.assign(n, {});
  g.connected.assign(n, std::vector<char>(n, 0));
  for (const auto& [a, b] : lig.bonds) {
    if (a >= local.size() || b >= local.size()) continue;
    const std::uint32_t ia = local[a], ib = local[b];
    if (ia == 0xFFFFFFFFu || ib == 0xFFFFFFFFu) continue;
    g.adj[ia].push_back(ib);
    g.adj[ib].push_back(ia);
    g.connected[ia][ib] = g.connected[ib][ia] = 1;
  }
  g.invariant = heavy_invariants(g.element, g.adj);
  return g;
}

// The automorphisms of the heavy-atom graph, as permutations of the heavy-atom list; the identity comes
// first. Used to tell two poses of the same molecule apart when they differ only by a symmetry of its own
// graph. Capped, so a molecule with many equivalent groups cannot blow the list up.
std::vector<std::vector<std::uint32_t>> heavy_automorphisms(const DockLigand& lig, std::size_t cap) {
  const HeavyGraph g = heavy_graph(lig);
  const std::size_t n = g.atoms.size();
  std::vector<std::vector<std::uint32_t>> out;
  std::vector<std::uint32_t> identity(n);
  for (std::size_t k = 0; k < n; ++k) identity[k] = static_cast<std::uint32_t>(k);
  out.push_back(identity);
  if (n == 0) return out;

  std::vector<std::uint32_t> map(n, 0xFFFFFFFFu);
  std::vector<char> used(n, 0);
  // Depth-first over the atoms in index order; a candidate must agree with every assignment made so far on
  // adjacency, which makes the search cheap for the molecules docking sees.
  struct Search {
    const HeavyGraph& g;
    std::vector<std::uint32_t>& map;
    std::vector<char>& used;
    std::vector<std::vector<std::uint32_t>>& out;
    std::size_t cap;

    void go(std::size_t depth) {
      if (out.size() >= cap) return;
      const std::size_t n = g.atoms.size();
      if (depth == n) {
        bool identity_map = true;
        for (std::size_t k = 0; k < n && identity_map; ++k) identity_map = map[k] == k;
        if (!identity_map) out.push_back(map);
        return;
      }
      for (std::uint32_t c = 0; c < n; ++c) {
        if (used[c] || g.element[depth] != g.element[c] || g.invariant[depth] != g.invariant[c]) continue;
        bool ok = true;
        for (std::size_t d = 0; d < depth && ok; ++d) ok = g.connected[depth][d] == g.connected[c][map[d]];
        if (!ok) continue;
        map[depth] = c;
        used[c] = 1;
        go(depth + 1);
        used[c] = 0;
        map[depth] = 0xFFFFFFFFu;
        if (out.size() >= cap) return;
      }
    }
  };
  Search s{g, map, used, out, cap};
  s.go(0);
  return out;
}

}  // namespace

double symmetry_rmsd(const DockLigand& a, const PointCloud& a_pose, const DockLigand& b, const PointCloud& b_pose,
                     std::size_t* matchings) {
  if (matchings) *matchings = 0;
  const HeavyGraph ga = heavy_graph(a), gb = heavy_graph(b);
  const std::size_t n = ga.atoms.size();
  if (n == 0 || n != gb.atoms.size()) return -1;
  if (a_pose.size() < a.coordinates.size() || b_pose.size() < b.coordinates.size()) return -1;

  // Visit a's atoms so that every atom after the first touches one already placed: the adjacency test then
  // prunes as early as possible.
  std::vector<std::uint32_t> order;
  std::vector<char> queued(n, 0);
  for (std::size_t s = 0; s < n; ++s) {
    if (queued[s]) continue;
    std::deque<std::uint32_t> q{static_cast<std::uint32_t>(s)};
    queued[s] = 1;
    while (!q.empty()) {
      const std::uint32_t v = q.front();
      q.pop_front();
      order.push_back(v);
      for (std::uint32_t w : ga.adj[v])
        if (!queued[w]) {
          queued[w] = 1;
          q.push_back(w);
        }
    }
  }

  std::vector<Vec3> pa(n), pb(n);
  for (std::size_t k = 0; k < n; ++k) {
    const auto& p = a_pose[ga.atoms[k]];
    const auto& q = b_pose[gb.atoms[k]];
    pa[k] = {p[0], p[1], p[2]};
    pb[k] = {q[0], q[1], q[2]};
  }

  // Depth-first assignment carrying the partial sum of squares, so a branch already worse than the best
  // complete matching is abandoned. The recursion is as deep as the molecule has heavy atoms.
  struct Search {
    const HeavyGraph& ga;
    const HeavyGraph& gb;
    const std::vector<std::uint32_t>& order;
    const std::vector<Vec3>& pa;
    const std::vector<Vec3>& pb;
    std::vector<std::uint32_t> map;
    std::vector<char> used;
    double best = -1;
    std::size_t tried = 0;

    void go(std::size_t depth, double partial) {
      if (tried >= kMaxMatchings) return;
      if (depth == order.size()) {
        if (best < 0 || partial < best) best = partial;
        ++tried;
        return;
      }
      const std::uint32_t v = order[depth];
      for (std::uint32_t c = 0; c < gb.element.size(); ++c) {
        if (used[c] || ga.element[v] != gb.element[c] || ga.invariant[v] != gb.invariant[c]) continue;
        bool ok = true;
        for (std::size_t d = 0; d < depth && ok; ++d) {
          const std::uint32_t u = order[d];
          ok = ga.connected[v][u] == gb.connected[c][map[u]];
        }
        if (!ok) continue;
        const Vec3 diff = pa[v] - pb[c];
        const double next = partial + dot(diff, diff);
        if (best >= 0 && next >= best) continue;
        map[v] = c;
        used[c] = 1;
        go(depth + 1, next);
        used[c] = 0;
        map[v] = 0xFFFFFFFFu;
        if (tried >= kMaxMatchings) return;
      }
    }
  };
  Search s{ga, gb, order, pa, pb, std::vector<std::uint32_t>(n, 0xFFFFFFFFu), std::vector<char>(n, 0), -1, 0};
  s.go(0, 0.0);
  if (matchings) *matchings = s.tried;
  if (s.best < 0) return -1;
  return std::sqrt(s.best / static_cast<double>(n));
}

// ---- scoring a pose as it stands ----------------------------------------------------------------

namespace {

// Exact sum over the receptor atoms within the cutoff (no grid), so the numbers reported do not carry the
// interpolation error of the maps.
DockedPose exact_score(const DockReceptor& receptor, const CellList& cells, const DockLigand& lig,
                       const std::vector<Vec3>& x, const Model* model, const ScoringWeights& w) {
  DockedPose p;
  for (std::uint32_t i : lig.heavy) {
    cells.for_each_near(x[i], w.cutoff, [&](std::uint32_t j, double r2) {
      const ScoringTerms t = pair_terms(lig.type[i], receptor.atoms[j].type, std::sqrt(r2), w);
      p.terms.gauss1 += t.gauss1;
      p.terms.gauss2 += t.gauss2;
      p.terms.repulsion += t.repulsion;
      p.terms.hydrophobic += t.hydrophobic;
      p.terms.hbond += t.hbond;
    });
  }
  p.intermolecular = p.terms.weighted(w);
  if (model)
    for (const auto& [i, j] : model->intra) {
      const Vec3 d = x[j] - x[i];
      p.internal += pair_energy(lig.type[i], lig.type[j], length(d), model->intra_w);
    }
  p.affinity = p.intermolecular / (1.0 + w.rotatable * static_cast<double>(lig.torsions.size()));
  p.coordinates.resize(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) p.coordinates[i] = {x[i][0], x[i][1], x[i][2]};
  return p;
}

}  // namespace

std::optional<DockedPose> score_pose(const DockReceptor& receptor, const DockLigand& ligand, const PointCloud& pose,
                                     const ScoringWeights& w) {
  if (pose.size() != ligand.coordinates.size()) return std::nullopt;
  std::vector<Vec3> x(pose.size());
  for (std::size_t i = 0; i < pose.size(); ++i) x[i] = {pose[i][0], pose[i][1], pose[i][2]};
  std::vector<Vec3> rp(receptor.atoms.size());
  for (std::size_t i = 0; i < receptor.atoms.size(); ++i)
    rp[i] = {receptor.atoms[i].x, receptor.atoms[i].y, receptor.atoms[i].z};
  CellList cells;
  cells.build(rp, w.cutoff);
  const Model m = build_model(ligand, nullptr, w);  // only its intramolecular pair list is used here
  return exact_score(receptor, cells, ligand, x, &m, w);
}

std::size_t pose_variable_count(const DockLigand& ligand) { return 6 + ligand.torsions.size(); }

double pose_variable_energy(const DockReceptor& receptor, const DockLigand& ligand,
                            const std::vector<double>& variables, const ScoringWeights& w,
                            std::vector<double>* gradient, PointCloud* coordinates) {
  if (variables.size() != pose_variable_count(ligand)) {
    if (gradient) gradient->assign(variables.size(), 0.0);
    return 0.0;
  }
  std::vector<Vec3> rp(receptor.atoms.size());
  for (std::size_t i = 0; i < receptor.atoms.size(); ++i)
    rp[i] = {receptor.atoms[i].x, receptor.atoms[i].y, receptor.atoms[i].z};
  CellList cells;
  cells.build(rp, w.cutoff);
  Model m = build_model(ligand, nullptr, w);
  m.receptor = &receptor;
  m.cells = &cells;
  Pose base;
  base.torsions.assign(ligand.torsions.size(), 0.0);
  PoseObjective obj(m, base);
  std::vector<double> g;
  const double e = obj(variables, g);
  if (gradient) *gradient = g;
  if (coordinates) {
    coordinates->resize(obj.world().size());
    for (std::size_t i = 0; i < obj.world().size(); ++i)
      (*coordinates)[i] = {obj.world()[i][0], obj.world()[i][1], obj.world()[i][2]};
  }
  return e;
}

// ---- the search ---------------------------------------------------------------------------------

std::optional<DockResult> dock(const DockReceptor& receptor, const DockLigand& ligand, const DockOptions& options,
                               Diagnostics& diag, const DockLigand* reference) {
  const auto started = std::chrono::steady_clock::now();
  if (receptor.atoms.empty()) {
    diag.error(DiagCode::InvalidInput, "the receptor has no atoms");
    return std::nullopt;
  }
  if (!options.box.valid()) {
    diag.error(DiagCode::InvalidInput, "the search box has no size");
    return std::nullopt;
  }
  if (options.runs <= 0 || options.modes <= 0) {
    diag.error(DiagCode::InvalidInput, "runs and modes must be positive");
    return std::nullopt;
  }
  std::vector<Vec3> rp(receptor.atoms.size());
  for (std::size_t i = 0; i < receptor.atoms.size(); ++i)
    rp[i] = {receptor.atoms[i].x, receptor.atoms[i].y, receptor.atoms[i].z};
  CellList cells;
  cells.build(rp, options.weights.cutoff);

  // One grid of the function as it is; the hunting model reads the same maps and squashes what it reads.
  const GridMaps grid = build_grid(receptor, ligand, options, cells, reference);
  const Model model = build_model(ligand, &grid, options.weights);
  Model hunt_model = model;
  hunt_model.inter_cap = kHuntCapIntermolecular;
  hunt_model.intra_cap = kHuntCapInternal;
  // Poses that differ only by a symmetry of the ligand's own graph are the same pose, so clustering uses
  // these permutations rather than the atom order.
  const std::vector<std::vector<std::uint32_t>> symmetry = heavy_automorphisms(ligand, 128);
  const int steps = options.steps > 0 ? options.steps
                                      : 200 + 100 * static_cast<int>(std::min<std::size_t>(ligand.torsions.size(), 20));

  DockResult result;
  result.heavy_atoms = ligand.heavy.size();
  result.rotatable_bonds = ligand.torsions.size();
  result.receptor_atoms = receptor.atoms.size();
  result.box = options.box;
  result.grid_spacing = grid.spacing;
  result.runs = options.runs;
  result.steps = options.refine_only ? 0 : steps;

  const double factor = 1.0 + options.weights.rotatable * static_cast<double>(ligand.torsions.size());
  auto rank_of = [&](const Energy& e) { return (e.intermolecular + e.internal + e.penalty) / factor; };

  // The input pose, minimised where it stands. Also the whole job when refine_only is set.
  auto refine_input = [&](std::vector<Vec3>& coordinates) {
    Pose p;
    p.torsions.assign(ligand.torsions.size(), 0.0);
    p.translation = model.centre;
    const Energy e = minimize_pose(model, p, options.local, coordinates);
    return e;
  };

  std::vector<std::vector<Candidate>> per_run(static_cast<std::size_t>(options.runs));
  std::size_t evaluations = 0;

  if (options.refine_only) {
    std::vector<Vec3> x;
    const Energy e = refine_input(x);
    per_run[0].push_back({rank_of(e), e.intermolecular, e.internal, x});
    evaluations = model.evaluations;
  } else {
    // Independent Monte Carlo runs. Each writes only its own slot and draws from its own generator, so the
    // result is the same whatever the thread count.
    std::vector<std::size_t> run_evaluations(static_cast<std::size_t>(options.runs), 0);
    parallel_for(static_cast<std::size_t>(options.runs), [&](std::size_t run) {
      Rng rng(options.seed * 0x9E3779B97F4A7C15ull + run);
      Model local_model = model;  // its own evaluation counter
      local_model.evaluations = 0;
      Model local_hunt = hunt_model;
      local_hunt.evaluations = 0;
      std::vector<Candidate>& out = per_run[run];
      std::vector<Vec3> x, refined_x;
      Pose current;
      current.torsions.assign(ligand.torsions.size(), 0.0);
      auto randomize = [&](Pose& p) {
        p.translation = {rng.uniform(options.box.cx - 0.5 * options.box.sx, options.box.cx + 0.5 * options.box.sx),
                         rng.uniform(options.box.cy - 0.5 * options.box.sy, options.box.cy + 0.5 * options.box.sy),
                         rng.uniform(options.box.cz - 0.5 * options.box.sz, options.box.cz + 0.5 * options.box.sz)};
        p.orientation = quat_normalize(rng.orientation());
        for (double& t : p.torsions) t = rng.uniform(-kPi, kPi);
      };
      // A Monte Carlo step moves the centre inside a ball of 2 A, or turns the ligand by a rotation vector
      // inside a ball of 2 A / radius of gyration, or sets one torsion at random. Scaling the rotation by
      // the ligand's size keeps the step a comparable displacement for a small and a large molecule: an
      // unscaled 2 rad kick throws a 45-atom ligand clean out of its site every step.
      double gyration = 0;
      for (std::uint32_t i : ligand.heavy) gyration += dot(local_model.base[i], local_model.base[i]);
      gyration = std::max(std::sqrt(gyration / static_cast<double>(ligand.heavy.size())), 1.0);
      // The walk itself compares total energies, not energies over the N_rot factor: dividing by it would
      // make the Metropolis temperature rise with the number of rotatable bonds, which is backwards.
      randomize(current);
      double current_rank = minimize_pose(local_hunt, current, options.local, x).total();
      auto remember = [&](double rank, const Energy& en, const std::vector<Vec3>& pose) {
        // Keep a short, diverse list: a new pose replaces a nearby one only when it is better.
        for (Candidate& c : out)
          if (heavy_rmsd(ligand, symmetry, c.coordinates, pose) < options.min_rmsd) {
            if (rank < c.rank) c = {rank, en.intermolecular, en.internal, pose};
            return;
          }
        out.push_back({rank, en.intermolecular, en.internal, pose});
        if (out.size() > 24) {
          std::sort(out.begin(), out.end(), [](const Candidate& a, const Candidate& b) { return a.rank < b.rank; });
          out.resize(24);
        }
      };
      // A pose worth keeping is minimised again on the function as it is, which is the energy it is then
      // ranked and reported by; the hunting energy only decides where the walk goes next.
      auto keep = [&](Pose p) {
        const Energy fe = minimize_pose(local_model, p, options.local, refined_x);
        remember(rank_of(fe), fe, refined_x);
      };
      keep(current);
      double best_hunt = current_rank;
      for (int step = 0; step < steps; ++step) {
        Pose candidate = current;
        const std::size_t which = rng.index(2 + ligand.torsions.size());
        if (which == 0) {
          candidate.translation = candidate.translation + 2.0 * rng.in_unit_ball();
        } else if (which == 1 && rng.uniform() < 0.25) {
          // Now and then a wholly new orientation. The small turn below is scaled by the ligand's size, so
          // for a large one it is at most some 25 degrees and a flip of the molecule end over end is out of
          // its reach: in a site with a two-fold symmetry, which protein dimers have, the walk would then
          // never leave whichever of the two orientations it happened to start in.
          candidate.orientation = quat_normalize(rng.orientation());
        } else if (which == 1) {
          candidate.orientation = quat_normalize(quat_multiply(
              quat_from_rotation_vector((2.0 / gyration) * rng.in_unit_ball()), candidate.orientation));
        } else if (rng.uniform() < 0.5) {
          candidate.torsions[which - 2] = rng.uniform(-kPi, kPi);
        } else {
          // Half the torsion moves nudge the angle instead of redrawing it. With a dozen torsions, redrawing
          // one from scratch undoes a conformation that is nearly right and is almost always rejected, so the
          // walk can reach a good conformation but never polish it.
          candidate.torsions[which - 2] += 0.35 * rng.normal();
        }
        const double cr = minimize_pose(local_hunt, candidate, options.local, x).total();
        if (cr < current_rank || rng.uniform() < std::exp(-(cr - current_rank) / std::max(options.temperature, 1e-6))) {
          current = candidate;
          current_rank = cr;
          if (cr < best_hunt || out.size() < 24) {
            keep(current);
            best_hunt = std::min(best_hunt, cr);
          }
        }
      }
      std::sort(out.begin(), out.end(), [](const Candidate& a, const Candidate& b) { return a.rank < b.rank; });
      run_evaluations[run] = local_model.evaluations + local_hunt.evaluations;
    });
    for (std::size_t k : run_evaluations) evaluations += k;
  }

  // Pool the runs in run order, then cluster: a pose is reported only if it is more than min_rmsd from every
  // better pose already reported.
  std::vector<Candidate> pool;
  for (const auto& v : per_run) pool.insert(pool.end(), v.begin(), v.end());
  std::stable_sort(pool.begin(), pool.end(), [](const Candidate& a, const Candidate& b) { return a.rank < b.rank; });
  std::vector<Candidate> kept;
  for (const Candidate& c : pool) {
    if (!kept.empty() && c.rank > kept.front().rank + options.energy_range) break;
    bool near = false;
    for (const Candidate& k : kept)
      if (heavy_rmsd(ligand, symmetry, k.coordinates, c.coordinates) < options.min_rmsd) {
        near = true;
        break;
      }
    if (near) continue;
    kept.push_back(c);
    if (kept.size() >= static_cast<std::size_t>(options.modes)) break;
  }
  if (kept.empty()) {
    diag.error(DiagCode::NotConverged, "the search found no pose inside the box");
    return std::nullopt;
  }

  // Report exact scores, not interpolated ones.
  PointCloud reference_pose;
  if (reference) {
    reference_pose.resize(reference->coordinates.size());
    for (std::size_t i = 0; i < reference->coordinates.size(); ++i)
      reference_pose[i] = {reference->coordinates[i][0], reference->coordinates[i][1], reference->coordinates[i][2]};
  }
  for (const Candidate& c : kept) {
    DockedPose p = exact_score(receptor, cells, ligand, c.coordinates, &model, options.weights);
    if (reference) {
      PointCloud pose(c.coordinates.size());
      for (std::size_t i = 0; i < c.coordinates.size(); ++i)
        pose[i] = {c.coordinates[i][0], c.coordinates[i][1], c.coordinates[i][2]};
      p.rmsd = symmetry_rmsd(ligand, pose, *reference, reference_pose);
    }
    result.poses.push_back(std::move(p));
  }

  // The reference pose minimised where it stands: the yardstick that separates a search failure from a
  // scoring failure.
  if (reference && !options.refine_only) {
    DockOptions ro = options;
    ro.refine_only = true;
    Model rm = build_model(*reference, &grid, options.weights);
    Pose p;
    p.torsions.assign(reference->torsions.size(), 0.0);
    p.translation = rm.centre;
    std::vector<Vec3> x;
    minimize_pose(rm, p, options.local, x);
    const DockedPose rpose = exact_score(receptor, cells, *reference, x, &rm, options.weights);
    result.reference_affinity = rpose.affinity;
    PointCloud refined(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) refined[i] = {x[i][0], x[i][1], x[i][2]};
    result.reference_rmsd = symmetry_rmsd(*reference, refined, *reference, reference_pose);
    result.has_reference = true;
  }

  result.evaluations = evaluations;
  result.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
  return result;
}

// ---- output -------------------------------------------------------------------------------------

std::string write_poses_sdf(const DockLigand& ligand, const std::vector<DockedPose>& poses) {
  std::string out;
  char buf[96];
  for (std::size_t k = 0; k < poses.size(); ++k) {
    const DockedPose& p = poses[k];
    if (p.coordinates.size() != ligand.coordinates.size()) continue;
    ChemicalSystemState st = ligand.state;
    const Species& s = st.species[ligand.species];
    for (std::uint32_t i = 0; i < s.atom_count; ++i) {
      st.atoms.x[s.atom_begin + i] = p.coordinates[i][0];
      st.atoms.y[s.atom_begin + i] = p.coordinates[i][1];
      st.atoms.z[s.atom_begin + i] = p.coordinates[i][2];
    }
    std::string record = write_sdf(st, ligand.species, "pose " + std::to_string(k + 1));
    const std::string tail = "$$$$\n";
    if (record.size() >= tail.size() && record.compare(record.size() - tail.size(), tail.size(), tail) == 0)
      record.resize(record.size() - tail.size());
    std::snprintf(buf, sizeof buf, "%.3f", p.affinity);
    record += "> <affinity_kcal_mol>\n";
    record += buf;
    record += "\n\n";
    if (p.rmsd >= 0) {
      std::snprintf(buf, sizeof buf, "%.3f", p.rmsd);
      record += "> <rmsd_reference>\n";
      record += buf;
      record += "\n\n";
    }
    record += "$$$$\n";
    out += record;
  }
  return out;
}

// ---- JSON ---------------------------------------------------------------------------------------

namespace {

std::string jnum(double v) {
  char b[40];
  std::snprintf(b, sizeof b, "%.6g", v);
  return std::isfinite(v) ? b : "null";
}

std::string fail_json(const Diagnostics& d) {
  return std::string("{\"ok\":false,\"diagnostics\":") + diagnostics_to_json(d) + "}";
}

// A molfile is at least four lines (title, program, comment, counts); a SMILES string is one, give or take
// a trailing newline, and a residue name is a word.
bool looks_like_molfile(std::string_view s) { return std::count(s.begin(), s.end(), '\n') >= 3; }

}  // namespace

std::string dock_json(std::string_view receptor_text, std::string_view ligand_text, const DockJsonOptions& options) {
  Diagnostics diag;
  const std::optional<Protein> protein = read_protein(receptor_text, diag);
  if (!protein) {
    if (!diag.has_errors()) diag.error(DiagCode::UnknownFileFormat, "the receptor is not a PDB or mmCIF file");
    return fail_json(diag);
  }

  // The reference ligand: a crystal pose as MOL/SDF, or the name of a hetero residue of the file. Either
  // places the box and clears its own site; the first is also the pose the RMSD is measured from.
  PointCloud reference_atoms;
  std::optional<DockLigand> reference;
  if (!options.reference.empty()) {
    if (looks_like_molfile(options.reference)) {
      reference = ligand_from_molfile(options.reference, diag, options.dock.seed);
      if (!reference) return fail_json(diag);
      for (std::uint32_t i : reference->heavy) reference_atoms.push_back(reference->coordinates[i]);
    } else {
      reference_atoms = protein_ligand_atoms(*protein, options.reference);
      if (reference_atoms.empty()) {
        diag.error(DiagCode::InvalidInput, "the receptor file has no residue named " + options.reference);
        return fail_json(diag);
      }
    }
  }

  DockOptions dock_options = options.dock;
  if (!dock_options.box.valid()) {
    if (reference_atoms.empty()) {
      diag.error(DiagCode::InvalidInput, "give a search box, or a reference ligand to place one around");
      return fail_json(diag);
    }
    dock_options.box = box_around(reference_atoms, options.box_margin);
  }

  const DockReceptor receptor = prepare_receptor(*protein, options.receptor, reference_atoms, diag);
  if (receptor.atoms.empty()) return fail_json(diag);

  const std::optional<DockLigand> ligand = looks_like_molfile(ligand_text)
                                               ? ligand_from_molfile(ligand_text, diag, dock_options.seed)
                                               : ligand_from_smiles(ligand_text, diag, dock_options.seed);
  if (!ligand) return fail_json(diag);

  const std::optional<DockResult> result =
      dock(receptor, *ligand, dock_options, diag, reference ? &*reference : nullptr);
  if (!result) return fail_json(diag);

  std::string out = "{\"ok\":true,\"smiles\":" + json_escape(ligand->smiles) +
                    ",\"heavy_atoms\":" + std::to_string(result->heavy_atoms) +
                    ",\"rotatable_bonds\":" + std::to_string(result->rotatable_bonds) +
                    ",\"receptor_atoms\":" + std::to_string(result->receptor_atoms) + ",\"box\":{\"centre\":[" +
                    jnum(result->box.cx) + "," + jnum(result->box.cy) + "," + jnum(result->box.cz) + "],\"size\":[" +
                    jnum(result->box.sx) + "," + jnum(result->box.sy) + "," + jnum(result->box.sz) + "]}";
  out += ",\"grid_spacing\":" + jnum(result->grid_spacing) + ",\"runs\":" + std::to_string(result->runs) +
         ",\"steps\":" + std::to_string(result->steps) +
         ",\"evaluations\":" + std::to_string(result->evaluations) + ",\"seconds\":" + jnum(result->seconds);
  out += ",\"receptor\":{\"waters_removed\":" + std::to_string(receptor.waters_removed) +
         ",\"hetero_removed\":" + std::to_string(receptor.hetero_removed) + ",\"donors\":" +
         std::to_string(receptor.donors) + ",\"acceptors\":" + std::to_string(receptor.acceptors) +
         ",\"hydrophobic\":" + std::to_string(receptor.hydrophobic) + ",\"kept_hetero\":[";
  for (std::size_t i = 0; i < receptor.kept_hetero.size(); ++i)
    out += (i ? "," : "") + json_escape(receptor.kept_hetero[i]);
  out += "]}";
  out += ",\"poses\":[";
  for (std::size_t k = 0; k < result->poses.size(); ++k) {
    const DockedPose& p = result->poses[k];
    out += (k ? ",{" : "{");
    out += "\"affinity\":" + jnum(p.affinity) + ",\"intermolecular\":" + jnum(p.intermolecular) +
           ",\"internal\":" + jnum(p.internal) + ",\"rmsd\":" + (p.rmsd >= 0 ? jnum(p.rmsd) : "null") +
           ",\"terms\":{\"gauss1\":" + jnum(p.terms.gauss1) + ",\"gauss2\":" + jnum(p.terms.gauss2) +
           ",\"repulsion\":" + jnum(p.terms.repulsion) + ",\"hydrophobic\":" + jnum(p.terms.hydrophobic) +
           ",\"hbond\":" + jnum(p.terms.hbond) + "},\"coordinates\":[";
    for (std::size_t i = 0; i < p.coordinates.size(); ++i)
      out += (i ? ",[" : "[") + jnum(p.coordinates[i][0]) + "," + jnum(p.coordinates[i][1]) + "," +
             jnum(p.coordinates[i][2]) + "]";
    out += "]}";
  }
  out += "]";
  if (!options.reference.empty()) {
    out += ",\"reference\":{\"atoms\":" + std::to_string(reference_atoms.size());
    if (reference) out += ",\"smiles\":" + json_escape(reference->smiles);
    if (result->has_reference)
      out += ",\"affinity\":" + jnum(result->reference_affinity) + ",\"rmsd\":" + jnum(result->reference_rmsd);
    out += "}";
  }
  if (options.poses_sdf) out += ",\"sdf\":" + json_escape(write_poses_sdf(*ligand, result->poses));
  out += ",\"diagnostics\":" + diagnostics_to_json(diag) + "}";
  return out;
}

}  // namespace chemlarp
