#include "chemlarp/conformer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "chemlarp/canonical.hpp"
#include "chemlarp/elements.hpp"
#include "chemlarp/molio.hpp"
#include "chemlarp/rings.hpp"
#include "chemlarp/smiles.hpp"

namespace chemlarp {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;

// Bondi van der Waals radii (Å).
double vdw_radius(std::uint8_t z) {
  switch (z) {
    case 1: return 1.20;
    case 5: return 1.92;
    case 6: return 1.70;
    case 7: return 1.55;
    case 8: return 1.52;
    case 9: return 1.47;
    case 14: return 2.10;
    case 15: return 1.80;
    case 16: return 1.80;
    case 17: return 1.75;
    case 35: return 1.85;
    case 53: return 1.98;
    default: return 2.0;
  }
}

// Platform-independent uniform [0, 1): std::mt19937_64 is fully specified by the
// standard, std::uniform_real_distribution is not.
double uniform01(std::mt19937_64& rng) { return static_cast<double>(rng() >> 11) * 0x1.0p-53; }

struct Matrix {
  std::size_t n = 0;
  std::vector<double> a;
  explicit Matrix(std::size_t size = 0, double fill = 0.0) : n(size), a(size * size, fill) {}
  double& operator()(std::size_t i, std::size_t j) { return a[i * n + j]; }
  double operator()(std::size_t i, std::size_t j) const { return a[i * n + j]; }
};

struct Bounds {
  Matrix lo, up;
};

double cis_trans_distance(double r1, double r2, double r3, double th1, double th2, double phi) {
  // j at origin, k on +x; i in the xy-plane at angle th1 from j->k; l at angle th2 from k->j, torsion phi.
  const double ix = r1 * std::cos(th1), iy = r1 * std::sin(th1);
  const double lx = r2 - r3 * std::cos(th2);
  const double ly = r3 * std::sin(th2) * std::cos(phi), lz = r3 * std::sin(th2) * std::sin(phi);
  return std::sqrt((lx - ix) * (lx - ix) + (ly - iy) * (ly - iy) + lz * lz);
}

struct Topology {
  std::size_t n = 0;
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> adj;  // (neighbour, global bond)
  std::vector<std::vector<std::uint32_t>> rings;                          // SSSR, local atom indices
};

Topology topology(const ChemicalSystemState& st, const Species& s, std::uint32_t sp) {
  Topology t;
  t.n = s.atom_count;
  t.adj.resize(t.n);
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    const std::uint32_t i = st.bonds.begin[b] - s.atom_begin, j = st.bonds.end[b] - s.atom_begin;
    t.adj[i].push_back({j, b});
    t.adj[j].push_back({i, b});
  }
  for (const auto& r : find_sssr(st, sp)) {
    std::vector<std::uint32_t> loc;
    for (std::uint32_t a : r.atoms) loc.push_back(a - s.atom_begin);
    t.rings.push_back(std::move(loc));
  }
  return t;
}

// Smallest ring containing all of the given (consecutive) atoms, 0 if none.
std::size_t ring_size_with(const Topology& t, std::initializer_list<std::uint32_t> atoms) {
  std::size_t best = 0;
  for (const auto& r : t.rings) {
    bool all = true;
    for (std::uint32_t a : atoms) all = all && std::find(r.begin(), r.end(), a) != r.end();
    if (all && (best == 0 || r.size() < best)) best = r.size();
  }
  return best;
}

bool build_bounds(const ChemicalSystemState& st, std::uint32_t sp, const std::vector<int>& types, const Topology& t,
                  Bounds& B, Diagnostics& diag) {
  const Species& s = st.species[sp];
  const std::size_t n = t.n;
  B.lo = Matrix(n, 0.0);
  B.up = Matrix(n, 1000.0);
  Matrix set(n, 0.0);  // 1 = 1-2, 2 = 1-3, 3 = 1-4 bound already fixed

  auto r0 = [&](std::uint32_t i, std::uint32_t j) {
    for (auto [k, b] : t.adj[i])
      if (k == j) return uff_rest_length(types[i], types[j], uff_bond_order(st, b));
    return 1.5;
  };
  auto angle = [&](std::uint32_t i, std::uint32_t j, std::uint32_t k) -> std::pair<double, double> {
    if (t.adj[j].size() > 4) return {85.0 * kDeg, 180.0 * kDeg};  // hypervalent centre: 90 or 180
    const std::size_t ring = ring_size_with(t, {i, j, k});
    if (ring >= 3 && ring <= 5) {
      const double th = kPi * static_cast<double>(ring - 2) / static_cast<double>(ring);
      return {th, th};
    }
    double th = uff_parameters(types[j]).theta0 * kDeg;
    if (t.adj[j].size() == 3 && uff_parameters(types[j]).theta0 > 115.0 && uff_parameters(types[j]).theta0 < 125.0)
      th = 120.0 * kDeg;
    return {th, th};
  };
  auto setb = [&](std::uint32_t i, std::uint32_t j, double lo, double up, int level) {
    if (set(i, j) != 0 && set(i, j) <= level) return;
    B.lo(i, j) = B.lo(j, i) = lo;
    B.up(i, j) = B.up(j, i) = up;
    set(i, j) = set(j, i) = level;
  };

  // 1-2
  for (std::uint32_t i = 0; i < n; ++i)
    for (auto [j, b] : t.adj[i]) {
      const double d = uff_rest_length(types[i], types[j], uff_bond_order(st, b));
      setb(i, j, d - 0.01, d + 0.01, 1);
    }
  // 1-3
  for (std::uint32_t j = 0; j < n; ++j) {
    for (std::size_t a = 0; a < t.adj[j].size(); ++a) {
      for (std::size_t c = a + 1; c < t.adj[j].size(); ++c) {
        const std::uint32_t i = t.adj[j][a].first, k = t.adj[j][c].first;
        const double rij = r0(i, j), rjk = r0(j, k);
        auto [tmin, tmax] = angle(i, j, k);
        const double dmin = std::sqrt(rij * rij + rjk * rjk - 2 * rij * rjk * std::cos(tmin));
        const double dmax = std::sqrt(rij * rij + rjk * rjk - 2 * rij * rjk * std::cos(tmax));
        setb(i, k, dmin - 0.04, dmax + 0.04, 2);
      }
    }
  }
  // 1-4, with E/Z from the stereo table and cis enforced inside planar rings.
  std::vector<const DoubleBondStereo*> dbs(s.bond_count, nullptr);
  for (const auto& d : st.double_bonds)
    if (st.bonds.species[d.bond] == sp) dbs[d.bond - s.bond_begin] = &d;
  for (std::uint32_t j = 0; j < n; ++j) {
    for (auto [k, bjk] : t.adj[j]) {
      if (k < j) continue;
      const bool planar = (st.bonds.flags[bjk] & kBondAromatic) || st.bonds.order[bjk] == 2;
      for (auto [i, b1] : t.adj[j]) {
        if (i == k) continue;
        for (auto [l, b2] : t.adj[k]) {
          if (l == j || l == i) continue;
          const double r1 = r0(i, j), r2 = r0(j, k), r3 = r0(k, l);
          const double th1 = angle(i, j, k).first, th2 = angle(j, k, l).first;
          const double dc = cis_trans_distance(r1, r2, r3, th1, th2, 0.0);
          const double dt = cis_trans_distance(r1, r2, r3, th1, th2, kPi);
          double lo = std::min(dc, dt) - 0.05, up = std::max(dc, dt) + 0.05;
          if (const DoubleBondStereo* d = dbs[bjk - s.bond_begin]) {
            const std::uint32_t gb = st.bonds.begin[bjk] - s.atom_begin;
            const std::uint32_t ri = (gb == j ? d->ref_begin : d->ref_end) - s.atom_begin;
            const std::uint32_t rl = (gb == j ? d->ref_end : d->ref_begin) - s.atom_begin;
            const bool cis = d->cis != (i != ri) != (l != rl);
            const double target = cis ? dc : dt;
            lo = target - 0.05;
            up = target + 0.05;
          } else if (planar && ring_size_with(t, {i, j, k, l}) > 0) {
            lo = dc - 0.05;
            up = dc + 0.05;
          }
          setb(i, l, lo, up, 3);
        }
      }
    }
  }
  // van der Waals lower bounds for everything else.
  for (std::uint32_t i = 0; i < n; ++i)
    for (std::uint32_t j = i + 1; j < n; ++j)
      if (set(i, j) == 0) {
        const double lo = 0.7 * (vdw_radius(st.atoms.element[s.atom_begin + i]) + vdw_radius(st.atoms.element[s.atom_begin + j]));
        B.lo(i, j) = B.lo(j, i) = lo;
      }

  // Triangle smoothing (Floyd-Warshall on upper bounds, then lower bounds).
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      if (i == k) continue;
      for (std::size_t j = i + 1; j < n; ++j) {
        if (j == k) continue;
        double& u = B.up(i, j);
        const double uk = B.up(i, k) + B.up(k, j);
        if (u > uk) u = B.up(j, i) = uk;
        double& l = B.lo(i, j);
        const double l1 = B.lo(i, k) - B.up(k, j), l2 = B.lo(j, k) - B.up(k, i);
        if (l < l1) l = B.lo(j, i) = l1;
        if (l < l2) l = B.lo(j, i) = l2;
        if (B.lo(i, j) > B.up(i, j) + 1e-6) {
          diag.error(DiagCode::InconsistentBounds, "distance bounds between atoms " + std::to_string(i) + " and " +
                                                       std::to_string(j) + " are inconsistent after triangle smoothing")
              .on_species(sp);
          return false;
        }
      }
    }
  }
  // Disconnected fragments: keep them a few Å apart.
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j)
      if (B.up(i, j) > 999.0) B.up(i, j) = B.up(j, i) = B.lo(i, j) + 6.0;
  return true;
}

struct ChiralTerm {
  std::array<std::uint32_t, 4> idx;  // local atoms (lone pair replaced by the centre)
  double lo, up;
};

struct PlanarTerm {
  std::array<std::uint32_t, 4> idx;  // centre and its three neighbours
};

double volume(const double* x, std::size_t dim, const std::array<std::uint32_t, 4>& p, std::array<std::array<double, 3>, 4>* grad) {
  auto P = [&](int k, int c) { return x[p[static_cast<std::size_t>(k)] * dim + static_cast<std::size_t>(c)]; };
  double v1[3], v2[3], v3[3];
  for (int c = 0; c < 3; ++c) {
    v1[c] = P(1, c) - P(0, c);
    v2[c] = P(2, c) - P(0, c);
    v3[c] = P(3, c) - P(0, c);
  }
  const double cx[3] = {v2[1] * v3[2] - v2[2] * v3[1], v2[2] * v3[0] - v2[0] * v3[2], v2[0] * v3[1] - v2[1] * v3[0]};
  const double V = v1[0] * cx[0] + v1[1] * cx[1] + v1[2] * cx[2];
  if (grad) {
    const double g2[3] = {v3[1] * v1[2] - v3[2] * v1[1], v3[2] * v1[0] - v3[0] * v1[2], v3[0] * v1[1] - v3[1] * v1[0]};
    const double g3[3] = {v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]};
    for (int c = 0; c < 3; ++c) {
      (*grad)[1][static_cast<std::size_t>(c)] = cx[c];
      (*grad)[2][static_cast<std::size_t>(c)] = g2[c];
      (*grad)[3][static_cast<std::size_t>(c)] = g3[c];
      (*grad)[0][static_cast<std::size_t>(c)] = -(cx[c] + g2[c] + g3[c]);
    }
  }
  return V;
}

// Distance-geometry refinement energy in `dim` dimensions.
double dg_energy(const std::vector<double>& x, std::vector<double>& g, std::size_t n, std::size_t dim, const Bounds& B,
                 const std::vector<ChiralTerm>& chiral, const std::vector<PlanarTerm>& planar, double w4) {
  g.assign(x.size(), 0.0);
  double E = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      double d2 = 0.0;
      for (std::size_t c = 0; c < dim; ++c) {
        const double t = x[i * dim + c] - x[j * dim + c];
        d2 += t * t;
      }
      const double u = B.up(i, j), l = B.lo(i, j);
      double dEdd2 = 0.0;
      if (d2 > u * u) {
        const double f = d2 / (u * u) - 1.0;
        E += f * f;
        dEdd2 = 2.0 * f / (u * u);
      } else if (d2 < l * l) {
        const double den = l * l + d2;
        const double f = 2.0 * l * l / den - 1.0;
        E += f * f;
        dEdd2 = 2.0 * f * (-2.0 * l * l / (den * den));
      } else {
        continue;
      }
      for (std::size_t c = 0; c < dim; ++c) {
        const double t = 2.0 * dEdd2 * (x[i * dim + c] - x[j * dim + c]);
        g[i * dim + c] += t;
        g[j * dim + c] -= t;
      }
    }
  }
  std::array<std::array<double, 3>, 4> vg{};
  for (const auto& ct : chiral) {
    const double V = volume(x.data(), dim, ct.idx, &vg);
    double dv = 0.0;
    if (V < ct.lo) dv = V - ct.lo;
    else if (V > ct.up) dv = V - ct.up;
    else continue;
    E += dv * dv;
    for (int k = 0; k < 4; ++k)
      for (int c = 0; c < 3; ++c) g[ct.idx[static_cast<std::size_t>(k)] * dim + static_cast<std::size_t>(c)] += 2.0 * dv * vg[static_cast<std::size_t>(k)][static_cast<std::size_t>(c)];
  }
  for (const auto& pt : planar) {
    const double V = volume(x.data(), dim, pt.idx, &vg);
    const double w = 0.1;
    E += w * V * V;
    for (int k = 0; k < 4; ++k)
      for (int c = 0; c < 3; ++c) g[pt.idx[static_cast<std::size_t>(k)] * dim + static_cast<std::size_t>(c)] += 2.0 * w * V * vg[static_cast<std::size_t>(k)][static_cast<std::size_t>(c)];
  }
  if (dim == 4 && w4 > 0.0) {
    for (std::size_t i = 0; i < n; ++i) {
      const double z4 = x[i * 4 + 3];
      E += w4 * z4 * z4;
      g[i * 4 + 3] += 2.0 * w4 * z4;
    }
  }
  return E;
}

// Top-k eigenpairs of a symmetric matrix by shifted power iteration with deflation.
bool top_eigen(const Matrix& G, std::size_t k, std::mt19937_64& rng, std::vector<double>& vals,
               std::vector<std::vector<double>>& vecs) {
  const std::size_t n = G.n;
  double shift = 0.0;  // Gershgorin: make the matrix positive semidefinite
  for (std::size_t i = 0; i < n; ++i) {
    double r = 0.0;
    for (std::size_t j = 0; j < n; ++j)
      if (j != i) r += std::fabs(G(i, j));
    shift = std::max(shift, r - G(i, i));
  }
  vals.clear();
  vecs.clear();
  std::vector<double> v(n), w(n);
  for (std::size_t e = 0; e < std::min(k, n); ++e) {
    for (double& x : v) x = uniform01(rng) - 0.5;
    double lambda = 0.0;
    for (int it = 0; it < 2000; ++it) {
      for (std::size_t i = 0; i < n; ++i) {
        double s = shift * v[i];
        for (std::size_t j = 0; j < n; ++j) s += G(i, j) * v[j];
        w[i] = s;
      }
      for (std::size_t p = 0; p < vecs.size(); ++p) {  // deflate found eigenvectors
        double d = 0.0;
        for (std::size_t i = 0; i < n; ++i) d += vecs[p][i] * v[i];
        for (std::size_t i = 0; i < n; ++i) w[i] -= (vals[p] + shift) * d * vecs[p][i];
      }
      double nw = 0.0;
      for (double x : w) nw += x * x;
      nw = std::sqrt(nw);
      if (nw < 1e-300) break;
      double change = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        const double nv = w[i] / nw;
        change = std::max(change, std::fabs(nv - v[i]));
        v[i] = nv;
      }
      lambda = nw - shift;
      if (change < 1e-9) break;
    }
    vals.push_back(lambda);
    vecs.push_back(v);
  }
  return !vals.empty() && vals[0] > 0.0;
}

struct Attempt {
  std::vector<double> xyz;  // 3n
  double energy = 0.0;
  bool chirality_ok = false;
};

Attempt embed_once(const Bounds& B, std::size_t n, const std::vector<ChiralTerm>& chiral,
                   const std::vector<PlanarTerm>& planar, std::mt19937_64& rng) {
  Attempt out;
  // Random distances within the bounds, then the metric matrix from centroid distances.
  Matrix D2(n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j) {
      const double d = B.lo(i, j) + uniform01(rng) * (B.up(i, j) - B.lo(i, j));
      D2(i, j) = D2(j, i) = d * d;
    }
  double sum_all = 0.0;
  for (std::size_t j = 0; j < n; ++j)
    for (std::size_t k = j + 1; k < n; ++k) sum_all += D2(j, k);
  std::vector<double> d0(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    double s = 0.0;
    for (std::size_t j = 0; j < n; ++j) s += D2(i, j);
    d0[i] = s / static_cast<double>(n) - sum_all / static_cast<double>(n * n);
  }
  Matrix G(n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) G(i, j) = 0.5 * (d0[i] + d0[j] - D2(i, j));
  std::vector<double> vals;
  std::vector<std::vector<double>> vecs;
  const std::size_t dim = 4;
  std::vector<double> x(n * dim, 0.0);
  if (n == 1) {
    out.xyz.assign(3, 0.0);
    out.chirality_ok = true;
    return out;
  }
  if (!top_eigen(G, dim, rng, vals, vecs)) {
    out.energy = 1e30;
    return out;
  }
  for (std::size_t e = 0; e < vals.size(); ++e) {
    const double s = std::sqrt(std::max(vals[e], 0.0));
    for (std::size_t i = 0; i < n; ++i) x[i * dim + e] = s * vecs[e][i];
  }
  // Small random jitter breaks exact degeneracies (e.g. all-planar starts).
  for (double& c : x) c += 0.05 * (uniform01(rng) - 0.5);

  LbfgsOptions lo;
  lo.max_iterations = 1500;
  lo.grad_tol = 1e-5;
  lo.max_step = 0.5;
  // Stage 1: bounds + chirality in 4D. Stage 2: squeeze out the 4th dimension.
  minimize_lbfgs([&](const std::vector<double>& xx, std::vector<double>& g) { return dg_energy(xx, g, n, 4, B, chiral, planar, 0.0); }, x, lo);
  minimize_lbfgs([&](const std::vector<double>& xx, std::vector<double>& g) { return dg_energy(xx, g, n, 4, B, chiral, planar, 1.0); }, x, lo);
  std::vector<double> x3(n * 3);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t c = 0; c < 3; ++c) x3[i * 3 + c] = x[i * dim + c];
  // Final 3D polish of the bounds.
  auto r = minimize_lbfgs([&](const std::vector<double>& xx, std::vector<double>& g) { return dg_energy(xx, g, n, 3, B, chiral, planar, 0.0); }, x3, lo);
  out.energy = r.f;
  out.chirality_ok = true;
  for (const auto& ct : chiral) {
    const double V = volume(x3.data(), 3, ct.idx, nullptr);
    if ((ct.lo > 0 && V <= 0) || (ct.up < 0 && V >= 0)) out.chirality_ok = false;
  }
  out.xyz = std::move(x3);
  return out;
}

void copy_coordinates(ChemicalSystemState& st, const Species& s, const std::vector<double>& xyz) {
  for (std::uint32_t i = 0; i < s.atom_count; ++i) {
    st.atoms.x[s.atom_begin + i] = xyz[3 * i];
    st.atoms.y[s.atom_begin + i] = xyz[3 * i + 1];
    st.atoms.z[s.atom_begin + i] = xyz[3 * i + 2];
  }
  ++st.revision;
}

}  // namespace

// ---- explicit hydrogens ------------------------------------------------------

std::int32_t add_explicit_hydrogens(ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag) {
  const Species old = st.species[sp];
  const std::uint32_t new_sp = static_cast<std::uint32_t>(st.species.size());
  const std::uint32_t atom_base = static_cast<std::uint32_t>(st.atoms.size());
  const std::uint32_t bond_base = static_cast<std::uint32_t>(st.bonds.size());
  std::vector<std::uint32_t> first_h(old.atom_count, UINT32_MAX);
  std::uint32_t added = 0;

  for (std::uint32_t i = 0; i < old.atom_count; ++i) {
    Atom a = st.atoms.get(old.atom_begin + i);
    a.implicit_h = 0;
    a.species = new_sp;
    st.atoms.push(a);
  }
  for (std::uint32_t b = 0; b < old.bond_count; ++b) {
    Bond bd = st.bonds.get(old.bond_begin + b);
    bd.begin = atom_base + (bd.begin - old.atom_begin);
    bd.end = atom_base + (bd.end - old.atom_begin);
    bd.species = new_sp;
    st.bonds.push(bd);
  }
  for (std::uint32_t i = 0; i < old.atom_count; ++i) {
    const int h = st.atoms.implicit_h[old.atom_begin + i];
    for (int k = 0; k < h; ++k) {
      Atom H;
      H.element = 1;
      H.species = new_sp;
      const std::uint32_t idx = atom_base + old.atom_count + added++;
      if (k == 0) first_h[i] = idx;
      st.atoms.push(H);
      st.bonds.push(Bond{atom_base + i, idx, 1, 0, new_sp});
    }
  }
  auto remap = [&](std::uint32_t ref, std::uint32_t centre_local) -> std::uint32_t {
    if (ref == kImplicitHRef) return first_h[centre_local];
    if (ref == kLonePairRef) return kLonePairRef;
    return atom_base + (ref - old.atom_begin);
  };
  const std::size_t ntet = st.tetrahedral.size();
  for (std::size_t t = 0; t < ntet; ++t) {
    const TetrahedralStereo ts = st.tetrahedral[t];
    if (st.atoms.species[ts.atom] != sp) continue;
    TetrahedralStereo n2 = ts;
    const std::uint32_t c = ts.atom - old.atom_begin;
    n2.atom = atom_base + c;
    for (auto& r : n2.refs) r = remap(r, c);
    st.tetrahedral.push_back(n2);
  }
  const std::size_t nax = st.axial.size();
  for (std::size_t t = 0; t < nax; ++t) {
    const AxialStereo a = st.axial[t];
    if (st.atoms.species[a.begin] != sp) continue;
    AxialStereo n2 = a;
    const std::uint32_t b0 = a.begin - old.atom_begin, e0 = a.end - old.atom_begin;
    n2.begin = atom_base + b0;
    n2.end = atom_base + e0;
    if (a.centre != kNoAtom) n2.centre = atom_base + (a.centre - old.atom_begin);
    for (std::size_t q = 0; q < 4; ++q) n2.refs[q] = remap(a.refs[q], q < 2 ? b0 : e0);
    st.axial.push_back(n2);
  }
  const std::size_t npl = st.planar.size();
  for (std::size_t t = 0; t < npl; ++t) {
    const PlanarStereo u = st.planar[t];
    if (u.atoms.empty() || st.atoms.species[u.atoms[0]] != sp) continue;
    PlanarStereo n2 = u;
    for (auto& a : n2.atoms) a = atom_base + (a - old.atom_begin);
    for (auto& a : n2.pilots) a = atom_base + (a - old.atom_begin);
    st.planar.push_back(n2);
  }
  const std::size_t ndb = st.double_bonds.size();
  for (std::size_t t = 0; t < ndb; ++t) {
    const DoubleBondStereo d = st.double_bonds[t];
    if (st.bonds.species[d.bond] != sp) continue;
    st.double_bonds.push_back({bond_base + (d.bond - old.bond_begin), atom_base + (d.ref_begin - old.atom_begin),
                               atom_base + (d.ref_end - old.atom_begin), d.cis});
  }
  Species s = old;
  s.atom_begin = atom_base;
  s.atom_count = old.atom_count + added;
  s.bond_begin = bond_base;
  s.bond_count = old.bond_count + added;
  st.append_species(std::move(s));
  st.rebuild_adjacency();
  st.species[new_sp].canonical_smiles = canonical_smiles(st, new_sp, diag);
  st.moles[new_sp] = st.moles[sp];
  return static_cast<std::int32_t>(new_sp);
}

// ---- embedding ---------------------------------------------------------------

namespace {
// Largest bond angle (degrees) at any 4-coordinate carbon. Real sp3 carbons stay
// far below 150 degrees even in strained cages; above it the structure is
// inverted (a bridge pulled through its ring), which happens when the requested
// stereo cannot exist, e.g. an impossible bridgehead pair in norbornane.
double worst_sp3_carbon_angle(const ChemicalSystemState& st, const Species& s) {
  double worst = 0.0;
  for (std::uint32_t g = s.atom_begin; g < s.atom_begin + s.atom_count; ++g) {
    if (st.atoms.element[g] != 6) continue;
    const auto nb = st.adjacency.neighbors(g);
    if (nb.size() != 4) continue;
    for (std::size_t i = 0; i < 4; ++i)
      for (std::size_t j = i + 1; j < 4; ++j) {
        const double ux = st.atoms.x[nb[i]] - st.atoms.x[g], uy = st.atoms.y[nb[i]] - st.atoms.y[g], uz = st.atoms.z[nb[i]] - st.atoms.z[g];
        const double vx = st.atoms.x[nb[j]] - st.atoms.x[g], vy = st.atoms.y[nb[j]] - st.atoms.y[g], vz = st.atoms.z[nb[j]] - st.atoms.z[g];
        const double c = (ux * vx + uy * vy + uz * vz) /
                         std::sqrt((ux * ux + uy * uy + uz * uz) * (vx * vx + vy * vy + vz * vz) + 1e-300);
        worst = std::max(worst, std::acos(std::clamp(c, -1.0, 1.0)) / kDeg);
      }
  }
  return worst;
}
constexpr double kMaxSp3CarbonAngle = 150.0;
}  // namespace

EmbedResult embed_species(ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag, const EmbedOptions& opt) {
  EmbedResult res;
  const Species s = st.species[sp];
  for (std::uint32_t g = s.atom_begin; g < s.atom_begin + s.atom_count; ++g)
    if (st.atoms.implicit_h[g] > 0) {
      diag.error(DiagCode::ImplicitHydrogens, "embedding needs explicit hydrogens; call add_explicit_hydrogens first").on_species(sp);
      return res;
    }
  const std::vector<int> types = assign_uff_types(st, sp, diag);
  for (int t : types)
    if (t < 0) return res;
  const Topology topo = topology(st, s, sp);
  Bounds B;
  if (!build_bounds(st, sp, types, topo, B, diag)) return res;
  const std::size_t n = s.atom_count;

  std::vector<ChiralTerm> chiral;
  for (const auto& t : st.tetrahedral) {
    if (st.atoms.species[t.atom] != sp) continue;
    ChiralTerm ct;
    bool lone_pair = false;
    for (std::size_t k = 0; k < 4; ++k) {
      const std::uint32_t r = t.refs[k];
      lone_pair = lone_pair || r == kLonePairRef;
      ct.idx[k] = (r == kLonePairRef ? t.atom : r) - s.atom_begin;
    }
    const double vmin = lone_pair ? 0.3 : 1.0;
    if (t.parity == Chirality::CounterClockwise) { ct.lo = -100.0; ct.up = -vmin; }
    else { ct.lo = vmin; ct.up = 100.0; }
    chiral.push_back(ct);
  }
  // Axial units: the same signed volume over the four ligands of the elongated tetrahedron.
  for (const auto& a : st.axial) {
    if (st.atoms.species[a.begin] != sp) continue;
    ChiralTerm ct;
    for (std::size_t k = 0; k < 4; ++k) ct.idx[k] = a.refs[k] - s.atom_begin;
    if (a.parity == Chirality::CounterClockwise) { ct.lo = -100.0; ct.up = -1.0; }
    else { ct.lo = 1.0; ct.up = 100.0; }
    chiral.push_back(ct);
  }
  std::vector<PlanarTerm> planar;
  for (std::uint32_t i = 0; i < n; ++i) {
    if (topo.adj[i].size() != 3) continue;
    const std::string lab = uff_parameters(types[i]).label;
    if (lab == "C_2" || lab == "C_R" || lab == "N_2" || lab == "N_R" || lab == "B_2")
      planar.push_back({{i, topo.adj[i][0].first, topo.adj[i][1].first, topo.adj[i][2].first}});
  }

  // Each attempt: distance geometry -> UFF -> read stereo back from the 3D structure.
  // Collect up to `candidates` attempts that reproduce all stereo and keep the lowest
  // energy one (a single random start can land in e.g. a twist-boat); if none
  // reproduces everything, keep the best (most stereo correct, then lowest energy).
  auto score_stereo = [&](EmbedResult& r) {
    r.stereo_centres = r.stereo_centres_ok = r.stereo_bonds = r.stereo_bonds_ok = 0;
    for (const auto& t : st.tetrahedral) {
      if (st.atoms.species[t.atom] != sp) continue;
      ++r.stereo_centres;
      r.stereo_centres_ok += parity_from_coordinates(st, t) == t.parity;
    }
    for (const auto& a : st.axial) {
      if (st.atoms.species[a.begin] != sp) continue;
      ++r.stereo_centres;
      r.stereo_centres_ok += axial_parity_from_coordinates(st, a) == a.parity;
    }
    for (const auto& d : st.double_bonds) {
      if (st.bonds.species[d.bond] != sp) continue;
      ++r.stereo_bonds;
      r.stereo_bonds_ok += cis_from_coordinates(st, d) == d.cis;
    }
  };
  std::vector<double> best_xyz;
  EmbedResult best;
  std::size_t best_ok = 0;
  double best_energy = 1e300;
  int complete = 0, inverted = 0;
  const int wanted = opt.minimize ? std::max(1, opt.candidates) : 1;
  for (int a = 0; a < opt.max_attempts && complete < wanted; ++a) {
    std::mt19937_64 rng(opt.seed * 1000003ULL + static_cast<std::uint64_t>(a) * 7919ULL);
    Attempt at = embed_once(B, n, chiral, planar, rng);
    res.attempts = a + 1;
    if (at.xyz.empty() || !at.chirality_ok) continue;
    EmbedResult cur;
    cur.dg_energy = at.energy;
    copy_coordinates(st, s, at.xyz);
    if (opt.minimize) {
      Diagnostics quiet;  // per-attempt warnings are summarised below
      cur.uff = uff_minimize(st, sp, quiet, opt.uff);
    }
    if (opt.minimize && worst_sp3_carbon_angle(st, s) > kMaxSp3CarbonAngle) {
      ++inverted;
      continue;  // inverted geometry: not a real structure
    }
    score_stereo(cur);
    const std::size_t ok = cur.stereo_centres_ok + cur.stereo_bonds_ok;
    const double energy = opt.minimize ? cur.uff.energy_final : cur.dg_energy;
    if (best_xyz.empty() || ok > best_ok || (ok == best_ok && energy < best_energy)) {
      best = cur;
      best_ok = ok;
      best_energy = energy;
      best_xyz.resize(3 * n);
      for (std::uint32_t i = 0; i < n; ++i) {
        best_xyz[3 * i] = st.atoms.x[s.atom_begin + i];
        best_xyz[3 * i + 1] = st.atoms.y[s.atom_begin + i];
        best_xyz[3 * i + 2] = st.atoms.z[s.atom_begin + i];
      }
    }
    if (ok == cur.stereo_centres + cur.stereo_bonds) ++complete;
  }
  if (best_xyz.empty()) {
    if (inverted > 0)
      diag.error(DiagCode::EmbeddingFailed,
                 "the requested stereochemistry cannot exist in 3D: every structure that has it is inverted (an sp3 carbon "
                 "bent past 150 degrees), as with an impossible bridgehead configuration; check the @/@@ marks")
          .on_species(sp);
    else
      diag.error(DiagCode::EmbeddingFailed, "distance geometry found no coordinates with the requested chirality after " +
                                                std::to_string(res.attempts) + " attempts")
          .on_species(sp);
    return res;
  }
  const int attempts = res.attempts;
  res = best;
  res.attempts = attempts;
  copy_coordinates(st, s, best_xyz);
  if (opt.minimize && !res.uff.converged)
    diag.warn(DiagCode::MinimizationNotConverged, "UFF minimisation stopped after " + std::to_string(res.uff.iterations) +
                                                      " iterations (max |gradient| " + std::to_string(res.uff.grad_max) + ")")
        .on_species(sp);
  if (res.stereo_centres_ok != res.stereo_centres || res.stereo_bonds_ok != res.stereo_bonds)
    diag.warn(DiagCode::StereoLost, "3D structure reproduces " + std::to_string(res.stereo_centres_ok) + "/" +
                                        std::to_string(res.stereo_centres) + " stereocentres and " +
                                        std::to_string(res.stereo_bonds_ok) + "/" + std::to_string(res.stereo_bonds) +
                                        " stereo double bonds" +
                                        (inverted > 0 ? "; the full specification only gave inverted (impossible) "
                                                        "structures, so the closest real one is returned"
                                                      : ""))
        .on_species(sp);
  res.ok = true;
  return res;
}

std::int32_t embed_smiles(ChemicalSystemState& st, std::string_view smiles, Diagnostics& diag, EmbedResult* result,
                          const EmbedOptions& opt) {
  const std::int32_t sp = parse_smiles(st, smiles, diag);
  if (sp < 0) return -1;
  const std::int32_t sph = add_explicit_hydrogens(st, static_cast<std::uint32_t>(sp), diag);
  EmbedResult r = embed_species(st, static_cast<std::uint32_t>(sph), diag, opt);
  if (result) *result = r;
  return r.ok ? sph : -1;
}

// ---- stereo from coordinates -------------------------------------------------

Chirality parity_from_coordinates(const ChemicalSystemState& st, const TetrahedralStereo& t) {
  std::array<double, 12> p{};
  for (std::size_t k = 0; k < 4; ++k) {
    const std::uint32_t a = t.refs[k] == kLonePairRef ? t.atom : t.refs[k];
    if (a == kImplicitHRef) return Chirality::None;
    p[3 * k] = st.atoms.x[a];
    p[3 * k + 1] = st.atoms.y[a];
    p[3 * k + 2] = st.atoms.z[a];
  }
  const double V = volume(p.data(), 3, {0, 1, 2, 3}, nullptr);
  if (std::fabs(V) < 1e-6) return Chirality::None;
  return V < 0 ? Chirality::CounterClockwise : Chirality::Clockwise;
}

Chirality axial_parity_from_coordinates(const ChemicalSystemState& st, const AxialStereo& a) {
  TetrahedralStereo t;
  t.atom = a.begin;
  t.refs = a.refs;
  return parity_from_coordinates(st, t);
}

double dihedral_degrees(const ChemicalSystemState& st, std::uint32_t a, std::uint32_t b, std::uint32_t c, std::uint32_t d) {
  auto P = [&](std::uint32_t i) { return std::array<double, 3>{st.atoms.x[i], st.atoms.y[i], st.atoms.z[i]}; };
  auto sub = [](std::array<double, 3> u, std::array<double, 3> v) { return std::array<double, 3>{u[0] - v[0], u[1] - v[1], u[2] - v[2]}; };
  auto cr = [](std::array<double, 3> u, std::array<double, 3> v) {
    return std::array<double, 3>{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
  };
  auto dt = [](std::array<double, 3> u, std::array<double, 3> v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
  const auto b1 = sub(P(b), P(a)), b2 = sub(P(c), P(b)), b3 = sub(P(d), P(c));
  const auto m = cr(b1, b2), nn = cr(b2, b3);
  return std::atan2(std::sqrt(dt(b2, b2)) * dt(b1, nn), dt(m, nn)) / kDeg;
}

bool cis_from_coordinates(const ChemicalSystemState& st, const DoubleBondStereo& d) {
  return std::fabs(dihedral_degrees(st, d.ref_begin, st.bonds.begin[d.bond], st.bonds.end[d.bond], d.ref_end)) < 90.0;
}

// ---- output ------------------------------------------------------------------

std::string write_sdf(const ChemicalSystemState& st, std::uint32_t sp, const std::string& title) {
  const Species& s = st.species[sp];
  std::string out = (title.empty() ? s.canonical_smiles : title) + "\n  ChemLarp          3D\n\n";
  char buf[128];
  const bool chiral = std::any_of(st.tetrahedral.begin(), st.tetrahedral.end(),
                                  [&](const TetrahedralStereo& t) { return st.atoms.species[t.atom] == sp; });
  std::snprintf(buf, sizeof buf, "%3u%3u  0  0%3d  0  0  0  0  0999 V2000\n", s.atom_count, s.bond_count, chiral ? 1 : 0);
  out += buf;
  std::vector<std::pair<std::uint32_t, int>> charges, isotopes;
  std::vector<int> order_sum(s.atom_count, 0);
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    order_sum[st.bonds.begin[b] - s.atom_begin] += st.bonds.order[b];
    order_sum[st.bonds.end[b] - s.atom_begin] += st.bonds.order[b];
  }
  for (std::uint32_t i = 0; i < s.atom_count; ++i) {
    const std::uint32_t g = s.atom_begin + i;
    const int vvv = mdl_valence_field(st.atoms.element[g], st.atoms.formal_charge[g], order_sum[i], st.atoms.implicit_h[g]);
    std::snprintf(buf, sizeof buf, "%10.4f%10.4f%10.4f %-3.3s 0  0  0  0  0%3d  0  0  0  0  0  0\n", st.atoms.x[g],
                  st.atoms.y[g], st.atoms.z[g], std::string(element_symbol(st.atoms.element[g])).c_str(), vvv);
    out += buf;
    if (st.atoms.formal_charge[g] != 0) charges.push_back({i + 1, st.atoms.formal_charge[g]});
    if (st.atoms.isotope[g] != 0) isotopes.push_back({i + 1, st.atoms.isotope[g]});
  }
  for (std::uint32_t b = s.bond_begin; b < s.bond_begin + s.bond_count; ++b) {
    std::snprintf(buf, sizeof buf, "%3u%3u%3u  0\n", st.bonds.begin[b] - s.atom_begin + 1, st.bonds.end[b] - s.atom_begin + 1,
                  static_cast<unsigned>(st.bonds.order[b]));
    out += buf;
  }
  auto emit = [&](const char* tag, const std::vector<std::pair<std::uint32_t, int>>& v) {
    for (std::size_t k = 0; k < v.size(); k += 8) {
      const std::size_t m = std::min<std::size_t>(8, v.size() - k);
      std::snprintf(buf, sizeof buf, "M  %s%3zu", tag, m);
      out += buf;
      for (std::size_t q = 0; q < m; ++q) {
        std::snprintf(buf, sizeof buf, " %3u %3d", v[k + q].first, v[k + q].second);
        out += buf;
      }
      out += "\n";
    }
  };
  emit("CHG", charges);
  emit("ISO", isotopes);
  out += "M  END\n$$$$\n";
  return out;
}

std::string write_xyz(const ChemicalSystemState& st, std::uint32_t sp, const std::string& comment) {
  const Species& s = st.species[sp];
  std::string out = std::to_string(s.atom_count) + "\n" + (comment.empty() ? s.canonical_smiles : comment) + "\n";
  char buf[96];
  for (std::uint32_t i = 0; i < s.atom_count; ++i) {
    const std::uint32_t g = s.atom_begin + i;
    std::snprintf(buf, sizeof buf, "%-2s %14.8f %14.8f %14.8f\n", std::string(element_symbol(st.atoms.element[g])).c_str(),
                  st.atoms.x[g], st.atoms.y[g], st.atoms.z[g]);
    out += buf;
  }
  return out;
}

}  // namespace chemlarp
