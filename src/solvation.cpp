// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
#include "chemlarp/solvation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <numbers>
#include <random>

#include "chemlarp/conformer.hpp"
#include "chemlarp/elements.hpp"
#include "chemlarp/ewald.hpp"
#include "chemlarp/json.hpp"
#include "chemlarp/mmff.hpp"
#include "chemlarp/parallel.hpp"
#include "chemlarp/uff.hpp"

namespace chemlarp {

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kAcc = 4.184e-4;         // (kcal/mol/Å)/amu -> Å/fs²
constexpr double kBoltz = 0.0019872043;   // kcal/(mol K)
constexpr double kAtm = 1.4583972e-5;     // 1 atm in kcal/(mol Å³)
constexpr double kAmuPerA3 = 1.66053907;  // g/cm³ per (amu/Å³)

// TIP3P (Jorgensen et al. 1983).
constexpr double kOH = 0.9572, kHOH = 104.52 * kPi / 180.0;
constexpr double kQO = -0.834, kQH = 0.417;
const double kSigO = std::pow(582000.0 / 595.0, 1.0 / 6.0);  // 3.15066 Å
const double kEpsO = 595.0 * 595.0 / (4.0 * 582000.0);       // 0.15207 kcal/mol
constexpr double kMO = 15.9994, kMH = 1.008;

using V3 = std::array<double, 3>;

struct System {
  std::size_t n = 0, ns = 0;  // atoms; the first ns are the solute
  std::vector<double> x, v, m, q, sig, eps;
  std::vector<int> mol;                                // molecule of each atom (solute: 0)
  std::vector<std::vector<std::uint32_t>> molecules;   // atoms of each molecule
  std::vector<std::array<std::uint32_t, 3>> waters;   // O, H, H (the solute too when it is a water)
  bool solute_water = false;
  std::shared_ptr<const MmffForceField> mmff;          // solute intramolecular
  V3 L{};
  std::size_t n_oxygen = 0;                            // water oxygens other than the solute
};

inline double wrap(double d, double L) { return d - L * std::round(d / L); }

// ---- Energies and forces ---------------------------------------------------------------------------------

struct Parts {
  double coul = 0.0;    // real + reciprocal + self + excluded corrections
  double lj_ww = 0.0, lj_sw = 0.0, lrc = 0.0, intra = 0.0;
  double total() const { return coul + lj_ww + lj_sw + lrc + intra; }
};

// erfc(βr)/r tabulated as a cubic Hermite spline in r (values and slopes exact at the nodes), so that the
// force is the exact derivative of the interpolated energy and energy conservation is kept.
struct ErfcTable {
  double r0 = 0.6, h = 0.005;
  std::vector<double> f, df;  // erfc(βr)/r and its r-derivative at r0 + k h
  void build(double beta, double rmax) {
    const std::size_t n = static_cast<std::size_t>((rmax - r0) / h) + 3;
    f.resize(n);
    df.resize(n);
    const double twob = 2.0 * beta / std::sqrt(kPi);
    for (std::size_t k = 0; k < n; ++k) {
      const double r = r0 + double(k) * h, e = std::erfc(beta * r);
      f[k] = e / r;
      df[k] = -(e / r + twob * std::exp(-beta * beta * r * r)) / r;
    }
  }
  // Value and derivative at r (r >= r0).
  void eval(double r, double& v, double& dv) const {
    const double t = (r - r0) / h;
    const std::size_t k = static_cast<std::size_t>(t);
    const double u = t - double(k), u2 = u * u, u3 = u2 * u;
    const double h00 = 2 * u3 - 3 * u2 + 1, h10 = u3 - 2 * u2 + u, h01 = -2 * u3 + 3 * u2, h11 = u3 - u2;
    v = h00 * f[k] + h10 * h * df[k] + h01 * f[k + 1] + h11 * h * df[k + 1];
    const double d00 = 6 * u2 - 6 * u, d10 = 3 * u2 - 4 * u + 1, d01 = -6 * u2 + 6 * u, d11 = 3 * u2 - 2 * u;
    dv = (d00 * f[k] + d01 * f[k + 1]) / h + d10 * df[k] + d11 * df[k + 1];
  }
};

struct Engine {
  const SolvationOptions& o;
  double beta, rc2, rl2;
  Pme pme;
  ErfcTable table;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> pairs;
  std::vector<double> x_list;  // positions at the last list build
  V3 L_list{};

  Engine(const SolvationOptions& opt, const System& s)
      : o(opt), beta(ewald_beta(opt.cutoff, opt.ewald_tolerance)), rc2(opt.cutoff * opt.cutoff),
        rl2((opt.cutoff + opt.skin) * (opt.cutoff + opt.skin)), pme(s.L, beta, opt.pme_order, opt.pme_spacing) {
    table.build(beta, opt.cutoff + 0.1);
  }

  void build_list(const System& s) {
    pairs.clear();
    for (std::uint32_t i = 0; i < s.n; ++i)
      for (std::uint32_t j = i + 1; j < s.n; ++j) {
        if (s.mol[i] == s.mol[j]) continue;
        if (s.q[i] == 0.0 && s.q[j] == 0.0 && (s.eps[i] == 0.0 || s.eps[j] == 0.0)) continue;
        double r2 = 0.0;
        for (int d = 0; d < 3; ++d) {
          const double dx = wrap(s.x[3 * i + d] - s.x[3 * j + d], s.L[d]);
          r2 += dx * dx;
        }
        if (r2 < rl2) pairs.push_back({i, j});
      }
    x_list = s.x;
    L_list = s.L;
  }

  void maybe_rebuild(const System& s) {
    if (x_list.size() != s.x.size() || s.L != L_list) return build_list(s);
    double m2 = 0.0;
    for (std::size_t i = 0; i < s.n; ++i) {
      double d2 = 0.0;
      for (int d = 0; d < 3; ++d) {
        const double dx = s.x[3 * i + d] - x_list[3 * i + d];
        d2 += dx * dx;
      }
      m2 = std::max(m2, d2);
    }
    if (4.0 * m2 > o.skin * o.skin) build_list(s);
  }

  // Soft-core LJ (Beutler et al. 1994) with σ, ε; returns U, sets ff (= −dU/dr / r) and dU/dλ.
  static double softcore(double r2, double sig, double eps, double lam, double alpha, double& ff, double& dudl) {
    const double s2 = sig * sig;
    const double R6 = (r2 / s2) * (r2 / s2) * (r2 / s2);
    const double D = std::max(alpha * (1.0 - lam) + R6, 1e-300);
    const double iD = 1.0 / D, iD2 = iD * iD;
    ff = 4.0 * eps * lam * (2.0 * iD2 * iD - iD2) * 6.0 * R6 / r2;
    dudl = 4.0 * eps * (iD2 - iD) + 4.0 * eps * lam * alpha * (2.0 * iD2 * iD - iD2);
    return 4.0 * eps * lam * (iD2 - iD);
  }

  // Long-range dispersion correction: water–water at full strength, solute–water scaled by λ_vdw.
  double lrc(const System& s, double lv, double* d_lv = nullptr) const {
    const double V = s.L[0] * s.L[1] * s.L[2], rc = o.cutoff;
    const double rho = double(s.n_oxygen) / V;
    auto tail = [&](double sig, double eps) {  // ∫_rc^∞ 4πr² 4ε[(σ/r)^12 − (σ/r)^6] dr per unit density
      const double sr3 = std::pow(sig / rc, 3), sr9 = sr3 * sr3 * sr3;
      return 16.0 * kPi * eps * sig * sig * sig * (sr9 / 9.0 - sr3 / 3.0);
    };
    // Water–water: N_O(N_O − 1)/2 pairs ≈ N_O ρ / 2 per molecule.
    const double nww = double(s.n_oxygen);
    double e = 0.5 * nww * (nww / V) * tail(kSigO, kEpsO);
    double esw = 0.0;
    for (std::size_t i = 0; i < s.ns; ++i)
      if (s.eps[i] > 0.0) esw += rho * tail(0.5 * (s.sig[i] + kSigO), std::sqrt(s.eps[i] * kEpsO));
    if (d_lv) *d_lv = esw;
    return e + lv * esw;
  }

  // Total energy (and forces when F) at (λ_elec, λ_vdw). Deterministic: fixed partitions summed in order.
  double compute(System& s, double le, double lv, std::vector<double>* F, Parts* parts) {
    maybe_rebuild(s);
    const std::size_t N = s.n;
    std::vector<double> q(s.q);
    for (std::size_t i = 0; i < s.ns; ++i) q[i] *= le;
    const double alpha = o.softcore_alpha;
    const std::size_t P = 32;
    const std::size_t extra = 4;  // coul, lj_ww, lj_sw, (unused)
    const std::size_t width = (F ? 3 * N : 0) + extra;
    const double b = beta, twob = 2.0 * beta / std::sqrt(kPi);
    std::vector<double> acc = parallel_sum(P, width, [&](std::size_t part, std::vector<double>& out) {
      const std::size_t a0 = partition_begin(pairs.size(), P, part), a1 = partition_begin(pairs.size(), P, part + 1);
      double* f = F ? out.data() : nullptr;
      double& ec = out[width - 4];
      double& eww = out[width - 3];
      double& esw = out[width - 2];
      for (std::size_t k = a0; k < a1; ++k) {
        const auto [i, j] = pairs[k];
        double d[3], r2 = 0.0;
        for (int c = 0; c < 3; ++c) {
          d[c] = wrap(s.x[3 * i + c] - s.x[3 * j + c], s.L[c]);
          r2 += d[c] * d[c];
        }
        if (r2 >= rc2) continue;
        double ff = 0.0;
        const double qq = kCoulomb * q[i] * q[j];
        if (qq != 0.0) {
          const double r = std::sqrt(r2);
          if (r >= table.r0) {
            double v, dv;
            table.eval(r, v, dv);
            ec += qq * v;
            ff -= qq * dv / r;
          } else {
            const double erfcv = std::erfc(b * r);
            ec += qq * erfcv / r;
            ff += qq * (erfcv / r + twob * std::exp(-b * b * r2)) / r2;
          }
        }
        if (s.eps[i] > 0.0 && s.eps[j] > 0.0) {
          const bool sw = (i < s.ns) != (j < s.ns);
          const double sig = 0.5 * (s.sig[i] + s.sig[j]), eps = std::sqrt(s.eps[i] * s.eps[j]);
          if (sw) {
            double f2, dl;
            esw += softcore(r2, sig, eps, lv, alpha, f2, dl);
            ff += f2;
          } else {
            const double s2 = sig * sig / r2, s6 = s2 * s2 * s2, s12 = s6 * s6;
            eww += 4.0 * eps * (s12 - s6);
            ff += 24.0 * eps * (2.0 * s12 - s6) / r2;
          }
        }
        if (f)
          for (int c = 0; c < 3; ++c) {
            f[3 * i + c] += ff * d[c];
            f[3 * j + c] -= ff * d[c];
          }
      }
    });
    Parts p;
    p.coul = acc[width - 4];
    p.lj_ww = acc[width - 3];
    p.lj_sw = acc[width - 2];
    if (F) F->assign(acc.begin(), acc.begin() + static_cast<std::ptrdiff_t>(3 * N));
    // Excluded pairs inside molecules: remove the reciprocal-space part.
    for (const auto& mol : s.molecules)
      for (std::size_t a = 0; a < mol.size(); ++a)
        for (std::size_t c = a + 1; c < mol.size(); ++c) {
          const std::uint32_t i = mol[a], j = mol[c];
          const double qq = kCoulomb * q[i] * q[j];
          if (qq == 0.0) continue;
          double d[3], r2 = 0.0;
          for (int k = 0; k < 3; ++k) {
            d[k] = wrap(s.x[3 * i + k] - s.x[3 * j + k], s.L[k]);
            r2 += d[k] * d[k];
          }
          const double r = std::sqrt(r2), erfv = std::erf(b * r);
          p.coul -= qq * erfv / r;
          if (F) {
            // −d/dr(−qq erf/r) / r
            const double ff = -qq * (erfv / r - twob * std::exp(-b * b * r2)) / r2;
            for (int k = 0; k < 3; ++k) {
              (*F)[3 * i + k] += ff * d[k];
              (*F)[3 * j + k] -= ff * d[k];
            }
          }
        }
    p.coul += ewald_self(b, q);
    pme.set_box(s.L);
    p.coul += pme.energy(s.x, q, F);
    p.lrc = lrc(s, lv);
    if (s.mmff) {
      std::vector<double> xs(s.x.begin(), s.x.begin() + static_cast<std::ptrdiff_t>(3 * s.ns)), g;
      p.intra = F ? s.mmff->energy_gradient(xs, g) : s.mmff->energy(xs);
      if (F)
        for (std::size_t k = 0; k < 3 * s.ns; ++k) (*F)[k] -= g[k];
    }
    if (parts) *parts = p;
    return p.total();
  }

  // Electrostatic energy only, at λ_elec (for the quadratic fit).
  double coulomb(System& s, double le) {
    Parts p;
    compute(s, le, 1.0, nullptr, &p);
    return p.coul;
  }

  // Solute–water soft-core energy at λ_vdw (and its λ-derivative), including the long-range part.
  double vdw_sw(const System& s, double lv, double& dudl) const {
    double e = 0.0;
    dudl = 0.0;
    for (const auto& [i, j] : pairs) {
      if ((i < s.ns) == (j < s.ns)) continue;
      if (s.eps[i] == 0.0 || s.eps[j] == 0.0) continue;
      double r2 = 0.0;
      for (int c = 0; c < 3; ++c) {
        const double dx = wrap(s.x[3 * i + c] - s.x[3 * j + c], s.L[c]);
        r2 += dx * dx;
      }
      if (r2 >= rc2) continue;
      double ff, dl;
      e += softcore(r2, 0.5 * (s.sig[i] + s.sig[j]), std::sqrt(s.eps[i] * s.eps[j]), lv, o.softcore_alpha, ff, dl);
      dudl += dl;
    }
    double dl_lrc = 0.0;
    lrc(s, lv, &dl_lrc);
    dudl += dl_lrc;  // the long-range part is linear in λ_vdw
    return e + lv * dl_lrc;
  }
};

// ---- Constraints ------------------------------------------------------------------------------------------

// SETTLE (Miyamoto & Kollman 1992): new positions of a rigid TIP3P water from its constrained old positions
// (b) and unconstrained new ones (a), in the form used by GROMACS (csettle).
bool settle(const double* bO, const double* bH1, const double* bH2, double* aO, double* aH1, double* aH2, double mO,
            double mH) {
  const double dOH = kOH, dHH = 2.0 * kOH * std::sin(0.5 * kHOH);
  const double wohh = mO + 2.0 * mH, rc = 0.5 * dHH;
  const double ra = 2.0 * mH * std::sqrt(dOH * dOH - rc * rc) / wohh, rb = std::sqrt(dOH * dOH - rc * rc) - ra;
  V3 xb0, xc0, com, xa1, xb1, xc1;
  for (int d = 0; d < 3; ++d) {
    xb0[d] = bH1[d] - bO[d];
    xc0[d] = bH2[d] - bO[d];
    com[d] = (aO[d] * mO + (aH1[d] + aH2[d]) * mH) / wohh;
    xa1[d] = aO[d] - com[d];
    xb1[d] = aH1[d] - com[d];
    xc1[d] = aH2[d] - com[d];
  }
  auto cross = [](const V3& a, const V3& b) { return V3{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; };
  auto dot = [](const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
  V3 ez = cross(xb0, xc0), ex = cross(xa1, ez), ey = cross(ez, ex);
  for (V3* e : {&ex, &ey, &ez}) {
    const double l = std::sqrt(dot(*e, *e));
    if (!(l > 0)) return false;
    for (double& c : *e) c /= l;
  }
  const double xb0d = dot(ex, xb0), yb0d = dot(ey, xb0), xc0d = dot(ex, xc0), yc0d = dot(ey, xc0);
  const double za1d = dot(ez, xa1);
  const double xb1d = dot(ex, xb1), yb1d = dot(ey, xb1), zb1d = dot(ez, xb1);
  const double xc1d = dot(ex, xc1), yc1d = dot(ey, xc1), zc1d = dot(ez, xc1);
  const double sinphi = za1d / ra;
  const double t = 1.0 - sinphi * sinphi;
  if (t <= 0.0) return false;
  const double cosphi = std::sqrt(t);
  const double sinpsi = (zb1d - zc1d) / (dHH * cosphi);
  const double t2 = 1.0 - sinpsi * sinpsi;
  if (t2 <= 0.0) return false;
  const double cospsi = std::sqrt(t2);
  const double ya2d = ra * cosphi, xb2d = -rc * cospsi;
  const double tt1 = -rb * cosphi, tt2 = rc * sinpsi * sinphi;
  const double yb2d = tt1 - tt2, yc2d = tt1 + tt2;
  const double alpha = xb2d * (xb0d - xc0d) + yb0d * yb2d + yc0d * yc2d;
  const double beta = xb2d * (yc0d - yb0d) + xb0d * yb2d + xc0d * yc2d;
  const double gamma = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;
  const double al2be2 = alpha * alpha + beta * beta;
  const double t3 = al2be2 - gamma * gamma;
  if (t3 < 0.0 || al2be2 <= 0.0) return false;
  const double sinthe = (alpha * gamma - beta * std::sqrt(t3)) / al2be2;
  const double costhe = std::sqrt(std::max(0.0, 1.0 - sinthe * sinthe));
  const double xa3d = -ya2d * sinthe, ya3d = ya2d * costhe, za3d = za1d;
  const double xb3d = xb2d * costhe - yb2d * sinthe, yb3d = xb2d * sinthe + yb2d * costhe, zb3d = zb1d;
  const double xc3d = -xb2d * costhe - yc2d * sinthe, yc3d = -xb2d * sinthe + yc2d * costhe, zc3d = zc1d;
  for (int d = 0; d < 3; ++d) {
    aO[d] = com[d] + xa3d * ex[d] + ya3d * ey[d] + za3d * ez[d];
    aH1[d] = com[d] + xb3d * ex[d] + yb3d * ey[d] + zb3d * ez[d];
    aH2[d] = com[d] + xc3d * ex[d] + yc3d * ey[d] + zc3d * ez[d];
  }
  return true;
}

// Removes the velocity components that would change the three distances of a rigid water (exact 3×3 solve).
void constrain_velocities(System& s, const std::array<std::uint32_t, 3>& w) {
  const std::uint32_t at[3][2] = {{w[0], w[1]}, {w[0], w[2]}, {w[1], w[2]}};
  double dd[3][3], rhs[3];
  std::array<V3, 3> dvec;
  for (int k = 0; k < 3; ++k)
    for (int c = 0; c < 3; ++c) dvec[k][c] = wrap(s.x[3 * at[k][1] + c] - s.x[3 * at[k][0] + c], s.L[c]);
  auto sgn = [&](std::uint32_t atom, int l) { return atom == at[l][1] ? 1.0 : atom == at[l][0] ? -1.0 : 0.0; };
  for (int k = 0; k < 3; ++k) {
    const std::uint32_t a = at[k][0], b = at[k][1];
    rhs[k] = 0.0;
    for (int c = 0; c < 3; ++c) rhs[k] -= (s.v[3 * b + c] - s.v[3 * a + c]) * dvec[k][c];
    for (int l = 0; l < 3; ++l) {
      double dl = 0.0;
      for (int c = 0; c < 3; ++c) dl += dvec[l][c] * dvec[k][c];
      dd[k][l] = dl * (sgn(b, l) / s.m[b] - sgn(a, l) / s.m[a]);
    }
  }
  // Solve dd τ = rhs (Cramer).
  auto det3 = [](const double m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  };
  const double D = det3(dd);
  if (std::fabs(D) < 1e-300) return;
  double tau[3];
  for (int c = 0; c < 3; ++c) {
    double mm[3][3];
    for (int r = 0; r < 3; ++r)
      for (int q = 0; q < 3; ++q) mm[r][q] = q == c ? rhs[r] : dd[r][q];
    tau[c] = det3(mm) / D;
  }
  for (int l = 0; l < 3; ++l)
    for (int c = 0; c < 3; ++c) {
      s.v[3 * at[l][0] + c] -= tau[l] * dvec[l][c] / s.m[at[l][0]];
      s.v[3 * at[l][1] + c] += tau[l] * dvec[l][c] / s.m[at[l][1]];
    }
}

// ---- Building the box -------------------------------------------------------------------------------------

struct Rng {
  std::mt19937_64 g;
  explicit Rng(std::uint64_t seed) : g(seed) {}
  double uniform() { return (static_cast<double>(g() >> 11) + 0.5) * 0x1.0p-53; }
  double normal() {
    const double u = uniform(), v = uniform();
    return std::sqrt(-2.0 * std::log(u)) * std::cos(2.0 * kPi * v);
  }
};

// TIP3P in its own frame: O at the origin, the H's in the xy plane.
std::array<V3, 3> water_geometry() {
  const double hx = kOH * std::sin(0.5 * kHOH), hy = kOH * std::cos(0.5 * kHOH);
  return {V3{0, 0, 0}, V3{hx, hy, 0}, V3{-hx, hy, 0}};
}

V3 rotate(const V3& p, const double R[3][3]) {
  return {R[0][0] * p[0] + R[0][1] * p[1] + R[0][2] * p[2], R[1][0] * p[0] + R[1][1] * p[1] + R[1][2] * p[2],
          R[2][0] * p[0] + R[2][1] * p[1] + R[2][2] * p[2]};
}

void random_rotation(Rng& rng, double R[3][3]) {
  // Uniform quaternion.
  double q[4];
  double n = 0.0;
  for (double& c : q) c = rng.normal(), n += c * c;
  n = std::sqrt(n);
  for (double& c : q) c /= n;
  const double w = q[0], x = q[1], y = q[2], z = q[3];
  const double M[3][3] = {{1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)},
                          {2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)},
                          {2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)}};
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) R[a][b] = M[a][b];
}

void add_water(System& s, const V3& centre, Rng& rng, bool solute) {
  double R[3][3];
  random_rotation(rng, R);
  const auto g = water_geometry();
  const int mol = static_cast<int>(s.molecules.size());
  std::array<std::uint32_t, 3> idx{};
  std::vector<std::uint32_t> atoms;
  for (int k = 0; k < 3; ++k) {
    const V3 p = rotate(g[static_cast<std::size_t>(k)], R);
    for (int d = 0; d < 3; ++d) s.x.push_back(centre[d] + p[d]);
    s.m.push_back(k == 0 ? kMO : kMH);
    s.q.push_back(k == 0 ? kQO : kQH);
    s.sig.push_back(k == 0 ? kSigO : 0.0);
    s.eps.push_back(k == 0 ? kEpsO : 0.0);
    s.mol.push_back(mol);
    idx[static_cast<std::size_t>(k)] = static_cast<std::uint32_t>(s.n);
    atoms.push_back(static_cast<std::uint32_t>(s.n));
    ++s.n;
  }
  s.waters.push_back(idx);
  s.molecules.push_back(atoms);
  if (solute) s.ns = 3, s.solute_water = true;
  else ++s.n_oxygen;
}

// Solute from SMILES: 3D structure (MMFF minimum), MMFF charges and intramolecular terms, UFF Lennard-Jones.
bool add_solute(System& s, std::string_view smiles, Diagnostics& diag) {
  ChemicalSystemState st;
  const std::int32_t sp = embed_smiles(st, smiles, diag);
  if (sp < 0) return false;
  const auto spi = static_cast<std::uint32_t>(sp);
  mmff_minimize(st, spi, diag, MmffVariant::Mmff94, LbfgsOptions{5000, 10, 1e-4, 1e-15, 0.3});
  auto ff = MmffForceField::build(st, spi, diag, MmffVariant::Mmff94);
  if (!ff) return false;
  const auto types = assign_uff_types(st, spi, diag);
  const Species& spc = st.species[spi];
  if (types.size() != spc.atom_count) {
    diag.error(DiagCode::InvalidInput, "no UFF types for the solute");
    return false;
  }
  double net = 0.0;
  for (double c : ff->charges()) net += c;
  if (std::fabs(net) > 1e-6) {
    diag.error(DiagCode::InvalidInput, "charged solutes are not supported (net charge " + std::to_string(net) + ")");
    return false;
  }
  V3 c{};
  for (std::uint32_t a = 0; a < spc.atom_count; ++a) {
    c[0] += st.atoms.x[spc.atom_begin + a] / spc.atom_count;
    c[1] += st.atoms.y[spc.atom_begin + a] / spc.atom_count;
    c[2] += st.atoms.z[spc.atom_begin + a] / spc.atom_count;
  }
  std::vector<std::uint32_t> atoms;
  for (std::uint32_t a = 0; a < spc.atom_count; ++a) {
    const std::uint32_t g = spc.atom_begin + a;
    s.x.push_back(st.atoms.x[g] - c[0] + 0.5 * s.L[0]);
    s.x.push_back(st.atoms.y[g] - c[1] + 0.5 * s.L[1]);
    s.x.push_back(st.atoms.z[g] - c[2] + 0.5 * s.L[2]);
    const ElementInfo* e = find_element_by_z(st.atoms.element[g]);
    s.m.push_back(e ? e->mass : 12.0);
    s.q.push_back(ff->charges()[a]);
    const auto& up = uff_parameters(types[a]);
    s.sig.push_back(up.x1 / std::pow(2.0, 1.0 / 6.0));
    s.eps.push_back(up.D1);
    s.mol.push_back(0);
    atoms.push_back(static_cast<std::uint32_t>(s.n));
    ++s.n;
  }
  // Hydrogen-mass repartitioning: hydrogens 3.024 amu, the difference taken from their bonded atom.
  for (std::uint32_t a = 0; a < spc.atom_count; ++a) {
    const std::uint32_t g = spc.atom_begin + a;
    if (st.atoms.element[g] != 1) continue;
    const auto nb = st.adjacency.neighbors(g);
    if (nb.empty()) continue;
    const std::uint32_t heavy = nb[0] - spc.atom_begin;
    const double dm = 3.024 - s.m[a];
    if (s.m[heavy] - dm < 1.0) continue;
    s.m[a] += dm;
    s.m[heavy] -= dm;
  }
  s.molecules.push_back(atoms);
  s.ns = spc.atom_count;
  s.mmff = std::make_shared<const MmffForceField>(std::move(*ff));
  return true;
}

// Box of `waters` TIP3P molecules (plus the solute, already placed at the centre when ns > 0).
void fill_water(System& s, int waters, Rng& rng) {
  const int n_side = static_cast<int>(std::ceil(std::cbrt(double(waters) * 1.3)));
  const double a = s.L[0] / n_side;
  std::vector<V3> sites;
  for (int i = 0; i < n_side; ++i)
    for (int j = 0; j < n_side; ++j)
      for (int k = 0; k < n_side; ++k) {
        const V3 p{(i + 0.5) * a, (j + 0.5) * a, (k + 0.5) * a};
        bool clash = false;
        for (std::size_t t = 0; t < s.ns && !clash; ++t) {
          double r2 = 0.0;
          for (int d = 0; d < 3; ++d) {
            const double dx = wrap(p[d] - s.x[3 * t + d], s.L[d]);
            r2 += dx * dx;
          }
          clash = r2 < 2.8 * 2.8;
        }
        if (!clash) sites.push_back(p);
      }
  // A random subset of the free sites.
  for (std::size_t i = sites.size(); i > 1; --i) std::swap(sites[i - 1], sites[static_cast<std::size_t>(rng.uniform() * double(i))]);
  sites.resize(std::min<std::size_t>(sites.size(), static_cast<std::size_t>(waters)));
  for (const auto& p : sites) add_water(s, p, rng, false);
  s.v.assign(3 * s.n, 0.0);
}

// ---- Dynamics ---------------------------------------------------------------------------------------------

struct Dynamics {
  System& s;
  Engine& e;
  const SolvationOptions& o;
  Rng rng;
  double kT;
  double le = 1.0, lv = 1.0;
  std::vector<double> F;
  double E = 0.0;
  Parts parts;
  int moves = 0, accepted = 0;
  double vmax;  // barostat step, Å³

  Dynamics(System& sys, Engine& eng, const SolvationOptions& opt, std::uint64_t seed)
      : s(sys), e(eng), o(opt), rng(seed), kT(kBoltz * opt.temperature_K),
        vmax(0.01 * sys.L[0] * sys.L[1] * sys.L[2]) {}

  void forces() { E = e.compute(s, le, lv, &F, &parts); }

  void thermalize() {
    for (std::size_t i = 0; i < s.n; ++i)
      for (int c = 0; c < 3; ++c) s.v[3 * i + c] = std::sqrt(kT * kAcc / s.m[i]) * rng.normal();
    for (const auto& w : s.waters) constrain_velocities(s, w);
  }

  void kick(double h) {
    for (std::size_t i = 0; i < s.n; ++i)
      for (int c = 0; c < 3; ++c) s.v[3 * i + c] += h * F[3 * i + c] * kAcc / s.m[i];
    for (const auto& w : s.waters) constrain_velocities(s, w);
  }

  void drift(double h) {
    const std::vector<double> old = s.x;
    for (std::size_t k = 0; k < 3 * s.n; ++k) s.x[k] += h * s.v[k];
    for (const auto& w : s.waters) {
      double* a = &s.x[3 * w[0]];
      double* b = &s.x[3 * w[1]];
      double* c = &s.x[3 * w[2]];
      if (!settle(&old[3 * w[0]], &old[3 * w[1]], &old[3 * w[2]], a, b, c, kMO, kMH)) continue;
      for (std::uint32_t at : w)
        for (int d = 0; d < 3; ++d) s.v[3 * at + d] = (s.x[3 * at + d] - old[3 * at + d]) / h;
    }
  }

  void thermostat(double h) {
    const double c1 = std::exp(-o.friction_ps * 1e-3 * h), c2 = std::sqrt(1.0 - c1 * c1);
    for (std::size_t i = 0; i < s.n; ++i)
      for (int c = 0; c < 3; ++c) s.v[3 * i + c] = c1 * s.v[3 * i + c] + c2 * std::sqrt(kT * kAcc / s.m[i]) * rng.normal();
    for (const auto& w : s.waters) constrain_velocities(s, w);
  }

  bool step(double dt) {
    kick(0.5 * dt);
    drift(0.5 * dt);
    thermostat(dt);
    drift(0.5 * dt);
    forces();
    kick(0.5 * dt);
    return std::isfinite(E);
  }

  // Monte Carlo volume move (isotropic; molecule centres scaled, molecules kept rigid).
  void barostat() {
    if (o.pressure_atm <= 0.0) return;
    const double V = s.L[0] * s.L[1] * s.L[2];
    const double Vn = V + vmax * (2.0 * rng.uniform() - 1.0);
    if (!(Vn > 0)) return;
    const double f = std::cbrt(Vn / V);
    if (0.5 * s.L[0] * f <= o.cutoff) return;
    const std::vector<double> xold = s.x;
    const V3 Lold = s.L;
    for (const auto& mol : s.molecules) {
      V3 c{};
      double mt = 0.0;
      for (std::uint32_t a : mol) {
        for (int d = 0; d < 3; ++d) c[d] += s.m[a] * s.x[3 * a + d];
        mt += s.m[a];
      }
      for (int d = 0; d < 3; ++d) {
        const double shift = (f - 1.0) * c[d] / mt;
        for (std::uint32_t a : mol) s.x[3 * a + d] += shift;
      }
    }
    for (double& l : s.L) l *= f;
    Parts pn;
    const double En = e.compute(s, le, lv, nullptr, &pn);
    const double nmol = double(s.molecules.size());
    const double w = -(En - E + o.pressure_atm * kAtm * (Vn - V)) / kT + nmol * std::log(Vn / V);
    ++moves;
    if (std::log(rng.uniform()) < w) {
      ++accepted;
      forces();
    } else {
      s.x = xold;
      s.L = Lold;
      e.build_list(s);
    }
    if (moves % 10 == 0) {
      const double rate = double(accepted) / double(moves);
      if (rate < 0.25) vmax *= 0.9;
      else if (rate > 0.75) vmax *= 1.1;
    }
  }

  // Steepest descent with the waters kept rigid, to remove the clashes of the starting lattice.
  void minimize(int steps) {
    forces();
    double h = 0.02;
    for (int it = 0; it < steps; ++it) {
      double fmax = 0.0;
      for (double f : F) fmax = std::max(fmax, std::fabs(f));
      if (fmax < 1e-6) break;
      const std::vector<double> old = s.x;
      const double Eold = E;
      for (std::size_t k = 0; k < 3 * s.n; ++k) s.x[k] += h * F[k] / fmax;
      for (const auto& w : s.waters)
        settle(&old[3 * w[0]], &old[3 * w[1]], &old[3 * w[2]], &s.x[3 * w[0]], &s.x[3 * w[1]], &s.x[3 * w[2]], kMO, kMH);
      forces();
      if (E < Eold) h = std::min(0.2, h * 1.2);
      else {
        s.x = old;
        forces();
        h *= 0.5;
      }
    }
  }

  double temperature() const {
    double k = 0.0;
    for (std::size_t i = 0; i < s.n; ++i)
      for (int c = 0; c < 3; ++c) k += 0.5 * s.m[i] * s.v[3 * i + c] * s.v[3 * i + c] / kAcc;
    const double dof = 3.0 * double(s.n) - 3.0 * double(s.waters.size());
    return 2.0 * k / (dof * kBoltz);
  }

  // Equilibration path: minimisation, a short run at 0.5 fs, then the full step.
  bool equilibrate(double ps) {
    minimize(300);
    thermalize();
    forces();
    for (int i = 0; i < 400; ++i)
      if (!step(0.5)) return false;
    const int n = static_cast<int>(ps * 1000.0 / o.dt_fs);
    for (int i = 0; i < n; ++i) {
      if (!step(o.dt_fs)) return false;
      if (o.barostat_every > 0 && (i + 1) % o.barostat_every == 0) barostat();
    }
    return true;
  }
};

System make_box(const SolvationOptions& o, int waters, double extra_mass) {
  System s;
  const double mass = waters * (kMO + 2 * kMH) + extra_mass;
  const double Lb = std::cbrt(mass * kAmuPerA3 / o.density);
  s.L = {Lb, Lb, Lb};
  return s;
}

std::string jn(double v) {
  if (!std::isfinite(v)) return "null";
  char b[32];
  std::snprintf(b, sizeof b, "%.8g", v);
  return b;
}

}  // namespace

// ---- Pure water --------------------------------------------------------------------------------------------

std::optional<WaterBoxResult> simulate_water(const SolvationOptions& o, double ps, Diagnostics& diag) {
  if (o.waters < 64 || o.waters > 4000 || !(ps > 0)) {
    diag.error(DiagCode::InvalidInput, "water box: 64-4000 molecules and a positive run length");
    return std::nullopt;
  }
  System s = make_box(o, o.waters, 0.0);
  if (0.5 * s.L[0] <= o.cutoff) {
    diag.error(DiagCode::InvalidInput, "the box is too small for the cutoff: use more water");
    return std::nullopt;
  }
  Rng rng(o.seed);
  fill_water(s, o.waters, rng);
  Engine e(o, s);
  e.build_list(s);
  Dynamics dyn(s, e, o, o.seed * 7919 + 1);
  if (!dyn.equilibrate(o.equilibration_ps)) {
    diag.error(DiagCode::NotConverged, "the water box became unstable during equilibration");
    return std::nullopt;
  }
  WaterBoxResult r;
  r.waters = static_cast<int>(s.waters.size());
  const int n = static_cast<int>(ps * 1000.0 / o.dt_fs), every = std::max(1, static_cast<int>(0.1 * 1000.0 / o.dt_fs));
  std::vector<double> rho, epot, temps, boxes;
  const double mass = double(s.waters.size()) * (kMO + 2 * kMH);
  const double rmax = 0.5 * s.L[0] * 0.98, dr = 0.05;
  std::vector<double> hist(static_cast<std::size_t>(rmax / dr) + 1, 0.0);
  int frames = 0;
  double vol_acc = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!dyn.step(o.dt_fs)) {
      diag.error(DiagCode::NotConverged, "the water box became unstable");
      return std::nullopt;
    }
    if (o.barostat_every > 0 && (i + 1) % o.barostat_every == 0) dyn.barostat();
    if ((i + 1) % every == 0) {
      const double V = s.L[0] * s.L[1] * s.L[2];
      rho.push_back(mass * kAmuPerA3 / V);
      epot.push_back(dyn.E / double(s.waters.size()));
      temps.push_back(dyn.temperature());
      boxes.push_back(s.L[0]);
      // O–O pairs for g(r).
      for (std::size_t a = 0; a < s.waters.size(); ++a)
        for (std::size_t b2 = a + 1; b2 < s.waters.size(); ++b2) {
          double r2 = 0.0;
          for (int d = 0; d < 3; ++d) {
            const double dx = wrap(s.x[3 * s.waters[a][0] + d] - s.x[3 * s.waters[b2][0] + d], s.L[d]);
            r2 += dx * dx;
          }
          const double rr = std::sqrt(r2);
          if (rr < rmax) hist[static_cast<std::size_t>(rr / dr)] += 2.0;
        }
      vol_acc += V;
      ++frames;
    }
  }
  auto mean_err = [](const std::vector<double>& y, double& err) {
    double m = 0.0;
    for (double q : y) m += q;
    m /= double(y.size());
    const double g = statistical_inefficiency(y);
    double var = 0.0;
    for (double q : y) var += (q - m) * (q - m);
    var /= double(y.size() > 1 ? y.size() - 1 : 1);
    err = std::sqrt(var * g / double(y.size()));
    return m;
  };
  double terr = 0.0, berr = 0.0;
  r.density = mean_err(rho, r.density_error);
  r.potential_per_molecule = mean_err(epot, r.potential_error);
  r.temperature = mean_err(temps, terr);
  r.box = mean_err(boxes, berr);
  const double Nw = double(s.waters.size()), Vavg = vol_acc / frames;
  for (std::size_t k = 0; k < hist.size(); ++k) {
    const double r0 = k * dr, r1 = r0 + dr;
    const double shell = 4.0 / 3.0 * kPi * (r1 * r1 * r1 - r0 * r0 * r0);
    const double g = hist[k] / (double(frames) * Nw * shell * (Nw - 1) / Vavg);
    r.g_oo.push_back({r0 + 0.5 * dr, g});
  }
  r.acceptance = dyn.moves ? double(dyn.accepted) / dyn.moves : 0.0;
  return r;
}

// ---- Hydration free energy --------------------------------------------------------------------------------

std::optional<HydrationResult> hydration_free_energy(std::string_view solute, const SolvationOptions& o, Diagnostics& diag) {
  std::vector<std::pair<double, double>> states = o.states;
  if (states.empty())
    states = {{1, 1},   {0.75, 1}, {0.5, 1},  {0.25, 1}, {0, 1},    {0, 0.9},  {0, 0.8},  {0, 0.7},  {0, 0.6},  {0, 0.5},
              {0, 0.4}, {0, 0.35}, {0, 0.3},  {0, 0.25}, {0, 0.2},  {0, 0.15}, {0, 0.1},  {0, 0.05}, {0, 0}};
  if (states.front() != std::pair<double, double>{1, 1} || states.back() != std::pair<double, double>{0, 0}) {
    diag.error(DiagCode::InvalidInput, "the windows must run from (1, 1) coupled to (0, 0) decoupled");
    return std::nullopt;
  }
  System s = make_box(o, o.waters, 0.0);
  Rng rng(o.seed);
  const bool is_water = solute == "water" || solute == "TIP3P";
  if (is_water) {
    add_water(s, {0.5 * s.L[0], 0.5 * s.L[1], 0.5 * s.L[2]}, rng, true);
  } else if (!add_solute(s, solute, diag)) {
    return std::nullopt;
  }
  if (0.5 * s.L[0] <= o.cutoff + 1.0) {
    diag.error(DiagCode::InvalidInput, "the box is too small for the cutoff: use more water");
    return std::nullopt;
  }
  fill_water(s, o.waters, rng);
  Engine e0(o, s);
  e0.build_list(s);
  {
    Dynamics dyn(s, e0, o, o.seed * 7919 + 3);
    if (!dyn.equilibrate(o.equilibration_ps)) {
      diag.error(DiagCode::NotConverged, "the solvated system became unstable during equilibration");
      return std::nullopt;
    }
  }
  const std::size_t K = states.size();
  const double kT = kBoltz * o.temperature_K;
  struct WindowData {
    std::vector<std::vector<double>> u;  // [state][sample], reduced
    std::vector<double> dudl_e, dudl_v;
    bool ok = true;
  };
  std::vector<WindowData> data(K);
  const int prod = static_cast<int>(o.production_ps * 1000.0 / o.dt_fs);
  const int eq = static_cast<int>(o.equilibration_ps * 1000.0 / o.dt_fs);
  const int every = std::max(1, static_cast<int>(o.sample_every_ps * 1000.0 / o.dt_fs));
  parallel_for(K, [&](std::size_t w) {
    System sw = s;
    Engine e(o, sw);
    e.build_list(sw);
    Dynamics dyn(sw, e, o, o.seed * 104729ULL + w * 7919ULL + 11);
    dyn.le = states[w].first;
    dyn.lv = states[w].second;
    dyn.thermalize();
    dyn.forces();
    WindowData& D = data[w];
    D.u.assign(K, {});
    for (int i = 0; i < eq + prod; ++i) {
      if (!dyn.step(o.dt_fs)) {
        D.ok = false;
        return;
      }
      if (o.barostat_every > 0 && (i + 1) % o.barostat_every == 0) dyn.barostat();
      if (i >= eq && (i - eq + 1) % every == 0) {
        // Electrostatics: E(λ) = a + bλ + cλ², from λ = 0, 1, −1.
        const double c0 = e.coulomb(sw, 0.0), c1 = e.coulomb(sw, 1.0), cm = e.coulomb(sw, -1.0);
        const double bq = 0.5 * (c1 - cm), cq = 0.5 * (c1 + cm) - c0;
        for (std::size_t k = 0; k < K; ++k) {
          double dl;
          const double ev = e.vdw_sw(sw, states[k].second, dl);
          const double le = states[k].first;
          D.u[k].push_back((c0 + bq * le + cq * le * le + ev) / kT);
        }
        double dlv;
        e.vdw_sw(sw, states[w].second, dlv);
        D.dudl_e.push_back(bq + 2.0 * cq * states[w].first);
        D.dudl_v.push_back(dlv);
      }
    }
  });
  for (std::size_t w = 0; w < K; ++w)
    if (!data[w].ok) {
      diag.error(DiagCode::NotConverged, "a window became unstable (λ_elec " + std::to_string(states[w].first) + ", λ_vdw " +
                                             std::to_string(states[w].second) + ")");
      return std::nullopt;
    }
  // Subsample, then MBAR, BAR and TI.
  HydrationResult res;
  res.solute = std::string(solute);
  res.solute_atoms = static_cast<int>(s.ns);
  res.waters = static_cast<int>(s.n_oxygen);
  res.box = s.L[0];
  std::vector<std::vector<std::size_t>> keep(K);
  std::vector<std::size_t> nk(K);
  for (std::size_t w = 0; w < K; ++w) {
    const WindowData& D = data[w];
    const std::size_t n = D.dudl_e.size();
    double g = 1.0;
    for (std::size_t nb : {w > 0 ? w - 1 : w, w + 1 < K ? w + 1 : w}) {
      if (nb == w) continue;
      std::vector<double> du(n);
      for (std::size_t i = 0; i < n; ++i) du[i] = D.u[nb][i] - D.u[w][i];
      g = std::max(g, statistical_inefficiency(du));
    }
    keep[w] = subsample_indices(n, g);
    nk[w] = keep[w].size();
    HydrationWindow hw;
    hw.lambda_elec = states[w].first;
    hw.lambda_vdw = states[w].second;
    hw.samples = keep[w].size();
    hw.inefficiency = g;
    for (std::size_t i : keep[w]) hw.dudl_elec += D.dudl_e[i] / double(keep[w].size()), hw.dudl_vdw += D.dudl_v[i] / double(keep[w].size());
    res.windows.push_back(hw);
  }
  std::vector<std::vector<double>> u(K);
  for (std::size_t w = 0; w < K; ++w)
    for (std::size_t k = 0; k < K; ++k)
      for (std::size_t i : keep[w]) u[k].push_back(data[w].u[k][i]);
  const MbarResult m = mbar(u, nk);
  // ΔG_hyd = G(coupled) − G(decoupled) = −(f_last − f_first) kT.
  res.dg_mbar = {-m.df[0][K - 1] * kT, m.ddf[0][K - 1] * kT};
  std::size_t mid = 0;
  for (std::size_t k = 0; k < K; ++k)
    if (states[k].first == 0.0 && states[k].second == 1.0) mid = k;
  res.dg_elec = -m.df[0][mid] * kT;
  res.dg_vdw = -m.df[mid][K - 1] * kT;
  res.min_neighbour_overlap = 1.0;
  for (std::size_t k = 0; k + 1 < K; ++k) res.min_neighbour_overlap = std::min(res.min_neighbour_overlap, m.overlap[k][k + 1]);
  double bs = 0.0, bv = 0.0;
  for (std::size_t k = 0; k + 1 < K; ++k) {
    std::vector<double> wf, wr;
    for (std::size_t i : keep[k]) wf.push_back(data[k].u[k + 1][i] - data[k].u[k][i]);
    for (std::size_t i : keep[k + 1]) wr.push_back(data[k + 1].u[k][i] - data[k + 1].u[k + 1][i]);
    const FreeEnergyEstimate b = bar(wf, wr);
    bs += b.df;
    bv += b.error * b.error;
  }
  res.dg_bar = {-bs * kT, std::sqrt(bv) * kT};
  // TI along the path: d G = ⟨∂U/∂λ⟩ dλ for each segment (λ_elec on the first part, λ_vdw on the second).
  double ti = 0.0;
  for (std::size_t k = 0; k + 1 < K; ++k) {
    const double dle = states[k + 1].first - states[k].first, dlv = states[k + 1].second - states[k].second;
    ti += 0.5 * dle * (res.windows[k].dudl_elec + res.windows[k + 1].dudl_elec) +
          0.5 * dlv * (res.windows[k].dudl_vdw + res.windows[k + 1].dudl_vdw);
  }
  res.dg_ti = {-ti, 0.0};
  if (res.min_neighbour_overlap < 0.03)
    diag.warn(DiagCode::NotConverged, "neighbouring windows overlap poorly; add windows where the overlap is small");
  return res;
}

std::string water_box_json(const SolvationOptions& o, double ps) {
  Diagnostics d;
  const auto r = simulate_water(o, ps, d);
  std::string s = std::string("{\"ok\":") + (r ? "true" : "false");
  if (r) {
    s += ",\"waters\":" + std::to_string(r->waters) + ",\"box\":" + jn(r->box) + ",\"density\":" + jn(r->density) +
         ",\"density_error\":" + jn(r->density_error) + ",\"potential_per_molecule\":" + jn(r->potential_per_molecule) +
         ",\"potential_error\":" + jn(r->potential_error) + ",\"temperature\":" + jn(r->temperature) +
         ",\"barostat_acceptance\":" + jn(r->acceptance) + ",\"g_oo\":[";
    for (std::size_t i = 0; i < r->g_oo.size(); ++i) s += std::string(i ? "," : "") + "[" + jn(r->g_oo[i].first) + "," + jn(r->g_oo[i].second) + "]";
    s += "]";
  }
  return s + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

std::string hydration_json(std::string_view solute, const SolvationOptions& o) {
  Diagnostics d;
  const auto r = hydration_free_energy(solute, o, d);
  std::string s = std::string("{\"ok\":") + (r ? "true" : "false");
  if (r) {
    auto est = [](const FreeEnergyEstimate& e) { return "{\"dg\":" + jn(e.df) + ",\"error\":" + jn(e.error) + "}"; };
    s += ",\"solute\":" + json_escape(r->solute) + ",\"solute_atoms\":" + std::to_string(r->solute_atoms) + ",\"waters\":" +
         std::to_string(r->waters) + ",\"box\":" + jn(r->box) + ",\"dg_hyd\":{\"mbar\":" + est(r->dg_mbar) + ",\"bar\":" + est(r->dg_bar) +
         ",\"ti\":" + est(r->dg_ti) + "},\"dg_elec\":" + jn(r->dg_elec) + ",\"dg_vdw\":" + jn(r->dg_vdw) +
         ",\"min_neighbour_overlap\":" + jn(r->min_neighbour_overlap) + ",\"windows\":[";
    for (std::size_t i = 0; i < r->windows.size(); ++i) {
      const auto& w = r->windows[i];
      s += std::string(i ? "," : "") + "{\"lambda_elec\":" + jn(w.lambda_elec) + ",\"lambda_vdw\":" + jn(w.lambda_vdw) +
           ",\"samples\":" + std::to_string(w.samples) + ",\"inefficiency\":" + jn(w.inefficiency) + ",\"dudl_elec\":" +
           jn(w.dudl_elec) + ",\"dudl_vdw\":" + jn(w.dudl_vdw) + "}";
    }
    s += "]";
  }
  return s + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

}  // namespace chemlarp
