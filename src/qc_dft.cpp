// Kohn-Sham exchange-correlation: molecular grid, functionals, and the XC potential.
//
// Functionals are written once as energy densities e(ra, rb, saa, sab, sbb) on
// dual numbers carrying the five first derivatives, so the potentials are exact
// derivatives of the published energy expressions.
#include <algorithm>
#include <array>
#include <cmath>
#include <map>

#include "chemlarp/parallel.hpp"
#include "qc_internal.hpp"

namespace chemlarp::qc {
namespace {

constexpr double kPi = 3.14159265358979323846;

// ---- dual numbers (value + 5 partial derivatives) -------------------------------------------

struct D {
  double v = 0.0;
  std::array<double, 5> d{};
  D() = default;
  D(double x) : v(x) {}  // NOLINT: constants promote implicitly
  static D var(double x, int i) {
    D r(x);
    r.d[static_cast<std::size_t>(i)] = 1.0;
    return r;
  }
};
inline D chain(const D& a, double f, double df) {
  D r(f);
  for (std::size_t i = 0; i < 5; ++i) r.d[i] = df * a.d[i];
  return r;
}
inline D operator+(const D& a, const D& b) { D r(a.v + b.v); for (std::size_t i = 0; i < 5; ++i) r.d[i] = a.d[i] + b.d[i]; return r; }
inline D operator-(const D& a, const D& b) { D r(a.v - b.v); for (std::size_t i = 0; i < 5; ++i) r.d[i] = a.d[i] - b.d[i]; return r; }
inline D operator-(const D& a) { D r(-a.v); for (std::size_t i = 0; i < 5; ++i) r.d[i] = -a.d[i]; return r; }
inline D operator*(const D& a, const D& b) { D r(a.v * b.v); for (std::size_t i = 0; i < 5; ++i) r.d[i] = a.d[i] * b.v + a.v * b.d[i]; return r; }
inline D operator/(const D& a, const D& b) {
  D r(a.v / b.v);
  const double ib2 = 1.0 / (b.v * b.v);
  for (std::size_t i = 0; i < 5; ++i) r.d[i] = (a.d[i] * b.v - a.v * b.d[i]) * ib2;
  return r;
}
inline D pow(const D& a, double p) { const double f = std::pow(a.v, p); return chain(a, f, p * f / a.v); }
inline D sqrt(const D& a) { const double f = std::sqrt(a.v); return chain(a, f, 0.5 / f); }
inline D exp(const D& a) { const double f = std::exp(a.v); return chain(a, f, f); }
inline D log(const D& a) { return chain(a, std::log(a.v), 1.0 / a.v); }
inline D atan(const D& a) { return chain(a, std::atan(a.v), 1.0 / (1.0 + a.v * a.v)); }
inline D asinh(const D& a) { return chain(a, std::asinh(a.v), 1.0 / std::sqrt(1.0 + a.v * a.v)); }

// ---- functionals ---------------------------------------------------------------------------

const double kCx = 0.75 * std::cbrt(6.0 / kPi);  // spin-resolved Slater: e = -kCx sum rho_s^(4/3)

D slater(const D& ra, const D& rb, bool a_on, bool b_on) {
  D e(0.0);
  if (a_on) e = e - kCx * pow(ra, 4.0 / 3.0);
  if (b_on) e = e - kCx * pow(rb, 4.0 / 3.0);
  return e;
}

// Becke 88 gradient correction (without the Slater part), per spin.
D b88_correction(const D& r, const D& s) {
  const double beta = 0.0042;
  const D r43 = pow(r, 4.0 / 3.0);
  const D x = sqrt(s) / r43;
  return -beta * r43 * x * x / (1.0 + 6.0 * beta * x * asinh(x));
}

// VWN-type correlation per particle at x = sqrt(rs) with parameters (A, x0, b, c).
D vwn_g(const D& x, double A, double x0, double b, double c) {
  const double Q = std::sqrt(4.0 * c - b * b);
  const D X = x * x + b * x + c;
  const double X0 = x0 * x0 + b * x0 + c;
  const D at = atan(Q / (2.0 * x + b));
  return A * (log(x * x / X) + (2.0 * b / Q) * at -
              (b * x0 / X0) * (log((x - x0) * (x - x0) / X) + (2.0 * (b + 2.0 * x0) / Q) * at));
}

const double kFz20 = 1.709920934161365617563962776245;  // f''(0)

D spin_f(const D& z) {
  return (pow(1.0 + z, 4.0 / 3.0) + pow(1.0 - z, 4.0 / 3.0) - 2.0) / (std::pow(2.0, 4.0 / 3.0) - 2.0);
}

struct VwnSet { double para[4], ferro[4], alpha[4]; };
// VWN5 (paper's fit to Ceperley-Alder) and VWN-RPA (fit to RPA), as (A, x0, b, c), hartree.
const VwnSet kVwn5 = {{0.0310907, -0.10498, 3.72744, 12.9352},
                      {0.01554535, -0.32500, 7.06042, 18.0578},
                      {-1.0 / (6.0 * kPi * kPi), -0.0047584, 1.13107, 13.0045}};
const VwnSet kVwnRpa = {{0.0310907, -0.409286, 13.0720, 42.7198},
                        {0.01554535, -0.743294, 20.1231, 101.578},
                        {-1.0 / (6.0 * kPi * kPi), -0.228344, 1.06835, 11.4813}};

// VWN correlation energy per volume. VWN5 uses the paper's spin interpolation
// with the spin stiffness alpha_c; the RPA fit is interpolated with f(zeta)
// alone (the "VWN1" form), as in B3LYP implementations (libxc, PySCF, Gaussian).
D vwn(const D& ra, const D& rb, const VwnSet& p, bool polarized) {
  const D rho = ra + rb;
  const D rs = pow(3.0 / (4.0 * kPi) / rho, 1.0 / 3.0);
  const D x = sqrt(rs);
  const D ep = vwn_g(x, p.para[0], p.para[1], p.para[2], p.para[3]);
  if (!polarized) return rho * ep;
  const D z = (ra - rb) / rho;
  const D ef = vwn_g(x, p.ferro[0], p.ferro[1], p.ferro[2], p.ferro[3]);
  const D f = spin_f(z);
  if (&p == &kVwnRpa) return rho * (ep + (ef - ep) * f);
  const D ac = vwn_g(x, p.alpha[0], p.alpha[1], p.alpha[2], p.alpha[3]);
  const D z4 = z * z * z * z;
  return rho * (ep + ac * f / kFz20 * (1.0 - z4) + (ef - ep) * f * z4);
}

// PW92 correlation per particle (modified-precision A values, as used with PBE).
D pw92_g(const D& rs, double A, double a1, double b1, double b2, double b3, double b4) {
  const D srs = sqrt(rs);
  return -2.0 * A * (1.0 + a1 * rs) * log(1.0 + 1.0 / (2.0 * A * (b1 * srs + b2 * rs + b3 * rs * srs + b4 * rs * rs)));
}
D pw92_eps(const D& rs, const D& z, bool polarized) {
  const D ep = pw92_g(rs, 0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294);
  if (!polarized) return ep;
  const D ef = pw92_g(rs, 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517);
  const D ac = -pw92_g(rs, 0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671);
  const D f = spin_f(z);
  const D z4 = z * z * z * z;
  return ep + ac * f / kFz20 * (1.0 - z4) + (ef - ep) * f * z4;
}

// PBE exchange for one spin channel (spin scaling: E_x[ra, rb] = (E_x[2ra] + E_x[2rb]) / 2).
D pbe_x_spin(const D& r, const D& s) {
  const double kappa = 0.804, mu = 0.2195149727645171;
  const D rho = 2.0 * r, sig = 4.0 * s;
  const D kf = pow(3.0 * kPi * kPi * rho, 1.0 / 3.0);
  const D ss = sig / (4.0 * kf * kf * rho * rho);  // s^2
  const D fx = 1.0 + kappa - kappa / (1.0 + mu * ss / kappa);
  const D ex_unif = -0.75 * std::cbrt(3.0 / kPi) * pow(rho, 4.0 / 3.0);
  return 0.5 * ex_unif * fx;
}

D pbe_c(const D& ra, const D& rb, const D& sigma, bool polarized) {
  const double beta = 0.06672455060314922, gamma = (1.0 - std::log(2.0)) / (kPi * kPi);
  const D rho = ra + rb;
  const D rs = pow(3.0 / (4.0 * kPi) / rho, 1.0 / 3.0);
  const D z = polarized ? (ra - rb) / rho : D(0.0);
  const D ec = pw92_eps(rs, z, polarized);
  const D phi = polarized ? 0.5 * (pow(1.0 + z, 2.0 / 3.0) + pow(1.0 - z, 2.0 / 3.0)) : D(1.0);
  const D kf = pow(3.0 * kPi * kPi * rho, 1.0 / 3.0);
  const D ks = sqrt(4.0 * kf / kPi);
  const D t2 = sigma / (4.0 * phi * phi * ks * ks * rho * rho);
  const D phi3 = phi * phi * phi;
  const D A = (beta / gamma) / (exp(-ec / (gamma * phi3)) - 1.0);
  const D at2 = A * t2;
  const D H = gamma * phi3 * log(1.0 + (beta / gamma) * t2 * (1.0 + at2) / (1.0 + at2 + at2 * at2));
  return rho * (ec + H);
}

// LYP correlation in the form of Miehlich, Savin, Stoll & Preuss (no Laplacian).
D lyp(const D& ra, const D& rb, const D& saa, const D& sab, const D& sbb) {
  const double a = 0.04918, b = 0.132, c = 0.2533, d = 0.349;
  const double cf = 0.3 * std::pow(3.0 * kPi * kPi, 2.0 / 3.0);
  const D rho = ra + rb;
  const D rm13 = pow(rho, -1.0 / 3.0);
  const D den = 1.0 + d * rm13;
  const D omega = exp(-c * rm13) / den * pow(rho, -11.0 / 3.0);
  const D delta = c * rm13 + d * rm13 / den;
  const D grad2 = saa + 2.0 * sab + sbb;
  const D t1 = -a * 4.0 / den * ra * rb / rho;
  const D inner = ra * rb * (std::pow(2.0, 11.0 / 3.0) * cf * (pow(ra, 8.0 / 3.0) + pow(rb, 8.0 / 3.0)) +
                             (47.0 / 18.0 - 7.0 / 18.0 * delta) * grad2 - (2.5 - delta / 18.0) * (saa + sbb) -
                             (delta - 11.0) / 9.0 * (ra / rho * saa + rb / rho * sbb)) -
                  2.0 / 3.0 * rho * rho * grad2 + (2.0 / 3.0 * rho * rho - ra * ra) * sbb +
                  (2.0 / 3.0 * rho * rho - rb * rb) * saa;
  return t1 - a * b * omega * inner;
}

D energy_density(Functional f, const D& ra, const D& rb, const D& saa, const D& sab, const D& sbb, bool polarized,
                 bool a_on, bool b_on) {
  switch (f) {
    case Functional::SVWN5:
      return slater(ra, rb, a_on, b_on) + vwn(ra, rb, kVwn5, polarized);
    case Functional::PBE: {
      D ex(0.0);
      if (a_on) ex = ex + pbe_x_spin(ra, saa);
      if (b_on) ex = ex + pbe_x_spin(rb, sbb);
      return ex + pbe_c(ra, rb, saa + 2.0 * sab + sbb, polarized);
    }
    case Functional::B3LYP: {
      D ex = 0.8 * slater(ra, rb, a_on, b_on);
      if (a_on) ex = ex + 0.72 * b88_correction(ra, saa);
      if (b_on) ex = ex + 0.72 * b88_correction(rb, sbb);
      return ex + 0.19 * vwn(ra, rb, kVwnRpa, polarized) + 0.81 * lyp(ra, rb, saa, sab, sbb);
    }
  }
  return D(0.0);
}

// ---- grid ------------------------------------------------------------------------------------

// Bragg-Slater radii (angstrom; H set to 0.35 as Becke recommends) for the partition size adjustment.
double bragg_radius(int z) {
  static const double r[] = {0, 0.35, 0.35, 1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50, 0.45,
                             1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 1.00};
  return z >= 1 && z <= 18 ? r[z] : 1.5;
}
// Treutler-Ahlrichs radial scaling xi (JCP 102 (1995) 346, table 1).
double ta_xi(int z) {
  static const double x[] = {0, 0.8, 0.9, 1.8, 1.4, 1.3, 1.1, 0.9, 0.9, 0.9, 0.9, 1.4, 1.3, 1.3, 1.2, 1.1, 1.0, 1.0, 1.0};
  return z >= 1 && z <= 18 ? x[z] : 1.0;
}

void gauss_legendre(int n, std::vector<double>& x, std::vector<double>& w) {
  x.resize(static_cast<std::size_t>(n));
  w.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    double z = std::cos(kPi * (i + 0.75) / (n + 0.5)), pp = 0.0;
    for (int it = 0; it < 100; ++it) {
      double p1 = 1.0, p2 = 0.0;
      for (int j = 1; j <= n; ++j) {
        const double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      const double dz = p1 / pp;
      z -= dz;
      if (std::fabs(dz) < 1e-15) break;
    }
    x[static_cast<std::size_t>(i)] = z;
    w[static_cast<std::size_t>(i)] = 2.0 / ((1.0 - z * z) * pp * pp);
  }
}

struct Point {
  double x, y, z, w;
};

std::vector<Point> molecular_grid(const Molecule& mol, int level, std::vector<std::pair<std::size_t, std::size_t>>& batches) {
  const int nat = static_cast<int>(mol.z.size());
  const int lv = std::clamp(level, 1, 3) - 1;
  const int ntheta = std::array<int, 3>{10, 18, 26}[static_cast<std::size_t>(lv)];
  // Angular product grids (Gauss-Legendre in cos(theta) x trapezoid in phi, exact for
  // spherical harmonics through degree 2*nt-1); smaller near the nuclei, where the
  // density is nearly spherical (pruning in the spirit of Treutler & Ahlrichs).
  auto make_sphere = [](int nt) {
    std::vector<double> ct, wt;
    gauss_legendre(nt, ct, wt);
    const int nphi = 2 * nt;
    std::vector<std::array<double, 4>> sp;
    for (int i = 0; i < nt; ++i)
      for (int j = 0; j < nphi; ++j) {
        const double ph = 2.0 * kPi * (j + 0.5) / nphi, st = std::sqrt(1.0 - ct[static_cast<std::size_t>(i)] * ct[static_cast<std::size_t>(i)]);
        sp.push_back({st * std::cos(ph), st * std::sin(ph), ct[static_cast<std::size_t>(i)], wt[static_cast<std::size_t>(i)] * 2.0 * kPi / nphi});
      }
    return sp;
  };
  const std::vector<std::array<double, 4>> sphere_full = make_sphere(ntheta);
  const std::vector<std::array<double, 4>> sphere_mid = make_sphere(std::max(6, (2 * ntheta) / 3));
  const std::vector<std::array<double, 4>> sphere_core = make_sphere(std::max(5, ntheta / 3));
  // inverse distances and Becke size adjustments
  std::vector<double> inv_r(static_cast<std::size_t>(nat) * nat, 0.0), adj(static_cast<std::size_t>(nat) * nat, 0.0);
  for (int a = 0; a < nat; ++a)
    for (int b = 0; b < nat; ++b) {
      if (a == b) continue;
      const double dx = mol.xyz[static_cast<std::size_t>(a)][0] - mol.xyz[static_cast<std::size_t>(b)][0];
      const double dy = mol.xyz[static_cast<std::size_t>(a)][1] - mol.xyz[static_cast<std::size_t>(b)][1];
      const double dz = mol.xyz[static_cast<std::size_t>(a)][2] - mol.xyz[static_cast<std::size_t>(b)][2];
      inv_r[static_cast<std::size_t>(a) * nat + b] = 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz);
      const double chi = bragg_radius(mol.z[static_cast<std::size_t>(a)]) / bragg_radius(mol.z[static_cast<std::size_t>(b)]);
      const double u = (chi - 1.0) / (chi + 1.0);
      adj[static_cast<std::size_t>(a) * nat + b] = std::clamp(u / (u * u - 1.0), -0.5, 0.5);
    }
  std::vector<Point> pts;
  std::vector<double> dist(static_cast<std::size_t>(nat)), cell(static_cast<std::size_t>(nat));
  for (int a = 0; a < nat; ++a) {
    const int z = mol.z[static_cast<std::size_t>(a)];
    const int row = z <= 2 ? 0 : z <= 10 ? 1 : 2;
    static const int nrad_tab[3][3] = {{40, 60, 90}, {55, 75, 110}, {65, 85, 125}};
    const int nr = nrad_tab[lv][row];
    const double xi = ta_xi(z), alpha = 0.6;
    for (int i = 1; i <= nr; ++i) {
      // Treutler-Ahlrichs M4 mapping of Chebyshev (second kind) nodes
      const double th = i * kPi / (nr + 1), x = std::cos(th);
      const double r = xi / std::log(2.0) * std::pow(1.0 + x, alpha) * std::log(2.0 / (1.0 - x));
      const double dr = xi / std::log(2.0) * (alpha * std::pow(1.0 + x, alpha - 1.0) * std::log(2.0 / (1.0 - x)) + std::pow(1.0 + x, alpha) / (1.0 - x));
      const double wr = kPi / (nr + 1) * std::sin(th) * dr * r * r;
      const auto& sphere = lv == 2 ? sphere_full : r < 0.25 * xi ? sphere_core : r < 0.6 * xi ? sphere_mid : sphere_full;
      for (const auto& s : sphere) {
        Point p{mol.xyz[static_cast<std::size_t>(a)][0] + r * s[0], mol.xyz[static_cast<std::size_t>(a)][1] + r * s[1],
                mol.xyz[static_cast<std::size_t>(a)][2] + r * s[2], wr * s[3]};
        // Becke partition weight of atom a at p
        for (int b = 0; b < nat; ++b) {
          const double dx = p.x - mol.xyz[static_cast<std::size_t>(b)][0], dy = p.y - mol.xyz[static_cast<std::size_t>(b)][1], dz = p.z - mol.xyz[static_cast<std::size_t>(b)][2];
          dist[static_cast<std::size_t>(b)] = std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        double total = 0.0;
        for (int b = 0; b < nat; ++b) {
          double pb = 1.0;
          for (int c = 0; c < nat && pb > 0.0; ++c) {
            if (c == b) continue;
            double mu = (dist[static_cast<std::size_t>(b)] - dist[static_cast<std::size_t>(c)]) * inv_r[static_cast<std::size_t>(b) * nat + c];
            mu += adj[static_cast<std::size_t>(b) * nat + c] * (1.0 - mu * mu);
            for (int k = 0; k < 3; ++k) mu = 1.5 * mu - 0.5 * mu * mu * mu;
            pb *= 0.5 * (1.0 - mu);
          }
          cell[static_cast<std::size_t>(b)] = pb;
          total += pb;
        }
        p.w *= total > 0.0 ? cell[static_cast<std::size_t>(a)] / total : 0.0;
        if (p.w > 1e-15) pts.push_back(p);
      }
    }
  }
  // Batches: points grouped by 2-bohr boxes (at most 256 points each), so a batch
  // only touches basis functions near it.
  const double box = 2.0;
  std::map<std::array<long, 3>, std::vector<std::size_t>> cells;
  for (std::size_t i = 0; i < pts.size(); ++i)
    cells[{static_cast<long>(std::floor(pts[i].x / box)), static_cast<long>(std::floor(pts[i].y / box)),
           static_cast<long>(std::floor(pts[i].z / box))}].push_back(i);
  std::vector<Point> sorted;
  sorted.reserve(pts.size());
  for (const auto& [key, idx] : cells)
    for (std::size_t k = 0; k < idx.size(); k += 256) {
      const std::size_t start = sorted.size();
      for (std::size_t q = k; q < std::min(idx.size(), k + 256); ++q) sorted.push_back(pts[idx[q]]);
      batches.push_back({start, sorted.size()});
    }
  return sorted;
}

}  // namespace

double exact_exchange_fraction(Functional f) { return f == Functional::B3LYP ? 0.2 : 0.0; }

struct XcIntegrator::Impl {
  const Basis* basis = nullptr;
  Functional f = Functional::SVWN5;
  std::vector<Point> pts;
  std::vector<std::pair<std::size_t, std::size_t>> batches;
  std::vector<double> shell_cut;  // radius beyond which a shell is below 1e-11
};

XcIntegrator::XcIntegrator(const Molecule& mol, const Basis& basis, Functional f, int grid_level) : impl_(new Impl) {
  impl_->basis = &basis;
  impl_->f = f;
  impl_->pts = molecular_grid(mol, grid_level, impl_->batches);
  for (const Shell& s : basis.shells) {
    double rc = 0.0;
    for (std::size_t k = 0; k < s.exps.size(); ++k) {
      const double c = std::fabs(s.coefs[k]) * 10.0;  // margin for the polynomial factor
      if (c > 1e-11) rc = std::max(rc, std::sqrt(std::log(c / 1e-11) / s.exps[k]));
    }
    impl_->shell_cut.push_back(rc + 1.0);
  }
}

XcIntegrator::~XcIntegrator() { delete impl_; }
std::size_t XcIntegrator::points() const { return impl_->pts.size(); }

XcResult XcIntegrator::evaluate(const Mat& pa, const Mat& pb, bool restricted) const {
  const Basis& B = *impl_->basis;
  const std::size_t n = static_cast<std::size_t>(B.n);
  const bool gga = impl_->f != Functional::SVWN5;
  XcResult out;
  out.va.assign(n * n, 0.0);
  if (!restricted) out.vb.assign(n * n, 0.0);
  // Batches are split into a fixed number of contiguous parts (it depends on the grid and the basis size,
  // not on the thread count); each part accumulates its own potential and energy, and the parts are added
  // in order, so the result is the same for any number of threads.
  struct Acc {
    Mat va, vb;
    double energy = 0.0, electrons = 0.0;
  };
  const std::size_t nb = impl_->batches.size();
  const std::size_t cap = std::max<std::size_t>(1, static_cast<std::size_t>(2e7 / static_cast<double>(n * n + 1)));
  const std::size_t parts = std::max<std::size_t>(1, std::min<std::size_t>({nb, 32, cap}));
  std::vector<Acc> acc(parts);
  parallel_for(parts, [&](std::size_t part) {
    Acc& A = acc[part];
    A.va.assign(n * n, 0.0);
    if (!restricted) A.vb.assign(n * n, 0.0);
    std::vector<int> fn;
    // function-major arrays: value[f * np + p]
    std::vector<double> phi, gx, gy, gz, xa, xb, wa, wb, sub;
    auto batch = [&](std::size_t b0, std::size_t b1) {
    const std::size_t np = b1 - b0;
    double cx = 0, cy = 0, cz = 0;
    for (std::size_t i = b0; i < b1; ++i) { cx += impl_->pts[i].x; cy += impl_->pts[i].y; cz += impl_->pts[i].z; }
    cx /= static_cast<double>(np); cy /= static_cast<double>(np); cz /= static_cast<double>(np);
    double rad = 0.0;
    for (std::size_t i = b0; i < b1; ++i)
      rad = std::max(rad, std::hypot(impl_->pts[i].x - cx, impl_->pts[i].y - cy, impl_->pts[i].z - cz));
    fn.clear();
    std::vector<std::pair<const Shell*, std::size_t>> shells;
    for (std::size_t s = 0; s < B.shells.size(); ++s) {
      const Shell& sh = B.shells[s];
      if (std::hypot(sh.center[0] - cx, sh.center[1] - cy, sh.center[2] - cz) > impl_->shell_cut[s] + rad) continue;
      shells.push_back({&sh, fn.size()});
      const int ncomp = (sh.L + 1) * (sh.L + 2) / 2;
      for (int k = 0; k < ncomp; ++k) fn.push_back(sh.first + k);
    }
    if (fn.empty()) return;
    std::size_t nf = fn.size();
    phi.assign(nf * np, 0.0);
    if (gga) { gx.assign(nf * np, 0.0); gy.assign(nf * np, 0.0); gz.assign(nf * np, 0.0); }
    for (const auto& [sh, off] : shells) {
      const auto comps = cartesian_components(sh->L);
      for (std::size_t ip = 0; ip < np; ++ip) {
        const Point& P = impl_->pts[b0 + ip];
        const double dx = P.x - sh->center[0], dy = P.y - sh->center[1], dz = P.z - sh->center[2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        double rad0 = 0.0, rad1 = 0.0;  // sum c e^{-a r^2}, sum -2 a c e^{-a r^2}
        for (std::size_t k = 0; k < sh->exps.size(); ++k) {
          const double e = sh->coefs[k] * std::exp(-sh->exps[k] * r2);
          rad0 += e;
          rad1 -= 2.0 * sh->exps[k] * e;
        }
        if (std::fabs(rad0) < 1e-16 && std::fabs(rad1) < 1e-16) continue;
        for (std::size_t c = 0; c < comps.size(); ++c) {
          const auto& pw = comps[c];
          const double sc = B.scale[static_cast<std::size_t>(sh->first) + c];
          const double xl = pw[0] ? std::pow(dx, pw[0]) : 1.0, yl = pw[1] ? std::pow(dy, pw[1]) : 1.0, zl = pw[2] ? std::pow(dz, pw[2]) : 1.0;
          const double ang = xl * yl * zl;
          const std::size_t idx = (off + c) * np + ip;
          phi[idx] = sc * ang * rad0;
          if (gga) {
            const double dxl = pw[0] ? pw[0] * std::pow(dx, pw[0] - 1) : 0.0;
            const double dyl = pw[1] ? pw[1] * std::pow(dy, pw[1] - 1) : 0.0;
            const double dzl = pw[2] ? pw[2] * std::pow(dz, pw[2] - 1) : 0.0;
            gx[idx] = sc * (dxl * yl * zl * rad0 + ang * rad1 * dx);
            gy[idx] = sc * (xl * dyl * zl * rad0 + ang * rad1 * dy);
            gz[idx] = sc * (xl * yl * dzl * rad0 + ang * rad1 * dz);
          }
        }
      }
    }
    // Keep only functions that matter somewhere in this batch.
    {
      std::size_t keep = 0;
      for (std::size_t f = 0; f < nf; ++f) {
        double mx = 0.0;
        for (std::size_t ip = 0; ip < np; ++ip) {
          const std::size_t idx = f * np + ip;
          mx = std::max(mx, std::fabs(phi[idx]));
          if (gga) mx = std::max({mx, std::fabs(gx[idx]), std::fabs(gy[idx]), std::fabs(gz[idx])});
        }
        if (mx < 1e-10) continue;
        if (keep != f) {
          fn[keep] = fn[f];
          std::copy(phi.begin() + static_cast<std::ptrdiff_t>(f * np), phi.begin() + static_cast<std::ptrdiff_t>((f + 1) * np), phi.begin() + static_cast<std::ptrdiff_t>(keep * np));
          if (gga) {
            std::copy(gx.begin() + static_cast<std::ptrdiff_t>(f * np), gx.begin() + static_cast<std::ptrdiff_t>((f + 1) * np), gx.begin() + static_cast<std::ptrdiff_t>(keep * np));
            std::copy(gy.begin() + static_cast<std::ptrdiff_t>(f * np), gy.begin() + static_cast<std::ptrdiff_t>((f + 1) * np), gy.begin() + static_cast<std::ptrdiff_t>(keep * np));
            std::copy(gz.begin() + static_cast<std::ptrdiff_t>(f * np), gz.begin() + static_cast<std::ptrdiff_t>((f + 1) * np), gz.begin() + static_cast<std::ptrdiff_t>(keep * np));
          }
        }
        ++keep;
      }
      if (keep == 0) return;
      fn.resize(keep);
      nf = keep;
    }
    // X = P_sub phi (nf x np), then rho = sum_f phi X and grad rho = 2 sum_f grad(phi) X
    auto halfdens = [&](const Mat& p, std::vector<double>& x) {
      x.assign(nf * np, 0.0);
      for (std::size_t f = 0; f < nf; ++f) {
        const double* prow = &p[static_cast<std::size_t>(fn[f]) * n];
        double* xf = &x[f * np];
        for (std::size_t g = 0; g < nf; ++g) {
          const double pfg = prow[fn[g]];
          if (pfg == 0.0) continue;
          const double* pg = &phi[g * np];
          for (std::size_t ip = 0; ip < np; ++ip) xf[ip] += pfg * pg[ip];
        }
      }
    };
    halfdens(pa, xa);
    if (!restricted) halfdens(pb, xb);
    std::vector<double> ra(np, 0.0), rb(np, 0.0);
    std::vector<std::array<double, 3>> ga(np, {0, 0, 0}), gb(np, {0, 0, 0});
    for (std::size_t f = 0; f < nf; ++f)
      for (std::size_t ip = 0; ip < np; ++ip) {
        const std::size_t idx = f * np + ip;
        ra[ip] += phi[idx] * xa[idx];
        if (gga) { ga[ip][0] += 2.0 * gx[idx] * xa[idx]; ga[ip][1] += 2.0 * gy[idx] * xa[idx]; ga[ip][2] += 2.0 * gz[idx] * xa[idx]; }
        if (!restricted) {
          rb[ip] += phi[idx] * xb[idx];
          if (gga) { gb[ip][0] += 2.0 * gx[idx] * xb[idx]; gb[ip][1] += 2.0 * gy[idx] * xb[idx]; gb[ip][2] += 2.0 * gz[idx] * xb[idx]; }
        }
      }
    if (restricted) { rb = ra; gb = ga; }
    // functional at each point: v_rho and the gradient factors
    std::vector<double> vra(np, 0.0), vrb(np, 0.0);
    std::vector<std::array<double, 3>> fa(np, {0, 0, 0}), fb(np, {0, 0, 0});
    for (std::size_t ip = 0; ip < np; ++ip) {
      const double w = impl_->pts[b0 + ip].w;
      A.electrons += w * (ra[ip] + rb[ip]);
      const bool a_on = ra[ip] > 1e-14, b_on = rb[ip] > 1e-14;
      if (!a_on && !b_on) continue;
      const double saa = ga[ip][0] * ga[ip][0] + ga[ip][1] * ga[ip][1] + ga[ip][2] * ga[ip][2];
      const double sbb = gb[ip][0] * gb[ip][0] + gb[ip][1] * gb[ip][1] + gb[ip][2] * gb[ip][2];
      const double sab = ga[ip][0] * gb[ip][0] + ga[ip][1] * gb[ip][1] + ga[ip][2] * gb[ip][2];
      const D e = energy_density(impl_->f, D::var(std::max(ra[ip], 1e-30), 0), D::var(std::max(rb[ip], 1e-30), 1),
                                 D::var(std::max(saa, 1e-60), 2), D::var(sab, 3), D::var(std::max(sbb, 1e-60), 4),
                                 !restricted, a_on, b_on);
      if (!std::isfinite(e.v)) continue;
      A.energy += w * e.v;
      vra[ip] = w * e.d[0];
      vrb[ip] = w * e.d[1];
      if (gga)
        for (int k = 0; k < 3; ++k) {
          fa[ip][static_cast<std::size_t>(k)] = w * (2.0 * e.d[2] * ga[ip][static_cast<std::size_t>(k)] + e.d[3] * gb[ip][static_cast<std::size_t>(k)]);
          fb[ip][static_cast<std::size_t>(k)] = w * (2.0 * e.d[4] * gb[ip][static_cast<std::size_t>(k)] + e.d[3] * ga[ip][static_cast<std::size_t>(k)]);
        }
    }
    // W[f][p] = 0.5 v_rho phi + f . grad(phi); V_sub = W phi^T + phi W^T
    auto potential = [&](const std::vector<double>& vr, const std::vector<std::array<double, 3>>& fv, std::vector<double>& wm, Mat& v) {
      wm.assign(nf * np, 0.0);
      for (std::size_t f = 0; f < nf; ++f)
        for (std::size_t ip = 0; ip < np; ++ip) {
          const std::size_t idx = f * np + ip;
          double t = 0.5 * vr[ip] * phi[idx];
          if (gga) t += fv[ip][0] * gx[idx] + fv[ip][1] * gy[idx] + fv[ip][2] * gz[idx];
          wm[idx] = t;
        }
      sub.assign(nf * nf, 0.0);
      for (std::size_t f = 0; f < nf; ++f) {
        const double* wf = &wm[f * np];
        for (std::size_t g = 0; g < nf; ++g) {
          const double* pg = &phi[g * np];
          double acc = 0.0;
          for (std::size_t ip = 0; ip < np; ++ip) acc += wf[ip] * pg[ip];
          sub[f * nf + g] = acc;
        }
      }
      for (std::size_t f = 0; f < nf; ++f)
        for (std::size_t g = 0; g < nf; ++g)
          v[static_cast<std::size_t>(fn[f]) * n + static_cast<std::size_t>(fn[g])] += sub[f * nf + g] + sub[g * nf + f];
    };
    potential(vra, fa, wa, A.va);
    if (!restricted) potential(vrb, fb, wb, A.vb);
      };
    for (std::size_t bi = partition_begin(nb, parts, part); bi < partition_begin(nb, parts, part + 1); ++bi)
      batch(impl_->batches[bi].first, impl_->batches[bi].second);
  });
  for (const Acc& A : acc) {
    out.energy += A.energy;
    out.electrons += A.electrons;
    for (std::size_t k = 0; k < n * n; ++k) out.va[k] += A.va[k];
    if (!restricted)
      for (std::size_t k = 0; k < n * n; ++k) out.vb[k] += A.vb[k];
  }
  if (restricted) out.vb = out.va;
  return out;
}

}  // namespace chemlarp::qc
