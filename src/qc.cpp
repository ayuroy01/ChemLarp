#include "chemlarp/parallel.hpp"
#include "chemlarp/qc.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>

#include "chemlarp/elements.hpp"
#include "qc_internal.hpp"

namespace chemlarp {
namespace qc {
namespace {
#include "qc_basis_data.inc"

constexpr double kPi = 3.14159265358979323846;

double double_factorial(int n) {  // n!! for n >= -1
  double r = 1.0;
  for (int k = n; k > 1; k -= 2) r *= k;
  return r;
}

// ---- Boys function ---------------------------------------------------------------------

constexpr int kBoysMax = 24;        // highest order tabulated (+6 Taylor terms needed)
constexpr double kBoysStep = 0.1;
constexpr double kBoysTMax = 40.0;
constexpr int kBoysGrid = 401;

struct BoysTable {
  std::vector<double> f;  // f[i * (kBoysMax + 7) + n] = F_n(i * step)
  BoysTable() : f(static_cast<std::size_t>(kBoysGrid) * (kBoysMax + 7)) {
    for (int i = 0; i < kBoysGrid; ++i) {
      const double t = i * kBoysStep;
      for (int n = 0; n < kBoysMax + 7; ++n) {
        // F_n(T) = exp(-T) sum_k (2T)^k / ((2n+1)(2n+3)...(2n+2k+1))
        double term = 1.0 / (2 * n + 1), sum = term;
        for (int k = 1; k < 2000; ++k) {
          term *= 2.0 * t / (2 * n + 2 * k + 1);
          sum += term;
          if (term < 1e-17 * sum) break;
        }
        f[static_cast<std::size_t>(i) * (kBoysMax + 7) + n] = std::exp(-t) * sum;
      }
    }
  }
};

const BoysTable& boys_table() {
  static const BoysTable table;
  return table;
}

// F_0..F_nmax at T.
void boys(int nmax, double t, double* out) {
  if (nmax > kBoysMax) throw std::logic_error("Boys order too high");
  if (t < kBoysTMax) {
    const BoysTable& tb = boys_table();
    const int i = static_cast<int>(std::lround(t / kBoysStep));
    const double d = t - i * kBoysStep;
    const double* row = &tb.f[static_cast<std::size_t>(i) * (kBoysMax + 7)];
    double v = 0.0, fact = 1.0, pw = 1.0;
    for (int k = 0; k < 7; ++k) {  // Taylor: F_n(T_i + d) = sum_k F_{n+k}(T_i) (-d)^k / k!
      v += row[nmax + k] * pw / fact;
      pw *= -d;
      fact *= (k + 1);
    }
    out[nmax] = v;
    const double e = std::exp(-t);
    for (int n = nmax; n > 0; --n) out[n - 1] = (2.0 * t * out[n] + e) / (2 * n - 1);
  } else {
    out[0] = 0.5 * std::sqrt(kPi / t);
    const double e = std::exp(-t);
    for (int n = 0; n < nmax; ++n) out[n + 1] = ((2 * n + 1) * out[n] - e) / (2.0 * t);
  }
}

// ---- Hermite expansion coefficients -----------------------------------------------------

// E[(i * (lb + 1) + j) * (la + lb + 1) + t] for one Cartesian direction.
struct Hermite1D {
  int la = 0, lb = 0;
  std::vector<double> e;
  double at(int i, int j, int t) const {
    if (t < 0 || t > i + j) return 0.0;
    return e[static_cast<std::size_t>((i * (lb + 1) + j) * (la + lb + 1) + t)];
  }
};

Hermite1D hermite(int la, int lb, double a, double b, double xa, double xb) {
  Hermite1D h;
  h.la = la;
  h.lb = lb;
  const int T = la + lb + 1;
  h.e.assign(static_cast<std::size_t>((la + 1) * (lb + 1) * T), 0.0);
  const double p = a + b, mu = a * b / p, x = xa - xb, xp = (a * xa + b * xb) / p;
  const double xpa = xp - xa, xpb = xp - xb, inv2p = 0.5 / p;
  auto E = [&](int i, int j, int t) -> double& { return h.e[static_cast<std::size_t>((i * (lb + 1) + j) * T + t)]; };
  E(0, 0, 0) = std::exp(-mu * x * x);
  for (int i = 0; i <= la; ++i)
    for (int j = 0; j <= lb; ++j) {
      if (i == 0 && j == 0) continue;
      for (int t = 0; t <= i + j; ++t) {
        if (i > 0) {
          E(i, j, t) = (t > 0 ? inv2p * h.at(i - 1, j, t - 1) : 0.0) + xpa * h.at(i - 1, j, t) + (t + 1) * h.at(i - 1, j, t + 1);
        } else {
          E(i, j, t) = (t > 0 ? inv2p * h.at(i, j - 1, t - 1) : 0.0) + xpb * h.at(i, j - 1, t) + (t + 1) * h.at(i, j - 1, t + 1);
        }
      }
    }
  return h;
}

// Hermite Coulomb integrals R_{tuv} (n = 0) for t + u + v <= L. The work array
// holds all orders n; R views its n = 0 layer, so nothing is copied or reallocated
// once the buffer has grown to the largest L.
struct HermiteR {
  int L = 0;
  const double* r = nullptr;  // index (t * (L+1) + u) * (L+1) + v
  double at(int t, int u, int v) const { return r[(t * (L + 1) + u) * (L + 1) + v]; }
};

void hermite_r(int L, double alpha, double x, double y, double z, HermiteR& out, std::vector<double>& work) {
  out.L = L;
  const int D = L + 1;
  const std::size_t layer = static_cast<std::size_t>(D) * D * D;
  if (work.size() < layer * D) work.resize(layer * D);
  double* w = work.data();
  auto W = [&](int n, int t, int u, int v) -> double& {
    return w[static_cast<std::size_t>(n) * layer + static_cast<std::size_t>((t * D + u) * D + v)];
  };
  double f[kBoysMax + 1];
  boys(L, alpha * (x * x + y * y + z * z), f);
  double pw = 1.0;
  for (int n = 0; n <= L; ++n) {
    W(n, 0, 0, 0) = pw * f[n];
    pw *= -2.0 * alpha;
  }
  for (int order = 1; order <= L; ++order)
    for (int n = 0; n + order <= L; ++n)
      for (int t = 0; t <= order; ++t)
        for (int u = 0; t + u <= order; ++u) {
          const int v = order - t - u;
          double val;
          if (t > 0) val = (t > 1 ? (t - 1) * W(n + 1, t - 2, u, v) : 0.0) + x * W(n + 1, t - 1, u, v);
          else if (u > 0) val = (u > 1 ? (u - 1) * W(n + 1, t, u - 2, v) : 0.0) + y * W(n + 1, t, u - 1, v);
          else val = (v > 1 ? (v - 1) * W(n + 1, t, u, v - 2) : 0.0) + z * W(n + 1, t, u, v - 1);
          W(n, t, u, v) = val;
        }
  out.r = w;
}

std::size_t pair_index(std::size_t i, std::size_t j) { return i >= j ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i; }

}  // namespace

std::vector<std::array<int, 3>> cartesian_components(int L) {
  std::vector<std::array<int, 3>> out;
  for (int x = L; x >= 0; --x)
    for (int y = L - x; y >= 0; --y) out.push_back({x, y, L - x - y});
  return out;
}

namespace {
const std::vector<std::array<int, 3>>& components(int L) {  // cached; L <= 4
  static const std::array<std::vector<std::array<int, 3>>, 5> table = {cartesian_components(0), cartesian_components(1),
                                                                       cartesian_components(2), cartesian_components(3),
                                                                       cartesian_components(4)};
  return table[static_cast<std::size_t>(L)];
}
}  // namespace

bool build_basis(const Molecule& mol, const std::string& name, Basis& out, std::string& error) {
  std::string key;
  for (char c : name) key += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  if (key == "6-31G(D)") key = "6-31G*";
  const double* table = nullptr;
  for (const auto& t : k_basis_tables)
    if (key == t.name) table = t.data;
  if (!table) {
    error = "unknown basis set '" + name + "' (STO-3G, 3-21G, 6-31G, 6-31G*)";
    return false;
  }
  out = Basis{};
  for (std::size_t a = 0; a < mol.z.size(); ++a) {
    bool found = false;
    for (const double* p = table; *p != 0.0;) {
      const int z = static_cast<int>(p[0]), L = static_cast<int>(p[1]), np = static_cast<int>(p[2]);
      if (z == mol.z[a]) {
        found = true;
        Shell s;
        s.center = mol.xyz[a];
        s.L = L;
        s.atom = static_cast<int>(a);
        s.first = out.n;
        for (int k = 0; k < np; ++k) {
          const double alpha = p[3 + k], c = p[3 + np + k];
          const double norm = std::pow(2.0 * alpha / kPi, 0.75) * std::pow(4.0 * alpha, 0.5 * L) /
                              std::sqrt(double_factorial(2 * L - 1));
          s.exps.push_back(alpha);
          s.coefs.push_back(c * norm);
        }
        for (const auto& pw : cartesian_components(L)) {
          // self-overlap of this contracted component, from the closed form of
          // integral x^(2l) exp(-p x^2) dx = (2l-1)!! / (2p)^l sqrt(pi/p)
          double ss = 0.0;
          for (int i = 0; i < np; ++i)
            for (int j = 0; j < np; ++j) {
              const double pp = s.exps[static_cast<std::size_t>(i)] + s.exps[static_cast<std::size_t>(j)];
              double v = std::pow(kPi / pp, 1.5);
              for (int d = 0; d < 3; ++d) v *= double_factorial(2 * pw[static_cast<std::size_t>(d)] - 1) / std::pow(2.0 * pp, pw[static_cast<std::size_t>(d)]);
              ss += s.coefs[static_cast<std::size_t>(i)] * s.coefs[static_cast<std::size_t>(j)] * v;
            }
          out.scale.push_back(1.0 / std::sqrt(ss));
          out.atom_of.push_back(static_cast<int>(a));
          out.powers.push_back(pw);
          ++out.n;
        }
        out.max_l = std::max(out.max_l, L);
        out.shells.push_back(std::move(s));
      }
      p += 3 + 2 * np;
    }
    if (!found) {
      error = std::string("basis set ") + name + " has no functions for " + std::string(element_symbol(static_cast<std::uint8_t>(mol.z[a])));
      return false;
    }
  }
  return true;
}

void symmetric_eigen(const Mat& a_in, int n, std::vector<double>& values, Mat& vectors) {
  Mat a = a_in;
  vectors.assign(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) vectors[static_cast<std::size_t>(i) * n + i] = 1.0;
  auto A = [&](int i, int j) -> double& { return a[static_cast<std::size_t>(i) * n + j]; };
  auto V = [&](int i, int j) -> double& { return vectors[static_cast<std::size_t>(i) * n + j]; };
  double scale = 0.0;
  for (double x : a) scale = std::max(scale, std::fabs(x));
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0;
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j) off += A(i, j) * A(i, j);
    if (off <= 1e-30 * (scale * scale + 1e-300) * n * n) break;
    for (int p = 0; p < n; ++p)
      for (int q = p + 1; q < n; ++q) {
        const double apq = A(p, q);
        if (std::fabs(apq) < 1e-300) continue;
        const double theta = (A(q, q) - A(p, p)) / (2.0 * apq);
        const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
        for (int k = 0; k < n; ++k) {  // A <- J^T A J
          const double akp = A(k, p), akq = A(k, q);
          A(k, p) = c * akp - s * akq;
          A(k, q) = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          const double apk = A(p, k), aqk = A(q, k);
          A(p, k) = c * apk - s * aqk;
          A(q, k) = s * apk + c * aqk;
        }
        for (int k = 0; k < n; ++k) {
          const double vkp = V(k, p), vkq = V(k, q);
          V(k, p) = c * vkp - s * vkq;
          V(k, q) = s * vkp + c * vkq;
        }
      }
  }
  std::vector<int> order(static_cast<std::size_t>(n));
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int x, int y) { return A(x, x) < A(y, y); });
  values.resize(static_cast<std::size_t>(n));
  Mat sorted(static_cast<std::size_t>(n) * n);
  for (int k = 0; k < n; ++k) {
    values[static_cast<std::size_t>(k)] = A(order[static_cast<std::size_t>(k)], order[static_cast<std::size_t>(k)]);
    for (int i = 0; i < n; ++i) sorted[static_cast<std::size_t>(i) * n + k] = V(i, order[static_cast<std::size_t>(k)]);
  }
  vectors = std::move(sorted);
}

namespace {

// One-electron integrals: overlap, kinetic, nuclear attraction, dipole (origin 0).
struct OneElectron {
  Mat s, t, v;
  std::array<Mat, 3> dipole;
};

OneElectron one_electron(const Molecule& mol, const Basis& b) {
  const int n = b.n;
  OneElectron o;
  for (Mat* m : {&o.s, &o.t, &o.v, &o.dipole[0], &o.dipole[1], &o.dipole[2]}) m->assign(static_cast<std::size_t>(n) * n, 0.0);
  std::vector<std::vector<double>> nuc_work;
  for (const Shell& A : b.shells)
    for (const Shell& B : b.shells) {
      const auto ca = cartesian_components(A.L), cb = cartesian_components(B.L);
      for (std::size_t ka = 0; ka < A.exps.size(); ++ka)
        for (std::size_t kb = 0; kb < B.exps.size(); ++kb) {
          const double a = A.exps[ka], bb = B.exps[kb], p = a + bb;
          const double dd = A.coefs[ka] * B.coefs[kb];
          std::array<Hermite1D, 3> E;
          std::array<double, 3> P{};
          for (int d = 0; d < 3; ++d) {
            E[static_cast<std::size_t>(d)] = hermite(A.L, B.L + 2, a, bb, A.center[static_cast<std::size_t>(d)], B.center[static_cast<std::size_t>(d)]);
            P[static_cast<std::size_t>(d)] = (a * A.center[static_cast<std::size_t>(d)] + bb * B.center[static_cast<std::size_t>(d)]) / p;
          }
          const double sq = std::sqrt(kPi / p);
          std::vector<HermiteR> nuc(mol.z.size());
          nuc_work.resize(mol.z.size());
          for (std::size_t c = 0; c < mol.z.size(); ++c)
            hermite_r(A.L + B.L, p, P[0] - mol.xyz[c][0], P[1] - mol.xyz[c][1], P[2] - mol.xyz[c][2], nuc[c], nuc_work[c]);
          for (std::size_t ia = 0; ia < ca.size(); ++ia)
            for (std::size_t ib = 0; ib < cb.size(); ++ib) {
              const int mu = A.first + static_cast<int>(ia), nu = B.first + static_cast<int>(ib);
              const auto& pa = ca[ia];
              const auto& pb = cb[ib];
              std::array<double, 3> S{}, T{}, D{};
              for (std::size_t d = 0; d < 3; ++d) {
                const int i = pa[d], j = pb[d];
                const Hermite1D& e = E[d];
                S[d] = e.at(i, j, 0) * sq;
                const double sm2 = j >= 2 ? e.at(i, j - 2, 0) * sq : 0.0;
                T[d] = -0.5 * (j * (j - 1) * sm2 - 2.0 * bb * (2 * j + 1) * S[d] + 4.0 * bb * bb * e.at(i, j + 2, 0) * sq);
                D[d] = (e.at(i, j, 1) + P[d] * e.at(i, j, 0)) * sq;
              }
              const std::size_t idx = static_cast<std::size_t>(mu) * n + nu;
              o.s[idx] += dd * S[0] * S[1] * S[2];
              o.t[idx] += dd * (T[0] * S[1] * S[2] + S[0] * T[1] * S[2] + S[0] * S[1] * T[2]);
              o.dipole[0][idx] += dd * D[0] * S[1] * S[2];
              o.dipole[1][idx] += dd * S[0] * D[1] * S[2];
              o.dipole[2][idx] += dd * S[0] * S[1] * D[2];
              double vsum = 0.0;
              for (std::size_t c = 0; c < mol.z.size(); ++c) {
                double acc = 0.0;
                for (int t = 0; t <= pa[0] + pb[0]; ++t)
                  for (int u = 0; u <= pa[1] + pb[1]; ++u)
                    for (int w = 0; w <= pa[2] + pb[2]; ++w)
                      acc += E[0].at(pa[0], pb[0], t) * E[1].at(pa[1], pb[1], u) * E[2].at(pa[2], pb[2], w) * nuc[c].at(t, u, w);
                vsum -= mol.z[c] * acc;
              }
              o.v[idx] += dd * 2.0 * kPi / p * vsum;
            }
        }
    }
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      const double f = b.scale[static_cast<std::size_t>(i)] * b.scale[static_cast<std::size_t>(j)];
      const std::size_t idx = static_cast<std::size_t>(i) * n + j;
      o.s[idx] *= f;
      o.t[idx] *= f;
      o.v[idx] *= f;
      for (auto& m : o.dipole) m[idx] *= f;
    }
  return o;
}

// Primitive-pair data for the two-electron integrals.
struct PrimPair {
  double p, coef;
  std::array<double, 3> P;
  std::array<Hermite1D, 3> E;
};

struct ShellPair {
  int a, b;  // shell indices, a >= b
  std::vector<PrimPair> prims;
};

ShellPair make_pair(const Basis& basis, int ia, int ib) {
  const Shell& A = basis.shells[static_cast<std::size_t>(ia)];
  const Shell& B = basis.shells[static_cast<std::size_t>(ib)];
  ShellPair sp{ia, ib, {}};
  for (std::size_t ka = 0; ka < A.exps.size(); ++ka)
    for (std::size_t kb = 0; kb < B.exps.size(); ++kb) {
      PrimPair pp;
      const double a = A.exps[ka], b = B.exps[kb];
      pp.p = a + b;
      pp.coef = A.coefs[ka] * B.coefs[kb];
      double ab2 = 0.0;
      for (std::size_t d = 0; d < 3; ++d) ab2 += (A.center[d] - B.center[d]) * (A.center[d] - B.center[d]);
      // Gaussian product prefactor: a pair below 1e-15 (with its pi/p factor) contributes nothing measurable
      if (std::fabs(pp.coef) * std::exp(-a * b / pp.p * ab2) * std::pow(kPi / pp.p, 1.5) < 1e-15) continue;
      for (std::size_t d = 0; d < 3; ++d) {
        pp.E[d] = hermite(A.L, B.L, a, b, A.center[d], B.center[d]);
        pp.P[d] = (a * A.center[d] + b * B.center[d]) / pp.p;
      }
      sp.prims.push_back(std::move(pp));
    }
  return sp;
}

// (ab|cd) for every component of the four shells, written into `out` (component-major).
void shell_quartet(const Basis& basis, const ShellPair& ab, const ShellPair& cd, std::vector<double>& out,
                   std::vector<double>& work, std::vector<double>& g) {
  const Shell& A = basis.shells[static_cast<std::size_t>(ab.a)];
  const Shell& B = basis.shells[static_cast<std::size_t>(ab.b)];
  const Shell& C = basis.shells[static_cast<std::size_t>(cd.a)];
  const Shell& D = basis.shells[static_cast<std::size_t>(cd.b)];
  const auto& ca = components(A.L);
  const auto& cb = components(B.L);
  const auto& cc = components(C.L);
  const auto& cdc = components(D.L);
  const int Lab = A.L + B.L, Lcd = C.L + D.L, L = Lab + Lcd;
  const std::size_t nab = ca.size() * cb.size(), ncd = cc.size() * cdc.size();
  out.assign(nab * ncd, 0.0);
  const int Dab = Lab + 1;
  HermiteR R;
  for (const PrimPair& p1 : ab.prims)
    for (const PrimPair& p2 : cd.prims) {
      const double p = p1.p, q = p2.p, alpha = p * q / (p + q);
      hermite_r(L, alpha, p1.P[0] - p2.P[0], p1.P[1] - p2.P[1], p1.P[2] - p2.P[2], R, work);
      const double pref = p1.coef * p2.coef * 2.0 * std::pow(kPi, 2.5) / (p * q * std::sqrt(p + q));
      for (std::size_t ic = 0; ic < cc.size(); ++ic)
        for (std::size_t id = 0; id < cdc.size(); ++id) {
          const auto& pc = cc[ic];
          const auto& pd = cdc[id];
          // g_{tuv} = sum (-1)^(tau+nu+phi) E^cd_tau E^cd_nu E^cd_phi R_{t+tau, u+nu, v+phi}
          g.assign(static_cast<std::size_t>(Dab) * Dab * Dab, 0.0);
          for (int tau = 0; tau <= pc[0] + pd[0]; ++tau) {
            const double ex = p2.E[0].at(pc[0], pd[0], tau);
            if (ex == 0.0) continue;
            for (int nu = 0; nu <= pc[1] + pd[1]; ++nu) {
              const double ey = p2.E[1].at(pc[1], pd[1], nu);
              if (ey == 0.0) continue;
              for (int phi = 0; phi <= pc[2] + pd[2]; ++phi) {
                const double ez = p2.E[2].at(pc[2], pd[2], phi);
                const double f = ((tau + nu + phi) & 1 ? -1.0 : 1.0) * ex * ey * ez;
                if (f == 0.0) continue;
                for (int t = 0; t <= Lab; ++t)
                  for (int u = 0; t + u <= Lab; ++u)
                    for (int v = 0; t + u + v <= Lab; ++v)
                      g[static_cast<std::size_t>((t * Dab + u) * Dab + v)] += f * R.at(t + tau, u + nu, v + phi);
              }
            }
          }
          const std::size_t col = ic * cdc.size() + id;
          for (std::size_t ia = 0; ia < ca.size(); ++ia)
            for (std::size_t ib = 0; ib < cb.size(); ++ib) {
              const auto& pa = ca[ia];
              const auto& pb = cb[ib];
              double acc = 0.0;
              for (int t = 0; t <= pa[0] + pb[0]; ++t) {
                const double ex = p1.E[0].at(pa[0], pb[0], t);
                for (int u = 0; u <= pa[1] + pb[1]; ++u) {
                  const double ey = p1.E[1].at(pa[1], pb[1], u);
                  for (int v = 0; v <= pa[2] + pb[2]; ++v)
                    acc += ex * ey * p1.E[2].at(pa[2], pb[2], v) * g[static_cast<std::size_t>((t * Dab + u) * Dab + v)];
                }
              }
              out[(ia * cb.size() + ib) * ncd + col] += pref * acc;
            }
        }
    }
}

// All unique two-electron integrals, packed (ij|kl) with i>=j, k>=l, ij>=kl.
struct Eri {
  int n = 0;
  std::vector<double> v;
  double at(int i, int j, int k, int l) const { return v[pair_index(pair_index(static_cast<std::size_t>(i), static_cast<std::size_t>(j)), pair_index(static_cast<std::size_t>(k), static_cast<std::size_t>(l)))]; }
};

bool two_electron(const Basis& basis, double limit_mb, Eri& eri, std::string& error) {
  const std::size_t n = static_cast<std::size_t>(basis.n);
  const std::size_t m = n * (n + 1) / 2;
  const double mb = static_cast<double>(m) * static_cast<double>(m + 1) / 2.0 * 8.0 / 1048576.0;
  if (mb > limit_mb) {
    error = "two-electron integrals need " + std::to_string(static_cast<long long>(mb)) + " MB (" + std::to_string(n) +
            " basis functions); the limit is " + std::to_string(static_cast<long long>(limit_mb)) + " MB: use a smaller basis or molecule";
    return false;
  }
  eri.n = basis.n;
  eri.v.assign(m * (m + 1) / 2, 0.0);
  const int ns = static_cast<int>(basis.shells.size());
  std::vector<ShellPair> pairs;
  for (int a = 0; a < ns; ++a)
    for (int b = 0; b <= a; ++b) pairs.push_back(make_pair(basis, a, b));
  // Schwarz bounds per shell pair. Each shell pair (and each row x below) writes only its own outputs, so
  // the loops run in parallel with results independent of the thread count.
  std::vector<double> q(pairs.size(), 0.0);
  parallel_for(pairs.size(), [&](std::size_t k) {
    std::vector<double> buf, work, g;
    shell_quartet(basis, pairs[k], pairs[k], buf, work, g);
    const Shell& A = basis.shells[static_cast<std::size_t>(pairs[k].a)];
    const Shell& B = basis.shells[static_cast<std::size_t>(pairs[k].b)];
    const std::size_t na = cartesian_components(A.L).size(), nb = cartesian_components(B.L).size(), nab = na * nb;
    double mx = 0.0;
    for (std::size_t i = 0; i < na; ++i)
      for (std::size_t j = 0; j < nb; ++j) {
        const double s = basis.scale[static_cast<std::size_t>(A.first) + i] * basis.scale[static_cast<std::size_t>(B.first) + j];
        mx = std::max(mx, std::fabs(buf[(i * nb + j) * nab + i * nb + j]) * s * s);
      }
    q[k] = std::sqrt(mx);
  });
  parallel_for(pairs.size(), [&](std::size_t xr) {
    const std::size_t x = pairs.size() - 1 - xr;  // longest rows first
    std::vector<double> buf, work, g;
    for (std::size_t y = 0; y <= x; ++y) {
      if (q[x] * q[y] < 1e-12) continue;
      shell_quartet(basis, pairs[x], pairs[y], buf, work, g);
      const Shell& A = basis.shells[static_cast<std::size_t>(pairs[x].a)];
      const Shell& B = basis.shells[static_cast<std::size_t>(pairs[x].b)];
      const Shell& C = basis.shells[static_cast<std::size_t>(pairs[y].a)];
      const Shell& D = basis.shells[static_cast<std::size_t>(pairs[y].b)];
      const std::size_t na = cartesian_components(A.L).size(), nb = cartesian_components(B.L).size();
      const std::size_t nc = cartesian_components(C.L).size(), nd = cartesian_components(D.L).size();
      for (std::size_t i = 0; i < na; ++i)
        for (std::size_t j = 0; j < nb; ++j)
          for (std::size_t k = 0; k < nc; ++k)
            for (std::size_t l = 0; l < nd; ++l) {
              const std::size_t fi = static_cast<std::size_t>(A.first) + i, fj = static_cast<std::size_t>(B.first) + j;
              const std::size_t fk = static_cast<std::size_t>(C.first) + k, fl = static_cast<std::size_t>(D.first) + l;
              const double s = basis.scale[fi] * basis.scale[fj] * basis.scale[fk] * basis.scale[fl];
              eri.v[pair_index(pair_index(fi, fj), pair_index(fk, fl))] = buf[(i * nb + j) * (nc * nd) + k * nd + l] * s;
            }
    }
  });
  return true;
}

// J (from pt) and K (from each of pa, pb) using every unique integral once. Each
// unique (ij|kl) is scaled by its degeneracy (halved for i=j, k=l and ij=kl) and
// applied in all eight permutations, which needs no duplicate tests.
void build_jk(const Eri& eri, const Mat& pt, const Mat* pa, const Mat* pb, Mat& j, Mat* ka, Mat* kb) {
  const std::size_t n = static_cast<std::size_t>(eri.n);
  const double* v = eri.v.data();
  // pairs in packed order: ij = i(i+1)/2 + j
  std::vector<std::pair<std::size_t, std::size_t>> pairs;
  pairs.reserve(n * (n + 1) / 2);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t jj = 0; jj <= i; ++jj) pairs.push_back({i, jj});
  const std::size_t np = pairs.size();
  const double* P = pt.data();
  // Rows ij are split into a fixed number of parts with about the same number of integrals (row ij has
  // ij+1); each part accumulates its own J and K, and the parts are added in order. The number of parts
  // depends on the basis size only, so the result is the same for any thread count.
  const std::size_t nmat = 1 + (ka ? 1 : 0) + (kb ? 1 : 0);
  const std::size_t cap = std::max<std::size_t>(1, static_cast<std::size_t>(3e7 / static_cast<double>(nmat * n * n + 1)));
  const std::size_t parts = std::max<std::size_t>(1, std::min<std::size_t>({np, 32, cap}));
  std::vector<std::size_t> start(parts + 1, np);
  {
    const double total = 0.5 * static_cast<double>(np) * static_cast<double>(np + 1);
    std::size_t row = 0;
    for (std::size_t k = 0; k < parts; ++k) {
      const double target = total * static_cast<double>(k) / static_cast<double>(parts);
      while (row < np && 0.5 * static_cast<double>(row) * static_cast<double>(row + 1) < target) ++row;
      start[k] = row;
    }
  }
  std::vector<Mat> pj(parts), pka(parts), pkb(parts);
  parallel_for(parts, [&](std::size_t part) {
    Mat& Jm = pj[part];
    Jm.assign(n * n, 0.0);
    if (ka) pka[part].assign(n * n, 0.0);
    if (kb) pkb[part].assign(n * n, 0.0);
    double* J = Jm.data();
    auto exchange = [&](const Mat& pm, Mat& km, std::size_t i, std::size_t jj, std::size_t k, std::size_t l, double x) {
      const double* Q = pm.data();
      double* K = km.data();
      K[i * n + k] += Q[jj * n + l] * x;
      K[jj * n + k] += Q[i * n + l] * x;
      K[i * n + l] += Q[jj * n + k] * x;
      K[jj * n + l] += Q[i * n + k] * x;
      K[k * n + i] += Q[l * n + jj] * x;
      K[l * n + i] += Q[k * n + jj] * x;
      K[k * n + jj] += Q[l * n + i] * x;
      K[l * n + jj] += Q[k * n + i] * x;
    };
    for (std::size_t ij = start[part]; ij < start[part + 1]; ++ij) {
      const auto [i, jj] = pairs[ij];
      std::size_t idx = ij * (ij + 1) / 2;
      for (std::size_t kl = 0; kl <= ij; ++kl, ++idx) {
        double x = v[idx];
        if (x == 0.0) continue;
        const auto [k, l] = pairs[kl];
        if (i == jj) x *= 0.5;
        if (k == l) x *= 0.5;
        if (ij == kl) x *= 0.5;
        const double pij = P[i * n + jj] + P[jj * n + i], pkl = P[k * n + l] + P[l * n + k];
        J[i * n + jj] += pkl * x;
        J[jj * n + i] += pkl * x;
        J[k * n + l] += pij * x;
        J[l * n + k] += pij * x;
        if (ka) exchange(*pa, pka[part], i, jj, k, l, x);
        if (kb) exchange(*pb, pkb[part], i, jj, k, l, x);
      }
    }
  });
  j.assign(n * n, 0.0);
  if (ka) ka->assign(n * n, 0.0);
  if (kb) kb->assign(n * n, 0.0);
  for (std::size_t part = 0; part < parts; ++part)
    for (std::size_t t = 0; t < n * n; ++t) {
      j[t] += pj[part][t];
      if (ka) (*ka)[t] += pka[part][t];
      if (kb) (*kb)[t] += pkb[part][t];
    }
}

Mat mul(const Mat& a, const Mat& b, int n, bool ta = false, bool tb = false) {
  Mat c(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < n; ++k) {
      const double aik = ta ? a[static_cast<std::size_t>(k) * n + i] : a[static_cast<std::size_t>(i) * n + k];
      if (aik == 0.0) continue;
      for (int j = 0; j < n; ++j) c[static_cast<std::size_t>(i) * n + j] += aik * (tb ? b[static_cast<std::size_t>(j) * n + k] : b[static_cast<std::size_t>(k) * n + j]);
    }
  return c;
}

double dot(const Mat& a, const Mat& b) {
  double s = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// Canonical orthogonaliser X = U s^-1/2 over overlap eigenvalues above 1e-7.
void orthogonaliser(const Mat& s, int n, Mat& X, int& m) {
  std::vector<double> val;
  Mat vec;
  symmetric_eigen(s, n, val, vec);
  std::vector<int> keep;
  for (int k = 0; k < n; ++k)
    if (val[static_cast<std::size_t>(k)] > 1e-7) keep.push_back(k);
  m = static_cast<int>(keep.size());
  X.assign(static_cast<std::size_t>(n) * m, 0.0);
  for (int i = 0; i < n; ++i)
    for (int c = 0; c < m; ++c) {
      const int k = keep[static_cast<std::size_t>(c)];
      X[static_cast<std::size_t>(i) * m + c] = vec[static_cast<std::size_t>(i) * n + k] / std::sqrt(val[static_cast<std::size_t>(k)]);
    }
}

// Solves F C = S C e in the orthogonal basis X; C is n x m.
void solve_fock(const Mat& X, int n, int m, const Mat& f, std::vector<double>& eps, Mat& c) {
  Mat tmp(static_cast<std::size_t>(n) * m, 0.0), fx(static_cast<std::size_t>(m) * m, 0.0);
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < n; ++k) {
      const double fik = f[static_cast<std::size_t>(i) * n + k];
      if (fik == 0.0) continue;
      for (int b = 0; b < m; ++b) tmp[static_cast<std::size_t>(i) * m + b] += fik * X[static_cast<std::size_t>(k) * m + b];
    }
  for (int a = 0; a < m; ++a)
    for (int i = 0; i < n; ++i) {
      const double xia = X[static_cast<std::size_t>(i) * m + a];
      for (int b = 0; b < m; ++b) fx[static_cast<std::size_t>(a) * m + b] += xia * tmp[static_cast<std::size_t>(i) * m + b];
    }
  Mat v;
  symmetric_eigen(fx, m, eps, v);
  c.assign(static_cast<std::size_t>(n) * m, 0.0);
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < m; ++k) {
      const double xik = X[static_cast<std::size_t>(i) * m + k];
      for (int b = 0; b < m; ++b) c[static_cast<std::size_t>(i) * m + b] += xik * v[static_cast<std::size_t>(k) * m + b];
    }
}

// Spherically averaged neutral-atom density: restricted SCF in which electrons
// beyond the closed shells are spread evenly over the degenerate frontier
// orbitals (e.g. O 2p^4 as 4/3 per p orbital), with simple damping.
Mat atom_density(int z, const std::string& basis_name) {
  Molecule m1;
  m1.z = {z};
  m1.xyz = {{0.0, 0.0, 0.0}};
  Basis b;
  std::string err;
  if (!build_basis(m1, basis_name, b, err)) return {};
  const int n = b.n;
  const OneElectron o = one_electron(m1, b);
  Mat h(o.t.size());
  for (std::size_t i = 0; i < h.size(); ++i) h[i] = o.t[i] + o.v[i];
  Eri e;
  if (!two_electron(b, 1e9, e, err)) return {};
  Mat X;
  int m = 0;
  orthogonaliser(o.s, n, X, m);
  Mat p(static_cast<std::size_t>(n) * n, 0.0);
  for (int iter = 0; iter < 300; ++iter) {
    Mat half(p), J, K;
    for (double& x : half) x *= 0.5;
    build_jk(e, p, &half, nullptr, J, &K, nullptr);
    Mat f(h);
    for (std::size_t i = 0; i < f.size(); ++i) f[i] += J[i] - K[i];
    std::vector<double> eps;
    Mat c;
    solve_fock(X, n, m, f, eps, c);
    std::vector<double> occ(static_cast<std::size_t>(m), 0.0);
    double left = z;
    for (int i = 0; i < m && left > 1e-12;) {
      int j = i;
      while (j < m && eps[static_cast<std::size_t>(j)] - eps[static_cast<std::size_t>(i)] < 1e-3) ++j;
      const double put = std::min(left, 2.0 * (j - i));
      for (int k = i; k < j; ++k) occ[static_cast<std::size_t>(k)] = put / (j - i);
      left -= put;
      i = j;
    }
    Mat pn(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        double sum = 0.0;
        for (int k = 0; k < m; ++k) sum += occ[static_cast<std::size_t>(k)] * c[static_cast<std::size_t>(i) * m + k] * c[static_cast<std::size_t>(j) * m + k];
        pn[static_cast<std::size_t>(i) * n + j] = sum;
      }
    double change = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) change = std::max(change, std::fabs(pn[i] - p[i]));
    const double mix = iter == 0 ? 1.0 : 0.5;
    for (std::size_t i = 0; i < p.size(); ++i) p[i] = (1.0 - mix) * p[i] + mix * pn[i];
    if (change < 1e-8) break;
  }
  return p;
}

// Superposition of atomic densities: atomic blocks on the diagonal of the total density.
Mat sad_guess(const Molecule& mol, const Basis& basis, const std::string& basis_name) {
  const int n = basis.n;
  Mat p(static_cast<std::size_t>(n) * n, 0.0);
  std::map<int, Mat> cache;
  int offset = 0;
  for (std::size_t a = 0; a < mol.z.size(); ++a) {
    auto it = cache.find(mol.z[a]);
    if (it == cache.end()) it = cache.emplace(mol.z[a], atom_density(mol.z[a], basis_name)).first;
    const Mat& pa = it->second;
    int na = 0;
    while (offset + na < n && basis.atom_of[static_cast<std::size_t>(offset + na)] == static_cast<int>(a)) ++na;
    if (static_cast<std::size_t>(na) * na == pa.size())
      for (int i = 0; i < na; ++i)
        for (int j = 0; j < na; ++j) p[static_cast<std::size_t>(offset + i) * n + offset + j] = pa[static_cast<std::size_t>(i) * na + j];
    offset += na;
  }
  return p;
}

}  // namespace
}  // namespace qc

// ---- public ------------------------------------------------------------------------------

const char* to_string(QcMethod m) noexcept {
  switch (m) {
    case QcMethod::HF: return "HF";
    case QcMethod::SVWN5: return "SVWN5";
    case QcMethod::PBE: return "PBE";
    case QcMethod::B3LYP: return "B3LYP";
  }
  return "HF";
}

bool parse_qc_method(const std::string& name, QcMethod& out) {
  std::string k;
  for (char c : name) k += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  if (k == "hf" || k == "rhf" || k == "uhf" || k == "scf") out = QcMethod::HF;
  else if (k == "svwn" || k == "svwn5" || k == "lda" || k == "lsda") out = QcMethod::SVWN5;
  else if (k == "pbe") out = QcMethod::PBE;
  else if (k == "b3lyp") out = QcMethod::B3LYP;
  else return false;
  return true;
}

QcResult run_qc(const ChemicalSystemState& st, std::uint32_t sp, Diagnostics& diag, const QcOptions& opt) {
  using namespace qc;
  QcResult res;
  res.method = to_string(opt.method);
  const Species& S = st.species[sp];
  Molecule mol;
  int electrons = 0;
  for (std::uint32_t g = S.atom_begin; g < S.atom_begin + S.atom_count; ++g) {
    if (st.atoms.implicit_h[g] > 0) {
      diag.error(DiagCode::ImplicitHydrogens, "quantum chemistry needs explicit hydrogens and 3D coordinates (embed first)").on_species(sp);
      return res;
    }
    mol.z.push_back(st.atoms.element[g]);
    mol.xyz.push_back({st.atoms.x[g] * kAngstromToBohr, st.atoms.y[g] * kAngstromToBohr, st.atoms.z[g] * kAngstromToBohr});
    electrons += st.atoms.element[g];
  }
  mol.charge = S.net_charge;
  electrons -= mol.charge;
  res.charge = mol.charge;
  Basis basis;
  std::string err;
  if (!build_basis(mol, opt.basis, basis, err)) {
    diag.error(DiagCode::InvalidInput, err).on_species(sp);
    return res;
  }
  res.basis = opt.basis;
  int mult = opt.multiplicity > 0 ? opt.multiplicity : (electrons % 2 ? 2 : 1);
  if ((electrons + mult - 1) % 2 != 0 || mult - 1 > electrons) {
    diag.error(DiagCode::InvalidInput, "multiplicity " + std::to_string(mult) + " is impossible with " + std::to_string(electrons) + " electrons").on_species(sp);
    return res;
  }
  const int na = (electrons + mult - 1) / 2, nb = electrons - na;
  const int n = basis.n;
  res.n_basis = n;
  res.n_alpha = na;
  res.n_beta = nb;
  res.multiplicity = mult;
  const bool uhf = mult > 1;
  res.unrestricted = uhf;

  for (std::size_t a = 0; a < mol.z.size(); ++a)
    for (std::size_t b = 0; b < a; ++b) {
      const double dx = mol.xyz[a][0] - mol.xyz[b][0], dy = mol.xyz[a][1] - mol.xyz[b][1], dz = mol.xyz[a][2] - mol.xyz[b][2];
      const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (r < 1e-3) {
        diag.error(DiagCode::InvalidInput, "two atoms are at the same position").on_species(sp);
        return res;
      }
      res.nuclear_repulsion += mol.z[a] * mol.z[b] / r;
    }

  const OneElectron one = one_electron(mol, basis);
  Mat h(one.t.size());
  for (std::size_t i = 0; i < h.size(); ++i) h[i] = one.t[i] + one.v[i];
  Eri eri;
  if (!two_electron(basis, opt.memory_limit_mb, eri, err)) {
    diag.error(DiagCode::InvalidInput, err).on_species(sp);
    return res;
  }
  // Canonical orthogonalisation X = U s^-1/2 over the kept eigenvalues.
  std::vector<double> sval;
  Mat svec;
  symmetric_eigen(one.s, n, sval, svec);
  std::vector<int> keep;
  for (int k = 0; k < n; ++k)
    if (sval[static_cast<std::size_t>(k)] > 1e-7) keep.push_back(k);
  const int m = static_cast<int>(keep.size());
  if (m < na) {
    diag.error(DiagCode::InvalidInput, "basis is linearly dependent").on_species(sp);
    return res;
  }
  if (m < n) diag.warn(DiagCode::InvalidInput, std::to_string(n - m) + " near-linearly-dependent basis combinations removed").on_species(sp);
  Mat X(static_cast<std::size_t>(n) * m, 0.0);
  for (int i = 0; i < n; ++i)
    for (int c = 0; c < m; ++c) {
      const int k = keep[static_cast<std::size_t>(c)];
      X[static_cast<std::size_t>(i) * m + c] = svec[static_cast<std::size_t>(i) * n + k] / std::sqrt(sval[static_cast<std::size_t>(k)]);
    }
  // F (n x n) -> eigenvalues and AO coefficients C (n x m)
  auto diagonalise = [&](const Mat& f, std::vector<double>& eps, Mat& c) {
    Mat fx(static_cast<std::size_t>(m) * m, 0.0), tmp(static_cast<std::size_t>(n) * m, 0.0);
    for (int i = 0; i < n; ++i)
      for (int k = 0; k < n; ++k) {
        const double fik = f[static_cast<std::size_t>(i) * n + k];
        if (fik == 0.0) continue;
        for (int b = 0; b < m; ++b) tmp[static_cast<std::size_t>(i) * m + b] += fik * X[static_cast<std::size_t>(k) * m + b];
      }
    for (int a = 0; a < m; ++a)
      for (int i = 0; i < n; ++i) {
        const double xia = X[static_cast<std::size_t>(i) * m + a];
        for (int b = 0; b < m; ++b) fx[static_cast<std::size_t>(a) * m + b] += xia * tmp[static_cast<std::size_t>(i) * m + b];
      }
    Mat v;
    symmetric_eigen(fx, m, eps, v);
    c.assign(static_cast<std::size_t>(n) * m, 0.0);
    for (int i = 0; i < n; ++i)
      for (int k = 0; k < m; ++k)
        for (int b = 0; b < m; ++b) c[static_cast<std::size_t>(i) * m + b] += X[static_cast<std::size_t>(i) * m + k] * v[static_cast<std::size_t>(k) * m + b];
  };
  auto density = [&](const Mat& c, int occ) {
    Mat p(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        double s = 0.0;
        for (int k = 0; k < occ; ++k) s += c[static_cast<std::size_t>(i) * m + k] * c[static_cast<std::size_t>(j) * m + k];
        p[static_cast<std::size_t>(i) * n + j] = s;
      }
    return p;
  };

  const bool dft = opt.method != QcMethod::HF;
  const Functional fn = opt.method == QcMethod::SVWN5 ? Functional::SVWN5 : opt.method == QcMethod::PBE ? Functional::PBE : Functional::B3LYP;
  const double hfx = dft ? exact_exchange_fraction(fn) : 1.0;
  std::unique_ptr<XcIntegrator> xc;
  if (dft) xc = std::make_unique<XcIntegrator>(mol, basis, fn, opt.grid_level);

  // Guess: orbitals of the Hartree-Fock Fock matrix built from a superposition of
  // atomic densities; occupying na alpha and nb beta of them starts open shells.
  std::vector<double> ea, eb;
  Mat ca, cb;
  {
    const Mat p0 = sad_guess(mol, basis, opt.basis);
    Mat half(p0), J, K;
    for (double& x : half) x *= 0.5;
    build_jk(eri, p0, &half, nullptr, J, &K, nullptr);
    Mat f(h);
    for (std::size_t i = 0; i < f.size(); ++i) f[i] += J[i] - K[i];
    diagonalise(f, ea, ca);
  }
  cb = ca;
  eb = ea;
  Mat pa = density(ca, na), pb = density(cb, nb);
  std::vector<Mat> diis_f, diis_e;  // stacked alpha|beta
  double e_old = 0.0;
  bool converged = false;
  int it = 0;
  for (it = 1; it <= opt.max_iterations; ++it) {
    Mat pt(pa.size());
    for (std::size_t i = 0; i < pt.size(); ++i) pt[i] = pa[i] + pb[i];
    Mat J, Ka, Kb;
    if (hfx != 0.0) build_jk(eri, pt, &pa, &pb, J, &Ka, &Kb);
    else build_jk(eri, pt, nullptr, nullptr, J, nullptr, nullptr);
    Mat fa(h), fb(h);
    for (std::size_t i = 0; i < fa.size(); ++i) {
      fa[i] += J[i] - (hfx != 0.0 ? hfx * Ka[i] : 0.0);
      fb[i] += J[i] - (hfx != 0.0 ? hfx * Kb[i] : 0.0);
    }
    double exc = 0.0;
    if (dft) {
      const XcResult x = xc->evaluate(pa, pb, !uhf);
      exc = x.energy;
      res.electrons_on_grid = x.electrons;
      for (std::size_t i = 0; i < fa.size(); ++i) {
        fa[i] += x.va[i];
        fb[i] += x.vb[i];
      }
    }
    const double e1 = dot(pt, h), ej = 0.5 * dot(pt, J);
    const double ek = hfx != 0.0 ? -0.5 * hfx * (dot(pa, Ka) + dot(pb, Kb)) : 0.0;
    const double energy = e1 + ej + ek + exc + res.nuclear_repulsion;
    res.one_electron = e1;
    res.two_electron = ej + ek;
    res.xc = exc;
    // DIIS error X^T (F P S - S P F) X per spin
    auto err_of = [&](const Mat& f, const Mat& p) {
      const Mat fps = mul(mul(f, p, n), one.s, n), spf = mul(mul(one.s, p, n), f, n);
      Mat e(static_cast<std::size_t>(m) * m, 0.0), tmp(static_cast<std::size_t>(n) * m, 0.0);
      for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k) {
          const double d = fps[static_cast<std::size_t>(i) * n + k] - spf[static_cast<std::size_t>(i) * n + k];
          if (d == 0.0) continue;
          for (int b = 0; b < m; ++b) tmp[static_cast<std::size_t>(i) * m + b] += d * X[static_cast<std::size_t>(k) * m + b];
        }
      for (int a = 0; a < m; ++a)
        for (int i = 0; i < n; ++i)
          for (int b = 0; b < m; ++b) e[static_cast<std::size_t>(a) * m + b] += X[static_cast<std::size_t>(i) * m + a] * tmp[static_cast<std::size_t>(i) * m + b];
      return e;
    };
    Mat e = err_of(fa, pa);
    if (uhf) {
      const Mat e2 = err_of(fb, pb);
      e.insert(e.end(), e2.begin(), e2.end());
    }
    double emax = 0.0;
    for (double x : e) emax = std::max(emax, std::fabs(x));
    if (it > 1 && std::fabs(energy - e_old) < opt.energy_tolerance && emax < opt.density_tolerance) {
      res.energy = energy;
      converged = true;
      break;
    }
    e_old = energy;
    res.energy = energy;
    Mat stacked = fa;
    if (uhf) stacked.insert(stacked.end(), fb.begin(), fb.end());
    diis_f.push_back(stacked);
    diis_e.push_back(e);
    if (diis_f.size() > 8) {
      diis_f.erase(diis_f.begin());
      diis_e.erase(diis_e.begin());
    }
    Mat fuse = stacked;
    const int k = static_cast<int>(diis_f.size());
    if (k >= 2) {
      // B c = rhs with the Lagrange row; solve by Gaussian elimination with pivoting
      const int dim = k + 1;
      std::vector<double> B(static_cast<std::size_t>(dim) * dim, 0.0), rhs(static_cast<std::size_t>(dim), 0.0);
      for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j) B[static_cast<std::size_t>(i) * dim + j] = dot(diis_e[static_cast<std::size_t>(i)], diis_e[static_cast<std::size_t>(j)]);
      for (int i = 0; i < k; ++i) B[static_cast<std::size_t>(i) * dim + k] = B[static_cast<std::size_t>(k) * dim + i] = -1.0;
      rhs[static_cast<std::size_t>(k)] = -1.0;
      bool ok = true;
      for (int c = 0; c < dim && ok; ++c) {
        int piv = c;
        for (int r = c + 1; r < dim; ++r)
          if (std::fabs(B[static_cast<std::size_t>(r) * dim + c]) > std::fabs(B[static_cast<std::size_t>(piv) * dim + c])) piv = r;
        if (std::fabs(B[static_cast<std::size_t>(piv) * dim + c]) < 1e-30) { ok = false; break; }
        for (int q = 0; q < dim; ++q) std::swap(B[static_cast<std::size_t>(c) * dim + q], B[static_cast<std::size_t>(piv) * dim + q]);
        std::swap(rhs[static_cast<std::size_t>(c)], rhs[static_cast<std::size_t>(piv)]);
        for (int r = 0; r < dim; ++r) {
          if (r == c) continue;
          const double f = B[static_cast<std::size_t>(r) * dim + c] / B[static_cast<std::size_t>(c) * dim + c];
          for (int q = c; q < dim; ++q) B[static_cast<std::size_t>(r) * dim + q] -= f * B[static_cast<std::size_t>(c) * dim + q];
          rhs[static_cast<std::size_t>(r)] -= f * rhs[static_cast<std::size_t>(c)];
        }
      }
      if (ok) {
        fuse.assign(stacked.size(), 0.0);
        for (int i = 0; i < k; ++i) {
          const double ci = rhs[static_cast<std::size_t>(i)] / B[static_cast<std::size_t>(i) * dim + i];
          for (std::size_t q = 0; q < fuse.size(); ++q) fuse[q] += ci * diis_f[static_cast<std::size_t>(i)][q];
        }
      }
    }
    Mat use_a(fuse.begin(), fuse.begin() + static_cast<std::ptrdiff_t>(static_cast<std::size_t>(n) * n));
    Mat use_b;
    if (uhf) use_b.assign(fuse.begin() + static_cast<std::ptrdiff_t>(static_cast<std::size_t>(n) * n), fuse.end());
    // Level shift (Saunders & Hillier, IJQC 7 (1973) 699) once DIIS has had 25
    // iterations: F + b (S - S P S) raises the virtual orbitals by b and damps
    // occupied-virtual oscillation; the converged solution is unchanged.
    if (it > 25) {
      const double shift = 0.5;
      auto apply = [&](Mat& f, const Mat& p) {
        const Mat sps = mul(mul(one.s, p, n), one.s, n);
        for (std::size_t q = 0; q < f.size(); ++q) f[q] += shift * (one.s[q] - sps[q]);
      };
      apply(use_a, uhf ? pa : Mat(pa));
      if (uhf) apply(use_b, pb);
    }
    diagonalise(use_a, ea, ca);
    if (uhf) {
      diagonalise(use_b, eb, cb);
    } else {
      cb = ca;
      eb = ea;
    }
    pa = density(ca, na);
    pb = density(cb, nb);
  }
  res.iterations = std::min(it, opt.max_iterations);
  res.converged = converged;
  res.orbital_energies_alpha = ea;
  if (uhf) res.orbital_energies_beta = eb;
  // Mulliken charges and dipole from the final density
  Mat pt(pa.size());
  for (std::size_t i = 0; i < pt.size(); ++i) pt[i] = pa[i] + pb[i];
  const Mat ps = mul(pt, one.s, n);
  res.mulliken.assign(mol.z.size(), 0.0);
  for (std::size_t a = 0; a < mol.z.size(); ++a) res.mulliken[a] = mol.z[a];
  for (int i = 0; i < n; ++i) res.mulliken[static_cast<std::size_t>(basis.atom_of[static_cast<std::size_t>(i)])] -= ps[static_cast<std::size_t>(i) * n + i];
  for (int d = 0; d < 3; ++d) {
    double mu = -dot(pt, one.dipole[static_cast<std::size_t>(d)]);
    for (std::size_t a = 0; a < mol.z.size(); ++a) mu += mol.z[a] * mol.xyz[a][static_cast<std::size_t>(d)];
    res.dipole[static_cast<std::size_t>(d)] = mu * kAuToDebye;
  }
  if (uhf) {  // <S^2> = Sz(Sz+1) + Nb - sum |<a|b>|^2
    const double sz = 0.5 * (na - nb);
    double overlap = 0.0;
    for (int i = 0; i < na; ++i)
      for (int j = 0; j < nb; ++j) {
        double s = 0.0;
        for (int mu = 0; mu < n; ++mu)
          for (int nu = 0; nu < n; ++nu)
            s += ca[static_cast<std::size_t>(mu) * m + i] * one.s[static_cast<std::size_t>(mu) * n + nu] * cb[static_cast<std::size_t>(nu) * m + j];
        overlap += s * s;
      }
    res.s_squared = sz * (sz + 1) + nb - overlap;
  }
  if (!converged)
    diag.warn(DiagCode::NotConverged, "SCF did not converge in " + std::to_string(opt.max_iterations) + " iterations").on_species(sp);
  res.ok = true;
  return res;
}

}  // namespace chemlarp
