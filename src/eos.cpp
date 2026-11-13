#include "chemlarp/eos.hpp"
#include "chemlarp/properties.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "chemlarp/smiles.hpp"
#include "chemlarp/system_state.hpp"
#include "chemlarp/thermo.hpp"

namespace chemlarp {
namespace {

constexpr double kSqrt2 = 1.4142135623730951;

struct Constants {
  double omega_a, omega_b, d1, d2;  // Z + d1 B, Z + d2 B in the attractive term
};

Constants constants(CubicEos e) {
  if (e == CubicEos::PengRobinson) return {0.45723553, 0.07779607, 1.0 + kSqrt2, 1.0 - kSqrt2};
  return {0.42748023, 0.08664035, 1.0, 0.0};
}

double kappa(CubicEos e, double w) {
  if (e == CubicEos::SoaveRedlichKwong) return 0.480 + 1.574 * w - 0.176 * w * w;
  if (w <= 0.491) return 0.37464 + 1.54226 * w - 0.26992 * w * w;
  return 0.379642 + 1.48503 * w - 0.164423 * w * w + 0.016666 * w * w * w;  // PR 1978
}

// Pure-component a_i(T) and b_i.
void pure_ab(const EosMixture& m, double T, std::vector<double>& a, std::vector<double>& b) {
  const Constants c = constants(m.eos);
  const std::size_t n = m.size();
  a.resize(n);
  b.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto& p = m.components[i];
    const double tr = T / p.Tc_K;
    const double s = 1.0 + kappa(m.eos, p.omega) * (1.0 - std::sqrt(tr));
    a[i] = c.omega_a * kGasConstant * kGasConstant * p.Tc_K * p.Tc_K / p.Pc_Pa * s * s;
    b[i] = c.omega_b * kGasConstant * p.Tc_K / p.Pc_Pa;
  }
}

struct MixAB {
  double a = 0, b = 0;
  std::vector<double> sum_xa;  // Σ_j x_j a_ij
  std::vector<double> bi;
};

MixAB mix_ab(const EosMixture& m, double T, const std::vector<double>& x) {
  std::vector<double> a, b;
  pure_ab(m, T, a, b);
  const std::size_t n = m.size();
  MixAB r;
  r.sum_xa.assign(n, 0.0);
  r.bi = b;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      const double kij = m.kij.size() == n * n ? m.kij[i * n + j] : 0.0;
      const double aij = std::sqrt(a[i] * a[j]) * (1.0 - kij);
      r.sum_xa[i] += x[j] * aij;
    }
    r.a += x[i] * r.sum_xa[i];
    r.b += x[i] * b[i];
  }
  return r;
}

// Real roots of z³ + p2 z² + p1 z + p0, ascending.
std::vector<double> cubic_roots(double p2, double p1, double p0) {
  const double q = (3 * p1 - p2 * p2) / 9.0;
  const double r = (9 * p2 * p1 - 27 * p0 - 2 * p2 * p2 * p2) / 54.0;
  const double disc = q * q * q + r * r;
  std::vector<double> z;
  if (disc > 0) {
    const double s = std::cbrt(r + std::sqrt(disc)), t = std::cbrt(r - std::sqrt(disc));
    z.push_back(s + t - p2 / 3.0);
  } else {
    const double th = std::acos(std::clamp(r / std::sqrt(-q * q * q), -1.0, 1.0));
    const double m = 2.0 * std::sqrt(-q);
    for (int k = 0; k < 3; ++k) z.push_back(m * std::cos((th + 2.0 * 3.14159265358979323846 * k) / 3.0) - p2 / 3.0);
  }
  // One Newton polish per root.
  for (double& v : z) {
    const double f = ((v + p2) * v + p1) * v + p0, df = (3 * v + 2 * p2) * v + p1;
    if (df != 0) v -= f / df;
  }
  std::sort(z.begin(), z.end());
  return z;
}

std::vector<double> normalised(std::vector<double> v) {
  double s = 0;
  for (double x : v) s += x;
  for (double& x : v) x /= s;
  return v;
}

double rachford_rice(const std::vector<double>& z, const std::vector<double>& K, double& beta) {
  // Bracket of physical solutions keeps every denominator positive.
  double kmax = *std::max_element(K.begin(), K.end()), kmin = *std::min_element(K.begin(), K.end());
  double lo = 1.0 / (1.0 - kmax), hi = 1.0 / (1.0 - kmin);
  auto f = [&](double b) {
    double s = 0;
    for (std::size_t i = 0; i < z.size(); ++i) s += z[i] * (K[i] - 1) / (1 + b * (K[i] - 1));
    return s;
  };
  lo += 1e-12 * std::fabs(lo);
  hi -= 1e-12 * std::fabs(hi);
  beta = std::clamp(beta, lo, hi);
  for (int it = 0; it < 100; ++it) {
    const double fv = f(beta);
    double df = 0;
    for (std::size_t i = 0; i < z.size(); ++i) {
      const double d = 1 + beta * (K[i] - 1);
      df -= z[i] * (K[i] - 1) * (K[i] - 1) / (d * d);
    }
    if (fv > 0) lo = beta; else hi = beta;
    double nb = beta - fv / df;
    if (!(nb > lo && nb < hi)) nb = 0.5 * (lo + hi);
    if (std::fabs(nb - beta) < 1e-14) { beta = nb; break; }
    beta = nb;
  }
  return f(beta);
}

}  // namespace

std::vector<double> compressibility_roots(const EosMixture& m, double T, double P, const std::vector<double>& x) {
  const Constants c = constants(m.eos);
  const MixAB ab = mix_ab(m, T, x);
  const double A = ab.a * P / (kGasConstant * kGasConstant * T * T), B = ab.b * P / (kGasConstant * T);
  const double u = c.d1 + c.d2, w = c.d1 * c.d2;
  // Z³ + (uB - B - 1) Z² + (A + wB² - uB - uB²) Z - (AB + wB² + wB³) = 0
  auto roots = cubic_roots(u * B - B - 1.0, A + w * B * B - u * B - u * B * B, -(A * B + w * B * B + w * B * B * B));
  std::vector<double> out;
  for (double z : roots)
    if (z > B) out.push_back(z);
  return out;
}

PhaseProperties phase_properties(const EosMixture& m, double T, double P, const std::vector<double>& x, Root root) {
  const Constants c = constants(m.eos);
  const MixAB ab = mix_ab(m, T, x);
  const double RT = kGasConstant * T;
  const double A = ab.a * P / (RT * RT), B = ab.b * P / RT;
  const auto roots = compressibility_roots(m, T, P, x);
  auto props = [&](double Z) {
    PhaseProperties p;
    p.Z = Z;
    p.molar_volume = Z * RT / P;
    p.ln_phi.resize(m.size());
    const double L = std::log((Z + c.d1 * B) / (Z + c.d2 * B));
    for (std::size_t i = 0; i < m.size(); ++i) {
      const double bb = ab.bi[i] / ab.b;
      p.ln_phi[i] = bb * (Z - 1) - std::log(Z - B) - A / ((c.d1 - c.d2) * B) * (2 * ab.sum_xa[i] / ab.a - bb) * L;
      p.gibbs_residual_RT += x[i] * p.ln_phi[i];
    }
    return p;
  };
  if (roots.empty()) return props(1.0);
  if (root == Root::Liquid) return props(roots.front());
  if (root == Root::Vapor) return props(roots.back());
  PhaseProperties best = props(roots.front());
  if (roots.size() > 1) {
    PhaseProperties v = props(roots.back());
    if (v.gibbs_residual_RT < best.gibbs_residual_RT) best = v;
  }
  return best;
}

std::optional<double> saturation_pressure(const CriticalProperties& comp, CubicEos eos, double T, Diagnostics& diag) {
  if (!(T < comp.Tc_K)) {
    diag.error(DiagCode::SupercriticalTemperature, comp.name + ": no vapour pressure above Tc = " + std::to_string(comp.Tc_K) + " K");
    return std::nullopt;
  }
  EosMixture m;
  m.eos = eos;
  m.components = {comp};
  const std::vector<double> x{1.0};
  double P = comp.Pc_Pa * std::exp(5.373 * (1 + comp.omega) * (1 - comp.Tc_K / T));  // Wilson estimate
  for (int it = 0; it < 200; ++it) {
    auto roots = compressibility_roots(m, T, P, x);
    if (roots.size() < 2) {  // only one root: step towards the two-root region
      const double Z = roots.empty() ? 1.0 : roots.front();
      P *= (Z > 0.3) ? 1.1 : 0.9;
      continue;
    }
    const double lphiL = phase_properties(m, T, P, x, Root::Liquid).ln_phi[0];
    const double lphiV = phase_properties(m, T, P, x, Root::Vapor).ln_phi[0];
    const double d = lphiL - lphiV;
    P *= std::exp(d);
    if (std::fabs(d) < 1e-12) return P;
  }
  diag.error(DiagCode::FlashNotConverged, comp.name + ": vapour pressure iteration did not converge at " + std::to_string(T) + " K");
  return std::nullopt;
}

const char* to_string(PhaseRegion r) noexcept {
  switch (r) {
    case PhaseRegion::Liquid: return "liquid";
    case PhaseRegion::Vapor: return "vapor";
    case PhaseRegion::Supercritical: return "supercritical";
    case PhaseRegion::TwoPhase: return "two-phase";
  }
  return "?";
}

FlashResult flash_tp(const EosMixture& m, double T, double P, const std::vector<double>& z_in, Diagnostics& diag) {
  FlashResult res;
  const std::size_t n = m.size();
  if (n == 0 || z_in.size() != n || !(T > 0) || !(P > 0)) {
    diag.error(DiagCode::InvalidInput, "flash needs one mole fraction per component, T > 0 and P > 0");
    return res;
  }
  for (double v : z_in)
    if (!(v >= 0) || !std::isfinite(v)) {
      diag.error(DiagCode::NegativeAmount, "feed mole fractions must be finite and >= 0");
      return res;
    }
  const std::vector<double> z = normalised(z_in);
  res.ok = true;

  // Michelsen stability test: trial phases from Wilson K (vapour-like and liquid-like).
  const PhaseProperties feed = phase_properties(m, T, P, z, Root::Stable);
  std::vector<double> d(n), Kw(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto& c = m.components[i];
    Kw[i] = c.Pc_Pa / P * std::exp(5.373 * (1 + c.omega) * (1 - c.Tc_K / T));
    d[i] = (z[i] > 0 ? std::log(z[i]) : -700.0) + feed.ln_phi[i];
  }
  std::vector<double> best_K;
  double best_tpd = 0.0;
  for (int trial = 0; trial < 2; ++trial) {
    std::vector<double> W(n);
    for (std::size_t i = 0; i < n; ++i) W[i] = trial == 0 ? z[i] * Kw[i] : z[i] / Kw[i];
    for (int it = 0; it < 300; ++it) {
      const auto ph = phase_properties(m, T, P, normalised(W), Root::Stable);
      double change = 0;
      for (std::size_t i = 0; i < n; ++i) {
        const double nw = std::exp(d[i] - ph.ln_phi[i]);
        change = std::max(change, std::fabs(std::log(std::max(nw, 1e-300)) - std::log(std::max(W[i], 1e-300))));
        W[i] = nw;
      }
      if (change < 1e-10) break;
    }
    double sw = 0;
    for (double v : W) sw += v;
    const double tpd = 1.0 - sw;  // < 0 means the feed is unstable
    bool trivial = true;
    const auto wn = normalised(W);
    for (std::size_t i = 0; i < n; ++i) trivial = trivial && std::fabs(wn[i] - z[i]) < 1e-6;
    if (!trivial && tpd < best_tpd - 1e-10) {
      best_tpd = tpd;
      best_K.resize(n);
      for (std::size_t i = 0; i < n; ++i) best_K[i] = trial == 0 ? std::max(wn[i], 1e-300) / std::max(z[i], 1e-300)
                                                                 : std::max(z[i], 1e-300) / std::max(wn[i], 1e-300);
    }
  }

  auto single = [&]() {
    res.vapor_fraction = 0;
    res.x = res.y = z;
    res.K.assign(n, 1.0);
    res.liquid = res.vapor = feed;
    // Single phase: supercritical above the pseudo-critical point (Kay's rule),
    // otherwise liquid-like when v/b < 1.75 (Pedersen et al.), else vapour.
    double Tpc = 0, Ppc = 0;
    for (std::size_t i = 0; i < n; ++i) {
      Tpc += z[i] * m.components[i].Tc_K;
      Ppc += z[i] * m.components[i].Pc_Pa;
    }
    const double b = mix_ab(m, T, z).b;
    if (T > Tpc && P > Ppc) res.region = PhaseRegion::Supercritical;
    else res.region = feed.molar_volume / b < 1.75 ? PhaseRegion::Liquid : PhaseRegion::Vapor;
    if (res.region == PhaseRegion::Vapor) res.vapor_fraction = 1.0;
    res.converged = true;
    return res;
  };
  if (best_K.empty()) return single();

  // Two-phase flash by successive substitution.
  std::vector<double> K = best_K, x(n), y(n);
  double beta = 0.5;
  for (res.iterations = 1; res.iterations <= 500; ++res.iterations) {
    rachford_rice(z, K, beta);
    for (std::size_t i = 0; i < n; ++i) {
      x[i] = z[i] / (1 + beta * (K[i] - 1));
      y[i] = K[i] * x[i];
    }
    x = normalised(x);
    y = normalised(y);
    res.liquid = phase_properties(m, T, P, x, Root::Liquid);
    res.vapor = phase_properties(m, T, P, y, Root::Vapor);
    double change = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const double nk = std::exp(res.liquid.ln_phi[i] - res.vapor.ln_phi[i]);
      change = std::max(change, std::fabs(std::log(nk / K[i])));
      K[i] = nk;
    }
    if (change < 1e-11) { res.converged = true; break; }
  }
  if (!res.converged) diag.warn(DiagCode::FlashNotConverged, "two-phase flash did not converge in 500 iterations (near-critical?)");
  if (beta <= 0.0 || beta >= 1.0 || !std::isfinite(beta)) return single();
  bool trivial = true;
  for (std::size_t i = 0; i < n; ++i) trivial = trivial && std::fabs(K[i] - 1.0) < 1e-6;
  if (trivial) return single();
  res.region = PhaseRegion::TwoPhase;
  res.vapor_fraction = beta;
  res.x = x;
  res.y = y;
  res.K = K;
  return res;
}

// ---- critical constants ------------------------------------------------------

namespace {
struct Row {
  const char* smiles;
  const char* name;
  double Tc, Pc_bar, omega;
};

// Poling, Prausnitz & O'Connell (2001), Appendix A.
constexpr Row kCritical[] = {
    {"C", "methane", 190.56, 45.99, 0.011},
    {"CC", "ethane", 305.32, 48.72, 0.099},
    {"CCC", "propane", 369.83, 42.48, 0.152},
    {"CCCC", "n-butane", 425.12, 37.96, 0.200},
    {"CC(C)C", "isobutane", 407.85, 36.40, 0.186},
    {"CCCCC", "n-pentane", 469.70, 33.70, 0.252},
    {"CCCCCC", "n-hexane", 507.60, 30.25, 0.301},
    {"CCCCCCC", "n-heptane", 540.20, 27.40, 0.350},
    {"CCCCCCCC", "n-octane", 568.70, 24.90, 0.399},
    {"C=C", "ethylene", 282.34, 50.41, 0.087},
    {"CC=C", "propylene", 364.90, 46.00, 0.142},
    {"c1ccccc1", "benzene", 562.05, 48.95, 0.210},
    {"Cc1ccccc1", "toluene", 591.75, 41.08, 0.264},
    {"C1CCCCC1", "cyclohexane", 553.50, 40.73, 0.211},
    {"CO", "methanol", 512.64, 80.97, 0.565},
    {"CCO", "ethanol", 513.92, 61.48, 0.649},
    {"CC(C)=O", "acetone", 508.20, 47.01, 0.307},
    {"O", "water", 647.14, 220.64, 0.344},
    {"N", "ammonia", 405.65, 112.80, 0.253},
    {"O=C=O", "carbon dioxide", 304.12, 73.74, 0.225},
    {"[C-]#[O+]", "carbon monoxide", 132.85, 34.94, 0.045},
    {"N#N", "nitrogen", 126.20, 33.98, 0.037},
    {"O=O", "oxygen", 154.58, 50.43, 0.022},
    {"[H][H]", "hydrogen", 33.19, 13.13, -0.216},
    {"[Ar]", "argon", 150.86, 48.98, -0.002},
    {"[He]", "helium", 5.19, 2.27, -0.390},
    {"S", "hydrogen sulfide", 373.40, 89.63, 0.090},
    {"O=S=O", "sulfur dioxide", 430.75, 78.84, 0.245},
    {"ClCl", "chlorine", 416.90, 79.80, 0.090},
};

const std::map<std::string, CriticalProperties>& table() {
  static const std::map<std::string, CriticalProperties> t = [] {
    std::map<std::string, CriticalProperties> m;
    for (const Row& r : kCritical) {
      ChemicalSystemState st;
      Diagnostics d;
      const std::int32_t sp = parse_smiles(st, r.smiles, d);
      if (sp < 0) continue;
      m[st.species[static_cast<std::size_t>(sp)].canonical_smiles] = {r.name, r.Tc, r.Pc_bar * 1e5, r.omega};
    }
    return m;
  }();
  return t;
}
}  // namespace

std::optional<CriticalProperties> critical_properties(std::string_view smiles) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smiles, d);
  if (sp < 0) return std::nullopt;
  auto it = table().find(st.species[static_cast<std::size_t>(sp)].canonical_smiles);
  if (it != table().end()) return it->second;
  // Otherwise the pure-component database, when it has all three constants.
  if (auto c = pure_component(smiles); c && std::isfinite(c->Tc_K) && std::isfinite(c->Pc_Pa) && std::isfinite(c->omega))
    return CriticalProperties{c->name, c->Tc_K, c->Pc_Pa, c->omega};
  return std::nullopt;
}

std::optional<EosMixture> mixture_from_smiles(const std::vector<std::string>& smiles, CubicEos eos, Diagnostics& diag) {
  EosMixture m;
  m.eos = eos;
  bool ok = true;
  for (const auto& s : smiles) {
    ChemicalSystemState st;
    const std::int32_t sp = parse_smiles(st, s, diag);
    if (sp < 0) { ok = false; continue; }
    auto c = critical_properties(s);
    if (!c) {
      diag.error(DiagCode::NoCriticalData, "no critical constants for " + s + " (" +
                                               st.species[static_cast<std::size_t>(sp)].formula + ")");
      ok = false;
      continue;
    }
    m.components.push_back(*c);
  }
  if (!ok) return std::nullopt;
  return m;
}

}  // namespace chemlarp
