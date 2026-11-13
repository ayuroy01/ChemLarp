#include "chemlarp/process.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace chemlarp {
namespace {

std::string num(double v) {
  char b[32];
  std::snprintf(b, sizeof b, "%g", v);
  return b;
}

// Copies the items of `from` into `to`, skipping messages already present.
void merge(const Diagnostics& from, Diagnostics& to) {
  for (const auto& it : from.items()) {
    bool seen = false;
    for (const auto& t : to.items()) seen = seen || (t.code == it.code && t.message == it.message);
    if (seen) continue;
    if (it.severity == Severity::Error) to.error(it.code, it.message);
    else to.warn(it.code, it.message);
  }
}

// Molar enthalpy of component i as ideal-gas vapour.
std::optional<double> h_vapor(const PureComponent& p, double T, Diagnostics& d) {
  auto h = ideal_gas_enthalpy(p, T, d);
  if (!h) return std::nullopt;
  return (std::isfinite(p.Hf_gas_Jmol) ? p.Hf_gas_Jmol : 0.0) + *h;
}

// Molar enthalpy of component i as liquid (see process.hpp).
std::optional<double> h_liquid(const PureComponent& p, double T, Diagnostics& d) {
  if (std::isfinite(p.Tc_K) && T >= p.Tc_K) {
    d.warn(DiagCode::PropertyExtrapolated,
           p.name + " is above its critical temperature in the liquid; its vapour enthalpy is used");
    return h_vapor(p, T, d);
  }
  if (p.liquid_heat_capacity.present() && std::isfinite(p.Tb_K) && std::isfinite(p.Hvap_Tb_Jmol)) {
    Diagnostics q;
    auto hv = h_vapor(p, p.Tb_K, q);
    auto dh = liquid_enthalpy_change(p, p.Tb_K, T, q);
    if (hv && dh) {
      merge(q, d);
      return *hv - p.Hvap_Tb_Jmol + *dh;
    }
  }
  auto hv = h_vapor(p, T, d);
  auto lv = enthalpy_of_vaporization(p, T, d);
  if (!hv || !lv) return std::nullopt;
  d.warn(DiagCode::PropertyExtrapolated, "liquid enthalpy of " + p.name + " from Watson's ΔHvap (no liquid Cp data)");
  return *hv - *lv;
}

// Enthalpy flow of a phase with component flows nf at T.
std::optional<double> phase_enthalpy(const ProcessComponents& c, double T, const std::vector<double>& nf, bool vapor,
                                     Diagnostics& d) {
  double H = 0.0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    if (nf[i] == 0.0) continue;
    auto h = vapor ? h_vapor(c.pure[i], T, d) : h_liquid(c.pure[i], T, d);
    if (!h) return std::nullopt;
    H += nf[i] * *h;
  }
  return H;
}

bool check_flows(const ProcessComponents& c, const std::vector<double>& n, double& F, Diagnostics& d) {
  F = 0.0;
  if (n.size() != c.size()) {
    d.error(DiagCode::InvalidMixtureInput, "expected " + std::to_string(c.size()) + " component flows, got " +
                                               std::to_string(n.size()));
    return false;
  }
  for (double v : n) {
    if (!std::isfinite(v) || v < 0.0) {
      d.error(DiagCode::InvalidMixtureInput, "component flows must be finite and non-negative");
      return false;
    }
    F += v;
  }
  if (F <= 0.0) {
    d.error(DiagCode::InvalidMixtureInput, "the stream has no flow");
    return false;
  }
  return true;
}

bool check_tp(double T, double P, Diagnostics& d) {
  if (std::isfinite(T) && T > 0.0 && std::isfinite(P) && P > 0.0) return true;
  d.error(DiagCode::InvalidMixtureInput, "temperature and pressure must be positive");
  return false;
}

// Rachford–Rice: Σ z_i (K_i − 1) / (1 + β (K_i − 1)) = 0 for β in (0, 1), by bisection on the monotone function.
double rachford_rice(const std::vector<double>& z, const std::vector<double>& K) {
  auto f = [&](double b) {
    double s = 0.0;
    for (std::size_t i = 0; i < z.size(); ++i) s += z[i] * (K[i] - 1.0) / (1.0 + b * (K[i] - 1.0));
    return s;
  };
  double lo = 0.0, hi = 1.0;
  if (f(lo) <= 0.0) return 0.0;
  if (f(hi) >= 0.0) return 1.0;
  for (int k = 0; k < 200 && hi - lo > 1e-15; ++k) {
    const double m = 0.5 * (lo + hi);
    (f(m) > 0.0 ? lo : hi) = m;
  }
  return 0.5 * (lo + hi);
}

// TP phase split without the enthalpy; fills vapor_fraction, x, y.
bool split_tp(const ProcessComponents& c, double T, double P, const std::vector<double>& z, Stream& s, Diagnostics& d) {
  const std::size_t n = c.size();
  s.x = z;
  s.y = z;
  auto bub = bubble_pressure(c.vle, T, z, d);
  if (!bub) return false;
  if (P >= bub->P_Pa) { s.vapor_fraction = 0.0; return true; }
  auto dew = dew_pressure(c.vle, T, z, d);
  if (!dew) return false;
  if (P <= dew->P_Pa) { s.vapor_fraction = 1.0; return true; }
  std::vector<double> psat(n), K(n), x = z;
  for (std::size_t i = 0; i < n; ++i) psat[i] = c.vle.antoine[i].pressure(T);
  for (int it = 0; it < 1000; ++it) {
    auto lg = activity_ln_gamma(c.vle, T, x, d);
    if (!lg) return false;
    for (std::size_t i = 0; i < n; ++i) K[i] = std::exp((*lg)[i]) * psat[i] / P;
    const double b = rachford_rice(z, K);
    std::vector<double> xn(n), yn(n);
    double sx = 0.0, sy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      xn[i] = z[i] / (1.0 + b * (K[i] - 1.0));
      yn[i] = K[i] * xn[i];
      sx += xn[i];
      sy += yn[i];
    }
    double change = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      xn[i] /= sx;
      yn[i] /= sy;
      change = std::max(change, std::fabs(xn[i] - x[i]));
    }
    x = xn;
    s.x = xn;
    s.y = yn;
    s.vapor_fraction = b;
    if (change < 1e-12) return true;
  }
  d.error(DiagCode::VleNotConverged, "two-phase flash at " + num(T) + " K, " + num(P) + " Pa did not converge");
  return false;
}

std::optional<Stream> tp(const ProcessComponents& c, double T, double P, const std::vector<double>& n, Diagnostics& d) {
  double F = 0.0;
  if (!check_flows(c, n, F, d) || !check_tp(T, P, d)) return std::nullopt;
  std::vector<double> z(n.size());
  for (std::size_t i = 0; i < n.size(); ++i) z[i] = n[i] / F;
  Stream s;
  s.T_K = T;
  s.P_Pa = P;
  s.n = n;
  if (!split_tp(c, T, P, z, s, d)) return std::nullopt;
  std::vector<double> nl(n.size()), nv(n.size());
  for (std::size_t i = 0; i < n.size(); ++i) {
    nv[i] = s.vapor_fraction * F * s.y[i];
    nl[i] = (1.0 - s.vapor_fraction) * F * s.x[i];
  }
  auto hl = phase_enthalpy(c, T, nl, false, d);
  auto hv = phase_enthalpy(c, T, nv, true, d);
  if (!hl || !hv) return std::nullopt;
  s.H_W = *hl + *hv;
  return s;
}

std::vector<Stream> phases_of(const Stream& s, const ProcessComponents& c, Diagnostics& d) {
  const double F = s.total();
  Stream v = s, l = s;
  for (std::size_t i = 0; i < s.n.size(); ++i) {
    v.n[i] = s.vapor_fraction * F * s.y[i];
    l.n[i] = (1.0 - s.vapor_fraction) * F * s.x[i];
  }
  v.vapor_fraction = 1.0;
  v.x = v.y;
  l.vapor_fraction = 0.0;
  l.y = l.x;
  v.H_W = phase_enthalpy(c, s.T_K, v.n, true, d).value_or(0.0);
  l.H_W = s.H_W - v.H_W;  // exact split of the flash enthalpy
  return {v, l};
}

}  // namespace

double Stream::total() const noexcept {
  double F = 0.0;
  for (double v : n) F += v;
  return F;
}

std::optional<ProcessComponents> process_components(const std::vector<std::string>& smiles, LiquidModel model,
                                                    Diagnostics& diag) {
  ProcessComponents c;
  c.smiles = smiles;
  bool ok = !smiles.empty();
  if (!ok) diag.error(DiagCode::InvalidMixtureInput, "no components");
  for (const auto& s : smiles) {
    auto p = pure_component(s);
    if (!p) {
      diag.error(DiagCode::NoPropertyData, "no pure-component data for " + s);
      ok = false;
      continue;
    }
    for (double a : p->cp_a)
      if (!std::isfinite(a)) {
        diag.error(DiagCode::NoPropertyData, "no ideal-gas heat capacity for " + p->name);
        ok = false;
        break;
      }
    if (!std::isfinite(p->Hvap_Tb_Jmol) || !std::isfinite(p->Tb_K)) {
      diag.error(DiagCode::NoPropertyData, "no enthalpy of vaporization at the normal boiling point for " + p->name);
      ok = false;
    }
    c.pure.push_back(*p);
  }
  if (!ok) return std::nullopt;
  if (model == LiquidModel::Unifac) {
    auto v = unifac_system_from_smiles(smiles, diag);
    if (!v) return std::nullopt;
    c.vle = *v;
  } else {
    c.vle.model = ActivityModel::Nrtl;  // NRTL with τ = 0: γ = 1
    c.vle.nrtl.n = smiles.size();
    c.vle.nrtl.alpha_c.assign(smiles.size() * smiles.size(), 0.3);
    for (const auto& s : smiles) {
      auto a = antoine_coefficients(s);
      if (!a) {
        diag.error(DiagCode::NoVaporPressureData, "no Antoine vapour-pressure coefficients for " + s);
        return std::nullopt;
      }
      c.vle.antoine.push_back(*a);
    }
  }
  return c;
}

std::optional<Stream> stream_tp(const ProcessComponents& c, double T, double P, const std::vector<double>& n,
                                Diagnostics& diag) {
  Diagnostics q;
  auto s = tp(c, T, P, n, q);
  merge(q, diag);
  return s;
}

std::optional<Stream> stream_ph(const ProcessComponents& c, double P, double H_W, const std::vector<double>& n,
                                Diagnostics& diag, double T_guess) {
  double F = 0.0;
  if (!check_flows(c, n, F, diag) || !check_tp(300.0, P, diag)) return std::nullopt;
  if (!std::isfinite(H_W)) {
    diag.error(DiagCode::InvalidMixtureInput, "enthalpy flow must be finite");
    return std::nullopt;
  }
  // H(T) at fixed P increases with T; find the bracket, then Illinois false position.
  Diagnostics quiet;
  auto g = [&](double T) -> std::optional<double> {
    auto s = tp(c, T, P, n, quiet);
    if (!s) return std::nullopt;
    return s->H_W - H_W;
  };
  double lo = T_guess > 0.0 ? T_guess : 300.0, hi = lo;
  auto glo = g(lo);
  if (!glo) { merge(quiet, diag); return std::nullopt; }
  double flo = *glo, fhi = flo;
  for (int k = 0; k < 200 && flo * fhi > 0.0; ++k) {
    if (flo > 0.0) { hi = lo; fhi = flo; lo = std::max(20.0, lo * 0.9 - 1.0); auto v = g(lo); if (!v) break; flo = *v; }
    else { lo = hi; flo = fhi; hi = hi * 1.1 + 1.0; auto v = g(hi); if (!v) break; fhi = *v; }
    if (lo <= 20.0 && flo > 0.0) break;
  }
  if (flo * fhi > 0.0) {
    merge(quiet, diag);
    diag.error(DiagCode::VleNotConverged, "no temperature found for the requested enthalpy at " + num(P) + " Pa");
    return std::nullopt;
  }
  const double scale = 1e-9 * (std::fabs(H_W) + 1e3 * F);
  int side = 0;
  double T = lo;
  for (int k = 0; k < 300; ++k) {
    T = (lo * fhi - hi * flo) / (fhi - flo);
    auto v = g(T);
    if (!v) { merge(quiet, diag); return std::nullopt; }
    if (std::fabs(*v) <= scale || hi - lo < 1e-10 * T) break;
    if (*v * fhi > 0.0) { hi = T; fhi = *v; if (side == -1) flo /= 2; side = -1; }
    else { lo = T; flo = *v; if (side == 1) fhi /= 2; side = 1; }
  }
  Diagnostics fin;
  auto s = tp(c, T, P, n, fin);
  if (!s) { merge(fin, diag); return std::nullopt; }
  if (std::fabs(s->H_W - H_W) > 1e-6 * (std::fabs(H_W) + 1e3 * F)) {
    // H jumps at T: a pure component or an azeotrope boils at one temperature. Split by the lever rule.
    std::vector<double> z(n.size());
    for (std::size_t i = 0; i < n.size(); ++i) z[i] = n[i] / F;
    auto hl = phase_enthalpy(c, T, n, false, fin), hv = phase_enthalpy(c, T, n, true, fin);
    if (!hl || !hv || *hv <= *hl) { merge(fin, diag); diag.error(DiagCode::VleNotConverged, "PH flash did not converge"); return std::nullopt; }
    s->vapor_fraction = std::clamp((H_W - *hl) / (*hv - *hl), 0.0, 1.0);
    s->x = z;
    s->y = z;
    s->H_W = *hl + s->vapor_fraction * (*hv - *hl);
  }
  merge(fin, diag);
  return s;
}

std::optional<Stream> mixer(const ProcessComponents& c, const std::vector<Stream>& inlets, double P_out, Diagnostics& diag) {
  if (inlets.empty()) {
    diag.error(DiagCode::InvalidMixtureInput, "a mixer needs at least one inlet");
    return std::nullopt;
  }
  std::vector<double> n(c.size(), 0.0);
  double H = 0.0, P = inlets[0].P_Pa, TF = 0.0, F = 0.0;
  for (const Stream& s : inlets) {
    if (s.n.size() != c.size()) {
      diag.error(DiagCode::InvalidMixtureInput, "inlet streams must have one flow per component");
      return std::nullopt;
    }
    for (std::size_t i = 0; i < n.size(); ++i) n[i] += s.n[i];
    H += s.H_W;
    P = std::min(P, s.P_Pa);
    TF += s.T_K * s.total();
    F += s.total();
  }
  return stream_ph(c, P_out > 0.0 ? P_out : P, H, n, diag, F > 0.0 ? TF / F : 300.0);
}

std::optional<std::vector<Stream>> splitter(const ProcessComponents& c, const Stream& in, const std::vector<double>& fractions,
                                            Diagnostics& diag) {
  double sum = 0.0;
  for (double f : fractions) {
    if (!std::isfinite(f) || f < 0.0) {
      diag.error(DiagCode::InvalidMixtureInput, "split fractions must be non-negative");
      return std::nullopt;
    }
    sum += f;
  }
  if (fractions.empty() || std::fabs(sum - 1.0) > 1e-9 || in.n.size() != c.size()) {
    diag.error(DiagCode::InvalidMixtureInput, "split fractions must sum to 1");
    return std::nullopt;
  }
  std::vector<Stream> out;
  for (double f : fractions) {
    Stream s = in;
    for (double& v : s.n) v *= f;
    s.H_W *= f;
    out.push_back(s);
  }
  return out;
}

std::optional<HeaterResult> heater_to_temperature(const ProcessComponents& c, const Stream& in, double T_out, double P_out,
                                                  Diagnostics& diag) {
  auto s = stream_tp(c, T_out, P_out > 0.0 ? P_out : in.P_Pa, in.n, diag);
  if (!s) return std::nullopt;
  return HeaterResult{*s, s->H_W - in.H_W};
}

std::optional<HeaterResult> heater_with_duty(const ProcessComponents& c, const Stream& in, double duty_W, double P_out,
                                             Diagnostics& diag) {
  auto s = stream_ph(c, P_out > 0.0 ? P_out : in.P_Pa, in.H_W + duty_W, in.n, diag, in.T_K);
  if (!s) return std::nullopt;
  return HeaterResult{*s, duty_W};
}

std::optional<Stream> valve(const ProcessComponents& c, const Stream& in, double P_out, Diagnostics& diag) {
  if (!(P_out > 0.0) || P_out > in.P_Pa) {
    diag.error(DiagCode::InvalidMixtureInput, "a valve needs an outlet pressure between 0 and the inlet pressure");
    return std::nullopt;
  }
  return stream_ph(c, P_out, in.H_W, in.n, diag, in.T_K);
}

std::optional<FlashDrumResult> flash_drum(const ProcessComponents& c, const Stream& in, double T, double P, Diagnostics& diag) {
  auto s = stream_tp(c, T, P, in.n, diag);
  if (!s) return std::nullopt;
  auto ph = phases_of(*s, c, diag);
  return FlashDrumResult{ph[0], ph[1], s->H_W - in.H_W};
}

std::optional<FlashDrumResult> flash_drum_adiabatic(const ProcessComponents& c, const Stream& in, double P, Diagnostics& diag) {
  auto s = stream_ph(c, P, in.H_W, in.n, diag, in.T_K);
  if (!s) return std::nullopt;
  auto ph = phases_of(*s, c, diag);
  return FlashDrumResult{ph[0], ph[1], 0.0};
}

}  // namespace chemlarp
