#include "chemlarp/process.hpp"

#include "chemlarp/formula.hpp"
#include "chemlarp/smiles.hpp"
#include "chemlarp/system_state.hpp"

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
  if (!hv) return std::nullopt;
  Diagnostics q;
  if (auto lv = enthalpy_of_vaporization(p, T, q)) {
    d.warn(DiagCode::PropertyExtrapolated, "liquid enthalpy of " + p.name + " from Watson's ΔHvap (no liquid Cp data)");
    return *hv - *lv;
  }
  // No enthalpy of vaporization (a gas such as CO2 without a normal boiling point): only traces dissolve
  // under the modified-Raoult model, and they are given their vapour enthalpy.
  d.warn(DiagCode::PropertyExtrapolated, "no enthalpy of vaporization for " + p.name +
                                             "; in the liquid it is given its vapour enthalpy");
  return hv;
}

// Pressure term of the liquid enthalpy, v_L (P − P_sat): small, and left out (silently) where the liquid
// volume or the vapour pressure is not defined (above Tc).
double liquid_pressure_term(const ProcessComponents& c, std::size_t i, double T, double P) {
  const PureComponent& p = c.pure[i];
  if (!(std::isfinite(p.Tc_K) && T < p.Tc_K)) return 0.0;
  Diagnostics q;
  auto v = saturated_liquid_volume(p, T, q);
  if (!v) return 0.0;
  for (std::size_t k = 0; k < c.condensable.size(); ++k)
    if (c.condensable[k] == i) return *v * (P - c.vle.antoine[k].pressure(T));
  return 0.0;
}

// Enthalpy flow of a phase with component flows nf at T and P.
std::optional<double> phase_enthalpy(const ProcessComponents& c, double T, double P, const std::vector<double>& nf,
                                     bool vapor, Diagnostics& d) {
  double H = 0.0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    if (nf[i] == 0.0) continue;
    auto h = vapor ? h_vapor(c.pure[i], T, d) : h_liquid(c.pure[i], T, d);
    if (!h) return std::nullopt;
    H += nf[i] * (*h + (vapor ? 0.0 : liquid_pressure_term(c, i, T, P)));
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
  const std::size_t n = c.size(), m = c.condensable.size();
  s.x = z;
  s.y = z;
  double znc = 1.0;
  std::vector<double> zc(m);
  for (std::size_t k = 0; k < m; ++k) {
    zc[k] = z[c.condensable[k]];
    znc -= zc[k];
  }
  znc = std::max(znc, 0.0);
  const bool gas = znc > 1e-14;
  if (m == 0 || znc >= 1.0 - 1e-14) { s.vapor_fraction = 1.0; return true; }  // nothing can condense
  if (!gas) {
    // Only condensable components: one liquid below the bubble pressure, one vapour above the dew pressure.
    auto bub = bubble_pressure(c.vle, T, zc, d);
    if (!bub) return false;
    if (P >= bub->P_Pa) { s.vapor_fraction = 0.0; return true; }
    auto dew = dew_pressure(c.vle, T, zc, d);
    if (!dew) return false;
    if (P <= dew->P_Pa) { s.vapor_fraction = 1.0; return true; }
  }
  // Rachford–Rice over the condensable components (K from γ(x) P_sat / P), plus non-condensables with K = ∞:
  // f(β) = Σ_c z (K − 1)/(1 + β (K − 1)) + z_nc / β, decreasing in β.
  std::vector<double> psat(m), K(m), x(m);
  for (std::size_t k = 0; k < m; ++k) {
    psat[k] = c.vle.antoine[k].pressure(T);
    x[k] = zc[k] / (1.0 - znc);
  }
  for (int it = 0; it < 1000; ++it) {
    auto lg = activity_ln_gamma(c.vle, T, x, d);
    if (!lg) return false;
    for (std::size_t k = 0; k < m; ++k) K[k] = std::exp((*lg)[k]) * psat[k] / P;
    double b;
    if (!gas) {
      b = rachford_rice(zc, K);
    } else {
      auto f = [&](double bb) {
        double v = znc / bb;
        for (std::size_t k = 0; k < m; ++k) v += zc[k] * (K[k] - 1.0) / (1.0 + bb * (K[k] - 1.0));
        return v;
      };
      if (f(1.0) >= 0.0) {
        b = 1.0;
      } else {
        double lo = 0.0, hi = 1.0;
        for (int k = 0; k < 200 && hi - lo > 1e-15; ++k) {
          const double mid = 0.5 * (lo + hi);
          (f(mid) > 0.0 ? lo : hi) = mid;
        }
        b = 0.5 * (lo + hi);
      }
    }
    std::vector<double> xn(m), yc(m);
    double sx = 0.0;
    for (std::size_t k = 0; k < m; ++k) {
      xn[k] = zc[k] / (1.0 + b * (K[k] - 1.0));
      yc[k] = K[k] * xn[k];
      sx += xn[k];
    }
    double change = 0.0;
    for (std::size_t k = 0; k < m; ++k) {
      xn[k] /= sx;
      change = std::max(change, std::fabs(xn[k] - x[k]));
    }
    x = xn;
    // Full-length compositions: no non-condensable in the liquid; vapour from the material balance.
    std::fill(s.x.begin(), s.x.end(), 0.0);
    for (std::size_t k = 0; k < m; ++k) s.x[c.condensable[k]] = xn[k];
    double sy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      s.y[i] = b > 0.0 ? (z[i] - (1.0 - b) * s.x[i]) / b : 0.0;
      s.y[i] = std::max(s.y[i], 0.0);
      sy += s.y[i];
    }
    for (double& v : s.y) v /= sy;
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
  auto hl = phase_enthalpy(c, T, P, nl, false, d);
  auto hv = phase_enthalpy(c, T, P, nv, true, d);
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
  v.H_W = phase_enthalpy(c, s.T_K, s.P_Pa, v.n, true, d).value_or(0.0);
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
    c.pure.push_back(*p);
  }
  if (!ok) return std::nullopt;
  // Components without a vapour-pressure equation are non-condensable: they stay in the vapour.
  std::vector<std::string> cond;
  for (std::size_t i = 0; i < smiles.size(); ++i) {
    if (antoine_coefficients(smiles[i])) {
      c.condensable.push_back(i);
      cond.push_back(smiles[i]);
    } else {
      diag.warn(DiagCode::NoVaporPressureData,
                "no vapour-pressure equation for " + c.pure[i].name + "; it is treated as non-condensable");
    }
  }
  if (cond.empty()) return c;
  if (model == LiquidModel::Unifac) {
    auto v = unifac_system_from_smiles(cond, diag);
    if (!v) return std::nullopt;
    c.vle = *v;
  } else {
    c.vle.model = ActivityModel::Nrtl;  // NRTL with τ = 0: γ = 1
    c.vle.nrtl.n = cond.size();
    c.vle.nrtl.alpha_c.assign(cond.size() * cond.size(), 0.3);
    for (const auto& s : cond) c.vle.antoine.push_back(*antoine_coefficients(s));
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
  const double scale = 1e-9 * (std::fabs(H_W) + 1e3 * F);
  double lo = T_guess > 0.0 ? T_guess : 300.0, hi = lo;
  auto glo = g(lo);
  if (!glo) { merge(quiet, diag); return std::nullopt; }
  double flo = *glo, fhi = flo;
  // The guess may already be the answer (a mixer fed with the phases of one flash); otherwise widen the
  // bracket until the residual changes sign. A collapsed bracket would make false position divide 0 by 0.
  double T = lo;
  bool done = std::fabs(flo) <= scale;
  for (int k = 0; k < 200 && !done && flo * fhi > 0.0; ++k) {
    std::optional<double> v;
    if (flo > 0.0) {
      hi = lo; fhi = flo;
      if (lo <= 20.0) break;
      lo = std::max(20.0, lo * 0.9 - 1.0);
      v = g(lo);
      if (!v) break;
      flo = *v;
      if (std::fabs(flo) <= scale) { T = lo; done = true; }
    } else {
      lo = hi; flo = fhi;
      hi = hi * 1.1 + 1.0;
      v = g(hi);
      if (!v) break;
      fhi = *v;
      if (std::fabs(fhi) <= scale) { T = hi; done = true; }
    }
  }
  if (!done && !(flo < 0.0 && fhi > 0.0)) {
    merge(quiet, diag);
    diag.error(DiagCode::VleNotConverged, "no temperature found for the requested enthalpy at " + num(P) + " Pa");
    return std::nullopt;
  }
  // Illinois false position on the bracket lo < hi with g(lo) < 0 < g(hi).
  for (int k = 0, side = 0; k < 300 && !done; ++k) {
    T = (lo * fhi - hi * flo) / (fhi - flo);
    auto v = g(T);
    if (!v) { merge(quiet, diag); return std::nullopt; }
    if (std::fabs(*v) <= scale || hi - lo < 1e-10 * T) break;
    if (*v > 0.0) { hi = T; fhi = *v; if (side == -1) flo /= 2; side = -1; }
    else { lo = T; flo = *v; if (side == 1) fhi /= 2; side = 1; }
  }
  Diagnostics fin;
  auto s = tp(c, T, P, n, fin);
  if (!s) { merge(fin, diag); return std::nullopt; }
  if (std::fabs(s->H_W - H_W) > 1e-6 * (std::fabs(H_W) + 1e3 * F)) {
    // H jumps at T: a pure component or an azeotrope boils at one temperature. Split by the lever rule.
    std::vector<double> z(n.size());
    for (std::size_t i = 0; i < n.size(); ++i) z[i] = n[i] / F;
    auto hl = phase_enthalpy(c, T, P, n, false, fin), hv = phase_enthalpy(c, T, P, n, true, fin);
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

namespace {

// ∫ Cp_ig/T dT antiderivative (J/(mol K)) of the Poling polynomial.
double cp_over_t(const PureComponent& p, double T) {
  const double* a = p.cp_a;
  return 8.314462618 * (a[0] * std::log(T) + T * (a[1] + T * (a[2] / 2 + T * (a[3] / 3 + T * a[4] / 4))));
}

std::optional<double> ph_temperature(const ProcessComponents& c, double P, double H, const std::vector<double>& n,
                                     double guess, Diagnostics& d) {
  auto s = stream_ph(c, P, H, n, d, guess);
  if (!s) return std::nullopt;
  return s->T_K;
}

}  // namespace

std::optional<HeatExchangerResult> heat_exchanger(const ProcessComponents& hot, const Stream& hot_in,
                                                  const ProcessComponents& cold, const Stream& cold_in,
                                                  ExchangerSpec spec, double value, Diagnostics& diag) {
  if (hot_in.T_K <= cold_in.T_K) {
    diag.error(DiagCode::TemperatureCross, "the hot inlet (" + num(hot_in.T_K) + " K) is not hotter than the cold inlet (" +
                                               num(cold_in.T_K) + " K)");
    return std::nullopt;
  }
  double Q = 0.0;
  if (spec == ExchangerSpec::HotOutletTemperature) {
    auto o = stream_tp(hot, value, hot_in.P_Pa, hot_in.n, diag);
    if (!o) return std::nullopt;
    Q = hot_in.H_W - o->H_W;
  } else if (spec == ExchangerSpec::ColdOutletTemperature) {
    auto o = stream_tp(cold, value, cold_in.P_Pa, cold_in.n, diag);
    if (!o) return std::nullopt;
    Q = o->H_W - cold_in.H_W;
  } else {
    Q = value;
  }
  if (!(std::isfinite(Q) && Q > 0.0)) {
    diag.error(DiagCode::InvalidMixtureInput, "the specification gives no heat flow from the hot to the cold stream");
    return std::nullopt;
  }
  HeatExchangerResult r;
  r.duty_W = Q;
  auto ho = stream_ph(hot, hot_in.P_Pa, hot_in.H_W - Q, hot_in.n, diag, hot_in.T_K);
  auto co = stream_ph(cold, cold_in.P_Pa, cold_in.H_W + Q, cold_in.n, diag, cold_in.T_K);
  if (!ho || !co) return std::nullopt;
  r.hot_out = *ho;
  r.cold_out = *co;
  // Counter-current profile: position f = 0 is the hot inlet and the cold outlet.
  constexpr int kZones = 20;
  std::vector<double> dT(kZones + 1);
  for (int k = 0; k <= kZones; ++k) {
    const double f = static_cast<double>(k) / kZones;
    Diagnostics q;
    auto th = k == 0 ? std::optional<double>(hot_in.T_K) : k == kZones ? std::optional<double>(ho->T_K)
                     : ph_temperature(hot, hot_in.P_Pa, hot_in.H_W - f * Q, hot_in.n, hot_in.T_K, q);
    auto tc = k == 0 ? std::optional<double>(co->T_K) : k == kZones ? std::optional<double>(cold_in.T_K)
                     : ph_temperature(cold, cold_in.P_Pa, cold_in.H_W + (1.0 - f) * Q, cold_in.n, cold_in.T_K, q);
    if (!th || !tc) {
      merge(q, diag);
      return std::nullopt;
    }
    dT[k] = *th - *tc;
  }
  r.min_approach_K = *std::min_element(dT.begin(), dT.end());
  if (r.min_approach_K <= 0.0) {
    diag.error(DiagCode::TemperatureCross, "temperature cross inside the exchanger: the hot stream falls " +
                                               num(-r.min_approach_K) + " K below the cold stream");
    return std::nullopt;
  }
  for (int k = 0; k < kZones; ++k) {
    const double a = dT[k], b = dT[k + 1];
    const double lm = std::fabs(a - b) < 1e-9 * (a + b) ? 0.5 * (a + b) : (a - b) / std::log(a / b);
    r.UA_WK += (Q / kZones) / lm;
  }
  return r;
}

std::optional<WorkResult> pump(const ProcessComponents& c, const Stream& in, double P_out, double efficiency, Diagnostics& diag) {
  if (in.vapor_fraction > 0.0) {
    diag.error(DiagCode::InvalidMixtureInput, "a pump needs a liquid feed (vapour fraction " + num(in.vapor_fraction) + ")");
    return std::nullopt;
  }
  if (!(P_out >= in.P_Pa) || !(efficiency > 0.0 && efficiency <= 1.0)) {
    diag.error(DiagCode::InvalidMixtureInput, "a pump needs an outlet pressure at or above the inlet and an efficiency in (0, 1]");
    return std::nullopt;
  }
  double V = 0.0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    if (in.n[i] == 0.0) continue;
    auto v = saturated_liquid_volume(c.pure[i], in.T_K, diag);
    if (!v) return std::nullopt;
    V += in.n[i] * *v;
  }
  WorkResult r;
  r.isentropic_work_W = V * (P_out - in.P_Pa);
  r.work_W = r.isentropic_work_W / efficiency;
  auto o = stream_ph(c, P_out, in.H_W + r.work_W, in.n, diag, in.T_K);
  if (!o) return std::nullopt;
  r.out = *o;
  return r;
}

std::optional<WorkResult> compressor(const ProcessComponents& c, const Stream& in, double P_out, double efficiency,
                                     Diagnostics& diag) {
  if (in.vapor_fraction < 1.0) {
    diag.error(DiagCode::InvalidMixtureInput, "a compressor needs a vapour feed (vapour fraction " + num(in.vapor_fraction) + ")");
    return std::nullopt;
  }
  if (!(P_out > 0.0) || !(efficiency > 0.0 && efficiency <= 1.0)) {
    diag.error(DiagCode::InvalidMixtureInput, "a compressor needs a positive outlet pressure and an efficiency in (0, 1]");
    return std::nullopt;
  }
  const double F = in.total(), T1 = in.T_K;
  // Isentropic ideal-gas path: Σ n_i [∫Cp_i/T dT] = F R ln(P2/P1); the mixing entropy does not change.
  const double target = F * 8.314462618 * std::log(P_out / in.P_Pa);
  auto ds = [&](double T) {
    double s = 0.0;
    for (std::size_t i = 0; i < c.size(); ++i) s += in.n[i] * (cp_over_t(c.pure[i], T) - cp_over_t(c.pure[i], T1));
    return s - target;
  };
  double lo = 10.0, hi = 6000.0;
  if (ds(lo) > 0.0 || ds(hi) < 0.0) {
    diag.error(DiagCode::VleNotConverged, "no isentropic outlet temperature between 10 and 6000 K");
    return std::nullopt;
  }
  for (int k = 0; k < 200 && hi - lo > 1e-10 * hi; ++k) {
    const double m = 0.5 * (lo + hi);
    (ds(m) > 0.0 ? hi : lo) = m;
  }
  const double T2s = 0.5 * (lo + hi);
  WorkResult r;
  for (std::size_t i = 0; i < c.size(); ++i) {
    if (in.n[i] == 0.0) continue;
    auto h1 = ideal_gas_enthalpy(c.pure[i], T1, diag), h2 = ideal_gas_enthalpy(c.pure[i], T2s, diag);
    if (!h1 || !h2) return std::nullopt;
    r.isentropic_work_W += in.n[i] * (*h2 - *h1);
  }
  r.work_W = P_out >= in.P_Pa ? r.isentropic_work_W / efficiency : r.isentropic_work_W * efficiency;
  auto o = stream_ph(c, P_out, in.H_W + r.work_W, in.n, diag, T2s);
  if (!o) return std::nullopt;
  r.out = *o;
  return r;
}

std::optional<ReactorResult> conversion_reactor(const ProcessComponents& c, const Stream& in, const std::vector<double>& nu,
                                                std::size_t key, double conversion, double T_out, double P_out,
                                                Diagnostics& diag) {
  const std::size_t n = c.size();
  if (nu.size() != n || key >= n || !(nu[key] < 0.0) || !(conversion >= 0.0 && conversion <= 1.0)) {
    diag.error(DiagCode::InvalidMixtureInput,
               "a conversion reactor needs one coefficient per component, a key reactant (negative coefficient) and a "
               "conversion in [0, 1]");
    return std::nullopt;
  }
  // Element balance from the structures.
  Composition sum{};
  for (std::size_t i = 0; i < n; ++i) {
    if (nu[i] == 0.0) continue;
    ChemicalSystemState st;
    Diagnostics q;
    const std::int32_t sp = parse_smiles(st, c.smiles[i], q);
    if (sp < 0) { merge(q, diag); return std::nullopt; }
    const Composition comp = compute_composition(st, static_cast<std::uint32_t>(sp));
    for (std::size_t e = 0; e < sum.size(); ++e) sum[e] += static_cast<std::int32_t>(std::lround(nu[i] * 1000.0)) * comp[e];
    if (!std::isfinite(c.pure[i].Hf_gas_Jmol)) {
      diag.error(DiagCode::NoPropertyData, "no enthalpy of formation for " + c.pure[i].name + ", which takes part in the reaction");
      return std::nullopt;
    }
  }
  for (std::int32_t e : sum)
    if (e != 0) {
      diag.error(DiagCode::UnbalancedReaction, "the stoichiometry does not conserve the elements");
      return std::nullopt;
    }
  ReactorResult r;
  r.extent_mols = conversion * in.n[key] / -nu[key];
  std::vector<double> out(n);
  const double F = in.total();
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = in.n[i] + nu[i] * r.extent_mols;
    if (out[i] < -1e-12 * F) {
      diag.error(DiagCode::InvalidMixtureInput, c.pure[i].name + " would be used up before the key component reaches the conversion");
      return std::nullopt;
    }
    out[i] = std::max(out[i], 0.0);
    r.heat_of_reaction_Jmol += nu[i] * (nu[i] == 0.0 ? 0.0 : c.pure[i].Hf_gas_Jmol);
  }
  const double P = P_out > 0.0 ? P_out : in.P_Pa;
  auto o = T_out > 0.0 ? stream_tp(c, T_out, P, out, diag) : stream_ph(c, P, in.H_W, out, diag, in.T_K);
  if (!o) return std::nullopt;
  r.out = *o;
  r.duty_W = o->H_W - in.H_W;
  return r;
}

}  // namespace chemlarp
