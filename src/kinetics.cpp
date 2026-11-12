#include "chemlarp/kinetics.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <optional>
#include <utility>
#include <vector>

#include "chemlarp/elements.hpp"

namespace chemlarp {

double ArrheniusRate::k(double T) const noexcept {
  return A * std::pow(T, b) * std::exp(-Ea_J / (kGasConstant * T));
}

std::int32_t Mechanism::find(std::string_view name) const {
  for (std::size_t i = 0; i < names.size(); ++i) {
    if (names[i].size() != name.size()) continue;
    bool same = true;
    for (std::size_t c = 0; c < name.size() && same; ++c)
      same = std::toupper(static_cast<unsigned char>(names[i][c])) == std::toupper(static_cast<unsigned char>(name[c]));
    if (same) return static_cast<std::int32_t>(i);
  }
  return -1;
}

namespace {

double int_pow(double c, double nu, bool& integral) {
  const double r = std::round(nu);
  integral = std::fabs(nu - r) < 1e-12 && r >= 0;
  if (integral) {
    double p = 1.0;
    for (int i = 0; i < static_cast<int>(r); ++i) p *= c;
    return p;
  }
  return std::pow(std::max(c, 0.0), nu);
}

double power(double c, double nu) {
  bool integral;
  return int_pow(c, nu, integral);
}

// d/dc c^nu
double dpower(double c, double nu) {
  bool integral;
  int_pow(c, nu, integral);
  if (integral) return nu == 0.0 ? 0.0 : nu * power(c, nu - 1.0);
  return c > 0.0 ? nu * std::pow(c, nu - 1.0) : 0.0;
}

constexpr double kLn10 = 2.302585092994046;

// Troe broadening factor F(T, Pr) and dF/dPr.
double troe_F(const TroeParameters& t, double T, double Pr, double* dF_dPr) {
  const double Fcent = std::max(t.F_cent(T), 1e-300);
  const double lFc = std::log10(Fcent);
  const double c = -0.4 - 0.67 * lFc;
  const double n = 0.75 - 1.27 * lFc;
  const double d = 0.14;
  const double x = std::log10(std::max(Pr, 1e-300));
  const double u = x + c;
  const double den = n - d * u;
  const double f = u / den;
  const double L = lFc / (1.0 + f * f);
  const double F = std::pow(10.0, L);
  if (dF_dPr) {
    if (Pr <= 1e-300) {
      *dF_dPr = 0.0;
    } else {
      const double dL_df = -2.0 * f * lFc / ((1.0 + f * f) * (1.0 + f * f));
      const double df_du = n / (den * den);
      *dF_dPr = F * kLn10 * dL_df * df_du / (Pr * kLn10);  // du/dPr = 1/(Pr ln10), dF/dL = F ln10
    }
  }
  return F;
}

}  // namespace

double TroeParameters::F_cent(double T) const noexcept {
  auto expo = [](double x) { return std::isfinite(x) ? std::exp(x) : 0.0; };
  double f = 0.0;
  if (std::fabs(T3) > 1e-300) f += (1.0 - a) * expo(-T / T3);
  if (std::fabs(T1) > 1e-300) f += a * expo(-T / T1);
  if (has_T2) f += expo(-T2 / T);
  return f;
}

double falloff_rate_constant(const KineticReaction& r, double T, double M, double* dk_dM) {
  const bool activated = r.kind == RateKind::ChemicallyActivated;
  const double k0 = activated ? r.forward.k(T) : r.low_or_high.k(T);
  const double kinf = activated ? r.low_or_high.k(T) : r.forward.k(T);
  if (!(kinf > 0.0)) {  // no high-pressure limit: pure low-pressure behaviour
    if (dk_dM) *dk_dM = activated ? 0.0 : k0;
    return activated ? k0 : k0 * M;
  }
  const double Pr = k0 * std::max(M, 0.0) / kinf;
  double dF = 0.0;
  const double F = r.troe ? troe_F(*r.troe, T, Pr, &dF) : 1.0;
  double k, dk_dPr;
  if (activated) {
    k = k0 / (1.0 + Pr) * F;
    dk_dPr = k0 * (-F / ((1.0 + Pr) * (1.0 + Pr)) + dF / (1.0 + Pr));
  } else {
    k = kinf * Pr / (1.0 + Pr) * F;
    dk_dPr = kinf * (F / ((1.0 + Pr) * (1.0 + Pr)) + Pr / (1.0 + Pr) * dF);
  }
  if (dk_dM) *dk_dM = dk_dPr * k0 / kinf;
  return k;
}

double plog_rate_constant(const KineticReaction& r, double T, double P, double* dk_dP) {
  if (dk_dP) *dk_dP = 0.0;
  if (r.plog.empty()) return 0.0;
  auto level_k = [&](const PlogLevel& l) {
    double k = 0.0;
    for (const auto& a : l.rates) k += a.k(T);
    return std::max(k, 1e-300);
  };
  if (P <= r.plog.front().P_Pa) return level_k(r.plog.front());
  if (P >= r.plog.back().P_Pa) return level_k(r.plog.back());
  std::size_t i = 0;
  while (i + 1 < r.plog.size() && r.plog[i + 1].P_Pa <= P) ++i;
  const double k1 = level_k(r.plog[i]), k2 = level_k(r.plog[i + 1]);
  const double slope = (std::log(k2) - std::log(k1)) / (std::log(r.plog[i + 1].P_Pa) - std::log(r.plog[i].P_Pa));
  const double k = std::exp(std::log(k1) + slope * (std::log(P) - std::log(r.plog[i].P_Pa)));
  if (dk_dP) *dk_dP = k * slope / P;
  return k;
}

// ---- rates -------------------------------------------------------------------

namespace {

double log_kc(const Mechanism& m, const KineticReaction& r, double T) {
  double dg = 0.0, dnu = 0.0;
  for (auto& [s, nu] : r.products) { dg += nu * m.thermo[s]->g_RT(T); dnu += nu; }
  for (auto& [s, nu] : r.reactants) { dg -= nu * m.thermo[s]->g_RT(T); dnu -= nu; }
  // K_c = K_p (P°/RT)^Δν
  return -dg + dnu * std::log(kStandardPressure / (kGasConstant * T));
}

// q = M (kf fwd - kr rev). For fall-off and PLOG, kr = ratio * kf with ratio = 1/K_c,
// so their pressure dependence enters through dkf_dM / dkf_dP only.
struct RateTerms {
  double kf = 0.0, kr = 0.0, M = 1.0, fwd = 1.0, rev = 1.0;
  double ratio = 0.0;  // kr / kf under detailed balance (0 if irreversible or explicit REV)
  double dkf_dM = 0.0, dkf_dP = 0.0;
};

double collider(const Mechanism& m, const KineticReaction& r, const double* C) {
  double M = 0.0;
  for (std::size_t i = 0; i < m.species_count(); ++i) M += r.efficiency[i] * C[i];
  return M;
}

RateTerms terms(const Mechanism& m, const KineticReaction& r, double T, const double* C, double P) {
  RateTerms t;
  switch (r.kind) {
    case RateKind::Elementary:
      t.kf = r.forward.k(T);
      break;
    case RateKind::ThirdBody:
      t.kf = r.forward.k(T);
      t.M = collider(m, r, C);
      break;
    case RateKind::FallOff:
    case RateKind::ChemicallyActivated:
      t.kf = falloff_rate_constant(r, T, collider(m, r, C), &t.dkf_dM);
      break;
    case RateKind::PressureLog:
      t.kf = plog_rate_constant(r, T, P, &t.dkf_dP);
      break;
  }
  if (r.reversible) {
    if (r.reverse) {
      t.kr = r.reverse->k(T);
    } else {
      t.ratio = std::exp(std::min(-log_kc(m, r, T), 700.0));
      t.kr = t.kf * t.ratio;
    }
  }
  for (auto& [s, nu] : r.reactants) t.fwd *= power(C[s], nu);
  if (r.reversible)
    for (auto& [s, nu] : r.products) t.rev *= power(C[s], nu);
  else
    t.rev = 0.0;
  return t;
}

double total_pressure(const Mechanism& m, double T, const double* C) {
  double c = 0.0;
  for (std::size_t i = 0; i < m.species_count(); ++i) c += C[i];
  return c * kGasConstant * T;
}

}  // namespace

double equilibrium_constant_c(const Mechanism& m, std::size_t r, double T) {
  return std::exp(log_kc(m, m.reactions[r], T));
}

void rates_of_progress(const Mechanism& m, double T, const double* C, double* q) {
  const double P = total_pressure(m, T, C);
  for (std::size_t r = 0; r < m.reactions.size(); ++r) {
    const RateTerms t = terms(m, m.reactions[r], T, C, P);
    q[r] = t.M * (t.kf * t.fwd - t.kr * t.rev);
  }
}

void production_rates(const Mechanism& m, double T, const double* C, double* wdot) {
  std::fill(wdot, wdot + m.species_count(), 0.0);
  const double P = total_pressure(m, T, C);
  for (const auto& r : m.reactions) {
    const RateTerms t = terms(m, r, T, C, P);
    const double q = t.M * (t.kf * t.fwd - t.kr * t.rev);
    for (auto& [s, nu] : r.reactants) wdot[s] -= nu * q;
    for (auto& [s, nu] : r.products) wdot[s] += nu * q;
  }
}

void production_rates_jacobian(const Mechanism& m, double T, const double* C, double* J) {
  const std::size_t S = m.species_count();
  std::fill(J, J + S * S, 0.0);
  std::vector<double> dq(S);
  const double P = total_pressure(m, T, C);
  for (const auto& r : m.reactions) {
    const RateTerms t = terms(m, r, T, C, P);
    std::fill(dq.begin(), dq.end(), 0.0);
    // Mass-action part, product rule over distinct species.
    for (auto& [j, nuj] : r.reactants) {
      double d = dpower(C[j], nuj);
      for (auto& [k, nuk] : r.reactants)
        if (k != j) d *= power(C[k], nuk);
      dq[j] += t.M * t.kf * d;
    }
    if (r.reversible) {
      for (auto& [j, nuj] : r.products) {
        double d = dpower(C[j], nuj);
        for (auto& [k, nuk] : r.products)
          if (k != j) d *= power(C[k], nuk);
        dq[j] -= t.M * t.kr * d;
      }
    }
    if (r.kind == RateKind::ThirdBody) {
      const double base = t.kf * t.fwd - t.kr * t.rev;
      for (std::size_t j = 0; j < S; ++j) dq[j] += r.efficiency[j] * base;
    } else if (r.kind == RateKind::FallOff || r.kind == RateKind::ChemicallyActivated) {
      const double base = t.dkf_dM * (t.fwd - t.ratio * t.rev);
      for (std::size_t j = 0; j < S; ++j) dq[j] += r.efficiency[j] * base;
    } else if (r.kind == RateKind::PressureLog) {
      const double base = t.dkf_dP * kGasConstant * T * (t.fwd - t.ratio * t.rev);  // dP/dC_j = RT
      for (std::size_t j = 0; j < S; ++j) dq[j] += base;
    }
    for (auto& [i, nu] : r.reactants)
      for (std::size_t j = 0; j < S; ++j) J[i * S + j] -= nu * dq[j];
    for (auto& [i, nu] : r.products)
      for (std::size_t j = 0; j < S; ++j) J[i * S + j] += nu * dq[j];
  }
}

// ---- reactors ----------------------------------------------------------------

const char* to_string(ReactorModel model) noexcept {
  switch (model) {
    case ReactorModel::IsothermalConstantVolume: return "iso-v";
    case ReactorModel::IsothermalConstantPressure: return "iso-p";
    case ReactorModel::AdiabaticConstantVolume: return "adiabatic-v";
    case ReactorModel::AdiabaticConstantPressure: return "adiabatic-p";
  }
  return "?";
}

std::optional<ReactorModel> reactor_model_from_string(std::string_view s) noexcept {
  if (s == "iso-v") return ReactorModel::IsothermalConstantVolume;
  if (s == "iso-p") return ReactorModel::IsothermalConstantPressure;
  if (s == "adiabatic-v") return ReactorModel::AdiabaticConstantVolume;
  if (s == "adiabatic-p") return ReactorModel::AdiabaticConstantPressure;
  return std::nullopt;
}

ReactorResult run_reactor(const Mechanism& m, const ReactorSetup& su, Diagnostics& diag) {
  ReactorResult res;
  const std::size_t S = m.species_count();
  const bool adiabatic =
      su.model == ReactorModel::AdiabaticConstantVolume || su.model == ReactorModel::AdiabaticConstantPressure;
  const bool const_v =
      su.model == ReactorModel::IsothermalConstantVolume || su.model == ReactorModel::AdiabaticConstantVolume;

  bool ok = true;
  if (su.moles.size() != S) {
    diag.error(DiagCode::InvalidReactorSetup, "initial moles must have one entry per mechanism species (" +
                                                  std::to_string(S) + ")");
    return res;
  }
  double N0 = 0.0;
  for (std::size_t i = 0; i < S; ++i) {
    if (!(su.moles[i] >= 0.0) || !std::isfinite(su.moles[i])) {
      diag.error(DiagCode::NegativeAmount, "initial amount of " + m.names[i] + " must be finite and >= 0").on_species(i);
      ok = false;
    }
    N0 += su.moles[i];
  }
  if (!(su.temperature_K > 0) || !(su.pressure_Pa > 0) || !(su.t_end > 0) || !(N0 > 0)) {
    diag.error(DiagCode::InvalidReactorSetup, "temperature, pressure, t_end and total moles must all be positive");
    ok = false;
  }
  if (adiabatic) {
    for (std::size_t i = 0; i < S; ++i) {
      if (!m.thermo[i] || !std::isfinite(m.thermo[i]->cp_R(su.temperature_K))) {
        diag.error(DiagCode::MissingHeatCapacity,
                   m.names[i] + " has no heat-capacity data; adiabatic reactors need NASA-7 thermo for every species")
            .on_species(i);
        ok = false;
      }
    }
  }
  if (!ok) return res;

  const double V0 = N0 * kGasConstant * su.temperature_K / su.pressure_Pa;
  const std::size_t n = S + (adiabatic ? 1 : 0);
  auto volume = [&](const std::vector<double>& y, double T) {
    if (const_v) return V0;
    double N = 0.0;
    for (std::size_t i = 0; i < S; ++i) N += y[i];
    return std::max(N, 1e-300) * kGasConstant * T / su.pressure_Pa;
  };
  std::vector<double> C(S), wdot(S);
  OdeRhs rhs = [&](double, const std::vector<double>& y, std::vector<double>& dy) {
    const double T = adiabatic ? y[S] : su.temperature_K;
    const double V = volume(y, T);
    for (std::size_t i = 0; i < S; ++i) C[i] = y[i] / V;
    production_rates(m, T, C.data(), wdot.data());
    for (std::size_t i = 0; i < S; ++i) dy[i] = V * wdot[i];
    if (adiabatic) {
      double num = 0.0, den = 0.0;
      for (std::size_t i = 0; i < S; ++i) {
        const auto& th = *m.thermo[i];
        if (const_v) {
          num += (th.h_RT(T) - 1.0) * kGasConstant * T * dy[i];    // u_i dn_i
          den += y[i] * (th.cp_R(T) - 1.0) * kGasConstant;          // n_i cv_i
        } else {
          num += th.h_RT(T) * kGasConstant * T * dy[i];             // h_i dn_i
          den += y[i] * th.cp_R(T) * kGasConstant;                  // n_i cp_i
        }
      }
      dy[S] = -num / den;
    }
  };
  std::vector<double> Jc(S * S), f0(n), f1(n), yT;
  OdeJacobian jac = [&](double t, const std::vector<double>& y, std::vector<double>& J) {
    // Species block: exact for constant volume (dn/dt = V w(n/V) => d/dn = dw/dC),
    // a good modified-Newton approximation at constant pressure.
    const double T = adiabatic ? y[S] : su.temperature_K;
    const double V = volume(y, T);
    for (std::size_t i = 0; i < S; ++i) C[i] = y[i] / V;
    production_rates_jacobian(m, T, C.data(), Jc.data());
    J.assign(n * n, 0.0);
    for (std::size_t i = 0; i < S; ++i)
      for (std::size_t j = 0; j < S; ++j) J[i * n + j] = Jc[i * S + j];
    if (!adiabatic) return;
    // Temperature column: one well-scaled forward difference (captures k(T), K_c(T), V(T)).
    rhs(t, y, f0);
    yT = y;
    const double dT = 1e-6 * std::max(T, 1.0);
    yT[S] += dT;
    rhs(t, yT, f1);
    for (std::size_t i = 0; i < n; ++i) J[i * n + S] = (f1[i] - f0[i]) / dT;
    // Temperature row, analytic: Tdot = -num/den, num = sum e_i ndot_i, den = sum n_i c_i,
    // with (e, c) = (u, cv) at constant volume and (h, cp) at constant pressure.
    double num = 0.0, den = 0.0;
    std::vector<double> e(S), cap(S);
    for (std::size_t i = 0; i < S; ++i) {
      const auto& th = *m.thermo[i];
      e[i] = (th.h_RT(T) - (const_v ? 1.0 : 0.0)) * kGasConstant * T;
      cap[i] = (th.cp_R(T) - (const_v ? 1.0 : 0.0)) * kGasConstant;
      num += e[i] * f0[i];
      den += y[i] * cap[i];
    }
    for (std::size_t j = 0; j < S; ++j) {
      double dnum = 0.0;
      for (std::size_t i = 0; i < S; ++i) dnum += e[i] * J[i * n + j];
      J[S * n + j] = -dnum / den + num * cap[j] / (den * den);
    }
  };

  auto snapshot = [&](double t, const std::vector<double>& y) {
    ReactorState s;
    s.t = t;
    s.T = adiabatic ? y[S] : su.temperature_K;
    s.V = volume(y, s.T);
    s.n.assign(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(S));
    double N = 0.0;
    for (double v : s.n) N += v;
    s.P = const_v ? N * kGasConstant * s.T / V0 : su.pressure_Pa;
    return s;
  };

  std::vector<double> y(su.moles);
  if (adiabatic) y.push_back(su.temperature_K);
  BdfOptions opt;
  opt.rtol = su.rtol;
  opt.atol.assign(n, su.atol * N0);
  if (adiabatic) opt.atol[S] = 1e-6;
  opt.max_steps = su.max_steps;

  res.samples.push_back(snapshot(0.0, y));
  double t_prev = 0.0, T_prev = su.temperature_K;
  OdeObserver on_step = [&](double t, const std::vector<double>& yy) {
    if (adiabatic) {
      const double rate = (yy[S] - T_prev) / (t - t_prev);
      if (rate > res.max_dTdt) {
        res.max_dTdt = rate;
        res.ignition_time = 0.5 * (t + t_prev);
      }
      T_prev = yy[S];
    }
    t_prev = t;
    if (su.output_times.empty()) res.samples.push_back(snapshot(t, yy));
  };
  OdeObserver on_output = [&](double t, const std::vector<double>& yy) { res.samples.push_back(snapshot(t, yy)); };

  res.ok = integrate_bdf(rhs, &jac, 0.0, y, su.t_end, su.output_times, on_output, &on_step, opt,
                         res.stats, diag);
  res.final_state = snapshot(su.t_end, y);
  if (adiabatic && res.ok) {
    double tmin = su.temperature_K, tmax = su.temperature_K;
    for (const auto& s : res.samples) { tmin = std::min(tmin, s.T); tmax = std::max(tmax, s.T); }
    tmax = std::max(tmax, res.final_state.T);
    for (std::size_t i = 0; i < S; ++i) {
      if (!m.thermo[i]->in_range(tmin) || !m.thermo[i]->in_range(tmax)) {
        diag.warn(DiagCode::TemperatureOutOfRange, m.names[i] + ": temperature left the thermo data range (" +
                                                       std::to_string(tmin) + "-" + std::to_string(tmax) + " K)")
            .on_species(i);
      }
    }
  }
  if (res.max_dTdt <= 0.0) res.ignition_time = -1.0;
  return res;
}

}  // namespace chemlarp
