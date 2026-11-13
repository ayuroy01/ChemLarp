#include "chemlarp/properties.hpp"

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <map>
#include <string>

#include "chemlarp/smiles.hpp"
#include "chemlarp/system_state.hpp"
#include "chemlarp/thermo.hpp"

namespace chemlarp {
namespace {

struct PropertyRow {
  const char *smiles, *name, *cas;
  double mw, Tc, Pc, Vc;
  char crit_source;
  double omega, Tm, Tb, Hvap_Tb, a0, a1, a2, a3, a4, cp_Tmin, cp_Tmax, cp_liquid, Z_RA, Hf, Gf, S;
};

#include "property_data.inc"

struct TransportRow {
  const char* smiles;
  char property, equation;
  double c[7];
  double Tmin, Tmax;
};

#include "transport_data.inc"

std::string canonical(std::string_view smiles) {
  ChemicalSystemState st;
  Diagnostics d;
  const std::int32_t sp = parse_smiles(st, smiles, d);
  if (sp < 0 || d.has_errors()) return {};
  return st.species[static_cast<std::size_t>(sp)].canonical_smiles;
}

const std::map<std::string, PureComponent>& table() {
  static const std::map<std::string, PureComponent> t = [] {
    std::map<std::string, PureComponent> m;
    for (const PropertyRow& r : kProperties) {
      const std::string key = canonical(r.smiles);
      if (key.empty() || m.count(key)) continue;
      PureComponent c;
      c.smiles = key;
      c.name = r.name;
      c.cas = r.cas;
      c.molar_mass = r.mw;
      c.Tc_K = r.Tc;
      c.Pc_Pa = r.Pc;
      c.Vc_m3mol = r.Vc;
      c.critical_source = r.crit_source;
      c.omega = r.omega;
      c.Tm_K = r.Tm;
      c.Tb_K = r.Tb;
      c.Hvap_Tb_Jmol = r.Hvap_Tb;
      c.cp_a[0] = r.a0; c.cp_a[1] = r.a1; c.cp_a[2] = r.a2; c.cp_a[3] = r.a3; c.cp_a[4] = r.a4;
      c.cp_Tmin_K = r.cp_Tmin;
      c.cp_Tmax_K = r.cp_Tmax;
      c.cp_liquid_298_JmolK = r.cp_liquid;
      c.Z_RA = r.Z_RA;
      c.Hf_gas_Jmol = r.Hf;
      c.Gf_gas_Jmol = r.Gf;
      c.S_gas_JmolK = r.S;
      m.emplace(key, std::move(c));
    }
    for (const TransportRow& r : kTransport) {
      auto it = m.find(canonical(r.smiles));
      if (it == m.end()) continue;
      PureComponent& c = it->second;
      TransportCorrelation* corr = r.property == 'L' ? &c.liquid_viscosity
                                : r.property == 'G' ? &c.gas_viscosity
                                : r.property == 'l' ? &c.liquid_conductivity
                                : r.property == 'g' ? &c.gas_conductivity
                                : r.property == 's' ? &c.surface_tension
                                : r.property == 'C' ? &c.liquid_heat_capacity
                                : r.property == 'V' ? &c.liquid_density : nullptr;
      if (!corr || corr->present()) continue;
      corr->equation = r.equation;
      for (int i = 0; i < 7; ++i) corr->c[i] = r.c[i];
      corr->Tmin_K = r.Tmin;
      corr->Tmax_K = r.Tmax;
    }
    return m;
  }();
  return t;
}

bool known(double v) { return std::isfinite(v); }

std::string num(double v) {
  char b[32];
  std::snprintf(b, sizeof b, "%g", v);
  return b;
}

bool valid_T(double T, Diagnostics& diag) {
  if (std::isfinite(T) && T > 0.0) return true;
  diag.error(DiagCode::InvalidMixtureInput, "temperature must be a positive number of kelvin");
  return false;
}

bool missing(const PureComponent& c, const char* what, Diagnostics& diag) {
  diag.error(DiagCode::NoPropertyData, std::string("no ") + what + " for " + c.name);
  return false;
}

bool has_cp(const PureComponent& c, double T, Diagnostics& diag) {
  if (!valid_T(T, diag)) return false;
  for (double a : c.cp_a)
    if (!known(a)) return missing(c, "ideal-gas heat-capacity coefficients", diag);
  if (known(c.cp_Tmin_K) && known(c.cp_Tmax_K) && (T < c.cp_Tmin_K || T > c.cp_Tmax_K))
    diag.warn(DiagCode::PropertyExtrapolated, "ideal-gas Cp of " + c.name + " extrapolated to " + num(T) + " K (fitted " +
                                                  num(c.cp_Tmin_K) + "-" + num(c.cp_Tmax_K) + " K)");
  return true;
}

// The liquid correlations need T below Tc; below the melting point they are extrapolations.
bool liquid_T(const PureComponent& c, double T, const char* what, Diagnostics& diag) {
  if (!valid_T(T, diag)) return false;
  if (!known(c.Tc_K)) return missing(c, "critical temperature", diag);
  if (T >= c.Tc_K) {
    diag.error(DiagCode::SupercriticalTemperature,
               std::string(what) + " of " + c.name + " asked at or above Tc = " + num(c.Tc_K) + " K");
    return false;
  }
  if (known(c.Tm_K) && T < c.Tm_K)
    diag.warn(DiagCode::PropertyExtrapolated, std::string(what) + " of " + c.name + " extrapolated below the melting point");
  return true;
}

constexpr double kT0 = kTref;

}  // namespace

std::optional<PureComponent> pure_component(std::string_view smiles) {
  const std::string key = canonical(smiles);
  if (key.empty()) return std::nullopt;
  auto it = table().find(key);
  if (it == table().end()) return std::nullopt;
  return it->second;
}

std::size_t pure_component_count() noexcept { return table().size(); }

std::optional<double> ideal_gas_cp(const PureComponent& c, double T, Diagnostics& diag) {
  if (!has_cp(c, T, diag)) return std::nullopt;
  const double* a = c.cp_a;
  return kGasConstant * (a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4]))));
}

std::optional<double> ideal_gas_enthalpy(const PureComponent& c, double T, Diagnostics& diag) {
  if (!has_cp(c, T, diag)) return std::nullopt;
  const double* a = c.cp_a;
  auto H = [&](double t) { return t * (a[0] + t * (a[1] / 2 + t * (a[2] / 3 + t * (a[3] / 4 + t * a[4] / 5)))); };
  return kGasConstant * (H(T) - H(kT0));
}

std::optional<double> ideal_gas_entropy(const PureComponent& c, double T, double P, Diagnostics& diag) {
  if (!(std::isfinite(P) && P > 0.0)) {
    diag.error(DiagCode::InvalidMixtureInput, "pressure must be a positive number of pascals");
    return std::nullopt;
  }
  if (!has_cp(c, T, diag)) return std::nullopt;
  if (!known(c.S_gas_JmolK)) {
    missing(c, "standard entropy", diag);
    return std::nullopt;
  }
  const double* a = c.cp_a;
  auto S = [&](double t) { return a[0] * std::log(t) + t * (a[1] + t * (a[2] / 2 + t * (a[3] / 3 + t * a[4] / 4))); };
  return c.S_gas_JmolK + kGasConstant * (S(T) - S(kT0) - std::log(P / kStandardPressure));
}

std::optional<double> enthalpy_of_vaporization(const PureComponent& c, double T, Diagnostics& diag) {
  if (!liquid_T(c, T, "enthalpy of vaporization", diag)) return std::nullopt;
  if (!known(c.Tb_K) || !known(c.Hvap_Tb_Jmol)) {
    missing(c, "enthalpy of vaporization at the normal boiling point", diag);
    return std::nullopt;
  }
  return c.Hvap_Tb_Jmol * std::pow((1.0 - T / c.Tc_K) / (1.0 - c.Tb_K / c.Tc_K), 0.38);
}

std::optional<double> saturated_liquid_volume(const PureComponent& c, double T, Diagnostics& diag) {
  if (!liquid_T(c, T, "saturated liquid volume", diag)) return std::nullopt;
  if (c.liquid_density.present()) {
    const double* k = c.liquid_density.c;
    const double tau = std::max(0.0, 1.0 - T / k[5]);
    const double t3 = std::cbrt(tau);
    const double rho = k[0] + k[1] * std::pow(tau, 0.35) + k[2] * t3 * t3 + tau * (k[3] + k[4] * t3);  // kg/m³
    if (rho > 0.0) return 1e-3 * k[6] / rho;
  }
  if (!known(c.Pc_Pa)) {
    missing(c, "critical pressure", diag);
    return std::nullopt;
  }
  double Z = c.Z_RA;
  if (!known(Z)) {
    if (!known(c.Vc_m3mol)) {
      missing(c, "Rackett parameter or critical volume", diag);
      return std::nullopt;
    }
    Z = c.Pc_Pa * c.Vc_m3mol / (kGasConstant * c.Tc_K);
    diag.warn(DiagCode::PropertyExtrapolated, "liquid volume of " + c.name +
                                                  " from the original Rackett equation with Zc (no fitted Z_RA)");
  }
  return kGasConstant * c.Tc_K / c.Pc_Pa * std::pow(Z, 1.0 + std::pow(1.0 - T / c.Tc_K, 2.0 / 7.0));
}

namespace {

// Evaluates a transport correlation; `liquid` properties are refused at or above Tc.
std::optional<double> transport(const PureComponent& c, const TransportCorrelation& t, double T, bool liquid,
                                const char* what, Diagnostics& diag) {
  if (!valid_T(T, diag)) return std::nullopt;
  if (!t.present()) {
    missing(c, what, diag);
    return std::nullopt;
  }
  const double Tc = t.equation == 'M' ? t.c[6] : c.Tc_K;
  if (liquid && known(Tc) && T >= Tc) {
    diag.error(DiagCode::SupercriticalTemperature,
               std::string(what) + " of " + c.name + " asked at or above Tc = " + num(Tc) + " K");
    return std::nullopt;
  }
  if (T < t.Tmin_K || T > t.Tmax_K)
    diag.warn(DiagCode::PropertyExtrapolated, std::string(what) + " of " + c.name + " extrapolated to " + num(T) +
                                                  " K (fitted " + num(t.Tmin_K) + "-" + num(t.Tmax_K) + " K)");
  const double* k = t.c;
  double v = 0.0;
  switch (t.equation) {
    case '0': v = k[0] + T * (k[1] + T * (k[2] + T * (k[3] + T * k[4]))); break;
    case '1': v = std::exp(k[0] + k[1] / T + k[2] * std::log(T) + k[3] * std::pow(T, k[4])); break;
    case '2': v = k[0] * std::pow(T, k[1]) / (1.0 + k[2] / T + k[3] / (T * T)); break;
    case 'M': {
      const double tau = 1.0 - T / k[6];
      for (int i = 0; i < 3; ++i)
        if (k[2 * i] != 0.0) v += k[2 * i] * std::pow(tau, k[2 * i + 1]);
      break;
    }
    case 'J': v = (k[0] - k[1] * (T - 273.15)) * 1e-3; break;
    default: break;
  }
  if (!std::isfinite(v) || v <= 0.0) {
    diag.error(DiagCode::NoPropertyData, std::string(what) + " of " + c.name + " is not defined at " + num(T) + " K");
    return std::nullopt;
  }
  return v;
}

}  // namespace

std::optional<double> liquid_viscosity(const PureComponent& c, double T, Diagnostics& diag) {
  return transport(c, c.liquid_viscosity, T, true, "liquid viscosity", diag);
}
std::optional<double> gas_viscosity(const PureComponent& c, double T, Diagnostics& diag) {
  return transport(c, c.gas_viscosity, T, false, "vapour viscosity", diag);
}
std::optional<double> liquid_thermal_conductivity(const PureComponent& c, double T, Diagnostics& diag) {
  return transport(c, c.liquid_conductivity, T, true, "liquid thermal conductivity", diag);
}
std::optional<double> gas_thermal_conductivity(const PureComponent& c, double T, Diagnostics& diag) {
  return transport(c, c.gas_conductivity, T, false, "vapour thermal conductivity", diag);
}
std::optional<double> surface_tension(const PureComponent& c, double T, Diagnostics& diag) {
  return transport(c, c.surface_tension, T, true, "surface tension", diag);
}

std::optional<double> liquid_cp(const PureComponent& c, double T, Diagnostics& diag) {
  return transport(c, c.liquid_heat_capacity, T, true, "liquid heat capacity", diag);
}

std::optional<double> liquid_enthalpy_change(const PureComponent& c, double T0, double T1, Diagnostics& diag) {
  // Both ends are checked (range, Tc) by evaluating Cp there.
  if (!liquid_cp(c, T0, diag) || !liquid_cp(c, T1, diag)) return std::nullopt;
  const double* k = c.liquid_heat_capacity.c;
  auto F = [&](double t) { return t * (k[0] + t * (k[1] / 2 + t * (k[2] / 3 + t * (k[3] / 4 + t * k[4] / 5)))); };
  return F(T1) - F(T0);
}

}  // namespace chemlarp
