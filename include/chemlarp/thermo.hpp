#pragma once
// Standard-state thermochemistry. Convention (NASA / JANAF / Cantera):
//   H°(T) referenced so that elements in their standard state have H = 0 at 298.15 K
//   (hence H°(298.15) = ΔfH°), S°(T) absolute (third-law), G° = H° - T S°.
// All species in one equilibrium problem must share this convention.

#include <array>
#include <optional>
#include <string>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/system_state.hpp"

namespace chemlarp {

inline constexpr double kGasConstant = 8.314462618;  // J/(mol K)
inline constexpr double kStandardPressure = 1.0e5;   // Pa
inline constexpr double kTref = 298.15;               // K

// NASA 7-coefficient polynomial (two temperature ranges):
//   Cp/R = a0 + a1 T + a2 T^2 + a3 T^3 + a4 T^4
//   H/RT = a0 + a1 T/2 + a2 T^2/3 + a3 T^3/4 + a4 T^4/5 + a5/T
//   S/R  = a0 ln T + a1 T + a2 T^2/2 + a3 T^3/3 + a4 T^4/4 + a6
struct Nasa7 {
  double t_low = 0, t_mid = 0, t_high = 0;
  std::array<double, 7> low{}, high{};

  const std::array<double, 7>& coeffs(double T) const noexcept { return T < t_mid ? low : high; }
  double cp_R(double T) const noexcept;
  double h_RT(double T) const noexcept;
  double s_R(double T) const noexcept;
  double g_RT(double T) const noexcept { return h_RT(T) - s_R(T); }
};

enum class ThermoSource { Nasa7Database, JobackEstimate };

struct SpeciesThermo {
  ThermoSource source = ThermoSource::Nasa7Database;
  Nasa7 nasa{};
  double h298_J = 0.0;  // Joback: ΔfH°(298.15) in J/mol (held constant)
  double s298_J = 0.0;  // Joback: absolute S°(298.15) in J/(mol K) (held constant)
  std::string label;

  double h_RT(double T) const noexcept;
  double s_R(double T) const noexcept;
  double g_RT(double T) const noexcept { return h_RT(T) - s_R(T); }
  double cp_R(double T) const noexcept;  // NaN for Joback estimates (no heat-capacity model)
  bool in_range(double T) const noexcept;
};

// Built-in NASA-7 data (GRI-Mech 3.0 thermo set), keyed by canonical SMILES.
std::optional<SpeciesThermo> builtin_thermo_canonical(const std::string& canonical_smiles);
// Formula lookup (first entry with that Hill formula); prefer builtin_thermo_canonical.
std::optional<SpeciesThermo> builtin_thermo(const std::string& hill_formula);

// Database lookup first; otherwise Joback group-additivity estimate from the graph.
// Appends NoThermoData / NoJobackGroup errors or ThermoApproximation / TemperatureOutOfRange warnings.
std::optional<SpeciesThermo> resolve_thermo(const ChemicalSystemState& state, std::uint32_t species,
                                            double temperature_K, Diagnostics& diag);

// Fill state.g_RT for every species at state.conditions. Returns false on any missing data.
bool bind_thermo(ChemicalSystemState& state, Diagnostics& diag);

}  // namespace chemlarp
