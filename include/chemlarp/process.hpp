#pragma once
// Steady-state material streams and unit operations.
//
// Thermodynamic model (an "ideal"-type package with activity coefficients), per component i:
//   vapour   H_V,i(T) = ΔfH°_i + ∫_298.15^T Cp_ig,i dT          (ideal gas; ΔfH° taken as 0 when not tabulated,
//                                                                which cancels in every non-reacting balance)
//   liquid   H_L,i(T) = H_V,i(Tb) − ΔHvap,i(Tb) + ∫_Tb^T Cp_L,i dT   (measured ΔHvap at the normal boiling
//            point and liquid Cp; no excess enthalpy and no pressure effect). Without a liquid Cp the
//            liquid is H_V,i(T) − ΔHvap,i(T) (Watson) and a warning says so; a component above its Tc in
//            the liquid is given its vapour enthalpy, also with a warning.
//   phases   modified Raoult's law y_i P = x_i γ_i P_i^sat with γ from UNIFAC or γ = 1 (ideal solution),
//            Antoine vapour pressures (activity.hpp); one liquid phase only.
// Flows are in mol/s, enthalpy flows in W (J/s), heat duties in W (positive = heat added).

#include <optional>
#include <string>
#include <vector>

#include "chemlarp/activity.hpp"
#include "chemlarp/diagnostics.hpp"
#include "chemlarp/properties.hpp"

namespace chemlarp {

enum class LiquidModel { Ideal, Unifac };

struct ProcessComponents {
  std::vector<std::string> smiles;
  std::vector<PureComponent> pure;
  ActivitySystem vle;  // vapour pressures and the activity model
  std::size_t size() const noexcept { return pure.size(); }
};

std::optional<ProcessComponents> process_components(const std::vector<std::string>& smiles, LiquidModel model,
                                                    Diagnostics& diag);

struct Stream {
  double T_K = 0.0, P_Pa = 0.0;
  std::vector<double> n;               // component flows, mol/s
  double vapor_fraction = 0.0;         // molar
  std::vector<double> x, y;            // liquid and vapour mole fractions (equal to z for a single phase)
  double H_W = 0.0;                    // enthalpy flow
  double total() const noexcept;
};

// Phase split and enthalpy at given T and P (TP flash), and at given P and enthalpy flow (PH flash).
std::optional<Stream> stream_tp(const ProcessComponents& c, double T, double P, const std::vector<double>& n,
                                Diagnostics& diag);
std::optional<Stream> stream_ph(const ProcessComponents& c, double P, double H_W, const std::vector<double>& n,
                                Diagnostics& diag, double T_guess = 0.0);

// Unit operations. Every unit conserves each component flow and closes its energy balance exactly
// (inlet enthalpy + duty = outlet enthalpy) up to the flash tolerance.
// Adiabatic mixer; outlet pressure P_out, or the lowest inlet pressure when P_out <= 0.
std::optional<Stream> mixer(const ProcessComponents& c, const std::vector<Stream>& inlets, double P_out, Diagnostics& diag);
// Splitter into outlets with the given fractions of the inlet (summing to 1); same T, P and phase split.
std::optional<std::vector<Stream>> splitter(const ProcessComponents& c, const Stream& in, const std::vector<double>& fractions,
                                            Diagnostics& diag);

struct HeaterResult {
  Stream out;
  double duty_W = 0.0;
};
// Heater or cooler to a given outlet temperature, or with a given duty; P_out <= 0 keeps the inlet pressure.
std::optional<HeaterResult> heater_to_temperature(const ProcessComponents& c, const Stream& in, double T_out, double P_out,
                                                  Diagnostics& diag);
std::optional<HeaterResult> heater_with_duty(const ProcessComponents& c, const Stream& in, double duty_W, double P_out,
                                             Diagnostics& diag);
// Valve: isenthalpic expansion to P_out.
std::optional<Stream> valve(const ProcessComponents& c, const Stream& in, double P_out, Diagnostics& diag);

struct FlashDrumResult {
  Stream vapor, liquid;
  double duty_W = 0.0;
};
// Flash drum at given T and P (duty reported), or adiabatic at P (duty zero).
std::optional<FlashDrumResult> flash_drum(const ProcessComponents& c, const Stream& in, double T, double P, Diagnostics& diag);
std::optional<FlashDrumResult> flash_drum_adiabatic(const ProcessComponents& c, const Stream& in, double P, Diagnostics& diag);

}  // namespace chemlarp
