#pragma once
// Gas-phase chemical kinetics: mechanisms, mass-action rates, batch reactors.
//
// Mechanism text format (.ocm). Species are bound to structures by SMILES, so
// element balance, thermo and identity come from the molecular graph:
//
//   species
//     H2   [H][H]
//     O2   O=O
//     ...
//   end
//   reactions  cal/mol  cm3          # Ea unit: cal/mol|kcal/mol|J/mol|kJ/mol|K ; volume unit: cm3|m3
//     H + O2 <=> O + OH      3.547e15  -0.406  16599     # A  b  Ea   (k = A T^b exp(-Ea/RT))
//     H2 + M <=> H + H + M   4.577e19  -1.40   104380
//       H2/2.5/ H2O/12.0/                                 # third-body efficiencies (CHEMKIN style)
//     2 B => B + C           3.0e7      0       0
//   end
//
// Arrows: "<=>" or "=" reversible, "=>" irreversible. Reverse rate constants
// come from detailed balance with the species thermo, k_r = k_f / K_c, so a
// reversible mechanism relaxes exactly to the Gibbs equilibrium of gibbs.hpp.
// Not supported yet (reported as UnsupportedReactionType): pressure fall-off
// "(+M)", PLOG, explicit REV parameters, non-elementary orders (FORD/RORD).

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "omnichem/diagnostics.hpp"
#include "omnichem/ode.hpp"
#include "omnichem/system_state.hpp"
#include "omnichem/thermo.hpp"

namespace omnichem {

struct ArrheniusRate {
  double A = 0.0;     // SI: (m^3/mol)^(order-1) / s
  double b = 0.0;     // temperature exponent
  double Ea_J = 0.0;  // activation energy, J/mol
  double k(double T) const noexcept;
};

struct KineticReaction {
  std::vector<std::pair<std::uint32_t, double>> reactants;  // (species, stoichiometric coefficient)
  std::vector<std::pair<std::uint32_t, double>> products;
  ArrheniusRate forward;
  bool reversible = false;
  bool third_body = false;
  std::vector<double> efficiency;  // per species, used when third_body (default 1)
  std::string equation;
  std::size_t line = 0;
};

struct Mechanism {
  ChemicalSystemState state;  // one species per mechanism species, in declaration order
  std::vector<std::string> names;
  std::vector<KineticReaction> reactions;
  std::vector<std::optional<SpeciesThermo>> thermo;  // per species (nullopt if unavailable)

  std::size_t species_count() const noexcept { return names.size(); }
  std::int32_t find(std::string_view name) const;  // case-insensitive, -1 if absent
};

std::optional<Mechanism> parse_mechanism(std::string_view text, Diagnostics& diag);
std::optional<Mechanism> load_mechanism(const std::string& path, Diagnostics& diag);

// C: concentrations, mol/m^3 (one per species). Outputs in mol/(m^3 s).
void rates_of_progress(const Mechanism& m, double T, const double* C, double* q);
void production_rates(const Mechanism& m, double T, const double* C, double* wdot);
// Row-major S x S, J[i*S + j] = d wdot_i / d C_j (exact, analytic).
void production_rates_jacobian(const Mechanism& m, double T, const double* C, double* J);
// K_c of reaction r in (mol/m^3)^(Δν).
double equilibrium_constant_c(const Mechanism& m, std::size_t r, double T);

enum class ReactorModel {
  IsothermalConstantVolume,
  IsothermalConstantPressure,
  AdiabaticConstantVolume,
  AdiabaticConstantPressure,
};

struct ReactorSetup {
  ReactorModel model = ReactorModel::IsothermalConstantVolume;
  double temperature_K = 1000.0;
  double pressure_Pa = 1.0e5;      // initial pressure; fixes V for constant-volume models
  std::vector<double> moles;       // initial moles per species
  double t_end = 1.0;              // s
  std::vector<double> output_times;  // empty: record every accepted step
  double rtol = 1e-8;
  double atol = 1e-14;             // absolute tolerance on moles, relative to total initial moles
  std::size_t max_steps = 500000;
};

struct ReactorState {
  double t = 0.0, T = 0.0, P = 0.0, V = 0.0;
  std::vector<double> n;  // moles
};

struct ReactorResult {
  bool ok = false;
  std::vector<ReactorState> samples;
  ReactorState final_state;
  BdfStats stats;
  double ignition_time = -1.0;  // time of maximum dT/dt (adiabatic models), -1 if none
  double max_dTdt = 0.0;        // K/s
};

ReactorResult run_reactor(const Mechanism& m, const ReactorSetup& setup, Diagnostics& diag);

const char* to_string(ReactorModel model) noexcept;
std::optional<ReactorModel> reactor_model_from_string(std::string_view s) noexcept;  // "iso-v", "iso-p", "adiabatic-v", "adiabatic-p"

}  // namespace omnichem
