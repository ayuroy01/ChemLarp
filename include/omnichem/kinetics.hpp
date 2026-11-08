#pragma once
// Gas-phase chemical kinetics: mechanisms, mass-action rates, batch reactors.
//
// Three input formats share one reaction model (see mechanism_io.cpp):
//   * .ocm          OmniChem native: species bound to SMILES, CHEMKIN-style reactions
//   * CHEMKIN-II    chem.inp + therm.dat (NASA-7 fixed-format cards), optional NAME->SMILES map
//   * Cantera YAML  species/thermo/reactions as written by ck2yaml or by hand
//
// .ocm format. Species are bound to structures by SMILES, so element balance,
// thermo and identity come from the molecular graph:
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
//     H + O2 (+M) <=> HO2 (+M)   1.475e12  0.60  0          # fall-off: main line is k_inf
//       LOW / 6.366e20 -1.72 524.8 /                        # k_0
//       TROE / 0.8 1e-30 1e30 /                              # a T3 T1 [T2]; omit for Lindemann
//       H2/2.0/ H2O/11.0/ O2/0.78/
//     CH3 + CH3 <=> C2H6     1.0e13 0 0
//       PLOG / 0.1  1.2e12 0 0 /                             # P (atm)  A  b  Ea
//       PLOG / 10.0 5.4e12 0 0 /
//     A + B <=> C + D        1e12 0 1000
//       REV / 3e11 0 500 /                                  # explicit reverse instead of detailed balance
//       DUPLICATE
//   end
//   thermo                                                    # optional NASA-7 cards (CHEMKIN fixed format)
//   ...
//   end
//
// Arrows: "<=>" or "=" reversible, "=>" irreversible. Unless REV is given,
// reverse rate constants come from detailed balance with the species thermo,
// k_r = k_f / K_c, so such mechanisms relax exactly to the Gibbs equilibrium.
// Not supported (UnsupportedReactionType): SRI, Chebyshev, FORD/RORD, LT, JAN,
// FIT1, ionic/electron species, and REV on fall-off or PLOG reactions.

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

enum class RateKind : std::uint8_t {
  Elementary,           // k(T)
  ThirdBody,            // [M] k(T)
  FallOff,              // k_inf Pr/(1+Pr) F,  Pr = k_0 [M]/k_inf   (Lindemann or Troe)
  ChemicallyActivated,  // k_0 1/(1+Pr) F
  PressureLog,          // PLOG: ln k interpolated linearly in ln P
};

// Troe broadening: Fcent = (1-a) exp(-T/T3) + a exp(-T/T1) + exp(-T2/T)
struct TroeParameters {
  double a = 0.0, T3 = 0.0, T1 = 0.0, T2 = 0.0;
  bool has_T2 = false;
  double F_cent(double T) const noexcept;
};

struct PlogLevel {
  double P_Pa = 0.0;
  std::vector<ArrheniusRate> rates;  // summed (duplicate entries at one pressure)
};

struct KineticReaction {
  std::vector<std::pair<std::uint32_t, double>> reactants;  // (species, stoichiometric coefficient)
  std::vector<std::pair<std::uint32_t, double>> products;
  RateKind kind = RateKind::Elementary;
  ArrheniusRate forward;      // Elementary/ThirdBody: k;  FallOff: k_inf;  ChemicallyActivated: k_0
  ArrheniusRate low_or_high;  // FallOff: k_0 (LOW);  ChemicallyActivated: k_inf (HIGH)
  std::optional<TroeParameters> troe;
  std::vector<PlogLevel> plog;             // PressureLog, sorted by pressure
  std::optional<ArrheniusRate> reverse;    // explicit REV (otherwise detailed balance)
  bool reversible = false;
  bool third_body = false;         // rate depends on a collider concentration [M]
  std::vector<double> efficiency;  // per species, used when third_body (default 1)
  bool duplicate = false;
  std::string equation;
  std::size_t line = 0;
};

// Effective forward rate constants (SI), exposed for testing and inspection.
double falloff_rate_constant(const KineticReaction& r, double T, double M_conc, double* dk_dM = nullptr);
double plog_rate_constant(const KineticReaction& r, double T, double P_Pa, double* dk_dP = nullptr);

struct Mechanism {
  ChemicalSystemState state;  // one species per mechanism species, in declaration order
  std::vector<std::string> names;
  std::vector<KineticReaction> reactions;
  std::vector<std::optional<SpeciesThermo>> thermo;  // per species (nullopt if unavailable)

  std::size_t species_count() const noexcept { return names.size(); }
  std::int32_t find(std::string_view name) const;  // case-insensitive, -1 if absent
};

enum class MechanismFormat { Auto, OmniChem, Chemkin, CanteraYaml };

MechanismFormat detect_mechanism_format(std::string_view text) noexcept;

// .ocm text (optionally with a `thermo` section of NASA-7 cards).
std::optional<Mechanism> parse_mechanism(std::string_view text, Diagnostics& diag);

// CHEMKIN-II. `thermo` may be empty when the mechanism has an inline THERMO
// section. `smiles_map` ("NAME SMILES" lines) binds species to structures; species
// without an entry take their composition from their thermo card.
std::optional<Mechanism> import_chemkin(std::string_view mechanism, std::string_view thermo,
                                        std::string_view smiles_map, Diagnostics& diag);

// Cantera YAML (subset: NASA7 thermo; elementary, three-body, falloff (Lindemann/Troe),
// chemically-activated and pressure-dependent-Arrhenius reactions).
std::optional<Mechanism> import_cantera_yaml(std::string_view yaml, std::string_view smiles_map, Diagnostics& diag);

// Any of the above, chosen by content (`format` = Auto) or explicitly.
std::optional<Mechanism> parse_any_mechanism(std::string_view text, std::string_view thermo, std::string_view smiles_map,
                                             MechanismFormat format, Diagnostics& diag);

// Reads files and dispatches on content; thermo_path / smiles_path may be empty.
std::optional<Mechanism> load_mechanism(const std::string& path, Diagnostics& diag, const std::string& thermo_path = "",
                                        const std::string& smiles_path = "");

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
