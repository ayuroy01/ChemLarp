#pragma once
// Aqueous solution chemistry: speciation, pH, ionic strength, saturation indices and equilibrium with
// minerals and gases, for dilute to moderately concentrated water (ionic strength up to about 0.5–1 mol/kg).
//
// Thermodynamic data: PHREEQC's phreeqc.dat (USGS, public domain; Parkhurst & Appelo 2013), converted by
// scripts/gen_aqueous_data.py. Every valence state of a redox element is its own component (Fe(2) and
// Fe(3), C(4) and C(-4), N(5) and N(-3), ...); no redox reactions couple them.
//
// Model, as in PHREEQC's ion-association approach:
//   m_j = K_j(T) Π_k a_k^ν_jk · a_H+^ν_jH · a_w^ν_jw / γ_j   for every species j formed from the master
//   species k of the components; mass balance Σ_j ν_jk m_j = total_k; charge balance Σ_j z_j m_j = 0 (pH)
//   or a fixed pH. log K(T) from the analytic expressions or van 't Hoff with ΔH. Activity coefficients:
//   WATEQ Debye–Hückel −A z² √I / (1 + B å √I) + b I where the database gives å and b, else Davies
//   −A z² (√I/(1 + √I) − 0.3 I); 0.1 I for uncharged species; a_w = 1 − 0.017 Σ m. A and B from the
//   dielectric constant (Malmberg–Maryott) and density (Kell) of water at 1 atm, 0–100 °C.
//   The mass of water is held at 1 kg (PHREEQC lets reactions add or remove water).
// Equilibrium phases (minerals or gases): dissolve or precipitate until the saturation index equals the
// target (for a gas, the target is log10 of its partial pressure in atm), limited by the amount present.
// Solution: Newton–Raphson in log10 activities of the master species (with the phase amounts), activity
// coefficients updated each iteration, active set for the phases.
//
// Text format (one statement per line, '#' comment):
//   temperature <°C>                       (default 25)
//   units mol/kgw | mmol/kgw               (default mmol/kgw)
//   pH <value> [charge]                    (charge: pH from the charge balance, value is the start;
//                                           default: charge balance)
//   <component> <total>                    (Ca 1.2, C(4) 2, N(-3) 0.5; an element without a valence
//                                           state takes its usual one: C -> C(4), S -> S(6), Fe -> Fe(2))
//   add <formula> <amount>                 (NaCl 10, CaCO3 1, (NH4)2SO4 0.5, FeCl3 0.1: elements counted,
//                                           valence states from the formula; H and O come from the water
//                                           and the charge balance)
//   equilibrium <phase> [<target SI or log P(atm)>] [<moles available>]   (default 0 and 10)

#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "chemlarp/diagnostics.hpp"

namespace chemlarp {

struct AqueousPhaseTarget {
  std::string name;
  double target_si = 0.0;
  double amount_mol = 10.0;
};

struct AqueousInput {
  double T_C = 25.0;
  std::vector<std::pair<std::string, double>> totals;  // component -> mol/kg water
  bool charge_balance = true;                           // pH from the charge balance; else fixed
  double pH = 7.0;
  std::vector<AqueousPhaseTarget> phases;
};

std::optional<AqueousInput> parse_aqueous(std::string_view text, Diagnostics& diag);

struct AqueousSpecies {
  std::string name;
  int charge = 0;
  double molality = 0.0, activity = 0.0, log_gamma = 0.0;
  std::vector<std::pair<std::string, double>> components;  // component -> moles per mole of species
};

struct SaturationIndex {
  std::string phase, formula;
  double si = 0.0, log_iap = 0.0, log_k = 0.0;
};

struct PhaseChange {
  std::string phase;
  double target_si = 0.0;
  double dissolved_mol = 0.0;  // per kg water; negative: precipitated
  double remaining_mol = 0.0;
};

struct AqueousResult {
  double T_K = 298.15, pH = 7.0, ionic_strength = 0.0, water_activity = 1.0;
  double A_dh = 0.0, B_dh = 0.0;           // Debye–Hückel A (kg^0.5 mol^-0.5) and B (kg^0.5 mol^-0.5 Å^-1)
  double charge_balance_eq = 0.0;          // Σ z m, eq/kg water
  double alkalinity_eq = 0.0;              // eq/kg water (PHREEQC's definition)
  std::vector<AqueousSpecies> species;     // by decreasing molality
  std::vector<std::pair<std::string, double>> totals;  // component totals after phase reactions
  std::vector<SaturationIndex> saturation;  // every phase whose components are all present
  std::vector<PhaseChange> phases;
  int iterations = 0;
};

std::optional<AqueousResult> aqueous_equilibrium(const AqueousInput& input, Diagnostics& diag);

// Parses and solves a solution text; JSON {"ok","T_K","pH","ionic_strength","water_activity","A_dh","B_dh",
// "charge_balance_eq","alkalinity_eq","iterations","totals":[{component,molality}],
// "species":[{name,charge,molality,activity,log_gamma,components:{...}}],"saturation":[{phase,formula,si,log_iap,log_k}],
// "phases":[{phase,target_si,dissolved_mol,remaining_mol}],"diagnostics"}.
std::string aqueous_json(std::string_view text);

// Names of the components and phases in the database.
std::vector<std::string> aqueous_components();
std::vector<std::string> aqueous_phases();

}  // namespace chemlarp
