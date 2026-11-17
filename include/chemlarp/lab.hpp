// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
#pragma once
// Virtual laboratory: mixing aqueous solutions, titrations, and the safety checks that go with them.
//
// Chemistry: the aqueous model (aqueous.hpp; PHREEQC's phreeqc.dat). Each portion is a volume of solution with
// solutes at given molarities; amounts are c·V and are put into the model as molalities n / (V · 0.99705 kg/L),
// i.e. the mass of water is taken as the volume times the density of water at 25 °C (exact for dilute
// solutions; at 1 mol/L the solute's own volume makes it a few per cent high). Volumes add on mixing.
//   Precipitation: phases that form quickly in the laboratory (carbonates, sulfates, fluorite, metal hydroxides
//     and sulfides, halite, sylvite) precipitate when supersaturated; slow minerals (silicates, dolomite,
//     quartz, oxides) do not.
//   Gases: in an open vessel CO2, H2S and NH3 leave as bubbles once their partial pressure would exceed 1 atm.
//   Heat: the reaction enthalpy is the change of Σ n_j ΔH_j between the separate portions and the mixture
//     (Hess's law; ΔH_j of each species and phase from the temperature dependence of its log K), including
//     what precipitates or bubbles out. The temperature change assumes 4.184 J g⁻¹ K⁻¹ and 1 g/mL for the whole
//     volume, no heat loss, and the mixture is recomputed at the final temperature. Excess (non-ideal) enthalpies
//     and heats of dilution are not included.
// Substances the model does not contain (organic compounds other than carbonate species, hypochlorite,
// permanganate, peroxides, cyanide, metals, ...) are listed as not modelled: they take part in the safety
// rules below but not in the chemistry.
//
// Safety, three kinds, each with its source:
//   computed       from the result: pH ≤ 2 or ≥ 11.5 (UN GHS presumes skin corrosion, chapter 3.2), boiling,
//                  and gases given off (with that gas's GHS classification);
//   incompatibility reactive groups of the reagents present, from the rules stated in National Research
//                  Council, "Prudent Practices in the Laboratory" (updated version, 2011), each with its section;
//   ghs            the GHS classification of each reagent (signal word, pictograms, H-statements) from PubChem
//                  (src/ghs_data.inc, scripts/gen_ghs_data.py).
// The checks are not a risk assessment and their absence does not mean a mixture is safe.

#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "chemlarp/diagnostics.hpp"

namespace chemlarp {

struct LabReagent {
  std::string id, name;
  std::string formula;  // for the aqueous model ("" when not modelled)
  std::string smiles;   // structure (GHS lookup)
  std::vector<std::string> groups;  // reactive groups used by the incompatibility rules
};

// The reagent catalogue (hydrochloric acid, sodium hydroxide, bleach, ...).
const std::vector<LabReagent>& lab_reagents();
const LabReagent* find_lab_reagent(std::string_view id);

struct LabComponent {
  std::string reagent;     // catalogue id, or a formula for the aqueous model ("CaCl2")
  double molarity = 0.0;   // mol/L in its portion
};

struct LabPortion {
  std::string label;
  double volume_mL = 0.0;
  double T_C = 25.0;
  std::vector<LabComponent> solutes;
};

struct LabHazard {
  std::string kind;      // computed | incompatibility | ghs
  std::string severity;  // danger | warning | info
  std::string message;
  std::string source;
  std::vector<std::string> substances;
};

struct LabPrecipitate {
  std::string phase, formula;
  double mol = 0.0, grams = 0.0;
};

struct LabGas {
  std::string phase;
  double mol = 0.0, mL = 0.0;  // volume at 25 °C and 1 atm (ideal gas)
};

struct LabMixResult {
  double volume_mL = 0.0;
  double T_initial_C = 25.0, T_final_C = 25.0;
  double heat_released_kJ = 0.0;  // > 0: exothermic
  double pH = 7.0, ionic_strength = 0.0;
  std::vector<LabPrecipitate> precipitates;
  std::vector<LabGas> gases;
  std::vector<std::pair<std::string, double>> major_species;  // mol/kg water, largest first (up to 10)
  std::vector<std::string> not_modelled;
  std::vector<LabHazard> hazards;
};

std::optional<LabMixResult> lab_mix(const std::vector<LabPortion>& portions, Diagnostics& diag);

struct LabTitrationPoint {
  double volume_mL = 0.0, pH = 7.0;
};

struct LabTitrationResult {
  std::vector<LabTitrationPoint> curve;
  std::vector<LabTitrationPoint> equivalence;  // steepest points of the curve (dpH/dV maxima)
  std::vector<std::string> not_modelled;
  std::vector<LabHazard> hazards;  // of the final mixture
};

// Adds the titrant (one reagent at its molarity) to the analyte in `points` equal steps up to max_mL.
std::optional<LabTitrationResult> lab_titrate(const LabPortion& analyte, const LabComponent& titrant, double max_mL,
                                              int points, Diagnostics& diag);

// Incompatibility rules over catalogue reagents (and water when `aqueous`).
std::vector<LabHazard> lab_incompatibilities(const std::vector<std::string>& reagent_ids, bool aqueous);

struct GhsEntry {
  std::string smiles, name;
  int cid = 0;
  std::string signal;                                   // Danger | Warning | ""
  std::vector<std::string> pictograms;                  // GHS01 ... GHS09
  std::vector<std::pair<std::string, double>> hazards;  // H-code, per cent of notifications
};

// By canonical SMILES (ChemLarp's); nullopt when the substance is not in the bundled table.
std::optional<GhsEntry> ghs_lookup(std::string_view canonical_smiles);
// The UN GHS hazard statement for an H-code ("H314" -> "Causes severe skin burns and eye damage"); "" if unknown.
std::string ghs_phrase(std::string_view code);
std::size_t ghs_table_size();

// JSON for the C API. Text formats:
//   mix:        portion <label> <volume mL> [<T °C>]   then lines   <reagent id or formula> <mol/L>
//   titration:  analyte <volume mL> [<T °C>]  then solute lines;  titrant <reagent> <mol/L>;  max <mL>;  points <n>
std::string lab_mix_json(std::string_view text);
std::string lab_titration_json(std::string_view text);
std::string lab_reagents_json();
std::string ghs_json(std::string_view smiles);

}  // namespace chemlarp
