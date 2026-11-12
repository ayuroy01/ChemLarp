#pragma once
// Reaction equations and stoichiometry validation.
//
// Syntax (tokens separated by whitespace):
//   [coef] SMILES + [coef] SMILES ... -> [coef] SMILES + ...
// '>>' is accepted in place of '->'. Coefficients are positive reals (default 1).
// Example: "1 C + 2 O=O -> 1 O=C=O + 2 O"

#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/system_state.hpp"

namespace chemlarp {

struct ReactionTerm {
  double coefficient = 1.0;
  std::uint32_t species = 0;
  std::size_t position = 0;  // offset of the SMILES token in the equation string
};

struct Reaction {
  std::vector<ReactionTerm> reactants;
  std::vector<ReactionTerm> products;
};

// Parses each SMILES into `state` (one species per term). Returns nullopt if any
// term fails to parse or the equation is malformed; all errors go to `diag`.
std::optional<Reaction> parse_reaction(ChemicalSystemState& state, std::string_view equation, Diagnostics& diag);

// Element-by-element and charge balance. Emits one ElementImbalance per unbalanced
// element (delta = products - reactants) and ChargeImbalance. Returns true if balanced.
bool check_balance(const ChemicalSystemState& state, const Reaction& rxn, Diagnostics& diag,
                   double tolerance = 1e-9);

struct ReactionThermo {
  double dG_RT = 0.0;  // Δr G° / RT
  double dH_J = 0.0;   // Δr H°, J/mol
  double dS_J = 0.0;   // Δr S°, J/(mol K)
  double K = 0.0;      // equilibrium constant, standard state 1 bar
};

// Requires state.g_RT bound at `temperature_K` (call bind_thermo first) — or binds it itself.
std::optional<ReactionThermo> reaction_thermo(const ChemicalSystemState& state, const Reaction& rxn,
                                              double temperature_K, Diagnostics& diag);

}  // namespace chemlarp
