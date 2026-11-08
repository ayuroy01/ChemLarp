#pragma once

#include <string>

#include "omnichem/system_state.hpp"

namespace omnichem {

// Element counts of one species: explicit atoms plus implicit hydrogens.
Composition compute_composition(const ChemicalSystemState& state, std::uint32_t species);

// Hill notation (C, then H, then alphabetical; alphabetical if no carbon), e.g. "C2H6O", "H4N+".
std::string hill_formula(const Composition& c, int net_charge = 0);

double molar_mass(const Composition& c);  // g/mol

}  // namespace omnichem
