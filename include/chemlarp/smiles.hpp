#pragma once
// OpenSMILES parser -> ChemicalSystemState.
//
// Supported: organic subset, bracket atoms (isotope, element, @/@@, H count,
// charge, atom class), bonds - = # $ : / \, branches, ring closures (0-9, %nn),
// '.' disconnected components, aromatic input with kekulization.
// Not yet supported (reported as UnsupportedFeature): '*' wildcard, @TH/@AL/@SP/@TB/@OH.
//
// Parsing is transactional: on any error the state is left exactly as it was
// and every problem found is appended to `diag`.

#include <cstdint>
#include <string_view>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/system_state.hpp"

namespace chemlarp {

// Returns the new species index, or -1 on error. `base_offset` is added to all
// reported character positions (used when the SMILES is embedded in a larger
// string such as a reaction equation).
std::int32_t parse_smiles(ChemicalSystemState& state, std::string_view smiles, Diagnostics& diag,
                          std::size_t base_offset = 0);

}  // namespace chemlarp
