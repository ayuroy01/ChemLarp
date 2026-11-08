#pragma once
// JSON emitters for the .ocx payloads (see docs/OCX_SPEC.md).

#include <string>
#include <string_view>

#include "omnichem/diagnostics.hpp"
#include "omnichem/gibbs.hpp"
#include "omnichem/system_state.hpp"

namespace omnichem {

std::string json_escape(std::string_view s);
std::string diagnostics_to_json(const Diagnostics& d);
std::string structure_to_json(const ChemicalSystemState& state);  // ocx "structure" payload
std::string thermo_to_json(const ChemicalSystemState& state);     // ocx "thermo" payload
std::string equilibrium_to_json(const ChemicalSystemState& state, const EquilibriumResult& r);

}  // namespace omnichem
