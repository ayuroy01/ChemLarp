#pragma once
// JSON emitters for the .clx payloads (see docs/CLX_SPEC.md).

#include <string>
#include <string_view>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/gibbs.hpp"
#include "chemlarp/system_state.hpp"

namespace chemlarp {

std::string json_escape(std::string_view s);
std::string diagnostics_to_json(const Diagnostics& d);
std::string structure_to_json(const ChemicalSystemState& state);  // clx "structure" payload
std::string thermo_to_json(const ChemicalSystemState& state);     // clx "thermo" payload
std::string equilibrium_to_json(const ChemicalSystemState& state, const EquilibriumResult& r);

}  // namespace chemlarp
