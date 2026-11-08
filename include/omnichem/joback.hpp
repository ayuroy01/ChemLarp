#pragma once
// Joback group-contribution estimate of ideal-gas ΔfH°(298) and ΔfG°(298),
// derived directly from the molecular graph. Reference: K. G. Joback & R. C. Reid,
// Chem. Eng. Commun. 57 (1987) 233; values as tabulated in Poling, Prausnitz &
// O'Connell, "The Properties of Gases and Liquids", 5th ed., Appendix C.
// Typical accuracy: a few kJ/mol for simple organics, worse for strained or
// highly functionalised molecules.

#include <cstdint>
#include <utility>
#include <vector>

#include "omnichem/diagnostics.hpp"
#include "omnichem/system_state.hpp"

namespace omnichem {

enum class JobackGroup : std::uint8_t {
  CH3, CH2, CH, C, EqCH2, EqCH, EqC, EqCEq, TrCH, TrC,
  RingCH2, RingCH, RingC, RingEqCH, RingEqC,
  F, Cl, Br, I,
  OHAlcohol, OHPhenol, OEther, ORingEther, Ketone, RingKetone, Aldehyde, Acid, Ester, OxoOther,
  NH2, NH, RingNH, N, RingNEq, Nitrile, Nitro,
  Count
};

const char* to_string(JobackGroup g) noexcept;

struct JobackResult {
  bool ok = false;
  double dHf_kJ = 0.0;  // ΔfH°(298.15 K), ideal gas, kJ/mol
  double dGf_kJ = 0.0;  // ΔfG°(298.15 K), ideal gas, kJ/mol
  std::vector<std::pair<JobackGroup, int>> groups;
};

JobackResult joback_estimate(const ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag);

}  // namespace omnichem
