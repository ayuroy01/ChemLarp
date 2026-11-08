#pragma once
// Ring perception and aromaticity.
//
// SSSR: candidate cycles are the shortest cycle through each ring bond; a GF(2)
// Gaussian elimination over bond-incidence vectors keeps the smallest
// independent set (size = cyclomatic number when the candidates span the cycle
// space — true for ordinary organic ring systems, not guaranteed for cages).
//
// Aromaticity: a cycle is aromatic when every atom is sp2-eligible and the
// pi-electron count is 4n+2 (Hückel). Fused systems are handled by also testing
// the envelopes (symmetric differences) of pairs of fused rings, e.g. azulene.

#include <cstdint>
#include <vector>

#include "omnichem/system_state.hpp"

namespace omnichem {

struct Ring {
  std::vector<std::uint32_t> atoms;  // global atom ids, in ring order
  std::vector<std::uint32_t> bonds;  // global bond ids, sorted
};

std::vector<Ring> find_sssr(const ChemicalSystemState& state, std::uint32_t species);

// Pi electrons atom `a` contributes to a ring, or -1 if it cannot be aromatic.
int aromatic_pi_electrons(const ChemicalSystemState& state, std::uint32_t atom);

// Recomputes kAtomAromatic / kBondAromatic for the species from its Kekulé
// structure, so equivalent Kekulé and aromatic inputs end up identical.
void perceive_aromaticity(ChemicalSystemState& state, std::uint32_t species);

}  // namespace omnichem
