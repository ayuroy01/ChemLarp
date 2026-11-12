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

#include "chemlarp/system_state.hpp"

namespace chemlarp {

struct Ring {
  std::vector<std::uint32_t> atoms;  // global atom ids, in ring order
  std::vector<std::uint32_t> bonds;  // global bond ids, sorted
};

std::vector<Ring> find_sssr(const ChemicalSystemState& state, std::uint32_t species);

// Relevant cycles (Vismara, Electron. J. Combin. 4 (1997) R9): every cycle that
// is not a sum of strictly shorter cycles. Unlike an SSSR this set is unique, so
// it does not depend on atom order; in bridged systems it can hold more rings
// than the cyclomatic number (all three six-rings of bicyclo[2.2.2]octane).
// Candidates are all shortest cycles through each ring bond (up to 64 per bond).
std::vector<Ring> find_relevant_rings(const ChemicalSystemState& state, std::uint32_t species);

// Pi electrons atom `a` contributes to a ring, or -1 if it cannot be aromatic.
int aromatic_pi_electrons(const ChemicalSystemState& state, std::uint32_t atom);

// Recomputes kAtomAromatic / kBondAromatic for the species from its Kekulé
// structure, so equivalent Kekulé and aromatic inputs end up identical.
void perceive_aromaticity(ChemicalSystemState& state, std::uint32_t species);

}  // namespace chemlarp
