#pragma once
// Cahn-Ingold-Prelog descriptors (IUPAC 2013 recommendations, P-92).
//
// Priorities come from a breadth-first exploration of the hierarchical digraph
// rooted at the stereocentre: multiple bonds and ring closures produce duplicate
// atoms; implicit hydrogens are leaves. Implemented: Rule 1a (atomic number) and
// Rule 2 (mass number; unlabelled atoms use the standard atomic weight).
// Not implemented: Rule 1b, Rules 3-5 (like/unlike, r/s pseudo-asymmetry,
// seqcis/seqtrans). When rules 1-2 leave a tie the descriptor is None and a
// CipUndetermined warning is emitted.

#include <cstdint>
#include <vector>

#include "omnichem/diagnostics.hpp"
#include "omnichem/system_state.hpp"

namespace omnichem {

enum class CipLabel : std::uint8_t { None = 0, R, S, E, Z };

const char* to_string(CipLabel l) noexcept;

// Returns >0 if substituent `a` outranks `b` on `root`, <0 if lower, 0 if tied.
// a and b are neighbour atom ids of root, or kImplicitHRef / kLonePairRef.
// `exclude` is a neighbour of root that is not a substituent (the far end of a double bond), or UINT32_MAX.
int cip_compare(const ChemicalSystemState& state, std::uint32_t root, std::uint32_t a, std::uint32_t b,
                std::uint32_t exclude = UINT32_MAX);

CipLabel cip_label(const ChemicalSystemState& state, const TetrahedralStereo& centre, Diagnostics* diag = nullptr);
CipLabel cip_label(const ChemicalSystemState& state, const DoubleBondStereo& bond, Diagnostics* diag = nullptr);

struct CipAssignment {
  std::uint32_t index;  // atom for centres, bond for double bonds
  CipLabel label;
};

// Labels for every stereo element of one species.
std::vector<CipAssignment> cip_centres(const ChemicalSystemState& state, std::uint32_t species, Diagnostics* diag = nullptr);
std::vector<CipAssignment> cip_double_bonds(const ChemicalSystemState& state, std::uint32_t species,
                                            Diagnostics* diag = nullptr);

}  // namespace omnichem
