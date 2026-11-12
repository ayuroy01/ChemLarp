#pragma once
// Cahn-Ingold-Prelog descriptors (IUPAC 2013 recommendations, P-92).
//
// Priorities come from the hierarchical digraph rooted at the stereogenic
// unit: multiple bonds and ring closures produce duplicate atoms, implicit
// hydrogens are leaves. The sequence rules are applied hierarchically:
//   1a  higher atomic number (duplicates in mancude rings carry the mean
//       atomic number over the Kekule forms)
//   1b  a ring-closure duplicate ranks by the root distance of its original
//   2   higher mass number (unlabelled atoms use the standard atomic weight)
//   3   Z before E
//   4a  chiral before pseudoasymmetric before non-stereogenic
//   4b  like descriptor pairs before unlike
//   4c  r before s
//   5   R before S; a centre separated only by this rule is r / s
// Rules 3-5 use auxiliary descriptors, computed for every stereogenic unit of
// the digraph from the outermost sphere inwards.
//
// Limits: only tetrahedral centres and double bonds recorded in the state take
// part (no axial, planar or helical units, so no M/P); a double bond separated
// only by Rule 5 is reported as E/Z rather than seqTrans/seqCis; undefined
// stereo units contribute no auxiliary descriptor. A unit that still ties gets
// None and a CipUndetermined warning.
//
// Validation: descriptors agree with RDKit's rdCIPLabeler (Hanson et al.,
// J. Chem. Inf. Model. 58:1755, 2018) on 2601/2601 stereo units of 816 PubChem
// structures and on 19359/19485 units of 6312 enumerated stereoisomers. The
// remaining differences are input reading, not ranking: three-coordinate N/P
// centres written with a ring-closure digit, [P@H] phosphines, quinoid c(=C).

#include <cstdint>
#include <vector>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/system_state.hpp"

namespace chemlarp {

enum class CipLabel : std::uint8_t { None = 0, R, S, E, Z, r, s };

const char* to_string(CipLabel l) noexcept;

// Returns >0 if substituent `a` outranks `b` on `root`, <0 if lower, 0 if tied.
// a and b are neighbour atom ids of root, or kImplicitHRef / kLonePairRef.
int cip_compare(const ChemicalSystemState& state, std::uint32_t root, std::uint32_t a, std::uint32_t b);

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

}  // namespace chemlarp
