#pragma once
// Canonical SMILES.
//
// 1. Normalised view: plain [H] atoms attached to one heavy atom are folded
//    into hydrogen counts; aromaticity is the perceived one (rings.hpp), so
//    Kekulé and aromatic spellings coincide.
// 2. Ranking (after Weininger, Weininger & Weininger, JCICS 29 (1989) 97):
//    atom invariants (degree, Z, isotope, charge, H, ring, aromatic, CIP label)
//    are refined by iterated neighbour-rank multisets until stable; remaining
//    ties are broken by splitting the lowest tied class and refining again.
//    With stereo present, each stereo element's parity relative to the ranks
//    (once its deciding neighbours have distinct ranks) is refined in as an
//    invariant, and while any element is still undecided every member of the
//    tied class is tried and the split giving the largest rank-ordered parity
//    vector is kept. Mirror-image spellings of an achiral molecule (cis/trans
//    rings such as C[C@H]1CC[C@@H](O)CC1) therefore give one string.
// 3. Writing: depth-first from the lowest-ranked atom, neighbours visited in
//    rank order, ring closures assigned the lowest free digit. Tetrahedral
//    parity and '/' '\' marks are recomputed for the output neighbour order.
//    Only stereo elements with a defined CIP descriptor are written.
//
// Known limitations: tie-breaking is canonical whenever the remaining tied
// atoms are symmetry-equivalent, which refinement guarantees for ordinary
// molecules but not for some highly regular graphs (cages). Centres without a
// CIP descriptor are not written, so ring stereo the CIP rules leave
// undescribed (e.g. the outer centres of 1,3,5-trisubstituted cyclohexanes) is
// lost and such molecules can canonicalise by spelling. E/Z marks in some
// partly aromatic macrocycles (chlorins) also depend on the spelling.

#include <cstdint>
#include <string>
#include <vector>

#include "omnichem/diagnostics.hpp"
#include "omnichem/system_state.hpp"

namespace omnichem {

struct SmilesWriteOptions {
  bool canonical = true;          // false: random atom order (see random_seed)
  bool kekule = false;            // write upper-case atoms and explicit '=' instead of aromatic notation
  bool stereo = true;             // write @/@@ and '/' '\'
  std::uint64_t random_seed = 0;  // used when canonical == false
};

std::string write_smiles(const ChemicalSystemState& state, std::uint32_t species, const SmilesWriteOptions& options,
                         Diagnostics* diag = nullptr);

std::string canonical_smiles(const ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag);

// Canonical rank per species atom (index = global atom - atom_begin); folded hydrogens get UINT32_MAX.
std::vector<std::uint32_t> canonical_ranks(const ChemicalSystemState& state, std::uint32_t species);

}  // namespace omnichem
