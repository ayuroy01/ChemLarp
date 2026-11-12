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
//    Stereo elements with a CIP descriptor are written, and so are ring
//    para-stereocentres (see stereogenic_centres).
//
// Known limitations: tie-breaking is canonical whenever the remaining tied
// atoms are symmetry-equivalent, which refinement guarantees for ordinary
// molecules but not for some highly regular graphs (cages). E/Z marks in some
// partly aromatic macrocycles (chlorins) depend on the spelling.

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

// Which entries of state.tetrahedral (all species; others 0) are stereogenic in
// `species` and are written: centres with a CIP descriptor, and ring
// para-stereocentres that CIP leaves undescribed — two symmetry-equivalent ring
// branches, two other different substituents, and one ring through the
// centre, both equivalent branches and another such centre; cage bridgeheads
// (three ring bonds) excluded. This is graph-based (potential) stereo: a mark
// that is redundant for one particular isomer is kept and written canonically. (the outer carbons of 1,3,5-trimethylcyclohexane).
std::vector<char> stereogenic_centres(const ChemicalSystemState& state, std::uint32_t species);

// Canonical rank per species atom (index = global atom - atom_begin); folded hydrogens get UINT32_MAX.
std::vector<std::uint32_t> canonical_ranks(const ChemicalSystemState& state, std::uint32_t species);

}  // namespace omnichem
