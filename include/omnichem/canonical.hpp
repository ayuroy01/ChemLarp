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
// 3. Writing: depth-first from the lowest-ranked atom, neighbours visited in
//    rank order, ring closures assigned the lowest free digit. Tetrahedral
//    parity and '/' '\' marks are recomputed for the output neighbour order.
//    Only stereo elements with a defined CIP descriptor are written.
//
// Known limitation: tie-breaking by lowest index is canonical whenever the
// remaining tied atoms are symmetry-equivalent, which refinement guarantees for
// ordinary molecules but not for some highly regular graphs (cages).

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
