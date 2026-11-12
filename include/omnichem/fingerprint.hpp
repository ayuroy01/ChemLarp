#pragma once
// Circular (Morgan / ECFP) fingerprints and similarity.
//
// Rogers & Hahn, J. Chem. Inf. Model. 50 (2010) 742. Hydrogens are folded into
// counts. Each heavy atom starts from the Daylight-style invariants of the
// paper: heavy neighbours, valence minus hydrogens, atomic number, mass number
// (isotope, else 0), formal charge, hydrogens, ring membership. Iteration k
// combines an atom's identifier with the sorted (bond type, neighbour
// identifier) pairs; aromatic bonds are a type of their own. An environment that
// covers exactly the same bonds as one already kept is dropped (structural
// duplicate removal), so radius 2 gives ECFP4. Identifiers are 32-bit (FNV-1a
// over 32-bit words, finalised with MurmurHash3's fmix32), and bits are the
// identifier modulo the length. The hash differs from other toolkits, so bit
// positions are not comparable with theirs; similarities are.
//
// use_chirality adds the CIP label of stereocentres to their initial invariant.

#include <cstdint>
#include <utility>
#include <vector>

#include "omnichem/system_state.hpp"

namespace omnichem {

struct MorganOptions {
  int radius = 2;           // 2 = ECFP4, 3 = ECFP6
  int bits = 2048;          // folded length
  bool use_chirality = false;
};

// Unfolded environments: (identifier, count), sorted by identifier.
std::vector<std::pair<std::uint32_t, std::uint32_t>> morgan_environments(const ChemicalSystemState& state,
                                                                         std::uint32_t species,
                                                                         const MorganOptions& options = {});

struct Fingerprint {
  int bits = 0;
  std::vector<std::uint64_t> words;  // bit i is words[i / 64] >> (i % 64)
  int count() const;                 // bits set
  bool test(int i) const { return (words[static_cast<std::size_t>(i) / 64] >> (i % 64)) & 1u; }
};

Fingerprint morgan_fingerprint(const ChemicalSystemState& state, std::uint32_t species, const MorganOptions& options = {});

// Both 0 when neither fingerprint has a bit set. Fingerprints must have equal length.
double tanimoto(const Fingerprint& a, const Fingerprint& b);
double dice(const Fingerprint& a, const Fingerprint& b);

}  // namespace omnichem
