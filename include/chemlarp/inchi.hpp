#pragma once
// IUPAC International Chemical Identifier (standard InChI, version 1) and
// InChIKey, generated from scratch.
//
// Clean-room implementation after the published descriptions (Heller et al.,
// J. Cheminform. 7:23, 2015; InChI Technical Manual): no code from the IUPAC
// reference library is used. Agreement with the reference software is checked
// against a table of reference strings (tests/data/inchi_reference.tsv).
//
// Pipeline per connected component:
//   1. normalisation: explicit hydrogens folded into counts, adjacent +/- pairs
//      neutralised, bonds to metals disconnected, (de)protonation to the neutral
//      form (reported in /p);
//   2. mobile-hydrogen groups (1,3 shifts X(H)-Q=Y and ring 1,5 shifts, merged
//      transitively);
//   3. canonical numbering: colour refinement on (element, connections), then
//      a search over individualisations that keeps the numberings with the
//      smallest connection table, hydrogen vector, isotope and stereo vectors;
//   4. layers: formula, /c, /h, /q, /p, /b, /t, /m, /s, /i.
// The InChIKey is the standard 27-character hash (SHA-256 based).
//
// Agreement with the reference software: 1483 of 1500 (98.9%) random PubChem
// compounds that were never used for tuning give the identical string and key.
// Treat the output as an identifier candidate, not as a certified InChI.
//
// Not covered: polymers, reconnected-metal layer (/r), fixed-H layer (/f),
// isotopic stereo and mobile isotopic H, allene / cumulene stereo, and several
// normalisation rules of the reference software (diazonium, N-oxides, charged
// mobile-H groups in cations, metal complexes of macrocycles).

#include <array>
#include <cstdint>
#include <string>
#include <string_view>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/system_state.hpp"

namespace chemlarp {

struct InchiResult {
  bool ok = false;
  std::string inchi;  // "InChI=1S/..."
  std::string key;    // 27 characters
};

InchiResult inchi(const ChemicalSystemState& state, std::uint32_t species, Diagnostics* diag = nullptr);

// Convenience: parse and identify. On a parse error the result has ok == false.
InchiResult inchi_from_smiles(std::string_view smiles, Diagnostics& diag);

// InChIKey of an InChI string; empty if the string is not an InChI.
std::string inchi_key(std::string_view inchi);

// SHA-256 digest (FIPS 180-4), exposed for tests.
std::array<std::uint8_t, 32> sha256(std::string_view data);

}  // namespace chemlarp
