#pragma once
// SMARTS substructure queries (Daylight Theory Manual, ch. 4).
//
// Atoms: bracket expressions with ! & , ; (high to low precedence, '&' may be
// left out), and the primitives * a A #n element symbols (upper case aliphatic,
// lower case aromatic) Dn Hn hn Rn rn vn Xn xn +n -n mass $(recursive) @ @@
// :map. Unbracketed atoms: B C N O P S F Cl Br I b c n o p s * A a.
// Bonds: - = # : ~ @ / \ with the same logic operators; no symbol means
// "single or aromatic". Branches, ring closures (0-9, %nn) and '.' between
// fragments (fragments may match the same target component).
//
// Target conventions (the same as typed SMILES): plain [H] atoms are folded into
// hydrogen counts; aromaticity is OmniChem's perception (rings.hpp); D counts
// explicit connections, H all hydrogens, h hydrogens not written as atoms,
// v the Kekulé valence including hydrogens, X = D + H; R and r use the relevant cycles
// (unique, unlike an SSSR; see rings.hpp), x counts ring bonds. Chirality (@, @@) is checked only when
// SmartsMatchOptions::use_chirality is set, as in most toolkits.
//
// Not supported: component-level grouping '(...)', reaction SMARTS '>',
// @? / @TH / @AL classes, ring-closure-only bonds across '.'.

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string_view>
#include <vector>

#include "omnichem/diagnostics.hpp"
#include "omnichem/system_state.hpp"

namespace omnichem {

struct SmartsQueryData;  // compiled form (smarts.cpp)

class SmartsQuery {
 public:
  bool valid() const noexcept { return data_ != nullptr; }
  std::size_t atom_count() const noexcept;
  const SmartsQueryData& data() const { return *data_; }

 private:
  friend SmartsQuery parse_smarts(std::string_view, Diagnostics&);
  std::shared_ptr<const SmartsQueryData> data_;
};

// Compiles a SMARTS string. On error the result is invalid and `diag` says why
// (with positions in the string).
SmartsQuery parse_smarts(std::string_view smarts, Diagnostics& diag);

struct SmartsMatchOptions {
  bool unique = true;            // one match per set of target atoms
  bool use_chirality = false;    // check @ / @@ against the target's stereo
  std::size_t max_matches = 1000;  // 0 = no limit
};

// Every embedding of the query in one species. Each match lists the global
// target atom for each query atom, in query order.
std::vector<std::vector<std::uint32_t>> match_smarts(const ChemicalSystemState& state, std::uint32_t species,
                                                     const SmartsQuery& query, const SmartsMatchOptions& options = {});

bool has_match(const ChemicalSystemState& state, std::uint32_t species, const SmartsQuery& query,
               bool use_chirality = false);

}  // namespace omnichem
