#pragma once
// 3D structure generation: explicit hydrogens, distance-geometry embedding,
// UFF refinement, stereo read-back from coordinates, SDF / XYZ output.
//
// Embedding follows the classical metric-matrix distance geometry approach
// (Crippen & Havel 1988; Blaney & Dixon 1994): bounds from ideal bond lengths
// and angles (UFF), cis/trans 1-4 limits, van der Waals lower bounds;
// Floyd-Warshall triangle smoothing; random distances within the bounds;
// coordinates from the top eigenvectors of the metric matrix in 4D; then a
// bound-violation + chiral-volume + planarity refinement that squeezes out the
// 4th dimension (as in RDKit's embedder), and finally UFF minimisation.
//
// Chirality convention: for refs (a, b, c, d) of a TetrahedralStereo,
// V = (b-a)·[(c-a)×(d-a)] is negative for '@' and positive for '@@'.

#include <cstdint>
#include <string>
#include <string_view>

#include "omnichem/diagnostics.hpp"
#include "omnichem/optimize.hpp"
#include "omnichem/system_state.hpp"
#include "omnichem/uff.hpp"

namespace omnichem {

// Appends a copy of `species` with every implicit hydrogen as an explicit [H]
// atom (stereo references remapped) and returns the new species index.
// Its canonical SMILES equals the original's.
std::int32_t add_explicit_hydrogens(ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag);

struct EmbedOptions {
  std::uint64_t seed = 42;
  int max_attempts = 30;
  int candidates = 5;    // stereo-correct conformers to generate; the lowest UFF energy is kept
  bool minimize = true;  // UFF refinement after embedding
  LbfgsOptions uff{5000, 10, 1e-3, 1e-15, 0.3};
};

struct EmbedResult {
  bool ok = false;
  int attempts = 0;
  double dg_energy = 0.0;  // final distance-geometry violation energy
  UffMinimizeResult uff;
  std::size_t stereo_centres = 0, stereo_centres_ok = 0;
  std::size_t stereo_bonds = 0, stereo_bonds_ok = 0;
};

// Embeds a species with explicit hydrogens; writes x, y, z in `state`.
EmbedResult embed_species(ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag,
                          const EmbedOptions& options = {});

// Convenience: parse, add hydrogens, embed. Returns the 3D species index or -1.
std::int32_t embed_smiles(ChemicalSystemState& state, std::string_view smiles, Diagnostics& diag,
                          EmbedResult* result = nullptr, const EmbedOptions& options = {});

// Stereo read back from coordinates.
Chirality parity_from_coordinates(const ChemicalSystemState& state, const TetrahedralStereo& centre);
bool cis_from_coordinates(const ChemicalSystemState& state, const DoubleBondStereo& bond);
// Dihedral angle a-b-c-d in degrees (-180, 180].
double dihedral_degrees(const ChemicalSystemState& state, std::uint32_t a, std::uint32_t b, std::uint32_t c,
                        std::uint32_t d);

// MDL molfile V2000 (one record ending in $$$$) and XYZ.
std::string write_sdf(const ChemicalSystemState& state, std::uint32_t species, const std::string& title = "");
std::string write_xyz(const ChemicalSystemState& state, std::uint32_t species, const std::string& comment = "");

}  // namespace omnichem
