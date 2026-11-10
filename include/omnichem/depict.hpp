#pragma once
// 2D structure depiction: coordinates, wedge/hash stereo bonds, SVG and a 2D
// MDL molfile.
//
// Layout follows the usual structure-diagram-generation scheme (Helson,
// "Structure diagram generation", Rev. Comput. Chem. 13, 1999; Clark et al.,
// J. Chem. Inf. Model. 46, 2006):
//   1. ring systems are laid out on their own: the first ring as a regular
//      polygon, every further ring by closing the arc between its already
//      placed atoms (fused, spiro and bridged rings alike), on the less crowded
//      side;
//   2. chains grow breadth-first with 120° zig-zags (180° at sp centres);
//      substituents go into the largest gap that is not inside a ring;
//      attached ring systems are rotated so their exocyclic bond lines up;
//   3. E/Z double bonds are fixed during growth by mirroring the new branch;
//   4. overlaps and bond crossings are relieved by reflecting subtrees across
//      acyclic single bonds;
//   5. each fragment is turned so its long axis is horizontal (snapped to 60°,
//      which keeps ring bonds on the hexagonal lattice).
// Wedges are chosen so the 2D drawing encodes every tetrahedral centre
// (IUPAC 2006 graphical-representation recommendations, ST-1): the narrow end
// sits on the stereocentre, non-ring bonds to non-stereo atoms are preferred,
// and ring-fusion centres get an explicit wedged hydrogen.
// Drawing proportions follow the ACS 1996 document style (line width 4.2% and
// double-bond spacing 18% of the bond length).

#include <array>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "omnichem/diagnostics.hpp"
#include "omnichem/system_state.hpp"

namespace omnichem {

struct DepictOptions {
  double bond_length = 30.0;  // px
  bool color = true;          // colour heteroatom labels (CPK-like); bonds use currentColor
  bool cip_labels = false;    // annotate (R)/(S) and (E)/(Z)
  bool atom_numbers = false;  // annotate 1-based atom numbers (as in diagnostics)
};

struct DepictAtom {
  std::int32_t atom = -1;  // global atom id, or -1 for a hydrogen drawn only to carry a wedge
  double x = 0.0, y = 0.0; // SVG px, y down
  std::string label;       // empty for skeletal carbon
};

enum class BondStereo : std::uint8_t { None = 0, Wedge = 1, Hash = 6 };  // MDL codes

struct DepictBond {
  std::uint32_t a = 0, b = 0;  // indices into Depiction::atoms; for wedges `a` is the narrow end
  std::uint8_t order = 1;
  BondStereo stereo = BondStereo::None;
};

struct Depiction {
  bool ok = false;
  std::vector<DepictAtom> atoms;
  std::vector<DepictBond> bonds;
  double width = 0.0, height = 0.0;  // px
  std::size_t stereo_centres = 0, stereo_centres_ok = 0;
  std::size_t stereo_bonds = 0, stereo_bonds_ok = 0;
  std::string svg;
  std::string molfile;  // MDL V2000, 2D, wedges as bond stereo 1/6
};

// Unit-bond-length coordinates (y up) for every atom of the species, indexed by
// atom - atom_begin.
std::vector<std::array<double, 2>> layout_2d(const ChemicalSystemState& state, std::uint32_t species);

Depiction depict(const ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag,
                 const DepictOptions& options = {});

// Convenience: parse and depict. On a parse error the result has ok == false.
Depiction depict_smiles(std::string_view smiles, Diagnostics& diag, const DepictOptions& options = {});

}  // namespace omnichem
