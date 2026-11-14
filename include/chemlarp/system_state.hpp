#pragma once
// ChemicalSystemState: the single in-memory representation shared by every stage
// of the pipeline (topology -> thermodynamics -> force field -> process).
//
// Storage is structure-of-arrays (SoA). Each column is a 64-byte aligned vector,
// so per-atom loops (coordinates, forces, masses) compile to packed SIMD. The
// `Atom` / `Bond` structs are value *views* materialised on demand for
// readability; they are not the storage format.
//
// Every species owns a contiguous [atom_begin, atom_begin + atom_count) slice and
// a contiguous bond slice, so a species can be handed to a worker as a pair of
// spans without copying.

#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "chemlarp/aligned.hpp"
#include "chemlarp/elements.hpp"

namespace chemlarp {

enum AtomFlag : std::uint8_t {
  kAtomAromatic = 1u << 0,  // aromatic by perception (Hückel 4n+2), independent of input notation
  kAtomBracket = 1u << 1,   // written as a bracket atom: hydrogen count is explicit
  kAtomInRing = 1u << 2,    // member of at least one ring (set by ring perception)
};

enum BondFlag : std::uint8_t {
  kBondAromatic = 1u << 0,     // aromatic by perception; `order` holds the Kekulé order
  kBondRingClosure = 1u << 1,  // created by a ring-closure digit
  kBondUp = 1u << 2,           // '/' directional single bond (E/Z input)
  kBondDown = 1u << 3,         // '\' directional single bond
  kBondInRing = 1u << 4,       // not a bridge in the molecular graph
};

enum class Chirality : std::int8_t {
  None = 0,
  CounterClockwise = 1,  // '@'  — neighbours viewed from the first, in input order
  Clockwise = 2,         // '@@'
};

enum class Phase : std::uint8_t { Gas, Liquid, Solid, Aqueous };

// --- value views ---------------------------------------------------------------

struct Atom {
  std::uint8_t element = 0;
  std::int8_t formal_charge = 0;
  std::uint8_t implicit_h = 0;
  std::uint8_t flags = 0;
  std::uint16_t isotope = 0;
  Chirality chirality = Chirality::None;
  std::uint16_t atom_map = 0;
  std::uint32_t species = 0;
  double x = 0.0, y = 0.0, z = 0.0;
};

struct Bond {
  std::uint32_t begin = 0;
  std::uint32_t end = 0;
  std::uint8_t order = 1;  // 1..4, Kekulé form (aromatic bonds are kekulized)
  std::uint8_t flags = 0;
  std::uint32_t species = 0;
};

// --- SoA tables -----------------------------------------------------------------

struct AtomTable {
  AlignedVector<double> x, y, z;
  AlignedVector<std::uint8_t> element;
  AlignedVector<std::int8_t> formal_charge;
  AlignedVector<std::uint8_t> implicit_h;
  AlignedVector<std::uint8_t> flags;
  AlignedVector<std::uint16_t> isotope;
  AlignedVector<std::int8_t> chirality;
  AlignedVector<std::uint16_t> atom_map;
  AlignedVector<std::uint32_t> species;

  std::size_t size() const noexcept { return element.size(); }
  void reserve(std::size_t n);
  void push(const Atom& a);
  Atom get(std::size_t i) const;
  void truncate(std::size_t n);
};

struct BondTable {
  AlignedVector<std::uint32_t> begin, end;
  AlignedVector<std::uint8_t> order;
  AlignedVector<std::uint8_t> flags;
  AlignedVector<std::uint32_t> species;

  std::size_t size() const noexcept { return begin.size(); }
  void reserve(std::size_t n);
  void push(const Bond& b);
  Bond get(std::size_t i) const;
  void truncate(std::size_t n);
};

// Compressed-sparse-row adjacency over all atoms. neighbors(i) and bonds_of(i)
// are parallel spans. Rebuilt in O(atoms + bonds) after topology changes.
struct Adjacency {
  AlignedVector<std::uint32_t> offsets;   // size atoms + 1
  AlignedVector<std::uint32_t> neighbor;  // size 2 * bonds
  AlignedVector<std::uint32_t> bond;      // size 2 * bonds

  std::span<const std::uint32_t> neighbors(std::size_t atom) const noexcept {
    return {neighbor.data() + offsets[atom], offsets[atom + 1] - offsets[atom]};
  }
  std::span<const std::uint32_t> bonds_of(std::size_t atom) const noexcept {
    return {bond.data() + offsets[atom], offsets[atom + 1] - offsets[atom]};
  }
};

using Composition = std::array<std::int32_t, kElementSlots>;  // atom count per atomic number

// Stereo is stored relative to explicit neighbour lists, independent of output order.
inline constexpr std::uint32_t kImplicitHRef = 0xFFFFFFFFu;  // the implicit hydrogen of a [C@H] centre
inline constexpr std::uint32_t kLonePairRef = 0xFFFFFFFEu;   // lone pair of a 3-coordinate centre, e.g. [S@](=O)

struct TetrahedralStereo {
  std::uint32_t atom = 0;
  std::array<std::uint32_t, 4> refs{};  // neighbours (global atom ids or k*Ref markers)
  Chirality parity = Chirality::None;   // looking from refs[0], refs[1..3] are CCW (@) or CW (@@)
};

struct DoubleBondStereo {
  std::uint32_t bond = 0;       // the double bond
  std::uint32_t ref_begin = 0;  // a neighbour of bonds.begin[bond]
  std::uint32_t ref_end = 0;    // a neighbour of bonds.end[bond]
  bool cis = false;             // ref_begin and ref_end on the same side
};

// Axial chirality: an allene or longer cumulene with an even number of double bonds (`centre` is its middle
// atom; `begin` and `end` are the terminal atoms), or a hindered single bond between two ring atoms
// (atropisomer; `centre` is kNoAtom and begin–end is the axis bond). refs[0] and refs[1] are bonded to
// `begin`, refs[2] and refs[3] to `end` (kImplicitHRef for an implicit hydrogen of that end). The parity
// uses the tetrahedral convention for the four refs: looking from refs[0], refs[1..3] counterclockwise (@)
// or clockwise (@@).
inline constexpr std::uint32_t kNoAtom = 0xFFFFFFFDu;
struct AxialStereo {
  std::uint32_t begin = 0, end = 0, centre = kNoAtom;
  std::array<std::uint32_t, 4> refs{};
  Chirality parity = Chirality::None;
};

struct Species {
  std::string smiles;
  std::string canonical_smiles;  // unique per structure incl. stereo (see canonical.hpp)
  std::string formula;  // Hill notation, charge suffix if any
  std::uint32_t atom_begin = 0, atom_count = 0;
  std::uint32_t bond_begin = 0, bond_count = 0;
  std::int32_t net_charge = 0;
  double molar_mass = 0.0;  // g/mol
  Phase phase = Phase::Gas;
  Composition composition{};
};

struct ThermoConditions {
  double temperature_K = 298.15;
  double pressure_Pa = 1.0e5;
};

struct ChemicalSystemState {
  AtomTable atoms;
  BondTable bonds;
  Adjacency adjacency;
  std::vector<Species> species;
  std::vector<TetrahedralStereo> tetrahedral;  // all species; filter by atoms.species[atom]
  std::vector<DoubleBondStereo> double_bonds;  // all species; filter by bonds.species[bond]
  std::vector<AxialStereo> axial;              // all species; filter by atoms.species[begin]

  ThermoConditions conditions;
  AlignedVector<double> moles;  // per species, mol
  AlignedVector<double> g_RT;   // per species, standard chemical potential mu°/RT at `conditions`

  std::uint64_t revision = 0;  // bumped on every topology or state mutation

  std::uint32_t append_species(Species s);
  void rebuild_adjacency();
};

}  // namespace chemlarp
