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
// Axial units (allenes, cumulenes with an even number of double bonds, atropisomeric biaryls) get Ra or Sa
// (P-93.5.1.1): the two ligands of each end are ranked by the rules above with that end as the root, the
// near pair precedes the far pair, and the four are then read like a tetrahedral centre (Ra ≡ M, Sa ≡ P);
// either end may be taken as near.
//
// Planar units (cyclophanes) get Rp or Sp by the pilot-atom rule (Cahn, Ingold & Prelog, Angew. Chem. Int.
// Ed. 5 (1966) 385; Prelog & Helmchen, ibid. 21 (1982) 567): the pilot is the highest-ranked out-of-plane atom
// bonded to the plane; from the plane atom a it bonds to, go to the higher-ranked in-plane neighbour b, then to
// the higher-ranked in-plane neighbour c of b (ring atom or substituent); a -> b -> c clockwise seen from the
// pilot is Rp. Metallocenes use Schlögl's convention (Top. Stereochem. 1 (1967) 39), the one in general use:
// seen from the side of the ring away from the metal, the shorter way round the ring from the highest-ranked
// substituent to the second is clockwise for Rp. Helicenes get P or M (helicity rule, P-93.5.1.3): P for a
// right-handed helix, i.e. positive torsion angles along the inner rim. Planar units are read from 3D
// coordinates only; a plane counts as stereogenic only when it carries a substituent or a ring heteroatom.
//
// Limits: axial units give no auxiliary descriptors to other units (rules 4-5 see only tetrahedral centres
// and double bonds); a double bond separated
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

enum class CipLabel : std::uint8_t { None = 0, R, S, E, Z, r, s, Ra, Sa, Rp, Sp, M, P };

const char* to_string(CipLabel l) noexcept;

// Returns >0 if substituent `a` outranks `b` on `root`, <0 if lower, 0 if tied.
// a and b are neighbour atom ids of root, or kImplicitHRef / kLonePairRef.
int cip_compare(const ChemicalSystemState& state, std::uint32_t root, std::uint32_t a, std::uint32_t b);

CipLabel cip_label(const ChemicalSystemState& state, const TetrahedralStereo& centre, Diagnostics* diag = nullptr);
CipLabel cip_label(const ChemicalSystemState& state, const DoubleBondStereo& bond, Diagnostics* diag = nullptr);
CipLabel cip_label(const ChemicalSystemState& state, const AxialStereo& axis, Diagnostics* diag = nullptr);
CipLabel cip_label(const ChemicalSystemState& state, const PlanarStereo& unit, Diagnostics* diag = nullptr);

// Compares two atoms anywhere in a species by their hierarchical digraphs (as if both were ligands of one
// common root): >0 if `a` ranks higher. Used to choose pilot atoms and rank metallocene substituents.
int cip_compare_atoms(const ChemicalSystemState& state, std::uint32_t a, std::uint32_t b);

struct CipAssignment {
  std::uint32_t index;  // atom for centres, bond for double bonds
  CipLabel label;
};

// Labels for every stereo element of one species.
std::vector<CipAssignment> cip_centres(const ChemicalSystemState& state, std::uint32_t species, Diagnostics* diag = nullptr);
std::vector<CipAssignment> cip_double_bonds(const ChemicalSystemState& state, std::uint32_t species,
                                            Diagnostics* diag = nullptr);
// index: the centre atom of an allene, or the axis bond of an atropisomer.
std::vector<CipAssignment> cip_axes(const ChemicalSystemState& state, std::uint32_t species, Diagnostics* diag = nullptr);
// index: position of the unit in state.planar.
std::vector<CipAssignment> cip_planes(const ChemicalSystemState& state, std::uint32_t species, Diagnostics* diag = nullptr);

}  // namespace chemlarp
