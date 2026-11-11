#pragma once
// Structure-file readers: MDL molfile (V2000 and V3000), SDF, PDB and CML.
//
// Every record becomes one species (disconnected parts such as counter-ions or
// waters stay in it, joined by '.'). Atoms keep the file's order and
// coordinates. Each record goes through the SMILES parser, so it gets the same
// valence checking, kekulization, aromaticity and canonical SMILES as typed
// input; Species::smiles is set to the canonical SMILES.
//
// Hydrogens. Plain hydrogen atoms (one single bond, no charge or isotope) are
// folded into hydrogen counts unless ReadOptions::keep_hydrogens is set.
// Implicit hydrogens follow the MDL convention: B C N O P S F Cl Br I are
// filled up to their lowest default valence for the given charge and radical;
// other elements get none. A V2000 valence field or V3000 VAL= overrides this,
// and CML hydrogenCount gives the total directly. PDB adds implicit hydrogens
// only when the file lists no hydrogen atoms at all.
//
// Stereo. 3D records: every candidate centre (4-coordinate C, Si, Ge, Sn, B,
// N+, P, As; 3-coordinate S, Se, P, As) and every double bond gets its
// configuration from the coordinates. 2D records: centres come from wedge (1,
// V3000 CFG=1) and hash (6, CFG=3) bonds at their narrow end; double bonds from
// the drawing unless marked "either" (V2000 stereo 3, a wavy bond on a
// neighbour, CFG=2). Candidates that the CIP rules do not find stereogenic are
// dropped, as with SMILES input. Double bonds in rings smaller than 8 and
// cumulated double bonds (allenes, azides) carry no E/Z. In 3D a centre's
// single hydrogen is placed opposite its heavy neighbours rather than taken
// from the file, so a badly placed H cannot invert it; nearly planar centres get
// no configuration. Atom parity fields and the chiral flag are ignored, as in
// most toolkits.
//
// Validation (RDKit 2026.09 reading the same files, structures compared by
// OmniChem canonical SMILES): PubChem 2D SDF 999/999; PubChem 3D SDF 865/866;
// OmniChem 2D molfiles 1498/1499 and 3D SDF 1465/1468 (1500 random PubChem
// compounds). The differences are stereo on cage atoms (adamantane) or =NH
// imines that the CIP model does not keep, hypervalent halogen normalisation,
// and one distorted 3D input. PDB: crambin
// (1CRN) and ubiquitin (1UBQ) give RDKit's InChIKey; 4HHB, which RDKit rejects,
// reads through the residue templates.
//
// PDB. ATOM/HETATM of the first model; the first alternate location of each
// atom. Bonds come from CONECT (a pair listed twice is a double bond) and, for
// atoms without CONECT, from distance (covalent radii + 0.45 A; no bonds to
// metals). Double bonds and charges of the standard amino acids come from
// residue templates; every other residue keeps single bonds and is reported
// with BondOrdersGuessed.
//
// Not read: query features (atom lists, A/Q/* atoms, query bonds), R-groups,
// S-group expansion, V2000 mass differences other than D/T, enhanced stereo,
// XYZ files (no connectivity).

#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "omnichem/diagnostics.hpp"
#include "omnichem/system_state.hpp"

namespace omnichem {

enum class MolFormat : std::uint8_t { Auto, Molfile, Pdb, Cml };

const char* to_string(MolFormat f) noexcept;

// From the file name extension when it is known, otherwise from the content.
MolFormat detect_format(std::string_view text, std::string_view filename = {});

struct ReadOptions {
  MolFormat format = MolFormat::Auto;
  bool keep_hydrogens = false;
  std::size_t max_records = 0;  // 0 = all
};

struct MolRecord {
  std::int32_t species = -1;  // -1 when the record could not be read
  std::string title;
  std::vector<std::pair<std::string, std::string>> properties;  // SDF data items, in file order
  int dimensions = 0;                                           // 0 (no coordinates), 2 or 3
  std::size_t line = 0;                                         // first line of the record, 1-based
};

// The V2000 valence field (columns 49-51) a writer must give an atom so that
// MDL readers (this one included) give it `hydrogens` implicit hydrogens:
// 0 when the default rule already does for B C N O P S F Cl Br I, 15 for
// "valence zero", otherwise bond_order_sum + hydrogens. Other elements always
// get an explicit value, since readers' valence tables differ for them.
int mdl_valence_field(std::uint8_t z, int charge, int bond_order_sum, int hydrogens) noexcept;

// Reads every record. A record that fails leaves the state untouched, keeps
// species = -1 and adds its errors to `diag`; the following records are still
// read. Diagnostic messages start with "line N:".
std::vector<MolRecord> read_molecules(ChemicalSystemState& state, std::string_view text, Diagnostics& diag,
                                      const ReadOptions& options = {});

}  // namespace omnichem
