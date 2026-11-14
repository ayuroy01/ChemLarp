#pragma once
// Proteins: PDB chains and residues, protonation at a given pH, hydrogen placement, sequences, isoelectric
// point and molecular weight.
//
// Records used: ATOM/HETATM of the first model; for atoms with alternate locations, the location with
// the highest occupancy (the first on a tie). SSBOND records and S–S distances below 2.5 Å give disulfides;
// a C(i)–N(i+1) distance above 2.0 Å is a chain break.
//
// Reading: PDB or mmCIF (detected). Missing side-chain atoms and missing residues can be built (repair()).
//
// Protonation: every titratable group takes the state favoured at the pH by its pKa, predicted from the
// structure (predict_pka(); reference values in water are the PROPKA 3 model values: Asp 3.80, Glu 4.50,
// His 6.50, Cys 9.00, Tyr 10.00, Lys 10.50, Arg 12.50, N-terminus 8.00, C-terminus 3.20; Olsson et al.,
// J. Chem. Theory Comput. 7 (2011) 525). Cysteines in disulfides or bound to a metal ion have no HG; a ring
// nitrogen of histidine bound to a metal ion stays unprotonated. Existing hydrogens are removed and all
// hydrogens are rebuilt; a missing C-terminal OXT is added.
//
// Hydrogen-bond network: the groups whose hydrogens are not fixed by the heavy atoms are optimised together:
// Asn and Gln amides and His rings may flip by 180° (N and O, or N and C, are not told apart by X-ray
// data), neutral His takes the Nδ1 (HID) or Nε2 (HIE) tautomer, Ser/Thr/Cys hydroxyl and thiol hydrogens
// rotate in 30° steps, Tyr OH lies in the ring plane on either side, NH3+ groups of Lys and the N-terminus
// rotate, and a protonated carboxyl group carries its hydrogen on either oxygen. Each state is scored by its
// hydrogen bonds (H···acceptor distance and donor–H···acceptor angle) and overlaps with every other atom;
// groups in contact are solved together (exhaustively for networks of up to 10^5 combinations, by iterated
// improvement beyond). A flip must gain more than 0.5 of a hydrogen bond over the deposited orientation.
// Groups for which the best alternative lies within 0.3 of a hydrogen bond are listed as undetermined
// (for example a histidine with no hydrogen-bond partner, which is given the more common HIE form).
//
// Hydrogen geometry (PDB version 3 names): C–H 1.09 Å (aromatic 1.08), N–H 1.01, O–H 0.96, S–H 1.34;
// tetrahedral centres from their heavy neighbours (methyl staggered), planar sp2 groups in plane, carboxylic
// acid hydrogens syn. Methylene pairs follow the PDB prochirality convention (checked against NMR structures).
//
// Sequence properties: net charge and isoelectric point from the composition with the Bjellqvist pKa set
// used by ExPASy Compute pI/Mw (Bjellqvist et al., Electrophoresis 14 (1993) 1023; 15 (1994) 529), or the
// PROPKA model set; average molecular weight from residue formulas and standard atomic weights.

#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "chemlarp/diagnostics.hpp"

namespace chemlarp {

struct ProteinAtom {
  std::string name, element;
  double x = 0, y = 0, z = 0, occupancy = 1.0, bfactor = 0.0;
  int charge = 0;
  bool added = false;     // hydrogen (or OXT) placed by protonate()
  bool modelled = false;  // heavy atom placed by repair(): not from the file
};

struct ProteinResidue {
  std::string name;           // as read (HIS, CYS, ...)
  std::string chain;
  int number = 0;
  char icode = ' ';
  bool hetero = false;        // HETATM residue
  std::string state;          // protonation state after protonate(): ASP/ASH, GLU/GLH, HID/HIE/HIP, LYS/LYN, CYS/CYX/CYM, TYR/TYM ...
  int charge = 0;             // formal charge of the residue in that state (termini included)
  bool n_terminal = false, c_terminal = false;
  bool predicted = false;     // whole residue built by repair() from the sequence: a prediction, not data
  std::vector<ProteinAtom> atoms;
  const ProteinAtom* atom(std::string_view n) const;
  ProteinAtom* atom(std::string_view n);
};

// A residue of the sequence that has no coordinates (REMARK 465, or mmCIF _pdbx_poly_seq_scheme).
struct MissingResidue {
  std::string name, chain;
  int number = 0;
  char icode = ' ';
};

struct Protein {
  std::string title;
  std::vector<ProteinResidue> residues;
  std::vector<std::pair<std::size_t, std::size_t>> disulfides;  // residue indices
  std::vector<std::size_t> chain_breaks;                        // residue index after each break
  std::vector<MissingResidue> missing_residues;                 // listed by the file, without coordinates
  std::size_t atom_count() const;
};

std::optional<Protein> read_protein_pdb(std::string_view text, Diagnostics& diag);
// mmCIF (PDBx): _atom_site with the author chain and residue numbering, so the result matches the PDB file
// of the same entry; chain identifiers may have several characters.
std::optional<Protein> read_protein_cif(std::string_view text, Diagnostics& diag);
// Either format, detected from the content.
std::optional<Protein> read_protein(std::string_view text, Diagnostics& diag);

// Missing heavy atoms and residues.
//
// Side chains: every standard residue with N, CA and C gets its missing side-chain atoms (and O) from
// internal coordinates and a rotamer library, both derived from X-ray structures by
// scripts/gen_protein_geometry.py. Existing atoms are kept; among the library rotamers compatible with them
// the one with the lowest score (−ln frequency + overlaps with all other atoms + a small contact term) is
// taken, and residues rebuilt next to each other are refined together.
//
// Loops: residues the file lists as missing (REMARK 465 / _pdbx_poly_seq_scheme) between two residues that
// have coordinates are built with ideal backbone geometry, closed onto both ends by cyclic coordinate descent
// (Canutescu & Dunbrack, Protein Sci. 12 (2003) 963) from several starting conformations, keeping the one
// with the fewest overlaps. Such residues are marked `predicted`: the chain is connected and free of
// overlaps, but the conformation is not experimental data. Missing residues at chain ends are built the
// same way (without closure) only when `termini` is set.
struct RepairOptions {
  bool side_chains = true;
  bool loops = true;
  bool termini = false;
  int loop_trials = 40;
};

struct RepairReport {
  std::size_t atoms_added = 0;
  std::vector<std::string> completed_residues;  // "A 42 LYS": residues that had atoms missing
  std::vector<std::string> built_residues;      // predicted residues
  std::vector<std::string> not_built;           // with the reason
};

RepairReport repair(Protein& p, const RepairOptions& options, Diagnostics& diag);

// pKa of every titratable group from the structure (heavy atoms only; hydrogens and waters are ignored).
//
//   pKa = reference + desolvation + hydrogen bonds + charge–charge
//
// The reference is the group's value in water; the three terms are linear in geometric features of the
// structure: how buried the group is (heavy atoms within 15.5 Å and within 6 Å), hydrogen-bond partners within
// 3–4 Å of its titratable atoms (backbone NH, backbone C=O, side-chain NH, hydroxyl, amide O), and the
// screened Coulomb sum over the other titratable groups and metal ions within 10 Å, separately for exposed
// and buried pairs, plus their direct contacts (salt bridges). The weights were fitted to measured pKa values (PKAD-2, Ancona et al., Database 2023)
// by scripts/fit_protein_pka.py; docs/ARCHITECTURE.md gives the cross-validated error next to PROPKA 3's on
// the same groups. Arginine keeps its reference value (no measurements to fit).
constexpr int kPkaFeatures = 10;
struct PkaSite {
  std::size_t residue = 0;  // index into Protein::residues
  std::string group;        // ASP GLU HIS CYS TYR LYS ARG N+ C-
  double reference = 0, pka = 0;
  double desolvation = 0, hydrogen_bonds = 0, charges = 0;  // the three terms (pKa units)
  double buried = 0;                                        // 0 exposed … 1 buried
  double features[kPkaFeatures] = {};
};
std::vector<PkaSite> predict_pka(const Protein& p);

struct ProtonationOptions {
  double pH = 7.4;
  bool add_oxt = true;
  bool structure_pka = true;   // pKa of each group from predict_pka(); false: reference values
  bool hydrogen_bonds = true;  // optimise the hydrogen-bond network (below); false: each group on its own
  bool flips = true;           // allow Asn/Gln/His side-chain flips
  double flip_penalty = 0.5;   // what a flip must gain over the deposited orientation (hydrogen-bond units)
};

struct ProtonationReport {
  std::size_t hydrogens_removed = 0, hydrogens_added = 0;
  int net_charge = 0;
  std::vector<std::string> missing_atoms;  // "A 42 LYS NZ"
  std::vector<std::string> unknown_residues;
  std::vector<PkaSite> pka;                // the pKa values the states were chosen from
  std::vector<std::string> flipped;        // "A 42 ASN": amide or imidazole turned by 180°
  std::vector<std::string> undetermined;   // groups whose orientation or tautomer the structure does not decide
};

ProtonationReport protonate(Protein& p, const ProtonationOptions& options, Diagnostics& diag);

// One-letter sequence per chain of the 20 standard amino acids (and their protonation-state names such as
// HID, CYX), in file order; ligands, waters and nonstandard residues are skipped.
std::vector<std::pair<std::string, std::string>> chain_sequences(const Protein& p);



// PDB format; residues built by repair() are written with occupancy 0.00 and counted in a REMARK.
std::string write_protein_pdb(const Protein& p);
// mmCIF (_atom_site only): for structures the PDB format cannot hold (more than 99,999 atoms, chain
// identifiers of several characters, residue numbers beyond 9999).
std::string write_protein_cif(const Protein& p);
bool fits_pdb_format(const Protein& p);

enum class PkaSet { Bjellqvist, PropkaModel };
// Net charge of a sequence (one-letter codes) at a pH, and the pH where it is zero (bisection to 1e-4).
double sequence_net_charge(std::string_view sequence, double pH, PkaSet set = PkaSet::Bjellqvist);
double isoelectric_point(std::string_view sequence, PkaSet set = PkaSet::Bjellqvist);
// Average molecular weight (g/mol) of the chain: residues plus one water.
double sequence_molecular_weight(std::string_view sequence);

// JSON for the C API: reads the text (PDB or mmCIF), repairs, protonates at the pH and returns {"ok","title",
// "chains":[{id,sequence,residues,pI,molecular_weight}],"residues":[{chain,number,icode,name,state,charge,
// n_terminal,c_terminal,predicted}] (titratable residues and termini), "pka":[{chain,number,name,group,
// reference,pKa,desolvation,hydrogen_bonds,charges,buried,features}],"disulfides":[[chain,num,chain,num]],
// "chain_breaks","missing_atoms","unknown_residues","flipped","undetermined","repair":{atoms_added,
// completed_residues,built_residues,not_built,missing_residues_listed,missing_residues_left},"atoms",
// "hydrogens_added","hydrogens_removed","net_charge","pI","titration":[[pH,charge]] (from the sequence),
// "titration_structure":[[pH,charge]] (from the predicted pKa values),"format":"pdb"|"cif","pdb","cif"
// (the structure in the format named; the other is empty),"diagnostics"}.
struct ProteinJsonOptions {
  double pH = 7.4;
  bool repair = true;         // complete side chains
  bool loops = true;          // build missing residues inside chains
  bool termini = false;       // build missing residues at chain ends
  bool flips = true;          // allow Asn/Gln/His flips
  bool structure_pka = true;  // false: reference pKa values
};
std::string protein_json(std::string_view text, const ProteinJsonOptions& options);
std::string protein_json(std::string_view text, double pH);

}  // namespace chemlarp
