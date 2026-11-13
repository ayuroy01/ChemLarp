#pragma once
// Proteins: PDB chains and residues, protonation at a given pH, hydrogen placement, sequences, isoelectric
// point and molecular weight.
//
// Reading: ATOM/HETATM records of the first model; for atoms with alternate locations, the location with
// the highest occupancy (the first on a tie). SSBOND records and S–S distances below 2.5 Å give disulfides;
// a C(i)–N(i+1) distance above 2.0 Å is a chain break.
//
// Protonation: every titratable group takes the state favoured at the pH by its model pKa (PROPKA 3 model
// values: Asp 3.80, Glu 4.50, His 6.50, Cys 9.00, Tyr 10.00, Lys 10.50, Arg 12.50, N-terminus 8.00,
// C-terminus 3.20; Olsson et al., J. Chem. Theory Comput. 7 (2011) 525), with no shifts from the protein
// environment. Cysteines in disulfides lose HG; a neutral histidine is protonated on the ring nitrogen that
// donates a hydrogen bond to an O or N acceptor (Nδ1: HID, Nε2: HIE; HIE when neither or both do), and a
// nitrogen bound to a metal ion stays unprotonated. Existing hydrogens are removed and all hydrogens are
// rebuilt; a missing C-terminal OXT is added.
//
// Hydrogen geometry (PDB version 3 names): C–H 1.09 Å (aromatic 1.08), N–H 1.01, O–H 0.96, S–H 1.34;
// tetrahedral centres from their heavy neighbours (methyl and NH3+ staggered), planar sp2 groups in plane,
// hydroxyl and thiol hydrogens turned toward the best hydrogen-bond acceptor (else anti), carboxylic acid
// hydrogens syn. Methylene pairs follow the PDB prochirality convention (checked against NMR structures).
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
  bool added = false;  // placed by protonate()
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
  std::vector<ProteinAtom> atoms;
  const ProteinAtom* atom(std::string_view n) const;
  ProteinAtom* atom(std::string_view n);
};

struct Protein {
  std::string title;
  std::vector<ProteinResidue> residues;
  std::vector<std::pair<std::size_t, std::size_t>> disulfides;  // residue indices
  std::vector<std::size_t> chain_breaks;                        // residue index after each break
  std::size_t atom_count() const;
};

std::optional<Protein> read_protein_pdb(std::string_view text, Diagnostics& diag);

// One-letter sequence per chain of the 20 standard amino acids (and their protonation-state names such as
// HID, CYX), in file order; ligands, waters and nonstandard residues are skipped.
std::vector<std::pair<std::string, std::string>> chain_sequences(const Protein& p);

struct ProtonationOptions {
  double pH = 7.4;
  bool add_oxt = true;
};

struct ProtonationReport {
  std::size_t hydrogens_removed = 0, hydrogens_added = 0;
  int net_charge = 0;
  std::vector<std::string> missing_atoms;  // "A 42 LYS NZ"
  std::vector<std::string> unknown_residues;
};

ProtonationReport protonate(Protein& p, const ProtonationOptions& options, Diagnostics& diag);

std::string write_protein_pdb(const Protein& p);

enum class PkaSet { Bjellqvist, PropkaModel };
// Net charge of a sequence (one-letter codes) at a pH, and the pH where it is zero (bisection to 1e-4).
double sequence_net_charge(std::string_view sequence, double pH, PkaSet set = PkaSet::Bjellqvist);
double isoelectric_point(std::string_view sequence, PkaSet set = PkaSet::Bjellqvist);
// Average molecular weight (g/mol) of the chain: residues plus one water.
double sequence_molecular_weight(std::string_view sequence);

// JSON for the C API: reads the PDB text, protonates at the pH and returns {"ok","title","chains":[{id,
// sequence,residues,pI,molecular_weight}],"residues":[{chain,number,icode,name,state,charge}] (titratable
// residues and termini only), "disulfides":[[chain,num,chain,num]],"chain_breaks","missing_atoms",
// "unknown_residues","atoms","hydrogens_added","hydrogens_removed","net_charge","titration":[[pH,charge]],
// "pdb","diagnostics"}.
std::string protein_json(std::string_view pdb_text, double pH);

}  // namespace chemlarp
