// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
#pragma once
// Protein–ligand docking: a flexible ligand placed into a rigid receptor.
//
// Scoring. An empirical function of the form published for AutoDock Vina (Trott & Olson, J. Comput. Chem. 31
// (2010) 455), written here from the paper rather than from that program's code. Heavy atoms only; hydrogens
// decide the atom types and are then ignored. Every term is a function of the *surface* distance
//
//     d = r − R_i − R_j
//
// between two atoms, with R the atom's X-Score radius, and the sum runs over pairs within 8 Å (of r):
//
//     gauss1(d)      = exp(−(d / 0.5)²)
//     gauss2(d)      = exp(−((d − 3) / 2)²)
//     repulsion(d)   = d²        for d < 0, else 0
//     hydrophobic(d) = 1 for d < 0.5, 0 for d > 1.5, linear between   (both atoms hydrophobic)
//     hbond(d)       = 1 for d < −0.7, 0 for d > 0, linear between    (one donor, one acceptor)
//
// weighted by −0.035579, −0.005156, +0.840245, −0.035069 and −0.587439 (the published fit). The reported
// affinity divides the intermolecular sum by the conformational-entropy factor of the same fit,
//
//     affinity = intermolecular / (1 + 0.058459 · N_rot)
//
// where N_rot is the number of active rotatable bonds. The search minimises intermolecular + intramolecular
// energy (the same five terms over ligand atom pairs more than three bonds apart that lie in different rigid
// fragments), so that a pose cannot win by folding the ligand into itself, and poses are ranked by that total
// over the same factor; the affinity reported is the intermolecular part alone, as in the paper, with the
// ligand's own energy reported next to it so that a strained pose is visible.
//
// Atom types (X-Score, as the paper describes them): carbon is hydrophobic unless it is bonded to a nitrogen
// or an oxygen; nitrogen and oxygen are always acceptors and are donors as well when they carry a hydrogen;
// sulfur and phosphorus are polar; the halogens are hydrophobic; metal ions count as donors. Radii (Å):
// C 1.9, N 1.8, O 1.7, S 2.0, P 2.1, F 1.5, Cl 1.8, Br 2.0, I 2.2, metal 1.2.
//
// Speed. The receptor is fixed, so its contribution to each ligand atom type is precomputed once on a grid
// (0.375 Å by default) over the search box and read back by trilinear interpolation; a ligand atom outside
// the grid gets a quadratic penalty that pushes it back in. Only the ligand's own pairs are summed directly.
//
// Ligand. From SMILES (embedded by conformer.hpp, explicit hydrogens) or from MOL/SDF text. The embedding
// fixes bond lengths, angles and ring conformations; docking varies only the position, the orientation and
// the torsions of the rotatable bonds — single, acyclic, not an amide C–N, and with at least two heavy
// neighbours at each end. Cutting those bonds splits the molecule into rigid fragments; the largest is the
// root of a torsion tree and every other fragment hangs off it through one torsion.
//
// Search. Monte Carlo over the pose variables (3 translations, a rotation vector, one angle per torsion),
// each step followed by a local L-BFGS minimisation and accepted or rejected by the Metropolis criterion at
// 1.2 kcal/mol on the total energy. While the walk hunts, repulsion is squashed smoothly towards a ceiling —
// 1.5 kcal/mol for a ligand atom's whole interaction with the receptor, 10 for each of the ligand's own
// pairs — so that a pose can slide out of the protein or thread into a buried site instead of being thrown
// back by a gradient of hundreds; a pose worth keeping is then minimised again on the function as it is,
// which is what every reported number comes from. Several independent runs start from random poses and run
// in parallel; their results are pooled, clustered by heavy-atom RMSD and ranked. A run is deterministic in
// the seed and independent of the thread count.
//
// RMSD. Heavy-atom RMSD to a reference pose with the atoms matched by graph isomorphism over elements and
// connectivity, minimised over every matching, so that symmetric groups (a phenyl ring, a carboxylate, a
// tert-butyl) do not inflate the number. Bond orders are ignored, so a reference drawn with another Kekulé
// structure or protonation still matches.
//
// Receptor. read_protein() then protonate() (protein.hpp), which decides the donors. Hetero groups that
// overlap a reference ligand are dropped when one is given; haem, metals, other cofactors, modified chain
// residues and the waters the structure resolves are kept (ReceptorOptions::remove_waters drops the waters).
// Connectivity for the typing comes from distances (covalent radii), since a PDB file carries none.
//
// Limits: the receptor is rigid, there is no explicit solvent or entropy model beyond the N_rot factor, metal
// coordination is only the generic donor term, and tautomers and protonation states are the caller's choice.
// docs/ARCHITECTURE.md reports redocking RMSD on public crystal complexes.

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/optimize.hpp"
#include "chemlarp/protein.hpp"
#include "chemlarp/system_state.hpp"

namespace chemlarp {

// --- scoring function ---------------------------------------------------------------------------

// X-Score atom types: `_H` hydrophobic, `_P` polar (neither donor nor acceptor), `_D` donor, `_A` acceptor,
// `_DA` both.
enum class XsType : std::uint8_t {
  None = 0,
  C_H,    // carbon bonded to no nitrogen or oxygen
  C_P,    // carbon bonded to a nitrogen or an oxygen
  N_P,    // nitrogen with neither a hydrogen nor a free lone pair (nitro, tertiary amine, proline)
  N_D,    // nitrogen with a hydrogen and no free lone pair (amide, ammonium, amidine, pyrrole)
  N_A,    // nitrogen with a free lone pair and no hydrogen (pyridine, imidazole, nitrile)
  N_DA,   // nitrogen with both
  O_A,    // oxygen without hydrogens
  O_DA,   // oxygen carrying at least one hydrogen
  S_P,
  P_P,
  F_H,
  Cl_H,
  Br_H,
  I_H,
  Met_D,  // metal ion
};

const char* to_string(XsType t) noexcept;
double xs_radius(XsType t) noexcept;  // Å; 0 for None
bool xs_hydrophobic(XsType t) noexcept;
bool xs_donor(XsType t) noexcept;
bool xs_acceptor(XsType t) noexcept;

// The type of atom `z` from what it is bonded to: `heteroatom_neighbour` is true when a carbon has a
// nitrogen or oxygen neighbour, `hydrogens` is how many hydrogens it carries, and `connections` is how many
// atoms it is bonded to in all (hydrogens included). Elements the function does not know give None and are
// left out of the score.
XsType xs_type(std::uint8_t z, bool heteroatom_neighbour, int hydrogens, int connections) noexcept;

struct ScoringWeights {
  double gauss1 = -0.035579;
  double gauss2 = -0.005156;
  double repulsion = 0.840245;
  double hydrophobic = -0.035069;
  double hbond = -0.587439;
  double rotatable = 0.058459;  // the N_rot factor of the affinity
  double cutoff = 8.0;          // Å, centre-to-centre
  // How much of the attraction (gauss1, gauss2, hydrophobic, hydrogen bond) the ligand's own pairs keep;
  // repulsion is always kept in full. 1 is the function as published, and is the default: the published
  // function applies itself unchanged inside the ligand. 0 leaves the ligand's own pairs only keeping it
  // from passing through itself, which removes any reward for folding it onto itself — redocking the
  // complexes in docs/ARCHITECTURE.md is no better either way, so the published form stands.
  double internal_attraction = 1.0;
};

struct ScoringTerms {
  double gauss1 = 0, gauss2 = 0, repulsion = 0, hydrophobic = 0, hbond = 0;
  double weighted(const ScoringWeights& w) const noexcept;
};

// The five unweighted terms for a pair of types at centre-to-centre distance r (Å). Beyond the cutoff every
// term is zero. `derivative` (optional) receives d/dr of the weighted sum.
ScoringTerms pair_terms(XsType a, XsType b, double r, const ScoringWeights& w = {}, double* derivative = nullptr);
double pair_energy(XsType a, XsType b, double r, const ScoringWeights& w = {}, double* derivative = nullptr);

// --- receptor -----------------------------------------------------------------------------------

struct ReceptorAtom {
  double x = 0, y = 0, z = 0;
  XsType type = XsType::None;
  std::string residue;  // "A 102 HIS"
  std::string name;     // atom name, e.g. "ND1"
};

struct DockReceptor {
  std::vector<ReceptorAtom> atoms;  // heavy atoms with a type
  std::size_t waters_removed = 0;
  std::size_t hetero_removed = 0;      // groups overlapping the reference ligand
  std::vector<std::string> kept_hetero;  // cofactors and metals kept, as residue labels
  std::size_t donors = 0, acceptors = 0, hydrophobic = 0;
};

struct ReceptorOptions {
  double pH = 7.4;
  bool protonate = true;  // run protonate() to find the donors; false: donors from the file's hydrogens
  bool repair = false;    // complete missing side chains first
  // Waters are kept by default: the ones a crystal structure resolves are part of the site, and deleting
  // them opens cavities that the diffuse part of the score sinks a ligand into. Over the twelve complexes in
  // docs/ARCHITECTURE.md, dropping them costs three of the eight that otherwise redock correctly. Set this
  // when docking into a structure whose waters are not wanted, such as one modelled without them. A water
  // the reference ligand sits on top of is removed either way: the ligand displaces it.
  bool remove_waters = false;
  // Å: a hetero group (a water included) with an atom this close to a reference-ligand atom goes with it.
  double ligand_overlap = 3.0;
};

// Reference-ligand atoms (heavy atoms, Å) used to clear overlapping hetero groups and to place the box.
using PointCloud = std::vector<std::array<double, 3>>;

DockReceptor prepare_receptor(const Protein& p, const ReceptorOptions& options, const PointCloud& reference,
                              Diagnostics& diag);
// The ligand of a residue name in the file ("MK1"), as heavy-atom coordinates; empty when it is not there.
PointCloud protein_ligand_atoms(const Protein& p, std::string_view residue_name);

// --- ligand -------------------------------------------------------------------------------------

// One rotatable bond: turning `angle` about begin→end moves `moves` (the atoms on the `end` side).
struct LigandTorsion {
  std::uint32_t begin = 0, end = 0;  // local atom indices
  std::vector<std::uint32_t> moves;
};

struct DockLigand {
  ChemicalSystemState state;        // the molecule with explicit hydrogens; coordinates are the input pose
  std::uint32_t species = 0;
  std::string smiles;               // canonical
  std::vector<std::uint8_t> element;     // per local atom
  std::vector<XsType> type;              // None for hydrogens and unknown elements
  std::vector<std::uint32_t> heavy;      // local indices of the scored atoms
  std::vector<std::array<double, 3>> coordinates;  // all atoms, input pose
  std::vector<std::pair<std::uint32_t, std::uint32_t>> bonds;
  std::vector<LigandTorsion> torsions;
  std::vector<std::uint32_t> fragment;   // rigid fragment of each atom
  std::size_t heavy_atoms() const noexcept { return heavy.size(); }
  std::size_t rotatable_bonds() const noexcept { return torsions.size(); }
};

// SMILES: parsed, hydrogens added, embedded (conformer.hpp). The embedding's seed is `seed`.
std::optional<DockLigand> ligand_from_smiles(std::string_view smiles, Diagnostics& diag, std::uint64_t seed = 42);
// MOL/SDF text: the first record, which must carry 3D coordinates, kept as they are (not re-embedded), so a
// crystal conformation survives. Hydrogens are used only to decide the atom types, so a record without them
// is read as it stands and the implicit hydrogen counts of the file serve instead.
std::optional<DockLigand> ligand_from_molfile(std::string_view text, Diagnostics& diag, std::uint64_t seed = 42);

// --- docking ------------------------------------------------------------------------------------

struct DockBox {
  double cx = 0, cy = 0, cz = 0;  // centre, Å
  double sx = 0, sy = 0, sz = 0;  // edge lengths, Å
  bool valid() const noexcept { return sx > 0 && sy > 0 && sz > 0; }
};

// The smallest box around `points` with `margin` added on every side.
DockBox box_around(const PointCloud& points, double margin = 4.0);

struct DockOptions {
  DockBox box;
  double grid_spacing = 0.375;  // Å
  int runs = 8;                 // independent Monte Carlo runs (Vina's "exhaustiveness")
  int steps = 0;                // Monte Carlo steps per run; 0: 200 + 100 per torsion
  int modes = 9;                // poses reported at most
  double energy_range = 3.0;    // kcal/mol: poses worse than best + this are dropped
  double min_rmsd = 1.0;        // Å: cluster separation of the reported poses
  double temperature = 1.2;     // kcal/mol, Metropolis
  std::uint64_t seed = 42;
  ScoringWeights weights;
  // The local minimisation after every Monte Carlo step. It runs tens of thousands of times in a job, so the
  // limits are loose: a pose only has to reach its own basin, and the next step perturbs it again.
  LbfgsOptions local{80, 8, 1e-2, 1e-12, 0.5};
  bool refine_only = false;     // skip the Monte Carlo: minimise the input pose where it stands
};

struct DockedPose {
  double affinity = 0;        // kcal/mol, intermolecular / (1 + w·N_rot)
  double intermolecular = 0;  // kcal/mol
  double internal = 0;        // kcal/mol, ligand's own pairs
  ScoringTerms terms;         // unweighted intermolecular terms
  std::vector<std::array<double, 3>> coordinates;  // every ligand atom, in the ligand's atom order
  double rmsd = -1;           // heavy-atom RMSD to the reference pose, or -1 when there is none
};

struct DockResult {
  // Ordered by what the search minimises: the intermolecular and intramolecular energy plus the penalty
  // that keeps the ligand in the box, over the same N_rot factor. So a pose with a strained conformation,
  // or one hanging over the edge of the box, can carry a slightly better `affinity` and still rank lower.
  std::vector<DockedPose> poses;
  std::size_t heavy_atoms = 0, rotatable_bonds = 0;
  std::size_t receptor_atoms = 0;
  DockBox box;
  // The spacing the maps were actually built at: the option's value, or coarser when the box was so large
  // that the maps would not have fitted in memory at it.
  double grid_spacing = 0;
  int runs = 0, steps = 0;
  std::size_t evaluations = 0;  // scoring-function evaluations over all runs
  double seconds = 0;
  // The reference pose minimised where it stands, when one was given: a docked pose worse than this is a
  // search failure, a docked pose better than this with a large RMSD is a scoring failure.
  double reference_affinity = 0, reference_rmsd = -1;
  bool has_reference = false;
};

// `reference` (optional) is a crystal pose of the same molecule: heavy atoms in any order, matched by graph
// isomorphism. It is used only for the RMSD of the result, never by the search.
std::optional<DockResult> dock(const DockReceptor& receptor, const DockLigand& ligand, const DockOptions& options,
                               Diagnostics& diag, const DockLigand* reference = nullptr);

// Heavy-atom RMSD between two poses of the same molecule, minimised over the graph isomorphisms between them
// (elements and connectivity; bond orders ignored). Negative when the two graphs do not match. `matchings`
// receives how many isomorphisms were tried.
double symmetry_rmsd(const DockLigand& a, const PointCloud& a_pose, const DockLigand& b, const PointCloud& b_pose,
                     std::size_t* matchings = nullptr);

// The score of a pose as it stands (no search): coordinates in the ligand's atom order.
std::optional<DockedPose> score_pose(const DockReceptor& receptor, const DockLigand& ligand,
                                     const PointCloud& pose, const ScoringWeights& w = {});

// The pose variables the search works in: three translations of the ligand's heavy-atom centroid, a
// three-component rotation vector about it (axis times angle), and one angle per torsion, in the order of
// DockLigand::torsions. All zero reproduces the input conformation at the origin, so
// `{centroid, 0, 0, 0, 0 ...}` reproduces the input pose.
std::size_t pose_variable_count(const DockLigand& ligand);

// Energy of the pose those variables describe: the intermolecular sum (exact, not read from a grid) plus the
// ligand's own pairs. With `gradient`, also the derivative with respect to each variable; with
// `coordinates`, the atom positions. This is the function the search minimises, minus the penalty that keeps
// the ligand inside the box, and it is exposed so that its gradient can be checked against its energy.
double pose_variable_energy(const DockReceptor& receptor, const DockLigand& ligand,
                            const std::vector<double>& variables, const ScoringWeights& w = {},
                            std::vector<double>* gradient = nullptr, PointCloud* coordinates = nullptr);

// Poses of `ligand` as an SDF file, one record each, with the affinity as a data item.
std::string write_poses_sdf(const DockLigand& ligand, const std::vector<DockedPose>& poses);

// --- JSON for the C API -------------------------------------------------------------------------

struct DockJsonOptions {
  DockOptions dock;
  ReceptorOptions receptor;
  // The reference ligand: either the name of a hetero residue of the receptor file ("MK1"), or MOL/SDF text
  // holding the crystal pose (detected from the content). Either one places the box when none is given and
  // clears the hetero groups sharing that site; the MOL/SDF form is also the pose the RMSD is measured from.
  std::string reference;
  double box_margin = 4.0;  // around the reference ligand when no box is given
  bool poses_sdf = true;
};

// `receptor_text` is a PDB or mmCIF file; `ligand` is a SMILES string or MOL/SDF text (detected).
// {"ok","smiles","heavy_atoms","rotatable_bonds","receptor_atoms","box":{"centre","size"},"runs","steps",
//  "evaluations","seconds","poses":[{affinity,intermolecular,internal,rmsd,terms:{...},
//  coordinates:[[x,y,z]...]}],"reference":{"atoms","smiles","affinity","rmsd"},
//  "receptor":{"waters_removed","hetero_removed","kept_hetero","donors","acceptors","hydrophobic"},
//  "sdf","diagnostics"}.
std::string dock_json(std::string_view receptor_text, std::string_view ligand, const DockJsonOptions& options);

}  // namespace chemlarp
