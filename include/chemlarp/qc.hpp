#pragma once
// Quantum chemistry: Hartree-Fock and Kohn-Sham DFT single points.
//
// Basis sets: STO-3G, 3-21G, 6-31G, 6-31G* (H-Ar; Cartesian d functions, as in
// Pople's definition of 6-31G*). Integrals over contracted Cartesian Gaussians by
// the McMurchie-Davidson scheme (Helgaker, Jorgensen & Olsen, Molecular
// Electronic-Structure Theory, Wiley 2000, ch. 9); the Boys function from a
// Taylor table (6 terms, step 0.1) below T = 40 and its asymptotic form above.
// Two-electron integrals are stored once with 8-fold symmetry and Schwarz
// screening (1e-12).
//
// SCF: canonical orthogonalisation (overlap eigenvalues below 1e-7 dropped),
// core-Hamiltonian guess, DIIS on the commutator FDS - SDF (Pulay, CPL 73 (1980)
// 393). RHF for closed shells; UHF / UKS when the multiplicity is above 1.
//
// DFT: Becke multicentre partition (JCP 88 (1988) 2547) with Bragg-Slater size
// adjustment; Treutler-Ahlrichs M4 radial grid (JCP 102 (1995) 346);
// Gauss-Legendre x trapezoid angular product grid. Functionals: SVWN5 (Slater +
// VWN5, Vosko, Wilk & Nusair, CJP 58 (1980) 1200), PBE (PRL 77 (1996) 3865, with
// PW92 correlation, PRB 45 (1992) 13244), B3LYP (Stephens et al., JPC 98 (1994)
// 11623: 0.2 HF + 0.08 Slater + 0.72 B88 + 0.19 VWN-RPA + 0.81 LYP; B88: Becke,
// PRA 38 (1988) 3098; LYP in the form of Miehlich et al., CPL 157 (1989) 200).
// Potentials come from forward-mode automatic differentiation of each
// functional, so only the energy expressions are written by hand.
//
// Units: energies in hartree, coordinates in angstrom in the state (converted to
// bohr), dipole in debye.

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/system_state.hpp"

namespace chemlarp {

enum class QcMethod : std::uint8_t { HF, SVWN5, PBE, B3LYP };

const char* to_string(QcMethod m) noexcept;
bool parse_qc_method(const std::string& name, QcMethod& out);  // "hf", "rhf", "uhf", "svwn", "lda", "pbe", "b3lyp"

struct QcOptions {
  QcMethod method = QcMethod::HF;
  std::string basis = "STO-3G";  // STO-3G, 3-21G, 6-31G, 6-31G* (case-insensitive; 6-31G(d) accepted)
  int multiplicity = 0;          // 0: lowest consistent with the electron count (singlet or doublet)
  int max_iterations = 128;
  double energy_tolerance = 1e-9;   // hartree, change between iterations
  double density_tolerance = 1e-7;  // max |DIIS error|
  int grid_level = 2;               // 1 coarse, 2 default, 3 fine (DFT only)
  double memory_limit_mb = 2048;    // two-electron integral storage
};

struct QcResult {
  bool ok = false;
  bool converged = false;
  bool unrestricted = false;
  int iterations = 0;
  std::string method, basis;
  int n_basis = 0, n_alpha = 0, n_beta = 0, multiplicity = 1, charge = 0;
  double energy = 0.0;          // total, hartree
  double nuclear_repulsion = 0.0;
  double one_electron = 0.0;
  double two_electron = 0.0;    // Coulomb + (fraction of) exact exchange
  double xc = 0.0;              // DFT exchange-correlation
  double electrons_on_grid = 0.0;
  std::vector<double> orbital_energies_alpha, orbital_energies_beta;  // hartree, ascending
  std::vector<double> mulliken;         // per species atom (state order)
  std::array<double, 3> dipole{};        // debye, origin at 0
  double s_squared = 0.0;                // <S^2> (UHF/UKS)
};

// Single point for one species with 3D coordinates and explicit hydrogens.
QcResult run_qc(const ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag, const QcOptions& options = {});

}  // namespace chemlarp
