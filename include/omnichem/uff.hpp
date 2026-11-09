#pragma once
// Universal Force Field (UFF): A. K. Rappé, C. J. Casewit, K. S. Colwell,
// W. A. Goddard III, W. M. Skiff, J. Am. Chem. Soc. 114 (1992) 10024.
//
// Terms (kcal/mol, Å): harmonic bond stretch with bond-order and
// electronegativity corrections; cosine-Fourier angle bend (linear special
// case); torsion by hybridisation class, divided over the torsions about each
// bond; Wilson inversion for sp2 C and N; Lennard-Jones 12-6 van der Waals for
// all pairs beyond 1-3. Electrostatics are not included (as in common UFF use).
// Covered elements: H, B, C, N, O, F, Si, P, S, Cl, Br, I. Hydrogens must be
// explicit (see add_explicit_hydrogens in conformer.hpp).

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "omnichem/diagnostics.hpp"
#include "omnichem/optimize.hpp"
#include "omnichem/system_state.hpp"

namespace omnichem {

struct UffParameters {
  const char* label;
  double r1;      // valence bond radius, Å
  double theta0;  // natural angle, degrees
  double x1;      // vdW distance, Å
  double D1;      // vdW well depth, kcal/mol
  double zeta;    // vdW scale (unused by the 12-6 form)
  double Z1;      // effective charge
  double Vi;      // sp3 torsion barrier, kcal/mol
  double Uj;      // sp2 torsion barrier, kcal/mol
  double Xi;      // GMP electronegativity
};

const UffParameters& uff_parameters(int type);
int uff_type_by_label(const std::string& label);  // -1 if unknown

// Per-atom UFF type (indexed by atom - atom_begin), -1 where typing failed (with an error in diag).
std::vector<int> assign_uff_types(const ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag);

// UFF bond order of bond b: 1, 2, 3, 1.5 aromatic, 1.41 amide C-N.
double uff_bond_order(const ChemicalSystemState& state, std::uint32_t bond);

// Natural bond length r_IJ = r_I + r_J + r_BO - r_EN (Å).
double uff_rest_length(int type_i, int type_j, double bond_order);

struct UffEnergyBreakdown {
  double bond = 0, angle = 0, torsion = 0, inversion = 0, vdw = 0;
  double total() const { return bond + angle + torsion + inversion + vdw; }
};

class UffForceField {
 public:
  struct Bond { std::uint32_t i, j; double k, r0; };
  struct Angle { std::uint32_t i, j, k; double K, c0, c1, c2; bool linear; };
  struct Torsion { std::uint32_t i, j, k, l; double V, cos_nphi0; int n; };
  struct Inversion { std::uint32_t center, j, k, l; double K, c0, c1, c2; };
  struct Vdw { std::uint32_t i, j; double D, x; };

  // Builds terms for one species. Atom indices are local (atom - atom_begin);
  // coordinate vectors are 3 * atom_count long, x0 y0 z0 x1 ...
  static std::optional<UffForceField> build(const ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag);

  double energy(const std::vector<double>& xyz) const;
  double energy_gradient(const std::vector<double>& xyz, std::vector<double>& grad) const;
  UffEnergyBreakdown breakdown(const std::vector<double>& xyz) const;

  std::size_t atom_count() const noexcept { return n_; }
  const std::vector<int>& types() const noexcept { return types_; }
  const std::vector<Bond>& bonds() const noexcept { return bonds_; }
  const std::vector<Angle>& angles() const noexcept { return angles_; }
  const std::vector<Torsion>& torsions() const noexcept { return torsions_; }
  const std::vector<Inversion>& inversions() const noexcept { return inversions_; }
  const std::vector<Vdw>& vdw_pairs() const noexcept { return vdw_; }

 private:
  double eval(const std::vector<double>& xyz, std::vector<double>* grad, UffEnergyBreakdown* parts) const;
  std::size_t n_ = 0;
  std::vector<int> types_;
  std::vector<Bond> bonds_;
  std::vector<Angle> angles_;
  std::vector<Torsion> torsions_;
  std::vector<Inversion> inversions_;
  std::vector<Vdw> vdw_;
};

struct UffMinimizeResult {
  bool ok = false;
  bool converged = false;
  int iterations = 0;
  double energy_initial = 0.0;  // kcal/mol
  double energy_final = 0.0;    // kcal/mol
  double grad_max = 0.0;        // kcal/mol/Å
};

// Minimises the species' coordinates in `state` (x, y, z columns) in place.
UffMinimizeResult uff_minimize(ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag,
                               const LbfgsOptions& options = {});

}  // namespace omnichem
