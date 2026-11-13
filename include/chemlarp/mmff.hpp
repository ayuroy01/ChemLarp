#pragma once
// Merck Molecular Force Field, MMFF94 and MMFF94s (T. A. Halgren, J. Comput. Chem. 17 (1996) 490, 520, 553,
// 587, 616; 20 (1999) 720, 730).
//
// Atom typing: the 99 numeric MMFF types from the structure (explicit hydrogens, Kekulé bond orders and
// formal charges as written: sulfonyl and phosphoryl groups with S=O / P=O double bonds or as dative
// S+–O− bonds), with MMFF's own aromaticity for 5- and 6-membered rings. Partial charges: formal charges
// shared over resonance-equivalent atoms, then q_i = (1 − M_i u_i) q0_i + Σ_k u_k q0_k + Σ_k ω_ki with
// the bond charge increments (or partial bond charge increments when a pair is not tabulated).
// Energy (kcal/mol, Å): quartic bond stretch, cubic angle bend (linear form for linear centres),
// stretch–bend, Wilson out-of-plane, threefold torsion, buffered 14-7 van der Waals and buffered Coulomb
// (dielectric 1, 1-4 interactions scaled by 0.75), with no cutoffs. Parameters come from Halgren's tables
// with MMFF's step-down equivalences; missing bond, angle, stretch–bend and torsion parameters use the
// empirical rules.

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/optimize.hpp"
#include "chemlarp/system_state.hpp"

namespace chemlarp {

enum class MmffVariant { Mmff94, Mmff94s };

struct MmffEnergyBreakdown {
  double bond = 0, angle = 0, stretch_bend = 0, out_of_plane = 0, torsion = 0, vdw = 0, electrostatic = 0;
  double total() const { return bond + angle + stretch_bend + out_of_plane + torsion + vdw + electrostatic; }
};

class MmffForceField {
 public:
  struct Bond { std::uint32_t i, j; double kb, r0; };
  struct Angle { std::uint32_t i, j, k; double ka, theta0; bool linear; };
  struct StretchBend { std::uint32_t i, j, k; double kijk, kkji, r0ij, r0kj, theta0; };
  struct OutOfPlane { std::uint32_t i, j, k, l; double koop; };  // j centre; l the bond out of the i-j-k plane
  struct Torsion { std::uint32_t i, j, k, l; double V1, V2, V3; };
  struct Pair { std::uint32_t i, j; double R, eps, qq; };  // vdW R*, ε; qq = 332.0716 q_i q_j × scale

  // Types, charges and parameters for one species (all hydrogens explicit). Atom indices are local
  // (atom - atom_begin); coordinates are 3N, x0 y0 z0 x1 ...
  static std::optional<MmffForceField> build(const ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag,
                                             MmffVariant variant = MmffVariant::Mmff94, double dielectric = 1.0);

  double energy(const std::vector<double>& xyz) const;
  double energy_gradient(const std::vector<double>& xyz, std::vector<double>& grad) const;
  MmffEnergyBreakdown breakdown(const std::vector<double>& xyz) const;

  std::size_t atom_count() const noexcept { return types_.size(); }
  const std::vector<int>& types() const noexcept { return types_; }
  const std::vector<std::string>& symbols() const noexcept { return symbols_; }
  const std::vector<double>& charges() const noexcept { return charges_; }
  const std::vector<double>& formal_charges() const noexcept { return q0_; }  // MMFF's shared formal charges
  const std::vector<Bond>& bonds() const noexcept { return bonds_; }
  const std::vector<Angle>& angles() const noexcept { return angles_; }
  const std::vector<StretchBend>& stretch_bends() const noexcept { return stbn_; }
  const std::vector<OutOfPlane>& out_of_planes() const noexcept { return oop_; }
  const std::vector<Torsion>& torsions() const noexcept { return torsions_; }
  std::size_t empirical_parameters() const noexcept { return empirical_; }

 private:
  double eval(const std::vector<double>& xyz, std::vector<double>* grad, MmffEnergyBreakdown* parts) const;
  std::vector<int> types_;
  std::vector<std::string> symbols_;
  std::vector<double> charges_, q0_;
  std::vector<Bond> bonds_;
  std::vector<Angle> angles_;
  std::vector<StretchBend> stbn_;
  std::vector<OutOfPlane> oop_;
  std::vector<Torsion> torsions_;
  std::vector<Pair> pairs_;
  std::size_t empirical_ = 0;
};

// MMFF atom types (numeric) for a species, or nullopt with an error naming an atom MMFF does not cover.
std::optional<std::vector<int>> mmff_atom_types(const ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag);

struct MmffMinimizeResult {
  bool ok = false, converged = false;
  int iterations = 0;
  double energy_initial = 0.0, energy_final = 0.0, grad_max = 0.0;
};

// Embeds a SMILES in 3D (explicit hydrogens), optionally minimizes it with MMFF94/94s, and returns JSON:
// {"ok","variant","elements","bonds":[[a,b,order]],"types","symbols","charges","formal_charges",
//  "energy":{total,bond,angle,stretch_bend,out_of_plane,torsion,vdw,electrostatic},"minimized",
//  "converged","iterations","empirical_parameters","sdf","xyz","diagnostics"} (kcal/mol, e).
std::string mmff_json(std::string_view smiles, MmffVariant variant, bool minimize);

// Minimizes the species' coordinates in `state` in place.
MmffMinimizeResult mmff_minimize(ChemicalSystemState& state, std::uint32_t species, Diagnostics& diag,
                                 MmffVariant variant = MmffVariant::Mmff94, const LbfgsOptions& options = {});

}  // namespace chemlarp
