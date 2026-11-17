// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
#pragma once
// Explicit-water simulation and hydration free energies.
//
// Water: TIP3P (Jorgensen, Chandrasekhar, Madura, Impey & Klein, J. Chem. Phys. 79 (1983) 926): rigid, r(OH)
// 0.9572 Å, HOH 104.52°, q(O) −0.834, q(H) +0.417, Lennard-Jones on O only (A = 582000 kcal Å¹²/mol,
// B = 595 kcal Å⁶/mol, i.e. σ 3.15066 Å, ε 0.152 kcal/mol). Rigid geometry by SETTLE (Miyamoto & Kollman,
// J. Comput. Chem. 13 (1992) 952) for positions and an exact three-constraint solve for velocities.
// Electrostatics: Ewald with smooth PME for the reciprocal part (ewald.hpp), real space within the cutoff.
// Lennard-Jones: truncated at the cutoff with the analytic long-range dispersion correction to energy and
// pressure (uniform fluid beyond the cutoff).
// Solute: an organic molecule from SMILES, flexible, with MMFF94 for everything inside the molecule
// (bonded terms, its own van der Waals and electrostatics) and MMFF94 partial charges; its Lennard-Jones
// interaction with water uses UFF's 12-6 parameters (x_i as r_min, D_i as ε; Lorentz–Berthelot mixing with
// TIP3P oxygen). This combination is not a published force field: the free energies it gives are that model's.
// The solute can also be a TIP3P water ("water"), for the model-consistent check of water in water.
// Dynamics: Langevin, BAOAB splitting with constraints (Leimkuhler & Matthews, Proc. R. Soc. A 472 (2016)
// 20160138), 2 fs; the solute's hydrogens weigh 3 amu (hydrogen-mass repartitioning: the mass comes from the
// atom they are bonded to), which changes the dynamics, not the configurational ensemble the free energies
// depend on. Pressure: Monte Carlo barostat (molecule centres scaled; Åqvist, Wennerström, Nervall, Bjelic &
// Brandsdal, Chem. Phys. Lett. 384 (2004) 288).
//
// Hydration free energy: the solute is decoupled from the water — charges first (λ_elec 1 → 0, linear), then
// van der Waals (λ_vdw 1 → 0, soft-core of Beutler et al. 1994, α = 0.5) — while everything inside the
// molecule stays as it is, so the decoupled state is the molecule in the gas phase. ΔG_hyd = G(coupled) −
// G(decoupled), from MBAR over all windows (free_energy.hpp), with BAR and TI as checks. The electrostatic
// energy is quadratic in λ_elec, so three evaluations give it at every λ.

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "chemlarp/diagnostics.hpp"
#include "chemlarp/free_energy.hpp"

namespace chemlarp {

struct SolvationOptions {
  int waters = 400;
  double temperature_K = 298.15;
  double pressure_atm = 1.0;      // <= 0: constant volume
  double density = 0.997;         // g/cm³, starting box
  double cutoff = 9.0;            // Å, real-space Coulomb and Lennard-Jones
  double skin = 1.5;              // Å, neighbour list
  double ewald_tolerance = 1e-5;  // erfc(β cutoff)
  double pme_spacing = 1.0;       // Å
  int pme_order = 6;
  double dt_fs = 2.0;
  double friction_ps = 1.0;       // Langevin γ, 1/ps
  int barostat_every = 25;        // steps between volume moves
  std::uint64_t seed = 1;
  // Decoupling windows (λ_elec, λ_vdw) from coupled to decoupled; empty: the default schedule.
  std::vector<std::pair<double, double>> states;
  double equilibration_ps = 20.0, production_ps = 100.0, sample_every_ps = 0.5;
  double softcore_alpha = 0.5;
};

struct WaterBoxResult {
  int waters = 0;
  double box = 0.0;                   // Å (cube edge), average
  double density = 0.0, density_error = 0.0;           // g/cm³
  double potential_per_molecule = 0.0, potential_error = 0.0;  // kcal/mol
  double temperature = 0.0;           // K, kinetic
  std::vector<std::pair<double, double>> g_oo;  // O–O radial distribution (r Å, g)
  double acceptance = 0.0;            // barostat
};

// Pure TIP3P water, `ps` picoseconds of production after the options' equilibration.
std::optional<WaterBoxResult> simulate_water(const SolvationOptions& options, double ps, Diagnostics& diag);

struct HydrationWindow {
  double lambda_elec = 1.0, lambda_vdw = 1.0;
  std::size_t samples = 0;
  double inefficiency = 1.0;
  double dudl_elec = 0.0, dudl_vdw = 0.0;  // ⟨∂U/∂λ⟩, kcal/mol
};

struct HydrationResult {
  std::string solute;
  int solute_atoms = 0, waters = 0;
  double box = 0.0;
  FreeEnergyEstimate dg_mbar, dg_bar, dg_ti;  // ΔG_hyd in kcal/mol
  double dg_elec = 0.0, dg_vdw = 0.0;         // MBAR, kcal/mol (charging and van der Waals parts)
  double min_neighbour_overlap = 0.0;
  std::vector<HydrationWindow> windows;
};

// `solute`: a SMILES string, or "water" for a TIP3P molecule.
std::optional<HydrationResult> hydration_free_energy(std::string_view solute, const SolvationOptions& options,
                                                     Diagnostics& diag);

// JSON for the C API.
std::string water_box_json(const SolvationOptions& options, double ps);
std::string hydration_json(std::string_view solute, const SolvationOptions& options);

}  // namespace chemlarp
