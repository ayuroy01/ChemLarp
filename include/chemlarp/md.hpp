#pragma once
// Molecular dynamics and normal-mode analysis.
//
// Units: length Å, time fs, mass amu (g/mol), energy kcal/mol (the units of the UFF force field).
// Integrators:
//   NVE       velocity Verlet (Swope et al., J. Chem. Phys. 76 (1982) 637).
//   Bussi     velocity Verlet with stochastic velocity rescaling, a canonical thermostat whose conserved
//             "effective energy" (total energy minus the heat the thermostat has added) checks the
//             integration (Bussi, Donadio & Parrinello, J. Chem. Phys. 126 (2007) 014101).
//   Langevin  BAOAB splitting (Leimkuhler & Matthews, Appl. Math. Res. Express 2013, 34).
// For a molecule in the gas phase the centre-of-mass motion and the overall rotation are removed at the
// start (NVE and Bussi conserve them), so it has 3N − 6 (3N − 5 if linear) kinetic degrees of freedom;
// Langevin dynamics does not conserve momentum and has 3N.
// Random numbers come from a 64-bit Mersenne Twister with the normal and gamma deviates computed here,
// so a run is reproducible across compilers and platforms for a given seed.
//
// Normal modes: Hessian by central differences of the analytic gradient, mass-weighted, with the
// translations and rotations projected out; frequencies in cm⁻¹ (negative: imaginary).
// Vibrational spectrum: Fourier transform of the mass-weighted velocity autocorrelation function.
//
// Lennard-Jones fluid (reduced units ε = σ = m = kB = 1): N particles from an fcc lattice in a cubic
// periodic box, the potential truncated at r_c and shifted to zero there; minimum-image pairs.

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "chemlarp/diagnostics.hpp"

namespace chemlarp {

// Energy (kcal/mol) and gradient (kcal/mol/Å, 3N) at coordinates x (Å, 3N: x0 y0 z0 x1 ...).
using PotentialFunction = std::function<double(const std::vector<double>& x, std::vector<double>& gradient)>;

enum class Thermostat { None, Bussi, Langevin };

struct MdOptions {
  double dt_fs = 0.5;
  std::size_t steps = 2000;
  std::size_t sample_every = 10;   // frames and energies recorded every this many steps
  double temperature_K = 300.0;    // target (thermostats) and initial velocities
  Thermostat thermostat = Thermostat::Bussi;
  double tau_fs = 100.0;           // Bussi coupling time, or 1/γ for Langevin
  std::uint64_t seed = 1;
  bool keep_positions = true;      // store coordinates in the frames
  bool keep_velocities = false;    // store velocities in the frames (for spectra)
};

struct MdFrame {
  double time_fs = 0.0;
  double kinetic = 0.0, potential = 0.0, total = 0.0;  // kcal/mol
  double conserved = 0.0;   // total energy minus the thermostat's accumulated heat (Bussi); total otherwise
  double temperature_K = 0.0;
  std::vector<double> x, v;  // Å, Å/fs
};

struct MdResult {
  std::vector<MdFrame> frames;
  std::vector<double> x, v;  // final state
  std::size_t degrees_of_freedom = 0;
  double mean_temperature_K = 0.0;
  double mean_potential = 0.0;
  double drift_kcal_mol_ps = 0.0;  // slope of the conserved energy (NVE, Bussi)
  double fluctuation = 0.0;        // standard deviation of the conserved energy, kcal/mol
};

// Runs dynamics from positions x (Å) with masses (amu). Initial velocities are drawn at the temperature
// (Maxwell–Boltzmann) unless `v0` is given.
std::optional<MdResult> run_md(const std::vector<double>& masses, const std::vector<double>& x,
                               const PotentialFunction& potential, const MdOptions& options, Diagnostics& diag,
                               const std::vector<double>* v0 = nullptr);

struct NormalModes {
  std::vector<double> frequencies_cm1;           // ascending; the 3N − 6 (or 3N − 5) vibrations
  std::vector<std::vector<double>> displacements;  // Cartesian displacement of each mode, normalized (3N)
  std::vector<double> reduced_masses;            // amu
};

std::optional<NormalModes> normal_modes(const std::vector<double>& masses, const std::vector<double>& x,
                                        const PotentialFunction& potential, Diagnostics& diag);

// Power spectrum of the mass-weighted velocities stored in frames (keep_velocities), on a wavenumber grid
// from 0 to max_cm1; returns intensity per point (arbitrary units, max 1).
std::vector<double> velocity_spectrum(const MdResult& result, const std::vector<double>& masses, double sample_dt_fs,
                                      double max_cm1, std::size_t points);

// ---- Molecules with UFF ----

struct MoleculeMdResult {
  std::vector<std::string> elements;  // symbol per atom
  std::vector<double> masses;
  struct BondRef {
    std::size_t a, b;
    int order;
  };
  std::vector<BondRef> bonds;  // atom indices and Kekulé bond order
  double minimized_energy = 0.0;      // kcal/mol
  MdResult md;
  std::optional<NormalModes> modes;
};

enum class ForceFieldChoice { Uff, Mmff94, Mmff94s };

// Embeds the SMILES in 3D (explicit hydrogens), minimizes with the force field (UFF, or MMFF94 / MMFF94s
// starting from the UFF structure), optionally computes the normal modes, and runs dynamics.
std::optional<MoleculeMdResult> molecule_md(std::string_view smiles, const MdOptions& options, bool compute_modes,
                                            Diagnostics& diag, ForceFieldChoice force_field = ForceFieldChoice::Uff);

// JSON for the C API: {"ok","elements","bonds":[[a,b,order]],"minimized_energy","dof","mean_temperature_K",
// "drift_kcal_mol_ps","fluctuation","frames":[{t,K,U,E,conserved,T}],"positions":[[...3N] per frame],
// "modes":[cm-1],"spectrum":{"max_cm1","intensity":[...]},"xyz","diagnostics"}.
std::string molecule_md_json(std::string_view smiles, const MdOptions& options, bool compute_modes,
                             ForceFieldChoice force_field = ForceFieldChoice::Uff);

// ---- Lennard-Jones fluid ----

struct LjOptions {
  int cells = 3;             // N = 4 cells³
  double density = 0.8;      // ρ* = N σ³ / V
  double temperature = 1.0;  // T* = kB T / ε
  double cutoff = 2.5;       // r_c / σ
  double dt = 0.004;         // τ = σ sqrt(m/ε)
  std::size_t equilibration = 2000, production = 10000;
  Thermostat thermostat = Thermostat::Bussi;
  double tau = 0.5;
  std::uint64_t seed = 1;
};

struct LjResult {
  std::size_t particles = 0;
  double box = 0.0;
  double energy_per_particle = 0.0, energy_error = 0.0;  // potential energy U/N (shifted potential), block error
  double pressure = 0.0, pressure_error = 0.0;           // virial pressure (truncated forces)
  double temperature = 0.0;                              // mean kinetic temperature
  double drift = 0.0;                                    // conserved-energy drift per particle per unit time
};

std::optional<LjResult> lennard_jones_fluid(const LjOptions& options, Diagnostics& diag);

}  // namespace chemlarp
