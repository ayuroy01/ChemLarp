// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
#pragma once
// Free-energy differences from equilibrium samples, and the alchemical coupling of a Lennard-Jones particle.
//
// All estimators work on reduced potentials u = U / kB T, so every free energy here is in units of kB T.
//
//   Exponential averaging (Zwanzig, J. Chem. Phys. 22 (1954) 1420): ΔF = −ln⟨exp(−Δu)⟩₀ from samples of
//     state 0 only; the error is the delta-method estimate, which is optimistic when the overlap is poor.
//   BAR (Bennett, J. Comput. Phys. 22 (1976) 245), in the maximum-likelihood form and with the asymptotic
//     variance of Shirts, Bair, Hooker & Pande, Phys. Rev. Lett. 91 (2003) 140601.
//   MBAR (Shirts & Chodera, J. Chem. Phys. 129 (2008) 124105): the reduced free energies f_k of K states
//     from samples drawn from any of them. Solved by self-consistent iteration followed by Newton steps on
//     the convex objective; uncertainties from the asymptotic covariance Θ = Wᵀ(I − W N Wᵀ)⁺W, evaluated in
//     its K × K form through the eigen-decomposition of N^½ WᵀW N^½. The overlap matrix O = WᵀW N
//     (Klimovich, Shirts & Mobley, J. Comput.-Aided Mol. Des. 29 (2015) 397) is reported with it: an
//     off-diagonal element between neighbouring states below about 0.03 means the result is not reliable.
//   Thermodynamic integration: trapezoidal rule over λ of ⟨∂u/∂λ⟩, errors propagated from the per-window
//     standard errors.
//   Statistical inefficiency g of a correlated time series (Chodera, Swope, Pitera, Seok & Dill, J. Chem.
//     Theory Comput. 3 (2007) 26): g = 1 + 2 Σ (1 − t/N) C(t)/C(0), the sum stopped at the first
//     non-positive autocorrelation; subsampling every ⌈g⌉ samples gives nearly independent ones.
//
// The samples must be (nearly) uncorrelated for the uncertainties to mean anything: subsample first.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "chemlarp/diagnostics.hpp"

namespace chemlarp {

struct FreeEnergyEstimate {
  double df = 0.0;     // kB T
  double error = 0.0;  // one standard error, kB T
};

// w[n] = u_1(x_n) − u_0(x_n) for samples x_n of state 0. ΔF = F_1 − F_0.
FreeEnergyEstimate exp_average(const std::vector<double>& w);

// w_forward[n] = u_1 − u_0 on samples of state 0; w_reverse[m] = u_0 − u_1 on samples of state 1. ΔF = F_1 − F_0.
FreeEnergyEstimate bar(const std::vector<double>& w_forward, const std::vector<double>& w_reverse);

struct MbarResult {
  bool converged = false;
  int iterations = 0;
  std::vector<double> f;                     // reduced free energies, f[0] = 0
  std::vector<std::vector<double>> df, ddf;  // df[i][j] = f_j − f_i and its standard error
  std::vector<std::vector<double>> overlap;  // O[i][j]; each row sums to 1
};

// u[k][n]: reduced potential of sample n evaluated in state k. The samples are grouped by the state they
// were drawn from, n_k[0] of state 0 first, then n_k[1] of state 1, and so on; a state with n_k = 0 is
// evaluated but not sampled (its error is reported as NaN).
MbarResult mbar(const std::vector<std::vector<double>>& u, const std::vector<std::size_t>& n_k, double tolerance = 1e-10,
                int max_iterations = 100000);

// ΔF = F(λ_last) − F(λ_first) from ⟨∂u/∂λ⟩ at each λ (ascending) and its standard error.
FreeEnergyEstimate thermodynamic_integration(const std::vector<double>& lambda, const std::vector<double>& mean,
                                             const std::vector<double>& error);

double statistical_inefficiency(const std::vector<double>& series);
// Indices 0, ⌈g⌉, 2⌈g⌉, ... below n.
std::vector<std::size_t> subsample_indices(std::size_t n, double g);

// ---- A Lennard-Jones particle coupled into the Lennard-Jones fluid ----------------------------------------
//
// Reduced units (ε = σ = m = kB = 1). N particles in a cubic periodic box at the given density; particle 0
// is the solute, whose interaction with every other particle is the soft-core potential of Beutler, Mark,
// van Schaik, Gerber & van Gunsteren (Chem. Phys. Lett. 222 (1994) 529):
//     U(r, λ) = 4λ [ 1/(α(1−λ) + r⁶)² − 1/(α(1−λ) + r⁶) ],   α = 0.5,
// truncated at the cutoff and shifted to zero there, like the solvent–solvent potential. λ = 1 is a full
// solvent particle, λ = 0 a non-interacting one. Each λ window is an independent Langevin (BAOAB) run from
// the fcc lattice. The excess chemical potential μ_ex = F(λ = 1) − F(λ = 0) is reported from MBAR, BAR
// between neighbouring windows, thermodynamic integration and exponential averaging, all on subsampled data.
// As an independent check, Widom test-particle insertion is run in a separate simulation of the pure fluid
// of N − 1 particles in the same box: μ_ex = −kT ln⟨exp(−ΔU/kT)⟩ over random insertions.
struct LjSolvationOptions {
  int cells = 3;  // N = 4 cells³
  double density = 0.8;
  double temperature = 1.2;
  double cutoff = 2.5;
  double dt = 0.004;
  double gamma = 2.0;  // Langevin friction, 1/τ
  // Denser where the soft-core particle opens its cavity (∂u/∂λ peaks near λ = 0.2 and changes sign near 0.4).
  std::vector<double> lambdas = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  std::size_t equilibration = 4000, production = 40000, sample_every = 10;
  std::size_t insertions = 200;  // Widom trial insertions per stored frame; 0 skips the Widom run
  double alpha = 0.5;
  std::uint64_t seed = 1;
  bool keep_matrix = false;  // return the subsampled reduced potentials (u_kn, n_k) for other MBAR programs
};

struct LjSolvationWindow {
  double lambda = 0.0;
  std::size_t samples = 0;      // after subsampling
  double inefficiency = 1.0;    // of the reduced potential difference to the neighbours
  double dudl = 0.0, dudl_error = 0.0;  // ⟨∂u/∂λ⟩ (reduced) and its standard error
};

struct LjSolvationResult {
  std::size_t particles = 0;
  double box = 0.0;
  std::vector<LjSolvationWindow> windows;
  // μ_ex / kT by each estimator.
  FreeEnergyEstimate mbar, bar, ti, exp_forward, widom;
  double min_neighbour_overlap = 0.0;
  bool mbar_converged = false;
  MbarResult mbar_detail;
  std::vector<std::vector<double>> u_kn;  // with keep_matrix: the matrix MBAR was given
  std::vector<std::size_t> n_k;
};

std::optional<LjSolvationResult> lj_solvation(const LjSolvationOptions& options, Diagnostics& diag);

// JSON for the C API.
// mbar_json: text with the sample counts on the first data line (K integers), then K lines of reduced
// potentials (u_kn, N values each); blank lines and lines starting with '#' are skipped.
std::string mbar_json(std::string_view text);
std::string lj_solvation_json(const LjSolvationOptions& options);

}  // namespace chemlarp
