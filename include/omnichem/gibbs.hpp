#pragma once
// Multi-species chemical equilibrium by Gibbs free-energy minimisation
// (ideal-gas phase), solved in the element-potential (dual) formulation.
// See docs/MATH.md §4 for the derivation of the Newton system.

#include <cstdint>
#include <vector>

#include "omnichem/diagnostics.hpp"
#include "omnichem/system_state.hpp"

namespace omnichem {

struct EquilibriumInput {
  std::vector<Composition> composition;  // column a_i of the element matrix, per species
  std::vector<double> g_RT;              // mu°_i / RT at the system temperature
  std::vector<double> n0;                // initial moles (defines element totals b = A n0)
  double pressure_Pa = 1.0e5;
};

struct EquilibriumOptions {
  int max_iterations = 200;
  double tolerance = 1e-12;   // on the scaled residual (relative element balance, sum x - 1)
  double max_step = 5.0;      // Newton step clamp in log units
};

struct EquilibriumResult {
  bool converged = false;
  int iterations = 0;
  std::vector<double> n;                   // equilibrium moles
  std::vector<double> mole_fraction;
  std::vector<std::uint8_t> elements;      // atomic numbers of the independent element rows
  std::vector<double> element_potentials;  // lambda_j / RT for each row in `elements`
  double total_moles = 0.0;
  double G_RT = 0.0;                       // sum_i n_i mu_i / RT at equilibrium
  double max_relative_element_error = 0.0;
};

EquilibriumResult minimize_gibbs_ideal_gas(const EquilibriumInput& input, Diagnostics& diag,
                                           const EquilibriumOptions& options = {});

// Convenience: bind thermo for every species at state.conditions, equilibrate from
// state.moles, and write the result back into state.moles (only on convergence).
EquilibriumResult equilibrate(ChemicalSystemState& state, Diagnostics& diag, const EquilibriumOptions& options = {});

// Dense Gaussian elimination with partial pivoting; `m` is row-major n x n.
// Returns false if the matrix is numerically singular.
bool solve_dense(std::vector<double>& m, std::vector<double>& rhs, std::size_t n);

}  // namespace omnichem
