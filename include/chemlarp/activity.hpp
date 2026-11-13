#pragma once
// Activity-coefficient models for liquid mixtures.
//
// NRTL (Renon & Prausnitz, AIChE J. 14 (1968) 135), multicomponent form:
//   ln γ_i = Σ_j x_j τ_ji G_ji / Σ_k x_k G_ki
//          + Σ_j [x_j G_ij / Σ_k x_k G_kj] (τ_ij − Σ_m x_m τ_mj G_mj / Σ_k x_k G_kj),
//   G_ij = exp(−α_ij τ_ij), τ_ii = 0,
//   τ_ij = a_ij + b_ij/T + e_ij ln T + f_ij T,  α_ij = c_ij + d_ij (T − 273.15 K).
//
// UNIFAC, original version (Fredenslund, Jones & Prausnitz, AIChE J. 21 (1975) 1086), with the subgroup
// volumes/areas and main-group interaction parameters published by DDBST (Hansen et al., Ind. Eng. Chem.
// Res. 30 (1991) 2352, and later revisions; see scripts/gen_unifac_data.py):
//   ln γ_i = ln γ_i^C + ln γ_i^R,
//   ln γ_i^C = ln(Φ_i/x_i) + 5 q_i ln(θ_i/Φ_i) + l_i − (Φ_i/x_i) Σ_j x_j l_j,  l_i = 5(r_i − q_i) − (r_i − 1),
//   ln γ_i^R = Σ_k ν_k^(i) (ln Γ_k − ln Γ_k^(i)),
//   ln Γ_k = Q_k [1 − ln(Σ_m Θ_m Ψ_mk) − Σ_m Θ_m Ψ_km / Σ_n Θ_n Ψ_nm],  Ψ_mn = exp(−a_mn / T).
// A pair of main groups without a published a_mn is reported (NoInteractionParameter), never taken as zero.

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "chemlarp/diagnostics.hpp"

namespace chemlarp {

// NRTL parameters for n components; every matrix is n*n, row-major (i*n + j), diagonal ignored.
// Empty matrices mean zero; alpha_c must be given (typically 0.2-0.47 for i != j).
struct NrtlParameters {
  std::size_t n = 0;
  std::vector<double> a, b, e, f;  // τ_ij = a + b/T + e ln T + f T
  std::vector<double> alpha_c, alpha_d;  // α_ij = c + d (T − 273.15)
};

// ln γ_i at temperature T (K) and liquid mole fractions x (normalised internally).
std::optional<std::vector<double>> nrtl_ln_gamma(const NrtlParameters& p, double T, const std::vector<double>& x,
                                                 Diagnostics& diag);

// One UNIFAC subgroup of the original UNIFAC table (DDBST numbering: 1 CH3, 2 CH2, ..., 16 H2O, ...).
struct UnifacSubgroupInfo {
  int id = 0, main_group = 0;
  double R = 0.0, Q = 0.0;
  std::string name;
};
std::optional<UnifacSubgroupInfo> unifac_subgroup(int id);
// Main-group interaction parameter a_mn (K), if published.
std::optional<double> unifac_interaction(int m, int n);

// A component as subgroup counts: {subgroup id, count}.
using UnifacGroups = std::vector<std::pair<int, int>>;

// ln γ_i at temperature T (K) and liquid mole fractions x (normalised internally).
std::optional<std::vector<double>> unifac_ln_gamma(const std::vector<UnifacGroups>& components, double T,
                                                   const std::vector<double>& x, Diagnostics& diag);

// Original UNIFAC subgroups of a compound from its structure (SMILES); nullopt with a diagnostic when an
// atom is not covered by a supported subgroup.
std::optional<UnifacGroups> unifac_groups_from_smiles(const std::string& smiles, Diagnostics& diag);

// ---- vapour pressure and vapour-liquid equilibrium ------------------------------------------------------

// Antoine equation log10(P/Pa) = A − B/(T + C), fitted for Tmin ≤ T ≤ Tmax (Poling, Prausnitz & O'Connell,
// The Properties of Gases and Liquids, 5th ed.; see scripts/gen_antoine_data.py).
struct AntoineCoefficients {
  std::string name;
  double A = 0.0, B = 0.0, C = 0.0, Tmin_K = 0.0, Tmax_K = 0.0;
  double pressure(double T) const;  // Pa
};
std::optional<AntoineCoefficients> antoine_coefficients(const std::string& smiles);

enum class ActivityModel { Unifac, Nrtl };

// A liquid mixture for modified-Raoult VLE: y_i P = x_i γ_i(T, x) P_i^sat(T) (ideal vapour).
struct ActivitySystem {
  ActivityModel model = ActivityModel::Unifac;
  std::vector<UnifacGroups> groups;            // UNIFAC
  NrtlParameters nrtl;                         // NRTL
  std::vector<AntoineCoefficients> antoine;    // one per component
  std::size_t size() const noexcept { return antoine.size(); }
};

// UNIFAC system from SMILES: groups by unifac_groups_from_smiles, vapour pressures from the Antoine table.
std::optional<ActivitySystem> unifac_system_from_smiles(const std::vector<std::string>& smiles, Diagnostics& diag);

std::optional<std::vector<double>> activity_ln_gamma(const ActivitySystem& sys, double T, const std::vector<double>& x,
                                                     Diagnostics& diag);

struct VleResult {
  double T_K = 0.0, P_Pa = 0.0;
  std::vector<double> x, y, gamma;  // liquid and vapour mole fractions, activity coefficients
  int iterations = 0;
  bool extrapolated = false;        // some P^sat used outside its Antoine range (also a warning)
};

// Bubble point: liquid x given; dew point: vapour y given.
std::optional<VleResult> bubble_pressure(const ActivitySystem& sys, double T, const std::vector<double>& x, Diagnostics& diag);
std::optional<VleResult> bubble_temperature(const ActivitySystem& sys, double P, const std::vector<double>& x, Diagnostics& diag);
std::optional<VleResult> dew_pressure(const ActivitySystem& sys, double T, const std::vector<double>& y, Diagnostics& diag);
std::optional<VleResult> dew_temperature(const ActivitySystem& sys, double P, const std::vector<double>& y, Diagnostics& diag);

}  // namespace chemlarp
