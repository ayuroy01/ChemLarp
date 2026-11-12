#pragma once
// Cubic equations of state and vapour-liquid equilibrium.
//
//   Peng-Robinson (1976):  P = RT/(v-b) - a(T) / (v² + 2bv - b²)
//   Soave-Redlich-Kwong (1972):  P = RT/(v-b) - a(T) / (v(v+b))
// with a_i = Ω_a R²Tc²/Pc · α(T), α = [1 + κ(1-√Tr)]², b_i = Ω_b R Tc/Pc,
// van der Waals one-fluid mixing a = ΣΣ x_i x_j √(a_i a_j)(1-k_ij), b = Σ x_i b_i.
//
// Flash: Michelsen tangent-plane stability test of the feed, then
// Rachford-Rice with successive substitution of K_i = φ_i^L / φ_i^V.
// References: Peng & Robinson, Ind. Eng. Chem. Fundam. 15 (1976) 59;
// Soave, Chem. Eng. Sci. 27 (1972) 1197; Michelsen, Fluid Phase Equil. 9 (1982) 1, 21;
// Rachford & Rice, J. Pet. Technol. 4 (1952) 19. Critical constants: Poling,
// Prausnitz & O'Connell, The Properties of Gases and Liquids, 5th ed., App. A.

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "chemlarp/diagnostics.hpp"

namespace chemlarp {

struct CriticalProperties {
  std::string name;
  double Tc_K = 0.0;
  double Pc_Pa = 0.0;
  double omega = 0.0;  // acentric factor
};

// Built-in table, looked up by structure (any SMILES spelling of the compound).
std::optional<CriticalProperties> critical_properties(std::string_view smiles);

enum class CubicEos { PengRobinson, SoaveRedlichKwong };
enum class Root { Liquid, Vapor, Stable };  // Stable: the root with the lower Gibbs energy

struct EosMixture {
  CubicEos eos = CubicEos::PengRobinson;
  std::vector<CriticalProperties> components;
  std::vector<double> kij;  // n*n binary interaction parameters (default 0)
  std::size_t size() const noexcept { return components.size(); }
};

struct PhaseProperties {
  double Z = 0.0;                // compressibility factor
  double molar_volume = 0.0;     // m^3/mol
  std::vector<double> ln_phi;    // ln fugacity coefficients
  double gibbs_residual_RT = 0;  // G^R / (nRT) = Σ x_i ln φ_i
};

// Real roots of the cubic in Z (ascending), only those with Z > B.
std::vector<double> compressibility_roots(const EosMixture& mix, double T, double P, const std::vector<double>& x);

PhaseProperties phase_properties(const EosMixture& mix, double T, double P, const std::vector<double>& x, Root root);

// Pure-component vapour pressure (equal liquid and vapour fugacity).
std::optional<double> saturation_pressure(const CriticalProperties& c, CubicEos eos, double T, Diagnostics& diag);

enum class PhaseRegion { Liquid, Vapor, Supercritical, TwoPhase };
const char* to_string(PhaseRegion r) noexcept;

struct FlashResult {
  bool ok = false;
  bool converged = false;
  int iterations = 0;
  PhaseRegion region = PhaseRegion::Liquid;
  double vapor_fraction = 0.0;  // β, moles of vapour per mole of feed
  std::vector<double> x, y, K;  // liquid and vapour mole fractions, K = y/x
  PhaseProperties liquid, vapor;
};

// Isothermal flash of feed z (mole fractions, normalised internally) at T (K), P (Pa).
FlashResult flash_tp(const EosMixture& mix, double T, double P, const std::vector<double>& z, Diagnostics& diag);

// Builds a mixture from SMILES; reports NoCriticalData for unknown components.
std::optional<EosMixture> mixture_from_smiles(const std::vector<std::string>& smiles, CubicEos eos, Diagnostics& diag);

}  // namespace chemlarp
