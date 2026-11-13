#pragma once
// Pure-component property database and the correlations built on it.
//
// Data: published compilations (critical constants from the IUPAC series, CRC Handbook or the PSRK appendix;
// acentric factors; CRC melting and boiling points, enthalpies of vaporization and ideal-gas formation
// properties; Poling et al. ideal-gas heat-capacity polynomials; Rackett parameters). See
// scripts/gen_property_data.py for the source of every field. A missing value is NaN.
//
// Correlations:
//   ideal-gas Cp/R = a0 + a1 T + a2 T² + a3 T³ + a4 T⁴ (Poling et al.), H and S by exact integration;
//   ΔHvap(T) = ΔHvap(Tb) [(1 − T/Tc)/(1 − Tb/Tc)]^0.38 (Watson, Ind. Eng. Chem. 35 (1943) 398);
//   saturated liquid volume from the VDI Heat Atlas (PPDS) density equation where tabulated, otherwise
//   V = (R Tc/Pc) Z_RA^(1 + (1 − T/Tc)^(2/7)) (Rackett, J. Chem. Eng. Data 15 (1970) 514, with Spencer &
//   Danner's Z_RA; Zc when no Z_RA is tabulated, reported as a warning).
//
// Transport properties (scripts/gen_transport_data.py), each a published correlation with its fitted range:
//   viscosity and thermal conductivity of the liquid and of the low-pressure vapour: Perry's Chemical
//   Engineers' Handbook, 8th ed., Tables 2-312 to 2-315 (DIPPR equations 100, 101, 102);
//   liquid heat capacity: Perry Table 2-153 (DIPPR equation 100);
//   surface tension: Mulero & Cachadiña, J. Phys. Chem. Ref. Data 41 (2012) 043105, else Jasper,
//   J. Phys. Chem. Ref. Data 1 (1972) 841.

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <string_view>

#include "chemlarp/diagnostics.hpp"

namespace chemlarp {

// One temperature correlation. equation: '0' DIPPR 100, '1' DIPPR 101, '2' DIPPR 102,
// 'M' Mulero–Cachadiña (c = σ0, n0, σ1, n1, σ2, n2, Tc), 'J' Jasper (c = a, b; mN/m, °C); '\0' none.
struct TransportCorrelation {
  char equation = '\0';
  double c[7] = {0, 0, 0, 0, 0, 0, 0};
  double Tmin_K = 0.0, Tmax_K = 0.0;
  bool present() const noexcept { return equation != '\0'; }
};

struct PureComponent {
  static constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
  std::string smiles, name, cas;
  double molar_mass = NaN;                       // g/mol
  double Tc_K = NaN, Pc_Pa = NaN, Vc_m3mol = NaN;  // critical constants
  char critical_source = '-';                    // 'I' IUPAC, 'C' CRC, 'P' PSRK appendix, '-' none
  double omega = NaN;                            // acentric factor
  double Tm_K = NaN, Tb_K = NaN;                 // normal melting and boiling points
  double Hvap_Tb_Jmol = NaN;                     // enthalpy of vaporization at Tb
  double cp_a[5] = {NaN, NaN, NaN, NaN, NaN};    // ideal-gas Cp/R polynomial
  double cp_Tmin_K = NaN, cp_Tmax_K = NaN;
  double cp_liquid_298_JmolK = NaN;
  double Z_RA = NaN;                             // Rackett parameter
  double Hf_gas_Jmol = NaN, Gf_gas_Jmol = NaN, S_gas_JmolK = NaN;  // ideal gas, 298.15 K, 1 bar
  TransportCorrelation liquid_viscosity, gas_viscosity, liquid_conductivity, gas_conductivity, surface_tension;
  TransportCorrelation liquid_heat_capacity;
  TransportCorrelation liquid_density;  // 'P': VDI PPDS saturated liquid density (c = ρc, a, b, c, d, Tc, M)
};

// Looked up by structure (any SMILES spelling); nullopt if the compound is not in the table.
std::optional<PureComponent> pure_component(std::string_view smiles);
std::size_t pure_component_count() noexcept;

// Each correlation returns nullopt with NoPropertyData when an input it needs is missing and warns
// PropertyExtrapolated when T is outside the range of the data behind it.
std::optional<double> ideal_gas_cp(const PureComponent& c, double T, Diagnostics& diag);        // J/(mol K)
std::optional<double> ideal_gas_enthalpy(const PureComponent& c, double T, Diagnostics& diag);  // H(T) − H(298.15 K), J/mol
// Absolute ideal-gas entropy at T and P (Pa): S°(298.15 K) + ∫Cp/T dT − R ln(P/1 bar), J/(mol K).
std::optional<double> ideal_gas_entropy(const PureComponent& c, double T, double P, Diagnostics& diag);
std::optional<double> enthalpy_of_vaporization(const PureComponent& c, double T, Diagnostics& diag);  // J/mol
std::optional<double> saturated_liquid_volume(const PureComponent& c, double T, Diagnostics& diag);   // m³/mol

// Transport properties; same error and warning conventions. Gas values are for low pressure (ideal gas).
std::optional<double> liquid_viscosity(const PureComponent& c, double T, Diagnostics& diag);              // Pa s
std::optional<double> gas_viscosity(const PureComponent& c, double T, Diagnostics& diag);                 // Pa s
std::optional<double> liquid_thermal_conductivity(const PureComponent& c, double T, Diagnostics& diag);   // W/(m K)
std::optional<double> gas_thermal_conductivity(const PureComponent& c, double T, Diagnostics& diag);      // W/(m K)
std::optional<double> surface_tension(const PureComponent& c, double T, Diagnostics& diag);               // N/m
std::optional<double> liquid_cp(const PureComponent& c, double T, Diagnostics& diag);                     // J/(mol K)
// ∫ Cp,liquid dT from T0 to T1 (J/mol), exact for the DIPPR polynomial; warns outside the fitted range.
std::optional<double> liquid_enthalpy_change(const PureComponent& c, double T0, double T1, Diagnostics& diag);

}  // namespace chemlarp
