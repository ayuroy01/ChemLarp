// Aqueous speciation, pH and solubility.
//
// Reference values are from PHREEQC 3 (IPhreeqc through phreeqpython, phreeqc.dat, solutions defined with the
// pH from the charge balance) unless a measurement is named; scripts/validate_aqueous.py compares many more
// cases. Known differences, documented in docs/ARCHITECTURE.md: PHREEQC adds a correction to the CO2(g) and
// CH4(g) solubility constants, and uses another constant for barite, than phreeqc.dat lists.
#include <cmath>
#include <cstdlib>
#include <string>

#include "chemlarp/aqueous.hpp"
#include "chemlarp/chemlarp_c.h"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

AqueousResult solve(const std::string& text) {
  Diagnostics d;
  auto in = parse_aqueous(text, d);
  CHECK(in.has_value());
  auto r = aqueous_equilibrium(*in, d);
  CHECK(r.has_value());
  return *r;
}

double molality(const AqueousResult& r, const std::string& name) {
  for (const auto& s : r.species)
    if (s.name == name) return s.molality;
  return 0.0;
}

double si(const AqueousResult& r, const std::string& phase) {
  for (const auto& s : r.saturation)
    if (s.phase == phase) return s.si;
  throw ::cltest::Failure("no saturation index for " + phase);
}

double total(const AqueousResult& r, const std::string& c) {
  for (const auto& [k, v] : r.totals)
    if (k == c) return v;
  return 0.0;
}

}  // namespace

TEST(aqueous_pure_water_against_phreeqc_and_iupac) {
  CHECK_NEAR(solve("temperature 0\n").pH, 7.46923, 1e-4);
  const AqueousResult w25 = solve("");
  CHECK_NEAR(w25.pH, 6.99738, 1e-4);
  CHECK_NEAR(2 * w25.pH, 13.995, 0.01);  // IUPAC pKw at 25 °C
  CHECK_NEAR(solve("temperature 100\n").pH, 6.11916, 1e-4);
  CHECK_NEAR(w25.A_dh, 0.5108, 0.002);  // Debye–Hückel A at 25 °C (0.509–0.511 in the literature)
}

TEST(aqueous_carbonate_solution_against_phreeqc) {
  // 1 mmol/kg CaCl2 + 2 mmol/kg NaHCO3.
  const AqueousResult r = solve("add CaCl2 1\nadd NaHCO3 2\n");
  CHECK_NEAR(r.pH, 8.12277, 2e-4);
  CHECK_NEAR(r.ionic_strength, 0.00493685, 1e-7);
  CHECK_NEAR(molality(r, "CO3-2"), 1.48984e-5, 0.002 * 1.48984e-5);
  CHECK_NEAR(molality(r, "CaHCO3+"), 1.76845e-5, 0.002 * 1.76845e-5);
  CHECK_NEAR(si(r, "Calcite"), 0.384066, 0.002);
  CHECK_NEAR(r.alkalinity_eq, 0.002, 1e-9);  // bicarbonate alkalinity is conserved
  CHECK_NEAR(r.charge_balance_eq, 0.0, 1e-14);
  CHECK_NEAR(total(r, "Ca"), 0.001, 1e-15);
  CHECK_NEAR(total(r, "C(4)"), 0.002, 1e-15);
  // Mass balance over the species.
  double c = 0.0;
  for (const char* n : {"HCO3-", "CO3-2", "CO2", "CaHCO3+", "CaCO3", "NaHCO3", "NaCO3-"}) c += molality(r, n);
  c += 2 * molality(r, "(CO2)2");
  CHECK_NEAR(c, 0.002, 1e-12);
}

TEST(aqueous_acids_bases_and_hydrolysis_against_phreeqc) {
  CHECK_NEAR(solve("add HCl 1\n").pH, 3.01475, 2e-4);
  CHECK_NEAR(solve("add NH4Cl 10\n").pH, 5.66516, 1e-3);
  CHECK_NEAR(solve("N(-3) 10\n").pH, 10.6102, 2e-4);  // ammonia: the charge balance makes it a base
  CHECK_NEAR(solve("add NaH2PO4 10\nadd Na2HPO4 10\n").pH, 6.96416, 1e-3);
  const AqueousResult al = solve("add AlCl3 0.1\n");
  CHECK_NEAR(al.pH, 4.51541, 2e-4);
  CHECK_NEAR(molality(al, "AlOH+2"), 2.11814e-5, 0.002 * 2.11814e-5);
  CHECK_NEAR(solve("add FeCl3 0.1\n").pH, 3.79312, 2e-4);
  // Fixed pH: the charge balance is then reported, not enforced.
  const AqueousResult f = solve("pH 4\nadd NaCl 1\n");
  CHECK_NEAR(f.pH, 4.0, 1e-12);
  CHECK(f.charge_balance_eq > 0.0);
}

TEST(aqueous_minerals_and_gases) {
  // Calcite in a closed system: PHREEQC pH 9.9068.
  CHECK_NEAR(solve("equilibrium Calcite 0\n").pH, 9.90681, 2e-3);
  // Rain water with atmospheric CO2 (10^-3.5 atm): pH 5.6–5.7 (textbook 5.6).
  CHECK_NEAR(solve("equilibrium CO2(g) -3.5\n").pH, 5.66, 0.02);
  // Calcite with atmospheric CO2: PHREEQC pH 8.2748 and 0.496 mmol/kg Ca with its own CO2 correction; the
  // database constant (Plummer & Busenberg 1982) gives 8.279 and 0.493.
  const AqueousResult cc = solve("equilibrium Calcite 0\nequilibrium CO2(g) -3.5\n");
  CHECK_NEAR(cc.pH, 8.2748, 0.01);
  CHECK_NEAR(total(cc, "Ca"), 0.000496, 0.01 * 0.000496);
  CHECK_NEAR(si(cc, "Calcite"), 0.0, 1e-9);
  CHECK_NEAR(si(cc, "CO2(g)"), -3.5, 1e-9);
  // Gypsum: PHREEQC 15.085 mmol/kg; measured 15.45 (Marshall & Slusher 1966).
  const AqueousResult g = solve("equilibrium Gypsum 0\n");
  CHECK_NEAR(g.phases[0].dissolved_mol, 0.015085, 0.002 * 0.015085);
  CHECK_NEAR(g.phases[0].dissolved_mol, 0.01545, 0.03 * 0.01545);
  // Fluorite and H2S(g), NH3(g).
  CHECK_NEAR(solve("equilibrium Fluorite 0\n").phases[0].dissolved_mol, 1.94944e-4, 0.002 * 1.94944e-4);
  CHECK_NEAR(solve("equilibrium H2S(g) -3\n").pH, 5.49596, 2e-4);
  CHECK_NEAR(solve("equilibrium NH3(g) -5\n").pH, 10.0176, 2e-4);
  // A limited amount dissolves completely when the solution stays undersaturated.
  const AqueousResult lim = solve("equilibrium Gypsum 0 0.001\n");
  CHECK_NEAR(lim.phases[0].dissolved_mol, 0.001, 1e-15);
  CHECK_NEAR(lim.phases[0].remaining_mol, 0.0, 1e-15);
  CHECK(si(lim, "Gypsum") < 0.0);
  // Supersaturated solution precipitates.
  const AqueousResult pr = solve("add CaCl2 20\nadd Na2SO4 20\nequilibrium Gypsum 0 0\n");
  CHECK(pr.phases[0].dissolved_mol < 0.0);
  CHECK_NEAR(si(pr, "Gypsum"), 0.0, 1e-9);
}

TEST(aqueous_seawater_against_phreeqc) {
  const AqueousResult r =
      solve("units mmol/kgw\npH 8.22 charge\nCa 10.6\nMg 54.7\nNa 485.8\nK 10.6\nCl 565.8\nS(6) 29.3\nC(4) 2.4\n");
  CHECK_NEAR(r.pH, 8.14231, 0.005);
  CHECK_NEAR(r.ionic_strength, 0.674326, 0.001);
  CHECK_NEAR(molality(r, "MgSO4"), 0.00714145, 0.01 * 0.00714145);
  CHECK_NEAR(si(r, "Calcite"), 0.72218, 0.01);
}

TEST(aqueous_temperature_dependence_against_phreeqc) {
  const AqueousResult r = solve("temperature 60\nadd CaCl2 1\nadd NaHCO3 2\n");
  CHECK_NEAR(r.T_K, 333.15, 1e-12);
  CHECK_NEAR(molality(r, "HCO3-"), 0.00188698, 0.002 * 0.00188698);
  CHECK_NEAR(si(solve("temperature 25\nadd CaCl2 1\nadd NaHCO3 2\n"), "Calcite") < si(r, "Calcite"), 1.0, 0.0);
}

TEST(aqueous_formulas_and_input_errors) {
  Diagnostics d;
  auto in = parse_aqueous("units mol/kgw\nadd Fe2(SO4)3 0.001\nadd NH4NO3 0.002\nadd CaSO4:2H2O 0.001\n", d);
  CHECK(in.has_value());
  auto get = [&](const std::string& c) {
    for (const auto& [k, v] : in->totals)
      if (k == c) return v;
    return 0.0;
  };
  CHECK_NEAR(get("Fe(3)"), 0.002, 1e-15);
  CHECK_NEAR(get("S(6)"), 0.004, 1e-15);
  CHECK_NEAR(get("N(5)"), 0.002, 1e-15);
  CHECK_NEAR(get("N(-3)"), 0.002, 1e-15);
  CHECK_NEAR(get("Ca"), 0.001, 1e-15);
  for (const char* bad : {"Unobtainium 1\n", "equilibrium Kryptonite 0\n", "temperature 150\n", "add Fe3O4 1\n",
                          "add NaClO4 1\n", "add 1\n", "Ca -1\n", "units ppm\n"}) {
    Diagnostics e;
    CHECK(!parse_aqueous(bad, e).has_value());
    CHECK(e.has_errors());
  }
}

TEST(aqueous_c_api_json) {
  char* out = cl_aqueous_json("add NaCl 10\nequilibrium Calcite 0\n");
  CHECK(out != nullptr);
  const std::string s(out);
  cl_free_string(out);
  CHECK(s.rfind("{\"ok\":true", 0) == 0);
  CHECK(s.find("\"species\":[") != std::string::npos);
  CHECK(s.find("\"phase\":\"Calcite\"") != std::string::npos);
  char* bad = cl_aqueous_json("Xx 1\n");
  CHECK(std::string(bad).find("\"ok\":false") != std::string::npos);
  cl_free_string(bad);
}
