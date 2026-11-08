#include "omnichem/thermo.hpp"

#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <unordered_map>

#include "omnichem/joback.hpp"
#include "omnichem/smiles.hpp"

namespace omnichem {

double Nasa7::cp_R(double T) const noexcept {
  const auto& a = coeffs(T);
  return a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4])));
}

double Nasa7::h_RT(double T) const noexcept {
  const auto& a = coeffs(T);
  return a[0] + T * (a[1] / 2 + T * (a[2] / 3 + T * (a[3] / 4 + T * a[4] / 5))) + a[5] / T;
}

double Nasa7::s_R(double T) const noexcept {
  const auto& a = coeffs(T);
  return a[0] * std::log(T) + T * (a[1] + T * (a[2] / 2 + T * (a[3] / 3 + T * a[4] / 4))) + a[6];
}

double SpeciesThermo::h_RT(double T) const noexcept {
  if (source == ThermoSource::Nasa7Database) return nasa.h_RT(T);
  return h298_J / (kGasConstant * T);
}

double SpeciesThermo::s_R(double T) const noexcept {
  if (source == ThermoSource::Nasa7Database) return nasa.s_R(T);
  return s298_J / kGasConstant;
}

double SpeciesThermo::cp_R(double T) const noexcept {
  if (source == ThermoSource::Nasa7Database) return nasa.cp_R(T);
  return std::nan("");
}

bool SpeciesThermo::in_range(double T) const noexcept {
  if (source == ThermoSource::Nasa7Database) return T >= nasa.t_low && T <= nasa.t_high;
  return T >= 250.0 && T <= 400.0;  // constant-ΔH/S estimate is only sound near 298 K
}

namespace {

struct Database {
  std::unordered_map<std::string, SpeciesThermo> by_canonical;  // primary key: canonical SMILES
  std::unordered_map<std::string, std::string> canonical_by_formula;
};

// GRI-Mech 3.0 thermodynamic data (Smith et al., http://combustion.berkeley.edu/gri-mech/).
// Entries are declared with a SMILES; keys are that SMILES canonicalised by our own
// pipeline, so lookups match any equivalent spelling of the same structure.
const Database& database() {
  static const Database db = [] {
    Database m;
    auto add = [&](const char* smiles, const char* label, double tl, double tm, double th,
                   std::array<double, 7> lo, std::array<double, 7> hi) {
      SpeciesThermo t;
      t.source = ThermoSource::Nasa7Database;
      t.nasa = Nasa7{tl, tm, th, lo, hi};
      t.label = label;
      ChemicalSystemState scratch;
      Diagnostics d;
      const std::int32_t sp = parse_smiles(scratch, smiles, d);
      if (sp < 0) return;  // unreachable for the literals below
      const Species& s = scratch.species[static_cast<std::size_t>(sp)];
      m.by_canonical.emplace(s.canonical_smiles, t);
      m.canonical_by_formula.emplace(s.formula, s.canonical_smiles);
    };
    add("[H][H]", "H2 (GRI-Mech 3.0)", 200, 1000, 3500,
        {2.34433112E+00, 7.98052075E-03, -1.94781510E-05, 2.01572094E-08, -7.37611761E-12, -9.17935173E+02, 6.83010238E-01},
        {3.33727920E+00, -4.94024731E-05, 4.99456778E-07, -1.79566394E-10, 2.00255376E-14, -9.50158922E+02, -3.20502331E+00});
    add("O=O", "O2 (GRI-Mech 3.0)", 200, 1000, 3500,
        {3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00},
        {3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00});
    add("O", "H2O (GRI-Mech 3.0)", 200, 1000, 3500,
        {4.19864056E+00, -2.03643410E-03, 6.52040211E-06, -5.48797062E-09, 1.77197817E-12, -3.02937267E+04, -8.49032208E-01},
        {3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00});
    add("C", "CH4 (GRI-Mech 3.0)", 200, 1000, 3500,
        {5.14987613E+00, -1.36709788E-02, 4.91800599E-05, -4.84743026E-08, 1.66693956E-11, -1.02466476E+04, -4.64130376E+00},
        {7.48514950E-02, 1.33909467E-02, -5.73285809E-06, 1.22292535E-09, -1.01815230E-13, -9.46834459E+03, 1.84373180E+01});
    add("O=C=O", "CO2 (GRI-Mech 3.0)", 200, 1000, 3500,
        {2.35677352E+00, 8.98459677E-03, -7.12356269E-06, 2.45919022E-09, -1.43699548E-13, -4.83719697E+04, 9.90105222E+00},
        {3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00});
    add("[C-]#[O+]", "CO (GRI-Mech 3.0)", 200, 1000, 3500,
        {3.57953347E+00, -6.10353680E-04, 1.01681433E-06, 9.07005884E-10, -9.04424499E-13, -1.43440860E+04, 3.50840928E+00},
        {2.71518561E+00, 2.06252743E-03, -9.98825771E-07, 2.30053008E-10, -2.03647716E-14, -1.41518724E+04, 7.81868772E+00});
    add("N#N", "N2 (GRI-Mech 3.0)", 300, 1000, 5000,
        {3.298677E+00, 1.4082404E-03, -3.963222E-06, 5.641515E-09, -2.444854E-12, -1.0208999E+03, 3.950372E+00},
        {2.92664E+00, 1.4879768E-03, -5.68476E-07, 1.0097038E-10, -6.753351E-15, -9.227977E+02, 5.980528E+00});
    add("[H]", "H atom (GRI-Mech 3.0)", 200, 1000, 3500,
        {2.50000000E+00, 7.05332819E-13, -1.99591964E-15, 2.30081632E-18, -9.27732332E-22, 2.54736599E+04, -4.46682853E-01},
        {2.50000001E+00, -2.30842973E-11, 1.61561948E-14, -4.73515235E-18, 4.98197357E-22, 2.54736599E+04, -4.46682914E-01});
    add("[O]", "O atom (GRI-Mech 3.0)", 200, 1000, 3500,
        {3.16826710E+00, -3.27931884E-03, 6.64306396E-06, -6.12806624E-09, 2.11265971E-12, 2.91222592E+04, 2.05193346E+00},
        {2.56942078E+00, -8.59741137E-05, 4.19484589E-08, -1.00177799E-11, 1.22833691E-15, 2.92175791E+04, 4.78433864E+00});
    add("[O]O", "HO2 (GRI-Mech 3.0)", 200, 1000, 3500,
        {4.30179801E+00, -4.74912051E-03, 2.11582891E-05, -2.42763894E-08, 9.29225124E-12, 2.94808040E+02, 3.71666245E+00},
        {4.01721090E+00, 2.23982013E-03, -6.33658150E-07, 1.14246370E-10, -1.07908535E-14, 1.11856713E+02, 3.78510215E+00});
    add("OO", "H2O2 (GRI-Mech 3.0)", 200, 1000, 3500,
        {4.27611269E+00, -5.42822417E-04, 1.67335701E-05, -2.15770813E-08, 8.62454363E-12, -1.77025821E+04, 3.43505074E+00},
        {4.16500285E+00, 4.90831694E-03, -1.90139225E-06, 3.71185986E-10, -2.87908305E-14, -1.78617877E+04, 2.91615662E+00});
    add("[Ar]", "Ar (GRI-Mech 3.0)", 300, 1000, 5000,
        {2.5, 0.0, 0.0, 0.0, 0.0, -7.45375E+02, 4.366},
        {2.5, 0.0, 0.0, 0.0, 0.0, -7.45375E+02, 4.366});
    add("[OH]", "OH radical (GRI-Mech 3.0)", 200, 1000, 3500,
        {3.99201543E+00, -2.40131752E-03, 4.61793841E-06, -3.88113333E-09, 1.36411470E-12, 3.61508056E+03, -1.03925458E-01},
        {3.09288767E+00, 5.48429716E-04, 1.26505228E-07, -8.79461556E-11, 1.17412376E-14, 3.85865700E+03, 4.47669610E+00});
    return m;
  }();
  return db;
}

// Absolute S°(298.15 K) of each element's reference state, per atom, J/(mol K).
// (H2, N2, O2, F2, Cl2 gases; C graphite; S rhombic; Br2 liquid; I2 solid.) Source: CODATA key values.
double element_reference_entropy(std::size_t z) {
  switch (z) {
    case 1: return 130.680 / 2;
    case 6: return 5.74;
    case 7: return 191.609 / 2;
    case 8: return 205.152 / 2;
    case 9: return 202.791 / 2;
    case 16: return 32.054;
    case 17: return 223.081 / 2;
    case 35: return 152.21 / 2;
    case 53: return 116.14 / 2;
    default: return std::nan("");
  }
}

}  // namespace

std::optional<SpeciesThermo> builtin_thermo(const std::string& hill_formula) {
  const auto& db = database();
  auto key = db.canonical_by_formula.find(hill_formula);
  if (key == db.canonical_by_formula.end()) return std::nullopt;
  return db.by_canonical.at(key->second);
}

std::optional<SpeciesThermo> builtin_thermo_canonical(const std::string& canonical_smiles) {
  const auto& db = database();
  auto it = db.by_canonical.find(canonical_smiles);
  if (it == db.by_canonical.end()) return std::nullopt;
  return it->second;
}

std::optional<SpeciesThermo> resolve_thermo(const ChemicalSystemState& state, std::uint32_t species,
                                            double temperature_K, Diagnostics& diag) {
  const Species& s = state.species[species];
  std::optional<SpeciesThermo> out;
  // Exact structural match first. A formula-only match is accepted with a
  // warning: it catches alternative valence-bond spellings of the same species
  // (e.g. [C]=O for carbon monoxide) but would be wrong for a true isomer.
  out = builtin_thermo_canonical(s.canonical_smiles);
  if (!out && s.net_charge == 0) {
    out = builtin_thermo(s.formula);
    if (out) {
      diag.warn(DiagCode::ThermoFormulaMatch,
                s.smiles + " matched database entry " + database().canonical_by_formula.at(s.formula) +
                    " by formula only; check that it is the same species")
          .on_species(species);
    }
  }

  if (!out) {
    Diagnostics local;
    JobackResult jb = joback_estimate(state, species, local);
    for (const auto& d : local.items()) {
      Diagnostic& copy = d.severity == Severity::Error ? diag.error(d.code, d.message) : diag.warn(d.code, d.message);
      copy.on_atom(d.atom).on_species(d.species);
    }
    if (!jb.ok) {
      diag.error(DiagCode::NoThermoData, "no thermodynamic data for " + s.formula + " (" + s.smiles +
                                             "): not in the built-in database and Joback groups are incomplete")
          .on_species(species);
      return std::nullopt;
    }
    double s_elements = 0.0;
    for (std::size_t z = 1; z < kElementSlots; ++z) {
      if (s.composition[z] == 0) continue;
      const double sz = element_reference_entropy(z);
      if (std::isnan(sz)) {
        diag.error(DiagCode::NoThermoData, "no element reference entropy for " +
                                               std::string(element_symbol(static_cast<std::uint8_t>(z))))
            .on_species(species);
        return std::nullopt;
      }
      s_elements += s.composition[z] * sz;
    }
    SpeciesThermo t;
    t.source = ThermoSource::JobackEstimate;
    t.h298_J = jb.dHf_kJ * 1000.0;
    const double dS_formation = (jb.dHf_kJ - jb.dGf_kJ) * 1000.0 / kTref;  // ΔfS° = (ΔfH° - ΔfG°)/T
    t.s298_J = dS_formation + s_elements;
    t.label = "Joback estimate";
    diag.warn(DiagCode::ThermoApproximation,
              s.formula + " (" + s.smiles + "): Joback group-additivity estimate, ΔH° and S° held at 298 K values")
        .on_species(species);
    out = t;
  }
  if (!out->in_range(temperature_K)) {
    diag.warn(DiagCode::TemperatureOutOfRange,
              s.formula + ": T = " + std::to_string(temperature_K) + " K is outside the data's valid range; extrapolating")
        .on_species(species);
  }
  return out;
}

bool bind_thermo(ChemicalSystemState& state, Diagnostics& diag) {
  const double T = state.conditions.temperature_K;
  if (!(T > 0.0) || !std::isfinite(T)) {
    diag.error(DiagCode::InvalidInput, "temperature must be a positive finite number of kelvin");
    return false;
  }
  bool ok = true;
  for (std::uint32_t i = 0; i < state.species.size(); ++i) {
    auto t = resolve_thermo(state, i, T, diag);
    if (!t) { ok = false; continue; }
    state.g_RT[i] = t->g_RT(T);
  }
  ++state.revision;
  return ok;
}

}  // namespace omnichem
