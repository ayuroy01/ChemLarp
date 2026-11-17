// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
#include "chemlarp/lab.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <set>
#include <sstream>

#include "chemlarp/aqueous.hpp"
#include "chemlarp/elements.hpp"
#include "chemlarp/json.hpp"

namespace chemlarp {

namespace {

// ---- GHS table (generated) ------------------------------------------------------------------------------
struct GhsRow {
  const char* smiles;
  const char* name;
  int cid;
  const char* signal;
  const char* pictograms;
  const char* hazards;  // "H314:100;H290:20.5"
};
const GhsRow kGhs[] = {
#include "ghs_data.inc"
};
struct GhsPhraseRow {
  const char* code;
  const char* text;
};
const GhsPhraseRow kGhsPhrases[] = {
#include "ghs_phrases.inc"
};

GhsEntry entry_of(const GhsRow& r) {
  GhsEntry e;
  e.smiles = r.smiles;
  e.name = r.name;
  e.cid = r.cid;
  e.signal = r.signal;
  std::istringstream p(r.pictograms);
  for (std::string w; p >> w;) e.pictograms.push_back(w);
  std::string h = r.hazards;
  std::size_t pos = 0;
  while (pos < h.size()) {
    std::size_t end = h.find(';', pos);
    if (end == std::string::npos) end = h.size();
    const std::string item = h.substr(pos, end - pos);
    const std::size_t colon = item.find(':');
    if (colon != std::string::npos) e.hazards.push_back({item.substr(0, colon), std::atof(item.c_str() + colon + 1)});
    pos = end + 1;
  }
  return e;
}

std::optional<GhsEntry> ghs_by_name(std::string_view name) {
  for (const auto& r : kGhs)
    if (name == r.name) return entry_of(r);
  return std::nullopt;
}

// ---- reagent catalogue ----------------------------------------------------------------------------------
struct ReagentRow {
  const char* id;
  const char* name;
  const char* formula;   // aqueous model; "" = not modelled
  const char* smiles;
  const char* groups;    // space-separated
  const char* ghs_name;  // name in the GHS table
};
// Groups: acid base oxidizer strong_oxidizer reducing organic alcohol halocarbon alkali_metal water_reactive
// cyanide ammonia hypochlorite iodine peroxide_inorganic alkali_peroxide permanganate sulfuric_acid nitric_acid.
// Oxidizers follow the examples of Prudent Practices Table 4.6; "strong" those named in §4.D.1.4.1 (nitric acid,
// chlorates, permanganates, peroxides, persulfates) and perchloric acid (§6.G.1).
const ReagentRow kReagents[] = {
    {"hcl", "Hydrochloric acid", "HCl", "Cl", "acid", "hydrochloric acid"},
    {"h2so4", "Sulfuric acid", "H2SO4", "OS(=O)(=O)O", "acid oxidizer sulfuric_acid", "sulfuric acid"},
    {"hno3", "Nitric acid", "HNO3", "O[N+](=O)[O-]", "acid oxidizer strong_oxidizer nitric_acid", "nitric acid"},
    {"h3po4", "Phosphoric acid", "H3PO4", "OP(=O)(O)O", "acid", "phosphoric acid"},
    {"hf", "Hydrofluoric acid", "HF", "F", "acid", "hydrofluoric acid"},
    {"naoh", "Sodium hydroxide", "NaOH", "[Na+].[OH-]", "base", "sodium hydroxide"},
    {"koh", "Potassium hydroxide", "KOH", "[K+].[OH-]", "base", "potassium hydroxide"},
    {"nh3", "Ammonia (aqueous)", "NH3", "N", "base ammonia", "ammonia"},
    {"nh4cl", "Ammonium chloride", "NH4Cl", "[NH4+].[Cl-]", "", "ammonium chloride"},
    {"nh4no3", "Ammonium nitrate", "NH4NO3", "[NH4+].[O-][N+](=O)[O-]", "oxidizer", "ammonium nitrate"},
    {"na2co3", "Sodium carbonate", "Na2CO3", "C(=O)([O-])[O-].[Na+].[Na+]", "base", "sodium carbonate"},
    {"nahco3", "Sodium bicarbonate", "NaHCO3", "C(=O)(O)[O-].[Na+]", "", "sodium bicarbonate"},
    {"nacl", "Sodium chloride", "NaCl", "[Na+].[Cl-]", "", "sodium chloride"},
    {"kcl", "Potassium chloride", "KCl", "[K+].[Cl-]", "", "potassium chloride"},
    {"kbr", "Potassium bromide", "KBr", "[K+].[Br-]", "", "potassium bromide"},
    {"naf", "Sodium fluoride", "NaF", "[Na+].[F-]", "", ""},
    {"cacl2", "Calcium chloride", "CaCl2", "[Ca+2].[Cl-].[Cl-]", "", "calcium chloride"},
    {"mgso4", "Magnesium sulfate", "MgSO4", "[O-]S(=O)(=O)[O-].[Mg+2]", "", "magnesium sulfate"},
    {"na2so4", "Sodium sulfate", "Na2SO4", "[O-]S(=O)(=O)[O-].[Na+].[Na+]", "", "sodium sulfate"},
    {"bacl2", "Barium chloride", "BaCl2", "[Ba+2].[Cl-].[Cl-]", "", "barium chloride"},
    {"cuso4", "Copper(II) sulfate", "CuSO4", "[O-]S(=O)(=O)[O-].[Cu+2]", "", "copper(II) sulfate"},
    {"feso4", "Iron(II) sulfate", "FeSO4", "[O-]S(=O)(=O)[O-].[Fe+2]", "", "iron(II) sulfate"},
    {"fecl3", "Iron(III) chloride", "FeCl3", "[Cl-].[Cl-].[Cl-].[Fe+3]", "", "iron(III) chloride"},
    {"zncl2", "Zinc chloride", "ZnCl2", "[Cl-].[Cl-].[Zn+2]", "", "zinc chloride"},
    {"pbno3", "Lead(II) nitrate", "Pb(NO3)2", "[O-][N+](=O)[O-].[O-][N+](=O)[O-].[Pb+2]", "oxidizer", "lead(II) nitrate"},
    {"na2s", "Sodium sulfide", "Na2S", "[Na+].[Na+].[S-2]", "", "sodium sulfide"},
    // Not in the solution model: safety rules and GHS only.
    {"agno3", "Silver nitrate", "", "[Ag+].[O-][N+](=O)[O-]", "oxidizer", "silver nitrate"},
    {"naocl", "Sodium hypochlorite (bleach)", "", "[Na+].[O-]Cl", "oxidizer hypochlorite", "sodium hypochlorite"},
    {"h2o2", "Hydrogen peroxide", "", "OO", "oxidizer strong_oxidizer", "hydrogen peroxide"},
    {"kmno4", "Potassium permanganate", "", "[K+].[O-][Mn](=O)(=O)=O", "oxidizer strong_oxidizer permanganate", "potassium permanganate"},
    {"k2cr2o7", "Potassium dichromate", "", "[K+].[K+].[O-][Cr](=O)(=O)O[Cr](=O)(=O)[O-]", "oxidizer", "potassium dichromate"},
    {"kclo3", "Potassium chlorate", "", "[K+].[O-]Cl(=O)=O", "oxidizer strong_oxidizer", "potassium chlorate"},
    {"na2o2", "Sodium peroxide", "", "[Na+].[Na+].[O-][O-]", "oxidizer strong_oxidizer peroxide_inorganic alkali_peroxide", "sodium peroxide"},
    {"kcn", "Potassium cyanide", "", "[K+].[C-]#N", "cyanide", "potassium cyanide"},
    {"nacn", "Sodium cyanide", "", "[Na+].[C-]#N", "cyanide", "sodium cyanide"},
    {"na", "Sodium (metal)", "", "[Na]", "alkali_metal water_reactive reducing", "sodium"},
    {"k", "Potassium (metal)", "", "[K]", "alkali_metal water_reactive reducing", "potassium"},
    {"li", "Lithium (metal)", "", "[Li]", "alkali_metal water_reactive reducing", "lithium"},
    {"nabh4", "Sodium borohydride", "", "[BH4-].[Na+]", "water_reactive reducing", "sodium borohydride"},
    {"lialh4", "Lithium aluminium hydride", "", "[AlH4-].[Li+]", "water_reactive reducing", "lithium aluminium hydride"},
    {"cao", "Calcium oxide", "", "[Ca+2].[O-2]", "water_reactive", "calcium oxide"},
    {"i2", "Iodine", "", "II", "iodine oxidizer", "iodine"},
    {"acetic", "Acetic acid", "", "CC(=O)O", "acid organic", "acetic acid"},
    {"ethanol", "Ethanol", "", "CCO", "organic alcohol", "ethanol"},
    {"methanol", "Methanol", "", "CO", "organic alcohol", "methanol"},
    {"ipa", "2-Propanol", "", "CC(C)O", "organic alcohol", "2-propanol"},
    {"acetone", "Acetone", "", "CC(C)=O", "organic", "acetone"},
    {"hexane", "Hexane", "", "CCCCCC", "organic", "hexane"},
    {"toluene", "Toluene", "", "Cc1ccccc1", "organic", "toluene"},
    {"dcm", "Dichloromethane", "", "ClCCl", "organic halocarbon", "dichloromethane"},
    {"chloroform", "Chloroform", "", "ClC(Cl)Cl", "organic halocarbon", "chloroform"},
    {"glucose", "Glucose", "", "OCC1OC(O)C(O)C(O)C1O", "organic", "glucose"},
    {"sucrose", "Sucrose", "", "OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O", "organic", "sucrose"},
};

const ReagentRow* reagent_row(std::string_view id) {
  for (const auto& r : kReagents)
    if (id == r.id) return &r;
  return nullptr;
}

std::vector<std::string> split_words(const char* s) {
  std::istringstream in(s);
  std::vector<std::string> out;
  for (std::string w; in >> w;) out.push_back(w);
  return out;
}

// ---- incompatibility rules (NRC, Prudent Practices in the Laboratory, 2011) ------------------------------
struct Rule {
  const char* a;
  const char* b;  // "water": any aqueous portion
  const char* severity;
  const char* message;
  const char* source;
};
const char* kBook = "National Research Council, Prudent Practices in the Laboratory (updated version, 2011), ";
const Rule kRules[] = {
    {"cyanide", "acid", "danger", "Cyanide with acid forms hydrogen cyanide (HCN), a highly toxic gas.", "section 4.D.2.3"},
    {"oxidizer", "reducing", "danger",
     "Oxidizing agents are incompatible with reducing agents (above all when concentrated): fire or explosion.",
     "sections 4.D.2.3 and 5.E.7"},
    {"strong_oxidizer", "organic", "danger",
     "Strong oxidizing agents with organic material can ignite spontaneously or react violently.",
     "sections 4.D.1.4.1 and 6.G.1"},
    {"permanganate", "sulfuric_acid", "danger",
     "Potassium permanganate with sulfuric acid has caused explosions (mixed cleaning baths).", "section 4.D.3.3"},
    {"nitric_acid", "alcohol", "danger", "Nitric acid with alcohols has caused explosions.", "section 4.D.3.3"},
    {"alkali_metal", "halocarbon", "danger", "Alkali metals are extremely reactive toward halogenated solvents.",
     "section 4.D.2.3"},
    {"water_reactive", "water", "danger", "Reacts violently with water, giving heat and flammable hydrogen gas.",
     "sections 4.D.2.1 and 6.G.6"},
    {"ammonia", "hypochlorite", "danger",
     "Ammonia reacts with hypochlorite (bleach), giving toxic gases (chlorine; chloramine is among the gases that "
     "damage the lungs).",
     "sections 6.G.6 and 4.C.3.1.2"},
    {"ammonia", "iodine", "danger", "Ammonia with iodine forms nitrogen triiodide, which explodes on touch.",
     "section 6.G.6"},
    {"peroxide_inorganic", "organic", "warning",
     "Inorganic peroxides can form organic peroxides and hydroperoxides with organic compounds.", "section 6.G.1"},
    {"alkali_peroxide", "water", "danger", "Alkali-metal peroxides react violently with water.", "section 6.G.1"},
};

// ---- formula masses -------------------------------------------------------------------------------------
// Molar mass of a mineral formula such as "CaSO4:2H2O", "Ca5(PO4)3OH", "Mg2Si3O7.5OH:3H2O"; 0 if unreadable.
double formula_mass(const std::string& f) {
  double total = 0.0;
  std::size_t start = 0;
  while (start <= f.size()) {
    std::size_t end = f.find(':', start);
    if (end == std::string::npos) end = f.size();
    std::string part = f.substr(start, end - start);
    double mult = 1.0;
    std::size_t i = 0;
    while (i < part.size() && (std::isdigit(static_cast<unsigned char>(part[i])) || part[i] == '.')) ++i;
    if (i > 0) mult = std::atof(part.substr(0, i).c_str()), part = part.substr(i);
    // Recursive descent over groups.
    std::size_t p = 0;
    std::function<double()> group = [&]() -> double {
      double m = 0.0;
      while (p < part.size() && part[p] != ')') {
        double unit = 0.0;
        if (part[p] == '(') {
          ++p;
          unit = group();
          if (p < part.size() && part[p] == ')') ++p;
        } else if (std::isupper(static_cast<unsigned char>(part[p]))) {
          std::string sym(1, part[p++]);
          while (p < part.size() && std::islower(static_cast<unsigned char>(part[p]))) sym += part[p++];
          const ElementInfo* e = find_element_by_symbol(sym);
          if (!e) return std::nan("");
          unit = e->mass;
        } else {
          return std::nan("");
        }
        std::size_t q = p;
        while (q < part.size() && (std::isdigit(static_cast<unsigned char>(part[q])) || part[q] == '.')) ++q;
        const double n = q > p ? std::atof(part.substr(p, q - p).c_str()) : 1.0;
        p = q;
        m += unit * n;
      }
      return m;
    };
    const double m = group();
    if (!std::isfinite(m)) return 0.0;
    total += mult * m;
    start = end + 1;
  }
  return total;
}

// Phases that form within a laboratory experiment (precipitates) and gases that leave an open vessel.
const std::set<std::string> kFastPhases = {
    "Calcite", "Witherite", "Strontianite", "Gypsum", "Celestite", "Barite", "Fluorite", "Al(OH)3(a)", "Fe(OH)3(a)",
    "FeS(ppt)", "Pyrochroite", "Zn(OH)2(e)", "Smithsonite", "Sphalerite", "Cd(OH)2", "Otavite", "Cerussite",
    "Anglesite", "Pb(OH)2", "Halite", "Sylvite", "Vivianite", "Melanterite"};
const std::set<std::string> kGases = {"CO2(g)", "H2S(g)", "NH3(g)"};
const std::map<std::string, std::string> kGasGhsName = {
    {"CO2(g)", "carbon dioxide"}, {"H2S(g)", "hydrogen sulfide"}, {"NH3(g)", "ammonia"}};

constexpr double kWaterDensity = 0.99705;  // kg/L at 25 °C
constexpr double kCp = 4.184;              // J g⁻¹ K⁻¹
constexpr double kMolarVolume = 24465.0;   // mL/mol, ideal gas at 25 °C, 1 atm

std::string num(double v, int prec = 6) {
  char b[48];
  std::snprintf(b, sizeof b, "%.*g", prec, v);
  return b;
}

struct Solved {
  AqueousResult r;
  std::vector<std::string> phases;
};

// Solves a solution of the given moles (formula -> mol) in kgw of water at T; with `react`, supersaturated
// fast phases precipitate and gases above 1 atm leave (iterated until the set is stable).
std::optional<Solved> solve(const std::vector<std::pair<std::string, double>>& moles, double kgw, double T_C, bool react,
                            const std::vector<std::string>& start_phases, Diagnostics& diag) {
  std::vector<std::string> phases = start_phases;
  for (int round = 0; round < 8; ++round) {
    std::string text = "temperature " + num(T_C, 10) + "\nunits mol/kgw\n";
    for (const auto& [f, n] : moles)
      if (n > 0) text += "add " + f + " " + num(n / kgw, 12) + "\n";
    for (const auto& p : phases) text += "equilibrium " + p + " 0 0\n";
    Diagnostics d;
    auto in = parse_aqueous(text, d);
    std::optional<AqueousResult> r;
    if (in) r = aqueous_equilibrium(*in, d);
    if (!r) {
      for (const auto& x : d.items()) (x.severity == Severity::Error ? diag.error(x.code, x.message) : diag.warn(x.code, x.message));
      return std::nullopt;
    }
    if (!react) return Solved{*r, phases};
    bool added = false;
    for (const auto& s : r->saturation) {
      const bool fast = kFastPhases.count(s.phase) > 0, gas = kGases.count(s.phase) > 0;
      if ((fast || gas) && s.si > 1e-6 && std::find(phases.begin(), phases.end(), s.phase) == phases.end()) {
        phases.push_back(s.phase);
        added = true;
      }
    }
    if (!added) return Solved{*r, phases};
  }
  diag.error(DiagCode::NotConverged, "the set of precipitates and gases did not settle");
  return std::nullopt;
}

double enthalpy_kJ(const AqueousResult& r, double kgw) {
  double h = r.enthalpy_kJ * kgw;
  for (const auto& p : r.phases) h += p.dissolved_mol * kgw * p.dissolution_enthalpy_kJ;
  return h;
}

LabHazard ghs_hazard(const std::string& who, const GhsEntry& g) {
  LabHazard h;
  h.kind = "ghs";
  h.severity = g.signal == "Danger" ? "danger" : g.signal == "Warning" ? "warning" : "info";
  // The classification is that of the pure substance (hydrogen chloride gas for hydrochloric acid, ammonia gas
  // for aqueous ammonia); a dilute solution is usually classified lower.
  h.message = who + " (GHS classification of the pure substance): " + (g.signal.empty() ? std::string("no signal word") : g.signal);
  for (const auto& [code, pct] : g.hazards) {
    const std::string phrase = ghs_phrase(code);
    h.message += "; " + code + (phrase.empty() ? "" : " " + phrase) + (pct < 100.0 ? " (" + num(pct, 3) + "% of notifications)" : "");
  }
  if (!g.pictograms.empty()) {
    h.message += "; pictograms";
    for (const auto& p : g.pictograms) h.message += " " + p;
  }
  h.source = "PubChem CID " + std::to_string(g.cid) + ", GHS classification (first source listed)";
  h.substances = {who};
  return h;
}

void computed_hazards(LabMixResult& out) {
  if (out.pH <= 2.0 || out.pH >= 11.5) {
    LabHazard h;
    h.kind = "computed";
    h.severity = "danger";
    h.message = "pH " + num(out.pH, 3) + ": corrosive to skin and eyes (GHS presumes skin corrosion at pH ≤ 2 or ≥ 11.5).";
    h.source = "UN GHS chapter 3.2 (skin corrosion/irritation), pH presumption; as a waste, US EPA 40 CFR 261.22 "
               "(pH ≤ 2 or ≥ 12.5)";
    out.hazards.push_back(h);
  }
  if (out.T_final_C >= 100.0) {
    LabHazard h;
    h.kind = "computed";
    h.severity = "danger";
    h.message = "The mixture would reach " + num(out.T_final_C, 3) +
                " °C and boil (spattering); the model assumes liquid water, so the result above 100 °C is not valid.";
    h.source = "energy balance of the mixture";
    out.hazards.push_back(h);
  }
  if (std::fabs(out.heat_released_kJ) > 1e-3) {
    LabHazard h;
    h.kind = "computed";
    h.severity = "info";
    h.message = std::string(out.heat_released_kJ > 0 ? "Releases " : "Absorbs ") + num(std::fabs(out.heat_released_kJ), 3) +
                " kJ: " + num(out.T_initial_C, 3) + " → " + num(out.T_final_C, 3) + " °C without heat loss.";
    h.source = "reaction enthalpies of the solution model (Hess's law)";
    out.hazards.push_back(h);
  }
  for (const auto& g : out.gases) {
    LabHazard h;
    h.kind = "computed";
    h.severity = "info";
    h.message = g.phase.substr(0, g.phase.size() - 3) + " is given off: " + num(g.mol, 3) + " mol, about " + num(g.mL, 3) +
                " mL at 25 °C and 1 atm.";
    h.source = "gas exceeds 1 atm in the solution model (open vessel)";
    auto it = kGasGhsName.find(g.phase);
    if (it != kGasGhsName.end())
      if (auto e = ghs_by_name(it->second)) {
        h.severity = e->signal == "Danger" ? "danger" : e->signal == "Warning" ? "warning" : "info";
        h.message += " Its GHS classification: " + (e->signal.empty() ? std::string("none listed") : e->signal);
        for (const auto& [code, pct] : e->hazards) h.message += " " + code;
        h.message += ".";
        h.source += "; PubChem CID " + std::to_string(e->cid);
      }
    out.hazards.push_back(h);
  }
}

// Collects modelled moles (formula -> mol) and the catalogue ids of a set of portions.
bool gather(const std::vector<LabPortion>& portions, std::vector<std::pair<std::string, double>>& moles,
            std::vector<std::string>& ids, std::vector<std::string>& not_modelled, double& volume_mL, Diagnostics& diag) {
  std::map<std::string, double> acc;
  volume_mL = 0.0;
  for (const auto& p : portions) {
    if (!(p.volume_mL > 0.0)) {
      diag.error(DiagCode::InvalidInput, "portion '" + p.label + "' needs a volume above 0 mL");
      return false;
    }
    volume_mL += p.volume_mL;
    for (const auto& c : p.solutes) {
      if (!(c.molarity >= 0.0)) {
        diag.error(DiagCode::InvalidInput, "'" + c.reagent + "' needs a concentration of 0 mol/L or more");
        return false;
      }
      std::string formula = c.reagent;
      if (const ReagentRow* r = reagent_row(c.reagent)) {
        if (std::find(ids.begin(), ids.end(), r->id) == ids.end()) ids.push_back(r->id);
        formula = r->formula;
        if (formula.empty()) {
          if (std::find(not_modelled.begin(), not_modelled.end(), r->name) == not_modelled.end()) not_modelled.push_back(r->name);
          continue;
        }
      }
      acc[formula] += c.molarity * p.volume_mL / 1000.0;
    }
  }
  moles.assign(acc.begin(), acc.end());
  return true;
}

}  // namespace

// ---- public ---------------------------------------------------------------------------------------------

const std::vector<LabReagent>& lab_reagents() {
  static const std::vector<LabReagent> list = [] {
    std::vector<LabReagent> v;
    for (const auto& r : kReagents) v.push_back({r.id, r.name, r.formula, r.smiles, split_words(r.groups)});
    return v;
  }();
  return list;
}

const LabReagent* find_lab_reagent(std::string_view id) {
  for (const auto& r : lab_reagents())
    if (r.id == id) return &r;
  return nullptr;
}

std::optional<GhsEntry> ghs_lookup(std::string_view smiles) {
  for (const auto& r : kGhs)
    if (smiles == r.smiles) return entry_of(r);
  return std::nullopt;
}

std::string ghs_phrase(std::string_view code) {
  for (const auto& p : kGhsPhrases)
    if (code == p.code) return p.text;
  return "";
}

std::size_t ghs_table_size() { return sizeof(kGhs) / sizeof(kGhs[0]); }

std::vector<LabHazard> lab_incompatibilities(const std::vector<std::string>& ids, bool aqueous) {
  std::vector<LabHazard> out;
  auto members = [&](const char* group) {
    std::vector<std::string> who;
    if (std::string(group) == "water") {
      if (aqueous) who.push_back("water");
      return who;
    }
    for (const auto& id : ids)
      if (const ReagentRow* r = reagent_row(id)) {
        const auto g = split_words(r->groups);
        if (std::find(g.begin(), g.end(), group) != g.end()) who.push_back(r->name);
      }
    return who;
  };
  for (const auto& rule : kRules) {
    const auto a = members(rule.a), b = members(rule.b);
    if (a.empty() || b.empty()) continue;
    // A pair needs two different substances.
    std::vector<std::string> both;
    for (const auto& x : a) both.push_back(x);
    for (const auto& y : b)
      if (std::find(both.begin(), both.end(), y) == both.end()) both.push_back(y);
    if (both.size() < 2) continue;
    LabHazard h;
    h.kind = "incompatibility";
    h.severity = rule.severity;
    h.message = rule.message;
    h.source = std::string(kBook) + rule.source;
    h.substances = both;
    out.push_back(h);
  }
  return out;
}

std::optional<LabMixResult> lab_mix(const std::vector<LabPortion>& portions, Diagnostics& diag) {
  if (portions.empty()) {
    diag.error(DiagCode::InvalidInput, "nothing to mix");
    return std::nullopt;
  }
  LabMixResult out;
  std::vector<std::pair<std::string, double>> moles;
  std::vector<std::string> ids;
  if (!gather(portions, moles, ids, out.not_modelled, out.volume_mL, diag)) return std::nullopt;
  // Initial temperature: volume-weighted (equal heat capacity per volume).
  double tv = 0.0;
  for (const auto& p : portions) tv += p.T_C * p.volume_mL;
  out.T_initial_C = tv / out.volume_mL;
  // Enthalpy of each portion on its own (no reactions with the others).
  double h0 = 0.0;
  for (const auto& p : portions) {
    std::vector<std::pair<std::string, double>> m1;
    std::vector<std::string> i1, nm1;
    double v1 = 0.0;
    Diagnostics quiet;
    if (!gather({p}, m1, i1, nm1, v1, quiet)) return std::nullopt;
    const double kgw = p.volume_mL / 1000.0 * kWaterDensity;
    auto s = solve(m1, kgw, out.T_initial_C, false, {}, diag);
    if (!s) return std::nullopt;
    h0 += enthalpy_kJ(s->r, kgw);
  }
  const double kgw = out.volume_mL / 1000.0 * kWaterDensity;
  auto mixed = solve(moles, kgw, out.T_initial_C, true, {}, diag);
  if (!mixed) return std::nullopt;
  out.heat_released_kJ = -(enthalpy_kJ(mixed->r, kgw) - h0);
  out.T_final_C = out.T_initial_C + out.heat_released_kJ * 1000.0 / (kCp * out.volume_mL);
  // The mixture at its final temperature (the model covers 0–100 °C).
  const double T_solve = std::clamp(out.T_final_C, 0.0, 100.0);
  auto final_state = std::fabs(T_solve - out.T_initial_C) > 0.01 ? solve(moles, kgw, T_solve, true, mixed->phases, diag) : mixed;
  if (!final_state) return std::nullopt;
  const AqueousResult& r = final_state->r;
  out.pH = r.pH;
  out.ionic_strength = r.ionic_strength;
  for (const auto& p : r.phases) {
    if (p.dissolved_mol >= -1e-12) continue;
    const double n = -p.dissolved_mol * kgw;
    if (kGases.count(p.phase)) {
      out.gases.push_back({p.phase, n, n * kMolarVolume});
    } else {
      std::string formula;
      for (const auto& s : r.saturation)
        if (s.phase == p.phase) formula = s.formula;
      out.precipitates.push_back({p.phase, formula, n, n * formula_mass(formula)});
    }
  }
  for (const auto& s : r.species) {
    if (out.major_species.size() >= 10) break;
    out.major_species.push_back({s.name, s.molality});
  }
  computed_hazards(out);
  for (auto& h : lab_incompatibilities(ids, true)) out.hazards.push_back(std::move(h));
  for (const auto& id : ids) {
    const ReagentRow* row = reagent_row(id);
    if (row && *row->ghs_name)
      if (auto g = ghs_by_name(row->ghs_name)) out.hazards.push_back(ghs_hazard(row->name, *g));
  }
  if (!out.not_modelled.empty()) {
    std::string list;
    for (const auto& n : out.not_modelled) list += (list.empty() ? "" : ", ") + n;
    diag.warn(DiagCode::ThermoApproximation, "not in the solution model, so left out of the chemistry (safety checks still apply): " + list);
  }
  return out;
}

std::optional<LabTitrationResult> lab_titrate(const LabPortion& analyte, const LabComponent& titrant, double max_mL, int points,
                                              Diagnostics& diag) {
  if (!(max_mL > 0.0) || points < 2 || points > 2000) {
    diag.error(DiagCode::InvalidInput, "titration: the maximum volume must be above 0 mL and the points between 2 and 2000");
    return std::nullopt;
  }
  LabTitrationResult out;
  std::vector<std::string> ids;
  for (int i = 0; i <= points; ++i) {
    const double v = max_mL * i / points;
    std::vector<LabPortion> ps = {analyte};
    if (v > 0) ps.push_back({"titrant", v, analyte.T_C, {titrant}});
    std::vector<std::pair<std::string, double>> moles;
    std::vector<std::string> nm;
    double vol = 0.0;
    ids.clear();
    if (!gather(ps, moles, ids, nm, vol, diag)) return std::nullopt;
    if (i == points) out.not_modelled = nm;
    auto s = solve(moles, vol / 1000.0 * kWaterDensity, analyte.T_C, true, {}, diag);
    if (!s) return std::nullopt;
    out.curve.push_back({v, s->r.pH});
  }
  // Equivalence points: local maxima of dpH/dV above a fifth of the largest, refined by a parabola.
  const std::size_t n = out.curve.size();
  std::vector<double> slope(n, 0.0);
  for (std::size_t i = 1; i + 1 < n; ++i)
    slope[i] = (out.curve[i + 1].pH - out.curve[i - 1].pH) / (out.curve[i + 1].volume_mL - out.curve[i - 1].volume_mL);
  double smax = 0.0;
  for (double s : slope) smax = std::max(smax, std::fabs(s));
  for (std::size_t i = 2; i + 2 < n; ++i) {
    const double a = std::fabs(slope[i - 1]), b = std::fabs(slope[i]), c = std::fabs(slope[i + 1]);
    if (b < 0.2 * smax || b < a || b <= c) continue;
    const double denom = a - 2 * b + c;
    const double off = denom != 0.0 ? 0.5 * (a - c) / denom : 0.0;
    const double h = out.curve[i + 1].volume_mL - out.curve[i].volume_mL;
    const double v = out.curve[i].volume_mL + off * h;
    // pH at v by linear interpolation.
    std::size_t j = i;
    while (j + 1 < n && out.curve[j + 1].volume_mL < v) ++j;
    while (j > 0 && out.curve[j].volume_mL > v) --j;
    const auto& p0 = out.curve[j];
    const auto& p1 = out.curve[std::min(j + 1, n - 1)];
    const double t = p1.volume_mL > p0.volume_mL ? (v - p0.volume_mL) / (p1.volume_mL - p0.volume_mL) : 0.0;
    out.equivalence.push_back({v, p0.pH + t * (p1.pH - p0.pH)});
  }
  // Hazards of the final mixture.
  std::vector<LabPortion> fin = {analyte, {"titrant", max_mL, analyte.T_C, {titrant}}};
  Diagnostics quiet;
  if (auto m = lab_mix(fin, quiet)) out.hazards = m->hazards;
  return out;
}

// ---- text input and JSON -------------------------------------------------------------------------------

namespace {

std::string jstr(const std::string& s) { return json_escape(s); }
std::string jnum(double v) { return std::isfinite(v) ? num(v, 10) : "null"; }

std::string hazards_json(const std::vector<LabHazard>& hz) {
  std::string s = "[";
  for (std::size_t i = 0; i < hz.size(); ++i) {
    const auto& h = hz[i];
    s += std::string(i ? "," : "") + "{\"kind\":" + jstr(h.kind) + ",\"severity\":" + jstr(h.severity) + ",\"message\":" +
         jstr(h.message) + ",\"source\":" + jstr(h.source) + ",\"substances\":[";
    for (std::size_t k = 0; k < h.substances.size(); ++k) s += (k ? "," : "") + jstr(h.substances[k]);
    s += "]}";
  }
  return s + "]";
}

std::string strings_json(const std::vector<std::string>& v) {
  std::string s = "[";
  for (std::size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + jstr(v[i]);
  return s + "]";
}

// Lines of words, '#' comments removed.
std::vector<std::vector<std::string>> lines_of(std::string_view text) {
  std::vector<std::vector<std::string>> out;
  std::istringstream in{std::string(text)};
  for (std::string line; std::getline(in, line);) {
    const auto hash = line.find('#');
    if (hash != std::string::npos) line.resize(hash);
    std::istringstream ls(line);
    std::vector<std::string> w;
    for (std::string x; ls >> x;) w.push_back(x);
    if (!w.empty()) out.push_back(w);
  }
  return out;
}

bool to_double(const std::string& s, double& v) {
  char* end = nullptr;
  v = std::strtod(s.c_str(), &end);
  return end != s.c_str() && *end == '\0' && std::isfinite(v);
}

}  // namespace

std::string lab_mix_json(std::string_view text) {
  Diagnostics d;
  std::vector<LabPortion> portions;
  for (const auto& w : lines_of(text)) {
    double v = 0.0;
    if (w[0] == "portion") {
      LabPortion p;
      if (w.size() < 3 || !to_double(w[2], p.volume_mL) || (w.size() > 3 && !to_double(w[3], p.T_C))) {
        d.error(DiagCode::InvalidInput, "write: portion <label> <volume mL> [<T °C>]");
        break;
      }
      p.label = w[1];
      portions.push_back(p);
    } else if (w.size() == 2 && to_double(w[1], v)) {
      if (portions.empty()) {
        d.error(DiagCode::InvalidInput, "a solute line before the first 'portion' line");
        break;
      }
      portions.back().solutes.push_back({w[0], v});
    } else {
      d.error(DiagCode::InvalidInput, "cannot read the line starting '" + w[0] + "'");
      break;
    }
  }
  std::optional<LabMixResult> r;
  if (!d.has_errors()) r = lab_mix(portions, d);
  std::string s = std::string("{\"ok\":") + (r ? "true" : "false");
  if (r) {
    s += ",\"volume_mL\":" + jnum(r->volume_mL) + ",\"T_initial_C\":" + jnum(r->T_initial_C) + ",\"T_final_C\":" + jnum(r->T_final_C) +
         ",\"heat_released_kJ\":" + jnum(r->heat_released_kJ) + ",\"pH\":" + jnum(r->pH) + ",\"ionic_strength\":" +
         jnum(r->ionic_strength) + ",\"precipitates\":[";
    for (std::size_t i = 0; i < r->precipitates.size(); ++i) {
      const auto& p = r->precipitates[i];
      s += std::string(i ? "," : "") + "{\"phase\":" + jstr(p.phase) + ",\"formula\":" + jstr(p.formula) + ",\"mol\":" + jnum(p.mol) +
           ",\"grams\":" + jnum(p.grams) + "}";
    }
    s += "],\"gases\":[";
    for (std::size_t i = 0; i < r->gases.size(); ++i)
      s += std::string(i ? "," : "") + "{\"phase\":" + jstr(r->gases[i].phase) + ",\"mol\":" + jnum(r->gases[i].mol) + ",\"mL\":" +
           jnum(r->gases[i].mL) + "}";
    s += "],\"major_species\":[";
    for (std::size_t i = 0; i < r->major_species.size(); ++i)
      s += std::string(i ? "," : "") + "{\"name\":" + jstr(r->major_species[i].first) + ",\"molality\":" + jnum(r->major_species[i].second) + "}";
    s += "],\"not_modelled\":" + strings_json(r->not_modelled) + ",\"hazards\":" + hazards_json(r->hazards);
  }
  return s + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

std::string lab_titration_json(std::string_view text) {
  Diagnostics d;
  LabPortion analyte;
  LabComponent titrant;
  double max_mL = 0.0, pts = 100.0;
  bool have_analyte = false, have_titrant = false;
  for (const auto& w : lines_of(text)) {
    double v = 0.0;
    if (w[0] == "analyte" && w.size() >= 2 && to_double(w[1], analyte.volume_mL)) {
      have_analyte = true;
      analyte.label = "analyte";
      if (w.size() > 2 && !to_double(w[2], analyte.T_C)) d.error(DiagCode::InvalidInput, "write: analyte <volume mL> [<T °C>]");
    } else if (w[0] == "titrant" && w.size() == 3 && to_double(w[2], titrant.molarity)) {
      titrant.reagent = w[1];
      have_titrant = true;
    } else if (w[0] == "max" && w.size() == 2 && to_double(w[1], max_mL)) {
    } else if (w[0] == "points" && w.size() == 2 && to_double(w[1], pts)) {
    } else if (w.size() == 2 && to_double(w[1], v) && have_analyte) {
      analyte.solutes.push_back({w[0], v});
    } else {
      d.error(DiagCode::InvalidInput, "cannot read the line starting '" + w[0] + "'");
    }
  }
  if (!have_analyte || !have_titrant) d.error(DiagCode::InvalidInput, "need an 'analyte <mL>' line and a 'titrant <reagent> <mol/L>' line");
  std::optional<LabTitrationResult> r;
  if (!d.has_errors()) r = lab_titrate(analyte, titrant, max_mL, static_cast<int>(pts), d);
  std::string s = std::string("{\"ok\":") + (r ? "true" : "false");
  if (r) {
    s += ",\"curve\":[";
    for (std::size_t i = 0; i < r->curve.size(); ++i)
      s += std::string(i ? "," : "") + "[" + jnum(r->curve[i].volume_mL) + "," + jnum(r->curve[i].pH) + "]";
    s += "],\"equivalence\":[";
    for (std::size_t i = 0; i < r->equivalence.size(); ++i)
      s += std::string(i ? "," : "") + "{\"volume_mL\":" + jnum(r->equivalence[i].volume_mL) + ",\"pH\":" + jnum(r->equivalence[i].pH) + "}";
    s += "],\"not_modelled\":" + strings_json(r->not_modelled) + ",\"hazards\":" + hazards_json(r->hazards);
  }
  return s + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

std::string lab_reagents_json() {
  std::string s = "{\"ok\":true,\"reagents\":[";
  bool first = true;
  for (const auto& r : kReagents) {
    s += std::string(first ? "" : ",") + "{\"id\":" + jstr(r.id) + ",\"name\":" + jstr(r.name) + ",\"formula\":" + jstr(r.formula) +
         ",\"modelled\":" + (*r.formula ? "true" : "false") + ",\"smiles\":" + jstr(r.smiles) + ",\"groups\":" +
         strings_json(split_words(r.groups));
    if (*r.ghs_name)
      if (auto g = ghs_by_name(r.ghs_name)) {
        s += ",\"ghs\":{\"signal\":" + jstr(g->signal) + ",\"cid\":" + std::to_string(g->cid) + ",\"pictograms\":" + strings_json(g->pictograms) +
             ",\"hazards\":[";
        for (std::size_t i = 0; i < g->hazards.size(); ++i)
          s += std::string(i ? "," : "") + "{\"code\":" + jstr(g->hazards[i].first) + ",\"percent\":" + jnum(g->hazards[i].second) +
               ",\"text\":" + jstr(ghs_phrase(g->hazards[i].first)) + "}";
        s += "]}";
      }
    s += "}";
    first = false;
  }
  return s + "],\"ghs_compounds\":" + std::to_string(ghs_table_size()) + "}";
}

std::string ghs_json(std::string_view smiles) {
  auto g = ghs_lookup(smiles);
  std::string s = std::string("{\"ok\":true,\"found\":") + (g ? "true" : "false");
  if (g) {
    s += ",\"name\":" + jstr(g->name) + ",\"cid\":" + std::to_string(g->cid) + ",\"signal\":" + jstr(g->signal) +
         ",\"pictograms\":" + strings_json(g->pictograms) + ",\"hazards\":[";
    for (std::size_t i = 0; i < g->hazards.size(); ++i)
      s += std::string(i ? "," : "") + "{\"code\":" + jstr(g->hazards[i].first) + ",\"percent\":" + jnum(g->hazards[i].second) +
           ",\"text\":" + jstr(ghs_phrase(g->hazards[i].first)) + "}";
    s += "]";
  }
  return s + "}";
}

}  // namespace chemlarp
