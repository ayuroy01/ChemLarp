// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
// Virtual lab: mixing, titrations, safety rules and GHS data. scripts/validate_lab.py compares the chemistry with
// PHREEQC; these tests pin the behaviour.
#include <cmath>
#include <string>

#include "chemlarp/aqueous.hpp"
#include "chemlarp/chemlarp_c.h"
#include "chemlarp/lab.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

bool has_hazard(const std::vector<LabHazard>& hz, const std::string& kind, const std::string& text) {
  for (const auto& h : hz)
    if (h.kind == kind && h.message.find(text) != std::string::npos) return true;
  return false;
}

}  // namespace

TEST(lab_neutralisation_heat_and_ph) {
  Diagnostics d;
  const auto r = lab_mix({{"acid", 50, 25, {{"hcl", 1.0}}}, {"base", 50, 25, {{"naoh", 1.0}}}}, d);
  CHECK(r.has_value());
  if (!r) return;
  // 0.05 mol x ~56 kJ/mol over 100 g x 4.184 J/(g K): about 6.7 K.
  CHECK(std::fabs(r->T_final_C - r->T_initial_C - 6.74) < 0.05);
  CHECK(std::fabs(r->heat_released_kJ - 2.818) < 0.01);
  CHECK(std::fabs(r->pH - 6.85) < 0.02);  // neutral at 31.7 °C
  CHECK(has_hazard(r->hazards, "computed", "Releases"));
  CHECK(has_hazard(r->hazards, "ghs", "Sodium hydroxide"));
  CHECK(!has_hazard(r->hazards, "computed", "corrosive"));
}

TEST(lab_precipitates_gases_and_corrosive_ph) {
  Diagnostics d;
  auto r = lab_mix({{"a", 50, 25, {{"bacl2", 0.1}}}, {"b", 50, 25, {{"na2so4", 0.1}}}}, d);
  CHECK(r && r->precipitates.size() == 1);
  if (r && !r->precipitates.empty()) {
    CHECK_EQ(r->precipitates[0].phase, std::string("Barite"));
    CHECK(std::fabs(r->precipitates[0].mol - 0.005) < 1e-5);
    CHECK(std::fabs(r->precipitates[0].grams - 0.005 * 233.39) < 0.01);  // BaSO4
  }
  r = lab_mix({{"a", 50, 25, {{"na2s", 0.5}}}, {"b", 50, 25, {{"hcl", 1.0}}}}, d);
  CHECK(r && r->gases.size() == 1 && r->gases[0].phase == "H2S(g)");
  if (r) CHECK(has_hazard(r->hazards, "computed", "H2S is given off"));
  // The danger rating of the gas comes from its GHS classification (hydrogen sulfide: Danger).
  bool danger = false;
  if (r)
    for (const auto& h : r->hazards) danger = danger || (h.kind == "computed" && h.message.find("H2S") != std::string::npos && h.severity == "danger");
  CHECK(danger);
  r = lab_mix({{"a", 10, 25, {{"naoh", 1.0}}}, {"b", 90, 25, {}}}, d);
  CHECK(r && r->pH > 12.5);
  if (r) CHECK(has_hazard(r->hazards, "computed", "corrosive"));
  // Calcite: nearly all the calcium precipitates (the case that used to stall the solver short of convergence).
  r = lab_mix({{"a", 40, 25, {{"cacl2", 0.05}}}, {"b", 60, 25, {{"na2co3", 0.05}}}}, d);
  CHECK(r && r->precipitates.size() == 1 && std::fabs(r->precipitates[0].mol - 0.002) < 1e-5);
}

TEST(lab_incompatibility_rules) {
  auto fires = [](std::vector<std::string> ids, bool aqueous, const std::string& text) {
    for (const auto& h : lab_incompatibilities(ids, aqueous))
      if (h.message.find(text) != std::string::npos && h.source.find("Prudent Practices") != std::string::npos) return true;
    return false;
  };
  CHECK(fires({"nh3", "naocl"}, true, "hypochlorite"));
  CHECK(fires({"kcn", "hcl"}, true, "hydrogen cyanide"));
  CHECK(fires({"na", "dcm"}, false, "halogenated"));
  CHECK(fires({"na"}, true, "violently with water"));
  CHECK(!fires({"na"}, false, "violently with water"));  // no water, no reaction
  CHECK(fires({"kmno4", "h2so4"}, false, "permanganate"));
  CHECK(fires({"hno3", "ethanol"}, false, "alcohols"));
  CHECK(fires({"h2o2", "acetone"}, false, "organic material"));
  CHECK(fires({"nh3", "i2"}, true, "triiodide"));
  CHECK(fires({"kclo3", "nabh4"}, false, "reducing agents"));
  // Benign pairs: no rule.
  CHECK(lab_incompatibilities({"nacl", "kcl"}, true).empty());
  CHECK(lab_incompatibilities({"hcl", "naoh"}, true).empty());  // heat is computed, not a rule
  // A substance does not react with itself through two of its own groups (nitric acid is acid and oxidizer).
  CHECK(lab_incompatibilities({"hno3"}, false).empty());
}

TEST(lab_not_modelled_and_ghs) {
  Diagnostics d;
  const auto r = lab_mix({{"a", 50, 25, {{"naocl", 0.7}}}, {"b", 50, 25, {{"nh3", 1.0}}}}, d);
  CHECK(r.has_value());
  if (r) {
    CHECK_EQ(r->not_modelled.size(), std::size_t{1});
    CHECK(has_hazard(r->hazards, "incompatibility", "hypochlorite"));
  }
  CHECK(d.contains(DiagCode::ThermoApproximation));
  // ghs_lookup takes ChemLarp's canonical SMILES; the C API canonicalizes first, so any spelling works there.
  const auto g = ghs_lookup("[OH-].[Na+]");
  for (const char* spelling : {"[Na+].[OH-]", "[OH-].[Na+]"}) {
    char* j = cl_ghs_json(spelling);
    CHECK(std::string(j).find("\"cid\":14798") != std::string::npos);
    cl_free_string(j);
  }
  CHECK(g && g->signal == "Danger" && g->cid == 14798);
  if (g) CHECK(!g->hazards.empty() && g->hazards[0].first == "H314");
  CHECK_EQ(ghs_phrase("H314"), std::string("Causes severe skin burns and eye damage"));
  CHECK(!ghs_lookup("CCCCCCCCCCCCCCCCCCCCCCCCC"));
  CHECK(ghs_table_size() > 150);
  // Every catalogue reagent with a GHS name finds its entry (lab_reagents_json lists them with "ghs").
  const std::string cat = lab_reagents_json();
  CHECK(cat.find("\"id\":\"naocl\"") != std::string::npos);
  CHECK(cat.find("\"ghs\":{\"signal\":\"Danger\"") != std::string::npos);
}

TEST(lab_titrations_find_equivalence_points) {
  Diagnostics d;
  auto r = lab_titrate({"a", 25, 25, {{"hcl", 0.1}}}, {"naoh", 0.1}, 50, 200, d);
  CHECK(r && r->equivalence.size() == 1);
  if (r && !r->equivalence.empty()) {
    CHECK(std::fabs(r->equivalence[0].volume_mL - 25.0) < 0.1);
    CHECK(std::fabs(r->equivalence[0].pH - 7.0) < 0.1);
  }
  r = lab_titrate({"a", 25, 25, {{"h3po4", 0.1}}}, {"naoh", 0.1}, 80, 320, d);
  CHECK(r && r->equivalence.size() == 2);
  if (r && r->equivalence.size() == 2) {
    CHECK(std::fabs(r->equivalence[0].volume_mL - 25.0) < 0.1);
    CHECK(std::fabs(r->equivalence[1].volume_mL - 50.0) < 0.1);
  }
  CHECK(!lab_titrate({"a", 25, 25, {{"hcl", 0.1}}}, {"naoh", 0.1}, 0, 10, d));
}

TEST(lab_text_input_errors) {
  CHECK(lab_mix_json("portion a 50\nhcl 1\nportion b 50\nnaoh 1\n").rfind("{\"ok\":true", 0) == 0);
  CHECK(lab_mix_json("hcl 1\n").find("before the first 'portion'") != std::string::npos);
  CHECK(lab_mix_json("portion a\n").find("portion <label>") != std::string::npos);
  CHECK(lab_mix_json("portion a 0\nhcl 1\n").find("volume above 0") != std::string::npos);
  CHECK(lab_titration_json("analyte 25\nhcl 0.1\n").find("titrant") != std::string::npos);
}
