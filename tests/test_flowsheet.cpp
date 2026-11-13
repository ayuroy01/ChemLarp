// Flowsheet parsing and the sequential-modular solver with recycle.
//
// Reference solutions are independent of the solver: closed-form balances for the propane dehydrogenation
// recycle of Felder & Rousseau, Elementary Principles of Chemical Processes (3rd ed.), Example 4.7-2, and the
// reactor-feed basis of their Example 4.7-3 (CO2 hydrogenation with purge), solved by hand algebra.
// Equation-oriented checks against thermo's correlations are in "Flowsheet solver: validation" in
// docs/ARCHITECTURE.md.
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>

#include "chemlarp/chemlarp_c.h"
#include "chemlarp/flowsheet.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

FlowsheetResult solve(const std::string& text, Diagnostics& d) {
  auto spec = parse_flowsheet(text, d);
  CHECK(spec.has_value());
  auto r = solve_flowsheet(*spec, d);
  CHECK(r.has_value());
  return *r;
}

const Stream& stream(const FlowsheetResult& r, const std::string& name) {
  for (const auto& s : r.streams)
    if (s.name == name) return s.stream;
  throw ::cltest::Failure("no stream " + name);
}

bool has(const Diagnostics& d, DiagCode code) {
  for (const auto& it : d.items())
    if (it.code == code) return true;
  return false;
}

std::string num(double v) {
  char b[40];
  std::snprintf(b, sizeof b, "%.17g", v);
  return b;
}

// Felder & Rousseau Example 4.7-2: 95 % overall conversion of propane with 0.555 % of the reactor-outlet
// propane and 95 % of the propylene leaving in the product.
std::string propane(const char* method) {
  const double X = 95.0 / (95.0 + 5.0 / 0.00555);
  return std::string("components CCC C=CC [HH]\n"
                     "feed F T=300 P=1bar flows CCC=100\n"
                     "mixer M F R -> Q1\n"
                     "reactor RX Q1 -> Q2 key=CCC X=") +
         num(X) +
         " T=800 nu CCC=-1 C=CC=1 [HH]=1\n"
         "separator S Q2 -> P R T=300 split CCC=0.00555 C=CC=0.95 [HH]=1  # product, recycle\n"
         "tolerance 1e-10\nmax-iterations 1000\nmethod " +
         method + "\n";
}

}  // namespace

TEST(flowsheet_propane_recycle_closed_form) {
  const double a = 1 - 0.00555, b = 0.05, X = 95.0 / (95.0 + 5.0 / 0.00555);
  const double Pin = 100.0 / (1 - a * (1 - X)), Yout = X * Pin / (1 - b);
  int iters[2] = {0, 0};
  int k = 0;
  for (const char* m : {"wegstein", "direct"}) {
    Diagnostics d;
    const FlowsheetResult r = solve(propane(m), d);
    CHECK(r.converged);
    CHECK(!d.has_errors());
    CHECK_EQ(r.tears.size(), std::size_t{1});
    CHECK_EQ(r.tears[0], std::string("R"));
    iters[k++] = r.iterations;
    const Stream& q1 = stream(r, "Q1");
    CHECK_NEAR(q1.n[0], Pin, 1e-6 * Pin);
    CHECK_NEAR(q1.n[1], b * Yout, 1e-6);
    const Stream& p = stream(r, "P");
    CHECK_NEAR(p.n[0], 5.0, 1e-6);
    CHECK_NEAR(p.n[1], 95.0, 1e-6);
    CHECK_NEAR(p.n[2], 95.0, 1e-6);
    // Published: 995.9 mol propane to the reactor per 100 mol fresh feed; recycle ratio 9.0.
    CHECK_NEAR(q1.n[0], 995.9, 0.05);
    const Stream& rr = stream(r, "R");
    CHECK_NEAR((rr.n[0] + rr.n[1]) / 100.0, 9.0, 0.01);
    CHECK_NEAR(r.mass_out_kgs, r.mass_in_kgs, 1e-8 * r.mass_in_kgs);
    CHECK_NEAR(r.energy_imbalance_W, 0.0, 1e-6 * std::fabs(q1.H_W));
  }
  CHECK(iters[0] < iters[1] / 4);  // Wegstein accelerates the slow (factor 0.9) direct substitution
}

TEST(flowsheet_purge_loop_reactor_feed) {
  // Felder & Rousseau Example 4.7-3 statement: reactor feed 28 % CO2, 70 % H2, 2 % inerts; 60 % single-pass
  // H2 conversion; condenser removes all methanol and water; fresh feed 0.4 % inerts. Basis 100 mol/s
  // reactor feed: condenser gas 14 CO2, 28 H2, 2 I; fresh feed n0 from the inert balance.
  const double yI = 2.0 / 44.0, n0 = (2.0 - 100.0 * yI) / (0.004 - yI), rec = 100.0 - n0;
  const double fC = 28.0 - rec * 14.0 / 44.0, fH = 70.0 - rec * 28.0 / 44.0, fI = 2.0 - rec * yI;
  const double purge = 1.0 - rec / 44.0;
  const std::string text = "components O=C=O [HH] [Ar] CO O\nfeed F T=300 P=50bar flows O=C=O=" + num(fC) +
                           " [HH]=" + num(fH) + " [Ar]=" + num(fI) +
                           "\nmixer M F R -> Q1\n"
                           "reactor RX Q1 -> Q2 key=[HH] X=0.6 T=500 nu O=C=O=-1 [HH]=-3 CO=1 O=1\n"
                           "separator C Q2 -> G L T=300 split O=C=O=1 [HH]=1 [Ar]=1\n"
                           "splitter SP G -> R PURGE fractions=" +
                           num(1 - purge) + "," + num(purge) + "\ntolerance 1e-10\n";
  Diagnostics d;
  const FlowsheetResult r = solve(text, d);
  CHECK(r.converged);
  CHECK_NEAR(fI / n0, 0.004, 1e-12);
  const Stream& q1 = stream(r, "Q1");
  CHECK_NEAR(q1.n[0], 28.0, 1e-6);
  CHECK_NEAR(q1.n[1], 70.0, 1e-6);
  CHECK_NEAR(q1.n[2], 2.0, 1e-6);
  CHECK_NEAR(stream(r, "L").n[3], 14.0, 1e-6);
  CHECK_NEAR(stream(r, "PURGE").total(), 44.0 - rec, 1e-6);
  CHECK_NEAR(r.mass_out_kgs, r.mass_in_kgs, 1e-8 * r.mass_in_kgs);
}

TEST(flowsheet_flash_recycle_is_a_fixed_point) {
  const std::string text =
      "components c1ccccc1 Cc1ccccc1\n"
      "feed F T=300 P=1atm flows c1ccccc1=50 Cc1ccccc1=50\n"
      "mixer M F R -> S1\nheater H S1 -> S2 T=368\nflash D S2 -> V L T=368 P=1atm\n"
      "splitter SP L -> R B fractions=0.6,0.4\n";
  Diagnostics d;
  const FlowsheetResult r = solve(text, d);
  CHECK(r.converged);
  CHECK(r.residuals.back() < 1e-7);
  // Re-running the loop's units from the converged mixer outlet reproduces the recycle.
  Diagnostics q;
  auto again = flash_drum(r.components, stream(r, "S2"), 368.0, 101325.0, q);
  CHECK(again.has_value());
  const Stream& R = stream(r, "R");
  for (std::size_t i = 0; i < 2; ++i) {
    CHECK_NEAR(0.6 * again->liquid.n[i], R.n[i], 1e-6);
    CHECK_NEAR(stream(r, "S1").n[i], 50.0 + R.n[i], 1e-9);
    // Component balance over the whole flowsheet: feed = vapour + bottoms.
    CHECK_NEAR(stream(r, "V").n[i] + stream(r, "B").n[i], 50.0, 1e-6);
  }
  CHECK_NEAR(r.energy_imbalance_W, 0.0, 1e-3);
  // Direct substitution reaches the same point.
  Diagnostics d2;
  const FlowsheetResult r2 = solve(text + "method direct\ntolerance 1e-10\nmax-iterations 500\n", d2);
  CHECK(r2.converged);
  CHECK_NEAR(stream(r2, "V").n[0], stream(r, "V").n[0], 1e-5);
}

TEST(flowsheet_energy_recycle) {
  // Gas loop with a heater: at steady state the product carries the duty, F (h(T) − h(300 K)) = Q.
  const std::string text =
      "components N#N\nfeed F T=300 P=1bar flows N#N=100\nmixer M F R -> S1\nheater H S1 -> S2 Q=0.5MW\n"
      "splitter SP S2 -> R OUT fractions=0.5,0.5\n";
  Diagnostics d;
  const FlowsheetResult r = solve(text, d);
  CHECK(r.converged);
  CHECK_NEAR(stream(r, "OUT").H_W - stream(r, "F").H_W, 5e5, 1e-3);
  CHECK_NEAR(r.energy_imbalance_W, 0.0, 1e-3);
  // thermo (HeatCapacityGas, nitrogen) gives 470.968 K; ChemLarp's Cp polynomial differs by 2e-4.
  CHECK_NEAR(stream(r, "OUT").T_K, 470.968, 0.2);
  CHECK_NEAR(stream(r, "S2").T_K, stream(r, "OUT").T_K, 1e-9);
}

TEST(flowsheet_feed_effluent_exchanger) {
  // Toluene hydrodealkylation: the cold feed is preheated by the adiabatic reactor effluent, cooled to 500 K.
  const std::string text =
      "components Cc1ccccc1 [HH] c1ccccc1 C\n"
      "feed F T=300 P=25bar flows Cc1ccccc1=10 [HH]=50\n"
      "exchanger E Q2 F -> PROD Q0 Thot=500\n"
      "heater H Q0 -> Q1 T=890\n"
      "reactor RX Q1 -> Q2 key=Cc1ccccc1 X=0.75 nu Cc1ccccc1=-1 [HH]=-1 c1ccccc1=1 C=1\n";
  // From an empty recycle the exchanger passes the cold feed through untouched; the solution follows.
  Diagnostics d;
  const FlowsheetResult r = solve(text, d);
  CHECK(r.converged);
  CHECK_NEAR(stream(r, "PROD").T_K, 500.0, 1e-6);
  // Overall: product enthalpy = feed enthalpy + heater duty (exchanger internal, reactor adiabatic).
  double heater = 0, ex = 0;
  for (const auto& u : r.units) {
    if (u.name == "H") heater = u.duty_W;
    if (u.name == "E")
      for (const auto& [k, v] : u.values)
        if (k == "exchanger_duty_W") ex = v;
  }
  CHECK(heater > 0 && ex > heater);
  CHECK_NEAR(stream(r, "PROD").H_W, stream(r, "F").H_W + heater, 1e-3);
  CHECK_NEAR(r.energy_imbalance_W, 0.0, 1e-3);
  CHECK_NEAR(stream(r, "Q0").H_W - stream(r, "F").H_W, ex, 1e-6 * ex);
  // ΔrH = 82.9 − 74.6 − 50.5 = −42.2 kJ/mol (tabulated ΔfH° of benzene, methane, toluene) heats the reactor.
  CHECK(stream(r, "Q2").T_K > 890.0);
}

// A heater duty and a hot-outlet temperature together fix the product enthalpy, so the loop has a solution
// only for one duty; otherwise the tear stream drifts and the solver must not report convergence.
TEST(flowsheet_ill_posed_loop_does_not_converge) {
  Diagnostics d;
  auto spec = parse_flowsheet(
      "components Cc1ccccc1 [HH] c1ccccc1 C\nfeed F T=300 P=25bar flows Cc1ccccc1=10 [HH]=50\n"
      "exchanger E Q2 F -> PROD Q0 Thot=500\nheater H Q0 -> Q1 Q=650kW\n"
      "reactor RX Q1 -> Q2 key=Cc1ccccc1 X=0.75 nu Cc1ccccc1=-1 [HH]=-1 c1ccccc1=1 C=1\n"
      "guess Q2 T=900 P=25bar flows Cc1ccccc1=2.5 [HH]=42.5 c1ccccc1=7.5 C=7.5\nmax-iterations 40\n",
      d);
  CHECK(spec.has_value());
  auto r = solve_flowsheet(*spec, d);
  CHECK(!r || !r->converged);
  CHECK(d.has_errors());
}

TEST(flowsheet_without_recycle_is_one_pass) {
  Diagnostics d;
  const FlowsheetResult r = solve(
      "components O\nfeed W T=25C P=1bar flows O=1\npump P1 W -> W2 P=10bar eff=0.75\n"
      "heater H W2 -> S T=500 # steam\nvalve V S -> S2 P=2bar\n",
      d);
  CHECK(r.converged);
  CHECK_EQ(r.iterations, 1);
  CHECK(r.tears.empty());
  CHECK_NEAR(stream(r, "W").T_K, 298.15, 1e-12);
  CHECK_NEAR(stream(r, "W2").P_Pa, 1e6, 1e-6);
  CHECK(r.units[0].work_W > 0);
  CHECK_NEAR(r.energy_imbalance_W, 0.0, 1e-6);
}

TEST(flowsheet_user_tear_and_guess) {
  // Tearing at the reactor feed instead of the recycle, with a guess, gives the same solution.
  Diagnostics d;
  std::string text = propane("wegstein");
  text += "guess Q1 T=300 P=1bar flows CCC=900 C=CC=4\n";
  const FlowsheetResult r = solve(text, d);
  CHECK(r.converged);
  CHECK_EQ(r.tears.size(), std::size_t{1});
  CHECK_EQ(r.tears[0], std::string("Q1"));
  CHECK_EQ(r.order[0], std::string("RX"));
  CHECK_NEAR(stream(r, "P").n[0], 5.0, 1e-6);
}

TEST(flowsheet_reports_non_convergence) {
  Diagnostics d;
  auto spec = parse_flowsheet(propane("direct") + "max-iterations 3\n", d);
  CHECK(spec.has_value());
  auto r = solve_flowsheet(*spec, d);
  CHECK(r.has_value());
  CHECK(!r->converged);
  CHECK_EQ(r->iterations, 3);
  CHECK(has(d, DiagCode::RecycleNotConverged));
}

TEST(flowsheet_closed_loop_without_feed_warns) {
  // The feed is wired to the hot side, so the loop through the cold side never receives material.
  Diagnostics d;
  const FlowsheetResult r = solve(
      "components CCC\nfeed F T=900 P=2bar flows CCC=1\nexchanger E F Q2 -> PROD Q0 Q=1kW\n"
      "heater H Q0 -> Q2 T=300\n",
      d);
  CHECK(r.converged);
  CHECK(has(d, DiagCode::FlowsheetTopology));
}

TEST(flowsheet_syntax_errors) {
  const char* bad[] = {
      "components O\nfeed F T=300 P=1bar flows X=1\nheater H F -> G T=350\n",        // unknown component
      "components O\nfeed F T=300 P=1bar flows O=1\nheater H F -> G\n",               // heater without T or Q
      "components O\nfeed F T=300 P=1bar flows O=1\nboiler H F -> G T=350\n",         // unknown unit
      "components O\nfeed F T=300 P=1bar flows O=1\nheater H F -> G T=350 Z=1\n",     // unknown parameter
      "components O\nfeed F T=300 P=1furlong flows O=1\nheater H F -> G T=350\n",     // unknown unit of measure
      "components O\nfeed F T=300 P=1bar flows O=1\nsplitter S F -> A B fractions=1\n",  // one fraction for two outlets
      "components O\nfeed F T=300 P=1bar flows O=1\n",                                 // no units
      "feed F T=300 P=1bar flows O=1\nheater H F -> G T=350\n",                        // no components
  };
  for (const char* t : bad) {
    Diagnostics d;
    CHECK(!parse_flowsheet(t, d).has_value());
    CHECK(has(d, DiagCode::FlowsheetSyntax));
  }
  // '#' inside a SMILES is not a comment; units of measure convert.
  Diagnostics d;
  auto s = parse_flowsheet("components N#N # nitrogen\nfeed F T=80C P=1.5bar flows N#N=2\nheater H F -> G Q=-2kW\n", d);
  CHECK(s.has_value());
  CHECK_EQ(s->components.size(), std::size_t{1});
  CHECK_EQ(s->components[0], std::string("N#N"));
  CHECK_NEAR(s->feeds[0].T_K, 353.15, 1e-12);
  CHECK_NEAR(s->feeds[0].P_Pa, 1.5e5, 1e-9);
  CHECK_NEAR(s->units[0].duty_W, -2000.0, 1e-12);
}

TEST(flowsheet_topology_errors) {
  const char* bad[] = {
      "components O\nfeed F T=300 P=1bar flows O=1\nheater H X -> G T=350\n",                      // no source
      "components O\nfeed F T=300 P=1bar flows O=1\nheater H F -> G T=350\nheater H2 F -> K T=350\n",  // two destinations
      "components O\nfeed F T=300 P=1bar flows O=1\nheater H F -> F T=350\n",                      // two sources
  };
  for (const char* t : bad) {
    Diagnostics d;
    auto s = parse_flowsheet(t, d);
    CHECK(s.has_value());
    CHECK(!solve_flowsheet(*s, d).has_value());
    CHECK(has(d, DiagCode::FlowsheetTopology));
  }
}

TEST(flowsheet_c_api_json) {
  char* out = cl_flowsheet_json(
      "components c1ccccc1 Cc1ccccc1\nfeed F T=300 P=1atm flows c1ccccc1=50 Cc1ccccc1=50\n"
      "mixer M F R -> S1\nflash D S1 -> V L T=368 P=1atm\nsplitter SP L -> R B fractions=0.5,0.5\n");
  CHECK(out != nullptr);
  const std::string s(out);
  cl_free_string(out);
  CHECK(s.rfind("{\"ok\":true", 0) == 0);
  CHECK(s.find("\"tears\":[\"R\"]") != std::string::npos);
  CHECK(s.find("\"streams\":[{\"name\":\"F\"") != std::string::npos);
  CHECK(s.find("\"energy_imbalance_W\"") != std::string::npos);
  char* bad = cl_flowsheet_json("components O\nheater H A -> B T=1\n");
  CHECK(std::string(bad).find("FlowsheetTopology") != std::string::npos);
  cl_free_string(bad);
}
