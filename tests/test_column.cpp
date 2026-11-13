// Rigorous distillation column.
//
// Reference values for benzene/toluene come from an independent equation-oriented solve of the same
// stage model with thermo's correlations (scripts/validate_column.py); the stage-by-stage check below
// re-flashes every stage's inflows with the separate PH flash of process.cpp, which isolates the column
// solver from the property data.
#include <cmath>
#include <string>
#include <vector>

#include "chemlarp/column.hpp"
#include "chemlarp/flowsheet.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {

ProcessComponents comps(const std::vector<std::string>& smiles, LiquidModel model = LiquidModel::Ideal) {
  Diagnostics d;
  auto c = process_components(smiles, model, d);
  CHECK(c.has_value());
  return *c;
}

Stream feed(const ProcessComponents& c, double T, std::vector<double> n) {
  Diagnostics d;
  auto s = stream_tp(c, T, 101325.0, n, d);
  CHECK(s.has_value());
  return *s;
}

ColumnSpec bt_spec(const ProcessComponents& c) {
  ColumnSpec s;
  s.stages = 15;
  s.feeds = {{8, feed(c, 360.0, {50.0, 50.0})}};
  s.P_top_Pa = 101325.0;
  s.reflux_ratio = 2.0;
  s.distillate_rate = 50.0;
  return s;
}

bool has(const Diagnostics& d, DiagCode code) {
  for (const auto& it : d.items())
    if (it.code == code) return true;
  return false;
}

}  // namespace

TEST(column_benzene_toluene_against_independent_solver) {
  const ProcessComponents c = comps({"c1ccccc1", "Cc1ccccc1"});
  const ColumnSpec s = bt_spec(c);
  Diagnostics d;
  auto r = distillation_column(c, s, d);
  CHECK(r.has_value());
  CHECK(r->residual < 1e-9);
  // Independent solve (thermo vapour pressures and enthalpies): 48.721 mol/s benzene overhead,
  // 353.73 K top, 382.55 K bottom, condenser −4636 kW, reboiler 4798 kW.
  CHECK_NEAR(r->distillate.n[0], 48.721, 0.5);
  CHECK_NEAR(r->stages.front().T_K, 353.73, 0.3);
  CHECK_NEAR(r->stages.back().T_K, 382.55, 0.5);
  CHECK_NEAR(r->condenser_duty_W, -4.636e6, 0.02 * 4.636e6);
  CHECK_NEAR(r->reboiler_duty_W, 4.798e6, 0.02 * 4.798e6);
  // Specifications and balances.
  CHECK_NEAR(r->distillate.total(), 50.0, 1e-8);
  CHECK_NEAR(r->stages[0].L, 2.0 * 50.0, 1e-8);
  for (std::size_t i = 0; i < 2; ++i) CHECK_NEAR(r->distillate.n[i] + r->bottoms.n[i], 50.0, 1e-8);
  const double E = s.feeds[0].stream.H_W + r->condenser_duty_W + r->reboiler_duty_W - r->distillate.H_W - r->bottoms.H_W;
  CHECK_NEAR(E, 0.0, 1e-3);
  // Temperatures rise down the column.
  for (std::size_t j = 1; j < r->stages.size(); ++j) CHECK(r->stages[j].T_K > r->stages[j - 1].T_K);
}

TEST(column_every_stage_is_an_adiabatic_flash) {
  // The inflows of each intermediate stage, flashed at the stage pressure with their enthalpy, give the
  // stage's temperature and its vapour and liquid flows.
  const ProcessComponents c = comps({"c1ccccc1", "Cc1ccccc1"});
  const ColumnSpec s = bt_spec(c);
  Diagnostics d;
  auto r = distillation_column(c, s, d);
  CHECK(r.has_value());
  const auto& st = r->stages;
  for (std::size_t j = 1; j + 1 < st.size(); ++j) {
    std::vector<double> n(2);
    double H = 0.0;
    auto hl = component_enthalpies(c, st[j - 1].T_K, st[j - 1].P_Pa, false, d);
    auto hv = component_enthalpies(c, st[j + 1].T_K, st[j + 1].P_Pa, true, d);
    CHECK(hl && hv);
    for (std::size_t i = 0; i < 2; ++i) {
      // The liquid from the condenser is the reflux; the vapour comes from the stage below.
      const double li = st[j - 1].L * st[j - 1].x[i], vi = st[j + 1].V * st[j + 1].y[i];
      n[i] = li + vi;
      H += li * (*hl)[i] + vi * (*hv)[i];
    }
    if (j + 1 == s.feeds[0].stage) {
      for (std::size_t i = 0; i < 2; ++i) n[i] += s.feeds[0].stream.n[i];
      H += s.feeds[0].stream.H_W;
    }
    auto f = stream_ph(c, st[j].P_Pa, H, n, d, st[j].T_K);
    CHECK(f.has_value());
    CHECK_NEAR(f->T_K, st[j].T_K, 1e-5);
    CHECK_NEAR(f->vapor_fraction * f->total(), st[j].V, 1e-5);
    CHECK_NEAR(f->y[0], st[j].y[0], 1e-7);
    CHECK_NEAR(f->x[0], st[j].x[0], 1e-7);
  }
}

TEST(column_partial_condenser_and_reflux_trend) {
  const ProcessComponents c = comps({"c1ccccc1", "Cc1ccccc1"});
  ColumnSpec s = bt_spec(c);
  s.condenser = Condenser::Partial;
  Diagnostics d;
  auto r = distillation_column(c, s, d);
  CHECK(r.has_value());
  CHECK_EQ(r->distillate.vapor_fraction, 1.0);
  // The vapour distillate is in equilibrium with the reflux: y = K(T, x) x.
  auto K = k_values(c, r->stages[0].T_K, r->stages[0].P_Pa, r->stages[0].x, d);
  CHECK(K.has_value());
  for (std::size_t i = 0; i < 2; ++i) CHECK_NEAR(r->distillate.x[i], (*K)[i] * r->stages[0].x[i], 1e-9);
  // A partial condenser is an extra stage: purer distillate than with a total condenser.
  auto t = distillation_column(c, bt_spec(c), d);
  CHECK(t.has_value());
  CHECK(r->distillate.n[0] > t->distillate.n[0]);
  // More reflux, purer distillate, more reboiler duty.
  ColumnSpec s3 = bt_spec(c);
  s3.reflux_ratio = 4.0;
  auto r3 = distillation_column(c, s3, d);
  CHECK(r3.has_value());
  CHECK(r3->distillate.n[0] > t->distillate.n[0]);
  CHECK(r3->reboiler_duty_W > t->reboiler_duty_W);
  // Bottoms specification instead of distillate: same column.
  ColumnSpec sb = bt_spec(c);
  sb.distillate_rate = 0.0;
  sb.bottoms_rate = 50.0;
  auto rb = distillation_column(c, sb, d);
  CHECK(rb.has_value());
  CHECK_NEAR(rb->distillate.n[0], t->distillate.n[0], 1e-6);
}

TEST(column_ethanol_water_stays_below_azeotrope) {
  // UNIFAC: the ethanol–water azeotrope (about 0.89 ethanol, 351.1 K at 1 atm) bounds the distillate.
  const ProcessComponents c = comps({"CCO", "O"}, LiquidModel::Unifac);
  ColumnSpec s;
  s.stages = 20;
  s.feeds = {{12, feed(c, 350.0, {10.0, 90.0})}};
  s.P_top_Pa = 101325.0;
  s.reflux_ratio = 3.0;
  s.distillate_rate = 11.0;
  Diagnostics d;
  auto r = distillation_column(c, s, d);
  CHECK(r.has_value());
  const double xD = r->distillate.x[0];
  // Independent solve with thermo's UNIFAC and vapour pressures: 0.8364.
  CHECK_NEAR(xD, 0.8364, 0.01);
  CHECK(xD < 0.894);
  CHECK(r->stages.front().T_K > 351.0);
  for (std::size_t i = 0; i < 2; ++i) CHECK_NEAR(r->distillate.n[i] + r->bottoms.n[i], s.feeds[0].stream.n[i], 1e-8);
}

TEST(column_input_errors) {
  const ProcessComponents c = comps({"c1ccccc1", "Cc1ccccc1"});
  auto expect_error = [&](ColumnSpec s) {
    Diagnostics d;
    CHECK(!distillation_column(c, s, d).has_value());
    CHECK(d.has_errors());
  };
  ColumnSpec s = bt_spec(c);
  s.bottoms_rate = 50.0;  // both D and B
  expect_error(s);
  s = bt_spec(c);
  s.distillate_rate = 150.0;  // more than the feed
  expect_error(s);
  s = bt_spec(c);
  s.stages = 2;
  expect_error(s);
  s = bt_spec(c);
  s.feeds[0].stage = 16;
  expect_error(s);
  s = bt_spec(c);
  s.reflux_ratio = 0.0;
  expect_error(s);
  // A non-condensable component (no vapour-pressure equation) is outside the column model.
  const ProcessComponents g = comps({"c1ccccc1", "O=C=O"});
  ColumnSpec sg = bt_spec(g);
  sg.feeds[0].stream = feed(g, 360.0, {50.0, 1.0});
  Diagnostics d;
  CHECK(!distillation_column(g, sg, d).has_value());
  CHECK(has(d, DiagCode::InvalidMixtureInput));
}

TEST(column_in_a_flowsheet_recycle) {
  // Part of the bottoms is recycled to the feed; the column converges inside the recycle loop.
  Diagnostics d;
  auto spec = parse_flowsheet(
      "components c1ccccc1 Cc1ccccc1\n"
      "feed F T=360 P=1atm flows c1ccccc1=40 Cc1ccccc1=60\n"
      "mixer M F R -> S\n"
      "column T1 S -> D B stages=12 feed=6 R=2.5 D=40 P=1atm\n"
      "splitter SP B -> R W fractions=0.2,0.8\n",
      d);
  CHECK(spec.has_value());
  auto r = solve_flowsheet(*spec, d);
  CHECK(r.has_value());
  CHECK(r->converged);
  double ben_out = 0.0;
  for (const auto& s : r->streams)
    if (s.destination.empty()) ben_out += s.stream.n[0];
  CHECK_NEAR(ben_out, 40.0, 1e-5);
  CHECK_NEAR(r->energy_imbalance_W, 0.0, 1.0);
  bool profile = false;
  for (const auto& u : r->units) profile = profile || (u.kind == UnitKind::Column && u.profile.size() == 12);
  CHECK(profile);
}
