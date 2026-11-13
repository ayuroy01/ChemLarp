// Material streams and unit operations.
//
// Water reference values are from the IAPWS-95 formulation (Wagner & Pruß, J. Phys. Chem. Ref. Data 31
// (2002) 387, as implemented in the `chemicals` package), which is independent of the property data the
// model uses. Tolerances reflect the model (ideal gas, Antoine vapour pressure, DIPPR liquid Cp), not
// numerical precision; see "Unit operations: validation" in docs/ARCHITECTURE.md.
#include <cmath>
#include <vector>

#include "chemlarp/process.hpp"
#include "test_framework.hpp"

using namespace chemlarp;

namespace {
ProcessComponents comps(const std::vector<std::string>& smiles, LiquidModel m = LiquidModel::Ideal) {
  Diagnostics d;
  auto c = process_components(smiles, m, d);
  CHECK(c.has_value());
  return *c;
}
Stream tp(const ProcessComponents& c, double T, double P, const std::vector<double>& n) {
  Diagnostics d;
  auto s = stream_tp(c, T, P, n, d);
  CHECK(s.has_value());
  return *s;
}
}  // namespace

TEST(water_throttling_against_iapws95) {
  // Saturated liquid at P1 expanded to 1 atm: vapour fraction from IAPWS-95 enthalpies.
  const ProcessComponents c = comps({"O"});
  const double cases[][2] = {{2e5, 0.037956}, {5e5, 0.097953}, {10e5, 0.15221}, {20e5, 0.216905}};
  for (const auto& k : cases) {
    Diagnostics d;
    auto b = bubble_temperature(c.vle, k[0], {1.0}, d);
    CHECK(b.has_value());
    const Stream in = tp(c, b->T_K - 1e-6, k[0], {1.0});
    CHECK_EQ(in.vapor_fraction, 0.0);
    auto out = valve(c, in, 101325.0, d);
    CHECK(out.has_value());
    CHECK_NEAR(out->vapor_fraction, k[1], 0.002);
    CHECK_NEAR(out->H_W, in.H_W, 1e-6 * std::fabs(in.H_W));  // isenthalpic
  }
}

TEST(water_heating_and_mixing_against_iapws95) {
  const ProcessComponents c = comps({"O"});
  Diagnostics d;
  auto h = heater_to_temperature(c, tp(c, 300.0, 101325.0, {1.0}), 360.0, 0.0, d);
  CHECK(h.has_value());
  CHECK_NEAR(h->duty_W, 4524.84, 0.005 * 4524.84);  // liquid, J/mol
  auto m = mixer(c, {tp(c, 290.0, 101325.0, {1.0}), tp(c, 350.0, 101325.0, {1.0})}, 0.0, d);
  CHECK(m.has_value());
  CHECK_NEAR(m->T_K, 320.020, 0.05);
  auto s = heater_to_temperature(c, tp(c, 400.0, 101325.0, {1.0}), 600.0, 0.0, d);
  CHECK(s.has_value());
  CHECK_NEAR(s->duty_W, 7177.77, 0.02 * 7177.77);  // steam: ideal-gas Cp, −1.4%
  // Boiling a pure liquid with a given duty: lever rule at the boiling point.
  auto b = bubble_temperature(c.vle, 101325.0, {1.0}, d);
  auto boil = heater_with_duty(c, tp(c, b->T_K - 1e-6, 101325.0, {1.0}), 20000.0, 0.0, d);
  CHECK(boil.has_value());
  CHECK_NEAR(boil->out.vapor_fraction, 20000.0 / 40650.95, 0.002);  // IAPWS-95 ΔHvap at 1 atm
  CHECK_NEAR(boil->out.T_K, b->T_K, 1e-6);
}

TEST(raoult_flash_matches_rachford_rice) {
  // Benzene + toluene, ideal solution: the split at 368 K and 1 atm from an independent Rachford–Rice solve
  // with the same Antoine coefficients.
  const ProcessComponents c = comps({"c1ccccc1", "Cc1ccccc1"});
  const Stream s = tp(c, 368.0, 101325.0, {0.5, 0.5});
  CHECK_NEAR(s.vapor_fraction, 0.4168863776, 1e-9);
  CHECK_NEAR(s.x[0], 0.4073956743, 1e-9);
  CHECK_NEAR(s.y[0], 0.6295289238, 1e-9);
  CHECK_EQ(tp(c, 360.0, 101325.0, {0.5, 0.5}).vapor_fraction, 0.0);  // below the bubble point
  CHECK_EQ(tp(c, 375.0, 101325.0, {0.5, 0.5}).vapor_fraction, 1.0);  // above the dew point
}

TEST(unit_operations_conserve_mass_and_energy) {
  const ProcessComponents c = comps({"CCO", "O"}, LiquidModel::Unifac);
  Diagnostics d;
  const Stream feed = tp(c, 300.0, 101325.0, {0.3, 0.7});
  auto h = heater_to_temperature(c, feed, 360.0, 0.0, d);
  CHECK(h.has_value());
  CHECK(h->out.vapor_fraction > 0.0 && h->out.vapor_fraction < 1.0);
  CHECK_NEAR(feed.H_W + h->duty_W, h->out.H_W, 1e-9 * std::fabs(feed.H_W));
  // Heating by the same duty returns the same state (PH flash inverts the TP flash).
  auto h2 = heater_with_duty(c, feed, h->duty_W, 0.0, d);
  CHECK(h2.has_value());
  CHECK_NEAR(h2->out.T_K, 360.0, 1e-6);
  CHECK_NEAR(h2->out.vapor_fraction, h->out.vapor_fraction, 1e-8);
  // Adiabatic flash: phases carry all of each component and all of the enthalpy.
  auto fl = flash_drum_adiabatic(c, h->out, 50662.5, d);
  CHECK(fl.has_value());
  for (int i = 0; i < 2; ++i) CHECK_NEAR(fl->vapor.n[i] + fl->liquid.n[i], h->out.n[i], 1e-12);
  CHECK_NEAR(fl->vapor.H_W + fl->liquid.H_W, h->out.H_W, 1e-6 * std::fabs(h->out.H_W));
  CHECK(fl->vapor.y[0] > fl->liquid.x[0]);  // ethanol is enriched in the vapour
  // Remixing the two phases gives back the flashed state.
  auto m = mixer(c, {fl->vapor, fl->liquid}, 0.0, d);
  CHECK(m.has_value());
  CHECK_NEAR(m->T_K, fl->vapor.T_K, 1e-6);
  // Splitter: fractions of everything.
  auto sp = splitter(c, feed, {0.25, 0.75}, d);
  CHECK(sp.has_value());
  CHECK_NEAR((*sp)[0].n[0] + (*sp)[1].n[0], feed.n[0], 1e-15);
  CHECK_NEAR((*sp)[0].H_W, 0.25 * feed.H_W, 1e-9 * std::fabs(feed.H_W));
  // Flash drum at fixed T: duty closes the balance.
  auto fd = flash_drum(c, feed, 355.0, 101325.0, d);
  CHECK(fd.has_value());
  CHECK_NEAR(feed.H_W + fd->duty_W, fd->vapor.H_W + fd->liquid.H_W, 1e-6 * std::fabs(feed.H_W));
  CHECK(!d.has_errors());
}

TEST(process_input_errors) {
  Diagnostics d;
  CHECK(!process_components({"CCO", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"}, LiquidModel::Ideal, d));
  CHECK(d.contains(DiagCode::NoPropertyData));
  const ProcessComponents c = comps({"O"});
  Diagnostics e;
  CHECK(!stream_tp(c, 300.0, 101325.0, {0.0}, e));
  CHECK(!stream_tp(c, -1.0, 101325.0, {1.0}, e));
  CHECK(!stream_tp(c, 300.0, 101325.0, {1.0, 1.0}, e));
  CHECK(!valve(c, tp(c, 300.0, 1e5, {1.0}), 2e5, e));  // a valve cannot raise the pressure
  CHECK(!splitter(c, tp(c, 300.0, 1e5, {1.0}), {0.5, 0.6}, e));
  CHECK_EQ(e.count(DiagCode::InvalidMixtureInput), std::size_t{5});
}
