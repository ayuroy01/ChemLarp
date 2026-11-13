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

TEST(ph_flash_from_the_exact_temperature) {
  // Regression: a PH flash started at its own solution (as a mixer fed with the phases of one flash) left a
  // collapsed bracket and divided 0 by 0 on platforms where the residual rounded to exactly zero.
  const ProcessComponents c = comps({"CCO", "O"}, LiquidModel::Unifac);
  for (double T : {300.0, 355.0, 380.0}) {  // liquid, two-phase, vapour
    const Stream s = tp(c, T, 101325.0, {0.3, 0.7});
    for (double guess : {T, T - 50.0, T + 50.0, 0.0}) {
      Diagnostics d;
      auto r = stream_ph(c, 101325.0, s.H_W, s.n, d, guess);
      CHECK(r.has_value());
      CHECK_NEAR(r->T_K, T, 1e-6);
      CHECK_NEAR(r->vapor_fraction, s.vapor_fraction, 1e-8);
    }
  }
}

TEST(pump_and_compressor) {
  const ProcessComponents w = comps({"O"});
  Diagnostics d;
  // Water 300 K, 1 → 100 bar: isentropic work 178.58 J/mol and, at 75% efficiency, 300.996 K (IAPWS-95).
  auto p = pump(w, tp(w, 300.0, 1e5, {1.0}), 100e5, 1.0, d);
  CHECK(p.has_value());
  CHECK_NEAR(p->work_W, 178.58, 0.01 * 178.58);
  auto p75 = pump(w, tp(w, 300.0, 1e5, {1.0}), 100e5, 0.75, d);
  CHECK_NEAR(p75->work_W, p75->isentropic_work_W / 0.75, 1e-9);
  CHECK_NEAR(p75->out.T_K, 300.996, 0.25);  // the model omits the 0.2 K isentropic temperature rise
  // Argon (Cp = 5R/2): T2 = T1 (P2/P1)^(2/5) exactly for an ideal gas.
  const ProcessComponents ar = comps({"[Ar]"});
  auto a = compressor(ar, tp(ar, 300.0, 1e5, {1.0}), 5e5, 1.0, d);
  CHECK(a.has_value());
  CHECK_NEAR(a->out.T_K, 300.0 * std::pow(5.0, 0.4), 0.05);
  // Steam 450 K, 0.1 → 0.5 bar: IAPWS-95 isentropic work 7313.3 J/mol, outlet 654.40 K; 80%: 703.42 K.
  auto s = compressor(w, tp(w, 450.0, 0.1e5, {1.0}), 0.5e5, 1.0, d);
  CHECK(s.has_value());
  CHECK_NEAR(s->work_W, 7313.3, 0.005 * 7313.3);
  CHECK_NEAR(s->out.T_K, 654.40, 1.5);
  auto s80 = compressor(w, tp(w, 450.0, 0.1e5, {1.0}), 0.5e5, 0.8, d);
  CHECK_NEAR(s80->out.T_K, 703.42, 1.5);
  // Expansion back: work is negative and the efficiency reduces what is recovered.
  auto x = compressor(w, s->out, 0.1e5, 1.0, d);
  CHECK(x.has_value());
  CHECK_NEAR(x->work_W, -s->work_W, 1e-6 * s->work_W);
  CHECK(!d.has_errors());
  Diagnostics e;
  CHECK(!pump(w, tp(w, 400.0, 1e5, {1.0}), 2e5, 1.0, e));        // vapour feed
  CHECK(!compressor(w, tp(w, 300.0, 1e5, {1.0}), 2e5, 1.0, e));  // liquid feed
  CHECK_EQ(e.count(DiagCode::InvalidMixtureInput), std::size_t{2});
}

TEST(heat_exchanger_profiles) {
  const ProcessComponents w = comps({"O"});
  const Stream hot = tp(w, 360.0, 101325.0, {1.0}), cold = tp(w, 290.0, 101325.0, {1.0});
  Diagnostics d;
  // Balanced liquid water: IAPWS-95 gives Q = 3018.9 J/mol and a cold outlet of 330.079 K; the approach is
  // nearly constant, so UA ≈ Q / 30 K.
  auto hx = heat_exchanger(w, hot, w, cold, ExchangerSpec::HotOutletTemperature, 320.0, d);
  CHECK(hx.has_value());
  CHECK_NEAR(hx->duty_W, 3018.9, 0.003 * 3018.9);
  CHECK_NEAR(hx->cold_out.T_K, 330.079, 0.1);
  CHECK_NEAR(hx->UA_WK, hx->duty_W / 30.0, 0.01 * hx->duty_W / 30.0);
  CHECK_NEAR(hot.H_W + cold.H_W, hx->hot_out.H_W + hx->cold_out.H_W, 1e-6 * std::fabs(hot.H_W));
  // The same exchanger by duty and by cold outlet temperature.
  auto byq = heat_exchanger(w, hot, w, cold, ExchangerSpec::Duty, hx->duty_W, d);
  CHECK_NEAR(byq->hot_out.T_K, 320.0, 1e-6);
  auto byc = heat_exchanger(w, hot, w, cold, ExchangerSpec::ColdOutletTemperature, hx->cold_out.T_K, d);
  CHECK_NEAR(byc->hot_out.T_K, 320.0, 1e-6);
  // Condensing steam against water: the pinch is inside the exchanger, at the dew point.
  const Stream steam = tp(w, 400.0, 101325.0, {1.0});
  auto cond = heat_exchanger(w, steam, w, tp(w, 290.0, 101325.0, {10.0}), ExchangerSpec::HotOutletTemperature, 360.0, d);
  CHECK(cond.has_value());
  CHECK(cond->min_approach_K > 0.0 && cond->min_approach_K < 400.0 - cond->cold_out.T_K);
  CHECK(!d.has_errors());
  // Asking for a cold outlet above the hot inlet is a temperature cross.
  Diagnostics e;
  CHECK(!heat_exchanger(w, hot, w, cold, ExchangerSpec::ColdOutletTemperature, 365.0, e));
  CHECK(e.contains(DiagCode::TemperatureCross));
}

TEST(conversion_reactor_heats_of_reaction) {
  // CH4 + 2 O2 -> CO2 + 2 H2O at 298.15 K (CO2 has no vapour-pressure equation: non-condensable).
  const ProcessComponents g = comps({"C", "O=O", "O=C=O", "O"});
  Diagnostics d;
  // All gas (low pressure): Σ ν ΔfH° = −802.5 kJ/mol from the CRC values (CODATA −802.3).
  auto lo = conversion_reactor(g, tp(g, 298.15, 1000.0, {1, 2, 0, 0}), {-1, -2, 1, 2}, 0, 1.0, 298.15, 0.0, d);
  CHECK(lo.has_value());
  CHECK_EQ(lo->out.vapor_fraction, 1.0);
  CHECK_NEAR(lo->duty_W, -802.3e3, 0.001 * 802.3e3);
  CHECK_NEAR(lo->duty_W, lo->heat_of_reaction_Jmol, 1e-6);
  // At 10 bar the water condenses (but for 0.32% in the CO2): higher heating value 890.56 kJ/mol (NIST),
  // corrected for that vapour.
  auto hi = conversion_reactor(g, tp(g, 298.15, 10e5, {1, 2, 0, 0}), {-1, -2, 1, 2}, 0, 1.0, 298.15, 0.0, d);
  CHECK(hi.has_value());
  CHECK_NEAR(hi->out.x[2], 0.0, 1e-15);  // CO2 stays in the vapour
  const double vapour_water = hi->out.vapor_fraction * hi->out.total() * hi->out.y[3];
  CHECK_NEAR(hi->duty_W - vapour_water * 43990.0, -890.56e3, 0.001 * 890.56e3);
  // Adiabatic: the duty is zero and the outlet hot.
  auto ad = conversion_reactor(g, tp(g, 298.15, 1e5, {1, 10, 0, 0}), {-1, -2, 1, 2}, 0, 1.0, 0.0, 0.0, d);
  CHECK(ad.has_value());
  CHECK_NEAR(ad->duty_W, 0.0, 1e-6);
  CHECK(ad->out.T_K > 1500.0);
  CHECK(!d.has_errors());
  Diagnostics e;
  CHECK(!conversion_reactor(g, tp(g, 298.15, 1e5, {1, 2, 0, 0}), {-1, -1, 1, 2}, 0, 1.0, 298.15, 0.0, e));
  CHECK(e.contains(DiagCode::UnbalancedReaction));
  Diagnostics f;  // not enough oxygen for full conversion
  CHECK(!conversion_reactor(g, tp(g, 298.15, 1e5, {1, 1, 0, 0}), {-1, -2, 1, 2}, 0, 1.0, 298.15, 0.0, f));
  CHECK(f.contains(DiagCode::InvalidMixtureInput));
}
