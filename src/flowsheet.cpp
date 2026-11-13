#include "chemlarp/flowsheet.hpp"

#include "chemlarp/json.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <set>
#include <sstream>

namespace chemlarp {

const char* to_string(UnitKind kind) noexcept {
  switch (kind) {
    case UnitKind::Mixer: return "mixer";
    case UnitKind::Splitter: return "splitter";
    case UnitKind::Heater: return "heater";
    case UnitKind::Valve: return "valve";
    case UnitKind::Flash: return "flash";
    case UnitKind::Separator: return "separator";
    case UnitKind::Pump: return "pump";
    case UnitKind::Compressor: return "compressor";
    case UnitKind::Exchanger: return "exchanger";
    case UnitKind::Reactor: return "reactor";
  }
  return "unit";
}

namespace {

// ---------------------------------------------------------------------------------------------------------
// Parsing

std::string strip_comment(const std::string& line) {
  for (std::size_t i = 0; i < line.size(); ++i)
    if (line[i] == '#' && (i == 0 || std::isspace(static_cast<unsigned char>(line[i - 1])))) return line.substr(0, i);
  return line;
}

std::vector<std::string> split_ws(const std::string& s) {
  std::istringstream in(s);
  std::vector<std::string> out;
  for (std::string t; in >> t;) out.push_back(t);
  return out;
}

bool to_number(const std::string& s, double& v) {
  if (s.empty()) return false;
  char* end = nullptr;
  v = std::strtod(s.c_str(), &end);
  return end == s.c_str() + s.size() && std::isfinite(v);
}

// A number with an optional unit: dim 'T' (K, C), 'P' (Pa, kPa, MPa, bar, atm), 'Q' (W, kW, MW).
bool to_quantity(const std::string& s, char dim, double& v) {
  if (s.empty()) return false;
  char* end = nullptr;
  v = std::strtod(s.c_str(), &end);
  if (end == s.c_str() || !std::isfinite(v)) return false;
  const std::string u(end);
  if (u.empty()) return true;
  if (dim == 'T') {
    if (u == "K") return true;
    if (u == "C") return v += 273.15, true;
  } else if (dim == 'P') {
    if (u == "Pa") return true;
    if (u == "kPa") return v *= 1e3, true;
    if (u == "MPa") return v *= 1e6, true;
    if (u == "bar") return v *= 1e5, true;
    if (u == "atm") return v *= 101325.0, true;
  } else if (dim == 'Q') {
    if (u == "W") return true;
    if (u == "kW") return v *= 1e3, true;
    if (u == "MW") return v *= 1e6, true;
  }
  return false;
}

// key=value, splitting at the last '=' (SMILES keys may contain '=').
bool key_value(const std::string& tok, std::string& k, std::string& v) {
  const auto p = tok.rfind('=');
  if (p == std::string::npos || p == 0 || p + 1 == tok.size()) return false;
  k = tok.substr(0, p);
  v = tok.substr(p + 1);
  return true;
}

struct Parser {
  FlowsheetSpec spec;
  Diagnostics& diag;
  int line = 0;
  bool ok = true;

  bool fail(const std::string& msg) {
    diag.error(DiagCode::FlowsheetSyntax, "line " + std::to_string(line) + ": " + msg);
    ok = false;
    return false;
  }

  std::optional<std::size_t> component(const std::string& s) {
    for (std::size_t i = 0; i < spec.components.size(); ++i)
      if (spec.components[i] == s) return i;
    if (spec.components.empty()) fail("declare the components before using them");
    else fail("'" + s + "' is not one of the components");
    return std::nullopt;
  }

  // Component values after a section keyword (flows, split, nu); unlisted components are 0.
  bool component_values(const std::vector<std::string>& t, std::size_t from, std::vector<double>& out) {
    out.assign(spec.components.size(), 0.0);
    std::vector<bool> seen(out.size(), false);
    for (std::size_t i = from; i < t.size(); ++i) {
      std::string k, v;
      double x = 0.0;
      if (!key_value(t[i], k, v) || !to_number(v, x)) return fail("expected <SMILES>=<number>, got '" + t[i] + "'");
      auto c = component(k);
      if (!c) return false;
      if (seen[*c]) return fail(k + " is given twice");
      seen[*c] = true;
      out[*c] = x;
    }
    return true;
  }

  bool stream_spec(const std::vector<std::string>& t, FeedSpec& f) {
    if (t.size() < 2) return fail(t[0] + " needs a stream name");
    f.name = t[1];
    std::size_t i = 2;
    for (; i < t.size() && t[i] != "flows"; ++i) {
      std::string k, v;
      if (!key_value(t[i], k, v)) return fail("expected T=, P= or 'flows', got '" + t[i] + "'");
      if (k == "T") {
        if (!to_quantity(v, 'T', f.T_K) || !(f.T_K > 0)) return fail("bad temperature '" + v + "'");
      } else if (k == "P") {
        if (!to_quantity(v, 'P', f.P_Pa) || !(f.P_Pa > 0)) return fail("bad pressure '" + v + "'");
      } else {
        return fail("unknown parameter '" + k + "' (expected T, P, then 'flows')");
      }
    }
    if (!(f.T_K > 0) || !(f.P_Pa > 0)) return fail(t[0] + " " + f.name + " needs T= and P=");
    if (i == t.size()) return fail(t[0] + " " + f.name + " needs 'flows <SMILES>=<mol/s> ...'");
    if (!component_values(t, i + 1, f.n)) return false;
    for (double v : f.n)
      if (v < 0) return fail("flows must be non-negative");
    return true;
  }

  bool unit(const std::vector<std::string>& t, UnitKind kind) {
    UnitSpec u;
    u.kind = kind;
    u.line = line;
    if (t.size() < 2) return fail(t[0] + " needs a name");
    u.name = t[1];
    for (const auto& o : spec.units)
      if (o.name == u.name) return fail("a unit named '" + u.name + "' already exists");
    std::size_t i = 2;
    for (; i < t.size() && t[i] != "->"; ++i) {
      if (t[i].find('=') != std::string::npos) return fail("expected inlet streams, then '->'");
      u.inlets.push_back(t[i]);
    }
    if (i == t.size()) return fail("expected '<inlets> -> <outlets>'");
    for (++i; i < t.size() && t[i].find('=') == std::string::npos && t[i] != "split" && t[i] != "nu"; ++i)
      u.outlets.push_back(t[i]);
    bool hasT = false, hasP = false, hasX = false, hasKey = false, hasEff = false, hasTh = false, hasTc = false;
    std::string key_name;
    for (; i < t.size() && t[i] != "split" && t[i] != "nu"; ++i) {
      std::string k, v;
      if (!key_value(t[i], k, v)) return fail("expected key=value, got '" + t[i] + "'");
      bool good = true;
      if (k == "T") good = hasT = to_quantity(v, 'T', u.T_K) && u.T_K > 0;
      else if (k == "P") good = hasP = to_quantity(v, 'P', u.P_Pa) && u.P_Pa > 0;
      else if (k == "Q") good = u.has_duty = to_quantity(v, 'Q', u.duty_W);
      else if (k == "Thot") good = hasTh = to_quantity(v, 'T', u.hot_T_K) && u.hot_T_K > 0;
      else if (k == "Tcold") good = hasTc = to_quantity(v, 'T', u.cold_T_K) && u.cold_T_K > 0;
      else if (k == "eff") good = hasEff = to_number(v, u.efficiency) && u.efficiency > 0 && u.efficiency <= 1;
      else if (k == "X") good = hasX = to_number(v, u.conversion) && u.conversion >= 0 && u.conversion <= 1;
      else if (k == "key") key_name = v, good = hasKey = true;
      else if (k == "fractions") {
        std::string item;
        std::istringstream in(v);
        while (std::getline(in, item, ',')) {
          double f = 0;
          if (!to_number(item, f) || f < 0) return fail("bad split fraction '" + item + "'");
          u.fractions.push_back(f);
        }
      } else {
        return fail("unknown parameter '" + k + "' for " + to_string(kind) + " " + u.name);
      }
      if (!good) return fail("bad value for " + k + ": '" + v + "'");
    }
    if (i < t.size()) {
      const bool split = t[i] == "split";
      if (split != (kind == UnitKind::Separator)) return fail("'" + t[i] + "' does not belong to a " + to_string(kind));
      if (!component_values(t, i + 1, u.per_component)) return false;
    }
    const std::string what = std::string(to_string(kind)) + " " + u.name;
    auto ports = [&](std::size_t nin, std::size_t nout) {
      // nin/nout == 0: at least one / two
      const bool in_ok = nin ? u.inlets.size() == nin : !u.inlets.empty();
      const bool out_ok = nout ? u.outlets.size() == nout : u.outlets.size() >= 2;
      if (in_ok && out_ok) return true;
      return fail(what + " has the wrong number of streams (" + std::to_string(u.inlets.size()) + " in, " +
                  std::to_string(u.outlets.size()) + " out)");
    };
    switch (kind) {
      case UnitKind::Mixer:
        if (!ports(0, 1)) return false;
        break;
      case UnitKind::Splitter:
        if (!ports(1, 0)) return false;
        if (u.fractions.size() != u.outlets.size())
          return fail(what + " needs fractions= with one value per outlet");
        break;
      case UnitKind::Heater:
        if (!ports(1, 1)) return false;
        if (hasT == u.has_duty) return fail(what + " needs either T= or Q=");
        break;
      case UnitKind::Valve:
        if (!ports(1, 1)) return false;
        if (!hasP) return fail(what + " needs P=");
        break;
      case UnitKind::Flash:
        if (!ports(1, 2)) return false;
        if (!hasP) return fail(what + " needs P= (and T=, or none for an adiabatic flash)");
        break;
      case UnitKind::Separator:
        if (!ports(1, 2)) return false;
        if (u.per_component.empty()) return fail(what + " needs 'split <SMILES>=<fraction to first outlet> ...'");
        for (double f : u.per_component)
          if (f < 0 || f > 1) return fail("separator split fractions must be in [0, 1]");
        break;
      case UnitKind::Pump:
      case UnitKind::Compressor:
        if (!ports(1, 1)) return false;
        if (!hasP) return fail(what + " needs P=");
        if (!hasEff) u.efficiency = 1.0;
        break;
      case UnitKind::Exchanger:
        if (!ports(2, 2)) return false;
        if (int(hasTh) + int(hasTc) + int(u.has_duty) != 1) return fail(what + " needs one of Thot=, Tcold= or Q=");
        break;
      case UnitKind::Reactor: {
        if (!ports(1, 1)) return false;
        if (!hasKey || !hasX || u.per_component.empty())
          return fail(what + " needs key=<SMILES>, X=<conversion> and 'nu <SMILES>=<coefficient> ...'");
        auto k = component(key_name);
        if (!k) return false;
        u.key = *k;
        if (!(u.per_component[u.key] < 0)) return fail("the key component of " + what + " must be a reactant (nu < 0)");
        break;
      }
    }
    spec.units.push_back(std::move(u));
    return true;
  }

  void statement(const std::vector<std::string>& t) {
    std::string kw = t[0];
    std::transform(kw.begin(), kw.end(), kw.begin(), [](unsigned char ch) { return std::tolower(ch); });
    static const std::map<std::string, UnitKind> kinds = {
        {"mixer", UnitKind::Mixer},         {"splitter", UnitKind::Splitter},   {"heater", UnitKind::Heater},
        {"cooler", UnitKind::Heater},       {"valve", UnitKind::Valve},         {"flash", UnitKind::Flash},
        {"separator", UnitKind::Separator}, {"sep", UnitKind::Separator},       {"pump", UnitKind::Pump},
        {"compressor", UnitKind::Compressor}, {"expander", UnitKind::Compressor}, {"turbine", UnitKind::Compressor},
        {"exchanger", UnitKind::Exchanger}, {"reactor", UnitKind::Reactor}};
    if (kw == "components") {
      if (!spec.components.empty()) fail("components are declared twice");
      else if (t.size() < 2) fail("components needs at least one SMILES");
      else spec.components.assign(t.begin() + 1, t.end());
      for (std::size_t i = 0; i < spec.components.size(); ++i)
        for (std::size_t j = 0; j < i; ++j)
          if (spec.components[i] == spec.components[j]) fail(spec.components[i] + " is listed twice");
    } else if (kw == "model") {
      if (t.size() != 2 || (t[1] != "ideal" && t[1] != "unifac")) fail("model must be 'ideal' or 'unifac'");
      else spec.model = t[1] == "unifac" ? LiquidModel::Unifac : LiquidModel::Ideal;
    } else if (kw == "feed" || kw == "guess") {
      FeedSpec f;
      if (stream_spec(t, f)) (kw == "feed" ? spec.feeds : spec.guesses).push_back(std::move(f));
    } else if (kw == "tear") {
      if (t.size() < 2) fail("tear needs stream names");
      spec.tears.insert(spec.tears.end(), t.begin() + 1, t.end());
    } else if (kw == "method") {
      if (t.size() != 2 || (t[1] != "wegstein" && t[1] != "direct")) fail("method must be 'wegstein' or 'direct'");
      else spec.method = t[1] == "direct" ? TearMethod::Direct : TearMethod::Wegstein;
    } else if (kw == "tolerance") {
      if (t.size() != 2 || !to_number(t[1], spec.tolerance) || !(spec.tolerance > 0 && spec.tolerance < 0.1))
        fail("tolerance must be a number between 0 and 0.1");
    } else if (kw == "max-iterations") {
      double v = 0;
      if (t.size() != 2 || !to_number(t[1], v) || v < 1 || v > 10000 || v != std::floor(v))
        fail("max-iterations must be a whole number from 1 to 10000");
      else spec.max_iterations = static_cast<int>(v);
    } else if (auto it = kinds.find(kw); it != kinds.end()) {
      unit(t, it->second);
    } else {
      fail("unknown statement '" + t[0] + "'");
    }
  }
};

// ---------------------------------------------------------------------------------------------------------
// Solving

Stream empty_stream(std::size_t nc, double T, double P) {
  Stream s;
  s.T_K = T;
  s.P_Pa = P;
  s.n.assign(nc, 0.0);
  s.x = s.y = s.n;
  return s;
}

bool is_empty(const Stream& s) { return !(s.total() > 0.0); }

void merge_prefixed(const Diagnostics& from, Diagnostics& to, const std::string& prefix) {
  for (const auto& it : from.items()) {
    const std::string msg = prefix + it.message;
    bool seen = false;
    for (const auto& t : to.items()) seen = seen || (t.code == it.code && t.message == msg);
    if (seen) continue;
    if (it.severity == Severity::Error) to.error(it.code, msg);
    else to.warn(it.code, msg);
  }
}

// One unit calculation. Inlets without flow (a tear stream before the first pass) give empty outlets.
std::optional<std::vector<Stream>> run_unit(const ProcessComponents& c, const UnitSpec& u, const std::vector<Stream>& in,
                                            UnitReport& rep, Diagnostics& d) {
  const std::size_t nc = c.size();
  rep = UnitReport{u.name, u.kind, 0.0, 0.0, {}};
  const Stream& a = in[0];
  switch (u.kind) {
    case UnitKind::Mixer: {
      std::vector<Stream> full;
      for (const auto& s : in)
        if (!is_empty(s)) full.push_back(s);
      if (full.empty()) return std::vector<Stream>{empty_stream(nc, a.T_K, u.P_Pa > 0 ? u.P_Pa : a.P_Pa)};
      auto o = mixer(c, full, u.P_Pa, d);
      if (!o) return std::nullopt;
      return std::vector<Stream>{*o};
    }
    case UnitKind::Splitter: {
      if (is_empty(a)) return std::vector<Stream>(u.outlets.size(), a);
      double sum = 0;
      for (double f : u.fractions) sum += f;
      std::vector<double> f = u.fractions;
      if (std::fabs(sum - 1.0) > 1e-6) {
        d.error(DiagCode::InvalidMixtureInput, "split fractions must sum to 1");
        return std::nullopt;
      }
      for (double& v : f) v /= sum;
      return splitter(c, a, f, d);
    }
    case UnitKind::Heater: {
      if (is_empty(a)) return std::vector<Stream>{empty_stream(nc, u.T_K > 0 ? u.T_K : a.T_K, u.P_Pa > 0 ? u.P_Pa : a.P_Pa)};
      auto r = u.has_duty ? heater_with_duty(c, a, u.duty_W, u.P_Pa, d) : heater_to_temperature(c, a, u.T_K, u.P_Pa, d);
      if (!r) return std::nullopt;
      rep.duty_W = r->duty_W;
      rep.values.push_back({"outlet_T_K", r->out.T_K});
      return std::vector<Stream>{r->out};
    }
    case UnitKind::Valve: {
      if (is_empty(a)) return std::vector<Stream>{empty_stream(nc, a.T_K, u.P_Pa)};
      auto o = valve(c, a, u.P_Pa, d);
      if (!o) return std::nullopt;
      rep.values.push_back({"outlet_T_K", o->T_K});
      return std::vector<Stream>{*o};
    }
    case UnitKind::Flash: {
      if (is_empty(a)) {
        Stream e = empty_stream(nc, u.T_K > 0 ? u.T_K : a.T_K, u.P_Pa);
        return std::vector<Stream>{e, e};
      }
      auto r = u.T_K > 0 ? flash_drum(c, a, u.T_K, u.P_Pa, d) : flash_drum_adiabatic(c, a, u.P_Pa, d);
      if (!r) return std::nullopt;
      rep.duty_W = r->duty_W;
      const double F = a.total();
      rep.values.push_back({"T_K", r->vapor.T_K});
      rep.values.push_back({"vapor_fraction", F > 0 ? r->vapor.total() / F : 0.0});
      // An outlet without flow is reported empty.
      Stream v = r->vapor, l = r->liquid;
      if (!(v.total() > 0)) v = empty_stream(nc, v.T_K, v.P_Pa);
      if (!(l.total() > 0)) l = empty_stream(nc, l.T_K, l.P_Pa);
      return std::vector<Stream>{v, l};
    }
    case UnitKind::Separator: {
      auto r = component_separator(c, a, u.per_component, u.T_K, u.P_Pa, d);
      if (!r) return std::nullopt;
      rep.duty_W = r->duty_W;
      return std::vector<Stream>{r->first, r->second};
    }
    case UnitKind::Pump:
    case UnitKind::Compressor: {
      if (is_empty(a)) return std::vector<Stream>{empty_stream(nc, a.T_K, u.P_Pa)};
      auto r = u.kind == UnitKind::Pump ? pump(c, a, u.P_Pa, u.efficiency, d) : compressor(c, a, u.P_Pa, u.efficiency, d);
      if (!r) return std::nullopt;
      rep.work_W = r->work_W;
      rep.values.push_back({"isentropic_work_W", r->isentropic_work_W});
      rep.values.push_back({"outlet_T_K", r->out.T_K});
      return std::vector<Stream>{r->out};
    }
    case UnitKind::Exchanger: {
      const Stream& b = in[1];
      if (is_empty(a) || is_empty(b)) {
        rep.values.push_back({"exchanger_duty_W", 0.0});
        return std::vector<Stream>{a, b};
      }
      ExchangerSpec spec = ExchangerSpec::Duty;
      double value = u.duty_W;
      if (u.hot_T_K > 0) spec = ExchangerSpec::HotOutletTemperature, value = u.hot_T_K;
      if (u.cold_T_K > 0) spec = ExchangerSpec::ColdOutletTemperature, value = u.cold_T_K;
      auto r = heat_exchanger(c, a, c, b, spec, value, d);
      if (!r) return std::nullopt;
      rep.values.push_back({"exchanger_duty_W", r->duty_W});
      rep.values.push_back({"min_approach_K", r->min_approach_K});
      rep.values.push_back({"UA_WK", r->UA_WK});
      return std::vector<Stream>{r->hot_out, r->cold_out};
    }
    case UnitKind::Reactor: {
      if (is_empty(a)) return std::vector<Stream>{empty_stream(nc, u.T_K > 0 ? u.T_K : a.T_K, u.P_Pa > 0 ? u.P_Pa : a.P_Pa)};
      auto r = conversion_reactor(c, a, u.per_component, u.key, u.conversion, u.T_K, u.P_Pa, d);
      if (!r) return std::nullopt;
      rep.duty_W = r->duty_W;
      rep.values.push_back({"extent_mols", r->extent_mols});
      rep.values.push_back({"heat_of_reaction_Jmol", r->heat_of_reaction_Jmol});
      rep.values.push_back({"outlet_T_K", r->out.T_K});
      return std::vector<Stream>{r->out};
    }
  }
  return std::nullopt;
}

// Relative change between two values of a tear stream: component flows against the total flow, enthalpy
// against max(|H|, F × 10 kJ/mol).
double tear_change(const Stream& x, const Stream& g) {
  const double F = std::max(x.total(), g.total());
  if (!(F > 0)) return 0.0;
  double r = 0.0;
  for (std::size_t i = 0; i < x.n.size(); ++i) r = std::max(r, std::fabs(g.n[i] - x.n[i]) / F);
  const double Hs = std::max({std::fabs(x.H_W), std::fabs(g.H_W), F * 1e4});
  return std::max(r, std::fabs(g.H_W - x.H_W) / Hs);
}

}  // namespace

std::optional<FlowsheetSpec> parse_flowsheet(std::string_view text, Diagnostics& diag) {
  Parser p{{}, diag};
  std::string all(text), raw;
  std::istringstream in(all);
  while (std::getline(in, raw)) {
    ++p.line;
    auto t = split_ws(strip_comment(raw));
    if (!t.empty()) p.statement(t);
  }
  if (p.ok && p.spec.components.empty()) p.fail("no 'components' line");
  if (p.ok && p.spec.units.empty()) p.fail("no units");
  if (!p.ok) return std::nullopt;
  return p.spec;
}

std::optional<FlowsheetResult> solve_flowsheet(const FlowsheetSpec& spec, Diagnostics& diag) {
  FlowsheetResult res;
  {
    auto c = process_components(spec.components, spec.model, diag);
    if (!c) return std::nullopt;
    res.components = std::move(*c);
  }
  const ProcessComponents& c = res.components;
  const std::size_t nc = c.size(), nu = spec.units.size();
  auto topo = [&](const std::string& m) {
    diag.error(DiagCode::FlowsheetTopology, m);
    return std::nullopt;
  };

  // Streams: source and destination.
  struct Port {
    int unit = -1;  // -1: feed
    std::size_t index = 0;
  };
  std::map<std::string, Port> source;
  std::map<std::string, Port> dest;
  std::vector<std::string> names;  // feeds, then unit outlets
  for (std::size_t f = 0; f < spec.feeds.size(); ++f) {
    const auto& fd = spec.feeds[f];
    if (fd.n.size() != nc) return topo("feed " + fd.name + " needs one flow per component");
    if (!source.emplace(fd.name, Port{-1, f}).second) return topo("stream " + fd.name + " is defined twice");
    names.push_back(fd.name);
  }
  for (std::size_t u = 0; u < nu; ++u)
    for (std::size_t k = 0; k < spec.units[u].outlets.size(); ++k) {
      const auto& s = spec.units[u].outlets[k];
      if (!source.emplace(s, Port{int(u), k}).second)
        return topo("stream " + s + " has two sources (" + spec.units[u].name + " and another)");
      names.push_back(s);
    }
  for (std::size_t u = 0; u < nu; ++u)
    for (std::size_t k = 0; k < spec.units[u].inlets.size(); ++k) {
      const auto& s = spec.units[u].inlets[k];
      if (!source.count(s)) return topo("stream " + s + " (inlet of " + spec.units[u].name + ") has no source: add a feed or a unit outlet");
      if (!dest.emplace(s, Port{int(u), k}).second) return topo("stream " + s + " goes to two units");
    }
  for (const auto& f : spec.feeds)
    if (!dest.count(f.name)) diag.warn(DiagCode::FlowsheetTopology, "feed " + f.name + " goes to no unit");

  // Tear streams: the user's, then back edges of a depth-first search.
  std::set<std::string> tears;
  auto check_tear = [&](const std::string& s) {
    auto it = source.find(s);
    if (it == source.end() || it->second.unit < 0 || !dest.count(s)) {
      diag.error(DiagCode::FlowsheetTopology, "tear stream " + s + " must run from one unit to another");
      return false;
    }
    return true;
  };
  for (const auto& s : spec.tears) {
    if (!check_tear(s)) return std::nullopt;
    if (tears.insert(s).second) res.tears.push_back(s);
  }
  for (const auto& g : spec.guesses) {
    if (!check_tear(g.name)) return std::nullopt;
    if (g.n.size() != nc) return topo("guess for " + g.name + " needs one flow per component");
    if (tears.insert(g.name).second) res.tears.push_back(g.name);
  }
  std::vector<int> color(nu, 0);
  std::function<void(std::size_t)> dfs = [&](std::size_t u) {
    color[u] = 1;
    for (const auto& s : spec.units[u].outlets) {
      if (tears.count(s) || !dest.count(s)) continue;
      const std::size_t v = static_cast<std::size_t>(dest[s].unit);
      if (color[v] == 1) {
        tears.insert(s);
        res.tears.push_back(s);
      } else if (color[v] == 0) {
        dfs(v);
      }
    }
    color[u] = 2;
  };
  for (std::size_t u = 0; u < nu; ++u)
    if (!color[u]) dfs(u);

  // Calculation order: topological order of the graph without tear streams, earliest-declared first.
  std::vector<int> indeg(nu, 0);
  for (std::size_t u = 0; u < nu; ++u)
    for (const auto& s : spec.units[u].inlets)
      if (source[s].unit >= 0 && !tears.count(s)) ++indeg[u];
  std::vector<std::size_t> order;
  std::vector<bool> done(nu, false);
  while (order.size() < nu) {
    std::size_t pick = nu;
    for (std::size_t u = 0; u < nu && pick == nu; ++u)
      if (!done[u] && indeg[u] == 0) pick = u;
    if (pick == nu) return topo("the units could not be ordered");  // not reached: back edges are torn
    done[pick] = true;
    order.push_back(pick);
    for (const auto& s : spec.units[pick].outlets)
      if (dest.count(s) && !tears.count(s)) --indeg[static_cast<std::size_t>(dest[s].unit)];
  }
  for (std::size_t u : order) res.order.push_back(spec.units[u].name);

  // Stream values.
  std::map<std::string, Stream> val;
  for (const auto& f : spec.feeds) {
    double F = 0;
    for (double v : f.n) F += v;
    if (!(F > 0)) {
      val[f.name] = empty_stream(nc, f.T_K, f.P_Pa);
      continue;
    }
    Diagnostics q;
    auto s = stream_tp(c, f.T_K, f.P_Pa, f.n, q);
    merge_prefixed(q, diag, "feed " + f.name + ": ");
    if (!s) return std::nullopt;
    val[f.name] = *s;
  }
  std::vector<Stream> x;  // tear values
  std::vector<bool> guessed;
  for (const auto& t : res.tears) {
    bool g_found = false;
    for (const auto& g : spec.guesses) g_found = g_found || g.name == t;
    guessed.push_back(g_found);
    Stream s = empty_stream(nc, 298.15, 101325.0);
    for (const auto& g : spec.guesses)
      if (g.name == t) {
        Diagnostics q;
        double F = 0;
        for (double v : g.n) F += v;
        if (F > 0) {
          auto gs = stream_tp(c, g.T_K, g.P_Pa, g.n, q);
          merge_prefixed(q, diag, "guess " + g.name + ": ");
          if (!gs) return std::nullopt;
          s = *gs;
        }
      }
    x.push_back(s);
  }

  std::vector<UnitReport> reports(nu);
  Diagnostics last;  // warnings of the latest pass
  auto pass = [&]() -> bool {
    last = Diagnostics{};
    for (std::size_t k = 0; k < res.tears.size(); ++k) val[res.tears[k]] = x[k];
    for (std::size_t u : order) {
      const UnitSpec& us = spec.units[u];
      std::vector<Stream> in;
      for (const auto& s : us.inlets) in.push_back(val[s]);
      Diagnostics q;
      auto out = run_unit(c, us, in, reports[u], q);
      merge_prefixed(q, last, us.name + ": ");
      if (!out) return false;
      for (std::size_t k = 0; k < us.outlets.size(); ++k) {
        // A tear stream keeps the value the downstream units used; its computed value is collected below.
        if (tears.count(us.outlets[k])) val["\x01" + us.outlets[k]] = (*out)[k];
        else val[us.outlets[k]] = (*out)[k];
      }
    }
    return true;
  };

  const std::size_t nt = res.tears.size();
  std::vector<std::vector<double>> xprev(nt), gprev(nt);
  for (int it = 1; it <= std::max(1, spec.max_iterations); ++it) {
    res.iterations = it;
    if (!pass()) {
      merge_prefixed(last, diag, "");
      std::string unguessed;
      for (std::size_t k = 0; k < nt; ++k)
        if (!guessed[k]) unguessed += (unguessed.empty() ? "" : ", ") + res.tears[k];
      if (!unguessed.empty())
        diag.warn(DiagCode::FlowsheetTopology, "the iteration started from empty tear streams (" + unguessed +
                                                   "); a 'guess' line for them may give the units a feasible start");
      return std::nullopt;
    }
    if (nt == 0) {
      res.converged = true;
      break;
    }
    std::vector<Stream> g;
    double r = 0.0;
    for (std::size_t k = 0; k < nt; ++k) {
      g.push_back(val["\x01" + res.tears[k]]);
      r = std::max(r, tear_change(x[k], g[k]));
    }
    res.residuals.push_back(r);
    if (r < spec.tolerance) {
      res.converged = true;
      break;
    }
    if (it == spec.max_iterations) break;
    for (std::size_t k = 0; k < nt; ++k) {
      std::vector<double> xv = x[k].n, gv = g[k].n;
      xv.push_back(x[k].H_W);
      gv.push_back(g[k].H_W);
      std::vector<double> nv = gv;
      const bool accelerate = spec.method == TearMethod::Wegstein && it > 2 && xprev[k].size() == xv.size();
      if (accelerate) {
        for (std::size_t j = 0; j < xv.size(); ++j) {
          const double dx = xv[j] - xprev[k][j], dg = gv[j] - gprev[k][j];
          double q = 0.0;
          if (std::fabs(dx) > 1e-12 * std::max({std::fabs(xv[j]), std::fabs(gv[j]), 1e-300})) {
            const double s = dg / dx;
            if (s != 1.0) q = std::clamp(s / (s - 1.0), -5.0, 0.0);
          }
          nv[j] = q * xv[j] + (1.0 - q) * gv[j];
        }
      }
      xprev[k] = xv;
      gprev[k] = gv;
      if (!accelerate) {
        x[k] = g[k];
        continue;
      }
      std::vector<double> n(nv.begin(), nv.end() - 1);
      double F = 0;
      for (double& v : n) F += (v = std::max(v, 0.0));
      if (!(F > 0)) {
        x[k] = empty_stream(nc, g[k].T_K, g[k].P_Pa);
        continue;
      }
      Diagnostics q;
      auto s = stream_ph(c, g[k].P_Pa, nv.back(), n, q, g[k].T_K);
      x[k] = s ? *s : g[k];  // an accelerated point outside the model's range falls back to direct substitution
    }
  }
  merge_prefixed(last, diag, "");
  for (std::size_t k = 0; k < nt; ++k)
    if (res.converged && is_empty(x[k]))
      diag.warn(DiagCode::FlowsheetTopology, "tear stream " + res.tears[k] +
                                                 " carries no flow: no feed reaches its loop (check the stream order of each unit)");
  if (!res.converged)
    diag.error(DiagCode::RecycleNotConverged,
               "the tear streams did not converge in " + std::to_string(res.iterations) + " iterations (last change " +
                   std::to_string(res.residuals.empty() ? 0.0 : res.residuals.back()) + ")");

  // Results.
  for (const auto& nm : names) {
    NamedStream ns;
    ns.name = nm;
    const Port src = source[nm];
    ns.source = src.unit < 0 ? "feed" : spec.units[static_cast<std::size_t>(src.unit)].name;
    if (dest.count(nm)) ns.destination = spec.units[static_cast<std::size_t>(dest[nm].unit)].name;
    ns.stream = val[nm];
    res.streams.push_back(std::move(ns));
  }
  res.units = reports;
  auto mass = [&](const Stream& s) {
    double m = 0;
    for (std::size_t i = 0; i < nc; ++i) m += s.n[i] * c.pure[i].molar_mass * 1e-3;
    return m;
  };
  double E = 0;
  for (const auto& ns : res.streams) {
    if (ns.source == "feed") res.mass_in_kgs += mass(ns.stream), E += ns.stream.H_W;
    if (ns.destination.empty()) res.mass_out_kgs += mass(ns.stream), E -= ns.stream.H_W;
  }
  for (const auto& u : res.units) E += u.duty_W + u.work_W;
  res.energy_imbalance_W = E;
  return res;
}

namespace {
std::string jnum(double v) {
  char b[40];
  std::snprintf(b, sizeof b, "%.12g", v);
  return std::isfinite(v) ? b : "null";
}
std::string jarr(const std::vector<double>& v) {
  std::string o = "[";
  for (std::size_t i = 0; i < v.size(); ++i) o += (i ? "," : "") + jnum(v[i]);
  return o + "]";
}
std::string jstrs(const std::vector<std::string>& v) {
  std::string o = "[";
  for (std::size_t i = 0; i < v.size(); ++i) o += (i ? "," : "") + json_escape(v[i]);
  return o + "]";
}
}  // namespace

std::string flowsheet_json(std::string_view text) {
  Diagnostics d;
  auto spec = parse_flowsheet(text, d);
  std::optional<FlowsheetResult> r;
  if (spec) r = solve_flowsheet(*spec, d);
  std::string o = std::string("{\"ok\":") + (r && r->converged && !d.has_errors() ? "true" : "false");
  if (r) {
    const auto& c = r->components;
    o += ",\"converged\":" + std::string(r->converged ? "true" : "false") + ",\"iterations\":" +
         std::to_string(r->iterations) + ",\"method\":" +
         json_escape(spec->method == TearMethod::Direct ? "direct" : "wegstein") + ",\"tears\":" + jstrs(r->tears) +
         ",\"order\":" + jstrs(r->order) + ",\"residuals\":" + jarr(r->residuals) + ",\"components\":[";
    for (std::size_t i = 0; i < c.size(); ++i)
      o += (i ? "," : "") + std::string("{\"smiles\":") + json_escape(c.smiles[i]) + ",\"name\":" +
           json_escape(c.pure[i].name) + ",\"molar_mass\":" + jnum(c.pure[i].molar_mass) + "}";
    o += "],\"streams\":[";
    for (std::size_t k = 0; k < r->streams.size(); ++k) {
      const auto& ns = r->streams[k];
      const Stream& s = ns.stream;
      double m = 0;
      for (std::size_t i = 0; i < c.size(); ++i) m += s.n[i] * c.pure[i].molar_mass * 1e-3;
      o += (k ? "," : "") + std::string("{\"name\":") + json_escape(ns.name) + ",\"from\":" + json_escape(ns.source) +
           ",\"to\":" + json_escape(ns.destination) + ",\"T_K\":" + jnum(s.T_K) + ",\"P_Pa\":" + jnum(s.P_Pa) +
           ",\"vapor_fraction\":" + jnum(s.vapor_fraction) + ",\"H_W\":" + jnum(s.H_W) + ",\"mass_kgs\":" + jnum(m) +
           ",\"flows\":" + jarr(s.n) + ",\"x\":" + jarr(s.x) + ",\"y\":" + jarr(s.y) + "}";
    }
    o += "],\"units\":[";
    for (std::size_t k = 0; k < r->units.size(); ++k) {
      const auto& u = r->units[k];
      o += (k ? "," : "") + std::string("{\"name\":") + json_escape(u.name) + ",\"type\":" + json_escape(to_string(u.kind)) +
           ",\"duty_W\":" + jnum(u.duty_W) + ",\"work_W\":" + jnum(u.work_W) + ",\"values\":{";
      for (std::size_t j = 0; j < u.values.size(); ++j)
        o += (j ? "," : "") + json_escape(u.values[j].first) + ":" + jnum(u.values[j].second);
      o += "}}";
    }
    o += "],\"balance\":{\"mass_in_kgs\":" + jnum(r->mass_in_kgs) + ",\"mass_out_kgs\":" + jnum(r->mass_out_kgs) +
         ",\"energy_imbalance_W\":" + jnum(r->energy_imbalance_W) + "}";
  }
  return o + ",\"diagnostics\":" + diagnostics_to_json(d) + "}";
}

}  // namespace chemlarp
