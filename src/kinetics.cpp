#include "omnichem/kinetics.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <sstream>
#include <utility>

#include "omnichem/elements.hpp"
#include "omnichem/smiles.hpp"

namespace omnichem {

double ArrheniusRate::k(double T) const noexcept {
  return A * std::pow(T, b) * std::exp(-Ea_J / (kGasConstant * T));
}

std::int32_t Mechanism::find(std::string_view name) const {
  for (std::size_t i = 0; i < names.size(); ++i) {
    if (names[i].size() != name.size()) continue;
    bool same = true;
    for (std::size_t c = 0; c < name.size() && same; ++c)
      same = std::toupper(static_cast<unsigned char>(names[i][c])) == std::toupper(static_cast<unsigned char>(name[c]));
    if (same) return static_cast<std::int32_t>(i);
  }
  return -1;
}

namespace {

// ---- parsing -----------------------------------------------------------------

struct Tok {
  std::string text;
  std::size_t pos;  // offset in the whole file
};

std::vector<Tok> split_ws(std::string_view line, std::size_t base) {
  std::vector<Tok> out;
  std::size_t i = 0;
  while (i < line.size()) {
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
    if (i >= line.size()) break;
    const std::size_t s = i;
    while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i]))) ++i;
    out.push_back({std::string(line.substr(s, i - s)), base + s});
  }
  return out;
}

std::string lower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

bool parse_number(const std::string& s, double& v) {
  if (s.empty()) return false;
  char* end = nullptr;
  v = std::strtod(s.c_str(), &end);
  return end == s.c_str() + s.size() && std::isfinite(v);
}

void merge(std::vector<std::pair<std::uint32_t, double>>& side) {
  std::sort(side.begin(), side.end());
  std::vector<std::pair<std::uint32_t, double>> out;
  for (auto& p : side) {
    if (!out.empty() && out.back().first == p.first) out.back().second += p.second;
    else out.push_back(p);
  }
  side = std::move(out);
}

double int_pow(double c, double nu, bool& integral) {
  const double r = std::round(nu);
  integral = std::fabs(nu - r) < 1e-12 && r >= 0;
  if (integral) {
    double p = 1.0;
    for (int i = 0; i < static_cast<int>(r); ++i) p *= c;
    return p;
  }
  return std::pow(std::max(c, 0.0), nu);
}

double power(double c, double nu) {
  bool integral;
  return int_pow(c, nu, integral);
}

// d/dc c^nu
double dpower(double c, double nu) {
  bool integral;
  int_pow(c, nu, integral);
  if (integral) return nu == 0.0 ? 0.0 : nu * power(c, nu - 1.0);
  return c > 0.0 ? nu * std::pow(c, nu - 1.0) : 0.0;
}

}  // namespace

std::optional<Mechanism> parse_mechanism(std::string_view text, Diagnostics& diag) {
  Mechanism m;
  enum class Section { None, Species, Reactions } section = Section::None;
  double ea_factor = 4.184;     // cal/mol default (CHEMKIN)
  double volume_factor = 1e-6;  // cm3 -> m3 default (CHEMKIN)
  const std::size_t errors_before = diag.error_count();
  std::size_t line_no = 0, offset = 0;

  while (offset <= text.size()) {
    std::size_t eol = text.find('\n', offset);
    if (eol == std::string_view::npos) eol = text.size();
    std::string_view line = text.substr(offset, eol - offset);
    const std::size_t line_base = offset;
    offset = eol + 1;
    ++line_no;
    // Comments: '!' anywhere, '#' only at line start or after whitespace ('#' is the SMILES triple bond).
    std::size_t cut = line.find('!');
    for (std::size_t c = 0; c < line.size() && c < cut; ++c)
      if (line[c] == '#' && (c == 0 || std::isspace(static_cast<unsigned char>(line[c - 1])))) { cut = c; break; }
    if (cut != std::string_view::npos) line = line.substr(0, cut);
    auto toks = split_ws(line, line_base);
    if (toks.empty()) {
      if (eol == text.size()) break;
      continue;
    }
    const std::string head = lower(toks[0].text);
    const std::string where = "line " + std::to_string(line_no) + ": ";

    if (head == "end") {
      section = Section::None;
      continue;
    }
    if (section == Section::None) {
      if (head == "species") {
        section = Section::Species;
        if (toks.size() > 1) diag.error(DiagCode::MechanismSyntax, where + "put species on their own lines after 'species'").at(toks[1].pos);
      } else if (head == "reactions") {
        section = Section::Reactions;
        for (std::size_t i = 1; i < toks.size(); ++i) {
          const std::string u = lower(toks[i].text);
          if (u == "cal/mol") ea_factor = 4.184;
          else if (u == "kcal/mol") ea_factor = 4184.0;
          else if (u == "j/mol") ea_factor = 1.0;
          else if (u == "kj/mol") ea_factor = 1000.0;
          else if (u == "k" || u == "kelvins") ea_factor = kGasConstant;
          else if (u == "cm3" || u == "moles") volume_factor = 1e-6;
          else if (u == "m3") volume_factor = 1.0;
          else diag.error(DiagCode::MechanismSyntax, where + "unknown unit '" + toks[i].text + "'").at(toks[i].pos);
        }
      } else {
        diag.error(DiagCode::MechanismSyntax, where + "expected 'species' or 'reactions' section").at(toks[0].pos);
      }
      continue;
    }

    if (section == Section::Species) {
      if (toks.size() != 2) {
        diag.error(DiagCode::MechanismSyntax, where + "species lines are 'NAME SMILES'").at(toks[0].pos);
        continue;
      }
      if (m.find(toks[0].text) >= 0) {
        diag.error(DiagCode::DuplicateSpecies, where + "species '" + toks[0].text + "' is declared twice").at(toks[0].pos);
        continue;
      }
      const std::int32_t sp = parse_smiles(m.state, toks[1].text, diag, toks[1].pos);
      if (sp < 0) continue;
      for (std::size_t i = 0; i + 1 < m.state.species.size(); ++i) {
        if (m.state.species[i].canonical_smiles == m.state.species.back().canonical_smiles) {
          diag.error(DiagCode::DuplicateSpecies, where + "'" + toks[0].text + "' has the same structure as '" + m.names[i] +
                                                     "' (" + m.state.species[i].canonical_smiles + ")")
              .at(toks[1].pos);
        }
      }
      m.names.push_back(toks[0].text);
      continue;
    }

    // ---- reactions section ----
    const bool has_arrow = line.find('=') != std::string_view::npos;
    if (!has_arrow) {
      // Third-body efficiencies: "H2/2.5/ H2O/12/" or "efficiencies: H2:2.5 H2O:12"
      if (m.reactions.empty() || !m.reactions.back().third_body) {
        diag.error(DiagCode::MechanismSyntax, where + "efficiencies must follow a third-body (+ M) reaction").at(toks[0].pos);
        continue;
      }
      std::string body(line);
      std::size_t start = 0;
      if (lower(toks[0].text) == "efficiencies:") start = body.find(':') + 1;
      for (char& ch : body) if (ch == '/' || ch == ':') ch = ' ';
      auto parts = split_ws(std::string_view(body).substr(start), line_base + start);
      if (parts.size() % 2 != 0) {
        diag.error(DiagCode::MechanismSyntax, where + "efficiencies come in NAME value pairs").at(toks[0].pos);
        continue;
      }
      for (std::size_t i = 0; i < parts.size(); i += 2) {
        const std::int32_t sp = m.find(parts[i].text);
        double v = 0;
        if (sp < 0) {
          diag.error(DiagCode::UnknownSpecies, where + "unknown species '" + parts[i].text + "' in efficiencies").at(parts[i].pos);
        } else if (!parse_number(parts[i + 1].text, v) || v < 0) {
          diag.error(DiagCode::InvalidRateParameter, where + "efficiency must be a number >= 0").at(parts[i + 1].pos);
        } else {
          m.reactions.back().efficiency[static_cast<std::size_t>(sp)] = v;
        }
      }
      continue;
    }

    std::size_t numeric_tail = 0;
    for (std::size_t i = toks.size(); i-- > 0;) {
      double v;
      if (!parse_number(toks[i].text, v)) break;
      ++numeric_tail;
    }
    if (toks.size() < 4 || numeric_tail < 3) {
      diag.error(DiagCode::MechanismSyntax, where + "reaction lines are 'EQUATION A b Ea' (three rate parameters)").at(toks[0].pos);
      continue;
    }
    double A = 0, b = 0, Ea = 0;
    const Tok& tA = toks[toks.size() - 3];
    const Tok& tb = toks[toks.size() - 2];
    const Tok& tE = toks[toks.size() - 1];
    bool params_ok = true;
    if (!parse_number(tA.text, A) || A < 0) {
      diag.error(DiagCode::InvalidRateParameter, where + "pre-exponential factor A must be a finite number >= 0").at(tA.pos);
      params_ok = false;
    }
    if (!parse_number(tb.text, b)) {
      diag.error(DiagCode::InvalidRateParameter, where + "temperature exponent b must be a finite number").at(tb.pos);
      params_ok = false;
    }
    if (!parse_number(tE.text, Ea)) {
      diag.error(DiagCode::InvalidRateParameter, where + "activation energy Ea must be a finite number").at(tE.pos);
      params_ok = false;
    }

    std::string eq;
    for (std::size_t i = 0; i + 3 < toks.size(); ++i) eq += toks[i].text;
    const std::size_t eq_pos = toks[0].pos;
    KineticReaction r;
    r.line = line_no;
    if (eq.find("(+") != std::string::npos) {
      diag.error(DiagCode::UnsupportedReactionType, where + "pressure-dependent fall-off '(+M)' reactions are not supported yet").at(eq_pos);
      continue;
    }
    std::size_t arrow = eq.find("<=>");
    std::size_t arrow_len = 3;
    r.reversible = true;
    if (arrow == std::string::npos) {
      arrow = eq.find("=>");
      arrow_len = 2;
      r.reversible = false;
      if (arrow == std::string::npos) {
        arrow = eq.find('=');
        arrow_len = 1;
        r.reversible = true;
      }
    }
    const std::string lhs = eq.substr(0, arrow), rhs = eq.substr(arrow + arrow_len);
    bool ok = params_ok;
    int third_lhs = 0, third_rhs = 0;
    auto parse_side = [&](const std::string& side, std::vector<std::pair<std::uint32_t, double>>& out, int& third) {
      if (side.empty()) {
        diag.error(DiagCode::EmptyReactionSide, where + "reaction side is empty").at(eq_pos);
        ok = false;
        return;
      }
      std::size_t s = 0;
      while (s <= side.size()) {
        std::size_t e = side.find('+', s);
        if (e == std::string::npos) e = side.size();
        std::string term = side.substr(s, e - s);
        s = e + 1;
        if (term.empty()) {
          diag.error(DiagCode::MechanismSyntax, where + "empty term (stray '+')").at(eq_pos);
          ok = false;
          if (e == side.size()) break;
          continue;
        }
        std::size_t k = 0;
        while (k < term.size() && (std::isdigit(static_cast<unsigned char>(term[k])) || term[k] == '.')) ++k;
        double coef = 1.0;
        if (k > 0 && (!parse_number(term.substr(0, k), coef) || coef <= 0)) {
          diag.error(DiagCode::InvalidCoefficient, where + "invalid coefficient in '" + term + "'").at(eq_pos);
          ok = false;
          coef = 1.0;
        }
        const std::string name = term.substr(k);
        if (name == "M" || name == "m") {
          ++third;
        } else {
          const std::int32_t sp = m.find(name);
          if (sp < 0) {
            diag.error(DiagCode::UnknownSpecies, where + "unknown species '" + name + "'").at(eq_pos);
            ok = false;
          } else {
            out.push_back({static_cast<std::uint32_t>(sp), coef});
          }
        }
        if (e == side.size()) break;
      }
    };
    parse_side(lhs, r.reactants, third_lhs);
    parse_side(rhs, r.products, third_rhs);
    if (third_lhs != third_rhs || third_lhs > 1) {
      diag.error(DiagCode::MechanismSyntax, where + "third body M must appear exactly once on each side").at(eq_pos);
      ok = false;
    }
    if (!ok) continue;
    merge(r.reactants);
    merge(r.products);
    r.third_body = third_lhs == 1;
    if (r.third_body) r.efficiency.assign(m.species_count(), 1.0);
    double order = r.third_body ? 1.0 : 0.0;
    for (auto& [s2, nu] : r.reactants) order += nu;
    r.forward = {A * std::pow(volume_factor, order - 1.0), b, Ea * ea_factor};
    r.equation = eq;

    // Element and charge balance from the SMILES-derived compositions.
    std::array<double, kElementSlots> delta{};
    double dq = 0;
    for (auto& [s2, nu] : r.reactants) {
      for (std::size_t z = 1; z < kElementSlots; ++z) delta[z] -= nu * m.state.species[s2].composition[z];
      dq -= nu * m.state.species[s2].net_charge;
    }
    for (auto& [s2, nu] : r.products) {
      for (std::size_t z = 1; z < kElementSlots; ++z) delta[z] += nu * m.state.species[s2].composition[z];
      dq += nu * m.state.species[s2].net_charge;
    }
    bool balanced = true;
    for (std::size_t z = 1; z < kElementSlots; ++z) {
      if (std::fabs(delta[z]) < 1e-9) continue;
      balanced = false;
      char buf[48];
      std::snprintf(buf, sizeof buf, "%+g", delta[z]);
      diag.error(DiagCode::ElementImbalance, where + std::string(element_symbol(static_cast<std::uint8_t>(z))) +
                                                 " is not conserved in " + eq + " (products - reactants = " + buf + ")")
          .at(eq_pos)
          .with_element(static_cast<std::uint8_t>(z), delta[z]);
    }
    if (std::fabs(dq) > 1e-9) {
      balanced = false;
      diag.error(DiagCode::ChargeImbalance, where + "charge is not conserved in " + eq).at(eq_pos).with_element(0, dq);
    }
    if (balanced) m.reactions.push_back(std::move(r));
    if (eol == text.size()) break;
  }
  if (section != Section::None)
    diag.error(DiagCode::MechanismSyntax, "missing 'end' after the last section").at(text.size());
  if (m.names.empty()) diag.error(DiagCode::MechanismSyntax, "mechanism declares no species").at(0);

  // Thermo: required for every species taking part in a reversible reaction.
  m.thermo.resize(m.species_count());
  std::vector<char> needs(m.species_count(), 0);
  for (const auto& r : m.reactions) {
    if (!r.reversible) continue;
    for (auto& [s2, nu] : r.reactants) needs[s2] = 1;
    for (auto& [s2, nu] : r.products) needs[s2] = 1;
  }
  for (std::uint32_t i = 0; i < m.species_count(); ++i) {
    Diagnostics local;
    m.thermo[i] = resolve_thermo(m.state, i, kTref, local);
    if (!needs[i]) continue;
    for (const auto& d : local.items()) {
      if (d.code == DiagCode::NoJobackGroup) continue;          // summarised by the NoThermoData error
      if (d.code == DiagCode::TemperatureOutOfRange) continue;  // checked against the actual run temperatures
      (d.severity == Severity::Error ? diag.error(d.code, m.names[i] + ": " + d.message)
                                     : diag.warn(d.code, m.names[i] + ": " + d.message))
          .on_species(i);
    }
  }
  if (diag.error_count() != errors_before) return std::nullopt;
  return m;
}

std::optional<Mechanism> load_mechanism(const std::string& path, Diagnostics& diag) {
  std::ifstream in(path);
  if (!in) {
    diag.error(DiagCode::InvalidInput, "cannot open mechanism file '" + path + "'");
    return std::nullopt;
  }
  std::stringstream ss;
  ss << in.rdbuf();
  return parse_mechanism(ss.str(), diag);
}

// ---- rates -------------------------------------------------------------------

namespace {

struct RateTerms {
  double kf, kr, M, fwd, rev;  // q = M (kf fwd - kr rev)
};

double log_kc(const Mechanism& m, const KineticReaction& r, double T) {
  double dg = 0.0, dnu = 0.0;
  for (auto& [s, nu] : r.products) { dg += nu * m.thermo[s]->g_RT(T); dnu += nu; }
  for (auto& [s, nu] : r.reactants) { dg -= nu * m.thermo[s]->g_RT(T); dnu -= nu; }
  // K_c = K_p (P°/RT)^Δν
  return -dg + dnu * std::log(kStandardPressure / (kGasConstant * T));
}

RateTerms terms(const Mechanism& m, const KineticReaction& r, double T, const double* C) {
  RateTerms t{};
  t.kf = r.forward.k(T);
  t.kr = r.reversible ? t.kf * std::exp(std::min(-log_kc(m, r, T), 700.0)) : 0.0;
  t.M = 1.0;
  if (r.third_body) {
    t.M = 0.0;
    for (std::size_t i = 0; i < m.species_count(); ++i) t.M += r.efficiency[i] * C[i];
  }
  t.fwd = 1.0;
  for (auto& [s, nu] : r.reactants) t.fwd *= power(C[s], nu);
  t.rev = 1.0;
  if (r.reversible)
    for (auto& [s, nu] : r.products) t.rev *= power(C[s], nu);
  return t;
}

}  // namespace

double equilibrium_constant_c(const Mechanism& m, std::size_t r, double T) {
  return std::exp(log_kc(m, m.reactions[r], T));
}

void rates_of_progress(const Mechanism& m, double T, const double* C, double* q) {
  for (std::size_t r = 0; r < m.reactions.size(); ++r) {
    const RateTerms t = terms(m, m.reactions[r], T, C);
    q[r] = t.M * (t.kf * t.fwd - t.kr * t.rev);
  }
}

void production_rates(const Mechanism& m, double T, const double* C, double* wdot) {
  std::fill(wdot, wdot + m.species_count(), 0.0);
  for (const auto& r : m.reactions) {
    const RateTerms t = terms(m, r, T, C);
    const double q = t.M * (t.kf * t.fwd - t.kr * t.rev);
    for (auto& [s, nu] : r.reactants) wdot[s] -= nu * q;
    for (auto& [s, nu] : r.products) wdot[s] += nu * q;
  }
}

void production_rates_jacobian(const Mechanism& m, double T, const double* C, double* J) {
  const std::size_t S = m.species_count();
  std::fill(J, J + S * S, 0.0);
  std::vector<double> dq(S);
  for (const auto& r : m.reactions) {
    const RateTerms t = terms(m, r, T, C);
    std::fill(dq.begin(), dq.end(), 0.0);
    // d(fwd)/dC_j and d(rev)/dC_j by the product rule over distinct species.
    for (auto& [j, nuj] : r.reactants) {
      double d = dpower(C[j], nuj);
      for (auto& [k, nuk] : r.reactants)
        if (k != j) d *= power(C[k], nuk);
      dq[j] += t.M * t.kf * d;
    }
    if (r.reversible) {
      for (auto& [j, nuj] : r.products) {
        double d = dpower(C[j], nuj);
        for (auto& [k, nuk] : r.products)
          if (k != j) d *= power(C[k], nuk);
        dq[j] -= t.M * t.kr * d;
      }
    }
    if (r.third_body) {
      const double base = t.kf * t.fwd - t.kr * t.rev;
      for (std::size_t j = 0; j < S; ++j) dq[j] += r.efficiency[j] * base;
    }
    for (auto& [i, nu] : r.reactants)
      for (std::size_t j = 0; j < S; ++j) J[i * S + j] -= nu * dq[j];
    for (auto& [i, nu] : r.products)
      for (std::size_t j = 0; j < S; ++j) J[i * S + j] += nu * dq[j];
  }
}

// ---- reactors ----------------------------------------------------------------

const char* to_string(ReactorModel model) noexcept {
  switch (model) {
    case ReactorModel::IsothermalConstantVolume: return "iso-v";
    case ReactorModel::IsothermalConstantPressure: return "iso-p";
    case ReactorModel::AdiabaticConstantVolume: return "adiabatic-v";
    case ReactorModel::AdiabaticConstantPressure: return "adiabatic-p";
  }
  return "?";
}

std::optional<ReactorModel> reactor_model_from_string(std::string_view s) noexcept {
  if (s == "iso-v") return ReactorModel::IsothermalConstantVolume;
  if (s == "iso-p") return ReactorModel::IsothermalConstantPressure;
  if (s == "adiabatic-v") return ReactorModel::AdiabaticConstantVolume;
  if (s == "adiabatic-p") return ReactorModel::AdiabaticConstantPressure;
  return std::nullopt;
}

ReactorResult run_reactor(const Mechanism& m, const ReactorSetup& su, Diagnostics& diag) {
  ReactorResult res;
  const std::size_t S = m.species_count();
  const bool adiabatic =
      su.model == ReactorModel::AdiabaticConstantVolume || su.model == ReactorModel::AdiabaticConstantPressure;
  const bool const_v =
      su.model == ReactorModel::IsothermalConstantVolume || su.model == ReactorModel::AdiabaticConstantVolume;

  bool ok = true;
  if (su.moles.size() != S) {
    diag.error(DiagCode::InvalidReactorSetup, "initial moles must have one entry per mechanism species (" +
                                                  std::to_string(S) + ")");
    return res;
  }
  double N0 = 0.0;
  for (std::size_t i = 0; i < S; ++i) {
    if (!(su.moles[i] >= 0.0) || !std::isfinite(su.moles[i])) {
      diag.error(DiagCode::NegativeAmount, "initial amount of " + m.names[i] + " must be finite and >= 0").on_species(i);
      ok = false;
    }
    N0 += su.moles[i];
  }
  if (!(su.temperature_K > 0) || !(su.pressure_Pa > 0) || !(su.t_end > 0) || !(N0 > 0)) {
    diag.error(DiagCode::InvalidReactorSetup, "temperature, pressure, t_end and total moles must all be positive");
    ok = false;
  }
  if (adiabatic) {
    for (std::size_t i = 0; i < S; ++i) {
      if (!m.thermo[i] || !std::isfinite(m.thermo[i]->cp_R(su.temperature_K))) {
        diag.error(DiagCode::MissingHeatCapacity,
                   m.names[i] + " has no heat-capacity data; adiabatic reactors need NASA-7 thermo for every species")
            .on_species(i);
        ok = false;
      }
    }
  }
  if (!ok) return res;

  const double V0 = N0 * kGasConstant * su.temperature_K / su.pressure_Pa;
  const std::size_t n = S + (adiabatic ? 1 : 0);
  auto volume = [&](const std::vector<double>& y, double T) {
    if (const_v) return V0;
    double N = 0.0;
    for (std::size_t i = 0; i < S; ++i) N += y[i];
    return std::max(N, 1e-300) * kGasConstant * T / su.pressure_Pa;
  };
  std::vector<double> C(S), wdot(S);
  OdeRhs rhs = [&](double, const std::vector<double>& y, std::vector<double>& dy) {
    const double T = adiabatic ? y[S] : su.temperature_K;
    const double V = volume(y, T);
    for (std::size_t i = 0; i < S; ++i) C[i] = y[i] / V;
    production_rates(m, T, C.data(), wdot.data());
    for (std::size_t i = 0; i < S; ++i) dy[i] = V * wdot[i];
    if (adiabatic) {
      double num = 0.0, den = 0.0;
      for (std::size_t i = 0; i < S; ++i) {
        const auto& th = *m.thermo[i];
        if (const_v) {
          num += (th.h_RT(T) - 1.0) * kGasConstant * T * dy[i];    // u_i dn_i
          den += y[i] * (th.cp_R(T) - 1.0) * kGasConstant;          // n_i cv_i
        } else {
          num += th.h_RT(T) * kGasConstant * T * dy[i];             // h_i dn_i
          den += y[i] * th.cp_R(T) * kGasConstant;                  // n_i cp_i
        }
      }
      dy[S] = -num / den;
    }
  };
  std::vector<double> Jc(S * S), f0(n), f1(n), yT;
  OdeJacobian jac = [&](double t, const std::vector<double>& y, std::vector<double>& J) {
    // Species block: exact for constant volume (dn/dt = V w(n/V) => d/dn = dw/dC),
    // a good modified-Newton approximation at constant pressure.
    const double T = adiabatic ? y[S] : su.temperature_K;
    const double V = volume(y, T);
    for (std::size_t i = 0; i < S; ++i) C[i] = y[i] / V;
    production_rates_jacobian(m, T, C.data(), Jc.data());
    J.assign(n * n, 0.0);
    for (std::size_t i = 0; i < S; ++i)
      for (std::size_t j = 0; j < S; ++j) J[i * n + j] = Jc[i * S + j];
    if (!adiabatic) return;
    // Temperature column: one well-scaled forward difference (captures k(T), K_c(T), V(T)).
    rhs(t, y, f0);
    yT = y;
    const double dT = 1e-6 * std::max(T, 1.0);
    yT[S] += dT;
    rhs(t, yT, f1);
    for (std::size_t i = 0; i < n; ++i) J[i * n + S] = (f1[i] - f0[i]) / dT;
    // Temperature row, analytic: Tdot = -num/den, num = sum e_i ndot_i, den = sum n_i c_i,
    // with (e, c) = (u, cv) at constant volume and (h, cp) at constant pressure.
    double num = 0.0, den = 0.0;
    std::vector<double> e(S), cap(S);
    for (std::size_t i = 0; i < S; ++i) {
      const auto& th = *m.thermo[i];
      e[i] = (th.h_RT(T) - (const_v ? 1.0 : 0.0)) * kGasConstant * T;
      cap[i] = (th.cp_R(T) - (const_v ? 1.0 : 0.0)) * kGasConstant;
      num += e[i] * f0[i];
      den += y[i] * cap[i];
    }
    for (std::size_t j = 0; j < S; ++j) {
      double dnum = 0.0;
      for (std::size_t i = 0; i < S; ++i) dnum += e[i] * J[i * n + j];
      J[S * n + j] = -dnum / den + num * cap[j] / (den * den);
    }
  };

  auto snapshot = [&](double t, const std::vector<double>& y) {
    ReactorState s;
    s.t = t;
    s.T = adiabatic ? y[S] : su.temperature_K;
    s.V = volume(y, s.T);
    s.n.assign(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(S));
    double N = 0.0;
    for (double v : s.n) N += v;
    s.P = const_v ? N * kGasConstant * s.T / V0 : su.pressure_Pa;
    return s;
  };

  std::vector<double> y(su.moles);
  if (adiabatic) y.push_back(su.temperature_K);
  BdfOptions opt;
  opt.rtol = su.rtol;
  opt.atol.assign(n, su.atol * N0);
  if (adiabatic) opt.atol[S] = 1e-6;
  opt.max_steps = su.max_steps;

  res.samples.push_back(snapshot(0.0, y));
  double t_prev = 0.0, T_prev = su.temperature_K;
  OdeObserver on_step = [&](double t, const std::vector<double>& yy) {
    if (adiabatic) {
      const double rate = (yy[S] - T_prev) / (t - t_prev);
      if (rate > res.max_dTdt) {
        res.max_dTdt = rate;
        res.ignition_time = 0.5 * (t + t_prev);
      }
      T_prev = yy[S];
    }
    t_prev = t;
    if (su.output_times.empty()) res.samples.push_back(snapshot(t, yy));
  };
  OdeObserver on_output = [&](double t, const std::vector<double>& yy) { res.samples.push_back(snapshot(t, yy)); };

  res.ok = integrate_bdf(rhs, &jac, 0.0, y, su.t_end, su.output_times, on_output, &on_step, opt,
                         res.stats, diag);
  res.final_state = snapshot(su.t_end, y);
  if (adiabatic && res.ok) {
    double tmin = su.temperature_K, tmax = su.temperature_K;
    for (const auto& s : res.samples) { tmin = std::min(tmin, s.T); tmax = std::max(tmax, s.T); }
    tmax = std::max(tmax, res.final_state.T);
    for (std::size_t i = 0; i < S; ++i) {
      if (!m.thermo[i]->in_range(tmin) || !m.thermo[i]->in_range(tmax)) {
        diag.warn(DiagCode::TemperatureOutOfRange, m.names[i] + ": temperature left the thermo data range (" +
                                                       std::to_string(tmin) + "-" + std::to_string(tmax) + " K)")
            .on_species(i);
      }
    }
  }
  if (res.max_dTdt <= 0.0) res.ignition_time = -1.0;
  return res;
}

}  // namespace omnichem
